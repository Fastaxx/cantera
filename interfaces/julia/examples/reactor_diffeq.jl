#=
Custom Reactor Model with DifferentialEquations.jl
==================================================

Writes the equations of a constant-pressure, adiabatic reactor in Julia, takes the
thermodynamic and kinetic properties from Cantera, and integrates the model with
OrdinaryDiffEq.jl. The result is checked against Cantera's own reactor network, and
the model is then extended with a heat-loss term.

Requires OrdinaryDiffEq, installed in an environment of its own so that it is not
added to the Cantera package::

    julia> using Pkg
    julia> Pkg.activate(temp=true)
    julia> Pkg.develop(path="interfaces/julia")
    julia> Pkg.add("OrdinaryDiffEq")
    julia> include("interfaces/julia/examples/reactor_diffeq.jl")

.. tags:: Julia, combustion, reactor network, ignition delay, user-defined model
=#

using Cantera
using OrdinaryDiffEq

# ## Model
#
# The state is `u = [T; Y]`, temperature and mass fractions. At constant pressure,
#
#     dY_k/dt = ω̇_k W_k / ρ
#     dT/dt   = -(Σ_k h_k ω̇_k) / (ρ c_p) - q_loss / (ρ c_p)
#
# with ω̇_k the net molar production rates, W_k the molecular weights and h_k the
# partial molar enthalpies. `q_loss` is a volumetric heat loss, zero for an
# adiabatic reactor.
#
# Everything the right-hand side needs is allocated once, in the parameters, and
# filled in place: `set_TPY!` reads the mass fractions straight from a view of `u`,
# and the `!` getters write into preallocated buffers. The function therefore does
# not allocate, which matters when the solver calls it thousands of times.

struct ReactorModel{G}
    gas::G                  # Cantera phase, used to evaluate properties
    P::Float64              # pressure [Pa]
    W::Vector{Float64}      # molecular weights [kg/kmol]
    wdot::Vector{Float64}   # buffer: net production rates [kmol/m^3/s]
    h::Vector{Float64}      # buffer: partial molar enthalpies [J/kmol]
    hA_V::Float64           # heat-loss coefficient per unit volume [W/m^3/K]
    T_wall::Float64         # wall temperature [K]
end

function ReactorModel(gas; hA_V=0.0, T_wall=300.0)
    ns = n_species(gas)
    return ReactorModel(gas, pressure(gas), molecular_weights(gas), zeros(ns), zeros(ns),
                        Float64(hA_V), Float64(T_wall))
end

function reactor_rhs!(du, u, m::ReactorModel, t)
    gas = m.gas
    set_TPY!(gas, u[1], m.P, @view u[2:end])
    net_production_rates!(m.wdot, gas)
    partial_molar_enthalpies!(m.h, gas)
    rho_cp = density(gas) * cp_mass(gas)

    @views du[2:end] .= m.wdot .* m.W ./ density(gas)
    heat_release = sum(m.h[k] * m.wdot[k] for k in eachindex(m.h))
    du[1] = -(heat_release + m.hA_V * (u[1] - m.T_wall)) / rho_cp
    return nothing
end

# ## Initial state

gas = Solution("gri30.yaml")
set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")
T0 = temperature(gas)
u0 = [T0; mass_fractions(gas)]
model = ReactorModel(gas)

du = similar(u0)
reactor_rhs!(du, u0, model, 0.0)
println("Allocations per right-hand-side call: ",
        @allocated(reactor_rhs!(du, u0, model, 0.0)), " bytes")

# ## Integration
#
# The problem is stiff, so it calls for an implicit solver such as `FBDF`. Such
# solvers need the Jacobian of the right-hand side. OrdinaryDiffEq computes it by
# automatic differentiation by default, which cannot see through the calls to
# Cantera's compiled library, so finite differences are requested instead.

tspan = (0.0, 1e-3)
tols = (reltol=1e-10, abstol=1e-14)
prob = ODEProblem(reactor_rhs!, u0, tspan, model)
alg = FBDF(autodiff=AutoFiniteDiff())
sol = solve(prob, alg; tols...)
println("OrdinaryDiffEq: ", length(sol.t), " steps, final T = ",
        round(sol.u[end][1], digits=2), " K")

# Ignition delay: time at which the temperature first rises 400 K above its
# initial value, located on the dense solution.
function ignition_delay(t, T, T0)
    i = findfirst(>(T0 + 400.0), T)
    i === nothing && return NaN
    # linear interpolation between the two samples that bracket the crossing
    return t[i-1] + (T0 + 400.0 - T[i-1]) * (t[i] - t[i-1]) / (T[i] - T[i-1])
end
t_fine = range(tspan...; length=20_001)
tau_diffeq = ignition_delay(t_fine, [sol(t)[1] for t in t_fine], T0)

# ## Comparison with Cantera's reactor network
#
# The same constant-pressure reactor, integrated by Cantera itself with the same
# tolerances. Both converge to the same ignition delay as the tolerances are
# tightened; at `rtol = 1e-8`, the two results would still differ by about 0.3%.

set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")
reactor = IdealGasConstPressureReactor(gas)
net = ReactorNet(reactor)
set_tolerances!(net; rtol=tols.reltol, atol=tols.abstol)
T_net = similar(t_fine)
for (i, t) in enumerate(t_fine)
    i > 1 && advance!(net, t)
    T_net[i] = temperature(reactor)
end
tau_net = ignition_delay(t_fine, T_net, T0)

println("Ignition delay: OrdinaryDiffEq ", round(tau_diffeq * 1e6, digits=2),
        " µs, ReactorNet ", round(tau_net * 1e6, digits=2), " µs")
println("Final temperature: OrdinaryDiffEq ", round(sol.u[end][1], digits=2),
        " K, ReactorNet ", round(T_net[end], digits=2), " K")

# ## Extending the model: heat losses
#
# Since the equations are plain Julia, changing the physics only takes a change to
# the right-hand side or its parameters. Here, a heat loss to a wall at 300 K
# delays ignition and lowers the temperature reached. With about three times more
# loss, `hA_V = 1e5`, the mixture would not ignite within 1 ms at all.

set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")
lossy = ReactorModel(gas; hA_V=3e4, T_wall=300.0)
sol_loss = solve(ODEProblem(reactor_rhs!, u0, tspan, lossy), alg; tols...)
tau_loss = ignition_delay(t_fine, [sol_loss(t)[1] for t in t_fine], T0)
println("With heat losses: ignition delay ", round(tau_loss * 1e6, digits=2),
        " µs, final T = ", round(sol_loss.u[end][1], digits=2), " K")
