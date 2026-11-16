# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Interfaces (surface phases), their adjacent bulk phases, surface coverages and
# multi-phase kinetics.

"""
    Interface(infile, name; adjacent=Solution[]) -> Solution

Construct an interface (surface) phase `name` from the input file `infile`, for
example `Interface("ptcombust.yaml", "Pt_surf")`.

`adjacent` lists the bulk phases the interface borders.  When it is empty, the
adjacent phases named in the input file are created along with the interface;
pass existing `Solution`s to share them, e.g. with a reactor.  The result is a
`Solution` whose kinetics manager covers the species of all these phases.
"""
function Interface(infile::AbstractString, name::AbstractString;
                   adjacent::AbstractVector{Solution}=Solution[])
    handles = Int32[s.handle for s in adjacent]
    h = check(LibCantera.sol_newInterface(infile, name, Int32(length(handles)), handles))
    return Solution(h)
end

"Number of bulk phases adjacent to an interface."
n_adjacent(s::Solution) = Int(check(LibCantera.sol_nAdjacent(s.handle)))

"""
    adjacent(surf, i) -> Solution
    adjacent(surf, name) -> Solution

Adjacent phase `i` (1-based) or `name` of the interface `surf`.  The result
refers to the same Cantera object as the phase used by the interface.
"""
adjacent(s::Solution, i::Integer) =
    Solution(check(LibCantera.sol_adjacent(s.handle, Int32(i - 1))))
adjacent(s::Solution, name::AbstractString) =
    Solution(check(LibCantera.sol_adjacentByName(s.handle, name)))

"Name of adjacent phase `i` (1-based) of an interface."
adjacent_name(s::Solution, i::Integer) =
    get_string((n, b) -> LibCantera.sol_adjacentName(s.handle, Int32(i - 1), n, b))

# ---- coverages --------------------------------------------------------------

"Fractional site coverages of the surface species (they sum to 1)."
coverages(g::ThermoLike) =
    get_array(n_species(g), (n, b) -> LibCantera.surf_getCoverages(_thermo_handle(g), n, b))

coverages!(out, g::ThermoLike) =
    get_array!(out, n_species(g),
        (n, b) -> LibCantera.surf_getCoverages(_thermo_handle(g), n, b))
@doc inplace_doc(:coverages, "surf") coverages!

"""
    set_coverages!(surf, θ)

Set the surface coverages from a numeric vector or a composition string
(e.g. `"PT(S):0.5, O(S):0.5"`).  The values are normalized to sum to 1.
"""
function set_coverages!(g::ThermoLike, theta::AbstractVector{<:Real})
    x = as_f64(theta)
    check(LibCantera.surf_setCoverages(_thermo_handle(g), Int32(length(x)), x))
    return g
end
function set_coverages!(g::ThermoLike, theta::AbstractString)
    check(LibCantera.surf_setCoveragesByName(_thermo_handle(g), theta))
    return g
end

"""
    set_unnormalized_coverages!(surf, θ)

Set the surface coverages without normalizing them, e.g. to evaluate rates at a
perturbed state when computing a Jacobian.
"""
function set_unnormalized_coverages!(g::ThermoLike, theta::AbstractVector{<:Real})
    x = as_f64(theta)
    check(LibCantera.surf_setCoveragesNoNorm(_thermo_handle(g), Int32(length(x)), x))
    return g
end

"Site density of the surface [kmol/m^2]."
site_density(g::ThermoLike) = checkd(LibCantera.surf_siteDensity(_thermo_handle(g)))

"Set the site density of the surface [kmol/m^2]."
function set_site_density!(g::ThermoLike, n)
    check(LibCantera.surf_setSiteDensity(_thermo_handle(g), Float64(n)))
    return g
end

"""
    advance_coverages!(surf, dt) -> surf

Integrate the surface coverages over a time `dt` [s], holding the adjacent bulk
phases at their current state.
"""
function advance_coverages!(g::KineticsLike, dt)
    check(LibCantera.kin_advanceCoverages(_kinetics_handle(g), Float64(dt)))
    return g
end

# ---- multi-phase kinetics ---------------------------------------------------

"Number of phases spanned by a kinetics manager (1 for a homogeneous phase)."
n_phases(g::KineticsLike) = Int(check(LibCantera.kin_nPhases(_kinetics_handle(g))))

"""
    kinetics_species_index(gas, name) -> Int

1-based index of species `name` in the kinetics manager's species list, which
spans all its phases, or `0` if it is not present.  This is the index into
arrays such as [`net_production_rates`](@ref) of an interface.
"""
function kinetics_species_index(g::KineticsLike, name::AbstractString)
    idx = LibCantera.kin_kineticsSpeciesIndex(_kinetics_handle(g), name)
    idx == -1 && return 0
    return Int(check(idx)) + 1
end

export Interface, n_adjacent, adjacent, adjacent_name,
       coverages, coverages!, set_coverages!, set_unnormalized_coverages!,
       site_density, set_site_density!, advance_coverages!, kinetics_species_index
