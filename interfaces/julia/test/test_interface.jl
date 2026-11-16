# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test

const SURF_MECH = "ptcombust.yaml"
const LEAN = "CH4:0.095, O2:0.21, AR:0.695"

@testset "Interface and adjacent phases" begin
    surf = Interface(SURF_MECH, "Pt_surf")
    @test n_adjacent(surf) == 1
    @test adjacent_name(surf, 1) == "gas"
    gas = adjacent(surf, 1)
    @test n_species(gas) == n_species(adjacent(surf, "gas")) == 32
    @test n_phases(surf) == 2
    @test n_phases(gas) == 1
    @test n_total_species(surf) == n_species(surf) + n_species(gas)

    # Kinetics species span the surface species, then the gas species.
    ns = n_species(surf)
    @test [kinetics_species_index(surf, nm) for nm in species_names(surf)] == 1:ns
    @test kinetics_species_index(surf, "CH4") == ns + species_index(gas, "CH4")
    @test kinetics_species_index(surf, "not-a-species") == 0
    @test length(net_production_rates(surf)) == n_total_species(surf)

    # An existing bulk phase can be shared with the interface.
    g = Solution(SURF_MECH, "gas")
    s2 = Interface(SURF_MECH, "Pt_surf"; adjacent=[g])
    set_TPX!(g, 1000.0, one_atm, "CH4:1")
    @test temperature(adjacent(s2, 1)) ≈ 1000.0
    close!(s2); close!(g); close!(surf)
end

@testset "Coverages" begin
    surf = Interface(SURF_MECH, "Pt_surf")
    ns = n_species(surf)
    iO = species_index(surf, "O(S)")
    iPt = species_index(surf, "PT(S)")
    @test sum(coverages(surf)) ≈ 1.0

    θ = zeros(ns); θ[iPt] = 2.0; θ[iO] = 2.0
    @test set_coverages!(surf, θ) === surf
    @test coverages(surf)[iO] ≈ 0.5               # normalized
    set_coverages!(surf, "PT(S):3, O(S):1")
    @test coverages(surf)[iO] ≈ 0.25
    set_unnormalized_coverages!(surf, θ)
    @test coverages(surf)[iO] ≈ 2.0                # not normalized

    # In-place getter: same values, no allocation, contiguous views accepted.
    out = zeros(ns)
    @test coverages!(out, surf) == coverages(surf)
    f!(o, s) = coverages!(o, s)
    f!(out, surf)
    @test (@allocated f!(out, surf)) == 0
    M = zeros(ns, 2)
    coverages!(view(M, :, 2), surf)
    @test M[:, 2] == coverages(surf)
    @test_throws DimensionMismatch coverages!(zeros(ns + 1), surf)

    n0 = site_density(surf)
    @test n0 > 0
    set_site_density!(surf, 2n0)
    @test site_density(surf) ≈ 2n0

    # Coverages are only defined for surface phases.
    gas = adjacent(surf, 1)
    @test_throws CanteraError coverages(gas)
    close!(surf)
end

@testset "advance_coverages!" begin
    surf = Interface(SURF_MECH, "Pt_surf")
    gas = adjacent(surf, 1)
    set_TPX!(gas, 900.0, one_atm, LEAN)
    set_TP!(surf, 900.0, one_atm)
    set_coverages!(surf, "PT(S):1")
    @test advance_coverages!(surf, 1.0) === surf
    θ = coverages(surf)
    @test sum(θ) ≈ 1.0
    # A lean mixture leaves platinum mostly covered by oxygen.
    @test argmax(θ) == species_index(surf, "O(S)")
    @test θ[species_index(surf, "O(S)")] > 0.5
    close!(surf)
end
