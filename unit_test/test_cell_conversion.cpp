#include "test_support.h"
#include "core/cell.h"
#include "core/mesh.h"
#include "hydro/flux.h"
#include "config/config.h"
#include "microphysics/eos/eos.h"
#include "microphysics/composition/species.h"
#include "physcon/constants.h"
#include <array>
#include <filesystem>
#include <fstream>
#include <limits>

namespace
{
using namespace hd_core;
using hd_microphysics::Equation_of_State;

std::array<double, 9> metadata(const Cell& cell)
{
    return {cell.gamma, cell.temperature, cell.n_baryon, cell.Y_e_free,
            cell.Y_Ion, cell.n_e_free, cell.n_Ion, cell.mean_molecular_weight,
            cell.dual_energy_correction};
}

void unchanged(const Cell& cell, const Cell& before)
{
    CHECK(cell.U.data == before.U.data);
    CHECK(metadata(cell) == metadata(before));
}

Cell context()
{
    Cell cell;
    cell.U[RHO] = 7.0;
    cell.U[INTERNAL_ENERGY] = 11.0;
    cell.U[ENERGY] = 11.0;
    if constexpr (NCOMP > 0) cell.U[COMP(0)] = 7.0;
    cell.temperature = 17.0;
    cell.n_baryon = 23.0;
    cell.Y_e_free = 0.12;
    cell.Y_Ion = 0.34;
    cell.n_e_free = 29.0;
    cell.n_Ion = 31.0;
    cell.mean_molecular_weight = 4.0;
    cell.dual_energy_correction = -37.0;
    for (int k = 0; k < NUSER; ++k) cell.U[CONUSR(k)] = 41.0 + k;
    return cell;
}

void check_constructor_density()
{
    const auto saved_model = hd_microphysics::eos_type;
    hd_microphysics::eos_type = Equation_of_State::IdealGas_Default;
    Primitives input{};
    input[RHO] = 1e-10;
    input[PRE] = 3.0;
    input[VEL(0)] = 2.0;
    if constexpr (NCOMP == 2)
    {
        input[COMP(0)] = 0.65;
        input[COMP(1)] = 0.35;
    }
    const Cell primitive_cell(input, 1.4);
    const Cell conserved_cell(primitive_cell.U, 1.4);
    CHECK(conserved_cell.U.data == primitive_cell.U.data);
    CHECK(near(primitive_cell.U[INTERNAL_ENERGY], 7.5));
    CHECK(near(primitive_cell.U[ENERGY], 7.5 + 2e-10));
    for (const Cell* cell : {&primitive_cell, &conserved_cell})
    {
        CHECK(cell->gamma == 1.4);
        CHECK(cell->temperature == 1.0);
        CHECK(cell->dual_energy_correction == 0.0);
        if constexpr (NCOMP > 0)
        {
            CHECK(near(cell->mean_molecular_weight, 1.0 / 1.35));
            CHECK(near(cell->Y_e_free, 0.35));
            CHECK(near(cell->Y_Ion, 1.0));
            CHECK(near(cell->n_baryon, input[RHO] / physcon::m_u));
            CHECK(near(cell->n_e_free, 0.35 * input[RHO] / physcon::m_u));
            CHECK(near(cell->n_Ion, input[RHO] / physcon::m_u));
        }
        else
        {
            Cell expected;
            expected.gamma = 1.4;
            CHECK(metadata(*cell) == metadata(expected));
        }
    }
    // A mesh's empty cells and explicit vacuum states cannot enter the density evaluator.
    for (const Cell& vacuum : {Cell{}, Cell(Primitives{}), Cell(Conserved{})})
    {
        CHECK(vacuum.U.data == Conserved{}.data);
        CHECK(metadata(vacuum) == metadata(Cell{}));
    }
    hd_microphysics::eos_type = saved_model;
}

void check_conversion(Equation_of_State model, double temperature, double ion_fraction)
{
    const Cell cell = context();
    const Cell before = cell;
    Primitives input{};
    input[RHO] = 1e-10;
    double kinetic = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        input[VEL(d)] = (d % 2 == 0 ? 1 : -1) * (d + 1) * 1e5;
        kinetic += 0.5 * input[RHO] * input[VEL(d)] * input[VEL(d)];
    }
    if constexpr (NCOMP > 0) input[COMP(0)] = 1.0;
    if constexpr (NCOMP == 2)
    {
        input[COMP(0)] = 1.0 - ion_fraction;
        input[COMP(1)] = ion_fraction;
    }
    for (int k = 0; k < NUSER; ++k) input[CONUSR(k)] = -43.0 - k;
    // RhoP must ignore the stale auxiliary primitive energy, even when it is NaN.
    input[INTERNAL_ENERGY] = std::numeric_limits<double>::quiet_NaN();
    const auto input_before = input.data;

    const double mu = model == Equation_of_State::IdealGas_Default
        ? cell.mean_molecular_weight : 1.0 / (1.0 + ion_fraction);
    const double gas_pressure = input[RHO] / mu * physcon::R_g * temperature;
    double ion_energy = 0.0;
    if (model == Equation_of_State::Ionized || model == Equation_of_State::Ionized_Radiative)
        ion_energy = input[RHO] * ion_fraction / physcon::m_u
            * hd_microphysics::species_table.species_active[1].energy_eV * physcon::eV;
    const double radiation = model == Equation_of_State::Ionized_Radiative
        ? physcon::a_Rad * std::pow(temperature, 4) : 0.0;
    input[PRE] = gas_pressure + radiation / 3.0;
    const double expected_internal = gas_pressure / (cell.gamma - 1.0) + ion_energy + radiation;

    const auto saved_model = hd_microphysics::eos_type;
    hd_microphysics::eos_type = static_cast<Equation_of_State>(-1);
    // Explicit Default conversions stay ideal-gas even with an invalid global model.
    const auto default_primitive = input.to_Conserved_Default(cell.gamma)
        .to_Primitives_Default(cell.gamma);
    CHECK(near(default_primitive[PRE], input[PRE]));
    CHECK(near(default_primitive[INTERNAL_ENERGY],
               input[PRE] / (cell.gamma - 1.0) / input[RHO]));
    const auto converted = cell.to_Conserved(input, model);
    unchanged(cell, before);
    CHECK(converted[RHO] == input[RHO]);
    CHECK(near(converted[INTERNAL_ENERGY], expected_internal, 1e-9));
    CHECK(near(converted[ENERGY], expected_internal + kinetic, 1e-9));
    for (int d = 0; d < NDIM; ++d)
        CHECK(converted[MOMDEN(d)] == input[RHO] * input[VEL(d)]);
    for (int k = 0; k < NCOMP; ++k)
        CHECK(converted[COMP(k)] == input[RHO] * input[COMP(k)]);
    for (int k = 0; k < NUSER; ++k)
        CHECK(converted[CONUSR(k)] == input[CONUSR(k)]);

    // Reuse the converted state even with an invalid global EoS and stale primitive e.
    // The flux must not solve EoS again or derive thermal energy from pressure/gamma.
    const auto converted_before = converted.data;
    for (int d = 0; d < NDIM; ++d)
    {
        const auto flux = hydro::physical_flux(input, converted, d);
        const double velocity = input[VEL(d)];
        CHECK(near(flux[RHO], input[RHO] * velocity));
        CHECK(near(flux[INTERNAL_ENERGY], expected_internal * velocity, 1e-9));
        CHECK(near(flux[ENERGY], (expected_internal + kinetic + input[PRE]) * velocity, 1e-9));
        for (int j = 0; j < NDIM; ++j)
            CHECK(near(flux[MOMDEN(j)], input[RHO] * velocity * input[VEL(j)]
                       + (j == d ? input[PRE] : 0.0)));
        for (int k = 0; k < NCOMP; ++k)
            CHECK(near(flux[COMP(k)], input[RHO] * velocity * input[COMP(k)]));
        for (int k = 0; k < NUSER; ++k) CHECK(flux[CONUSR(k)] == 0.0);
    }
    CHECK(converted.data == converted_before);

    Cell evolved = cell;
    evolved.U = converted;
    const Cell evolved_before = evolved;
    const auto primitive = evolved.Get_Primitives(model);
    unchanged(evolved, evolved_before);
    CHECK(near(primitive[PRE], input[PRE], 1e-9));
    CHECK(near(primitive[INTERNAL_ENERGY], expected_internal / input[RHO], 1e-9));
    CHECK(primitive[RHO] == input[RHO]);
    for (int d = 0; d < NDIM; ++d)
        CHECK(near(primitive[VEL(d)], input[VEL(d)]));
    for (int k = 0; k < NCOMP; ++k)
        CHECK(near(primitive[COMP(k)], input[COMP(k)]));
    for (int k = 0; k < NUSER; ++k)
        CHECK(primitive[CONUSR(k)] == input[CONUSR(k)]);

    // Refresh the actual Cell once, then recover the same primitives from its cache.
    Cell refreshed = evolved;
    double refreshed_pressure = 0.0;
    hd_microphysics::evaluate_eos(refreshed, hd_microphysics::EoSInput::RhoE,
                                 model, refreshed_pressure);
    CHECK(refreshed.U.data == evolved.U.data);
    const Cell refreshed_before = refreshed;
    const auto fast_primitive = refreshed.Get_Primitives_From_Updated_State(model);
    CHECK(refreshed.Get_Pressure_From_Updated_State(model) == fast_primitive[PRE]);
    unchanged(refreshed, refreshed_before);
    CHECK(near(fast_primitive[PRE], refreshed_pressure, 1e-9));
    for (int k = 0; k < NVAR; ++k)
        CHECK(near(fast_primitive[k], primitive[k], 1e-9));

    Cell destination = cell;
    destination.Set_Conserved_From_Primitives(input, model);
    CHECK(destination.U.data == converted.data);
    CHECK(metadata(destination) == metadata(before));

    // Default arguments read the existing global selection at the time of the call.
    hd_microphysics::eos_type = model;
    auto face = input;
    face.update_internal_energy(cell.gamma);
    CHECK(near(face[INTERNAL_ENERGY], expected_internal / input[RHO], 1e-9));
    for (int k = 0; k < NVAR; ++k)
        if (k != INTERNAL_ENERGY) CHECK(face[k] == input[k]);
    for (int d = 0; d < NDIM; ++d)
        CHECK(hydro::physical_flux(face, d, cell.gamma).data
              == hydro::physical_flux(face, converted, d).data);
    CHECK(cell.to_Conserved(input).data == converted.data);
    CHECK(evolved.Get_Primitives().data == primitive.data);
    CHECK(refreshed.Get_Primitives_From_Updated_State().data == fast_primitive.data);
    CHECK(refreshed.Get_Pressure_From_Updated_State() == fast_primitive[PRE]);

    // A refreshed stage includes ghost storage and leaves every conserved field intact.
    MeshBlock stage;
    stage.nactive.fill(1);
    stage.nghost.fill(1);
    std::size_t count = 1;
    for (int d = 0; d < NDIM; ++d) count *= 3;
    stage.cells.assign(count, evolved);
    stage.cells.front() = Cell{};
    const auto vacuum_before = stage.cells.front();
    hd_microphysics::refresh_thermodynamics(stage);
    unchanged(stage.cells.front(), vacuum_before);
    for (std::size_t flat = 1; flat < stage.cells.size(); ++flat)
    {
        CHECK(stage.cells[flat].U.data == evolved.U.data);
        CHECK(near(stage.cells[flat].temperature, temperature, 1e-9));
        CHECK(near(stage.cells[flat].mean_molecular_weight, mu, 1e-9));
        CHECK(near(stage.cells[flat].Get_Pressure_From_Updated_State(), input[PRE], 1e-9));
    }
    Cell single = evolved;
    hd_microphysics::eos_type = static_cast<Equation_of_State>(-1);
    hd_microphysics::refresh_thermodynamics(single, model);
    unchanged(single, stage.cells.back());
    hd_microphysics::eos_type = model;

    destination = cell;
    destination.Set_Conserved_From_Primitives(input);
    CHECK(destination.U.data == converted.data);
    CHECK(metadata(destination) == metadata(before));
    hd_microphysics::eos_type = saved_model;

    CHECK(std::isnan(input[INTERNAL_ENERGY]));
    for (int k = 0; k < NVAR; ++k)
        if (k != PRE && k != INTERNAL_ENERGY) CHECK(input[k] == input_before[k]);
}

void check_energy_selection_and_failure()
{
    const auto model = Equation_of_State::IdealGas_Default;
    Cell cell = context();
    cell.U[RHO] = 1.0;
    cell.U[ENERGY] = 10.0;
    cell.U[INTERNAL_ENERGY] = 1.0;
    auto before = cell;
    auto primitive = cell.Get_Primitives(model);
    CHECK(near(primitive[PRE], 10.0 * (cell.gamma - 1.0)));
    CHECK(primitive[INTERNAL_ENERGY] == 10.0);
    unchanged(cell, before);
    double pressure = 0.0;
    hd_microphysics::evaluate_eos(cell, hd_microphysics::EoSInput::RhoE, model, pressure);
    before = cell;
    const auto fast_total = cell.Get_Primitives_From_Updated_State(model);
    CHECK(near(fast_total[PRE], primitive[PRE]));
    CHECK(fast_total[INTERNAL_ENERGY] == 10.0);
    unchanged(cell, before);

    cell.U[MOMDEN(0)] = 1e10;
    cell.U[ENERGY] = cell.U.kinetic_energy();
    cell.U[INTERNAL_ENERGY] = 3.0;
    before = cell;
    primitive = cell.Get_Primitives(model);
    CHECK(near(primitive[PRE], 3.0 * (cell.gamma - 1.0)));
    CHECK(primitive[INTERNAL_ENERGY] == 3.0);
    unchanged(cell, before);
    hd_microphysics::evaluate_eos(cell, hd_microphysics::EoSInput::RhoE, model, pressure);
    before = cell;
    const auto fast = cell.Get_Primitives_From_Updated_State(model);
    CHECK(near(fast[PRE], primitive[PRE]));
    CHECK(fast[INTERNAL_ENERGY] == 3.0);
    unchanged(cell, before);

    const auto invalid_model = static_cast<Equation_of_State>(-1);
    expect_throw<std::invalid_argument>([&] { cell.Get_Primitives(invalid_model); });
    expect_throw<std::invalid_argument>([&] { cell.Get_Primitives_From_Updated_State(invalid_model); });
    expect_throw<std::invalid_argument>([&] { cell.Get_Pressure_From_Updated_State(invalid_model); });
    expect_throw<std::invalid_argument>([&] { cell.to_Conserved(primitive, invalid_model); });
    expect_throw<std::invalid_argument>([&] { cell.Set_Conserved_From_Primitives(primitive, invalid_model); });
    unchanged(cell, before);

    for (double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()})
    {
        Cell invalid_cache = cell;
        invalid_cache.temperature = invalid;
        expect_throw<std::invalid_argument>([&] {
            invalid_cache.Get_Primitives_From_Updated_State(model);
        });
        expect_throw<std::invalid_argument>([&] {
            invalid_cache.Get_Pressure_From_Updated_State(model);
        });
        invalid_cache = cell;
        invalid_cache.mean_molecular_weight = invalid;
        expect_throw<std::invalid_argument>([&] {
            invalid_cache.Get_Primitives_From_Updated_State(model);
        });
        expect_throw<std::invalid_argument>([&] {
            invalid_cache.Get_Pressure_From_Updated_State(model);
        });
        if (invalid != 0.0)
        {
            invalid_cache = cell;
            invalid_cache.U[RHO] = invalid;
            expect_throw<std::invalid_argument>([&] {
                invalid_cache.Get_Pressure_From_Updated_State(model);
            });
        }
    }
    Cell overflow = cell;
    overflow.temperature = std::numeric_limits<double>::max();
    expect_throw<std::invalid_argument>([&] {
        overflow.Get_Pressure_From_Updated_State(model);
    });

    primitive[PRE] = -1.0;
    const auto invalid_face = primitive.data;
    expect_throw<std::invalid_argument>([&] { primitive.update_internal_energy(cell.gamma); });
    CHECK(primitive.data == invalid_face);
    expect_throw<std::invalid_argument>([&] { cell.Set_Conserved_From_Primitives(primitive, model); });
    unchanged(cell, before);

    if constexpr (NCOMP == 2)
    {
        // Density caches are updated on the scratch Cell before Newton rejects T0.
        primitive = Primitives{};
        primitive[RHO] = 1e-10;
        primitive[PRE] = std::numeric_limits<double>::max();
        primitive[COMP(0)] = 1.0;
        expect_throw<std::invalid_argument>([&] {
            cell.Set_Conserved_From_Primitives(primitive, Equation_of_State::Ionized_Radiative);
        });
        unchanged(cell, before);
    }

    Cell vacuum = context();
    vacuum.U = Conserved{};
    const auto vacuum_before = vacuum;
    hd_microphysics::refresh_thermodynamics(vacuum, model);
    unchanged(vacuum, vacuum_before);
    CHECK(vacuum.Get_Primitives(model).data == Primitives{}.data);
    CHECK(vacuum.Get_Primitives_From_Updated_State(model).data == Primitives{}.data);
    CHECK(vacuum.Get_Pressure_From_Updated_State(model) == 0.0);
    unchanged(vacuum, vacuum_before);
    // Vacuum needs no thermodynamic cache; extra conserved fields survive.
    vacuum.temperature = vacuum.mean_molecular_weight = 0.0;
    CHECK(vacuum.Get_Pressure_From_Updated_State(model) == 0.0);
    for (int k = 0; k < NUSER; ++k) vacuum.U[CONUSR(k)] = 53.0 + k;
    const auto empty_primitive = vacuum.Get_Primitives_From_Updated_State(model);
    CHECK(empty_primitive[RHO] == 0.0);
    CHECK(empty_primitive[PRE] == 0.0);
    for (int k = 0; k < NUSER; ++k)
        CHECK(empty_primitive[CONUSR(k)] == vacuum.U[CONUSR(k)]);
    Primitives empty{};
    empty[INTERNAL_ENERGY] = std::numeric_limits<double>::quiet_NaN();
    empty.update_internal_energy(cell.gamma);
    CHECK(empty.data == Primitives{}.data);
    CHECK(hydro::physical_flux(empty, Conserved{}, 0).data == Conserved{}.data);
    expect_throw<std::out_of_range>([&] {
        hydro::physical_flux(empty, Conserved{}, NDIM);
    });
    CHECK(cell.to_Conserved(empty, model).data == Conserved{}.data);
    unchanged(cell, before);
    cell.Set_Conserved_From_Primitives(empty, model);
    CHECK(cell.U.data == Conserved{}.data);
    CHECK(metadata(cell) == metadata(before));
}
}

int main(int argc, char** argv)
{
    return run_test("Cell EoS conversion isolation and U-only setter", [&] {
        if constexpr (hd_core::NCOMP == 2)
        {
            CHECK(argc == 2);
            const auto directory = std::filesystem::absolute(argv[1]);
            std::filesystem::create_directories(directory);
            const auto active_file = directory / "active.txt";
            std::ofstream(active_file) << "h1 0 1\n";
            config::SimulationConfig cfg;
            cfg.composition_activate_file = active_file.string();
            hd_microphysics::initialize_species_table(cfg);
            CHECK(hd_microphysics::species_table.species_active[0].charge == 0);
            CHECK(hd_microphysics::species_table.species_active[1].charge == 1);
        }
        else
            hd_microphysics::initialize_species_table(config::SimulationConfig{});
        check_constructor_density();
        for (double temperature : {100.0, 1e4, 1e7})
        {
            check_conversion(Equation_of_State::IdealGas_Default, temperature, 0.0);
            if constexpr (hd_core::NCOMP == 2)
                for (double ion_fraction : {0.0, 0.35, 1.0})
                    for (auto model : {Equation_of_State::IdealGas_Mixed,
                         Equation_of_State::Ionized, Equation_of_State::Ionized_Radiative})
                        check_conversion(model, temperature, ion_fraction);
        }
        check_energy_selection_and_failure();
    });
}
