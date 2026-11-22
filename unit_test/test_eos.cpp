#include "test_support.h"
#include "core/cell.h"
#include "config/config.h"
#include "microphysics/composition/species.h"
#include "microphysics/eos/eos.h"
#include "physcon/constants.h"
#include <filesystem>
#include <fstream>
#include <limits>

namespace
{
using namespace hd_core;
using namespace hd_microphysics;

Cell state()
{
    Cell cell;
    cell.U[RHO] = 1e-10;
    cell.U[MOMDEN(0)] = 2e-5;
    cell.mean_molecular_weight = 0.7;
    if constexpr (NCOMP == 2)
    {
        cell.U[COMP(0)] = 0.65 * cell.U[RHO];
        cell.U[COMP(1)] = 0.35 * cell.U[RHO];
    }
    return cell;
}

void roundtrip(Equation_of_State model)
{
    for (double gamma : {1.2, 4.0/3.0, 5.0/3.0})
    for (double temperature : {100.0, 1e4, 1e7})
    {
        auto cell = state();
        cell.gamma = gamma;
        cell.temperature = temperature;
        const double mu = model == Equation_of_State::IdealGas_Default ? 0.7 : 1.0/1.35;
        const double pg = cell.U[RHO] * physcon::R_g * temperature / mu;
        const double rad = model == Equation_of_State::Ionized_Radiative
            ? physcon::a_Rad * std::pow(temperature, 4) : 0.0;
        const double ion = model == Equation_of_State::Ionized || model == Equation_of_State::Ionized_Radiative
            ? cell.U[COMP(1)] / physcon::m_u * species_table.species_active[1].energy_eV * physcon::eV : 0.0;
        const double expected_energy = pg/(gamma-1) + ion + rad;
        double pressure = -1;
        evaluate_eos(cell, EoSInput::RhoT, model, pressure);
        CHECK(near(pressure, pg + rad/3, 2e-10));
        CHECK(near(cell.U[INTERNAL_ENERGY], expected_energy, 2e-10));
        CHECK(near(cell.U[ENERGY], expected_energy + cell.U.kinetic_energy(), 2e-10));
        const auto initial = cell.U.data;
        cell.temperature = std::numeric_limits<double>::quiet_NaN();
        evaluate_eos(cell, EoSInput::RhoE, model, pressure);
        CHECK(cell.U.data == initial);
        CHECK(near(cell.temperature, temperature, 2e-10));
        CHECK(near(cell.Get_Pressure_From_Updated_State(model), pg + rad/3, 2e-10));
        const double input_pressure = pressure;
        evaluate_eos(cell, EoSInput::RhoP, model, pressure);
        CHECK(pressure == input_pressure);
        CHECK(near(cell.temperature, temperature, 2e-10));
        CHECK(near(cell.U[INTERNAL_ENERGY], expected_energy, 2e-10));

        // The explicit model must win over a conflicting global selection.
        eos_type = Equation_of_State::Ionized_Radiative;
        const auto before = cell.U.data;
        const double sound = cell.get_sound_speed(model);
        const double A = pg + 4*rad/3;
        const double B = pg/(gamma-1) + 4*rad;
        const double expected_c2 = model == Equation_of_State::Ionized_Radiative
            ? (pg + A*A/B)/cell.U[RHO] : gamma*pg/cell.U[RHO];
        CHECK(near(sound*sound, expected_c2, 3e-10));
        CHECK(cell.U.data == before);
    }
}

void failures(Equation_of_State model)
{
    double pressure = 1;
    auto cell = state();
    cell.temperature = 1e4;
    evaluate_eos(cell, EoSInput::RhoT, model, pressure);
    const auto valid = cell;
    for (double value : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()})
    {
        cell = valid; cell.U[RHO] = value;
        expect_throw<std::invalid_argument>([&] { evaluate_eos(cell, EoSInput::RhoE, model, pressure); });
        cell = valid; cell.temperature = value;
        expect_throw<std::invalid_argument>([&] { evaluate_eos(cell, EoSInput::RhoT, model, pressure); });
        expect_throw<std::invalid_argument>([&] { cell.get_sound_speed(model); });
        cell = valid; pressure = value;
        expect_throw<std::invalid_argument>([&] { evaluate_eos(cell, EoSInput::RhoP, model, pressure); });
        cell = valid; cell.mean_molecular_weight = value;
        expect_throw<std::invalid_argument>([&] { cell.get_sound_speed(model); });
        if (model == Equation_of_State::IdealGas_Default)
            expect_throw<std::invalid_argument>([&] { evaluate_eos(cell, EoSInput::RhoT, model, pressure); });
    }
    for (double gamma : {1.0, 0.5, std::numeric_limits<double>::quiet_NaN()})
    {
        cell = valid; cell.gamma = gamma;
        expect_throw<std::invalid_argument>([&] { evaluate_eos(cell, EoSInput::RhoE, model, pressure); });
        expect_throw<std::invalid_argument>([&] { cell.get_sound_speed(model); });
    }
    if (model == Equation_of_State::Ionized || model == Equation_of_State::Ionized_Radiative)
    {
        cell = valid;
        cell.U[MOMDEN(0)] = 0;
        const double ion = cell.U[COMP(1)] / physcon::m_u * species_table.species_active[1].energy_eV * physcon::eV;
        for (double fraction : {0.5, 1.0})
        {
            cell.U[ENERGY] = cell.U[INTERNAL_ENERGY] = fraction*ion;
            expect_throw<std::invalid_argument>([&] { refresh_thermodynamics(cell, model); });
        }
    }
    CHECK(Cell{}.get_sound_speed(model) == 0.0);
}
}

int main(int argc, char** argv)
{
    return run_test("EoS roundtrips, explicit sound-speed model and admissibility", [&] {
        CHECK(NCOMP == 0 || NCOMP == 2);
        config::SimulationConfig cfg;
        if constexpr (NCOMP == 2)
        {
            CHECK(argc == 2);
            const auto directory = std::filesystem::absolute(argv[1]);
            std::filesystem::create_directories(directory);
            cfg.composition_activate_file = (directory / "active.txt").string();
            std::ofstream(cfg.composition_activate_file) << "h1 0 1\n";
        }
        initialize_species_table(cfg);
        CHECK(parse_equation_of_state("idealgas_default") == Equation_of_State::IdealGas_Default);
        CHECK(parse_equation_of_state("IdealGas_Mixed") == Equation_of_State::IdealGas_Mixed);
        CHECK(parse_equation_of_state("Ionized") == Equation_of_State::Ionized);
        CHECK(parse_equation_of_state("iOnIzEd_rAdIaTiVe") == Equation_of_State::Ionized_Radiative);
        for (const char* name : {"", "IdealGas", "Ionized_Unknown"})
            expect_throw<std::invalid_argument>([&] { parse_equation_of_state(name); });
        for (auto model : {Equation_of_State::IdealGas_Default, Equation_of_State::IdealGas_Mixed,
                           Equation_of_State::Ionized, Equation_of_State::Ionized_Radiative})
        {
            if (NCOMP == 0 && model != Equation_of_State::IdealGas_Default) continue;
            roundtrip(model);
            failures(model);
        }
        auto cell = state();
        double pressure = 1;
        const auto before = cell.U.data;
        expect_throw<std::invalid_argument>([&] {
            evaluate_eos(cell, EoSInput::RhoT, static_cast<Equation_of_State>(-1), pressure);
        });
        expect_throw<std::invalid_argument>([&] {
            evaluate_eos(cell, static_cast<EoSInput>(-1), Equation_of_State::IdealGas_Default, pressure);
        });
        CHECK(cell.U.data == before && pressure == 1);
    });
}
