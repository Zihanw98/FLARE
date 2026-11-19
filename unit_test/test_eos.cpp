#include "test_support.h"
#include "core/cell.h"
#include "microphysics/eos/eos.h"
#include "physcon/constants.h"
#include <array>

int main()
{
    using hd_microphysics::Equation_of_State;
    using hd_microphysics::EoSInput;
    using hd_microphysics::evaluate_eos;
    using hd_microphysics::parse_equation_of_state;
    return run_test("EoS selection and unified evaluation", [] {
        struct Selection
        {
            const char* name;
            Equation_of_State type;
        };
        const std::array<Selection, 4> selections{{
            {"IdealGas_Default", Equation_of_State::IdealGas_Default},
            {"IdealGas_Mixed", Equation_of_State::IdealGas_Mixed},
            {"Ionized", Equation_of_State::Ionized},
            {"Ionized_Radiative", Equation_of_State::Ionized_Radiative},
        }};
        hd_core::Cell cell{};
        cell.U[hd_core::RHO] = 2.0;
        cell.U[hd_core::MOMDEN(0)] = 6.0;
        cell.mean_molecular_weight = 0.5;
        cell.temperature = 1.0 / (4.0 * physcon::R_g);
        CHECK(hd_microphysics::eos_type == Equation_of_State::IdealGas_Default);
        for (const auto& selection : selections)
        {
            hd_microphysics::eos_type = parse_equation_of_state(selection.name);
            CHECK(hd_microphysics::eos_type == selection.type);
            for (const auto input : {EoSInput::RhoP, EoSInput::RhoE})
            {
                const auto initial_state = cell.U.data;
                const auto initial_temperature = cell.temperature;
                double pressure = 3.0;
                bool unavailable = false;
                try { evaluate_eos(cell, input, hd_microphysics::eos_type, pressure); }
                catch (const std::logic_error& error)
                {
                    unavailable = true;
                    CHECK(std::string(error.what()) ==
                          (input == EoSInput::RhoP ? "EoS RhoP is not implemented"
                                                  : "EoS RhoE is not implemented"));
                }
                CHECK(unavailable);
                CHECK(cell.U.data == initial_state);
                CHECK(cell.temperature == initial_temperature && pressure == 3.0);
            }
        }

        // The explicit model argument controls evaluation, independently of the global selection.
        double pressure = -1.0;
        evaluate_eos(cell, EoSInput::RhoT, Equation_of_State::IdealGas_Default, pressure);
        CHECK(near(pressure, 1.0));
        CHECK(near(cell.U[hd_core::INTERNAL_ENERGY], 1.5));
        CHECK(near(cell.U[hd_core::ENERGY], 10.5));

        hd_microphysics::eos_type = parse_equation_of_state("IdealGas_Default");
        pressure = -1.0;
        evaluate_eos(cell, EoSInput::RhoT, hd_microphysics::eos_type, pressure);
        CHECK(near(pressure, 1.0));
        CHECK(parse_equation_of_state("idealgas_default") == Equation_of_State::IdealGas_Default);
        CHECK(parse_equation_of_state("iOnIzEd_rAdIaTiVe") == Equation_of_State::Ionized_Radiative);
        for (const char* name : {"", "IdealGas", "Ionized_Unknown"})
            expect_throw<std::invalid_argument>([&] { parse_equation_of_state(name); });
        const auto initial_state = cell.U.data;
        expect_throw<std::invalid_argument>([&] {
            evaluate_eos(cell, EoSInput::RhoT, static_cast<Equation_of_State>(-1), pressure);
        });
        expect_throw<std::invalid_argument>([&] {
            evaluate_eos(cell, static_cast<EoSInput>(-1), hd_microphysics::eos_type, pressure);
        });
        CHECK(cell.U.data == initial_state && near(pressure, 1.0));
    });
}
