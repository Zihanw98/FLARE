#include "test_support.h"
#include "core/cell.h"
#include "microphysics/composition/species.h"
#include "microphysics/eos/eos.h"
#include <array>

int main()
{
    using hd_microphysics::Equation_of_State;
    using hd_microphysics::get_eos_function;
    using hd_microphysics::parse_equation_of_state;
    return run_test("EoS parsing and function dispatch", [] {
        struct Selection
        {
            const char* name;
            Equation_of_State type;
        };
        const std::array<Selection, 5> selections{{
            {"IdealGas_Default", Equation_of_State::IdealGas_Default},
            {"IdealGas_Single", Equation_of_State::IdealGas_Single},
            {"IdealGas_Mixed", Equation_of_State::IdealGas_Mixed},
            {"Ionized", Equation_of_State::Ionized},
            {"Ionized_Radiative", Equation_of_State::Ionized_Radiative},
        }};
        hd_core::Cell cell{};
        const hd_microphysics::SpeciesTable species{};
        const auto initial_state = cell.U.data;
        const auto initial_gamma = cell.gamma;
        const auto initial_temperature = cell.temperature;
        for (const auto& selection : selections)
        {
            const auto type = parse_equation_of_state(selection.name);
            CHECK(type == selection.type);
            const auto eos = get_eos_function(type);
            CHECK(eos != nullptr);

            // 调用公共指针必须到达选中的占位实现，并明确拒绝未实现的物理计算。
            bool unavailable = false;
            try { eos(cell, species); }
            catch (const std::logic_error& error)
            {
                unavailable = true;
                CHECK(std::string(error.what()) ==
                      std::string("EoS ") + selection.name + " is not implemented");
            }
            CHECK(unavailable);
            CHECK(cell.U.data == initial_state);
            CHECK(cell.gamma == initial_gamma && cell.temperature == initial_temperature);
        }
        CHECK(parse_equation_of_state("idealgas_default") == Equation_of_State::IdealGas_Default);
        CHECK(parse_equation_of_state("iOnIzEd_rAdIaTiVe") == Equation_of_State::Ionized_Radiative);
        for (const char* name : {"", "IdealGas", "Ionized_Unknown"})
            expect_throw<std::invalid_argument>([&] { parse_equation_of_state(name); });
        expect_throw<std::invalid_argument>([] {
            get_eos_function(static_cast<Equation_of_State>(-1));
        });
    });
}
