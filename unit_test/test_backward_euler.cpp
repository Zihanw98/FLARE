#include "test_support.h"
#include "user/user_api.h"
#include "hydro/riemann/exact.h"
#include "hydro/wave_timestep.h"
// Exercise the private Jacobian without adding an unfinished solver to the API.
// The archive's backward_euler object is not pulled into this executable.
#include "../src/time/backward_euler.cpp"

namespace
{
class ZeroFluxSolver : public hydro::RiemannSolver
{
public:
    hd_core::Conserved flux(const hd_core::Primitives&, const hd_core::Primitives&,
                           int, double, double) const override
    {
        hydro::wave_timestep.face_speed = 0.0;
        return {};
    }
};

void jacobian_tests()
{
    hd_core::MeshBlock mesh;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        mesh.nactive[d] = d == 0 ? 3 : 1;
        mesh.nghost[d] = 1;
        mesh.set_linear_uniform_faces(d, 0.0, 1.0);
        mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
    }
    mesh.allocate();
    mesh.build_geometry();
    std::vector<std::size_t> active;
    for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
    {
        auto& cell = mesh.cells[flat];
        cell.U[hd_core::RHO] = 1.0;
        cell.U[hd_core::ENERGY] = 4.0;
        if (mesh.is_active(flat))
        {
            active.push_back(flat);
            cell.U[hd_core::RHO] += 0.1 * active.size();
            cell.U[hd_core::MOMDEN(0)] = 0.1 * active.size();
        }
    }
    std::reverse(active.begin(), active.end()); // Ordering need not follow storage.
    const auto old = mesh;
    const config::SimulationConfig cfg;
    const hd_user::UserHooks hooks;
    const ZeroFluxSolver zero_solver;
    constexpr auto scheme = hydro::SpatialScheme::Godunov;
    constexpr double time = 0.75;
    constexpr double cfl = 0.4;
    int source_calls = 0;
    const hd_source::SourceList sources{{"linear periodic coupling", hd_source::SourceKind::General,
        [&](const hd_core::MeshBlock& state, double t, const config::SimulationConfig&,
            std::vector<hd_core::Conserved>& rhs, double&) {
            CHECK(t == time);
            ++source_calls;
            for (auto flat : active)
            {
                auto left = state.coordinates(flat);
                --left[0]; // At the lower boundary this reads a periodic ghost.
                for (int k = 0; k < hd_core::NVAR; ++k)
                    rhs[flat][k] += -2.0 * state.cells[flat].U[k]
                        + t * state.at(left).U[(k + 1) % hd_core::NVAR];
            }
        }}};
    const auto check_unchanged = [&](const hd_core::MeshBlock& before) {
        CHECK(mesh.face == before.face);
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            CHECK(mesh.cells[flat].U.data == before.cells[flat].U.data);
            CHECK(mesh.cells[flat].gamma == before.cells[flat].gamma);
        }
    };

    for (double dt : {0.0, 0.2})
    {
        const auto R = hd_time::Residual(mesh, old, time, dt, scheme,
                                         sources, zero_solver, cfg, hooks, cfl);
        const auto before = mesh;
        source_calls = 0;
        const auto J = hd_time::finite_difference_jacobian(mesh, old, active, R,
            time, dt, scheme, sources, zero_solver, cfg, hooks, cfl);
        const auto n = active.size() * hd_core::NVAR;
        CHECK(J.size() == n);
        CHECK(source_calls == static_cast<int>(n));
        for (std::size_t i = 0; i < active.size(); ++i)
        {
            auto left = mesh.coordinates_active(active[i]);
            left[0] = (left[0] + mesh.nactive[0] - 1) % mesh.nactive[0];
            const auto left_flat = mesh.flat_index(left, true);
            for (int k = 0; k < hd_core::NVAR; ++k)
                for (std::size_t column = 0; column < n; ++column)
                {
                    const auto row = i * hd_core::NVAR + k;
                    double expected = row == column ? 1.0 + 2.0 * dt : 0.0;
                    if (active[column / hd_core::NVAR] == left_flat
                        && static_cast<int>(column % hd_core::NVAR) == (k + 1) % hd_core::NVAR)
                        expected -= dt * time;
                    CHECK(near(J[row][column], expected, 2e-7));
                }
        }
        check_unchanged(before);

        auto duplicate = active;
        duplicate.back() = duplicate.front();
        expect_throw<std::invalid_argument>([&] {
            hd_time::finite_difference_jacobian(mesh, old, duplicate, R,
                time, dt, scheme, sources, zero_solver, cfg, hooks, cfl);
        });
        const hd_source::SourceList failing{{"failure", hd_source::SourceKind::General,
            [](const hd_core::MeshBlock&, double, const config::SimulationConfig&,
               std::vector<hd_core::Conserved>&, double&) {
                throw std::runtime_error("source failure");
            }}};
        expect_throw<std::runtime_error>([&] {
            hd_time::finite_difference_jacobian(mesh, old, active, R,
                time, dt, scheme, failing, zero_solver, cfg, hooks, cfl);
        });
        check_unchanged(before);
    }

    // Independent central directional difference of the real hydrodynamic RHS.
    const hydro::ExactSolver exact_solver;
    constexpr double dt = 0.03;
    const auto R = hd_time::Residual(mesh, old, time, dt, scheme,
                                     {}, exact_solver, cfg, hooks, cfl);
    const auto J = hd_time::finite_difference_jacobian(mesh, old, active, R,
        time, dt, scheme, {}, exact_solver, cfg, hooks, cfl);
    const auto n = J.size();
    std::vector<double> direction(n);
    auto plus = mesh;
    auto minus = mesh;
    constexpr double h = 1e-5;
    for (std::size_t j = 0; j < n; ++j)
    {
        direction[j] = std::sin(static_cast<double>(j + 1));
        plus.cells[active[j / hd_core::NVAR]].U[j % hd_core::NVAR] += h * direction[j];
        minus.cells[active[j / hd_core::NVAR]].U[j % hd_core::NVAR] -= h * direction[j];
    }
    const auto Rp = hd_time::Residual(plus, old, time, dt, scheme,
                                      {}, exact_solver, cfg, hooks, cfl);
    const auto Rm = hd_time::Residual(minus, old, time, dt, scheme,
                                      {}, exact_solver, cfg, hooks, cfl);
    for (std::size_t i = 0; i < n; ++i)
    {
        double Jv = 0.0;
        for (std::size_t j = 0; j < n; ++j) Jv += J[i][j] * direction[j];
        const auto flat = active[i / hd_core::NVAR];
        const int k = static_cast<int>(i % hd_core::NVAR);
        CHECK(near(Jv, (Rp[flat][k] - Rm[flat][k]) / (2.0 * h), 2e-5));
    }

    // Positive momentum probes would make this cold moving gas unphysical.
    for (auto flat : active)
    {
        mesh.cells[flat].U = {};
        mesh.cells[flat].U[hd_core::RHO] = 1.0;
        mesh.cells[flat].U[hd_core::MOMDEN(0)] = 1.0;
        mesh.cells[flat].U[hd_core::ENERGY] = 0.5 + 1e-12;
    }
    const auto cold_R = hd_time::Residual(mesh, old, time, 0.0, scheme,
                                          {}, zero_solver, cfg, hooks, cfl);
    const auto cold_J = hd_time::finite_difference_jacobian(mesh, old, active, cold_R,
        time, 0.0, scheme, {}, zero_solver, cfg, hooks, cfl);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            CHECK(near(cold_J[i][j], i == j ? 1.0 : 0.0, 2e-7));
}
}

int main()
{
    return run_test("backward Euler residual Jacobian", jacobian_tests);
}
