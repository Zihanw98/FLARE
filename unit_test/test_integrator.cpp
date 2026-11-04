#include "test_support.h"
#include "time/integrator.h"
#include <limits>

int main()
{
    return run_test("integrator active update and rollback", [] {
        hd_core::MeshBlock mesh{};
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            mesh.nactive[d] = 3 + d;
            mesh.nghost[d] = 1 + d;
            mesh.set_uniform_faces(d, 0.0, 1.0);
            mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
        }
        mesh.allocate();
        std::vector<bool> active(mesh.cells.size(), true);
        std::vector<hd_core::Conserved> rhs(mesh.cells.size());
        std::size_t last_active = 0;
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            mesh.cells[flat].U[hd_core::RHO] = 1.0;
            mesh.cells[flat].U[hd_core::ENERGY] = 10.0;
            mesh.cells[flat].gamma = flat % 2 ? 1.4 : 5.0 / 3.0;
            const auto index = mesh.coordinates(flat);
            for (int d = 0; d < hd_core::NDIM; ++d)
                active[flat] = active[flat] && index[d] >= mesh.nghost[d]
                    && index[d] < mesh.nghost[d] + mesh.nactive[d];
            rhs[flat].data.fill(active[flat] ? 0.25
                : std::numeric_limits<double>::quiet_NaN());
            if (active[flat]) last_active = flat;
        }
        const auto original = mesh.cells;
        auto trial = mesh;
        hd_time::add_scaled_rhs(trial, rhs, 0.5);
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            CHECK(mesh.cells[flat].U.data == original[flat].U.data);
            for (int k = 0; k < hd_core::NVAR; ++k)
                CHECK(trial.cells[flat].U[k] == original[flat].U[k]
                      + (active[flat] ? 0.125 : 0.0));
            CHECK(trial.cells[flat].gamma == original[flat].gamma);
        }
        CHECK(trial.face == mesh.face);
        CHECK(trial.bc_lower == mesh.bc_lower);
        CHECK(trial.bc_upper == mesh.bc_upper);

        // Fail at the last active cell: preceding candidate updates must not leak.
        auto bad_rhs = rhs;
        bad_rhs[last_active][hd_core::ENERGY] = -1e6;
        auto check_unchanged = [&] {
            for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
            {
                CHECK(mesh.cells[flat].U.data == original[flat].U.data);
                CHECK(mesh.cells[flat].gamma == original[flat].gamma);
            }
        };
        expect_throw<std::invalid_argument>([&] {
            hd_time::add_scaled_rhs(mesh, bad_rhs, 0.5);
        });
        check_unchanged();
        bad_rhs[last_active][hd_core::ENERGY] = std::numeric_limits<double>::quiet_NaN();
        expect_throw<std::invalid_argument>([&] {
            hd_time::add_scaled_rhs(mesh, bad_rhs, 0.5);
        });
        check_unchanged();
        expect_throw<std::invalid_argument>([&] { hd_time::add_scaled_rhs(mesh, {}, 0.5); });
        check_unchanged();
        expect_throw<std::invalid_argument>([&] {
            hd_time::add_scaled_rhs(mesh, rhs, std::numeric_limits<double>::infinity());
        });
        check_unchanged();
        hd_time::add_scaled_rhs(mesh, bad_rhs, 0.0); // Zero scale is a no-op.
        check_unchanged();
        hd_time::add_scaled_rhs(trial, rhs, -0.5);
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
            CHECK(trial.cells[flat].U.data == original[flat].U.data);
        mesh.cells.clear();
        expect_throw<std::logic_error>([&] { hd_time::add_scaled_rhs(mesh, {}, 0.5); });
    });
}
