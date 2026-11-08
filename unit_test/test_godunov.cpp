#include "test_support.h"
#include "core/mesh.h"
#include "hydro/riemann/riemann.h"
#include "hydro/spatial.h"
#include "hydro/wave_timestep.h"
#include <vector>

namespace
{
struct FaceCall { double left, right, gamma_left, gamma_right; int direction; };

// A prescribed nonconstant flux gives an independently known divergence.
struct ProbeSolver final : hydro::RiemannSolver
{
    mutable std::vector<FaceCall> calls;
    hd_core::Conserved flux(const hd_core::Primitives& left, const hd_core::Primitives& right,
                            int d, double gl, double gr) const override
    {
        const double l = left[hd_core::RHO], r = right[hd_core::RHO];
        calls.push_back({l, r, gl, gr, d});
        hydro::wave_timestep.face_speed = 1.0 + d + 0.1 * l;
        hd_core::Conserved result{};
        for (int k = 0; k < hd_core::NVAR; ++k)
            result[k] = (k + 1.0) * (d + 1.0) * (l + 2.0 * r);
        return result;
    }
};
}

int main()
{
    return run_test("Godunov face neighbors, divergence sign and CFL reduction", [] {
        hd_core::MeshBlock mesh{};
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            mesh.nactive[d] = 2 + d;
            mesh.nghost[d] = 2;
            mesh.face[d].push_back(-1.0);
            for (int i = 0; i < mesh.ntotal(d); ++i)
                mesh.face[d].push_back(mesh.face[d].back() + 0.2 * (d + 1) * (1 + 0.05 * i));
        }
        mesh.allocate();
        mesh.build_geometry();
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            hd_core::Primitives w{};
            w[hd_core::RHO] = 2.0 + 0.01 * flat;
            w[hd_core::PRE] = 2.0;
            mesh.cells[flat] = hd_core::Cell(w, 1.4 + 0.0001 * flat);
        }
        const auto original = mesh.cells;
        ProbeSolver solver;
        std::vector<hd_core::Conserved> rhs;
        hydro::compute_hydro_rhs(mesh, hydro::SpatialScheme::Godunov, solver, rhs);
        CHECK(rhs.size() == mesh.cells.size());

        // Independently enumerate face-local coordinates and translate to storage.
        std::size_t call_number = 0;
        for (int normal = 0; normal < hd_core::NDIM; ++normal)
        {
            auto shape = mesh.nactive;
            ++shape[normal];
            std::size_t face_count = 1;
            for (int extent : shape) face_count *= extent;
            for (std::size_t face = 0; face < face_count; ++face)
            {
                auto right = mesh.nghost;
                auto remainder = face;
                for (int d = 0; d < hd_core::NDIM; ++d)
                {
                    right[d] += remainder % shape[d];
                    remainder /= shape[d];
                }
                auto left = right;
                --left[normal];
                CHECK(call_number < solver.calls.size());
                const auto& call = solver.calls[call_number++];
                CHECK(call.direction == normal);
                CHECK(call.left == mesh.at(left).U[hd_core::RHO]);
                CHECK(call.right == mesh.at(right).U[hd_core::RHO]);
                CHECK(call.gamma_left == mesh.at(left).gamma);
                CHECK(call.gamma_right == mesh.at(right).gamma);
            }
        }
        CHECK(call_number == solver.calls.size());

        double expected_max_rate = 0.0;
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            const auto index = mesh.coordinates(flat);
            bool active = true;
            for (int d = 0; d < hd_core::NDIM; ++d)
                active &= index[d] >= mesh.nghost[d] && index[d] < mesh.nghost[d] + mesh.nactive[d];
            double mass_rhs = 0.0, rate = 0.0;
            std::size_t stride = 1;
            if (active)
            {
                for (int d = 0; d < hd_core::NDIM; ++d)
                {
                    const double width = mesh.face[d][index[d] + 1] - mesh.face[d][index[d]];
                    // rho rises by 0.01 * stride between neighboring cells.
                    // F_right - F_left = (d+1) * 3 * 0.01 * stride > 0.
                    mass_rhs -= (d + 1.0) * 0.03 * stride / width;
                    rate += (1.0 + d + 0.1 * (2.0 + 0.01 * flat)) / width;
                    stride *= mesh.nactive[d] + 2 * mesh.nghost[d];
                }
                expected_max_rate = std::max(expected_max_rate, rate);
            }
            for (int k = 0; k < hd_core::NVAR; ++k)
                CHECK(near(rhs[flat][k], (k + 1.0) * mass_rhs));
            CHECK(mesh.cells[flat].U.data == original[flat].U.data);
            CHECK(mesh.cells[flat].gamma == original[flat].gamma);
        }
        CHECK(near(hydro::wave_timestep.max_inverse_dt, expected_max_rate));
    });
}
