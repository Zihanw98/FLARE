#include "test_support.h"
#include "core/mesh.h"
#include "config/config.h"
#include "hydro/facefield.h"
#include <limits>

int main()
{
    return run_test("mesh indexing, geometry and active-state validation", [] {
        hd_core::MeshBlock mesh{};
        std::size_t count = 1;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            mesh.nactive[d] = 3 + d;
            mesh.nghost[d] = 1 + d;
            mesh.set_linear_uniform_faces(d, -1.0, 2.0 + d);
            count *= mesh.nactive[d] + 2 * mesh.nghost[d];
            CHECK(near(mesh.cell_center(d, mesh.nghost[d]), -0.5));
            CHECK(near(mesh.cell_width(d, 0), 1.0));
            CHECK(near(mesh.active_length(d), 3.0 + d));
        }
        mesh.allocate();
        CHECK(mesh.cells.size() == count);
        const auto& const_mesh = mesh;
        for (std::size_t flat = 0; flat < count; ++flat)
        {
            const auto storage = mesh.coordinates(flat);
            const auto local = mesh.coordinates_active(flat);
            CHECK(mesh.flat_index(storage) == flat);
            CHECK(&const_mesh.at(storage) == &mesh.cells[flat]);
            CHECK(&mesh.at(storage) == &mesh.cells[flat]);
            bool active = true;
            std::size_t expected = 0, stride = 1;
            for (int d = 0; d < hd_core::NDIM; ++d)
            {
                expected += storage[d] * stride;
                stride *= mesh.nactive[d] + 2 * mesh.nghost[d];
                CHECK(local[d] == storage[d] - mesh.nghost[d]);
                CHECK(near(mesh.cell_width(d, storage), 1.0));
                active &= local[d] >= 0 && local[d] < mesh.nactive[d];
            }
            CHECK(expected == flat); // Axis 0 has unit stride.
            if (active)
            {
                CHECK(mesh.flat_index(local, true) == flat);
                CHECK(&mesh.at_active(local) == &mesh.cells[flat]);
                CHECK(&const_mesh.at_active(local) == &mesh.cells[flat]);
            }
        }
        // Along an axis widths may vary; the other coordinates do not enter.
        mesh.face[0][1] = mesh.face[0][0] + 0.25;
        CHECK(near(mesh.cell_width(0, 0), 0.25));
        CHECK(near(mesh.cell_width(0, 1), 1.75));
        CHECK(near(mesh.cell_center(0, 0), mesh.face[0][0] + 0.125));
        mesh.validate_geometry();

        auto index = hd_core::Index{};
        index[0] = -1;
        expect_throw<std::out_of_range>([&] { mesh.at(index); });
        index[0] = mesh.nactive[0];
        expect_throw<std::out_of_range>([&] { mesh.at_active(index); });
        expect_throw<std::out_of_range>([&] { mesh.coordinates(count); });
        expect_throw<std::out_of_range>([&] { mesh.cell_width(hd_core::NDIM, 0); });
        expect_throw<std::out_of_range>([&] { mesh.cell_center(0, mesh.ntotal(0)); });
        expect_throw<std::invalid_argument>([&] { mesh.at({}); });
        auto invalid = mesh;
        invalid.face[0][1] = invalid.face[0][0];
        expect_throw<std::invalid_argument>([&] { invalid.allocate(); });
        invalid = mesh;
        invalid.face[0].pop_back();
        expect_throw<std::invalid_argument>([&] { invalid.validate_geometry(); });
        invalid = mesh;
        invalid.cells.pop_back();
        expect_throw<std::logic_error>([&] { invalid.coordinates(0); });
        invalid = mesh;
        invalid.nghost[0] = std::numeric_limits<int>::max();
        expect_throw<std::overflow_error>([&] { invalid.ntotal(0); });

        // Poison every ghost, including edges/corners, and populate active cells.
        auto states = mesh;
        std::vector<std::size_t> active_cells;
        for (std::size_t flat = 0; flat < states.cells.size(); ++flat)
        {
            auto& cell = states.cells[flat];
            cell.U[hd_core::RHO] = -1.0;
            cell.gamma = std::numeric_limits<double>::quiet_NaN();
            const auto local = states.coordinates_active(flat);
            bool active = true;
            for (int d = 0; d < hd_core::NDIM; ++d)
                active &= local[d] >= 0 && local[d] < states.nactive[d];
            if (!active) continue;
            cell.U[hd_core::RHO] = 1.0;
            cell.U[hd_core::ENERGY] = 2.0;
            cell.gamma = flat % 2 ? 1.4 : 5.0 / 3.0;
            active_cells.push_back(flat);
        }
        const auto& const_states = states;
        const_states.validate_active(); // Invalid ghosts must not be examined.
        for (auto flat : active_cells)
        {
            states.cells[flat].U[hd_core::RHO] = -1.0;
            expect_throw<std::invalid_argument>([&] { states.validate_active(); });
            states.cells[flat].U[hd_core::RHO] = 1.0;
        }
        auto& last_active = states.cells[active_cells.back()];
        const auto saved = last_active;
        last_active.U[hd_core::ENERGY] = -1.0;
        expect_throw<std::invalid_argument>([&] { states.validate_active(); });
        last_active = saved;
        last_active.U[hd_core::MOMDEN(0)] = std::numeric_limits<double>::quiet_NaN();
        expect_throw<std::invalid_argument>([&] { states.validate_active(); });
        last_active = saved;
        last_active.gamma = 1.0;
        expect_throw<std::invalid_argument>([&] { states.validate_active(); });
        last_active = saved;
        last_active.U = {}; // Preserve Conserved::validate's vacuum convention.
        states.validate_active();

        invalid = states;
        invalid.cells.clear();
        expect_throw<std::logic_error>([&] { invalid.validate_active(); });
        invalid = states;
        invalid.cells.pop_back();
        expect_throw<std::logic_error>([&] { invalid.validate_active(); });

        // End-to-end builder: coordinates must exist before allocate() validates.
        config::SimulationConfig cfg;
        double domain_volume = 1.0;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            cfg.nactive[d] = 2 + d;
            cfg.nghost[d] = 1 + d;
            cfg.lower[d] = -1.0;
            cfg.upper[d] = cfg.lower[d] + (2 + d) * (1.5 + d);
            cfg.boundary_lower[d] = cfg.boundary_upper[d] = "periodic";
            domain_volume *= cfg.upper[d] - cfg.lower[d];
        }
        hd_core::MeshBlock built;
        built.build_mesh(cfg);
        CHECK(built.volume_cell.size() == built.cells.size());
        double active_volume = 0.0;
        for (std::size_t flat = 0; flat < built.cells.size(); ++flat)
        {
            const auto local = built.coordinates_active(flat);
            bool active = true;
            for (int d = 0; d < hd_core::NDIM; ++d)
                active &= local[d] >= 0 && local[d] < built.nactive[d];
            if (active) active_volume += built.volume_cell[flat];
        }
        CHECK(near(active_volume, domain_volume));

        // Nonuniform grid, different widths per axis. Rebuilding must use the
        // current faces and preserve physical states, not recreate the cfg grid.
        built.cells[0].U[hd_core::RHO] = 7.0;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            built.face[d][0] = -2.0;
            for (int i = 0; i < built.ntotal(d); ++i)
                built.face[d][i + 1] = built.face[d][i] + (d + 1) * (i + 1);
        }
        built.build_geometry();
        CHECK(built.cells[0].U[hd_core::RHO] == 7.0);
        const auto shape = built.ntotal();
        for (std::size_t flat = 0; flat < built.volume_cell.size(); ++flat)
        {
            const auto index = built.coordinates(flat);
            double expected_volume = 1.0;
            for (int d = 0; d < hd_core::NDIM; ++d)
                expected_volume *= (d + 1) * (index[d] + 1);
            CHECK(near(built.volume_cell[flat], expected_volume));
        }
        for (int normal = 0; normal < hd_core::NDIM; ++normal)
        {
            auto face_shape = shape;
            ++face_shape[normal];
            std::size_t face_count = 1;
            for (const int n : face_shape) face_count *= n;
            const auto& areas = built.face_area[normal];
            CHECK(areas.shape == face_shape);
            CHECK(areas.faces.size() == face_count);
            for (std::size_t flat = 0; flat < face_count; ++flat)
            {
                auto remainder = flat;
                double expected_area = 1.0;
                for (int d = 0; d < hd_core::NDIM; ++d)
                {
                    const auto i = remainder % face_shape[d];
                    remainder /= face_shape[d];
                    if (d != normal) expected_area *= (d + 1) * (i + 1);
                }
                CHECK(near(areas.faces[flat], expected_area));
                CHECK(near(areas.at(areas.coordinates(flat)), expected_area));
            }

            // Fluxes remain active-only. Areas include ghosts, so the same
            // physical face is found by adding the ghost offset in every axis.
            auto flux_shape = built.nactive;
            ++flux_shape[normal];
            hydro::FaceFlux flux;
            flux.allocate(flux_shape);
            for (std::size_t f = 0; f < flux.faces.size(); ++f)
            {
                auto storage_face = flux.coordinates(f);
                double expected_area = 1.0;
                for (int d = 0; d < hd_core::NDIM; ++d)
                {
                    storage_face[d] += built.nghost[d];
                    if (d != normal) expected_area *= (d + 1) * (storage_face[d] + 1);
                }
                CHECK(near(areas.at(storage_face), expected_area));
            }

            // For a cell I, I is its lower face and I + e_normal its upper
            // face; the latter is also the next cell's lower face.
            for (std::size_t flat = 0; flat < built.cells.size(); ++flat)
            {
                auto upper = built.coordinates(flat);
                const auto lower_flat = areas.flat_index(upper);
                ++upper[normal];
                const auto upper_flat = areas.flat_index(upper);
                std::size_t normal_stride = 1;
                for (int d = 0; d < normal; ++d) normal_stride *= face_shape[d];
                CHECK(upper_flat == lower_flat + normal_stride);
                CHECK(near(areas.faces[lower_flat], areas.at(upper)));
            }
        }

        // A failed rebuild does not replace either cache.
        const auto saved_volumes = built.volume_cell;
        const auto saved_areas = built.face_area;
        built.face[0][1] = built.face[0][0];
        expect_throw<std::invalid_argument>([&] { built.build_geometry(); });
        CHECK(built.volume_cell == saved_volumes);
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            CHECK(built.face_area[d].shape == saved_areas[d].shape);
            CHECK(built.face_area[d].faces == saved_areas[d].faces);
        }

        cfg.geometry = "spherical";
        cfg.lower[0] = 0.0;
        cfg.upper[0] = 2.0;
        cfg.nactive[0] = 2;
        cfg.boundary_lower[0] = "reflective";
        hd_core::MeshBlock sphere;
        if constexpr (hd_core::NDIM == 1)
        {
            sphere.build_mesh(cfg); // Faces: -1, 0, 1, 2, 3.
            const double four_pi = 4.0 * std::acos(-1.0);
            CHECK(sphere.volume_cell.size() == 4);
            CHECK(sphere.face_area[0].faces.size() == 5);
            CHECK(sphere.face_area[0].at({1}) == 0.0);
            CHECK(near(sphere.face_area[0].at({0}), four_pi));
            CHECK(near(sphere.face_area[0].at({3}), 4.0 * four_pi));
            CHECK(near(sphere.face_area[0].at({4}), 9.0 * four_pi));
            CHECK(near(sphere.volume_cell[0], four_pi / 3.0));
            CHECK(near(sphere.volume_cell[1], sphere.volume_cell[0]));
            CHECK(near(sphere.volume_cell[2], 7.0 * four_pi / 3.0));
            CHECK(near(sphere.volume_cell[1] + sphere.volume_cell[2],
                       8.0 * four_pi / 3.0));

            // Geometry rebuild works without allocated physical states and
            // remains accurate for thin shells far from the origin.
            sphere.cells.clear();
            const double radius = 1e8;
            sphere.face[0] = {radius - 1, radius, radius + 1, radius + 2, radius + 3};
            sphere.build_geometry();
            CHECK(sphere.cells.empty());
            CHECK(near(sphere.volume_cell[1],
                       four_pi * (radius * radius + radius + 1.0 / 3.0)));
            const auto shell_volumes = sphere.volume_cell;
            const auto shell_areas = sphere.face_area;
            sphere.set_linear_uniform_faces(0, 0.0, 1e200);
            expect_throw<std::invalid_argument>([&] { sphere.build_geometry(); });
            CHECK(sphere.volume_cell == shell_volumes);
            CHECK(sphere.face_area[0].shape == shell_areas[0].shape);
            CHECK(sphere.face_area[0].faces == shell_areas[0].faces);

            cfg.lower[0] = -1.0;
            expect_throw<std::invalid_argument>([&] { sphere.build_mesh(cfg); });
        }
        else
        {
            expect_throw<std::invalid_argument>([&] { sphere.build_mesh(cfg); });
        }
        built.geometry = hd_core::Geometry::Cylindrical;
        built.set_linear_uniform_faces(0, cfg.lower[0], cfg.upper[0]);
        expect_throw<std::invalid_argument>([&] { built.build_geometry(); });
    });
}
