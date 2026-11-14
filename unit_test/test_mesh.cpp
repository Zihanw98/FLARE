#include "test_support.h"
#include "core/mesh.h"
#include "config/config.h"
#include "hydro/facefield.h"
#include <limits>

namespace
{
void check_topology(const hd_core::MeshBlock& mesh)
{
    const auto invalid = hd_core::MeshBlock::invalid_active_index;
    const auto active_count = mesh.active_to_storage.size();
    CHECK(mesh.total_count == mesh.cells.size());
    CHECK(mesh.active_count == active_count);
    CHECK(mesh.cell_faces.size() == active_count);
    std::size_t expected_faces = 0;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        const bool periodic = mesh.bc_lower[d] == boundary::BoundaryCondition::Periodic;
        expected_faces += active_count + (periodic ? 0 : active_count / mesh.nactive[d]);
    }
    CHECK(mesh.face_info.size() == expected_faces);
    std::vector<int> references(mesh.face_info.size(), 0);
    for (std::size_t cell = 0; cell < active_count; ++cell)
    {
        const auto storage = mesh.coordinates(mesh.active_to_storage[cell]);
        const auto local = mesh.coordinates_active(mesh.active_to_storage[cell]);
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            const bool periodic = mesh.bc_lower[d] == boundary::BoundaryCondition::Periodic;
            const auto& pair = mesh.cell_faces[cell][d];
            const auto& lower = mesh.face_info.at(pair.lower);
            const auto& upper = mesh.face_info.at(pair.upper);
            CHECK(lower.direction == d && upper.direction == d);
            CHECK(lower.right_active == cell);
            CHECK(upper.left_active == cell);
            CHECK(lower.storage_coordinates == storage);
            auto upper_coordinate = storage;
            ++upper_coordinate[d];
            if (periodic && local[d] == mesh.nactive[d] - 1)
                upper_coordinate[d] = mesh.nghost[d];
            CHECK(upper.storage_coordinates == upper_coordinate);

            for (int side : {-1, 1})
            {
                auto neighbor_coordinate = local;
                neighbor_coordinate[d] += side;
                const bool outside = neighbor_coordinate[d] < 0 ||
                    neighbor_coordinate[d] >= mesh.nactive[d];
                std::size_t neighbor = invalid;
                if (periodic || !outside)
                {
                    if (neighbor_coordinate[d] < 0) neighbor_coordinate[d] = mesh.nactive[d] - 1;
                    if (neighbor_coordinate[d] == mesh.nactive[d]) neighbor_coordinate[d] = 0;
                    neighbor = mesh.storage_to_active.at(mesh.flat_index(neighbor_coordinate, true));
                    CHECK(neighbor != invalid);
                    const auto& other = mesh.cell_faces.at(neighbor)[d];
                    CHECK((side < 0 ? other.upper : other.lower) ==
                          (side < 0 ? pair.lower : pair.upper));
                }
                CHECK((side < 0 ? lower.left_active : upper.right_active) == neighbor);
            }
            ++references[pair.lower];
            ++references[pair.upper];
            if constexpr (hd_core::NDIM == 1)
            {
                // Preserve the old LTS face numbering, including periodic face 0.
                CHECK(pair.lower == cell);
                CHECK(pair.upper == (periodic && cell + 1 == active_count ? 0 : cell + 1));
            }
        }
    }
    for (std::size_t face = 0; face < mesh.face_info.size(); ++face)
    {
        const auto& info = mesh.face_info[face];
        const int sides = (info.left_active != invalid) + (info.right_active != invalid);
        CHECK(sides >= 1);
        CHECK(references[face] == sides);
    }
}

void check_same_topology(const hd_core::MeshBlock& a, const hd_core::MeshBlock& b)
{
    CHECK(a.total_count == b.total_count);
    CHECK(a.active_count == b.active_count);
    CHECK(a.active_to_storage == b.active_to_storage);
    CHECK(a.storage_to_active == b.storage_to_active);
    CHECK(a.face_info.size() == b.face_info.size());
    CHECK(a.cell_faces.size() == b.cell_faces.size());
    for (std::size_t f = 0; f < a.face_info.size(); ++f)
    {
        CHECK(a.face_info[f].direction == b.face_info[f].direction);
        CHECK(a.face_info[f].storage_coordinates == b.face_info[f].storage_coordinates);
        CHECK(a.face_info[f].left_active == b.face_info[f].left_active);
        CHECK(a.face_info[f].right_active == b.face_info[f].right_active);
    }
    for (std::size_t c = 0; c < a.cell_faces.size(); ++c)
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            CHECK(a.cell_faces[c][d].lower == b.cell_faces[c][d].lower);
            CHECK(a.cell_faces[c][d].upper == b.cell_faces[c][d].upper);
        }
}
}

int main()
{
    return run_test("mesh indexing, geometry and active-state validation", [] {
        hd_core::MeshBlock mesh{};
        CHECK(mesh.total_count == 0);
        CHECK(mesh.active_count == 0);
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
        check_topology(mesh);
        auto bounded = mesh;
        bounded.bc_lower.fill(boundary::BoundaryCondition::Reflective);
        bounded.bc_upper.fill(boundary::BoundaryCondition::Reflective);
        bounded.build_topology();
        check_topology(bounded);
        CHECK(mesh.cells.size() == count);
        CHECK(mesh.storage_to_active.size() == count);
        std::size_t active_count = 0;
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
                CHECK(mesh.active_to_storage.at(active_count) == flat);
                CHECK(mesh.storage_to_active.at(flat) == active_count);
                CHECK(mesh.active_to_storage.at(mesh.storage_to_active[flat]) == flat);
                ++active_count;
                CHECK(mesh.flat_index(local, true) == flat);
                CHECK(&mesh.at_active(local) == &mesh.cells[flat]);
                CHECK(&const_mesh.at_active(local) == &mesh.cells[flat]);
            }
            else
                CHECK(mesh.storage_to_active.at(flat) == hd_core::MeshBlock::invalid_active_index);
        }
        CHECK(mesh.active_to_storage.size() == active_count);

        // The builder recreates all connectivity without reallocating cell storage.
        auto rebuilt = mesh;
        rebuilt.active_to_storage.clear();
        rebuilt.storage_to_active.clear();
        rebuilt.face_info.clear();
        rebuilt.cell_faces.clear();
        auto* const cell_storage = rebuilt.cells.data();
        rebuilt.cells[0].U[hd_core::RHO] = 7.0;
        rebuilt.build_topology();
        check_same_topology(rebuilt, mesh);
        CHECK(rebuilt.cells.data() == cell_storage);
        CHECK(rebuilt.cells[0].U[hd_core::RHO] == 7.0);

        // Reallocating a copied grid must rebuild its maps even when the
        // storage size is unchanged, without modifying the original grid.
        auto resized = mesh;
        resized.nactive[0] -= 2;
        ++resized.nghost[0];
        resized.set_linear_uniform_faces(0, -1.0, 2.0);
        resized.allocate();
        check_topology(resized);
        CHECK(resized.cells.size() == count);
        std::vector<std::size_t> expected_active;
        std::vector<std::size_t> expected_inverse(count, hd_core::MeshBlock::invalid_active_index);
        for (std::size_t flat = 0; flat < resized.cells.size(); ++flat)
            if (resized.is_active(flat))
            {
                expected_inverse[flat] = expected_active.size();
                expected_active.push_back(flat);
            }
        CHECK(resized.active_to_storage == expected_active);
        CHECK(resized.storage_to_active == expected_inverse);
        CHECK(mesh.active_to_storage.size() == active_count);
        CHECK(resized.active_to_storage != mesh.active_to_storage);
        CHECK(resized.storage_to_active != mesh.storage_to_active);
        CHECK(mesh.storage_to_active == rebuilt.storage_to_active);

        // Without ghosts both maps are the identity in every dimension.
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            resized.nactive[d] = 2 + d;
            resized.nghost[d] = 0;
            resized.set_linear_uniform_faces(d, -1.0, 2.0 + d);
        }
        resized.allocate();
        check_topology(resized);
        CHECK(resized.active_to_storage.size() == resized.cells.size());
        CHECK(resized.storage_to_active.size() == resized.cells.size());
        for (std::size_t flat = 0; flat < resized.cells.size(); ++flat)
        {
            CHECK(resized.active_to_storage[flat] == flat);
            CHECK(resized.storage_to_active[flat] == flat);
        }

        // Every combination of periodic/physical axes, with asymmetric shapes,
        // no ghosts or unequal ghost widths, and a one-cell periodic direction.
        for (int mask = 0; mask < (1 << hd_core::NDIM); ++mask)
        {
            for (int ghosts : {0, 1})
            {
                hd_core::MeshBlock topology;
                for (int d = 0; d < hd_core::NDIM; ++d)
                {
                    topology.nactive[d] = 1 + d;
                    topology.nghost[d] = ghosts * (d + 1);
                    topology.bc_lower[d] = topology.bc_upper[d] = (mask & (1 << d))
                        ? boundary::BoundaryCondition::Periodic
                        : boundary::BoundaryCondition::Reflective;
                    topology.set_linear_uniform_faces(d, -1.0, 2.0 + d);
                }
                topology.allocate();
                check_topology(topology);
                topology.build_geometry();
                for (const auto& info : topology.face_info)
                    CHECK(topology.face_area[info.direction].at(info.storage_coordinates) > 0.0);

                // A boundary-only change rebuilds faces without changing cells.
                const auto* storage = topology.cells.data();
                const auto active_map = topology.active_to_storage;
                const auto reverse_map = topology.storage_to_active;
                const auto old_face_count = topology.face_info.size();
                const bool was_periodic = topology.bc_lower[0] == boundary::BoundaryCondition::Periodic;
                topology.bc_lower[0] = topology.bc_upper[0] = was_periodic
                    ? boundary::BoundaryCondition::Reflective : boundary::BoundaryCondition::Periodic;
                topology.build_topology();
                check_topology(topology);
                CHECK(topology.cells.data() == storage);
                CHECK(topology.active_to_storage == active_map);
                CHECK(topology.storage_to_active == reverse_map);
                const auto boundary_plane = active_map.size() / topology.nactive[0];
                CHECK(topology.face_info.size() == (was_periodic
                    ? old_face_count + boundary_plane : old_face_count - boundary_plane));
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
        check_same_topology(invalid, mesh);
        invalid = mesh;
        invalid.face[0].pop_back();
        expect_throw<std::invalid_argument>([&] { invalid.validate_geometry(); });
        invalid = mesh;
        invalid.cells.pop_back();
        expect_throw<std::logic_error>([&] { invalid.coordinates(0); });
        expect_throw<std::logic_error>([&] { invalid.build_topology(); });
        check_same_topology(invalid, mesh);
        invalid = mesh;
        invalid.bc_lower[0] = boundary::BoundaryCondition::Reflective;
        expect_throw<std::invalid_argument>([&] { invalid.build_topology(); });
        check_same_topology(invalid, mesh);
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
        check_topology(built);
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
        cfg.boundary_lower[0] = cfg.boundary_upper[0] = "reflective";
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
