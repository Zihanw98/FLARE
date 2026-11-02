#include "test_support.h"
#include "core/mesh.h"
#include <limits>

int main()
{
    return run_test("mesh indexing and nonuniform geometry", [] {
        hd_core::MeshBlock mesh{};
        std::size_t count = 1;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            mesh.nactive[d] = 3 + d;
            mesh.nghost[d] = 1 + d;
            mesh.set_uniform_faces(d, -1.0, 2.0 + d);
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
    });
}
