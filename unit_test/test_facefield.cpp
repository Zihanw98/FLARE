#include "test_support.h"
#include "hydro/facefield.h"

int main()
{
    return run_test("face-field shapes, indexing and allocation", [] {
        hd_core::Index active{};
        for (int d = 0; d < hd_core::NDIM; ++d) active[d] = 2 + d;
        for (int normal = 0; normal < hd_core::NDIM; ++normal)
        {
            auto shape = active;
            ++shape[normal]; // Only the normal direction has one extra face.
            hydro::FaceField field;
            field.allocate(shape);
            std::size_t count = 1;
            for (int extent : shape) count *= extent;
            CHECK(field.faces.size() == count);
            const auto& const_field = field;
            for (std::size_t flat = 0; flat < count; ++flat)
            {
                const auto index = field.coordinates(flat);
                std::size_t expected = 0, stride = 1;
                for (int d = 0; d < hd_core::NDIM; ++d)
                {
                    CHECK(index[d] >= 0 && index[d] < shape[d]);
                    expected += index[d] * stride;
                    stride *= shape[d];
                }
                CHECK(expected == flat);
                CHECK(field.flat_index(index) == flat);
                CHECK(&field.at(index) == &field.faces[flat]);
                CHECK(&const_field.at(index) == &field.faces[flat]);
                for (double value : field.faces[flat].data) CHECK(value == 0.0);
                field.at(index)[hd_core::RHO] = flat + 1.0;
            }
            auto invalid = shape;
            invalid[0] = 0;
            expect_throw<std::invalid_argument>([&] { field.allocate(invalid); });
            CHECK(field.shape == shape && field.faces.size() == count);
            CHECK(field.faces.back()[hd_core::RHO] == count);
            expect_throw<std::out_of_range>([&] { field.coordinates(count); });
            expect_throw<std::out_of_range>([&] { field.at(shape); });
            expect_throw<std::invalid_argument>([&] { field.at({}); });
            invalid.fill(0);
            invalid[0] = -1;
            expect_throw<std::out_of_range>([&] { field.at(invalid); });
            field.allocate(shape); // Reallocation clears previous fluxes.
            for (const auto& face : field.faces)
                for (double value : face.data) CHECK(value == 0.0);
            field.faces.pop_back();
            expect_throw<std::logic_error>([&] { field.coordinates(0); });
        }
    });
}
