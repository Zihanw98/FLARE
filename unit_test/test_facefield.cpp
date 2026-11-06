#include "test_support.h"
#include "hydro/facefield.h"
#include <limits>
#include <type_traits>

static_assert(std::is_same_v<hydro::FaceFlux,
                             hd_core::FaceField<hd_core::Conserved>>);

template<typename T>
double& marker(T& value)
{
    if constexpr (std::is_same_v<T, double>) return value;
    else return value[hd_core::RHO];
}

template<typename T>
void check_zero(const T& value)
{
    if constexpr (std::is_same_v<T, double>) CHECK(value == 0.0);
    else for (const double component : value.data) CHECK(component == 0.0);
}

template<typename T>
void check_field()
{
        hd_core::Index active{};
        for (int d = 0; d < hd_core::NDIM; ++d) active[d] = 2 + d;
        for (int normal = 0; normal < hd_core::NDIM; ++normal)
        {
            auto shape = active;
            ++shape[normal]; // Only the normal direction has one extra face.
            hd_core::FaceField<T> field;
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
                check_zero(field.faces[flat]);
                marker(field.at(index)) = flat + 1.0;
            }
            auto invalid = shape;
            invalid[0] = 0;
            expect_throw<std::invalid_argument>([&] { field.allocate(invalid); });
            CHECK(field.shape == shape && field.faces.size() == count);
            CHECK(marker(field.faces.back()) == count);
            expect_throw<std::out_of_range>([&] { field.coordinates(count); });
            expect_throw<std::out_of_range>([&] { field.at(shape); });
            expect_throw<std::invalid_argument>([&] { field.at({}); });
            invalid.fill(0);
            invalid[0] = -1;
            expect_throw<std::out_of_range>([&] { field.at(invalid); });
            // Exercise both initializer-list accessors and flat_index overload.
            const auto check_list = [&](std::initializer_list<int> index) {
                const auto flat = field.flat_index(index);
                CHECK(&field.at(index) == &field.faces[flat]);
                CHECK(&const_field.at(index) == &field.faces[flat]);
                marker(field.at(index)) = 12.5;
                CHECK(marker(field.faces[flat]) == 12.5);
            };
            if constexpr (hd_core::NDIM == 1) check_list({1});
            else if constexpr (hd_core::NDIM == 2) check_list({1, 1});
            else check_list({1, 1, 1});
            if constexpr (hd_core::NDIM > 1)
            {
                auto huge = shape;
                huge.fill(std::numeric_limits<int>::max());
                expect_throw<std::length_error>([&] { field.allocate(huge); });
                CHECK(field.shape == shape && field.faces.size() == count);
            }
            field.allocate(shape); // Reallocation clears previous values.
            for (const auto& face : field.faces)
                check_zero(face);
            field.faces.pop_back();
            expect_throw<std::logic_error>([&] { field.coordinates(0); });
            expect_throw<std::logic_error>([&] { field.flat_index(hd_core::Index{}); });
            field.allocate(shape);
            field.shape[0] = 0;
            expect_throw<std::logic_error>([&] { field.coordinates(0); });
            expect_throw<std::logic_error>([&] { field.at(hd_core::Index{}); });
        }
}

int main()
{
    return run_test("face-field shapes, indexing and allocation", [] {
        check_field<hd_core::Conserved>();
        check_field<double>();
    });
}
