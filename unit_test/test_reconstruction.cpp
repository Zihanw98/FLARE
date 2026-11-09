#include "test_support.h"
#include "hydro/reconstruction.h"

namespace
{
using hydro::Reconstructor;
using hydro::Limiter;

hd_core::MeshBlock make_mesh()
{
    hd_core::MeshBlock mesh;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        mesh.nactive[d] = 4;
        mesh.nghost[d] = 2;
        mesh.face[d].push_back(0);
        for (double width : {0.2, 0.3, 0.4, 2.0, 0.1, 0.5, 0.6, 0.7})
            mesh.face[d].push_back(mesh.face[d].back() + width);
    }
    mesh.allocate();
    mesh.build_geometry();
    return mesh;
}

// A large positive gradient on one side makes unlimited PLM negative on the
// other face. Test density and pressure independently, including tiny scales.
void constrained(int direction, int component, double scale)
{
    auto mesh = make_mesh();
    auto index = mesh.nghost;
    index[direction] = 3;
    const double width = mesh.cell_width(direction, index);
    const double spacing = mesh.cell_center(direction, 4) - mesh.cell_center(direction, 3);
    const auto fill = [&](double overshoot) {
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            const int i = mesh.coordinates(flat)[direction];
            hd_core::Primitives w{};
            w[hd_core::RHO] = w[hd_core::PRE] = scale;
            if (i >= 4) w[component] *= 1 + 4 * spacing / width * (1 + overshoot);
            // All velocity components must share the pressure/density factor.
            for (int d = 0; d < hd_core::NDIM; ++d)
                w[hd_core::VEL(d)] = 0.01 * (d + 1) * (i - 3);
            mesh.cells[flat] = hd_core::Cell(w, 2.0);
        }
    };
    const auto faces = [&](Reconstructor method) {
        hydro::ReconstructionCache cache(mesh, method, Limiter::none, 1.5);
        auto right = index;
        ++right[direction];
        const auto a = cache.face(direction, index);
        const auto b = cache.face(direction, right);
        const auto direct_a = hydro::Reconstruct(method, mesh, direction, index, Limiter::none);
        const auto direct_b = hydro::Reconstruct(method, mesh, direction, right, Limiter::none);
        CHECK(a.first.data == direct_a.first.data && a.second.data == direct_a.second.data);
        CHECK(b.first.data == direct_b.first.data && b.second.data == direct_b.second.data);
        CHECK(cache.face(direction, index).second.data == a.second.data);
        return std::make_pair(a.second, b.first); // Both faces of the same cell.
    };

    fill(0.5);
    const auto old_cells = mesh.cells;
    const auto center = mesh.at(index).U.to_Primtivies(2.0);
    const auto old_faces = faces(Reconstructor::PLM);
    CHECK(old_faces.first.data == center.data); // Old whole-state fallback.
    const auto fixed = faces(Reconstructor::PLM_corrected);
    fixed.first.validate(2.0);
    fixed.second.validate(2.0);
    for (int k : {hd_core::RHO, hd_core::PRE})
    {
        CHECK(fixed.first[k] >= 0.99e-10 * center[k]);
        CHECK(fixed.second[k] >= 0.99e-10 * center[k]);
        CHECK(near((fixed.first[k] + fixed.second[k]) / center[k], 2.0));
    }
    const double theta = (fixed.second[component] - center[component])
                       / (old_faces.second[component] - center[component]);
    CHECK(theta > 0 && theta < 1);
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        const int k = hd_core::VEL(d);
        CHECK(near(fixed.second[k], theta * old_faces.second[k]));
        CHECK(near(fixed.first[k], -fixed.second[k]));
        CHECK(fixed.first[k] != center[k]); // Must not silently fall back.
    }
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
        CHECK(mesh.cells[i].U.data == old_cells[i].U.data);

    fill(-1e-7);
    const auto below = faces(Reconstructor::PLM_corrected);
    const auto old_below = faces(Reconstructor::PLM);
    CHECK(below.first.data == old_below.first.data);
    CHECK(below.second.data == old_below.second.data);
    fill(1e-7);
    const auto above = faces(Reconstructor::PLM_corrected);
    CHECK(std::abs(above.first[component] - below.first[component]) / scale < 2e-7);
    for (int d = 0; d < hd_core::NDIM; ++d)
        CHECK(std::abs(above.first[hd_core::VEL(d)] - below.first[hd_core::VEL(d)]) < 1e-7);
}

void cache_and_unchanged()
{
    auto mesh = make_mesh();
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
    {
        hd_core::Primitives w{};
        w[hd_core::RHO] = 1 + (i % 5) * 0.01;
        w[hd_core::PRE] = 2 + (i % 7) * 0.01;
        mesh.cells[i] = hd_core::Cell(w, 1.4);
    }
    for (auto limiter : {Limiter::none, Limiter::minmod, Limiter::MC,
                        Limiter::van_Leer, Limiter::MCbeta, Limiter::superbee})
    {
        hydro::ReconstructionCache cache(mesh, Reconstructor::PLM_corrected, limiter, 1.7);
        for (auto flat : mesh.get_all_active_flats())
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            auto index = mesh.coordinates(flat);
            const auto old = hydro::Reconstruct(Reconstructor::PLM, mesh, d, index, limiter, 1.7);
            const auto fixed = hydro::Reconstruct(Reconstructor::PLM_corrected, mesh, d, index, limiter, 1.7);
            const auto cached = cache.face(d, index);
            CHECK(old.first.data == fixed.first.data && old.second.data == fixed.second.data);
            CHECK(cached.first.data == fixed.first.data && cached.second.data == fixed.second.data);
        }
    }
}

void monotone_nonuniform(int direction)
{
    auto mesh = make_mesh();
    for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
    {
        const int i = mesh.coordinates(flat)[direction];
        hd_core::Primitives w{};
        w[hd_core::RHO] = w[hd_core::PRE] = i < 3 ? 0.01 : i == 3 ? 1.0 : 10.0;
        mesh.cells[flat] = hd_core::Cell(w, 2.0);
    }
    auto index = mesh.nghost;
    index[direction] = 3;
    const auto original = hydro::Reconstruct(Reconstructor::PLM, mesh, direction, index, Limiter::MC);
    CHECK(original.second.data == mesh.at(index).U.to_Primtivies(2.0).data);
    for (auto limiter : {Limiter::none, Limiter::minmod, Limiter::MC,
                        Limiter::van_Leer, Limiter::MCbeta, Limiter::superbee})
    {
        hydro::ReconstructionCache cache(mesh, Reconstructor::PLM_corrected, limiter, 1.7);
        const auto direct = hydro::Reconstruct(Reconstructor::PLM_corrected, mesh, direction, index, limiter, 1.7);
        const auto fixed = cache.face(direction, index);
        CHECK(fixed.first.data == direct.first.data && fixed.second.data == direct.second.data);
        fixed.first.validate(2.0);
        fixed.second.validate(2.0);
        CHECK(fixed.second[hd_core::PRE] < 1.0); // Preserve the gradient.
    }
}
}

int main()
{
    return run_test("PLM_corrected positivity, continuity and compatibility", [] {
        CHECK(hydro::parse_reconstructor("pLm_CoRrEcTeD") == Reconstructor::PLM_corrected);
        CHECK(hydro::parse_reconstructor("PLM") == Reconstructor::PLM);
        CHECK(hydro::parse_reconstructor("PCM") == Reconstructor::PCM);
        for (int d = 0; d < hd_core::NDIM; ++d)
        for (int component : {hd_core::RHO, hd_core::PRE})
        for (double scale : {1.0, 1e-25})
            constrained(d, component, scale);
        for (int d = 0; d < hd_core::NDIM; ++d)
            monotone_nonuniform(d);
        cache_and_unchanged();
    });
}
