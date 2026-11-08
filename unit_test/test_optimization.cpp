#include "test_support.h"
#include "core/mesh.h"
#include "hydro/reconstruction.h"
#include "hydro/riemann/exact.h"
#include "hydro/wave_timestep.h"
#include "utils/performance.h"
#include <cstdlib>
#include <new>
#include <random>
#include <numeric>

namespace { bool track_allocations = false; std::size_t allocations = 0; }
void* operator new(std::size_t n)
{
    if (track_allocations) ++allocations;
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc();
}
[[gnu::noinline]] void operator delete(void* p) noexcept { std::free(p); }
[[gnu::noinline]] void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void* operator new[](std::size_t n) { return ::operator new(n); }
[[gnu::noinline]] void operator delete[](void* p) noexcept { ::operator delete(p); }
[[gnu::noinline]] void operator delete[](void* p, std::size_t) noexcept { ::operator delete(p); }

hd_core::Primitives state(double rho, double p, double u = 0, int d = 0)
{
    hd_core::Primitives w{};
    w[hd_core::RHO] = rho; w[hd_core::PRE] = p; w[hd_core::VEL(d)] = u;
    return w;
}

void conversions()
{
    std::mt19937 rng(5);
    for (int sample = 0; sample < 1000; ++sample)
    {
        const double gamma = 1.1 + (rng() % 800) * 0.001;
        auto w = state(0.01 + (rng() % 1000), 100 + (rng() % 1000));
        for (int d = 0; d < hd_core::NDIM; ++d) w[hd_core::VEL(d)] = (static_cast<int>(rng() % 100) - 50) * 0.01;
        track_allocations = true;
        const auto u = w.to_Conserved(gamma);
        const auto recovered = u.to_Primtivies(gamma);
        track_allocations = false;
        CHECK(allocations == 0);
        // Original vector formulas provide an independent byte-level reference.
        CHECK(u[hd_core::ENERGY] == w[hd_core::PRE] / (gamma - 1) + 0.5 * w[hd_core::RHO] * w.Velocity().squared_norm());
        CHECK(recovered[hd_core::PRE] == (gamma - 1) * (u[hd_core::ENERGY] - 0.5 * u[hd_core::RHO] * u.Velocity().squared_norm()));
    }
}

void reconstruction()
{
    hd_core::MeshBlock mesh;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        mesh.nactive[d] = 4; mesh.nghost[d] = 2;
        mesh.face[d].push_back(0);
        for (int i = 0; i < mesh.ntotal(d); ++i)
            mesh.face[d].push_back(mesh.face[d].back() + 0.1 * (1 + i * i));
    }
    mesh.allocate(); mesh.build_geometry();
    for (int pass = 0; pass < 2; ++pass)
    {
        for (std::size_t i = 0; i < mesh.cells.size(); ++i)
        {
            const auto index = mesh.coordinates(i);
            int ghost_axes = 0;
            for (int d = 0; d < hd_core::NDIM; ++d)
                ghost_axes += index[d] < 2 || index[d] >= 6;
            if (ghost_axes > 1) { mesh.cells[i].U[hd_core::RHO] = -1; continue; }
            // Discontinuous positive states exercise independent positivity fallback.
            auto w = state((i % 3 == 0 ? 0.01 : 10.0) + pass, (i % 5 == 0 ? 0.001 : 2.0) + pass);
            mesh.cells[i] = hd_core::Cell(w, 1.3 + 0.01 * (i % 20));
        }
        for (auto method : {hydro::Reconstructor::PCM, hydro::Reconstructor::PLM})
        for (auto limiter : {hydro::Limiter::none, hydro::Limiter::minmod, hydro::Limiter::MC,
                            hydro::Limiter::van_Leer, hydro::Limiter::MCbeta, hydro::Limiter::superbee})
        {
            hydro::ReconstructionCache cache(mesh, method, limiter, 1.7);
            for (int d = 0; d < hd_core::NDIM; ++d)
            {
                auto shape = mesh.nactive; ++shape[d];
                const auto count = std::accumulate(shape.begin(), shape.end(), std::size_t{1}, std::multiplies<std::size_t>());
                for (std::size_t f = 0; f < count; ++f)
                {
                    auto right = mesh.nghost; auto offset = f;
                    for (int n = 0; n < hd_core::NDIM; ++n) { right[n] += offset % shape[n]; offset /= shape[n]; }
                    const auto original = hydro::Reconstruct(method, mesh, d, right, limiter, 1.7);
                    for (int repeat = 0; repeat < 2; ++repeat)
                    {
                        const auto cached = cache.face(d, right);
                        CHECK(original.first.data == cached.first.data);
                        CHECK(original.second.data == cached.second.data);
                    }
                }
            }
        }
    }
}

void roots()
{
    using hydro::ExactPressureMethod;
    hydro::ExactSolver solver(1e-12, 512), reference(1e-12, 512, ExactPressureMethod::Bisection);
    std::mt19937 rng(37);
    std::uniform_real_distribution<double> unit(0, 1);
    std::uint64_t newton_iterations = 0, old_iterations = 0;
    double max_error = 0;
    for (int d = 0; d < hd_core::NDIM; ++d)
    for (int i = 0; i < 2000; ++i)
    {
        const double gl = 1.1 + 0.8 * unit(rng), gr = 1.1 + 0.8 * unit(rng);
        auto l = state(std::pow(10., -5 + 10 * unit(rng)), std::pow(10., -5 + 10 * unit(rng)));
        auto r = state(std::pow(10., -5 + 10 * unit(rng)), std::pow(10., -5 + 10 * unit(rng)));
        const double cl = std::sqrt(gl * l[hd_core::PRE] / l[hd_core::RHO]);
        const double cr = std::sqrt(gr * r[hd_core::PRE] / r[hd_core::RHO]);
        l[hd_core::VEL(d)] = (unit(rng) - 0.5) * 8 * cl;
        r[hd_core::VEL(d)] = (unit(rng) - 0.5) * 8 * cr;
        for (int t = 0; t < hd_core::NDIM; ++t) if (t != d)
        { l[hd_core::VEL(t)] = 0.3 * cl; r[hd_core::VEL(t)] = -0.7 * cr; }
        performance::configure(performance::Mode::Summary);
        const auto a = solver.flux(l, r, d, gl, gr);
        const auto speed_a = hydro::wave_timestep.face_speed;
        const auto stats = performance::exact_roots;
        CHECK(stats.solves == 1 && stats.failures == 0);
        CHECK(stats.initial_guesses + stats.newton_steps + stats.bisection_steps == stats.iterations);
        CHECK(std::accumulate(stats.histogram.begin(), stats.histogram.end(), std::uint64_t{0}) == 1);
        newton_iterations += stats.iterations;
        performance::configure(performance::Mode::Summary);
        const auto b = reference.flux(l, r, d, gl, gr);
        old_iterations += performance::exact_roots.iterations;
        CHECK(near(speed_a, hydro::wave_timestep.face_speed, 2e-10));
        // Flux components can nearly cancel (especially contact/energy flux).
        // Scale error by the physical flux scale, not by that tiny difference.
        const double speed = std::max(std::abs(l[hd_core::VEL(d)]) + cl, std::abs(r[hd_core::VEL(d)]) + cr);
        const double mass_scale = std::max(l[hd_core::RHO], r[hd_core::RHO]) * speed;
        const double momentum_scale = mass_scale * speed + std::max(l[hd_core::PRE], r[hd_core::PRE]);
        const double energy_scale = std::max(l.to_Conserved(gl)[hd_core::ENERGY] + l[hd_core::PRE],
            r.to_Conserved(gr)[hd_core::ENERGY] + r[hd_core::PRE]) * speed;
        for (int k = 0; k < hd_core::NHYDRO; ++k)
        {
            const double scale = k == hd_core::RHO ? mass_scale : k == hd_core::ENERGY ? energy_scale : momentum_scale;
            const double error = std::abs(a[k] - b[k]) / scale;
            max_error = std::max(error, max_error);
            CHECK(error < 2e-10);
        }
    }
    CHECK(newton_iterations < old_iterations);
    // Near-vacuum rarefactions; the pressure is many orders smaller than the input.
    for (double fraction : {0.5, 0.99, 0.999999, 1.01})
    {
        const double velocity = fraction * 2 * std::sqrt(1.4) / 0.4;
        const auto l = state(1, 1, -velocity), r = state(1, 1, velocity);
        const auto a = solver.flux(l, r, 0, 1.4, 1.4);
        const auto b = reference.flux(l, r, 0, 1.4, 1.4);
        for (int k = 0; k < hd_core::NHYDRO; ++k) CHECK(near(a[k], b[k], 1e-10));
    }
    performance::configure(performance::Mode::Summary);
    hydro::ExactSolver fail(1e-12, 1);
    expect_throw<std::runtime_error>([&] { fail.flux(state(1, 1000), state(1, 0.01), 0, 1.4, 1.4); });
    CHECK(performance::exact_roots.failures == 1);
    performance::configure(performance::Mode::Off);
    solver.flux(state(1, 1), state(1, 1), 0, 1.4, 1.4);
    CHECK(performance::exact_roots.solves == 0);
    std::cout << "Random roots: Newton iterations=" << newton_iterations << ", bisection="
              << old_iterations << ", max scaled flux error=" << max_error << '\n';
}

int main()
{
    return run_test("allocation-free conversions, cached reconstruction, safeguarded exact roots", [] {
        conversions(); reconstruction(); roots();
    });
}
