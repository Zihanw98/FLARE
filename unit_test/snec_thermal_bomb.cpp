// Case-only adapter: SNEC shells and a Lagrangian-mass thermal source on FLARE's
// fixed Eulerian grid. No framework source or physical constant is modified.
#include "user/user_api.h"
#include "source/source.h"
#include "microphysics/eos/eos.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
struct Shell { double left, right, rho, pressure, velocity, mass_left, mass_right; };
std::vector<Shell> shells;
std::vector<double> bomb_edges, bomb_specific_weight;
double bomb_energy, bomb_end, bomb_dt;

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg) {
    using namespace hd_core;
    if (mesh.geometry != Geometry::Spherical || cfg.gamma != 1.4 ||
        hd_microphysics::eos_type != hd_microphysics::Equation_of_State::IdealGas_Default ||
        !cfg.spherical_gravity_enabled || !cfg.spherical_self_gravity_enabled ||
        cfg.constant_heating_enabled || cfg.constant_gravity_enabled ||
        mesh.bc_lower[0] != boundary::BoundaryCondition::Reflective ||
        mesh.bc_upper[0] != boundary::BoundaryCondition::User || cfg.time_integrator != "rk2")
        throw std::invalid_argument("SNEC comparison requires spherical ideal-gas RK2, self gravity, reflective/user boundaries");
    std::ifstream input("initial_shells.dat");
    int n;
    if (!(input >> n) || n != mesh.nactive[0]) throw std::runtime_error("Invalid initial_shells.dat count");
    shells.resize(n);
    for (auto& s : shells) {
        if (!(input >> s.left >> s.right >> s.rho >> s.pressure >> s.velocity >> s.mass_left >> s.mass_right) ||
            !(s.right > s.left && s.rho > 0 && s.pressure > 0 && s.mass_right > s.mass_left))
            throw std::runtime_error("Invalid SNEC shell");
    }
    if (std::abs(cfg.spherical_central_mass / shells.front().mass_left - 1) > 1e-12)
        throw std::runtime_error("Central mass differs from SNEC mass cut");
    const int ng = mesh.nghost[0];
    for (int i = 0; i < n; ++i) mesh.face[0][ng+i] = shells[i].left;
    mesh.face[0][ng+n] = shells.back().right;
    for (int j = 1; j <= ng; ++j) {
        mesh.face[0][ng-j] = shells.front().left-j*(shells.front().right-shells.front().left);
        mesh.face[0][ng+n+j] = shells.back().right+j*(shells.back().right-shells.back().left);
    }
    mesh.allocate();
    mesh.build_geometry();
    std::filesystem::create_directories(cfg.output_directory);
    std::ofstream grid(cfg.output_directory+"/grid.csv");
    grid << std::setprecision(17) << "r_left,r_right,volume,mass_left,mass_right\n";
    for (int i = 0; i < n; ++i) {
        const auto& s = shells[i];
        Primitives w{};
        w[RHO] = s.rho; w[PRE] = s.pressure; w[VEL(0)] = s.velocity;
        mesh.at_active({i}) = Cell(w, cfg.gamma);
        const double vol = mesh.volume_cell[mesh.flat_index({i},true)];
        if (std::abs(s.rho*vol/(s.mass_right-s.mass_left)-1) > 2e-9)
            throw std::runtime_error("SNEC mass/volume mapping mismatch");
        grid << s.left << ',' << s.right << ',' << vol << ',' << s.mass_left << ',' << s.mass_right << '\n';
    }
    std::ifstream bomb("bomb.dat");
    int nb;
    if (!(bomb >> bomb_energy >> bomb_end >> bomb_dt >> nb) || nb < 1 || bomb_energy < 0 || bomb_end <= 0 || bomb_dt <= 0)
        throw std::runtime_error("Invalid bomb.dat");
    bomb_edges.resize(nb+1); bomb_specific_weight.resize(nb);
    double norm = 0;
    for (int j = 0; j < nb; ++j) {
        double lo, hi, q;
        if (!(bomb >> lo >> hi >> q) || hi <= lo || q <= 0 || (j && lo != bomb_edges[j]))
            throw std::runtime_error("Invalid bomb mass support");
        bomb_edges[j] = lo; bomb_edges[j+1] = hi; bomb_specific_weight[j] = q;
        norm += (hi-lo)*q;
    }
    if (std::abs(norm-1) > 1e-12) throw std::runtime_error("Bomb weights not normalized");
}

void register_sources(const config::SimulationConfig& cfg, const hd_core::MeshBlock& mesh, hd_source::SourceList& sources) {
    // Integrate SNEC's piecewise-constant specific heating over the *current*
    // enclosed-mass interval of each Eulerian cell. Fixed radial heating would
    // heat the wrong material as soon as the ejecta start moving.
    sources.push_back({"snec_mass_bomb", hd_source::SourceKind::Thermal,
        [](const hd_core::MeshBlock& mesh, double time, const config::SimulationConfig& cfg,
           std::vector<hd_core::Conserved>& rhs, double& dt_limit) {
            if (time > bomb_end) return;
            dt_limit = std::min(dt_limit, bomb_dt);
            const double c = std::log(100.0)/bomb_end;
            const double luminosity = bomb_energy*c*std::exp(-c*time)/0.99;
            double mass = cfg.spherical_central_mass;
            int j = 0;
            for (auto flat : mesh.active_to_storage) {
                const double dm = mesh.cells[flat].U[hd_core::RHO]*mesh.volume_cell[flat];
                const double end = mass+dm;
                while (j < static_cast<int>(bomb_specific_weight.size()) && bomb_edges[j+1] <= mass) ++j;
                double fraction = 0;
                for (int k = j; k < static_cast<int>(bomb_specific_weight.size()) && bomb_edges[k] < end; ++k)
                    fraction += std::max(0.0, std::min(end,bomb_edges[k+1])-std::max(mass,bomb_edges[k]))*bomb_specific_weight[k];
                const double rate = luminosity*fraction/mesh.volume_cell[flat];
                rhs[flat][hd_core::ENERGY] += rate;
                rhs[flat][hd_core::INTERNAL_ENERGY] += rate;
                if (rate > 0) dt_limit = std::min(dt_limit, 0.05*mesh.cells[flat].U.internal_energy()/rate);
                mass = end;
                if (mass >= bomb_edges.back()) break;
            }
        }, true});
    // Exercise the actual source on both the initial grid and a deliberately
    // redistributed mass field. Integral power must be independent of where
    // Lagrangian heating-bin boundaries fall inside Eulerian cells.
    for (int rearrange = 0; rearrange < 2; ++rearrange) {
        auto probe = mesh;
        if (rearrange) {
            const int count = std::min(40,mesh.nactive[0]);
            std::vector<double> dm(count);
            for (int i = 0; i < count; ++i) {
                auto f = mesh.active_to_storage[i];
                dm[i] = mesh.cells[f].U[hd_core::RHO]*mesh.volume_cell[f];
            }
            for (int i = 0; i < count; ++i) {
                auto f = probe.active_to_storage[i];
                probe.cells[f].U[hd_core::RHO] = dm[count-1-i]/probe.volume_cell[f];
            }
        }
        for (double t : {0.0,0.5*bomb_end,1.1*bomb_end}) {
            std::vector<hd_core::Conserved> rates(probe.cells.size());
            double limit = std::numeric_limits<double>::infinity();
            sources.back().func(probe,t,cfg,rates,limit);
            double power = 0;
            for (auto f : probe.active_to_storage) power += rates[f][hd_core::ENERGY]*probe.volume_cell[f];
            const double c = std::log(100.0)/bomb_end;
            const double expected = t > bomb_end ? 0 : bomb_energy*c*std::exp(-c*t)/0.99;
            if (std::abs(power-expected) > 2e-12*std::max(1.0,expected))
                throw std::runtime_error("Bomb source integral failed on redistributed mass");
        }
    }
    std::ofstream checks(cfg.output_directory+"/source_checks.json");
    checks << "{\"initial_and_redistributed_mass_power\":true,\"after_bomb_zero\":true}\n";
}

void outer_boundary(hd_core::MeshBlock& mesh, double, int direction, boundary::BoundarySide side) {
    if (direction != 0 || side != boundary::BoundarySide::Upper) throw std::runtime_error("Unexpected boundary");
    // Free outflow, far outside the early shock. SNEC has a moving free surface.
    const int end = mesh.nghost[0]+mesh.nactive[0];
    for (int i = end; i < mesh.ntotal(0); ++i) mesh.at({i}) = mesh.at({end-1});
}
}

hd_user::UserHooks hd_user::make_user_hooks() {
    UserHooks hooks;
    hooks.initialize = initialize;
    hooks.register_sources = register_sources;
    hooks.userbc = outer_boundary;
    return hooks;
}
