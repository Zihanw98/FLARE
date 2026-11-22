#include "test_support.h"
#include "time/implicit.h"
#include "core/mesh.h"
#include "microphysics/eos/eos.h"
#include "microphysics/composition/species.h"
#include "physcon/constants.h"
#include "user/user_api.h"
#include <fstream>
#include <limits>

namespace
{
using namespace hd_core;
using hd_time::TimeIntegrator;
constexpr double pi = 3.14159265358979323846;
constexpr double infinity = std::numeric_limits<double>::infinity();

// Isolate the transport operator from acoustic/advection errors.
class NoHydro : public hydro::RiemannSolver
{
public:
    double speed = 0.0;
    Conserved flux(const Primitives&, const Primitives&, int, double, double,
                   double& wave, double& velocity) const override
    {
        wave = speed;
        velocity = 0.0;
        return {};
    }
};

MeshBlock make_mesh(int n = 24, bool stretched = false, bool spherical = false)
{
    MeshBlock mesh;
    for (int d = 0; d < NDIM; ++d)
    {
        mesh.nactive[d] = d == 0 ? n : 2;
        mesh.nghost[d] = 1;
        mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
        mesh.set_linear_uniform_faces(d, 0.0, 1.0);
    }
    if (stretched || spherical)
    {
        auto& f = mesh.face[0];
        f[1] = spherical ? 1.0 : 0.0;
        for (int i = 0; i < n; ++i)
            f[i+2] = f[i+1] + ((stretched && i == n/2) ? 0.2 : 1.0)/n;
        f[0] = f[1] - (spherical ? f[2]-f[1] : f[n+1]-f[n]);
        f[n+2] = f[n+1] + (spherical ? f[n+1]-f[n] : f[2]-f[1]);
    }
    if (spherical)
    {
        mesh.geometry = Geometry::Spherical;
        mesh.bc_lower[0] = mesh.bc_upper[0] = boundary::BoundaryCondition::Reflective;
    }
    mesh.allocate();
    mesh.build_geometry();
    for (auto flat : mesh.active_to_storage)
    {
        auto& cell = mesh.cells[flat];
        const auto index = mesh.coordinates(flat);
        const double x = (mesh.cell_center(0, index)-mesh.face[0][1])/mesh.active_length(0);
        cell.gamma = 5.0/3.0;
        cell.U[RHO] = 1.0;
        cell.U[ENERGY] = cell.U[INTERNAL_ENERGY] = 2.0 + 0.1*std::cos(2*pi*x);
        if constexpr (NCOMP > 0) cell.U[COMP(0)] = 1.0;
    }
    boundary::FillBoundary(mesh, 0.0, {});
    hd_microphysics::refresh_thermodynamics(mesh);
    return mesh;
}

// Constant-coefficient energy diffusion, with insulating physical boundaries.
// Used as a scalar diffusion verification, not as a radiation physics closure.
hd_source::TransportFunction diffusion(double D, bool composition = false)
{
    return [=](const MeshBlock& mesh, std::size_t face, double,
               const config::SimulationConfig&, Conserved& flux, std::array<double,2>& rates) {
        const auto& f = mesh.face_info[face];
        if (f.direction != 0 || f.left_active == MeshBlock::invalid_active_index
            || f.right_active == MeshBlock::invalid_active_index) return;
        auto right = f.storage_coordinates, left = right;
        --left[0];
        const double distance = mesh.cell_center(0, right)-mesh.cell_center(0, left);
        const double area = mesh.face_area[0].at(right);
        rates[0] = area*D/(distance*mesh.volume_cell[mesh.active_to_storage[f.left_active]]);
        rates[1] = area*D/(distance*mesh.volume_cell[mesh.active_to_storage[f.right_active]]);
        const auto& L = mesh.at(left);
        const auto& R = mesh.at(right);
        if (composition)
        {
            if constexpr (NCOMP == 2)
            {
                const double J = -D*(R.U[COMP(0)]-L.U[COMP(0)])/distance;
                flux[COMP(0)] = J;
                flux[COMP(1)] = -J;
                const double T = 0.5*(L.temperature + R.temperature);
                const auto& species = hd_microphysics::species_table.species_active;
                const double h_difference = 2.5*physcon::R_g*T*(1.0/species[0].A-1.0/species[1].A);
                flux[ENERGY] = flux[INTERNAL_ENERGY] = h_difference*J;
                // This test's coupled energy flux is weaker than this conservative bound.
                rates[0] *= 8.0;
                rates[1] *= 8.0;
            }
        }
        else flux[ENERGY] = flux[INTERNAL_ENERGY] =
            -D*(R.U[INTERNAL_ENERGY]-L.U[INTERNAL_ENERGY])/distance;
    };
}

Conserved integral(const MeshBlock& mesh)
{
    Conserved value{};
    for (auto flat : mesh.active_to_storage) value += mesh.cells[flat].U*mesh.volume_cell[flat];
    return value;
}

config::SimulationConfig settings()
{
    config::SimulationConfig cfg;
    cfg.cfl = 0.35;
    cfg.cfl_max = 0.5;
    cfg.lts_max_levels = 3;
    cfg.reconstructor = "PCM";
    cfg.implicit_adaptive = false;
    cfg.implicit_tolerance = 1e-10;
    cfg.implicit_report_file.clear();
    return cfg;
}

double step(MeshBlock& mesh, TimeIntegrator method, double time, double limit,
            const hd_source::SourceTerm& transport, config::SimulationConfig cfg = settings(),
            double hydro_speed = 0.0)
{
    NoHydro solver;
    solver.speed = hydro_speed;
    double next = infinity;
    return hd_time::advance_one_step(mesh, time, limit, method, hydro::SpatialScheme::Godunov,
                                     solver, {transport}, cfg, {}, cfg.cfl, next);
}

void compare(const MeshBlock& a, const MeshBlock& b, double tolerance = 1e-12)
{
    for (auto flat : a.active_to_storage)
        for (int k = 0; k < NVAR; ++k)
            CHECK(near(a.cells[flat].U[k], b.cells[flat].U[k], tolerance));
}

void global_and_single_level()
{
    for (auto methods : {std::array<TimeIntegrator,2>{TimeIntegrator::Euler,TimeIntegrator::LTSEuler},
                         std::array<TimeIntegrator,2>{TimeIntegrator::RK2,TimeIntegrator::LTSRK2}})
    {
        auto a = make_mesh(), b = a;
        const auto before = integral(a);
        auto cfg = settings(); cfg.lts_max_levels = 1;
        hd_source::SourceTerm transport{"diffusion",hd_source::SourceKind::Transport,{},false,{}}; transport.transport_flux = diffusion(0.1);
        const double da = step(a, methods[0], 0.0, 1.0, transport, cfg);
        const double db = step(b, methods[1], 0.0, 1.0, transport, cfg);
        CHECK(near(da, db, 1e-14));
        CHECK(near(da, cfg.cfl/(2*0.1*24*24), 1e-14));
        compare(a,b);
        CHECK(near(integral(a)[ENERGY], before[ENERGY], 2e-14));
        CHECK(a.cells[a.active_to_storage[0]].U[ENERGY] < make_mesh().cells[a.active_to_storage[0]].U[ENERGY]);
    }
}

void multilevel_conservation(bool spherical = false)
{
    for (auto method : {TimeIntegrator::LTSEuler, TimeIntegrator::LTSRK2})
    {
        auto mesh = make_mesh(24, true, spherical);
        const auto before = integral(mesh);
        std::vector<int> visits(mesh.face_info.size());
        const auto physical = diffusion(0.1);
        hd_source::SourceTerm transport{"diffusion",hd_source::SourceKind::Transport,{},false,{}};
        transport.transport_flux = [&](const MeshBlock& stage, std::size_t face, double t,
                                   const config::SimulationConfig& cfg, Conserved& f,
                                   std::array<double,2>& rates) {
            ++visits[face];
            physical(stage,face,t,cfg,f,rates);
        };
        double t = 0.0;
        for (int n = 0; n < 5; ++n) t += step(mesh,method,t,1.0,transport);
        const auto after = integral(mesh);
        CHECK(near(after[ENERGY],before[ENERGY],2e-13));
        CHECK(near(after[RHO],before[RHO],2e-14));
        int least = 1000000, most = 0;
        for (std::size_t f = 0; f < visits.size(); ++f)
            if (mesh.face_info[f].direction == 0) {
                least = std::min(least,visits[f]); most = std::max(most,visits[f]);
            }
        CHECK(most > least); // Exercise actual coarse/fine time interfaces.
        for (auto flat : mesh.active_to_storage) CHECK(mesh.cells[flat].U[ENERGY] > 1.8);
    }
}

void eigenmode_and_implicit()
{
    // Exact decay for the semi-discrete centered diffusion operator.
    constexpr double D = 0.1, end = 0.01;
    auto mesh = make_mesh();
    const auto old = mesh;
    hd_source::SourceTerm transport{"diffusion",hd_source::SourceKind::Transport,{},false,{}}; transport.transport_flux = diffusion(D);
    double t = 0;
    while (t < end) t += step(mesh,TimeIntegrator::RK2,t,std::min(0.0001,end-t),transport);
    const double lambda = 4*D*24*24*std::pow(std::sin(pi/24),2);
    for (auto flat : mesh.active_to_storage)
        CHECK(near(mesh.cells[flat].U[ENERGY],
                   2+(old.cells[flat].U[ENERGY]-2)*std::exp(-lambda*end),1e-9));

    // BE must include transport in the residual and may exceed the explicit limit.
    mesh = old;
    auto cfg = settings(); cfg.dt_implicit = 0.02;
    hd_implicit::Integrator implicit(mesh,cfg);
    NoHydro solver;
    const double dt = implicit.advance(mesh,0.0,cfg.dt_implicit,hydro::SpatialScheme::Godunov,
                                       solver,{transport},cfg,{});
    CHECK(near(dt,cfg.dt_implicit));
    CHECK(dt > cfg.cfl/(2*D*24*24));
    for (auto flat : mesh.active_to_storage)
        CHECK(near(mesh.cells[flat].U[ENERGY],
                   2+(old.cells[flat].U[ENERGY]-2)/(1+lambda*dt),2e-9));
    CHECK(near(integral(mesh)[ENERGY],integral(old)[ENERGY],2e-10));
}

void time_limits_and_rollback()
{
    for (auto method : {TimeIntegrator::Euler, TimeIntegrator::RK2,
                        TimeIntegrator::LTSEuler, TimeIntegrator::LTSRK2})
    {
        auto mesh = make_mesh(8);
        auto cfg = settings(); cfg.lts_max_levels = 1;
        hd_source::SourceTerm transport{"diffusion",hd_source::SourceKind::Transport,{},false,{}}; transport.transport_flux = diffusion(0.1);
        const double dt = step(mesh,method,0.0,1.0,transport,cfg,1.0);
        const double hydro_rate = 8.0 + (NDIM-1)*2.0;
        CHECK(near(dt,cfg.cfl/(hydro_rate+2*0.1*8*8),1e-14));
        const auto before = mesh;
        transport.transport_flux = [](const MeshBlock&,std::size_t,double,
                                  const config::SimulationConfig&,Conserved&,std::array<double,2>& r) {
            r[0] = -1.0;
        };
        expect_throw<std::invalid_argument>([&] { step(mesh,method,dt,1.0,transport,cfg); });
        compare(mesh,before,0.0);
    }
    // Stage-dependent bounds must be checked again and rejected caches discarded.
    for (auto method : {TimeIntegrator::RK2,TimeIntegrator::LTSRK2})
    {
        auto mesh = make_mesh(8);
        const auto physical = diffusion(0.1);
        hd_source::SourceTerm transport{"diffusion",hd_source::SourceKind::Transport,{},false,{}};
        transport.transport_flux = [=](const MeshBlock& m,std::size_t f,double t,
                                   const config::SimulationConfig& c,Conserved& F,std::array<double,2>& r) {
            physical(m,f,t,c,F,r);
            if (t > 0.0) { F *= 16.0; r[0] *= 16.0; r[1] *= 16.0; }
        };
        auto cfg = settings(); cfg.lts_max_levels = 1;
        const double dt = step(mesh,method,0.0,1.0,transport,cfg);
        CHECK(dt <= cfg.cfl_max/(16*2*0.1*8*8));
        CHECK(near(integral(mesh)[ENERGY],2.0,1e-13));
    }
}

void boundary_energy_loss()
{
    // Prescribed outward flux is an exchange with the exterior, not a volume source.
    auto mesh = make_mesh(12);
    mesh.bc_lower[0] = mesh.bc_upper[0] = boundary::BoundaryCondition::Reflective;
    mesh.build_topology();
    boundary::FillBoundary(mesh,0.0,{});
    const auto initial = mesh;
    hd_source::SourceTerm transport{"diffusion",hd_source::SourceKind::Transport,{},false,{}};
    transport.transport_flux = [](const MeshBlock& m,std::size_t face,double t,
                              const config::SimulationConfig&,Conserved& F,std::array<double,2>&) {
        const auto& f = m.face_info[face];
        if (f.direction == 0 && f.right_active == MeshBlock::invalid_active_index)
            F[ENERGY] = F[INTERNAL_ENERGY] = 0.1*(1+t);
    };
    for (auto method : {TimeIntegrator::RK2,TimeIntegrator::LTSRK2})
    {
        mesh = initial;
        const double dt = step(mesh,method,0.0,0.01,transport);
        CHECK(near(integral(mesh)[ENERGY]-integral(initial)[ENERGY],-0.1*(dt+dt*dt/2),1e-13));
    }
}

void composition_conservation()
{
    if constexpr (NCOMP == 2)
    {
        for (auto method : {TimeIntegrator::RK2,TimeIntegrator::LTSRK2})
        {
            auto mesh = make_mesh(24,true);
            for (auto flat : mesh.active_to_storage)
            {
                const double x = mesh.cell_center(0,mesh.coordinates(flat))/mesh.active_length(0);
                mesh.cells[flat].U[COMP(0)] = 0.5+0.2*std::cos(2*pi*x);
                mesh.cells[flat].U[COMP(1)] = 1-mesh.cells[flat].U[COMP(0)];
            }
            boundary::FillBoundary(mesh,0.0,{});
            hd_microphysics::refresh_thermodynamics(mesh);
            const auto before = integral(mesh);
            hd_source::SourceTerm transport{"diffusion",hd_source::SourceKind::Transport,{},false,{}}; transport.transport_flux = diffusion(0.1,true);
            step(mesh,method,0.0,1.0,transport);
            const auto after = integral(mesh);
            for (int k : {COMP(0),COMP(1),ENERGY}) CHECK(near(before[k],after[k],2e-13));
            for (auto flat : mesh.active_to_storage)
                hd_microphysics::validate_cell_composition(mesh.cells[flat]);
        }
    }
}

void register_transport(const config::SimulationConfig&, const MeshBlock&,
                        hd_source::SourceList& sources)
{
    for (double D : {0.03,0.07})
    {
        hd_source::SourceTerm term{"diffusion",hd_source::SourceKind::Transport,{},false,{}};
        term.transport_flux = diffusion(D);
        // The type tag must prevent this local callback from running.
        term.func = [](const MeshBlock&,double,const config::SimulationConfig&,
                       std::vector<Conserved>&,double&) { throw std::runtime_error("wrong source path"); };
        sources.push_back(term);
    }
}

void source_kind_dispatch()
{
    const auto mesh = make_mesh(8);
    auto cfg = settings();
    cfg.constant_gravity_enabled = cfg.constant_heating_enabled = true;
    cfg.gravity[0] = 0.2;
    cfg.heating_rate = 0.01;
    hd_user::UserHooks hooks; hooks.register_sources = register_transport;
    auto sources = hd_source::make_sources(cfg,mesh,hooks);
    CHECK(sources.size() == 4);
    // A normal source's face callback must not be used.
    sources[0].transport_flux = [](const MeshBlock&,std::size_t,double,const config::SimulationConfig&,
                                   Conserved&,std::array<double,2>&) { throw std::runtime_error("wrong face path"); };
    NoHydro solver;
    double rate,limit;
    const auto rhs = hd_time::evaluate_RHS(mesh,0.0,hydro::SpatialScheme::Godunov,sources,
                                          solver,cfg,cfg.cfl,rate,limit).second;
    CHECK(near(rate,2*0.1*8*8,1e-14));
    Conserved total{};
    for (auto flat : mesh.active_to_storage)
    {
        CHECK(near(rhs[flat][MOMDEN(0)],0.2));
        total += rhs[flat]*mesh.volume_cell[flat];
    }
    CHECK(near(total[ENERGY],0.01,1e-13));
    CHECK(near(total[INTERNAL_ENERGY],0.01,1e-13));
    for (auto method : {TimeIntegrator::Euler,TimeIntegrator::RK2,
                        TimeIntegrator::LTSEuler,TimeIntegrator::LTSRK2})
    {
        auto trial = mesh;
        double next = infinity;
        const double dt = hd_time::advance_one_step(trial,0.0,0.001,method,
            hydro::SpatialScheme::Godunov,solver,sources,cfg,{},cfg.cfl,next);
        CHECK(near(integral(trial)[MOMDEN(0)],0.2*dt,1e-13));
    }
    sources.back().transport_flux = {};
    expect_throw<std::invalid_argument>([&] {
        hd_time::evaluate_RHS(mesh,0.0,hydro::SpatialScheme::Godunov,sources,
                             solver,cfg,cfg.cfl,rate,limit);
    });
}
}

int main(int argc,char** argv)
{
    return run_test("face transport: diffusion, LTS conservation, limits, implicit and boundaries",[&] {
        hd_microphysics::eos_type = hd_microphysics::Equation_of_State::IdealGas_Default;
        if constexpr (NCOMP == 2)
        {
            CHECK(argc > 1);
            const std::string path = std::string(argv[1])+"-species.txt";
            { std::ofstream file(path); file << "h1 0 0\nhe4 0 0\n"; }
            auto cfg = settings(); cfg.composition_activate_file = path;
            hd_microphysics::initialize_species_table(cfg);
            hd_microphysics::eos_type = hd_microphysics::Equation_of_State::IdealGas_Mixed;
        }
        global_and_single_level();
        multilevel_conservation();
        if constexpr (NDIM == 1) multilevel_conservation(true);
        eigenmode_and_implicit();
        time_limits_and_rollback();
        boundary_energy_loss();
        composition_conservation();
        source_kind_dispatch();
    });
}
