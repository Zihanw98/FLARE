#include "test_support.h"
#include "time/integrator.h"
#include "time/implicit.h"
#include <limits>
#include <iomanip>
#include "hydro/riemann/exact.h"
#include "user/user_api.h"

double expansion_pressure = 1.e-6;
hd_core::Cell expansion_state(const hd_core::MeshBlock& m, std::size_t k, double t)
{
        using namespace hd_core;
        const int dimension = m.geometry == Geometry::Spherical ? 3 : NDIM;
        const double gamma = 5.0/3.0;
        Primitives w;
        w[RHO] = std::pow(1+t, -dimension);
        w[PRE] = expansion_pressure*std::pow(1+t, -dimension*gamma);
        const auto index = m.coordinates(k);
        for (int d=0; d<NDIM; ++d) w[VEL(d)] = m.cell_center(d,index)/(1+t);
        return Cell(w,gamma);
}


// Optional physics regressions, kept in the existing integrator test executable.
// Usage: test_integrator --physics sedov|expansion n cartesian|spherical method [p0] [grid]
int dual_energy_physics(int argc, char** argv)
{
    using namespace hd_core;
    const bool blast = std::string(argv[2]) == "sedov";
    const int n = std::stoi(argv[3]);
    const bool spherical = std::string(argv[4]) == "spherical";
    const std::string method = argv[5];
    const double p0 = argc > 6 ? std::stod(argv[6]) : 1.e-6;
    expansion_pressure = p0;
    const bool narrow = argc > 7 && std::string(argv[7]) == "narrow";
    const double gamma = 5.0/3.0, end = blast ? 0.02 : 0.5;
    const int dimension = spherical ? 3 : NDIM;
    MeshBlock mesh;
    mesh.geometry = spherical ? Geometry::Spherical : Geometry::Cartesian;
    for (int d = 0; d < NDIM; ++d)
    {
        mesh.nactive[d] = n; mesh.nghost[d] = 2;
        mesh.set_linear_uniform_faces(d, blast ? (spherical ? 0.0 : -0.5) : 1.0,
                                         blast ? (spherical ? 1.0 : 0.5) : 2.0);
        mesh.bc_lower[d] = mesh.bc_upper[d] = blast ? boundary::BoundaryCondition::Reflective
                                                   : boundary::BoundaryCondition::User;
    }
    if (narrow)
    {
        CHECK(NDIM == 1 && !blast);
        auto& f = mesh.face[0];
        const double width = 1.0/(n - 1 + 1.0/64);
        f[2] = 1;
        for (int i = 0; i < n; ++i) f[i+3] = f[i+2] + width*(i == n/2 ? 1.0/64 : 1.0);
        f[1] = f[2]-width; f[0] = f[1]-width;
        f[n+3] = f[n+2]+width; f[n+4] = f[n+3]+width;
    }
    mesh.allocate(); mesh.build_geometry();
    auto radius = [&](std::size_t k) {
        const auto index = mesh.coordinates(k);
        double r2 = 0;
        for (int d = 0; d < NDIM; ++d) r2 += std::pow(mesh.cell_center(d,index),2);
        return std::sqrt(r2);
    };
    const double injection_radius = spherical ? 2.0/n : std::max(0.04, 1.0/n);
    double injection_volume = 0;
    for (std::size_t k=0; k<mesh.cells.size(); ++k)
        if (mesh.is_active(k) && radius(k)<injection_radius) injection_volume += mesh.volume_cell[k];
    if (blast) CHECK(injection_volume > 0);
    for (std::size_t k=0; k<mesh.cells.size(); ++k)
    {
        if (!blast) mesh.cells[k] = expansion_state(mesh,k,0);
        else
        {
            Primitives w; w[RHO]=1; w[PRE]=p0;
            if (mesh.is_active(k) && radius(k)<injection_radius) w[PRE] += (gamma-1)/injection_volume;
            mesh.cells[k] = Cell(w,gamma);
        }
    }
    config::SimulationConfig cfg;
    cfg.cfl = 0.4; cfg.cfl_max=0.5; cfg.lts_max_levels=8;
    cfg.reconstructor="PLM"; cfg.limiter="MC";
    hd_user::UserHooks hooks;
    if (!blast) hooks.userbc = [](MeshBlock& m,double t,int,boundary::BoundarySide) {
        for (std::size_t k=0; k<m.cells.size(); ++k) if (!m.is_active(k)) m.cells[k]=expansion_state(m,k,t);
    };
    hd_source::SourceList sources;
    if (spherical) sources.push_back(hd_source::make_geometric_source());
    hydro::ExactSolver solver(1.e-10,200);
    hydro::StarPressureCache hints;
    hd_time::LTSEulerWorkspace workspace;
    const auto integrator=hd_time::parse_time_integrator(method);
    std::unique_ptr<hd_implicit::Integrator> implicit;
    if (method == "implicit")
    {
        cfg.dt_implicit=0.005; cfg.implicit_tolerance=1.e-8;
        implicit=std::make_unique<hd_implicit::Integrator>(mesh,cfg);
    }
    auto integral = [&](int v) {
        long double sum=0;
        for (std::size_t k=0; k<mesh.cells.size(); ++k)
            if (mesh.is_active(k)) sum+=mesh.cells[k].U[v]*mesh.volume_cell[k];
        return double(sum);
    };
    const double mass0=integral(RHO), energy0=integral(ENERGY);
    double t=0, next=std::numeric_limits<double>::infinity(); int steps=0;
    boundary::FillBoundary(mesh,0,hooks);
    while (t<end && steps<200000)
    {
        t += hd_time::advance_one_step(mesh,t,end-t,integrator,hydro::SpatialScheme::Godunov,
            solver,sources,cfg,hooks,cfg.cfl,next,&workspace,&hints,implicit.get());
        ++steps;
    }
    CHECK(near(t,end,1.e-13));
    double error=0, volume=0, correction=0, pmin=1.e300, shock_radius=0;
    for (std::size_t k=0; k<mesh.cells.size(); ++k) if(mesh.is_active(k))
    {
        const auto& cell=mesh.cells[k]; cell.U.validate(cell.gamma);
        CHECK(cell.U[INTERNAL_ENERGY]>0);
        const auto w=cell.U.to_Primtivies(cell.gamma);
        pmin=std::min(pmin,w[PRE]); correction+=cell.dual_energy_correction*mesh.volume_cell[k];
        if (blast && w[RHO]>1.1) shock_radius=std::max(shock_radius,radius(k));
        bool interior=true;
        const auto index=mesh.coordinates(k);
        for (int d=0; d<NDIM; ++d)
            interior &= mesh.cell_center(d,index)>=1.3 && mesh.cell_center(d,index)<=1.9;
        if (!blast && interior)
        {
            error+=std::abs(w[PRE]/(p0*std::pow(1+t,-dimension*gamma))-1)*mesh.volume_cell[k];
            volume+=mesh.volume_cell[k];
        }
    }
    const double energy_error=(integral(ENERGY)-energy0)/energy0;
    if (blast)
    {
        CHECK(shock_radius>injection_radius);
        CHECK(std::abs(integral(RHO)/mass0-1)<1.e-10);
        CHECK(std::abs((integral(ENERGY)-energy0-correction)/energy0)<1.e-10);
        CHECK(std::abs(energy_error)<0.01);
    }
    else CHECK(error/volume<0.15); // Coarse Euler is first order; RK2 convergence is reported separately.
    std::cout<<std::setprecision(12)<<"physics,"<<(blast?"sedov":"expansion")<<','<<NDIM<<','
        <<(spherical?"spherical":"cartesian")<<','<<n<<','<<method<<','<<p0<<','
        <<(narrow?"narrow":"uniform")<<','<<t<<','<<steps<<','<<(blast?0:error/volume)<<','
        <<pmin<<','<<shock_radius<<','<<energy_error<<','<<correction/energy0<<'\n';
    return 0;
}

int main(int argc, char** argv)
{
    if (argc >= 6 && std::string(argv[1]) == "--physics")
        return run_test("dual energy physics", [&] { dual_energy_physics(argc, argv); });
    return run_test("integrator active update and rollback", [] {
        hd_core::MeshBlock mesh{};
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            mesh.nactive[d] = 3 + d;
            mesh.nghost[d] = 1 + d;
            mesh.set_linear_uniform_faces(d, 0.0, 1.0);
            mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
        }
        mesh.allocate();
        mesh.build_geometry();
        std::vector<bool> active(mesh.cells.size(), true);
        std::vector<hd_core::Conserved> rhs(mesh.cells.size());
        std::size_t last_active = 0;
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            mesh.cells[flat].U[hd_core::RHO] = 1.0;
            mesh.cells[flat].U[hd_core::ENERGY] = 10.0;
            mesh.cells[flat].gamma = flat % 2 ? 1.4 : 5.0 / 3.0;
            const auto index = mesh.coordinates(flat);
            for (int d = 0; d < hd_core::NDIM; ++d)
                active[flat] = active[flat] && index[d] >= mesh.nghost[d]
                    && index[d] < mesh.nghost[d] + mesh.nactive[d];
            rhs[flat].data.fill(active[flat] ? 0.25
                : std::numeric_limits<double>::quiet_NaN());
            if (active[flat]) last_active = flat;
        }
        const auto original = mesh.cells;
        auto trial = mesh;
        hd_time::add_scaled_rhs(trial, rhs, 0.5);
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            CHECK(mesh.cells[flat].U.data == original[flat].U.data);
            for (int k = 0; k < hd_core::NVAR; ++k)
                CHECK(trial.cells[flat].U[k] == original[flat].U[k]
                      + (active[flat] ? 0.125 : 0.0));
            CHECK(trial.cells[flat].gamma == original[flat].gamma);
        }
        CHECK(trial.face == mesh.face);
        CHECK(trial.bc_lower == mesh.bc_lower);
        CHECK(trial.bc_upper == mesh.bc_upper);

        // Fail at the last active cell: preceding candidate updates must not leak.
        auto bad_rhs = rhs;
        bad_rhs[last_active][hd_core::ENERGY] = -1e6;
        auto check_unchanged = [&] {
            for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
            {
                CHECK(mesh.cells[flat].U.data == original[flat].U.data);
                CHECK(mesh.cells[flat].gamma == original[flat].gamma);
            }
        };
        expect_throw<std::invalid_argument>([&] {
            hd_time::add_scaled_rhs(mesh, bad_rhs, 0.5);
        });
        check_unchanged();
        bad_rhs[last_active][hd_core::ENERGY] = std::numeric_limits<double>::quiet_NaN();
        expect_throw<std::invalid_argument>([&] {
            hd_time::add_scaled_rhs(mesh, bad_rhs, 0.5);
        });
        check_unchanged();
        expect_throw<std::invalid_argument>([&] { hd_time::add_scaled_rhs(mesh, {}, 0.5); });
        check_unchanged();
        expect_throw<std::invalid_argument>([&] {
            hd_time::add_scaled_rhs(mesh, rhs, std::numeric_limits<double>::infinity());
        });
        check_unchanged();
        hd_time::add_scaled_rhs(mesh, bad_rhs, 0.0); // Zero scale is a no-op.
        check_unchanged();
        hd_time::add_scaled_rhs(trial, rhs, -0.5);
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
            CHECK(trial.cells[flat].U.data == original[flat].U.data);
        mesh.cells.clear();
        expect_throw<std::logic_error>([&] { hd_time::add_scaled_rhs(mesh, {}, 0.5); });
    });
}
