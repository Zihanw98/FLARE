/**
 * @file reconstruction.cpp
 * @brief Hydrodynamic spatial operators.
 */
#include "reconstruction.h"
#include "utils.h"
#include <stdexcept>
namespace hydro
{
    ReconstructionCache::ReconstructionCache(const hd_core::MeshBlock& mesh,
        Reconstructor method, Limiter limiter, double beta)
        : mesh_(mesh), method_(method), limiter_(limiter), beta_(beta),
          primitives_(mesh.cells.size()), primitive_ready_(mesh.cells.size(), 0)
    {
        if (method != Reconstructor::PCM && method != Reconstructor::PLM)
            throw std::runtime_error("Invalid Reconstructor");
        std::size_t stride = 1;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            strides_[d] = stride;
            stride *= mesh.ntotal(d);
        }
        if (method == Reconstructor::PCM) return;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            if (mesh.nghost[d] < 2)
                throw std::runtime_error("PLM Reconstructor requires at least 2 ghost zone for each side");
            slopes_[d].resize(mesh.cells.size());
            slope_ready_[d].assign(mesh.cells.size(), 0);
            centers_[d].resize(mesh.ntotal(d));
            widths_[d].resize(mesh.ntotal(d));
            for (int i = 0; i < mesh.ntotal(d); ++i)
            {
                centers_[d][i] = mesh.cell_center(d, i);
                widths_[d][i] = mesh.cell_width(d, i);
            }
        }
    }

    const hd_core::Primitives& ReconstructionCache::primitive(std::size_t flat)
    {
        if (!primitive_ready_[flat])
        {
            const auto& cell = mesh_.cells[flat];
            primitives_[flat] = cell.U.to_Primtivies(cell.gamma);
            primitive_ready_[flat] = 1;
        }
        return primitives_[flat];
    }

    const hd_core::Primitives& ReconstructionCache::slope(int d, std::size_t flat, int coordinate)
    {
        if (!slope_ready_[d][flat])
        {
            const auto& q = primitive(flat);
            const auto& ql = primitive(flat - strides_[d]);
            const auto& qr = primitive(flat + strides_[d]);
            const double dxl = centers_[d][coordinate] - centers_[d][coordinate - 1];
            const double dxr = centers_[d][coordinate + 1] - centers_[d][coordinate];
            for (int k = 0; k < hd_core::NHYDRO; ++k)
                slopes_[d][flat][k] = limited_slope((q[k]-ql[k])/dxl,
                    (qr[k]-q[k])/dxr, limiter_, beta_);
            slope_ready_[d][flat] = 1;
        }
        return slopes_[d][flat];
    }

    std::pair<hd_core::Primitives, hd_core::Primitives> ReconstructionCache::face(
        int d, const hd_core::Index& right)
    {
        if (d < 0 || d >= hd_core::NDIM) throw std::out_of_range("Invalid reconstruction direction");
        const auto rf = mesh_.flat_index(right);
        if (right[d] < 1) throw std::out_of_range("Invalid reconstruction left cell");
        const auto lf = rf - strides_[d];
        auto wl = primitive(lf), wr = primitive(rf);
        if (method_ == Reconstructor::PCM) return {wl, wr};
        if (right[d] < 2 || right[d] + 1 >= mesh_.ntotal(d))
            throw std::out_of_range("Invalid PLM stencil");
        const auto& sl = slope(d, lf, right[d] - 1);
        const auto& sr = slope(d, rf, right[d]);
        for (int k = 0; k < hd_core::NHYDRO; ++k)
        {
            // Preserve the original multiply/divide order for byte compatibility.
            wl[k] += sl[k] * widths_[d][right[d] - 1] / 2;
            wr[k] -= sr[k] * widths_[d][right[d]] / 2;
        }
        try { wl.validate(mesh_.cells[lf].gamma); }
        catch (const std::invalid_argument&) { wl = primitive(lf); }
        try { wr.validate(mesh_.cells[rf].gamma); }
        catch (const std::invalid_argument&) { wr = primitive(rf); }
        return {wl, wr};
    }

    Reconstructor parse_reconstructor(const std::string& name)
    {
        if (utils::iequals(name, "PLM"))
            return Reconstructor::PLM;
        return Reconstructor::PCM;
    }
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct(const Reconstructor& reconstructor, const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter, double limiter_beta)
    {
        switch(reconstructor)
        {
            case Reconstructor::PCM:
                return Reconstruct_PCM(mesh,direction,right_index,limiter,limiter_beta);
            case Reconstructor::PLM:
                return Reconstruct_PLM(mesh,direction,right_index,limiter,limiter_beta);
            default:
                throw std::runtime_error("Invalid Reconstructor");

        }
    }
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PCM(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter,double limiter_beta)
    {
        hd_core::Index left_index = right_index;
        left_index[direction]--;
        auto WL = mesh.at(left_index).U.to_Primtivies(mesh.at(left_index).gamma);
        auto WR = mesh.at(right_index).U.to_Primtivies(mesh.at(right_index).gamma);
        return {WL,WR};

        
    }
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PLM(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter, double limiter_beta)
    {
        if(mesh.nghost[direction]<2)
            throw std::runtime_error("PLM Reconstructor requires at least 2 ghost zone for each side");
        hd_core::Index left_index = right_index;
        left_index[direction]--;
        hd_core::Index left_left_index = left_index;
        left_left_index[direction]--;
        
        hd_core::Index right_right_index = right_index;
        right_right_index[direction]++;
        
        auto gamma_L = mesh.at(left_index).gamma;
        auto gamma_R = mesh.at(right_index).gamma;
        auto gamma_LL = mesh.at(left_left_index).gamma;
        auto gamma_RR = mesh.at(right_right_index).gamma;
        auto WL = mesh.at(left_index).U.to_Primtivies(gamma_L);
        auto WR = mesh.at(right_index).U.to_Primtivies(gamma_R);

        hd_core::Primitives slope_left;
        hd_core::Primitives slope_right;
        
        const hd_core::Primitives q_l = WL;
        const hd_core::Primitives q_r = WR;
        const hd_core::Primitives q_ll = mesh.at(left_left_index).U.to_Primtivies(gamma_LL);
        const hd_core::Primitives q_rr = mesh.at(right_right_index).U.to_Primtivies(gamma_RR);
        const double dx_l = mesh.cell_center(direction,left_index) - mesh.cell_center(direction,left_left_index);
        const double dx_m = mesh.cell_center(direction,right_index) - mesh.cell_center(direction,left_index);
        const double dx_r = mesh.cell_center(direction,right_right_index) - mesh.cell_center(direction,right_index);
        for(int i=0 ; i<hd_core::NHYDRO; i++)
        {
            slope_left[i] = limited_slope((q_l[i]-q_ll[i])/dx_l,(q_r[i]-q_l[i])/dx_m,limiter,limiter_beta);
            slope_right[i] = limited_slope((q_r[i]-q_l[i])/dx_m,(q_rr[i]-q_r[i])/dx_r,limiter,limiter_beta);
            WL[i] += slope_left[i] * mesh.cell_width(direction,left_index)/2;
            WR[i] -= slope_right[i] * mesh.cell_width(direction,right_index)/2;
        }
        try
        {
            WL.validate(gamma_L);
        }
        catch (const std::invalid_argument&)
        {
            WL = q_l;
        }

        try
        {
            WR.validate(gamma_R);
        }
        catch (const std::invalid_argument&)
        {
            WR = q_r;
        }
        return {WL,WR};
    }

}
