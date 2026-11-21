/**
 * @file reconstruction.cpp
 * @brief Hydrodynamic spatial operators.
 */
#include "reconstruction.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace hydro
{
    namespace
    {
        void compute_species_slopes(const hd_core::Primitives& left, const hd_core::Primitives& center,
                                    const hd_core::Primitives& right, double left_spacing, double right_spacing,
                                    Limiter limiter, double beta, hd_core::Primitives& slope)
        {
            for (int k = 0; k < hd_core::NCOMP; ++k)
            {
                const int j = hd_core::COMP(k);
                slope[j] = limited_slope((center[j] - left[j]) / left_spacing,
                                         (right[j] - center[j]) / right_spacing, limiter, beta);
            }
        }

        void reconstruct_species_face(const hd_core::Primitives& center, const hd_core::Primitives& slope,
                                      double offset, hd_core::Primitives& face)
        {
            if constexpr (hd_core::NCOMP == 0) return;
            if (center[hd_core::RHO] == 0.0)
            {
                // A vacuum cell has no composition to supply to a gas face.
                face = center;
                return;
            }

            double sum = 0.0;
            for (int k = 0; k < hd_core::NCOMP; ++k)
            {
                const int j = hd_core::COMP(k);
                face[j] = std::clamp(center[j] + offset * slope[j], 0.0, 1.0);
                sum += face[j];
            }
            // Independent limiters can change the sum. Normalize this side of the face
            // so its species fluxes add up to the Riemann solver's mass flux.
            if (sum == 1.0) return;
            for (int k = 0; k < hd_core::NCOMP; ++k)
            {
                const int j = hd_core::COMP(k);
                face[j] = sum > 0.0 ? face[j] / sum : center[j];
            }
        }

        // Limit both extrapolations W +/- slope*width/2 with the same factor.
        // Apply after the user's slope limiter; leave unconstrained slopes intact.
        void correct_positive_slope(const hd_core::Primitives& center,
                                    hd_core::Primitives& slope, double width)
        {
            constexpr double retained_fraction = 1e-10;
            double theta = 1.0;
            for (int k : {hd_core::RHO, hd_core::PRE})
            {
                const double change = std::abs(slope[k] * width / 2);
                const double allowed = (1 - retained_fraction) * center[k];
                if (change > allowed)
                    theta = std::min(theta, allowed / change);
            }
            if (theta < 1.0)
                for (int k = 0; k < hd_core::INTERNAL_ENERGY; ++k)
                    slope[k] *= theta;
        }

        std::pair<hd_core::Primitives, hd_core::Primitives> reconstruct_plm(
            const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index,
            const Limiter& limiter, double limiter_beta, bool corrected);
    }

    ReconstructionCache::ReconstructionCache(const hd_core::MeshBlock& mesh,
        Reconstructor method, Limiter limiter, double beta)
    {
        reset(mesh, method, limiter, beta);
    }

    void ReconstructionCache::reset(const hd_core::MeshBlock& mesh,
        Reconstructor method, Limiter limiter, double beta)
    {
        if (method != Reconstructor::PCM && method != Reconstructor::PLM
            && method != Reconstructor::PLM_corrected)
            throw std::runtime_error("Invalid Reconstructor");
        mesh_ = &mesh;
        method_ = method;
        limiter_ = limiter;
        beta_ = beta;
        primitives_.resize(mesh.cells.size());
        primitive_epoch_.resize(mesh.cells.size(), 0);
        if (++epoch_ == 0)
        {
            std::fill(primitive_epoch_.begin(), primitive_epoch_.end(), 0);
            for (auto& stamps : slope_epoch_) std::fill(stamps.begin(), stamps.end(), 0);
            epoch_ = 1;
        }
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
            slope_epoch_[d].resize(mesh.cells.size(), 0);
            centers_[d].resize(mesh.ntotal(d));
            widths_[d].resize(mesh.ntotal(d));
            for (int i = 0; i < mesh.ntotal(d); ++i)
            {
                centers_[d][i] = mesh.cell_center(d, i);
                widths_[d][i] = mesh.cell_width(d, i);
            }
        }
    }

    const hd_core::Primitives& ReconstructionCache::primitive(std::size_t cell_flat)
    {
        if (primitive_epoch_[cell_flat] != epoch_)
        {
            const auto& cell = mesh_->cells[cell_flat];
            primitives_[cell_flat] = cell.Get_Primitives();
            primitive_epoch_[cell_flat] = epoch_;
        }
        return primitives_[cell_flat];
    }

    const hd_core::Primitives& ReconstructionCache::slope(
        int direction, std::size_t cell_flat, int axis_index)
    {
        if (slope_epoch_[direction][cell_flat] != epoch_)
        {
            const auto& center_state = primitive(cell_flat);
            const auto& left_state = primitive(cell_flat - strides_[direction]);
            const auto& right_state = primitive(cell_flat + strides_[direction]);
            const double left_spacing = centers_[direction][axis_index] - centers_[direction][axis_index - 1];
            const double right_spacing = centers_[direction][axis_index + 1] - centers_[direction][axis_index];
            // Reconstruct rho, v, p and mass fractions. Derive e/rho at the face.
            for (int k = 0; k < hd_core::INTERNAL_ENERGY; ++k)
                slopes_[direction][cell_flat][k] = limited_slope((center_state[k]-left_state[k])/left_spacing,
                    (right_state[k]-center_state[k])/right_spacing, limiter_, beta_);
            compute_species_slopes(left_state, center_state, right_state, left_spacing, right_spacing,
                                   limiter_, beta_, slopes_[direction][cell_flat]);
            if (method_ == Reconstructor::PLM_corrected)
                correct_positive_slope(center_state, slopes_[direction][cell_flat],
                                       widths_[direction][axis_index]);
            slope_epoch_[direction][cell_flat] = epoch_;
        }
        return slopes_[direction][cell_flat];
    }

    std::pair<hd_core::Primitives, hd_core::Primitives> ReconstructionCache::face(
        int direction, const hd_core::Index& right_index)
    {
        if (direction < 0 || direction >= hd_core::NDIM)
            throw std::out_of_range("Invalid reconstruction direction");
        const auto right_flat = mesh_->flat_index(right_index);
        if (right_index[direction] < 1)
            throw std::out_of_range("Invalid reconstruction left cell");
        const auto left_flat = right_flat - strides_[direction];
        // Face states are copies: extrapolation must not change cached cell centers.
        const auto& left_center = primitive(left_flat);
        const auto& right_center = primitive(right_flat);
        auto left_state = left_center;
        auto right_state = right_center;
        if (method_ == Reconstructor::PCM) return {left_state, right_state};
        if (right_index[direction] < 2 || right_index[direction] + 1 >= mesh_->ntotal(direction))
            throw std::out_of_range("Invalid PLM stencil");
        const auto& left_slope = slope(direction, left_flat, right_index[direction] - 1);
        const auto& right_slope = slope(direction, right_flat, right_index[direction]);
        const double left_width = widths_[direction][right_index[direction] - 1];
        const double right_width = widths_[direction][right_index[direction]];
        for (int k = 0; k < hd_core::INTERNAL_ENERGY; ++k)
        {
            // Preserve the original multiply/divide order for byte compatibility.
            left_state[k] += left_slope[k] * left_width / 2;
            right_state[k] -= right_slope[k] * right_width / 2;
        }
        reconstruct_species_face(left_center, left_slope, left_width / 2, left_state);
        reconstruct_species_face(right_center, right_slope, -right_width / 2, right_state);
        try
        {
            left_state.update_internal_energy(mesh_->cells[left_flat].gamma);
            left_state.validate(mesh_->cells[left_flat].gamma);
        }
        catch (const std::invalid_argument&) { left_state = left_center; }
        try
        {
            right_state.update_internal_energy(mesh_->cells[right_flat].gamma);
            right_state.validate(mesh_->cells[right_flat].gamma);
        }
        catch (const std::invalid_argument&) { right_state = right_center; }
        return {left_state, right_state};
    }

    Reconstructor parse_reconstructor(const std::string& name)
    {
        if (utils::iequals(name, "PLM"))
            return Reconstructor::PLM;
        if (utils::iequals(name, "PLM_corrected"))
            return Reconstructor::PLM_corrected;
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
            case Reconstructor::PLM_corrected:
                return Reconstruct_PLM_corrected(mesh,direction,right_index,limiter,limiter_beta);
            default:
                throw std::runtime_error("Invalid Reconstructor");

        }
    }
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PCM(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& /*limiter*/,double /*limiter_beta*/)
    {
        hd_core::Index left_index = right_index;
        left_index[direction]--;
        auto WL = mesh.at(left_index).Get_Primitives();
        auto WR = mesh.at(right_index).Get_Primitives();
        return {WL,WR};

        
    }
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PLM(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter, double limiter_beta)
    {
        return reconstruct_plm(mesh, direction, right_index, limiter, limiter_beta, false);
    }

    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PLM_corrected(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter, double limiter_beta)
    {
        return reconstruct_plm(mesh, direction, right_index, limiter, limiter_beta, true);
    }

    namespace
    {
    std::pair<hd_core::Primitives,hd_core::Primitives> reconstruct_plm(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter, double limiter_beta, bool corrected)
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
        auto WL = mesh.at(left_index).Get_Primitives();
        auto WR = mesh.at(right_index).Get_Primitives();

        hd_core::Primitives slope_left;
        hd_core::Primitives slope_right;
        
        const hd_core::Primitives q_l = WL;
        const hd_core::Primitives q_r = WR;
        const hd_core::Primitives q_ll = mesh.at(left_left_index).Get_Primitives();
        const hd_core::Primitives q_rr = mesh.at(right_right_index).Get_Primitives();
        const double dx_l = mesh.cell_center(direction,left_index) - mesh.cell_center(direction,left_left_index);
        const double dx_m = mesh.cell_center(direction,right_index) - mesh.cell_center(direction,left_index);
        const double dx_r = mesh.cell_center(direction,right_right_index) - mesh.cell_center(direction,right_index);
        for(int i=0 ; i<hd_core::INTERNAL_ENERGY; i++)
        {
            slope_left[i] = limited_slope((q_l[i]-q_ll[i])/dx_l,(q_r[i]-q_l[i])/dx_m,limiter,limiter_beta);
            slope_right[i] = limited_slope((q_r[i]-q_l[i])/dx_m,(q_rr[i]-q_r[i])/dx_r,limiter,limiter_beta);
        }
        compute_species_slopes(q_ll, q_l, q_r, dx_l, dx_m, limiter, limiter_beta, slope_left);
        compute_species_slopes(q_l, q_r, q_rr, dx_m, dx_r, limiter, limiter_beta, slope_right);
        const double left_width = mesh.cell_width(direction, left_index);
        const double right_width = mesh.cell_width(direction, right_index);
        if (corrected)
        {
            correct_positive_slope(q_l, slope_left, left_width);
            correct_positive_slope(q_r, slope_right, right_width);
        }
        for(int i=0 ; i<hd_core::INTERNAL_ENERGY; i++)
        {
            WL[i] += slope_left[i] * left_width/2;
            WR[i] -= slope_right[i] * right_width/2;
        }
        reconstruct_species_face(q_l, slope_left, left_width / 2, WL);
        reconstruct_species_face(q_r, slope_right, -right_width / 2, WR);
        try
        {
            WL.update_internal_energy(gamma_L);
            WL.validate(gamma_L);
        }
        catch (const std::invalid_argument&)
        {
            WL = q_l;
        }

        try
        {
            WR.update_internal_energy(gamma_R);
            WR.validate(gamma_R);
        }
        catch (const std::invalid_argument&)
        {
            WR = q_r;
        }
        return {WL,WR};
    }
    } // namespace

}
