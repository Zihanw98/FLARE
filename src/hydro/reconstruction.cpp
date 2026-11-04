#include "reconstruction.h"
#include "utils.h"
#include <stdexcept>
namespace hydro
{
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
                return Reconstruct_PCM(mesh,direction,right_index,limiter);
            case Reconstructor::PLM:
                return Reconstruct_PLM(mesh,direction,right_index,limiter);
            default:
                throw std::runtime_error("Invalid Reconstrcutor");

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
        
        const hd_core::Primitives& q_l = WL;
        const hd_core::Primitives& q_r = WR;
        const hd_core::Primitives& q_ll = mesh.at(left_left_index).U.to_Primtivies(gamma_LL);
        const hd_core::Primitives& q_rr = mesh.at(right_right_index).U.to_Primtivies(gamma_RR);
        const double dx_l = mesh.cell_center(direction,left_index) - mesh.cell_center(direction,left_left_index);
        const double dx_m = mesh.cell_center(direction,right_index) - mesh.cell_center(direction,left_index);
        const double dx_r = mesh.cell_center(direction,right_right_index) - mesh.cell_center(direction,right_index);
        for(int i=0 ; i<hd_core::NHYDRO; i++)
        {
            slope_left[i] = limited_slope((q_l[i]-q_ll[i])/dx_l,(q_r[i]-q_l[i])/dx_m,limiter);
            slope_right[i] = limited_slope((q_r[i]-q_l[i])/dx_m,(q_rr[i]-q_r[i])/dx_r,limiter);
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