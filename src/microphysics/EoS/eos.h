#pragma once

namespace hd_microphysics
{

    enum struct EoS
    {
        IdealGas_Default,       //Default HD, gamma = 5/3, mean_molecular_weight = 1
        IdealGas_Single,        //Default HD, gamma and mean_molecular_weight can be modified but fixed and same value for every cell
        IdealGas_Mixed,         //Different kinds of ideal gas, with no electron/photon pressure
        Ionized,
        Ionized_Radiation
    };
}
