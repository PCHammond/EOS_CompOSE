#include <cctk.h>
#include <cctk_Arguments.h>

CCTK_INT EOS_CompOSE_GetHandle_(CCTK_STRING name)
{
    if (CCTK_EQUALS(name, "2D_Polytrope"))
        return 1;
    if (CCTK_EQUALS(name, "Ideal_Fluid"))
        return 2;
    if (CCTK_EQUALS(name, "Hybrid"))
        return 3;
    if (CCTK_EQUALS(name, "nuc_eos"))
        return 4;
    if (CCTK_EQUALS(name, "cold_tabulated"))
        return 5;
    if (CCTK_EQUALS(name, "barotropic_tabulated"))
        return 6;
    if (CCTK_EQUALS(name, "CompOSE_tabulated"))
        return 7;
    return 0;
}
