#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EOS_CompOSE_Startup(CCTK_ARGUMENTS)

    use EOS_CompOSE_Module
    use EOS_CompOSE_InterpCoeffs_Module
    implicit none

    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_ARGUMENTS

    if(poly_gamma_initial .gt. 0d0) then
        poly_gamma_ini = poly_gamma_initial
    else
        poly_gamma_ini = poly_gamma
    end if

    poly_k_cgs = poly_k * rho_gf**poly_gamma_ini / press_gf
    
    call EOS_CompOSE_SetInterpCoeffs(EOS_CompOSE_interpCoeffs)

end subroutine EOS_CompOSE_Startup
