module EOS_CompOSE_Module
    implicit none

    !Physical constants
    real(kind=8), parameter :: const_Msol   = 1.98847d30       !kg
    real(kind=8), parameter :: const_c      = 2.99792458d8     !m * s-1
    real(kind=8), parameter :: const_G      = 6.67408d-11      !m3 * kg-1 * s-2
    real(kind=8), parameter :: const_Mn_kg  = 1.6726d-27       !kg
    real(kind=8), parameter :: const_fm     = 1.0d-15          !m
    real(kind=8), parameter :: const_qe     = 1.6021796208d-19 !C
    real(kind=8), parameter :: const_Mn_Mev = 9.3828d2         !MeV/c**2

    !Unit conversion factors calculated by hand
    real(kind=8), parameter :: unit_rho_tabToCode     = 2.7081965813424477d-03
    real(kind=8), parameter :: unit_press_tabToCode   = 2.8864091366088927d-06
    real(kind=8), parameter :: unit_dP_dnb_tabToCode  = 1.0658048815563104d-03
    real(kind=8), parameter :: unit_dP_deps_tabToCode = 2.7082599646973916d-03

    integer, dimension(64,64) :: EOS_CompOSE_interpCoeffs

    real(kind=8), allocatable :: EOS_CompOSE_table(:,:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_yeArray(:)
    real(kind=8), allocatable :: EOS_CompOSE_rhoArray(:)
    real(kind=8), allocatable :: EOS_CompOSE_tempArray(:)

    real(kind=8), allocatable :: EOS_CompOSE_pressTable(:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_epsTable(:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_dP_dRhoTable(:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_dP_depsTable(:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_c_s2Table(:,:,:)

    real(kind=8), allocatable :: EOS_CompOSE_pressInterp(:,:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_epsInterp(:,:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_dP_dRhoInterp(:,:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_dP_depsInterp(:,:,:,:)
    real(kind=8), allocatable :: EOS_CompOSE_c_s2Interp(:,:,:,:)

    integer :: EOS_CompOSE_yeCount
    integer :: EOS_CompOSE_rhoCount
    integer :: EOS_CompOSE_tempCount
    integer :: EOS_CompOSE_varCount

    real(kind=8) :: EOS_CompOSE_d_ye
    real(kind=8) :: EOS_CompOSE_d_rho
    real(kind=8) :: EOS_CompOSE_d_temp

    !EOS_Omni constants
    ! conversion factors between cgs and M_Sun = c = G = 1
    ! see EOS_Omni/doc/units.py
    real(kind=8), parameter :: rho_gf = 1.61887093132742d-18
    real(kind=8), parameter :: press_gf = 1.80123683248503d-39
    real(kind=8), parameter :: eps_gf = 1.11265005605362d-21
    real(kind=8), parameter :: time_gf = 2.03040204956746d05
    real(kind=8), parameter :: mass_gf =  5.02916918125126d-34
    real(kind=8), parameter :: length_gf = 6.77269222552442d-06

    ! Inverses of the numbers above, calculated manually instead of by
    ! the compiler

    real(kind=8), parameter :: inv_rho_gf = 6.17714470405638d17
    real(kind=8), parameter :: inv_press_gf = 5.55174079257738d38
    real(kind=8), parameter :: inv_eps_gf = 8.98755178736818d20
    real(kind=8), parameter :: inv_time_gf = 4.92513293223396d-6
    real(kind=8), parameter :: inv_mass_gf = 1.98840000000000d33
    real(kind=8), parameter :: inv_length_gf = 1.47651770773117d05

    real(kind=8), parameter :: clite = 2.99792458d10
    real(kind=8), parameter :: cliteinv2 = 1.11265005605362d-21

    ! These values are set by EOS_Omni_Startup
    real(kind=8) :: poly_k_cgs = 0.0d0
    real(kind=8) :: poly_gamma_ini

end module EOS_CompOSE_Module
