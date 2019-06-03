#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

! eoskey:
! 1 --- polytropic EOS
! 7 --- CompOSE tabulated

subroutine EOS_CompOSE_Press(eoskey,keytemp,rf_precision,npoints, &
                             rho,eps,temp,ye,press,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(out)   :: press(npoints)

    CCTK_REAL :: xpress

    ! local vars
    integer          :: i
    character(256)   :: warnstring

    anyerr    = 0
    keyerr(:) = 0

    select case (eoskey)
        case (1)
            ! polytropic EOS
            if(keytemp.eq.1) then
                do i=1,npoints
                    eps(i) = press_gf * poly_k_cgs * &
                        (rho(i)*inv_rho_gf)**(poly_gamma) / &
                        (poly_gamma - 1.0d0) / rho(i)
                end do
            end if

            do i=1,npoints
                press(i) = press_gf * poly_k_cgs * &
                    (rho(i)*inv_rho_gf)**poly_gamma
            end do

        case (7)
            ! CompOSE tabulated
            if (keytemp.eq.1) then
            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, eps(i))
            end do
            end if

            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_pressInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, xpress)
                press(i) = 10**xpress
            end do

        case DEFAULT
            write(warnstring,*) "eoskey ",eoskey," not implemented!"
            call CCTK_ERROR(warnstring)
            STOP
    end select

end subroutine EOS_CompOSE_Press

subroutine EOS_CompOSE_Press_cs2(eoskey,keytemp,rf_precision,npoints, &
                                 rho,eps,temp,ye,press,cs2,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(out)   :: press(npoints)
    CCTK_REAL, intent(out)   :: cs2(npoints)

    CCTK_REAL :: xpress

    ! local vars
    integer          :: i
    character(256)   :: warnstring

    anyerr    = 0
    keyerr(:) = 0

    select case (eoskey)
        case (1)
            ! polytropic EOS
            if(keytemp.eq.1) then
                do i=1,npoints
                    eps(i) = press_gf * poly_k_cgs * &
                        (rho(i)*inv_rho_gf)**(poly_gamma) / &
                        (poly_gamma - 1.0d0) / rho(i)
                end do
            end if

            do i=1,npoints
                press(i) = press_gf * poly_k_cgs * &
                    (rho(i)*inv_rho_gf)**poly_gamma
                cs2(i) = poly_gamma * press(i) / rho(i) / &
                    (1 + eps(i) + press(i)/rho(i))
            end do

        case (7)
            ! CompOSE tabulated
            if (keytemp.eq.1) then
            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, eps(i))
            end do
            end if

            do i=1,npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_pressInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, xpress)
                press(i) = 10**xpress

                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_c_s2Interp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, cs2(i))
                if (cs2(i).lt.0.0d0) then
                    cs2(i)=0.0d0
                end if
            end do

        case DEFAULT
            write(warnstring,*) "eoskey ",eoskey," not implemented!"
            call CCTK_ERROR(warnstring)
            STOP
    end select

end subroutine EOS_CompOSE_Press_cs2

subroutine EOS_CompOSE_PressOMP(eoskey,keytemp,rf_precision,npoints, &
                                rho,eps,temp,ye,press,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module

    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(out)   :: press(npoints)

    CCTK_REAL :: xpress

    ! local vars
    integer          :: i
    character(256)   :: warnstring

    CCTK_INT :: my_anyerr

    anyerr    = 0
    keyerr(:) = 0

    select case (eoskey)
        case (1)
            ! polytropic EOS
            if(keytemp.eq.1) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i=1,npoints
                    eps(i) = press_gf * poly_k_cgs * &
                        (rho(i)*inv_rho_gf)**(poly_gamma) / &
                        (poly_gamma - 1.0d0) / rho(i)
                end do
                !$OMP END PARALLEL DO
            end if
            !$OMP PARALLEL DO PRIVATE(i)
            do i=1,npoints
                press(i) = press_gf * poly_k_cgs * &
                    (rho(i)*inv_rho_gf)**poly_gamma
            end do
            !$OMP END PARALLEL DO

        case (7)
            ! CompOSE tabulated
            if (keytemp.eq.1) then
            !$OMP PARALLEL DO PRIVATE(i)
            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, eps(i))
            end do
            !$OMP END PARALLEL DO
            end if

            !$OMP PARALLEL DO PRIVATE(i,xpress)
            do i=1,npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_pressInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, xpress)
                press(i) = 10**xpress
            end do
            !$OMP END PARALLEL DO
        case DEFAULT
            write(warnstring,*) "eoskey ",eoskey," not implemented!"
            call CCTK_ERROR(warnstring)
            STOP
    end select

end subroutine EOS_CompOSE_PressOMP

subroutine EOS_CompOSE_DPressByDEps(eoskey,keytemp,rf_precision,npoints, &
                                    rho,eps,temp,ye,dpdepsrho,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(out)   :: dpdepsrho(npoints)

    CCTK_REAL :: xdpde

    ! local vars
    integer          :: i
    character(256)   :: warnstring

    anyerr    = 0
    keyerr(:) = 0

    select case (eoskey)
        case (1)
            ! polytropic EOS
            if(keytemp.eq.1) then
                do i=1,npoints
                    eps(i) = press_gf * poly_k_cgs * &
                        (rho(i)*inv_rho_gf)**(poly_gamma) / &
                        (poly_gamma - 1.0d0) / rho(i)
                end do
            end if
            do i=1,npoints
                dpdepsrho(i) = 0.0d0
            end do
        case (7)
            ! CompOSE tabulated
            if (keytemp.eq.1) then
            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, eps(i))
            end do
            end if

            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_dP_depsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, xdpde)
                dpdepsrho(i) = 10**xdpde
            end do

        case DEFAULT
            write(warnstring,*) "eoskey ",eoskey," not implemented!"
            call CCTK_ERROR(warnstring)
            STOP
    end select

end subroutine EOS_CompOSE_DPressByDEps

subroutine EOS_CompOSE_DPressByDRho(eoskey,keytemp,rf_precision,npoints, &
                                    rho,eps,temp,ye,dpdrhoe,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(out)   :: dpdrhoe(npoints)

    CCTK_REAL :: xdpdr

    ! local vars
    integer          :: i
    character(256)   :: warnstring

    anyerr    = 0
    keyerr(:) = 0

    select case (eoskey)
        case (1)
            ! polytropic EOS
            if(keytemp.eq.1) then
                do i=1,npoints
                    eps(i) = press_gf * poly_k_cgs * &
                        (rho(i)*inv_rho_gf)**(poly_gamma) / &
                        (poly_gamma - 1.0d0) / rho(i)
                end do
            end if
            do i=1,npoints
                dpdrhoe(i) = press_gf * poly_k_cgs *  &
                    poly_gamma * inv_rho_gf *        &
                    (rho(i)*inv_rho_gf) ** (poly_gamma - 1.d0)
            enddo

        case (7)
            ! CompOSE tabulated
            if (keytemp.eq.1) then
            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, eps(i))
            end do
            end if

            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_dP_dRhoInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, xdpdr)
                dpdrhoe(i) = 10**xdpdr
            end do

        case DEFAULT
            write(warnstring,*) "eoskey ",eoskey," not implemented!"
            call CCTK_ERROR(warnstring)
            STOP
    end select

end subroutine EOS_CompOSE_DPressByDRho

subroutine EOS_CompOSE_dpderho_dpdrhoe(eoskey,keytemp,rf_precision,npoints,&
                                  rho,eps,temp,ye,dpderho,dpdrhoe,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(out)   :: dpderho(npoints)
    CCTK_REAL, intent(out)   :: dpdrhoe(npoints)

    CCTK_REAL :: xdpde, xdpdr

    ! local vars
    integer          :: i
    character(256)   :: warnstring

    anyerr    = 0
    keyerr(:) = 0

    select case(eoskey)
        case (1)
            ! polytropic EOS
            if(keytemp.eq.1) then
                do i=1,npoints
                    eps(i) = press_gf * poly_k_cgs * &
                        (rho(i)*inv_rho_gf)**(poly_gamma) / &
                        (poly_gamma - 1.0d0) / rho(i)
                end do
            end if
            do i=1,npoints
                dpdrhoe(i) = press_gf * poly_k_cgs *  &
                    poly_gamma * inv_rho_gf *        &
                    (rho(i)*inv_rho_gf) ** (poly_gamma - 1.d0)
                dpderho(i) = 0.0d0
            end do
        case (7)
            ! CompOSE tabulated
            if (keytemp.eq.1) then
            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, eps(i))
            end do
            end if

            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_dP_depsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, xdpde)
                dpderho(i) = 10**xdpde

                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_dP_dRhoInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, xdpdr)
                dpdrhoe(i) = 10**xdpdr
            end do

        case DEFAULT
            write(warnstring,*) "eoskey ",eoskey," not implemented!"
            call CCTK_ERROR(warnstring)
            STOP
    end select

end subroutine EOS_CompOSE_dpderho_dpdrhoe

subroutine EOS_CompOSE_cs2(eoskey,keytemp,rf_precision,npoints,&
                           rho,eps,temp,ye,cs2,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(out)   :: cs2(npoints)

    ! local vars
    integer          :: i
    character(256)   :: warnstring
    real(kind=8) :: xpress

    anyerr    = 0
    keyerr(:) = 0

    select case (eoskey)
        case (1)
            ! polytropic EOS
            if(keytemp.eq.1) then
                do i=1,npoints
                    eps(i) = press_gf * poly_k_cgs * &
                        (rho(i)*inv_rho_gf)**(poly_gamma) / &
                        (poly_gamma - 1.0d0) / rho(i)
                end do
            end if
            do i=1,npoints
                xpress = press_gf*poly_k_cgs * &
                    (rho(i)*inv_rho_gf)**(poly_gamma)
                cs2(i) = poly_gamma * xpress / rho(i) / &
                    (1 + eps(i) + xpress/rho(i))
            end do
        case (7)
            ! CompOSE tabulated
            if (keytemp.eq.1) then
            do i=1, npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, eps(i))
            end do
            end if

            do i=1,npoints
                call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
                    EOS_CompOSE_c_s2Interp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, cs2(i))
                if (cs2(i).lt.0.0d0) then
                    cs2(i)=0.0d0
                end if
            end do

        case DEFAULT
            write(warnstring,*) "eoskey ",eoskey," not implemented!"
            call CCTK_ERROR(warnstring)
            STOP
    end select
end subroutine EOS_CompOSE_cs2

subroutine EOS_CompOSE_eps_from_press(eoskey,keytemp,rf_precision,npoints,&
                                      rho,eps,temp,ye,press,xeps,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints),press(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(out)   :: xeps(npoints)

    ! local vars
    integer        :: i
    character(256) :: warnstring
    CCTK_REAL :: xpress(npoints)
    real(kind=8) :: temp_current,temp_min,temp_max
    real(kind=8) :: inv_grp1 = 0.38196601125010515d0
    character(len=256) :: message

    if(keytemp.eq.1) then
        anyerr = 1
        keyerr(:) = -1
    else
        anyerr    = 0
        keyerr(:) = 0
    end if

    select case (eoskey)
        case (1)
            ! polytropic EOS
            do i=1,npoints
                xeps(i) = press(i) / (poly_gamma - 1.0d0) / rho(i)
            end do

        case (7)
            ! CompOSE tabulated
            do i=1, npoints
                call EOS_CompOSE_Press(eoskey,0,rf_precision,1, &
                    rho(i),xeps,eos_compose_temp_max,ye(i),xpress(i),keyerr(i),anyerr)
                if (press(i).gt.xpress(i)) then
                    call CCTK_ERROR("Pressure too high for EOS")
                    STOP
                end if
        
                call EOS_CompOSE_Press(eoskey,0,rf_precision,1, &
                    rho(i),xeps,eos_compose_temp_min,ye(i),xpress(i),keyerr(i),anyerr)
                if (press(i).lt.xpress(i)) then
                    !if (abs(press(i)-xpress(i)).lt.rf_precision) then
                    !    temp_current = eos_compose_temp_min
                    !    cycle
                    !end if
                    !call CCTK_INFO("rho, temp_min, ye:")
                    !write(message,*) rho(i), eos_compose_temp_min, ye(i)
                    !call CCTK_INFO(message)
                    !call CCTK_INFO("Pressure requested is:")
                    !write(message,*) press(i)
                    !call CCTK_INFO(message)
                    !call CCTK_INFO("Pressure at minimum temperature is:")
                    !write(message,*) xpress(i)
                    !call CCTK_INFO(message)
                    !call CCTK_ERROR("Pressure too low for EOS")
                    !STOP
                    temp_current = eos_compose_temp_min
                    cycle
                end if

                call EOS_CompOSE_Press(eoskey,0,rf_precision,1, &
                    rho(i),xeps,temp(i),ye(i),xpress(i),keyerr(i),anyerr)
                if (press(i).gt.xpress(i)) then
                    temp_current = temp(i)
                    temp_min = temp(i)
                    temp_max = eos_compose_temp_max
                else if (press(i).lt.xpress(i)) then
                    temp_current = temp(i)
                    temp_min = eos_compose_temp_min
                    temp_max = temp(i)
                else
                    temp_current = temp(i)
                    cycle
                end if

                do while (abs(log10(press(i))-log10(xpress(i)))/abs(log10(press(i))) .gt. rf_precision)
                    temp_current = temp_min + (temp_max-temp_min)*inv_grp1
                    call EOS_CompOSE_Press(eoskey,0,rf_precision,1, &
                        rho(i),xeps,temp_current,ye(i),xpress(i),keyerr(i),anyerr)
                    if (xpress(i).gt.press(i)) then
                        temp_max = temp_current
                    else
                        temp_min = temp_current
                    end if
                end do

                temp(i) = temp_current
                call EOS_CompOSE_Press(eoskey,1,rf_precision,1, &
                        rho(i),eps(i),temp(i),ye(i),xpress(i),keyerr(i),anyerr)
            end do
        case DEFAULT
            write(warnstring,*) "eoskey ",eoskey," not implemented!"
            call CCTK_ERROR(warnstring)
            STOP
    end select
end subroutine EOS_CompOSE_eps_from_press

subroutine EOS_CompOSE_Eps(eoskey,keytemp,rf_precision,npoints, &
    rho,eps,temp,ye,keyerr,anyerr)

    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)

    ! local vars
    integer          :: i
    character(256)   :: warnstring

    anyerr    = 0
    keyerr(:) = 0

    select case (eoskey)
    case (7)
    ! CompOSE tabulated
    do i=1, npoints
        call EOS_CompOSE_Eval(ye(i), log10(rho(i)), log10(temp(i)), &
            EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
            EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
            EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
            EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
            EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, eps(i))
    end do

    case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
    end select
end subroutine EOS_CompOSE_Eps

subroutine EOS_CompOSE_temp_from_eps(eoskey,keytemp,rf_precision,npoints, &
                                     rho,eps,temp,ye,keyerr,anyerr,epstol)
    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
    CCTK_INT, intent(out)    :: keyerr(npoints)
    CCTK_INT, intent(out)    :: anyerr
    CCTK_REAL, intent(in)    :: rf_precision
    CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
    CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
    CCTK_REAL, intent(in) :: epstol

    !local vars
    real(kind=8) :: temp_current,temp_min,temp_max,xeps
    integer :: i
    real(kind=8) :: inv_grp1 = 0.38196601125010515d0

    anyerr    = 0
    keyerr(:) = 0

    do i=1,npoints
        call EOS_CompOSE_Eps(eoskey,keytemp,rf_precision,1, &
        rho(i),xeps,eos_compose_temp_max,ye(i),keyerr(i),anyerr)
        if (eps(i).gt.xeps) then
            call CCTK_ERROR("Epsilon too high for EOS")
            STOP
        end if

        call EOS_CompOSE_Eps(eoskey,keytemp,rf_precision,1, &
        rho(i),xeps,eos_compose_temp_min,ye(i),keyerr(i),anyerr)
        if (eps(i).lt.xeps) then
            !call CCTK_ERROR("Epsilon too low for EOS")
            !STOP
            temp(i) = eos_compose_temp_min
            eps(i)  = xeps
        end if

        call EOS_CompOSE_Eps(eoskey,keytemp,rf_precision,1, &
            rho(i),xeps,temp(i),ye(i),keyerr(i),anyerr)
        if (eps(i).gt.xeps) then
            temp_min = temp(i)
            temp_max = eos_compose_temp_max
        else
            temp_min = eos_compose_temp_min
            temp_max = temp(i)
        end if

        do while (abs(eps(i)-xeps)/abs(eps(i)) .gt. epstol)
            temp_current = temp_min + (temp_max-temp_min)*inv_grp1
            call EOS_CompOSE_Eps(eoskey,keytemp,rf_precision,1, &
                rho(i),xeps,temp_current,ye(i),keyerr(i),anyerr)
            if (xeps.gt.eps(i)) then
                temp_max = temp_current
            else
                temp_min = temp_current
            end if
        end do
        temp(i) = temp_current
    end do

end subroutine EOS_CompOSE_temp_from_eps