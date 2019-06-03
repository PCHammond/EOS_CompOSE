#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EOS_CompOSE_Test_Table(CCTK_ARGUMENTS)
    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    character(len=512) :: message

    integer, allocatable :: test_ye_bin(:)
    integer, allocatable :: test_rho_bin(:)
    integer, allocatable :: test_temp_bin(:)

    real(kind=8), allocatable :: test_ye_pos(:)
    real(kind=8), allocatable :: test_rho_pos(:)
    real(kind=8), allocatable :: test_temp_pos(:)

    real :: rand
    integer :: i, j, k, offset(3)
    real(kind=8) :: results_current(8,5)

    real(kind=8) :: ye_offset, rho_offset, temp_offset
    real(kind=8) :: ye_current, rho_current, temp_current

    integer, dimension(3,8) :: offset_matrix = reshape( (/ &
        0,0,0, 1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1, 1,1,1 /), &
        shape(offset_matrix))

    integer :: x1_offset, x2_offset, x3_offset
    real(kind=8) :: x1_temp, x2_temp, x3_temp

    allocate(test_ye_bin(test_eos_table_count))
    allocate(test_rho_bin(test_eos_table_count))
    allocate(test_temp_bin(test_eos_table_count))
    allocate(test_ye_pos(test_eos_table_count))
    allocate(test_rho_pos(test_eos_table_count))
    allocate(test_temp_pos(test_eos_table_count))

    call random_seed()

    call CCTK_INFO("ye Array:")
    do i=1,EOS_CompOSE_yeCount
        write(message,"(I3,A3,E17.8E3)") i," : ", EOS_CompOSE_yeArray(i)
        call CCTK_INFO(message)
    end do

    call CCTK_INFO("log10(rho) Array:")
    do i=1,EOS_CompOSE_rhoCount
        write(message,"(I3,A3,E17.8E3)") i," : ",EOS_CompOSE_rhoArray(i)
        call CCTK_INFO(message)
    end do

    call CCTK_INFO("log10(temp) Array:")
    do i=1,EOS_CompOSE_tempCount
        write(message,"(I3,A3,E17.8E3)") i," : ", EOS_CompOSE_tempArray(i)
        call CCTK_INFO(message)
    end do

    write(message,*) "d_ye: ", EOS_CompOSE_d_ye
    call CCTK_INFO(message)

    write(message,*) "d_rho: ", EOS_CompOSE_d_rho
    call CCTK_INFO(message)

    write(message,*) "d_temp: ", EOS_CompOSE_d_temp
    call CCTK_INFO(message)

    call CCTK_INFO("Begin Interp Coeff Test")
    do i=1,64
        write(message,"(I2,A3,E17.8E3)") i," : ",real(EOS_CompOSE_interpCoeffs(i,1))
        call CCTK_INFO(message)
    end do

    call CCTK_INFO("Testing point ")
    call EOS_CompOSE_Eval_get_offset(0.1827458d0,-8.779401645324358d0,-5.0d0, &
        EOS_CompOSE_d_ye,EOS_CompOSE_d_rho,EOS_CompOSE_d_temp, &
        EOS_CompOSE_yeArray,EOS_CompOSE_rhoArray,EOS_CompOSE_tempArray, &
        EOS_CompOSE_yeCount,EOS_CompOSE_rhoCount,EOS_CompOSE_tempCount, &
        x1_offset,x2_offset,x3_offset,x1_temp,x2_temp,x3_temp)
    write(message,*) x1_offset, x2_offset, x3_offset
    call CCTK_INFO(message)
    write(message,*) x1_temp, x2_temp, x3_temp
    call CCTK_INFO(message)

    call EOS_CompOSE_Eval(0.1827458d0,-8.779401645324358d0,-5.0d0, &
        EOS_CompOSE_epsInterp, &
        EOS_CompOSE_d_ye,EOS_CompOSE_d_rho,EOS_CompOSE_d_temp, &
        EOS_CompOSE_yeArray,EOS_CompOSE_rhoArray,EOS_CompOSE_tempArray, &
        EOS_CompOSE_yeCount,EOS_CompOSE_rhoCount,EOS_CompOSE_tempCount, &
        x1_temp)
    write(message,*) x1_temp
    call CCTK_INFO(message)

    do i=1, test_eos_table_count
        call random_number(rand)
        test_ye_bin(i) = int(rand*(EOS_CompOSE_yeCount-6))+10
        test_ye_pos(i) = EOS_CompOSE_yeArray(test_ye_bin(i))
        call random_number(rand)
        test_rho_bin(i) = int(rand*(EOS_CompOSE_rhoCount-6))+10
        test_rho_pos(i) = EOS_CompOSE_rhoArray(test_rho_bin(i))
        call random_number(rand)
        test_temp_bin(i) = int(rand*(EOS_CompOSE_tempCount-6))+10
        test_temp_pos(i) = EOS_CompOSE_tempArray(test_temp_bin(i))
    end do

    ye_offset = test_eos_table_offset*EOS_CompOSE_d_ye
    rho_offset = test_eos_table_offset*EOS_CompOSE_d_rho
    temp_offset = test_eos_table_offset*EOS_CompOSE_d_temp

    do i=1, test_eos_table_count
        write(message,"(A14,I3,A2,I3,A2,I3)") "Testing point ", &
            test_ye_bin(i), ", ", test_rho_bin(i), ", ", test_temp_bin(i)
        call CCTK_INFO(message)
        write(message,"(A14,E17.8E3,A2,E17.8E3,A2,E17.8E3)") "Testing point ", &
            test_ye_pos(i), ", ", 10.0d0**test_rho_pos(i), ", ", 10.0d0**test_temp_pos(i)
        call CCTK_INFO(message)

        write(message,"(A8,E17.8E3)") "press = ", &
            EOS_CompOSE_pressTable(test_ye_bin,test_rho_bin,test_temp_bin)
        call CCTK_INFO(message)

        write(message,"(A6,E17.8E3)") "eps = ", &
            EOS_CompOSE_epsTable(test_ye_bin,test_rho_bin,test_temp_bin)
        call CCTK_INFO(message)

        write(message,"(A10,E17.8E3)") "dP_dRho = ", &
            EOS_CompOSE_dP_dRhoTable(test_ye_bin,test_rho_bin,test_temp_bin)
        call CCTK_INFO(message)

        write(message,"(A10,E17.8E3)") "dP_deps = ", &
            EOS_CompOSE_dP_depsTable(test_ye_bin,test_rho_bin,test_temp_bin)
        call CCTK_INFO(message)

        write(message,"(A7,E17.8E3)") "c_s2 = ", &
            EOS_CompOSE_c_s2Table(test_ye_bin,test_rho_bin,test_temp_bin)
        call CCTK_INFO(message)

        do j=1,8
            offset = offset_matrix(:,j)
            ye_current = test_ye_pos(i) + (2.0d0 * offset(1) - 1.0d0)*ye_offset
            rho_current = test_rho_pos(i) + (2.0d0 * offset(2) - 1.0d0)*rho_offset
            temp_current = test_temp_pos(i) + (2.0d0 * offset(3) - 1.0d0)*temp_offset

            write(message,"(A11,E17.8E3,A2,E17.8E3,A2,E17.8E3)") "Working on ", &
                ye_current, ", ", 10.0d0**rho_current, ", ", 10.0d0**temp_current
            call CCTK_INFO(message)

            call EOS_CompOSE_Eval(ye_current, rho_current, temp_current, &
                EOS_CompOSE_pressInterp, EOS_CompOSE_d_ye, &
                EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, results_current(j,1))

            call EOS_CompOSE_Eval(ye_current, rho_current, temp_current, &
                EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, results_current(j,2))

            call EOS_CompOSE_Eval(ye_current, rho_current, temp_current, &
                EOS_CompOSE_dP_dRhoInterp, EOS_CompOSE_d_ye, &
                EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, results_current(j,3))

            call EOS_CompOSE_Eval(ye_current, rho_current, temp_current, &
                EOS_CompOSE_dP_depsInterp, EOS_CompOSE_d_ye, &
                EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, results_current(j,4))

            call EOS_CompOSE_Eval(ye_current, rho_current, temp_current, &
                EOS_CompOSE_c_s2Interp, EOS_CompOSE_d_ye, &
                EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, results_current(j,5))
        end do
        do j=1,5
            do k=1,8
                write(message,"(I2,A1,I2,A1,I2,A1,E17.8E3)"), i, ",", j, ",", k, ":", results_current(k,j)
                call CCTK_INFO(message)
            end do
        end do
    end do
end subroutine EOS_CompOSE_Test_Table

subroutine EOS_CompOSE_Create_Test_Table(CCTK_ARGUMENTS)
    use EOS_CompOSE_Module
    use EOS_CompOSE_Eval_Module
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    real(kind=8) :: interp_d_ye, interp_d_rho, interp_d_temp
    real(kind=8), allocatable :: interp_ye_array(:), interp_rho_array(:), interp_temp_array(:)
    integer :: interp_count_ye, interp_count_rho, interp_count_temp

    real(kind=8), dimension(8) :: current_data
    integer :: i,j,k,l

    character(len=256) :: eosTestTableFilename
    integer :: slength
    logical :: table_exists
    character(len=512) :: message

    interp_d_ye   = 0.25d0*EOS_CompOSE_d_ye
    interp_d_rho  = 0.25d0*EOS_CompOSE_d_rho
    interp_d_temp = 0.25d0*EOS_CompOSE_d_temp

    interp_count_ye   = 4*(EOS_CompOSE_yeCount-1)   + 1
    interp_count_rho  = 4*(EOS_CompOSE_rhoCount-1)  + 1
    interp_count_temp = 4*(EOS_CompOSE_tempCount-1) + 1

    allocate(interp_temp_array(interp_count_temp))
    do i=1,interp_count_temp
        interp_temp_array(i) = EOS_CompOSE_tempArray(((i-1)/4) + 1) + mod(i-1,4)*interp_d_temp
    end do

    allocate(interp_rho_array(interp_count_rho))
    do j=1,interp_count_rho
        interp_rho_array(j) = EOS_CompOSE_rhoArray(((j-1)/4) + 1) + mod(j-1,4)*interp_d_rho
    end do

    allocate(interp_ye_array(interp_count_ye))
    do k=1,interp_count_ye
        interp_ye_array(k) = EOS_CompOSE_yeArray(((k-1)/4) + 1) + mod(k-1,4)*interp_d_ye
    end do

    call CCTK_FortranString(slength, eos_compose_test_table_name, eosTestTableFilename)
    call CCTK_INFO("Writing 4x Interpolated Table to:")
    call CCTK_INFO(eosTestTableFilename)
    call CCTK_INFO("Table will have dimensions:")
    write(message,*) interp_count_temp, interp_count_rho, interp_count_ye, 8
    call CCTK_INFO(message)
    if (CCTK_MyProc(cctkGH)==0) then
    open (unit=669,file=eosTestTableFilename,status="new",access="stream")

    do i=1,interp_count_temp
        write(message,*) "Current temp index: ", i
        call CCTK_INFO(message)
        current_data(1) = 10**interp_temp_array(i)
        do j=1,interp_count_rho
            current_data(2) = 10**interp_rho_array(j)
            do k=1,interp_count_ye
                current_data(3) = interp_ye_array(k)
                call EOS_CompOSE_Eval(interp_ye_array(k), interp_rho_array(j), interp_temp_array(i), &
                    EOS_CompOSE_pressInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, current_data(4))
                call EOS_CompOSE_Eval(interp_ye_array(k), interp_rho_array(j), interp_temp_array(i), &
                    EOS_CompOSE_epsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, current_data(5))
                call EOS_CompOSE_Eval(interp_ye_array(k), interp_rho_array(j), interp_temp_array(i), &
                    EOS_CompOSE_dP_dRhoInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, current_data(6))
                call EOS_CompOSE_Eval(interp_ye_array(k), interp_rho_array(j), interp_temp_array(i), &
                    EOS_CompOSE_dP_depsInterp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, current_data(7))
                call EOS_CompOSE_Eval(interp_ye_array(k), interp_rho_array(j), interp_temp_array(i), &
                    EOS_CompOSE_c_s2Interp, EOS_CompOSE_d_ye, &
                    EOS_CompOSE_d_rho, EOS_CompOSE_d_temp, &
                    EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, &
                    EOS_CompOSE_tempArray, EOS_CompOSE_yeCount, &
                    EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, current_data(8))
                do l=1,8
                    write (669) current_data(l)
                end do
            end do
        end do
    end do

    close(669)
    end if
    call CCTK_INFO("Interpolated table written")

end subroutine EOS_CompOSE_Create_Test_Table

subroutine EOS_CompOSE_DumpInterpCoeffs(CCTK_ARGUMENTS)
    use EOS_CompOSE_Module
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    character(len=512) :: message

    integer :: slength
    character(len=256) :: eosInterpCoeffFilename_base

    character(len=512) :: eosInterpCoeffFilename_press, eosInterpCoeffFilename_eps
    character(len=512) :: eosInterpCoeffFilename_dpdr, eosInterpCoeffFilename_dpde
    character(len=512) :: eosInterpCoeffFilename_cs2

    integer :: i,j,k,l

    call CCTK_FortranString(slength, eos_compose_interp_coeffs_file, eosInterpCoeffFilename_base)

    write(eosInterpCoeffFilename_press,"(A,A)") trim(adjustl(eosInterpCoeffFilename_base)), ".press.dat"
    write(eosInterpCoeffFilename_eps,"(A,A)") trim(adjustl(eosInterpCoeffFilename_base)), ".eps.dat"
    write(eosInterpCoeffFilename_dpdr,"(A,A)") trim(adjustl(eosInterpCoeffFilename_base)), ".dpdr.dat"
    write(eosInterpCoeffFilename_dpde,"(A,A)") trim(adjustl(eosInterpCoeffFilename_base)), ".dpde.dat"
    write(eosInterpCoeffFilename_cs2,"(A,A)") trim(adjustl(eosInterpCoeffFilename_base)), ".cs2.dat"

    call CCTK_INFO("Writing pressure interpolator to :")
    call CCTK_INFO(eosInterpCoeffFilename_press)
    call CCTK_INFO("Writing eps interpolator to :")
    call CCTK_INFO(eosInterpCoeffFilename_eps)
    call CCTK_INFO("Writing dpdr interpolator to :")
    call CCTK_INFO(eosInterpCoeffFilename_dpdr)
    call CCTK_INFO("Writing dpde interpolator to :")
    call CCTK_INFO(eosInterpCoeffFilename_dpde)
    call CCTK_INFO("Writing cs2 interpolator to :")
    call CCTK_INFO(eosInterpCoeffFilename_cs2)

    if (CCTK_MyProc(cctkGH)==0) then
    open (unit=700,file=eosInterpCoeffFilename_press,status="replace",access="stream")
    open (unit=701,file=eosInterpCoeffFilename_eps,status="replace",access="stream")
    open (unit=702,file=eosInterpCoeffFilename_dpdr,status="replace",access="stream")
    open (unit=703,file=eosInterpCoeffFilename_dpde,status="replace",access="stream")
    open (unit=704,file=eosInterpCoeffFilename_cs2,status="replace",access="stream")

    do i=1, EOS_CompOSE_tempCount-1
        write(message,*) "Current temp index: ", i
        call CCTK_INFO(message)
        do j=1, EOS_CompOSE_rhoCount-1
            do k=1, EOS_CompOSE_yeCount-1
                do l=1,64
                    write(700) EOS_CompOSE_pressInterp(l,k,j,i)
                    write(701) EOS_CompOSE_epsInterp(l,k,j,i)
                    write(702) EOS_CompOSE_dP_dRhoInterp(l,k,j,i)
                    write(703) EOS_CompOSE_dP_depsInterp(l,k,j,i)
                    write(704) EOS_CompOSE_c_s2Interp(l,k,j,i)
                end do
            end do
        end do
    end do

end if

    call CCTK_INFO("Interpolation tables written")
end subroutine EOS_CompOSE_DumpInterpCoeffs
