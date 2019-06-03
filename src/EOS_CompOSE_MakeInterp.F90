#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EOS_CompOSE_LoadTable(CCTK_ARGUMENTS)
    use EOS_CompOSE_Module
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    character(len=512) :: warnline
    character(len=512) :: message
    character(len=256) :: eosTableFilename
    integer :: slength
    logical :: table_exists

    integer :: rowCount, colCount, row
    integer :: tempIdx, rhoIdx, yeIdx
    real(kind=8), allocatable :: tableEOSRaw(:,:)

    integer, dimension(4) :: table_shape
    integer, dimension(3) :: variable_shape

    EOS_CompOSE_yeCount = eos_compose_tableshape(1)
    EOS_CompOSE_rhoCount = eos_compose_tableshape(2)
    EOS_CompOSE_tempCount = eos_compose_tableshape(3)
    EOS_CompOSE_varCount = eos_compose_tableshape(4)

    call CCTK_FortranString(slength, eos_compose_table_name, eosTableFilename)

    inquire(file=trim(adjustl(eosTableFilename)), exist=table_exists)

    if(.not.table_exists) then
        write(warnline,"(A10,A,A15)") "EOS file ", &
            trim(adjustl(eosTableFilename)), " does not exist!"
        call CCTK_ERROR(warnline)
        stop
    endif

    rowCount = EOS_CompOSE_tempCount*EOS_CompOSE_rhoCount*EOS_CompOSE_yeCount
    colCount = EOS_CompOSE_varCount

    write(message,"(A14,A)") "Opening table ", trim(adjustl(eosTableFilename))
    call CCTK_INFO(message)

    open(667, file=eosTableFilename, status="old", action="read")
    allocate(tableEOSRaw(colCount, rowCount))
    do row=1, rowCount
        read(667,*) tableEOSRaw(:,row)
    end do

    write(message,"(A19,I4,A2,I4,A2,I4,A2,I4)") "Reshaping table to ", &
        EOS_CompOSE_yeCount, " x", EOS_CompOSE_rhoCount, " x", &
        EOS_CompOSE_tempCount, " x", EOS_CompOSE_varCount
    call CCTK_INFO(message)

    allocate(EOS_CompOSE_table(EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, &
        EOS_CompOSE_tempCount, EOS_CompOSE_varCount))
    do row=1, rowCount
        yeIdx = mod(row-1,EOS_CompOSE_yeCount)+1
        rhoIdx = mod((row-1)/EOS_CompOSE_yeCount,EOS_CompOSE_rhoCount)+1
        tempIdx = ((row-1)/(EOS_CompOSE_yeCount*EOS_CompOSE_rhoCount))+1

        !Ye
        EOS_CompOSE_table(yeIdx,rhoIdx,tempIdx,1) = tableEOSRaw(3,row)

        !Density
        EOS_CompOSE_table(yeIdx,rhoIdx,tempIdx,2) = &
            log10(tableEOSRaw(2,row)*unit_rho_tabToCode)

        !Temperature
        EOS_CompOSE_table(yeIdx,rhoIdx,tempIdx,3) = log10(tableEOSRaw(1,row))

        !Pressure
        EOS_CompOSE_table(yeIdx,rhoIdx,tempIdx,4) = &
            log10(tableEOSRaw(4,row)*unit_press_tabToCode)

        !epsilon
        EOS_CompOSE_table(yeIdx,rhoIdx,tempIdx,5) = tableEOSRaw(5,row)

        !dP / drho
        EOS_CompOSE_table(yeIdx,rhoIdx,tempIdx,6) = &
            log10(tableEOSRaw(6,row)*unit_dP_dnb_tabToCode)

        !dP / deps
        EOS_CompOSE_table(yeIdx,rhoIdx,tempIdx,7) = &
            log10(tableEOSRaw(7,row)*unit_dP_deps_tabToCode)

        !c_s ** 2
        EOS_CompOSE_table(yeIdx,rhoIdx,tempIdx,8) = tableEOSRaw(8,row)
    end do

    allocate(EOS_CompOSE_yeArray(EOS_CompOSE_yeCount))
    allocate(EOS_CompOSE_rhoArray(EOS_CompOSE_rhoCount))
    allocate(EOS_CompOSE_tempArray(EOS_CompOSE_tempCount))
    EOS_CompOSE_yeArray = EOS_CompOSE_table(:,1,1,1)
    EOS_CompOSE_rhoArray = EOS_CompOSE_table(1,:,1,2)
    EOS_CompOSE_tempArray = EOS_CompOSE_table(1,1,:,3)
    EOS_CompOSE_d_ye = EOS_CompOSE_yeArray(2) - EOS_CompOSE_yeArray(1)
    EOS_CompOSE_d_rho = EOS_CompOSE_rhoArray(2) - EOS_CompOSE_rhoArray(1)
    EOS_CompOSE_d_temp = EOS_CompOSE_tempArray(2) - EOS_CompOSE_tempArray(1)

    call CCTK_INFO("Building pressure interpolator")
    !Pressure interpolator
    allocate(EOS_CompOSE_pressTable(EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, &
        EOS_CompOSE_tempCount))
    EOS_CompOSE_pressTable = EOS_CompOSE_table(:,:,:,4)
    allocate(EOS_CompOSE_pressInterp(64, EOS_CompOSE_yeCount-1, &
        EOS_CompOSE_rhoCount-1, EOS_CompOSE_tempCount-1))
    call EOS_CompOSE_BuildInterpolator(EOS_CompOSE_PressTable, &
        EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, EOS_CompOSE_tempArray, &
        EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, &
        EOS_CompOSE_d_ye,    EOS_CompOSE_d_rho,    EOS_CompOSE_d_temp, &
        EOS_CompOSE_pressInterp)

    write(message,"(A28,E17.8E3)") "Maximum Pressure in table = ", &
        maxval(EOS_CompOSE_pressInterp(1,:,:,:))
    call CCTK_INFO(message)
    write(message,"(A28,E17.8E3)") "Minimum Pressure in table = ", &
        minval(EOS_CompOSE_pressInterp(1,:,:,:))
    call CCTK_INFO(message)

    call CCTK_INFO("Building eps interpolator")
    !epsilon interpolator
    allocate(EOS_CompOSE_epsTable(EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, &
        EOS_CompOSE_tempCount))
    EOS_CompOSE_epsTable = EOS_CompOSE_table(:,:,:,5)
    allocate(EOS_CompOSE_epsInterp(64, EOS_CompOSE_yeCount-1, &
        EOS_CompOSE_rhoCount-1, EOS_CompOSE_tempCount-1))
    call EOS_CompOSE_BuildInterpolator(EOS_CompOSE_epsTable, &
        EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, EOS_CompOSE_tempArray, &
        EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, &
        EOS_CompOSE_d_ye,    EOS_CompOSE_d_rho,    EOS_CompOSE_d_temp, &
        EOS_CompOSE_epsInterp)

    write(message,"(A27,E17.8E3)") "Maximum epsilon in table = ", &
        maxval(EOS_CompOSE_epsInterp(1,:,:,:))
    call CCTK_INFO(message)
    write(message,"(A27,E17.8E3)") "Minimum epsilon in table = ", &
        minval(EOS_CompOSE_epsInterp(1,:,:,:))
    call CCTK_INFO(message)
    write(message,*) EOS_CompOSE_table(18,81,1,5)
    call CCTK_INFO(message)
    write(message,*) EOS_CompOSE_epsInterp(1,18,81,1)
    call CCTK_INFO(message)

    call CCTK_INFO("Building dP/drho interpolator")
    !dP / drho interpolator
    allocate(EOS_CompOSE_dP_dRhoTable(EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, &
        EOS_CompOSE_tempCount))
    EOS_CompOSE_dP_dRhoTable = EOS_CompOSE_table(:,:,:,6)
    allocate(EOS_CompOSE_dP_dRhoInterp(64, EOS_CompOSE_yeCount-1, &
        EOS_CompOSE_rhoCount-1, EOS_CompOSE_tempCount-1))
    call EOS_CompOSE_BuildInterpolator(EOS_CompOSE_dP_dRhoTable, &
        EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, EOS_CompOSE_tempArray, &
        EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, &
        EOS_CompOSE_d_ye,    EOS_CompOSE_d_rho,    EOS_CompOSE_d_temp, &
        EOS_CompOSE_dP_dRhoInterp)

    write(message,"(A27,E17.8E3)") "Maximum dP/drho in table = ", &
        maxval(EOS_CompOSE_dP_dRhoInterp(1,:,:,:))
    call CCTK_INFO(message)
    write(message,"(A27,E17.8E3)") "Minimum dP/drho in table = ", &
        minval(EOS_CompOSE_dP_dRhoInterp(1,:,:,:))
    call CCTK_INFO(message)

    call CCTK_INFO("Building dP/deps interpolator")
    !dP / deps interpolator
    allocate(EOS_CompOSE_dP_depsTable(EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, &
        EOS_CompOSE_tempCount))
    EOS_CompOSE_dP_depsTable = EOS_CompOSE_table(:,:,:,7)
    allocate(EOS_CompOSE_dP_depsInterp(64, EOS_CompOSE_yeCount-1, &
        EOS_CompOSE_rhoCount-1, EOS_CompOSE_tempCount-1))
    call EOS_CompOSE_BuildInterpolator(EOS_CompOSE_dP_depsTable, &
        EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, EOS_CompOSE_tempArray, &
        EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, &
        EOS_CompOSE_d_ye,    EOS_CompOSE_d_rho,    EOS_CompOSE_d_temp, &
        EOS_CompOSE_dP_depsInterp)

    write(message,"(A27,E17.8E3)") "Maximum dP/deps in table = ", &
        maxval(EOS_CompOSE_dP_depsInterp(1,:,:,:))
    call CCTK_INFO(message)
    write(message,"(A27,E17.8E3)") "Minimum dP/deps in table = ", &
        minval(EOS_CompOSE_dP_depsInterp(1,:,:,:))
    call CCTK_INFO(message)

    call CCTK_INFO("Building c_s^2 interpolator")
    ! c_s ** 2 interpolator
    allocate(EOS_CompOSE_c_s2Table(EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, &
        EOS_CompOSE_tempCount))
    EOS_CompOSE_c_s2Table = EOS_CompOSE_table(:,:,:,8)
    allocate(EOS_CompOSE_c_s2Interp(64, EOS_CompOSE_yeCount-1, &
        EOS_CompOSE_rhoCount-1, EOS_CompOSE_tempCount-1))
    call EOS_CompOSE_BuildInterpolator(EOS_CompOSE_c_s2Table, &
        EOS_CompOSE_yeArray, EOS_CompOSE_rhoArray, EOS_CompOSE_tempArray, &
        EOS_CompOSE_yeCount, EOS_CompOSE_rhoCount, EOS_CompOSE_tempCount, &
        EOS_CompOSE_d_ye,    EOS_CompOSE_d_rho,    EOS_CompOSE_d_temp, &
        EOS_CompOSE_c_s2Interp)

    write(message,"(A23,E17.8E3)") "Maximum cs2 in table = ", &
        maxval(EOS_CompOSE_c_s2Interp(1,:,:,:))
    call CCTK_INFO(message)
    write(message,"(A23,E17.8E3)") "Minimum cs2 in table = ", &
        minval(EOS_CompOSE_c_s2Interp(1,:,:,:))
    call CCTK_INFO(message)

    deallocate(tableEOSRaw)
end subroutine EOS_CompOSE_LoadTable

subroutine EOS_CompOSE_BuildInterpolator(variable, &
                                         yeArray, rhoArray, tempArray, &
                                         yeCount, rhoCount, tempCount, &
                                         d_ye,    d_rho,    d_temp, &
                                         poly_coeffs_final)
    use EOS_CompOSE_Module
    use EOS_CompOSE_Differences_Module
    use EOS_CompOSE_InterpCoeffs_Module
    implicit none

    DECLARE_CCTK_FUNCTIONS

    !Input variables
    real(kind=8), intent(in) :: d_ye, d_rho, d_temp
    integer, intent(in) :: yeCount, rhoCount, tempCount
    real(kind=8), intent(in) :: variable(yeCount,rhoCount,tempCount)
    real(kind=8), intent(in) :: yeArray(yeCount), rhoArray(rhoCount), tempArray(tempCount)

    !Temporary variables
    real(kind=8), allocatable :: diff_array(:,:,:,:)
    real(kind=8), allocatable :: interp_matrix(:,:,:,:)
    real(kind=8), allocatable :: poly_coeffs_temp(:,:,:,:)
    !integer, dimension(3,8) :: offset_matrix = reshape( (/ &
    !    0,0,0, 0,0,1, 0,1,0, 0,1,1, 1,0,0, 1,0,1, 1,1,0, 1,1,1 /), &
    !    shape(offset_matrix))
    integer, dimension(3,8) :: offset_matrix = reshape( (/ &
        0,0,0, 1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1, 1,1,1 /), &
        shape(offset_matrix))
    integer, dimension(3) :: offset
    integer :: tempIdx, rhoIdx, yeIdx
    integer :: i, j

    integer, dimension(3) :: variable_shape
    character(len=512) :: message

    !Output variables
    real(kind=8), intent(out) :: poly_coeffs_final(64,yeCount-1,rhoCount-1,tempCount-1)

    !Calculate finite differences
    allocate(diff_array(yeCount,rhoCount,tempCount,8))
    call EOS_CompOSE_DiffArray(variable,yeCount,rhoCount,tempCount,diff_array)

    !Get cube coefficients
    allocate(interp_matrix(yeCount-1,rhoCount-1,tempCount-1,64))
    !loop over corners
    do i=1,8
        offset = offset_matrix(:,i)
        !loop over diffs
        do j=0,7
            interp_matrix(:,:,:,i+8*j) = &
                diff_array(offset(1)+1:yeCount   + offset(1)-1, &
                           offset(2)+1:rhoCount  + offset(2)-1, &
                           offset(3)+1:tempCount + offset(3)-1,j+1)
        end do
    end do

    !Calculate polynomial coefficients
    allocate(poly_coeffs_temp(yeCount-1,rhoCount-1,tempCount-1,64))
    do j=1,64
        poly_coeffs_temp(:,:,:,j) = 0.0
        do i=1,64
            poly_coeffs_temp(:,:,:,j) = poly_coeffs_temp(:,:,:,j) + &
                EOS_CompOSE_interpCoeffs(i,j)*interp_matrix(:,:,:,i)
        end do
    end do

    !Change of variable
    do tempIdx = 1, tempCount-1
        do rhoIdx = 1, rhoCount-1
            do yeIdx = 1, yeCount-1
                do i=1,64
                    poly_coeffs_final(i, yeIdx, rhoIdx, tempIdx) = &
                        poly_coeffs_temp(yeIdx, rhoIdx, tempIdx, i)
                end do
            end do
        end do
    end do
end subroutine EOS_CompOSE_BuildInterpolator
