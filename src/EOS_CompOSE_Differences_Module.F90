#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

module EOS_CompOSE_Differences_Module
    implicit none
contains
    subroutine dF_dx1(F_array,dF_dx1_array,len_x1,len_x2,len_x3)
        real(kind=8), intent(in) :: F_array(len_x1,len_x2,len_x3)
        real(kind=8), intent(out) :: dF_dx1_array(len_x1,len_x2,len_x3)
        integer, intent(in) :: len_x1,len_x2,len_x3
        dF_dx1_array(1,:,:) = -3.0*F_array(1,:,:) + 4.0*F_array(2,:,:) - F_array(3,:,:)
        dF_dx1_array(len_x1,:,:) = F_array(len_x1-2,:,:) - 4.0*F_array(len_x1-1,:,:) + F_array(len_x1,:,:)
        dF_dx1_array(2:len_x1-1,:,:) = F_array(3:len_x1,:,:) - F_array(1:len_x1-2,:,:)
    end subroutine dF_dx1

    subroutine dF_dx2(F_array,dF_dx2_array,len_x1,len_x2,len_x3)
        real(kind=8), intent(in) :: F_array(len_x1,len_x2,len_x3)
        real(kind=8), intent(out) :: dF_dx2_array(len_x1,len_x2,len_x3)
        integer, intent(in) :: len_x1,len_x2,len_x3
        dF_dx2_array(:,1,:) = -3.0*F_array(:,1,:) + 4.0*F_array(:,2,:) - F_array(:,3,:)
        dF_dx2_array(:,len_x2,:) = F_array(:,len_x2-2,:) - 4.0*F_array(:,len_x2-1,:) + F_array(:,len_x2,:)
        dF_dx2_array(:,2:len_x2-1,:) = F_array(:,3:len_x2,:) - F_array(:,1:len_x2-2,:)
    end subroutine dF_dx2

    subroutine dF_dx3(F_array,dF_dx3_array,len_x1,len_x2,len_x3)
        real(kind=8), intent(in) :: F_array(len_x1,len_x2,len_x3)
        real(kind=8), intent(out) :: dF_dx3_array(len_x1,len_x2,len_x3)
        integer, intent(in) :: len_x1,len_x2,len_x3
        dF_dx3_array(:,:,1) = -3.0*F_array(:,:,1) + 4.0*F_array(:,:,2) - F_array(:,:,3)
        dF_dx3_array(:,:,len_x3) = F_array(:,:,len_x3-2) - 4.0*F_array(:,:,len_x3-1) + F_array(:,:,len_x3)
        dF_dx3_array(:,:,2:len_x3-1) = F_array(:,:,3:len_x3) - F_array(:,:,1:len_x3-2)
    end subroutine dF_dx3

    subroutine EOS_CompOSE_DiffArray(variable, yeCount, rhoCount, tempCount, &
                                     diff_array)
        !Input variables
        integer, intent(in) :: yeCount, rhoCount, tempCount
        real(kind=8), intent(in) :: variable(yeCount,rhoCount,tempCount)

        !Output variables
        real(kind=8), intent(out) :: diff_array(yeCount,rhoCount,tempCount,8)

        DECLARE_CCTK_FUNCTIONS

        diff_array(:,:,:,1) = variable(:,:,:)

        call dF_dx1(diff_array(:,:,:,1),diff_array(:,:,:,2),yeCount,rhoCount,tempCount) !d/dx
        call dF_dx2(diff_array(:,:,:,1),diff_array(:,:,:,3),yeCount,rhoCount,tempCount) !d/dy
        call dF_dx3(diff_array(:,:,:,1),diff_array(:,:,:,4),yeCount,rhoCount,tempCount) !d/dz

        call dF_dx2(diff_array(:,:,:,2),diff_array(:,:,:,5),yeCount,rhoCount,tempCount)  !d/dy(d/dx) = dd/dxdy
        call dF_dx3(diff_array(:,:,:,2),diff_array(:,:,:,6),yeCount,rhoCount,tempCount) !d/dz(d/dx) = dd/dxdz
        call dF_dx3(diff_array(:,:,:,3),diff_array(:,:,:,7),yeCount,rhoCount,tempCount) !d/dz(d/dy) = dd/dydz

        call dF_dx3(diff_array(:,:,:,5),diff_array(:,:,:,8),yeCount,rhoCount,tempCount) !d/dz(dd/dxdy) = ddd/dxdydz

    end subroutine EOS_CompOSE_DiffArray
end module EOS_CompOSE_Differences_Module
