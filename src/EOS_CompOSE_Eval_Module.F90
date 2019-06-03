module EOS_CompOSE_Eval_Module
    implicit none
    contains
    subroutine EOS_CompOSE_Eval(x1,x2,x3,polycoeffs,d_x1,d_x2,d_x3,x1_array,x2_array,x3_array,x1_count,x2_count,x3_count,result)
        integer, intent(in) :: x1_count, x2_count, x3_count
        real(kind=8), intent(in) :: x1, x2, x3, d_x1, d_x2, d_x3
        real(kind=8), intent(in) :: polycoeffs(64,x1_count-1,x2_count-1,x3_count-1)
        real(kind=8), intent(in) :: x1_array(x1_count), x2_array(x2_count), x3_array(x3_count)
        real(kind=8), intent(out) :: result
        integer :: i, x1_offset, x2_offset, x3_offset

        real(kind=8) :: x1_temp, x2_temp, x3_temp

        x1_offset = int((x1-x1_array(1))/d_x1)+1
        x2_offset = int((x2-x2_array(1))/d_x2)+1
        x3_offset = int((x3-x3_array(1))/d_x3)+1

        if (x1_offset.lt.1) then
            x1_offset = 1
            x1_temp = 0.0d0
            !print *, "x1 under bounds"
        else if (x1_offset.gt.x1_count-1) then
            x1_offset = x1_count-1
            x1_temp = 1.0d0
            !print *, "x1 above bounds"
        else
            x1_temp = (x1 - x1_array(x1_offset))/d_x1
        end if

        !if (x1_temp.lt.0.0d0) then
        !    x1_temp = 0.0d0
        !else if (x1_temp.gt.1.0d0) then
        !    x1_temp = 1.0d0
        !end if

        if (x2_offset.lt.1) then
            x2_offset = 1
            x2_temp = 0.0d0
            !print *, "x2 under bounds"
        else if (x2_offset.gt.x2_count-1) then
            x2_offset = x2_count-1
            x2_temp = 1.0d0
            !print *, "x2 above bounds"
        else
            x2_temp = (x2 - x2_array(x2_offset))/d_x2
        end if

        !if (x2_temp.lt.0.0d0) then
        !    x2_temp = 0.0d0
        !else if (x2_temp.gt.1.0d0) then
        !    x2_temp = 1.0d0
        !end if

        if (x3_offset.lt.1) then
            x3_offset = 1
            x3_temp = 0.0d0
            !print *, "x3 under bounds"
        else if (x3_offset.gt.x3_count-1) then
            x3_offset = x3_count-1
            x3_temp = 1.0d0
            !print *, "x3 above bounds"
        else
            x3_temp = (x3 - x3_array(x3_offset))/d_x3
        end if

        !if (x3_temp.lt.0.0d0) then
        !    x3_temp = 0.0d0
        !else if (x3_temp.gt.1.0d0) then
        !    x3_temp = 1.0d0
        !end if

        result = 0.0d0
        do i=0,63
            result = result + (polycoeffs(i+1,x1_offset,x2_offset,x3_offset)* &
                                  (x1_temp**(mod(i,4)))* &
                                  (x2_temp**(mod(i/4,4)))* &
                                  (x3_temp**(mod(i/16,4))))
        end do
    end subroutine EOS_CompOSE_Eval

    subroutine EOS_CompOSE_Eval_get_offset(x1,x2,x3,d_x1,d_x2,d_x3,x1_array,x2_array,x3_array, &
            x1_count,x2_count,x3_count,x1_offset,x2_offset,x3_offset,x1_temp,x2_temp,x3_temp)
        integer, intent(in) :: x1_count, x2_count, x3_count
        real(kind=8), intent(in) :: x1, x2, x3, d_x1, d_x2, d_x3
        real(kind=8), intent(in) :: x1_array(x1_count), x2_array(x2_count), x3_array(x3_count)
        integer, intent(out) :: x1_offset, x2_offset, x3_offset
        real(kind=8), intent(out) :: x1_temp, x2_temp, x3_temp

        x1_offset = int((x1-x1_array(1))/d_x1)+1
        x2_offset = int((x2-x2_array(1))/d_x2)+1
        x3_offset = int((x3-x3_array(1))/d_x3)+1

        if (x1_offset.lt.1) then
            x1_offset = 1
            x1_temp = 0.0d0
            !print *, "x1 under bounds"
        else if (x1_offset.gt.x1_count-1) then
            x1_offset = x1_count-1
            x1_temp = 1.0d0
            !print *, "x1 above bounds"
        else
            x1_temp = (x1 - x1_array(x1_offset))/d_x1
        end if

        !if (x1_temp.lt.0.0d0) then
        !    x1_temp = 0.0d0
        !else if (x1_temp.gt.1.0d0) then
        !    x1_temp = 1.0d0
        !end if

        if (x2_offset.lt.1) then
            x2_offset = 1
            x2_temp = 0.0d0
            !print *, "x2 under bounds"
        else if (x2_offset.gt.x2_count-1) then
            x2_offset = x2_count-1
            x2_temp = 1.0d0
            !print *, "x2 above bounds"
        else
            x2_temp = (x2 - x2_array(x2_offset))/d_x2
        end if

        !if (x2_temp.lt.0.0d0) then
        !    x2_temp = 0.0d0
        !else if (x2_temp.gt.1.0d0) then
        !    x2_temp = 1.0d0
        !end if

        if (x3_offset.lt.1) then
            x3_offset = 1
            x3_temp = 0.0d0
            !print *, "x3 under bounds"
        else if (x3_offset.gt.x3_count-1) then
            x3_offset = x3_count-1
            x3_temp = 1.0d0
            !print *, "x3 above bounds"
        else
            x3_temp = (x3 - x3_array(x3_offset))/d_x3
        end if

        !if (x3_temp.lt.0.0d0) then
        !    x3_temp = 0.0d0
        !else if (x3_temp.gt.1.0d0) then
        !    x3_temp = 1.0d0
        !end if
    end subroutine EOS_CompOSE_Eval_get_offset
end module EOS_CompOSE_Eval_Module
