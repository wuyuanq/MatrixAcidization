
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_constructMat

    use DBF_model
    use DBF_globalData
    implicit none

Contains

    ! '|' equation
    ! '|' element
    ! laplace operation
    subroutine dctz_xxlap(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: field_ind, equ_ind
        integer :: indexl, indexr, indexd, indexu, indextemp
        integer :: i, j

        ! the field index
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols + 1
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        if(pcol /= pncols) then
            indextemp = localncols
        else
            indextemp = localncols + 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(field(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call coordiToGlobalInd(myid-1, ind_kind, localncols, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, 1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, 1, j)
                    elseif((pcol/=pncols).and.(i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call coordiToGlobalInd(myid+1, ind_kind, 1, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, localncols, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, localncols, j)
                    elseif((j==0).and.(i/=0).and.(i/=indextemp+1).and.(.not.((pcol==1).and.(i==1).and. &!
                        (isDiriX0_p(ylower+j)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                        (isDiriX1_p(ylower+j)==0)))) then
                        call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, 1)
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=indextemp+1).and.(.not.((pcol==1).and. &!
                        (i==1).and.(isDiriX0_p(ylower+j-2)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                        (isDiriX1_p(ylower+j-2)==0)))) then
                        call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, localnrows)
                    elseif((i>=1).and.(i<=indextemp).and.(j>=1).and.(j<=localnrows)) then
                        call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                        if(((j/=1).and.(i>=2).and.(i<=indextemp-1)).or.((j/=1).and.(i==1).and.(pcol/=1)).or. &!
                            ((j/=1).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-2)/=0)).or.((j/=1).and.(i==indextemp) &!
                            .and.(pcol/=pncols)).or.((j/=1).and.(i==indextemp).and.(pcol==pncols).and. &!
                            (isDiriX1_p(ylower+j-2)/=0))) then
                            call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        end if
                        if((i>=3).or.((i==2).and.(pcol/=1)).or.((i==2).and.(pcol==1).and.(isDiriX0_p(ylower+j-1)/=0))) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        end if
                        call setMatValue(field_ind, field_ind, resi(i,j), fi_kind, i, j)
                        if((i<=indextemp-2).or.((i==indextemp-1).and.(pcol/=pncols)).or.((i==indextemp-1).and. &!
                            (pcol==pncols).and.(isDiriX1_p(ylower+j-1)/=0))) then
                            call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                        end if
                        if(((j/=localnrows).and.(i>=2).and.(i<=indextemp-1)).or.((j/=localnrows).and.(i==1).and. &!
                            (pcol/=1)).or.((j/=localnrows).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j)/=0)).or. &!
                            ((j/=localnrows).and.(i==indextemp).and.(pcol/=pncols)).or.((j/=localnrows).and. &!
                            (i==indextemp).and.(pcol==pncols).and.(isDiriX1_p(ylower+j)/=0))) then
                            call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine dctz_xxlap

    ! '|' equation
    ! '.' element
    ! gradient operation
    subroutine dctz_xpgra(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind

        integer :: indexl, indexr, indexd, indexu
        integer :: field_ind, equ_ind
        integer :: i, j

        ! the field index
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols
        indexd = 1
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                if(field(i,j) == 1) then
                    if(i == 0) then
                        call coordiToGlobalInd(myid-1, f_ind_kind, localncols, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, 1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, 1, j)
                    elseif((i==1).and.(pcol==1)) then
                        call coordiToGlobalInd(myid, f_ind_kind, i, j, field_ind)
                        if(isDiriX0_p(ylower+j-1) /= 0) then
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        end if
                        call coordiToGlobalInd(myid, e_ind_kind, i+1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                    elseif((i==localncols).and.(pcol/=pncols)) then
                        call coordiToGlobalInd(myid, f_ind_kind, localncols, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, localncols, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, localncols, j)
                    elseif((i==localncols).and.(pcol==pncols)) then
                        call coordiToGlobalInd(myid, f_ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        if(isDiriX1_p(ylower+j-1) /= 0) then
                            call coordiToGlobalInd(myid, e_ind_kind, i+1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                        end if
                    else
                        call coordiToGlobalInd(myid, f_ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        call coordiToGlobalInd(myid, e_ind_kind, i+1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                    end if
                end if

            end do
        end do

    end subroutine dctz_xpgra

    ! '-' equation
    ! '-' element
    ! laplace operation
    subroutine dctz_yylap(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: indexl, indexr, indexd, indexu, indextemp
        integer :: field_ind, equ_ind
        integer :: i, j

        ! the field index
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows + 1

        if(prow /= pnrows) then
            indextemp = localnrows
        else
            indextemp = localnrows + 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(field(i,j) == 1) then

                    if((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, 1)
                    elseif((prow/=pnrows).and.(j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, localnrows)
                    elseif((i==0).and.(j/=0).and.(j/=indextemp+1).and.(.not.((prow==1).and.(j==1).and. &!
                        (isDiriY0_p(xlower+i)==0))).and.(.not.((prow==pnrows).and.(j==indextemp).and. &!
                        (isDiriY1_p(xlower+i)==0)))) then
                        call coordiToGlobalInd(myid-1, ind_kind, localncols, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, 1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, 1, j)
                    elseif((i==localncols+1).and.(j/=0).and.(j/=indextemp+1).and.(.not.((prow==1).and.(j==1).and. &!
                        (isDiriY0_p(xlower+i-2)==0))).and.(.not.((prow==pnrows).and.(j==indextemp).and. &!
                        (isDiriY1_p(xlower+i-2)==0)))) then
                        call coordiToGlobalInd(myid+1, ind_kind, 1, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, localncols, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, localncols, j)
                    elseif((j>=1).and.(j<=indextemp).and.(i>=1).and.(i<=localncols)) then
                        call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                        if((j>=3).or.((j==2).and.(prow/=1)).or.((j==2).and.(prow==1).and.(isDiriY0_p(xlower+i-1)/=0))) then
                            call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        end if
                        if(((i/=1).and.(j>=2).and.(j<=indextemp-1)).or.((i/=1).and.(j==1).and.(prow/=1)).or.((i/=1).and. &!
                            (j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-2)/=0)).or.((i/=1).and.(j==indextemp).and. &!
                            (prow/=pnrows)).or.((i/=1).and.(j==indextemp).and.(prow==pnrows).and. &!
                            (isDiriY1_p(xlower+i-2)/=0))) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        end if
                        call setMatValue(field_ind, field_ind, resi(i,j), fi_kind, i, j)
                        if(((i/=localncols).and.(j>=2).and.(j<=indextemp-1)).or.((i/=localncols).and.(j==1).and.(prow/=1)).or. &!
                            ((i/=localncols).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i)/=0)).or.((i/=localncols).and. &!
                            (j==indextemp).and.(prow/=pnrows)).or.((i/=localncols).and.(j==indextemp).and.(prow==pnrows) &!
                            .and.(isDiriY1_p(xlower+i)/=0))) then
                            call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                        end if
                        if((j<=indextemp-2).or.((j==indextemp-1).and.(prow/=pnrows)).or.((j==indextemp-1).and. &!
                            (prow==pnrows).and.(isDiriY1_p(xlower+i-1)/=0))) then
                            call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine dctz_yylap

    ! '-' equation
    ! '.' element
    ! gradient operation
    subroutine dctz_ypgra(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind

        integer :: indexl, indexr, indexd, indexu
        integer :: field_ind, equ_ind
        integer :: i, j

        ! the field index
        indexl = 1
        indexr = localncols
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                if(field(i,j) == 1) then
                    if(j == 0) then
                        call coordiToGlobalInd(myid-pncols, f_ind_kind, i, localnrows, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, 1)
                    elseif((j==1).and.(prow==1)) then
                        call coordiToGlobalInd(myid, f_ind_kind, i, j, field_ind)
                        if(isDiriY0_p(xlower+i-1) /= 0) then
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        end if
                        call coordiToGlobalInd(myid, e_ind_kind, i, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                    elseif((j==localnrows).and.(prow/=pnrows)) then
                        call coordiToGlobalInd(myid, f_ind_kind, i, localnrows, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, localnrows)
                    elseif((j==localnrows).and.(prow==pnrows)) then
                        call coordiToGlobalInd(myid, f_ind_kind, i, localnrows, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, localnrows)
                        if(isDiriY1_p(xlower+i-1) /= 0) then
                            call coordiToGlobalInd(myid, e_ind_kind, i, localnrows+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, localnrows+1)
                        end if
                    else
                        call coordiToGlobalInd(myid, f_ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        call coordiToGlobalInd(myid, e_ind_kind, i, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                    end if
                end if

            end do
        end do

    end subroutine dctz_ypgra

    ! '.' equation
    ! '|' element
    ! divergence operation
    subroutine dctz_pxdiv(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind

        integer :: indexl, indexr, indexd, indexu
        integer :: field_ind, equ_ind
        integer :: i, j

        ! the field index
        indexl = 1
        indexr = localncols+1
        indexd = 1
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                if(field(i,j) == 1) then
                    if(i == 1) then
                        call coordiToGlobalInd(myid, f_ind_kind, 1, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, 1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, 1, j)
                    elseif((i==localncols+1).and.(pcol/=pncols)) then
                        call coordiToGlobalInd(myid+1, f_ind_kind, 1, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, localncols, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, localncols, j)
                    elseif((i==localncols+1).and.(pcol==pncols)) then
                        call coordiToGlobalInd(myid, f_ind_kind, localncols+1, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, localncols, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, localncols, j)
                    else
                        call coordiToGlobalInd(myid, f_ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i-1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        call coordiToGlobalInd(myid, e_ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                    end if
                end if

            end do
        end do

    end subroutine dctz_pxdiv

    ! '.' equation
    ! '-' element
    ! divergence operation
    subroutine dctz_pydiv(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind

        integer :: indexl, indexr, indexd, indexu
        integer :: field_ind, equ_ind
        integer :: i, j

        ! the field index
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1

        do j = indexd, indexu
            do i = indexl, indexr

                if(field(i,j) == 1) then
                    if(j == 1) then
                        call coordiToGlobalInd(myid, f_ind_kind, i, 1, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, 1)
                    elseif((j==localnrows+1).and.(prow/=pnrows)) then
                        call coordiToGlobalInd(myid+pncols, f_ind_kind, i, 1, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, localnrows)
                    elseif((j==localnrows+1).and.(prow==pnrows)) then
                        call coordiToGlobalInd(myid, f_ind_kind, i, localnrows+1, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, localnrows)
                    else
                        call coordiToGlobalInd(myid, f_ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, e_ind_kind, i, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        call coordiToGlobalInd(myid, e_ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                    end if
                end if

            end do
        end do

    end subroutine dctz_pydiv

    ! '.' equation
    ! '.' element
    ! derivative operation
    subroutine dctz_ppder(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: field_ind
        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                if(field(i,j) == 1) then
                    call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                    call setMatValue(field_ind, field_ind, resi(i,j), fi_kind, i, j)
                end if

            end do
        end do

    end subroutine dctz_ppder

    ! '.' equation
    ! '.' element
    ! laplace operation - 9 elements
    subroutine dctz_pplap9(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: indexl, indexr, indexd, indexu
        integer :: field_ind, equ_ind
        integer :: i, j

        ! the field index
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(field(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call coordiToGlobalInd(myid-1, ind_kind, localncols, j, field_ind)
                        if((j>=2).and.(j<=localnrows-1)) then
                            call coordiToGlobalInd(myid, ind_kind, i+1, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j-1), fi_kind, i+1, j-1)
                            call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                            call coordiToGlobalInd(myid, ind_kind, i+1, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j+1), fi_kind, i+1, j+1)
                        elseif(j == 1) then
                            call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                            call coordiToGlobalInd(myid, ind_kind, i+1, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j+1), fi_kind, i+1, j+1)
                        elseif(j == localnrows) then
                            call coordiToGlobalInd(myid, ind_kind, i+1, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j-1), fi_kind, i+1, j-1)
                            call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                        end if
                    elseif((i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call coordiToGlobalInd(myid+1, ind_kind, 1, j, field_ind)
                        if((j>=2).and.(j<=localnrows-1)) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j-1), fi_kind, i-1, j-1)
                            call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                            call coordiToGlobalInd(myid, ind_kind, i-1, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j+1), fi_kind, i-1, j+1)
                        elseif(j == 1) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                            call coordiToGlobalInd(myid, ind_kind, i-1, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j+1), fi_kind, i-1, j+1)
                        elseif(j == localnrows) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j-1), fi_kind, i-1, j-1)
                            call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        end if
                    elseif((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, field_ind)
                        if((i>=2).and.(i<=localncols-1)) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j+1), fi_kind, i-1, j+1)
                            call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                            call coordiToGlobalInd(myid, ind_kind, i+1, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j+1), fi_kind, i+1, j+1)
                        elseif(i == 1) then
                            call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                            call coordiToGlobalInd(myid, ind_kind, i+1, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j+1), fi_kind, i+1, j+1)
                        elseif(i == localncols) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j+1), fi_kind, i-1, j+1)
                            call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                        end if
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, field_ind)
                        if((i>=2).and.(i<=localncols-1)) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j-1), fi_kind, i-1, j-1)
                            call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                            call coordiToGlobalInd(myid, ind_kind, i+1, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j-1), fi_kind, i+1, j-1)
                        elseif(i == 1) then
                            call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                            call coordiToGlobalInd(myid, ind_kind, i+1, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j-1), fi_kind, i+1, j-1)
                        elseif(i == localncols) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j-1), fi_kind, i-1, j-1)
                            call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        end if
                    elseif((i==0).and.(j==0)) then
                        call coordiToGlobalInd(myid-1-pncols, ind_kind, localncols, localnrows, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, 1, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(1,1), fi_kind, 1, 1)
                    elseif((i==0).and.(j==localnrows+1)) then
                        call coordiToGlobalInd(myid-1+pncols, ind_kind, localncols, 1, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, 1, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(1,localnrows), fi_kind, 1, localnrows)
                    elseif((i==localncols+1).and.(j==0)) then
                        call coordiToGlobalInd(myid+1-pncols, ind_kind, 1, localnrows, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, localncols, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols,1), fi_kind, localncols, 1)
                    elseif((i==localncols+1).and.(j==localnrows+1)) then
                        call coordiToGlobalInd(myid+1+pncols, ind_kind, 1, 1, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, localncols, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols,localnrows), fi_kind, localncols, localnrows)
                    elseif((i==1).and.(j==1)) then
                        call coordiToGlobalInd(myid, ind_kind, 1, 1, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, 1, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(1,1), fi_kind, 1, 1)
                        call coordiToGlobalInd(myid, ind_kind, 2, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(2,1), fi_kind, 2, 1)
                        call coordiToGlobalInd(myid, ind_kind, 1, 2, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(1,2), fi_kind, 1, 2)
                        call coordiToGlobalInd(myid, ind_kind, 2, 2, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(2,2), fi_kind, 2, 2)
                    elseif((i==1).and.(j==localnrows)) then
                        call coordiToGlobalInd(myid, ind_kind, 1, localnrows, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, 1, localnrows-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(1,localnrows-1), fi_kind, 1, localnrows-1)
                        call coordiToGlobalInd(myid, ind_kind, 2, localnrows-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(2,localnrows-1), fi_kind, 2, localnrows-1)
                        call coordiToGlobalInd(myid, ind_kind, 1, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(1,localnrows), fi_kind, 1, localnrows)
                        call coordiToGlobalInd(myid, ind_kind, 2, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(2,localnrows), fi_kind, 2, localnrows)
                    elseif((i==localncols).and.(j==1)) then
                        call coordiToGlobalInd(myid, ind_kind, localncols, 1, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, localncols-1, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols-1,1), fi_kind, localncols-1, 1)
                        call coordiToGlobalInd(myid, ind_kind, localncols, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols,1), fi_kind, localncols, 1)
                        call coordiToGlobalInd(myid, ind_kind, localncols-1, 2, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols-1,2), fi_kind, localncols-1, 2)
                        call coordiToGlobalInd(myid, ind_kind, localncols, 2, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols, 2), fi_kind, localncols, 2)
                    elseif((i==localncols).and.(j==localnrows)) then
                        call coordiToGlobalInd(myid, ind_kind, localncols, localnrows, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, localncols-1, localnrows-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols-1,localnrows-1), fi_kind, localncols-1, localnrows-1)
                        call coordiToGlobalInd(myid, ind_kind, localncols, localnrows-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols,localnrows-1), fi_kind, localncols, localnrows-1)
                        call coordiToGlobalInd(myid, ind_kind, localncols-1, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols-1,localnrows), fi_kind, localncols-1, localnrows)
                        call coordiToGlobalInd(myid, ind_kind, localncols, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(localncols, localnrows), fi_kind, localncols, localnrows)
                    elseif((i==1).and.(j>=2).and.(j<=localnrows-1)) then
                        call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j-1), fi_kind, i+1, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                        call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j+1), fi_kind, i+1, j+1)
                    elseif((i==localncols).and.(j>=2).and.(j<=localnrows-1)) then
                        call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j-1), fi_kind, i-1, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        call coordiToGlobalInd(myid, ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j+1), fi_kind, i-1, j+1)
                        call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                    elseif((j==1).and.(i>=2).and.(i<=localncols-1)) then
                        call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        call coordiToGlobalInd(myid, ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j+1), fi_kind, i-1, j+1)
                        call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j+1), fi_kind, i+1, j+1)
                    elseif((j==localnrows).and.(i>=2).and.(i<=localncols-1)) then
                        call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j-1), fi_kind, i-1, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j-1), fi_kind, i+1, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        call coordiToGlobalInd(myid, ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                    else
                        call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j-1), fi_kind, i-1, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j-1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j-1), fi_kind, i+1, j-1)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        call coordiToGlobalInd(myid, ind_kind, i, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j), fi_kind, i, j)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                        call coordiToGlobalInd(myid, ind_kind, i-1, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j+1), fi_kind, i-1, j+1)
                        call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                        call coordiToGlobalInd(myid, ind_kind, i+1, j+1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j+1), fi_kind, i+1, j+1)
                    end if

                end if

            end do
        end do

    end subroutine dctz_pplap9

    ! '.' equation
    ! '.' element
    ! laplace operation - 5 elements
    subroutine dctz_pplap5(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: field_ind, equ_ind
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        ! the field index
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(field(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call coordiToGlobalInd(myid-1, ind_kind, localncols, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, 1, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, 1, j)
                    elseif((i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call coordiToGlobalInd(myid+1, ind_kind, 1, j, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, localncols, j, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, localncols, j)
                    elseif((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i, 1, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, 1)
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, field_ind)
                        call coordiToGlobalInd(myid, ind_kind, i, localnrows, equ_ind)
                        call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, localnrows)
                    elseif((i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows)) then
                        call coordiToGlobalInd(myid, ind_kind, i, j, field_ind)
                        call setMatValue(field_ind, field_ind, resi(i,j), fi_kind, i, j)
                        if(i /= 1) then
                            call coordiToGlobalInd(myid, ind_kind, i-1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j), fi_kind, i-1, j)
                        end if
                        if(i /= localncols) then
                            call coordiToGlobalInd(myid, ind_kind, i+1, j, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j), fi_kind, i+1, j)
                        end if
                        if(j /= 1) then
                            call coordiToGlobalInd(myid, ind_kind, i, j-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1), fi_kind, i, j-1)
                        end if
                        if(j /= localnrows) then
                            call coordiToGlobalInd(myid, ind_kind, i, j+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1), fi_kind, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine dctz_pplap5

    subroutine setMatValue(col, row, value, fi_kind, eq_i, eq_j)

        integer, intent(in) :: col ! global column index
        integer, intent(in) :: row ! global row index
        real(kind=8), intent(in) :: value
        integer, intent(in) :: fi_kind ! field kind
        integer, intent(in) :: eq_i ! equation x-direction coordinate
        integer, intent(in) :: eq_j ! equation y-direction coordinate

        integer, dimension(:), pointer :: Acols
        integer, dimension(:), pointer :: Arows
        real(kind=8), dimension(:), pointer :: Avalues
        integer, dimension(:), pointer :: AEntryBase
        integer, dimension(:), pointer :: AEntryNum

        integer :: indexr, pos, base, tail, shend, m ,n
        integer :: left, right, mid

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! fi_kind has three digits.
        ! The first digit represents the kind of element:
        ! 1 means the '|' element.
        ! 2 means the '-' element.
        ! 3 means the '.' element.
        ! The last two digits represent the kind of equation:
        ! 51 means the static x-velocity equation.
        ! 52 means the dynamic x-velocity equation.
        ! 53 means the static y-velocity equation.
        ! 54 means the dynamic y-velocity equation.
        ! 55 means the pressure equation.
        ! 56 means the concentration equation.
        ! 57 means the temperature equation.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(fi_kind == 151) then ! AxxStatic
            Acols => AxxCols
            Arows => AxxRows
            Avalues => AxxStaticValues
            AEntryBase => AxxEntryBase
            AEntryNum => AxxEntryNum
        elseif(fi_kind == 152) then ! AxxDyn
            Acols => AxxCols
            Arows => AxxRows
            Avalues => AxxDynValues
            AEntryBase => AxxEntryBase
            AEntryNum => AxxEntryNum
        elseif(fi_kind == 253) then ! AyyStatic
            Acols => AyyCols
            Arows => AyyRows
            Avalues => AyyStaticValues
            AEntryBase => AyyEntryBase
            AEntryNum => AyyEntryNum
        elseif(fi_kind == 254) then ! AyyDyn
            Acols => AyyCols
            Arows => AyyRows
            Avalues => AyyDynValues
            AEntryBase => AyyEntryBase
            AEntryNum => AyyEntryNum
        elseif(fi_kind == 355) then ! Ap
            Acols => ApCols
            Arows => ApRows
            Avalues => ApValues
            AEntryBase => ApEntryBase
            AEntryNum => ApEntryNum
        elseif(fi_kind == 356) then ! Acf
            Acols => AcfCols
            Arows => AcfRows
            Avalues => AcfValues
            AEntryBase => AcfEntryBase
            AEntryNum => AcfEntryNum
        elseif(fi_kind == 357) then ! Atem
            Acols => AtemCols
            Arows => AtemRows
            Avalues => AtemValues
            AEntryBase => AtemEntryBase
            AEntryNum => AtemEntryNum
        end if

        if((fi_kind==151).or.(fi_kind==152)) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            pos = indexr*(eq_j-1) + eq_i
        else
            pos = localncols*(eq_j-1) + eq_i
        end if

        base = AEntryBase(pos)
        tail = base + AEntryNum(pos) - 1

        if(t == 2) then

            do n = base, tail
                if(Acols(n) > col) then
                    do m = n+1, tail
                        if(Acols(m) == 0) then
                            shend = m - 1
                            exit
                        end if
                    end do
                    do m = shend, n, -1
                        Acols(m+1) = Acols(m)
                        Arows(m+1) = Arows(m)
                        Avalues(m+1) = Avalues(m)
                    end do
                    Acols(n) = col
                    Arows(n) = row
                    Avalues(n) = value
                    exit
                elseif(Acols(n) == col) then
                    Avalues(n) = value
                    exit
                elseif(Acols(n) == 0) then
                    Acols(n) = col
                    Arows(n) = row
                    Avalues(n) = value
                    exit
                end if
            end do

        else

            left = base
            right = tail
            mid = (left+right)/2
            do n = 1, AEntryNum(pos)
                if(Acols(mid) == col) then
                    Avalues(mid) = value
                    exit
                elseif(Acols(mid) < col) then
                    left = mid
                    mid = (left+right)/2
                    if((right-left) == 1) then
                        if(Acols(left) == col) then
                            Avalues(left) = value
                            exit
                        elseif(Acols(right) == col) then
                            Avalues(right) = value
                            exit
                        end if
                    end if
                elseif(Acols(mid) > col) then
                    right = mid
                    mid = (left+right)/2
                    if((right-left) == 1) then
                        if(Acols(left) == col) then
                            Avalues(left) = value
                            exit
                        elseif(Acols(right) == col) then
                            Avalues(right) = value
                            exit
                        end if
                    end if
                end if
            end do

        end if

    end subroutine setMatValue

    ! change from the local index to the global index
    subroutine coordiToGlobalInd(pid, ind_kind, local_i, local_j, global_ind)

        integer, intent(in) :: pid
        integer, intent(in) :: ind_kind ! the kind of the index
        integer, intent(in) :: local_i
        integer, intent(in) :: local_j
        integer, intent(out) :: global_ind

        integer :: p_pcol, p_prow
        integer :: xBase, yBase, pBase, uBase

        p_pcol = mod(pid,pncols) + 1
        p_prow = pid/pncols + 1

        ! the base index '|' of the processor for the x-edge
        xBase = (p_prow-1) * localnrows * (nx+1)
        xBase = xBase + (p_pcol-1) * localncols * localnrows

        ! the base index '-' of the processor for the y-edge
        yBase = (p_prow-1) * localnrows * nx
        yBase = yBase + (p_pcol-1) * localncols * localnrows
        if(p_prow == pnrows) then
            yBase = yBase + (p_pcol-1) * localncols
        end if

        ! the base index '.' of the processor for the point at the cell center
        pBase = (p_prow-1) * localnrows * nx
        pBase = pBase + (p_pcol-1) * localncols * localnrows

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ind_kind has three digits.
        ! The first digit represents the kind of element:
        ! 1 means the '|' element.
        ! 2 means the '-' element.
        ! 3 means the '.' element.
        ! The last two digits represent the kind of matrix:
        ! 11 means the velocity matrix.
        ! 12 means the pressure matrix.
        ! 13 means the concentration matrix.
        ! 14 means the temperature matrix.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! the '|' index of the velocity matrix
        if(ind_kind == 111) then
            if(p_pcol == pncols) then
                global_ind = xBase + yBase + (local_j-1)*(localncols+1) + local_i
            else
                global_ind = xBase + yBase + (local_j-1)*localncols + local_i
            end if
        ! the '-' index of the velocity matrix
        elseif(ind_kind == 211) then
            if(p_pcol == pncols) then
                uBase = localnrows * (localncols+1)
            else
                uBase = localnrows * localncols
            end if
            global_ind = xBase + yBase + uBase + (local_j-1)*localncols + local_i
        ! the '.' index of the pressure matrix
        elseif(ind_kind == 312) then
            global_ind = pBase + (local_j-1)*localncols + local_i
        ! the '.' index of the concentration matrix
        elseif(ind_kind == 313) then
            global_ind = pBase + (local_j-1)*localncols + local_i
        ! the '.' index of the temperature matrix
        elseif(ind_kind == 314) then
            global_ind = pBase + (local_j-1)*localncols + local_i
        end if

    end subroutine coordiToGlobalInd

    subroutine genExpField(bx, by, local_nx, local_ny, field, isField)

        integer, intent(in) :: bx, by, local_nx, local_ny
        integer, dimension(:,:), pointer, intent(in out) :: field
        logical, intent(out) :: isField

        field(:,:) = 0

        if((bx>local_nx).or.(by>local_ny)) then
            isField = .false.
        else
            field(bx:local_nx:3, by:local_ny:3) = 1
            isField = .true.
        end if

    end subroutine genExpField

end module DBF_constructMat
