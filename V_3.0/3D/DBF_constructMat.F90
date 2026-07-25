
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

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: field_ind, equ_ind
        integer :: indexl, indexr, indexd, indexu, indexf, indexb, indextemp
        integer :: i, j, k

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
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        if(play /= pnlays) then
            indexb = localnlays + 1
        else
            indexb = localnlays
        end if
        if(pcol /= pncols) then
            indextemp = localncols
        else
            indextemp = localncols + 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then

                        if((i==0).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid-1, ind_kind, localncols, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, 1, j, k)
                        elseif((pcol/=pncols).and.(i==localncols+1).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and. &!
                            (k/=localnlays+1)) then
                            call coordiToGlobalInd(myid+1, ind_kind, 1, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, localncols, j, k)
                        elseif((j==0).and.(i/=0).and.(i/=indextemp+1).and.(k/=0).and.(k/=localnlays+1).and.(.not.((pcol==1).and. &!
                            (i==1).and.(isDiriX0_p(ylower+j,zlower+k-1)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                            (isDiriX1_p(ylower+j,zlower+k-1)==0)))) then
                            call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, 1, k)
                        elseif((j==localnrows+1).and.(i/=0).and.(i/=indextemp+1).and.(k/=0).and.(k/=localnlays+1).and.(.not. &!
                            ((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-2,zlower+k-1)==0))).and.(.not.((pcol==pncols).and. &!
                            (i==indextemp).and.(isDiriX1_p(ylower+j-2,zlower+k-1)==0)))) then
                            call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, localnrows, k)
                        elseif((k==0).and.(i/=0).and.(i/=indextemp+1).and.(j/=0).and.(j/=localnrows+1).and.(.not.((pcol==1).and. &!
                            (i==1).and.(isDiriX0_p(ylower+j-1,zlower+k)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                            (isDiriX1_p(ylower+j-1,zlower+k)==0)))) then
                            call coordiToGlobalInd(myid-pncols*pnrows, ind_kind, i, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, 1)
                        elseif((k==localnlays+1).and.(i/=0).and.(i/=indextemp+1).and.(j/=0).and.(j/=localnrows+1).and.(.not. &!
                            ((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-2)==0))).and.(.not.((pcol==pncols).and. &!
                            (i==indextemp).and.(isDiriX1_p(ylower+j-1,zlower+k-2)==0)))) then
                            call coordiToGlobalInd(myid+pncols*pnrows, ind_kind, i, j, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, localnlays)
                        elseif((i>=1).and.(i<=indextemp).and.(j>=1).and.(j<=localnrows).and.(k>=1).and.(k<=localnlays)) then
                            call coordiToGlobalInd(myid, ind_kind, i, j, k, field_ind)
                            if((i>=3).or.((i==2).and.(pcol/=1)).or.((i==2).and.(pcol==1).and. &!
                                (isDiriX0_p(ylower+j-1,zlower+k-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, i-1, j, k)
                            end if
                            if((i<=indextemp-2).or.((i==indextemp-1).and.(pcol/=pncols)).or.((i==indextemp-1).and. &!
                                (pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, i+1, j, k)
                            end if
                            if(((j/=1).and.(i>=2).and.(i<=indextemp-1)).or.((j/=1).and.(i==1).and.(pcol/=1)).or. &!
                                ((j/=1).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-2,zlower+k-1)/=0)).or.((j/=1).and. &!
                                (i==indextemp).and.(pcol/=pncols)).or.((j/=1).and.(i==indextemp).and.(pcol==pncols).and. &!
                                (isDiriX1_p(ylower+j-2,zlower+k-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, j-1, k)
                            end if
                            if(((j/=localnrows).and.(i>=2).and.(i<=indextemp-1)).or.((j/=localnrows).and.(i==1).and. &!
                                (pcol/=1)).or.((j/=localnrows).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j,zlower+k-1)/=0)) &!
                                .or.((j/=localnrows).and.(i==indextemp).and.(pcol/=pncols)).or.((j/=localnrows).and. &!
                                (i==indextemp).and.(pcol==pncols).and.(isDiriX1_p(ylower+j,zlower+k-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, j+1, k)
                            end if
                            if(((k/=1).and.(i>=2).and.(i<=indextemp-1)).or.((k/=1).and.(i==1).and.(pcol/=1)).or. &!
                                ((k/=1).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k-2)/=0)).or.((k/=1).and. &!
                                (i==indextemp).and.(pcol/=pncols)).or.((k/=1).and.(i==indextemp).and.(pcol==pncols).and. &!
                                (isDiriX1_p(ylower+j-1,zlower+k-2)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, k-1)
                            end if
                            if(((k/=localnlays).and.(i>=2).and.(i<=indextemp-1)).or.((k/=localnlays).and.(i==1).and. &!
                                (pcol/=1)).or.((k/=localnlays).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k)/=0)) &!
                                .or.((k/=localnlays).and.(i==indextemp).and.(pcol/=pncols)).or.((k/=localnlays).and. &!
                                (i==indextemp).and.(pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, k+1)
                            end if
                            call setMatValue(field_ind, field_ind, resi(i,j,k), fi_kind, i, j, k)
                        end if

                    end if

                end do
            end do
        end do

    end subroutine dctz_xxlap

    ! '|' equation
    ! '.' element
    ! gradient operation
    subroutine dctz_xpgra(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        ! notice that the dimensions of field are different from the dimensions of resi
        ! in the function
        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: field_ind, equ_ind
        integer :: i, j, k

        ! the field index
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then
                        if(i == 0) then
                            call coordiToGlobalInd(myid-1, f_ind_kind, localncols, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, 1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, 1, j, k)
                        elseif((i==1).and.(pcol==1)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            if(isDiriX0_p(ylower+j-1,zlower+k-1) /= 0) then
                                call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                            end if
                            call coordiToGlobalInd(myid, e_ind_kind, i+1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, i+1, j, k)
                        elseif((i==localncols).and.(pcol/=pncols)) then
                            call coordiToGlobalInd(myid, f_ind_kind, localncols, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, localncols, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, localncols, j, k)
                        elseif((i==localncols).and.(pcol==pncols)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                            if(isDiriX1_p(ylower+j-1,zlower+k-1) /= 0) then
                                call coordiToGlobalInd(myid, e_ind_kind, i+1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, i+1, j, k)
                            end if
                        else
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                            call coordiToGlobalInd(myid, e_ind_kind, i+1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, i+1, j, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine dctz_xpgra

    ! '-' equation
    ! '-' element
    ! laplace operation
    subroutine dctz_yylap(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: indexl, indexr, indexd, indexu, indexf, indexb, indextemp
        integer :: field_ind, equ_ind
        integer :: i, j, k

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
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        if(play /= pnlays) then
            indexb = localnlays + 1
        else
            indexb = localnlays
        end if
        if(prow /= pnrows) then
            indextemp = localnrows
        else
            indextemp = localnrows + 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then

                        if((i==0).and.(j/=0).and.(j/=indextemp+1).and.(k/=0).and.(k/=localnlays+1).and.(.not.((prow==1) &!
                            .and.(j==1).and.(isDiriY0_p(xlower+i,zlower+k-1)==0))).and.(.not.((prow==pnrows).and.(j==indextemp) &!
                            .and.(isDiriY1_p(xlower+i,zlower+k-1)==0)))) then
                            call coordiToGlobalInd(myid-1, ind_kind, localncols, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, 1, j, k)
                        elseif((i==localncols+1).and.(j/=0).and.(j/=indextemp+1).and.(k/=0).and.(k/=localnlays+1).and. &!
                            (.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-2,zlower+k-1)==0))).and.(.not.((prow==pnrows) &!
                            .and.(j==indextemp).and.(isDiriY1_p(xlower+i-2,zlower+k-1)==0)))) then
                            call coordiToGlobalInd(myid+1, ind_kind, 1, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, localncols, j, k)
                        elseif((j==0).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, 1, k)
                        elseif((prow/=pnrows).and.(j==localnrows+1).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and. &!
                            (k/=localnlays+1)) then
                            call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, localnrows, k)
                        elseif((k==0).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=indextemp+1).and.(.not.((prow==1) &!
                            .and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k)==0))).and.(.not.((prow==pnrows).and. &!
                            (j==indextemp).and. (isDiriY1_p(xlower+i-1,zlower+k)==0)))) then
                            call coordiToGlobalInd(myid-pncols*pnrows, ind_kind, i, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, 1)
                        elseif((k==localnlays+1).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=indextemp+1).and. &!
                            (.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-2)==0))).and.(.not. &!
                            ((prow==pnrows).and.(j==indextemp).and.(isDiriY1_p(xlower+i-1,zlower+k-2)==0)))) then
                            call coordiToGlobalInd(myid+pncols*pnrows, ind_kind, i, j, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, localnlays)
                        elseif((j>=1).and.(j<=indextemp).and.(i>=1).and.(i<=localncols).and.(k>=1).and.(k<=localnlays)) then
                            call coordiToGlobalInd(myid, ind_kind, i, j, k, field_ind)
                            if(((i/=1).and.(j>=2).and.(j<=indextemp-1)).or.((i/=1).and.(j==1).and.(prow/=1)).or.((i/=1) &!
                                .and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-2,zlower+k-1)/=0)).or.((i/=1).and. &!
                                (j==indextemp).and.(prow/=pnrows)).or.((i/=1).and.(j==indextemp).and.(prow==pnrows).and. &!
                                (isDiriY1_p(xlower+i-2,zlower+k-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, i-1, j, k)
                            end if
                            if(((i/=localncols).and.(j>=2).and.(j<=indextemp-1)).or.((i/=localncols).and.(j==1).and. &!
                                (prow/=1)).or.((i/=localncols).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i,zlower+k-1)/=0)) &!
                                .or.((i/=localncols).and.(j==indextemp).and.(prow/=pnrows)).or.((i/=localncols).and.(j==indextemp) &!
                                .and.(prow==pnrows).and.(isDiriY1_p(xlower+i,zlower+k-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, i+1, j, k)
                            end if
                            if((j>=3).or.((j==2).and.(prow/=1)).or.((j==2).and.(prow==1).and. &!
                                (isDiriY0_p(xlower+i-1,zlower+k-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, j-1, k)
                            end if
                            if((j<=indextemp-2).or.((j==indextemp-1).and.(prow/=pnrows)).or.((j==indextemp-1).and. &!
                                (prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, j+1, k)
                            end if
                            if(((k/=1).and.(j>=2).and.(j<=indextemp-1)).or.((k/=1).and.(j==1).and.(prow/=1)).or. &!
                                ((k/=1).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k-2)/=0)).or.((k/=1).and. &!
                                (j==indextemp).and.(prow/=pnrows)).or.((k/=1).and.(j==indextemp).and.(prow==pnrows).and. &!
                                (isDiriY1_p(xlower+i-1,zlower+k-2)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, k-1)
                            end if
                            if(((k/=localnlays).and.(j>=2).and.(j<=indextemp-1)).or.((k/=localnlays).and.(j==1).and. &!
                                (prow/=1)).or.((k/=localnlays).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k)/=0)) &!
                                .or.((k/=localnlays).and.(j==indextemp).and.(prow/=pnrows)).or.((k/=localnlays).and. &!
                                (j==indextemp).and.(prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, k+1)
                            end if
                            call setMatValue(field_ind, field_ind, resi(i,j,k), fi_kind, i, j, k)
                        end if

                    end if

                end do
            end do
        end do

    end subroutine dctz_yylap

    ! '-' equation
    ! '.' element
    ! gradient operation
    subroutine dctz_ypgra(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: field_ind, equ_ind
        integer :: i, j, k

        ! the field index
        indexl = 1
        indexr = localncols
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows
        indexf = 1
        indexb = localnlays

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then
                        if(j == 0) then
                            call coordiToGlobalInd(myid-pncols, f_ind_kind, i, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, 1, k)
                        elseif((j==1).and.(prow==1)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            if(isDiriY0_p(xlower+i-1,zlower+k-1) /= 0) then
                                call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                            end if
                            call coordiToGlobalInd(myid, e_ind_kind, i, j+1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, j+1, k)
                        elseif((j==localnrows).and.(prow/=pnrows)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, localnrows, k)
                        elseif((j==localnrows).and.(prow==pnrows)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, localnrows, k)
                            if(isDiriY1_p(xlower+i-1,zlower+k-1) /= 0) then
                                call coordiToGlobalInd(myid, e_ind_kind, i, localnrows+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, localnrows+1, k)
                            end if
                        else
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j+1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, j+1, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine dctz_ypgra

    ! '|_' equation
    ! '|_' element
    ! laplace operation
    subroutine dctz_zzlap(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: field_ind, equ_ind
        integer :: indexl, indexr, indexd, indexu, indexf, indexb, indextemp
        integer :: i, j, k

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
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        indexb = localnlays + 1
        
        if(play /= pnlays) then
            indextemp = localnlays
        else
            indextemp = localnlays + 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then

                        if((i==0).and.(k/=0).and.(k/=indextemp+1).and.(j/=0).and.(j/=localnrows+1).and.(.not.((play==1).and. &!
                            (k==1).and.(isDiriZ0_p(xlower+i,ylower+j-1)==0))).and.(.not.((play==pnlays).and.(k==indextemp).and. &!
                            (isDiriZ1_p(xlower+i,ylower+j-1)==0)))) then
                            call coordiToGlobalInd(myid-1, ind_kind, localncols, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, 1, j, k)
                        elseif((i==localncols+1).and.(k/=0).and.(k/=indextemp+1).and.(j/=0).and.(j/=localnrows+1).and. &!
                            (.not.((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-2,ylower+j-1)==0))).and.(.not.((play==pnlays) &!
                            .and.(k==indextemp).and.(isDiriZ1_p(xlower+i-2,ylower+j-1)==0)))) then
                            call coordiToGlobalInd(myid+1, ind_kind, 1, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, localncols, j, k)
                        elseif((j==0).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=indextemp+1).and.(.not.((play==1) &!
                            .and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j)==0))).and.(.not.((play==pnlays).and.(k==indextemp) &!
                            .and.(isDiriZ1_p(xlower+i-1,ylower+j)==0)))) then
                            call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, 1, k)
                        elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=indextemp+1).and. &!
                            (.not.((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-2)==0))).and.(.not.((play==pnlays) &!
                            .and.(k==indextemp).and.(isDiriZ1_p(xlower+i-1,ylower+j-2)==0)))) then
                            call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, localnrows, k)
                        elseif((k==0).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call coordiToGlobalInd(myid-pncols*pnrows, ind_kind, i, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, 1)
                        elseif((play/=pnlays).and.(k==localnlays+1).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and. &!
                            (j/=localnrows+1)) then
                            call coordiToGlobalInd(myid+pncols*pnrows, ind_kind, i, j, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, localnlays)
                        elseif((k>=1).and.(k<=indextemp).and.(i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows)) then
                            call coordiToGlobalInd(myid, ind_kind, i, j, k, field_ind)
                            if(((i/=1).and.(k>=2).and.(k<=indextemp-1)).or.((i/=1).and.(k==1).and.(play/=1)).or.((i/=1).and. &!
                                (k==1).and.(play==1).and.(isDiriZ0_p(xlower+i-2,ylower+j-1)/=0)).or.((i/=1).and.(k==indextemp) &!
                                .and.(play/=pnlays)).or.((i/=1).and.(k==indextemp).and.(play==pnlays).and. &!
                                (isDiriZ1_p(xlower+i-2,ylower+j-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, i-1, j, k)
                            end if
                            if(((i/=localncols).and.(k>=2).and.(k<=indextemp-1)).or.((i/=localncols).and.(k==1).and.(play/=1)) &!
                                .or.((i/=localncols).and.(k==1).and.(play==1).and.(isDiriZ0_p(xlower+i,ylower+j-1)/=0)).or. &!
                                ((i/=localncols).and.(k==indextemp).and.(play/=pnlays)).or.((i/=localncols).and.(k==indextemp) &!
                                .and.(play==pnlays).and.(isDiriZ1_p(xlower+i,ylower+j-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, i+1, j, k)
                            end if
                            if(((j/=1).and.(k>=2).and.(k<=indextemp-1)).or.((j/=1).and.(k==1).and.(play/=1)).or. &!
                                ((j/=1).and.(k==1).and.(play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-2)/=0)).or.((j/=1).and. &!
                                (k==indextemp).and.(play/=pnlays)).or.((j/=1).and.(k==indextemp).and.(play==pnlays).and. &!
                                (isDiriZ1_p(xlower+i-1,ylower+j-2)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, j-1, k)
                            end if
                            if(((j/=localnrows).and.(k>=2).and.(k<=indextemp-1)).or.((j/=localnrows).and.(k==1).and. &!
                                (play/=1)).or.((j/=localnrows).and.(k==1).and.(play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j)/=0)) &!
                                .or.((j/=localnrows).and.(k==indextemp).and.(play/=pnlays)).or.((j/=localnrows).and. &!
                                (k==indextemp).and.(play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, j+1, k)
                            end if
                            if((k>=3).or.((k==2).and.(play/=1)).or.((k==2).and.(play==1).and. &!
                                (isDiriZ0_p(xlower+i-1,ylower+j-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, k-1)
                            end if
                            if((k<=indextemp-2).or.((k==indextemp-1).and.(play/=pnlays)).or.((k==indextemp-1).and. &!
                                (play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)/=0))) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, k+1)
                            end if
                            call setMatValue(field_ind, field_ind, resi(i,j,k), fi_kind, i, j, k)
                        end if

                    end if

                end do
            end do
        end do

    end subroutine dctz_zzlap

    ! '|_' equation
    ! '.' element
    ! gradient operation
    subroutine dctz_zpgra(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: field_ind, equ_ind
        integer :: i, j, k

        ! the field index
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        indexb = localnlays
        
        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then
                        if(k == 0) then
                            call coordiToGlobalInd(myid-pncols*pnrows, f_ind_kind, i, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, 1)
                        elseif((k==1).and.(play==1)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            if(isDiriZ0_p(xlower+i-1,ylower+j-1) /= 0) then
                                call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                            end if
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, k+1)
                        elseif((k==localnlays).and.(play/=pnlays)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, localnlays)
                        elseif((k==localnlays).and.(play==pnlays)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, localnlays)
                            if(isDiriZ1_p(xlower+i-1,ylower+j-1) /= 0) then
                                call coordiToGlobalInd(myid, e_ind_kind, i, j, localnlays+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, localnlays+1)
                            end if
                        else
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k+1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, k+1)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine dctz_zpgra

    ! '.' equation
    ! '|' element
    ! divergence operation
    subroutine dctz_pxdiv(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: field_ind, equ_ind
        integer :: i, j, k

        ! the field index
        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then
                        if(i == 1) then
                            call coordiToGlobalInd(myid, f_ind_kind, 1, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, 1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, 1, j, k)
                        elseif((i==localncols+1).and.(pcol/=pncols)) then
                            call coordiToGlobalInd(myid+1, f_ind_kind, 1, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, localncols, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, localncols, j, k)
                        elseif((i==localncols+1).and.(pcol==pncols)) then
                            call coordiToGlobalInd(myid, f_ind_kind, localncols+1, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, localncols, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, localncols, j, k)
                        else
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i-1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, i-1, j, k)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine dctz_pxdiv

    ! '.' equation
    ! '-' element
    ! divergence operation
    subroutine dctz_pydiv(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: field_ind, equ_ind
        integer :: i, j, k

        ! the field index
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1
        indexf = 1
        indexb = localnlays

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then
                        if(j == 1) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, 1, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, 1, k)
                        elseif((j==localnrows+1).and.(prow/=pnrows)) then
                            call coordiToGlobalInd(myid+pncols, f_ind_kind, i, 1, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, localnrows, k)
                        elseif((j==localnrows+1).and.(prow==pnrows)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, localnrows+1, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, localnrows, k)
                        else
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j-1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, j-1, k)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine dctz_pydiv

    ! '.' equation
    ! '|_' element
    ! divergence operation
    subroutine dctz_pzdiv(field, resi, f_ind_kind, e_ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: f_ind_kind, e_ind_kind, fi_kind
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: field_ind, equ_ind
        integer :: i, j, k

        ! the field index
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays + 1

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then
                        if(k == 1) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, 1, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, 1)
                        elseif((k==localnlays+1).and.(play/=pnlays)) then
                            call coordiToGlobalInd(myid+pncols*pnrows, f_ind_kind, i, j, 1, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, localnlays)
                        elseif((k==localnlays+1).and.(play==pnlays)) then
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, localnlays+1, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, localnlays)
                        else
                            call coordiToGlobalInd(myid, f_ind_kind, i, j, k, field_ind)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k-1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, k-1)
                            call coordiToGlobalInd(myid, e_ind_kind, i, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k), fi_kind, i, j, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine dctz_pzdiv

    ! '.' equation
    ! '.' element
    ! derivative operation
    subroutine dctz_ppder(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind
        
        integer :: field_ind
        integer :: i, j, k

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    if(field(i,j,k) == 1) then
                        call coordiToGlobalInd(myid, ind_kind, i, j, k, field_ind)
                        call setMatValue(field_ind, field_ind, resi(i,j,k), fi_kind, i, j, k)
                    end if

                end do
            end do
        end do

    end subroutine dctz_ppder

    ! '.' equation
    ! '.' element
    ! laplace operation - 19 elements
    subroutine dctz_pplap19(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: field_ind, equ_ind
        integer :: i, j, k

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
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        if(play /= pnlays) then
            indexb = localnlays + 1
        else
            indexb = localnlays
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then

                        if((i==0).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid-1, ind_kind, localncols, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(1,j,k), fi_kind, 1, j, k)
                            if(j/=1) then
                                call coordiToGlobalInd(myid, ind_kind, 1, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(1,j-1,k), fi_kind, 1, j-1, k)
                            end if
                            if(j/=localnrows) then
                                call coordiToGlobalInd(myid, ind_kind, 1, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(1,j+1,k), fi_kind, 1, j+1, k)
                            end if
                            if(k/=1) then
                                call coordiToGlobalInd(myid, ind_kind, 1, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(1,j,k-1), fi_kind, 1, j, k-1)
                            end if
                            if(k/=localnlays) then
                                call coordiToGlobalInd(myid, ind_kind, 1, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(1,j,k+1), fi_kind, 1, j, k+1)
                            end if
                        end if

                        if((i==localncols+1).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid+1, ind_kind, 1, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(localncols,j,k), fi_kind, localncols, j, k)
                            if(j/=1) then
                                call coordiToGlobalInd(myid, ind_kind, localncols, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(localncols,j-1,k), fi_kind, localncols, j-1, k)
                            end if
                            if(j/=localnrows) then
                                call coordiToGlobalInd(myid, ind_kind, localncols, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(localncols,j+1,k), fi_kind, localncols, j+1, k)
                            end if
                            if(k/=1) then
                                call coordiToGlobalInd(myid, ind_kind, localncols, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(localncols,j,k-1), fi_kind, localncols, j, k-1)
                            end if
                            if(k/=localnlays) then
                                call coordiToGlobalInd(myid, ind_kind, localncols, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(localncols,j,k+1), fi_kind, localncols, j, k+1)
                            end if
                        end if

                        if((j==0).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,1,k), fi_kind, i, 1, k)
                            if(i/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, 1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,1,k), fi_kind, i-1, 1, k)
                            end if
                            if(i/=localncols) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, 1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,1,k), fi_kind, i+1, 1, k)
                            end if
                            if(k/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i, 1, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,1,k-1), fi_kind, i, 1, k-1)
                            end if
                            if(k/=localnlays) then
                                call coordiToGlobalInd(myid, ind_kind, i, 1, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,1,k+1), fi_kind, i, 1, k+1)
                            end if
                        end if

                        if((j==localnrows+1).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,localnrows,k), fi_kind, i, localnrows, k)
                            if(i/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, localnrows, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,localnrows,k), fi_kind, i-1, localnrows, k)
                            end if
                            if(i/=localncols) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, localnrows, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,localnrows,k), fi_kind, i+1, localnrows, k)
                            end if
                            if(k/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i, localnrows, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,localnrows,k-1), fi_kind, i, localnrows, k-1)
                            end if
                            if(k/=localnlays) then
                                call coordiToGlobalInd(myid, ind_kind, i, localnrows, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,localnrows,k+1), fi_kind, i, localnrows, k+1)
                            end if
                        end if

                        if((k==0).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call coordiToGlobalInd(myid-pncols*pnrows, ind_kind, i, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,1), fi_kind, i, j, 1)
                            if(i/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, 1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,1), fi_kind, i-1, j, 1)
                            end if
                            if(i/=localncols) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, 1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,1), fi_kind, i+1, j, 1)
                            end if
                            if(j/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, 1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,1), fi_kind, i, j-1, 1)
                            end if
                            if(j/=localnrows) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, 1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,1), fi_kind, i, j+1, 1)
                            end if
                        end if

                        if((k==localnlays+1).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call coordiToGlobalInd(myid+pncols*pnrows, ind_kind, i, j, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,localnlays), fi_kind, i, j, localnlays)
                            if(i/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, localnlays, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,localnlays), fi_kind, i-1, j, localnlays)
                            end if
                            if(i/=localncols) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, localnlays, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,localnlays), fi_kind, i+1, j, localnlays)
                            end if
                            if(j/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, localnlays, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,localnlays), fi_kind, i, j-1, localnlays)
                            end if
                            if(j/=localnrows) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, localnlays, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,localnlays), fi_kind, i, j+1, localnlays)
                            end if
                        end if

                        if((i==0).and.(j==0).and.(k>=1).and.(k<=localnlays)) then
                            call coordiToGlobalInd(myid-1-pncols, ind_kind, localncols, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(1,1,k), fi_kind, 1, 1, k)
                        end if

                        if((i==0).and.(j==localnrows+1).and.(k>=1).and.(k<=localnlays)) then
                            call coordiToGlobalInd(myid-1+pncols, ind_kind, localncols, 1, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(1,localnrows,k), fi_kind, 1, localnrows, k)
                        end if

                        if((i==0).and.(k==0).and.(j>=1).and.(j<=localnrows)) then
                            call coordiToGlobalInd(myid-1-pncols*pnrows, ind_kind, localncols, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(1,j,1), fi_kind, 1, j, 1)
                        end if

                        if((i==0).and.(k==localnlays+1).and.(j>=1).and.(j<=localnrows)) then
                            call coordiToGlobalInd(myid-1+pncols*pnrows, ind_kind, localncols, j, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(1,j,localnlays), fi_kind, 1, j, localnlays)
                        end if

                        if((i==localncols+1).and.(j==0).and.(k>=1).and.(k<=localnlays)) then
                            call coordiToGlobalInd(myid+1-pncols, ind_kind, 1, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(localncols,1,k), fi_kind, localncols, 1, k)
                        end if

                        if((i==localncols+1).and.(j==localnrows+1).and.(k>=1).and.(k<=localnlays)) then
                            call coordiToGlobalInd(myid+1+pncols, ind_kind, 1, 1, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(localncols,localnrows,k), fi_kind, localncols, localnrows, k)
                        end if

                        if((i==localncols+1).and.(k==0).and.(j>=1).and.(j<=localnrows)) then
                            call coordiToGlobalInd(myid+1-pncols*pnrows, ind_kind, 1, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(localncols,j,1), fi_kind, localncols, j, 1)
                        end if

                        if((i==localncols+1).and.(k==localnlays+1).and.(j>=1).and.(j<=localnrows)) then
                            call coordiToGlobalInd(myid+1+pncols*pnrows, ind_kind, 1, j, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(localncols,j,localnlays), fi_kind, localncols, j, localnlays)
                        end if

                        if((j==0).and.(k==0).and.(i>=1).and.(i<=localncols)) then
                            call coordiToGlobalInd(myid-pncols-pncols*pnrows, ind_kind, i, localnrows, &!
                                localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, 1, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,1,1), fi_kind, i, 1, 1)
                        end if

                        if((j==0).and.(k==localnlays+1).and.(i>=1).and.(i<=localncols)) then
                            call coordiToGlobalInd(myid-pncols+pncols*pnrows, ind_kind, i, localnrows, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, 1, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,1,localnlays), fi_kind, i, 1, localnlays)
                        end if

                        if((j==localnrows+1).and.(k==0).and.(i>=1).and.(i<=localncols)) then
                            call coordiToGlobalInd(myid+pncols-pncols*pnrows, ind_kind, i, 1, localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, localnrows, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,localnrows,1), fi_kind, i, localnrows, 1)
                        end if

                        if((j==localnrows+1).and.(k==localnlays+1).and.(i>=1).and.(i<=localncols)) then
                            call coordiToGlobalInd(myid+pncols+pncols*pnrows, ind_kind, i, 1, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, localnrows, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,localnrows,localnlays), fi_kind, i, localnrows, localnlays)
                        end if

                        if((i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows).and.(k>=1).and.(k<=localnlays)) then
                            call coordiToGlobalInd(myid, ind_kind, i, j, k, field_ind)
                            call setMatValue(field_ind, field_ind, resi(i,j,k), fi_kind, i, j, k)
                            if(i/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, i-1, j, k)
                            end if
                            if(i/=localncols) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, i+1, j, k)
                            end if
                            if(j/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, j-1, k)
                            end if
                            if(j/=localnrows) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, j+1, k)
                            end if
                            if(k/=1) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, k-1)
                            end if
                            if(k/=localnlays) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, k+1)
                            end if
                            if((i/=1).and.(j/=1)) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j-1,k), fi_kind, i-1, j-1, k)
                            end if
                            if((i/=1).and.(j/=localnrows)) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j+1,k), fi_kind, i-1, j+1, k)
                            end if
                            if((i/=1).and.(k/=1)) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,k-1), fi_kind, i-1, j, k-1)
                            end if
                            if((i/=1).and.(k/=localnlays)) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,k+1), fi_kind, i-1, j, k+1)
                            end if
                            if((i/=localncols).and.(j/=1)) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j-1,k), fi_kind, i+1, j-1, k)
                            end if
                            if((i/=localncols).and.(j/=localnrows)) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j+1,k), fi_kind, i+1, j+1, k)
                            end if
                            if((i/=localncols).and.(k/=1)) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,k-1), fi_kind, i+1, j, k-1)
                            end if
                            if((i/=localncols).and.(k/=localnlays)) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,k+1), fi_kind, i+1, j, k+1)
                            end if
                            if((j/=1).and.(k/=1)) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,k-1), fi_kind, i, j-1, k-1)
                            end if
                            if((j/=1).and.(k/=localnlays)) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,k+1), fi_kind, i, j-1, k+1)
                            end if
                            if((j/=localnrows).and.(k/=1)) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,k-1), fi_kind, i, j+1, k-1)
                            end if
                            if((j/=localnrows).and.(k/=localnlays)) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,k+1), fi_kind, i, j+1, k+1)
                            end if
                        end if

                    end if

                end do
            end do
        end do

    end subroutine dctz_pplap19

    ! '.' equation
    ! '.' element
    ! laplace operation - 7 elements
    subroutine dctz_pplap7(field, resi, ind_kind, fi_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: field
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: ind_kind, fi_kind

        integer :: field_ind, equ_ind
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: i, j, k

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
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        if(play /= pnlays) then
            indexb = localnlays + 1
        else
            indexb = localnlays
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if(field(i,j,k) == 1) then

                        if((i==0).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid-1, ind_kind, localncols, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, 1, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, 1, j, k)
                        elseif((i==localncols+1).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and. &!
                            (k/=localnlays+1)) then
                            call coordiToGlobalInd(myid+1, ind_kind, 1, j, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, localncols, j, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, localncols, j, k)
                        elseif((j==0).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid-pncols, ind_kind, i, localnrows, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, 1, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, 1, k)
                        elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call coordiToGlobalInd(myid+pncols, ind_kind, i, 1, k, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, localnrows, k, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, localnrows, k)
                        elseif((k==0).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call coordiToGlobalInd(myid-pncols*pnrows, ind_kind, i, j, localnlays, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, 1, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, 1)
                        elseif((k==localnlays+1).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call coordiToGlobalInd(myid+pncols*pnrows, ind_kind, i, j, 1, field_ind)
                            call coordiToGlobalInd(myid, ind_kind, i, j, localnlays, equ_ind)
                            call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, localnlays)
                        elseif((i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows).and.(k>=1).and.(k<=localnlays)) then
                            call coordiToGlobalInd(myid, ind_kind, i, j, k, field_ind)
                            call setMatValue(field_ind, field_ind, resi(i,j,k), fi_kind, i, j, k)
                            if(i /= 1) then
                                call coordiToGlobalInd(myid, ind_kind, i-1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i-1,j,k), fi_kind, i-1, j, k)
                            end if
                            if(i /= localncols) then
                                call coordiToGlobalInd(myid, ind_kind, i+1, j, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i+1,j,k), fi_kind, i+1, j, k)
                            end if
                            if(j /= 1) then
                                call coordiToGlobalInd(myid, ind_kind, i, j-1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j-1,k), fi_kind, i, j-1, k)
                            end if
                            if(j /= localnrows) then
                                call coordiToGlobalInd(myid, ind_kind, i, j+1, k, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j+1,k), fi_kind, i, j+1, k)
                            end if
                            if(k /= 1) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k-1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k-1), fi_kind, i, j, k-1)
                            end if
                            if(k /= localnlays) then
                                call coordiToGlobalInd(myid, ind_kind, i, j, k+1, equ_ind)
                                call setMatValue(field_ind, equ_ind, resi(i,j,k+1), fi_kind, i, j, k+1)
                            end if
                        end if

                    end if

                end do
            end do
        end do

    end subroutine dctz_pplap7

    subroutine setMatValue(col, row, value, fi_kind, eq_i, eq_j, eq_k)

        integer, intent(in) :: col ! global column index
        integer, intent(in) :: row ! global row index
        real(kind=8), intent(in) :: value
        integer, intent(in) :: fi_kind ! field kind
        integer, intent(in) :: eq_i ! equation x-direction coordinate
        integer, intent(in) :: eq_j ! equation y-direction coordinate
        integer, intent(in) :: eq_k ! equation z-direction coordinate

        integer, dimension(:), pointer :: Acols
        integer, dimension(:), pointer :: Arows
        real(kind=8), dimension(:), pointer :: Avalues
        integer, dimension(:), pointer :: AEntryBase
        integer, dimension(:), pointer :: AEntryNum

        integer :: pos, base, tail, shend, m ,n
        integer :: left, right, mid

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! fi_kind has three digits.
        ! The first digit represents the kind of element:
        ! 1 means the '|' element.
        ! 2 means the '-' element.
        ! 3 means the '|_' element.
        ! 4 means the '.' element.
        ! The last two digits represent the kind of equation:
        ! 51 means the static x-velocity equation.
        ! 52 means the dynamic x-velocity equation.
        ! 53 means the static y-velocity equation.
        ! 54 means the dynamic y-velocity equation.
        ! 55 means the static z-velocity equation.
        ! 56 means the dynamic z-velocity equation.
        ! 57 means the pressure equation.
        ! 58 means the concentration equation.
        ! 59 means the temperature equation.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! matrix kind
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
        elseif(fi_kind == 355) then ! AzzStatic
            Acols => AzzCols
            Arows => AzzRows
            Avalues => AzzStaticValues
            AEntryBase => AzzEntryBase
            AEntryNum => AzzEntryNum
        elseif(fi_kind == 356) then ! AzzDyn
            Acols => AzzCols
            Arows => AzzRows
            Avalues => AzzDynValues
            AEntryBase => AzzEntryBase
            AEntryNum => AzzEntryNum
        elseif(fi_kind == 457) then ! Ap
            Acols => ApCols
            Arows => ApRows
            Avalues => ApValues
            AEntryBase => ApEntryBase
            AEntryNum => ApEntryNum
        elseif(fi_kind == 458) then ! Acf
            Acols => AcfCols
            Arows => AcfRows
            Avalues => AcfValues
            AEntryBase => AcfEntryBase
            AEntryNum => AcfEntryNum
        elseif(fi_kind == 459) then ! Atem
            Acols => AtemCols
            Arows => AtemRows
            Avalues => AtemValues
            AEntryBase => AtemEntryBase
            AEntryNum => AtemEntryNum
        end if

        call coordiToLocalInd(fi_kind, eq_i, eq_j, eq_k, pos)

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

    ! change from the local coordinate to the global index
    subroutine coordiToGlobalInd(pid, ind_kind, local_i, local_j, local_k, global_ind)

        integer, intent(in) :: pid
        integer, intent(in) :: ind_kind
        integer, intent(in) :: local_i, local_j, local_k
        integer, intent(out) :: global_ind

        integer :: p_pcol, p_prow, p_play
        integer :: xbase, yBase, zBase, pBase, uBase, vBase

        p_play = pid/(pnrows*pncols)+1
        p_prow = (pid-(p_play-1)*pnrows*pncols)/pncols+1
        p_pcol = (pid-(p_play-1)*pnrows*pncols)-(p_prow-1)*pncols+1

        ! the base index '|' of the processor for the x-edge
        xBase = (p_play-1)*(nx+1)*ny*localnlays
        xBase = xBase + (p_prow-1)*(nx+1)*localnrows*localnlays
        xBase = xBase + (p_pcol-1)*localncols*localnrows*localnlays

        ! the base index '-' of the processor for the y-edge
        yBase = (p_play-1)*nx*(ny+1)*localnlays
        yBase = yBase + (p_prow-1)*nx*localnrows*localnlays
        yBase = yBase + (p_pcol-1)*localncols*localnrows*localnlays
        if(p_prow == pnrows) then
            yBase = yBase + (p_pcol-1)*localncols*localnlays
        end if

        ! the base index '|_' of the processor for the z-edge
        zBase = (p_play-1)*nx*ny*localnlays
        zBase = zBase + (p_prow-1)*nx*localnrows*localnlays
        zBase = zBase + (p_pcol-1)*localncols*localnrows*localnlays
        if(p_play == pnlays) then
            zBase = zBase + (p_prow-1)*nx*localnrows + (p_pcol-1)*localncols*localnrows
        end if

        ! the base index '.' of the processor for the point at the cell center
        pBase = (p_play-1)*nx*ny*localnlays
        pBase = pBase + (p_prow-1)*nx*localnrows*localnlays
        pBase = pBase + (p_pcol-1)*localncols*localnrows*localnlays

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ind_kind has three digits.
        ! The first digit represents the kind of element:
        ! 1 means the '|' element.
        ! 2 means the '-' element.
        ! 3 means the '|_' element.
        ! 4 means the '.' element.
        ! The last two digits represent the kind of matrix:
        ! 11 means the velocity matrix.
        ! 12 means the pressure matrix.
        ! 13 means the concentration matrix.
        ! 14 means the temperature matrix.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! the '|' index of the velocity matrix
        if(ind_kind == 111) then
            if(p_pcol == pncols) then
                global_ind = xBase+yBase+zBase+(local_k-1)*(localncols+1)*localnrows+ &!
                    (local_j-1)*(localncols+1)+local_i
            else
                global_ind = xBase+yBase+zBase+(local_k-1)*localncols*localnrows+ &!
                    (local_j-1)*localncols+local_i
            end if
        ! the '-' index of the velocity matrix
        elseif(ind_kind == 211) then
            if(p_pcol == pncols) then
                uBase = (localncols+1)*localnrows*localnlays
            else
                uBase = localncols*localnrows*localnlays
            end if
            if(p_prow == pnrows) then
                global_ind = xBase+yBase+zBase+uBase+(local_k-1)*localncols*(localnrows+1)+ &!
                    (local_j-1)*localncols+local_i
            else
                global_ind = xBase+yBase+zBase+uBase+(local_k-1)*localncols*localnrows+ &!
                    (local_j-1)*localncols+local_i
            end if
        ! the '|_' index of the velocity matrix
        elseif(ind_kind == 311) then
            if(p_pcol == pncols) then
                uBase = (localncols+1)*localnrows*localnlays
            else
                uBase = localncols*localnrows*localnlays
            end if
            if(p_prow == pnrows) then
                vBase = localncols*(localnrows+1)*localnlays
            else
                vBase = localncols*localnrows*localnlays
            end if
            global_ind = xBase+yBase+zBase+uBase+vBase+(local_k-1)*localncols*localnrows+ &!
                (local_j-1)*localncols + local_i
        ! the '.' index of the pressure, concentration, temperature matrix
        else
            global_ind = pBase + (local_k-1)*localncols*localnrows + (local_j-1)*localncols + local_i
        end if

    end subroutine coordiToGlobalInd

    subroutine coordiToLocalInd(ind_kind, local_i, local_j, local_k, local_ind)

        integer, intent(in) :: ind_kind
        integer, intent(in) :: local_i
        integer, intent(in) :: local_j
        integer, intent(in) :: local_k
        integer, intent(out) :: local_ind

        integer :: indexr, indexu

        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if

        ! x-momentum equation
        if((ind_kind==151).or.(ind_kind==152)) then
            local_ind = (local_k-1)*indexr*localnrows + (local_j-1)*indexr + local_i
        ! y-momentum equation
        elseif((ind_kind==253).or.(ind_kind==254)) then
            local_ind = (local_k-1)*localncols*indexu + (local_j-1)*localncols + local_i
        ! z-momentum equation, continuity equation,
        ! concentration equation, temperature equation
        else
            local_ind = (local_k-1)*localncols*localnrows + (local_j-1)*localncols + local_i
        end if

    end subroutine coordiToLocalInd

    subroutine genExpField(bx, by, bz, local_nx, local_ny, local_nz, field, isField)

        integer, intent(in) :: bx, by, bz, local_nx, local_ny, local_nz
        integer, dimension(:,:,:), pointer, intent(in out) :: field
        logical, intent(out) :: isField

        field(:,:,:) = 0

        if((bx>local_nx).or.(by>local_ny).or.(bz>local_nz)) then
            isField = .false.
        else
            field(bx:local_nx:3, by:local_ny:3, bz:local_nz:3) = 1
            isField = .true.
        end if

    end subroutine genExpField

end module DBF_constructMat


