
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

    subroutine constructAxx(velx, resi, sd_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: sd_kind ! static values or dynamic values

        integer :: m_kind
        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu, indexf, indexb, indextemp
        integer :: i, j, k

        if(sd_kind == 1) then
            m_kind = 1
        else
            m_kind = 2
        end if

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

                    if(velx(i,j,k) == 1) then

                        if((i==0).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid-1, 1, localncols, j, k, fieldInd)
                            call index_convert_local_global(myid, 1, 1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j,k), m_kind, 1, j, k)
                        elseif((pcol/=pncols).and.(i==localncols+1).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and. &!
                            (k/=localnlays+1)) then
                            call index_convert_local_global(myid+1, 1, 1, j, k, fieldInd)
                            call index_convert_local_global(myid, 1, localncols, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j,k), m_kind, localncols, j, k)
                        elseif((j==0).and.(i/=0).and.(i/=indextemp+1).and.(k/=0).and.(k/=localnlays+1).and.(.not.((pcol==1).and. &!
                            (i==1).and.(isDiriX0_p(ylower+j,zlower+k-1)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                            (isDiriX1_p(ylower+j,zlower+k-1)==0)))) then
                            call index_convert_local_global(myid-pncols, 1, i, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 1, i, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1,k), m_kind, i, 1, k)
                        elseif((j==localnrows+1).and.(i/=0).and.(i/=indextemp+1).and.(k/=0).and.(k/=localnlays+1).and.(.not. &!
                            ((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-2,zlower+k-1)==0))).and.(.not.((pcol==pncols).and. &!
                            (i==indextemp).and.(isDiriX1_p(ylower+j-2,zlower+k-1)==0)))) then
                            call index_convert_local_global(myid+pncols, 1, i, 1, k, fieldInd)
                            call index_convert_local_global(myid, 1, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1,k), m_kind, i, localnrows, k)
                        elseif((k==0).and.(i/=0).and.(i/=indextemp+1).and.(j/=0).and.(j/=localnrows+1).and.(.not.((pcol==1).and. &!
                            (i==1).and.(isDiriX0_p(ylower+j-1,zlower+k)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                            (isDiriX1_p(ylower+j-1,zlower+k)==0)))) then
                            call index_convert_local_global(myid-pncols*pnrows, 1, i, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 1, i, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k+1), m_kind, i, j, 1)
                        elseif((k==localnlays+1).and.(i/=0).and.(i/=indextemp+1).and.(j/=0).and.(j/=localnrows+1).and.(.not. &!
                            ((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-2)==0))).and.(.not.((pcol==pncols).and. &!
                            (i==indextemp).and.(isDiriX1_p(ylower+j-1,zlower+k-2)==0)))) then
                            call index_convert_local_global(myid+pncols*pnrows, 1, i, j, 1, fieldInd)
                            call index_convert_local_global(myid, 1, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k-1), m_kind, i, j, localnlays)
                        elseif((i>=1).and.(i<=indextemp).and.(j>=1).and.(j<=localnrows).and.(k>=1).and.(k<=localnlays)) then
                            call index_convert_local_global(myid, 1, i, j, k, fieldInd)
                            if((i>=3).or.((i==2).and.(pcol/=1)).or.((i==2).and.(pcol==1).and. &!
                                (isDiriX0_p(ylower+j-1,zlower+k-1)/=0))) then
                                call index_convert_local_global(myid, 1, i-1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,k), m_kind, i-1, j, k)
                            end if
                            if((i<=indextemp-2).or.((i==indextemp-1).and.(pcol/=pncols)).or.((i==indextemp-1).and. &!
                                (pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)/=0))) then
                                call index_convert_local_global(myid, 1, i+1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,k), m_kind, i+1, j, k)
                            end if
                            if(((j/=1).and.(i>=2).and.(i<=indextemp-1)).or.((j/=1).and.(i==1).and.(pcol/=1)).or. &!
                                ((j/=1).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-2,zlower+k-1)/=0)).or.((j/=1).and. &!
                                (i==indextemp).and.(pcol/=pncols)).or.((j/=1).and.(i==indextemp).and.(pcol==pncols).and. &!
                                (isDiriX1_p(ylower+j-2,zlower+k-1)/=0))) then
                                call index_convert_local_global(myid, 1, i, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,k), m_kind, i, j-1, k)
                            end if
                            if(((j/=localnrows).and.(i>=2).and.(i<=indextemp-1)).or.((j/=localnrows).and.(i==1).and. &!
                                (pcol/=1)).or.((j/=localnrows).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j,zlower+k-1)/=0)) &!
                                .or.((j/=localnrows).and.(i==indextemp).and.(pcol/=pncols)).or.((j/=localnrows).and. &!
                                (i==indextemp).and.(pcol==pncols).and.(isDiriX1_p(ylower+j,zlower+k-1)/=0))) then
                                call index_convert_local_global(myid, 1, i, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,k), m_kind, i, j+1, k)
                            end if
                            if(((k/=1).and.(i>=2).and.(i<=indextemp-1)).or.((k/=1).and.(i==1).and.(pcol/=1)).or. &!
                                ((k/=1).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k-2)/=0)).or.((k/=1).and. &!
                                (i==indextemp).and.(pcol/=pncols)).or.((k/=1).and.(i==indextemp).and.(pcol==pncols).and. &!
                                (isDiriX1_p(ylower+j-1,zlower+k-2)/=0))) then
                                call index_convert_local_global(myid, 1, i, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k-1), m_kind, i, j, k-1)
                            end if
                            if(((k/=localnlays).and.(i>=2).and.(i<=indextemp-1)).or.((k/=localnlays).and.(i==1).and. &!
                                (pcol/=1)).or.((k/=localnlays).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k)/=0)) &!
                                .or.((k/=localnlays).and.(i==indextemp).and.(pcol/=pncols)).or.((k/=localnlays).and. &!
                                (i==indextemp).and.(pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k)/=0))) then
                                call index_convert_local_global(myid, 1, i, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k+1), m_kind, i, j, k+1)
                            end if
                            call setMatValue(fieldInd, fieldInd, resi(i,j,k), m_kind, i, j, k)
                        end if

                    end if

                end do
            end do
        end do

    end subroutine constructAxx

    subroutine constructAxp(pres, resi)

        ! notice that the dimensions of field are different from the dimensions of resi
        ! in the function
        integer, dimension(:,:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: fieldInd, equInd
        integer :: i, j, k

        ! the field index (the pressure)
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

                    if(pres(i,j,k) == 1) then
                        if(i == 0) then
                            call index_convert_local_global(myid-1, 4, localncols, j, k, fieldInd)
                            call index_convert_local_global(myid, 1, 1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j,k), 3, 1, j, k)
                        elseif((i==1).and.(pcol==1)) then
                            call index_convert_local_global(myid, 4, i, j, k, fieldInd)
                            if(isDiriX0_p(ylower+j-1,zlower+k-1) /= 0) then
                                call index_convert_local_global(myid, 1, i, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k), 3, i, j, k)
                            end if
                            call index_convert_local_global(myid, 1, i+1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j,k), 3, i+1, j, k)
                        elseif((i==localncols).and.(pcol/=pncols)) then
                            call index_convert_local_global(myid, 4, localncols, j, k, fieldInd)
                            call index_convert_local_global(myid, 1, localncols, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 3, localncols, j, k)
                        elseif((i==localncols).and.(pcol==pncols)) then
                            call index_convert_local_global(myid, 4, i, j, k, fieldInd)
                            call index_convert_local_global(myid, 1, i, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 3, i, j, k)
                            if(isDiriX1_p(ylower+j-1,zlower+k-1) /= 0) then
                                call index_convert_local_global(myid, 1, i+1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,k), 3, i+1, j, k)
                            end if
                        else
                            call index_convert_local_global(myid, 4, i, j, k, fieldInd)
                            call index_convert_local_global(myid, 1, i, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 3, i, j, k)
                            call index_convert_local_global(myid, 1, i+1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j,k), 3, i+1, j, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine constructAxp

    subroutine constructAyy(vely, resi, sd_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: sd_kind

        integer :: m_kind
        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu, indexf, indexb, indextemp
        integer :: i, j, k

        if(sd_kind == 1) then
            m_kind = 4
        else
            m_kind = 5
        end if

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

                    if(vely(i,j,k) == 1) then

                        if((i==0).and.(j/=0).and.(j/=indextemp+1).and.(k/=0).and.(k/=localnlays+1).and.(.not.((prow==1) &!
                            .and.(j==1).and.(isDiriY0_p(xlower+i,zlower+k-1)==0))).and.(.not.((prow==pnrows).and.(j==indextemp) &!
                            .and.(isDiriY1_p(xlower+i,zlower+k-1)==0)))) then
                            call index_convert_local_global(myid-1, 2, localncols, j, k, fieldInd)
                            call index_convert_local_global(myid, 2, 1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j,k), m_kind, 1, j, k)
                        elseif((i==localncols+1).and.(j/=0).and.(j/=indextemp+1).and.(k/=0).and.(k/=localnlays+1).and. &!
                            (.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-2,zlower+k-1)==0))).and.(.not.((prow==pnrows) &!
                            .and.(j==indextemp).and.(isDiriY1_p(xlower+i-2,zlower+k-1)==0)))) then
                            call index_convert_local_global(myid+1, 2, 1, j, k, fieldInd)
                            call index_convert_local_global(myid, 2, localncols, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j,k), m_kind, localncols, j, k)
                        elseif((j==0).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid-pncols, 2, i, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 2, i, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1,k), m_kind, i, 1, k)
                        elseif((prow/=pnrows).and.(j==localnrows+1).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and. &!
                            (k/=localnlays+1)) then
                            call index_convert_local_global(myid+pncols, 2, i, 1, k, fieldInd)
                            call index_convert_local_global(myid, 2, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1,k), m_kind, i, localnrows, k)
                        elseif((k==0).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=indextemp+1).and.(.not.((prow==1) &!
                            .and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k)==0))).and.(.not.((prow==pnrows).and. &!
                            (j==indextemp).and. (isDiriY1_p(xlower+i-1,zlower+k)==0)))) then
                            call index_convert_local_global(myid-pncols*pnrows, 2, i, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 2, i, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k+1), m_kind, i, j, 1)
                        elseif((k==localnlays+1).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=indextemp+1).and. &!
                            (.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-2)==0))).and.(.not. &!
                            ((prow==pnrows).and.(j==indextemp).and.(isDiriY1_p(xlower+i-1,zlower+k-2)==0)))) then
                            call index_convert_local_global(myid+pncols*pnrows, 2, i, j, 1, fieldInd)
                            call index_convert_local_global(myid, 2, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k-1), m_kind, i, j, localnlays)
                        elseif((j>=1).and.(j<=indextemp).and.(i>=1).and.(i<=localncols).and.(k>=1).and.(k<=localnlays)) then
                            call index_convert_local_global(myid, 2, i, j, k, fieldInd)
                            if(((i/=1).and.(j>=2).and.(j<=indextemp-1)).or.((i/=1).and.(j==1).and.(prow/=1)).or.((i/=1) &!
                                .and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-2,zlower+k-1)/=0)).or.((i/=1).and. &!
                                (j==indextemp).and.(prow/=pnrows)).or.((i/=1).and.(j==indextemp).and.(prow==pnrows).and. &!
                                (isDiriY1_p(xlower+i-2,zlower+k-1)/=0))) then
                                call index_convert_local_global(myid, 2, i-1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,k), m_kind, i-1, j, k)
                            end if
                            if(((i/=localncols).and.(j>=2).and.(j<=indextemp-1)).or.((i/=localncols).and.(j==1).and. &!
                                (prow/=1)).or.((i/=localncols).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i,zlower+k-1)/=0)) &!
                                .or.((i/=localncols).and.(j==indextemp).and.(prow/=pnrows)).or.((i/=localncols).and.(j==indextemp) &!
                                .and.(prow==pnrows).and.(isDiriY1_p(xlower+i,zlower+k-1)/=0))) then
                                call index_convert_local_global(myid, 2, i+1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,k), m_kind, i+1, j, k)
                            end if
                            if((j>=3).or.((j==2).and.(prow/=1)).or.((j==2).and.(prow==1).and. &!
                                (isDiriY0_p(xlower+i-1,zlower+k-1)/=0))) then
                                call index_convert_local_global(myid, 2, i, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,k), m_kind, i, j-1, k)
                            end if
                            if((j<=indextemp-2).or.((j==indextemp-1).and.(prow/=pnrows)).or.((j==indextemp-1).and. &!
                                (prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)/=0))) then
                                call index_convert_local_global(myid, 2, i, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,k), m_kind, i, j+1, k)
                            end if
                            if(((k/=1).and.(j>=2).and.(j<=indextemp-1)).or.((k/=1).and.(j==1).and.(prow/=1)).or. &!
                                ((k/=1).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k-2)/=0)).or.((k/=1).and. &!
                                (j==indextemp).and.(prow/=pnrows)).or.((k/=1).and.(j==indextemp).and.(prow==pnrows).and. &!
                                (isDiriY1_p(xlower+i-1,zlower+k-2)/=0))) then
                                call index_convert_local_global(myid, 2, i, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k-1), m_kind, i, j, k-1)
                            end if
                            if(((k/=localnlays).and.(j>=2).and.(j<=indextemp-1)).or.((k/=localnlays).and.(j==1).and. &!
                                (prow/=1)).or.((k/=localnlays).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k)/=0)) &!
                                .or.((k/=localnlays).and.(j==indextemp).and.(prow/=pnrows)).or.((k/=localnlays).and. &!
                                (j==indextemp).and.(prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k)/=0))) then
                                call index_convert_local_global(myid, 2, i, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k+1), m_kind, i, j, k+1)
                            end if
                            call setMatValue(fieldInd, fieldInd, resi(i,j,k), m_kind, i, j, k)
                        end if

                    end if

                end do
            end do
        end do

    end subroutine constructAyy

    subroutine constructAyp(pres, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: fieldInd, equInd
        integer :: i, j, k

        ! the field index (the pressure)
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

                    if(pres(i,j,k) == 1) then
                        if(j == 0) then
                            call index_convert_local_global(myid-pncols, 4, i, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 2, i, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1,k), 6, i, 1, k)
                        elseif((j==1).and.(prow==1)) then
                            call index_convert_local_global(myid, 4, i, j, k, fieldInd)
                            if(isDiriY0_p(xlower+i-1,zlower+k-1) /= 0) then
                                call index_convert_local_global(myid, 2, i, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k), 6, i, j, k)
                            end if
                            call index_convert_local_global(myid, 2, i, j+1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1,k), 6, i, j+1, k)
                        elseif((j==localnrows).and.(prow/=pnrows)) then
                            call index_convert_local_global(myid, 4, i, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 2, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 6, i, localnrows, k)
                        elseif((j==localnrows).and.(prow==pnrows)) then
                            call index_convert_local_global(myid, 4, i, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 2, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 6, i, localnrows, k)
                            if(isDiriY1_p(xlower+i-1,zlower+k-1) /= 0) then
                                call index_convert_local_global(myid, 2, i, localnrows+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,k), 6, i, localnrows+1, k)
                            end if
                        else
                            call index_convert_local_global(myid, 4, i, j, k, fieldInd)
                            call index_convert_local_global(myid, 2, i, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 6, i, j, k)
                            call index_convert_local_global(myid, 2, i, j+1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1,k), 6, i, j+1, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine constructAyp

    subroutine constructAzz(velz, resi, sd_kind)

        integer, dimension(:,:,:), pointer, intent(in) :: velz
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        integer, intent(in) :: sd_kind

        integer :: m_kind
        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu, indexf, indexb, indextemp
        integer :: i, j, k

        if(sd_kind == 1) then
            m_kind = 7
        else
            m_kind = 8
        end if

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

                    if(velz(i,j,k) == 1) then

                        if((i==0).and.(k/=0).and.(k/=indextemp+1).and.(j/=0).and.(j/=localnrows+1).and.(.not.((play==1).and. &!
                            (k==1).and.(isDiriZ0_p(xlower+i,ylower+j-1)==0))).and.(.not.((play==pnlays).and.(k==indextemp).and. &!
                            (isDiriZ1_p(xlower+i,ylower+j-1)==0)))) then
                            call index_convert_local_global(myid-1, 3, localncols, j, k, fieldInd)
                            call index_convert_local_global(myid, 3, 1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j,k), m_kind, 1, j, k)
                        elseif((i==localncols+1).and.(k/=0).and.(k/=indextemp+1).and.(j/=0).and.(j/=localnrows+1).and. &!
                            (.not.((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-2,ylower+j-1)==0))).and.(.not.((play==pnlays) &!
                            .and.(k==indextemp).and.(isDiriZ1_p(xlower+i-2,ylower+j-1)==0)))) then
                            call index_convert_local_global(myid+1, 3, 1, j, k, fieldInd)
                            call index_convert_local_global(myid, 3, localncols, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j,k), m_kind, localncols, j, k)
                        elseif((j==0).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=indextemp+1).and.(.not.((play==1) &!
                            .and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j)==0))).and.(.not.((play==pnlays).and.(k==indextemp) &!
                            .and.(isDiriZ1_p(xlower+i-1,ylower+j)==0)))) then
                            call index_convert_local_global(myid-pncols, 3, i, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 3, i, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1,k), m_kind, i, 1, k)
                        elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=indextemp+1).and. &!
                            (.not.((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-2)==0))).and.(.not.((play==pnlays) &!
                            .and.(k==indextemp).and.(isDiriZ1_p(xlower+i-1,ylower+j-2)==0)))) then
                            call index_convert_local_global(myid+pncols, 3, i, 1, k, fieldInd)
                            call index_convert_local_global(myid, 3, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1,k), m_kind, i, localnrows, k)
                        elseif((k==0).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call index_convert_local_global(myid-pncols*pnrows, 3, i, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 3, i, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k+1), m_kind, i, j, 1)
                        elseif((play/=pnlays).and.(k==localnlays+1).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and. &!
                            (j/=localnrows+1)) then
                            call index_convert_local_global(myid+pncols*pnrows, 3, i, j, 1, fieldInd)
                            call index_convert_local_global(myid, 3, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k-1), m_kind, i, j, localnlays)
                        elseif((k>=1).and.(k<=indextemp).and.(i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows)) then
                            call index_convert_local_global(myid, 3, i, j, k, fieldInd)
                            if(((i/=1).and.(k>=2).and.(k<=indextemp-1)).or.((i/=1).and.(k==1).and.(play/=1)).or.((i/=1).and. &!
                                (k==1).and.(play==1).and.(isDiriZ0_p(xlower+i-2,ylower+j-1)/=0)).or.((i/=1).and.(k==indextemp) &!
                                .and.(play/=pnlays)).or.((i/=1).and.(k==indextemp).and.(play==pnlays).and. &!
                                (isDiriZ1_p(xlower+i-2,ylower+j-1)/=0))) then
                                call index_convert_local_global(myid, 3, i-1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,k), m_kind, i-1, j, k)
                            end if
                            if(((i/=localncols).and.(k>=2).and.(k<=indextemp-1)).or.((i/=localncols).and.(k==1).and.(play/=1)) &!
                                .or.((i/=localncols).and.(k==1).and.(play==1).and.(isDiriZ0_p(xlower+i,ylower+j-1)/=0)).or. &!
                                ((i/=localncols).and.(k==indextemp).and.(play/=pnlays)).or.((i/=localncols).and.(k==indextemp) &!
                                .and.(play==pnlays).and.(isDiriZ1_p(xlower+i,ylower+j-1)/=0))) then
                                call index_convert_local_global(myid, 3, i+1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,k), m_kind, i+1, j, k)
                            end if
                            if(((j/=1).and.(k>=2).and.(k<=indextemp-1)).or.((j/=1).and.(k==1).and.(play/=1)).or. &!
                                ((j/=1).and.(k==1).and.(play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-2)/=0)).or.((j/=1).and. &!
                                (k==indextemp).and.(play/=pnlays)).or.((j/=1).and.(k==indextemp).and.(play==pnlays).and. &!
                                (isDiriZ1_p(xlower+i-1,ylower+j-2)/=0))) then
                                call index_convert_local_global(myid, 3, i, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,k), m_kind, i, j-1, k)
                            end if
                            if(((j/=localnrows).and.(k>=2).and.(k<=indextemp-1)).or.((j/=localnrows).and.(k==1).and. &!
                                (play/=1)).or.((j/=localnrows).and.(k==1).and.(play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j)/=0)) &!
                                .or.((j/=localnrows).and.(k==indextemp).and.(play/=pnlays)).or.((j/=localnrows).and. &!
                                (k==indextemp).and.(play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j)/=0))) then
                                call index_convert_local_global(myid, 3, i, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,k), m_kind, i, j+1, k)
                            end if
                            if((k>=3).or.((k==2).and.(play/=1)).or.((k==2).and.(play==1).and. &!
                                (isDiriZ0_p(xlower+i-1,ylower+j-1)/=0))) then
                                call index_convert_local_global(myid, 3, i, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k-1), m_kind, i, j, k-1)
                            end if
                            if((k<=indextemp-2).or.((k==indextemp-1).and.(play/=pnlays)).or.((k==indextemp-1).and. &!
                                (play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)/=0))) then
                                call index_convert_local_global(myid, 3, i, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k+1), m_kind, i, j, k+1)
                            end if
                            call setMatValue(fieldInd, fieldInd, resi(i,j,k), m_kind, i, j, k)
                        end if

                    end if

                end do
            end do
        end do

    end subroutine constructAzz

    subroutine constructAzp(pres, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: fieldInd, equInd
        integer :: i, j, k

        ! the field index (the pressure)
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

                    if(pres(i,j,k) == 1) then
                        if(k == 0) then
                            call index_convert_local_global(myid-pncols*pnrows, 4, i, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 3, i, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k+1), 9, i, j, 1)
                        elseif((k==1).and.(play==1)) then
                            call index_convert_local_global(myid, 4, i, j, k, fieldInd)
                            if(isDiriZ0_p(xlower+i-1,ylower+j-1) /= 0) then
                                call index_convert_local_global(myid, 3, i, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k), 9, i, j, k)
                            end if
                            call index_convert_local_global(myid, 3, i, j, k+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k+1), 9, i, j, k+1)
                        elseif((k==localnlays).and.(play/=pnlays)) then
                            call index_convert_local_global(myid, 4, i, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 3, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 9, i, j, localnlays)
                        elseif((k==localnlays).and.(play==pnlays)) then
                            call index_convert_local_global(myid, 4, i, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 3, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 9, i, j, localnlays)
                            if(isDiriZ1_p(xlower+i-1,ylower+j-1) /= 0) then
                                call index_convert_local_global(myid, 3, i, j, localnlays+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k+1), 9, i, j, localnlays+1)
                            end if
                        else
                            call index_convert_local_global(myid, 4, i, j, k, fieldInd)
                            call index_convert_local_global(myid, 3, i, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 9, i, j, k)
                            call index_convert_local_global(myid, 3, i, j, k+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k+1), 9, i, j, k+1)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine constructAzp

    subroutine constructAcx(velx, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: fieldInd, equInd
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

                    if(velx(i,j,k) == 1) then
                        if(i == 1) then
                            call index_convert_local_global(myid, 1, 1, j, k, fieldInd)
                            call index_convert_local_global(myid, 4, 1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 10, 1, j, k)
                        elseif((i==localncols+1).and.(pcol/=pncols)) then
                            call index_convert_local_global(myid+1, 1, 1, j, k, fieldInd)
                            call index_convert_local_global(myid, 4, localncols, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j,k), 10, localncols, j, k)
                        elseif((i==localncols+1).and.(pcol==pncols)) then
                            call index_convert_local_global(myid, 1, localncols+1, j, k, fieldInd)
                            call index_convert_local_global(myid, 4, localncols, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j,k), 10, localncols, j, k)
                        else
                            call index_convert_local_global(myid, 1, i, j, k, fieldInd)
                            call index_convert_local_global(myid, 4, i-1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j,k), 10, i-1, j, k)
                            call index_convert_local_global(myid, 4, i, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 10, i, j, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine constructAcx

    subroutine constructAcy(vely, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: fieldInd, equInd
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

                    if(vely(i,j,k) == 1) then
                        if(j == 1) then
                            call index_convert_local_global(myid, 2, i, 1, k, fieldInd)
                            call index_convert_local_global(myid, 4, i, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 11, i, 1, k)
                        elseif((j==localnrows+1).and.(prow/=pnrows)) then
                            call index_convert_local_global(myid+pncols, 2, i, 1, k, fieldInd)
                            call index_convert_local_global(myid, 4, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1,k), 11, i, localnrows, k)
                        elseif((j==localnrows+1).and.(prow==pnrows)) then
                            call index_convert_local_global(myid, 2, i, localnrows+1, k, fieldInd)
                            call index_convert_local_global(myid, 4, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1,k), 11, i, localnrows, k)
                        else
                            call index_convert_local_global(myid, 2, i, j, k, fieldInd)
                            call index_convert_local_global(myid, 4, i, j-1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1,k), 11, i, j-1, k)
                            call index_convert_local_global(myid, 4, i, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 11, i, j, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine constructAcy

    subroutine constructAcz(velz, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: velz
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: fieldInd, equInd
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

                    if(velz(i,j,k) == 1) then
                        if(k == 1) then
                            call index_convert_local_global(myid, 3, i, j, 1, fieldInd)
                            call index_convert_local_global(myid, 4, i, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 12, i, j, 1)
                        elseif((k==localnlays+1).and.(play/=pnlays)) then
                            call index_convert_local_global(myid+pncols*pnrows, 3, i, j, 1, fieldInd)
                            call index_convert_local_global(myid, 4, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k-1), 12, i, j, localnlays)
                        elseif((k==localnlays+1).and.(play==pnlays)) then
                            call index_convert_local_global(myid, 3, i, j, localnlays+1, fieldInd)
                            call index_convert_local_global(myid, 4, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k-1), 12, i, j, localnlays)
                        else
                            call index_convert_local_global(myid, 3, i, j, k, fieldInd)
                            call index_convert_local_global(myid, 4, i, j, k-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k-1), 12, i, j, k-1)
                            call index_convert_local_global(myid, 4, i, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k), 12, i, j, k)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine constructAcz

    subroutine constructAcp(pres, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        
        integer :: fieldInd
        integer :: i, j, k

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    if(pres(i,j,k) == 1) then
                        call index_convert_local_global(myid, 4, i, j, k, fieldInd)
                        call setMatValue(fieldInd, fieldInd, resi(i,j,k), 13, i, j, k)
                    end if

                end do
            end do
        end do

    end subroutine constructAcp

    subroutine constructAcf(conc, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: conc
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi
        
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: fieldInd, equInd
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

                    if(conc(i,j,k) == 1) then

                        if((i==0).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid-1, 5, localncols, j, k, fieldInd)
                            call index_convert_local_global(myid, 5, 1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(1,j,k), 14, 1, j, k)
                            if(j/=1) then
                                call index_convert_local_global(myid, 5, 1, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(1,j-1,k), 14, 1, j-1, k)
                            end if
                            if(j/=localnrows) then
                                call index_convert_local_global(myid, 5, 1, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(1,j+1,k), 14, 1, j+1, k)
                            end if
                            if(k/=1) then
                                call index_convert_local_global(myid, 5, 1, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(1,j,k-1), 14, 1, j, k-1)
                            end if
                            if(k/=localnlays) then
                                call index_convert_local_global(myid, 5, 1, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(1,j,k+1), 14, 1, j, k+1)
                            end if
                        end if

                        if((i==localncols+1).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid+1, 5, 1, j, k, fieldInd)
                            call index_convert_local_global(myid, 5, localncols, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(localncols,j,k), 14, localncols, j, k)
                            if(j/=1) then
                                call index_convert_local_global(myid, 5, localncols, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(localncols,j-1,k), 14, localncols, j-1, k)
                            end if
                            if(j/=localnrows) then
                                call index_convert_local_global(myid, 5, localncols, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(localncols,j+1,k), 14, localncols, j+1, k)
                            end if
                            if(k/=1) then
                                call index_convert_local_global(myid, 5, localncols, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(localncols,j,k-1), 14, localncols, j, k-1)
                            end if
                            if(k/=localnlays) then
                                call index_convert_local_global(myid, 5, localncols, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(localncols,j,k+1), 14, localncols, j, k+1)
                            end if
                        end if

                        if((j==0).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid-pncols, 5, i, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 5, i, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,1,k), 14, i, 1, k)
                            if(i/=1) then
                                call index_convert_local_global(myid, 5, i-1, 1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,1,k), 14, i-1, 1, k)
                            end if
                            if(i/=localncols) then
                                call index_convert_local_global(myid, 5, i+1, 1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,1,k), 14, i+1, 1, k)
                            end if
                            if(k/=1) then
                                call index_convert_local_global(myid, 5, i, 1, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,1,k-1), 14, i, 1, k-1)
                            end if
                            if(k/=localnlays) then
                                call index_convert_local_global(myid, 5, i, 1, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,1,k+1), 14, i, 1, k+1)
                            end if
                        end if

                        if((j==localnrows+1).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid+pncols, 5, i, 1, k, fieldInd)
                            call index_convert_local_global(myid, 5, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,localnrows,k), 14, i, localnrows, k)
                            if(i/=1) then
                                call index_convert_local_global(myid, 5, i-1, localnrows, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,localnrows,k), 14, i-1, localnrows, k)
                            end if
                            if(i/=localncols) then
                                call index_convert_local_global(myid, 5, i+1, localnrows, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,localnrows,k), 14, i+1, localnrows, k)
                            end if
                            if(k/=1) then
                                call index_convert_local_global(myid, 5, i, localnrows, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,localnrows,k-1), 14, i, localnrows, k-1)
                            end if
                            if(k/=localnlays) then
                                call index_convert_local_global(myid, 5, i, localnrows, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,localnrows,k+1), 14, i, localnrows, k+1)
                            end if
                        end if

                        if((k==0).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call index_convert_local_global(myid-pncols*pnrows, 5, i, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 5, i, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,1), 14, i, j, 1)
                            if(i/=1) then
                                call index_convert_local_global(myid, 5, i-1, j, 1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,1), 14, i-1, j, 1)
                            end if
                            if(i/=localncols) then
                                call index_convert_local_global(myid, 5, i+1, j, 1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,1), 14, i+1, j, 1)
                            end if
                            if(j/=1) then
                                call index_convert_local_global(myid, 5, i, j-1, 1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,1), 14, i, j-1, 1)
                            end if
                            if(j/=localnrows) then
                                call index_convert_local_global(myid, 5, i, j+1, 1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,1), 14, i, j+1, 1)
                            end if
                        end if

                        if((k==localnlays+1).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call index_convert_local_global(myid+pncols*pnrows, 5, i, j, 1, fieldInd)
                            call index_convert_local_global(myid, 5, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,localnlays), 14, i, j, localnlays)
                            if(i/=1) then
                                call index_convert_local_global(myid, 5, i-1, j, localnlays, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,localnlays), 14, i-1, j, localnlays)
                            end if
                            if(i/=localncols) then
                                call index_convert_local_global(myid, 5, i+1, j, localnlays, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,localnlays), 14, i+1, j, localnlays)
                            end if
                            if(j/=1) then
                                call index_convert_local_global(myid, 5, i, j-1, localnlays, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,localnlays), 14, i, j-1, localnlays)
                            end if
                            if(j/=localnrows) then
                                call index_convert_local_global(myid, 5, i, j+1, localnlays, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,localnlays), 14, i, j+1, localnlays)
                            end if
                        end if

                        if((i==0).and.(j==0).and.(k>=1).and.(k<=localnlays)) then
                            call index_convert_local_global(myid-1-pncols, 5, localncols, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 5, 1, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(1,1,k), 14, 1, 1, k)
                        end if

                        if((i==0).and.(j==localnrows+1).and.(k>=1).and.(k<=localnlays)) then
                            call index_convert_local_global(myid-1+pncols, 5, localncols, 1, k, fieldInd)
                            call index_convert_local_global(myid, 5, 1, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(1,localnrows,k), 14, 1, localnrows, k)
                        end if

                        if((i==0).and.(k==0).and.(j>=1).and.(j<=localnrows)) then
                            call index_convert_local_global(myid-1-pncols*pnrows, 5, localncols, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 5, 1, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(1,j,1), 14, 1, j, 1)
                        end if

                        if((i==0).and.(k==localnlays+1).and.(j>=1).and.(j<=localnrows)) then
                            call index_convert_local_global(myid-1+pncols*pnrows, 5, localncols, j, 1, fieldInd)
                            call index_convert_local_global(myid, 5, 1, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(1,j,localnlays), 14, 1, j, localnlays)
                        end if

                        if((i==localncols+1).and.(j==0).and.(k>=1).and.(k<=localnlays)) then
                            call index_convert_local_global(myid+1-pncols, 5, 1, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 5, localncols, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(localncols,1,k), 14, localncols, 1, k)
                        end if

                        if((i==localncols+1).and.(j==localnrows+1).and.(k>=1).and.(k<=localnlays)) then
                            call index_convert_local_global(myid+1+pncols, 5, 1, 1, k, fieldInd)
                            call index_convert_local_global(myid, 5, localncols, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(localncols,localnrows,k), 14, localncols, localnrows, k)
                        end if

                        if((i==localncols+1).and.(k==0).and.(j>=1).and.(j<=localnrows)) then
                            call index_convert_local_global(myid+1-pncols*pnrows, 5, 1, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 5, localncols, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(localncols,j,1), 14, localncols, j, 1)
                        end if

                        if((i==localncols+1).and.(k==localnlays+1).and.(j>=1).and.(j<=localnrows)) then
                            call index_convert_local_global(myid+1+pncols*pnrows, 5, 1, j, 1, fieldInd)
                            call index_convert_local_global(myid, 5, localncols, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(localncols,j,localnlays), 14, localncols, j, localnlays)
                        end if

                        if((j==0).and.(k==0).and.(i>=1).and.(i<=localncols)) then
                            call index_convert_local_global(myid-pncols-pncols*pnrows, 5, i, localnrows, &!
                                localnlays, fieldInd)
                            call index_convert_local_global(myid, 5, i, 1, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,1,1), 14, i, 1, 1)
                        end if

                        if((j==0).and.(k==localnlays+1).and.(i>=1).and.(i<=localncols)) then
                            call index_convert_local_global(myid-pncols+pncols*pnrows, 5, i, localnrows, 1, fieldInd)
                            call index_convert_local_global(myid, 5, i, 1, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,1,localnlays), 14, i, 1, localnlays)
                        end if

                        if((j==localnrows+1).and.(k==0).and.(i>=1).and.(i<=localncols)) then
                            call index_convert_local_global(myid+pncols-pncols*pnrows, 5, i, 1, localnlays, fieldInd)
                            call index_convert_local_global(myid, 5, i, localnrows, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,localnrows,1), 14, i, localnrows, 1)
                        end if

                        if((j==localnrows+1).and.(k==localnlays+1).and.(i>=1).and.(i<=localncols)) then
                            call index_convert_local_global(myid+pncols+pncols*pnrows, 5, i, 1, 1, fieldInd)
                            call index_convert_local_global(myid, 5, i, localnrows, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,localnrows,localnlays), 14, i, localnrows, localnlays)
                        end if

                        if((i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows).and.(k>=1).and.(k<=localnlays)) then
                            call index_convert_local_global(myid, 5, i, j, k, fieldInd)
                            call setMatValue(fieldInd, fieldInd, resi(i,j,k), 14, i, j, k)
                            if(i/=1) then
                                call index_convert_local_global(myid, 5, i-1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,k), 14, i-1, j, k)
                            end if
                            if(i/=localncols) then
                                call index_convert_local_global(myid, 5, i+1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,k), 14, i+1, j, k)
                            end if
                            if(j/=1) then
                                call index_convert_local_global(myid, 5, i, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,k), 14, i, j-1, k)
                            end if
                            if(j/=localnrows) then
                                call index_convert_local_global(myid, 5, i, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,k), 14, i, j+1, k)
                            end if
                            if(k/=1) then
                                call index_convert_local_global(myid, 5, i, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k-1), 14, i, j, k-1)
                            end if
                            if(k/=localnlays) then
                                call index_convert_local_global(myid, 5, i, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k+1), 14, i, j, k+1)
                            end if
                            if((i/=1).and.(j/=1)) then
                                call index_convert_local_global(myid, 5, i-1, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j-1,k), 14, i-1, j-1, k)
                            end if
                            if((i/=1).and.(j/=localnrows)) then
                                call index_convert_local_global(myid, 5, i-1, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j+1,k), 14, i-1, j+1, k)
                            end if
                            if((i/=1).and.(k/=1)) then
                                call index_convert_local_global(myid, 5, i-1, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,k-1), 14, i-1, j, k-1)
                            end if
                            if((i/=1).and.(k/=localnlays)) then
                                call index_convert_local_global(myid, 5, i-1, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,k+1), 14, i-1, j, k+1)
                            end if
                            if((i/=localncols).and.(j/=1)) then
                                call index_convert_local_global(myid, 5, i+1, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j-1,k), 14, i+1, j-1, k)
                            end if
                            if((i/=localncols).and.(j/=localnrows)) then
                                call index_convert_local_global(myid, 5, i+1, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j+1,k), 14, i+1, j+1, k)
                            end if
                            if((i/=localncols).and.(k/=1)) then
                                call index_convert_local_global(myid, 5, i+1, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,k-1), 14, i+1, j, k-1)
                            end if
                            if((i/=localncols).and.(k/=localnlays)) then
                                call index_convert_local_global(myid, 5, i+1, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,k+1), 14, i+1, j, k+1)
                            end if
                            if((j/=1).and.(k/=1)) then
                                call index_convert_local_global(myid, 5, i, j-1, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,k-1), 14, i, j-1, k-1)
                            end if
                            if((j/=1).and.(k/=localnlays)) then
                                call index_convert_local_global(myid, 5, i, j-1, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,k+1), 14, i, j-1, k+1)
                            end if
                            if((j/=localnrows).and.(k/=1)) then
                                call index_convert_local_global(myid, 5, i, j+1, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,k-1), 14, i, j+1, k-1)
                            end if
                            if((j/=localnrows).and.(k/=localnlays)) then
                                call index_convert_local_global(myid, 5, i, j+1, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,k+1), 14, i, j+1, k+1)
                            end if
                        end if

                    end if

                end do
            end do
        end do

    end subroutine constructAcf

    subroutine constructAtem(tempe, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: tempe
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: resi

        integer :: fieldInd, equInd
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

                    if(tempe(i,j,k) == 1) then

                        if((i==0).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid-1, 6, localncols, j, k, fieldInd)
                            call index_convert_local_global(myid, 6, 1, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j,k), 15, 1, j, k)
                        elseif((i==localncols+1).and.(j/=0).and.(j/=localnrows+1).and.(k/=0).and. &!
                            (k/=localnlays+1)) then
                            call index_convert_local_global(myid+1, 6, 1, j, k, fieldInd)
                            call index_convert_local_global(myid, 6, localncols, j, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j,k), 15, localncols, j, k)
                        elseif((j==0).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid-pncols, 6, i, localnrows, k, fieldInd)
                            call index_convert_local_global(myid, 6, i, 1, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1,k), 15, i, 1, k)
                        elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1).and.(k/=0).and.(k/=localnlays+1)) then
                            call index_convert_local_global(myid+pncols, 6, i, 1, k, fieldInd)
                            call index_convert_local_global(myid, 6, i, localnrows, k, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1,k), 15, i, localnrows, k)
                        elseif((k==0).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call index_convert_local_global(myid-pncols*pnrows, 6, i, j, localnlays, fieldInd)
                            call index_convert_local_global(myid, 6, i, j, 1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k+1), 15, i, j, 1)
                        elseif((k==localnlays+1).and.(i/=0).and.(i/=localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                            call index_convert_local_global(myid+pncols*pnrows, 6, i, j, 1, fieldInd)
                            call index_convert_local_global(myid, 6, i, j, localnlays, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j,k-1), 15, i, j, localnlays)
                        elseif((i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows).and.(k>=1).and.(k<=localnlays)) then
                            call index_convert_local_global(myid, 6, i, j, k, fieldInd)
                            call setMatValue(fieldInd, fieldInd, resi(i,j,k), 15, i, j, k)
                            if(i /= 1) then
                                call index_convert_local_global(myid, 6, i-1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i-1,j,k), 15, i-1, j, k)
                            end if
                            if(i /= localncols) then
                                call index_convert_local_global(myid, 6, i+1, j, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i+1,j,k), 15, i+1, j, k)
                            end if
                            if(j /= 1) then
                                call index_convert_local_global(myid, 6, i, j-1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j-1,k), 15, i, j-1, k)
                            end if
                            if(j /= localnrows) then
                                call index_convert_local_global(myid, 6, i, j+1, k, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j+1,k), 15, i, j+1, k)
                            end if
                            if(k /= 1) then
                                call index_convert_local_global(myid, 6, i, j, k-1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k-1), 15, i, j, k-1)
                            end if
                            if(k /= localnlays) then
                                call index_convert_local_global(myid, 6, i, j, k+1, equInd)
                                call setMatValue(fieldInd, equInd, resi(i,j,k+1), 15, i, j, k+1)
                            end if
                        end if

                    end if

                end do
            end do
        end do

    end subroutine constructAtem

    subroutine setMatValue(col, row, value, m_kind, eq_i, eq_j, eq_k)

        integer, intent(in) :: col
        integer, intent(in) :: row
        real(kind=8), intent(in) :: value
        integer, intent(in) :: m_kind ! matrix kind
        integer, intent(in) :: eq_i ! equation x-direction coordinate
        integer, intent(in) :: eq_j ! equation y-direction coordinate
        integer, intent(in) :: eq_k ! equation z-direction coordinate

        integer, dimension(:), pointer :: Acols
        integer, dimension(:), pointer :: Arows
        real(kind=8), dimension(:), pointer :: Avalues
        integer, dimension(:), pointer :: AEntryBase
        integer, dimension(:), pointer :: AEntryNum

        integer :: eq_kind ! equation kind
        integer :: pos, base, tail, shend, m ,n
        integer :: left, right, mid

        ! matrix kind
        if(m_kind == 1) then ! AxxStatic
            Acols => AxxCols
            Arows => AxxRows
            Avalues => AxxStaticValues
            AEntryBase => AxxEntryBase
            AEntryNum => AxxEntryNum
            eq_kind = 1
        elseif(m_kind == 2) then ! AxxDyn
            Acols => AxxCols
            Arows => AxxRows
            Avalues => AxxDynValues
            AEntryBase => AxxEntryBase
            AEntryNum => AxxEntryNum
            eq_kind = 1
        elseif(m_kind == 3) then ! Axp
            Acols => AxpCols
            Arows => AxpRows
            Avalues => AxpValues
            AEntryBase => AxpEntryBase
            AEntryNum => AxpEntryNum
            eq_kind = 1
        elseif(m_kind == 4) then ! AyyStatic
            Acols => AyyCols
            Arows => AyyRows
            Avalues => AyyStaticValues
            AEntryBase => AyyEntryBase
            AEntryNum => AyyEntryNum
            eq_kind = 2
        elseif(m_kind == 5) then ! AyyDyn
            Acols => AyyCols
            Arows => AyyRows
            Avalues => AyyDynValues
            AEntryBase => AyyEntryBase
            AEntryNum => AyyEntryNum
            eq_kind = 2
        elseif(m_kind == 6) then ! Ayp
            Acols => AypCols
            Arows => AypRows
            Avalues => AypValues
            AEntryBase => AypEntryBase
            AEntryNum => AypEntryNum
            eq_kind = 2
        elseif(m_kind == 7) then ! AzzStatic
            Acols => AzzCols
            Arows => AzzRows
            Avalues => AzzStaticValues
            AEntryBase => AzzEntryBase
            AEntryNum => AzzEntryNum
            eq_kind = 3
        elseif(m_kind == 8) then ! AzzDyn
            Acols => AzzCols
            Arows => AzzRows
            Avalues => AzzDynValues
            AEntryBase => AzzEntryBase
            AEntryNum => AzzEntryNum
            eq_kind = 3
        elseif(m_kind == 9) then ! Azp
            Acols => AzpCols
            Arows => AzpRows
            Avalues => AzpValues
            AEntryBase => AzpEntryBase
            AEntryNum => AzpEntryNum
            eq_kind = 3
        elseif(m_kind == 10) then ! Acx
            Acols => AcxCols
            Arows => AcxRows
            Avalues => AcxValues
            AEntryBase => AcxEntryBase
            AEntryNum => AcxEntryNum
            eq_kind = 4
        elseif(m_kind == 11) then ! Acy
            Acols => AcyCols
            Arows => AcyRows
            Avalues => AcyValues
            AEntryBase => AcyEntryBase
            AEntryNum => AcyEntryNum
            eq_kind = 4
        elseif(m_kind == 12) then ! Acz
            Acols => AczCols
            Arows => AczRows
            Avalues => AczValues
            AEntryBase => AczEntryBase
            AEntryNum => AczEntryNum
            eq_kind = 4
        elseif(m_kind == 13) then ! Acp
            Acols => AcpCols
            Arows => AcpRows
            Avalues => AcpValues
            AEntryBase => AcpEntryBase
            AEntryNum => AcpEntryNum
            eq_kind = 4
        elseif(m_kind == 14) then ! Acf
            Acols => AcfCols
            Arows => AcfRows
            Avalues => AcfValues
            AEntryBase => AcfEntryBase
            AEntryNum => AcfEntryNum
            eq_kind = 5
        elseif(m_kind == 15) then ! Atem
            Acols => AtemCols
            Arows => AtemRows
            Avalues => AtemValues
            AEntryBase => AtemEntryBase
            AEntryNum => AtemEntryNum
            eq_kind = 6
        end if

        call equCoorditoInd(eq_kind, eq_i, eq_j, eq_k, pos)

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

    ! change the index of the unknowns from the local index to the global index
    subroutine index_convert_local_global(pid, kind, local_i, local_j, local_k, global_ind)

        integer, intent(in) :: pid
        integer, intent(in) :: kind
        integer, intent(in) :: local_i, local_j, local_k
        integer, intent(out) :: global_ind

        integer :: p_pcol, p_prow, p_play
        integer :: base, ubase, vbase, wbase, Cfbase, Tembase

        p_play = pid/(pnrows*pncols)+1
        p_prow = (pid-(p_play-1)*pnrows*pncols)/pncols+1
        p_pcol = (pid-(p_play-1)*pnrows*pncols)-(p_prow-1)*pncols+1

        base = (p_play-1)*localnlays*((nx+1)*ny+nx*(ny+1)+2*nx*ny)
        base = base + (p_prow-1)*localnrows*localnlays*((nx+1)+3*nx)
        if(p_play == pnlays) then
            base = base + nx*(p_prow-1)*localnrows
        end if
        base = base + (p_pcol-1)*localncols*(4*localnrows*localnlays)
        if(p_prow == pnrows) then
            base = base + (p_pcol-1)*localncols*localnlays
        end if
        if(p_play == pnlays) then
            base = base + (p_pcol-1)*localncols*localnrows
        end if

        if(p_pcol == pncols) then
            ubase = (localncols+1)*localnrows*localnlays
        else
            ubase = localncols*localnrows*localnlays
        end if
        if(p_prow == pnrows) then
            vbase = localncols*(localnrows+1)*localnlays
        else
            vbase = localncols*localnrows*localnlays
        end if
        if(p_play == pnlays) then
            wbase = localncols*localnrows*(localnlays+1)
        else
            wbase = localncols*localnrows*localnlays
        end if

        Cfbase = (p_play-1)*nx*ny*localnlays + (p_prow-1)*nx*localnrows*localnlays + &!
            (p_pcol-1)*localncols*localnrows*localnlays

        Tembase = (p_play-1)*nx*ny*localnlays + (p_prow-1)*nx*localnrows*localnlays + &!
            (p_pcol-1)*localncols*localnrows*localnlays

        ! the 'x' index
        if(kind == 1) then
            if(p_pcol == pncols) then
                global_ind = base + (local_k-1)*(localncols+1)*localnrows + (local_j-1)*(localncols+1) + local_i
            else
                global_ind = base + (local_k-1)*localncols*localnrows + (local_j-1)*localncols + local_i
            end if
        ! the 'y' index
        elseif(kind == 2) then
            if(p_prow == pnrows) then
                global_ind = base + ubase + (local_k-1)*localncols*(localnrows+1) + (local_j-1)*localncols + local_i
            else
                global_ind = base + ubase + (local_k-1)*localncols*localnrows + (local_j-1)*localncols + local_i
            end if
        ! the 'z' index
        elseif(kind == 3) then
            global_ind = base + ubase + vbase + (local_k-1)*localncols*localnrows + (local_j-1)*localncols + local_i
        ! the 'p' index
        elseif(kind == 4) then
            global_ind = base + ubase + vbase + wbase + (local_k-1)*localncols*localnrows + &!
                (local_j-1)*localncols+local_i
        ! the 'Cf' index
        elseif(kind == 5) then
            global_ind = Cfbase + (local_k-1)*localncols*localnrows + (local_j-1)*localncols + local_i
        ! the 'Tem' index
        elseif(kind == 6) then
            global_ind = Tembase + (local_k-1)*localncols*localnrows + (local_j-1)*localncols + local_i
        end if

    end subroutine index_convert_local_global

    subroutine equCoorditoInd(eq_kind, eq_i, eq_j, eq_k, eq_ind)

        integer, intent(in) :: eq_kind
        integer, intent(in) :: eq_i
        integer, intent(in) :: eq_j
        integer, intent(in) :: eq_k
        integer, intent(out) :: eq_ind

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
        if(eq_kind == 1) then
            eq_ind = (eq_k-1)*indexr*localnrows + (eq_j-1)*indexr + eq_i
        ! y-momentum equation
        elseif(eq_kind == 2) then
            eq_ind = (eq_k-1)*localncols*indexu + (eq_j-1)*localncols + eq_i
        ! z-momentum equation, continuity equation,
        ! concentration equation, temperature equation
        else
            eq_ind = (eq_k-1)*localncols*localnrows + (eq_j-1)*localncols + eq_i
        end if

    end subroutine equCoorditoInd

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


