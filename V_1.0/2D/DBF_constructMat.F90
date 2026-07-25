
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

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: sd_kind ! static values or dynamic values

        integer :: m_kind, fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu, indextemp
        integer :: i, j

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
        indexr = localncols+1
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

                if(velx(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid-1, 1, localncols, j, fieldInd)
                        call index_convert_local_global(myid, 1, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, 1, j)
                    elseif((pcol/=pncols).and.(i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid+1, 1, 1, j, fieldInd)
                        call index_convert_local_global(myid, 1, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, localncols, j)
                    elseif((j==0).and.(i/=0).and.(i/=indextemp+1).and.(.not.((pcol==1).and.(i==1).and. &!
                        (isDiriX0_p(ylower+j)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                        (isDiriX1_p(ylower+j)==0)))) then
                        call index_convert_local_global(myid-pncols, 1, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 1, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, 1)
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=indextemp+1).and.(.not.((pcol==1).and. &!
                        (i==1).and.(isDiriX0_p(ylower+j-2)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                        (isDiriX1_p(ylower+j-2)==0)))) then
                        call index_convert_local_global(myid+pncols, 1, i, 1, fieldInd)
                        call index_convert_local_global(myid, 1, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, localnrows)
                    elseif((i>=1).and.(i<=indextemp).and.(j>=1).and.(j<=localnrows)) then
                        call index_convert_local_global(myid, 1, i, j, fieldInd)
                        if(((j/=1).and.(i>=2).and.(i<=indextemp-1)).or.((j/=1).and.(i==1).and.(pcol/=1)).or. &!
                            ((j/=1).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-2)/=0)).or.((j/=1).and.(i==indextemp) &!
                            .and.(pcol/=pncols)).or.((j/=1).and.(i==indextemp).and.(pcol==pncols).and. &!
                            (isDiriX1_p(ylower+j-2)/=0))) then
                            call index_convert_local_global(myid, 1, i, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, j-1)
                        end if
                        if((i>=3).or.((i==2).and.(pcol/=1)).or.((i==2).and.(pcol==1).and.(isDiriX0_p(ylower+j-1)/=0))) then
                            call index_convert_local_global(myid, 1, i-1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, i-1, j)
                        end if
                        call setMatValue(fieldInd, fieldInd, resi(i,j), m_kind, i, j)
                        if((i<=indextemp-2).or.((i==indextemp-1).and.(pcol/=pncols)).or.((i==indextemp-1).and. &!
                            (pcol==pncols).and.(isDiriX1_p(ylower+j-1)/=0))) then
                            call index_convert_local_global(myid, 1, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, i+1, j)
                        end if
                        if(((j/=localnrows).and.(i>=2).and.(i<=indextemp-1)).or.((j/=localnrows).and.(i==1).and. &!
                            (pcol/=1)).or.((j/=localnrows).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j)/=0)).or. &!
                            ((j/=localnrows).and.(i==indextemp).and.(pcol/=pncols)).or.((j/=localnrows).and. &!
                            (i==indextemp).and.(pcol==pncols).and.(isDiriX1_p(ylower+j)/=0))) then
                            call index_convert_local_global(myid, 1, i, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine constructAxx

    subroutine constructAxp(pres, resi)

        ! notice that the dimensions of field are different from the dimensions of resi
        ! in the function
        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi

        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
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

                if(pres(i,j) == 1) then
                    if(i == 0) then
                        call index_convert_local_global(myid-1, 3, localncols, j, fieldInd)
                        call index_convert_local_global(myid, 1, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), 3, 1, j)
                    elseif((i==1).and.(pcol==1)) then
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        if(isDiriX0_p(ylower+j-1) /= 0) then
                            call index_convert_local_global(myid, 1, i, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j), 3, i, j)
                        end if
                        call index_convert_local_global(myid, 1, i+1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), 3, i+1, j)
                    elseif((i==localncols).and.(pcol/=pncols)) then
                        call index_convert_local_global(myid, 3, localncols, j, fieldInd)
                        call index_convert_local_global(myid, 1, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 3, localncols, j)
                    elseif((i==localncols).and.(pcol==pncols)) then
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        call index_convert_local_global(myid, 1, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 3, i, j)
                        if(isDiriX1_p(ylower+j-1) /= 0) then
                            call index_convert_local_global(myid, 1, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), 3, i+1, j)
                        end if
                    else
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        call index_convert_local_global(myid, 1, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 3, i, j)
                        call index_convert_local_global(myid, 1, i+1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), 3, i+1, j)
                    end if
                end if

            end do
        end do

    end subroutine constructAxp

    subroutine constructAyy(vely, resi, sd_kind)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: sd_kind

        integer :: m_kind, fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu, indextemp
        integer :: i, j

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

        if(prow /= pnrows) then
            indextemp = localnrows
        else
            indextemp = localnrows + 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(vely(i,j) == 1) then

                    if((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid-pncols, 2, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 2, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, 1)
                    elseif((prow/=pnrows).and.(j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid+pncols, 2, i, 1, fieldInd)
                        call index_convert_local_global(myid, 2, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, localnrows)
                    elseif((i==0).and.(j/=0).and.(j/=indextemp+1).and.(.not.((prow==1).and.(j==1).and. &!
                        (isDiriY0_p(xlower+i)==0))).and.(.not.((prow==pnrows).and.(j==indextemp).and. &!
                        (isDiriY1_p(xlower+i)==0)))) then
                        call index_convert_local_global(myid-1, 2, localncols, j, fieldInd)
                        call index_convert_local_global(myid, 2, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, 1, j)
                    elseif((i==localncols+1).and.(j/=0).and.(j/=indextemp+1).and.(.not.((prow==1).and.(j==1).and. &!
                        (isDiriY0_p(xlower+i-2)==0))).and.(.not.((prow==pnrows).and.(j==indextemp).and. &!
                        (isDiriY1_p(xlower+i-2)==0)))) then
                        call index_convert_local_global(myid+1, 2, 1, j, fieldInd)
                        call index_convert_local_global(myid, 2, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, localncols, j)
                    elseif((j>=1).and.(j<=indextemp).and.(i>=1).and.(i<=localncols)) then
                        call index_convert_local_global(myid, 2, i, j, fieldInd)
                        if((j>=3).or.((j==2).and.(prow/=1)).or.((j==2).and.(prow==1).and.(isDiriY0_p(xlower+i-1)/=0))) then
                            call index_convert_local_global(myid, 2, i, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, j-1)
                        end if
                        if(((i/=1).and.(j>=2).and.(j<=indextemp-1)).or.((i/=1).and.(j==1).and.(prow/=1)).or.((i/=1).and. &!
                            (j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-2)/=0)).or.((i/=1).and.(j==indextemp).and. &!
                            (prow/=pnrows)).or.((i/=1).and.(j==indextemp).and.(prow==pnrows).and. &!
                            (isDiriY1_p(xlower+i-2)/=0))) then
                            call index_convert_local_global(myid, 2, i-1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, i-1, j)
                        end if
                        call setMatValue(fieldInd, fieldInd, resi(i,j), m_kind, i, j)
                        if(((i/=localncols).and.(j>=2).and.(j<=indextemp-1)).or.((i/=localncols).and.(j==1).and.(prow/=1)).or. &!
                            ((i/=localncols).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i)/=0)).or.((i/=localncols).and. &!
                            (j==indextemp).and.(prow/=pnrows)).or.((i/=localncols).and.(j==indextemp).and.(prow==pnrows) &!
                            .and.(isDiriY1_p(xlower+i)/=0))) then
                            call index_convert_local_global(myid, 2, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, i+1, j)
                        end if
                        if((j<=indextemp-2).or.((j==indextemp-1).and.(prow/=pnrows)).or.((j==indextemp-1).and. &!
                            (prow==pnrows).and.(isDiriY1_p(xlower+i-1)/=0))) then
                            call index_convert_local_global(myid, 2, i, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine constructAyy

    subroutine constructAyp(pres, resi)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        
        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
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

                if(pres(i,j) == 1) then
                    if(j == 0) then
                        call index_convert_local_global(myid-pncols, 3, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 2, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), 6, i, 1)
                    elseif((j==1).and.(prow==1)) then
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        if(isDiriY0_p(xlower+i-1) /= 0) then
                            call index_convert_local_global(myid, 2, i, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j), 6, i, j)
                        end if
                        call index_convert_local_global(myid, 2, i, j+1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), 6, i, j+1)
                    elseif((j==localnrows).and.(prow/=pnrows)) then
                        call index_convert_local_global(myid, 3, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 2, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 6, i, localnrows)
                    elseif((j==localnrows).and.(prow==pnrows)) then
                        call index_convert_local_global(myid, 3, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 2, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 6, i, localnrows)
                        if(isDiriY1_p(xlower+i-1) /= 0) then
                            call index_convert_local_global(myid, 2, i, localnrows+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), 6, i, localnrows+1)
                        end if
                    else
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        call index_convert_local_global(myid, 2, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 6, i, j)
                        call index_convert_local_global(myid, 2, i, j+1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), 6, i, j+1)
                    end if
                end if

            end do
        end do

    end subroutine constructAyp

    subroutine constructAcx(velx, resi)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        
        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        ! the field index
        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                if(velx(i,j) == 1) then
                    if(i == 1) then
                        call index_convert_local_global(myid, 1, 1, j, fieldInd)
                        call index_convert_local_global(myid, 3, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 7, 1, j)
                    elseif((i==localncols+1).and.(pcol/=pncols)) then
                        call index_convert_local_global(myid+1, 1, 1, j, fieldInd)
                        call index_convert_local_global(myid, 3, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), 7, localncols, j)
                    elseif((i==localncols+1).and.(pcol==pncols)) then
                        call index_convert_local_global(myid, 1, localncols+1, j, fieldInd)
                        call index_convert_local_global(myid, 3, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), 7, localncols, j)
                    else
                        call index_convert_local_global(myid, 1, i, j, fieldInd)
                        call index_convert_local_global(myid, 3, i-1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), 7, i-1, j)
                        call index_convert_local_global(myid, 3, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 7, i, j)
                    end if
                end if

            end do
        end do

    end subroutine constructAcx

    subroutine constructAcy(vely, resi)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        
        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        ! the field index
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1

        do j = indexd, indexu
            do i = indexl, indexr

                if(vely(i,j) == 1) then
                    if(j == 1) then
                        call index_convert_local_global(myid, 2, i, 1, fieldInd)
                        call index_convert_local_global(myid, 3, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 8, i, 1)
                    elseif((j==localnrows+1).and.(prow/=pnrows)) then
                        call index_convert_local_global(myid+pncols, 2, i, 1, fieldInd)
                        call index_convert_local_global(myid, 3, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), 8, i, localnrows)
                    elseif((j==localnrows+1).and.(prow==pnrows)) then
                        call index_convert_local_global(myid, 2, i, localnrows+1, fieldInd)
                        call index_convert_local_global(myid, 3, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), 8, i, localnrows)
                    else
                        call index_convert_local_global(myid, 2, i, j, fieldInd)
                        call index_convert_local_global(myid, 3, i, j-1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), 8, i, j-1)
                        call index_convert_local_global(myid, 3, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 8, i, j)
                    end if
                end if

            end do
        end do

    end subroutine constructAcy

    subroutine constructAcp(pres, resi)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        
        integer :: fieldInd
        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                if(pres(i,j) == 1) then
                    call index_convert_local_global(myid, 3, i, j, fieldInd)
                    call setMatValue(fieldInd, fieldInd, resi(i,j), 9, i, j)
                end if

            end do
        end do

    end subroutine constructAcp

    subroutine constructAcf(conc, resi)

        integer, dimension(:,:), pointer, intent(in) :: conc
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        
        integer :: fieldInd, equInd
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

                if(conc(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid-1, 4, localncols, j, fieldInd)
                        call index_convert_local_global(myid, 4, i+1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), 10, i+1, j)
                        if(j-1 >= 1) then
                            call index_convert_local_global(myid, 4, i+1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j-1), 10, i+1, j-1)
                        end if
                        if(j+1 <= localnrows) then
                            call index_convert_local_global(myid, 4, i+1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j+1), 10, i+1, j+1)
                        end if
                    elseif((i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid+1, 4, 1, j, fieldInd)
                        call index_convert_local_global(myid, 4, i-1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), 10, i-1, j)
                        if(j-1 >= 1) then
                            call index_convert_local_global(myid, 4, i-1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j-1), 10, i-1, j-1)
                        end if
                        if(j+1 <= localnrows) then
                            call index_convert_local_global(myid, 4, i-1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j+1), 10, i-1, j+1)
                        end if
                    elseif((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid-pncols, 4, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 4, i, j+1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), 10, i, j+1)
                        if(i-1 >= 1) then
                            call index_convert_local_global(myid, 4, i-1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j+1), 10, i-1, j+1)
                        end if
                        if(i+1 <= localncols) then
                            call index_convert_local_global(myid, 4, i+1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j+1), 10, i+1, j+1)
                        end if
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid+pncols, 4, i, 1, fieldInd)
                        call index_convert_local_global(myid, 4, i, j-1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), 10, i, j-1)
                        if(i-1 >= 1) then
                            call index_convert_local_global(myid, 4, i-1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j-1), 10, i-1, j-1)
                        end if
                        if(i+1 <= localncols) then
                            call index_convert_local_global(myid, 4, i+1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j-1), 10, i+1, j-1)
                        end if
                    elseif((i==0).and.(j==0)) then
                        call index_convert_local_global(myid-1-pncols, 4, localncols, localnrows, fieldInd)
                        call index_convert_local_global(myid, 4, 1, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(1,1), 10, 1, 1)
                    elseif((i==0).and.(j==localnrows+1)) then
                        call index_convert_local_global(myid-1+pncols, 4, localncols, 1, fieldInd)
                        call index_convert_local_global(myid, 4, 1, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(1,localnrows), 10, 1, localnrows)
                    elseif((i==localncols+1).and.(j==0)) then
                        call index_convert_local_global(myid+1-pncols, 4, 1, localnrows, fieldInd)
                        call index_convert_local_global(myid, 4, localncols, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(localncols,1), 10, localncols, 1)
                    elseif((i==localncols+1).and.(j==localnrows+1)) then
                        call index_convert_local_global(myid+1+pncols, 4, 1, 1, fieldInd)
                        call index_convert_local_global(myid, 4, localncols, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(localncols,localnrows), 10, localncols, localnrows)
                    else
                        call index_convert_local_global(myid, 4, i, j, fieldInd)
                        call index_convert_local_global(myid, 4, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 10, i, j)
                        if((i-1>=1).and.(j-1>=1)) then
                            call index_convert_local_global(myid, 4, i-1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j-1), 10, i-1, j-1)
                        end if
                        if(j-1 >= 1) then
                            call index_convert_local_global(myid, 4, i, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1), 10, i, j-1)
                        end if
                        if((i+1<=localncols).and.(j-1>=1)) then
                            call index_convert_local_global(myid, 4, i+1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j-1), 10, i+1, j-1)
                        end if
                        if(i-1 >= 1) then
                            call index_convert_local_global(myid, 4, i-1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j), 10, i-1, j)
                        end if
                        if(i+1 <= localncols) then
                            call index_convert_local_global(myid, 4, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), 10, i+1, j)
                        end if
                        if((i-1>=1).and.(j+1<=localnrows)) then
                            call index_convert_local_global(myid, 4, i-1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j+1), 10, i-1, j+1)
                        end if
                        if(j+1 <= localnrows) then
                            call index_convert_local_global(myid, 4, i, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), 10, i, j+1)
                        end if
                        if((i+1<=localncols).and.(j+1<=localnrows)) then
                            call index_convert_local_global(myid, 4, i+1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j+1), 10, i+1, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine constructAcf

    subroutine setMatValue(col, row, value, m_kind, eq_i, eq_j)

        integer, intent(in) :: col ! global column index
        integer, intent(in) :: row ! global row index
        real(kind=8), intent(in) :: value
        integer, intent(in) :: m_kind ! matrix kind
        integer, intent(in) :: eq_i ! equation x-direction coordinate
        integer, intent(in) :: eq_j ! equation y-direction coordinate

        integer, dimension(:), pointer :: Acols
        integer, dimension(:), pointer :: Arows
        real(kind=8), dimension(:), pointer :: Avalues
        integer, dimension(:), pointer :: AEntryBase
        integer, dimension(:), pointer :: AEntryNum

        integer :: eq_kind ! equation kind
        integer :: indexr, pos, base, tail, shend, m ,n
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
        elseif(m_kind == 7) then ! Acx
            Acols => AcxCols
            Arows => AcxRows
            Avalues => AcxValues
            AEntryBase => AcxEntryBase
            AEntryNum => AcxEntryNum
            eq_kind = 3
        elseif(m_kind == 8) then ! Acy
            Acols => AcyCols
            Arows => AcyRows
            Avalues => AcyValues
            AEntryBase => AcyEntryBase
            AEntryNum => AcyEntryNum
            eq_kind = 3
        elseif(m_kind == 9) then ! Acp
            Acols => AcpCols
            Arows => AcpRows
            Avalues => AcpValues
            AEntryBase => AcpEntryBase
            AEntryNum => AcpEntryNum
            eq_kind = 3
        elseif(m_kind == 10) then ! Acf
            Acols => AcfCols
            Arows => AcfRows
            Avalues => AcfValues
            AEntryBase => AcfEntryBase
            AEntryNum => AcfEntryNum
            eq_kind = 4
        end if

        ! equation kind
        if(eq_kind == 1) then
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


      !  isExist = .false.
       ! base = AEntryBase(pos)
        !tail = base + AEntryNum(pos) - 1
       ! do n = base, tail
        !    if(Acols(n) == col) then
         !       Avalues(n) = value
          !      isExist = .true.
           !     exit
         !   end if
       ! end do

        !if(.not.isExist) then
         !   do n = base, tail
          !      if(Acols(n) == 0) then
           !         Acols(n) = col
            !        Arows(n) = row
             !       Avalues(n) = value
              !      exit
               ! end if
           ! end do
        !end if

    end subroutine setMatValue

    ! change the index of the unknowns from the local index to the global index
    subroutine index_convert_local_global(pid, kind, local_i, local_j, global_ind)

        integer, intent(in) :: pid
        integer, intent(in) :: kind
        integer, intent(in) :: local_i
        integer, intent(in) :: local_j
        integer, intent(out) :: global_ind

        integer :: p_pcol, p_prow
        integer :: base, ubase, vbase, Cfbase

        p_pcol = mod(pid,pncols)+1
        p_prow = pid/pncols+1

        base = (p_prow-1)*((nx+1)*localnrows+nx*localnrows+nx*localnrows)
        base = base + (p_pcol-1)*(localncols*localnrows*3)
        if(p_prow == pnrows) then
            base = base + (p_pcol-1)*localncols
        end if

        if(p_pcol == pncols) then
            ubase = localnrows*(localncols+1)
        else
            ubase = localnrows*localncols
        end if

        if(p_prow == pnrows) then
            vbase = (localnrows+1)*localncols
        else
            vbase = localnrows*localncols
        end if

        Cfbase = (p_prow-1)*nx*localnrows + (p_pcol-1)*localncols*localnrows

        ! the '|' index
        if(kind == 1) then
            if(p_pcol == pncols) then
                global_ind = base + (local_j-1)*(localncols+1) + local_i
            else
                global_ind = base + (local_j-1)*localncols + local_i
            end if
        ! the '-' index
        elseif(kind == 2) then
            global_ind = base + ubase + (local_j-1)*localncols + local_i
        ! the '.' index
        elseif(kind == 3) then
            global_ind = base + ubase + vbase + (local_j-1)*localncols+local_i
        ! the 'Cf' index
        elseif(kind == 4) then
            global_ind = Cfbase + (local_j-1)*localncols + local_i
        end if

    end subroutine index_convert_local_global

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
