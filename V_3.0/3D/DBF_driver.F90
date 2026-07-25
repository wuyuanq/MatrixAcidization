
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_driver

    use DBF_model
    use DBF_globalData
    use DBF_resi
    use DBF_constructMat
    use DBF_exportResults

    implicit none
    include 'mpif.h'

contains

    subroutine genRandomNum()

        real(kind=8) :: rand
        real(kind=8), dimension(:), allocatable :: random
        integer :: ierr, i

        open(unit=10, file=trim(adjustl(FRANDOMTXT)), status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file error. ', ierr
            stop
        end if

        allocate(random(RANDOMSIZE))
        i = 1
        do while (i <= RANDOMSIZE)
            call random_number(rand)
            rand = 6.D-2*rand + 1.5D-1
            random(i) = rand
            i = i + 1
        end do
        write(10, fmt="(f8.6)") random(:)
        deallocate(random)
        close(10)

    end subroutine genRandomNum

    subroutine computeMatEntryNum()

        integer :: indexr, indexu, indexb
        integer :: eq_ind
        integer :: i, j, k

        AxxEntryNum_standard = 7
        AyyEntryNum_standard = 7
        AzzEntryNum_standard = 7
        ApEntryNum_standard = 7
        AcfEntryNum_standard = 19
        AtemEntryNum_standard = 7

        AxxEntryNum(:) = AxxEntryNum_standard
        AyyEntryNum(:) = AyyEntryNum_standard
        AzzEntryNum(:) = AzzEntryNum_standard
        ApEntryNum(:) = ApEntryNum_standard
        AcfEntryNum(:) = AcfEntryNum_standard
        AtemEntryNum(:) = AtemEntryNum_standard

        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call coordiToLocalInd(151, 1, j, k, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call coordiToLocalInd(151, localncols+1, j, k, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, indexr
                    call coordiToLocalInd(151, i, 1, k, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, indexr
                    call coordiToLocalInd(151, i, localnrows, k, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, indexr
                    call coordiToLocalInd(151, i, j, 1, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, indexr
                    call coordiToLocalInd(151, i, j, localnlays, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        do k = 1, localnlays
            do j = 1, localnrows
                if((pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                    call coordiToLocalInd(151, 1, j, k, eq_ind)
                    AxxEntryNum(eq_ind) = 1
                end if
                if((pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                    call coordiToLocalInd(151, localncols+1, j, k, eq_ind)
                    AxxEntryNum(eq_ind) = 1
                end if
            end do
        end do
       
        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if
        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, indexu
                    call coordiToLocalInd(253, 1, j, k, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, indexu
                    call coordiToLocalInd(253, localncols, j, k, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    call coordiToLocalInd(253, i, 1, k, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    call coordiToLocalInd(253, i, localnrows+1, k, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, indexu
                do i = 1, localncols
                    call coordiToLocalInd(253, i, j, 1, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, indexu
                do i = 1, localncols
                    call coordiToLocalInd(253, i, j, localnlays, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        do k = 1, localnlays
            do i = 1, localncols
                if((prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                    call coordiToLocalInd(253, i, 1, k, eq_ind)
                    AyyEntryNum(eq_ind) = 1
                end if
                if((prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                    call coordiToLocalInd(253, i, localnrows+1, k, eq_ind)
                    AyyEntryNum(eq_ind) = 1
                end if
            end do
        end do

        if(play /= pnlays) then
            indexb = localnlays
        else
            indexb = localnlays + 1
        end if
        if(pcol == 1) then
            do k = 1, indexb
                do j = 1, localnrows
                    call coordiToLocalInd(355, 1, j, k, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, indexb
                do j = 1, localnrows
                    call coordiToLocalInd(355, localncols, j, k, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, indexb
                do i = 1, localncols
                    call coordiToLocalInd(355, i, 1, k, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, indexb
                do i = 1, localncols
                    call coordiToLocalInd(355, i, localnrows, k, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    call coordiToLocalInd(355, i, j, 1, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    call coordiToLocalInd(355, i, j, localnlays+1, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        do j = 1, localnrows
            do i = 1, localncols
                if((play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                    call coordiToLocalInd(355, i, j, 1, eq_ind)
                    AzzEntryNum(eq_ind) = 1
                end if
                if((play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                    call coordiToLocalInd(355, i, j, localnlays+1, eq_ind)
                    AzzEntryNum(eq_ind) = 1
                end if
            end do
        end do

        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call coordiToLocalInd(457, 1, j, k, eq_ind)
                    ApEntryNum(eq_ind) = ApEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call coordiToLocalInd(457, localncols, j, k, eq_ind)
                    ApEntryNum(eq_ind) = ApEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    call coordiToLocalInd(457, i, 1, k, eq_ind)
                    ApEntryNum(eq_ind) = ApEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    call coordiToLocalInd(457, i, localnrows, k, eq_ind)
                    ApEntryNum(eq_ind) = ApEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    call coordiToLocalInd(457, i, j, 1, eq_ind)
                    ApEntryNum(eq_ind) = ApEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    call coordiToLocalInd(457, i, j, localnlays, eq_ind)
                    ApEntryNum(eq_ind) = ApEntryNum(eq_ind) - 1
                end do
            end do
        end if

        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call coordiToLocalInd(458, 1, j, k, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call coordiToLocalInd(458, localncols, j, k, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    call coordiToLocalInd(458, i, 1, k, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    call coordiToLocalInd(458, i, localnrows, k, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    call coordiToLocalInd(458, i, j, 1, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    call coordiToLocalInd(458, i, j, localnlays, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if((pcol==1).and.(prow==1)) then
            do k = 1, localnlays
                call coordiToLocalInd(458, 1, 1, k, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==1).and.(prow==pnrows)) then
            do k = 1, localnlays
                call coordiToLocalInd(458, 1, localnrows, k, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==pncols).and.(prow==1)) then
            do k = 1, localnlays
                call coordiToLocalInd(458, localncols, 1, k, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==pncols).and.(prow==pnrows)) then
            do k = 1, localnlays
                call coordiToLocalInd(458, localncols, localnrows, k, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==1).and.(play==1)) then
            do j = 1, localnrows
                call coordiToLocalInd(458, 1, j, 1, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==1).and.(play==pnlays)) then
            do j = 1, localnrows
                call coordiToLocalInd(458, 1, j, localnlays, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==pncols).and.(play==1)) then
            do j = 1, localnrows
                call coordiToLocalInd(458, localncols, j, 1, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==pncols).and.(play==pnlays)) then
            do j = 1, localnrows
                call coordiToLocalInd(458, localncols, j, localnlays, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((prow==1).and.(play==1)) then
            do i = 1, localncols
                call coordiToLocalInd(458, i, 1, 1, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((prow==1).and.(play==pnlays)) then
            do i = 1, localncols
                call coordiToLocalInd(458, i, 1, localnlays, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((prow==pnrows).and.(play==1)) then
            do i = 1, localncols
                call coordiToLocalInd(458, i, localnrows, 1, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((prow==pnrows).and.(play==pnlays)) then
            do i = 1, localncols
                call coordiToLocalInd(458, i, localnrows, localnlays, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if

        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call coordiToLocalInd(459, 1, j, k, eq_ind)
                    AtemEntryNum(eq_ind) = AtemEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call coordiToLocalInd(459, localncols, j, k, eq_ind)
                    AtemEntryNum(eq_ind) = AtemEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    call coordiToLocalInd(459, i, 1, k, eq_ind)
                    AtemEntryNum(eq_ind) = AtemEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    call coordiToLocalInd(459, i, localnrows, k, eq_ind)
                    AtemEntryNum(eq_ind) = AtemEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    call coordiToLocalInd(459, i, j, 1, eq_ind)
                    AtemEntryNum(eq_ind) = AtemEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    call coordiToLocalInd(459, i, j, localnlays, eq_ind)
                    AtemEntryNum(eq_ind) = AtemEntryNum(eq_ind) - 1
                end do
            end do
        end if

    end subroutine computeMatEntryNum

    subroutine communicateP(scalN)

        ! 1 means dm, 2 means porosity, 3 means permeability, 4 means pressure
        integer, intent(in) :: scalN

        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu, requestf, requestb, requestlu, &!
            requestlb, requestrd, requestru, requestrf, requestrb, requestdb, requestuf, requestub
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, k, c, ierr

        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(pcol /= 1) then
                if(scalN == 1) then
                    sentSize = localnrows*localnlays
                elseif(scalN == 2) then
                    sentSize = localnrows*localnlays*2
                elseif(scalN == 3) then
                    sentSize = localnrows*localnlays*3
                end if
                allocate(sent(sentSize))
             c = 0
                do k = 1, localnlays
                    do j = 1, localnrows
                        if(scalN == 1) then
                            c = c + 1
                            sent(c) = dm(1,j,k)
                        elseif(scalN == 2) then
                            c = c + 1
                            sent(c) = poro(1,j,k)
                            c = c + 1
                            sent(c) = poro_old(1,j,k)
                        elseif(scalN == 3) then
                            c = c + 1
                            sent(c) = Kxx(1,j,k)
                            c = c + 1
                            sent(c) = Kyy(1,j,k)
                            c = c + 1
                            sent(c) = Kzz(1,j,k)
                        end if
                    end do
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
                deallocate(sent)
            end if
        end if

        if(pcol /= pncols) then
            if(scalN == 1) then
                sentSize = localnrows*localnlays
            elseif(scalN == 2) then
                sentSize = localnrows*localnlays*4
            elseif(scalN == 3) then
                sentSize = localnrows*localnlays*3
            elseif(scalN == 4) then
                sentSize = localnrows*localnlays
            end if
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    if(scalN == 1) then
                        c = c + 1
                        sent(c) = dm(localncols,j,k)
                    elseif(scalN == 2) then
                        c = c + 1
                        sent(c) = poro(localncols-1,j,k)
                        c = c + 1
                        sent(c) = poro_old(localncols-1,j,k)
                        c = c + 1
                        sent(c) = poro(localncols,j,k)
                        c = c + 1
                        sent(c) = poro_old(localncols,j,k)
                    elseif(scalN == 3) then
                        c = c + 1
                        sent(c) = Kxx(localncols,j,k)
                        c = c + 1
                        sent(c) = Kyy(localncols,j,k)
                        c = c + 1
                        sent(c) = Kzz(localncols,j,k)
                    elseif(scalN == 4) then
                        c = c + 1
                        sent(c) = p(localncols,j,k)
                    end if
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(prow /= 1) then
                if(scalN == 1) then
                    sentSize = localncols*localnlays
                elseif(scalN == 2) then
                    sentSize = localncols*localnlays*2
                elseif(scalN == 3) then
                    sentSize = localncols*localnlays*3
                end if
                allocate(sent(sentSize))
                c = 0
                do k = 1, localnlays
                    do i = 1, localncols
                        if(scalN == 1) then
                            c = c + 1
                            sent(c) = dm(i,1,k)
                        elseif(scalN == 2) then
                            c = c + 1
                            sent(c) = poro(i,1,k)
                            c = c + 1
                            sent(c) = poro_old(i,1,k)
                        elseif(scalN == 3) then
                            c = c + 1
                            sent(c) = Kxx(i,1,k)
                            c = c + 1
                            sent(c) = Kyy(i,1,k)
                            c = c + 1
                            sent(c) = Kzz(i,1,k)
                        end if
                    end do
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
                deallocate(sent)
            end if
        end if

        if(prow /= pnrows) then
            if(scalN == 1) then
                sentSize = localncols*localnlays
            elseif(scalN == 2) then
                sentSize = localncols*localnlays*4
            elseif(scalN == 3) then
                sentSize = localncols*localnlays*3
            elseif(scalN == 4) then
                sentSize = localncols*localnlays
            end if
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do i = 1, localncols
                    if(scalN == 1) then
                        c = c + 1
                        sent(c) = dm(i,localnrows,k)
                    elseif(scalN == 2) then
                        c = c + 1
                        sent(c) = poro(i,localnrows-1,k)
                        c = c + 1
                        sent(c) = poro_old(i,localnrows-1,k)
                        c = c + 1
                        sent(c) = poro(i,localnrows,k)
                        c = c + 1
                        sent(c) = poro_old(i,localnrows,k)
                    elseif(scalN == 3) then
                        c = c + 1
                        sent(c) = Kxx(i,localnrows,k)
                        c = c + 1
                        sent(c) = Kyy(i,localnrows,k)
                        c = c + 1
                        sent(c) = Kzz(i,localnrows,k)
                    elseif(scalN == 4) then
                        c = c + 1
                        sent(c) = p(i,localnrows,k)
                    end if
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(play /= 1) then
                if(scalN == 1) then
                    sentSize = localncols*localnrows
                elseif(scalN == 2) then
                    sentSize = localncols*localnrows*2
                elseif(scalN == 3) then
                    sentSize = localncols*localnrows*3
                end if
                allocate(sent(sentSize))
                c = 0
                do j = 1, localnrows
                    do i = 1, localncols
                        if(scalN == 1) then
                            c = c + 1
                            sent(c) = dm(i,j,1)
                        elseif(scalN == 2) then
                            c = c + 1
                            sent(c) = poro(i,j,1)
                            c = c + 1
                            sent(c) = poro_old(i,j,1)
                        elseif(scalN == 3) then
                            c = c + 1
                            sent(c) = Kxx(i,j,1)
                            c = c + 1
                            sent(c) = Kyy(i,j,1)
                            c = c + 1
                            sent(c) = Kzz(i,j,1)
                        end if
                    end do
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid-pncols*pnrows, myid, MPI_COMM_WORLD, requestf, ierr)
                deallocate(sent)
            end if
        end if

        if(play /= pnlays) then
            if(scalN == 1) then
                sentSize = localncols*localnrows
            elseif(scalN == 2) then
                sentSize = localncols*localnrows*4
            elseif(scalN == 3) then
                sentSize = localncols*localnrows*3
            elseif(scalN == 4) then
                sentSize = localncols*localnrows
            end if
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                do i = 1, localncols
                    if(scalN == 1) then
                        c = c + 1
                        sent(c) = dm(i,j,localnlays)
                    elseif(scalN == 2) then
                        c = c + 1
                        sent(c) = poro(i,j,localnlays-1)
                        c = c + 1
                        sent(c) = poro_old(i,j,localnlays-1)
                        c = c + 1
                        sent(c) = poro(i,j,localnlays)
                        c = c + 1
                        sent(c) = poro_old(i,j,localnlays)
                    elseif(scalN == 3) then
                        c = c + 1
                        sent(c) = Kxx(i,j,localnlays)
                        c = c + 1
                        sent(c) = Kyy(i,j,localnlays)
                        c = c + 1
                        sent(c) = Kzz(i,j,localnlays)
                    elseif(scalN == 4) then
                        c = c + 1
                        sent(c) = p(i,j,localnlays)
                    end if
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols*pnrows, myid, MPI_COMM_WORLD, requestb, ierr)
            deallocate(sent)
        end if

        if(scalN==2) then
            if((pcol/=1).and.(prow/=pnrows)) then
                sentSize = localnlays*2
                allocate(sent(sentSize))
                c = 0
                do k = 1, localnlays
                    c = c + 1
                    sent(c) = poro(1,localnrows,k)
                    c = c + 1
                    sent(c) = poro_old(1,localnrows,k)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid-1+pncols, myid, MPI_COMM_WORLD, requestlu, ierr)
                deallocate(sent)
            end if
        end if

        if(scalN==2) then
            if((pcol/=1).and.(play/=pnlays)) then
                sentSize = localnrows*2
                allocate(sent(sentSize))
                c = 0
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = poro(1,j,localnlays)
                    c = c + 1
                    sent(c) = poro_old(1,j,localnlays)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid-1+pncols*pnrows, myid, MPI_COMM_WORLD, requestlb, ierr)
                deallocate(sent)
            end if
        end if

        if(scalN==2) then
            if((pcol/=pncols).and.(prow/=1)) then
                sentSize = localnlays*2
                allocate(sent(sentSize))
                c = 0
                do k = 1, localnlays
                    c = c + 1
                    sent(c) = poro(localncols,1,k)
                    c = c + 1
                    sent(c) = poro_old(localncols,1,k)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid+1-pncols, myid, MPI_COMM_WORLD, requestrd, ierr)
                deallocate(sent)
            end if
        end if

        if(scalN==2) then
            if((pcol/=pncols).and.(prow/=pnrows)) then
                sentSize = localnlays*2
                allocate(sent(sentSize))
                c = 0
                do k = 1, localnlays
                    c = c + 1
                    sent(c) = poro(localncols,localnrows,k)
                    c = c + 1
                    sent(c) = poro_old(localncols,localnrows,k)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid+1+pncols, myid, MPI_COMM_WORLD, requestru, ierr)
                deallocate(sent)
            end if
        end if

        if(scalN==2) then
            if((pcol/=pncols).and.(play/=1)) then
                sentSize = localnrows*2
                allocate(sent(sentSize))
                c = 0
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = poro(localncols,j,1)
                    c = c + 1
                    sent(c) = poro_old(localncols,j,1)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid+1-pncols*pnrows, myid, MPI_COMM_WORLD, requestrf, ierr)
                deallocate(sent)
            end if
        end if

        if(scalN==2) then
            if((pcol/=pncols).and.(play/=pnlays)) then
                sentSize = localnrows*2
                allocate(sent(sentSize))
                c = 0
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = poro(localncols,j,localnlays)
                    c = c + 1
                    sent(c) = poro_old(localncols,j,localnlays)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid+1+pncols*pnrows, myid, MPI_COMM_WORLD, requestrb, ierr)
                deallocate(sent)
            end if
        end if

        if(scalN==2) then
            if((prow/=1).and.(play/=pnlays)) then
                sentSize = localncols*2
                allocate(sent(sentSize))
                c = 0
                do i = 1, localncols
                    c = c + 1
                    sent(c) = poro(i,1,localnlays)
                    c = c + 1
                    sent(c) = poro_old(i,1,localnlays)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid-pncols+pncols*pnrows, myid, MPI_COMM_WORLD, requestdb, ierr)
                deallocate(sent)
            end if
        end if

        if(scalN==2) then
            if((prow/=pnrows).and.(play/=1)) then
                sentSize = localncols*2
                allocate(sent(sentSize))
                c = 0
                do i = 1, localncols
                    c = c + 1
                    sent(c) = poro(i,localnrows,1)
                    c = c + 1
                    sent(c) = poro_old(i,localnrows,1)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid+pncols-pncols*pnrows, myid, MPI_COMM_WORLD, requestuf, ierr)
                deallocate(sent)
            end if
        end if

        if(scalN==2) then
            if((prow/=pnrows).and.(play/=pnlays)) then
                sentSize = localncols*2
                allocate(sent(sentSize))
                c = 0
                do i = 1, localncols
                    c = c + 1
                    sent(c) = poro(i,localnrows,localnlays)
                    c = c + 1
                    sent(c) = poro_old(i,localnrows,localnlays)
                end do
                call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                    myid+pncols+pncols*pnrows, myid, MPI_COMM_WORLD, requestub, ierr)
                deallocate(sent)
            end if
        end if

        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(pcol /= pncols) then
                if(scalN == 1) then
                    recvSize = localnrows*localnlays
                elseif(scalN == 2) then
                    recvSize = localnrows*localnlays*2
                elseif(scalN == 3) then
                    recvSize = localnrows*localnlays*3
                end if
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
                c = 0
                do k = 1, localnlays
                    do j = 1, localnrows
                        if(scalN == 1) then
                            c = c + 1
                            dm(localncols+1,j,k) = recv(c)
                        elseif(scalN == 2) then
                            c = c + 1
                            poro(localncols+1,j,k) = recv(c)
                            c = c + 1
                            poro_old(localncols+1,j,k) = recv(c)
                        elseif(scalN == 3) then
                            c = c + 1
                            Kxx(localncols+1,j,k) = recv(c)
                            c = c + 1
                            Kyy(localncols+1,j,k) = recv(c)
                            c = c + 1
                            Kzz(localncols+1,j,k) = recv(c)
                        end if
                    end do
                end do
                deallocate(recv)
            end if
        end if

        if(pcol /= 1) then
            if(scalN == 1) then
                recvSize = localnrows*localnlays
            elseif(scalN == 2) then
                recvSize = localnrows*localnlays*4
            elseif(scalN == 3) then
                recvSize = localnrows*localnlays*3
            elseif(scalN == 4) then
                recvSize = localnrows*localnlays
            end if
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    if(scalN == 1) then
                        c = c + 1
                        dm(0,j,k) = recv(c)
                    elseif(scalN == 2) then
                        c = c + 1
                        poro(-1,j,k) = recv(c)
                        c = c + 1
                        poro_old(-1,j,k) = recv(c)
                        c = c + 1
                        poro(0,j,k) = recv(c)
                        c = c + 1
                        poro_old(0,j,k) = recv(c)
                    elseif(scalN == 3) then
                        c = c + 1
                        Kxx(0,j,k) = recv(c)
                        c = c + 1
                        Kyy(0,j,k) = recv(c)
                        c = c + 1
                        Kzz(0,j,k) = recv(c)
                    elseif(scalN == 4) then
                        c = c + 1
                        p(0,j,k) = recv(c)
                    end if
                end do
            end do
            deallocate(recv)
        end if

        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(prow /= pnrows) then
                if(scalN == 1) then
                    recvSize = localncols*localnlays
                elseif(scalN == 2) then
                    recvSize = localncols*localnlays*2
                elseif(scalN == 3) then
                    recvSize = localncols*localnlays*3
                end if
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
                c = 0
                do k = 1, localnlays
                    do i = 1, localncols
                        if(scalN == 1) then
                            c = c + 1
                            dm(i,localnrows+1,k) = recv(c)
                        elseif(scalN == 2) then
                            c = c + 1
                            poro(i,localnrows+1,k) = recv(c)
                            c = c + 1
                            poro_old(i,localnrows+1,k) = recv(c)
                        elseif(scalN == 3) then
                            c = c + 1
                            Kxx(i,localnrows+1,k) = recv(c)
                            c = c + 1
                            Kyy(i,localnrows+1,k) = recv(c)
                            c = c + 1
                            Kzz(i,localnrows+1,k) = recv(c)
                        end if
                    end do
                end do
                deallocate(recv)
            end if
        end if

        if(prow /= 1) then
            if(scalN == 1) then
                recvSize = localncols*localnlays
            elseif(scalN == 2) then
                recvSize = localncols*localnlays*4
            elseif(scalN == 3) then
                recvSize = localncols*localnlays*3
            elseif(scalN == 4) then
                recvSize = localncols*localnlays
            end if
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do i = 1, localncols
                    if(scalN == 1) then
                        c = c + 1
                        dm(i,0,k) = recv(c)
                    elseif(scalN == 2) then
                        c = c + 1
                        poro(i,-1,k) = recv(c)
                        c = c + 1
                        poro_old(i,-1,k) = recv(c)
                        c = c + 1
                        poro(i,0,k) = recv(c)
                        c = c + 1
                        poro_old(i,0,k) = recv(c)
                    elseif(scalN == 3) then
                        c = c + 1
                        Kxx(i,0,k) = recv(c)
                        c = c + 1
                        Kyy(i,0,k) = recv(c)
                        c = c + 1
                        Kzz(i,0,k) = recv(c)
                    elseif(scalN == 4) then
                        c = c + 1
                        p(i,0,k) = recv(c)
                    end if
                end do
            end do
            deallocate(recv)
        end if

        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(play /= pnlays) then
                if(scalN == 1) then
                    recvSize = localncols*localnrows
                elseif(scalN == 2) then
                    recvSize = localncols*localnrows*2
                elseif(scalN == 3) then
                    recvSize = localncols*localnrows*3
                end if
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid+pncols*pnrows, myid+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
                c = 0
                do j = 1, localnrows
                    do i = 1, localncols
                        if(scalN == 1) then
                            c = c + 1
                            dm(i,j,localnlays+1) = recv(c)
                        elseif(scalN == 2) then
                            c = c + 1
                            poro(i,j,localnlays+1) = recv(c)
                            c = c + 1
                            poro_old(i,j,localnlays+1) = recv(c)
                        elseif(scalN == 3) then
                            c = c + 1
                            Kxx(i,j,localnlays+1) = recv(c)
                            c = c + 1
                            Kyy(i,j,localnlays+1) = recv(c)
                            c = c + 1
                            Kzz(i,j,localnlays+1) = recv(c)
                        end if
                    end do
                end do
                deallocate(recv)
            end if
        end if

        if(play /= 1) then
            if(scalN == 1) then
                recvSize = localncols*localnrows
            elseif(scalN == 2) then
                recvSize = localncols*localnrows*4
            elseif(scalN == 3) then
                recvSize = localncols*localnrows*3
            elseif(scalN == 4) then
                recvSize = localncols*localnrows
            end if
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols*pnrows, myid-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                do i = 1, localncols
                    if(scalN == 1) then
                        c = c + 1
                        dm(i,j,0) = recv(c)
                    elseif(scalN == 2) then
                        c = c + 1
                        poro(i,j,-1) = recv(c)
                        c = c + 1
                        poro_old(i,j,-1) = recv(c)
                        c = c + 1
                        poro(i,j,0) = recv(c)
                        c = c + 1
                        poro_old(i,j,0) = recv(c)
                    elseif(scalN == 3) then
                        c = c + 1
                        Kxx(i,j,0) = recv(c)
                        c = c + 1
                        Kyy(i,j,0) = recv(c)
                        c = c + 1
                        Kzz(i,j,0) = recv(c)
                    elseif(scalN == 4) then
                        c = c + 1
                        p(i,j,0) = recv(c)
                    end if
                end do
            end do
            deallocate(recv)
        end if

        if(scalN==2) then
            if((pcol/=pncols).and.(prow/=1)) then
                recvSize = localnlays*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid+1-pncols, myid+1-pncols, MPI_COMM_WORLD, status, ierr)
                c = 0
                do k = 1, localnlays
                    c = c + 1
                    poro(localncols+1,0,k) = recv(c)
                    c = c + 1
                    poro_old(localncols+1,0,k) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if(scalN==2) then
            if((pcol/=pncols).and.(play/=1)) then
                recvSize = localnrows*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid+1-pncols*pnrows, myid+1-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
                c = 0
                do j = 1, localnrows
                    c = c + 1
                    poro(localncols+1,j,0) = recv(c)
                    c = c + 1
                    poro_old(localncols+1,j,0) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if(scalN==2) then
            if((pcol/=1).and.(prow/=pnrows)) then
                recvSize = localnlays*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid-1+pncols, myid-1+pncols, MPI_COMM_WORLD, status, ierr)
                c = 0
                do k = 1, localnlays
                    c = c + 1
                    poro(0,localnrows+1,k) = recv(c)
                    c = c + 1
                    poro_old(0,localnrows+1,k) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if(scalN==2) then
            if((pcol/=1).and.(prow/=1)) then
                recvSize = localnlays*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid-1-pncols, myid-1-pncols, MPI_COMM_WORLD, status, ierr)
                c = 0
                do k = 1, localnlays
                    c = c + 1
                    poro(0,0,k) = recv(c)
                    c = c + 1
                    poro_old(0,0,k) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if(scalN==2) then
            if((pcol/=1).and.(play/=pnlays)) then
                recvSize = localnrows*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid-1+pncols*pnrows, myid-1+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
                c = 0
                do j = 1, localnrows
                    c = c + 1
                    poro(0,j,localnlays+1) = recv(c)
                    c = c + 1
                    poro_old(0,j,localnlays+1) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if(scalN==2) then
            if((pcol/=1).and.(play/=1)) then
                recvSize = localnrows*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid-1-pncols*pnrows, myid-1-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
                c = 0
                do j = 1, localnrows
                    c = c + 1
                    poro(0,j,0) = recv(c)
                    c = c + 1
                    poro_old(0,j,0) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if(scalN==2) then
            if((prow/=pnrows).and.(play/=1)) then
                recvSize = localncols*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid+pncols-pncols*pnrows, myid+pncols-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
                c = 0
                do i = 1, localncols
                    c = c + 1
                    poro(i,localnrows+1,0) = recv(c)
                    c = c + 1
                    poro_old(i,localnrows+1,0) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if(scalN==2) then
            if((prow/=1).and.(play/=pnlays)) then
                recvSize = localncols*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid-pncols+pncols*pnrows, myid-pncols+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
                c = 0
                do i = 1, localncols
                    c = c + 1
                    poro(i,0,localnlays+1) = recv(c)
                    c = c + 1
                    poro_old(i,0,localnlays+1) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if(scalN==2) then
            if((prow/=1).and.(play/=1)) then
                recvSize = localncols*2
                allocate(recv(recvSize))
                call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                    myid-pncols-pncols*pnrows, myid-pncols-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
                c = 0
                do i = 1, localncols
                    c = c + 1
                    poro(i,0,0) = recv(c)
                    c = c + 1
                    poro_old(i,0,0) = recv(c)
                end do
                deallocate(recv)
            end if
        end if

        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(pcol /= 1) then
                call MPI_WAIT(requestl, status, ierr)
            end if
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(prow /= 1) then
                call MPI_WAIT(requestd, status, ierr)
            end if
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if((scalN==1).or.(scalN==2).or.(scalN==3)) then
            if(play /= 1) then
                call MPI_WAIT(requestf, status, ierr)
            end if
        end if
        if(play /= pnlays) then
            call MPI_WAIT(requestb, status, ierr)
        end if

        if(scalN==2) then
            if((pcol/=1).and.(prow/=pnrows)) then
                call MPI_WAIT(requestlu, status, ierr)
            end if
            if((pcol/=1).and.(play/=pnlays)) then
                call MPI_WAIT(requestlb, status, ierr)
            end if
            if((pcol/=pncols).and.(prow/=1)) then
                call MPI_WAIT(requestrd, status, ierr)
            end if
            if((pcol/=pncols).and.(prow/=pnrows)) then
                call MPI_WAIT(requestru, status, ierr)
            end if
            if((pcol/=pncols).and.(play/=1)) then
                call MPI_WAIT(requestrf, status, ierr)
            end if
            if((pcol/=pncols).and.(play/=pnlays)) then
                call MPI_WAIT(requestrb, status, ierr)
            end if
            if((prow/=1).and.(play/=pnlays)) then
                call MPI_WAIT(requestdb, status, ierr)
            end if
            if((prow/=pnrows).and.(play/=1)) then
                call MPI_WAIT(requestuf, status, ierr)
            end if
            if((prow/=pnrows).and.(play/=pnlays)) then
                call MPI_WAIT(requestub, status, ierr)
            end if
        end if

    end subroutine communicateP

    subroutine communicateV()

        integer :: indexr, indexu, indexb
        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu, requestf, requestb, &!
            requestld, requestlu, requestlf, requestlb, requestrd, requestrf, &!
            requestdf, requestdb, requestuf
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, k, c, ierr

        if(pcol /= 1) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            sentSize = localnrows*localnlays + indexu*localnlays + localnrows*indexb
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = vx(1,j,k)
                end do
            end do
            do k = 1, localnlays
                do j = 1, indexu
                    c = c + 1
                    sent(c) = vy(1,j,k)
                end do
            end do
            do k = 1, indexb
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = vz(1,j,k)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            sentSize = indexu*localnlays + localnrows*indexb
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do j = 1, indexu
                    c = c + 1
                    sent(c) = vy(localncols,j,k)
                end do
            end do
            do k = 1, indexb
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = vz(localncols,j,k)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            sentSize = indexr*localnlays + localncols*localnlays + localncols*indexb
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do i = 1, indexr
                    c = c + 1
                    sent(c) = vx(i,1,k)
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    sent(c) = vy(i,1,k)
                end do
            end do
            do k = 1, indexb
                do i = 1, localncols
                    c = c + 1
                    sent(c) = vz(i,1,k)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            sentSize = indexr*localnlays + localncols*indexb
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do i = 1, indexr
                    c = c + 1
                    sent(c) = vx(i,localnrows,k)
                end do
            end do
            do k = 1, indexb
                do i = 1, localncols
                    c = c + 1
                    sent(c) = vz(i,localnrows,k)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if(play /= 1) then
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
            sentSize = indexr*localnrows + localncols*indexu + localncols*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                do i = 1, indexr
                    c = c + 1
                    sent(c) = vx(i,j,1)
                end do
            end do
            do j = 1, indexu
                do i = 1, localncols
                    c = c + 1
                    sent(c) = vy(i,j,1)
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    sent(c) = vz(i,j,1)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols*pnrows, myid, MPI_COMM_WORLD, requestf, ierr)
            deallocate(sent)
        end if

        if(play /= pnlays) then
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
            sentSize = indexr*localnrows + localncols*indexu
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                do i = 1, indexr
                    c = c + 1
                    sent(c) = vx(i,j,localnlays)
                end do
            end do
            do j = 1, indexu
                do i = 1, localncols
                     c = c + 1
                    sent(c) = vy(i,j,localnlays)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols*pnrows, myid, MPI_COMM_WORLD, requestb, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=1)) then
            sentSize = 2*localnlays
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                c = c + 1
                sent(c) = vx(1,1,k)
            end do
            do k = 1, localnlays
                c = c + 1
                sent(c) = vy(1,1,k)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-1-pncols, myid, MPI_COMM_WORLD, requestld, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            sentSize = localnlays
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                c = c + 1
                sent(c) = vx(1,localnrows,k)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-1+pncols, myid, MPI_COMM_WORLD, requestlu, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(play/=1)) then
            sentSize = 2*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = vx(1,j,1)
            end do
            do j = 1, localnrows
                c = c + 1
                sent(c) = vz(1,j,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-1-pncols*pnrows, myid, MPI_COMM_WORLD, requestlf, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(play/=pnlays)) then
            sentSize = localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = vx(1,j,localnlays)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-1+pncols*pnrows, myid, MPI_COMM_WORLD, requestlb, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            sentSize = localnlays
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                c = c + 1
                sent(c) = vy(localncols,1,k)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+1-pncols, myid, MPI_COMM_WORLD, requestrd, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(play/=1)) then
            sentSize = localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = vz(localncols,j,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+1-pncols*pnrows, myid, MPI_COMM_WORLD, requestrf, ierr)
            deallocate(sent)
        end if

        if((prow/=1).and.(play/=1)) then
            sentSize = 2*localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = vy(i,1,1)
            end do
            do i = 1, localncols
                c = c + 1
                sent(c) = vz(i,1,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols-pncols*pnrows, myid, MPI_COMM_WORLD, requestdf, ierr)
            deallocate(sent)
        end if

        if((prow/=1).and.(play/=pnlays)) then
            sentSize = localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = vy(i,1,localnlays)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols+pncols*pnrows, myid, MPI_COMM_WORLD, requestdb, ierr)
            deallocate(sent)
        end if

        if((prow/=pnrows).and.(play/=1)) then
            sentSize = localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = vz(i,localnrows,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols-pncols*pnrows, myid, MPI_COMM_WORLD, requestuf, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            recvSize = localnrows*localnlays + indexu*localnlays + localnrows*indexb
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    vx(localncols+1,j,k) = recv(c)
                end do
            end do
            do k = 1, localnlays
                do j = 1, indexu
                    c = c + 1
                    vy(localncols+1,j,k) = recv(c)
                end do
            end do
            do k = 1, indexb
                do j = 1, localnrows
                    c = c + 1
                    vz(localncols+1,j,k) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            recvSize = indexu*localnlays + localnrows*indexb
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do j = 1, indexu
                    c = c + 1
                    vy(0,j,k) = recv(c)
                end do
            end do
            do k = 1, indexb
                do j = 1, localnrows
                    c = c + 1
                    vz(0,j,k) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            recvSize = indexr*localnlays + localncols*localnlays + localncols*indexb
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do i = 1, indexr
                    c = c + 1
                    vx(i,localnrows+1,k) = recv(c)
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    vy(i,localnrows+1,k) = recv(c)
                end do
            end do
            do k = 1, indexb
                do i = 1, localncols
                    c = c + 1
                    vz(i,localnrows+1,k) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            recvSize = indexr*localnlays + localncols*indexb
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do i = 1, indexr
                    c = c + 1
                    vx(i,0,k) = recv(c)
                end do
            end do
            do k = 1, indexb
                do i = 1, localncols
                    c = c + 1
                    vz(i,0,k) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(play /= pnlays) then
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
            recvSize = indexr*localnrows + localncols*indexu + localncols*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols*pnrows, myid+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                do i = 1, indexr
                    c = c + 1
                    vx(i,j,localnlays+1) = recv(c)
                end do
            end do
            do j = 1, indexu
                do i = 1, localncols
                    c = c + 1
                    vy(i,j,localnlays+1) = recv(c)
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    vz(i,j,localnlays+1) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(play /= 1) then
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
            recvSize = indexr*localnrows + localncols*indexu
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols*pnrows, myid-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                do i = 1, indexr
                    c = c + 1
                    vx(i,j,0) = recv(c)
                end do
            end do
            do j = 1, indexu
                do i = 1, localncols
                    c = c + 1
                    vy(i,j,0) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=pnrows)) then
            recvSize = 2*localnlays
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+1+pncols, myid+1+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                c = c + 1
                vx(localncols+1,localnrows+1,k) = recv(c)
            end do
            do k = 1, localnlays
                c = c + 1
                vy(localncols+1,localnrows+1,k) = recv(c)
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            recvSize = localnlays
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+1-pncols, myid+1-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                c = c + 1
                vx(localncols+1,0,k) = recv(c)
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(play/=pnlays)) then
            recvSize = 2*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+1+pncols*pnrows, myid+1+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                vx(localncols+1,j,localnlays+1) = recv(c)
            end do
            do j = 1, localnrows
                c = c + 1
                vz(localncols+1,j,localnlays+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(play/=1)) then
            recvSize = localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+1-pncols*pnrows, myid+1-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                vx(localncols+1,j,0) = recv(c)
            end do
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            recvSize = localnlays
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-1+pncols, myid-1+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                c = c + 1
                vy(0,localnrows+1,k) = recv(c)
            end do
            deallocate(recv)
        end if

        if((pcol/=1).and.(play/=pnlays)) then
            recvSize = localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-1+pncols*pnrows, myid-1+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                vz(0,j,localnlays+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if((prow/=pnrows).and.(play/=pnlays)) then
            recvSize = 2*localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols+pncols*pnrows, myid+pncols+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                vy(i,localnrows+1,localnlays+1) = recv(c)
            end do
            do i = 1, localncols
                c = c + 1
                vz(i,localnrows+1,localnlays+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if((prow/=pnrows).and.(play/=1)) then
            recvSize = localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols-pncols*pnrows, myid+pncols-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                vy(i,localnrows+1,0) = recv(c)
            end do
            deallocate(recv)
        end if

        if((prow/=1).and.(play/=pnlays)) then
            recvSize = localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols+pncols*pnrows, myid-pncols+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                vz(i,0,localnlays+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if(play /= 1) then
            call MPI_WAIT(requestf, status, ierr)
        end if
        if(play /= pnlays) then
            call MPI_WAIT(requestb, status, ierr)
        end if
        if((pcol/=1).and.(prow/=1)) then
            call MPI_WAIT(requestld, status, ierr)
        end if
        if((pcol/=1).and.(prow/=pnrows)) then
            call MPI_WAIT(requestlu, status, ierr)
        end if
        if((pcol/=1).and.(play/=1)) then
            call MPI_WAIT(requestlf, status, ierr)
        end if
        if((pcol/=1).and.(play/=pnlays)) then
            call MPI_WAIT(requestlb, status, ierr)
        end if
        if((prow/=1).and.(pcol/=pncols)) then
            call MPI_WAIT(requestrd, status, ierr)
        end if
        if((play/=1).and.(pcol/=pncols)) then
            call MPI_WAIT(requestrf, status, ierr)
        end if
        if((prow/=1).and.(play/=1)) then
            call MPI_WAIT(requestdf, status, ierr)
        end if
        if((prow/=1).and.(play/=pnlays)) then
            call MPI_WAIT(requestdb, status, ierr)
        end if
        if((play/=1).and.(prow/=pnrows)) then
            call MPI_WAIT(requestuf, status, ierr)
        end if

    end subroutine communicateV

    subroutine initialize(modelCase)

        type(model), intent(in out) :: modelCase

        integer :: xmomSize, ymomSize, zmomSize, massSize, concenSize, enerSize
        ! The arrays that will communicate between different functions must use pointer type instead of allocatalbe
        ! type. Pointer type can make sure that the subscripts of the arrays keep the same in the calling and called
        ! functions. However, if the subscripts of the arrays begin with 1, using allocatable type is also OK.
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: p_pcol, p_prow, p_play, global_ind_b, global_ind_e
        logical :: alive
        integer :: i, j, k, n, ierr

        call MPI_Init(ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nProcs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)

#if defined(MUMPS) || defined(HYPRE)
        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)
#endif

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timestart = MPI_Wtime()
        solvertime = 0.D0

        pncols = modelcase%pncols
        pnrows = modelcase%pnrows
        pnlays = modelCase%pnlays
        Lx = modelCase%Lx
        Ly = modelCase%Ly
        Lz = modelCase%Lz
        nx = modelCase%nx
        ny = modelCase%ny
        nz = modelCase%nz
        timeEnd = modelCase%timeEnd
        nt = modelCase%nt
        visc = modelCase%visc
        rhof = modelCase%rhof
        rhos = modelCase%rhos
        alpha = modelCase%alpha
        dmRef = modelCase%dmRef
        TemRef_dm = modelCase%TemRef_dm
        ksRef = modelCase%ksRef
        TemRef_ks = modelCase%TemRef_ks
        epslon = modelCase%epslon
        alphaOS = modelCase%alphaOS
        lamdaX = modelCase%lamdaX
        lamdaT = modelCase%lamdaT
        lamdaf = modelCase%lamdaf
        lamdas = modelCase%lamdas
        thetaf = modelCase%thetaf
        thetas = modelCase%thetas
        radiusInit = modelCase%radiusInit
        ShInfinity = modelCase%ShInfinity
        al = modelCase%al
        gravX = modelCase%gravX
        gravY = modelCase%gravY
        gravZ = modelCase%gravZ

        allocate(xs(nx+1))
        allocate(ys(ny+1))
        allocate(zs(nz+1))
        allocate(ts(nt+1))
        allocate(src(nx,ny,nz))
        allocate(poroInit(nx,ny,nz))
        allocate(KxxInit(nx,ny,nz))
        allocate(KyyInit(nx,ny,nz))
        allocate(KzzInit(nx,ny,nz))
        allocate(avInit(nx,ny,nz))
        allocate(vxBdryX0(ny,nz))
        allocate(vxBdryX1(ny,nz))
        allocate(vxBdryY0(nx+1,nz))
        allocate(vxBdryY1(nx+1,nz))
        allocate(vxBdryZ0(nx+1,ny))
        allocate(vxBdryZ1(nx+1,ny))
        allocate(vyBdryX0(ny+1,nz))
        allocate(vyBdryX1(ny+1,nz))
        allocate(vyBdryY0(nx,nz))
        allocate(vyBdryY1(nx,nz))
        allocate(vyBdryZ0(nx,ny+1))
        allocate(vyBdryZ1(nx,ny+1))
        allocate(vzBdryX0(ny,nz+1))
        allocate(vzBdryX1(ny,nz+1))
        allocate(vzBdryY0(nx,nz+1))
        allocate(vzBdryY1(nx,nz+1))
        allocate(vzBdryZ0(nx,ny))
        allocate(vzBdryZ1(nx,ny))
        allocate(isDiriX0_p(ny,nz))
        allocate(isDiriX1_p(ny,nz))
        allocate(isDiriY0_p(nx,nz))
        allocate(isDiriY1_p(nx,nz))
        allocate(isDiriZ0_p(nx,ny))
        allocate(isDiriZ1_p(nx,ny))
        allocate(pBdryX0(ny,nz))
        allocate(pBdryX1(ny,nz))
        allocate(pBdryY0(nx,nz))
        allocate(pBdryY1(nx,nz))
        allocate(pBdryZ0(nx,ny))
        allocate(pBdryZ1(nx,ny))
        allocate(pInit(nx,ny,nz))
        allocate(isDiriX0_Cf(ny,nz))
        allocate(isDiriX1_Cf(ny,nz))
        allocate(isDiriY0_Cf(nx,nz))
        allocate(isDiriY1_Cf(nx,nz))
        allocate(isDiriZ0_Cf(nx,ny))
        allocate(isDiriZ1_Cf(nx,ny))
        allocate(CfBdryX0(ny,nz))
        allocate(CfBdryX1(ny,nz))
        allocate(CfBdryY0(nx,nz))
        allocate(CfBdryY1(nx,nz))
        allocate(CfBdryZ0(nx,ny))
        allocate(CfBdryZ1(nx,ny))
        allocate(CfInit(nx,ny,nz))
        allocate(TemBdryX0(ny,nz))
        allocate(TemBdryX1(ny,nz))
        allocate(TemBdryY0(nx,nz))
        allocate(TemBdryY1(nx,nz))
        allocate(TemBdryZ0(nx,ny))
        allocate(TemBdryZ1(nx,ny))
        allocate(TemInit(nx,ny,nz))

        xs = modelCase%xs
        ys = modelCase%ys
        zs = modelCase%zs
        ts = modelCase%ts
        src = modelCase%src
        poroInit = modelCase%poroInit
        KxxInit = modelCase%KxxInit
        KyyInit = modelCase%KyyInit
        KzzInit = modelCase%KzzInit
        avInit = modelCase%avInit
        vxBdryX0 = modelCase%vxBdryX0
        vxBdryX1 = modelCase%vxBdryX1
        vxBdryY0 = modelCase%vxBdryY0
        vxBdryY1 = modelCase%vxBdryY1
        vxBdryZ0 = modelCase%vxBdryZ0
        vxBdryZ1 = modelCase%vxBdryZ1
        vyBdryX0 = modelCase%vyBdryX0
        vyBdryX1 = modelCase%vyBdryX1
        vyBdryY0 = modelCase%vyBdryY0
        vyBdryY1 = modelCase%vyBdryY1
        vyBdryZ0 = modelCase%vyBdryZ0
        vyBdryZ1 = modelCase%vyBdryZ1
        vzBdryX0 = modelCase%vzBdryX0
        vzBdryX1 = modelCase%vzBdryX1
        vzBdryY0 = modelCase%vzBdryY0
        vzBdryY1 = modelCase%vzBdryY1
        vzBdryZ0 = modelCase%vzBdryZ0
        vzBdryZ1 = modelCase%vzBdryZ1
        isDiriX0_p = modelCase%isDiriX0_p
        isDiriX1_p = modelCase%isDiriX1_p
        isDiriY0_p = modelCase%isDiriY0_p
        isDiriY1_p = modelCase%isDiriY1_p
        isDiriZ0_p = modelCase%isDiriZ0_p
        isDiriZ1_p = modelCase%isDiriZ1_p
        pBdryX0 = modelCase%pBdryX0
        pBdryX1 = modelCase%pBdryX1
        pBdryY0 = modelCase%pBdryY0
        pBdryY1 = modelCase%pBdryY1
        pBdryZ0 = modelCase%pBdryZ0
        pBdryZ1 = modelCase%pBdryZ1
        pInit = modelCase%pInit
        isDiriX0_Cf = modelCase%isDiriX0_Cf
        isDiriX1_Cf = modelCase%isDiriX1_Cf
        isDiriY0_Cf = modelCase%isDiriY0_Cf
        isDiriY1_Cf = modelCase%isDiriY1_Cf
        isDiriZ0_Cf = modelCase%isDiriZ0_Cf
        isDiriZ1_Cf = modelCase%isDiriZ1_Cf
        CfBdryX0 = modelCase%CfBdryX0
        CfBdryX1 = modelCase%CfBdryX1
        CfBdryY0 = modelCase%CfBdryY0
        CfBdryY1 = modelCase%CfBdryY1
        CfBdryZ0 = modelCase%CfBdryZ0
        CfBdryZ1 = modelCase%CfBdryZ1
        CfInit = modelCase%CfInit
        isDiriX0_Tem = modelCase%isDiriX0_Tem
        isDiriX1_Tem = modelCase%isDiriX1_Tem
        isDiriY0_Tem = modelCase%isDiriY0_Tem
        isDiriY1_Tem = modelCase%isDiriY1_Tem
        isDiriZ0_Tem = modelCase%isDiriZ0_Tem
        isDiriZ1_Tem = modelCase%isDiriZ1_Tem
        TemBdryX0 = modelCase%TemBdryX0
        TemBdryX1 = modelCase%TemBdryX1
        TemBdryY0 = modelCase%TemBdryY0
        TemBdryY1 = modelCase%TemBdryY1
        TemBdryZ0 = modelCase%TemBdryZ0
        TemBdryZ1 = modelCase%TemBdryZ1
        TemInit = modelCase%TemInit
        soludoc = modelCase%soludoc

        localncols = nx/pncols
        localnrows = ny/pnrows
        localnlays = nz/pnlays

        play = myid/(pnrows*pncols)+1
        prow = (myid-(play-1)*pnrows*pncols)/pncols+1
        pcol = (myid-(play-1)*pnrows*pncols)-(prow-1)*pncols+1

        xlower = (pcol-1)*localncols+1
        xupper = pcol*localncols
        ylower = (prow-1)*localnrows+1
        yupper = prow*localnrows
        zlower = (play-1)*localnlays+1
        zupper = play*localnlays

        call coordiToGlobalInd(myid, 111, 1, 1, 1, ilower_v)
        if(play /= pnlays) then
            call coordiToGlobalInd(myid, 311, localncols, localnrows, localnlays, iupper_v)
        else
            call coordiToGlobalInd(myid, 311, localncols, localnrows, localnlays+1, iupper_v)
        end if
        local_x_size_v = iupper_v - ilower_v + 1

        if((nProcs>1).and.(myid==0)) then
            allocate(slave_v_data_size(nProcs-1))
            do i = 1, nProcs-1
                call coordiToGlobalInd(i, 111, 1, 1, 1, global_ind_b)
                p_play = i/(pnrows*pncols)+1
                if(p_play /= pnlays) then
                    call coordiToGlobalInd(i, 311, localncols, localnrows, localnlays, global_ind_e)
                else
                    call coordiToGlobalInd(i, 311, localncols, localnrows, localnlays+1, global_ind_e)
                end if
                slave_v_data_size(i) = global_ind_e - global_ind_b + 1
            end do
        end if

        call coordiToGlobalInd(myid, 412, 1, 1, 1, ilower_p)
        call coordiToGlobalInd(myid, 412, localncols, localnrows, localnlays, iupper_p)
        local_x_size_p = iupper_p - ilower_p + 1

        call coordiToGlobalInd(myid, 413, 1, 1, 1, ilower_Cf)
        call coordiToGlobalInd(myid, 413, localncols, localnrows, localnlays, iupper_Cf)
        local_x_size_Cf = iupper_Cf - ilower_Cf + 1

        call coordiToGlobalInd(myid, 414, 1, 1, 1, ilower_Tem)
        call coordiToGlobalInd(myid, 414, localncols, localnrows, localnlays, iupper_Tem)
        local_x_size_Tem = iupper_Tem - ilower_Tem + 1

        t = 2

        ! initialize dm
        allocate(dm(0:localncols+1, 0:localnrows+1, 0:localnlays+1))

        ! initialize kc
        allocate(kc(1:localncols, 1:localnrows, 1:localnlays))

        ! initialize ks
        allocate(ks(1:localncols, 1:localnrows, 1:localnlays))

        ! initialize hx, hy, hz, poro
        if(pcol /= 1) then
            indexl = -1
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        if(prow /= 1) then
            indexd = -1
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if
        if(play /= 1) then
            indexf = -1
        else
            indexf = 1
        end if
        if(play /= pnlays) then
            indexb = localnlays + 1
        else
            indexb = localnlays
        end if

        allocate(hx(indexl:indexr))
        do i = indexl, indexr
            hx(i) = xs(xlower+i) - xs(xlower+i-1)
        end do

        allocate(hy(indexd:indexu))
        do j = indexd, indexu
            hy(j) = ys(ylower+j) - ys(ylower+j-1)
        end do
       
        allocate(hz(indexf:indexb))
        do k = indexf, indexb
            hz(k) = zs(zlower+k) - zs(zlower+k-1)
        end do

        allocate(poro(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(poro_old(indexl:indexr,indexd:indexu,indexf:indexb))
        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    poro(i,j,k) = poroInit(xlower+i-1,ylower+j-1,zlower+k-1)
                    poro_old(i,j,k) = poro(i,j,k)
                end do
            end do
        end do

        ! initialize poroHarmX
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
        allocate(poroHarmX(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(poroHarmX_old(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(poroHarmXInit(indexl:indexr,indexd:indexu,indexf:indexb))

        ! initialize poroHarmY
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
        allocate(poroHarmY(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(poroHarmY_old(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(poroHarmYInit(indexl:indexr,indexd:indexu,indexf:indexb))

        ! initialize poroHarmZ
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
        allocate(poroHarmZ(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(poroHarmZ_old(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(poroHarmZInit(indexl:indexr,indexd:indexu,indexf:indexb))

        ! initialize K
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
        allocate(Kxx(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(Kyy(indexl:indexr,indexd:indexu,indexf:indexb))
        allocate(Kzz(indexl:indexr,indexd:indexu,indexf:indexb))
        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    Kxx(i,j,k) = KxxInit(xlower+i-1,ylower+j-1,zlower+k-1)
                    Kyy(i,j,k) = KyyInit(xlower+i-1,ylower+j-1,zlower+k-1)
                    Kzz(i,j,k) = KzzInit(xlower+i-1,ylower+j-1,zlower+k-1)
                end do
            end do
        end do

        ! initialize KxxHarm
        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays
        allocate(KxxHarm(indexl:indexr,indexd:indexu,indexf:indexb))

        ! initialize KyyHarm
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1
        indexf = 1
        indexb = localnlays
        allocate(KyyHarm(indexl:indexr,indexd:indexu,indexf:indexb))

        ! initialize KzzHarm
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays + 1
        allocate(KzzHarm(indexl:indexr,indexd:indexu,indexf:indexb))

        ! initialize av
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays
        allocate(av(indexl:indexr,indexd:indexu,indexf:indexb))
        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    av(i,j,k) = avInit(xlower+i-1,ylower+j-1,zlower+k-1)
                end do
            end do
        end do

        ! initialize p
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        indexb = localnlays
        allocate(p(indexl:indexr,indexd:indexu,indexf:indexb))

        ! initialize vx, vy, vz
        indexl = 1
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
        allocate(vx(indexl:indexr,indexd:indexu,indexf:indexb))
        vx(:,:,:) = 0.D0
        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    if(isDiriX0_p(ylower+j-1,zlower+k-1) == 0) then
                        vx(1,j,k) = vxBdryX0(ylower+j-1,zlower+k-1)
                    end if
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    if(isDiriX1_p(ylower+j-1,zlower+k-1) == 0) then
                        vx(localncols+1,j,k) = vxBdryX1(ylower+j-1,zlower+k-1)
                    end if
                end do
            end do
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
        indexd = 1
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
        allocate(vy(indexl:indexr,indexd:indexu,indexf:indexb))
        vy(:,:,:) = 0.D0
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    if(isDiriY0_p(xlower+i-1,zlower+k-1) == 0) then
                        vy(i,1,k) = vyBdryY0(xlower+i-1,zlower+k-1)
                    end if
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    if(isDiriY1_p(xlower+i-1,zlower+k-1) == 0) then
                        vy(i,localnrows+1,k) = vyBdryY1(xlower+i-1,zlower+k-1)
                    end if
                end do
            end do
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
        indexf = 1
        indexb = localnlays + 1
        allocate(vz(indexl:indexr,indexd:indexu,indexf:indexb))
        vz(:,:,:) = 0.D0
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    if(isDiriZ0_p(xlower+i-1,ylower+j-1) == 0) then
                        vz(i,j,1) = vzBdryZ0(xlower+i-1,ylower+j-1)
                    end if
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    if(isDiriZ1_p(xlower+i-1,ylower+j-1) == 0) then
                        vz(i,j,localnlays+1) = vzBdryZ1(xlower+i-1,ylower+j-1)
                    end if
                end do
            end do
        end if

        ! initialize Cf
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays
        allocate(Cf(indexl:indexr,indexd:indexu,indexf:indexb))
        Cf(indexl:indexr,indexd:indexu,indexf:indexb) = CfInit(xlower:xupper,ylower:yupper,zlower:zupper)

        ! initialize Tem
        allocate(Tem(indexl:indexr,indexd:indexu,indexf:indexb))
        Tem(indexl:indexr,indexd:indexu,indexf:indexb) = TemInit(xlower:xupper,ylower:yupper,zlower:zupper)

        ! compute the size of the equations
        if(pcol /= pncols) then
            xmomSize = localncols*localnrows*localnlays
        else
            xmomSize = (localncols+1)*localnrows*localnlays
        end if
        if(prow /= pnrows) then
            ymomSize = localncols*localnrows*localnlays
        else
            ymomSize = localncols*(localnrows+1)*localnlays
        end if
        if(play /= pnlays) then
            zmomSize = localncols*localnrows*localnlays
        else
            zmomSize = localncols*localnrows*(localnlays+1)
        end if
        massSize = localncols*localnrows*localnlays
        concenSize = localncols*localnrows*localnlays
        enerSize = localncols*localnrows*localnlays

        allocate(AxxEntryNum(xmomSize))
        allocate(AyyEntryNum(ymomSize))
        allocate(AzzEntryNum(zmomSize))
        allocate(ApEntryNum(massSize))
        allocate(AcfEntryNum(concenSize))
        allocate(AtemEntryNum(enerSize))

        call computeMatEntryNum()

        allocate(AxxEntryBase(xmomSize))
        allocate(AyyEntryBase(ymomSize))
        allocate(AzzEntryBase(zmomSize))
        allocate(ApEntryBase(massSize))
        allocate(AcfEntryBase(concenSize))
        allocate(AtemEntryBase(enerSize))

        AxxEntryBase(1) = 1
        do n = 2, xmomSize
            AxxEntryBase(n) = AxxEntryBase(n-1) + AxxEntryNum(n-1)
        end do

        AyyEntryBase(1) = 1
        do n = 2, ymomSize
            AyyEntryBase(n) = AyyEntryBase(n-1) + AyyEntryNum(n-1)
        end do

        AzzEntryBase(1) = 1
        do n = 2, zmomSize
            AzzEntryBase(n) = AzzEntryBase(n-1) + AzzEntryNum(n-1)
        end do

        ApEntryBase(1) = 1
        do n = 2, massSize
            ApEntryBase(n) = ApEntryBase(n-1) + ApEntryNum(n-1)
        end do

        AcfEntryBase(1) = 1
        do n = 2, concenSize
            AcfEntryBase(n) = AcfEntryBase(n-1) + AcfEntryNum(n-1)
        end do

        AtemEntryBase(1) = 1
        do n = 2, enerSize
            AtemEntryBase(n) = AtemEntryBase(n-1) + AtemEntryNum(n-1)
        end do

        ! compute matrix size
        AxxSize = 0
        do n = 1, xmomSize
            AxxSize = AxxSize + AxxEntryNum(n)
        end do
        AyySize = 0
        do n = 1, ymomSize
            AyySize = AyySize + AyyEntryNum(n)
        end do
        AzzSize = 0
        do n = 1, zmomSize
            AzzSize = AzzSize + AzzEntryNum(n)
        end do
        ApSize = 0
        do n = 1, massSize
            ApSize = ApSize + ApEntryNum(n)
        end do
        AcfSize = 0
        do n = 1, concenSize
            AcfSize = AcfSize + AcfEntryNum(n)
        end do
        AtemSize = 0
        do n = 1, enerSize
            AtemSize = AtemSize + AtemEntryNum(n)
        end do

        ! initialize matrix
        allocate(local_rhs_static_v(local_x_size_v))
        allocate(local_rhs_v(local_x_size_v))
        allocate(local_rhs_p(local_x_size_p))
        allocate(local_rhs_Cf(local_x_size_Cf))
        allocate(local_rhs_Tem(local_x_size_Tem))

        allocate(AxxCols(AxxSize))
        allocate(AxxRows(AxxSize))
        allocate(AxxStaticValues(AxxSize))
        allocate(AxxDynValues(AxxSize))
        allocate(AxxValues(AxxSize))
        allocate(AyyCols(AyySize))
        allocate(AyyRows(AyySize))
        allocate(AyyStaticValues(AyySize))
        allocate(AyyDynValues(AyySize))
        allocate(AyyValues(AyySize))
        allocate(AzzCols(AzzSize))
        allocate(AzzRows(AzzSize))
        allocate(AzzStaticValues(AzzSize))
        allocate(AzzDynValues(AzzSize))
        allocate(AzzValues(AzzSize))
        allocate(ApCols(ApSize))
        allocate(ApRows(ApSize))
        allocate(ApValues(ApSize))
        allocate(AcfCols(AcfSize))
        allocate(AcfRows(AcfSize))
        allocate(AcfValues(AcfSize))
        allocate(AtemCols(AtemSize))
        allocate(AtemRows(AtemSize))
        allocate(AtemValues(AtemSize))

        AxxCols(:) = 0
        AxxRows(:) = 0
        AxxStaticValues(:) = 0.D0
        AxxDynValues(:) = 0.D0
        AxxValues(:) = 0.D0
        AyyCols(:) = 0
        AyyRows(:) = 0
        AyyStaticValues(:) = 0.D0
        AyyDynValues(:) = 0.D0
        AyyValues(:) = 0.D0
        AzzCols(:) = 0
        AzzRows(:) = 0
        AzzStaticValues(:) = 0.D0
        AzzDynValues(:) = 0.D0
        AzzValues(:) = 0.D0
        ApCols(:) = 0
        ApRows(:) = 0
        ApValues(:) = 0.D0
        AcfCols(:) = 0
        AcfRows(:) = 0
        AcfValues(:) = 0.D0
        AtemCols(:) = 0
        AtemRows(:) = 0
        AtemValues(:) = 0.D0

        presDropInit = 0.D0
        isFindPresDropInit = .false.

        ! open the data output files
        if(myid == 0) then
            inquire(file = trim(adjustl(soludoc)), exist = alive)
            if(.not.alive) then
                call system("mkdir "//trim(adjustl(soludoc)))
            end if
            open(unit=40, file=trim(adjustl(soludoc))//'/his_poro_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_poro_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=41, file=trim(adjustl(soludoc))//'/his_Kxx_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_Kxx_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=42, file=trim(adjustl(soludoc))//'/his_av_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_av_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=43, file=trim(adjustl(soludoc))//'/his_p_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_p_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=44, file=trim(adjustl(soludoc))//'/his_Cf_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_Cf_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=45, file=trim(adjustl(soludoc))//'/his_Tem_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_Tem_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=46, file=trim(adjustl(soludoc))//'/his_q_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_q_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=47, file=trim(adjustl(soludoc))//'/his_lp_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_lp_avg.txt', ' error. ', ierr
                stop
            end if
        end if

#ifdef HYPRE

        allocate(initial_x_guess_v(local_x_size_v))
        allocate(initial_x_guess_p(local_x_size_p))
        allocate(initial_x_guess_Cf(local_x_size_Cf))
        allocate(initial_x_guess_Tem(local_x_size_Tem))

#endif

        if(mod(nt, NUMFRAME) /= 0) then
            if(myid == 0) then
                print *, 'The number of time steps must be divided by the number of frames.'
            end if
            stop
        end if

        deallocate(modelCase%xs)
        deallocate(modelCase%ys)
        deallocate(modelCase%zs)
        deallocate(modelCase%ts)
        deallocate(modelCase%src)
        deallocate(modelCase%poroInit)
        deallocate(modelCase%KxxInit)
        deallocate(modelCase%KyyInit)
        deallocate(modelCase%KzzInit)
        deallocate(modelCase%avInit)
        deallocate(modelCase%vxBdryX0)
        deallocate(modelCase%vxBdryX1)
        deallocate(modelCase%vxBdryY0)
        deallocate(modelCase%vxBdryY1)
        deallocate(modelCase%vxBdryZ0)
        deallocate(modelCase%vxBdryZ1)
        deallocate(modelCase%vyBdryX0)
        deallocate(modelCase%vyBdryX1)
        deallocate(modelCase%vyBdryY0)
        deallocate(modelCase%vyBdryY1)
        deallocate(modelCase%vyBdryZ0)
        deallocate(modelCase%vyBdryZ1)
        deallocate(modelCase%vzBdryX0)
        deallocate(modelCase%vzBdryX1)
        deallocate(modelCase%vzBdryY0)
        deallocate(modelCase%vzBdryY1)
        deallocate(modelCase%vzBdryZ0)
        deallocate(modelCase%vzBdryZ1)
        deallocate(modelCase%isDiriX0_p)
        deallocate(modelCase%isDiriX1_p)
        deallocate(modelCase%isDiriY0_p)
        deallocate(modelCase%isDiriY1_p)
        deallocate(modelCase%isDiriZ0_p)
        deallocate(modelCase%isDiriZ1_p)
        deallocate(modelCase%pBdryX0)
        deallocate(modelCase%pBdryX1)
        deallocate(modelCase%pBdryY0)
        deallocate(modelCase%pBdryY1)
        deallocate(modelCase%pBdryZ0)
        deallocate(modelCase%pBdryZ1)
        deallocate(modelCase%pInit)
        deallocate(modelCase%isDiriX0_Cf)
        deallocate(modelCase%isDiriX1_Cf)
        deallocate(modelCase%isDiriY0_Cf)
        deallocate(modelCase%isDiriY1_Cf)
        deallocate(modelCase%isDiriZ0_Cf)
        deallocate(modelCase%isDiriZ1_Cf)
        deallocate(modelCase%CfBdryX0)
        deallocate(modelCase%CfBdryX1)
        deallocate(modelCase%CfBdryY0)
        deallocate(modelCase%CfBdryY1)
        deallocate(modelCase%CfBdryZ0)
        deallocate(modelCase%CfBdryZ1)
        deallocate(modelCase%CfInit)
        deallocate(modelCase%isDiriX0_Tem)
        deallocate(modelCase%isDiriX1_Tem)
        deallocate(modelCase%isDiriY0_Tem)
        deallocate(modelCase%isDiriY1_Tem)
        deallocate(modelCase%isDiriZ0_Tem)
        deallocate(modelCase%isDiriZ1_Tem)
        deallocate(modelCase%TemBdryX0)
        deallocate(modelCase%TemBdryX1)
        deallocate(modelCase%TemBdryY0)
        deallocate(modelCase%TemBdryY1)
        deallocate(modelCase%TemBdryZ0)
        deallocate(modelCase%TemBdryZ1)
        deallocate(modelCase%TemInit)

    end subroutine initialize

    ! Generate the values on the right-hand side and the coefficients of the matrix A,
    ! and the subroutine will generate the values that will not change with the time iteration.
    subroutine genStaticPara(fi_kind)

        integer, intent(in) :: fi_kind

        integer :: findexl, findexr, findexd, findexu, findexf, findexb ! field index
        integer :: eindexl, eindexr, eindexd, eindexu, eindexf, eindexb ! equation index
        integer, dimension(:,:,:), pointer :: field
        logical :: isField
        real(kind=8), dimension(:,:,:), pointer :: rhs
        real(kind=8), dimension(:,:,:), pointer :: resiA, resitemp
        integer :: global_ind
        integer :: i, j, k, n, c

        if(fi_kind == 151) then
            if(pcol /= 1) then
                findexl = 0
            else
                findexl = 1
            end if
            findexr = localncols + 1
            if(prow /= 1) then
                findexd = 0
            else
                findexd = 1
            end if
            if(prow /= pnrows) then
                findexu = localnrows + 1
            else
                findexu = localnrows
            end if
            if(play /= 1) then
                findexf = 0
            else
                findexf = 1
            end if
            if(play /= pnlays) then
                findexb = localnlays + 1
            else
                findexb = localnlays
            end if
            eindexl = 1
            if(pcol /= pncols) then
                eindexr = localncols
            else
                eindexr = localncols + 1
            end if
            eindexd = 1
            eindexu = localnrows
            eindexf = 1
            eindexb = localnlays
        elseif(fi_kind == 253) then
            if(pcol /= 1) then
                findexl = 0
            else
                findexl = 1
            end if
            if(pcol /= pncols) then
                findexr = localncols + 1
            else
                findexr = localncols
            end if
            if(prow /= 1) then
                findexd = 0
            else
                findexd = 1
            end if
            findexu = localnrows + 1
            if(play /= 1) then
                findexf = 0
            else
                findexf = 1
            end if
            if(play /= pnlays) then
                findexb = localnlays + 1
            else
                findexb = localnlays
            end if
            eindexl = 1
            eindexr = localncols
            eindexd = 1
            if(prow /= pnrows) then
                eindexu = localnrows
            else
                eindexu = localnrows + 1
            end if
            eindexf = 1
            eindexb = localnlays
        elseif(fi_kind == 355) then
            if(pcol /= 1) then
                findexl = 0
            else
                findexl = 1
            end if
            if(pcol /= pncols) then
                findexr = localncols + 1
            else
                findexr = localncols
            end if
            if(prow /= 1) then
                findexd = 0
            else
                findexd = 1
            end if
            if(prow /= pnrows) then
                findexu = localnrows + 1
            else
                findexu = localnrows
            end if
            if(play /= 1) then
                findexf = 0
            else
                findexf = 1
            end if
            findexb = localnlays + 1
            eindexl = 1
            eindexr = localncols
            eindexd = 1
            eindexu = localnrows
            eindexf = 1
            if(play /= pnlays) then
                eindexb = localnlays
            else
                eindexb = localnlays + 1
            end if
        end if

        allocate(field(findexl:findexr,findexd:findexu,findexf:findexb))
        field(:,:,:) = 0
        allocate(rhs(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        if(fi_kind == 151) then
            call Resi_xmom_vx_b(field, rhs)
        elseif(fi_kind == 253) then
            call Resi_ymom_vy_b(field, rhs)
        elseif(fi_kind == 355) then
            call Resi_zmom_vz_b(field, rhs)
        end if
        allocate(resiA(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, field, isField)
                    if(isField) then
                        if(fi_kind == 151) then
                            call Resi_xmom_vx_b(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_xxlap(field, resitemp, 111, 151)
                        elseif(fi_kind == 253) then
                            call Resi_ymom_vy_b(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_yylap(field, resitemp, 211, 253)
                        elseif(fi_kind == 355) then
                            call Resi_zmom_vz_b(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_zzlap(field, resitemp, 311, 355)
                        end if
                    end if
                end do
            end do
        end do

        if(fi_kind == 151) then

            do k = 1, localnlays
                do j = 1, localnrows
                    if((pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                        call coordiToGlobalInd(myid, 111, 1, j, k, global_ind)
                        do n = 1, AxxSize
                            if(AxxRows(n)==global_ind) then
                                AxxValues(n) = AxxStaticValues(n)
                            end if
                        end do
                    end if
                    if((pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                        call coordiToGlobalInd(myid, 111, localncols+1, j, k, global_ind)
                        do n = 1, AxxSize
                            if(AxxRows(n)==global_ind) then
                                AxxValues(n) = AxxStaticValues(n)
                            end if
                        end do
                    end if
                end do
            end do

            c = 0
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        c = c + 1
                        local_rhs_static_v(c) = -rhs(i,j,k)
                        if(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)).or.((pcol==pncols).and. &!
                            (i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0))) then
                            local_rhs_v(c) = local_rhs_static_v(c)
                        end if
                    end do
                end do
            end do

        elseif(fi_kind == 253) then

            do k = 1, localnlays
                do i = 1, localncols
                    if((prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                        call coordiToGlobalInd(myid, 211, i, 1, k, global_ind)
                        do n = 1, AyySize
                            if(AyyRows(n)==global_ind) then
                                AyyValues(n) = AyyStaticValues(n)
                            end if
                        end do
                    end if
                    if((prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                        call coordiToGlobalInd(myid, 211, i, localnrows+1, k, global_ind)
                        do n = 1, AyySize
                            if(AyyRows(n)==global_ind) then
                                AyyValues(n) = AyyStaticValues(n)
                            end if
                        end do
                    end if
                end do
            end do

            if(pcol /= pncols) then
                c = localncols * localnrows * localnlays
            else
                c = (localncols+1) * localnrows * localnlays
            end if
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        c = c + 1
                        local_rhs_static_v(c) = -rhs(i,j,k)
                        if(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)).or. &!
                            ((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0))) then
                            local_rhs_v(c) = local_rhs_static_v(c)
                        end if
                    end do
                end do
            end do

        elseif(fi_kind == 355) then

            do j = 1, localnrows
                do i = 1, localncols
                    if((play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                        call coordiToGlobalInd(myid, 311, i, j, 1, global_ind)
                        do n = 1, AzzSize
                            if(AzzRows(n)==global_ind) then
                                AzzValues(n) = AzzStaticValues(n)
                            end if
                        end do
                    end if
                    if((play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                        call coordiToGlobalInd(myid, 311, i, j, localnlays+1, global_ind)
                        do n = 1, AzzSize
                            if(AzzRows(n)==global_ind) then
                                AzzValues(n) = AzzStaticValues(n)
                            end if
                        end do
                    end if
                end do
            end do

            if(pcol /= pncols) then
                c = localncols * localnrows * localnlays
            else
                c = (localncols+1) * localnrows * localnlays
            end if
            if(prow /= pnrows) then
                c = c + localncols * localnrows * localnlays
            else
                c = c + localncols * (localnrows+1) * localnlays
            end if
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        c = c + 1
                        local_rhs_static_v(c) = -rhs(i,j,k)
                        if(((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)).or. &!
                            ((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0))) then
                            local_rhs_v(c) = local_rhs_static_v(c)
                        end if
                    end do
                end do
            end do

        end if

        deallocate(field)
        deallocate(rhs)
        deallocate(resiA)
        deallocate(resitemp)

    end subroutine genStaticPara

    subroutine computedm()

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu, requestf, requestb
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, k, c, ierr

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    dm(i,j,k) = dmRef * exp(Eg/Rg*(1/TemRef_dm-1/Tem(i,j,k)))
                end do
            end do
        end do

        call communicateP(1)

    end subroutine computedm

    subroutine computeks()

        integer :: i, j, k

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    ks(i,j,k) = ksRef * exp(Eg/Rg*(1/TemRef_ks-1/Tem(i,j,k)))

                end do
            end do
        end do

    end subroutine computeks

    subroutine computekc()

        real(kind=8) :: Sc, vmodulus, radius, Rep
        integer :: i, j, k

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    vmodulus = dsqrt((vx(i+1,j,k)-vx(i,j,k))**2.D0+(vy(i,j+1,k)-vy(i,j,k))**2.D0+ &!
                        (vz(i,j,k+1)-vz(i,j,k))**2.D0)

                    radius = radiusInit*(1.D0-poroInit(xlower+i-1,ylower+j-1,zlower+k-1))/ &!
                        (poroInit(xlower+i-1,ylower+j-1,zlower+k-1)*(1.D0-poro(i,j,k)))*poro(i,j,k)

                    Rep = 2.D0*vmodulus*radius/(visc/rhof)

                    Sc = visc/rhof*dm(i,j,k)

                    kc(i,j,k) = (ShInfinity+7.D-1*Rep**(1.D0/2.D0)*Sc**(1.D0/3.D0))*dm(i,j,k)/2.D0/radius

                end do
            end do
        end do

    end subroutine computekc

    subroutine computePoro()

        real(kind=8) :: coe
        integer :: i, j, k

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    coe = avInit(xlower+i-1,ylower+j-1,zlower+k-1)*al*Cf(i,j,k)*kc(i,j,k)*ks(i,j,k)*(ts(t)-ts(t-1))/ &!
                        (rhos*(kc(i,j,k)+ks(i,j,k))*(1.D0-poroInit(xlower+i-1,ylower+j-1,zlower+k-1)))
                    poro_old(i,j,k) = poro(i,j,k)
                    poro(i,j,k) = (coe+poro_old(i,j,k))/(1+coe)
                end do
            end do
        end do

        call communicateP(2)

    end subroutine computePoro

    subroutine computePoroHarm(var)

        integer, intent(in) :: var

        real(kind=8), dimension(:,:,:), pointer :: HarmX, HarmY, HarmZ
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: i, j, k

        if(var == 1) then
            HarmX => poroHarmX_old
            HarmY => poroHarmY_old
            HarmZ => poroHarmZ_old
        elseif(var == 2) then
            HarmX => poroHarmX
            HarmY => poroHarmY
            HarmZ => poroHarmZ
        else
            print *, 'The var value in computePoroHarm is wrong!'
            stop
        end if

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

        if(pcol == 1) then
            HarmX(1,indexd:indexu,indexf:indexb) = poro(1,indexd:indexu,indexf:indexb)
            indexl = indexl + 1
        end if

        if(pcol == pncols) then
            HarmX(localncols+1,indexd:indexu,indexf:indexb) = poro(localncols,indexd:indexu,indexf:indexb)
            indexr = indexr - 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    HarmX(i,j,k) = (hx(i-1)+hx(i)) / (hx(i-1)/poro(i-1,j,k)+hx(i)/poro(i,j,k))
                end do
            end do
        end do

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

        if(prow == 1) then
            HarmY(indexl:indexr,1,indexf:indexb) = poro(indexl:indexr,1,indexf:indexb)
            indexd = indexd + 1
        end if

        if(prow == pnrows) then
            HarmY(indexl:indexr,localnrows+1,indexf:indexb) = poro(indexl:indexr,localnrows,indexf:indexb)
            indexu = indexu - 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    HarmY(i,j,k) = (hy(j-1)+hy(j)) / (hy(j-1)/poro(i,j-1,k)+hy(j)/poro(i,j,k))
                end do
            end do
        end do

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

        if(play == 1) then
            HarmZ(indexl:indexr,indexd:indexu,1) = poro(indexl:indexr,indexd:indexu,1)
            indexf = indexf + 1
        end if

        if(play == pnlays) then
            HarmZ(indexl:indexr,indexd:indexu,localnlays+1) = poro(indexl:indexr,indexd:indexu,localnlays)
            indexb = indexb - 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    HarmZ(i,j,k) = (hz(k-1)+hz(k)) / (hz(k-1)/poro(i,j,k-1)+hz(k)/poro(i,j,k))
                end do
            end do
        end do

        if((t == 2).and.(var == 1)) then
            poroHarmXInit = HarmX
            poroHarmYInit = HarmY
            poroHarmZInit = HarmZ
        end if

    end subroutine computePoroHarm

    subroutine computeK()

        integer :: i, j, k

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    Kxx(i,j,k) = poro(i,j,k)/poroInit(xlower+i-1,ylower+j-1,zlower+k-1)*(poro(i,j,k)* &!
                        (1.D0-poroInit(xlower+i-1,ylower+j-1,zlower+k-1))/poroInit(xlower+i-1,ylower+j-1,zlower+k-1) &!
                        /(1.D0-poro(i,j,k)))**2 * KxxInit(xlower+i-1,ylower+j-1,zlower+k-1)
                    Kyy(i,j,k) = poro(i,j,k)/poroInit(xlower+i-1,ylower+j-1,zlower+k-1)*(poro(i,j,k)* &!
                        (1.D0-poroInit(xlower+i-1,ylower+j-1,zlower+k-1))/poroInit(xlower+i-1,ylower+j-1,zlower+k-1) &!
                        /(1.D0-poro(i,j,k)))**2 * KyyInit(xlower+i-1,ylower+j-1,zlower+k-1)
                    Kzz(i,j,k) = poro(i,j,k)/poroInit(xlower+i-1,ylower+j-1,zlower+k-1)*(poro(i,j,k)* &!
                        (1.D0-poroInit(xlower+i-1,ylower+j-1,zlower+k-1))/poroInit(xlower+i-1,ylower+j-1,zlower+k-1) &!
                        /(1.D0-poro(i,j,k)))**2 * KzzInit(xlower+i-1,ylower+j-1,zlower+k-1)
                end do
            end do
        end do

        call communicateP(3)

    end subroutine computeK

    subroutine computeKHarm()

        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: i, j, k

        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays

        if(pcol == 1) then
            KxxHarm(1,indexd:indexu,indexf:indexb) = Kxx(1,indexd:indexu,indexf:indexb)
            indexl = indexl + 1
        end if

        if(pcol == pncols) then
            KxxHarm(localncols+1,indexd:indexu,indexf:indexb) = Kxx(localncols,indexd:indexu,indexf:indexb)
            indexr = indexr - 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    KxxHarm(i,j,k) = (hx(i-1)+hx(i)) / (hx(i-1)/Kxx(i-1,j,k)+hx(i)/Kxx(i,j,k))
                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1
        indexf = 1
        indexb = localnlays

        if(prow == 1) then
            KyyHarm(indexl:indexr,1,indexf:indexb) = Kyy(indexl:indexr,1,indexf:indexb)
            indexd = indexd + 1
        end if

        if(prow == pnrows) then
            KyyHarm(indexl:indexr,localnrows+1,indexf:indexb) = Kyy(indexl:indexr,localnrows,indexf:indexb)
            indexu = indexu - 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    KyyHarm(i,j,k) = (hy(j-1)+hy(j)) / (hy(j-1)/Kyy(i,j-1,k)+hy(j)/Kyy(i,j,k))
                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays + 1

        if(play == 1) then
            KzzHarm(indexl:indexr,indexd:indexu,1) = Kzz(indexl:indexr,indexd:indexu,1)
            indexf = indexf + 1
        end if

        if(play == pnlays) then
            KzzHarm(indexl:indexr,indexd:indexu,localnlays+1) = Kzz(indexl:indexr,indexd:indexu,localnlays)
            indexb = indexb - 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    KzzHarm(i,j,k) = (hz(k-1)+hz(k)) / (hz(k-1)/Kzz(i,j,k-1)+hz(k)/Kzz(i,j,k))
                end do
            end do
        end do

    end subroutine computeKHarm

    subroutine computeav()

        integer :: i, j, k

        ! suppose Kxx=Kyy=Kzz. if you find the anisotropic permeability equation to compute av, you can change the
        ! computation here.
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    av(i,j,k) = avInit(xlower+i-1,ylower+j-1,zlower+k-1)*(1.D0-poro(i,j,k))/ &!
                        (1.D0-poroInit(xlower+i-1,ylower+j-1,zlower+k-1))

                end do
            end do
        end do

    end subroutine computeav

    ! Generate the values on the right-hand side and the coefficients of the matrix A,
    ! and the subroutine will generate the values that will change with the time iteration.
    subroutine genDynPara(fi_kind)

        integer, intent(in) :: fi_kind

        integer :: findexl, findexr, findexd, findexu, findexf, findexb
        integer :: eindexl, eindexr, eindexd, eindexu, eindexf, eindexb
        integer, dimension(:,:,:), pointer :: field
        logical :: isField
        real(kind=8), dimension(:,:,:), pointer :: rhs
        real(kind=8), dimension(:,:,:), pointer :: resiA, resitemp
        integer :: AxxBeInd, AyyBeInd, AzzBeInd, global_ind
        integer :: i, j, k, n, c

        if(fi_kind == 152) then
            if(pcol /= 1) then
                findexl = 0
            else
                findexl = 1
            end if
            findexr = localncols + 1
            if(prow /= 1) then
                findexd = 0
            else
                findexd = 1
            end if
            if(prow /= pnrows) then
                findexu = localnrows + 1
            else
                findexu = localnrows
            end if
            if(play /= 1) then
                findexf = 0
            else
                findexf = 1
            end if
            if(play /= pnlays) then
                findexb = localnlays + 1
            else
                findexb = localnlays
            end if
            eindexl = 1
            if(pcol /= pncols) then
                eindexr = localncols
            else
                eindexr = localncols + 1
            end if
            eindexd = 1
            eindexu = localnrows
            eindexf = 1
            eindexb = localnlays
        elseif(fi_kind == 254) then
            if(pcol /= 1) then
                findexl = 0
            else
                findexl = 1
            end if
            if(pcol /= pncols) then
                findexr = localncols + 1
            else
                findexr = localncols
            end if
            if(prow /= 1) then
                findexd = 0
            else
                findexd = 1
            end if
            findexu = localnrows + 1
            if(play /= 1) then
                findexf = 0
            else
                findexf = 1
            end if
            if(play /= pnlays) then
                findexb = localnlays + 1
            else
                findexb = localnlays
            end if
            eindexl = 1
            eindexr = localncols
            eindexd = 1
            if(prow /= pnrows) then
                eindexu = localnrows
            else
                eindexu = localnrows + 1
            end if
            eindexf = 1
            eindexb = localnlays
        elseif(fi_kind == 356) then
            if(pcol /= 1) then
                findexl = 0
            else
                findexl = 1
            end if
            if(pcol /= pncols) then
                findexr = localncols + 1
            else
                findexr = localncols
            end if
            if(prow /= 1) then
                findexd = 0
            else
                findexd = 1
            end if
            if(prow /= pnrows) then
                findexu = localnrows + 1
            else
                findexu = localnrows
            end if
            if(play /= 1) then
                findexf = 0
            else
                findexf = 1
            end if
            findexb = localnlays + 1
            eindexl = 1
            eindexr = localncols
            eindexd = 1
            eindexu = localnrows
            eindexf = 1
            if(play /= pnlays) then
                eindexb = localnlays
            else
                eindexb = localnlays + 1
            end if
        else
            findexl = 0
            findexr = localncols + 1
            findexd = 0
            findexu = localnrows + 1
            findexf = 0
            findexb = localnlays + 1
            eindexl = 1
            eindexr = localncols
            eindexd = 1
            eindexu = localnrows
            eindexf = 1
            eindexb = localnlays
        end if
        allocate(field(findexl:findexr,findexd:findexu,findexf:findexb))
        field(:,:,:) = 0
        allocate(rhs(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        rhs(:,:,:) = 0.D0

        if(fi_kind == 152) then
            call Resi_xmom_vx(field, rhs)
        elseif(fi_kind == 254) then
            call Resi_ymom_vy(field, rhs)
        elseif(fi_kind == 356) then
            call Resi_zmom_vz(field, rhs)
        elseif(fi_kind == 457) then
            call Resi_mass_p(field, rhs)
        elseif(fi_kind == 458) then
            call Resi_concen_cf(field, rhs)
        elseif(fi_kind == 459) then
            call Resi_ener_tem(field, rhs)
        end if

        allocate(resiA(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, field, isField)
                    if(isField) then
                        if(fi_kind == 152) then
                            call Resi_xmom_vx(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_xxlap(field, resitemp, 111, 152)
                        elseif(fi_kind == 254) then
                            call Resi_ymom_vy(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_yylap(field, resitemp, 211, 254)
                        elseif(fi_kind == 356) then
                            call Resi_zmom_vz(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_zzlap(field, resitemp, 311, 356)
                        elseif(fi_kind == 457) then
                            call Resi_mass_p(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_pplap7(field, resitemp, 412, 457)
                        elseif(fi_kind == 458) then
                            call Resi_concen_cf(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_pplap19(field, resitemp, 413, 458)
                        elseif(fi_kind == 459) then
                            call Resi_ener_tem(field, resiA)
                            resitemp = resiA - rhs
                            call dctz_pplap7(field, resitemp, 414, 459)
                        end if
                    end if
                end do
            end do
        end do

        if(fi_kind == 152) then

            AxxBeInd = 1
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        if(.not.(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)).or.((pcol==pncols).and. &!
                            (i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)))) then
                            call coordiToGlobalInd(myid, 111, i, j, k, global_ind)
                            do n = AxxBeInd, AxxSize
                                if(AxxRows(n)==global_ind) then
                                    AxxValues(n) = AxxStaticValues(n) + AxxDynValues(n)
                                else
                                    AxxBeInd = n
                                    exit
                                end if
                            end do
                        elseif(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)).or.((pcol==pncols).and. &!
                            (i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0))) then
                            AxxBeInd = AxxBeInd + 1
                        end if
                    end do
                end do
            end do

            c = 0
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        c = c + 1
                        if(.not.(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)).or.((pcol==pncols).and. &!
                            (i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)))) then
                            local_rhs_v(c) = local_rhs_static_v(c) - rhs(i,j,k)
                        end if
                    end do
                end do
            end do

        elseif(fi_kind == 254) then

            AyyBeInd = 1
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        if(.not.(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)).or.((prow==pnrows).and. &!
                            (j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)))) then
                            call coordiToGlobalInd(myid, 211, i, j, k, global_ind)
                            do n = AyyBeInd, AyySize
                                if(AyyRows(n)==global_ind) then
                                    AyyValues(n) = AyyStaticValues(n) + AyyDynValues(n)
                                else
                                    AyyBeInd = n
                                    exit
                                end if
                            end do
                        elseif(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)).or.((prow==pnrows).and. &!
                            (j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0))) then
                            AyyBeInd = AyyBeInd + 1
                        end if
                    end do
                end do
            end do

            if(pcol /= pncols) then
                c = localncols * localnrows * localnlays
            else
                c = (localncols+1) * localnrows * localnlays
            end if
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        c = c + 1
                        if(.not.(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)).or.((prow==pnrows).and. &!
                            (j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)))) then
                            local_rhs_v(c) = local_rhs_static_v(c) - rhs(i,j,k)
                        end if
                    end do
                end do
            end do

        elseif(fi_kind == 356) then

            AzzBeInd = 1
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        if(.not.(((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)).or.((play==pnlays).and. &!
                            (k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)))) then
                            call coordiToGlobalInd(myid, 311, i, j, k, global_ind)
                            do n = AzzBeInd, AzzSize
                                if(AzzRows(n)==global_ind) then
                                    AzzValues(n) = AzzStaticValues(n) + AzzDynValues(n)
                                else
                                    AzzBeInd = n
                                    exit
                                end if
                            end do
                        elseif(((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)).or.((play==pnlays).and. &!
                            (k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0))) then
                            AzzBeInd = AzzBeInd + 1
                        end if
                    end do
                end do
            end do

            if(pcol /= pncols) then
                c = localncols * localnrows * localnlays
            else
                c = (localncols+1) * localnrows * localnlays
            end if
            if(prow /= pnrows) then
                c = c + localncols * localnrows * localnlays
            else
                c = c + localncols * (localnrows+1) * localnlays
            end if
            do k = eindexf, eindexb
                do j = eindexd, eindexu
                    do i = eindexl, eindexr
                        c = c + 1
                        if(.not.(((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)).or.((play==pnlays).and. &!
                            (k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)))) then
                            local_rhs_v(c) = local_rhs_static_v(c) - rhs(i,j,k)
                        end if
                    end do
                end do
            end do

        end if

        c = 0
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    if(fi_kind == 457) then
                        local_rhs_p(c) = -rhs(i,j,k)
                    elseif(fi_kind == 458) then
                        local_rhs_Cf(c) = -rhs(i,j,k)
                    elseif(fi_kind == 459) then
                        local_rhs_Tem(c) = -rhs(i,j,k)
                    end if
                end do
            end do
        end do

        deallocate(field)
        deallocate(rhs)
        deallocate(resiA)
        deallocate(resitemp)

    end subroutine genDynPara

    subroutine solve(mat_kind)

        integer, intent(in) :: mat_kind

        integer :: AEntryNum_standard, local_x_size
        real(kind=8), dimension(:), allocatable :: local_rhs
        integer :: ASize, ilower, iupper
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer, dimension(:), allocatable :: ACols, ARows
        real(kind=8), dimension(:), allocatable :: AValues
        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), pointer :: local_x
        integer :: ABeInd, nVelx, nVely
        real(kind=8) :: solvertimestart, solvertimefinish
        integer :: i, j, k, l, n, c, ierr

#ifdef LAPACK

        ! the variables used in LAPACK
        real(kind=8), dimension(:,:), pointer :: A_lapack
        real(kind=8), dimension(:), pointer :: b_lapack
        integer :: LAPACKINFO
        integer, dimension(:), pointer :: IPIV

#elif defined(UMFPACK)

        ! the variables used in UMFPACK
        integer, dimension(:), allocatable :: Ap, Ai
        real(kind=8), dimension(:), allocatable :: Ax
        integer(kind=8) :: symbolic, numeric
        real(kind=8) :: control(20), umfinfo(90)

#elif defined(MUMPS)

        ! the variables used in MUMPS
        type(DMUMPS_STRUC) mumps_par
        integer(kind=8) :: mumps_NNZ_loc
        integer, dimension(:), allocatable ::  mumps_IRN_loc, mumps_JCN_loc
        real(kind=8), dimension(:), allocatable :: mumps_A_loc
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data

#elif defined(HYPRE)

        ! the parameters and variables used in HYPRE
        integer, parameter :: HYPRE_PARCSR = 5555
        integer(kind=8) :: A
        integer(kind=8) :: b
        integer(kind=8) :: x
        integer(kind=8) :: parcsr_A
        integer(kind=8) :: par_b
        integer(kind=8) :: par_x
        integer(kind=8) :: precond, solver
        integer :: jlower, jupper
        integer, dimension(:), pointer :: rows
        real(kind=8), dimension(:), allocatable :: initial_x_guess

#endif

        if(mat_kind == 11) then
            AEntryNum_standard = AxxEntryNum_standard
            local_x_size = local_x_size_v
            allocate(local_rhs(1:local_x_size))
            local_rhs(1:local_x_size) = local_rhs_v(1:local_x_size)
            ilower = ilower_v
            iupper = iupper_v
#ifdef HYPRE
            if(play /= 1) then
                call coordiToGlobalInd(myid-pncols*pnrows, 111, 1, 1, localnlays, jlower)
            elseif(prow /= 1) then
                call coordiToGlobalInd(myid-pncols, 111, 1, localnrows, 1, jlower)
            elseif(pcol /= 1) then
                call coordiToGlobalInd(myid-1, 111, localncols, 1, 1, jlower)
            else
                jlower = ilower
            end if
            if(play /= pnlays) then
                call coordiToGlobalInd(myid+pncols*pnrows, 311, localncols, localnrows, 1, jupper)
            elseif(prow /= pnrows) then
                call coordiToGlobalInd(myid+pncols, 311, localncols, 1, localnlays+1, jupper)
            elseif(pcol /= pncols) then
                call coordiToGlobalInd(myid+1, 311, 1, localnrows, localnlays+1, jupper)
            else
                jupper = iupper
            end if
#endif
        elseif(mat_kind == 12) then
            AEntryNum_standard = ApEntryNum_standard
            local_x_size = local_x_size_p
            allocate(local_rhs(1:local_x_size))
            local_rhs(1:local_x_size) = local_rhs_p(1:local_x_size)
            ASize = ApSize
            allocate(ACols(ASize))
            allocate(ARows(ASize))
            allocate(AValues(ASize))
            ACols(1:ASize) = ApCols(1:ASize)
            ARows(1:ASize) = ApRows(1:ASize)
            AValues(1:ASize) = ApValues(1:ASize)
            ilower = ilower_p
            iupper = iupper_p
#ifdef HYPRE
            if(play /= 1) then
                call coordiToGlobalInd(myid-pncols*pnrows, 412, 1, 1, localnlays, jlower)
            elseif(prow /= 1) then
                call coordiToGlobalInd(myid-pncols, 412, 1, localnrows, 1, jlower)
            elseif(pcol /= 1) then
                call coordiToGlobalInd(myid-1, 412, localncols, 1, 1, jlower)
            else
                jlower = ilower
            end if
            if(play /= pnlays) then
                call coordiToGlobalInd(myid+pncols*pnrows, 412, localncols, localnrows, 1, jupper)
            elseif(prow /= pnrows) then
                call coordiToGlobalInd(myid+pncols, 412, localncols, 1, localnlays, jupper)
            elseif(pcol /= pncols) then
                call coordiToGlobalInd(myid+1, 412, 1, localnrows, localnlays, jupper)
            else
                jupper = iupper
            end if
#endif
        elseif(mat_kind == 13) then
            AEntryNum_standard = AcfEntryNum_standard
            local_x_size = local_x_size_Cf
            allocate(local_rhs(1:local_x_size))
            local_rhs(1:local_x_size) = local_rhs_Cf(1:local_x_size)
            ASize = AcfSize
            allocate(ACols(ASize))
            allocate(ARows(ASize))
            allocate(AValues(ASize))
            ACols(1:ASize) = AcfCols(1:ASize)
            ARows(1:ASize) = AcfRows(1:ASize)
            AValues(1:ASize) = AcfValues(1:ASize)
            ilower = ilower_Cf
            iupper = iupper_Cf
#ifdef HYPRE
            if((prow /= 1).and.(play /= 1)) then
                call coordiToGlobalInd(myid-pncols-pncols*pnrows, 413, 1, localnrows, localnlays, jlower)
            elseif((pcol /= 1).and.(play /= 1)) then
                call coordiToGlobalInd(myid-1-pncols*pnrows, 413, localncols, 1, localnlays, jlower)
            elseif((pcol /= 1).and.(prow /= 1)) then
                call coordiToGlobalInd(myid-1-pncols, 413, localncols, localnrows, 1, jlower)
            elseif(play /= 1) then
                call coordiToGlobalInd(myid-pncols*pnrows, 413, 1, 1, localnlays, jlower)
            elseif(prow /= 1) then
                call coordiToGlobalInd(myid-pncols, 413, 1, localnrows, 1, jlower)
            elseif(pcol /= 1) then
                call coordiToGlobalInd(myid-1, 413, localncols, 1, 1, jlower)
            else
                jlower = ilower
            end if
            if((prow /= pnrows).and.(play /= pnlays)) then
                call coordiToGlobalInd(myid+pncols+pncols*pnrows, 413, localncols, 1, 1, jupper)
            elseif((pcol /= pncols).and.(play /= pnlays)) then
                call coordiToGlobalInd(myid+1+pncols*pnrows, 413, 1, localnrows, 1, jupper)
            elseif((pcol /= pncols).and.(prow /= pnrows)) then
                call coordiToGlobalInd(myid+1+pncols, 413, 1, 1, localnlays, jupper)
            elseif(play /= pnlays) then
                call coordiToGlobalInd(myid+pncols*pnrows, 413, localncols, localnrows, 1, jupper)
            elseif(prow /= pnrows) then
                call coordiToGlobalInd(myid+pncols, 413, localncols, 1, localnlays, jupper)
            elseif(pcol /= pncols) then
                call coordiToGlobalInd(myid+1, 413, 1, localnrows, localnlays, jupper)
            else
                jupper = iupper
            end if
#endif
        elseif(mat_kind == 14) then
            AEntryNum_standard = AtemEntryNum_standard
            local_x_size = local_x_size_Tem
            allocate(local_rhs(1:local_x_size))
            local_rhs(1:local_x_size) = local_rhs_Tem(1:local_x_size)
            ASize = AtemSize
            allocate(ACols(ASize))
            allocate(ARows(ASize))
            allocate(AValues(ASize))
            ACols(1:ASize) = AtemCols(1:ASize)
            ARows(1:ASize) = AtemRows(1:ASize)
            AValues(1:ASize) = AtemValues(1:ASize)
            ilower = ilower_Tem
            iupper = iupper_Tem
#ifdef HYPRE
            if(play /= 1) then
                call coordiToGlobalInd(myid-pncols*pnrows, 414, 1, 1, localnlays, jlower)
            elseif(prow /= 1) then
                call coordiToGlobalInd(myid-pncols, 414, 1, localnrows, 1, jlower)
            elseif(pcol /= 1) then
                call coordiToGlobalInd(myid-1, 414, localncols, 1, 1, jlower)
            else
                jlower = ilower
            end if
            if(play /= pnlays) then
                call coordiToGlobalInd(myid+pncols*pnrows, 414, localncols, localnrows, 1, jupper)
            elseif(prow /= pnrows) then
                call coordiToGlobalInd(myid+pncols, 414, localncols, 1, localnlays, jupper)
            elseif(pcol /= pncols) then
                call coordiToGlobalInd(myid+1, 414, 1, localnrows, localnlays, jupper)
            else
                jupper = iupper
            end if
#endif
        end if

        allocate(cols(AEntryNum_standard))
        allocate(values(AEntryNum_standard))
        allocate(local_x(local_x_size))

#ifdef LAPACK

        allocate(A_lapack(local_x_size,local_x_size))
        allocate(b_lapack(local_x_size))
        allocate(IPIV(local_x_size))

        A_lapack(:,:) = 0.D0
        b_lapack(:) = local_rhs(:)
        IPIV(:) = 0.D0

#elif defined(UMFPACK)

        allocate(Ap(1:local_x_size+1))
        allocate(Ai(1:AEntryNum_standard*local_x_size))
        allocate(Ax(1:AEntryNum_standard*local_x_size))

        Ap(1) = 0
        Ai(:) = 0
        Ax(:) = 0.D0

#elif defined(MUMPS)

        mumps_par%COMM = MPI_COMM_WORLD
        mumps_par%SYM = 0
        mumps_par%PAR = 1
        mumps_par%JOB = -1
        call DMUMPS(mumps_par)

        mumps_par%ICNTL(4) = 0
        mumps_par%ICNTL(7) = 15
        mumps_par%ICNTL(18) = 3
        if(mumps_par%MYID == 0) then
            if(mat_kind == 11) then
                mumps_par%N = (nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1)
            elseif(mat_kind == 12) then
                mumps_par%N = nx*ny*nz
            elseif(mat_kind == 13) then
                mumps_par%N = nx*ny*nz
            elseif(mat_kind == 14) then
                mumps_par%N = nx*ny*nz
            end if
        end if
        allocate(mumps_par%RHS(mumps_par%N))

        allocate(mumps_IRN_loc(local_x_size*AEntryNum_standard))
        allocate(mumps_JCN_loc(local_x_size*AEntryNum_standard))
        allocate(mumps_A_loc(local_x_size*AEntryNum_standard))

        mumps_NNZ_loc = 0
        mumps_IRN_loc(:) = 0
        mumps_JCN_loc(:) = 0
        mumps_A_loc(:) = 0.D0

        if(mumps_par%MYID /= 0) then
            call MPI_IBSEND(local_rhs, local_x_size, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(mumps_par%MYID == 0) then
            mumps_par%RHS(1:local_x_size) = local_rhs(1:local_x_size)
            c = local_x_size + 1
            do i = 1, nProcs-1
                if(mat_kind == 11) then
                    allocate(slave_data(slave_v_data_size(i)))
                    call MPI_RECV(slave_data, slave_v_data_size(i), MPI_DOUBLE_PRECISION, &!
                        i, i, MPI_COMM_WORLD, status, ierr)
                    mumps_par%RHS(c:c+slave_v_data_size(i)-1) = slave_data(:)
                    c = c + slave_v_data_size(i)
                    deallocate(slave_data)
                else
                    allocate(slave_data(local_x_size))
                    call MPI_RECV(slave_data, local_x_size, MPI_DOUBLE_PRECISION, &!
                        i, i, MPI_COMM_WORLD, status, ierr)
                    mumps_par%RHS(c:c+local_x_size-1) = slave_data(:)
                    c = c + local_x_size
                    deallocate(slave_data)
                end if
            end do
        end if

        if(mumps_par%MYID /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

#elif defined(HYPRE)

        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, A, ierr)
        call HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR, ierr)
        call HYPRE_IJMatrixInitialize(A, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, b, ierr)
        call HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(b, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, x, ierr)
        call HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(x, ierr)

        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver, ierr)
        call HYPRE_ParCSRGMRESSetMaxIter(solver, 50000, ierr)
        !call HYPRE_ParCSRGMRESSetPrintLevel(solver, 2, ierr)
        !call HYPRE_ParCSRGMRESSetLogging(solver, 1, ierr)

        call HYPRE_ParaSailsCreate(MPI_COMM_WORLD, precond,ierr)
        call HYPRE_ParaSailsSetParams(precond, -9D-1, 2, ierr)
        call HYPRE_ParaSailsSetFilter(precond, -9D-1, ierr)
        ! Because the matrix A is nonsymmetric and indefinite, you must choose the parameter as 0.
        call HYPRE_ParaSailsSetSym(precond, 0)
        !call HYPRE_ParaSailsSetLogging(precond, 1, ierr)

        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver, 1, precond, ierr)

        allocate(rows(local_x_size))
        do n = 1, local_x_size
            rows(n) = ilower + n - 1
        end do

        allocate(initial_x_guess(local_x_size))
        if(mat_kind == 11) then
            if(t == 2) then
                c = 0
                indexl = 1
                if(pcol /= pncols) then
                    indexr = localncols
                else
                    indexr = localncols + 1
                end if
                indexd = 1
                indexu = localnrows
                indexf = 1
                indexb = localnlays
                do k = indexf, indexb
                    do j = indexd, indexu
                        do i = indexl, indexr
                            c = c + 1
                            initial_x_guess(c) = vx(i,j,k)
                        end do
                    end do
                end do
                indexl = 1
                indexr = localncols
                indexd = 1
                if(prow /= pnrows) then
                    indexu = localnrows
                else
                    indexu = localnrows + 1
                end if
                indexf = 1
                indexb = localnlays
                do k = indexf, indexb
                    do j = indexd, indexu
                        do i = indexl, indexr
                            c = c + 1
                            initial_x_guess(c) = vy(i,j,k)
                        end do
                    end do
                end do
                indexl = 1
                indexr = localncols
                indexd = 1
                indexu = localnrows
                indexf = 1
                if(play /= pnlays) then
                    indexb = localnlays
                else
                    indexb = localnlays + 1
                end if
                do k = indexf, indexb
                    do j = indexd, indexu
                        do i = indexl, indexr
                            c = c + 1
                            initial_x_guess(c) = vz(i,j,k)
                        end do
                    end do
                end do
            else
                initial_x_guess(1:local_x_size) = initial_x_guess_v(1:local_x_size)
            end if
        elseif(mat_kind == 12) then
            if(t == 2) then
                c = 0
                do k = 1, localnlays
                    do j = 1, localnrows
                        do i = 1, localncols
                            c = c + 1
                            initial_x_guess(c) = pInit(xlower+i-1,ylower+j-1,zlower+k-1)
                        end do
                    end do
                end do
            else
                initial_x_guess(1:local_x_size) = initial_x_guess_p(1:local_x_size)
            end if
        elseif(mat_kind == 13) then
            if(t == 2) then
                c = 0
                do k = 1, localnlays
                    do j = 1, localnrows
                        do i = 1, localncols
                            c = c + 1
                            initial_x_guess(c) = CfInit(xlower+i-1,ylower+j-1,zlower+k-1)
                        end do
                    end do
                end do
            else
                initial_x_guess(1:local_x_size) = initial_x_guess_Cf(1:local_x_size)
            end if
        elseif(mat_kind == 14) then
            if(t == 2) then
                c = 0
                do k = 1, localnlays
                    do j = 1, localnrows
                        do i = 1, localncols
                            c = c + 1
                            initial_x_guess(c) = TemInit(xlower+i-1,ylower+j-1,zlower+k-1)
                        end do
                    end do
                end do
            else
                initial_x_guess(1:local_x_size) = initial_x_guess_Tem(1:local_x_size)
            end if
        end if

        call HYPRE_IJVectorSetValues(b, local_x_size, rows, local_rhs, ierr)
        call HYPRE_IJVectorSetValues(x, local_x_size, rows, initial_x_guess, ierr)

        call HYPRE_IJVectorAssemble(b, ierr)
        call HYPRE_IJVectorAssemble(x, ierr)

        call HYPRE_IJVectorGetObject(b, par_b, ierr)
        call HYPRE_IJVectorGetObject(x, par_x, ierr)

        deallocate(initial_x_guess)

#endif

        if(mat_kind == 11) then
            if(pcol /= pncols) then
                nVelx = localncols*localnrows*localnlays
            else
                nVelx = (localncols+1)*localnrows*localnlays
            end if
            if(prow /= pnrows) then
                nVely = localncols*localnrows*localnlays
            else
                nVely = localncols*(localnrows+1)*localnlays
            end if
        end if

        ABeInd = 1
        do n = ilower, iupper

            cols(:) = 0
            values(:) = 0.D0
            c = 0

            if(mat_kind == 11) then
                ! the line is in the x-velocity part
                if(n <= ilower+nVelx-1) then
                    do l = ABeInd, AxxSize
                        if(AxxRows(l) == n) then
                            c = c + 1
                            cols(c) = AxxCols(l)
                            values(c) = AxxValues(l)
                        else
                            ABeInd = l
                            exit
                        end if
                    end do
                ! the line is in the y-velocity part
                elseif((n >= ilower+nVelx).and.(n <= ilower+nVelx+nVely-1)) then
                    if(n == ilower+nVelx) then
                        ABeInd = 1
                    end if
                    do l = ABeInd, AyySize
                        if(AyyRows(l) == n) then
                            c = c + 1
                            cols(c) = AyyCols(l)
                            values(c) = AyyValues(l)
                        else
                            ABeInd = l
                            exit
                        end if
                    end do
                ! the line is in the z-velocity part
                elseif(n >= ilower+nVelx+nVely) then
                    if(n == ilower+nVelx+nVely) then
                        ABeInd = 1
                    end if
                    do l = ABeInd, AzzSize
                        if(AzzRows(l) == n) then
                            c = c + 1
                            cols(c) = AzzCols(l)
                            values(c) = AzzValues(l)
                        else
                            ABeInd = l
                            exit
                        end if
                    end do
                end if
            else
                do l = ABeInd, ASize
                    if(ARows(l) == n) then
                        c = c + 1
                        cols(c) = ACols(l)
                        values(c) = AValues(l)
                    else
                        ABeInd = l
                        exit
                    end if
                end do
            end if

#ifdef LAPACK

            do l = 1, c
                A_lapack(n,cols(l)) = values(l)
            end do

#elif defined(UMFPACK)

            Ap(n+1) = Ap(n) + c
            Ai(Ap(n)+1:Ap(n+1)) = cols(1:c) - 1
            Ax(Ap(n)+1:Ap(n+1)) = values(1:c)

#elif defined(MUMPS)

            mumps_IRN_loc(mumps_NNZ_loc+1:mumps_NNZ_loc+c) = n
            mumps_JCN_loc(mumps_NNZ_loc+1:mumps_NNZ_loc+c) = cols(1:c)
            mumps_A_loc(mumps_NNZ_loc+1:mumps_NNZ_loc+c) = values(1:c)
            mumps_NNZ_loc = mumps_NNZ_loc + c

#elif defined(HYPRE)

            call HYPRE_IJMatrixSetValues(A, 1, c, n, cols, values, ierr)

#endif

        end do

#ifdef LAPACK

        solvertimestart = MPI_Wtime()
        call dgesv(local_x_size, 1, A_lapack, local_x_size, IPIV, b_lapack, local_x_size, LAPACKINFO)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        local_x(:) = b_lapack(:)
        if(LAPACKINFO /= 0) then
            print *, 'LAPACK solver error. Info = ', LAPACKINFO
            stop
        end if

#elif defined(UMFPACK)

        solvertimestart = MPI_Wtime()
        call umf4def(control)
        call umf4sym(local_x_size, local_x_size, Ap, Ai, Ax, symbolic, control, umfinfo)
        call umf4num(Ap, Ai, Ax, symbolic, numeric, control, umfinfo)
        call umf4fsym(symbolic)
        call umf4sol(1, local_x, local_rhs, numeric, control, umfinfo) ! 1 means A'x=b
        if(umfinfo(1) < 0) then
            print *, 'UMFPACK solver error. Info: ', umfinfo(1)
            stop
        end if
        call umf4fnum(numeric)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart

#elif defined(MUMPS)

        mumps_par%NNZ_loc = mumps_NNZ_loc
        allocate(mumps_par%IRN_loc(mumps_par%NNZ_loc))
        allocate(mumps_par%JCN_loc(mumps_par%NNZ_loc))
        allocate(mumps_par%A_loc(mumps_par%NNZ_loc))
        mumps_par%IRN_loc(1:mumps_par%NNZ_loc) = mumps_IRN_loc(1:mumps_NNZ_loc)
        mumps_par%JCN_loc(1:mumps_par%NNZ_loc) = mumps_JCN_loc(1:mumps_NNZ_loc)
        mumps_par%A_loc(1:mumps_par%NNZ_loc) = mumps_A_loc(1:mumps_NNZ_loc)

        mumps_par%JOB = 6
        solvertimestart = MPI_Wtime()
        call DMUMPS(mumps_par)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(mumps_par%INFOG(1) < 0) then
            print *, 'MUMPS solver error. Info = ', mumps_par%INFOG(1)
            stop
        end if

        if(mumps_par%MYID == 0) then
            local_x(:) = mumps_par%RHS(1:local_x_size)
            allocate(requestarray(nProcs-1))
            c = local_x_size + 1
            do i = 1, nProcs-1
                if(mat_kind == 11) then
                    allocate(slave_data(slave_v_data_size(i)))
                    slave_data(:) = mumps_par%RHS(c:c+slave_v_data_size(i)-1)
                    call MPI_IBSEND(slave_data, slave_v_data_size(i), MPI_DOUBLE_PRECISION, i, myid, &!
                        MPI_COMM_WORLD, requestarray(i), ierr)
                    c = c + slave_v_data_size(i)
                    deallocate(slave_data)
                else
                    allocate(slave_data(local_x_size))
                    slave_data(:) = mumps_par%RHS(c:c+local_x_size-1)
                    call MPI_IBSEND(slave_data, local_x_size, MPI_DOUBLE_PRECISION, i, myid, &!
                        MPI_COMM_WORLD, requestarray(i), ierr)
                    c = c + local_x_size
                    deallocate(slave_data)
                end if
            end do
        end if

        if(mumps_par%MYID /= 0) then
            allocate(slave_data(local_x_size))
            call MPI_RECV(slave_data, local_x_size, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
            local_x(:) = slave_data(:)
            deallocate(slave_data)
        end if

        if(mumps_par%MYID == 0) then
            do i = 1, nProcs-1
                call MPI_WAIT(requestarray(i), status, ierr)
            end do
            deallocate(requestarray)
        end if

        deallocate(mumps_par%IRN_loc)
        deallocate(mumps_par%JCN_loc)
        deallocate(mumps_par%A_loc)

#elif defined(HYPRE)

        call HYPRE_IJMatrixAssemble(A, ierr)
        call HYPRE_IJMatrixGetObject(A, parcsr_A, ierr)

        solvertimestart = MPI_Wtime()
        call HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x, ierr)
        call HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b, par_x, ierr)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(ierr /= 0) then
            if(myid == 0) then
                print *, 'HYPRE solver error. ierr = ', ierr, mat_kind
                stop
            end if
        end if

        call HYPRE_IJVectorGetValues(x, local_x_size, rows, local_x, ierr)

        ! let the solution of this time step be the initial x guess in the next time step.
        ! by this way, the number of solver iteration steps can be reduced greatly.
        if(mat_kind == 11) then
            initial_x_guess_v(:) = local_x(:)
        elseif(mat_kind == 12) then
            initial_x_guess_p(:) = local_x(:)
        elseif(mat_kind == 13) then
            initial_x_guess_Cf(:) = local_x(:)
        elseif(mat_kind == 14) then
            initial_x_guess_Tem(:) = local_x(:)
        end if

        deallocate(rows)

#endif

        if(mat_kind == 11) then

            c = 0
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, indexr
                        c = c + 1
                        vx(i,j,k) = local_x(c)
                    end do
                end do
            end do

            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            do k = 1, localnlays
                do j = 1, indexu
                    do i = 1, localncols
                        c = c + 1
                        vy(i,j,k) = local_x(c)
                    end do
                end do
            end do

            if(play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            do k = 1, indexb
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        vz(i,j,k) = local_x(c)
                    end do
                end do
            end do

            call communicateV()

        elseif(mat_kind == 12) then

            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        p(i,j,k) = local_x(c)
                    end do
                end do
            end do

            call communicateP(4)

        elseif(mat_kind == 13) then

            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        Cf(i,j,k) = local_x(c)
                    end do
                end do
            end do

        elseif(mat_kind == 14) then

            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        Tem(i,j,k) = local_x(c)
                    end do
                end do
            end do

        end if

        deallocate(local_rhs)
        if(mat_kind /= 11) then
            deallocate(ACols)
            deallocate(ARows)
            deallocate(AValues)
        end if
        deallocate(cols)
        deallocate(values)
        deallocate(local_x)

#ifdef LAPACK

        deallocate(A_lapack)
        deallocate(b_lapack)
        deallocate(IPIV)

#elif defined(UMFPACK)

        deallocate(Ap)
        deallocate(Ai)
        deallocate(Ax)

#elif defined(MUMPS)

        deallocate(mumps_par%RHS)
        mumps_par%JOB = -2
        call DMUMPS(mumps_par)

        deallocate(mumps_IRN_loc)
        deallocate(mumps_JCN_loc)
        deallocate(mumps_A_loc)

#elif defined(HYPRE)

        call HYPRE_IJMatrixDestroy(A, ierr)
        call HYPRE_IJVectorDestroy(b, ierr)
        call HYPRE_IJVectorDestroy(x, ierr)

        call HYPRE_ParCSRGMRESDestroy(solver, ierr)
        call HYPRE_ParaSailsDestroy(precond, ierr)

#endif

    end subroutine solve

    subroutine outputHisData()

        real(kind=8), dimension(:), allocatable :: local_data, recv
        real(kind=8), dimension(:), pointer :: global_data
        integer :: local_data_size
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        real(kind=8) :: poroavg, Kxxavg, avavg, pavg, Cfavg, Temavg, qavg
        real(kind=8) :: localqsum, localpsum
        integer :: i, j, k, c, ierr

        local_data_size = 6
        allocate(local_data(local_data_size))
        local_data(:) = 0.D0
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    local_data(1) = local_data(1) + poro(i,j,k)
                    local_data(2) = local_data(2) + Kxx(i,j,k)
                    local_data(3) = local_data(3) + av(i,j,k)
                    local_data(4) = local_data(4) + p(i,j,k)
                    local_data(5) = local_data(5) + Cf(i,j,k)
                    local_data(6) = local_data(6) + Tem(i,j,k)
                end do
            end do
        end do

        ! output the average values of poro, Kxx, av, p, Cf, Tem
        if(myid /= 0) then
            call MPI_IBSEND(local_data, local_data_size, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
        end if

        if(myid == 0) then

            allocate(global_data(local_data_size*nProcs))
            allocate(recv(local_data_size))

            global_data(1:local_data_size) = local_data(1:local_data_size)
            c = local_data_size + 1
            do i = 1, nProcs-1
                call MPI_RECV(recv, local_data_size, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                global_data(c:c+local_data_size-1) = recv(:)
                c = c + local_data_size
            end do
            
            poroavg = 0.D0
            do c = 1, local_data_size*nProcs, local_data_size
                poroavg = poroavg + global_data(c)
            end do
            poroavg = poroavg/(nx*ny*nz)

            Kxxavg = 0.D0
            do c = 2, local_data_size*nProcs, local_data_size
                Kxxavg = Kxxavg + global_data(c)
            end do
            Kxxavg = Kxxavg/(nx*ny*nz)

            avavg = 0.D0
            do c = 3, local_data_size*nProcs, local_data_size
                avavg = avavg + global_data(c)
            end do
            avavg = avavg/(nx*ny*nz)

            pavg = 0.D0
            do c = 4, local_data_size*nProcs, local_data_size
                pavg = pavg + global_data(c)
            end do
            pavg = pavg/(nx*ny*nz)

            Cfavg = 0.D0
            do c = 5, local_data_size*nProcs, local_data_size
                Cfavg = Cfavg + global_data(c)
            end do
            Cfavg = Cfavg/(nx*ny*nz)

            Temavg = 0.D0
            do c = 6, local_data_size*nProcs, local_data_size
                Temavg = Temavg + global_data(c)
            end do
            Temavg = Temavg/(nx*ny*nz)

            write(40, fmt='(es24.16)', iostat=ierr) poroavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(41, fmt='(es24.16)', iostat=ierr) Kxxavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(42, fmt='(es24.16)', iostat=ierr) avavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(43, fmt='(es24.16)', iostat=ierr) pavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(44, fmt='(es24.16)', iostat=ierr) Cfavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(45, fmt='(es24.16)', iostat=ierr) Temavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if

            deallocate(global_data)
            deallocate(recv)

        end if

        if(myid /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

        deallocate(local_data)

        ! output Q at the exit
        if(play == 1) then

            localqsum = 0.D0
            do j = 1, localnrows
                do i = 1, localncols
                    localqsum = localqsum + vz(i,j,1)
                end do
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localqsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pncols*pnrows))
            allocate(recv(1))

            global_data(1) = localqsum
            c = 1
            do i = 1, pncols*pnrows-1
                call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                c = c + 1
                global_data(c) = recv(1)
            end do

            qavg = 0.D0
            do c = 1, pncols*pnrows
                qavg = qavg + global_data(c)
            end do
            qavg = qavg/(nx*ny)

            write(46, fmt='(es24.16)', iostat=ierr) qavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if

            deallocate(global_data)
            deallocate(recv)

        end if

        if((play==1).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

        ! output p at the entry
        if(play == pnlays) then

            localpsum = 0.D0
            do j = 1, localnrows
                do i = 1, localncols
                    localpsum = localpsum + p(i,j,localnlays)
                end do
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pncols*pnrows))
            allocate(recv(1))

            if(play == pnlays) then
                global_data(1) = localpsum
                c = 1
                do i = nProcs-pncols*pnrows, nProcs-1
                    if(i > 0) then
                        call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                        c = c + 1
                        global_data(c) = recv(1)
                    end if
                end do
            else
                c = 0
                do i = nProcs-pncols*pnrows, nProcs-1
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end do
            end if

            pavg = 0.D0
            do c = 1, pncols*pnrows
                pavg = pavg + global_data(c)
            end do
            pavg = pavg/(nx*ny)

            write(47, fmt='(es24.16)', iostat=ierr) pavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if

            deallocate(global_data)
            deallocate(recv)

        end if

        if((play==pnlays).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

    end subroutine outputHisData

    subroutine computePresDrop(isBT)

        logical, intent(out) :: isBT
        real(kind=8) :: inpavg, outpavg, localpsum, presDrop
        real(kind=8), dimension(:), allocatable :: global_data, recv
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer :: i, j, k, c, ierr

        if(play == 1) then

            localpsum = 0.D0
            do j = 1, localnrows
                do i = 1, localncols
                    localpsum = localpsum + p(i,j,1)
                end do
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pncols*pnrows))
            allocate(recv(1))

            global_data(1) = localpsum
            c = 1
            do i = 1, pncols*pnrows-1
                call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                c = c + 1
                global_data(c) = recv(1)
            end do

            outpavg = 0.D0
            do c = 1, pncols*pnrows
                outpavg = outpavg + global_data(c)
            end do
            outpavg = outpavg/(nx*ny)

            deallocate(global_data)
            deallocate(recv)

        end if

        if((play==1).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

        if(play == pnlays) then

            localpsum = 0.D0
            do j = 1, localnrows
                do i = 1, localncols
                    localpsum = localpsum + p(i,j,localnlays)
                end do
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pncols*pnrows))
            allocate(recv(1))

            if(play == pnlays) then
                global_data(1) = localpsum
                c = 1
                do i = nProcs-pncols*pnrows, nProcs-1
                    if(i > 0) then
                        call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                        c = c + 1
                        global_data(c) = recv(1)
                    end if
                end do
            else
                c = 0
                do i = nProcs-pncols*pnrows, nProcs-1
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end do
            end if

            inpavg = 0.D0
            do c = 1, pncols*pnrows
                inpavg = inpavg + global_data(c)
            end do
            inpavg = inpavg/(nx*ny)

            deallocate(global_data)
            deallocate(recv)

        end if

        if((play==pnlays).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

        if(myid == 0) then
            isBT = .false.
            presDrop = abs(outpavg-inpavg)
            if(.not.isFindPresDropInit) then
                presDropInit = presDrop
                isFindPresDropInit = .true.
            end if
            if(presDrop/presDropInit < 1.D-2) then
                isBT = .true.
                print *, 'Breakthrough has been achieved! Program stops now.'
                print *, 'The breakthrough time is ', (t-1)*timeEnd/nt, ' seconds.'
                print *, 'The pore volume to breakthrough is ', (t-1)*timeEnd/nt*abs(vzBdryZ1(1,1))
            elseif(mod(t,100) == 0) then
                print *, 'The normalized pressure drop = ', presDrop/presDropInit*1.D2, '%'
            end if
        end if
        if(nProcs > 1) then
            call MPI_BCAST(isBT, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        end if

    end subroutine computePresDrop

    ! output the raw results
    subroutine outputRawData()

        real(kind=8), dimension(:), allocatable :: local_data, recv
        real(kind=8), dimension(:), pointer :: global_data
        integer :: local_data_size, slave_data_size
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer :: indexr, indexu, indexb
        integer :: i, j, k, c, ierr

        local_data_size = local_x_size_v+5*localncols*localnrows*localnlays
        allocate(local_data(local_data_size))

        c = 0
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_data(c) = poro(i,j,k)
                end do
            end do
        end do

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_data(c) = Kxx(i,j,k)
                end do
            end do
        end do

        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, indexr
                    c = c + 1
                    local_data(c) = vx(i,j,k)
                end do
            end do
        end do

        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if
        do k = 1, localnlays
            do j = 1, indexu
                do i = 1, localncols
                    c = c + 1
                    local_data(c) = vy(i,j,k)
                end do
            end do
        end do

        if(play /= pnlays) then
            indexb = localnlays
        else
            indexb = localnlays + 1
        end if
        do k = 1, indexb
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_data(c) = vz(i,j,k)
                end do
            end do
        end do

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_data(c) = p(i,j,k)
                end do
            end do
        end do

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_data(c) = Cf(i,j,k)
                end do
            end do
        end do

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_data(c) = Tem(i,j,k)
                end do
            end do
        end do

        if(myid /= 0) then
            call MPI_IBSEND(local_data, local_data_size, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(myid == 0) then

            allocate(global_data((nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1)+5*nx*ny*nz))

            global_data(1:local_data_size) = local_data(1:local_data_size)
            c = local_data_size + 1
            do i = 1, nProcs-1
                slave_data_size = slave_v_data_size(i)+5*localncols*localnrows*localnlays
                allocate(recv(slave_data_size))
                call MPI_RECV(recv, slave_data_size, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                global_data(c:c+slave_data_size-1) = recv(:)
                c = c + slave_data_size
                deallocate(recv)
            end do

            call exportResults(global_data)

            deallocate(global_data)

        end if

        if(myid /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

        deallocate(local_data)

    end subroutine outputRawData

    subroutine finalize()

        real(kind=8) :: timefinish
        integer :: ierr

        deallocate(xs)
        deallocate(ys)
        deallocate(zs)
        deallocate(ts)
        deallocate(src)
        deallocate(poroInit)
        deallocate(KxxInit)
        deallocate(KyyInit)
        deallocate(KzzInit)
        deallocate(avInit)
        deallocate(vxBdryX0)
        deallocate(vxBdryX1)
        deallocate(vxBdryY0)
        deallocate(vxBdryY1)
        deallocate(vxBdryZ0)
        deallocate(vxBdryZ1)
        deallocate(vyBdryX0)
        deallocate(vyBdryX1)
        deallocate(vyBdryY0)
        deallocate(vyBdryY1)
        deallocate(vyBdryZ0)
        deallocate(vyBdryZ1)
        deallocate(vzBdryX0)
        deallocate(vzBdryX1)
        deallocate(vzBdryY0)
        deallocate(vzBdryY1)
        deallocate(vzBdryZ0)
        deallocate(vzBdryZ1)
        deallocate(isDiriX0_p)
        deallocate(isDiriX1_p)
        deallocate(isDiriY0_p)
        deallocate(isDiriY1_p)
        deallocate(isDiriZ0_p)
        deallocate(isDiriZ1_p)
        deallocate(pBdryX0)
        deallocate(pBdryX1)
        deallocate(pBdryY0)
        deallocate(pBdryY1)
        deallocate(pBdryZ0)
        deallocate(pBdryZ1)
        deallocate(pInit)
        deallocate(isDiriX0_Cf)
        deallocate(isDiriX1_Cf)
        deallocate(isDiriY0_Cf)
        deallocate(isDiriY1_Cf)
        deallocate(isDiriZ0_Cf)
        deallocate(isDiriZ1_Cf)
        deallocate(CfBdryX0)
        deallocate(CfBdryX1)
        deallocate(CfBdryY0)
        deallocate(CfBdryY1)
        deallocate(CfBdryZ0)
        deallocate(CfBdryZ1)
        deallocate(CfInit)
        deallocate(isDiriX0_Tem)
        deallocate(isDiriX1_Tem)
        deallocate(isDiriY0_Tem)
        deallocate(isDiriY1_Tem)
        deallocate(isDiriZ0_Tem)
        deallocate(isDiriZ1_Tem)
        deallocate(TemBdryX0)
        deallocate(TemBdryX1)
        deallocate(TemBdryY0)
        deallocate(TemBdryY1)
        deallocate(TemBdryZ0)
        deallocate(TemBdryZ1)
        deallocate(TemInit)

        deallocate(dm)
        deallocate(kc)
        deallocate(ks)
        deallocate(hx)
        deallocate(hy)
        deallocate(hz)
        deallocate(poro)
        deallocate(poro_old)
        deallocate(poroHarmX)
        deallocate(poroHarmX_old)
        deallocate(poroHarmXInit)
        deallocate(poroHarmY)
        deallocate(poroHarmY_old)
        deallocate(poroHarmYInit)
        deallocate(poroHarmZ)
        deallocate(poroHarmZ_old)
        deallocate(poroHarmZInit)
        deallocate(Kxx)
        deallocate(Kyy)
        deallocate(Kzz)
        deallocate(KxxHarm)
        deallocate(KyyHarm)
        deallocate(KzzHarm)
        deallocate(av)
        deallocate(vx)
        deallocate(vy)
        deallocate(vz)
        deallocate(p)
        deallocate(Cf)
        deallocate(Tem)

        deallocate(local_rhs_static_v)
        deallocate(local_rhs_v)
        deallocate(local_rhs_p)
        deallocate(local_rhs_Cf)
        deallocate(local_rhs_Tem)
        deallocate(AxxCols)
        deallocate(AxxRows)
        deallocate(AxxStaticValues)
        deallocate(AxxDynValues)
        deallocate(AxxValues)
        deallocate(AyyCols)
        deallocate(AyyRows)
        deallocate(AyyStaticValues)
        deallocate(AyyDynValues)
        deallocate(AyyValues)
        deallocate(AzzCols)
        deallocate(AzzRows)
        deallocate(AzzStaticValues)
        deallocate(AzzDynValues)
        deallocate(AzzValues)
        deallocate(ApCols)
        deallocate(ApRows)
        deallocate(ApValues)
        deallocate(AcfCols)
        deallocate(AcfRows)
        deallocate(AcfValues)
        deallocate(AtemCols)
        deallocate(AtemRows)
        deallocate(AtemValues)
        deallocate(AxxEntryNum)
        deallocate(AyyEntryNum)
        deallocate(AzzEntryNum)
        deallocate(ApEntryNum)
        deallocate(AcfEntryNum)
        deallocate(AtemEntryNum)
        deallocate(AxxEntryBase)
        deallocate(AyyEntryBase)
        deallocate(AzzEntryBase)
        deallocate(ApEntryBase)
        deallocate(AcfEntryBase)
        deallocate(AtemEntryBase)

        if((nProcs>1).and.(myid == 0)) then
            deallocate(slave_v_data_size)
        end if

        if(myid == 0) then
            close(40)
            close(41)
            close(42)
            close(43)
            close(44)
            close(45)
            close(46)
            close(47)
        end if

#ifdef HYPRE

        deallocate(initial_x_guess_v)
        deallocate(initial_x_guess_p)
        deallocate(initial_x_guess_Cf)
        deallocate(initial_x_guess_Tem)

#endif

#if defined(MUMPS) || defined(HYPRE)
        call MPI_BUFFER_DETACH(buffer, buffer_size, ierr)
#endif

        if(myid == 0) then
            print *, 'Solver time = ', solvertime, ' seconds.'
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timefinish = MPI_Wtime()
        if(myid == 0) then
            print *, 'Elapsed time = ', timefinish-timestart, ' seconds.'
        end if

        call MPI_Finalize(ierr)

    end subroutine finalize

    subroutine driver(modelCase)

        type(model), intent(in out) :: modelCase
        logical :: isBT

        call initialize(modelCase)

        call genStaticPara(151)
        call genStaticPara(253)
        call genStaticPara(355)

        ! time iteration
        do t = 2, nt+1

            call computedm()
            call computeks()
            if(t == 2) then
                call computekc()
            end if
            call computePoroHarm(1) ! old
            call computePoro()
            call computePoroHarm(2) ! new

            call computeK()
            call computeKHarm()

            call computeav()

            call genDynPara(457)
            call solve(12)

            call genDynPara(152)
            call genDynPara(254)
            call genDynPara(356)
            call solve(11)

            call computekc()
            call genDynPara(458)
            call solve(13)

            call genDynPara(459)
            call solve(14)

            call outputHisData()
            call computePresDrop(isBT)
            if(isBT) then
                call outputRawData()
                exit
            elseif(t == nt+1) then
                if(myid == 0) then
                    print *, 'Program stops now, but breakthrough has NOT been achieved.'
                    print *, 'More simulation time is needed to achieve breakthrough!'
                end if
            end if
            if((t==2).or.(mod(t-1,nt/NUMFRAME)==0)) then
                call outputRawData()
            end if

        end do

        call finalize()

    end subroutine driver

end module DBF_driver

