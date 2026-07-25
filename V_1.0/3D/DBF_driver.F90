
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

        AxxEntryNum(:) = 7
        AxpEntryNum(:) = 2
        AyyEntryNum(:) = 7
        AypEntryNum(:) = 2
        AzzEntryNum(:) = 7
        AzpEntryNum(:) = 2
        AcxEntryNum(:) = 2
        AcyEntryNum(:) = 2
        AczEntryNum(:) = 2
        AcpEntryNum(:) = 1
        AcfEntryNum(:) = 19

        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if

        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call equCoorditoInd(1, 1, j, k, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call equCoorditoInd(1, localncols+1, j, k, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, indexr
                    call equCoorditoInd(1, i, 1, k, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, indexr
                    call equCoorditoInd(1, i, localnrows, k, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, indexr
                    call equCoorditoInd(1, i, j, 1, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, indexr
                    call equCoorditoInd(1, i, j, localnlays, eq_ind)
                    AxxEntryNum(eq_ind) = AxxEntryNum(eq_ind) - 1
                end do
            end do
        end if
        do k = 1, localnlays
            do j = 1, localnrows
                if((pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                    call equCoorditoInd(1, 1, j, k, eq_ind)
                    AxxEntryNum(eq_ind) = 1
                end if
                if((pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                    call equCoorditoInd(1, localncols+1, j, k, eq_ind)
                    AxxEntryNum(eq_ind) = 1
                end if
            end do
        end do

        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call equCoorditoInd(1, 1, j, k, eq_ind)
                    if((pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                        AxpEntryNum(eq_ind) = 0
                    else
                        AxpEntryNum(eq_ind) = 1
                    end if
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call equCoorditoInd(1, localncols+1, j, k, eq_ind)
                    if((pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                        AxpEntryNum(eq_ind) = 0
                    else
                        AxpEntryNum(eq_ind) = 1
                    end if
                end do
            end do
        end if
       
        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if

        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, indexu
                    call equCoorditoInd(2, 1, j, k, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, indexu
                    call equCoorditoInd(2, localncols, j, k, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    call equCoorditoInd(2, i, 1, k, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    call equCoorditoInd(2, i, localnrows+1, k, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, indexu
                do i = 1, localncols
                    call equCoorditoInd(2, i, j, 1, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, indexu
                do i = 1, localncols
                    call equCoorditoInd(2, i, j, localnlays, eq_ind)
                    AyyEntryNum(eq_ind) = AyyEntryNum(eq_ind) - 1
                end do
            end do
        end if
        do k = 1, localnlays
            do i = 1, localncols
                if((prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                    call equCoorditoInd(2, i, 1, k, eq_ind)
                    AyyEntryNum(eq_ind) = 1
                end if
                if((prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                    call equCoorditoInd(2, i, localnrows+1, k, eq_ind)
                    AyyEntryNum(eq_ind) = 1
                end if
            end do
        end do

        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    call equCoorditoInd(2, i, 1, k, eq_ind)
                    if((prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                        AypEntryNum(eq_ind) = 0
                    else
                        AypEntryNum(eq_ind) = 1
                    end if
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    call equCoorditoInd(2, i, localnrows+1, k, eq_ind)
                    if((prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                        AypEntryNum(eq_ind) = 0
                    else
                        AypEntryNum(eq_ind) = 1
                    end if
                end do
            end do
        end if

        if(play /= pnlays) then
            indexb = localnlays
        else
            indexb = localnlays + 1
        end if

        if(pcol == 1) then
            do k = 1, indexb
                do j = 1, localnrows
                    call equCoorditoInd(3, 1, j, k, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, indexb
                do j = 1, localnrows
                    call equCoorditoInd(3, localncols, j, k, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, indexb
                do i = 1, localncols
                    call equCoorditoInd(3, i, 1, k, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, indexb
                do i = 1, localncols
                    call equCoorditoInd(3, i, localnrows, k, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    call equCoorditoInd(3, i, j, 1, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    call equCoorditoInd(3, i, j, localnlays+1, eq_ind)
                    AzzEntryNum(eq_ind) = AzzEntryNum(eq_ind) - 1
                end do
            end do
        end if
        do j = 1, localnrows
            do i = 1, localncols
                if((play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                    call equCoorditoInd(3, i, j, 1, eq_ind)
                    AzzEntryNum(eq_ind) = 1
                end if
                if((play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                    call equCoorditoInd(3, i, j, localnlays+1, eq_ind)
                    AzzEntryNum(eq_ind) = 1
                end if
            end do
        end do

        if(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    call equCoorditoInd(3, i, j, 1, eq_ind)
                    if((play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                        AzpEntryNum(eq_ind) = 0
                    else
                        AzpEntryNum(eq_ind) = 1
                    end if
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    call equCoorditoInd(3, i, j, localnlays+1, eq_ind)
                    if((play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                        AzpEntryNum(eq_ind) = 0
                    else
                        AzpEntryNum(eq_ind) = 1
                    end if
                end do
            end do
        end if

        if(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call equCoorditoInd(4, 1, j, k, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    call equCoorditoInd(4, localncols, j, k, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    call equCoorditoInd(4, i, 1, k, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    call equCoorditoInd(4, i, localnrows, k, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    call equCoorditoInd(4, i, j, 1, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    call equCoorditoInd(4, i, j, localnlays, eq_ind)
                    AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) - 5
                end do
            end do
        end if
        if((pcol==1).and.(prow==1)) then
            do k = 1, localnlays
                call equCoorditoInd(4, 1, 1, k, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==1).and.(prow==pnrows)) then
            do k = 1, localnlays
                call equCoorditoInd(4, 1, localnrows, k, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==pncols).and.(prow==1)) then
            do k = 1, localnlays
                call equCoorditoInd(4, localncols, 1, k, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==pncols).and.(prow==pnrows)) then
            do k = 1, localnlays
                call equCoorditoInd(4, localncols, localnrows, k, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==1).and.(play==1)) then
            do j = 1, localnrows
                call equCoorditoInd(4, 1, j, 1, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==1).and.(play==pnlays)) then
            do j = 1, localnrows
                call equCoorditoInd(4, 1, j, localnlays, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==pncols).and.(play==1)) then
            do j = 1, localnrows
                call equCoorditoInd(4, localncols, j, 1, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((pcol==pncols).and.(play==pnlays)) then
            do j = 1, localnrows
                call equCoorditoInd(4, localncols, j, localnlays, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((prow==1).and.(play==1)) then
            do i = 1, localncols
                call equCoorditoInd(4, i, 1, 1, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((prow==1).and.(play==pnlays)) then
            do i = 1, localncols
                call equCoorditoInd(4, i, 1, localnlays, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((prow==pnrows).and.(play==1)) then
            do i = 1, localncols
                call equCoorditoInd(4, i, localnrows, 1, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if
        if((prow==pnrows).and.(play==pnlays)) then
            do i = 1, localncols
                call equCoorditoInd(4, i, localnrows, localnlays, eq_ind)
                AcfEntryNum(eq_ind) = AcfEntryNum(eq_ind) + 1
            end do
        end if

    end subroutine computeMatEntryNum

    subroutine initialize(modelCase)

        type(model), intent(in out) :: modelCase

        integer :: xmomeSize, ymomeSize, zmomeSize, contiSize, cfSize
        ! The arrays that will communicate between different functions must use pointer type instead of allocatalbe
        ! type. Pointer type can make sure that the subscripts of the arrays keep the same in the calling and called
        ! functions. However, if the subscripts of the arrays begin with 1, using allocatable type is also OK.
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: global_ind_b, global_ind_e
        logical :: alive
        integer :: i, j, k, n, c, ierr

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
        epslon = modelCase%epslon
        dm = modelCase%dm
        alphaOS = modelCase%alphaOS
        lamdaX = modelCase%lamdaX
        lamdaT = modelCase%lamdaT
        radiusInit = modelCase%radiusInit
        ShInfinity = modelCase%ShInfinity
        ks = modelCase%ks
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

        call index_convert_local_global(myid, 1, 1, 1, 1, ilower)
        call index_convert_local_global(myid, 4, localncols, localnrows, localnlays, iupper)
        local_x_size = iupper - ilower + 1

        if((nProcs>1).and.(myid==0)) then
            allocate(slave_vp_data_size(nProcs-1))
            do i = 1, nProcs-1
                call index_convert_local_global(i, 1, 1, 1, 1, global_ind_b)
                call index_convert_local_global(i, 4, localncols, localnrows, localnlays, global_ind_e)
                slave_vp_data_size(i) = global_ind_e-global_ind_b+1
            end do
        end if

        call index_convert_local_global(myid, 5, 1, 1, 1, ilower_Cf)
        call index_convert_local_global(myid, 5, localncols, localnrows, localnlays, iupper_Cf)
        local_x_size_Cf = iupper_Cf - ilower_Cf + 1

        t = 2

        ! initialize kc
        allocate(kc(1:localncols, 1:localnrows, 1:localnlays))

        ! initialize hx, hy, hz
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
        allocate(hx(indexl:indexr))
        do i = indexl, indexr
            hx(i) = xs(xlower+i) - xs(xlower+i-1)
        end do

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
        allocate(hy(indexd:indexu))
        do j = indexd, indexu
            hy(j) = ys(ylower+j) - ys(ylower+j-1)
        end do

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
        allocate(hz(indexf:indexb))
        do k = indexf, indexb
            hz(k) = zs(zlower+k) - zs(zlower+k-1)
        end do
        
        ! initialize poro
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

        ! initialize poroHarm
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

        ! initialize KHarm
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
        if(play /= pnlays) then
            indexb = localnlays
        else
            indexb = localnlays + 1
        end if

        allocate(KxxHarm(1:indexr,1:localnrows,1:localnlays))
        allocate(KyyHarm(1:localncols,1:indexu,1:localnlays))
        allocate(KzzHarm(1:localncols,1:localnrows,1:indexb))

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
                    if(isDiriY0_p(i,k) == 0) then
                        vy(i,1,k) = vyBdryY0(xlower+i-1,zlower+k-1)
                    end if
                end do
            end do
        end if
        if(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    if(isDiriY1_p(i,k) == 0) then
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
                    if(isDiriZ0_p(i,j) == 0) then
                        vz(i,j,1) = vzBdryZ0(xlower+i-1,ylower+j-1)
                    end if
                end do
            end do
        end if
        if(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    if(isDiriZ1_p(i,j) == 0) then
                        vz(i,j,localnlays+1) = vzBdryZ1(xlower+i-1,ylower+j-1)
                    end if
                end do
            end do
        end if

        ! initialize p
        allocate(p(1:localncols,1:localnrows,1:localnlays))
        p(1:localncols,1:localnrows,1:localnlays) = pInit(xlower:xupper,ylower:yupper,zlower:zupper)

        ! initialize Cf
        allocate(Cf(1:localncols,1:localnrows,1:localnlays))
        Cf(1:localncols,1:localnrows,1:localnlays) = CfInit(xlower:xupper,ylower:yupper,zlower:zupper)

        ! compute the size of the equations
        if(pcol /= pncols) then
            xmomeSize = localncols*localnrows*localnlays
        else
            xmomeSize = (localncols+1)*localnrows*localnlays
        end if
        if(prow /= pnrows) then
            ymomeSize = localncols*localnrows*localnlays
        else
            ymomeSize = localncols*(localnrows+1)*localnlays
        end if
        if(play /= pnlays) then
            zmomeSize = localncols*localnrows*localnlays
        else
            zmomeSize = localncols*localnrows*(localnlays+1)
        end if
        contiSize = localncols*localnrows*localnlays
        cfSize = localncols*localnrows*localnlays

        allocate(AxxEntryNum(xmomeSize))
        allocate(AxpEntryNum(xmomeSize))
        allocate(AyyEntryNum(ymomeSize))
        allocate(AypEntryNum(ymomeSize))
        allocate(AzzEntryNum(zmomeSize))
        allocate(AzpEntryNum(zmomeSize))
        allocate(AcxEntryNum(contiSize))
        allocate(AcyEntryNum(contiSize))
        allocate(AczEntryNum(contiSize))
        allocate(AcpEntryNum(contiSize))
        allocate(AcfEntryNum(cfSize))

        call computeMatEntryNum()

        allocate(AxxEntryBase(xmomeSize))
        allocate(AxpEntryBase(xmomeSize))
        allocate(AyyEntryBase(ymomeSize))
        allocate(AypEntryBase(ymomeSize))
        allocate(AzzEntryBase(zmomeSize))
        allocate(AzpEntryBase(zmomeSize))
        allocate(AcxEntryBase(contiSize))
        allocate(AcyEntryBase(contiSize))
        allocate(AczEntryBase(contiSize))
        allocate(AcpEntryBase(contiSize))
        allocate(AcfEntryBase(cfSize))

        AxxEntryBase(1) = 1
        AxpEntryBase(1) = 1
        do n = 2, xmomeSize
            AxxEntryBase(n) = AxxEntryBase(n-1) + AxxEntryNum(n-1)
            AxpEntryBase(n) = AxpEntryBase(n-1) + AxpEntryNum(n-1)
        end do

        AyyEntryBase(1) = 1
        AypEntryBase(1) = 1
        do n = 2, ymomeSize
            AyyEntryBase(n) = AyyEntryBase(n-1) + AyyEntryNum(n-1)
            AypEntryBase(n) = AypEntryBase(n-1) + AypEntryNum(n-1)
        end do

        AzzEntryBase(1) = 1
        AzpEntryBase(1) = 1
        do n = 2, zmomeSize
            AzzEntryBase(n) = AzzEntryBase(n-1) + AzzEntryNum(n-1)
            AzpEntryBase(n) = AzpEntryBase(n-1) + AzpEntryNum(n-1)
        end do

        AcxEntryBase(1) = 1
        AcyEntryBase(1) = 1
        AczEntryBase(1) = 1
        AcpEntryBase(1) = 1
        do n = 2, contiSize
            AcxEntryBase(n) = AcxEntryBase(n-1) + AcxEntryNum(n-1)
            AcyEntryBase(n) = AcyEntryBase(n-1) + AcyEntryNum(n-1)
            AczEntryBase(n) = AczEntryBase(n-1) + AczEntryNum(n-1)
            AcpEntryBase(n) = AcpEntryBase(n-1) + AcpEntryNum(n-1)
        end do

        AcfEntryBase(1) = 1
        do n = 2, cfSize
            AcfEntryBase(n) = AcfEntryBase(n-1) + AcfEntryNum(n-1)
        end do

        ! compute matrix size
        AxxSize = 0
        AxpSize = 0
        do n = 1, xmomeSize
            AxxSize = AxxSize + AxxEntryNum(n)
            AxpSize = AxpSize + AxpEntryNum(n)
        end do

        AyySize = 0
        AypSize = 0
        do n = 1, ymomeSize
            AyySize = AyySize + AyyEntryNum(n)
            AypSize = AypSize + AypEntryNum(n)
        end do
        
        AzzSize = 0
        AzpSize = 0
        do n = 1, zmomeSize
            AzzSize = AzzSize + AzzEntryNum(n)
            AzpSize = AzpSize + AzpEntryNum(n)
        end do

        AcxSize = 2*localncols*localnrows*localnlays
        AcySize = 2*localncols*localnrows*localnlays
        AczSize = 2*localncols*localnrows*localnlays
        AcpSize = localncols*localnrows*localnlays

        AcfSize = 0
        do n = 1, cfSize
            AcfSize = AcfSize + AcfEntryNum(n)
        end do

        ! initialize matrix
        allocate(local_rhs(local_x_size))
        allocate(local_rhs_static(local_x_size))
        allocate(local_rhs_Cf(local_x_size_Cf))

        allocate(AxxCols(AxxSize))
        allocate(AxxRows(AxxSize))
        allocate(AxxStaticValues(AxxSize))
        allocate(AxxDynValues(AxxSize))
        allocate(AxxValues(AxxSize))
        allocate(AxpCols(AxpSize))
        allocate(AxpRows(AxpSize))
        allocate(AxpValues(AxpSize))
        allocate(AyyCols(AyySize))
        allocate(AyyRows(AyySize))
        allocate(AyyStaticValues(AyySize))
        allocate(AyyDynValues(AyySize))
        allocate(AyyValues(AyySize))
        allocate(AypCols(AypSize))
        allocate(AypRows(AypSize))
        allocate(AypValues(AypSize))
        allocate(AzzCols(AzzSize))
        allocate(AzzRows(AzzSize))
        allocate(AzzStaticValues(AzzSize))
        allocate(AzzDynValues(AzzSize))
        allocate(AzzValues(AzzSize))
        allocate(AzpCols(AzpSize))
        allocate(AzpRows(AzpSize))
        allocate(AzpValues(AzpSize))
        allocate(AcxCols(AcxSize))
        allocate(AcxRows(AcxSize))
        allocate(AcxValues(AcxSize))
        allocate(AcyCols(AcySize))
        allocate(AcyRows(AcySize))
        allocate(AcyValues(AcySize))
        allocate(AczCols(AczSize))
        allocate(AczRows(AczSize))
        allocate(AczValues(AczSize))
        allocate(AcpCols(AcpSize))
        allocate(AcpRows(AcpSize))
        allocate(AcpValues(AcpSize))
        allocate(AcfCols(AcfSize))
        allocate(AcfRows(AcfSize))
        allocate(AcfValues(AcfSize))

        AxxCols(:) = 0
        AxxRows(:) = 0
        AxxStaticValues(:) = 0.D0
        AxxDynValues(:) = 0.D0
        AxxValues(:) = 0.D0
        AxpCols(:) = 0
        AxpRows(:) = 0
        AxpValues(:) = 0.D0
        AyyCols(:) = 0
        AyyRows(:) = 0
        AyyStaticValues(:) = 0.D0
        AyyDynValues(:) = 0.D0
        AyyValues(:) = 0.D0
        AypCols(:) = 0
        AypRows(:) = 0
        AypValues(:) = 0.D0
        AzzCols(:) = 0
        AzzRows(:) = 0
        AzzStaticValues(:) = 0.D0
        AzzDynValues(:) = 0.D0
        AzzValues(:) = 0.D0
        AzpCols(:) = 0
        AzpRows(:) = 0
        AzpValues(:) = 0.D0
        AcxCols(:) = 0
        AcxRows(:) = 0
        AcxValues(:) = 0.D0
        AcyCols(:) = 0
        AcyRows(:) = 0
        AcyValues(:) = 0.D0
        AczCols(:) = 0
        AczRows(:) = 0
        AczValues(:) = 0.D0
        AcpCols(:) = 0
        AcpRows(:) = 0
        AcpValues(:) = 0.D0
        AcfCols(:) = 0
        AcfRows(:) = 0
        AcfValues(:) = 0.D0

        presDropInit = 0.D0
        isFindPresDropInit = .false.

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
            open(unit=45, file=trim(adjustl(soludoc))//'/his_q_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_q_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=46, file=trim(adjustl(soludoc))//'/his_lp_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_lp_avg.txt', ' error. ', ierr
                stop
            end if
        end if

#ifdef LAPACK

        allocate(A_lapack(local_x_size,local_x_size))
        allocate(b_lapack(local_x_size))
        allocate(IPIV(local_x_size))

        allocate(A_lapack_Cf(local_x_size_Cf,local_x_size_Cf))
        allocate(b_lapack_Cf(local_x_size_Cf))
        allocate(IPIV_Cf(local_x_size_Cf))

#elif defined(UMFPACK)

        allocate(Ap(1:local_x_size+1))
        allocate(Ai(1:9*local_x_size))
        allocate(Ax(1:9*local_x_size))

        allocate(Ap_Cf(1:local_x_size_Cf+1))
        allocate(Ai_Cf(1:19*local_x_size_Cf))
        allocate(Ax_Cf(1:19*local_x_size_Cf))

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
            mumps_par%N = (nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1)+nx*ny*nz
        end if
        allocate(mumps_par%RHS(mumps_par%N))

        allocate(mumps_IRN_loc(local_x_size*9))
        allocate(mumps_JCN_loc(local_x_size*9))
        allocate(mumps_A_loc(local_x_size*9))

        mumps_par_Cf%COMM = MPI_COMM_WORLD
        mumps_par_Cf%SYM = 0
        mumps_par_Cf%PAR = 1
        mumps_par_Cf%JOB = -1
        call DMUMPS(mumps_par_Cf)

        mumps_par_Cf%ICNTL(4) = 0
        mumps_par_Cf%ICNTL(7) = 15
        mumps_par_Cf%ICNTL(18) = 3
        if(mumps_par_Cf%MYID == 0) then
            mumps_par_Cf%N = nx*ny*nz
        end if
        allocate(mumps_par_Cf%RHS(mumps_par_Cf%N))

        allocate(mumps_IRN_loc_Cf(local_x_size_Cf*19))
        allocate(mumps_JCN_loc_Cf(local_x_size_Cf*19))
        allocate(mumps_A_loc_Cf(local_x_size_Cf*19))

#elif defined(HYPRE)

        if(play /= 1) then
            call index_convert_local_global(myid-pncols*pnrows, 1, 1, 1, localnlays, jlower)
        elseif(prow /= 1) then
            call index_convert_local_global(myid-pncols, 1, 1, localnrows, 1, jlower)
        elseif(pcol /= 1) then
            call index_convert_local_global(myid-1, 1, localncols, 1, 1, jlower)
        else
            jlower = ilower
        end if
        if(play /= pnlays) then
            call index_convert_local_global(myid+pncols*pnrows, 4, localncols, localnrows, 1, jupper)
        elseif(prow /= pnrows) then
            call index_convert_local_global(myid+pncols, 4, localncols, 1, localnlays, jupper)
        elseif(pcol /= pncols) then
            call index_convert_local_global(myid+1, 4, 1, localnrows, localnlays, jupper)
        else
            jupper = iupper
        end if

        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, jlower, jupper, A, ierr)
        call HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR, ierr)
        call HYPRE_IJMatrixInitialize(A, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, b, ierr)
        call HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(b, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, x, ierr)
        call HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(x, ierr)

        if((prow /= 1).and.(play /= 1)) then
            call index_convert_local_global(myid-pncols-pncols*pnrows, 5, 1, localnrows, localnlays, jlower_Cf)
        elseif((pcol /= 1).and.(play /= 1)) then
            call index_convert_local_global(myid-1-pncols*pnrows, 5, localncols, 1, localnlays, jlower_Cf)
        elseif((pcol /= 1).and.(prow /= 1)) then
            call index_convert_local_global(myid-1-pncols, 5, localncols, localnrows, 1, jlower_Cf)
        elseif(play /= 1) then
            call index_convert_local_global(myid-pncols*pnrows, 5, 1, 1, localnlays, jlower_Cf)
        elseif(prow /= 1) then
            call index_convert_local_global(myid-pncols, 5, 1, localnrows, 1, jlower_Cf)
        elseif(pcol /= 1) then
            call index_convert_local_global(myid-1, 5, localncols, 1, 1, jlower_Cf)
        else
            jlower_Cf = ilower_Cf
        end if
        if((prow /= pnrows).and.(play /= pnlays)) then
            call index_convert_local_global(myid+pncols+pncols*pnrows, 5, localncols, 1, 1, jupper_Cf)
        elseif((pcol /= pncols).and.(play /= pnlays)) then
            call index_convert_local_global(myid+1+pncols*pnrows, 5, 1, localnrows, 1, jupper_Cf)
        elseif((pcol /= pncols).and.(prow /= pnrows)) then
            call index_convert_local_global(myid+1+pncols, 5, 1, 1, localnlays, jupper_Cf)
        elseif(play /= pnlays) then
            call index_convert_local_global(myid+pncols*pnrows, 5, localncols, localnrows, 1, jupper_Cf)
        elseif(prow /= pnrows) then
            call index_convert_local_global(myid+pncols, 5, localncols, 1, localnlays, jupper_Cf)
        elseif(pcol /= pncols) then
            call index_convert_local_global(myid+1, 5, 1, localnrows, localnlays, jupper_Cf)
        else
            jupper_Cf = iupper_Cf
        end if

        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower_Cf, iupper_Cf, jlower_Cf, jupper_Cf, A_Cf, ierr)
        call HYPRE_IJMatrixSetObjectType(A_Cf, HYPRE_PARCSR, ierr)
        call HYPRE_IJMatrixInitialize(A_Cf, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Cf, iupper_Cf, b_Cf, ierr)
        call HYPRE_IJVectorSetObjectType(b_Cf, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(b_Cf, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Cf, iupper_Cf, x_Cf, ierr)
        call HYPRE_IJVectorSetObjectType(x_Cf, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(x_Cf, ierr)

        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver, ierr)
        call HYPRE_ParCSRGMRESSetMaxIter(solver, 1000, ierr)
        call HYPRE_ParCSRGMRESSetTol(solver, 1.0d-7, ierr)
        !call HYPRE_ParCSRGMRESSetPrintLevel(solver, 2, ierr)
        !call HYPRE_ParCSRGMRESSetLogging(solver, 1, ierr)

        call HYPRE_ParaSailsCreate(MPI_COMM_WORLD, precond, ierr)
        call HYPRE_ParaSailsSetParams(precond, -0.9, 2, ierr)
        call HYPRE_ParaSailsSetFilter(precond, -0.9, ierr)
        ! Because the matrix A is nonsymmetric and indefinite, you must choose the parameter as 0.
        call HYPRE_ParaSailsSetSym(precond, 0)
        call HYPRE_ParaSailsSetLogging(precond, 1, ierr)

        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver, 1, precond, ierr)

        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver_Cf, ierr)
        call HYPRE_ParCSRGMRESSetMaxIter(solver_Cf, 1000, ierr)
        call HYPRE_ParCSRGMRESSetTol(solver_Cf, 1.0d-7, ierr)
        !call HYPRE_ParCSRGMRESSetPrintLevel(solver_Cf, 2, ierr)
        !call HYPRE_ParCSRGMRESSetLogging(solver_Cf, 1, ierr)

        call HYPRE_ParaSailsCreate(MPI_COMM_WORLD, precond_Cf,ierr)
        call HYPRE_ParaSailsSetParams(precond_Cf, -0.9, 2, ierr)
        call HYPRE_ParaSailsSetFilter(precond_Cf, -0.9, ierr)
        call HYPRE_ParaSailsSetLoadbal(precond_Cf, 0.9, ierr)
        ! Because the matrix A is nonsymmetric and indefinite, you must choose the parameter as 0.
        call HYPRE_ParaSailsSetSym(precond_Cf, 0)
        call HYPRE_ParaSailsSetLogging(precond_Cf, 1, ierr)

        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver_Cf, 1, precond_Cf, ierr)

        allocate(rows(local_x_size))
        do n = 1, local_x_size
            rows(n) = ilower + n - 1
        end do

        allocate(rows_Cf(local_x_size_Cf))
        do n = 1, local_x_size_Cf
            rows_Cf(n) = ilower_Cf + n - 1
        end do

        ! the initial guess of x in the solver iteration
        allocate(initial_x_guess(local_x_size))
        initial_x_guess(:) = 0.D0
        c = local_x_size - localncols*localnrows*localnlays
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    initial_x_guess(c) = pInit(xlower+i-1,ylower+j-1,zlower+k-1)
                end do
            end do
        end do

        allocate(initial_x_guess_Cf(local_x_size_Cf))
        c = 0
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    initial_x_guess_Cf(c) = CfInit(xlower+i-1,ylower+j-1,zlower+k-1)
                end do
            end do
        end do

#endif

        if(mod(nt, NUMFRAME) /= 0) then
            print *, 'The number of time steps must be divided by the number of frames.'
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

    end subroutine initialize

    ! Generate the values on the right-hand side and the coefficients of the matrix A,
    ! and the subroutine will generate the values that will not change with the time iteration.
    subroutine genStaticPara_vp()

        integer :: findexl, findexr, findexd, findexu, findexf, findexb ! field index
        integer :: eindexl, eindexr, eindexd, eindexu, eindexf, eindexb ! equation index
        integer :: global_ind
        integer, dimension(:,:,:), pointer :: velx, vely, velz, pres
        logical :: isField
        real(kind=8), dimension(:,:,:), pointer :: rhs_velx_b, rhs_dpdx, rhs_vely_b, rhs_dpdy, rhs_velz_b, rhs_dpdz, rhs_dudx, &!
            rhs_dvdy, rhs_dwdz
        real(kind=8), dimension(:,:,:), pointer :: resiAxx_b, resiAxp, resiAyy_b, resiAyp, resiAzz_b, resiAzp, resiAcx, resiAcy, &!
            resiAcz, resitemp
        integer :: i, j, k, n, c

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
        allocate(velx(findexl:findexr,findexd:findexu,findexf:findexb))
        velx(:,:,:) = 0

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
        allocate(rhs_velx_b(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_velx_b(velx, rhs_velx_b)
        allocate(resiAxx_b(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, velx, isField)
                    if(isField) then
                        call Resi_velx_b(velx, resiAxx_b)
                        resitemp = resiAxx_b - rhs_velx_b
                        call constructAxx(velx, resitemp, 1)
                    end if
                end do
            end do
        end do
        deallocate(velx)
        deallocate(resiAxx_b)
        deallocate(resitemp)

        do k = 1, localnlays
            do j = 1, localnrows
                if((pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                    call index_convert_local_global(myid, 1, 1, j, k, global_ind)
                    do n = 1, AxxSize
                        if(AxxRows(n)==global_ind) then
                            AxxValues(n) = AxxStaticValues(n)
                        end if
                    end do
                end if
                if((pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                    call index_convert_local_global(myid, 1, localncols+1, j, k, global_ind)
                    do n = 1, AxxSize
                        if(AxxRows(n)==global_ind) then
                            AxxValues(n) = AxxStaticValues(n)
                        end if
                    end do
                end if
            end do
        end do

        if(pcol /= 1) then
            findexl = 0
        else
            findexl = 1
        end if
        findexr = localncols
        findexd = 1
        findexu = localnrows
        findexf = 1
        findexb = localnlays
        allocate(pres(findexl:findexr,findexd:findexu,findexf:findexb))
        pres(:,:,:) = 0
        allocate(rhs_dpdx(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_dpdx(pres, rhs_dpdx)
        allocate(resiAxp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, pres, isField)
                    if(isField) then
                        call Resi_dpdx(pres, resiAxp)
                        resitemp = resiAxp - rhs_dpdx
                        call constructAxp(pres, resitemp)
                    end if
                end do
            end do
        end do
        deallocate(pres)
        deallocate(resiAxp)
        deallocate(resitemp)

        c = 0
        do k = eindexf, eindexb
            do j = eindexd, eindexu
                do i = eindexl, eindexr
                    c = c + 1
                    local_rhs_static(c) = -(rhs_velx_b(i,j,k) + rhs_dpdx(i,j,k) - rhof*gravX)
                    if(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)).or.((pcol==pncols).and. &!
                        (i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0))) then
                        local_rhs(c) = local_rhs_static(c) - rhof*gravX
                    end if
                end do
            end do
        end do
        deallocate(rhs_velx_b)
        deallocate(rhs_dpdx)

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
        allocate(vely(findexl:findexr,findexd:findexu, findexf:findexb))
        vely(:,:,:) = 0

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
        allocate(rhs_vely_b(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_vely_b(vely, rhs_vely_b)
        allocate(resiAyy_b(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, vely, isField)
                    if(isField) then
                        call Resi_vely_b(vely, resiAyy_b)
                        resitemp = resiAyy_b - rhs_vely_b
                        call constructAyy(vely, resitemp, 1)
                    end if
                end do
            end do
        end do
        deallocate(vely)
        deallocate(resiAyy_b)
        deallocate(resitemp)

        do k = 1, localnlays
            do i = 1, localncols
                if((prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                    call index_convert_local_global(myid, 2, i, 1, k, global_ind)
                    do n = 1, AyySize
                        if(AyyRows(n)==global_ind) then
                            AyyValues(n) = AyyStaticValues(n)
                        end if
                    end do
                end if
                if((prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                    call index_convert_local_global(myid, 2, i, localnrows+1, k, global_ind)
                    do n = 1, AyySize
                        if(AyyRows(n)==global_ind) then
                            AyyValues(n) = AyyStaticValues(n)
                        end if
                    end do
                end if
            end do
        end do

        findexl = 1
        findexr = localncols
        if(prow /= 1) then
            findexd = 0
        else
            findexd = 1
        end if
        findexu = localnrows
        findexf = 1
        findexb = localnlays
        allocate(pres(findexl:findexr,findexd:findexu,findexf:findexb))
        pres(:,:,:) = 0
        allocate(rhs_dpdy(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_dpdy(pres, rhs_dpdy)
        allocate(resiAyp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, pres, isField)
                    if(isField) then
                        call Resi_dpdy(pres, resiAyp)
                        resitemp = resiAyp - rhs_dpdy
                        call constructAyp(pres, resitemp)
                    end if
                end do
            end do
        end do
        deallocate(pres)
        deallocate(resiAyp)
        deallocate(resitemp)

        do k = eindexf, eindexb
            do j = eindexd, eindexu
                do i = eindexl, eindexr
                    c = c + 1
                    local_rhs_static(c) = -(rhs_vely_b(i,j,k) + rhs_dpdy(i,j,k) - rhof*gravY)
                    if(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)).or.((prow==pnrows).and.(j==localnrows+1) &!
                        .and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0))) then
                        local_rhs(c) = local_rhs_static(c) - rhof*gravY
                    end if
                end do
            end do
        end do
        deallocate(rhs_vely_b)
        deallocate(rhs_dpdy)

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
        allocate(velz(findexl:findexr,findexd:findexu, findexf:findexb))
        velz(:,:,:) = 0

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
        allocate(rhs_velz_b(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_velz_b(velz, rhs_velz_b)
        allocate(resiAzz_b(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, velz, isField)
                    if(isField) then
                        call Resi_velz_b(velz, resiAzz_b)
                        resitemp = resiAzz_b - rhs_velz_b
                        call constructAzz(velz, resitemp, 1)
                    end if
                end do
            end do
        end do
        deallocate(velz)
        deallocate(resiAzz_b)
        deallocate(resitemp)

        do j = 1, localnrows
            do i = 1, localncols
                if((play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                    call index_convert_local_global(myid, 3, i, j, 1, global_ind)
                    do n = 1, AzzSize
                        if(AzzRows(n)==global_ind) then
                            AzzValues(n) = AzzStaticValues(n)
                        end if
                    end do
                end if
                if((play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                    call index_convert_local_global(myid, 3, i, j, localnlays+1, global_ind)
                    do n = 1, AzzSize
                        if(AzzRows(n)==global_ind) then
                            AzzValues(n) = AzzStaticValues(n)
                        end if
                    end do
                end if
            end do
        end do

        findexl = 1
        findexr = localncols
        findexd = 1
        findexu = localnrows
        if(play /= 1) then
            findexf = 0
        else
            findexf = 1
        end if
        findexb = localnlays
        allocate(pres(findexl:findexr,findexd:findexu,findexf:findexb))
        pres(:,:,:) = 0
        allocate(rhs_dpdz(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_dpdz(pres, rhs_dpdz)
        allocate(resiAzp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, pres, isField)
                    if(isField) then
                        call Resi_dpdz(pres, resiAzp)
                        resitemp = resiAzp - rhs_dpdz
                        call constructAzp(pres, resitemp)
                    end if
                end do
            end do
        end do
        deallocate(pres)
        deallocate(resiAzp)
        deallocate(resitemp)

        do k = eindexf, eindexb
            do j = eindexd, eindexu
                do i = eindexl, eindexr
                    c = c + 1
                    local_rhs_static(c) = -(rhs_velz_b(i,j,k) + rhs_dpdz(i,j,k) - rhof*gravZ)
                    if(((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)).or. &!
                        ((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0))) then
                        local_rhs(c) = local_rhs_static(c) - rhof*gravZ
                    end if
                end do
            end do
        end do
        deallocate(rhs_velz_b)
        deallocate(rhs_dpdz)

        findexl = 1
        findexr = localncols + 1
        findexd = 1
        findexu = localnrows
        findexf = 1
        findexb = localnlays
        allocate(velx(findexl:findexr,findexd:findexu,findexf:findexb))
        velx(:,:,:) = 0
        allocate(rhs_dudx(1:localncols,1:localnrows,1:localnlays))
        call Resi_dudx(velx, rhs_dudx)
        allocate(resiAcx(1:localncols,1:localnrows,1:localnlays))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, velx, isField)
                    if(isField) then
                        call Resi_dudx(velx, resiAcx)
                        call constructAcx(velx, resiAcx)
                    end if
                end do
            end do
        end do
        deallocate(velx)
        deallocate(resiAcx)

        findexl = 1
        findexr = localncols
        findexd = 1
        findexu = localnrows+1
        findexf = 1
        findexb = localnlays
        allocate(vely(findexl:findexr,findexd:findexu,findexf:findexb))
        vely(:,:,:) = 0
        allocate(rhs_dvdy(1:localncols,1:localnrows,1:localnlays))
        call Resi_dvdy(vely, rhs_dvdy)
        allocate(resiAcy(1:localncols,1:localnrows,1:localnlays))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, vely, isField)
                    if(isField) then
                        call Resi_dvdy(vely, resiAcy)
                        call constructAcy(vely, resiAcy)
                    end if
                end do
            end do
        end do
        deallocate(vely)
        deallocate(resiAcy)

        findexl = 1
        findexr = localncols
        findexd = 1
        findexu = localnrows
        findexf = 1
        findexb = localnlays + 1
        allocate(velz(findexl:findexr,findexd:findexu,findexf:findexb))
        velz(:,:,:) = 0
        allocate(rhs_dwdz(1:localncols,1:localnrows,1:localnlays))
        call Resi_dwdz(velz, rhs_dwdz)
        allocate(resiAcz(1:localncols,1:localnrows,1:localnlays))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, velz, isField)
                    if(isField) then
                        call Resi_dwdz(velz, resiAcz)
                        call constructAcz(velz, resiAcz)
                    end if
                end do
            end do
        end do
        deallocate(velz)
        deallocate(resiAcz)

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_rhs_static(c) = -(rhs_dudx(i,j,k) + rhs_dvdy(i,j,k) + rhs_dwdz(i,j,k))
                end do
            end do
        end do
        deallocate(rhs_dudx)
        deallocate(rhs_dvdy)
        deallocate(rhs_dwdz)

    end subroutine genStaticPara_vp

    subroutine computekc()

        real(kind=8) :: Sc, vmodulus, radius, Rep
        integer :: i, j, k

        Sc = visc/rhof*dm

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    vmodulus = dsqrt((vx(i+1,j,k)-vx(i,j,k))**2.D0+(vy(i,j+1,k)-vy(i,j,k))**2.D0+ &!
                        (vz(i,j,k+1)-vz(i,j,k))**2.D0)

                    radius = radiusInit*(1.D0-poroInit(xlower+i-1,ylower+j-1,zlower+k-1))/ &!
                        (poroInit(xlower+i-1,ylower+j-1,zlower+k-1)*(1.D0-poro(i,j,k)))*poro(i,j,k)

                    Rep = 2.D0*vmodulus*radius/(visc/rhof)

                    kc(i,j,k) = (ShInfinity+7.D-1*Rep**(1.D0/2.D0)*Sc**(1.D0/3.D0))*dm/2.D0/radius

                end do
            end do
        end do

    end subroutine computekc

    subroutine computePoro()

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu, requestf, requestb, requestlu, &!
            requestlb, requestrd, requestru, requestrf, requestrb, requestdb, requestuf, requestub
        real(kind=8), dimension(:), allocatable :: sent, recv
        real(kind=8) :: coe
        integer :: sentSize, recvSize
        integer :: i, j, k, c, ierr

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    coe = avInit(xlower+i-1,ylower+j-1,zlower+k-1)*al*Cf(i,j,k)*kc(i,j,k)*ks*(ts(t)-ts(t-1))/ &!
                        (rhos*(kc(i,j,k)+ks)*(1.D0-poroInit(xlower+i-1,ylower+j-1,zlower+k-1)))
                    poro_old(i,j,k) = poro(i,j,k)
                    poro(i,j,k) = (coe+poro_old(i,j,k))/(1+coe)
                end do
            end do
        end do

        ! send
        ! the 6 faces
        if(pcol /= 1) then
            sentSize = 2*localnrows*localnlays
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = poro(1,j,k)
                    c = c + 1
                    sent(c) = poro_old(1,j,k)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            sentSize = 4*localnrows*localnlays
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = localncols-1, localncols
                        c = c + 1
                        sent(c) = poro(i,j,k)
                        c = c + 1
                        sent(c) = poro_old(i,j,k)
                    end do
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            sentSize = localncols*2*localnlays
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    sent(c) = poro(i,1,k)
                    c = c + 1
                    sent(c) = poro_old(i,1,k)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = localncols*4*localnlays
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do j = localnrows-1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        sent(c) = poro(i,j,k)
                        c = c + 1
                        sent(c) = poro_old(i,j,k)
                    end do
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if(play /= 1) then
            sentSize = localncols*localnrows*2
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    sent(c) = poro(i,j,1)
                    c = c + 1
                    sent(c) = poro_old(i,j,1)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols*pnrows, myid, MPI_COMM_WORLD, requestf, ierr)
            deallocate(sent)
        end if

        if(play /= pnlays) then
            sentSize = localncols*localnrows*4
            allocate(sent(sentSize))
            c = 0
            do k = localnlays-1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        sent(c) = poro(i,j,k)
                        c = c + 1
                        sent(c) = poro_old(i,j,k)
                    end do
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols*pnrows, myid, MPI_COMM_WORLD, requestb, ierr)
            deallocate(sent)
        end if

        ! the 9 edges
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

        ! receive
        ! the 6 faces
        if(pcol /= pncols) then
            recvSize = 2*localnrows*localnlays
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    poro(localncols+1,j,k) = recv(c)
                    c = c + 1
                    poro_old(localncols+1,j,k) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            recvSize = 4*localnrows*localnlays
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = -1, 0
                        c = c + 1
                        poro(i,j,k) = recv(c)
                        c = c + 1
                        poro_old(i,j,k) = recv(c)
                    end do
                end do
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            recvSize = localncols*2*localnlays
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    poro(i,localnrows+1,k) = recv(c)
                    c = c + 1
                    poro_old(i,localnrows+1,k) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = localncols*4*localnlays
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do j = -1, 0
                    do i = 1, localncols
                        c = c + 1
                        poro(i,j,k) = recv(c)
                        c = c + 1
                        poro_old(i,j,k) = recv(c)
                    end do
                end do
            end do
            deallocate(recv)
        end if

        if(play /= pnlays) then
            recvSize = localncols*localnrows*2
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols*pnrows, myid+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    poro(i,j,localnlays+1) = recv(c)
                    c = c + 1
                    poro_old(i,j,localnlays+1) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(play /= 1) then
            recvSize = localncols*localnrows*4
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols*pnrows, myid-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = -1, 0
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        poro(i,j,k) = recv(c)
                        c = c + 1
                        poro_old(i,j,k) = recv(c)
                    end do
                end do
            end do
            deallocate(recv)
        end if

        ! the 9 edges
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

        ! wait
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

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestr, requestu, requestb
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, k, c, ierr

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

        ! send
        ! the 3 faces
        if(pcol /= pncols) then
            sentSize = localnrows*localnlays*3
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = Kxx(localncols,j,k)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = Kyy(localncols,j,k)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    sent(c) = Kzz(localncols,j,k)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = localncols*localnlays*3
            allocate(sent(sentSize))
            c = 0
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    sent(c) = Kxx(i,localnrows,k)
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    sent(c) = Kyy(i,localnrows,k)
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    sent(c) = Kzz(i,localnrows,k)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if(play /= pnlays) then
            sentSize = localncols*localnrows*3
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    sent(c) = Kxx(i,j,localnlays)
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    sent(c) = Kyy(i,j,localnlays)
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    sent(c) = Kzz(i,j,localnlays)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, &!
                myid+pncols*pnrows, myid, MPI_COMM_WORLD, requestb, ierr)
            deallocate(sent)
        end if

        ! receive
        ! the 3 faces
        if(pcol /= 1) then
            recvSize = localnrows*localnlays*3
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    Kxx(0,j,k) = recv(c)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    Kyy(0,j,k) = recv(c)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    c = c + 1
                    Kzz(0,j,k) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = localncols*localnlays*3
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    Kxx(i,0,k) = recv(c)
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    Kyy(i,0,k) = recv(c)
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    c = c + 1
                    Kzz(i,0,k) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(play /= 1) then
            recvSize = localncols*localnrows*3
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, &!
                myid-pncols*pnrows, myid-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    Kxx(i,j,0) = recv(c)
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    Kyy(i,j,0) = recv(c)
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    Kzz(i,j,0) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        ! wait
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if(play /= pnlays) then
            call MPI_WAIT(requestb, status, ierr)
        end if

    end subroutine computeK

    subroutine computeKHarm()

        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: i, j, k

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
        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if
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
        if(play /= pnlays) then
            indexb = localnlays
        else
            indexb = localnlays + 1
        end if

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
                    av(i,j,k) = avInit(xlower+i-1,ylower+j-1,zlower+k-1)*poro(i,j,k)/ &!
                        poroInit(xlower+i-1,ylower+j-1,zlower+k-1)*dsqrt(1.D0*KxxInit(xlower+i-1,ylower+j-1,zlower+k-1) &!
                        *poro(i,j,k)/Kxx(i,j,k)/poroInit(xlower+i-1,ylower+j-1,zlower+k-1))
                end do
            end do
        end do

    end subroutine computeav

    ! Generate the values on the right-hand side and the coefficients of the matrix A,
    ! and the subroutine will generate the values that will change with the time iteration.
    subroutine genDynPara_vp()

        integer :: findexl, findexr, findexd, findexu, findexf, findexb
        integer :: eindexl, eindexr, eindexd, eindexu, eindexf, eindexb
        integer :: global_ind
        integer, dimension(:,:,:), pointer :: velx, vely, velz, pres
        logical :: isField
        real(kind=8), dimension(:,:,:), pointer :: rhs_velx, rhs_vely, rhs_velz, rhs_dpdt
        real(kind=8), dimension(:,:,:), pointer :: resiAxx, resiAyy, resiAzz, resiAcp, resitemp
        integer :: AxxBeInd, AyyBeInd, AzzBeInd
        integer :: i, j, k, n, c

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
        allocate(velx(findexl:findexr,findexd:findexu,findexf:findexb))
        velx(:,:,:) = 0

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
        allocate(rhs_velx(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_velx(velx, rhs_velx)
        allocate(resiAxx(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, velx, isField)
                    if(isField) then
                        call Resi_velx(velx, resiAxx)
                        ! generate the dynamic coefficients of Axx. The scheme is just like the DNA copy according
                        ! to a template, here the AxxCols and the AxxRows are just the template.
                        resitemp = resiAxx - rhs_velx
                        call constructAxx(velx, resitemp, 2)
                    end if
                end do
            end do
        end do
        deallocate(resiAxx)
        deallocate(resitemp)
        deallocate(velx)

        AxxBeInd = 1
        do k = eindexf, eindexb
            do j = eindexd, eindexu
                do i = eindexl, eindexr
                    if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0))).and.(.not.((pcol==pncols) &!
                        .and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)))) then
                        call index_convert_local_global(myid, 1, i, j, k, global_ind)
                        do n = AxxBeInd, AxxSize
                            if(AxxRows(n)==global_ind) then
                                AxxValues(n) = AxxStaticValues(n) + AxxDynValues(n)
                            else
                                AxxBeInd = n
                                exit
                            end if
                        end do
                    else
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
                        local_rhs(c) = local_rhs_static(c) - rhs_velx(i,j,k)
                    end if
                end do
            end do
        end do
        deallocate(rhs_velx)

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
        allocate(vely(findexl:findexr,findexd:findexu,findexf:findexb))
        vely(:,:,:) = 0

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
        allocate(rhs_vely(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_vely(vely, rhs_vely)
        allocate(resiAyy(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, vely, isField)
                    if(isField) then
                        call Resi_vely(vely, resiAyy)
                        resitemp = resiAyy - rhs_vely
                        call constructAyy(vely, resitemp, 2)
                    end if
                end do
            end do
        end do
        deallocate(resiAyy)
        deallocate(resitemp)
        deallocate(vely)

        AyyBeInd = 1
        do k = eindexf, eindexb
            do j = eindexd, eindexu
                do i = eindexl, eindexr
                    if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0))).and.(.not.((prow==pnrows).and. &!
                        (j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)))) then
                        call index_convert_local_global(myid, 2, i, j, k, global_ind)
                        do n = AyyBeInd, AyySize
                            if(AyyRows(n)==global_ind) then
                                AyyValues(n) = AyyStaticValues(n) + AyyDynValues(n)
                            else
                                AyyBeInd = n
                                exit
                            end if
                        end do
                    else
                        AyyBeInd = AyyBeInd + 1
                    end if
                end do
            end do
        end do

        do k = eindexf, eindexb
            do j = eindexd, eindexu
                do i = eindexl, eindexr
                    c = c + 1
                    if(.not.(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)).or.((prow==pnrows).and. &!
                        (j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)))) then
                        local_rhs(c) = local_rhs_static(c) - rhs_vely(i,j,k)
                    end if
                end do
            end do
        end do
        deallocate(rhs_vely)

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
        allocate(velz(findexl:findexr,findexd:findexu,findexf:findexb))
        velz(:,:,:) = 0

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
        allocate(rhs_velz(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        call Resi_velz(velz, rhs_velz)
        allocate(resiAzz(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, velz, isField)
                    if(isField) then
                        call Resi_velz(velz, resiAzz)
                        resitemp = resiAzz - rhs_velz
                        call constructAzz(velz, resitemp, 2)
                    end if
                end do
            end do
        end do
        deallocate(resiAzz)
        deallocate(resitemp)
        deallocate(velz)

        AzzBeInd = 1
        do k = eindexf, eindexb
            do j = eindexd, eindexu
                do i = eindexl, eindexr
                    if((.not.((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0))).and.(.not.((play==pnlays).and. &!
                        (k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)))) then
                        call index_convert_local_global(myid, 3, i, j, k, global_ind)
                        do n = AzzBeInd, AzzSize
                            if(AzzRows(n)==global_ind) then
                                AzzValues(n) = AzzStaticValues(n) + AzzDynValues(n)
                            else
                                AzzBeInd = n
                                exit
                            end if
                        end do
                    else
                        AzzBeInd = AzzBeInd + 1
                    end if
                end do
            end do
        end do

        do k = eindexf, eindexb
            do j = eindexd, eindexu
                do i = eindexl, eindexr
                    c = c + 1
                    if(.not.(((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)).or.((play==pnlays).and. &!
                        (k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)))) then
                        local_rhs(c) = local_rhs_static(c) - rhs_velz(i,j,k)
                    end if
                end do
            end do
        end do
        deallocate(rhs_velz)

        allocate(pres(1:localncols,1:localnrows,1:localnlays))
        pres(:,:,:) = 0
        allocate(rhs_dpdt(1:localncols,1:localnrows,1:localnlays))
        call Resi_dpdt(pres, rhs_dpdt)
        allocate(resiAcp(1:localncols,1:localnrows,1:localnlays))
        allocate(resitemp(1:localncols,1:localnrows,1:localnlays))
        do k = 1, 3
            do j = 1, 3
                do i = 1, 3
                    call genExpField(i, j, k, localncols, localnrows, localnlays, pres, isField)
                    if(isField) then
                        call Resi_dpdt(pres, resiAcp)
                        resitemp = resiAcp - rhs_dpdt
                        call constructAcp(pres, resitemp)
                    end if
                end do
            end do
        end do
        deallocate(resiAcp)
        deallocate(resitemp)
        deallocate(pres)

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_rhs(c) = local_rhs_static(c) - rhs_dpdt(i,j,k) - &!
                        kc(i,j,k)*ks*Cf(i,j,k)*av(i,j,k)*al/rhos/(kc(i,j,k)+ks)
                end do
            end do
        end do
        deallocate(rhs_dpdt)

    end subroutine genDynPara_vp

    subroutine computevp()
       
        integer :: indexr, indexu, indexb
        integer :: nVelx, nVely, nVelz, nPres
        integer :: AxxBeInd, AxpBeInd, AyyBeInd, AypBeInd, AzzBeInd, AzpBeInd, AcxBeInd, AcyBeInd, AczBeInd, AcpBeInd
        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), allocatable :: local_x
        integer :: status(MPI_STATUS_SIZE)
        integer :: request, requestl, requestr, requestd, requestu, requestf, requestb, requestld, requestlu, &!
            requestlf, requestlb, requestrd, requestrf, requestdf, requestdb, requestuf
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        real(kind=8) :: solvertimestart, solvertimefinish
        integer :: i, j, k, l, n, c, num_iter, ierr

        allocate(cols(9))
        allocate(values(9)) ! Each row of A has at most 9 nonzero values.
        allocate(local_x(local_x_size))

        cols(:) = 0
        values(:) = 0.D0

#ifdef LAPACK

        A_lapack(:,:) = 0.D0
        b_lapack(:) = local_rhs(:)
        IPIV(:) = 0.D0

#elif defined(UMFPACK)

        Ap(1) = 0
        Ai(:) = 0
        Ax(:) = 0.D0

#elif defined(MUMPS)

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
                allocate(slave_data(slave_vp_data_size(i)))
                call MPI_RECV(slave_data, slave_vp_data_size(i), MPI_DOUBLE_PRECISION, &!
                    i, i, MPI_COMM_WORLD, status, ierr)
                mumps_par%RHS(c:c+slave_vp_data_size(i)-1) = slave_data(:)
                c = c + slave_vp_data_size(i)
                deallocate(slave_data)
            end do
        end if

        if(mumps_par%MYID /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

#elif defined(HYPRE)

        call HYPRE_IJVectorSetValues(b, local_x_size, rows, local_rhs, ierr)
        call HYPRE_IJVectorSetValues(x, local_x_size, rows, initial_x_guess, ierr)

        call HYPRE_IJVectorAssemble(b, ierr)
        call HYPRE_IJVectorAssemble(x, ierr)

        call HYPRE_IJVectorGetObject(b, par_b, ierr)
        call HYPRE_IJVectorGetObject(x, par_x, ierr)

#endif

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
        if(play /= pnlays) then
            nVelz = localncols*localnrows*localnlays
        else
            nVelz = localncols*localnrows*(localnlays+1)
        end if
        nPres = localncols*localnrows*localnlays

        AxxBeInd = 1
        AxpBeInd = 1
        AyyBeInd = 1
        AypBeInd = 1
        AzzBeInd = 1
        AzpBeInd = 1
        AcxBeInd = 1
        AcyBeInd = 1
        AczBeInd = 1
        AcpBeInd = 1

        do n = ilower, iupper

            cols(:) = 0
            values(:) = 0
            c = 0

            ! the line is in the x-momentum part
            if(n <= ilower+nVelx-1) then

                do l = AxxBeInd, AxxSize
                    if(AxxRows(l) == n) then
                        c = c + 1
                        cols(c) = AxxCols(l)
                        values(c) = AxxValues(l)
                    else
                        AxxBeInd = l
                        exit
                    end if
                end do
                do l = AxpBeInd, AxpSize
                    if(AxpRows(l) == n) then
                        c = c + 1
                        cols(c) = AxpCols(l)
                        values(c) = AxpValues(l)
                    else
                        AxpBeInd = l
                        exit
                    end if
                end do

            ! the line is in the y-momentum part
            elseif((n >= ilower+nVelx).and.(n <= ilower+nVelx+nVely-1)) then

                do l = AyyBeInd, AyySize
                    if(AyyRows(l) == n) then
                        c = c + 1
                        cols(c) = AyyCols(l)
                        values(c) = AyyValues(l)
                    else
                        AyyBeInd = l
                        exit
                    end if
                end do
                do l = AypBeInd, AypSize
                    if(AypRows(l) == n) then
                        c = c + 1
                        cols(c) = AypCols(l)
                        values(c) = AypValues(l)
                    else
                        AypBeInd = l
                        exit
                    end if
                end do

            ! the line is in the z-momentum part
            elseif((n >= ilower+nVelx+nVely).and.(n <= ilower+nVelx+nVely+nVelz-1)) then

                do l = AzzBeInd, AzzSize
                    if(AzzRows(l) == n) then
                        c = c + 1
                        cols(c) = AzzCols(l)
                        values(c) = AzzValues(l)
                    else
                        AzzBeInd = l
                        exit
                    end if
                end do
                do l = AzpBeInd, AzpSize
                    if(AzpRows(l) == n) then
                        c = c + 1
                        cols(c) = AzpCols(l)
                        values(c) = AzpValues(l)
                    else
                        AzpBeInd = l
                        exit
                    end if
                end do

            ! the line is in the continuity part
            elseif(n >= ilower+nVelx+nVely+nVelz) then

                do l = AcxBeInd, AcxSize
                    if(AcxRows(l) == n) then
                        c = c + 1
                        cols(c) = AcxCols(l)
                        values(c) = AcxValues(l)
                    else
                        AcxBeInd = l
                        exit
                    end if
                end do
                do l = AcyBeInd, AcySize
                    if(AcyRows(l) == n) then
                        c = c + 1
                        cols(c) = AcyCols(l)
                        values(c) = AcyValues(l)
                    else
                        AcyBeInd = l
                        exit
                    end if
                end do
                do l = AczBeInd, AczSize
                    if(AczRows(l) == n) then
                        c = c + 1
                        cols(c) = AczCols(l)
                        values(c) = AczValues(l)
                    else
                        AczBeInd = l
                        exit
                    end if
                end do
                do l = AcpBeInd, AcpSize
                    if(AcpRows(l) == n) then
                        c = c + 1
                        cols(c) = AcpCols(l)
                        values(c) = AcpValues(l)
                    else
                        AcpBeInd = l
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

        ! you have to make sure that the number of processors is set to 1 when using such method.
        solvertimestart = MPI_Wtime()
        call dgesv(local_x_size, 1, A_lapack, local_x_size, IPIV, b_lapack, local_x_size, LAPACKINFO)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart

        local_x(:) = b_lapack(:)
        if(LAPACKINFO /= 0) then
            print *, 'LAPACK solver error. INFO = ', LAPACKINFO
            stop
        end if

#elif defined(UMFPACK)

        solvertimestart = MPI_Wtime()
        call umf4def(control)
        call umf4sym(local_x_size, local_x_size, Ap, Ai, Ax, symbolic, control, umfinfo)
        call umf4num(Ap, Ai, Ax, symbolic, numeric, control, umfinfo)
        call umf4fsym(symbolic)
        call umf4sol(1, local_x, local_rhs, numeric, control, umfinfo) ! 1 means A'x=b
        if(umfinfo(1) .lt. 0) then
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

        mumps_par%ICNTL(14) = 100
        mumps_par%JOB = 6
        solvertimestart = MPI_Wtime()
        call DMUMPS(mumps_par)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(mumps_par%INFOG(1) < 0) then
            stop
        end if

        if(mumps_par%MYID == 0) then
            local_x(:) = mumps_par%RHS(1:local_x_size)
            allocate(requestarray(nProcs-1))
            c = local_x_size + 1
            do i = 1, nProcs-1
                allocate(slave_data(slave_vp_data_size(i)))
                slave_data(:) = mumps_par%RHS(c:c+slave_vp_data_size(i)-1)
                call MPI_IBSEND(slave_data, slave_vp_data_size(i), MPI_DOUBLE_PRECISION, i, myid, &!
                    MPI_COMM_WORLD, requestarray(i), ierr)
                c = c + slave_vp_data_size(i)
                deallocate(slave_data)
            end do
        end if

        if(mumps_par%MYID /= 0) then
            allocate(slave_data(local_x_size))
            call MPI_RECV(slave_data, local_x_size, MPI_DOUBLE_PRECISION, &!
                0, 0, MPI_COMM_WORLD, status, ierr)
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
        call HYPRE_ParCSRGMRESGetNumIteratio(solver, num_iter, ierr)

        if(ierr /= 0) then
            if(myid == 0) then
                print *, 'HYPRE solver error. ierr = ', ierr
                stop
            end if
        end if

        call HYPRE_IJVectorGetValues(x, local_x_size, rows, local_x, ierr)

        ! let the solution of this time step be the initial x guess in the next time step.
        ! by this way, the number of solver iteration steps can be reduced greatly.
        initial_x_guess(:) = local_x(:)

#endif

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

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    p(i,j,k) = local_x(c)
                end do
            end do
        end do

        ! communicate, send
        ! the 6 faces
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

        ! communicate, receive
        ! the 6 faces
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

        ! wait for the completion of the communication
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

        deallocate(cols)
        deallocate(values)
        deallocate(local_x)

    end subroutine computevp

    subroutine genDynPara_Cf()

        integer :: findexl, findexr, findexd, findexu, findexf, findexb
        integer :: eindexl, eindexr, eindexd, eindexu, eindexf, eindexb
        integer, dimension(:,:,:), pointer :: conc
        logical :: isField
        real(kind=8), dimension(:,:,:), pointer :: rhs_Cf
        real(kind=8), dimension(:,:,:), pointer :: resiAcf, resitemp
        integer :: i, j, k, c

        findexl = 0
        findexr = localncols + 1
        findexd = 0
        findexu = localnrows + 1
        findexf = 0
        findexb = localnlays + 1
        allocate(conc(findexl:findexr,findexd:findexu,findexf:findexb))
        conc(:,:,:) = 0

        eindexl = 1
        eindexr = localncols
        eindexd = 1
        eindexu = localnrows
        eindexf = 1
        eindexb = localnlays
        allocate(rhs_Cf(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        rhs_Cf(:,:,:) = 0.D0

        call Resi_Cf(conc, rhs_Cf)
        allocate(resiAcf(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu,eindexf:eindexb))
        do k = findexf, findexf+2
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, k, findexr, findexu, findexb, conc, isField)
                    if(isField) then
                        call Resi_Cf(conc, resiAcf)
                        resitemp = resiAcf - rhs_Cf
                        call constructAcf(conc, resitemp)
                    end if
                end do
            end do
        end do

        c = 0
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    local_rhs_Cf(c) = -rhs_Cf(i,j,k)
                end do
            end do
        end do

        deallocate(conc)
        deallocate(rhs_Cf)
        deallocate(resiAcf)
        deallocate(resitemp)

    end subroutine genDynPara_Cf

    subroutine computeCf()

        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), pointer :: local_x
        integer :: AcfBeInd
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data
        real(kind=8) :: solvertimestart, solvertimefinish
        integer :: i, j, k, l, n, c, num_iter, ierr

        allocate(cols(19))
        allocate(values(19)) ! Each row of A has at most 19 nonzero values.
        allocate(local_x(local_x_size_Cf))

#ifdef LAPACK

        A_lapack_Cf(:,:) = 0.D0
        b_lapack_Cf(:) = local_rhs_Cf(:)
        IPIV_Cf(:) = 0.D0

#elif defined(UMFPACK)

        Ap_Cf(1) = 0
        Ai_Cf(:) = 0
        Ax_Cf(:) = 0.D0

#elif defined(MUMPS)

        mumps_NNZ_loc_Cf = 0
        mumps_IRN_loc_Cf(:) = 0
        mumps_JCN_loc_Cf(:) = 0
        mumps_A_loc_Cf(:) = 0.D0

        if(mumps_par_Cf%MYID /= 0) then
            call MPI_IBSEND(local_rhs_Cf, local_x_size_Cf, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(mumps_par_Cf%MYID == 0) then
            mumps_par_Cf%RHS(1:local_x_size_Cf) = local_rhs_Cf(1:local_x_size_Cf)
            c = local_x_size_Cf + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Cf))
                call MPI_RECV(slave_data, local_x_size_Cf, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                mumps_par_Cf%RHS(c:c+local_x_size_Cf-1) = slave_data(:)
                c = c + local_x_size_Cf
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Cf%MYID /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

#elif defined(HYPRE)

        call HYPRE_IJVectorSetValues(b_Cf, local_x_size_Cf, rows_Cf, local_rhs_Cf, ierr)
        call HYPRE_IJVectorSetValues(x_Cf, local_x_size_Cf, rows_Cf, initial_x_guess_Cf, ierr)

        call HYPRE_IJVectorAssemble(b_Cf, ierr)
        call HYPRE_IJVectorAssemble(x_Cf, ierr)

        call HYPRE_IJVectorGetObject(b_Cf, par_b_Cf, ierr)
        call HYPRE_IJVectorGetObject(x_Cf, par_x_Cf, ierr)

#endif

        AcfBeInd = 1
        do n = ilower_Cf, iupper_Cf

            cols(:) = 0
            values(:) = 0.D0

            c = 0
            do l = AcfBeInd, AcfSize
                if(AcfRows(l) == n) then
                    c = c + 1
                    cols(c) = AcfCols(l)
                    values(c) = AcfValues(l)
                else
                    AcfBeInd = l
                    exit
                end if
            end do

#ifdef LAPACK

            do l = 1, c
                A_lapack_Cf(n,cols(l)) = values(l)
            end do

#elif defined(UMFPACK)

            Ap_Cf(n+1) = Ap_Cf(n) + c
            Ai_Cf(Ap_Cf(n)+1:Ap_Cf(n+1)) = cols(1:c) - 1
            Ax_Cf(Ap_Cf(n)+1:Ap_Cf(n+1)) = values(1:c)

#elif defined(MUMPS)

            mumps_IRN_loc_Cf(mumps_NNZ_loc_Cf+1:mumps_NNZ_loc_Cf+c) = n
            mumps_JCN_loc_Cf(mumps_NNZ_loc_Cf+1:mumps_NNZ_loc_Cf+c) = cols(1:c)
            mumps_A_loc_Cf(mumps_NNZ_loc_Cf+1:mumps_NNZ_loc_Cf+c) = values(1:c)
            mumps_NNZ_loc_Cf = mumps_NNZ_loc_Cf + c

#elif defined(HYPRE)

            call HYPRE_IJMatrixSetValues(A_Cf, 1, c, n, cols, values, ierr)

#endif

        end do

#ifdef LAPACK

        ! you have to make sure that the number of processors is set to 1 when using such method.
        solvertimestart = MPI_Wtime()
        call dgesv(local_x_size_Cf, 1, A_lapack_Cf, local_x_size_Cf, IPIV_Cf, b_lapack_Cf, local_x_size_Cf, LAPACKINFO)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        local_x(:) = b_lapack_Cf(:)
        if(LAPACKINFO /= 0) then
            print *, 'LAPACK solver error when solving Cf. INFO = ', LAPACKINFO
            stop
        end if

#elif defined(UMFPACK)

        solvertimestart = MPI_Wtime()
        call umf4def(control)
        call umf4sym(local_x_size_Cf, local_x_size_Cf, Ap_Cf, Ai_Cf, Ax_Cf, symbolic, control, umfinfo)
        call umf4num(Ap_Cf, Ai_Cf, Ax_Cf, symbolic, numeric, control, umfinfo)
        call umf4fsym(symbolic)
        call umf4sol(1, local_x, local_rhs_Cf, numeric, control, umfinfo) ! 1 means A'x=b
        if(umfinfo(1) < 0) then
            print *, 'UMFPACK solver error when solving Cf. Info: ', umfinfo(1)
            stop
        end if
        call umf4fnum(numeric)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart

#elif defined(MUMPS)

        mumps_par_Cf%NNZ_loc = mumps_NNZ_loc_Cf
        allocate(mumps_par_Cf%IRN_loc(mumps_par_Cf%NNZ_loc))
        allocate(mumps_par_Cf%JCN_loc(mumps_par_Cf%NNZ_loc))
        allocate(mumps_par_Cf%A_loc(mumps_par_Cf%NNZ_loc))
        mumps_par_Cf%IRN_loc(1:mumps_par_Cf%NNZ_loc) = mumps_IRN_loc_Cf(1:mumps_NNZ_loc_Cf)
        mumps_par_Cf%JCN_loc(1:mumps_par_Cf%NNZ_loc) = mumps_JCN_loc_Cf(1:mumps_NNZ_loc_Cf)
        mumps_par_Cf%A_loc(1:mumps_par_Cf%NNZ_loc) = mumps_A_loc_Cf(1:mumps_NNZ_loc_Cf)

        mumps_par_Cf%ICNTL(14) = 100
        mumps_par_Cf%JOB = 6
        solvertimestart = MPI_Wtime()
        call DMUMPS(mumps_par_Cf)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(mumps_par_Cf%INFOG(1) < 0) then
            stop
        end if

        if(mumps_par_Cf%MYID == 0) then
            local_x(:) = mumps_par_Cf%RHS(1:local_x_size_Cf)
            allocate(requestarray(nProcs-1))
            c = local_x_size_Cf + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Cf))
                slave_data(:) = mumps_par_Cf%RHS(c:c+local_x_size_Cf-1)
                call MPI_IBSEND(slave_data, local_x_size_Cf, MPI_DOUBLE_PRECISION, i, myid, &!
                    MPI_COMM_WORLD, requestarray(i), ierr)
                c = c + local_x_size_Cf
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Cf%MYID /= 0) then
            allocate(slave_data(local_x_size_Cf))
            call MPI_RECV(slave_data, local_x_size_Cf, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
            local_x(:) = slave_data(:)
            deallocate(slave_data)
        end if

        if(mumps_par_Cf%MYID == 0) then
            do i = 1, nProcs-1
                call MPI_WAIT(requestarray(i), status, ierr)
            end do
            deallocate(requestarray)
        end if

        deallocate(mumps_par_Cf%IRN_loc)
        deallocate(mumps_par_Cf%JCN_loc)
        deallocate(mumps_par_Cf%A_loc)

#elif defined(HYPRE)

        call HYPRE_IJMatrixAssemble(A_Cf, ierr)
        call HYPRE_IJMatrixGetObject(A_Cf, parcsr_A_Cf, ierr)

        solvertimestart = MPI_Wtime()
        call HYPRE_ParCSRGMRESSetup(solver_Cf, parcsr_A_Cf, par_b_Cf, par_x_Cf, ierr)
        call HYPRE_ParCSRGMRESSolve(solver_Cf, parcsr_A_Cf, par_b_Cf, par_x_Cf, ierr)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        call HYPRE_ParCSRGMRESGetNumIteratio(solver_Cf, num_iter, ierr)

        if(ierr /= 0) then
            if(myid == 0) then
                print *, 'HYPRE solver error when solving Cf. ierr = ', ierr
                stop
            end if
        end if

        call HYPRE_IJVectorGetValues(x_Cf, local_x_size_Cf, rows_Cf, local_x, ierr)

        ! let the solution of this time step be the initial x guess in the next time step.
        ! by this way, the number of solver iteration steps can be reduced greatly.
        initial_x_guess_Cf(:) = local_x(:)

#endif

        c = 0
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    Cf(i,j,k) = local_x(c)
                end do
            end do
        end do

        deallocate(cols)
        deallocate(values)
        deallocate(local_x)

    end subroutine computeCf

    subroutine outputHisData()

        real(kind=8), dimension(:), allocatable :: local_data, recv
        real(kind=8), dimension(:), pointer :: global_data
        integer :: local_data_size
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        real(kind=8) :: poroavg, Kxxavg, avavg, Cfavg, qavg, pavg
        real(kind=8) :: localqsum, localpsum
        integer :: i, j, k, c, ierr

        local_data_size = 5
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
                end do
            end do
        end do

        ! output the average values of poro, Kxx, av, p, Cf
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

            write(45, fmt='(es24.16)', iostat=ierr) qavg
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

            write(46, fmt='(es24.16)', iostat=ierr) pavg
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

        local_data_size = local_x_size+3*localncols*localnrows*localnlays
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

        if(myid /= 0) then

            call MPI_IBSEND(local_data, local_data_size, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(myid == 0) then

            allocate(global_data((nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1)+4*nx*ny*nz))

            global_data(1:local_data_size) = local_data(1:local_data_size)
            c = local_data_size + 1
            do i = 1, nProcs-1
                slave_data_size = slave_vp_data_size(i)+3*localncols*localnrows*localnlays
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

        deallocate(hx)
        deallocate(hy)
        deallocate(hz)
        deallocate(kc)
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

        deallocate(local_rhs)
        deallocate(local_rhs_static)
        deallocate(local_rhs_Cf)
        deallocate(AxxCols)
        deallocate(AxxRows)
        deallocate(AxxStaticValues)
        deallocate(AxxDynValues)
        deallocate(AxxValues)
        deallocate(AxpCols)
        deallocate(AxpRows)
        deallocate(AxpValues)
        deallocate(AyyCols)
        deallocate(AyyRows)
        deallocate(AyyStaticValues)
        deallocate(AyyDynValues)
        deallocate(AyyValues)
        deallocate(AypCols)
        deallocate(AypRows)
        deallocate(AypValues)
        deallocate(AzzCols)
        deallocate(AzzRows)
        deallocate(AzzStaticValues)
        deallocate(AzzDynValues)
        deallocate(AzzValues)
        deallocate(AzpCols)
        deallocate(AzpRows)
        deallocate(AzpValues)
        deallocate(AcxCols)
        deallocate(AcxRows)
        deallocate(AcxValues)
        deallocate(AcyCols)
        deallocate(AcyRows)
        deallocate(AcyValues)
        deallocate(AczCols)
        deallocate(AczRows)
        deallocate(AczValues)
        deallocate(AcpCols)
        deallocate(AcpRows)
        deallocate(AcpValues)
        deallocate(AcfCols)
        deallocate(AcfRows)
        deallocate(AcfValues)
        deallocate(AxxEntryNum)
        deallocate(AxpEntryNum)
        deallocate(AyyEntryNum)
        deallocate(AypEntryNum)
        deallocate(AzzEntryNum)
        deallocate(AzpEntryNum)
        deallocate(AcxEntryNum)
        deallocate(AcyEntryNum)
        deallocate(AczEntryNum)
        deallocate(AcpEntryNum)
        deallocate(AcfEntryNum)
        deallocate(AxxEntryBase)
        deallocate(AxpEntryBase)
        deallocate(AyyEntryBase)
        deallocate(AypEntryBase)
        deallocate(AzzEntryBase)
        deallocate(AzpEntryBase)
        deallocate(AcxEntryBase)
        deallocate(AcyEntryBase)
        deallocate(AczEntryBase)
        deallocate(AcpEntryBase)
        deallocate(AcfEntryBase)

        if((nProcs>1).and.(myid == 0)) then
            deallocate(slave_vp_data_size)
        end if

        if(myid == 0) then
            close(40)
            close(41)
            close(42)
            close(43)
            close(44)
            close(45)
            close(46)
        end if

#ifdef LAPACK

        deallocate(A_lapack)
        deallocate(b_lapack)
        deallocate(IPIV)
        deallocate(A_lapack_Cf)
        deallocate(b_lapack_Cf)
        deallocate(IPIV_Cf)

#elif defined(UMFPACK)

        deallocate(Ap)
        deallocate(Ai)
        deallocate(Ax)
        deallocate(Ap_Cf)
        deallocate(Ai_Cf)
        deallocate(Ax_Cf)

#elif defined(MUMPS)

        deallocate(mumps_par%RHS)
        deallocate(mumps_par_Cf%RHS)

        mumps_par%JOB = -2
        call DMUMPS(mumps_par)
        mumps_par_Cf%JOB = -2
        call DMUMPS(mumps_par_Cf)

        deallocate(mumps_IRN_loc)
        deallocate(mumps_JCN_loc)
        deallocate(mumps_A_loc)
        deallocate(mumps_IRN_loc_Cf)
        deallocate(mumps_JCN_loc_Cf)
        deallocate(mumps_A_loc_Cf)

#elif defined(HYPRE)

        call HYPRE_ParaSailsDestroy(precond, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver, ierr)
        call HYPRE_ParaSailsDestroy(precond_Cf, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver_Cf, ierr)

        call HYPRE_IJMatrixDestroy(A, ierr)
        call HYPRE_IJVectorDestroy(b, ierr)
        call HYPRE_IJVectorDestroy(x, ierr)
        call HYPRE_IJMatrixDestroy(A_Cf, ierr)
        call HYPRE_IJVectorDestroy(b_Cf, ierr)
        call HYPRE_IJVectorDestroy(x_Cf, ierr)

        deallocate(rows)
        deallocate(initial_x_guess)
        deallocate(rows_Cf)
        deallocate(initial_x_guess_Cf)

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

        call genStaticPara_vp()

        call computekc()

        ! time iteration
        do t = 2, nt+1

            call computePoroHarm(1) ! old
            call computePoro()
            call computePoroHarm(2) ! new

            call computeK()
            call computeKHarm()

            call computeav()

            call computekc()
            call genDynPara_vp()
            call computevp()

            call computekc()
            call genDynPara_Cf()
            call computeCf()

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

