
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

    subroutine initialize(modelCase)

        type(model), intent(in out) :: modelCase

        integer :: xmomeSize, ymomeSize, contiSize, cfSize, temSize
        ! The arrays that will communicate between different functions must use pointer type instead of allocatalbe
        ! type. Pointer type can make sure that the subscripts of the arrays keep the same in the calling and called
        ! functions. However, if the subscripts of the arrays begin with 1, using allocatable type is also OK.
        integer :: indexl, indexr, indexd, indexu
        integer :: global_ind_b, global_ind_e
        logical :: alive
        integer :: i, j, n, c, ierr

        call MPI_Init(ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nProcs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)

#if defined(MUMPS) || defined(HYPRE)
        call MPI_BUFFER_ATTACH(buffer, buffer_size, ierr)
#endif

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timestart = MPI_Wtime()
        solvertime = 0.D0

        pncols = modelcase%pncols
        pnrows = modelcase%pnrows
        Lx = modelCase%Lx
        Ly = modelCase%Ly
        nx = modelCase%nx
        ny = modelCase%ny
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

        allocate(xs(nx+1))
        allocate(ys(ny+1))
        allocate(ts(nt+1))
        allocate(src(nx,ny))
        allocate(poroInit(nx,ny))
        allocate(KxxInit(nx,ny))
        allocate(KyyInit(nx,ny))
        allocate(avInit(nx,ny))
        allocate(vxBdryX0(ny))
        allocate(vxBdryX1(ny))
        allocate(vxBdryY0(nx+1))
        allocate(vxBdryY1(nx+1))
        allocate(vyBdryX0(ny+1))
        allocate(vyBdryX1(ny+1))
        allocate(vyBdryY0(nx))
        allocate(vyBdryY1(nx))
        allocate(isDiriX0_p(ny))
        allocate(isDiriX1_p(ny))
        allocate(isDiriY0_p(nx))
        allocate(isDiriY1_p(nx))
        allocate(pBdryX0(ny))
        allocate(pBdryX1(ny))
        allocate(pBdryY0(nx))
        allocate(pBdryY1(nx))
        allocate(pInit(nx,ny))
        allocate(isDiriX0_Cf(ny))
        allocate(isDiriX1_Cf(ny))
        allocate(isDiriY0_Cf(nx))
        allocate(isDiriY1_Cf(nx))
        allocate(CfBdryX0(ny))
        allocate(CfBdryX1(ny))
        allocate(CfBdryY0(nx))
        allocate(CfBdryY1(nx))
        allocate(CfInit(nx,ny))
        allocate(isDiriX0_Tem(ny))
        allocate(isDiriX1_Tem(ny))
        allocate(isDiriY0_Tem(nx))
        allocate(isDiriY1_Tem(nx))
        allocate(TemBdryX0(ny))
        allocate(TemBdryX1(ny))
        allocate(TemBdryY0(nx))
        allocate(TemBdryY1(nx))
        allocate(TemInit(nx,ny))

        xs = modelCase%xs
        ys = modelCase%ys
        ts = modelCase%ts
        src = modelCase%src
        poroInit = modelCase%poroInit
        KxxInit = modelCase%KxxInit
        KyyInit = modelCase%KyyInit
        avInit = modelCase%avInit
        vxBdryX0 = modelCase%vxBdryX0
        vxBdryX1 = modelCase%vxBdryX1
        vxBdryY0 = modelCase%vxBdryY0
        vxBdryY1 = modelCase%vxBdryY1
        vyBdryX0 = modelCase%vyBdryX0
        vyBdryX1 = modelCase%vyBdryX1
        vyBdryY0 = modelCase%vyBdryY0
        vyBdryY1 = modelCase%vyBdryY1
        isDiriX0_p = modelCase%isDiriX0_p
        isDiriX1_p = modelCase%isDiriX1_p
        isDiriY0_p = modelCase%isDiriY0_p
        isDiriY1_p = modelCase%isDiriY1_p
        pBdryX0 = modelCase%pBdryX0
        pBdryX1 = modelCase%pBdryX1
        pBdryY0 = modelCase%pBdryY0
        pBdryY1 = modelCase%pBdryY1
        pInit = modelCase%pInit
        isDiriX0_Cf = modelCase%isDiriX0_Cf
        isDiriX1_Cf = modelCase%isDiriX1_Cf
        isDiriY0_Cf = modelCase%isDiriY0_Cf
        isDiriY1_Cf = modelCase%isDiriY1_Cf
        CfBdryX0 = modelCase%CfBdryX0
        CfBdryX1 = modelCase%CfBdryX1
        CfBdryY0 = modelCase%CfBdryY0
        CfBdryY1 = modelCase%CfBdryY1
        CfInit = modelCase%CfInit
        isDiriX0_Tem = modelCase%isDiriX0_Tem
        isDiriX1_Tem = modelCase%isDiriX1_Tem
        isDiriY0_Tem = modelCase%isDiriY0_Tem
        isDiriY1_Tem = modelCase%isDiriY1_Tem
        TemBdryX0 = modelCase%TemBdryX0
        TemBdryX1 = modelCase%TemBdryX1
        TemBdryY0 = modelCase%TemBdryY0
        TemBdryY1 = modelCase%TemBdryY1
        TemInit = modelCase%TemInit
        soludoc = modelCase%soludoc

        localncols = nx/pncols
        localnrows = ny/pnrows

        pcol = mod(myid,pncols)+1
        prow = myid/pncols+1

        xlower = (pcol-1)*localncols+1
        xupper = pcol*localncols
        ylower = (prow-1)*localnrows+1
        yupper = prow*localnrows

        call index_convert_local_global(myid, 1, 1, 1, ilower)
        call index_convert_local_global(myid, 3, localncols, localnrows, iupper)
        local_x_size = iupper - ilower + 1

        if((nProcs>1).and.(myid==0)) then
            allocate(slave_vp_data_size(nProcs-1))
            do i = 1, nProcs-1
                call index_convert_local_global(i, 1, 1, 1, global_ind_b)
                call index_convert_local_global(i, 3, localncols, localnrows, global_ind_e)
                slave_vp_data_size(i) = global_ind_e - global_ind_b + 1
            end do
        end if

        call index_convert_local_global(myid, 4, 1, 1, ilower_Cf)
        call index_convert_local_global(myid, 4, localncols, localnrows, iupper_Cf)
        local_x_size_Cf = iupper_Cf - ilower_Cf + 1

        call index_convert_local_global(myid, 5, 1, 1, ilower_Tem)
        call index_convert_local_global(myid, 5, localncols, localnrows, iupper_Tem)
        local_x_size_Tem = iupper_Tem - ilower_Tem + 1

        t = 2

        ! initialize dm, kc, ks
        allocate(dm(0:localncols+1, 0:localnrows+1))
        allocate(kc(1:localncols, 1:localnrows))
        allocate(ks(1:localncols, 1:localnrows))

        ! initialize hx, hy
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
        
        ! initialize poro
        allocate(poro(indexl:indexr,indexd:indexu))
        allocate(poro_old(indexl:indexr,indexd:indexu))
        do j = indexd, indexu
            do i = indexl, indexr
                poro(i,j) = poroInit(xlower+i-1,ylower+j-1)
                poro_old(i,j) = poro(i,j)
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
        allocate(poroHarmX(indexl:indexr,indexd:indexu))
        allocate(poroHarmX_old(indexl:indexr,indexd:indexu))
        allocate(poroHarmXInit(indexl:indexr,indexd:indexu))

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
        allocate(poroHarmY(indexl:indexr,indexd:indexu))
        allocate(poroHarmY_old(indexl:indexr,indexd:indexu))
        allocate(poroHarmYInit(indexl:indexr,indexd:indexu))

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
        allocate(Kxx(indexl:indexr,indexd:indexu))
        allocate(Kyy(indexl:indexr,indexd:indexu))
        do j = indexd, indexu
            do i = indexl, indexr
                Kxx(i,j) = KxxInit(xlower+i-1,ylower+j-1)
                Kyy(i,j) = KyyInit(xlower+i-1,ylower+j-1)
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

        allocate(KxxHarm(1:indexr,1:localnrows))
        allocate(KyyHarm(1:localncols,1:indexu))

        ! initialize av
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        allocate(av(indexl:indexr,indexd:indexu))
        do j = indexd, indexu
            do i = indexl, indexr
                av(i,j) = avInit(xlower+i-1,ylower+j-1)
            end do
        end do

        ! initialize vx, vy
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
        allocate(vx(indexl:indexr,indexd:indexu))
        vx(:,:) = 0.D0
        if(pcol == 1) then
            do j = 1, localnrows
                if(isDiriX0_p(ylower+j-1) == 0) then
                    vx(1,j) = vxBdryX0(ylower+j-1)
                end if
            end do
        end if
        if(pcol == pncols) then
            do j = 1, localnrows
                if(isDiriX1_p(ylower+j-1) == 0) then
                    vx(localncols+1,j) = vxBdryX1(ylower+j-1)
                end if
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
        allocate(vy(indexl:indexr,indexd:indexu))
        vy(:,:) = 0.D0
        if(prow == 1) then
            do i = 1, localncols
                if(isDiriY0_p(xlower+i-1) == 0) then
                    vy(i,1) = vyBdryY0(xlower+i-1)
                end if
            end do
        end if
        if(prow == pnrows) then
            do i = 1, localncols
                if(isDiriY1_p(xlower+i-1) == 0) then
                    vy(i,localnrows+1) = vyBdryY1(xlower+i-1)
                end if
            end do
        end if

        ! initialize p
        allocate(p(1:localncols,1:localnrows))
        p(1:localncols,1:localnrows) = pInit(xlower:xupper,ylower:yupper)

        ! initialize Cf
        allocate(Cf(1:localncols,1:localnrows))
        Cf(1:localncols,1:localnrows) = CfInit(xlower:xupper,ylower:yupper)

        ! initialize Tem
        allocate(Tem(1:localncols,1:localnrows))
        Tem(1:localncols,1:localnrows) = TemInit(xlower:xupper,ylower:yupper)

        ! compute the size of the equations
        if(pcol /= pncols) then
            xmomeSize = localncols*localnrows
        else
            xmomeSize = (localncols+1)*localnrows
        end if
        if(prow /= pnrows) then
            ymomeSize = localncols*localnrows
        else
            ymomeSize = localncols*(localnrows+1)
        end if
        contiSize = localncols*localnrows
        cfSize = localncols*localnrows
        temSize = localncols*localnrows

        allocate(AxxEntryNum(xmomeSize))
        allocate(AxpEntryNum(xmomeSize))
        allocate(AyyEntryNum(ymomeSize))
        allocate(AypEntryNum(ymomeSize))
        allocate(AcxEntryNum(contiSize))
        allocate(AcyEntryNum(contiSize))
        allocate(AcpEntryNum(contiSize))
        allocate(AcfEntryNum(cfSize))
        allocate(AtemEntryNum(temSize))

        AxxEntryNum(:) = 5
        AxpEntryNum(:) = 2
        AyyEntryNum(:) = 5
        AypEntryNum(:) = 2
        AcxEntryNum(:) = 2
        AcyEntryNum(:) = 2
        AcpEntryNum(:) = 1
        AcfEntryNum(:) = 9
        AtemEntryNum(:) = 5

        if(pcol == pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if

        if(pcol == 1) then
            do n = 0, localnrows-1
                AxxEntryNum(1+indexr*n) = AxxEntryNum(1+indexr*n) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                AxxEntryNum(indexr*n) = AxxEntryNum(indexr*n) - 1
            end do
        end if
        if(prow == 1) then
            do n = 1, indexr
                AxxEntryNum(n) = AxxEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, indexr
                AxxEntryNum(indexr*(localnrows-1)+n) = AxxEntryNum(indexr*(localnrows-1)+n) - 1
            end do
        end if
        do j = 1, localnrows
            if((isDiriX0_p(ylower+j-1)==0).and.(pcol==1)) then
                if(((j==1).and.(prow==1)).or.((j==localnrows).and.(prow==pnrows))) then
                    AxxEntryNum(1+indexr*(j-1)) = AxxEntryNum(1+indexr*(j-1)) - 2
                else
                    AxxEntryNum(1+indexr*(j-1)) = AxxEntryNum(1+indexr*(j-1)) - 3
                end if
            end if
            if((isDiriX1_p(ylower+j-1)==0).and.(pcol==pncols)) then
                if(((j==1).and.(prow==1)).or.((j==localnrows).and.(prow==pnrows))) then
                    AxxEntryNum(indexr*j) = AxxEntryNum(indexr*j) - 2
                else
                    AxxEntryNum(indexr*j) = AxxEntryNum(indexr*j) - 3
                end if
            end if
        end do

        if(pcol == 1) then
            do n = 0, localnrows-1
                AxpEntryNum(1+indexr*n) = AxpEntryNum(1+indexr*n) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                AxpEntryNum(indexr*n) = AxpEntryNum(indexr*n) - 1
            end do
        end if
        do j = 1, localnrows
            if((pcol==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                AxpEntryNum(1+indexr*(j-1)) = AxpEntryNum(1+indexr*(j-1)) - 1
            end if
            if((pcol==pncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                AxpEntryNum(indexr*j) = AxpEntryNum(indexr*j) - 1
            end if
        end do

        if(prow == pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        if(prow == 1) then
            do n = 1, localncols
                AyyEntryNum(n) = AyyEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                AyyEntryNum(localncols*localnrows+n) = AyyEntryNum(localncols*localnrows+n) - 1
            end do
        end if
        if(pcol == 1) then
            do n = 1, indexu
                AyyEntryNum(1+(n-1)*localncols) = AyyEntryNum(1+(n-1)*localncols) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, indexu
                AyyEntryNum(n*localncols) = AyyEntryNum(n*localncols) - 1
            end do
        end if
        do i = 1, localncols
            if((isDiriY0_p(xlower+i-1)==0).and.(prow==1)) then
                if(((i==1).and.(pcol==1)).or.((i==localncols).and.(pcol==pncols))) then
                    AyyEntryNum(i) = AyyEntryNum(i) - 2
                else
                    AyyEntryNum(i) = AyyEntryNum(i) - 3
                end if
            end if
            if((isDiriY1_p(xlower+i-1)==0).and.(prow==pnrows)) then
                if(((i==1).and.(pcol==1)).or.((i==localncols).and.(pcol==pncols))) then
                    AyyEntryNum(localncols*localnrows+i) = AyyEntryNum(localncols*localnrows+i) - 2
                else
                    AyyEntryNum(localncols*localnrows+i) = AyyEntryNum(localncols*localnrows+i) - 3
                end if
            end if
        end do

        if(prow == 1) then
            do n = 1, localncols
                AypEntryNum(n) = AypEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                AypEntryNum(localncols*localnrows+n) = AypEntryNum(localncols*localnrows+n) - 1
            end do
        end if
        do i = 1, localncols
            if((prow==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                AypEntryNum(i) = AypEntryNum(i) - 1
            end if
            if((prow==pnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                AypEntryNum(localncols*localnrows+i) = AypEntryNum(localncols*localnrows+i) - 1
            end if
        end do

        if(pcol == 1) then
            do n = 0, localnrows-1
                AcfEntryNum(1+localncols*n) = AcfEntryNum(1+localncols*n) - 3
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                AcfEntryNum(localncols*n) = AcfEntryNum(localncols*n) - 3
            end do
        end if
        if(prow == 1) then
            do n = 1, localncols
                AcfEntryNum(n) = AcfEntryNum(n) - 3
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                AcfEntryNum(localncols*(localnrows-1)+n) = AcfEntryNum(localncols*(localnrows-1)+n) - 3
            end do
        end if
        if((pcol==1).and.(prow==1)) then
            AcfEntryNum(1) = AcfEntryNum(1) + 1
        end if
        if((pcol==1).and.(prow==pnrows)) then
            AcfEntryNum(1+localncols*(localnrows-1)) = AcfEntryNum(1+localncols*(localnrows-1)) + 1
        end if
        if((pcol==pncols).and.(prow==1)) then
            AcfEntryNum(localncols) = AcfEntryNum(localncols) + 1
        end if
        if((pcol==pncols).and.(prow==pnrows)) then
            AcfEntryNum(localncols*localnrows) = AcfEntryNum(localncols*localnrows) + 1
        end if

        if(pcol == 1) then
            do n = 0, localnrows-1
                AtemEntryNum(1+localncols*n) = AtemEntryNum(1+localncols*n) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                AtemEntryNum(localncols*n) = AtemEntryNum(localncols*n) - 1
            end do
        end if
        if(prow == 1) then
            do n = 1, localncols
                AtemEntryNum(n) = AtemEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                AtemEntryNum(localncols*(localnrows-1)+n) = AtemEntryNum(localncols*(localnrows-1)+n) - 1
            end do
        end if

        allocate(AxxEntryBase(xmomeSize))
        allocate(AxpEntryBase(xmomeSize))
        allocate(AyyEntryBase(ymomeSize))
        allocate(AypEntryBase(ymomeSize))
        allocate(AcxEntryBase(contiSize))
        allocate(AcyEntryBase(contiSize))
        allocate(AcpEntryBase(contiSize))
        allocate(AcfEntryBase(cfSize))
        allocate(AtemEntryBase(temSize))

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

        AcxEntryBase(1) = 1
        AcyEntryBase(1) = 1
        AcpEntryBase(1) = 1
        do n = 2, contiSize
            AcxEntryBase(n) = AcxEntryBase(n-1) + AcxEntryNum(n-1)
            AcyEntryBase(n) = AcyEntryBase(n-1) + AcyEntryNum(n-1)
            AcpEntryBase(n) = AcpEntryBase(n-1) + AcpEntryNum(n-1)
        end do

        AcfEntryBase(1) = 1
        do n = 2, cfSize
            AcfEntryBase(n) = AcfEntryBase(n-1) + AcfEntryNum(n-1)
        end do

        AtemEntryBase(1) = 1
        do n = 2, temSize
            AtemEntryBase(n) = AtemEntryBase(n-1) + AtemEntryNum(n-1)
        end do

        ! compute matrix size
        AxxSize = 0
        do n = 1, xmomeSize
            AxxSize = AxxSize + AxxEntryNum(n)
        end do
        AxpSize = 0
        do n = 1, xmomeSize
            AxpSize = AxpSize + AxpEntryNum(n)
        end do
        AyySize = 0
        do n = 1, ymomeSize
            AyySize = AyySize + AyyEntryNum(n)
        end do
        AypSize = 0
        do n = 1, ymomeSize
            AypSize = AypSize + AypEntryNum(n)
        end do
        AcxSize = 2*localncols*localnrows
        AcySize = 2*localncols*localnrows
        AcpSize = localncols*localnrows
        AcfSize = 0
        do n = 1, cfSize
            AcfSize = AcfSize + AcfEntryNum(n)
        end do
        AtemSize = 0
        do n = 1, temSize
            AtemSize = AtemSize + AtemEntryNum(n)
        end do

        ! initialize matrix
        allocate(local_rhs(local_x_size))
        allocate(local_rhs_static(local_x_size))
        allocate(local_rhs_Cf(local_x_size_Cf))
        allocate(local_rhs_Tem(local_x_size_Tem))

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
        allocate(AcxCols(AcxSize))
        allocate(AcxRows(AcxSize))
        allocate(AcxValues(AcxSize))
        allocate(AcyCols(AcySize))
        allocate(AcyRows(AcySize))
        allocate(AcyValues(AcySize))
        allocate(AcpCols(AcpSize))
        allocate(AcpRows(AcpSize))
        allocate(AcpValues(AcpSize))
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
        AcxCols(:) = 0
        AcxRows(:) = 0
        AcxValues(:) = 0.D0
        AcyCols(:) = 0
        AcyRows(:) = 0
        AcyValues(:) = 0.D0
        AcpCols(:) = 0
        AcpRows(:) = 0
        AcpValues(:) = 0.D0
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

#ifdef LAPACK

        allocate(A_lapack(local_x_size,local_x_size))
        allocate(b_lapack(local_x_size))
        allocate(IPIV(local_x_size))

        allocate(A_lapack_Cf(local_x_size_Cf,local_x_size_Cf))
        allocate(b_lapack_Cf(local_x_size_Cf))
        allocate(IPIV_Cf(local_x_size_Cf))

        allocate(A_lapack_Tem(local_x_size_Tem,local_x_size_Tem))
        allocate(b_lapack_Tem(local_x_size_Tem))
        allocate(IPIV_Tem(local_x_size_Tem))

#elif defined(UMFPACK)

        allocate(Ap(1:local_x_size+1))
        allocate(Ai(1:7*local_x_size))
        allocate(Ax(1:7*local_x_size))

        allocate(Ap_Cf(1:local_x_size_Cf+1))
        allocate(Ai_Cf(1:9*local_x_size_Cf))
        allocate(Ax_Cf(1:9*local_x_size_Cf))

        allocate(Ap_Tem(1:local_x_size_Tem+1))
        allocate(Ai_Tem(1:5*local_x_size_Tem))
        allocate(Ax_Tem(1:5*local_x_size_Tem))

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
            mumps_par%N = (nx+1)*ny+nx*(ny+1)+nx*ny
        end if
        allocate(mumps_par%RHS(mumps_par%N))

        allocate(mumps_IRN_loc(local_x_size*7))
        allocate(mumps_JCN_loc(local_x_size*7))
        allocate(mumps_A_loc(local_x_size*7))

        mumps_par_Cf%COMM = MPI_COMM_WORLD
        mumps_par_Cf%SYM = 0
        mumps_par_Cf%PAR = 1
        mumps_par_Cf%JOB = -1
        call DMUMPS(mumps_par_Cf)

        mumps_par_Cf%ICNTL(4) = 0
        mumps_par_Cf%ICNTL(7) = 15
        mumps_par_Cf%ICNTL(18) = 3
        if(mumps_par_Cf%MYID == 0) then
            mumps_par_Cf%N = nx*ny
        end if
        allocate(mumps_par_Cf%RHS(mumps_par_Cf%N))

        allocate(mumps_IRN_loc_Cf(local_x_size_Cf*9))
        allocate(mumps_JCN_loc_Cf(local_x_size_Cf*9))
        allocate(mumps_A_loc_Cf(local_x_size_Cf*9))

        mumps_par_Tem%COMM = MPI_COMM_WORLD
        mumps_par_Tem%SYM = 0
        mumps_par_Tem%PAR = 1
        mumps_par_Tem%JOB = -1
        call DMUMPS(mumps_par_Tem)

        mumps_par_Tem%ICNTL(4) = 0
        mumps_par_Tem%ICNTL(7) = 15
        mumps_par_Tem%ICNTL(18) = 3
        if(mumps_par_Tem%MYID == 0) then
            mumps_par_Tem%N = nx*ny
        end if
        allocate(mumps_par_Tem%RHS(mumps_par_Tem%N))

        allocate(mumps_IRN_loc_Tem(local_x_size_Tem*5))
        allocate(mumps_JCN_loc_Tem(local_x_size_Tem*5))
        allocate(mumps_A_loc_Tem(local_x_size_Tem*5))

#elif defined(HYPRE)

        if(prow /= 1) then
            call index_convert_local_global(myid-pncols, 1, 1, localnrows, jlower)
        elseif(pcol /= 1) then
            call index_convert_local_global(myid-1, 1, localncols, 1, jlower)
        else
            jlower = ilower
        end if
        if(prow /= pnrows) then
            call index_convert_local_global(myid+pncols, 3, localncols, 1, jupper)
        elseif(pcol /= pncols) then
            call index_convert_local_global(myid+1, 3, 1, localnrows, jupper)
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

        if((pcol /= 1).and.(prow /= 1)) then
            call index_convert_local_global(myid-1-pncols, 4, localncols, localnrows, jlower_Cf)
        elseif((pcol /= 1).and.(prow == 1)) then
            call index_convert_local_global(myid-1, 4, localncols, 1, jlower_Cf)
        elseif((pcol == 1).and.(prow /= 1)) then
            call index_convert_local_global(myid-pncols, 4, 1, localnrows, jlower_Cf)
        else
            jlower_Cf = ilower_Cf
        end if
        if((pcol /= pncols).and.(prow /= pnrows)) then
            call index_convert_local_global(myid+1+pncols, 4, 1, 1, jupper_Cf)
        elseif((pcol /= pncols).and.(prow == pnrows)) then
            call index_convert_local_global(myid+1, 4, 1, localnrows, jupper_Cf)
        elseif((pcol == pncols).and.(prow /= pnrows)) then
            call index_convert_local_global(myid+pncols, 4, localncols, 1, jupper_Cf)
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

        if(prow /= 1) then
            call index_convert_local_global(myid-pncols, 5, 1, localnrows, jlower_Tem)
        elseif(pcol /= 1) then
            call index_convert_local_global(myid-1, 5, localncols, 1, jlower_Tem)
        else
            jlower_Tem = ilower_Tem
        end if
        if(prow /= pnrows) then
            call index_convert_local_global(myid+pncols, 5, localncols, 1, jupper_Tem)
        elseif(pcol /= pncols) then
            call index_convert_local_global(myid+1, 5, 1, localnrows, jupper_Tem)
        else
            jupper_Tem = iupper_Tem
        end if

        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower_Tem, iupper_Tem, jlower_Tem, jupper_Tem, A_Tem, ierr)
        call HYPRE_IJMatrixSetObjectType(A_Tem, HYPRE_PARCSR, ierr)
        call HYPRE_IJMatrixInitialize(A_Tem, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Tem, iupper_Tem, b_Tem, ierr)
        call HYPRE_IJVectorSetObjectType(b_Tem, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(b_Tem, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Tem, iupper_Tem, x_Tem, ierr)
        call HYPRE_IJVectorSetObjectType(x_Tem, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(x_Tem, ierr)

        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver, ierr)
        call HYPRE_ParCSRGMRESSetPrintLevel(solver, 2, ierr)
        call HYPRE_ParCSRGMRESSetLogging(solver, 1, ierr)

        call HYPRE_ParaSailsCreate(MPI_COMM_WORLD, precond,ierr)
        call HYPRE_ParaSailsSetParams(precond, -9D-1, 2, ierr)
        call HYPRE_ParaSailsSetFilter(precond, -9D-1, ierr)
        ! Because the matrix A is nonsymmetric and indefinite, you must choose the parameter as 0.
        call HYPRE_ParaSailsSetSym(precond, 0)
        call HYPRE_ParaSailsSetLogging(precond, 1, ierr)

        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver, 1, precond, ierr)

        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver_Cf, ierr)
        call HYPRE_ParCSRGMRESSetPrintLevel(solver_Cf, 2, ierr)
        call HYPRE_ParCSRGMRESSetLogging(solver_Cf, 1, ierr)

        call HYPRE_ParaSailsCreate(MPI_COMM_WORLD, precond_Cf,ierr)
        call HYPRE_ParaSailsSetParams(precond_Cf, -9D-1, 2, ierr)
        call HYPRE_ParaSailsSetFilter(precond_Cf, -9D-1, ierr)
        ! Because the matrix A is nonsymmetric and indefinite, you must choose the parameter as 0.
        call HYPRE_ParaSailsSetSym(precond_Cf, 0)
        call HYPRE_ParaSailsSetLogging(precond_Cf, 1, ierr)

        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver_Cf, 1, precond_Cf, ierr)

        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver_Tem, ierr)
        call HYPRE_ParCSRGMRESSetPrintLevel(solver_Tem, 2, ierr)
        call HYPRE_ParCSRGMRESSetLogging(solver_Tem, 1, ierr)

        call HYPRE_ParaSailsCreate(MPI_COMM_WORLD, precond_Tem,ierr)
        call HYPRE_ParaSailsSetParams(precond_Tem, -9D-1, 2, ierr)
        call HYPRE_ParaSailsSetFilter(precond_Tem, -9D-1, ierr)
        ! Because the matrix A is nonsymmetric and indefinite, you must choose the parameter as 0.
        call HYPRE_ParaSailsSetSym(precond_Tem, 0)
        call HYPRE_ParaSailsSetLogging(precond_Tem, 1, ierr)

        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver_Tem, 1, precond_Tem, ierr)

        allocate(rows(local_x_size))
        do n = 1, local_x_size
            rows(n) = ilower + n - 1
        end do

        allocate(rows_Cf(local_x_size_Cf))
        do n = 1, local_x_size_Cf
            rows_Cf(n) = ilower_Cf + n - 1
        end do

        allocate(rows_Tem(local_x_size_Tem))
        do n = 1, local_x_size_Tem
            rows_Tem(n) = ilower_Tem + n - 1
        end do

        ! the initial guess of x in the solver iteration
        allocate(initial_x_guess(local_x_size))
        initial_x_guess(:) = 0.D0
        c = local_x_size - localncols*localnrows
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                initial_x_guess(c) = pInit(xlower+i-1,ylower+j-1)
            end do
        end do

        allocate(initial_x_guess_Cf(local_x_size_Cf))
        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                initial_x_guess_Cf(c) = CfInit(xlower+i-1,ylower+j-1)
            end do
        end do

        allocate(initial_x_guess_Tem(local_x_size_Tem))
        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                initial_x_guess_Tem(c) = TemInit(xlower+i-1,ylower+j-1)
            end do
        end do

#endif

        if(mod(nt, NUMFRAME) /= 0) then
            if(myid == 0) then
                print *, 'The number of time steps must be divided by the number of frames.'
            end if
            stop
        end if

        deallocate(modelCase%xs)
        deallocate(modelCase%ys)
        deallocate(modelCase%ts)
        deallocate(modelCase%src)
        deallocate(modelCase%poroInit)
        deallocate(modelCase%KxxInit)
        deallocate(modelCase%KyyInit)
        deallocate(modelCase%avInit)
        deallocate(modelCase%vxBdryX0)
        deallocate(modelCase%vxBdryX1)
        deallocate(modelCase%vxBdryY0)
        deallocate(modelCase%vxBdryY1)
        deallocate(modelCase%vyBdryX0)
        deallocate(modelCase%vyBdryX1)
        deallocate(modelCase%vyBdryY0)
        deallocate(modelCase%vyBdryY1)
        deallocate(modelCase%isDiriX0_p)
        deallocate(modelCase%isDiriX1_p)
        deallocate(modelCase%isDiriY0_p)
        deallocate(modelCase%isDiriY1_p)
        deallocate(modelCase%pBdryX0)
        deallocate(modelCase%pBdryX1)
        deallocate(modelCase%pBdryY0)
        deallocate(modelCase%pBdryY1)
        deallocate(modelCase%pInit)
        deallocate(modelCase%isDiriX0_Cf)
        deallocate(modelCase%isDiriX1_Cf)
        deallocate(modelCase%isDiriY0_Cf)
        deallocate(modelCase%isDiriY1_Cf)
        deallocate(modelCase%CfBdryX0)
        deallocate(modelCase%CfBdryX1)
        deallocate(modelCase%CfBdryY0)
        deallocate(modelCase%CfBdryY1)
        deallocate(modelCase%CfInit)
        deallocate(modelCase%isDiriX0_Tem)
        deallocate(modelCase%isDiriX1_Tem)
        deallocate(modelCase%isDiriY0_Tem)
        deallocate(modelCase%isDiriY1_Tem)
        deallocate(modelCase%TemBdryX0)
        deallocate(modelCase%TemBdryX1)
        deallocate(modelCase%TemBdryY0)
        deallocate(modelCase%TemBdryY1)
        deallocate(modelCase%TemInit)

    end subroutine initialize

    ! Generate the values on the right-hand side and the coefficients of the matrix A,
    ! and the subroutine will generate the values that will not change with the time iteration.
    subroutine genStaticPara_vp()

        integer :: findexl, findexr, findexd, findexu ! field index
        integer :: eindexl, eindexr, eindexd, eindexu ! equation index
        integer :: global_ind
        integer, dimension(:,:), pointer :: velx, vely, pres
        logical :: isField
        real(kind=8), dimension(:,:), pointer :: rhs_velx_b, rhs_dpdx, rhs_vely_b, rhs_dpdy, rhs_dudx, rhs_dvdy
        real(kind=8), dimension(:,:), pointer :: resiAxx_b, resiAcx, resiAyy_b, resiAcy, resiAxp, resiAyp, resitemp
        integer :: i, j, n, c

        if(pcol == 1) then
            findexl = 1
        else
            findexl = 0
        end if
        findexr = localncols + 1
        if(prow == 1) then
            findexd = 1
        else
            findexd = 0
        end if
        if(prow == pnrows) then
            findexu = localnrows
        else
            findexu = localnrows + 1
        end if
        allocate(velx(findexl:findexr,findexd:findexu))
        velx(:,:) = 0

        eindexl = 1
        if(pcol /= pncols) then
            eindexr = localncols
        else
            eindexr = localncols + 1
        end if
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_velx_b(eindexl:eindexr,eindexd:eindexu))
        call Resi_velx_b(velx, rhs_velx_b)
        allocate(resiAxx_b(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, velx, isField)
                if(isField) then
                    call Resi_velx_b(velx, resiAxx_b)
                    resitemp = resiAxx_b - rhs_velx_b
                    call constructAxx(velx, resitemp, 1)
                end if
            end do
        end do
        deallocate(velx)
        deallocate(resiAxx_b)
        deallocate(resitemp)

        do j = 1, localnrows
            if((pcol==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                call index_convert_local_global(myid, 1, 1, j, global_ind)
                do n = 1, AxxSize
                    if(AxxRows(n)==global_ind) then
                        AxxValues(n) = AxxStaticValues(n)
                    end if
                end do
            end if
            if((pcol==pncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                call index_convert_local_global(myid, 1, localncols+1, j, global_ind)
                do n = 1, AxxSize
                    if(AxxRows(n)==global_ind) then
                        AxxValues(n) = AxxStaticValues(n)
                    end if
                end do
            end if
        end do

        if(pcol == 1) then
            findexl = 1
        else
            findexl = 0
        end if
        findexr = localncols
        findexd = 1
        findexu = localnrows
        allocate(pres(findexl:findexr,findexd:findexu))
        pres(:,:) = 0
        allocate(rhs_dpdx(eindexl:eindexr,eindexd:eindexu))
        call Resi_dpdx(pres, rhs_dpdx)
        allocate(resiAxp(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, pres, isField)
                if(isField) then
                    call Resi_dpdx(pres, resiAxp)
                    resitemp = resiAxp - rhs_dpdx
                    call constructAxp(pres, resitemp)
                end if
            end do
        end do
        deallocate(pres)
        deallocate(resiAxp)
        deallocate(resitemp)

        c = 0
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                local_rhs_static(c) = -(rhs_velx_b(i,j) + rhs_dpdx(i,j) - rhof*gravX)
                if(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and. &!
                    (i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0))) then
                    local_rhs(c) = local_rhs_static(c) - rhof*gravX
                end if
            end do
        end do
        deallocate(rhs_velx_b)
        deallocate(rhs_dpdx)

        if(pcol == 1) then
            findexl = 1
        else
            findexl = 0
        end if
        if(pcol == pncols) then
            findexr = localncols
        else
            findexr = localncols + 1
        end if
        if(prow == 1) then
            findexd = 1
        else
            findexd = 0
        end if
        findexu = localnrows + 1
        allocate(vely(findexl:findexr,findexd:findexu))
        vely(:,:) = 0
        eindexl = 1
        eindexr = localncols
        eindexd = 1
        if(prow /= pnrows) then
            eindexu = localnrows
        else
            eindexu = localnrows + 1
        end if
        allocate(rhs_vely_b(eindexl:eindexr,eindexd:eindexu))
        call Resi_vely_b(vely, rhs_vely_b)
        allocate(resiAyy_b(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, vely, isField)
                if(isField) then
                    call Resi_vely_b(vely, resiAyy_b)
                    resitemp = resiAyy_b - rhs_vely_b
                    call constructAyy(vely, resitemp, 1)
                end if
            end do
        end do
        deallocate(vely)
        deallocate(resiAyy_b)
        deallocate(resitemp)

        do i = 1, localncols
            if((prow==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                call index_convert_local_global(myid, 2, i, 1, global_ind)
                do n = 1, AyySize
                    if(AyyRows(n)==global_ind) then
                        AyyValues(n) = AyyStaticValues(n)
                    end if
                end do
            end if
            if((prow==pnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                call index_convert_local_global(myid, 2, i, localnrows+1, global_ind)
                do n = 1, AyySize
                    if(AyyRows(n)==global_ind) then
                        AyyValues(n) = AyyStaticValues(n)
                    end if
                end do
            end if
        end do

        findexl = 1
        findexr = localncols
        if(prow == 1) then
            findexd = 1
        else
            findexd = 0
        end if
        findexu = localnrows
        allocate(pres(findexl:findexr,findexd:findexu))
        pres(:,:) = 0
        allocate(rhs_dpdy(eindexl:eindexr,eindexd:eindexu))
        call Resi_dpdy(pres, rhs_dpdy)
        allocate(resiAyp(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, pres, isField)
                if(isField) then
                    call Resi_dpdy(pres, resiAyp)
                    resitemp = resiAyp - rhs_dpdy
                    call constructAyp(pres, resitemp)
                end if
            end do
        end do
        deallocate(pres)
        deallocate(resiAyp)
        deallocate(resitemp)

        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                local_rhs_static(c) = -(rhs_vely_b(i,j) + rhs_dpdy(i,j) - rhof*gravY)
                if(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and.(j==localnrows+1) &!
                    .and.(isDiriY1_p(xlower+i-1)==0))) then
                    local_rhs(c) = local_rhs_static(c) - rhof*gravY
                end if
            end do
        end do
        deallocate(rhs_vely_b)
        deallocate(rhs_dpdy)

        findexl = 1
        findexr = localncols+1
        findexd = 1
        findexu = localnrows
        allocate(velx(findexl:findexr,findexd:findexu))
        velx(:,:) = 0
        allocate(rhs_dudx(1:localncols,1:localnrows))
        call Resi_dudx(velx, rhs_dudx)
        allocate(resiAcx(1:localncols,1:localnrows))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, velx, isField)
                if(isField) then
                    call Resi_dudx(velx, resiAcx)
                    call constructAcx(velx, resiAcx)
                end if
            end do
        end do
        deallocate(velx)
        deallocate(resiAcx)

        findexl = 1
        findexr = localncols
        findexd = 1
        findexu = localnrows+1
        allocate(vely(findexl:findexr,findexd:findexu))
        vely(:,:) = 0
        allocate(rhs_dvdy(1:localncols,1:localnrows))
        call Resi_dvdy(vely, rhs_dvdy)
        allocate(resiAcy(1:localncols,1:localnrows))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, vely, isField)
                if(isField) then
                    call Resi_dvdy(vely, resiAcy)
                    call constructAcy(vely, resiAcy)
                end if
            end do
        end do
        deallocate(vely)
        deallocate(resiAcy)

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_rhs_static(c) = -(rhs_dudx(i,j) + rhs_dvdy(i,j))
            end do
        end do
        deallocate(rhs_dudx)
        deallocate(rhs_dvdy)

    end subroutine genStaticPara_vp

    subroutine computedm()

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, c, ierr

        do j = 1, localnrows
            do i = 1, localncols

                dm(i,j) = dmRef * exp(Eg/Rg*(1.D0/TemRef_dm-1.D0/Tem(i,j)))

            end do
        end do

        ! send
        ! the 4 edges
        if(pcol /= 1) then
            sentSize = localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = dm(1,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            sentSize = localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = dm(localncols,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            sentSize = localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = dm(i,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = dm(i,localnrows)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        ! receive
        ! the 4 edges
        if(pcol /= pncols) then
            recvSize = localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                dm(localncols+1,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            recvSize = localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                dm(0,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            recvSize = localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                dm(i,localnrows+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                dm(i,0) = recv(c)
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

    end subroutine computedm

    subroutine computeks()

        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                ks(i,j) = ksRef * exp(Eg/Rg*(1.D0/TemRef_ks-1.D0/Tem(i,j)))

            end do
        end do

    end subroutine computeks

    subroutine computekc()

        real(kind=8) :: Sc, vmodulus, radius, Rep
        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                vmodulus = dsqrt((vx(i+1,j)-vx(i,j))**2.D0+(vy(i,j+1)-vy(i,j))**2.D0)

                radius = radiusInit*(1.D0-poroInit(xlower+i-1,ylower+j-1))/ &!
                    (poroInit(xlower+i-1,ylower+j-1)*(1.D0-poro(i,j)))*poro(i,j)

                Rep = 2.D0*vmodulus*radius/(visc/rhof)

                Sc = visc/rhof*dm(i,j)

                kc(i,j) = (ShInfinity+7.D-1*Rep**(1.D0/2.D0)*Sc**(1.D0/3.D0))*dm(i,j)/2.D0/radius

            end do
        end do

    end subroutine computekc

    subroutine computePoro()

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu, requestlu, requestrd, requestru
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        real(kind=8) :: coe
        integer :: i, j, c, ierr

        do j = 1, localnrows
            do i = 1, localncols
                coe = avInit(xlower+i-1,ylower+j-1)*al*Cf(i,j)*kc(i,j)*ks(i,j)*(ts(t)-ts(t-1))/ &!
                    (rhos*(kc(i,j)+ks(i,j))*(1.D0-poroInit(xlower+i-1,ylower+j-1)))
                poro_old(i,j) = poro(i,j)
                poro(i,j) = (coe+poro_old(i,j))/(1+coe)
            end do
        end do

        if(pcol /= 1) then
            sentSize = 2*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = poro(1,j)
                c = c + 1
                sent(c) = poro_old(1,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            sentSize = 4*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                do i = localncols-1, localncols
                    c = c + 1
                    sent(c) = poro(i,j)
                    c = c + 1
                    sent(c) = poro_old(i,j)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            sentSize = 2*localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = poro(i,1)
                c = c + 1
                sent(c) = poro_old(i,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = localncols*4
            allocate(sent(sentSize))
            c = 0
            do j = localnrows-1, localnrows
                do i = 1, localncols
                    c = c + 1
                    sent(c) = poro(i,j)
                    c = c + 1
                    sent(c) = poro_old(i,j)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(sent(2))
            sent(1) = poro(1,localnrows)
            sent(2) = poro_old(1,localnrows)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid-1+pncols, myid, MPI_COMM_WORLD, requestlu, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(sent(2))
            sent(1) = poro(localncols,1)
            sent(2) = poro_old(localncols,1)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid+1-pncols, myid, MPI_COMM_WORLD, requestrd, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=pnrows)) then
            allocate(sent(2))
            sent(1) = poro(localncols,localnrows)
            sent(2) = poro_old(localncols,localnrows)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid+1+pncols, myid, MPI_COMM_WORLD, requestru, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            recvSize = 2*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                poro(localncols+1,j) = recv(c)
                c = c + 1
                poro_old(localncols+1,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            recvSize = 4*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                do i = -1, 0
                    c = c + 1
                    poro(i,j) = recv(c)
                    c = c + 1
                    poro_old(i,j) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            recvSize = 2*localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                poro(i,localnrows+1) = recv(c)
                c = c + 1
                poro_old(i,localnrows+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = localncols*4
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = -1, 0
                do i = 1, localncols
                    c = c + 1
                    poro(i,j) = recv(c)
                    c = c + 1
                    poro_old(i,j) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid+1-pncols, myid+1-pncols, MPI_COMM_WORLD, status, ierr)
            poro(localncols+1,0) = recv(1)
            poro_old(localncols+1,0) = recv(2)
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid-1+pncols, myid-1+pncols, MPI_COMM_WORLD, status, ierr)
            poro(0,localnrows+1) = recv(1)
            poro_old(0,localnrows+1) = recv(2)
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=1)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid-1-pncols, myid-1-pncols, MPI_COMM_WORLD, status, ierr)
            poro(0,0) = recv(1)
            poro_old(0,0) = recv(2)
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
        if((pcol/=1).and.(prow/=pnrows)) then
            call MPI_WAIT(requestlu, status, ierr)
        end if
        if((pcol/=pncols).and.(prow/=1)) then
            call MPI_WAIT(requestrd, status, ierr)
        end if
        if((pcol/=pncols).and.(prow/=pnrows)) then
            call MPI_WAIT(requestru, status, ierr)
        end if

    end subroutine computePoro

    subroutine computePoroHarm(var)

        integer, intent(in) :: var

        real(kind=8), dimension(:,:), pointer :: HarmX, HarmY
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        if(var == 1) then
            HarmX => poroHarmX_old
            HarmY => poroHarmY_old
        elseif(var == 2) then
            HarmX => poroHarmX
            HarmY => poroHarmY
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

        if(pcol == 1) then
            HarmX(1,indexd:indexu) = poro(1,indexd:indexu)
            indexl = indexl + 1
        end if

        if(pcol == pncols) then
            HarmX(localncols+1,indexd:indexu) = poro(localncols,indexd:indexu)
            indexr = indexr - 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr
                HarmX(i,j) = (hx(i-1)+hx(i)) / (hx(i-1)/poro(i-1,j)+hx(i)/poro(i,j))
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

        if(prow == 1) then
            HarmY(indexl:indexr,1) = poro(indexl:indexr,1)
            indexd = indexd + 1
        end if

        if(prow == pnrows) then
            HarmY(indexl:indexr,localnrows+1) = poro(indexl:indexr,localnrows)
            indexu = indexu - 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr
                HarmY(i,j) = (hy(j-1)+hy(j)) / (hy(j-1)/poro(i,j-1)+hy(j)/poro(i,j))
            end do
        end do

        if((t == 2).and.(var == 1)) then
            poroHarmXInit = HarmX
            poroHarmYInit = HarmY
        end if

    end subroutine computePoroHarm

    subroutine computeK()

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestr, requestu
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, c, ierr

        do j = 1, localnrows
            do i = 1, localncols
                Kxx(i,j) = poro(i,j)/poroInit(xlower+i-1,ylower+j-1)*(poro(i,j)*(1.D0-poroInit(xlower+i-1,ylower+j-1)) &!
                    /poroInit(xlower+i-1,ylower+j-1)/(1.D0-poro(i,j)))**2 * KxxInit(xlower+i-1,ylower+j-1)
                Kyy(i,j) = poro(i,j)/poroInit(xlower+i-1,ylower+j-1)*(poro(i,j)*(1.D0-poroInit(xlower+i-1,ylower+j-1)) &!
                    /poroInit(xlower+i-1,ylower+j-1)/(1.D0-poro(i,j)))**2 * KyyInit(xlower+i-1,ylower+j-1)
            end do
        end do

        if(pcol /= pncols) then
            sentSize = 2*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = Kxx(localncols,j)
            end do
            do j = 1, localnrows
                c = c + 1
                sent(c) = Kyy(localncols,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = 2*localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = Kxx(i,localnrows)
            end do
            do i = 1, localncols
                c = c + 1
                sent(c) = Kyy(i,localnrows)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if(pcol /= 1) then
            recvSize = 2*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                Kxx(0,j) = recv(c)
            end do
            do j = 1, localnrows
                c = c + 1
                Kyy(0,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = 2*localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                Kxx(i,0) = recv(c)
            end do
            do i = 1, localncols
                c = c + 1
                Kyy(i,0) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if

    end subroutine computeK

    subroutine computeKHarm()

        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        indexl = 1
        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        indexd = 1
        indexu = localnrows

        if(pcol == 1) then
            KxxHarm(1,indexd:indexu) = Kxx(1,indexd:indexu)
            indexl = indexl + 1
        end if

        if(pcol == pncols) then
            KxxHarm(localncols+1,indexd:indexu) = Kxx(localncols,indexd:indexu)
            indexr = indexr - 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr
                KxxHarm(i,j) = (hx(i-1)+hx(i)) / (hx(i-1)/Kxx(i-1,j)+hx(i)/Kxx(i,j))
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        if(prow /= pnrows) then
            indexu = indexu
        else
            indexu = indexu + 1
        end if

        if(prow == 1) then
            KyyHarm(indexl:indexr,1) = Kyy(indexl:indexr,1)
            indexd = indexd + 1
        end if

        if(prow == pnrows) then
            KyyHarm(indexl:indexr,localnrows+1) = Kyy(indexl:indexr,localnrows)
            indexu = indexu - 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr
                KyyHarm(i,j) = (hy(j-1)+hy(j)) / (hy(j-1)/Kyy(i,j-1)+hy(j)/Kyy(i,j))
            end do
        end do

    end subroutine computeKHarm

    subroutine computeav()

        integer :: i, j

        ! suppose kxx=kyy. if you find the anisotropic permeability equation to compute av, you can change the
        ! computation here.
        do j = 1, localnrows
            do i = 1, localncols
                av(i,j) = avInit(xlower+i-1,ylower+j-1)*poro(i,j)/poroInit(xlower+i-1,ylower+j-1)* &!
                    dsqrt(1.D0*KxxInit(xlower+i-1,ylower+j-1)*poro(i,j)/Kxx(i,j)/poroInit(xlower+i-1,ylower+j-1))
            end do
        end do

    end subroutine computeav

    ! Generate the values on the right-hand side and the coefficients of the matrix A,
    ! and the subroutine will generate the values that will change with the time iteration.
    subroutine genDynPara_vp()
       
        integer :: findexl, findexr, findexd, findexu
        integer :: eindexl, eindexr, eindexd, eindexu
        integer :: global_ind
        integer, dimension(:,:), pointer :: velx, vely, pres
        logical :: isField
        real(kind=8), dimension(:,:), pointer :: rhs_velx, rhs_vely, rhs_dpdt
        real(kind=8), dimension(:,:), pointer :: resiAxx, resiAyy, resiAcp, resitemp
        integer :: AxxBeInd, AyyBeInd
        integer :: i, j, n, c
        
        if(pcol == 1) then
            findexl = 1
        else
            findexl = 0
        end if
        findexr = localncols + 1
        if(prow == 1) then
            findexd = 1
        else
            findexd = 0
        end if
        if(prow == pnrows) then
            findexu = localnrows
        else
            findexu = localnrows + 1
        end if
        allocate(velx(findexl:findexr,findexd:findexu))
        velx(:,:) = 0

        eindexl = 1
        if(pcol /= pncols) then
            eindexr = localncols
        else
            eindexr = localncols + 1
        end if
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_velx(eindexl:eindexr,eindexd:eindexu))
        call Resi_velx(velx, rhs_velx)
        allocate(resiAxx(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, velx, isField)
                if(isField) then
                    call Resi_velx(velx, resiAxx)
                    ! generate the dynamic coefficients of Axx. The scheme is just like the DNA copy according
                    ! to a template, here the AxxCols and the AxxRows are just the template.
                    resitemp = resiAxx - rhs_velx
                    call constructAxx(velx, resitemp, 2)
                end if
            end do
        end do
        deallocate(velx)
        deallocate(resiAxx)
        deallocate(resitemp)

        AxxBeInd = 1
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0))).and.(.not.((pcol==pncols) &!
                    .and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                    call index_convert_local_global(myid, 1, i, j, global_ind)
                    do n = AxxBeInd, AxxSize
                        if(AxxRows(n)==global_ind) then
                            AxxValues(n) = AxxStaticValues(n) + AxxDynValues(n)
                        else
                            AxxBeInd = n
                            exit
                        end if
                    end do
                elseif(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and. &!
                    (i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0))) then
                    AxxBeInd = AxxBeInd + 1
                end if
            end do
        end do

        c = 0
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                if(.not.(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and. &!
                    (i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                    local_rhs(c) = local_rhs_static(c) - rhs_velx(i,j)
                end if
            end do
        end do
        deallocate(rhs_velx)

        if(pcol == 1) then
            findexl = 1
        else
            findexl = 0
        end if
        if(pcol == pncols) then
            findexr = localncols
        else
            findexr = localncols + 1
        end if
        if(prow == 1) then
            findexd = 1
        else
            findexd = 0
        end if
        findexu = localnrows + 1
        allocate(vely(findexl:findexr,findexd:findexu))
        vely(:,:) = 0

        eindexl = 1
        eindexr = localncols
        eindexd = 1
        if(prow /= pnrows) then
            eindexu = localnrows
        else
            eindexu = localnrows + 1
        end if
        allocate(rhs_vely(eindexl:eindexr,eindexd:eindexu))
        call Resi_vely(vely, rhs_vely)
        allocate(resiAyy(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, vely, isField)
                if(isField) then
                    call Resi_vely(vely, resiAyy)
                    resitemp = resiAyy - rhs_vely
                    call constructAyy(vely, resitemp, 2)
                end if
            end do
        end do
        deallocate(vely)
        deallocate(resiAyy)
        deallocate(resitemp)

        AyyBeInd = 1
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0))).and.(.not.((prow==pnrows).and. &!
                    (j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                    call index_convert_local_global(myid, 2, i, j, global_ind)
                    do n = AyyBeInd, AyySize
                        if(AyyRows(n)==global_ind) then
                            AyyValues(n) = AyyStaticValues(n) + AyyDynValues(n)
                        else
                            AyyBeInd = n
                            exit
                        end if
                    end do
                elseif(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and. &!
                    (j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0))) then
                    AyyBeInd = AyyBeInd + 1
                end if
            end do
        end do

        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                if(.not.(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and. &!
                    (j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                    local_rhs(c) = local_rhs_static(c) - rhs_vely(i,j)
                end if
            end do
        end do
        deallocate(rhs_vely)

        allocate(pres(1:localncols,1:localnrows))
        pres(:,:) = 0
        allocate(rhs_dpdt(1:localncols,1:localnrows))
        call Resi_dpdt(pres, rhs_dpdt)
        allocate(resiAcp(1:localncols,1:localnrows))
        allocate(resitemp(1:localncols,1:localnrows))
        do j = 1, 3
            do i = 1, 3
                call genExpField(i, j, localncols, localnrows, pres, isField)
                if(isField) then
                    call Resi_dpdt(pres, resiAcp)
                    resitemp = resiAcp - rhs_dpdt
                    call constructAcp(pres, resitemp)
                end if
            end do
        end do
        deallocate(pres)
        deallocate(resiAcp)
        deallocate(resitemp)

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_rhs(c) = local_rhs_static(c) - rhs_dpdt(i,j) - &!
                    kc(i,j)*ks(i,j)*Cf(i,j)*av(i,j)*al/rhos/(kc(i,j)+ks(i,j))
            end do
        end do
        deallocate(rhs_dpdt)

    end subroutine genDynPara_vp

    subroutine computevp()

        integer :: indexr, indexu
        integer :: nVelx, nVely, nPres
        integer :: AxxBeInd, AxpBeInd, AyyBeInd, AypBeInd, AcxBeInd, AcyBeInd, AcpBeInd
        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), pointer :: local_x
        integer :: status(MPI_STATUS_SIZE)
        integer :: request, requestl, requestr, requestd, requestu, requestld, requestlu, requestrd
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        real(kind=8) :: solvertimestart, solvertimefinish
        integer :: i, j, l, n, c, num_iter, ierr

        allocate(cols(7))
        allocate(values(7)) ! Each row of A has at most 7 nonzero values.
        allocate(local_x(local_x_size))

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
                call MPI_RECV(slave_data, slave_vp_data_size(i), MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
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

        if(pcol == pncols) then
            nVelx = (localncols+1)*localnrows
        else
            nVelx = localncols*localnrows
        end if
        if(prow == pnrows) then
            nVely = localncols*(localnrows+1)
        else
            nVely = localncols*localnrows
        end if
        nPres = localncols*localnrows

        AxxBeInd = 1
        AxpBeInd = 1
        AyyBeInd = 1
        AypBeInd = 1
        AcxBeInd = 1
        AcyBeInd = 1
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

            ! the line is in the continuity part
            elseif(n >= ilower+nVelx+nVely) then

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
        indexr = localncols
        if(pcol == pncols) then
            indexr = localncols + 1
        end if
        do j = 1, localnrows
            do i = 1, indexr
                c = c + 1
                vx(i,j) = local_x(c)
            end do
        end do

        indexu = localnrows
        if(prow == pnrows) then
            indexu = localnrows + 1
        end if
        do j = 1, indexu
            do i = 1, localncols
                c = c + 1
                vy(i,j) = local_x(c)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                p(i,j) = local_x(c)
            end do
        end do

        if(pcol /= 1) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            sentSize = localnrows + indexu
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = vx(1,j)
            end do
            do j = 1, indexu
                c = c + 1
                sent(c) = vy(1,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            sentSize = indexu
            allocate(sent(sentSize))
            c = 0
            do j = 1, indexu
                c = c + 1
                sent(c) = vy(localncols,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            sentSize = indexr + localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, indexr
                c = c + 1
                sent(c) = vx(i,1)
            end do
            do i = 1, localncols
                c = c + 1
                sent(c) = vy(i,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            sentSize = indexr
            allocate(sent(sentSize))
            c = 0
            do i = 1, indexr
                c = c + 1
                sent(c) = vx(i,localnrows)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=1)) then
            allocate(sent(2))
            sent(1) = vx(1,1)
            sent(2) = vy(1,1)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid-1-pncols, myid, MPI_COMM_WORLD, requestld, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(sent(1))
            sent(1) = vx(1,localnrows)
            call MPI_IBSEND(sent, 1, MPI_DOUBLE_PRECISION, myid-1+pncols, myid, MPI_COMM_WORLD, requestlu, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(sent(1))
            sent(1) = vy(localncols,1)
            call MPI_IBSEND(sent, 1, MPI_DOUBLE_PRECISION, myid+1-pncols, myid, MPI_COMM_WORLD, requestrd, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            recvSize = localnrows + indexu
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                vx(localncols+1,j) = recv(c)
            end do
            do j = 1, indexu
                c = c + 1
                vy(localncols+1,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            recvSize = indexu
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, indexu
                c = c + 1
                vy(0,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            recvSize = indexr + localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, indexr
                c = c + 1
                vx(i,localnrows+1) = recv(c)
            end do
            do i = 1, localncols
                c = c + 1
                vy(i,localnrows+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            recvSize = indexr
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, indexr
                c = c + 1
                vx(i,0) = recv(c)
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=pnrows)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid+1+pncols, myid+1+pncols, MPI_COMM_WORLD, status, ierr)
            vx(localncols+1,localnrows+1) = recv(1)
            vy(localncols+1,localnrows+1) = recv(2)
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(recv(1))
            call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, myid+1-pncols, myid+1-pncols, MPI_COMM_WORLD, status, ierr)
            vx(localncols+1,0) = recv(1)
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(recv(1))
            call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, myid-1+pncols, myid-1+pncols, MPI_COMM_WORLD, status, ierr)
            vy(0,localnrows+1) = recv(1)
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
        if((pcol/=1).and.(prow/=1)) then
            call MPI_WAIT(requestld, status, ierr)
        end if
        if((pcol/=1).and.(prow/=pnrows)) then
            call MPI_WAIT(requestlu, status, ierr)
        end if
        if((pcol/=pncols).and.(prow/=1)) then
            call MPI_WAIT(requestrd, status, ierr)
        end if

        deallocate(cols)
        deallocate(values)
        deallocate(local_x)

    end subroutine computevp

    subroutine genDynPara_Cf()

        integer :: findexl, findexr, findexd, findexu
        integer :: eindexl, eindexr, eindexd, eindexu
        integer, dimension(:,:), pointer :: conc
        logical :: isField
        real(kind=8), dimension(:,:), pointer :: rhs_Cf
        real(kind=8), dimension(:,:), pointer :: resiAcf, resitemp
        integer :: i, j, c

        findexl = 0
        findexr = localncols + 1
        findexd = 0
        findexu = localnrows + 1
        allocate(conc(findexl:findexr,findexd:findexu))
        conc(:,:) = 0

        eindexl = 1
        eindexr = localncols
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_Cf(eindexl:eindexr,eindexd:eindexu))
        rhs_Cf(:,:) = 0.D0

        call Resi_Cf(conc, rhs_Cf)

        allocate(resiAcf(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, conc, isField)
                if(isField) then
                    call Resi_Cf(conc, resiAcf)
                    resitemp = resiAcf - rhs_Cf
                    call constructAcf(conc, resitemp)
                end if
            end do
        end do

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_rhs_Cf(c) = -rhs_Cf(i,j)
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
        integer :: i, j, l, n, c, num_iter, ierr

        allocate(cols(9))
        allocate(values(9)) ! Each row of A has at most 9 nonzero values.
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
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                Cf(i,j) = local_x(c)
            end do
        end do

        deallocate(cols)
        deallocate(values) 
        deallocate(local_x)

    end subroutine computeCf

    subroutine genDynPara_Tem()

        integer :: findexl, findexr, findexd, findexu
        integer :: eindexl, eindexr, eindexd, eindexu
        integer, dimension(:,:), pointer :: tempe
        logical :: isField
        real(kind=8), dimension(:,:), pointer :: rhs_Tem
        real(kind=8), dimension(:,:), pointer :: resiAtem, resitemp
        integer :: i, j, c

        findexl = 0
        findexr = localncols + 1
        findexd = 0
        findexu = localnrows + 1
        allocate(tempe(findexl:findexr,findexd:findexu))
        tempe(:,:) = 0

        eindexl = 1
        eindexr = localncols
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_Tem(eindexl:eindexr,eindexd:eindexu))
        rhs_Tem(:,:) = 0.D0

        call Resi_Tem(tempe, rhs_Tem)

        allocate(resiAtem(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, tempe, isField)
                if(isField) then
                    call Resi_Tem(tempe, resiAtem)
                    resitemp = resiAtem - rhs_Tem
                    call constructAtem(tempe, resitemp)
                end if
            end do
        end do

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_rhs_Tem(c) = -rhs_Tem(i,j)
            end do
        end do

        deallocate(tempe)
        deallocate(rhs_Tem)
        deallocate(resiAtem)
        deallocate(resitemp)

    end subroutine genDynPara_Tem

    subroutine computeTem()

        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), pointer :: local_x
        integer :: AtemBeInd
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data
        real(kind=8) :: solvertimestart, solvertimefinish
        integer :: i, j, l, n, c, num_iter, ierr

        allocate(cols(5))
        allocate(values(5)) ! Each row of A has at most 5 nonzero values.
        allocate(local_x(local_x_size_Tem))

#ifdef LAPACK

        A_lapack_Tem(:,:) = 0.D0
        b_lapack_Tem(:) = local_rhs_Tem(:)
        IPIV_Tem(:) = 0.D0

#elif defined(UMFPACK)

        Ap_Tem(1) = 0
        Ai_Tem(:) = 0
        Ax_Tem(:) = 0.D0

#elif defined(MUMPS)

        mumps_NNZ_loc_Tem= 0
        mumps_IRN_loc_Tem(:) = 0
        mumps_JCN_loc_Tem(:) = 0
        mumps_A_loc_Tem(:) = 0.D0

        if(mumps_par_Tem%MYID /= 0) then
            call MPI_IBSEND(local_rhs_Tem, local_x_size_Tem, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(mumps_par_Tem%MYID == 0) then
            mumps_par_Tem%RHS(1:local_x_size_Tem) = local_rhs_Tem(1:local_x_size_Tem)
            c = local_x_size_Tem + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Tem))
                call MPI_RECV(slave_data, local_x_size_Tem, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    mumps_par_Tem%RHS(c:c+local_x_size_Tem-1) = slave_data(:)
                c = c + local_x_size_Tem
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Tem%MYID /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

#elif defined(HYPRE)

        call HYPRE_IJVectorSetValues(b_Tem, local_x_size_Tem, rows_Tem, local_rhs_Tem, ierr)
        call HYPRE_IJVectorSetValues(x_Tem, local_x_size_Tem, rows_Tem, initial_x_guess_Tem, ierr)

        call HYPRE_IJVectorAssemble(b_Tem, ierr)
        call HYPRE_IJVectorAssemble(x_Tem, ierr)

        call HYPRE_IJVectorGetObject(b_Tem, par_b_Tem, ierr)
        call HYPRE_IJVectorGetObject(x_Tem, par_x_Tem, ierr)

#endif

        AtemBeInd = 1
        do n = ilower_Tem, iupper_Tem

            cols(:) = 0
            values(:) = 0.D0

            c = 0
            do l = AtemBeInd, AtemSize
                if(AtemRows(l) == n) then
                    c = c + 1
                    cols(c) = AtemCols(l)
                    values(c) = AtemValues(l)
                else
                    AtemBeInd = l
                    exit
                end if
            end do

#ifdef LAPACK
            do l = 1, c
                A_lapack_Tem(n,cols(l)) = values(l)
            end do
#elif defined(UMFPACK)
            Ap_Tem(n+1) = Ap_Tem(n) + c
            Ai_Tem(Ap_Tem(n)+1:Ap_Tem(n+1)) = cols(1:c) - 1
            Ax_Tem(Ap_Tem(n)+1:Ap_Tem(n+1)) = values(1:c)
#elif defined(MUMPS)
            mumps_IRN_loc_Tem(mumps_NNZ_loc_Tem+1:mumps_NNZ_loc_Tem+c) = n
            mumps_JCN_loc_Tem(mumps_NNZ_loc_Tem+1:mumps_NNZ_loc_Tem+c) = cols(1:c)
            mumps_A_loc_Tem(mumps_NNZ_loc_Tem+1:mumps_NNZ_loc_Tem+c) = values(1:c)
            mumps_NNZ_loc_Tem = mumps_NNZ_loc_Tem + c
#elif defined(HYPRE)
            call HYPRE_IJMatrixSetValues(A_Tem, 1, c, n, cols, values, ierr)
#endif

        end do

#ifdef LAPACK

        ! you have to make sure that the number of processors is set to 1 when using such method.
        solvertimestart = MPI_Wtime()
        call dgesv(local_x_size_Tem, 1, A_lapack_Tem, local_x_size_Tem, IPIV_Tem, b_lapack_Tem, local_x_size_Tem, LAPACKINFO)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        local_x(:) = b_lapack_Tem(:)
        if(LAPACKINFO /= 0) then
            print *, 'LAPACK solver error when solving Tem. INFO = ', LAPACKINFO
            stop
        end if

#elif defined(UMFPACK)

        solvertimestart = MPI_Wtime()
        call umf4def(control)
        call umf4sym(local_x_size_Tem, local_x_size_Tem, Ap_Tem, Ai_Tem, Ax_Tem, symbolic, control, umfinfo)
        call umf4num(Ap_Tem, Ai_Tem, Ax_Tem, symbolic, numeric, control, umfinfo)
        call umf4fsym(symbolic)
        call umf4sol(1, local_x, local_rhs_Tem, numeric, control, umfinfo) ! 1 means A'x=b
        if(umfinfo(1) < 0) then
            print *, 'UMFPACK solver error when solving Tem. Info: ', umfinfo(1)
            stop
        end if
        call umf4fnum(numeric)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart

#elif defined(MUMPS)

        mumps_par_Tem%NNZ_loc = mumps_NNZ_loc_Tem
        allocate(mumps_par_Tem%IRN_loc(mumps_par_Tem%NNZ_loc))
        allocate(mumps_par_Tem%JCN_loc(mumps_par_Tem%NNZ_loc))
        allocate(mumps_par_Tem%A_loc(mumps_par_Tem%NNZ_loc))
        mumps_par_Tem%IRN_loc(1:mumps_par_Tem%NNZ_loc) = mumps_IRN_loc_Tem(1:mumps_NNZ_loc_Tem)
        mumps_par_Tem%JCN_loc(1:mumps_par_Tem%NNZ_loc) = mumps_JCN_loc_Tem(1:mumps_NNZ_loc_Tem)
        mumps_par_Tem%A_loc(1:mumps_par_Tem%NNZ_loc) = mumps_A_loc_Tem(1:mumps_NNZ_loc_Tem)

        mumps_par_Tem%JOB = 6
        solvertimestart = MPI_Wtime()
        call DMUMPS(mumps_par_Tem)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(mumps_par_Tem%INFOG(1) < 0) then
            stop
        end if

        if(mumps_par_Tem%MYID == 0) then
            local_x(:) = mumps_par_Tem%RHS(1:local_x_size_Tem)
            allocate(requestarray(nProcs-1))
            c = local_x_size_Tem + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Tem))
                slave_data(:) = mumps_par_Tem%RHS(c:c+local_x_size_Tem-1)
                call MPI_IBSEND(slave_data, local_x_size_Tem, MPI_DOUBLE_PRECISION, i, myid, &!
                    MPI_COMM_WORLD, requestarray(i), ierr)
                c = c + local_x_size_Tem
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Tem%MYID /= 0) then
            allocate(slave_data(local_x_size_Tem))
            call MPI_RECV(slave_data, local_x_size_Tem, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
            local_x(:) = slave_data(:)
            deallocate(slave_data)
        end if

        if(mumps_par_Tem%MYID == 0) then
            do i = 1, nProcs-1
                call MPI_WAIT(requestarray(i), status, ierr)
            end do
            deallocate(requestarray)
        end if

        deallocate(mumps_par_Tem%IRN_loc)
        deallocate(mumps_par_Tem%JCN_loc)
        deallocate(mumps_par_Tem%A_loc)

#elif defined(HYPRE)

        call HYPRE_IJMatrixAssemble(A_Tem, ierr)
        call HYPRE_IJMatrixGetObject(A_Tem, parcsr_A_Tem, ierr)

        solvertimestart = MPI_Wtime()
        call HYPRE_ParCSRGMRESSetup(solver_Tem, parcsr_A_Tem, par_b_Tem, par_x_Tem, ierr)
        call HYPRE_ParCSRGMRESSolve(solver_Tem, parcsr_A_Tem, par_b_Tem, par_x_Tem, ierr)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        call HYPRE_ParCSRGMRESGetNumIteratio(solver_Tem, num_iter, ierr)
        if(ierr /= 0) then
            if(myid == 0) then
                print *, 'HYPRE solver error when solving Tem. ierr = ', ierr
                stop
            end if
        end if

        call HYPRE_IJVectorGetValues(x_Tem, local_x_size_Tem, rows_Tem, local_x, ierr)

        ! let the solution of this time step be the initial x guess in the next time step.
        ! by this way, the number of solver iteration steps can be reduced greatly.
        initial_x_guess_Tem(:) = local_x(:)

#endif

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                Tem(i,j) = local_x(c)
            end do
        end do

        deallocate(cols)
        deallocate(values)
        deallocate(local_x)

    end subroutine computeTem

    subroutine outputHisData()

        real(kind=8), dimension(:), allocatable :: local_data, recv
        real(kind=8), dimension(:), pointer :: global_data
        integer :: local_data_size
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        real(kind=8) :: poroavg, Kxxavg, avavg, pavg, Cfavg, Temavg, qavg
        real(kind=8) :: localqsum, localpsum
        integer :: i, j, c, ierr

        local_data_size = 6
        allocate(local_data(local_data_size))
        local_data(:) = 0.D0
        do j = 1, localnrows
            do i = 1, localncols
                local_data(1) = local_data(1) + poro(i,j)
                local_data(2) = local_data(2) + Kxx(i,j)
                local_data(3) = local_data(3) + av(i,j)
                local_data(4) = local_data(4) + p(i,j)
                local_data(5) = local_data(5) + Cf(i,j)
                local_data(6) = local_data(6) + Tem(i,j)
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
            poroavg = poroavg/(nx*ny)

            Kxxavg = 0.D0
            do c = 2, local_data_size*nProcs, local_data_size
                Kxxavg = Kxxavg + global_data(c)
            end do
            Kxxavg = Kxxavg/(nx*ny)

            avavg = 0.D0
            do c = 3, local_data_size*nProcs, local_data_size
                avavg = avavg + global_data(c)
            end do
            avavg = avavg/(nx*ny)

            pavg = 0.D0
            do c = 4, local_data_size*nProcs, local_data_size
                pavg = pavg + global_data(c)
            end do
            pavg = pavg/(nx*ny)

            Cfavg = 0.D0
            do c = 5, local_data_size*nProcs, local_data_size
                Cfavg = Cfavg + global_data(c)
            end do
            Cfavg = Cfavg/(nx*ny)

            Temavg = 0.D0
            do c = 6, local_data_size*nProcs, local_data_size
                Temavg = Temavg + global_data(c)
            end do
            Temavg = Temavg/(nx*ny)

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
        if(pcol == pncols) then

            localqsum = 0.D0
            do j = 1, localnrows
                localqsum = localqsum + vx(localncols+1,j)
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localqsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pnrows))
            allocate(recv(1))

            if(pcol == pncols) then
                global_data(1) = localqsum
                c = 1
                do i = pncols-1, nProcs-1, pncols
                    if(i > 0) then
                        call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                        c = c + 1
                        global_data(c) = recv(1)
                    end if
                end do
            else
                c = 0
                do i = pncols-1, nProcs-1, pncols
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end do
            end if

            qavg = 0.D0
            do c = 1, pnrows
                qavg = qavg + global_data(c)
            end do
            qavg = qavg/ny

            write(46, fmt='(es24.16)', iostat=ierr) qavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if

            deallocate(global_data)
            deallocate(recv)

        end if

        if((pcol==pncols).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

        ! output p at the entry
        if(pcol == 1) then

            localpsum = 0.D0
            do j = 1, localnrows
                localpsum = localpsum + p(1,j)
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pnrows))
            allocate(recv(1))

            global_data(1) = localpsum
            c = 1
            do i = 0, nProcs-1, pncols
                if(i > 0) then
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end if
            end do

            pavg = 0.D0
            do c = 1, pnrows
                pavg = pavg + global_data(c)
            end do
            pavg = pavg/ny

            write(47, fmt='(es24.16)', iostat=ierr) pavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if

            deallocate(global_data)
            deallocate(recv)

        end if

        if((pcol==1).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

    end subroutine outputHisData

    subroutine computePresDrop(isBT)

        logical, intent(out) :: isBT
        real(kind=8) :: inpavg, outpavg, localpsum, presDrop
        real(kind=8), dimension(:), allocatable :: global_data, recv
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer :: i, j, c, ierr

        if(pcol == 1) then

            localpsum = 0.D0
            do j = 1, localnrows
                localpsum = localpsum + p(1,j)
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pnrows))
            allocate(recv(1))

            global_data(1) = localpsum
            c = 1
            do i = 0, nProcs-1, pncols
                if(i > 0) then
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end if
            end do

            inpavg = 0.D0
            do c = 1, pnrows
                inpavg = inpavg + global_data(c)
            end do
            inpavg = inpavg/ny

            deallocate(global_data)
            deallocate(recv)

        end if

        if((pcol==1).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

        if(pcol == pncols) then

            localpsum = 0.D0
            do j = 1, localnrows
                localpsum = localpsum + p(localncols,j)
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pnrows))
            allocate(recv(1))

            if(pcol == pncols) then
                global_data(1) = localpsum
                c = 1
                do i = pncols-1, nProcs-1, pncols
                    if(i > 0) then
                        call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                        c = c + 1
                        global_data(c) = recv(1)
                    end if
                end do
            else
                c = 0
                do i = pncols-1, nProcs-1, pncols
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end do
            end if

            outpavg = 0.D0
            do c = 1, pnrows
                outpavg = outpavg + global_data(c)
            end do
            outpavg = outpavg/ny

            deallocate(global_data)
            deallocate(recv)

        end if

        if((pcol==pncols).and.(myid/=0)) then
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
                print *, 'The pore volume to breakthrough is ', (t-1)*timeEnd/nt*abs(vxBdryX0(1))
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
        integer :: indexr, indexu
        integer :: i, j, c, ierr

        local_data_size = local_x_size+4*localncols*localnrows
        allocate(local_data(local_data_size))

        c = 0

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = poro(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = Kxx(i,j)
            end do
        end do

        if(pcol == pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        do j = 1, localnrows
            do i = 1, indexr
                c = c + 1
                local_data(c) = vx(i,j)
            end do
        end do

        if(prow == pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if
        do j = 1, indexu
            do i = 1, localncols
                c = c + 1
                local_data(c) = vy(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = p(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = Cf(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = Tem(i,j)
            end do
        end do

        if(myid /= 0) then
            call MPI_IBSEND(local_data, local_data_size, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(myid == 0) then

            allocate(global_data((nx+1)*ny+nx*(ny+1)+5*nx*ny))

            global_data(1:local_data_size) = local_data(1:local_data_size)
            c = local_data_size + 1
            do i = 1, nProcs-1
                slave_data_size = slave_vp_data_size(i)+4*localncols*localnrows
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
        deallocate(ts)
        deallocate(src)
        deallocate(poroInit)
        deallocate(KxxInit)
        deallocate(KyyInit)
        deallocate(avInit)
        deallocate(vxBdryX0)
        deallocate(vxBdryX1)
        deallocate(vxBdryY0)
        deallocate(vxBdryY1)
        deallocate(vyBdryX0)
        deallocate(vyBdryX1)
        deallocate(vyBdryY0)
        deallocate(vyBdryY1)
        deallocate(isDiriX0_p)
        deallocate(isDiriX1_p)
        deallocate(isDiriY0_p)
        deallocate(isDiriY1_p)
        deallocate(pBdryX0)
        deallocate(pBdryX1)
        deallocate(pBdryY0)
        deallocate(pBdryY1)
        deallocate(pInit)
        deallocate(isDiriX0_Cf)
        deallocate(isDiriX1_Cf)
        deallocate(isDiriY0_Cf)
        deallocate(isDiriY1_Cf)
        deallocate(CfBdryX0)
        deallocate(CfBdryX1)
        deallocate(CfBdryY0)
        deallocate(CfBdryY1)
        deallocate(CfInit)
        deallocate(isDiriX0_Tem)
        deallocate(isDiriX1_Tem)
        deallocate(isDiriY0_Tem)
        deallocate(isDiriY1_Tem)
        deallocate(TemBdryX0)
        deallocate(TemBdryX1)
        deallocate(TemBdryY0)
        deallocate(TemBdryY1)
        deallocate(TemInit)

        deallocate(dm)
        deallocate(kc)
        deallocate(ks)
        deallocate(hx)
        deallocate(hy)
        deallocate(poro)
        deallocate(poro_old)
        deallocate(poroHarmX)
        deallocate(poroHarmX_old)
        deallocate(poroHarmXInit)
        deallocate(poroHarmY)
        deallocate(poroHarmY_old)
        deallocate(poroHarmYInit)
        deallocate(Kxx)
        deallocate(Kyy)
        deallocate(KxxHarm)
        deallocate(KyyHarm)
        deallocate(av)
        deallocate(vx)
        deallocate(vy)
        deallocate(p)
        deallocate(Cf)
        deallocate(Tem)
        
        deallocate(local_rhs)
        deallocate(local_rhs_static)
        deallocate(local_rhs_Cf)
        deallocate(local_rhs_Tem)
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
        deallocate(AcxCols)
        deallocate(AcxRows)
        deallocate(AcxValues)
        deallocate(AcyCols)
        deallocate(AcyRows)
        deallocate(AcyValues)
        deallocate(AcpCols)
        deallocate(AcpRows)
        deallocate(AcpValues)
        deallocate(AcfCols)
        deallocate(AcfRows)
        deallocate(AcfValues)
        deallocate(AtemCols)
        deallocate(AtemRows)
        deallocate(AtemValues)
        deallocate(AxxEntryNum)
        deallocate(AxpEntryNum)
        deallocate(AyyEntryNum)
        deallocate(AypEntryNum)
        deallocate(AcxEntryNum)
        deallocate(AcyEntryNum)
        deallocate(AcpEntryNum)
        deallocate(AcfEntryNum)
        deallocate(AtemEntryNum)
        deallocate(AxxEntryBase)
        deallocate(AxpEntryBase)
        deallocate(AyyEntryBase)
        deallocate(AypEntryBase)
        deallocate(AcxEntryBase)
        deallocate(AcyEntryBase)
        deallocate(AcpEntryBase)
        deallocate(AcfEntryBase)
        deallocate(AtemEntryBase)

        if((nProcs>1).and.(myid==0)) then
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
            close(47)
        end if

#ifdef LAPACK

        deallocate(A_lapack)
        deallocate(b_lapack)
        deallocate(IPIV)
        deallocate(A_lapack_Cf)
        deallocate(b_lapack_Cf)
        deallocate(IPIV_Cf)
        deallocate(A_lapack_Tem)
        deallocate(b_lapack_Tem)
        deallocate(IPIV_Tem)

#elif defined(UMFPACK)

        deallocate(Ap)
        deallocate(Ai)
        deallocate(Ax)
        deallocate(Ap_Cf)
        deallocate(Ai_Cf)
        deallocate(Ax_Cf)
        deallocate(Ap_Tem)
        deallocate(Ai_Tem)
        deallocate(Ax_Tem)

#elif defined(MUMPS)

        deallocate(mumps_par%RHS)
        deallocate(mumps_par_Cf%RHS)
        deallocate(mumps_par_Tem%RHS)

        mumps_par%JOB = -2
        call DMUMPS(mumps_par)
        mumps_par_Cf%JOB = -2
        call DMUMPS(mumps_par_Cf)
        mumps_par_Tem%JOB = -2
        call DMUMPS(mumps_par_Tem)

        deallocate(mumps_IRN_loc)
        deallocate(mumps_JCN_loc)
        deallocate(mumps_A_loc)
        deallocate(mumps_IRN_loc_Cf)
        deallocate(mumps_JCN_loc_Cf)
        deallocate(mumps_A_loc_Cf)
        deallocate(mumps_IRN_loc_Tem)
        deallocate(mumps_JCN_loc_Tem)
        deallocate(mumps_A_loc_Tem)

#elif defined(HYPRE)

        call HYPRE_ParaSailsDestroy(precond, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver, ierr)
        call HYPRE_ParaSailsDestroy(precond_Cf, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver_Cf, ierr)
        call HYPRE_ParaSailsDestroy(precond_Tem, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver_Tem, ierr)

        call HYPRE_IJMatrixDestroy(A, ierr)
        call HYPRE_IJVectorDestroy(b, ierr)
        call HYPRE_IJVectorDestroy(x, ierr)
        call HYPRE_IJMatrixDestroy(A_Cf, ierr)
        call HYPRE_IJVectorDestroy(b_Cf, ierr)
        call HYPRE_IJVectorDestroy(x_Cf, ierr)
        call HYPRE_IJMatrixDestroy(A_Tem, ierr)
        call HYPRE_IJVectorDestroy(b_Tem, ierr)
        call HYPRE_IJVectorDestroy(x_Tem, ierr)

        deallocate(rows)
        deallocate(rows_Cf)
        deallocate(rows_Tem)
        deallocate(initial_x_guess)
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

        call genStaticPara_vp()

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

            call computekc()
            call genDynPara_vp()
            call computevp()

            call computekc()
            call genDynPara_Cf()
            call computeCf()

            call genDynPara_Tem()
            call computeTem()

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
