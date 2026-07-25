
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

program DBF_infile

    use DBF_model
    use DBF_driver
    use RST_proceAlloc

    implicit none
    type(model) :: modelCase
    real(kind=8) :: rand
    real(kind=8), dimension(:), allocatable :: random
    integer :: i, j, k, times, ierr

    times = 1  !!!!
    modelCase%Lx = 1.D-1
    modelCase%Ly = 4.D-2
    modelCase%nx = 100!180*times
    modelCase%ny = 40!72*times
    modelCase%timeEnd = 9.D2 !!!!
    modelCase%nt = 9.D2 !!!!
    modelCase%visc = 1.D-3
    modelCase%rhof = 1.D3!1.01D3
    modelCase%rhos = 2.7D3!2.71D3
    !modelCase%alpha = 1.D0
    modelCase%dmRef = 3.6D-9
    modelCase%TemRef_dm = 2.98D2
    modelCase%ksRef = 2.D-3
    modelCase%TemRef_ks = 2.98D2
    modelCase%epslon = 0.D0
    modelCase%alphaOS = 5.D-1
    modelCase%lamdaX = 5.D-1
    modelCase%lamdaT = 1.D-1
    modelCase%lamdaf = 5.8D-1
    modelCase%lamdas = 5.526D0
    modelCase%thetaf = 4.184D3
    modelCase%thetas = 2.D2
    modelCase%radiusInit = 1.D-6
    modelCase%ShInfinity = 3.D0!3.66D0
    modelCase%al = 1.3D0!5.D-2
    modelCase%gravX = 0.D0
    modelCase%gravY = 0.D0!-9.807

    allocate(modelCase%xs(modelCase%nx+1))
    do i = 1, modelCase%nx+1
        modelCase%xs(i) = (i-1)*modelCase%Lx/modelCase%nx
    end do
    allocate(modelCase%ys(modelCase%ny+1))
    do j = 1, modelCase%ny+1
        modelCase%ys(j) = (j-1)*modelCase%Ly/modelCase%ny
    end do
    allocate(modelCase%ts(modelCase%nt+1))
    do k = 1, modelCase%nt+1
        modelCase%ts(k) = (k-1)*modelCase%timeEnd/modelCase%nt
    end do

    allocate(modelCase%src(modelCase%nx,modelCase%ny))
    modelCase%src(:,:) = 0.D0

    !call genRandomNum()
    allocate(random(RANDOMSIZE))
    open(unit=10, file=trim(adjustl(FRANDOMTXT)), iostat=ierr)
    if(ierr /= 0) then
        print *, 'open file error. ', ierr
        stop
    end if
    read(10, fmt="(f8.6)") random(:)
    close(10)
    allocate(modelCase%poroInit(modelCase%nx,modelCase%ny))
    k = 0
    do j = 1, modelCase%ny
        do i = 1, modelCase%nx
            k = k + 1
            if(k <= RANDOMSIZE) then
                modelCase%poroInit(i,j) = random(k)
            else
                print *, 'The random numbers are not enough!'
                print *, 'Please regenerate the random number file.'
                stop
            end if
        end do
    end do
    deallocate(random)
    ! add the fracture
    do i = 11, 20
        j = 52 - 2*i
        modelCase%poroInit(i,j) = 0.7
        modelCase%poroInit(i,j-1) = 0.7
    end do
    do i = 31, 40
        j = 2*i - 51
        modelCase%poroInit(i,j) = 0.7
        modelCase%poroInit(i,j+1) = 0.7
    end do
    do i = 61, 70
        j = 152 - 2*i
        modelCase%poroInit(i,j) = 0.7
        modelCase%poroInit(i,j-1) = 0.7
    end do
    do i = 81, 90
        j = 2*i - 151
        modelCase%poroInit(i,j) = 0.7
        modelCase%poroInit(i,j+1) = 0.7
    end do

    allocate(modelCase%KxxInit(modelCase%nx,modelCase%ny))
    modelCase%KxxInit(:,:) = 1.D-15!9.869233D-16
    allocate(modelCase%KyyInit(modelCase%nx,modelCase%ny))
    modelCase%KyyInit(:,:) = 1.D-15!9.869233D-16
    ! add the fracture
    do i = 11, 20
        j = 52 - 2*i
        modelCase%KxxInit(i,j) = 1.D-12
        modelCase%KxxInit(i,j-1) = 1.D-12
        modelCase%KyyInit(i,j) = 1.D-12
        modelCase%KyyInit(i,j-1) = 1.D-12
    end do
    do i = 31, 40
        j = 2*i - 51
        modelCase%KxxInit(i,j) = 1.D-12
        modelCase%KxxInit(i,j+1) = 1.D-12
        modelCase%KyyInit(i,j) = 1.D-12
        modelCase%KyyInit(i,j+1) = 1.D-12
    end do
    do i = 61, 70
        j = 152 - 2*i
        modelCase%KxxInit(i,j) = 1.D-12
        modelCase%KxxInit(i,j-1) = 1.D-12
        modelCase%KyyInit(i,j) = 1.D-12
        modelCase%KyyInit(i,j-1) = 1.D-12
    end do
    do i = 81, 90
        j = 2*i - 151
        modelCase%KxxInit(i,j) = 1.D-12
        modelCase%KxxInit(i,j+1) = 1.D-12
        modelCase%KyyInit(i,j) = 1.D-12
        modelCase%KyyInit(i,j+1) = 1.D-12
    end do

    allocate(modelCase%avInit(modelCase%nx,modelCase%ny))
    modelCase%avInit(:,:) = 1.8D2!5.D-1

    allocate(modelCase%vxBdryX0(modelCase%ny))
    modelCase%vxBdryX0(:) = -1.D-3!-4.17D-6!!!!
    allocate(modelCase%vxBdryX1(modelCase%ny))
    modelCase%vxBdryX1(:) = 0.D0
    allocate(modelCase%vyBdryX0(modelCase%ny+1))
    modelCase%vyBdryX0(:) = 0.D0
    allocate(modelCase%vyBdryX1(modelCase%ny+1))
    modelCase%vyBdryX1(:) = 0.D0
    allocate(modelCase%vxBdryY0(modelCase%nx+1))
    modelCase%vxBdryY0(:) = 0.D0
    allocate(modelCase%vxBdryY1(modelCase%nx+1))
    modelCase%vxBdryY1(:) = 0.D0
    allocate(modelCase%vyBdryY0(modelCase%nx))
    modelCase%vyBdryY0(:) = 0.D0
    allocate(modelCase%vyBdryY1(modelCase%nx))
    modelCase%vyBdryY1(:) = 0.D0

    allocate(modelCase%isDiriX0_p(modelCase%ny))
    modelCase%isDiriX0_p(:) = 0
    allocate(modelCase%isDiriX1_p(modelCase%ny))
    modelCase%isDiriX1_p(:) = 1
    allocate(modelCase%isDiriY0_p(modelCase%nx))
    modelCase%isDiriY0_p(:) = 0
    allocate(modelCase%isDiriY1_p(modelCase%nx))
    modelCase%isDiriY1_p(:) = 0

    allocate(modelCase%pBdryX0(modelCase%ny))
    modelCase%pBdryX0(:) = 0.D0
    allocate(modelCase%pBdryX1(modelCase%ny))
    modelCase%pBdryX1(:) = 1.D5!1.52D7
    allocate(modelCase%pBdryY0(modelCase%nx))
    modelCase%pBdryY0(:) = 0.D0
    allocate(modelCase%pBdryY1(modelCase%nx))
    modelCase%pBdryY1(:) = 0.D0
    allocate(modelCase%pInit(modelCase%nx,modelCase%ny))
    modelCase%pInit(:,:) = 1.D5!1.52D7

    allocate(modelCase%isDiriX0_Cf(modelCase%ny))
    modelCase%isDiriX0_Cf(:) = 1
    allocate(modelCase%isDiriX1_Cf(modelCase%ny))
    modelCase%isDiriX1_Cf(:) = 0
    allocate(modelCase%isDiriY0_Cf(modelCase%nx))
    modelCase%isDiriY0_Cf(:) = 0
    allocate(modelCase%isDiriY1_Cf(modelCase%nx))
    modelCase%isDiriY1_Cf(:) = 0

    allocate(modelCase%CfBdryX0(modelCase%ny))
    modelCase%CfBdryX0(:) = 2.17D2!5.D2
    allocate(modelCase%CfBdryX1(modelCase%ny))
    modelCase%CfBdryX1(:) = 0.D0
    allocate(modelCase%CfBdryY0(modelCase%nx))
    modelCase%CfBdryY0(:) = 0.D0
    allocate(modelCase%CfBdryY1(modelCase%nx))
    modelCase%CfBdryY1(:) = 0.D0
    allocate(modelCase%CfInit(modelCase%nx,modelCase%ny))
    modelCase%CfInit(:,:) = 0.D0

    allocate(modelCase%isDiriX0_Tem(modelCase%ny))
    modelCase%isDiriX0_Tem(:) = 1
    allocate(modelCase%isDiriX1_Tem(modelCase%ny))
    modelCase%isDiriX1_Tem(:) = 0
    allocate(modelCase%isDiriY0_Tem(modelCase%nx))
    modelCase%isDiriY0_Tem(:) = 0
    allocate(modelCase%isDiriY1_Tem(modelCase%nx))
    modelCase%isDiriY1_Tem(:) = 0

    allocate(modelCase%TemBdryX0(modelCase%ny))
    modelCase%TemBdryX0(:) = 2.95D2
    allocate(modelCase%TemBdryX1(modelCase%ny))
    modelCase%TemBdryX1(:) = 0.D0
    allocate(modelCase%TemBdryY0(modelCase%nx))
    modelCase%TemBdryY0(:) = 0.D0
    allocate(modelCase%TemBdryY1(modelCase%nx))
    modelCase%TemBdryY1(:) = 0.D0
    allocate(modelCase%TemInit(modelCase%nx,modelCase%ny))
    modelCase%TemInit(:,:) = 2.95D2

    modelCase%soludoc = 'case1'

    call proceAlloc(4, modelCase%nx, modelCase%ny, modelCase%pncols, modelCase%pnrows)

    call driver(modelCase)

end program DBF_infile


The breakthrough time is    278.00000000000000       seconds.
The pore volume to breakthrough is   0.27800000000000002
Solver time =    95.807069157599472       seconds.
Elapsed time =    100.11413053295109       seconds.
