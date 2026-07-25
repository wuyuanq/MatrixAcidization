
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
    modelCase%nx = 100*times
    modelCase%ny = 40*times
    modelCase%timeEnd = 2.D6 !!!!
    modelCase%nt = 2.D4  !!!!
    modelCase%visc = 1.D-3
    modelCase%rhof = 1.01D3
    modelCase%rhos = 2.71D3
    modelCase%alpha = 1.D0
    modelCase%epslon = 0.D0
    modelCase%dm = 3.6D-9
    modelCase%alphaOS = 5.D-1
    modelCase%lamdaX = 5.D-1
    modelCase%lamdaT = 1.D-1
    modelCase%radiusInit = 1.D-6
    modelCase%ShInfinity = 3.66D0
    modelCase%ks = 2.D-3
    modelCase%al = 5.D-2
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
    modelCase%poroInit(50,11:30) = 0.7

    allocate(modelCase%KxxInit(modelCase%nx,modelCase%ny))
    modelCase%KxxInit(:,:) = 1.D-15
    allocate(modelCase%KyyInit(modelCase%nx,modelCase%ny))
    modelCase%KyyInit(:,:) = 1.D-15
    ! add the fracture
    modelCase%KxxInit(50,11:30) = 1.D-12
    modelCase%KyyInit(50,11:30) = 1.D-12

    allocate(modelCase%avInit(modelCase%nx,modelCase%ny))
    modelCase%avInit(:,:) = 5.D-1!!!!!!

    allocate(modelCase%vxBdryX0(modelCase%ny))
    modelCase%vxBdryX0(:) = 4.D-6
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
    modelCase%pBdryX1(:) = 0.D0
    allocate(modelCase%pBdryY0(modelCase%nx))
    modelCase%pBdryY0(:) = 0.D0
    allocate(modelCase%pBdryY1(modelCase%nx))
    modelCase%pBdryY1(:) = 0.D0
    allocate(modelCase%pInit(modelCase%nx,modelCase%ny))
    modelCase%pInit(:,:) = 0.D0

    allocate(modelCase%isDiriX0_Cf(modelCase%ny))
    modelCase%isDiriX0_Cf(:) = 1
    allocate(modelCase%isDiriX1_Cf(modelCase%ny))
    modelCase%isDiriX1_Cf(:) = 0
    allocate(modelCase%isDiriY0_Cf(modelCase%nx))
    modelCase%isDiriY0_Cf(:) = 0
    allocate(modelCase%isDiriY1_Cf(modelCase%nx))
    modelCase%isDiriY1_Cf(:) = 0

    allocate(modelCase%CfBdryX0(modelCase%ny))
    modelCase%CfBdryX0(:) = 5.D2
    allocate(modelCase%CfBdryX1(modelCase%ny))
    modelCase%CfBdryX1(:) = 0.D0
    allocate(modelCase%CfBdryY0(modelCase%nx))
    modelCase%CfBdryY0(:) = 0.D0
    allocate(modelCase%CfBdryY1(modelCase%nx))
    modelCase%CfBdryY1(:) = 0.D0
    allocate(modelCase%CfInit(modelCase%nx,modelCase%ny))
    modelCase%CfInit(:,:) = 0.D0

    modelCase%soludoc = 'case1'

    call proceAlloc(1, modelCase%nx, modelCase%ny, modelCase%pncols, modelCase%pnrows)

    call driver(modelCase)

end program DBF_infile


Breakthrough has been achieved! Program stops now.
 The breakthrough time is    1301700.0000000000       seconds.
 The pore volume to breakthrough is    5.2067999999999994
 Solver time =    2615.4185700000280       seconds.
 Elapsed time =    3030.3381640000002       seconds.
