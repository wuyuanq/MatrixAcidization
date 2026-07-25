
!!$ Author:
!!$   Yuanqing Wu, DGUT, P.R.China
!!$
!!$ History:
!!$   2025-5-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

program DBF_infile

    use DBF_globalData
    use DBF_driver
    use RST_proceAlloc

    implicit none

    real(kind=8) :: rand
    real(kind=8), dimension(:), allocatable :: random
    integer :: i, j, k, times, ierr

    times = 1  !!!!
    Lx = 1.D-1
    Ly = 4.D-2
    nx = 180*times/2
    ny = 72*times/2
    timeEnd = 2.D6 !!!!
    nt = 2.D5  !!!!
    viscw = 1.D-3
    viscn = 1.D-2
    gammaw = 0.D0!5.8D2
    gamman = 0.D0!5.3D4
    gammawn = 0.D0!3.71D5
    rhofw = 1.01D3
    rhofn = 9.0D2
    rhos = 2.71D3
    dm = 3.6D-9
    dl = 1.D-3
    dt = 1.D-4
    alphaOS = 5.D-1
    radiusInit = 1.D-6
    ShInfinity = 3.66D0
    ks = 2.D-3
    al = 5.D-2
    gravX = 0.D0
    gravY = 0.D0!-9.807

    allocate(xs(nx+1))
    do i = 1, nx+1
        xs(i) = (i-1)*Lx/nx
    end do
    allocate(ys(ny+1))
    do j = 1, ny+1
        ys(j) = (j-1)*Ly/ny
    end do
    allocate(ts(nt+1))
    do k = 1, nt+1
        ts(k) = (k-1)*timeEnd/nt
    end do

    allocate(src(nx,ny))
    src(:,:) = 0.D0

    !call genRandomNum()
    allocate(random(RANDOMSIZE))
    open(unit=10, file=trim(adjustl(FRANDOMTXT)), iostat=ierr)
    if(ierr /= 0) then
        print *, 'open file error. ', ierr
        stop
    end if
    read(10, fmt="(f8.6)") random(:)
    close(10)
    allocate(poroInit(nx,ny))
    k = 0
    do j = 1, ny
        do i = 1, nx
            k = k + 1
            if(k <= RANDOMSIZE) then
                poroInit(i,j) = random(k)
            else
                print *, 'The random numbers are not enough!'
                print *, 'Please regenerate the random number file.'
                stop
            end if
        end do
    end do
    deallocate(random)
poroInit(:,ny/2) = 0.6
poroInit(:,ny/2+1) = 0.6

    allocate(isDiriX0_Sw(ny))
    isDiriX0_Sw(:) = 1
    allocate(isDiriX1_Sw(ny))
    isDiriX1_Sw(:) = 0
    allocate(isDiriY0_Sw(nx))
    isDiriY0_Sw(:) = 0
    allocate(isDiriY1_Sw(nx))
    isDiriY1_Sw(:) = 0

    allocate(SwBdryX0(ny))
    SwBdryX0(:) = 1.D0
    allocate(SwBdryX1(ny))
    SwBdryX1(:) = 0.D0
    allocate(SwBdryY0(nx))
    SwBdryY0(:) = 0.D0
    allocate(SwBdryY1(nx))
    SwBdryY1(:) = 0.D0
    allocate(SwInit(nx,ny))
    SwInit(:,:) = 1.D-2

    allocate(KxxInit(nx,ny))
    KxxInit(:,:) = 1.D-15
    allocate(KyyInit(nx,ny))
    KyyInit(:,:) = 1.D-15
KxxInit(:,ny/2) = 1.D-10
KyyInit(:,ny/2) = 1.D-10
KxxInit(:,ny/2+1) = 1.D-10
KyyInit(:,ny/2+1) = 1.D-10

    allocate(avInit(nx,ny))
    avInit(:,:) = 5.D-1

    allocate(vxwBdryX0(ny))
    vxwBdryX0(:) = 4.D-6
    allocate(vxwBdryX1(ny))
    vxwBdryX1(:) = 0.D0
    allocate(vywBdryX0(ny+1))
    vywBdryX0(:) = 0.D0
    allocate(vywBdryX1(ny+1))
    vywBdryX1(:) = 0.D0
    allocate(vxwBdryY0(nx+1))
    vxwBdryY0(:) = 0.D0
    allocate(vxwBdryY1(nx+1))
    vxwBdryY1(:) = 0.D0
    allocate(vywBdryY0(nx))
    vywBdryY0(:) = 0.D0
    allocate(vywBdryY1(nx))
    vywBdryY1(:) = 0.D0

    allocate(vxnBdryX0(ny))
    vxnBdryX0(:) = 0.D0
    allocate(vxnBdryX1(ny))
    vxnBdryX1(:) = 0.D0
    allocate(vynBdryX0(ny+1))
    vynBdryX0(:) = 0.D0
    allocate(vynBdryX1(ny+1))
    vynBdryX1(:) = 0.D0
    allocate(vxnBdryY0(nx+1))
    vxnBdryY0(:) = 0.D0
    allocate(vxnBdryY1(nx+1))
    vxnBdryY1(:) = 0.D0
    allocate(vynBdryY0(nx))
    vynBdryY0(:) = 0.D0
    allocate(vynBdryY1(nx))
    vynBdryY1(:) = 0.D0

    allocate(isDiriX0_p(ny))
    isDiriX0_p(:) = 0
    allocate(isDiriX1_p(ny))
    isDiriX1_p(:) = 1
    allocate(isDiriY0_p(nx))
    isDiriY0_p(:) = 0
    allocate(isDiriY1_p(nx))
    isDiriY1_p(:) = 0

    allocate(pBdryX0(ny))
    pBdryX0(:) = 0.D0
    allocate(pBdryX1(ny))
    pBdryX1(:) = 0.D0
    allocate(pBdryY0(nx))
    pBdryY0(:) = 0.D0
    allocate(pBdryY1(nx))
    pBdryY1(:) = 0.D0
    allocate(pInit(nx,ny))
    pInit(:,:) = 0.D0

    allocate(isDiriX0_Cf(ny))
    isDiriX0_Cf(:) = 1
    allocate(isDiriX1_Cf(ny))
    isDiriX1_Cf(:) = 0
    allocate(isDiriY0_Cf(nx))
    isDiriY0_Cf(:) = 0
    allocate(isDiriY1_Cf(nx))
    isDiriY1_Cf(:) = 0

    allocate(CfBdryX0(ny))
    CfBdryX0(:) = 5.D2
    allocate(CfBdryX1(ny))
    CfBdryX1(:) = 0.D0
    allocate(CfBdryY0(nx))
    CfBdryY0(:) = 0.D0
    allocate(CfBdryY1(nx))
    CfBdryY1(:) = 0.D0
    allocate(CfInit(nx,ny))
    CfInit(:,:) = 0.D0

    soludoc = 'case'

    call proceAlloc(1, nx, ny, pncols, pnrows)

    call driver()

end program DBF_infile

Breakthrough has been achieved! Program stops now.
The breakthrough time is    226830.00000000000       seconds.
The pore volume to breakthrough is   0.90731999999999990
Solver time =    4712.9814759996498       seconds.
Elapsed time =    5735.5821210000004       seconds.
