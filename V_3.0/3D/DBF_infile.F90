
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
    integer :: i, j, k, c, ierr

    ! 0.5 M HCl, CaCO3
    modelCase%Lx = (4.D-2)/1
    modelCase%Ly = (4.D-2)/1
    modelCase%Lz = (1.D-1)/1
    modelCase%nx = 4!36/1
    modelCase%ny = 4!36/1
    modelCase%nz = 4!90/1
    modelCase%timeEnd = 3.44D6
    modelCase%nt = 5.8D3
    modelCase%visc = 1.D-3
    modelCase%rhof = 1.01D3
    modelCase%rhos = 2.71D3
    modelCase%alpha = 1.D0
    modelCase%dmRef = 3.6D-9
    modelCase%TemRef_dm = 2.98D2
    modelCase%ksRef = 2.D-3
    modelCase%TemRef_ks = 2.98D2
    modelCase%epslon = 0.D0
    modelCase%alphaOS = 5.D-1
    modelCase%lamdaX = 5.D-1
    modelCase%lamdaT = 1.D-1
    modelCase%lamdaf = 5.8D-1 !6.D-1
    modelCase%lamdas = 5.526D0 !1.4D0
    modelCase%thetaf = 4.184D3 !2.51D3
    modelCase%thetas = 2.D2 !8.1D2
    modelCase%radiusInit = 1.D-6
    modelCase%ShInfinity = 3.66D0
    modelCase%al = 5.D-2
    modelCase%gravX = 0.D0
    modelCase%gravY = 0.D0
    modelCase%gravZ = 0.D0

    allocate(modelCase%xs(modelCase%nx+1))
    do i = 1, modelCase%nx+1
        modelCase%xs(i) = (i-1)*modelCase%Lx/modelCase%nx
    end do
    allocate(modelCase%ys(modelCase%ny+1))
    do j = 1, modelCase%ny+1
        modelCase%ys(j) = (j-1)*modelCase%Ly/modelCase%ny
    end do
    allocate(modelCase%zs(modelCase%nz+1))
    do k = 1, modelCase%nz+1
        modelCase%zs(k) = (k-1)*modelCase%Lz/modelCase%nz
    end do
    allocate(modelCase%ts(modelCase%nt+1))
    do k = 1, modelCase%nt+1
        modelCase%ts(k) = (k-1)*modelCase%timeEnd/modelCase%nt
    end do

    allocate(modelCase%src(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%src(:,:,:) = 0.D0

    !call genRandomNum()
    allocate(random(RANDOMSIZE))
    open(unit=10, file=trim(adjustl(FRANDOMTXT)), iostat=ierr)
    if(ierr /= 0) then
        print *, 'open file error. ', ierr
        stop
    end if
    read(10, fmt="(f8.6)") random(:)
    close(10)
    allocate(modelCase%poroInit(modelCase%nx,modelCase%ny,modelCase%nz))
    c = 0
    do k = 1, modelCase%nz
        do j = 1, modelCase%ny
            do i = 1, modelCase%nx
                c = c + 1
                if(c <= RANDOMSIZE) then
                    modelCase%poroInit(i,j,k) = random(c)
                else
                    print *, 'The random numbers are not enough!'
                    print *, 'Please regenerate the random number file.'
                    stop
                end if
            end do
        end do
    end do
    deallocate(random)

    allocate(modelCase%KxxInit(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%KxxInit(:,:,:) = 9.869233D-16
    allocate(modelCase%KyyInit(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%KyyInit(:,:,:) = 9.869233D-16
    allocate(modelCase%KzzInit(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%KzzInit(:,:,:) = 9.869233D-16

    allocate(modelCase%avInit(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%avInit(:,:,:) = 5.D-1!5.D3

    allocate(modelCase%vxBdryX0(modelCase%ny,modelCase%nz))
    modelCase%vxBdryX0(:,:) = 0.D0
    allocate(modelCase%vxBdryX1(modelCase%ny,modelCase%nz))
    modelCase%vxBdryX1(:,:) = 0.D0
    allocate(modelCase%vyBdryX0(modelCase%ny+1,modelCase%nz))
    modelCase%vyBdryX0(:,:) = 0.D0
    allocate(modelCase%vyBdryX1(modelCase%ny+1,modelCase%nz))
    modelCase%vyBdryX1(:,:) = 0.D0
    allocate(modelCase%vzBdryX0(modelCase%ny,modelCase%nz+1))
    modelCase%vzBdryX0(:,:) = 0.D0
    allocate(modelCase%vzBdryX1(modelCase%ny,modelCase%nz+1))
    modelCase%vzBdryX1(:,:) = 0.D0

    allocate(modelCase%vxBdryY0(modelCase%nx+1,modelCase%nz))
    modelCase%vxBdryY0(:,:) = 0.D0
    allocate(modelCase%vxBdryY1(modelCase%nx+1,modelCase%nz))
    modelCase%vxBdryY1(:,:) = 0.D0
    allocate(modelCase%vyBdryY0(modelCase%nx,modelCase%nz))
    modelCase%vyBdryY0(:,:) = 0.D0
    allocate(modelCase%vyBdryY1(modelCase%nx,modelCase%nz))
    modelCase%vyBdryY1(:,:) = 0.D0
    allocate(modelCase%vzBdryY0(modelCase%nx,modelCase%nz+1))
    modelCase%vzBdryY0(:,:) = 0.D0
    allocate(modelCase%vzBdryY1(modelCase%nx,modelCase%nz+1))
    modelCase%vzBdryY1(:,:) = 0.D0

    allocate(modelCase%vxBdryZ0(modelCase%nx+1,modelCase%ny))
    modelCase%vxBdryZ0(:,:) = 0.D0
    allocate(modelCase%vxBdryZ1(modelCase%nx+1,modelCase%ny))
    modelCase%vxBdryZ1(:,:) = 0.D0
    allocate(modelCase%vyBdryZ0(modelCase%nx,modelCase%ny+1))
    modelCase%vyBdryZ0(:,:) = 0.D0
    allocate(modelCase%vyBdryZ1(modelCase%nx,modelCase%ny+1))
    modelCase%vyBdryZ1(:,:) = 0.D0
    allocate(modelCase%vzBdryZ0(modelCase%nx,modelCase%ny))
    modelCase%vzBdryZ0(:,:) = 0.D0
    allocate(modelCase%vzBdryZ1(modelCase%nx,modelCase%ny))
    modelCase%vzBdryZ1(:,:) = 1.04D-6

    allocate(modelCase%isDiriX0_p(modelCase%ny,modelCase%nz))
    modelCase%isDiriX0_p(:,:) = 0
    allocate(modelCase%isDiriX1_p(modelCase%ny,modelCase%nz))
    modelCase%isDiriX1_p(:,:) = 0
    allocate(modelCase%isDiriY0_p(modelCase%nx,modelCase%nz))
    modelCase%isDiriY0_p(:,:) = 0
    allocate(modelCase%isDiriY1_p(modelCase%nx,modelCase%nz))
    modelCase%isDiriY1_p(:,:) = 0
    allocate(modelCase%isDiriZ0_p(modelCase%nx,modelCase%ny))
    modelCase%isDiriZ0_p(:,:) = 1
    allocate(modelCase%isDiriZ1_p(modelCase%nx,modelCase%ny))
    modelCase%isDiriZ1_p(:,:) = 0

    allocate(modelCase%pBdryX0(modelCase%ny,modelCase%nz))
    modelCase%pBdryX0(:,:) = 0.D0
    allocate(modelCase%pBdryX1(modelCase%ny,modelCase%nz))
    modelCase%pBdryX1(:,:) = 0.D0
    allocate(modelCase%pBdryY0(modelCase%nx,modelCase%nz))
    modelCase%pBdryY0(:,:) = 0.D0
    allocate(modelCase%pBdryY1(modelCase%nx,modelCase%nz))
    modelCase%pBdryY1(:,:) = 0.D0
    allocate(modelCase%pBdryZ0(modelCase%nx,modelCase%ny))
    modelCase%pBdryZ0(:,:) = 1.52D7
    allocate(modelCase%pBdryZ1(modelCase%nx,modelCase%ny))
    modelCase%pBdryZ1(:,:) = 0.D0

    allocate(modelCase%pInit(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%pInit(:,:,:) = 1.52D7

    allocate(modelCase%isDiriX0_Cf(modelCase%ny,modelCase%nz))
    modelCase%isDiriX0_Cf(:,:) = 0
    allocate(modelCase%isDiriX1_Cf(modelCase%ny,modelCase%nz))
    modelCase%isDiriX1_Cf(:,:) = 0
    allocate(modelCase%isDiriY0_Cf(modelCase%nx,modelCase%nz))
    modelCase%isDiriY0_Cf(:,:) = 0
    allocate(modelCase%isDiriY1_Cf(modelCase%nx,modelCase%nz))
    modelCase%isDiriY1_Cf(:,:) = 0
    allocate(modelCase%isDiriZ0_Cf(modelCase%nx,modelCase%ny))
    modelCase%isDiriZ0_Cf(:,:) = 0
    allocate(modelCase%isDiriZ1_Cf(modelCase%nx,modelCase%ny))
    modelCase%isDiriZ1_Cf(:,:) = 1

    allocate(modelCase%CfBdryX0(modelCase%ny,modelCase%nz))
    modelCase%CfBdryX0(:,:) = 0.D0
    allocate(modelCase%CfBdryX1(modelCase%ny,modelCase%nz))
    modelCase%CfBdryX1(:,:) = 0.D0
    allocate(modelCase%CfBdryY0(modelCase%nx,modelCase%nz))
    modelCase%CfBdryY0(:,:) = 0.D0
    allocate(modelCase%CfBdryY1(modelCase%nx,modelCase%nz))
    modelCase%CfBdryY1(:,:) = 0.D0 ! mole/m^3
    allocate(modelCase%CfBdryZ0(modelCase%nx,modelCase%ny))
    modelCase%CfBdryZ0(:,:) = 0.D0
    allocate(modelCase%CfBdryZ1(modelCase%nx,modelCase%ny))
    modelCase%CfBdryZ1(:,:) = 5.D2

    allocate(modelCase%CfInit(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%CfInit(:,:,:) = 0.D0

    allocate(modelCase%isDiriX0_Tem(modelCase%ny,modelCase%nz))
    modelCase%isDiriX0_Tem(:,:) = 0
    allocate(modelCase%isDiriX1_Tem(modelCase%ny,modelCase%nz))
    modelCase%isDiriX1_Tem(:,:) = 0
    allocate(modelCase%isDiriY0_Tem(modelCase%nx,modelCase%nz))
    modelCase%isDiriY0_Tem(:,:) = 0
    allocate(modelCase%isDiriY1_Tem(modelCase%nx,modelCase%nz))
    modelCase%isDiriY1_Tem(:,:) = 0
    allocate(modelCase%isDiriZ0_Tem(modelCase%nx,modelCase%ny))
    modelCase%isDiriZ0_Tem(:,:) = 0
    allocate(modelCase%isDiriZ1_Tem(modelCase%nx,modelCase%ny))
    modelCase%isDiriZ1_Tem(:,:) = 1

    allocate(modelCase%TemBdryX0(modelCase%ny,modelCase%nz))
    modelCase%TemBdryX0(:,:) = 0.D0
    allocate(modelCase%TemBdryX1(modelCase%ny,modelCase%nz))
    modelCase%TemBdryX1(:,:) = 0.D0
    allocate(modelCase%TemBdryY0(modelCase%nx,modelCase%nz))
    modelCase%TemBdryY0(:,:) = 0.D0
    allocate(modelCase%TemBdryY1(modelCase%nx,modelCase%nz))
    modelCase%TemBdryY1(:,:) = 0.D0
    allocate(modelCase%TemBdryZ0(modelCase%nx,modelCase%ny))
    modelCase%TemBdryZ0(:,:) = 0.D0
    allocate(modelCase%TemBdryZ1(modelCase%nx,modelCase%ny))
    modelCase%TemBdryZ1(:,:) = 2.95D2

    allocate(modelCase%TemInit(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%TemInit(:,:,:) = 2.95D2

    modelCase%soludoc = 'case1'

    call proceAlloc(8, modelCase%nx, modelCase%ny, modelCase%nz, modelCase%pncols, &!
        modelCase%pnrows, modelCase%pnlays)  !!!!!!!

    call driver(modelCase)

end program DBF_infile
