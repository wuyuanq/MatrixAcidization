
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_globalData

    implicit none
#ifdef MUMPS
    include 'dmumps_struc.h'
#endif

    ! the program parameters
    integer, parameter :: MAX_BUF_SIZE = 1.D8
    integer, parameter :: NUMFRAME = 100 ! number of frames in the movie
    integer, parameter :: RANDOMSIZE = 1.D7
    character(len=20), parameter :: FRANDOMTXT = '../../random.txt'
    logical, parameter :: isDarcy = .true.
    logical, parameter :: isBrinkman = .true.
    logical, parameter :: isForchheimer = .true. ! it should be always true in this code

    ! the model variables
    integer :: pncols, pnrows
    real(kind=8) :: Lx, Ly
    integer :: nx, ny
    real(kind=8) :: timeEnd
    integer :: nt
    real(kind=8) :: visc
    real(kind=8) :: rhof
    real(kind=8) :: rhos
    real(kind=8) :: alpha
    real(kind=8) :: dmRef
    real(kind=8) :: TemRef_dm
    real(kind=8) :: ksRef
    real(kind=8) :: TemRef_ks
    real(kind=8) :: epslon
    real(kind=8) :: alphaOS
    real(kind=8) :: lamdaX
    real(kind=8) :: lamdaT
    real(kind=8) :: lamdaf
    real(kind=8) :: lamdas
    real(kind=8) :: thetaf
    real(kind=8) :: thetas
    real(kind=8) :: radiusInit
    real(kind=8) :: ShInfinity
    real(kind=8) :: al
    real(kind=8) :: gravX, gravY
    real(kind=8), dimension(:), allocatable :: xs, ys
    real(kind=8), dimension(:), allocatable :: ts
    real(kind=8), dimension(:,:), allocatable :: src
    real(kind=8), dimension(:,:), allocatable :: poroInit
    real(kind=8), dimension(:,:), allocatable :: KxxInit, KyyInit
    real(kind=8), dimension(:,:), allocatable :: avInit
    real(kind=8), dimension(:), allocatable :: vxBdryX0, vxBdryX1, vxBdryY0, vxBdryY1, &!
        vyBdryX0, vyBdryX1, vyBdryY0, vyBdryY1
    integer, dimension(:), allocatable :: isDiriX0_p, isDiriX1_p, isDiriY0_p, isDiriY1_p
    real(kind=8), dimension(:), allocatable :: pBdryX0, pBdryX1, pBdryY0, pBdryY1
    real(kind=8), dimension(:,:), allocatable :: pInit
    integer, dimension(:), allocatable :: isDiriX0_Cf, isDiriX1_Cf, isDiriY0_Cf, isDiriY1_Cf
    real(kind=8), dimension(:), allocatable :: CfBdryX0, CfBdryX1, CfBdryY0, CfBdryY1
    real(kind=8), dimension(:,:), allocatable :: CfInit
    integer, dimension(:), allocatable :: isDiriX0_Tem, isDiriX1_Tem, isDiriY0_Tem, isDiriY1_Tem
    real(kind=8), dimension(:), allocatable :: TemBdryX0, TemBdryX1, TemBdryY0, TemBdryY1
    real(kind=8), dimension(:,:), allocatable :: TemInit
    character(len=20) :: soludoc

    ! the global parameters and variables
    real(kind=8), parameter :: Rg = 8.314
    real(kind=8), parameter :: Eg = 5.02416D4  !!!!!!!!!!!!!!!!
    real(kind=8), dimension(:,:), pointer :: dm
    real(kind=8), dimension(:,:), pointer :: kc
    real(kind=8), dimension(:,:), pointer :: ks
    real(kind=8), dimension(:), pointer :: hx, hy
    real(kind=8), dimension(:,:), pointer :: poro, poro_old
    real(kind=8), dimension(:,:), pointer :: poroHarmX, poroHarmY
    real(kind=8), dimension(:,:), pointer :: poroHarmX_old, poroHarmY_old
    real(kind=8), dimension(:,:), pointer :: poroHarmXInit, poroHarmYInit
    real(kind=8), dimension(:,:), allocatable :: Kxx, Kyy
    real(kind=8), dimension(:,:), pointer :: KxxHarm, KyyHarm
    real(kind=8), dimension(:,:), pointer :: av
    real(kind=8), dimension(:,:), pointer :: vx, vy
    real(kind=8), dimension(:,:), pointer :: p
    real(kind=8), dimension(:,:), pointer :: Cf
    real(kind=8), dimension(:,:), pointer :: Tem

    real(kind=8), dimension(:), pointer :: local_rhs_static_v, local_rhs_v, local_rhs_p, &!
        local_rhs_Cf, local_rhs_Tem
    integer, dimension(:), pointer :: AxxCols, AyyCols, ApCols, AcfCols, AtemCols
    integer, dimension(:), pointer :: AxxRows, AyyRows, ApRows, AcfRows, AtemRows
    real(kind=8), dimension(:), pointer :: AxxValues, AyyValues, AxxStaticValues, AxxDynValues, &!
        AyyStaticValues, AyyDynValues, ApValues, AcfValues, AtemValues
    integer :: AxxSize, AyySize, ApSize, AcfSize, AtemSize
    integer :: AxxEntryNum_standard, AyyEntryNum_standard, ApEntryNum_standard, AcfEntryNum_standard, AtemEntryNum_standard
    integer, dimension(:), pointer :: AxxEntryNum, AyyEntryNum, ApEntryNum, AcfEntryNum, AtemEntryNum
    integer, dimension(:), pointer :: AxxEntryBase, AyyEntryBase, ApEntryBase, AcfEntryBase, AtemEntryBase

    integer :: nProcs, myid
    integer :: buffer_size = MAX_BUF_SIZE
    real(kind=8) :: buffer(MAX_BUF_SIZE)
    integer :: localncols, localnrows
    integer :: pcol, prow
    integer :: xlower, xupper, ylower, yupper
    integer :: ilower_v, iupper_v, ilower_p, iupper_p, ilower_Cf, iupper_Cf, ilower_Tem, iupper_Tem
    integer :: local_x_size_v, local_x_size_p, local_x_size_Cf, local_x_size_Tem
    integer, dimension(:), allocatable :: slave_v_data_size
    real(kind=8) :: timestart, solvertime

    real(kind=8) :: presDropInit
    logical :: isFindPresDropInit
    integer :: t

#ifdef HYPRE

    real(kind=8), dimension(:), allocatable :: initial_x_guess_v, initial_x_guess_p, &!
        initial_x_guess_Cf, initial_x_guess_Tem

#endif

end module DBF_globalData

