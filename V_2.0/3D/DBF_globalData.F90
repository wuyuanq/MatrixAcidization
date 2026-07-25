
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
    logical, parameter :: isForchheimer = .true.

    ! the model variables
    integer :: pncols, pnrows, pnlays
    real(kind=8) :: Lx, Ly, Lz
    integer :: nx, ny, nz
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
    real(kind=8) :: gravX, gravY, gravZ
    real(kind=8), dimension(:), allocatable :: xs, ys, zs
    real(kind=8), dimension(:), allocatable :: ts
    real(kind=8), dimension(:,:,:), allocatable :: src
    real(kind=8), dimension(:,:,:), allocatable :: poroInit
    real(kind=8), dimension(:,:,:), allocatable :: KxxInit, KyyInit, KzzInit
    real(kind=8), dimension(:,:,:), allocatable :: avInit
    real(kind=8), dimension(:,:), allocatable :: vxBdryX0, vxBdryX1, vxBdryY0, vxBdryY1, vxBdryZ0, &!
        vxBdryZ1, vyBdryX0, vyBdryX1, vyBdryY0, vyBdryY1, vyBdryZ0, vyBdryZ1, vzBdryX0, vzBdryX1, &!
        vzBdryY0, vzBdryY1, vzBdryZ0, vzBdryZ1
    integer, dimension(:,:), allocatable :: isDiriX0_p, isDiriX1_p, isDiriY0_p, isDiriY1_p, isDiriZ0_p, isDiriZ1_p
    real(kind=8), dimension(:,:), allocatable :: pBdryX0, pBdryX1, pBdryY0, pBdryY1, pBdryZ0, pBdryZ1
    real(kind=8), dimension(:,:,:), allocatable :: pInit
    integer, dimension(:,:), allocatable :: isDiriX0_Cf, isDiriX1_Cf, isDiriY0_Cf, isDiriY1_Cf, &!
        isDiriZ0_Cf, isDiriZ1_Cf
    real(kind=8), dimension(:,:), allocatable :: CfBdryX0, CfBdryX1, CfBdryY0, CfBdryY1, CfBdryZ0, CfBdryZ1
    real(kind=8), dimension(:,:,:), allocatable :: CfInit
    integer, dimension(:,:), allocatable :: isDiriX0_Tem, isDiriX1_Tem, isDiriY0_Tem, isDiriY1_Tem, &!
        isDiriZ0_Tem, isDiriZ1_Tem
    real(kind=8), dimension(:,:), allocatable :: TemBdryX0, TemBdryX1, TemBdryY0, TemBdryY1, TemBdryZ0, TemBdryZ1
    real(kind=8), dimension(:,:,:), allocatable :: TemInit
    character(len=20) :: soludoc

    ! the global parameters and variables
    real(kind=8), parameter :: Rg = 8.314
    real(kind=8), parameter :: Eg = 5.02416D4  !!!!!!!!!!!!!!!!
    real(kind=8), dimension(:,:,:), pointer :: dm
    real(kind=8), dimension(:,:,:), pointer :: kc
    real(kind=8), dimension(:,:,:), pointer :: ks
    real(kind=8), dimension(:), pointer :: hx, hy, hz
    real(kind=8), dimension(:,:,:), pointer :: poro, poro_old
    real(kind=8), dimension(:,:,:), pointer :: poroHarmX, poroHarmY, poroHarmZ
    real(kind=8), dimension(:,:,:), pointer :: poroHarmX_old, poroHarmY_old, poroHarmZ_old
    real(kind=8), dimension(:,:,:), pointer :: poroHarmXInit, poroHarmYInit, poroHarmZInit
    real(kind=8), dimension(:,:,:), allocatable :: Kxx, Kyy, Kzz
    real(kind=8), dimension(:,:,:), pointer :: KxxHarm, KyyHarm, KzzHarm
    real(kind=8), dimension(:,:,:), pointer :: av
    real(kind=8), dimension(:,:,:), pointer :: vx, vy, vz
    real(kind=8), dimension(:,:,:), pointer :: p
    real(kind=8), dimension(:,:,:), pointer :: Cf
    real(kind=8), dimension(:,:,:), pointer :: Tem
    
    real(kind=8), dimension(:), allocatable :: local_rhs, local_rhs_static, local_rhs_Cf, local_rhs_Tem
    integer, dimension(:), pointer :: AxxCols, AxpCols, AyyCols, AypCols, AzzCols, AzpCols, AcxCols, AcyCols, &!
        AczCols, AcpCols, AcfCols, AtemCols
    integer, dimension(:), pointer :: AxxRows, AxpRows, AyyRows, AypRows, AzzRows, AzpRows, AcxRows, AcyRows, &!
        AczRows, AcpRows, AcfRows, AtemRows
    real(kind=8), dimension(:), pointer :: AxxValues, AxpValues, AyyValues, AypValues, AzzValues, AzpValues, &!
        AcxValues, AcyValues, AczValues, AcpValues, AxxStaticValues, AxxDynValues, AyyStaticValues, AyyDynValues, &!
        AzzStaticValues, AzzDynValues, AcfValues, AtemValues
    integer :: AxxSize, AxpSize, AyySize, AypSize, AzzSize, AzpSize, AcxSize, AcySize, AczSize, AcpSize, AcfSize, AtemSize
    integer, dimension(:), pointer :: AxxEntryNum, AxpEntryNum, AyyEntryNum, AypEntryNum, AzzEntryNum, AzpEntryNum, &!
        AcxEntryNum, AcyEntryNum, AczEntryNum, AcpEntryNum, AcfEntryNum, AtemEntryNum
    integer, dimension(:), pointer :: AxxEntryBase, AxpEntryBase, AyyEntryBase, AypEntryBase, AzzEntryBase, &!
        AzpEntryBase, AcxEntryBase, AcyEntryBase, AczEntryBase, AcpEntryBase, AcfEntryBase, AtemEntryBase

    integer :: nProcs, myid
    integer :: buffer_size = MAX_BUF_SIZE
    real(kind=8) :: buffer(MAX_BUF_SIZE)
    integer :: localncols, localnrows, localnlays
    integer :: pcol, prow, play
    integer :: xlower, xupper, ylower, yupper, zlower, zupper
    integer :: ilower, iupper, ilower_Cf, iupper_Cf, ilower_Tem, iupper_Tem
    integer :: local_x_size, local_x_size_Cf, local_x_size_Tem
    integer, dimension(:), allocatable :: slave_vp_data_size
    real(kind=8) :: timestart, solvertime

    real(kind=8) :: presDropInit
    logical :: isFindPresDropInit
    integer :: t

#ifdef LAPACK

    ! the variables used in LAPACK
    real(kind=8), dimension(:,:), pointer :: A_lapack, A_lapack_Cf, A_lapack_Tem
    real(kind=8), dimension(:), pointer :: b_lapack, b_lapack_Cf, b_lapack_Tem
    integer :: LAPACKINFO
    integer, dimension(:), pointer :: IPIV, IPIV_Cf, IPIV_Tem

#elif defined(UMFPACK)

    ! the variables used in UMFPACK
    integer, dimension(:), allocatable :: Ap, Ai, Ap_Cf, Ai_Cf, Ap_Tem, Ai_Tem
    real(kind=8), dimension(:), allocatable :: Ax, Ax_Cf, Ax_Tem
    integer(kind=8) :: symbolic, numeric
    real(kind=8) :: control(20), umfinfo(90)

#elif defined(MUMPS)

    ! the variables used in MUMPS
    type(DMUMPS_STRUC) mumps_par, mumps_par_Cf, mumps_par_Tem
    integer(kind=8) :: mumps_NNZ_loc, mumps_NNZ_loc_Cf, mumps_NNZ_loc_Tem
    integer, dimension(:), allocatable ::  mumps_IRN_loc, mumps_JCN_loc, mumps_IRN_loc_Cf, mumps_JCN_loc_Cf, &!
        mumps_IRN_loc_Tem, mumps_JCN_loc_Tem
    real(kind=8), dimension(:), allocatable :: mumps_A_loc, mumps_A_loc_Cf, mumps_A_loc_Tem

#elif defined(HYPRE)

    ! the parameters and variables used in HYPRE
    integer, parameter :: HYPRE_PARCSR = 5555
    integer(kind=8) :: A, A_Cf, A_Tem
    integer(kind=8) :: b, b_Cf, b_Tem
    integer(kind=8) :: x, x_Cf, x_Tem
    integer(kind=8) :: parcsr_A, parcsr_A_Cf, parcsr_A_Tem
    integer(kind=8) :: par_b, par_b_Cf, par_b_Tem
    integer(kind=8) :: par_x, par_x_Cf, par_x_Tem
    integer(kind=8) :: precond, precond_Cf, precond_Tem
    integer(kind=8) :: solver, solver_Cf, solver_Tem
    integer :: jlower, jlower_Cf, jlower_Tem
    integer :: jupper, jupper_Cf, jupper_Tem
    integer, dimension(:), pointer :: rows, rows_Cf, rows_Tem
    real(kind=8), dimension(:), allocatable :: initial_x_guess, initial_x_guess_Cf, initial_x_guess_Tem

#endif

end module DBF_globalData

