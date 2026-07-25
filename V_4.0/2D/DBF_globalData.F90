
!!$ Author:
!!$   Yuanqing Wu, DGUT, P.R.China
!!$
!!$ History:
!!$   2025-5-9 by Yuanqing Wu
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
    integer :: pncols
    integer :: pnrows
    real(kind=8) :: Lx
    real(kind=8) :: Ly
    integer :: nx
    integer :: ny
    real(kind=8) :: timeEnd
    integer :: nt
    real(kind=8), target :: viscw, viscn
    real(kind=8) :: gammaw, gamman, gammawn
    real(kind=8), target :: rhofw, rhofn
    real(kind=8) :: rhos
    real(kind=8) :: dm, dl, dt
    real(kind=8) :: alphaOS
    real(kind=8) :: lamdaX
    real(kind=8) :: lamdaT
    real(kind=8) :: radiusInit
    real(kind=8) :: ShInfinity
    real(kind=8) :: ks
    real(kind=8) :: al
    real(kind=8) :: gravX
    real(kind=8) :: gravY
    real(kind=8), dimension(:), allocatable :: xs
    real(kind=8), dimension(:), allocatable :: ys
    real(kind=8), dimension(:), allocatable :: ts
    real(kind=8), dimension(:,:), allocatable :: src
    real(kind=8), dimension(:,:), allocatable :: poroInit
    real(kind=8), dimension(:,:), allocatable :: KxxInit
    real(kind=8), dimension(:,:), allocatable :: KyyInit
    real(kind=8), dimension(:,:), allocatable :: avInit
    integer, dimension(:), allocatable :: isDiriX0_Sw
    integer, dimension(:), allocatable :: isDiriX1_Sw
    integer, dimension(:), allocatable :: isDiriY0_Sw
    integer, dimension(:), allocatable :: isDiriY1_Sw
    real(kind=8), dimension(:), allocatable, target :: SwBdryX0
    real(kind=8), dimension(:), allocatable, target :: SwBdryX1
    real(kind=8), dimension(:), allocatable, target :: SwBdryY0
    real(kind=8), dimension(:), allocatable, target :: SwBdryY1
    real(kind=8), dimension(:,:), allocatable :: SwInit
    real(kind=8), dimension(:), allocatable, target :: vxwBdryX0, vxnBdryX0
    real(kind=8), dimension(:), allocatable, target :: vxwBdryX1, vxnBdryX1
    real(kind=8), dimension(:), allocatable, target :: vxwBdryY0, vxnBdryY0
    real(kind=8), dimension(:), allocatable, target :: vxwBdryY1, vxnBdryY1
    real(kind=8), dimension(:), allocatable, target :: vywBdryX0, vynBdryX0
    real(kind=8), dimension(:), allocatable, target :: vywBdryX1, vynBdryX1
    real(kind=8), dimension(:), allocatable, target :: vywBdryY0, vynBdryY0
    real(kind=8), dimension(:), allocatable, target :: vywBdryY1, vynBdryY1
    integer, dimension(:), allocatable :: isDiriX0_p
    integer, dimension(:), allocatable :: isDiriX1_p
    integer, dimension(:), allocatable :: isDiriY0_p
    integer, dimension(:), allocatable :: isDiriY1_p
    real(kind=8), dimension(:), allocatable :: pBdryX0
    real(kind=8), dimension(:), allocatable :: pBdryX1
    real(kind=8), dimension(:), allocatable :: pBdryY0
    real(kind=8), dimension(:), allocatable :: pBdryY1
    real(kind=8), dimension(:,:), allocatable :: pInit
    integer, dimension(:), allocatable :: isDiriX0_Cf
    integer, dimension(:), allocatable :: isDiriX1_Cf
    integer, dimension(:), allocatable :: isDiriY0_Cf
    integer, dimension(:), allocatable :: isDiriY1_Cf
    real(kind=8), dimension(:), allocatable :: CfBdryX0
    real(kind=8), dimension(:), allocatable :: CfBdryX1
    real(kind=8), dimension(:), allocatable :: CfBdryY0
    real(kind=8), dimension(:), allocatable :: CfBdryY1
    real(kind=8), dimension(:,:), allocatable :: CfInit
    character(len = 10) :: soludoc

    ! the global variables
    real(kind=8), dimension(:,:), pointer :: kc
    real(kind=8), dimension(:), pointer :: hx, hy
    real(kind=8), dimension(:,:), pointer :: poro, poro_old
    real(kind=8), dimension(:,:), pointer :: poroEdgeX, poroEdgeY
    real(kind=8), dimension(:,:), pointer :: poroEdgeX_old, poroEdgeY_old
    real(kind=8), dimension(:,:), pointer :: poroEdgeXInit, poroEdgeYInit
    real(kind=8), dimension(:,:), pointer :: Sw, Sw_old
    real(kind=8), dimension(:,:), pointer :: SwEdgeX, SwEdgeY
    real(kind=8), dimension(:,:), pointer :: SwEdgeX_old, SwEdgeY_old
    real(kind=8), dimension(:,:), pointer :: Kxx, Kyy
    real(kind=8), dimension(:,:), pointer :: KxxEdge, KyyEdge
    real(kind=8), dimension(:,:), pointer :: av
    real(kind=8), dimension(:,:), pointer :: vxw, vyw, vxn, vyn
    real(kind=8), dimension(:,:), pointer :: p
    real(kind=8), dimension(:,:), pointer :: Cf

    real(kind=8), dimension(:), pointer :: local_rhs_vp, local_rhs_vp_static, local_rhs_Cf, local_rhs_Sw
    integer, dimension(:), pointer :: AxxwCols, AxpwCols, AyywCols, AypwCols, AcxwCols, AcywCols
    integer, dimension(:), pointer :: AxxnCols, AxpnCols, AyynCols, AypnCols, AcxnCols, AcynCols
    integer, dimension(:), pointer :: ACfCols, ASwCols
    integer, dimension(:), pointer :: AxxwRows, AxpwRows, AyywRows, AypwRows, AcxwRows, AcywRows
    integer, dimension(:), pointer :: AxxnRows, AxpnRows, AyynRows, AypnRows, AcxnRows, AcynRows
    integer, dimension(:), pointer :: ACfRows, ASwRows
    real(kind=8), dimension(:), pointer :: AxxwValues, AxpwValues, AyywValues, AypwValues, AcxwValues, AcywValues
    real(kind=8), dimension(:), pointer :: AxxwStaticValues, AxxwDynValues, AyywStaticValues, AyywDynValues
    real(kind=8), dimension(:), pointer :: AxxnValues, AxpnValues, AyynValues, AypnValues, AcxnValues, AcynValues
    real(kind=8), dimension(:), pointer :: AxxnStaticValues, AxxnDynValues, AyynStaticValues, AyynDynValues
    real(kind=8), dimension(:), pointer :: ACfValues, ASwValues
    integer :: AxxSize, AxpSize, AyySize, AypSize, AcxSize, AcySize, ACfSize, ASwSize
    integer, dimension(:), pointer :: AxxEntryNum, AxpEntryNum, AyyEntryNum, AypEntryNum, &!
        AcxEntryNum, AcyEntryNum, ACfEntryNum, ASwEntryNum
    integer, dimension(:), pointer :: AxxEntryBase, AxpEntryBase, AyyEntryBase, AypEntryBase, &!
        AcxEntryBase, AcyEntryBase, ACfEntryBase, ASwEntryBase

    integer :: nProcs, myid
    integer :: buffer_size = MAX_BUF_SIZE
    real(kind=8) :: buffer(MAX_BUF_SIZE)
    integer :: localncols, localnrows
    integer :: pcol, prow
    integer :: xlower, xupper, ylower, yupper
    integer :: ilower_vp, ilower_Cf, ilower_Sw
    integer :: iupper_vp, iupper_Cf, iupper_Sw
    integer :: local_x_size_vp, local_x_size_Cf, local_x_size_Sw
    integer, dimension(:), allocatable :: slave_vp_data_size
    real(kind=8) :: timestart, solvertime

    real(kind=8) :: presDropInit
    logical :: isFindPresDropInit
    integer :: t

#ifdef LAPACK

    ! the variables used in LAPACK
    real(kind=8), dimension(:,:), pointer :: A_lapack_vp, A_lapack_Cf, A_lapack_Sw
    real(kind=8), dimension(:), pointer :: b_lapack_vp, b_lapack_Cf, b_lapack_Sw
    integer :: LAPACKINFO
    integer, dimension(:), pointer :: IPIV_vp, IPIV_Cf, IPIV_Sw

#elif defined(UMFPACK)

    ! the variables used in UMFPACK
    integer, dimension(:), allocatable :: Ap_vp, Ap_Cf, Ap_Sw
    integer, dimension(:), allocatable :: Ai_vp, Ai_Cf, Ai_Sw
    real(kind=8), dimension(:), allocatable :: Ax_vp, Ax_Cf, Ax_Sw
    integer(kind=8) :: symbolic, numeric
    real(kind=8) :: control(20), umfinfo(90)

#elif defined(MUMPS)

    ! the variables used in MUMPS
    type(DMUMPS_STRUC) mumps_par_vp, mumps_par_Cf, mumps_par_Sw
    integer(kind=8) :: mumps_NNZ_loc_vp, mumps_NNZ_loc_Cf, mumps_NNZ_loc_Sw
    integer, dimension(:), allocatable :: mumps_IRN_loc_vp, mumps_IRN_loc_Cf, mumps_IRN_loc_Sw
    integer, dimension(:), allocatable :: mumps_JCN_loc_vp, mumps_JCN_loc_Cf, mumps_JCN_loc_Sw
    real(kind=8), dimension(:), allocatable :: mumps_A_loc_vp, mumps_A_loc_Cf, mumps_A_loc_Sw

#elif defined(HYPRE)

    ! the parameters and variables used in HYPRE
    integer, parameter :: HYPRE_PARCSR = 5555
    integer(kind=8) :: A_vp, A_Cf, A_Sw
    integer(kind=8) :: b_vp, b_Cf, b_Sw
    integer(kind=8) :: x_vp, x_Cf, x_Sw
    integer(kind=8) :: parcsr_A_vp, parcsr_A_Cf, parcsr_A_Sw
    integer(kind=8) :: par_b_vp, par_b_Cf, par_b_Sw
    integer(kind=8) :: par_x_vp, par_x_Cf, par_x_Sw
    integer(kind=8) :: precond_vp, precond_Cf, precond_Sw
    integer(kind=8) :: solver_vp, solver_Cf, solver_Sw
    integer :: jlower_vp, jlower_Cf, jlower_Sw
    integer :: jupper_vp, jupper_Cf, jupper_Sw
    integer, dimension(:), pointer :: rows_vp, rows_Cf, rows_Sw
    real(kind=8), dimension(:), allocatable :: initial_x_guess_vp, initial_x_guess_Cf, initial_x_guess_Sw

#endif

end module DBF_globalData

