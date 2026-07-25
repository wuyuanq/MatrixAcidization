
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_model

    implicit none

    type :: model

        integer :: pncols
        integer :: pnrows
        integer :: pnlays
        real(kind=8) :: Lx
        real(kind=8) :: Ly
        real(kind=8) :: Lz
        integer :: nx
        integer :: ny
        integer :: nz
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
        real(kind=8) :: gravX
        real(kind=8) :: gravY
        real(kind=8) :: gravZ
        real(kind=8), dimension(:), allocatable :: xs
        real(kind=8), dimension(:), allocatable :: ys
        real(kind=8), dimension(:), allocatable :: zs
        real(kind=8), dimension(:), allocatable :: ts
        real(kind=8), dimension(:,:,:), allocatable :: src
        real(kind=8), dimension(:,:,:), allocatable :: poroInit
        real(kind=8), dimension(:,:,:), allocatable :: KxxInit
        real(kind=8), dimension(:,:,:), allocatable :: KyyInit
        real(kind=8), dimension(:,:,:), allocatable :: KzzInit
        real(kind=8), dimension(:,:,:), allocatable :: avInit
        real(kind=8), dimension(:,:), allocatable :: vxBdryX0
        real(kind=8), dimension(:,:), allocatable :: vxBdryX1
        real(kind=8), dimension(:,:), allocatable :: vxBdryY0
        real(kind=8), dimension(:,:), allocatable :: vxBdryY1
        real(kind=8), dimension(:,:), allocatable :: vxBdryZ0
        real(kind=8), dimension(:,:), allocatable :: vxBdryZ1
        real(kind=8), dimension(:,:), allocatable :: vyBdryX0
        real(kind=8), dimension(:,:), allocatable :: vyBdryX1
        real(kind=8), dimension(:,:), allocatable :: vyBdryY0
        real(kind=8), dimension(:,:), allocatable :: vyBdryY1
        real(kind=8), dimension(:,:), allocatable :: vyBdryZ0
        real(kind=8), dimension(:,:), allocatable :: vyBdryZ1
        real(kind=8), dimension(:,:), allocatable :: vzBdryX0
        real(kind=8), dimension(:,:), allocatable :: vzBdryX1
        real(kind=8), dimension(:,:), allocatable :: vzBdryY0
        real(kind=8), dimension(:,:), allocatable :: vzBdryY1
        real(kind=8), dimension(:,:), allocatable :: vzBdryZ0
        real(kind=8), dimension(:,:), allocatable :: vzBdryZ1
        integer, dimension(:,:), allocatable :: isDiriX0_p
        integer, dimension(:,:), allocatable :: isDiriX1_p
        integer, dimension(:,:), allocatable :: isDiriY0_p
        integer, dimension(:,:), allocatable :: isDiriY1_p
        integer, dimension(:,:), allocatable :: isDiriZ0_p
        integer, dimension(:,:), allocatable :: isDiriZ1_p
        real(kind=8), dimension(:,:), allocatable :: pBdryX0
        real(kind=8), dimension(:,:), allocatable :: pBdryX1
        real(kind=8), dimension(:,:), allocatable :: pBdryY0
        real(kind=8), dimension(:,:), allocatable :: pBdryY1
        real(kind=8), dimension(:,:), allocatable :: pBdryZ0
        real(kind=8), dimension(:,:), allocatable :: pBdryZ1
        real(kind=8), dimension(:,:,:), allocatable :: pInit
        integer, dimension(:,:), allocatable :: isDiriX0_Cf
        integer, dimension(:,:), allocatable :: isDiriX1_Cf
        integer, dimension(:,:), allocatable :: isDiriY0_Cf
        integer, dimension(:,:), allocatable :: isDiriY1_Cf
        integer, dimension(:,:), allocatable :: isDiriZ0_Cf
        integer, dimension(:,:), allocatable :: isDiriZ1_Cf
        real(kind=8), dimension(:,:), allocatable :: CfBdryX0
        real(kind=8), dimension(:,:), allocatable :: CfBdryX1
        real(kind=8), dimension(:,:), allocatable :: CfBdryY0
        real(kind=8), dimension(:,:), allocatable :: CfBdryY1
        real(kind=8), dimension(:,:), allocatable :: CfBdryZ0
        real(kind=8), dimension(:,:), allocatable :: CfBdryZ1
        real(kind=8), dimension(:,:,:), allocatable :: CfInit
        integer, dimension(:,:), allocatable :: isDiriX0_Tem
        integer, dimension(:,:), allocatable :: isDiriX1_Tem        
        integer, dimension(:,:), allocatable :: isDiriY0_Tem
        integer, dimension(:,:), allocatable :: isDiriY1_Tem
        integer, dimension(:,:), allocatable :: isDiriZ0_Tem
        integer, dimension(:,:), allocatable :: isDiriZ1_Tem
        real(kind=8), dimension(:,:), allocatable :: TemBdryX0
        real(kind=8), dimension(:,:), allocatable :: TemBdryX1
        real(kind=8), dimension(:,:), allocatable :: TemBdryY0
        real(kind=8), dimension(:,:), allocatable :: TemBdryY1
        real(kind=8), dimension(:,:), allocatable :: TemBdryZ0
        real(kind=8), dimension(:,:), allocatable :: TemBdryZ1
        real(kind=8), dimension(:,:,:), allocatable :: TemInit
        character(len = 10) :: soludoc 

    end type model

end module DBF_model
