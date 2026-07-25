
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_export2tecplot

    use DBF_model
    use DBF_globalData
    implicit none

contains

    subroutine exportPointCenteredIJK(data, varName, zoneName, title, fileName)

        real(kind=8), dimension(:,:,:), pointer, intent(in) :: data
        character(len=*), intent(in) :: varName, zoneName, title, fileName

        integer :: ix, iy, iz, ierr

        open(unit=11, file=fileName, status="replace", action="write", iostat=ierr)
        if(ierr /= 0) then
            print *, "Failed to open '", fileName, "'"
            stop
        end if
        write(11, fmt=*) 'TITLE = "', trim(adjustl(title)), '"'
        write(11, fmt=*) 'VARIABLES = "X", "Y", "Z", "', trim(adjustl(varName)), '"'
        write(11,*) 'ZONE T="', trim(adjustl(zoneName)), '" DATAPACKING=POINT, I=', nx+1, ', J=', ny+1, ', K=', nz+1
        do iz = 1, nz+1
            do iy = 1, ny+1
                do ix = 1, nx+1
                    write(11,*) xs(ix), ys(iy), zs(iz), data(ix, iy, iz)
                end do
            end do
        end do
        close(11)

    end subroutine exportPointCenteredIJK

    subroutine exportPointCenteredIJK_vars(dats, varNames, zoneName, title, fileName)

        real(kind=8), dimension(:,:,:,:), pointer, intent(in) :: dats
        character(len=*), dimension(:), intent(in) :: varNames
        character(len=*), intent(in) :: zoneName, title, fileName

        integer :: ix, iy, iz, iVar, nVar, ierr

        nVar = size(varNames)
        open(unit=11, file=fileName, status="replace", action="write", iostat=ierr)
        if(ierr /= 0) then
            print *, "Failed to open '", fileName, "'" 
            stop
        end if
        write(11,fmt='(a,a,a)') 'TITLE = "', trim(adjustl(title)), '"'
        write(11,fmt='(a)', advance='no') 'VARIABLES = "X", "Y", "Z"'
        do iVar = 1, nVar 
            write(11,fmt='(a,a,a)', advance='no') ', "', trim(adjustl(varNames(iVar))), '"'
        end do 
        write(11,*)
        nVar = size(dats, 4)
        write(11,*) 'ZONE T="', trim(adjustl(zoneName)), '" DATAPACKING=POINT, I=', nx+1, ', J=', ny+1, ', K=', nz+1
        do iz = 1, nz+1
            do iy = 1, ny+1
                do ix = 1, nx+1
                    write(11, fmt='(g15.5,g15.5,g15.5)', advance='no') xs(ix), ys(iy), zs(iz)
                    do iVar = 1, nVar
                        write(11, fmt='(g15.5)', advance='no') dats(ix, iy, iz, iVar)
                    end do
                    write(11,*)
                end do
            end do
        end do
        close(11)

    end subroutine exportPointCenteredIJK_vars

    subroutine xFaceCtr2nodCtr(xFaceCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:,:), pointer, intent(in) :: xFaceCtrDat
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: nodCtrDat

        integer :: iy, iz

        nodCtrDat(:,:,:) = 0.D0
        do iz = 0, 1
            do iy = 0, 1
                nodCtrDat(:,1+iy:ny+iy,1+iz:nz+iz) = nodCtrDat(:,1+iy:ny+iy,1+iz:nz+iz) + xFaceCtrDat
            end do
        end do
        nodCtrDat(:,2:ny,:) = nodCtrDat(:,2:ny,:) / 2.0
        nodCtrDat(:,:,2:nz) = nodCtrDat(:,:,2:nz) / 2.0

    end subroutine xFaceCtr2nodCtr

    subroutine yFaceCtr2nodCtr(yFaceCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:,:), pointer, intent(in) :: yFaceCtrDat
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: nodCtrDat

        integer :: ix, iz

        nodCtrDat(:,:,:) = 0.D0
        do iz = 0, 1
            do ix = 0, 1
                nodCtrDat(1+ix:nx+ix,:,1+iz:nz+iz) = nodCtrDat(1+ix:nx+ix,:,1+iz:nz+iz) + yFaceCtrDat
            end do
        end do
        nodCtrDat(2:nx,:,:) = nodCtrDat(2:nx,:,:) / 2.D0
        nodCtrDat(:,:,2:nz) = nodCtrDat(:,:,2:nz) / 2.D0

    end subroutine yFaceCtr2nodCtr

    subroutine zFaceCtr2nodCtr(zFaceCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:,:), pointer, intent(in) :: zFaceCtrDat
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: nodCtrDat

        integer :: ix, iy

        nodCtrDat(:,:,:) = 0.D0
        do iy = 0, 1
            do ix = 0, 1
                nodCtrDat(1+ix:nx+ix,1+iy:ny+iy,:) = nodCtrDat(1+ix:nx+ix,1+iy:ny+iy,:) + zFaceCtrDat
            end do
        end do
        nodCtrDat(2:nx,:,:) = nodCtrDat(2:nx,:,:) / 2.D0
        nodCtrDat(:,2:ny,:) = nodCtrDat(:,2:ny,:) / 2.D0

    end subroutine zFaceCtr2nodCtr

    subroutine cellCtr2nodCtr(cellCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:,:), pointer, intent(in) :: cellCtrDat
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: nodCtrDat

        integer :: ix, iy, iz

        nodCtrDat(:,:,:) = 0.D0
        do iz = 0, 1
            do iy = 0, 1
                do ix = 0, 1
                    nodCtrDat(1+ix:nx+ix,1+iy:ny+iy,1+iz:nz+iz) = nodCtrDat(1+ix:nx+ix,1+iy:ny+iy,1+iz:nz+iz) + cellCtrDat
                end do
            end do
        end do
        nodCtrDat(2:nx,:,:) = nodCtrDat(2:nx,:,:) / 2.D0
        nodCtrDat(:,2:ny,:) = nodCtrDat(:,2:ny,:) / 2.D0
        nodCtrDat(:,:,2:nz) = nodCtrDat(:,:,2:nz) / 2.D0

    end subroutine cellCtr2nodCtr

    subroutine export2tecplot(g_poro, g_Kxx, g_vx, g_vy, g_vz, g_p, g_Cf)

        real(kind=8), dimension(:,:,:), pointer, intent(in) :: g_poro
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: g_Kxx
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: g_vx, g_vy, g_vz
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: g_p
        real(kind=8), dimension(:,:,:), pointer, intent(in) :: g_Cf

        character(len=10) :: chart
        real(kind=8), dimension(:,:,:), pointer :: poro_a
        real(kind=8), dimension(:,:,:), pointer :: Kxx_a
        real(kind=8), dimension(:,:,:), pointer :: vx_a, vy_a, vz_a
        real(kind=8), dimension(:,:,:), pointer :: p_a
        real(kind=8), dimension(:,:,:), pointer :: Cf_a
        real(kind=8), dimension(:,:,:,:), pointer :: vars_a

        allocate(poro_a(nx+1,ny+1,nz+1))
        allocate(Kxx_a(nx+1,ny+1,nz+1))
        allocate(vx_a(nx+1,ny+1,nz+1))
        allocate(vy_a(nx+1,ny+1,nz+1))
        allocate(vz_a(nx+1,ny+1,nz+1))
        allocate(p_a(nx+1,ny+1,nz+1))
        allocate(Cf_a(nx+1,ny+1,nz+1))
        allocate(vars_a(nx+1,ny+1,nz+1,7))

        call cellCtr2nodCtr(g_poro, poro_a)
        call cellCtr2nodCtr(g_Kxx, Kxx_a)
        call xFaceCtr2nodCtr(g_vx, vx_a)
        call yFaceCtr2nodCtr(g_vy, vy_a)
        call zFaceCtr2nodCtr(g_vz, vz_a)
        call cellCtr2nodCtr(g_p, p_a)
        call cellCtr2nodCtr(g_Cf, Cf_a)

        vars_a(:,:,:,1) = poro_a
        vars_a(:,:,:,2) = Kxx_a
        vars_a(:,:,:,3) = vx_a
        vars_a(:,:,:,4) = vy_a
        vars_a(:,:,:,5) = vz_a
        vars_a(:,:,:,6) = p_a
        vars_a(:,:,:,7) = Cf_a

        write(chart,'(i10)') t

        call exportPointCenteredIJK_vars(vars_a, (/"poro", "Kxx ", "vx  ", "vy  ", "vz  ", "p   ", "Cf  "/), &!
            "zone1", "result from single phase flow", trim(adjustl(soludoc))//"/out_singlePhase_"// &!
            trim(adjustl(chart))//".plt")

        deallocate(poro_a)
        deallocate(Kxx_a)
        deallocate(vx_a)
        deallocate(vy_a)
        deallocate(vz_a)
        deallocate(p_a)
        deallocate(Cf_a)
        deallocate(vars_a)

    end subroutine export2tecplot

end module DBF_export2tecplot

