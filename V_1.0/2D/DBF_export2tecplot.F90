
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

    subroutine exportPointCenteredIJ(data, varName, zoneName, title, fileName)

        real(kind=8), dimension(:,:), pointer, intent(in) :: data
        character(len=*), intent(in) :: varName, zoneName, title, fileName

        integer :: ix, iy, ierr

        open(unit=10, file=fileName, status="replace", action="write", iostat=ierr)
        if(ierr /= 0) then
            print *, "Failed to open '", fileName, "'"
            stop
        end if
        write(10, fmt=*) 'TITLE = "', trim(adjustl(title)), '"'
        write(10, fmt=*) 'VARIABLES = "X", "Y", "', trim(adjustl(varName)), '"'
        write(10,*) 'ZONE T="', trim(adjustl(zoneName)), '" DATAPACKING=POINT, I=', nx+1, ', J=', ny+1
        do iy = 1, ny+1
            do ix = 1, nx+1
                write(10,*) xs(ix), ys(iy), data(ix, iy)
            end do
        end do
        close(10)

    end subroutine exportPointCenteredIJ

    subroutine exportPointCenteredIJ_vars(dats, varNames, zoneName, title, fileName)

        real(kind=8), dimension(:,:,:), pointer, intent(in) :: dats
        character(len=*), dimension(:), intent(in) :: varNames
        character(len=*), intent(in) :: zoneName, title, fileName

        integer :: ix, iy, iVar, nVar, ierr

        nVar = size(varNames)
        open(unit=10, file=fileName, status="replace", action="write", iostat=ierr)
        if(ierr /= 0) then
            print *, "Failed to open '", fileName, "'" 
            stop
        end if
        write(10,fmt='(a,a,a)') 'TITLE = "', trim(adjustl(title)), '"'
        write(10,fmt='(a)', advance='no') 'VARIABLES = "X", "Y"'
        do iVar = 1, nVar 
            write(10,fmt='(a,a,a)', advance='no') ', "', trim(adjustl(varNames(iVar))), '"'
        end do 
        write(10,*)
        nVar = size(dats, 3)  
        write(10,*) 'ZONE T="', trim(adjustl(zoneName)), '" DATAPACKING=POINT, I=', nx+1, ', J=', ny+1
        do iy = 1, ny+1
            do ix = 1, nx+1
                write(10, fmt='(g15.5,g15.5)', advance='no') xs(ix), ys(iy)
                do iVar = 1, nVar
                    write(10, fmt='(g15.5)', advance='no') dats(ix, iy, iVar)
                end do
                write(10,*)
            end do
        end do
        close(10)

    end subroutine exportPointCenteredIJ_vars

    subroutine xEdgCtr2nodCtr(xEdgCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:), pointer, intent(in) :: xEdgCtrDat
        real(kind=8), dimension(:,:), pointer, intent(in out) :: nodCtrDat

        integer :: iy

        nodCtrDat(:,:) = 0
        do iy = 0, 1
            nodCtrDat(:,1+iy:ny+iy) = nodCtrDat(:,1+iy:ny+iy) + xEdgCtrDat
        end do
        nodCtrDat(:,2:ny) = nodCtrDat(:,2:ny) / 2.D0

    end subroutine xEdgCtr2nodCtr

    subroutine yEdgCtr2nodCtr(yEdgCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:), pointer, intent(in) :: yEdgCtrDat
        real(kind=8), dimension(:,:), pointer, intent(in out) :: nodCtrDat

        integer :: ix

        nodCtrDat(:,:) = 0
        do ix = 0, 1
            nodCtrDat(1+ix:nx+ix,:) = nodCtrDat(1+ix:nx+ix,:) + yEdgCtrDat
        end do
        nodCtrDat(2:nx,:) = nodCtrDat(2:nx,:) / 2.D0

    end subroutine yEdgCtr2nodCtr

    subroutine cellCtr2nodCtr(cellCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:), pointer, intent(in) :: cellCtrDat
        real(kind=8), dimension(:,:), pointer, intent(in out) :: nodCtrDat

        integer :: ix, iy

        nodCtrDat(:,:) = 0
        do ix = 0, 1
            do iy = 0, 1
                nodCtrDat(1+ix:nx+ix,1+iy:ny+iy) = nodCtrDat(1+ix:nx+ix,1+iy:ny+iy) + cellCtrDat
            end do
        end do
        nodCtrDat(2:nx,:) = nodCtrDat(2:nx,:) / 2.D0
        nodCtrDat(:,2:ny) = nodCtrDat(:,2:ny) / 2.D0

    end subroutine cellCtr2nodCtr

    subroutine export2tecplot(g_poro, g_Kxx, g_vx, g_vy, g_p, g_Cf)

        real(kind=8), dimension(:,:), pointer, intent(in) :: g_poro
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_Kxx
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_vx, g_vy
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_p
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_Cf

        character(len=10) :: chart
        real(kind=8), dimension(:,:), pointer :: poro_a
        real(kind=8), dimension(:,:), pointer :: Kxx_a
        real(kind=8), dimension(:,:), pointer :: vx_a, vy_a
        real(kind=8), dimension(:,:), pointer :: p_a
        real(kind=8), dimension(:,:), pointer :: Cf_a
        real(kind=8), dimension(:,:,:), pointer :: vars_a

        allocate(poro_a(nx+1,ny+1))
        allocate(Kxx_a(nx+1,ny+1))
        allocate(vx_a(nx+1,ny+1))
        allocate(vy_a(nx+1,ny+1))
        allocate(p_a(nx+1,ny+1))
        allocate(Cf_a(nx+1,ny+1))
        allocate(vars_a(nx+1,ny+1,6))

        call cellCtr2nodCtr(g_poro, poro_a)
        call cellCtr2nodCtr(g_Kxx, Kxx_a)
        call xEdgCtr2nodCtr(g_vx, vx_a)
        call yEdgCtr2nodCtr(g_vy, vy_a)
        call cellCtr2nodCtr(g_p, p_a)
        call cellCtr2nodCtr(g_Cf, Cf_a)

        vars_a(:,:,1) = poro_a
        vars_a(:,:,2) = Kxx_a
        vars_a(:,:,3) = vx_a
        vars_a(:,:,4) = vy_a
        vars_a(:,:,5) = p_a
        vars_a(:,:,6) = Cf_a

        write(chart,'(i10)') t

        call exportPointCenteredIJ_vars(vars_a, (/"poro", "Kxx ", "vx  ", "vy  ", "p   ", "Cf  "/), &!
            "zone1", "result from single phase flow", trim(adjustl(soludoc))//"/out_singlePhase_"// &!
            trim(adjustl(chart))//".plt")

        deallocate(poro_a)
        deallocate(Kxx_a)
        deallocate(vx_a)
        deallocate(vy_a)
        deallocate(p_a)
        deallocate(Cf_a)
        deallocate(vars_a)

    end subroutine export2tecplot

end module DBF_export2tecplot

