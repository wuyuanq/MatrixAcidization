
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_exportResults
  
    use DBF_model
    use DBF_globalData
    use DBF_export2Matlab
    use DBF_export2tecplot
    implicit none

contains

    subroutine exportResults(global_data)

        real(kind=8), dimension(:), pointer, intent(in) :: global_data

        character(len=10) :: chart
        integer :: p_pcol, p_prow
        integer :: p_xlower, p_ylower
        integer :: indexl, indexr, indexd, indexu
        real(kind=8), dimension(:,:), pointer :: g_poro
        real(kind=8), dimension(:,:), pointer :: g_Kxx
        real(kind=8), dimension(:,:), pointer :: g_vx
        real(kind=8), dimension(:,:), pointer :: g_vy
        real(kind=8), dimension(:,:), pointer :: g_p
        real(kind=8), dimension(:,:), pointer :: g_Cf
        real(kind=8), dimension(:,:), pointer :: g_Tem
        character(len=40) :: fporotxt, fKxxtxt, fvxtxt, fvytxt, fptxt, fCftxt, fTemtxt, &!
            fvxmidtxt, fvymidtxt
        integer :: pid, ierr
        integer :: i, j, c

        write(chart,'(i10)') t
        fporotxt = trim(adjustl(soludoc))//'/soln_poro_raw_'//trim(adjustl(chart))//'.txt'
        fKxxtxt = trim(adjustl(soludoc))//'/soln_Kxx_raw_'//trim(adjustl(chart))//'.txt'
        fvxtxt = trim(adjustl(soludoc))//'/soln_vx_raw_'//trim(adjustl(chart))//'.txt'
        fvytxt = trim(adjustl(soludoc))//'/soln_vy_raw_'//trim(adjustl(chart))//'.txt'
        fptxt = trim(adjustl(soludoc))//'/soln_p_raw_'//trim(adjustl(chart))//'.txt'
        fCftxt = trim(adjustl(soludoc))//'/soln_Cf_raw_'//trim(adjustl(chart))//'.txt'
        fTemtxt = trim(adjustl(soludoc))//'/soln_Tem_raw_'//trim(adjustl(chart))//'.txt'
        fvxmidtxt = trim(adjustl(soludoc))//'/soln_vxmid_raw_'//trim(adjustl(chart))//'.txt'
        fvymidtxt = trim(adjustl(soludoc))//'/soln_vymid_raw_'//trim(adjustl(chart))//'.txt'

        allocate(g_poro(nx,ny))
        allocate(g_Kxx(nx,ny))
        allocate(g_vx(nx+1,ny))
        allocate(g_vy(nx,ny+1))
        allocate(g_p(nx,ny))
        allocate(g_Cf(nx,ny))
        allocate(g_Tem(nx,ny))

        c = 0
        do pid = 0, nProcs-1

            p_pcol = mod(pid,pncols)+1
            p_prow = pid/pncols+1

            p_xlower = (p_pcol-1)*localncols+1
            p_ylower = (p_prow-1)*localnrows+1

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_poro(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_Kxx(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            indexl = 1
            if(p_pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            indexd = 1
            indexu = localnrows

            do j = indexd, indexu
                do i = indexl, indexr
                    c = c + 1
                    g_vx(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            indexl = 1
            indexr = localncols
            indexd = 1
            if(p_prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if

            do j = indexd, indexu
                do i = indexl, indexr
                    c = c + 1
                    g_vy(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_p(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_Cf(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols            
                    c = c + 1
                    g_Tem(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

        end do

        ! export results to the raw txt files
        open(unit=10, file=fporotxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fporotxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_poro(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fporotxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fKxxtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fKxxtxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_Kxx(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fKxxtxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fvxtxt, status='replace', iostat=ierr)
        open(unit=11, file=fvxmidtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fvxtxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx+1
                write(10, fmt='(es24.16)', iostat=ierr) g_vx(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fvxtxt, ' error. ', ierr
                    stop
                end if
                if(i == nx/2+1) then
                    write(11, fmt='(es12.5)', iostat=ierr) (ys(j)+ys(j+1))/2.D0
                    write(11, fmt='(es12.5)', iostat=ierr) g_vx(i,j)
                end if
            end do
        end do
        close(10)
        close(11)

        open(unit=10, file=fvytxt, status='replace', iostat=ierr)
        open(unit=11, file=fvymidtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fvytxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny+1
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_vy(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fvytxt, ' error. ', ierr
                    stop
                end if
                if(j == ny/2+1) then
                    write(11, fmt='(es12.5)', iostat=ierr) (xs(i)+xs(i+1))/2.D0
                    write(11, fmt='(es12.5)', iostat=ierr) g_vy(i,j)
                end if
            end do
        end do
        close(10)
        close(11)

        open(unit=10, file=fptxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fptxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_p(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fptxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fCftxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fCftxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16e3)', iostat=ierr) g_Cf(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fCftxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fTemtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fTemtxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16e3)', iostat=ierr) g_Tem(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fTemtxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        ! export results to Matlab file
        call export2Matlab()

        ! export results to tecplot file
        call export2tecplot(g_poro, g_Kxx, g_vx, g_vy, g_p, g_Cf, g_Tem)

        deallocate(g_poro)
        deallocate(g_Kxx)
        deallocate(g_vx)
        deallocate(g_vy)
        deallocate(g_p)
        deallocate(g_Cf)
        deallocate(g_Tem)

    end subroutine exportResults

end module DBF_exportResults
