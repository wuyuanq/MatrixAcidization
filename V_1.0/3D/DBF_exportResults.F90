
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
        integer :: p_pcol, p_prow, p_play
        integer :: p_xlower, p_ylower, p_zlower
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        real(kind=8), dimension(:,:,:), pointer :: g_poro
        real(kind=8), dimension(:,:,:), pointer :: g_Kxx
        real(kind=8), dimension(:,:,:), pointer :: g_vx, g_vy, g_vz
        real(kind=8), dimension(:,:,:), pointer :: g_p
        real(kind=8), dimension(:,:,:), pointer :: g_Cf
        character(len=40) :: fporotxt, fKxxtxt, fvxtxt, fvytxt, fvztxt, fptxt, fCftxt, fvxmidtxt, fvymidtxt
        integer :: pid, ierr
        integer :: i, j, k, c

        write(chart,'(i10)') t
        fporotxt = trim(adjustl(soludoc))//'/soln_poro_raw_'//trim(adjustl(chart))//'.txt'
        fKxxtxt = trim(adjustl(soludoc))//'/soln_Kxx_raw_'//trim(adjustl(chart))//'.txt'
        fvxtxt = trim(adjustl(soludoc))//'/soln_vx_raw_'//trim(adjustl(chart))//'.txt'
        fvytxt = trim(adjustl(soludoc))//'/soln_vy_raw_'//trim(adjustl(chart))//'.txt'
        fvztxt = trim(adjustl(soludoc))//'/soln_vz_raw_'//trim(adjustl(chart))//'.txt'
        fptxt = trim(adjustl(soludoc))//'/soln_p_raw_'//trim(adjustl(chart))//'.txt'
        fCftxt = trim(adjustl(soludoc))//'/soln_Cf_raw_'//trim(adjustl(chart))//'.txt'
        fvxmidtxt = trim(adjustl(soludoc))//'/soln_vxmid_raw_'//trim(adjustl(chart))//'.txt'
        fvymidtxt = trim(adjustl(soludoc))//'/soln_vymid_raw_'//trim(adjustl(chart))//'.txt'

        allocate(g_poro(nx,ny,nz))
        allocate(g_Kxx(nx,ny,nz))
        allocate(g_vx(nx+1,ny,nz))
        allocate(g_vy(nx,ny+1,nz))
        allocate(g_vz(nx,ny,nz+1))
        allocate(g_p(nx,ny,nz))
        allocate(g_Cf(nx,ny,nz))

        c = 0
        do pid = 0, nProcs-1

            p_play = pid/(pnrows*pncols)+1
            p_prow = (pid-(p_play-1)*pnrows*pncols)/pncols+1
            p_pcol = (pid-(p_play-1)*pnrows*pncols)-(p_prow-1)*pncols+1

            p_xlower = (p_pcol-1)*localncols+1
            p_ylower = (p_prow-1)*localnrows+1
            p_zlower = (p_play-1)*localnlays+1

            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        g_poro(p_xlower+i-1,p_ylower+j-1,p_zlower+k-1) = global_data(c)
                    end do
                end do
            end do

            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        g_Kxx(p_xlower+i-1,p_ylower+j-1,p_zlower+k-1) = global_data(c)
                    end do
                end do
            end do

            if(p_pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, indexr
                        c = c + 1
                        g_vx(p_xlower+i-1,p_ylower+j-1,p_zlower+k-1) = global_data(c)
                    end do
                end do
            end do

            if(p_prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            do k = 1, localnlays
                do j = 1, indexu
                    do i = 1, localncols
                        c = c + 1
                        g_vy(p_xlower+i-1,p_ylower+j-1,p_zlower+k-1) = global_data(c)
                    end do
                end do
            end do

            if(p_play /= pnlays) then
                indexb = localnlays
            else
                indexb = localnlays + 1
            end if
            do k = 1, indexb
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        g_vz(p_xlower+i-1,p_ylower+j-1,p_zlower+k-1) = global_data(c)
                    end do
                end do
            end do

            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        g_p(p_xlower+i-1,p_ylower+j-1,p_zlower+k-1) = global_data(c)
                    end do
                end do
            end do

            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        c = c + 1
                        g_Cf(p_xlower+i-1,p_ylower+j-1,p_zlower+k-1) = global_data(c)
                    end do
                end do
            end do

        end do

        ! export results to the raw txt files
        open(unit=10, file=fporotxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fporotxt, ' error. ', ierr
            stop
        end if
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(10, fmt='(es24.16)', iostat=ierr) g_poro(i,j,k)
                    if(ierr /= 0) then
                        print *, 'write file ', fporotxt, ' error. ', ierr
                        stop
                    end if
                end do
            end do
        end do
        close(10)

        open(unit=10, file=fKxxtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fKxxtxt, ' error. ', ierr
            stop
        end if
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(10, fmt='(es24.16)', iostat=ierr) g_Kxx(i,j,k)
                    if(ierr /= 0) then
                        print *, 'write file ', fKxxtxt, ' error. ', ierr
                        stop
                    end if
                end do
            end do
        end do
        close(10)

        open(unit=10, file=fvxtxt, status='replace', iostat=ierr)
        open(unit=11, file=fvxmidtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fvxtxt, ' error. ', ierr
            stop
        end if
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx+1
                    write(10, fmt='(es24.16)', iostat=ierr) g_vx(i,j,k)
                    if(ierr /= 0) then
                        print *, 'write file ', fvxtxt, ' error. ', ierr
                        stop
                    end if
                    if((i==nx/2+1).and.(k==nz/2)) then
                        write(11, fmt='(es12.5)', iostat=ierr) (ys(j)+ys(j+1))/2.D0
                        write(11, fmt='(es12.5)', iostat=ierr) (g_vx(i,j,k)+g_vx(i,j,k+1))/2.D0
                    end if
                end do
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
        do k = 1, nz
            do j = 1, ny+1
                do i = 1, nx
                    write(10, fmt='(es24.16)', iostat=ierr) g_vy(i,j,k)
                    if(ierr /= 0) then
                        print *, 'write file ', fvytxt, ' error. ', ierr
                        stop
                    end if
                    if((j==ny/2+1).and.(k==nz/2)) then
                        write(11, fmt='(es12.5)', iostat=ierr) (xs(i)+xs(i+1))/2.D0
                        write(11, fmt='(es12.5)', iostat=ierr) (g_vy(i,j,k)+g_vy(i,j,k+1))/2.D0
                    end if
                end do
            end do
        end do
        close(10)
        close(11)

        open(unit=10, file=fvztxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fvztxt, ' error. ', ierr
            stop
        end if
        do k = 1, nz+1
            do j = 1, ny
                do i = 1, nx
                    write(10, fmt='(es24.16)', iostat=ierr) g_vz(i,j,k)
                    if(ierr /= 0) then
                        print *, 'write file ', fvztxt, ' error. ', ierr
                        stop
                    end if
                end do
            end do
        end do
        close(10)

        open(unit=10, file=fptxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fptxt, ' error. ', ierr
            stop
        end if
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(10, fmt='(es24.16)', iostat=ierr) g_p(i,j,k)
                    if(ierr /= 0) then
                        print *, 'write file ', fptxt, ' error. ', ierr
                        stop
                    end if
                end do
            end do
        end do
        close(10)

        open(unit=10, file=fCftxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fCftxt, ' error. ', ierr
            stop
        end if
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(10, fmt='(es24.16e3)', iostat=ierr) g_Cf(i,j,k)
                    if(ierr /= 0) then
                        print *, 'write file ', fCftxt, ' error. ', ierr
                        stop
                    end if
                end do
            end do
        end do
        close(10)

        ! export results to Matlab file
        call export2Matlab()

        ! export results to tecplot file
        call export2tecplot(g_poro, g_Kxx, g_vx, g_vy, g_vz, g_p, g_Cf)

        deallocate(g_poro)
        deallocate(g_Kxx)
        deallocate(g_vx)
        deallocate(g_vy)
        deallocate(g_vz)
        deallocate(g_p)
        deallocate(g_Cf)

    end subroutine exportResults

end module DBF_exportResults
