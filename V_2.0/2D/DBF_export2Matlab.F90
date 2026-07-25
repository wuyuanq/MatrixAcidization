
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_export2Matlab
    
    use DBF_model
    use DBF_globalData
    implicit none

contains

    subroutine export2Matlab()

        character(len=30) :: fpmatlabm
        character(len=8) :: charLx, charLy
        character(len=5) :: charnx, charny
        character(len=11) :: chartimeEnd, charnt, chart, charnframe
        logical :: alive
        integer :: ierr

        fpmatlabm = trim(adjustl(soludoc))//"/matlabplot.m"

        open(unit=10, file=fpmatlabm, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fpmatlabm, ' error. ', ierr
            stop
        end if

        write(10, fmt="(a)") "path('/Users/yuanqingwu/research/MatrixAcidization/V_2.0/2D', path);"
        write(charLx,'(f8.4)') Lx
        write(10, fmt="(a)") "Lx = "//trim(adjustl(charLx))//";"
        write(charLy,'(f8.4)') Ly
        write(10, fmt="(a)") "Ly = "//trim(adjustl(charLy))//";"
        write(chartimeEnd,'(f11.1)') timeEnd
        write(10, fmt="(a)") "timeEnd = "//trim(adjustl(chartimeEnd))//";"
        write(charnx,'(i5)') nx
        write(10, fmt="(a)") "nx = "//trim(adjustl(charnx))//";"
        write(charny,'(i5)') ny
        write(10, fmt="(a)") "ny = "//trim(adjustl(charny))//";"
        write(charnt,'(i10)') nt
        write(10, fmt="(a)") "nt = "//trim(adjustl(charnt))//";"
        write(10, fmt="(a)") "model.nx = nx;"
        write(10, fmt="(a)") "model.ny = ny;"
        write(10, fmt="(a)") "model.nt = nt;"
        write(10, fmt="(a)") "model.xs = (0:nx)*Lx/nx;"
        write(10, fmt="(a)") "model.ys = (0:ny)*Ly/ny;"
        write(10, fmt="(a)") "model.ts = (0:nt)*timeEnd/nt;"
        write(10, fmt="(a)") "model.timeEnd = timeEnd;"
        write(charnframe,'(i10)') NUMFRAME
        write(10, fmt="(a)") "NUMFRAME = "//trim(adjustl(charnframe))//";"
        write(chart,'(i10)') t
        write(10, fmt="(a)") "fporotxt = 'soln_poro_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "fKxxtxt = 'soln_Kxx_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "fvxtxt = 'soln_vx_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "fvytxt = 'soln_vy_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "fptxt = 'soln_p_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "fCftxt = 'soln_Cf_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "fTemtxt = 'soln_Tem_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "fporohistxt = 'his_poro_avg.txt';"
        write(10, fmt="(a)") "fkxxhistxt = 'his_Kxx_avg.txt';"
        write(10, fmt="(a)") "favhistxt = 'his_av_avg.txt';"
        write(10, fmt="(a)") "fphistxt = 'his_p_avg.txt';"
        write(10, fmt="(a)") "fcfhistxt = 'his_Cf_avg.txt';"
        write(10, fmt="(a)") "fTemhistxt = 'his_Tem_avg.txt';"
        write(10, fmt="(a)") "fqhistxt = 'his_q_avg.txt';"
        write(10, fmt="(a)") "flphistxt = 'his_lp_avg.txt';"
        write(10, fmt="(a)") "model.soludoc = 'matlabplots';"
        inquire(file = trim(adjustl(soludoc))//'/matlabplots', exist = alive)
        if(.not.alive) then
            call system("mkdir "//trim(adjustl(soludoc))//"/matlabplots")
        end if
        write(10, fmt="(a)") "DBF_plot(model, fporotxt, fKxxtxt, fvxtxt, fvytxt, fptxt, fCftxt, fTemtxt, NUMFRAME);"
        write(10, fmt="(a)") &!
            "DBF_his(model, fporohistxt, fKxxhistxt, favhistxt, fphistxt, fCfhistxt, fTemhistxt, fqhistxt, flphistxt);"
        write(10, fmt="(a)") "rmpath('/Users/yuanqingwu/research/MatrixAcidization/V_2.0/2D');"

        close(10)

    end subroutine export2Matlab

end module DBF_export2Matlab
