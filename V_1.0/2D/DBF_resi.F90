
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_resi

    use DBF_model
    use DBF_globalData
    implicit none

Contains

    ! Calc Resi x-momentum
    ! The derivatives are evaluated on edges

    subroutine Resi_velx_b(velx, resi)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi        

        integer :: j

        resi(:,:) = 0.D0

        ! Applying Neumann BC
        do j = 1, localnrows
            if((pcol==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                resi(1,j) = velx(1,j) - vxBdryX0(ylower+j-1)
            end if
            if((pcol==pncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                resi(localncols+1,j) = velx(localncols+1,j) - vxBdryX1(ylower+j-1)
            end if
        end do

    end subroutine Resi_velx_b

    subroutine Resi_velx(velx, resi)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        real(kind=8) :: velxLaplace, Forchh
        real(kind=8) :: vyAv, vAbs, DVxDx, DVxDy
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        indexl = 1
        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        indexd = 1
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                ! Determine the advection term
                if((pcol==1).and.(i==1)) then
                    vyAv = 2.5D-1 * (vyBdryX0(ylower+j-1)+vyBdryX0(ylower+j)+vy(i,j)+vy(i,j+1))
                elseif((pcol==pncols).and.(i == localncols+1)) then
                    vyAv = 2.5D-1 * (vy(i-1,j)+vy(i-1,j+1)+vyBdryX1(ylower+j-1)+vyBdryX1(ylower+j))
                else
                    vyAv = 2.5D-1 * (vy(i-1,j)+vy(i-1,j+1)+vy(i,j)+vy(i,j+1))
                end if

                vAbs = dsqrt(vyAv**2.D0 + vx(i,j)**2.D0)

                if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                    DVxDx = (velx(i+1,j)/poroHarmX(i+1,j)-vxBdryX0(ylower+j-1)/poroHarmX(i,j))/hx(1)
                elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==1)) then
                    DVxDx = (velx(i+1,j)/poroHarmX(i+1,j)-velx(i,j)/poroHarmX(i,j))/hx(1)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)) then
                    DVxDx = (vxBdryX1(ylower+j-1)/poroHarmX(i,j)-velx(i-1,j)/poroHarmX(i-1,j))/hx(localncols)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==1)) then
                    DVxDx = (velx(i,j)/poroHarmX(i,j)-velx(i-1,j)/poroHarmX(i-1,j))/hx(localncols)
                else
                    if(vx(i,j)>0.D0) then
                        DVxDx = (velx(i,j)/poroHarmX(i,j)-velx(i-1,j)/poroHarmX(i-1,j))/hx(i-1)
                    else
                        DVxDx = (velx(i+1,j)/poroHarmX(i+1,j)-velx(i,j)/poroHarmX(i,j))/hx(i)
                    end if
                end if

                if((prow==1).and.(j==1)) then
                    DVxDy = (velx(i,j)/poroHarmX(i,j)-vxBdryY0(xlower+i-1)/poroHarmX(i,j))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows)) then
                    DVxDy = (vxBdryY1(xlower+i-1)/poroHarmX(i,j)-velx(i,j)/poroHarmX(i,j))/hy(j)
                else
                    if(vyAv>0.D0) then
                        DVxDy = (velx(i,j)/poroHarmX(i,j)-velx(i,j-1)/poroHarmX(i,j-1))/((hy(j-1)+hy(j))/2.D0)
                    else
                        DVxDy = (velx(i,j+1)/poroHarmX(i,j+1)-velx(i,j)/poroHarmX(i,j))/((hy(j)+hy(j+1))/2.D0)
                    end if
                end if

                resi(i,j) = - rhof/poroHarmX_old(i,j) * (vx(i,j)*DVxDx+vyAv*DVxDy)

                ! Determine the Darcian term
                if(isDarcy) then
                    if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0))).and. &!
                        (.not.((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                            resi(i,j) = resi(i,j) - visc/KxxHarm(i,j)*velx(i,j)
                    end if
                end if

                ! Determine the Brinkman term
                if(isBrinkman) then
                    if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                        velxLaplace = (velx(i+1,j)/poroHarmX(i+1,j)-vxBdryX0(ylower+j-1)/poroHarmX(i,j))/(hx(1)**2.D0)
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==1)) then
                        velxLaplace = (velx(i+1,j)/poroHarmX(i+1,j)-velx(i,j)/poroHarmX(i,j))/(hx(1)**2.D0)
                    elseif((pcol==pncols).and.(i == localncols+1).and.(isDiriX1_p(ylower+j-1)==0)) then
                        velxLaplace = (vxBdryX1(ylower+j-1)/poroHarmX(i,j)-velx(i-1,j)/poroHarmX(i-1,j))/(hx(localncols)**2.D0)
                    elseif((pcol==pncols).and.(i == localncols+1).and.(isDiriX1_p(ylower+j-1)==1)) then
                        velxLaplace = (velx(i,j)/poroHarmX(i,j)-velx(i-1,j)/poroHarmX(i-1,j))/(hx(localncols)**2.D0)
                    else
                        velxLaplace = ((velx(i+1,j)/poroHarmX(i+1,j)-velx(i,j)/poroHarmX(i,j))/hx(i) - &!
                            (velx(i,j)/poroHarmX(i,j)-velx(i-1,j)/poroHarmX(i-1,j))/hx(i-1)) / ((hx(i-1)+hx(i))/2.D0)
                    end if

                    if((prow==1).and.(j == 1)) then
                        velxLaplace = velxLaplace + ((velx(i,j+1)/poroHarmX(i,j+1)-velx(i,j)/poroHarmX(i,j))/ &!
                            ((hy(j+1)+hy(j))/2.D0) - (velx(i,j)/poroHarmX(i,j)-vxBdryY0(xlower+i-1)/ &!
                            poroHarmX(i,j))/hy(1)) / hy(j)
                    elseif((prow==pnrows).and.(j == localnrows)) then
                        velxLaplace = velxLaplace + ((vxBdryY1(xlower+i-1)/poroHarmX(i,j)-velx(i,j)/ &!
                            poroHarmX(i,j))/hy(localnrows) - (velx(i,j)/poroHarmX(i,j)-velx(i,j-1)/ &!
                            poroHarmX(i,j-1))/((hy(j)+hy(j-1))/2.D0)) / hy(j)
                    else
                        velxLaplace = velxLaplace + ((velx(i,j+1)/poroHarmX(i,j+1)-velx(i,j)/ &!
                            poroHarmX(i,j))/((hy(j+1)+hy(j))/2.D0) - (velx(i,j)/poroHarmX(i,j)- &!
                            velx(i,j-1)/poroHarmX(i,j-1))/((hy(j)+hy(j-1))/2.D0)) / hy(j)
                    end if

                    resi(i,j) = resi(i,j) + visc*velxLaplace
                end if

                ! Determine the Forchhimer term
                if(isForchheimer) then
                    if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0))).and. &!
                        (.not.((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                        Forchh = 1.75D0/dsqrt(1.5D2*poroHarmX(i,j)**3.D0)
                        resi(i,j) = resi(i,j) - (Forchh*rhof/dsqrt(KxxHarm(i,j)))*(vAbs**alpha)*velx(i,j)
                    end if
                end if

                ! Determine the time term
                resi(i,j) = resi(i,j) - rhof*(velx(i,j)/poroHarmX(i,j)-vx(i,j)/poroHarmX_old(i,j))/(timeEnd/nt)

            end do
        end do

    end subroutine Resi_velx

    subroutine Resi_dpdx(pres, resi)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        real(kind=8), dimension(:,:), allocatable :: pPad
        real(kind=8), dimension(:), allocatable :: xc
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        ! the index for equation
        indexl = 1
        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        indexd = 1
        indexu = localnrows
        
        allocate(pPad(0:localncols+1,1:localnrows))
        allocate(xc(0:localncols+1))

        xc(1:localncols) = (xs(xlower:xlower+localncols-1)+xs(xlower+1:xlower+localncols))/2.D0
        if(pcol /= 1) then
            xc(0) = (xs(xlower-1)+xs(xlower))/2.D0
        else
            xc(0) = xs(xlower)
        end if
        if(pcol /= pncols) then
            xc(localncols+1) = (xs(xlower+localncols)+xs(xlower+localncols+1))/2.D0
        else
            xc(localncols+1) = xs(xlower+localncols)
        end if

        pPad(1:localncols,1:localnrows) = pres(1:localncols,1:localnrows)
        ! insert BC
        if(pcol /= 1) then
            pPad(0,1:localnrows) = pres(0,1:localnrows)
        else
            pPad(0,1:localnrows) = pBdryX0(ylower:ylower+localnrows-1) * isDiriX0_p(ylower:ylower+localnrows-1)
        end if
        if(pcol /= pncols) then
            pPad(localncols+1,1:localnrows) = pres(localncols+1,1:localnrows)
        else
            pPad(localncols+1,1:localnrows) = pBdryX1(ylower:ylower+localnrows-1) * isDiriX1_p(ylower:ylower+localnrows-1)
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                ! Applying Neumann BC
                if(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and.(i==localncols+1) &!
                    .and.(isDiriX1_p(ylower+j-1)==0))) then
                    resi(i,j) = 0.D0
                else
                    resi(i,j) = -(pPad(i,j) - pPad(i-1,j))/(xc(i) - xc(i-1))
                end if

            end do
        end do

        deallocate(pPad)
        deallocate(xc)

    end subroutine Resi_dpdx

    ! End Calc Resi x-momentum

    ! Calc Resi y-momentum
    ! The derivatives are evaluated on edges

    subroutine Resi_vely_b(vely, resi)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi

        integer :: i

        resi(:,:) = 0.D0

        ! Applying Neumann BC
        do i = 1, localncols
            if((prow==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                resi(i,1) = vely(i,1) - vyBdryY0(xlower+i-1)
            end if
            if((prow==pnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                resi(i,localnrows+1) = vely(i,localnrows+1) - vyBdryY1(xlower+i-1)
            end if
        end do
        
    end subroutine Resi_vely_b

    subroutine Resi_vely(vely, resi)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        real(kind=8) :: velyLaplace, Forchh
        real(kind=8) :: vxAv, vAbs, DVyDx, DVyDy
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        indexl = 1
        indexr = localncols
        indexd = 1
        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                ! Determine the advection term
                if((prow==1).and.(j==1)) then
                    vxAv = 2.5D-1 * (vxBdryY0(xlower+i-1)+vxBdryY0(xlower+i)+vx(i,j)+vx(i+1,j))
                elseif((prow==pnrows).and.(j==localnrows+1)) then
                    vxAv = 2.5D-1 * (vx(i,j-1)+vx(i+1,j-1)+vxBdryY1(xlower+i-1)+vxBdryY1(xlower+i))
                else
                    vxAv = 2.5D-1 * (vx(i,j-1)+vx(i+1,j-1)+vx(i,j)+vx(i+1,j))
                end if

                vAbs = dsqrt(vxAv**2.D0 + vy(i,j)**2.D0)

                if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                    DVyDy = (vely(i,j+1)/poroHarmY(i,j+1)-vyBdryY0(xlower+i-1)/poroHarmY(i,j))/hy(j)
                elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==1)) then
                    DVyDy = (vely(i,j+1)/poroHarmY(i,j+1)-vely(i,j)/poroHarmY(i,j))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)) then
                    DVyDy = (vyBdryY1(xlower+i-1)/poroHarmY(i,j)-vely(i,j-1)/poroHarmY(i,j-1))/hy(j-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==1)) then
                    DVyDy = (vely(i,j)/poroHarmY(i,j)-vely(i,j-1)/poroHarmY(i,j-1))/hy(j-1)
                else
                    if(vy(i,j)>0.D0) then
                        DVyDy = (vely(i,j)/poroHarmY(i,j)-vely(i,j-1)/poroHarmY(i,j-1))/hy(j-1)
                    else
                        DVyDy = (vely(i,j+1)/poroHarmY(i,j+1)-vely(i,j)/poroHarmY(i,j))/hy(j)
                    end if
                end if
                if((pcol==1).and.(i==1)) then
                    DVyDx = (vely(i,j)/poroHarmY(i,j)-vyBdryX0(ylower+j-1)/poroHarmY(i,j))/hx(i)
                elseif((pcol==pncols).and.(i==localncols)) then
                    DVyDx = (vyBdryX1(ylower+j-1)/poroHarmY(i,j)-vely(i,j)/poroHarmY(i,j))/hx(i)
                else
                    if(vxAv>0.D0) then
                        DVyDx = (vely(i,j)/poroHarmY(i,j)-vely(i-1,j)/poroHarmY(i-1,j))/((hx(i-1)+hx(i))/2.D0)
                    else
                        DVyDx = (vely(i+1,j)/poroHarmY(i+1,j)-vely(i,j)/poroHarmY(i,j))/((hx(i)+hx(i+1))/2.D0)
                    end if
                end if

                resi(i,j) = - rhof/poroHarmY_old(i,j) * (vy(i,j)*DVyDy+vxAv*DVyDx)

                ! Determine the Darcian term
                if(isDarcy) then
                    if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0))).and. &!
                        (.not.((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                        resi(i,j) = resi(i,j) - visc/KyyHarm(i,j)*vely(i,j)
                    end if
                end if

                ! Determine the Brinkman term
                if(isBrinkman) then
                    if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                        velyLaplace = (vely(i,j+1)/poroHarmY(i,j+1)-vyBdryY0(xlower+i-1)/poroHarmY(i,j))/(hy(j)**2.D0)
                    elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==1)) then
                        velyLaplace = (vely(i,j+1)/poroHarmY(i,j+1)-vely(i,j)/poroHarmY(i,j))/(hy(j)**2.D0)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)) then
                        velyLaplace = (vyBdryY1(xlower+i-1)/poroHarmY(i,j)-vely(i,j-1)/poroHarmY(i,j-1))/(hy(j-1)**2.D0)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==1)) then
                        velyLaplace = (vely(i,j)/poroHarmY(i,j)-vely(i,j-1)/poroHarmY(i,j-1))/(hy(j-1)**2.D0)
                    else
                        velyLaplace = ((vely(i,j+1)/poroHarmY(i,j+1)-vely(i,j)/poroHarmY(i,j))/hy(j) - &!
                            (vely(i,j)/poroHarmY(i,j)-vely(i,j-1)/poroHarmY(i,j-1))/hy(j-1)) / ((hy(j-1)+hy(j))/2.D0)
                    end if
                    if((pcol==1).and.(i==1)) then
                        velyLaplace = velyLaplace + ((vely(i+1,j)/poroHarmY(i+1,j)-vely(i,j)/poroHarmY(i,j))/ &!
                            ((hx(i+1)+hx(i))/2.D0) - (vely(i,j)/poroHarmY(i,j)-vyBdryX0(ylower+j-1)/ &!
                            poroHarmY(i,j))/hx(i)) / hx(i)
                    elseif((pcol==pncols).and.(i==localncols)) then
                        velyLaplace = velyLaplace + ((vyBdryX1(ylower+j-1)/poroHarmY(i,j)-vely(i,j)/ &!
                            poroHarmY(i,j))/hx(i) - (vely(i,j)/poroHarmY(i,j)-vely(i-1,j)/ &!
                            poroHarmY(i-1,j))/((hx(i)+hx(i-1))/2.D0)) / hx(i)
                    else
                        velyLaplace = velyLaplace + ((vely(i+1,j)/poroHarmY(i+1,j)-vely(i,j)/ &!
                            poroHarmY(i,j))/((hx(i+1)+hx(i))/2.D0) - (vely(i,j)/poroHarmY(i,j)-vely(i-1,j)/ &!
                            poroHarmY(i-1,j))/((hx(i)+hx(i-1))/2.D0)) / hx(i)
                    end if

                    resi(i,j) = resi(i,j) + visc*velyLaplace
                end if

                ! Determine the Forchheimer term
                if(isForchheimer) then
                    if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0))).and. &!
                        (.not.((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                        Forchh = 1.75D0/dsqrt(1.5D2*poroHarmY(i,j)**3.D0)
                        resi(i,j) = resi(i,j) - (Forchh*rhof/dsqrt(KyyHarm(i,j)))*(vAbs**alpha)*vely(i,j)
                    end if
                end if

                ! Determine the time term
                resi(i,j) = resi(i,j) - rhof*(vely(i,j)/poroHarmY(i,j)-vy(i,j)/poroHarmY_old(i,j))/(timeEnd/nt)

            end do
        end do

    end subroutine Resi_vely

    subroutine Resi_dpdy(pres, resi)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        real(kind=8), dimension(:,:), allocatable :: pPad
        real(kind=8), dimension(:), allocatable :: yc
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j
 
        indexl = 1
        indexr = localncols
        indexd = 1
        if(prow /= pnrows) then
            indexu = localnrows 
        else
            indexu = localnrows + 1
        end if

        allocate(pPad(1:localncols,0:localnrows+1))
        allocate(yc(0:localnrows+1))

        yc(1:localnrows) = (ys(ylower:ylower+localnrows-1)+ys(ylower+1:ylower+localnrows))/2.D0
        if(prow /= 1) then
            yc(0) = (ys(ylower-1)+ys(ylower))/2.D0
        else
            yc(0) = ys(1)
        end if
        if(prow /= pnrows) then
            yc(localnrows+1) = (ys(ylower+localnrows)+ys(ylower+localnrows+1))/2.D0
        else
            yc(localnrows+1) = ys(ylower+localnrows)
        end if

        pPad(1:localncols,1:localnrows) = pres(1:localncols,1:localnrows)
        ! insert BC
        if(prow /= 1) then
            pPad(1:localncols,0) = pres(1:localncols,0)
        else
            pPad(1:localncols,0) = pBdryY0(xlower:xlower+localncols-1) * isDiriY0_p(xlower:xlower+localncols-1)
        end if
        if(prow /= pnrows) then
            pPad(1:localncols,localnrows+1) = pres(1:localncols,localnrows+1)
        else
            pPad(1:localncols,localnrows+1) = pBdryY1(xlower:xlower+localncols-1) * isDiriY1_p(xlower:xlower+localncols-1)
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                ! Applying Neumann BC
                if(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and.(j==localnrows+1) &!
                    .and.(isDiriY1_p(xlower+i-1)==0))) then
                    resi(i,j) = 0
                else
                    resi(i,j) = -(pPad(i,j) - pPad(i,j-1))/(yc(j) - yc(j-1))
                end if

            end do
        end do

        deallocate(pPad)
        deallocate(yc)

    end subroutine Resi_dpdy

    ! End Calc Resi y-momentum

    ! Calc Resi Continuity

    subroutine Resi_dudx(velx, resi)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                    resi(i,j) = - (velx(i+1,j) - vxBdryX0(ylower+j-1)) / hx(i)
                elseif((pcol==pncols).and.(i==localncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                    resi(i,j) = - (vxBdryX1(ylower+j-1) - velx(i,j)) / hx(i)
                else
                    resi(i,j) = - (velx(i+1,j) - velx(i,j)) / hx(i)
                end if

            end do
        end do

    end subroutine Resi_dudx

    subroutine Resi_dvdy(vely, resi)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real (kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                    resi(i,j) = - (vely(i,j+1) - vyBdryY0(xlower+i-1)) / hy(j)
                elseif((prow==pnrows).and.(j==localnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                    resi(i,j) = - (vyBdryY1(xlower+i-1) - vely(i,j)) / hy(j)
                else
                    resi(i,j) = - (vely(i,j+1) - vely(i,j)) / hy(j)
                end if 

            end do
        end do

    end subroutine Resi_dvdy

    subroutine Resi_dpdt(pres, resi)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi

        resi(:,:) = - epslon * (pres(:,:)-p(:,:))/(timeEnd/nt)

    end subroutine Resi_dpdt

    ! End Calc Resi Continuity

    ! Calc Resi Concentration

    subroutine Resi_Cf(conc, resi)

        integer, dimension(:,:), pointer, intent(in) :: conc
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        integer :: indexl, indexr, indexd, indexu
        real(kind=8) :: dl, dt
        real(kind=8), dimension(:,:,:), allocatable :: Ex, Ey, Dx, Dy
        real(kind=8), dimension(:,:), allocatable :: Cfbarx, Cfbary
        real(kind=8) :: vxAv, vyAv, vmodulus
        real(kind=8) :: div1, div2
        real(kind=8), dimension(:,:), allocatable :: dCfdx, dCfdy
        real(kind=8) :: dCfdxleft, dCfdxright, dCfdxdown, dCfdxup, dCfdyleft, dCfdyright, dCfdydown, dCfdyup
        real(kind=8) :: reaction
        integer :: i, j

        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows

        allocate(Ex(indexl:indexr,indexd:indexu,2))
        allocate(Dx(indexl:indexr,indexd:indexu,2))
        allocate(Cfbarx(indexl:indexr,indexd:indexu))

        do j = indexd, indexu
            do i = indexl, indexr

                if((pcol==1).and.(i==1)) then
                    vyAv = 2.5D-1 * (vyBdryX0(ylower+j-1)+vyBdryX0(ylower+j)+vy(i,j)+vy(i,j+1))
                elseif((pcol==pncols).and.(i == localncols+1)) then
                    vyAv = 2.5D-1 * (vy(i-1,j)+vy(i-1,j+1)+vyBdryX1(ylower+j-1)+vyBdryX1(ylower+j))
                else
                    vyAv = 2.5D-1 * (vy(i-1,j)+vy(i-1,j+1)+vy(i,j)+vy(i,j+1))
                end if

                vmodulus = dsqrt((vx(i,j)**2.D0+vyAv**2.D0)*1.D0)

                if(vmodulus /= 0.D0) then
                    Ex(i,j,1) = vx(i,j)**2.D0 / vmodulus**2.D0
                    Ex(i,j,2) = vx(i,j)*vyAv / vmodulus**2.D0  ! the second column and the first row of the matrix
                else
                    Ex(i,j,1) = 0.D0
                    Ex(i,j,2) = 0.D0
                end if

                dl = alphaOS*dm + 2.D0*lamdaX*vmodulus*radiusInit*(1-poroHarmXInit(i,j))/ &!
                    (poroHarmXInit(i,j)*(1-poroHarmX(i,j)))

                dt = alphaOS*dm + 2.D0*lamdaT*vmodulus*radiusInit*(1-poroHarmXInit(i,j))/ &!
                    (poroHarmXInit(i,j)*(1-poroHarmX(i,j)))

                Dx(i,j,1) = (dm+vmodulus*dt) + vmodulus*(dl-dt)*Ex(i,j,1)
                Dx(i,j,2) = vmodulus*(dl-dt)*Ex(i,j,2)

            end do
        end do

        do j = indexd, indexu
            do i = indexl, indexr

                if((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1)==0)) then
                    Cfbarx(i,j) = conc(i,j)
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1)==1)) then
                    Cfbarx(i,j) = CfBdryX0(ylower+j-1)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1)==0)) then
                    Cfbarx(i,j) = conc(i-1,j)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1)==1)) then
                    Cfbarx(i,j) = CfBdryX1(ylower+j-1)
                elseif(vx(i,j) > 0.D0) then
                    Cfbarx(i,j) = conc(i-1,j)
                else
                    Cfbarx(i,j) = conc(i,j)
                end if

            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1

        allocate(Ey(indexl:indexr,indexd:indexu,2))
        allocate(Dy(indexl:indexr,indexd:indexu,2))
        allocate(Cfbary(indexl:indexr,indexd:indexu))

        do j = indexd, indexu
            do i = indexl, indexr

                if((prow==1).and.(j==1)) then
                    vxAv = 2.5D-1 * (vxBdryY0(xlower+i-1)+vxBdryY0(xlower+i)+vx(i,j)+vx(i+1,j))
                elseif((prow==pnrows).and.(j==localnrows+1)) then
                    vxAv = 2.5D-1 * (vx(i,j-1)+vx(i+1,j-1)+vxBdryY1(xlower+i-1)+vxBdryY1(xlower+i))
                else
                    vxAv = 2.5D-1 * (vx(i,j-1)+vx(i+1,j-1)+vx(i,j)+vx(i+1,j))
                end if

                vmodulus = dsqrt((vy(i,j)**2.D0+vxAv**2.D0)*1.D0)

                if(vmodulus /= 0.D0) then
                    Ey(i,j,1) = vy(i,j)*vxAv / vmodulus**2.D0
                    Ey(i,j,2) = vy(i,j)**2.D0 / vmodulus**2.D0
                else
                    Ey(i,j,1) = 0.D0
                    Ey(i,j,2) = 0.D0
                end if

                dl = alphaOS*dm + 2.D0*lamdaX*vmodulus*radiusInit*(1-poroHarmYInit(i,j))/ &!
                    (poroHarmYInit(i,j)*(1-poroHarmY(i,j)))

                dt = alphaOS*dm + 2.D0*lamdaT*vmodulus*radiusInit*(1-poroHarmYInit(i,j))/ &!
                    (poroHarmYInit(i,j)*(1-poroHarmY(i,j)))

                Dy(i,j,1) = vmodulus*(dl-dt)*Ey(i,j,1)
                Dy(i,j,2) = (dm+vmodulus*dt) + vmodulus*(dl-dt)*Ey(i,j,2)

            end do
        end do

        do j = indexd, indexu
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1)==0)) then
                    Cfbary(i,j) = conc(i,j)
                elseif((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1)==1)) then
                    Cfbary(i,j) = CfBdryY0(xlower+i-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1)==0)) then
                    Cfbary(i,j) = conc(i,j-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1)==1)) then
                    Cfbary(i,j) = CfBdryY1(xlower+i-1)
                elseif(vy(i,j) > 0.D0) then
                    Cfbary(i,j) = conc(i,j-1)
                else
                    Cfbary(i,j) = conc(i,j)
                end if

            end do
        end do

        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if
        allocate(dCfdx(1:localncols+1,indexd:indexu))

        do j = indexd, indexu
            do i = 1, localncols+1

                if((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1)==0)) then
                    dCfdx(i,j) = 0.D0
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1)==1)) then
                    dCfdx(i,j) = (conc(i,j)-CfBdryX0(ylower+j-1))/hx(i)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1)==0)) then
                    dCfdx(i,j) = 0.D0
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1)==1)) then
                    dCfdx(i,j) = (CfBdryX1(ylower+j-1)-conc(i-1,j))/hx(i-1)
                else
                    dCfdx(i,j) = (conc(i,j)-conc(i-1,j))/((hx(i)+hx(i-1))/2.D0)
                end if

            end do
        end do

        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if

        allocate(dCfdy(indexl:indexr,1:localnrows+1))

        do j = 1, localnrows+1
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1)==0)) then
                    dCfdy(i,j) = 0.D0
                elseif((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1)==1)) then
                    dCfdy(i,j) = (conc(i,j)-CfBdryY0(xlower+i-1))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1)==0)) then
                    dCfdy(i,j) = 0.D0
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1)==1)) then
                    dCfdy(i,j) = (CfBdryY1(xlower+i-1)-conc(i,j-1))/hy(j-1)
                else
                    dCfdy(i,j) = (conc(i,j)-conc(i,j-1))/((hy(j)+hy(j-1))/2.D0)
                end if

            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols

                div1 = (vx(i+1,j)*Cfbarx(i+1,j) - vx(i,j)*Cfbarx(i,j))/hx(i) + &!
                    (vy(i,j+1)*Cfbary(i,j+1) - vy(i,j)*Cfbary(i,j))/hy(j)

                dCfdxleft = dCfdx(i,j)
                dCfdxright = dCfdx(i+1,j)
                if((prow==1).and.(j==1)) then
                    dCfdxdown = (dCfdx(i,j)+dCfdx(i+1,j)) / 2.D0
                else
                    dCfdxdown = (dCfdx(i,j)+dCfdx(i+1,j)+dCfdx(i,j-1)+dCfdx(i+1,j-1)) / 4.D0
                end if
                if((prow==pnrows).and.(j==localnrows)) then
                    dCfdxup = (dCfdx(i,j)+dCfdx(i+1,j)) / 2.D0
                else
                    dCfdxup = (dCfdx(i,j)+dCfdx(i+1,j)+dCfdx(i,j+1)+dCfdx(i+1,j+1)) / 4.D0
                end if

                if((pcol==1).and.(i==1)) then
                    dCfdyleft = (dCfdy(i,j)+dCfdy(i,j+1)) / 2.D0
                else
                    dCfdyleft = (dCfdy(i-1,j)+dCfdy(i-1,j+1)+dCfdy(i,j)+dCfdy(i,j+1)) / 4.D0
                end if
                if((pcol==pncols).and.(i==localncols)) then
                    dCfdyright = (dCfdy(i,j)+dCfdy(i,j+1)) / 2.D0
                else
                    dCfdyright = (dCfdy(i,j)+dCfdy(i,j+1)+dCfdy(i+1,j)+dCfdy(i+1,j+1)) / 4.D0
                end if
                dCfdydown = dCfdy(i,j)
                dCfdyup = dCfdy(i,j+1)

                div2 = (dCfdxright*poroHarmX(i+1,j)*Dx(i+1,j,1) - dCfdxleft*poroHarmX(i,j)*Dx(i,j,1))/hx(i) + &!
                    (dCfdyright*poroHarmX(i+1,j)*Dx(i+1,j,2) - dCfdyleft*poroHarmX(i,j)*Dx(i,j,2))/hx(i) + &!
                    (dCfdxup*poroHarmY(i,j+1)*Dy(i,j+1,1) - dCfdxdown*poroHarmY(i,j)*Dy(i,j,1))/hy(j) + &!
                    (dCfdyup*poroHarmY(i,j+1)*Dy(i,j+1,2) - dCfdydown*poroHarmY(i,j)*Dy(i,j,2))/hy(j)

                reaction = av(i,j)*conc(i,j)*kc(i,j)*ks/(kc(i,j)+ks)

                resi(i,j) = src(xlower+i-1,ylower+j-1)-reaction+div2-div1-(conc(i,j)*poro(i,j)-Cf(i,j)*poro_old(i,j))/ &!
                    (timeEnd/nt)

            end do
        end do

        deallocate(Ex)
        deallocate(Ey)
        deallocate(Dx)
        deallocate(Dy)
        deallocate(Cfbarx)
        deallocate(Cfbary)
        deallocate(dCfdx)
        deallocate(dCfdy)

    end subroutine Resi_Cf

    ! End Calc Resi Concentration

end module DBF_resi
