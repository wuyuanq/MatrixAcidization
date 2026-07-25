
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

    ! Calc Resi x-velocity
    ! The derivatives are evaluated on edges

    subroutine Resi_xmom_vx_b(velx, resi)

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

    end subroutine Resi_xmom_vx_b

    subroutine Resi_xmom_vx(velx, resi)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi

        real(kind=8) :: dpdx, vyAv, vAbs, Forchh, ForchhTerm, vxD
        real(kind=8) :: DVxDx, DVxDy, velxLaplace
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

                if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                    dpdx = 0.D0
                elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==1)) then
                    dpdx = (p(i,j)-pBdryX0(ylower+j-1)) / hx(i)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)) then
                    dpdx = 0.D0
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==1)) then
                    dpdx = (pBdryX1(ylower+j-1)-p(i-1,j)) / hx(i-1)
                else
                    dpdx = (p(i,j)-p(i-1,j))/((hx(i)+hx(i-1))/2.D0)
                end if

                if((pcol==1).and.(i==1)) then
                    vyAv = (vyBdryX0(ylower+j-1)+vyBdryX0(ylower+j)+vy(i,j)+vy(i,j+1)) / 4.D0
                elseif((pcol==pncols).and.(i == localncols+1)) then
                    vyAv = (vy(i-1,j)+vy(i-1,j+1)+vyBdryX1(ylower+j-1)+vyBdryX1(ylower+j)) / 4.D0
                else
                    vyAv = (vy(i-1,j)+vy(i-1,j+1)+vy(i,j)+vy(i,j+1)) / 4.D0
                end if

                vAbs = dsqrt(vx(i,j)**2.D0 + vyAv**2.D0)

                Forchh = 1.75D0/dsqrt(1.5D2*poroHarmX(i,j)**3.D0)
                ForchhTerm = rhof*Forchh/dsqrt(KxxHarm(i,j))*vAbs
                vxD = - (dpdx-rhof*gravX)/(ForchhTerm+visc/KxxHarm(i,j))

                ! Determine the time term
                resi(i,j) = - rhof*(velx(i,j)/poroHarmX(i,j)-vx(i,j)/poroHarmX_old(i,j))/(timeEnd/nt)

                ! Determine the Forchhimer term
                if(isForchheimer) then
                    if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0))).and. &!
                        (.not.((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                        resi(i,j) = resi(i,j) - &!
                            (visc/KxxHarm(i,j)+Forchh*rhof/dsqrt(KxxHarm(i,j))*vAbs)*(velx(i,j)-vxD)
                    end if
                end if

                ! Determine the advection term
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

                resi(i,j) = resi(i,j) - rhof/poroHarmX_old(i,j) * (vx(i,j)*DVxDx+vyAv*DVxDy)

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

            end do
        end do

    end subroutine Resi_xmom_vx

    ! End Calc Resi x-velocity

    ! Calc Resi y-velocity
    ! The derivatives are evaluated on edges

    subroutine Resi_ymom_vy_b(vely, resi)

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
        
    end subroutine Resi_ymom_vy_b

    subroutine Resi_ymom_vy(vely, resi)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        real(kind=8) :: dpdy, vxAv, vAbs, Forchh, ForchhTerm, vyD
        real(kind=8) :: DVyDx, DVyDy, velyLaplace
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

                if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                    dpdy = 0.D0
                elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==1)) then
                    dpdy = (p(i,j)-pBdryY0(xlower+i-1)) / hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)) then
                    dpdy = 0.D0
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==1)) then
                    dpdy = (pBdryY1(xlower+i-1)-p(i,j-1)) / hy(j-1)
                else
                    dpdy = (p(i,j)-p(i,j-1))/((hy(j)+hy(j-1))/2.D0)
                end if

                if((prow==1).and.(j==1)) then
                    vxAv = (vxBdryY0(xlower+i-1)+vxBdryY0(xlower+i)+vx(i,j)+vx(i+1,j)) / 4.D0
                elseif((prow==pnrows).and.(j==localnrows+1)) then
                    vxAv = (vx(i,j-1)+vx(i+1,j-1)+vxBdryY1(xlower+i-1)+vxBdryY1(xlower+i)) / 4.D0
                else
                    vxAv = (vx(i,j-1)+vx(i+1,j-1)+vx(i,j)+vx(i+1,j)) / 4.D0
                end if

                vAbs = dsqrt(vxAv**2.D0 + vy(i,j)**2.D0)

                Forchh = 1.75D0/dsqrt(1.5D2*poroHarmY(i,j)**3.D0)
                ForchhTerm = rhof*Forchh/dsqrt(KyyHarm(i,j))*vAbs
                vyD = - (dpdy-rhof*gravY)/(ForchhTerm+visc/KyyHarm(i,j))

                ! Determine the time term
                resi(i,j) = - rhof*(vely(i,j)/poroHarmY(i,j)-vy(i,j)/poroHarmY_old(i,j))/(timeEnd/nt)

                ! Determine the Forchhimer term
                if(isForchheimer) then
                    if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0))).and. &!
                        (.not.((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                        resi(i,j) = resi(i,j) - &!
                            (visc/KyyHarm(i,j)+Forchh*rhof/dsqrt(KyyHarm(i,j))*vAbs)*(vely(i,j)-vyD)
                    end if
                end if

                ! Determine the advection term
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

                resi(i,j) = resi(i,j) - rhof/poroHarmY_old(i,j) * (vy(i,j)*DVyDy+vxAv*DVyDx)

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

            end do
        end do

    end subroutine Resi_ymom_vy

    ! End Calc Resi y-velocity

    ! Calc Resi Mass

    subroutine Resi_mass_p(pres, resi)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi

        integer :: indexl, indexr, indexd, indexu
        real(kind=8), dimension(:,:), allocatable :: dpdx, dpdy
        real(kind=8) :: vyAv_l, vyAv_r, vxAv_d, vxAv_u
        real(kind=8) :: vAbs_l, vAbs_r, vAbs_d, vAbs_u
        real(kind=8) :: div_l, div_r, div_d, div_u
        real(kind=8) :: ForchhTerm_l, ForchhTerm_r, ForchhTerm_d, ForchhTerm_u, div
        integer :: i, j

        indexd = 1
        indexu = localnrows
        allocate(dpdx(1:localncols+1,indexd:indexu))

        do j = indexd, indexu
            do i = 1, localncols+1

                if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                    dpdx(i,j) = 0.D0
                elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==1)) then
                    dpdx(i,j) = (pres(i,j)-pBdryX0(ylower+j-1))/hx(i)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)) then
                    dpdx(i,j) = 0.D0
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==1)) then
                    dpdx(i,j) = (pBdryX1(ylower+j-1)-pres(i-1,j))/hx(i-1)
                else
                    dpdx(i,j) = (pres(i,j)-pres(i-1,j))/((hx(i)+hx(i-1))/2.D0)
                end if

            end do
        end do

        indexl = 1
        indexr = localncols
        allocate(dpdy(indexl:indexr,1:localnrows+1))

        do j = 1, localnrows+1
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                    dpdy(i,j) = 0.D0
                elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==1)) then
                    dpdy(i,j) = (pres(i,j)-pBdryY0(xlower+i-1))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)) then
                    dpdy(i,j) = 0.D0
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==1)) then
                    dpdy(i,j) = (pBdryY1(xlower+i-1)-pres(i,j-1))/hy(j-1)
                else
                    dpdy(i,j) = (pres(i,j)-pres(i,j-1))/((hy(j)+hy(j-1))/2.D0)
                end if

            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols

                if((pcol==1).and.(i==1)) then
                    vyAv_l = (vyBdryX0(ylower+j-1)+vyBdryX0(ylower+j)+vy(i,j)+vy(i,j+1)) / 4.D0
                else
                    vyAv_l = (vy(i-1,j)+vy(i-1,j+1)+vy(i,j)+vy(i,j+1)) / 4.D0
                end if
                if((pcol==pncols).and.(i == localncols)) then
                    vyAv_r = (vy(i,j)+vy(i,j+1)+vyBdryX1(ylower+j-1)+vyBdryX1(ylower+j)) / 4.D0
                else
                    vyAv_r = (vy(i+1,j)+vy(i+1,j+1)+vy(i,j)+vy(i,j+1)) / 4.D0
                end if
                if((prow==1).and.(j==1)) then
                    vxAv_d = (vxBdryY0(xlower+i-1)+vxBdryY0(xlower+i)+vx(i,j)+vx(i+1,j)) / 4.D0
                else
                    vxAv_d = (vx(i,j-1)+vx(i+1,j-1)+vx(i,j)+vx(i+1,j)) / 4.D0
                end if
                if((prow==pnrows).and.(j==localnrows)) then
                    vxAv_u = (vx(i,j)+vx(i+1,j)+vxBdryY1(xlower+i-1)+vxBdryY1(xlower+i)) / 4.D0
                else
                    vxAv_u = (vx(i,j+1)+vx(i+1,j+1)+vx(i,j)+vx(i+1,j)) / 4.D0
                end if
 
                vAbs_l = dsqrt(vx(i,j)**2.D0 + vyAv_l**2.D0)
                vAbs_r = dsqrt(vx(i+1,j)**2.D0 + vyAv_r**2.D0)
                vAbs_d = dsqrt(vxAv_d**2.D0 + vy(i,j)**2.D0)
                vAbs_u = dsqrt(vxAv_u**2.D0 + vy(i,j+1)**2.D0)

                ForchhTerm_l = rhof*1.75D0/dsqrt(1.5D2*poroHarmX(i,j)**3.D0)/dsqrt(KxxHarm(i,j))*vAbs_l
                ForchhTerm_r = rhof*1.75D0/dsqrt(1.5D2*poroHarmX(i+1,j)**3.D0)/dsqrt(KxxHarm(i+1,j))*vAbs_r
                ForchhTerm_d = rhof*1.75D0/dsqrt(1.5D2*poroHarmY(i,j)**3.D0)/dsqrt(KyyHarm(i,j))*vAbs_d
                ForchhTerm_u = rhof*1.75D0/dsqrt(1.5D2*poroHarmY(i,j+1)**3.D0)/dsqrt(KyyHarm(i,j+1))*vAbs_u

                if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                    div_l = vxBdryX0(ylower+j-1)
                else
                    div_l = (-dpdx(i,j)+rhof*gravX)/(ForchhTerm_l+visc/KxxHarm(i,j))
                end if
                if((pcol==pncols).and.(i==localncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                    div_r = vxBdryX1(ylower+j-1)
                else
                    div_r = (-dpdx(i+1,j)+rhof*gravX)/(ForchhTerm_r+visc/KxxHarm(i+1,j))
                end if
                if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                    div_d = vyBdryY0(xlower+i-1)
                else
                    div_d = (-dpdy(i,j)+rhof*gravY)/(ForchhTerm_d+visc/KyyHarm(i,j))
                end if
                if((prow==pnrows).and.(j==localnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                    div_u = vyBdryY1(xlower+i-1)
                else
                    div_u = (-dpdy(i,j+1)+rhof*gravY)/(ForchhTerm_u+visc/KyyHarm(i,j+1))
                end if

                div = (div_r-div_l)/hx(i) + (div_u-div_d)/hy(j)

                resi(i,j) = src(xlower+i-1,ylower+j-1)-div-(poro(i,j)-poro_old(i,j))/(timeEnd/nt)

            end do
        end do

        deallocate(dpdx)
        deallocate(dpdy)

    end subroutine Resi_mass_p

    ! End Calc Resi Mass

    ! Calc Resi Concentration

    subroutine Resi_concen_cf(conc, resi)

        integer, dimension(:,:), pointer, intent(in) :: conc
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        integer :: indexl, indexr, indexd, indexu
        real(kind=8) :: dl, dt
        real(kind=8), dimension(:,:,:), allocatable :: Ex, Ey, Dx, Dy
        real(kind=8), dimension(:,:), allocatable :: Cfbarx, Cfbary
        real(kind=8) :: vxAv, vyAv, vmodulus
        real(kind=8) :: dmbarx, dmbary
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
                    vyAv = (vyBdryX0(ylower+j-1)+vyBdryX0(ylower+j)+vy(i,j)+vy(i,j+1)) / 4.D0
                elseif((pcol==pncols).and.(i == localncols+1)) then
                    vyAv = (vy(i-1,j)+vy(i-1,j+1)+vyBdryX1(ylower+j-1)+vyBdryX1(ylower+j)) / 4.D0
                else
                    vyAv = (vy(i-1,j)+vy(i-1,j+1)+vy(i,j)+vy(i,j+1)) / 4.D0
                end if

                vmodulus = dsqrt((vx(i,j)**2.D0+vyAv**2.D0)*1.D0)

                if(vmodulus /= 0.D0) then
                    Ex(i,j,1) = vx(i,j)**2.D0 / vmodulus**2.D0
                    Ex(i,j,2) = vx(i,j)*vyAv / vmodulus**2.D0  ! the second column and the first row of the matrix
                else
                    Ex(i,j,1) = 0.D0
                    Ex(i,j,2) = 0.D0
                end if

                if((pcol==1).and.(i==1)) then
                    dmbarx = dm(i,j)
                elseif((pcol==pncols).and.(i==localncols+1)) then
                    dmbarx = dm(i-1,j)
                elseif(vx(i,j) > 0.D0) then
                    dmbarx = dm(i-1,j)
                else
                    dmbarx = dm(i,j)
                end if

                dl = alphaOS*dmbarx + 2.D0*lamdaX*vmodulus*radiusInit*(1.D0-poroHarmXInit(i,j))/ &!
                    (poroHarmXInit(i,j)*(1.D0-poroHarmX(i,j)))
                dt = alphaOS*dmbarx + 2.D0*lamdaT*vmodulus*radiusInit*(1.D0-poroHarmXInit(i,j))/ &!
                    (poroHarmXInit(i,j)*(1.D0-poroHarmX(i,j)))
                Dx(i,j,1) = (dmbarx+vmodulus*dt) + vmodulus*(dl-dt)*Ex(i,j,1)
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
                    vxAv = (vxBdryY0(xlower+i-1)+vxBdryY0(xlower+i)+vx(i,j)+vx(i+1,j)) / 4.D0
                elseif((prow==pnrows).and.(j==localnrows+1)) then
                    vxAv = (vx(i,j-1)+vx(i+1,j-1)+vxBdryY1(xlower+i-1)+vxBdryY1(xlower+i)) / 4.D0
                else
                    vxAv = (vx(i,j-1)+vx(i+1,j-1)+vx(i,j)+vx(i+1,j)) / 4.D0
                end if

                vmodulus = dsqrt((vy(i,j)**2.D0+vxAv**2.D0)*1.D0)

                if(vmodulus /= 0.D0) then
                    Ey(i,j,1) = vy(i,j)*vxAv / vmodulus**2.D0
                    Ey(i,j,2) = vy(i,j)**2.D0 / vmodulus**2.D0
                else
                    Ey(i,j,1) = 0.D0
                    Ey(i,j,2) = 0.D0
                end if

                if((prow==1).and.(j==1)) then
                    dmbary = dm(i,j)
                elseif((prow==pnrows).and.(j==localnrows+1)) then
                    dmbary = dm(i,j-1)
                elseif(vy(i,j) > 0.D0) then
                    dmbary = dm(i,j-1)
                else
                    dmbary = dm(i,j)
                end if

                dl = alphaOS*dmbary + 2.D0*lamdaX*vmodulus*radiusInit*(1.D0-poroHarmYInit(i,j))/ &!
                    (poroHarmYInit(i,j)*(1.D0-poroHarmY(i,j)))
                dt = alphaOS*dmbary + 2.D0*lamdaT*vmodulus*radiusInit*(1.D0-poroHarmYInit(i,j))/ &!
                    (poroHarmYInit(i,j)*(1.D0-poroHarmY(i,j)))
                Dy(i,j,1) = vmodulus*(dl-dt)*Ey(i,j,1)
                Dy(i,j,2) = (dmbary+vmodulus*dt) + vmodulus*(dl-dt)*Ey(i,j,2)

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

                reaction = av(i,j)*conc(i,j)*kc(i,j)*ks(i,j)/(kc(i,j)+ks(i,j))

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

    end subroutine Resi_concen_cf

    ! End Calc Resi Concentration

    ! Calc Resi Temperature

    subroutine Resi_ener_tem(tempe, resi)

        integer, dimension(:,:), pointer, intent(in) :: tempe
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi

        integer :: indexl, indexr, indexd, indexu
        real(kind=8), dimension(:,:), allocatable :: Tembarx, Tembary
        real(kind=8), dimension(:,:), allocatable :: dTemdx, dTemdy
        real(kind=8) :: dvdt, tp, vmodulus, v1, v2, Forchh, HrT, arh
        real(kind=8) :: div1, div2, div3
        integer :: i, j

        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows
        allocate(Tembarx(indexl:indexr,indexd:indexu))

        do j = indexd, indexu
            do i = indexl, indexr

                if((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1)==0)) then
                    Tembarx(i,j) = tempe(i,j)
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1)==1)) then
                    Tembarx(i,j) = TemBdryX0(ylower+j-1)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1)==0)) then
                    Tembarx(i,j) = tempe(i-1,j)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1)==1)) then
                    Tembarx(i,j) = TemBdryX1(ylower+j-1)
                elseif(vx(i,j) > 0.D0) then
                    Tembarx(i,j) = tempe(i-1,j)
                else
                    Tembarx(i,j) = tempe(i,j)
                end if

            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1
        allocate(Tembary(indexl:indexr,indexd:indexu))

        do j = indexd, indexu
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1)==0)) then
                    Tembary(i,j) = tempe(i,j)
                elseif((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1)==1)) then
                    Tembary(i,j) = TemBdryY0(xlower+i-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1)==0)) then
                    Tembary(i,j) = tempe(i,j-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1)==1)) then
                    Tembary(i,j) = TemBdryY1(xlower+i-1)
                elseif(vy(i,j) > 0.D0) then
                    Tembary(i,j) = tempe(i,j-1)
                else
                    Tembary(i,j) = tempe(i,j)
                end if

            end do
        end do

        indexd = 1
        indexu = localnrows
        allocate(dTemdx(1:localncols+1,indexd:indexu))

        do j = indexd, indexu
            do i = 1, localncols+1

                if((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1)==0)) then
                    dTemdx(i,j) = 0.D0
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1)==1)) then
                    dTemdx(i,j) = (tempe(i,j)-TemBdryX0(ylower+j-1))/hx(i)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1)==0)) then
                    dTemdx(i,j) = 0.D0
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1)==1)) then
                    dTemdx(i,j) = (TemBdryX1(ylower+j-1)-tempe(i-1,j))/hx(i-1)
                else
                    dTemdx(i,j) = (tempe(i,j)-tempe(i-1,j))/((hx(i)+hx(i-1))/2.D0)
                end if

            end do
        end do

        indexl = 1
        indexr = localncols
        allocate(dTemdy(indexl:indexr,1:localnrows+1))

        do j = 1, localnrows+1
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1)==0)) then
                    dTemdy(i,j) = 0.D0
                elseif((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1)==1)) then
                    dTemdy(i,j) = (tempe(i,j)-TemBdryY0(xlower+i-1))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1)==0)) then
                    dTemdy(i,j) = 0.D0
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1)==1)) then
                    dTemdy(i,j) = (TemBdryY1(xlower+i-1)-tempe(i,j-1))/hy(j-1)
                else
                    dTemdy(i,j) = (tempe(i,j)-tempe(i,j-1))/((hy(j)+hy(j-1))/2.D0)
                end if

            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols

                dvdt = (poro(i,j)*rhof*thetaf*tempe(i,j)+(1-poro(i,j))*rhos*thetas*tempe(i,j) - &!
                    poro_old(i,j)*rhof*thetaf*Tem(i,j)-(1-poro_old(i,j))*rhos*thetas*Tem(i,j))/(timeEnd/nt)

                div1 = rhof*thetaf * ((vx(i+1,j)*Tembarx(i+1,j) - vx(i,j)*Tembarx(i,j))/hx(i) + &!
                    (vy(i,j+1)*Tembary(i,j+1) - vy(i,j)*Tembary(i,j))/hy(j))

                div2 = (dTemdx(i+1,j)*(poroHarmX(i+1,j)*lamdaf+(1-poroHarmX(i+1,j))*lamdas) - &!
                    dTemdx(i,j)*(poroHarmX(i,j)*lamdaf+(1-poroHarmX(i,j))*lamdas))/hx(i) + &!
                    (dTemdy(i,j+1)*(poroHarmY(i,j+1)*lamdaf+(1-poroHarmY(i,j+1))*lamdas) - &!
                    dTemdy(i,j)*(poroHarmY(i,j)*lamdaf+(1-poroHarmY(i,j))*lamdas))/hy(j)

                div3 = p(i,j) * ((vx(i+1,j)-vx(i,j))/hx(i) + (vy(i,j+1)-vy(i,j))/hy(j))

                tp = visc * (((vx(i+1,j)/poroHarmX(i+1,j)-vx(i,j)/poroHarmX(i,j))/hx(i))**2.D0 + &!
                    ((vy(i,j+1)/poroHarmY(i,j+1)-vy(i,j)/poroHarmY(i,j))/hy(j))**2.D0)

                vmodulus = dsqrt(((vx(i,j)+vx(i+1,j))/2.D0)**2.D0 + ((vy(i,j)+vy(i,j+1))/2.D0)**2.D0)
                v1 = visc/Kxx(i,j)*vmodulus**2.D0

                Forchh = 1.75D0/dsqrt(1.5D2*poro(i,j)**3.D0)
                v2 = rhof*Forchh/dsqrt(Kxx(i,j))*vmodulus**3.D0

                HrT = abs(-9.702D3+1.697D1*Tem(i,j)-2.34D-3*Tem(i,j)**2.D0)
                arh = av(i,j) * kc(i,j)*ks(i,j)/(kc(i,j)+ks(i,j))*Cf(i,j) * HrT

                resi(i,j) = div2 - div3 + tp + v1 + v2 + arh - div1 - dvdt

            end do
        end do

        deallocate(Tembarx)
        deallocate(Tembary)
        deallocate(dTemdx)
        deallocate(dTemdy)

    end subroutine Resi_ener_tem

    ! End Calc Resi Temperature

end module DBF_resi
