
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

        integer, dimension(:,:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi
        
        integer :: j,k

        resi(:,:,:) = 0.D0

        ! Applying Neumann BC
        do k = 1, localnlays
            do j = 1, localnrows
                if((pcol==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                    resi(1,j,k) = velx(1,j,k) - vxBdryX0(ylower+j-1,zlower+k-1)
                end if
                if((pcol==pncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                    resi(localncols+1,j,k) = velx(localncols+1,j,k) - vxBdryX1(ylower+j-1,zlower+k-1)
                end if
            end do
        end do

    end subroutine Resi_xmom_vx_b

    subroutine Resi_xmom_vx(velx, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi

        real(kind=8) :: dpdx, vyAv, vzAv, vAbs, Forchh, ForchhTerm, vxD
        real(kind=8) :: DVxDx, DVxDy, DVxDz, velxLaplace
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: i, j, k

        indexl = 1
        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                        dpdx = 0.D0
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==1)) then
                        dpdx = (p(i,j,k)-pBdryX0(ylower+j-1,zlower+k-1)) / hx(i)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                        dpdx = 0.D0
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==1)) then
                        dpdx = (pBdryX1(ylower+j-1,zlower+k-1)-p(i-1,j,k)) / hx(i-1)
                    else
                        dpdx = (p(i,j,k)-p(i-1,j,k))/((hx(i)+hx(i-1))/2.D0)
                    end if

                    if((pcol==1).and.(i==1)) then
                        vyAv = 2.5D-1 * (vyBdryX0(ylower+j-1,zlower+k-1)+vyBdryX0(ylower+j,zlower+k-1)+vy(i,j,k)+vy(i,j+1,k))
                        vzAv = 2.5D-1 * (vzBdryX0(ylower+j-1,zlower+k-1)+vzBdryX0(ylower+j-1,zlower+k)+vz(i,j,k)+vz(i,j,k+1))
                    elseif((pcol==pncols).and.(i == localncols+1)) then
                        vyAv = 2.5D-1 * (vy(i-1,j,k)+vy(i-1,j+1,k)+vyBdryX1(ylower+j-1,zlower+k-1)+vyBdryX1(ylower+j,zlower+k-1))
                        vzAv = 2.5D-1 * (vz(i-1,j,k)+vz(i-1,j,k+1)+vzBdryX1(ylower+j-1,zlower+k-1)+vzBdryX1(ylower+j-1,zlower+k))
                    else
                        vyAv = 2.5D-1 * (vy(i-1,j,k)+vy(i-1,j+1,k)+vy(i,j,k)+vy(i,j+1,k))
                        vzAv = 2.5D-1 * (vz(i-1,j,k)+vz(i-1,j,k+1)+vz(i,j,k)+vz(i,j,k+1))
                    end if

                    vAbs = dsqrt(vx(i,j,k)**2.D0 + vyAv**2.D0 + vzAv**2.D0)

                    Forchh = 1.75D0/dsqrt(1.5D2*poroHarmX(i,j,k)**3.D0)
                    ForchhTerm = rhof*Forchh/dsqrt(KxxHarm(i,j,k))*vAbs
                    vxD = - (dpdx-rhof*gravX)/(ForchhTerm+visc/KxxHarm(i,j,k))

                    ! Determine the time term
                    resi(i,j,k) = - rhof*(velx(i,j,k)/poroHarmX(i,j,k)-vx(i,j,k)/poroHarmX_old(i,j,k))/(timeEnd/nt)

                    ! Determine the Forchhimer term
                    if(isForchheimer) then
                        if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0))).and. &!
                            (.not.((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)))) then
                            resi(i,j,k) = resi(i,j,k) - (visc/KxxHarm(i,j,k)+Forchh*rhof/dsqrt(KxxHarm(i,j,k))*vAbs)*(velx(i,j,k)-vxD)
                        end if
                    end if

                    ! Determine the advection term
                    if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                        DVxDx = (velx(i+1,j,k)/poroHarmX(i+1,j,k)-vxBdryX0(ylower+j-1,zlower+k-1)/poroHarmX(i,j,k))/hx(1)
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==1)) then
                        DVxDx = (velx(i+1,j,k)/poroHarmX(i+1,j,k)-velx(i,j,k)/poroHarmX(i,j,k))/hx(1)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                        DVxDx = (vxBdryX1(ylower+j-1,zlower+k-1)/poroHarmX(i,j,k)-velx(i-1,j,k)/poroHarmX(i-1,j,k))/hx(localncols)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==1)) then
                        DVxDx = (velx(i,j,k)/poroHarmX(i,j,k)-velx(i-1,j,k)/poroHarmX(i-1,j,k))/hx(localncols)
                    else
                        if(vx(i,j,k)>0.D0) then
                            DVxDx = (velx(i,j,k)/poroHarmX(i,j,k)-velx(i-1,j,k)/poroHarmX(i-1,j,k))/hx(i-1)
                        else
                            DVxDx = (velx(i+1,j,k)/poroHarmX(i+1,j,k)-velx(i,j,k)/poroHarmX(i,j,k))/hx(i)
                        end if
                    end if

                    if((prow==1).and.(j==1)) then
                        DVxDy = (velx(i,j,k)/poroHarmX(i,j,k)-vxBdryY0(xlower+i-1,zlower+k-1)/poroHarmX(i,j,k))/hy(j)
                    elseif((prow==pnrows).and.(j==localnrows)) then
                        DVxDy = (vxBdryY1(xlower+i-1,zlower+k-1)/poroHarmX(i,j,k)-velx(i,j,k)/poroHarmX(i,j,k))/hy(j)
                    else
                        if(vyAv>0.D0) then
                            DVxDy = (velx(i,j,k)/poroHarmX(i,j,k)-velx(i,j-1,k)/poroHarmX(i,j-1,k))/((hy(j-1)+hy(j))/2.D0)
                        else
                            DVxDy = (velx(i,j+1,k)/poroHarmX(i,j+1,k)-velx(i,j,k)/poroHarmX(i,j,k))/((hy(j)+hy(j+1))/2.D0)
                        end if
                    end if

                    if((play==1).and.(k==1)) then
                        DVxDz = (velx(i,j,k)/poroHarmX(i,j,k)-vxBdryZ0(xlower+i-1,ylower+j-1)/poroHarmX(i,j,k))/hz(k)
                    elseif((play==pnlays).and.(k==localnlays)) then
                        DVxDz = (vxBdryZ1(xlower+i-1,ylower+j-1)/poroHarmX(i,j,k)-velx(i,j,k)/poroHarmX(i,j,k))/hz(k)
                    else
                        if(vzAv>0.D0) then
                            DVxDz = (velx(i,j,k)/poroHarmX(i,j,k)-velx(i,j,k-1)/poroHarmX(i,j,k-1))/((hz(k-1)+hz(k))/2.D0)
                        else
                            DVxDz = (velx(i,j,k+1)/poroHarmX(i,j,k+1)-velx(i,j,k)/poroHarmX(i,j,k))/((hz(k)+hz(k+1))/2.D0)
                        end if
                    end if

                    resi(i,j,k) = resi(i,j,k) - rhof/poroHarmX_old(i,j,k) * (vx(i,j,k)*DVxDx+vyAv*DVxDy+vzAv*DVxDz)

                    ! Determine the Brinkman term
                    if(isBrinkman) then
                        if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                            velxLaplace = (velx(i+1,j,k)/poroHarmX(i+1,j,k)-vxBdryX0(ylower+j-1,zlower+k-1)/poroHarmX(i,j,k))/(hx(1)**2.D0)
                        elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==1)) then
                            velxLaplace = (velx(i+1,j,k)/poroHarmX(i+1,j,k)-velx(i,j,k)/poroHarmX(i,j,k))/(hx(1)**2.D0)
                        elseif((pcol==pncols).and.(i == localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                            velxLaplace = (vxBdryX1(ylower+j-1,zlower+k-1)/poroHarmX(i,j,k)-velx(i-1,j,k)/poroHarmX(i-1,j,k))/(hx(localncols)**2.D0)
                        elseif((pcol==pncols).and.(i == localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==1)) then
                            velxLaplace = (velx(i,j,k)/poroHarmX(i,j,k)-velx(i-1,j,k)/poroHarmX(i-1,j,k))/(hx(localncols)**2.D0)
                        else
                            velxLaplace = ((velx(i+1,j,k)/poroHarmX(i+1,j,k)-velx(i,j,k)/poroHarmX(i,j,k))/hx(i) - &!
                                (velx(i,j,k)/poroHarmX(i,j,k)-velx(i-1,j,k)/poroHarmX(i-1,j,k))/hx(i-1)) / ((hx(i-1)+hx(i))/2.D0)
                        end if

                        if((prow==1).and.(j == 1)) then
                            velxLaplace = velxLaplace + ((velx(i,j+1,k)/poroHarmX(i,j+1,k)-velx(i,j,k)/poroHarmX(i,j,k))/ &!
                                ((hy(j+1)+hy(j))/2.D0) - (velx(i,j,k)/poroHarmX(i,j,k)-vxBdryY0(xlower+i-1,zlower+k-1)/ &!
                                poroHarmX(i,j,k))/hy(1)) / hy(j)
                        elseif((prow==pnrows).and.(j == localnrows)) then
                            velxLaplace = velxLaplace + ((vxBdryY1(xlower+i-1,zlower+k-1)/poroHarmX(i,j,k)-velx(i,j,k)/ &!
                                poroHarmX(i,j,k))/hy(localnrows) - (velx(i,j,k)/poroHarmX(i,j,k)-velx(i,j-1,k)/ &!
                                poroHarmX(i,j-1,k))/((hy(j)+hy(j-1))/2.D0)) / hy(j)
                        else
                            velxLaplace = velxLaplace + ((velx(i,j+1,k)/poroHarmX(i,j+1,k)-velx(i,j,k)/poroHarmX(i,j,k))/ &!
                                ((hy(j+1)+hy(j))/2.D0) - (velx(i,j,k)/poroHarmX(i,j,k)-velx(i,j-1,k)/poroHarmX(i,j-1,k))/ &!
                                ((hy(j)+hy(j-1))/2.D0)) / hy(j)
                        end if

                        if((play==1).and.(k == 1)) then
                            velxLaplace = velxLaplace + ((velx(i,j,k+1)/poroHarmX(i,j,k+1)-velx(i,j,k)/poroHarmX(i,j,k))/ &!
                                ((hz(k+1)+hz(k))/2.D0) - (velx(i,j,k)/poroHarmX(i,j,k)-vxBdryZ0(xlower+i-1,ylower+j-1)/ &!
                                poroHarmX(i,j,k))/hz(1)) / hz(k)
                        elseif((play==pnlays).and.(k == localnlays)) then
                            velxLaplace = velxLaplace + ((vxBdryZ1(xlower+i-1,ylower+j-1)/poroHarmX(i,j,k)-velx(i,j,k)/ &!
                                poroHarmX(i,j,k))/hz(localnlays) - (velx(i,j,k)/poroHarmX(i,j,k)-velx(i,j,k-1)/ &!
                                poroHarmX(i,j,k-1))/((hz(k)+hz(k-1))/2.D0)) / hz(k)
                        else
                            velxLaplace = velxLaplace + ((velx(i,j,k+1)/poroHarmX(i,j,k+1)-velx(i,j,k)/ &!
                                poroHarmX(i,j,k))/((hz(k+1)+hz(k))/2.D0) - (velx(i,j,k)/poroHarmX(i,j,k)-velx(i,j,k-1)/ &!
                                poroHarmX(i,j,k-1))/((hz(k)+hz(k-1))/2.D0)) / hz(k)
                        end if

                        resi(i,j,k) = resi(i,j,k) + visc*velxLaplace

                    end if
                    
                end do
            end do
        end do

    end subroutine Resi_xmom_vx

    ! End Calc Resi x-velocity

    ! Calc Resi y-velocity
    ! The derivatives are evaluated on edges

    subroutine Resi_ymom_vy_b(vely, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi

        integer :: i, k

        resi(:,:,:) = 0.D0

        ! Applying Neumann BC
        do k = 1, localnlays
            do i = 1, localncols
                if((prow==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                    resi(i,1,k) = vely(i,1,k) - vyBdryY0(xlower+i-1,zlower+k-1)
                end if
                if((prow==pnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                    resi(i,localnrows+1,k) = vely(i,localnrows+1,k) - vyBdryY1(xlower+i-1,zlower+k-1)
                end if
            end do
        end do
        
    end subroutine Resi_ymom_vy_b

    subroutine Resi_ymom_vy(vely, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi
       
        real(kind=8) :: dpdy, vxAv, vzAv, vAbs, Forchh, ForchhTerm, vyD
        real(kind=8) :: DVyDx, DVyDy, DVyDz, velyLaplace
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: i, j, k
        
        indexl = 1
        indexr = localncols
        indexd = 1
        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if
        indexf = 1
        indexb = localnlays

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                        dpdy = 0.D0
                    elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==1)) then
                        dpdy = (p(i,j,k)-pBdryY0(xlower+i-1,zlower+k-1)) / hy(j)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                        dpdy = 0.D0
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==1)) then
                        dpdy = (pBdryY1(xlower+i-1,zlower+k-1)-p(i,j-1,k)) / hy(j-1)
                    else
                        dpdy = (p(i,j,k)-p(i,j-1,k))/((hy(j)+hy(j-1))/2.D0)
                    end if

                    ! Determine the advection term
                    if((prow==1).and.(j==1)) then
                        vxAv = 2.5D-1 * (vxBdryY0(xlower+i-1,zlower+k-1)+vxBdryY0(xlower+i,zlower+k-1)+vx(i,j,k)+vx(i+1,j,k))
                        vzAv = 2.5D-1 * (vzBdryY0(xlower+i-1,zlower+k-1)+vzBdryY0(xlower+i-1,zlower+k)+vz(i,j,k)+vz(i,j,k+1))
                    elseif((prow==pnrows).and.(j==localnrows+1)) then
                        vxAv = 2.5D-1 * (vx(i,j-1,k)+vx(i+1,j-1,k)+vxBdryY1(xlower+i-1,zlower+k-1)+vxBdryY1(xlower+i,zlower+k-1))
                        vzAv = 2.5D-1 * (vz(i,j-1,k)+vz(i,j-1,k+1)+vzBdryY1(xlower+i-1,zlower+k-1)+vzBdryY1(xlower+i-1,zlower+k))
                    else
                        vxAv = 2.5D-1 * (vx(i,j-1,k)+vx(i+1,j-1,k)+vx(i,j,k)+vx(i+1,j,k))
                        vzAv = 2.5D-1 * (vz(i,j-1,k)+vz(i,j-1,k+1)+vz(i,j,k)+vz(i,j,k+1))
                    end if

                    vAbs = dsqrt(vxAv**2.D0 + vy(i,j,k)**2.D0 + vzAv**2.D0)

                    Forchh = 1.75D0/dsqrt(1.5D2*poroHarmY(i,j,k)**3.D0)
                    ForchhTerm = rhof*Forchh/dsqrt(KyyHarm(i,j,k))*vAbs
                    vyD = - (dpdy-rhof*gravY)/(ForchhTerm+visc/KyyHarm(i,j,k))

                    ! Determine the time term
                    resi(i,j,k) = - rhof*(vely(i,j,k)/poroHarmY(i,j,k)-vy(i,j,k)/poroHarmY_old(i,j,k))/(timeEnd/nt)

                    ! Determine the Forchheimer term
                    if(isForchheimer) then
                        if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0))).and. &!
                            (.not.((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)))) then
                            resi(i,j,k) = resi(i,j,k) - (visc/KyyHarm(i,j,k)+Forchh*rhof/dsqrt(KyyHarm(i,j,k))*vAbs)*(vely(i,j,k)-vyD)
                        end if
                    end if

                    ! Determine the advection term
                    if((pcol==1).and.(i==1)) then
                        DVyDx = (vely(i,j,k)/poroHarmY(i,j,k)-vyBdryX0(ylower+j-1,zlower+k-1)/poroHarmY(i,j,k))/hx(i)
                    elseif((pcol==pncols).and.(i==localncols)) then
                        DVyDx = (vyBdryX1(ylower+j-1,zlower+k-1)/poroHarmY(i,j,k)-vely(i,j,k)/poroHarmY(i,j,k))/hx(i)
                    else
                        if(vxAv>0.D0) then
                            DVyDx = (vely(i,j,k)/poroHarmY(i,j,k)-vely(i-1,j,k)/poroHarmY(i-1,j,k))/((hx(i-1)+hx(i))/2.D0)
                        else
                            DVyDx = (vely(i+1,j,k)/poroHarmY(i+1,j,k)-vely(i,j,k)/poroHarmY(i,j,k))/((hx(i)+hx(i+1))/2.D0)
                        end if
                    end if
                    if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                        DVyDy = (vely(i,j+1,k)/poroHarmY(i,j+1,k)-vyBdryY0(xlower+i-1,zlower+k-1)/poroHarmY(i,j,k))/hy(j)
                    elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==1)) then
                        DVyDy = (vely(i,j+1,k)/poroHarmY(i,j+1,k)-vely(i,j,k)/poroHarmY(i,j,k))/hy(j)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                        DVyDy = (vyBdryY1(xlower+i-1,zlower+k-1)/poroHarmY(i,j,k)-vely(i,j-1,k)/poroHarmY(i,j-1,k))/hy(j-1)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==1)) then
                        DVyDy = (vely(i,j,k)/poroHarmY(i,j,k)-vely(i,j-1,k)/poroHarmY(i,j-1,k))/hy(j-1)
                    else
                        if(vy(i,j,k)>0.D0) then
                            DVyDy = (vely(i,j,k)/poroHarmY(i,j,k)-vely(i,j-1,k)/poroHarmY(i,j-1,k))/hy(j-1)
                        else
                            DVyDy = (vely(i,j+1,k)/poroHarmY(i,j+1,k)-vely(i,j,k)/poroHarmY(i,j,k))/hy(j)
                        end if
                    end if
                    if((play==1).and.(k==1)) then
                        DVyDz = (vely(i,j,k)/poroHarmY(i,j,k)-vyBdryZ0(xlower+i-1,ylower+j-1)/poroHarmY(i,j,k))/hz(k)
                    elseif((play==pnlays).and.(k==localnlays)) then
                        DVyDz = (vyBdryZ1(xlower+i-1,ylower+j-1)/poroHarmY(i,j,k)-vely(i,j,k)/poroHarmY(i,j,k))/hz(k)
                    else
                        if(vzAv>0.D0) then
                            DVyDz = (vely(i,j,k)/poroHarmY(i,j,k)-vely(i,j,k-1)/poroHarmY(i,j,k-1))/((hz(k-1)+hz(k))/2.D0)
                        else
                            DVyDz = (vely(i,j,k+1)/poroHarmY(i,j,k+1)-vely(i,j,k)/poroHarmY(i,j,k))/((hz(k)+hz(k+1))/2.D0)
                        end if
                    end if

                    resi(i,j,k) = resi(i,j,k) - rhof/poroHarmY_old(i,j,k) * (vxAv*DVyDx+vy(i,j,k)*DVyDy+vzAv*DVyDz)

                    ! Determine the Brinkman term
                    if(isBrinkman) then
                        if((pcol==1).and.(i==1)) then
                            velyLaplace = ((vely(i+1,j,k)/poroHarmY(i+1,j,k)-vely(i,j,k)/poroHarmY(i,j,k))/((hx(i+1)+hx(i))/2.D0) &!
                                - (vely(i,j,k)/poroHarmY(i,j,k)-vyBdryX0(ylower+j-1,zlower+k-1)/poroHarmY(i,j,k))/hx(i)) / hx(i)
                        elseif((pcol==pncols).and.(i==localncols)) then
                            velyLaplace = ((vyBdryX1(ylower+j-1,zlower+k-1)/poroHarmY(i,j,k)-vely(i,j,k)/poroHarmY(i,j,k))/hx(i) - &!
                                (vely(i,j,k)/poroHarmY(i,j,k)-vely(i-1,j,k)/poroHarmY(i-1,j,k))/((hx(i)+hx(i-1))/2.D0)) / hx(i)
                        else
                            velyLaplace = ((vely(i+1,j,k)/poroHarmY(i+1,j,k)-vely(i,j,k)/poroHarmY(i,j,k))/ &!
                                ((hx(i+1)+hx(i))/2.D0) - (vely(i,j,k)/poroHarmY(i,j,k)-vely(i-1,j,k)/ &!
                                poroHarmY(i-1,j,k))/((hx(i)+hx(i-1))/2.D0)) / hx(i)
                        end if
                        if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                            velyLaplace = velyLaplace + (vely(i,j+1,k)/poroHarmY(i,j+1,k)-vyBdryY0(xlower+i-1,zlower+k-1)/ &!
                                poroHarmY(i,j,k))/(hy(j)**2.D0)
                        elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==1)) then
                            velyLaplace = velyLaplace + (vely(i,j+1,k)/poroHarmY(i,j+1,k)-vely(i,j,k)/ &!
                                poroHarmY(i,j,k))/(hy(j)**2.D0)
                        elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                            velyLaplace = velyLaplace + (vyBdryY1(xlower+i-1,zlower+k-1)/poroHarmY(i,j,k)-vely(i,j-1,k)/ &!
                                poroHarmY(i,j-1,k))/(hy(j-1)**2.D0)
                        elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==1)) then
                            velyLaplace = velyLaplace + (vely(i,j,k)/poroHarmY(i,j,k)-vely(i,j-1,k)/ &!
                                poroHarmY(i,j-1,k))/(hy(j-1)**2.D0)
                        else
                            velyLaplace = velyLaplace + ((vely(i,j+1,k)/poroHarmY(i,j+1,k)-vely(i,j,k)/ &!
                                poroHarmY(i,j,k))/hy(j) - (vely(i,j,k)/poroHarmY(i,j,k)-vely(i,j-1,k)/ &!
                                poroHarmY(i,j-1,k))/hy(j-1)) / ((hy(j-1)+hy(j))/2.D0)
                        end if
                        if((play==1).and.(k==1)) then
                            velyLaplace = velyLaplace + ((vely(i,j,k+1)/poroHarmY(i,j,k+1)-vely(i,j,k)/ &!
                                poroHarmY(i,j,k))/((hz(k+1)+hz(k))/2.D0) - (vely(i,j,k)/poroHarmY(i,j,k)- &!
                                vyBdryZ0(xlower+i-1,ylower+j-1)/poroHarmY(i,j,k))/hz(k)) / hz(k)
                        elseif((play==pnlays).and.(k==localnlays)) then
                            velyLaplace = velyLaplace + ((vyBdryZ1(xlower+i-1,ylower+j-1)/poroHarmY(i,j,k)-vely(i,j,k)/ &!
                                poroHarmY(i,j,k))/hz(k) - (vely(i,j,k)/poroHarmY(i,j,k)-vely(i,j,k-1)/ &!
                                poroHarmY(i,j,k-1))/((hz(k)+hz(k-1))/2.D0)) / hz(k)
                        else
                            velyLaplace = velyLaplace + ((vely(i,j,k+1)/poroHarmY(i,j,k+1)-vely(i,j,k)/ &!
                                poroHarmY(i,j,k))/((hz(k+1)+hz(k))/2.D0) - (vely(i,j,k)/poroHarmY(i,j,k)- &!
                                vely(i,j,k-1)/poroHarmY(i,j,k-1))/((hz(k)+hz(k-1))/2.D0)) / hz(k)
                        end if

                        resi(i,j,k) = resi(i,j,k) + visc*velyLaplace

                    end if

                end do
            end do
        end do

    end subroutine Resi_ymom_vy

    ! End Calc Resi y-velocity

    ! Calc Resi z-velocity
    ! The derivatives are evaluated on edges

    subroutine Resi_zmom_vz_b(velz, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: velz
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi

        integer :: i, j

        resi(:,:,:) = 0.D0

        ! Applying Neumann BC
        do j = 1, localnrows
            do i = 1, localncols
                if((play==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                    resi(i,j,1) = velz(i,j,1) - vzBdryZ0(xlower+i-1,ylower+j-1)
                end if
                if((play==pnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                    resi(i,j,localnlays+1) = velz(i,j,localnlays+1) - vzBdryZ1(xlower+i-1,ylower+j-1)
                end if
            end do
        end do

    end subroutine Resi_zmom_vz_b

    subroutine Resi_zmom_vz(velz, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: velz
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi
        
        real(kind=8) :: dpdz, vxAv, vyAv, vAbs, Forchh, ForchhTerm, vzD
        real(kind=8) :: DVzDx, DVzDy, DVzDz, velzLaplace
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: i, j, k

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        if(play /= pnlays) then
            indexb = localnlays
        else
            indexb = localnlays + 1
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                        dpdz = 0.D0
                    elseif((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==1)) then
                        dpdz = (p(i,j,k)-pBdryZ0(xlower+i-1,ylower+j-1)) / hz(k)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                        dpdz = 0.D0
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==1)) then
                        dpdz = (pBdryZ1(xlower+i-1,ylower+j-1)-p(i,j,k-1)) / hz(k-1)
                    else
                        dpdz = (p(i,j,k)-p(i,j,k-1))/((hz(k)+hz(k-1))/2.D0)
                    end if

                    if((play==1).and.(k==1)) then
                        vxAv = 2.5D-1 * (vxBdryZ0(xlower+i-1,ylower+j-1)+vxBdryZ0(xlower+i,ylower+j-1)+vx(i,j,k)+vx(i+1,j,k))
                        vyAv = 2.5D-1 * (vyBdryZ0(xlower+i-1,ylower+j-1)+vyBdryZ0(xlower+i-1,ylower+j)+vy(i,j,k)+vy(i,j+1,k))
                    elseif((play==pnlays).and.(k==localnlays+1)) then
                        vxAv = 2.5D-1 * (vx(i,j,k-1)+vx(i+1,j,k-1)+vxBdryZ1(xlower+i-1,ylower+j-1)+vxBdryZ1(xlower+i,ylower+j-1))
                        vyAv = 2.5D-1 * (vy(i,j,k-1)+vy(i,j+1,k-1)+vyBdryZ1(xlower+i-1,ylower+j-1)+vyBdryZ1(xlower+i-1,ylower+j))
                    else
                        vxAv = 2.5D-1 * (vx(i,j,k-1)+vx(i+1,j,k-1)+vx(i,j,k)+vx(i+1,j,k))
                        vyAv = 2.5D-1 * (vy(i,j,k-1)+vy(i,j+1,k-1)+vy(i,j,k)+vy(i,j+1,k))
                    end if

                    vAbs = dsqrt(vxAv**2.D0 + vyAv**2.D0 + vz(i,j,k)**2.D0)

                    Forchh = 1.75D0/dsqrt(1.5D2*poroHarmZ(i,j,k)**3.D0)
                    ForchhTerm = rhof*Forchh/dsqrt(KzzHarm(i,j,k))*vAbs
                    vzD = - (dpdz-rhof*gravZ)/(ForchhTerm+visc/KzzHarm(i,j,k))

                    ! Determine the time term
                    resi(i,j,k) = - rhof*(velz(i,j,k)/poroHarmZ(i,j,k)-vz(i,j,k)/poroHarmZ_old(i,j,k))/(timeEnd/nt)

                    ! Determine the Forchheimer term
                    if(isForchheimer) then
                        if((.not.((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0))).and. &!
                            (.not.((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)))) then
                            resi(i,j,k) = resi(i,j,k) - (visc/KzzHarm(i,j,k)+Forchh*rhof/dsqrt(KzzHarm(i,j,k))*vAbs)*(velz(i,j,k)-vzD)
                        end if
                    end if

                    ! Determine the advection term
                    if((pcol==1).and.(i==1)) then
                        DVzDx = (velz(i,j,k)/poroHarmZ(i,j,k)-vzBdryX0(ylower+j-1,zlower+k-1)/poroHarmZ(i,j,k))/hx(i)
                    elseif((pcol==pncols).and.(i==localncols)) then
                        DVzDx = (vzBdryX1(ylower+j-1,zlower+k-1)/poroHarmZ(i,j,k)-velz(i,j,k)/poroHarmZ(i,j,k))/hx(i)
                    else
                        if(vxAv>0.D0) then
                            DVzDx = (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i-1,j,k)/poroHarmZ(i-1,j,k))/((hx(i-1)+hx(i))/2.D0)
                        else
                            DVzDx = (velz(i+1,j,k)/poroHarmZ(i+1,j,k)-velz(i,j,k)/poroHarmZ(i,j,k))/((hx(i)+hx(i+1))/2.D0)
                        end if
                    end if
                    if((prow==1).and.(j==1)) then
                        DVzDy = (velz(i,j,k)/poroHarmZ(i,j,k)-vzBdryY0(xlower+i-1,zlower+k-1)/poroHarmZ(i,j,k))/hy(j)
                    elseif((prow==pnrows).and.(j==localnrows)) then
                        DVzDy = (vzBdryY1(xlower+i-1,zlower+k-1)/poroHarmZ(i,j,k)-velz(i,j,k)/poroHarmZ(i,j,k))/hy(j)
                    else
                        if(vyAv>0.D0) then
                            DVzDy = (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i,j-1,k)/poroHarmZ(i,j-1,k))/((hy(j-1)+hy(j))/2.D0)
                        else
                            DVzDy = (velz(i,j+1,k)/poroHarmZ(i,j+1,k)-velz(i,j,k)/poroHarmZ(i,j,k))/((hy(j)+hy(j+1))/2.D0)
                        end if
                    end if
                    if((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                        DVzDz = (velz(i,j,k+1)/poroHarmZ(i,j,k+1)-vzBdryZ0(xlower+i-1,ylower+j-1)/poroHarmZ(i,j,k))/hz(k)
                    elseif((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==1)) then
                        DVzDz = (velz(i,j,k+1)/poroHarmZ(i,j,k+1)-velz(i,j,k)/poroHarmZ(i,j,k))/hz(k)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                        DVzDz = (vzBdryZ1(xlower+i-1,ylower+j-1)/poroHarmZ(i,j,k)-velz(i,j,k-1)/poroHarmZ(i,j,k-1))/hz(k-1)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==1)) then
                        DVzDz = (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i,j,k-1)/poroHarmZ(i,j,k-1))/hz(k-1)
                    else
                        if(vz(i,j,k)>0.D0) then
                            DVzDz = (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i,j,k-1)/poroHarmZ(i,j,k-1))/hz(k-1)
                        else
                            DVzDz = (velz(i,j,k+1)/poroHarmZ(i,j,k+1)-velz(i,j,k)/poroHarmZ(i,j,k))/hz(k)
                        end if
                    end if

                    resi(i,j,k) = resi(i,j,k) - rhof/poroHarmZ_old(i,j,k) * (vxAv*DVzDx+vyAv*DVzDy+vz(i,j,k)*DVzDz)

                    ! Determine the Brinkman term
                    if(isBrinkman) then
                        if((pcol==1).and.(i==1)) then
                            velzLaplace = ((velz(i+1,j,k)/poroHarmZ(i+1,j,k)-velz(i,j,k)/poroHarmZ(i,j,k))/ &!
                                ((hx(i+1)+hx(i))/2.D0) - (velz(i,j,k)/poroHarmZ(i,j,k)-vzBdryX0(ylower+j-1,zlower+k-1)/ &!
                                poroHarmZ(i,j,k))/hx(i)) / hx(i)
                        elseif((pcol==pncols).and.(i==localncols)) then
                            velzLaplace = ((vzBdryX1(ylower+j-1,zlower+k-1)/poroHarmZ(i,j,k)-velz(i,j,k)/ &!
                                poroHarmZ(i,j,k))/hx(i) - (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i-1,j,k)/ &!
                                poroHarmZ(i-1,j,k))/((hx(i)+hx(i-1))/2.D0)) / hx(i)
                        else
                            velzLaplace = ((velz(i+1,j,k)/poroHarmZ(i+1,j,k)-velz(i,j,k)/poroHarmZ(i,j,k))/ &!
                                ((hx(i+1)+hx(i))/2.D0) - (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i-1,j,k)/ &!
                                poroHarmZ(i-1,j,k))/((hx(i)+hx(i-1))/2.D0)) / hx(i)
                        end if
                        if((prow==1).and.(j==1)) then
                            velzLaplace = velzLaplace + ((velz(i,j+1,k)/poroHarmZ(i,j+1,k)-velz(i,j,k)/ &!
                                poroHarmZ(i,j,k))/((hy(j+1)+hy(j))/2.D0) - (velz(i,j,k)/ &!
                                poroHarmZ(i,j,k)-vzBdryY0(xlower+i-1,zlower+k-1)/poroHarmZ(i,j,k))/hy(j)) / hy(j)
                        elseif((prow==pnrows).and.(j==localnrows)) then
                            velzLaplace = velzLaplace + ((vzBdryY1(xlower+i-1,zlower+k-1)/poroHarmZ(i,j,k)- &!
                                velz(i,j,k)/poroHarmZ(i,j,k))/hy(j) - (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i,j-1,k)/ &!
                                poroHarmZ(i,j-1,k))/((hy(j)+hy(j-1))/2.D0)) / hy(j)
                        else
                            velzLaplace = velzLaplace + ((velz(i,j+1,k)/poroHarmZ(i,j+1,k)-velz(i,j,k)/ &!
                                poroHarmZ(i,j,k))/((hy(j+1)+hy(j))/2.D0) - (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i,j-1,k)/ &!
                                poroHarmZ(i,j-1,k))/((hy(j)+hy(j-1))/2.D0)) / hy(j)
                        end if
                        if((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                            velzLaplace = velzLaplace + (velz(i,j,k+1)/poroHarmZ(i,j,k+1)-vzBdryZ0(xlower+i-1,ylower+j-1)/ &!
                                poroHarmZ(i,j,k))/(hz(k)**2.D0)
                        elseif((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==1)) then
                            velzLaplace = velzLaplace + (velz(i,j,k+1)/poroHarmZ(i,j,k+1)-velz(i,j,k)/ &!
                                poroHarmZ(i,j,k))/(hz(k)**2.D0)
                        elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                            velzLaplace = velzLaplace + (vzBdryZ1(xlower+i-1,ylower+j-1)/poroHarmZ(i,j,k)-velz(i,j,k-1)/ &!
                                poroHarmZ(i,j,k-1))/(hz(k-1)**2.D0)
                        elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==1)) then
                            velzLaplace = velzLaplace + (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i,j,k-1)/ &!
                                poroHarmZ(i,j,k-1))/(hz(k-1)**2.D0)
                        else
                            velzLaplace = velzLaplace + ((velz(i,j,k+1)/poroHarmZ(i,j,k+1)-velz(i,j,k)/ &!
                                poroHarmZ(i,j,k))/hz(k) - (velz(i,j,k)/poroHarmZ(i,j,k)-velz(i,j,k-1)/ &!
                                poroHarmZ(i,j,k-1))/hz(k-1)) / ((hz(k-1)+hz(k))/2.D0)
                        end if
                        
                        resi(i,j,k) = resi(i,j,k) + visc*velzLaplace
                    end if
                   
                end do
            end do
        end do

    end subroutine Resi_zmom_vz

    ! End Calc Resi z-velocity

    ! Calc Resi Mass

    subroutine Resi_mass_p(pres, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi

        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        real(kind=8), dimension(:,:,:), allocatable :: dpdx, dpdy, dpdz
        real(kind=8) :: vyAv_l, vzAv_l, vyAv_r, vzAv_r, vxAv_d, vzAv_d, vxAv_u, &!
            vzAv_u, vxAv_f, vyAv_f, vxAv_b, vyAv_b
        real(kind=8) :: vAbs_l, vAbs_r, vAbs_d, vAbs_u, vAbs_f, vAbs_b
        real(kind=8) :: div_l, div_r, div_d, div_u, div_f, div_b
        real(kind=8) :: ForchhTerm_l, ForchhTerm_r, ForchhTerm_d, ForchhTerm_u, &!
            ForchhTerm_f, ForchhTerm_b, div
        integer :: i, j, k

        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays
        allocate(dpdx(1:localncols+1,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = 1, localncols+1

                    if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                        dpdx(i,j,k) = 0.D0
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==1)) then
                        dpdx(i,j,k) = (pres(i,j,k)-pBdryX0(ylower+j-1,zlower+k-1))/hx(i)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                        dpdx(i,j,k) = 0.D0
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==1)) then
                        dpdx(i,j,k) = (pBdryX1(ylower+j-1,zlower+k-1)-pres(i-1,j,k))/hx(i-1)
                    else
                        dpdx(i,j,k) = (pres(i,j,k)-pres(i-1,j,k))/((hx(i)+hx(i-1))/2.D0)
                    end if

                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexf = 1
        indexb = localnlays
        allocate(dpdy(indexl:indexr,1:localnrows+1,indexf:indexb))

        do k = indexf, indexb
            do j = 1, localnrows+1
                do i = indexl, indexr

                    if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                        dpdy(i,j,k) = 0.D0
                    elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==1)) then
                        dpdy(i,j,k) = (pres(i,j,k)-pBdryY0(xlower+i-1,zlower+k-1))/hy(j)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                        dpdy(i,j,k) = 0.D0
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==1)) then
                        dpdy(i,j,k) = (pBdryY1(xlower+i-1,zlower+k-1)-pres(i,j-1,k))/hy(j-1)
                    else
                        dpdy(i,j,k) = (pres(i,j,k)-pres(i,j-1,k))/((hy(j)+hy(j-1))/2.D0)
                    end if

                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        allocate(dpdz(indexl:indexr,indexd:indexu,1:localnlays+1))

        do k = 1, localnlays+1
            do j = indexd, indexu
                do i = indexl, indexr

                    if((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                        dpdz(i,j,k) = 0.D0
                    elseif((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==1)) then
                        dpdz(i,j,k) = (pres(i,j,k)-pBdryZ0(xlower+i-1,ylower+j-1))/hz(k)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                        dpdz(i,j,k) = 0.D0
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==1)) then
                        dpdz(i,j,k) = (pBdryZ1(xlower+i-1,ylower+j-1)-pres(i,j,k-1))/hz(k-1)
                    else
                        dpdz(i,j,k) = (pres(i,j,k)-pres(i,j,k-1))/((hz(k)+hz(k-1))/2.D0)
                    end if

                end do
            end do
        end do

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    if((pcol==1).and.(i==1)) then
                        vyAv_l = (vyBdryX0(ylower+j-1,zlower+k-1)+vyBdryX0(ylower+j,zlower+k-1)+vy(i,j,k)+vy(i,j+1,k)) / 4.D0
                        vzAv_l = (vzBdryX0(ylower+j-1,zlower+k-1)+vzBdryX0(ylower+j-1,zlower+k)+vz(i,j,k)+vz(i,j,k+1)) / 4.D0
                    else
                        vyAv_l = (vy(i-1,j,k)+vy(i-1,j+1,k)+vy(i,j,k)+vy(i,j+1,k)) / 4.D0
                        vzAv_l = (vz(i-1,j,k)+vz(i-1,j,k+1)+vz(i,j,k)+vz(i,j,k+1)) / 4.D0
                    end if
                    if((pcol==pncols).and.(i==localncols)) then
                        vyAv_r = (vy(i,j,k)+vy(i,j+1,k)+vyBdryX1(ylower+j-1,zlower+k-1)+vyBdryX1(ylower+j,zlower+k-1)) / 4.D0
                        vzAv_r = (vz(i,j,k)+vz(i,j,k+1)+vzBdryX1(ylower+j-1,zlower+k-1)+vzBdryX1(ylower+j-1,zlower+k)) / 4.D0
                    else
                        vyAv_r = (vy(i+1,j,k)+vy(i+1,j+1,k)+vy(i,j,k)+vy(i,j+1,k)) / 4.D0
                        vzAv_r = (vz(i+1,j,k)+vz(i+1,j,k+1)+vz(i,j,k)+vz(i,j,k+1)) / 4.D0
                    end if

                    if((prow==1).and.(j==1)) then
                        vxAv_d = (vxBdryY0(xlower+i-1,zlower+k-1)+vxBdryY0(xlower+i,zlower+k-1)+vx(i,j,k)+vx(i+1,j,k)) / 4.D0
                        vzAv_d = (vzBdryY0(xlower+i-1,zlower+k-1)+vzBdryY0(xlower+i-1,zlower+k)+vz(i,j,k)+vz(i,j,k+1)) / 4.D0
                    else
                        vxAv_d = (vx(i,j-1,k)+vx(i+1,j-1,k)+vx(i,j,k)+vx(i+1,j,k)) / 4.D0
                        vzAv_d = (vz(i,j-1,k)+vz(i,j-1,k+1)+vz(i,j,k)+vz(i,j,k+1)) / 4.D0
                    end if
                    if((prow==pnrows).and.(j==localnrows)) then
                        vxAv_u = (vx(i,j,k)+vx(i+1,j,k)+vxBdryY1(xlower+i-1,zlower+k-1)+vxBdryY1(xlower+i,zlower+k-1)) / 4.D0
                        vzAv_u = (vz(i,j,k)+vz(i,j,k+1)+vzBdryY1(xlower+i-1,zlower+k-1)+vzBdryY1(xlower+i-1,zlower+k)) / 4.D0
                    else
                        vxAv_u = (vx(i,j+1,k)+vx(i+1,j+1,k)+vx(i,j,k)+vx(i+1,j,k)) / 4.D0
                        vzAv_u = (vz(i,j+1,k)+vz(i,j+1,k+1)+vz(i,j,k)+vz(i,j,k+1)) / 4.D0
                    end if

                    if((play==1).and.(k==1)) then
                        vxAv_f = (vxBdryZ0(xlower+i-1,ylower+j-1)+vxBdryZ0(xlower+i,ylower+j-1)+vx(i,j,k)+vx(i+1,j,k)) / 4.D0
                        vyAv_f = (vyBdryZ0(xlower+i-1,ylower+j-1)+vyBdryZ0(xlower+i-1,ylower+j)+vy(i,j,k)+vy(i,j+1,k)) / 4.D0
                    else
                        vxAv_f = (vx(i,j,k-1)+vx(i+1,j,k-1)+vx(i,j,k)+vx(i+1,j,k)) / 4.D0
                        vyAv_f = (vy(i,j,k-1)+vy(i,j+1,k-1)+vy(i,j,k)+vy(i,j+1,k)) / 4.D0
                    end if
                    if((play==pnlays).and.(k==localnlays)) then
                        vxAv_b = (vxBdryZ1(xlower+i-1,ylower+j-1)+vxBdryZ1(xlower+i,ylower+j-1)+vx(i,j,k)+vx(i+1,j,k)) / 4.D0
                        vyAv_b = (vyBdryZ1(xlower+i-1,ylower+j-1)+vyBdryZ1(xlower+i-1,ylower+j)+vy(i,j,k)+vy(i,j+1,k)) / 4.D0
                    else
                        vxAv_b = (vx(i,j,k)+vx(i+1,j,k)+vx(i,j,k+1)+vx(i+1,j,k+1)) / 4.D0
                        vyAv_b = (vy(i,j,k)+vy(i,j+1,k)+vy(i,j,k+1)+vy(i,j+1,k+1)) / 4.D0
                    end if

                    vAbs_l = dsqrt(vx(i,j,k)**2.D0 + vyAv_l**2.D0 + vzAv_l**2.D0)
                    vAbs_r = dsqrt(vx(i+1,j,k)**2.D0 + vyAv_r**2.D0 + vzAv_r**2.D0)
                    vAbs_d = dsqrt(vxAv_d**2.D0 + vy(i,j,k)**2.D0 + vzAv_d**2.D0)
                    vAbs_u = dsqrt(vxAv_u**2.D0 + vy(i,j+1,k)**2.D0 + vzAv_u**2.D0)
                    vAbs_f = dsqrt(vxAv_f**2.D0 + vyAv_f**2.D0 + vz(i,j,k)**2.D0)
                    vAbs_b = dsqrt(vxAv_b**2.D0 + vyAv_b**2.D0 + vz(i,j,k+1)**2.D0)

                    ForchhTerm_l = rhof*1.75D0/dsqrt(1.5D2*poroHarmX(i,j,k)**3.D0)/dsqrt(KxxHarm(i,j,k))*vAbs_l
                    ForchhTerm_r = rhof*1.75D0/dsqrt(1.5D2*poroHarmX(i+1,j,k)**3.D0)/dsqrt(KxxHarm(i+1,j,k))*vAbs_r
                    ForchhTerm_d = rhof*1.75D0/dsqrt(1.5D2*poroHarmY(i,j,k)**3.D0)/dsqrt(KyyHarm(i,j,k))*vAbs_d
                    ForchhTerm_u = rhof*1.75D0/dsqrt(1.5D2*poroHarmY(i,j+1,k)**3.D0)/dsqrt(KyyHarm(i,j+1,k))*vAbs_u
                    ForchhTerm_f = rhof*1.75D0/dsqrt(1.5D2*poroHarmZ(i,j,k)**3.D0)/dsqrt(KzzHarm(i,j,k))*vAbs_f
                    ForchhTerm_b = rhof*1.75D0/dsqrt(1.5D2*poroHarmZ(i,j,k+1)**3.D0)/dsqrt(KzzHarm(i,j,k+1))*vAbs_b

                    if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1,zlower+k-1)==0)) then
                        div_l = vxBdryX0(ylower+j-1,zlower+k-1)
                    else
                        div_l = (-dpdx(i,j,k)+rhof*gravX)/(ForchhTerm_l+visc/KxxHarm(i,j,k))
                    end if
                    if((pcol==pncols).and.(i==localncols).and.(isDiriX1_p(ylower+j-1,zlower+k-1)==0)) then
                        div_r = vxBdryX1(ylower+j-1,zlower+k-1)
                    else
                        div_r = (-dpdx(i+1,j,k)+rhof*gravX)/(ForchhTerm_r+visc/KxxHarm(i+1,j,k))
                    end if
                    if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1,zlower+k-1)==0)) then
                        div_d = vyBdryY0(xlower+i-1,zlower+k-1)
                    else
                        div_d = (-dpdy(i,j,k)+rhof*gravY)/(ForchhTerm_d+visc/KyyHarm(i,j,k))
                    end if
                    if((prow==pnrows).and.(j==localnrows).and.(isDiriY1_p(xlower+i-1,zlower+k-1)==0)) then
                        div_u = vyBdryY1(xlower+i-1,zlower+k-1)
                    else
                        div_u = (-dpdy(i,j+1,k)+rhof*gravY)/(ForchhTerm_u+visc/KyyHarm(i,j+1,k))
                    end if
                    if((play==1).and.(k==1).and.(isDiriZ0_p(xlower+i-1,ylower+j-1)==0)) then
                        div_f = vzBdryZ0(xlower+i-1,ylower+j-1)
                    else
                        div_f = (-dpdz(i,j,k)+rhof*gravZ)/(ForchhTerm_f+visc/KzzHarm(i,j,k))
                    end if
                    if((play==pnlays).and.(k==localnlays).and.(isDiriZ1_p(xlower+i-1,ylower+j-1)==0)) then
                        div_b = vzBdryZ1(xlower+i-1,ylower+j-1)
                    else
                        div_b = (-dpdz(i,j,k+1)+rhof*gravZ)/(ForchhTerm_b+visc/KzzHarm(i,j,k+1))
                    end if

                    div = (div_r-div_l)/hx(i) + (div_u-div_d)/hy(j) + (div_b-div_f)/hz(k)

                    resi(i,j,k) = src(xlower+i-1,ylower+j-1,zlower+k-1)-div-(poro(i,j,k)-poro_old(i,j,k))/(timeEnd/nt)
                
                end do
            end do
        end do

        deallocate(dpdx)
        deallocate(dpdy)
        deallocate(dpdz)

    end subroutine Resi_mass_p

    ! End Calc Resi Mass

    ! Calc Resi Concentration

    subroutine Resi_concen_cf(conc, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: conc
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi

        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        real(kind=8) :: dl, dt
        real(kind=8), dimension(:,:,:,:), allocatable :: Ex, Ey, Ez, Dx, Dy, Dz
        real(kind=8), dimension(:,:,:), allocatable :: Cfbarx, Cfbary, Cfbarz
        real(kind=8) :: vxAv, vyAv, vzAv, vmodulus
        real(kind=8) :: dmbarx, dmbary, dmbarz
        real(kind=8) :: div1, div2
        real(kind=8), dimension(:,:,:), allocatable :: dCfdx, dCfdy, dCfdz
        real(kind=8) :: dCfdxleft, dCfdxright, dCfdxdown, dCfdxup, dCfdxfront, dCfdxback, dCfdyleft, dCfdyright, dCfdydown, &!
            dCfdyup, dCfdyfront, dCfdyback, dCfdzleft, dCfdzright, dCfdzdown, dCfdzup, dCfdzfront, dCfdzback
        real(kind=8) :: reaction
        integer :: i, j, k

        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays

        allocate(Ex(indexl:indexr,indexd:indexu,indexf:indexb,3))
        allocate(Dx(indexl:indexr,indexd:indexu,indexf:indexb,3))
        allocate(Cfbarx(indexl:indexr,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((pcol==1).and.(i==1)) then
                        vyAv = 2.5D-1 * (vyBdryX0(ylower+j-1,zlower+k-1)+vyBdryX0(ylower+j,zlower+k-1)+vy(i,j,k)+vy(i,j+1,k))
                        vzAv = 2.5D-1 * (vzBdryX0(ylower+j-1,zlower+k-1)+vzBdryX0(ylower+j-1,zlower+k)+vz(i,j,k)+vz(i,j,k+1))
                    elseif((pcol==pncols).and.(i == localncols+1)) then
                        vyAv = 2.5D-1 * (vy(i-1,j,k)+vy(i-1,j+1,k)+vyBdryX1(ylower+j-1,zlower+k-1)+vyBdryX1(ylower+j,zlower+k-1))
                        vzAv = 2.5D-1 * (vz(i-1,j,k)+vz(i-1,j,k+1)+vzBdryX1(ylower+j-1,zlower+k-1)+vzBdryX1(ylower+j-1,zlower+k))
                    else
                        vyAv = 2.5D-1 * (vy(i-1,j,k)+vy(i-1,j+1,k)+vy(i,j,k)+vy(i,j+1,k))
                        vzAv = 2.5D-1 * (vz(i-1,j,k)+vz(i-1,j,k+1)+vz(i,j,k)+vz(i,j,k+1))
                    end if

                    vmodulus = dsqrt((vx(i,j,k)**2.D0+vyAv**2.D0+vzAv**2.D0)*1.D0)

                    if(vmodulus /= 0.D0) then
                        Ex(i,j,k,1) = vx(i,j,k)**2.D0 / vmodulus**2.D0
                        Ex(i,j,k,2) = vx(i,j,k)*vyAv / vmodulus**2.D0  ! the second column and the first row of the matrix
                        Ex(i,j,k,3) = vx(i,j,k)*vzAv / vmodulus**2.D0
                    else
                        Ex(i,j,k,1) = 0.D0
                        Ex(i,j,k,2) = 0.D0
                        Ex(i,j,k,3) = 0.D0
                    end if

                    if((pcol==1).and.(i==1)) then
                        dmbarx = dm(i,j,k)
                    elseif((pcol==pncols).and.(i==localncols+1)) then
                        dmbarx = dm(i-1,j,k)
                    elseif(vx(i,j,k) > 0.D0) then
                        dmbarx = dm(i-1,j,k)
                    else
                        dmbarx = dm(i,j,k)
                    end if

                    dl = alphaOS*dmbarx + 2.D0*lamdaX*vmodulus*radiusInit*(1.D0-poroHarmXInit(i,j,k))/ &!
                        (poroHarmXInit(i,j,k)*(1.D0-poroHarmX(i,j,k)))
                    dt = alphaOS*dmbarx + 2.D0*lamdaT*vmodulus*radiusInit*(1.D0-poroHarmXInit(i,j,k))/ &!
                        (poroHarmXInit(i,j,k)*(1.D0-poroHarmX(i,j,k)))

                    Dx(i,j,k,1) = (dmbarx+vmodulus*dt) + vmodulus*(dl-dt)*Ex(i,j,k,1)
                    Dx(i,j,k,2) = vmodulus*(dl-dt)*Ex(i,j,k,2)
                    Dx(i,j,k,3) = vmodulus*(dl-dt)*Ex(i,j,k,3)

                end do
            end do
        end do

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1,zlower+k-1)==0)) then
                        Cfbarx(i,j,k) = conc(i,j,k)
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1,zlower+k-1)==1)) then
                        Cfbarx(i,j,k) = CfBdryX0(ylower+j-1,zlower+k-1)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1,zlower+k-1)==0)) then
                        Cfbarx(i,j,k) = conc(i-1,j,k)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1,zlower+k-1)==1)) then
                        Cfbarx(i,j,k) = CfBdryX1(ylower+j-1,zlower+k-1)
                    elseif(vx(i,j,k) > 0.D0) then
                        Cfbarx(i,j,k) = conc(i-1,j,k)
                    else
                        Cfbarx(i,j,k) = conc(i,j,k)
                    end if

                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1
        indexf = 1
        indexb = localnlays

        allocate(Ey(indexl:indexr,indexd:indexu,indexf:indexb,3))
        allocate(Dy(indexl:indexr,indexd:indexu,indexf:indexb,3))
        allocate(Cfbary(indexl:indexr,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((prow==1).and.(j==1)) then
                        vxAv = 2.5D-1 * (vxBdryY0(xlower+i-1,zlower+k-1)+vxBdryY0(xlower+i,zlower+k-1)+vx(i,j,k)+vx(i+1,j,k))
                        vzAv = 2.5D-1 * (vzBdryY0(xlower+i-1,zlower+k-1)+vzBdryY0(xlower+i-1,zlower+k)+vz(i,j,k)+vz(i,j,k+1))
                    elseif((prow==pnrows).and.(j==localnrows+1)) then
                        vxAv = 2.5D-1 * (vx(i,j-1,k)+vx(i+1,j-1,k)+vxBdryY1(xlower+i-1,zlower+k-1)+vxBdryY1(xlower+i,zlower+k-1))
                        vzAv = 2.5D-1 * (vz(i,j-1,k)+vz(i,j-1,k+1)+vzBdryY1(xlower+i-1,zlower+k-1)+vzBdryY1(xlower+i-1,zlower+k))
                    else
                        vxAv = 2.5D-1 * (vx(i,j-1,k)+vx(i+1,j-1,k)+vx(i,j,k)+vx(i+1,j,k))
                        vzAv = 2.5D-1 * (vz(i,j-1,k)+vz(i,j-1,k+1)+vz(i,j,k)+vz(i,j,k+1))
                    end if

                    vmodulus = dsqrt((vxAv**2.D0+vy(i,j,k)**2.D0+vzAv**2.D0)*1.D0)

                    if(vmodulus /= 0.D0) then
                        Ey(i,j,k,1) = vy(i,j,k)*vxAv / vmodulus**2.D0
                        Ey(i,j,k,2) = vy(i,j,k)**2.D0 / vmodulus**2.D0
                        Ey(i,j,k,3) = vy(i,j,k)*vzAv / vmodulus**2.D0
                    else
                        Ey(i,j,k,1) = 0.D0
                        Ey(i,j,k,2) = 0.D0
                        Ey(i,j,k,3) = 0.D0
                    end if

                    if((prow==1).and.(j==1)) then
                        dmbary = dm(i,j,k)
                    elseif((prow==pnrows).and.(j==localnrows+1)) then
                        dmbary = dm(i,j-1,k)
                    elseif(vy(i,j,k) > 0.D0) then
                        dmbary = dm(i,j-1,k)
                    else
                        dmbary = dm(i,j,k)
                    end if

                    dl = alphaOS*dmbary + 2.D0*lamdaX*vmodulus*radiusInit*(1.D0-poroHarmYInit(i,j,k))/ &!
                        (poroHarmYInit(i,j,k)*(1.D0-poroHarmY(i,j,k)))
                    dt = alphaOS*dmbary + 2.D0*lamdaT*vmodulus*radiusInit*(1.D0-poroHarmYInit(i,j,k))/ &!
                        (poroHarmYInit(i,j,k)*(1.D0-poroHarmY(i,j,k)))

                    Dy(i,j,k,1) = vmodulus*(dl-dt)*Ey(i,j,k,1)
                    Dy(i,j,k,2) = (dmbary+vmodulus*dt) + vmodulus*(dl-dt)*Ey(i,j,k,2)
                    Dy(i,j,k,3) = vmodulus*(dl-dt)*Ey(i,j,k,3)

                end do
            end do
        end do

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1,zlower+k-1)==0)) then
                        Cfbary(i,j,k) = conc(i,j,k)
                    elseif((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1,zlower+k-1)==1)) then
                        Cfbary(i,j,k) = CfBdryY0(xlower+i-1,zlower+k-1)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1,zlower+k-1)==0)) then
                        Cfbary(i,j,k) = conc(i,j-1,k)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1,zlower+k-1)==1)) then
                        Cfbary(i,j,k) = CfBdryY1(xlower+i-1,zlower+k-1)
                    elseif(vy(i,j,k) > 0.D0) then
                        Cfbary(i,j,k) = conc(i,j-1,k)
                    else
                        Cfbary(i,j,k) = conc(i,j,k)
                    end if

                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays + 1

        allocate(Ez(indexl:indexr,indexd:indexu,indexf:indexb,3))
        allocate(Dz(indexl:indexr,indexd:indexu,indexf:indexb,3))
        allocate(Cfbarz(indexl:indexr,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((play==1).and.(k==1)) then
                        vxAv = 2.5D-1 * (vxBdryZ0(xlower+i-1,ylower+j-1)+vxBdryZ0(xlower+i,ylower+j-1)+vx(i,j,k)+vx(i+1,j,k))
                        vyAv = 2.5D-1 * (vyBdryZ0(xlower+i-1,ylower+j-1)+vyBdryZ0(xlower+i-1,ylower+j)+vy(i,j,k)+vy(i,j+1,k))
                    elseif((play==pnlays).and.(k==localnlays+1)) then
                        vxAv = 2.5D-1 * (vx(i,j,k-1)+vx(i+1,j,k-1)+vxBdryZ1(xlower+i-1,ylower+j-1)+vxBdryZ1(xlower+i,ylower+j-1))
                        vyAv = 2.5D-1 * (vy(i,j,k-1)+vy(i,j+1,k-1)+vyBdryZ1(xlower+i-1,ylower+j-1)+vyBdryZ1(xlower+i-1,ylower+j))
                    else
                        vxAv = 2.5D-1 * (vx(i,j,k-1)+vx(i+1,j,k-1)+vx(i,j,k)+vx(i+1,j,k))
                        vyAv = 2.5D-1 * (vy(i,j,k-1)+vy(i,j+1,k-1)+vy(i,j,k)+vy(i,j+1,k))
                    end if

                    vmodulus = dsqrt((vxAv**2.D0+vyAv**2.D0+vz(i,j,k)**2.D0)*1.D0)

                    if(vmodulus /= 0.D0) then
                        Ez(i,j,k,1) = vz(i,j,k)*vxAv / vmodulus**2.D0
                        Ez(i,j,k,2) = vz(i,j,k)*vyAv / vmodulus**2.D0
                        Ez(i,j,k,3) = vz(i,j,k)**2.D0 / vmodulus**2.D0
                    else
                        Ez(i,j,k,1) = 0.D0
                        Ez(i,j,k,2) = 0.D0
                        Ez(i,j,k,3) = 0.D0
                    end if

                    if((play==1).and.(k==1)) then
                        dmbarz = dm(i,j,k)
                    elseif((play==pnlays).and.(k==localnlays+1)) then
                        dmbarz = dm(i,j,k-1)
                    elseif(vz(i,j,k) > 0.D0) then
                        dmbarz = dm(i,j,k-1)
                    else
                        dmbarz = dm(i,j,k)
                    end if

                    dl = alphaOS*dmbarz + 2.D0*lamdaX*vmodulus*radiusInit*(1.D0-poroHarmZInit(i,j,k))/ &!
                        (poroHarmZInit(i,j,k)*(1.D0-poroHarmZ(i,j,k)))
                    dt = alphaOS*dmbarz + 2.D0*lamdaT*vmodulus*radiusInit*(1.D0-poroHarmZInit(i,j,k))/ &!
                        (poroHarmZInit(i,j,k)*(1.D0-poroHarmZ(i,j,k)))

                    Dz(i,j,k,1) = vmodulus*(dl-dt)*Ez(i,j,k,1)
                    Dz(i,j,k,2) = vmodulus*(dl-dt)*Ez(i,j,k,2)
                    Dz(i,j,k,3) = (dmbarz+vmodulus*dt) + vmodulus*(dl-dt)*Ez(i,j,k,3)

                end do
            end do
        end do

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((play==1).and.(k==1).and.(isDiriZ0_Cf(xlower+i-1,ylower+j-1)==0)) then
                        Cfbarz(i,j,k) = conc(i,j,k)
                    elseif((play==1).and.(k==1).and.(isDiriZ0_Cf(xlower+i-1,ylower+j-1)==1)) then
                        Cfbarz(i,j,k) = CfBdryZ0(xlower+i-1,ylower+j-1)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_Cf(xlower+i-1,ylower+j-1)==0)) then
                        Cfbarz(i,j,k) = conc(i,j,k-1)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_Cf(xlower+i-1,ylower+j-1)==1)) then
                        Cfbarz(i,j,k) = CfBdryZ1(xlower+i-1,ylower+j-1)
                    elseif(vz(i,j,k) > 0.D0) then
                        Cfbarz(i,j,k) = conc(i,j,k-1)
                    else
                        Cfbarz(i,j,k) = conc(i,j,k)
                    end if

                end do
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
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        if(play /= pnlays) then
            indexb = localnlays + 1
        else
            indexb = localnlays
        end if
        allocate(dCfdx(1:localncols+1,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = 1, localncols+1

                    if((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1,zlower+k-1)==0)) then
                        dCfdx(i,j,k) = 0.D0
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1,zlower+k-1)==1)) then
                        dCfdx(i,j,k) = (conc(i,j,k)-CfBdryX0(ylower+j-1,zlower+k-1))/hx(i)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1,zlower+k-1)==0)) then
                        dCfdx(i,j,k) = 0.D0
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1,zlower+k-1)==1)) then
                        dCfdx(i,j,k) = (CfBdryX1(ylower+j-1,zlower+k-1)-conc(i-1,j,k))/hx(i-1)
                    else
                        dCfdx(i,j,k) = (conc(i,j,k)-conc(i-1,j,k))/((hx(i)+hx(i-1))/2.D0)
                    end if

                end do
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
        if(play /= 1) then
            indexf = 0
        else
            indexf = 1
        end if
        if(play /= pnlays) then
            indexb = localnlays + 1
        else
            indexb = localnlays
        end if
        allocate(dCfdy(indexl:indexr,1:localnrows+1,indexf:indexb))

        do k = indexf, indexb
            do j = 1, localnrows+1
                do i = indexl, indexr

                    if((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1,zlower+k-1)==0)) then
                        dCfdy(i,j,k) = 0.D0
                    elseif((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1,zlower+k-1)==1)) then
                        dCfdy(i,j,k) = (conc(i,j,k)-CfBdryY0(xlower+i-1,zlower+k-1))/hy(j)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1,zlower+k-1)==0)) then
                        dCfdy(i,j,k) = 0.D0
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1,zlower+k-1)==1)) then
                        dCfdy(i,j,k) = (CfBdryY1(xlower+i-1,zlower+k-1)-conc(i,j-1,k))/hy(j-1)
                    else
                        dCfdy(i,j,k) = (conc(i,j,k)-conc(i,j-1,k))/((hy(j)+hy(j-1))/2.D0)
                    end if

                end do
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
        allocate(dCfdz(indexl:indexr,indexd:indexu,1:localnlays+1))

        do k = 1, localnlays+1
            do j = indexd, indexu
                do i = indexl, indexr

                    if((play==1).and.(k==1).and.(isDiriZ0_Cf(xlower+i-1,ylower+j-1)==0)) then
                        dCfdz(i,j,k) = 0.D0
                    elseif((play==1).and.(k==1).and.(isDiriZ0_Cf(xlower+i-1,ylower+j-1)==1)) then
                        dCfdz(i,j,k) = (conc(i,j,k)-CfBdryZ0(xlower+i-1,ylower+j-1))/hz(k)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_Cf(xlower+i-1,ylower+j-1)==0)) then
                        dCfdz(i,j,k) = 0.D0
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_Cf(xlower+i-1,ylower+j-1)==1)) then
                        dCfdz(i,j,k) = (CfBdryZ1(xlower+i-1,ylower+j-1)-conc(i,j,k-1))/hz(k-1)
                    else
                        dCfdz(i,j,k) = (conc(i,j,k)-conc(i,j,k-1))/((hz(k)+hz(k-1))/2.D0)
                    end if

                end do
            end do
        end do

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    div1 = (vx(i+1,j,k)*Cfbarx(i+1,j,k) - vx(i,j,k)*Cfbarx(i,j,k))/hx(i) + &!
                        (vy(i,j+1,k)*Cfbary(i,j+1,k) - vy(i,j,k)*Cfbary(i,j,k))/hy(j) + &!
                        (vz(i,j,k+1)*Cfbarz(i,j,k+1) - vz(i,j,k)*Cfbarz(i,j,k))/hz(k)

                    dCfdxleft = dCfdx(i,j,k)
                    dCfdxright = dCfdx(i+1,j,k)
                    if((prow==1).and.(j==1)) then
                        dCfdxdown = (dCfdx(i,j,k)+dCfdx(i+1,j,k)) / 2.D0
                    else
                        dCfdxdown = (dCfdx(i,j,k)+dCfdx(i+1,j,k)+dCfdx(i,j-1,k)+dCfdx(i+1,j-1,k)) / 4.D0
                    end if
                    if((prow==pnrows).and.(j==localnrows)) then
                        dCfdxup = (dCfdx(i,j,k)+dCfdx(i+1,j,k)) / 2.D0
                    else
                        dCfdxup = (dCfdx(i,j,k)+dCfdx(i+1,j,k)+dCfdx(i,j+1,k)+dCfdx(i+1,j+1,k)) / 4.D0
                    end if
                    if((play==1).and.(k==1)) then
                        dCfdxfront = (dCfdx(i,j,k)+dCfdx(i+1,j,k)) / 2.D0
                    else
                        dCfdxfront = (dCfdx(i,j,k)+dCfdx(i+1,j,k)+dCfdx(i,j,k-1)+dCfdx(i+1,j,k-1)) / 4.D0
                    end if
                    if((play==pnlays).and.(k==localnlays)) then
                        dCfdxback = (dCfdx(i,j,k)+dCfdx(i+1,j,k)) / 2.D0
                    else
                        dCfdxback = (dCfdx(i,j,k)+dCfdx(i+1,j,k)+dCfdx(i,j,k+1)+dCfdx(i+1,j,k+1)) / 4.D0
                    end if

                    if((pcol==1).and.(i==1)) then
                        dCfdyleft = (dCfdy(i,j,k)+dCfdy(i,j+1,k)) / 2.D0
                    else
                        dCfdyleft = (dCfdy(i-1,j,k)+dCfdy(i-1,j+1,k)+dCfdy(i,j,k)+dCfdy(i,j+1,k)) / 4.D0
                    end if
                    if((pcol==pncols).and.(i==localncols)) then
                        dCfdyright = (dCfdy(i,j,k)+dCfdy(i,j+1,k)) / 2.D0
                    else
                        dCfdyright = (dCfdy(i,j,k)+dCfdy(i,j+1,k)+dCfdy(i+1,j,k)+dCfdy(i+1,j+1,k)) / 4.D0
                    end if
                    dCfdydown = dCfdy(i,j,k)
                    dCfdyup = dCfdy(i,j+1,k)
                    if((play==1).and.(k==1)) then
                        dCfdyfront = (dCfdy(i,j,k)+dCfdy(i,j+1,k)) / 2.D0
                    else
                        dCfdyfront = (dCfdy(i,j,k-1)+dCfdy(i,j+1,k-1)+dCfdy(i,j,k)+dCfdy(i,j+1,k)) / 4.D0
                    end if
                    if((play==pnlays).and.(k==localnlays)) then
                        dCfdyback = (dCfdy(i,j,k)+dCfdy(i,j+1,k)) / 2.D0
                    else
                        dCfdyback = (dCfdy(i,j,k)+dCfdy(i,j+1,k)+dCfdy(i,j,k+1)+dCfdy(i,j+1,k+1)) / 4.D0
                    end if

                    if((pcol==1).and.(i==1)) then
                        dCfdzleft = (dCfdz(i,j,k)+dCfdz(i,j,k+1)) / 2.D0
                    else
                        dCfdzleft = (dCfdz(i-1,j,k)+dCfdz(i-1,j,k+1)+dCfdz(i,j,k)+dCfdz(i,j,k+1)) / 4.D0
                    end if
                    if((pcol==pncols).and.(i==localncols)) then
                        dCfdzright = (dCfdz(i,j,k)+dCfdz(i,j,k+1)) / 2.D0
                    else
                        dCfdzright = (dCfdz(i,j,k)+dCfdz(i,j,k+1)+dCfdz(i+1,j,k)+dCfdz(i+1,j,k+1)) / 4.D0
                    end if
                    if((prow==1).and.(j==1)) then
                        dCfdzdown = (dCfdz(i,j,k)+dCfdz(i,j,k+1)) / 2.D0
                    else
                        dCfdzdown = (dCfdz(i,j-1,k)+dCfdz(i,j-1,k+1)+dCfdz(i,j,k)+dCfdz(i,j,k+1)) / 4.D0
                    end if
                    if((prow==pnrows).and.(j==localnrows)) then
                        dCfdzup = (dCfdz(i,j,k)+dCfdz(i,j,k+1)) / 2.D0
                    else
                        dCfdzup = (dCfdz(i,j,k)+dCfdz(i,j,k+1)+dCfdz(i,j+1,k)+dCfdz(i,j+1,k+1)) / 4.D0
                    end if
                    dCfdzfront = dCfdz(i,j,k)
                    dCfdzback = dCfdz(i,j,k+1)

                    div2 = (dCfdxright*poroHarmX(i+1,j,k)*Dx(i+1,j,k,1) - dCfdxleft*poroHarmX(i,j,k)*Dx(i,j,k,1))/hx(i) + &!
                        (dCfdyright*poroHarmX(i+1,j,k)*Dx(i+1,j,k,2) - dCfdyleft*poroHarmX(i,j,k)*Dx(i,j,k,2))/hx(i) + &!
                        (dCfdzright*poroHarmX(i+1,j,k)*Dx(i+1,j,k,3) - dCfdzleft*poroHarmX(i,j,k)*Dx(i,j,k,3))/hx(i) + &!
                        (dCfdxup*poroHarmY(i,j+1,k)*Dy(i,j+1,k,1) - dCfdxdown*poroHarmY(i,j,k)*Dy(i,j,k,1))/hy(j) + &!
                        (dCfdyup*poroHarmY(i,j+1,k)*Dy(i,j+1,k,2) - dCfdydown*poroHarmY(i,j,k)*Dy(i,j,k,2))/hy(j) + &!
                        (dCfdzup*poroHarmY(i,j+1,k)*Dy(i,j+1,k,3) - dCfdzdown*poroHarmY(i,j,k)*Dy(i,j,k,3))/hy(j) + &!
                        (dCfdxback*poroHarmZ(i,j,k+1)*Dz(i,j,k+1,1) - dCfdxfront*poroHarmZ(i,j,k)*Dz(i,j,k,1))/hz(k) + &!
                        (dCfdyback*poroHarmZ(i,j,k+1)*Dz(i,j,k+1,2) - dCfdyfront*poroHarmZ(i,j,k)*Dz(i,j,k,2))/hz(k) + &!
                        (dCfdzback*poroHarmZ(i,j,k+1)*Dz(i,j,k+1,3) - dCfdzfront*poroHarmZ(i,j,k)*Dz(i,j,k,3))/hz(k)

                    reaction = av(i,j,k)*conc(i,j,k)*kc(i,j,k)*ks(i,j,k)/(kc(i,j,k)+ks(i,j,k))

                    resi(i,j,k) = src(xlower+i-1,ylower+j-1,zlower+k-1)-reaction+div2-div1-(conc(i,j,k)*poro(i,j,k)-Cf(i,j,k)* &!
                        poro_old(i,j,k))/(timeEnd/nt)

                end do
            end do
        end do

        deallocate(Ex)
        deallocate(Ey)
        deallocate(Ez)
        deallocate(Dx)
        deallocate(Dy)
        deallocate(Dz)
        deallocate(Cfbarx)
        deallocate(Cfbary)
        deallocate(Cfbarz)
        deallocate(dCfdx)
        deallocate(dCfdy)
        deallocate(dCfdz)

    end subroutine Resi_concen_cf

    ! End Calc Resi Concentration

    ! Calc Resi Temperature

    subroutine Resi_ener_tem(tempe, resi)

        integer, dimension(:,:,:), pointer, intent(in) :: tempe
        real(kind=8), dimension(:,:,:), pointer, intent(in out) :: resi

        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        real(kind=8), dimension(:,:,:), allocatable :: Tembarx, Tembary, Tembarz
        real(kind=8), dimension(:,:,:), allocatable :: dTemdx, dTemdy, dTemdz
        real(kind=8) :: dvdt, tp, vmodulus, v1, v2, Forchh, HrT, arh
        real(kind=8) :: div1, div2, div3
        integer :: i, j, k

        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays
        allocate(Tembarx(indexl:indexr,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1,zlower+k-1)==0)) then
                        Tembarx(i,j,k) = tempe(i,j,k)
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1,zlower+k-1)==1)) then
                        Tembarx(i,j,k) = TemBdryX0(ylower+j-1,zlower+k-1)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1,zlower+k-1)==0)) then
                        Tembarx(i,j,k) = tempe(i-1,j,k)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1,zlower+k-1)==1)) then
                        Tembarx(i,j,k) = TemBdryX1(ylower+j-1,zlower+k-1)
                    elseif(vx(i,j,k) > 0.D0) then
                        Tembarx(i,j,k) = tempe(i-1,j,k)
                    else
                        Tembarx(i,j,k) = tempe(i,j,k)
                    end if

                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1
        indexf = 1
        indexb = localnlays
        allocate(Tembary(indexl:indexr,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1,zlower+k-1)==0)) then
                        Tembary(i,j,k) = tempe(i,j,k)
                    elseif((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1,zlower+k-1)==1)) then
                        Tembary(i,j,k) = TemBdryY0(xlower+i-1,zlower+k-1)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1,zlower+k-1)==0)) then
                        Tembary(i,j,k) = tempe(i,j-1,k)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1,zlower+k-1)==1)) then
                        Tembary(i,j,k) = TemBdryY1(xlower+i-1,zlower+k-1)
                    elseif(vy(i,j,k) > 0.D0) then
                        Tembary(i,j,k) = tempe(i,j-1,k)
                    else
                        Tembary(i,j,k) = tempe(i,j,k)
                    end if

                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays + 1
        allocate(Tembarz(indexl:indexr,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((play==1).and.(k==1).and.(isDiriZ0_Tem(xlower+i-1,ylower+j-1)==0)) then
                        Tembarz(i,j,k) = tempe(i,j,k)
                    elseif((play==1).and.(k==1).and.(isDiriZ0_Tem(xlower+i-1,ylower+j-1)==1)) then
                        Tembarz(i,j,k) = TemBdryZ0(xlower+i-1,ylower+j-1)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_Tem(xlower+i-1,ylower+j-1)==0)) then
                        Tembarz(i,j,k) = tempe(i,j,k-1)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_Tem(xlower+i-1,ylower+j-1)==1)) then
                        Tembarz(i,j,k) = TemBdryZ1(xlower+i-1,ylower+j-1)
                    elseif(vz(i,j,k) > 0.D0) then
                        Tembarz(i,j,k) = tempe(i,j,k-1)
                    else
                        Tembarz(i,j,k) = tempe(i,j,k)
                    end if

                end do
            end do
        end do


        indexd = 1
        indexu = localnrows
        indexf = 1
        indexb = localnlays
        allocate(dTemdx(1:localncols+1,indexd:indexu,indexf:indexb))

        do k = indexf, indexb
            do j = indexd, indexu
                do i = 1, localncols+1

                    if((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1,zlower+k-1)==0)) then
                        dTemdx(i,j,k) = 0.D0
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1,zlower+k-1)==1)) then
                        dTemdx(i,j,k) = (tempe(i,j,k)-TemBdryX0(ylower+j-1,zlower+k-1))/hx(i)
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1,zlower+k-1)==0)) then
                        dTemdx(i,j,k) = 0.D0
                    elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1,zlower+k-1)==1)) then
                        dTemdx(i,j,k) = (TemBdryX1(ylower+j-1,zlower+k-1)-tempe(i-1,j,k))/hx(i-1)
                    else
                        dTemdx(i,j,k) = (tempe(i,j,k)-tempe(i-1,j,k))/((hx(i)+hx(i-1))/2.D0)
                    end if

                end do
            end do
        end do


        indexl = 1
        indexr = localncols
        indexf = 1
        indexb = localnlays
        allocate(dTemdy(indexl:indexr,1:localnrows+1,indexf:indexb))

        do k = indexf, indexb
            do j = 1, localnrows+1
                do i = indexl, indexr

                    if((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1,zlower+k-1)==0)) then
                        dTemdy(i,j,k) = 0.D0
                    elseif((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1,zlower+k-1)==1)) then
                        dTemdy(i,j,k) = (tempe(i,j,k)-TemBdryY0(xlower+i-1,zlower+k-1))/hy(j)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1,zlower+k-1)==0)) then
                        dTemdy(i,j,k) = 0.D0
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1,zlower+k-1)==1)) then
                        dTemdy(i,j,k) = (TemBdryY1(xlower+i-1,zlower+k-1)-tempe(i,j-1,k))/hy(j-1)
                    else
                        dTemdy(i,j,k) = (tempe(i,j,k)-tempe(i,j-1,k))/((hy(j)+hy(j-1))/2.D0)
                    end if

                end do
            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        allocate(dTemdz(indexl:indexr,indexd:indexu,1:localnlays+1))

        do k = 1, localnlays+1
            do j = indexd, indexu
                do i = indexl, indexr

                    if((play==1).and.(k==1).and.(isDiriZ0_Tem(xlower+i-1,ylower+j-1)==0)) then
                        dTemdz(i,j,k) = 0.D0
                    elseif((play==1).and.(k==1).and.(isDiriZ0_Tem(xlower+i-1,ylower+j-1)==1)) then
                        dTemdz(i,j,k) = (tempe(i,j,k)-TemBdryZ0(xlower+i-1,ylower+j-1))/hz(k)
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_Tem(xlower+i-1,ylower+j-1)==0)) then
                        dTemdz(i,j,k) = 0.D0
                    elseif((play==pnlays).and.(k==localnlays+1).and.(isDiriZ1_Tem(xlower+i-1,ylower+j-1)==1)) then
                        dTemdz(i,j,k) = (TemBdryZ1(xlower+i-1,ylower+j-1)-tempe(i,j,k-1))/hz(k-1)
                    else
                        dTemdz(i,j,k) = (tempe(i,j,k)-tempe(i,j,k-1))/((hz(k)+hz(k-1))/2.D0)
                    end if

                end do
            end do
        end do

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    dvdt = (poro(i,j,k)*rhof*thetaf*tempe(i,j,k)+(1-poro(i,j,k))*rhos*thetas*tempe(i,j,k) - &!
                        poro_old(i,j,k)*rhof*thetaf*Tem(i,j,k)-(1-poro_old(i,j,k))*rhos*thetas*Tem(i,j,k))/(timeEnd/nt)

                    div1 = rhof*thetaf * ((vx(i+1,j,k)*Tembarx(i+1,j,k) - vx(i,j,k)*Tembarx(i,j,k))/hx(i) + &!
                        (vy(i,j+1,k)*Tembary(i,j+1,k) - vy(i,j,k)*Tembary(i,j,k))/hy(j) + &!
                        (vz(i,j,k+1)*Tembarz(i,j,k+1) - vz(i,j,k)*Tembarz(i,j,k))/hz(k))

                    div2 = (dTemdx(i+1,j,k)*(poroHarmX(i+1,j,k)*lamdaf+(1-poroHarmX(i+1,j,k))*lamdas) - &!
                        dTemdx(i,j,k)*(poroHarmX(i,j,k)*lamdaf+(1-poroHarmX(i,j,k))*lamdas))/hx(i) + &!

                        (dTemdy(i,j+1,k)*(poroHarmY(i,j+1,k)*lamdaf+(1-poroHarmY(i,j+1,k))*lamdas) - &!
                        dTemdy(i,j,k)*(poroHarmY(i,j,k)*lamdaf+(1-poroHarmY(i,j,k))*lamdas))/hy(j) + &!

                        (dTemdz(i,j,k+1)*(poroHarmZ(i,j,k+1)*lamdaf+(1-poroHarmZ(i,j,k+1))*lamdas) - &!
                        dTemdz(i,j,k)*(poroHarmZ(i,j,k)*lamdaf+(1-poroHarmZ(i,j,k))*lamdas))/hz(k)

                    div3 = p(i,j,k) * ((vx(i+1,j,k)-vx(i,j,k))/hx(i) + (vy(i,j+1,k)-vy(i,j,k))/hy(j) + &!
                        (vz(i,j,k+1)-vz(i,j,k))/hz(k))

                    tp = visc * (((vx(i+1,j,k)/poroHarmX(i+1,j,k)-vx(i,j,k)/poroHarmX(i,j,k))/hx(i))**2.D0 + &!
                        ((vy(i,j+1,k)/poroHarmY(i,j+1,k)-vy(i,j,k)/poroHarmY(i,j,k))/hy(j))**2.D0 + &!
                        ((vz(i,j,k+1)/poroHarmZ(i,j,k+1)-vz(i,j,k)/poroHarmZ(i,j,k))/hz(k))**2.D0)

                    vmodulus = dsqrt(((vx(i,j,k)+vx(i+1,j,k))/2.D0)**2.D0 + &!
                        ((vy(i,j,k)+vy(i,j+1,k))/2.D0)**2.D0  + ((vz(i,j,k+1)+vz(i,j,k))/2.D0)**2.D0)
                    v1 = visc/Kxx(i,j,k)*vmodulus**2.D0

                    Forchh = 1.75D0/dsqrt(1.5D2*poro(i,j,k)**3.D0)
                    v2 = rhof*Forchh/dsqrt(Kxx(i,j,k))*vmodulus**3.D0

                    HrT = abs(-9.702D3+1.697D1*Tem(i,j,k)-2.34D-3*Tem(i,j,k)**2.D0)
                    arh = av(i,j,k) * kc(i,j,k)*ks(i,j,k)/(kc(i,j,k)+ks(i,j,k))*Cf(i,j,k) * HrT

                    resi(i,j,k) = div2 - div3 + tp + v1 + v2 + arh - div1 - dvdt

                end do
            end do
        end do

        deallocate(Tembarx)
        deallocate(Tembary)
        deallocate(Tembarz)
        deallocate(dTemdx)
        deallocate(dTemdy)
        deallocate(dTemdz)

    end subroutine Resi_ener_tem

    ! End Calc Resi Temperature

end module DBF_resi
