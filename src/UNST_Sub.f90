! Calculation Subroutines
! Coded by K.Kawaike and TK Labo
! Released on July 7th 2025

!-------------------------------
! Calculate flux
!-------------------------------
subroutine flux
    use unst_globals_mod
    implicit none
    real(8) hhe, hhw, hhep, hhwp, hhan, sgn, hh1 !Discontinuous(Drop formula, Overflow formula)
    real(8) u11, v11, u13, v13, sqx, ram !Continuous
    real(8) outa2    !Add storage model_1(Paddy field dam)
    real(8) rem    !Add storage model_2(School ground & Park)
    real(8) f, sumf, dvp1, dvp2, dv    !Equation of Continuity
    real(8) sadj, dvc, dv_cextre !tracer
    real(8) rt
    real(8) plant_force                 ! add plant force
    integer li, k, me1, me2
    integer me, it, ilt, nn, ct, nt
    integer iflg0, iflg
    real(8), allocatable :: rr(:), q(:), gamma(:)
    real(8) blink, h1, h2, uvmn, vol1, vol2
    integer ii
    real(8) temp_dsdepth

    ! equation of motion

    !$omp parallel do default(shared) &
    !$omp private(hhe, hhw, hhep, hhwp, hhan, sgn, hh1, u11, v11, u13, v13, sqx, ram, li, k, me1, me2) &
    !$omp private(plant_force, blink, h1, h2, uvmn, vol1, vol2)
  do li = 1, link
    if(limesh(li, 2) == 0) cycle
    if(inf(limesh(li, 2)) == 0) goto 300
    if(inf(limesh(li, 1)) == 0) goto 300
    if(rbeta(li) == 0.0d0) goto 300
    if(unsth(limesh(li, 1)) <= th .and. unsth(limesh(li, 2)) <= th) goto 300
        hhe = unsth(limesh(li, 2)) + baseo(limesh(li, 2))
        hhw = unsth(limesh(li, 1)) + baseo(limesh(li, 1))
        hhep = unsth(limesh(li, 2)) - th
        hhwp = unsth(limesh(li, 1)) - th

    ! Embankment added by CTI r.nishizawa
    if (morid == 1) then
    if(infl(li) == 1) then
        blink= sqrt((dnox(lnode1(li)) - dnox(lnode2(li))) **2 + (dnox(lnode1(li)) - dnox(lnode2(li)) ** 2))
        h1 = max(hhe, hhw) - max(baseo(limesh(li, 1)), baseo(limesh(li,2)), zbbk(li))
        h2 = min(hhe, hhw) - max(baseo(limesh(li, 1)), baseo(limesh(li,2)), zbbk(li))

        hhan = hhe - hhw
        sgn = hhan/abs(hhan)

        if (h1 > th) then
            if((h2/h1) <= 2.0d0/3.0d0 .and. h1 > th) then
                um(li) = -sgn * 0.35d0 * h1 * sqrt(2.0d0*gg*h1)*ux(li)
                vn(li) = -sgn * 0.35d0 * h1 * sqrt(2.0d0*gg*h1)*uy(li)
            elseif((h2/h1) > 2.0d0/3.0d0 .and. h1 > th) then
                um(li) = -sgn * 0.91d0 * h2 * sqrt(2.0d0*gg*(h1-h2))*ux(li)
                vn(li) = -sgn * 0.91d0 * h2 * sqrt(2.0d0*gg*(h1-h2))*uy(li)
            endif

            uvmn = sqrt(um(li)**2 + vn(li)**2) * blink * dt2
            vol1 = unsth(limesh(li,1)) * smesh(limesh(li,1)) * (1.0d0 - lambda(limesh(li,1)))
            vol2 = unsth(limesh(li,2)) * smesh(limesh(li,2)) * (1.0d0 - lambda(limesh(li,2)))

            if (hhw > hhe .and. uvmn > vol1) then
                um(li) = um(li) * vol1/uvmn
                vn(li) = vn(li) * vol1/uvmn
            elseif (hhe > hhw .and. uvmn > vol2) then
                um(li) = um(li) * vol2/uvmn
                vn(li) = vn(li) * vol2/uvmn
            endif
        else
            um(li) = 0.0d0
            vn(li) = 0.0d0
        endif
        lhan(li) = 1
        cycle
    endif
    endif

    !Paddy field dam by k.kawaike
    if(inf(limesh(li,1))==61 .or. inf(limesh(li,1))==63) then
        if(inf(limesh(li,2))>20 .and. inf(limesh(li,2))<26) then
            if(unsth(limesh(li, 1)) > th .or. hhe > baseo(limesh(li, 1))) then
                outa2 = 0.25d0*pi*phi(inf(limesh(li,1)))**2   !オリフィス(孔)の断面積
                hhan = hhe - hhw
                if(abs(hhan)==0.0d0) goto 300
                sgn = hhan/abs(hhan)
                um(li) = - sgn*0.6d0*outa2*sqrt(2.0d0*gg*abs(hhan))*ux(li)/20.0d0   !流量係数を0.6とおく．また，オリフィスが20m間隔で並んでいるとする．
                vn(li) = - sgn*0.6d0*outa2*sqrt(2.0d0*gg*abs(hhan))*uy(li)/20.0d0   !Let Coefficient of discharge be 0.6. Considered that the holes are lined up every 20m.
                if(sgn>=0.0d0) count2 = count2 + 1
                if(sgn<0.0d0) count1 = count1 + 1
                lhan(li) = 1
            else
                count3 = count3 + 1
                um(li) = 0.0d0
                vn(li) = 0.0d0
                lhan(li) = 1
            endif
        else
            count4 = count4 + 1
            goto 311
        endif
        cycle
    elseif(inf(limesh(li,2))==61 .or. inf(limesh(li,2))==63) then
        if(inf(limesh(li,1))>20 .and. inf(limesh(li,1))<26) then
            if(unsth(limesh(li, 2)) > th .or. hhw > baseo(limesh(li, 2))) then
                outa2 = 0.25d0*pi*phi(inf(limesh(li,2)))**2   !オリフィス(孔)の断面積
                hhan = hhe - hhw
                if(abs(hhan)==0.0d0) goto 300
                sgn = hhan/abs(hhan)
                um(li) = - sgn*0.6d0*outa2*sqrt(2.0d0*gg*abs(hhan))*ux(li)/20.0d0   !流量係数を0.6とおく．また，オリフィスが20m間隔で並んでいるとする．
                vn(li) = - sgn*0.6d0*outa2*sqrt(2.0d0*gg*abs(hhan))*uy(li)/20.0d0   !Let Coefficient of discharge be 0.6. Considered that the holes are lined up every 20m.
                if(sgn>0.0d0) count1 = count1 + 1
                if(sgn<=0.0d0) count2 = count2 + 1
                lhan(li) = 1
            else
                count3 = count3 + 1
                um(li) = 0.0d0
                vn(li) = 0.0d0
                lhan(li) = 1
            endif
        else
            count4 = count4 + 1
            goto 311
        endif
        cycle
311 endif

    ! when the water surface is not continuous
    ! 段落ちの式
    if(hhe < baseo(limesh(li, 1))) then
        if ((inf(limesh(li,1))==71 .and. inf(limesh(li,2))<71) .or. &
            (inf(limesh(li,2))==71 .and. inf(limesh(li,1))<71)) then
            !田んぼダム適用 by k.yamamura(inf=71)
            if(unsth(limesh(li,1)) > lh + th) then
                um(li) = 0.544d0*(unsth(limesh(li, 1)) - lh)*sqrt(gg*(unsth(limesh(li, 1)) - lh))*ux(li)
                vn(li) = 0.544d0*(unsth(limesh(li, 1)) - lh)*sqrt(gg*(unsth(limesh(li, 1)) - lh))*uy(li)
            else
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        else
            ! 田んぼダム適用場所以外
            if(unsth(limesh(li, 1)) > th) then
                if(baseo(limesh(li,2))==-9999.0d0 .and. baseo(limesh(li,1))/=-9999.0d0) then
                    !自由流出
                    um(li) = unsth(limesh(li,1))*sqrt(2.0d0*gg*unsth(limesh(li,1)))*ux(li)
                    vn(li) = unsth(limesh(li,1))*sqrt(2.0d0*gg*unsth(limesh(li,1)))*uy(li)
                else
                    um(li) = 0.544d0*unsth(limesh(li, 1))*sqrt(gg*unsth(limesh(li, 1)))*ux(li)
                    vn(li) = 0.544d0*unsth(limesh(li, 1))*sqrt(gg*unsth(limesh(li, 1)))*uy(li)
                endif
            else
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        endif
        lhan(li) = 1
        cycle
    elseif (hhw < baseo(limesh(li, 2))) then
        if ((inf(limesh(li,1))==71 .and. inf(limesh(li,2))<71) .or. &
            (inf(limesh(li,2))==71 .and. inf(limesh(li,1))<71)) then
            ! 田んぼダム適用 by k.yamamura(inf=71)
            if(unsth(limesh(li,2)) > lh + th) then
                um(li) = -0.544d0*(unsth(limesh(li, 2)) - lh)*sqrt(gg*(unsth(limesh(li, 2)) - lh))*ux(li)
                vn(li) = -0.544d0*(unsth(limesh(li, 2)) - lh)*sqrt(gg*(unsth(limesh(li, 2)) - lh))*uy(li)
            else
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        else
            !　田んぼダム適用以外
            if(unsth(limesh(li,2)) > th) then
                if(baseo(limesh(li,1))==-9999.0d0 .and. baseo(limesh(li,2))/=-9999.0d0) then
                    !自由流出
                    um(li) = -unsth(limesh(li,2))*sqrt(2.0d0*gg*unsth(limesh(li,2)))*ux(li)
                    vn(li) = -unsth(limesh(li,2))*sqrt(2.0d0*gg*unsth(limesh(li,2)))*uy(li)
                else
                    um(li) = -0.544d0*unsth(limesh(li, 2))*sqrt(gg*unsth(limesh(li, 2)))*ux(li)
                    vn(li) = -0.544d0*unsth(limesh(li, 2))*sqrt(gg*unsth(limesh(li, 2)))*uy(li)
                endif
            else
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        endif
        lhan(li) = 1
        cycle

    ! 完全越流の式
    elseif(hhep*hhwp < 0.0d0) then
        if(unsth(limesh(li, 2)) > 0.0d0 .or. unsth(limesh(li, 1)) > 0.0d0) then
            hhan = hhep - hhwp
            sgn = hhan/abs(hhan)
            hh1 =   max(unsth(limesh(li, 2))+baseo(limesh(li, 2)), unsth(limesh(li, 1)) + baseo(limesh(li, 1))) &
                    - max(baseo(limesh(li, 2)), baseo(limesh(li, 1)))
            if ((inf(limesh(li,1))==71 .and. inf(limesh(li,2))<71) .or. &
                (inf(limesh(li,2))==71 .and. inf(limesh(li,1))<71)) then
                ! 田んぼダム適用
                if(hh1 > lh + th) then
                    um(li) = - sgn*0.35d0*(hh1-lh)*sqrt(2.0d0*gg*(hh1-lh))*ux(li)
                    vn(li) = - sgn*0.35d0*(hh1-lh)*sqrt(2.0d0*gg*(hh1-lh))*uy(li)
                else
                    um(li) = 0.0d0
                    vn(li) = 0.0d0
                endif
            else
                ! 田んぼダム適用以外
                um(li) = - sgn*0.35d0*hh1*sqrt(2.0d0*gg*hh1)*ux(li)
                vn(li) = - sgn*0.35d0*hh1*sqrt(2.0d0*gg*hh1)*uy(li)
            endif
        else
            um(li) = 0.0d0
            vn(li) = 0.0d0
        endif
        lhan(li) = 1
        cycle
    endif

    ! when the water surface is continuous
    ! convective terms
    u11 = 0.0d0
    v11 = 0.0d0
    do k = 1, ko(limesh(li, 1))
    if(lhano(melink(limesh(li, 1), k)) == 1) goto 321
    if(melink(limesh(li, 1), k) == li) goto 321
    ! if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 1) .and. limesh(melink(limesh(li, 1), k), 2)==0) goto 321
    if(uu(melink(limesh(li, 1), k))==0.0d0 .and. vv(melink(limesh(li, 1), k))==0.0d0) goto 321
    if(uu(melink(limesh(li, 1), k))*node_dy(limesh(li,1),k) > 0.0d0) then
        me1 = limesh(li, 1)
    else
        if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 1)) me1 = limesh(melink(limesh(li, 1), k), 2)
        if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 2)) me1 = limesh(melink(limesh(li, 1), k), 1)
        ! water level
        if(me1 == 0 .and. dsmesh == 1) then
            if (dsinf(limesh(li,1))== 2 .or. dsinf(limesh(li,1))== 3) me1 = limesh(li, 1)
        endif
    endif

    if(vv(melink(limesh(li, 1), k))*node_dx(limesh(li,1),k) <0.0d0) then
        me2 = limesh(li, 1)
    else
        if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 1)) me2 = limesh(melink(limesh(li, 1), k), 2)
        if(limesh(li, 1) == limesh(melink(limesh(li, 1), k), 2)) me2 = limesh(melink(limesh(li, 1), k), 1)
        ! water level
        if(me2 == 0 .and. dsmesh == 1) then
            if (dsinf(limesh(li,1))== 2 .or. dsinf(limesh(li,1))== 3) me2 = limesh(li, 1)
        endif
    endif
    if (me1==0 .and. me2==0) goto 321
    u11 = u11 + uu(melink(limesh(li, 1), k))*umm(me1)*node_dy(limesh(li,1),k) - vv(melink(limesh(li, 1), k))*umm(me2)*node_dx(limesh(li,1),k)
    v11 = v11 + uu(melink(limesh(li, 1), k))*vnm(me1)*node_dy(limesh(li,1),k) - vv(melink(limesh(li, 1), k))*vnm(me2)*node_dx(limesh(li,1),k)
321 enddo

    do k = 1, ko(limesh(li, 2))
    if(lhano(melink(limesh(li, 2), k)) == 1) goto 322
    if(melink(limesh(li, 2), k) == li) goto 322
    ! if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 1) .and. limesh(melink(limesh(li, 2), k), 2)==0) goto 322
    if(uu(melink(limesh(li, 2), k))==0.0d0 .and. vv(melink(limesh(li, 2), k))==0.0d0) goto 322
    if(uu(melink(limesh(li, 2), k))*node_dy(limesh(li,2),k) > 0.0d0) then
        me1 = limesh(li, 2)
    else
        if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 1)) me1 = limesh(melink(limesh(li, 2), k), 2)
        if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 2)) me1 = limesh(melink(limesh(li, 2), k), 1)
        ! water level
        if(me1 == 0 .and. dsmesh == 1) then
            if (dsinf(limesh(li,2)) == 2 .or. dsinf(limesh(li,2)) == 3) me1 = limesh(li, 2)
        endif
    endif

    if(vv(melink(limesh(li, 2), k))*node_dx(limesh(li,2),k) < 0.0d0) then
        me2 = limesh(li, 2)
    else
        if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 1)) me2 = limesh(melink(limesh(li, 2), k), 2)
        if(limesh(li, 2) == limesh(melink(limesh(li, 2), k), 2)) me2 = limesh(melink(limesh(li, 2), k), 1)
        ! water level
        if(me2 == 0 .and. dsmesh == 1) then
            if (dsinf(limesh(li,2)) == 2 .or. dsinf(limesh(li,2)) == 3) me2 = limesh(li, 2)
        endif
    endif
    if (me1==0 .and. me2==0) goto 322
    u11 = u11 + uu(melink(limesh(li, 2), k))*umm(me1)*node_dy(limesh(li,2),k) - vv(melink(limesh(li, 2), k))*umm(me2)*node_dx(limesh(li,2),k)
    v11 = v11 + uu(melink(limesh(li, 2), k))*vnm(me1)*node_dy(limesh(li,2),k) - vv(melink(limesh(li, 2), k))*vnm(me2)*node_dx(limesh(li,2),k)
322  enddo
    u11 = u11*dt2/scv(li)
    v11 = v11*dt2/scv(li)

    ! gravitational term
    if(dl(li) /= 0.0d0) then
    u13 = gg*hl(li)*dt2/dl(li)*ux(li)*(unsth(limesh(li, 2)) + baseo(limesh(li, 2)) - unsth(limesh(li, 1)) - baseo(limesh(li, 1)))
    v13 = gg*hl(li)*dt2/dl(li)*uy(li)*(unsth(limesh(li, 2)) + baseo(limesh(li, 2)) - unsth(limesh(li, 1)) - baseo(limesh(li, 1)))
    else
    u13 = 0.0d0
    v13 = 0.0d0
    endif

    ! shear term
    sqx = sqrt(uu(li)**2 + vv(li)**2)
    ram = gg*rnx(li)**2*sqx/hl(li)**1.333333
    if(hl(li) <= th) ram = 0.0d0

    ! added by k.yamamura
    ! defended forest term
    if(plantDa==1) then
        if(dk_val(li)>0.0d0) then
            plant_force = gg*(sqx**2.0d0)/(dk_val(li)**2.0d0)
            ram = ram + plant_force
        endif
    endif

    ! um, vn calculation
    um(li) = ((1.0d0 - dt2*ram*(1.0d0 - fita))*umo(li) - u11 - u13)/(1.0d0 + dt2*ram*fita)
    vn(li) = ((1.0d0 - dt2*ram*(1.0d0 - fita))*vno(li) - v11 - v13)/(1.0d0 + dt2*ram*fita)
    lhan(li) = 0
    cycle


300 um(li) = 0.0d0
    vn(li) = 0.0d0
    lhan(li) = 1
! 301 if(um(li) >  10.0d0) um(li) =  10.0d0
!     if(um(li) < -10.0d0) um(li) = -10.0d0
!     if(vn(li) >  10.0d0) vn(li) =  10.0d0
!     if(vn(li) < -10.0d0) vn(li) = -10.0d0
    enddo
    !$omp end parallel do

    return

    ! equation of continuity
    entry suisin

    allocate(rr(mesh), q(mesh), gamma(link))
    if (paddydam==1) paddy_q = 0.0d0

    it = int(unsttime/dtrain) + 1

    if(drainarea==1) vol = 0.0d0
    if(drainarea==1) vol_dr = 0.0d0
    ct = int(unsttime / (dt2+0.000001d0))
    nt = int((unsttime + (dt2+0.000001d0)) / (dt2+0.000001d0))
    !$omp parallel do default(shared),private(li)
    do li = 1, link
        gamma(li) = 1.0d0
    enddo
    !$omp end parallel do
    iflg0 = 0

400 iflg0 = iflg0 + 1
    iflg = 0
    dvp1 = 0.0d0
    dvp2 = 0.0d0
    dvc = 0.0d0
    dv_cextre = 0.0d0

    !$omp parallel do default(shared) &
    !$omp private(f,sumf,dv,rem,rt,k,me,sadj,ilt,nn) &
    !$omp reduction(+:dvp1, dvp2, dvc, dv_cextre)
    do me = 1, mesh
        if(inf(me) == 0) cycle
        if(dsmesh==1) then
            ! water level(dstype:2) added by d.baba
            if (dsinf(me)==2) then
                ii = int( unsttime / dsdt(me)) + 1
                temp_dsdepth = dsdepth(me, ii) + (unsttime - dsdt(me)*dble(ii - 1))/dsdt(me)*(dsdepth(me, ii + 1) - dsdepth(me, ii))
                unsth(me) = max(temp_dsdepth - baseo(me), 0.000d0)
                cycle
            ! uniform flow(dstype:3) added by d.baba
            elseif (dsinf(me)==3) then
                do k = 1, ko(me)
                    if (limesh(melink(me, k),1) /= dsupper(me) .and. limesh(melink(me, k),2) /= dsupper(me)) then
                        um(melink(me, k)) = um(melink(dsupper(me), k))
                        vn(melink(me, k)) = vn(melink(dsupper(me), k))
                    endif
                enddo
                unsth(me) = unsth(dsupper(me))
                cycle
            endif
        endif
        ! added by k.yamamura
        rr(me) = rain(it, urain_i(me), urain_j(me))

        ! inflow by flux
        sumf = 0.0d0
        sadj = 0.0d0
        do k = 1, ko(me)
            umbeta(melink(me, k)) = um(melink(me, k))*rbeta(melink(me, k))
            vnbeta(melink(me, k)) = vn(melink(me, k))*rbeta(melink(me, k))
            f = umbeta(melink(me, k))*node_dy(me,k) - vnbeta(melink(me, k))*node_dx(me,k)
            sumf = sumf + f
            if(str_type==1) then
            if(f<0 .and. limesh(melink(me,k),1)==me) sadj = sadj + co(limesh(melink(me,k),2))*f
            if(f<0 .and. limesh(melink(me,k),2)==me) sadj = sadj + co(limesh(melink(me,k),1))*f
            if(f>=0) sadj = sadj + co(me)*f
            endif
        enddo

        ! calculate depth
        q(me) = dt2*(-sumf/(1.0d0 - lambda(me)) + rr(me)*rnof(me)*smesh(me))
        rem = strmx(me) - stro(me)
        if(rem <= 0.0d0 .or. strmx(me) <= 0.0d0) goto 403
        str(me) = stro(me) + min(q(me),rem)
        dvc = dvc - unstc(me)*min(q(me),rem)
        q(me) = q(me) - min(q(me),rem)

        ! 水深の更新
    403 unsth(me) = ho(me) + q(me)/smesh(me)
        if(paddydam==1 .and. inf(me) == 71) then
            if(unsth(me) >= etp(paddyid(me))) then
                unsth(me) = unsth(me) - etp(paddyid(me))
                q(me) = q(me) - (etp(paddyid(me)) * smesh(me))
            endif
        endif
        if(baseo(me)==-9999.0d0) unsth(me) = 0.0d0

        ! 浸透モデル
        if( unst_soildepth(me) .gt. 0.d0 .and. unst_ksv(me) .gt. 0.d0 ) call unst_infilt(me)

        if(str_type==0) then
            unsth(me) = max(unsth(me), 0.0d0)
        elseif(str_type==1) then
            unstc(me) = ((co(me)*ho(me)+cr(me)*rr(me)*rnof(me)*dt2)*smesh(me) - sadj/(1.0d0 - lambda(me))*dt2)/(unsth(me)*smesh(me))
            if(unsth(me)<=1.0d-4) goto 500
        endif
        goto 501

        !負の水深について
    500 dv_cextre = dv_cextre + (co(me)*ho(me)+cr(me)*rr(me)*rnof(me)*dt2)*smesh(me) - sadj/(1.0d0 - lambda(me))*dt2
        unstc(me) = 0.0d0
        if(unsth(me)<0.0d0) then
            iflg = 1
            rt = (ho(me) - th)/(ho(me) - unsth(me))
            do k = 1, ko(me)
                gamma(melink(me, k)) = min(gamma(melink(me, k)), rt)
            enddo
        endif

        !　calculate sewerage and fields depth
        ! added by k.yamamura
    501 if(drainarea==1) then
        if(inf_dr(me)>0) then
                if(unsth(me)*smesh(me)<=vol_dr(me)) then
                    vol(me) = 0.0d0
                else
                    vol(me) = vol_dr(me)
                    unsth(me) = (unsth(me)*smesh(me) - vol(me)) / smesh(me)
                endif
                !Time-delayed discharge to the outlet point
                !$omp critical
                ilt = tripTime(me)
                if(ilt>0) then
                    do nn = 1, dr_no
                        dhj(ilt,nn) = dhj(ilt,nn) + vol(me)
                    enddo
                endif
                !$omp end critical
        endif
        endif

        ! calculate paddy field flow
        if(paddydam==1) then
            if(paddyid(me)>0) then
                !$omp atomic
                paddy_q(paddyid(me)) = paddy_q(paddyid(me)) + (unsth(me) * smesh(me))
            endif
        endif

        ! check water balance
        if(paddydam==1) then
            if(inf(me) /= 71) qr_sum(me) = qr_sum(me) + q(me)
        elseif(drainarea==1) then
            if(inf_dr(me)>0) qr_sum(me) = qr_sum(me) + (q(me) - vol(me))
            if(inf_dr(me)==0) qr_sum(me) = qr_sum(me) + q(me)
        else
            if(baseo(me)/=-9999.0d0) qr_sum(me) = qr_sum(me) + q(me)
        endif
    enddo
    !$omp end parallel do

    if(iflg0 <= 2 .and. iflg == 1) goto 400
    svc = svc + dvc
    v_cextre = v_cextre +dv_cextre

    if(str_type == 1) then
    !$omp parallel do default(shared),private(me) &
    !$omp reduction(+:v_minus_all, vrain, svc, v_cminus, v_cplus)
    do me = 1, mesh
        if(unsth(me)<0.0d0) then
            v_minus_all = v_minus_all + unsth(me)*smesh(me)
            !$omp critical
            v_minus(inf(me)) = v_minus(inf(me)) + unsth(me)*smesh(me)
            !$omp end critical
            v_cminus = v_cminus + unstc(me)*unsth(me)*smesh(me)
            unstc(me) = 0.0d0
        endif
        unsth(me) = max(unsth(me), 0.0d0)
        if(unstc(me)<0.0d0) v_cminus = v_cminus + unstc(me)*unsth(me)*smesh(me)
        unstc(me) = max(unstc(me), 0.0d0)
        if(unstc(me)>1.0d0) v_cplus = v_cplus + (unstc(me)-1.0d0)*unsth(me)*smesh(me)
        unstc(me) = min(unstc(me), 1.0d0)
        vrain = vrain + rr(me)*rnof(me)*smesh(me)*dt2
        svc = svc + cr(me)*rr(me)*rnof(me)*smesh(me)*dt2
    enddo
    !$omp end parallel do
    endif

    ! calculate paddy field dam outflow
    if (paddydam==1) call paddyflow
    if (ct /= nt .and. paddydam==1) call timedelay_paddyflow

    ! calculate sewerage and fields outflow
    if (ct /= nt .and. drainarea==1)  call drainflow

    qout = qout + dvp1 + dvp2

    deallocate(rr, q, gamma)
end subroutine flux

!----------------------------
! Calculate velocity
!----------------------------
subroutine velocity
    use unst_globals_mod
    implicit none
    integer li, me, k, idsf2

    ! hl calculation
    !$omp parallel do default(shared),private(li)
    do li = 1, link
        if(limesh(li, 2) == 0) then
        hl(li) = unsth(limesh(li, 1))
        else
        hl(li) = unsth(limesh(li, 1))*rthl(li, 1) + unsth(limesh(li, 2))*rthl(li, 2)
        endif
    enddo
    !$omp end parallel do

    ! umm, vnm calculation
    !$omp parallel do default(shared),private(me, k)
    do me = 1, mesh
        umm(me) = 0.0d0
        vnm(me) = 0.0d0
    do k = 1, ko(me)
        umm(me) = umm(me) + um(melink(me, k))*rtuv_x(me, k)
        vnm(me) = vnm(me) + vn(melink(me, k))*rtuv_y(me, k)
    enddo
    enddo
    !$omp end parallel do

    if(dsmesh==1) then
    !$omp parallel do default(shared),private(idsf2, k)
    do idsf2 = 1, idsfilter2
        do k = 1, ko(dsfilter2(idsf2))
            if(limesh(melink(dsfilter2(idsf2), k),2)==0) then
                um(melink(dsfilter2(idsf2), k)) = umm(dsfilter2(idsf2))
                vn(melink(dsfilter2(idsf2), k)) = vnm(dsfilter2(idsf2))
                if(unsth(dsfilter2(idsf2)) < th) then
                if((umm(dsfilter2(idsf2))*node_dy(dsfilter2(idsf2),k) - &
                     vnm(dsfilter2(idsf2))*node_dx(dsfilter2(idsf2),k)) > 0.0d0) then
                        um(melink(dsfilter2(idsf2), k)) = 0.0d0
                        vn(melink(dsfilter2(idsf2), k)) = 0.0d0
                endif
                endif
            endif
        enddo
    enddo
    !$omp end parallel do
    endif

    ! uu, vv calculation
    !$omp parallel do default(shared),private(me, k)
    do li = 1, link
        if(hl(li) > th) then
        uu(li) = um(li)/hl(li)
        vv(li) = vn(li)/hl(li)
        else
        uu(li) = 0.0d0
        vv(li) = 0.0d0
        endif
    enddo
    !$omp end parallel do

    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        if(unsth(me) < th) then
        uum(me) = 0.0d0
        vvm(me) = 0.0d0
        else
        uum(me) = umm(me)/unsth(me)
        vvm(me) = vnm(me)/unsth(me)
        endif
    enddo
    !$omp end parallel do

end subroutine velocity

!-------------------------
! Replace (Update) data
!-------------------------
subroutine replace
    use unst_globals_mod
    implicit none
    integer li, me

    !$omp parallel do default(shared),private(li)
    do li = 1, link
        umo(li) = um(li)
        vno(li) = vn(li)
        lhano(li) = lhan(li)
    enddo
    !$omp end parallel do

    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        ho(me) = unsth(me)
    !   stro(me) = str(me)
    !   co(me) = c(me)
    enddo
    !$omp end parallel do

end subroutine replace

!--------------------------------------
! Calculate inflow on the boundaries
!--------------------------------------
subroutine lkyokai
    use unst_globals_mod
    implicit none
    real(8) yqin, xqin, rqin
    integer j, ii

    ii = int(unsttime/dtq) + 1
    !$omp parallel do default(shared),private(j, yqin, xqin, rqin)
    do j = 1, iqnum
        yqin = qinv(j, ii) + (unsttime - dtq*dble(ii - 1))/dtq*(qinv(j, ii + 1) - qinv(j, ii))
        xqin = qinu(j, ii) + (unsttime - dtq*dble(ii - 1))/dtq*(qinu(j, ii + 1) - qinu(j, ii))
        rqin = qin(j, ii) + (unsttime - dtq*dble(ii - 1))/dtq*(qin(j, ii + 1) - qin(j, ii))

        if(lkyokai_dx(j) - lkyokai_dy(j) >= 0) then
            yqin = yqin - rqin
        else
            xqin = xqin - rqin
        endif

        if(lkyokai_dx(j) > 0) vn(inl(j)) = yqin/lkyokai_dx(j)
        if(lkyokai_dy(j) > 0) um(inl(j)) = xqin/lkyokai_dy(j)
    enddo
    !$omp end parallel do

end subroutine lkyokai

!---------------------------------------------
! Aggregate inundation statistics
!---------------------------------------------
subroutine sumqa(sv, sa, saj)
    use unst_globals_mod
    implicit none
    real(8) sv, sa, saj
    integer me

    ! inundation water level
    sv = 0.0d0
    sa = 0.0d0
    saj = 0.0d0
    !$omp parallel do default(shared),private(me),reduction(+:sv, sa, saj)
    do me = 1, mesh
        if(inf(me) == 0) goto 100
        sv = sv + unsth(me)*smesh(me)*(1.0d0 - lambda(me))
        if(unsth(me) > 1.0d-3) sa  = sa  + smesh(me)
        if(unsth(me) > 1.0d-3) saj = saj + smesh(me)*(1.0d0 - lambda(me))  ! flooded area
100 enddo
    !$omp end parallel do

end subroutine sumqa

!----------------------------------
! Calculate infiltration process
!----------------------------------
subroutine unst_infilt(me)
    use unst_globals_mod
    implicit none

    real(8) unst_gampt_ff_temp
    integer me

    unst_gampt_f(me) = 0.d0
    unst_gampt_ff_temp = unst_gampt_ff(me)
    if( unst_gampt_ff_temp .le. 0.01d0 ) unst_gampt_ff_temp = 0.01d0

    ! unst_gampt_f(me) : infiltration capacity [m/s]
    ! unst_gampt_ff : accumulated infiltration depth [m]
    unst_gampt_f(me) = unst_ksv(me) * (1.d0 + unst_faif(me) * unst_gammaa(me) / unst_gampt_ff_temp)

    ! unst_gampt_f(me) : infiltration capacity -> infiltration rate [m/s]
    if( unst_gampt_f(me) .ge. unsth(me) / unstdt ) unst_gampt_f(me) = unsth(me) / unstdt

    ! unst_gampt_ff should not exceeds a certain level
    if( unst_infilt_limit(me) .ge. 0.d0 .and. unst_gampt_ff(me) .ge. unst_infilt_limit(me) ) unst_gampt_f(me) = 0.d0

    ! update unst_gampt_ff [m]
    unst_gampt_ff(me) = unst_gampt_ff(me) + unst_gampt_f(me) * unstdt

    ! hs : hs - infiltration rate * dt [m]
    unsth(me) = unsth(me) - unst_gampt_f(me) * unstdt
    if( unsth(me) .le. 0.d0 ) unsth(me) = 0.d0

end subroutine unst_infilt
