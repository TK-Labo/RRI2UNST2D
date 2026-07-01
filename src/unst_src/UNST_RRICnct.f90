! Subroutines for data exchange between RRI and UNST models
! Coded by TK Labo
! Released on July 7th 2025

!-----------------------------------------------------------------------
! Transfer of flow data from RRI to UNST
! In:
!   ny: RRIの y 軸方向の格子数/RRI num of rows
!   nx: RRIの x 軸方向の格子数/RRI num of cols
!   i4: RRIの流量配列の数 (4)/RRI num of q-array
!   qr_ave: RRIの河川流量/RRI river-q
!   qs_ave: RRIの斜面流量/RRI slope-q
!   hr: RRIの河川水深/RRI river depth(h)
!   hs: RRIの斜面水深/RRI slope depth(h)
!   area: RRIの単一の格子の面積/RRI cell area
!   time: RRIの時間/RRI cal time
! In/Out:
!   uflg: UNSTの実行フラグ/UNST2D RUN FLAG
!-----------------------------------------------------------------------
subroutine unst_qin(ny, nx, i4, &
     qr_ave, qs_ave, hr, hs, area, time, uflg)
    use unst_globals_mod
    use d1riv_globals_mod
    implicit none
    integer, intent(in) :: ny, nx, i4
    real(8), intent(in) :: qr_ave(ny,nx), qs_ave(i4,ny,nx)
    real(8), intent(in) :: hs(ny,nx), hr(ny,nx)
    real(8), intent(in) :: area, time
    integer, intent(inout) :: uflg

    real(8) qu(iqnum), qv(iqnum), unst_qr(iqnum)
    integer i, ii, id, n
    integer uflg_1d

    ! 現在の時間ステップのインデックスを計算
    ii = int(time/dtq) + 1

    if(qin_type==0) then
        !$omp parallel do default(shared),private(i, id)
        do i = 1, iqnum
            ! target UNST2D mesh id
            id = limesh(inl(i),1)
            if(lkyokai_dir(i)==0) cycle
            ! x方向の流量計算 (RRI斜面方向からUNST方向への変換)
            ! convert x direction q from RRI q
            qu(i) = (qs_ave(1, rsetsu_i(id), rsetsu_j(id)) &
                + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) &
                    - qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

            ! y方向の流量計算 (RRI斜面方向からUNST方向への変換)
            ! convert y direction q from RRI q
            qv(i) =  (qs_ave(2, rsetsu_i(id), rsetsu_j(id)) &
                + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) &
                    + qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

            ! 流入流量の計算
            ! set qin array
            qin(i, ii) = abs(qr_ave(rsetsu_i(id), rsetsu_j(id)))
            qinu(i, ii) = qu(i)
            qinv(i, ii) = qv(i)
        end do
        !$omp end parallel do
    elseif(qin_type==1) then
        !$omp parallel do default(shared),private(i, id)
        do i = 1, iqnum
            ! target UNST2D mesh id
            id = limesh(inl(i),1)
            if(lkyokai_dir(i)==0) cycle
            ! x方向の流量計算 (RRI斜面方向からUNST方向への変換)
            ! convert x direction q from RRI q
            qu(i) = (qs_ave(1, rsetsu_i(id), rsetsu_j(id)) &
                + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) &
                    - qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

            ! y方向の流量計算 (RRI斜面方向からUNST方向への変換)
            ! convert y direction q from RRI q
            qv(i) =  (qs_ave(2, rsetsu_i(id), rsetsu_j(id)) &
                + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) &
                    + qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

            unst_qr(i) = qr_ave(rsetsu_i(id), rsetsu_j(id))

            select case(qin_inf(id))
            case(1)
                qv(i) = qv(i) + &
                        (qs_ave(2, rsetsu_i(id), rsetsu_j(id)+1) &
                        + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)+1) &
                        + qs_ave(4, rsetsu_i(id), rsetsu_j(id)+1)) / 2.d0) * area
                unst_qr(i) = unst_qr(i) + qr_ave(rsetsu_i(id), rsetsu_j(id)+1) * 0.5d0
            case(2)
                qu(i) = qu(i) + &
                        (qs_ave(1, rsetsu_i(id)+1, rsetsu_j(id)) &
                        + (qs_ave(3, rsetsu_i(id)+1, rsetsu_j(id)) &
                        - qs_ave(4, rsetsu_i(id)+1, rsetsu_j(id))) / 2.d0) * area
                unst_qr(i) = unst_qr(i) + qr_ave(rsetsu_i(id)+1, rsetsu_j(id)) * 0.5d0
            case(3)
                qv(i) = qv(i) + &
                        (qs_ave(2, rsetsu_i(id), rsetsu_j(id)+1) &
                        + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)+1) &
                        + qs_ave(4, rsetsu_i(id), rsetsu_j(id)+1)) / 2.d0) * area
                unst_qr(i) = unst_qr(i) + qr_ave(rsetsu_i(id), rsetsu_j(id)+1) * 0.5d0
            case(4)
                qu(i) = qu(i) + &
                        (qs_ave(1, rsetsu_i(id)-1, rsetsu_j(id)) &
                        + (qs_ave(3, rsetsu_i(id)-1, rsetsu_j(id)) &
                        - qs_ave(4, rsetsu_i(id)-1, rsetsu_j(id))) / 2.d0) * area
                unst_qr(i) = unst_qr(i) + qr_ave(rsetsu_i(id)-1, rsetsu_j(id)) * 0.5d0
            case(5)
                qv(i) = qv(i) + &
                        (qs_ave(2, rsetsu_i(id), rsetsu_j(id)-1) &
                        + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)-1) &
                        + qs_ave(4, rsetsu_i(id), rsetsu_j(id)-1)) / 2.d0) * area
                unst_qr(i) = unst_qr(i) + qr_ave(rsetsu_i(id), rsetsu_j(id)-1) * 0.5d0
            case(6)
                qu(i) = qu(i) + &
                        (qs_ave(1, rsetsu_i(id)-1, rsetsu_j(id)) &
                        + (qs_ave(3, rsetsu_i(id)-1, rsetsu_j(id)) &
                        - qs_ave(4, rsetsu_i(id)-1, rsetsu_j(id))) / 2.d0) * area
                unst_qr(i) = unst_qr(i) + qr_ave(rsetsu_i(id)-1, rsetsu_j(id)) * 0.5d0
            case(7)
                qv(i) = qv(i) + &
                        (qs_ave(2, rsetsu_i(id), rsetsu_j(id)-1) &
                        + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)-1) &
                        + qs_ave(4, rsetsu_i(id), rsetsu_j(id)-1)) / 2.d0) * area
                unst_qr(i) = unst_qr(i) + qr_ave(rsetsu_i(id), rsetsu_j(id)-1) * 0.5d0
            case(8)
                qu(i) = qu(i) + &
                        (qs_ave(1, rsetsu_i(id)+1, rsetsu_j(id)) &
                        + (qs_ave(3, rsetsu_i(id)+1, rsetsu_j(id)) &
                        - qs_ave(4, rsetsu_i(id)+1, rsetsu_j(id))) / 2.d0) * area
                unst_qr(i) = unst_qr(i) + qr_ave(rsetsu_i(id)+1, rsetsu_j(id)) * 0.5d0
            end select

            ! 流入流量の計算
            ! set qin array
            qin(i, ii) = unst_qr(i)  ! update v.1.0.5
            qinu(i, ii) = qu(i)
            qinv(i, ii) = qv(i)
        end do
        !$omp end parallel do
    endif
    ! safety
    qin(:, min(ii+1, iqin+1)) = qin(:, ii)
    qinu(:, min(ii+1, iqin+1)) = qinu(:, ii)
    qinv(:, min(ii+1, iqin+1)) = qinv(:, ii)

    ! 流入判定フラグの設定
    ! 流入がある場合　uflg = 1
    ! 流入がない場合　uflg = 0
    if (any(qin > 0.0d0) .or. any(hr > 0.0d0) .or. any(hs > 0.0d0)) then
        uflg = 1
    else
        uflg = 0
    endif

    ! 結果をファイルに出力
    write(fkyokaiq_unit, '(A8, f8.1)') 'time = ', time
    write(fkyokaiq_unit, '(10f14.5)') (qin(i, ii), i = 1, iqnum)

    if(d1riv==1) then
        uflg_1d = 0
        do n = 1, nbound_1d
            i = b_idx_1d(n)
            if(ntype_1d(i)/=-100) cycle
            b_data_1d(ii,n) = qr_ave(rsetsu_i_1d(i), rsetsu_j_1d(i))
            b_data_1d(min(ii+1,max_nb),n) = b_data_1d(ii,n)
            if(b_data_1d(ii,n)/=0.0d0) uflg_1d = 1
        enddo

        if(uflg==0 .and. uflg_1d==1) then
            uflg = 1
        endif
    endif

    end subroutine unst_qin


!-----------------------------------------------------------------------
! Transfer of depth data from UNST to RRI
! In:
!   ny: RRIの y 軸方向の格子数/RRI num of rows
!   nx: RRIの x 軸方向の格子数/RRI num of cols
!   domain: RRIの領域設定/RRI active domain
!   riv: RRIの河川の設定/RRI river domain
!   area: RRIの単一の格子の面積/RRI cell area
! Out:
!   hs: RRIの斜面水深（更新される）/RRI slope depth(update)
!   hr: RRIの河川水深（更新される）/RRI river depth(update)
!-----------------------------------------------------------------------
subroutine unst2rri(ny, nx, domain, riv, hs, hr)
    use unst_globals_mod
    use d1riv_globals_mod
    implicit none

    integer, intent(in) :: ny, nx
    integer, intent(in) :: domain(ny, nx), riv(ny, nx)
    real(8), intent(inout) :: hs(ny,nx), hr(ny,nx)

    real(8), allocatable :: hsmxdif(:, :), hsmindif(:, :), hrmxdif(:, :), hrmindif(:, :)
    real(8), allocatable :: hrr(:, :)
    integer i, j, me, n
    integer, allocatable :: counthr(:, :), counths(:, :), counthrr(:, :)
    real(8), allocatable :: minbaseo(:, :)
    integer, allocatable :: u_to_r_update(:, :)

    allocate(hsmxdif(ny,nx), hsmindif(ny,nx), hrmxdif(ny,nx), hrmindif(ny,nx))
    allocate(hrr(ny,nx), counthr(ny,nx), counths(ny,nx), counthrr(ny,nx), minbaseo(ny, nx))
    allocate(u_to_r_update(ny,nx))

    ! 配列の初期化
    counthr = 0
    counths = 0
    counthrr = 0
    hsmxdif = 0.0d0
    hsmindif = 0.0d0
    hrmxdif = 0.0d0
    hrmindif = 0.0d0
    hrr = 0.0d0
    minbaseo = 99999999.0d0
    u_to_r_update = 0

    ! UNST各メッシュからRRIメッシュへの水深データ転送（集計）
    do me = 1, mesh
        i = rsetsu_i(me)
        j = rsetsu_j(me)

        ! 計算対象外セルはスキップ
        if(i==0 .or. j==0) cycle
        if (domain(i,j) == 0) cycle
        if(baseo(me) == -9999.0d0) cycle

        if (riv(i,j) /= 0) then
            ! 河川セルの場合
            if(counthr(i,j)==0) hr(i,j) = 0.0d0

            if(inf(me) == 24) then
                ! 地盤高が最も低いメッシュを選んで水深を設定
                if(baseo(me) < minbaseo(i,j)) then
                    minbaseo(i,j) = baseo(me)
                    hr(i,j) = unsth(me)
                    counthr(i,j) = 1
                    minbaseo(i,j) = min(minbaseo(i,j), baseo(me))
                elseif(baseo(me) == minbaseo(i,j)) then
                    counthr(i,j) = counthr(i,j) + 1
                    hr(i,j) = hr(i,j) + unsth(me)
                endif
                hrmxdif(i,j) = max(hrmxdif(i,j), unsth(me))
                hrmindif(i,j) = min(hrmindif(i,j), unsth(me))
            else
                counthrr(i,j) = counthrr(i,j) + 1
                hrr(i,j) = hrr(i,j) + unsth(me)
            endif
        else
            ! 斜面セルの場合
            if(counths(i,j)==0) hs(i,j) = 0.0d0
            counths(i,j) = counths(i,j) + 1
            hs(i,j) = hs(i,j) + unsth(me)
            hsmxdif(i,j) = max(hsmxdif(i,j), unsth(me))
            hsmindif(i,j) = min(hsmindif(i,j), unsth(me))
        endif
    enddo

    ! 平均水深の計算
    do me = 1, mesh
        i = rsetsu_i(me)
        j = rsetsu_j(me)
        if(i==0 .or. j==0) cycle
        if (domain(i,j) == 0) cycle
        if(baseo(me) == -9999.0d0) cycle
        if(u_to_r_update(i,j) == 1) cycle  ! fixed v.1.0.5

        if (riv(i,j) /= 0) then
            ! 河川セルの場合の平均化
            if(counthr(i,j) > 0) then
                hr(i,j) = hr(i,j) / counthr(i,j)
                u_to_r_update(i,j) = 1  ! fixed v.1.0.5
            elseif(counthr(i,j) == 0) then
                hr(i,j) = hrr(i,j) / counthrr(i,j)
                u_to_r_update(i,j) = 1  ! fixed v.1.0.5
            endif
        else
            ! 斜面セルの場合の平均化
            hs(i,j) = hs(i,j) / counths(i,j)
        endif
    enddo

    if(d1riv==1) then
        if(rsetsu_1d==1) then
            do n = 1, ndan
                i = rsetsu_i_1d(n)
                j = rsetsu_j_1d(n)
                if(i==0 .or. j==0) cycle
                if (domain(i,j) == 0) cycle
                if (riv(i,j) /= 0) hr(i,j) = h_1d(n)
            enddo
        endif
    endif

    deallocate(hsmxdif, hsmindif, hrmxdif, hrmindif)
    deallocate(hrr, counthr, counths, counthrr, minbaseo, u_to_r_update)

end subroutine unst2rri
