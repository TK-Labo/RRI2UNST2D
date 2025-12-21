! Subroutines for data exchange between RRI and UNST models
! Coded by TK Labo
! Released on July 7th 2025

!-----------------------------------------------------------------------
! Transfer of flow data from RRI to UNST
! In:
!   ny: RRIの y 軸方向の格子数
!   nx: RRIの x 軸方向の格子数
!   i4: RRIの流量配列の数 (4)
!   qr_ave: RRIの河川流量
!   qs_ave: RRIの斜面流量
!   hr: RRIの河川水深
!   hs: RRIの斜面水深
!   area: RRIの単一の格子の面積
!   time: RRIの時間
! In/Out:
!   uflg: UNSTの実行フラグ
!-----------------------------------------------------------------------
subroutine unst_qin(ny, nx, i4, &
     qr_ave, qs_ave, hr, hs, area, time, uflg)
use unst_globals_mod
implicit none

! 外部から受け取る配列
integer, intent(in) :: ny, nx, i4
real(8), intent(in) :: qr_ave(ny,nx), qs_ave(i4,ny,nx)
real(8), intent(in) :: hs(ny,nx), hr(ny,nx)
real(8), intent(in) :: area, time
integer, intent(inout) :: uflg

real(8) qu(iqnum), qv(iqnum), hr_and_hs(iqnum)
integer i, ii, id

! 現在の時間ステップのインデックスを計算
ii = int(time/dtq) + 1

if(qin_type==0) then
    !$omp parallel do default(shared),private(i, id)
    do i = 1, iqnum
        ! 対応するRRIメッシュのIDを取得
        id = limesh(inl(i),1)

        ! x方向の流量計算 (RRI斜面方向からUNST方向への変換)
        qu(i) = (qs_ave(1, rsetsu_i(id), rsetsu_j(id)) &
            + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) &
                - qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

        ! y方向の流量計算 (RRI斜面方向からUNST方向への変換)
        qv(i) =  (qs_ave(2, rsetsu_i(id), rsetsu_j(id)) &
            + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) &
                + qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

        ! 流入流量の計算
        qin(i, ii) = abs(qr_ave(rsetsu_i(id), rsetsu_j(id)))
        qinu(i, ii) = qu(i)
        qinv(i, ii) = qv(i)
    end do
    !$omp end parallel do
elseif(qin_type==1) then
    !$omp parallel do default(shared),private(i, id)
    do i = 1, iqnum
        ! 対応するRRIメッシュのIDを取得
        id = limesh(inl(i),1)

        ! x方向の流量計算 (RRI斜面方向からUNST方向への変換)
        qu(i) = (qs_ave(1, rsetsu_i(id), rsetsu_j(id)) &
            + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) &
                - qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

        ! y方向の流量計算 (RRI斜面方向からUNST方向への変換)
        qv(i) =  (qs_ave(2, rsetsu_i(id), rsetsu_j(id)) &
            + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) &
                + qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

        ! Cal Volume(m3)  v.1.0.5
        hr_and_hs(i) = hs(rsetsu_i(id), rsetsu_j(id))*area &
            + hr(rsetsu_i(id), rsetsu_j(id))*vin_coef(i)
        
        ! 流入流量の計算
        qin(i, ii) = qr_ave(rsetsu_i(id), rsetsu_j(id))  ! update v.1.0.5
        qinu(i, ii) = qu(i)
        qinv(i, ii) = qv(i)
        vin(i, ii) = hr_and_hs(i)  ! update v.1.0.5
    end do
    !$omp end parallel do
endif

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
if(qin_type==0) then
    write(fkyokaiq_unit, '(10f14.5)') (qin(i, ii), i = 1, iqnum)
elseif(qin_type==1) then
    write(fkyokaiq_unit, '(10f14.5)') (vin(i, ii), i = 1, iqnum)
endif

end subroutine unst_qin


!-----------------------------------------------------------------------
! Transfer of depth data from UNST to RRI
! In:
!   ny: RRIの y 軸方向の格子数
!   nx: RRIの x 軸方向の格子数
!   domain: RRIの領域設定
!   riv: RRIの河川の設定
!   area: RRIの単一の格子の面積
! Out:
!   hs: RRIの斜面水深（更新される）
!   hr: RRIの河川水深（更新される）
!-----------------------------------------------------------------------
subroutine unst2rri(ny, nx, domain, riv, hs, hr)
    use unst_globals_mod
    implicit none

    integer, intent(in) :: ny, nx
    integer, intent(in) :: domain(ny, nx), riv(ny, nx)
    real(8), intent(inout) :: hs(ny,nx), hr(ny,nx)

    real(8), allocatable :: hsmxdif(:, :), hsmindif(:, :), hrmxdif(:, :), hrmindif(:, :)
    real(8), allocatable :: hrr(:, :)
    integer i, j, me
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
        if (domain(i,j) == 0) cycle
        if(baseo(me) == -9999.0d0) cycle
        if(u_to_r_update(i,j) == 1) cycle

        if (riv(i,j) /= 0) then
            ! 河川セルの場合の平均化
            if(counthr(i,j) > 0) then
                hr(i,j) = hr(i,j) / counthr(i,j)
                u_to_r_update(i,j) = 1
            elseif(counthr(i,j) == 0) then
                hr(i,j) = hrr(i,j) / counthrr(i,j)
                u_to_r_update(i,j) = 1
            endif
        else
            ! 斜面セルの場合の平均化
            hs(i,j) = hs(i,j) / counths(i,j)
        endif
    enddo

    deallocate(hsmxdif, hsmindif, hrmxdif, hrmindif)
    deallocate(hrr, counthr, counths, counthrr, minbaseo, u_to_r_update)

end subroutine unst2rri
