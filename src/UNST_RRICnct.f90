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

real(8) qu(iqnum), qv(iqnum)
integer i, ii, id

! 現在の時間ステップのインデックスを計算
ii = int(time/dtq) + 1

!$omp parallel do default(shared),private(i, id)
do i = 1, iqnum
    ! 対応するRRIメッシュのIDを取得
    id = limesh(inl(i),1)

    ! x方向の流量計算 (RRI斜面方向からUNST方向への変換)
    qu(i) = (qs_ave(1, rsetsu_i(id), rsetsu_j(id)) + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) - qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

    ! y方向の流量計算 (RRI斜面方向からUNST方向への変換)
    qv(i) =  (qs_ave(2, rsetsu_i(id), rsetsu_j(id)) + (qs_ave(3, rsetsu_i(id), rsetsu_j(id)) + qs_ave(4, rsetsu_i(id), rsetsu_j(id))) / 2.d0) * area

    ! 流入流量の計算
    qin(i, ii) = abs(qr_ave(rsetsu_i(id), rsetsu_j(id)))
    qinu(i, ii) = qu(i)
    qinv(i, ii) = qv(i)

    ! 流入判定フラグの設定
    ! 流入がある場合　uflg = 1
    ! 流入がない場合　uflg = 0
    !$omp critical
    if(uflg == 0) then
        if(any(qin > 0.0d0)) uflg = 1
        if(any(hr > 0.0d0)) uflg = 1
        if(any(hs > 0.0d0)) uflg = 1
    end if
    !$omp end critical

end do
!$omp end parallel do

! 結果をファイルに出力
write(100103, '(A8, f8.1)') 'time = ', time
write(100103, '(10f14.5)') (qin(i, ii), i = 1, iqnum)

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
subroutine unst2rri(ny, nx, domain, riv, area, hs, hr)
    use unst_globals_mod
    implicit none

    integer, intent(in) :: ny, nx
    integer, intent(in) :: domain(ny, nx), riv(ny, nx)
    real(8), intent(in) :: area
    real(8), intent(inout) :: hs(ny,nx), hr(ny,nx)

    real(8) hsmxdif(ny,nx), hsmindif(ny,nx), hrmxdif(ny,nx), hrmindif(ny,nx)
    real(8) hrr(ny,nx)
    integer i, j, me
    integer counthr(ny,nx), counths(ny,nx), counthrr(ny,nx)
    real(8) minbaseo(ny, nx)

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

        if (riv(i,j) /= 0) then
            ! 河川セルの場合の平均化
            if(counthr(i,j) > 0) then
                hr(i,j) = hr(i,j) / counthr(i,j)
            elseif(counthr(i,j) == 0) then
                hr(i,j) = hrr(i,j) / counthrr(i,j)
            endif
        else
            ! 斜面セルの場合の平均化
            hs(i,j) = hs(i,j) / counths(i,j)
        endif
    enddo

    ! 河川と斜面の相互作用を計算
    call funcrs(hr, hs)
    write(*,*) 'UNST >>> RRI h replaced'

end subroutine unst2rri
