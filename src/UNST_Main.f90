!***********************************************************************
! UNST_Main.f90
! Coupled model with RRI model
! Coded by K.Kawaike, Added by TK Labo
! Released on July 7th 2025
!***********************************************************************

!-----------------------------------------------------------------------
! UNSTメインサブルーチン
! 入力
!   ny: RRIの y 軸方向の格子数
!   nx: RRIの x 軸方向の格子数
!   domain: RRIの領域設定
!   riv: RRIの河川の設定
!   area: RRIの単一の格子の面積
!   time: RRIの時間
! 入力/出力:
!   hs: RRIの斜面水深（更新される）
!   hr: RRIの河川水深（更新される）
!-----------------------------------------------------------------------
subroutine UNST(ny, nx, domain, riv, area, time, hs, hr)
    use unst_globals_mod
    implicit none
    integer, intent(in) :: ny, nx
    integer, intent(in) :: domain(ny, nx), riv(ny, nx)
    real(8), intent(in) :: area, time
    real(8), intent(inout) :: hs(ny,nx), hr(ny,nx)
    integer me, li, k

    ! 計算開始メッセージの出力
    write(*,*) 'UNST timestep ===  ',int(time-timmax), '  >>>>>  ', int(time)

    ! UNSTの時間設定
    unsttime = time - timmax
    mstep = int(unsttime / unstdt)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                  時間計算ループ開始
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !  運動方程式の計算（フラックス計算）
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  1 call flux   ! セル間のフラックスを計算
    call lkyokai ! 境界条件の適用

    ! 先端処理（水深が閾値以下のメッシュにおける流速の設定）
    !$omp parallel do default(shared),private(me,li,k)
    do me = 1, mesh
        if(unsth(me) >= th) goto 200
        do k = 1, ko(me)
            li = melink(me, k)
            if((um(li)*node_dy(me,k) - vn(li)*node_dx(me,k)) > 0.0d0) then
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        enddo
200 enddo
    !$omp end parallel do

    ! 時間ステップの更新
    unsttime = unsttime + unstdt
    mstep = mstep + 1

    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !  連続式の計算（水位計算）
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    call suisin  ! 水位の更新計算

    ! 流速計算
    call velocity

    ! 最大水深・流速の更新
    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        hmax(me) = max(hmax(me), unsth(me))
        uummax(me) = max(uummax(me), abs(uum(me)))
        vvmmax(me) = max(vvmmax(me), abs(vvm(me)))
    enddo
    !$omp end parallel do

    ! データの更新（新データ >>> 旧データ）
    call replace

    ! 時間ステップの更新
    unsttime = unsttime + unstdt
    mstep = mstep + 1

    ! 結果出力
    if(mod(mstep, lkout) == 0) call diskwrite   ! ディスク出力
    if(mod(mstep, lpout) == 0) call dispwrite   ! 画面出力
    if(paddydam==1 .and. mod(mstep, lkout) == 0) call paddywrite  ! 田んぼダム出力

    ! 時間判定（次の計算を行うかどうか）
    if(unsttime + unstdt <= time) goto 1
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                  時間計算ループ終了
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! UNSTからRRIへのデータ転送
    call unst2rri(ny, nx, domain, riv, area, hs, hr)

end subroutine UNST
