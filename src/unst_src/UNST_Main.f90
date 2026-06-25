!***********************************************************************
! UNST_Main.f90
! Coupled model with RRI model
! Coded by K.Kawaike, Added by TK Labo
! Released on July 7th 2025
!***********************************************************************

!-----------------------------------------------------------------------
! UNSTメインサブルーチン/UNST2D main subroutine
! 入力/Input:
!   ny: RRIの y 軸方向の格子数/num of meshes y axis direction
!   nx: RRIの x 軸方向の格子数/num of meshes x axis direction 
!   domain: RRIの領域設定/RRI active mesh switch
!   riv: RRIの河川の設定/RRI river mesh switch
!   time: RRIの時間/RRI time
! 入力-出力/input-output:
!   hs: RRIの斜面水深（更新される）/slope water depth(updated)
!   hr: RRIの河川水深（更新される）/river water depth(updated)
!-----------------------------------------------------------------------
subroutine UNST2D(ny, nx, domain, riv, time, hs, hr)
    use unst_globals_mod
    use unst_cal_sub
    use d1riv_globals_mod
    use unst_d1_main
    use unst_cnct_1d2d
    use unst_wrfile
    implicit none
    interface
        subroutine unst2rri(ny, nx, domain, riv, hs, hr)
            integer, intent(in) :: ny, nx
            integer, intent(in) :: domain(ny,nx), riv(ny,nx)
            real(8), intent(inout) :: hs(ny,nx), hr(ny,nx)
        end subroutine unst2rri
    end interface
    integer, intent(in) :: ny, nx
    integer, intent(in) :: domain(ny, nx), riv(ny, nx)
    real(8), intent(in) :: time
    real(8), intent(inout) :: hs(ny,nx), hr(ny,nx)
    integer me, li, k, i

    ! write message - start UNST2D
    write(*,*) 'UNST2D timestep ===  ',int(time-timmax), '  >>>>>  ', int(time)

    ! UNST2D time set
    unsttime = time - timmax
    if(disk_flag) next_disk_t = next_disk_t + dkout
    if(disp_flag) next_disp_t = next_disp_t + dpout

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                    Loop Start
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    do while (unsttime + dt2 <= time)
        if(d1riv==1) then
            q_n_1d = q_1d
            a_n_1d = a_1d
            call calc_2d_to_1d_inflow
            call calc_1d_to_2d_outflow
            call d1riv_main(1)
        endif
        ! U^n -> U^*
        ! == Equation of motion (cal flux) ==
        um_n = um
        vn_n = vn
        unsth_n = unsth
        call flux   ! cal link flux
        call lkyokai ! apply boundary condition
        call limit_front
        ! == Equation of continuity (cal waterdepth) ==
        call suisin  ! cal mesh water depth
        ! cal velocity
        call velocity(1)
        ! update variable（new >>> old）
        call replace
        riv_eq = 0.0d0

        ! U^* -> U^n+1
        if(d1riv==1) then
            call calc_2d_to_1d_inflow
            call calc_1d_to_2d_outflow
            call d1riv_main(2)
        endif
        call flux   ! cal link flux
        um = 0.5d0 * (um_n + um)
        vn = 0.5d0 * (vn_n + vn)
        call lkyokai ! apply boundary condition
        call limit_front
        ! == Equation of continuity (cal waterdepth) ==
        call suisin  ! cal mesh water depth
        unsth = 0.5d0 * (unsth_n + unsth)
        ! cal velocity
        call velocity(2)

        ! update max record
        !$omp parallel do default(shared),private(me)
        do me = 1, mesh
            hmax(me) = max(hmax(me), unsth(me))
            uummax(me) = max(uummax(me), abs(uum(me)))
            vvmmax(me) = max(vvmmax(me), abs(vvm(me)))
        enddo
        !$omp end parallel do

        ! update variable（new >>> old）
        call replace

        ! reset subflow
        if(d1riv==1) then
            !$omp parallel do default(shared),private(i)
            do i = 1, ndan
                h_1dmax(i) = max(h_1dmax(i), h_1d(i))
                vv_1dmax(i) = max((vv_1dmax(i)), abs(vv_1d(i)))
            enddo
            !$omp end parallel do
        endif

        ! update time(2)
        unsttime = unsttime + dt2

        ! output result
        if(disk_flag) then
            call diskwrite   ! disk
            ! call dispwrite   ! display
            if(paddydam==1) call paddywrite  ! paddydam
        endif

        if(disp_flag) call dispwrite   ! display

        riv_eq = 0.0d0

        if(d1riv==1) call unstdt_adapt_1d
        call unstdt_adapt(time)

    enddo
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                      Loop End
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! UNST to RRI
    if(cnct_mode==1) call unst2rri(ny, nx, domain, riv, hs, hr)  ! select v.1.0.5

end subroutine UNST2D