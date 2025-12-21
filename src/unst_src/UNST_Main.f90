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
    use unst_1d_main
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
    mstep = int(unsttime / unstdt)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                    Loop Start
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    do while (unsttime + unstdt <= time)
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !  Equation of motion (cal flux)
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    call flux   ! cal link flux
    call lkyokai ! apply boundary condition

    ! flow flux(velocity) ahead of the flood inundation front set
    !$omp parallel do default(shared),private(me,li,k)
    do me = 1, mesh
        if(unsth(me) >= th) cycle
        do k = 1, ko(me)
            li = melink(me, k)
            if((um(li)*node_dy(me,k) - vn(li)*node_dx(me,k)) > 0.0d0) then
                um(li) = 0.0d0
                vn(li) = 0.0d0
            endif
        enddo
    enddo
    !$omp end parallel do

    if(d1riv==1) then
        call calc_2d_to_1d_inflow
        call calc_1d_to_2d_outflow
        call d1riv_main(unsttime)
    endif

    ! update time(1)
    unsttime = unsttime + unstdt
    mstep = mstep + 1

    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !  Equation of continuity (cal waterdepth)
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    call suisin  ! cal mesh water depth

    ! cal velocity
    call velocity

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
        riv_eq = 0.0d0
        !$omp parallel do default(shared),private(i)
        do i = 1, ndan
            h_1dmax(i) = max(h_1dmax(i), h_1d(i))
        enddo
        !$omp end parallel do
    endif

    ! update time(2)
    unsttime = unsttime + unstdt
    mstep = mstep + 1

    ! output result
    if(mod(mstep, lkout) == 0) call diskwrite   ! disk
    if(mod(mstep, lpout) == 0) call dispwrite   ! display
    if(paddydam==1 .and. mod(mstep, lkout) == 0) call paddywrite  ! paddydam

    enddo
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                      Loop End
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! UNST to RRI
    if(cnct_mode==1) call unst2rri(ny, nx, domain, riv, hs, hr)

end subroutine UNST2D
