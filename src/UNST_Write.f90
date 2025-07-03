! Write data
! Coded by K.Kawaike and TK Labo
! Release on July 7th 2025

module write_procedures
    use unst_globals_mod
    contains
    ! 配列全体書き込み用サブルーチン
    subroutine write_array_data(unit_num, data, fmt_data)
        implicit none
        integer, intent(in) :: unit_num
        real(8), intent(in) :: data(:)
        character(*), intent(in) :: fmt_data
        integer :: me

        write(unit_num, '(a,f8.0,a)') ' time=', unsttime, '(s)'
        write(unit_num, fmt_data) (data(me), me = 1, mesh)  ! meshもグローバル変数
    end subroutine

    subroutine write_paddyarray_data(unit_num, data, fmt_data)
        implicit none
        integer, intent(in) :: unit_num
        real(8), intent(in) :: data(:)
        character(*), intent(in) :: fmt_data
        integer :: i

        write(unit_num, '(a,f8.0,a)') ' time=', unsttime, '(s)'
        write(unit_num, fmt_data) (data(i), i = 1, paddy)  ! meshもグローバル変数
    end subroutine

    ! paddy書き込み用サブルーチン
    subroutine write_paddy_data(unit_num)
        implicit none
        integer, intent(in) :: unit_num
        integer :: me, i

        write(unit_num, 1040) unsttime
        do me = 1, nhp
            do i = 1, 72000
                if(dhp(i,me) > 0.0d0) write(unit_num, 1041) i, dhp(i,me), me
            enddo
        enddo
    1040 format(' time=', f8.0, '(s)')
    1041 format(i8, f10.5, i8)
    end subroutine

    !RRI_UNST専用
    subroutine prewrite_array_data(unit_num, data, time, fmt_data)
        implicit none
        integer, intent(in) :: unit_num
        real(8), intent(in) :: data(:)
        real(8), intent(in) :: time
        character(*), intent(in) :: fmt_data
        integer :: me

        write(unit_num, '(a,f8.0,a)') ' time=', time - (timmax / 2.0d0), '(s)'
        write(unit_num, fmt_data) (data(me), me = 1, mesh)
        write(unit_num, '(a,f8.0,a)') ' time=', time, '(s)'
        write(unit_num, fmt_data) (data(me), me = 1, mesh)
    end subroutine

    subroutine prewrite_paddyarray_data(unit_num, data, time, fmt_data)
        implicit none
        integer, intent(in) :: unit_num
        real(8), intent(in) :: data(:)
        real(8), intent(in) :: time
        character(*), intent(in) :: fmt_data
        integer :: i

        write(unit_num, '(a,f8.0,a)') ' time=', time - (timmax / 2.0d0), '(s)'
        write(unit_num, fmt_data) (data(i), i = 1, paddy)  ! meshもグローバル変数
        write(unit_num, '(a,f8.0,a)') ' time=', time, '(s)'
        write(unit_num, fmt_data) (data(i), i = 1, paddy)  ! meshもグローバル変数
    end subroutine
end module write_procedures

subroutine wrfile
    use unst_globals_mod
    use write_procedures
    implicit none
    real(8) sv, sa, saj

    ! display indication
    entry dispwrite
    call sumqa(sv, sa, saj)
    write(*, 1000) unsttime, sv
    if(str_type == 1) then
    !write(10071,1000) unsttime, sv, sa, saj
    ! write(10071,1000) unsttime, sv, sv3, sv6, sv26, sumstr45
    ! write(10073,1000) unsttime, sv3, v_minus(31)+v_minus(32)
    ! write(10076,1005) unsttime, sv6, v_minus(61), count1, count2, count3, count4
    ! write(10072,1000) unsttime, sv26, v_minus(26)
    ! write(10074,1000) unsttime, sumstr45, sumstr4, sumstr5
    write(10078,1000) unsttime, v_minus_all, v_minus(21), v_minus(22), v_minus(23), v_minus(24), v_minus(25)
    write(10077,1099) unsttime, svc, v_cminus, v_cplus, v_cextre
    endif
1000 format('UNST------   unsttime=', f8.0, '(s)', 10f18.4)
! 1005 format('   unsttime=', f8.0, '(s)', 2f18.4, 4i10)
1099 format('   unsttime=', f8.0, '(s)', 10f21.7)
    return

    ! writing data to a file
    entry diskwrite

    !$omp parallel
    !$omp sections
    !$omp section
    call write_array_data(10091, unsth, '(10f15.3)')

    !$omp section
    call write_array_data(10094, uum, '(10f15.3)')

    !$omp section
    call write_array_data(10095, vvm, '(10f15.3)')

    !$omp section
    call write_array_data(100101, ((unsth*smesh) - qr_sum), '(10e15.7)')

    !$omp section
    call write_array_data(100102, qr_sum, '(10f15.3)')

    !$omp section
    if(str_type == 1) call write_array_data(100110, unstc, '(10f7.4)')

    !$omp end sections
    !$omp end parallel
    return

    ! Maximum depth of inundation written
    ! Maximum water level written
    entry wrhmax

    !$omp parallel
    !$omp sections
    !$omp section
    call write_array_data(10093, hmax,  '(10f15.3)')

    !$omp section
    call write_array_data(10096, uummax,  '(10f15.3)')

    !$omp section
    call write_array_data(10097, vvmmax,  '(10f15.3)')
    !$omp end sections
    !$omp end parallel
    return

    entry paddywrite

    !$omp parallel
    !$omp sections
    !$omp section
    call write_paddy_data(10098)

    !$omp section
    call write_paddyarray_data(10099, pqh,  '(10f15.9)')

    !$omp section
    call write_paddyarray_data(100100, totalqp,  '(10e15.7)')
    !$omp end sections
    !$omp end parallel
    return

end subroutine wrfile

subroutine prewrfile(time)
    use unst_globals_mod
    use write_procedures
    implicit none

    real(8), intent(in) :: time

    ! display indication
    entry predispwrite

    write(*, 1000) time - (timmax / 2.0d0), 0.0d0
    if(str_type == 1) then
    !write(10071,1000) time - (timmax / 2.0d0), sv, sa, saj
    ! write(10071,1000) time - (timmax / 2.0d0), sv, sv3, sv6, sv26, sumstr45
    ! write(10073,1000) time - (timmax / 2.0d0), sv3, v_minus(31)+v_minus(32)
    ! write(10076,1005) time - (timmax / 2.0d0), sv6, v_minus(61), count1, count2, count3, count4
    ! write(10072,1000) time - (timmax / 2.0d0), sv26, v_minus(26)
    ! write(10074,1000) time - (timmax / 2.0d0), sumstr45, sumstr4, sumstr5
    write(10078,1000) time - (timmax / 2.0d0), v_minus_all, v_minus(21), v_minus(22), v_minus(23), v_minus(24), v_minus(25)
    write(10077,1099) time - (timmax / 2.0d0), svc, v_cminus, v_cplus, v_cextre
    endif

    write(*, 1000) time, 0.0d0
    if(str_type == 1) then
    !write(10071,1000) time, sv, sa, saj
    ! write(10071,1000) time, sv, sv3, sv6, sv26, sumstr45
    ! write(10073,1000) time, sv3, v_minus(31)+v_minus(32)
    ! write(10076,1005) time, sv6, v_minus(61), count1, count2, count3, count4
    ! write(10072,1000) time, sv26, v_minus(26)
    ! write(10074,1000) time, sumstr45, sumstr4, sumstr5
    write(10078,1000) time, v_minus_all, v_minus(21), v_minus(22), v_minus(23), v_minus(24), v_minus(25)
    write(10077,1099) time, svc, v_cminus, v_cplus, v_cextre
    endif
1000 format('UNST------   unsttime=', f8.0, '(s)', 10f18.4)
! 1005 format('   unsttime=', f8.0, '(s)', 2f18.4, 4i10)
1099 format('   unsttime=', f8.0, '(s)', 10f21.7)
    return

    ! writing data to a file
    entry prediskwrite

    !$omp parallel
    !$omp sections
    !$omp section
    call prewrite_array_data(10091, unsth, time, '(10f15.3)')

    !$omp section
    call prewrite_array_data(10094, uum, time, '(10f15.3)')

    !$omp section
    call prewrite_array_data(10095, vvm, time, '(10f15.3)')

    !$omp section
    call prewrite_array_data(100101, ((unsth*smesh) - qr_sum), time, '(10e15.7)')

    !$omp section
    call prewrite_array_data(100102, qr_sum, time, '(10f15.3)')

    !$omp section
    if(str_type == 1) call prewrite_array_data(100110, unstc, time, '(10f7.4)')

    !$omp section
    if(paddydam == 1) then
        write(10098, 7010) time - (timmax / 2.0d0)
        write(10098, 7010) time
    7010 format(' time=', f8.0, '(s)')

        call prewrite_paddyarray_data(10099, pqh, time, '(10f15.9)')
        call prewrite_paddyarray_data(100100, totalqp, time, '(10e15.7)')
    endif
    !$omp end sections
    !$omp end parallel

    return

end subroutine prewrfile
