! Data read
! Coded by K.Kawaike and TK Labo
! Released on July 7th 2025

subroutine unst_rdat(ny_rain, temprain, nx_rain, tt_max_rain,&
     ny, nx, lasth, dt, &
     xllcorner_rain, yllcorner_rain, cellsize_rain_x, cellsize_rain_y,&
     xllcorner, yllcorner, cellsize)
    use unst_globals_mod
    implicit none

    integer, intent(in) :: tt_max_rain, nx_rain, ny_rain
    integer, intent(in) :: ny, nx, lasth, dt
    real(8), intent(in) :: temprain(0:tt_max_rain, ny_rain, nx_rain)
    real(8), intent(in) :: xllcorner_rain, yllcorner_rain, &
         cellsize_rain_x, cellsize_rain_y, &
         xllcorner, yllcorner, cellsize

    integer no, me, li, k, i, j, tt, dummy, read_count, ios
    character*50  fnode, flink, fmesh, finf, fbs, fqin, fmesh2ij
    character*50  fplantF, fplantN, fplantD, fplanta
    character*50  fpaddy, fpqout, fpaddy_param, fdhp, fpaddyh, fpaddyq
    character*50  finf_dr, fdrain, fmorid, fdsmesh, frn
    character*50  fh, fhmx, fuum, fvvm, fuumx, fvvmx, fstorage, fq
    character*50  fga, fpaddy_cluster
    character*50  fn(1)
    character*100 line
    integer, allocatable :: infmn(:)
    real(8), allocatable :: temp_mn(:)
    real(8), allocatable :: temp_soildepth(:), temp_gammaa(:), temp_ksv(:), temp_faif(:)
    integer n, numrn
    integer sep_rtuv
    integer, allocatable :: dsme(:), dstype(:), dsdepthnum(:)

    real(8) dkout, dpout

    data fn/ 'cntl.dat'/

    !parameter set
    open(10, file = fn(1), action = 'read')
    read(10, 1002) unstdt
    read(10, 1001) calldt
    read(10, 1002) dkout
    read(10, 1002) dpout
    read(10, 1002) ocpy
    read(10, 1002) unstbeta
    ! read(10, 1003) str_type  2507 廃止予定
    read(10, *)

    ! filename set
    read(10, 1009) fnode
    read(10, 1009) flink
    read(10, 1009) fmesh
    read(10, 1009) finf
    read(10, 1009) fbs
    read(10, 1009) fqin
    read(10, 1009) fmesh2ij
    read(10, 1009) frn
    read(10, *)
    read(10, *)

    read(10, 1003) dsmesh
    read(10, 1009) fdsmesh
    read(10, *)
    read(10, *)

    ! plant force filename set
    read(10, 1003) plantFN
    !filename
    read(10, 1009) fplantF
    read(10, 1009) fplantN
    read(10, *)
    read(10, *)

    read(10, 1003) plantDa
    !filename
    read(10, 1009) fplantD
    read(10, 1009) fplanta
    read(10, *)
    read(10, *)

    ! paddydam filename set
    read(10, 1003) paddydam
    ! filename
    read(10, 1009) fpaddy
    read(10, 1009) fpqout
    read(10, 1009) fpaddy_param
    read(10, *)
    read(10, 1009) fdhp
    read(10, 1009) fpaddyh
    read(10, 1009) fpaddyq
    read(10, *)
    read(10, *)

    ! drainarea filename set
    read(10, 1003) drainarea
    !filename
    read(10, 1009) finf_dr
    read(10, 1009) fdrain
    read(10, *)
    read(10, *)

    ! morid set
    read(10, 1003) morid
    ! filename
    read(10, 1009) fmorid
    read(10, *)
    read(10, *)

    ! 遺伝的アルゴリズム set
    ! read(10, 1003) ga
    ! filename
    ! read(10, 1009) fga
    ! read(10, 1009) fpaddy_cluster
    ! read(10,*) paddyclass
    ! read(10, *)
    ! read(10, *)

    read(10, 1009) fh
    read(10, 1009) fhmx
    read(10, 1009) fuum
    read(10, 1009) fvvm
    read(10, 1009) fuumx
    read(10, 1009) fvvmx
    read(10, 1009) fstorage
    read(10, 1009) fq

1001 format(7x, i11)
1002 format(7x, f11.2)
1003 format(17x, i1)
1009 format(a50)
    close(10)
!     ===============================================
    dtq = dble(dt)
    dt2 = 2.0d0*unstdt
    th = 1.0d-3
    gg = 9.8d0
    fita = 0.5d0
    pi = 3.14d0

    !write time interval
    lpout = int(dpout/unstdt)
    lkout = int(dkout/unstdt)
!     =====================================================
    !data files open
    ! 2D input
    open(21, file = fnode, action = 'read')
    open(22, file = flink, action = 'read')
    open(23, file = fmesh, action = 'read')
    open(24, file = finf, action = 'read')
    open(25, file = fbs, action = 'read')
    open(26, file = fqin, action = 'read')

    if(plantFN==1) then
        open(28, file = fplantF, action = 'read')
        open(29, file = fplantN, action = 'read')
    endif
    if(plantDa==1) then
        open(30, file = fplantD, action = 'read')
        open(31, file = fplanta, action = 'read')
    endif

    if(paddydam==1) then
        open(32, file = fpaddy, action = 'read')
        open(33, file = fpqout, action = 'read')
        open(36, file = fpaddy_param, action = 'read')

        ! output
        open(10098, file = fdhp, action = 'write')
        open(10099, file = fpaddyh, action = 'write')
        open(100100, file = fpaddyq, action = 'write')
    endif

    if(drainarea==1) then
        open(34, file = finf_dr, action = 'read')
        open(35, file = fdrain, action = 'read')
    endif

    open(37, file = fmesh2ij, action = 'read')
    if(dsmesh==1) open(38, file = fdsmesh, action = 'read')
    if(morid==1) open(39, file = fmorid, action = 'read')
    open(40, file = frn, action = 'read')

    ! 遺伝的アルゴリズム add
    ! if(ga==1) then
    !     open(41, file = fga, action = 'read')
    !     open(42, file = fpaddy_cluster, action = 'read')
    ! endif

    ! output file set
    if(str_type==1) then
    ! open(10071, file = 'out/unst/vl.dat', action = 'write')
    ! open(10072, file = 'out/unst/vl26.dat', action = 'write')
    ! open(10073, file = 'out/unst/vl3.dat', action = 'write')
    ! open(10074, file = 'out/unst/str45.dat', action = 'write')
    ! open(10076, file = 'out/unst/vl6.dat', action = 'write')
    open(10077, file = 'out/unst/tracer.dat', action = 'write')
    open(10078, file = 'out/unst/v_minus.dat', action = 'write')
    endif

    open(10091, file = fh, action = 'write')
    open(10093, file = fhmx, action = 'write')
    open(10094, file = fuum, action = 'write')
    open(10095, file = fvvm, action = 'write')
    open(10096, file = fuumx, action = 'write')
    open(10097, file = fvvmx, action = 'write')
    open(100101, file = fstorage, action = 'write')
    open(100102, file = fq, action = 'write')
    open(100103, file = 'out/unst/kyokaiq.dat', action = 'write')
    if(str_type==1) open(100110, file = 'out/unst/c.dat', action = 'write')

    unsttime = 0.0d0
    mstep = 0

    !node data
    read(21, *) node
    allocate(dnox(node), dnoy(node))
    do no = 1, node
        read(21, 1212) dnox(no), dnoy(no)
    enddo
1212 format(8x, 2f10.2)
    close(21)

    !link data
    read(22, *) link
    allocate(limesh(link, 2), linode(link, 2))
    allocate(scv(link), rthl(link, 2))
    allocate(ux(link), uy(link))
    do li = 1, link
        read(22, 1222) limesh(li, 1), limesh(li, 2), linode(li, 1), linode(li, 2)
        read(22, 1223) scv(li), rthl(li, 1), rthl(li, 2)
        read(22, 1224) ux(li), uy(li)
    enddo
1222 format(8x, 2i8, 2x, 2i8)
1223 format(8x, f10.2, 2f10.5)
1224 format(18x, 2f10.5)
    close(22)

    ! mesh data
    sep_rtuv = 0
    read(23, '(A)', iostat = ios) line

    ! count data
    line = adjustl(trim(line))
    read_count = 0
    do j = 1, len_trim(line)
        if (line(j:j) == ' ' .and. line(j+1:j+1) /= ' ') then
            read_count = read_count + 1
        endif
    enddo

    if (read_count == 0) then
        read(line, *, iostat = ios) mesh
        write(*,*) 'mesh - ', mesh
        write(*,*) '    interpolation type:', sep_rtuv
    elseif(read_count == 1) then
        read(line, *, iostat = ios) mesh, sep_rtuv
        write(*,*) 'mesh - ', mesh
        write(*,*) '    interpolation type:', sep_rtuv
    else
        write(*,*) 'ERROR : mesh.dat'
    endif

    allocate(ko(mesh), menode(mesh, 6), melink(mesh, 6), smesh(mesh), xmesh(mesh), ymesh(mesh))
    allocate(rtuv_x(mesh, 6), rtuv_y(mesh, 6))
    allocate(qr_sum(mesh))
    do me = 1, mesh
        read(23, 1232) ko(me), (menode(me, k), k = 1, ko(me))
        read(23, 1233) (melink(me, k), k = 1, ko(me))
        read(23, 1234) smesh(me), xmesh(me), ymesh(me)
        if (sep_rtuv==0) then
            ! weight distance only
            read(23, 1235) (rtuv_x(me, k), k = 1, ko(me))
            rtuv_y(me, 1:ko(me)) = rtuv_x(me, 1:ko(me))
        elseif (sep_rtuv==1) then
            ! weight distance and angle added by d.baba
            read(23, 1235) (rtuv_x(me, k), k = 1, ko(me))
            read(23, 1235) (rtuv_y(me, k), k = 1, ko(me))
        endif
    enddo
1232 format(8x, i5, 5x, 20i8)
1233 format(13x, 20i8)
1234 format(6x, 3f15.2)
1235 format(13x, 20f10.5)
    close(23)

    ! inf data read
    allocate(inf(mesh))
    do me = 1, mesh
        read(24, *) inf(me)
    enddo
    close(24)

    ! 標高データ
    allocate(baseo(mesh))
    do me = 1, mesh
        read(25, *) dummy, baseo(me)
    enddo
    close(25)

    ! 占有率
    allocate(lambda(mesh))
    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        lambda(me) = 0.0d0
    enddo
    !$omp end parallel do

    ! 通過率
    allocate(rbeta(link))
    !$omp parallel do default(shared),private(li)
    do li = 1, link
    rbeta(li) = unstbeta
    if(limesh(li,2)==0) then
        if(inf(limesh(li,1))==63) then
            rbeta(li) = 0.0d0
            goto 3001
        endif
    endif
3001 enddo
    !$omp end parallel do

    ! 流入データ
    iqin = lasth * 3600 / dt
    print *, iqin
    read(26, *) iqnum
    allocate(inl(iqnum), qin(iqnum, iqin+1), qinu(iqnum, iqin+1), qinv(iqnum, iqin+1) )
    allocate(lkyokai_dx(iqnum), lkyokai_dy(iqnum))

    qin = 0.0d0
    qinu = 0.0d0
    qinv = 0.0d0

    do i = 1, iqnum
        read(26, *) inl(i)
    enddo
    close(26)

    ! rain data read
    allocate(rain(0:tt_max_rain, ny_rain, nx_rain))
    do tt = 0, tt_max_rain
        do i = 1, ny_rain
            do j = 1, nx_rain
                rain(tt, i, j) = temprain(tt, i, j)
            enddo
        enddo
    enddo

    ! メッシュの緯度経度を読み込む
    allocate(urain_j(mesh), urain_i(mesh))
    allocate(rri_x(mesh), rri_y(mesh))
    do me = 1, mesh
     read(37, *) rri_x(me), rri_y(me)
    end do

    do me = 1, mesh
        urain_j(me) = int( (rri_x(me) - xllcorner_rain) / cellsize_rain_x) + 1
        urain_i(me) = ny_rain - int( (rri_y(me) - yllcorner_rain) / cellsize_rain_y)
        if(urain_i(me) <=0 .or. urain_j(me) <= 0) write(*,*) 'error rain file or mesh4326.dat'
    enddo
    close(37)

    !rsetsuzoku set
    allocate( rsetsu_i(mesh), rsetsu_j(mesh) )
    do me = 1, mesh
        rsetsu_j(me) = int( (rri_x(me) - xllcorner) / cellsize) + 1
        rsetsu_i(me) = ny - int( (rri_y(me) - yllcorner) / cellsize)
        if(rsetsu_i(me) <=0 .or. rsetsu_j(me) <= 0) write(*,*) 'error qinfile'
    enddo

    ! downstream data read
    if(dsmesh==1) then
        read(38, *) dsmenum
        allocate( dsme(dsmenum), dstype(dsmenum), dsdepthnum(dsmenum))
        allocate( dsdt(mesh), dsinf(mesh), dsupper(mesh), dsfilter2(dsmenum))
        dstype = 1
        dsinf = 0
        dsfilter2 = 0
        idsfilter2 = 0
        dsdepthnum = 0
        dsdt = lasth * 3600
        read(38, *)
        do i = 1, dsmenum
            read(38, '(A)', iostat = ios) line
            if (ios /= 0) exit

            ! count data
            line = adjustl(trim(line))
            read_count = 0
            do j = 1, len_trim(line)
                if (line(j:j) == ' ' .and. line(j+1:j+1) /= ' ') then
                    read_count = read_count + 1
                endif
            enddo

            ! 読み込んだデータ数に応じて処理
            if (read_count == 0) then
                read(line, *) dsme(i)
            elseif (read_count == 1) then
                read(line, *) dsme(i), dstype(i)
            elseif (read_count == 2) then
                read(line, *) dsme(i), dstype(i), dsdepthnum(i)
            endif

            dsinf(dsme(i)) = dstype(i)
            idsdepth = max(dsdepthnum(i), idsdepth, 0)
            ! cal time step
            if(dstype(i)==2) then
                idsfilter2 = idsfilter2 + 1
                dsfilter2(idsfilter2) = dsme(i)
                if (dsdepthnum(i)<=0) then
                    write(*,*) 'ERROR : meshid' , dsme(i) , 'is illegal value, No discharge.'
                    dstype(i) = 0
                    dsinf(dsme(i)) = 0
                    idsfilter2 = idsfilter2 - 1
                else
                    dsdt(dsme(i)) = lasth * 3600/dsdepthnum(i)
                endif
            endif
        enddo

        ! dstype2: water level(m)
        if (any(2 == dstype)) then
            read(38, *)
            allocate(dsdepth(mesh,idsdepth+1))

            dsdepth = 0.0d0

            ! downstream mesh rtuv set
            allocate(sum_rtuv_x(mesh), sum_rtuv_y(mesh))
        endif

        if (any(3 == dstype)) write(*,*) 'UNIFORM FLOW: target meshid - upper meshid'

        do i = 1, dsmenum
            ! dstype1: 自由流出する下流端メッシュのセット
            if (dstype(i)==1) then
                baseo(dsme(i)) = -9999.0d0

            ! dstype2: water level(m)
            elseif (dstype(i)==2) then
                ! read water level(m)
                read(38, *) (dsdepth(dsme(i),j),j=1,dsdepthnum(i))
                dsdepth(dsme(i),dsdepthnum(i)+1) = dsdepth(dsme(i),dsdepthnum(i))
                ! cal downstream rtuv set
                sum_rtuv_x(dsme(i)) = 0.0d0
                sum_rtuv_y(dsme(i)) = 0.0d0
                do k = 1, ko(dsme(i))
                    if (limesh(melink(dsme(i),k),2)/=0) then
                        sum_rtuv_x(dsme(i)) = sum_rtuv_x(dsme(i)) + rtuv_x(dsme(i),k)
                        sum_rtuv_y(dsme(i)) = sum_rtuv_y(dsme(i)) + rtuv_y(dsme(i),k)
                    endif
                enddo
                do k = 1, ko(dsme(i))
                    if (limesh(melink(dsme(i),k),2)/=0) then
                        rtuv_x(dsme(i),k) = rtuv_x(dsme(i), k)/sum_rtuv_x(dsme(i))
                        rtuv_y(dsme(i),k) = rtuv_y(dsme(i), k)/sum_rtuv_y(dsme(i))
                    else
                        rtuv_x(dsme(i),k) = 0.0d0
                        rtuv_y(dsme(i),k) = 0.0d0
                    endif
                enddo

            ! dstype3: uniform flow(beta)
            ! only square mesh -- copy upper mesh env.
            elseif (dstype(i)==3) then
                if (ko(dsme(i))==4) then
                    do k = 1, ko(dsme(i))
                        ! set upper meshid
                        if (limesh(melink(dsme(i), k),2) == 0) then
                            if (limesh(melink(dsme(i), mod(k+1,4)+1),1) == dsme(i)) then
                                dsupper(dsme(i)) = limesh(melink(dsme(i), mod(k+1,4)+1),2)
                            else
                                dsupper(dsme(i)) = limesh(melink(dsme(i), mod(k+1,4)+1),1)
                            endif
                            write(*,*) dsme(i) , '-' , dsupper(dsme(i))
                        endif
                    enddo
                else
                    write(*,*) ' ERROR : meshid', dsme(i), 'is not square, No discharge.'
                    dsinf(dsme(i)) = 0
                endif
            endif
        enddo
        close(38)
        if (any(2 == dstype)) deallocate(sum_rtuv_x, sum_rtuv_y)
        deallocate(dsme, dstype, dsdepthnum)
    endif


    ! 雨水トレーサー濃度
    allocate( unstc(mesh), co(mesh), cr(mesh), str(mesh), stro(mesh), strmx(mesh) )
    !$omp parallel do default(shared),private(me)
    do me = 1, mesh
        cr(me) = 0.0d0
        if(inf(me)==61) cr(me) = 1.0d0
    enddo
    !$omp end parallel do

    ! 粗度係数データ read
    allocate( mn(mesh) )
    allocate( unst_soildepth(mesh), unst_gammaa(mesh), unst_ksv(mesh), unst_faif(mesh), unst_infilt_limit(mesh) )
    allocate( unst_gampt_ff(mesh), unst_gampt_f(mesh) )
    read(40, *) numrn
    read(40, *)

    ! 一旦必要な配列を確保
    allocate( infmn(numrn), temp_mn(numrn) )
    allocate( temp_soildepth(numrn), temp_gammaa(numrn), temp_ksv(numrn), temp_faif(numrn) )

    ! 初期化
    temp_soildepth = 0.0d0
    temp_gammaa = 0.0d0
    temp_ksv = 0.0d0
    temp_faif = 0.0d0

    ! ファイルから読み込み
    do n = 1, numrn
        ! 1行読み込み
        read(40, '(A)', iostat = ios) line
        if (ios /= 0) exit

        ! 行内のデータ数をカウント
        read_count = 0
        do i = 1, len_trim(line)
            if (line(i:i) == ' ' .and. line(i+1:i+1) /= ' ') then
                read_count = read_count + 1
            endif
        enddo
        read_count = read_count + 1

        ! 読み込んだデータ数に応じて処理
        if (read_count >= 2) then
            read(line, *) infmn(n), temp_mn(n)
        endif

        if (read_count >= 7) then
            read(line, *) infmn(n), temp_mn(n), temp_soildepth(n), temp_gammaa(n), temp_ksv(n), temp_faif(n)
        endif
    enddo
    close(40)

    !$omp parallel do default(shared),private(me, n)
    do me = 1, mesh
        ! 初期値を0に設定
        unst_soildepth(me) = 0.0d0
        unst_gammaa(me) = 0.0d0
        unst_ksv(me) = 0.0d0
        unst_faif(me) = 0.0d0
        unst_infilt_limit(me) = 0.0d0
        unst_gampt_ff(me) = 0.0d0
        unst_gampt_f(me) = 0.0d0

        do n = 1, numrn
            if(inf(me) == infmn(n)) then
                mn(me) = temp_mn(n)
                unst_soildepth(me) = temp_soildepth(n)
                unst_gammaa(me) = temp_gammaa(n)
                unst_ksv(me) = temp_ksv(n)
                unst_faif(me) = temp_faif(n)
                if( unst_soildepth(me) .gt. 0.d0 .and. unst_ksv(me) .gt. 0.d0 ) unst_infilt_limit(me) = unst_soildepth(me) * unst_gammaa(me)
            endif
        enddo
    enddo
    !$omp end parallel do

    ! 一時配列を解放
    deallocate( infmn, temp_mn, temp_soildepth, temp_gammaa, temp_ksv, temp_faif )

    !allocate variables
    allocate(unsth(mesh), ho(mesh), hl(link), hmax(mesh), uummax(mesh), vvmmax(mesh))
    allocate(um(link), umo(link), umm(0:mesh), uu(link), vn(link), vno(link), vnm(0:mesh), vv(link))
    allocate(rnof(mesh), umbeta(link), vnbeta(link))
    allocate(uum(mesh), vvm(mesh))
    allocate(lhan(link), lhano(link))
    allocate(rnx(link), dl(link))
    allocate(node_dx(mesh,6), node_dy(mesh,6))

end subroutine unst_rdat


! ヨシ帯抵抗（倒伏過程考慮）data read
! 開発中 July 7th 2025
subroutine plantFNdat
    use unst_globals_mod
    implicit none
    integer me

    ! plant resistance / 1 (syokusei teiko / hon)
    allocate(plantF_array(mesh), plantN_array(mesh))
    do me = 1, mesh
          read(28, *) plantF_array(me)
    enddo

    ! number of plants(honsu)
    do me = 1, mesh
        read(29, *) plantN_array(me)
    enddo
    close(29)
!
end subroutine plantFNdat


! 樹林帯抵抗 data read
subroutine plantDadat
    use unst_globals_mod
    implicit none
    integer me

    ! 幹の直径　(m)
    allocate(plant_D_array(mesh), plant_a_array(mesh), plant_lambda(mesh))
    allocate(dk_val(link))
    do me = 1, mesh
        read(30, *) plant_D_array(me)
    enddo
    close(30)

    ! 樹木密度(本 / m^2)
    do me = 1, mesh
        read(31, *) plant_a_array(me)
    enddo
    close(31)
!
end subroutine plantDadat


! 田んぼダム data read
subroutine paddydat
    use unst_globals_mod
    implicit none
    integer me, pa, ii

    ! 水田id read
    allocate(paddyid(mesh))
    do me = 1, mesh
        read(32,*) paddyid(me)
    enddo
    close(32)

    ! 水路までの距離
    read(33, *) paddy
    read(33, *)
    allocate(pqout_idx(paddy), pdrain(paddy), min_pmeshid(paddy), device(paddy))
    allocate(orifice_num(paddy), min_dist(paddy), etp(paddy))
    allocate(dr_dist(paddy), pqh(paddy), paddy_q(paddy), totalqp(paddy))
    if(paddy>0) then
    do pa = 1, paddy
          read(33, *) pqout_idx(pa), min_pmeshid(pa), min_dist(pa), orifice_num(pa), &
                         pdrain(pa), dr_dist(pa), device(pa), etp(pa)
    enddo
    read(33, *)
    read(33, *)
    read(33, *) nhp
    allocate(dhp(72000, nhp), phid(nhp))
    do ii = 1, nhp
          read(33, *) phid(ii)
    enddo
    endif
    close(33)

    ! 遺伝的アルゴリズム　add  2507 廃止予定
    ! if(ga==1) then
    !     allocate( genes(paddyclass) )
    !     read(41, *) (genes(pa), pa = 1, paddyclass)
    !     close(41)

    !     allocate( paddycluster(mesh) )
    !     do me = 1, mesh
    !         read(42,*) paddycluster(me)
    !         if(paddycluster(me) > 0) then
    !             device(paddyid(me)) = genes(paddycluster(me))
    !         endif
    !     enddo
    !     close(42)
    !     deallocate(paddycluster)
    ! endif

    ! paddy field parameters
    read(36, 1361) lh           ! 畦畔高
    read(36, 1363) wh1          ! 分離型のセキ板の高さ
    read(36, 1361) wh2          ! 一体型のセキ板の高さ
    read(36, 1361) ww1          ! 落水口幅
    read(36, 1361) ww2          ! 器具の切欠幅
    read(36, 1362) ca           ! 器具の中心角
    read(36, 1362) wtyp         ! セキのタイプ（1:四角セキ, 2:三角セキ）
    read(36, 1361) dld          ! 畦畔天端と器具上端の高さの差
    read(36, 1361) dd           ! 器具の穴の直径
    read(36, 1361) unstdh       ! 田面から器具の穴中心までの高さ
    read(36, 1361) p_data       ! 排水溝パイプ直径
    read(36, 1361) ph           ! 田面からパイプ中心までの高さ
1361 format(7x, f11.2)
1362 format(7x, i11)
1363 format(7x, f11.3)
!
    !allocate variables
    allocate(ttp(mesh))

end subroutine paddydat


! sewerage and fields are data read
subroutine draindat
    use unst_globals_mod
    implicit none
    integer me, nn, drn

    ! inf_dr(sewerage and fields id) data read
    ! (0:not 1~:sewerage and fields are)
    allocate(inf_dr(mesh))
    do me = 1, mesh
    read(34, *) inf_dr(me)
    enddo
    close(34)

    ! drainage system planning　data read
    ! (keikaku haisui noryoku)
    read(35, *) dr_no
    if(dr_no>0) then
    read(35, *)
    allocate(drp(dr_no), drc(dr_no))
    allocate(drr(dr_no), drr_dist(dr_no))
    allocate(vol_dr(mesh), vol(mesh))
    allocate(dhj(72000, dr_no))
    do nn = 1, dr_no
          read(35, *)  drn, drr(nn), drp(nn), drc(nn), drr_dist(nn)
          if (nn == dr_no) exit
          if (drr(nn) < 0.0) then
                print *, 'Enter the planned treatment rainfall.'
          endif
    enddo
    endif

    !allocate variables
    allocate(tripTime(mesh))

end subroutine draindat


! morido are data read
subroutine moriddat
    use unst_globals_mod
    implicit none
    integer mmo, li

    read(39, *) mmorid
    read(39, *)

    allocate( nmorili(mmorid), infl(link) )
    allocate( zbbk(link) )

    do mmo = 1, mmorid
        read(39, *) nmorili(mmo), zbbk(nmorili(mmo))
        do li = 1, link
            if(li == nmorili(mmo)) infl(li)=1
        enddo
    enddo

end subroutine moriddat
