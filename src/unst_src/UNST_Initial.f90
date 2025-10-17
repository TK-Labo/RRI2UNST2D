! Set initiald conditions
! Coded by K.Kawaike and TK Labo
! Released on July 7th 2025

subroutine unst_initiald
    use unst_globals_mod
    implicit none
    integer me, li, k, k2, j
    real(8) plant_lambda_link, mesh_dx, mesh_dy

    !-----------------------
    ! Initialize time param
    !-----------------------
    unsttime = 0.0d0
    mstep = 0

    !---------------------
    ! Initialize variable
    !---------------------
    unst_error_v = 0.0d0  ! v1.0.5
    if(dsmesh==1) unst_dis_v = 0.0d0  ! v1.0.5

    ! coefficient of roughness by inf
    !$omp parallel do default(shared),private(me, k, k2)
    do me = 1, mesh
        unsth(me) = 0.0d0

        if (dsmesh==1) then
            if (ds_inf(me)==2) unsth(me) = ds_wl(me,1)
        endif
        ho(me) = unsth(me)
        hmax(me) = unsth(me)
        qr_sum(me) = 0.0d0
        rnof(me) = 1.00d0

        if(plantDa==1) then
            if (plant_a_array(me) /= 0.0d0 .and. plant_D_array(me) /= 0.0d0) then
                plant_lambda(me) = plant_a_array(me)*plant_D_array(me)
            else
                plant_lambda(me) = 0.0d0
            endif
        endif

        do k = 1, ko(me)
            k2 = mod(k, ko(me)) + 1
            node_dx(me, k) = dnox(menode(me, k2)) - dnox(menode(me, k))
            node_dy(me, k) = dnoy(menode(me, k2)) - dnoy(menode(me, k))
        enddo

    enddo
    !$omp end parallel do

    !$omp parallel do default(shared),private(li, plant_lambda_link)
    do li = 1, link
        um(li)  = 0.0d0
        umo(li) = 0.0d0
        vn(li)  = 0.0d0
        vno(li) = 0.0d0
!     ---------------------------------- 補間水深
        if(limesh(li, 2) /= 0) then
        hl(li) = unsth(limesh(li, 1))*rthl(li, 1) + unsth(limesh(li, 2))*rthl(li, 2)
        else
        hl(li) = unsth(limesh(li, 1))
        endif
        lhan(li) = 0
        lhano(li) = lhan(li)

        dl(li) = 0.0d0
        rnx(li) = 0.0d0
        if(limesh(li,2) /= 0) then
            mesh_dx = xmesh(limesh(li, 2)) - xmesh(limesh(li, 1))
            mesh_dy = ymesh(limesh(li, 2)) - ymesh(limesh(li, 1))
            dl(li) = sqrt(mesh_dx**2 + mesh_dy**2)
            rnx(li) = 0.5d0*(mn(limesh(li, 1)) + mn(limesh(li, 2)))
        endif

        if(plantDa==1) then
            if(limesh(li,2)/=0) then
                plant_lambda_link = 0.5d0*(plant_lambda(limesh(li, 1))+plant_lambda(limesh(li, 2)))
            else
                plant_lambda_link = plant_lambda(limesh(li, 1))
            endif
            dk_val(li) = 0.0d0
            ! 抗力係数CD = 1.2d0
            if(plant_lambda_link>0.0d0) dk_val(li) = 1.0d0/sqrt(plant_lambda_link*1.2d0/(2.0d0*gg))
        endif
    enddo
    !$omp end parallel do

    ! 上流端のdx,dy計算
    lkyokai_dx = 0.0d0
    lkyokai_dy = 0.0d0
    !$omp parallel do default(shared),private(j)
    do j = 1, iqnum
        lkyokai_dx(j) = dnox(linode(inl(j), 1)) - dnox(linode(inl(j), 2))
        lkyokai_dy(j) = dnoy(linode(inl(j), 1)) - dnoy(linode(inl(j), 2))
    enddo
    !$omp end parallel do

    ! 境界参照用
    umm(0) = 0.0d0
    vnm(0) = 0.0d0

end subroutine unst_initiald
