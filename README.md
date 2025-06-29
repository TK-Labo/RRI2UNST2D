# RRI2UNST2D

## Citation
このコードを利用した計算結果の公表・頒布に際しては、以下の論文を引用してください。\
Please cite the following paper when publishing or distributing calculation results using this code.\
 \
川池 健司, 井上 和也, 戸田 圭一（2000）非構造格子の都市氾濫解析への適用，水工学論文集，44：461-466.\
Kenji KAWAIKE, Kazuya INOUE, Kei-ichi TODA (2000) Applications of unstructured meshes to inundation flow analysis in urban area, Annual journal of Hydraulic Engineering, JSCE, 44, 461-466.\
https://doi.org/10.2208/prohe.44.461 \
 \
山村 孝輝，西野 駿治，山田 真史，佐山 敬洋，川池健司，瀧健太郎（2025），流域治水計画検討のための降雨流出氾濫(RRI)モデルと非構造格子二次元不定流(UNST-2D)モデルの連成解析法の検討 河川技術論文集，31：469–474.\
Koki YAMAMURA, Shunji NISHINO, Masafumi YAMADA, Takahiro SAYAMA, Kenji KAWAIKE, Kentaro TAKI (2025) Proposal of coupled analysis of Rainfall-Runoff-Inundation (RRI) model and unstructured flood management planning, Advances in river engineering, JSCE, 31, 469-474.\

## Introduction
このプロジェクトは降雨流出氾濫（RRI: Rainfall-Runoff-Inundation）モデルと非構造格子二次元不定流モデル（UNST2D: Unstructured grid 2D unsteady flow model）を連携させることによって、降雨による流出・氾濫現象をより精度高くシミュレーションするためのモデルです。

## Coupled calculation with the RRI model
 \
RRIモデルとの連成計算を行う場合、国立研究開発法人土木研究所よりRRIモデル（ver.1.4.2.7）を入手し、RRI.f90に以下を追加します。\
https://www.pwri.go.jp/icharm/research/rri/index_j.html  
 \
### 変数追加（68行目に挿入）
  
```fortran
integer uflg
uflg = 0
```

### 主要なコード追加

#### 1. UNST-2D解析関連の追加（574行目に挿入）

```fortran
call unst_rdat(ny_rain, qp, nx_rain, tt_max_rain)
if(plantFN==1) call plantFNdat
if(plantDa==1) call plantDadat
if(paddydam==1) call paddydat
if(drainarea==1) call draindat
if(morid==1) call moriddat

call unst_initiald
if(paddydam==1) call paddyinitiald
if(drainarea==1) call draininitiald

if(str_type == 0) then
    phi(61) = 0.20d0
    phi(63) = 0.20d0
elseif(str_type == 1) then
    phi(61) = 0.05d0
    phi(63) = 0.20d0
endif

call diskwrite
call dispwrite
```

### 2. 降雨データ処理（525行目に挿入）

```fortran
dtrain = t_rain(1) - t_rain(0)
```

### 3. 出力タイムステップ設定（592~593行目の間に挿入）

```fortran
timmax = dble(dt) * dble(calldt)
```

### 4. UNST-2Dの処理ロジックの追加（1192行目に挿入）

```fortran
call unst_qin(qr_ave, qs_ave, uflg, hr, hs)
if(mod(time, timmax) == 0 .and. uflg==0) call predispwrite  
if(mod(time, timmax) == 0 .and. uflg==0) call prediskwrite
if(mod(time, timmax) == 0 .and. uflg==1) call UNST(hs, hr)
```

### 5. ファイルクローズ処理とメモリ解放処理（1203行目に挿入）

```fortran
call wrhmax

write(*, 1999) time
1999 format('      - normal end -  time=', f8.0)

if(str_type==1) then
    close(10077)
    close(10078)
endif
close(10091)
close(10093)
close(10094)
close(10095)
close(10096)
close(10097)
close(100101)
close(100102)
close(100103)

if(paddydam==1) then
    close(10098)
    close(10099)
    close(100100)
endif

deallocate(baseo, dnox, dnoy, smesh, scv, rthl, ux, uy, xmesh, ymesh, rtuv)
deallocate(limesh, linode, inf, ko, menode, melink, inl, qin, lkyokai_dx, lkyokai_dy)
deallocate(unsth, ho, hl, hmax, uummax, vvmmax)
deallocate(um, umo, umm, uu, vn, vno, vnm, vv)
deallocate(mn, rnof, lambda, rbeta, umbeta, vnbeta)
deallocate(uum, vvm, lhan, lhano, qr_sum, rnx, dl)
deallocate(unstc, co, str, stro, strmx)
if(plantFN==1) deallocate(plantF_array, plantN_array)
if(plantDa==1) deallocate(plant_D_array, plant_a_array, dk_val)
if(paddydam==1) deallocate(paddyid, pqout_idx, pdrain, min_pmeshid, device)
if(paddydam==1) deallocate(orifice_num, min_dist, psmesh, dr_dist, dhp, phid)
if(paddydam==1) deallocate(paddy_q, pqh, drain2phidx)
if(drainarea==1) deallocate(inf_dr, drp, drr, drr_dist, dhj)
if(dsmesh==1) deallocate(dsme)
if(ga==1) deallocate(genes)
```
