# RRI2UNST2D

https://github.com/user-attachments/assets/d922fac3-b09c-4861-9c57-251e41bb8a58

## Introduction
このプロジェクトでは、降雨流出氾濫 (RRI: Rainfall-Runoff-Inundation) モデルと非構造格子二次元不定流モデル (UNST2D: Unstructured grid 2D unsteady flow model) を連携させ、降雨による流出・氾濫現象をより精度高くシミュレーションするためのモデルを開発しています。  
In this project, we are developing a model to more accurately simulate runoff and flooding phenomena caused by rainfall by linking the Rainfall-Runoff-Inundation (RRI) model with the Unstructured grid 2D unsteady flow model (UNST2D).  

降雨流出氾濫（RRI）モデルと非構造格子二次元不定流（UNST2D）モデルを結合し，連成計算が可能な解析法を提案しています．RRIモデルをベースに任意領域でUNST-2Dの分担範囲を選択し，RRIモデルで計算された流量フラックスをUNST-2Dモデルの外縁に境界条件として与え，UNST-2Dモデルで計算された水深を重複領域上にあるRRIモデルの各メッシュに水位として逐次返すことで，両モデルの計算結果を整合させながら連成計算を進めます．

## Citation
このコードを利用した計算結果の公表・頒布に際しては、以下の論文を引用してください。  
Please cite the following paper when publishing or distributing calculation results using this code.

川池健司, 井上和也, 戸田圭一 (2000) 非構造格子の都市氾濫解析への適用, 水工学論文集, 44: 461-466.  
Kenji KAWAIKE, Kazuya INOUE, Kei-ichi TODA (2000) Applications of unstructured meshes to inundation flow analysis in urban area, Annual journal of Hydraulic Engineering, JSCE, 44, 461-466.
https://doi.org/10.2208/prohe.44.461

山村孝輝, 西野駿治, 山田真史, 佐山敬洋, 川池健司, 瀧健太郎 (2025) 流域治水計画検討のための降雨流出氾濫 (RRI) モデルと非構造格子二次元不定流 (UNST-2D) モデルの連成解析法の検討，河川技術論文集, 31: 469–474.  
Koki YAMAMURA, Shunji NISHINO, Masafumi YAMADA, Takahiro SAYAMA, Kenji KAWAIKE, Kentaro TAKI (2025) Proposal of coupled analysis of Rainfall-Runoff-Inundation (RRI) model and unstructured flood management planning, Advances in river engineering, JSCE, 31, 469-474.

RRIモデルとの連成計算を行う場合は、あわせてRRIモデルのプログラム利用許諾規約にしたがってください。  
If you are performing coupled calculations with the RRI model, please also quote the information specified in the license for the RRI model.
https://www.pwri.go.jp/icharm/research/rri/index_j.html

## Contents

プロジェクトは次の主要なディレクトリとファイルで構成されています:

- `src/`: Fortranソースコードが格納されたディレクトリ
  - UNST2Dモデル関連のファイル（`UNST*.f90`）
  - Makefile

The project consists of the following main directories and files:

- `src/`: Directory where Fortran source code is stored
  - UNST2D model related files（`UNST*.f90`）
  - Makefile


## Compile
RRIとUNST2Dによる連成計算を行うプログラム`RRI_UNST.exe`を生成する手順は以下の通りです。

1. RRI関連コードの内、`RRI-CUI/source/1.4.2.7/*.f90`を本プロジェクトの`src`ディレクトリにコピーしてください。
2. RRIのソースコードに連成計算に必要なコードを追加してください。コードの追加方法については [Coupled calculation with the RRI model](#Coupled calculation with the RRI model) を参照してください。
3. `src`ディレクトリで以下のコマンドを実行すると、`RRI_UNST.exe`が生成されます。

The procedure for generating the program RRI_UNST.exe that performs coupled calculations using RRI and UNST2D is as follows.  
  
1. Among the RRI related codes, copy `RRI-CUI/source/1.4.2.7/*.f90` to the `src` directory of this project.
2. Please add the code required for coupled calculation to the RRI source code. For how to add the code, see [Coupled calculation with the RRI model](#Coupled calculation with the RRI model).
3. Execute the following command in the `src` directory to generate `RRI_UNST.exe`.


```bash
make
```


## Run
必要な入力ファイルを準備した後、以下のコマンドを実行すると計算が始まります:  
After preparing the necessary input files, run the following command to start the calculation:

```bash
./RRI_UNST.exe
```

## Coupled calculation with the RRI model
RRIモデルとの連成計算を行う場合、国立研究開発法人土木研究所よりRRIモデル（ver.1.4.2.7）を入手し、以下の手順に従って`RRI.f90`にコードを追加してください。以下の手順の (Line. 行数) で記載の行数は、`RRI.f90` (ver.1.4.2.7) の行数を表します。
GNU patch コマンドが利用可能であれば、パッチ ファイル`modify_rri.patch`を用いて自動的に修正を適用できます。下記のコマンドで、`RRI.f90`を上書きします。    
  
If you want to perform a coupled calculation with the RRI model, obtain the RRI model (ver.1.4.2.7) from the Public Works Research Institute and add the code to `RRI.f90` according to the following procedure. The line numbers in (Line.) in the following procedure refer to the line numbers in `RRI.f90` (ver.1.4.2.7). If you have the GNU patch command available, you can apply the modifications automatically using the patch file `modify_rri.patch`. Overwrite `RRI.f90` with the command below.
  
  
```bash
patch RRI.f90 < modify_rri.patch
```
  

RRIモデルは以下のURLよりダウンロードできます。  
The RRI model can be downloaded from the following URL.   
https://www.pwri.go.jp/icharm/research/rri/index_j.html
 
  
    
### Additional Code for RRI.f90

#### 1. Load Module after line.8 of RRI.f90 (ver.1.4.2.7).
```diff
 use globals
+use unst_globals_mod
```

#### 2. Add Variables after line.67 of RRI.f90 (ver.1.4.2.7).

```diff
 character*256 ofile_ro, outdir1, outdir2
 integer kk, l
+integer uflg
+uflg = 0
```

#### 3. End of file condition after line.508 of RRI.f90 (ver.1.4.2.7).
This modification is for compiliation with gfortran.

```diff
  if( ios.lt.0 ) exit
+ if( ios.ne.0 ) exit
  tt = tt + 1
```


#### 4. Rain Setting after line.527 of RRI.f90 (ver.1.4.2.7).

```diff
 qp = qp / 3600.d0 / 1000.d0
+dtrain = t_rain(1) - t_rain(0)
 do j = 1, nx
```


#### 5. Call UNST2D  after line.577 of RRI.f90 (ver.1.4.2.7).

```diff
  close(11)
 endif

+call unst_rdat(ny_rain, qp, nx_rain, tt_max_rain, &
+    ny, nx, lasth, dt, &
+    xllcorner_rain, yllcorner_rain, cellsize_rain_x, cellsize_rain_y, &
+    xllcorner, yllcorner, cellsize)
+if(plantFN==1) call plantFNdat
+if(plantDa==1) call plantDadat
+if(paddydam==1) call paddydat
+if(drainarea==1) call draindat
+if(morid==1) call moriddat
+
+call unst_initiald
+if(paddydam==1) call paddyinitiald
+if(drainarea==1) call draininitiald
+
+if(str_type == 0) then
+  phi(61) = 0.20d0
+  phi(63) = 0.20d0
+elseif(str_type == 1) then
+  phi(61) = 0.05d0
+  phi(63) = 0.20d0
+endif
+
+call diskwrite
+call dispwrite
+
 ! For TSAS Output (Initial Condition)
 call sub_slo_ij2idx( hs, hs_idx )
 call sub_riv_ij2idx( hr, hr_idx )
```

#### 6. Timestep Setting after line.595 of RRI.f90 (ver.1.4.2.7).

```diff
 out_next = nint(out_dt)
+timmax = dble(dt) * dble(calldt)
 tt = 0
```

#### 7. UNST2D > RRI after line.1195 of RRI.f90 (ver.1.4.2.7)

```diff
   endif
  endif

+ call unst_qin(ny, nx, i4, qr_ave, qs_ave, hr, hs, area, time, uflg)
+ if(mod(time, timmax) == 0 .and. uflg==0) call predispwrite(time)
+ if(mod(time, timmax) == 0 .and. uflg==0) call prediskwrite(time)
+ if(mod(time, timmax) == 0 .and. uflg==1) call UNST(ny, nx, domain, riv, area, time, hs, hr)
+
  ! check water balance
  if(mod(t, 1).eq.0) then
```

#### 8. File closing and memory release processing after line.1205 of RRI.f90 (ver.1.4.2.7)

```diff
 enddo

+call wrhmax
+
+write(*, 1999) time
+1999 format('      - normal end -  time=', f8.0)
+
+if(str_type==1) then
+  close(10077)
+  close(10078)
+endif
+close(10091)
+close(10093)
+close(10094)
+close(10095)
+close(10096)
+close(10097)
+close(100101)
+close(100102)
+close(100103)
+
+if(paddydam==1) then
+  close(10098)
+  close(10099)
+  close(100100)
+endif
+
+deallocate(baseo, dnox, dnoy, smesh, scv, rthl, ux, uy, xmesh, ymesh, rtuv_x, rtuv_y)
+deallocate(limesh, linode, inf, ko, menode, melink, inl, qin, lkyokai_dx, lkyokai_dy)
+deallocate(unsth, ho, hl, hmax, uummax, vvmmax)
+deallocate(um, umo, umm, uu, vn, vno, vnm, vv)
+deallocate(mn, rnof, lambda, rbeta, umbeta, vnbeta)
+deallocate(uum, vvm, lhan, lhano, qr_sum, rnx, dl)
+deallocate(unstc, co, str, stro, strmx)
+if(plantFN==1) deallocate(plantF_array, plantN_array)
+if(plantDa==1) deallocate(plant_D_array, plant_a_array, dk_val)
+if(paddydam==1) deallocate(paddyid, pqout_idx, pdrain, min_pmeshid, device)
+if(paddydam==1) deallocate(orifice_num, min_dist, psmesh, dr_dist, dhp, phid)
+if(paddydam==1) deallocate(paddy_q, pqh, drain2phidx)
+if(drainarea==1) deallocate(inf_dr, drp, drr, drr_dist, dhj)
+if(dsmesh==1) deallocate(dsdt, dsinf, dsupper, dsfilter2)
+if(ga==1) deallocate(genes)
+
 !pause

 end program RRI
```

## License

Licensed under the [MIT](https://github.com/TK-Labo/RRI2UNST2D/blob/main/LICENSE) license.  
  
Copyright (c) 2025 K.Kawaike & TK Labo

## Coded by

Kenji Kawaike & TK Labo

TK Labo Members:

- Koki Yamamura 2021-2025
- Daiki Baba 2022-
- Shunji Nishino 2023-
- Kentaro Taki 2017-
