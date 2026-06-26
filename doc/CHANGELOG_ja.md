# Changelog

## [1.0.6] - 2026-07-07
### Fixed
- (UNST_***.f90共通) 一部の変数の配列アクセスを列優先に最適化
- (UNST_Read.f90) RRI->UNST2Dの接続方法(メッシュ指定)の場合の処理設定の一部を修正
- (UNST_Read.f90) 植生抵抗のインプットを追加，格子ごとに設定可能に変更
- (UNST_Initial.f90) 辺の長さの算出式の誤りを修正
- (UNST_Sub.f90) 植生抵抗(樹木)の正負の誤りを修正
- (UNST_Riv.f90) 一次元河道のロジック修正(参照流量の修正)

### Add
- 開発中
  - (UNST_***.f90共通) 倒伏を考慮した植生抵抗の雛形を記述
- 時間発展スキームの変更 (開発中)
  - Leap-frogからRK2(Heun)に変更
    - 計算初期段階等の高速化を実現
    - マニュアル未更新（随時更新中）

### Change
- Leap-frogスキームモデルをversions/v1_0_5_leapfrogへ移動
- (UNST_Read2.f90) d1riv_cntl.datのフォーマットを変更
  

## [1.0.5] - 2025-12-31
### Add
- 変更記録ファイル(CHANGELOG_ja.md / CHANGELOG_en.md)の追加
- 機能追加
  - (UNST_***.f90共通) 負の水深補正に使用した累計流量の出力を追加
  - (UNST_***.f90共通) 排水された累計流量の出力を追加
  - (UNST_***.f90共通) RRI->UNST2Dの接続方法の追加(メッシュ指定)
  - (UNST_***.f90共通) RRI<-UNST2Dの水深フィードバック有無の選択機能追加
  - (UNST_***.f90共通) 一次元河道(フラクショナルステップ)モデルに対応
  - (UNST_Read2.f90) 一次元河道(フラクショナルステップ)モデルの追加
  - (UNST_Riv.f90) 一次元河道(フラクショナルステップ)モデルの追加
  - (UNST_Mod2.f90) 一次元河道(フラクショナルステップ)モデルの追加
- Windows用ビルドバッチファイル(make_win.bat)の追加
- Windows用patch反映pythonスクリプト(apply_patch.py)の追加
- 可視化用サンプルpythonスクリプト(point.py, movie.py, out_to_csv.py)の追加

### Fixed
- (UNST_Sub.f90) 田んぼダム有効時(paddydam==1)の水収支処理の不具合を修正
- (UNST_Sub.f90) 境界流入処理の不具合を修正
- (UNST_Sub.f90) 倒伏を考慮しない植生抵抗の不正確な処理を修正
- (UNST_Cnct.f90) UNST2DからRRIへ受け渡す平均水深処理の不具合を修正

### Change
- cntl.datをUNST2D_cntl.datに名称変更(UNST2D単体モデルのcntl.datとの差別化)
- UNST2D_cntl.datのフォーマットの変更
  - beta行の削除
  - ocpy行の削除
  - kyokaiq.datの出力先を指定するよう変更
  - 倒伏を考慮しない植生抵抗のインプットを追加
- ソースコードの変更
  - (UNST_***.f90共通) コメントの充実化
  - (UNST_***.f90共通) str(貯留)機能の削除
  - (UNST_***.f90共通) 排水処理の各種変数名を変更
  - (UNST_***.f90共通) 倒伏を考慮しない植生抵抗のインプットを追加(抗力係数，有効な樹高)
  - (UNST_Main.f90) gotoによるループをdo while文によるループに変更
  - (UNST_Main.f90) module化に応じたuseの追加
  - (UNST_Main.f90) interfaceによりrri2unstを明示的に
  - (UNST_Read.f90) 全体をモジュール化(module unst_read)
  - (UNST_Read.f90) 装置番号をnewunitに変更
  - (UNST_Read.f90) UNST2Dのインプット諸元をターミナル上に出力するよう変更
  - (UNST_Read.f90) inf.datで占有率lambdaの指定を可能に
  - (UNST_Read.f90) UNST2D領域外の降雨エリアに対応
  - (UNST_Read.f90) UNST2D領域外のRRIエリアに対応
  - (UNST_Read.f90) 排水処理(dsmesh)の読み取りをsubroutine dsmeshに分離
  - (UNST_Read.f90) morido.datで通過率rbetaの指定を可能に
  - (UNST_Initial.f90) 初期条件反映処理の分割(可読性の向上)
  - (UNST_Sub.f90) gotoの削減
  - (UNST_Sub.f90) 全体ををモジュール化(module unst_cal_sub)
  - (UNST_Sub.f90) subroutine sumqaを削除
  - (UNST_Sub.f90) entry suisinをsubroutine化(subroutine suisin)
  - (UNST_Write.f90) 全体をmodule化(unst_wrfile, unst_prewfile)
  - (UNST_Write.f90) UNST_Sub.f90のsubroutine sumqaの内容をsubroutine dispwriteに移植
  - (UNST_Write.f90) entry paddywriteをsubroutine化(subroutine paddywrite)
  - (UNST_Write.f90) 出力フォーマットの一部を動的に変更
  - (UNST_Mod.f90) 各種追加機能に応じた変数の追加・削除
  - (UNST_Mod.f90) 変数の整理
  - (modify_rri.patch/RRI_UNST.f90) module化に応じたuseの追加
  - (modify_rri.patch/RRI_UNST.f90) 各種変更に応じた引数の追加
  - (modify_rri.patch/RRI_UNST.f90) dsmeshの読み込み処理の追加(call dsmesh)
  - (modify_rri.patch/RRI_UNST.f90) 一次元河道関連の処理を追加
- Makefileの更新
  - ビルド対象からUNST_Break.f90を除外
  - 一次元河道関連のソースを追加
- マニュアルの更新(ver.1.0.5)

## [1.0.0] - 2025-07-07
### Publish
- RRI-UNST2Dをgithub上に公開