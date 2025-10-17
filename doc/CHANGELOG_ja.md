# Changelog

## [1.1.0] - comming soon
### Add
- RRIとUNST2Dの接続方法の追加
  - 水深接続の追加
- テストデータの公開
- Windows用ビルドバッチファイル(make.bat)の追加

### Fixed
- UNST2Dの呼び出しメソッドの修正(接続メッシュ下流からの降雨イベントに対応)
  - RRI-UNST2D接続メッシュに流量がない場合UNST2Dが呼び出されない現象を修正

### Change
- cntl.datをUNST2D_cntl.datに名称変更(UNST2D単体モデルのcntl.datとの差別化)
- cntl.datのフォーマットの変更
  - ocpy(占有率), beta(通過率)の削除
  - RRI-UNST2Dの接続タイプ指定(RRIへのフィードバックの有無)の追加
  - lkyokaiq.datの出力先をcntl.datで指定するように変更
- 現状未反映となっている占有率の指定をinf.datの横に記載するよう変更
- 通過率の指定を線盛土オプションと統合
- マニュアルの更新(ver.1.1.0)

## [1.0.5] - 2025-10-17
### Add
- 変更記録ファイル(CHANGELOG_ja.md / CHANGELOG_en.md)の追加
- (UNST_***.f90共通) 負の水深補正に使用した累計流量の出力を追加
- (UNST_***.f90共通) 排水された累計流量の出力を追加

### Fixed
- (UNST_Sub.f90) 田んぼダム有効時(paddydam==1)の水収支処理の不具合を修正
- (UNST_Elements.f90) UNST2DからRRIへ受け渡す平均水深処理の不具合を修正

### Change
- ソースコードの変更
  - (UNST_***.f90共通) コメントの充実化
  - (UNST_***.f90共通) str(貯留)機能の削除
  - (UNST_***.f90共通) 排水処理の各種変数名を変更
  - (UNST_Main.f90) gotoによるループをdo while文によるループに変更
  - (UNST_Main.f90) module化に応じたuseの追加
  - (UNST_Main.f90) interfaceによりrri2unstを明示的に
  - (UNST_Read.f90) 全体をモジュール化(module unst_read)
  - (UNST_Read.f90) 装置番号をnewunitに変更
  - (UNST_Read.f90) UNST2Dのインプット諸元をターミナル上に出力するよう変更
  - (UNST_Read.f90) 排水処理(dsmesh)の読み取りをsubroutine dsmesh.datに分離
  - (UNST_Sub.f90) gotoの削減
  - (UNST_Sub.f90) 全体ををモジュール化(module unst_cal_sub)
  - (UNST_Sub.f90) subroutine sumqaを削除
  - (UNST_Sub.f90) entry suisinをsubroutine化(subroutine suisin)
  - (UNST_Write.f90) 全体をmodule化(unst_wrfile, unst_prewfile)
  - (UNST_Write.f90) UNST_Sub.f90のsubroutine sumqaの内容をsubroutine dispwriteに移植
  - (UNST_Write.f90) entry paddywriteをsubroutine化(subroutine paddywrite)
  - (UNST_Write.f90) 出力フォーマットの一部を動的に変更
  - (UNST_Mod.f90) 変数の整理
  - (modify_rri.patch/RRI_UNST.f90) module化に応じたuseの追加

- Makefileのビルド対象からUNST_Break.f90を除外
- マニュアルの更新(ver.1.0.5)

## [1.0.0] - 2025-07-07
### Publish
- RRI-UNST2Dをgithub上に公開