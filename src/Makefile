#==============================================================================
# RRI-UNST連携モデル用Makefile
#
# 使い方:
#   make        : プログラムをコンパイルして実行ファイルを生成
#   make clean  : 中間ファイルと実行ファイルを削除
#
# 最終更新: 2023-12-01
#==============================================================================

#------------------------------------------------------------------------------
# コンパイラと最適化オプションの設定
#------------------------------------------------------------------------------
# Fortranコンパイラの指定
FC = gfortran

# 通常のコンパイルオプション（最適化レベル2、OpenMP有効）
FFLAGS = -g -O2 -cpp -ffree-line-length-none -fopenmp

# デバッグモード用オプション
# アンコメントして使用する場合は上のFFLAGSをコメントアウト
# FFLAGS = -g -O0 -cpp -ffree-line-length-none -fopenmp -fbacktrace -fcheck=all

# リンクオプションの指定
LDFLAGS = -fopenmp

#------------------------------------------------------------------------------
# ソースファイルとオブジェクトファイルの設定
#------------------------------------------------------------------------------
# ソースファイルの拡張子
FSUFFIX = .f90

# ソースファイルの取得（ワイルドカードで自動検出）
EXCLUDE = RRI_Break.f90  # RRI_Break.f90 は除きます。
ALLSOURCES = $(wildcard *$(FSUFFIX))
SOURCES := $(filter-out $(EXCLUDE), $(ALLSOURCES))

# オブジェクトファイルの拡張子
OSUFFIX = .o

# オブジェクトファイルの取得（ソースファイルの拡張子を置換）
OBJECTS = $(patsubst %$(FSUFFIX),%$(OSUFFIX),$(SOURCES))

# 実行ファイル名
TARGET = RRI_UNST.exe

#------------------------------------------------------------------------------
# ビルドルール
#------------------------------------------------------------------------------
# デフォルトターゲット：クリーン→ビルド→後処理
all: clean $(TARGET) post_process

# 実行ファイルの生成
$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(FC) $(LDFLAGS) -o $@ $^ || (rm -f $(OBJECTS) $(SOURCES:.f90=.d); exit 1)
	@echo "Build successful!"

# オブジェクトファイルの生成
%$(OSUFFIX): %$(FSUFFIX)
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -c $< || (rm -f $@ $(patsubst %$(FSUFFIX),%.d,$<); exit 1)

# 依存関係ファイルの生成
%.d: %$(FSUFFIX)
	@echo "Generating dependencies for $<..."
	$(FC) -cpp -MM $< > $@ || (rm -f $@; exit 1)

# 依存関係ファイルのインクルード
-include $(SOURCES:.f90=.d)

#------------------------------------------------------------------------------
# クリーンアップルール
#------------------------------------------------------------------------------
# 中間ファイルと実行ファイルの削除
.PHONY: clean
clean:
	@echo "Cleaning build files..."
	rm -f $(OBJECTS) $(TARGET) $(SOURCES:.f90=.d) $(SOURCES:.f90=.mod)
	@echo "Clean completed."

# ビルド後の中間ファイル削除（モジュールファイルとオブジェクトファイル）
.PHONY: post_process
post_process:
	@echo "Removing intermediate files..."
	rm -f $(OBJECTS) $(SOURCES:.f90=.d) $(SOURCES:.f90=.mod)
	@echo "Post-processing completed."

#------------------------------------------------------------------------------
# ヘルプ
#------------------------------------------------------------------------------
.PHONY: help
help:
	@echo "RRI-UNST連携モデル Makefile"
	@echo ""
	@echo "使用可能なターゲット:"
	@echo "  make       - プログラムをビルドして実行ファイルを生成"
	@echo "  make clean - 中間ファイルと実行ファイルを削除"
	@echo "  make help  - このヘルプメッセージを表示"
	@echo ""
	@echo "コンパイルオプションを変更する場合はMakefileを編集してください。"