#==============================================================================
# RRI-UNST連携モデル用Makefile
#
# 使い方:
#   make        : プログラムをコンパイルして実行ファイルを生成
#   make clean  : 中間ファイルと実行ファイルを削除
#
# 最終更新:
#==============================================================================

#------------------------------------------------------------------------------
# コンパイラと最適化オプションの設定
#------------------------------------------------------------------------------
# Fortranコンパイラの指定
FC = gfortran
#FC = ifx

# 通常のコンパイルオプション（最適化レベル2、OpenMP有効）
# gfortran 用
FFLAGS = -O2 -fopenmp -cpp -ffree-line-length-none
# ifx 用
#FFLAGS = -O2 -qopenmp

# # デバッグモード用オプション
# # gfortran 用
# FFLAGS = -g -O0 -cpp -ffree-line-length-none -fopenmp -fbacktrace -fcheck=all --warn-all
# # ifx 用
# FFLAGS = -g -traceback -O0 -check all -traceback -fpe0 -warn all --extend-source

# UNST と RRI のソースファイルのディレクトリ
UNST_DIR = unst_src
RRI_DIR = 1.4.2.7

#------------------------------------------------------------------------------
# ソースファイルとオブジェクトファイルの設定
#------------------------------------------------------------------------------
# ソースファイルの拡張子
FSUFFIX = .f90

# オブジェクトファイルの拡張子
OSUFFIX = .o

# ビルド用ディレクトリ (中間ファイルの出力先)
BUILD_DIR = build

# UNST のオブジェクトファイルの取得
# モジュール オブジェクト
UNST_MODULE = ${BUILD_DIR}/UNST_Mod.o ${BUILD_DIR}/UNST_Write.o
# モジュール以外
UNST_SOURCES = $(filter-out $(UNST_DIR)/UNST_Mod.f90 $(UNST_DIR)/UNST_Write.f90, \
 $(wildcard $(UNST_DIR)/*$(FSUFFIX)))
UNST_OBJECTS = $(UNST_MODULE)
UNST_OBJECTS += $(patsubst $(UNST_DIR)/%$(FSUFFIX), $(BUILD_DIR)/%$(OSUFFIX), $(UNST_SOURCES))

# RRI のオブジェクトファイルの取得
# モジュールオブジェクト
RRI_MODULE = ${BUILD_DIR}/RRI_Mod.o ${BUILD_DIR}/RRI_Mod2.o \
 ${BUILD_DIR}/RRI_Mod_Dam.o ${BUILD_DIR}/RRI_Mod_Tecout.o
# モジュール以外
RRI_SOURCES = $(filter-out $(RRI_DIR)/RRI_Mod.f90 $(RRI_DIR)/RRI_Mod2.f90 \
 $(RRI_DIR)/RRI_Mod_Dam.f90 $(RRI_DIR)/RRI_Mod_Tecout.f90 \
 $(RRI_DIR)/RRI_Break.f90 $(RRI_DIR)/RRI.f90 $(RRI_DIR)/RRI_UNST.f90, $(wildcard $(RRI_DIR)/*.f90))
RRI_OBJECTS = $(RRI_MODULE)
RRI_OBJECTS += $(patsubst $(RRI_DIR)/%$(FSUFFIX), $(BUILD_DIR)/%$(OSUFFIX), $(RRI_SOURCES))

# メインプログラムのオブジェクト
MAIN_OBJECT = ${BUILD_DIR}/RRI_UNST$(OSUFFIX)

# オブジェクト全体
OBJECTS = $(UNST_OBJECTS) $(RRI_OBJECTS) $(MAIN_OBJECT)

# ソースコード全体
SOURCES = $(wildcard $(UNST_DIR)/*$(FSUFFIX)) $(wildcard $(RRI_DIR)/*$(FSUFFIX))

# 実行ファイル名
TARGET = RRI_UNST.exe

#------------------------------------------------------------------------------
# ビルドルール
#------------------------------------------------------------------------------
all: $(BUILD_DIR) $(TARGET)

# リンクの指定
LDFLAGS = -L$(BUILD_DIR)

# インクルードの指定
INCLDIR = $(BUILD_DIR)

# .mod の出力先の指定
ifeq ($(FC),ifx)
  MODFLAG = -module $(BUILD_DIR)
else
  MODFLAG = -J$(BUILD_DIR)
endif

# 依存関係の生成ルール (gfortran と ifx に対応)
ifeq ($(findstring gfortran,$(FC)),gfortran)
    DEPRULE = $(FC) $(FFLAGS) -c -I$(INCLDIR) $(MODFLAG) -o $@ $< ; \
              $(FC) $(FFLAGS) -cpp -M $< -I$(INCLDIR) $(MODFLAG) > $(basename $@).d
else ifeq ($(findstring ifx,$(FC)),ifx)
    DEPRULE = $(FC) $(FFLAGS) -c -cpp -o $(BUILD_DIR)/$*$(OSUFFIX) \
	-gen-dep=$@ -gen-depformat=make -I$(INCLDIR) $(MODFLAG) $<
endif

# 依存関係ファイルの生成
$(BUILD_DIR)/%.d: $(UNST_DIR)/%.f90 | $(BUILD_DIR)
	@echo "Generating dependencies for $<..."
	$(DEPRULE)

$(BUILD_DIR)/%.d: $(RRI_DIR)/%.f90 | $(BUILD_DIR)
	@echo "Generating dependencies for $<..."
	$(DEPRULE)

# オブジェクトファイルの生成
$(BUILD_DIR)/%.o: $(UNST_DIR)/%.f90 $(BUILD_DIR)/%.d | $(BUILD_DIR)
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -c $< -I$(INCLDIR) $(MODFLAG) -o $@

$(BUILD_DIR)/%.o: $(RRI_DIR)/%.f90 $(BUILD_DIR)/%.d | $(BUILD_DIR)
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -c $< -I$(INCLDIR) $(MODFLAG) -o $@

# ビルドディレクトリの作成
$(BUILD_DIR):
	@echo "Making build directory ..."
	mkdir -p $(BUILD_DIR)

# 実行ファイルの生成
$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(FC) -o $@ $^ $(FFLAGS) $(LDFLAGS)
	@echo "Build successful!"

#------------------------------------------------------------------------------
# クリーンアップルール
#------------------------------------------------------------------------------
# 中間ファイルと実行ファイルの削除
clean:
	test -n "$(BUILD_DIR)" && test "$(BUILD_DIR)" != "/" && rm -rf $(BUILD_DIR)
	test -n "$(TARGET)" && test "$(TARGET)" != "/" && rm -rf $(TARGET)

# 中間ファイルの削除
post_process:
	test -n "$(BUILD_DIR)" && test "$(BUILD_DIR)" != "/" && rm -rf $(BUILD_DIR)

#------------------------------------------------------------------------------
# RRI_UNST.f90 の作成
#------------------------------------------------------------------------------
RRI_ORIGINAL = $(RRI_DIR)/RRI.f90
RRI_UNST = $(RRI_DIR)/RRI_UNST.f90
PATCH_FILE = modify_rri.patch

rri_unst_main:
	@echo "Applying patch to $(RRI_ORIGINAL) to create $(RRI_UNST)."
	@cat $(RRI_ORIGINAL) > $(RRI_UNST).tmp
	@patch $(RRI_UNST).tmp < $(PATCH_FILE)
	@mv $(RRI_UNST).tmp $(RRI_UNST)
	@echo "$(RRI_UNST) is created successfully."


#------------------------------------------------------------------------------
# ヘルプ
#------------------------------------------------------------------------------
help:
	@echo "RRI-UNST連携モデル Makefile"
	@echo ""
	@echo "使用可能なターゲット:"
	@echo "  make       - プログラムをビルドして実行ファイルを生成"
	@echo "  make clean - 中間ファイルと実行ファイルを削除"
	@echo "  make post_process - 中間ファイルを削除"
	@echo "  make rri_unst_main - RRI_UNST.f90 を作成"
	@echo "  make help  - このヘルプメッセージを表示"
	@echo ""
	@echo "コンパイルオプションを変更する場合はMakefileを編集してください。"