# result_to_csv by d.baba  TK Labo
import pandas as pd
import re
import numpy as np

# input
hdat = "../out/unst/h.dat"
udat = "../out/unst/uum.dat"
vdat = "../out/unst/vvm.dat"
hmaxdat = "../out/unst/hmax.dat"
uummaxdat = "../out/unst/uummax.dat"
vvmmaxdat = "../out/unst/vvmmax.dat"
epsg = 2448
# output
## csv
out_hcsv = "../out/unst/h.csv"
out_ucsv = "../out/unst/uum.csv"
out_vcsv = "../out/unst/vvm.csv"
speed_csv = "../out/unst/speed.csv"
hydrodynamicforce_csv = "../out/unst/hdf.csv"
## max csv
hmax_csv = "../out/unst/hmax.csv"
uummax_csv = "../out/unst/uummax.csv"
vvmmax_csv = "../out/unst/vvmmax.csv"
speedmax_csv = "../out/unst/speedmax.csv"
hydrodynamicforcemax_csv = "../out/unst/hdfmax.csv"

# process
def process_line(line:str, width:int=15):
    """
    Processes a line of data:
    - Splits the line into fixed-width segments.
    - Replaces `***************` or `Infinity` or `NaN` with 9999.999.
    - Converts segments to floats.
    
    Parameters
    ----------
    line: str
        Input line containing the data.
    width: int
        Width of each data segment(default:15)

    Returns
    ------
    list: list
        List of processed float values.
    """
    if len(line) != 150:
        padded_line = " " * (width - (len(line) % width)) + line
    else:
        padded_line = line
    processed_data = []
    for i in range(0, len(padded_line), width):
        segment = padded_line[i:i+width]  # Extract fixed-width segment
        if segment == "***************" or segment == "       Infinity" or segment =="            Nan":
            processed_data.append(9999.999)  # Replace `***************` with 9999.999
        else:
            processed_data.append(float(segment))  # Convert to float
    return processed_data

def unst_output(file_path:str):
    '''
    DAT to CSV
    
    Parameters
    ----------
    file_path: str
        output path(h.dat/uum.dat/vvm.dat)

    Returns
    ------
    df: DataFrame
        output df
    
    '''

    with open(file_path, 'r') as file:
        lines = file.readlines()

    current_label = None
    data_dict = {}
    current_data = []

    for line in lines:
        line = line.strip()
        # ラベル行を検出
        if line.startswith('time='):
            # 新しいラベルが検出されたら、前のラベルのデータを保存
            if current_label is not None and current_data:
                data_dict[current_label] = current_data
                current_data = []
            # ラベルを取得
            current_label = re.search(r"time=\s*([\d.]+\(s\))", line).group(1)
        elif current_label:
            # データ行を読み込む
            #current_data.extend(map(float, line.split()))
            current_data.extend(process_line(line))
    
    # 最後のラベルのデータを保存
    if current_label is not None and current_data:
        data_dict[current_label] = current_data

    # dfへ変換
    df = pd.DataFrame(data_dict)

    # mesh id列の追加
    # df.insert(0, 'id', df.index+1)

    return df

# データの読み込み
# 時系列データの読み込み
print("reading h.dat ...")
hdat_df = unst_output(hdat)
hdat_df2 = hdat_df.copy()
hdat_df2.insert(0, 'id', hdat_df2.index+1)
hdat_df2.to_csv(out_hcsv, index=False)

print("reading uum.dat ...")
udat_df = unst_output(udat)
udat_df2 = udat_df.copy()
udat_df2.insert(0, 'id', udat_df2.index+1)
udat_df2.to_csv(out_ucsv, index=False)

print("reading vvm.dat ...")
vdat_df = unst_output(vdat)
vdat_df2 = vdat_df.copy()
vdat_df2.insert(0, 'id', vdat_df2.index+1)
vdat_df2.to_csv(out_vcsv, index=False)

print("meke speed ...")
speed_df = (udat_df ** 2 + vdat_df ** 2).map(np.sqrt)
speed_df2 = speed_df.copy()
speed_df2.insert(0, 'id', speed_df2.index+1)
speed_df2.to_csv(speed_csv, index=False)

print("make Hydrodynamic force ...")
hdf_df = (udat_df ** 2 + vdat_df ** 2) * hdat_df
hdf_df.to_csv(hydrodynamicforce_csv, index=False)

# 最大値データの読み込み
##
print("reading hmax.dat ...")
hmax_df = unst_output(hmaxdat)
hmax_df.insert(0, "id", hmax_df.index+1)
hmax_df.to_csv(hmax_csv, index=False)
print("reading uummax.dat ...")
uummax_df = unst_output(uummaxdat)
uummax_df.insert(0, "id", uummax_df.index+1)
uummax_df.to_csv(uummax_csv, index=False)
print("reading vvmmax.dat ...")
vvmmax_df = unst_output(vvmmaxdat)
vvmmax_df.insert(0, "id", vvmmax_df.index+1)
vvmmax_df.to_csv(vvmmax_csv, index=False)
print("calculate speed max ...")
speed_max = speed_df.max(axis=1)
speed_max_df = speed_max.to_frame(name="Maxspeed")
speed_max_df.insert(0, "id", speed_max_df.index+1)
speed_max_df.to_csv(speedmax_csv, index=False)
print("calculate hydrodynamic force max ...")
hdf_max = hdf_df.max(axis=1)
hdf_max_df = hdf_max.to_frame(name="Maxfp")
hdf_max_df.insert(0, "id", hdf_max_df.index+1)
hdf_max_df.to_csv(hydrodynamicforcemax_csv, index=False)