# movie.py  by d.baba TK Labo
# UNST2Dの結果を面的可視化するツール
import pandas as pd
import re
import geopandas as gpd
from shapely import wkt
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# ffmegが動かない場合は有効にしてffmegのbinフォルダーのパスを指定
#plt.rcParams["animation.ffmpeg_path"] = r"C:\ffmpeg\bin\ffmpeg.exe" 
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
#import tqdm
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap

# 読み込み
skip_dattocsv = 0  # h.csvが既にある場合有効(1)
dat = "../out/unst/h.dat"  # 可視化したい出力結果(h.dat or uum.dat or vvm.dat)
shp_file = "../MK_mesh/pred/CSV/wkt.csv"  # 格子データ(shp/csv(wkt) など)
epsg = 6674  # 座標系
# output
## dat → csv
out_csv = "../out/unst/h.csv"  # h.csvのパス
## movie option
# 色の指定など
set_colortype = 0  # 0:select 1:浸水想定区域図マニュアル(連続) 2:浸水想定区域図マニュアル(離散)
camp_c = "Blues"  # select時の色調
n_bins = 100  # 一度に100階調を作成
vmin = 0.0  # レンジ(最小値)
vmax = 0.4  # レンジ(最大値)
data_label = "depth(m)"
interval = 100
output_name = "unst_animation.mp4"  # 出力パス

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
        unst出力結果(h.datなど)

    Returns
    ------
    df: DataFrame
        csv形式で整理されたunst出力結果
    
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
            print(f'read data : {line}')
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
if skip_dattocsv != 1:
    print('== UNST RESULT MOVIE ver1.0.0 ==')
    dat_df = unst_output(dat)
else:
    dat_df = pd.read_csv(out_csv)

gdf = gpd.read_file(shp_file)
gdf = gdf.sort_values(by='meid', key=lambda x: x.astype(int)).reset_index(drop=True)

dataset = pd.concat([gdf, dat_df], axis=1)#横方向に結合
dataset['geometry'] = dataset['WKT'].apply(wkt.loads)
dataset = dataset.drop(columns='WKT')
dataset = gpd.GeoDataFrame(dataset, geometry='geometry')
dataset.set_crs(epsg=epsg, inplace=True)

# 時系列データ列を取得（"(s)"が含まれる列を選択）
time_columns = [col for col in dataset.columns if "(s)" in col]

# 時系列データの値を取得
time_values = np.array(dataset[time_columns])

# カラーマップ設定
if set_colortype == 0:
    custom_cmap = plt.get_cmap(camp_c)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
elif set_colortype == 1:
    colors = [(1, 1, 0), (0, 1, 1), (1, 0, 1)]  # Yellow, Cyan, Magenta
    custom_cmap = LinearSegmentedColormap.from_list('custom', colors, N=n_bins)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
elif set_colortype == 2:
    colors = [(247/255,245/255,169/255),
              (255/255,216/255,192/255),
              (255/255,183/255,183/255), 
              (255/255,145/255,145/255),
              (242/255,133/255,201/255),
              (220/255,122/255,220/255)
              ]
    select_boundary = [0, 1.75, 4.0, 7.5, 15.0, 20]  # 0〜5の範囲で境界値を設定
    positions = [sb/max(select_boundary) for sb in select_boundary]
    custom_cmap = LinearSegmentedColormap.from_list('custom', list(zip(positions, colors)), N=n_bins)
    norm = plt.Normalize(vmin=0, vmax=20)  # Normalizeの範囲を0〜5に設定
elif set_colortype == 3:
    colors = [(247/255,245/255,169/255),
              (255/255,216/255,192/255),
              (255/255,183/255,183/255), 
              (255/255,145/255,145/255),
              (242/255,133/255,201/255),
              (220/255,122/255,220/255)
              ]  # 各区間の色
    # BoundaryNormを使用して離散的な色分けを実現
    select_boundary = [0.5, 3.0, 5.0, 10.0, 20.0, 9999.0]  # 水深の境界値
    norm = BoundaryNorm(select_boundary, len(colors))
    custom_cmap = LinearSegmentedColormap.from_list('custom', colors, N=len(colors))

# アニメーション作成
fig, ax = plt.subplots(figsize=(8, 8))
# dataset.plot(ax=ax, color="white", edgecolor="black")  # 初期プロット
minx, miny, maxx, maxy = dataset.total_bounds  # ジオメトリの最小・最大範囲
buffer = 0.1  # 5%の余裕を持たせる
ax.set_xlim(minx - (maxx - minx)*buffer , maxx + (maxx - minx)*buffer)
ax.set_ylim(miny - (maxy - miny)*buffer, maxy + (maxy - miny)*buffer)
ax.set_aspect('equal')  # アスペクト比固定

# 初期パッチを設定
polygons = [Polygon(np.array(geom.exterior.coords)) for geom in dataset.geometry]
patches = PatchCollection(polygons, cmap=custom_cmap, norm=norm)
patches.set_array(time_values[:, 0])  # 初期値を設定
ax.add_collection(patches)

# 凡例を追加
sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax)
cbar.set_label(data_label, rotation=270, labelpad=20)

#progress_bar = tqdm.tqdm(total=len(time_columns), desc="Progress", unit="frame")

# アニメーションの更新関数
def update(frame):
    patches.set_array(time_values[:, frame])
    ax.set_title(f"Time: {time_columns[frame]}", loc='center')
    #progress_bar.update(1)  # 進捗バーを更新

ani = FuncAnimation(fig, update, frames=len(time_columns), interval=interval)
# アニメーション保存または表示

ani.save(output_name, writer="ffmpeg")  # MP4に保存