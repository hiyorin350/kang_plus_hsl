import os
import numpy as np
import cv2
from conversions import *

def angle_to_normal_vector(angle):
    """
    2次元空間において、指定された角度での直線の法線ベクトルを計算する。
    
    :param angle: 直線がX軸と成す角度（度単位）
    :return: 法線ベクトル（numpy配列）
    """
    # 角度をラジアンに変換
    angle_rad = np.radians(angle)
    
    # 直線の法線ベクトルの計算
    # 直線が成す角度に90度を加える（垂直な方向）
    nx = np.cos(angle_rad + np.pi / 2)
    ny = np.sin(angle_rad + np.pi / 2)
    
    return np.array([nx, ny, 0])

def project_pixels_to_color_plane(image, u):
    """
    射影された画像を返す関数。

    :param image: 入力画像（直交HSL 色空間）
    :param u: 色平面の法線ベクトル
    :return: 射影された画像
    """
    # 画像の形状を取得
    height, width, _ = image.shape

    # 射影された画像を格納するための配列を初期化
    projected_image = np.zeros_like(image)

    # 各画素に対して射影を行う
    for i in range(height):
        for j in range(width):
            # 画素の色ベクトルを取得
            color_vector = image[i, j, :]

            # 色ベクトルを色平面に射影
            projected_vector = color_vector - np.dot(color_vector, u) * u

            # 射影された色ベクトルを保存
            projected_image[i, j, :] = projected_vector

    return projected_image

def process_images_in_folder(input_folder, output_folder):
    # フォルダ内の全ての画像ファイルを取得
    for filename in os.listdir(input_folder):
        if filename.endswith('.ppm'):
            # 画像の読み込み
            image_path = os.path.join(input_folder, filename)
            image = cv2.imread(image_path)

            # HSL変換と直交HSL変換
            hsl_image = rgb_to_hsl(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))

            hsl_cartesian = hsl_to_cartesian(hsl_image)

            # 画像の形状を取得
            height, width, _ = image.shape

            # 2色覚平面の法線ベクトル
            u = angle_to_normal_vector(50.19789)

            # 射影処理
            projected_image = project_pixels_to_color_plane(hsl_cartesian, u)

            # 射影結果をRGBに変換
            img_out_hsl = cartesian_to_hsl(projected_image)
            img_out_rgb = hsl_to_rgb(img_out_hsl)

            # BGRに変換して保存
            img_out_bgr = cv2.cvtColor(img_out_rgb, cv2.COLOR_RGB2BGR)
            output_path = os.path.join(output_folder, f'projected_{filename}')
            cv2.imwrite(output_path, img_out_bgr)

            print(f"Processed and saved: {output_path}")

# 入力フォルダと出力フォルダの指定
input_folder = 'images/chart26/rotates/'
output_folder = 'images/chart26/projections/'

# 出力フォルダが存在しない場合は作成
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# フォルダ内の全ての画像に対して処理を実行
process_images_in_folder(input_folder, output_folder)
