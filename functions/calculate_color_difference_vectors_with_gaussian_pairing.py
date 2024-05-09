import numpy as np
import cv2
from conversions import *

def calculate_color_difference_vectors_with_gaussian_pairing(image):#TODO この段階で重みを計算
    """
    ガウシアンペアリングを使用して画像内の各ピクセルの色差ベクトルを計算します。

    :param image: CIE L*a*b* 色空間の入力画像。
    :param sigma: ペアリングに使用されるガウス分布の標準偏差。
    :return: 色差ベクトルの配列。重み関数w
    sigmaは分散ではなく標準偏差であることに注意
    """
    height, width, _ = image.shape
    N = height * width
    Xl = np.zeros((N, 3))  # 色差行列の初期化
    w = np.zeros(N)

    sigma = np.sqrt((2 / np.pi) * np.sqrt(2 * min(height, width)))# 近隣ピクセルの参照に使用

    sigma_weight = 5 #明度差重みのパラメータ

    hsl_image = rgb_to_hsl(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))

    mhsl_image = hsl_to_mhsl(hsl_image)

    cartesian_mhsl = hsl_to_cartesian(mhsl_image)

    flat_mhsl_image = cartesian_mhsl.reshape(N, 3)

    # ガウス分布を使用してランダムなペアを生成
    for i in range(N):
        # ピクセルの座標を計算
        x, y = i % width, i // width

        # ガウス分布に基づいてランダムなオフセットを生成
        offset_x, offset_y = np.random.normal(0, sigma, 2)
        nx, ny = int(x + offset_x), int(y + offset_y)

        # 新しい座標が画像の範囲内に収まるようにする
        nx = max(0, min(nx, width - 1))
        ny = max(0, min(ny, height - 1))

        # 色差ベクトルを計算
        neighbor_index = ny * width + nx
        Xl[i, 0] = flat_mhsl_image[i, 0] - flat_mhsl_image[neighbor_index, 0]#直交座標HSLのx
        Xl[i, 1] = flat_mhsl_image[i, 1] - flat_mhsl_image[neighbor_index, 1]#y
        Xl[i, 2] = flat_mhsl_image[i, 2] - flat_mhsl_image[neighbor_index, 2]#z(=L)

        w[i] = np.exp(-np.square(flat_mhsl_image[i, 2] - flat_mhsl_image[neighbor_index, 2]) / 2 * np.square(sigma_weight))
        # print(w[i])
    
    return Xl, w
