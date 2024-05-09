import numpy as np
import cv2

def rgb_to_hsl(image):
    # imageは(高さ, 幅, 3)の形状のNumPy配列と仮定
    # dtypeをfloatに変換して計算を行う
    image = image.astype(np.float32) / 255.0

    R = image[:, :, 0]
    G = image[:, :, 1]
    B = image[:, :, 2]

    max_color = np.maximum(np.maximum(R, G), B)
    min_color = np.minimum(np.minimum(R, G), B)

    L = (max_color + min_color) / 2

    delta = max_color - min_color
    S = np.zeros_like(L)
    
    # 彩度の計算
    S[delta != 0] = delta[delta != 0] / (1 - np.abs(2 * L[delta != 0] - 1))

    H = np.zeros_like(L)
    # 色相の計算
    # Rが最大値
    idx = (max_color == R) & (delta != 0)
    H[idx] = 60 * (((G[idx] - B[idx]) / delta[idx]) % 6)

    # Gが最大値
    idx = (max_color == G) & (delta != 0)
    H[idx] = 60 * (((B[idx] - R[idx]) / delta[idx]) + 2)

    # Bが最大値
    idx = (max_color == B) & (delta != 0)
    H[idx] = 60 * (((R[idx] - G[idx]) / delta[idx]) + 4)

    # 彩度と輝度をパーセンテージに変換
    S = S * 100
    L = L * 100

    return np.stack([H, S, L], axis=-1)

def hsl_to_rgb(hsl_image):
    H, S, L = hsl_image[:, :, 0], hsl_image[:, :, 1], hsl_image[:, :, 2]
    H /= 360  # Hを0から1の範囲に正規化
    S /= 100  # Sを0から1の範囲に正規化
    L /= 100  # Lを0から1の範囲に正規化

    def hue_to_rgb(p, q, t):
        # tが0より小さい場合、1を加算
        t[t < 0] += 1
        # tが1より大きい場合、1を減算
        t[t > 1] -= 1
        # t < 1/6の場合
        r = np.copy(p)
        r[t < 1/6] = p[t < 1/6] + (q[t < 1/6] - p[t < 1/6]) * 6 * t[t < 1/6]
        # 1/6 <= t < 1/2の場合
        r[(t >= 1/6) & (t < 1/2)] = q[(t >= 1/6) & (t < 1/2)]
        # 1/2 <= t < 2/3の場合
        r[(t >= 1/2) & (t < 2/3)] = p[(t >= 1/2) & (t < 2/3)] + (q[(t >= 1/2) & (t < 2/3)] - p[(t >= 1/2) & (t < 2/3)]) * (2/3 - t[(t >= 1/2) & (t < 2/3)]) * 6
        # t >= 2/3の場合、rは変更なし（pの値を保持）
        
        return r

    rgb_image = np.zeros_like(hsl_image)
    q = np.where(L < 0.5, L * (1 + S), L + S - L * S)
    p = 2 * L - q

    rgb_image[:, :, 0] = hue_to_rgb(p, q, H + 1/3)  # R
    rgb_image[:, :, 1] = hue_to_rgb(p, q, H)        # G
    rgb_image[:, :, 2] = hue_to_rgb(p, q, H - 1/3)  # B

    return np.clip(rgb_image * 255, 0, 255).astype(np.uint8)

# 画像データのダミー例を作成してテストする
# image = cv2.imread('/Users/hiyori/kang_hsl/images/Lena.ppm')
# hsl_image = rgb_to_hsl(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
# re_rgb_image = hsl_to_rgb(re_hsl_image)
# cv2.imshow('test', cv2.cvtColor(re_rgb_image, cv2.COLOR_RGB2BGR))
# cv2.waitKey(0)
# cv2.destroyAllWindows()

def hsl_to_cartesian(image):
    """
    HSL色空間で表された画像を直交座標系に変換する。
    :param image: HSL色空間の画像 (高さ x 幅 x 3のNumPy配列)
    :return: 直交座標系に変換された画像 (高さ x 幅 x 3のNumPy配列)
    """
    height, width, _ = image.shape
    cartesian_image = np.zeros_like(image, dtype=float)
    
    h_rad = np.deg2rad(image[:, :, 0])  # 色相をラジアンに変換
    s = image[:, :, 1]  # 彩度
    l = image[:, :, 2]  # 輝度
    
    cartesian_image[:, :, 0] = s * np.cos(h_rad)  # x
    cartesian_image[:, :, 1] = s * np.sin(h_rad)  # y
    cartesian_image[:, :, 2] = l  # z
    
    return cartesian_image

def cartesian_to_hsl(cartesian_image):
    """
    直交座標系に変換されたHSL色空間の画像を通常のHSLに戻す。
    :param cartesian_image: 直交座標系に変換された画像 (高さ x 幅 x 3のNumPy配列)
    :return: HSL色空間の画像 (高さ x 幅 x 3のNumPy配列)
    """
    height, width, _ = cartesian_image.shape
    hsl_image = np.zeros_like(cartesian_image, dtype=float)
    
    # xとy座標から色相(H)と彩度(S)を計算
    x = cartesian_image[:, :, 0]
    y = cartesian_image[:, :, 1]
    hsl_image[:, :, 0] = (np.arctan2(y, x) * (180 / np.pi)) % 360  # 色相H
    hsl_image[:, :, 1] = np.sqrt(x**2 + y**2)  # 彩度S
    
    # z座標は輝度(L)に直接対応
    hsl_image[:, :, 2] = cartesian_image[:, :, 2]  # 輝度L
    
    return hsl_image