import numpy as np
import cv2
from conversions import *
#FIXME 何も起きていない

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

def project_pixels_to_color_plane(image, u):#TODO HSLでも射影　labから直っていない？
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
            # 画素の色ベクトルを取得 TODO
            color_vector = image[i, j, :]
            # print(image[i, j, :])

            # 色ベクトルを色平面に射影
            projected_vector = color_vector - np.dot(color_vector, u) * u

            # 射影された色ベクトルを保存
            projected_image[i, j, :] = projected_vector

            # print(projected_vector)
            # print()

    return projected_image

# 画像の読み込み
image = cv2.imread('/Users/hiyori/kang_mhsl/images/Chart26_kang_rotate_mhsl.ppm')
print(image[120,110,:])
hsl_image = rgb_to_hsl(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
mhsl_image = hsl_to_mhsl(hsl_image)
mhsl_cartesian = hsl_to_cartesian(mhsl_image)

height, width, _ = image.shape
N = height * width

u = angle_to_normal_vector(50.19789)#2色覚平面(EB平面)

projected_image = project_pixels_to_color_plane(mhsl_cartesian, u)

img_out_rgb = hsl_to_rgb(mhsl_to_hsl(cartesian_to_hsl(projected_image)))

# if img_out_rgb.dtype == np.float32 or img_out_rgb.dtype == np.float64:
#     # 最大値が1.0を超えない場合、255を掛ける　入っていない
#     if img_out_rgb.max() <= 1.0:
#         img_out_rgb = (img_out_rgb * 255).astype(np.uint8)

# 射影された画像を表示
print("done!")

img_out_bgr = cv2.cvtColor(img_out_rgb, cv2.COLOR_RGB2BGR)

cv2.imwrite('/Users/hiyori/kang_mhsl/images/Chart26_kang_pjt_mhsl.ppm',img_out_bgr)
cv2.imshow('hsl_projected', img_out_bgr)
cv2.waitKey(0)
cv2.destroyAllWindows()
