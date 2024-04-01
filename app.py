import numpy as np
import cv2
from functions import *
from scipy.optimize import minimize

image = cv2.imread('/Users/hiyori/kang/images/chart26.ppm')
assert image is not None, "読み込みに失敗しました"

height, width, _ = image.shape
N = height * width

result_image = image

# G(u)を計算する関数、今回は未搭載
# g = global_contrast.global_contrast(cv2.cvtColor(image, cv2.COLOR_BGR2LAB))

#TODO J(u)を計算する関数、今回は未搭載
# j = alpha * l + (1 - alpha) * g

#Xlを計算する関数
Xl = calculate_color_difference_vectors_with_gaussian_pairing.calculate_color_difference_vectors_with_gaussian_pairing(image)
Al = Xl.T @ Xl

# 目的関数と制約条件の定義
def objective(u):
    return ((u.T @ Al @ u) / N)

# L軸に垂直な制約
def constraint_perpendicular_to_L(u):
    e = np.array([0, 0, 1])  # L軸を指す標準基底ベクトル
    return np.dot(u, e)

# 単位ベクトル制約
def constraint_unit_vector(u):
    return np.dot(u, u) - 1

# 初期値の設定
u0 = np.random.rand(3)
u0 /= np.linalg.norm(u0) #単位ベクトル化

# 制約の設定
constraints = [{'type': 'eq', 'fun': constraint_perpendicular_to_L},
               {'type': 'eq', 'fun': constraint_unit_vector}]

# 逐次二次計画法による最適化問題の解決
res = minimize(objective, u0, constraints=constraints, method = "SLSQP")

# 最適化された u の値
optimized_u = res.x
optimized_u = np.reshape(optimized_u, (3,1))

# 最適化された結果の u^T Al u の値を計算
optimized_value = (optimized_u.T @ Al @ optimized_u) / N

# 最適色平面と、二色覚平面の差分だけ回す
img_out = cycle.cycle(image, optimized_u)#TODO cycle_imageに変更

#最終的な画像を出す
print("done!")
print(optimized_u)
# cv2.imshow('result', img_out)
# cv2.waitKey()
# cv2.destroyAllWindows()