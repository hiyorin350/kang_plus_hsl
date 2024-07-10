import numpy as np
import cv2
import conversions

# 画像の読み込み
image = cv2.cvtColor(cv2.imread('images/chart26/chart26.ppm'), cv2.COLOR_BGR2RGB)
height, width, channels = image.shape
N = height * width

hsl_image = conversions.rgb_to_hsl(image)

flat_hsl = hsl_image.reshape(N, 3)

hsl_out = np.zeros_like(image)

u = np.array([ 0.02581419, 0.99966676, 0.]) 
# u = np.array([ 0.0028714, 0.99999588, 0.])
optimum_theta = (np.arctan2(u[1], u[0]))
print(np.rad2deg(optimum_theta))
dichromatic_theta = 50.19789 #ls=50での2色覚平面

# 色相を調整
# print(hsl_image[100, 100, 0])
hsl_image[:, :, 0] = (hsl_image[:, :, 0] + (dichromatic_theta - (np.rad2deg(optimum_theta)) + 90)) % 360
# print(hsl_image[100, 100, 0])
img_out = cv2.cvtColor(conversions.hsl_to_rgb(hsl_image), cv2.COLOR_RGB2BGR)

# 回転された画像を表示
cv2.imwrite('images/chart26/chart26_kang_plus_rotate_hsl.ppm',img_out)
cv2.imshow('cycle_image', img_out)
cv2.waitKey(0)
cv2.destroyAllWindows()