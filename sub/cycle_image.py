import numpy as np
import cv2
import conversions

# 画像の読み込み
image = cv2.imread('/Users/hiyori/kang_plus/images/Chart26.ppm')
height, width, channels = image.shape
N = height * width

hsl_image = conversions.rgb_to_hsl(image)
mhsl_image = conversions.hsl_to_mhsl(hsl_image)

flat_hsl = mhsl_image.reshape(N, 3)

hsl_out = np.zeros_like(image)

u = np.array([ -0.21965709, 0.97557734, 0.])
optimum_theta = (np.arctan2(u[1], u[0]))
dichromatic_theta = 50.19789 #ls=50での2色覚平面

# 色相を調整
mhsl_image[:, :, 0] = (mhsl_image[:, :, 0] + (dichromatic_theta - (np.rad2deg(optimum_theta)) + 90)) % 360

img_out = cv2.cvtColor(conversions.hsl_to_rgb(conversions.mhsl_to_hsl(mhsl_image)), cv2.COLOR_RGB2BGR)

# 回転された画像を表示
cv2.imwrite('/Users/hiyori/kang_plus/images/Chart26_kang_plus_rotate_mhsl.ppm',img_out)
cv2.imshow('cycle_image', img_out)
cv2.waitKey(0)
cv2.destroyAllWindows()