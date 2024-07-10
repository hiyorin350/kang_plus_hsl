import numpy as np
import cv2
import conversions

# 画像の読み込み
image = cv2.cvtColor(cv2.imread('images/chart26/chart26.ppm'), cv2.COLOR_BGR2RGB)
height, width, channels = image.shape
N = height * width

# RGBからHSLへの変換
hsl_image = conversions.rgb_to_hsl(image)

# 画像のフラット化
flat_hsl = hsl_image.reshape(N, 3)

# 出力用の配列
hsl_out = np.zeros_like(image)

u = np.array([0.0728749, 0.99734109, 0.])

# 色相を15度ずつ変えて画像を生成
for i in range(0, 360, 15):
    adjusted_hsl_image = hsl_image.copy()
    adjusted_hsl_image[:, :, 0] = (adjusted_hsl_image[:, :, 0] + i) % 360
    img_out = cv2.cvtColor(conversions.hsl_to_rgb(adjusted_hsl_image), cv2.COLOR_RGB2BGR)
    
    output_filename = f'images/chart26/rotates/chart26_kang_plus_rotate_mhsl_{i}.ppm'
    cv2.imwrite(output_filename, img_out)
    print(f'Saved: {output_filename}')

print("All images have been saved.")
