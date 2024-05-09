import numpy as np
import cv2

def clip_lab_within_rgb_gamut(lab_color, step=1):
    """
    Lab色空間の色をRGB色空間に収まるように調整する。この際、L（明度）とh（色相）を保持し、C（彩度）のみを調整する。

    引数:
    lab_color (numpy.ndarray): Lab色を表す配列。
    step (int): Cを減少させる際のステップサイズ。

    戻り値:
    numpy.ndarray: 調整されたLab色。
    """
    L, a, b = lab_color

    # 初期のC（彩度）とh（色相）を計算
    C = np.sqrt(a**2 + b**2)
    h = np.degrees(np.arctan2(b, a))

    # LabからRGBへ変換
    lab_color_reshaped = lab_color.reshape(1, 1, 3)
    rgb_color = cv2.cvtColor(lab_color_reshaped, cv2.COLOR_Lab2BGR)

    # Cを減少させながらRGB色空間内に収まるかチェック
    while not np.all((rgb_color >= 0) & (rgb_color <= 255)):
        C -= step
        if C < 0:
            C = 0
            break

        a = C * np.cos(np.radians(h))
        b = C * np.sin(np.radians(h))
        lab_color = np.array([L, a, b], dtype=np.float32)
        lab_color_reshaped = lab_color.reshape(1, 1, 3)
        rgb_color = cv2.cvtColor(lab_color_reshaped, cv2.COLOR_Lab2BGR)

    return np.array([L, a, b])

def rotate_hue_lab(lab_color, angle):
    """
    LAB色空間において、指定した角度だけ色相（h）を回転させる。

    引数:
    lab_color (numpy.ndarray): Lab色を表す配列。
    angle (float): 回転させる角度（度単位）。

    戻り値:
    numpy.ndarray: 色相が回転されたLab色。
    """
    n, _ = lab_color.shape
    l = np.zeros(n)
    a = np.zeros(n)
    b = np.zeros(n)

    for i in range(n):
        l[i], a[i], b[i] = lab_color[i][0], lab_color[i][1], lab_color[i][2]
    
    h = np.degrees(np.arctan2(b, a))  # 現在の色相を計算
    h_rotated = (h + angle) % 360  # 角度を加え、360度で割った余りを取る

    # 回転した色相で新しいaとbを計算
    C = np.sqrt(a**2 + b**2)
    a_rotated = C * np.cos(np.radians(h_rotated))
    b_rotated = C * np.sin(np.radians(h_rotated))

    return np.array([l, a_rotated, b_rotated])

def cycle(image, u):
    """
    LAB色空間において、指定した角度だけ色相（h）を回転させる。

    引数:
    image (numpy.ndarray): BGR色空間を表す配列。
    u (NDArray[Any]): 最適色空間の法線ベクトル。

    戻り値:
    img_out: BGR色空間での回転後の画像を表す配列
    """
    # 画像の読み込み
    # image = cv2.imread('/Users/hiyori/kang/images/Chart26.ppm')

    lab_image = cv2.cvtColor(image, cv2.COLOR_BGR2LAB)

    height, width, _ = image.shape
    N = height * width

    flat_lab_image = lab_image.reshape(N, 3)

    lab_out = np.zeros_like(image)

    lab_process = np.zeros((N,3))#8bit扱いから32bit扱いに変換
    for i in range(N):
        lab_process[i][0] = 100 * flat_lab_image[i][0] / 255
        lab_process[i][1] = flat_lab_image[i][1] - 128
        lab_process[i][2] = flat_lab_image[i][2] - 128

    #lab_processは(縦、横、チャンネル)形式だが、rotateは(l,a,b)形式を引数としている

    flat_lab = lab_process.reshape(N, 3)

    # u = np.array([ 0. ,-0.07, 0.99])
    optimum_theta = (np.arctan2(u[2], u[1]))
    dichromatic_theta = 90 + 11.48

    rotated_lab = rotate_hue_lab(flat_lab, dichromatic_theta - (np.rad2deg(optimum_theta)) + 90)

    l_out = np.zeros(N)
    a_out = np.zeros(N)
    b_out = np.zeros(N)

    for i in range(N):
        l_out[i] = 255 * rotated_lab[0,i] / 100
        a_out[i] = rotated_lab[1,i] + 128
        b_out[i] = rotated_lab[2,i] + 128

    count = 0
    for i in range(height):
        for j in range(width):
            lab_out[i][j] = (l_out[count], a_out[count], b_out[count])
            count += 1

    img_out = cv2.cvtColor(lab_out, cv2.COLOR_LAB2BGR)

    return img_out