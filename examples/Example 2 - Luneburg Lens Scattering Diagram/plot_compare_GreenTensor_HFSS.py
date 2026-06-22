#DATA FROM HFSS // need upload Ephi=90_from_HFSS.csv ||| compare with GreenTensor
######################################################################

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Загружаем данные из файла (HFSS) — положите CSV рядом со скриптом
file_path = "Ephi=90_from_HFSS.csv"
print(os.path.exists(file_path))  # Проверяем, существует ли файл
df = pd.read_csv(file_path)

# Извлекаем необходимые данные
theta_deg = df["Theta [deg]"]
r_db = df["dB20normalize(rETotal) []"]

# Преобразуем углы в радианы
theta_rad = np.radians(theta_deg)

# Создаем полярный график
fig, axs = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})

# График из GreenTensor (доступно только при запуске в одном окружении с расчётным скриптом)
if "teta" in dir() and "DN_NORM_phi" in dir():
    axs.plot(teta, DN_NORM_phi, color='orange', linestyle='-', linewidth=1, label='GreenTensor')

# График из HFSS
axs.plot(theta_rad, r_db, color='blue', linestyle='--', linewidth=1, label="HFSS")

# Настройки графика
axs.set_theta_zero_location("N")  # 0° сверху
axs.set_theta_direction(-1)  # Углы идут против часовой стрелки
axs.legend()

# Отображаем график
plt.show()