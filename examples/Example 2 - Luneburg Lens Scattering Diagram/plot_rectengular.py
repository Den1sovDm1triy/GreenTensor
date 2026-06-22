#DATA FROM HFSS // need upload Ephi=0.csv
#########################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Загружаем данные из файла — положите CSV рядом со скриптом
file_path = "Ephi=0_from_HFSS.csv"
print(os.path.exists(file_path))  # Проверяем, существует ли файл

df = pd.read_csv(file_path)

# Извлекаем необходимые данные
theta_deg = df["Theta [deg]"]
r_db = df["dB20normalize(rETotal) []"]

# Преобразуем углы в радианы
theta_rad = np.radians(theta_deg)

# Строим график
plt.figure(figsize=(8, 6))
plt.plot(theta_rad, r_db, color='blue', linestyle='-', linewidth=1, label="Experimental Data")

# Настройки графика
plt.xlabel("Theta (radians)")
plt.ylabel("dB20normalize(rETotal)")
plt.title("Experimental Data Plot")
plt.legend()
plt.grid(True)

# Отображаем график
plt.show()