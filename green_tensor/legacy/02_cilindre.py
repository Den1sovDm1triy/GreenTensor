# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt
import logging
from scipy.sparse.linalg import gmres
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class CylinderDiffraction:
    def __init__(self, k0a, k0L, layer_thicknesses, layer_eps, wave_orientation='along_axis'):
        self.k0a = k0a
        self.k0L = k0L
        self.layer_thicknesses = np.array(layer_thicknesses)
        self.layer_eps = np.array(layer_eps, dtype=complex)
        self.wave_orientation = wave_orientation
        
        if not np.isclose(np.sum(self.layer_thicknesses), 1.0, atol=1e-3):
            raise ValueError("Сумма относительных толщин слоев должна быть равна 1")
        
        if len(self.layer_thicknesses) != len(self.layer_eps):
            raise ValueError("Количество слоев в thicknesses и eps должно совпадать")
        
        # Волновое число в безразмерных единицах
        self.k0 = 1.0
        
        # Параметры дискретизации (уменьшено для производительности)
        self.N_rho = 10
        self.N_phi = 10
        self.N_z = 10
        
        # Сетки координат
        self.rho_grid = None
        self.phi_grid = None
        self.z_grid = None
        
        # Поля и результаты
        self.internal_field = None
        self.far_field_theta = None
        self.far_field_phi = None
        
        # Убрана поправка Лоренца как нефизичная в данном контексте
        self.diagonal_correction = 0.0
    
    def _build_radial_grid(self):
        layer_thicknesses_abs = self.layer_thicknesses * self.k0a
        boundaries = [0.0]
        for thickness in layer_thicknesses_abs:
            boundaries.append(boundaries[-1] + thickness)
        
        rho_segments = []
        for i, thickness in enumerate(layer_thicknesses_abs):
            n_in_layer = max(3, int(self.N_rho * self.layer_thicknesses[i]))
            start = boundaries[i] + 1e-12
            end = boundaries[i+1]
            segment = np.linspace(start, end, n_in_layer, endpoint=False)
            rho_segments.append(segment)
        
        self.rho_grid = np.concatenate(rho_segments)
        self.N_rho = len(self.rho_grid)
        self.phi_grid = np.linspace(0, 2*np.pi, self.N_phi, endpoint=False)
        self.z_grid = np.linspace(-self.k0L/2, self.k0L/2, self.N_z)
        logging.info(f"Сетка построена: rho={len(self.rho_grid)}, phi={len(self.phi_grid)}, z={len(self.z_grid)}")
    
    def _epsilon_r(self, rho_index):
        rho = self.rho_grid[rho_index]
        cumulative = 0.0
        for i, rel_thickness in enumerate(self.layer_thicknesses):
            abs_thickness = rel_thickness * self.k0a
            cumulative += abs_thickness
            if rho <= cumulative:
                return self.layer_eps[i]
        return self.layer_eps[-1]
    
    def _g_scalar(self, r, r_prime):
        dx = r[0] - r_prime[0]
        dy = r[1] - r_prime[1]
        dz = r[2] - r_prime[2]
        R = np.sqrt(dx**2 + dy**2 + dz**2 + 1e-16)
        return np.exp(1j * self.k0 * R) / (4 * np.pi * R)
    
    def _G_tensor(self, r, r_prime):
        dx = r[0] - r_prime[0]
        dy = r[1] - r_prime[1]
        dz = r[2] - r_prime[2]
        R = np.sqrt(dx**2 + dy**2 + dz**2 + 1e-16)
        
        if R > 1e-8:
            R_hat = np.array([dx/R, dy/R, dz/R])
        else:
            R_hat = np.zeros(3)
        
        g = np.exp(1j * self.k0 * R) / (4 * np.pi * R)
        G = np.zeros((3, 3), dtype=complex)
        
        for i in range(3):
            for j in range(3):
                delta_ij = 1 if i == j else 0
                term1 = (delta_ij - R_hat[i] * R_hat[j]) * g
                term2 = (delta_ij - 3 * R_hat[i] * R_hat[j]) * (1 - 1j*self.k0*R) * g / (self.k0**2 * R**2)
                G[i, j] = term1 + term2
        
        return G
    
    def _cylindrical_to_cartesian(self, rho, phi, z):
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return np.array([x, y, z])
    
    def _incident_field(self, rho, phi, z, component):
        # Преобразуем цилиндрические координаты в декартовы
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        
        if self.wave_orientation == 'along_axis':
            # Волна распространяется вдоль оси Z
            if component == 0:  # x-компонента
                return 0.0
            elif component == 1:  # y-компонента
                return 0.0
            elif component == 2:  # z-компонента
                return np.exp(-1j * self.k0 * z)
        else:
            # Волна распространяется вдоль оси X (перпендикулярно оси цилиндра)
            if component == 0:  # x-компонента
                return 0.0
            elif component == 1:  # y-компонента
                return 0.0
            elif component == 2:  # z-компонента
                return np.exp(-1j * self.k0 * x)
    
    def solve_internal_field(self):
        self._build_radial_grid()
        
        drho = self.rho_grid[1] - self.rho_grid[0] if len(self.rho_grid) > 1 else self.k0a/self.N_rho
        dphi = self.phi_grid[1] - self.phi_grid[0] if len(self.phi_grid) > 1 else 2*np.pi/self.N_phi
        dz = self.z_grid[1] - self.z_grid[0] if len(self.z_grid) > 1 else self.k0L/self.N_z
        
        N_points = self.N_rho * self.N_phi * self.N_z  # Теперь включает все азимуты
        N_total = 3 * N_points
        logging.info(f"Решение для {N_total} неизвестных...")
        
        A = np.zeros((N_total, N_total), dtype=complex)
        b = np.zeros(N_total, dtype=complex)
        
        # Создаем список всех точек с их координатами и индексами
        points = []
        for i_rho in range(self.N_rho):
            for i_phi in range(self.N_phi):  # Добавлен цикл по phi
                for i_z in range(self.N_z):
                    rho = self.rho_grid[i_rho]
                    phi = self.phi_grid[i_phi]
                    z = self.z_grid[i_z]
                    cart_coord = self._cylindrical_to_cartesian(rho, phi, z)
                    points.append({
                        'rho_idx': i_rho,
                        'phi_idx': i_phi,
                        'z_idx': i_z,
                        'cart_coord': cart_coord
                    })
        
        # Рассчитываем объемный элемент (корректно для цилиндрических координат)
        for i, p_i in enumerate(points):
            idx_x, idx_y, idx_z = 3*i, 3*i+1, 3*i+2
            
            # Падающее поле в декартовых компонентах
            b[idx_x] = self._incident_field(p_i['cart_coord'][0], p_i['cart_coord'][1], p_i['cart_coord'][2], 0)
            b[idx_y] = self._incident_field(p_i['cart_coord'][0], p_i['cart_coord'][1], p_i['cart_coord'][2], 1)
            b[idx_z] = self._incident_field(p_i['cart_coord'][0], p_i['cart_coord'][1], p_i['cart_coord'][2], 2)
            
            for j, p_j in enumerate(points):
                jdx_x, jdx_y, jdx_z = 3*j, 3*j+1, 3*j+2
                G_tensor = self._G_tensor(p_i['cart_coord'], p_j['cart_coord'])
                delta_eps = self._epsilon_r(p_j['rho_idx']) - 1.0
                
                # Корректный объемный элемент для цилиндрических координат
                rho_j = self.rho_grid[p_j['rho_idx']]
                dV = rho_j * drho * dphi * dz
                
                coeff = -self.k0**2 * delta_eps * dV
                
                if i == j:
                    # Диагональные элементы (без поправки Лоренца)
                    A[idx_x, jdx_x] = 1.0
                    A[idx_y, jdx_y] = 1.0
                    A[idx_z, jdx_z] = 1.0
                else:
                    # Внедиагональные элементы
                    A[idx_x, jdx_x] += coeff * G_tensor[0, 0]
                    A[idx_x, jdx_y] += coeff * G_tensor[0, 1]
                    A[idx_x, jdx_z] += coeff * G_tensor[0, 2]
                    
                    A[idx_y, jdx_x] += coeff * G_tensor[1, 0]
                    A[idx_y, jdx_y] += coeff * G_tensor[1, 1]
                    A[idx_y, jdx_z] += coeff * G_tensor[1, 2]
                    
                    A[idx_z, jdx_x] += coeff * G_tensor[2, 0]
                    A[idx_z, jdx_y] += coeff * G_tensor[2, 1]
                    A[idx_z, jdx_z] += coeff * G_tensor[2, 2]
        
        # Обработка численных ошибок
        A = np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0)
        b = np.nan_to_num(b, nan=0.0, posinf=0.0, neginf=0.0)
        
        cond_number = np.linalg.cond(A)
        logging.info(f"Число обусловленности матрицы: {cond_number:.2e}")
        
        try:
            E_solution = solve(A, b)
            logging.info("Система успешно решена прямым методом")
        except np.linalg.LinAlgError as e:
            logging.error(f"Ошибка решения СЛАУ: {e}")
            E_solution, info = gmres(A, b, atol=1e-5)
            if info == 0:
                logging.info("GMRES успешно сошелся")
            else:
                logging.warning("GMRES не сошелся, используем псевдообратную матрицу")
                E_solution = np.linalg.pinv(A) @ b
        
        # Изменение формы массива внутреннего поля (rho, phi, z, компоненты)
        self.internal_field = E_solution.reshape((self.N_rho, self.N_phi, self.N_z, 3))
        
        max_field = np.max(np.abs(self.internal_field))
        min_field = np.min(np.abs(self.internal_field))
        logging.info(f"Внутреннее поле: max={max_field:.4e}, min={min_field:.4e}")
        
        return self.internal_field
    
    def compute_far_field(self, theta, phi_obs=0):
        if self.internal_field is None:
            raise RuntimeError("Сначала выполните solve_internal_field")
        
        drho = self.rho_grid[1] - self.rho_grid[0] if len(self.rho_grid) > 1 else self.k0a/self.N_rho
        dphi = self.phi_grid[1] - self.phi_grid[0] if len(self.phi_grid) > 1 else 2*np.pi/self.N_phi
        dz = self.z_grid[1] - self.z_grid[0] if len(self.z_grid) > 1 else self.k0L/self.N_z
        
        n_theta = len(theta)
        E_theta = np.zeros(n_theta, dtype=complex)
        E_phi = np.zeros(n_theta, dtype=complex)
        const = 1j * self.k0**2 / (4 * np.pi)
        
        # Направление наблюдения в декартовых координатах
        r_hat_x = np.sin(theta) * np.cos(phi_obs)
        r_hat_y = np.sin(theta) * np.sin(phi_obs)
        r_hat_z = np.cos(theta)
        
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        cos_phi_obs = np.cos(phi_obs)
        sin_phi_obs = np.sin(phi_obs)
        
        # Матрицы преобразования для компонент theta и phi
        T_theta = np.zeros((n_theta, 3))
        T_theta[:, 0] = cos_theta * cos_phi_obs
        T_theta[:, 1] = cos_theta * sin_phi_obs
        T_theta[:, 2] = -sin_theta
        
        T_phi = np.zeros((n_theta, 3))
        T_phi[:, 0] = -sin_phi_obs
        T_phi[:, 1] = cos_phi_obs
        T_phi[:, 2] = 0.0  # Важно: z-компонента = 0
        
        # Тензор Грина для дальней зоны
        I = np.eye(3)
        r_outer = np.zeros((3, 3))
        r_outer[0, 0] = r_hat_x[0] * r_hat_x[0] if n_theta == 1 else 1  # Заглушка для скалярного случая
        
        # Суммирование по всем точкам источника
        for i_rho in range(self.N_rho):
            rho = self.rho_grid[i_rho]
            delta_eps = self._epsilon_r(i_rho) - 1.0
            
            for i_phi in range(self.N_phi):
                phi_src = self.phi_grid[i_phi]
                for i_z in range(self.N_z):
                    z_val = self.z_grid[i_z]
                    
                    # Декартовы координаты источника
                    x_src, y_src, z_src = self._cylindrical_to_cartesian(rho, phi_src, z_val)
                    
                    # Фаза для направления наблюдения
                    k_dot_r = self.k0 * (r_hat_x * x_src + r_hat_y * y_src + r_hat_z * z_src)
                    phase = np.exp(-1j * k_dot_r)
                    
                    # Поле в точке источника
                    E_src = self.internal_field[i_rho, i_phi, i_z]
                    
                    # Вычисление тензора Грина для дальней зоны
                    G_ee = np.zeros((n_theta, 3, 3), dtype=complex)
                    for i in range(n_theta):
                        # Тензорная структура: I - r_hat ⊗ r_hat
                        r_hat = np.array([r_hat_x[i], r_hat_y[i], r_hat_z[i]])
                        G_ee[i] = (I - np.outer(r_hat, r_hat)) * phase[i]
                    
                    # Рассеянное поле от элемента объема
                    E_sc = const * delta_eps * np.einsum('kij,j->ki', G_ee, E_src)
                    
                    # Проекция на компоненты theta и phi
                    E_theta_i = np.sum(T_theta * E_sc, axis=1)
                    E_phi_i = np.sum(T_phi * E_sc, axis=1)
                    
                    # Объемный элемент и суммирование
                    dV = rho * drho * dphi * dz
                    E_theta += E_theta_i * dV
                    E_phi += E_phi_i * dV
        
        return E_theta, E_phi
    
    def compute_rcs(self, theta, phi_obs=0):
        """
        Вычисление ЭПР (RCS) с нормировкой на длину цилиндра
        
        :param theta: Угол от оси z в радианах (массив)
        :param phi_obs: Азимутальный угол наблюдения в радианах (по умолчанию 0)
        :return: Нормированный ЭПР (σ/λ)
        """
        E_theta, E_phi = self.compute_far_field(theta, phi_obs)
        rcs = 4 * np.pi * (np.abs(E_theta)**2 + np.abs(E_phi)**2)
        return rcs / self.k0L  # Нормировка на длину цилиндра
    
    def plot_radiation_pattern(self, theta_steps=181):
        """
        Построение диаграммы направленности в полярных координатах
        для двух плоскостей (XZ и YZ) при перпендикулярном падении
        с полным кругом 360 градусов
        """
        # Углы theta от 0 до 360 градусов (полный круг)
        theta_deg = np.linspace(0, 360, theta_steps)
        theta_rad = np.deg2rad(theta_deg)
        
        # Вычисление поля в дальней зоне
        E_theta, E_phi = self.compute_far_field(theta_rad)
        
        # Мощность в дБ
        P_theta = 20 * np.log10(np.abs(E_theta) + 1e-12)
        P_phi = 20 * np.log10(np.abs(E_phi) + 1e-12)
        
        # Нормировка
        P_theta -= np.max(P_theta)
        P_phi -= np.max(P_phi)
        
        plt.figure(figsize=(10, 8))
        ax = plt.subplot(111, projection='polar')
        
        # Построение на полярной диаграмме (полный круг)
        ax.plot(theta_rad, P_theta, 'b-', linewidth=2, label='Eθ')
        ax.plot(theta_rad, P_phi, 'r-', linewidth=2, label='Eφ')
        ax.set_theta_zero_location('E')  # 0 градусов справа
        ax.set_theta_direction(-1)       # По часовой стрелке
        ax.set_rlabel_position(90)       # Положение меток радиальной оси
        ax.set_ylim(-40, 0)
        ax.legend(loc='lower right')
        ax.set_title(f'Диаграмма направленности (360°)\nПадение: {self.wave_orientation}, k0a={self.k0a:.2f}, k0L={self.k0L:.2f}')
        
        plt.tight_layout()
        plt.show()
    
    def plot_xy_plane_pattern(self):
        """
        Построение диаграммы направленности в плоскости XY (горизонтальной)
        для случая поперечного падения волны
        """
        if self.wave_orientation != 'perpendicular':
            logging.warning("Плоскость XY имеет смысл только для поперечного падения волны")
        
        # Фиксируем theta = 90° (горизонтальная плоскость)
        theta = np.array([np.pi/2])
        
        # Углы phi от 0 до 360 градусов
        phi_deg = np.linspace(0, 360, 361)
        phi_rad = np.deg2rad(phi_deg)
        
        # Вычисление RCS с нормировкой на длину
        rcs_normalized = np.zeros_like(phi_deg)
        for i, phi in enumerate(phi_rad):
            rcs_normalized[i] = self.compute_rcs(theta, phi)
        
        # Переводим в дБ
        rcs_db = 10 * np.log10(rcs_normalized + 1e-12)
        
        plt.figure(figsize=(10, 8))
        ax = plt.subplot(111, projection='polar')
        
        # Построение на полярной диаграмме
        ax.plot(phi_rad, rcs_db, 'g-', linewidth=2)
        ax.set_theta_zero_location('E')  # 0 градусов справа
        ax.set_theta_direction(-1)       # По часовой стрелке
        ax.set_rlabel_position(90)       # Положение меток радиальной оси
        ax.set_title(f'Диаграмма направленности в плоскости XY (θ=90°)\nПадение: перпендикулярно, k0a={self.k0a:.2f}, k0L={self.k0L:.2f}')
        
        plt.tight_layout()
        plt.show()
    
    def plot_internal_field(self, component='magnitude', phi_index=0):
        """
        Визуализация внутреннего поля в безразмерных координатах
        :param component: 'x', 'y', 'z', 'magnitude' или 'phase'
        :param phi_index: Индекс азимута для визуализации
        """
        if self.internal_field is None:
            raise RuntimeError("Сначала выполните solve_internal_field")
        
        if phi_index >= self.N_phi:
            raise ValueError("Недопустимый индекс азимута")
        
        plt.figure(figsize=(10, 6))
        
        # Выбор компоненты для визуализации
        if component == 'x':
            field_data = np.real(self.internal_field[:, phi_index, :, 0])
            label = 'Re(Ex)'
        elif component == 'y':
            field_data = np.real(self.internal_field[:, phi_index, :, 1])
            label = 'Re(Ey)'
        elif component == 'z':
            field_data = np.real(self.internal_field[:, phi_index, :, 2])
            label = 'Re(Ez)'
        elif component == 'magnitude':
            field_slice = self.internal_field[:, phi_index, :, :]
            field_data = np.sqrt(np.abs(field_slice[:, :, 0])**2 + 
                               np.abs(field_slice[:, :, 1])**2 + 
                               np.abs(field_slice[:, :, 2])**2)
            label = '|E|'
        elif component == 'phase':
            field_data = np.angle(self.internal_field[:, phi_index, :, 0])
            label = 'Phase(Ex)'
        else:
            raise ValueError("Недопустимое значение component. Допустимы: 'x', 'y', 'z', 'magnitude', 'phase'")
        
        # Определение типа поля для заголовка
        if self.wave_orientation == 'along_axis':
            wave_type = "вдоль оси"
        else:
            wave_type = "перпендикулярно оси"
        
        plt.imshow(field_data, 
                   extent=[self.z_grid[0], self.z_grid[-1], 
                           self.rho_grid[0], self.rho_grid[-1]],
                   aspect='auto', cmap='RdBu', origin='lower')
        plt.colorbar(label=label)
        plt.xlabel('k₀·z')
        plt.ylabel('k₀·ρ')
        plt.title(f'Внутреннее поле ({label}) при φ={np.rad2deg(self.phi_grid[phi_index]):.1f}°\nПадение волны {wave_type}, k₀a={self.k0a:.2f}, k₀L={self.k0L:.2f}')
        plt.show()

class CylinderVisualization:
    def __init__(self, cylinder):
        self.cylinder = cylinder
        if self.cylinder.internal_field is None:
            raise RuntimeError("Сначала выполните solve_internal_field для цилиндра")
    
    def plot_rcs_vs_radius(self, radii, fixed_k0L=3.0):
        """
        Радиолокационный коэффициент рассеяния (RCS) в зависимости от радиуса
        с нормировкой на длину цилиндра
        
        :param radii: Массив безразмерных радиусов (k0a)
        :param fixed_k0L: Фиксированная безразмерная длина цилиндра
        """
        from tqdm import tqdm
        rcs_values = []
        
        logging.info("Расчет RCS в зависимости от радиуса...")
        for k0a in tqdm(radii):
            # Создаем временный цилиндр с текущим радиусом
            temp_cylinder = CylinderDiffraction(
                k0a, fixed_k0L,
                self.cylinder.layer_thicknesses,
                self.cylinder.layer_eps,
                self.cylinder.wave_orientation
            )
            temp_cylinder.solve_internal_field()
            
            # Вычисляем нормированный RCS в обратном направлении (θ=180°)
            theta = np.array([np.pi])  # 180 градусов
            rcs_norm = temp_cylinder.compute_rcs(theta)
            rcs_values.append(rcs_norm[0])
        
        # Визуализация
        plt.figure(figsize=(10, 6))
        plt.plot(radii, 10 * np.log10(rcs_values), 'b-o', linewidth=2)
        plt.xlabel('Безразмерный радиус (k₀a)')
        plt.ylabel('Нормированный RCS (дБ)')
        plt.title(f'Нормированный RCS (σ/λ) в зависимости от радиуса\nФиксированная длина k₀L={fixed_k0L:.1f}')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        
        return radii, rcs_values
    
    def plot_rcs_vs_length(self, lengths, fixed_k0a=0.5):
        """
        Радиолокационный коэффициент рассеяния (RCS) в зависимости от длины
        с нормировкой на длину цилиндра
        
        :param lengths: Массив безразмерных длин (k0L)
        :param fixed_k0a: Фиксированный безразмерный радиус цилиндра
        """
        from tqdm import tqdm
        rcs_values = []
        
        logging.info("Расчет RCS в зависимости от длины...")
        for k0L in tqdm(lengths):
            # Создаем временный цилиндр с текущей длиной
            temp_cylinder = CylinderDiffraction(
                fixed_k0a, k0L,
                self.cylinder.layer_thicknesses,
                self.cylinder.layer_eps,
                self.cylinder.wave_orientation
            )
            temp_cylinder.solve_internal_field()
            
            # Вычисляем нормированный RCS в обратном направлении (θ=180°)
            theta = np.array([np.pi])  # 180 градусов
            rcs_norm = temp_cylinder.compute_rcs(theta)
            rcs_values.append(rcs_norm[0])
        
        # Визуализация
        plt.figure(figsize=(10, 6))
        plt.plot(lengths, 10 * np.log10(rcs_values), 'r-o', linewidth=2)
        plt.xlabel('Безразмерная длина (k₀L)')
        plt.ylabel('Нормированный RCS (дБ)')
        plt.title(f'Нормированный RCS (σ/λ) в зависимости от длины\nФиксированный радиус k₀a={fixed_k0a:.1f}')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        
        return lengths, rcs_values

# Пример использования
if __name__ == "__main__":
    # Безразмерные параметры цилиндра (уменьшены для скорости)
    k0a = 0.1*np.pi  # Безразмерный радиус
    k0L = 100.0*np.pi  # Безразмерная длина
    
    # Многослойная структура (ядро + оболочка)
    layer_thicknesses = [1]  # 70% ядро, 30% оболочка
    layer_eps = [1]  # Комплексные ε
    
    # Случай 1: Волна вдоль оси цилиндра
    print("Расчет для волны вдоль оси...")
    cylinder_along = CylinderDiffraction(k0a, k0L, layer_thicknesses, layer_eps, 'along_axis')
    cylinder_along.solve_internal_field()
    cylinder_along.plot_radiation_pattern()  # ДН 360°
    cylinder_along.plot_internal_field(component='magnitude', phi_index=0)
    
    # Случай 2: Волна перпендикулярна оси цилиндра
    print("\nРасчет для волны перпендикулярно оси...")
    cylinder_perp = CylinderDiffraction(k0a, k0L, layer_thicknesses, layer_eps, 'perpendicular')
    cylinder_perp.solve_internal_field()
    cylinder_perp.plot_radiation_pattern()  # ДН 360°
    cylinder_perp.plot_xy_plane_pattern()   # ДН в плоскости XY
    cylinder_perp.plot_internal_field(component='magnitude', phi_index=0)
    
    # Визуализация RCS с нормировкой
    vis = CylinderVisualization(cylinder_perp)
    radii = np.linspace(0.1, 30.0, 90)  # Уменьшено количество точек
    vis.plot_rcs_vs_radius(radii, fixed_k0L=3.0)
    
    lengths = np.linspace(1.0, 5.0, 5)  # Уменьшено количество точек
    vis.plot_rcs_vs_length(lengths, fixed_k0a=1.0)