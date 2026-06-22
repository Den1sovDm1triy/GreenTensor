import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import LinearOperator, gmres 
import time
from tqdm import tqdm
from scipy.special import jv, hankel1
from numba import jit, prange
import concurrent.futures
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

class MultilayerCuboidVIESolver: 
    def __init__(self, freq_hz, kappa_a, kappa_b, kappa_L, layer_fracs, layer_epsr, 
                 grid_shape=(24, 24, 28), padding_factor=2, eta_reg_frac=1e-3, verbose=True):
        
        # Проверка входных параметров
        assert kappa_a > 0 and kappa_b > 0 and kappa_L > 0, "Размеры должны быть положительными"
        assert all(np.real(epsr) > 0 for epsr in layer_epsr), "Действительная часть εr должна быть положительной"
        assert all(np.imag(epsr) <= 0 for epsr in layer_epsr), "Мнимая часть εr должна быть неположительной"
        assert padding_factor >= 1, "Коэффициент дополнения должен быть >=1"
        assert len(layer_fracs) == len(layer_epsr), "Количество слоев не совпадает"
        
        self.c0 = 299792458.0
        self.eps0 = 8.854187817e-12
        self.mu0 = 4e-7 * np.pi
        self.eta0 = np.sqrt(self.mu0 / self.eps0)

        # Частота, волновое число
        self.freq = freq_hz
        self.omega = 2 * np.pi * self.freq
        self.k0 = self.omega * np.sqrt(self.eps0 * self.mu0)
        self.lam = 2 * np.pi / self.k0

        # Геометрия
        self.a = kappa_a / self.k0
        self.b = kappa_b / self.k0
        self.L = kappa_L / self.k0

        # Слои
        S = sum(layer_fracs)
        self.layer_fracs = [f / S for f in layer_fracs]
        self.layer_epsr = [complex(epsr) for epsr in layer_epsr]
        self.N_layers = len(self.layer_fracs)

        # Сетка
        self.Nx, self.Ny, self.Nz = grid_shape
        self.x_min, self.x_max = -self.a / 2, self.a / 2
        self.y_min, self.y_max = -self.b / 2, self.b / 2
        self.z_min, self.z_max = -self.L / 2, self.L / 2

        self.x = np.linspace(self.x_min, self.x_max, self.Nx, endpoint=False)
        self.y = np.linspace(self.y_min, self.y_max, self.Ny, endpoint=False)
        self.z = np.linspace(self.z_min, self.z_max, self.Nz, endpoint=False)
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
        self.dV = self.dx * self.dy * self.dz

        self.X, self.Y, self.Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # Характеристическая функция
        self.chi = self._build_chi()

        # Спектральный оператор
        self.padding_factor = int(padding_factor)
        self.Npx = self.padding_factor * self.Nx
        self.Npy = self.padding_factor * self.Ny
        self.Npz = self.padding_factor * self.Nz
        self.eta_reg = eta_reg_frac * (self.k0 ** 2)
        self._build_spectral_green()

        # Плейсхолдеры
        self._Einc = None
        self._Esol = None
        self._Esca = None
        self.verbose = verbose
        self._s_matrix_cache = {}

    def _build_chi(self):
        chi = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.complex128)
        
        # Границы слоёв по z (должны быть в пределах z_min..z_max)
        z_bounds = [self.z_min]
        acc = self.z_min
        for f in self.layer_fracs:
            acc += f * self.L
            z_bounds.append(min(acc, self.z_max))  # Гарантируем не выход за границы
        
        # Проверка границ
        print(f"Z bounds: {z_bounds}")
        print(f"Z range: {self.z_min} to {self.z_max}")
        
        # Создание 3D-маски с правильными границами
        mask_xy = (np.abs(self.X) <= self.a/2) & (np.abs(self.Y) <= self.b/2)
        
        for m in range(self.N_layers):
            zl, zr = z_bounds[m], z_bounds[m+1]
            print(f"Layer {m}: z={zl:.3f} to {zr:.3f}, epsr={self.layer_epsr[m]}")
            
            # Строгое условие для слоев (последний слой включает правую границу)
            if m == self.N_layers - 1:
                mask_z = (self.Z >= zl) & (self.Z <= zr)
            else:
                mask_z = (self.Z >= zl) & (self.Z < zr)
            
            mask_layer = mask_xy & mask_z
            chi[mask_layer] = (self.layer_epsr[m] - 1.0)
        
        # Проверка результата
        print(f"Chi non-zero: {np.count_nonzero(chi)}/{chi.size}")
        print(f"Chi values: {np.unique(chi)}")
        
        return chi

    @property
    def Einc(self):
        return self._Einc
    
    @property
    def Esol(self):
        return self._Esol
    
    @property
    def Esca(self):
        if self._Esca is None and self._Esol is not None and self._Einc is not None:
            self._Esca = self._Esol - self._Einc
        return self._Esca

    def _s_matrix_layer(self, k_rho, k0, eps1, eps2):
        """Вычисляет S-матрицу для перехода между двумя слоями."""
        kz1 = np.lib.scimath.sqrt(eps1 * k0**2 - k_rho**2)
        kz2 = np.lib.scimath.sqrt(eps2 * k0**2 - k_rho**2)
        
        # Коэффициенты отражения и прохождения
        r_te = (kz1 - kz2) / (kz1 + kz2)
        t_te = 2 * kz1 / (kz1 + kz2)
        r_tm = (eps2 * kz1 - eps1 * kz2) / (eps2 * kz1 + eps1 * kz2)
        t_tm = (2 * np.sqrt(eps1) * kz1) / (eps2 * kz1 + eps1 * kz2)  # Исправленная формула
        
        return r_te, t_te, r_tm, t_tm

    def _s_matrix_multilayer(self, k_rho, k0, eps_list, d_list):
        """Вычисляет общую S-матрицу для многослойной структуры, возвращая коэффициенты отражения и прохождения."""
        # Кеширование результатов
        cache_key = (k_rho, tuple(eps_list), tuple(d_list))
        if cache_key in self._s_matrix_cache:
            return self._s_matrix_cache[cache_key]
        
        # Инициализация: начинаем с последнего интерфейса (подложка)
        R_te, R_tm = 0j, 0j
        T_te, T_tm = 1.0, 1.0
        
        # Идем от предпоследнего слоя к первому (снизу вверх)
        for i in range(len(eps_list)-2, -1, -1):
            eps1 = eps_list[i]
            eps2 = eps_list[i+1]
            d = d_list[i]
            
            # Волновые числа в текущем слое
            kz1 = np.lib.scimath.sqrt(eps1 * k0**2 - k_rho**2)
            kz2 = np.lib.scimath.sqrt(eps2 * k0**2 - k_rho**2)
            
            # Локальные коэффициенты на интерфейсе
            r_te12, t_te12, r_tm12, t_tm12 = self._s_matrix_layer(k_rho, k0, eps1, eps2)
            r_te21, t_te21, r_tm21, t_tm21 = self._s_matrix_layer(k_rho, k0, eps2, eps1)
            
            # Фазовый множитель в слое
            phase = np.exp(1j * kz2 * d)
            phase_sq = phase * phase
            
            # Обновление коэффициентов для TE поляризации
            denom_te = 1 - r_te21 * R_te * phase_sq
            R_te = r_te12 + (t_te12 * R_te * t_te21 * phase_sq) / denom_te
            T_te = (T_te * t_te21 * phase) / denom_te
            
            # Обновление коэффициентов для TM поляризации
            denom_tm = 1 - r_tm21 * R_tm * phase_sq
            R_tm = r_tm12 + (t_tm12 * R_tm * t_tm21 * phase_sq) / denom_tm
            T_tm = (T_tm * t_tm21 * phase) / denom_tm
            
        # Сохраняем в кеше
        self._s_matrix_cache[cache_key] = (R_te, R_tm, T_te, T_tm)
        return R_te, R_tm, T_te, T_tm

    @staticmethod
    @jit(nopython=True, fastmath=True, parallel=True)
    def _free_space_einc_kernel(X, Y, Z, rs, p_vec, m_vec, omega, mu0, k0, eta0):
        Nx, Ny, Nz = X.shape
        Einc_x = np.zeros((Nx, Ny, Nz), dtype=np.complex128)
        Einc_y = np.zeros((Nx, Ny, Nz), dtype=np.complex128)
        Einc_z = np.zeros((Nx, Ny, Nz), dtype=np.complex128)
        
        rx = X - rs[0]
        ry = Y - rs[1]
        rz = Z - rs[2]
        
        for i in prange(Nx):
            for j in prange(Ny):
                for k in prange(Nz):
                    r_val = max(1e-12, np.sqrt(rx[i,j,k]**2 + ry[i,j,k]**2 + rz[i,j,k]**2))
                    
                    # Векторное произведение и тензор Грина
                    rhat_x = rx[i,j,k] / r_val
                    rhat_y = ry[i,j,k] / r_val
                    rhat_z = rz[i,j,k] / r_val
                    
                    g = np.exp(1j * k0 * r_val) / (4 * np.pi * r_val)
                    gp = g * (1j * k0 - 1.0 / r_val)
                    gpp = g * (-k0**2 - 2j * k0 / r_val + 2.0 / (r_val**2))
                    
                    nnxx = rhat_x * rhat_x
                    nnyy = rhat_y * rhat_y
                    nnzz = rhat_z * rhat_z
                    nnxy = rhat_x * rhat_y
                    nnxz = rhat_x * rhat_z
                    nnyz = rhat_y * rhat_z
                    gp_over_r = gp / r_val
                    
                    d2g_xx = gpp * nnxx + gp_over_r * (1 - nnxx)
                    d2g_yy = gpp * nnyy + gp_over_r * (1 - nnyy)
                    d2g_zz = gpp * nnzz + gp_over_r * (1 - nnzz)
                    d2g_xy = gpp * nnxy - gp_over_r * nnxy
                    d2g_xz = gpp * nnxz - gp_over_r * nnxz
                    d2g_yz = gpp * nnyz - gp_over_r * nnyz
                    
                    invk2 = 1.0 / (k0**2)
                    Gxx = g + invk2 * d2g_xx
                    Gyy = g + invk2 * d2g_yy
                    Gzz = g + invk2 * d2g_zz
                    Gxy = invk2 * d2g_xy
                    Gxz = invk2 * d2g_xz
                    Gyz = invk2 * d2g_yz
                    
                    Gp_x = Gxx * p_vec[0] + Gxy * p_vec[1] + Gxz * p_vec[2]
                    Gp_y = Gxy * p_vec[0] + Gyy * p_vec[1] + Gyz * p_vec[2]
                    Gp_z = Gxz * p_vec[0] + Gyz * p_vec[1] + Gzz * p_vec[2]
                    
                    term1_x = 1j * omega * mu0 * Gp_x
                    term1_y = 1j * omega * mu0 * Gp_y
                    term1_z = 1j * omega * mu0 * Gp_z
                    
                    gradg_x = gp * rhat_x
                    gradg_y = gp * rhat_y
                    gradg_z = gp * rhat_z
                    
                    curl_gm_x = gradg_y * m_vec[2] - gradg_z * m_vec[1]
                    curl_gm_y = gradg_z * m_vec[0] - gradg_x * m_vec[2]
                    curl_gm_z = gradg_x * m_vec[1] - gradg_y * m_vec[0]
                    
                    Einc_x[i,j,k] = term1_x - curl_gm_x
                    Einc_y[i,j,k] = term1_y - curl_gm_y
                    Einc_z[i,j,k] = term1_z - curl_gm_z
        
        return Einc_x, Einc_y, Einc_z

    def _dipole_field_layered(self, rs, p_vec, obs_pos):
        """Вычисляет поле диполя в слоистой среде с учетом проходящей и отраженной волн."""
        # Параметры слоев: добавляем воздух сверху
        eps_list = [1.0] + self.layer_epsr
        d_list = [frac * self.L for frac in self.layer_fracs]
        
        # Волновое число
        k0 = self.k0
        
        # Координаты источника и точки наблюдения
        xs, ys, zs = rs
        xo, yo, zo = obs_pos
        
        # Разность координат в плоскости XY
        rho = np.sqrt((xo - xs)**2 + (yo - ys)**2)
        phi = np.arctan2(yo - ys, xo - xs)
        
        # Преобразование дипольного момента
        p_rho = p_vec[0] * np.cos(phi) + p_vec[1] * np.sin(phi)
        p_phi = -p_vec[0] * np.sin(phi) + p_vec[1] * np.cos(phi)
        p_z = p_vec[2]
        
        # Границы структуры
        z_top = self.z_max  # Верхняя граница структуры
        z_bottom = self.z_min  # Нижняя граница структуры
        
        # Параметры интегрирования
        k_rho_max = 10 * k0
        n_points = 1000
        k_rho_arr = np.linspace(0, k_rho_max, n_points)
        dk_rho = k_rho_arr[1] - k_rho_arr[0]
        
        E_rho, E_phi, E_z = 0j, 0j, 0j
        
        for k_rho in k_rho_arr:
            # Вычисляем S-матрицу для всей структуры
            R_te, R_tm, T_te, T_tm = self._s_matrix_multilayer(k_rho, k0, eps_list, d_list)
            
            # Волновое число в воздухе
            kz0 = np.lib.scimath.sqrt(k0**2 - k_rho**2)
            
            # Волновое число в подложке (последний слой)
            kz_sub = np.lib.scimath.sqrt(eps_list[-1] * k0**2 - k_rho**2)
            
            # Определяем положение точки наблюдения относительно структуры
            if zo >= z_top:  # Точка наблюдения выше структуры
                # Отраженная волна
                phase = np.exp(1j * kz0 * (zs + zo - 2*z_top))
                E_ref_te = R_te * p_phi * phase
                E_ref_tm = R_tm * p_rho * phase
                
                # Прямая волна (только если источник выше точки наблюдения)
                if zs > zo:
                    phase_dir = np.exp(1j * kz0 * np.abs(zs - zo))
                    E_ref_te += p_phi * phase_dir
                    E_ref_tm += p_rho * phase_dir
                    
                E_ref_x = E_ref_tm * np.cos(phi) - E_ref_te * np.sin(phi)
                E_ref_y = E_ref_tm * np.sin(phi) + E_ref_te * np.cos(phi)
                E_field = np.array([E_ref_x, E_ref_y, 0])
                
            elif zo <= z_bottom:  # Точка наблюдения ниже структуры
                # Прошедшая волна
                phase = np.exp(1j * kz0 * (zs - z_top)) * np.exp(1j * kz_sub * (zo - z_bottom))
                E_trans_te = T_te * p_phi * phase
                E_trans_tm = T_tm * p_rho * phase
                
                E_trans_x = E_trans_tm * np.cos(phi) - E_trans_te * np.sin(phi)
                E_trans_y = E_trans_tm * np.sin(phi) + E_trans_te * np.cos(phi)
                E_field = np.array([E_trans_x, E_trans_y, 0])
                
            else:  # Точка наблюдения внутри структуры
                # Упрощенный подход: используем только отраженную компоненту
                phase = np.exp(1j * kz0 * (zs + zo - 2*z_top))
                E_ref_te = R_te * p_phi * phase
                E_ref_tm = R_tm * p_rho * phase
                
                E_ref_x = E_ref_tm * np.cos(phi) - E_ref_te * np.sin(phi)
                E_ref_y = E_ref_tm * np.sin(phi) + E_ref_te * np.cos(phi)
                E_field = np.array([E_ref_x, E_ref_y, 0])
            
            # Интеграл Зоммерфельда
            weight = k_rho * jv(0, k_rho * rho) * dk_rho
            E_rho += E_field[0] * weight
            E_phi += E_field[1] * weight
        
        # Преобразование в декартовы координаты
        Ex = E_rho * np.cos(phi) - E_phi * np.sin(phi)
        Ey = E_rho * np.sin(phi) + E_phi * np.cos(phi)
        Ez = 0
        
        return np.array([Ex, Ey, Ez])
    def visualize_field_3d(self, E=None, component='amplitude', level=0.5, opacity=0.8, 
                          colormap='viridis', title='3D Field Visualization'):
        """
        Визуализирует поле в 3D с использованием изоповерхностей.
        
        Параметры:
            E (ndarray): Поле для визуализации (по умолчанию Esol)
            component (str): Компонента для визуализации ('x','y','z','amplitude')
            level (float): Уровень изоповерхности (относительный)
            opacity (float): Прозрачность поверхности (0-1)
            colormap (str): Цветовая карта
            title (str): Заголовок графика
        """
        if E is None:
            E = self._Esol
        assert E is not None, "Сначала решите задачу (solve)."
        
        # Выбор компоненты
        if component == 'x':
            field_data = np.abs(E[..., 0])
        elif component == 'y':
            field_data = np.abs(E[..., 1])
        elif component == 'z':
            field_data = np.abs(E[..., 2])
        else:  # amplitude
            field_data = np.sqrt(np.abs(E[..., 0])**2 + 
                                np.abs(E[..., 1])**2 + 
                                np.abs(E[..., 2])**2)
        
        # Проверка наличия данных
        if np.all(field_data == 0):
            print("Предупреждение: все значения поля равны нулю. Визуализация невозможна.")
            return None
        
        # Нормализация
        max_val = np.max(field_data)
        min_val = np.min(field_data)
        field_data_normalized = field_data / max_val
        
        # Автоматическая корректировка уровня
        if level < 0 or level > 1:
            level = 0.5  # Значение по умолчанию при некорректном вводе
        
        iso_level = level
        if iso_level < min_val/max_val or iso_level > 1.0:
            # Автоматическая установка уровня как 50% между min и max
            iso_level = 0.5 * (min_val/max_val + 1.0)
            print(f"Предупреждение: скорректирован уровень изоповерхности до {iso_level:.2f}")
        
        try:
            # Построение изоповерхности
            verts, faces, _, _ = measure.marching_cubes(
                field_data_normalized, 
                level=iso_level,
                spacing=(self.dx, self.dy, self.dz)
            )
        except Exception as e:
            print(f"Ошибка при построении изоповерхности: {e}")
            return None
        
        # Создание 3D фигуры
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Создание полигональной поверхности
        mesh = Poly3DCollection(verts[faces], alpha=opacity)
        mesh.set_edgecolor('k')
        mesh.set_facecolor(plt.cm.get_cmap(colormap)(0.5))
        ax.add_collection3d(mesh)
        
        # Настройка осей
        ax.set_xlim(0, field_data_normalized.shape[0])
        ax.set_ylim(0, field_data_normalized.shape[1])
        ax.set_zlim(0, field_data_normalized.shape[2])
        
        # Преобразование координат в физические единицы
        x_ticks = np.linspace(0, self.Nx-1, 5)
        y_ticks = np.linspace(0, self.Ny-1, 5)
        z_ticks = np.linspace(0, self.Nz-1, 5)
        
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.set_zticks(z_ticks)
        
        ax.set_xticklabels([f"{self.x_min + t*self.dx:.2e}" for t in x_ticks])
        ax.set_yticklabels([f"{self.y_min + t*self.dy:.2e}" for t in y_ticks])
        ax.set_zticklabels([f"{self.z_min + t*self.dz:.2e}" for t in z_ticks])
        
        ax.set_xlabel('X (м)')
        ax.set_ylabel('Y (м)')
        ax.set_zlabel('Z (м)')
        ax.set_title(f"{title} | Изоповерхность {level*100:.0f}%")
        
        # Добавление цветовой шкалы
        mappable = plt.cm.ScalarMappable(cmap=colormap)
        mappable.set_array(np.linspace(0, max_val, 100))
        cbar = fig.colorbar(mappable, ax=ax, shrink=0.5, aspect=10)
        cbar.set_label(f'|E| ({component})')
        
        return fig

    def visualize_field_slices(self, E=None, component='amplitude', slices=None, 
                             cmap='viridis', title='Field Slices'):
        """
        Визуализирует срезы поля в трех плоскостях.
        
        Параметры:
            E (ndarray): Поле для визуализации (по умолчанию Esol)
            component (str): Компонента для визуализации ('x','y','z','amplitude')
            slices (tuple): Координаты срезов (x0, y0, z0) или None для середины
            cmap (str): Цветовая карта
            title (str): Заголовок графика
        """
        if E is None:
            E = self._Esol
        assert E is not None, "Сначала решите задачу (solve)."
        
        # Выбор компоненты
        if component == 'x':
            field_data = np.abs(E[..., 0])
        elif component == 'y':
            field_data = np.abs(E[..., 1])
        elif component == 'z':
            field_data = np.abs(E[..., 2])
        else:  # amplitude
            field_data = np.sqrt(np.abs(E[..., 0])**2 + 
                                np.abs(E[..., 1])**2 + 
                                np.abs(E[..., 2])**2)
        
        # Определение срезов
        if slices is None:
            x0 = self.Nx // 2
            y0 = self.Ny // 2
            z0 = self.Nz // 2
        else:
            x0, y0, z0 = slices
        
        # Создание фигуры
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle(title, fontsize=16)
        
        # XY срез
        im0 = axes[0].imshow(field_data[:, :, z0].T, origin='lower', 
                           extent=[self.x_min, self.x_max, self.y_min, self.y_max],
                           cmap=cmap, aspect='equal')
        axes[0].set_title(f'XY срез (z = {self.z[z0]:.2e} м)')
        axes[0].set_xlabel('X (м)')
        axes[0].set_ylabel('Y (м)')
        plt.colorbar(im0, ax=axes[0], label=f'|E| ({component})')
        
        # XZ срез
        im1 = axes[1].imshow(field_data[:, y0, :].T, origin='lower',
                           extent=[self.x_min, self.x_max, self.z_min, self.z_max],
                           cmap=cmap, aspect='auto')
        axes[1].set_title(f'XZ срез (y = {self.y[y0]:.2e} м)')
        axes[1].set_xlabel('X (м)')
        axes[1].set_ylabel('Z (м)')
        plt.colorbar(im1, ax=axes[1], label=f'|E| ({component})')
        
        # YZ срез
        im2 = axes[2].imshow(field_data[x0, :, :].T, origin='lower',
                           extent=[self.y_min, self.y_max, self.z_min, self.z_max],
                           cmap=cmap, aspect='auto')
        axes[2].set_title(f'YZ срез (x = {self.x[x0]:.2e} м)')
        axes[2].set_xlabel('Y (м)')
        axes[2].set_ylabel('Z (м)')
        plt.colorbar(im2, ax=axes[2], label=f'|E| ({component})')
        
        plt.tight_layout()
        return fig
    
    def huygens_incident(self, rs, p_vec, m_vec, use_layered=False):
        if use_layered:
            E_layered = np.zeros((self.Nx, self.Ny, self.Nz, 3), dtype=np.complex128)
            
            if self.verbose:
                print("Calculating layered field (slow operation)...")
                iterable = tqdm(np.ndindex(self.Nx, self.Ny, self.Nz), 
                                total=self.Nx*self.Ny*self.Nz)
            else:
                iterable = np.ndindex(self.Nx, self.Ny, self.Nz)
                
            for i, j, k in iterable:
                pos = (self.X[i,j,k], self.Y[i,j,k], self.Z[i,j,k])
                E_layered[i,j,k] = self._dipole_field_layered(rs, p_vec, pos)
            
            self._Einc = E_layered
            return self._Einc
        
        # Оптимизированный расчет для свободного пространства с JIT
        p_vec = np.asarray(p_vec, dtype=np.complex128)
        m_vec = np.asarray(m_vec, dtype=np.complex128)
        
        Einc_x, Einc_y, Einc_z = self._free_space_einc_kernel(
            self.X, self.Y, self.Z, rs, p_vec, m_vec, 
            self.omega, self.mu0, self.k0, self.eta0
        )
        
        self._Einc = np.stack([Einc_x, Einc_y, Einc_z], axis=-1)
        return self._Einc

    def _build_spectral_green(self):
        kx = 2 * np.pi * np.fft.fftfreq(self.Npx, d=self.dx)
        ky = 2 * np.pi * np.fft.fftfreq(self.Npy, d=self.dy)
        kz = 2 * np.pi * np.fft.fftfreq(self.Npz, d=self.dz)
        self.KX, self.KY, self.KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        K2 = self.KX**2 + self.KY**2 + self.KZ**2

        den = (self.k0**2 - K2) + 1j * self.eta_reg
        
        # Исправление для k=0
        den = np.where(np.abs(den) < 1e-12, 1j * self.eta_reg, den)
        
        k_outer = np.zeros((self.Npx, self.Npy, self.Npz, 3, 3), dtype=np.complex128)
        k_outer[..., 0, 0] = self.KX * self.KX
        k_outer[..., 0, 1] = self.KX * self.KY
        k_outer[..., 0, 2] = self.KX * self.KZ
        k_outer[..., 1, 0] = self.KY * self.KX
        k_outer[..., 1, 1] = self.KY * self.KY
        k_outer[..., 1, 2] = self.KY * self.KZ
        k_outer[..., 2, 0] = self.KZ * self.KX
        k_outer[..., 2, 1] = self.KZ * self.KY
        k_outer[..., 2, 2] = self.KZ * self.KZ

        I3 = np.eye(3, dtype=np.complex128).reshape(1, 1, 1, 3, 3)
        proj = I3 - (1.0 / (self.k0**2 + 0j)) * k_outer
        self.Ghat = proj / den[..., None, None]
        
        # Корректная обработка k=0
        k0_mask = (np.abs(self.KX) < 1e-12) & (np.abs(self.KY) < 1e-12) & (np.abs(self.KZ) < 1e-12)
        self.Ghat[k0_mask] = np.eye(3) * (1/(3*self.eta_reg))

    def _pad_field(self, F):
        window_x = np.hamming(self.Nx)[:, None, None]
        window_y = np.hamming(self.Ny)[None, :, None]
        window_z = np.hamming(self.Nz)[None, None, :]
        window = window_x * window_y * window_z
        
        P = np.zeros((self.Npx, self.Npy, self.Npz, 3), dtype=np.complex128)
        P[:self.Nx, :self.Ny, :self.Nz, :] = F * window[..., None]
        return P

    def _unpad_field(self, P):
        return P[:self.Nx, :self.Ny, :self.Nz, :]

    def green_conv_fft(self, Jvec):
        P = self._pad_field(Jvec)
        
        # Векторизованное преобразование Фурье
        Jhat = np.fft.fftn(P, axes=(0,1,2))
        
        # Векторизованное умножение с использованием einsum
        Ehat = np.einsum('ijklm, ijkm->ijkl', self.Ghat, Jhat)
        
        # Обратное преобразование
        E = np.fft.ifftn(Ehat, axes=(0,1,2))
        return self._unpad_field(E) * self.dV

    def _apply_A(self, E_vec):
        E = E_vec.reshape(self.Nx, self.Ny, self.Nz, 3)
        P = self.chi[..., None] * E
        GE = self.green_conv_fft(P)
        out = E - (self.k0**2) * GE
        return out.ravel()

    def solve(self, tol=5e-4, restart=40, maxiter=100, method='gmres', n_born=3, born_tol=1e-3):
        if self._Einc is None:
            rs = np.array([0.0, 0.0, self.z_max + 0.6 * self.lam])
            p0 = 1.0
            p_vec = np.array([1.0, 0.0, 0.0]) * p0
            m_vec = np.array([0.0, -self.eta0 * p0, 0.0])
            
            self.huygens_incident(rs, p_vec, m_vec, use_layered=True)

        if method == 'born':
            if self.verbose:  
                contrast = np.max(np.abs(self.chi))
                print(f"Born iterations: n={n_born}, tol={born_tol}, max_contrast={contrast:.2f}")
            E = self._Einc.copy()
            for it in range(n_born):
                t0 = time.time()
                GE = self.green_conv_fft(self.chi[..., None] * E)
                E_new = self._Einc + (self.k0**2) * GE
                
                diff = np.linalg.norm(E_new - E) / np.linalg.norm(E)
                if self.verbose:
                    print(f"  Born {it+1}/{n_born}: |ΔE|={diff:.3e}, dt={time.time()-t0:.2f}s")
                if diff < born_tol:
                    if self.verbose:
                        print(f"  Born converged at iteration {it+1}")
                    break
                E = E_new
            self._Esol = E
            self._Esca = self._Esol - self._Einc
            return self._Esol, self.Esca, 0

        N_unknowns = 3 * self.Nx * self.Ny * self.Nz
        Aop = LinearOperator((N_unknowns, N_unknowns), matvec=self._apply_A, dtype=np.complex128)
        b = self._Einc.ravel()

        residuals = []

        def cb(res):
            residuals.append(res)
            if self.verbose and (len(residuals) == 1 or len(residuals) % 5 == 0):
                print(f"  GMRES iter {len(residuals)}: |r|={res:.3e}")

        if self.verbose:
            print("GMRES: solve A E = Einc ...")
        t0 = time.time()
        E_sol_vec, info = gmres(Aop, b, tol=tol, restart=restart, maxiter=maxiter,
                                callback=cb, callback_type='pr_norm')
        if self.verbose:
            print(f"GMRES done in {time.time()-t0:.2f}s, info={info} (0=success)")

        self._Esol = E_sol_vec.reshape(self.Nx, self.Ny, self.Nz, 3)
        self._Esca = self._Esol - self._Einc
        return self._Esol, self.Esca, info
    def farfield(self, E=None, n_theta=361, plane='xz'):
        if E is None:
            E = self._Esol
        assert E is not None, "Сначала решите задачу (solve)."
        
        theta = np.linspace(0.00001, 2*np.pi, n_theta)
        
        if plane.lower() == 'xz':
            shat = np.stack([np.sin(theta), np.zeros_like(theta), np.cos(theta)], axis=1)
        elif plane.lower() == 'yz':
            shat = np.stack([np.zeros_like(theta), np.sin(theta), np.cos(theta)], axis=1)
        else:
            raise ValueError("Недопустимая плоскость. Используйте 'xz' или 'yz'")

        # Используем оптимизированный kernel
        F_theta, F_phi = self._farfield_kernel(
            self.X, self.Y, self.Z, self.chi, E, self.k0, self.dV, shat
        )
        
        Fabs = np.sqrt(np.abs(F_theta)**2 + np.abs(F_phi)**2)
        return theta, F_theta, F_phi, Fabs
    @staticmethod
    @jit(nopython=True, fastmath=True, parallel=True)
    def _farfield_kernel(X, Y, Z, chi, E, k0, dV, shat):
        Nx, Ny, Nz = X.shape
        n_theta = shat.shape[0]
        
        F_theta = np.zeros(n_theta, dtype=np.complex128)
        F_phi = np.zeros(n_theta, dtype=np.complex128)
        
        #  константа (добавлен множитель 1j)
        const = dV * (k0**2) / (4 * np.pi) * 1j
        
        for i_theta in prange(n_theta):
            s = shat[i_theta].astype(np.complex128)
            F_vec = np.zeros(3, dtype=np.complex128)
            
            for i in range(Nx):
                for j in range(Ny):
                    for k in range(Nz):
                        if np.abs(chi[i, j, k]) > 1e-12:
                            rx = X[i, j, k]
                            ry = Y[i, j, k]
                            rz = Z[i, j, k]
                            
                            # скалярное произведение (s·r)
                            dot_rs = s[0]*rx + s[1]*ry + s[2]*rz
                            phase = np.exp(-1j * k0 * dot_rs)
                            
                            Px = chi[i, j, k] * E[i, j, k, 0]
                            Py = chi[i, j, k] * E[i, j, k, 1]
                            Pz = chi[i, j, k] * E[i, j, k, 2]
                            
                            F_vec[0] += phase * Px
                            F_vec[1] += phase * Py
                            F_vec[2] += phase * Pz
            
            # Проекция: F_vec = [I - ss^T] * F_vec
            s_dot_F = s[0]*F_vec[0] + s[1]*F_vec[1] + s[2]*F_vec[2]
            F_vec[0] -= s[0] * s_dot_F
            F_vec[1] -= s[1] * s_dot_F
            F_vec[2] -= s[2] * s_dot_F
            
            
            F_vec = const * F_vec
            
            # Вычисление ортов сферической системы
            s_real = np.array([np.real(s[0]), np.real(s[1]), np.real(s[2])])
            norm_s = np.sqrt(s_real[0]**2 + s_real[1]**2 + s_real[2]**2)
            if norm_s < 1e-12:
                continue
                
            cos_theta = s_real[2] / norm_s
            sin_theta = np.sqrt(1 - cos_theta**2)
            phi_i = np.arctan2(s_real[1], s_real[0])
            
            hat_theta = np.array([
                cos_theta * np.cos(phi_i),
                cos_theta * np.sin(phi_i),
                -sin_theta
            ])
            
            hat_phi = np.array([
                -np.sin(phi_i),
                np.cos(phi_i),
                0.0
            ])
            
            # Векторные проекции
            F_theta[i_theta] = F_vec[0]*hat_theta[0] + F_vec[1]*hat_theta[1] + F_vec[2]*hat_theta[2]
            F_phi[i_theta] = F_vec[0]*hat_phi[0] + F_vec[1]*hat_phi[1] + F_vec[2]*hat_phi[2]
        
        return F_theta, F_phi

    def plot_polarization_components(self, theta, F_theta, F_phi, plane='xz', 
                                    title='Компоненты поляризации', log_scale=False):
        plt.figure(figsize=(10, 6))
        
        abs_theta = np.abs(F_theta)
        abs_phi = np.abs(F_phi)
        max_val = max(np.max(abs_theta), np.max(abs_phi))
        norm_theta = abs_theta / max_val
        norm_phi = abs_phi / max_val
        
        theta_deg = np.degrees(theta)
        
        if log_scale:
            plt.semilogy(theta_deg, norm_theta, 'b-', linewidth=2, label=r'$F_\theta$ (Вертикальная)')
            plt.semilogy(theta_deg, norm_phi, 'r-', linewidth=2, label=r'$F_\phi$ (Горизонтальная)')
            plt.ylabel('Нормированная амплитуда (логарифм)')
        else:
            plt.plot(theta_deg, norm_theta, 'b-', linewidth=2, label=r'$F_\theta$ (Вертикальная)')
            plt.plot(theta_deg, norm_phi, 'r-', linewidth=2, label=r'$F_\phi$ (Горизонтальная)')
            plt.ylabel('Нормированная амплитуда')
        plt.xlabel(f'Угол θ в плоскости {plane.upper()} [градусы]')
        plt.title(f'{title} | Плоскость {plane.upper()}')
        plt.grid(True, which='both', linestyle='--', alpha=0.6)
        plt.legend(loc='best')
        plt.xlim(0, 360)
        plt.xticks(np.arange(0, 361, 45))
        plt.ylim(bottom=1e-4 if log_scale else 0)
        plt.tight_layout()
        
    def plot_slice_z(self, E=None, z0=0.0, cmap='magma', title='|E| срез по z'):
        if E is None:
            E = self._Esol
        assert E is not None, "Сначала решите задачу (solve)."
        iz = np.argmin(np.abs(self.z - z0))
        Eabs = np.sqrt(np.abs(E[:, :, iz, 0])**2 + 
                       np.abs(E[:, :, iz, 1])**2 + 
                       np.abs(E[:, :, iz, 2])**2)
        plt.figure(figsize=(6, 5))
        plt.title(f'{title}, z≈{self.z[iz]:.3g} м (index {iz})')
        plt.imshow(Eabs.T, origin='lower',
                extent=[self.x_min, self.x_max, self.y_min, self.y_max],
                aspect='equal', cmap=cmap)
        plt.colorbar(label='|E|')
        plt.xlabel('x, м')
        plt.ylabel('y, м')
        plt.tight_layout()

    def plot_farfield_polar(self, theta, Fabs_norm, title='Нормированная диаграмма (x–z)'):
        plt.figure(figsize=(6, 6))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, Fabs_norm, lw=2)
        ax.set_title(title)
        ax.set_rticks([0.25, 0.5, 0.75, 1.0])
        ax.grid(True)
        plt.tight_layout()
        
    def plot_polar_components(self, theta, F_theta, F_phi, plane='xz', title=None):
        max_val = max(np.max(np.abs(F_theta)), np.max(np.abs(F_phi)))
        F_theta_norm = np.abs(F_theta) / max_val
        F_phi_norm = np.abs(F_phi) / max_val
        E_theta_db = 20 * np.log10(F_theta_norm)
        E_phi_db = 20 * np.log10(F_phi_norm)
        
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
        ax.set_ylim(-40, 0)
        ax.grid(True, linestyle='--', alpha=0.7)
        
        ax.plot(theta, E_theta_db, 
                color='#1f77b4', linestyle='--', linewidth=1.5, 
                alpha=0.9, label=r'F$_\theta$($\theta$), дБ')
        ax.plot(theta, E_phi_db, 
                color='#ff7f0e', linestyle='-', linewidth=2.5, 
                alpha=0.9, label=r'F$_\varphi$($\theta$), дБ')
        
        ax.set_theta_zero_location('E')
        ax.set_theta_direction(-1)
        
        if title is None:
            title = f'Компоненты поляризации в плоскости {plane.upper()}'
        ax.set_title(title, fontsize=15, pad=20)
        
        ax.text(0, ax.get_rmax() + 7, r'$\theta$$\degree$', fontsize=14, ha='right', va='center')
        ax.legend(loc='upper right', fontsize=14, frameon=True, framealpha=0.95)
        
        plt.tight_layout()
        return fig
class ZCoordinateSweep:
    def __init__(self, solver_instance, z_start, z_stop, z_step, phi0=0.0, ref_db=0.0):
        assert z_stop > z_start, "z_stop должно быть больше z_start"
        assert z_step > 0, "z_step должен быть положительным"
        
        self.solver = solver_instance
        self.z_start = z_start
        self.z_stop = z_stop
        self.z_step = z_step
        self.phi0 = phi0
        self.ref_db = ref_db
        self.z_values = None
        self.results_theta = None
        self.results_phi = None
        self.theta_grid = None
        self.ref_amplitude = None

    def run_sweep(self, plane='xz', n_theta=361, method='born', n_born=3, born_tol=1e-3):
        n_points = int(np.ceil((self.z_stop - self.z_start) / self.z_step)) + 1
        self.z_values = np.linspace(self.z_start, self.z_stop, n_points)
        self.results_theta = []
        self.results_phi = []
        
        # Сохраняем сетку углов
        _, self.theta_grid, _, _ = self.solver.farfield(plane=plane, n_theta=n_theta)
        
        # Параллельное выполнение
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = []
            for z in self.z_values:
                futures.append(executor.submit(self._compute_for_z, z, plane, n_theta, method, n_born, born_tol))
            
            for future in tqdm(concurrent.futures.as_completed(futures), 
                              total=len(futures), desc="Sweeping z-coordinate"):
                F_theta, F_phi = future.result()
                self.results_theta.append(np.abs(F_theta))
                self.results_phi.append(np.abs(F_phi))
        
        self.results_theta = np.array(self.results_theta)
        self.results_phi = np.array(self.results_phi)
        combined = np.concatenate([self.results_theta.ravel(), self.results_phi.ravel()])
        self.ref_amplitude = np.max(combined) if np.max(combined) > 0 else 1.0

    def _compute_for_z(self, z, plane, n_theta, method, n_born, born_tol):
        rs = [0.0, 0.0, z]
        p_vec = [1.0, 0.0, 0.0]
        m_vec = [0.0, -self.solver.eta0, 0.0]
        self.solver.huygens_incident(rs, p_vec, m_vec, use_layered=True)
        self.solver.solve(method=method, n_born=n_born, born_tol=born_tol)
        _, F_theta, F_phi, _ = self.solver.farfield(plane=plane, n_theta=n_theta)
        return F_theta, F_phi

    def to_db(self, data):
        data = np.where(data <= 0, 1e-20, data)
        return 20 * np.log10(data / self.ref_amplitude) + self.ref_db

    def plot_heatmap(self, component='both', output_file='heatmap_z_db.png', 
                    cmap='turbo', vmin=-40, vmax=0, figsize=(12, 8)):
        if self.results_theta is None:
            raise RuntimeError("Сначала выполните run_sweep()")
        
        if component == 'theta':
            data = self.to_db(self.results_theta)
            title_suffix = r'$E_\theta$'
        elif component == 'phi':
            data = self.to_db(self.results_phi)
            title_suffix = r'$E_\phi$'
        else:
            combined = np.sqrt(self.results_theta**2 + self.results_phi**2)
            data = self.to_db(combined)
            title_suffix = r'$\sqrt{E_\theta^2 + E_\phi^2}$'
        
        theta_deg = np.degrees(self.theta_grid)
        Z, Theta = np.meshgrid(self.z_values, theta_deg)
        data = data.T
        
        fig, ax = plt.subplots(figsize=figsize)
        pc = ax.pcolormesh(Theta, Z, data, shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        
        ax.set_xlabel('Угол θ [градусы]', fontsize=12)
        ax.set_ylabel('Координата z [м]', fontsize=12)
        ax.set_title(f'Зависимость поля от z и θ | {title_suffix} [дБ]', fontsize=14)
        ax.grid(True, linestyle='--', alpha=0.7)
        
        contour = ax.contour(Theta, Z, data, levels=np.arange(vmin, vmax+1, 5),
                            colors='k', linewidths=0.5, alpha=0.5)
        ax.clabel(contour, inline=True, fontsize=8, fmt='%d dB')
        
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('Амплитуда поля, дБ', fontsize=12)
        
        ax.xaxis.set_major_locator(plt.MultipleLocator(30))
        ax.yaxis.set_major_locator(plt.AutoLocator())
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return fig

    def plot_polar_heatmap(self, component='both', output_file='polar_heatmap_z_db.png',
                          cmap='turbo', vmin=-40, vmax=0, figsize=(10, 8)):
        if self.results_theta is None:
            raise RuntimeError("Сначала выполните run_sweep()")
            
        if component == 'theta':
            data = self.to_db(self.results_theta)
            title_suffix = r'$E_\theta$'
        elif component == 'phi':
            data = self.to_db(self.results_phi)
            title_suffix = r'$E_\phi$'
        else:
            combined = np.sqrt(self.results_theta**2 + self.results_phi**2)
            data = self.to_db(combined)
            title_suffix = r'$\sqrt{E_\theta^2 + E_\phi^2}$'
        
        theta_rad = self.theta_grid
        R, Theta = np.meshgrid(self.z_values, theta_rad)
        data = data.T
        
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'polar'})
        pc = ax.pcolormesh(Theta, R, data, shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        
        ax.set_theta_zero_location('E')
        ax.set_theta_direction(-1)
        ax.set_title(f'Полярная карта поля | {title_suffix} [дБ]', fontsize=14, pad=20)
        
        cbar = fig.colorbar(pc, ax=ax, pad=0.1)
        cbar.set_label('Амплитуда поля, дБ', fontsize=12)
        
        ax.text(0, np.max(R)*1.1, 'z [м]', fontsize=12, ha='center')
        
        contour = ax.contour(Theta, R, data, levels=np.arange(vmin, vmax+1, 5),
                           colors='w', linewidths=0.5, alpha=0.7)
        ax.clabel(contour, inline=True, fontsize=8, fmt='%d dB')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return fig

if __name__ == "__main__":
    # Параметры задачи
    freq = 2.4e9  # Гц
    kappa_a =1.0 * np.pi
    kappa_b = 1.0 * np.pi
    kappa_L = 1.0 * np.pi

    # Физически осмысленные параметры слоев
    layer_fracs = [1]  # доли толщины
    layer_epsr = [5]  # Комплексные εr
    grid_shape = (20, 20, 10)  # узлы по (x,y,z)
    padding_factor = 2

    solver = MultilayerCuboidVIESolver(
        freq_hz=freq,
        kappa_a=kappa_a, 
        kappa_b=kappa_b, 
        kappa_L=kappa_L,
        layer_fracs=layer_fracs, 
        layer_epsr=layer_epsr,
        grid_shape=grid_shape,
        padding_factor=padding_factor,
        eta_reg_frac=1e-3,
        verbose=True  
    )
    
    # Решаем VIE
    Esol, Esca, info = solver.solve(method='born', n_born=3, born_tol=1e-3)

    # Визуализация в ближней зоне
    solver.plot_slice_z(E=Esol, z0=0.1)
    plt.show()

    # Дальняя зона в плоскости XZ
    theta_xz, F_theta_xz, F_phi_xz, Fabs_xz = solver.farfield(plane='xz')
    
    # Дальняя зона в плоскости YZ
    theta_yz, F_theta_yz, F_phi_yz, Fabs_yz = solver.farfield(plane='yz')

    # Визуализация в полярных координатах
    fig_xz = solver.plot_polar_components(theta_xz, F_theta_xz, F_phi_xz, 
                                         plane='xz', 
                                         title='Поляризационные компоненты в плоскости XZ')
    
    fig_yz = solver.plot_polar_components(theta_yz, F_theta_yz, F_phi_yz, 
                                         plane='yz', 
                                         title='Поляризационные компоненты в плоскости YZ')
    print(F_theta_xz, F_phi_xz)
    plt.show()
    # После решения задачи
    solver.solve(method='born')
    solver.visualize_field_3d(component='amplitude', level=0.3, opacity=0.7)
    plt.show()
    solver.visualize_field_slices(component='z', cmap='hot')
    plt.show()
    # Анализ по z-координате
    z_sweep = ZCoordinateSweep(
        solver_instance=solver,
        z_start=0*np.pi,
        z_stop=10*np.pi,
        z_step=0.5*np.pi,
        ref_db=0.0
    )

    #z_sweep.run_sweep(plane='xz', method='born')
    #z_sweep.plot_heatmap()
    #plt.show()