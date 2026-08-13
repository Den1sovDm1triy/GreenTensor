# SPDX-License-Identifier: MIT

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

class PlaneWaveSolver: 
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
        
        # Границы слоёв по z
        z_bounds = [self.z_min]
        acc = self.z_min
        for f in self.layer_fracs:
            acc += f * self.L
            z_bounds.append(min(acc, self.z_max))
        
        # Создание 3D-маски
        mask_xy = (np.abs(self.X) <= self.a/2) & (np.abs(self.Y) <= self.b/2)
        
        for m in range(self.N_layers):
            zl, zr = z_bounds[m], z_bounds[m+1]
            if m == self.N_layers - 1:
                mask_z = (self.Z >= zl) & (self.Z <= zr)
            else:
                mask_z = (self.Z >= zl) & (self.Z < zr)
            
            mask_layer = mask_xy & mask_z
            chi[mask_layer] = (self.layer_epsr[m] - 1.0)
        
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
        
        r_te = (kz1 - kz2) / (kz1 + kz2)
        t_te = 2 * kz1 / (kz1 + kz2)
        r_tm = (eps2 * kz1 - eps1 * kz2) / (eps2 * kz1 + eps1 * kz2)
        t_tm = (2 * np.sqrt(eps1) * kz1) / (eps2 * kz1 + eps1 * kz2)
        
        return r_te, t_te, r_tm, t_tm

    def _s_matrix_multilayer(self, k_rho, k0, eps_list, d_list):
        """Вычисляет общую S-матрицу для многослойной структуры."""
        cache_key = (k_rho, tuple(eps_list), tuple(d_list))
        if cache_key in self._s_matrix_cache:
            return self._s_matrix_cache[cache_key]
        
        R_te, R_tm = 0j, 0j
        T_te, T_tm = 1.0, 1.0
        
        for i in range(len(eps_list)-2, -1, -1):
            eps1 = eps_list[i]
            eps2 = eps_list[i+1]
            d = d_list[i]
            
            kz1 = np.lib.scimath.sqrt(eps1 * k0**2 - k_rho**2)
            kz2 = np.lib.scimath.sqrt(eps2 * k0**2 - k_rho**2)
            
            r_te12, t_te12, r_tm12, t_tm12 = self._s_matrix_layer(k_rho, k0, eps1, eps2)
            r_te21, t_te21, r_tm21, t_tm21 = self._s_matrix_layer(k_rho, k0, eps2, eps1)
            
            phase = np.exp(1j * kz2 * d)
            phase_sq = phase * phase
            
            denom_te = 1 - r_te21 * R_te * phase_sq
            R_te = r_te12 + (t_te12 * R_te * t_te21 * phase_sq) / denom_te
            T_te = (T_te * t_te21 * phase) / denom_te
            
            denom_tm = 1 - r_tm21 * R_tm * phase_sq
            R_tm = r_tm12 + (t_tm12 * R_tm * t_tm21 * phase_sq) / denom_tm
            T_tm = (T_tm * t_tm21 * phase) / denom_tm
            
        self._s_matrix_cache[cache_key] = (R_te, R_tm, T_te, T_tm)
        return R_te, R_tm, T_te, T_tm

    def set_plane_wave(self, theta_inc=0.0, phi_inc=0.0, E0_mag=1.0, polarization='x'):
        """Устанавливает падающую плоскую волну с заданными параметрами."""
        # Рассчитываем волновой вектор
        kx = self.k0 * np.sin(theta_inc) * np.cos(phi_inc)
        ky = self.k0 * np.sin(theta_inc) * np.sin(phi_inc)
        kz = -self.k0 * np.cos(theta_inc)  # отрицательный, так как волна идет вниз

        # Вектор электрического поля (перпендикулярен волновому вектору)
        if polarization == 'x':
            E0 = np.array([E0_mag, 0, 0], dtype=np.complex128)
        elif polarization == 'y':
            E0 = np.array([0, E0_mag, 0], dtype=np.complex128)
        else:
            raise ValueError("Поляризация должна быть 'x' или 'y'")

        # Корректировка вектора E0 для перпендикулярности
        k_inc_vec = np.array([kx, ky, kz])
        k_norm = np.linalg.norm(k_inc_vec)
        if k_norm > 1e-10:
            k_hat = k_inc_vec / k_norm
            E0_parallel = np.dot(E0, k_hat) * k_hat
            E0 = E0 - E0_parallel
            if self.verbose:
                print("Скорректирован E0 для перпендикулярности k_inc.")

        # Фазовый множитель
        phase = np.exp(-1j * (kx * self.X + ky * self.Y + kz * self.Z))

        # Падающее поле
        E_inc_x = E0[0] * phase
        E_inc_y = E0[1] * phase
        E_inc_z = E0[2] * phase

        self._Einc = np.stack([E_inc_x, E_inc_y, E_inc_z], axis=-1)

    def _build_spectral_green(self):
        kx = 2 * np.pi * np.fft.fftfreq(self.Npx, d=self.dx)
        ky = 2 * np.pi * np.fft.fftfreq(self.Npy, d=self.dy)
        kz = 2 * np.pi * np.fft.fftfreq(self.Npz, d=self.dz)
        self.KX, self.KY, self.KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        K2 = self.KX**2 + self.KY**2 + self.KZ**2

        den = (self.k0**2 - K2) + 1j * self.eta_reg
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
        Jhat = np.fft.fftn(P, axes=(0,1,2))
        Ehat = np.einsum('ijklm, ijkm->ijkl', self.Ghat, Jhat)
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
            # Установка плоской волны по умолчанию при отсутствии падающего поля
            self.set_plane_wave()

        if method == 'born':
            contrast = np.max(np.abs(self.chi))
            if self.verbose:  
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
        
        theta = np.linspace(0.0, 2*np.pi, n_theta)
        
        if plane.lower() == 'xz':
            shat = np.stack([np.sin(theta), np.zeros_like(theta), np.cos(theta)], axis=1)
        elif plane.lower() == 'yz':
            shat = np.stack([np.zeros_like(theta), np.sin(theta), np.cos(theta)], axis=1)
        else:
            raise ValueError("Недопустимая плоскость. Используйте 'xz' или 'yz'")

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
        
        const = dV * (k0**2) / (4*np.pi)
        
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
                            
                            dot_rs = rx*s[0] + ry*s[1] + rz*s[2]
                            phase = np.exp(-1j * k0 * dot_rs)
                            
                            Px = chi[i, j, k] * E[i, j, k, 0]
                            Py = chi[i, j, k] * E[i, j, k, 1]
                            Pz = chi[i, j, k] * E[i, j, k, 2]
                            
                            F_vec[0] += phase * Px
                            F_vec[1] += phase * Py
                            F_vec[2] += phase * Pz
            
            s_dot_F = s[0]*F_vec[0] + s[1]*F_vec[1] + s[2]*F_vec[2]
            F_vec[0] = F_vec[0] - s[0] * s_dot_F
            F_vec[1] = F_vec[1] - s[1] * s_dot_F
            F_vec[2] = F_vec[2] - s[2] * s_dot_F
            
            F_vec = const * 1j * k0 * F_vec
            
            s_real = np.array([np.real(s[0]), np.real(s[1]), np.real(s[2])])
            norm_s = np.sqrt(s_real[0]**2 + s_real[1]**2 + s_real[2]**2)
            cos_theta = s_real[2] / norm_s
            sin_theta = np.sqrt(1 - cos_theta**2)
            phi_i = np.arctan2(s_real[1], s_real[0])
            
            hat_theta_x = cos_theta * np.cos(phi_i)
            hat_theta_y = cos_theta * np.sin(phi_i)
            hat_theta_z = -sin_theta
            
            hat_phi_x = -np.sin(phi_i)
            hat_phi_y = np.cos(phi_i)
            hat_phi_z = 0.0
            
            F_theta[i_theta] = F_vec[0]*hat_theta_x + F_vec[1]*hat_theta_y + F_vec[2]*hat_theta_z
            F_phi[i_theta] = F_vec[0]*hat_phi_x + F_vec[1]*hat_phi_y + F_vec[2]*hat_phi_z
        
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
            plt.semilogy(theta_deg, norm_theta, 'b-', linewidth=2, label=r'$E_\theta$ (Вертикальная)')
            plt.semilogy(theta_deg, norm_phi, 'r-', linewidth=2, label=r'$E_\phi$ (Горизонтальная)')
            plt.ylabel('Нормированная амплитуда (логарифм)')
        else:
            plt.plot(theta_deg, norm_theta, 'b-', linewidth=2, label=r'$E_\theta$ (Вертикальная)')
            plt.plot(theta_deg, norm_phi, 'r-', linewidth=2, label=r'$E_\phi$ (Горизонтальная)')
            plt.ylabel('Нормированная амплитуда')
        plt.xlabel(f'Угол θ в плоскости {plane.upper()} [градусы]')
        plt.title(f'{title} | Плоскость {plane.upper()}')
        plt.grid(True, which='both', linestyle='--', alpha=0.6)
        plt.legend(loc='best')
        plt.xlim(0, 360)
        plt.xticks(np.arange(0, 361, 45))
        plt.ylim(bottom=1e-4 if log_scale else 0)
        plt.tight_layout()
        plt.show()
    
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
                alpha=0.9, label=r'E$_\theta$($\theta$), дБ')
        ax.plot(theta, E_phi_db, 
                color='#ff7f0e', linestyle='-', linewidth=2.5, 
                alpha=0.9, label=r'E$_\varphi$($\theta$), дБ')
        
        ax.set_theta_zero_location('E')
        ax.set_theta_direction(1)
        
        if title is None:
            title = f'Компоненты поляризации в плоскости {plane.upper()}'
        ax.set_title(title, fontsize=15, pad=20)
        
        ax.text(0, ax.get_rmax() + 7, r'$\theta$$\degree$', fontsize=14, ha='right', va='center')
        ax.legend(loc='upper right', fontsize=14, frameon=True, framealpha=0.95)
        
        plt.tight_layout()
        return fig
    
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
        plt.show()
    # Методы визуализации (plot_slice_z, plot_polar_components и др. остаются без изменений)
    # ... (реализовать по аналогии с оригинальным классом)

# Пример использования
if __name__ == "__main__":
    freq = 2.4e9  # Гц
    kappa_a = 10 * np.pi
    kappa_b = 10 * np.pi
    kappa_L = 10 * np.pi
    layer_fracs = [1] 
    layer_epsr = [5] 
    grid_shape = (20, 20, 20)

    # Создаем один решатель для обеих компонент
    solver = PlaneWaveSolver(
        freq_hz=freq,
        kappa_a=kappa_a, 
        kappa_b=kappa_b, 
        kappa_L=kappa_L,
        layer_fracs=layer_fracs, 
        layer_epsr=layer_epsr,
        grid_shape=grid_shape,
        verbose=True  
    )
    
    # Первый расчет для X-поляризации
    solver.set_plane_wave(theta_inc=np.radians(0), phi_inc=0, E0_mag=1.0, polarization='x')
    solver.solve(method='born', n_born=4, born_tol=1e-3)
    theta, F_theta_x, F_phi_x, _ = solver.farfield(plane='xz')
    
    # Второй расчет для Y-поляризации
    solver.set_plane_wave(theta_inc=np.radians(0), phi_inc=0, E0_mag=1.0, polarization='y')
    solver.solve(method='born', n_born=4, born_tol=1e-3)
    _, F_theta_y, F_phi_y, _ = solver.farfield(plane='xz')
    
    # Суммируем компоненты от обеих поляризаций
    F_theta_total = F_theta_x + F_theta_y
    F_phi_total = F_phi_x + F_phi_y
    
    # Строим суммарную диаграмму
    solver.plot_polar_components(theta, F_theta_total, F_phi_total, 
                                plane='xz', 
                                title='Суммарное поле в XZ плоскости')


    plt.show()