import math
import cmath
import scipy
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

class RCSCalculator:
    def __init__(self, k0=0.5, toch=10, n=3, phi=math.pi/2, a=None, eps=None, miy=None):
        """
        Инициализация калькулятора RCS
        
        Параметры:
        k0 - волновое число
        toch - количество учитываемых членов ряда
        n - число слоев (последний слой - воздух)
        phi - азимутальный угол (по умолчанию pi)
        a, eps, miy - параметры слоев (радиусы, диэлектрические проницаемости, магнитные проницаемости)
        """
        self.k0 = k0
        self.toch = toch
        self.n = n
        self.phi = phi
        
        # Параметры материалов по умолчанию (3 слоя)
        self.a = a if a is not None else [0.1, 0.2, 1]
        self.eps = eps if eps is not None else [-1.7e7j, -1.7e7j, 1]
        self.miy = miy if miy is not None else [1, 1, 1]
        
        # Инициализация всех необходимых массивов
        self._initialize_arrays()
        
    def _initialize_arrays(self):
        """Инициализация всех массивов и переменных"""
        # Углы для расчета
        self.teta_start = 0.01
        self.teta_stop = 360
        self.step = math.pi/180
        self.teta_diap = abs(self.teta_stop) - abs(self.teta_start)
        self.steps = int(((self.teta_diap * (math.pi/180)) / self.step))
        self.teta = np.zeros(int(self.steps))
        self.cos_teta = np.zeros(int(self.steps))
        
        # Переменные среды
        self.alfa = np.zeros(self.n, dtype=complex)
        self.beta = np.zeros(self.n, dtype=complex)
        self.sigma = np.zeros(self.n, dtype=complex)
        self.k = np.zeros((self.n, self.n), dtype=complex)
        
        # Функции Бесселя, Неймана и др.
        self.J = np.zeros(self.toch, dtype=complex)
        self.Jpr = np.zeros(self.toch, dtype=complex)
        self.N = np.zeros(self.toch, dtype=complex)
        self.Npr = np.zeros(self.toch, dtype=complex)
        self.C = np.zeros((self.toch, len(self.sigma)-1), dtype=complex)
        self.Cpr = np.zeros((self.toch, len(self.sigma)-1), dtype=complex)
        self.S = np.zeros((self.toch, len(self.sigma)-1), dtype=complex)
        self.Spr = np.zeros((self.toch, len(self.sigma)-1), dtype=complex)
        
        # Импедансы и адмитансы
        self.Z = np.zeros((self.toch, len(self.a)), dtype=complex)
        self.Y = np.zeros((self.toch, len(self.a)), dtype=complex)
        
        # Модифицированные функции
        self.mJ = np.zeros(self.toch, dtype=complex)
        self.mJpr = np.zeros(self.toch, dtype=complex)
        self.mH = np.zeros(self.toch, dtype=complex)
        self.mHpr = np.zeros(self.toch, dtype=complex)
        
        # Коэффициенты рассеяния
        self.Mn = np.zeros(self.toch, dtype=complex)
        self.Nn = np.zeros(self.toch, dtype=complex)
        
        # Поля
        self.E_kp = np.zeros((self.toch, self.steps), dtype=complex)
        self.E_op = np.zeros((self.toch, self.steps), dtype=complex)
        self.S_teta = np.zeros((self.toch, self.steps), dtype=complex)
        self.S_phi = np.zeros((self.toch, self.steps), dtype=complex)
        self.E_teta = np.zeros((1, self.steps), dtype=complex)
        self.E_phi = np.zeros((1, self.steps), dtype=complex)
        
    def calculate_medium_parameters(self):
        """Расчет параметров среды"""
        for i in range(self.n):
            self.alfa[i] = cmath.atan(self.eps[i].imag / self.eps[i].real) if self.eps[i].real != 0 else cmath.pi/2
            self.beta[i] = math.atan(self.miy[i].imag / self.miy[i].real)
            self.sigma[i] = cmath.sqrt(abs(self.eps[i]) * abs(self.miy[i]))
    
    def calculate_k_coefficients(self):
        """Расчет коэффициентов k"""
        j = 0
        for i in range(self.n):
            self.k[i][j] = self.k0 * self.a[i] * self.sigma[j]
            if j < self.n - 1:
                j += 1
                self.k[i][j] = self.k0 * self.a[i] * self.sigma[j]
    
    def Jfunc(self, i, j1, j2):
        """Функция Бесселя первого рода"""
        nu = i + 1
        return scipy.special.jv(nu + 0.5, self.k[j1][j2]) * cmath.sqrt(self.k[j1][j2] * math.pi/2)
    
    def Jprfunc(self, i, j1, j2, tie):
        """Производная функции Бесселя первого рода"""
        nu = i + 1
        if not tie:
            return ((nu / (2 * nu + 1)) * 
               scipy.special.jv(nu - 0.5, self.k[j1][j2]) * 
               cmath.sqrt(self.k[j1][j2] * math.pi/2) -
               ((nu + 1) / (2 * nu + 1)) * 
               scipy.special.jv(nu + 1.5, self.k[j1][j2]) * 
               cmath.sqrt(self.k[j1][j2] * math.pi/2) +
               (self.J[i] / self.k[j1][j2]))
        else:
           return ((nu / (2 * nu + 1)) * 
               (scipy.special.jv(nu - 0.5, self.k[j1][j2]) * 
               cmath.sqrt(self.k[j1][j2] * math.pi/2) / self.k[j1][j2]) * self.k[j1][j2] -
               ((nu + 1) / (2 * nu + 1)) * 
               (scipy.special.jv(nu + 1.5, self.k[j1][j2]) * 
               cmath.sqrt(self.k[j1][j2] * math.pi/2) / self.k[j1][j2]) * self.k[j1][j2] +
               (scipy.special.jv(nu + 0.5, self.k[j1][j2]) * 
               cmath.sqrt(self.k[j1][j2] * math.pi/2) / self.k[j1][j2]))
    
    def Nfunc(self, i, j1, j2):
        """Функция Неймана"""
        return scipy.special.yv((i+1) + 0.5, self.k[j1][j2]) * cmath.sqrt(self.k[j1][j2] * math.pi/2)
    
    def Nprfunc(self, i, j1, j2, tie):
        """Производная функции Неймана"""
        nu = i + 1
        if not tie:
            return ((nu / (2 * nu + 1)) * scipy.special.yv(nu - 0.5, self.k[j1][j2]) * cmath.sqrt(self.k[j1][j2] * math.pi/2) -
                   ((nu + 1) / (2 * nu + 1)) * scipy.special.yv(nu + 1.5, self.k[j1][j2]) * cmath.sqrt(self.k[j1][j2] * math.pi/2) +
                   (self.Nfunc(i, j1, j2) / self.k[j1][j2]))
        else:
            return ((nu / (2 * nu + 1)) * (scipy.special.yv(nu - 0.5, self.k[j1][j2]) * cmath.sqrt(self.k[j1][j2] * math.pi/2) / self.k[j1][j2]) * self.k[j1][j2] -
                   ((nu + 1) / (2 * nu + 1)) * (scipy.special.yv(nu + 1.5, self.k[j1][j2]) * cmath.sqrt(self.k[j1][j2] * math.pi/2) / self.k[j1][j2]) * self.k[j1][j2] +
                   (self.Nfunc(i, j1, j2) / self.k[j1][j2]))
    
    def calculate_bessel_functions(self):
        """Расчет функций Бесселя и Неймана"""
        for i in range(self.toch):
            self.J[i] = self.Jfunc(i, 0, 0)
            self.Jpr[i] = self.Jprfunc(i, 0, 0, False)
            self.N[i] = self.Nfunc(i, 0, 0)
            self.Npr[i] = self.Nprfunc(i, 0, 0, False)
    
    def calculate_CS_functions(self):
        """Расчет функций C, S и их производных"""
        for i in range(self.toch-1):
            for j in range(len(self.sigma)-1):
                self.C[i][j] = (self.Jfunc(i, j+1, j+1) * self.Nprfunc(i, j, j+1, True) - 
                               self.Nfunc(i, j+1, j+1) * self.Jprfunc(i, j, j+1, True))
                self.Cpr[i][j] = (self.Jprfunc(i, j+1, j+1, True) * self.Nprfunc(i, j, j+1, True) - 
                                 self.Nprfunc(i, j+1, j+1, True) * self.Jprfunc(i, j, j+1, True))
                self.S[i][j] = (self.Nfunc(i, j+1, j+1) * self.Jfunc(i, j, j+1) - 
                               self.Jfunc(i, j+1, j+1) * self.Nfunc(i, j, j+1))
                self.Spr[i][j] = (self.Nprfunc(i, j+1, j+1, True) * self.Jfunc(i, j, j+1) - 
                                 self.Jprfunc(i, j+1, j+1, True) * self.Nfunc(i, j, j+1))
    
    def calculate_impedances(self):
        """Расчет импедансов и адмитансов"""
        for i in range(self.toch - 1):
            for h in range(len(self.a)):
                if h == 0:
                    numerator = cmath.exp(self.alfa[1] * 1j) * abs(self.eps[1])
                    denominator = cmath.exp(self.alfa[0] * 1j) * abs(self.eps[0])
                    self.Z[i][h] = cmath.sqrt(numerator / denominator) * (self.Jpr[i] / self.J[i])
                
                    numerator = cmath.exp(self.alfa[0] * 1j) * abs(self.eps[0])
                    denominator = cmath.exp(self.alfa[1] * 1j) * abs(self.eps[1])
                    self.Y[i][h] = cmath.sqrt(numerator / denominator) * (self.Jpr[i] / self.J[i])
                
                elif h == (len(self.a) - 1):
                    numerator = cmath.exp(self.alfa[h+1] * 1j) * abs(self.eps[h+1])
                    denominator = cmath.exp(self.alfa[h] * 1j) * abs(self.eps[h])
                    sqrt_part = cmath.sqrt(numerator / denominator)
                
                    term1 = self.Cpr[i][h-1] + self.Z[i][h-1] * self.Spr[i][h-1]
                    term2 = self.C[i][h-1] + self.Z[i][h-1] * self.S[i][h-1]
                    self.Z[i][h] = sqrt_part * (term1 / term2) / 2
                
                    numerator = cmath.exp(self.alfa[h] * 1j) * abs(self.eps[h])
                    denominator = cmath.exp(self.alfa[h+1] * 1j) * abs(self.eps[h+1])
                    sqrt_part = cmath.sqrt(numerator / denominator)
                    term1 = self.Cpr[i][h-1] + self.Y[i][h-1] * self.Spr[i][h-1]
                    term2 = self.C[i][h-1] + self.Y[i][h-1] * self.S[i][h-1]
                    self.Y[i][h] = sqrt_part * (term1 / term2) * 2
                
                else:
                    numerator = cmath.exp(self.alfa[h+1] * 1j) * abs(self.eps[h+1])
                    denominator = cmath.exp(self.alfa[h] * 1j) * abs(self.eps[h])
                    sqrt_part = cmath.sqrt(numerator / denominator)
                
                    term1 = self.Cpr[i][h-1] + self.Z[i][h-1] * self.Spr[i][h-1]
                    term2 = self.C[i][h-1] + self.Z[i][h-1] * self.S[i][h-1]
                    self.Z[i][h] = sqrt_part * (term1 / term2)

                    numerator = cmath.exp(self.alfa[h] * 1j) * abs(self.eps[h])
                    denominator = cmath.exp(self.alfa[h+1] * 1j) * abs(self.eps[h+1])
                    sqrt_part = cmath.sqrt(numerator / denominator)

                    term1 = self.Cpr[i][h-1] + self.Y[i][h-1] * self.Spr[i][h-1]
                    term2 = self.C[i][h-1] + self.Y[i][h-1] * self.S[i][h-1]
                    self.Y[i][h] = sqrt_part * (term1 / term2)

    def Hfunc(self, i, k0):
        """Функция Ханкеля второго рода"""
        nu = i + 1
        return scipy.special.hankel1(nu + 0.5, k0) * math.sqrt(k0 * math.pi/2)
    
    def Hprfunc(self, i, k0):
        """Производная функции Ханкеля второго рода"""
        nu = i + 1
        term1 = (nu / (2 * nu + 1)) * (scipy.special.hankel1(nu - 0.5, k0) * cmath.sqrt(k0 * math.pi/2)) / k0
        term2 = ((nu + 1) / (2 * nu + 1)) * (scipy.special.hankel1(nu + 1.5, k0) * cmath.sqrt(k0 * math.pi/2)) / k0
        term3 = (scipy.special.hankel1(nu + 0.5, k0) * cmath.sqrt(k0 * math.pi/2)) / k0
    
        return (term1 * k0 - term2 * k0 + term3)
    
    def calculate_modified_functions(self):
        """Расчет модифицированных функций"""
        k1 = self.k0
        k00 = self.k[0][0]
        self.k[0][0] = self.k0
        
        for i in range(self.toch):
            self.mJ[i] = self.Jfunc(i, 0, 0)
            self.mJpr[i] = self.Jprfunc(i, 0, 0, True)
            self.mH[i] = self.Hfunc(i, self.k0)
            self.mHpr[i] = self.Hprfunc(i, self.k0)
            
        self.k0 = k1
        self.k[0][0] = k00
    
    def calculate_scattering_coefficients(self):
        """Расчет коэффициентов рассеяния"""
        for i in range(self.toch):
            nu = i + 1
            self.Mn[i] = (self.Z[i][self.n-2] * self.mJ[i] - self.mJpr[i]) / (self.Z[i][self.n-2] * self.mH[i] - self.mHpr[i])
            self.Mn[i] = self.Mn[i].real - self.Mn[i].imag * 1j
            self.Nn[i] = (self.Y[i][self.n-2] * self.mJ[i] - self.mJpr[i]) / (self.Y[i][self.n-2] * self.mH[i] - self.mHpr[i])
            self.Nn[i] = self.Nn[i].real - self.Nn[i].imag * 1j
    
    def calculate_angles(self):
        """Расчет углов для диаграммы направленности"""
        for i in range(self.steps):
            if i == 0:
                self.teta[i] = self.teta_start * (math.pi/180)
            else:
                self.teta[i] = self.teta[i-1] + self.step
            self.cos_teta[i] = math.cos(self.teta[i])
    
    def calculate_legendre_functions(self):
        """Расчет функций Лежандра"""
        self.pii = np.zeros((self.toch+1, 2*self.steps+1))
        self.tay = np.zeros((self.toch+1, 2*self.steps+1))
        
        for i in range(self.toch):
            m = i + 1
            Lm0 = scipy.special.lpmv(0, m, self.cos_teta)
            Lm1 = scipy.special.lpmv(1, m, self.cos_teta)
            
            if m < 2:
                Lm2 = 0
            else:
                Lm2 = scipy.special.lpmv(2, m, self.cos_teta)
            
            for z in range(len(self.teta)):
                if (self.teta[z] > 0) and (self.teta[z] < math.pi):
                    self.pii[i][z] = (1 * Lm1[z]) / math.sin(self.teta[z])
                elif (self.teta[z] > math.pi) and (self.teta[z] < 2*math.pi):
                    self.pii[i][z] = (-1 * Lm1[z]) / math.sin(self.teta[z])
            
            for z in range(len(self.teta)):
                if m < 2:
                    self.tay[i][z] = 0.5 * (-m * (m + 1) * Lm0[z])
                else:
                    self.tay[i][z] = 0.5 * (Lm2[z] - m * (m + 1) * Lm0[z])
    
    def calculate_circular_polarization(self):
        """Расчет полей для круговой поляризации"""
        for p in range(1, self.toch):
            for z in range(len(self.teta)):
                self.E_op[p][z] = ((2*p + 1) / (p*(p + 1))) * ((-1)**p) * (self.tay[p][z] - self.pii[p][z]) * (self.Mn[p] + self.Nn[p])
                self.E_kp[p][z] = ((2*p + 1) / (p*(p + 1))) * ((-1)**p) * (self.tay[p][z] + self.pii[p][z]) * (self.Mn[p] - self.Nn[p])
        
        self.P1 = np.sum(self.E_op, axis=0)
        self.P2 = np.sum(self.E_kp, axis=0)
        self.Pab1 = np.abs(self.P1)
        self.Pab2 = np.abs(self.P2)
    
    def calculate_linear_polarization(self):
        """Расчет полей для линейной поляризации"""
        for z in range(len(self.teta)):
            for p in range(self.toch):
                y = p + 1
                self.S_teta[p][z] = ((2*y + 1)/(y*(y + 1))) * ((-1)**y) * (-1 * (self.tay[p][z] * self.Mn[p] - self.pii[p][z] * self.Nn[p]) * 
                                     math.cos(self.teta[z]) * math.cos(self.phi)**2 - (self.pii[p][z] * self.Mn[p] - self.tay[p][z] * self.Nn[p]) * 
                                     math.sin(self.phi)**2)
                self.S_phi[p][z] = ((2*y + 1)/(y*(y + 1))) * ((-1)**y) * ((self.tay[p][z] * self.Mn[p] - self.pii[p][z] * self.Nn[p]) * 
                                    math.cos(self.phi) * math.sin(self.phi) - (self.pii[p][z] * self.Mn[p] - self.tay[p][z] * self.Nn[p]) * 
                                    math.cos(self.teta[z])**2 * math.sin(self.phi)*math.cos(self.phi))
                self.E_teta[0][z] += self.S_teta[p][z]
                self.E_phi[0][z] += self.S_phi[p][z]
            
            self.E_teta[0][z] = (1 - (math.sin(self.teta[z]) * math.cos(self.phi))**2)**(-0.5) * self.E_teta[0][z]
            self.E_phi[0][z] = (1 - (math.sin(self.teta[z]) * math.cos(self.phi))**2)**(-0.5) * self.E_phi[0][z]
            
            for p in range(self.toch):
                self.E_teta[0][z] = abs(self.E_teta[0][z])
                self.E_phi[0][z] = abs(self.E_phi[0][z])
    
    def normalize_results(self):
        """Нормализация результатов"""
        self.tetay = np.zeros(self.steps)
        
        for i in range(self.steps):
            self.teta[i] = self.teta[i] - math.pi
        
        for i in range(len(self.teta)):
            self.tetay[i] = (self.teta[i] * (self.steps / (2 * math.pi)))
        
        # Нормированные E к k0a
        self.DN_NORM_lin_k0a_teta = self.E_teta[0] / self.k0
        self.DN_NORM_lin_k0a_phi = self.E_phi[0] / self.k0
        self.DN_NORM_circle_k0a_op = self.Pab1 / self.k0
        self.DN_NORM_circle_k0a_kp = self.Pab2 / self.k0
        
        # Нормированные E dB
        E_teta_max = np.max(self.E_teta[0])
        E_phi_max = np.max(self.E_phi[0])
        Pab1_max = np.max(self.Pab1)
        Pab2_max = np.max(self.Pab2)
        
        self.DN_NORM_lin_dB_teta = 20 * np.log10(self.E_teta[0] / E_teta_max)
        self.DN_NORM_lin_dB_phi = 20 * np.log10(self.E_phi[0] / E_phi_max)
        self.DN_NORM_circle_dB_op = 20 * np.log10(self.Pab1 / Pab1_max)
        self.DN_NORM_circle_dB_kp = 20 * np.log10(self.Pab2 / Pab2_max)
    
    def plot_results(self):
        """Построение графиков результатов"""
        # График в декартовых координатах
        fig1, ax1 = plt.subplots(figsize=(6, 6))
        ax1.set(xlim=(-180, 180))
        
        ax1.plot(self.tetay, self.DN_NORM_lin_dB_phi, color='red', linestyle='-', linewidth=1, label=r'$E_{\phi}$, dB')
        ax1.plot(self.tetay, self.DN_NORM_lin_dB_teta, color='green', linestyle='-', linewidth=1, label=r'$E_{\theta}$, dB')
        ax1.set_ylim(-60, 0)
        ax1.xaxis.set_major_locator(ticker.MaxNLocator(5))
        
        def format_degrees(x, pos):
            return f'{int(x)}°'
        
        ax1.xaxis.set_major_formatter(ticker.FuncFormatter(format_degrees))
        ax1.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
        ax1.set_ylabel(r'$E_{norm}$, dB', fontsize=14)
        ax1.set_xlabel(r'$\theta$$^{\circ}$', fontsize=14)
        ax1.legend()
        ax1.grid(True)
        
        # Графики в полярных координатах
        fig2, axs = plt.subplots(1, 2, figsize=(10, 5), subplot_kw={'projection': 'polar'})
        
        axs[0].plot(self.teta, self.DN_NORM_lin_dB_phi, color='black', linestyle='-', linewidth=1, label='GreenTensor')
        axs[0].set_title(r'$E_{\phi}$, dB', pad=30)
        axs[0].legend(loc='best')
        axs[0].set_ylim(-60, 0)
        
        axs[1].plot(self.teta, self.DN_NORM_lin_dB_teta, color='blue', linestyle='-', linewidth=1, label='GreenTensor')
        axs[1].set_title(r'$E_{\theta}$, dB', pad=30)
        axs[1].legend(loc='best')
        axs[1].set_ylim(-60, 0)
        
        plt.tight_layout()
        plt.show()
    
    def run_calculation(self):
        """Выполнение полного расчета"""
        #print('Параметры материалов:')
        #print(f'\na = {self.a}\neps = {self.eps}\nmiy = {self.miy}')
        
        
        self.calculate_medium_parameters()
        self.calculate_k_coefficients()
        self.calculate_bessel_functions()
        self.calculate_CS_functions()
        
        # Добавляем последний элемент в alfa и eps если нужно
        if self.eps[len(self.eps)-1] != (len(self.eps)-1):
            self.alfa = np.append(self.alfa, 0)
            self.eps = np.append(self.eps, len(self.eps))
        
        self.calculate_impedances()
        self.calculate_modified_functions()
        self.calculate_scattering_coefficients()
        self.calculate_angles()
        self.calculate_legendre_functions()
        self.calculate_circular_polarization()
        self.calculate_linear_polarization()
        self.normalize_results()
        #self.plot_results()



class ScatteringCalculator:
    def __init__(self, toch=10, k0a_start=0.25, k0a_stop=5.0, k0a_step=0.05):
        """
        Параметры:
            toch: количество членов в разложении
            k0a_start: начальное значение k0a
            k0a_stop: конечное значение k0a
            k0a_step: шаг изменения k0a
        """
        self.toch = toch
        self.k0a_start = k0a_start
        self.k0a_stop = k0a_stop
        self.k0a_step = k0a_step
        
        # Инициализация массивов коэффициентов
        self.i_arr = np.arange(1, self.toch + 1)  # начинаем с 1 до toch включительно
        self.coeffs_1 = 2 * self.i_arr + 1       # (2n + 1)
        self.coeffs_2 = 2 * (self.i_arr**2) - 1  # (2n² - 1)
        self.coeffs_3 = self.i_arr * (self.i_arr + 1)  # n(n+1)
        self.minus_1_pow = (-1)**self.i_arr
        
        # Создаем экземпляр RCSCalculator
        self.rcs_calculator = RCSCalculator()
    
    def calculate(self):
        """Основной метод расчета"""
        # Генерация массива k0a
        self.k0a = np.arange(self.k0a_start, self.k0a_stop + self.k0a_step/2, self.k0a_step)
        k0a_steps = len(self.k0a)
        
        # Инициализация результирующих массивов
        self.sigma_s = np.zeros(k0a_steps)
        self.sigma_theta = np.zeros(k0a_steps)  
        self.sigma_phi = np.zeros(k0a_steps)
        self.sigma_r = np.zeros(k0a_steps)
        self.sigma_p = np.zeros(k0a_steps)
        
        for k, current_k0a in enumerate(self.k0a):
            self._calculate_for_k(k, current_k0a)
        
        return self.k0a, self.sigma_s, self.sigma_r, self.sigma_p, self.sigma_theta, self.sigma_phi
    
    def _calculate_for_k(self, k, current_k0a):
        """Вычисления для конкретного значения k0a"""
        self.rcs_calculator.k0 = current_k0a
        self.rcs_calculator.run_calculation()
        
        Mn = self.rcs_calculator.Mn
        Nn = self.rcs_calculator.Nn
        
        # Векторизованные вычисления
        abs_Mn_sq = np.abs(Mn)**2
        abs_Nn_sq = np.abs(Nn)**2
        
        # Вычисление сумм
        sum_sigma_1 = np.sum(self.coeffs_1 * (abs_Mn_sq + abs_Nn_sq))
        sum_sigma_2 = np.sum(self.coeffs_1 * self.minus_1_pow * (Mn - Nn))
        sum_sigma_3 = np.sum(self.coeffs_1 * (Mn + Nn))
        sum_sigma_4 = np.sum((self.coeffs_1/self.coeffs_3)*(self.coeffs_2 * abs_Mn_sq + self.coeffs_1 * abs_Nn_sq))
        sum_sigma_5 = np.sum((self.coeffs_1/self.coeffs_3) * (self.coeffs_2 * abs_Nn_sq + self.coeffs_1 * abs_Mn_sq))
        
        # Сохранение результатов
        k0a_sq = current_k0a**2
        self.sigma_s[k] = (2 / k0a_sq) * sum_sigma_1
        self.sigma_r[k] = (1 / k0a_sq) * np.abs(sum_sigma_2)**2
        self.sigma_p[k] = (1 / k0a_sq) * np.abs(sum_sigma_3)**2
        self.sigma_theta[k] = (1 / k0a_sq) * sum_sigma_4
        self.sigma_phi[k] = (1 / k0a_sq) * sum_sigma_5
    #  r'$E_{\theta}$'
    def plot_results(self):
        """Визуализация результатов в трех отдельных графиках"""
        plt.figure(figsize=(12, 10))
    
        # Первый график - σ_s и компоненты
        plt.subplot(3, 1, 1)
        plt.plot(self.k0a, self.sigma_s, 'b-', label= r'${\sigma}_{s}$')
        plt.plot(self.k0a, self.sigma_theta, 'r--', label= r'${\sigma}_{\theta}$')
        plt.plot(self.k0a, self.sigma_phi, 'g--', label= r'${\sigma}_{\phi}$')
        plt.title('Полный коэффициент рассеяния')
        plt.xlabel('k0a')
        plt.ylabel(r'${\sigma}$')
        plt.grid(True)
        plt.legend()
        
        # Второй график - σ_r
        plt.subplot(3, 1, 2)
        plt.plot(self.k0a, self.sigma_r, 'm-', label= r'${\sigma}_{r}$')
        plt.title('Радиолокационный коэффициент рассеяния')
        plt.xlabel('k0a')
        plt.grid(True)
        plt.legend()
        
        # Третий график - σ_p
        plt.subplot(3, 1, 3)
        plt.plot(self.k0a, self.sigma_p, 'c-', label= r'${\sigma}_{p}$')
        plt.title('Коэффициент рассеяния в попутном направлении')
        plt.xlabel('k0a')
        plt.grid(True)
        plt.legend()
        
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    # Создание экземпляров
    scattering_calc = ScatteringCalculator()
    rcs_calc = RCSCalculator()
    
    # Выполнение расчетов
    results = scattering_calc.calculate()
    
    # Визуализация результатов
    scattering_calc.plot_results()
    
    # Выполнение расчетов для RCSCalculator
    rcs_calc.run_calculation()
    if hasattr(rcs_calc, 'plot_results'):
        rcs_calc.plot_results()