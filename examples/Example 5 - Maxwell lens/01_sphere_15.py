# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

import math
import cmath
import scipy
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from scipy.stats import pearsonr
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import Normalize, LogNorm
import time


DN_Ansys_FEM_Ephi_0 = [-23.990185338254, -24.0135923600654, -24.1097027639874, -24.257619493249, -24.4291406463063, -24.5967269728459, -24.7439243824073, -24.8741533685856, -25.0131052270013, -25.2038131152549, -25.4983783539484, -25.9513428151712, -26.6168744888576, -27.5484688568982, -28.7961340491898, -30.3858955731699, -32.2316292612247, -33.8760293753118, -34.3490550672958, -33.3576855426192, -31.8344704314439, -30.4942185872922, -29.5550604985933, -29.053109285659, -28.9908733493641, -29.379185096221, -30.243439949591, -31.5981136936146, -33.310854222579, -34.6546807655953, -34.2633160164529, -32.3398104761707, -30.247598677324, -28.5664153629155, -27.3842561991549, -26.683679921049, -26.4491615540821, -26.6929957828345, -27.4705752775184, -28.9112394698999, -31.3058356669695, -35.4292387041212, -44.5613708860959, -43.0407402681627, -34.8322320159633, -30.9477423209318, -28.6717257390174, -27.3210697967237, -26.6382163953661, -26.5172173083207, -26.9303318741155, -27.9092172175643, -29.550077093066, -32.0201649079446, -35.3533996617561, -37.5508705535069, -35.2355986125763, -32.0576022853971, -29.7490993468099, -28.229758102401, -27.3186929090256, -26.8956375433646, -26.8814812087636, -27.2109890003418, -27.8022814251344, -28.5159548197616, -29.1194372663322, -29.3372158851287, -29.046492964254, -28.3937011786946, -27.6345651651292, -26.953342858472, -26.4302099149692, -26.0730384013674, -25.8452233446344, -25.6821781689163, -25.5077208248937, -25.2586551229644, -24.91141473473, -24.4919310799077, -24.0600161442685, -23.6826066652708, -23.4141851875989, -23.2886699928631, -23.3169265150741, -23.483339884811, -23.7383895683674, -23.9907646744162, -24.1147036025925, -23.9949491884892, -23.597028214994, -22.9930758711761, -22.31085979185, -21.6679719188039, -21.143930155267, -20.7829323896842, -20.6047956334724, -20.6131099755458, -20.7983032917987, -21.1359188608689, -21.5810826829939, -22.0622220140756, -22.4825607776379, -22.7426926046244, -22.7851024626754, -22.6279451183183, -22.3516521750195, -22.0533142320625, -21.8112042612727, -21.6741353694052, -21.6651398335667, -21.7881361769809, -22.0330563295135, -22.3792301214488, -22.7984722815495, -23.2596276978219, -23.7355271852595, -24.2112133835338, -24.6899918409804, -25.193861823911, -25.7580655456269, -26.4228817608734, -27.2260661125423, -28.1968679194956, -29.3495684598977, -30.6715681983283, -32.0970791154573, -33.455839238663, -34.4188846729934, -34.5867784561243, -33.8288516107737, -32.4215532137592, -30.7404644443394, -29.0378554635525, -27.4466758611745, -26.0324626460877, -24.8259728970625, -23.8390428402163, -23.0718070289627, -22.5151983406991, -22.1502299873843, -21.9448288367657, -21.8496696336554, -21.7969658060001, -21.7095588030505, -21.5257467214371, -21.2304148476408, -20.8675712384846, -20.5227127518728, -20.2938251133088, -20.2744661285, -20.5548730994885, -21.2374674625025, -22.4690625737095, -24.5152376053668, -27.9991789734617, -35.2467369131077, -42.1398080812442, -29.942191998991, -25.2933307090527, -22.7271659726898, -21.3306773470247, -20.853059485628, -21.3075584900352, -22.9463824965053, -26.2371444803495, -27.9476722549083, -22.1102652599143, -16.9405797206346, -13.1614692582796, -10.2551807120067, -7.93546344138625, -6.04642050904107, -4.49600829868165, -3.22604331749419, -2.19790071466837, -1.38513754761066, -0.769426582103044, -0.338204180026377, -0.08327569356358, -0.0868600049644348, -0.345326589557174, -0.779988014966816, -1.39897108079458, -2.21474425931448, -3.24549210463763, -4.51743193126563, -6.06881412944438, -7.95716914219601, -10.2733392426086, -13.1709047068484, -16.9318705779056, -22.0756642274177, -28.0787848440817, -26.6248519227815, -23.2711856543618, -21.6053526684193, -21.1526544581756, -21.6534829472752, -23.0993272643919, -25.7658761057742, -30.679771433204, -45.0806098540013, -34.4941836505203, -27.7592909937593, -24.4278449157321, -22.4583800002688, -21.2757533322534, -20.6293509950343, -20.3785166322184, -20.4235879812102, -20.6750757261901, -21.0383157885018, -21.4124580173747, -21.7085342446879, -21.880836501627, -21.9460970848755, -21.9697708652168, -22.0324834193885, -22.2039410609849, -22.5339531143875, -23.0541829716955, -23.782713489844, -24.7273863473767, -25.8865614373729, -27.2466699197544, -28.7753056107566, -30.4065396587983, -32.0118281454678, -33.3556978632508, -34.0945906531001, -33.9687169391645, -33.0783907456512, -31.7867031521275, -30.413127474398, -29.1323175916629, -28.0156140594431, -27.0791405762893, -26.3100440465458, -25.6786461882345, -25.1450963306376, -24.6658319692499, -24.2022403281819, -23.7304045788961, -23.2476354184044, -22.7718097511525, -22.3339817270692, -21.9687583150929, -21.7066538683259, -21.5695823019003, -21.5681863396114, -21.6989602496768, -21.939579143795, -22.242491262565, -22.5309371980233, -22.7082518081583, -22.6914843866101, -22.4570796852574, -22.0573597592036, -21.5894771560877, -21.1506257922315, -20.8138114339364, -20.6252637681911, -20.6104385467241, -20.779694876194, -21.1303452941205, -21.6440623191779, -22.2793209011776, -22.9608252204713, -23.5756723747076, -23.9973527780341, -24.147555439695, -24.050085528991, -23.8131113839328, -23.5618082638653, -23.390513160124, -23.3524626938754, -23.4667059070456, -23.7256470385474, -24.0987093723195, -24.5345114379246, -24.9676753985171, -25.3365574558116, -25.6096011620859, -25.8035859159798, -25.9774000098189, -26.2066108459247, -26.5574929916698, -27.0686834316938, -27.7332146906758, -28.4695588179737, -29.0896038110368, -29.3349590102212, -29.064434348225, -28.410527705866, -27.652460964177, -27.0204304450518, -26.6499801639437, -26.619289132967, -26.9895872411759, -27.8345360183478, -29.2663168737355, -31.45881597195, -34.5517841000342, -37.3230575372772, -35.7338845142099, -32.3352317641627, -29.733033088893, -27.9933653622772, -26.941569891376, -26.4715301288025, -26.5452317091529, -27.1860871933674, -28.4959948665986, -30.7257217096227, -34.5356400853707, -42.4355726482404, -44.3334641771887, -35.2766084866084, -31.1260290415924, -28.7155436891011, -27.2680984361912, -26.4921799236979, -26.2595961466797, -26.5180367531602, -27.2613544429054, -28.5170103802687, -30.3249159772457, -32.6258650543074, -34.7684441139172, -35.0354804016477, -33.3541368609218, -31.4509938357518, -30.0242794033804, -29.1419472622821, -28.760212617375, -28.8415424641008, -29.3700820728178, -30.3425748453381, -31.7285770035141, -33.3301384778654, -34.4603920330316, -34.1356435257382, -32.5656644940854, -30.7500307327966, -29.1816472611157, -27.9552564495533, -27.0437768653467, -26.3936682408238, -25.9470714774113, -25.645256013503, -25.4299776666035, -25.2480216265103, -25.0594235848957, -24.8456419170974, -24.6119582844169, -24.3815931930591, -24.1848242486951, -24.0486580934673, -23.990185338254]
class RCSCalculator:
    def __init__(self, k0=5, toch=20, n=7, phi=0, a=None, eps=None, miy=None, k1=None):
        """
        Инициализация калькулятора RCS
        max k0=100000, toch=93
        Параметры:
        k0 - волновое число
        toch - количество учитываемых членов ряда
        n - число слоев (последний слой - воздух)
        phi - азимутальный угол (по умолчанию pi)
        a, eps, miy - параметры слоев (радиусы, диэлектрические проницаемости, магнитные проницаемости)
        """
        if k1 == None: self.k1=k0 
        else: self.k1 = k1 
        self.k0 = k0
        self.toch = toch
        self.n = n
        self.phi = phi
        
        # Параметры материалов по умолчанию (3 слоя)
        self.a = a if a is not None else [0.15, 0.30, 0.45, 0.63, 0.82, 1, 1]
        self.eps = eps if eps is not None else [3.90, 3.52, 2.98, 2.48, 1.76, 1.50, 1]
        # 3.48 -0.0037j 0.55
        #a = [0.53, 0.75, 0.93, 1] [0.001,  0.05, 1]
        #eps = [1.86, 1.57, 1.28, 1]  [10e20-10e10j, 3.48 -0.0037j, 1]
        self.miy = miy if miy is not None else [1, 1, 1, 1, 1, 1, 1]
        
        # Инициализация всех необходимых массивов
        self._initialize_arrays()
        
    def _initialize_arrays(self):
        """Инициализация всех массивов и переменных"""
        # Углы для расчета
        self.teta_start = 0.01
        self.teta_stop = 360
        self.step = math.pi/180
        self.teta_diap = abs(self.teta_stop) - abs(self.teta_start)
        self.steps = int(((self.teta_diap * (math.pi/180)) / self.step)) + 1
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
            self.alfa[i] = cmath.atan(self.eps[i].imag / self.eps[i].real) if self.eps[i].real != 0 else math.pi/2
            self.beta[i] = math.atan(self.miy[i].imag / self.miy[i].real)
            self.sigma[i] = cmath.sqrt(abs(self.eps[i]) * abs(self.miy[i]))
            if self.eps[-1] != len(self.eps) - 1:
                self.alfa = np.append(self.alfa, 0)
                self.eps = np.append(self.eps, len(self.eps))
    
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
        J = (scipy.special.jv(nu + 0.5, self.k[j1][j2])) * (cmath.sqrt(self.k[j1][j2] * math.pi/2))
        return J
    
    def Jprfunc(self, i, j1, j2, tie):
        """Производная функции Бесселя первого рода"""
        nu = i + 1

        if tie == False:
            Jpr = ((nu / (2 * nu + 1)) *  (scipy.special.jv(nu - 0.5, self.k[j1][j2]) * cmath.sqrt(self.k[j1][j2] * math.pi/2)) - \
            ((nu + 1) / (2 * nu + 1)) *  (scipy.special.jv(nu + 1.5, self.k[j1][j2]) * cmath.sqrt(self.k[j1][j2] * math.pi/2)) + \
            (self.J[i] / self.k[j1][j2]))
        else:
            Jpr = ((nu / (2 * nu + 1)) * ((scipy.special.jv(nu - 0.5,self.k[j1][j2]) * (cmath.sqrt(self.k[j1][j2] * math.pi/2))) / self.k[j1][j2])) * self.k[j1][j2] - \
            (((nu + 1) / (2 * nu + 1)) * ((scipy.special.jv(nu + 1.5,self.k[j1][j2])) * (cmath.sqrt(self.k[j1][j2] * math.pi/2))) / self.k[j1][j2]) * self.k[j1][j2] + \
            ((scipy.special.jv(nu + 0.5,self.k[j1][j2])) * (cmath.sqrt(self.k[j1][j2] * math.pi/2))) / self.k[j1][j2]
        return Jpr
    
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
        for i in range(self.toch):
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
                    self.Z[i][h] = sqrt_part * (term1 / term2) /2
                
                    numerator = cmath.exp(self.alfa[h] * 1j) * abs(self.eps[h])
                    denominator = cmath.exp(self.alfa[h+1] * 1j) * abs(self.eps[h+1])
                    sqrt_part = cmath.sqrt(numerator / denominator)
                    term1 = self.Cpr[i][h-1] + self.Y[i][h-1] * self.Spr[i][h-1]
                    term2 = self.C[i][h-1] + self.Y[i][h-1] * self.S[i][h-1]
                    self.Y[i][h] = sqrt_part * (term1 / term2) *2
                
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
    def Hfunc(self, i, k1):
        """Функция Ханкеля второго рода"""
        nu = i + 1
        H = (scipy.special.hankel1(nu + 0.5,self.k1)) * (cmath.sqrt(self.k1 * math.pi/2))
        return H
    
    def Hprfunc(self, i, k1):
        """Производная функции Ханкеля второго рода"""
        nu = i + 1
        Hpr = ((nu / (2 * nu + 1)) * (((scipy.special.hankel1(nu - 0.5,self.k1) * (cmath.sqrt(self.k1 * math.pi/2))) / self.k1)) * self.k1 - \
        (((nu + 1) / (2 * nu + 1)) * ((scipy.special.hankel1(nu + 1.5,self.k1)) * (cmath.sqrt(self.k1 * math.pi/2))) / self.k1) * self.k1 + \
        ((scipy.special.hankel1(nu + 0.5,self.k1)) * (cmath.sqrt(self.k1 * math.pi/2))) / self.k1)
        return Hpr

    
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
            n = i + 1
            #Элемент Гюгенса дифракция
            self.Mn[i] = (self.Z[i][self.n-1] * self.mJ[i] - self.mJpr[i]) / (self.Z[i][self.n-1] * self.mH[i] - self.mHpr[i])
            self.Nn[i] = (self.Y[i][self.n-1] * self.mJ[i] - self.mJpr[i]) / (self.Y[i][self.n-1] * self.mH[i] - self.mHpr[i])
            #Рупор на поверхности сферы
            #self.Mn[i] = (self.Z[i][self.n-1]  - 1j) / (self.Z[i][self.n-1] * self.mH[i] - self.mHpr[i])             
            #self.Nn[i] = (self.Y[i][self.n-1] - 1j) / (self.Y[i][self.n-1] * self.mH[i] - self.mHpr[i])
            self.Mn[i] = self.Mn[i].real - self.Mn[i].imag * 1j
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
            M = scipy.special.lpmv(0, m, self.cos_teta)
            Lm0 = M
            M = scipy.special.lpmv(1, m, self.cos_teta)
            Lm1 = M
            
            if m < 2:
                Lm2 = 0
            else:
                M = scipy.special.lpmv(2, m, self.cos_teta)
                Lm2 = M
            
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
        for z in range(len(self.teta)):
            for p in range(self.toch):
                y=p+1
                self.E_op[p][z] = ((2*y + 1) / (y*(y + 1))) * ((-1)**y) * (self.tay[p][z] - self.pii[p][z]) * (self.Mn[p] + self.Nn[p])
                self.E_kp[p][z] = ((2*y + 1) / (y*(y + 1))) * ((-1)**y) * (self.tay[p][z] + self.pii[p][z]) * (self.Mn[p] - self.Nn[p])
        
        self.P1 = np.sum(self.E_op, axis=0)
        self.P2 = np.sum(self.E_kp, axis=0)
        self.Pab1 = np.real(self.P1)
        self.Pab2 = np.real(self.P2)
    
    def calculate_linear_polarization(self):
        """Расчет полей для линейной поляризации"""
        for z in range(len(self.teta)):
            for p in range(self.toch):
                y = p + 1
                #self.S_teta[p][z] = ((2*y + 1)/(y*(y + 1))) * ((-1)**y) * (self.tay[p][z] * self.Mn[p] - self.pii[p][z] * self.Nn[p])
                self.S_teta[p][z] = ((2*y + 1)/(y*(y + 1))) * ((-1)**y) * (-1 * (self.tay[p][z] * self.Mn[p] - self.pii[p][z] * self.Nn[p]) * 
                                     math.cos(self.teta[z]) * math.cos(self.phi)**2 - (self.pii[p][z] * self.Mn[p] - self.tay[p][z] * self.Nn[p]) * 
                                     math.sin(self.phi)**2)
                #self.S_phi[p][z] = ((2*y + 1)/(y*(y + 1))) * ((-1)**y) * (self.tay[p][z] * self.Nn[p] - self.pii[p][z] * self.Mn[p])
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
            self.teta[i] = self.teta[i] 
        
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

#####################
# СРАВНЕНИЕ С ANSYS #
#####################

class ResultPlotter:
    def __init__(self, comp, theta, E1, E2):
        """
        Инициализация класса.
        :param theta: Массив углов в градусах (0-360).
        :param E1: Массив значений E1.
        :param E2: Массив значений E2.
        """
        if comp == 1*math.pi/2:
            self.comp = r'$_\phi$'
        else:
            self.comp = r'$_\theta$'
            
        self.theta = np.array(theta)
        self.E1 = (np.array(E1))
        self.E2 = (np.array(E2))
        self.E3 = (np.array(E3))
        if np.min(self.E1+self.E2) > -60: self.Ymin = np.min(self.E1+self.E2)
        else: self.Ymin = -1
        #self.Ymin = np.min(self.E1+self.E2)
        self.Ymax = np.max(self.E1)
        # Сдвигаем массивы на 180 элементов
        #self.E1_shifted = np.roll(self.E1, 180)
        self.E1_shifted = (self.E1)
        #self.E2_shifted = np.roll(self.E2, 180)
        self.E2_shifted = (self.E2)
        # Создаем маску для области в пределах -3 дБ от максимума E1
        E1_max = np.max(self.E1_shifted)
        self.valid_mask = self.E1_shifted >= (E1_max - 3)
        self.squared_errors = np.zeros_like(self.E1_shifted)
        self.squared_errors3 = np.zeros_like(self.E1_shifted)
        self.squared_errors[self.valid_mask] = np.sqrt((self.E1_shifted[self.valid_mask] - self.E2_shifted[self.valid_mask]) ** 2)
        self.squared_errors3 = np.sqrt(abs((self.E1_shifted - self.E2_shifted) ** 2))
        # Среднеквадратичная ошибка только для валидной области
        self.mseMax = np.mean(self.squared_errors3)
        self.mse = np.mean(self.squared_errors[self.valid_mask]) if np.any(self.valid_mask) else 0
        # Вычисляем корреляцию Пирсона только для валидной области
        self.pearson_corr = 0
        self.regression_line = None
        

    def calc_pearson_correlation(self):
        """Вычисление корреляции Пирсона между E1 и E2 в валидной области."""
        if not np.any(self.valid_mask):
            return
        
        # Выбираем данные из валидной области
        E1_valid = np.real(self.E1_shifted[self.valid_mask])
        E2_valid = np.real(self.E2_shifted[self.valid_mask])
        
        # Вычисляем корреляцию Пирсона
        self.pearson_corr, _ = pearsonr(E1_valid, E2_valid)
        
        # Рассчитываем линию регрессии
        coeffs = np.polyfit(E1_valid, E2_valid, 1)
        self.regression_line = np.poly1d(coeffs)
    def plot_polar(self):
        """Построение полярного графика E1 и E2 (после сдвига)."""
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
        theta_rad = (self.theta)
        ax.set_ylim(-50, 1)
        # Графики со сдвинутыми данными
        ax.plot(theta_rad, self.E1_shifted, 
                color='#1f77b4', linestyle='--', linewidth=1.5, 
                alpha=0.9, label=r'E$_\varphi$($\theta$), дБ (Ansys HFSS)')
        ax.plot(theta_rad, self.E2_shifted, 
                color='#ff7f0e', linestyle='-', linewidth=2.5, 
                alpha=0.9, label=r'E$_\varphi$($\theta$), дБ (GreenTensor)')
        #ax.set_title(f'E{self.comp}, дБ', fontsize=15, pad=20)
        #ax.set_title(r'$k_0a = 0.6, \varepsilon = 78.73 - i \cdot 12.28$', fontsize=14, pad=20)
        ax.legend(loc='best', fontsize=14, frameon=True, framealpha=0.95)
        ax.set_theta_zero_location('E')  # 0° сверху
        ax.set_theta_direction(1)       # По часовой стрелке
        plt.legend(loc='upper right')
        plt.tight_layout()
        ax.text(0,  # 90° 
        ax.get_rmax() +12, r'$\theta$$\degree$', fontsize=14, ha='right', va='center')
        return fig

    def plot_cartesian(self):
        """Декартов график E1 и E2 (после сдвига) для углов 0-180 градусов."""
        fig, ax = plt.subplots(figsize=(8, 6))
        # Маска для углов 0-180
        mask = (self.theta >= 0) & (self.theta <= 180)
        theta_subset = np.degrees(self.theta[mask])
        E1_subset = self.E1_shifted[mask]
        E2_subset = self.E2_shifted[mask]
        
        # Построение графиков
        ax.plot(theta_subset, E1_subset, 
                color='#1f77b4', linestyle='--', linewidth=1.5, 
                alpha=0.9, label='Ansys')
        ax.plot(theta_subset, E2_subset, 
                color='#ff7f0e', linestyle='-', linewidth=2.5, 
                alpha=0.9, label='GreenTensor')
        # Выделение валидной области
        angle_mask = (self.theta >= 0) & (self.theta <= 180)

        valid_mask = self.valid_mask[angle_mask]
        ax.fill_between(theta_subset, -70, 5, 
                        where=valid_mask, color='green', 
                        alpha=0.15, label='Область по уровню -3дБ')
        ax.set_xlim(0, 180)
        ax.set_ylim(self.Ymin, 0)
        ax.set_xlabel(r'$\theta$$\degree$', fontsize=14)
        ax.set_ylabel(f'E{self.comp}, дБ', fontsize=14)
        ax.legend(fontsize=12)
        ax.grid(True)
        # Форматирование осей
        ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(15))
        plt.tight_layout()
        plt.legend(loc='upper right')
        return fig

    def plot_error(self):
        fig, ax = plt.subplots(figsize=(8, 6))
        # Маска для углов 0-180
        angle_mask = (self.theta >= 0) & (self.theta <= 180)
        theta_subset = np.degrees(self.theta[angle_mask])
        errors_subset = self.squared_errors[angle_mask]
        
        # Отдельно выделяем валидную область
        valid_mask = self.valid_mask[angle_mask]
        ax.fill_between(theta_subset, np.min(errors_subset), np.max(self.squared_errors3), 
                        where=valid_mask, color='green', 
                        alpha=0.15, label='Область по уровню -3дБ')
        ax.plot(theta_subset, self.squared_errors3, 
                color='#2ca02c', linestyle='-', linewidth=2.5, 
                alpha=0.9, label='Cреднеквадратическая ошибка')
        
        ax.set_xlim(0, 360)
        ax.set_ylim(0)
        ax.set_xlabel(r'$\theta$$\degree$', fontsize=14)
        ax.set_ylabel('Cреднеквадратическая ошибка, дБ', fontsize=14)
        ax.legend(fontsize=12)
        ax.grid(True)
        
        # Форматирование осей
        ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(15))
        plt.tight_layout()
        plt.legend(loc='upper right')
        print(f'RMSE(-3dB): {self.mse:.4f}')
        print(f'RMSE: {self.mseMax:.4f}')
        return fig
    
    def plot_pearson_correlation(self):
        self.calc_pearson_correlation()
        """Визуализация корреляции Пирсона в валидной области."""
        if not np.any(self.valid_mask):
            print("No valid data points for correlation analysis")
            return None
        
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Выбираем данные из валидной области
        E1_valid = self.E1_shifted[self.valid_mask]
        E2_valid = self.E2_shifted[self.valid_mask]
        
        # Создаем scatter plot
        ax.scatter(E1_valid, E2_valid, 
                   color='#1f77b4', alpha=0.7, s=60,
                   label='Точки данных')
        
        # Добавляем линию регрессии
        x_vals = np.linspace(min(E1_valid), max(E1_valid), 100)
        y_vals = self.regression_line(x_vals)
        ax.plot(x_vals, y_vals, 
                color='#d62728', linestyle='--', linewidth=2.5,
                label=f'Линия регрессии: y = {self.regression_line.coeffs[0]:.4f}x + {self.regression_line.coeffs[1]:.4f}')
        
        # Добавляем идеальную линию
        min_val = min(min(E1_valid), min(E2_valid))
        max_val = max(max(E1_valid), max(E2_valid))
        ax.plot([min_val, max_val], [min_val, max_val], 
                color='#2ca02c', linestyle=':', linewidth=2,
                label='Эталонная корреляция (y = x)')
        
        # Настраиваем график
        ax.set_xlabel(f'E{self.comp}, дБ (Ansys)', fontsize=14)
        ax.set_ylabel(f'E{self.comp}, дБ (GreenTensor)', fontsize=14)
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)
        
        # Добавляем аннотацию с коэффициентом корреляции
        ax.annotate(f'r = {self.pearson_corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    fontsize=14, bbox=dict(boxstyle="round,pad=0.3", 
                                           fc="white", ec="gray", alpha=0.8))
        
        # Устанавливаем равные масштабы осей
        ax.set_aspect('equal', adjustable='box')
        plt.legend(loc='upper right')
        plt.tight_layout()
        print(f'Корреляция Пирсона: {self.pearson_corr:.4f}')
        return fig

    def plot_all(self):
        """Построение всех графиков."""
        figs = []
        figs.append(self.plot_polar())
        #figs.append(self.plot_cartesian())
        #figs.append(self.plot_error())
        #figs.append(self.plot_pearson_correlation())
        plt.show()
        return figs


#######################################################
#  #
#######################################################
    
class ScatteringCalculator:
    def __init__(self, tochl=10, k0a_start=0.25, k0a_stop=5.0, k0a_step=0.05):
        """
        Параметры:
            toch: количество членов в разложении
            k0a_start: начальное значение k0a
            k0a_stop: конечное значение k0a
            k0a_step: шаг изменения k0a
        """
        self.toch = tochl
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
        
        abs_Mn_sq = np.abs(Mn)**2
        abs_Nn_sq = np.abs(Nn)**2
        
        # Вычисление сумм
        sum_sigma_1 = np.sum(self.coeffs_1 * (abs_Mn_sq + abs_Nn_sq))
        sum_sigma_2 = np.sum(self.coeffs_1 * self.minus_1_pow * (Mn - Nn))
        sum_sigma_3 = np.sum(self.coeffs_1 * (Mn + Nn))
        sum_sigma_4 = np.sum((self.coeffs_1/self.coeffs_3)*(self.coeffs_2 * abs_Mn_sq + self.coeffs_1 * abs_Nn_sq))
        sum_sigma_5 = np.sum((self.coeffs_1/self.coeffs_3) * (self.coeffs_2 * abs_Nn_sq + self.coeffs_1 * abs_Mn_sq))
        
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


############################
# ПОСТРОЕНИЕ ТЕПЛОВЫХ КАРТ #
############################
        

class ParameterSweep:
    def __init__(self, simulation_class, k0_start, k0_stop, k0_step, phi0):
        """
        Класс для параметрического анализа по волновому числу k0.
        
        :param simulation_class: Класс для моделирования (должен иметь метод run и поле E)
        :param k0_start: Начальное значение k0
        :param k0_stop: Конечное значение k0
        :param k0_step: Шаг изменения k0
        :param phi0: Фиксированное значение угла phi0
        """
        self.simulation_class = simulation_class
        self.k0_start = k0_start
        self.k0_stop = k0_stop
        self.k0_step = k0_step
        self.phi0 = (phi0)
        self.k0_values = None
        self.results = None

    def run_sweep(self):
        """Выполнить цикл по k0 и собрать результаты."""
        n_points = int(np.ceil((self.k0_stop - self.k0_start) / self.k0_step)) + 1
        self.k0_values = np.linspace(self.k0_start, self.k0_stop, n_points)
        self.results = []
        
        for k0 in self.k0_values:
            # Создаем новый экземпляр класса для каждого k0
            sim = self.simulation_class(k1=k0, phi=self.phi0)
            sim.run_calculation()  # Запускаем расчет
            
            # Получаем и проверяем поле E
            e_field = (np.real(sim.DN_NORM_lin_k0a_teta))
            self.results.append(e_field)
            print(self.phi0)
        
        # Преобразуем результаты в массив (k0 x углы)
        self.results = np.array(self.results)

    def plot_heatmap(self, output_file='heatmap.png'):
        """Построить тепловую карту в полярных координатах."""
        if self.results is None:
            raise RuntimeError("Сначала выполните run_sweep()")
        
        # Создаем полярную сетку
        theta = np.radians(np.linspace(0, 359, 360))  # Углы в радианах (361 точка для замкнутой сетки)
        r = self.k0_values  # Радиусы = значения k0
        
        # Создаем 2D сетку для pcolormesh
        R, Theta = np.meshgrid(r, theta)
        
        # Транспонируем результаты для соответствия размерности сетки
        data = self.results.T  # Теперь форма: (углы x k0)
        
        # Создаем график
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 8))
        
        # Построение тепловой карты
        pc = ax.pcolormesh(Theta, R, data, shading='auto', cmap='turbo', vmin=0, vmax=1.4)
        
        # Настройка внешнего вида
        ax.set_theta_zero_location('E')  # 0° наверху (север)
        ax.set_theta_direction(1)       # По часовой стрелке
        #plt.title(f'Тепловая карта (φ={np.degrees(self.phi0):f} град.)')
        fig.colorbar(pc, ax=ax, label=r'E$_\varphi$, нормированное к $k_0a$')
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        #plt.close(fig)

#######################################################
# СРАВНЕНИЕ ДИАГРАММ НАПРАВЛЕННОСТЕЙ С МЕНЯЮЩИМСЯ k0a #
#######################################################
        
class RCSVariation:
    def __init__(self, k0_start, k0_stop, step, toch=10, n=3, phi=0.0, a=None, eps=None, miy=None, k1=None):
        """
        Инициализация параметров для вариационного анализа.
        
        :param k0_start: Начальное значение k0
        :param k0_stop: Конечное значение k0
        :param step: Шаг изменения k0
        :param toch: Размер сетки (по умолчанию 20)
        :param n: Порядок аппроксимации (по умолчанию 4)
        :param phi: Угол падения в радианах (по умолчанию 0.0)
        :param a, eps, miy, k1: Дополнительные параметры для RCSCalculator
        """
        self.k0_start = k0_start
        self.k0_stop = k0_stop
        self.step = step
        self.toch = toch
        self.n = n
        self.phi = phi
        self.a = a
        self.eps = eps
        self.miy = miy
        self.k1 = None
        self.results = []  # Список для хранения результатов 
        self.k0_values = []  # Список использованных значений k0
        
    def run_calculations(self):
        """Выполняет цикл расчетов для диапазона k0."""
        k0_range = np.arange(self.k0_start, self.k0_stop + self.step, self.step)
        
        for k0 in k0_range:
            # Создаем экземпляр калькулятора с текущим k0
            calculator = RCSCalculator(
                k1=k0,
                toch=self.toch,
                n=self.n,
                phi=self.phi,
                a=self.a,
                eps=self.eps,
                miy=self.miy
            )
            # Сохраняем значение k0
            self.k0_values.append(k0)
            calculator.run_calculation()
            # Сохраняем результат расчета 
            self.results.append(np.roll(np.abs(calculator.DN_NORM_lin_k0a_teta), 180))

    def plot_polar(self):
        """Строит полярный график для всех результатов."""
        plt.figure(figsize=(10, 10))
        ax = plt.subplot(111, projection='polar')
        # Углы в градусах (0-359) -> радианы
        theta = np.deg2rad(np.arange(0, 360))
        #plt.ylim(0, 0.4)
        # Рисуем каждый результат
        for i, result in enumerate(self.results):
            ax.plot(theta, result, label=r'$k_0a$'+f'={self.k0_values[i]:.1f}')
        
        # Настройка графика
        ax.set_theta_zero_location('N')  # 0° наверху
        ax.set_theta_direction(-1)      # По часовой стрелке
        ax.set_title(r'F$_\varphi$, $\varphi$ = 0$\degree$', pad=20, fontsize=20)
        ax.set_theta_zero_location('E')
        ax.set_theta_direction(1)
        ax.text(0,  # 90° 
        ax.get_rmax() +.01, r'$\theta$$\degree$', fontsize=20, ha='right', va='center')
        ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.1), fontsize=18)
        
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    # Создание экземпляров
    #scattering_calc = ScatteringCalculator()
    k0=15
    toch=30
    n=7
    phi2=1*math.pi/2
    phi1=1*math.pi/2
    a=None
    eps=None
    miy=None
    k1=k0
    rcs_calc1 = RCSCalculator(k0, toch, n, phi1, a, eps, miy, k1=1*k0)
    rcs_calc2 = RCSCalculator(k0, toch, n, phi2, a, eps, miy, k1=1*k0)
    
    # Выполнение расчетов
    #results = scattering_calc.calculate()
    
    # Визуализация результатов
    #scattering_calc.plot_results()
    # Выполнение расчетов для RCSCalculator

    rcs_calc1.run_calculation()
    rcs_calc2.run_calculation()

    theta = rcs_calc1.teta
    labl = rcs_calc1.phi
    E1 = rcs_calc1.DN_NORM_lin_dB_teta
    E2 = np.roll(DN_Ansys_FEM_Ephi_0, 0)
    E3 = np.roll(DN_Ansys_FEM_Ephi_0, 0)
    # Создание экземпляра класса и построение графиков
    #1 - phi, 0 - theta
    print(E1)
    plotter = ResultPlotter(labl, theta, E2, E1)
    plotter.plot_all()
    