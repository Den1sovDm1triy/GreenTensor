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


DN_Ansys_FEM_Ephi_0 = [-32.536693327507, -33.3557083621987, -33.1376082787953, -29.7603624673751, -26.6371248065763, -24.9087488883262, -24.5107526750701, -25.3384492900005, -27.2722905942612, -29.7123668384369, -30.8667947899504, -30.2937111475207, -29.725725553496, -29.6515270590993, -29.5285251773541, -28.9906631817453, -28.5384632062469, -28.8644918682823, -30.3254507853649, -32.3863457148081, -32.2384994522856, -30.0078953533655, -28.640643389518, -28.9660059759851, -31.2816072744557, -34.5939210635546, -32.3274628314807, -28.7308292474648, -27.1678900808789, -27.571316157101, -30.1902706433169, -36.3834856418195, -39.8086691355848, -33.2257785243471, -30.9299072980424, -31.3693870527476, -34.4799278284724, -42.3199830841176, -43.6245681940575, -36.9517370775809, -35.2443211299139, -36.3883351323681, -40.9371933430205, -60.0382283685429, -42.5591903606359, -37.2422352237608, -35.2671664273485, -34.9717813828178, -35.880448859093, -37.8829336218933, -41.2505960311666, -47.1680076670882, -52.2433080360066, -44.5903460144919, -40.1271833788271, -37.2630298456857, -35.1426942513487, -33.6778723797238, -33.1016852987786, -33.7653845207159, -36.2220844610054, -41.2948795996646, -42.2598574431503, -37.3083022485879, -35.2298436563153, -35.6874154060343, -38.6212627510913, -42.0015398970117, -39.0321079361072, -36.5603859074145, -36.7672582636487, -39.0428448066898, -38.0483968839596, -33.821234540194, -31.3592292993931, -30.9438477465295, -32.5078118849291, -34.9877308644983, -33.0683091248234, -29.4195440738327, -27.4190385525748, -27.1073345846672, -28.3301783535296, -30.4332709830804, -30.4769456827038, -28.0634751128271, -26.2444949108465, -25.7260984976396, -26.357592832859, -27.5645979273204, -28.0448173186436, -27.2416397383076, -26.3307243245955, -26.0715386022051, -26.3805731206749, -26.6620368109919, -26.3976778115279, -25.8915937935525, -25.7585990730026, -26.2293721646602, -27.0137092158546, -27.3825758438416, -27.047063390316, -26.6861160000483, -26.9887324031148, -28.2607063307316, -30.5211849354094, -33.1987579622046, -34.7087799343592, -34.7906751555668, -34.8660965374188, -35.159024065851, -35.1568147866359, -34.9445603361472, -34.8433171785953, -34.1900660938815, -32.263487194219, -30.0426271339498, -28.5129198828322, -27.9001417903191, -28.0116530915424, -28.2535305020965, -27.8714938643235, -26.932616302465, -26.175948645531, -26.0664944522182, -26.6736397864316, -27.7216695195587, -28.5594139779201, -28.7278013770236, -28.6916749730425, -29.1761618904487, -30.5785134167765, -33.0431305996312, -36.4744273315364, -40.1402094071796, -41.0650694981128, -37.1897086743634, -32.9240774711391, -29.9248570047638, -28.23198044035, -27.7502331035378, -28.4139764748301, -30.0160002537845, -31.4221660328443, -30.6344608894293, -28.6279632709019, -27.1227753898914, -26.5417763176091, -26.9803647684602, -28.6313271121107, -31.895660326848, -36.3312280421973, -35.3236045492401, -32.2372332645883, -31.2291147911532, -31.6782040954471, -31.2018783305435, -28.6613818260013, -26.585443332349, -26.2317149502695, -28.2268386387003, -32.3025754387445, -28.036557837025, -22.8044119621531, -20.1547566700109, -19.4968756412083, -20.7714986207283, -23.9815644305925, -24.5450186844149, -20.5219204524289, -18.347424667206, -19.0635246412245, -24.8362062730051, -20.6352493465365, -11.7587688324725, -6.76313793417833, -3.5669573601608, -1.52466310023244, -0.373690399224927, -0.370897219272283, -1.51931083162391, -3.55962653219378, -6.75517737217892, -11.7541229918702, -20.6517834294917, -24.7614476998208, -18.991556651139, -18.2617926245523, -20.3971983963322, -24.3908760256795, -24.036015517592, -20.8808296635539, -19.6117781522875, -20.2837675253022, -22.9721364007499, -28.2941437907464, -32.2729387149467, -28.0146833951606, -26.0085686068769, -26.3024386063291, -28.2634588480497, -30.7691952844115, -31.5413968860885, -31.3120735118915, -32.3945576999533, -35.5105105240296, -36.4465684901424, -32.0022715557632, -28.7750177293323, -27.1607636101731, -26.7446942645244, -27.3120449848989, -28.7324157735611, -30.5473755040877, -31.1770083163434, -29.8267907139686, -28.2781423817459, -27.6133468656883, -28.0718239482735, -29.7392125918531, -32.7022289357605, -36.7859964516517, -40.0163853842105, -39.2184349727605, -36.1880758334363, -33.0113098322621, -30.60172808598, -29.21300609793, -28.750342647801, -28.8217951274182, -28.6786205623012, -27.8400088643578, -26.7949125412592, -26.2191237260691, -26.3957846119702, -27.2501367642522, -28.250145285899, -28.521535966933, -28.0503711339306, -27.7454470908645, -28.2229146065467, -29.6553334268218, -31.8222236367382, -33.7701141018458, -34.4667121854475, -34.5851310548162, -34.8843453261012, -35.0909325865206, -34.9940760051889, -35.0458306149428, -35.0755186175432, -33.5774702962234, -30.8106143186843, -28.4886593901626, -27.1991393676818, -26.9072963284956, -27.2741277083539, -27.5416616039379, -27.0250905565886, -26.1248297342684, -25.6051821468037, -25.7311594208735, -26.2470705816876, -26.5075627603161, -26.1929003626106, -25.8448875094738, -26.0848753179913, -27.0311310366496, -27.9769001666084, -27.6469819248857, -26.4583882188084, -25.7938036252523, -26.2882482251401, -28.116089212482, -30.5901294034141, -30.5424599378506, -28.3670763794861, -27.1133469596444, -27.428515553215, -29.4749328003059, -33.2772534689269, -35.3145236091839, -32.6445909951993, -31.0192051930736, -31.4509626052464, -34.0406755144892, -38.9033491426652, -40.5064745928916, -37.5892523357563, -37.015895687776, -38.8530645101168, -40.2556567283522, -37.4381441614417, -35.0502312147411, -34.6753999140833, -36.5010295487156, -41.0266449116638, -43.0011531402978, -37.9425527780383, -35.1063533910429, -34.1787932079242, -34.5479458526585, -35.9007687610108, -38.177742186324, -41.8161967804636, -48.5666772850863, -50.8793699011802, -43.053322847797, -39.0402481011453, -36.6545794849097, -35.1711759070258, -34.4847346816824, -34.8243220722156, -36.7587603730571, -42.0170179400654, -62.8244989408897, -40.3829586033005, -35.7864877492512, -34.6193723850613, -36.2033848120384, -41.5273317265042, -40.3222740737137, -33.8680307036094, -31.016326736398, -30.6622422445709, -32.8576151461031, -38.3175891580369, -35.9502518071274, -30.4115263401492, -27.9619097887265, -27.6632425219038, -29.2794943190768, -32.6091995809926, -33.860686000125, -30.8776904958584, -28.9488988085832, -28.8506565946164, -30.3544287351328, -32.4418882813339, -32.0499696706473, -29.8624114093188, -28.4772253997449, -28.2743504054858, -28.9464011568029, -29.8201796826583, -30.162786016042, -30.111875375281, -30.2058760244535, -30.0875413991666, -28.7746177948354, -26.7336950362385, -25.1511464206057, -24.5509512108983, -25.107821090218, -26.9525707829424, -30.0813828784782, -33.0592059362057, -33.0451488144363, -32.536693327507]
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
                    self.Z[i][h] = sqrt_part * (term1 / term2)  /2 
                
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
        ax.set_ylim(-80, 1)
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
    k0=31.4
    toch=60
    n=7
    phi2=1*math.pi/2
    phi1=0*math.pi/2
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
    plotter = ResultPlotter(labl, theta, E2, E1)
    plotter.plot_all()
    