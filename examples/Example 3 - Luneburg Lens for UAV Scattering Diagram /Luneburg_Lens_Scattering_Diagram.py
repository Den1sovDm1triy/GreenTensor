import math
import cmath
import scipy
import matplotlib.pyplot as plt
import json
import matplotlib.ticker as ticker
import numpy as np

#Задаем радиус Линзы
k0 = 4.33 * (math.pi)
print ('Радиус сферы k0 =', k0)

#Точность расчетов берется в два раза больше чем значение k0
toch = 28

print ('Точность расчетов toch =', toch)



#Параметры материала Линзы
#Для этого примера не допускаются параметры с реактивной (внимной) частью параметров следы eps / miy

n = 5
eps = [0 * n for i in range(n)]
miy = [0 * n for i in range(n)]
a = [0 * n for i in range(n)]


a = [0.47, 0.67, 0.82, 0.94, 1]
eps = [1.88, 1.67, 1.44, 1.22, 1]
miy = [1, 1, 1, 1, 1]


    #Параметры материалов
print ('Параметры материалов')
print ('\n a =', a, '\n eps =', eps,'\n miy =', miy)

#Диапазон расчетных углов
teta_start = 0.01
teta_stop = 360
step = math.pi/180
teta_diap = abs(teta_stop)-abs(teta_start)
steps = int(((teta_diap)*(math.pi/180)) / step)
teta = [0 * n for i in range(steps)]
cos_teta = [0 * n for i in range(steps)]

  #Создаем массив для переменных среды
alfa = [0 * n for i in range(n)]
beta = [0 * n for i in range(n)]
etta = [0 * n for i in range(n)]
k= [[0] * n for i in range(n)]

    #Расчет переменных, входящих в коэффициенты k
for i in range(n):
  alfa[i] = math.atan((eps[i]).imag / (eps[i]).real)
  beta[i] = math.atan((miy[i]).imag / (miy[i]).real)
  etta[i] = math.sqrt(math.fabs(eps[i]) * math.fabs(miy[i]))

    #Conrol Point
    #print ('alfa:', alfa)
    #print ('beta:', beta)
    #print ('etta:', etta)

    #Расчет коэффициентов среды k

j = 0; #Индекс для рассчета коэффициентов k

for i in range (n):
  k[i][j] = k0 * a[i] * etta[j]
  if j < n - 1:
    j = j + 1
    k[i][j] = k0 * a[i] * etta[j]

    #Conrol Point
    #print ('k:', k)

    ####################################
    # Определяем переменные для функции Бесселя, Неймана, их производных, C, S, их производных
    ####################################

J = [0 * n for i in range(toch)]
Jpr = [0 * n for i in range(toch)]
N = [0 * n for i in range(toch)]
Npr = [0 * n for i in range(toch)]
C = [[0] * (len(etta)-1) for i in range(toch)]
Cpr = [[0] * (len(etta)-1) for i in range(toch)]
S = [[0] * (len(etta)-1) for i in range(toch)]
Spr = [[0] * (len(etta)-1) for i in range(toch)]

    ####################################
    # Определяем модифицированную функцию Бесселя первого рода как функцию Jfunc(i, j1, j2)
    # где: i -- порядок, j1&j2 -- координаты коэффициента k в массиве k[j1][j2]
    ####################################

def Jfunc(i, j1, j2):
  nu = i + 1
  J = (scipy.special.jv(nu + 0.5, k[j1][j2])) * (math.sqrt(k[j1][j2] * math.pi/2))
  return J

    ####################################
    # Определяем производную функции Бесселя первого рода как функцию Jprfunc(i, j1, j2, tie)
    # где: i -- порядок, j1&j2 -- координаты коэффициента k в массиве k[j1][j2], tie (True | False) -- наличие связи между слоями
    ####################################

def Jprfunc(i, j1, j2, tie):

  nu = i + 1

  if tie == False:
    Jpr = ((nu / (2 * nu + 1)) *  (scipy.special.jv(nu - 0.5, k[j1][j2]) * math.sqrt(k[j1][j2] * math.pi/2)) - \
    ((nu + 1) / (2 * nu + 1)) *  (scipy.special.jv(nu + 1.5, k[j1][j2]) * math.sqrt(k[j1][j2] * math.pi/2)) + \
    (J[i] / k[j1][j2]))
  else:
    Jpr = ((nu / (2 * nu + 1)) * ((scipy.special.jv(nu - 0.5,k[j1][j2]) * (math.sqrt(k[j1][j2] * math.pi/2))) / k[j1][j2])) * k[j1][j2] - \
     (((nu + 1) / (2 * nu + 1)) * ((scipy.special.jv(nu + 1.5,k[j1][j2])) * (math.sqrt(k[j1][j2] * math.pi/2))) / k[j1][j2]) * k[j1][j2] + \
     ((scipy.special.jv(nu + 0.5,k[j1][j2])) * (math.sqrt(k[j1][j2] * math.pi/2))) / k[j1][j2]
  return Jpr

    ####################################
    # Определеям модифицированную фунцию Бесселя второго рода (функцию Неймана) как функцию Nfunc(i, j1, j2)
    # где: i -- порядок, j1&j2 -- координаты коэффициента k в массиве k[j1][j2]
    ####################################

def Nfunc(i, j1, j2):
  N = scipy.special.yv((i+1) + 0.5, k[j1][j2]) * math.sqrt(k[j1][j2]* math.pi/2)
  return N

    ####################################
    # Производная функции Неймана как функцию Nprfunc(i, j1, j2, tie)
    # где: i -- порядок, j1&j2 -- координаты коэффициента k в массиве k[j1][j2], tie (True | False) -- наличие связи между слоями
    ####################################

def Nprfunc(i, j1, j2, tie):

  nu = i + 1

  if tie == False:
    Npr = ((nu / (2 * nu + 1)) *  (scipy.special.yv(nu - 0.5, k[j1][j2]) * math.sqrt(k[j1][j2] * math.pi/2)) - \
    ((nu + 1) / (2 * nu + 1)) *  (scipy.special.yv(nu + 1.5, k[j1][j2]) * math.sqrt(k[j1][j2] * math.pi/2)) + \
    (Nfunc(i, j1, j2) / k[j1][j2]))
  else:
    Npr = (((nu/(2 * nu + 1)) * (((scipy.special.yv(nu - 0.5,k[j1][j2])) * (math.sqrt(k[j1][j2] * math.pi/2))) / k[j1][j2])) * k[j1][j2] - \
    (((nu + 1) / (2 * nu + 1)) * ((scipy.special.yv(nu + 1.5,k[j1][j2])) * (math.sqrt(k[j1][j2] * math.pi/2))) / k[j1][j2]) * k[j1][j2] + \
    (Nfunc(i, j1, j2)) / k[j1][j2])
  return Npr

    ####################################
    # Вычисляем массивы значений функции Бесселя первого рода, функции Неймана и их производных для i = 1...toch (значение точности расчетов)
    # связь между слоями не учитываем (tie == False), для k[0][0]  (j1, j2 == 0)
    ####################################

for i in range(toch):
  J[i] = Jfunc(i, 0, 0)
  Jpr[i] = Jprfunc(i, 0, 0, False)
  N[i] = Nfunc(i, 0, 0)
  Npr[i] = Nprfunc(i, 0, 0, False)


    ####################################
    # Вычисляем массивы значений функций C, S и их производных (Cpr и Spr)
    # в Npr и Jpr учитываем связь между слоями (tie == True)
    # C = J * Npr - N * Jpr
    # Cpr(n) = Jpr * Npr - Npr * Jpr
    # S(n) = N * J - J * N
    # Spr(n) = Npr * J - Jpr * N
    ####################################

for i in range(toch-1):
  for j in range(len(etta)-1):
    C[i][j] = (Jfunc(i, (j+1), (j+1)) * Nprfunc(i, j, (j+1), True)) - (Nfunc(i, (j+1), (j+1)) * Jprfunc(i, j, (j+1), True))
    Cpr[i][j] = Jprfunc(i, (j+1), (j+1), True) * Nprfunc(i, j, (j+1), True) - Nprfunc(i, (j+1), (j+1), True) * Jprfunc(i, j, (j+1), True)
    S[i][j] = Nfunc(i,(j+1),(j+1)) * Jfunc(i,j,(j+1)) - Jfunc(i,(j+1),(j+1)) * Nfunc(i,j,(j+1))
    Spr[i][j] = Nprfunc(i,(j+1),(j+1), True) * Jfunc(i, j, (j+1)) - Jprfunc(i,(j+1),(j+1), True) * Nfunc(i, j, (j+1))

    ####################################

    #Conrol Point
    #print('N:', N)
    #print('Npr:', Npr)
    #print('J:', J)
    #print('Jpr:', Jpr)
    #print('C:', C)
    #print('Cpr:', Cpr)
    #print('S:', S)
    #print('Spr:', Spr)

    ####################################
    # Добавляем в конец массива alfa элемент со значением 0
    # Добавляем в конец массива eps элемент со значением длинны массива eps
    ####################################

if eps[len(eps)-1] != (len(eps)-1):
    alfa.append(0)
    eps.append(len(eps))

    #Conrol Point
    #print('alfa[', len(alfa)-1,']:', alfa[len(alfa)-1])
    #print('eps[', len(eps)-1,']:', eps[len(eps)-1])

    ####################################
    # Определяем импедансы (Z) и адмитансы (Y)
    ####################################

Z = [[0] * (len(a)) for i in range(toch)]
Y = [[0] * (len(a)) for i in range(toch)]

for i in range(toch - 1):
  for h in range(len(a)):

    ####################################
    #Задаем импедансы (Z) и адмитансы (Y) для первого слоя
    ####################################

    if h == 0:
      Z[i][h] = (cmath.sqrt((cmath.exp(alfa[1] * 1j) * abs(eps[1])) / ((cmath.exp(alfa[0] * 1j) * abs(eps[0]))))) * ((Jpr[i])/(J[i]))
      Y[i][h] = (cmath.sqrt((cmath.exp(alfa[0] * 1j) * abs(eps[0])) / ((cmath.exp(alfa[1] * 1j) * abs(eps[1]))))) * ((Jpr[i])/(J[i]))

    ####################################
    #Задаем импедансы (Z) и адмитансы (Y) для последнего слоя
    ####################################

    else:
      if h == (len(a) - 1):
        Z[i][h] = (cmath.sqrt((cmath.exp(alfa[h+1] * 1j) * abs(eps[h+1])) / ((cmath.exp(alfa[h] * 1j) * abs(eps[h]))))) * \
                  (Cpr[i][h-1] + Z[i][h-1] * Spr[i][h-1]) / (C[i][h-1] + Z[i][h-1] * S[i][h-1]) / 2
        Y[i][h] = (cmath.sqrt((cmath.exp(alfa[h] * 1j) * abs(eps[h])) / ((cmath.exp(alfa[h+1] * 1j) * abs(eps[h+1]))))) * \
                  (Cpr[i][h-1] + Y[i][h-1] * Spr[i][h-1]) / (C[i][h-1] + Y[i][h-1] * S[i][h-1]) * 2

    ####################################
    #Задаем импедансы (Z) и адмитансы (Y) для промежуточных слоёв
    ####################################

      else:
        Z[i][h] = (cmath.sqrt((cmath.exp(alfa[h+1] * 1j) * abs(eps[h+1])) / ((cmath.exp(alfa[h] * 1j) * abs(eps[h]))))) * \
                (Cpr[i][h-1] + Z[i][h-1] * Spr[i][h-1]) / (C[i][h-1] + Z[i][h-1] * S[i][h-1])
        Y[i][h] = (cmath.sqrt((cmath.exp(alfa[h] * 1j) * abs(eps[h])) / ((cmath.exp(alfa[h+1] * 1j) * abs(eps[h+1]))))) * \
                (Cpr[i][h-1] + Y[i][h-1] * Spr[i][h-1]) / (C[i][h-1] + Y[i][h-1] * S[i][h-1])

    #Conrol Point
    #print('Z:', Z)
    #print('Y:', Y)

    ####################################
    # Определяем массивы для mJ, mJpr, mH, mHpr
    ####################################

mJ = [0 * n for i in range(toch)]
mJpr = [0 * n for i in range(toch)]
mH = [0 * n for i in range(toch)]
mHpr = [0 * n for i in range(toch)]

    ####################################
    #Функция Ханкеоя второго рода как функцию Hfunc(i)
    #где: i -- порядок
    ####################################

def Hfunc(i):

  nu = i + 1
  H = (scipy.special.hankel1(nu + 0.5,k0)) * (math.sqrt(k0 * math.pi/2))
  return H

    ####################################
    #Производная функции Ханкеоя второго рода как функцию Hprfunc(i)
    #где: i -- порядок
    ####################################

def Hprfunc(i):
  nu = i + 1
  Hpr = ((nu / (2 * nu + 1)) * (((scipy.special.hankel1(nu - 0.5,k0) * (cmath.sqrt(k0 * math.pi/2))) / k0)) * k0 - \
  (((nu + 1) / (2 * nu + 1)) * ((scipy.special.hankel1(nu + 1.5,k0)) * (cmath.sqrt(k0 * math.pi/2))) / k0) * k0 + \
  ((scipy.special.hankel1(nu + 0.5,k0)) * (cmath.sqrt(k0 * math.pi/2))) / k0)
  return Hpr

    ####################################
    #Заполняем массивы mJ, mJpr, mH, mHpr
    ####################################

k1 = k0
k00 = k[0][0]
k[0][0] = k0
for i in range(toch):
  mJ[i] = Jfunc(i, 0, 0)
  mJpr[i] = Jprfunc(i, 0, 0, True)
  mH[i] = Hfunc(i)
  mHpr[i] = Hprfunc(i)
k0 = k1
k[0][0] = k00

    #Conrol Point
    #print('mJ:', mJ)
    #print('mJpr:', mJpr)
    #print('mH:', mH)
    #print('mHpr:', mHpr)

    ####################################
    # Определяем массивы для Mn и Nn
    ####################################

Mn = [0 * n for i in range(toch)]
Nn = [0 * n for i in range(toch)]

    ####################################
    # Задаём Mn и Nn
    ####################################

for i in range(toch):
  nu = i + 1
  Mn[i] = (Z[i][h] * mJ[i] - mJpr[i]) / (Z[i][h] * mH[i] - mHpr[i])
  Mn[i] = Mn[i].real - Mn[i].imag * 1j
  Nn[i] = (Y[i][h] * mJ[i] - mJpr[i]) / (Y[i][h] * mH[i] - mHpr[i])
  Nn[i] = Nn[i].real - Nn[i].imag * 1j

    #Conrol Point
    #print('Mn:', Mn)
    #print('Nn:', Nn)



    #Conrol Point
    #print(step, steps, math.pi)

for i in range(steps):
  if i == 0:
    teta[i] = (teta_start)*(math.pi/180)
  else:
    teta[i] = teta[i-1] + step

    #Conrol Point
    #print(teta)

for i in range(steps):
  cos_teta[i] = math.cos(teta[i])

    #Conrol Point
    #print(cos_teta)

M = [0 * n for i in range(steps)]
Lm0=[0 * n for i in range(steps)]
Lm1=[0 * n for i in range(steps)]
Lm2=[0 * n for i in range(steps)]
pii = [[0] * ((2*steps+1)) for i in range(toch+1)]
tay = [[0] * ((2*steps+1)) for i in range(toch+1)]
m=0

for i in range(toch):
  m = m+1
  M = scipy.special.lpmv(0, m, cos_teta)
  Lm0 = M
  M = scipy.special.lpmv(1, m, cos_teta)
  Lm1 = M

  if m<2:
    Lm2 = 0
  else:
    M = scipy.special.lpmv(2, m, cos_teta)
    Lm2 = M

  for z in range(len(teta)):
    if (teta[z] > 0) & (teta[z] < math.pi):
      pii[i][z] = ((1)*Lm1[z])/(math.sin(teta[z]))
    else:
      if (teta[z] > math.pi) & (teta[z] < 2*math.pi):
        pii[i][z] = ((-1)*Lm1[z])/(math.sin(teta[z]))

  for z in range(len(teta)):
    if m<2:
      tay[i][z] = (1/2)*(-m*(m+1)*Lm0[z])
    else:
      tay[i][z] = (1/2)*(Lm2[z]-m*(m+1)*Lm0[z])

E_teta= [[0] * (len(teta)) for i in range(toch)]
P_teta= [[0] * (len(teta)) for i in range(toch)]
E_kp= [[0] * (len(teta)) for i in range(toch)]
E_op= [[0] * (len(teta)) for i in range(toch)]
y=0





for p in range(1, toch): 
    for z in range(len(teta)):
        E_op[p][z] = ((2*p + 1)/(p*(p + 1))) * ((-1)**p) * (tay[p][z] - pii[p][z]) * (Mn[p] + Nn[p])
        E_kp[p][z] = ((2*p + 1)/(p*(p + 1))) * ((-1)**p) * (tay[p][z] + pii[p][z]) * (Mn[p] - Nn[p])

P1 = np.sum(E_op, axis=0)
P2 = np.sum(E_kp, axis=0)

Pab1 = np.abs(P1)
Pab2 = np.abs(P2)



for z in range(len(teta)):
  for p in range(toch):
    y=p+1
    E_teta[p][z]=((((2*y+1)/(y*(y+1)))*((-1)**y))*(tay[p][z]*Mn[p]-pii[p][z]*Nn[p]))

  for p in range(toch):
    y=p+1
    P_teta[0][z]=(P_teta[0][z]+E_teta[p][z])
  for p in range(toch):
    P_teta[0][z]=abs(P_teta[0][z])

    #Conrol Point
    '''
    print('P_teta:', P_teta)
    for z in range((toch)):
        print('E_teta:', E_teta[z][180])
    '''

    #Normalize E diagram
tetay = [0 * n for i in range(steps)]
DN_NORM = [0 * n for i in range(len(teta))]

P_teta_max = 0

for i in range(len(teta)):
  if P_teta[0][i] > P_teta_max:
    P_teta_max = P_teta[0][i]

    #Conrol Point
    #print(P_teta_max)
for i in range(steps):
  teta[i] = teta[i] - math.pi
for i in range(len(teta)):
  tetay[i] = (teta[i]*(steps/(2*math.pi)))

#tetay.reverse()

    #Conrol Point
    #print(tetay)
for i in range(len(teta)):
  DN_NORM[i] = 20*math.log10(P_teta[0][i]/P_teta_max)

    #Conrol Point
    #print(DN_NORM)
fig, ax = plt.subplots(figsize=(6, 6))
ax.set(xlim=(-180, 180))

# Создание графика
plt.plot(tetay, DN_NORM, color='blue', linestyle='-', linewidth=1, label='Greentensor')

# Настройка шагов сетки на оси X
plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(5))  # Ограничение до 5 шагов

# Добавление знака градуса к подписям сетки
def format_degrees(x, pos):
    return f'{int(x)}°'

plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(format_degrees))

# Включение подписей сетки по оси X
plt.tick_params(axis='x', which='both', bottom=True, labelbottom=True)

# Подписи осей
plt.ylabel(r'$E_{norm}$, dB', fontsize=14)  # Ось Y: Enorm в dB
plt.xlabel(r'$\theta$$^{\circ}$', fontsize=14)  # Ось X: Θ в градусах
plt.tight_layout()

# Легенда
plt.legend()

    #Show grid
plt.grid(True)
    #Polar Plot
fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': 'polar'})

ax.plot(teta, DN_NORM, color='blue', linestyle='-', linewidth=1, label='GreenTensor')
#ax.plot(teta, P1, color='red', linestyle='-', linewidth=1, label='Ansys1')
# Подписи осей
plt.ylabel(r'$E_{norm}$, dB', fontsize=14, labelpad=25)  # Ось Y: Enorm в dB
plt.xlabel(r'$\theta$$^{\circ}$', fontsize=14)  # Ось X: Θ в градусах
plt.tight_layout()

    #Legend
ax.legend(loc='upper right')
#ax.set_yticklabels([])
#ax.set_xticklabels([])
  #Titleplt.title('Polar Scatterplot')
  #Print
print('E0_Norm = ' + str(DN_NORM))
plt.show()