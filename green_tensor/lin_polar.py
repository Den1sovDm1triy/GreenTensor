import math
import cmath
import scipy.special

#Задаем радиус Линзы
k0 = 6 * (math.pi)
print ('Радиус сферы k0 =', k0)

#Точность расчетов берется в два раза больше чем значение k0
#Функция ceil() определяет, какая из границ интервала наибольшая и записывает её в результат округления.
#toch = math.ceil(k0*2)
toch = 50
print ('Точность расчетов toch =', toch)

#Параметры материала Линзы
#Для этого примера не допускаются параметры с реактивной (внимной) частью параметров следы eps / miy
#так как из формулы расчета etta был исключен первый множитель (экспонента)
#возникала ошибка can't convert complex to float (нужно решить чуть позднее)
#a - нормированные радиусы слоев
a = [0.53, 0.75, 0.93, 1]
#a = [0.34, 0.49, 0.59, 0.69, 0.77, 0.84, 0.91, 0.97]
#eps - диэлектрическая проницаемость материала
eps = [1.86, 1.57, 1.28, 1]
#eps = [1, 1, 1, 1, 1, 1, 1, 1]
#miy = [1.94, 1.82, 1.71, 1.59, 1.47, 1.35, 1.24, 1.12]
miy = [1, 1, 1, 1]
#miy - магнитная проницаемость материала
print ('\n a =', a, '\n eps =', eps,'\n miy =', miy)

#Расчет коэффициентов k, которые связывают слои
#Например для 4-х слойной линзы это: k0a1, k1a1, k2a1, k2a2, k3a2, k3a3, k4a3, k4a4

#Создаем массив для переменных среды
n = len(eps)
alfa = [0 * n for i in range(n)]
beta = [0 * n for i in range(n)]
etta = [0 * n for i in range(n)]
k= [[0] * n for i in range(n)]


#Расчет переменных, входящих в коэффициенты k
print('Расчет переменных, входящих в коэффициенты k')

for i in range(n):
  alfa[i] = math.atan((eps[i]).imag / (eps[i]).real)
  beta[i] = math.atan((miy[i]).imag / (miy[i]).real)
  etta[i] = math.sqrt(math.fabs(eps[i]) * math.fabs(miy[i]))


print ('alfa:', alfa)
print ('beta:', beta)
print ('etta:', etta)
print('###-###-###')

#Расчет коэффициентов k
print('Расчет переменных, входящих в коэффициенты k')

j = 0; #Индекс для рассчета коэффициентов k

for i in range (n):
  k[i][j] = k0 * a[i] * etta[j]
  if j < n - 1:
    j = j + 1
    k[i][j] = k0 * a[i] * etta[j]
print ('k:', k)

print('###-###-###')

from traitlets.traitlets import ForwardDeclaredInstance
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
#
# C = J * Npr - N * Jpr
# Cpr(n) = Jpr * Npr - Npr * Jpr
# S(n) = N * J - J * N
# Spr(n) = Npr * J - Jpr * N
#
####################################

for i in range(toch-1):
  for j in range(len(etta)-1):
    C[i][j] = (Jfunc(i, (j+1), (j+1)) * Nprfunc(i, j, (j+1), True)) - (Nfunc(i, (j+1), (j+1)) * Jprfunc(i, j, (j+1), True))
    Cpr[i][j] = Jprfunc(i, (j+1), (j+1), True) * Nprfunc(i, j, (j+1), True) - Nprfunc(i, (j+1), (j+1), True) * Jprfunc(i, j, (j+1), True)
    S[i][j] = Nfunc(i,(j+1),(j+1)) * Jfunc(i,j,(j+1)) - Jfunc(i,(j+1),(j+1)) * Nfunc(i,j,(j+1))
    Spr[i][j] = Nprfunc(i,(j+1),(j+1), True) * Jfunc(i, j, (j+1)) - Jprfunc(i,(j+1),(j+1), True) * Nfunc(i, j, (j+1))

####################################

print('N:', N)
print('Npr:', Npr)
print('J:', J)
print('Jpr:', Jpr)
print('C:', C)
print('Cpr:', Cpr)
print('S:', S)
print('Spr:', Spr)

####################################
# Добавляем в конец массива alfa элемент со значением 0
# Добавляем в конец массива eps элемент со значением длинны массива eps
####################################

if eps[len(eps)-1] != (len(eps)-1):
    alfa.append(0)
    eps.append(len(eps))

print('alfa[', len(alfa)-1,']:', alfa[len(alfa)-1])
print('eps[', len(eps)-1,']:', eps[len(eps)-1])

####################################
# Определяем импедансы (Z) и адмитансы (Y)
####################################

Z = [[0] * (len(a)) for i in range(toch)]
Y = [[0] * (len(a)) for i in range(toch)]

for i in range(toch - 1):
  for h in range(len(a)):

####################################
# Задаем импедансы (Z) и адмитансы (Y) для первого слоя
####################################

    if h == 0:
     Z[i][h] = (cmath.sqrt((cmath.exp(alfa[1] * 1j) * abs(eps[1])) / ((cmath.exp(alfa[0] * 1j) * abs(eps[0]))))) * ((Jpr[i])/(J[i]))
     Y[i][h] = (cmath.sqrt((cmath.exp(alfa[0] * 1j) * abs(eps[0])) / ((cmath.exp(alfa[1] * 1j) * abs(eps[1]))))) * ((Jpr[i])/(J[i]))

####################################
# Задаем импедансы (Z) и адмитансы (Y) для последнего слоя
####################################

    else:
      if h == (len(a) - 1):
       Z[i][h] = (cmath.sqrt((cmath.exp(alfa[h+1] * 1j) * abs(eps[h+1])) / ((cmath.exp(alfa[h] * 1j) * abs(eps[h]))))) * \
                 (Cpr[i][h-1] + Z[i][h-1] * Spr[i][h-1]) / (C[i][h-1] + Z[i][h-1] * S[i][h-1]) / 2
       Y[i][h] = (cmath.sqrt((cmath.exp(alfa[h] * 1j) * abs(eps[h])) / ((cmath.exp(alfa[h+1] * 1j) * abs(eps[h+1]))))) * \
                 (Cpr[i][h-1] + Y[i][h-1] * Spr[i][h-1]) / (C[i][h-1] + Y[i][h-1] * S[i][h-1]) * 2

####################################
# Задаем импедансы (Z) и адмитансы (Y) для промежуточных слоёв
####################################

      else:
         Z[i][h] = (cmath.sqrt((cmath.exp(alfa[h+1] * 1j) * abs(eps[h+1])) / ((cmath.exp(alfa[h] * 1j) * abs(eps[h]))))) * \
                 (Cpr[i][h-1] + Z[i][h-1] * Spr[i][h-1]) / (C[i][h-1] + Z[i][h-1] * S[i][h-1])
         Y[i][h] = (cmath.sqrt((cmath.exp(alfa[h] * 1j) * abs(eps[h])) / ((cmath.exp(alfa[h+1] * 1j) * abs(eps[h+1]))))) * \
                 (Cpr[i][h-1] + Y[i][h-1] * Spr[i][h-1]) / (C[i][h-1] + Y[i][h-1] * S[i][h-1])

print('Z:', Z)
print('Y:', Y)

####################################
# Определяем массивы для mJ, mJpr, mH, mHpr
####################################

mJ = [0 * n for i in range(toch)]
mJpr = [0 * n for i in range(toch)]
mH = [0 * n for i in range(toch)]
mHpr = [0 * n for i in range(toch)]

####################################
# Функция Ханкеоя второго рода как функцию Hfunc(i)
# где: i -- порядок
####################################

def Hfunc(i):

  nu = i + 1
  H = (scipy.special.hankel1(nu + 0.5,k0)) * (math.sqrt(k0 * math.pi/2))
  return H

####################################
# Производная функции Ханкеоя второго рода как функцию Hprfunc(i)
# где: i -- порядок
####################################

def Hprfunc(i):
  nu = i + 1
  Hpr = ((nu / (2 * nu + 1)) * (((scipy.special.hankel1(nu - 0.5,k0) * (cmath.sqrt(k0 * math.pi/2))) / k0)) * k0 - \
  (((nu + 1) / (2 * nu + 1)) * ((scipy.special.hankel1(nu + 1.5,k0)) * (cmath.sqrt(k0 * math.pi/2))) / k0) * k0 + \
  ((scipy.special.hankel1(nu + 0.5,k0)) * (cmath.sqrt(k0 * math.pi/2))) / k0)
  return Hpr

####################################
# Заполняем массивы mJ, mJpr, mH, mHpr
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

print('mJ:', mJ)
print('mJpr:', mJpr)
print('mH:', mH)
print('mHpr:', mHpr)

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

print('Mn:', Mn)
print('Nn:', Nn)

#teta=0.01:pi/180:2*pi; %диапазон изменения углов для декартовой системы координат

i = [1, 2, 3]
teta_start = 0.01
teta_stop = 360
teta_diap = abs(teta_stop)-abs(teta_start)
step = math.pi/360
steps = int(((teta_diap)*(math.pi/180)) / step)
teta = [0 * n for i in range(steps)]
cos_teta = [0 * n for i in range(steps)]
print(step, steps, math.pi)
for i in range(steps):
  if i == 0:
    teta[i] = (teta_start)*(math.pi/180)
  else:
    teta[i] = teta[i-1] + step
print(teta)

for i in range(steps):
  cos_teta[i] = math.cos(teta[i])
print(cos_teta)

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
y=0

for z in range(len(teta)):
  for p in range(toch):
    y=p+1
    E_teta[p][z]=((((2*y+1)/(y*(y+1)))*((-1)**y))*(tay[p][z]*Mn[p]-pii[p][z]*Nn[p]))

  for p in range(toch):
    y=p+1
    P_teta[0][z]=(P_teta[0][z]+E_teta[p][z])
  for p in range(toch):
    P_teta[0][z]=abs(P_teta[0][z])


#print('P_teta:', P_teta)
for z in range((toch)):
    print('E_teta:', E_teta[z][180])

tetay = [0 * n for i in range(steps)]
DN_NORM = [0 * n for i in range(len(teta))]

P_teta_max = 0

for i in range(len(teta)):
  if P_teta[0][i] > P_teta_max:
    P_teta_max = P_teta[0][i]
print(P_teta_max)

for i in range(len(teta)):
  tetay[i] = teta[i]*(steps/math.pi)

tetay.reverse()
print(tetay)

for i in range(len(teta)):
  DN_NORM[i] = 20*math.log10(P_teta[0][i]/P_teta_max)

print(DN_NORM)

Zlay = [0 * n for i in range(len(a))]
Ylay = [0 * n for i in range(len(a))]
Mlay = [[0] * (len(a)) for i in range(toch)]
Nlay = [[0] * (len(a)) for i in range(toch)]
Mnlay = [0 * n for i in range(len(a))]
Nnlay = [0 * n for i in range(len(a))]

for i in range(toch):
  for h in range(len(a)):
    Zlay[h] = (Zlay[h] + (Z[i][h]))
    Ylay[h] = (Ylay[h] + (Y[i][h]))

for i in range(toch):
  for h in range(len(a)):
      nu = i + 1
      Mlay[i][h] = (Z[i][2] * mJ[i] - mJpr[i]) / (Z[i][2] * mH[i] - mHpr[i])
      Mlay[i][h] = Mn[i].real - Mn[i].imag * 1j
      Nlay[i][h] = (Y[i][h] * mJ[i] - mJpr[i]) / (Y[i][h] * mH[i] - mHpr[i])
      Nlay[i][h] = Nn[i].real - Nn[i].imag * 1j
      #print(mJ[i], mJpr[i], Z[i][h], Y[i][h], Mlay[i][h], Nlay[i][h])

for i in range(toch - 1):
  for h in range(len(a)):
     Mnlay[h] = (Mnlay[h] + Mlay[i][h])
     Nnlay[h] = abs(Nnlay[h] + Mn[i])

print(Zlay)
print(Ylay)

import matplotlib.pyplot as plt

plt.plot(tetay,DN_NORM,color='blue', linestyle='-', linewidth=2, label='Sphere Scatterplot')



# Показать сетку
plt.grid(True)

plt.show()

# Создайте поларный график
fig = plt.figure(figsize=(6, 6))  # Размеры фигуры
ax = fig.add_subplot(111, projection='polar')  # Создание полярной системы координат

# Polar Plot
ax.plot(teta, DN_NORM, color='blue', linestyle='-', linewidth=2, label='Sphere Scatterplot')

# Legend
ax.legend(loc='upper right')

plt.title('Polar Scatterplot')
plt.show()