import numpy as np
import sys

if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    filename = "matrix.txt"

# Считывание матрицы A из файла
file = open("matrix.txt")
n = int(file.readline())                    # В первой строке стоит порядок матрицы n
data = file.read()                          # Далее считываем в строку все числа
A = np.fromstring(string=data, dtype=np.float32, sep=" \n")     # Переводим в одномерный массив
A = np.reshape(A, (n, n))                   # Перегруппируем в двумерный массив n * n

x0 = np.linspace(1, n, n).reshape((n, 1))   # Генерим тестовый столбец x0: ( 1 ... n )T
b = np.dot(A, x0).reshape((n, 1))           # Прогоняем его через матрицу A. Решением системы должен быть наш x0

Ab = np.hstack((A, b))                      # Ab = ( A | b )

for i in range(n):
    if Ab[i, i] == 0:                       # Если наткнулись на ноль на главной диагонали, можно попробовать
        nzs = np.nonzero(Ab[i:, i])         # Поменять местами эту строку и другую пониже, в которой в этом столбце
        nzs = nzs[0]                        # Не ноль
        if len(nzs) == 0:                   # Если такого нет, то в матрице есть нудевой столбец и она вырождена
            print("Матрица вырождена.")
            sys.exit(1)
        nz = nzs[0] + i
        tmp = Ab[nz,:].copy()
        Ab[nz,:] = Ab[i,:]
        Ab[i,:] = tmp

    Ab[i,:] /= Ab[i,i]                      # Делим строку на Aii
    for i1 in range(n):                     # Зануляем все эл-ты столбца, кроме i-того вычитанием домноженных строк
        if i1 == i:
            continue
        Ab[i1,:] -= Ab[i, :] * Ab[i1, i]

x = Ab[:,n]                                 # Матрица слева приведена к единичной, справа приписан столбец решения

print(x)                                    # Печатаем x
print(np.linalg.norm(x - x0.reshape((n))))  # Печатаем |x - x0|
