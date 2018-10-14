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

# Метод ПВР:
# D = диагональная, Dii = Aii
# B = D^-1 * A
# B = L + E + U         L - н.треуг, U - верх.треуг, на гл. диаг. нули
#
# x(s+1) = x(s) - w * (L * x(s+1) + (E + U) * x(s) - c)
#
# Суть: сразу же после того, как посчитали i-тую компоненту x, далее подставляем уже ее

aii = A.diagonal().reshape((n, 1))
B = A / aii
c = b / aii

x = b.copy()
r = np.dot(A, x) - b

w = 1.5
eps = 1e-9

while np.linalg.norm(r) > eps:
    for i in range(n):
        x[i] = x[i] - w * (np.dot(B[i, :], x) - c[i])

    np.dot(A, x, r)
    r -= b
    print("|r|\t\t\t= ", np.linalg.norm(r))

print()

print(" x\t\t\t=", x.reshape((n)))
print("|x - x0|\t= ", np.linalg.norm(x - x0))