# -*- coding: utf-8 -*-

import numpy as np

def f(x):
    return np.sin(x)

def F(x):
    return -np.cos(x)

def integrateSimpson(x0, x1, n):
    I = 0
    grid = np.linspace(x0, x1, n + 1)
    for i in range(n):
        a = grid[i]
        b = grid[i + 1]
        m = (a + b) / 2

        # Формула Симпсона
        I += (b - a) / 6.0 * (f(a) + 4 * f(m) + f(b))

    return I


#Корни полинома Лежандра 3 степени
t1 = 0
t2 = np.sqrt(3.0 / 5.0)
t3 = -np.sqrt(3.0 / 5.0)

#Веса квадратуры
c1 = 2.0 / 3.0 * (1 + 3 * t2 * t3) / (t1 - t2) / (t1 - t3)
c2 = 2.0 / 3.0 * (1 + 3 * t1 * t3) / (t2 - t1) / (t2 - t3)
c3 = 2.0 / 3.0 * (1 + 3 * t1 * t2) / (t3 - t1) / (t3 - t2)

def GaussSum(a, b):
    x1 = (b - a) / 2 * t1 + (a + b) / 2
    x2 = (b - a) / 2 * t2 + (a + b) / 2
    x3 = (b - a) / 2 * t3 + (a + b) / 2

    s = f(x1) * c1 + f(x2) * c2 + f(x3) * c3

    return (b - a) / 2 * s

def integrateGauss(x0, x1, n):
    I = 0
    grid = np.linspace(x0, x1, n + 1)
    for i in range(n):
        a = grid[i]
        b = grid[i + 1]
        
        # Считаем квадратуру Гаусса на каждом отрезочке и складываем
        I += GaussSum(a, b)

    return I

def calcIntegrals(x0, x1, n):
    Int1 = F(x1) - F(x0)
    Int2 = integrateSimpson(x0, x1, n)
    Int3 = integrateGauss(x0, x1, n)

    print "----- n =", n, "-----"
    print "Int1 (первообр.) =", Int1
    print "Int2 (Симпсон)   =", Int2
    print "Int3 (Гаусс)     =", Int3
    print
    print "|Int1 - Int2| =", np.abs(Int1 - Int2)
    print "|Int1 - Int3| =", np.abs(Int1 - Int3)
    print

       
Ns = [10, 20, 40, 80]
x0 = 0
x1 = 1
#x1 = np.pi

for n in Ns:
    calcIntegrals(x0, x1, n)
