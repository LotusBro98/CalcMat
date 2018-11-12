# -*- coding: utf-8 -*-

import numpy as np
import sys


def gauss(A, b):
    n = b.shape[0]
    Ab = np.hstack((A, b.reshape(n, 1)))        # Ab = ( A | b )

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
    return x

def J(x):
    return np.vstack([
            [4 - x[3],  -1,         1,          -x[0]],
            [1,         -2,         3 - x[3],   -x[2]],
            [-1,        3 - x[3],   -2,         -x[1]],
            [2 * x[0],  2 * x[1],   2 * x[2],   0    ]
        ])

def F(x):
    return np.asarray([
            4*x[0]      - x[1]      + x[2]      - x[0] * x[3],
            x[0]        - 2 * x[1]  + 3 * x[2]  - x[2] * x[3],
            -x[0]       + 3 * x[1]  - 2 * x[2]  - x[1] * x[3],
            x[0]*x[0]   + x[1]*x[1] + x[2]*x[2] - 1
        ])


def Newton(F, J, x0, eps = 1e-4):
    x = x0
    while True:
        x -= gauss(J(x), F(x))                  # x = x - J^-1(x) * F(x)
        
        df = np.linalg.norm(F(x))
        print (df)
        if (df < eps):
            break

    return x


x0 = np.asarray([0, 0.1, 0.2, 0.3])
print(Newton(F, J, x0, 1e-13))                  # x = (0, (0.5)^(0.5), (0.5)^(0.5), 1)
