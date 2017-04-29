# CSD optimointialgoritmi
# Taneli Leppänen

from scipy.optimize import minimize
from math import sqrt
import numpy as np


def laske_d_pituus(d):
    d_pituus = 0

    for alkio in d:
        d_pituus += alkio ** 2

    d_pituus = sqrt(d_pituus)

    return d_pituus


def pysayta(d, v, eps1, eps2):

    d_pituus = laske_d_pituus(d)

    if d_pituus < eps1 and v < eps2:
        return True

    else:
        return False


def laske_v(g_lista):

    v = 0

    for g in g_lista:

        if g > v:
            v = 0 + g

    return v


def laske_askel(R, x, gamma, d, v):

    d_pituus = laske_d_pituus(d)
    beta = gamma * (d_pituus ** 2)

    fitness0 = laske_fitness(x, R, v)

    j = 0
    while True:
        t = (0.5**j)
        x1 = x + t*d
        print("x1 = " + str(x1))
        fitness1 = laske_fitness(x1, R, v)
        print("Fitness0 = " + str(fitness0))
        print("Fitness1 = " + str(fitness1))
        print("Fitness1 + t*beta = " + str(fitness1 + t * beta))
        input("Continue>")
        print()

        if fitness1 + t*beta < fitness0:
            return t

        j += 1


def laske_d(x, g_arvo):

    g1 = -212 / (x[0]**2)
    g2 = -164 / (x[1]**2)
    g3 = -116 / (x[2]**2)
    g4 = -68 / (x[3]**2)
    g5 = -1670400 / (x[0]**3)
    g6 = -993600 / (x[1]**3)
    g7 = -489600 / (x[2]**3)
    g8 = -158400 / (x[3]**3)

    fun = lambda d: 7500*d[0] + 7500*d[1] + 7500*d[2] + 7500*x[3] * 0.5*(d[0]**2 + d[1]**2 + d[2]**2 + d[3]**2)

    cons = ({'type': 'ineq', 'fun': lambda d: g1*d[0] + g1*d[1] + g1*d[2] + g1*d[3] + g_arvo[0]},
            {'type': 'ineq', 'fun': lambda d: g2*d[0] + g2*d[1] + g2*d[2] + g2*d[3] + g_arvo[1]},
            {'type': 'ineq', 'fun': lambda d: g3*d[0] + g3*d[1] + g3*d[2] + g3*d[3] + g_arvo[2]},
            {'type': 'ineq', 'fun': lambda d: g4*d[0] + g4*d[1] + g4*d[2] + g4*d[3] + g_arvo[3]},
            {'type': 'ineq', 'fun': lambda d: g5*d[0] + g5*d[1] + g5*d[2] + g5*d[3] + g_arvo[4]},
            {'type': 'ineq', 'fun': lambda d: g6*d[0] + g6*d[1] + g6*d[2] + g6*d[3] + g_arvo[5]},
            {'type': 'ineq', 'fun': lambda d: g7*d[0] + g7*d[1] + g7*d[2] + g7*d[3] + g_arvo[6]},
            {'type': 'ineq', 'fun': lambda d: g8*d[0] + g8*d[1] + g8*d[2] + g8*d[3] + g_arvo[7]})

    res = minimize(fun, (0, 0, 0), method='SLSQP', constraints=cons)

    return res.x


def paivita(x, d, alph):
    x = x + alph*d
    return x

def laske_rajoitusrikkomat(x):

    cons = [(5300 / (25 * x[0])) - 8,
            (4100 / (25 * x[1])) - 8,
            (2900 / (25 * x[2])) - 8,
            (1700 / (25 * x[3])) - 8,
            ((12 * 3480000) / (50 * x[0]**2)) - 40,
            ((12 * 2070000) / (50 * x[1]**2)) - 40,
            ((12 * 1020000) / (50 * x[2]**2)) - 40,
            ((12 * 330000) / (50 * x[3]**2)) - 40]

    return cons


def laske_fitness(x, R, v):

    fitness = 25 * 300 * (x[0] + x[1] + x[3] + x[4])

    fitness += R*v
    print(v)
    print(fitness)

    return fitness


def main():

    x0 = [150, 150, 150, 150]
    eps1 = 0.001
    eps2 = 0.001
    gamma = 0.05
    R = 10000

    iter_max = 10

    x = np.array(x0)

    for k in range(iter_max):

        g_arvot = laske_rajoitusrikkomat(x)
        v = laske_v(g_arvot)
        d = np.array(laske_d(x, g_arvot))

        if pysayta(d, v, eps1, eps2):
            print("Optimi pisteessa: " + str(x) + "\n" + "Kierroksella: " + str(k))
            break

        alph = laske_askel(R, x, gamma, d, v)
        x = paivita(x, d, alph)


main()