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


def laske_askel(R, x, d, v, delta, i_max):

    d_pituus = laske_d_pituus(d)
    d = d / d_pituus
    fitness0 = laske_fitness(x, R, v)
    m = 1
    print("Suunta: " + str(d))

    while True:
        alpha = delta*(1.618**m)
        x1 = x + alpha*d
        fitness1 = laske_fitness(x1, R, v)

        if fitness1 > fitness0:

            if m < 2:
                m = 2

            break

        print("Piste = " + str(x1))
        print("Fitness0 = " + str(fitness0))
        print("Fitness1 = " + str(fitness1))
        print("m = " + str(m))
        print("Alpha = " + str(alpha))

        fitness0 = 0 + fitness1
        m += 1

        input("Continue>")

    alpha_y = delta*(1.618**m)
    alpha_a = delta*(1.618**(m-2))

    while True:
        vali = alpha_y - alpha_a

        if vali < i_max:
            return alpha_a

        a = vali*(1 - 0.618)
        b = vali*0.618

        fitness_a = laske_fitness((x + d*(alpha_a + a)), R, v)
        fitness_b = laske_fitness((x + d*(alpha_a + b)), R, v)

        if fitness_a < fitness_b:
            alpha_y = alpha_a + b

        else:
            alpha_a = alpha_a + a


def laske_d(x, g_arvot, nvars):

    g1 = [(1/(4*x[0]*x[1])) - ((x[0] + x[1])/(4*(x[0]**2)*x[1])),
          (1 / (4 * x[0] * x[1])) - ((x[0] + x[1]) / (4 * (x[1] ** 2) * x[0])),
          ]

    g2 = [(1 / (4 * x[0] * x[1])) - ((x[0] - x[1]) / (4 * (x[0] ** 2) * x[1])),
          (- 1 / (4 * x[0] * x[1])) - (x[0] - x[1]) / (4 * x[0] * (x[1] ** 2))
          ]

    c = np.array([3109888511975475/2199023255552, 3109888511975475/2199023255552])
    b = np.array(g_arvot)
    A = np.array([])


def paivita(x, d, alph):
    x = x + alph*d
    return x

def laske_rajoitusrikkomat(x):

    cons = [((2 * (2 ** 0.5) * 1000 * (x[0] + x[1]) * 10000) /
            (210000 * (4 * x[0] * x[1]) + 4 * x[0] * x[1])) - 0.5,
            ((2 * (2 ** 0.5) * 1000 * (x[0] - x[1]) * 10000) /
            (210000 * (4 * x[0] * x[1]))) - 1
            ]

    print(cons)

    return cons

def laske_fitness(x, R, v):

    fitness = (2**0.5)*1000*(x[0] + x[1])

    fitness += R*v

    return fitness


def main():

    nvars = 2
    x0 = [180, 180]
    eps1 = 0.01
    eps2 = 0.01
    R = 1000000
    delta = 10
    i_max = 0.01

    iter_max = 10

    x = x0

    for k in range(iter_max):

        g_arvot = laske_rajoitusrikkomat(x)
        v = laske_v(g_arvot)
        d = laske_d(x, g_arvot, nvars)

        if pysayta(d, v, eps1, eps2):
            print("Optimi pisteessa: " + str(x) + "\n" + "Kierroksella: " + str(k))
            break

        alpha = laske_askel(R, x, d, v, delta, i_max)
        x = paivita(x, d, alpha)


main()