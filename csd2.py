# CSD optimointialgoritmi
# Taneli Leppänen


from scipy.optimize import minimize
from math import sqrt
import numpy as np
import time

# Tulostaa laskennan ohessa erinäisiä tietoja jos arvo on True
TULOSTA = False


# Laskee suuntavektorin pituuden
def laske_d_pituus(d):
    d_pituus = 0

    for alkio in d:
        d_pituus += alkio ** 2

    d_pituus = sqrt(d_pituus)

    return d_pituus


# Tarkastelee toteutuvatko pyäsytysehdot
def pysayta(d, x, eps1, eps2):

    d_pituus = laske_d_pituus(d)
    v = laske_v(x)

    if TULOSTA is True:
        print("Pysäytysehdot:")
        print("d_pituus = " + str(d_pituus) + ", eps1 = " + str(eps1))
        print("v = " + str(v) + ", eps2 = " + str(eps2) + "\n")

    if d_pituus < eps1 and v < eps2:
        return True

    else:
        return False


# Määrittää sakon suuruuden
def laske_v(x):

    g_arvot = laske_rajoitusrikkomat(x)
    v = 0

    for g in g_arvot:

        if g > v:
            v = 0 + g

    return v


# Laskee askeleen pituuden hakusuuntaan d käyttäen epätarkkaa viivaltahakua
def laske_askel(R, x, gamma, d):

    d_pituus = laske_d_pituus(d)
    beta = gamma * (d_pituus ** 2)

    fitness0 = laske_fitness(x, R)
    j = 0

    while True:
        t = (0.5**j)
        x1 = x + t*d
        fitness1 = laske_fitness(x1, R)

        if TULOSTA is True:
            print("x1 = " + str(x1))
            print("Fitness0 = " + str(fitness0))
            print("Fitness1 = " + str(fitness1))
            print("Fitness1 + t*beta = " + str(fitness1 + t*beta))
            input("Continue> \n")

        if fitness1 + t*beta < fitness0 and fitness1 > 0:
            return t

        j += 1


# Laskee hakusuunnan d QP -tehtävästä
# Funktioon määritellään rajoitusehotojen derivaatat sekä määritellään QP -alitehtävän
# kohdefunktio sekä rajoitusehdot
def laske_d(x, g_arvo):

    g1 = [(1 / (4 * x[0] * x[1])) - ((x[0] + x[1]) / (4 * (x[0] ** 2) * x[1])),
          (1 / (4 * x[0] * x[1])) - ((x[0] + x[1]) / (4 * (x[1] ** 2) * x[0])),
          ]

    g2 = [(1 / (4 * x[0] * x[1])) - ((x[0] - x[1]) / (4 * (x[0] ** 2) * x[1])),
          (- 1 / (4 * x[0] * x[1])) - (x[0] - x[1]) / (4 * x[0] * (x[1] ** 2))
          ]

    fun = lambda d: (3109888511975475 / 2199023255552) * d[0] + (3109888511975475 / 2199023255552) * d[1] + \
                    0.5 * (d[0] ** 2 + d[1] ** 2)

    cons = ({'type': 'ineq', 'fun': lambda d: -(g1[0] * d[0] + g1[1] * d[1]) - g_arvo[0]},
            {'type': 'ineq', 'fun': lambda d: -(g2[0] * d[0] + g2[1] * d[1]) - g_arvo[1]},
            )

    res = minimize(fun, x, method='SLSQP', constraints=cons)

    if TULOSTA is True:
        print("Suunnanhaku:")
        print(res.success)
        print(res.status)
        print("Suunta: " + str(res.x) + "\n")

    return res.x


# Päivittää pisteen hakusuunnan ja askelpituuden perusteella
def paivita(x, d, alph):
    x = x + alph*d
    return x


# Laskee rajoitusrikkomat
# Tähän funktioon määritellään rajoitusehdot
def laske_rajoitusrikkomat(x):

    cons = [((2 * (2 ** 0.5) * 1000 * (x[0] + x[1]) * 10000) /
             (210000 * (4 * x[0] * x[1]) + 4 * x[0] * x[1])) - 0.5,
            ((2 * (2 ** 0.5) * 1000 * (x[0] - x[1]) * 10000) /
             (210000 * (4 * x[0] * x[1]))) - 1
            ]

    return cons


# Laskee pisteen sopivuuden
# Tähän funktioon määritellään kohdefunktio
def laske_fitness(x, R):

    v = laske_v(x)

    fitness = (2 ** 0.5) * 1000 * (x[0] + x[1])

    fitness += R * v

    return fitness



def main():

    t0 = time.time()

    # Käyttäjä asettaa
    x = [120, 120]
    eps1 = 0.5
    eps2 = 0.5
    gamma = 0.4
    R = 1000000
    iter_max = 10

    for k in range(iter_max):

        g_arvot = laske_rajoitusrikkomat(x)
        d = np.array(laske_d(x, g_arvot))

        if pysayta(d, x, eps1, eps2):

            fval = laske_fitness(x, R)
            print("Ratkaisu saavutettu")
            print("Kierros: " + str(k))
            print("x = " + str(x))
            print("fval = " + str(fval))

            for g in range(len(g_arvot)):
                print("g" + str(g + 1) + ": " + str(g_arvot[g]))
            t1 = time.time()

            print("Ajoaika: " + str(t1 - t0))
            exit()

        alph = laske_askel(R, x, gamma, d)
        x = paivita(x, d, alph)

    print("Optimipistettä ei löydettty")
    print("Päätöspiste: " + str(x))

    t1 = time.time()

    print("Ajoaika: " + str(t1 - t0))

main()