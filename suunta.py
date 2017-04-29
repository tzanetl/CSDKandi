import numpy as np
from scipy import optimize

def main(c, b, nvars, x0):

    H = np.identity(nvars) * 2

    def fun(d, c):
        return np.dot(c.T, d) + 0.5 * (np.dot(d.T, d))


    def jac(d):
        return c + d


    cons = {'type': 'ineq',
            'fun': lambda d: b - np.dot(A, d),
            'jac': lambda d: -A}


    def ratkaiseQP():

        res = optimize.minimize(fun, x0, jac=jac, constraints=cons, method='SLSQP')

