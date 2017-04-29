import numpy as np
from scipy import optimize

H = np.array([[2., 0.],
              [0., 8.]])

c = np.array([0, -32])

c0 = 64

A = np.array([[1., 1.],
              [-1., 2.],
              [-1., 0.],
              [0., -1.],
              [0., 1.]])

b = np.array([7., 4., 0., 0., 4.])

x0 = np.random.randn(2)


def loss(x, sign=1.):
    return sign * (0.5 * np.dot(x.T, np.dot(H, x)) + np.dot(c, x) + c0)


def jac(x, sign=1.):
    return sign * (np.dot(x.T, H) + c)


cons = {'type': 'ineq',
        'fun': lambda x: b - np.dot(A, x),
        'jac': lambda x: -A}

opt = {'disp': False}


def solve():
    res_cons = optimize.minimize(loss, x0, jac=jac, constraints=cons,
                                 method='SLSQP', options=opt)

    res_uncons = optimize.minimize(loss, x0, jac=jac, method='SLSQP',
                                   options=opt)

    print('\nConstrained:')
    print(res_cons)

    print('\nUnconstrained:')
    print(res_uncons)

solve()