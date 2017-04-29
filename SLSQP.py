from scipy.optimize import minimize

def main():

    x = [3,3]

    fun = lambda d: (-d[0] - d[1]) + 0.5*((d[0]**2) + (d[1]**2))

    cons = ({'type': 'ineq', 'fun': lambda d: (1/3)*d[0] + (1/3)*d[1] - (2/3)},
            {'type': 'ineq', 'fun': lambda d: -d[0] - 1},
            {'type': 'ineq', 'fun': lambda d: -d[1] - 1})


    res = minimize(fun, x, method='SLSQP', constraints=cons )

    print(res.x)
    print(res.success)
    print(res.fun)

main()