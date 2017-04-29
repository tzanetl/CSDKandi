import numpy as np

c = np.array([1, 2])
d = np.array([2, 4])

fun = np.dot(c.T, d) + 0.5 * (np.dot(d.T, d))

print(fun)