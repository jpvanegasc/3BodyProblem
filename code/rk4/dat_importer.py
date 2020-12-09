import numpy as np


def get_initial(filename, *args, **kwargs):
    data = np.loadtxt(filename, *args, **kwargs)

    R = data[:, -1]
    m = data[:, -2]
    x0 = data[:, 0:3]
    v0 = data[:, 3:6]

    return R, m, x0, v0
