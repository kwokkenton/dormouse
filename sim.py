import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Want to create an array of wavelets that can propagate in the z direction

grid = np.mgrid[0:5, 0:5]

def phase_func(x):
    return 1
    return 2*np.pi*x**2

def psi(r1,t):
    """ Calculates the field at coordinate r1 at time t

    Args:
        r1 (np.array): [description]
        t1 (float): [description]
    """

    w = 1
    c = 1 # wave speed
    delta_r = 0.01

    def v(r2,t):
        return np.cos((w*t+ phase_func(r2)))
    
    integral = 0
    # Integrate over space
    for r2 in np.arange(-1,1,delta_r):
        dist = np.linalg.norm((r1-r2))
        t_ret = dist/c
        if t - t_ret < 0:
            continue
        else:
            integral += v(r2, t-t_ret)/ (2 * np.pi * dist) * delta_r
    
    return integral

r1 = 2
##plt.plot(psi(r1,20))
#plt.show()
print(psi([1,2,3],22))
