import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Want to create an array of wavelets that can propagate in the z direction

grid = np.mgrid[0:5, 0:5]

def phase_func(x):
    return 2*np.pi*x@x

def psi(r1,t):
    """ Calculates the field at coordinate r1 at time t

    Args:
        r1 (np.array): position where you want to measure
        t1 (float): [description]
    """

    w = 1
    c = 1 # wave speed
    delta_r = 0.01

    def v(r2, t):
        return np.cos((w*t+ phase_func(r2)))
    
    integral = 0
    # Integrate over flat transducer
    for r2 in np.arange(-1,1,delta_r):
        r2 = np.array([r2, 0]) # At z = 0
        dist = np.linalg.norm((r1-r2))
        
        t_ret = dist/c

        integral += (v(r2, t-t_ret)/ (2 * np.pi * dist) * delta_r 
                * np.heaviside(t-t_ret, 1))
    
    return integral

if __name__ == '__main__':
    #psi = np.vectorize(psi)
    print(psi(np.array([5,0]),5.5))
    