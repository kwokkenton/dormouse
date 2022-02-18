#%%
# Import relevant packages
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np

#%%
# Define Class to handle plotting and basic slicing
class slate:
    def __init__(self, path):
        self.data = np.array(sio.loadmat(path)['data'])

    def slice(self, N):
        self.sliced = self.data[:,:,N]
        return self.sliced

    def plot_slice(self, N):
        plt.figure(figsize=(8,8))
        Z = self.slice(N)
        X,Y = np.meshgrid(range(Z.shape[0]+1),range(Z.shape[1]+1))
        plt.pcolormesh(X,Y, Z.T, vmin=-1.5, vmax=1.5)
        pass

    def shape(self):
        return self.data.shape
    
    def rotate(self, N, angle):
        self.rotated = rotate(self.slice(N).T, angle=angle, reshape= 1)
        return self.rotated

    def plot_rotated(self):
        plt.figure()
        rotated = self.rotated
        X_new, Y_new  = np.meshgrid(range(rotated.shape[0]+1),range(rotated.shape[1]+1))
        plt.pcolormesh(X_new,Y_new, rotated.T, vmin=-1.5, vmax=1.5)
        return
        

#%%
# Load data from .mat files from Matlab
wrap = slate('../raw_data/p_wrap.mat')
unwrap = slate('../raw_data/p_unwrap.mat')

#%%
# Plot the wrapped and unwrapped files side by side
# Import functions and set backend
from matplotlib.widgets import Slider, Button
plt.set_cmap('RdBu')
%matplotlib

size = 256
N = 800

# Initialise figure
fig, ax = plt.subplots(1,2, figsize = (8,4))
X,Y = np.meshgrid(range(size+1),range(size+1))
quad0 = ax[0].pcolormesh(X,Y, wrap.slice(N).T, vmin=-1.5, vmax=1.5)
quad1 = ax[1].pcolormesh(X,Y, unwrap.slice(N).T, vmin=-1.5, vmax=1.5)

# Set slider for timestep N and adjust spacing 
allowed_N = np.arange(0, 1300, 50)
plt.subplots_adjust(bottom=0.25)

axfreq = plt.axes([0.25, 0.1, 0.65, 0.03])
N_slider = Slider(
    ax=axfreq,
    label='Timestep (N)',
    valmin= 100,
    valmax = 1200,
    valstep=allowed_N,
    valinit= N,
)

# The function to be called anytime a slider's value changes
def update(val):
    quad0.set_array(wrap.slice(val).T)
    quad1.set_array(unwrap.slice(val).T)
    fig.canvas.draw_idle()


# register the update function with each slider
N_slider.on_changed(update)

plt.show()
# %%

from scipy.ndimage.interpolation import rotate
wrap.rotate(N=600, angle = 30)
wrap.plot_rotated()
# %%
unwrap.rotate(N=600, angle = 30)
unwrap.plot_rotated()

# %% 
# Plot along axes
w = wrap.rotate(N=600, angle = 30)

def plot_maxes(w):
    ind = np.unravel_index(np.argmax(w, axis=None), w.shape)


    plt.figure()
    plt.title('Along x')
    plt.plot(w[ind[0],:])
    plt.show()

    plt.figure()
    plt.title('Along y')
    plt.plot(w[:,ind[1]])
    plt.show()

plot_maxes(wrap.rotate(N=600, angle = 30))
plot_maxes(unwrap.rotate(N=600, angle = 30))
# %%
