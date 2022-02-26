# %%
# Import relevant packages
from scipy.signal import correlate
from scipy.ndimage.interpolation import rotate
from matplotlib.widgets import Slider, Button
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np

# %%
# Define Class to handle plotting and basic slicing


class slate:
    def __init__(self, path):
        data = sio.loadmat(path)['data']

        # Import Matlab file and store in array
        self.p = np.array(data['p'][0][0])
        self.ux = np.array(data['ux'][0][0])
        self.uy = np.array(data['uy'][0][0])
        self.data = {'p': self.p, 'ux': self.ux, 'uy': self.uy}

    def slice(self, N):
        self.sliced = {'p': 0, 'ux': 0, 'uy': 0}
        for i in self.data.keys():
            self.sliced[i] = self.data[i][:, :, N]
        return self.sliced

    def plot_slice(self, N, var='p'):
        plt.figure(figsize=(8, 8))
        Z = self.slice(N)[var]
        X, Y = np.meshgrid(range(Z.shape[0]+1), range(Z.shape[1]+1))
        plt.pcolormesh(X, Y, Z.T, vmin=-1.5, vmax=1.5)
        pass

    def shape(self):
        return self.data.shape

    def rotate(self, N, angle, var='p', mode='single'):
        if mode == 'single':
            self.rotated = rotate(self.slice(N)[var].T, angle=angle, reshape=1)
        else:
            self.p_rot = rotate(wrap.data['p'], angle, reshape=True)
            self.ux_rot = rotate(wrap.data['ux'], angle, reshape=True)
            self.uy_rot = rotate(wrap.data['uy'], angle, reshape=True)
            self.rotated = {'p': self.p_rot,
                            'ux': self.ux_rot, 'uy': self.uy_rot}

        return self.rotated

    def plot_rotated(self):
        plt.figure()
        rotated = self.rotated
        X_new, Y_new = np.meshgrid(
            range(rotated.shape[0]+1), range(rotated.shape[1]+1))
        plt.pcolormesh(X_new, Y_new, rotated.T, vmin=-1.5, vmax=1.5)
        return


# %%
# Load data from .mat files from Matlab
date = '20220221'
wrap = slate(f'../raw_data/{date}_focus_wrap.mat')
unwrap = slate(f'../raw_data/{date}_focus.mat')

# %%
# Plot the wrapped and unwrapped files side by side
# Import functions and set backend
mode = 'ux'

plt.set_cmap('RdBu')
%matplotlib

size = 256
N = 800

# Initialise figure
fig, ax = plt.subplots(1, 2, figsize=(8, 4))
X, Y = np.meshgrid(range(size+1), range(size+1))
quad0 = ax[0].pcolormesh(X, Y, wrap.slice(N)[mode].T)
quad1 = ax[0].pcolormesh(X, Y, unwrap.slice(N)[mode].T)

# Set slider for timestep N and adjust spacing
allowed_N = np.arange(0, 1300, 50)
plt.subplots_adjust(bottom=0.25)

axfreq = plt.axes([0.25, 0.1, 0.65, 0.03])
N_slider = Slider(
    ax=axfreq,
    label='Timestep (N)',
    valmin=100,
    valmax=1200,
    valstep=allowed_N,
    valinit=N,
)

# The function to be called anytime a slider's value changes


def update(val):
    quad0.set_array(wrap.slice(val)[mode].T)
    quad1.set_array(unwrap.slice(val)[mode].T)
    fig.canvas.draw_idle()


# register the update function with each slider
N_slider.on_changed(update)
plt.show()
# %%

wrap.rotate(N=600, angle=30)
wrap.plot_rotated()
# %%
unwrap.rotate(N=600, angle=30)
unwrap.plot_rotated()

# %%
# Plot along axes, to make picture clearer
w = wrap.rotate(N=600, angle=30)


def plot_maxes(w):
    ind = np.unravel_index(np.argmax(w, axis=None), w.shape)

    plt.figure()
    plt.title('Along x')
    plt.plot(w[ind[0], :])
    plt.show()

    plt.figure()
    plt.title('Along y')
    plt.plot(w[:, ind[1]])
    plt.show()


plot_maxes(wrap.rotate(N=600, angle=30))
plot_maxes(unwrap.rotate(N=600, angle=30))

# %% Implement the Cross Correlation
# Cross correlation calculates the locations where two signals are
# similar
plt.close('all')
N = 300
a1 = wrap.rotate(N, 30)
a2 = unwrap.rotate(N, 30)

fig, ax = plt.subplots(1, 3, figsize=(12, 4))

X, Y = np.meshgrid(range(a1.shape[0]+1), range(a1.shape[0]+1))
ax[0].pcolormesh(X, Y, correlate(a1, a2, 'same'))
ax[1].pcolormesh(X, Y, a1)
ax[2].pcolormesh(X, Y, a2)

# Comparing the cross correlation and the auto-correlation
fig, ax = plt.subplots(1, 2, figsize=(8, 4))
X, Y = np.meshgrid(range(a1.shape[0]+1), range(a1.shape[0]+1))
ax[0].pcolormesh(X, Y, correlate(a1, a1, 'same'))
ax[1].pcolormesh(X, Y, correlate(a2, a2, 'same'))

# TODO Try and quantify the spread

# %% Implement focus tightness measurement
rotate(wrap.data['p'][:, :, 0:2], 30, reshape=True).shape

# %%
(np.average(np.sqrt(wrap.data['p']*wrap.data['p']), axis=2)).shape

plt.imshow((np.average(np.sqrt(unwrap.data['p']*unwrap.data['p']), axis=2)))
plt.imshow((np.average(np.sqrt(wrap.data['p']*wrap.data['p']), axis=2)))
# %%
