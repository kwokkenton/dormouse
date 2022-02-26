""" This file preprocesses data collected from Matlab

.mat file -> .npy file
rotated according to the angle the beam was formed in
this 
"""

# Normalisation
from scipy.ndimage.interpolation import rotate
import scipy.io as sio
import numpy as np


def import_mat(path):
    """imports .mat file in the saved kwave format
        .mat file should have 3 saved variables as a time series

    Args:
        path (str): relative filepath to .mat file

    Returns:
        data (np.array): 4 dimensional Numpy array of time series data
    """
    data = sio.loadmat(path)['data'][0][0]
    data = np.array([data[0], data[1], data[2]])
    return data


def rotate_all(data, angle):
    return rotate(data, -angle, axes=(2, 1), reshape=False)


for mode in ['focus', 'focus_wrap']:
    for angle in [0, 10, 20, ]:
        date = '0222'
        name = '_'.join([date, mode, str(angle)])
        fpath = '../raw_data/' + name + '.mat'

        data = import_mat(fpath)
        rotated = rotate(data, -angle, axes=(2, 1), reshape=False)

        savepath = '../processed_data/' + name + '.npy'
        with open(savepath, 'wb') as f:
            np.save(f, rotated)

print('Done')
