""" This file preprocesses data collected from Matlab

.mat file -> .npy file
rotated according to the angle the beam was formed in
this configuration

"""

# Normalisation
from scipy.ndimage.interpolation import rotate
import scipy.io as sio
import numpy as np


def import_mat(path):
    """Imports .mat file in the saved kwave format
        .mat file should have 3 saved variables as a time series

    Args:
        path (str): relative filepath to .mat file

    Returns:
        data (np.array): 4 dimensional Numpy array of time series data
    """
    data = sio.loadmat(path)['data'][0][0]
    data = np.array([data[i] for i in range(len(data))])
    return data


def rotate_all(data, angle):
    """Rotates data

    Args:
        data (np.array): 4 dimensional Numpy array of time series data
        angle (_type_): angle to normal to which beam is steered

    Returns:
        _type_: _description_
    """
    return rotate(data, -angle, axes=(2, 1), reshape=False)


if __name__ == '__main__':
    for mode in ['focus', 'focus_wrap']:
        for angle in [0, 10, 20, 30]:
            date = '0302'
            name = '_'.join([date, mode, str(angle)])
            fpath = '../raw_data/' + name + '.mat'

            data = import_mat(fpath)
            rotated = rotate(data, -angle, axes=(2, 1), reshape=False)

            savepath = '../processed_data/' + name + '.npy'
            with open(savepath, 'wb') as f:
                np.save(f, rotated)

    print('Done')
