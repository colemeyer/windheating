
"""Fit temperature, electron heat conduction flux, and density datasets to polynomials.

This script fits input datasets to polynomials and outputs resulting coefficients and fit data.

Example Usages
--------------
python fit.py --bins_dir ../bins/
"""


import os
import numpy as np
import argparse


def fit_data(param_type, data):
    """
    Fits input dataset to polynomials.

    :param param_type: string; options: ['q','n_e','Tp','Te']
    :param data: (n x 2) array; (distance, parameter)
    :return data: (n x 2) array; (distance, parameter)
    """

    r = np.linspace(0.063, 5.44, 1001)

    if param_type == 'q':
        coeffs = np.polyfit(np.log(data[:, 0]), np.log(data[:, 1] / 0.01), 2)
        fit = 0.01 * np.exp(coeffs[0] * np.log(r) ** 2 + coeffs[1] * np.log(r) + coeffs[2])
    elif param_type == 'n_e':
        coeffs = np.polyfit(np.log(data[:, 0]), np.log(data[:, 1]), 1)
        fit = np.exp(coeffs[0] * np.log(r) + coeffs[1])
    else:
        coeffs = np.polyfit(np.log(data[:, 0]), np.log(data[:, 1] / 1e5), 2)
        fit = 1e5 * np.exp(coeffs[0] * np.log(r) ** 2 + coeffs[1] * np.log(r) + coeffs[2])

    return fit, coeffs


def main():

    parser = argparse.ArgumentParser(
        description='''Fit temperature, electron heat conduction flux, and density datasets to polynomials.'''
    )

    parser.add_argument('--bins_dir',
                        type=str, default='../bins/',
                        help='origin directory for binned datasets (default: \'../bins/\')')

    parser.add_argument('--fits_dir',
                        type=str, default='../fits/',
                        help='destination directory for fits datasets (default: \'../fits/\')')

    args = parser.parse_args()

    if args.bins_dir[-1] != '/':
        args.bins_dir += '/'
    if args.fits_dir[-1] != '/':
        args.fits_dir += '/'

    if not os.path.isdir(args.bins_dir):
        raise Exception("bins_dir directory does not exist.")
    if not os.path.isdir(args.fits_dir):
        create_dir = input('fits_dir directory does not exist. Would you like to create it? (y/n)')
        if create_dir == ("y" or "Y" or "yes" or "Yes" or "YES"):
            os.mkdir(args.fits_dir)
        elif create_dir == ("n" or "N" or "no" or "No" or "NO"):
            raise Exception("Program terminated.")
        else:
            raise Exception("User didn't provide y/n response.")

    print('Fitting data...')

    for filename in os.listdir(args.bins_dir):
        if filename[-4:] == '.csv':
            data = np.genfromtxt(args.bins_dir + filename, delimiter=',')

            assert np.shape(data)[1] == 2, 'Data array wrong size.'

            fit, coeffs = fit_data(filename[:-4], data)
            print(filename[:-4]+': '+str(coeffs))
            np.savetxt(args.fits_dir + filename, fit, delimiter=",")

    np.savetxt(args.fits_dir + 'r.csv', np.array(np.linspace(0.063, 5.44, 1001), dtype='float'), delimiter=",")

    print('Done.')


if __name__ == '__main__':
    main()
