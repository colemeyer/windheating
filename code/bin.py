
"""Bin filtered temperature, electron heat conduction flux, and density datasets.

This script separates data points within each dataset into 100 averaged bins.

Example Usages
--------------
python bin.py --filt_dir ../filter/
"""


import os
import numpy as np
import argparse
import warnings


def bin_data(data):
    """
    Bins input dataset into 100 points.

    :param data: (n x 2) array; (distance, parameter)
    :return data: (n x 2) array; (distance, parameter)
    """

    num_bins = 100

    # Identify bin edges
    bin_edges = np.logspace(np.log10(0.063), np.log10(5.44), num_bins)

    # Separate data points into bins and return corresponding indices
    bin_indices = np.digitize(data[:, 0], bin_edges[:-1]) - 1

    binned_data = np.ones((num_bins, 2))
    for i in range(num_bins):
        # Set distance value of each bin
        binned_data[i, 0] = 0.5 * (bin_edges[1] - bin_edges[0]) + bin_edges[i]

        # Average across bin to find parameter at bin i
        binned_data[i, 1] = np.mean(data[np.where(bin_indices == i)[0], 1])

    return binned_data[~np.isnan(binned_data[:, 1]), :]


def main():

    parser = argparse.ArgumentParser(
        description='''Bin filtered temperature, electron heat conduction flux, and density datasets.'''
    )

    parser.add_argument('--filt_dir',
                        type=str, default='../filter/',
                        help='origin directory for filtered datasets (default: \'../filter/\')')

    parser.add_argument('--bins_dir',
                        type=str, default='../bins/',
                        help='destination directory for binned datasets (default: \'../bins/\')')

    args = parser.parse_args()

    if args.filt_dir[-1] != '/':
        args.filt_dir += '/'
    if args.bins_dir[-1] != '/':
        args.bins_dir += '/'

    if not os.path.isdir(args.filt_dir):
        raise Exception("filt_dir does not exist.")
    if not os.path.isdir(args.bins_dir):
        create_dir = input('bins_dir directory does not exist. Would you like to create it? (y/n)')
        if create_dir == ("y" or "Y" or "yes" or "Yes" or "YES"):
            os.mkdir(args.bins_dir)
        elif create_dir == ("n" or "N" or "no" or "No" or "NO"):
            raise Exception("Program terminated.")
        else:
            raise Exception("User didn't provide y/n response.")

    print('Binning data...')

    # Hide warnings resulting from empty bins
    warnings.filterwarnings("ignore")
    for filename in os.listdir(args.filt_dir):
        if filename[-4:] == '.csv':
            data = np.genfromtxt(args.filt_dir + filename, delimiter=',')

            assert np.shape(data)[1] == 2, 'Data array wrong size.'

            np.savetxt(args.bins_dir + filename, bin_data(data), delimiter=",")

    print('Done.')


if __name__ == '__main__':
    main()
