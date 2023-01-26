
"""Apply velocity cutoff to raw temperature, electron heat conduction flux, and density datasets.

This script removes all data below the specified velocity cutoff and removes the velocity data from each array.

Example Usages
--------------
python filter.py --vel_cutoff 650 --raw_dir ../raw/

python filter.py --filt_dir ../filter/
"""


import os
import numpy as np
import argparse


def filter_by_vel(data, vel_lim):
    """
    Filters input dataset by corresponding velocity limit.

    :param data: (n x 3) array; (distance, velocities, parameter)
    :param vel_lim: scalar; velocity cutoff
    :return data: (n x 2) array; (distance, parameter)
    """

    data = data[np.where(data[:, 1] > vel_lim)[0], :]
    data = np.column_stack((data[:, 0], data[:, 2]))

    return data


def main():
    parser = argparse.ArgumentParser(
        description='''Apply velocity cutoff to raw temperature, electron heat conduction flux, and density datasets.'''
    )

    parser.add_argument('--vel_cutoff',
                        type=int, default=600,
                        help='velocity cutoff in km/s (default: 600)')

    parser.add_argument('--raw_dir',
                        type=str, default='../raw/',
                        help='origin directory for unfiltered datasets (default: \'../raw/\')')

    parser.add_argument('--filt_dir',
                        type=str, default='../filter/',
                        help='destination directory for filtered datasets (default: \'../filter/\')')

    args = parser.parse_args()

    if args.raw_dir[-1] != '/':
        args.raw_dir += '/'
    if args.filt_dir[-1] != '/':
        args.filt_dir += '/'

    if not os.path.isdir(args.raw_dir):
        raise Exception("raw_dir directory does not exist.")
    if not os.path.isdir(args.filt_dir):
        create_dir = input('filt_dir directory does not exist. Would you like to create it? (y/n)')
        if create_dir == ("y" or "Y" or "yes" or "Yes" or "YES"):
            os.mkdir(args.filt_dir)
        elif create_dir == ("n" or "N" or "no" or "No" or "NO"):
            raise Exception("Program terminated.")
        else:
            raise Exception("User didn't provide y/n response.")

    print('Applying velocity cutoff...')

    for filename in os.listdir(args.raw_dir+'PSP/'):
        if filename[-4:] == '.csv':
            data = np.concatenate((np.genfromtxt(args.raw_dir+'PSP/'+filename, delimiter=','),
                                   np.genfromtxt(args.raw_dir+'Helios/'+filename, delimiter=','),
                                   np.genfromtxt(args.raw_dir+'Ulysses/'+filename, delimiter=',')))
            data = filter_by_vel(data, args.vel_cutoff)

            assert np.shape(data) != (0, 2), 'Velocity cutoff too large.'

            np.savetxt(args.filt_dir+filename, data, delimiter=",")

    print('Done.')


if __name__ == '__main__':
    main()
