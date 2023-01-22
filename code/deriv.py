
"""Take numerical derivative of temperature, electron heat conduction flux, and density fits.

This script takes a numerical derivative of fits using the standard centered-difference approximation.

Example Usages
--------------
python deriv.py --fits_dir ../fits/

python deriv.py --deriv_dir ../deriv/
"""


import os
import numpy as np
import astropy.units as u
import argparse


def calculate_last_term(r, q):
    """
    Calculate the \"last term\" in the electron internal energy conservation equation from Cranmer et al. 2009.

    :param r: (n x 1) array: distance
    :param q: (n x 1) array: electron heat conduction flux
    :return last_term: (n x 1) array: "last term" of electron internal energy conservation equation from Cranmer et
    al. 2009
    """

    colat = (15 * u.deg).to(u.rad)
    rot_freq = 2.7e-6 * u.s ** -1

    dr = r[1] - r[0]

    long_bit = np.cos(np.arctan(rot_freq * r * np.sin(colat) / (700 * u.km * u.s ** -1))) ** 2
    last_term = (np.gradient((q * r ** 2 * long_bit).to(u.erg / u.s), dr) / r ** 2).to(
        u.erg * u.cm ** -3 * u.s ** -1)

    return last_term


def main():
    parser = argparse.ArgumentParser(
        description='''Take numerical derivative of temperature, electron heat conduction flux, and density fits.'''
    )

    parser.add_argument('--fits_dir',
                        type=str, default='../fits/',
                        help='origin directory for fits (default: \'../fits/\')')

    parser.add_argument('--deriv_dir',
                        type=str, default='../deriv/',
                        help='destination directory for derivatives (default: \'../deriv/\')')

    args = parser.parse_args()

    if args.fits_dir[-1] != '/':
        args.fits_dir += '/'
    if args.deriv_dir[-1] != '/':
        args.deriv_dir += '/'

    if not os.path.isdir(args.fits_dir):
        raise Exception("fits_dir directory does not exist.")
    if not os.path.isdir(args.deriv_dir):
        create_dir = input('deriv_dir directory does not exist. Would you like to create it? (y/n)')
        if create_dir == ("y" or "Y" or "yes" or "Yes" or "YES"):
            os.mkdir(args.deriv_dir)
        elif create_dir == ("n" or "N" or "no" or "No" or "NO"):
            raise Exception("Program terminated.")
        else:
            raise Exception("User didn't provide y/n response.")

    print('Taking derivative...')

    r = np.genfromtxt(args.fits_dir + 'r.csv', delimiter=',')
    dr = r[1] - r[0]

    for filename in os.listdir(args.fits_dir):
        if filename[-4:] == '.csv' and filename != 'r.csv':
            if filename == 'q.csv':
                # Take derivative of q
                q = np.genfromtxt(args.fits_dir + filename, delimiter=',')
                dq = np.array(np.gradient(q, dr), dtype='float')
                np.savetxt(args.deriv_dir + 'd' + filename, dq, delimiter=",")

                last_term = calculate_last_term(r * u.AU, q * u.erg * u.cm ** -2 * u.s ** -1)
                np.savetxt(args.deriv_dir + 'last_term.csv', last_term.value, delimiter=",")

            else:
                data = np.genfromtxt(args.fits_dir + filename, delimiter=',')
                ddata = np.array(np.gradient(data, dr), dtype='float')
                np.savetxt(args.deriv_dir + 'd' + filename, ddata, delimiter=",")

    print('Done.')


if __name__ == '__main__':
    main()
