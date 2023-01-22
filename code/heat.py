
"""Calculate heating rates and error corresponding to standard deviation of residuals of
        temperature, electron heat conduction flux, and density.

This script uses internal energy conservation equations described by Cranmer et al. 2009
to derive proton and electron heating rates. Additionally, the script varies proton and
electron temperature, electron heat conduction flux, and electron density by the standard
deviations of their respective residuals to identify error envelopes of each heating rate.

Example Usages
--------------
python heat.py --fits_dir ../fits/ --heat_dir ../heat/

python heat.py --deriv_dir ../deriv/
"""

from deriv import *
import os
import numpy as np
import astropy.units as u
import argparse


def calculate_err(r, n_e, dn_e, Tp, dTp, Te, dTe, last_term):
    """
    Calculate heating rate errors by varying temperature, electron heat conduction flux, and density by their
    respective standard deviation of residuals.

    :param r: (n x 1) array; distance
    :param n_e: (n x 1) array; electron density
    :param dn_e: (n x 1) array; electron density derivative
    :param Tp: (n x 1) array; proton temperature
    :param dTp: (n x 1) array; proton temperature derivative
    :param Te: (n x 1) array; electron temperature
    :param dTe: (n x 1) array; electron temperature derivative
    :param last_term: (n x 1) array; last term of electron heating rate
    :return Qp_err: (n x 2) array; proton heating rate error (high, low)
    :return Qe_err: (n x 2) array; electron heating rate error (high, low)
    """

    Qp_high = 0
    Qp_low = 0
    Qe_high = 0
    Qe_low = 0
    for i in [-1, 1, 0]:
        for j in [-1, 1, 0]:
            for k in [-1, 1, 0]:
                for f in [-1, 1, 0]:

                    # PSP+Helios+Ulysses standard deviation of residuals
                    Tp_const = (1 + 0.133) ** i
                    Te_const = (1 + 0.122) ** j
                    q_const = (1 + 0.387) ** k
                    n_const = (1 + 0.125) ** f

                    Tp_ = Tp * Tp_const
                    dTp_ = dTp * Tp_const
                    Te_ = Te * Te_const
                    dTe_ = dTe * Te_const
                    n_e_ = n_e * n_const
                    dn_e_ = dn_e * n_const
                    last_term_ = last_term * q_const

                    Qp, Qe = calculate_heat(r, n_e_, dn_e_, Tp_, dTp_, Te_, dTe_, last_term_)

                    if i == -1 and j == -1 and k == -1 and f == -1:
                        Qp_high = Qp.copy()
                        Qp_low = Qp.copy()
                        Qe_high = Qe.copy()
                        Qe_low = Qe.copy()

                    else:
                        if np.sum(Qp) > np.sum(Qp_high):
                            Qp_high = Qp.copy()
                        elif np.sum(Qp) < np.sum(Qp_low):
                            Qp_low = Qp.copy()

                        if np.sum(Qe) > np.sum(Qe_high):
                            Qe_high = Qe.copy()
                        elif np.sum(Qe) < np.sum(Qe_low):
                            Qe_low = Qe.copy()

    Qp_err = np.column_stack((Qp_high, Qp_low))
    Qe_err = np.column_stack((Qe_high, Qe_low))

    return Qp_err, Qe_err


def calculate_heat(r, n_e, dn_e, Tp, dTp, Te, dTe, last_term):
    """
    Calculates heating rates using internal energy conservation equations described by Cranmer et al. 2009.

    :param r: (n x 1) array; distance
    :param n_e: (n x 1) array; electron density
    :param dn_e: (n x 1) array; electron density derivative
    :param Tp: (n x 1) array; proton temperature
    :param dTp: (n x 1) array; proton temperature derivative
    :param Te: (n x 1) array; electron temperature
    :param dTe: (n x 1) array; electron temperature derivative
    :param last_term: (n x 1) array; last term of electron heating rate
    :return Qp: (n x 1) array; proton heating rate
    :return Qe: (n x 1) array; electron heating rate
    """

    n_p = n_e / 1.1
    dn_p = dn_e / 1.1
    kb = 1.381e-23 * u.m ** 2 * u.kg * u.s ** -2 * u.K ** -1
    v = 700 * np.ones(np.shape(n_e)[0]) * u.km / u.s
    vpe = 8.4e-9 * (n_e / (2.5 * u.cm ** -3)) * (Te / (1e5 * u.K)) ** -1.5 / u.s
    vep = 8.4e-9 * (n_p / (2.5 * u.cm ** -3)) * (Tp / (1e5 * u.K)) ** -1.5 / u.s

    Qp_terms = np.ones((np.shape(r)[0], 3))
    Qp_terms[:, 0] = (1.5 * n_p * v * kb * dTp).to(u.erg * u.cm ** -3 * u.s ** -1)
    Qp_terms[:, 1] = (-v * kb * Tp * dn_p).to(u.erg * u.cm ** -3 * u.s ** -1)
    Qp_terms[:, 2] = (-1.5 * n_p * kb * vpe * (Te - Tp)).to(u.erg * u.cm ** -3 * u.s ** -1)
    Qp = Qp_terms[:, 0] + Qp_terms[:, 1] + Qp_terms[:, 2]

    Qe_terms = np.ones((np.shape(r)[0], 4))
    Qe_terms[:, 0] = (1.5 * n_e * v * kb * dTe).to(u.erg * u.cm ** -3 * u.s ** -1)
    Qe_terms[:, 1] = (-v * kb * Te * dn_e).to(u.erg * u.cm ** -3 * u.s ** -1)
    Qe_terms[:, 2] = (-1.5 * n_e * kb * vep * (Tp - Te)).to(u.erg * u.cm ** -3 * u.s ** -1)
    Qe_terms[:, 3] = last_term.to(u.erg * u.cm ** -3 * u.s ** -1)
    Qe = Qe_terms[:, 0] + Qe_terms[:, 1] + Qe_terms[:, 2] + Qe_terms[:, 3]

    return Qp, Qe


def main():
    parser = argparse.ArgumentParser(
        description='''Calculate heating rates and error corresponding to standard deviation of residuals of
        temperature, electron heat conduction flux, and density.'''
    )

    parser.add_argument('--fits_dir',
                        type=str, default='../fits/',
                        help='origin directory for fits (default: \'../fits/\')')

    parser.add_argument('--deriv_dir',
                        type=str, default='../deriv/',
                        help='origin directory for derivatives (default: \'../deriv/\')')

    parser.add_argument('--heat_dir',
                        type=str, default='../heat/',
                        help='destination directory for heating rates (default: \'../heat/\')')

    args = parser.parse_args()

    if args.fits_dir[-1] != '/':
        args.fits_dir += '/'
    if args.deriv_dir[-1] != '/':
        args.deriv_dir += '/'
    if args.heat_dir[-1] != '/':
        args.heat_dir += '/'

    if not os.path.isdir(args.fits_dir):
        raise Exception("fits_dir directory does not exist.")
    if not os.path.isdir(args.deriv_dir):
        raise Exception("deriv_dir directory does not exist.")
    if not os.path.isdir(args.heat_dir):
        create_dir = input('heat_dir directory does not exist. Would you like to create it? (y/n)')
        if create_dir == ("y" or "Y" or "yes" or "Yes" or "YES"):
            os.mkdir(args.heat_dir)
        elif create_dir == ("n" or "N" or "no" or "No" or "NO"):
            raise Exception("Program terminated.")
        else:
            raise Exception("User didn't provide y/n response.")

    print('Loading data...')

    r = np.genfromtxt(args.fits_dir + 'r.csv', delimiter=',')[1:] * u.AU
    n_e = np.genfromtxt(args.fits_dir + 'n_e.csv', delimiter=',')[1:] * u.cm ** -3
    dn_e = np.genfromtxt(args.deriv_dir + 'dn_e.csv', delimiter=',')[1:] * u.cm ** -3 * u.AU ** -1
    Tp = np.genfromtxt(args.fits_dir + 'Tp.csv', delimiter=',')[1:] * u.K
    dTp = np.genfromtxt(args.deriv_dir + 'dTp.csv', delimiter=',')[1:] * u.K * u.AU ** -1
    Te = np.genfromtxt(args.fits_dir + 'Te.csv', delimiter=',')[1:] * u.K
    dTe = np.genfromtxt(args.deriv_dir + 'dTe.csv', delimiter=',')[1:] * u.K * u.AU ** -1
    last_term = np.genfromtxt(args.deriv_dir + 'last_term.csv', delimiter=',')[1:] * u.erg * u.cm ** -3 * u.s ** -1

    print('Calculating heating rates...')
    Qp, Qe = calculate_heat(r, n_e, dn_e, Tp, dTp, Te, dTe, last_term)

    print('Calculating errors...')
    Qp_err, Qe_err = calculate_err(r, n_e, dn_e, Tp, dTp, Te, dTe, last_term)

    print('Calculating Cranmer et al. 2009 heating rates...')
    c_r = r[34:].value  # distances corresponding to Helios and Ulysses data
    c_q = 0.01 * np.exp(- 0.7032 - 2.115 * np.log(c_r) - 0.2545 * np.log(c_r) ** 2) * u.erg * u.cm ** -2 * u.s ** -1
    c_n_e = 2.75 * c_r ** -2 * u.cm ** -3
    c_dn_e = np.gradient(c_n_e, c_r[1] - c_r[0]) * u.AU ** -1
    c_Tp = 1e5 * np.exp(0.9711 - 0.7988 * np.log(c_r) + 0.07062 * np.log(c_r) ** 2) * u.K
    c_dTp = np.gradient(c_Tp, c_r[1] - c_r[0]) * u.AU ** -1
    c_Te = 1e5 * np.exp(0.03460 - 0.4333 * np.log(c_r) + 0.08383 * np.log(c_r) ** 2) * u.K
    c_dTe = np.gradient(c_Te, c_r[1] - c_r[0]) * u.AU ** -1
    c_r *= u.AU
    c_last_term = calculate_last_term(c_r, c_q)

    c_Qp, c_Qe = calculate_heat(c_r, c_n_e, c_dn_e, c_Tp, c_dTp, c_Te, c_dTe, c_last_term)

    np.savetxt(args.heat_dir + 'r.csv', r.value, delimiter=",")
    np.savetxt(args.heat_dir+'Qp.csv', Qp, delimiter=",")
    np.savetxt(args.heat_dir + 'Qp_err.csv', Qp_err, delimiter=",")
    np.savetxt(args.heat_dir+'Qe.csv', Qe, delimiter=",")
    np.savetxt(args.heat_dir+'Qe_err.csv', Qe_err, delimiter=",")
    np.savetxt(args.heat_dir+'Qtot.csv', Qp+Qe, delimiter=",")
    np.savetxt(args.heat_dir+'Qtot_err.csv', Qp_err+Qe_err, delimiter=",")
    np.savetxt(args.heat_dir + 'c_r.csv', c_r.value, delimiter=",")
    np.savetxt(args.heat_dir+'c_Qtot.csv', c_Qp+c_Qe, delimiter=",")

    print('Done.')


if __name__ == '__main__':
    main()
