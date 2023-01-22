
"""Plot fits + binned data, heating rates, heating rate ratio, and/or total heating rate.

Example Usages
--------------
python plot.py --bins_dir ../bins/ --plot_fits --save_plots

python plot.py --plot_fits
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import matplotlib as mpl
import argparse

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'serif'


def plot_fits(args, r, q, q_fit, n_e, n_e_fit, Tp, Tp_fit, Te, Te_fit):
    """
    Plot and show binned data + fits plots.

    :param args: (1 x 10) array: (bins_dir, fits_dir, heat_dir, plots_dir, plot_fits, plot_heats, plot_heat_ratio,
    plot_total_heat, plot_all, save_plots)
    :param r: (n x 1) array: distance
    :param q: (n x 2) array: (distance, binned electron heat conduction flux)
    :param q_fit: (n x 1) array: electron heat conduction flux fits
    :param n_e: (n x 2) array: (distance, binned electron density)
    :param n_e_fit: (n x 1) array: electron density fits
    :param Tp: (n x 2) array: (distance, binned proton temperature)
    :param Tp_fit: (n x 1) array: proton temperature fits
    :param Te: (n x 2) array: (distance, binned electron temperature)
    :param Te_fit: (n x 1) array: electron temperature fits
    """

    # Initial set-up
    fig, axs = plt.subplots(3, 1, figsize=(7, 14))
    plt.setp(axs, xscale='log', yscale='log', xlim=[0.055, 6])
    for i in range(3):
        axs[i].tick_params(axis='both', which='major', direction='in', top=True, right=True, length=6, width=1,
                           labelsize=26, pad=10)
        axs[i].tick_params(axis='both', which='minor', direction='in', top=True, right=True, length=3, width=1,
                           labelsize=26)
    axs[0].set_ylabel(r'$\rm T\;(K)$', size=26)
    axs[1].set_ylabel(r'$\rm q_{||,e}\;(erg\;cm^{-2}\;s^{-1})$', size=26)
    axs[2].set_ylabel(r'$\rm n_e\;(cm^{-3})$', size=26)
    axs[2].set_xlabel(r"$\rm r\;(AU)$", size=26)
    axs[0].set_ylim([1.5e4, 5e6])

    # Dotted lines, text, and legends
    axs[0].text(0.07, 85694.8, r"$\rm PSP\;data$", c='k', fontsize=24, family="serif")
    axs[0].text(0.35, 1564567.3, r"$\rm Helios+Ulysses\;data$", c='k', fontsize=24, family="serif")
    axs[0].text(0.065, 26815.0, r"$\rm (a)$", c='k', fontsize=24, family="serif")
    axs[0].text(0.115, 2.2e6, r"$\rm T_p$", c='r', fontsize=30, family="serif")
    axs[0].text(0.85, 4e4, r"$\rm T_e$", c='b', fontsize=30, family="serif")
    axs[1].text(0.07, 0.000649, r"$\rm PSP\;data$", c='k', fontsize=24, family="serif")
    axs[1].text(0.35, 0.604, r"$\rm Helios+Ulysses\;data$", c='k', fontsize=24, family="serif")
    axs[1].text(0.065, 4.216e-5, r"$\rm (b)$", c='k', fontsize=24, family="serif")
    axs[2].text(0.07, 0.773, r"$\rm PSP\;data$", c='k', fontsize=24, family="serif")
    axs[2].text(0.35, 338.436, r"$\rm Helios+Ulysses\;data$", c='k', fontsize=24, family="serif")
    axs[2].text(0.065, 0.0679, r"$\rm (c)$", c='k', fontsize=24, family="serif")

    dotted_line = np.column_stack((0.26 * np.ones(5), np.linspace(17855.7, 4200332.2, 5)))
    axs[0].plot(dotted_line[:, 0], dotted_line[:, 1], c='k', linestyle=':')

    dotted_line = np.column_stack((0.26 * np.ones(5), np.linspace(1.619e-5, 6.176, 5)))
    axs[1].plot(dotted_line[:, 0], dotted_line[:, 1], c='k', linestyle=':')

    dotted_line = np.column_stack((0.26 * np.ones(5), np.linspace(0.029, 2675.43, 5)))
    axs[2].plot(dotted_line[:, 0], dotted_line[:, 1], c='k', linestyle=':')

    # Tick labels
    axs[0].set_xticks([0.1, 1])
    axs[0].set_xticklabels(['', ''], family="serif")
    axs[0].set_yticks([1e5, 1e6])
    axs[0].set_yticklabels([r'$10^{5}$', r'$10^{6}$'], family="serif")
    axs[1].set_xticks([0.1, 1])
    axs[1].set_xticklabels(['', ''], family="serif")
    axs[1].set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    axs[1].set_yticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$'], family="serif")
    axs[2].set_xticks([0.1, 1])
    axs[2].set_xticklabels([r'$10^{-1}$', r'$10^{0}$'], family="serif")
    axs[2].set_yticks([1e-1, 1e0, 1e1, 1e2, 1e3])
    axs[2].set_yticklabels([r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$', r'$10^{3}$'], family="serif")

    # Plot and show data
    axs[0].scatter(Tp[np.where(Tp[:, 0] < 0.25)[0], 0], Tp[np.where(Tp[:, 0] < 0.25)[0], 1], c='r', s=35, marker='*')
    axs[0].scatter(Te[np.where(Te[:, 0] < 0.25)[0], 0], Te[np.where(Te[:, 0] < 0.25)[0], 1], c='b', s=35, marker='*')
    axs[0].scatter(Tp[np.where((Tp[:, 0] > 0.25) & (Tp[:, 0] < 1))[0], 0],
                   Tp[np.where((Tp[:, 0] > 0.25) & (Tp[:, 0] < 1))[0], 1], c='r', s=25, marker='s')
    axs[0].scatter(Te[np.where((Te[:, 0] > 0.25) & (Te[:, 0] < 1))[0], 0],
                   Te[np.where((Te[:, 0] > 0.25) & (Te[:, 0] < 1))[0], 1], c='b', s=25, marker='s')
    axs[0].scatter(Tp[np.where(Tp[:, 0] > 1)[0], 0], Tp[np.where(Tp[:, 0] > 1)[0], 1], c='r', s=30, marker='v')
    axs[0].scatter(Te[np.where(Te[:, 0] > 1)[0], 0], Te[np.where(Te[:, 0] > 1)[0], 1], c='b', s=30, marker='v')
    axs[1].scatter(q[np.where(q[:, 0] < 0.25)[0], 0], q[np.where(q[:, 0] < 0.25)[0], 1], c='k', s=35, marker='*')
    axs[1].scatter(q[np.where((q[:, 0] > 0.25) & (q[:, 0] < 1))[0], 0],
                   q[np.where((q[:, 0] > 0.25) & (q[:, 0] < 1))[0], 1], c='k', s=25, marker='s')
    axs[1].scatter(q[np.where(q[:, 0] > 1)[0], 0], q[np.where(q[:, 0] > 1)[0], 1], c='k', s=30, marker='v')
    axs[2].scatter(n_e[np.where(n_e[:, 0] < 0.25)[0], 0], n_e[np.where(n_e[:, 0] < 0.25)[0], 1], c='k', s=35,
                   marker='*')
    axs[2].scatter(n_e[np.where((n_e[:, 0] > 0.25) & (n_e[:, 0] < 1))[0], 0],
                   n_e[np.where((n_e[:, 0] > 0.25) & (n_e[:, 0] < 1))[0], 1], c='k', s=25, marker='s')
    axs[2].scatter(n_e[np.where(n_e[:, 0] > 1)[0], 0], n_e[np.where(n_e[:, 0] > 1)[0], 1], c='k', s=30, marker='v')

    axs[0].plot(r, Tp_fit, c='k', linestyle='-', linewidth=2)
    axs[0].plot(r, Te_fit, c='k', linestyle='-', linewidth=2)
    axs[1].plot(r, q_fit, c='k', linestyle='-', linewidth=2)
    axs[2].plot(r, n_e_fit, c='k', linestyle='-', linewidth=2)

    plt.tight_layout()
    if args.save_plots:
        plt.savefig(args.plots_dir + 'RawData.pdf', facecolor='white', dpi=300)
    plt.show()


def plot_heat(args, r, Qp, Qp_err, Qe, Qe_err):
    """
    Plot and show heating rate plot.

    :param args: (1 x 10) array: (bins_dir, fits_dir, heat_dir, plots_dir, plot_fits, plot_heats, plot_heat_ratio,
    plot_total_heat, plot_all, save_plots)
    :param r: (n x 1) array: distance
    :param Qp: (n x 1) array: proton heating rate
    :param Qp_err: (n x 2) array: proton heating rate error (high, low)
    :param Qe: (n x 1) array: electron heating rate
    :param Qe_err: (n x 2) array: electron heating rate error (high, low)
    """

    # Initial set-up
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.setp(ax, xscale='log', yscale='log', xlim=[0.06, 6], ylim=[1e-19, 7e-11])
    ax.tick_params(axis='both', which='major', direction='in', top=True, right=True, length=6, width=1, labelsize=26,
                   pad=10)
    ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True, length=3, width=1, labelsize=26)
    ax.set_ylabel(r'$\rm Q\;(erg\;cm^{-3}\;s^{-1})$', size=26)
    ax.set_xlabel(r"$\rm r\;(AU)$", size=26)

    # Tick labels
    ax.set_xticks([0.1, 1])
    ax.set_xticklabels([r'$10^{-1}$', r'$10^{0}$'], family="serif")
    ax.set_yticks([1e-19, 1e-17, 1e-15, 1e-13, 1e-11])
    ax.set_yticklabels([r'$10^{-19}$', r'$10^{-17}$', r'$10^{-15}$', r'$10^{-13}$', r'$10^{-11}$'], family="serif")

    # Plot data
    ax.plot(r, Qp, c='r', linewidth=2, linestyle='-')
    ax.plot(r, Qe, c='b', linewidth=2, linestyle='-')
    ax.fill_between(r, Qp_err[:, 0], Qp_err[:, 1], facecolor='r', alpha=0.2, interpolate=True)
    ax.fill_between(r, Qe_err[:, 0], Qe_err[:, 1], facecolor='b', alpha=0.2, interpolate=True)

    # Configure legend
    p1, = ax.plot([], [], c='r', linewidth=3)
    p2, = ax.plot([], [], c='b', linewidth=3)
    ax.legend([p1, p2], [r'Proton Heating', r'Electron Heating'], handler_map={tuple: HandlerTuple(ndivide=0)},
              fontsize=22, frameon=False, loc=1, bbox_to_anchor=(0, 0, 1, 1), labelspacing=1)

    plt.tight_layout()
    if args.save_plots:
        plt.savefig(args.plots_dir + 'Heat.pdf', facecolor='white', dpi=300)
    plt.show()


def plot_heat_ratio(args, r, Qp, Qp_err, Qe, Qe_err):
    """
    Plot and show heating rate plot.

    :param args: (1 x 10) array: (bins_dir, fits_dir, heat_dir, plots_dir, plot_fits, plot_heats, plot_heat_ratio,
    plot_total_heat, plot_all, save_plots)
    :param r: (n x 1) array: distance
    :param Qp: (n x 1) array: proton heating rate
    :param Qp_err: (n x 2) array: proton heating rate error (high, low)
    :param Qe: (n x 1) array: electron heating rate
    :param Qe_err: (n x 2) array: electron heating rate error (high, low)
    """

    # Initial set-up
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.setp(ax, xscale='log', xlim=[0.06, 6], ylim=[0.35, 0.95])
    ax.tick_params(axis='both', which='major', direction='in', top=True, right=True, length=6, width=1, labelsize=26,
                   pad=10)
    ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True, length=3, width=1, labelsize=26)
    ax.set_ylabel(r'$\rm Q_p/(Q_p + Q_e)$', size=26)
    ax.set_xlabel(r"$\rm r\;(AU)$", size=26)

    # Tick labels
    ax.set_xticks([0.1, 1])
    ax.set_xticklabels([r'$10^{-1}$', r'$10^{0}$'], family="serif")
    ax.set_yticks([0.4, 0.6, 0.8])
    ax.set_yticklabels([r'$0.4$', r'$0.6$', r'$0.8$'], family="serif")

    # Plot data
    ax.plot(r, Qp / (Qp + Qe), c='k', linewidth=3, linestyle='-')
    ax.fill_between(r, Qp_err[:, 0] / (Qp_err[:, 0] + Qe_err[:, 1]), Qp_err[:, 1] / (Qp_err[:, 1] + Qe_err[:, 0]),
                    facecolor='r', alpha=0.1, interpolate=True)

    # Configure legend
    p1, = ax.plot([], [], c='k', linewidth=3, linestyle='-')
    ax.legend([p1], [r'Derived Heating Rates'], handler_map={tuple: HandlerTuple(ndivide=0)},
              fontsize=20, frameon=False, loc=3, bbox_to_anchor=(-0.01, 0.13, 1, 1), labelspacing=1)

    plt.tight_layout()
    if args.save_plots:
        plt.savefig(args.plots_dir + 'HeatRatio.pdf', facecolor='white', dpi=300)
    plt.show()


def plot_total_heat(args, r, c_r, Qtot, Qtot_err, c_Qtot):
    """
    Plot and show heating rate plot.

    :param args: (1 x 10) array: (bins_dir, fits_dir, heat_dir, plots_dir, plot_fits, plot_heats, plot_heat_ratio,
    plot_total_heat, plot_all, save_plots)
    :param r: (n x 1) array: distance
    :param c_r: (n x 1) array: Cranmer et al. 2009 heating rate distance
    :param Qtot: (n x 1) array: total heating rate (Qp + Qe)
    :param Qtot_err: (n x 2) array: total heating rate error (high, low)
    :param c_Qtot: (n x 1) array: Cranmer et al. 2009 total heating rate
    """

    # Initial set-up
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.setp(ax, xscale='log', yscale='log', xlim=[0.06, 6], ylim=[7e-19, 3e-10])
    ax.tick_params(axis='both', which='major', direction='in', top=True, right=True, length=6, width=1, labelsize=26,
                   pad=10)
    ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True, length=3, width=1, labelsize=26)
    ax.set_ylabel(r'$\rm Q\;(erg\;cm^{-3}\;s^{-1})$', size=26)
    ax.set_xlabel(r"$\rm r\;(AU)$", size=26)

    # Tick labels
    ax.set_xticks([0.1, 1])
    ax.set_xticklabels([r'$10^{-1}$', r'$10^{0}$'], family="serif")
    ax.set_yticks([1e-17, 1e-15, 1e-13, 1e-11])
    ax.set_yticklabels([r'$10^{-17}$', r'$10^{-15}$', r'$10^{-13}$', r'$10^{-11}$'], family="serif")

    # Plot data
    ax.plot(r, Qtot, c='k', linewidth=2, linestyle='-')
    ax.plot(c_r, c_Qtot, c='k', linewidth=2, linestyle='--')
    ax.fill_between(r, Qtot_err[:, 0], Qtot_err[:, 1], facecolor='r', alpha=0.2, interpolate=True)

    # Configure legend
    p1, = ax.plot([], [], c='k', linewidth=3, linestyle='-')
    p2, = ax.plot([], [], c='k', linewidth=3, linestyle='--')
    p3, _, _ = ax.errorbar([0.167, 0.251], [1.2e-12, 9.4e-14], fmt='o', c='k', marker='*', markersize=17,
                           markerfacecolor='fuchsia')
    first_legend = ax.legend([p1, p2], [r'Derived Heating', r'Cranmer et al. 2009'],
                             handler_map={tuple: HandlerTuple(ndivide=0)}, fontsize=20, frameon=False, loc=1,
                             bbox_to_anchor=(0, -0.02, 1.02, 1), labelspacing=1)
    ax.add_artist(first_legend)
    ax.legend([p3], ["Bandyopadhyay\net al. 2020"], handler_map={tuple: HandlerTuple(ndivide=0)},
              fontsize=22, frameon=False, loc=3, bbox_to_anchor=(-0.05, 0.07, 1.02, 1), labelspacing=1)

    plt.tight_layout()
    if args.save_plots:
        plt.savefig(args.plots_dir + 'TotalHeat.pdf', facecolor='white', dpi=300)
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='''Plot fits + binned data, heating rates, heating rate ratio, and/or total heating rate.'''
    )

    parser.add_argument('--bins_dir',
                        type=str, default='../bins/',
                        help='origin directory for binned datasets (default: \'../bins/\')')

    parser.add_argument('--fits_dir',
                        type=str, default='../fits/',
                        help='origin directory for fits datasets (default: \'../fits/\')')

    parser.add_argument('--heat_dir',
                        type=str, default='../heat/',
                        help='origin directory for heat datasets (default: \'../heat/\')')

    parser.add_argument('--plots_dir',
                        type=str, default='../plots/',
                        help='destination directory for plots (default: \'../plots/\')')

    parser.add_argument('--plot_fits',
                        action='store_true', default=False,
                        help='plot fits? (default: False)')

    parser.add_argument('--plot_heat',
                        action='store_true', default=False,
                        help='plot heating rates? (default: False)')

    parser.add_argument('--plot_heat_ratio',
                        action='store_true', default=False,
                        help='plot heating rate ratio? (default: False)')

    parser.add_argument('--plot_total_heat',
                        action='store_true', default=False,
                        help='plot total heating rate? (default: False)')

    parser.add_argument('--plot_all',
                        action='store_true', default=False,
                        help='plot all options? (default: False)')

    parser.add_argument('--save_plots',
                        action='store_true', default=False,
                        help='save plots? (default: False)')

    args = parser.parse_args()

    if args.bins_dir[-1] != '/':
        args.bins_dir += '/'
    if args.fits_dir[-1] != '/':
        args.fits_dir += '/'
    if args.heat_dir[-1] != '/':
        args.heat_dir += '/'
    if args.plots_dir[-1] != '/':
        args.plots_dir += '/'

    if not os.path.isdir(args.bins_dir):
        raise Exception("bins_dir directory does not exist.")
    if not os.path.isdir(args.fits_dir):
        raise Exception("fits_dir directory does not exist.")
    if not os.path.isdir(args.heat_dir):
        raise Exception("heat_dir directory does not exist.")
    if not os.path.isdir(args.plots_dir) and args.save_plots:
        create_dir = input('plots_dir directory does not exist. Would you like to create it? (y/n)')
        if create_dir == ("y" or "Y" or "yes" or "Yes" or "YES"):
            os.mkdir(args.plots_dir)
        elif create_dir == ("n" or "N" or "no" or "No" or "NO"):
            raise Exception("Program terminated.")
        else:
            raise Exception("User didn't provide y/n response.")

    if args.plot_fits or args.plot_heat or args.plot_heat_ratio or args.plot_total_heat or args.plot_all:
        print('Plotting data...')

    # Load datasets
    q = np.genfromtxt(args.bins_dir + 'q.csv', delimiter=',')
    n_e = np.genfromtxt(args.bins_dir + 'n_e.csv', delimiter=',')
    Tp = np.genfromtxt(args.bins_dir + 'Tp.csv', delimiter=',')
    Te = np.genfromtxt(args.bins_dir + 'Te.csv', delimiter=',')

    r = np.genfromtxt(args.fits_dir + 'r.csv', delimiter=',')
    q_fit = np.genfromtxt(args.fits_dir + 'q.csv', delimiter=',')
    n_e_fit = np.genfromtxt(args.fits_dir + 'n_e.csv', delimiter=',')
    Tp_fit = np.genfromtxt(args.fits_dir + 'Tp.csv', delimiter=',')
    Te_fit = np.genfromtxt(args.fits_dir + 'Te.csv', delimiter=',')

    Qp = np.genfromtxt(args.heat_dir + 'Qp.csv', delimiter=',')
    Qp_err = np.genfromtxt(args.heat_dir + 'Qp_err.csv', delimiter=',')
    Qe = np.genfromtxt(args.heat_dir + 'Qe.csv', delimiter=',')
    Qe_err = np.genfromtxt(args.heat_dir + 'Qe_err.csv', delimiter=',')
    Qtot = np.genfromtxt(args.heat_dir + 'Qtot.csv', delimiter=',')
    Qtot_err = np.genfromtxt(args.heat_dir + 'Qtot_err.csv', delimiter=',')
    c_r = np.genfromtxt(args.heat_dir + 'c_r.csv', delimiter=',')
    c_Qtot = np.genfromtxt(args.heat_dir + 'c_Qtot.csv', delimiter=',')

    if args.plot_all:
        plot_fits(args, r, q, q_fit, n_e, n_e_fit, Tp, Tp_fit, Te, Te_fit)
        plot_heat(args, r[1:], Qp, Qp_err, Qe, Qe_err)
        plot_heat_ratio(args, r[1:], Qp, Qp_err, Qe, Qe_err)
        plot_total_heat(args, r[1:], c_r, Qtot, Qtot_err, c_Qtot)
    else:
        if args.plot_fits:
            plot_fits(args, r, q, q_fit, n_e, n_e_fit, Tp, Tp_fit, Te, Te_fit)
        if args.plot_heat:
            plot_heat(args, r[1:], Qp, Qp_err, Qe, Qe_err)
        if args.plot_heat_ratio:
            plot_heat_ratio(args, r[1:], Qp, Qp_err, Qe, Qe_err)
        if args.plot_total_heat:
            plot_total_heat(args, r[1:], c_r, Qtot, Qtot_err, c_Qtot)

    if args.plot_fits or args.plot_heat or args.plot_heat_ratio or args.plot_total_heat or args.plot_all:
        print('Done.')
    else:
        print('Did not choose any plots. See \'python plot.py -h\'.')


if __name__ == '__main__':
    main()
