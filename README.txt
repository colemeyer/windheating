This package allows the user to reproduce the plots from Bandyopadhyay et al. 2023 using data from Parker Solar Probe (PSP), Helios, and Ulysses spacecraft between heliocentric distances of 0.063 and 5.44 AU. The specific plots are the following:

- Figure 1: In situ measurements for the fast solar wind: (a) plasma temperature, (b) electron heat conduction flux,
and (c) electron density. Binned PSP (stars), Helios (squares), and Ulysses (triangles) data are shown in all plots,
with protons in red and electrons in blue.
- Figure 2: Empirically derived heating rates vs. heliocentric distance. Also shown are the error envelopes described
in the text.
- Figure 3: Total heating rates (proton plus electron) vs. heliocentric distance. Error envelopes are described in text.
- Figure 4: Proton heating to total heating ratio vs. heliocentric distance. Error envelopes are described in text.

The package includes 'code' and 'filter'. The 'filter' directory consists of all the velocity-filtered PSP, Helios, and Ulysses data needed to reproduce the plots described above. Scripts in 'code' should be run in the following order when deriving and plotting heating rates for the first time, given their dependencies:

bin.py –> fit.py –> deriv.py –> heat.py –> plot.py