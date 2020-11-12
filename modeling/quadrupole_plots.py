from glob import glob
import os

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np


def make_gifs(root):

    os.system('convert -delay 50 {}*_diff.pdf {}diff_freq.gif'.format(root, root))

    return


def parse_args():

    parser = ArgumentParser(
        description="Make plots of all 40 beams diff xx-yy.",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('-c', '--calibname', default='Cyg A',
                        help="Specify the calibrator. (default: '%(default)s').")
    parser.add_argument('-t', "--taskid", default="190531207",
                        help="The first taskid in the set. (default: '%(default)s').")
    parser.add_argument('-o', '--root', default='/Users/hess/apertif/scheduling/aperdrift/modeling/CygA_190607/',
                        help="Specify the root directory. \n(default: '%(default)s').")

    args = parser.parse_args()
    return args


def main():

    args = parse_args()
    if args.root[-1] != '/':
        args.root += '/'

    # Read in difference files and beam files

    xxyy_diff = glob(args.root + 'CygA_*diff.fits')
    stokesI = glob(args.root + 'CygA_*I.fits')

    # Plot difference with beam contours for all beams at one frequency

    row_0 = ['00']
    rows = [['39', '38', '37', '36', '35', '34', '33'],
            ['32', '31', '30', '29', '28', '27', 'aa'],
            ['26', '25', '24', '23', '22', '21', 'aa'],
            ['20', '19', '18', '17', '16', '15', 'aa'],
            ['14', '13', '12', '11', '10', '09', '08'],
            ['07', '06', '05', '04', '03', '02', '01']]

    fig, axs = plt.subplots(6, 7, figsize=(12, 10))

    for ax, bm in zip(np.array(axs).flatten(), np.array(rows).flatten()):
        if bm != 'aa':
            im = [s for s in xxyy_diff if "_"+bm+"_" in s]
            hdu = fits.open(im[0])
            im2 = [s for s in stokesI if "_"+bm+"_" in s]
            hduc = fits.open(im2[0])
            c1, c2 = hduc[0].header['CRPIX1'], hduc[0].header['CRPIX2']
            cax = ax.imshow(hdu[0].data[6, max(0, int(c2-25)):int(c2+25), max(0, int(c1-25)):int(c1+25)],
                            interpolation='none', vmin=-0.1, vmax=0.1, cmap='viridis')
            ax.contour(hduc[0].data[6, max(0, int(c2-25)):int(c2+25), max(0, int(c1-25)):int(c1+25)],
                       levels=[0.25, 0.5, 0.75], colors='k')
            ax.set_title("Beam {}".format(bm), size=8)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            # FITS pixels indexed at 1!
            ax.scatter(hdu[0].header['CRPIX1']-max(0, int(c1-24)), hdu[0].header['CRPIX2']-max(0, int(c2-24)),
                       marker='x', c='k', s=3)
        else:
            ax.axis('off')
    # cbar = fig.colorbar(cax, ax=axs.ravel().tolist(), shrink=0.95, aspect=25)
    fig.suptitle("Polarization Difference, Chan 6: {}".format(args.root[-7:-1]))
    cbaxes = fig.add_axes([0.915, 0.125, 0.03, 0.745])
    cb = plt.colorbar(cax, cax = cbaxes)

    plt.savefig(args.root + '/CygA_{}_diff.pdf'.format(args.root[-7:-1]))

    fig, axs = plt.subplots(3, 4, figsize=(12, 10))

    beam = "00"
    im = [s for s in xxyy_diff if "_" + beam + "_" in s]
    hdu = fits.open(im[0])
    im2 = [s for s in stokesI if "_" + beam + "_" in s]
    hduc = fits.open(im2[0])
    c1, c2 = hduc[0].header['CRPIX1'], hduc[0].header['CRPIX2']
    for ax, chan in zip(np.array(axs).flatten()[:10], range(10)):
        cax = ax.imshow(hdu[0].data[chan, max(0, int(c2 - 25)):int(c2 + 25), max(0, int(c1 - 25)):int(c1 + 25)],
                        interpolation='none', vmin=-0.1, vmax=0.1, cmap='viridis')
        ax.contour(hduc[0].data[6, max(chan, int(c2 - 25)):int(c2 + 25), max(0, int(c1 - 25)):int(c1 + 25)],
                   levels=[0.25, 0.5, 0.75], colors='k')
        ax.set_title("Channel {}".format(chan), size=8)
        ax.set_xlabel("RA")
        ax.set_ylabel("Declination")
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)
        # FITS pixels indexed at 1!
        ax.scatter(hdu[0].header['CRPIX1'] - max(0, int(c1 - 24)), hdu[0].header['CRPIX2'] - max(0, int(c2 - 24)),
                   marker='x', c='k', s=3)
    for ax in np.array(axs).flatten()[10:12]:
        ax.axis('off')

    fig.suptitle("Polarization Diff as functn of Freq, Beam {}: {}".format(beam, args.root[-7:-1]))
    cbaxes = fig.add_axes([0.915, 0.125, 0.03, 0.745])
    cb = plt.colorbar(cax, cax = cbaxes)

    plt.savefig(args.root + '/CygA_{}_{}_diff_freq.pdf'.format(args.root[-7:-1], beam))


if __name__ == '__main__':
    main()
