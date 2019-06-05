# scan2fits: Create XX & YY beam models from drift scans
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)
__author__ = "Kelley M. Hess"
__date__ = "$04-jun-2019 16:00:00$"
__version__ = "0.1"

from glob import glob
import os

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.coordinates import SkyCoord, EarthLocation, FK5
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

from modules.telescope_params import westerbork

def taskid2equinox(taskid):

    # Automatically take the date of the observaitons from the taskid to calculate apparent coordinates of calibrator
    year = 2000 + int(str(taskid)[0:2])
    month = str(taskid)[2:4]
    day = str(taskid)[4:6]
    equinox = Time('{}-{}-{}'.format(year, month, day))

    return equinox.decimalyear


def make_gifs(root):

    os.system('convert -delay 50 {}*db0_reconstructed.png {}all_beams0.gif'.format(root, root))
    os.system('convert -delay 50 {}*_difference.png {}diff_xx-yy.gif'.format(root, root))

    return


def parse_args():

    parser = ArgumentParser(
        description="Make a model of all 40 beams from drift scans.",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('-c', '--calibname', default='Cyg A',
                        help="Specify the calibrator. (default: '%(default)s').")
    parser.add_argument('-t', "--taskid", default="190531207",
                        help="The first taskid in the set. (default: '%(default)s').")
    parser.add_argument('-o', '--root', default='/Users/hess/apertif/scheduling/aperdrift/modeling/CygA_190531/',
                        help="Specify the root directory. \n(default: '%(default)s').")
    parser.add_argument('-g', '--make_gifs',
                        help="(Re)Make gifs of figures? (default is False).",
                        action='store_true')
    # parser.add_argument('-v', "--verbose",
    #                     help="If option is included, print time estimate for several drift combos.",
    #                     action='store_true')

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    np.warnings.filterwarnings('ignore')

    # Find calibrator position
    calib = SkyCoord.from_name(args.calibname)

    # Put all the output from drift_scan_auto_corr.ipynb in a unique folder per source, per set of drift scans.
    datafiles = glob(args.root + '*exported_data.csv')
    datafiles.sort()
    posfiles = glob(args.root + '*hadec.csv')
    posfiles.sort()

    # Put calibrator into apparent coordinates (because that is what the telescope observes it in.)
    test = calib.transform_to('fk5')
    calibnow = test.transform_to(FK5(equinox='J{}'.format(taskid2equinox(args.taskid))))

    corr_im = []
    diff_im = []
    for beam in range(40):
        print(beam, end=' ')
        # Create the vectors which contain all data from all scans for a given beam which has been specified above.
        x, y, z_xx, z_yy = [], [], [], []
        for file, pos in zip(datafiles, posfiles):
            data = Table.read(file, format='csv')
            hadec = Table.read(pos, format='csv')
            hadec_start = SkyCoord(ra=hadec['ha'], dec=hadec['dec'], unit=(u.rad, u.rad))  # From ALTA (same as above)

            time_mjd = Time(data['time'] / (3600 * 24), format='mjd')
            lst = time_mjd.sidereal_time('apparent', westerbork().lon)

            HAcal = lst - calibnow.ra  # in sky coords
            dHAsky = HAcal - hadec_start[beam].ra + (24 * u.hourangle)  # in sky coords in hours
            dHAsky.wrap_at('180d', inplace=True)
            dHAphys = dHAsky * np.cos(hadec_start[beam].dec.deg * u.deg)  # physical offset in hours

            x = np.append(x, dHAphys.deg)
            y = np.append(y, np.full(len(dHAphys.deg), hadec_start[beam].dec.deg))
            z_xx = np.append(z_xx, data['auto_corr_beam_' + str(beam) + '_xx'] / np.median(
                data['auto_corr_beam_' + str(beam) + '_xx']))
            z_yy = np.append(z_yy, data['auto_corr_beam_' + str(beam) + '_yy'] / np.median(
                data['auto_corr_beam_' + str(beam) + '_yy']))

        # # Add a fake drift that goes to zero power at 1 deg above last scan
        # x=np.append(x,dHAphys.deg)
        # y=np.append(y,np.full(len(dHAphys.deg),max(y)+1.0))
        # z_xx=np.append(z_xx,np.full(len(dHAphys.deg),1))
        # z_yy=np.append(z_yy,np.full(len(dHAphys.deg),1))

        # # Add a fake drift that goes to zero power at 1 deg below first scan
        # x=np.append(x,dHAphys.deg)
        # y=np.append(y,np.full(len(dHAphys.deg),min(y)-1.0))
        # z_xx=np.append(z_xx,np.full(len(dHAphys.deg),1))
        # z_yy=np.append(z_yy,np.full(len(dHAphys.deg),1))

        # Create the 2D grid and do a cubic interpolation
        cell_size = 105. / 3600.
        tx = np.arange(min(x), max(x), cell_size)
        ty = np.arange(min(y), max(y), cell_size)
        XI, YI = np.meshgrid(tx, ty)
        gridcubx = interpolate.griddata((x, y), z_xx, (XI, YI), method='cubic')
        zero_gridcubx = gridcubx - np.nanmedian(gridcubx)
        gridcuby = interpolate.griddata((x, y), z_yy, (XI, YI), method='cubic')
        zero_gridcuby = gridcuby - np.nanmedian(gridcuby)

        # Find the reference pixel at the apparent coordinates of the calibrator
        ref_pixy = (calibnow.dec.deg - min(y)) / cell_size
        ref_pixx = (-min(x)) / cell_size

        # Find the peak of the primary beam to normalize
        norm_xx = np.max(gridcubx[int(ref_pixy)-3:int(ref_pixy)+4,int(ref_pixx)-3:int(ref_pixx)+4])
        norm_yy = np.max(gridcuby[int(ref_pixy) - 3:int(ref_pixy) + 4, int(ref_pixx) - 3:int(ref_pixx) + 4])
        zero_norm_xx = np.max(zero_gridcubx[int(ref_pixy)-3:int(ref_pixy)+4,int(ref_pixx)-3:int(ref_pixx)+4])
        zero_norm_yy = np.max(zero_gridcuby[int(ref_pixy) - 3:int(ref_pixy) + 4, int(ref_pixx) - 3:int(ref_pixx) + 4])
        if beam == 0:
            norm0_xx = np.max(gridcubx[int(ref_pixy) - 3:int(ref_pixy) + 4, int(ref_pixx) - 3:int(ref_pixx) + 4])
            norm0_yy = np.max(gridcuby[int(ref_pixy) - 3:int(ref_pixy) + 4, int(ref_pixx) - 3:int(ref_pixx) + 4])

        # Convert to decibels
        db_xx = np.log10(gridcubx/norm_xx) * 10.
        db_yy = np.log10(gridcuby/norm_yy) * 10.
        db0_xx = np.log10(gridcubx/norm0_xx) * 10.
        db0_yy = np.log10(gridcuby/norm0_yy) * 10.

        wcs = WCS(naxis=2)
        wcs.wcs.cdelt = np.array([-cell_size, cell_size])
        wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        wcs.wcs.crval = [calib.ra.to_value(u.deg), calib.dec.to_value(u.deg)]
        wcs.wcs.crpix = [ref_pixx, ref_pixy]
        header = wcs.to_header()

        hdux = fits.PrimaryHDU(db_xx, header=header)
        hduy = fits.PrimaryHDU(db_yy, header=header)
        # hdulx = fits.HDUList([hdux])
        # hduly = fits.HDUList([hduy])

        # Save the FITS files
        hdux.writeto(args.root + '{}_{}_{:02}xx.fits'.format(args.calibname.replace(" ", ""), args.taskid[:-3], beam),
                     overwrite=True)
        hduy.writeto(args.root + '{}_{}_{:02}yy.fits'.format(args.calibname.replace(" ", ""), args.taskid[:-3], beam),
                     overwrite=True)

        fig1 = plt.figure(figsize=(6, 9))
        ax1 = fig1.add_subplot(2, 1, 1, projection=wcs.celestial)
        ax1.grid(lw=1, color='white')
        ax1.set_title("Beam {:02} - XX Correlation - Cubic".format(beam))
        ax1.set_ylabel("Declination [J2000]")
        ax1.set_xlabel("Right Ascension [J2000]")
        # im1 = ax1.imshow(gridcubx, norm=colors.LogNorm(vmin=0.93, vmax=1.2), cmap='magma', animated=True)
        im1 = ax1.imshow(db0_xx, cmap='magma', vmax=-4.5, vmin=-6.3)
        ax2 = fig1.add_subplot(2, 1, 2, projection=wcs.celestial)
        ax2.grid(lw=1, color='white')
        ax2.set_title("Beam {:02} - YY Correlation - Cubic".format(beam))
        ax2.set_ylabel("Declination [J2000]")
        ax2.set_xlabel("Right Ascension [J2000]")
        # im2 = ax2.imshow(gridcuby, norm=colors.LogNorm(vmin=0.93, vmax=1.2), cmap='magma', animated=True)
        im2 = ax2.imshow(db0_yy, cmap='magma', vmax=-4.5, vmin=-6.3)
        corr_im.append([im1, im2])
        plt.savefig(args.root + '{}_{}_{:02}db0_reconstructed.png'.format(args.calibname.replace(" ", ""),
                                                                       args.taskid, beam))
        plt.close('all')

        # Plot the difference between XX and YY for every beam
        diffcub = zero_gridcubx/zero_norm_xx - zero_gridcuby/zero_norm_yy

        fig2 = plt.figure(figsize=(10, 9))
        ax1 = fig2.add_subplot(1, 1, 1, projection=wcs.celestial)
        ax1.grid(lw=1, color='white')
        ax1.set_title("Beam {:02} - Difference (XX$-$YY)".format(beam))
        ax1.set_ylabel("Declination [J2000]")
        ax1.set_xlabel("Right Ascension [J2000]")
        ax1.scatter(ref_pixx, ref_pixy, marker='x', color='black')
        im3 = ax1.imshow(diffcub, vmin=-0.1, vmax=0.1)
        plt.colorbar(im3)
        diff_im.append([im3])
        plt.savefig(args.root + '{}_{}_{:02}_difference.png'.format(args.calibname.replace(" ", ""),
                                                                    args.taskid, beam))
        plt.close('all')

    if args.make_gifs:
        make_gifs(args.root)

if __name__ == '__main__':
    main()
