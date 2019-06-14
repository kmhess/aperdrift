# scan2fits: Create XX & YY beam models from drift scans
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)
__author__ = "Kelley M. Hess"
__date__ = "$04-jun-2019 16:00:00$"
__version__ = "0.2"

from glob import glob
import os

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.coordinates import SkyCoord, FK5
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS, Wcsprm
# import matplotlib.colors as colors
# import matplotlib.pyplot as plt
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
        description="Make cubes of all 40 beams from drift scans.",
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

    cell_size = 105. / 3600.
    freqchunks = 10

    # Put all the output from drift_scan_auto_corr.ipynb in a unique folder per source, per set of drift scans.
    datafiles = glob(args.root + '*_*_exported_data.csv')
    datafiles.sort()
    posfiles = glob(args.root + '*hadec.csv')
    posfiles.sort()

    # Put calibrator into apparent coordinates (because that is what the telescope observes it in.)
    test = calib.transform_to('fk5')
    calibnow = test.transform_to(FK5(equinox='J{}'.format(taskid2equinox(args.taskid))))

    for beam in range(15, 20):
        print(beam, end=' ')

        # Create the vectors which contain all data from all scans for a given beam which has been specified above.
        # cube_xx, cube_y, db_xx, db_yy = [], [], [], []
        for f in range(freqchunks):
            x, y, z_xx, z_yy = [], [], [], []
            for file, pos in zip(datafiles[f::10], posfiles):
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
                z_xx = np.append(z_xx, data['auto_corr_beam_' + str(beam) + '_xx'] - np.median(
                    data['auto_corr_beam_' + str(beam) + '_xx']))
                z_yy = np.append(z_yy, data['auto_corr_beam_' + str(beam) + '_yy'] - np.median(
                    data['auto_corr_beam_' + str(beam) + '_yy']))

            # Create the 2D plane, do a cubic interpolation, and append it to the cube.
            tx = np.arange(min(x), max(x), cell_size)
            ty = np.arange(min(y), max(y), cell_size)
            XI, YI = np.meshgrid(tx, ty)
            gridcubx = interpolate.griddata((x, y), z_xx, (XI, YI), method='cubic')  # median already subtracted
            gridcuby = interpolate.griddata((x, y), z_yy, (XI, YI), method='cubic')

            # Find the reference pixel at the apparent coordinates of the calibrator
            ref_pixy = (calibnow.dec.deg - min(y)) / cell_size
            ref_pixx = (-min(x)) / cell_size
            # ref_pixz =

            # Find the peak of the primary beam to normalize
            norm_xx = np.max(gridcubx[int(ref_pixy)-3:int(ref_pixy)+4, int(ref_pixx)-3:int(ref_pixx)+4])
            norm_yy = np.max(gridcuby[int(ref_pixy) - 3:int(ref_pixy) + 4, int(ref_pixx) - 3:int(ref_pixx) + 4])
            # if beam == 0:
            #     norm0_xx = np.max(gridcubx[int(ref_pixy) - 3:int(ref_pixy) + 4, int(ref_pixx) - 3:int(ref_pixx) + 4])
            #     norm0_yy = np.max(gridcuby[int(ref_pixy) - 3:int(ref_pixy) + 4, int(ref_pixx) - 3:int(ref_pixx) + 4])

            # Create 3D array with proper size for given scan set to save data as a cube
            if f == 0:
                cube_xx = np.zeros((freqchunks, gridcubx.shape[0], gridcubx.shape[1]))
                cube_yy = np.zeros((freqchunks, gridcuby.shape[0], gridcuby.shape[1]))
                db_xx = np.zeros((freqchunks, gridcubx.shape[0], gridcubx.shape[1]))
                db_yy = np.zeros((freqchunks, gridcuby.shape[0], gridcuby.shape[1]))

            cube_xx[f, :, :] = gridcubx/norm_xx
            cube_yy[f, :, :] = gridcuby/norm_yy

            # Convert to decibels
            db_xx[f, :, :] = np.log10(gridcubx/norm_xx) * 10.
            db_yy[f, :, :] = np.log10(gridcuby/norm_yy) * 10.

        wcs = WCS(naxis=3)
        wcs.wcs.cdelt = np.array([-cell_size, cell_size, 12.207e3*1500])
        wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'FREQ']
        wcs.wcs.crval = [calib.ra.to_value(u.deg), calib.dec.to_value(u.deg), 1219.609e6+(12.207e3*(500+1500/2))]
        wcs.wcs.crpix = [ref_pixx, ref_pixy, 1]
        wcs.wcs.specsys = 'TOPOCENT'
        wcs.wcs.restfrq = 1.420405752e+9
        header = wcs.to_header()

        hdux_db = fits.PrimaryHDU(db_xx, header=header)
        hduy_db = fits.PrimaryHDU(db_yy, header=header)
        hdux = fits.PrimaryHDU(cube_xx, header=header)
        hduy = fits.PrimaryHDU(cube_yy, header=header)

        # Save the FITS files
        hdux_db.writeto(args.root + '{}_{}_{:02}xx_db_cube.fits'.format(args.calibname.replace(" ", ""), args.taskid[:-3],
                                                                        beam), overwrite=True)
        hduy_db.writeto(args.root + '{}_{}_{:02}yy_db_cube.fits'.format(args.calibname.replace(" ", ""), args.taskid[:-3],
                                                                        beam), overwrite=True)
        hdux.writeto(args.root + '{}_{}_{:02}xx_cube.fits'.format(args.calibname.replace(" ", ""), args.taskid[:-3],
                                                                  beam), overwrite=True)
        hduy.writeto(args.root + '{}_{}_{:02}yy_cube.fits'.format(args.calibname.replace(" ", ""), args.taskid[:-3],
                                                                  beam), overwrite=True)


if __name__ == '__main__':
    main()
