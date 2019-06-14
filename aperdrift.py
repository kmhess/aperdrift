# aperdrift: Schedule a drift scan on a calibrator
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)
__author__ = "Kelley M. Hess"
__date__ = "$28-apr-2019 16:00:00$"
__version__ = "0.4.1"

import datetime

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.coordinates import SkyCoord, FK5
from astropy.table import Table, unique
from astropy.time import Time,TimeDecimalYear
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

from modules.calibrators import *
from modules.telescope_params import westerbork

###################################################################
# Define some specific functions for drift scans

def read_beams():
    beams = Table.read('ancillary/pattern39+1.txt', format='ascii')
    return beams


def do_drift(obstime_utc, drift_time):
    return obstime_utc + datetime.timedelta(minutes=drift_time)


def parse_args():

    parser = ArgumentParser(
        description="Make a driftscan schedule for the Apertif imaging surveys.",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('-d', '--drifts_per_beam', default=1,
                        help='Specify the number of drifts per beam line (default: %(default)s).\n 1 & even numbers drift through bm00.')
    parser.add_argument('-c', '--calib', default='3C147',
                        help='Specify the calibrator. (default: %(default)s).')
    parser.add_argument('-s', "--starttime_utc", default="2019-04-25 16:00:00", #"2019-04-25 07:35:00",
                        help="The start time in ** UTC ** ! - format 'YYYY-MM-DD HH:MM:SS' (default: '%(default)s').",
                        type=datetime.datetime.fromisoformat)
    parser.add_argument('-o', '--output', default='_temp',
                        help="Specify the root of output csv and png files (default: CALname_driftDATE%(default)s.csv.)")
    parser.add_argument('-v', "--verbose",
                        help="If option is included, print time estimate for several drift combos.",
                        action='store_true')

    args = parser.parse_args()
    return args
###################################################################

def main():

    args = parse_args()

    # User supplied **UTC** start time:
    start_obstime_utc = args.starttime_utc
    current_lst = Time(start_obstime_utc).sidereal_time('apparent', westerbork().lon)

    # User supplied calibrator:
    calib_name = args.calib
    drift_cal = SkyCoord.from_name(calib_name)

    # Convert to apparent coordinates for setting up the observations
    test = drift_cal.transform_to('fk5')
    drift_cal = test.transform_to(FK5(equinox='J'+str(Time(start_obstime_utc).decimalyear)))

    print("\n##################################################################")
    print("Calibrator is: {}".format(calib_name))
    print("Calibrator position is: {} in J{:7.2f} coordinates".format(drift_cal.to_string('hmsdms'), Time(start_obstime_utc).decimalyear))
    print("\t in degrees: {} {}".format(drift_cal.ra.deg, drift_cal.dec.deg))

    print("Starting LST is :", current_lst)

    wrap = 0 * u.hourangle
    if (current_lst-drift_cal.ra).value > 12.0:
        wrap = 24 * u.hourangle
    if (current_lst-drift_cal.ra).value < -12.0:
        wrap = -24 * u.hourangle

    print("Starting HA of calibrator is: {}".format(current_lst-drift_cal.ra-wrap))
    if np.abs((current_lst-drift_cal.ra-wrap).hourangle) > 5.:
        print("CALIBRATOR IS NOT ACTUALLY UP (but calculations are still right).")
    print("")

    beams = read_beams()

    rows = np.array(unique(beams, keys='dDec'))

    diff = rows[0][1] - rows[1][1]

    print("Starting observations! UTC: {}".format(args.starttime_utc))

    # Calculate the starting position of the beams and the length of the drift
    # (Do for multiple options which will print if ags.verbose = True):
    for n in range(1,5):
        # Only include explicit drift through beam 0 if only drifting through peak.
        if n > 1:
            r_nozero = np.delete(rows, 3)
        else:
            r_nozero = rows
        dec_cen = []
        dec_row = []
        ra_start = []
        ha_start = []
        ha_end = []
        drift_time = []
        # Add drifts at the start
        if n > 1:
            for i in range(n + 0, 0, -1):  # changed from n-1
                dec_row.append(drift_cal.dec.deg + rows[0][1] + (diff) / (n-1) * i)
                dec_cen.append(drift_cal.dec.deg - rows[0][1] - (diff) / (n-1) * i)
                # Drift starts at high RA, low HA
                ha_start.append(np.min(beams[np.where(beams['dDec'] == rows[0][1])]['dHA']))
                ha_end.append(np.max(beams[np.where(beams['dDec'] == rows[0][1])]['dHA']))
                ra_start.append(drift_cal.ra.deg + (ha_start[-1] - 0.75) / np.cos(dec_row[-1] * u.deg))
                drift_time.append((np.abs(ha_end[-1] - ha_start[-1] + 1.5) / np.cos((dec_row[-1]) * u.deg)) * 12. / 180. * 60. * u.min)
        # Do drifts for all the beams
        for r in r_nozero[:-1]:
            for i in range(n):
                dec_row.append(drift_cal.dec.deg + r[1] - (diff) / n * i)
                dec_cen.append(drift_cal.dec.deg - r[1] + (diff) / n * i)
                ha_start.append(np.min(beams[np.where(beams['dDec'] == r[1])]['dHA']))
                ha_end.append(np.max(beams[np.where(beams['dDec'] == r[1])]['dHA']))
                ra_start.append(drift_cal.ra.deg + (ha_start[-1] - 0.9) / np.cos(dec_row[-1] * u.deg))
                drift_time.append((np.abs(ha_end[-1] - ha_start[-1] + 1.8) / np.cos((dec_row[-1]) * u.deg)) * 12. / 180. * 60. * u.min)
        # Add drifts at the end
        if n > 1:
            for i in range(0, n+1):  # changed from n-1
                dec_row.append(drift_cal.dec.deg + rows[-1][1] - (diff) / (n-1) * i)
                dec_cen.append(drift_cal.dec.deg - rows[-1][1] + (diff) / (n-1) * i)
                # Drift starts at high RA, low HA
                ha_start.append(np.min(beams[np.where(beams['dDec'] == rows[-1][1])]['dHA']))
                ha_end.append(np.max(beams[np.where(beams['dDec'] == rows[-1][1])]['dHA']))
                ra_start.append(drift_cal.ra.deg + (ha_start[-1] - 0.75) / np.cos(dec_row[-1] * u.deg))
                drift_time.append((np.abs(ha_end[-1] - ha_start[-1] + 1.5) / np.cos((dec_row[-1]) * u.deg)) * 12. / 180. * 60. * u.min)

        start_pos = SkyCoord(ra=np.array(ra_start), dec=dec_cen, unit='deg')

        if args.verbose and (int(args.drifts_per_beam) != n):
            print("Drift time: {}".format([drift_time[i].value for i in range(len(drift_time))]))
            print("Number of drifts/beam (approx): {}.".format(n))
            print("Total time: {}.\n".format(np.ceil(sum(drift_time) + 2. * (len(drift_time)-1) * u.min)))

        if (args.verbose and int(args.drifts_per_beam) == n) or (int(args.drifts_per_beam) == n):
            print("Drift time: {}".format([drift_time[i].value for i in range(len(drift_time))]))
            print("Number of drifts/beam (approx): {}.".format(n))
            print("Total time: {}.\n".format(np.ceil(sum(drift_time) + 2. * (len(drift_time)-1) * u.min)))

            # Plot and save the drifts for the one actually requested:
            fig,ax=plt.subplots(figsize=(9, 6))
            ax.scatter(drift_cal.ra.deg, drift_cal.dec.deg, c='red', s=20)
            ax.scatter(drift_cal.ra.deg + beams['dHA'] / np.cos((drift_cal.dec.deg + beams['dDec']) * u.deg),
                   drift_cal.dec.deg + beams['dDec'], s=10, marker='o', facecolor='black')
            ax.scatter(ra_start, dec_row, s=10, marker='*', facecolor='brown')
            for i in range(len(dec_cen)):
                # ax.plot([ra_start[i], ra_start[i] + (ha_end[i] - ha_start[i] + 1.2) / np.cos((dec_row[i]) * u.deg)],
                ax.plot([ra_start[i], ra_start[i] + drift_time[i] / u.min / 60. /12. * 180.],
                        [dec_row[i], dec_row[i]])
            xlims = ax.get_xlim()
            plt.xlim(xlims[1]+0.25,xlims[0]-0.25)
            figfilename = calib_name.replace(' ', '') + '_driftscan.png'
            plt.savefig(figfilename)

            # Open & prepare CSV file to write parset parameters to, in format given by V.M. Moss.
            # Don't worry about slew time because 2 minute wait will always be longer.
            with open(calib_name.replace(' ','') + "_drift" + args.starttime_utc.strftime("%Y%m%d") + args.output + ".csv", "w") as csvfile:
                csvfile.write('source,ra,ha,dec,date1,time1,date2,time2,int,type,weight,beam,switch_type,freqmode,centfreq\n')
                for i in range(len(dec_cen)):
                    sidereal_t = Time(start_obstime_utc).sidereal_time('apparent', westerbork().lon)
                    wrap = 0 * u.hourangle
                    if (sidereal_t - drift_cal.ra).value > 12.0 :
                        wrap = 24 * u.hourangle
                    if (sidereal_t - drift_cal.ra).value < -12.0 :
                        wrap = -24 * u.hourangle
                    telescope_position_hadec = sidereal_t - start_pos[i].ra - wrap
                    end_obstime_utc = do_drift(start_obstime_utc, drift_time[i].value)
                    date1, time1 = start_obstime_utc.strftime('%Y-%m-%d'), start_obstime_utc.strftime('%H:%M:%S')
                    date2, time2 = end_obstime_utc.strftime('%Y-%m-%d'), end_obstime_utc.strftime('%H:%M:%S')
                    offset = (drift_cal.dec.deg - dec_cen[i]) * 60.     # units in arcmins
                    sign = '+' if int(offset) >= 0 else ''
                    csvfile.write('{}drift{}{:02},,{:.6f},{:.6f},{},{},{},{},10,T,compound,0,system,300,1370\n'.format(calib_name.replace(' ',''), sign, int(offset),
                                                                        telescope_position_hadec.deg, dec_cen[i], date1, time1, date2, time2))
                    start_obstime_utc = end_obstime_utc + datetime.timedelta(minutes=2.0)
            print(end_obstime_utc)

    # Don't add the last 2 minutes to write out the end time because we don't care about the delay for the data writer.
    print("Assuming {} drifts, ending observations! UTC: {}".format(args.drifts_per_beam,end_obstime_utc))

    print("CSV file written to {}".format(csvfile.name))
    print("PNG file written to {}".format(figfilename))
    print("##################################################################\n")


if __name__ == '__main__':
    # import timeit
    # print(timeit.timeit("main()", setup="from __main__ import main"))
    main()
