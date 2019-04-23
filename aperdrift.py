# aperdrift: Schedule a drift scan on a calibrator
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)
__author__ = "Kelley M. Hess"
__date__ = "$07-mar-2019 16:00:00$"
__version__ = "0.1"

import datetime

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.coordinates import SkyCoord
from astropy.table import Table, unique
from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

from modules.calc_slewtime import calc_slewtime  # Wants [ra, dec] start/end positions in radians; outputs seconds.
from modules.calibrators import *
from modules.telescope_params import westerbork

###################################################################
# Define some specific functions for drift scans

def read_beams():
    beams = Table.read('pattern39+1.txt', format='ascii')
    return beams


def do_drift(obstime_utc, drift_time):
    return obstime_utc + datetime.timedelta(minutes=np.abs(drift_time))


def parse_args():
    # *** THIS IS A WORK IN PROGRESS AND NOT YET IMPLEMENTED: ***
    parser = ArgumentParser(
        description="Make a driftscan schedule for the Apertif imaging surveys.",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('-d', '--drifts_per_beam', default=1,
                        help='Specify the number of drifts per beam line. Should be an odd integer. (default: %(default)s).')
    parser.add_argument('-s', "--starttime_utc", default="2019-03-25 20:00:00",
                        help="The start time in ** UTC ** ! - format 'YYYY-MM-DD HH:MM:SS' (default: '%(default)s').",
                        type=datetime.datetime.fromisoformat)
    parser.add_argument('-o', '--output', default='temp',
                        help='Specify the root of output csv and png files (default: imaging_sched_%(default)s.csv.)')
    parser.add_argument('-v', "--verbose",
                        help="If option is included, print updated UTC times after each scan.",
                        action='store_true')

    args = parser.parse_args()
    return args
###################################################################

def main():

    # User supplied **UTC** start time:
    starttime_utc = Time('2019-03-16 09:20:00')  # Start observing at 5 pm Local (4 UTC) on 01 Mar 2019.

    obstime_utc = starttime_utc
    current_lst = Time(obstime_utc).sidereal_time('apparent', westerbork().lon)

    # User supplied number of rows:
    nrows = 7

    calib_name = '3C48'
    drift_cal = SkyCoord.from_name(calib_name)

    print("\n##################################################################")
    print("Calibrator position is: {}".format(drift_cal.to_string('hmsdms')))
    print("\t in degrees: {} {}".format(drift_cal.ra.deg, drift_cal.dec.deg))

    print("Starting LST is :", current_lst)

    wrap = 0 * u.hourangle
    if (current_lst-drift_cal.ra).value > 12.0:
        wrap = 24*u.hourangle

    print("Starting HA of calibrator is: {}".format(current_lst-drift_cal.ra-wrap))

    beams = read_beams()

    rows = np.array(unique(beams, keys='dDec'))

    dec_cen = []
    dec_row = []
    ra_start = []
    ha_start = []
    ha_end = []
    drift_time = []
    for i in range(len(rows)):
        dec_cen.append(drift_cal.dec.deg - rows[i][1])
        dec_row.append(drift_cal.dec.deg + rows[i][1])
        ha_start.append(np.max(beams[np.where(beams['dDec'] == rows[i][1])]['dHA']))
        ha_end.append(np.min(beams[np.where(beams['dDec'] == rows[i][1])]['dHA']))
        ra_start.append(drift_cal.ra.deg + (ha_start[i] + 0.5) / np.cos(dec_row[i] * u.deg))
        drift_time.append(((ha_end[i] - ha_start[i] - 1.0) / np.cos((dec_row[i]) * u.deg)) * 12. / 180. * 60. * u.min)
    start_pos = SkyCoord(ra=np.array(ra_start), dec=dec_row, unit='deg')
    print("Drift time: {}".format(drift_time))

    plt.figure(figsize=(10, 6))
    for i in range(len(rows)):
        plt.scatter(drift_cal.ra.deg, drift_cal.dec.deg, c='red', s=20)
        plt.scatter(drift_cal.ra.deg + beams['dHA'] / np.cos((drift_cal.dec.deg + beams['dDec']) * u.deg),
                    drift_cal.dec.deg + beams['dDec'], s=10, marker='o', facecolor='black')
        plt.scatter(ra_start, dec_row, s=10, marker='*', facecolor='brown')
        plt.plot([ra_start[i], ra_start[i] + (ha_end[i] - ha_start[i] - 1.0) / np.cos((dec_row[i]) * u.deg)],
                 [dec_row[i], dec_row[i]])
        plt.xlim(np.ceil(np.min(ra_start)+1.), np.ceil(np.min(ra_start)-4.))
    # plt.savefig('3c147_driftscan.png'.format(datetime.datetime.now().strftime("%m%d-%H%M%S")))
    plt.savefig('3c48_driftscan.png')

    print('\n Total time (not counting slew): {}'.format(sum(drift_time) + 2. * len(drift_time) * u.min))

    telescope_position = start_pos[0]
    telescope_position = SkyCoord('15h09m38s 60d41m51s')

    with open("3c48_drift" + starttime_utc.strftime("%Y%m%d") + ".sh", "w") as file:
        for i in range(len(rows)):
            slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
                                        [start_pos[i].ra.radian, start_pos[i].dec.radian])

            start_obstime_utc = obstime_utc + datetime.timedelta(seconds=slew_seconds)
            current_lst = Time(start_obstime_utc).sidereal_time('apparent', westerbork().lon)
            telescope_position_hadec = current_lst - start_pos[i].ra - wrap
            end_obstime_utc = do_drift(start_obstime_utc, drift_time[i].value)

            print("atdb_service --field_name=3C48_drift --field_ha={:.6f} --field_dec={:.6f} --starttime='{}' "
                  "--endtime='{}' --parset_only --parset_location=/home/apertif/hess/parset_start_observation_driftscan_atdb.template "
                  "--pattern=square_39p1 --observing_mode=imaging --integration_factor=10 --telescopes=2345679ABCD "
                  "--central_frequency=1400 --data_dir=/data/apertif/ --operation=specification --atdb_host=prod".
                  format(telescope_position_hadec.deg, start_pos[i].dec.deg,
                         start_obstime_utc.strftime("%Y-%m-%d %H:%M:%S"), end_obstime_utc.strftime("%Y-%m-%d %H:%M:%S")))
            file.write(
                "atdb_service --field_name=3C48_drift --field_ha={:.6f} --field_dec={:.6f} --starttime='{}' "
                "--endtime='{}' --parset_only --parset_location=/home/apertif/hess/parset_start_observation_driftscan_atdb.template "
                "--pattern=square_39p1 --observing_mode=imaging --integration_factor=10 --telescopes=2345679ABCD "
                "--central_frequency=1400 --data_dir=/data/apertif/ --operation=specification --atdb_host=prod\n".
                format(telescope_position_hadec.deg, start_pos[i].dec.deg,
                       start_obstime_utc.strftime("%Y-%m-%d %H:%M:%S"), end_obstime_utc.strftime("%Y-%m-%d %H:%M:%S")))

            obstime_utc = end_obstime_utc + datetime.timedelta(minutes=2.0)
            telescope_position = start_pos[i]

    # Subtract the last 2 minutes to write out the end time because we don't care about the delay for the data writer.
    print("\nEnding observations! UTC: " + str(obstime_utc - datetime.timedelta(minutes=2.0)))
    print("ATDB commands written to {}".format(file.name))
    print("##################################################################\n")


if __name__ == '__main__':
    main()
