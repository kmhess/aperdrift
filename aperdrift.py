# aperdrift: Schedule a drift scan on a calibrator
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)
__author__ = "Kelley M. Hess"
__date__ = "$07-mar-2019 16:00:00$"
__version__ = "0.1"

import csv
import datetime

from astropy.coordinates import Angle, Longitude, SkyCoord
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


def do_drift(obstimeUTC,drift_time):
    return obstimeUTC + datetime.timedelta(minutes=drift_time)

###################################################################

def main():

    # User supplied **UTC** start time:
    starttimeUTC = datetime.datetime(2019, 3, 8, 20, 30, 0)  # Start observing at 5 pm Local (4 UTC) on 01 Mar 2019.
    endtimeUTC = datetime.datetime(2019, 3, 8, 23, 30, 0)  # Start observing at 5 pm Local (4 UTC) on 01 Mar 2019.

    # User supplied filename:
    csv_filename = 'drift_scan_3c147.csv'

    # User supplied number of rows:
    nrows = 7

    calib_name = '3C147'
    drift_cal = SkyCoord.from_name(calib_name)
    print("Calibrator position is: ", drift_cal.to_string('hmsdms'))
    print("\t in degrees: ", drift_cal.ra.deg, drift_cal.dec.deg)

    beams = read_beams()

    rows = np.array(unique(beams, keys='dDec'))

    declination = []
    ra_start = []
    hour_angle_degs = []
    drift_degs = []
    drift_ra = []
    drift_time = []
    for i in range(len(rows)):
        declination.append(drift_cal.dec.deg + rows[i][1])
        hour_angle_degs.append([np.max(beams[np.where(beams['dDec'] == rows[i][1])]['dHA']),
                                np.min(beams[np.where(beams['dDec'] == rows[i][1])]['dHA'])])
        drift_degs.append(1.0 + np.max(beams[np.where(beams['dDec'] == rows[i][1])]['dHA']) -
                          np.min(beams[np.where(beams['dDec'] == rows[i][1])]['dHA']))
        drift_ra.append(drift_degs[i] / np.cos(declination[i] * u.deg))
        drift_time.append(drift_ra[i] * 12. / 180. * 60. * u.min)
        ra_start.append((drift_cal.ra.deg - (-0.0 + (hour_angle_degs[i][0])) / np.cos(declination[i] * u.deg)))
        # print(ra_start[i], declination[i], drift_degs[i], drift_ra[i], drift_time[i])
    start_pos = SkyCoord(ra=np.array(ra_start), dec=declination, unit='deg')

    plt.figure()
    for i in range(len(rows)):
        plt.scatter(drift_cal.ra.deg - beams['dHA'] / np.cos(beams['dDec'] * u.deg), drift_cal.dec.deg + beams['dDec'], s=10, marker='o', facecolor='black')
        #plt.scatter(drift_cal.ra.deg - beams['dHA'] + 1, drift_cal.dec.deg + beams['dDec'], s=10, marker='o', facecolor='brown')
        plt.plot([ra_start[i], ra_start[i] + drift_degs[i]], [declination[i], declination[i]])
    plt.savefig('test.png')

    print('\n Total time (not counting slew): {}'.format(sum(drift_time) + 2. * len(drift_time) * u.min))

    header = ['source', 'ra', 'dec', 'date1', 'time1', 'date2', 'time2', 'int', 'type', 'weight', 'beam', 'switch_type']

    # with open(csv_filename, 'w') as csvfile:
    #     writer = csv.writer(csvfile)
    #     writer.writerow(header)

    telescope_position = start_pos[0]
    obstimeUTC = starttimeUTC
    currentLST = Time(obstimeUTC).sidereal_time('apparent', westerbork().lon)
    print("Starting LST is :", currentLST)

    for i in range(len(rows)):
        slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
                                    [start_pos[i].ra.radian, start_pos[i].dec.radian])

        start_obstimeUTC = obstimeUTC + datetime.timedelta(seconds=slew_seconds)
        currentLST = Time(start_obstimeUTC).sidereal_time('apparent', westerbork().lon)
        telescope_position_hadec = currentLST - start_pos[i].ra
        # print("HA? ",telescope_position_hadec)

        end_obstimeUTC = do_drift(start_obstimeUTC,drift_time[i].value)
            # write_to_csv(csvfile, calib_name, telescope_position, new_obstimeUTC, after_target)
        print("atdb_service --field_name=3C147_drift --field_ha={:.6f} --field_dec={:.6f} --starttime='{}' --endtime='{}' --parset_only --parset_location=/home/apertif/hess/parset_start_observation_atdb.template --pattern=square_39p1 --observing_mode=imaging --integration_factor=10 --telescopes=2345679ABCD --central_frequency=1400 --data_dir=/data/apertif/ --operation=specification --atdb_host=prod".
              format(telescope_position_hadec.deg, start_pos[i].dec.deg, start_obstimeUTC.strftime("%Y-%m-%d %H:%M:%S"), end_obstimeUTC.strftime("%Y-%m-%d %H:%M:%S")))
        obstimeUTC = end_obstimeUTC
        telescope_position = start_pos[i]
         # print("Scan {} observed {}.".format(i, telescope_position))

    print("Ending observations! UTC: " + str(obstimeUTC))
    print("##################################################################\n")


if __name__ == '__main__':
    main()