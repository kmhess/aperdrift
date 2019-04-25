# aperdrift: Schedule a drift scan on a calibrator
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)
__author__ = "Kelley M. Hess"
__date__ = "$24-apr-2019 16:00:00$"
__version__ = "0.2"

import csv
import datetime

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.coordinates import SkyCoord
from astropy.table import Table, unique
from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

from modules.calibrators import *
from modules.telescope_params import westerbork

###################################################################
# Define some specific functions for drift scans

def read_beams():
    beams = Table.read('pattern39+1.txt', format='ascii')
    return beams


# def whatsup(obstime_utc = args.starttime_utc):
#     cals = flux_cal + pol_cal
#     print(cals)
#     return

def do_drift(obstime_utc, drift_time):
    return obstime_utc + datetime.timedelta(minutes=drift_time)


def parse_args():
    # *** THIS IS A WORK IN PROGRESS AND NOT YET IMPLEMENTED: ***
    parser = ArgumentParser(
        description="Make a driftscan schedule for the Apertif imaging surveys.",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('-d', '--drifts_per_beam', default=1,
                        help='Specify the number of drifts per beam line. Should be an odd integer. (default: %(default)s).')
    parser.add_argument('-s', "--starttime_utc", default="2019-04-24 16:00:00", #"2019-03-16 09:20:00",
                        help="The start time in ** UTC ** ! - format 'YYYY-MM-DD HH:MM:SS' (default: '%(default)s').",
                        type=datetime.datetime.fromisoformat)
    parser.add_argument('-o', '--output', default='temp',
                        help="Specify the root of output csv and png files (default: imaging_sched_%(default)s.csv.)")
    parser.add_argument('-c', '--command_line',
                        help="If option is included, create a .sh file to run atdbspec.",
                        action='store_true')
    parser.add_argument('-v', "--verbose",
                        help="If option is included, print updated UTC times after each scan.",
                        action='store_true')

    args = parser.parse_args()
    return args
###################################################################

def main():

    args = parse_args()

    # User supplied **UTC** start time:
    start_obstime_utc = args.starttime_utc
    current_lst = Time(start_obstime_utc).sidereal_time('apparent', westerbork().lon)

    calib_name = '3C286'
    # calib_name = whatsup(args.starttime_utc)
    drift_cal = SkyCoord.from_name(calib_name)

    print("\n##################################################################")
    print("Calibrator is: {}".format(calib_name))
    print("Calibrator position is: {}".format(drift_cal.to_string('hmsdms')))
    print("\t in degrees: {} {}".format(drift_cal.ra.deg, drift_cal.dec.deg))

    print("Starting LST is :", current_lst)

    wrap = 0 * u.hourangle
    if (current_lst-drift_cal.ra).value > 12.0:
        wrap = 24 * u.hourangle

    print("Starting HA of calibrator is: {}".format(current_lst-drift_cal.ra-wrap))

    beams = read_beams()

    rows = np.array(unique(beams, keys='dDec'))

# Calculate the starting position of the beams and the length of the drift:
    dec_cen = []
    dec_row = []
    ra_start = []
    ha_start = []
    ha_end = []
    drift_time = []
    for i in range(len(rows)):
        # This is for plotting & dec correction purposes!
        dec_row.append(drift_cal.dec.deg + rows[i][1])
        # The center of the PAF should decrease in declination throughout beam calculation.
        dec_cen.append(drift_cal.dec.deg - rows[i][1])
        ha_start.append(np.min(beams[np.where(beams['dDec'] == rows[i][1])]['dHA']))
        ha_end.append(np.max(beams[np.where(beams['dDec'] == rows[i][1])]['dHA']))
        ra_start.append(drift_cal.ra.deg + (ha_start[i] - 0.5) / np.cos(dec_row[i] * u.deg))
        drift_time.append((np.abs(ha_end[i] - ha_start[i] + 1.0) / np.cos((dec_row[i]) * u.deg)) * 12. / 180. * 60. * u.min)
    start_pos = SkyCoord(ra=np.array(ra_start), dec=dec_cen, unit='deg')
    print("Drift time: {}".format([drift_time[i].value for i in range(len(drift_time))]))

# Plot the drift scans:
    fig,ax=plt.subplots(figsize=(9, 6))
    ax.scatter(drift_cal.ra.deg, drift_cal.dec.deg, c='red', s=20)
    ax.scatter(drift_cal.ra.deg + beams['dHA'] / np.cos((drift_cal.dec.deg + beams['dDec']) * u.deg),
               drift_cal.dec.deg + beams['dDec'], s=10, marker='o', facecolor='black')
    ax.scatter(ra_start, dec_row, s=10, marker='*', facecolor='brown')
    for i in range(len(rows)):
        ax.plot([ra_start[i], ra_start[i] + (ha_end[i] - ha_start[i] + 1.0) / np.cos((dec_row[i]) * u.deg)],
                 [dec_cen[i], dec_cen[i]])
    xlims = ax.get_xlim()
    plt.xlim(xlims[1]+0.25,xlims[0]-0.25)
    # plt.savefig('3c147_driftscan.png'.format(datetime.datetime.now().strftime("%m%d-%H%M%S")))
    plt.savefig(calib_name + '_driftscan.png')

    print('\n Total time (not counting slew): {}'.format(sum(drift_time) + 2. * len(drift_time) * u.min))

    # Open & prepare CSV file to write parset parameters to, in format given by V.M. Moss.
    # Don't worry about slew time because 2 minute wait will always be longer.
    with open(calib_name + "_drift" + args.starttime_utc.strftime("%Y%m%d") + ".csv", "w") as csvfile:
        csvfile.write('source,ha,dec,date1,time1,date,time2,int,type,weight,beam,switch_type\n')
        for i in range(len(rows)):
            current_lst = Time(start_obstime_utc).sidereal_time('apparent', westerbork().lon)
            wrap = 0 * u.hourangle
            if (current_lst - drift_cal.ra).value > 12.0:
                wrap = 24 * u.hourangle
            telescope_position_hadec = current_lst - start_pos[i].ra - wrap
            end_obstime_utc = do_drift(start_obstime_utc, drift_time[i].value)
            date1, time1 = start_obstime_utc.strftime('%Y-%m-%d'), start_obstime_utc.strftime('%H:%M:%S')
            date2, time2 = end_obstime_utc.strftime('%Y-%m-%d'), end_obstime_utc.strftime('%H:%M:%S')
            csvfile.write('{}_drift,{:.6f},{:.6f},{},{},{},{},10,S*,compound,0,system,300\n'.format(calib_name, telescope_position_hadec.deg,
                                                                                                    dec_cen[i], date1, time1, date2, time2))
            start_obstime_utc = end_obstime_utc + datetime.timedelta(minutes=2.0)

    # Open & prepare SH file to run atdbspec, in format given by V.M. Moss.
    if args.command_line:
        print(" *** NOTE THE sh FORMAT FOR ATDBSPEC MIGHT BE DEPRECATED! *** ")
        start_obstime_utc = args.starttime_utc

        with open(calib_name + "_drift" + args.starttime_utc.strftime("%Y%m%d") + ".sh", "w") as file:
            for i in range(len(rows)):
                current_lst = Time(start_obstime_utc).sidereal_time('apparent', westerbork().lon)
                wrap = 0 * u.hourangle
                if (current_lst - drift_cal.ra).value > 12.0:
                    wrap = 24 * u.hourangle
                telescope_position_hadec = current_lst - start_pos[i].ra - wrap
                end_obstime_utc = do_drift(start_obstime_utc, drift_time[i].value)

                print("atdb_service --field_name=3C48_drift --field_ha={:.6f} --field_dec={:.6f} --starttime='{}' "
                      "--endtime='{}' --parset_only --parset_location=/home/apertif/hess/parset_start_observation_driftscan_atdb.template "
                      "--pattern=square_39p1 --observing_mode=imaging --integration_factor=10 --telescopes=2345679ABCD "
                      "--central_frequency=1400 --data_dir=/data/apertif/ --operation=specification --atdb_host=prod".
                      format(telescope_position_hadec.deg, dec_cen[i],
                             start_obstime_utc.strftime("%Y-%m-%d %H:%M:%S"), end_obstime_utc.strftime("%Y-%m-%d %H:%M:%S")))
                file.write(
                    "atdb_service --field_name=3C48_drift --field_ha={:.6f} --field_dec={:.6f} --starttime='{}' "
                    "--endtime='{}' --parset_only --parset_location=/home/apertif/hess/parset_start_observation_driftscan_atdb.template "
                    "--pattern=square_39p1 --observing_mode=imaging --integration_factor=10 --telescopes=2345679ABCD "
                    "--central_frequency=1400 --data_dir=/data/apertif/ --operation=specification --atdb_host=prod\n".
                    format(telescope_position_hadec.deg, dec_cen[i],
                           start_obstime_utc.strftime("%Y-%m-%d %H:%M:%S"), end_obstime_utc.strftime("%Y-%m-%d %H:%M:%S")))

                start_obstime_utc = end_obstime_utc + datetime.timedelta(minutes=2.0)

    # Don't add the last 2 minutes to write out the end time because we don't care about the delay for the data writer.
    print("\nEnding observations! UTC: " + str(end_obstime_utc))

    print("CSV file written to {}".format(csvfile.name))
    if args.command_line:
        print("ATDB commands written to {}".format(file.name))
    print("##################################################################\n")


if __name__ == '__main__':
    main()
