# drift_scan_auto_corr: visualise the auto correlations of the drift scans and
# extract auto correlation data from measurement sets.
# Helga Denes 05/06/2019
# Edited by K.M.Hess (hess@astro.rug.nl)
__author__ = "Helga Denes"
__date__ = "$05-jun-2019 16:00:00$"
__version__ = "0.2"

from argparse import ArgumentParser, RawTextHelpFormatter
import casacore.tables as pt
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd


def plot_row(df, row, pol, args):
    fig, axes = plt.subplots(1, len(row), figsize=(30, 3), subplot_kw={'xticks': [], 'yticks': []})
    fig.subplots_adjust(hspace=0.3, wspace=0.2)

    for ax, beam in zip(axes.flat, row):
        for i in range(12):
            ax.plot(df['time'] - df['time'][0], df['auto_corr_beam_{}_{}_antenna_{}'.format(beam, pol, i)] /
                    np.max(df['auto_corr_beam_{}_{}_antenna_{}'.format(beam, pol, i)]), '-', label='ANT'+str(i))
        ax.plot(df['time'] - df['time'][0],  df['auto_corr_beam_{}_{}'.format(beam, pol)] /
                np.max(df['auto_corr_beam_{}_{}'.format(beam, pol)]), '.-', mec='k', mfc='k', label='average')
        ax.text(50, 0.995, str(beam))

        ax.grid(True, alpha=0.3)
        # ax.set_ylim(600,2800)
        ax.set_ylabel(str(pol)+' Amplitude')
        ax.set_xlabel('Time')

    fig.savefig(str(args.output_path)+str(args.task_id)+'_row_from_'+str(row[0])+'_'+str(pol)+'.png', bbox_inches='tight')


def extract_data(antennas, beams, channel_range, exclude, args):
    df = pd.DataFrame()

    for j in beams:
        t = pt.taql(
            'select TIME, gmeans(means(abs(DATA[' + channel_range + ',0]),0)) as XXPOL, gmeans(means(abs(DATA[5000:20000,3]),0)) as YYPOL from ' + str(
                args.data_location) + str(args.task_id) + '/WSRTA{}_B0{:02}.MS where ANTENNA1==ANTENNA2 {} GROUP BY TIME'.format(args.task_id, j, exclude))
        # t = pt.taql('select TIME, gmeans(means(abs(DATA[{},0]),0)) as XXPOL, gmeans(means(abs(DATA[5000:20000,3]),0))' +
        #             ' as YYPOL from {}{}/WSRTA{}_B0{:02}.MS where ANTENNA1==ANTENNA2 {}' +
        #             ' GROUP BY TIME'.format(channel_range, args.data_location, args.task_id, args.task_id, j, exclude))
        times = np.array(pt.tablecolumn(t, 'TIME'))
        data_beam_xx = np.array(pt.tablecolumn(t, 'XXPOL'))
        auto_corr_xx = np.reshape(data_beam_xx, len(data_beam_xx))
        data_beam_yy = np.array(pt.tablecolumn(t, 'YYPOL'))
        auto_corr_yy = np.reshape(data_beam_yy, len(data_beam_yy))
        df['time'] = times
        df['auto_corr_beam_' + str(j) + '_xx'] = auto_corr_xx
        df['auto_corr_beam_' + str(j) + '_yy'] = auto_corr_yy
        # df = pd.DataFrame.from_dict({'time':times, 'auto_corr_beam_'+str(j):auto_corr})

        for i in antennas:
            t_ant = pt.taql(
                'select TIME, means(abs(DATA[' + channel_range + ',0]),0) as XXPOL, means(abs(DATA[5000:20000,3]),0) as YYPOL from ' + str(
                    args.data_location) + str(args.task_id) + '/WSRTA{}_B0{:02}.MS where ANTENNA1==ANTENNA2 AND ANTENNA1={} GROUP BY TIME'.format(args.task_id, j, i))
            # t_ant = pt.taql(
            #     'select TIME, means(abs(DATA[{},0]),0) as XXPOL, means(abs(DATA[5000:20000,3]),0) as' +
            #     ' YYPOL from {}{}/WSRTA{}_B0{:02}.MS where ANTENNA1==ANTENNA2 AND ANTENNA1={}' +
            #     ' GROUP BY TIME'.format(channel_range, args.data_location, args.task_id, args.task_id, j, i))
            data_ant_xx = np.array(pt.tablecolumn(t_ant, 'XXPOL'))
            data_ant_yy = np.array(pt.tablecolumn(t_ant, 'YYPOL'))
            auto_corr_xx_ant = np.reshape(data_ant_xx, len(data_ant_xx))
            auto_corr_yy_ant = np.reshape(data_ant_yy, len(data_ant_yy))

            df['auto_corr_beam_' + str(j) + '_xx_antenna_' + str(i)] = auto_corr_xx_ant
            df['auto_corr_beam_' + str(j) + '_yy_antenna_' + str(i)] = auto_corr_yy_ant

    return df


def parse_args():

    parser = ArgumentParser(
        description="Visualize and extract auto correlation data from measurment sets.\n"
                    "Call as: > python ../drift_scan_auto_corr.py -t 190531210 -d /data/hess/apertif/CygA_190531/",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('-t', "--task_id", default="190531207",
                        help="The measurement set taskid to extract. (default: '%(default)s').")
    parser.add_argument('-d', '--data_location', default='/Users/hess/apertif/scheduling/aperdrift/modeling/CygA_190531/',
                        help="Specify the root directory. \n(default: '%(default)s').")
    parser.add_argument('-o', '--output_path', default='',
                        help="NOT IMPLEMENTED. Specify the output directory. \n(default: Same as --data_location).")
    parser.add_argument('-a', "--antennas", default="[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]",
                        help="NOT IMPLEMENTED. Array of antennas to include. (default: '%(default)s').")
    # parser.add_argument('-v', "--verbose",
    #                     help="If option is included, print time estimate for several drift combos.",
    #                     action='store_true')

    args = parser.parse_args()
    return args


def main():

    args = parse_args()
    output_path = args.data_location + args.task_id +'/'

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # task_id = '190531208'
    # data_location = '/data/hess/apertif/CygA_190531/'
    # output_path = '/data/hess/apertif/CygA_190531/' + str(args.task_id) + '/'

    # First read HA & Dec from MS file header and save into a csv file for later use.
    # These are the best calculated RA/Dec positions of the centers of the beams.

    beams = range(0, 40)

    hadec = pd.DataFrame()
    ha = np.zeros(40)
    dec = np.zeros(40)
    for j in beams:
        field_table = pt.table('{}{}/WSRTA{}_B0{:02}.MS::FIELD'.format(args.data_location, args.task_id,
                                                                       args.task_id, j))
        (ha[j], dec[j]) = field_table.getcol("DELAY_DIR")[0, 0]
    hadec['ha'] = ha
    hadec['dec'] = dec

    hadec.to_csv(str(output_path) + str(args.task_id) + '_hadec.csv')

    # Extract data and export data into a csv file

    antennas = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    exclude = ''  # to exclude bad antennas
    # exclude = 'AND ANTENNA1!= 11 AND ANTENNA2!=11 AND ANTENNA1!= 10 AND ANTENNA2!=10'
    # channel_range = ['5000:20000'] #RFI free channel range (more or less)
    channel_range = []
    for i in range(5000, 20000, 1500):
        channel_range.append(str(i) + ':' + str(i + 1500))

    for i in range(len(channel_range)):
        print('Extracting data for bin: ', i)
        df_1 = extract_data(antennas, beams, channel_range[i], exclude, args)

        df_1.to_csv(str(output_path) + str(args.task_id) + '_' + str(i) + '_exported_data.csv')

    # # to get antenna names from a measurment set
    # taql_antnames = "SELECT NAME FROM {0}::ANTENNA".format(
    #     str(data_location) + str(task_id) + '/WSRTA' + str(task_id) + '_B000.MS')
    # t = pt.taql(taql_antnames)
    # ant_names = t.getcol("NAME")


if __name__ == '__main__':
    main()
