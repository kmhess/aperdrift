# beam_fit.py
# This program will read in an Apertif FITS file containing one beam map and
# fit it with interp2d.
# 6/13/19 DJP

from glob import glob

from scipy.interpolate import RectBivariateSpline
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from numpy import arange
import pandas as pd


def main():
    files=glob('CygA_190607/CygA_190607_*_I.fits')
    files.sort()

    hdu=fits.open(files[0])
    beam=hdu[0].data
    f0=hdu[0].header['CRVAL3']
    fdelt=hdu[0].header['CDELT3']
    hdu.close()

    chan = 9
    freq = f0 + fdelt * chan
    bmap, fbeam, model = [], [], []
    # for i,t in enumerate(test[:33]):
    for i,t in enumerate(files):
        hdu=fits.open(t)
        beam=hdu[0].data
        hdu.close()
        bmap.append(beam[chan,int(hdu[0].header['CRPIX2'])-20:int(hdu[0].header['CRPIX2'])+20,
                     int(hdu[0].header['CRPIX1'])-20:int(hdu[0].header['CRPIX1'])+20])
        x=arange(0,bmap[i].shape[1])
        y=arange(0,bmap[i].shape[0])
        fbeam.append(RectBivariateSpline(y,x,bmap[i]))
        model.append(fbeam[i](y,x))

    fig, axes = plt.subplots(5, 8, figsize=(15, 12))
    fig.suptitle('Beams for f = {} Hz'.format(freq))
    for i, ax in enumerate(np.array(axes).flatten()):
        #     ax.imshow(bmap[i],origin='lower',vmax=0.35)
        ax.imshow(bmap[i], origin='lower', vmax=300000)
        #     ax.contour(model[i],levels=np.linspace(0.01,1.,1),colors='white')
        ax.contour(model[i], levels=np.linspace(0, 100000, 700000), colors='white')
        ax.set_title('Beam {:02}'.format(i))
    plt.savefig('models_190607.png')

    df = pd.DataFrame()
    maxlen = len(model[0].flatten())
    for b,m in enumerate(model):
        col = list(m.flatten())
        if maxlen != len(col):
            col.extend(['']*(maxlen-len(col)))
        df['B{:02}_{}'.format(b,int(freq))] = col

    df.to_csv('/Users/hess/apertif/scheduling/aperdrift/modeling/models/beam_models_190607_unnorm.csv')


if __name__ == '__main__':
    main()