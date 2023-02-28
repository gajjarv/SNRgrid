#!/usr/bin/env python
# To plot SNR vs RA vs DEC
# - Vishal (27 feb 23)
from SNRline import SNR
from blimpy import Waterfall as wt
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from astropy.coordinates import SkyCoord
import astropy.units as u
import random
from astropy.coordinates import Angle
from scipy.interpolate import griddata
import glob
import argparse


class Point:
    def __init__(self, ra, dec, snr):
        self.ra = ra
        self.dec = dec
        self.snr = snr

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some arguments.')
    parser.add_argument('fil', type=str, help='Path to the .fil file')
    parser.add_argument('--f_start_on', type=float, default=6668.37, help='f_start for on-line')
    parser.add_argument('--f_stop_on', type=float, default=6668.40, help='f_stop for on-line')
    parser.add_argument('--f_start_off', type=float, default=6668.65, help='f_start for off-line')
    parser.add_argument('--f_stop_off', type=float, default=6668.75, help='f_stop for off-line')
    args = parser.parse_args()

    #On Line
    f_start_on = args.f_start_on
    f_stop_on = args.f_stop_on

    #Off line
    f_start_off = args.f_start_off
    f_stop_off = args.f_stop_off
    beams = []


    ''' 
    fil_pattern = args.fil
    fil_list = glob.glob(fil_pattern)

    for fil in fil_list:
    	print(fil)
	fb = wt(fil)
	snr = SNR(fb,f_start_on,f_stop_on,f_start_off,f_stop_off)
	beams.append(Point(fb.header['src_raj']+raoffset,fb.header['src_dej']+decoffset,snr))
    '''

    fb = wt(args.fil)
    snr = SNR(fb,f_start_on,f_stop_on,f_start_off,f_stop_off)
    raoffset = 0.0*u.deg
    snr1 = snr
    for i in range(10):
        decoffset = 0.0*u.deg
        for j in range(10):
            beams.append(Point(fb.header['src_raj']+raoffset,fb.header['src_dej']+decoffset,snr1))
            decoffset+=0.5*u.deg
            snr1=snr+0.1*snr*random.uniform(-1,1)
        raoffset+=0.5*u.deg

    # Convert RA and DEC coordinates of the Point objects into arrays
    ra_hourangle = np.array([beam.ra.value for beam in beams])
    dec_deg = np.array([beam.dec.to('degree').value for beam in beams])
    snrs = np.array([beam.snr for beam in beams])
    
    # interpolate data onto a regular grid
    xi = np.linspace(ra_hourangle.min(), ra_hourangle.max(), 100)
    yi = np.linspace(dec_deg.min(), dec_deg.max(), 100)
    zi = griddata((ra_hourangle, dec_deg), snrs, (xi[None,:], yi[:,None]), method='cubic')
    
    # create contour plot
    fig, ax = plt.subplots(figsize=(10, 10))
    cmap = 'viridis'
    levels = np.linspace(snrs.min(), snrs.max(), 20)
    contour_set = ax.contour(xi, yi, zi, levels=levels, cmap=cmap)
    ax.scatter(ra_hourangle,dec_deg)   
 
    # create colorbar
    cbar = plt.colorbar(contour_set, ax=ax)
    cbar.ax.set_ylabel('SNR')
    
    # set axis labels
    ax.set_xlabel('RA [h]')
    ax.set_ylabel('DEC [deg]')
    
    plt.show()

	

