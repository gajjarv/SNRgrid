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

def get_onspec(fb,f_start_on,f_stop_on,f_start_off,f_stop_off): #Pass filterbank object name 
        onfb = fb.grab_data(f_start=f_start_on,f_stop=f_stop_on)
        offfb = fb.grab_data(f_start=f_start_off,f_stop=f_stop_off)

        onspec = np.mean(onfb[1],axis=0)
        offspec = np.mean(offfb[1],axis=0)

        onspec = (savgol_filter(onspec,11,4))
        offspec = (savgol_filter(offspec,11,4))

        onspec = onspec-np.mean(offspec)
        onspec = onspec/np.std(offspec)

        return onspec


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

    fil_pattern = args.fil
    fil_pattern = fil_pattern + "*.fil"
    fil_list = glob.glob(fil_pattern)

    for fil in fil_list:
        print(fil)
        fb = wt(fil)
        snr = SNR(fb,f_start_on,f_stop_on,f_start_off,f_stop_off)
        beams.append(Point(fb.header['src_raj'],fb.header['src_dej'],snr))

    max_snr = max([point.snr for point in beams])
    min_snr = min([point.snr for point in beams])
    half_max_snr = 0.5 * max_snr

    max_snr_fil = fil_list[np.argmax(snrs)]
    min_snr_fil = fil_list[np.argmin(snrs)]
    half_max_snr_fil = fil_list[np.argmin(np.abs(snrs - half_max_snr))]

    onspec_max = get_onspec(wt(max_snr_fil),f_start_on,f_stop_on,f_start_off,f_stop_off)
    onsepc_min = get_onspec(wt(min_snr_fil),f_start_on,f_stop_on,f_start_off,f_stop_off)	
    onspec_hp  = get_onspec(wt(half_max_snr_fil),f_start_on,f_stop_on,f_start_off,f_stop_off)

    # Convert RA and DEC coordinates of the Point objects into arrays
    ra_hourangle = np.array([beam.ra.value for beam in beams])
    dec_deg = np.array([beam.dec.to('degree').value for beam in beams])
    snrs = np.array([beam.snr for beam in beams])
	
    # interpolate data onto a regular grid
    xi = np.linspace(ra_hourangle.min(), ra_hourangle.max(), 100)
    yi = np.linspace(dec_deg.min(), dec_deg.max(), 100)
    zi = griddata((ra_hourangle, dec_deg), snrs, (xi[None,:], yi[:,None]), method='cubic')
	
    # create contour plot and frequency vs power subplot
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    cmap = 'viridis'
    levels = np.linspace(snrs.min(), snrs.max(), 20)
    
    # plot RA vs DEC contour
    contour_set = axs[0].contour(xi, yi, zi, levels=levels, cmap=cmap)
    axs[0].scatter(ra_hourangle, dec_deg)
    cbar = plt.colorbar(contour_set, ax=axs[0])
    cbar.ax.set_ylabel('SNR')
    axs[0].set_xlabel('RA [h]')
    axs[0].set_ylabel('DEC [deg]')
    
    # plot Frequency vs Power
    axs[1].plot(onspec_max)
    axs[1].plot(onspec_min)
    axs[1].plot(onspec_hp)
    axs[1].set_xlabel('Frequency [MHz]')
    axs[1].set_ylabel('Power')
    
    plt.savefig('Grid_pattern_RA_DEC_Freq_Power.png')

