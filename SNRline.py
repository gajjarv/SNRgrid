#!/usr/bin/env python
# To calculate SNR of spectral line
# - Vishal (23 feb 23)

from blimpy import Waterfall as wt
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


def SNR(fb,f_start_on,f_stop_on,f_start_off,f_stop_off): #Pass blimpy filterbnk object with ON and OFF line
	onfb = fb.grab_data(f_start=f_start_on,f_stop=f_stop_on)
	offfb = fb.grab_data(f_start=f_start_off,f_stop=f_stop_off)	

	onspec = np.mean(onfb[1],axis=0)
	offspec = np.mean(offfb[1],axis=0)

	onspec = (savgol_filter(onspec,11,4))
	offspec = (savgol_filter(offspec,11,4))

	onspec = onspec-np.mean(offspec)
	onspec = onspec/np.std(offspec)

	offspec = offspec-np.mean(offspec)
	offspec = offspec/np.std(offspec)
	
	#plt.plot(onspec)
	#plt.plot(offspec[0:len(onspec)])	
	#plt.show()

	snr = np.sum(onspec)
	return snr

if __name__ == "__main__":
	#On Line
	f_start_on = 6668.37
	f_stop_on = 6668.40
	
	#Off line
	f_start_off = 6668.65
	f_stop_off = 6668.75

	fil = sys.argv[1]
	#print fil
	fb = wt(fil)
	snr = SNR(fb,f_start_on,f_stop_on,f_start_off,f_stop_off)
	#print snr


