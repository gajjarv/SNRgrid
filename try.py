import numpy as np
import matplotlib.pyplot as plt

class Point:
    def __init__(self, ra, dec, snr):
        self.ra = ra
        self.dec = dec
        self.snr = snr

# create an array of Point objects
points = [Point(ra, dec, snr) for ra, dec, snr in zip(np.random.rand(100), np.random.rand(100), np.random.rand(100))]

# convert the array of Point objects into a 2D array
ra_values = np.array([point.ra for point in points])
dec_values = np.array([point.dec for point in points])
snr_values = np.array([point.snr for point in points])

# specify the RA and DEC ranges and bin sizes
ra_range = [0, 1] # in degrees
dec_range = [0, 1] # in degrees
ra_bin_size = 0.1 # in degrees
dec_bin_size = 0.1 # in degrees

# convert RA and DEC values to pixel coordinates
ra_bins = np.arange(ra_range[0], ra_range[1]+ra_bin_size, ra_bin_size)
dec_bins = np.arange(dec_range[0], dec_range[1]+dec_bin_size, dec_bin_size)

ra_pix = np.digitize(ra_values, ra_bins) - 1
dec_pix = np.digitize(dec_values, dec_bins) - 1

# create a 2D array of SNR values binned in RA and DEC
counts = np.zeros((len(dec_bins)-1, len(ra_bins)-1))
for i, (ra, dec) in enumerate(zip(ra_pix, dec_pix)):
    counts[dec, ra] += snr_values[i]

# plot the 2D array as an image using imshow
fig, ax = plt.subplots()
im = ax.imshow(counts, extent=[ra_bins[0], ra_bins[-1], dec_bins[0], dec_bins[-1]], origin='lower', cmap='viridis')
cbar = plt.colorbar(im)
ax.set_xlabel('RA (deg)')
ax.set_ylabel('DEC (deg)')
plt.show()

