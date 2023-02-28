This is a code that takes a number of filterbank files as an input and measure SNR of a spectral line in it. 
usage: SNR_RA_DEC.py [-h] [--f_start_on F_START_ON] [--f_stop_on F_STOP_ON] [--f_start_off F_START_OFF] [--f_stop_off F_STOP_OFF] fil

Process some arguments.

positional arguments:
  fil                   Path to the .fil file

optional arguments:
  -h, --help            show this help message and exit
  --f_start_on F_START_ON
                        f_start for on-line
  --f_stop_on F_STOP_ON
                        f_stop for on-line
  --f_start_off F_START_OFF
                        f_start for off-line
  --f_stop_off F_STOP_OFF
                        f_stop for off-line


It creates a grid of measured SNR as a function of RA and DEC. 
