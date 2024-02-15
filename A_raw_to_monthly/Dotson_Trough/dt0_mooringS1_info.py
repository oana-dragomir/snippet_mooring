"""
Mooring S1 (check readme from Tiago)

!! permission to use the data from Anna Wahlin (Uni. of Gothenburg)

DATA
> there are 3 files for this mooring: S110, S111, S112 (2010-2014)
> in mat: u,v, matbv (metres above bottom) - hourly data
> one file for the period 2014-2016

Tiago's PROCESSING:
>> regrid vertically (300-550, every 10m)
>> apply a 40-h Hanning filter to remove tides
>> compute monthly averages

Last modified: 7 Jan 2020

"""

# Import modules
from netCDF4 import num2date
import pandas as pd
import xarray as xr
import sys

#-------------------------------------------------------------------
workdir = '/Volumes/SamT5/PhD/data/'
mdir = workdir + 'moorings/SWE_Dotson/'

localdir = '/Users/ocd1n16/PhD_git/'
auxscriptdir = localdir + 'aux_func/'
sys.path.append(auxscriptdir)
import aux_matlab_func as matf
#-------------------------------------------------------------------

# time reference
time_units='days since 0001-01-01 00:00:00.0'


# S1    
fnames = ['S110', 'S111', 'S112', 'S211', 'S212', 'S412']
files = [n + '.mat' for n in fnames]

k = int(input("which file ("0":S110, "1":S111, "2":S112, "3":S211, "4":S212, "5":S412): "))
dm = matf.loadmat(mdir + files[k])
sm = dm[fnames[k]]

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~    
print(">> MOORING: ", sm['name'])
print(sm.keys(), '\n')
mat = sm['mat']
print("hourly-averaged data: ", sm['mat'].keys(), '\n')
print("raw (ADCP) data: ", sm['ins'][0].keys(), '\n')
print('lon: %.3f | lat: %.3f' % (sm['lon'], sm['lat']))
print('depth: %.2f' % sm['depth'])
#print(sm['note'])

# matlab references the number of days to 0000-00-00, 
# whereas python to 0001-01-01 00:00.0 
# -> subtract a year (365 days) form matlab time to account for that
date_in = num2date(sm['in']-365, units=time_units, calendar='gregorian')
date_out = num2date(sm['out']-365, units=time_units, calendar='gregorian')
print("date in:", date_in)
print("date out: \n", date_out)
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~    

# S1 2014-2016
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~    
mdir3 = workdir + 'moorings/SWE_Dotson/S114_adcp/'
dm = matf.loadmat(mdir3 + 'open_S1_ADCP_20142016.mat')
print(dm.keys())

time = dm['Time14'][:, 0]-365 # time is on a grid and in matlab datenum format
date = num2date(time, units=time_units, calendar='gregorian')

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~    
