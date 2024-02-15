"""
other moorings in the Amundsen from Pierre Dutrieux
> extract data, compute monthly means and save in separate files
> save hourly data

Last modified: 8 Feb 2022
"""
import numpy as np
import xarray as xr
import sys
#-------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------
voldir = '/Volumes/SamT5/PhD/data/'
mdir = voldir + 'moorings/'
mdir_new = mdir + 'Currentmeters_PD_monthly/'

localdir = '/Users/ocd1n16/PhD_local/'

auxscriptdir = localdir + 'scripts/aux_func/'
sys.path.append(auxscriptdir)
import aux_func_trend as fc
import aux_matlab_func as matf

#-----------------------------------------------------------------
#-------------------------------------------------------------------
#            ~ ~ ~     MOORING DATA      ~ ~ ~   
#------------------------------------------------------------------
dm = matf.loadmat(mdir + 'Currentmeters_PierreDutrieux.mat')

# currentmeters available
print(dm['Currentmeters'].keys())
cmeters = dm['Currentmeters']

#------------------------------------------------------------------
def extract_time(ndays_ts):
    """
    Convert from matlab format number of days to xarray data format.
    """
    attrs = {'units': 'days since 0001-01-01 00:00:00.0'}
    xr_ndays = xr.Dataset({"time": ("time", ndays_ts, attrs)})
    xr_time = xr.decode_cf(xr_ndays)
    return xr_time

def extract_dataset(moor_name, instr_name, moor_time):
    """
    Extract u and v from a mooring and store them in an xarray dataset
    Store lat, lon, depth and time.
    """
    # some moorings have variable depth
    # ideally can do a more complex analysis bu I'll just take the average here for simplicity
    
    depth_moor = np.asarray(cmeters[moor_name][instr_name]['depth'])
    if depth_moor.size > 1:
        depth_avg = abs(np.mean(depth_moor))
    else:
        depth_avg = abs(depth_moor)
    print("depth: ", depth_avg)

    ds = xr.Dataset({'u' : ('time', cmeters[moor_name][instr_name]['u']),
        'v' : ('time', cmeters[moor_name][instr_name]['v'])},
        coords={'time' : moor_time.time.values,
        'lat' : cmeters[moor_name]['lat'],
        'lon' : cmeters[moor_name]['lon'],
        'depth' : depth_avg})
    return ds

#------------------------------------------------------------------
# BSR14 (entrance to Dotson trough) depth: 388m
#------------------------------------------------------------------
bsr14_time = extract_time(cmeters['bsr14']['AEM_0033']['date']-365)
bsr14 = extract_dataset('bsr14', 'AEM_0033', bsr14_time)

# >> save file
bsr14.to_netcdf(mdir_new + 'bsr14_raw.nc')

# daily and monthly averages
bsr14_crop = bsr14.sel(time=slice('2009-03-01', '2010-11-01'))
bsr14_monthly = bsr14_crop.resample(time='1MS').mean()

# >> save file
bsr14_monthly.to_netcdf(mdir_new + 'bsr14.nc')
#------------------------------------------------------------------
# bsr12 (north of Dotson)
#------------------------------------------------------------------
bsr12_time0 = extract_time(cmeters['bsr12']['AEM_0036']['date']-365)
bsr12_time1 = extract_time(cmeters['bsr12']['Aquadopp_5390']['date']-365)
# depth 400
bsr12_0 = extract_dataset('bsr12', 'AEM_0036', bsr12_time0)
# depth 553
bsr12_1 = extract_dataset('bsr12', 'Aquadopp_5390', bsr12_time1)

# >> save file
bsr12_0.to_netcdf(mdir_new + 'bsr12_0_raw.nc')
bsr12_1.to_netcdf(mdir_new + 'bsr12_1_raw.nc')

bsr12_0_crop = bsr12_0.sel(time=slice('2009-03-01', '2010-12-01'))
bsr12_0_monthly = bsr12_0_crop.resample(time='1MS').mean()

bsr12_1_crop = bsr12_1.sel(time=slice('2009-03-01','2010-12-01'))
bsr12_1_monthly = bsr12_1_crop.resample(time='1MS').mean()

# >> save file
bsr12_0_monthly.to_netcdf(mdir_new + 'bsr12_0.nc')
bsr12_1_monthly.to_netcdf(mdir_new + 'bsr12_1.nc')
#------------------------------------------------------------------
# bsr13a (north of Dotson) depth: 338m
#------------------------------------------------------------------
bsr13a_time = extract_time(cmeters['bsr13a']['AEM_0033']['date']-365)
bsr13a = extract_dataset('bsr13a', 'AEM_0033', bsr13a_time)

# >> save file
bsr13a.to_netcdf(mdir_new + 'bsr13a_raw.nc')

bsr13a_crop = bsr13a.sel(time=slice('2011-01-01', '2012-02-01'))
bsr13a_monthly = bsr13a_crop.resample(time='1MS').mean()

# >> save file
bsr13a_monthly.to_netcdf(mdir_new + 'bsr13a.nc')
#---------------------------
# istar 1, depth 581
#---------------------------
istar1_time = extract_time(cmeters['istar1']['aqd_9392']['date']-365)
istar1 = extract_dataset('istar1', 'aqd_9392', istar1_time)

# >> save file
istar1.to_netcdf(mdir_new + 'istar1_raw.nc')

istar1_crop = istar1.sel(time=slice('2012-03-01', '2014-02-01'))
istar1_monthly = istar1_crop.resample(time='1MS').mean()

# >> save file
istar1_monthly.to_netcdf(mdir_new + 'istar1.nc')
#------------------------------------------------------------------
# istar 4, depth 485.8
#------------------------------------------------------------------
istar4_time = extract_time(cmeters['istar4']['aqd_9396']['date']-365)
istar4 = extract_dataset('istar4', 'aqd_9396', istar4_time)

# >> save file
istar4.to_netcdf(mdir_new + 'istar4_raw.nc')

# crop from 2012-05 instead of 2012-03
istar4_crop = istar4.sel(time=slice('2012-05-01', '2014-02-01'))
istar4_monthly = istar4_crop.resample(time='1MS').mean()

# >> save file
istar4_monthly.to_netcdf(mdir_new + 'istar4.nc')
#---------------------------
# istar 5
#---------------------------
istar5_time0 = extract_time(cmeters['istar5']['aqd_9375']['date']-365)
istar5_time1 = extract_time(cmeters['istar5']['aqd_9386']['date']-365)
istar5_time2 = extract_time(cmeters['istar5']['aqd_9368']['date']-365)
# depth 558
istar5_0 = extract_dataset('istar5', 'aqd_9375', istar5_time0)
# depth 398
istar5_1 = extract_dataset('istar5', 'aqd_9386', istar5_time1)
# depth 717
istar5_2 = extract_dataset('istar5', 'aqd_9368', istar5_time2)

# >> save file
istar5_0.to_netcdf(mdir_new + 'istar5_0_raw.nc')
istar5_1.to_netcdf(mdir_new + 'istar5_1_raw.nc')
istar5_2.to_netcdf(mdir_new + 'istar5_2_raw.nc')


istar5_0_crop = istar5_0.sel(time=slice('2012-03-01', '2014-02-01'))
istar5_0_monthly = istar5_0_crop.resample(time='1MS').mean()

istar5_1_crop = istar5_1.sel(time=slice('2012-03-01', '2014-02-01'))
istar5_1_monthly = istar5_1_crop.resample(time='1MS').mean()

istar5_2_crop = istar5_2.sel(time=slice('2012-03-01','2014-02-01'))
istar5_2_monthly = istar5_2_crop.resample(time='1MS').mean()

# >> save file
istar5_0_monthly.to_netcdf(mdir_new + 'istar5_0.nc')
istar5_1_monthly.to_netcdf(mdir_new + 'istar5_1.nc')
istar5_2_monthly.to_netcdf(mdir_new + 'istar5_2.nc')
#------------------------------------------------------------------
# trough_e (easternmost mooring) 
#------------------------------------------------------------------
troughE_time0 = extract_time(cmeters['trough_e']['aqd_9375']['date']-365)
troughE_time1 = extract_time(cmeters['trough_e']['aqd_9396']['date']-365)
# depth 467
troughE0 = extract_dataset('trough_e', 'aqd_9375', troughE_time0)
# depth 595
troughE1 = extract_dataset('trough_e', 'aqd_9396', troughE_time1)

# >> save file
troughE0.to_netcdf(mdir_new + 'troughE0_raw.nc')
troughE1.to_netcdf(mdir_new + 'troughE1_raw.nc')

troughE0_crop = troughE0.sel(time=slice('2014-03-01', '2016-01-01'))
troughE0_monthly = troughE0_crop.resample(time='1MS').mean()

troughE1_crop = troughE1.sel(time=slice('2014-03-01', '2016-01-01'))
troughE1_monthly = troughE1_crop.resample(time='1MS').mean()

# >> save file
troughE0_monthly.to_netcdf(mdir_new + 'troughE0.nc')
troughE1_monthly.to_netcdf(mdir_new + 'troughE1.nc')
#------------------------------------------------------------------
# trough_w 
#------------------------------------------------------------------
troughW_time0 = extract_time(cmeters['trough_w']['aqd_9368']['date']-365)
troughW_time1 = extract_time(cmeters['trough_w']['aqd_9386']['date']-365)
# depth 423
troughW0 = extract_dataset('trough_w', 'aqd_9368', troughW_time0)
# depth 555
troughW1 = extract_dataset('trough_w', 'aqd_9386', troughW_time1)

# >> save file
troughW0.to_netcdf(mdir_new + 'troughW0_raw.nc')
troughW1.to_netcdf(mdir_new + 'troughW1_raw.nc')

troughW0_crop = troughW0.sel(time=slice('2014-03-01', '2016-01-01'))
troughW0_monthly = troughW0_crop.resample(time='1MS').mean()

troughW1_crop = troughW1.sel(time=slice('2014-03-01', '2016-01-01'))
troughW1_monthly = troughW1_crop.resample(time='1MS').mean()

# >> save file
troughW0_monthly.to_netcdf(mdir_new + 'troughW0.nc')
troughW1_monthly.to_netcdf(mdir_new + 'troughW1.nc')

print("Done.")
