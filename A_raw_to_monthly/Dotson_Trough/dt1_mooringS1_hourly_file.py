"""
Extract velocity from ADCP and T,S,P from microcats
for S1 mooring.

Resample data hourly - although it is already at that resolution but 
resampling removes duplicates when merging

[skip the Hanning filter because monthly averaging removes the effect of tides] 
Filter tides (but do some interpolation beforehand to have less nans)

Save all in a file.

Last modified: 14 Dec 2023
"""
import numpy as np
from numpy import ma
#from netCDF4 import num2date

import xarray as xr
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.dates as mdates

from netCDF4 import num2date

from pycurrents.num import interp1

import gsw
import sys

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
# functions
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
# time is on a grid and in matlab datenum format
time_units = "days since 0001-01-01 00:00:00.0"

s1_lon = -116.358
s1_lat = -72.468

def decode_time(var, var_units):
  """
  Use xr.decode_cf to convert time as a number with a reference
  to datetime64[ns]. Uses a units attribute like 'days since 0001-01-01'

  > returns an array of datetime64[ns]
  """
  attrs = {"units": var_units}
  t = xr.Dataset({'time':('time', var, attrs)})
  t_cf = xr.decode_cf(t)
  return t_cf.time.values

def plot_vel_check(xg, yg, var, vlims, title):
  vmin = vlims[0]
  vmax = vlims[1]

  plt.ion()
  fig, ax = plt.subplots(figsize=(15, 4))
  axins1 = inset_axes(ax, width="30%",
                      height="3.5%",
                      loc='upper right')
  cs=ax.pcolormesh(xg, yg, var, 
                cmap=cm.viridis_r,#cm.RdBu_r, 
                vmin=vmin, vmax=vmax)
  ax.annotate(title, xy=(0.01, .93), 
             xycoords='axes fraction', 
             bbox=dict(fc='lavender', alpha=.7, ec='lavender'))
  fig.colorbar(cs, cax=axins1, orientation='horizontal')
  ax.grid(True, c='dimgrey', lw=.5)
  ax.invert_yaxis()
  plt.tight_layout()


################################################################
# restrict data to 526-390 m [Dotto, 2020, Dotson]
# changed to 404, 532 [Oct 2021]
step = 8 
dep_min, dep_max = 404, 548 #340, 532#382, 526  
nbins = int(1 + (dep_max - dep_min)/step)
gdepth = np.linspace(dep_min, dep_max, nbins) 

################################################################


# Directories
# --------------------------------------------------------
workdir = '/Volumes/SamT5//PhD/data/'
mdir = workdir + 'moorings/SWE_Dotson/'
mdir14 = workdir + 'moorings/SWE_Dotson/S114_adcp/'
mcdir14 = workdir + 'moorings/SWE_Dotson/S114_microcat/'

localdir = '/Users/ocd1n16/PhD_local/'
scriptdir = localdir + 'scripts/'

auxscriptdir = localdir + 'scripts/aux_func/'
sys.path.append(auxscriptdir)
import aux_func_trend as ft
import aux_matlab_func as matf
# --------------------------------------------------------
# 2010 - 2014
# averaged every hour (hh:00)
fnames = ['S110', 'S111', 'S112']
files = [n + '.mat' for n in fnames]

u0, v0, t0, s0, tim0 = [], [], [], [], []

for k in range(len(fnames)):
    dm = matf.loadmat(mdir + files[k])
    sm = dm[fnames[k]]
    mat = sm['mat']
    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~    
    print("\n >> MOORING: ", sm['name'])
    print(sm.keys(), '\n')
    print("hourly-averaged data: ", sm['mat'].keys(), '\n')
    print("raw (ADCP) data: ", sm['ins'][0].keys(), '\n')
    print('lon: %.3f | lat: %.3f' % (sm['lon'], sm['lat']))
    print('depth: %.2f' % sm['depth'])
    #print(sm['note'])

    # matlab references the number of days to 0000-00-00, 
    # whereas python to 0001-01-01 00:00.0 
    # -> subtract a year (365 days) from matlab time to account for that
    date_in = num2date(sm['in']-365, units=time_units, calendar='gregorian')
    date_out = num2date(sm['out']-365, units=time_units, calendar='gregorian')
    print("date in: ", date_in)
    print("date out: ", date_out)
    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~    

    # velocity has the shape time x depth levels (m above bottom)
    u = mat['u']
    v = mat['v']
    t = mat['t']
    s = mat['s']
    time = mat['time']-365

    # depth = depth_of_mooring - height_above_bottom
    vdepth = sm['depth'] - mat['mabv']
    tdepth = sm['depth'] - mat['mabt']
    sdepth = sm['depth'] - mat['mabs']

    u0.append(interp1(vdepth, u, gdepth, axis=1))
    v0.append(interp1(vdepth, v, gdepth, axis=1))
    t0.append(interp1(tdepth, t, gdepth, axis=1))
    s0.append(interp1(sdepth, s, gdepth, axis=1))
    tim0.append(time)

u_part1 = np.asarray([item for sublist in u0 for item in sublist])
v_part1 = np.asarray([item for sublist in v0 for item in sublist])

t_part1 = np.asarray([item for sublist in t0 for item in sublist])
s_part1 = np.asarray([item for sublist in s0 for item in sublist])

tim_part1 = ma.asarray([item for sublist in tim0 for item in sublist])
mdate_part1 = decode_time(tim_part1, time_units)

# VELOCITY
# - - - - - - - - - - 
xg, yg = np.meshgrid(mdate_part1, gdepth)
v_limits = [-5, 5]
#title_v = 'v rotated by 42 deg clockwise'
title_v = 'u [cm/s]'
plot_vel_check(xg, yg, u_part1.T*100, v_limits, title_v)


# TEMPERATURE
# - - - - - - - - - - 
fig, ax = plt.subplots(figsize=(15, 4))
#ax.pcolormesh(xg, yg, t10.T, 
#              cmap=cm.plasma, 
#              vmin=-1.5, vmax=1)
axins1 = inset_axes(ax, width="30%",
                    height="5%",
                    loc='upper right')
levels = np.linspace(-2, 1.5,8) 
cs = ax.contourf(xg, yg, t_part1.T, 
            cmap=cm.plasma,
            levels=levels,
            extend='both') 
ax.contour(xg, yg, t_part1.T, levels=[0],
           colors='k', linewidths=1)
ax.invert_yaxis()
fig.colorbar(cs, cax=axins1, orientation='horizontal')

ax.annotate('temperature (in-situ?)', xy=(0.01, .93), 
           xycoords='axes fraction', 
           bbox=dict(fc='lavender', alpha=.7, ec='lavender'))
ax.grid(True, c='dimgrey',ls=':', lw=.5)
plt.tight_layout()

# SALINITY
# - - - - - - - - - - 
fig, ax = plt.subplots(figsize=(15, 4))
axins1 = inset_axes(ax, width="30%",
                    height="5%",
                    loc='upper right')
levels = np.linspace(34., 34.8, 9) 
cs = ax.contourf(xg, yg, s_part1.T, 
            cmap=cm.viridis,
            levels=levels,
            extend='both') 
ax.contour(xg, yg, t_part1.T, levels=[0],
           colors='k', linewidths=1)
ax.invert_yaxis()
fig.colorbar(cs, cax=axins1, orientation='horizontal')

ax.annotate('salinity (in-situ?)', xy=(0.01, .93), 
           xycoords='axes fraction', 
           bbox=dict(fc='lavender', alpha=.7, ec='lavender'))
ax.grid(True, c='dimgrey',ls=':', lw=.5)
plt.tight_layout()
   

# --------------------------------------------------------
# velocity: 2014-2016
# --------------------------------------------------------  
dm = matf.loadmat(mdir14 + 'open_S1_ADCP_20142016_quarter.mat')
print("fields in file for S1 2014-2016: \n", dm.keys())

# time is on a grid and in matlab datenum format
tim_part2 = dm['Time_ADCP_2016'][:, 0]-365   
# datetime64[ns] format
mdate_part2 = decode_time(tim_part2, time_units)
#date14 = num2date(var, units=time_units, calendar='gregorian')
#mdate14 = mdates.num2date(var-1)

ds_part2 = xr.Dataset({'u': (('time', 'depth'), dm['u_adcp_S1_2016']),
                   'v': (('time', 'depth'), dm['v_adcp_S1_2016'])},
                   coords={'time': mdate_part2,
                           'depth': dm['Depth_ADCP_2016'][0, :]})
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  
# 1. hourly resampling
hds_part2 = ds_part2.resample(time="1H").mean()

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  
# 2. regrid on the same depth grid as above 
u_part2 = interp1(hds_part2.depth.values,
                hds_part2.u.values, gdepth, axis=1)
v_part2 = interp1(hds_part2.depth.values,
                hds_part2.v.values, gdepth, axis=1)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  
# 3. combine time series
u_all = np.concatenate((u_part1, u_part2), axis=0)
v_all = np.concatenate((v_part1, v_part2), axis=0)
tim_all = np.hstack((mdate_part1, hds_part2.time.values))

vs = xr.Dataset({'u': (('time', 'depth'), u_all),
                'v': (('time', 'depth'), v_all)},
                coords={'time': tim_all,
                        'depth': gdepth})


vs.u.attrs['long_name'] = 'hourly_zonal_velocity_adcp'
vs.u.attrs['units'] = 'm/s'
vs.v.attrs['long_name'] = 'hourly_meridional_velocity_adcp'
vs.v.attrs['units'] = 'm/s'

vs.attrs['lon'] = s1_lon
vs.attrs['lat'] = s1_lat

# compute depth average
vs_avg = vs.mean('depth')
u_avg = vs_avg.u.values

# discard outliers outside |0.2| m/s (some measurement anomalies in 2010-12)
cond_outliers = np.logical_and(u_avg>0.2, u_avg<-0.2)
vs.u[cond_outliers] = np.nan
vs.v[cond_outliers] = np.nan

# save vel file
vs.to_netcdf(mdir + '../dragomir_phd/s1_raw_uv_hourly_adcp_404m_548m_2010_2016.nc')


#---------------------------------------------------------
# temp, sal 2014-2016
# --------------------------------------------------------  

dm = matf.loadmat(mcdir14 + 'S1_SBE37_20142016_quarter.mat')
print("\n fields in file for S1 2014-2016: \n", dm.keys())

# 1. variables from microcat sensors
ts_tim_part2 = dm['time'][:,0]-365 # time is on a grid and in matlab datenum format
ts_mdate_part2 = decode_time(ts_tim_part2, time_units)

# pressure (sensors drift in time!)
# temperature (in-situ?) and salinity (practical salinity?)
p0 = ma.masked_invalid(dm['p'])
t0 = ma.masked_invalid(dm['t'])
s0 = ma.masked_invalid(dm['s'])

# 2.1 find and remove linear trend in pressure sensors (account for drift)
# 2.2 average the detrended pressure time series to find the depths
r, c = p0.shape
avp = np.ones((c,))
for i in range(c):
    pc = p0[:, i]
    tc = ts_tim_part2[~pc.mask]
    pc = pc[~pc.mask]
    trend, _ = fc.trend_ci_np(tc, pc, 0.95)
    dt = np.hstack((0, tc[1:] - tc[:-1]))
    avp[i] = ma.mean(pc - trend.slope.values*dt)

# 3. re-grid T, S to match previous record
t_part2 = interp1(avp, t0, gdepth, axis=1)
s_part2 = interp1(avp, s0, gdepth, axis=1)

ts_part2 = xr.Dataset({'temp': (('time', 'depth'), t_part2),
                'sal': (('time', 'depth'), s_part2)},
                coords={'time': ts_mdate_part2,
                        'depth': gdepth})
# 3. resample hourly
hts_part2 = ts_part2.resample(time='1H').mean()

# 4. combine time series and resample hourly
t_all = np.concatenate((t_part1, hts_part2.temp.values), axis=0)
s_all = np.concatenate((s_part1, hts_part2.sal.values), axis=0)
date_all = ma.hstack((mdate_part1, hts_part2.time.values))

# convert to CT, SA
sa = gsw.SA_from_SP(s_all, gdepth, s1_lon, s1_lat)
ct = gsw.CT_from_t(sa, t_all, gdepth)

ts = xr.Dataset({'temp': (('time', 'depth'), t_all),
                'sal': (('time', 'depth'), s_all),
                'ct' : (('time', 'depth'), ct),
                'sa' : (('time', 'depth'), sa)},
                coords={'time': date_all,
                        'depth': gdepth})
ts.temp.attrs['long_name'] = 'temperature_microcats'
ts.ct.attrs['long_name'] = 'conservative_temperature'
ts.temp.attrs['units'] = 'degrees Celsius'
ts.sal.attrs['long_name'] = 'salinity_microcats'
ts.sa.attrs['long_name'] = 'absolute_salinity'
ts.sal.attrs['units'] = 'g/kg'
ts.sal.attrs['units'] = 'g/kg'
ts.depth.attrs['units'] = 'm'
ts.attrs['lon'] = s1_lon
ts.attrs['lat'] = s1_lat

# SAVE temp-sal file
ts.to_netcdf(mdir + '../dragomir_phd/s1_raw_hourly_microcat_404m_532m_2010_2016.nc')

