import numpy as np 
import pandas as pd
import glob,sys,os
import xarray as xr
import datetime
import matplotlib.pyplot as plt
from wavelets import *

# dilation coefficients for daily data, used for onset
daya=np.arange(70,320,2)
# dilation coefficients for calculation of MSD
dayo=np.arange(50,100,2)
# flag to determine if midsummer drought or double peak structures are to be determined
analyze_MSD=True


# this where one defines the path and or file where precip data is stored
mainpath='/xpt/berimbau.local/data3/jorgegf/obs_precip/'
data_filename=mainpath+'CMORPH/'
flist=glob.glob(data_filename+'*.nc')
#load data
precipitation_data=xr.open_mfdataset(flist)['rain']
print(precipitation_data)
# slice data for a grid point
series=precipitation_data.sel(lat=16,lon=-91,method='nearest')
print(series)
#select a year we want to analyze 
year=2014
print(series.time)
# following conditions select the year analyzed and the two years right before and right after it
year_series=series.where((series.time.dt.year>=year-1)&(series.time.dt.year<=year+1),drop=True)
# the method is best when longer time series are analyzed, so what we do here is append a few months of the year before (Nov-Dec) and a few months of the following year (Jan-Apr)
year_series=year_series.where((year_series.time.dt.year==year)|((series.time.dt.year==year+1)&(series.time.dt.month<5))|
		((series.time.dt.year==year-1)&(series.time.dt.month>10)),drop=True)
newt=pd.to_datetime(year_series.time)
haarsum=haarcovtransfm(year_series.values,daya,summed=True)
onset_i=np.where(haarsum==np.max(haarsum))[0]
onset_f=np.where(haarsum==np.min(haarsum))[0]
print(year,'Onset date ',newt[onset_i],'Retreat date',newt[onset_f])
if analyze_MSD:
	# the threshold days parameter is meant to determine where the precipitation time-series is cutoff to analyze the MSD timing only.
	# 2-10 values are typically good enough and the results should not be sensitive to this choice when in this parameter space.
	threshold_days=5
	indices=np.arange(onset_i-threshold_days,onset_f+threshold_days)
	msd_precipitation=year_series.isel(time=indices)
	msd_dates=pd.to_datetime(msd_precipitation.time)
	haarsum2=haarcovtransfm(msd_precipitation.values,dayo,summed=True)
	msde=np.where(haarsum2==np.max(haarsum2))[0]
	msdo=np.where(haarsum2==np.min(haarsum2))[0]
	msdo_date=msd_dates[msdo]
	msde_date=msd_dates[msde]
	print(year,'MSD Onset date ',msdo_date,'End of MSD date',msde_date)
# write results to txt file
f=open('MSD_stats.txt','a')
init = datetime.datetime(year,4,7)
end = datetime.datetime(year,11,30)
fig=plt.figure(figsize=(10.1,6.1),dpi=170)
plt.subplot(211)
#plt.subplot(221+idd)	
#plt.plot(msd_dates,msd_precipitation,c='k',linewidth=2,label='All year precip')
plt.plot(newt,year_series,c='navy',linewidth=2.5,linestyle='--',label='Precipitation')
plt.axvline(x=newt[onset_i],c='dodgerblue',label='Onset')
plt.axvline(x=newt[onset_f],c='orangered',label='Demise')
plt.xlim([init,end])
plt.grid()
plt.legend()
plt.ylabel(r'pr mm day$^{-1}$',fontsize=12)
ax=plt.subplot(212)
plt.plot(newt,haarsum,c='k',label='Long',linewidth=2)
plt.xlim([init,end])
plt.axvline(x=newt[onset_i],c='dodgerblue',label='Onset')
plt.axvline(x=newt[onset_f],c='orangered',label='Demise')
plt.ylabel(r'$\sum W_f$',fontsize=12)
plt.grid()
plt.legend()
plt.tight_layout(h_pad=0.13,w_pad=0.3,pad=0.3)
plt.savefig('/home/jorgegf/plots/msd_example.png')

