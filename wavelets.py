#!/usr/bin/python/env
# -*- coding: utf8 -*-

###########  Wavelet toolbox ############################
### Following Brooks (2003) and Grabon (2010) ###########
### by Jorge L. García Franco on May-Aug 2016 ###########
### contact: jgcaspark@ciencias.unam.mx or hotmail.com ##
#########################################################
# Import numpy for math processing
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
import xarray as xr
from typing import List, Union
import logging
#### Main Function haarcovtransform #####################
# 6 inputs:
# allprf= Retrodispersion profile matrix[m,n] where
# m = len(height array) and n = len(time array)
# z = height array (vector)
# i = index for current time value. (int)
# a = Dilation type: Automated or Standard (str)
# f = Skip between each height step in resolution (int)
def haarcovtransfm(allprf,a,meaned=True,summed=False):
#### Declare global variables, with a0 being the minimum dilation posible.
####
	global fi
	global z
	global a0
	global iz
	global b
	global top
## Bottom is 80 m due to noise generated artifacts by surface, newmlh is initially this value to jumpstart.
## wf = wavelet transform coefficients array, initially  empty.
	b=range(len(allprf))
# Selecting current profile for current time.
#	prf=allprf[:,i]
# a = Automated implies following recursive algorithm for finding transition zone by algorithm of Brooks.
#	if a=='Auto':
	detail=True
	wf=wfab(allprf,a)
	if summed:
		return np.sum(wf,axis=1)
	elif meaned:
		return np.mean(wf,axis=1)
	else:
		return wf
# if not Automated, Standard implies using standard dilation a = 60m (found by Grabon and useful for UNAM profiles.
#	elif a=='Standard':
#		a=60
#		detail=False

#######################################################################################################################
### haarval = function to compute wavelet coefficient Wf(a,b) for every a,b.
def haarval(prf,a,b0):
	global z
# wnlen is the window size to compute the positive and negative pulse.
	wnlen=a/2
	fun=0
# Loop through z. Assigning weight according to Haar.

	return fun/a

#######################################################################################################################
### Function findtops, finds recursively the top and bottom of the transition zone, and mlh by
### finding the suitable dilation.
def findtops(prf,wf,newmlh,a):
### First round of coefficients to find top (c1) and bottom (c2) heights.
	c1=0.6
	c2=0.4
	global a0
	bt=b
	while a>a0:
		maxi=np.max(wf)
		imaxi=np.argmax(wf)
		#Top index retrieval
		topindex=0
		wf6=wf[imaxi]
		while wf6 > c1*maxi and imaxi+topindex != len(wf)-1:
			wf6=wf[imaxi+topindex]
			topindex+=1
		#Bottom index retrieval
		botindex=1
		wf4=wf[imaxi-botindex]
		while wf4 > c2*maxi and imaxi-botindex!=0:
			botindex+=1
			wf4=wf[imaxi-botindex]
		a=a-20
		if bt[imaxi+topindex-1]-bt[imaxi]<=a0 or bt[imaxi]-bt[imaxi-botindex]:
			break
		bt=bt[imaxi-botindex:imaxi+topindex]
		c1=c1-0.02
		c2=c2+0.02
		wf=[]
		### Find wavelet transform coefficients given current dilation.
		for n,b0 in enumerate(bt):
			covtransform=haarval(prf,a,b0)
			wf.append(covtransform)
	return bt[imaxi-botindex],bt[imaxi],bt[imaxi+topindex-1]
####################################################################################
### firstmlh: function to obtain first approximation to mlh given the first dilation observed.
### It is written to avoid ceiling or floor mlh values being floor = 200 and top =3000.
def wfab(prf,a):
	index=0
	wf=np.zeros((len(prf),len(a)))
	for index,ai in enumerate(a):
	#Loop until newmlh is not current bottomo or top
#	while newmlh<=bottom+10 or newmlh >= top-10:
		b=np.arange(ai/2.,len(prf)-ai/2.)
		#print(b)
		for ib,bi in enumerate(b): 
			fun=0
	#		print(bi-ai/2.,bi,bi+ai/2.+1,ai)
	#		print(np.arange(bi-ai/2.,bi+ai/2+1))
	#		print(prf[int(bi-ai/2.):int(bi+ai/2.)+1],int(ai))
			apost=prf[int(bi+1):int(bi+ai/2.)+1]
			aprior=prf[int(bi-ai/2.):int(bi)]
			wf[int(bi),index]=(np.nansum(apost)-np.nansum(aprior))/float(ai)
			continue			
#			shortprf=prf[bi-	
			for i,z0 in enumerate(z):
				if z0 < b0-wnlen or z0 == b0:
					continue
				elif z0 > b0+wnlen:
					break
				elif z0 >= b0-wnlen and z0 < b0:
					fun=fun+prf[i]
				elif z0 > b0 and z0 <= b0+wnlen:
					fun=fun-prf[i]
			
#		break	

	return wf
def calc_msd(
    dataset_name: str,
    year: int,
    ds: xr.DataArray,
    output_file: str = "MSD_obs_table.txt"
) -> List[Union[int, datetime.date, float]]:
    """
    Calculate Midsummer Drought (MSD) timing metrics for a given year.

    This function computes the onset and end dates of the rainy season and the midsummer drought
    using a Haar wavelet covariance transform on daily precipitation data. Results are appended
    to a CSV-style text file.

    Parameters
    ----------
    dataset_name : str
        Identifier for the dataset, used in logging output.
    year : int
        Year for which the MSD metrics are calculated.
    ds : xr.DataArray
        Time-indexed precipitation data. Must have a 'time' coordinate.
    output_file : str, optional
        Path to the output CSV file where results are appended. Default is 'MSD_obs_table.txt'.

    Returns
    -------
    results : list
        A list containing:
        - year (int)
        - rainy season onset_date (datetime.date)
        - rainy season end_date (datetime.date)
        - MSD onset within season (datetime.date)
        - MSD end within season (datetime.date)
        - coef3 (float): minimum wavelet coefficient prior to MSD
        - coef1 (float): maximum wavelet coefficient during MSD
        - coef2 (float): amplitude (coef1 - coef3)
        - jjas_mean (float): mean precipitation between onset and end

    Notes
    -----
    - Uses `haarcovtransfm` to compute the Haar covariance transform.
    - Onset and end are determined by the max/min of the seasonal transform.
    - Within-season MSD metrics use a shorter window around the drought period.
    """
    # Configure logging
    output_file=dataset_name+output_file
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Define Haar transform windows based on dataset length
    if ds.size < 400:
        window_season = np.arange(70, 120, 2)
    else:
        window_season = np.arange(50, 140, 2)
    window_msd = np.arange(34, 94, 2)
    threshold = 60

    # Define analysis period: Oct 15 previous year to May 10 next year
    date_start = datetime.date(year - 1, 10, 15)
    date_end = datetime.date(year + 1, 5, 10)
    print(date_start,date_end)
    print(ds)
    # Subset the precipitation series
    precip = ds.sel(
        time=slice(pd.to_datetime(date_start), pd.to_datetime(date_end))
    )
    
    logger.info(
        "Calculating MSD for %s, year %d: %d time points",
        dataset_name,
        year,
        precip.size,
    )
    print(precip)
    # Compute Haar covariance transform over the full season
    times = pd.to_datetime(precip.time.values)
    season_transform = haarcovtransfm(precip, window_season)

    # Identify rainy season onset and end by max/min of transform
    onset_idx = int(np.argmax(season_transform))
    end_idx = int(np.argmin(season_transform))
    onset_date = times[onset_idx].date()
    end_date = times[end_idx].date()
    print(onset_date,end_date)
    logger.info(
        "Year %d rainy season: onset %s, end %s",
        year,
        onset_date,
        end_date,
    )

    # Extract data between onset and end for MSD analysis
    rainy_period = precip.sel(
        time=slice(pd.to_datetime(onset_date), pd.to_datetime(end_date))
    )
    rainy_times = pd.to_datetime(rainy_period.time.values)

    # Compute Haar transform for midsummer drought period
    msd_transform = haarcovtransfm(rainy_period, window_msd)

    # Identify MSD onset and end indices within rainy period
    msd_onset_idx = int(np.argmin(msd_transform[:-threshold]))
    msd_end_idx = int(np.argmax(msd_transform[threshold:])) + threshold
    msd_onset_date = rainy_times[msd_onset_idx].date()
    msd_end_date = rainy_times[msd_end_idx].date()

    # Compute coefficients and mean precipitation
    coef3 = float(np.min(msd_transform[:-threshold]))
    coef1 = float(np.max(msd_transform[threshold:]))
    coef2 = coef1 - coef3
    jjas_mean = float(rainy_period.mean().item())

    # Prepare results list
    results = [
        year,
        onset_date,
        end_date,
        msd_onset_date,
        msd_end_date,
        coef3,
        coef1,
        coef2,
        jjas_mean,
    ]

    # Append to output file
    with open(output_file, "a") as f:
        f.write(
            ",".join(str(x) for x in results) + "\n"
        )

    return results
