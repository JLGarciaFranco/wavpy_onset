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
#### Main Function haarcovtransform #####################
# 6 inputs:
# allprf= Retrodispersion profile matrix[m,n] where
# m = len(height array) and n = len(time array)
# z = height array (vector)
# i = index for current time value. (int)
# a = Dilation type: Automated or Standard (str)
# f = Skip between each height step in resolution (int)
def haarcovtransfm(allprf,a):
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
	return wf
# if not Automated, Standard implies using standard dilation a = 60m (found by Grabon and useful for UNAM profiles.
#	elif a=='Standard':
#		a=60
#		detail=False
# Call recursive function to find mlh (clearing out bottom values, i.e., no floor mlh value is permitted.
# Inspect firstmlh function below if needed.
	newmlh,wf=firstmlh(prf,a,bottom)
# If Automated, find top, bot and mlh using recursive algorithm.
	if detail:
		bot,newmlh,top=findtops(prf,wf,newmlh,a)
	else:
		bot=newmlh
		top=newmlh
### RETURN MLH ####
	return bot,newmlh,top


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
