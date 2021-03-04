import math, sys, os, shutil, string, cmath
import numpy as np
import bottleneck as bn
import pandas as pd

import astropy
from astropy.coordinates import SkyCoord
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
from astropy.io import ascii
import astropy
from astropy import stats 
from astropy.io import fits
from astropy import table
from astropy import convolution
from numpy import mean,median

from decimal import *
from progressbar import *
from os.path import expanduser
from itertools import tee, islice, chain, izip

from termcolor import colored
from pandas import DataFrame

import scipy
from scipy import stats
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import gaussian_filter1d
from scipy import signal
from scipy.constants import physical_constants
import scipy.optimize as opt
import warnings
from scipy.optimize import OptimizeWarning
import scipy.integrate as integrate
from scipy.special import wofz


import logging
import itertools
####Fnc_Stk_Mth####
def q_2_incl(b_a_ratio):
	q = 0.22
	angle = math.degrees(math.acos(np.sqrt(((b_a_ratio**2) - (q**2)) /(1-(q**2)))))
	return angle

def fwhm2sigma(fwhm):
	return fwhm/(2*np.sqrt(2*np.log(2)))

def sigma2fwhm(sigma):
	return sigma*(2*np.sqrt(2*np.log(2)))

def astpy_conv_gaus_1dkernel(kernel):
	#http://docs.astropy.org/en/stable/api/astropy.convolution.Gaussian1DKernel.html#astropy.convolution.Gaussian1DKernel
	return astropy.convolution.Gaussian1DKernel(stddev=kernel)#,mode='integrate')#center, linear_interp,oversample,integrate

def split_list(alist, wanted_parts=1):
	length = len(alist)
	return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
			for i in range(wanted_parts) ]
def Lmbd_2_Vel(lambda_reference,lambda_observed,args,**kwargs):
	velocity = c_speed * (1- lambda_reference/lambda_observed)
	return velocity

def lw_sgma2fwhm(linewidth):
    return np.sqrt(8*np.log(2))*linewidth

def func_Linear(X,slope,b):
    return (slope*X) + b

def func_1D_Gaussian_Emm(X,X_0,A,SIGMA):
	if X_0 == 999999.99999 and A == 999999.99999 and SIGMA == 999999.99999:
		X_0   = 0
		A     = 0
		SIGMA = 0
		k     = 1
	else:
		pass
		k = 1
	return k + (A*np.exp(-(X-X_0)**2/(2*SIGMA**2)))

def func_1D_Gaussian(X,X_0,A,SIGMA):
	if X_0 == 999999.99999 and A == 999999.99999 and SIGMA == 999999.99999:
		X_0   = 0
		A     = 0
		SIGMA = 0
		k     = 1
	else:
		pass
		k = 1
	return k + (A*np.exp(-(X-X_0)**2/(2*SIGMA**2)))

def func_1D_Gaussian_O(X,X_0,A,SIGMA,OFFSET):
	if X_0 == 999999.99999 and A == 999999.99999 and SIGMA == 999999.99999 and OFFSET == 999999.99999:
		X_0    = 0
		A      = 0
		SIGMA  = 0
		OFFSET = 0
		k      = 1
	else:
		k = 1
	return OFFSET + k + (A*np.exp(-((X-X_0)**2)/(2*SIGMA**2)))

def func_2_1D_Gaussian(X,X_0_1,A_1,SIGMA_1,X_0_2,A_2,SIGMA_2):
    return 1 + (A_1*np.exp(-(X-X_0_1)**2/(2*SIGMA_1**2))) + (A_2*np.exp(-(X-X_0_2)**2/(2*SIGMA_2**2)))

def func_2_1D_Gaussian_O(X,X_0_1,A_1,SIGMA_1,X_0_2,A_2,SIGMA_2, OFFSET):
    return OFFSET + 1 + (A_1*np.exp(-(X-X_0_1)**2/(2*SIGMA_1**2))) + (A_2*np.exp(-(X-X_0_2)**2/(2*SIGMA_2**2)))

def func_Lorentzian(X,X_0,GAMMA):
    return 1 + (GAMMA / np.pi / ((X-X_0)**2 + GAMMA **2))

def func_Voigt(X,X_0,ALPHA,GAMMA):
    SIGMA = (ALPHA / np.sqrt(2 * np.log(2)))
    return 1 + np.real(wofz(((X-X_0)+1j*GAMMA) / SIGMA / np.sqrt(2))) / SIGMA / np.sqrt(2*np.pi)

def Cosmo_Scale_Par(z_cosmo_parm,Uni_Cosmo_Pars,*args, **kwargs):
	#Info obtained via:
	#http://www.astro.ucla.edu/~wright/CosmoCalc.html
	#For Ho = 70, OmegaM = 0.270, Omegavac = 0.730, z = 0.020
	#
	#It is now 13.861 Gyr since the Big Bang.
	#The age at redshift z was 13.585 Gyr.
	#The light travel time was 0.275 Gyr.
	#The comoving radial distance, which goes into Hubble's law, is 85.3 Mpc or 0.278 Gly.
	#The comoving volume within redshift z is 0.003 Gpc3.
	#The angular size distance DA is 83.634 Mpc or 0.272778 Gly.
	#This gives a scale of 0.405 kpc/".
	#The luminosity distance DL is 87.0 Mpc or 0.284 Gly.
	#1 Gly = 1,000,000,000 light years or 9.461*1026 cm.
	#1 Gyr = 1,000,000,000 years.
	#1 Mpc = 1,000,000 parsecs = 3.08568*1024 cm, or 3,261,566 light years.

	#astropy.cosmology.FLRW.arcsec_per_kpc_comoving
	#print cosmo.H(0)
	#print cosmo.Om(0)
	#print cosmo.Tcmb(0)
	#print cosmo.Neff
	#print cosmo.m_nu
	scale    = (1/Uni_Cosmo_Pars.arcsec_per_kpc_proper(z_cosmo_parm)) #kpc/arcsec
	scaleL   = Uni_Cosmo_Pars.luminosity_distance(z_cosmo_parm) #dlGrp Mpc
	c_speed  = physical_constants['speed of light in vacuum'][0]/1000
	DELTAZ   = rad_vel_sep/c_speed
	sepkpc   = rad_sep[1]*scale

	cosmo_check = kwargs.get('cosmo_check', False)

	if cosmo_check == True:

		print ''
		print 'Cosmology:'
		print
		print 'Cosmological parameters:'
		print Uni_Cosmo_Pars
		print 
		print 'At redshift                                      : ', z_cosmo_parm 
		print 'The angular separations [arc sec]                : ', rad_sep[1]
		print 'Angular distance scale                           : ', scale #,'[kpc/arcsec]'
		print 'The angular separations [kpc]                    : ', sepkpc.value
		print 'Luminosity distance(D_L)                         : ', scaleL#',[Mpc]'
		print ''
		print 'A radial velocity difference (DeltaV_r)          : ',rad_vel_sep
		print 'Corresponds to a redshift difference (Delta_z) of: ',DELTAZ,DELTAZ
		print 'c                                                : ',c_speed
		print
	elif cosmo_check == False:
		pass

	return scale.value,scaleL.value,sepkpc
####Fnc_Stk_Mth####
