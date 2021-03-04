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
from progressbar import *               # just a simple progress bar
from os.path import expanduser
from itertools import tee, islice, chain, izip
from pysynphot import observation
from pysynphot import spectrum
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

import pyraf
from pyraf import iraf

import logging
import itertools

from Fnc_Stk_Fts import *
from Fnc_Stk_Mth import *
from Fnc_Stk_Dir import *

import Lines_Dictionary
LINES = Lines_Dictionary.LINES_STK #FOR MASKING##

####Fnc_Stk_Spc####
def rebin_spec(wave, specin, wavnew,*args, **kwargs):
	#http://pysynphot.readthedocs.io/en/latest/ref_api.html#module-pysynphot.observation
	#http://pysynphot.readthedocs.io/en/latest/ref_api.html#pysynphot.observation.Observation.validate_overlap
	#http://pysynphot.readthedocs.io/en/latest/ref_api.html#pysynphot.observation.validate_overlap
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew, force=None)#extrap taper None
    return obs.binflux

def new_wave_range(new_lambda_zero,new_lambda_end,new_lambda_step,*args, **kwargs):
	new_zs_lambda = ((new_lambda_end-new_lambda_zero)/new_lambda_step)+1
	new_sp_lambda = np.arange(new_lambda_zero,new_lambda_zero+(new_zs_lambda*new_lambda_step),new_lambda_step)
	new_zs_lambda = len(new_sp_lambda)
	return new_sp_lambda,new_lambda_step,new_zs_lambda,new_lambda_end,new_lambda_zero

def SNR(flux):
	"""
	flux = array(flux)
	Values that are exactly zero (padded) are skipped
	flux = array(flux[where(flux != 0.0)])
	https://stdatu.stsci.edu/vodocs/der_snr.pdf
	http://www.stecf.org/software/ASTROsoft/DER_SNR/
	"""
	flux = flux[np.where(flux != 0.0)]
	n    = len(flux)      
	# For spectra shorter than this, no value can be returned
	if (n>4):
	  signal = bn.nanmedian(flux)
	  noise  = 0.6052697 * bn.nanmedian(abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
	  return float(signal / noise),signal,noise 
	else:
	  return 666.6,666.6,666.6

def Spectra_x_y(specfile,*args, **kwargs):
	hdulist_sp  = fits.open(specfile)
	sp_intens   = hdulist_sp[0].data
	sp_lambda0  = Header_Get(specfile,'CRVAL1')
	sp_step     = Header_Get(specfile,'CDELT1')
	sp_header   = hdulist_sp[0].header
	header_lbl  = kwargs.get('header_lbl','CDELT1')
	sp_hdr      = Header_Get(specfile,header_lbl)
	
	if sp_intens is None:
		print 'None Type'
		np.array(sp_intens, dtype=np.float)
		sp_intens    = np.empty(1117)
		sp_intens[:] = np.nan
	else:
		pass

	if len(sp_intens)==1:
		sp_lambda   = np.arange(sp_lambda0,sp_lambda0+(len(sp_intens[0])*sp_step),sp_step)
		sp_intens   = sp_intens[0][:]
	elif len(sp_intens) != 1:
		sp_lambda   = np.arange(sp_lambda0,sp_lambda0+(len(sp_intens)*sp_step),sp_step)

	if len(sp_lambda) == len(sp_intens) + 1:
		sp_lambda   = sp_lambda[0:len(sp_lambda) - 1]
	elif len(sp_lambda) == len(sp_intens) - 1:
		sp_lambda   = sp_lambda[0:len(sp_lambda) + 1]
	return sp_lambda,sp_intens,sp_lambda0,sp_step,sp_hdr,specfile

def Spectra_x_y_Updt(xyupdt_ifn,msk_typ,x_min,x_max,x_type,*args, **kwargs):
	res_sp_s   = Spectra_x_y(xyupdt_ifn)

	if x_type == 'lambda':
		n_idx_in   = int(abs(x_min - res_sp_s[2]) / res_sp_s[3])
		n_idx_fn   = int(abs(x_max - res_sp_s[2]) / res_sp_s[3]) + 1
	elif x_type == 'pixels':
		n_idx_in   = int(abs(x_min))
		n_idx_fn   = int(abs(x_max))
	elif x_type == 'all':
		n_idx_in   = 0
		n_idx_fn   = -1

	hdulist_sp  = fits.open(xyupdt_ifn, mode='update')
	hdulist_sp[0].data  = hdulist_sp[0].data.flatten()

	if msk_typ == 'constant':
		cnt_val = kwargs.get('cnt_val', None)
		hdulist_sp[0].data[n_idx_in:n_idx_fn] = cnt_val
	elif msk_typ == 'NaN':
		hdulist_sp[0].data[n_idx_in:n_idx_fn] = np.nan
	elif msk_typ == 'continuum':
		spfn_i_2 = kwargs.get('spfn_i_2', None)
		res_sp_c = Spectra_x_y(spfn_i_2)
		hdulist_sp[0].data[n_idx_in:n_idx_fn] = res_sp_c[1][n_idx_in:n_idx_fn]
	elif msk_typ == 'all':
		new_val_array = kwargs.get('new_val_array',None)
		if len(new_val_array) != len(res_sp_s[1]):
			print 'Error'
		elif len(new_val_array) == len(res_sp_s[1]):
			pass
		hdulist_sp[0].data[n_idx_in:n_idx_fn] = new_val_array

	hdulist_sp.flush()
	hdulist_sp.close()
	res_sp_m   = Spectra_x_y(xyupdt_ifn)
	return res_sp_m

def Spectra_x_y_Get(xyupdt_ifn,x_min,x_max,x_type,*args, **kwargs):
	res_sp_s   = Spectra_x_y(xyupdt_ifn)

	if x_type == 'lambda':
		n_idx_in   = int(abs(x_min - res_sp_s[2]) / res_sp_s[3])
		n_idx_fn   = int(abs(x_max - res_sp_s[2]) / res_sp_s[3]) + 1
	elif x_type == 'pixels':
		n_idx_in   = int(abs(x_min))
		n_idx_fn   = int(abs(x_max))
	elif x_type == 'all':
		n_idx_in   = 0
		n_idx_fn   = -1

	hdulist_sp  = fits.open(xyupdt_ifn, mode='update')
	hdulist_sp[0].data  = hdulist_sp[0].data.flatten()
	intens_values = hdulist_sp[0].data[n_idx_in:n_idx_fn]
	hdulist_sp.flush()
	hdulist_sp.close()
	return intens_values

def Spectra_Masking(mask_ifn,msk_typ,*args, **kwargs):
	rshft_corr        = kwargs.get('rshft_corr', 0)
	rshft_corr_direct = kwargs.get('rshft_corr_direct',False)

	msk_abs_lne       = kwargs.get('msk_abs_lne'  ,False)
	msk_blu_rgn       = kwargs.get('msk_blu_rgn'  ,False)
	blu_lmb_min       = kwargs.get('blu_lmb_min'  ,500)
	blu_lmb_max       = kwargs.get('blu_lmb_max'  ,1200)
	rwt_file          = kwargs.get('rwt_file'     ,False)
	msk_cte_val       = kwargs.get('msk_cte_val'  ,0)

	if rwt_file == False:
		mask_ofn   = str(mask_ifn.split('.fits',1)[0]) + '-m.fits'
		org_spec   = Spectra_x_y(mask_ifn)
		os.system('cp ' + mask_ifn   + ' ' + mask_ofn)
		Header_History_Step(mask_ifn,mask_ofn)
	elif rwt_file==True:
		mask_ofn = mask_ifn

	if msk_abs_lne == True:
		for lines in range(len(LINES[0])):
			if rshft_corr_direct == True:
				lambda_cen = LINES[0][lines] * (rshft_corr)
				lmb_min = (lambda_cen - (LINES[1][lines]*(rshft_corr)))
				lmb_max = (lambda_cen + (LINES[1][lines]*(rshft_corr)))
			elif rshft_corr_direct == False:
		 		lambda_cen = LINES[0][lines] * (1+rshft_corr)
				lmb_min = (lambda_cen - (LINES[1][lines]*(1+rshft_corr)))
				lmb_max = (lambda_cen + (LINES[1][lines]*(1+rshft_corr)))

			if (org_spec[0][0]<lambda_cen<org_spec[0][-1]) and (org_spec[0][0]<lmb_min<org_spec[0][-1]) and (org_spec[0][0]<lmb_max<org_spec[0][-1]):
				Spectra_x_y_Updt(mask_ofn,msk_typ,lmb_min,lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
			elif (org_spec[0][0]<lambda_cen<org_spec[0][-1]) and (org_spec[0][0]>lmb_min<org_spec[0][-1]) and (org_spec[0][0]<lmb_max<org_spec[0][-1]):
				lmb_min = org_spec[0][0]
				Spectra_x_y_Updt(mask_ofn,msk_typ,lmb_min,lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
			elif (org_spec[0][0]<lambda_cen<org_spec[0][-1]) and (org_spec[0][0]<lmb_min<org_spec[0][-1]) and (org_spec[0][0]<lmb_max>org_spec[0][-1]):
				lmb_max = org_spec[0][-1]
				Spectra_x_y_Updt(mask_ofn,msk_typ,lmb_min,lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
			else:
				pass
	elif msk_abs_lne == False:
		pass
	if msk_blu_rgn == True:
		if rshft_corr_direct == True:
			blu_lmb_min = blu_lmb_min * (rshft_corr)
			blu_lmb_max = blu_lmb_max * (rshft_corr)
		elif rshft_corr_direct == False:
			blu_lmb_min = blu_lmb_min * (1+rshft_corr)
			blu_lmb_max = blu_lmb_max * (1+rshft_corr)
		if (org_spec[0][0]<blu_lmb_min<org_spec[0][-1]) and (org_spec[0][0]<blu_lmb_max<org_spec[0][-1]):
			Spectra_x_y_Updt(mask_ofn,msk_typ,blu_lmb_min,blu_lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
		elif (org_spec[0][0]>blu_lmb_min<org_spec[0][-1]) and (org_spec[0][0]<blu_lmb_max<org_spec[0][-1]):
			blu_lmb_min = org_spec[0][0]
			Spectra_x_y_Updt(mask_ofn,msk_typ,blu_lmb_min,blu_lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
		elif (org_spec[0][0]<blu_lmb_min<org_spec[0][-1]) and (org_spec[0][0]<blu_lmb_max>org_spec[0][-1]):
			blu_lmb_max = org_spec[0][-1]
			Spectra_x_y_Updt(mask_ofn,msk_typ,blu_lmb_min,blu_lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
		else:
			pass		
	elif msk_blu_rgn == False:
		pass
	return mask_ofn

def Spectra_Cont_IRAF(Cont_ip_sfn_IRAF,Cont_log_IRAF,*args, **kwargs):
	Cont_type_IRAF     = kwargs.get('Cont_type_IRAF'    ,'ratio'  )  # Continuum fitting type fit,ratio,difference
	Cont_lines_IRAF    = kwargs.get('Cont_lines_IRAF'   ,'*'      )  # Image lines to be fit
	Cont_funct_IRAF    = kwargs.get('Cont_funct_IRAF'   ,'spline3')  # Fitting function: legendre, chebyshev, spline1, spline3
	Cont_order_IRAF    = kwargs.get('Cont_order_IRAF'   ,49       )  # Order Polynomial / num pieces spline
	Cont_override_IRAF = kwargs.get('Cont_override_IRAF','yes'    )  # Override previous norm spec
	Cont_replace_IRAF  = kwargs.get('Cont_replace_IRAF' ,'no'     )  # Replace rejected points by fit?
	Cont_low_rej_IRAF  = kwargs.get('Cont_low_rej_IRAF' ,3        )  # Low rejection in sigma of fit
	Cont_high_rej_IRAF = kwargs.get('Cont_high_rej_IRAF',3        )  # High rejection in sigma of fit

	Cont_op_sfn_IRAF_1      = str(Cont_ip_sfn_IRAF.split('.fits',1)[0]) + '-c.fits'
	if os.path.exists(Cont_op_sfn_IRAF_1)==True:
		os.system('rm ' + str(Cont_op_sfn_IRAF_1))
	else:
		pass
	iraf.continuum.input    = Cont_ip_sfn_IRAF
	iraf.continuum.output   = Cont_op_sfn_IRAF_1
	iraf.continuum.type     = Cont_type_IRAF
	iraf.continuum.lines    = Cont_lines_IRAF
	iraf.continuum.band     = 1
	iraf.continuum.logfile  = Cont_log_IRAF
	iraf.continuum.function = Cont_funct_IRAF
	iraf.continuum.order    = Cont_order_IRAF
	iraf.continuum.replace  = Cont_replace_IRAF
	iraf.continuum.low_rej  = Cont_low_rej_IRAF
	iraf.continuum.high_rej = Cont_high_rej_IRAF
	iraf.continuum.override = Cont_override_IRAF
	iraf.continuum.interac  = 'no'
	iraf.continuum.mode     = 'h'
	iraf.continuum()

	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_TYP',Cont_type_IRAF    ,header_comment='Continuum IRAF type')
	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_FIT',Cont_funct_IRAF   ,header_comment='Continuum IRAF function')
	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_ORD',Cont_order_IRAF   ,header_comment='Continuum IRAF order')
	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_RPC',Cont_replace_IRAF ,header_comment='Continuum IRAF replacement')
	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_LRJ',Cont_low_rej_IRAF ,header_comment='Continuum IRAF low sigma rejection')
	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_HRJ',Cont_high_rej_IRAF,header_comment='Continuum IRAF high sigma rejection')

	try:
		Header_Copy(Cont_op_sfn_IRAF_1,Cont_ip_sfn_IRAF,'h_s_c')
	except KeyError:
		Header_Get_Add(Cont_op_sfn_IRAF_1,'h_s_c',0,header_comment='Step 0')
		pass
	try:
		Header_Copy(Cont_op_sfn_IRAF_1,Cont_ip_sfn_IRAF,'h_s_0')
	except KeyError:
		Header_Get_Add(Cont_op_sfn_IRAF_1,'h_s_0',Cont_ip_sfn_IRAF,header_comment='Step 0')
		pass
	Header_History_Step(Cont_ip_sfn_IRAF,Cont_op_sfn_IRAF_1)

	Cont_op_sfn_IRAF_2      = str(Cont_ip_sfn_IRAF.split('.fits',1)[0]) + '-c-f.fits'
	if os.path.exists(Cont_op_sfn_IRAF_2)==True:
		os.system('rm ' + str(Cont_op_sfn_IRAF_2))
	else:
		pass

	iraf.continuum.output   = Cont_op_sfn_IRAF_2
	iraf.continuum.type     = 'fit'
	iraf.continuum()

	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_TYP',Cont_type_IRAF    ,header_comment='Continuum IRAF type')
	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_FIT',Cont_funct_IRAF   ,header_comment='Continuum IRAF function')
	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_ORD',Cont_order_IRAF   ,header_comment='Continuum IRAF order')
	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_RPC',Cont_replace_IRAF ,header_comment='Continuum IRAF replacement')
	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_LRJ',Cont_low_rej_IRAF ,header_comment='Continuum IRAF low sigma rejection')
	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_HRJ',Cont_high_rej_IRAF,header_comment='Continuum IRAF high sigma rejection')

	try:
		Header_Copy(Cont_op_sfn_IRAF_2,Cont_ip_sfn_IRAF,'h_s_c')
	except KeyError:
		Header_Get_Add(Cont_op_sfn_IRAF_2,'h_s_c',0,header_comment='Step 0')
		pass
	try:
		Header_Copy(Cont_op_sfn_IRAF_2,Cont_ip_sfn_IRAF,'h_s_0')
	except KeyError:
		Header_Get_Add(Cont_op_sfn_IRAF_2,'h_s_0',Cont_ip_sfn_IRAF,header_comment='Step 0')
		pass

	Header_History_Step(Cont_ip_sfn_IRAF,Cont_op_sfn_IRAF_2)

	Header_Get_Add(Cont_ip_sfn_IRAF,'CNT_FN0',str((Cont_ip_sfn_IRAF.rsplit('/',1)[1]).split('.',1)[0])   ,header_comment='Continuum IRAF Spec file pre-cont')
	Header_Get_Add(Cont_ip_sfn_IRAF,'CNT_FNM',str((Cont_op_sfn_IRAF_1.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Continuum IRAF Spec file pst-cont')
	Header_Get_Add(Cont_ip_sfn_IRAF,'CNT_FNF',str((Cont_op_sfn_IRAF_2.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Continuum IRAF Spec file cnt-fit')

	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_FN0',str((Cont_ip_sfn_IRAF.rsplit('/',1)[1]).split('.',1)[0])   ,header_comment='Continuum IRAF Spec file pre-cont')
	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_FNM',str((Cont_op_sfn_IRAF_1.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Continuum IRAF Spec file pst-cont')
	Header_Get_Add(Cont_op_sfn_IRAF_1,'CNT_FNF',str((Cont_op_sfn_IRAF_2.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Continuum IRAF Spec file cnt-fit')

	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_FN0',str((Cont_ip_sfn_IRAF.rsplit('/',1)[1]).split('.',1)[0])   ,header_comment='Continuum IRAF Spec file pre-cont')
	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_FNM',str((Cont_op_sfn_IRAF_1.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Continuum IRAF Spec file pst-cont')
	Header_Get_Add(Cont_op_sfn_IRAF_2,'CNT_FNF',str((Cont_op_sfn_IRAF_2.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Continuum IRAF Spec file cnt-fit')

	return Cont_op_sfn_IRAF_1,Cont_op_sfn_IRAF_2

def Spectra_Smooth(spc_ipf_smt,krnl_typ_smt,krnl_sz_smt,*args,**kwargs):
	spc_opf_smt = str(spc_ipf_smt.split('.fits',1)[0]) + '-smt.fits'
	spc2bsmth   = Spectra_x_y(spc_ipf_smt)

	if krnl_typ_smt == 'gaussian':
		grid_krnl = astpy_conv_gaus_1dkernel(krnl_sz_smt)
	elif krnl_typ_smt == 'boxcar':
		grid_krnl = astropy.convolution.Box1DKernel(krnl_sz_smt)
	elif pkrnl_typ_smt == 'mexican':
		grid_krnl = astropy.convolution.MexicanHat1DKernel(krnl_sz_smt)

	#http://docs.astropy.org/en/stable/api/astropy.convolution.convolve.html
	spc_pst_smt = astropy.convolution.convolve(spc2bsmth[1], grid_krnl)#,boundary='extend') #fill wrap extend

	Wrt_FITS_File(spc_pst_smt,spc_opf_smt)

	Header_Get_Add(spc_opf_smt,'CRVAL1',spc2bsmth[2],header_comment='Smoothing-CRVAL1')
	Header_Get_Add(spc_opf_smt,'CDELT1',spc2bsmth[3],header_comment='Smoothing-CDELT1')
	Header_Get_Add(spc_opf_smt,'CD1_1' ,spc2bsmth[3],header_comment='Smoothing-CD1_1')
	Header_Get_Add(spc_opf_smt,'CD2_2' ,spc2bsmth[3],header_comment='Smoothing-CD2_2')

	Header_Copy(spc_opf_smt,spc_ipf_smt,'h_s_c')
	Header_Copy(spc_opf_smt,spc_ipf_smt,'h_s_0')
	Header_History_Step(spc_ipf_smt,spc_opf_smt)
	return spc_opf_smt

def Spectra_Cont_GetVal(spc_2b_wght,*args,**kwargs):
	gcv_lmbd_i = kwargs.get('gcv_lmbd_i',1400)
	gcv_lmbd_f = kwargs.get('gcv_lmbd_f',1500)
	x_type     = kwargs.get('x_type','lambda')

	dest_dir_fr      = par_frg_dir + '/' + str(rad_sep[0][0]) + '-' + str(rad_sep[0][-1]) + '/'
	dest_dir_bk      = par_bkg_dir + '/' + str(rad_sep[0][0]) + '-' + str(rad_sep[0][-1]) + '/'

	try:
		specfile_step = fts_rd1_res + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = fts_rd2_res + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = fts_rd3_res + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = fts_rd4_res + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = fts_rd1_lst + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = fts_rd2_lst + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = fts_rd3_lst + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = fts_rd4_lst + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = str_rd1_stk + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = str_rd2_stk + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = str_rd3_stk + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = str_rd4_stk + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = spc_2b_wght 
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = res_stk_res + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = ind_stk_res + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = dest_dir_fr + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = dest_dir_bk + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	try:
		specfile_step = str_bst_stk + spc_2b_wght + '.fits'
		values = Spectra_x_y_Get(specfile_step,gcv_lmbd_i,gcv_lmbd_f,x_type=x_type)
	except:
		pass
	return values,len(values)

def Spectra_Stats(spc_stt_ifn,*args,**kwargs):
	spc_stt_lmbd_i = kwargs.get('spc_stt_lmbd_i',1200)
	spc_stt_lmbd_f = kwargs.get('spc_stt_lmbd_f',1750)
	x_type         = kwargs.get('x_type','lambda')

	spc_val_stt = Spectra_Cont_GetVal(spc_stt_ifn,gcv_lmbd_i=spc_stt_lmbd_i,gcv_lmbd_f=spc_stt_lmbd_f,x_type=x_type)
	SPC_NUM  = spc_val_stt[1]
	SPC_SUM  = bn.nansum(spc_val_stt[0])
	SPC_AVG  = bn.nanmean(spc_val_stt[0])
	SPC_MED  = bn.nanmedian(spc_val_stt[0])
	SPC_STD  = bn.nanstd(spc_val_stt[0])
	SPC_VAR  = bn.nanvar(spc_val_stt[0])
	SPC_MIN  = bn.nanmax(spc_val_stt[0])
	SPC_MAX  = bn.nanmin(spc_val_stt[0])
	SPC_SNR  = SNR(spc_val_stt[0])
	SNR_SNR  = float(SPC_SNR[0])
	SNR_SGN  = float(SPC_SNR[1])
	SNR_NSE  = float(SPC_SNR[2])

	if SNR_NSE == 0:
		SNR_SNR = 0
	else:
		pass
	if np.isnan(SPC_SUM) == True:
		SPC_SUM = 11111
	else:
		pass
	if np.isnan(SPC_AVG) == True:
		SPC_AVG = 22222
	else:
		pass
	if np.isnan(SPC_MED) == True:
		SPC_MED = 33333
	else:
		pass
	if np.isnan(SPC_STD) == True:
		SPC_STD = 44444
	else:
		pass
	if np.isnan(SPC_VAR) == True:
		SPC_VAR = 55555
	else:
		pass
	if np.isnan(SPC_MIN) == True:
		SPC_MIN = 66666
	else:
		pass
	if np.isnan(SPC_MAX) == True:
		SPC_MAX = 77777
	else:
		pass
	if np.isnan(SNR_SGN) == True:
		SNR_SGN = 99999
		SNR_NSE = 99999
		SNR_SNR = 99999
	else:
		pass
	if np.isnan(SNR_NSE) == True:
		SNR_SGN = 88888
		SNR_NSE = 88888
		SNR_SNR = 88888
	else:
		pass
	Header_Get_Add(spc_stt_ifn,'SPC_NUM',SPC_NUM       ,header_comment='Spec Stats Number')
	Header_Get_Add(spc_stt_ifn,'SPC_BOT',spc_stt_lmbd_i,header_comment='Spec Stats lower limit (A)')
	Header_Get_Add(spc_stt_ifn,'SPC_TOP',spc_stt_lmbd_f,header_comment='Spec Stats upper limit (A)')
	Header_Get_Add(spc_stt_ifn,'SPC_SUM',SPC_SUM       ,header_comment='Spec Stats sum')
	Header_Get_Add(spc_stt_ifn,'SPC_AVG',SPC_AVG       ,header_comment='Spec Stats mean')
	Header_Get_Add(spc_stt_ifn,'SPC_MED',SPC_MED       ,header_comment='Spec Stats median')
	Header_Get_Add(spc_stt_ifn,'SPC_STD',SPC_STD       ,header_comment='Spec Stats standard dev')
	Header_Get_Add(spc_stt_ifn,'SPC_VAR',SPC_VAR       ,header_comment='Spec Stats variance')
	Header_Get_Add(spc_stt_ifn,'SPC_MIN',SPC_MIN       ,header_comment='Spec Stats min')
	Header_Get_Add(spc_stt_ifn,'SPC_MAX',SPC_MAX       ,header_comment='Spec Stats max')
	Header_Get_Add(spc_stt_ifn,'SNR_SNR',SNR_SNR       ,header_comment='Spec Stats SNR Stoehr+08 snr')
	Header_Get_Add(spc_stt_ifn,'SNR_SGN',SNR_SGN       ,header_comment='Spec Stats SNR Stoehr+08 signal (median)')
	Header_Get_Add(spc_stt_ifn,'SNR_NSE',SNR_NSE       ,header_comment='Spec Stats SNR Stoehr+08 noise')

	spc_val_stt = Spectra_Cont_GetVal(spc_stt_ifn,gcv_lmbd_i=spc_stt_lmbd_i,gcv_lmbd_f=spc_stt_lmbd_f,x_type='all')
	SPC_NUA  = spc_val_stt[1]
	SPC_SUA  = bn.nansum(spc_val_stt[0])
	SPC_AVA  = bn.nanmean(spc_val_stt[0])
	SPC_MEA  = bn.nanmedian(spc_val_stt[0])
	SPC_STA  = bn.nanstd(spc_val_stt[0])
	SPC_VAA  = bn.nanvar(spc_val_stt[0])
	SPC_MIA  = bn.nanmax(spc_val_stt[0])
	SPC_MAA  = bn.nanmin(spc_val_stt[0])
	SPC_SNA  = SNR(spc_val_stt[0])
	SNR_SNA  = float(SPC_SNA[0])
	SNR_SGA  = float(SPC_SNA[1])
	SNR_NSA  = float(SPC_SNA[2])
	if SNR_NSA == 0:
		SNR_SNA = 0
	else:
		pass
	if np.isnan(SPC_SUA) == True:
		SPC_SUA = 11111
	else:
		pass
	if np.isnan(SPC_AVA) == True:
		SPC_AVA = 22222
	else:
		pass
	if np.isnan(SPC_MEA) == True:
		SPC_MEA = 33333
	else:
		pass
	if np.isnan(SPC_STA) == True:
		SPC_STA = 44444
	else:
		pass
	if np.isnan(SPC_VAA) == True:
		SPC_VAA = 55555
	else:
		pass
	if np.isnan(SPC_MIA) == True:
		SPC_MIA = 66666
	else:
		pass
	if np.isnan(SPC_MAA) == True:
		SPC_MAA = 77777
	else:
		pass
	if np.isnan(SNR_SGA) == True:
		SNR_SGA = 99999
		SNR_NSA = 99999
		SNR_SNA = 99999
	else:
		pass
	if np.isnan(SNR_NSA) == True:
		SNR_SGA = 88888
		SNR_NSA = 88888
		SNR_SNA = 88888
	else:
		pass
		
	Header_Get_Add(spc_stt_ifn,'SPC_NUA',SPC_NUA       ,header_comment='Spec Stats Number ALL')
	Header_Get_Add(spc_stt_ifn,'SPC_BOA','ALL'         ,header_comment='Spec Stats lower limit ALL (A)')
	Header_Get_Add(spc_stt_ifn,'SPC_TOA','ALL'         ,header_comment='Spec Stats upper limit ALL (A)')
	Header_Get_Add(spc_stt_ifn,'SPC_SUA',SPC_SUA       ,header_comment='Spec Stats sum ALL')
	Header_Get_Add(spc_stt_ifn,'SPC_AVA',SPC_AVA       ,header_comment='Spec Stats mean ALL')
	Header_Get_Add(spc_stt_ifn,'SPC_MEA',SPC_MEA       ,header_comment='Spec Stats median ALL')
	Header_Get_Add(spc_stt_ifn,'SPC_STA',SPC_STA       ,header_comment='Spec Stats standard dev ALL')
	Header_Get_Add(spc_stt_ifn,'SPC_VAA',SPC_VAA       ,header_comment='Spec Stats variance ALL')
	Header_Get_Add(spc_stt_ifn,'SPC_MIA',SPC_MIA       ,header_comment='Spec Stats min ALL')
	Header_Get_Add(spc_stt_ifn,'SPC_MAA',SPC_MAA       ,header_comment='Spec Stats max ALL')
	Header_Get_Add(spc_stt_ifn,'SNR_SNA',SNR_SNA       ,header_comment='Spec Stats SNR Stoehr+08 snr ALL')
	Header_Get_Add(spc_stt_ifn,'SNR_SGA',SNR_SGA       ,header_comment='Spec Stats SNR Stoehr+08 signal (median) ALL')
	Header_Get_Add(spc_stt_ifn,'SNR_NSA',SNR_NSA       ,header_comment='Spec Stats SNR Stoehr+08 noise ALL')

def get_last_lambda_item(spc_ifn_lastlambdaitem,*args, **kwargs):
	info = Spectra_x_y(spc_ifn_lastlambdaitem)
	return info[0][-1]
####Fnc_Stk_Spc####

