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

from termcolor import colored
from Fnc_Stk_Dir import *

#from Fnc_Stk_Spc import *
#from Fnc_Stk_Utl import *

####Fnc_Stk_Tbl####
def readtable_fg_bg_glx(table_name,format_tbl,*args, **kwargs):
	bs_func     = kwargs.get('bs_func',None)
	ref_cat_2br = kwargs.get('ref_cat_2br',False)
	#Random_Vars = kwargs.get('Random_Vars',False)

	if 'BS_MST' in table_name:
		print
		print colored('Reading table: ','green')
		print colored(table_name,'green')
		print
		print colored('BS_MST table!','yellow')
		print
		print colored('Function: ' + bs_func,'yellow')
		print colored('Reference catalogue table: ' + str(ref_cat_2br),'yellow')
		print
		print format_tbl
		ftbl = astropy.table.Table.read(table_name, format=format_tbl)
		if bs_func == '_med':
			print '1: ' + bs_func
			c1  = ftbl['bs_spc_file_med']
		elif bs_func == '_avg':
			print '2: ' + bs_func
			c1  = ftbl['bs_spc_file_avg']
		elif bs_func == '_avw':
			print '3: ' + bs_func
			c1  = ftbl['bs_spc_file_avw']
		elif bs_func == '_med-c':
			print '4: ' + bs_func
			c1  = ftbl['bs_spc_file_med-c']
		elif bs_func == '_avg-c':
			print '5: ' + bs_func
			c1  = ftbl['bs_spc_file_avg-c']
		elif bs_func == '_avw-c':
			print '6: ' + bs_func
			c1  = ftbl['bs_spc_file_avw-c']
		elif bs_func == '_avw-c-smt':
			print '7: ' + bs_func
			c1  = ftbl['bs_spc_file_avw-c-s']
		elif bs_func == '_avg-c-smt':
			print '8: ' + bs_func
			c1  = ftbl['bs_spc_file_avg-c-s']
		elif bs_func == '_med-c-smt':
			print '9: ' + bs_func
			c1  = ftbl['bs_spc_file_med-c-s']
		elif bs_func == None:
			print 'Function does not exist!'
			print bs_func
			quit()
		else:
			print 'Function does not exist!'
			print bs_func
			quit()			
		return(ftbl,c1)
	elif '-PRP-' in table_name:
		print
		print colored('Reading table: ','green')
		print colored(table_name,'green')
		print 'PRP'
		print colored(ref_cat_2br,'green')
		print
		ftbl = astropy.table.Table.read(table_name, format=format_tbl)
		c1   = ftbl['id_F']
		c2   = ftbl['z_F']
		c3   = ftbl['zf_F']
		c4   = ftbl['id_B']
		c5   = ftbl['z_B']
		c6   = ftbl['zf_B']
		c7   = ftbl['DELTAZ']
		c8   = ftbl['SEP_arcsec']
		c9   = ftbl['arcsec/kpc']
		c10  = ftbl['SEP_kpc']
		c11  = ftbl['spc_f_F']
		c12  = ftbl['spc_f_n_F']
		c13  = ftbl['spc_f_B']
		c14  = ftbl['spc_f_n_B']

		c15  = ftbl['mass_B']
		c16  = ftbl['Age_B']
		c17  = ftbl['SFR_B']
		c18  = ftbl['sSFR_B']
		c19  = ftbl['Lnuv_B']
		c20  = ftbl['mass_F']
		c21  = ftbl['Age_F']
		c22  = ftbl['SFR_F']
		c23  = ftbl['sSFR_F']
		c24  = ftbl['Lnuv_F']
		c25  = ftbl['magi_F']

		return(ftbl,
				c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
				c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
				c21,c22,c23,c24,c25)
	elif '_PRP_MRP' in table_name and ref_cat_2br == False:
		print
		print colored('Reading table: ','green')
		print colored(table_name,'green')
		print '_PRP_MRP'
		print colored(ref_cat_2br,'green')
		print
		ftbl = astropy.table.Table.read(table_name, format=format_tbl)
		c1   = ftbl['id_F']
		c2   = ftbl['z_F']
		c3   = ftbl['zf_F']
		c4   = ftbl['id_B']
		c5   = ftbl['z_B']
		c6   = ftbl['zf_B']
		c7   = ftbl['DELTAZ']
		c8   = ftbl['SEP_arcsec']
		c9   = ftbl['arcsec/kpc']
		c10  = ftbl['SEP_kpc']
		c11  = ftbl['spc_f_F']
		c12  = ftbl['spc_f_n_F']
		c13  = ftbl['spc_f_B']
		c14  = ftbl['spc_f_n_B']

		c15  = ftbl['mass_B']
		c16  = ftbl['Age_B']
		c17  = ftbl['SFR_B']
		c18  = ftbl['sSFR_B']
		c19  = ftbl['Lnuv_B']
		c20  = ftbl['mass_F']
		c21  = ftbl['Age_F']
		c22  = ftbl['SFR_F']
		c23  = ftbl['sSFR_F']
		c24  = ftbl['Lnuv_F']
		c25  = ftbl['magi_F']

		c26  = abs(ftbl['PHI'])
		c27  = ftbl['PHI_ABS']
		c28  = ftbl['re_B']
		c29  = ftbl['n_B']
		c30  = ftbl['q_B']
		c31  = ftbl['re_F_kpc']
		c32  = ftbl['n_F']
		c33  = ftbl['q_F']

		PHI_slice_new = []
		for j,angle_item in enumerate(c26):
			#print angle_item
			if angle_item>90:
				#print j,angle_item,180-angle_item
				PHI_slice_new.append(abs((180-angle_item)-90))
			else:
				PHI_slice_new.append(abs(angle_item-90))
		c26 = np.asarray(PHI_slice_new)
		add_phi_crc_col = kwargs.get('add_phi_crc_col',False)
		if add_phi_crc_col == True:
			print
			print colored('Adding PHI Corrected Column!','yellow')
			print
			ftbl_Mdf					= astropy.table.Table()

			ftbl_Mdf['id_F']			= c1 
			ftbl_Mdf['z_F']				= c2 
			ftbl_Mdf['zf_F']			= c3 
			ftbl_Mdf['id_B']			= c4 
			ftbl_Mdf['z_B']				= c5 
			ftbl_Mdf['zf_B']			= c6 
			ftbl_Mdf['DELTAZ']			= c7 
			ftbl_Mdf['SEP_arcsec']		= c8 
			ftbl_Mdf['arcsec/kpc']		= c9 
			ftbl_Mdf['SEP_kpc']			= c10
			ftbl_Mdf['spc_f_F']			= c11
			ftbl_Mdf['spc_f_n_F']		= c12
			ftbl_Mdf['spc_f_B']			= c13
			ftbl_Mdf['spc_f_n_B']		= c14

			ftbl_Mdf['mass_B']			= c15
			ftbl_Mdf['Age_B']			= c16
			ftbl_Mdf['SFR_B']			= c17
			ftbl_Mdf['sSFR_B']			= c18
			ftbl_Mdf['Lnuv_B']			= c19
			ftbl_Mdf['mass_F']			= c20
			ftbl_Mdf['Age_F']			= c21
			ftbl_Mdf['SFR_F']			= c22
			ftbl_Mdf['sSFR_F']			= c23
			ftbl_Mdf['Lnuv_F']			= c24
			ftbl_Mdf['magi_F']			= c25

			ftbl_Mdf['PHI']				= c26
			ftbl_Mdf['PHI_ABS']			= c27
			ftbl_Mdf['PHI_CRC']			= PHI_slice_new
			ftbl_Mdf['re_B']			= c28
			ftbl_Mdf['n_B']				= c29
			ftbl_Mdf['q_B']				= c30
			ftbl_Mdf['re_F_kpc']		= c31
			ftbl_Mdf['n_F']				= c32
			ftbl_Mdf['q_F']				= c33

			ftbl_Mdf.write(table_name, format=tbl_format_opt,overwrite=True)#'ascii.fixed_width_two_line')	
			print
			print colored('Modified Table: ','green')
			print colored(table_name,'green')
			print
		else:
			pass
		return(ftbl,
				c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
				c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
				c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,
				c31,c32,c33)
	elif '-PRP.' in table_name and ref_cat_2br == True:
		print
		print colored('Reading table: ','green')
		print colored(table_name,'green')
		print 'PRP'
		print colored(ref_cat_2br,'green')
		print
		ftbl = astropy.table.Table.read(table_name, format=format_tbl)
		c1   = ftbl['ident']
		c2   = ftbl['z_spec']
		c3   = ftbl['zflags']
		c4   = ftbl['alpha']
		c5   = ftbl['delta']
		c6   = ftbl['id_iau']
		c7   = ftbl['epoch']
		c8   = ftbl['spec1d']
		c9   = ftbl['spec1dnoise']
		c10  = ftbl['pointing']
		c11  = ftbl['quad']
		c12  = ftbl['slit']
		c13  = ftbl['obj']
		c14  = ftbl['obj']

		c15  = ftbl['mass']
		c16  = ftbl['Age']
		c17  = ftbl['sfr']
		c18  = ftbl['ssfr']
		c19  = ftbl['L_nuv_rest']
		c20  = ftbl['mass']
		c21  = ftbl['Age']
		c22  = ftbl['sfr']
		c23  = ftbl['ssfr']
		c24  = ftbl['L_nuv_rest']
		c25  = ftbl['magi']

		return(ftbl,
				c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
				c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
				c21,c22,c23,c24,c25)
	elif '_PRP_MRP.' in table_name and ref_cat_2br == True:
		print
		print colored('Reading table: ','green')
		print colored(table_name,'green')
		print '_PRP_MRP'
		print colored(ref_cat_2br,'green')
		print
		ftbl = astropy.table.Table.read(table_name, format=format_tbl)
		c1   = ftbl['ident_1']
		c2   = ftbl['z_spec']
		c3   = ftbl['zflags']
		c4   = ftbl['alpha']
		c5   = ftbl['delta']
		c6   = ftbl['id_iau']
		c7   = ftbl['epoch']
		c8   = ftbl['spec1d']
		c9   = ftbl['spec1dnoise']
		c10  = ftbl['pointing']
		c11  = ftbl['quad']
		c12  = ftbl['slit']
		c13  = ftbl['obj']
		c14  = ftbl['obj']

		c15  = ftbl['mass']
		c16  = ftbl['Age']
		c17  = ftbl['sfr']
		c18  = ftbl['ssfr']
		c19  = ftbl['L_nuv_rest']
		c20  = ftbl['mass']
		c21  = ftbl['Age']
		c22  = ftbl['sfr']
		c23  = ftbl['ssfr']
		c24  = ftbl['L_nuv_rest']
		c25  = ftbl['magi']

		c26  = c25#abs(ftbl['PHI'])
		c27  = c25#ftbl['PHI_ABS']
		c28  = ftbl['re']
		c29  = ftbl['n']
		c30  = ftbl['q']
		c31  = ftbl['re']
		c32  = ftbl['n']
		c33  = ftbl['q']

		PHI_slice_new = []
		for j,angle_item in enumerate(c26):
			#print angle_item
			if angle_item>90:
				#print j,angle_item,180-angle_item
				PHI_slice_new.append(abs((180-angle_item)-90))
			else:
				PHI_slice_new.append(abs(angle_item-90))
		c26 = np.asarray(PHI_slice_new)
		return(ftbl,
				c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
				c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
				c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,
				c31,c32,c33)
	else:
		print
		print colored('Reading table: ','green')
		print colored(table_name,'green')
		print colored(ref_cat_2br,'green')
		print
		ftbl = astropy.table.Table.read(table_name, format=format_tbl)
		c1  = ftbl['id_F']
		c2  = ftbl['z_F']
		c3  = ftbl['zf_F']
		c4  = ftbl['id_B']
		c5  = ftbl['z_B']
		c6  = ftbl['zf_B']
		c7  = ftbl['DELTAZ']
		c8  = ftbl['SEP_arcsec']
		c9  = ftbl['arcsec/kpc']
		c10 = ftbl['SEP_kpc']
		c11 = ftbl['spc_f_F']
		c12 = ftbl['spc_f_n_F']
		c13 = ftbl['spc_f_B']
		c14 = ftbl['spc_f_n_B']

		#if Random_Vars == False:
		 #pass
		#elif Random_Vars == 'z_F':
			#c2  = ftbl['RND_z_F']
			#c7  = ftbl['RND_DELTAZ']
		#elif Random_Vars == 'SEP_arcsec':
			#c8  = ftbl['RND_SEP_arcsec']
		#elif Random_Vars == 'Both':
			#c2  = ftbl['RND_z_F']
			#c7  = ftbl['RND_DELTAZ']
			#c8  = ftbl['RND_SEP_arcsec']
		return(ftbl,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)

def readtable_fg_bg_glx_cdf(table_name,format_tbl,*args, **kwargs):
	cdf2bobt_lne = kwargs.get('cdf2bobt_lne',None)
	#Random_Vars = kwargs.get('Random_Vars',False)
	print
	print colored('Reading table: ','green')
	print colored(table_name,'green')
	print
	ftbl = astropy.table.Table.read(table_name, format=format_tbl)
	if 'gaussM.csv' in table_name:
		try:
			c1    = ftbl[cdf2bobt_lne + '_CGL1LC']
			c2    = ftbl[cdf2bobt_lne + '_AGL1LC']
			c3    = ftbl[cdf2bobt_lne + '_SGL1LC']
			c4    = ftbl[cdf2bobt_lne + '_CLCMLC']
			c5    = ftbl[cdf2bobt_lne + '_ALCMLC']
			c6    = ftbl[cdf2bobt_lne + '_SLCMLC']
			c7    = ftbl[cdf2bobt_lne + '_CGL2LC']
			c8    = ftbl[cdf2bobt_lne + '_AGL2LC']
			c9    = ftbl[cdf2bobt_lne + '_SGL2LC']
			c10   = ftbl[cdf2bobt_lne + '_WLCML']
			c11   = ftbl[cdf2bobt_lne + '_ELCML']
			c12   = ftbl[cdf2bobt_lne + '_WGL1L']
			c13   = ftbl[cdf2bobt_lne + '_EGL1L']
			c14   = ftbl[cdf2bobt_lne + '_WGL2L']
			c15   = ftbl[cdf2bobt_lne + '_EGL2L']
			c16   = ftbl[cdf2bobt_lne + '_WPLPL']
			c17   = ftbl[cdf2bobt_lne + '_WEPLL']
			c18   = ftbl[cdf2bobt_lne + '_WTOTL']
			c19   = ftbl[cdf2bobt_lne + '_WETTL']
		except KeyError:
			print 'line 308'
			print'Fnc_Stk_TTbl.py'
			print
			quit()
		return(ftbl,
				c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
				c11,c12,c13,c14,c15,c16,c17,c18,c19
				)


	else:
		try:
			c1   = ftbl[cdf2bobt_lne + '_CGLC']
			c2   = ftbl[cdf2bobt_lne + '_AGLC']
			c3   = ftbl[cdf2bobt_lne + '_SGLC']
			c4   = ftbl[cdf2bobt_lne + '_FGLC']
			c5   = ftbl[cdf2bobt_lne + '_WGLC']
			c6   = ftbl[cdf2bobt_lne + '_EGLC']
			c7   = ftbl[cdf2bobt_lne + '_CGLEC']
			c8   = ftbl[cdf2bobt_lne + '_AGLEC']
			c9   = ftbl[cdf2bobt_lne + '_SGLEC']
			c10  = ftbl[cdf2bobt_lne + '_CH2GL']
			c11  = ftbl[cdf2bobt_lne + '_CHRGL']
		except KeyError:
			c1   = ftbl[cdf2bobt_lne + '_CGAC']
			c2   = ftbl[cdf2bobt_lne + '_AGAC']
			c3   = ftbl[cdf2bobt_lne + '_SGAC']
			c4   = ftbl[cdf2bobt_lne + '_FGAC']
			c5   = ftbl[cdf2bobt_lne + '_WGAC']
			c6   = ftbl[cdf2bobt_lne + '_EGAC']
			c7   = ftbl[cdf2bobt_lne + '_CGAEC']
			c8   = ftbl[cdf2bobt_lne + '_AGAEC']
			c9   = ftbl[cdf2bobt_lne + '_SGAEC']		
			c10  = ftbl[cdf2bobt_lne + '_CH2GA']
			c11  = ftbl[cdf2bobt_lne + '_CHRGA']			

		return(ftbl,
				c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
				c11)

def readtable_Lit(table_name,format_tbl,*args, **kwargs):
	#Random_Vars = kwargs.get('Random_Vars',False)
	print table_name
	print
	print colored('Reading Literature table: ','yellow')
	print colored(table_name,'yellow')
	print	
	ftbl = astropy.table.Table.read(table_name, format=format_tbl)
	return ftbl

def Confident_Intervals(table_name,format_tbl,*args, **kwargs):
	cdf2bobt_lne = kwargs.get('cdf2bobt_lne',None)
	cdf_2b_rtrn  = kwargs.get('cdf_2b_rtrn','1sgm')
	rmv_uft_spc  = kwargs.get('rmv_uft_spc',False)
	ems_lne_ctb  = kwargs.get('ems_lne_ctb',False)
	cdf_lne      = readtable_fg_bg_glx_cdf(table_name,format_tbl,**kwargs)
	if 'gaussM.csv' in table_name and ems_lne_ctb == False:
		ctr_cdf      = cdf_lne[1]
		amp_cdf      = cdf_lne[2]
		sgm_cdf      = cdf_lne[3]

		ctr_cdf      = cdf_lne[7]
		amp_cdf      = cdf_lne[8]
		sgm_cdf      = cdf_lne[9]

		ewd1_cdf     = cdf_lne[12]
		ewd2_cdf     = cdf_lne[14]
		ewd_cdf      = cdf_lne[10]		
	elif 'gaussM.csv' in table_name and ems_lne_ctb == True:
		ctr_cdf      = cdf_lne[1]
		amp_cdf      = cdf_lne[2]
		sgm_cdf      = cdf_lne[3]

		ctr_cdf      = cdf_lne[7]
		amp_cdf      = cdf_lne[8]
		sgm_cdf      = cdf_lne[9]

		ewd1_cdf     = cdf_lne[12]
		ewd2_cdf     = cdf_lne[14]
		ewd_cdf      = ewd1_cdf + ewd2_cdf #cdf_lne[10]		
	else:
		ctr_cdf      = cdf_lne[1]
		amp_cdf      = cdf_lne[2]
		sgm_cdf      = cdf_lne[3]
		fwhm_cdf     = cdf_lne[4]
		ewd_cdf      = cdf_lne[5]
		ewd_cdf_err  = cdf_lne[6]
		ctr_cdf_err  = cdf_lne[7]
		amp_cdf_err  = cdf_lne[8]
		sgm_cdf_err  = cdf_lne[9]
		ch2_cdf      = cdf_lne[10]
		ch2r_cdf     = cdf_lne[11]

    #bn.nansum(np.array(img_stat)     , axis =0)#np.nansum(np.array(img_stat)   , axis=0)#
	#bn.nanmean(np.array(img_stat)    , axis=0)#np.nanmean(np.array(img_stat)  , axis=0) #
	#bn.nanmedian(np.array(img_stat)  , axis=0)#np.nanmedian(np.array(img_stat), axis=0) #
	#bn.nanstd(np.array(img_stat)     , axis=0)#np.nanstd(np.array(img_stat)   , axis=0) #
	if rmv_uft_spc == True:
		indx_msk_uf = (np.where(np.asarray(ewd_cdf)==999999.99999)[0])
		indx_msk_uf = np.asarray(indx_msk_uf)
		ewd_cdf     = np.delete(ewd_cdf,indx_msk_uf)
		print
		print colored('Removing unfitted values from BS distribution!','yellow')
		print colored('N: '+str(len(indx_msk_uf)),'yellow')
		print
	else:
		pass

	ewd_cdf_1sl = np.nanpercentile(np.array(ewd_cdf), 15.9, axis=0)
	ewd_cdf_1sh = np.nanpercentile(np.array(ewd_cdf), 84.1, axis=0)

	ewd_cdf_2sl = np.nanpercentile(np.array(ewd_cdf), 2.30, axis=0)
	ewd_cdf_2sh = np.nanpercentile(np.array(ewd_cdf), 97.7, axis=0)

	ewd_cdf_3sl = np.nanpercentile(np.array(ewd_cdf), 0.20, axis=0)
	ewd_cdf_3sh = np.nanpercentile(np.array(ewd_cdf), 99.8, axis=0)

	ewd_cdf_p25 = np.nanpercentile(np.array(ewd_cdf), 25, axis=0)
	ewd_cdf_p75 = np.nanpercentile(np.array(ewd_cdf), 75, axis=0)

	if cdf_2b_rtrn == '1sgm':
		print
		print colored('1-Sigma errors','yellow')
		print colored(str(ewd_cdf_1sl)+', '+str(ewd_cdf_1sh),'magenta')
		return ewd_cdf_1sl,ewd_cdf_1sh
	elif cdf_2b_rtrn == '2sgm':
		print
		print colored('2-Sigma errors','yellow')
		print colored(str(ewd_cdf_2sl)+', '+str(ewd_cdf_2sh),'magenta')
		return ewd_cdf_2sl,ewd_cdf_2sh
	elif cdf_2b_rtrn == '3sgm':
		print
		print colored('3-Sigma errors','yellow')
		print colored(str(ewd_cdf_3sl)+', '+str(ewd_cdf_3sh),'magenta')
		return ewd_cdf_3sl,ewd_cdf_3sh
	else:
		print
		print colored('1-Sigma errors','yellow')
		print colored(str(ewd_cdf_1sl)+', '+str(ewd_cdf_1sh),'magenta')
		return ewd_cdf_1sl,ewd_cdf_1sh

def stats_table(tbl_stt_in,format_tbl,*args, **kwargs):
	df = pd.read_csv(tbl_stt_in)
	col_nms = df.head()

	ftbl_stt = (df.describe(percentiles=[.159,.841,.0230,.977,0.0020,.998], include='all'))

	tbl_stt_out = stt_dir_res + (str(tbl_stt_in.split('.csv',1)[0])).split('/')[-1] + '-stt.csv'
	ftbl_stt.to_csv(tbl_stt_out)

	ftbl_stt_c    = astropy.table.Table.read(tbl_stt_out, format='csv')
	tbl_stt_out_c = stt_dir_res + (str(tbl_stt_in.split('.csv',1)[0])).split('/')[-1] + '-stt.dat'
	ftbl_stt_c.write(tbl_stt_out_c,  format='ascii.fixed_width_two_line',overwrite=True)

	print colored('Stat file: '+str(tbl_stt_out_c),'green')

	return tbl_stt_out_c

####Fnc_Stk_Utl####