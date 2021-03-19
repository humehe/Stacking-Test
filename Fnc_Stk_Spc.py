import sys, os
import bottleneck as bn
#import astropy
from astropy import convolution as apcvl

#import scipy
#from scipy import stats
from pysynphot import observation
from pysynphot import spectrum

import pyraf
from pyraf import iraf
from pyraf.iraf import continuum as pyfcnt

from Fnc_Stk_Fts import *
from Fnc_Stk_Mth import *
from Fnc_Stk_Dir import *

import Lines_Dictionary
LINES = Lines_Dictionary.LINES_STK #FOR MASKING##

####Fnc_Stk_Spc####
def Shift_Spectra(spec_ifn,redshift_ref,id_ref,*args, **kwargs):

	shft_ref   = kwargs.get('shft_ref', 'rf')
	redshift_0 = kwargs.get('redshift_0',redshift_ref)
	id_0       = kwargs.get('id_0',id_ref)

	sel_pre_shf = kwargs.get('sel_pre_shf',True)
	sel_pre_cnt = kwargs.get('sel_pre_cnt',True)
	sel_pre_msk = kwargs.get('sel_pre_msk',False)

	Header_Get_Add(spec_ifn,'CD1_1'  ,5.355)
	Header_Get_Add(spec_ifn,'CD2_2'  ,5.355)
	Header_Get_Add(spec_ifn,'CDELT1' ,5.355)

	lmbstp0fg  = Header_Get_Add(spec_ifn,'CDELT1',5.355)
	lambda0fg  = Header_Get_Add(spec_ifn,'CRVAL1',3510.71)

	pre_shft_spec = Spectra_x_y(spec_ifn)
	pre_shft_lmbd = pre_shft_spec[0]

	spec_ofn   = str(spec_ifn.split('.fits',1)[0])   +'-s2-' + str(shft_ref) + '.fits'
	hdulist_fg = fits.open(spec_ifn)
	intens_fg  = hdulist_fg[0].data

	lambda_fg2fg   = np.around(lambdashifted(lambda0fg ,redshift_ref),decimals=2)
	lmbstp_fg2fg   = np.around(lambdashifted(lmbstp0fg ,redshift_ref),decimals=3)

	os.system('cp ' + spec_ifn   + ' ' + spec_ofn)
	Header_History_Step(spec_ifn,spec_ofn)

	Header_Updt(spec_ofn,'CRVAL1',lambda_fg2fg)
	Header_Updt(spec_ofn,'CDELT1',lmbstp_fg2fg)
	Header_Updt(spec_ofn,'CD2_2' ,lmbstp_fg2fg)
	Header_Updt(spec_ofn,'CD1_1' ,lmbstp_fg2fg)

	Header_Get_Add(spec_ofn,'ID_0'   ,int(id_0)   ,header_comment='ID')
	Header_Get_Add(spec_ofn,'ID_REF' ,int(id_ref) ,header_comment='ID Reference')
 	Header_Get_Add(spec_ofn,'Z_0'    ,redshift_0  ,header_comment='Redshift')
	Header_Get_Add(spec_ofn,'Z_REF'  ,redshift_ref,header_comment='Redshift Reference')

	return spec_ofn

def Shift_Spectra_Cat(spc_cat_itn,*args, **kwargs):
	from Fnc_Stk_Spc import Spectra_Cont_IRAF
	sel_pre_shf = kwargs.get('sel_pre_shf',True)
	sel_pre_cnt = kwargs.get('sel_pre_cnt',True)
	sel_pre_msk = kwargs.get('sel_pre_msk',False)
	sfx_tbl_otp = kwargs.get('sfx_tbl_otp','')
	verbose     = kwargs.get('verbose'    ,False)
	

	print
	print colored('Spectra shifted      :'+str(sel_pre_shf),'yellow')
	print colored('Spectra cont fitted  :'+str(sel_pre_cnt),'yellow')
	print colored('Spectra masked       :'+str(sel_pre_msk),'yellow')
	print

	sel_cnt_typ     = kwargs.get('sel_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	sel_cnt_lns     = kwargs.get('sel_cnt_lns'     ,'*')       # Image lines to be fit
	sel_cnt_fnc     = kwargs.get('sel_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	sel_cnt_ord     = kwargs.get('sel_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	sel_cnt_ovr     = kwargs.get('sel_cnt_ovr'     ,'yes')     # Override previous norm spec
	sel_cnt_rpl     = kwargs.get('sel_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	sel_cnt_lrj     = kwargs.get('sel_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	sel_cnt_hrj     = kwargs.get('sel_cnt_hrj'     ,3)         # High rejection in sigma of fit

	sel_msk_type    = kwargs.get('sel_msk_type'    ,'NaN')
	sel_msk_cte_val = kwargs.get('sel_msk_cte_val' ,1)
	sel_msk_abs_lne = kwargs.get('sel_msk_abs_lne' ,False)
	sel_msk_spc_rgn = kwargs.get('sel_msk_spc_rgn' ,False)
	sel_msk_lmb_min = kwargs.get('sel_msk_lmb_min' ,500)
	sel_msk_lmb_max = kwargs.get('sel_msk_lmb_max' ,1210)

	print
	print colored('Shifting spectra files.','yellow')
	print colored('Reading table : ' + spc_cat_itn,'green')
	print

	cat          = readtable_fg_bg_glx(spc_cat_itn,tbl_format_ipt,**kwargs)

	if 'PRP' in spc_cat_itn:
		tbl_cat_M        = cat[0]
		ident_fg_glx_M   = cat[5]
		z_fg_glx_M       = cat[6]
		z_f_fg_glx_M     = cat[7]
		magi_fg_glx_M    = cat[8]
		deltaz_fg_glx_M  = cat[9]
		sepas_fg_glx_M   = cat[10]
	else:
		tbl_cat_M        = cat[0]
		ident_fg_glx_M   = cat[1]
		z_fg_glx_M       = cat[2]
		z_f_fg_glx_M     = cat[3]
		magi_fg_glx_M    = cat[4]
		deltaz_fg_glx_M  = cat[9]
		sepas_fg_glx_M   = cat[10]

	tbl_cat_M.sort('z_F')
	name_sorted = (spc_cat_itn.split('.'+tbl_format_ipt,1)[0])+'-sorted'+'.'+tbl_format_ipt
	print
	print colored('Input Catalogue sorted by z_F','yellow')
	print colored(name_sorted,'green')
	print
	tbl_cat_M.write(name_sorted, format=tbl_format_opt, overwrite=True)

	widgets  = ['Evaluating pairs for '+ str(len(ident_fg_glx_M)) + ' foreground galaxies: ',
	Percentage(), ' ', Bar(marker='>',left='[',right=']'),
	' ', ETA(), ' ', FileTransferSpeed()]
	pbar1    = ProgressBar(widgets=widgets, maxval=len(ident_fg_glx_M))
	pbar1.start()

	for foregorund_galaxy in range(len(ident_fg_glx_M)):
		print
		pbar1.update(foregorund_galaxy)

		id_fg           = ident_fg_glx_M[foregorund_galaxy]
		z_fg            = z_fg_glx_M[foregorund_galaxy]
		z_f_fg          = z_f_fg_glx_M[foregorund_galaxy]
		magi_fg         = magi_fg_glx_M[foregorund_galaxy]
		
		if verbose == True:
			print ''
			print 'Foreground galaxy: ',foregorund_galaxy+1,'/',len(ident_fg_glx_M)
			print
			print ident_fg_glx_M[foregorund_galaxy],z_fg_glx_M[foregorund_galaxy],z_f_fg_glx_M[foregorund_galaxy],magi_fg_glx_M[foregorund_galaxy]
			print id_fg,z_fg,z_f_fg,magi_fg
		elif verbose == False:
			pass

		print
		print colored(par_dir_res+'/' +str(id_fg)+'/'+str(id_fg)+'_pairs' + sfx_tbl_otp + tbl_ext_ipt,'magenta')
		print 
		fg_bg_glx      = readtable_fg_bg_glx(par_dir_res+'/' +str(id_fg)+'/'+str(id_fg)+'_pairs' + sfx_tbl_otp + tbl_ext_ipt,tbl_format_ipt)
		
		tbl_fg_bg_glx  = fg_bg_glx[0]

		ident_bg_glx_p = fg_bg_glx[1]
		z_bg_glx       = fg_bg_glx[2]
		z_f_bg_glx     = fg_bg_glx[3]
		magi_bg_glx    = fg_bg_glx[4]

		ident_fg_glx_p = fg_bg_glx[5]
		z_fg_glx       = fg_bg_glx[6]
		z_f_fg_glx     = fg_bg_glx[7]
		magi_fg_glx    = fg_bg_glx[8]

		deltaz_bg_glx  = fg_bg_glx[9]
		sepas_bg_glx   = fg_bg_glx[10]
		sepkpc_bg_glx  = fg_bg_glx[12]

		spc_f_fg_glx   = fg_bg_glx[15]
		spc_f_n_fg_glx = fg_bg_glx[16]
		spc_f_bg_glx   = fg_bg_glx[17]
		spc_f_n_bg_glx = fg_bg_glx[18]

		if (len(ident_bg_glx_p) >= 2):
			subdir        = par_dir_res + '/' + str(id_fg) +'/'
			dest_dir_fr   = par_frg_dir + str(rad_sep[0][0]) + '-' + str(rad_sep[0][-1]) + '/'

			specfile_fg   = subdir + str(spc_f_fg_glx[0])
			specfile_n_fg = subdir + str(spc_f_n_fg_glx[0])

			if sel_pre_cnt == True:
				spc_ipt_cnt = specfile_fg
			elif sel_pre_cnt == False:
				pass

			if sel_pre_msk == True and sel_pre_cnt == True:
				spc_ipt_msk = str(specfile_fg.split('.fits',1)[0]) + '-c.fits'
 	 		elif sel_pre_msk == True and sel_pre_cnt == False:
				spc_ipt_msk = specfile_fg
			else:
				pass

 			if sel_pre_shf == True and sel_pre_msk == True and sel_pre_cnt == True:
				spc_ipt_shf     = str(specfile_fg.split('.fits',1)[0]) + '-c-m.fits'
				spc_ipt_shf_fit = str(specfile_fg.split('.fits',1)[0]) + '-c-f.fits'

			elif sel_pre_shf == True and sel_pre_msk == True and sel_pre_cnt == False:
				spc_ipt_shf     = str(specfile_fg.split('.fits',1)[0]) + '-m.fits'
				spc_ipt_shf_fit = spc_ipt_shf

			elif sel_pre_shf == True and sel_pre_msk == False and sel_pre_cnt == True:
				spc_ipt_shf     = str(specfile_fg.split('.fits',1)[0]) + '-c.fits'
				spc_ipt_shf_fit = str(specfile_fg.split('.fits',1)[0]) + '-c-f.fits'

			elif sel_pre_shf == True and sel_pre_msk == False and sel_pre_cnt == False:
				spc_ipt_shf     = specfile_fg
				spc_ipt_shf_fit = spc_ipt_shf 

			else:
				pass

			if sel_pre_cnt == True:
				spc_opt_cnt  = Spectra_Cont_IRAF(spc_ipt_cnt,subdir + 'log_cont_' + str(id_fg),
								Cont_type_IRAF     = sel_cnt_typ,Cont_lines_IRAF    = sel_cnt_lns,
								Cont_funct_IRAF    = sel_cnt_fnc,Cont_order_IRAF    = sel_cnt_ord,
								Cont_override_IRAF = sel_cnt_ovr,Cont_replace_IRAF  = sel_cnt_rpl,
								Cont_low_rej_IRAF  = sel_cnt_lrj,Cont_high_rej_IRAF = sel_cnt_hrj)

				Header_Add(spc_opt_cnt[0],'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				Header_Add(spc_opt_cnt[1],'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				Header_Add(specfile_fg   ,'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				os.system ('cp ' + str(spc_opt_cnt[0]) + ' ' + dest_dir_fr)
				os.system ('cp ' + str(spc_opt_cnt[1]) + ' ' + dest_dir_fr)
				os.system ('cp ' + specfile_fg         + ' ' + dest_dir_fr)

			elif sel_pre_cnt == True and spectra_noise == True:
				Header_Add(specfile_n_fg ,'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				os.system ('cp ' + specfile_n_fg       + ' ' + dest_dir_fr)
			elif sel_pre_cnt == False:
				pass

			if sel_pre_msk == True:
				spc_opt_msk = Spectra_Masking(spc_ipt_msk,msk_typ=sel_msk_type,rshft_corr=z_fg,
								rshft_corr_direct=False,msk_abs_lne=sel_msk_abs_lne,
								msk_spc_rgn=sel_msk_spc_rgn,msk_lmb_min=sel_msk_lmb_min,msk_lmb_max=sel_msk_lmb_max)
				Header_Add(spc_opt_msk  ,'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				Header_Add(specfile_fg  ,'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				os.system ('cp ' + spc_opt_msk	 + ' ' + dest_dir_fr)
				os.system ('cp ' + specfile_fg   + ' ' + dest_dir_fr)
			elif sel_pre_msk == True and spectra_noise == True:
				Header_Add(specfile_n_fg,'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				os.system ('cp ' + specfile_n_fg + ' ' + dest_dir_fr)
			elif sel_pre_msk == False:
				pass

			if sel_pre_shf == True:
				spc_opt_shf     = Shift_Spectra(spc_ipt_shf    ,z_fg,id_fg)
				spc_opt_shf_fit = Shift_Spectra(spc_ipt_shf_fit,z_fg,id_fg)

				Header_Add(spc_opt_shf    ,'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				Header_Add(spc_opt_shf    ,'SHT_FN0',str((spc_ipt_shf_fit.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Spec file pre-shift')
				Header_Add(spc_opt_shf    ,'SHT_FNS',str((spc_opt_shf_fit.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Spec file post-shift')


				Header_Add(spc_opt_shf_fit,'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				Header_Add(spc_opt_shf_fit,'SHT_FN0',str((spc_ipt_shf_fit.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Spec file pre-shift')
				Header_Add(spc_opt_shf_fit,'SHT_FNS',str((spc_opt_shf_fit.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Spec file post-shift')

				Header_Add(specfile_fg    ,'MAG_I',magi_fg,header_comment='i-band mag VUDS')

				os.system ('cp ' + spc_opt_shf     + ' ' + dest_dir_fr)
				os.system ('cp ' + spc_opt_shf_fit + ' ' + dest_dir_fr)
				os.system ('cp ' + spc_ipt_shf     + ' ' + dest_dir_fr)
				os.system ('cp ' + spc_ipt_shf_fit + ' ' + dest_dir_fr)

				os.system ('cp ' + specfile_fg     + ' ' + dest_dir_fr)
			elif sel_pre_shf == True and spectra_noise == True:
				specfile_n_fg_s = Shift_Spectra(specfile_n_fg  ,z_fg,id_fg)

				Header_Add(specfile_n_fg  ,'MAG_I',magi_fg,header_comment='i-band mag VUDS')
				Header_Add(specfile_n_fg_s,'MAG_I',magi_fg,header_comment='i-band mag VUDS')

				os.system ('cp ' + specfile_n_fg   + ' ' + dest_dir_fr)
				os.system ('cp ' + specfile_n_fg_s + ' ' + dest_dir_fr)	
			elif sel_pre_shf == False:	
				pass

			print
			widgets = ['Shifting background galaxies spectra: ' + str(len(ident_fg_glx_p)-1) + ', paired with foreground galaxy: '+ str(id_fg) + ' ', 
			Percentage(), ' ', Bar(marker='*',left='[',right=']'),
			' ', ETA(), ' ', FileTransferSpeed()]
			pbar    = ProgressBar(widgets=widgets, maxval=len(ident_bg_glx_p))
			pbar.start()
			
			for background_galaxy in range(1,len(ident_bg_glx_p)):
				pbar.update(background_galaxy)

				dest_dir_bk      = par_bkg_dir + str(rad_sep[0][0]) + '-' + str(rad_sep[0][-1]) + '/'

				specfile_bg      = subdir + str(spc_f_bg_glx[background_galaxy])
				specfile_n_bg    = subdir + str(spc_f_n_bg_glx[background_galaxy])

				id_bg           = ident_bg_glx_p[background_galaxy]
				z_bg            = z_bg_glx[background_galaxy]
				z_f_bg          = z_f_bg_glx[background_galaxy]
				magi_bg         = magi_bg_glx[background_galaxy]

				if verbose == True:
					print 
					print colored('Background galaxy to be shifted: ','yellow')
					print background_galaxy,ident_bg_glx_p[background_galaxy],z_bg,sepas_bg_glx[background_galaxy],deltaz_bg_glx[background_galaxy],magi_bg_glx[background_galaxy]
					print 
				elif verbose == False:
					pass

				if os.path.isfile(specfile_bg) or os.path.isfile(specfile_n_bg): 

					if sel_pre_cnt == True:
						spc_ipt_cnt = specfile_bg
					elif sel_pre_cnt == False:
						pass

					if sel_pre_msk == True and sel_pre_cnt == True:
						spc_ipt_msk = str(specfile_bg.split('.fits',1)[0]) + '-c.fits'
		 	 		elif sel_pre_msk == True and sel_pre_cnt == False:
						spc_ipt_msk = specfile_bg
					else:
						pass

		 			if sel_pre_shf == True and sel_pre_msk == True and sel_pre_cnt == True:
						spc_ipt_shf     = str(specfile_bg.split('.fits',1)[0]) + '-c-m.fits'
						spc_ipt_shf_fit = str(specfile_bg.split('.fits',1)[0]) + '-c-f.fits'

					elif sel_pre_shf == True and sel_pre_msk == True and sel_pre_cnt == False:
						spc_ipt_shf     = str(specfile_bg.split('.fits',1)[0]) + '-m.fits'
						spc_ipt_shf_fit = spc_ipt_shf

					elif sel_pre_shf == True and sel_pre_msk == False and sel_pre_cnt == True:
						spc_ipt_shf     = str(specfile_bg.split('.fits',1)[0]) + '-c.fits'
						spc_ipt_shf_fit = str(specfile_bg.split('.fits',1)[0]) + '-c-f.fits'

					elif sel_pre_shf == True and sel_pre_msk == False and sel_pre_cnt == False:
						spc_ipt_shf     = specfile_bg
						spc_ipt_shf_fit = spc_ipt_shf

					else:
						pass

					if sel_pre_cnt == True:
						spc_opt_cnt = Spectra_Cont_IRAF(spc_ipt_cnt,subdir + 'log_cont_' + str(id_fg),
										Cont_type_IRAF     = sel_cnt_typ,Cont_lines_IRAF    = sel_cnt_lns,
										Cont_funct_IRAF    = sel_cnt_fnc,Cont_order_IRAF    = sel_cnt_ord,
										Cont_override_IRAF = sel_cnt_ovr,Cont_replace_IRAF  = sel_cnt_rpl,
										Cont_low_rej_IRAF  = sel_cnt_lrj,Cont_high_rej_IRAF = sel_cnt_hrj)

						os.system ('cp ' + str(spc_opt_cnt[0]) + ' ' + dest_dir_bk)
						Header_Add(spc_opt_cnt[0],'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(spc_opt_cnt[1],'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(specfile_bg   ,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						os.system ('cp ' + str(spc_opt_cnt[1]) + ' ' + dest_dir_bk)
						os.system ('cp ' + specfile_bg         + ' ' + dest_dir_bk)
					elif sel_pre_cnt == True and spectra_noise == True:

						Header_Add(specfile_n_bg ,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
					elif sel_pre_cnt == False:
						pass

					if sel_pre_msk == True:
						spc_opt_msk = Spectra_Masking(spc_ipt_msk,msk_typ=sel_msk_type,rshft_corr=z_bg,
										rshft_corr_direct=False,msk_abs_lne=sel_msk_abs_lne,
										msk_spc_rgn=sel_msk_spc_rgn,msk_lmb_min=sel_msk_lmb_min,msk_lmb_max=sel_msk_lmb_max)
						Header_Add(spc_opt_msk  ,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(specfile_bg  ,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						os.system ('cp ' + spc_opt_msk	 + ' ' + dest_dir_bk)
						os.system ('cp ' + specfile_bg   + ' ' + dest_dir_bk)
					elif sel_pre_msk == True and spectra_noise == True:
						Header_Add(specfile_n_bg,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						os.system ('cp ' + specfile_n_bg + ' ' + dest_dir_bk)
					elif sel_pre_msk == False:
						pass

					if sel_pre_shf == True:
						spc_opt_shf1     = Shift_Spectra(spc_ipt_shf    ,z_fg,id_fg,shft_ref=id_fg,id_0=id_bg,redshift_0=z_bg)
						spc_opt_shf1_fit = Shift_Spectra(spc_ipt_shf_fit,z_fg,id_fg,shft_ref=id_fg,id_0=id_bg,redshift_0=z_bg)

						spc_opt_shf2     = Shift_Spectra(spc_ipt_shf    ,z_bg,id_bg)
						spc_opt_shf2_fit = Shift_Spectra(spc_ipt_shf_fit,z_bg,id_bg)

						Header_Add(spc_opt_shf1    ,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(spc_opt_shf1    ,'SHT_FN0',str((spc_ipt_shf_fit.rsplit('/',1)[1]).rsplit('.',1)[0]) ,header_comment='Spec file pre-shift')
						Header_Add(spc_opt_shf1    ,'SHT_FNS',str((spc_opt_shf1_fit.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Spec file post-shift')

						Header_Add(spc_opt_shf1_fit,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(spc_opt_shf1_fit,'SHT_FN0',str((spc_ipt_shf_fit.rsplit('/',1)[1]).rsplit('.',1)[0]) ,header_comment='Spec file pre-shift')
						Header_Add(spc_opt_shf1_fit,'SHT_FNS',str((spc_opt_shf1_fit.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Spec file post-shift')

						Header_Add(spc_opt_shf2    ,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(spc_opt_shf2    ,'SHT_FN0',str((spc_ipt_shf_fit.rsplit('/',1)[1]).rsplit('.',1)[0]) ,header_comment='Spec file pre-shift')
						Header_Add(spc_opt_shf2    ,'SHT_FNS',str((spc_opt_shf2_fit.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Spec file post-shift')

						Header_Add(spc_opt_shf2_fit,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(spc_opt_shf2_fit,'SHT_FN0',str((spc_ipt_shf_fit.rsplit('/',1)[1]).rsplit('.',1)[0]) ,header_comment='Spec file pre-shift')
						Header_Add(spc_opt_shf2_fit,'SHT_FNS',str((spc_opt_shf2_fit.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='Spec file post-shift')

						Header_Add(specfile_bg  ,'MAG_I',magi_bg,header_comment='i-band mag VUDS')

						os.system ('cp ' + spc_opt_shf1      + ' ' + dest_dir_bk)
						os.system ('cp ' + spc_opt_shf1_fit  + ' ' + dest_dir_bk)
						os.system ('cp ' + spc_opt_shf2      + ' ' + dest_dir_bk)
						os.system ('cp ' + spc_opt_shf2_fit  + ' ' + dest_dir_bk)
						os.system ('cp ' + spc_ipt_shf       + ' ' + dest_dir_bk)
						os.system ('cp ' + spc_ipt_shf_fit   + ' ' + dest_dir_bk)

						os.system ('cp ' + specfile_bg       + ' ' + dest_dir_bk)

					if sel_pre_shf == True and spectra_noise == True:
						spc_opt_shf3     = Shift_Spectra(specfile_n_bg,z_fg,id_fg,shft_ref=id_fg,id_0=id_bg,redshift_0=z_bg)
						spc_opt_shf4     = Shift_Spectra(specfile_n_bg,z_bg,id_bg)

						Header_Add(spc_opt_shf3,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(spc_opt_shf4,'MAG_I',magi_bg,header_comment='i-band mag VUDS')
						Header_Add(specfile_n_bg,'MAG_I',magi_bg,header_comment='i-band mag VUDS')

						os.system ('cp ' + spc_opt_shf3      + ' ' + dest_dir_bk)
						os.system ('cp ' + spc_opt_shf4      + ' ' + dest_dir_bk)
						os.system ('cp ' + specfile_n_bg     + ' ' + dest_dir_bk)
 
					elif sel_pre_shf == False:	
						pass
				elif not os.path.isfile(specfile_bg):
					print
					print 'No Spectra File exists: ',specfile_bg

				elif not os.path.isfile(specfile_n_bg):
					print
					print 'No Spectra File (noise) exists: ',specfile_n_bg

			pbar.finish()

		elif len(ident_bg_glx_p) == 1:
			print
			print 'NO BACKGROUND GALAXIES NEARBY'

	pbar1.finish()

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
	msk_spc_rgn       = kwargs.get('msk_spc_rgn'  ,False)
	msk_lmb_min       = kwargs.get('msk_lmb_min'  ,500)
	msk_lmb_max       = kwargs.get('msk_lmb_max'  ,1200)
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
	if msk_spc_rgn == True:
		if rshft_corr_direct == True:
			msk_lmb_min = msk_lmb_min * (rshft_corr)
			msk_lmb_max = msk_lmb_max * (rshft_corr)
		elif rshft_corr_direct == False:
			msk_lmb_min = msk_lmb_min * (1+rshft_corr)
			msk_lmb_max = msk_lmb_max * (1+rshft_corr)
		if (org_spec[0][0]<msk_lmb_min<org_spec[0][-1]) and (org_spec[0][0]<msk_lmb_max<org_spec[0][-1]):
			Spectra_x_y_Updt(mask_ofn,msk_typ,msk_lmb_min,msk_lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
		elif (org_spec[0][0]>msk_lmb_min<org_spec[0][-1]) and (org_spec[0][0]<msk_lmb_max<org_spec[0][-1]):
			msk_lmb_min = org_spec[0][0]
			Spectra_x_y_Updt(mask_ofn,msk_typ,msk_lmb_min,msk_lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
		elif (org_spec[0][0]<msk_lmb_min<org_spec[0][-1]) and (org_spec[0][0]<msk_lmb_max>org_spec[0][-1]):
			msk_lmb_max = org_spec[0][-1]
			Spectra_x_y_Updt(mask_ofn,msk_typ,msk_lmb_min,msk_lmb_max,'lambda',cnt_val=msk_cte_val,spfn_i_2=str(mask_ifn .split('.fits',1)[0]) + '-c-f.fits')
		else:
			pass		
	elif msk_spc_rgn == False:
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
	pyfcnt.input    = Cont_ip_sfn_IRAF
	pyfcnt.output   = Cont_op_sfn_IRAF_1
	pyfcnt.type     = Cont_type_IRAF
	pyfcnt.lines    = Cont_lines_IRAF
	pyfcnt.band     = 1
	pyfcnt.logfile  = Cont_log_IRAF
	pyfcnt.function = Cont_funct_IRAF
	pyfcnt.order    = Cont_order_IRAF
	pyfcnt.replace  = Cont_replace_IRAF
	pyfcnt.low_rej  = Cont_low_rej_IRAF
	pyfcnt.high_rej = Cont_high_rej_IRAF
	pyfcnt.override = Cont_override_IRAF
	pyfcnt.interac  = 'no'
	pyfcnt.mode     = 'h'
	pyfcnt()

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

	pyfcnt.output   = Cont_op_sfn_IRAF_2
	pyfcnt.type     = 'fit'
	pyfcnt()

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
		grid_krnl = apcvl.Box1DKernel(krnl_sz_smt)
	elif pkrnl_typ_smt == 'mexican':
		grid_krnl = apcvl.MexicanHat1DKernel(krnl_sz_smt)

	#http://docs.astropy.org/en/stable/api/apcvl.convolve.html
	spc_pst_smt = apcvl.convolve(spc2bsmth[1], grid_krnl)#,boundary='extend') #fill wrap extend

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



