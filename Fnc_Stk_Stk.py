import bottleneck as bn
#import astropy
#from astropy import stats 
from astropy import table as aptbl
from astropy import stats as apsts
#import itertools
from itertools import chain as itchn
from itertools import product as itpdc

#from termcolor import colored
import scipy.integrate as integrate

from Fnc_Stk_Fts import *
from Fnc_Stk_Mth import *
from Fnc_Stk_Spc import *
from Fnc_Stk_Dir import *
from Fnc_Stk_Utl import *
from Fnc_Stk_Tbl import *
from Fnc_Stk_Plt import *

####Fnc_Stk_Stk####
def monotonic(x):
	dx = np.diff(x)
	return np.all(dx <= 0) or np.all(dx >= 0)

def Stk_Fit_Lines(specfile,*args,**kwargs):
	dest_dir      = kwargs.get('dest_dir',None)
	plt_fit       = kwargs.get('plt_fit',False)
	verbose       = kwargs.get('verbose',False)

	z_glx_Ps      = kwargs.get('z_glx_Ps'    ,0)
	lmb_min       = kwargs.get('lmb_min',1200)
	lmb_max       = kwargs.get('lmb_max',1700)
	
	stk_function  = kwargs.get('stk_function','med')

	fit_fnct      = kwargs.get('fit_fnct','gauss')
	fit_type      = kwargs.get('fit_type','lmfit')
	pre_off_plt   = kwargs.get('pre_off_plt',False)

	lmb_min_lim   = lmb_min
	lmb_max_lim   = lmb_max

	epssave       = kwargs.get('epssave',False)
	showplot      = kwargs.get('showplot',False)

	Cube2bPlot_1D_Err  = kwargs.get('Cube2bPlot_1D_Err', None)

	org_spc_fle   = kwargs.get('org_spc_fle',None)
	ivl_fts_hdr   = kwargs.get('ivl_fts_hdr',False)

	mke_lne_fit   = kwargs.get('mke_lne_fit',True)
	uft_lne_vls   = kwargs.get('uft_lne_vls',False)
	fit_vls_hdr   = kwargs.get('fit_vls_hdr',True)
	ofs_ctr_fit   = kwargs.get('ofs_ctr_fit',True)

	fix_ctr_gau    = kwargs.get('fix_ctr_gau',False)
	fix_pre_gau    = kwargs.get('fix_pre_gau',False)
	fix_pst_gau    = kwargs.get('fix_pst_gau',False)
	pre_shf_lim    = kwargs.get('pre_shf_lim',5 )
	pst_shf_lim    = kwargs.get('pst_shf_lim',15)
	pre_shf_ctr    = kwargs.get('pre_shf_ctr',2.5 )
	pst_shf_ctr    = kwargs.get('pst_shf_ctr',20)

	ivl_fts_hdr    = kwargs.get('ivl_fts_hdr',False)
	mke_lne_fit    = kwargs.get('mke_lne_fit',True)

	fix_ctr_gau_1  = kwargs.get('fix_ctr_gau_1',False)
	fix_ctr_gau_2  = kwargs.get('fix_ctr_gau_2',False)
	fix_mdl_gau    = kwargs.get('fix_mdl_gau',False)
	mdl_shf_ctr    = kwargs.get('mdl_shf_ctr',1 )
	mdl_shf_lim    = kwargs.get('mdl_shf_lim',5)

	int_vlf_hdr    = kwargs.get('int_vlf_hdr',True)
	fit_vls_hdr    = kwargs.get('fit_vls_hdr',True)
	uft_lne_vls    = kwargs.get('uft_lne_vls',False)


	import Lines_Dictionary
	LINES = Lines_Dictionary.LINES_PLT_BG
	print
	print colored('Fitting ' + str(len(LINES[0])) + ' lines in the range ' +str(lmb_min_lim) +'-'+str(lmb_max_lim),'yellow') 
	print "\n".join([lineinrange[0] + '-' +str(lineinrange[1]) for lineinrange in zip(LINES[4],LINES[0])])
	print

	if 'Fg' in specfile:
		glx_type = 'Fg'
		clr_plt  = 'red'
		clr_fit  = 'blue'
	elif 'Bg' in specfile:
		glx_type = 'Bg'
		clr_plt  = 'blue'
		clr_fit  = 'red'

	sep_label     = str((specfile.split('-stk',1)[0]).rsplit('_as-',2)[1])
	plt_sufix_fnm = ((specfile.rsplit('/')[-1])).split('.fits')[0]

	try:
		z_glx_zro = Header_Get(specfile,'Z_0')
	except KeyError:
		z_glx_zro = 0

	try:
		z_glx_ref = Header_Get(specfile,'Z_REF')
	except KeyError:
		z_glx_ref = 0

	z_glx_rsf=0

	lambda_sp,inten_sp,crval_sp,cdel1_sp = Spectra_x_y(specfile)[0], Spectra_x_y(specfile)[1],Spectra_x_y(specfile)[2],Spectra_x_y(specfile)[3]
	label_sp = glx_type +' (sep: ' + sep_label + ' arcsec) '  + ' ('+ str(Header_Get(specfile,'STK_OPR'))  + ' N = '+str(str(Header_Get(specfile,'stk_num')))  + ') STD: ' #+ str(np.round(std_glx_med,4))
	###
	MSK_NTMS=2.5
	for lines in range(len(LINES[0])):
		if lmb_min < LINES[0][lines]*(1+z_glx_Ps) < lmb_max :
			print
			print colored('Fitting line:','yellow')
			print LINES[4][lines]+'-'+str(LINES[0][lines])
			print
			########################################################PLOT PER LINE########################################################
			fxsize=11
			fysize=8
			f = plt.figure(num=None, figsize=(fxsize, fysize), dpi=180, facecolor='w',
				edgecolor='k')
			plt.subplots_adjust(
				left 	= (26/25.4)/fxsize, 
				bottom 	= (16/25.4)/fysize, 
				right 	= 1 - (4/25.4)/fxsize, 
				top 	= 1 - (4/25.4)/fysize)
			plt.subplots_adjust(hspace=0)

			#f.suptitle('An overall title', size=20)
			gs0 = gridspec.GridSpec(1, 1)

			#############################################################STACK###########################################################
			gs11 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[0])
				
			ax110 = plt.Subplot(f, gs11[0:3,0])
			f.add_subplot(ax110)

			ax110.set_rasterization_zorder(1)
			plt.autoscale(enable=True, axis='y', tight=False)
			ax110.xaxis.set_tick_params(labelsize=16)
			ax110.yaxis.set_tick_params(labelsize=16)
			#ax110.set_title(PLOT_TITLE)
			xticklabels = ax110.get_xticklabels()
			plt.setp(xticklabels, visible=True)
			yticklabels = ax110.get_yticklabels()
			plt.setp(yticklabels, visible=True)
			ax110.yaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)
			ax110.xaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)

			minorLocator_x   = plt.MultipleLocator(1)
			majorLocator_x   = plt.MultipleLocator(10)
			ax110.xaxis.set_minor_locator(minorLocator_x)
			ax110.xaxis.set_major_locator(majorLocator_x)
			plt.tick_params(which='both', width=0.7)
			plt.tick_params(which='major', length=5)
			plt.tick_params(which='minor', length=2)
			ax110.minorticks_on()

			plt.xlabel('$\lambda$',fontsize=16)
			plt.ylabel('F$_\lambda$ (ergs/sec/cm$^{2}/\mathrm{\AA}$)',fontsize=16)

			if 'Bg' in specfile:
				colors = "bgrcmykw"
			elif 'Fg' in specfile:
				colors = "rgbcmykw"
			FILES = [specfile]
			for index,specfile_glx in enumerate(FILES):
				glx = Spectra_x_y(specfile_glx)
				lambda_glx,inten_glx,crval_glx,cdel1_glx,cd1_glx = glx[0], glx[1], glx[2], glx[3], glx[4]
				stk_glx_nmb = Header_Get(specfile_glx,'STK_NUM')

				print
				print colored('Fitting stacked spectrum file: '+str(specfile_glx),'yellow')
				print colored('Fitting stacked spectrum original file: '+str(org_spc_fle),'yellow')
				print

				#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
				if ivl_fts_hdr == True:
					try:
						L1_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF_0')        #LINES-1 Wdt-Fit  1GF-IntVal      WIDTH-FIT
						L2_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP_0')        #LINES-2 Wdt-Plt  1GF-IntVal      WIDTH-PLT
						L7_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF_0')        #LINES-7 Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
						L8_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
						#L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')       #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
						print
						print colored('Initial fit variables from fits header!','yellow')
						print colored('Headers:','yellow')
						print
						print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','yellow')
						print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','yellow')
						print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','yellow')
						print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','yellow')
						#print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
					except ValueError:
						print
						print colored('Headers containing initial fit variables NOT found!','yellow')
						print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
						print colored('Headers:','yellow')
						print
						print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','yellow')
						print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','yellow')
						print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','yellow')
						print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','yellow')
						print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
						print
						quit()
					try:
						L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
						print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
						print
					except KeyError:
						print
						print colored('Header NOT found!','yellow')
						print colored('Adding Header with default valuee 0.001:','yellow')
						print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
						Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF_0',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
						L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
					################FOR CASES WHERE KLINE WIIDHT ==0################
					if L1_0  == 0:
						L1_0  = 1#LINES[1][lines]
					else:
						pass
					################FOR CASES WHERE KLINE WIIDHT ==0################
				elif ivl_fts_hdr == False:
					L1_0  = LINES[1][lines]
					L2_0  = LINES[2][lines]
					L7_0  = LINES[7][lines]
					L8_0  = LINES[8][lines]
					L10_0 = LINES[10][lines]
					print
					print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
					print					
				
				print str(LINES[5][lines])+'_WF_0' + ': ' + str(L1_0)
				print str(LINES[5][lines])+'_WP_0' + ': ' + str(L2_0)
				print str(LINES[5][lines])+'_CF_0' + ': ' + str(L7_0)
				print str(LINES[5][lines])+'_CO_0' + ': ' + str(L8_0)
				print str(LINES[5][lines])+'_AF_0' + ': ' + str(L10_0)

				#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################

				########################################################LINE-FIT#######################################################
				if   'Dbl' in   LINES[3][lines] and fit_fnct=='gauss' and fit_type == 'lmfit' and mke_lne_fit == True and uft_lne_vls == False:
					fit_typ = 'G'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					print 
					print colored('Double Line Fit (Ind)','yellow')
					print colored(LINES[3][lines-2]  + '-' + str(LINES[0][lines-2])   + '-' + str(LINES[1][lines-2]),'cyan')
					print LINES[3][lines]+ '-' + str(LINES[0][lines]) + '-'  + str(LINES[1][lines])
					print colored(LINES[3][lines-1]+ '-' + str(LINES[0][lines-1]) + '-' + str(LINES[1][lines-1]),'magenta')

					from lmfit import Model
	
					########################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#########################
					if ivl_fts_hdr == True:
						try:
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial fit variables from fits header!','yellow')
							print
							print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WF02' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WP02' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CF02' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CO02' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AF02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print
						except ValueError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print '*****'
							print
							quit()
						try:
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')							
							print colored(str(LINES[5][lines-2])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')							
							print colored(str(LINES[5][lines-1])+'_AF02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF01',0.001,header_comment = str(LINES[3][lines-2]) + str(LINES[0][lines-2]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF02',0.001,header_comment = str(LINES[3][lines-1]) + str(LINES[0][lines-1]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines-2])+'_AF01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines-1])+'_AF02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('*****Success!******','yellow')
							print
						################FOR CASES WHERE KLINE WIIDHT ==0################
						if L1_1  == 0:
							L1_1  = 1#LINES[1][lines]
						else:
							pass
						if L1_2  == 0:
							L1_2  = 1#LINES[1][lines]
						else:
							pass
						################FOR CASES WHERE KLINE WIIDHT ==0################
					elif ivl_fts_hdr == False:
						L1_1  = LINES[1][lines-2]
						L2_1  = LINES[2][lines-2]
						L7_1  = LINES[7][lines-2]
						L8_1  = LINES[8][lines-2]
						L10_1 = LINES[10][lines-2]

						L1_2  = LINES[1][lines-1]
						L2_2  = LINES[2][lines-1]
						L7_2  = LINES[7][lines-1]
						L8_2  = LINES[8][lines-1]
						L10_2 = LINES[10][lines-1]
						print
						print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
						print					
					print colored('Initial Values: ','cyan')
					print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')					
					print colored(str(LINES[5][lines])+'_WF01' + ': ' + str(L1_1),'cyan')
					print colored(str(LINES[5][lines])+'_WP01' + ': ' + str(L2_1),'cyan')
					print colored(str(LINES[5][lines])+'_CF01' + ': ' + str(L7_1),'cyan')
					print colored(str(LINES[5][lines])+'_CO01' + ': ' + str(L8_1),'cyan')
					print colored(str(LINES[5][lines])+'_AF01' + ': ' + str(L10_1),'cyan')
					print
					print colored('Initial Values: ','magenta')
					print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')					
					print colored(str(LINES[5][lines])+'_WF02' + ': ' + str(L1_2),'magenta')
					print colored(str(LINES[5][lines])+'_WP02' + ': ' + str(L2_2),'magenta')
					print colored(str(LINES[5][lines])+'_CF02' + ': ' + str(L7_2),'magenta')
					print colored(str(LINES[5][lines])+'_CO02' + ': ' + str(L8_2),'magenta')
					print colored(str(LINES[5][lines])+'_AF02' + ': ' + str(L10_2),'magenta')
					print
					########################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#########################
					##################################################CENTRAL GAUSSIAN-1###################################################
					lmb_min_lim_line_ft = (LINES[0][lines-2]+L8_1) - MSK_NTMS*L1_1 
					lmb_max_lim_line_ft = (LINES[0][lines-2]+L8_1) + MSK_NTMS*L1_1
					lmb_min_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_1#- 20#L2_1 - 10 
					lmb_max_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_1#+ 20#L2_1 + 10

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function

					idx_ctr_ft_reg = np.abs(lambda_glx[mask_ft] - (LINES[0][lines-2]+L8_1)).argmin()

					X0_f2DG    = (LINES[0][lines-2]+L8_1)
					SIGMA_f2DG = L1_1
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))
					A_f2DG     = -(1-inten_glx[mask_ft][idx_ctr_ft_reg])

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
				
					#################################################CENTRAL GAUSSIAN-1-C##################################################
					if fix_ctr_gau_1 == False:
						print
						print colored('Fitting 1st line','cyan')
						print colored('1-0-Fitting Central line','cyan')
						print
						#################################################CENTRAL GAUSSIAN-1-0##################################################
						try:
							gmodel_0           = Model(func_1D_Gaussian)
							gmodel_0.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG - (X0_f2DG*L7_1), max=X0_f2DG + (X0_f2DG*L7_1))
							gmodel_0.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*L10_1), max=A_f2DG  + (A_f2DG*L10_1))
							gmodel_0.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_0             = gmodel_0.make_params()							
							result_0_1         = gmodel_0.fit(inten_glx[mask_ft],pars_0,
													X=lambda_glx[mask_ft],
													X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_0_1         = result_0_1.params['X_0'].value
							AMPL_G_0_1         = result_0_1.params['A'].value
							SGMA_G_0_1         = abs(result_0_1.params['SIGMA'].value)
							FWHM_G_0_1         = lw_sgma2fwhm(SGMA_G_0_1)
							W_0_1              = integrate.quad(lambda x: AMPL_G_0_1*np.exp(-((x)**2)/(2*SGMA_G_0_1**2)), -np.inf, np.inf)
							EW_0_1             = np.round(abs(np.asarray(W_0_1[0])),10)
							EWE_0_1            = np.round(abs(np.asarray(W_0_1[1])),10)
							data_fitted_0      = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_0_1,AMPL_G_0_1,SGMA_G_0_1)

							CTRE_G_0_1_E       = result_0_1.params['X_0'].stderr
							AMPL_G_0_1_E       = result_0_1.params['A'].stderr
							SGMA_G_0_1_E       = result_0_1.params['SIGMA'].stderr

							CTRE_G_0_1_cor     = result_0_1.params['X_0'].correl
							AMPL_G_0_1_cor     = result_0_1.params['A'].correl
							SGMA_G_0_1_cor     = result_0_1.params['SIGMA'].correl

							chisqr_0_1         = result_0_1.chisqr
							redchi_0_1         = result_0_1.redchi
						except (RuntimeError,ValueError,TypeError):
							popt_0, pcov_0   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_0           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_0         = 999999.99999
							AMPL_G_0         = 999999.99999
							SGMA_G_0         = 999999.99999
							FWHM_G_0         = 999999.99999
							EW_0             = 999999.99999
							EWE_0            = 999999.99999

							CTRE_G_0_E      = 999999.99999
							AMPL_G_0_E      = 999999.99999
							SGMA_G_0_E      = 999999.99999

							CTRE_G_0_cor    = 999999.99999
							AMPL_G_0_cor    = 999999.99999
							SGMA_G_0_cor    = 999999.99999

							chisqr_0        = 999999.99999
							redchi_0        = 999999.99999
						if fit_vls_hdr == True and fix_ctr_gau_1 ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF01',float(CTRE_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF01',float(AMPL_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FF01',float(FWHM_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF01',float(EW_0_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EF01',float(EWE_0_1)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('The fit (CTR-1-0) values will be added to the fits headers!','cyan')
							print
						else:
							print
							print colored('The fit (CTR-1-0) values will NOT be added to the fits headers!','cyan')
							print
						#################################################CENTRAL GAUSSIAN-1-0##################################################
						#################################################CENTRAL GAUSSIAN-1-C##################################################
						try:
							gmodel_O           = Model(func_1D_Gaussian_O)
							gmodel_O.set_param_hint('X_0'   , value=X0_f2DG   , min=X0_f2DG-(X0_f2DG*L7_1), max=X0_f2DG+(X0_f2DG*L7_1))
							gmodel_O.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*L10_1) , max=A_f2DG  + (A_f2DG*L10_1))
							gmodel_O.set_param_hint('SIGMA' , value=SIGMA_f2DG)
							gmodel_O.set_param_hint('OFFSET', value=max(inten_glx[mask_ft])-1)
							pars_O             = gmodel_O.make_params()

							result_O_1         = gmodel_O.fit(inten_glx[mask_ft],pars_O,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,OFFSET=max(inten_glx[mask_ft])-1,
													nan_policy = 'omit')
							CTRE_G_O_1         = result_O_1.params['X_0'].value
							AMPL_G_O_1         = result_O_1.params['A'].value
							SGMA_G_O_1         = abs(result_O_1.params['SIGMA'].value)
							OFST_G_O_1         = abs(result_O_1.params['OFFSET'].value)
							FWHM_G_O_1         = lw_sgma2fwhm(SGMA_G_O_1)
							W_O_1              = integrate.quad(lambda x: AMPL_G_O_1*np.exp(-((x)**2)/(2*SGMA_G_O_1**2)), -np.inf, np.inf)
							EW_O_1             = np.round(abs(np.asarray(W_O_1[0])),10)
							EWE_O_1            = np.round(abs(np.asarray(W_O_1[1])),10)
							data_fitted_O_1    = func_1D_Gaussian_O((lambda_glx[mask_ft]),CTRE_G_O_1,AMPL_G_O_1,SGMA_G_O_1,OFST_G_O_1)

							CTRE_G_O_E         = result_O_1.params['X_0'].stderr
							AMPL_G_O_E         = result_O_1.params['A'].stderr
							SGMA_G_O_E         = result_O_1.params['SIGMA'].stderr

							CTRE_G_O_cor       = result_O_1.params['X_0'].correl
							AMPL_G_O_cor       = result_O_1.params['A'].correl
							SGMA_G_O_cor       = result_O_1.params['SIGMA'].correl

							chisqr_O_1         = result_O_1.chisqr
							redchi_O_1         = result_O_1.redchi
							
							#####################################################################################################################
							if ofs_ctr_fit == True:
								print
								print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
								print
								lmb_min_lim_line_ft = (CTRE_G_O_1-L8_1) - MSK_NTMS*L1_1 
								lmb_max_lim_line_ft = (CTRE_G_O_1+L8_1) + MSK_NTMS*L1_1
								mask_ft     = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
							else:
								pass
							#####################################################################################################################
							initial_guess_C    = (X0_f2DG,A_f2DG,SIGMA_f2DG)

							gmodel_C           = Model(func_1D_Gaussian)
							gmodel_C.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG-(X0_f2DG*L7_1), max=X0_f2DG+(X0_f2DG*L7_1))
							gmodel_C.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*L10_1) , max=A_f2DG  + (A_f2DG*L10_1))
							gmodel_C.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_C             = gmodel_C.make_params()
							result_C_1         = gmodel_C.fit(inten_glx[mask_ft],pars_C,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_C_1         = result_C_1.params['X_0'].value
							AMPL_G_C_1         = result_C_1.params['A'].value
							SGMA_G_C_1         = abs(result_C_1.params['SIGMA'].value)
							FWHM_G_C_1         = lw_sgma2fwhm(SGMA_G_C_1)

							W_C_1              = integrate.quad(lambda x: AMPL_G_C_1*np.exp(-((x)**2)/(2*SGMA_G_C_1**2)), -np.inf, np.inf)
							EW_C_1             = np.round(abs(np.asarray(W_C_1[0])),10)
							EWE_C_1            = np.round(abs(np.asarray(W_C_1[1])),10)
							data_fitted_C_1    = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_C_1,AMPL_G_C_1,SGMA_G_C_1)

							CTRE_G_C_1_E         = result_C_1.params['X_0'].stderr
							AMPL_G_C_1_E         = result_C_1.params['A'].stderr
							SGMA_G_C_1_E         = result_C_1.params['SIGMA'].stderr

							CTRE_G_C_1_cor       = result_C_1.params['X_0'].correl
							AMPL_G_C_1_cor       = result_C_1.params['A'].correl
							SGMA_G_C_1_cor       = result_C_1.params['SIGMA'].correl

							AMPL_SNR_1           = AMPL_G_C_1
							CTRE_SNR_1           = CTRE_G_C_1
							SGMA_SNR_1           = abs(SGMA_G_C_1)

							if CTRE_G_C_1_E == None:
								CTRE_G_C_1_E = 999999.99999
							else:
								pass
							if AMPL_G_C_1_E == None:
								AMPL_G_C_1_E = 999999.99999
							else:
								pass
							if SGMA_G_C_1_E == None:
								SGMA_G_C_1_E = 999999.99999
							else:
								pass
							if CTRE_G_C_1_cor == None:
								CTRE_G_C_1_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_1_cor == None:
								AMPL_G_C_1_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_1_cor == None:
								SGMA_G_C_1_cor = 999999.99999
							else:
								pass
							chisqr_C_1      = result_C_1.chisqr
							redchi_C_1      = result_C_1.redchi
							##inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O #OFFSET +
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C_1, pcov_C_1 = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_1           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_1         = 999999.99999
							AMPL_G_C_1         = 999999.99999
							SGMA_G_C_1         = 999999.99999
							FWHM_G_C_1         = 999999.99999
							EW_C_1             = 999999.99999
							EWE_C_1            = 999999.99999

							CTRE_G_C_1_E       = 999999.99999
							AMPL_G_C_1_E       = 999999.99999
							SGMA_G_C_1_E       = 999999.99999
							CTRE_G_C_1_cor     = 999999.99999
							AMPL_G_C_1_cor     = 999999.99999
							SGMA_G_C_1_cor     = 999999.99999
							chisqr_C_1         = 999999.99999
							redchi_C_1         = 999999.99999

							popt_O_1 ,pcov_O_1 = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
							perr_O_1           = [999999.99999,999999.99999,999999.99999,999999.99999]
							CTRE_G_O_1         = 999999.99999
							AMPL_G_O_1         = 999999.99999
							SGMA_G_O_1         = 999999.99999
							OFST_G_O_1         = 999999.99999
							FWHM_G_O_1         = 999999.99999
							EW_O_1             = 999999.99999
							EWE_O_1            = 999999.99999

							CTRE_G_O_1_E      = 999999.99999
							AMPL_G_O_1_E      = 999999.99999
							SGMA_G_O_1_E      = 999999.99999
							CTRE_G_O_1_cor    = 999999.99999
							AMPL_G_O_1_cor    = 999999.99999
							SGMA_G_O_1_cor    = 999999.99999
							OFST_G_O_1_cor    = 999999.99999
							chisqr_O_1        = 999999.99999
							redchi_O_1        = 999999.99999

							AMPL_SNR_1        = 999999.99999
							CTRE_SNR_1        = 999999.99999
							SGMA_SNR_1        = 999999.99999
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines-2])+'-'+str(LINES[3][lines-2]),'cyan')
						print
						print colored(str(LINES[5][lines-2])+'_CGLC: ' + str(CTRE_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_AGLC: ' + str(AMPL_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_SGLC: ' + str(SGMA_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_FGLC: ' + str(FWHM_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_WGLC: ' + str(EW_C_1)      ,'cyan')
						print colored(str(LINES[5][lines-2])+'_EGLC: ' + str(EWE_C_1)     ,'cyan')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O
						if fit_vls_hdr == True and fix_ctr_gau_1==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                         ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                         ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLO1',float(CTRE_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALO1',float(AMPL_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLO1',float(FWHM_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLO1',float(EW_O_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELO1',float(EWE_O_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OFO1',float(OFST_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLC1',float(CTRE_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALC1',float(AMPL_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLC1',float(SGMA_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLC1',float(FWHM_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLC1',float(EW_C_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELC1',float(EWE_C_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CEC1',float(CTRE_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AEC1',float(AMPL_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SEC1',float(SGMA_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CEC1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AEC1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SEC1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))								
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHL1',float(chisqr_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRL1',float(redchi_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines-2]) + ' Chi2 Reduced 1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('The fit (CTR-1-C & CTR-1-O) values will be added to the fits headers!','cyan')
							print
						else:
							print
							print colored('The fit (CTR-1-C & CTR-1-O) values will NOT be added to the fits headers!','cyan')
							print

						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values G-1 for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF01',float(L1_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-1-1 Wdt-Fit  1GF-IntVal')         #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WP01',float(L2_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-2-1 Wdt-Plt  1GF-IntVal')         #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF01',float(L7_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-7-1 Ctr Fit Bnds  1GF-IntVal')    #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CO01',float(L8_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-8-1 Ctr Fit Ofst  1GF-IntVal')    #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF01',float(L10_1),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10-1 Amp Fit Bnds  1GF-IntVal')   #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial Values: ','cyan')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')					
							print colored(str(LINES[5][lines])+'_WF01' + ': ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WP01' + ': ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CF01' + ': ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CO01' + ': ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AF01' + ': ' + str(L10_1),'cyan')
							print
						else:
							print
							print colored('Initial Guess Values G-1 for line Fitting will NOT be recorded!','yellow')
							print
							pass
					#################################################CENTRAL GAUSSIAN-1-C##################################################
					elif fix_ctr_gau_1 == True:
						print
						print colored('1-CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						try:
							CTRE_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF01')
							AMPL_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF01')
							FWHM_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF01')
							EW_0_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF01')
							EWE_0_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF01')

							CTRE_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLO1')
							AMPL_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALO1')
							FWHM_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLO1')
							EW_O_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLO1')
							EWE_O_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELO1')
							OFST_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFO1')

							CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC1')
							AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC1')
							SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC1')
							FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC1')
							EW_C_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLC1')
							EWE_C_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELC1')
							CTRE_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CEC1')
							AMPL_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AEC1')
							SGMA_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SEC1')

							chisqr_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHL1')
							redchi_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRL1')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines-2])+'_CF01','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau_1=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CGF0')
							AMPL_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_AGF0')
							FWHM_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_FGF0')
							EW_0_1        = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_WGF0')
							EWE_0_1       = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_EGF0')

							CTRE_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CGLO')
							AMPL_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_AGLO')
							FWHM_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_FGLO')
							EW_O_1        = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_WGLO')
							EWE_O_1       = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_EGLO')
							OFST_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_OFSO')

							CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CGLC')
							AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_AGLC')
							SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_SGLC')
							FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_FGLC')
							EW_C_1        = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_WGLC')
							EWE_C_1       = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_EGLC')
							CTRE_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CLEC')
							AMPL_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_ALEC')
							SGMA_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_SLEC')

							chisqr_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CHGL')
							redchi_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CRGL')
					print
					print colored(str(LINES[0][lines-2])+'-'+str(LINES[3][lines-2])+'-CTR','cyan')
					print
					print colored(str(LINES[5][lines-2])+'_CLC1: ' + str(CTRE_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines-2])+'_ALC1: ' + str(AMPL_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines-2])+'_SLC1: ' + str(SGMA_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines-2])+'_FLC1: ' + str(FWHM_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines-2])+'_WLC1: ' + str(EW_C_1)      ,'cyan')
					print colored(str(LINES[5][lines-2])+'_ELC1: ' + str(EWE_C_1),'cyan')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','cyan')
					print colored(str(CTRE_G_C_1)+', '+str(AMPL_G_C_1)+', '+str(SGMA_G_C_1),'cyan')
					print
					##################################################CENTRAL GAUSSIAN-1###################################################
					##################################################CENTRAL GAUSSIAN-2###################################################
					lmb_min_lim_line_ft = (LINES[0][lines-1]+L8_2) - MSK_NTMS*LINES[1][lines-1] 
					lmb_max_lim_line_ft = (LINES[0][lines-1]+L8_2) + MSK_NTMS*LINES[1][lines-1]
					lmb_min_lim_line    = (LINES[0][lines-1]+L8_2)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_2#- 20#L2_2 - 10 
					lmb_max_lim_line    = (LINES[0][lines-1]+L8_2)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_2#+ 20#L2_2 + 10

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function
					idx_ctr_ft_reg = np.abs(lambda_glx[mask_ft] - (LINES[0][lines-1]+L8_2)).argmin()

					X0_f2DG    = (LINES[0][lines-1]+L8_2)
					SIGMA_f2DG = LINES[1][lines-1]
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))
					A_f2DG     = -(1-inten_glx[mask_ft][idx_ctr_ft_reg])

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
					##################################################CENTRAL GAUSSIAN-2###################################################
					if fix_ctr_gau_2 == False:
						print
						print colored('Fitting 2nd line','magenta')
						print colored('2-0-Fitting Central line','magenta')
						print
						#################################################CENTRAL GAUSSIAN-2-0##################################################
						try:
							gmodel_0           = Model(func_1D_Gaussian)
							gmodel_0.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG - (X0_f2DG*L7_2), max=X0_f2DG+(X0_f2DG*L7_2))
							gmodel_0.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*L10_2), max=A_f2DG  + (A_f2DG*L10_2))
							gmodel_0.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_0             = gmodel_0.make_params()							
							result_0_2         = gmodel_0.fit(inten_glx[mask_ft],pars_0,
													X=lambda_glx[mask_ft],
													X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_0_2         = result_0_2.params['X_0'].value
							AMPL_G_0_2         = result_0_2.params['A'].value
							SGMA_G_0_2         = abs(result_0_2.params['SIGMA'].value)
							FWHM_G_0           = lw_sgma2fwhm(SGMA_G_0_2)
							W_0_2              = integrate.quad(lambda x: AMPL_G_0_2*np.exp(-((x)**2)/(2*SGMA_G_0_2**2)), -np.inf, np.inf)
							EW_0_2             = np.round(abs(np.asarray(W_0_2[0])),10)
							EWE_0_2            = np.round(abs(np.asarray(W_0_2[1])),10)
							data_fitted_0_2    = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_0_2,AMPL_G_0_2,SGMA_G_0_2)

							CTRE_G_0_2_E       = result_0_2.params['X_0'].stderr
							AMPL_G_0_2_E       = result_0_2.params['A'].stderr
							SGMA_G_0_2_E       = result_0_2.params['SIGMA'].stderr

							CTRE_G_0_2_cor     = result_0_2.params['X_0'].correl
							AMPL_G_0_2_cor     = result_0_2.params['A'].correl
							SGMA_G_0_2_cor     = result_0_2.params['SIGMA'].correl

							chisqr_0           = result_0_2.chisqr
							redchi_0           = result_0_2.redchi
						except (RuntimeError,ValueError,TypeError):
							popt_0, pcov_0   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_0           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_0_2         = 999999.99999
							AMPL_G_0_2         = 999999.99999
							SGMA_G_0_2         = 999999.99999
							FWHM_G_0         = 999999.99999
							EW_0_2             = 999999.99999
							EWE_0_2            = 999999.99999

							CTRE_G_0_2_E      = 999999.99999
							AMPL_G_0_2_E      = 999999.99999
							SGMA_G_0_2_E      = 999999.99999

							CTRE_G_0_2_cor    = 999999.99999
							AMPL_G_0_2_cor    = 999999.99999
							SGMA_G_0_2_cor    = 999999.99999

							chisqr_0        = 999999.99999
							redchi_0        = 999999.99999
						if fit_vls_hdr == True and fix_ctr_gau_2 ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF02',float(CTRE_G_0_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF02',float(AMPL_G_0_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FF02',float(FWHM_G_0)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF02',float(EW_0_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EF02',float(EWE_0_2)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('The fit (CTR-0) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-0) values will NOT be added to the fits headers!','magenta')
							print
						#################################################CENTRAL GAUSSIAN-2-0##################################################
						#################################################CENTRAL GAUSSIAN-2-C##################################################
						try:
							gmodel_O           = Model(func_1D_Gaussian_O)
							gmodel_O.set_param_hint('X_0'   , value=X0_f2DG  , min=X0_f2DG - (X0_f2DG*L7_2), max=X0_f2DG + (X0_f2DG*L7_2))
							gmodel_O.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*L10_2), max=A_f2DG  + (A_f2DG*L10_2))
							gmodel_O.set_param_hint('SIGMA' , value=SIGMA_f2DG)
							gmodel_O.set_param_hint('OFFSET', value=max(inten_glx[mask_ft])-1)
							pars_O             = gmodel_O.make_params()

							result_O_2         = gmodel_O.fit(inten_glx[mask_ft],pars_O,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,OFFSET=max(inten_glx[mask_ft])-1,
													nan_policy = 'omit')
							CTRE_G_O_2         = result_O_2.params['X_0'].value
							AMPL_G_O_2         = result_O_2.params['A'].value
							SGMA_G_O_2         = abs(result_O_2.params['SIGMA'].value)
							OFST_G_O_2         = abs(result_O_2.params['OFFSET'].value)
							FWHM_G_O_2         = lw_sgma2fwhm(SGMA_G_O_2)
							W_O_2              = integrate.quad(lambda x: AMPL_G_O_2*np.exp(-((x)**2)/(2*SGMA_G_O_2**2)), -np.inf, np.inf)
							EW_O_2             = np.round(abs(np.asarray(W_O_2[0])),10)
							EWE_O_2            = np.round(abs(np.asarray(W_O_2[1])),10)
							data_fitted_O_2    = func_1D_Gaussian_O((lambda_glx[mask_ft]),CTRE_G_O_2,AMPL_G_O_2,SGMA_G_O_2,OFST_G_O_2)

							CTRE_G_O_2_E       = result_O_2.params['X_0'].stderr
							AMPL_G_O_2_E       = result_O_2.params['A'].stderr
							SGMA_G_O_E         = result_O_2.params['SIGMA'].stderr

							CTRE_G_O_2_cor     = result_O_2.params['X_0'].correl
							AMPL_G_O_2_cor     = result_O_2.params['A'].correl
							SGMA_G_O_cor       = result_O_2.params['SIGMA'].correl

							chisqr_O_2         = result_O_2.chisqr
							redchi_O_2         = result_O_2.redchi
							
							#####################################################################################################################
							if ofs_ctr_fit == True:
								print
								print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
								print
								lmb_min_lim_line_ft = (CTRE_G_O_2-L8_2) - MSK_NTMS*LINES[1][lines-1] 
								lmb_max_lim_line_ft = (CTRE_G_O_2+L8_2) + MSK_NTMS*LINES[1][lines-1]
								#lmb_min_lim_line    = (CTRE_G_O-L8_2)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_2#-20#L2_2 - 10
								#lmb_max_lim_line    = (CTRE_G_O+L8_2)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_2#+20#L2_2 + 10
								#mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
								mask_ft     = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
							else:
								pass
							#####################################################################################################################
							initial_guess_C    = (X0_f2DG,A_f2DG,SIGMA_f2DG)

							gmodel_C           = Model(func_1D_Gaussian)
							gmodel_C.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG - (X0_f2DG*(X0_f2DG*L7_2)), max=X0_f2DG + (X0_f2DG*(X0_f2DG*L7_2)))
							gmodel_C.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG *(A_f2DG*L10_2)), max=A_f2DG  + (A_f2DG*(A_f2DG*L10_2)))
							gmodel_C.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_C             = gmodel_C.make_params()
							result_C_2         = gmodel_C.fit(inten_glx[mask_ft],pars_C,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_C_2         = result_C_2.params['X_0'].value
							AMPL_G_C_2         = result_C_2.params['A'].value
							SGMA_G_C_2         = abs(result_C_2.params['SIGMA'].value)
							FWHM_G_C_2         = lw_sgma2fwhm(SGMA_G_C_2)

							W_C_2              = integrate.quad(lambda x: AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2)), -np.inf, np.inf)
							EW_C_2             = np.round(abs(np.asarray(W_C_2[0])),10)
							EWE_C_2            = np.round(abs(np.asarray(W_C_2[1])),10)
							data_fitted_C      = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_C_2,AMPL_G_C_2,SGMA_G_C_2)

							CTRE_G_C_2_E       = result_C_2.params['X_0'].stderr
							AMPL_G_C_2_E       = result_C_2.params['A'].stderr
							SGMA_G_C_2_E       = result_C_2.params['SIGMA'].stderr

							CTRE_G_C_2_cor     = result_C_2.params['X_0'].correl
							AMPL_G_C_2_cor     = result_C_2.params['A'].correl
							SGMA_G_C_2_cor     = result_C_2.params['SIGMA'].correl

							AMPL_SNR_2           = AMPL_G_C_2
							CTRE_SNR_2           = CTRE_G_C_2
							SGMA_SNR_2           = abs(SGMA_G_C_2)

							if CTRE_G_C_2_E == None:
								CTRE_G_C_2_E = 999999.99999
							else:
								pass
							if AMPL_G_C_2_E == None:
								AMPL_G_C_2_E = 999999.99999
							else:
								pass
							if SGMA_G_C_2_E == None:
								SGMA_G_C_2_E = 999999.99999
							else:
								pass
							if CTRE_G_C_2_cor == None:
								CTRE_G_C_2_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_2_cor == None:
								AMPL_G_C_2_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_2_cor == None:
								SGMA_G_C_2_cor = 999999.99999
							else:
								pass
							chisqr_C_2        = result_C_2.chisqr
							redchi_C_2        = result_C_2.redchi
							##inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O #OFFSET +
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C_2, pcov_C_2 = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_2           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_2         = 999999.99999
							AMPL_G_C_2         = 999999.99999
							SGMA_G_C_2         = 999999.99999
							FWHM_G_C_2         = 999999.99999
							EW_C_2             = 999999.99999
							EWE_C_2            = 999999.99999

							CTRE_G_C_2_E       = 999999.99999
							AMPL_G_C_2_E       = 999999.99999
							SGMA_G_C_2_E       = 999999.99999
							CTRE_G_C_2_cor     = 999999.99999
							AMPL_G_C_2_cor     = 999999.99999
							SGMA_G_C_2_cor     = 999999.99999
							chisqr_C_2         = 999999.99999
							redchi_C_2         = 999999.99999

							popt_O_2 ,pcov_O_2 = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
							perr_O_2           = [999999.99999,999999.99999,999999.99999,999999.99999]
							CTRE_G_O_2         = 999999.99999
							AMPL_G_O_2         = 999999.99999
							SGMA_G_O_2         = 999999.99999
							OFST_G_O_2         = 999999.99999
							FWHM_G_O_2         = 999999.99999
							EW_O_2             = 999999.99999
							EWE_O_2            = 999999.99999

							CTRE_G_O_2_E       = 999999.99999
							AMPL_G_O_2_E       = 999999.99999
							SGMA_G_O_2_E       = 999999.99999
							CTRE_G_O_2_cor     = 999999.99999
							AMPL_G_O_2_cor     = 999999.99999
							SGMA_G_O_2_cor     = 999999.99999
							OFST_G_O_2_cor     = 999999.99999
							chisqr_O_2         = 999999.99999
							redchi_O_2         = 999999.99999

							AMPL_SNR_2         = 999999.99999
							CTRE_SNR_2         = 999999.99999
							SGMA_SNR_2         = 999999.99999
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines-1])+'-'+str(LINES[3][lines-1]),'yellow')
						print
						print colored(str(LINES[5][lines-1])+'_CGLC: ' + str(CTRE_G_C_1)  ,'yellow')
						print colored(str(LINES[5][lines-1])+'_AGLC: ' + str(AMPL_G_C_1)  ,'yellow')
						print colored(str(LINES[5][lines-1])+'_SGLC: ' + str(SGMA_G_C_1)  ,'yellow')
						print colored(str(LINES[5][lines-1])+'_FGLC: ' + str(FWHM_G_C_1)  ,'yellow')
						print colored(str(LINES[5][lines-1])+'_WGLC: ' + str(EW_C_1)      ,'yellow')
						print colored(str(LINES[5][lines-1])+'_EGLC: ' + str(EWE_C_1),'yellow')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O

						if fit_vls_hdr == True and fix_ctr_gau_2==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                         ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                         ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLO2',float(CTRE_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALO2',float(AMPL_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLO2',float(FWHM_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLO2',float(EW_O_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELO2',float(EWE_O_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OFO2',float(OFST_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLC2',float(CTRE_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALC2',float(AMPL_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLC2',float(SGMA_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLC2',float(FWHM_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLC2',float(EW_C_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELC2',float(EWE_C_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CEC2',float(CTRE_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AEC2',float(AMPL_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SEC2',float(SGMA_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CEC2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AEC2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SEC2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))


							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHL2',float(chisqr_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 2-1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRL2',float(redchi_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines-1]) + ' Chi2 Reduced 2-1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('The fit (CTR-C & CTR_O) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-C & CTR_O) values will NOT be added to the fits headers!','magenta')
							print

						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values G-2 for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF02',float(L1_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-1-2 Wdt-Fit  1GF-IntVal')         #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WP02',float(L2_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-2-2 Wdt-Plt  1GF-IntVal')         #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF02',float(L7_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-7-2 Ctr Fit Bnds  1GF-IntVal')    #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CO02',float(L8_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-8-2 Ctr Fit Ofst  1GF-IntVal')    #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF02',float(L10_2),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10-2 Amp Fit Bnds  1GF-IntVal')   #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial Values: ','cyan')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')					
							print colored(str(LINES[5][lines])+'_WF02' + ': ' + str(L1_2),'cyan')
							print colored(str(LINES[5][lines])+'_WP02' + ': ' + str(L2_2),'cyan')
							print colored(str(LINES[5][lines])+'_CF02' + ': ' + str(L7_2),'cyan')
							print colored(str(LINES[5][lines])+'_CO02' + ': ' + str(L8_2),'cyan')
							print colored(str(LINES[5][lines])+'_AF02' + ': ' + str(L10_2),'cyan')
							print
						else:
							print
							print colored('Initial Guess Values G-1 for line Fitting will NOT be recorded!','yellow')
							print
							pass							
						#################################################CENTRAL GAUSSIAN-2-C##################################################					
					elif fix_ctr_gau_2 == True:
						print
						print colored('0 CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						try:
							CTRE_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF02')
							AMPL_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF02')
							FWHM_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF02')
							EW_0_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF02')
							EWE_0_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF02')

							CTRE_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLO2')
							AMPL_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALO2')
							FWHM_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLO2')
							EW_O_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLO2')
							EWE_O_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELO2')
							OFST_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFO2')

							CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC2')
							AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC2')
							SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC2')
							FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC2')
							EW_C_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLC2')
							EWE_C_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELC2')
							CTRE_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CEC2')
							AMPL_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AEC2')
							SGMA_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SEC2')

							chisqr_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHL2')
							redchi_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRL2')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines-1])+'_CF01','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau_2=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CGF0')
							AMPL_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_AGF0')
							FWHM_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_FGF0')
							EW_0_2        = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_WGF0')
							EWE_0_2       = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_EGF0')

							CTRE_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CGLO')
							AMPL_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_AGLO')
							FWHM_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_FGLO')
							EW_O_2        = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_WGLO')
							EWE_O_2       = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_EGLO')
							OFST_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_OFSO')

							CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CGLC')
							AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_AGLC')
							SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_SGLC')
							FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_FGLC')
							EW_C_2        = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_WGLC')
							EWE_C_2       = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_EGLC')
							CTRE_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CLEC')
							AMPL_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_ALEC')
							SGMA_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_SLEC')

							chisqr_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CHGL')
							redchi_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CRGL')
					print
					print colored(str(LINES[0][lines-1])+'-'+str(LINES[3][lines-1])+'-CTR','magenta')
					print
					print colored(str(LINES[5][lines-1])+'_CLC2: ' + str(CTRE_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines-1])+'_ALC2: ' + str(AMPL_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines-1])+'_SLC2: ' + str(SGMA_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines-1])+'_FLC2: ' + str(FWHM_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines-1])+'_WLC2: ' + str(EW_C_2)      ,'magenta')
					print colored(str(LINES[5][lines-1])+'_ELC2: ' + str(EWE_C_2),'magenta')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','magenta')
					print colored(str(CTRE_G_C_2)+', '+str(AMPL_G_C_2)+', '+str(SGMA_G_C_2),'magenta')
					print
					##################################################CENTRAL GAUSSIAN-2###################################################
					###############################################COMPUTING TOTAL AREA###############################################
					print colored('Computing Flux Area','yellow')
					###############################################COMPUTING TOTAL AREA###############################################
					CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC1')
					AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC1')
					SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC1')

					CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC2')
					AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC2')
					SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC2')

					print
					print colored('Computing Areas using info from fits headers.','yellow')
					print 'AMPL_G_C_1: ',str(AMPL_G_C_1),'SGMA_G_C_1: ',str(SGMA_G_C_1)
					print 'AMPL_G_C_2: ',str(AMPL_G_C_2),'SGMA_G_C_2: ',str(SGMA_G_C_2)					
					print

					W_C_1     = integrate.quad(lambda x: AMPL_G_C_1*np.exp(-((x)**2)/(2*SGMA_G_C_1**2)), -np.inf, np.inf)
					EW_C_1    = np.round(abs(np.asarray(W_C_1[0])),10)
					EWE_C_1   = np.round(abs(np.asarray(W_C_1[1])),10)

					W_C_2     = integrate.quad(lambda x: AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2)), -np.inf, np.inf)
					EW_C_2    = np.round(abs(np.asarray(W_C_2[0])),10)
					EWE_C_2   = np.round(abs(np.asarray(W_C_2[1])),10)

					W_C       = integrate.quad(lambda x:  AMPL_G_C_1*np.exp(-((x)**2)/(2*SGMA_G_C_1**2)) + 
												AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2)), 
												-np.inf, np.inf
												)
					EW_C      = np.round((np.asarray(W_C[0])),10)
					EWE_C     = np.round(abs(np.asarray(W_C[1])),10)

					if EW_C_1 == 999999.99999:
						EW_C_1 = 0
					else:
						pass
					if EW_C_2 == 999999.99999:
						EW_C_2 = 0
					else:
						pass
						
					print
					print colored('Areas     :','yellow')
					print colored('Area CTR-1: ' + str(EW_C_1),'yellow')
					print colored('Area CTR-2: ' + str(EW_C_2),'yellow')					
					print colored('Area CTR-B: ' + str(EW_C),'yellow')
					print					
					###############################################COMPUTING TOTAL AREA###############################################
					#############################################ADDING AREA TO FTIS HEADER#############################################
					if fit_vls_hdr == True:
						print
						print colored('The Areas values will be updated to the fits headers!','magenta')
						print
						print colored('Area CTR-1: '                  + str(EW_C_1)   + '-' +str(LINES[5][lines])+'_WMC1','yellow')
						print colored('Area CTR-2: '                  + str(EW_C_2)   + '-' +str(LINES[5][lines])+'_WMC2','yellow')
						print colored('Area CTR-B: '                  + str(EW_C)     + '-' +str(LINES[5][lines])+'_WGM1','yellow')
						print						
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC1',float(EW_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C1-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC1',float(EWE_C_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C1-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC2',float(EW_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C2-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC2',float(EWE_C_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C2-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM1',float(EW_C)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C1-C2 Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM1',float(EWE_C)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C1-C2 Crct' + str(fit_type))
					else:
						print
						print colored('The Areas values will NOT be updated to the fits headers!','magenta')
						print
					#############################################ADDING AREA TO FTIS HEADER#############################################
				elif 'Dbl' in LINES[3][lines] and fit_fnct=='gaussM' and fit_type == 'lmfit' and mke_lne_fit == True and uft_lne_vls == False:
					fit_typ = 'GM'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					print 
					print colored('Double Line Fit (Ind)','yellow')
					print colored(LINES[3][lines-2]  + '-' + str(LINES[0][lines-2])   + '-' + str(LINES[1][lines-2]),'cyan')
					print LINES[3][lines]+ '-' + str(LINES[0][lines]) + '-'  + str(LINES[1][lines])
					print colored(LINES[3][lines-1]+ '-' + str(LINES[0][lines-1]) + '-' + str(LINES[1][lines-1]),'magenta')

					from lmfit import Model

					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
					if ivl_fts_hdr == True:
						try:
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial fit variables from fits header!','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WM01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WM01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WM02' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WM02' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print '*****'
							print
							print
							print colored('Initial fit variables from fits header!','yellow')
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print							
							#quit()
						try:
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')							
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')							
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AM01',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AM02',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('*****Success!******','yellow')
							print
						################FOR CASES WHERE KLINE WIIDHT ==0################
						if L1_1  == 0:
							L1_1  = 1
						else:
							pass
						if L1_2  == 0:
							L1_2  = 1
						else:
							pass
						################FOR CASES WHERE KLINE WIIDHT ==0################
					elif ivl_fts_hdr == False:
						L1_1  = LINES[1][lines-2]
						L2_1  = LINES[2][lines-2]
						L7_1  = LINES[7][lines-2]
						L8_1  = LINES[8][lines-2]
						L10_1 = LINES[10][lines-2]

						L1_2  = LINES[1][lines-1]
						L2_2  = LINES[2][lines-1]
						L7_2  = LINES[7][lines-1]
						L8_2  = LINES[8][lines-1]
						L10_2 = LINES[10][lines-1]
						print
						print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
						print					
					print colored('Initial Values: ','cyan')
					print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')					
					print colored(str(LINES[5][lines-2])+'_WF01' + ': ' + str(L1_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_WP01' + ': ' + str(L2_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_CF01' + ': ' + str(L7_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_CO01' + ': ' + str(L8_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_AF01' + ': ' + str(L10_1),'cyan')
					print
					print colored('Initial Values: ','magenta')
					print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')					
					print colored(str(LINES[5][lines-1])+'_WF02' + ': ' + str(L1_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_WP02' + ': ' + str(L2_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_CF02' + ': ' + str(L7_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_CO02' + ': ' + str(L8_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_AF02' + ': ' + str(L10_2),'magenta')
					print
					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
					lmb_min_lim_line_ft = (LINES[0][lines-2]+L8_1) - MSK_NTMS*L1_1 
					lmb_max_lim_line_ft = (LINES[0][lines-2]+L8_1) + MSK_NTMS*L1_1
					lmb_min_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_1#- 20#L2_1 - 10 
					lmb_max_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_1#+ 20#L2_1 + 10

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function

					idx_ctr_ft_reg = np.abs(lambda_glx[mask_ft] - (LINES[0][lines-2]+L8_1)).argmin()

					X0_f2DG    = (LINES[0][lines-2]+L8_1)
					SIGMA_f2DG = L1_1
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))
					A_f2DG     = -(1-inten_glx[mask_ft][idx_ctr_ft_reg])

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
				
					##################################################CENTRAL GAUSSIAN-1###################################################
					if fix_ctr_gau_1 == False:
						print
						print colored('Fitting 1st line','cyan')
						print colored('1-0-Fitting Central line','cyan')
						print
						#################################################CENTRAL GAUSSIAN-1-0##################################################
						try:
							gmodel_0           = Model(func_1D_Gaussian)
							gmodel_0.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG - (X0_f2DG*L7_1), max=X0_f2DG + (X0_f2DG*L7_1))
							gmodel_0.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*L10_1), max=A_f2DG  + (A_f2DG*L10_1))
							gmodel_0.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_0             = gmodel_0.make_params()							
							result_0_1         = gmodel_0.fit(inten_glx[mask_ft],pars_0,
													X=lambda_glx[mask_ft],
													X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_0_1         = result_0_1.params['X_0'].value
							AMPL_G_0_1         = result_0_1.params['A'].value
							SGMA_G_0_1         = abs(result_0_1.params['SIGMA'].value)
							FWHM_G_0_1         = lw_sgma2fwhm(SGMA_G_0_1)
							W_0_1              = integrate.quad(lambda x: AMPL_G_0_1*np.exp(-((x)**2)/(2*SGMA_G_0_1**2)), -np.inf, np.inf)
							EW_0_1             = np.round(abs(np.asarray(W_0_1[0])),10)
							EWE_0_1            = np.round(abs(np.asarray(W_0_1[1])),10)
							data_fitted_0      = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_0_1,AMPL_G_0_1,SGMA_G_0_1)

							CTRE_G_0_1_E       = result_0_1.params['X_0'].stderr
							AMPL_G_0_1_E       = result_0_1.params['A'].stderr
							SGMA_G_0_1_E       = result_0_1.params['SIGMA'].stderr

							CTRE_G_0_1_cor     = result_0_1.params['X_0'].correl
							AMPL_G_0_1_cor     = result_0_1.params['A'].correl
							SGMA_G_0_1_cor     = result_0_1.params['SIGMA'].correl

							chisqr_0_1         = result_0_1.chisqr
							redchi_0_1         = result_0_1.redchi
						except (RuntimeError,ValueError,TypeError):
							popt_0_1, pcov_0_1   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_0_1           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_0_1         = 999999.99999
							AMPL_G_0_1         = 999999.99999
							SGMA_G_0_1         = 999999.99999
							FWHM_G_0_1         = 999999.99999
							EW_0_1             = 999999.99999
							EWE_0_1            = 999999.99999

							CTRE_G_0_1_E      = 999999.99999
							AMPL_G_0_1_E      = 999999.99999
							SGMA_G_0_1_E      = 999999.99999

							CTRE_G_0_1_cor    = 999999.99999
							AMPL_G_0_1_cor    = 999999.99999
							SGMA_G_0_1_cor    = 999999.99999

							chisqr_0_1        = 999999.99999
							redchi_0_1        = 999999.99999
						if fit_vls_hdr == True and fix_ctr_gau_1 ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM01',float(CTRE_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM01',float(AMPL_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FM01',float(FWHM_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM01',float(EW_0_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EM01',float(EWE_0_1)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('The fit (CTR-1-0) values will be added to the fits headers!','cyan')
							print
						else:
							print
							print colored('The fit (CTR-1-0) values will NOT be added to the fits headers!','cyan')
							print
						#################################################CENTRAL GAUSSIAN-1-0##################################################
						#################################################CENTRAL GAUSSIAN-1-C##################################################
						try:
							gmodel_O           = Model(func_1D_Gaussian_O)
							gmodel_O.set_param_hint('X_0'   , value=X0_f2DG  , min=X0_f2DG - (X0_f2DG*L7_1), max=X0_f2DG + (X0_f2DG*L7_1))
							gmodel_O.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*L10_1), max=A_f2DG  + (A_f2DG*L10_1))
							gmodel_O.set_param_hint('SIGMA' , value=SIGMA_f2DG)
							gmodel_O.set_param_hint('OFFSET', value=max(inten_glx[mask_ft])-1)
							pars_O             = gmodel_O.make_params()

							result_O_1         = gmodel_O.fit(inten_glx[mask_ft],pars_O,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,OFFSET=max(inten_glx[mask_ft])-1,
													nan_policy = 'omit')
							CTRE_G_O_1         = result_O_1.params['X_0'].value
							AMPL_G_O_1         = result_O_1.params['A'].value
							SGMA_G_O_1         = abs(result_O_1.params['SIGMA'].value)
							OFST_G_O_1         = abs(result_O_1.params['OFFSET'].value)
							FWHM_G_O_1         = lw_sgma2fwhm(SGMA_G_O_1)
							W_O_1              = integrate.quad(lambda x: AMPL_G_O_1*np.exp(-((x)**2)/(2*SGMA_G_O_1**2)), -np.inf, np.inf)
							EW_O_1             = np.round(abs(np.asarray(W_O_1[0])),10)
							EWE_O_1            = np.round(abs(np.asarray(W_O_1[1])),10)
							data_fitted_O_1    = func_1D_Gaussian_O((lambda_glx[mask_ft]),CTRE_G_O_1,AMPL_G_O_1,SGMA_G_O_1,OFST_G_O_1)

							CTRE_G_O_E         = result_O_1.params['X_0'].stderr
							AMPL_G_O_E         = result_O_1.params['A'].stderr
							SGMA_G_O_E         = result_O_1.params['SIGMA'].stderr

							CTRE_G_O_cor       = result_O_1.params['X_0'].correl
							AMPL_G_O_cor       = result_O_1.params['A'].correl
							SGMA_G_O_cor       = result_O_1.params['SIGMA'].correl

							chisqr_O_1         = result_O_1.chisqr
							redchi_O_1         = result_O_1.redchi
							
							#####################################################################################################################
							if ofs_ctr_fit == True:
								print
								print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
								print
								lmb_min_lim_line_ft = (CTRE_G_O_1-L8_1) - MSK_NTMS*L1_1 
								lmb_max_lim_line_ft = (CTRE_G_O_1+L8_1) + MSK_NTMS*L1_1
								mask_ft     = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
							else:
								pass
							#####################################################################################################################
							initial_guess_C    = (X0_f2DG,A_f2DG,SIGMA_f2DG)

							gmodel_C           = Model(func_1D_Gaussian)
							gmodel_C.set_param_hint('X_0'  , value=X0_f2DG, min=X0_f2DG - (X0_f2DG*(X0_f2DG*L7_1)), max=X0_f2DG + (X0_f2DG*(X0_f2DG*L7_1)))
							gmodel_C.set_param_hint('A'    , value=A_f2DG , min=A_f2DG  - (A_f2DG*(A_f2DG*L10_1)) , max=A_f2DG  + (A_f2DG*(A_f2DG*L10_1)))
							gmodel_C.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_C             = gmodel_C.make_params()
							result_C_1         = gmodel_C.fit(inten_glx[mask_ft],pars_C,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_C_1         = result_C_1.params['X_0'].value
							AMPL_G_C_1         = result_C_1.params['A'].value
							SGMA_G_C_1         = abs(result_C_1.params['SIGMA'].value)
							FWHM_G_C_1         = lw_sgma2fwhm(SGMA_G_C_1)
							W_C_1              = integrate.quad(lambda x: AMPL_G_C_1*np.exp(-((x)**2)/(2*SGMA_G_C_1**2)), -np.inf, np.inf)
							EW_C_1             = np.round(abs(np.asarray(W_C_1[0])),10)
							EWE_C_1            = np.round(abs(np.asarray(W_C_1[1])),10)
							data_fitted_C_1    = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_C_1,AMPL_G_C_1,SGMA_G_C_1)

							CTRE_G_C_1_E         = result_C_1.params['X_0'].stderr
							AMPL_G_C_1_E         = result_C_1.params['A'].stderr
							SGMA_G_C_1_E         = result_C_1.params['SIGMA'].stderr

							CTRE_G_C_1_cor       = result_C_1.params['X_0'].correl
							AMPL_G_C_1_cor       = result_C_1.params['A'].correl
							SGMA_G_C_1_cor       = result_C_1.params['SIGMA'].correl

							AMPL_SNR_1           = AMPL_G_C_1
							CTRE_SNR_1           = CTRE_G_C_1
							SGMA_SNR_1           = abs(SGMA_G_C_1)

							if CTRE_G_C_1_E == None:
								CTRE_G_C_1_E = 999999.99999
							else:
								pass
							if AMPL_G_C_1_E == None:
								AMPL_G_C_1_E = 999999.99999
							else:
								pass
							if SGMA_G_C_1_E == None:
								SGMA_G_C_1_E = 999999.99999
							else:
								pass
							if CTRE_G_C_1_cor == None:
								CTRE_G_C_1_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_1_cor == None:
								AMPL_G_C_1_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_1_cor == None:
								SGMA_G_C_1_cor = 999999.99999
							else:
								pass
							chisqr_C_1      = result_C_1.chisqr
							redchi_C_1      = result_C_1.redchi
							##inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O #OFFSET +
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C_1, pcov_C_1 = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_1           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_1         = 999999.99999
							AMPL_G_C_1         = 999999.99999
							SGMA_G_C_1         = 999999.99999
							FWHM_G_C_1         = 999999.99999
							EW_C_1             = 999999.99999
							EWE_C_1            = 999999.99999

							CTRE_G_C_1_E       = 999999.99999
							AMPL_G_C_1_E       = 999999.99999
							SGMA_G_C_1_E       = 999999.99999
							CTRE_G_C_1_cor     = 999999.99999
							AMPL_G_C_1_cor     = 999999.99999
							SGMA_G_C_1_cor     = 999999.99999
							chisqr_C_1         = 999999.99999
							redchi_C_1         = 999999.99999

							popt_O_1 ,pcov_O_1 = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
							perr_O_1           = [999999.99999,999999.99999,999999.99999,999999.99999]
							CTRE_G_O_1         = 999999.99999
							AMPL_G_O_1         = 999999.99999
							SGMA_G_O_1         = 999999.99999
							OFST_G_O_1         = 999999.99999
							FWHM_G_O_1         = 999999.99999
							EW_O_1             = 999999.99999
							EWE_O_1            = 999999.99999

							CTRE_G_O_1_E      = 999999.99999
							AMPL_G_O_1_E      = 999999.99999
							SGMA_G_O_1_E      = 999999.99999
							CTRE_G_O_1_cor    = 999999.99999
							AMPL_G_O_1_cor    = 999999.99999
							SGMA_G_O_1_cor    = 999999.99999
							OFST_G_O_1_cor    = 999999.99999
							chisqr_O_1        = 999999.99999
							redchi_O_1        = 999999.99999

							AMPL_SNR_1        = 999999.99999
							CTRE_SNR_1        = 999999.99999
							SGMA_SNR_1        = 999999.99999
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines-2])+'-'+str(LINES[3][lines-2]),'cyan')
						print
						print colored(str(LINES[5][lines-2])+'_CGLC: ' + str(CTRE_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_AGLC: ' + str(AMPL_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_SGLC: ' + str(SGMA_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_FGLC: ' + str(FWHM_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_WGLC: ' + str(EW_C_1)      ,'cyan')
						print colored(str(LINES[5][lines-2])+'_EGLC: ' + str(EWE_C_1)     ,'cyan')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O
						if fit_vls_hdr == True and fix_ctr_gau_1==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                         ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                         ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CMO1',float(CTRE_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AMO1',float(AMPL_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FMO1',float(FWHM_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WMO1',float(EW_O_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EMO1',float(EWE_O_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OMO1',float(OFST_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CMC1',float(CTRE_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AMC1',float(AMPL_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SMC1',float(SGMA_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FMC1',float(FWHM_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC1',float(EW_C_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC1',float(EWE_C_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CM1E',float(CTRE_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AM1E',float(AMPL_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SM1E',float(SGMA_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CM1E',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AM1E',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SM1E',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))


							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHM1',float(chisqr_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRM1',float(redchi_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines-2]) + ' Chi2 Reduced 1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('The fit (CTR-1-C & CTR-1-O) values will be added to the fits headers!','cyan')
							print
						else:
							print
							print colored('The fit (CTR-1-C & CTR-1-O) values will NOT be added to the fits headers!','cyan')
							print

						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values G-1 for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM01',float(L1_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-1-1 Wdt-Fit  1GF-IntVal')         #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM01',float(L2_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-2-1 Wdt-Plt  1GF-IntVal')         #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM01',float(L7_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-7-1 Ctr Fit Bnds  1GF-IntVal')    #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM01',float(L8_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-8-1 Ctr Fit Ofst  1GF-IntVal')    #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM01',float(L10_1),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10-1 Amp Fit Bnds  1GF-IntVal')   #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial Values: ','cyan')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')					
							print colored(str(LINES[5][lines])+'_WM01' + ': ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WM01' + ': ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ': ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ': ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AM01' + ': ' + str(L10_1),'cyan')
							print
						else:
							print
							print colored('Initial Guess Values G-1 for line Fitting will NOT be recorded!','yellow')
							print
							pass
					#################################################CENTRAL GAUSSIAN-1-C##################################################					
					elif fix_ctr_gau_1 == True:
						print
						print colored('1-CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						try:
							CTRE_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CM01')
							AMPL_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AM01')
							FWHM_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FM01')
							EW_0_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WM01')
							EWE_0_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EM01')

							CTRE_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMO1')
							AMPL_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMO1')
							FWHM_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMO1')
							EW_O_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WMO1')
							EWE_O_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EMO1')
							OFST_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OMO1')

							CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC1')
							AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC1')
							SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC1')
							FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMC1')
							EW_C_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WMC1')
							EWE_C_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EMC1')
							CTRE_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC1')
							AMPL_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC1')
							SGMA_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC1')

							chisqr_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHM1')
							redchi_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRM1')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines-2])+'_CF01','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau_1=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF01')
							AMPL_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF01')
							FWHM_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF01')
							EW_0_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF01')
							EWE_0_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF01')

							CTRE_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLO1')
							AMPL_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALO1')
							FWHM_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLO1')
							EW_O_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLO1')
							EWE_O_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELO1')
							OFST_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFO1')

							CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC1')
							AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC1')
							SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC1')
							FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC1')
							EW_C_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLC1')
							EWE_C_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELC1')
							CTRE_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CEC1')
							AMPL_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AEC1')
							SGMA_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SEC1')

							chisqr_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHL1')
							redchi_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRL1')
					print
					print colored(str(LINES[0][lines-2])+'-'+str(LINES[3][lines-2])+'-CTR','cyan')
					print colored('From: '+str(LINES[3][lines])+'-CTR','cyan')
					print
					print colored(str(LINES[5][lines])+'_CMC1: ' + str(CTRE_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines])+'_AMC1: ' + str(AMPL_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines])+'_SMC1: ' + str(SGMA_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines])+'_FMC1: ' + str(FWHM_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines])+'_WMC1: ' + str(EW_C_1)      ,'cyan')
					print colored(str(LINES[5][lines])+'_EMC1: ' + str(EWE_C_1),'cyan')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','cyan')
					print colored(str(CTRE_G_C_1)+', '+str(AMPL_G_C_1)+', '+str(SGMA_G_C_1),'cyan')
					print
					##################################################CENTRAL GAUSSIAN-1###################################################
					lmb_min_lim_line_ft = (LINES[0][lines-1]+L8_2) - MSK_NTMS*LINES[1][lines-1] 
					lmb_max_lim_line_ft = (LINES[0][lines-1]+L8_2) + MSK_NTMS*LINES[1][lines-1]
					lmb_min_lim_line    = (LINES[0][lines-1]+L8_2)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_2
					lmb_max_lim_line    = (LINES[0][lines-1]+L8_2)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_2

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function
					idx_ctr_ft_reg = np.abs(lambda_glx[mask_ft] - (LINES[0][lines-1]+L8_2)).argmin()

					X0_f2DG    = (LINES[0][lines-1]+L8_2)
					SIGMA_f2DG = LINES[1][lines-1]
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))
					A_f2DG     = -(1-inten_glx[mask_ft][idx_ctr_ft_reg])

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
					##################################################CENTRAL GAUSSIAN-2###################################################
					if fix_ctr_gau_2 == False:
						print
						print colored('Fitting 2nd line','magenta')
						print colored('2-0-Fitting Central line','magenta')
						print
						#################################################CENTRAL GAUSSIAN-2-0##################################################
						try:
							gmodel_0           = Model(func_1D_Gaussian)
							gmodel_0.set_param_hint('X_0'  , value=X0_f2DG , min=X0_f2DG - (X0_f2DG*L7_2), max=X0_f2DG + (X0_f2DG*L7_2))
							gmodel_0.set_param_hint('A'    , value=A_f2DG  , min=A_f2DG  - (A_f2DG*L10_2), max=A_f2DG  + (A_f2DG*L10_2))
							gmodel_0.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_0             = gmodel_0.make_params()							
							result_0_2         = gmodel_0.fit(inten_glx[mask_ft],pars_0,
													X=lambda_glx[mask_ft],
													X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_0_2         = result_0_2.params['X_0'].value
							AMPL_G_0_2         = result_0_2.params['A'].value
							SGMA_G_0_2         = abs(result_0_2.params['SIGMA'].value)
							FWHM_G_0           = lw_sgma2fwhm(SGMA_G_0_2)
							W_0_2              = integrate.quad(lambda x: AMPL_G_0_2*np.exp(-((x)**2)/(2*SGMA_G_0_2**2)), -np.inf, np.inf)
							EW_0_2             = np.round(abs(np.asarray(W_0_2[0])),10)
							EWE_0_2            = np.round(abs(np.asarray(W_0_2[1])),10)
							data_fitted_0_2    = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_0_2,AMPL_G_0_2,SGMA_G_0_2)

							CTRE_G_0_2_E       = result_0_2.params['X_0'].stderr
							AMPL_G_0_2_E       = result_0_2.params['A'].stderr
							SGMA_G_0_2_E       = result_0_2.params['SIGMA'].stderr

							CTRE_G_0_2_cor     = result_0_2.params['X_0'].correl
							AMPL_G_0_2_cor     = result_0_2.params['A'].correl
							SGMA_G_0_2_cor     = result_0_2.params['SIGMA'].correl

							chisqr_0           = result_0_2.chisqr
							redchi_0           = result_0_2.redchi
						except (RuntimeError,ValueError,TypeError):
							popt_0, pcov_0   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_0           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_0_2         = 999999.99999
							AMPL_G_0_2         = 999999.99999
							SGMA_G_0_2         = 999999.99999
							FWHM_G_0         = 999999.99999
							EW_0_2             = 999999.99999
							EWE_0_2            = 999999.99999

							CTRE_G_0_2_E      = 999999.99999
							AMPL_G_0_2_E      = 999999.99999
							SGMA_G_0_2_E      = 999999.99999

							CTRE_G_0_2_cor    = 999999.99999
							AMPL_G_0_2_cor    = 999999.99999
							SGMA_G_0_2_cor    = 999999.99999

							chisqr_0        = 999999.99999
							redchi_0        = 999999.99999
						if fit_vls_hdr == True and fix_ctr_gau_2 ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM02',float(CTRE_G_0_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM02',float(AMPL_G_0_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FM02',float(FWHM_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM02',float(EW_0_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EM02',float(EWE_0_2)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('The fit (CTR-2-0) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-2-0) values will NOT be added to the fits headers!','magenta')
							print
						#################################################CENTRAL GAUSSIAN-2-0##################################################
						#################################################CENTRAL GAUSSIAN-2-C##################################################
						try:
							gmodel_O           = Model(func_1D_Gaussian_O)
							gmodel_O.set_param_hint('X_0'   , value=X0_f2DG , min=X0_f2DG - (X0_f2DG*L7_2), max=X0_f2DG + (X0_f2DG*L7_2))
							gmodel_O.set_param_hint('A'     , value=A_f2DG  , min=A_f2DG  - (A_f2DG*L10_2), max=A_f2DG  + (A_f2DG*L10_2))
							gmodel_O.set_param_hint('SIGMA' , value=SIGMA_f2DG)
							gmodel_O.set_param_hint('OFFSET', value=max(inten_glx[mask_ft])-1)
							pars_O             = gmodel_O.make_params()

							result_O_2         = gmodel_O.fit(inten_glx[mask_ft],pars_O,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,OFFSET=max(inten_glx[mask_ft])-1,
													nan_policy = 'omit')
							CTRE_G_O_2         = result_O_2.params['X_0'].value
							AMPL_G_O_2         = result_O_2.params['A'].value
							SGMA_G_O_2         = abs(result_O_2.params['SIGMA'].value)
							OFST_G_O_2         = abs(result_O_2.params['OFFSET'].value)
							FWHM_G_O_2         = lw_sgma2fwhm(SGMA_G_O_2)
							W_O_2              = integrate.quad(lambda x: AMPL_G_O_2*np.exp(-((x)**2)/(2*SGMA_G_O_2**2)), -np.inf, np.inf)
							EW_O_2             = np.round(abs(np.asarray(W_O_2[0])),10)
							EWE_O_2            = np.round(abs(np.asarray(W_O_2[1])),10)
							data_fitted_O_2    = func_1D_Gaussian_O((lambda_glx[mask_ft]),CTRE_G_O_2,AMPL_G_O_2,SGMA_G_O_2,OFST_G_O_2)

							CTRE_G_O_2_E       = result_O_2.params['X_0'].stderr
							AMPL_G_O_2_E       = result_O_2.params['A'].stderr
							SGMA_G_O_E         = result_O_2.params['SIGMA'].stderr

							CTRE_G_O_2_cor     = result_O_2.params['X_0'].correl
							AMPL_G_O_2_cor     = result_O_2.params['A'].correl
							SGMA_G_O_cor       = result_O_2.params['SIGMA'].correl

							chisqr_O_2         = result_O_2.chisqr
							redchi_O_2         = result_O_2.redchi
							
							#####################################################################################################################
							if ofs_ctr_fit == True:
								print
								print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
								print
								lmb_min_lim_line_ft = (CTRE_G_O_2-L8_2) - MSK_NTMS*LINES[1][lines-1] 
								lmb_max_lim_line_ft = (CTRE_G_O_2+L8_2) + MSK_NTMS*LINES[1][lines-1]
								mask_ft     = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
							else:
								pass
							#####################################################################################################################
							initial_guess_C    = (X0_f2DG,A_f2DG,SIGMA_f2DG)

							gmodel_C           = Model(func_1D_Gaussian)
							gmodel_C.set_param_hint('X_0'  , value=X0_f2DG, min=X0_f2DG - (X0_f2DG*(X0_f2DG*L7_2)), max=X0_f2DG + (X0_f2DG*(X0_f2DG*L7_2)))
							gmodel_C.set_param_hint('A'    , value=A_f2DG , min=A_f2DG  - (A_f2DG*(A_f2DG*L10_2)) , max=A_f2DG  + (A_f2DG*(A_f2DG*L10_2)))
							gmodel_C.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_C             = gmodel_C.make_params()
							result_C_2         = gmodel_C.fit(inten_glx[mask_ft],pars_C,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_C_2         = result_C_2.params['X_0'].value
							AMPL_G_C_2         = result_C_2.params['A'].value
							SGMA_G_C_2         = abs(result_C_2.params['SIGMA'].value)
							FWHM_G_C_2         = lw_sgma2fwhm(SGMA_G_C_2)

							W_C_2              = integrate.quad(lambda x: AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2)), -np.inf, np.inf)
							EW_C_2             = np.round(abs(np.asarray(W_C_2[0])),10)
							EWE_C_2            = np.round(abs(np.asarray(W_C_2[1])),10)
							data_fitted_C      = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_C_2,AMPL_G_C_2,SGMA_G_C_2)

							CTRE_G_C_2_E       = result_C_2.params['X_0'].stderr
							AMPL_G_C_2_E       = result_C_2.params['A'].stderr
							SGMA_G_C_2_E       = result_C_2.params['SIGMA'].stderr

							CTRE_G_C_2_cor     = result_C_2.params['X_0'].correl
							AMPL_G_C_2_cor     = result_C_2.params['A'].correl
							SGMA_G_C_2_cor     = result_C_2.params['SIGMA'].correl

							AMPL_SNR_2           = AMPL_G_C_2
							CTRE_SNR_2           = CTRE_G_C_2
							SGMA_SNR_2           = abs(SGMA_G_C_2)

							if CTRE_G_C_2_E == None:
								CTRE_G_C_2_E = 999999.99999
							else:
								pass
							if AMPL_G_C_2_E == None:
								AMPL_G_C_2_E = 999999.99999
							else:
								pass
							if SGMA_G_C_2_E == None:
								SGMA_G_C_2_E = 999999.99999
							else:
								pass
							if CTRE_G_C_2_cor == None:
								CTRE_G_C_2_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_2_cor == None:
								AMPL_G_C_2_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_2_cor == None:
								SGMA_G_C_2_cor = 999999.99999
							else:
								pass
							chisqr_C_2        = result_C_2.chisqr
							redchi_C_2        = result_C_2.redchi
							##inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O #OFFSET +
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C_2, pcov_C_2 = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_2           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_2         = 999999.99999
							AMPL_G_C_2         = 999999.99999
							SGMA_G_C_2         = 999999.99999
							FWHM_G_C_2         = 999999.99999
							EW_C_2             = 999999.99999
							EWE_C_2            = 999999.99999

							CTRE_G_C_2_E       = 999999.99999
							AMPL_G_C_2_E       = 999999.99999
							SGMA_G_C_2_E       = 999999.99999
							CTRE_G_C_2_cor     = 999999.99999
							AMPL_G_C_2_cor     = 999999.99999
							SGMA_G_C_2_cor     = 999999.99999
							chisqr_C_2         = 999999.99999
							redchi_C_2         = 999999.99999

							popt_O_2 ,pcov_O_2 = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
							perr_O_2           = [999999.99999,999999.99999,999999.99999,999999.99999]
							CTRE_G_O_2         = 999999.99999
							AMPL_G_O_2         = 999999.99999
							SGMA_G_O_2         = 999999.99999
							OFST_G_O_2         = 999999.99999
							FWHM_G_O_2         = 999999.99999
							EW_O_2             = 999999.99999
							EWE_O_2            = 999999.99999

							CTRE_G_O_2_E       = 999999.99999
							AMPL_G_O_2_E       = 999999.99999
							SGMA_G_O_2_E       = 999999.99999
							CTRE_G_O_2_cor     = 999999.99999
							AMPL_G_O_2_cor     = 999999.99999
							SGMA_G_O_2_cor     = 999999.99999
							OFST_G_O_2_cor     = 999999.99999
							chisqr_O_2         = 999999.99999
							redchi_O_2         = 999999.99999

							AMPL_SNR_2         = 999999.99999
							CTRE_SNR_2         = 999999.99999
							SGMA_SNR_2         = 999999.99999
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines-1])+'-'+str(LINES[3][lines-1]),'yellow')
						print
						print colored(str(LINES[5][lines])+'_CGLC: ' + str(CTRE_G_C_2)  ,'yellow')
						print colored(str(LINES[5][lines])+'_AGLC: ' + str(AMPL_G_C_2)  ,'yellow')
						print colored(str(LINES[5][lines])+'_SGLC: ' + str(SGMA_G_C_2)  ,'yellow')
						print colored(str(LINES[5][lines])+'_FGLC: ' + str(FWHM_G_C_2)  ,'yellow')
						print colored(str(LINES[5][lines])+'_WGLC: ' + str(EW_C_2)      ,'yellow')
						print colored(str(LINES[5][lines])+'_EGLC: ' + str(EWE_C_2),'yellow')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O

						if fit_vls_hdr == True and fix_ctr_gau_2==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                         ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                         ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CMO2',float(CTRE_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AMO2',float(AMPL_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FMO2',float(FWHM_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WMO2',float(EW_O_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EMO2',float(EWE_O_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OMO2',float(OFST_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CMC2',float(CTRE_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AMC2',float(AMPL_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SMC2',float(SGMA_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FMC2',float(FWHM_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC2',float(EW_C_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC2',float(EWE_C_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CM2E',float(CTRE_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AM2E',float(AMPL_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SM2E',float(SGMA_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CM2E',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AM2E',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SM2E',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHM2',float(chisqr_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 2-1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRM2',float(redchi_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines-1]) + ' Chi2 Reduced 2-1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('The fit (CTR-2-C & CTR-2-O) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-2-C & CTR-2-O) values will NOT be added to the fits headers!','magenta')
							print

						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values G-2 for line Fitting will be recorded!','magenta')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','magenta')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM02',float(L1_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-1-2 Wdt-Fit  1GF-IntVal')         #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM02',float(L2_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-2-2 Wdt-Plt  1GF-IntVal')         #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM02',float(L7_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-7-2 Ctr Fit Bnds  1GF-IntVal')    #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM02',float(L8_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-8-2 Ctr Fit Ofst  1GF-IntVal')    #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM02',float(L10_2),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10-2 Amp Fit Bnds  1GF-IntVal')   #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial Values: ','magenta')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')					
							print colored(str(LINES[5][lines])+'_WM02' + ': ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WM02' + ': ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ': ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ': ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AM02' + ': ' + str(L10_2),'magenta')
							print
						else:
							print
							print colored('Initial Guess Values G-2 for line Fitting will NOT be recorded!','yellow')
							print
							pass							
						#################################################CENTRAL GAUSSIAN-2-C##################################################					
					elif fix_ctr_gau_2 == True:
						print
						print colored('2-CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						print colored(LINES[3][lines],'yellow')
						print
						try:
							CTRE_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CM02')
							AMPL_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AM02')
							FWHM_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FM02')
							EW_0_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WM02')
							EWE_0_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EM02')

							CTRE_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMO2')
							AMPL_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMO2')
							FWHM_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMO2')
							EW_O_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WMO2')
							EWE_O_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EMO2')
							OFST_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OMO2')

							CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC2')
							AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC2')
							SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC2')
							FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMC2')
							EW_C_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WMC2')
							EWE_C_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EMC2')
							CTRE_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC2')
							AMPL_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC2')
							SGMA_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC2')

							chisqr_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHM2')
							redchi_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRM2')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines-1])+'_CF01','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau_2=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF02')
							AMPL_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF02')
							FWHM_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF02')
							EW_0_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF02')
							EWE_0_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF02')

							CTRE_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLO2')
							AMPL_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALO2')
							FWHM_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLO2')
							EW_O_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLO2')
							EWE_O_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELO2')
							OFST_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFO2')

							CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC2')
							AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC2')
							SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC2')
							FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC2')
							EW_C_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLC2')
							EWE_C_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELC2')
							CTRE_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CEC2')
							AMPL_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AEC2')
							SGMA_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SEC2')

							chisqr_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHL2')
							redchi_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRL2')
					print
					print colored(str(LINES[0][lines-1])+'-'+str(LINES[3][lines-1])+'-CTR','magenta')
					print colored('From '+str(LINES[3][lines])+'-CTR','cyan')
					print
					print colored(str(LINES[5][lines])+'_CMC2: ' + str(CTRE_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines])+'_AMC2: ' + str(AMPL_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines])+'_SMC2: ' + str(SGMA_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines])+'_FMC2: ' + str(FWHM_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines])+'_WMC2: ' + str(EW_C_2)      ,'magenta')
					print colored(str(LINES[5][lines])+'_EMC2: ' + str(EWE_C_2),'magenta')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','magenta')
					print colored(str(CTRE_G_C_2)+', '+str(AMPL_G_C_2)+', '+str(SGMA_G_C_2),'magenta')
					print
					##################################################CENTRAL GAUSSIAN-2###################################################	
					CTRE_G_0 = CTRE_G_0_1
					#####################################################PRE GAUSSIAN#################################################
					if fix_pre_gau == False and pst_shf_lim>0:
						###########################################DEFINING PRE-PST REGIONS##################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
							lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						else:
							print
							print colored('Using Expected Line Center to Define Fitting Region!','yellow')
							print
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
							
						print
						print colored('Region limits for fitting.','yellow')
						print colored('Line Center    : ' + str(LINES[0][lines]),'cyan')
						print colored('pre_shf_ctr    : ' + str(pre_shf_ctr),'cyan')
						print colored('pre_shf_lim    : ' + str(pre_shf_lim),'cyan')
						print colored('Line Center    : ' + str(LINES[0][lines] - pre_shf_ctr),'cyan')
						print colored('Lower Limit    : ' + str(LINES[0][lines] - (pre_shf_ctr-pre_shf_lim)),'cyan')
						print colored('Upper Limit    : ' + str(LINES[0][lines] - (pre_shf_ctr+pre_shf_lim)),'cyan')
						print 
						print colored('Central    : ' + str(LINES[0][lines]),'cyan')
						print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'cyan')
						print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'cyan')
						print
						print colored('Limits:','cyan')
						print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
						print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr - pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr + pre_shf_lim
						print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr - pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr + pst_shf_lim
						print
						print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr - pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr + pst_shf_lim
						print
						###########################################DEFINING PRE-PST REGIONS##################################################
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						x_a = lambda_glx[mask_ft_pre][X0_f2DG_indx_PRE]
						y_a = inten_glx[mask_ft_pre][X0_f2DG_indx_PRE]
						try:
							print
							print colored('1-Fitting gaussian before line','cyan')
							print
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							A_f2DG_PRE             = (max(inten_glx[mask_ft_pre])-1)
							X0_f2DG_PRE            = LINES[0][lines] - pre_shf_ctr
							SIGMA_f2DG_PRE         = SIGMA_f2DG/2.5
							initial_guess_C_PRE    = (X0_f2DG_PRE,A_f2DG_PRE,SIGMA_f2DG_PRE)
							
							x_a = lambda_glx[mask_ft_pre][X0_f2DG_indx_PRE]
							y_a = inten_glx[mask_ft_pre][X0_f2DG_indx_PRE]
							
							print
							print colored('Initial Guess Values PRE : ','cyan')
							print colored(initial_guess_C_PRE,'cyan')
							print X0_f2DG,pre_shf_ctr,X0_f2DG-pre_shf_ctr
							print
							
							gmodel_C_PRE           = Model(func_1D_Gaussian_Emm)
							gmodel_C_PRE.set_param_hint('X_0'  , value=X0_f2DG_PRE , min=X0_f2DG_PRE-(X0_f2DG_PRE*LINES[7][lines]), max=X0_f2DG_PRE+(X0_f2DG_PRE*LINES[7][lines]))
							gmodel_C_PRE.set_param_hint('A'    , value=A_f2DG_PRE  , min=A_f2DG_PRE -(A_f2DG_PRE*LINES[10][lines]), max=A_f2DG_PRE +(A_f2DG_PRE*LINES[10][lines]))
							gmodel_C_PRE.set_param_hint('SIGMA', value=SIGMA_f2DG_PRE)
							pars_C_PRE             = gmodel_C_PRE.make_params()
							result_C_PRE           = gmodel_C_PRE.fit(inten_glx[mask_ft_pre],pars_C_PRE,
													X=lambda_glx[mask_ft_pre],X_0=X0_f2DG_PRE,A=A_f2DG_PRE,SIGMA=SIGMA_f2DG_PRE,
													nan_policy = 'omit')
							CTRE_G_C_PRE           = result_C_PRE.params['X_0'].value
							AMPL_G_C_PRE           = result_C_PRE.params['A'].value
							SGMA_G_C_PRE           = abs(result_C_PRE.params['SIGMA'].value)
							FWHM_G_C_PRE           = lw_sgma2fwhm(SGMA_G_C_PRE)

							W_C_PRE                = integrate.quad(lambda x: AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2)), -np.inf, np.inf)
							EW_C_PRE               = np.round(abs(np.asarray(W_C_PRE[0])),10)
							EWE_C_PRE              = np.round(abs(np.asarray(W_C_PRE[1])),10)
							data_fitted_C_PRE      = func_1D_Gaussian_Emm((lambda_glx[mask_ft_pre]), CTRE_G_C_PRE,AMPL_G_C_PRE,SGMA_G_C_PRE)
							
							CTRE_G_C_PRE_E         = result_C_PRE.params['X_0'].stderr
							AMPL_G_C_PRE_E         = result_C_PRE.params['A'].stderr
							SGMA_G_C_PRE_E         = result_C_PRE.params['SIGMA'].stderr
							
							CTRE_G_C_PRE_cor       = result_C_PRE.params['X_0'].correl
							AMPL_G_C_PRE_cor       = result_C_PRE.params['A'].correl
							SGMA_G_C_PRE_cor       = result_C_PRE.params['SIGMA'].correl
							
							AMPL_SNR               = AMPL_G_C_PRE
							CTRE_SNR               = CTRE_G_C_PRE
							SGMA_SNR               = abs(SGMA_G_C_PRE)
							
							if CTRE_G_C_PRE_E == None or np.isnan(CTRE_G_C_PRE_E):
								CTRE_G_C_PRE_E = 999999.99999
							else:
								pass
							if AMPL_G_C_PRE_E == None:
								AMPL_G_C_PRE_E = 999999.99999
							else:
								pass
							if SGMA_G_C_PRE_E == None:
								SGMA_G_C_PRE_E = 999999.99999
							else:
								pass
							if CTRE_G_C_PRE_cor == None:
								CTRE_G_C_PRE_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_PRE_cor == None:
								AMPL_G_C_PRE_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_PRE_cor == None:
								SGMA_G_C_PRE_cor = 999999.99999
							else:
								pass
							chisqr_C_PRE           = result_C_PRE.chisqr
							redchi_C_PRE           = result_C_PRE.redchi
							
							
							W_C_PR1    = integrate.quad(lambda x: AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2)), x_a, np.inf)
							EW_C_PR1   = np.round(abs(np.asarray(W_C_PR1[0])),10)
							EWE_C_PR1  = np.round(abs(np.asarray(W_C_PR1[1])),10)
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C_PRE, pcov_C_PRE  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_PRE          = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_PRE        = 999999.99999
							AMPL_G_C_PRE        = 999999.99999
							SGMA_G_C_PRE        = 999999.99999
							FWHM_G_C_PRE        = 999999.99999
							EW_C_PRE            = 999999.99999
							EWE_C_PRE           = 999999.99999
							
							CTRE_G_C_PRE_E      = 999999.99999
							AMPL_G_C_PRE_E      = 999999.99999
							SGMA_G_C_PRE_E      = 999999.99999
							CTRE_G_C_PRE_cor    = 999999.99999
							AMPL_G_C_PRE_cor    = 999999.99999
							SGMA_G_C_PRE_cor    = 999999.99999
							chisqr_C_PRE        = 999999.99999
							redchi_C_PRE        = 999999.99999
							
						if fit_vls_hdr == True and fix_pre_gau == False:						
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHML',float(chisqr_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRML',float(redchi_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 Reduced 1GF' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGM1',float(CTRE_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGM1',float(AMPL_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGM1',float(SGMA_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGM1',float(FWHM_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM1',float(EW_C_PRE)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM1',float(EWE_C_PRE)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CME1',float(CTRE_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AME1',float(AMPL_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SME1',float(SGMA_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							except:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CME1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AME1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SME1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))

							
							Header_Add(specfile_glx,str(LINES[5][lines])+'_XAM1',float(x_a),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PRE GAU-LNR X1 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YAM1',float(y_a),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PRE GAU-LNR Y1 COO')
							print
							print colored('The fit (PRE) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (PRE) values will NOT be added to the fits headers!','magenta')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GSM1',float(pre_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for PRE G-1')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GCM1',float(pre_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for PRE G-1')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass
					elif fix_pre_gau == True or (fix_pre_gau == False and pre_shf_lim<=0):
						try:
							print
							print colored('1 PRE-gaussian already fitted!','yellow')
							print colored('Values from fits header','yellow')
							CTRE_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM1')
							AMPL_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM1')
							SGMA_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM1')
							FWHM_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM1')
							EW_C_PRE       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGM1')
							EWE_C_PRE      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGM1')
							CTRE_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CME1')
							AMPL_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_AME1')
							SGMA_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SME1')
							
							x_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XAM1')
							y_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YAM1')
							
							pre_shf_lim    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GSM1')
							pre_shf_ctr    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GCM1')
							
						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!' ,'yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CGM1','yellow')
							print colored ('Gotta Fit first before fixing a component (PRE)!','yellow')
							print colored ('Or UnFix (fix_pre_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print 'line 8598'
							print
							quit()
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-PRE','cyan')
					print
					print '******************************************************************************'
					print colored('Boundaries for Gaussian Fitting (PRE-CTR): ' + str(pre_shf_ctr),'cyan')
					print colored('Boundaries for Gaussian Fitting (PRE-LIM): ' + str(pre_shf_lim),'cyan')
					print '******************************************************************************'					
					print
					print colored(str(LINES[5][lines])+'_CGM1: ' + str(CTRE_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_AGM1: ' + str(AMPL_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_SGM1: ' + str(SGMA_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_FGM1: ' + str(FWHM_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_WGM1: ' + str(EW_C_PRE)      ,'cyan')
					print colored(str(LINES[5][lines])+'_EGM1: ' + str(EWE_C_PRE),'cyan')
					print
					print colored(str(LINES[5][lines])+'_XAM1 : ' + str(x_a),'cyan')
					print colored(str(LINES[5][lines])+'_YAM1 : ' + str(y_a),'cyan')
					print colored(str(LINES[5][lines])+'_GSM1 : ' + str(pre_shf_lim),'cyan')
					print colored(str(LINES[5][lines])+'_GCM1 : ' + str(pre_shf_ctr),'cyan')
					print
					print colored('Fit Values (PRE) Center, Amplitude, Sigma ('+fit_type+'):','cyan')
					print colored(str(CTRE_G_C_PRE)+', '+str(AMPL_G_C_PRE)+', '+str(SGMA_G_C_PRE),'cyan')
					print
					#####################################################PRE GAUSSIAN#################################################
					###################################################POST GAUSSIAN##################################################
					if fix_pst_gau == False and pst_shf_lim>0:
						###########################################DEFINING PRE-PST REGIONS##################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
							lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						else:
							print
							print colored('Using Expected Line Center to Define Fitting Region!','yellow')
							print
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
							
						print
						print colored('Region limits for fitting.','yellow')
						print colored('Line Center    : ' + str(LINES[0][lines]),'cyan')
						print colored('pst_shf_ctr    : ' + str(pst_shf_ctr),'cyan')
						print colored('pst_shf_lim    : ' + str(pst_shf_lim),'cyan')
						print colored('Line Center    : ' + str(LINES[0][lines] + pst_shf_ctr),'cyan')
						print colored('Lower Limit    : ' + str(LINES[0][lines] - (pst_shf_ctr-pst_shf_lim)),'cyan')
						print colored('Upper Limit    : ' + str(LINES[0][lines] - (pst_shf_ctr+pst_shf_lim)),'cyan')
						print 
						print colored('Central    : ' + str(LINES[0][lines]),'cyan')
						print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'cyan')
						print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'cyan')
						print
						print colored('Limits:','yellow')
						print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
						print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim
						print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
						print
						print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
						print
						###########################################DEFINING PRE-PST REGIONS##################################################						
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						x_b = lambda_glx[mask_ft_pst][X0_f2DG_indx_PST]
						y_b = inten_glx[mask_ft_pst][X0_f2DG_indx_PST]
						try:
							print colored('2-Fitting gaussian after line','cyan')
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
							A_f2DG_PST             = (max(inten_glx[mask_ft_pst])-1)
							X0_f2DG_PST            = LINES[0][lines] + pst_shf_ctr#
							SIGMA_f2DG_PST         = SIGMA_f2DG/2
							initial_guess_C_PST    = (X0_f2DG_PST,A_f2DG_PST,SIGMA_f2DG_PST)
							
							print
							print colored('Initial Guess Values PST : ','cyan')
							print colored(initial_guess_C_PST,'cyan')
							print
							
							gmodel_C_PST           = Model(func_1D_Gaussian_Emm)
							gmodel_C_PST.set_param_hint('X_0'  , value=X0_f2DG_PST , min=X0_f2DG_PST-(X0_f2DG_PST*LINES[7][lines]), max=X0_f2DG_PST+(X0_f2DG_PST*LINES[7][lines]))
							gmodel_C_PST.set_param_hint('A'    , value=A_f2DG_PST  , min=A_f2DG_PST -(A_f2DG_PST*LINES[10][lines]), max=A_f2DG_PST +(A_f2DG_PST*LINES[10][lines]))#min=A_f2DG_PST-0.001, max=A_f2DG_PST)
							gmodel_C_PST.set_param_hint('SIGMA', value=SIGMA_f2DG_PST)
							pars_C_PST             = gmodel_C_PST.make_params()
							result_C_PST           = gmodel_C_PST.fit(inten_glx[mask_ft_pst],pars_C_PST,
													X=lambda_glx[mask_ft_pst],X_0=X0_f2DG_PST,A=A_f2DG_PST,SIGMA=SIGMA_f2DG_PST,
													nan_policy = 'omit')
							CTRE_G_C_PST           = result_C_PST.params['X_0'].value
							AMPL_G_C_PST           = result_C_PST.params['A'].value
							SGMA_G_C_PST           = abs(result_C_PST.params['SIGMA'].value)
							FWHM_G_C_PST           = lw_sgma2fwhm(SGMA_G_C_PST)

							W_C_PST                = integrate.quad(lambda x: AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)), -np.inf, np.inf)
							EW_C_PST               = np.round(abs(np.asarray(W_C_PST[0])),10)
							EWE_C_PST              = np.round(abs(np.asarray(W_C_PST[1])),10)
							data_fitted_C_PST      = func_1D_Gaussian_Emm((lambda_glx[mask_ft_pst]), CTRE_G_C_PST,AMPL_G_C_PST,SGMA_G_C_PST)
							
							CTRE_G_C_PST_E         = result_C_PST.params['X_0'].stderr
							AMPL_G_C_PST_E         = result_C_PST.params['A'].stderr
							SGMA_G_C_PST_E         = result_C_PST.params['SIGMA'].stderr
							
							CTRE_G_C_PST_cor       = result_C_PST.params['X_0'].correl
							AMPL_G_C_PST_cor       = result_C_PST.params['A'].correl
							SGMA_G_C_PST_cor       = result_C_PST.params['SIGMA'].correl
							
							AMPL_SNR           = AMPL_G_C_PST
							CTRE_SNR           = CTRE_G_C_PST
							SGMA_SNR           = abs(SGMA_G_C_PST)
							
							if CTRE_G_C_PST_E == None:
								CTRE_G_C_PST_E = 999999.99999
							else:
								pass
							if AMPL_G_C_PST_E == None:
								AMPL_G_C_PST_E = 999999.99999
							else:
								pass
							if SGMA_G_C_PST_E == None:
								SGMA_G_C_PST_E = 999999.99999
							else:
								pass
							if CTRE_G_C_PST_cor == None:
								CTRE_G_C_PST_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_PST_cor == None:
								AMPL_G_C_PST_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_PST_cor == None:
								SGMA_G_C_PST_cor = 999999.99999
							else:
								pass
							chisqr_C_PST           = result_C_PST.chisqr
							redchi_C_PST           = result_C_PST.redchi
							
							W_C_PS2    = integrate.quad(lambda x: AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)), -np.inf, x_b)
							EW_C_PS2   = np.round(abs(np.asarray(W_C_PS2[0])),10)
							EWE_C_PS2  = np.round(abs(np.asarray(W_C_PS2[1])),10)
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C_PST, pcov_C_PST  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_PST          = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_PST        = 999999.99999
							AMPL_G_C_PST        = 999999.99999
							SGMA_G_C_PST        = 999999.99999
							FWHM_G_C_PST        = 999999.99999
							EW_C_PST            = 999999.99999
							EWE_C_PST           = 999999.99999
							
														
							CTRE_G_C_PST_E      = 999999.99999
							AMPL_G_C_PST_E      = 999999.99999
							SGMA_G_C_PST_E      = 999999.99999
							CTRE_G_C_PST_cor    = 999999.99999
							AMPL_G_C_PST_cor    = 999999.99999
							SGMA_G_C_PST_cor    = 999999.99999
							chisqr_C_PST        = 999999.99999
							redchi_C_PST        = 999999.99999
														
						if fit_vls_hdr == True and fix_pst_gau == False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGM2',float(CTRE_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGM2',float(AMPL_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGM2',float(SGMA_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGM2',float(FWHM_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM2',float(EW_C_PST)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM2',float(EWE_C_PST)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CME2',float(CTRE_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AME2',float(AMPL_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SME2',float(SGMA_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CME2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AME2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SME2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_XAM2',float(x_b),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PST GAU-LNR X2 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YAM2',float(y_b),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PST GAU-LNR Y2 COO')
						else:
							print
							print colored('The fit (PST) values will NOT be added to the fits headers!','magenta')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GSM2',float(pst_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for PST G-2')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GCM2',float(pst_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for PST G-2')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass
					elif fix_pst_gau == True or (fix_pst_gau == False and pst_shf_lim<=0):
						try:
							print
							print colored('2 PST-gaussian already fitted!','yellow')
							print colored('Values from fits header','yellow')
							CTRE_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM2')
							AMPL_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM2')
							SGMA_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM2')
							FWHM_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM2')
							EW_C_PST       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGM2')
							EWE_C_PST      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGM2')
							CTRE_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CME2')
							AMPL_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_AME2')
							SGMA_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SME2')
							
							x_b            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XAM2')
							y_b            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YAM2')
							
							pst_shf_lim    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GSM2')
							pst_shf_ctr    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GCM2')
						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
							print colored ('Gotta Fit first before fixing a component (PST)!','yellow')
							print colored ('Or UnFix (fix_pst_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							print 'line 8828'
							quit()
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-PST','magenta')
					print
					print '******************************************************************************'
					print colored('Boundaries for Gaussian Fitting (PST-CTR): ' + str(pst_shf_ctr),'magenta')					
					print colored('Boundaries for Gaussian Fitting (PST-LIM): ' + str(pst_shf_lim),'magenta')
					print '******************************************************************************'
					print
					print colored(str(LINES[5][lines])+'_CGM2: ' + str(CTRE_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_AGM2: ' + str(AMPL_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_SGM2: ' + str(SGMA_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_FGM2: ' + str(FWHM_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_WGM2: ' + str(EW_C_PST)      ,'magenta')
					print colored(str(LINES[5][lines])+'_EGM2: ' + str(EWE_C_PST),'magenta')
					print
					print colored(str(LINES[5][lines])+'_GSM2 : ' + str(pst_shf_lim),'magenta')
					print colored(str(LINES[5][lines])+'_GCM2 : ' + str(pst_shf_ctr),'magenta')
					print colored(str(LINES[5][lines])+'_XAM2 : ' + str(x_b),'magenta')
					print colored(str(LINES[5][lines])+'_YAM2 : ' + str(y_b),'magenta')
					print
					print colored('Fit Values (PST) Center, Amplitude, Sigma ('+fit_type+'):','magenta')
					print colored(str(CTRE_G_C_PST)+', '+str(AMPL_G_C_PST)+', '+str(SGMA_G_C_PST),'magenta')
					print					
					###################################################POST GAUSSIAN##################################################										

					###################################################MDL GAUSSIAN###################################################
					if fix_mdl_gau == False:# and mdl_shf_lim>0:
						if mdl_shf_lim>0:
							#########################################DEFINING PRE-PST-MDL REGIONS################################################
							if ofs_ctr_fit == True:
								print
								print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
								print
								lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
								lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
								mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
								mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
								mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
								mask_ft_mdl = (lambda_glx >= LINES[0][lines] + mdl_shf_ctr - mdl_shf_lim) & (lambda_glx <= LINES[0][lines] + mdl_shf_ctr + mdl_shf_lim)
								X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
								X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
								X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
							else:
								print
								print colored('Using Expected Line Center to Define Fitting Region!','yellow')
								print
								mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
								mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
								mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
								mask_ft_mdl = (lambda_glx >= LINES[0][lines] + mdl_shf_ctr - mdl_shf_lim) & (lambda_glx <= LINES[0][lines] + mdl_shf_ctr + mdl_shf_lim)
								X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
								X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
								X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
							print
							print colored('Region limits for fitting.','yellow')
							print colored('Line Center    : ' + str(LINES[0][lines]),'cyan')
							print colored('pst_shf_ctr    : ' + str(pst_shf_ctr),'cyan')
							print colored('pst_shf_lim    : ' + str(pst_shf_lim),'cyan')
							print colored('Line Center    : ' + str(LINES[0][lines] + pst_shf_ctr),'cyan')
							print colored('Lower Limit    : ' + str(LINES[0][lines] - (pst_shf_ctr-pst_shf_lim)),'cyan')
							print colored('Upper Limit    : ' + str(LINES[0][lines] - (pst_shf_ctr+pst_shf_lim)),'cyan')
							print 
							print colored('Central    : ' + str(LINES[0][lines]),'cyan')
							print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'cyan')
							print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'cyan')
							print colored('Central+MDL: ' + str(mdl_shf_ctr) + '-' + str(LINES[0][lines]+mdl_shf_ctr),'cyan')
							print
							print colored('Limits:','yellow')
							print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
							print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim
							print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
							print 'MDL    : lambda_glx >= ',LINES[0][lines] + mdl_shf_ctr-mdl_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + mdl_shf_ctr+mdl_shf_lim
							print
							print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
							print
							#########################################DEFINING PRE-PST-MDL REGIONS################################################
							X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
							x_c = lambda_glx[mask_ft_mdl][X0_f2DG_indx_MDL]
							y_c = inten_glx[mask_ft_mdl][X0_f2DG_indx_MDL]
							try:
								print colored('3-Fitting gaussian between lines','green')
								X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
								A_f2DG_MDL             = (max(inten_glx[mask_ft_mdl])-1)
								X0_f2DG_MDL            = LINES[0][lines] + mdl_shf_ctr
								SIGMA_f2DG_MDL         = SIGMA_f2DG/2
								initial_guess_C_MDL    = (X0_f2DG_MDL,A_f2DG_MDL,SIGMA_f2DG_MDL)
								
								print
								print colored('Initial Guess Values MDL : ','green')
								print colored(initial_guess_C_MDL,'green')
								print
								
								gmodel_C_MDL           = Model(func_1D_Gaussian_Emm)
								gmodel_C_MDL.set_param_hint('X_0'  , value=X0_f2DG_MDL , min=X0_f2DG_MDL-(X0_f2DG_MDL*LINES[7][lines]), max=X0_f2DG_MDL+(X0_f2DG_MDL*LINES[7][lines]))
								gmodel_C_MDL.set_param_hint('A'    , value=A_f2DG_MDL  , min=A_f2DG_MDL -(A_f2DG_MDL*LINES[10][lines]), max=A_f2DG_MDL +(A_f2DG_MDL*LINES[10][lines]))#min=A_f2DG_MDL-0.001, max=A_f2DG_MDL)
								gmodel_C_MDL.set_param_hint('SIGMA', value=SIGMA_f2DG_MDL)
								pars_C_MDL             = gmodel_C_MDL.make_params()
								result_C_MDL           = gmodel_C_MDL.fit(inten_glx[mask_ft_mdl],pars_C_MDL,
														X=lambda_glx[mask_ft_mdl],X_0=X0_f2DG_MDL,A=A_f2DG_MDL,SIGMA=SIGMA_f2DG_MDL,
														nan_policy = 'omit')
								CTRE_G_C_MDL           = result_C_MDL.params['X_0'].value
								AMPL_G_C_MDL           = result_C_MDL.params['A'].value
								SGMA_G_C_MDL           = abs(result_C_MDL.params['SIGMA'].value)
								FWHM_G_C_MDL           = lw_sgma2fwhm(SGMA_G_C_MDL)

								W_C_MDL                = integrate.quad(lambda x: AMPL_G_C_MDL*np.exp(-((x)**2)/(2*SGMA_G_C_MDL**2)), -np.inf, np.inf)
								EW_C_MDL               = np.round(abs(np.asarray(W_C_MDL[0])),10)
								EWE_C_MDL              = np.round(abs(np.asarray(W_C_MDL[1])),10)
								data_fitted_C_MDL      = func_1D_Gaussian_Emm((lambda_glx[mask_ft_mdl]), CTRE_G_C_MDL,AMPL_G_C_MDL,SGMA_G_C_MDL)
								
								CTRE_G_C_MDL_E         = result_C_MDL.params['X_0'].stderr
								AMPL_G_C_MDL_E         = result_C_MDL.params['A'].stderr
								SGMA_G_C_MDL_E         = result_C_MDL.params['SIGMA'].stderr
								
								CTRE_G_C_MDL_cor       = result_C_MDL.params['X_0'].correl
								AMPL_G_C_MDL_cor       = result_C_MDL.params['A'].correl
								SGMA_G_C_MDL_cor       = result_C_MDL.params['SIGMA'].correl
								
								AMPL_SNR               = AMPL_G_C_MDL
								CTRE_SNR               = CTRE_G_C_MDL
								SGMA_SNR               = abs(SGMA_G_C_MDL)
								
								if CTRE_G_C_MDL_E == None:
									CTRE_G_C_MDL_E = 999999.99999
								else:
									pass
								if AMPL_G_C_MDL_E == None:
									AMPL_G_C_MDL_E = 999999.99999
								else:
									pass
								if SGMA_G_C_MDL_E == None:
									SGMA_G_C_MDL_E = 999999.99999
								else:
									pass
								if CTRE_G_C_MDL_cor == None:
									CTRE_G_C_MDL_cor = 999999.99999
								else:
									pass
								if AMPL_G_C_MDL_cor == None:
									AMPL_G_C_MDL_cor = 999999.99999
								else:
									pass
								if SGMA_G_C_MDL_cor == None:
									SGMA_G_C_MDL_cor = 999999.99999
								else:
									pass
								chisqr_C_MDL           = result_C_MDL.chisqr
								redchi_C_MDL           = result_C_MDL.redchi
								
								W_C_PS3    = integrate.quad(lambda x: AMPL_G_C_MDL*np.exp(-((x)**2)/(2*SGMA_G_C_MDL**2)), -np.inf, x_c)
								EW_C_PS3   = np.round(abs(np.asarray(W_C_PS3[0])),10)
								EWE_C_PS3  = np.round(abs(np.asarray(W_C_PS3[1])),10)
							except (RuntimeError,ValueError,TypeError):
								print colored('RuntimeError','green')
								popt_C_MDL, pcov_C_MDL  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
								perr_C_MDL          = [999999.99999,999999.99999,999999.99999]
								CTRE_G_C_MDL        = 999999.99999
								AMPL_G_C_MDL        = 999999.99999
								SGMA_G_C_MDL        = 999999.99999
								FWHM_G_C_MDL        = 999999.99999
								EW_C_MDL            = 999999.99999
								EWE_C_MDL           = 999999.99999
								
															
								CTRE_G_C_MDL_E      = 999999.99999
								AMPL_G_C_MDL_E      = 999999.99999
								SGMA_G_C_MDL_E      = 999999.99999
								CTRE_G_C_MDL_cor    = 999999.99999
								AMPL_G_C_MDL_cor    = 999999.99999
								SGMA_G_C_MDL_cor    = 999999.99999
								chisqr_C_MDL        = 999999.99999
								redchi_C_MDL        = 999999.99999
						else:
							popt_C_MDL, pcov_C_MDL  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_MDL          = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_MDL        = 999999.99999
							AMPL_G_C_MDL        = 999999.99999
							SGMA_G_C_MDL        = 999999.99999
							FWHM_G_C_MDL        = 999999.99999
							EW_C_MDL            = 999999.99999
							EWE_C_MDL           = 999999.99999
							
														
							CTRE_G_C_MDL_E      = 999999.99999
							AMPL_G_C_MDL_E      = 999999.99999
							SGMA_G_C_MDL_E      = 999999.99999
							CTRE_G_C_MDL_cor    = 999999.99999
							AMPL_G_C_MDL_cor    = 999999.99999
							SGMA_G_C_MDL_cor    = 999999.99999
							chisqr_C_MDL        = 999999.99999
							redchi_C_MDL        = 999999.99999
							
							x_c                 = 999999.99999
							y_c                 = 999999.99999							
						if fit_vls_hdr == True and fix_mdl_gau == False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGM3',float(CTRE_G_C_MDL)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGM3',float(AMPL_G_C_MDL)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGM3',float(SGMA_G_C_MDL)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGM3',float(FWHM_G_C_MDL)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM3',float(EW_C_MDL)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM3',float(EWE_C_MDL)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CME2',float(CTRE_G_C_MDL_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AME2',float(AMPL_G_C_MDL_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SME2',float(SGMA_G_C_MDL_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CME2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AME2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SME2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_XAM3',float(x_c),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' MDL GAU-LNR X3 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YAM3',float(y_c),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' MDL GAU-LNR Y3 COO')
						else:
							print
							print colored('The fit (MDL) values will NOT be added to the fits headers!','green')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GSM3',float(mdl_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for MDL G-3')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GCM3',float(mdl_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for MDL G-3')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass
					elif fix_mdl_gau == True:
						try:
							print
							print colored('3- MDL-gaussian already fitted!','yellow')
							print colored('Values from fits header','yellow')
							CTRE_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM3')
							AMPL_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM3')
							SGMA_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM3')
							FWHM_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM3')
							EW_C_MDL       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGM3')
							EWE_C_MDL      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGM3')
							CTRE_G_C_MDL_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CME2')
							AMPL_G_C_MDL_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_AME2')
							SGMA_G_C_MDL_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SME2')
							
							x_c            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XAM3')
							y_c            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YAM3')
							EW_C_PS3       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WRM3')
							EWE_C_PS3      = Header_Get(specfile_glx,str(LINES[5][lines])+'_ERM3')
							
							mdl_shf_lim    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GSM3')
							mdl_shf_ctr    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GCM3')
						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
							print colored ('Gotta Fit first before fixing a component (PST)!','yellow')
							print colored ('Or UnFix (fix_pst_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print 'line 9088'
							print
							quit()
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-MDL','green')
					print
					print '******************************************************************************'
					print colored('Boundaries for Gaussian Fitting (PST-CTR): ' + str(mdl_shf_ctr),'green')					
					print colored('Boundaries for Gaussian Fitting (PST-LIM): ' + str(mdl_shf_lim),'green')
					print '******************************************************************************'
					print
					print colored(str(LINES[5][lines])+'_CGM3: ' + str(CTRE_G_C_MDL)  ,'green')
					print colored(str(LINES[5][lines])+'_AGM3: ' + str(AMPL_G_C_MDL)  ,'green')
					print colored(str(LINES[5][lines])+'_SGM3: ' + str(SGMA_G_C_MDL)  ,'green')
					print colored(str(LINES[5][lines])+'_FGM3: ' + str(FWHM_G_C_MDL)  ,'green')
					print colored(str(LINES[5][lines])+'_WGM3: ' + str(EW_C_MDL)      ,'green')
					print colored(str(LINES[5][lines])+'_EGM3: ' + str(EWE_C_MDL),'green')
					print
					print colored(str(LINES[5][lines])+'_GSM3 : ' + str(mdl_shf_lim),'green')
					print colored(str(LINES[5][lines])+'_GCM3 : ' + str(mdl_shf_ctr),'green')
					print colored(str(LINES[5][lines])+'_XAM3 : ' + str(x_c),'green')
					print colored(str(LINES[5][lines])+'_YAM3 : ' + str(y_c),'green')
					print
					print colored('Fit Values (PST) Center, Amplitude, Sigma ('+fit_type+'):','green')
					print colored(str(CTRE_G_C_MDL)+', '+str(AMPL_G_C_MDL)+', '+str(SGMA_G_C_MDL),'green')
					print					
					###################################################MDL GAUSSIAN##################################################						
					###############################################COMPUTING TOTAL AREA###############################################
					print colored('Computing Flux Area','yellow')
					#############################################COMPUTING LINEAR AREA###################################################
					slope_line1 = (y_a-y_b)/(x_a-x_b)
					slope_line2 = (y_b-y_a)/(x_b-x_a)
					b1 = y_a - (slope_line1*x_a)
					b2 = y_b - (slope_line1*x_b)
					print
					print  colored('Computing Linear Area considering peak points:','yellow')
					print 'Point A: ',x_a,y_a
					print 'Point B: ',x_b,y_b
					print 'Slope: ',slope_line1
					print 'Slope: ',slope_line2
					print 'b: ',b1,b2
					print
					#############################################COMPUTING LINEAR AREA###################################################
					###############################################COMPUTING TOTAL AREA###############################################
					CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC1')
					AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC1')
					SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC1')

					CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC2')
					AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC2')
					SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC2')

					CTRE_G_C_PRE  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM1')
					AMPL_G_C_PRE  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM1')
					SGMA_G_C_PRE  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM1')

					CTRE_G_C_PST  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM2')
					AMPL_G_C_PST  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM2')
					SGMA_G_C_PST  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM2')

					CTRE_G_C_MDL  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM3')
					AMPL_G_C_MDL  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM3')
					SGMA_G_C_MDL  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM3')

					print
					print colored('Computing Areas using info from fits headers.','yellow')
					print 'CTRE_G_C_1: ',str(CTRE_G_C_1),' AMPL_G_C-1: ',str(AMPL_G_C_1),' SGMA_G_C-1: ',str(SGMA_G_C_1)
					print 'CTRE_G_C_2: ',str(CTRE_G_C_2),' AMPL_G_C-2: ',str(AMPL_G_C_2),' SGMA_G_C-2: ',str(SGMA_G_C_2)
					print 'CTRE_G_C_PRE: ',str(CTRE_G_C_PRE),' AMPL_G_C_PRE: ',str(AMPL_G_C_PRE),' SGMA_G_C_PRE: ',str(SGMA_G_C_PRE)
					print 'CTRE_G_C_PST: ',str(CTRE_G_C_PST),' AMPL_G_C_PST: ',str(AMPL_G_C_PST),' SGMA_G_C_PST: ',str(SGMA_G_C_PST)
					print 'CTRE_G_C_MDL: ',str(CTRE_G_C_MDL),' AMPL_G_C_MDL: ',str(AMPL_G_C_MDL),' SGMA_G_C_MDL: ',str(SGMA_G_C_PST)
					print

					W_C_1     = integrate.quad(lambda x: AMPL_G_C_1*np.exp(-((x)**2)/(2*SGMA_G_C_1**2)), -np.inf, np.inf)
					EW_C_1    = np.round(abs(np.asarray(W_C_1[0])),10)
					EWE_C_1   = np.round(abs(np.asarray(W_C_1[1])),10)

					W_C_2     = integrate.quad(lambda x: AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2)), -np.inf, np.inf)
					EW_C_2    = np.round(abs(np.asarray(W_C_2[0])),10)
					EWE_C_2   = np.round(abs(np.asarray(W_C_2[1])),10)

					W_C_PRE   = integrate.quad(lambda x: AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2)), -np.inf, np.inf)
					EW_C_PRE  = np.round(abs(np.asarray(W_C_PRE[0])),10)
					EWE_C_PRE = np.round(abs(np.asarray(W_C_PRE[1])),10)


					W_C_PST   = integrate.quad(lambda x: AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)), -np.inf, np.inf)
					EW_C_PST  = np.round(abs(np.asarray(W_C_PST[0])),10)
					EWE_C_PST = np.round(abs(np.asarray(W_C_PST[1])),10)

					W_C_MDL   = integrate.quad(lambda x: AMPL_G_C_MDL*np.exp(-((x)**2)/(2*SGMA_G_C_MDL**2)), -np.inf, np.inf)
					EW_C_MDL  = np.round(abs(np.asarray(W_C_MDL[0])),10)
					EWE_C_MDL = np.round(abs(np.asarray(W_C_MDL[1])),10)

					W_PLP     = integrate.quad(lambda x: AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2))    + 
												AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2))+ 
												AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2))+
												AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)),
												-np.inf, np.inf
												)
					EW_PLP    = np.round((np.asarray(W_PLP[0])),10)
					EWE_PLP   = np.round(abs(np.asarray(W_PLP[1])),10)

					if EW_C_1 == 999999.99999:
						EW_C_1 = 0
					else:
						pass
					if EW_C_2 == 999999.99999:
						EW_C_2 = 0
					else:
						pass
					if EW_C_PRE == 999999.99999:
						EW_C_PRE = 0
					else:
						pass
					if EW_C_PST == 999999.99999:
						EW_C_PST = 0
					else:
						pass						
					EWMT    = EW_C_1 + EW_C_2 + EW_C_PRE + EW_C_PST #+ EW_C_LNR

					print
					print colored('Areas     :','yellow')
					print colored('Area CTR-1: ' + str(EW_C_1),'yellow')
					print colored('Area CTR-2: ' + str(EW_C_2),'yellow')					
					print colored('Area PRE-G: ' + str(EW_C_PRE),'blue')
					print colored('Area PST-G: ' + str(EW_C_PST),'magenta')
					print colored('Area MDL-G: ' + str(EW_C_MDL),'green')
					print colored('Area PRE-CTR1-CTR2-PST: ' + str(EW_PLP),'yellow')
					print					
					###############################################COMPUTING TOTAL AREA###############################################
					#############################################ADDING AREA TO FTIS HEADER#############################################
					if fit_vls_hdr == True:
						print
						print colored('The Areas values will be updated to the fits headers!','magenta')
						print
						print colored('Area CTR-1: '                  + str(EW_C_1)   + '-' +str(LINES[5][lines])+'_WMC1','blue')
						print colored('Area CTR-2: '                  + str(EW_C_2)   + '-' +str(LINES[5][lines])+'_WMC2','red')
						print colored('Area PRE-G: '                  + str(EW_C_PRE) + '-' +str(LINES[5][lines])+'_WGM1','cyan')
						print colored('Area PST-G: '                  + str(EW_C_PST) + '-' +str(LINES[5][lines])+'_WGM2','magenta')
						print colored('Area MDL-G: '                  + str(EW_C_MDL) + '-' +str(LINES[5][lines])+'_WGM3','green')
						print colored('Area PLP-G: '                  + str(EW_PLP)   + '-' +str(LINES[5][lines])+'_WMPP','yellow')
						print colored('Area TOT=PRE-CTR1-CTR2-PST: '  + str(EWMT)     + '-' +str(LINES[5][lines])+'_WMPT','yellow')
						print						
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC1',float(EW_C_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C1-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC1',float(EWE_C_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C1-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC2',float(EW_C_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C2-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC2',float(EWE_C_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C2-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM1',float(EW_C_PRE)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   PRE Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM1',float(EWE_C_PRE)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  PRE Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM2',float(EW_C_PST)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   PST Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM2',float(EWE_C_PST)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  PST Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM3',float(EW_C_MDL)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   MDL Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM3',float(EWE_C_MDL)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  MDL Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMPP',float(EW_C_PST)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   PRE-C1-C2-PST Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMPP',float(EWE_C_PST)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  PRE-C1-C2-PST Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMPT',float(EWE_C_PST)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' WE   PLP Crct' + str(fit_type))

					else:
						print
						print colored('The Areas values will NOT be updated to the fits headers!','magenta')
						print
					#############################################ADDING AREA TO FTIS HEADER#############################################
					#################################DEFINE AREA FOR LATER ON SUBSTRACT OFFSET FOR PLOTING###############################
					###########################################DEFINING PRE-PST REGIONS##################################################
					if ofs_ctr_fit == True:
						print
						print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
						print
						lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
						lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
						mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
						mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_mdl = (lambda_glx >= LINES[0][lines] + mdl_shf_ctr - mdl_shf_lim) & (lambda_glx <= LINES[0][lines] + mdl_shf_ctr + mdl_shf_lim)
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
					else:
						print
						print colored('Using Expected Line Center to Define Fitting Region!','yellow')
						print
						mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
						mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_mdl = (lambda_glx >= LINES[0][lines] + mdl_shf_ctr - mdl_shf_lim) & (lambda_glx <= LINES[0][lines] + mdl_shf_ctr + mdl_shf_lim)
						try:
							X0_f2DG_indx_PRE = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						except ValueError:
							X0_f2DG_indx_PRE = LINES[0][lines]

						try:
							X0_f2DG_indx_PST = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						except ValueError:
							X0_f2DG_indx_PST = LINES[0][lines]

						try:
							X0_f2DG_indx_MDL = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
						except ValueError:
							X0_f2DG_indx_MDL = LINES[0][lines]

						
					print
					print colored('Region limits for fitting.','yellow')
					print colored('Central    : ' + str(LINES[0][lines]),'yellow')
					print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'cyan')
					print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'magenta')
					print colored('Central+MDL: ' + str(mdl_shf_ctr) + '-' + str(LINES[0][lines]+mdl_shf_ctr),'green')
					print
					print colored('Limits:','yellow')
					print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
					print 'PRE    : lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim,'-','lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim
					print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
					print 'MDL    : lambda_glx >= ',LINES[0][lines] + mdl_shf_ctr-mdl_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + mdl_shf_ctr+mdl_shf_lim
					print
					print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
					print
					############################################DEFINING PRE-PST REGIONS##################################################
					##################################DEFINE AREA FOR LATER ON SUBSTRACT OFFSET FOR PLOTING###############################
				elif 'Dbl' not in LINES[3][lines] and fit_fnct=='gauss'   and fit_type == 'lmfit' and mke_lne_fit == True and uft_lne_vls == False:
					fit_typ = 'G'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					from lmfit import Model

					lmb_min_lim_line_ft = (LINES[0][lines]-LINES[8][lines]) - MSK_NTMS*LINES[1][lines]  
					lmb_max_lim_line_ft = (LINES[0][lines]+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
					lmb_min_lim_line    = (LINES[0][lines]-LINES[8][lines])*(1+z_glx_Ps) - 1.5*MSK_NTMS*LINES[2][lines] 
					lmb_max_lim_line    = (LINES[0][lines]+LINES[8][lines])*(1+z_glx_Ps) + 1.5*MSK_NTMS*LINES[2][lines]

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function

					X0_f2DG    = (LINES[0][lines]+LINES[8][lines])
					SIGMA_f2DG = LINES[1][lines]
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
					try:
						gmodel_0           = Model(func_1D_Gaussian)
						gmodel_0.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG-(X0_f2DG*LINES[7][lines]), max=X0_f2DG+(X0_f2DG*LINES[7][lines]))
						gmodel_0.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*LINES[10][lines]) , max=A_f2DG  + (A_f2DG*LINES[10][lines]))
						gmodel_0.set_param_hint('SIGMA', value=SIGMA_f2DG)
						pars_0             = gmodel_0.make_params()							
						result_0           = gmodel_0.fit(inten_glx[mask_ft],pars_0,
												X=lambda_glx[mask_ft],
												X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
												nan_policy = 'omit')
						CTRE_G_0           = result_0.params['X_0'].value
						AMPL_G_0           = result_0.params['A'].value
						SGMA_G_0           = abs(result_0.params['SIGMA'].value)
						FWHM_G_0           = lw_sgma2fwhm(SGMA_G_0)
						W_0                = integrate.quad(lambda x: AMPL_G_0*np.exp(-((x)**2)/(2*SGMA_G_0**2)), -np.inf, np.inf)
						EW_0               = np.round(abs(np.asarray(W_0[0])),3)
						EWE_0              = np.round(abs(np.asarray(W_0[1])),10)
						data_fitted_0      = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_0,AMPL_G_0,SGMA_G_0)

						CTRE_G_0_E         = result_0.params['X_0'].stderr
						AMPL_G_0_E         = result_0.params['A'].stderr
						SGMA_G_0_E         = result_0.params['SIGMA'].stderr

						CTRE_G_0_cor       = result_0.params['X_0'].correl
						AMPL_G_0_cor       = result_0.params['A'].correl
						SGMA_G_0_cor       = result_0.params['SIGMA'].correl

						chisqr_0           = result_0.chisqr
						redchi_0           = result_0.redchi
					except (RuntimeError,ValueError,TypeError):
						popt_0, pcov_0     = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_0             = [999999.99999,999999.99999,999999.99999]
						CTRE_G_0           = 999999.99999
						AMPL_G_0           = 999999.99999
						SGMA_G_0           = 999999.99999
						FWHM_G_0           = 999999.99999
						EW_0               = 999999.99999
						EWE_0              = 999999.99999

						CTRE_G_0_E         = 999999.99999
						AMPL_G_0_E         = 999999.99999
						SGMA_G_0_E         = 999999.99999

						CTRE_G_0_cor       = 999999.99999
						AMPL_G_0_cor       = 999999.99999
						SGMA_G_0_cor       = 999999.99999

						chisqr_0           = 999999.99999
						redchi_0           = 999999.99999
					if fit_vls_hdr == True:
						Header_Add(specfile_glx,str(LINES[5][lines])+'_CGF0',float(CTRE_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_AGF0',float(AMPL_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_SGF0',float(SGMA_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'SGMA 1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_FGF0',float(FWHM_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGF0',float(EW_0)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGF0',float(EWE_0)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct')
					else:
						print
						print colored('The fit values will NOT be added to the fits headers!','magenta')
						print
					try:
						gmodel_O           = Model(func_1D_Gaussian_O)
						gmodel_O.set_param_hint('X_0'   , value=X0_f2DG   , min=X0_f2DG-(X0_f2DG*LINES[7][lines]), max=X0_f2DG+(X0_f2DG*LINES[7][lines]))
						gmodel_O.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*LINES[10][lines]) , max=A_f2DG  + (A_f2DG*LINES[10][lines]))
						gmodel_O.set_param_hint('SIGMA' , value=SIGMA_f2DG)
						gmodel_O.set_param_hint('OFFSET', value=max(inten_glx[mask_ft])-1)
						pars_O             = gmodel_O.make_params()

						result_O           = gmodel_O.fit(inten_glx[mask_ft],pars_O,
												X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,OFFSET=max(inten_glx[mask_ft])-1,
												nan_policy = 'omit')
						CTRE_G_O           = result_O.params['X_0'].value
						AMPL_G_O           = result_O.params['A'].value
						SGMA_G_O           = abs(result_O.params['SIGMA'].value)
						OFST_G_O           = abs(result_O.params['OFFSET'].value)
						FWHM_G_O           = lw_sgma2fwhm(SGMA_G_O)
						W_O                = integrate.quad(lambda x: AMPL_G_O*np.exp(-((x)**2)/(2*SGMA_G_O**2)), -np.inf, np.inf)
						EW_O               = np.round(abs(np.asarray(W_O[0])),3)
						EWE_O              = np.round(abs(np.asarray(W_O[1])),10)
						data_fitted_O      = func_1D_Gaussian_O((lambda_glx[mask_ft]),CTRE_G_O,AMPL_G_O,SGMA_G_O,OFST_G_O)

						CTRE_G_O_E         = result_O.params['X_0'].stderr
						AMPL_G_O_E         = result_O.params['A'].stderr
						SGMA_G_O_E         = result_O.params['SIGMA'].stderr

						CTRE_G_O_cor       = result_O.params['X_0'].correl
						AMPL_G_O_cor       = result_O.params['A'].correl
						SGMA_G_O_cor       = result_O.params['SIGMA'].correl

						chisqr_O           = result_O.chisqr
						redchi_O           = result_O.redchi
						
						#####################################################################################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_O-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
							lmb_max_lim_line_ft = (CTRE_G_O+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							#lmb_min_lim_line    = (CTRE_G_O-LINES[8][lines])*(1+z_glx_Ps) - 1.5*MSK_NTMS*LINES[2][lines]#-20#LINES[2][lines] - 10
							#lmb_max_lim_line    = (CTRE_G_O+LINES[8][lines])*(1+z_glx_Ps) + 1.5*MSK_NTMS*LINES[2][lines]#+20#LINES[2][lines] + 10
							#mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
							mask_ft     = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
						else:
							pass
						#####################################################################################################################

						initial_guess_C    = (X0_f2DG,A_f2DG,SIGMA_f2DG)

						gmodel_C           = Model(func_1D_Gaussian)
						gmodel_C.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG-(X0_f2DG*LINES[7][lines]), max=X0_f2DG+(X0_f2DG*LINES[7][lines]))
						gmodel_C.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*LINES[10][lines]) , max=A_f2DG  + (A_f2DG*LINES[10][lines]))
						gmodel_C.set_param_hint('SIGMA', value=SIGMA_f2DG)
						pars_C             = gmodel_C.make_params()
						result_C           = gmodel_C.fit(inten_glx[mask_ft],pars_C,
												X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
												nan_policy = 'omit')
						CTRE_G_C           = result_C.params['X_0'].value
						AMPL_G_C           = result_C.params['A'].value
						SGMA_G_C           = abs(result_C.params['SIGMA'].value)
						FWHM_G_C           = lw_sgma2fwhm(SGMA_G_C)
						W_C                = integrate.quad(lambda x: AMPL_G_C*np.exp(-((x)**2)/(2*SGMA_G_C**2)), -np.inf, np.inf)

						EW_C               = np.round(abs(np.asarray(W_C[0])),3)
						EWE_C              = np.round(abs(np.asarray(W_C[1])),10)
						data_fitted_C      = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_C,AMPL_G_C,SGMA_G_C)

						CTRE_G_C_E         = result_C.params['X_0'].stderr
						AMPL_G_C_E         = result_C.params['A'].stderr
						SGMA_G_C_E         = result_C.params['SIGMA'].stderr

						CTRE_G_C_cor       = result_C.params['X_0'].correl
						AMPL_G_C_cor       = result_C.params['A'].correl
						SGMA_G_C_cor       = result_C.params['SIGMA'].correl

						AMPL_SNR           = AMPL_G_C
						CTRE_SNR           = CTRE_G_C
						SGMA_SNR           = abs(SGMA_G_C)

						if CTRE_G_C_E == None:
							CTRE_G_C_E = 999999.99999
						else:
							pass
						if AMPL_G_C_E == None:
							AMPL_G_C_E = 999999.99999
						else:
							pass
						if SGMA_G_C_E == None:
							SGMA_G_C_E = 999999.99999
						else:
							pass
						if CTRE_G_C_cor == None:
							CTRE_G_C_cor = 999999.99999
						else:
							pass
						if AMPL_G_C_cor == None:
							AMPL_G_C_cor = 999999.99999
						else:
							pass
						if SGMA_G_C_cor == None:
							SGMA_G_C_cor = 999999.99999
						else:
							pass
						chisqr_C           = result_C.chisqr
						redchi_C           = result_C.redchi
					except (RuntimeError,ValueError,TypeError):
						print colored('RuntimeError','cyan')
						popt_C, pcov_C  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C          = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C        = 999999.99999
						AMPL_G_C        = 999999.99999
						SGMA_G_C        = 999999.99999
						FWHM_G_C        = 999999.99999
						EW_C            = 999999.99999
						EWE_C           = 999999.99999

						CTRE_G_C_E      = 999999.99999
						AMPL_G_C_E      = 999999.99999
						SGMA_G_C_E      = 999999.99999
						CTRE_G_C_cor    = 999999.99999
						AMPL_G_C_cor    = 999999.99999
						SGMA_G_C_cor    = 999999.99999
						chisqr_C        = 999999.99999
						redchi_C        = 999999.99999

						popt_O, pcov_O  = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
						perr_O          = [999999.99999,999999.99999,999999.99999,999999.99999]
						CTRE_G_O        = 999999.99999
						AMPL_G_O        = 999999.99999
						SGMA_G_O        = 999999.99999
						OFST_G_O        = 999999.99999
						FWHM_G_O        = 999999.99999
						EW_O            = 999999.99999
						EWE_O           = 999999.99999

						CTRE_G_O_E      = 999999.99999
						AMPL_G_O_E      = 999999.99999
						SGMA_G_O_E      = 999999.99999
						CTRE_G_O_cor    = 999999.99999
						AMPL_G_O_cor    = 999999.99999
						SGMA_G_O_cor    = 999999.99999
						OFST_G_O_cor    = 999999.99999
						chisqr_O        = 999999.99999
						redchi_O        = 999999.99999

						AMPL_SNR        = 999999.99999
						CTRE_SNR        = 999999.99999
						SGMA_SNR        = 999999.99999

					print
					print colored(specfile_glx,'cyan')
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines]),'yellow')
					print
					print colored(str(LINES[5][lines])+'_CGLC: ' + str(CTRE_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_AGLC: ' + str(AMPL_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_SGLC: ' + str(SGMA_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_FGLC: ' + str(FWHM_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_WGLC: ' + str(EW_C)      ,'yellow')
					print colored(str(LINES[5][lines])+'_EGLC: ' + str(EWE_C),'yellow')
					print
					print
					print
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','yellow')
					print colored(str(CTRE_G_C)+', '+str(AMPL_G_C)+', '+str(SGMA_G_C),'yellow')
					print
					print
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','yellow')
					print colored(str(CTRE_G_C)+', '+str(AMPL_G_C)+', '+str(SGMA_G_C),'yellow')
					print

					if fit_vls_hdr == True:
						Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                       ,header_comment = 'Fit function')
						Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                       ,header_comment = 'Fit method')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_CGLO',float(CTRE_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_AGLO',float(AMPL_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_SGLO',float(SGMA_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_FGLO',float(FWHM_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGLO',float(EW_O)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGLO',float(EWE_O)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_OFSO',float(OFST_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 1GF Offst' + str(fit_type))

						Header_Add(specfile_glx,str(LINES[5][lines])+'_CGLC',float(CTRE_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_AGLC',float(AMPL_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_SGLC',float(SGMA_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_FGLC',float(FWHM_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGLC',float(EW_C)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGLC',float(EWE_C)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))

						try:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLEC',float(CTRE_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALEC',float(AMPL_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLEC',float(SGMA_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
						except ValueError:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLEC',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALEC',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLEC',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))


						Header_Add(specfile_glx,str(LINES[5][lines])+'_CHGL',float(chisqr_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_CRGL',float(redchi_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 Reduced 1GF' + str(fit_type))
					else:
						print
						print colored('The fit values will NOT be added to the fits headers!','magenta')
						print
						pass

					###############################################COMPUTING TOTAL AREA###############################################
					Header_Get(specfile_glx,str(LINES[5][lines])+'_AGLC')
					Header_Get(specfile_glx,str(LINES[5][lines])+'_SGLC')

					print colored('Computing Flux Area','yellow')
					print
					print colored('Computing Area using info from fits headers.','yellow')
					print 'AMPL_G_C_CRC: ',str(AMPL_G_C),'SGMA_G_CRC: ',str(SGMA_G_C)
					print

					W_C   = integrate.quad(lambda x: AMPL_G_C*np.exp(-((x)**2)/(2*SGMA_G_C**2)), -np.inf, np.inf)
					EW_C  = np.round(abs(np.asarray(W_C[0])),10)
					EWE_C = np.round(abs(np.asarray(W_C[1])),10)

					print
					print colored('Areas     :','yellow')
					print colored('Area CRC-G: ' + str(EW_C),'yellow')
					print
					################################################COMPUTING TOTAL AREA################################################
					#############################################ADDING AREA TO FTIS HEADER#############################################
					if fit_vls_hdr == True:
						print
						print colored('The Areas values will be updated to the fits headers!','magenta')
						print
						print colored('Area CRC-G: ' + str(EW_C)   + '-' +str(LINES[5][lines])+'_WGLC','yellow')
						print

						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGLC',float(EW_C)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGLC',float(EWE_C)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
						print
						print colored('The LINEAR AND TOTAL Areas values will NOT be updated to the fits headers!','magenta')
						print
					else:
						pass
					#############################################ADDING AREA TO FTIS HEADER#############################################
				elif 'Dbl' not in LINES[3][lines] and fit_fnct=='gaussM'  and fit_type == 'lmfit' and mke_lne_fit == True and uft_lne_vls == False:
					fit_typ = 'GM'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					from lmfit import Model

	
					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
					if ivl_fts_hdr == True:
						try:
							L1_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF_0')        #LINES-1 Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP_0')        #LINES-2 Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF_0')        #LINES-7 Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print
							print colored('Initial fit variables from fits header!','yellow')
							print colored('Headers:','yellow')
							print
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_0),'yellow')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_0),'yellow')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_0),'yellow')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_0),'yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_0),'yellow')
							print colored('*****Success!******','magenta')
							print
						except ValueError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
							print '*****'
							print
							quit()
						try:
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
							print colored('*****Success!******','magenta')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF_0',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('*****Success!******','magenta')
							print
						################FOR CASES WHERE KLINE WIIDHT ==0################
						if L1_0  == 0:
							L1_0  = 1#LINES[1][lines]
						else:
							pass
						################FOR CASES WHERE KLINE WIIDHT ==0################
					elif ivl_fts_hdr == False:
						L1_0  = LINES[1][lines]
						L2_0  = LINES[2][lines]
						L7_0  = LINES[7][lines]
						L8_0  = LINES[8][lines]
						L10_0 = LINES[10][lines]
						print
						print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
						print					
					print str(LINES[5][lines])+'_WF_0' + ': ' + str(L1_0)
					print str(LINES[5][lines])+'_WP_0' + ': ' + str(L2_0)
					print str(LINES[5][lines])+'_CF_0' + ': ' + str(L7_0)
					print str(LINES[5][lines])+'_CO_0' + ': ' + str(L8_0)
					print str(LINES[5][lines])+'_AF_0' + ': ' + str(L10_0)
					#quit()

					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################

					lmb_min_lim_line_ft = (LINES[0][lines]-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
					lmb_max_lim_line_ft = (LINES[0][lines]+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
					lmb_min_lim_line    = (LINES[0][lines]-LINES[8][lines])*(1+z_glx_Ps) - 1.5*MSK_NTMS*LINES[2][lines]
					lmb_max_lim_line    = (LINES[0][lines]+LINES[8][lines])*(1+z_glx_Ps) + 1.5*MSK_NTMS*LINES[2][lines]

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function

					X0_f2DG    = (LINES[0][lines]+LINES[8][lines])
					SIGMA_f2DG = LINES[1][lines]
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
				
					##################################################CENTRAL GAUSSIAN###################################################
					if fix_ctr_gau == False:
						print
						print colored('0-Fitting Central line','cyan')
						print
						#################################################CENTRAL GAUSSIAN-0##################################################
						try:
							gmodel_0           = Model(func_1D_Gaussian)
							gmodel_0.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG-(X0_f2DG*LINES[7][lines]), max=X0_f2DG+(X0_f2DG*LINES[7][lines]))
							gmodel_0.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*LINES[10][lines]) , max=A_f2DG  + (A_f2DG*LINES[10][lines]))
							gmodel_0.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_0             = gmodel_0.make_params()							
							result_0           = gmodel_0.fit(inten_glx[mask_ft],pars_0,
													X=lambda_glx[mask_ft],
													X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_0           = result_0.params['X_0'].value
							AMPL_G_0           = result_0.params['A'].value
							SGMA_G_0           = abs(result_0.params['SIGMA'].value)
							FWHM_G_0           = lw_sgma2fwhm(SGMA_G_0)
							W_0                = integrate.quad(lambda x: AMPL_G_0*np.exp(-((x)**2)/(2*SGMA_G_0**2)), -np.inf, np.inf)
							EW_0               = np.round(abs(np.asarray(W_0[0])),10)
							EWE_0              = np.round(abs(np.asarray(W_0[1])),10)
							data_fitted_0      = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_0,AMPL_G_0,SGMA_G_0)

							CTRE_G_0_E         = result_0.params['X_0'].stderr
							AMPL_G_0_E         = result_0.params['A'].stderr
							SGMA_G_0_E         = result_0.params['SIGMA'].stderr

							CTRE_G_0_cor       = result_0.params['X_0'].correl
							AMPL_G_0_cor       = result_0.params['A'].correl
							SGMA_G_0_cor       = result_0.params['SIGMA'].correl

							chisqr_0           = result_0.chisqr
							redchi_0           = result_0.redchi
						except (RuntimeError,ValueError,TypeError):
							popt_0, pcov_0   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_0           = [999999.99999,999999.99999,999999.99999]
							CTRE_G_0         = 999999.99999
							AMPL_G_0         = 999999.99999
							SGMA_G_0         = 999999.99999
							FWHM_G_0         = 999999.99999
							EW_0             = 999999.99999
							EWE_0            = 999999.99999

							CTRE_G_0_E      = 999999.99999
							AMPL_G_0_E      = 999999.99999
							SGMA_G_0_E      = 999999.99999

							CTRE_G_0_cor    = 999999.99999
							AMPL_G_0_cor    = 999999.99999
							SGMA_G_0_cor    = 999999.99999

							chisqr_0        = 999999.99999
							redchi_0        = 999999.99999
						if fit_vls_hdr == True and fix_ctr_gau ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF0M',float(CTRE_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF0M',float(AMPL_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FF0M',float(FWHM_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF0M',float(EW_0)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EF0M',float(EWE_0)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('The fit (CTR-0) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-0) values will NOT be added to the fits headers!','magenta')
							print
						#################################################CENTRAL GAUSSIAN-0##################################################
						#################################################CENTRAL GAUSSIAN-C##################################################
						try:
							gmodel_O           = Model(func_1D_Gaussian_O)
							gmodel_O.set_param_hint('X_0'   , value=X0_f2DG   , min=X0_f2DG-(X0_f2DG*LINES[7][lines]), max=X0_f2DG+(X0_f2DG*LINES[7][lines]))
							gmodel_O.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*LINES[10][lines]) , max=A_f2DG  + (A_f2DG*LINES[10][lines]))
							gmodel_O.set_param_hint('SIGMA' , value=SIGMA_f2DG)
							gmodel_O.set_param_hint('OFFSET', value=max(inten_glx[mask_ft])-1)
							pars_O             = gmodel_O.make_params()

							result_O           = gmodel_O.fit(inten_glx[mask_ft],pars_O,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,OFFSET=max(inten_glx[mask_ft])-1,
													nan_policy = 'omit')
							CTRE_G_O           = result_O.params['X_0'].value
							AMPL_G_O           = result_O.params['A'].value
							SGMA_G_O           = abs(result_O.params['SIGMA'].value)
							OFST_G_O           = abs(result_O.params['OFFSET'].value)
							FWHM_G_O           = lw_sgma2fwhm(SGMA_G_O)
							W_O                = integrate.quad(lambda x: AMPL_G_O*np.exp(-((x)**2)/(2*SGMA_G_O**2)), -np.inf, np.inf)
							EW_O               = np.round(abs(np.asarray(W_O[0])),10)
							EWE_O              = np.round(abs(np.asarray(W_O[1])),10)
							data_fitted_O      = func_1D_Gaussian_O((lambda_glx[mask_ft]),CTRE_G_O,AMPL_G_O,SGMA_G_O,OFST_G_O)

							CTRE_G_O_E         = result_O.params['X_0'].stderr
							AMPL_G_O_E         = result_O.params['A'].stderr
							SGMA_G_O_E         = result_O.params['SIGMA'].stderr

							CTRE_G_O_cor       = result_O.params['X_0'].correl
							AMPL_G_O_cor       = result_O.params['A'].correl
							SGMA_G_O_cor       = result_O.params['SIGMA'].correl

							chisqr_O           = result_O.chisqr
							redchi_O           = result_O.redchi
							
							#####################################################################################################################
							if ofs_ctr_fit == True:
								print
								print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
								print
								lmb_min_lim_line_ft = (CTRE_G_O-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
								lmb_max_lim_line_ft = (CTRE_G_O+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
								#lmb_min_lim_line    = (CTRE_G_O-LINES[8][lines])*(1+z_glx_Ps) - 1.5*MSK_NTMS*LINES[2][lines]#-20#LINES[2][lines] - 10
								#lmb_max_lim_line    = (CTRE_G_O+LINES[8][lines])*(1+z_glx_Ps) + 1.5*MSK_NTMS*LINES[2][lines]#+20#LINES[2][lines] + 10
								#mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
								mask_ft     = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
							else:
								pass
							#####################################################################################################################

							
							initial_guess_C    = (X0_f2DG,A_f2DG,SIGMA_f2DG)

							gmodel_C           = Model(func_1D_Gaussian)
							gmodel_C.set_param_hint('X_0'  , value=X0_f2DG   , min=X0_f2DG-(X0_f2DG*LINES[7][lines]), max=X0_f2DG+(X0_f2DG*LINES[7][lines]))
							gmodel_C.set_param_hint('A'    , value=A_f2DG    , min=A_f2DG  - (A_f2DG*LINES[10][lines]) , max=A_f2DG  + (A_f2DG*LINES[10][lines]))
							gmodel_C.set_param_hint('SIGMA', value=SIGMA_f2DG)
							pars_C             = gmodel_C.make_params()
							result_C           = gmodel_C.fit(inten_glx[mask_ft],pars_C,
													X=lambda_glx[mask_ft],X_0=X0_f2DG,A=A_f2DG,SIGMA=SIGMA_f2DG,
													nan_policy = 'omit')
							CTRE_G_C           = result_C.params['X_0'].value
							AMPL_G_C           = result_C.params['A'].value
							SGMA_G_C           = abs(result_C.params['SIGMA'].value)
							FWHM_G_C           = lw_sgma2fwhm(SGMA_G_C)
							W_C                = integrate.quad(lambda x: AMPL_G_C*np.exp(-((x)**2)/(2*SGMA_G_C**2)), -np.inf, np.inf)

							EW_C               = np.round(abs(np.asarray(W_C[0])),10)
							EWE_C              = np.round(abs(np.asarray(W_C[1])),10)
							data_fitted_C      = func_1D_Gaussian((lambda_glx[mask_ft]), CTRE_G_C,AMPL_G_C,SGMA_G_C)

							CTRE_G_C_E         = result_C.params['X_0'].stderr
							AMPL_G_C_E         = result_C.params['A'].stderr
							SGMA_G_C_E         = result_C.params['SIGMA'].stderr

							CTRE_G_C_cor       = result_C.params['X_0'].correl
							AMPL_G_C_cor       = result_C.params['A'].correl
							SGMA_G_C_cor       = result_C.params['SIGMA'].correl

							AMPL_SNR           = AMPL_G_C
							CTRE_SNR           = CTRE_G_C
							SGMA_SNR           = abs(SGMA_G_C)

							if CTRE_G_C_E == None:
								CTRE_G_C_E = 999999.99999
							else:
								pass
							if AMPL_G_C_E == None:
								AMPL_G_C_E = 999999.99999
							else:
								pass
							if SGMA_G_C_E == None:
								SGMA_G_C_E = 999999.99999
							else:
								pass
							if CTRE_G_C_cor == None:
								CTRE_G_C_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_cor == None:
								AMPL_G_C_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_cor == None:
								SGMA_G_C_cor = 999999.99999
							else:
								pass
							chisqr_C        = result_C.chisqr
							redchi_C        = result_C.redchi
							##inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O #OFFSET +
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C, pcov_C = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C          = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C        = 999999.99999
							AMPL_G_C        = 999999.99999
							SGMA_G_C        = 999999.99999
							FWHM_G_C        = 999999.99999
							EW_C            = 999999.99999
							EWE_C           = 999999.99999

							CTRE_G_C_E      = 999999.99999
							AMPL_G_C_E      = 999999.99999
							SGMA_G_C_E      = 999999.99999
							CTRE_G_C_cor    = 999999.99999
							AMPL_G_C_cor    = 999999.99999
							SGMA_G_C_cor    = 999999.99999
							chisqr_C        = 999999.99999
							redchi_C        = 999999.99999

							popt_O ,pcov_O  = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
							perr_O          = [999999.99999,999999.99999,999999.99999,999999.99999]
							CTRE_G_O        = 999999.99999
							AMPL_G_O        = 999999.99999
							SGMA_G_O        = 999999.99999
							OFST_G_O        = 999999.99999
							FWHM_G_O        = 999999.99999
							EW_O            = 999999.99999
							EWE_O           = 999999.99999

							CTRE_G_O_E      = 999999.99999
							AMPL_G_O_E      = 999999.99999
							SGMA_G_O_E      = 999999.99999
							CTRE_G_O_cor    = 999999.99999
							AMPL_G_O_cor    = 999999.99999
							SGMA_G_O_cor    = 999999.99999
							OFST_G_O_cor    = 999999.99999
							chisqr_O        = 999999.99999
							redchi_O        = 999999.99999

							AMPL_SNR     = 999999.99999
							CTRE_SNR     = 999999.99999
							SGMA_SNR     = 999999.99999
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines]),'yellow')
						print
						print colored(str(LINES[5][lines])+'_CGLC: ' + str(CTRE_G_C)  ,'yellow')
						print colored(str(LINES[5][lines])+'_AGLC: ' + str(AMPL_G_C)  ,'yellow')
						print colored(str(LINES[5][lines])+'_SGLC: ' + str(SGMA_G_C)  ,'yellow')
						print colored(str(LINES[5][lines])+'_FGLC: ' + str(FWHM_G_C)  ,'yellow')
						print colored(str(LINES[5][lines])+'_WGLC: ' + str(EW_C)      ,'yellow')
						print colored(str(LINES[5][lines])+'_EGLC: ' + str(EWE_C),'yellow')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O

						if fit_vls_hdr == True and fix_ctr_gau==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                       ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                       ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLOM',float(CTRE_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALOM',float(AMPL_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLOM',float(FWHM_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLOM',float(EW_O)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELOM',float(EWE_O)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OFOM',float(OFST_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLCM',float(CTRE_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALCM',float(AMPL_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLCM',float(SGMA_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLCM',float(FWHM_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLCM',float(EW_C)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELCM',float(EWE_C)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))

							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CECM',float(CTRE_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AECM',float(AMPL_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SECM',float(SGMA_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CECM',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AECM',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SECM',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))


							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHLM',float(chisqr_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRLM',float(redchi_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 Reduced 1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('The fit (CTR-C & CTR_O) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-C & CTR_O) values will NOT be added to the fits headers!','magenta')
							print
						#################################################CENTRAL GAUSSIAN-C##################################################					
					elif fix_ctr_gau == True:
						print
						print colored('0 CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						try:
							CTRE_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF0M')
							AMPL_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF0M')
							FWHM_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF0M')
							EW_0        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF0M')
							EWE_0       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF0M')

							CTRE_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLOM')
							AMPL_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALOM')
							FWHM_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLOM')
							EW_O        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLOM')
							EWE_O       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELOM')
							OFST_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFOM')

							CTRE_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLCM')
							AMPL_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALCM')
							SGMA_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLCM')
							FWHM_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLCM')
							EW_C        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLCM')
							EWE_C       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELCM')
							CTRE_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CECM')
							AMPL_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AECM')
							SGMA_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SECM')

							chisqr_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHLM')
							redchi_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRLM')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGF0')
							AMPL_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGF0')
							FWHM_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGF0')
							EW_0        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGF0')
							EWE_0       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGF0')

							CTRE_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGLO')
							AMPL_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGLO')
							FWHM_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGLO')
							EW_O        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGLO')
							EWE_O       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGLO')
							OFST_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFSO')

							CTRE_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGLC')
							AMPL_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGLC')
							SGMA_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGLC')
							FWHM_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGLC')
							EW_C        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGLC')
							EWE_C       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGLC')
							CTRE_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLEC')
							AMPL_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALEC')
							SGMA_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLEC')

							chisqr_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHGL')
							redchi_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRGL')
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-CTR','green')
					print
					print colored(str(LINES[5][lines])+'_CLCM: ' + str(CTRE_G_C)  ,'green')
					print colored(str(LINES[5][lines])+'_ALCM: ' + str(AMPL_G_C)  ,'green')
					print colored(str(LINES[5][lines])+'_SLCM: ' + str(SGMA_G_C)  ,'green')
					print colored(str(LINES[5][lines])+'_FLCM: ' + str(FWHM_G_C)  ,'green')
					print colored(str(LINES[5][lines])+'_WLCM: ' + str(EW_C)      ,'green')
					print colored(str(LINES[5][lines])+'_ELCM: ' + str(EWE_C),'green')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','green')
					print colored(str(CTRE_G_C)+', '+str(AMPL_G_C)+', '+str(SGMA_G_C),'green')
					print
					##################################################CENTRAL GAUSSIAN###################################################
					#####################################################PRE GAUSSIAN####################################################
					if fix_pre_gau == False:
						###########################################DEFINING PRE-PST REGIONS##################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
							lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						else:
							print
							print colored('Using Expected Line Center to Define Fitting Region!','yellow')
							print
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]

						print
						print colored('Region limits for fitting.','yellow')
						print colored('Central    : ' + str(LINES[0][lines]),'yellow')
						print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'yellow')
						print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'yellow')
						print
						print colored('Limits:','yellow')
						print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
						print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim
						print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
						print
						print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
						print
						###########################################DEFINING PRE-PST REGIONS##################################################
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						x_a = lambda_glx[mask_ft_pre][X0_f2DG_indx_PRE]
						y_a = inten_glx[mask_ft_pre][X0_f2DG_indx_PRE]
						try:
							print
							print colored('1-Fitting gaussian before line','cyan')
							print
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							A_f2DG_PRE             = (max(inten_glx[mask_ft_pre])-1)
							X0_f2DG_PRE            = X0_f2DG-pre_shf_ctr#-2.5#lambda_glx[mask_ft_pre][X0_f2DG_indx_PRE]#X0_f2DG-pre_shf_lim#2.5#-2.5#SIGMA_f2DG/2#pre_shf_lim#lambda_glx[mask_ft_pre][X0_f2DG_indx_PRE]#X0_f2DG-2.5#
							SIGMA_f2DG_PRE         = SIGMA_f2DG/2.5
							initial_guess_C_PRE    = (X0_f2DG_PRE,A_f2DG_PRE,SIGMA_f2DG_PRE)

							x_a = lambda_glx[mask_ft_pre][X0_f2DG_indx_PRE]
							y_a = inten_glx[mask_ft_pre][X0_f2DG_indx_PRE]

							print
							print colored('Initial Guess Values PRE : ','cyan')
							print colored(initial_guess_C_PRE,'cyan')
							print

							gmodel_C_PRE           = Model(func_1D_Gaussian_Emm)
							gmodel_C_PRE.set_param_hint('X_0'  , value=X0_f2DG_PRE , min=X0_f2DG_PRE-(X0_f2DG_PRE*LINES[7][lines]), max=X0_f2DG_PRE+(X0_f2DG_PRE*LINES[7][lines]))
							gmodel_C_PRE.set_param_hint('A'    , value=A_f2DG_PRE  , min=A_f2DG_PRE -(A_f2DG_PRE*LINES[10][lines]), max=A_f2DG_PRE +(A_f2DG_PRE*LINES[10][lines]))
							gmodel_C_PRE.set_param_hint('SIGMA', value=SIGMA_f2DG_PRE)
							pars_C_PRE             = gmodel_C_PRE.make_params()
							result_C_PRE           = gmodel_C_PRE.fit(inten_glx[mask_ft_pre],pars_C_PRE,
													X=lambda_glx[mask_ft_pre],X_0=X0_f2DG_PRE,A=A_f2DG_PRE,SIGMA=SIGMA_f2DG_PRE,
													nan_policy = 'omit')
							CTRE_G_C_PRE           = result_C_PRE.params['X_0'].value
							AMPL_G_C_PRE           = result_C_PRE.params['A'].value
							SGMA_G_C_PRE           = abs(result_C_PRE.params['SIGMA'].value)
							FWHM_G_C_PRE           = lw_sgma2fwhm(SGMA_G_C_PRE)
							W_C_PRE                = integrate.quad(lambda x: AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2)), -np.inf, np.inf)

							EW_C_PRE               = np.round(abs(np.asarray(W_C_PRE[0])),10)
							EWE_C_PRE              = np.round(abs(np.asarray(W_C_PRE[1])),10)
							data_fitted_C_PRE      = func_1D_Gaussian_Emm((lambda_glx[mask_ft_pre]), CTRE_G_C_PRE,AMPL_G_C_PRE,SGMA_G_C_PRE)

							CTRE_G_C_PRE_E         = result_C_PRE.params['X_0'].stderr
							AMPL_G_C_PRE_E         = result_C_PRE.params['A'].stderr
							SGMA_G_C_PRE_E         = result_C_PRE.params['SIGMA'].stderr

							CTRE_G_C_PRE_cor       = result_C_PRE.params['X_0'].correl
							AMPL_G_C_PRE_cor       = result_C_PRE.params['A'].correl
							SGMA_G_C_PRE_cor       = result_C_PRE.params['SIGMA'].correl

							AMPL_SNR               = AMPL_G_C_PRE
							CTRE_SNR               = CTRE_G_C_PRE
							SGMA_SNR               = abs(SGMA_G_C_PRE)

							if CTRE_G_C_PRE_E == None or np.isnan(CTRE_G_C_PRE_E):
								CTRE_G_C_PRE_E = 999999.99999
							else:
								pass
							if AMPL_G_C_PRE_E == None:
								AMPL_G_C_PRE_E = 999999.99999
							else:
								pass
							if SGMA_G_C_PRE_E == None:
								SGMA_G_C_PRE_E = 999999.99999
							else:
								pass
							if CTRE_G_C_PRE_cor == None:
								CTRE_G_C_PRE_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_PRE_cor == None:
								AMPL_G_C_PRE_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_PRE_cor == None:
								SGMA_G_C_PRE_cor = 999999.99999
							else:
								pass
							chisqr_C_PRE           = result_C_PRE.chisqr
							redchi_C_PRE           = result_C_PRE.redchi


							W_C_PR1    = integrate.quad(lambda x: AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2)), x_a, np.inf)
							EW_C_PR1   = np.round(abs(np.asarray(W_C_PR1[0])),10)
							EWE_C_PR1  = np.round(abs(np.asarray(W_C_PR1[1])),10)
							#inten_glx[mask_ft_plp] = inten_glx[mask_ft_plp] + OFST_G_O #OFFSET +
						#except (RuntimeError,ValueError,TypeError):
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C_PRE, pcov_C_PRE  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_PRE          = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_PRE        = 999999.99999
							AMPL_G_C_PRE        = 999999.99999
							SGMA_G_C_PRE        = 999999.99999
							FWHM_G_C_PRE        = 999999.99999
							EW_C_PRE            = 999999.99999
							EWE_C_PRE           = 999999.99999

							CTRE_G_C_PRE_E      = 999999.99999
							AMPL_G_C_PRE_E      = 999999.99999
							SGMA_G_C_PRE_E      = 999999.99999
							CTRE_G_C_PRE_cor    = 999999.99999
							AMPL_G_C_PRE_cor    = 999999.99999
							SGMA_G_C_PRE_cor    = 999999.99999
							chisqr_C_PRE        = 999999.99999
							redchi_C_PRE        = 999999.99999

							EW_C_PR1            = 999999.99999
							EWE_C_PR1           = 999999.99999
						if fit_vls_hdr == True and fix_pre_gau == False:						
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHGL',float(chisqr_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRGL',float(redchi_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 Reduced 1GF' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGL1',float(CTRE_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGL1',float(AMPL_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGL1',float(SGMA_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGL1',float(FWHM_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGL1',float(EW_C_PRE)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGL1',float(EWE_C_PRE)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLE1',float(CTRE_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_ALE1',float(AMPL_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SLE1',float(SGMA_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_ALE1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SLE1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))


							Header_Add(specfile_glx,str(LINES[5][lines])+'_XA1',float(x_a),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PRE GAU-LNR X1 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YA1',float(y_a),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PRE GAU-LNR Y1 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WPR1',float(EW_C_PR1),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])   + ' EW   PRE' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EPR1',float(EWE_C_PR1),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE  PRE' + str(fit_type))
							print
							print colored('The fit (PRE) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (PRE) values will NOT be added to the fits headers!','magenta')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GS1',float(pre_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for PRE G-1')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GC1',float(pre_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for PRE G-1')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass
					elif fix_pre_gau == True:
						try:
							print
							print colored('1 PRE-gaussian already fitted!','yellow')
							print colored('Values from fits header','yellow')
							chisqr_C       = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHGL')
							redchi_C       = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRGL')
							CTRE_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGL1')
							AMPL_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGL1')
							SGMA_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGL1')
							FWHM_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGL1')
							EW_C_PRE       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGL1')
							EWE_C_PRE      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGL1')
							CTRE_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLE1')
							AMPL_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALE1')
							SGMA_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLE1')

							x_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XA1')
							y_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YA1')
							EW_C_PR1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WPR1')
							EWE_C_PR1      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EPR1')

							pre_shf_lim    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GS1')
							pre_shf_ctr    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GC1')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!' ,'yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
							print colored ('Gotta Fit first before fixing a component (PRE)!','yellow')
							print colored ('Or UnFix (fix_pre_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							quit()
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-PRE','cyan')
					print
					print '******************************************************************************'
					print colored('Boundaries for Gaussian Fitting (PRE-CTR): ' + str(pre_shf_ctr),'cyan')
					print colored('Boundaries for Gaussian Fitting (PRE-LIM): ' + str(pre_shf_lim),'cyan')
					print '******************************************************************************'					
					print
					print colored(str(LINES[5][lines])+'_CGL1: ' + str(CTRE_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_AGL1: ' + str(AMPL_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_SGL1: ' + str(SGMA_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_FGL1: ' + str(FWHM_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_WGL1: ' + str(EW_C_PRE)      ,'cyan')
					print colored(str(LINES[5][lines])+'_EGL1: ' + str(EWE_C_PRE),'cyan')
					print
					print colored(str(LINES[5][lines])+'_XA1 : ' + str(x_a),'cyan')
					print colored(str(LINES[5][lines])+'_YA1 : ' + str(y_a),'cyan')
					print colored(str(LINES[5][lines])+'_WPR1: ' + str(EW_C_PR1),'cyan')
					print colored(str(LINES[5][lines])+'_EPR1: ' + str(EWE_C_PR1),'cyan')
					print colored(str(LINES[5][lines])+'_GS1 : ' + str(pre_shf_lim),'cyan')
					print colored(str(LINES[5][lines])+'_GC1 : ' + str(pre_shf_ctr),'cyan')
					print
					print colored('Fit Values (PRE) Center, Amplitude, Sigma ('+fit_type+'):','cyan')
					print colored(str(CTRE_G_C_PRE)+', '+str(AMPL_G_C_PRE)+', '+str(SGMA_G_C_PRE),'cyan')
					print
					#####################################################PRE GAUSSIAN#################################################
					###################################################POST GAUSSIAN##################################################
					if fix_pst_gau == False:
						###########################################DEFINING PRE-PST REGIONS##################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
							lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						else:
							print
							print colored('Using Expected Line Center to Define Fitting Region!','yellow')
							print
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]

						print
						print colored('Region limits for fitting.','yellow')
						print colored('Central    : ' + str(LINES[0][lines]),'yellow')
						print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'yellow')
						print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'yellow')
						print
						print colored('Limits:','yellow')
						print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
						print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim
						print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
						print
						print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
						print
						###########################################DEFINING PRE-PST REGIONS##################################################						
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						x_b = lambda_glx[mask_ft_pst][X0_f2DG_indx_PST]
						y_b = inten_glx[mask_ft_pst][X0_f2DG_indx_PST]
						try:
							print colored('2-Fitting gaussian after line','cyan')
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
							A_f2DG_PST             = (max(inten_glx[mask_ft_pst])-1)
							X0_f2DG_PST            = X0_f2DG+pst_shf_ctr
							SIGMA_f2DG_PST         = SIGMA_f2DG/2
							initial_guess_C_PST    = (X0_f2DG_PST,A_f2DG_PST,SIGMA_f2DG_PST)

							print
							print colored('Initial Guess Values PST : ','cyan')
							print colored(initial_guess_C_PST,'cyan')
							print

							gmodel_C_PST           = Model(func_1D_Gaussian_Emm)
							gmodel_C_PST.set_param_hint('X_0'  , value=X0_f2DG_PST , min=X0_f2DG_PST-(X0_f2DG_PST*LINES[7][lines]), max=X0_f2DG_PST+(X0_f2DG_PST*LINES[7][lines]))
							gmodel_C_PST.set_param_hint('A'    , value=A_f2DG_PST  , min=A_f2DG_PST -(A_f2DG_PST*LINES[10][lines]), max=A_f2DG_PST +(A_f2DG_PST*LINES[10][lines]))#min=A_f2DG_PST-0.001, max=A_f2DG_PST)
							gmodel_C_PST.set_param_hint('SIGMA', value=SIGMA_f2DG_PST)
							pars_C_PST             = gmodel_C_PST.make_params()
							result_C_PST           = gmodel_C_PST.fit(inten_glx[mask_ft_pst],pars_C_PST,
													X=lambda_glx[mask_ft_pst],X_0=X0_f2DG_PST,A=A_f2DG_PST,SIGMA=SIGMA_f2DG_PST,
													nan_policy = 'omit')
							CTRE_G_C_PST           = result_C_PST.params['X_0'].value
							AMPL_G_C_PST           = result_C_PST.params['A'].value
							SGMA_G_C_PST           = abs(result_C_PST.params['SIGMA'].value)
							FWHM_G_C_PST           = lw_sgma2fwhm(SGMA_G_C_PST)
							W_C_PST                = integrate.quad(lambda x: AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)), -np.inf, np.inf)

							EW_C_PST               = np.round(abs(np.asarray(W_C_PST[0])),10)
							EWE_C_PST              = np.round(abs(np.asarray(W_C_PST[1])),10)
							data_fitted_C_PST      = func_1D_Gaussian_Emm((lambda_glx[mask_ft_pst]), CTRE_G_C_PST,AMPL_G_C_PST,SGMA_G_C_PST)

							CTRE_G_C_PST_E         = result_C_PST.params['X_0'].stderr
							AMPL_G_C_PST_E         = result_C_PST.params['A'].stderr
							SGMA_G_C_PST_E         = result_C_PST.params['SIGMA'].stderr

							CTRE_G_C_PST_cor       = result_C_PST.params['X_0'].correl
							AMPL_G_C_PST_cor       = result_C_PST.params['A'].correl
							SGMA_G_C_PST_cor       = result_C_PST.params['SIGMA'].correl

							AMPL_SNR           = AMPL_G_C_PST
							CTRE_SNR           = CTRE_G_C_PST
							SGMA_SNR           = abs(SGMA_G_C_PST)

							if CTRE_G_C_PST_E == None:
								CTRE_G_C_PST_E = 999999.99999
							else:
								pass
							if AMPL_G_C_PST_E == None:
								AMPL_G_C_PST_E = 999999.99999
							else:
								pass
							if SGMA_G_C_PST_E == None:
								SGMA_G_C_PST_E = 999999.99999
							else:
								pass
							if CTRE_G_C_PST_cor == None:
								CTRE_G_C_PST_cor = 999999.99999
							else:
								pass
							if AMPL_G_C_PST_cor == None:
								AMPL_G_C_PST_cor = 999999.99999
							else:
								pass
							if SGMA_G_C_PST_cor == None:
								SGMA_G_C_PST_cor = 999999.99999
							else:
								pass
							chisqr_C_PST           = result_C_PST.chisqr
							redchi_C_PST           = result_C_PST.redchi

							W_C_PS2    = integrate.quad(lambda x: AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)), -np.inf, x_b)
							EW_C_PS2   = np.round(abs(np.asarray(W_C_PS2[0])),10)
							EWE_C_PS2  = np.round(abs(np.asarray(W_C_PS2[1])),10)
						except (RuntimeError,ValueError,TypeError):
							print colored('RuntimeError','cyan')
							popt_C_PST, pcov_C_PST  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_PST          = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_PST        = 999999.99999
							AMPL_G_C_PST        = 999999.99999
							SGMA_G_C_PST        = 999999.99999
							FWHM_G_C_PST        = 999999.99999
							EW_C_PST            = 999999.99999
							EWE_C_PST           = 999999.99999

							
							CTRE_G_C_PST_E      = 999999.99999
							AMPL_G_C_PST_E      = 999999.99999
							SGMA_G_C_PST_E      = 999999.99999
							CTRE_G_C_PST_cor    = 999999.99999
							AMPL_G_C_PST_cor    = 999999.99999
							SGMA_G_C_PST_cor    = 999999.99999
							chisqr_C_PST        = 999999.99999
							redchi_C_PST        = 999999.99999

							EW_C_PR2            = 999999.99999
							EWE_C_PR2           = 999999.99999

						if fit_vls_hdr == True and fix_pst_gau == False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGL2',float(CTRE_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGL2',float(AMPL_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGL2',float(SGMA_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGL2',float(FWHM_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGL2',float(EW_C_PST)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGL2',float(EWE_C_PST)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CLE2',float(CTRE_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_ALE2',float(AMPL_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SLE2',float(SGMA_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CLE2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_ALE2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SLE2',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))


							Header_Add(specfile_glx,str(LINES[5][lines])+'_XA2',float(x_b),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PST GAU-LNR X2 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YA2',float(y_b),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PST GAU-LNR Y2 COO')
							#Header_Add(specfile_glx,str(LINES[5][lines])+'_WPS2',float(EW_C_PS2),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])   + ' EW   PST' + str(fit_type))
							#Header_Add(specfile_glx,str(LINES[5][lines])+'_EPS2',float(EWE_C_PS2),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE  PST' + str(fit_type))
						else:
							print
							print colored('The fit (PST) values will NOT be added to the fits headers!','magenta')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GS2',float(pst_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for PST G-2')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GC2',float(pst_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for PST G-2')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass
						#inten_glx[mask_ft_plp] = inten_glx[mask_ft_plp] + OFST_G_O #OFFSET -
					elif fix_pst_gau == True:
						try:
							print
							print colored('2 PST-gaussian already fitted!','yellow')
							print colored('Values from fits header','yellow')
							CTRE_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGL2')
							AMPL_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGL2')
							SGMA_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGL2')
							FWHM_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGL2')
							EW_C_PST       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGL2')
							EWE_C_PST      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGL2')
							CTRE_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLE2')
							AMPL_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALE2')
							SGMA_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLE2')

							x_b         = Header_Get(specfile_glx,str(LINES[5][lines])+'_XA2')
							y_b         = Header_Get(specfile_glx,str(LINES[5][lines])+'_YA2')
							EW_C_PS2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_WPS2')
							EWE_C_PS2   = Header_Get(specfile_glx,str(LINES[5][lines])+'_EPS2')

							pst_shf_lim    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GS2')
							pst_shf_ctr    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GC2')
						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
							print colored ('Gotta Fit first before fixing a component (PST)!','yellow')
							print colored ('Or UnFix (fix_pst_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							quit()
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-PST','magenta')
					print
					print '******************************************************************************'
					print colored('Boundaries for Gaussian Fitting (PST-CTR): ' + str(pst_shf_ctr),'magenta')					
					print colored('Boundaries for Gaussian Fitting (PST-LIM): ' + str(pst_shf_lim),'magenta')
					print '******************************************************************************'
					print
					print colored(str(LINES[5][lines])+'_CGL2: ' + str(CTRE_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_AGL2: ' + str(AMPL_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_SGL2: ' + str(SGMA_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_FGL2: ' + str(FWHM_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_WGL2: ' + str(EW_C_PST)      ,'magenta')
					print colored(str(LINES[5][lines])+'_EGL2: ' + str(EWE_C_PST),'magenta')
					print
					print colored(str(LINES[5][lines])+'_GS2 : ' + str(pst_shf_lim),'magenta')
					print colored(str(LINES[5][lines])+'_GC2 : ' + str(pst_shf_ctr),'magenta')
					print colored(str(LINES[5][lines])+'_XA2 : ' + str(x_b),'magenta')
					print colored(str(LINES[5][lines])+'_YA2 : ' + str(y_b),'magenta')
					#print colored(str(LINES[5][lines])+'_WPS2: ' + str(EW_C_PS2),'magenta')
					#print colored(str(LINES[5][lines])+'_EPS2: ' + str(EWE_C_PS2),'magenta')
					print
					print colored('Fit Values (PST) Center, Amplitude, Sigma ('+fit_type+'):','magenta')
					print colored(str(CTRE_G_C_PST)+', '+str(AMPL_G_C_PST)+', '+str(SGMA_G_C_PST),'magenta')
					print					
					###################################################POST GAUSSIAN##################################################										

					###############################################COMPUTING TOTAL AREA###############################################
					print colored('Computing Flux Area','yellow')
					#############################################COMPUTING LINEAR AREA##################################################
					slope_line1 = (y_a-y_b)/(x_a-x_b)
					slope_line2 = (y_b-y_a)/(x_b-x_a)
					b1 = y_a - (slope_line1*x_a)
					b2 = y_b - (slope_line1*x_b)
					print
					print colored('Linear Parameters considering peak points:','yellow')
					print 'Point A: ',x_a,y_a
					print 'Point B: ',x_b,y_b
					print 'Slope: ',slope_line1
					print 'Slope: ',slope_line2
					print 'b: ',b1,b2
					print
					#############################################COMPUTING LINEAR AREA##################################################
					################################################COMPUTING TOTAL AREA################################################
					CTRE_G_C_PRE = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGL1')
					AMPL_G_C_PRE = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGL1')
					SGMA_G_C_PRE = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGL1')

					CTRE_G_C     = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLCM')
					AMPL_G_C     = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALCM')
					SGMA_G_C     = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLCM')

					CTRE_G_C_PST = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGL2')
					AMPL_G_C_PST = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGL2')
					SGMA_G_C_PST = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGL2')

					print
					print colored('Computing Areas using info from fits headers.','yellow')
					print 'AMPL_G_C: ',str(AMPL_G_C),'SGMA_G_C: ',str(SGMA_G_C)
					print 'AMPL_G_C_PRE: ',str(AMPL_G_C_PRE),'SGMA_G_C_PRE: ',str(SGMA_G_C_PRE)
					print 'AMPL_G_C_PST: ',str(AMPL_G_C_PST),'SGMA_G_C_PST: ',str(SGMA_G_C_PST)
					print

					W_C       = integrate.quad(lambda x: AMPL_G_C*np.exp(-((x)**2)/(2*SGMA_G_C**2)), -np.inf, np.inf)
					EW_C      = np.round((np.asarray(W_C[0])),10)
					EWE_C     = np.round(abs(np.asarray(W_C[1])),10)

					W_C_PRE   = integrate.quad(lambda x: AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2)), -np.inf, np.inf)
					EW_C_PRE  = np.round((np.asarray(W_C_PRE[0])),10)
					EWE_C_PRE = np.round(abs(np.asarray(W_C_PRE[1])),10)

					W_C_PST   = integrate.quad(lambda x: AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)), -np.inf, np.inf)
					EW_C_PST  = np.round((np.asarray(W_C_PST[0])),10)
					EWE_C_PST = np.round(abs(np.asarray(W_C_PST[1])),10)

					###############ALTERNATTIVE##############
					#####ROOTS ARE NEEDED FOR INT LIMITS#####
					#W_PLP   = integrate.quad(lambda x: AMPL_G_C*np.exp(-((x)**2)/(2*SGMA_G_C**2)) +
												#(AMPL_G_C_PST*np.exp(-((x-(CTRE_G_C_PST-CTRE_G_C))**2)/(2*SGMA_G_C_PST**2)))+
												#(AMPL_G_C_PRE*np.exp(-((x-(CTRE_G_C_PRE-CTRE_G_C))**2)/(2*SGMA_G_C_PRE**2))), 
												##-np.inf, np.inf
												#lmb_min_lim_line_ft-CTRE_G_C,lmb_max_lim_line_ft-CTRE_G_C
												#)
					#EW_PLP  = np.round((np.asarray(W_PLP[0])),10)
					#EWE_PLP = np.round(abs(np.asarray(W_PLP[1])),10)
					#EWMT    = EW_C + EW_C_PRE + EW_C_PST #+ EW_C_LNR
					#print CTRE_G_C
					#print CTRE_G_C_PRE 
					#print CTRE_G_C_PST
					#print (CTRE_G_C_PST-CTRE_G_C)
					#print (CTRE_G_C_PRE-CTRE_G_C)
					#print lmb_min_lim_line_ft-CTRE_G_C
					#print lmb_max_lim_line_ft-CTRE_G_C
					###############ALTERNATTIVE##############
					W_PLP   = integrate.quad(lambda x: AMPL_G_C*np.exp(-((x)**2)/(2*SGMA_G_C**2)) +
												(AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)))+
												(AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2))), 
												-np.inf, np.inf
												)
					EW_PLP  = np.round((np.asarray(W_PLP[0])),10)
					EWE_PLP = np.round(abs(np.asarray(W_PLP[1])),10)

					if EW_C_PRE == 999999.99999:
						EW_C_PRE = 0
					else:
						pass
					if EWE_C_PST == 999999.99999:
						EWE_C_PST = 0
					else:
						pass
					if EW_C == 999999.99999:
						EW_C = 0
					else:
						pass

					EWMT    = EW_C  + EW_C_PRE  + EW_C_PST  #+ EW_C_LNR
					EWEMT   = EWE_C + EWE_C_PRE + EWE_C_PST #+ EW_C_LNR

					print
					print colored('Areas     :','yellow')
					print colored('Area CTR-G: ' + str(EW_C),'yellow')
					print colored('Area PRE-G: ' + str(EW_C_PRE),'blue')
					print colored('Area PST-G: ' + str(EW_C_PST),'magenta')
					print colored('Area PRE-CTR-PST: ' + str(EW_PLP),'yellow')
					print colored('Area TOT=PRE-CTR-PST: ' + str(EWMT),'yellow')
					print
					################################################COMPUTING TOTAL AREA################################################
					#############################################ADDING AREA TO FTIS HEADER#############################################
					if fit_vls_hdr == True:
						print
						print colored('The Areas values will be updated to the fits headers!','magenta')
						print
						print colored('Area CTR-G: '           + str(EW_C)     + '-' +str(LINES[5][lines])+'_WLCM','yellow')
						print colored('Area PRE-G: '           + str(EW_C_PRE) + '-' +str(LINES[5][lines])+'_ELCM','yellow')
						print colored('Area PST-G: '           + str(EW_C_PST) + '-' +str(LINES[5][lines])+'_WPST','yellow')
						print colored('Area PRE-CTR-PST: '     + str(EW_PLP)   + '-' +str(LINES[5][lines])+'_WPLP','yellow')
						print colored('Area TOT=PRE-CTR-PST: ' + str(EWMT)     + '-' +str(LINES[5][lines])+'_WTOT','yellow')
						print
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WLCM',float(EW_C)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EW CTR'  + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_ELCM',float(EWE_C)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE CTR' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGL1',float(EW_C_PRE) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EW PRE'  + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGL1',float(EWE_C_PRE),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE PRE' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGL2',float(EW_C_PST) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EW PST'  + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGL2',float(EWE_C_PST),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE PST' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WPLP',float(EW_PLP)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EW CTR-PRE-PST'  + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WEPL',float(EWE_PLP)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE CTR-PRE-PST' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WTOT',float(EWMT)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EW TOT'  + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WETT',float(EWMT)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE TOT' + str(fit_type))
					else:
						print
						print colored('The LINEAR AND TOTAL Areas values will NOT be updated to the fits headers!','magenta')
						print
					#############################################ADDING AREA TO FTIS HEADER#############################################						
					#################################DEFINE AREA FOR LATER ON SUBSTRACT OFFSET FOR PLOTING###############################
					###########################################DEFINING PRE-PST REGIONS##################################################
					if ofs_ctr_fit == True:
						print
						print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
						print
						lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
						lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
						mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
						mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
					else:
						print
						print colored('Using Expected Line Center to Define Fitting Region!','yellow')
						print
						mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
						mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]

					print
					print colored('Region limits for fitting.','yellow')
					print colored('Central    : ' + str(LINES[0][lines]),'yellow')
					print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'yellow')
					print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'yellow')
					print
					print colored('Limits:','yellow')
					print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
					print 'PRE    : lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim,'-','lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim
					print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
					print
					print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
					print
					###########################################DEFINING PRE-PST REGIONS##################################################
					#################################DEFINE AREA FOR LATER ON SUBSTRACT OFFSET FOR PLOTING###############################
				elif 'Dbl' in LINES[3][lines] and fit_fnct=='gauss' and fit_type == 'lmfit' and uft_lne_vls == True:
					print 'Line-fitting Cleaning Method. To be checked Fnc_Stk_Plt.py def(Plot_Idp_Spc_Lne) line 10731!'
					fit_typ = 'G'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					print 
					print colored('Double Line Fit (Ind)','yellow')
					print colored(LINES[3][lines-2]  + '-' + str(LINES[0][lines-2])   + '-' + str(LINES[1][lines-2]),'cyan')
					print LINES[3][lines]+ '-' + str(LINES[0][lines]) + '-'  + str(LINES[1][lines])
					print colored(LINES[3][lines-1]+ '-' + str(LINES[0][lines-1]) + '-' + str(LINES[1][lines-1]),'magenta')

					from lmfit import Model
	
					########################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#########################
					if ivl_fts_hdr == True:
						try:
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial fit variables from fits header!','yellow')
							print
							print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WF02' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WP02' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CF02' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CO02' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AF02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print
						except ValueError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print '*****'
							print
							quit()
						try:
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')							
							print colored(str(LINES[5][lines-2])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')							
							print colored(str(LINES[5][lines-1])+'_AF02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF01',0.001,header_comment = str(LINES[3][lines-2]) + str(LINES[0][lines-2]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF02',0.001,header_comment = str(LINES[3][lines-1]) + str(LINES[0][lines-1]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines-2])+'_AF01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines-1])+'_AF02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('*****Success!******','yellow')
							print
						################FOR CASES WHERE KLINE WIIDHT ==0################
						if L1_1  == 0:
							L1_1  = 1#LINES[1][lines]
						else:
							pass
						if L1_2  == 0:
							L1_2  = 1#LINES[1][lines]
						else:
							pass
						################FOR CASES WHERE KLINE WIIDHT ==0################
					elif ivl_fts_hdr == False:
						L1_1  = LINES[1][lines-2]
						L2_1  = LINES[2][lines-2]
						L7_1  = LINES[7][lines-2]
						L8_1  = LINES[8][lines-2]
						L10_1 = LINES[10][lines-2]

						L1_2  = LINES[1][lines-1]
						L2_2  = LINES[2][lines-1]
						L7_2  = LINES[7][lines-1]
						L8_2  = LINES[8][lines-1]
						L10_2 = LINES[10][lines-1]
						print
						print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
						print					
					print colored('Initial Values: ','cyan')
					print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')					
					print colored(str(LINES[5][lines])+'_WF01' + ': ' + str(L1_1),'cyan')
					print colored(str(LINES[5][lines])+'_WP01' + ': ' + str(L2_1),'cyan')
					print colored(str(LINES[5][lines])+'_CF01' + ': ' + str(L7_1),'cyan')
					print colored(str(LINES[5][lines])+'_CO01' + ': ' + str(L8_1),'cyan')
					print colored(str(LINES[5][lines])+'_AF01' + ': ' + str(L10_1),'cyan')
					print
					print colored('Initial Values: ','magenta')
					print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')					
					print colored(str(LINES[5][lines])+'_WF02' + ': ' + str(L1_2),'magenta')
					print colored(str(LINES[5][lines])+'_WP02' + ': ' + str(L2_2),'magenta')
					print colored(str(LINES[5][lines])+'_CF02' + ': ' + str(L7_2),'magenta')
					print colored(str(LINES[5][lines])+'_CO02' + ': ' + str(L8_2),'magenta')
					print colored(str(LINES[5][lines])+'_AF02' + ': ' + str(L10_2),'magenta')
					print
					########################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#########################
					##################################################CENTRAL GAUSSIAN-1###################################################
					lmb_min_lim_line_ft = (LINES[0][lines-2]+L8_1) - MSK_NTMS*L1_1 
					lmb_max_lim_line_ft = (LINES[0][lines-2]+L8_1) + MSK_NTMS*L1_1
					lmb_min_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_1#- 20#L2_1 - 10 
					lmb_max_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_1#+ 20#L2_1 + 10

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function

					idx_ctr_ft_reg = np.abs(lambda_glx[mask_ft] - (LINES[0][lines-2]+L8_1)).argmin()

					X0_f2DG    = (LINES[0][lines-2]+L8_1)
					SIGMA_f2DG = L1_1
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))
					A_f2DG     = -(1-inten_glx[mask_ft][idx_ctr_ft_reg])

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
				
					#################################################CENTRAL GAUSSIAN-1-C##################################################
					if fix_ctr_gau_1 == False:
						print
						print colored('Fitting 1st line','cyan')
						print colored('1-0-Fitting Central line','cyan')
						print
						#################################################CENTRAL GAUSSIAN-1-0##################################################
						popt_0_1, pcov_0_1   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_0_1             = [999999.99999,999999.99999,999999.99999]
						CTRE_G_0_1           = 999999.99999
						AMPL_G_0_1           = 999999.99999
						SGMA_G_0_1           = 999999.99999
						FWHM_G_0_1           = 999999.99999
						EW_0_1               = 999999.99999
						EWE_0_1              = 999999.99999

						CTRE_G_0_1_E    = 999999.99999
						AMPL_G_0_1_E    = 999999.99999
						SGMA_G_0_1_E    = 999999.99999

						CTRE_G_0_cor    = 999999.99999
						AMPL_G_0_cor    = 999999.99999
						SGMA_G_0_cor    = 999999.99999

						chisqr_0_1        = 999999.99999
						redchi_0_1        = 999999.99999
						if fit_vls_hdr == True and fix_ctr_gau_1 ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF01',float(CTRE_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF01',float(AMPL_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FF01',float(FWHM_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF01',float(EW_0_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EF01',float(EWE_0_1)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('**************************CLEANING**************************','yellow')
							print colored('The fit (CTR-1-0) values will be added to the fits headers!','cyan')
							print colored('**************************CLEANING**************************','yellow')
							print
						else:
							print
							print colored('The fit (CTR-1-0) values will NOT be added to the fits headers!','cyan')
							print
						#################################################CENTRAL GAUSSIAN-1-0##################################################
						#################################################CENTRAL GAUSSIAN-1-C##################################################
						print
						print colored('Clean Line','cyan')
						print
						popt_C_1, pcov_C_1 = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C_1           = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C_1         = 999999.99999
						AMPL_G_C_1         = 999999.99999
						SGMA_G_C_1         = 999999.99999
						FWHM_G_C_1         = 999999.99999
						EW_C_1             = 999999.99999
						EWE_C_1            = 999999.99999

						CTRE_G_C_1_E       = 999999.99999
						AMPL_G_C_1_E       = 999999.99999
						SGMA_G_C_1_E       = 999999.99999
						CTRE_G_C_1_cor     = 999999.99999
						AMPL_G_C_1_cor     = 999999.99999
						SGMA_G_C_1_cor     = 999999.99999
						chisqr_C_1         = 999999.99999
						redchi_C_1         = 999999.99999

						popt_O_1 ,pcov_O_1 = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
						perr_O_1           = [999999.99999,999999.99999,999999.99999,999999.99999]
						CTRE_G_O_1         = 999999.99999
						AMPL_G_O_1         = 999999.99999
						SGMA_G_O_1         = 999999.99999
						OFST_G_O_1         = 999999.99999
						FWHM_G_O_1         = 999999.99999
						EW_O_1             = 999999.99999
						EWE_O_1            = 999999.99999

						CTRE_G_O_1_E      = 999999.99999
						AMPL_G_O_1_E      = 999999.99999
						SGMA_G_O_1_E      = 999999.99999
						CTRE_G_O_1_cor    = 999999.99999
						AMPL_G_O_1_cor    = 999999.99999
						SGMA_G_O_1_cor    = 999999.99999
						OFST_G_O_1_cor    = 999999.99999
						chisqr_O_1        = 999999.99999
						redchi_O_1        = 999999.99999

						AMPL_SNR_1        = 999999.99999
						CTRE_SNR_1        = 999999.99999
						SGMA_SNR_1        = 999999.99999						
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines-2])+'-'+str(LINES[3][lines-2]),'cyan')
						print
						print colored(str(LINES[5][lines-2])+'_CGLC: ' + str(CTRE_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_AGLC: ' + str(AMPL_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_SGLC: ' + str(SGMA_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_FGLC: ' + str(FWHM_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_WGLC: ' + str(EW_C_1)      ,'cyan')
						print colored(str(LINES[5][lines-2])+'_EGLC: ' + str(EWE_C_1)     ,'cyan')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O
						if fit_vls_hdr == True and fix_ctr_gau_1==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                         ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                         ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLO1',float(CTRE_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALO1',float(AMPL_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLO1',float(FWHM_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLO1',float(EW_O_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELO1',float(EWE_O_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OFO1',float(OFST_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLC1',float(CTRE_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALC1',float(AMPL_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLC1',float(SGMA_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLC1',float(FWHM_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLC1',float(EW_C_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELC1',float(EWE_C_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CEC1',float(CTRE_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AEC1',float(AMPL_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SEC1',float(SGMA_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHL1',float(chisqr_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRL1',float(redchi_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines-2]) + ' Chi2 Reduced 1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('*******************************CLEANING*******************************','yellow')
							print colored('The fit (CTR-1-C & CTR-1-O) values will be added to the fits headers!','cyan')
							print colored('*******************************CLEANING*******************************','yellow')
							print
						else:
							print
							print colored('The fit (CTR-1-C & CTR-1-O) values will NOT be added to the fits headers!','cyan')
							print

						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values G-1 for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF01',float(L1_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-1-1 Wdt-Fit  1GF-IntVal')         #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WP01',float(L2_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-2-1 Wdt-Plt  1GF-IntVal')         #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF01',float(L7_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-7-1 Ctr Fit Bnds  1GF-IntVal')    #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CO01',float(L8_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-8-1 Ctr Fit Ofst  1GF-IntVal')    #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF01',float(L10_1),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10-1 Amp Fit Bnds  1GF-IntVal')   #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial Values: ','cyan')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')					
							print colored(str(LINES[5][lines])+'_WF01' + ': ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WP01' + ': ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CF01' + ': ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CO01' + ': ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AF01' + ': ' + str(L10_1),'cyan')
							print
						else:
							print
							print colored('Initial Guess Values G-1 for line Fitting will NOT be recorded!','yellow')
							print
							pass
					#################################################CENTRAL GAUSSIAN-1-C##################################################
					elif fix_ctr_gau_1 == True:
						print
						print colored('1st CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						try:
							CTRE_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF01')
							AMPL_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF01')
							FWHM_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF01')
							EW_0_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF01')
							EWE_0_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF01')

							CTRE_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLO1')
							AMPL_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALO1')
							FWHM_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLO1')
							EW_O_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLO1')
							EWE_O_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELO1')
							OFST_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFO1')

							CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC1')
							AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC1')
							SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC1')
							FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC1')
							EW_C_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLC1')
							EWE_C_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELC1')
							CTRE_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CEC1')
							AMPL_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AEC1')
							SGMA_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SEC1')

							chisqr_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHL1')
							redchi_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRL1')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines-2])+'_CF01','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau_1=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CGF0')
							AMPL_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_AGF0')
							FWHM_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_FGF0')
							EW_0_1        = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_WGF0')
							EWE_0_1       = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_EGF0')

							CTRE_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CGLO')
							AMPL_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_AGLO')
							FWHM_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_FGLO')
							EW_O_1        = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_WGLO')
							EWE_O_1       = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_EGLO')
							OFST_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_OFSO')

							CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CGLC')
							AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_AGLC')
							SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_SGLC')
							FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_FGLC')
							EW_C_1        = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_WGLC')
							EWE_C_1       = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_EGLC')
							CTRE_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CLEC')
							AMPL_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_ALEC')
							SGMA_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_SLEC')

							chisqr_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CHGL')
							redchi_C_1    = Header_Get(specfile_glx,str(LINES[5][lines-2])+'_CRGL')
					print
					print colored(str(LINES[0][lines-2])+'-'+str(LINES[3][lines-2])+'-CTR','cyan')
					print
					print colored(str(LINES[5][lines-2])+'_CLC1: ' + str(CTRE_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines-2])+'_ALC1: ' + str(AMPL_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines-2])+'_SLC1: ' + str(SGMA_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines-2])+'_FLC1: ' + str(FWHM_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines-2])+'_WLC1: ' + str(EW_C_1)      ,'cyan')
					print colored(str(LINES[5][lines-2])+'_ELC1: ' + str(EWE_C_1),'cyan')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','cyan')
					print colored(str(CTRE_G_C_1)+', '+str(AMPL_G_C_1)+', '+str(SGMA_G_C_1),'cyan')
					print
					##################################################CENTRAL GAUSSIAN-1###################################################
					##################################################CENTRAL GAUSSIAN-2###################################################
					lmb_min_lim_line_ft = (LINES[0][lines-1]+L8_2) - MSK_NTMS*LINES[1][lines-1] 
					lmb_max_lim_line_ft = (LINES[0][lines-1]+L8_2) + MSK_NTMS*LINES[1][lines-1]
					lmb_min_lim_line    = (LINES[0][lines-1]+L8_2)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_2#- 20#L2_2 - 10 
					lmb_max_lim_line    = (LINES[0][lines-1]+L8_2)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_2#+ 20#L2_2 + 10

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function
					idx_ctr_ft_reg = np.abs(lambda_glx[mask_ft] - (LINES[0][lines-1]+L8_2)).argmin()

					X0_f2DG    = (LINES[0][lines-1]+L8_2)
					SIGMA_f2DG = LINES[1][lines-1]
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))
					A_f2DG     = -(1-inten_glx[mask_ft][idx_ctr_ft_reg])

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass
					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
					##################################################CENTRAL GAUSSIAN-2###################################################
					if fix_ctr_gau_2 == False:
						print
						print colored('Fitting 2nd line','magenta')
						print colored('2-0-Fitting Central line','magenta')
						print
						#################################################CENTRAL GAUSSIAN-2-0##################################################
						popt_0, pcov_0   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_0           = [999999.99999,999999.99999,999999.99999]
						CTRE_G_0_2         = 999999.99999
						AMPL_G_0_2         = 999999.99999
						SGMA_G_0_2         = 999999.99999
						FWHM_G_0         = 999999.99999
						EW_0_2             = 999999.99999
						EWE_0_2            = 999999.99999

						CTRE_G_0_2_E      = 999999.99999
						AMPL_G_0_2_E      = 999999.99999
						SGMA_G_0_2_E      = 999999.99999

						CTRE_G_0_2_cor    = 999999.99999
						AMPL_G_0_2_cor    = 999999.99999
						SGMA_G_0_2_cor    = 999999.99999

						chisqr_0        = 999999.99999
						redchi_0        = 999999.99999
						if fit_vls_hdr == True and fix_ctr_gau_2 ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF02',float(CTRE_G_0_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF02',float(AMPL_G_0_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FF02',float(FWHM_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF02',float(EW_0_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EF02',float(EWE_0_2)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('**************************CLEANING**************************','yellow')
							print colored('The fit (CTR-0) values will be added to the fits headers!','magenta')
							print colored('**************************CLEANING**************************','yellow')
							print
						else:
							print
							print colored('The fit (CTR-0) values will NOT be added to the fits headers!','magenta')
							print
						#################################################CENTRAL GAUSSIAN-2-0##################################################
						#################################################CENTRAL GAUSSIAN-2-C##################################################
						print colored('Clean Line','cyan')
						popt_C_2, pcov_C_2 = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C_2           = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C_2         = 999999.99999
						AMPL_G_C_2         = 999999.99999
						SGMA_G_C_2         = 999999.99999
						FWHM_G_C_2         = 999999.99999
						EW_C_2             = 999999.99999
						EWE_C_2            = 999999.99999

						CTRE_G_C_2_E       = 999999.99999
						AMPL_G_C_2_E       = 999999.99999
						SGMA_G_C_2_E       = 999999.99999
						CTRE_G_C_2_cor     = 999999.99999
						AMPL_G_C_2_cor     = 999999.99999
						SGMA_G_C_2_cor     = 999999.99999
						chisqr_C_2         = 999999.99999
						redchi_C_2         = 999999.99999

						popt_O_2 ,pcov_O_2 = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
						perr_O_2           = [999999.99999,999999.99999,999999.99999,999999.99999]
						CTRE_G_O_2         = 999999.99999
						AMPL_G_O_2         = 999999.99999
						SGMA_G_O_2         = 999999.99999
						OFST_G_O_2         = 999999.99999
						FWHM_G_O_2         = 999999.99999
						EW_O_2             = 999999.99999
						EWE_O_2            = 999999.99999

						CTRE_G_O_2_E       = 999999.99999
						AMPL_G_O_2_E       = 999999.99999
						SGMA_G_O_2_E       = 999999.99999
						CTRE_G_O_2_cor     = 999999.99999
						AMPL_G_O_2_cor     = 999999.99999
						SGMA_G_O_2_cor     = 999999.99999
						OFST_G_O_2_cor     = 999999.99999
						chisqr_O_2         = 999999.99999
						redchi_O_2         = 999999.99999

						AMPL_SNR_2         = 999999.99999
						CTRE_SNR_2         = 999999.99999
						SGMA_SNR_2         = 999999.99999						
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines-1])+'-'+str(LINES[3][lines-1]),'yellow')
						print
						print colored(str(LINES[5][lines-1])+'_CGLC: ' + str(CTRE_G_C_1)  ,'yellow')
						print colored(str(LINES[5][lines-1])+'_AGLC: ' + str(AMPL_G_C_1)  ,'yellow')
						print colored(str(LINES[5][lines-1])+'_SGLC: ' + str(SGMA_G_C_1)  ,'yellow')
						print colored(str(LINES[5][lines-1])+'_FGLC: ' + str(FWHM_G_C_1)  ,'yellow')
						print colored(str(LINES[5][lines-1])+'_WGLC: ' + str(EW_C_1)      ,'yellow')
						print colored(str(LINES[5][lines-1])+'_EGLC: ' + str(EWE_C_1),'yellow')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O

						if fit_vls_hdr == True and fix_ctr_gau_2==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLO2',float(CTRE_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALO2',float(AMPL_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLO2',float(FWHM_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLO2',float(EW_O_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELO2',float(EWE_O_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OFO2',float(OFST_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLC2',float(CTRE_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALC2',float(AMPL_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLC2',float(SGMA_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLC2',float(FWHM_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLC2',float(EW_C_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELC2',float(EWE_C_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CEC2',float(CTRE_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AEC2',float(AMPL_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SEC2',float(SGMA_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHL2',float(chisqr_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 2-1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRL2',float(redchi_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines-1]) + ' Chi2 Reduced 2-1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('*****************************CLEANING*****************************','yellow')
							print colored('The fit (CTR-C & CTR_O) values will be added to the fits headers!','magenta')
							print colored('*****************************CLEANING*****************************','yellow')
							print
						else:
							print
							print colored('The fit (CTR-C & CTR_O) values will NOT be added to the fits headers!','magenta')
							print

						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values G-2 for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF02',float(L1_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-1-2 Wdt-Fit  1GF-IntVal')         #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WP02',float(L2_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-2-2 Wdt-Plt  1GF-IntVal')         #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF02',float(L7_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-7-2 Ctr Fit Bnds  1GF-IntVal')    #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CO02',float(L8_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-8-2 Ctr Fit Ofst  1GF-IntVal')    #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF02',float(L10_2),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10-2 Amp Fit Bnds  1GF-IntVal')   #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial Values: ','cyan')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')					
							print colored(str(LINES[5][lines])+'_WF02' + ': ' + str(L1_2),'cyan')
							print colored(str(LINES[5][lines])+'_WP02' + ': ' + str(L2_2),'cyan')
							print colored(str(LINES[5][lines])+'_CF02' + ': ' + str(L7_2),'cyan')
							print colored(str(LINES[5][lines])+'_CO02' + ': ' + str(L8_2),'cyan')
							print colored(str(LINES[5][lines])+'_AF02' + ': ' + str(L10_2),'cyan')
							print
						else:
							print
							print colored('Initial Guess Values G-2 for line Fitting will NOT be recorded!','yellow')
							print
							pass							
						#################################################CENTRAL GAUSSIAN-2-C##################################################					
					elif fix_ctr_gau_2 == True:
						print
						print colored('2nd CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						try:
							CTRE_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF02')
							AMPL_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF02')
							FWHM_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF02')
							EW_0_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF02')
							EWE_0_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF02')

							CTRE_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLO2')
							AMPL_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALO2')
							FWHM_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLO2')
							EW_O_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLO2')
							EWE_O_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELO2')
							OFST_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFO2')

							CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC2')
							AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC2')
							SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC2')
							FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC2')
							EW_C_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLC2')
							EWE_C_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELC2')
							CTRE_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CEC2')
							AMPL_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AEC2')
							SGMA_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SEC2')

							chisqr_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHL2')
							redchi_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRL2')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines-1])+'_CF01','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau_2=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CGF0')
							AMPL_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_AGF0')
							FWHM_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_FGF0')
							EW_0_2        = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_WGF0')
							EWE_0_2       = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_EGF0')

							CTRE_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CGLO')
							AMPL_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_AGLO')
							FWHM_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_FGLO')
							EW_O_2        = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_WGLO')
							EWE_O_2       = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_EGLO')
							OFST_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_OFSO')

							CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CGLC')
							AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_AGLC')
							SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_SGLC')
							FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_FGLC')
							EW_C_2        = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_WGLC')
							EWE_C_2       = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_EGLC')
							CTRE_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CLEC')
							AMPL_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_ALEC')
							SGMA_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_SLEC')

							chisqr_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CHGL')
							redchi_C_2    = Header_Get(specfile_glx,str(LINES[5][lines-1])+'_CRGL')
					print
					print colored(str(LINES[0][lines-1])+'-'+str(LINES[3][lines-1])+'-CTR','magenta')
					print
					print colored(str(LINES[5][lines-1])+'_CLC2: ' + str(CTRE_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines-1])+'_ALC2: ' + str(AMPL_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines-1])+'_SLC2: ' + str(SGMA_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines-1])+'_FLC2: ' + str(FWHM_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines-1])+'_WLC2: ' + str(EW_C_2)      ,'magenta')
					print colored(str(LINES[5][lines-1])+'_ELC2: ' + str(EWE_C_2),'magenta')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','magenta')
					print colored(str(CTRE_G_C_2)+', '+str(AMPL_G_C_2)+', '+str(SGMA_G_C_2),'magenta')
					print
					##################################################CENTRAL GAUSSIAN-2###################################################
					###############################################COMPUTING TOTAL AREA###############################################
					print colored('Computing Flux Area','yellow')
					###############################################COMPUTING TOTAL AREA###############################################
					CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC1')
					AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC1')
					SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC1')

					CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC2')
					AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC2')
					SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC2')

					print
					print colored('Computing Areas using info from fits headers.','yellow')
					print 'AMPL_G_C_1: ',str(AMPL_G_C_1),'SGMA_G_C_1: ',str(SGMA_G_C_1)
					print 'AMPL_G_C_2: ',str(AMPL_G_C_2),'SGMA_G_C_2: ',str(SGMA_G_C_2)					
					print

					W_C_1     = integrate.quad(lambda x: AMPL_G_C_1*np.exp(-((x)**2)/(2*SGMA_G_C_1**2)), -np.inf, np.inf)
					EW_C_1    = np.round(abs(np.asarray(W_C_1[0])),10)
					EWE_C_1   = np.round(abs(np.asarray(W_C_1[1])),10)

					W_C_2     = integrate.quad(lambda x: AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2)), -np.inf, np.inf)
					EW_C_2    = np.round(abs(np.asarray(W_C_2[0])),10)
					EWE_C_2   = np.round(abs(np.asarray(W_C_2[1])),10)

					W_C       = integrate.quad(lambda x:  AMPL_G_C_1*np.exp(-((x)**2)/(2*SGMA_G_C_1**2)) + 
												AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2)), 
												-np.inf, np.inf
												)
					EW_C      = np.round((np.asarray(W_C[0])),10)
					EWE_C     = np.round(abs(np.asarray(W_C[1])),10)

					if EW_C_1 == 999999.99999:
						EW_C_1 = 0
					else:
						pass
					if EW_C_2 == 999999.99999:
						EW_C_2 = 0
					else:
						pass
						
					print
					print colored('Areas     :','yellow')
					print colored('Area CTR-1: ' + str(EW_C_1),'yellow')
					print colored('Area CTR-2: ' + str(EW_C_2),'yellow')					
					print colored('Area CTR-B: ' + str(EW_C),'yellow')
					print					
					###############################################COMPUTING TOTAL AREA###############################################
					#############################################ADDING AREA TO FTIS HEADER#############################################
					if fit_vls_hdr == True:
						print
						print colored('The Areas values will be updated to the fits headers!','magenta')
						print
						print colored('Area CTR-1: '                  + str(EW_C_1)   + '-' +str(LINES[5][lines])+'_WMC1','yellow')
						print colored('Area CTR-2: '                  + str(EW_C_2)   + '-' +str(LINES[5][lines])+'_WMC2','yellow')
						print colored('Area CTR-B: '                  + str(EW_C)     + '-' +str(LINES[5][lines])+'_WGM1','yellow')
						print						
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC1',float(EW_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C1-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC1',float(EWE_C_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C1-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC2',float(EW_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C2-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC2',float(EWE_C_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C2-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM1',float(EW_C)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C1-C2 Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM1',float(EWE_C)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C1-C2 Crct' + str(fit_type))
					else:
						print
						print colored('The Areas values will NOT be updated to the fits headers!','magenta')
						print
					#############################################ADDING AREA TO FTIS HEADER#############################################
				elif 'Dbl' in LINES[3][lines] and fit_fnct=='gaussM' and fit_type == 'lmfit' and uft_lne_vls == True:
					print 'Line-fitting Cleaning Method. To be checked Fnc_Stk_Plt.py def(Plot_Idp_Spc_Lne) line 11920!'
					fit_typ = 'GM'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					print 
					print colored('Double Line Fit (Ind)','yellow')
					print colored(LINES[3][lines-2]  + '-' + str(LINES[0][lines-2])   + '-' + str(LINES[1][lines-2]),'cyan')
					print LINES[3][lines]+ '-' + str(LINES[0][lines]) + '-'  + str(LINES[1][lines])
					print colored(LINES[3][lines-1]+ '-' + str(LINES[0][lines-1]) + '-' + str(LINES[1][lines-1]),'magenta')

					from lmfit import Model

					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
					if ivl_fts_hdr == True:
						try:
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial fit variables from fits header!','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WM01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WM01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WM02' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WM02' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print '*****'
							print
							print
							print colored('Initial fit variables from fits header!','yellow')
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print							
							#quit()
						try:
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')							
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')							
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AM01',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AM02',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('*****Success!******','yellow')
							print
						################FOR CASES WHERE KLINE WIIDHT ==0################
						if L1_1  == 0:
							L1_1  = 1#LINES[1][lines]
						else:
							pass
						if L1_2  == 0:
							L1_2  = 1#LINES[1][lines]
						else:
							pass
						################FOR CASES WHERE KLINE WIIDHT ==0################
					elif ivl_fts_hdr == False:
						L1_1  = LINES[1][lines-2]
						L2_1  = LINES[2][lines-2]
						L7_1  = LINES[7][lines-2]
						L8_1  = LINES[8][lines-2]
						L10_1 = LINES[10][lines-2]

						L1_2  = LINES[1][lines-1]
						L2_2  = LINES[2][lines-1]
						L7_2  = LINES[7][lines-1]
						L8_2  = LINES[8][lines-1]
						L10_2 = LINES[10][lines-1]
						print
						print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
						print					
					print colored('Initial Values: ','cyan')
					print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')					
					print colored(str(LINES[5][lines-2])+'_WF01' + ': ' + str(L1_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_WP01' + ': ' + str(L2_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_CF01' + ': ' + str(L7_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_CO01' + ': ' + str(L8_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_AF01' + ': ' + str(L10_1),'cyan')
					print
					print colored('Initial Values: ','magenta')
					print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')					
					print colored(str(LINES[5][lines-1])+'_WF02' + ': ' + str(L1_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_WP02' + ': ' + str(L2_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_CF02' + ': ' + str(L7_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_CO02' + ': ' + str(L8_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_AF02' + ': ' + str(L10_2),'magenta')
					print
					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
					lmb_min_lim_line_ft = (LINES[0][lines-2]+L8_1) - MSK_NTMS*L1_1 
					lmb_max_lim_line_ft = (LINES[0][lines-2]+L8_1) + MSK_NTMS*L1_1
					lmb_min_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_1#- 20#L2_1 - 10 
					lmb_max_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_1#+ 20#L2_1 + 10

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function

					idx_ctr_ft_reg = np.abs(lambda_glx[mask_ft] - (LINES[0][lines-2]+L8_1)).argmin()

					X0_f2DG    = (LINES[0][lines-2]+L8_1)
					SIGMA_f2DG = L1_1
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))
					A_f2DG     = -(1-inten_glx[mask_ft][idx_ctr_ft_reg])

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
				
					##################################################CENTRAL GAUSSIAN-1###################################################
					if fix_ctr_gau_1 == False:
						print
						print colored('Fitting 1st line','cyan')
						print colored('1-0-Fitting Central line','cyan')
						print
						#################################################CENTRAL GAUSSIAN-1-0##################################################
						popt_0_1, pcov_0_1   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_0_1           = [999999.99999,999999.99999,999999.99999]
						CTRE_G_0_1         = 999999.99999
						AMPL_G_0_1         = 999999.99999
						SGMA_G_0_1         = 999999.99999
						FWHM_G_0_1         = 999999.99999
						EW_0_1             = 999999.99999
						EWE_0_1            = 999999.99999

						CTRE_G_0_1_E      = 999999.99999
						AMPL_G_0_1_E      = 999999.99999
						SGMA_G_0_1_E      = 999999.99999

						CTRE_G_0_1_cor    = 999999.99999
						AMPL_G_0_1_cor    = 999999.99999
						SGMA_G_0_1_cor    = 999999.99999

						chisqr_0        = 999999.99999						
						if fit_vls_hdr == True and fix_ctr_gau_1 ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM01',float(CTRE_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM01',float(AMPL_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FM01',float(FWHM_G_0_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM01',float(EW_0_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EM01',float(EWE_0_1)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct-1' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('The fit (CTR-1-0) values will be added to the fits headers!','cyan')
							print 'aaaaa'
							print
						else:
							print
							print colored('The fit (CTR-1-0) values will NOT be added to the fits headers!','cyan')
							print
						#################################################CENTRAL GAUSSIAN-1-0##################################################
						#################################################CENTRAL GAUSSIAN-1-C##################################################
						popt_C_1, pcov_C_1 = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C_1           = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C_1         = 999999.99999
						AMPL_G_C_1         = 999999.99999
						SGMA_G_C_1         = 999999.99999
						FWHM_G_C_1         = 999999.99999
						EW_C_1             = 999999.99999
						EWE_C_1            = 999999.99999

						CTRE_G_C_1_E       = 999999.99999
						AMPL_G_C_1_E       = 999999.99999
						SGMA_G_C_1_E       = 999999.99999
						CTRE_G_C_1_cor     = 999999.99999
						AMPL_G_C_1_cor     = 999999.99999
						SGMA_G_C_1_cor     = 999999.99999
						chisqr_C_1         = 999999.99999
						redchi_C_1         = 999999.99999

						popt_O_1 ,pcov_O_1 = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
						perr_O_1           = [999999.99999,999999.99999,999999.99999,999999.99999]
						CTRE_G_O_1         = 999999.99999
						AMPL_G_O_1         = 999999.99999
						SGMA_G_O_1         = 999999.99999
						OFST_G_O_1         = 999999.99999
						FWHM_G_O_1         = 999999.99999
						EW_O_1             = 999999.99999
						EWE_O_1            = 999999.99999

						CTRE_G_O_1_E      = 999999.99999
						AMPL_G_O_1_E      = 999999.99999
						SGMA_G_O_1_E      = 999999.99999
						CTRE_G_O_1_cor    = 999999.99999
						AMPL_G_O_1_cor    = 999999.99999
						SGMA_G_O_1_cor    = 999999.99999
						OFST_G_O_1_cor    = 999999.99999
						chisqr_O_1        = 999999.99999
						redchi_O_1        = 999999.99999

						AMPL_SNR_1        = 999999.99999
						CTRE_SNR_1        = 999999.99999
						SGMA_SNR_1        = 999999.99999						
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines-2])+'-'+str(LINES[3][lines-2]),'cyan')
						print
						print colored(str(LINES[5][lines-2])+'_CGLC: ' + str(CTRE_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_AGLC: ' + str(AMPL_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_SGLC: ' + str(SGMA_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_FGLC: ' + str(FWHM_G_C_1)  ,'cyan')
						print colored(str(LINES[5][lines-2])+'_WGLC: ' + str(EW_C_1)      ,'cyan')
						print colored(str(LINES[5][lines-2])+'_EGLC: ' + str(EWE_C_1)     ,'cyan')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O
						if fit_vls_hdr == True and fix_ctr_gau_1==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                         ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                         ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CMO1',float(CTRE_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AMO1',float(AMPL_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FMO1',float(FWHM_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WMO1',float(EW_O_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EMO1',float(EWE_O_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OMO1',float(OFST_G_O_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CMC1',float(CTRE_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AMC1',float(AMPL_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SMC1',float(SGMA_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FMC1',float(FWHM_G_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC1',float(EW_C_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC1',float(EWE_C_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM1E',float(CTRE_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM1E',float(AMPL_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SM1E',float(SGMA_G_C_1_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHM1',float(chisqr_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRM1',float(redchi_C_1)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines-2]) + ' Chi2 Reduced 1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('The fit (CTR-1-C & CTR-1-O) values will be added to the fits headers!','cyan')
							print
						else:
							print
							print colored('The fit (CTR-1-C & CTR-1-O) values will NOT be added to the fits headers!','cyan')
							print

						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values G-1 for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM01',float(L1_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-1-1 Wdt-Fit  1GF-IntVal')         #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM01',float(L2_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-2-1 Wdt-Plt  1GF-IntVal')         #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM01',float(L7_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-7-1 Ctr Fit Bnds  1GF-IntVal')    #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM01',float(L8_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-8-1 Ctr Fit Ofst  1GF-IntVal')    #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM01',float(L10_1),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10-1 Amp Fit Bnds  1GF-IntVal')   #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial Values: ','cyan')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')					
							print colored(str(LINES[5][lines])+'_WM01' + ': ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WM01' + ': ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ': ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ': ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AM01' + ': ' + str(L10_1),'cyan')
							print
						else:
							print
							print colored('Initial Guess Values G-1 for line Fitting will NOT be recorded!','yellow')
							print
							pass
					#################################################CENTRAL GAUSSIAN-1-C##################################################					
					elif fix_ctr_gau_1 == True:
						print
						print colored('1-CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						try:
							CTRE_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CM01')
							AMPL_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AM01')
							FWHM_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FM01')
							EW_0_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WM01')
							EWE_0_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EM01')

							CTRE_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMO1')
							AMPL_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMO1')
							FWHM_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMO1')
							EW_O_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WMO1')
							EWE_O_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EMO1')
							OFST_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OMO1')

							CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC1')
							AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC1')
							SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC1')
							FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMC1')
							EW_C_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WMC1')
							EWE_C_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EMC1')
							CTRE_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC1')
							AMPL_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC1')
							SGMA_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC1')

							chisqr_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHM1')
							redchi_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRM1')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines-2])+'_CF01','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau_1=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF01')
							AMPL_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF01')
							FWHM_G_0_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF01')
							EW_0_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF01')
							EWE_0_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF01')

							CTRE_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLO1')
							AMPL_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALO1')
							FWHM_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLO1')
							EW_O_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLO1')
							EWE_O_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELO1')
							OFST_G_O_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFO1')

							CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC1')
							AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC1')
							SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC1')
							FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC1')
							EW_C_1        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLC1')
							EWE_C_1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELC1')
							CTRE_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CEC1')
							AMPL_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AEC1')
							SGMA_G_C_1_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SEC1')

							chisqr_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHL1')
							redchi_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRL1')
					print
					print colored(str(LINES[0][lines-2])+'-'+str(LINES[3][lines-2])+'-CTR','cyan')
					print colored('From: '+str(LINES[3][lines])+'-CTR','cyan')
					print
					print colored(str(LINES[5][lines])+'_CMC1: ' + str(CTRE_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines])+'_AMC1: ' + str(AMPL_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines])+'_SMC1: ' + str(SGMA_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines])+'_FMC1: ' + str(FWHM_G_C_1)  ,'cyan')
					print colored(str(LINES[5][lines])+'_WMC1: ' + str(EW_C_1)      ,'cyan')
					print colored(str(LINES[5][lines])+'_EMC1: ' + str(EWE_C_1),'cyan')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','cyan')
					print colored(str(CTRE_G_C_1)+', '+str(AMPL_G_C_1)+', '+str(SGMA_G_C_1),'cyan')
					print
					##################################################CENTRAL GAUSSIAN-1###################################################
					lmb_min_lim_line_ft = (LINES[0][lines-1]+L8_2) - MSK_NTMS*LINES[1][lines-1] 
					lmb_max_lim_line_ft = (LINES[0][lines-1]+L8_2) + MSK_NTMS*LINES[1][lines-1]
					lmb_min_lim_line    = (LINES[0][lines-1]+L8_2)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_2
					lmb_max_lim_line    = (LINES[0][lines-1]+L8_2)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_2

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] 
					idx_ctr_ft_reg = np.abs(lambda_glx[mask_ft] - (LINES[0][lines-1]+L8_2)).argmin()

					X0_f2DG    = (LINES[0][lines-1]+L8_2)
					SIGMA_f2DG = LINES[1][lines-1]
					A_f2DG     = -(1-(min(inten_glx[mask_ft])))
					A_f2DG     = -(1-inten_glx[mask_ft][idx_ctr_ft_reg])

					max_val    = -(1-min(inten_glx[mask_ft]))
					lambda_max = (lambda_glx[np.where(inten_glx==min(inten_glx[mask_ft]))[0]][0])					

					if pre_off_plt == True:
						plt.step(lambda_glx[mask_pl], inten_glx[mask_pl],
								where='mid',lw=3.0,alpha=0.5,linestyle=':',color='gray',
								label='Original Spectrum')
					elif pre_off_plt == False:
						pass

					initial_guess_O   = (X0_f2DG,A_f2DG,SIGMA_f2DG,max(inten_glx[mask_ft])-1)
					initial_guess_0   = (X0_f2DG,A_f2DG,SIGMA_f2DG)
					##################################################CENTRAL GAUSSIAN-2###################################################
					if fix_ctr_gau_2 == False:
						print
						print colored('Fitting 2nd line','magenta')
						print colored('2-0-Fitting Central line','magenta')
						print
						#################################################CENTRAL GAUSSIAN-2-0##################################################

						popt_0, pcov_0   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_0           = [999999.99999,999999.99999,999999.99999]
						CTRE_G_0_2         = 999999.99999
						AMPL_G_0_2         = 999999.99999
						SGMA_G_0_2         = 999999.99999
						FWHM_G_0         = 999999.99999
						EW_0_2             = 999999.99999
						EWE_0_2            = 999999.99999

						CTRE_G_0_2_E      = 999999.99999
						AMPL_G_0_2_E      = 999999.99999
						SGMA_G_0_2_E      = 999999.99999

						CTRE_G_0_2_cor    = 999999.99999
						AMPL_G_0_2_cor    = 999999.99999
						SGMA_G_0_2_cor    = 999999.99999

						chisqr_0        = 999999.99999
						redchi_0        = 999999.99999						
						if fit_vls_hdr == True and fix_ctr_gau_2 ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM02',float(CTRE_G_0_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM02',float(AMPL_G_0_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FM02',float(FWHM_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM02',float(EW_0_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EM02',float(EWE_0_2)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct-2' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('The fit (CTR-2-0) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-2-0) values will NOT be added to the fits headers!','magenta')
							print
						#################################################CENTRAL GAUSSIAN-2-0##################################################
						#################################################CENTRAL GAUSSIAN-2-C##################################################

						print colored('Cleaning Method','cyan')
						popt_C_2, pcov_C_2 = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C_2           = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C_2         = 999999.99999
						AMPL_G_C_2         = 999999.99999
						SGMA_G_C_2         = 999999.99999
						FWHM_G_C_2         = 999999.99999
						EW_C_2             = 999999.99999
						EWE_C_2            = 999999.99999

						CTRE_G_C_2_E       = 999999.99999
						AMPL_G_C_2_E       = 999999.99999
						SGMA_G_C_2_E       = 999999.99999
						CTRE_G_C_2_cor     = 999999.99999
						AMPL_G_C_2_cor     = 999999.99999
						SGMA_G_C_2_cor     = 999999.99999
						chisqr_C_2         = 999999.99999
						redchi_C_2         = 999999.99999

						popt_O_2 ,pcov_O_2 = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
						perr_O_2           = [999999.99999,999999.99999,999999.99999,999999.99999]
						CTRE_G_O_2         = 999999.99999
						AMPL_G_O_2         = 999999.99999
						SGMA_G_O_2         = 999999.99999
						OFST_G_O_2         = 999999.99999
						FWHM_G_O_2         = 999999.99999
						EW_O_2             = 999999.99999
						EWE_O_2            = 999999.99999

						CTRE_G_O_2_E       = 999999.99999
						AMPL_G_O_2_E       = 999999.99999
						SGMA_G_O_2_E       = 999999.99999
						CTRE_G_O_2_cor     = 999999.99999
						AMPL_G_O_2_cor     = 999999.99999
						SGMA_G_O_2_cor     = 999999.99999
						OFST_G_O_2_cor     = 999999.99999
						chisqr_O_2         = 999999.99999
						redchi_O_2         = 999999.99999

						AMPL_SNR_2         = 999999.99999
						CTRE_SNR_2         = 999999.99999
						SGMA_SNR_2         = 999999.99999						
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines-1])+'-'+str(LINES[3][lines-1]),'yellow')
						print
						print colored(str(LINES[5][lines])+'_CGLC: ' + str(CTRE_G_C_2)  ,'yellow')
						print colored(str(LINES[5][lines])+'_AGLC: ' + str(AMPL_G_C_2)  ,'yellow')
						print colored(str(LINES[5][lines])+'_SGLC: ' + str(SGMA_G_C_2)  ,'yellow')
						print colored(str(LINES[5][lines])+'_FGLC: ' + str(FWHM_G_C_2)  ,'yellow')
						print colored(str(LINES[5][lines])+'_WGLC: ' + str(EW_C_2)      ,'yellow')
						print colored(str(LINES[5][lines])+'_EGLC: ' + str(EWE_C_2),'yellow')
						print
						print
						#inten_glx[mask_ft] = inten_glx[mask_ft] + OFST_G_O

						if fit_vls_hdr == True and fix_ctr_gau_2==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                         ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                         ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CMO2',float(CTRE_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AMO2',float(AMPL_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FMO2',float(FWHM_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WMO2',float(EW_O_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EMO2',float(EWE_O_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OMO2',float(OFST_G_O_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 2-1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CMC2',float(CTRE_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AMC2',float(AMPL_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SMC2',float(SGMA_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FMC2',float(FWHM_G_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC2',float(EW_C_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC2',float(EWE_C_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM2E',float(CTRE_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM2E',float(AMPL_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SM2E',float(SGMA_G_C_2_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 2-1GF Crct' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHM2',float(chisqr_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 2-1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRM2',float(redchi_C_2)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines-1]) + ' Chi2 Reduced 2-1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('The fit (CTR-2-C & CTR-2-O) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-2-C & CTR-2-O) values will NOT be added to the fits headers!','magenta')
							print

						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values G-2 for line Fitting will be recorded!','magenta')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','magenta')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM02',float(L1_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-1-2 Wdt-Fit  1GF-IntVal')         #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WM02',float(L2_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-2-2 Wdt-Plt  1GF-IntVal')         #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM02',float(L7_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-7-2 Ctr Fit Bnds  1GF-IntVal')    #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CM02',float(L8_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-8-2 Ctr Fit Ofst  1GF-IntVal')    #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AM02',float(L10_2),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10-2 Amp Fit Bnds  1GF-IntVal')   #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial Values: ','magenta')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')					
							print colored(str(LINES[5][lines])+'_WM02' + ': ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WM02' + ': ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ': ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ': ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AM02' + ': ' + str(L10_2),'magenta')
							print
						else:
							print
							print colored('Initial Guess Values G-2 for line Fitting will NOT be recorded!','yellow')
							print
							pass							
						#################################################CENTRAL GAUSSIAN-2-C##################################################					
					elif fix_ctr_gau_2 == True:
						print
						print colored('2-CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						print colored(LINES[3][lines],'yellow')
						print
						try:
							CTRE_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CM02')
							AMPL_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AM02')
							FWHM_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FM02')
							EW_0_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WM02')
							EWE_0_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EM02')

							CTRE_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMO2')
							AMPL_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMO2')
							FWHM_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMO2')
							EW_O_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WMO2')
							EWE_O_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EMO2')
							OFST_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OMO2')

							CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC2')
							AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC2')
							SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC2')
							FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMC2')
							EW_C_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WMC2')
							EWE_C_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EMC2')
							CTRE_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC2')
							AMPL_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC2')
							SGMA_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC2')

							chisqr_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHM2')
							redchi_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRM2')

						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines-1])+'_CF01','yellow')
							print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
							print colored ('Or UnFix (fix_ctr_gau_2=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							CTRE_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF02')
							AMPL_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF02')
							FWHM_G_0_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF02')
							EW_0_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF02')
							EWE_0_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF02')

							CTRE_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLO2')
							AMPL_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALO2')
							FWHM_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLO2')
							EW_O_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLO2')
							EWE_O_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELO2')
							OFST_G_O_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFO2')

							CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC2')
							AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC2')
							SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC2')
							FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC2')
							EW_C_2        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLC2')
							EWE_C_2       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELC2')
							CTRE_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CEC2')
							AMPL_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AEC2')
							SGMA_G_C_2_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SEC2')

							chisqr_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHL2')
							redchi_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRL2')
					print
					print colored(str(LINES[0][lines-1])+'-'+str(LINES[3][lines-1])+'-CTR','magenta')
					print colored('From '+str(LINES[3][lines])+'-CTR','cyan')
					print
					print colored(str(LINES[5][lines])+'_CMC2: ' + str(CTRE_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines])+'_AMC2: ' + str(AMPL_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines])+'_SMC2: ' + str(SGMA_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines])+'_FMC2: ' + str(FWHM_G_C_2)  ,'magenta')
					print colored(str(LINES[5][lines])+'_WMC2: ' + str(EW_C_2)      ,'magenta')
					print colored(str(LINES[5][lines])+'_EMC2: ' + str(EWE_C_2),'magenta')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','magenta')
					print colored(str(CTRE_G_C_2)+', '+str(AMPL_G_C_2)+', '+str(SGMA_G_C_2),'magenta')
					print
					##################################################CENTRAL GAUSSIAN-2###################################################	
					CTRE_G_0 = CTRE_G_0_1
					#####################################################PRE GAUSSIAN#################################################
					if fix_pre_gau == False and pst_shf_lim>0:
						###########################################DEFINING PRE-PST REGIONS##################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
							lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						else:
							print
							print colored('Using Expected Line Center to Define Fitting Region!','yellow')
							print
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
							
						print
						print colored('Region limits for fitting.','yellow')
						print colored('Line Center    : ' + str(LINES[0][lines]),'cyan')
						print colored('pre_shf_ctr    : ' + str(pre_shf_ctr),'cyan')
						print colored('pre_shf_lim    : ' + str(pre_shf_lim),'cyan')
						print colored('Line Center    : ' + str(LINES[0][lines] - pre_shf_ctr),'cyan')
						print colored('Lower Limit    : ' + str(LINES[0][lines] - (pre_shf_ctr-pre_shf_lim)),'cyan')
						print colored('Upper Limit    : ' + str(LINES[0][lines] - (pre_shf_ctr+pre_shf_lim)),'cyan')
						print 
						print colored('Central    : ' + str(LINES[0][lines]),'cyan')
						print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'cyan')
						print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'cyan')
						print
						print colored('Limits:','cyan')
						print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
						print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr - pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr + pre_shf_lim
						print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr - pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr + pst_shf_lim
						print
						print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr - pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr + pst_shf_lim
						print
						###########################################DEFINING PRE-PST REGIONS##################################################
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						x_a = lambda_glx[mask_ft_pre][X0_f2DG_indx_PRE]
						y_a = inten_glx[mask_ft_pre][X0_f2DG_indx_PRE]

						print colored('Cleaning Method','cyan')
						popt_C_PRE, pcov_C_PRE  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C_PRE          = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C_PRE        = 999999.99999
						AMPL_G_C_PRE        = 999999.99999
						SGMA_G_C_PRE        = 999999.99999
						FWHM_G_C_PRE        = 999999.99999
						EW_C_PRE            = 999999.99999
						EWE_C_PRE           = 999999.99999
						
						CTRE_G_C_PRE_E      = 999999.99999
						AMPL_G_C_PRE_E      = 999999.99999
						SGMA_G_C_PRE_E      = 999999.99999
						CTRE_G_C_PRE_cor    = 999999.99999
						AMPL_G_C_PRE_cor    = 999999.99999
						SGMA_G_C_PRE_cor    = 999999.99999
						chisqr_C_PRE        = 999999.99999
						redchi_C_PRE        = 999999.99999						
							
						if fit_vls_hdr == True and fix_pre_gau == False:						
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHML',float(chisqr_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRML',float(redchi_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 Reduced 1GF' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGM1',float(CTRE_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGM1',float(AMPL_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGM1',float(SGMA_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGM1',float(FWHM_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM1',float(EW_C_PRE)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM1',float(EWE_C_PRE)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CME1',float(CTRE_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AME1',float(AMPL_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SME1',float(SGMA_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							
							Header_Add(specfile_glx,str(LINES[5][lines])+'_XAM1',float(x_a),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PRE GAU-LNR X1 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YAM1',float(y_a),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PRE GAU-LNR Y1 COO')
							print
							print colored('The fit (PRE) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (PRE) values will NOT be added to the fits headers!','magenta')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GSM1',float(pre_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for PRE G-1')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GCM1',float(pre_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for PRE G-1')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass
					elif fix_pre_gau == True or (fix_pre_gau == False and pre_shf_lim<=0):
						try:
							print
							print colored('1 PRE-gaussian already fitted!','yellow')
							print colored('Values from fits header','yellow')
							chisqr_C       = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHML')
							redchi_C       = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRML')
							CTRE_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM1')
							AMPL_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM1')
							SGMA_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM1')
							FWHM_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM1')
							EW_C_PRE       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGM1')
							EWE_C_PRE      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGM1')
							CTRE_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CME1')
							AMPL_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_AME1')
							SGMA_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SME1')
							
							x_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XAM1')
							y_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YAM1')
						
							pre_shf_lim    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GSM1')
							pre_shf_ctr    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GCM1')
							
						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!' ,'yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
							print colored ('Gotta Fit first before fixing a component (PRE)!','yellow')
							print colored ('Or UnFix (fix_pre_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							quit()
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-PRE','cyan')
					print
					print '******************************************************************************'
					print colored('Boundaries for Gaussian Fitting (PRE-CTR): ' + str(pre_shf_ctr),'cyan')
					print colored('Boundaries for Gaussian Fitting (PRE-LIM): ' + str(pre_shf_lim),'cyan')
					print '******************************************************************************'					
					print
					print colored(str(LINES[5][lines])+'_CGM1: ' + str(CTRE_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_AGM1: ' + str(AMPL_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_SGM1: ' + str(SGMA_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_FGM1: ' + str(FWHM_G_C_PRE)  ,'cyan')
					print colored(str(LINES[5][lines])+'_WGM1: ' + str(EW_C_PRE)      ,'cyan')
					print colored(str(LINES[5][lines])+'_EGM1: ' + str(EWE_C_PRE),'cyan')
					print
					print colored(str(LINES[5][lines])+'_XAM1 : ' + str(x_a),'cyan')
					print colored(str(LINES[5][lines])+'_YAM1 : ' + str(y_a),'cyan')
					print colored(str(LINES[5][lines])+'_GSM1 : ' + str(pre_shf_lim),'cyan')
					print colored(str(LINES[5][lines])+'_GCM1 : ' + str(pre_shf_ctr),'cyan')
					print
					print colored('Fit Values (PRE) Center, Amplitude, Sigma ('+fit_type+'):','cyan')
					print colored(str(CTRE_G_C_PRE)+', '+str(AMPL_G_C_PRE)+', '+str(SGMA_G_C_PRE),'cyan')
					print
					#####################################################PRE GAUSSIAN#################################################
					###################################################POST GAUSSIAN##################################################
					if fix_pst_gau == False and pst_shf_lim>0:
						###########################################DEFINING PRE-PST REGIONS##################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
							lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						else:
							print
							print colored('Using Expected Line Center to Define Fitting Region!','yellow')
							print
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
							
						print
						print colored('Region limits for fitting.','yellow')
						print colored('Line Center    : ' + str(LINES[0][lines]),'cyan')
						print colored('pst_shf_ctr    : ' + str(pst_shf_ctr),'cyan')
						print colored('pst_shf_lim    : ' + str(pst_shf_lim),'cyan')
						print colored('Line Center    : ' + str(LINES[0][lines] + pst_shf_ctr),'cyan')
						print colored('Lower Limit    : ' + str(LINES[0][lines] - (pst_shf_ctr-pst_shf_lim)),'cyan')
						print colored('Upper Limit    : ' + str(LINES[0][lines] - (pst_shf_ctr+pst_shf_lim)),'cyan')
						print 
						print colored('Central    : ' + str(LINES[0][lines]),'cyan')
						print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'cyan')
						print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'cyan')
						print
						print colored('Limits:','yellow')
						print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
						print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim
						print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
						print
						print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
						print
						###########################################DEFINING PRE-PST REGIONS##################################################						
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						x_b = lambda_glx[mask_ft_pst][X0_f2DG_indx_PST]
						y_b = inten_glx[mask_ft_pst][X0_f2DG_indx_PST]

						print colored('Cleaning Method','cyan')
						popt_C_PST, pcov_C_PST  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C_PST          = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C_PST        = 999999.99999
						AMPL_G_C_PST        = 999999.99999
						SGMA_G_C_PST        = 999999.99999
						FWHM_G_C_PST        = 999999.99999
						EW_C_PST            = 999999.99999
						EWE_C_PST           = 999999.99999
						
													
						CTRE_G_C_PST_E      = 999999.99999
						AMPL_G_C_PST_E      = 999999.99999
						SGMA_G_C_PST_E      = 999999.99999
						CTRE_G_C_PST_cor    = 999999.99999
						AMPL_G_C_PST_cor    = 999999.99999
						SGMA_G_C_PST_cor    = 999999.99999
						chisqr_C_PST        = 999999.99999
						redchi_C_PST        = 999999.99999						
						if fit_vls_hdr == True and fix_pst_gau == False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGM2',float(CTRE_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGM2',float(AMPL_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGM2',float(SGMA_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGM2',float(FWHM_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM2',float(EW_C_PST)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM2',float(EWE_C_PST)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CME2',float(CTRE_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AME2',float(AMPL_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SME2',float(SGMA_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							
							Header_Add(specfile_glx,str(LINES[5][lines])+'_XAM2',float(x_b),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PST GAU-LNR X2 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YAM2',float(y_b),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PST GAU-LNR Y2 COO')
						else:
							print
							print colored('The fit (PST) values will NOT be added to the fits headers!','magenta')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GSM2',float(pst_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for PST G-2')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GCM2',float(pst_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for PST G-2')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass
						#inten_glx[mask_ft_plp] = inten_glx[mask_ft_plp] + OFST_G_O #OFFSET -
					elif fix_pst_gau == True or (fix_pst_gau == False and pst_shf_lim<=0):
						try:
							print
							print colored('2 PST-gaussian already fitted!','yellow')
							print colored('Values from fits header','yellow')
							CTRE_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM2')
							AMPL_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM2')
							SGMA_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM2')
							FWHM_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM2')
							EW_C_PST       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGM2')
							EWE_C_PST      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGM2')
							CTRE_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CME2')
							AMPL_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_AME2')
							SGMA_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SME2')
							
							x_b            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XAM2')
							y_b            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YAM2')
							
							pst_shf_lim    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GSM2')
							pst_shf_ctr    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GCM2')
						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
							print colored ('Gotta Fit first before fixing a component (PST)!','yellow')
							print colored ('Or UnFix (fix_pst_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							quit()
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-PST','magenta')
					print
					print '******************************************************************************'
					print colored('Boundaries for Gaussian Fitting (PST-CTR): ' + str(pst_shf_ctr),'magenta')					
					print colored('Boundaries for Gaussian Fitting (PST-LIM): ' + str(pst_shf_lim),'magenta')
					print '******************************************************************************'
					print
					print colored(str(LINES[5][lines])+'_CGM2: ' + str(CTRE_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_AGM2: ' + str(AMPL_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_SGM2: ' + str(SGMA_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_FGM2: ' + str(FWHM_G_C_PST)  ,'magenta')
					print colored(str(LINES[5][lines])+'_WGM2: ' + str(EW_C_PST)      ,'magenta')
					print colored(str(LINES[5][lines])+'_EGM2: ' + str(EWE_C_PST),'magenta')
					print
					print colored(str(LINES[5][lines])+'_GSM2 : ' + str(pst_shf_lim),'magenta')
					print colored(str(LINES[5][lines])+'_GCM2 : ' + str(pst_shf_ctr),'magenta')
					print colored(str(LINES[5][lines])+'_XAM2 : ' + str(x_b),'magenta')
					print colored(str(LINES[5][lines])+'_YAM2 : ' + str(y_b),'magenta')
					print
					print colored('Fit Values (PST) Center, Amplitude, Sigma ('+fit_type+'):','magenta')
					print colored(str(CTRE_G_C_PST)+', '+str(AMPL_G_C_PST)+', '+str(SGMA_G_C_PST),'magenta')
					print					
					###################################################POST GAUSSIAN##################################################										

					###################################################MDL GAUSSIAN###################################################
					if fix_mdl_gau == False:
						if mdl_shf_lim>0:
							#########################################DEFINING PRE-PST-MDL REGIONS################################################
							if ofs_ctr_fit == True:
								print
								print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
								print
								lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
								lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
								mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
								mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
								mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
								mask_ft_mdl = (lambda_glx >= LINES[0][lines] + mdl_shf_ctr - mdl_shf_lim) & (lambda_glx <= LINES[0][lines] + mdl_shf_ctr + mdl_shf_lim)
								X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
								X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
								X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
							else:
								print
								print colored('Using Expected Line Center to Define Fitting Region!','yellow')
								print
								mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
								mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
								mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
								mask_ft_mdl = (lambda_glx >= LINES[0][lines] + mdl_shf_ctr - mdl_shf_lim) & (lambda_glx <= LINES[0][lines] + mdl_shf_ctr + mdl_shf_lim)
								X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
								X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
								X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
							print
							print colored('Region limits for fitting.','yellow')
							print colored('Line Center    : ' + str(LINES[0][lines]),'cyan')
							print colored('pst_shf_ctr    : ' + str(pst_shf_ctr),'cyan')
							print colored('pst_shf_lim    : ' + str(pst_shf_lim),'cyan')
							print colored('Line Center    : ' + str(LINES[0][lines] + pst_shf_ctr),'cyan')
							print colored('Lower Limit    : ' + str(LINES[0][lines] - (pst_shf_ctr-pst_shf_lim)),'cyan')
							print colored('Upper Limit    : ' + str(LINES[0][lines] - (pst_shf_ctr+pst_shf_lim)),'cyan')
							print 
							print colored('Central    : ' + str(LINES[0][lines]),'cyan')
							print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'cyan')
							print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'cyan')
							print colored('Central+MDL: ' + str(mdl_shf_ctr) + '-' + str(LINES[0][lines]+mdl_shf_ctr),'cyan')
							print
							print colored('Limits:','yellow')
							print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
							print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim
							print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
							print 'MDL    : lambda_glx >= ',LINES[0][lines] + mdl_shf_ctr-mdl_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + mdl_shf_ctr+mdl_shf_lim
							print
							print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
							print
							#########################################DEFINING PRE-PST-MDL REGIONS################################################
							X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
							x_c = lambda_glx[mask_ft_mdl][X0_f2DG_indx_MDL]
							y_c = inten_glx[mask_ft_mdl][X0_f2DG_indx_MDL]

							print colored('Cleaning Method','green')
							popt_C_MDL, pcov_C_MDL  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_MDL          = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_MDL        = 999999.99999
							AMPL_G_C_MDL        = 999999.99999
							SGMA_G_C_MDL        = 999999.99999
							FWHM_G_C_MDL        = 999999.99999
							EW_C_MDL            = 999999.99999
							EWE_C_MDL           = 999999.99999
							
														
							CTRE_G_C_MDL_E      = 999999.99999
							AMPL_G_C_MDL_E      = 999999.99999
							SGMA_G_C_MDL_E      = 999999.99999
							CTRE_G_C_MDL_cor    = 999999.99999
							AMPL_G_C_MDL_cor    = 999999.99999
							SGMA_G_C_MDL_cor    = 999999.99999
							chisqr_C_MDL        = 999999.99999
							redchi_C_MDL        = 999999.99999							
						else:
							popt_C_MDL, pcov_C_MDL  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
							perr_C_MDL          = [999999.99999,999999.99999,999999.99999]
							CTRE_G_C_MDL        = 999999.99999
							AMPL_G_C_MDL        = 999999.99999
							SGMA_G_C_MDL        = 999999.99999
							FWHM_G_C_MDL        = 999999.99999
							EW_C_MDL            = 999999.99999
							EWE_C_MDL           = 999999.99999
							
														
							CTRE_G_C_MDL_E      = 999999.99999
							AMPL_G_C_MDL_E      = 999999.99999
							SGMA_G_C_MDL_E      = 999999.99999
							CTRE_G_C_MDL_cor    = 999999.99999
							AMPL_G_C_MDL_cor    = 999999.99999
							SGMA_G_C_MDL_cor    = 999999.99999
							chisqr_C_MDL        = 999999.99999
							redchi_C_MDL        = 999999.99999
							
							x_c                 = 999999.99999
							y_c                 = 999999.99999							
						if fit_vls_hdr == True and fix_mdl_gau == False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGM3',float(CTRE_G_C_MDL)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGM3',float(AMPL_G_C_MDL)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGM3',float(SGMA_G_C_MDL)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGM3',float(FWHM_G_C_MDL)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM3',float(EW_C_MDL)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM3',float(EWE_C_MDL)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CME2',float(CTRE_G_C_MDL_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AME2',float(AMPL_G_C_MDL_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SME2',float(SGMA_G_C_MDL_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							
							Header_Add(specfile_glx,str(LINES[5][lines])+'_XAM3',float(x_c),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' MDL GAU-LNR X3 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YAM3',float(y_c),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' MDL GAU-LNR Y3 COO')
						else:
							print
							print colored('The fit (MDL) values will NOT be added to the fits headers!','green')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GSM3',float(mdl_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for MDL G-3')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GCM3',float(mdl_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for MDL G-3')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass						
					elif fix_mdl_gau == True:
						try:
							print
							print colored('3- MDL-gaussian already fitted!','yellow')
							print colored('Values from fits header','yellow')
							CTRE_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM3')
							AMPL_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM3')
							SGMA_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM3')
							FWHM_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM3')
							EW_C_MDL       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGM3')
							EWE_C_MDL      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGM3')
							CTRE_G_C_MDL_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CME2')
							AMPL_G_C_MDL_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_AME2')
							SGMA_G_C_MDL_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SME2')
							
							x_c            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XAM3')
							y_c            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YAM3')
							EW_C_PS3       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WRM3')
							EWE_C_PS3      = Header_Get(specfile_glx,str(LINES[5][lines])+'_ERM3')
							
							mdl_shf_lim    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GSM3')
							mdl_shf_ctr    = Header_Get(specfile_glx,str(LINES[5][lines])+'_GCM3')
						except KeyError:
							print
							print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
							print colored ('File  : ' + specfile_glx,'yellow')
							print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
							print colored ('Gotta Fit first before fixing a component (PST)!','yellow')
							print colored ('Or UnFix (fix_pst_gau=False) For Plotting!','yellow')
							print colored ('Quitting!','yellow')
							print
							quit()
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-MDL','green')
					print
					print '******************************************************************************'
					print colored('Boundaries for Gaussian Fitting (PST-CTR): ' + str(mdl_shf_ctr),'green')					
					print colored('Boundaries for Gaussian Fitting (PST-LIM): ' + str(mdl_shf_lim),'green')
					print '******************************************************************************'
					print
					print colored(str(LINES[5][lines])+'_CGM3: ' + str(CTRE_G_C_MDL)  ,'green')
					print colored(str(LINES[5][lines])+'_AGM3: ' + str(AMPL_G_C_MDL)  ,'green')
					print colored(str(LINES[5][lines])+'_SGM3: ' + str(SGMA_G_C_MDL)  ,'green')
					print colored(str(LINES[5][lines])+'_FGM3: ' + str(FWHM_G_C_MDL)  ,'green')
					print colored(str(LINES[5][lines])+'_WGM3: ' + str(EW_C_MDL)      ,'green')
					print colored(str(LINES[5][lines])+'_EGM3: ' + str(EWE_C_MDL),'green')
					print
					print colored(str(LINES[5][lines])+'_GSM3 : ' + str(mdl_shf_lim),'green')
					print colored(str(LINES[5][lines])+'_GCM3 : ' + str(mdl_shf_ctr),'green')
					print colored(str(LINES[5][lines])+'_XAM3 : ' + str(x_c),'green')
					print colored(str(LINES[5][lines])+'_YAM3 : ' + str(y_c),'green')
					print
					print colored('Fit Values (PST) Center, Amplitude, Sigma ('+fit_type+'):','green')
					print colored(str(CTRE_G_C_MDL)+', '+str(AMPL_G_C_MDL)+', '+str(SGMA_G_C_MDL),'green')
					print					
					###################################################MDL GAUSSIAN##################################################						
					###############################################COMPUTING TOTAL AREA###############################################
					print colored('Computing Flux Area','yellow')
					#############################################COMPUTING LINEAR AREA###################################################
					slope_line1 = (y_a-y_b)/(x_a-x_b)
					slope_line2 = (y_b-y_a)/(x_b-x_a)
					b1 = y_a - (slope_line1*x_a)
					b2 = y_b - (slope_line1*x_b)
					print
					print  colored('Computing Linear Area considering peak points:','yellow')
					print 'Point A: ',x_a,y_a
					print 'Point B: ',x_b,y_b
					print 'Slope: ',slope_line1
					print 'Slope: ',slope_line2
					print 'b: ',b1,b2
					print
					#############################################COMPUTING LINEAR AREA###################################################
					###############################################COMPUTING TOTAL AREA###############################################
					CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC1')
					AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC1')
					SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC1')

					CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC2')
					AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC2')
					SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC2')

					CTRE_G_C_PRE  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM1')
					AMPL_G_C_PRE  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM1')
					SGMA_G_C_PRE  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM1')

					CTRE_G_C_PST  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM2')
					AMPL_G_C_PST  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM2')
					SGMA_G_C_PST  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM2')

					CTRE_G_C_MDL  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM3')
					AMPL_G_C_MDL  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM3')
					SGMA_G_C_MDL  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM3')

					print
					print colored('Computing Areas using info from fits headers.','yellow')
					print 'CTRE_G_C_1: ',str(CTRE_G_C_1),' AMPL_G_C-1: ',str(AMPL_G_C_1),' SGMA_G_C-1: ',str(SGMA_G_C_1)
					print 'CTRE_G_C_2: ',str(CTRE_G_C_2),' AMPL_G_C-2: ',str(AMPL_G_C_2),' SGMA_G_C-2: ',str(SGMA_G_C_2)
					print 'CTRE_G_C_PRE: ',str(CTRE_G_C_PRE),' AMPL_G_C_PRE: ',str(AMPL_G_C_PRE),' SGMA_G_C_PRE: ',str(SGMA_G_C_PRE)
					print 'CTRE_G_C_PST: ',str(CTRE_G_C_PST),' AMPL_G_C_PST: ',str(AMPL_G_C_PST),' SGMA_G_C_PST: ',str(SGMA_G_C_PST)
					print 'CTRE_G_C_MDL: ',str(CTRE_G_C_MDL),' AMPL_G_C_MDL: ',str(AMPL_G_C_MDL),' SGMA_G_C_MDL: ',str(SGMA_G_C_PST)
					print

					W_C_1     = integrate.quad(lambda x: AMPL_G_C_1*np.exp(-((x)**2)/(2*SGMA_G_C_1**2)), -np.inf, np.inf)
					EW_C_1    = np.round(abs(np.asarray(W_C_1[0])),10)
					EWE_C_1   = np.round(abs(np.asarray(W_C_1[1])),10)

					W_C_2     = integrate.quad(lambda x: AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2)), -np.inf, np.inf)
					EW_C_2    = np.round(abs(np.asarray(W_C_2[0])),10)
					EWE_C_2   = np.round(abs(np.asarray(W_C_2[1])),10)

					W_C_PRE   = integrate.quad(lambda x: AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2)), -np.inf, np.inf)
					EW_C_PRE  = np.round(abs(np.asarray(W_C_PRE[0])),10)
					EWE_C_PRE = np.round(abs(np.asarray(W_C_PRE[1])),10)


					W_C_PST   = integrate.quad(lambda x: AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)), -np.inf, np.inf)
					EW_C_PST  = np.round(abs(np.asarray(W_C_PST[0])),10)
					EWE_C_PST = np.round(abs(np.asarray(W_C_PST[1])),10)

					W_C_MDL   = integrate.quad(lambda x: AMPL_G_C_MDL*np.exp(-((x)**2)/(2*SGMA_G_C_MDL**2)), -np.inf, np.inf)
					EW_C_MDL  = np.round(abs(np.asarray(W_C_MDL[0])),10)
					EWE_C_MDL = np.round(abs(np.asarray(W_C_MDL[1])),10)

					W_PLP     = integrate.quad(lambda x: AMPL_G_C_2*np.exp(-((x)**2)/(2*SGMA_G_C_2**2))    + 
												AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2))+ 
												AMPL_G_C_PRE*np.exp(-((x)**2)/(2*SGMA_G_C_PRE**2))+
												AMPL_G_C_PST*np.exp(-((x)**2)/(2*SGMA_G_C_PST**2)),
												-np.inf, np.inf
												)
					EW_PLP    = np.round((np.asarray(W_PLP[0])),10)
					EWE_PLP   = np.round(abs(np.asarray(W_PLP[1])),10)

					if EW_C_1 == 999999.99999:
						EW_C_1 = 0
					else:
						pass
					if EW_C_2 == 999999.99999:
						EW_C_2 = 0
					else:
						pass
					if EW_C_PRE == 999999.99999:
						EW_C_PRE = 0
					else:
						pass
					if EW_C_PST == 999999.99999:
						EW_C_PST = 0
					else:
						pass						
					EWMT    = EW_C_1 + EW_C_2 + EW_C_PRE + EW_C_PST #+ EW_C_LNR

					print
					print colored('Areas     :','yellow')
					print colored('Area CTR-1: ' + str(EW_C_1),'yellow')
					print colored('Area CTR-2: ' + str(EW_C_2),'yellow')					
					print colored('Area PRE-G: ' + str(EW_C_PRE),'blue')
					print colored('Area PST-G: ' + str(EW_C_PST),'magenta')
					print colored('Area MDL-G: ' + str(EW_C_MDL),'green')
					print colored('Area PRE-CTR1-CTR2-PST: ' + str(EW_PLP),'yellow')
					print					
					###############################################COMPUTING TOTAL AREA###############################################
					#############################################ADDING AREA TO FTIS HEADER#############################################
					if fit_vls_hdr == True:
						print
						print colored('The Areas values will be updated to the fits headers!','magenta')
						print
						print colored('Area CTR-1: '                  + str(EW_C_1)   + '-' +str(LINES[5][lines])+'_WMC1','blue')
						print colored('Area CTR-2: '                  + str(EW_C_2)   + '-' +str(LINES[5][lines])+'_WMC2','red')
						print colored('Area PRE-G: '                  + str(EW_C_PRE) + '-' +str(LINES[5][lines])+'_WGM1','cyan')
						print colored('Area PST-G: '                  + str(EW_C_PST) + '-' +str(LINES[5][lines])+'_WGM2','magenta')
						print colored('Area MDL-G: '                  + str(EW_C_MDL) + '-' +str(LINES[5][lines])+'_WGM3','green')
						print colored('Area PLP-G: '                  + str(EW_PLP)   + '-' +str(LINES[5][lines])+'_WMPP','yellow')
						print colored('Area TOT=PRE-CTR1-CTR2-PST: '  + str(EWMT)     + '-' +str(LINES[5][lines])+'_WMPT','yellow')
						print						
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC1',float(EW_C_1)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C1-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC1',float(EWE_C_1)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C1-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMC2',float(EW_C_2)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   C2-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMC2',float(EWE_C_2)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  C2-GF Crct' + str(fit_fnct) + '-' +str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM1',float(EW_C_PRE)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   PRE Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM1',float(EWE_C_PRE)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  PRE Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM2',float(EW_C_PST)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   PST Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM2',float(EWE_C_PST)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  PST Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGM3',float(EW_C_MDL)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   MDL Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGM3',float(EWE_C_MDL)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  MDL Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMPP',float(EW_C_PST)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   PRE-C1-C2-PST Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EMPP',float(EWE_C_PST)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  PRE-C1-C2-PST Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WMPT',float(EWE_C_PST)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' WE   PLP Crct' + str(fit_type))

					else:
						print
						print colored('The Areas values will NOT be updated to the fits headers!','magenta')
						print
					#############################################ADDING AREA TO FTIS HEADER#############################################
					#################################DEFINE AREA FOR LATER ON SUBSTRACT OFFSET FOR PLOTING###############################
					###########################################DEFINING PRE-PST REGIONS##################################################
					if ofs_ctr_fit == True:
						print
						print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
						print
						lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
						lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
						mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
						mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_mdl = (lambda_glx >= LINES[0][lines] + mdl_shf_ctr - mdl_shf_lim) & (lambda_glx <= LINES[0][lines] + mdl_shf_ctr + mdl_shf_lim)
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						X0_f2DG_indx_MDL       = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
					else:
						print
						print colored('Using Expected Line Center to Define Fitting Region!','yellow')
						print
						mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
						mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_mdl = (lambda_glx >= LINES[0][lines] + mdl_shf_ctr - mdl_shf_lim) & (lambda_glx <= LINES[0][lines] + mdl_shf_ctr + mdl_shf_lim)
						try:
							X0_f2DG_indx_PRE = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						except ValueError:
							X0_f2DG_indx_PRE = LINES[0][lines]

						try:
							X0_f2DG_indx_PST = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						except ValueError:
							X0_f2DG_indx_PST = LINES[0][lines]

						try:
							X0_f2DG_indx_MDL = np.where(inten_glx[mask_ft_mdl]==(max(inten_glx[mask_ft_mdl])))[0]
						except ValueError:
							X0_f2DG_indx_MDL = LINES[0][lines]

						
					print
					print colored('Region limits for fitting.','yellow')
					print colored('Central    : ' + str(LINES[0][lines]),'yellow')
					print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'cyan')
					print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'magenta')
					print colored('Central+MDL: ' + str(mdl_shf_ctr) + '-' + str(LINES[0][lines]+mdl_shf_ctr),'green')
					print
					print colored('Limits:','yellow')
					print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
					print 'PRE    : lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim,'-','lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim
					print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
					print 'MDL    : lambda_glx >= ',LINES[0][lines] + mdl_shf_ctr-mdl_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + mdl_shf_ctr+mdl_shf_lim
					print
					print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
					print
					############################################DEFINING PRE-PST REGIONS##################################################
					##################################DEFINE AREA FOR LATER ON SUBSTRACT OFFSET FOR PLOTING###############################
				elif 'Dbl' not in LINES[3][lines] and fit_fnct=='gauss'   and fit_type == 'lmfit' and uft_lne_vls == True:
					print 'Line-fitting Cleaning Method. To be checked Fnc_Stk_Plt.py def(Plot_Idp_Spc_Lne) 7554!'
					
					fit_typ = 'G'
					#GAUSSIAN FIT
					MSK_NTMS=2.5

					lmb_min_lim_line_ft = (LINES[0][lines]-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
					lmb_max_lim_line_ft = (LINES[0][lines]+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
					lmb_min_lim_line    = (LINES[0][lines]-LINES[8][lines])*(1+z_glx_Ps) - 1.5*MSK_NTMS*LINES[2][lines]
					lmb_max_lim_line    = (LINES[0][lines]+LINES[8][lines])*(1+z_glx_Ps) + 1.5*MSK_NTMS*LINES[2][lines]

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)

					popt_0, pcov_0     = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
					perr_0             = [999999.99999,999999.99999,999999.99999]
					CTRE_G_0           = 999999.99999
					AMPL_G_0           = 999999.99999
					SGMA_G_0           = 999999.99999
					FWHM_G_0           = 999999.99999
					EW_0               = 999999.99999
					EWE_0              = 999999.99999

					CTRE_G_0_E         = 999999.99999
					AMPL_G_0_E         = 999999.99999
					SGMA_G_0_E         = 999999.99999

					CTRE_G_0_cor       = 999999.99999
					AMPL_G_0_cor       = 999999.99999
					SGMA_G_0_cor       = 999999.99999

					chisqr_0           = 999999.99999
					redchi_0           = 999999.99999
					if fit_vls_hdr == True:
						Header_Add(specfile_glx,str(LINES[5][lines])+'_CGF0',float(CTRE_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_AGF0',float(AMPL_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_FGF0',float(FWHM_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGF0',float(EW_0)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGF0',float(EWE_0)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct')
					else:
						print
						print colored('The fit values will NOT be added to the fits headers!','magenta')
						print

					perr_C          = [999999.99999,999999.99999,999999.99999]
					CTRE_G_C        = 999999.99999
					AMPL_G_C        = 999999.99999
					SGMA_G_C        = 999999.99999
					FWHM_G_C        = 999999.99999
					EW_C            = 999999.99999
					EWE_C           = 999999.99999

					CTRE_G_C_E      = 999999.99999
					AMPL_G_C_E      = 999999.99999
					SGMA_G_C_E      = 999999.99999
					CTRE_G_C_cor    = 999999.99999
					AMPL_G_C_cor    = 999999.99999
					SGMA_G_C_cor    = 999999.99999
					chisqr_C        = 999999.99999
					redchi_C        = 999999.99999

					popt_O, pcov_O  = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
					perr_O          = [999999.99999,999999.99999,999999.99999,999999.99999]
					CTRE_G_O        = 999999.99999
					AMPL_G_O        = 999999.99999
					SGMA_G_O        = 999999.99999
					OFST_G_O        = 999999.99999
					FWHM_G_O        = 999999.99999
					EW_O            = 999999.99999
					EWE_O           = 999999.99999

					CTRE_G_O_E      = 999999.99999
					AMPL_G_O_E      = 999999.99999
					SGMA_G_O_E      = 999999.99999
					CTRE_G_O_cor    = 999999.99999
					AMPL_G_O_cor    = 999999.99999
					SGMA_G_O_cor    = 999999.99999
					OFST_G_O_cor    = 999999.99999
					chisqr_O        = 999999.99999
					redchi_O        = 999999.99999

					AMPL_SNR        = 999999.99999
					CTRE_SNR        = 999999.99999
					SGMA_SNR        = 999999.99999
					print
					print colored(specfile_glx,'cyan')
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines]),'yellow')
					print
					print colored(str(LINES[5][lines])+'_CGLC: ' + str(CTRE_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_AGLC: ' + str(AMPL_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_SGLC: ' + str(SGMA_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_FGLC: ' + str(FWHM_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_WGLC: ' + str(EW_C)      ,'yellow')
					print colored(str(LINES[5][lines])+'_EGLC: ' + str(EWE_C),'yellow')
					print
					print
					if fit_vls_hdr == True:
						Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                       ,header_comment = 'Fit function')
						Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                       ,header_comment = 'Fit method')
						Header_Add(specfile_glx,str(LINES[5][lines])+'_CGLO',float(CTRE_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_AGLO',float(AMPL_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_FGLO',float(FWHM_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGLO',float(EW_O)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGLO',float(EWE_O)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Offst' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_OFSO',float(OFST_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 1GF Offst' + str(fit_type))

						Header_Add(specfile_glx,str(LINES[5][lines])+'_CGLC',float(CTRE_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_AGLC',float(AMPL_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_SGLC',float(SGMA_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_FGLC',float(FWHM_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WGLC',float(EW_C)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_EGLC',float(EWE_C)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))

						Header_Add(specfile_glx,str(LINES[5][lines])+'_CLEC',float(CTRE_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_ALEC',float(AMPL_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_SLEC',float(SGMA_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))

						Header_Add(specfile_glx,str(LINES[5][lines])+'_CHGL',float(chisqr_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_CRGL',float(redchi_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 Reduced 1GF' + str(fit_type))
					else:
						print
						print colored('The fit values will NOT be added to the fits headers!','magenta')
						print
						pass
				elif 'Dbl' not in LINES[3][lines] and fit_fnct=='gaussM'  and fit_type == 'lmfit' and uft_lne_vls == True:
					fit_typ = 'GM'
					#GAUSSIAN FIT
					MSK_NTMS=2.5

					lmb_min_lim_line_ft = (LINES[0][lines]-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
					lmb_max_lim_line_ft = (LINES[0][lines]+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
					lmb_min_lim_line    = (LINES[0][lines]-LINES[8][lines])*(1+z_glx_Ps) - 1.5*MSK_NTMS*LINES[2][lines]
					lmb_max_lim_line    = (LINES[0][lines]+LINES[8][lines])*(1+z_glx_Ps) + 1.5*MSK_NTMS*LINES[2][lines]

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function


									
					##################################################CENTRAL GAUSSIAN###################################################
					if fix_ctr_gau == False:
						print
						print colored('0-Fitting Central line','cyan')
						print
						#################################################CENTRAL GAUSSIAN-0##################################################

						popt_0, pcov_0   = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_0           = [999999.99999,999999.99999,999999.99999]
						CTRE_G_0         = 999999.99999
						AMPL_G_0         = 999999.99999
						SGMA_G_0         = 999999.99999
						FWHM_G_0         = 999999.99999
						EW_0             = 999999.99999
						EWE_0            = 999999.99999

						CTRE_G_0_E      = 999999.99999
						AMPL_G_0_E      = 999999.99999
						SGMA_G_0_E      = 999999.99999

						CTRE_G_0_cor    = 999999.99999
						AMPL_G_0_cor    = 999999.99999
						SGMA_G_0_cor    = 999999.99999

						chisqr_0        = 999999.99999
						redchi_0        = 999999.99999
						if fit_vls_hdr == True and fix_ctr_gau ==False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CF0M',float(CTRE_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Ctr  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AF0M',float(AMPL_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'Amp  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FF0M',float(FWHM_G_0) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'FWHM 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WF0M',float(EW_0)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EW   1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EF0M',float(EWE_0)    ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + 'EWE  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							print
							print colored('The fit (CTR-0) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-0) values will NOT be added to the fits headers!','magenta')
							print
						#################################################CENTRAL GAUSSIAN-0##################################################
						#################################################CENTRAL GAUSSIAN-C##################################################
						popt_C, pcov_C = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C          = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C        = 999999.99999
						AMPL_G_C        = 999999.99999
						SGMA_G_C        = 999999.99999
						FWHM_G_C        = 999999.99999
						EW_C            = 999999.99999
						EWE_C           = 999999.99999

						CTRE_G_C_E      = 999999.99999
						AMPL_G_C_E      = 999999.99999
						SGMA_G_C_E      = 999999.99999
						CTRE_G_C_cor    = 999999.99999
						AMPL_G_C_cor    = 999999.99999
						SGMA_G_C_cor    = 999999.99999
						chisqr_C        = 999999.99999
						redchi_C        = 999999.99999

						popt_O ,pcov_O  = [999999.99999,999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999,999999.99999]
						perr_O          = [999999.99999,999999.99999,999999.99999,999999.99999]
						CTRE_G_O        = 999999.99999
						AMPL_G_O        = 999999.99999
						SGMA_G_O        = 999999.99999
						OFST_G_O        = 999999.99999
						FWHM_G_O        = 999999.99999
						EW_O            = 999999.99999
						EWE_O           = 999999.99999

						CTRE_G_O_E      = 999999.99999
						AMPL_G_O_E      = 999999.99999
						SGMA_G_O_E      = 999999.99999
						CTRE_G_O_cor    = 999999.99999
						AMPL_G_O_cor    = 999999.99999
						SGMA_G_O_cor    = 999999.99999
						OFST_G_O_cor    = 999999.99999
						chisqr_O        = 999999.99999
						redchi_O        = 999999.99999

						AMPL_SNR     = 999999.99999
						CTRE_SNR     = 999999.99999
						SGMA_SNR     = 999999.99999
						print
						print colored(specfile_glx,'cyan')
						print
						print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines]),'yellow')
						print
						print colored(str(LINES[5][lines])+'_CGLC: ' + str(CTRE_G_C)  ,'yellow')
						print colored(str(LINES[5][lines])+'_AGLC: ' + str(AMPL_G_C)  ,'yellow')
						print colored(str(LINES[5][lines])+'_SGLC: ' + str(SGMA_G_C)  ,'yellow')
						print colored(str(LINES[5][lines])+'_FGLC: ' + str(FWHM_G_C)  ,'yellow')
						print colored(str(LINES[5][lines])+'_WGLC: ' + str(EW_C)      ,'yellow')
						print colored(str(LINES[5][lines])+'_EGLC: ' + str(EWE_C),'yellow')
						print
						print


						if fit_vls_hdr == True and fix_ctr_gau==False:
							Header_Add(specfile_glx,'BSF_FCT',str(fit_fnct)                       ,header_comment = 'Fit function')
							Header_Add(specfile_glx,'BSF_MTH',str(fit_type)                       ,header_comment = 'Fit method')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLOM',float(CTRE_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALOM',float(AMPL_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLOM',float(FWHM_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLOM',float(EW_O)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELOM',float(EWE_O)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Offst' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_OFOM',float(OFST_G_O)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ofst 1GF Offst' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLCM',float(CTRE_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALCM',float(AMPL_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLCM',float(SGMA_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FLCM',float(FWHM_G_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WLCM',float(EW_C)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ELCM',float(EWE_C)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CECM',float(CTRE_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AECM',float(AMPL_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SECM',float(SGMA_G_C_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))

							except ValueError:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CECM',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_AECM',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SECM',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_fnct) + '-' +str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHLM',float(chisqr_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_fnct) + '-' +str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRLM',float(redchi_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 Reduced 1GF' + str(fit_fnct) + '-' +str(fit_type))

							print
							print colored('The fit (CTR-C & CTR_O) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (CTR-C & CTR_O) values will NOT be added to the fits headers!','magenta')
							print
						#################################################CENTRAL GAUSSIAN-C##################################################					
					elif fix_ctr_gau == True:
						print
						print colored('0 CTR-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')

						CTRE_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CF0M')
						AMPL_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AF0M')
						FWHM_G_0    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FF0M')
						EW_0        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WF0M')
						EWE_0       = Header_Get(specfile_glx,str(LINES[5][lines])+'_EF0M')

						CTRE_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLOM')
						AMPL_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALOM')
						FWHM_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLOM')
						EW_O        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLOM')
						EWE_O       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELOM')
						OFST_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_OFOM')

						CTRE_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLCM')
						AMPL_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALCM')
						SGMA_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLCM')
						FWHM_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLCM')
						EW_C        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLCM')
						EWE_C       = Header_Get(specfile_glx,str(LINES[5][lines])+'_ELCM')
						CTRE_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_CECM')
						AMPL_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_AECM')
						SGMA_G_C_E  = Header_Get(specfile_glx,str(LINES[5][lines])+'_SECM')

						chisqr_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CHLM')
						redchi_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CRLM')

					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-CTR','yellow')
					print
					print colored(str(LINES[5][lines])+'_CLCM: ' + str(CTRE_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_ALCM: ' + str(AMPL_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_SLCM: ' + str(SGMA_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_FLCM: ' + str(FWHM_G_C)  ,'yellow')
					print colored(str(LINES[5][lines])+'_WLCM: ' + str(EW_C)      ,'yellow')
					print colored(str(LINES[5][lines])+'_ELCM: ' + str(EWE_C),'yellow')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','yellow')
					print colored(str(CTRE_G_C)+', '+str(AMPL_G_C)+', '+str(SGMA_G_C),'yellow')
					print
					##################################################CENTRAL GAUSSIAN###################################################

					#####################################################PRE GAUSSIAN####################################################
					if fix_pre_gau == False:
						###########################################DEFINING PRE-PST REGIONS##################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
							lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						else:
							print
							print colored('Using Expected Line Center to Define Fitting Region!','yellow')
							print
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]

						print
						print colored('Region limits for fitting.','yellow')
						print colored('Central    : ' + str(LINES[0][lines]),'yellow')
						print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'yellow')
						print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'yellow')
						print
						print colored('Limits:','yellow')
						print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
						print 'PRE    : lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim,'-','lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim
						print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
						print
						print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
						print
						###########################################DEFINING PRE-PST REGIONS##################################################

						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						x_a = lambda_glx[mask_ft_pre][X0_f2DG_indx_PRE]
						y_a = inten_glx[mask_ft_pre][X0_f2DG_indx_PRE]

						popt_C_PRE, pcov_C_PRE  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C_PRE          = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C_PRE        = 999999.99999
						AMPL_G_C_PRE        = 999999.99999
						SGMA_G_C_PRE        = 999999.99999
						FWHM_G_C_PRE        = 999999.99999
						EW_C_PRE            = 999999.99999
						EWE_C_PRE           = 999999.99999

						CTRE_G_C_PRE_E      = 999999.99999
						AMPL_G_C_PRE_E      = 999999.99999
						SGMA_G_C_PRE_E      = 999999.99999
						CTRE_G_C_PRE_cor    = 999999.99999
						AMPL_G_C_PRE_cor    = 999999.99999
						SGMA_G_C_PRE_cor    = 999999.99999
						chisqr_C_PRE        = 999999.99999
						redchi_C_PRE        = 999999.99999

						EW_C_PR1            = 999999.99999
						EWE_C_PR1           = 999999.99999
						if fit_vls_hdr == True and fix_pre_gau == False:						
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CHGL',float(chisqr_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 1GF' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CRGL',float(redchi_C)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Chi2 Reduced 1GF' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGL1',float(CTRE_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGL1',float(AMPL_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGL1',float(SGMA_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGL1',float(FWHM_G_C_PRE)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGL1',float(EW_C_PRE)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGL1',float(EWE_C_PRE)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							try:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CLE1',float(CTRE_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_ALE1',float(AMPL_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SLE1',float(SGMA_G_C_PRE_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))
							except:
								Header_Add(specfile_glx,str(LINES[5][lines])+'_CLE1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_ALE1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
								Header_Add(specfile_glx,str(LINES[5][lines])+'_SLE1',float(999999.99999),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))


							Header_Add(specfile_glx,str(LINES[5][lines])+'_XA1',float(x_a),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PRE GAU-LNR X1 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YA1',float(y_a),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PRE GAU-LNR Y1 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WPR1',float(EW_C_PR1),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])   + ' EW   PRE' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EPR1',float(EWE_C_PR1),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE  PRE' + str(fit_type))
							print
							print colored('The fit (PRE) values will be added to the fits headers!','magenta')
							print
						else:
							print
							print colored('The fit (PRE) values will NOT be added to the fits headers!','magenta')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GS1',float(pre_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for PRE G-1')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GC1',float(pre_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for PRE G-1')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass						
					elif fix_pre_gau == True:
						pass
						
						print
						print colored('1 PRE-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						CTRE_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGL1')
						AMPL_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGL1')
						SGMA_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGL1')
						FWHM_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGL1')
						EW_C_PRE       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGL1')
						EWE_C_PRE      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGL1')
						CTRE_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLE1')
						AMPL_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALE1')
						SGMA_G_C_PRE_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLE1')

						x_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XA1')
						y_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YA1')
						EW_C_PR1       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WPR1')
						EWE_C_PR1      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EPR1')

					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-PRE','yellow')
					print
					print colored(str(LINES[5][lines])+'_CGL1: ' + str(CTRE_G_C_PRE)  ,'yellow')
					print colored(str(LINES[5][lines])+'_AGL1: ' + str(AMPL_G_C_PRE)  ,'yellow')
					print colored(str(LINES[5][lines])+'_SGL1: ' + str(SGMA_G_C_PRE)  ,'yellow')
					print colored(str(LINES[5][lines])+'_FGL1: ' + str(FWHM_G_C_PRE)  ,'yellow')
					print colored(str(LINES[5][lines])+'_WGL1: ' + str(EW_C_PRE)      ,'yellow')
					print colored(str(LINES[5][lines])+'_EGL1: ' + str(EWE_C_PRE),'yellow')
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','yellow')
					print colored(str(CTRE_G_C_PRE)+', '+str(AMPL_G_C_PRE)+', '+str(SGMA_G_C_PRE),'yellow')
					print
					#####################################################PRE GAUSSIAN#################################################
					###################################################POST GAUSSIAN##################################################
					if fix_pst_gau == False:
						###########################################DEFINING PRE-PST REGIONS##################################################
						if ofs_ctr_fit == True:
							print
							print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
							print
							lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines]
							lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						else:
							print
							print colored('Using Expected Line Center to Define Fitting Region!','yellow')
							print
							mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
							mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
							X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
							X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]

						print
						print colored('Region limits for fitting.','yellow')
						print colored('Central    : ' + str(LINES[0][lines]),'yellow')
						print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'yellow')
						print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'yellow')
						print
						print colored('Limits:','yellow')
						print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
						print 'PRE    : lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim,'-','lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim
						print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
						print
						print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
						print
						###########################################DEFINING PRE-PST REGIONS##################################################

						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
						x_b = lambda_glx[mask_ft_pst][X0_f2DG_indx_PST]
						y_b = inten_glx[mask_ft_pst][X0_f2DG_indx_PST]

						popt_C_PST, pcov_C_PST  = [999999.99999,999999.99999,999999.99999],[999999.99999,999999.99999,999999.99999]
						perr_C_PST          = [999999.99999,999999.99999,999999.99999]
						CTRE_G_C_PST        = 999999.99999
						AMPL_G_C_PST        = 999999.99999
						SGMA_G_C_PST        = 999999.99999
						FWHM_G_C_PST        = 999999.99999
						EW_C_PST            = 999999.99999
						EWE_C_PST           = 999999.99999

						
						CTRE_G_C_PST_E      = 999999.99999
						AMPL_G_C_PST_E      = 999999.99999
						SGMA_G_C_PST_E      = 999999.99999
						CTRE_G_C_PST_cor    = 999999.99999
						AMPL_G_C_PST_cor    = 999999.99999
						SGMA_G_C_PST_cor    = 999999.99999
						chisqr_C_PST        = 999999.99999
						redchi_C_PST        = 999999.99999

						EW_C_PR2            = 999999.99999
						EWE_C_PR2           = 999999.99999

						if fit_vls_hdr == True and fix_pst_gau == False:
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CGL2',float(CTRE_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_AGL2',float(AMPL_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SGL2',float(SGMA_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_FGL2',float(FWHM_G_C_PST)  ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' FWHM 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WGL2',float(EW_C_PST)      ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EW   1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EGL2',float(EWE_C_PST)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_CLE2',float(CTRE_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Ctr  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_ALE2',float(AMPL_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Amp  err 1GF Crct' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_SLE2',float(SGMA_G_C_PST_E),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Sigm err 1GF Crct' + str(fit_type))

							Header_Add(specfile_glx,str(LINES[5][lines])+'_XA2',float(x_b),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PST GAU-LNR X2 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_YA2',float(y_b),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' PST GAU-LNR Y2 COO')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_WPS2',float(EW_C_PST),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])   + ' EW   PST' + str(fit_type))
							Header_Add(specfile_glx,str(LINES[5][lines])+'_EPS2',float(EWE_C_PST),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EWE  PST' + str(fit_type))
						else:
							print
							print colored('The fit (PST) values will NOT be added to the fits headers!','magenta')
							print
							pass
						if uft_lne_vls == False and fit_vls_hdr == True and int_vlf_hdr==True:
							print
							print colored('Initial Guess Values for line Fitting will be recorded!','yellow')
							print colored('Info input parameters (pre_gauss_shf, pre_gauss_ctr) in Spec_Stk_Stack_0.4.py!','yellow')
							print
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GS2',float(pst_shf_lim) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Shf from expt ctr wvl for PST G-2')
							Header_Add(specfile_glx,str(LINES[5][lines])+'_GC2',float(pst_shf_ctr) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' Gaussian Ctr Int Val for PST G-2')
						else:
							print
							print colored('Initial Guess Values for line Fitting will NOT be recorded!','yellow')
							print
							pass
					elif fix_pst_gau == True:
						pass
						
						print
						print colored('2 PST-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						CTRE_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGL2')
						AMPL_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGL2')
						SGMA_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGL2')
						FWHM_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGL2')
						EW_C_PST       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGL2')
						EWE_C_PST      = Header_Get(specfile_glx,str(LINES[5][lines])+'_EGL2')
						CTRE_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLE2')
						AMPL_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALE2')
						SGMA_G_C_PST_E = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLE2')

						x_b         = Header_Get(specfile_glx,str(LINES[5][lines])+'_XA2')
						y_b         = Header_Get(specfile_glx,str(LINES[5][lines])+'_YA2')
						EW_C_PS2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_WPS2')
						EWE_C_PS2   = Header_Get(specfile_glx,str(LINES[5][lines])+'_EPS2')
						
					print
					print colored(str(LINES[0][lines])+'-'+str(LINES[3][lines])+'-PST','yellow')
					print
					print colored(str(LINES[5][lines])+'_CGL2: ' + str(CTRE_G_C_PST)  ,'yellow')
					print colored(str(LINES[5][lines])+'_AGL2: ' + str(AMPL_G_C_PST)  ,'yellow')
					print colored(str(LINES[5][lines])+'_SGL2: ' + str(SGMA_G_C_PST)  ,'yellow')
					print colored(str(LINES[5][lines])+'_FGL2: ' + str(FWHM_G_C_PST)  ,'yellow')
					print colored(str(LINES[5][lines])+'_WGL2: ' + str(EW_C_PST)      ,'yellow')
					print colored(str(LINES[5][lines])+'_EGL2: ' + str(EWE_C_PST),'yellow')
					print
					print colored('Fit Values Center, Amplitude, Sigma ('+fit_type+'):','yellow')
					print colored(str(CTRE_G_C_PST)+', '+str(AMPL_G_C_PST)+', '+str(SGMA_G_C_PST),'yellow')
					print					
					###################################################POST GAUSSIAN##################################################										

					###############################################COMPUTING TOTAL AREA###############################################
					print colored('Computing Flux Area','yellow')
					#############################################COMPUTING LINEAR AREA###################################################
					slope_line1 = (y_a-y_b)/(x_a-x_b)
					slope_line2 = (y_b-y_a)/(x_b-x_a)
					b1 = y_a - (slope_line1*x_a)
					b2 = y_b - (slope_line1*x_b)
					print
					print  colored('Computing Linear Area considering peak points:','yellow')
					print 'Point A: ',x_a,y_a
					print 'Point B: ',x_b,y_b
					print 'Slope: ',slope_line1
					print 'Slope: ',slope_line2
					print 'b: ',b1,b2
					print
					W_C_LNR   = integrate.quad(lambda x: slope_line1*(x) + b1 - 1, x_a, x_b)
					EW_C_LNR  = np.round(abs(np.asarray(W_C_LNR[0])),3)
					EWE_C_LNR = np.round(abs(np.asarray(W_C_LNR[1])),10)

					print (slope_line1*(x_a**2))/2 + (b1*x_a)
					print (slope_line1*(x_b**2))/2 + (b1*x_b)
					print ((slope_line1*(x_a**2))/2 + (b1*x_a)) - ((slope_line1*(x_b**2))/2 + (b1*x_b))
					#############################################COMPUTING LINEAR AREA###################################################

					if EW_C_PRE == 999999.99999:
						EW_C_PRE = 0
					else:
						pass
					if EWE_C_PST == 999999.99999:
						EWE_C_PST = 0
					else:
						pass

					EWMT       = EW_C - EW_C_PRE - EWE_C_PST
					print
					print colored('Area CTR-G: ' + str(EW_C),'yellow')
					print colored('Area PRE-G: ' + str(EW_C_PR1),'blue')
					print colored('Area PST-G: ' + str(EW_C_PST),'magenta')
					print colored('Area LNR:   ' + str(EW_C_LNR),'cyan')
					print colored('Area TOT:   ' + str(EWMT),'yellow')
					print
					###############################################COMPUTING TOTAL AREA###############################################
					if fit_vls_hdr == True:
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WR1',float(EW_C_LNR),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])  + ' EW   LNR' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_ER1',float(EWE_C_LNR),header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' EWE  LNR' + str(fit_type))
						Header_Add(specfile_glx,str(LINES[5][lines])+'_WTOT',float(EWMT),header_comment = str(LINES[3][lines]) + str(LINES[0][lines])       + ' EW   TOT' + str(fit_type))
					else:
						print
						print colored('The LINEAR AND TOTAL Areas values will NOT be updated to the fits headers!','magenta')
						print
					#################################DEFINE AREA FOR LATER ON SUBSTRACT OFFSET FOR PLOTING###############################
					###########################################DEFINING PRE-PST REGIONS##################################################
					if ofs_ctr_fit == True:
						print
						print colored('Using Fitted Center From Offset Fit to Redefine Fitting Region!','yellow')
						print
						lmb_min_lim_line_ft = (CTRE_G_0-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
						lmb_max_lim_line_ft = (CTRE_G_0+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
						mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
						mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]
					else:
						print
						print colored('Using Expected Line Center to Define Fitting Region!','yellow')
						print
						mask_ft_pre = (lambda_glx <= LINES[0][lines] - pre_shf_ctr + pre_shf_lim) & (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim)
						mask_ft_pst = (lambda_glx >= LINES[0][lines] + pst_shf_ctr - pst_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						mask_ft_plp = (lambda_glx >= LINES[0][lines] - pre_shf_ctr - pre_shf_lim) & (lambda_glx <= LINES[0][lines] + pst_shf_ctr + pst_shf_lim)
						X0_f2DG_indx_PRE       = np.where(inten_glx[mask_ft_pre]==(max(inten_glx[mask_ft_pre])))[0]
						X0_f2DG_indx_PST       = np.where(inten_glx[mask_ft_pst]==(max(inten_glx[mask_ft_pst])))[0]

					print
					print colored('Region limits for fitting.','yellow')
					print colored('Central    : ' + str(LINES[0][lines]),'yellow')
					print colored('Central-PRE: ' + str(pre_shf_ctr) + '-' + str(LINES[0][lines]-pre_shf_ctr),'yellow')
					print colored('Central+PST: ' + str(pst_shf_ctr) + '-' + str(LINES[0][lines]+pst_shf_ctr),'yellow')
					print
					print colored('Limits:','yellow')
					print 'Central: lambda_glx >= ',lmb_min_lim_line_ft  ,'-','lambda_glx <= ',lmb_max_lim_line_ft
					print 'PRE    : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pre_shf_lim
					print 'PST    : lambda_glx >= ',LINES[0][lines] + pst_shf_ctr-pst_shf_lim,'-','lambda_glx <= ',LINES[0][lines] + pst_shf_ctr+pst_shf_lim
					print
					print 'Total  : lambda_glx >= ',LINES[0][lines] - pre_shf_ctr-pre_shf_lim,'-','lambda_glx <= ',LINES[0][lines] - pre_shf_ctr+pst_shf_lim
					print

					###########################################DEFINING PRE-PST REGIONS##################################################
					#################################DEFINE AREA FOR LATER ON SUBSTRACT OFFSET FOR PLOTING###############################
				elif 'Dbl' in LINES[3][lines] and fit_fnct=='gauss' and mke_lne_fit == False:
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					print 
					print colored('Double Line Fit (Ind)','yellow')
					print colored(LINES[3][lines-2]  + '-' + str(LINES[0][lines-2])   + '-' + str(LINES[1][lines-2]),'cyan')
					print LINES[3][lines]+ '-' + str(LINES[0][lines]) + '-'  + str(LINES[1][lines])
					print
					ivl_fts_hdr = True
					########################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#########################
					if ivl_fts_hdr == True:
						try:
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial fit variables from fits header!','yellow')
							print
							print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WF02' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WP02' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CF02' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CO02' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AF02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print
						except ValueError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print '*****'
							print
							quit()
						try:
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')							
							print colored(str(LINES[5][lines-2])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')							
							print colored(str(LINES[5][lines-1])+'_AF02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF01',0.001,header_comment = str(LINES[3][lines-2]) + str(LINES[0][lines-2]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF02',0.001,header_comment = str(LINES[3][lines-1]) + str(LINES[0][lines-1]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines-2])+'_AF01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines-1])+'_AF02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('*****Success!******','yellow')
							print
						################FOR CASES WHERE KLINE WIIDHT ==0################
						if L1_1  == 0:
							L1_1  = 1#LINES[1][lines]
						else:
							pass
						if L1_2  == 0:
							L1_2  = 1#LINES[1][lines]
						else:
							pass
						################FOR CASES WHERE KLINE WIIDHT ==0################
					elif ivl_fts_hdr == False:
						L1_1  = LINES[1][lines-2]
						L2_1  = LINES[2][lines-2]
						L7_1  = LINES[7][lines-2]
						L8_1  = LINES[8][lines-2]
						L10_1 = LINES[10][lines-2]

						L1_2  = LINES[1][lines-1]
						L2_2  = LINES[2][lines-1]
						L7_2  = LINES[7][lines-1]
						L8_2  = LINES[8][lines-1]
						L10_2 = LINES[10][lines-1]
						print
						print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
						print					
					print colored('Initial Values: ','cyan')
					print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')					
					print colored(str(LINES[5][lines])+'_WF01' + ': ' + str(L1_1),'cyan')
					print colored(str(LINES[5][lines])+'_WP01' + ': ' + str(L2_1),'cyan')
					print colored(str(LINES[5][lines])+'_CF01' + ': ' + str(L7_1),'cyan')
					print colored(str(LINES[5][lines])+'_CO01' + ': ' + str(L8_1),'cyan')
					print colored(str(LINES[5][lines])+'_AF01' + ': ' + str(L10_1),'cyan')
					print
					print colored('Initial Values: ','magenta')
					print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')					
					print colored(str(LINES[5][lines])+'_WF02' + ': ' + str(L1_2),'magenta')
					print colored(str(LINES[5][lines])+'_WP02' + ': ' + str(L2_2),'magenta')
					print colored(str(LINES[5][lines])+'_CF02' + ': ' + str(L7_2),'magenta')
					print colored(str(LINES[5][lines])+'_CO02' + ': ' + str(L8_2),'magenta')
					print colored(str(LINES[5][lines])+'_AF02' + ': ' + str(L10_2),'magenta')
					print
					########################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#########################
					##################################################CENTRAL GAUSSIAN-1###################################################
					lmb_min_lim_line_ft = (LINES[0][lines-2]+L8_1) - MSK_NTMS*L1_1 
					lmb_max_lim_line_ft = (LINES[0][lines-2]+L8_1) + MSK_NTMS*L1_1
					lmb_min_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_1
					lmb_max_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_1

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] 

					CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC1')
					AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC1')
					SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC1')
					FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC1')

					CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLC2')
					AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALC2')
					SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLC2')
					FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLC2')
				elif 'Dbl' in LINES[3][lines] and fit_fnct=='gaussM'and mke_lne_fit == False:
					fit_typ = 'GM'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					print 
					print colored('Double Line Fit (Ind)','yellow')
					print colored(LINES[3][lines-2]  + '-' + str(LINES[0][lines-2])   + '-' + str(LINES[1][lines-2]),'cyan')
					print LINES[3][lines]+ '-' + str(LINES[0][lines]) + '-'  + str(LINES[1][lines])
					print colored(LINES[3][lines-1]+ '-' + str(LINES[0][lines-1]) + '-' + str(LINES[1][lines-1]),'magenta')
					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
					ivl_fts_hdr = True
					if ivl_fts_hdr == True:
						try:
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WM02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CM02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('Initial fit variables from fits header!','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WM01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WM01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CM01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WM02' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WM02' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CM02' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','cyan')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','magenta')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print '*****'
							print
							print
							print colored('Initial fit variables from fits header!','yellow')
							L1_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF01')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP01')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF01')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_1  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO01')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF01')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							L1_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF02')        #LINES-1  Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP02')        #LINES-2  Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF02')        #LINES-7  Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_2  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO02')        #LINES-8  Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF02')        #LINES-10 Ctr Fit Amp   1GF-IntVal LNE_AMP_BNDS

							print
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored('Headers:','cyan')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_1),'cyan')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_1),'cyan')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_1),'cyan')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_1),'cyan')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_1),'cyan')
							print
							print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')
							print colored('Headers:','magenta')
							print colored(str(LINES[5][lines])+'_WF01' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_2),'magenta')
							print colored(str(LINES[5][lines])+'_WP01' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_2),'magenta')
							print colored(str(LINES[5][lines])+'_CF01' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_2),'magenta')
							print colored(str(LINES[5][lines])+'_CO01' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_2),'magenta')
							print colored(str(LINES[5][lines])+'_AF01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_2),'magenta')
							print colored('*****Success!******','yellow')
							print							
							#quit()
						try:
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_2 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							print colored('*****Success!******','yellow')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored('1: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'cyan')							
							print colored(str(LINES[5][lines])+'_AM01' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','cyan')
							print colored('2: ' + LINES[3][lines]  + '-' + str(LINES[0][lines]),'magenta')							
							print colored(str(LINES[5][lines])+'_AM02' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','magenta')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AM01',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AM02',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM01')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							L10_1 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AM02')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('*****Success!******','yellow')
							print
						################FOR CASES WHERE KLINE WIIDHT ==0################
						if L1_1  == 0:
							L1_1  = 1#LINES[1][lines]
						else:
							pass
						if L1_2  == 0:
							L1_2  = 1#LINES[1][lines]
						else:
							pass
						################FOR CASES WHERE KLINE WIIDHT ==0################
					elif ivl_fts_hdr == False:
						L1_1  = LINES[1][lines-2]
						L2_1  = LINES[2][lines-2]
						L7_1  = LINES[7][lines-2]
						L8_1  = LINES[8][lines-2]
						L10_1 = LINES[10][lines-2]

						L1_2  = LINES[1][lines-1]
						L2_2  = LINES[2][lines-1]
						L7_2  = LINES[7][lines-1]
						L8_2  = LINES[8][lines-1]
						L10_2 = LINES[10][lines-1]
						print
						print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
						print					
					print colored('Initial Values: ','cyan')
					print colored('1: ' + LINES[3][lines-2]  + '-' + str(LINES[0][lines-2]),'cyan')					
					print colored(str(LINES[5][lines-2])+'_WF01' + ': ' + str(L1_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_WP01' + ': ' + str(L2_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_CF01' + ': ' + str(L7_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_CO01' + ': ' + str(L8_1),'cyan')
					print colored(str(LINES[5][lines-2])+'_AF01' + ': ' + str(L10_1),'cyan')
					print
					print colored('Initial Values: ','magenta')
					print colored('2: ' + LINES[3][lines-1]  + '-' + str(LINES[0][lines-1]),'magenta')					
					print colored(str(LINES[5][lines-1])+'_WF02' + ': ' + str(L1_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_WP02' + ': ' + str(L2_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_CF02' + ': ' + str(L7_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_CO02' + ': ' + str(L8_2),'magenta')
					print colored(str(LINES[5][lines-1])+'_AF02' + ': ' + str(L10_2),'magenta')
					print
					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
					lmb_min_lim_line_ft = (LINES[0][lines-2]+L8_1) - MSK_NTMS*L1_1 
					lmb_max_lim_line_ft = (LINES[0][lines-2]+L8_1) + MSK_NTMS*L1_1
					lmb_min_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) - 1.5*MSK_NTMS*L2_1
					lmb_max_lim_line    = (LINES[0][lines-2]+L8_1)*(1+z_glx_Ps) + 1.5*MSK_NTMS*L2_1

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0]

					print
					print colored('1-CTR-gaussian already fitted!','yellow')
					print colored('Values from fits header','yellow')
					try:
						CTRE_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC1')
						AMPL_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC1')
						SGMA_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC1')
						FWHM_G_C_1    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMC1')
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines-2])+'_CF01','yellow')
						print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
						print colored ('Or UnFix (fix_ctr_gau_1=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print colored ('Line 15425!','yellow')
						print
						quit()
					print
					print colored('2-CTR-gaussian already fitted!','yellow')
					print colored('Values from fits header','yellow')
					print colored(LINES[3][lines],'yellow')
					print
					try:
						CTRE_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CMC2')
						AMPL_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_AMC2')
						SGMA_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SMC2')
						FWHM_G_C_2    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FMC2')
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines-1])+'_CF01','yellow')
						print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
						print colored ('Or UnFix (fix_ctr_gau_2=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print colored ('Line 15468!','yellow')
						print
						quit()
					try:
						print
						print colored('1 PRE-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						CTRE_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM1')
						AMPL_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM1')
						SGMA_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM1')
						FWHM_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM1')
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!' ,'yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines])+'_CGM1','yellow')
						print colored ('Gotta Fit first before fixing a component (PRE)!','yellow')
						print colored ('Or UnFix (fix_pre_gau=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print colored ('Line 15503!','yellow')
						print
						quit()
					try:
						print
						print colored('2 PST-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						CTRE_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM2')
						AMPL_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM2')
						SGMA_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM2')
						FWHM_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM2')
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
						print colored ('Gotta Fit first before fixing a component (PST)!','yellow')
						print colored ('Or UnFix (fix_pst_gau=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print
						print colored ('Line 15530!','yellow')
						quit()
					try:
						print
						print colored('3- MDL-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						CTRE_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGM3')
						AMPL_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGM3')
						SGMA_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGM3')
						FWHM_G_C_MDL   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGM3')
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
						print colored ('Gotta Fit first before fixing a component (PST)!','yellow')
						print colored ('Or UnFix (fix_pst_gau=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print colored ('Line 15561!','yellow')
						print
						quit()			
				elif 'Dbl' not in LINES[3][lines] and fit_fnct=='gauss'   and mke_lne_fit == False:
					fit_typ = 'G'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print

					lmb_min_lim_line_ft = (LINES[0][lines]-LINES[8][lines]) - MSK_NTMS*LINES[1][lines]  
					lmb_max_lim_line_ft = (LINES[0][lines]+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
					lmb_min_lim_line    = (LINES[0][lines]-LINES[8][lines])*(1+z_glx_Ps) - 1.5*MSK_NTMS*LINES[2][lines] 
					lmb_max_lim_line    = (LINES[0][lines]+LINES[8][lines])*(1+z_glx_Ps) + 1.5*MSK_NTMS*LINES[2][lines]

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function

					ivl_fts_hdr = True
					########################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#########################
					if ivl_fts_hdr == True:
						try:
							L1_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF_0')        #LINES-1 Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP_0')        #LINES-2 Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF_0')        #LINES-7 Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_0),'yellow')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_0),'yellow')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_0),'yellow')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_0),'yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_0),'yellow')
							print
						except KeyError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
							print '*****'
							print
							quit()
						try:
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_0),'yellow')
							print colored('*****Success-1!******','magenta')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF_0',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_0),'yellow')
							print colored('*****Success-2!******','magenta')
							print
					else:
						pass
					try:
						AMPL_G_O   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGLO')

						CTRE_G_C   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGLC')
						AMPL_G_C   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGLC')
						SGMA_G_C   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGLC')
						FWHM_G_C   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGLC')
						EW_C       = Header_Get(specfile_glx,str(LINES[5][lines])+'_WGLC')
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines-2])+'_CF01','yellow')
						print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
						print colored ('Or UnFix (fix_ctr_gau_1=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print colored ('Line 15654!','yellow')
						print
						quit()
				elif 'Dbl' not in LINES[3][lines] and fit_fnct=='gaussM'  and mke_lne_fit == False:
					fit_typ = 'GM'
					#GAUSSIAN FIT
					MSK_NTMS=2.5
					print
					print colored('1D Gaussian Fit Mode Choosen: Lmfit (Offset)','cyan')
					print
					ivl_fts_hdr = True

					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################
					if ivl_fts_hdr == True:
						try:
							L1_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WF_0')        #LINES-1 Wdt-Fit  1GF-IntVal      WIDTH-FIT
							L2_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_WP_0')        #LINES-2 Wdt-Plt  1GF-IntVal      WIDTH-PLT
							L7_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CF_0')        #LINES-7 Ctr Fit Bnds  1GF-IntVal CTR-FIT-BNDS
							L8_0  = Header_Get(org_spc_fle,str(LINES[5][lines])+'_CO_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal CTR-OFFSET
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')       #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print
							print colored('Initial fit variables from fits header!','yellow')
							print colored('Headers:','yellow')
							print
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal      : ' + str(L1_0),'yellow')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal      : ' + str(L2_0),'yellow')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal : ' + str(L7_0),'yellow')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal : ' + str(L8_0),'yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal: ' + str(L10_0),'yellow')
							print colored('*****Success!******','magenta')
							print
						except ValueError:
							print
							print colored('Headers containing initial fit variables NOT found!','yellow')
							print colored('Run Spec_Stk_Stack first or use value from Lines_Dictionary.py with ivl_fts_hdr = False:','yellow')
							print colored('Headers:','yellow')
							print
							print colored(str(LINES[5][lines])+'_WF_0' + ' LINES-1 Wdt-Fit  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_WP_0' + ' LINES-2 Wdt-Plt  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_CF_0' + ' LINES-7 Ctr Fit Bnds  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_CO_0' + ' LINES-8 Ctr Fit Ofst  1GF-IntVal','yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
							print '*****'
							print
							quit()
						try:
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
							print colored('*****Success!******','magenta')
							print
						except KeyError:
							print
							print colored('Header NOT found!','yellow')
							print colored('Adding Header with default valuee 0.001:','yellow')
							print colored(str(LINES[5][lines])+'_AF_0' + ' LINES-10 Amp Fit Bnds  1GF-IntVal','yellow')
							Header_Get_Add(org_spc_fle,str(LINES[5][lines])+'_AF_0',0.001,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' LINES-10 Amp Fit Bnds  1GF-IntVal')
							L10_0 = Header_Get(org_spc_fle,str(LINES[5][lines])+'_AF_0')        #LINES-8 Ctr Fit Ofst  1GF-IntVal LNE_AMP_BNDS
							print colored('*****Success!******','magenta')
							print
						################FOR CASES WHERE KLINE WIIDHT ==0################
						if L1_0  == 0:
							L1_0  = 1#LINES[1][lines]
						else:
							pass
						################FOR CASES WHERE KLINE WIIDHT ==0################
					elif ivl_fts_hdr == False:
						L1_0  = LINES[1][lines]
						L2_0  = LINES[2][lines]
						L7_0  = LINES[7][lines]
						L8_0  = LINES[8][lines]
						L10_0 = LINES[10][lines]
						print
						print colored('Initial fit variables from Lines_Dictionary.py!','yellow')
						print					
					print str(LINES[5][lines])+'_WF_0' + ': ' + str(L1_0)
					print str(LINES[5][lines])+'_WP_0' + ': ' + str(L2_0)
					print str(LINES[5][lines])+'_CF_0' + ': ' + str(L7_0)
					print str(LINES[5][lines])+'_CO_0' + ': ' + str(L8_0)
					print str(LINES[5][lines])+'_AF_0' + ': ' + str(L10_0)
					#quit()

					#####################################################Getting Line Fitting Initial Values From Statcked Spectrum Fits Header#####################################################

					lmb_min_lim_line_ft = (LINES[0][lines]-LINES[8][lines]) - MSK_NTMS*LINES[1][lines] 
					lmb_max_lim_line_ft = (LINES[0][lines]+LINES[8][lines]) + MSK_NTMS*LINES[1][lines]
					lmb_min_lim_line    = (LINES[0][lines]-LINES[8][lines])*(1+z_glx_Ps) - 1.5*MSK_NTMS*LINES[2][lines]
					lmb_max_lim_line    = (LINES[0][lines]+LINES[8][lines])*(1+z_glx_Ps) + 1.5*MSK_NTMS*LINES[2][lines]

					mask_pl    = (lambda_glx >= lmb_min_lim_line)    & (lambda_glx <= lmb_max_lim_line)
					mask_ft    = (lambda_glx >= lmb_min_lim_line_ft) & (lambda_glx <= lmb_max_lim_line_ft)
					label_glx  = (specfile_glx.split('sep_as-',1)[1]).split('-stk',1)[0] #+ ' ' + stk_function

					print
					print colored('0 CTR-gaussian already fitted!','yellow')
					print colored('Values from fits header','yellow')
					try:
						CTRE_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLOM')
						AMPL_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALOM')
						FWHM_G_O    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLOM')

						CTRE_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_CLCM')
						AMPL_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_ALCM')
						SGMA_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_SLCM')
						FWHM_G_C    = Header_Get(specfile_glx,str(LINES[5][lines])+'_FLCM')
						EW_C        = Header_Get(specfile_glx,str(LINES[5][lines])+'_WLCM')
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
						print colored ('Gotta Fit first before fixing a component (CTR)!','yellow')
						print colored ('Or UnFix (fix_ctr_gau=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print colored ('Line 15782!','yellow')
						print
						quit()
					try:
						print
						print colored('1 PRE-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						CTRE_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGL1')
						AMPL_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGL1')
						SGMA_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGL1')
						FWHM_G_C_PRE   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGL1')

						x_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_XA1')
						y_a            = Header_Get(specfile_glx,str(LINES[5][lines])+'_YA1')
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!' ,'yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
						print colored ('Gotta Fit first before fixing a component (PRE)!','yellow')
						print colored ('Or UnFix (fix_pre_gau=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print
						print colored ('Line 9761!','yellow')
						quit()
					try:
						print
						print colored('2 PST-gaussian already fitted!','yellow')
						print colored('Values from fits header','yellow')
						CTRE_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_CGL2')
						AMPL_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_AGL2')
						SGMA_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_SGL2')
						FWHM_G_C_PST   = Header_Get(specfile_glx,str(LINES[5][lines])+'_FGL2')

						x_b         = Header_Get(specfile_glx,str(LINES[5][lines])+'_XA2')
						y_b         = Header_Get(specfile_glx,str(LINES[5][lines])+'_YA2')
						
					except KeyError:
						print
						print colored ('Gaussian Fit Parameters are NOT found on Fits Header!','yellow')
						print colored ('File  : ' + specfile_glx,'yellow')
						print colored ('Header: ' + str(LINES[5][lines])+'_CF0M','yellow')
						print colored ('Gotta Fit first before fixing a component (PST)!','yellow')
						print colored ('Or UnFix (fix_pst_gau=False) For Plotting!','yellow')
						print colored ('Quitting!','yellow')
						print
						print colored ('Line 15849!','yellow')
						quit()

					slope_line1 = (y_a-y_b)/(x_a-x_b)
					slope_line2 = (y_b-y_a)/(x_b-x_a)
					b1 = y_a - (slope_line1*x_a)
					b2 = y_b - (slope_line1*x_b)
					print
					print  colored('Computing Linear Area considering peak points:','yellow')
					print 'Point A: ',x_a,y_a
					print 'Point B: ',x_b,y_b
					print 'Slope: ',slope_line1
					print 'Slope: ',slope_line2
					print 'b: ',b1,b2
					print				
				else:
					pass
				########################################################LINE-FIT#######################################################
				AMPL_G_O=A_f2DG
				########################################################LINE-SNR#######################################################
				spc_stt_lmbd_i = kwargs.get('spc_stt_lmbd_i',1200)
				spc_stt_lmbd_f = kwargs.get('spc_stt_lmbd_f',1750)
				x_type         = kwargs.get('x_type','lambda')
				
				spc_val_stt    = Spectra_Cont_GetVal(specfile_glx,gcv_lmbd_i=spc_stt_lmbd_i,gcv_lmbd_f=spc_stt_lmbd_f,x_type=x_type)
				
				SPC_SNR_FUNC_UB_1   = SNR(inten_glx)
				SPC_SNR_FUNC_UB_2   = SNR(spc_val_stt[0])
				SPC_SNR_FUNC_UB_3   = SNR(inten_glx[mask_ft])
				SPC_SNR_FUNC_UB_4   = SNR(inten_glx[mask_pl])
				
				
				SPC_SNR_UB_1 = abs(AMPL_SNR)/SPC_SNR_FUNC_UB_1[2]
				SPC_SNR_UB_2 = abs(AMPL_SNR)/SPC_SNR_FUNC_UB_2[2]
				SPC_SNR_UB_3 = abs(AMPL_SNR)/SPC_SNR_FUNC_UB_3[2]
				SPC_SNR_UB_4 = abs(AMPL_SNR)/SPC_SNR_FUNC_UB_4[2]
				
				###### SNR out the line plot region #####
				if AMPL_SNR == 999999.99999:
					SPC_SNR_UB_1 = 999999.99999
					SPC_SNR_UB_2 = 999999.99999
					SPC_SNR_UB_3 = 999999.99999
					SPC_SNR_UB_4 = 999999.99999
					
					bin_size     = 999999.99999
					bin_number   = 999999.99999
					SPC_NSE_BN_1 = 999999.99999
					SPC_NSE_BN_2 = 999999.99999
					SPC_NSE_BN_3 = 999999.99999
					SPC_SNR_BN_1 = 999999.99999
					SPC_SNR_BN_2 = 999999.99999
					SPC_SNR_BN_3 = 999999.99999
				else:
					lne_lim_inf_sgm,lne_lim_sup_sgm = CTRE_SNR-SGMA_SNR,CTRE_SNR+SGMA_SNR
					lne_lim_inf_fwh,lne_lim_sup_fwh = CTRE_SNR-(sigma2fwhm(SGMA_SNR)/2),CTRE_SNR+(sigma2fwhm(SGMA_SNR)/2)
					chn_lim_ctr   = np.digitize(1*abs(CTRE_SNR),lambda_glx[mask_pl],right=False)
					chn_lim_inf   = np.digitize(1*abs(lne_lim_inf_fwh),lambda_glx[mask_pl],right=False)
					chn_lim_sup   = np.digitize(1*abs(lne_lim_sup_fwh),lambda_glx[mask_pl],right=True)

					#####################SIGMA-FIT-TOO-WIDE-THAT-FALLS-OUT-OF-FIR-REGION####################
					#######################IF-IT-HAPPENS-DIGITIZE-RETURNS-LEN(ARRAY)/0######################

					if (chn_lim_inf == len(lambda_glx[mask_pl])) or (chn_lim_sup == len(lambda_glx[mask_pl])) or (lambda_glx[mask_pl][0]>CTRE_SNR) or (CTRE_SNR>lambda_glx[mask_pl][1]):
						print
						print 'Out of Bounds!-1'
						print
						chn_lim_ctr = 2						
						chn_lim_inf = chn_lim_ctr - 1
						chn_lim_sup = chn_lim_ctr + 1
					elif (chn_lim_inf == -1) or (chn_lim_sup == 1):
						print
						print 'Out of Bounds!-2'
						print
						chn_lim_inf = chn_lim_ctr - 1
						chn_lim_sup = chn_lim_ctr + 1

					else:
						pass
					#####################SIGMA-FIT-TOO-WIDE-THAT-FALLS-OUT-OF-FIR-REGION####################

					indx_inf      = np.arange(0,chn_lim_inf)
					indx_sup      = np.arange(chn_lim_sup+1,len(lambda_glx[mask_pl]))
					indx_tot      = np.hstack((indx_inf,indx_sup))
										
					flx_tot       = inten_glx[mask_pl][indx_tot]#[0,1,2,3,4,5,6,7,8,9]
					bin_size      = len(lambda_glx[mask_pl][chn_lim_inf:chn_lim_sup+2])
					bin_number    = np.ceil(len(inten_glx[mask_pl][indx_tot])/float(bin_size))
					binned_arrays = split_list(flx_tot, wanted_parts=int(bin_number))
					
					NST_STS_AVG = []
					NST_STS_MED = []
					NST_STS_SUM = []
					NST_STS_STD = []
					for array_iteration in binned_arrays:
						NST_STS_AVG.append(bn.nanmean(array_iteration))
						NST_STS_MED.append(bn.nanmedian(array_iteration))
						NST_STS_SUM.append(bn.nansum(array_iteration))
						NST_STS_STD.append(bn.nanstd(array_iteration))
					
					if bin_number == 1:
						SPC_NSE_BN_1 = NST_STS_STD[0]
						SPC_NSE_BN_2 = NST_STS_STD[0]
						SPC_NSE_BN_3 = NST_STS_STD[0]
					else:
						SPC_NSE_BN_1 = bn.nanstd(NST_STS_AVG)
						SPC_NSE_BN_2 = bn.nanstd(NST_STS_MED)
						SPC_NSE_BN_3 = bn.nanstd(NST_STS_SUM)
					
					SPC_SNR_BN_1 = abs(AMPL_SNR)/SPC_NSE_BN_1
					SPC_SNR_BN_2 = abs(AMPL_SNR)/SPC_NSE_BN_2
					SPC_SNR_BN_3 = abs(AMPL_SNR)/SPC_NSE_BN_3
					
					if verbose == True:
						pass
					else:
						pass
				###### SNR out the line plot region #####


				Header_Add(specfile_glx,str(LINES[5][lines])+'_NSU1',float(SPC_SNR_FUNC_UB_1[2]) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE UBIN SPC REG ALL')
				Header_Add(specfile_glx,str(LINES[5][lines])+'_NSU2',float(SPC_SNR_FUNC_UB_2[2]) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE UBIN SPC REG CNT')
				Header_Add(specfile_glx,str(LINES[5][lines])+'_NSU3',float(SPC_SNR_FUNC_UB_3[2]) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE UBIN SPC REG FIT')
				Header_Add(specfile_glx,str(LINES[5][lines])+'_NSU4',float(SPC_SNR_FUNC_UB_4[2]) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE UBIN SPC REG PLT')

				Header_Add(specfile_glx,str(LINES[5][lines])+'_SNU1',float(SPC_SNR_UB_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR UBIN SPC REG ALL')
				Header_Add(specfile_glx,str(LINES[5][lines])+'_SNU2',float(SPC_SNR_UB_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR UBIN SPC REG CNT')
				Header_Add(specfile_glx,str(LINES[5][lines])+'_SNU3',float(SPC_SNR_UB_3) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR UBIN SPC REG FIT')
				Header_Add(specfile_glx,str(LINES[5][lines])+'_SNU4',float(SPC_SNR_UB_4) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR UBIN SPC REG PLT')

				Header_Add(specfile_glx,str(LINES[5][lines])+'_SNBS',float(bin_size)     ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR BINNED Size')
				Header_Add(specfile_glx,str(LINES[5][lines])+'_SNBN',float(bin_number)   ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR BINNED # Bins')

				try:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_NSB1',float(SPC_NSE_BN_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE BINNED AVG')
				except ValueError:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_NSB1',float(999999.99999) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE BINNED AVG')
				try:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_NSB2',float(SPC_NSE_BN_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE BINNED MED')
				except ValueError:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_NSB2',float(999999.99999) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE BINNED MED')
				try:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_NSB3',float(SPC_NSE_BN_3) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE BINNED SUM')				
				except ValueError:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_NSB3',float(999999.99999) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' NSE BINNED SUM')				

				try:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_SNB1',float(SPC_SNR_BN_1) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR BINNED AVG')
				except ValueError:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_SNB1',float(999999.99999) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR BINNED AVG')
				try:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_SNB2',float(SPC_SNR_BN_2) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR BINNED MED')
				except ValueError:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_SNB2',float(999999.99999) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR BINNED MED')
				try:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_SNB3',float(SPC_SNR_BN_3) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR BINNED SUM')				
				except ValueError:
					Header_Add(specfile_glx,str(LINES[5][lines])+'_SNB3',float(999999.99999) ,header_comment = str(LINES[3][lines]) + str(LINES[0][lines]) + ' SNR BINNED SUM')				

				########################################################LINE-SNR#######################################################
				
				if fit_fnct is not None:
					if AMPL_G_O <0: #popt_O[1]
						lambda_glx_lne_fit = np.arange(lmb_min_lim_line,lmb_max_lim_line,0.1)
						if fit_fnct=='gauss' and not 'DblF' in LINES[3][lines]:
							if pre_off_plt == True:
								try:
									plt.plot(lambda_glx_lne_fit,
												func_1D_Gaussian(lambda_glx_lne_fit,*popt_0),
												color='gray'  ,ls=':',lw=1.5,
												label = label_glx + ' g fit Org: ' + 
												'EW: ' + str(EW_0), 
												alpha=0.5)
									plt.plot(lambda_glx_lne_fit,
												func_1D_Gaussian_O(lambda_glx_lne_fit,*popt_O),
												color='gray'  ,ls='--',lw=1.5,
												label = label_glx + ' g fit Off: ' + 
												'EW: ' + str(EW_O) + ', '  
												'OFS: ' + str(np.round(OFST_G_O,2)), 
												alpha=0.5)
									pass
								except UnboundLocalError:
									plt.plot(lambda_glx_lne_fit,
												func_1D_Gaussian(lambda_glx_lne_fit,CTRE_G_0,AMPL_G_0,SGMA_G_0)     ,
												color='gray'  ,ls=':',lw=1.5,
												label = label_glx + ' g fit Org: ' + 
												'EW: ' + str(EW_0), 
												alpha=0.5)
									plt.plot(lambda_glx_lne_fit,
												func_1D_Gaussian_O(lambda_glx_lne_fit,CTRE_G_O,AMPL_G_O,SGMA_G_O,OFST_G_O)  ,
												color='gray'  ,ls='--',lw=1.5,
												label = label_glx + ' g fit Off: ' + 
												'EW: ' + str(EW_O) + ', '  
												'OFS: ' + str(np.round(OFST_G_O,2)), 
												alpha=0.5)
							elif pre_off_plt == False:
								pass
							try:
								plt.plot(lambda_glx_lne_fit,
										func_1D_Gaussian(lambda_glx_lne_fit,*popt_C)    ,
										color=colors[index]   ,ls='-',lw=1.5,
										label = label_glx + ' '+ 
										"\n" +
										'$\lambda_{f}$: ' + str(np.round(CTRE_G_C,2)) +
										"\n" +
										'$\sigma_{f}$: ' + str(np.round(SGMA_G_C,2)) +
										"\n" +
										'EW: ' + str(np.round(EW_C,2)) + ', ' +
										"\n" +
										'N : ' + str(stk_glx_nmb) +  ', ' +
										'SNR (M) : ' + str(np.round(SPC_SNR_BN_2,2))+ ', ' 
										#'$\chi^{2}$: '+str(chisqr_C)+ ', ' +
										#'$\chi^{2}_{r}$:' + str(redchi_C)
										,
										alpha=1.0)
							except UnboundLocalError:
								plt.plot(lambda_glx_lne_fit,
										func_1D_Gaussian(lambda_glx_lne_fit,CTRE_G_C,AMPL_G_C,SGMA_G_C)    ,
										color=colors[index]   ,ls='-',lw=1.5,
										label = label_glx + ' '+
										"\n" +
										'$\lambda_{f}$: ' + str(np.round(CTRE_G_C,2)) +
										"\n" +
										'$\sigma_{f}$: ' + str(np.round(SGMA_G_C,2)) +
										"\n" +
										'EW: ' + str(np.round(EW_C,2)) + ', ' +
										"\n" +
										'N : ' + str(stk_glx_nmb) +  ', ' +
										'SNR (M) : ' + str(np.round(SPC_SNR_BN_2,2))+ ', ' 
										#'$\chi^{2}$: '+str(chisqr_C)+ ', ' +
										#'$\chi^{2}_{r}$:' + str(redchi_C)
										,
										alpha=1.0)
							plt.step(lambda_glx[mask_pl], inten_glx[mask_pl], colors[index],where='mid',lw=1.0,alpha=0.6,linestyle='-',color='blue')#label=label_glx,
						elif fit_fnct=='lorentz' and not 'DblF' in LINES[3][lines]:
							plt.plot(lambda_glx_lne_fit,func_Lorentzian(lambda_glx_lne_fit,*popt)   ,
									color=colors[index]  ,ls='--',lw=1.5, 
									label = label_glx + ' l fit ' + 
									'EW: ' + str(EW))# + ', EWE: '+str(EWE))
						elif fit_fnct=='voigt' and not 'DblF' in LINES[3][lines]:
							plt.plot(lambda_glx_lne_fit,func_Voigt(lambda_glx_lne_fit,X_0=CTRE_V_0,ALPHA=ALPH_V_0,GAMMA=GAMA_V_0)        ,
									color=colors[index]  ,ls='--',lw=1.5, 
									label = label_glx + ' v fit $\lambda_{0}$='+ str(np.round(CTRE_V_0,2))+',' + 
									'EW: ' + str(EW_V_0) + ', ' +
									"\n" +
									'N : ' + str(stk_glx_nmb) +  ', ' +
									'SNR (M) : ' + str(np.round(SPC_SNR_BN_2,2))+ ', ' #+
									#'$\chi^{2}$: '+str(chisqr_V_0)+ ', ' +
									#'$\chi^{2}_{r}$:' + str(redchi_V_0)
									,
									alpha=1.0)
							pass
					else:
						pass
				else:
					pass
				########################################################LINE-FIT#######################################################

			plt.xlim([lmb_min_lim_line-2,lmb_max_lim_line+2])
			xmin, xmax = plt.xlim()
			plt.xlim((xmin,xmax))

			#plt.ylim([min_y,max_y])
			#ymin, ymax = plt.ylim()
			#plt.ylim((ymin,ymax))
			min_y, max_y = ax110.get_ylim()

			left, width    = lmb_min_lim_line, (lmb_max_lim_line - lmb_min_lim_line)
			bottom, height = min_y, (max_y - min_y)
			right          = lmb_max_lim_line + width
			top            = max_y
			
			plt.plot([lmb_min_lim_line, lmb_max_lim_line], [1, 1],
						color='black', lw=1, alpha=0.6)			
			plt.plot([LINES[0][lines]+L8_0,LINES[0][lines]+L8_0],[min_y, max_y], #LINES[8][lines]
						color='black', lw=1, alpha=0.6,ls=':',
						label=LINES[9][lines]+' ' + str(LINES[0][lines])+' $\mathrm{\mathrm{\AA}}$')
			#plt.plot([lmb_min_lim_line_ft,lmb_min_lim_line_ft], [min_y, max_y], 
						#color='black', lw=1, alpha=0.6,ls=':')
			#plt.plot([lmb_max_lim_line_ft,lmb_max_lim_line_ft], [min_y, max_y], 
						#color='black', lw=1, alpha=0.6,ls=':')


			lg2=plt.legend(loc=3,prop={'size':12})
			lg2.draw_frame(False)

			################################################################SAVE###########################################################

			if dest_dir != None:
				PLOTFILENAME1 = str(dest_dir) + plt_sufix_fnm+'-'+str(int(LINES[0][lines]))+'-'+str(LINES[4][lines])+ '-' + fit_type + '-'+fit_typ + '.pdf'
				PLOTFILENAME2 = str(dest_dir) + plt_sufix_fnm+'-'+str(int(LINES[0][lines]))+'-'+str(LINES[4][lines])+ '-' + fit_type + '-'+fit_typ + '.eps'

			elif dest_dir == None:
				PLOTFILENAME1 = ind_bst_plt + plt_sufix_fnm+'-'+str(int(LINES[0][lines]))+'-'+str(LINES[4][lines])+ '-' + fit_type + '-' + fit_typ + '.pdf'
				PLOTFILENAME2 = ind_bst_plt + plt_sufix_fnm+'-'+str(int(LINES[0][lines]))+'-'+str(LINES[4][lines])+ '-' + fit_type + '-' + fit_typ + '.eps'

			plt.savefig(PLOTFILENAME1)
			if verbose == True:
				print
				print colored('Generated Plot: ' + str(PLOTFILENAME1),'green')
			elif verbose ==False:
				pass
			if epssave == True:
				plt.savefig(PLOTFILENAME2, rasterized=True)
				#os.system('open Spectra.eps')
			elif epssave == False:
				pass
			if showplot == True:
				#os.system('open '+str(PLOTFILENAME))
				pass
			elif showplot == False:
				pass	
			plt.close('all')
			########################################################PLOT PER LINE########################################################
		elif LINES[0][lines]*(1+z_glx_Ps) > 6700 :
			break

def Select_Subsamples(slct_smp_cat,col_name,slice_itv,*args, **kwargs):
	bs_func     = kwargs.get('bs_func'    ,None)
	#Random_Vars = kwargs.get('Random_Vars',False)
	sel_pre_shf = kwargs.get('sel_pre_shf',True)
	sel_pre_cnt = kwargs.get('sel_pre_cnt',True)
	sel_pre_msk = kwargs.get('sel_pre_msk',False)
	spc_nse     = kwargs.get('spc_nse'    ,False)

	slc_smp     = kwargs.get('slc_smp', True)
	slc_int     = kwargs.get('slc_int', True)

	test_fg     = kwargs.get('test_fg',False)
	test_bg     = kwargs.get('test_bg',False)
	add_phi_crc_col = kwargs.get('add_phi_crc_col',False)
	print
	print colored('Slicing Sample to define sub-samples.','cyan')
	print colored('From table: '+slct_smp_cat,'green')
	print
	print colored('Pre-Processed spectra files properties : ','yellow')
	print colored('Spectra shifted                        : '+str(sel_pre_shf),'yellow')
	print colored('Spectra cont fitted                    : '+str(sel_pre_cnt),'yellow')
	print colored('Spectra masked                         : '+str(sel_pre_msk),'yellow')
	print colored('Add & Save PHI CRC Column              : '+str(add_phi_crc_col),'yellow')
	print

	if  '-BS-' in slct_smp_cat:
		cpy_fts_dir = fts_bst_lst 
		nme_fts_dir = cpy_fts_dir
		print
		print 'step: ',bs_func
		print
	elif '-BS_MST' in slct_smp_cat:
		cpy_fts_dir = ind_bst_lst
		nme_fts_dir = str_bst_stk
		print
		print 'step: ',bs_func
		print
	else:
		cpy_fts_dir = ind_stk_res
		nme_fts_dir = cpy_fts_dir

	if 'BS_MST' in slct_smp_cat:

		print
		print colored('BS-MASTER!','cyan')
		print colored (slct_smp_cat,'cyan')
		print			

		cat_parent_tbl_mas_slice = readtable_fg_bg_glx(slct_smp_cat,tbl_format_ipt,**kwargs)

		cat_parent_tbl_all_slice = cat_parent_tbl_mas_slice[0]

		spc_f_BS_slice       = cat_parent_tbl_mas_slice[1]
		spc_f_n_BS_slice     = spc_f_BS_slice


	else:
		cat_parent_tbl_mas_slice = aptbl.Table.read(slct_smp_cat, format=tbl_format_ipt)

		cat_parent_tbl_all_slice = cat_parent_tbl_mas_slice #cat_parent_tbl_mas_slice[0]


		id_F_slice               = cat_parent_tbl_mas_slice['id_F']
		z_F_slice                = cat_parent_tbl_mas_slice['z_F']
		zf_F_slice               = cat_parent_tbl_mas_slice['zf_F']
		id_B_slice               = cat_parent_tbl_mas_slice['id_B']
		z_B_slice                = cat_parent_tbl_mas_slice['z_B']
		zf_B_slice               = cat_parent_tbl_mas_slice['zf_B']
		DELTAZ_slice             = cat_parent_tbl_mas_slice['DELTAZ']
		SEP_arcsec_slice         = cat_parent_tbl_mas_slice['SEP_arcsec']
		arcsec_kpc_slice         = cat_parent_tbl_mas_slice['arcsec/kpc']
		SEP_kpc_slice            = cat_parent_tbl_mas_slice['SEP_kpc']
		spc_f_F_slice            = cat_parent_tbl_mas_slice['spc_f_F']
		spc_f_n_F_slice          = cat_parent_tbl_mas_slice['spc_f_n_F']
		spc_f_B_slice            = cat_parent_tbl_mas_slice['spc_f_B']
		spc_f_n_B_slice          = cat_parent_tbl_mas_slice['spc_f_n_B']

		if '-PRP' in slct_smp_cat:
			print
			print colored('PRP!','cyan')
			print colored (slct_smp_cat,'cyan')
			print			
			mass_B_slice         = cat_parent_tbl_mas_slice['mass_B']
			Age_B_slice          = cat_parent_tbl_mas_slice['Age_B']
			SFR_B_slice          = cat_parent_tbl_mas_slice['SFR_B']
			sSFR_B_slice         = cat_parent_tbl_mas_slice['sSFR_B']
			Lnuv_B_slice         = cat_parent_tbl_mas_slice['Lnuv_B']
			mass_B_slice         = cat_parent_tbl_mas_slice['mass_F']
			Age_F_slice          = cat_parent_tbl_mas_slice['Age_F']
			SFR_F_slice          = cat_parent_tbl_mas_slice['SFR_F']
			sSFR_F_slice         = cat_parent_tbl_mas_slice['sSFR_F']
			Lnuv_F_slice         = cat_parent_tbl_mas_slice['Lnuv_F']
			magi_F_slice         = cat_parent_tbl_mas_slice['magi_F']
		elif '_PRP_MRP' in slct_smp_cat:
			print
			print colored('Morpho!','cyan')
			print colored (slct_smp_cat,'cyan')
			print
			mass_B_slice         = cat_parent_tbl_mas_slice['mass_B']
			Age_B_slice          = cat_parent_tbl_mas_slice['Age_B']
			SFR_B_slice          = cat_parent_tbl_mas_slice['SFR_B']
			sSFR_B_slice         = cat_parent_tbl_mas_slice['sSFR_B']
			Lnuv_B_slice         = cat_parent_tbl_mas_slice['Lnuv_B']
			mass_B_slice         = cat_parent_tbl_mas_slice['mass_F']
			Age_F_slice          = cat_parent_tbl_mas_slice['Age_F']
			SFR_F_slice          = cat_parent_tbl_mas_slice['SFR_F']
			sSFR_F_slice         = cat_parent_tbl_mas_slice['sSFR_F']
			Lnuv_F_slice         = cat_parent_tbl_mas_slice['Lnuv_F']
			magi_F_slice         = cat_parent_tbl_mas_slice['magi_F']
			PHI_slice            = abs(cat_parent_tbl_mas_slice['PHI'])
			re_F_slice           = cat_parent_tbl_mas_slice['re_F_kpc']
			n_F_slice            = cat_parent_tbl_mas_slice['n_F']
			q_F_slice            = cat_parent_tbl_mas_slice['q_F']
			PHI_slice_new=[]
			for j,angle_item in enumerate(PHI_slice):
				if angle_item>90:
					PHI_slice_new.append(180-angle_item)
				else:
					PHI_slice_new.append(angle_item)
			PHI_slice = np.asarray(PHI_slice_new)
		else:
			pass

	tbls_fnms                  = []
	subsample_stamps_fni       = []
	subsample_stamps_noise_fni = []

	subsample_crval1s = []
	subsample_cdelt1s = []
	subsample_crvalns = []
	subsample_magt_is = []

	subsample_spctype = []

	if   sel_pre_shf == True and sel_pre_msk == True and sel_pre_cnt == True and test_fg == True:
		spc_ifn_ext = '-c-m-s2-'
		print
		print colored('Using spectra files with cont fit, masked ext : '+ spc_ifn_ext,'yellow')
	elif sel_pre_shf == True and sel_pre_msk == True and sel_pre_cnt == True and test_fg == False:
		spc_ifn_ext = '-c-s2-'
		print
		print colored('Using spectra files no correction with ext    : '+ spc_ifn_ext,'yellow')

	elif sel_pre_shf == True and sel_pre_msk == True and sel_pre_cnt == False and test_fg == True:
		spc_ifn_ext = '-m-s2-'
		print
		print colored('Using spectra files with ext masked           : '+ spc_ifn_ext,'yellow')

	elif sel_pre_shf == True and sel_pre_msk == True and sel_pre_cnt == False and test_fg == False:
		spc_ifn_ext = '-s2-'
		print
		print colored('Using spectra files no correction with ext    : '+ spc_ifn_ext,'yellow')
	elif sel_pre_shf == True and sel_pre_msk == False and sel_pre_cnt == True and (test_fg == False or test_fg == True):
		spc_ifn_ext = '-c-s2-'
		print
		print colored('Using spectra files with ext cont fit         : '+ spc_ifn_ext,'yellow')

	elif sel_pre_shf == True and sel_pre_msk == False and sel_pre_cnt == False and (test_fg == False or test_fg == True):
		spc_ifn_ext   = '-s2-'
		print
		print colored('Using spectra files no correction with ext    : '+ spc_ifn_ext,'yellow')
	else:
		pass
	
	if slc_smp == True:
		if col_name == 'redshift_fg':
			col_slice = z_F_slice
			slicedby  = 'z_F'
			print 'Filtering table by: ', slicedby
		elif col_name == 'redshift_fg_flag':
			col_slice = zf_F_slice
			slicedby  = 'zf_F'
			print 'Filtering table by: ', slicedby
		elif col_name == 'redshift_bk':
			col_slice = z_B_slice
			slicedby  = 'z_B'
			print 'Filtering table by: ', slicedby
		elif col_name == 'redshift_bk_flag':
			col_slice = zf_B_slice
			slicedby  = 'zf_B'
			print 'Filtering table by: ', slicedby
		elif col_name == 'deltaz':
			col_slice = DELTAZ_slice
			slicedby  = 'DELTAZ'
			print 'Filtering table by: ', slicedby
		elif col_name == 'sep_as':
			col_slice = SEP_arcsec_slice
			slicedby  = 'sep_as'
			print 'Filtering table by: ', slicedby
		elif col_name == 'sep_kpc':
			col_slice = SEP_kpc_slice
			slicedby  = 'sep_kpc'
			print 'Filtering table by: ', slicedby

		elif col_name == 'mass_bg':
			col_slice = mass_B_slice
			slicedby  = 'mass_B'
			print 'Filtering table by: ', slicedby
			
		elif col_name == 'Age_bg':
			col_slice =Age_B_slice
			slicedby  = 'Age_B'
			print 'Filtering table by: ', slicedby

		elif col_name == 'SFR_bg':
			col_slice =SFR_B_slice
			slicedby  = 'SFR_B'
			print 'Filtering table by: ', slicedby

		elif col_name == 'sSFR_bg':
			col_slice = sSFR_B_slice
			slicedby  = 'sSFR_B'
			print 'Filtering table by: ', slicedby

		elif col_name == 'Lnuv_bg':
			col_slice = Lnuv_B_slice
			slicedby  = 'Lnuv_B'
			print 'Filtering table by: ', slicedby

		elif col_name == 'mass_fg':
			col_slice = mass_B_slice
			slicedby  = 'mass_F'
			print 'Filtering table by: ', slicedby
			
		elif col_name == 'Age_fg':
			col_slice =Age_F_slice
			slicedby  = 'Age_F'
			print 'Filtering table by: ', slicedby

		elif col_name == 'SFR_fg':
			col_slice =SFR_F_slice
			slicedby  = 'SFR_F'
			print 'Filtering table by: ', slicedby

		elif col_name == 'sSFR_fg':
			col_slice = sSFR_F_slice
			slicedby  = 'sSFR_F'
			print 'Filtering table by: ', slicedby

		elif col_name == 'Lnuv_fg':
			col_slice = Lnuv_F_slice
			slicedby  = 'Lnuv_F'
			print 'Filtering table by: ', slicedby
		elif col_name == 'magi_fg':
			col_slice = magi_F_slice
			slicedby  = 'magi_F'
			print 'Filtering table by: ', slicedby
		###################MORPHO###################			
		elif col_name == 'phi':
			col_slice = PHI_slice
			slicedby  = 'PHI'
			print 'Filtering table by: ', slicedby
		elif col_name == 'icl_fg':
			col_slice = q_F_slice
			slicedby  = 'q_F'
			print 'Filtering table by: ', slicedby
		elif col_name == 'n_sersic_fg':
			col_slice = n_F_slice
			slicedby  = 'n_F'
			print 'Filtering table by: ', slicedby
		elif col_name == 'r_eff_fg':
			col_slice = re_F_slice
			slicedby  = 're_F'
			print 'Filtering table by: ', slicedby
		else:
			print colored('Such col name does not exist! ' + col_name,'yellow')
			quit()
		###################MORPHO###################
		print 
		print colored('Values used for slicing: '+str(slice_itv),'yellow')
		for previous, item, nxt in Prev_Next(slice_itv):
			slc_int_bin = slice_itv[1]-slice_itv[0]
			if   item <= slice_itv[-1]:
			   if slc_int == True:
				print
				print 'Objects in the range: ',item,'<',col_name,'<=',nxt
				print
				subsampl_tbl_fn   = str(slct_smp_cat.split('.'+tbl_format_ipt)[0]) + '-ss-' + str(slicedby) +'-' + str(item)         + '-' + str(nxt)           + tbl_ext_opt
				subsampl_tbl_fn_A = str(slct_smp_cat.split('.'+tbl_format_ipt)[0]) + '-ss-' + str(slicedby) +'-' + str(slice_itv[0]) + '-' + str(slice_itv[-1]) + tbl_ext_opt
				subsample = np.where((item < col_slice) & (col_slice <= nxt))
				subsample = np.ravel(subsample)
			   elif slc_int == False:#:
				print
				print colored('Objects with values equal to: '+str(item),'cyan')
				subsampl_tbl_fn   = str(slct_smp_cat.split('.'+tbl_format_ipt)[0]) + '-ss-' + str(slicedby) + '-' + str(item)         + tbl_ext_opt
				subsampl_tbl_fn_A = str(slct_smp_cat.split('.'+tbl_format_ipt)[0]) + '-ss-' + str(slicedby) + '-' + str(slice_itv[0]) + '-' + str(slice_itv[-1]) + tbl_ext_opt
				subsample = np.where(item==col_slice)
				subsample = np.ravel(subsample)
			   if len(cat_parent_tbl_all_slice[subsample])!=0:
				print 'Non empty'
				print
				print 'Number of objects in this range: ',len(cat_parent_tbl_all_slice[subsample]),(len(subsample))

				rtsubsmpl = cat_parent_tbl_all_slice[subsample]
				rtsubsmpl.write(subsampl_tbl_fn,  format=tbl_format_opt, overwrite=True)

				print 'Results containing object list in table: '
				print subsampl_tbl_fn
				print colored('Selecting files ...','yellow')

				tbls_fnms.append(str(subsampl_tbl_fn))
				subsample_spctype.append(str(spc_ifn_ext))

				if  (('P_Fg_ALL' in slct_smp_cat) or ('P_Fg_ECDFS' in slct_smp_cat) or ('P_Fg_COSMOS' in slct_smp_cat) or ('P_Fg_VVDS2H' in slct_smp_cat) ):
					#[ os.system('cp ' + fg_sbdir[0] + str(spc_f_F_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf')           + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_F_slice))]
					[ os.system('cp ' + fg_sbdir[0] + str(spc_f_F_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf')           + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_F_slice))]
					if spc_nse == True:
						[ os.system('cp ' + fg_sbdir[0] + str(spc_f_n_F_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf')           + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_n_F_slice))]
					elif spc_nse == False:
						pass

					name       = [ cpy_fts_dir + str(spc_f_F_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf')           + '.fits' for j in subsample]
					name_noise = [ cpy_fts_dir + str(spc_f_n_F_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf')           + '.fits' for j in subsample]
				
				elif (('P_Bg_ALL' in slct_smp_cat) or ('P_Bg_ECDFS' in slct_smp_cat) or ('P_Bg_COSMOS' in slct_smp_cat) or ('P_Bg_VVDS2H' in slct_smp_cat) ) and test_bg == False:
					#[ os.system('cp ' + bg_sbdir[0] + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str(id_F_slice[j])  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_B_slice))]
					[ os.system('cp ' + bg_sbdir[0] + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str(id_F_slice[j])  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_B_slice))]
					if spc_nse == True:
						[ os.system('cp ' + bg_sbdir[0] + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str(id_F_slice[j])  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_n_B_slice))]
					elif spc_nse == False:
						pass
					
					name       = [ cpy_fts_dir + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str(id_F_slice[j]) + '.fits' for j in subsample]
					name_noise = [ cpy_fts_dir + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str(id_F_slice[j]) + '.fits' for j in subsample]


				elif (('P_Bg_ALL' in slct_smp_cat) or ('P_Bg_ECDFS' in slct_smp_cat) or ('P_Bg_COSMOS' in slct_smp_cat) or ('P_Bg_VVDS2H' in slct_smp_cat) ) and test_bg == True:
					print
					print colored('Warning: Background galaxies check mode selected!','green')
					#[ os.system('cp ' + bg_sbdir[0] + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf')           + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_B_slice))]
					[ os.system('cp ' + bg_sbdir[0] + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf')           + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_B_slice))]
					if spc_nse == True:
						[ os.system('cp ' + bg_sbdir[0] + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf')           + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_n_B_slice))]
					elif spc_nse == False:
						pass
					
					name       = [ cpy_fts_dir + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf')           + '.fits' for j in subsample]
					name_noise = [ cpy_fts_dir + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf')           + '.fits' for j in subsample]

				crval1s = [Header_Get(specfile,'CRVAL1') for specfile in name ]
				cdelt1s = [Header_Get(specfile,'CDELT1') for specfile in name ]
				crvalns = [get_last_lambda_item(specfile)     for specfile in name ]

				print
				print colored('Selected files: ' + str((len(name))) + ', ' + str(len(name_noise)) + ' like: ','yellow')
				print colored(name[0],'yellow')
				print colored(name_noise[0],'yellow')

				subsample_stamps_fni.append(name)
				subsample_stamps_noise_fni.append(name_noise)

				subsample_crval1s.append(crval1s)
				subsample_cdelt1s.append(cdelt1s)
				subsample_crvalns.append(crvalns)

			   elif len(cat_parent_tbl_all_slice[subsample])==0:
				print
				print 'Empty'
				print 'Number of objects in this range: ',len(cat_parent_tbl_all_slice[subsample])
			else :
				break

		#Def a new function for stacking tables
		template = readtable_fg_bg_glx(tbls_fnms[-1],tbl_format_opt)[0] 
		template.remove_rows(slice(0,-1))
		template.remove_row(0)	
		for previous, item, nxt in Prev_Next(tbls_fnms):
			if   item != slice_itv[-1]:
				a = readtable_fg_bg_glx(item,tbl_format_opt)[0]
				template = aptbl.vstack([template,a])
			else:
				break
		template.write(subsampl_tbl_fn_A,  format=tbl_format_opt, overwrite=True)
		#End a new function for stacking tables

		tbls_fnms.append(str(subsampl_tbl_fn_A))
		subsample_stamps_fni.append(list(chain.from_iterable(subsample_stamps_fni)))
		subsample_stamps_noise_fni.append(list(chain.from_iterable(subsample_stamps_noise_fni)))
		subsample_crval1s.append(list(chain.from_iterable(subsample_crval1s)))
		subsample_cdelt1s.append(list(chain.from_iterable(subsample_cdelt1s)))
		subsample_crvalns.append(list(chain.from_iterable(subsample_crvalns)))
		subsample_spctype.append(str(spc_ifn_ext))	

		print
		print colored('Number of subsamples created according to the '+ str(col_name) + ' criteria selection: ' + str(len(tbls_fnms)) + ': '+str(slice_itv),'cyan')
		print
		print colored('The results are contained in: ','green')
		print "\n".join([tblnm for tblnm in tbls_fnms])
		print

	elif slc_smp == False:
		print
		print colored('No slicing the sample, just reading the table for Stacking','yellow')
		print
		subsampl_tbl_fn   = slct_smp_cat
		tbls_fnms.append(str(subsampl_tbl_fn))
		subsample_spctype.append(str(spc_ifn_ext))
		print 
		print slct_smp_cat		

		if  'P_Fg_' + CAT_PARENT in slct_smp_cat  and ('BS_MST' not in slct_smp_cat) and ('-BS-' not in slct_smp_cat) :
			[ os.system('cp ' + fg_sbdir[0] + str(spc_f_F_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf') + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_F_slice))]
			if spc_nse == True:
				[ os.system('cp ' + fg_sbdir[0] + str(spc_f_n_F_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf') + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_n_F_slice))]
			elif spc_nse == False:
				pass
		
			name       = [ cpy_fts_dir + str(spc_f_F_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf') + '.fits' for j in range(len(spc_f_F_slice))]
			name_noise = [ cpy_fts_dir + str(spc_f_n_F_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf') + '.fits' for j in range(len(spc_f_n_F_slice))]

		elif 'P_Bg_' + CAT_PARENT in slct_smp_cat and ( 'BS_MST' not in slct_smp_cat)  and ('-BS-' not in slct_smp_cat) and test_bg == False:
			[ os.system('cp ' + bg_sbdir[0] + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str(id_F_slice[j])  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_B_slice))]
			if spc_nse == True:
				[ os.system('cp ' + bg_sbdir[0] + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str(id_F_slice[j])  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_n_B_slice))]
			elif spc_nse == False:
				pass
		
			name       = [ cpy_fts_dir + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str(id_F_slice[j]) + '.fits' for j in range(len(spc_f_B_slice))]
			name_noise = [ cpy_fts_dir + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str(id_F_slice[j]) + '.fits' for j in range(len(spc_f_n_B_slice))]

		elif 'P_Bg_' + CAT_PARENT in slct_smp_cat and ( 'BS_MST' not in slct_smp_cat)  and ('-BS-' not in slct_smp_cat) and test_bg == True:
			print
			print colored('Warning: Background galaxies check mode selected! (Background)','green')
			[ os.system('cp ' + bg_sbdir[0] + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf') + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_B_slice))]
			if spc_nse == True:
				[ os.system('cp ' + bg_sbdir[0] + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf') + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_n_B_slice))]
			elif spc_nse == False:
				pass
		
			name       = [ cpy_fts_dir + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf') + '.fits' for j in range(len(spc_f_B_slice))]
			name_noise = [ cpy_fts_dir + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf') + '.fits' for j in range(len(spc_f_n_B_slice))]

		elif (('P_Bg_' + CAT_PARENT in slct_smp_cat)) and ('-BS-' in slct_smp_cat):
			print
			print colored('Warning: Bootstrap galaxies mode selected! (Background)','green')
			[ os.system('cp ' + bg_sbdir[0] + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str(id_F_slice[j])  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_B_slice))]
			if spc_nse == True:
				[ os.system('cp ' + bg_sbdir[0] + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str(id_F_slice[j])  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_n_B_slice))]
			elif spc_nse == False:
				pass
		
			name       = [ nme_fts_dir + str(spc_f_B_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str(id_F_slice[j]) + '.fits' for j in range(len(spc_f_B_slice))]
			name_noise = [ nme_fts_dir + str(spc_f_n_B_slice[j].split('.fits',1)[0]) + '-s2-'      + str(id_F_slice[j]) + '.fits' for j in range(len(spc_f_n_B_slice))]

		elif (('P_Bg_' + CAT_PARENT in slct_smp_cat)) and ('BS_MST' in slct_smp_cat):
			print
			print colored('Warning: Bootstrap galaxies mode selected! (Background-MST)','green')

			name       = [ nme_fts_dir + str(spc_f_BS_slice[j]) for j in range(len(spc_f_BS_slice))]
			name_noise = name

		elif (('P_Fg_' + CAT_PARENT in slct_smp_cat)) and ('-BS-' in slct_smp_cat):
			print
			print colored('Warning: Bootstrap galaxies mode selected! (Foreground)','green')
			[ os.system('cp ' + fg_sbdir[0] + str(spc_f_F_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf')  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_F_slice))]
			if spc_nse == True:
				[ os.system('cp ' + fg_sbdir[0] + str(spc_f_n_F_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf')  + '.fits ' + cpy_fts_dir) for j in range(len(spc_f_n_F_slice))]
			elif spc_nse == False:
				pass
		
			name       = [ nme_fts_dir + str(spc_f_F_slice[j].split('.fits',1)[0])   + spc_ifn_ext + str('rf') + '.fits' for j in range(len(spc_f_F_slice))]
			name_noise = [ nme_fts_dir + str(spc_f_n_F_slice[j].split('.fits',1)[0]) + '-s2-'      + str('rf') + '.fits' for j in range(len(spc_f_n_F_slice))]

		elif (('P_Fg_' + CAT_PARENT in slct_smp_cat)) and ('BS_MST' in slct_smp_cat):
			print
			print colored('Warning: Bootstrap galaxies mode selected! (Foreground-MST)','green')

			name       = [ nme_fts_dir + str(spc_f_BS_slice[j]) for j in range(len(spc_f_BS_slice))]
			name_noise = name

		crval1s = [Header_Get(specfile,'CRVAL1')  for specfile in name ]
		cdelt1s = [Header_Get(specfile,'CDELT1')  for specfile in name ]
		crvalns = [get_last_lambda_item(specfile) for specfile in name ]
		magt_is = [Header_Get(specfile,'MAG_I')   for specfile in name ]

		print
		print colored('Selecting: ' + str((len(name))) + ', ' + str(len(name_noise)) + ' files like: ','yellow')
		print colored(name[0],'yellow')
		print colored(name_noise[0],'yellow')

		subsample_stamps_fni.append(name)
		subsample_stamps_noise_fni.append(name_noise)

		subsample_crval1s.append(crval1s)
		subsample_cdelt1s.append(cdelt1s)
		subsample_crvalns.append(crvalns)
		subsample_magt_is.append(magt_is)
	
	return tbls_fnms,subsample_stamps_fni,subsample_stamps_noise_fni,subsample_crval1s,subsample_cdelt1s,subsample_crvalns,subsample_spctype

def Interpolating_Spectra(spc_tbl_RS,spc_ifn_RS,new_lambda_bin,new_wavelength_array,*args,**kwargs):
	SS_indx         = kwargs.get('SS_indx'         ,None)

	sel_pre_shf     = kwargs.get('sel_pre_shf'     ,True)
	sel_pre_cnt     = kwargs.get('sel_pre_cnt'     ,True)
	sel_pre_msk     = kwargs.get('sel_pre_msk'     ,False)

	pre_cnt         = kwargs.get('pre_cnt'         ,False)
	pre_cnt_typ     = kwargs.get('pre_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pre_cnt_lns     = kwargs.get('pre_cnt_lns'     ,'*')       # Image lines to be fit
	pre_cnt_fnc     = kwargs.get('pre_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pre_cnt_ord     = kwargs.get('pre_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pre_cnt_ovr     = kwargs.get('pre_cnt_ovr'     ,'yes')     # Override previous norm spec
	pre_cnt_rpl     = kwargs.get('pre_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pre_cnt_lrj     = kwargs.get('pre_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pre_cnt_hrj     = kwargs.get('pre_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pre     = kwargs.get('smt_spc_pre'     ,True)
	smt_shp_pre     = kwargs.get('smt_shp_pre'     ,'gaussian')
	smt_sze_pre     = kwargs.get('smt_sze_pre'     ,1)

	pre_msk         = kwargs.get('pre_msk'         ,False)
	pre_msk_typ     = kwargs.get('pre_msk_typ'     ,'NaN')
	pre_msk_cte_val = kwargs.get('pre_msk_cte_val' ,1)
	pre_msk_abs_lne = kwargs.get('pre_msk_abs_lne' ,False)
	pre_msk_rgn     = kwargs.get('pre_msk_rgn' ,False)
	pre_msk_min     = kwargs.get('pre_msk_min' ,500)
	pre_msk_max     = kwargs.get('pre_msk_max' ,1210)

	sig_clp         = kwargs.get('sig_clp'         ,False)
	sig_cut         = kwargs.get('sig_cut'         ,3)
	sig_fct         = kwargs.get('sig_fct'         ,mean)
	sig_fll         = kwargs.get('sig_fll'         ,np.nan)

	wgt_typ         = kwargs.get('wgt_typ'         ,None)
	get_cont_flux   = kwargs.get('get_cont_flux'   ,True)
	gcv_lmbd_i      = kwargs.get('gcv_lmbd_i'      ,1430)
	gcv_lmbd_f      = kwargs.get('gcv_lmbd_f'      ,1480)

	spc_nse         = kwargs.get('spc_nse'         ,False)

	pst_cnt         = kwargs.get('pst_cnt'         ,True)
	pst_cnt_typ     = kwargs.get('pst_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pst_cnt_lns     = kwargs.get('pst_cnt_lns'     ,'*')       # Image lines to be fit
	pst_cnt_fnc     = kwargs.get('pst_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pst_cnt_ord     = kwargs.get('pst_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pst_cnt_ovr     = kwargs.get('pst_cnt_ovr'     ,'yes')     # Override previous norm spec
	pst_cnt_rpl     = kwargs.get('pst_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pst_cnt_lrj     = kwargs.get('pst_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pst_cnt_hrj     = kwargs.get('pst_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pst     = kwargs.get('smt_spc_pst'     ,True)
	smt_shp_pst     = kwargs.get('smt_shp_pst'     ,'gaussian')
	smt_sze_pst     = kwargs.get('smt_sze_pst'     ,1)

	if pre_cnt == True and not ('-BS-' in spc_ifn_RS):
		print
		print '1'
		print
		spc_ifn_RS_cnt_M = Spectra_Cont_IRAF(spc_ifn_RS,lst_stk_res + 'log_cont-pre-cnt',
						Cont_type_IRAF     = pre_cnt_typ,Cont_lines_IRAF    = pre_cnt_lns,
						Cont_funct_IRAF    = pre_cnt_fnc,Cont_order_IRAF    = pre_cnt_ord,
						Cont_override_IRAF = pre_cnt_ovr,Cont_replace_IRAF  = pre_cnt_rpl,
						Cont_low_rej_IRAF  = pre_cnt_lrj,Cont_high_rej_IRAF = pre_cnt_hrj)
		spc_ifn_RS_cnt = spc_ifn_RS_cnt_M[0]

		Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'ID_0'  ,header_comment='ID')
		Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'ID_REF',header_comment='ID Reference')
		Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'Z_0'   ,header_comment='Redshift')
		Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'Z_REF' ,header_comment='Redshift Reference')
		Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'MAG_I' ,header_comment='i-band magnitude')
		Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'h_s_0')
		Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'h_s_c')
		Header_History_Step(spc_ifn_RS,spc_ifn_RS_cnt)
		loop_back_step = int(Header_Get(spc_ifn_RS_cnt,'h_s_c'))
		while (loop_back_step > 0):
			loop_back_step = int(loop_back_step - 1)
			loop_back_step_header = 'h_s_' + str(loop_back_step)
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,loop_back_step_header)
		else:
			pass		
		if 	sel_pre_cnt == True and not 'noise' in spc_ifn_RS and '-BS_MST' not in spc_tbl_RS:
			print
			print '1-a'
			print
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_TYP',header_comment='Continuum IRAF type')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_FIT',header_comment='Continuum IRAF function')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_ORD',header_comment='Continuum IRAF order')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'MAG_I'  ,header_comment='i-band magnitude')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'SHT_FNS',header_comment='Spec file post-shift')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'SHT_FN0',header_comment='Spec file post-shift')
		elif sel_pre_cnt == True and not 'noise' in spc_ifn_RS and '-BS_MST' in spc_tbl_RS:
			print
			print '1-b'
			print
			Header_Get_Add(spc_ifn_RS_cnt,'CNT_TYP',cont_typ     ,header_comment='Continuum IRAF type')
			Header_Get_Add(spc_ifn_RS_cnt,'CNT_FIT',cont_funct   ,header_comment='Continuum IRAF function')
			Header_Get_Add(spc_ifn_RS_cnt,'CNT_ORD',cont_order   ,header_comment='Continuum IRAF order')
			Header_Get_Add(spc_ifn_RS_cnt,'CNT_RPC',cont_replace ,header_comment='Continuum IRAF replacement')
			Header_Get_Add(spc_ifn_RS_cnt,'CNT_LRJ',cont_low_rej ,header_comment='Continuum IRAF low sigma rejection')
			Header_Get_Add(spc_ifn_RS_cnt,'CNT_HRJ',cont_high_rej,header_comment='Continuum IRAF high sigma rejection')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'MAG_I',header_comment='i-band magnitude')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
			Header_Copy(spc_ifn_RS_cnt,spc_ifn_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
		else:
			pass
		spc_ifn_RS    = spc_ifn_RS_cnt
	elif pre_cnt == False:
		pass

	if smt_spc_pre == True and not ('-BS-' in spc_ifn_RS):
		spc_ifn_RS_smt = Spectra_Smooth(spc_ifn_RS,smt_shp_pre,smt_sze_pre,*args,**kwargs)
		Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'ID_0'  ,header_comment='ID')
		Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'ID_REF',header_comment='ID Reference')
		Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'Z_0'   ,header_comment='Redshift')
		Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'Z_REF' ,header_comment='Redshift Reference')
		Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'MAG_I' ,header_comment='i-band magnitude')
		Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'h_s_0')
		Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'h_s_c')
		Header_History_Step(spc_ifn_RS,spc_ifn_RS_smt)
		loop_back_step = int(Header_Get(spc_ifn_RS_smt,'h_s_c'))
		while (loop_back_step > 0):
			loop_back_step = int(loop_back_step - 1)
			loop_back_step_header = 'h_s_' + str(loop_back_step)
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,loop_back_step_header)
		else:
			pass
		if 	smt_spc_pre == True and 'noise' not in spc_ifn_RS and '-BS_MST' not in spc_tbl_RS:
			try:
				Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_TYP',header_comment='Continuum IRAF type')
				Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_FIT',header_comment='Continuum IRAF function')
				Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_ORD',header_comment='Continuum IRAF order')
				Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
				Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
				Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
			except KeyError: 
				pass
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'MAG_I'  ,header_comment='i-band magnitude'  )
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'SHT_FNS',header_comment='Spec file post-shift')
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'SHT_FN0',header_comment='Spec file post-shift')
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		elif smt_spc_pre == True and 'noise' not in spc_ifn_RS and '-BS_MST' in spc_tbl_RS:
			try:
				Header_Get_Add(spc_ifn_RS_smt,'CNT_TYP',cont_typ     ,header_comment='Continuum IRAF type')
				Header_Get_Add(spc_ifn_RS_smt,'CNT_FIT',cont_funct   ,header_comment='Continuum IRAF function')
				Header_Get_Add(spc_ifn_RS_smt,'CNT_ORD',cont_order   ,header_comment='Continuum IRAF order')
				Header_Get_Add(spc_ifn_RS_smt,'CNT_RPC',cont_replace ,header_comment='Continuum IRAF replacement')
				Header_Get_Add(spc_ifn_RS_smt,'CNT_LRJ',cont_low_rej ,header_comment='Continuum IRAF low sigma rejection')
				Header_Get_Add(spc_ifn_RS_smt,'CNT_HRJ',cont_high_rej,header_comment='Continuum IRAF high sigma rejection')
			except KeyError:
				pass
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'MAG_I'  ,header_comment='i-band magnitude')
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
			Header_Copy(spc_ifn_RS_smt,spc_ifn_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
		else:
			pass
		spc_ifn_RS     = spc_ifn_RS_smt
	elif smt_spc_pre == False:
		pass

	spec_RS      = Spectra_x_y(spc_ifn_RS)
	spc_ifn_RS   = spec_RS[5]

	spec_RS_big_lims = new_wave_range(new_wavelength_array[0],new_wavelength_array[-1],spec_RS[3])
	
	min_el       = float(spec_RS[0][0])
	max_el       = float(spec_RS[0][-1])
	elements     = np.array([min_el,max_el])
	fillin       = np.digitize(elements,spec_RS_big_lims[0])

	cae_abajo    = fillin[0]
	cae_arriba   = fillin[1]

	falta_abajo  = cae_abajo
	falta_arriba = abs(len(new_wavelength_array) - cae_arriba + 1)

	prev_lambda  = np.arange(spec_RS[0][0]  - (spec_RS[3] * (falta_abajo)),spec_RS[0][0]-spec_RS[3],spec_RS[3])
	post_lambda  = np.arange(spec_RS[0][-1] + spec_RS[3],spec_RS[0][-1]+(spec_RS[3]*falta_arriba),spec_RS[3])

	if len(prev_lambda)>0:
		if   np.isclose(prev_lambda[-1],spec_RS[0][0],rtol=1e-1,atol=0.0) == True:
			pass
		elif np.isclose(prev_lambda[-1],spec_RS[0][0],rtol=1e-1,atol=0.0) == False:
			print
			print 'Lower side overlap: No Match!'
			print abs((prev_lambda[-1]) - spec_RS[0][0])
			print (1e-1 * abs(spec_RS[0][0]))
			print abs((prev_lambda[-1]) - spec_RS[0][0]) < (1e-1 * abs(spec_RS[0][0]))
			print prev_lambda[-1],'+', spec_RS[3],'dif',spec_RS[0][0],1e-1
			print spec_RS[3]
			print prev_lambda[0],'-',prev_lambda[-1]
			print spec_RS[0][0],'-',spec_RS[0][-1]
			print 1E-1,1E-3,1E-4
			print 'Quitting 1'
			quit()
	else:
		pass

	if len(post_lambda)>0:
		if   np.isclose(post_lambda[0],spec_RS[0][-1],rtol=1e-1,atol=0.0)==True:
			pass
		elif np.isclose(post_lambda[0],spec_RS[0][-1],rtol=1e-1,atol=0.0)==False:
			print
			print 'Upper side overlap: No Match!'
			print abs((post_lambda[0]) - spec_RS[0][-1])
			print (1e-1 * abs(spec_RS[0][0]))
			print abs((post_lambda[0]) - spec_RS[0][-1]) < (1e-1 * abs(spec_RS[0][-1]))
			print post_lambda[0],'-',spec_RS[3],'dif',spec_RS[0][-1],1e-1
			print spec_RS[3]
			print spec_RS[0][0],'-',spec_RS[0][-1]
			print post_lambda[0],'-',post_lambda[-1]
			print 'Quitting 2'
			quit()
	else:
		pass
	new_lambda    = np.concatenate((prev_lambda,spec_RS[0],post_lambda),axis=0)

	if monotonic(new_lambda)==True:
		pass
	elif monotonic(new_lambda)==False:
		print
		print 'Generated Spectra: Not monotonic!'
		print 'Quitting 3'
		quit()

	if len(new_lambda) - 1 == len(new_wavelength_array):
		new_lambda = np.delete(new_lambda,-1,0)
	elif len(new_lambda)  + 1 == len(new_wavelength_array):
		missing = [np.around(new_lambda[-1] + spec_RS[3],2)]
		np.append(new_lambda,missing)
		new_lambda = np.concatenate((new_lambda,missing),0)
	elif len(new_lambda)  == len(new_wavelength_array):
		pass

	if len(new_lambda) == len(new_wavelength_array):
		pass
	elif len(new_lambda) != len(new_wavelength_array):
		len_dif = len(new_lambda) - len(new_wavelength_array)
		for deleting in range(len_dif):
			new_lambda    = np.delete(new_lambda,-1,0)

	prev_intens = np.empty(falta_abajo)
	post_intens = np.empty(falta_arriba)

	prev_intens[:] = np.nan
	post_intens[:] = np.nan

	new_intensity  = np.concatenate((prev_intens,spec_RS[1],post_intens),axis=0)

		
	if len(new_lambda) == len(new_intensity):
		pass
	elif len(new_lambda) != len(new_intensity):
		len_dif = (len(new_intensity) - len(new_lambda))
		for deleting in range(len_dif):
			new_intensity = np.delete(new_intensity,-1,0)
	spec_res_int = rebin_spec(new_lambda, new_intensity, new_wavelength_array)
	if '-BS-' in spc_tbl_RS:
		if SS_indx==None:
			Int_Spc_ofn  = str(spc_ifn_RS.split('.fits',1)[0]) + '-BS-'+ (spc_tbl_RS.split('-BS-')[-1]).split(tbl_ext_opt)[0] +'-int.fits'
		elif SS_indx!=None:
			Int_Spc_ofn  = str(spc_ifn_RS.split('.fits',1)[0]) + '-BS-'+ (spc_tbl_RS.split('-BS-')[-1]).split(tbl_ext_opt)[0] +'-int-'+str(SS_indx+1)+'.fits'
	elif '-BS_MST_' in spc_tbl_RS:
		if SS_indx==None:
			Int_Spc_ofn  = ind_bst_lst + str((str(spc_ifn_RS.split('.fits',1)[0]) + '-int.fits').rsplit('/',1)[1])
		elif SS_indx!=None:
			Int_Spc_ofn  = ind_bst_lst + str((str(spc_ifn_RS.split('.fits',1)[0]) + '-int-'+str(SS_indx+1)+'.fits').rsplit('/',1)[1])
	else:
		if SS_indx==None:
			Int_Spc_ofn  = str(spc_ifn_RS.split('.fits',1)[0]) + '-int.fits'
		elif SS_indx!=None:
			Int_Spc_ofn  = str(spc_ifn_RS.split('.fits',1)[0]) + '-int-'+str(SS_indx+1)+'.fits'


	if os.path.exists(Int_Spc_ofn):
		os.system('rm ' + str(Int_Spc_ofn))
	else:
		pass
	hdu          = fits.PrimaryHDU(spec_res_int)
	hdulist      = fits.HDUList([hdu])

	hdulist.writeto(Int_Spc_ofn,overwrite=True)

	Header_Updt(Int_Spc_ofn,'CRVAL1',spec_RS_big_lims[4])
	Header_Updt(Int_Spc_ofn,'CDELT1',new_lambda_bin)
	Header_Updt(Int_Spc_ofn,'CD1_1' ,new_lambda_bin)
	Header_Updt(Int_Spc_ofn,'CD2_2' ,new_lambda_bin)

	if '-BS_MST' in spc_tbl_RS:
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'Z_0'  ,header_comment='Redshift')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'Z_REF',header_comment='Redshift Reference')
	else:
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'ID_0'  ,header_comment='ID')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'ID_REF',header_comment='ID Reference')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'Z_0'   ,header_comment='Redshift')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'Z_REF' ,header_comment='Redshift Reference')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'MAG_I' ,header_comment='i-band magnitude')

	if 	sel_pre_cnt == True and 'noise' not in spc_ifn_RS and '-BS_MST' not in spc_tbl_RS:
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'MAG_I'  ,header_comment='i-band magnitude')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'SHT_FNS',header_comment='Spec file post-shift')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'SHT_FN0',header_comment='Spec file post-shift')
	elif sel_pre_cnt == True and 'noise' not in spc_ifn_RS and '-BS_MST' in spc_tbl_RS:
		Header_Get_Add(Int_Spc_ofn,'CNT_TYP',pre_cnt_typ,header_comment='Continuum IRAF type')
		Header_Get_Add(Int_Spc_ofn,'CNT_FIT',pre_cnt_fnc,header_comment='Continuum IRAF function')
		Header_Get_Add(Int_Spc_ofn,'CNT_ORD',pre_cnt_ord,header_comment='Continuum IRAF order')
		Header_Get_Add(Int_Spc_ofn,'CNT_RPC',pre_cnt_rpl,header_comment='Continuum IRAF replacement')
		Header_Get_Add(Int_Spc_ofn,'CNT_LRJ',pre_cnt_lrj,header_comment='Continuum IRAF low sigma rejection')
		Header_Get_Add(Int_Spc_ofn,'CNT_HRJ',pre_cnt_hrj,header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'MAG_I'  ,header_comment='i-band magnitude')
		try:
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
		except KeyError:
			print
			print '6-a'
			print
			print
			print colored('Header not Found: !','yellow')
			print colored('File: ' + spc_ifn_RS,'yellow')
			print colored('CNT_FNM','yellow')
			print colored('CNT_FN0','yellow')
			print colored('CNT_FNF','yellow')
			print
			print colored ('They will NOT be copied','yellow')
			print colored('Into file: '+Int_Spc_ofn,'yellow')
			print colored('From file: '+spc_ifn_RS,'yellow')
			print
			prev_smt_cnt_file = spc_ifn_RS.split('-smt')[0]+'.fits'
			if 'smt' in spc_ifn_RS and os.path.exists(prev_smt_cnt_file):
				print
				print '6-a-1'
				print
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
			elif not os.path.exists(prev_smt_cnt_file):
				print
				print '6-a-2'
				print				
				print
				print colored('File: '+ os.path.exists(prev_smt_cnt_file) + 'does not exists','yellow')
				print colored('From which headers would be copied','yellow')
				print colored('Quitting! (lne-2726 in Fnc_Stk_Stt.py','yellow')
				print
				quit()
			pass
			pass
	elif sel_pre_cnt == False:
		Header_Get_Add(Int_Spc_ofn,'CNT_TYP',pre_cnt     ,header_comment='Continuum IRAF PRE')
		pass

	if smt_spc_pre == True and 'noise' not in spc_ifn_RS and '-BS_MST' not in spc_tbl_RS:
		Header_Get_Add(Int_Spc_ofn,'PRE_SMT',str(smt_spc_pre) ,header_comment='Pre-Smooth')
		Header_Get_Add(Int_Spc_ofn,'PSM_KRN',str(smt_shp_pre) ,header_comment='Pre-Smooth kernel type')
		Header_Get_Add(Int_Spc_ofn,'PSM_SZE',int(smt_sze_pre) ,header_comment='Pre-Smooth kernel size')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'MAG_I'  ,header_comment='i-band magnitude')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'SHT_FNS',header_comment='Spec file post-shift')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'SHT_FN0',header_comment='Spec file post-shift')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'MAG_I'  ,header_comment='i-band magnitude')
		try:
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
		except KeyError: 
			print
			print '7-b'
			print
			print
			print colored('Header not Found: !','yellow')
			print colored('File: ' + Int_Spc_ofn,'yellow')
			print colored('CNT_FNM','yellow')
			print colored('CNT_FN0','yellow')
			print colored('CNT_FNF','yellow')
			print
			print colored ('They will NOT be copied','yellow')
			print colored('Into file: '+Int_Spc_ofn,'yellow')
			print colored('From file: '+spc_ifn_RS,'yellow')
			print
			prev_smt_cnt_file = spc_ifn_RS.split('-smt')[0]+'.fits'
			if 'smt' in spc_ifn_RS and os.path.exists(prev_smt_cnt_file):
				print
				print '7-b-1'
				print
				print
				print colored('Copying headers from pre-smooth spectra!','yellow')
				print colored('From File: ' + prev_smt_cnt_file,'yellow')
				print colored('Into File: ' + Int_Spc_ofn,'yellow')
				print colored('And also','yellow')
				print colored('Into File: ' + spc_ifn_RS,'yellow')
				print
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
			elif not os.path.exists(prev_smt_cnt_file):
				print
				print '7-b-2'
				print
				print
				print colored('File: '+ os.path.exists(prev_smt_cnt_file) + 'does not exists','yellow')
				print colored('From which headers would be copied','yellow')
				print colored('Quitting! (lne-2726 in Fnc_Stk_Stt.py','yellow')
				print
				print '7-$$$$$'
				quit()
	elif smt_spc_pre == True and 'noise' not in spc_ifn_RS and '-BS_MST' in spc_tbl_RS:
		Header_Get_Add(Int_Spc_ofn,'PRE_SMT',str(smt_spc_pre) ,header_comment='Pre-Smooth')
		Header_Get_Add(Int_Spc_ofn,'PSM_KRN',str(smt_shp_pre) ,header_comment='Pre-Smooth kernel type')
		Header_Get_Add(Int_Spc_ofn,'PSM_SZE',int(smt_sze_pre) ,header_comment='Pre-Smooth kernel size')
		Header_Copy(Int_Spc_ofn,spc_ifn_RS,'MAG_I'  ,header_comment='i-band magnitude')
		try:
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'SHT_FNS',header_comment='Spec file post-shift')
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'SHT_FN0',header_comment='Spec file post-shift')
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
		except KeyError:
			print
			print '8-b'
			print
			print
			print colored('Header not Found: !','yellow')
			print colored('File: ' + spc_ifn_RS,'yellow')
			print colored('CNT_FNM','yellow')
			print colored('CNT_FN0','yellow')
			print colored('CNT_FNF','yellow')
			print
			print colored ('They will NOT be copied','yellow')
			print colored('Into file: '+Int_Spc_ofn,'yellow')
			print colored('From file: '+spc_ifn_RS,'yellow')
			print
			prev_smt_cnt_file = spc_ifn_RS.split('-smt')[0]+'.fits'
			if 'smt' in spc_ifn_RS and os.path.exists(prev_smt_cnt_file):
				print
				print '8-b-1'
				print
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
				Header_Copy(Int_Spc_ofn,prev_smt_cnt_file,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
				Header_Copy(spc_ifn_RS,prev_smt_cnt_file,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')
			elif not os.path.exists(prev_smt_cnt_file):
				print
				print '8-b-2'
				print
				print
				print colored('File: '+ os.path.exists(prev_smt_cnt_file) + 'does not exists','yellow')
				print colored('From which headers would be copied','yellow')
				print colored('Quitting! (lne-2726 in Fnc_Stk_Stt.py','yellow')
				print
				print '8-$$$$$'
				quit()

	elif smt_spc_pre == False:
		print
		print '9'
		print
		Header_Get_Add(Int_Spc_ofn,'PRE_SMT',str(smt_spc_pre) ,header_comment='Pre-Smooth')

	Header_Copy(Int_Spc_ofn,spc_ifn_RS,'h_s_0')
	Header_Copy(Int_Spc_ofn,spc_ifn_RS,'h_s_c')
	Header_History_Step(spc_ifn_RS,Int_Spc_ofn)
	loop_back_step = int(Header_Get(Int_Spc_ofn,'h_s_c'))
	while (loop_back_step > 0):
		loop_back_step = int(loop_back_step - 1)
		loop_back_step_header = 'h_s_' + str(loop_back_step)
		try:
			Header_Copy(Int_Spc_ofn,spc_ifn_RS,loop_back_step_header)
		except KeyError:
			pass
	else:
		pass


	if pre_msk == True:
		z_ref_msk    = (1+ Header_Get(Int_Spc_ofn,'Z_0'))/ (1+ Header_Get(Int_Spc_ofn,'Z_REF'))
		spc_opt_msk  = Spectra_Masking(Int_Spc_ofn,msk_typ=pre_msk_typ,rshft_corr=z_ref_msk,rshft_corr_direct=True,
						msk_abs_lne=pre_msk_abs_lne,
						msk_blu_rgn=pre_msk_rgn,blu_lmb_min=pre_msk_min,blu_lmb_max=pre_msk_max)
		Int_Spc_ofn  = spc_opt_msk
		Header_Get_Add(Int_Spc_ofn,'SCP_PMK',str(pre_msk)        ,header_comment='Pre-Masking')
		Header_Get_Add(Int_Spc_ofn,'MSK_TYP',str(pre_msk_typ)    ,header_comment='Pre-Masking Type')
		Header_Get_Add(Int_Spc_ofn,'MSK_VAL',int(pre_msk_cte_val),header_comment='Pre-Masking Constant Values')
		Header_Get_Add(Int_Spc_ofn,'MSK_PXL',str(pre_msk_rgn)    ,header_comment='Pre-Masking Spec Region')
		Header_Get_Add(Int_Spc_ofn,'PIX_INT',int(pre_msk_min)    ,header_comment='Pre-Masking Lower lambda limit')
		Header_Get_Add(Int_Spc_ofn,'PIX_LST',int(pre_msk_max)    ,header_comment='Pre-Masking Upper lambda limit')
		Header_Get_Add(Int_Spc_ofn,'MSK_LNS',str(pre_msk_abs_lne),header_comment='Pre-Masking Absorption Lines')
	elif pre_msk == False:
		Header_Get_Add(Int_Spc_ofn,'SCP_PMK',str(pre_msk)        ,header_comment='Pre-Masking')
		pass

	if get_cont_flux == True:
		if '-BS_MST' not in Int_Spc_ofn:
			cnt_hdr = 'SHT_FNS'
		else:
			pass
		if '-BS_MST' in Int_Spc_ofn or ('-BS-' in Int_Spc_ofn and '-stk-' in Int_Spc_ofn):
			cnt_hdr = 'CNT_FNF'
		else:
			pass
		cont_val = Spectra_Cont_GetVal(str(Header_Get(spc_ifn_RS,cnt_hdr)),gcv_lmbd_i=gcv_lmbd_i,gcv_lmbd_f=gcv_lmbd_f)
		CFX_SUM  = bn.nansum(cont_val[0])
		CFX_MED  = bn.nanmedian(cont_val[0])
		CFX_AVG  = bn.nanmean(cont_val[0])
		CFX_NUM  = cont_val[1]

		if np.isnan(CFX_SUM) == True or CFX_SUM == 'NAN.0':
			CFX_SUM = 0.0
		else:
			pass
		if np.isnan(CFX_MED) == True or CFX_MED == 'NAN.0':
			CFX_MED = 0.0
		else:
			pass
		if np.isnan(CFX_AVG) == True or CFX_AVG == 'NAN.0':
			CFX_AVG = 0.0
		else:
			pass
		if np.isnan(CFX_NUM) == True or CFX_NUM == 'NAN.0':
			CFX_NUM = 0.0
		else:
			pass

		Header_Get_Add(Int_Spc_ofn,'CFX_TOP',gcv_lmbd_i,header_comment='Weight lambda lower limit')
		Header_Get_Add(Int_Spc_ofn,'CFX_BOT',gcv_lmbd_f,header_comment='Weight lambda upper limit')
		Header_Get_Add(Int_Spc_ofn,'CFX_NUM',CFX_NUM   ,header_comment='Weight number of pixels used')
		Header_Get_Add(Int_Spc_ofn,'CFX_SUM',CFX_SUM   ,header_comment='Weight sum')
		Header_Get_Add(Int_Spc_ofn,'CFX_MED',CFX_MED   ,header_comment='Weight median')
		Header_Get_Add(Int_Spc_ofn,'CFX_AVG',CFX_AVG   ,header_comment='Weight average')
		Header_Get_Add(Int_Spc_ofn,'WGT_UNT',1         ,header_comment='Weight unitary')
		Int_Spc_ofn
	elif get_cont_flux == False:
		Header_Get_Add(Int_Spc_ofn,'WGT_UNT',1         ,header_comment='Weight unitary')
	return Int_Spc_ofn
	
def Stack_Img_Op(stck_img_op_stck,name,img_2bstack,spc_nse,wght_img_2bstack,*args, **kwargs):
	wgt_var          = kwargs.get('wgt_var'        ,None)
	wgt_nrm_hdr      = kwargs.get('wgt_nrm_hdr'    ,None)
	wrt_fits         = kwargs.get('wrt_fits'       ,True)
	stack_ext        = kwargs.get('stack_ext'      ,None)
	new_CRVAL1_head  = kwargs.get('new_CRVAL1_head',None)
	new_CDELT1_head  = kwargs.get('new_CDELT1_head',None)
	bs_func          = kwargs.get('bs_func'        ,'')

	sel_pre_shf     = kwargs.get('sel_pre_shf'     ,True)
	sel_pre_cnt     = kwargs.get('sel_pre_cnt'     ,True)
	sel_pre_msk     = kwargs.get('sel_pre_msk'     ,False)

	pre_cnt         = kwargs.get('pre_cnt'         ,False)
	pre_cnt_typ     = kwargs.get('pre_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pre_cnt_lns     = kwargs.get('pre_cnt_lns'     ,'*')       # Image lines to be fit
	pre_cnt_fnc     = kwargs.get('pre_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pre_cnt_ord     = kwargs.get('pre_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pre_cnt_ovr     = kwargs.get('pre_cnt_ovr'     ,'yes')     # Override previous norm spec
	pre_cnt_rpl     = kwargs.get('pre_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pre_cnt_lrj     = kwargs.get('pre_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pre_cnt_hrj     = kwargs.get('pre_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pre     = kwargs.get('smt_spc_pre'     ,True)
	smt_shp_pre     = kwargs.get('smt_shp_pre'     ,'gaussian')
	smt_sze_pre     = kwargs.get('smt_sze_pre'     ,1)

	pre_msk         = kwargs.get('pre_msk'         ,False)
	pre_msk_typ     = kwargs.get('pre_msk_typ'     ,'NaN')
	pre_msk_cte_val = kwargs.get('pre_msk_cte_val' ,1)
	pre_msk_abs_lne = kwargs.get('pre_msk_abs_lne' ,False)
	pre_msk_rgn     = kwargs.get('pre_msk_rgn'     ,False)
	pre_msk_min     = kwargs.get('pre_msk_min'     ,500)
	pre_msk_max     = kwargs.get('pre_msk_max'     ,1210)

	sig_clp         = kwargs.get('sig_clp'         ,False)
	sig_cut         = kwargs.get('sig_cut'         ,3)
	sig_fct         = kwargs.get('sig_fct'         ,mean)
	sig_fll         = kwargs.get('sig_fll'         ,np.nan)

	wgt_typ         = kwargs.get('wgt_typ'         ,None)
	get_cont_flux   = kwargs.get('get_cont_flux'   ,True)
	gcv_lmbd_i      = kwargs.get('gcv_lmbd_i'      ,1430)
	gcv_lmbd_f      = kwargs.get('gcv_lmbd_f'      ,1480)

	spc_nse         = kwargs.get('spc_nse'         ,False)
	noise_imgs      = kwargs.get('noise_imgs'      ,None)

	pst_cnt         = kwargs.get('pst_cnt'         ,True)
	pst_cnt_typ     = kwargs.get('pst_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pst_cnt_lns     = kwargs.get('pst_cnt_lns'     ,'*')       # Image lines to be fit
	pst_cnt_fnc     = kwargs.get('pst_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pst_cnt_ord     = kwargs.get('pst_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pst_cnt_ovr     = kwargs.get('pst_cnt_ovr'     ,'yes')     # Override previous norm spec
	pst_cnt_rpl     = kwargs.get('pst_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pst_cnt_lrj     = kwargs.get('pst_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pst_cnt_hrj     = kwargs.get('pst_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pst     = kwargs.get('smt_spc_pst'     ,True)
	smt_shp_pst     = kwargs.get('smt_shp_pst'     ,'gaussian')
	smt_sze_pst     = kwargs.get('smt_sze_pst'     ,1)

	test_bg         = kwargs.get('test_bg'         ,False)
	test_fg         = kwargs.get('test_fg'         ,False)

	img_res         = np.zeros(shape=img_2bstack[0].shape)

	bs_nmb_itr      = kwargs.get('bs_nmb_itr'      ,999999)

	stk_pct_mde     = kwargs.get('stk_pct_mde'     ,False)
	stk_wgt_mde     = kwargs.get('stk_wgt_mde'     ,False)

	print
	print colored('Pre-Processed spectra files properties : ','yellow')
	print colored('Spectra shifted                        : '+str(sel_pre_shf),'yellow')
	print colored('Spectra cont fitted                    : '+str(sel_pre_cnt),'yellow')
	print colored('Spectra masked                         : '+str(sel_pre_msk),'yellow')
	print

	print 'Number of galaxies to be stacked (histogram): ',len(img_2bstack)
	print
	print colored('Pre-Stacking continuum       : '+ str(pre_cnt)    ,'yellow')
	print colored('Pre-Stacking smooth          : '+ str(smt_spc_pre),'yellow')
	print colored('Pre-Stacking mask            : '+ str(pre_msk)    ,'yellow')
	print colored('Sigma-clipping for stacking  : '+ str(sig_clp)    ,'yellow')
	print

	if sig_clp == True:
		img_flt       = apsts.sigma_clip(img_2bstack,sigma=sig_cut,axis=0,iters=None,cenfunc=sig_fct, copy=True)

		print
		print colored('Sigma-clipping for stacking!','yellow')
		print colored('Sigma Cut                    : ' + str(sig_cut),'yellow')
		print colored('Central function             : ' + str(sig_fct), 'yellow')
		print colored('Rejected values replaced with: ' + str(sig_fll),'yellow')
		print

		img_flt.set_fill_value(sig_fll)
		img_flt_filled = img_flt.filled()
		img_stat       = img_flt_filled
	elif sig_clp == False:
		img_stat   = img_2bstack

	wght_img_copy = wght_img_2bstack
	wght_img_stat = np.reshape(wght_img_2bstack,(len(img_stat),1))
	wght_img_stat = np.asarray(wght_img_stat)
	img_staw      = img_stat * wght_img_stat
	
	Transpose  = np.asarray(img_stat).T
	Transposw  = np.asarray(img_staw).T

	img_stat_hst = []
	img_stat_hsw = []
	img_stat_smw = []
	wght_img_stat = np.squeeze(wght_img_stat)

	for j in range(len(Transpose)):
		wght_img_stat_j =  Transposw[j]/Transposw[j]* wght_img_stat
		if np.isnan(sig_fll) == True:
			non_msk_num = int(np.count_nonzero(~np.isnan(Transpose[j])))
			msk_num     = int(np.count_nonzero(np.isnan(Transpose[j])))
			img_stat_hst.append(float(non_msk_num))

			non_msk_num_wghts = int(np.count_nonzero(~np.isnan(wght_img_stat_j)))
			msk_num_wghts     = int(np.count_nonzero(np.isnan(wght_img_stat_j)))
			img_stat_hsw.append(float(non_msk_num_wghts))
			img_stat_smw.append(float(bn.nansum(wght_img_stat_j)))

		elif np.isnan(sig_fll) == False:
			pass
			non_msk_num = int(np.count_nonzero(Transpose[j]!=sig_fll))
			img_stat_hst.append(float(non_msk_num))

			non_msk_num_wghts = int(np.count_nonzero(wght_img_stat_j!=sig_fll))
			img_stat_hsw.append(float(non_msk_num_wghts))
			img_stat_smw.append(float(bn.nansum(wght_img_stat_j)))

		else:
			pass


	img_sts_hst = np.asarray(img_stat_hst)
	img_sts_frc = img_sts_hst/len(img_2bstack)

	img_res_sum = bn.nansum(np.array(img_stat)     , axis=0)
	img_res_avg = bn.nanmean(np.array(img_stat)    , axis=0)
	img_res_med = bn.nanmedian(np.array(img_stat)  , axis=0)
	img_res_std = bn.nanstd(np.array(img_stat)     , axis=0)
	img_res_rms = np.sqrt((img_res_avg**2) + (img_res_std**2))
	
	img_sts_hsw = np.asarray(img_stat_hsw)
	img_sts_wsu = np.asarray(img_stat_smw)
	img_res_suw = bn.nansum(np.array(img_staw), axis=0)
	img_res_avw = img_res_suw/img_sts_wsu

	print colored('Stacked images through : sum, mean, median, + weighted stats: ','yellow')
	print
	print colored('With dimensions: ','yellow')
	print colored(str(img_res_sum.shape)+str(img_res_sum.dtype),'yellow')
	print		

	if stk_pct_mde == True:
		print colored('Stacked images through : percentiles: ','yellow')
		print colored('17., 83.0, (1 sigma)','yellow')
		print colored('2.5, 97.5, (2 sigma)','yellow')
		print colored('0.5, 99.5, (3 sigma)','yellow')
		print colored('25., 75.0, (interquantile)','yellow')
		print
		print colored('With dimensions: ','yellow')
		print colored(str(img_res_sum.shape)+str(img_res_sum.dtype),'yellow')
		print		
		img_res_1sl = np.nanpercentile(np.array(img_stat), 15.9, axis=0)
		img_res_1sh = np.nanpercentile(np.array(img_stat), 84.1, axis=0)

		img_res_2sl = np.nanpercentile(np.array(img_stat), 2.30, axis=0)
		img_res_2sh = np.nanpercentile(np.array(img_stat), 97.7, axis=0)

		img_res_3sl = np.nanpercentile(np.array(img_stat), 0.20, axis=0)
		img_res_3sh = np.nanpercentile(np.array(img_stat), 99.8, axis=0)

		img_res_p25 = np.nanpercentile(np.array(img_stat), 25, axis=0)
		img_res_p75 = np.nanpercentile(np.array(img_stat), 75, axis=0)
		
	else:
		pass
	if stk_wgt_mde == True:
		print colored('Stacked images through : weighted-auxiliary files ','yellow')
		print
		print colored('With dimensions: ','yellow')
		print colored(str(img_res_sum.shape)+str(img_res_sum.dtype),'yellow')
		print		
		wgt_sts_val = np.asarray(wght_img_copy)
		wgt_sts_vln = np.asarray(wght_img_copy) / bn.nansum(np.asarray(wght_img_copy))
		wgt_sts_hst = np.asarray(wght_img_copy) / bn.nansum(np.asarray(wght_img_copy))
	else:
		pass
	if wrt_fits==True:
		if  '-BS-' in name:
			print colored(name,'yellow')
			spc_dir_dst = str_bst_stk
		elif  '-BS_MST' in name:
			print colored(name,'yellow')
			spc_dir_dst = stt_bst_stk
		else:
			spc_dir_dst = res_stk_res

		if test_bg == True or test_fg == True:
			test_stk_sfx = '-rf'
			print
			print colored('Warning: Background galaxies check mode selected!','green')
			print colored('Adding: ' + test_stk_sfx + ' to the stacked files','green')
			print
		else:
			test_stk_sfx = ''

		spec_file_sum_ofn = spc_dir_dst + str(name) + bs_func + '-stk-sum' + test_stk_sfx + '.fits'
		spec_file_avg_ofn = spc_dir_dst + str(name) + bs_func + '-stk-avg' + test_stk_sfx + '.fits'
		spec_file_med_ofn = spc_dir_dst + str(name) + bs_func + '-stk-med' + test_stk_sfx + '.fits'
		spec_file_hst_ofn = spc_dir_dst + str(name) + bs_func + '-stk-hst' + test_stk_sfx + '.fits'
		spec_file_std_ofn = spc_dir_dst + str(name) + bs_func + '-stk-std' + test_stk_sfx + '.fits'
		spec_file_rms_ofn = spc_dir_dst + str(name) + bs_func + '-stk-rms' + test_stk_sfx + '.fits'

		spec_file_sum     = Wrt_FITS_File(img_res_sum,spec_file_sum_ofn)
		spec_file_avg     = Wrt_FITS_File(img_res_avg,spec_file_avg_ofn)
		spec_file_med     = Wrt_FITS_File(img_res_med,spec_file_med_ofn)
		spec_file_hst     = Wrt_FITS_File(img_sts_hst,spec_file_hst_ofn)
		spec_file_std     = Wrt_FITS_File(img_res_std,spec_file_std_ofn)
		spec_file_rms     = Wrt_FITS_File(img_res_rms,spec_file_rms_ofn)

		spec_file_avw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-avw' + test_stk_sfx + '.fits'
		spec_file_suw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-suw' + test_stk_sfx + '.fits'
		spec_file_wsu_ofn = spc_dir_dst + str(name) + bs_func + '-stk-wsu' + test_stk_sfx + '.fits'
		spec_file_hsw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-hsw' + test_stk_sfx + '.fits'

		spec_file_avw     = Wrt_FITS_File(img_res_avw,spec_file_avw_ofn)
		spec_file_suw     = Wrt_FITS_File(img_res_suw,spec_file_suw_ofn)
		spec_file_wsu     = Wrt_FITS_File(img_sts_wsu,spec_file_wsu_ofn)
		spec_file_hsw     = Wrt_FITS_File(img_sts_hsw,spec_file_hsw_ofn)

		OPT_STCK_CRE      = [
							spec_file_sum,spec_file_avg,spec_file_med,
							spec_file_hst,
							spec_file_std,spec_file_rms,
							spec_file_hsw,
							spec_file_wsu,spec_file_suw,spec_file_avw
							]

		[Header_Updt(spec_sts_res,'CRVAL1',new_CRVAL1_head) for spec_sts_res in OPT_STCK_CRE]
		[Header_Updt(spec_sts_res,'CDELT1',new_CDELT1_head) for spec_sts_res in OPT_STCK_CRE]
		[Header_Updt(spec_sts_res,'CD1_1' ,new_CDELT1_head) for spec_sts_res in OPT_STCK_CRE]
		[Header_Updt(spec_sts_res,'CD2_2' ,new_CDELT1_head) for spec_sts_res in OPT_STCK_CRE]

		[Header_Get_Add(spec_sts_res,'STK_NUM',str(len(img_2bstack)),header_comment='Number of galaxies used for Stack') for spec_sts_res in OPT_STCK_CRE]
		[Header_Get_Add(spec_sts_res,'h_s_c',0,header_comment='History Step Last')             for spec_sts_res in OPT_STCK_CRE]
		[Header_Get_Add(spec_sts_res,'h_s_0',str((spec_sts_res.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='History Step Init') for spec_sts_res in OPT_STCK_CRE]

		if stk_wgt_mde == True:		
			wght_file_val_ofn = spc_dir_dst + str(name) + bs_func + '-wgt-val' + test_stk_sfx + '.fits'
			wght_file_vln_ofn = spc_dir_dst + str(name) + bs_func + '-wgt-vln' + test_stk_sfx + '.fits'
			wght_file_hst_ofn = spc_dir_dst + str(name) + bs_func + '-wgt-hst' + test_stk_sfx + '.fits'

			wght_file_val     = Wrt_FITS_File(wgt_sts_val,wght_file_val_ofn)
			wght_file_vln     = Wrt_FITS_File(wgt_sts_vln,wght_file_vln_ofn)
			wght_file_hst     = Wrt_FITS_File(wgt_sts_hst,wght_file_hst_ofn)

			OPT_WGHT_FLS      = [wght_file_val,wght_file_vln,wght_file_hst]

			[Header_Updt(spec_sts_res,'CRVAL1',1) for spec_sts_res in OPT_WGHT_FLS]
			[Header_Updt(spec_sts_res,'CDELT1',1) for spec_sts_res in OPT_WGHT_FLS]
			[Header_Updt(spec_sts_res,'CD1_1' ,1) for spec_sts_res in OPT_WGHT_FLS]
			[Header_Updt(spec_sts_res,'CD2_2' ,1) for spec_sts_res in OPT_WGHT_FLS]

			OPT_WGT_OPR        = [
								'WVL',
								'WVN',
								'WHS',
								]

		else:
			pass

		if stk_pct_mde == True:
			spec_file_p25_ofn = spc_dir_dst + str(name) + bs_func + '-stk-p25' + test_stk_sfx + '.fits'
			spec_file_p75_ofn = spc_dir_dst + str(name) + bs_func + '-stk-p75' + test_stk_sfx + '.fits'
			spec_file_1sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-1sl' + test_stk_sfx + '.fits'
			spec_file_1sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-1sh' + test_stk_sfx + '.fits'
			spec_file_2sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-2sl' + test_stk_sfx + '.fits'
			spec_file_2sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-2sh' + test_stk_sfx + '.fits'
			spec_file_3sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-3sl' + test_stk_sfx + '.fits'
			spec_file_3sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-3sh' + test_stk_sfx + '.fits'

			spec_file_p25     = Wrt_FITS_File(img_res_p25,spec_file_p25_ofn)
			spec_file_p75     = Wrt_FITS_File(img_res_p75,spec_file_p75_ofn)
			spec_file_1sl     = Wrt_FITS_File(img_res_1sl,spec_file_1sl_ofn)
			spec_file_1sh     = Wrt_FITS_File(img_res_1sh,spec_file_1sh_ofn)
			spec_file_2sl     = Wrt_FITS_File(img_res_2sl,spec_file_2sl_ofn)
			spec_file_2sh     = Wrt_FITS_File(img_res_2sh,spec_file_2sh_ofn)
			spec_file_3sl     = Wrt_FITS_File(img_res_3sl,spec_file_3sl_ofn)
			spec_file_3sh     = Wrt_FITS_File(img_res_3sh,spec_file_3sh_ofn)

			OPT_STCK_PCT      = [
								spec_file_p25,spec_file_p75,
								spec_file_1sl,spec_file_1sh,
								spec_file_2sl,spec_file_2sh,
								spec_file_3sl,spec_file_3sh
								]

			OPT_PCT_OPR        = [
								'P25','P75',
								'1SL','1SH',
								'2SL','2SH',
								'3SL','3SH'
								]
			[Header_Get_Add(spec_sts_res,'STK_NUM',str(len(img_2bstack)),header_comment='Number of galaxies used for Stack') for spec_sts_res in OPT_STCK_PCT]
			[Header_Get_Add(spec_sts_res,'h_s_c',0,header_comment='History Step Last')             for spec_sts_res in OPT_STCK_PCT]
			[Header_Get_Add(spec_sts_res,'h_s_0',str((spec_sts_res.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='History Step Init') for spec_sts_res in OPT_STCK_PCT]

			[Header_Updt(spec_sts_res,'CRVAL1',new_CRVAL1_head) for spec_sts_res in OPT_STCK_PCT]
			[Header_Updt(spec_sts_res,'CDELT1',new_CDELT1_head) for spec_sts_res in OPT_STCK_PCT]
			[Header_Updt(spec_sts_res,'CD1_1' ,new_CDELT1_head) for spec_sts_res in OPT_STCK_PCT]
			[Header_Updt(spec_sts_res,'CD2_2' ,new_CDELT1_head) for spec_sts_res in OPT_STCK_PCT]

		else:
			pass
	else:
		pass

	#######################POST-PROCESSING########################
	#########################CORE-FILES###########################
	if pst_cnt==False and smt_spc_pst == False:
		print
		print colored('Generating CORE composite files','yellow')
		print colored('2-CRE','yellow')
		print
		OPT_STCK_CRE      = [
							spec_file_sum,spec_file_avg,spec_file_med,
							spec_file_hst,
							spec_file_std,spec_file_rms,
							spec_file_hsw,
							spec_file_wsu,spec_file_suw,spec_file_avw
							]

		OPT_STCK_OPR      = [
							'SUM','AVG','MED',
							'HST',
							'STD','RMS',
							'HSW',
							'WSU','SUW','AVW',
							]

		OPT_STCK_WGT      = [
							spec_file_hsw,
							spec_file_wsu,spec_file_suw,spec_file_avw,
							]							
		FNL_SPEC_RES      = OPT_STCK_CRE
	elif pst_cnt==True and smt_spc_pst == False:
		print
		print colored('Generating CORE composite files','yellow')
		print colored('3-CRE','yellow')
		print
		spec_file_sum_cnt  = Spectra_Cont_IRAF(spec_file_sum_ofn,ind_stk_res + 'log_cont-stk-sum',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_avg_cnt  = Spectra_Cont_IRAF(spec_file_avg_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_med_cnt  = Spectra_Cont_IRAF(spec_file_med_ofn,ind_stk_res + 'log_cont-stk-med',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_std_cnt  = Spectra_Cont_IRAF(spec_file_std_ofn,ind_stk_res + 'log_cont-stk-std',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_rms_cnt  = Spectra_Cont_IRAF(spec_file_rms_ofn,ind_stk_res + 'log_cont-stk-rms',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_avw_cnt  = Spectra_Cont_IRAF(spec_file_avw_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_suw_cnt  = Spectra_Cont_IRAF(spec_file_suw_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)

		spec_file_sum_cnt_RS = spec_file_sum_cnt[0]
		spec_file_avg_cnt_RS = spec_file_avg_cnt[0]
		spec_file_med_cnt_RS = spec_file_med_cnt[0]
		spec_file_std_cnt_RS = spec_file_std_cnt[0]
		spec_file_rms_cnt_RS = spec_file_rms_cnt[0]
		spec_file_avw_cnt_RS = spec_file_avw_cnt[0]
		spec_file_suw_cnt_RS = spec_file_suw_cnt[0]

		Header_Copy(spec_file_sum_cnt_RS,spec_file_sum_ofn,'h_s_0')
		Header_Copy(spec_file_sum_cnt_RS,spec_file_sum_ofn,'h_s_c')

		Header_Copy(spec_file_avg_cnt_RS,spec_file_avg_ofn,'h_s_0')
		Header_Copy(spec_file_avg_cnt_RS,spec_file_avg_ofn,'h_s_c')

		Header_Copy(spec_file_med_cnt_RS,spec_file_med_ofn,'h_s_0')
		Header_Copy(spec_file_med_cnt_RS,spec_file_med_ofn,'h_s_c')

		Header_Copy(spec_file_std_cnt_RS,spec_file_std_ofn,'h_s_0')
		Header_Copy(spec_file_std_cnt_RS,spec_file_std_ofn,'h_s_c')

		Header_Copy(spec_file_rms_cnt_RS,spec_file_rms_ofn,'h_s_0')
		Header_Copy(spec_file_rms_cnt_RS,spec_file_rms_ofn,'h_s_c')

		Header_Copy(spec_file_avw_cnt_RS,spec_file_avw_ofn,'h_s_0')
		Header_Copy(spec_file_avw_cnt_RS,spec_file_avw_ofn,'h_s_c')

		Header_Copy(spec_file_suw_cnt_RS,spec_file_suw_ofn,'h_s_0')
		Header_Copy(spec_file_suw_cnt_RS,spec_file_suw_ofn,'h_s_c')
		
		Header_History_Step(spec_file_sum_ofn,spec_file_sum_cnt_RS)
		Header_History_Step(spec_file_avg_ofn,spec_file_avg_cnt_RS)
		Header_History_Step(spec_file_med_ofn,spec_file_med_cnt_RS)
		Header_History_Step(spec_file_std_ofn,spec_file_std_cnt_RS)
		Header_History_Step(spec_file_rms_ofn,spec_file_rms_cnt_RS)
		Header_History_Step(spec_file_avw_ofn,spec_file_avw_cnt_RS)
		Header_History_Step(spec_file_suw_ofn,spec_file_suw_cnt_RS)

		Header_History_Step_Backwards_Loop(spec_file_sum_cnt_RS,spec_file_sum_ofn)
		Header_History_Step_Backwards_Loop(spec_file_avg_cnt_RS,spec_file_avg_ofn)
		Header_History_Step_Backwards_Loop(spec_file_med_cnt_RS,spec_file_med_ofn)
		Header_History_Step_Backwards_Loop(spec_file_std_cnt_RS,spec_file_std_ofn)
		Header_History_Step_Backwards_Loop(spec_file_rms_cnt_RS,spec_file_rms_ofn)
		Header_History_Step_Backwards_Loop(spec_file_avw_cnt_RS,spec_file_avw_ofn)
		Header_History_Step_Backwards_Loop(spec_file_suw_cnt_RS,spec_file_suw_ofn)

		OPT_STCK_CRE      = [
							spec_file_sum,spec_file_avg,spec_file_med,
							spec_file_hst,
							spec_file_std,spec_file_rms,
							spec_file_hsw,
							spec_file_wsu,spec_file_suw,spec_file_avw,								
							spec_file_sum_cnt[0],spec_file_sum_cnt[1],
							spec_file_avg_cnt[0],spec_file_avg_cnt[1],
							spec_file_med_cnt[0],spec_file_med_cnt[1],
							spec_file_std_cnt[0],spec_file_std_cnt[1],
							spec_file_rms_cnt[0],spec_file_rms_cnt[1],
							spec_file_avw_cnt[0],spec_file_avw_cnt[1],
							spec_file_suw_cnt[0],spec_file_suw_cnt[1],
							]		
		OPT_STCK_CNT      = [
							spec_file_sum_cnt[0],
							spec_file_avg_cnt[0],
							spec_file_med_cnt[0],
							spec_file_std_cnt[0],
							spec_file_rms_cnt[0],
							spec_file_avw_cnt[0],
							spec_file_suw_cnt[0]
							]
		OPT_STCK_OPR      = [
							'SUM','AVG','MED',
							'HST',
							'STD','RMS',
							'HSW',
							'WSU','SUW','AVW',
							'SUC','SUF',
							'AVC','AVF',
							'MEC','MEF',
							'STC','STF',
							'RMC','RMF',
							'AWC','AWF',
							'SWC','SWF'
							]
		OPT_STCK_WGT      = [
							spec_file_hsw,
							spec_file_wsu,spec_file_suw,spec_file_avw,
							spec_file_avw_cnt[0],spec_file_avw_cnt[1],
							spec_file_suw_cnt[0],spec_file_suw_cnt[1],
							]	

		FNL_SPEC_RES     = OPT_STCK_CRE
	elif pst_cnt==False and smt_spc_pst == True:
		print
		print colored('Generating CORE composite files','yellow')
		print colored('4-CRE','yellow')
		print
		spec_file_sum_smt     = Spectra_Smooth(spec_file_sum_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_avg_smt     = Spectra_Smooth(spec_file_avg_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_med_smt     = Spectra_Smooth(spec_file_med_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_avw_smt     = Spectra_Smooth(spec_file_avw_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_suw_smt     = Spectra_Smooth(spec_file_suw_ofn,smt_shp_pst,smt_sze_pst)

		try:
			print spec_file_sum_smt,Header_Get(spec_file_sum_smt,'STK_NUM')
			print spec_file_sum,Header_Get(spec_file_sum,'STK_NUM')
			print spec_file_avg_smt,Header_Get(spec_file_avg_smt,'STK_NUM')
			print spec_file_avg,Header_Get(spec_file_avg,'STK_NUM')
			print spec_file_med_smt,Header_Get(spec_file_med_smt,'STK_NUM')
			print spec_file_med,Header_Get(spec_file_med,'STK_NUM')
			print spec_file_avw_smt,Header_Get(spec_file_avw_smt,'STK_NUM')
			print spec_file_avw,Header_Get(spec_file_avw,'STK_NUM')
			print spec_file_suw_smt,Header_Get(spec_file_suw_smt,'STK_NUM')
			print spec_file_suw,Header_Get(spec_file_suw,'STK_NUM')
		except KeyError:
			pass

		Header_Copy(spec_file_sum_smt,spec_file_sum,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_avg_smt,spec_file_avg,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_med_smt,spec_file_med,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_avw_smt,spec_file_avw,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_suw_smt,spec_file_avw,'STK_NUM',header_comment='Number of galaxies used for Stack')

		OPT_STCK_CRE      = [
							spec_file_sum,spec_file_avg,spec_file_med,
							spec_file_hst,
							spec_file_std,spec_file_rms,
							spec_file_hsw,
							spec_file_wsu,spec_file_suw,spec_file_avw,
							spec_file_sum_smt,spec_file_avg_smt,spec_file_med_smt,
							spec_file_avw_smt,spec_file_suw_smt		
							]
		OPT_STCK_SMT      = [
							spec_file_sum_smt,spec_file_avg_smt,spec_file_med_smt,
							spec_file_avw_smt,spec_file_suw_smt
							]
		OPT_STCK_OPR      = [
							'SUM','AVG','MED',
							'HST',
							'STD','RMS',
							'HSW',
							'WSU','SUW','AVW',
							'SUMS','AVGS','MEDS',
							'AVWS','SUWS',
							]

		OPT_STCK_WGT      = [
							spec_file_hsw,
							spec_file_wsu,spec_file_suw,spec_file_avw,
							spec_file_avw_smt,spec_file_suw_smt,
							]								
		FNL_SPEC_RES      = OPT_STCK_CRE
	elif pst_cnt==True and smt_spc_pst == True:
		print
		print colored('Generating CORE composite files','yellow')
		print colored('5-CRE','yellow')
		print

		spec_file_sum_smt     = Spectra_Smooth(spec_file_sum_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_avg_smt     = Spectra_Smooth(spec_file_avg_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_med_smt     = Spectra_Smooth(spec_file_med_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_avw_smt     = Spectra_Smooth(spec_file_avw_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_suw_smt     = Spectra_Smooth(spec_file_suw_ofn,smt_shp_pst,smt_sze_pst)

		Header_Copy(spec_file_sum_smt,spec_file_sum,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_avg_smt,spec_file_avg,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_med_smt,spec_file_med,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_avw_smt,spec_file_avw,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_suw_smt,spec_file_avw,'STK_NUM',header_comment='Number of galaxies used for Stack')

		spec_file_sum_cnt  = Spectra_Cont_IRAF(spec_file_sum_ofn,ind_stk_res + 'log_cont-stk-sum',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_avg_cnt  = Spectra_Cont_IRAF(spec_file_avg_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_med_cnt  = Spectra_Cont_IRAF(spec_file_med_ofn,ind_stk_res + 'log_cont-stk-med',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_std_cnt  = Spectra_Cont_IRAF(spec_file_std_ofn,ind_stk_res + 'log_cont-stk-std',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_rms_cnt  = Spectra_Cont_IRAF(spec_file_rms_ofn,ind_stk_res + 'log_cont-stk-rms',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_avw_cnt  = Spectra_Cont_IRAF(spec_file_avw_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_suw_cnt  = Spectra_Cont_IRAF(spec_file_suw_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)

		spec_file_sum_cnt_RS = spec_file_sum_cnt[0]
		spec_file_avg_cnt_RS = spec_file_avg_cnt[0]
		spec_file_med_cnt_RS = spec_file_med_cnt[0]
		spec_file_std_cnt_RS = spec_file_std_cnt[0]
		spec_file_rms_cnt_RS = spec_file_rms_cnt[0]
		spec_file_avw_cnt_RS = spec_file_avw_cnt[0]
		spec_file_suw_cnt_RS = spec_file_suw_cnt[0]

		Header_Copy(spec_file_sum_cnt_RS,spec_file_sum_ofn,'h_s_0')
		Header_Copy(spec_file_sum_cnt_RS,spec_file_sum_ofn,'h_s_c')

		Header_Copy(spec_file_avg_cnt_RS,spec_file_avg_ofn,'h_s_0')
		Header_Copy(spec_file_avg_cnt_RS,spec_file_avg_ofn,'h_s_c')

		Header_Copy(spec_file_med_cnt_RS,spec_file_med_ofn,'h_s_0')
		Header_Copy(spec_file_med_cnt_RS,spec_file_med_ofn,'h_s_c')

		Header_Copy(spec_file_std_cnt_RS,spec_file_std_ofn,'h_s_0')
		Header_Copy(spec_file_std_cnt_RS,spec_file_std_ofn,'h_s_c')

		Header_Copy(spec_file_rms_cnt_RS,spec_file_rms_ofn,'h_s_0')
		Header_Copy(spec_file_rms_cnt_RS,spec_file_rms_ofn,'h_s_c')

		Header_Copy(spec_file_avw_cnt_RS,spec_file_avw_ofn,'h_s_0')
		Header_Copy(spec_file_avw_cnt_RS,spec_file_avw_ofn,'h_s_c')

		Header_Copy(spec_file_suw_cnt_RS,spec_file_suw_ofn,'h_s_0')
		Header_Copy(spec_file_suw_cnt_RS,spec_file_suw_ofn,'h_s_c')
		
		Header_History_Step(spec_file_sum_ofn,spec_file_sum_cnt_RS)
		Header_History_Step(spec_file_avg_ofn,spec_file_avg_cnt_RS)
		Header_History_Step(spec_file_med_ofn,spec_file_med_cnt_RS)
		Header_History_Step(spec_file_std_ofn,spec_file_std_cnt_RS)
		Header_History_Step(spec_file_rms_ofn,spec_file_rms_cnt_RS)
		Header_History_Step(spec_file_avw_ofn,spec_file_avw_cnt_RS)
		Header_History_Step(spec_file_suw_ofn,spec_file_suw_cnt_RS)

		Header_History_Step_Backwards_Loop(spec_file_sum_cnt_RS,spec_file_sum_ofn)
		Header_History_Step_Backwards_Loop(spec_file_avg_cnt_RS,spec_file_avg_ofn)
		Header_History_Step_Backwards_Loop(spec_file_med_cnt_RS,spec_file_med_ofn)
		Header_History_Step_Backwards_Loop(spec_file_std_cnt_RS,spec_file_std_ofn)
		Header_History_Step_Backwards_Loop(spec_file_rms_cnt_RS,spec_file_rms_ofn)
		Header_History_Step_Backwards_Loop(spec_file_avw_cnt_RS,spec_file_avw_ofn)
		Header_History_Step_Backwards_Loop(spec_file_suw_cnt_RS,spec_file_suw_ofn)

		spec_file_sum_cnt_smt = Spectra_Smooth(spec_file_sum_cnt_RS,smt_shp_pst,smt_sze_pst)
		spec_file_avg_cnt_smt = Spectra_Smooth(spec_file_avg_cnt_RS,smt_shp_pst,smt_sze_pst)
		spec_file_med_cnt_smt = Spectra_Smooth(spec_file_med_cnt_RS,smt_shp_pst,smt_sze_pst)
		spec_file_avw_cnt_smt = Spectra_Smooth(spec_file_avw_cnt_RS,smt_shp_pst,smt_sze_pst)
		spec_file_suw_cnt_smt = Spectra_Smooth(spec_file_suw_cnt_RS,smt_shp_pst,smt_sze_pst)


		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw,'STK_NUM',header_comment='Number of galaxies used for Stack')

		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_sum_cnt_smt,spec_file_sum_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_avg_cnt_smt,spec_file_avg_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_med_cnt_smt,spec_file_med_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_avw_cnt_smt,spec_file_avw_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_suw_cnt_smt,spec_file_avw_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')


		OPT_STCK_CRE = [
						spec_file_sum,spec_file_avg,spec_file_med,spec_file_hst,
						spec_file_std,
						spec_file_rms,
						spec_file_hsw,spec_file_wsu,spec_file_suw,spec_file_avw,
						spec_file_sum_cnt[0],spec_file_sum_cnt[1],
						spec_file_avg_cnt[0],spec_file_avg_cnt[1],
						spec_file_med_cnt[0],spec_file_med_cnt[1],
						spec_file_std_cnt[0],spec_file_std_cnt[1],
						spec_file_rms_cnt[0],spec_file_rms_cnt[1],
						spec_file_avw_cnt[0],spec_file_avw_cnt[1],
						spec_file_suw_cnt[0],spec_file_suw_cnt[1],
						spec_file_sum_smt,spec_file_avg_smt,spec_file_med_smt,
						spec_file_avw_smt,spec_file_suw_smt,
						spec_file_sum_cnt_smt,spec_file_avg_cnt_smt,spec_file_med_cnt_smt,
						spec_file_avw_cnt_smt,spec_file_suw_cnt_smt]


		OPT_STCK_CNT = [
						spec_file_sum_cnt[0],
						spec_file_avg_cnt[0],
						spec_file_med_cnt[0],
						spec_file_std_cnt[0],
						spec_file_rms_cnt[0],
						spec_file_avw_cnt[0],
						spec_file_suw_cnt[0]]

		OPT_STCK_SMT = [
						spec_file_sum_smt,spec_file_avg_smt,spec_file_med_smt,
						spec_file_avw_smt,spec_file_suw_smt,
						spec_file_sum_cnt_smt,spec_file_avg_cnt_smt,spec_file_med_cnt_smt,
						spec_file_avw_cnt_smt,spec_file_suw_cnt_smt]

		OPT_STCK_WGT = [
						spec_file_hsw,spec_file_wsu,spec_file_suw,spec_file_avw,
						spec_file_avw_cnt[0],spec_file_avw_cnt[1],
						spec_file_suw_cnt[0],spec_file_suw_cnt[1],
						spec_file_avw_smt,spec_file_suw_smt,
						spec_file_avw_cnt_smt,spec_file_suw_cnt_smt]

		OPT_STCK_OPR = [
						'SUM','AVG','MED','HST',
						'STD',
						'RMS',
						'HSW','WSU','SUW','AVW',
						'SUC','SUF',
						'AVC','AVF',
						'MEC','MEF',
						'STC','STF',
						'RMC','RMF',
						'AWC','AWF',
						'SWC','SWF',
						'SUMS','AVGS','MEDS',
						'AVWS','SUWS',
						'SUCS','AVCS','MECS',
						'AWCS','SWCS']
		FNL_SPEC_RES = OPT_STCK_CRE
	#########################CORE-FILES###########################
	#######################POST-PROCESSING########################

	#######################POST-PROCESSING########################
	##########################PCT-FILES###########################
	if stk_pct_mde == True and pst_cnt==False and smt_spc_pst == False:
		print
		print colored('Generating Percentile composite files','yellow')
		print colored('Generating Percentile composite files','yellow')
		print colored('2-PCT','yellow')
		print
		pass

		OPT_STCK_PCT      = [
							spec_file_p25,spec_file_p75,
							spec_file_1sl,spec_file_1sh,
							spec_file_2sl,spec_file_2sh,
							spec_file_3sl,spec_file_3sh
							]

		OPT_PCT_OPR        = [
							'P25','P75',
							'1SL','1SH',
							'2SL','2SH',
							'3SL','3SH'
							]
					
		FNL_SPEC_RES_PCT  = OPT_STCK_PCT
	elif stk_pct_mde == True and pst_cnt==True and smt_spc_pst == False:
		print colored('Generating Percentile composite files','yellow')
		print colored('3-PCT','yellow')
		print
		spec_file_p25_cnt  = Spectra_Cont_IRAF(spec_file_p25_ofn,ind_stk_res + 'log_cont-stk-sum',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_p75_cnt  = Spectra_Cont_IRAF(spec_file_p75_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_1sl_cnt  = Spectra_Cont_IRAF(spec_file_1sl_ofn,ind_stk_res + 'log_cont-stk-med',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_1sh_cnt  = Spectra_Cont_IRAF(spec_file_1sh_ofn,ind_stk_res + 'log_cont-stk-std',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_2sl_cnt  = Spectra_Cont_IRAF(spec_file_2sl_ofn,ind_stk_res + 'log_cont-stk-rms',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_2sh_cnt  = Spectra_Cont_IRAF(spec_file_2sh_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_3sl_cnt  = Spectra_Cont_IRAF(spec_file_3sl_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_3sh_cnt  = Spectra_Cont_IRAF(spec_file_3sh_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)

		spec_file_p25_cnt_RS = spec_file_p25_cnt[0]
		spec_file_p75_cnt_RS = spec_file_p75_cnt[0]
		spec_file_1sl_cnt_RS = spec_file_1sl_cnt[0]
		spec_file_1sh_cnt_RS = spec_file_1sh_cnt[0]
		spec_file_2sl_cnt_RS = spec_file_2sl_cnt[0]
		spec_file_2sh_cnt_RS = spec_file_2sh_cnt[0]
		spec_file_3sl_cnt_RS = spec_file_3sl_cnt[0]
		spec_file_3sh_cnt_RS = spec_file_3sh_cnt[0]

		Header_Copy(spec_file_p25_cnt_RS,spec_file_p25_ofn,'h_s_0')
		Header_Copy(spec_file_p25_cnt_RS,spec_file_p25_ofn,'h_s_c')

		Header_Copy(spec_file_p75_cnt_RS,spec_file_p75_ofn,'h_s_0')
		Header_Copy(spec_file_p75_cnt_RS,spec_file_p75_ofn,'h_s_c')

		Header_Copy(spec_file_1sl_cnt_RS,spec_file_1sl_ofn,'h_s_0')
		Header_Copy(spec_file_1sl_cnt_RS,spec_file_1sl_ofn,'h_s_c')

		Header_Copy(spec_file_1sh_cnt_RS,spec_file_1sh_ofn,'h_s_0')
		Header_Copy(spec_file_1sh_cnt_RS,spec_file_1sh_ofn,'h_s_c')

		Header_Copy(spec_file_2sl_cnt_RS,spec_file_2sl_ofn,'h_s_0')
		Header_Copy(spec_file_2sl_cnt_RS,spec_file_2sl_ofn,'h_s_c')

		Header_Copy(spec_file_2sh_cnt_RS,spec_file_2sh_ofn,'h_s_0')
		Header_Copy(spec_file_2sh_cnt_RS,spec_file_2sh_ofn,'h_s_c')

		Header_Copy(spec_file_3sl_cnt_RS,spec_file_3sl_ofn,'h_s_0')
		Header_Copy(spec_file_3sl_cnt_RS,spec_file_3sl_ofn,'h_s_c')

		Header_Copy(spec_file_3sh_cnt_RS,spec_file_3sh_ofn,'h_s_0')
		Header_Copy(spec_file_3sh_cnt_RS,spec_file_3sh_ofn,'h_s_c')
		
		Header_History_Step(spec_file_sum_ofn,spec_file_sum_cnt_RS)
		Header_History_Step(spec_file_avg_ofn,spec_file_avg_cnt_RS)
		Header_History_Step(spec_file_med_ofn,spec_file_med_cnt_RS)
		Header_History_Step(spec_file_std_ofn,spec_file_std_cnt_RS)
		Header_History_Step(spec_file_rms_ofn,spec_file_rms_cnt_RS)
		Header_History_Step(spec_file_avw_ofn,spec_file_avw_cnt_RS)
		Header_History_Step(spec_file_suw_ofn,spec_file_suw_cnt_RS)
		Header_History_Step(spec_file_suw_ofn,spec_file_suw_cnt_RS)

		Header_History_Step_Backwards_Loop(spec_file_p25_cnt_RS,spec_file_p25_ofn)
		Header_History_Step_Backwards_Loop(spec_file_p75_cnt_RS,spec_file_p75_ofn)
		Header_History_Step_Backwards_Loop(spec_file_1sl_cnt_RS,spec_file_1sl_ofn)
		Header_History_Step_Backwards_Loop(spec_file_1sh_cnt_RS,spec_file_1sh_ofn)
		Header_History_Step_Backwards_Loop(spec_file_2sl_cnt_RS,spec_file_2sl_ofn)
		Header_History_Step_Backwards_Loop(spec_file_2sh_cnt_RS,spec_file_2sh_ofn)
		Header_History_Step_Backwards_Loop(spec_file_3sl_cnt_RS,spec_file_3sl_ofn)
		Header_History_Step_Backwards_Loop(spec_file_3sh_cnt_RS,spec_file_3sh_ofn)

		OPT_STCK_PCT      = [
							spec_file_p25,spec_file_p75,
							spec_file_1sl,spec_file_1sh,
							spec_file_2sl,spec_file_2sh,
							spec_file_3sl,spec_file_3sh,
							spec_file_p25_cnt[0],spec_file_p25_cnt[1],
							spec_file_p75_cnt[0],spec_file_p75_cnt[1],
							spec_file_1sl_cnt[0],spec_file_1sl_cnt[1],
							spec_file_1sh_cnt[0],spec_file_1sh_cnt[1],
							spec_file_2sl_cnt[0],spec_file_2sl_cnt[1],
							spec_file_2sh_cnt[0],spec_file_2sh_cnt[1],
							spec_file_3sl_cnt[0],spec_file_3sl_cnt[1],
							spec_file_3sh_cnt[0],spec_file_3sh_cnt[1]
							]

		OPT_PCT_CNT      = [
							spec_file_p25_cnt[0],
							spec_file_p75_cnt[0],
							spec_file_1sl_cnt[0],
							spec_file_1sh_cnt[0],
							spec_file_2sl_cnt[0],
							spec_file_2sh_cnt[0],
							spec_file_3sl_cnt[0],
							spec_file_3sh_cnt[0]
							]							

		OPT_PCT_OPR        = [
							'P25','P75',
							'1SL','1SH',
							'2SL','2SH',
							'3SL','3SH',
							'P5C','P5F',
							'P5C','P5F',
							'1LC','1LF',
							'1HC','1HF',
							'2LC','2LF',
							'2HC','2HF',
							'3LC','3LF',
							'3HC','3HF',
							]
		FNL_SPEC_RES_PCT = OPT_STCK_PCT

	elif stk_pct_mde == True and pst_cnt==False and smt_spc_pst == True:
		print
		print colored('Generating Percentile composite files','yellow')
		print colored('4-PCT','yellow')
		print
		spec_file_p25_smt     = Spectra_Smooth(spec_file_p25_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_p75_smt     = Spectra_Smooth(spec_file_p75_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_1sl_smt     = Spectra_Smooth(spec_file_1sl_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_1sh_smt     = Spectra_Smooth(spec_file_1sh_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_2sl_smt     = Spectra_Smooth(spec_file_2sl_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_2sh_smt     = Spectra_Smooth(spec_file_2sh_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_3sl_smt     = Spectra_Smooth(spec_file_3sl_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_3sh_smt     = Spectra_Smooth(spec_file_3sh_ofn,smt_shp_pst,smt_sze_pst)						

		try:
			print spec_file_p25_smt,Header_Get(spec_file_p25_smt,'STK_NUM')
			print spec_file_p25,Header_Get(spec_file_p25,'STK_NUM')
			print spec_file_p75_smt,Header_Get(spec_file_p75_smt,'STK_NUM')
			print spec_file_p75,Header_Get(spec_file_p75,'STK_NUM')
			print spec_file_1sl_smt,Header_Get(spec_file_1sl_smt,'STK_NUM')
			print spec_file_1sl,Header_Get(spec_file_1sl,'STK_NUM')
			print spec_file_1sh_smt,Header_Get(spec_file_1sh_smt,'STK_NUM')
			print spec_file_1sh,Header_Get(spec_file_1sh,'STK_NUM')
			print spec_file_2sl_smt,Header_Get(spec_file_2sl_smt,'STK_NUM')
			print spec_file_2sl,Header_Get(spec_file_2sl,'STK_NUM')
			print spec_file_2sh_smt,Header_Get(spec_file_2sh_smt,'STK_NUM')
			print spec_file_2sh,Header_Get(spec_file_2sh,'STK_NUM')
			print spec_file_3sl_smt,Header_Get(spec_file_3sl_smt,'STK_NUM')
			print spec_file_3sl,Header_Get(spec_file_3sl,'STK_NUM')
			print spec_file_3sh_smt,Header_Get(spec_file_3sh_smt,'STK_NUM')
			print spec_file_3sh,Header_Get(spec_file_3sh,'STK_NUM')

		except KeyError:
			pass

		Header_Copy(spec_file_p25_smt,spec_file_p25,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_p75_smt,spec_file_p75,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_1sl_smt,spec_file_1sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_1sh_smt,spec_file_1sh,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_2sl_smt,spec_file_2sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_2sh_smt,spec_file_2sh,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_3sl_smt,spec_file_3sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_3sh_smt,spec_file_3sh,'STK_NUM',header_comment='Number of galaxies used for Stack')		

		print spec_file_p25_smt,Header_Get(spec_file_p25_smt,'STK_NUM')
		print spec_file_p25,Header_Get(spec_file_p25,'STK_NUM')
		print spec_file_p75_smt,Header_Get(spec_file_p75_smt,'STK_NUM')
		print spec_file_p75,Header_Get(spec_file_p75,'STK_NUM')
		print spec_file_1sl_smt,Header_Get(spec_file_1sl_smt,'STK_NUM')
		print spec_file_1sl,Header_Get(spec_file_1sl,'STK_NUM')
		print spec_file_1sh_smt,Header_Get(spec_file_1sh_smt,'STK_NUM')
		print spec_file_1sh,Header_Get(spec_file_1sh,'STK_NUM')
		print spec_file_2sl_smt,Header_Get(spec_file_2sl_smt,'STK_NUM')
		print spec_file_2sl,Header_Get(spec_file_2sl,'STK_NUM')
		print spec_file_2sh_smt,Header_Get(spec_file_2sh_smt,'STK_NUM')
		print spec_file_2sh,Header_Get(spec_file_2sh,'STK_NUM')
		print spec_file_3sl_smt,Header_Get(spec_file_3sl_smt,'STK_NUM')
		print spec_file_3sl,Header_Get(spec_file_3sl,'STK_NUM')
		print spec_file_3sh_smt,Header_Get(spec_file_3sh_smt,'STK_NUM')
		print spec_file_3sh,Header_Get(spec_file_3sh,'STK_NUM')

		OPT_STCK_PCT      = [
							spec_file_p25,spec_file_p75,
							spec_file_1sl,spec_file_1sh,
							spec_file_2sl,spec_file_2sh,
							spec_file_3sl,spec_file_3sh,
							spec_file_p25_smt,spec_file_p75_smt,
							spec_file_1sl_smt,spec_file_1sh_smt,
							spec_file_2sl_smt,spec_file_2sh_smt,
							spec_file_3sl_smt,spec_file_3sh_smt
							]

		OPT_PCT_OPR        = [
							'P25','P75',
							'1SL','1SH',
							'2SL','2SH',
							'3SL','3SH',
							'P2S','P7S',
							'1LS','1HS',
							'2LS','2HS',
							'3LS','3HS'

							]
		FNL_SPEC_RES_PCT  = OPT_STCK_PCT

		OPT_PCT_SMT      = [
							spec_file_p25_smt,spec_file_p75_smt,
							spec_file_1sl_smt,spec_file_1sh_smt,
							spec_file_2sl_smt,spec_file_2sh_smt,
							spec_file_3sl_smt,spec_file_3sh_smt,
							]
	elif stk_pct_mde == True and pst_cnt==True and smt_spc_pst == True:
		print
		print colored('Generating Percentile composite files','yellow')
		print colored('5-PCT','yellow')
		print

		spec_file_p25_smt     = Spectra_Smooth(spec_file_p25_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_p75_smt     = Spectra_Smooth(spec_file_p75_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_1sl_smt     = Spectra_Smooth(spec_file_1sl_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_1sh_smt     = Spectra_Smooth(spec_file_1sh_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_2sl_smt     = Spectra_Smooth(spec_file_2sl_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_2sh_smt     = Spectra_Smooth(spec_file_2sh_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_3sl_smt     = Spectra_Smooth(spec_file_3sl_ofn,smt_shp_pst,smt_sze_pst)
		spec_file_3sh_smt     = Spectra_Smooth(spec_file_3sh_ofn,smt_shp_pst,smt_sze_pst)						


		Header_Copy(spec_file_p25_smt,spec_file_p25,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_p75_smt,spec_file_p75,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_1sl_smt,spec_file_1sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_1sh_smt,spec_file_1sh,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_2sl_smt,spec_file_2sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_2sh_smt,spec_file_2sh,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_3sl_smt,spec_file_3sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_3sh_smt,spec_file_3sh,'STK_NUM',header_comment='Number of galaxies used for Stack')		

		try:
			print spec_file_p25_smt,Header_Get(spec_file_p25_smt,'STK_NUM')
		except KeyError:
			pass
			print colored('No Header! STK_NUM','yellow')
			print
		print spec_file_p25,Header_Get(spec_file_p25,'STK_NUM')
		try:
			print spec_file_p75_smt,Header_Get(spec_file_p75_smt,'STK_NUM')
		except KeyError:
			pass
			print colored('No Header! STK_NUM','yellow')
			print
		print spec_file_p75,Header_Get(spec_file_p75,'STK_NUM')
		try:
			print spec_file_1sl_smt,Header_Get(spec_file_1sl_smt,'STK_NUM')
		except KeyError:
			pass
			print colored('No Header! STK_NUM','yellow')
			print
		print spec_file_1sl,Header_Get(spec_file_1sl,'STK_NUM')
		try:
			print spec_file_1sh_smt,Header_Get(spec_file_1sh_smt,'STK_NUM')
		except KeyError:
			pass
			print colored('No Header! STK_NUM','yellow')
			print
		print spec_file_1sh,Header_Get(spec_file_1sh,'STK_NUM')
		try:
			print spec_file_2sl_smt,Header_Get(spec_file_2sl_smt,'STK_NUM')
		except KeyError:
			pass
			print colored('No Header! STK_NUM','yellow')
			print
		print spec_file_2sl,Header_Get(spec_file_2sl,'STK_NUM')

		try:
			print spec_file_2sh_smt,Header_Get(spec_file_2sh_smt,'STK_NUM')
		except KeyError:
			pass
			print colored('No Header! STK_NUM','yellow')
			print
		print spec_file_2sh,Header_Get(spec_file_2sh,'STK_NUM')
		try:
			print spec_file_3sl_smt,Header_Get(spec_file_3sl_smt,'STK_NUM')
		except KeyError:
			pass
			print colored('No Header! STK_NUM','yellow')
			print
		print spec_file_3sl,Header_Get(spec_file_3sl,'STK_NUM')
		try:
			print spec_file_3sh_smt,Header_Get(spec_file_3sh_smt,'STK_NUM')
		except KeyError:
			pass
			print colored('No Header! STK_NUM','yellow')
			print
		print spec_file_3sh,Header_Get(spec_file_3sh,'STK_NUM')	

		print spec_file_p25_smt,Header_Get(spec_file_p25_smt,'STK_NUM')
		print spec_file_p25,Header_Get(spec_file_p25,'STK_NUM')
		print spec_file_p75_smt,Header_Get(spec_file_p75_smt,'STK_NUM')
		print spec_file_p75,Header_Get(spec_file_p75,'STK_NUM')
		print spec_file_1sl_smt,Header_Get(spec_file_1sl_smt,'STK_NUM')
		print spec_file_1sl,Header_Get(spec_file_1sl,'STK_NUM')
		print spec_file_1sh_smt,Header_Get(spec_file_1sh_smt,'STK_NUM')
		print spec_file_1sh,Header_Get(spec_file_1sh,'STK_NUM')
		print spec_file_2sl_smt,Header_Get(spec_file_2sl_smt,'STK_NUM')
		print spec_file_2sl,Header_Get(spec_file_2sl,'STK_NUM')
		print spec_file_2sh_smt,Header_Get(spec_file_2sh_smt,'STK_NUM')
		print spec_file_2sh,Header_Get(spec_file_2sh,'STK_NUM')
		print spec_file_3sl_smt,Header_Get(spec_file_3sl_smt,'STK_NUM')
		print spec_file_3sl,Header_Get(spec_file_3sl,'STK_NUM')
		print spec_file_3sh_smt,Header_Get(spec_file_3sh_smt,'STK_NUM')
		print spec_file_3sh,Header_Get(spec_file_3sh,'STK_NUM')	

		spec_file_p25_cnt  = Spectra_Cont_IRAF(spec_file_p25_ofn,ind_stk_res + 'log_cont-stk-sum',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_p75_cnt  = Spectra_Cont_IRAF(spec_file_p75_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_1sl_cnt  = Spectra_Cont_IRAF(spec_file_1sl_ofn,ind_stk_res + 'log_cont-stk-med',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_1sh_cnt  = Spectra_Cont_IRAF(spec_file_1sh_ofn,ind_stk_res + 'log_cont-stk-std',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_2sl_cnt  = Spectra_Cont_IRAF(spec_file_2sl_ofn,ind_stk_res + 'log_cont-stk-rms',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_2sh_cnt  = Spectra_Cont_IRAF(spec_file_2sh_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_3sl_cnt  = Spectra_Cont_IRAF(spec_file_3sl_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)
		spec_file_3sh_cnt  = Spectra_Cont_IRAF(spec_file_3sh_ofn,ind_stk_res + 'log_cont-stk-avg',
							Cont_type_IRAF     = pst_cnt_typ,Cont_lines_IRAF    = pst_cnt_lns,
							Cont_funct_IRAF    = pst_cnt_fnc,Cont_order_IRAF    = pst_cnt_ord,
							Cont_override_IRAF = pst_cnt_ovr,Cont_replace_IRAF  = pst_cnt_rpl,
							Cont_low_rej_IRAF  = pst_cnt_lrj,Cont_high_rej_IRAF = pst_cnt_hrj)

		spec_file_p25_cnt_RS = spec_file_p25_cnt[0]
		spec_file_p75_cnt_RS = spec_file_p75_cnt[0]
		spec_file_1sl_cnt_RS = spec_file_1sl_cnt[0]
		spec_file_1sh_cnt_RS = spec_file_1sh_cnt[0]
		spec_file_2sl_cnt_RS = spec_file_2sl_cnt[0]
		spec_file_2sh_cnt_RS = spec_file_2sh_cnt[0]
		spec_file_3sl_cnt_RS = spec_file_3sl_cnt[0]
		spec_file_3sh_cnt_RS = spec_file_3sh_cnt[0]

		Header_Copy(spec_file_p25_cnt_RS,spec_file_p25_ofn,'h_s_0')
		Header_Copy(spec_file_p25_cnt_RS,spec_file_p25_ofn,'h_s_c')

		Header_Copy(spec_file_p75_cnt_RS,spec_file_p75_ofn,'h_s_0')
		Header_Copy(spec_file_p75_cnt_RS,spec_file_p75_ofn,'h_s_c')

		Header_Copy(spec_file_1sl_cnt_RS,spec_file_1sl_ofn,'h_s_0')
		Header_Copy(spec_file_1sl_cnt_RS,spec_file_1sl_ofn,'h_s_c')

		Header_Copy(spec_file_1sh_cnt_RS,spec_file_1sh_ofn,'h_s_0')
		Header_Copy(spec_file_1sh_cnt_RS,spec_file_1sh_ofn,'h_s_c')

		Header_Copy(spec_file_2sl_cnt_RS,spec_file_2sl_ofn,'h_s_0')
		Header_Copy(spec_file_2sl_cnt_RS,spec_file_2sl_ofn,'h_s_c')

		Header_Copy(spec_file_2sh_cnt_RS,spec_file_2sh_ofn,'h_s_0')
		Header_Copy(spec_file_2sh_cnt_RS,spec_file_2sh_ofn,'h_s_c')

		Header_Copy(spec_file_3sl_cnt_RS,spec_file_3sl_ofn,'h_s_0')
		Header_Copy(spec_file_3sl_cnt_RS,spec_file_3sl_ofn,'h_s_c')

		Header_Copy(spec_file_3sh_cnt_RS,spec_file_3sh_ofn,'h_s_0')
		Header_Copy(spec_file_3sh_cnt_RS,spec_file_3sh_ofn,'h_s_c')
		
		Header_History_Step(spec_file_sum_ofn,spec_file_sum_cnt_RS)
		Header_History_Step(spec_file_avg_ofn,spec_file_avg_cnt_RS)
		Header_History_Step(spec_file_med_ofn,spec_file_med_cnt_RS)
		Header_History_Step(spec_file_std_ofn,spec_file_std_cnt_RS)
		Header_History_Step(spec_file_rms_ofn,spec_file_rms_cnt_RS)
		Header_History_Step(spec_file_avw_ofn,spec_file_avw_cnt_RS)
		Header_History_Step(spec_file_suw_ofn,spec_file_suw_cnt_RS)
		Header_History_Step(spec_file_suw_ofn,spec_file_suw_cnt_RS)

		Header_History_Step_Backwards_Loop(spec_file_p25_cnt_RS,spec_file_p25_ofn)
		Header_History_Step_Backwards_Loop(spec_file_p75_cnt_RS,spec_file_p75_ofn)
		Header_History_Step_Backwards_Loop(spec_file_1sl_cnt_RS,spec_file_1sl_ofn)
		Header_History_Step_Backwards_Loop(spec_file_1sh_cnt_RS,spec_file_1sh_ofn)
		Header_History_Step_Backwards_Loop(spec_file_2sl_cnt_RS,spec_file_2sl_ofn)
		Header_History_Step_Backwards_Loop(spec_file_2sh_cnt_RS,spec_file_2sh_ofn)
		Header_History_Step_Backwards_Loop(spec_file_3sl_cnt_RS,spec_file_3sl_ofn)
		Header_History_Step_Backwards_Loop(spec_file_3sh_cnt_RS,spec_file_3sh_ofn)

		spec_file_p25_cnt_smt = Spectra_Smooth(spec_file_p25_cnt[0],smt_shp_pst,smt_sze_pst)
		spec_file_p75_cnt_smt = Spectra_Smooth(spec_file_p75_cnt[0],smt_shp_pst,smt_sze_pst)
		spec_file_1sl_cnt_smt = Spectra_Smooth(spec_file_1sl_cnt[0],smt_shp_pst,smt_sze_pst)
		spec_file_1sh_cnt_smt = Spectra_Smooth(spec_file_1sh_cnt[0],smt_shp_pst,smt_sze_pst)
		spec_file_2sl_cnt_smt = Spectra_Smooth(spec_file_2sl_cnt[0],smt_shp_pst,smt_sze_pst)
		spec_file_2sh_cnt_smt = Spectra_Smooth(spec_file_2sh_cnt[0],smt_shp_pst,smt_sze_pst)
		spec_file_3sl_cnt_smt = Spectra_Smooth(spec_file_3sl_cnt[0],smt_shp_pst,smt_sze_pst)
		spec_file_3sh_cnt_smt = Spectra_Smooth(spec_file_3sh_cnt[0],smt_shp_pst,smt_sze_pst)


		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl,'STK_NUM',header_comment='Number of galaxies used for Stack')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh,'STK_NUM',header_comment='Number of galaxies used for Stack')

		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_p25_cnt_smt,spec_file_p25_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_p75_cnt_smt,spec_file_p75_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_1sl_cnt_smt,spec_file_1sl_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_1sh_cnt_smt,spec_file_1sh_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_2sl_cnt_smt,spec_file_2sl_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_2sh_cnt_smt,spec_file_2sh_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_3sl_cnt_smt,spec_file_3sl_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_TYP',header_comment='Continuum IRAF type')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_FIT',header_comment='Continuum IRAF function')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_ORD',header_comment='Continuum IRAF order')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_RPC',header_comment='Continuum IRAF replacement')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_LRJ',header_comment='Continuum IRAF low sigma rejection')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_HRJ',header_comment='Continuum IRAF high sigma rejection')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_FNM',header_comment='Continuum IRAF Spec file pre-cont')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_FN0',header_comment='Continuum IRAF Spec file pst-cont')
		Header_Copy(spec_file_3sh_cnt_smt,spec_file_3sh_cnt_RS,'CNT_FNF',header_comment='Continuum IRAF Spec file cnt-fit')

		OPT_STCK_PCT      = [
							spec_file_p25,spec_file_p75,
							spec_file_1sl,spec_file_1sh,
							spec_file_2sl,spec_file_2sh,
							spec_file_3sl,spec_file_3sh,
							spec_file_p25_smt,spec_file_p75_smt,
							spec_file_1sl_smt,spec_file_1sh_smt,
							spec_file_2sl_smt,spec_file_2sh_smt,
							spec_file_3sl_smt,spec_file_3sh_smt,
							spec_file_p25_cnt[0],spec_file_p25_cnt[1],
							spec_file_p75_cnt[0],spec_file_p75_cnt[1],
							spec_file_1sl_cnt[0],spec_file_1sl_cnt[1],
							spec_file_1sh_cnt[0],spec_file_1sh_cnt[1],
							spec_file_2sl_cnt[0],spec_file_2sl_cnt[1],
							spec_file_2sh_cnt[0],spec_file_2sh_cnt[1],
							spec_file_3sl_cnt[0],spec_file_3sl_cnt[1],
							spec_file_3sh_cnt[0],spec_file_3sh_cnt[1],
							spec_file_p25_cnt_smt,spec_file_p75_cnt_smt,
							spec_file_1sl_cnt_smt,spec_file_1sh_cnt_smt,
							spec_file_2sl_cnt_smt,spec_file_2sh_cnt_smt,
							spec_file_3sl_cnt_smt,spec_file_3sh_cnt_smt
							]

		OPT_PCT_SMT      = [
							spec_file_p25_smt,spec_file_p75_smt,
							spec_file_1sl_smt,spec_file_1sh_smt,
							spec_file_2sl_smt,spec_file_2sh_smt,
							spec_file_3sl_smt,spec_file_3sh_smt,
							spec_file_p25_cnt_smt,spec_file_p75_cnt_smt,
							spec_file_1sl_cnt_smt,spec_file_1sh_cnt_smt,
							spec_file_2sl_cnt_smt,spec_file_2sh_cnt_smt,
							spec_file_3sl_cnt_smt,spec_file_3sh_cnt_smt
							]

		OPT_PCT_CNT      = [
							spec_file_p25_cnt[0],
							spec_file_p75_cnt[0],
							spec_file_1sl_cnt[0],
							spec_file_1sh_cnt[0],
							spec_file_2sl_cnt[0],
							spec_file_2sh_cnt[0],
							spec_file_3sl_cnt[0],
							spec_file_3sh_cnt[0]
							]

		OPT_PCT_OPR        = [
							'P25','P75',
							'1SL','1SH',
							'2SL','2SH',
							'3SL','3SH',
							'P2S','P7S',
							'1LS','1HS',
							'2LS','2HS',
							'3LS','3HS',
							'P5C','P5F',
							'P5C','P5F',
							'1LC','1LF',
							'1HC','1HF',
							'2LC','2LF',
							'2HC','2HF',
							'3LC','3LF',
							'3HC','3HF',
							'P5CS',
							'P5CS',
							'1LCS',
							'1HCS',
							'2LCS',
							'2HCS',
							'3LCS',
							'3HCS'
							]

		FNL_SPEC_RES_PCT  = OPT_STCK_PCT
	elif stk_pct_mde == False:
		pass
	##########################PCT-FILES###########################
	#######################POST-PROCESSING########################

	if wrt_fits == True:
		print
		print colored('Checking length! 2b commented on Fnc_Stk_Stk.py','yellow')
		print colored('OPT_STCK_CRE, OPT_STCK_OPR','yellow')
		print colored(str(len(OPT_STCK_CRE))+'-'+str(len(OPT_STCK_OPR)),'yellow')
		print
		##########ADDING IMPORTANT HEADERS FOR STACKED FILES##########
		#########################CORE-FILES###########################
		for stck_res_file ,stck_opr in zip(OPT_STCK_CRE,OPT_STCK_OPR):
			Header_Get_Add(stck_res_file,'STK_OPR',str(stck_opr),header_comment='Stack operation')
			Header_Get_Add(stck_res_file,'MAG_I'  ,0            ,header_comment='i-band magnitude')
			Header_Get_Add(stck_res_file,'Z_0'    ,0            ,header_comment='Redshift')
			Header_Get_Add(stck_res_file,'Z_REF'  ,0            ,header_comment='Redshift reference')
			Header_Get_Add(stck_res_file,'SEL_SHF',sel_pre_shf  ,header_comment='Redshift reference')
			Header_Get_Add(stck_res_file,'SEL_CNT',sel_pre_cnt  ,header_comment='Redshift reference')
			Header_Get_Add(stck_res_file,'SEL_MSK',sel_pre_msk  ,header_comment='Redshift reference')

			if pre_cnt == True:
				Header_Get_Add(stck_res_file,'PCT_PCT',str(pre_cnt)    ,header_comment='Pre-Continuum Fit')
				Header_Get_Add(stck_res_file,'PCT_TYP',str(pre_cnt_typ),header_comment='Continuum IRAF type')
				Header_Get_Add(stck_res_file,'PCT_FIT',str(pre_cnt_fnc),header_comment='Continuum IRAF function')
				Header_Get_Add(stck_res_file,'PCT_ORD',str(pre_cnt_ovr),header_comment='Continuum IRAF order')
				Header_Get_Add(stck_res_file,'PCT_RPC',str(pre_cnt_rpl),header_comment='Continuum IRAF replacement')
				Header_Get_Add(stck_res_file,'PCT_LRJ',str(pre_cnt_lrj),header_comment='Continuum IRAF low sigma rejection')
				Header_Get_Add(stck_res_file,'PCT_HRJ',str(pre_cnt_hrj),header_comment='Continuum IRAF high sigma rejection')
			elif pre_cnt == False:
				Header_Get_Add(stck_res_file,'PCT_PCT',str(pre_cnt),header_comment='Pre-Continuum Fit')
				pass

			if smt_spc_pre == True:
				Header_Get_Add(stck_res_file,'PSM_PSM',str(smt_spc_pre)    ,header_comment='Pre-Smooth')
				Header_Get_Add(stck_res_file,'PSM_KTP',str(smt_shp_pre)    ,header_comment='Pre-Smooth Kernel Type')
				Header_Get_Add(stck_res_file,'PSM_KSZ',str(smt_sze_pre)    ,header_comment='Pre-Smooth Kernel Size')
			elif smt_spc_pre == False:
				Header_Get_Add(stck_res_file,'PSM_PSM',str(smt_spc_pre)    ,header_comment='Pre-Smooth')

			if pre_msk == True:			
				Header_Get_Add(stck_res_file,'PMK_PMK',str(pre_msk)        ,header_comment='Pre-Smooth')
				Header_Get_Add(stck_res_file,'PMK_MTP',str(pre_msk_typ)    ,header_comment='Pre-Smooth Mask Type')
				Header_Get_Add(stck_res_file,'PMK_CVL',str(pre_msk_cte_val),header_comment='Pre-Smooth Constant Value')
				Header_Get_Add(stck_res_file,'PMK_ABL',str(pre_msk_abs_lne),header_comment='Pre-Smooth Mask Absorption Lines')
				Header_Get_Add(stck_res_file,'PMK_BRG',str(pre_msk_rgn)    ,header_comment='Pre-Smooth Mask Spec Region')
				Header_Get_Add(stck_res_file,'PMK_BMN',str(pre_msk_min)    ,header_comment='Pre-Smooth Mask Spec Region Min')
				Header_Get_Add(stck_res_file,'PMK_BMX',str(pre_msk_max)    ,header_comment='Pre-Smooth Mask Spec Region Max')
			elif pre_msk == False:
				Header_Get_Add(stck_res_file,'PMK_PMK',str(pre_msk)        ,header_comment='Pre-Smooth')

			if sig_clp == True:
				Header_Get_Add(stck_res_file,'SCP_SCP',str(sig_clp),header_comment='Sigma Clipping')
				Header_Get_Add(stck_res_file,'SCP_CUT',int(sig_cut),header_comment='Sigma Clipping low/high rejection Values')
				Header_Get_Add(stck_res_file,'SCP_CFN',str(sig_fct),header_comment='Sigma Clipping Central Function')
				Header_Get_Add(stck_res_file,'SCP_MFV',str(sig_fll),header_comment='Sigma Clipping Filling Value')
			elif sig_clp == False:
				Header_Get_Add(stck_res_file,'SCP_SCP',str(sig_clp),header_comment='Sigma Clipping')
				pass
			if  '-BS-' in name:
				pass
			elif  '-BS_MST' in name:
				Header_Get_Add(stck_res_file,'BST_NMB',str(bs_nmb_itr)  ,header_comment='Bootstrap Repetitions')
			else:
				pass
		##########ADDING IMPORTANT HEADERS FOR STACKED FILES##########
		#########################CORE-FILES###########################

		if stk_wgt_mde == True:
			print
			print colored('Checking length! 2b commented on Fnc_Stk_Stk.py','yellow')
			print colored('OPT_STCK_CRE, OPT_STCK_OPR','yellow')
			print colored(str(len(OPT_WGHT_FLS))+'-'+str(len(OPT_WGT_OPR)),'yellow')
			print
			for stck_res_file ,stck_opr in zip(OPT_WGHT_FLS,OPT_WGT_OPR):
				Header_Get_Add(stck_res_file,'STK_OPR',str(stck_opr),header_comment='Stack operation')
				Header_Get_Add(stck_res_file,'MAG_I'  ,0            ,header_comment='i-band magnitude')
				Header_Get_Add(stck_res_file,'Z_0'    ,0            ,header_comment='Redshift')
				Header_Get_Add(stck_res_file,'Z_REF'  ,0            ,header_comment='Redshift reference')
				Header_Get_Add(stck_res_file,'SEL_SHF',sel_pre_shf  ,header_comment='Redshift reference')
				Header_Get_Add(stck_res_file,'SEL_CNT',sel_pre_cnt  ,header_comment='Redshift reference')
				Header_Get_Add(stck_res_file,'SEL_MSK',sel_pre_msk  ,header_comment='Redshift reference')

				if pre_cnt == True:
					Header_Get_Add(stck_res_file,'PCT_PCT',str(pre_cnt)    ,header_comment='Pre-Continuum Fit')
					Header_Get_Add(stck_res_file,'PCT_TYP',str(pre_cnt_typ),header_comment='Continuum IRAF type')
					Header_Get_Add(stck_res_file,'PCT_FIT',str(pre_cnt_fnc),header_comment='Continuum IRAF function')
					Header_Get_Add(stck_res_file,'PCT_ORD',str(pre_cnt_ovr),header_comment='Continuum IRAF order')
					Header_Get_Add(stck_res_file,'PCT_RPC',str(pre_cnt_rpl),header_comment='Continuum IRAF replacement')
					Header_Get_Add(stck_res_file,'PCT_LRJ',str(pre_cnt_lrj),header_comment='Continuum IRAF low sigma rejection')
					Header_Get_Add(stck_res_file,'PCT_HRJ',str(pre_cnt_hrj),header_comment='Continuum IRAF high sigma rejection')
				elif pre_cnt == False:
					Header_Get_Add(stck_res_file,'PCT_PCT',str(pre_cnt),header_comment='Pre-Continuum Fit')
					pass

				if smt_spc_pre == True:
					Header_Get_Add(stck_res_file,'PSM_PSM',str(smt_spc_pre)    ,header_comment='Pre-Smooth')
					Header_Get_Add(stck_res_file,'PSM_KTP',str(smt_shp_pre)    ,header_comment='Pre-Smooth Kernel Type')
					Header_Get_Add(stck_res_file,'PSM_KSZ',str(smt_sze_pre)    ,header_comment='Pre-Smooth Kernel Size')
				elif smt_spc_pre == False:
					Header_Get_Add(stck_res_file,'PSM_PSM',str(smt_spc_pre)    ,header_comment='Pre-Smooth')

				if pre_msk == True:			
					Header_Get_Add(stck_res_file,'PMK_PMK',str(pre_msk)        ,header_comment='Pre-Smooth')
					Header_Get_Add(stck_res_file,'PMK_MTP',str(pre_msk_typ)    ,header_comment='Pre-Smooth Mask Type')
					Header_Get_Add(stck_res_file,'PMK_CVL',str(pre_msk_cte_val),header_comment='Pre-Smooth Constant Value')
					Header_Get_Add(stck_res_file,'PMK_ABL',str(pre_msk_abs_lne),header_comment='Pre-Smooth Mask Absorption Lines')
					Header_Get_Add(stck_res_file,'PMK_BRG',str(pre_msk_rgn)    ,header_comment='Pre-Smooth Mask Spec Region')
					Header_Get_Add(stck_res_file,'PMK_BMN',str(pre_msk_min)    ,header_comment='Pre-Smooth Mask Spec Region Min')
					Header_Get_Add(stck_res_file,'PMK_BMX',str(pre_msk_max)    ,header_comment='Pre-Smooth Mask Spec Region Max')
				elif pre_msk == False:
					Header_Get_Add(stck_res_file,'PMK_PMK',str(pre_msk)        ,header_comment='Pre-Smooth')

				if sig_clp == True:
					Header_Get_Add(stck_res_file,'SCP_SCP',str(sig_clp),header_comment='Sigma Clipping')
					Header_Get_Add(stck_res_file,'SCP_CUT',int(sig_cut),header_comment='Sigma Clipping low/high rejection Values')
					Header_Get_Add(stck_res_file,'SCP_CFN',str(sig_fct),header_comment='Sigma Clipping Central Function')
					Header_Get_Add(stck_res_file,'SCP_MFV',str(sig_fll),header_comment='Sigma Clipping Filling Value')
				elif sig_clp == False:
					Header_Get_Add(stck_res_file,'SCP_SCP',str(sig_clp),header_comment='Sigma Clipping')
					pass
				if  '-BS-' in name:
					pass
				elif  '-BS_MST' in name:
					Header_Get_Add(stck_res_file,'BST_NMB',str(bs_nmb_itr)  ,header_comment='Bootstrap Repetitions')
				else:
					pass
			##
		else:
			pass
		if stk_pct_mde == True:
			print
			print colored('Checking length! 2b commented on Fnc_Stk_Stk.py','yellow')
			print colored('OPT_STCK_CRE, OPT_STCK_OPR','yellow')
			print colored(str(len(OPT_STCK_PCT))+'-'+str(len(OPT_PCT_OPR)),'yellow')
			print
			for stck_res_file ,stck_opr in zip(OPT_STCK_PCT,OPT_PCT_OPR):
				Header_Get_Add(stck_res_file,'STK_OPR',str(stck_opr),header_comment='Stack operation')
				Header_Get_Add(stck_res_file,'MAG_I'  ,0            ,header_comment='i-band magnitude')
				Header_Get_Add(stck_res_file,'Z_0'    ,0            ,header_comment='Redshift')
				Header_Get_Add(stck_res_file,'Z_REF'  ,0            ,header_comment='Redshift reference')
				Header_Get_Add(stck_res_file,'SEL_SHF',sel_pre_shf  ,header_comment='Redshift reference')
				Header_Get_Add(stck_res_file,'SEL_CNT',sel_pre_cnt  ,header_comment='Redshift reference')
				Header_Get_Add(stck_res_file,'SEL_MSK',sel_pre_msk  ,header_comment='Redshift reference')

				if pre_cnt == True:
					Header_Get_Add(stck_res_file,'PCT_PCT',str(pre_cnt)    ,header_comment='Pre-Continuum Fit')
					Header_Get_Add(stck_res_file,'PCT_TYP',str(pre_cnt_typ),header_comment='Continuum IRAF type')
					Header_Get_Add(stck_res_file,'PCT_FIT',str(pre_cnt_fnc),header_comment='Continuum IRAF function')
					Header_Get_Add(stck_res_file,'PCT_ORD',str(pre_cnt_ovr),header_comment='Continuum IRAF order')
					Header_Get_Add(stck_res_file,'PCT_RPC',str(pre_cnt_rpl),header_comment='Continuum IRAF replacement')
					Header_Get_Add(stck_res_file,'PCT_LRJ',str(pre_cnt_lrj),header_comment='Continuum IRAF low sigma rejection')
					Header_Get_Add(stck_res_file,'PCT_HRJ',str(pre_cnt_hrj),header_comment='Continuum IRAF high sigma rejection')
				elif pre_cnt == False:
					Header_Get_Add(stck_res_file,'PCT_PCT',str(pre_cnt),header_comment='Pre-Continuum Fit')
					pass

				if smt_spc_pre == True:
					Header_Get_Add(stck_res_file,'PSM_PSM',str(smt_spc_pre)    ,header_comment='Pre-Smooth')
					Header_Get_Add(stck_res_file,'PSM_KTP',str(smt_shp_pre)    ,header_comment='Pre-Smooth Kernel Type')
					Header_Get_Add(stck_res_file,'PSM_KSZ',str(smt_sze_pre)    ,header_comment='Pre-Smooth Kernel Size')
				elif smt_spc_pre == False:
					Header_Get_Add(stck_res_file,'PSM_PSM',str(smt_spc_pre)    ,header_comment='Pre-Smooth')

				if pre_msk == True:			
					Header_Get_Add(stck_res_file,'PMK_PMK',str(pre_msk)        ,header_comment='Pre-Smooth')
					Header_Get_Add(stck_res_file,'PMK_MTP',str(pre_msk_typ)    ,header_comment='Pre-Smooth Mask Type')
					Header_Get_Add(stck_res_file,'PMK_CVL',str(pre_msk_cte_val),header_comment='Pre-Smooth Constant Value')
					Header_Get_Add(stck_res_file,'PMK_ABL',str(pre_msk_abs_lne),header_comment='Pre-Smooth Mask Absorption Lines')
					Header_Get_Add(stck_res_file,'PMK_BRG',str(pre_msk_rgn)    ,header_comment='Pre-Smooth Mask Spec Region')
					Header_Get_Add(stck_res_file,'PMK_BMN',str(pre_msk_min)    ,header_comment='Pre-Smooth Mask Spec Region Min')
					Header_Get_Add(stck_res_file,'PMK_BMX',str(pre_msk_max)    ,header_comment='Pre-Smooth Mask Spec Region Max')
				elif pre_msk == False:
					Header_Get_Add(stck_res_file,'PMK_PMK',str(pre_msk)        ,header_comment='Pre-Smooth')

				if sig_clp == True:
					Header_Get_Add(stck_res_file,'SCP_SCP',str(sig_clp),header_comment='Sigma Clipping')
					Header_Get_Add(stck_res_file,'SCP_CUT',int(sig_cut),header_comment='Sigma Clipping low/high rejection Values')
					Header_Get_Add(stck_res_file,'SCP_CFN',str(sig_fct),header_comment='Sigma Clipping Central Function')
					Header_Get_Add(stck_res_file,'SCP_MFV',str(sig_fll),header_comment='Sigma Clipping Filling Value')
				elif sig_clp == False:
					Header_Get_Add(stck_res_file,'SCP_SCP',str(sig_clp),header_comment='Sigma Clipping')
					pass
				if  '-BS-' in name:
					pass
				elif  '-BS_MST' in name:
					Header_Get_Add(stck_res_file,'BST_NMB',str(bs_nmb_itr)  ,header_comment='Bootstrap Repetitions')
				else:
					pass
			if pst_cnt == True:
				[Header_Get_Add(stck_res_file_cnt,'PCT_PCT',str(pst_cnt)      ,header_comment='Pst-Continuum Fit')                  for stck_res_file_cnt in OPT_PCT_CNT]
				[Header_Get_Add(stck_res_file_cnt,'PCT_TYP',str(pst_cnt_typ)  ,header_comment='Continuum IRAF type')                for stck_res_file_cnt in OPT_PCT_CNT]
				[Header_Get_Add(stck_res_file_cnt,'PCT_FIT',str(pst_cnt_fnc)  ,header_comment='Continuum IRAF function')            for stck_res_file_cnt in OPT_PCT_CNT]
				[Header_Get_Add(stck_res_file_cnt,'PCT_ORD',str(pst_cnt_ovr)  ,header_comment='Continuum IRAF order')               for stck_res_file_cnt in OPT_PCT_CNT]
				[Header_Get_Add(stck_res_file_cnt,'PCT_RPC',str(pst_cnt_rpl)  ,header_comment='Continuum IRAF replacement')         for stck_res_file_cnt in OPT_PCT_CNT]
				[Header_Get_Add(stck_res_file_cnt,'PCT_LRJ',str(pst_cnt_lrj)  ,header_comment='Continuum IRAF low sigma rejection') for stck_res_file_cnt in OPT_PCT_CNT]
				[Header_Get_Add(stck_res_file_cnt,'PCT_HRJ',str(pst_cnt_hrj)  ,header_comment='Continuum IRAF high sigma rejection')for stck_res_file_cnt in OPT_PCT_CNT]
			elif pst_cnt == False:
				pass
			if smt_spc_pst == True:
				[Header_Get_Add(stck_res_file_smth,'PST_SMT',str(smt_spc_pst) ,header_comment='Post-Smooth')                        for stck_res_file_smth in OPT_PCT_SMT]
				[Header_Get_Add(stck_res_file_smth,'PSS_KRN',str(smt_shp_pst) ,header_comment='Post-Smooth kernel type')            for stck_res_file_smth in OPT_PCT_SMT]
				[Header_Get_Add(stck_res_file_smth,'PSS_SZE',int(smt_sze_pst) ,header_comment='Post-Smooth kernel size')            for stck_res_file_smth in OPT_PCT_SMT]
			elif smt_spc_pst == False:
				pass
			else:
				pass
		else:
			pass
		print
		print colored('Stats on Stacked Spectra (Core)','yellow')
		print
		[Spectra_Stats(stck_res_file_cnt) for stck_res_file_cnt in OPT_STCK_CRE]
		[Header_Get_Add(stck_res_file_wgt,'WGT_TYP',str(wgt_typ)      ,header_comment='Weight Type')          for stck_res_file_wgt in OPT_STCK_CRE]
		[Header_Get_Add(stck_res_file_wgt,'SCP_SCP',str(get_cont_flux),header_comment='Weight Cont Flux')     for stck_res_file_wgt in OPT_STCK_CRE]
		[Header_Get_Add(stck_res_file_wgt,'SCP_SCP',str(gcv_lmbd_i)   ,header_comment='Weight Cont Flux min') for stck_res_file_wgt in OPT_STCK_CRE]
		[Header_Get_Add(stck_res_file_wgt,'SCP_SCP',str(gcv_lmbd_f)   ,header_comment='Weight Cont Flux max') for stck_res_file_wgt in OPT_STCK_CRE]
		if pst_cnt == True:
			print colored('Stats on Stacked Spectra (Continuum-Fit)','yellow')
			[Spectra_Stats(stck_res_file_cnt) for stck_res_file_cnt in OPT_STCK_CNT]
			[Header_Get_Add(stck_res_file_cnt,'PCT_PCT',str(pst_cnt)      ,header_comment='Pst-Continuum Fit')                  for stck_res_file_cnt in OPT_STCK_CNT]
			[Header_Get_Add(stck_res_file_cnt,'PCT_TYP',str(pst_cnt_typ)  ,header_comment='Continuum IRAF type')                for stck_res_file_cnt in OPT_STCK_CNT]
			[Header_Get_Add(stck_res_file_cnt,'PCT_FIT',str(pst_cnt_fnc)  ,header_comment='Continuum IRAF function')            for stck_res_file_cnt in OPT_STCK_CNT]
			[Header_Get_Add(stck_res_file_cnt,'PCT_ORD',str(pst_cnt_ovr)  ,header_comment='Continuum IRAF order')               for stck_res_file_cnt in OPT_STCK_CNT]
			[Header_Get_Add(stck_res_file_cnt,'PCT_RPC',str(pst_cnt_rpl)  ,header_comment='Continuum IRAF replacement')         for stck_res_file_cnt in OPT_STCK_CNT]
			[Header_Get_Add(stck_res_file_cnt,'PCT_LRJ',str(pst_cnt_lrj)  ,header_comment='Continuum IRAF low sigma rejection') for stck_res_file_cnt in OPT_STCK_CNT]
			[Header_Get_Add(stck_res_file_cnt,'PCT_HRJ',str(pst_cnt_hrj)  ,header_comment='Continuum IRAF high sigma rejection')for stck_res_file_cnt in OPT_STCK_CNT]
		elif pst_cnt == False:
			pass
		if smt_spc_pst == True:
			print colored('Stats on Stacked Spectra (Smoothed)','yellow')
			[Spectra_Stats(stck_res_file_cnt) for stck_res_file_cnt in OPT_STCK_SMT]
			[Header_Get_Add(stck_res_file_smth,'PST_SMT',str(smt_spc_pst) ,header_comment='Post-Smooth')                        for stck_res_file_smth in OPT_STCK_SMT]
			[Header_Get_Add(stck_res_file_smth,'PSS_KRN',str(smt_shp_pst) ,header_comment='Post-Smooth kernel type')            for stck_res_file_smth in OPT_STCK_SMT]
			[Header_Get_Add(stck_res_file_smth,'PSS_SZE',int(smt_sze_pst) ,header_comment='Post-Smooth kernel size')            for stck_res_file_smth in OPT_STCK_SMT]
		elif smt_spc_pst == False:
			pass

		print 'Imaged Stacked files names: '
		print colored("\n".join([result_stacked for result_stacked in OPT_STCK_CRE]),'cyan')
		print

		if pst_cnt==False and smt_spc_pst == False:
			pass
		elif pst_cnt==True and smt_spc_pst == False:
			pass
			print colored("\n".join([result_stacked for result_stacked in OPT_STCK_CNT]),'white')
		elif pst_cnt==False and smt_spc_pst == True:
			pass
			print colored("\n".join([result_stacked for result_stacked in OPT_STCK_SMT]),'green')
		elif pst_cnt==True and smt_spc_pst == True:
			pass
			print colored("\n".join([result_stacked for result_stacked in OPT_STCK_CNT]),'white')
			print colored("\n".join([result_stacked for result_stacked in OPT_STCK_SMT]),'green')
		if stk_wgt_mde == True:
			print colored(wght_file_val,'yellow')
			print colored(wght_file_vln,'yellow')
			print colored(wght_file_hst,'yellow')
			print colored("\n".join([result_stacked for result_stacked in OPT_WGHT_FLS]),'yellow')		
		else:
			pass
		if stk_pct_mde == True and pst_cnt==False and smt_spc_pst == False:
			pass
			print colored("\n".join([result_stacked for result_stacked in OPT_STCK_PCT]),'cyan')		
		elif stk_pct_mde == True and pst_cnt==True and smt_spc_pst == False:
			print colored("\n".join([result_stacked for result_stacked in OPT_PCT_CNT]),'white')
		elif stk_pct_mde == True and pst_cnt==False and smt_spc_pst == True:
			print colored("\n".join([result_stacked for result_stacked in OPT_PCT_SMT]),'green')			
		elif stk_pct_mde == True and pst_cnt==True and smt_spc_pst == True:
			print colored("\n".join([result_stacked for result_stacked in OPT_PCT_CNT]),'white')
			print colored("\n".join([result_stacked for result_stacked in OPT_PCT_SMT]),'green')

		else:
			pass
	else:
		FNL_SPEC_RES =[]
	return FNL_SPEC_RES

def Stack(Stck_cat,name,stamps_fni,stamps_noise_fni,*args, **kwargs):
	spc_type_2bstck  = kwargs.get('spc_type_2bstck',None)
	bs_func          = kwargs.get('bs_func'        ,'')
	upd_header_stck  = kwargs.get('upd_header_stck',False)
	new_CRVAL1_head  = kwargs.get('new_CRVAL1_head',None)
	new_CDELT1_head  = kwargs.get('new_CDELT1_head',None)
	wrt_fits         = kwargs.get('wrt_fits'       ,True)

	sel_pre_shf     = kwargs.get('sel_pre_shf'     ,True)
	sel_pre_cnt     = kwargs.get('sel_pre_cnt'     ,True)
	sel_pre_msk     = kwargs.get('sel_pre_msk'     ,False)

	pre_cnt         = kwargs.get('pre_cnt'         ,False)
	pre_cnt_typ     = kwargs.get('pre_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pre_cnt_lns     = kwargs.get('pre_cnt_lns'     ,'*')       # Image lines to be fit
	pre_cnt_fnc     = kwargs.get('pre_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pre_cnt_ord     = kwargs.get('pre_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pre_cnt_ovr     = kwargs.get('pre_cnt_ovr'     ,'yes')     # Override previous norm spec
	pre_cnt_rpl     = kwargs.get('pre_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pre_cnt_lrj     = kwargs.get('pre_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pre_cnt_hrj     = kwargs.get('pre_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pre     = kwargs.get('smt_spc_pre'     ,True)
	smt_shp_pre     = kwargs.get('smt_shp_pre'     ,'gaussian')
	smt_sze_pre     = kwargs.get('smt_sze_pre'     ,1)

	pre_msk         = kwargs.get('pre_msk'         ,False)
	pre_msk_typ     = kwargs.get('pre_msk_typ'     ,'NaN')
	pre_msk_cte_val = kwargs.get('pre_msk_cte_val' ,1)
	pre_msk_abs_lne = kwargs.get('pre_msk_abs_lne' ,False)
	pre_msk_rgn     = kwargs.get('pre_msk_rgn'     ,False)
	pre_msk_min     = kwargs.get('pre_msk_min'     ,500)
	pre_msk_max     = kwargs.get('pre_msk_max'     ,1210)

	sig_clp         = kwargs.get('sig_clp'         ,False)
	sig_cut         = kwargs.get('sig_cut'         ,3)
	sig_fct         = kwargs.get('sig_fct'         ,mean)
	sig_fll         = kwargs.get('sig_fll'         ,np.nan)

	wgt_typ         = kwargs.get('wgt_typ'         ,None)
	get_cont_flux   = kwargs.get('get_cont_flux'   ,True)
	gcv_lmbd_i      = kwargs.get('gcv_lmbd_i'      ,1430)
	gcv_lmbd_f      = kwargs.get('gcv_lmbd_f'      ,1480)

	spc_nse         = kwargs.get('spc_nse'         ,False)

	pst_cnt         = kwargs.get('pst_cnt'         ,True)
	pst_cnt_typ     = kwargs.get('pst_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pst_cnt_lns     = kwargs.get('pst_cnt_lns'     ,'*')       # Image lines to be fit
	pst_cnt_fnc     = kwargs.get('pst_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pst_cnt_ord     = kwargs.get('pst_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pst_cnt_ovr     = kwargs.get('pst_cnt_ovr'     ,'yes')     # Override previous norm spec
	pst_cnt_rpl     = kwargs.get('pst_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pst_cnt_lrj     = kwargs.get('pst_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pst_cnt_hrj     = kwargs.get('pst_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pst     = kwargs.get('smt_spc_pst'     ,True)
	smt_shp_pst     = kwargs.get('smt_shp_pst'     ,'gaussian')
	smt_sze_pst     = kwargs.get('smt_sze_pst'     ,1)

	test_bg         = kwargs.get('test_bg'         ,False)
	test_fg         = kwargs.get('test_fg'         ,False)

	prt_tbl_stk     = kwargs.get('prt_tbl_stk'     ,None)


	bs_nmb_itr      = kwargs.get('bs_nmb_itr'     ,999999)

	upd_header_stck = kwargs.get('upd_header_stck', True)

	stk_pct_mde     = kwargs.get('stk_pct_mde'     ,False)
	stk_wgt_mde     = kwargs.get('stk_wgt_mde'     ,False)

	print
	print 'Stacking galaxies: ',len(stamps_fni)

	if wgt_typ == 'i-band-mag':
		wgt_var     = 'MAG_I'
		wgt_nrm_hdr = 'MAG_IN'
	elif wgt_typ == 'cont-flux-sum':
		wgt_var     = 'CFX_SUM'
		wgt_nrm_hdr = 'CFX_SUN'
	elif wgt_typ == 'cont-flux-med':
		wgt_var     = 'CFX_MED'
		wgt_nrm_hdr = 'CFX_MEN'
	elif wgt_typ == 'cont-flux-avg':
		wgt_var     = 'CFX_AVG'
		wgt_nrm_hdr = 'CFX_AVN'
	elif wgt_typ == None:
		wgt_var     = 'WGT_UNT'
		wgt_nrm_hdr = 'WGT_UNN'

	img_stack       = []
	img_noise_stack = []
	img_wghts       = []

	img_stack       = [fits.getdata(img,memmap=False)       for img in stamps_fni]
	img_wghts       = [Header_Get(img,wgt_var)              for img in stamps_fni]

	if spc_nse == True:
		img_noise_stack = [fits.getdata(img_noise,memmap=False) for img_noise in stamps_noise_fni]
	elif spc_nse == False:
		img_noise_stack = []
		pass

	norm = sum(img_wghts)
	for img in stamps_fni:
		wgt_unorm = Header_Get(img,wgt_var)
		Header_Get_Add(img,wgt_nrm_hdr,wgt_unorm/norm   ,header_comment='Weight normalized')

	img_wghts       = [Header_Get(img,wgt_nrm_hdr) for img in stamps_fni]

	print
	print colored('Stacking galaxies: '+str(len(stamps_fni))+', '+str(len(img_stack)),'yellow')
	print colored('With a number of elements: '+str(len(img_stack[0])),'yellow')
	print colored('Selected weight: '+wgt_var,'yellow')
	print

	stck_stat   = Stack_Img_Op(Stck_cat,name,img_stack,spc_nse, img_wghts,
				test_bg         = test_bg           ,test_fg         = test_fg,
				stack_ext       = spc_type_2bstck   ,bs_func         = bs_func         ,new_CRVAL1_head = new_CRVAL1_head       ,new_CDELT1_head = new_CDELT1_head,
				wgt_var         = wgt_var           ,wgt_nrm_hdr     = wgt_nrm_hdr     ,wrt_fits        = wrt_fits              ,
				sel_pre_shf     = sel_pre_shf       ,sel_pre_cnt     = sel_pre_cnt     ,sel_pre_msk     = sel_pre_msk           ,
				pre_cnt         = pre_cnt           ,pre_cnt_typ     = pre_cnt_typ     ,pre_cnt_lns     = pre_cnt_lns           ,
				pre_cnt_fnc     = pre_cnt_fnc       ,pre_cnt_ord     = pre_cnt_ord     ,pre_cnt_ovr     = pre_cnt_ovr           ,
				pre_cnt_rpl     = pre_cnt_rpl       ,pre_cnt_lrj     = pre_cnt_lrj     ,pre_cnt_hrj     = pre_cnt_hrj           ,
				smt_spc_pre     = smt_spc_pre       ,smt_shp_pre     = smt_shp_pre     ,smt_sze_pre     = smt_sze_pre           ,
				pre_msk         = pre_msk           ,
				pre_msk_typ     = pre_msk_typ       ,pre_msk_cte_val = pre_msk_cte_val ,pre_msk_abs_lne = pre_msk_abs_lne,
				pre_msk_rgn     = pre_msk_rgn       ,pre_msk_min     = pre_msk_min     ,pre_msk_max     = pre_msk_max           ,										
				noise_imgs      = img_noise_stack   ,
				wgt_typ         = wgt_typ           ,
				get_cont_flux   = get_cont_flux     ,gcv_lmbd_i      = gcv_lmbd_i      ,gcv_lmbd_f      = gcv_lmbd_f,
				sig_clp         = sig_clp           ,
				sig_cut         = sig_cut           ,sig_fct         = sig_fct         ,sig_fll         = sig_fll,
				pst_cnt         = pst_cnt           ,pst_cnt_typ     = pst_cnt_typ     ,pst_cnt_lns     = pst_cnt_lns           ,
				pst_cnt_fnc     = pst_cnt_fnc       ,pst_cnt_ord     = pst_cnt_ord     ,pst_cnt_ovr     = pst_cnt_ovr           ,
				pst_cnt_rpl     = pst_cnt_rpl       ,pst_cnt_lrj     = pst_cnt_lrj     ,pst_cnt_hrj     = pst_cnt_hrj           ,
				smt_spc_pst     = smt_spc_pst       ,smt_shp_pst     = smt_shp_pst     ,smt_sze_pst     = smt_sze_pst           ,
				bs_nmb_itr      = bs_nmb_itr,
				stk_pct_mde     = stk_pct_mde       ,stk_wgt_mde     = stk_wgt_mde)

	if upd_header_stck == True and not  'BS_MST' in prt_tbl_stk:
		print
		print 'upd_header_stck',upd_header_stck
		print
		pass

		print
		print colored('Statistics from table '+prt_tbl_stk,'yellow')
		print
		stt_tbl_stk = readtable_fg_bg_glx(prt_tbl_stk,tbl_format_ipt,bs_func=bs_func)
		z_fg    = stt_tbl_stk[2]
		z_bg    = stt_tbl_stk[5]
		sep_as  = stt_tbl_stk[8]
		sep_kpc = stt_tbl_stk[10]

		z_fg_sample_avg     = np.mean(z_fg)
		z_fg_sample_med     = np.median(z_fg)
		z_fg_sample_1sl     = np.nanpercentile(z_fg, 15.9)
		z_fg_sample_1sh     = np.nanpercentile(z_fg, 84.1)
		z_fg_sample_2sl     = np.nanpercentile(z_fg, 2.30)
		z_fg_sample_2sh     = np.nanpercentile(z_fg, 97.7)
		z_fg_sample_3sl     = np.nanpercentile(z_fg, 0.20)
		z_fg_sample_3sh     = np.nanpercentile(z_fg, 99.8)
		z_fg_sample_p25     = np.nanpercentile(z_fg, 25.0)
		z_fg_sample_p75     = np.nanpercentile(z_fg, 75.0)

		z_bg_sample_avg     = np.mean(z_bg)
		z_bg_sample_med     = np.median(z_bg)
		z_bg_sample_1sl     = np.nanpercentile(z_bg, 15.9)
		z_bg_sample_1sh     = np.nanpercentile(z_bg, 84.1)
		z_bg_sample_2sl     = np.nanpercentile(z_bg, 2.30)
		z_bg_sample_2sh     = np.nanpercentile(z_bg, 97.7)
		z_bg_sample_3sl     = np.nanpercentile(z_bg, 0.20)
		z_bg_sample_3sh     = np.nanpercentile(z_bg, 99.8)
		z_bg_sample_p25     = np.nanpercentile(z_bg, 25.0)
		z_bg_sample_p75     = np.nanpercentile(z_bg, 75.0)

		sep_as_sample_avg     = np.mean(sep_as)
		sep_as_sample_med     = np.median(sep_as)
		sep_as_sample_1sl     = np.nanpercentile(sep_as, 15.9)
		sep_as_sample_1sh     = np.nanpercentile(sep_as, 84.1)
		sep_as_sample_2sl     = np.nanpercentile(sep_as, 2.30)
		sep_as_sample_2sh     = np.nanpercentile(sep_as, 97.7)
		sep_as_sample_3sl     = np.nanpercentile(sep_as, 0.20)
		sep_as_sample_3sh     = np.nanpercentile(sep_as, 99.8)
		sep_as_sample_p25     = np.nanpercentile(sep_as, 25.0)
		sep_as_sample_p75     = np.nanpercentile(sep_as, 75.0)

		sep_kpc_sample_avg     = np.mean(sep_kpc)
		sep_kpc_sample_med     = np.median(sep_kpc)
		sep_kpc_sample_1sl     = np.nanpercentile(sep_kpc, 15.9)
		sep_kpc_sample_1sh     = np.nanpercentile(sep_kpc, 84.1)
		sep_kpc_sample_2sl     = np.nanpercentile(sep_kpc, 2.30)
		sep_kpc_sample_2sh     = np.nanpercentile(sep_kpc, 97.7)
		sep_kpc_sample_3sl     = np.nanpercentile(sep_kpc, 0.20)
		sep_kpc_sample_3sh     = np.nanpercentile(sep_kpc, 99.8)
		sep_kpc_sample_p25     = np.nanpercentile(sep_kpc, 25.0)
		sep_kpc_sample_p75     = np.nanpercentile(sep_kpc, 75.0)

		z_fg_sample_avg_hdr = 'ZFG_AVG'
		z_fg_sample_med_hdr = 'ZFG_MED'
		z_fg_sample_1sl_hdr = 'ZFG_1SL'
		z_fg_sample_1sh_hdr = 'ZFG_1SH'
		z_fg_sample_2sl_hdr = 'ZFG_2SL'
		z_fg_sample_2sh_hdr = 'ZFG_2SH'
		z_fg_sample_3sl_hdr = 'ZFG_3SL'
		z_fg_sample_3sh_hdr = 'ZFG_3SH'
		z_fg_sample_p25_hdr = 'ZFG_P25'
		z_fg_sample_p75_hdr = 'ZFG_P75'

		z_bg_sample_avg_hdr = 'ZBG_AVG'
		z_bg_sample_med_hdr = 'ZBG_MED'
		z_bg_sample_1sl_hdr = 'ZBG_1SL'
		z_bg_sample_1sh_hdr = 'ZBG_1SH'
		z_bg_sample_2sl_hdr = 'ZBG_2SL'
		z_bg_sample_2sh_hdr = 'ZBG_2SH'
		z_bg_sample_3sl_hdr = 'ZBG_3SL'
		z_bg_sample_3sh_hdr = 'ZBG_3SH'
		z_bg_sample_p25_hdr = 'ZBG_P25'
		z_bg_sample_p75_hdr = 'ZBG_P75'

		sep_as_sample_avg_hdr = 'SAS_AVG'
		sep_as_sample_med_hdr = 'SAS_MED'
		sep_as_sample_1sl_hdr = 'SAS_1SL'
		sep_as_sample_1sh_hdr = 'SAS_1SH'
		sep_as_sample_2sl_hdr = 'SAS_2SL'
		sep_as_sample_2sh_hdr = 'SAS_2SH'
		sep_as_sample_3sl_hdr = 'SAS_3SL'
		sep_as_sample_3sh_hdr = 'SAS_3SH'
		sep_as_sample_p25_hdr = 'SAS_P25'
		sep_as_sample_p75_hdr = 'SAS_P75'

		sep_kpc_sample_avg_hdr = 'SKP_AVG'
		sep_kpc_sample_med_hdr = 'SKP_MED'
		sep_kpc_sample_1sl_hdr = 'SKP_1SL'
		sep_kpc_sample_1sh_hdr = 'SKP_1SH'
		sep_kpc_sample_2sl_hdr = 'SKP_2SL'
		sep_kpc_sample_2sh_hdr = 'SKP_2SH'
		sep_kpc_sample_3sl_hdr = 'SKP_3SL'
		sep_kpc_sample_3sh_hdr = 'SKP_3SH'
		sep_kpc_sample_p25_hdr = 'SKP_P25'
		sep_kpc_sample_p75_hdr = 'SKP_P75'

		z_fg_sample_avg_cmt = 'Stat z fg avg'
		z_fg_sample_med_cmt = 'Stat z fg med'
		z_fg_sample_1sl_cmt = 'Stat z fg 1sl'
		z_fg_sample_1sh_cmt = 'Stat z fg 1sh'
		z_fg_sample_2sl_cmt = 'Stat z fg 2sl'
		z_fg_sample_2sh_cmt = 'Stat z fg 2sh'
		z_fg_sample_3sl_cmt = 'Stat z fg 3sl'
		z_fg_sample_3sh_cmt = 'Stat z fg 3sh'
		z_fg_sample_p25_cmt = 'Stat z fg p25'
		z_fg_sample_p75_cmt = 'Stat z fg p75'

		z_bg_sample_avg_cmt = 'Stat z bg avg'
		z_bg_sample_med_cmt = 'Stat z bg med'
		z_bg_sample_1sl_cmt = 'Stat z bg 1sl'
		z_bg_sample_1sh_cmt = 'Stat z bg 1sh'
		z_bg_sample_2sl_cmt = 'Stat z bg 2sl'
		z_bg_sample_2sh_cmt = 'Stat z bg 2sh'
		z_bg_sample_3sl_cmt = 'Stat z bg 3sl'
		z_bg_sample_3sh_cmt = 'Stat z bg 3sh'
		z_bg_sample_p25_cmt = 'Stat z bg p25'
		z_bg_sample_p75_cmt = 'Stat z bg p75'

		sep_as_sample_avg_cmt = 'Stat sep arcsec avg'
		sep_as_sample_med_cmt = 'Stat sep arcsec med'
		sep_as_sample_1sl_cmt = 'Stat sep arcsec 1sl'
		sep_as_sample_1sh_cmt = 'Stat sep arcsec 1sh'
		sep_as_sample_2sl_cmt = 'Stat sep arcsec 2sl'
		sep_as_sample_2sh_cmt = 'Stat sep arcsec 2sh'
		sep_as_sample_3sl_cmt = 'Stat sep arcsec 3sl'
		sep_as_sample_3sh_cmt = 'Stat sep arcsec 3sh'
		sep_as_sample_p25_cmt = 'Stat sep arcsec p25'
		sep_as_sample_p75_cmt = 'Stat sep arcsec p75'

		sep_kpc_sample_avg_cmt = 'Stat sep kpc avg'
		sep_kpc_sample_med_cmt = 'Stat sep kpc med'
		sep_kpc_sample_1sl_cmt = 'Stat sep kpc 1sl'
		sep_kpc_sample_1sh_cmt = 'Stat sep kpc 1sh'
		sep_kpc_sample_2sl_cmt = 'Stat sep kpc 2sl'
		sep_kpc_sample_2sh_cmt = 'Stat sep kpc 2sh'
		sep_kpc_sample_3sl_cmt = 'Stat sep kpc 3sl'
		sep_kpc_sample_3sh_cmt = 'Stat sep kpc 3sh'
		sep_kpc_sample_p25_cmt = 'Stat sep kpc p25'
		sep_kpc_sample_p75_cmt = 'Stat sep kpc p75'

		Z_FG_HDR = [
					z_fg_sample_avg_hdr,
					z_fg_sample_med_hdr,
					z_fg_sample_1sl_hdr,
					z_fg_sample_1sh_hdr,
					z_fg_sample_2sl_hdr,
					z_fg_sample_2sh_hdr,
					z_fg_sample_3sl_hdr,
					z_fg_sample_3sh_hdr,
					z_fg_sample_p25_hdr,
					z_fg_sample_p75_hdr
				]
		Z_BG_HDR = [
					z_bg_sample_avg_hdr,
					z_bg_sample_med_hdr,
					z_bg_sample_1sl_hdr,
					z_bg_sample_1sh_hdr,
					z_bg_sample_2sl_hdr,
					z_bg_sample_2sh_hdr,
					z_bg_sample_3sl_hdr,
					z_bg_sample_3sh_hdr,
					z_bg_sample_p25_hdr,
					z_bg_sample_p75_hdr
				]				
		SEP_AS_HDR = [
					sep_as_sample_avg_hdr,
					sep_as_sample_med_hdr,
					sep_as_sample_1sl_hdr,
					sep_as_sample_1sh_hdr,
					sep_as_sample_2sl_hdr,
					sep_as_sample_2sh_hdr,
					sep_as_sample_3sl_hdr,
					sep_as_sample_3sh_hdr,
					sep_as_sample_p25_hdr,
					sep_as_sample_p75_hdr
				]
		SEP_KPC_HDR = [
					sep_kpc_sample_avg_hdr,
					sep_kpc_sample_med_hdr,
					sep_kpc_sample_1sl_hdr,
					sep_kpc_sample_1sh_hdr,
					sep_kpc_sample_2sl_hdr,
					sep_kpc_sample_2sh_hdr,
					sep_kpc_sample_3sl_hdr,
					sep_kpc_sample_3sh_hdr,
					sep_kpc_sample_p25_hdr,
					sep_kpc_sample_p75_hdr
				]

		Z_FG_CMT = [
					z_fg_sample_avg_cmt,
					z_fg_sample_med_cmt,
					z_fg_sample_1sl_cmt,
					z_fg_sample_1sh_cmt,
					z_fg_sample_2sl_cmt,
					z_fg_sample_2sh_cmt,
					z_fg_sample_3sl_cmt,
					z_fg_sample_3sh_cmt,
					z_fg_sample_p25_cmt,
					z_fg_sample_p75_cmt
				]
		Z_BG_CMT = [
					z_bg_sample_avg_cmt,
					z_bg_sample_med_cmt,
					z_bg_sample_1sl_cmt,
					z_bg_sample_1sh_cmt,
					z_bg_sample_2sl_cmt,
					z_bg_sample_2sh_cmt,
					z_bg_sample_3sl_cmt,
					z_bg_sample_3sh_cmt,
					z_bg_sample_p25_cmt,
					z_bg_sample_p75_cmt
				]
		SEP_AS_CMT = [				
					sep_as_sample_avg_cmt,
					sep_as_sample_med_cmt,
					sep_as_sample_1sl_cmt,
					sep_as_sample_1sh_cmt,
					sep_as_sample_2sl_cmt,
					sep_as_sample_2sh_cmt,
					sep_as_sample_3sl_cmt,
					sep_as_sample_3sh_cmt,
					sep_as_sample_p25_cmt,
					sep_as_sample_p75_cmt
				]			
		SEP_KPC_CMT = [				
					sep_kpc_sample_avg_cmt,
					sep_kpc_sample_med_cmt,
					sep_kpc_sample_1sl_cmt,
					sep_kpc_sample_1sh_cmt,
					sep_kpc_sample_2sl_cmt,
					sep_kpc_sample_2sh_cmt,
					sep_kpc_sample_3sl_cmt,
					sep_kpc_sample_3sh_cmt,
					sep_kpc_sample_p25_cmt,
					sep_kpc_sample_p75_cmt
				]
		Z_FG_VAL = [				
					z_fg_sample_avg,
					z_fg_sample_med,
					z_fg_sample_1sl,
					z_fg_sample_1sh,
					z_fg_sample_2sl,
					z_fg_sample_2sh,
					z_fg_sample_3sl,
					z_fg_sample_3sh,
					z_fg_sample_p25,
					z_fg_sample_p75
				]
		Z_BG_VAL = [				
					z_bg_sample_avg,
					z_bg_sample_med,
					z_bg_sample_1sl,
					z_bg_sample_1sh,
					z_bg_sample_2sl,
					z_bg_sample_2sh,
					z_bg_sample_3sl,
					z_bg_sample_3sh,
					z_bg_sample_p25,
					z_bg_sample_p75
				]
		SEP_AS_VAL = [				
					sep_as_sample_avg,
					sep_as_sample_med,
					sep_as_sample_1sl,
					sep_as_sample_1sh,
					sep_as_sample_2sl,
					sep_as_sample_2sh,
					sep_as_sample_3sl,
					sep_as_sample_3sh,
					sep_as_sample_p25,
					sep_as_sample_p75
				]
		SEP_KPC_VAL = [				
					sep_kpc_sample_avg,
					sep_kpc_sample_med,
					sep_kpc_sample_1sl,
					sep_kpc_sample_1sh,
					sep_kpc_sample_2sl,
					sep_kpc_sample_2sh,
					sep_kpc_sample_3sl,
					sep_kpc_sample_3sh,
					sep_kpc_sample_p25,
					sep_kpc_sample_p75
				]
		STK_SPC_RES = [stamps_fni,stamps_noise_fni]
		print
		print  "\n".join([(str(element)) for element in stck_stat])
		print
		print colored('Adding redshift, sep (as & kpc) stats (bg and fg galaxies) on stacked files.','yellow')
		print
		[Header_Get_Add(element[0],element[1][0],element[1][1],header_comment=element[1][2]) for element in itpdc(stck_stat,zip(Z_FG_HDR,Z_FG_VAL,Z_FG_CMT))]
		[Header_Get_Add(element[0],element[1][0],element[1][1],header_comment=element[1][2]) for element in itpdc(stck_stat,zip(Z_BG_HDR,Z_BG_VAL,Z_BG_CMT))]
		[Header_Get_Add(element[0],element[1][0],element[1][1],header_comment=element[1][2]) for element in itpdc(stck_stat,zip(SEP_AS_HDR,SEP_AS_VAL,SEP_AS_CMT))]
		[Header_Get_Add(element[0],element[1][0],element[1][1],header_comment=element[1][2]) for element in itpdc(stck_stat,zip(SEP_KPC_HDR,SEP_KPC_VAL,SEP_KPC_CMT))]

		if (('-PRP-' in prt_tbl_stk) and 'mass_F'in prt_tbl_stk) or (('_PRP_MRP' in prt_tbl_stk) and 'mass_F'in prt_tbl_stk):
			var_prp_slc = stt_tbl_stk[20]
			var_prp_hdr = 'MSF'
			var_prp_cmt = 'Mass Fg Glxs'
		elif (('-PRP-' in prt_tbl_stk) and 'Age_F'in prt_tbl_stk) or (('_PRP_MRP' in prt_tbl_stk) and 'Age_F'in prt_tbl_stk):
			var_prp_slc = stt_tbl_stk[21]
			var_prp_hdr = 'AGF'		
			var_prp_cmt = 'Age Fg Glxs'
		elif (('-PRP-' in prt_tbl_stk) and 'SFR_F'in prt_tbl_stk) or (('_PRP_MRP' in prt_tbl_stk) and 'SFR_F'in prt_tbl_stk):
			var_prp_slc = stt_tbl_stk[22]
			var_prp_hdr = 'SFR'		
			var_prp_cmt = 'SFR Fg Glxs'
		elif (('-PRP-' in prt_tbl_stk) and '-sSFR_F-'in prt_tbl_stk) or (('_PRP_MRP'in prt_tbl_stk) and '-sSFR_F-'in prt_tbl_stk):
			var_prp_slc = stt_tbl_stk[23]
			var_prp_hdr = 'sSF'
			var_prp_cmt = 'sSF Fg Glxs'
		elif (('-PRP-' in prt_tbl_stk) and 'Lnuv_F'in prt_tbl_stk) or (('_PRP_MRP' in prt_tbl_stk) and 'Lnuv_F'in prt_tbl_stk):
			var_prp_slc = stt_tbl_stk[24]
			var_prp_hdr = 'UVF'		
			var_prp_cmt = 'LNUV Fg Glxs'
		elif (('-PRP-' in prt_tbl_stk) and 'magi_F'in prt_tbl_stk) or (('_PRP_MRP' in prt_tbl_stk) and 'magi_F'in prt_tbl_stk):
			var_prp_slc = stt_tbl_stk[25]
			var_prp_hdr = 'magi'		
			var_prp_cmt = 'magi Fg Glxs'
		elif ('_PRP_MRP' in prt_tbl_stk) and 'PHI'in prt_tbl_stk:
			print
			print len(stt_tbl_stk)
			var_prp_slc = stt_tbl_stk[26]
			var_prp_hdr = 'PHI'		
			var_prp_cmt = 'Phi angle between Fg sma and Bg'
		elif ('_PRP_MRP' in prt_tbl_stk) and 're_F'in prt_tbl_stk:
			print
			print len(stt_tbl_stk)
			var_prp_slc = stt_tbl_stk[31]
			var_prp_hdr = 'ref'
			var_prp_cmt = 'Effective Radius Fg Glxs'
		elif ('_PRP_MRP' in prt_tbl_stk) and 'n_F'in prt_tbl_stk:
			print
			print len(stt_tbl_stk)
			var_prp_slc = stt_tbl_stk[32]
			var_prp_hdr = 'nsr'		
			var_prp_cmt = 'Sersic coefficient Fg Glxs'
		elif ('_PRP_MRP' in prt_tbl_stk) and 'q_F'in prt_tbl_stk:
			print
			print len(stt_tbl_stk)
			var_prp_slc = stt_tbl_stk[33]
			var_prp_hdr = 'q'		
			var_prp_cmt = 'Inclination q -> i Fg Glxs'
		else:
			print
			print len(stt_tbl_stk)
			var_prp_slc = stt_tbl_stk[10]
			var_prp_hdr = 'XXX'
			var_prp_cmt = 'SEP [kpc] Fg Glxs XXX ELSE CASE'

		var_prp_slc_fg_sample_avg     = np.mean(var_prp_slc)
		var_prp_slc_fg_sample_med     = np.median(var_prp_slc)
		var_prp_slc_fg_sample_1sl     = np.nanpercentile(var_prp_slc, 15.9)
		var_prp_slc_fg_sample_1sh     = np.nanpercentile(var_prp_slc, 84.1)
		var_prp_slc_fg_sample_2sl     = np.nanpercentile(var_prp_slc, 2.30)
		var_prp_slc_fg_sample_2sh     = np.nanpercentile(var_prp_slc, 97.7)
		var_prp_slc_fg_sample_3sl     = np.nanpercentile(var_prp_slc, 0.20)
		var_prp_slc_fg_sample_3sh     = np.nanpercentile(var_prp_slc, 99.8)
		var_prp_slc_fg_sample_p25     = np.nanpercentile(var_prp_slc, 25.0)
		var_prp_slc_fg_sample_p75     = np.nanpercentile(var_prp_slc, 75.0)

		var_prp_slc_fg_sample_avg_hdr = var_prp_hdr + '_AVG'
		var_prp_slc_fg_sample_med_hdr = var_prp_hdr + '_MED'
		var_prp_slc_fg_sample_1sl_hdr = var_prp_hdr + '_1SL'
		var_prp_slc_fg_sample_1sh_hdr = var_prp_hdr + '_1SH'
		var_prp_slc_fg_sample_2sl_hdr = var_prp_hdr + '_2SL'
		var_prp_slc_fg_sample_2sh_hdr = var_prp_hdr + '_2SH'
		var_prp_slc_fg_sample_3sl_hdr = var_prp_hdr + '_3SL'
		var_prp_slc_fg_sample_3sh_hdr = var_prp_hdr + '_3SH'
		var_prp_slc_fg_sample_p25_hdr = var_prp_hdr + '_P25'
		var_prp_slc_fg_sample_p75_hdr = var_prp_hdr + '_P75'

		var_prp_slc_fg_sample_avg_cmt = var_prp_cmt + ' avg'
		var_prp_slc_fg_sample_med_cmt = var_prp_cmt + ' med'
		var_prp_slc_fg_sample_1sl_cmt = var_prp_cmt + ' 1sl'
		var_prp_slc_fg_sample_1sh_cmt = var_prp_cmt + ' 1sh'
		var_prp_slc_fg_sample_2sl_cmt = var_prp_cmt + ' 2sl'
		var_prp_slc_fg_sample_2sh_cmt = var_prp_cmt + ' 2sh'
		var_prp_slc_fg_sample_3sl_cmt = var_prp_cmt + ' 3sl'
		var_prp_slc_fg_sample_3sh_cmt = var_prp_cmt + ' 3sh'
		var_prp_slc_fg_sample_p25_cmt = var_prp_cmt + ' p25'
		var_prp_slc_fg_sample_p75_cmt = var_prp_cmt + ' p75'

		VAR_PRP_SLC_FG_HDR = [
					var_prp_slc_fg_sample_avg_hdr,
					var_prp_slc_fg_sample_med_hdr,
					var_prp_slc_fg_sample_1sl_hdr,
					var_prp_slc_fg_sample_1sh_hdr,
					var_prp_slc_fg_sample_2sl_hdr,
					var_prp_slc_fg_sample_2sh_hdr,
					var_prp_slc_fg_sample_3sl_hdr,
					var_prp_slc_fg_sample_3sh_hdr,
					var_prp_slc_fg_sample_p25_hdr,
					var_prp_slc_fg_sample_p75_hdr
				]
		VAR_PRP_SLC_FG_CMT = [
					var_prp_slc_fg_sample_avg_cmt,
					var_prp_slc_fg_sample_med_cmt,
					var_prp_slc_fg_sample_1sl_cmt,
					var_prp_slc_fg_sample_1sh_cmt,
					var_prp_slc_fg_sample_2sl_cmt,
					var_prp_slc_fg_sample_2sh_cmt,
					var_prp_slc_fg_sample_3sl_cmt,
					var_prp_slc_fg_sample_3sh_cmt,
					var_prp_slc_fg_sample_p25_cmt,
					var_prp_slc_fg_sample_p75_cmt
				]
		VAR_PRP_SLC_FG_VAL = [				
					var_prp_slc_fg_sample_avg,
					var_prp_slc_fg_sample_med,
					var_prp_slc_fg_sample_1sl,
					var_prp_slc_fg_sample_1sh,
					var_prp_slc_fg_sample_2sl,
					var_prp_slc_fg_sample_2sh,
					var_prp_slc_fg_sample_3sl,
					var_prp_slc_fg_sample_3sh,
					var_prp_slc_fg_sample_p25,
					var_prp_slc_fg_sample_p75
				]
		print
		print colored('Adding Var ('+var_prp_cmt+') stats (bg and fg galaxies) on stacked files.','yellow')
		print colored('Adding Var ('+var_prp_cmt+') header type: ' + VAR_PRP_SLC_FG_HDR[1],'yellow')
		print 
		print 'VAR_PRP_SLC_FG_HDR','VAR_PRP_SLC_FG_VAL','VAR_PRP_SLC_FG_CMT'
		print len(VAR_PRP_SLC_FG_HDR),len(VAR_PRP_SLC_FG_VAL),len(VAR_PRP_SLC_FG_CMT)
		print 'stck_stat'
		print len(stck_stat)
		[Header_Get_Add(element[0],element[1][0],element[1][1],header_comment=element[1][2]) for element in itpdc(stck_stat,zip(VAR_PRP_SLC_FG_HDR,VAR_PRP_SLC_FG_VAL,VAR_PRP_SLC_FG_CMT))]
	elif 'BS_MST' in prt_tbl_stk:
		print
		print colored('BS_MST table!','yellow')
		print colored('The headers of the stacked fits files will NOT be updated with the corresponding variable stats','yellow')
		print colored('upd_header_stck '+ str(upd_header_stck),'yellow')
		print
		STK_SPC_RES = [stamps_fni,stamps_noise_fni]
	else:
		pass
		print
		print colored('The headers of the stacked fits files will NOT be updated with the corresponding variable stats','yellow')
		print colored('upd_header_stck '+ str(upd_header_stck),'yellow')
		print
		STK_SPC_RES = [stamps_fni,stamps_noise_fni]
	return STK_SPC_RES,stck_stat

def Stack_Subsample(SubSmpl,*args, **kwargs):
	bs_func          = kwargs.get('bs_func'        ,'')
	wrt_fits         = kwargs.get('wrt_fits'       ,True)

	sel_pre_shf     = kwargs.get('sel_pre_shf'     ,True)
	sel_pre_cnt     = kwargs.get('sel_pre_cnt'     ,True)
	sel_pre_msk     = kwargs.get('sel_pre_msk'     ,False)

	pre_cnt         = kwargs.get('pre_cnt'         ,False)
	pre_cnt_typ     = kwargs.get('pre_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pre_cnt_lns     = kwargs.get('pre_cnt_lns'     ,'*')       # Image lines to be fit
	pre_cnt_fnc     = kwargs.get('pre_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pre_cnt_ord     = kwargs.get('pre_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pre_cnt_ovr     = kwargs.get('pre_cnt_ovr'     ,'yes')     # Override previous norm spec
	pre_cnt_rpl     = kwargs.get('pre_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pre_cnt_lrj     = kwargs.get('pre_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pre_cnt_hrj     = kwargs.get('pre_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pre     = kwargs.get('smt_spc_pre'     ,True)
	smt_shp_pre     = kwargs.get('smt_shp_pre'     ,'gaussian')
	smt_sze_pre     = kwargs.get('smt_sze_pre'     ,1)

	pre_msk         = kwargs.get('pre_msk'         ,False)
	pre_msk_typ     = kwargs.get('pre_msk_typ'     ,'NaN')
	pre_msk_cte_val = kwargs.get('pre_msk_cte_val' ,1)
	pre_msk_abs_lne = kwargs.get('pre_msk_abs_lne' ,False)
	pre_msk_rgn     = kwargs.get('pre_msk_rgn' ,False)
	pre_msk_min     = kwargs.get('pre_msk_min' ,500)
	pre_msk_max     = kwargs.get('pre_msk_max' ,1210)

	sig_clp         = kwargs.get('sig_clp'         ,False)
	sig_cut         = kwargs.get('sig_cut'         ,3)
	sig_fct         = kwargs.get('sig_fct'         ,mean)
	sig_fll         = kwargs.get('sig_fll'         ,np.nan)

	wgt_typ         = kwargs.get('wgt_typ'         ,None)
	get_cont_flux   = kwargs.get('get_cont_flux'   ,True)
	gcv_lmbd_i      = kwargs.get('gcv_lmbd_i'      ,1430)
	gcv_lmbd_f      = kwargs.get('gcv_lmbd_f'      ,1480)

	spc_nse         = kwargs.get('spc_nse'         ,False)

	pst_cnt         = kwargs.get('pst_cnt'         ,True)
	pst_cnt_typ     = kwargs.get('pst_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pst_cnt_lns     = kwargs.get('pst_cnt_lns'     ,'*')       # Image lines to be fit
	pst_cnt_fnc     = kwargs.get('pst_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pst_cnt_ord     = kwargs.get('pst_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pst_cnt_ovr     = kwargs.get('pst_cnt_ovr'     ,'yes')     # Override previous norm spec
	pst_cnt_rpl     = kwargs.get('pst_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pst_cnt_lrj     = kwargs.get('pst_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pst_cnt_hrj     = kwargs.get('pst_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pst     = kwargs.get('smt_spc_pst'     ,True)
	smt_shp_pst     = kwargs.get('smt_shp_pst'     ,'gaussian')
	smt_sze_pst     = kwargs.get('smt_sze_pst'     ,1)

	test_bg         = kwargs.get('test_bg'         ,False)
	test_fg         = kwargs.get('test_fg'         ,False)

	bs_nmb_itr      = kwargs.get('bs_nmb_itr'      ,999999)

	upd_header_stck = kwargs.get('upd_header_stck' ,True)

	stk_pct_mde     = kwargs.get('stk_pct_mde'     ,False)
	stk_wgt_mde     = kwargs.get('stk_wgt_mde'     ,False)

	print
	print colored('Stack Percentile (PCT) mode: ' + str(stk_pct_mde),'yellow')
	print colored('Stack Weight     (WGT) mode: ' + str(stk_wgt_mde),'yellow')
	print colored('BS repetition number       : ' + str(bs_nmb_itr),'yellow')
	print

	Stack_Spec_Res = []

	for SubSet in range(len(SubSmpl[0])):
		print
		print 'Reading files from table: '
		print colored(str(SubSmpl[0][SubSet]),'cyan')
		print 'Using files as: '
		print colored(str(SubSmpl[1][0][0]),'yellow')
		print colored(str(SubSmpl[1][0][-1]),'yellow')

		SubSmpl_nm = str((SubSmpl[0][SubSet].split(tbl_ext_opt,1)[0]).split('/')[-1])
		print 

		BIG_SPECTRA  = new_wave_range(min(SubSmpl[3][SubSet]),max(SubSmpl[5][SubSet]),max(SubSmpl[4][SubSet]))
		print
		print 'The wavelength range for this subset: ',min(SubSmpl[3][SubSet]),'-',max(SubSmpl[5][SubSet])
		print 'With a wavelength binsize           : ',max(SubSmpl[4][SubSet])
		print
		print colored('Interpolating ' + str(len(SubSmpl[1][SubSet])) + ' spectra','cyan')
		print colored('Previous identical interpolated file(s) will be deleted (if any)!','yellow')

		if spc_nse == True:
			spec_n = [Interpolating_Spectra(str(SubSmpl[0][SubSet]),spc2bitp,max(SubSmpl[4][SubSet]),BIG_SPECTRA[0],
			SS_indx=SubSet,
			sel_pre_shf     = sel_pre_shf     , sel_pre_cnt     = sel_pre_cnt    , sel_pre_msk     = sel_pre_msk    ,
			pre_cnt         = pre_cnt         , pre_cnt_typ     = pre_cnt_typ    , pre_cnt_lns     = pre_cnt_lns    ,    
			pre_cnt_fnc     = pre_cnt_fnc     , pre_cnt_ord     = pre_cnt_ord    , pre_cnt_ovr     = pre_cnt_ovr    ,    
			pre_cnt_rpl     = pre_cnt_rpl     , pre_cnt_lrj     = pre_cnt_lrj    , pre_cnt_hrj     = pre_cnt_hrj    ,    
			smt_spc_pre     = smt_spc_pre     , smt_shp_pre     = smt_shp_pre    , smt_sze_pre     = smt_sze_pre    ,
			pre_msk         = pre_msk         , pre_msk_typ     = pre_msk_typ    , pre_msk_abs_lne = pre_msk_abs_lne,
			pre_msk_rgn     = pre_msk_rgn     , pre_msk_min     = pre_msk_min    , pre_msk_max     = pre_msk_max    ,
			get_cont_flux   = False           , gcv_lmbd_i      = gcv_lmbd_i     , gcv_lmbd_f      = gcv_lmbd_f) for spc2bitp in SubSmpl[2][SubSet]] 
		elif spc_nse == False:
			spec_n = [0]
			pass

		spec   = [Interpolating_Spectra(str(SubSmpl[0][SubSet]),spc2bitp,max(SubSmpl[4][SubSet]),BIG_SPECTRA[0],
			SS_indx         = SubSet,
			sel_pre_shf     = sel_pre_shf     , sel_pre_cnt     = sel_pre_cnt    , sel_pre_msk     = sel_pre_msk    ,
			pre_cnt         = pre_cnt         , pre_cnt_typ     = pre_cnt_typ    , pre_cnt_lns     = pre_cnt_lns    ,    
			pre_cnt_fnc     = pre_cnt_fnc     , pre_cnt_ord     = pre_cnt_ord    , pre_cnt_ovr     = pre_cnt_ovr    ,    
			pre_cnt_rpl     = pre_cnt_rpl     , pre_cnt_lrj     = pre_cnt_lrj    , pre_cnt_hrj     = pre_cnt_hrj    ,    
			smt_spc_pre     = smt_spc_pre     , smt_shp_pre     = smt_shp_pre    , smt_sze_pre     = smt_sze_pre    ,
			pre_msk         = pre_msk         , pre_msk_typ     = pre_msk_typ    , pre_msk_abs_lne = pre_msk_abs_lne,
			pre_msk_rgn     = pre_msk_rgn     , pre_msk_min     = pre_msk_min    , pre_msk_max     = pre_msk_max    ,
			get_cont_flux   = get_cont_flux   , gcv_lmbd_i      = gcv_lmbd_i     , gcv_lmbd_f      = gcv_lmbd_f) for spc2bitp in SubSmpl[1][SubSet]]
		if '-BS-' in SubSmpl_nm:
			cpy_fts_dir = fts_bst_lst
		elif  '-BS_MST' in SubSmpl_nm:
			cpy_fts_dir = ind_bst_lst
		else:
			cpy_fts_dir = ind_stk_res

		if len(spec)==1:
			spec.append(spec[0])
			spec_n.append(spec_n[0])
		elif len(spec)>1:
			pass

		last_stack_files_txt       = [str(cpy_fts_dir) + str(j).split('/')[-1] for j in spec]
		last_stack_files_txt_fnm   = str(cpy_fts_dir) + str(SubSmpl_nm) + '.txt'

		last_stack_files_txt_n     = [str(cpy_fts_dir) + str(j).split('/')[-1] for j in spec_n]
		last_stack_files_txt_fnm_n = str(cpy_fts_dir) + str(SubSmpl_nm) + '_n.txt'

		np.savetxt(last_stack_files_txt_fnm  , last_stack_files_txt  ,  delimiter=" ", fmt="%s", newline='\n')
		np.savetxt(last_stack_files_txt_fnm_n, last_stack_files_txt_n,  delimiter=" ", fmt="%s", newline='\n')

		print 'Stacked files in: '
		print colored(last_stack_files_txt_fnm,'cyan')
		if spc_nse == True:
			print 'Stacked files (noise) in: '
			print colored(last_stack_files_txt_fnm_n,'cyan')
		elif spc_nse == False:
			pass

		Spec_Res = Stack(cat_parent,SubSmpl_nm,spec,spc_nse,
			bs_func         = bs_func             ,
			test_bg         = test_bg             ,test_fg         = test_fg                ,
			upd_header_stck = upd_header_stck     ,new_CRVAL1_head = min(SubSmpl[3][SubSet]),new_CDELT1_head = max(SubSmpl[4][SubSet]),spc_type_2bstck = SubSmpl[6][SubSet],			
			sel_pre_shf     = sel_pre_shf         ,sel_pre_cnt     = sel_pre_cnt            ,sel_pre_msk     = sel_pre_msk            ,
			pre_cnt         = pre_cnt             ,pre_cnt_typ     = pre_cnt_typ            ,pre_cnt_lns     = pre_cnt_lns            ,
			pre_cnt_fnc     = pre_cnt_fnc         ,pre_cnt_ord     = pre_cnt_ord            ,pre_cnt_ovr     = pre_cnt_ovr            ,
			pre_cnt_rpl     = pre_cnt_rpl         ,pre_cnt_lrj     = pre_cnt_lrj            ,pre_cnt_hrj     = pre_cnt_hrj            ,
			smt_spc_pre     = smt_spc_pre         ,smt_shp_pre     = smt_shp_pre            ,smt_sze_pre     = smt_sze_pre            ,
			pre_msk         = pre_msk             ,
			pre_msk_typ     = pre_msk_typ         ,pre_msk_abs_lne = pre_msk_abs_lne        ,pre_msk_cte_val = pre_msk_cte_val,
			pre_msk_rgn     = pre_msk_rgn         ,pre_msk_min     = pre_msk_min            ,pre_msk_max     = pre_msk_max            ,
			wgt_typ         = wgt_typ             ,
			get_cont_flux   = get_cont_flux       ,gcv_lmbd_i      = gcv_lmbd_i             ,gcv_lmbd_f      = gcv_lmbd_f             ,
			sig_clp         = sig_clp             ,
			sig_cut         = sig_cut             ,sig_fct         = sig_fct                ,sig_fll         = sig_fll                ,
			pst_cnt         = pst_cnt             ,pst_cnt_typ     = pst_cnt_typ            ,pst_cnt_lns     = pst_cnt_lns            ,
			pst_cnt_fnc     = pst_cnt_fnc         ,pst_cnt_ord     = pst_cnt_ord            ,pst_cnt_ovr     = pst_cnt_ovr            ,
			pst_cnt_rpl     = pst_cnt_rpl         ,pst_cnt_lrj     = pst_cnt_lrj            ,pst_cnt_hrj     = pst_cnt_hrj            ,
			smt_spc_pst     = smt_spc_pst         ,smt_shp_pst     = smt_shp_pst            ,smt_sze_pst     = smt_sze_pst            ,
			bs_nmb_itr      = bs_nmb_itr          ,
			prt_tbl_stk     = (SubSmpl[0][SubSet]),
			stk_pct_mde     = stk_pct_mde         ,stk_wgt_mde     = stk_wgt_mde)

		if spc_nse == True:
			Spec_Res[0][0] = str(last_stack_files_txt_fnm)
			Spec_Res[0][1] = str(last_stack_files_txt_fnm_n)
			Stack_Spec_Res.append(Spec_Res)
		elif spc_nse == False:
			Spec_Res[0][0] = str(last_stack_files_txt_fnm)
			Stack_Spec_Res.append(Spec_Res[0][0])
	return Stack_Spec_Res

def Bootstrap(bs_tbl,bs_var,bs_i,bs_pct,*args, **kwargs):
	bs_func          = kwargs.get('bs_func'        ,'')
	wrt_fits         = kwargs.get('wrt_fits'       ,True)

	sel_pre_shf     = kwargs.get('sel_pre_shf'     ,True)
	sel_pre_cnt     = kwargs.get('sel_pre_cnt'     ,True)
	sel_pre_msk     = kwargs.get('sel_pre_msk'     ,False)

	pre_cnt         = kwargs.get('pre_cnt'         ,False)
	pre_cnt_typ     = kwargs.get('pre_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pre_cnt_lns     = kwargs.get('pre_cnt_lns'     ,'*')       # Image lines to be fit
	pre_cnt_fnc     = kwargs.get('pre_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pre_cnt_ord     = kwargs.get('pre_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pre_cnt_ovr     = kwargs.get('pre_cnt_ovr'     ,'yes')     # Override previous norm spec
	pre_cnt_rpl     = kwargs.get('pre_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pre_cnt_lrj     = kwargs.get('pre_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pre_cnt_hrj     = kwargs.get('pre_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pre     = kwargs.get('smt_spc_pre'     ,True)
	smt_shp_pre     = kwargs.get('smt_shp_pre'     ,'gaussian')
	smt_sze_pre     = kwargs.get('smt_sze_pre'     ,1)

	pre_msk         = kwargs.get('pre_msk'         ,False)
	pre_msk_typ     = kwargs.get('pre_msk_typ'     ,'NaN')
	pre_msk_cte_val = kwargs.get('pre_msk_cte_val' ,1)
	pre_msk_abs_lne = kwargs.get('pre_msk_abs_lne' ,False)
	pre_msk_rgn     = kwargs.get('pre_msk_rgn' ,False)
	pre_msk_min     = kwargs.get('pre_msk_min' ,500)
	pre_msk_max     = kwargs.get('pre_msk_max' ,1210)

	sig_clp         = kwargs.get('sig_clp'         ,False)
	sig_cut         = kwargs.get('sig_cut'         ,3)
	sig_fct         = kwargs.get('sig_fct'         ,mean)
	sig_fll         = kwargs.get('sig_fll'         ,np.nan)

	wgt_typ         = kwargs.get('wgt_typ'         ,None)
	get_cont_flux   = kwargs.get('get_cont_flux'   ,True)
	gcv_lmbd_i      = kwargs.get('gcv_lmbd_i'      ,1430)
	gcv_lmbd_f      = kwargs.get('gcv_lmbd_f'      ,1480)

	spc_nse         = kwargs.get('spc_nse'         ,False)

	pst_cnt         = kwargs.get('pst_cnt'         ,True)
	pst_cnt_typ     = kwargs.get('pst_cnt_typ'     ,'ratio')   # Continuum fitting type fit,ratio,difference
	pst_cnt_lns     = kwargs.get('pst_cnt_lns'     ,'*')       # Image lines to be fit
	pst_cnt_fnc     = kwargs.get('pst_cnt_fnc'     ,'spline3') # Fitting function: legendre, chebyshev, spline1, spline3
	pst_cnt_ord     = kwargs.get('pst_cnt_ord'     ,49)        # Order Polynomial / num pieces spline
	pst_cnt_ovr     = kwargs.get('pst_cnt_ovr'     ,'yes')     # Override previous norm spec
	pst_cnt_rpl     = kwargs.get('pst_cnt_rpl'     ,'no')      # Replace rejected points by fit?
	pst_cnt_lrj     = kwargs.get('pst_cnt_lrj'     ,3)         # Low rejection in sigma of fit
	pst_cnt_hrj     = kwargs.get('pst_cnt_hrj'     ,3)         # High rejection in sigma of fit

	smt_spc_pst     = kwargs.get('smt_spc_pst'     ,True)
	smt_shp_pst     = kwargs.get('smt_shp_pst'     ,'gaussian')
	smt_sze_pst     = kwargs.get('smt_sze_pst'     ,1)

	cmp_bst_cyc     = kwargs.get('cmp_bst_cyc'     ,False)          # Complete Prrevious Broken BS Cycle
	int_bst_brk     = kwargs.get('int_bst_brk'     ,0)              # Starting iteration point
	lst_bst_brk     = kwargs.get('lst_bst_brk'     ,10)             # Last iteration point

	bst_cyc_otb     = kwargs.get('bst_cyc_otb',False)


	tbl_bst_ipt = readtable_fg_bg_glx(bs_tbl,tbl_format_opt)
	bs_var      = tbl_bst_ipt[0][:]

	if bs_pct == 1:
		bs_s   = np.arange(0,int(len(bs_var)*bs_pct),1)
		bs_r   = 0
	elif bs_pct !=1:
		bs_s   = np.arange(0,int(len(bs_var)*bs_pct),1)
		bs_r   = np.arange(0,len(bs_var) - len(bs_s),1)


	print 
	print 'Bootstrap  : '
	print 'Sample size: ',len(bs_var)
	print 'Iterations : ',bs_i
	print 'Straps     : ',len(bs_s)
	if bs_pct == 1:
		print 'Repetions  : ',0
	elif bs_pct != 1:
		print 'Repetions  : ',len(bs_r)

	print 

	VAR_BS         = []
	TBL_BS_MED     = []
	TBL_BS_AVG     = []
	TBL_BS_STD     = []
	TBL_BS_AVW     = []

	TBL_BS_MED_C   = []
	TBL_BS_AVG_C   = []
	TBL_BS_AVW_C   = []

	TBL_BS_MED_C_S = []
	TBL_BS_AVG_C_S = []
	TBL_BS_AVW_C_S = []

	widgets = ['Bootstrap for galaxies to be stacked: ', Percentage(), ' ', Bar(marker='*',left='[',right=']'),
			   ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

	pbar = ProgressBar(widgets=widgets, maxval=bs_i)
	pbar.start()
	####CLASSIC BOOTSTRAP OR COMPLETE BROKEN BOOTSTRAP#####
	if cmp_bst_cyc == False:
		####CLASSIC BOOTSTRAP
		for main_it in range(bs_i):
			VAR_BS  = []
			for s_it in range(len(bs_s)) : [VAR_BS.append(np.random.choice(bs_s))] 
			VAR_BS_POOL = VAR_BS
			if bs_pct !=1:
				for r_it in range(len(bs_r)) : [VAR_BS.append(np.random.choice(VAR_BS_POOL))]
			elif bs_pct ==1:
				pass
			pbar.update(main_it)
			rtsubsmpl_bs  = tbl_bst_ipt[0][VAR_BS]
			tbl_bs_opt    = str_bst_tbl + (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'.csv').split('/')[-1]
			
			rtsubsmpl_bs.write(tbl_bs_opt,  format=tbl_format_opt, overwrite=True)

			stamps_bootstrap = Select_Subsamples(tbl_bs_opt,None,None,test_fg = True, test_bg = False, slc_int = False,slc_smp = False)
			stacks_bootstrap = np.array(Stack_Subsample(stamps_bootstrap,
								sel_pre_shf     = sel_pre_shf       ,sel_pre_cnt     = sel_pre_cnt           ,sel_pre_msk     = sel_pre_msk,
								pre_cnt         = pre_cnt           ,pre_cnt_typ     = pre_cnt_typ           ,pre_cnt_lns     = pre_cnt_lns,
								pre_cnt_fnc     = pre_cnt_fnc       ,pre_cnt_ord     = pre_cnt_ord           ,pre_cnt_ovr     = pre_cnt_ovr,
								pre_cnt_rpl     = pre_cnt_rpl       ,pre_cnt_lrj     = pre_cnt_lrj           ,pre_cnt_hrj     = pre_cnt_hrj,
								smt_spc_pre     = smt_spc_pre       ,smt_shp_pre     = smt_shp_pre           ,smt_sze_pre     = smt_sze_pre,
								pre_msk         = pre_msk           ,
								pre_msk_typ     = pre_msk_typ       ,pre_msk_abs_lne = False                 ,pre_msk_cte_val = pre_msk_cte_val,
								pre_msk_rgn     = pre_msk_rgn       ,
								pre_msk_min     = pre_msk_min       ,pre_msk_max     = pre_msk_max            ,
								sig_clp         = sig_clp           ,
								sig_cut         = sig_cut           ,sig_fct         = sig_fct               ,sig_fll         = sig_fll,
								wgt_typ         = wgt_typ           ,
								get_cont_flux   = get_cont_flux     ,gcv_lmbd_i      = gcv_lmbd_i            ,gcv_lmbd_f      = gcv_lmbd_f,							
								wrt_fits        = wrt_fits          ,spc_nse         = spc_nse               ,
								pst_cnt         = pst_cnt           ,pst_cnt_typ     = pst_cnt_typ           ,pst_cnt_lns     = pst_cnt_lns,
								pst_cnt_fnc     = pst_cnt_fnc       ,pst_cnt_ord     = pst_cnt_ord           ,pst_cnt_ovr     = pst_cnt_ovr,
								pst_cnt_rpl     = pst_cnt_rpl       ,pst_cnt_lrj     = pst_cnt_lrj           ,pst_cnt_hrj     = pst_cnt_hrj,
								smt_spc_pst     = smt_spc_pst       ,smt_shp_pst     = smt_shp_pst           ,smt_sze_pst     = smt_sze_pst,
								bs_nmb_itr      = bs_i))

			bs_file_med      = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med.fits'
			bs_file_avg      = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg.fits'
			bs_file_std      = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'std.fits'
			bs_file_avw      = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw.fits'

			bs_file_med_c    = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med-c.fits'
			bs_file_avg_c    = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg-c.fits'
			bs_file_avw_c    = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw-c.fits'

			bs_file_med_c_s  = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med-c-smt.fits'
			bs_file_avg_c_s  = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg-c-smt.fits'
			bs_file_avw_c_s  = (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw-c-smt.fits'

			TBL_BS_MED.append(bs_file_med)
			TBL_BS_AVG.append(bs_file_avg)
			TBL_BS_STD.append(bs_file_std)
			TBL_BS_AVW.append(bs_file_avw)

			TBL_BS_MED_C.append(bs_file_med_c)
			TBL_BS_AVG_C.append(bs_file_avg_c)
			TBL_BS_AVW_C.append(bs_file_avw_c)

			TBL_BS_MED_C_S.append(bs_file_med_c_s)
			TBL_BS_AVG_C_S.append(bs_file_avg_c_s)
			TBL_BS_AVW_C_S.append(bs_file_avw_c_s)

		rtM                        = aptbl.Table()
		rtM['bs_spc_file_med']     = TBL_BS_MED
		rtM['bs_spc_file_avg']     = TBL_BS_AVG
		rtM['bs_spc_file_std']     = TBL_BS_STD
		rtM['bs_spc_file_avw']     = TBL_BS_AVW
		rtM['bs_spc_file_med-c']   = TBL_BS_MED_C
		rtM['bs_spc_file_avg-c']   = TBL_BS_AVG_C
		rtM['bs_spc_file_avw-c']   = TBL_BS_AVW_C
		rtM['bs_spc_file_med-c-s'] = TBL_BS_MED_C_S
		rtM['bs_spc_file_avg-c-s'] = TBL_BS_AVG_C_S
		rtM['bs_spc_file_avw-c-s'] = TBL_BS_AVW_C_S

		tbl_bs_mst                 = str_bst_tbl + (str(bs_tbl.split('.csv')[0]) + '-BS_MST_' +str(bs_i) + '.' + tbl_format_opt).split('/')[-1]
		print
		rtM.write(tbl_bs_mst, format=tbl_format_opt,overwrite=True)
		print 'Results containing Bootstrap galaxies in table: '
		print tbl_bs_mst
		
		pbar.finish()

		print ' Tables located in directory: ',str_bst_tbl		

		return(tbl_bs_mst)
		####CLASSIC BOOTSTRAP
	####BROKEN BOOTSTRAP######
	elif cmp_bst_cyc == True and bst_cyc_otb == False:
		print 
		print colored('Bootstrap Restoring Mode! ' + str(cmp_bst_cyc),'magenta')
		print colored('Completing previous broken Bootstrap Cycle! ' + str(cmp_bst_cyc),'yellow')
		print
		print colored('Restartiing from initial Value : ' + str(int_bst_brk),'yellow')
		print colored('Last BS Value                  : ' + str(lst_bst_brk),'yellow')
		print
		
		for main_it in range(int_bst_brk-1,lst_bst_brk+1):
			print 
			print colored('Bootstrap Restoring Mode! ' + str(cmp_bst_cyc),'magenta')
			print colored('Completing previous broken Bootstrap Cycle! ' + str(cmp_bst_cyc),'yellow')
			print
			print colored('Restartiing from initial Value : ' + str(int_bst_brk),'yellow')
			print colored('Last BS Value                  : ' + str(lst_bst_brk),'yellow')
			print
			
			VAR_BS  = []
			for s_it in range(len(bs_s)) : [VAR_BS.append(np.random.choice(bs_s))] 
			VAR_BS_POOL = VAR_BS
			if bs_pct !=1:
				for r_it in range(len(bs_r)) : [VAR_BS.append(np.random.choice(VAR_BS_POOL))]
			elif bs_pct ==1:
				pass
			pbar.update(main_it)
			rtsubsmpl_bs  = tbl_bst_ipt[0][VAR_BS]
			tbl_bs_opt    = str_bst_tbl + (str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'.csv').split('/')[-1]
			
			rtsubsmpl_bs.write(tbl_bs_opt,  format=tbl_format_opt, overwrite=True)

			stamps_bootstrap = Select_Subsamples(tbl_bs_opt,None,None,test_fg = True, test_bg = False, slc_int = False,slc_smp = False)
			stacks_bootstrap = np.array(Stack_Subsample(stamps_bootstrap,
								sel_pre_shf     = sel_pre_shf       ,sel_pre_cnt     = sel_pre_cnt           ,sel_pre_msk     = sel_pre_msk,
								pre_cnt         = pre_cnt           ,pre_cnt_typ     = pre_cnt_typ           ,pre_cnt_lns     = pre_cnt_lns,
								pre_cnt_fnc     = pre_cnt_fnc       ,pre_cnt_ord     = pre_cnt_ord           ,pre_cnt_ovr     = pre_cnt_ovr,
								pre_cnt_rpl     = pre_cnt_rpl       ,pre_cnt_lrj     = pre_cnt_lrj           ,pre_cnt_hrj     = pre_cnt_hrj,
								smt_spc_pre     = smt_spc_pre       ,smt_shp_pre     = smt_shp_pre           ,smt_sze_pre     = smt_sze_pre,
								pre_msk         = pre_msk           ,
								pre_msk_typ     = pre_msk_typ       ,pre_msk_abs_lne = False                 ,pre_msk_cte_val = pre_msk_cte_val,
								pre_msk_rgn     = pre_msk_rgn       ,
								pre_msk_min     = pre_msk_min       ,pre_msk_max     = pre_msk_max           ,
								sig_clp         = sig_clp           ,
								sig_cut         = sig_cut           ,sig_fct         = sig_fct               ,sig_fll         = sig_fll,
								wgt_typ         = wgt_typ           ,
								get_cont_flux   = get_cont_flux     ,gcv_lmbd_i      = gcv_lmbd_i            ,gcv_lmbd_f      = gcv_lmbd_f,							
								wrt_fits        = wrt_fits          ,spc_nse         = spc_nse               ,
								pst_cnt         = pst_cnt           ,pst_cnt_typ     = pst_cnt_typ           ,pst_cnt_lns     = pst_cnt_lns,
								pst_cnt_fnc     = pst_cnt_fnc       ,pst_cnt_ord     = pst_cnt_ord           ,pst_cnt_ovr     = pst_cnt_ovr,
								pst_cnt_rpl     = pst_cnt_rpl       ,pst_cnt_lrj     = pst_cnt_lrj           ,pst_cnt_hrj     = pst_cnt_hrj,
								smt_spc_pst     = smt_spc_pst       ,smt_shp_pst     = smt_shp_pst           ,smt_sze_pst     = smt_sze_pst,
								bs_nmb_itr      = bs_i))
			
		[TBL_BS_MED.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med.fits') for main_it in range(bs_i)]#bs_file_med)
		[TBL_BS_AVG.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg.fits') for main_it in range(bs_i)]#bs_file_avg)
		[TBL_BS_STD.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'std.fits') for main_it in range(bs_i)]#bs_file_std)
		[TBL_BS_AVW.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw.fits') for main_it in range(bs_i)]#bs_file_avw)

		[TBL_BS_MED_C.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med-c.fits') for main_it in range(bs_i)]#bs_file_med_c)
		[TBL_BS_AVG_C.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg-c.fits') for main_it in range(bs_i)]#bs_file_avg_c)
		[TBL_BS_AVW_C.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw-c.fits') for main_it in range(bs_i)]#bs_file_avw_c)

		[TBL_BS_MED_C_S.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med-c-smt.fits') for main_it in range(bs_i)]#bs_file_med_c_s)
		[TBL_BS_AVG_C_S.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg-c-smt.fits') for main_it in range(bs_i)]#bs_file_avg_c_s)
		[TBL_BS_AVW_C_S.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw-c-smt.fits') for main_it in range(bs_i)]#bs_file_avw_c_s)
		
		rtM                        = aptbl.Table()
		rtM['bs_spc_file_med']     = TBL_BS_MED
		rtM['bs_spc_file_avg']     = TBL_BS_AVG
		rtM['bs_spc_file_std']     = TBL_BS_STD
		rtM['bs_spc_file_avw']     = TBL_BS_AVW
		rtM['bs_spc_file_med-c']   = TBL_BS_MED_C
		rtM['bs_spc_file_avg-c']   = TBL_BS_AVG_C
		rtM['bs_spc_file_avw-c']   = TBL_BS_AVW_C
		rtM['bs_spc_file_med-c-s'] = TBL_BS_MED_C_S
		rtM['bs_spc_file_avg-c-s'] = TBL_BS_AVG_C_S
		rtM['bs_spc_file_avw-c-s'] = TBL_BS_AVW_C_S

		
		tbl_bs_mst                 = str_bst_tbl + (str(bs_tbl.split('.csv')[0]) + '-BS_MST_' +str(bs_i) + '.' + tbl_format_opt).split('/')[-1]
		print
		rtM.write(tbl_bs_mst, format=tbl_format_opt,overwrite=True)
		print colored('No Output table in BS Restoring Mode!','yellow')
		print colored('Results containing Bootstrap galaxies in table: ','green')
		print colored(tbl_bs_mst,'green')
		
		
		pbar.finish()

		print ' Tables located in directory: ',str_bst_tbl		

		return(tbl_bs_mst)
	####BROKEN BOOTSTRAP######
	####BROKEN BOOTSTRAP-ONLY TABLE######
	elif cmp_bst_cyc == True and bst_cyc_otb == True:
		print 
		print colored('Bootstrap Restoring Mode! ' + str(cmp_bst_cyc),'magenta')
		print colored('Completing previous broken Bootstrap Cycle! ' + str(cmp_bst_cyc),'yellow')
		print
		print colored('Restartiing from initial Value : ' + str(int_bst_brk),'yellow')
		print colored('Last BS Value                  : ' + str(lst_bst_brk),'yellow')
		print
		
		[TBL_BS_MED.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med.fits') for main_it in range(bs_i)]#bs_file_med)
		[TBL_BS_AVG.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg.fits') for main_it in range(bs_i)]#bs_file_avg)
		[TBL_BS_STD.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'std.fits') for main_it in range(bs_i)]#bs_file_std)
		[TBL_BS_AVW.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw.fits') for main_it in range(bs_i)]#bs_file_avw)

		[TBL_BS_MED_C.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med-c.fits') for main_it in range(bs_i)]#bs_file_med_c)
		[TBL_BS_AVG_C.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg-c.fits') for main_it in range(bs_i)]#bs_file_avg_c)
		[TBL_BS_AVW_C.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw-c.fits') for main_it in range(bs_i)]#bs_file_avw_c)

		[TBL_BS_MED_C_S.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'med-c-smt.fits') for main_it in range(bs_i)]#bs_file_med_c_s)
		[TBL_BS_AVG_C_S.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avg-c-smt.fits') for main_it in range(bs_i)]#bs_file_avg_c_s)
		[TBL_BS_AVW_C_S.append((str(bs_tbl.split('.csv')[0]) + '-BS-'+str(main_it+1)+'-stk-').split('/')[-1] + 'avw-c-smt.fits') for main_it in range(bs_i)]#bs_file_avw_c_s)
			

		
		rtM                        = aptbl.Table()
		rtM['bs_spc_file_med']     = TBL_BS_MED
		rtM['bs_spc_file_avg']     = TBL_BS_AVG
		rtM['bs_spc_file_std']     = TBL_BS_STD
		rtM['bs_spc_file_avw']     = TBL_BS_AVW
		rtM['bs_spc_file_med-c']   = TBL_BS_MED_C
		rtM['bs_spc_file_avg-c']   = TBL_BS_AVG_C
		rtM['bs_spc_file_avw-c']   = TBL_BS_AVW_C
		rtM['bs_spc_file_med-c-s'] = TBL_BS_MED_C_S
		rtM['bs_spc_file_avg-c-s'] = TBL_BS_AVG_C_S
		rtM['bs_spc_file_avw-c-s'] = TBL_BS_AVW_C_S

		
		tbl_bs_mst                 = str_bst_tbl + (str(bs_tbl.split('.csv')[0]) + '-BS_MST_' +str(bs_i) + '.' + tbl_format_opt).split('/')[-1]
		print
		rtM.write(tbl_bs_mst, format=tbl_format_opt,overwrite=True)
		print colored('Output table in BS Restoring Mode!','yellow')
		print colored('Results containing Bootstrap galaxies in table: ','green')
		print colored(tbl_bs_mst,'green')
		
		
		pbar.finish()

		print ' Tables located in directory: ',str_bst_tbl		

		return(tbl_bs_mst)
		####BROKEN BOOTSTRAP######
	####BROKEN BOOTSTRAP-ONLY TABLE######
	####CLASSIC BOOTSTRAP OR COMPLETE BROKEN BOOTSTRAP#####
####Fnc_Stk_Stk####