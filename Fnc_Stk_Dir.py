import sys, os
import numpy as np
from numpy import mean,median
from progressbar import *
from termcolor import colored

from Fnc_Stk_Utl import *

home = os.path.expanduser("~") + '/Desktop/Example/'

############################################################################
###############################Input Catalogue##############################
#Catalogue
cat_parent             	= 'cosmos' #cosmos vvds2h ecdfs all_0
CAT_PARENT             	= cat_parent.upper()

#Spectra Directories
cats_dir               	= home     + 'Cats/COSMOS/'
cat_dir                	= cats_dir + CAT_PARENT + '/'
spc_dir                	= cat_dir  + CAT_PARENT + '/spec1d/'
spc_dirn               	= cat_dir  + CAT_PARENT + '/spec1dnoise/'

#Input Table
sfx_tbl_ext_opt        	= ''#'_PRP_MRP' #'-PRP' ''
cat_ipt_tbl            	= cat_dir + 'cesam_' + cat_parent + '_spectra_' + cat_parent + '_catalog' + sfx_tbl_ext_opt
cat_ipt_tbl            	= cat_dir + 'cesam_' + cat_parent + '_spectra_' + cat_parent + '_catalog'

#tables
tbl_format_ipt         	= 'csv' 										#'ascii.fixed_width_two_line'       #ascii,csv,fits,ascii.fixed_width_two_line
tbl_format_opt         	= 'csv' 										#'ascii.fixed_width_two_line'       #ascii,csv,fits,ascii.fixed_width_two_line
###############################Input Catalogue##############################
############################################################################

############################################################################
###############################Pair Selection###############################
slt_prs                	= True                                      	#Select pairs
unq_tbl_ipt            	= True											#Search unique elements input table
head_un_ipt            	= 'ident'										#Table header used for uniqueness 'spec1d' for select pairs
copy_spc_files         	= True                               		    #Copy spec files
copy_spn_files         	= False                               		    #Copy spec noise files
z_lim_inf              	= 1.5                                       	#lower redshift limit for pairs
zl1                    	= 0.01                                      	#lower redshift limit for selection
zl2                    	= 1.00                                      	#upper redshift limit for selection
rad_vel_sep            	= 3000                                      	#km/s
rad_sep                	= [[0,40],[0,5,10,15,20]]                   	#arcsec
rad_sep_23             	= [[0,40],[0,23]]                           	#arcsec
rad_sep_23_2           	= [[0,40],[0,16.5,23]]                      	#arcsec
rad_sep_23_3           	= [[0,40],[0,13.6,19.0,23]]                 	#arcsec
rad_sep_23_4           	= [[0,40],[0,11.8,16.5,20,23]]              	#arcsec Removed: 27,30
###############################Pair Selection###############################
############################################################################

############################################################################ 
###############################Cuts SubSamples############################## 
#Cuts: redshift_fg, redshift_fg_flag, redshift_bk, redshift_bk_flag, deltaz, sep_as, sep_kpc
z_itv                  	= [1.5,3.6]                                 	#[0,0.5,1,1.5,2,2.5,3,3.5,4,5]
z_flag_itv_fg          	= [3,4,33,34,43,44]                         	#33,34,43,44 Should be increasing
z_flag_itv_bg          	= [3,4,33,34,43,44]                         	#33,34,43,44 Should be increasing
dz_itv                 	= [0,0.2,0.4,0.6,0.8,1.0]                  		#
SEP_as_itv             	= rad_sep[1]                                	#
SEP_as_itv_23          	= rad_sep_23[1]
SEP_as_itv_23_4        	= rad_sep_23_4[1]
SEP_as_itv_23_2        	= rad_sep_23_2[1]
SEP_as_itv_23_3        	= rad_sep_23_3[1]
SEP_as_itv_23_4        	= rad_sep_23_4[1]
SEP_kpc_itv            	= []  
###############################Cuts SubSamples##############################
############################################################################

############################################################################
#############################Stacking Parameters############################
#SPECTRA SELECTION
selec_spec_shift       = True                                      # Shift spectra
selec_spec_contn       = True                                      # Continuum Fitting/Normalization
selec_spec_masks       = False                                     # Mask

#Stacks Pre-Processing 
#Pre-Processing Continuum
pre_continuum          = False                                     # Continuum Fitting/Normalization
pre_cont_typ           = 'ratio'                                   # Continuum fitting type fit,ratio,difference
pre_cont_lines         = '*'                                       # Image lines to be fit
pre_cont_funct         = 'spline3'                                 # Fitting function: legendre, chebyshev, spline1, spline3
pre_cont_order         = 49                                        # Order Polynomial / num pieces spline
pre_cont_override      = 'yes'                                     # Override previous norm spec
pre_cont_replace       = 'no'                                      # Replace rejected points by fit?
pre_cont_low_rej       = 3                                         # Low rejection in sigma of fit
pre_cont_high_rej      = 3                                         # High rejection in sigma of fit

#Pre-Processing Smoothing
pre_smooth             = True                                      # smooth after interpolation and before stacking
pre_smooth_shape       = 'gaussian'                                # gaussian,boxcar,mexican
pre_smooth_size        = 1                                         # kernel size pixel units

#Pre-Processing MASKING
pre_mask               = True                                      # mask spectra after smoothing (stacks)
pre_msk_abs_lines      = True                                      # mask IS absorptions lines
pre_mask_type          = 'NaN'                                     # continuum/constant/NaN
pre_mask_cte_val       = 0                                         # constant value for masking
pre_mask_lw            = 2                                         # line width (A)
pre_mask_regn          = True                                      # mask initial spectra pixels
pre_mask_regn_int      = 300                                       # intial pix
pre_mask_regn_fnl      = 912                                       # final pix

#Sigma-Clip
sigma_clipping         = True                                      # Sigma clipping
sigma_cut              = 3                                         # sigma cut
sigma_cen_fct          = mean                                      # median, mean
sigma_msk_fill_val     = np.nan                                    # np.nan, value

#Weighting
weight_type            = 'cont-flux-med'                           # i-band-mag,cont-flux-sum,cont-flux-med,cont-flux-avg None:
weight_cnt_flux_get    = True                                      # mask any given wavelength region 
weight_cnt_flux_lmb_0  = 1430                                      # initial lambda
weight_cnt_flux_lmb_n  = 1480                                      # final lambda

#Noise Files
spectra_noise          = False                                     #Include Noise files in the Stacks

#Stacks Post Processing 
post_continuum         = True	                                   # Fit Cont after stacking
post_cont_typ          = 'ratio'                                   # Continuum fitting type fit,ratio,difference
post_cont_lines        = '*'                                       # Image lines to be fit
post_cont_funct        = 'spline3'                                 # Fitting function: legendre, chebyshev, spline1, spline3
post_cont_order        = 49                                        # Order Polynomial / num pieces spline
post_cont_override     = 'yes'                                     # Override previous norm spec
post_cont_replace      = 'no'                                      # Replace rejected points by fit?
post_cont_low_rej      = 5                                         # Low rejection in sigma of fit
post_cont_high_rej     = 5                                         # High rejection in sigma of fit
post_smooth            = True                                      # smooth after stacking
post_smooth_shape      = 'gaussian'                                # smooth after stacking
post_smooth_size       = 1                                         # smooth after stacking
#############################Stacking Parameters############################
############################################################################

############################################################################ 
##################################Bootstrap#################################
bs_iteration_num       = 2                                       # Bootstrap iterations
bs_percentage          = 1                                         # Bootstrap fraction to resample #0.60
bs_function            = '_avw-c'                                  # '_avg','_avw','_med'

#BS broken cycle
comp_BS_run			   = False									   #Complete Broken BS Cycle iterattion
bst_brk_int			   = 2										   #Init iteration step
bst_brk_lst			   = bs_iteration_num						   #Last iteration step (default = bs_iteration_num)
crts_BS_otb			   = False									   #Creates BS-MST table Only For Completed BS Repetitions = bs_iteration_num

##################################Bootstrap#################################
############################################################################

###############################DIRECTORIES###############################
#MAIN DIRECTORIES#
stk_dir_res       		= home + 'Stack_Results/'
cat_dir_res       		= stk_dir_res + CAT_PARENT
par_dir_res       		= cat_dir_res + '/PAIRS' 
par_frg_dir       		= par_dir_res + '/FRGRD/'
par_bkg_dir       		= par_dir_res + '/BKGRD/'

#TABLES#
tbl_dir_res       		= cat_dir_res + '/TABLES'
frg_dir_res       		= tbl_dir_res + '/FRGRD/'
bkg_dir_res       		= tbl_dir_res + '/BKGRD/'
stt_dir_res       		= tbl_dir_res + '/STATS/'
irf_dir_res       		= tbl_dir_res + '/IRAF/'
std_dir_res       		= tbl_dir_res + '/Literature/'
plt_tbl_res       		= tbl_dir_res + '/PLOTS/'

#PLOTS#
plt_dir_res       		= cat_dir_res + '/PLOTS' 
ind_plt_res       		= plt_dir_res + '/IND/' 
frg_ind_plt       		= ind_plt_res + '/FRGRD/' 
bkg_ind_plt       		= ind_plt_res + '/BKGRD/' 
res_plt_res       		= plt_dir_res + '/RESULTS/' 
ewr_plt_res       		= plt_dir_res + '/EW/' 
fit_plt_res       		= plt_dir_res + '/FIT/' 

#STACKS#
spc_stk_res       		= cat_dir_res + '/STACKS'
res_stk_res       		= spc_stk_res + '/RESULTS/'
lst_stk_res       		= spc_stk_res + '/LAST-FITS'
ind_stk_res       		= lst_stk_res + '/RESULTS/'

#BOOTSRAPS#
bst_dir_res       		= cat_dir_res + '/BOOTSTRAP'
tbl_bst_res       		= bst_dir_res + '/TABLES/'
str_bst_tbl       		= tbl_bst_res + '/STRAPS/'
stt_bst_tbl       		= tbl_bst_res + '/STATS/'

#BOOTSRAPS-STACKS#
stk_bst_res       		= bst_dir_res + '/STACKS'
stt_bst_stk       		= stk_bst_res + '/STATS-BST/'
str_bst_stk       		= stk_bst_res + '/STATS-STR/'
lst_bst_res       		= stk_bst_res + '/LAST-FITS'
ind_bst_lst       		= lst_bst_res + '/STATS-BST/' 
fts_bst_lst       		= lst_bst_res + '/STATS-STR/'

#BOOTSRAPS-PLOTS#
plt_bst_res       		= bst_dir_res + '/PLOTS'
res_bst_plt       		= plt_bst_res + '/RESULTS/'
ind_bst_plt       		= plt_bst_res + '/INDIVIDUAL-BS/'
###############################DIRECTORIES###############################

#Output  Tables Selected Pairs
individual_table  		= True                               					#Create PAIRS gral table #PREV ind_tbl
general_table     		= True                               					#Create PAIRS indv table #PREV grl_tbl
grl_tbl_nmB       		= 'P_Bg_'+ CAT_PARENT 
grl_tbl_nmF       		= 'P_Fg_'+ CAT_PARENT 

#Output  Tables Selected Pairs Repeated elements
unq_tbl_opt       		= False													#Search unique elements output table
hed_un_opt_F      		= 'id_F'                                				#Table header used for uniqueness 'id_B' Select Pairs
hed_un_opt_B      		= 'id_B'                                				#Table header used for uniqueness
grl_tbl_nmB_U     		= 'P_Bg_'+ CAT_PARENT +'_U'			  					#+'_U_PRP' Select pairs
grl_tbl_nmF_U     		= 'P_Fg_'+ CAT_PARENT +'_U'			  					#+'_U_PRP' Select pairs

verbose_gral           = False
#############################################################################################################################
DIR_CAT_IPT = [cats_dir]
DIR_SPC_IPT = [spc_dir,spc_dirn]
DIR_RES     = [stk_dir_res,cat_dir_res,par_dir_res,par_frg_dir,par_bkg_dir,
				tbl_dir_res,frg_dir_res,bkg_dir_res,stt_dir_res,
				plt_dir_res,ind_plt_res,frg_ind_plt,bkg_ind_plt,res_plt_res,
				ewr_plt_res,fit_plt_res,irf_dir_res,std_dir_res,plt_tbl_res,
				spc_stk_res,res_stk_res,lst_stk_res,ind_stk_res,
				bst_dir_res,tbl_bst_res,stt_bst_tbl,str_bst_tbl,
				stk_bst_res,stt_bst_stk,str_bst_stk,
				lst_bst_res,ind_bst_lst,fts_bst_lst,
				plt_bst_res,res_bst_plt,ind_bst_plt]

if tbl_format_ipt == 'ascii' or tbl_format_ipt == 'ascii.fixed_width_two_line':
	tbl_ext_ipt = '.dat'
elif tbl_format_ipt == 'csv':	
	tbl_ext_ipt = '.csv'
elif tbl_format_ipt == 'fits':	
	tbl_ext_ipt = '.fits'

if tbl_format_opt == 'ascii' or tbl_format_opt == 'ascii.fixed_width_two_line':
	tbl_ext_opt = '.dat'
elif tbl_format_opt == 'csv':	
	tbl_ext_opt = '.csv'
elif tbl_format_opt == 'fits':	
	tbl_ext_opt = '.fits'

cat_tbl    = cat_ipt_tbl + tbl_ext_ipt
cat_tbl_U  = cat_ipt_tbl + tbl_ext_opt

op_tbl_B   = bkg_dir_res + grl_tbl_nmB   + '_' + str(rad_sep[0][0]) + '-' + str(rad_sep[0][-1]) + sfx_tbl_ext_opt + tbl_ext_opt      #header of column for uniqueness
op_tbl_F   = frg_dir_res + grl_tbl_nmF   + '_' + str(rad_sep[0][0]) + '-' + str(rad_sep[0][-1]) + sfx_tbl_ext_opt + tbl_ext_opt      #header of column for uniqueness

op_tbl_B_U = bkg_dir_res + grl_tbl_nmB   + '_' + str(rad_sep[0][0]) + '-' + str(rad_sep[0][-1]) + '_U' + sfx_tbl_ext_opt + tbl_ext_opt      #header of column for uniqueness
op_tbl_F_U = frg_dir_res + grl_tbl_nmF   + '_' + str(rad_sep[0][0]) + '-' + str(rad_sep[0][-1]) + '_U' + sfx_tbl_ext_opt + tbl_ext_opt      #header of column for uniqueness

#############################################################################################################################
def Check_directories(cat_tbl_chk,cat_chk,*args, **kwargs):
	verbose = kwargs.get('verbose',False)
	if os.path.exists(str(cat_tbl_chk)+str(tbl_ext_ipt))==True and os.path.exists(DIR_SPC_IPT[0])==True and os.path.exists(DIR_SPC_IPT[1])==True:
		if verbose == True:
			print
			print 'Checking input directories of the catalogue : ',cat_chk
			print
			print 'Catalogue table exists             : ', str(cat_tbl_chk)+str(tbl_ext_ipt)
			print 'Spectra directories exists         : ', str(DIR_SPC_IPT[0])
			print 'Spectra directories exists (noise) : ', str(DIR_SPC_IPT[1])
			print
			print 'Checking Result Directories.'
			print
		elif verbose == False:
			pass

		for tree in DIR_RES:
			if os.path.isdir(tree)==True:
				pass
				print 'Directory exists: ', tree
			elif os.path.isdir(tree)==False:
				print 'Directory does not exist, creating it: ', tree
				os.makedirs(tree)
	elif os.path.exists(str(cat_tbl_chk)+str(tbl_ext_ipt))==False or os.path.exists(DIR_SPC_IPT[0])==False or os.path.exists(DIR_SPC_IPT[1])==False:
		print
		print colored('Some of essential the directories does not exist.','yellow')
		print colored('Check input directories: ','yellow')
		print
		print colored(str(cat_tbl_chk)+str(tbl_ext_ipt),'yellow')
		print colored(DIR_SPC_IPT[0],'yellow')
		print colored(DIR_SPC_IPT[1],'yellow')
		print
#############################################################################################################################


fg_sbdir = Def_Sub_Dirs_Slice_xtr(par_frg_dir,rad_sep[0])[0]
bg_sbdir = Def_Sub_Dirs_Slice_xtr(par_bkg_dir,rad_sep[0])[0]

fg_sbdir_ind_plt = Def_Sub_Dirs_Slice_all(frg_ind_plt,rad_sep[1])[0]
bg_sbdir_ind_plt = Def_Sub_Dirs_Slice_all(bkg_ind_plt,rad_sep[1])[0]


DIR_RES.extend(fg_sbdir)
DIR_RES.extend(bg_sbdir)

DIR_RES.extend(fg_sbdir_ind_plt)
DIR_RES.extend(bg_sbdir_ind_plt)

Check_directories(cat_ipt_tbl,cat_parent,tbl_ext_ipt=tbl_ext_ipt)

