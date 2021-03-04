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

from Fnc_Stk_Utl import *

home = expanduser("~") + '/Desktop/VUDS/'
home = expanduser("~") + '/VUDS/'

#Catalogue
cat_parent             = 'all_0' #cosmos vvds2h ecdfs all
CAT_PARENT             = cat_parent.upper()

#Spectra Directories
cats_dir               = home     + 'Cats/Dez2016/'
cat_dir                = cats_dir + CAT_PARENT + '/'
spc_dir                = cat_dir  + CAT_PARENT + '/spec1d/'
spc_dirn               = cat_dir  + CAT_PARENT + '/spec1dnoise/'

#Input Table
sfx_tbl_ext_opt        = '-PRP'#'_PRP_MRP' #'-PRP' ''
cat_ipt_tbl            = cat_dir + 'cesam_vuds_spectra_' + cat_parent + '_catalog' + sfx_tbl_ext_opt
cat_ipt_tbl            = cat_dir + 'cesam_vuds_spectra_' + cat_parent + '_catalog'
unq_tbl_ipt            = True
head_un_ipt            = 'ident'#'spec1d'

#tables
tbl_format_ipt         = 'csv' #'ascii.fixed_width_two_line'       #ascii,csv,fits,ascii.fixed_width_two_line
tbl_format_opt         = 'csv' #'ascii.fixed_width_two_line'       #ascii,csv,fits,ascii.fixed_width_two_line

#Pair Selection 
slt_prs                = True                                      #Select pairs
z_lim_inf              = 1.5                                       #lower redshift limit for pairs
zl1                    = 0.01                                      #lower redshift limit for selection
zl2                    = 1.00                                      #upper redshift limit for selection
rad_vel_sep            = 3000                                      #km/s
rad_sep                = [[0,40],[0,5,10,15,20]]                   #arcsec
rad_sep_20             = [[0,40],[0,20]]                           #arcsec
rad_sep_23             = [[0,40],[0,23]]                           #arcsec
rad_sep_25             = [[0,40],[0,25]]                           #arcsec
rad_sep_20_2           = [[0,40],[0,14.5,20]]                      #arcsec
rad_sep_20_3           = [[0,40],[0,11.8,16.5,20]]                 #arcsec
rad_sep_20_4           = [[0,40],[0,5,10,15,20]]                   #arcsec
rad_sep_23_2           = [[0,40],[0,16.5,23]]                      #arcsec
rad_sep_23_3           = [[0,40],[0,13.6,19.0,23]]                 #arcsec
rad_sep_23_4           = [[0,40],[0,11.8,16.5,20,23]]              #arcsec Removed: 27,30
rad_sep_25_2           = [[0,40],[0,18.0,25]]                      #arcsec
rad_sep_25_3           = [[0,40],[0,15.0,21.0,25]]                 #arcsec
rad_sep_25_5           = [[0,40],[0,11.8,16.5,20,23,25]]           #arcsec Removed: 27,30
rad_sep_25_5_1         = [[0,40],[0,5,10,15,20,25]]                #arcsec Removed: 30
rad_sep_25_5_2         = [[0,40],[0,11.8,16.5,20,23,25]]           #arcsec Removed: 27,30
##rad_sep                = [[0,40],[20,25,30,35,40]]                #arcsec
##rad_sep                = [[0,40],[0,5,10,15,20,25,30,35,40]]      #arcsec
##rad_sep                = [[0,40],[20,40]]                         #arcsec
##rad_sep                = [[0,40],[0,20]]                          #arcsec
#cp_spc_fls             = True                                      #Copy spec files
#cp_spn_fls             = True                                      #Copy spec noise files

#Cuts: redshift_fg, redshift_fg_flag, redshift_bk, redshift_bk_flag, deltaz, sep_as, sep_kpc
z_itv                  = [1.5,3.6]                                 #[0,0.5,1,1.5,2,2.5,3,3.5,4,5]
z_flag_itv_fg          = [3,4,33,34,43,44]                         #33,34,43,44 Should be increasing
z_flag_itv_bg          = [3,4,33,34,43,44]                         #33,34,43,44 Should be increasing
#z_flag_itv_bg          = [3,4,9,33,34,39,43,44,49,93,94,99,349]    #
#z_flag_itv_bg          = [2,3,4,9,22,23,24,29,32,33,34,39,42,43,44,49,92,93,94,99,3492]                  #
#z_flag_itv_bg          = [1,2,3,4,9,11,12,13,14,19,21,22,23,24,29,31,32,33,34,39,41,42,43,44,49,91,92,93,94,99,34921]                  #
dz_itv                 = [0,0.2,0.4,0.6,0.8,1.0]                   #
SEP_as_itv             = rad_sep[1]                                #
SEP_as_itv_20          = rad_sep_20[1]
SEP_as_itv_23          = rad_sep_23[1]
SEP_as_itv_25          = rad_sep_25[1]
SEP_as_itv_20_2        = rad_sep_20_2[1]
SEP_as_itv_20_3        = rad_sep_20_3[1]
SEP_as_itv_20_4        = rad_sep_20_4[1]
SEP_as_itv_23_4        = rad_sep_23_4[1]
SEP_as_itv_23_2        = rad_sep_23_2[1]
SEP_as_itv_23_3        = rad_sep_23_3[1]
SEP_as_itv_23_4        = rad_sep_23_4[1]
SEP_as_itv_25_2        = rad_sep_25_2[1]
SEP_as_itv_25_3        = rad_sep_25_3[1]
SEP_as_itv_25_5        = rad_sep_25_5[1]
SEP_as_itv_25_5_1      = rad_sep_25_5_1[1]
SEP_as_itv_25_5_2      = rad_sep_25_5_2[1]
SEP_kpc_itv            = []  

mass_itv               = [8,9.50,9.72,9.9,14]                      #4 bins 45,46,46,44/per bin for z_flag 3,4,33,34,43,44
mass_itv_20_2          = [8,9.72,14]                               #2 bins 96,85/per bin for z_flag 3,4,33,34,43,44
mass_itv_20_3          = [8,9.57,9.90,14]                          #3 bins 58,61,62/per bin for z_flag 3,4,33,34,43,44
mass_itv_20_4          = [8,9.50,9.72,9.9,14]                      #4 bins 45,46,46,44/per bin for z_flag 3,4,33,34,43,44
mass_itv_23_2          = [8,9.72,14]                               #2 bins 96,85/per bin for z_flag 3,4,33,34,43,44
mass_itv_23_3          = [8,9.54,9.88,14]                          #3 bins 58,61,62/per bin for z_flag 3,4,33,34,43,44
mass_itv_23_4          = [8,9.50,9.72,10.02,14]                    #4 bins 45,46,46,44/per bin for z_flag 3,4,33,34,43,44
mass_itv_25_2          = [8,9.72,14]                               #2 bins 96,85/per bin for z_flag 3,4,33,34,43,44
mass_itv_25_3          = [8,9.53,9.90,14]                          #3 bins 58,61,62/per bin for z_flag 3,4,33,34,43,44
mass_itv_25_4          = [8,9.45,9.70,9.97,14]                     #4 bins 45,46,46,44/per bin for z_flag 3,4,33,34,43,44

Age_itv                = [0,1.70E8,3.50E8,5.50E8,2.31E9]           #4 bins 62,60,59,42/per bin for z_flag 3,4,33,34,43,44
Age_itv_20_2           = [0,3.40E8,2.31E9]                         #2 bins 88,93/per bin for z_flag 3,4,33,34,43,44
Age_itv_20_3           = [0,2.20E8,4.50E8,2.31E9]                  #3 bins 62,60,59/per bin for z_flag 3,4,33,34,43,44
Age_itv_20_4           = [0,1.70e8,3.50E8,5.50E8,2.31E9]           #4 bins 62,60,59,42/per bin for z_flag 3,4,33,34,43,44
Age_itv_23_2           = [0,3.41E8,2.31E9]                         #2 bins 88,93/per bin for z_flag 3,4,33,34,43,44
Age_itv_23_3           = [0,2.10E8,4.50E8,2.31E9]                  #3 bins 62,60,59/per bin for z_flag 3,4,33,34,43,44
Age_itv_23_4           = [0,1.74E8,3.43E8,5.32e8,2.31e9]            #4 bins 62,60,59,42/per bin for z_flag 3,4,33,34,43,44
Age_itv_25_2           = [0,3.40E8,2.31E9]                         #2 bins 88,93/per bin for z_flag 3,4,33,34,43,44
Age_itv_25_3           = [0,2.00E8,4.50E8,2.31E9]                  #3 bins 62,60,59/per bin for z_flag 3,4,33,34,43,44
Age_itv_25_4           = [0,1.70E8,3.32E8,5.5e8,2.31e9]            #4 bins 62,60,59,42/per bin for z_flag 3,4,33,34,43,44

SFR_itv                = [0,1.12,1.35,1.7,3.0]                     #4 bins 60,67,54/per bin for z_flag 3,4,33,34,43,44
SFR_itv_20_2           = [0,1.30,3.0]                              #2 bins 102,79/per bin for z_flag 3,4,33,34,43,44
SFR_itv_20_3           = [0,1.18,1.6,3.0]                          #3 bins 60,67,54/per bin for z_flag 3,4,33,34,43,44
SFR_itv_20_4           = [0,1.12,1.35,1.7,3.0]                     #4 bins 60,67,54/per bin for z_flag 3,4,33,34,43,44
SFR_itv_23_2           = [0,1.31,3.0]                              #2 bins 102,79/per bin for z_flag 3,4,33,34,43,44
SFR_itv_23_3           = [0,1.18,1.53,3.0]                          #3 bins 60,67,54/per bin for z_flag 3,4,33,34,43,44
SFR_itv_23_4           = [0,1.10,1.32,1.66,3.0]                     #4 bins 60,67,54/per bin for z_flag 3,4,33,34,43,44
SFR_itv_25_2           = [0,1.35,3.0]                              #2 bins 102,79/per bin for z_flag 3,4,33,34,43,44
SFR_itv_25_3           = [0,1.18,1.6,3.0]                          #3 bins 60,67,54/per bin for z_flag 3,4,33,34,43,44
SFR_itv_25_4           = [0,1.12,1.35,1.7,3.0]                     #4 bins 60,67,54/per bin for z_flag 3,4,33,34,43,44


sSFR_itv               = [-13,-8.55,-8.34,-8.06,-6]                #4 bins 45,44,46,45/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_20_2          = [-13,-8.35,-6]                            #2 bins 90,91/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_20_3          = [-13,-8.50,-8.14,-6]                      #3 bins 61,60,60/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_20_4          = [-13,-8.55,-8.34,-8.06,-6]                 #4 bins 45,44,46,45/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_23_2          = [-13,-8.33,-6]                            #2 bins 90,91/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_23_3          = [-13,-8.50,-8.18,-6]                      #3 bins 61,60,60/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_23_4          = [-13,-8.57,-8.34,-8.06,-6]                 #4 bins 45,44,46,45/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_25_2          = [-13,-8.35,-6]                            #2 bins 90,91/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_25_3          = [-13,-8.50,-8.14,-6]                      #3 bins 61,60,60/per bin for z_flag 3,4,33,34,43,44
sSFR_itv_25_4          = [-13,-8.55,-8.34,-8.06,-6]                 #4 bins 45,44,46,45/per bin for z_flag 3,4,33,34,43,44

Lnuv_itv               = [8,10.01,10.19,10.50,14]                  #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_20_2          = [8,10.19,12]                              #2 bins 89,92/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_20_3          = [8,10.01,10.42,14]                        #3 bins 61,61,59/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_20_4          = [8,9.94,10.21,10.55,14]                   #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_23_2          = [8,10.20,12]                              #2 bins 89,92/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_23_3          = [8,10.01,10.40,14]                        #3 bins 61,61,59/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_23_4          = [8,9.93,10.20,10.53,14]                   #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_25_2          = [8,10.19,12]                              #2 bins 89,92/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_25_3          = [8,10.01,10.42,14]                        #3 bins 61,61,59/per bin for z_flag 3,4,33,34,43,44
Lnuv_itv_25_4          = [8,9.94,10.21,10.55,14]                   #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44

magi_itv               = [22,10.01,10.19,10.50,26]                 #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
magi_itv_20_2          = [22,24.45,26]                             #2 bins 89,92/per bin for z_flag 3,4,33,34,43,44
magi_itv_20_3          = [22,24.23,24.67,26]                       #3 bins 61,60,60/per bin for z_flag 3,4,33,34,43,44
magi_itv_20_4          = [22,24.10,24.45,24.73,26]                 #4 bins 45,44,44,48/per bin for z_flag 3,4,33,34,43,44
magi_itv_23_2          = [22,24.45,26]                             #2 bins 89,92/per bin for z_flag 3,4,33,34,43,44
magi_itv_23_3          = [22,24.22,24.66,26]                       #3 bins 61,60,60/per bin for z_flag 3,4,33,34,43,44
magi_itv_23_4          = [22,24.10,24.45,24.75,26]                 #4 bins 45,44,44,48/per bin for z_flag 3,4,33,34,43,44
magi_itv_25_2          = [22,24.45,26]                             #2 bins 89,92/per bin for z_flag 3,4,33,34,43,44
magi_itv_25_3          = [22,24.23,24.67,26]                       #3 bins 61,60,60/per bin for z_flag 3,4,33,34,43,44
magi_itv_25_4          = [22,24.10,24.45,24.75,26]                 #4 bins 45,44,44,48/per bin for z_flag 3,4,33,34,43,44

phi_itv                = [0,45,90]                                 #2 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
phi_itv_23_2           = [0,45,90]                                 #2 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
phi_itv_23_3           = [0,30,60,90]                              #3 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44

icl_itv                = [0,0.22,0.535,0.8725,1]                   #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
icl_itv_23_2           = [0.22,0.47,1]                             #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
icl_itv_23_3           = [0.22,0.535,0.8725,1]                     #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
icl_itv_23_4           = [0,0.22,0.535,0.8725,1]                   #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44

n_srs_itv              = [0,2.5,3]                                 #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
n_srs_itv_23_2         = [0,2.5,3]                                 #4 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44

r_eff_itv              = [0,1,2.3,10]                              #3 bins 45,44,45,47/per bin for z_flag 3,4,33,34,43,44
r_eff_itv_23_2         = [0,1.7,10]                                #2 bins 89,92/per bin for z_flag 3,4,33,34,43,44
r_eff_itv_23_3         = [0,1,2.3,10]                              #3 bins 61,60,60/per bin for z_flag 3,4,33,34,43,44

#Results
stk_dir_res       = home + 'Stack_Results/'
cat_dir_res       = stk_dir_res + CAT_PARENT
par_dir_res       = cat_dir_res + '/PAIRS' 
par_frg_dir       = par_dir_res + '/FRGRD/'
par_bkg_dir       = par_dir_res + '/BKGRD/'

tbl_dir_res       = cat_dir_res + '/TABLES'
frg_dir_res       = tbl_dir_res + '/FRGRD/'
bkg_dir_res       = tbl_dir_res + '/BKGRD/'
stt_dir_res       = tbl_dir_res + '/STATS/'
irf_dir_res       = tbl_dir_res + '/IRAF/'
std_dir_res       = tbl_dir_res + '/Literature/'
plt_tbl_res       = tbl_dir_res + '/PLOTS/'

plt_dir_res       = cat_dir_res + '/PLOTS' 
ind_plt_res       = plt_dir_res + '/IND/' 
frg_ind_plt       = ind_plt_res + '/FRGRD/' 
bkg_ind_plt       = ind_plt_res + '/BKGRD/' 
res_plt_res       = plt_dir_res + '/RESULTS/' 
ewr_plt_res       = plt_dir_res + '/EW/' 
fit_plt_res       = plt_dir_res + '/FIT/' 

spc_stk_res       = cat_dir_res + '/STACKS'
res_stk_res       = spc_stk_res + '/RESULTS/'
lst_stk_res       = spc_stk_res + '/LAST-FITS'
ind_stk_res       = lst_stk_res + '/RESULTS/'

bst_dir_res       = cat_dir_res + '/BOOTSTRAP'
tbl_bst_res       = bst_dir_res + '/TABLES/'
str_bst_tbl       = tbl_bst_res + '/STRAPS/'
stt_bst_tbl       = tbl_bst_res + '/STATS/'

stk_bst_res       = bst_dir_res + '/STACKS'
stt_bst_stk       = stk_bst_res + '/STATS-BST/'
str_bst_stk       = stk_bst_res + '/STATS-STR/'

lst_bst_res       = stk_bst_res + '/LAST-FITS'
ind_bst_lst       = lst_bst_res + '/STATS-BST/' 
fts_bst_lst       = lst_bst_res + '/STATS-STR/'

plt_bst_res       = bst_dir_res + '/PLOTS'
res_bst_plt       = plt_bst_res + '/RESULTS/'
ind_bst_plt       = plt_bst_res + '/INDIVIDUAL-BS/'

#Output  Tables
ind_tbl           = True
grl_tbl           = True

grl_tbl_nmB       = 'P_Bg_'+ CAT_PARENT 
grl_tbl_nmF       = 'P_Fg_'+ CAT_PARENT 

unq_tbl_opt       = True
hed_un_opt_F      = 'id_F'                                #header of column for uniqueness
hed_un_opt_B      = 'id_B'                                #header of column for uniqueness
grl_tbl_nmB_U     = 'P_Bg_'+ CAT_PARENT +'_U'
grl_tbl_nmF_U     = 'P_Fg_'+ CAT_PARENT +'_U'

verbose           = False
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
		print 'Some of essential the directories does not exist.'
		print 'Check input directories: '
		print str(cat_tbl_chk)+str(tbl_ext_ipt)
		print DIR_SPC_IPT[0]
		print DIR_SPC_IPT[1]
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
