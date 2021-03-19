from Fnc_Stk_Fts import *
from Fnc_Stk_Spc import *
from Fnc_Stk_Mth import *
from Fnc_Stk_Dir import *
from Fnc_Stk_Utl import *
from Fnc_Stk_Tbl import *
from Fnc_Stk_Stk import *
from Fnc_Stk_Plt import *

os.system('clear')

##################################################STACK###################################################

prefixB = 'P_Bg_' + CAT_PARENT + '_0-40-PRP-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
prefixF = 'P_Fg_' + CAT_PARENT + '_0-40-PRP-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'


stamps_subsample_redshift_flag_fg     = Select_Subsamples(op_tbl_F,'redshift_bk_flag',z_flag_itv_bg, test_fg = False, test_bg = False, slc_int = False)
stamps_subsample_redshift_flag_fg_bg  = Select_Subsamples(stamps_subsample_redshift_flag_fg[0][-1],'redshift_fg_flag',z_flag_itv_fg, test_fg = False, 
										test_bg = False, slc_int = False)
fg_in                                 = str(stamps_subsample_redshift_flag_fg_bg[0][-1])
stamps_subsample_sep_fg               = Select_Subsamples(fg_in,'sep_as'     ,SEP_as_itv_23,z_flag_itv_fg, test_fg = False, test_bg = False)
f                                     = np.array(Stack_Subsample(stamps_subsample_sep_fg      ,
										sel_pre_shf     = selec_spec_shift   ,sel_pre_cnt     = selec_spec_contn      ,sel_pre_msk     = selec_spec_masks      ,
										pre_cnt         = pre_continuum      ,pre_cnt_typ     = pre_cont_typ          ,pre_cnt_lns     = pre_cont_lines        ,
										pre_cnt_fnc     = pre_cont_funct     ,pre_cnt_ord     = pre_cont_order        ,pre_cnt_ovr     = pre_cont_override     ,
										pre_cnt_rpl     = pre_cont_replace   ,pre_cnt_lrj     = pre_cont_low_rej      ,pre_cnt_hrj     = pre_cont_high_rej     ,
										smt_spc_pre     = pre_smooth         ,smt_shp_pre     = pre_smooth_shape      ,smt_sze_pre     = pre_smooth_size       ,
										pre_msk         = False              ,
										pre_msk_typ     = pre_mask_type      ,pre_msk_abs_lne = False                 ,pre_msk_cte_val = pre_mask_cte_val      ,
										pre_msk_rgn     = False              ,pre_msk_min     = pre_mask_regn_int     ,pre_msk_max     = pre_mask_regn_fnl     ,
										sig_clp         = sigma_clipping     ,
										sig_cut         = sigma_cut          ,sig_fct         = sigma_cen_fct         ,sig_fll         = sigma_msk_fill_val    ,
										wgt_typ         = weight_type        ,
										get_cont_flux   = weight_cnt_flux_get,gcv_lmbd_i      = weight_cnt_flux_lmb_0 ,gcv_lmbd_f      = weight_cnt_flux_lmb_n ,
										wrt_fits        = True               ,spc_nse         = spectra_noise         ,
										pst_cnt         = post_continuum     ,pst_cnt_typ     = post_cont_typ         ,pst_cnt_lns     = post_cont_lines       ,
										pst_cnt_fnc     = post_cont_funct    ,pst_cnt_ord     = post_cont_order       ,pst_cnt_ovr     = post_cont_override    ,
										pst_cnt_rpl     = post_cont_replace  ,pst_cnt_lrj     = post_cont_low_rej     ,pst_cnt_hrj     = post_cont_high_rej    ,
										smt_spc_pst     = post_smooth        ,smt_shp_pst     = post_smooth_shape     ,smt_sze_pst     = post_smooth_size))
##################################################STACK###################################################

###################################################STAT###################################################
print
[stats_table(tblnm,tbl_format_opt) for tblnm in stamps_subsample_sep_fg[0]]
[stats_table(tblnm,tbl_format_opt) for tblnm in stamps_subsample_sep_fg[0]]
print
###################################################STAT###################################################

###################################################PLOT###################################################
Plot_All_Spec_All_Int(
			frgrnd_plt     = True,
			bkgrnd_plt     = False,
			n_int_spt      = 4                , int_typ_spl     = 'sep_as',
			min_x_lim_Idp  = 1200             , max_x_lim_Idp   = 1900,
			plt_ind_spec   = False            , plt_cnt_stk_spc = False,
			wgt_typ        = weight_type,
			autoaxis_Idp   = False            , aaxs_Idp_ml_y   = True, 
			min_y_lim_Idp  = 0.5              , max_y_lim_Idp   = 2.0,
			lower_shift    = 0                , upper_shift     = 0, 
			only_stt_tbl   = False            , 
			SNR_lines      = None             , show_legends    = True,
			max_sep        = 23               ,
			mlt_stk_fct    = 'avg'            ,
			fpt_foreground = True             ,fpt_background   = False,
			plt_stk_med    = True             ,plt_stk_avg      = True  ,plt_stk_avw = False)
###################################################PLOT###################################################


####################################################FIT###################################################
Plot_Idp_Spc_Lne(
		int_typ_spl     = 'sep_as'        ,stk_function   = 'med-c-smt',
		lmb_min         = 1200            ,lmb_max        = 1900 ,
		fit_type        = 'lmfit'         ,fit_fnct       = 'gauss',
		verbose         = True            ,autoaxis       = True ,
		pre_off_plt     = False           ,ofs_ctr_fit    = False ,
		n_int_spt       = 4               ,
		lower_shift     = 4               ,upper_shift    = 0,
		max_sep         = 23              ,
		mlt_stk_fct     = 'med'           ,
		mke_lne_fit     = True            , 
		fit_vls_hdr     = True            ,
		int_vlf_hdr     = True            ,
		uft_lne_vls     = False           ,
		cnt_bnp_adj     = True            ,
		fpt_foreground  = True            ,fpt_background = False,
		fix_ctr_gau     = False,
		fix_pre_gau     = False,
		fix_pst_gau     = False,
		fix_ctr_gau_1   = False           ,
		fix_ctr_gau_2   = False     	  ,
		pre_shf_lim     = 0               ,pst_shf_lim    = 0,
		pre_shf_ctr     = 0               ,pst_shf_ctr    = 0,
		fix_mdl_gau     = False           ,
		mdl_shf_ctr     = 0               ,mdl_shf_lim    = 0,
		ivl_fts_hdr     = False
					)
####################################################FIT###################################################

##################################################PLOT-FIT################################################
Plot_Slc_Spc_Lne(
		int_typ_spl    = 'sep_as'       ,stk_function   = 'med-c-smt',
		lmb_min        = 1200           ,lmb_max        = 1900,
		stk_fct        = ['med'],
		fit_type       = 'lmfit'        ,fit_fnct       = 'gauss',
		verbose        = True           ,autoaxis       = True,
		pre_off_plt    = False          ,n_int_spt      = 4 ,
		lower_shift     = 4             ,upper_shift    = 0 ,   
		plt_ind_fit    = True           ,
		autoaxis_SSL   = True           ,
		lbl_col_idv    = True           ,
		fpt_foreground = True           ,fpt_background = False,
		max_sep        = 23             ,
		empty_plots    = 2              ,landscape_plt  = True,
		splt_ind_lns   = True)
##################################################PLOT-FIT################################################

#################################################BOOTSTRAP################################################

prefixB       = 'P_Bg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
prefixF       = 'P_Fg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
tables        = [frg_dir_res + prefixF + '0-23.csv']
bs_function_s = ['_med-c-smt']
print
print "\n".join(['Function: ' + function for function in bs_function_s])
print
print "\n".join(['Function: ' + function for function in tables])
print
print bs_iteration_num

for tbl_2b_btstr in tables:
	Boot_out = Bootstrap(tbl_2b_btstr,None,bs_iteration_num,bs_percentage,
				cmp_bst_cyc     = comp_BS_run        ,int_bst_brk     = bst_brk_int           ,lst_bst_brk     = bst_brk_lst,
				bst_cyc_otb     = crts_BS_otb        ,
				sel_pre_shf     = selec_spec_shift   ,sel_pre_cnt     = selec_spec_contn      ,sel_pre_msk     = selec_spec_masks,
				pre_cnt         = pre_continuum      ,pre_cnt_typ     = pre_cont_typ          ,pre_cnt_lns     = pre_cont_lines,
				pre_cnt_fnc     = pre_cont_funct     ,pre_cnt_ord     = pre_cont_order        ,pre_cnt_ovr     = pre_cont_override,
				pre_cnt_rpl     = pre_cont_replace   ,pre_cnt_lrj     = pre_cont_low_rej      ,pre_cnt_hrj     = pre_cont_high_rej,
				smt_spc_pre     = pre_smooth         ,smt_shp_pre     = pre_smooth_shape      ,smt_sze_pre     = pre_smooth_size,
				pre_msk         = False              ,
                pre_msk_typ     = pre_mask_type      ,pre_msk_abs_lne = False                 ,pre_msk_cte_val = pre_mask_cte_val,
				pre_msk_rgn     = False              ,pre_lmb_min     = pre_mask_regn_int     ,pre_lmb_max     = pre_mask_regn_fnl,
				sig_clp         = sigma_clipping     ,
                sig_cut         = sigma_cut          ,sig_fct         = sigma_cen_fct         ,sig_fll         = sigma_msk_fill_val,
				wgt_typ         = weight_type        ,
                get_cont_flux   = weight_cnt_flux_get,gcv_lmbd_i       = weight_cnt_flux_lmb_0 ,gcv_lmbd_f      = weight_cnt_flux_lmb_n,
				wrt_fits        = True               ,spc_nse         = spectra_noise         ,
				pst_cnt         = post_continuum     ,pst_cnt_typ     = post_cont_typ         ,pst_cnt_lns     = post_cont_lines,
				pst_cnt_fnc     = post_cont_funct    ,pst_cnt_ord     = post_cont_order       ,pst_cnt_ovr     = post_cont_override,
				pst_cnt_rpl     = post_cont_replace  ,pst_cnt_lrj     = post_cont_low_rej     ,pst_cnt_hrj     = post_cont_high_rej,
				smt_spc_pst     = post_smooth        ,smt_shp_pst     = post_smooth_shape     ,smt_sze_pst     = post_smooth_size)
	for bs_function in bs_function_s:
		Boot_out = str_bst_tbl + (tbl_2b_btstr.split('/')[-1]).split('.csv')[0] + '-BS_MST_'+ str(bs_iteration_num)+'.csv'
		print
		print 'Stacking all straps:'
		print 'From table:',Boot_out
		print 'Function: ',bs_function
		stamps_bootstrap_1 = Select_Subsamples(Boot_out,None,None,test_fg = True, test_bg = False, slc_int = False, slc_smp = False,bs_func = bs_function)
		stacks_bootstrap_1 = np.array(Stack_Subsample(stamps_bootstrap_1,bs_func = bs_function,
					sel_pre_shf     = selec_spec_shift   ,sel_pre_cnt     = selec_spec_contn      ,sel_pre_msk     = selec_spec_masks ,
					pre_cnt         = pre_continuum      ,pre_cnt_typ     = pre_cont_typ          ,pre_cnt_lns     = pre_cont_lines,
					pre_cnt_fnc     = pre_cont_funct     ,pre_cnt_ord     = pre_cont_order        ,pre_cnt_ovr     = pre_cont_override,
					pre_cnt_rpl     = pre_cont_replace   ,pre_cnt_lrj     = pre_cont_low_rej      ,pre_cnt_hrj     = pre_cont_high_rej,
					smt_spc_pre     = pre_smooth         ,smt_shp_pre     = pre_smooth_shape      ,smt_sze_pre     = pre_smooth_size,
					pre_msk         = False              ,pre_msk_typ     = pre_mask_type         ,
					pre_msk_abs_lne = False              ,pre_msk_cte_val = pre_mask_cte_val      ,
					pre_msk_rgn     = False              ,pre_lmb_min     = pre_mask_regn_int     ,pre_lmb_max     = pre_mask_regn_fnl,
					sig_clp         = sigma_clipping     ,
                    sig_cut         = sigma_cut          ,sig_fct         = sigma_cen_fct         ,sig_fll         = sigma_msk_fill_val,
					wgt_typ         = weight_type        ,
                    get_cont_flux   = weight_cnt_flux_get,gcv_lmbd_i      = weight_cnt_flux_lmb_0 ,gcv_lmbd_f      = weight_cnt_flux_lmb_n,
					wrt_fits        = True               ,spc_nse         = spectra_noise         ,
					pst_cnt         = post_continuum     ,pst_cnt_typ     = post_cont_typ         ,pst_cnt_lns     = post_cont_lines,
					pst_cnt_fnc     = post_cont_funct    ,pst_cnt_ord     = post_cont_order       ,pst_cnt_ovr     = post_cont_override,
					pst_cnt_rpl     = post_cont_replace  ,pst_cnt_lrj     = post_cont_low_rej     ,pst_cnt_hrj     = post_cont_high_rej,
					smt_spc_pst     = post_smooth        ,smt_shp_pst     = post_smooth_shape     ,smt_sze_pst     = post_smooth_size,
					stk_pct_mde     = True               ,stk_wgt_mde     = False))
#################################################BOOTSTRAP################################################

###############################################PLOT-BOOTSTRAP#############################################
prefixF 		= 'P_Fg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
separation      = ['0-23'] 			#sep_as
bs_function_s   = ['med-c-smt']

for element in itpdc(separation,bs_function_s):	
	sep         = element[0]
	bs_function = element[1]

	print
	print 'Generating BS plot (',str(bs_iteration_num),') for: sep: ',sep,', function: ',bs_function

	fits_file     = [
	                 res_stk_res + prefixF + sep+'-stk-'+bs_function +'.fits',
	                 res_stk_res + prefixF + sep+'-stk-hst.fits'
	                 ]                 
	fits_file_err = [
	                 stt_bst_stk + prefixF + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-med-c-smt.fits', 
	                 stt_bst_stk + prefixF + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-1sl-c-smt.fits',
	                 stt_bst_stk + prefixF + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-1sh-c-smt.fits',
	                 stt_bst_stk + prefixF + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-2sl-c-smt.fits',
	                 stt_bst_stk + prefixF + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-2sh-c-smt.fits',
	                 stt_bst_stk + prefixF + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-3sl-c-smt.fits',
	                 stt_bst_stk + prefixF + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-3sh-c-smt.fits',
	                 stt_bst_stk + prefixF + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-hst.fits'
	                 ]
	Plot_Idp_Spc_BS(fits_file,fits_file_err,bs_iteration_num,
					nsigma               = 2,
					min_x_lim_Idp        = 1150 , max_x_lim_Idp  = 1900,
					autoaxis_Idp         = False, aaxs_Idp_ml_y  = True, 
					min_y_lim_Idp        = 0.5  , max_y_lim_Idp  = 2.0,
					sep_lin_min          = 10 )

###############################################PLOT-BOOTSTRAP#############################################


##########################################FIT-BOOTSTRAP-REPETITIONS########################################
import Lines_Dictionary
LINES = Lines_Dictionary.LINES_PLT_BG

prefixF  = 'P_Fg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
prefix   = prefixF
function = 'med-c-smt'

fitting_m        = 'lmfit' 
fitting_f        = 'gauss' 

pb = ProgressBar(bs_iteration_num)
pb.start()
for strap in range(1,bs_iteration_num+1):
	pb.update(strap)
	spec_test   = str_bst_stk + prefix + '0-23' + '-BS-'  + str(strap) + '-stk-' + function +'.fits'		
	spec_test_0 = res_stk_res + prefix + '0-23' + '-stk-' + function   + '.fits'		
	Stk_Fit_Lines(spec_test,
		lmb_min      = 1200,
		lmb_max      = 1220,
		plt_fit      = True,
		verbose      = True,
		stk_function = function,
		fit_type     = 'lmfit',
		fit_fnct     = 'gauss',
		pre_off_plt  = False,
		org_spc_fle  = spec_test_0,
		ivl_fts_hdr  = True)
pb.finish()
##########################################FIT-BOOTSTRAP-REPETITIONS########################################
