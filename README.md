# VSAT-1D
The Valparaíso Stacking Analysis Tool (VSAT) provides a series of tools for selecting, stacking and anlysing 1D spectra. It is intended for stacking samples of sepctra belonging to large extragalactic catalogues by selecting subsamples of galaxies defined by their available properties (_e.g. redshift, stellar mass, star formatiion rate_) being possible to generate diverse (_e.g. median, average, weighted average, histogram_) composite spectra. However it is possible to also use VSAT on smaller datasets containing any type of astronomical object.

![Alt text](./Images/bootstrap.jpg?raw=true "Stacked spectra computed through median values including CIs.")

## Content

1. Fnc_Stk_Dir.py:
   - Location of the input catalogue and spectral data. 
   - Parameters for selecting subsamples of galaxies according to their physical properties.
   - Stacking and Boootstrap parameters.
   - Location of the resulting products of the stacking analyses _e.g. tables, plots, pre-stacking processed spectra and stacked spectra_.

2. Fnc_Stk_Mth.py:
   - Math functions (e.g. cosmological constants, gaussian, lorentzian and voigt profiles for line emmision/absorption fitting) needed throughout the stacking analysis.

3. Fnc_Stk_Tbl.py 
   - Functions to read, write and modify different tables. 

4. Fnc_Stk_Fts.py 
   - Funtions to access and modify (add, modify, delete) fits headers.

5. Fnc_Stk_Plt.py
   - Plot templates used throughout the stacking analysis. 

6. Fnc_Stk_Spc.py 
   - Functions (_e.g. smoothing, continuum fitting, masking_) for pre-processing the spectra prior to the stacking analysis

7. Fnc_Stk_Stk.py 
   - Core of the stacking tool.
   - Bootstrap function to compute the CIs of the stacked spectra. 
   - Function to select different galaxy subsamples for stacking.
   - Fitting tool for the different emmision/absorption lines through a simple/multiple component gaussian profiles.

8. Fnc_Stk_Utl.py 
   - Auxiliary functions for the stacking analysis.

9. Lines_Dictionary.py
   - Contains a list of identified emmission and absorption lines.

## Parameters
It is possible to pre-processes the spectra before combining them to generate a composite spectrum. This includes continuum substraction/normalization, smoothing, line masking and wavelength shift. The final composite spectra can be then re-processed (_i.e. continuum fitting or smoothing_). Param.py file contains all the parameters needed in each step of the stacking analysis and used in the Example section. 

###### "Pre-Processing Continuum"
**pre_continuum** (_bool, optional_) enables the continuum fitting prior to the stacking, 
**pre_cont_typ** sets the continuum fitting type: ```fit```,```ratio``` or ```difference```, **pre_cont_funct** sets the fitting function: ```legendre```, ```chebyshev```, ```spline1``` or ```spline3``` and **pre_cont_order**  sets the polynomial order. Two fits files will be generated, the continuum normalized/substracted spectra (```spectra-c.fits```) and the continuum fitted (```spectra-c-f.fits```) .

###### "Pre-Processing Smoothing"
**pre_smooth** (_bool, optional_) enables the spectral smoothing, **pre_smooth_shape** selects the smothing kernel (_i.e. gaussian,boxcar,mexican_) and **pre_smooth_size** sets the kernel size in pixel units.

###### "Pre-Processing MASKING"
**pre_mask**  (_bool, optional_) enables spectra masking (after smooothing), **pre_msk_abs_lines** (_bool, optional_) enables line masking from a list of lines included in Lines_Dictionary.py, **pre_mask_type** sets the replacement value for masking: ```NaN``` for numpy NaN, ```constant```  for a constant value or ```continuum``` to use continuum values extracted from the previously continuum fit spectrum, **pre_mask_cte_val**  sets the constant value if ```constant``` is selected, **pre_mask_lw** sets the width around the line center for line masking and **pre_mask_regn** (_bool, optional_) enables masking of a region of the sepctrum delimited by **pre_mask_regn_int** and **pre_mask_regn_fnl**. If ```rshft_corr_direct=True```, then the locatioon of the masks will be corrected by a factor defined by ```rshft_corr```.

###### "Sigma-Clip"
**sigma_clipping** (_bool, optional_) enables sigma cliipping for stacking, **sigma_cut** sets the _n-sigma_ parameter for clipping, **sigma_cen_fct** sets the central function for clipping: ```mean```or ```median``` and **sigma_msk_fill_val** sets the substitute value for clipped value: 
```np.nan``` or ```value```.

###### #Weighting"
**weight_type** 'cont-flux-med' sets the weight type to generate the average weigthing stacked spectra, (_e.g._ ```cont-flux-sum```,```cont-flux-med```, ```cont-flux-avg``` or ```None```), **weight_cnt_flux_get** (_bool_) will get the continuum in a region delimeted by **weight_cnt_flux_lmb_0** and **weight_cnt_flux_lmb_n**. If **weight_type** == ```None``` wights will be set to unity.

###### "Noise Files"
**spectra_noise** (_bool, optional_) includes noise files in the stacking analysis.

###### "Stacks Post Processing "
**post_continuum** and **post_smooth** (_bool, optional_) their parameters are similar to ***pre_continuum***  and ***pre_smooth*** parameters for spectra pre-processing.

![Alt text](./Images/step.jpg?raw=true "Pre-processing of stacked spetra.")
## Lines Dictionary
   - Contains a list of identified emmission and absorption lines in the UV spectral range [800-3000]AA from [IAC's OTELO spectral line summary](http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php), [Shapley+03](https://ui.adsabs.harvard.edu/abs/2003ApJ...588...65S/abstract), [Halliday+08](https://ui.adsabs.harvard.edu/abs/2008A%26A...479..417H/abstract) and [Le Fevre+15](https://ui.adsabs.harvard.edu/abs/2015A%26A...576A..79L/abstract). However this can de modified. The list contains different parameters for the Stacking procedure:
     - 0: Central wavelength.
     - 1: Line width.
     - 2: Line width region for plotting purposes.
     - 3: Line id for stacking routines. 
     - 4: Line id for labelling plots.
     - 5: Line id for fits headers.
     - 6: Line marker for plotting purposes.
     - 7: Line fitting: Central wavelength constrains.
     - 8: Line fitting: Central wavelength offset.
     - 9: Line id for stacking routines.
     - 10: Line fitting: Amplitude bound constrains.

## Plots

## Line fitting
VSAT uses lmfit for line fitting and by default a simple gaussian, however it is possible to use ```GM``` mode to fit multiple gaussian components for emmission before and after the central wavelength of some lines as CIV. The line fitting procedure is performed under a defined wavelength range which can contain a single line or multiple lines as defined in the Line_dictionary.py file. 

![Alt text](./Images/FitSingle.jpg?raw=true "Pre-processing of stacked spetra.")


![Alt text](./Images/FitMultiple.jpg?raw=true "Pre-processing of stacked spetra.")
## Example

Download the Example.zip file which include a catalogue and spectra. Then by simple running ```python Example.py``` will complete all the following steps below.

###### "Stacking"
The following snippet will stack galaxies from the COSMOS field. 

First we define a subsample of galaxies according to some redshift and separation constrains from the ```P_Fg_COSMOS_0-40.csv```  table. 

```python
stamps_subsample_redshift_flag_fg     = Select_Subsamples(op_tbl_F,'redshift_bk_flag',z_flag_itv_bg, test_fg = False, test_bg = False, slc_int = False)
stamps_subsample_redshift_flag_fg_bg  = Select_Subsamples(stamps_subsample_redshift_flag_fg[0][-1],'redshift_fg_flag',z_flag_itv_fg, test_fg = False, 
					test_bg = False, slc_int = False)
fg_in                                 = str(stamps_subsample_redshift_flag_fg_bg[0][-1])
stamps_subsample_sep_fg               = Select_Subsamples(fg_in,'sep_as'     ,SEP_as_itv_23,z_flag_itv_fg, test_fg = False, test_bg = False)

```
This will generate the following tables saved in the tables directory (```~/Example/Stack_Results/COSMOS/TABLES/FRGRD```):
- P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23.csv
- P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44.csv
- P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3.csv
- P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-34.csv
- P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-4.csv
- P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-43.csv
- P_Fg_COSMOS_0-40-ss-zf_B-3-44.csv
- P_Fg_COSMOS_0-40-ss-zf_B-3.csv
- P_Fg_COSMOS_0-40-ss-zf_B-33.csv
- P_Fg_COSMOS_0-40-ss-zf_B-4.csv
- P_Fg_COSMOS_0-40-ss-zf_B-43.csv
- 
Similarly ```Select_Subsamples``` can be used to define subsamples of objects restricted by any given object property. 

Next we stack the subsample. 

```python
f = np.array(Stack_Subsample(stamps_subsample_sep_fg ,
	sel_pre_shf     = selec_spec_shift   ,sel_pre_cnt     = selec_spec_contn      ,sel_pre_msk     = selec_spec_masks      ,
	pre_cnt         = pre_continuum      ,pre_cnt_typ     = pre_cont_typ          ,pre_cnt_lns     = pre_cont_lines        ,
	pre_cnt_fnc     = pre_cont_funct     ,pre_cnt_ord     = pre_cont_order        ,pre_cnt_ovr     = pre_cont_override     ,
	pre_cnt_rpl     = pre_cont_replace   ,pre_cnt_lrj     = pre_cont_low_rej      ,pre_cnt_hrj     = pre_cont_high_rej     ,
	smt_spc_pre     = pre_smooth         ,smt_shp_pre     = pre_smooth_shape      ,smt_sze_pre     = pre_smooth_size       ,
	pre_msk         = False              ,
	pre_msk_typ     = pre_mask_type      ,pre_msk_abs_lne = False                 ,pre_msk_cte_val = pre_mask_cte_val      ,
	pre_msk_blu_rgn = False              ,pre_lmb_min     = pre_mask_regn_int     ,pre_lmb_max     = pre_mask_regn_fnl,
	sig_clp         = sigma_clipping     ,
	sig_cut         = sigma_cut          ,sig_fct         = sigma_cen_fct         ,sig_fll         = sigma_msk_fill_val    ,
	wgt_typ         = weight_type        ,
	get_cont_flux   = weight_cnt_flux_get,gcv_lmbd_i      = weight_cnt_flux_lmb_0 ,gcv_lmbd_f      = weight_cnt_flux_lmb_n ,
	wrt_fits        = True               ,spc_nse         = spectra_noise         ,
	pst_cnt         = post_continuum     ,pst_cnt_typ     = post_cont_typ         ,pst_cnt_lns     = post_cont_lines       ,
	pst_cnt_fnc     = post_cont_funct    ,pst_cnt_ord     = post_cont_order       ,pst_cnt_ovr     = post_cont_override    ,
	pst_cnt_rpl     = post_cont_replace  ,pst_cnt_lrj     = post_cont_low_rej     ,pst_cnt_hrj     = post_cont_high_rej    ,
	smt_spc_pst     = post_smooth        ,smt_shp_pst     = post_smooth_shape     ,smt_sze_pst     = post_smooth_size))
```
This will generate the following fits files in the results directory (```~/Example/Stack_Results/COSMOS/STACKS/RESULTS/```):

```
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-sum.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-avg.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-med.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-hst.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-std.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-rms.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-hsw.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-wsu.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-suw.fits
 - P_Fg_COSMOS_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-0-23-stk-avw.fits
```

It will also generate ```*-c.fits```and ```*-smt.fits```files if ```pst_cnt``` and ```smt_spc_pst ``` are set to ```True```.

Auxiliary files generated during the stacking process are saved in ```~/Example/Stack_Results/COSMOS/STACKS/LAST-FITS/RESULTS```
and can be deleted. These files include a copy of the original spectra, continuum corrected, smoothed and interpolated spectra, continuum fit log files and .tex files with the spectra used in the stacking process. Notice that a same file can be used in differeent stackings and hence the -int.fits files include a number of the interpolated version.

###### "Stats"
Statistical values from the stacked galaxies can be obtained and saved in tables (```~/Example/Stack_Results/COSMOS/TABLES/STATS```) through:

```python
[stats_table(tblnm,tbl_format_opt) for tblnm in stamps_subsample_sep_fg[0]]
```

###### "Plots"
To plot the generated composite spectra (_e.g. average, median weighted average_)including post-processed versions (_continuum normalized and smoothed_) can bee generated with:

```python
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
			plt_stk_med    = True             ,plt_stk_avg      = True ,plt_stk_avw = True)
```
This will generate and save a pdf plot file (```~/Example/Stack_Results/COSMOS/PLOTS/RESULTS/```) including the median, average, and weighted average in the upper panel and an histogram of the number of spectra combined per wavelength element in the lower panel.

![Alt text](./Images/Stacked.jpg?raw=true "Stacked spectra computed COSMOS field.")

**plt_cnt_stk_spc** = ```True``` generates a plot of all the individual spectra used to generate the composite spectra (upper panel) with the composite specra in the bottom panel.

![Alt text](./Images/Stacked-Contribution.jpg?raw=true "Stacked spectra COSMOS field.")

**plt_ind_spec** = ```True``` will generate individual plots (```~/Example/Stack_Results/COSMOS/PLOTS/IND//FRGRD/0-23/```) of all the spectra used to generate the compsite spectra. 

![Alt text](./Images/Spec-Individual.jpg?raw=true "Stacked spectra COSMOS field.")

It will also generate a detailed plot of every step of the pre-processing procedure prior to the stacking process of all the combined spectra.
![Alt text](./Images/Spec-Step.jpg?raw=true "Stacked spectra COSMOS field.")

###### "Line Fitting"

```python
Plot_Idp_Spc_Lne(
		int_typ_spl     = 'sep_as'        ,stk_function   = 'avg-c-smt',
		lmb_min         = 1200            ,lmb_max        = 1700 ,
		fit_type        = 'lmfit'         ,fit_fnct       = 'gauss',
		verbose         = True            ,autoaxis       = True ,
		pre_off_plt     = False           ,ofs_ctr_fit    = False ,
		n_int_spt       = 4               ,
		lower_shift     = 4               ,upper_shift    = 0,
		max_sep         = 23              ,
		mlt_stk_fct     = 'avg'           ,
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
```

All the line initial and fitted parameters (line center, amplitude, line width) will be saved in the corresponding composite fits header. For example for Lyman alpha the following headers will be created:

 - L05_CGLC=    1207.832434368174 / Lya1215.7 Ctr  1GF Crctlmfit
 - L05_AGLC=  -0.4166563657612596 / Lya1215.7 Amp  1GF Crctlmfit
 - L05_SGLC=    3.475333003707025 / Lya1215.7 Sigm 1GF Crctlmfit
 - L05_FGLC=    8.183783820286921 / Lya1215.7 FWHM 1GF Crctlmfit
 - L05_WGLC=         3.6296469195 / Lya1215.7 EW   1GF Crctlmfit                   

L05 corresponds to the line identifier for fits headers purpooses and defined in the Line_Dictionarry.py file.

![Alt text](./Images/LINE-FIT-COSMOS-avg-c-smt-G-Ind-Splt.jpg?raw=true "Stacked spectra COSMOS field.")	
```python
Plot_Slc_Spc_Lne(
		int_typ_spl    = 'sep_as'       ,stk_function   = 'avg-c-smt',
		lmb_min        = 1200           ,lmb_max        = 1900,
		stk_fct        = ['avg'],
		fit_type       = 'lmfit'        ,fit_fnct       = 'gauss',
		verbose        = True           ,autoaxis       = True,
		pre_off_plt    = False          ,n_int_spt      = 4 ,
		lower_shift     = 4             ,upper_shift    = 0 ,   #Foregorund
		plt_ind_fit    = True           ,
		autoaxis_SSL   = True           ,
		lbl_col_idv    = True           ,
		fpt_foreground = True           ,fpt_background = False,
		max_sep        = 23             ,
		empty_plots    = 2              ,landscape_plt  = True,
		splt_ind_lns   = True)
```
![Alt text](./Images/LINE-FIT-COSMOS-avg-c-smt-G-Mlt-Splt.jpg?raw=true "Stacked spectra COSMOS field.")	

###### "Bootstrap"
To compute the Confidence Inteervals (CIs) of the composite spectra it is possible to bootstrap the spectra used in the stacking process. 

**bs_iteration_num** defines the number of bootstrap repetitionsm, **bs_percentage** defines the percentaje of elements to be resampled in each iteration (_default=1_), **bs_function**  defines the function used for the bootstrap.

It is possible to complete a broken BS process. 
**comp_BS_run**	(_bool_) enables the BS completion
**bst_brk_int** sets the reinizilation step at which the process broke, **bst_brk_lst** sets the BS number repetitions (_default=bs_iteration_num_)  #Last iteration step (default = bs_iteration_num) and **crts_BS_otb** (_bool_) creates the BS master table only for Completed BS Repetitions(_= bs_iteration_num_) but without master table. 


To generate the CIs, first we define:

```python
prefixF       = 'P_Fg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
tables        = [frg_dir_res + prefixF + '0-23.csv']
bs_function_s = [ '_med-c-smt']
```

Then we run the BS process:

```
for tbl_2b_btstr in tables:
	Boot_out = Bootstrap(tbl_2b_btstr,None,bs_iteration_num,bs_percentage,
				cmp_bst_cyc     = comp_BS_run          ,int_bst_brk     = bst_brk_int           ,lst_bst_brk     = bst_brk_lst           ,
				bst_cyc_otb     = crts_BS_otb          ,
				sel_pre_shf     = selec_spec_shift     ,sel_pre_cnt     = selec_spec_contn      ,sel_pre_msk     = selec_spec_masks      ,
				pre_cnt         = pre_continuum        ,pre_cnt_typ     = pre_cont_typ          ,pre_cnt_lns     = pre_cont_lines        ,
				pre_cnt_fnc     = pre_cont_funct       ,pre_cnt_ord     = pre_cont_order        ,pre_cnt_ovr     = pre_cont_override     ,
				pre_cnt_rpl     = pre_cont_replace     ,pre_cnt_lrj     = pre_cont_low_rej      ,pre_cnt_hrj     = pre_cont_high_rej     ,
				smt_spc_pre     = pre_smooth           ,smt_shp_pre     = pre_smooth_shape      ,smt_sze_pre     = pre_smooth_size       ,
				pre_msk         = False                ,pre_msk_typ     = pre_mask_type         ,pre_msk_abs_lne = False                 ,
				pre_msk_cte_val = pre_mask_cte_val     ,
				pre_msk_rgn     = False                ,pre_lmb_min.    = pre_mask_regn_int     ,pre_lmb_max     = pre_mask_regn_fnl,
				sig_clp         = sigma_clipping       ,sig_cut         = sigma_cut             ,sig_fct         = sigma_cen_fct         ,
				sig_fll         = sigma_msk_fill_val   ,
				wgt_typ         = weight_type          ,
				get_cont_flux   = weight_cnt_flux_get  ,gcv_lmbd_i      = weight_cnt_flux_lmb_0 ,
				gcv_lmbd_f      = weight_cnt_flux_lmb_n,
				wrt_fits        = True                 ,spc_nse         = spectra_noise         ,
				pst_cnt         = post_continuum       ,pst_cnt_typ     = post_cont_typ         ,pst_cnt_lns     = post_cont_lines       ,
				pst_cnt_fnc     = post_cont_funct      ,pst_cnt_ord     = post_cont_order       ,pst_cnt_ovr     = post_cont_override    ,
				pst_cnt_rpl     = post_cont_replace    ,pst_cnt_lrj     = post_cont_low_rej     ,pst_cnt_hrj     = post_cont_high_rej    ,
				smt_spc_pst     = post_smooth          ,smt_shp_pst     = post_smooth_shape     ,smt_sze_pst     = post_smooth_size)
	for bs_function in bs_function_s:
		Boot_out = str_bst_tbl + (tbl_2b_btstr.split('/')[-1]).split('.csv')[0] + '-BS_MST_'+ str(bs_iteration_num)+'.csv'
		print
		print 'Stacking all straps:'
		print 'From table:',Boot_out
		print 'Function: ',bs_function
		stamps_bootstrap_1 = Select_Subsamples(Boot_out,None,None,test_fg = True, test_bg = False, slc_int = False, slc_smp = False,bs_func = bs_function)
		stacks_bootstrap_1 = np.array(Stack_Subsample(stamps_bootstrap_1,bs_func = bs_function,
					sel_pre_shf     = selec_spec_shift     ,sel_pre_cnt     = selec_spec_contn      ,sel_pre_msk     = selec_spec_masks      ,
					pre_cnt         = pre_continuum        ,pre_cnt_typ     = pre_cont_typ          ,pre_cnt_lns     = pre_cont_lines        ,
					pre_cnt_fnc     = pre_cont_funct       ,pre_cnt_ord     = pre_cont_order        ,pre_cnt_ovr     = pre_cont_override     ,
					pre_cnt_rpl     = pre_cont_replace     ,pre_cnt_lrj     = pre_cont_low_rej      ,pre_cnt_hrj     = pre_cont_high_rej     ,
					smt_spc_pre     = pre_smooth           ,smt_shp_pre     = pre_smooth_shape      ,smt_sze_pre     = pre_smooth_size       ,
					pre_msk         = False                ,pre_msk_typ     = pre_mask_type         ,
					pre_msk_abs_lne = False                ,pre_msk_cte_val = pre_mask_cte_val      ,
					pre_msk_rgn     = False                ,pre_lmb_min     = pre_mask_regn_int     ,pre_lmb_max     = pre_mask_regn_fnl,
					sig_clp         = sigma_clipping       ,sig_cut         = sigma_cut             ,sig_fct         = sigma_cen_fct         ,
					sig_fll         = sigma_msk_fill_val   ,
					wgt_typ         = weight_type          ,
					get_cont_flux   = weight_cnt_flux_get  ,gcv_lmbd_i      = weight_cnt_flux_lmb_0 ,
					gcv_lmbd_f      = weight_cnt_flux_lmb_n,
					wrt_fits        = True                 ,spc_nse         = spectra_noise         ,
					pst_cnt         = post_continuum       ,pst_cnt_typ     = post_cont_typ         ,pst_cnt_lns     = post_cont_lines       ,
					pst_cnt_fnc     = post_cont_funct      ,pst_cnt_ord     = post_cont_order       ,pst_cnt_ovr     = post_cont_override    ,
					pst_cnt_rpl     = post_cont_replace    ,pst_cnt_lrj     = post_cont_low_rej     ,pst_cnt_hrj     = post_cont_high_rej    ,
					smt_spc_pst     = post_smooth          ,smt_shp_pst     = post_smooth_shape     ,smt_sze_pst     = post_smooth_size,
					stk_pct_mde     = True                 ,stk_wgt_mde     = False))

```

This process will create into the ```~/BOOTSTRAP/STACKS/``` directory three diferent subdirectories: 
 - ```~/BOOTSTRAP/STACKS/LAST-FITS``` contains all the files used to generate stacked boootstrap in each repetition with a similar structure similar to  ```~/Example/Stack_Results/COSMOS/STACKS/RESULTS/``` directory and described above in the **Stacking** section.
 - ```~/BOOTSTRAP/STACKS/STATS-STR``` contains all the stacked boootstrap repetitions (_e.g. FILE_NAME-BS-1-stk-avg.fits, FILE_NAME-BS-2-stk-avg.fits, ... FILE_NAME-BS-N-stk-avg.fits_) stacked boootstrap repetitions.
 - ```~/BOOTSTRAP/STACKS/STATS-BST``` contains the bootsrap stacked spectra considering all the _N_ stacked boootstrap repetitions.

Finally, to plot the Stacked spectra including the CIs generated thorugh the BS process; first we define:

```
prefixF 	= 'P_Fg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
separation      = ['0-23'] 			#sep_as
bs_function_s   = ['med-c-smt']
```

and then

```
for element in itertools.product(separation,bs_function_s):	
	sep         = element[0]
	bs_function = element[1]

	print
	print 'Generating BS plot (',str(bs_iteration_num),') for: sep: ',sep,', function: ',bs_function

	fits_file     = [
	                 res_stk_res + prefix + sep+'-stk-'+bs_function +'.fits',
	                 res_stk_res + prefix + sep+'-stk-hst.fits'
	                 ]                 
	fits_file_err = [
	                 stt_bst_stk + prefix + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-med-c-smt.fits', 
	                 stt_bst_stk + prefix + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-1sl-c-smt.fits',
	                 stt_bst_stk + prefix + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-1sh-c-smt.fits',
	                 stt_bst_stk + prefix + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-2sl-c-smt.fits',
	                 stt_bst_stk + prefix + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-2sh-c-smt.fits',
	                 stt_bst_stk + prefix + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-3sl-c-smt.fits',
	                 stt_bst_stk + prefix + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-3sh-c-smt.fits',
	                 stt_bst_stk + prefix + sep+'-BS_MST_' + str(bs_iteration_num) + '_' + bs_function +'-stk-hst.fits'
	                 ]
	Plot_Idp_Spc_BS(fits_file,fits_file_err,bs_iteration_num,
					nsigma               = 2,
					min_x_lim_Idp        = 1150 , max_x_lim_Idp  = 1900,
					autoaxis_Idp         = False, aaxs_Idp_ml_y  = True, 
					min_y_lim_Idp        = 0.5  , max_y_lim_Idp  = 2.0,
					sep_lin_min          = 10 )
```
![Alt text](./Images/P_Fg_COSMOS_BS_MST_100_med-c-smt-1150-1900.jpg?raw=true "Stacked spectra COSMOS field.")	

This plot will be saved in ```~/Example/Stack_Results/COSMOS/BOOTSTRAP/PLOTS/RESULTS/```.
**Notice that this plots is only useful to visualize the distribution of spectra generated through the bootstrap.** To properly generate CIs of the emission/absorption EW lines measured above, all the bootstrapped spectra should be fitted to obtain a distribution of EW. This can be easily done with:

```
import Lines_Dictionary
LINES = Lines_Dictionary.LINES_PLT_BG

prefixF  = 'P_Fg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
prefix   = prefixF
function = 'med-c-smt'

bs_iteration_num = 2
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
		lmb_max      = 1900,
		plt_fit      = True,
		verbose      = True,
		stk_function = function,
		fit_type     = 'lmfit',
		fit_fnct     = 'gauss',
		pre_off_plt  = True,
		org_spc_fle  = spec_test_0,
		ivl_fts_hdr  = True)
pb.finish()
```
If ```ivl_fts_hdr=True``` then the initial guess values will be read from the fits files headers, otherwise they will be set by the defualt line quantities defined in the Line_Dictionary.py file .
This will generate individual plots for each line profile fitted, located in ```~/Example/Stack_Results/COSMOS/BOOTSTRAP/PLOTS/INDIVIDUAL-BS/```.

As explained above all the fitted parameters (_i.e. line properties (amplitude, sigma, line center), fitted values_) are saved in the corresponding compoosite spectra.
## Dependencies
Currently VSAT works only with astropy 2.0 as it relies on pyraf continuum task for continuum normalization. However a new version will be released dropping this dependency.
 - [astropy](https://www.astropy.org)
 - [bottleneck](https://pypi.org/project/Bottleneck/)
 - [pandas](https://pandas.pydata.org)
 - [scipy](https://www.scipy.org)
 - [numpy](https://numpy.org)
 - [lmfit](https://lmfit.github.io/lmfit-py/)
 - [matplotlib](https://matplotlib.org)
 - [pysynphot](https://pysynphot.readthedocs.io/en/latest/)
 - [pyraf](https://astroconda.readthedocs.io/en/latest/installation.html)
 - [termcolor](https://pypi.org/project/termcolor/)
## License




