# VSAT-1D
Stacking 1D spectra tools

![Alt text](./Images/bootstrap.jpg?raw=true "Stacked spectra computed through median values including CIs.")

## Content

1. Fnc_Stk_Dir.py:
   - Location of the input catalogue and spectral data. 
   - Parameters for selecting galaxy pairs and subsamples of galaxies according to their physical properties. 
   - Location of the resulting products of the stacking analyses e.g. tables, plots, pre stacking processed spectra and stacked spectra.

2. Fnc_Stk_Mth.py:
   - Useful math functions (e.g. cosmological constants, gaussian, lorentzian and voigt profiles for line emmision/absorption fitting).

3. Fnc_Stk_Tbl.py 
   - Functions to access different tables. 

4. Fnc_Stk_Fts.py 
   - Funtions to access and modify (add, modify, delete) fits headers.

5. Fnc_Stk_Plt.py
   - Plot templates used through all the stacking analysis. 

6. Fnc_Stk_Spc.py 
   - Contains a series of functions for pre-processing the spectra needed for the stacking analysis

7. Fnc_Stk_Stk.py 
   - Contains the core of the stacking tool
   - Bootstrap function to compute the CIs of the stacked spectra. 
   - Function to select different galaxy subsamples for stacking.
   - Fitting tool for the different emmision/absorption lines through a simple/multiple component gaussian profiles.

8. Fnc_Stk_Utl.py 
   - Contains auxiliary functions for the stacking analysis


9. Lines_Dictionary.py
   - Contains a list of identified emmission and absorption lines.
## Parameters
It is possible to perform a pre-processing of the spectra before stacking them to create a composite spectrum. This includes continuum substraction/normalization, gaussian smoothing, line masking and wavelength shift. The final composite spectra can be processed to fit the continuum and smooth it. Param.py file contains all the parameters of each stage. 

###### "Pre-Processing Continuum"
   - pre_continuum          = False                                     # Continuum Fitting/Normalization
   - pre_cont_typ           = 'ratio'                                   # Continuum fitting type fit,ratio,difference
   - pre_cont_lines         = '*'                                       # Image lines to be fit
   - pre_cont_funct         = 'spline3'                                 # Fitting function: legendre, chebyshev, spline1, spline3
   - pre_cont_order         = 49                                        # Order Polynomial / num pieces spline
   - pre_cont_override      = 'yes'                                     # Override previous norm spec
   - pre_cont_replace       = 'no'                                      # Replace rejected points by fit?
   - pre_cont_low_rej       = 3                                         # Low rejection in sigma of fit
   - pre_cont_high_rej      = 3                                         # High rejection in sigma of fit

###### "Pre-Processing Smoothing"
   - pre_smooth             = True                                      # smooth after interpolation and before stacking
   - pre_smooth_shape       = 'gaussian'                                # gaussian,boxcar,mexican
   - pre_smooth_size        = 1                                         # kernel size

###### "Pre-Processing MASKING"
   - pre_mask               = True                                      # mask spectra after smoothing (stacks)
   - pre_msk_abs_lines      = True                                      # mask IS absorptions lines
   - pre_mask_type          = 'NaN'                                     # continuum/constant/NaN
   - pre_mask_cte_val       = 0                                         # constant value for masking
   - pre_mask_lw            = 2                                         # line width (A)
   - pre_mask_blue_regn     = True                                      # mask initial spectra pixels
   - pre_mask_blue_regn_int = 300                                       # intial pix
   - pre_mask_blue_regn_fnl = 912                                       # final pix

###### "Sigma-Clip"
   - sigma_clipping         = True                                      # Sigma clipping
   - sigma_cut              = 3                                         # sigma cut
   - sigma_cen_fct          = mean                                      # median, mean
   - sigma_msk_fill_val     = np.nan                                    # np.nan, value

###### #Weighting"
   - weight_type            = 'cont-flux-med'                           # i-band-mag,cont-flux-sum,cont-flux-med,cont-flux-avg None:
   - weight_cnt_flux_get    = True                                      # mask any given wavelength region 
   - weight_cnt_flux_lmb_0  = 1430                                      # initial lambda
   - weight_cnt_flux_lmb_n  = 1480                                      # final lambda


###### "Noise Files"
   - spectra_noise          = False                                     #Include Noise files in the Stacks

###### "Stacks Post Processing "
   - post_continuum         = False                                     # Fit Cont after stacking
   - post_cont_typ          = 'ratio'                                   # Continuum fitting type fit,ratio,difference
   - post_cont_lines        = '*'                                       # Image lines to be fit
   - post_cont_funct        = 'spline3'                                 # Fitting function: legendre, chebyshev, spline1, spline3
   - post_cont_order        = 9                                         # Order Polynomial / num pieces spline
   - post_cont_override     = 'yes'                                     # Override previous norm spec
   - post_cont_replace      = 'no'                                      # Replace rejected points by fit?
   - post_cont_low_rej      = 5                                         # Low rejection in sigma of fit
   - post_cont_high_rej     = 5                                         # High rejection in sigma of fit
   - post_smooth            = True                                      # smooth after stacking
   - post_smooth_shape      = 'gaussian'                                # smooth after stacking
   - post_smooth_size       = 1                                         # smooth after stacking


![Alt text](./Images/step.jpg?raw=true "Pre-processing of stacked spetra.")
## Lines Dictionary
   - Contains a list of identified emmission and absorption lines in the UV spectral range [800-3000]AA from [IAC's OTELO spectral line summary](http://research.iac.es/proyecto/otelo/pages/data-tools/spectral-line-summary.php), [Shapley+03](https://ui.adsabs.harvard.edu/abs/2003ApJ...588...65S/abstract), [Halliday+08](https://ui.adsabs.harvard.edu/abs/2008A%26A...479..417H/abstract) and [Le Fevre+15](https://ui.adsabs.harvard.edu/abs/2015A%26A...576A..79L/abstract). However this can de modified. The list contains different parameters for the Stacking procedure:
     - 0. Central wavelength.
     - 1. Line width.
     - 2. Line width region for plotting purposes.
     - 3. Line id for stacking routines. 
     - 4. Line id for labelling plots.
     - 5. Line id for fits headers.
     - 6. Line marker for plotting purposes.
     - 7. Line fitting: Central wavelength constrains.
     - 8. Line fitting: Central wavelength offset.
     - 9. Line id for stacking routines.
     - 10. Line fitting: Amplitude bound constrains.

## Plots

## Line fitting
VSAT uses lmfit for line fitting and by default a simple gaussian is used as a line profile. It is possible to use GM mode to fit multiple gaussian emmission before and after the central wavelength of some lines as CIV. The line fitting procedure is performed under a defined wavelength range which can contain a single line or multiple lines as defined in the Line_dictionary.py file. 

![Alt text](./Images/FitSingle.jpg?raw=true "Pre-processing of stacked spetra.")


![Alt text](./Images/FitMultiple.jpg?raw=true "Pre-processing of stacked spetra.")
## Examples
###### "Stacking"
The following snippet will stack galaxies from the COSMOS field. 

First we define a subsample of 102 galaxies according to redshift wuality flag and their impact parameter. 
```python
stamps_subsample_redshift_flag_fg     = Select_Subsamples(op_tbl_F                                ,'redshift_bk_flag',z_flag_itv_bg, test_fg = False, test_bg = False, slc_int = False)
stamps_subsample_redshift_flag_fg_bg  = Select_Subsamples(stamps_subsample_redshift_flag_fg[0][-1],'redshift_fg_flag',z_flag_itv_fg, test_fg = False, test_bg = False, slc_int = False)
fg_in                                 = str(stamps_subsample_redshift_flag_fg_bg[0][-1])
stamps_subsample_sep_fg               = Select_Subsamples(fg_in,'sep_as'     ,SEP_as_itv_23,z_flag_itv_fg, test_fg = False, test_bg = False)#, slc_smp=False, sel_pre_cnt = selec_spec_contn)
```
Next we stack the subsample. 
```python
f                                     = np.array(Stack_Subsample(stamps_subsample_sep_fg      ,
					sel_pre_shf     = selec_spec_shift   ,sel_pre_cnt     = selec_spec_contn      ,sel_pre_msk     = selec_spec_masks      ,
					pre_cnt         = pre_continuum      ,pre_cnt_typ     = pre_cont_typ          ,pre_cnt_lns     = pre_cont_lines        ,
					pre_cnt_fnc     = pre_cont_funct     ,pre_cnt_ord     = pre_cont_order        ,pre_cnt_ovr     = pre_cont_override     ,
					pre_cnt_rpl     = pre_cont_replace   ,pre_cnt_lrj     = pre_cont_low_rej      ,pre_cnt_hrj     = pre_cont_high_rej     ,
					smt_spc_pre     = pre_smooth         ,smt_shp_pre     = pre_smooth_shape      ,smt_sze_pre     = pre_smooth_size       ,
					pre_msk         = False              ,
					pre_msk_typ     = pre_mask_type      ,pre_msk_abs_lne = False                 ,pre_msk_cte_val = pre_mask_cte_val      ,
					pre_msk_blu_rgn = False              ,pre_blu_lmb_min = pre_mask_blue_regn_int,pre_blu_lmb_max = pre_mask_blue_regn_fnl,
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
###### "Stats"
Statistical values from the stacked galaxies can be obtained through:
```python
[stats_table(tblnm,tbl_format_opt) for tblnm in stamps_subsample_sep_fg[0]]
```

###### "Plots"
The following snippet will plot the stacked spectra continuum fit and smoothed, including
an histogram of the number of spectra combined per wavelength.
```python
Plot_All_Spec_All_Int(
			frgrnd_plt     = True,
			bkgrnd_plt     = False,
			n_int_spt      = n_int_splt_by    , int_typ_spl     = lst_spt_prp,
			min_x_lim_Idp  = lambda_low       , max_x_lim_Idp   = lambda_hgh,
			plt_ind_spec   = plot_spectra_idv , plt_cnt_stk_spc = plot_spectra_stk,
			wgt_typ        = weight_type,
			autoaxis_Idp   = False            , aaxs_Idp_ml_y   = True, 
			min_y_lim_Idp  = 0.5              , max_y_lim_Idp   = 2.0,
			lower_shift    = 0                , upper_shift     = 0, 
			only_stt_tbl   = False            , 
			SNR_lines      = 'None'           , show_legends    = True,
			max_sep        = max_red_sep      ,
			mlt_stk_fct    = plt_fit_stk_fct  ,
			fpt_foreground = fit_plt_fg       ,fpt_background   = fit_plt_bg,
			plt_stk_med    = plot_stack_med   ,plt_stk_avg      = plot_stack_avg  ,plt_stk_avw = plot_stack_avw)

```
![Alt text](./Images/Stacked.jpg?raw=true "Stacked spectra computed COSMOS field.")

It is also possible to plot the individual spectra, if plt_cnt_stk_spc is set True, generating:

![Alt text](./Images/Stacked-Contribution.jpg?raw=true "Stacked spectra COSMOS field.")

It is also possible to plot the individual spectra files
![Alt text](./Images/Spec-Individual.jpg?raw=true "Stacked spectra COSMOS field.")

And the pre-processing steps oof each spectra before combining them.
![Alt text](./Images/Spec-Step.jpg?raw=true "Stacked spectra COSMOS field.")
###### "Line Fitting"
```python
for function in plt_fit_stk_fct:
	lambda_low        = 1300
	lambda_hgh        = 1305
	Plot_Idp_Spc_Lne(
					int_typ_spl     = lst_spt_prp     ,stk_function   = function + fct_extra,
					lmb_min         = lambda_low      ,lmb_max        = lambda_hgh ,
					fit_type        = fitting_m       ,fit_fnct       = fitting_f ,
					verbose         = True            ,autoaxis       = True ,
					pre_off_plt     = False           ,ofs_ctr_fit    = False ,
					n_int_spt       = n_int_splt_by   ,
					lower_shift     = 0               ,upper_shift    = 1            ,   #0-1 All
					max_sep         = max_red_sep     ,
					mlt_stk_fct     = plt_fit_stk_fct ,
					mke_lne_fit     = mke_new_lne_fit , 
					fit_vls_hdr     = fit_upd_hdr     ,
					int_vlf_hdr     = fit_ivf_hdr     ,
					uft_lne_vls     = fit_luf_hdr     ,
					cnt_bnp_adj     = cnt_bnp_reg     ,
					fpt_foreground  = fit_plt_fg      ,fpt_background = fit_plt_bg,
					fix_ctr_gau     = fix_ctr_gaussian,
					fix_pre_gau     = fix_pre_gaussian,
					fix_pst_gau     = fix_pst_gaussian,
					fix_ctr_gau_1   = fix_gau_1       ,
					fix_ctr_gau_2   = fix_gau_2		  ,
					pre_shf_lim     = pre_gauss_shf   ,pst_shf_lim    = pst_gauss_shf,
					pre_shf_ctr     = pre_gauss_ctr   ,pst_shf_ctr    = pst_gauss_ctr,
					fix_mdl_gau     = fix_mdl_gaussian,
					mdl_shf_ctr     = mdl_gauss_ctr   ,mdl_shf_lim    = mdl_gauss_shf,
					ivl_fts_hdr     = int_vls_prv_fit
					)
	
	

![Alt text](./Images/LINE-FIT-COSMOS-avg-c-smt-G-Ind-Splt.jpg?raw=true "Stacked spectra COSMOS field.")	
```python
	Plot_Slc_Spc_Lne(
					int_typ_spl    = lst_spt_prp    ,stk_function   = function + fct_extra ,
					lmb_min        = lambda_low     ,lmb_max        = lambda_hgh ,
					stk_fct        = plt_fit_stk_fct,
					fit_type       = fitting_m      ,fit_fnct       = fitting_f     ,
					verbose        = True           ,autoaxis       = True          ,
					pre_off_plt    = False          ,n_int_spt      = n_int_splt_by ,
					lower_shift    = 0              ,upper_shift    = 1             ,
					plt_ind_fit    = slcs_plt_fit   ,
					autoaxis_SSL   = True           ,
					lbl_col_idv    = True           ,#nmb_cols      = 2,
					fpt_foreground = fit_plt_fg     ,fpt_background = fit_plt_bg,
					max_sep        = max_red_sep    ,
					empty_plots    = 2              ,landscape_plt  = True,
					splt_ind_lns   = slcs_plt_ind)
```
![Alt text](./Images/LINE-FIT-COSMOS-avg-c-smt-G-Mlt-Splt?raw=true "Stacked spectra COSMOS field.")	
## Dependencies
Currently VSAT works only with astropy 2.0 as it relies on pyraf continuum task for continuum normalization. However a new version will be released dropping this dependency.
 - [astropy](https://www.astropy.org)
 - [bottleneck](https://pypi.org/project/Bottleneck/)
 - [pandas](https://pandas.pydata.org)
 - [scipy](https://www.scipy.org)
 - [numpy](https://numpy.org)
 - [lmfit](https://lmfit.github.io/lmfit-py/)
 - [pysynphot](https://pysynphot.readthedocs.io/en/latest/)
 - [pyraf](https://astroconda.readthedocs.io/en/latest/installation.html)
 - [termcolor](https://pypi.org/project/termcolor/)
## License

###### "Scientist"



