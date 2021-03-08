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

###### "Plots"
```python

Plot_All_Spec_All_Int(n_int_spt      = n_int_splt_by    , int_typ_spl     = lst_spt_prp,
			min_x_lim_Idp  = lambda_low       , max_x_lim_Idp   = lambda_hgh,
			plt_ind_spec   = plot_spectra_idv , plt_cnt_stk_spc = plot_spectra_stk,
			wgt_typ        = weight_type,
			autoaxis_Idp   = False            , aaxs_Idp_ml_y   = True, 
			min_y_lim_Idp  = 0.5              , max_y_lim_Idp   = 2.0,
			lower_shift    = 0                , upper_shift     = 0, 
			only_stt_tbl   = False            , 
			SNR_lines      = 'None'           , show_legends    = True,
			max_sep        = max_red_sep      ,
			mlt_stk_fct    = plt_fit_stk_fct             ,
			fpt_foreground = fit_plt_fg       ,fpt_background   = fit_plt_bg,
			plt_stk_med    = plot_stack_med   ,plt_stk_avg      = plot_stack_avg  ,plt_stk_avw = plot_stack_avw)

```
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



