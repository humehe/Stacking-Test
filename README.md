# ValStack
Stacking 1D spectra tools

Contains:

Fnc_Stk_Dir.py Defines: 

Location of the input catalogue and spectral data. 

Parameters for selecting galaxy pairs and subsamples of galaxies according to their physical properties. 

Location of the resulting products of the stacking analyses e.g. tables, plots, pre stacking processed spectra and stacked spectra.


Fnc_Stk_Mth.py Contains useful math functions (e.g. cosmological constants, gaussian, lorentzian and voigt profiles for line emmision/absorption) for the stacking anaysis.


Fnc_Stk_Tbl.py Contains functions to access different tables. 


Fnc_Stk_Fts.py Contains funtions to access and modify (add, modify, delete) fits headers.


Fnc_Stk_Plt.py Contains different plot templates used through all the stacking analysis. 


Fnc_Stk_Spc.py Contains a series of functions for pre-processing the spectra needed for the stacking analysis


Fnc_Stk_Stk.py Contains the core of the stacking tool, a bootstrap function to compute the CIs of the stacked spectra, a function to select different galaxy subsamples for stacking and a tool to fit the different emmision/absorption lines through a simple/multiple component gaussian profiles.


Fnc_Stk_Utl.py Contains auxiliary functions for the stacking analysis


Lines_Dictionary.py Contains a list of identified emmission and absorption lines in the UV spectral range [800-3000]AA.


![Alt text](../VUDS-Z.pdf?raw=true "Title")
