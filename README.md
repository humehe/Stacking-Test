# ValStack
Stacking 1D spectra tools

1. First list item
   - First nested list item
     - Second nested list item
     - Third
      - Fourth
Contains:

1. Fnc_Stk_Dir.py:
   - Location of the input catalogue and spectral data. 
   - Parameters for selecting galaxy pairs and subsamples of galaxies according to their physical properties. 
   - Location of the resulting products of the stacking analyses e.g. tables, plots, pre stacking processed spectra and stacked spectra.


2. Fnc_Stk_Mth.py:
   - Contains useful math functions (e.g. cosmological constants, gaussian, lorentzian and voigt profiles for line emmision/absorption) for the stacking anaysis.


3. Fnc_Stk_Tbl.py 
   - Functions to access different tables. 


4. Fnc_Stk_Fts.py 
   -Funtions to access and modify (add, modify, delete) fits headers.

5. Fnc_Stk_Plt.py
   - Contains different plot templates used through all the stacking analysis. 

6. Fnc_Stk_Spc.py 
   - Contains a series of functions for pre-processing the spectra needed for the stacking analysis

7. Fnc_Stk_Stk.py 
   - Contains the core of the stacking tool
   - Bootstrap function to compute the CIs of the stacked spectra. 
   - Function to select different galaxy subsamples for stacking.
   - Fitting tool for the different emmision/absorption lines through a simple/multiple component gaussian profiles.


8. Fnc_Stk_Utl.py 
   - Contains auxiliary functions for the stacking analysis


Lines_Dictionary.py 
   - Contains a list of identified emmission and absorption lines in the UV spectral range [800-3000]AA.


![Alt text](./Images/VUDS-Z.jpg?raw=true "Title")

