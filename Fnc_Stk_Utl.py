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
from progressbar import *               # just a simple progress bar
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

from Fnc_Stk_Dir import *
####Fnc_Stk_Utl####
def Delete_Element_Array(array2bused,string2bdeleted):
	element2bdeleted = np.where(array2bused==string2bdeleted)
	array2bused = np.delete(array2bused, element2bdeleted)
	return array2bused

def Prev_Next(some_iterable,*args, **kwargs):
    prevs, items, nexts = tee(some_iterable, 3)
    prevs = chain([None], prevs)
    nexts = chain(islice(nexts, 1, None), [None])
    return izip(prevs, items, nexts)

def Def_Sub_Dirs_Slice_xtr(Prt_Dir,Slices,*args, **kwargs):
	sub_dir_slc  = []
	slices_split = []
	sub_dir_slc.append(Prt_Dir + str(Slices[0]) + '-' + str(Slices[-1]) + '/')
	slices_split.append(str(Slices[0]) + '-' + str(Slices[-1]))
	return sub_dir_slc,slices_split

def Def_Sub_Dirs_Slice_all(Prt_Dir,Slices,*args, **kwargs):
	sub_dir_slc  = []
	slices_split = []
	for previous, item, nxt in Prev_Next(Slices):
		slc_int_bin = Slices[1]-Slices[0]
		if item < Slices[-1]:
		   sub_dir_slc.append(Prt_Dir + str(item) + '-' + str(nxt) + '/')
		   slices_split.append(str(item) + '-' + str(nxt))
		else:
		    break
	sub_dir_slc.append(Prt_Dir + str(Slices[0]) + '-' + str(Slices[-1]) + '/')
	slices_split.append(str(Slices[0]) + '-' + str(Slices[-1]))
	return sub_dir_slc,slices_split
####Fnc_Stk_Utl####