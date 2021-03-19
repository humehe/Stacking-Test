import math

from astropy import coordinates as apcoo
from astropy import cosmology as apcosmo
#from astropy import cosmology
#from astropy.cosmology import apcosmo.FlatLambdaCDM
from astropy.io import ascii
#from astropy import convolution

from scipy.constants import physical_constants
from scipy.special import wofz
from astropy import units as u

from astropy import convolution as apcvl
from astropy import table as aptbl

from Fnc_Stk_Dir import *
from Fnc_Stk_Fts import *
from Fnc_Stk_Tbl import *

####Fnc_Stk_Mth####
def proj_sep(ra_o1,dec_o1,ra_o2,dec_o2):
	#http://cads.iiap.res.in/tools/angularSeparation
	#http://celestialwonders.com/tools/starAngleCalc.html
	#http://docs.astropy.org/en/stable/coordinates/matchsep.html
	#https://en.wikipedia.org/wiki/Great-circle_distance
	#http://www.gyes.eu/calculator/calculator_page1.htm

	ra_dif       = abs(ra_o1   - ra_o2)
	dec_dif      = abs(dec_o1  - dec_o2)
	dec_mean     = abs((dec_o1 + dec_o2)/2)

	theta        = math.degrees(math.atan(dec_dif/ra_dif))
	theta2       = math.degrees(math.atan2(dec_dif,ra_dif))

	cosdec_mean  = math.cos(math.radians(dec_mean))
	sep_deg      = math.sqrt(ra_dif**2+(dec_dif*cosdec_mean)**2)
	sep_am       = sep_deg*60
	sep_as       = sep_deg*3600

	c1=apcoo.SkyCoord(ra_o1,dec_o1,frame='icrs',unit='deg')
	c2=apcoo.SkyCoord(ra_o2,dec_o2,frame='icrs',unit='deg')
	c3=c1.separation(c2)
 	center = c1
	target = c2
	aframe = center.skyoffset_frame()
	new_coordintes = target.transform_to(aframe)

	new_lon_o2 = (new_coordintes.lon*u.deg).value
	new_lat_o2 = (new_coordintes.lat*u.deg).value

	return sep_deg,sep_am,sep_as,new_lon_o2,new_lat_o2,ra_dif,dec_dif,dec_mean,theta,theta2

def Select_Pairs(Slct_Pair_ip_tbl,*args, **kwargs):
	sfx_tbl_otp = kwargs.get('sfx_tbl_otp', '')
	verbose     = kwargs.get('verbose', False)
	cp_spc_fls  = kwargs.get('cp_spc_fls',False)
	cp_spn_fls  = kwargs.get('cp_spn_fls',False)
	ind_tbl     = kwargs.get('ind_tbl',True)
	grl_tbl     = kwargs.get('grl_tbl',True)

	cat      = aptbl.Table.read(Slct_Pair_ip_tbl, format=tbl_format_ipt)

	ident    = cat['ident']
	ar       = cat['alpha']
	dec      = cat['delta']
	spc_f    = cat['spec1d']
	spc_f_n  = cat['spec1dnoise']
	z        = cat['z_spec']
	zf       = cat['zflags']
	point    = cat['pointing']
	quad     = cat['quad']
	slit     = cat['slit']
	obj      = cat['obj']
	magi     = cat['magi']

	z_cat    = np.mean(z)


	if 'PRP' in Slct_Pair_ip_tbl and not 'MRP' in Slct_Pair_ip_tbl:
		Model        = cat['Model']
		ebv          = cat['ebv']
		chi_1        = cat['chi']
		Age          = cat['Age']
		age_low68    = cat['age_low68']
		age_high68   = cat['age_high68']
		mass         = cat['mass']
		mass_low68   = cat['mass_low68']
		mass_high68  = cat['mass_high68']
		sfr          = cat['sfr']
		sfr_low68    = cat['sfr_low68']
		sfr_high68   = cat['sfr_high68']
		ssfr         = cat['ssfr']
		ssfr_low68   = cat['ssfr_low68']
		ssfr_high68  = cat['ssfr_high68']
		L_nuv_rest   = cat['L_nuv_rest']
		M_FUV_Galex  = cat['M_FUV_Galex']
		M_NUV_Galex  = cat['M_NUV_Galex']
		M_u_CFHT     = cat['M_u_CFHT']
		M_B_Subaru   = cat['M_B_Subaru']
		M_V_Subaru   = cat['M_V_Subaru']
		M_g_Subaru   = cat['M_g_Subaru']
		M_r_Subaru   = cat['M_r_Subaru']
		M_i_Subaru   = cat['M_i_Subaru']
		M_z_Subaru   = cat['M_z_Subaru']
		M_J_VISTA    = cat['M_J_VISTA']
		M_K_VISTA    = cat['M_K_VISTA']
		M_i_CFHT     = cat['M_i_CFHT']
		M_IRAC1      = cat['M_IRAC1']
		M_IRAC2      = cat['M_IRAC2']
		M_IRAC3      = cat['M_IRAC3']
		M_IRAC4      = cat['M_IRAC4']
	elif 'MRP' in Slct_Pair_ip_tbl and not 'PRP' in Slct_Pair_ip_tbl:
		re           = cat['re']
		re_err       = cat['re_err']
		n            = cat['n']
		n_err        = cat['n_err']
		q            = cat['q']
		q_err        = cat['q_err']
		theta        = cat['theta']
		theta_err    = cat['theta_err']
		chi_2        = cat['chi']
		flag         = cat['flag']
	elif 'PRP' in Slct_Pair_ip_tbl and 'MRP' in Slct_Pair_ip_tbl:
		Model       = cat['Model']
		ebv         = cat['ebv']
		chi_1       = cat['chi_1']
		Age         = cat['Age']
		age_low68   = cat['age_low68']
		age_high68  = cat['age_high68']
		mass        = cat['mass']
		mass_low68  = cat['mass_low68']
		mass_high68 = cat['mass_high68']
		sfr         = cat['sfr']
		sfr_low68   = cat['sfr_low68']
		sfr_high68  = cat['sfr_high68']
		ssfr        = cat['ssfr']
		ssfr_low68  = cat['ssfr_low68']
		ssfr_high68 = cat['ssfr_high68']
		L_nuv_rest  = cat['L_nuv_rest']
		M_FUV_Galex = cat['M_FUV_Galex']
		M_NUV_Galex = cat['M_NUV_Galex']
		M_u_CFHT    = cat['M_u_CFHT']
		M_B_Subaru  = cat['M_B_Subaru']
		M_V_Subaru  = cat['M_V_Subaru']
		M_g_Subaru  = cat['M_g_Subaru']
		M_r_Subaru  = cat['M_r_Subaru']
		M_i_Subaru  = cat['M_i_Subaru']
		M_z_Subaru  = cat['M_z_Subaru']
		M_J_VISTA   = cat['M_J_VISTA']
		M_K_VISTA   = cat['M_K_VISTA']
		M_i_CFHT    = cat['M_i_CFHT']
		M_IRAC1     = cat['M_IRAC1']
		M_IRAC2     = cat['M_IRAC2']
		M_IRAC3     = cat['M_IRAC3']
		M_IRAC4     = cat['M_IRAC4']
		re          = cat['re']
		re_err      = cat['re_err']
		n           = cat['n']
		n_err       = cat['n_err']
		q           = cat['q']
		q_err       = cat['q_err']
		theta       = cat['theta']
		theta_err   = cat['theta_err']
		chi_2       = cat['chi_2']
		flag        = cat['flag']
	else:
		print
		print colored('No prefix added (PRP, MRP, PRP-MRP)! Reading regular table.','yellow')
		print

	cosmo = apcosmo.FlatLambdaCDM(70,0.3)

	Cosmo_Scale_Par(0.0279,cosmo)
	Cosmo_Scale_Par(2.280,cosmo)

	Cosmo_Scale_Par(z_cat,cosmo,cosmo_check=True)
	
	cat.sort('z_spec')

	DELTAZ_P    = []
	SEP_DG_P    = []
	SEP_AM_P    = []
	SEP_AS_P    = []
	SEP_KPC_P   = []
	SCALEKPC_P  = []
	NEW_RA_P    = []
	NEW_DEC_P   = []

	DIF_ALP_P   = []
	DIF_DLT_P   = []
	OMG1_P      = []
	OMG2_P      = []
	PHI_P       = []
	PHIA_P      = []


	indexes_B    = []
	PAIR_B       = []
	PAIR_Z_B     = []
	PAIR_ZF_B    = []
	PAIR_SFN_B   = []
	PAIR_SFN_N_B = []
	PAIR_MAGI_B  = []

	indexes_F    = []
	PAIR_F       = []
	PAIR_Z_F     = []
	PAIR_ZF_F    = []
	PAIR_SFN_F   = []
	PAIR_SFN_N_F = []
	PAIR_MAGI_F  = []

	PAIR_MASS_B  = []
	PAIR_Age_B   = []
	PAIR_sfr_B   = []
	PAIR_ssfr_B  = []
	PAIR_L_nuv_B = []

	PAIR_MASS_F  = []
	PAIR_Age_F   = []
	PAIR_sfr_F   = []
	PAIR_ssfr_F  = []
	PAIR_L_nuv_F = []

	PAIR_re_B        = []
	PAIR_re_err_B    = []
	PAIR_n_B         = []
	PAIR_n_err_B     = []
	PAIR_q_B         = []
	PAIR_q_err_B     = []
	PAIR_theta_B     = []
	PAIR_theta_err_B = []
	PAIR_chi_2_B     = []
	PAIR_flag_B      = []

	PAIR_re_F        = []
	PAIR_re_err_F    = []
	PAIR_n_F         = []
	PAIR_n_err_F     = []
	PAIR_q_F         = []
	PAIR_q_err_F     = []
	PAIR_theta_F     = []
	PAIR_theta_err_F = []
	PAIR_chi_2_F     = []
	PAIR_flag_F      = []

	if slt_prs == True:
		widgets = ['Evaluating possible pairs for '+ str(len(z)) + ' galaxies: ', 
		Percentage(), ' ', Bar(marker='>',left='[',right=']'),
		' ', ETA(), ' ', FileTransferSpeed()]
		pbar    = ProgressBar(widgets=widgets, maxval=len(z))
		pbar.start()
		print

		for j in range(len(z)):
			pbar.update(j)
			zref          = z[j]
			indexes_FB    = []

			DELTAZ_FB     = []
			SEP_DG_FB     = []
			SEP_AM_FB     = []
			SEP_AS_FB     = []
			SEP_KPC_FB    = []
			SCALEKPC_FB   = []
			NEW_RA_FB     = []
			NEW_DEC_FB    = []

			DIF_ALP_FB    = []
			DIF_DLT_FB    = []
			OMG1_FB       = []
			OMG2_FB       = []
			PHI_FB        = []
			PHIA_FB       = []

			PAIR_FB       = []
			PAIR_Z_FB     = []
			PAIR_ZF_FB    = []
			PAIR_MAGI_FB  = []

			PAIR_SFN_FB_F   = []
			PAIR_SFN_N_FB_F = []

			PAIR_SFN_FB_B   = []
			PAIR_SFN_N_FB_B = []

			if zref >= z_lim_inf:

				scale_pair  = Cosmo_Scale_Par(zref,cosmo)[0]
				scaleL_pair = Cosmo_Scale_Par(zref,cosmo)[1]
				
				print ''
				print colored('Foreground galaxy to be evaluated: ','yellow')
				print colored('index id ar dec z z_ref','yellow')
				print colored(str(j)+ ' ' + str(ident[j]) + ' ' + str(ar[j]) + ' ' + str(dec[j]) + ' ' + str(z[j]) + ' ' + str(zref),'yellow')
				print

				indexes_FB.append(j)
				PAIR_SFN_FB_F.append(str(spc_f[j]))
				PAIR_SFN_N_FB_F.append(str(spc_f_n[j]))
				DELTAZ_FB.append(0)
				SEP_DG_FB.append(0)
				SEP_AM_FB.append(0)
				SEP_AS_FB.append(0)
				SEP_KPC_FB.append(0)
				SCALEKPC_FB.append(0)
				NEW_RA_FB.append(0)
				NEW_DEC_FB.append(0)

				DIF_ALP_FB.append(0)
				DIF_DLT_FB.append(0)
				OMG1_FB.append(0)
				OMG2_FB.append(0)
				PHI_FB.append(0)
				PHIA_FB.append(0)


				PAIR_FB.append(0)
				PAIR_Z_FB.append(0)
				PAIR_ZF_FB.append(0)
				PAIR_MAGI_FB.append(0)

				PAIR_SFN_FB_B.append(0)
				PAIR_SFN_N_FB_B.append(0)

				glxdir     = par_dir_res + '/' + str(ident[j]) +'/'

				for w in range(len(z)-j):
					zdif      = z[w+j]-zref
					sep_pairs = proj_sep(ar[j],dec[j],ar[w+j],dec[w+j])
					dif_dg    = sep_pairs[0]
					dif_am    = sep_pairs[1]
					dif_as    = sep_pairs[2]
					new_ar    = sep_pairs[3]
					new_dc    = sep_pairs[4]
					dif_kpc   = dif_as*scale_pair
					dif_alp   = sep_pairs[5]
					dlt_avg   = sep_pairs[6]
					dif_dlt   = sep_pairs[7]
					omg1      = sep_pairs[8]
					omg2      = sep_pairs[9]

					if '_PRP_MRP' in Slct_Pair_ip_tbl:

						pos_agl   = theta[w+j]

						if pos_agl > 90:
							pos_agl = 180 - pos_agl
						elif pos_agl < 90:
							pos_agl = 180+pos_agl
					else:
						pass

					if (zl1 <= zdif <= zl2) and (dif_as <= rad_sep[0][-1]):

						print
						print colored('Foreground-Background galaxy pair candiate is within the limits!','cyan')
						print colored('Foreground galaxy: '+str(ident[j])+' '+str(w),'yellow')
						print colored('Background galaxy: '+str(w+j)+' '+str(ident[w+j]),'yellow')
						print colored('Redshift separatiion: '+ str(zl1)+'<='+str(zdif)+'<='+str(zl2),'yellow')
						print colored('Separation: '+str(dif_as),'yellow')
						print
						

						if   '_PRP_MRP' in Slct_Pair_ip_tbl:
							####################COMPUTE AZIMUTHAL ANGLE (PHI) BETWEEN#################### 
							################THE PROJECTED sma(ROTATION AXIS ) of FG DISK#################
							#############AND THE PROJECTED VECTOR FROM THE CENTER OF FG TO BG############
							####################AS DEFINDED BY BORDOLOI+11 APJ_743_10####################
							###############PHI>45 REPRESENTS LOS THAT PASS NEAR THE PLANE################
							############PHI<45 REPRESENTS LOS THAT PASS NEAR THE SYMMETRY AXIS###########
							#phi       = abs(omg1)+abs(pos_agl)
							if   ((dif_alp>0) and (dif_dlt>0)  and (pos_agl>0) and (omg1<(90-pos_agl))) or ((dif_alp<0) and (dif_dlt<0) and (pos_agl>0) and (omg1<(90-pos_agl))) or ((dif_alp<0) and (dif_dlt>0) and (pos_agl<0) and (omg1<(90-pos_agl))) or ((dif_alp>0) and (dif_dlt<0) and (pos_agl<0) and (omg1<(90-pos_agl))):
								phi = pos_agl + abs(omg1)
							elif ((dif_alp>0) and (dif_dlt>0)  and (pos_agl>0) and (omg1>(90-pos_agl))) or ((dif_alp<0) and (dif_dlt<0) and (pos_agl>0) and (omg1>(90-pos_agl))) or ((dif_alp<0) and (dif_dlt>0) and (pos_agl<0) and (omg1>(90-pos_agl))) or ((dif_alp>0) and (dif_dlt<0) and (pos_agl<0) and (omg1>(90-pos_agl))):
								phi = 180 - (pos_agl + abs(omg1))
							elif ((dif_alp<0) and (dif_dlt>0)  and (pos_agl>0) and (omg1>(90-pos_agl))) or ((dif_alp>0) and (dif_dlt<0) and (pos_agl>0) and (omg1>(90-pos_agl))) or ((dif_alp>0) and (dif_dlt>0) and (pos_agl<0) and (omg1>(pos_agl))) or ((dif_alp<0) and (dif_dlt<0) and (pos_agl<0) and (omg1>(pos_agl))):
								phi = abs(omg1) - pos_agl
							elif ((dif_alp<0) and (dif_dlt>0)  and (pos_agl>0) and (pos_agl>omg1)) or ((dif_alp>0) and (dif_dlt<0) and (pos_agl>0) and (pos_agl>omg1)) or ((dif_alp>0) and (dif_dlt>0) and (pos_agl<0) and (pos_agl>omg1)) or ((dif_alp<0) and (dif_dlt<0) and (pos_agl<0) and (pos_agl>omg1)):
								phi = pos_agl - abs(omg1)
							else:
								print 'Else'
								pass

							if '_PRP_MRP' in Slct_Pair_ip_tbl:
								phia = abs(phi)
							else:
								pass

							####################COMPUTE AZIMUTHAL ANGLE (PHI) BETWEEN#################### 
							################THE PROJECTED sma(ROTATION AXIS ) of FG DISK#################
							#############AND THE PROJECTED VECTOR FROM THE CENTER OF FG TO BG############
							####################AS DEFINDED BY BORDOLOI+11 APJ_743_10####################
							###############PHI>45 REPRESENTS LOS THAT PASS NEAR THE PLANE################
							############PHI<45 REPRESENTS LOS THAT PASS NEAR THE SYMMETRY AXIS###########
						else:
							pass

						if not os.path.exists(glxdir):
							os.mkdir(glxdir)
							if verbose == True:
								print 
								print 'Creating directory      : ',glxdir
							elif verbose == False:
								pass
						elif os.path.exists(glxdir):
							if verbose == True:
								print
								print 'Directory already exists: ',glxdir
							elif verbose == False:
								pass

						specfile_fg   = spc_dir  + str(spc_f[j])
						specfile_n_fg = spc_dirn + str(spc_f_n[j])

						sfn_cpt   = glxdir + str(spc_f[j])
						snfn_cpt  = glxdir + str(spc_f_n[j])

						if cp_spc_fls == True and not os.path.exists(sfn_cpt):
							Header_Get_Add(specfile_fg,'h_s_c',0,header_comment='History Step Last')
							Header_Get_Add(specfile_fg,'h_s_0',str(spc_f[j].rsplit('.',1)[0]),header_comment='History Step Init')
							os.system('cp ' + specfile_fg   + ' ' + glxdir)
							if verbose == True:
								print
								print 'Copying foreground galaxy spectra            : ',specfile_fg
								print 'To Foreground galaxy directory               : ',glxdir
							elif verbose == False:
								pass
						elif cp_spc_fls == True and os.path.exists(sfn_cpt):
							if verbose == True:
								print
								print 'Foreground galaxy spectra                    : ',specfile_fg
								print 'Already Exists in Foreground galaxy directory: ',glxdir
							elif verbose == False:
								pass
						elif cp_spc_fls == False:
							pass

						if cp_spn_fls == True and not os.path.exists(snfn_cpt):
							Header_Get_Add(specfile_n_fg,'h_s_c',0,header_comment='History Step Last')
							Header_Get_Add(specfile_n_fg,'h_s_0',str(spc_f_n[j].rsplit('.',1)[0]),header_comment='History Step Init')
							os.system('cp ' + specfile_n_fg  + ' ' + glxdir)
							if verbose == True:
								print
								print 'Copying foreground galaxy noise pectra       : ',specfile_n_fg
								print 'To Foreground galaxy directory               : ',glxdir
							elif verbose == False:
								pass
						elif cp_spn_fls == True and os.path.exists(snfn_cpt):
							if verbose == True:
								print
								print 'Foreground galaxy spectra  (noise)           : ',specfile_n_fg
								print 'Already Exists in Foreground galaxy directory: ',glxdir
							elif verbose == False:
								pass
						elif cp_spn_fls == False:
							pass

						specfile_bg   = spc_dir  + str(spc_f[w+j])
						specfile_n_bg = spc_dirn + str(spc_f_n[w+j])

						sfn_cpt       = glxdir + str(spc_f[w+j])
						snfn_cpt      = glxdir + str(spc_f_n[w+j])

						if cp_spc_fls == True and not os.path.exists(sfn_cpt):
							Header_Get_Add(specfile_bg,'h_s_c',0,header_comment='History Step Last')
							Header_Get_Add(specfile_bg,'h_s_0',str(spc_f[w+j].rsplit('.',1)[0]),header_comment='History Step Init')
							os.system('cp ' + specfile_bg   + ' ' + glxdir)
							if verbose == True:
								print 
								print 'Copying background galaxy spectra            : ',specfile_bg
								print 'To Foreground galaxy directory               : ',glxdir
							elif verbose == False:
								pass
						if cp_spc_fls == True and os.path.exists(sfn_cpt):
							if verbose == True:
								print
								print 'Foreground galaxy spectra                    : ',specfile_bg
								print 'Already Exists in Foreground galaxy directory: ',glxdir
							elif verbose == False:
								pass
						elif cp_spc_fls == False:
							pass
						if cp_spn_fls == True and not os.path.exists(snfn_cpt):
							Header_Get_Add(specfile_n_bg,'h_s_c',0,header_comment='History Step Last')
							Header_Get_Add(specfile_n_bg,'h_s_0',str(spc_f_n[w+j].rsplit('.',1)[0]),header_comment='History Step Init')							
							os.system('cp ' + specfile_n_bg  + ' ' + glxdir)
							if verbose == True:
								print
								print 'Copying background galaxy noise pectra       : ',specfile_n_bg
								print 'To Foreground galaxy directory               : ',glxdir
							elif verbose == False:
								pass
						if cp_spn_fls == True and  os.path.exists(snfn_cpt):
							if verbose == True:
								print
								print 'Foreground galaxy spectra  (noise)           : ',specfile_n_bg
								print 'Already Exists in Foreground galaxy directory: ',glxdir
							elif verbose == False:
								pass
						elif cp_spn_fls == False:
							pass

						##GENERAL TABLE##
						DELTAZ_P.append(str(zdif))
						SEP_DG_P.append(str(dif_dg))
						SEP_AM_P.append(str(dif_am))
						SEP_AS_P.append(str(dif_as))
						SEP_KPC_P.append(str(dif_kpc))
						SCALEKPC_P.append(str(scale_pair))
						NEW_RA_P.append(str(new_ar))
						NEW_DEC_P.append(str(new_dc))

						indexes_F.append(j)
						PAIR_F.append(str(ident[w+j]))
						PAIR_Z_F.append(str(z[w+j]))
						PAIR_ZF_F.append(str(zf[w+j]))
						PAIR_MAGI_F.append(str(magi[w+j]))

						PAIR_SFN_F.append(str(spc_f[w+j]))
						PAIR_SFN_N_F.append(str(spc_f_n[w+j]))

						indexes_B.append(w+j)
						PAIR_B.append(str(ident[j]))
						PAIR_Z_B.append(str(z[j]))
						PAIR_ZF_B.append(str(zf[j]))
						PAIR_MAGI_B.append(str(magi[j]))

						PAIR_SFN_B.append(str(spc_f[j]))
						PAIR_SFN_N_B.append(str(spc_f_n[j]))

						if 'PRP' in Slct_Pair_ip_tbl and not 'MRP' in Slct_Pair_ip_tbl:
							PAIR_MASS_F.append(str(mass[w+j]))
							PAIR_Age_F.append(str(Age[w+j]))
							PAIR_sfr_F.append(str(sfr[w+j]))
							PAIR_ssfr_F.append(str(ssfr[w+j]))
							PAIR_L_nuv_F.append(str(L_nuv_rest[w+j]))
				
							PAIR_MASS_B.append(str(mass[j]))
							PAIR_Age_B.append(str(Age[j]))
							PAIR_sfr_B.append(str(sfr[j]))
							PAIR_ssfr_B.append(str(ssfr[j]))
							PAIR_L_nuv_B.append(str(L_nuv_rest[j]))
						elif 'MRP' in Slct_Pair_ip_tbl and not 'PRP' in Slct_Pair_ip_tbl:
							PAIR_re_F.append(str(re[w+j]))
							PAIR_re_err_F.append(str(re_err[w+j]))
							PAIR_n_F.append(str(n[w+j]))
							PAIR_n_err_F.append(str(n_err[w+j]))
							PAIR_q_F.append(str(q[w+j]))
							PAIR_q_err_F.append(str(q_err[w+j]))
							PAIR_theta_F.append(str(theta[w+j]))
							PAIR_theta_err_F.append(str(theta_err[w+j]))
							PAIR_chi_2_F.append(str(chi_2[w+j]))
							PAIR_flag_F.append(str(flag[w+j]))

							PAIR_re_B.append(str(re[j]))
							PAIR_re_err_B.append(str(re_err[j]))
							PAIR_n_B.append(str(n[j]))
							PAIR_n_err_B.append(str(n_err[j]))
							PAIR_q_B.append(str(q[j]))
							PAIR_q_err_B.append(str(q_err[j]))
							PAIR_theta_B.append(str(theta[j]))
							PAIR_theta_err_B.append(str(theta_err[j]))
							PAIR_chi_2_B.append(str(chi_2[j]))
							PAIR_flag_B.append(str(flag[j]))

						elif 'PRP' in Slct_Pair_ip_tbl and 'MRP' in Slct_Pair_ip_tbl:
							print
							print colored('Indices for fg, bg:'+str(w)+','+str(j),'yellow')
							print colored('Index for bg: '+str(w+j),'yellow')
							print
							PAIR_MASS_F.append(str(mass[w+j]))
							PAIR_Age_F.append(str(Age[w+j]))
							PAIR_sfr_F.append(str(sfr[w+j]))
							PAIR_ssfr_F.append(str(ssfr[w+j]))
							PAIR_L_nuv_F.append(str(L_nuv_rest[w+j]))

							PAIR_MASS_B.append(str(mass[j]))
							PAIR_Age_B.append(str(Age[j]))
							PAIR_sfr_B.append(str(sfr[j]))
							PAIR_ssfr_B.append(str(ssfr[j]))
							PAIR_L_nuv_B.append(str(L_nuv_rest[j]))

							PAIR_re_F.append(str(re[w+j]))
							PAIR_re_err_F.append(str(re_err[w+j]))
							PAIR_n_F.append(str(n[w+j]))
							PAIR_n_err_F.append(str(n_err[w+j]))
							PAIR_q_F.append(str(q[w+j]))
							PAIR_q_err_F.append(str(q_err[w+j]))
							PAIR_theta_F.append(str(theta[w+j]))
							PAIR_theta_err_F.append(str(theta_err[w+j]))
							PAIR_chi_2_F.append(str(chi_2[w+j]))
							PAIR_flag_F.append(str(flag[w+j]))

							PAIR_re_B.append(str(re[j]))
							PAIR_re_err_B.append(str(re_err[j]))
							PAIR_n_B.append(str(n[j]))
							PAIR_n_err_B.append(str(n_err[j]))
							PAIR_q_B.append(str(q[j]))
							PAIR_q_err_B.append(str(q_err[j]))
							PAIR_theta_B.append(str(theta[j]))
							PAIR_theta_err_B.append(str(theta_err[j]))
							PAIR_chi_2_B.append(str(chi_2[j]))
							PAIR_flag_B.append(str(flag[j]))

							DIF_ALP_P.append(dif_alp)
							DIF_DLT_P.append(dif_dlt)
							OMG1_P.append(omg1)
							OMG2_P.append(omg2)
							PHI_P.append(phi)
							PHIA_P.append(phia)							

						else:
							pass

						##GENERAL TABLE##

						##INDIVIDUAL TABLE##
						PAIR_FB.append(str(ident[j]))
						PAIR_Z_FB.append(str(z[j]))
						PAIR_ZF_FB.append(str(zf[j]))
						PAIR_MAGI_FB.append(str(magi[j]))

						indexes_FB.append(w+j)
						DELTAZ_FB.append(str(zdif))
						SEP_DG_FB.append(str(dif_dg))
						SEP_AM_FB.append(str(dif_am))
						SEP_AS_FB.append(str(dif_as))
						SEP_KPC_FB.append(str(dif_kpc))
						SCALEKPC_FB.append(str(scale_pair))
						NEW_RA_FB.append(new_ar)
						NEW_DEC_FB.append(new_dc)

						if 'PRP' in Slct_Pair_ip_tbl and 'MRP' in Slct_Pair_ip_tbl:
							DIF_ALP_FB.append(dif_alp)
							DIF_DLT_FB.append(dif_dlt)
							OMG1_FB.append(omg1)
							OMG2_FB.append(omg2)
							PHI_FB.append(phi)
							PHIA_FB.append(phia)
						else:
							pass

						PAIR_SFN_FB_F.append(str(spc_f[j]))
						PAIR_SFN_N_FB_F.append(str(spc_f_n[j]))

						PAIR_SFN_FB_B.append(str(spc_f[w+j]))
						PAIR_SFN_N_FB_B.append(str(spc_f_n[w+j]))

						##INDIVIDUAL TABLE##
						print
						print colored('Creating individual pairs for foreground galaxy: ','yellow')
						print colored('id: '+str(ident[j]),'yellow')
						if ind_tbl == True:
							rt = cat[indexes_FB]

							rt                = aptbl.Table()
							rt['id_B']        = ident[indexes_FB]
							rt['z_B']         = z[indexes_FB]
							rt['zf_B']        = zf[indexes_FB]
							rt['magi_B']      = magi[indexes_FB]

							rt['id_F']        = PAIR_FB
							rt['z_F']         = PAIR_Z_FB
							rt['zf_F']        = PAIR_ZF_FB
							rt['magi_F']      = PAIR_MAGI_FB

							rt['DELTAZ']      = DELTAZ_FB
							rt['SEP_arcsec']  = SEP_AS_FB
							rt['arcsec/kpc']  = SCALEKPC_FB
							rt['SEP_kpc']     = SEP_KPC_FB
							rt['NEW_RA']      = NEW_RA_FB
							rt['NEW_DEC']     = NEW_DEC_FB

							if 'PRP' in Slct_Pair_ip_tbl and 'MRP' in Slct_Pair_ip_tbl:
								rt['DIF_ALP']     = DIF_ALP_FB
								rt['DIF_DLT']     = DIF_DLT_FB
								rt['OMG1']        = OMG1_FB
								rt['OMG2']        = OMG2_FB
								rt['PHI']         = PHI_FB
								rt['PHI_ABS']     = PHIA_FB
							else:
								pass

							rt['spc_f_F']     = PAIR_SFN_FB_F
							rt['spc_f_n_F']   = PAIR_SFN_N_FB_F

							rt['spc_f_B']     = PAIR_SFN_FB_B
							rt['spc_f_n_B']   = PAIR_SFN_N_FB_B


							rt.sort('SEP_arcsec')
							name = glxdir+str(ident[j])+'_pairs' + sfx_tbl_otp + tbl_ext_opt
							rt.write(name, format=tbl_format_opt, overwrite=True)
							if verbose == True:
								print ''
								print colored('Table with results by foreground galaxy','green')
								print colored(rt,'green')
								print
							elif verbose == False:
								pass
							print colored('Pairs in table: '+name,'green')
							print
						elif ind_tbl == False:
							pass
					elif zdif > zl2:
						pass
						break
					elif zl1>zdif:
						pass
					elif dif_as > rad_sep[0][-1]:
						pass

			elif zref < z_lim_inf:
				pass
		pbar.finish()
		if grl_tbl == True:
			rtF = aptbl.Table()
			rtB = aptbl.Table()

			rtF['id_F']       = ident[indexes_F]
			rtF['z_F']        = z[indexes_F]
			rtF['zf_F']       = zf[indexes_F]
			rtF['magi_F']     = magi[indexes_F]

			if 'PRP' in Slct_Pair_ip_tbl and not 'MRP' in Slct_Pair_ip_tbl:
				rtF['mass_F']     = mass[indexes_F]
				rtF['Age_F']      = Age[indexes_F]
				rtF['SFR_F']      = sfr[indexes_F]
				rtF['sSFR_F']     = ssfr[indexes_F]
				rtF['Lnuv_F']     = L_nuv_rest[indexes_F]
			elif 'MRP' in Slct_Pair_ip_tbl and not 'PRP' in Slct_Pair_ip_tbl:
				rtF['re_F']        = re[indexes_F]
				rtF['re_err_F']    = re_err[indexes_F]
				rtF['n_F']         = n[indexes_F]
				rtF['n_err _F']    = n_err [indexes_F]
				rtF['q_F']         = q[indexes_F]
				rtF['q_err_F']     = q_err[indexes_F]
				rtF['theta_F']     = theta[indexes_F]
				rtF['theta_err_F'] = theta_err[indexes_F]
				rtF['chi_2_F']     = chi_2[indexes_F]
				rtF['flag_F']      = flag[indexes_F]
			elif 'PRP' in Slct_Pair_ip_tbl and 'MRP' in Slct_Pair_ip_tbl:
				rtF['mass_F']      = mass[indexes_F]
				rtF['Age_F']       = Age[indexes_F]
				rtF['SFR_F']       = sfr[indexes_F]
				rtF['sSFR_F']      = ssfr[indexes_F]
				rtF['Lnuv_F']      = L_nuv_rest[indexes_F]
				rtF['re_F']        = re[indexes_F]
				rtF['re_err_F']    = re_err[indexes_F]
				rtF['n_F']         = n[indexes_F]
				rtF['n_err _F']    = n_err [indexes_F]
				rtF['q_F']         = q[indexes_F]
				rtF['q_err_F']     = q_err[indexes_F]
				rtF['theta_F']     = theta[indexes_F]
				rtF['theta_err_F'] = theta_err[indexes_F]
				rtF['chi_2_F']     = chi_2[indexes_F]
				rtF['flag_F']      = flag[indexes_F]
			else:
				pass
			rtF['id_B']       = PAIR_F
			rtF['z_B']        = PAIR_Z_F
			rtF['zf_B']       = PAIR_ZF_F
			rtF['magi_B']     = PAIR_MAGI_F

			if 'PRP' in Slct_Pair_ip_tbl and not 'MRP' in Slct_Pair_ip_tbl:
				rtF['mass_B']     = PAIR_MASS_F
				rtF['Age_B']      = PAIR_Age_F
				rtF['SFR_B']      = PAIR_sfr_F
				rtF['sSFR_B']     = PAIR_ssfr_F
				rtF['Lnuv_B']     = PAIR_L_nuv_F
			elif 'MRP' in Slct_Pair_ip_tbl and not 'PRP' in Slct_Pair_ip_tbl:
				rtF['re_B']        = PAIR_re_F
				rtF['re_err_B']    = PAIR_re_err_F
				rtF['n_B']         = PAIR_n_F
				rtF['n_err _B']    = PAIR_n_err_F
				rtF['q_B']         = PAIR_q_F
				rtF['q_err_B']     = PAIR_q_err_F
				rtF['theta_B']     = PAIR_theta_F
				rtF['theta_err_B'] = PAIR_theta_err_F
				rtF['chi_2_B']     = PAIR_chi_2_F
				rtF['flag_B']      = PAIR_flag_F

			elif 'PRP' in Slct_Pair_ip_tbl and  'MRP' in Slct_Pair_ip_tbl:
				rtF['mass_B']     = PAIR_MASS_F
				rtF['Age_B']      = PAIR_Age_F
				rtF['SFR_B']      = PAIR_sfr_F
				rtF['sSFR_B']     = PAIR_ssfr_F
				rtF['Lnuv_B']     = PAIR_L_nuv_F

				rtF['re_B']        = PAIR_re_F
				rtF['re_err_B']    = PAIR_re_err_F
				rtF['n_B']         = PAIR_n_F
				rtF['n_err _B']    = PAIR_n_err_F
				rtF['q_B']         = PAIR_q_F
				rtF['q_err_B']     = PAIR_q_err_F
				rtF['theta_B']     = PAIR_theta_F
				rtF['theta_err_B'] = PAIR_theta_err_F
				rtF['chi_2_B']     = PAIR_chi_2_F
				rtF['flag_B']      = PAIR_flag_F

				rtF['DIF_ALP']    = DIF_ALP_P
				rtF['DIF_DLT']    = DIF_DLT_P
				rtF['OMG1']       = OMG1_P
				rtF['OMG2']       = OMG2_P
				rtF['PHI']        = PHI_P
				rtF['PHI_ABS']    = PHIA_P

			else:
				pass
			rtF['DELTAZ']     = DELTAZ_P
			rtF['SEP_arcsec'] = SEP_AS_P
			rtF['arcsec/kpc'] = SCALEKPC_P
			rtF['SEP_kpc']    = SEP_KPC_P
			rtF['NEW_RA']     = NEW_RA_P
			rtF['NEW_DEC']    = NEW_DEC_P

			rtF['spc_f_F']    = spc_f[indexes_F]
			rtF['spc_f_n_F']  = spc_f_n[indexes_F]

			rtF['spc_f_B']    = PAIR_SFN_F
			rtF['spc_f_n_B']  = PAIR_SFN_N_F

			rtB['id_B']       = ident[indexes_B]
			rtB['z_B']        = z[indexes_B]
			rtB['zf_B']       = zf[indexes_B]
			rtB['magi_B']     = magi[indexes_B]

			if 'PRP' in Slct_Pair_ip_tbl and not 'MRP' in Slct_Pair_ip_tbl:
				rtB['mass_B']     = mass[indexes_B]
				rtB['Age_B']      = Age[indexes_B]
				rtB['SFR_B']      = sfr[indexes_B]
				rtB['sSFR_B']     = ssfr[indexes_B]
				rtB['Lnuv_B']     = L_nuv_rest[indexes_B]
			elif 'MRP' in Slct_Pair_ip_tbl and not 'PRP' in Slct_Pair_ip_tbl:
				rtB['re_B']        = re[indexes_B]
				rtB['re_err_B']    = re_err[indexes_B]
				rtB['n_B']         = n[indexes_B]
				rtB['n_err _B']    = n_err [indexes_B]
				rtB['q_B']         = q[indexes_B]
				rtB['q_err_B']     = q_err[indexes_B]
				rtB['theta_B']     = theta[indexes_B]
				rtB['theta_err_B'] = theta_err[indexes_B]
				rtB['chi_B']       = chi_2[indexes_B]
				rtB['flag_B']      = flag[indexes_B]

			elif 'PRP' in Slct_Pair_ip_tbl and  'MRP' in Slct_Pair_ip_tbl:
				rtB['mass_B']     = mass[indexes_B]
				rtB['Age_B']      = Age[indexes_B]
				rtB['SFR_B']      = sfr[indexes_B]
				rtB['sSFR_B']     = ssfr[indexes_B]
				rtB['Lnuv_B']     = L_nuv_rest[indexes_B]

				rtB['re_B']        = re[indexes_B]
				rtB['re_err_B']    = re_err[indexes_B]
				rtB['n_B']         = n[indexes_B]
				rtB['n_err _B']    = n_err [indexes_B]
				rtB['q_B']         = q[indexes_B]
				rtB['q_err_B']     = q_err[indexes_B]
				rtB['theta_B']     = theta[indexes_B]
				rtB['theta_err_B'] = theta_err[indexes_B]
				rtB['chi_B']       = chi_2[indexes_B]
				rtB['flag_B']      = flag[indexes_B]
			else:
				pass
			rtB['id_F']       = PAIR_B
			rtB['z_F']        = PAIR_Z_B
			rtB['zf_F']       = PAIR_ZF_B
			rtB['magi_F']     = PAIR_MAGI_B

			if 'PRP' in Slct_Pair_ip_tbl and not 'MRP' in Slct_Pair_ip_tbl:
				rtB['mass_F']     = PAIR_MASS_B
				rtB['Age_F']      = PAIR_Age_B
				rtB['SFR_F']      = PAIR_sfr_B
				rtB['sSFR_F']     = PAIR_ssfr_B
				rtB['Lnuv_F']     = PAIR_L_nuv_B
			elif 'MRP' in Slct_Pair_ip_tbl and not 'PRP' in Slct_Pair_ip_tbl:
				rtB['re_F']        = PAIR_re_B
				rtB['re_err_F']    = PAIR_re_err_B
				rtB['n_F']         = PAIR_n_B
				rtB['n_err _F']    = PAIR_n_err_B
				rtB['q_F']         = PAIR_q_B
				rtB['q_err_F']     = PAIR_q_err_B
				rtB['theta_F']     = PAIR_theta_B
				rtB['theta_err_F'] = PAIR_theta_err_B
				rtB['chi_F']       = PAIR_chi_2_B
				rtB['flag_F']      = PAIR_flag_B

			elif 'PRP' in Slct_Pair_ip_tbl and 'MRP' in Slct_Pair_ip_tbl:
				rtB['mass_F']     = PAIR_MASS_B
				rtB['Age_F']      = PAIR_Age_B
				rtB['SFR_F']      = PAIR_sfr_B
				rtB['sSFR_F']     = PAIR_ssfr_B
				rtB['Lnuv_F']     = PAIR_L_nuv_B

				rtB['re_F']        = PAIR_re_B
				rtB['re_err_F']    = PAIR_re_err_B
				rtB['n_F']         = PAIR_n_B
				rtB['n_err _F']    = PAIR_n_err_B
				rtB['q_F']         = PAIR_q_B
				rtB['q_err_F']     = PAIR_q_err_B
				rtB['theta_F']     = PAIR_theta_B
				rtB['theta_err_F'] = PAIR_theta_err_B
				rtB['chi_F']       = PAIR_chi_2_B
				rtB['flag_F']      = PAIR_flag_B

				rtB['DIF_ALP']     = DIF_ALP_P
				rtB['DIF_DLT']     = DIF_DLT_P
				rtB['OMG1']        = OMG1_P
				rtB['OMG2']        = OMG2_P
				rtB['PHI']         = PHI_P
				rtB['PHI_ABS']     = PHIA_P
			else:
				pass

			rtB['DELTAZ']      = DELTAZ_P
			rtB['SEP_arcsec']  = SEP_AS_P
			rtB['arcsec/kpc']  = SCALEKPC_P
			rtB['SEP_kpc']     = SEP_KPC_P
			rtB['NEW_RA_B_F']  = NEW_RA_P
			rtB['NEW_DEC_B_F'] = NEW_DEC_P

			rtB['spc_f_B']     = spc_f[indexes_B]
			rtB['spc_f_n_B']   = spc_f_n[indexes_B]

			rtB['spc_f_F']     = PAIR_SFN_B
			rtB['spc_f_n_F']   = PAIR_SFN_N_B

			print
			rtB.write(op_tbl_B, format=tbl_format_opt,overwrite=True)
			print 'Results containing Foreground galaxies in table: ',op_tbl_B

			rtF.write(op_tbl_F, format=tbl_format_opt,overwrite=True)
			print 'Results containing Background galaxies in table: ',op_tbl_F
		elif grl_tbl == False:
			pass
	elif slt_prs == False:
		pass	

def lambdashifted(lambdao,z):
	lambdae = lambdao/(1+z)
	return lambdae

def q_2_incl(b_a_ratio):
	q = 0.22
	angle = math.degrees(math.acos(np.sqrt(((b_a_ratio**2) - (q**2)) /(1-(q**2)))))
	return angle

def fwhm2sigma(fwhm):
	return fwhm/(2*np.sqrt(2*np.log(2)))

def sigma2fwhm(sigma):
	return sigma*(2*np.sqrt(2*np.log(2)))

def astpy_conv_gaus_1dkernel(kernel):
	#http://docs.astropy.org/en/stable/api/astropy.convolution.Gaussian1DKernel.html#astropy.convolution.Gaussian1DKernel
	return apcvl.Gaussian1DKernel(stddev=kernel)

def split_list(alist, wanted_parts=1):
	length = len(alist)
	return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
			for i in range(wanted_parts) ]
def Lmbd_2_Vel(lambda_reference,lambda_observed,args,**kwargs):
	velocity = c_speed * (1- lambda_reference/lambda_observed)
	return velocity

def lw_sgma2fwhm(linewidth):
    return np.sqrt(8*np.log(2))*linewidth

def func_Linear(X,slope,b):
    return (slope*X) + b

def func_1D_Gaussian_Emm(X,X_0,A,SIGMA):
	if X_0 == 999999.99999 and A == 999999.99999 and SIGMA == 999999.99999:
		X_0   = 0
		A     = 0
		SIGMA = 0
		k     = 1
	else:
		pass
		k = 1
	return k + (A*np.exp(-(X-X_0)**2/(2*SIGMA**2)))

def func_1D_Gaussian(X,X_0,A,SIGMA):
	if X_0 == 999999.99999 and A == 999999.99999 and SIGMA == 999999.99999:
		X_0   = 0
		A     = 0
		SIGMA = 0
		k     = 1
	else:
		pass
		k = 1
	return k + (A*np.exp(-(X-X_0)**2/(2*SIGMA**2)))

def func_1D_Gaussian_O(X,X_0,A,SIGMA,OFFSET):
	if X_0 == 999999.99999 and A == 999999.99999 and SIGMA == 999999.99999 and OFFSET == 999999.99999:
		X_0    = 0
		A      = 0
		SIGMA  = 0
		OFFSET = 0
		k      = 1
	else:
		k = 1
	return OFFSET + k + (A*np.exp(-((X-X_0)**2)/(2*SIGMA**2)))

def func_2_1D_Gaussian(X,X_0_1,A_1,SIGMA_1,X_0_2,A_2,SIGMA_2):
    return 1 + (A_1*np.exp(-(X-X_0_1)**2/(2*SIGMA_1**2))) + (A_2*np.exp(-(X-X_0_2)**2/(2*SIGMA_2**2)))

def func_2_1D_Gaussian_O(X,X_0_1,A_1,SIGMA_1,X_0_2,A_2,SIGMA_2, OFFSET):
    return OFFSET + 1 + (A_1*np.exp(-(X-X_0_1)**2/(2*SIGMA_1**2))) + (A_2*np.exp(-(X-X_0_2)**2/(2*SIGMA_2**2)))

def func_Lorentzian(X,X_0,GAMMA):
    return 1 + (GAMMA / np.pi / ((X-X_0)**2 + GAMMA **2))

def func_Voigt(X,X_0,ALPHA,GAMMA):
    SIGMA = (ALPHA / np.sqrt(2 * np.log(2)))
    return 1 + np.real(wofz(((X-X_0)+1j*GAMMA) / SIGMA / np.sqrt(2))) / SIGMA / np.sqrt(2*np.pi)

def Cosmo_Scale_Par(z_cosmo_parm,Uni_Cosmo_Pars,*args, **kwargs):
	#Info obtained via:
	#http://www.astro.ucla.edu/~wright/CosmoCalc.html
	#For Ho = 70, OmegaM = 0.270, Omegavac = 0.730, z = 0.020
	#
	#It is now 13.861 Gyr since the Big Bang.
	#The age at redshift z was 13.585 Gyr.
	#The light travel time was 0.275 Gyr.
	#The comoving radial distance, which goes into Hubble's law, is 85.3 Mpc or 0.278 Gly.
	#The comoving volume within redshift z is 0.003 Gpc3.
	#The angular size distance DA is 83.634 Mpc or 0.272778 Gly.
	#This gives a scale of 0.405 kpc/".
	#The luminosity distance DL is 87.0 Mpc or 0.284 Gly.
	#1 Gly = 1,000,000,000 light years or 9.461*1026 cm.
	#1 Gyr = 1,000,000,000 years.
	#1 Mpc = 1,000,000 parsecs = 3.08568*1024 cm, or 3,261,566 light years.

	#astropy.cosmology.FLRW.arcsec_per_kpc_comoving
	#print cosmo.H(0)
	#print cosmo.Om(0)
	#print cosmo.Tcmb(0)
	#print cosmo.Neff
	#print cosmo.m_nu
	scale    = (1/Uni_Cosmo_Pars.arcsec_per_kpc_proper(z_cosmo_parm)) #kpc/arcsec
	scaleL   = Uni_Cosmo_Pars.luminosity_distance(z_cosmo_parm) #dlGrp Mpc
	c_speed  = physical_constants['speed of light in vacuum'][0]/1000
	DELTAZ   = rad_vel_sep/c_speed
	sepkpc   = rad_sep[1]*scale

	cosmo_check = kwargs.get('cosmo_check', False)

	if cosmo_check == True:

		print ''
		print 'Cosmology:'
		print
		print 'Cosmological parameters:'
		print Uni_Cosmo_Pars
		print 
		print 'At redshift                                      : ', z_cosmo_parm 
		print 'The angular separations [arc sec]                : ', rad_sep[1]
		print 'Angular distance scale                           : ', scale #,'[kpc/arcsec]'
		print 'The angular separations [kpc]                    : ', sepkpc.value
		print 'Luminosity distance(D_L)                         : ', scaleL#',[Mpc]'
		print ''
		print 'A radial velocity difference (DeltaV_r)          : ',rad_vel_sep
		print 'Corresponds to a redshift difference (Delta_z) of: ',DELTAZ,DELTAZ
		print 'c                                                : ',c_speed
		print
	elif cosmo_check == False:
		pass

	return scale.value,scaleL.value,sepkpc
####Fnc_Stk_Mth####