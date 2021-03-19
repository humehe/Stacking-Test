from astropy.io import fits

####Fnc_Stk_Fts####
def Wrt_FITS_File(img_inpt_array,otp_img_fn,*args, **kwargs):
	hdu        = fits.PrimaryHDU(img_inpt_array)
	hdulist    = fits.HDUList([hdu])
	hdulist.writeto(otp_img_fn,overwrite=True)
	hdulist.close()
	return otp_img_fn
	
def Header_Get(image,field):
	hdulist = fits.open(image)
	try:
		header  = float(hdulist[0].header[field])
	except:
		header  = str(hdulist[0].header[field])
	hdulist.close()
	return header

def Header_Updt(fitsfile_in,field,value,*args, **kwargs):
	hdulist_afh  = fits.open(fitsfile_in, mode='update')
	prihdr       = hdulist_afh[0].header
	prihdr.set(field, value)
	hdulist_afh.flush()
	hdulist_afh.close()

def Header_Add(fitsfile_in,field,value,*args, **kwargs):
	header_comment = kwargs.get('header_comment',None)
	hdulist_afh  = fits.open(fitsfile_in, mode='update')
	prihdr       = hdulist_afh[0].header
	prihdr.set(field, value,comment=header_comment)
	hdulist_afh.flush()
	hdulist_afh.close()	

def Header_Get_Add(fits_ipfn,header,value,*args, **kwargs):
	header_comment = kwargs.get('header_comment',None)
	try:
		head_val      = (Header_Get(fits_ipfn,header))
	except KeyError:
		Header_Add(fits_ipfn,header,value,header_comment=header_comment)
		head_val      = (Header_Get(fits_ipfn,header))
	return head_val	

def Header_Copy(fits_to,fits_from,header,*args,**kwargs):
	header_comment = kwargs.get('header_comment',None)
	value = Header_Get(fits_from,header)
	Header_Get_Add(fits_to,header,value,header_comment=header_comment)

def Header_History_Step(fits_ipfn_prv,fits_ipfn_pst):
	try:
		head_val      = int(Header_Get(fits_ipfn_prv,'h_s_c'))
	except KeyError:
		head_val = 0
	c_step = int(head_val + 1)
	header = 'h_s_'+str(c_step)
	Header_Add(fits_ipfn_pst,header,str((fits_ipfn_pst.rsplit('/',1)[1]).rsplit('.',1)[0]),header_comment='History Step')
	Header_Updt(fits_ipfn_pst,'h_s_c',c_step)

def Header_History_Step_Backwards_Loop(fits_ipfn_crt,fits_ipfn_prv):
	loop_back_step = int(Header_Get(fits_ipfn_crt,'h_s_c'))
	while (loop_back_step > 0):
		loop_back_step = int(loop_back_step - 1)
		loop_back_step_header = 'h_s_' + str(loop_back_step)
		Header_Copy(fits_ipfn_crt,fits_ipfn_prv,loop_back_step_header)
	else:
		pass
####Fnc_Stk_Fts####