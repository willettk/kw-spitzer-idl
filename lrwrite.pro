; Quick routine to write the individual, calibrated fluxes into wavelength-sorted files for each target. -KW, Sep 07


;restore,'~/Astronomy/Research/Spitzer/OHM/fname.sav'
restore,'~/Astronomy/Research/Spitzer/control/conname.sav'

fname = conname
nfiles = n_elements(fname)

for i = 0, nfiles - 1 do begin

	;dir = '~/Astronomy/Research/Spitzer/OHM/data/idl_spectra/calibrated/coadd/'
	dir = '~/Astronomy/Research/Spitzer/control/data/idl_spectra/calibrated/coadd/'

	readcol, dir+fname(i)+'_ll1_cal.tbl', $
		ll1_order, ll1_wave, ll1_flux, ll1_err, ll1_bit, format = 'i,f,f,f,i',skipline=1,/silent

	readcol, dir+fname(i)+'_ll2_cal.tbl', $
		ll2_order, ll2_wave, ll2_flux, ll2_err, ll2_bit, format = 'i,f,f,f,i',skipline=1,/silent

	readcol, dir+fname(i)+'_ll3_cal.tbl', $
		ll3_order, ll3_wave, ll3_flux, ll3_err, ll3_bit, format = 'i,f,f,f,i',skipline=1,/silent

	readcol, dir+fname(i)+'_sl1_cal.tbl', $
		sl1_order, sl1_wave, sl1_flux, sl1_err, sl1_bit, format = 'i,f,f,f,i',skipline=1,/silent

	readcol, dir+fname(i)+'_sl2_cal.tbl', $
		sl2_order, sl2_wave, sl2_flux, sl2_err, sl2_bit, format = 'i,f,f,f,i',skipline=1,/silent

	readcol, dir+fname(i)+'_sl3_cal.tbl', $
		sl3_order, sl3_wave, sl3_flux, sl3_err, sl3_bit, format = 'i,f,f,f,i',skipline=1,/silent

	wavesort = sort([ll1_wave,ll2_wave,ll3_wave,sl1_wave,sl2_wave,sl3_wave])

	wave1 = [ll1_wave,ll2_wave,ll3_wave,sl1_wave,sl2_wave,sl3_wave]
	flux1 = [ll1_flux,ll2_flux,ll3_flux,sl1_flux,sl2_flux,sl3_flux]
	order1 = [ll1_order,ll2_order,ll3_order,sl1_order,sl2_order,sl3_order]
	bit1 = [ll1_bit,ll2_bit,ll3_bit,sl1_bit,sl2_bit,sl3_bit]
	err1 = [ll1_err,ll2_err,ll3_err,sl1_err,sl2_err,sl3_err]

	wave = wave1(wavesort)
	flux = flux1(wavesort)
	order = order1(wavesort)
	bit = bit1(wavesort)
	err = err1(wavesort)

	; Write out the sorted, combined data

	forprint, order, wave, flux, err, bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]      Uncert.        Bit type', $
		textout = dir+fname(i)+'_lr_cal.tbl', /silent


endfor



end
