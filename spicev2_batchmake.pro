; pro spicev2_batchmake

con_names = ['control004','control008','control013', $			; Objects in first version of CS (Darling)
	'control021','control022','control023','control024', $		; Baan
	'control025','control026','control027','control028', $		; Baan
	'control029','control030', $					; Baan
	'control031','control032','control033','control034', $		; Darling
	'control035','control036','control037','control038']		; Darling

n = n_elements(con_names)
path = '/Users/willettk/Astronomy/Research/Spitzer/spitzmed/optimal/'

openw, lun1, '~/Astronomy/Research/Spitzer/control/spicev2_con_lr.batch', /append, /get_lun
openw, lun2, '~/Astronomy/Research/Spitzer/control/spicev2_con_hr.batch', /append, /get_lun

for i = 0, n-1 do begin
	ch0_1o_1p_bcd  = path+con_names(i)+'_ch0_sub_medarr_1o_1p_bcd_clean.fits'
	ch0_1o_1p_mask = path+con_names(i)+'_ch0_sub_medarr_1o_1p_bmask_clean.fits'
	ch0_1o_1p_func = path+con_names(i)+'_ch0_sub_medarr_1o_1p_func_clean.fits'
	ch0_1o_2p_bcd  = path+con_names(i)+'_ch0_sub_medarr_1o_2p_bcd_clean.fits'
	ch0_1o_2p_mask = path+con_names(i)+'_ch0_sub_medarr_1o_2p_bmask_clean.fits'
	ch0_1o_2p_func = path+con_names(i)+'_ch0_sub_medarr_1o_2p_func_clean.fits'
	ch0_2o_1p_bcd  = path+con_names(i)+'_ch0_sub_medarr_2o_1p_bcd_clean.fits'
	ch0_2o_1p_mask = path+con_names(i)+'_ch0_sub_medarr_2o_1p_bmask_clean.fits'
	ch0_2o_1p_func = path+con_names(i)+'_ch0_sub_medarr_2o_1p_func_clean.fits'
	ch0_2o_2p_bcd  = path+con_names(i)+'_ch0_sub_medarr_2o_2p_bcd_clean.fits'
	ch0_2o_2p_mask = path+con_names(i)+'_ch0_sub_medarr_2o_2p_bmask_clean.fits'
	ch0_2o_2p_func = path+con_names(i)+'_ch0_sub_medarr_2o_2p_func_clean.fits'
	ch2_1o_1p_bcd  = path+con_names(i)+'_ch2_sub_medarr_1o_1p_bcd_clean.fits'
	ch2_1o_1p_mask = path+con_names(i)+'_ch2_sub_medarr_1o_1p_bmask_clean.fits'
	ch2_1o_1p_func = path+con_names(i)+'_ch2_sub_medarr_1o_1p_func_clean.fits'
	ch2_1o_2p_bcd  = path+con_names(i)+'_ch2_sub_medarr_1o_2p_bcd_clean.fits'
	ch2_1o_2p_mask = path+con_names(i)+'_ch2_sub_medarr_1o_2p_bmask_clean.fits'
	ch2_1o_2p_func = path+con_names(i)+'_ch2_sub_medarr_1o_2p_func_clean.fits'
	ch2_2o_1p_bcd  = path+con_names(i)+'_ch2_sub_medarr_2o_1p_bcd_clean.fits'
	ch2_2o_1p_mask = path+con_names(i)+'_ch2_sub_medarr_2o_1p_bmask_clean.fits'
	ch2_2o_1p_func = path+con_names(i)+'_ch2_sub_medarr_2o_1p_func_clean.fits'
	ch2_2o_2p_bcd  = path+con_names(i)+'_ch2_sub_medarr_2o_2p_bcd_clean.fits'
	ch2_2o_2p_mask = path+con_names(i)+'_ch2_sub_medarr_2o_2p_bmask_clean.fits'
	ch2_2o_2p_func = path+con_names(i)+'_ch2_sub_medarr_2o_2p_func_clean.fits'

	ch1_1p_bcd  = path+con_names(i)+'_ch1_medarr_1p_bcd_clean.fits'
	ch1_1p_mask = path+con_names(i)+'_ch1_medarr_1p_bmask_clean.fits'
	ch1_1p_func = path+con_names(i)+'_ch1_medarr_1p_func_clean.fits'
	ch1_2p_bcd  = path+con_names(i)+'_ch1_medarr_2p_bcd_clean.fits'
	ch1_2p_mask = path+con_names(i)+'_ch1_medarr_2p_bmask_clean.fits'
	ch1_2p_func = path+con_names(i)+'_ch1_medarr_2p_func_clean.fits'
	ch3_1p_bcd  = path+con_names(i)+'_ch3_sub_medarr_1p_bcd_clean.fits'
	ch3_1p_mask = path+con_names(i)+'_ch3_sub_medarr_1p_bmask_clean.fits'
	ch3_1p_func = path+con_names(i)+'_ch3_sub_medarr_1p_func_clean.fits'
	ch3_2p_bcd  = path+con_names(i)+'_ch3_sub_medarr_2p_bcd_clean.fits'
	ch3_2p_mask = path+con_names(i)+'_ch3_sub_medarr_2p_bmask_clean.fits'
	ch3_2p_func = path+con_names(i)+'_ch3_sub_medarr_2p_func_clean.fits'

	printf, lun1, ch0_1o_1p_bcd+', '+ch0_1o_1p_mask+', '+ch0_1o_1p_func
	printf, lun1, ch0_1o_2p_bcd+', '+ch0_1o_2p_mask+', '+ch0_1o_2p_func
	printf, lun1, ch0_2o_1p_bcd+', '+ch0_2o_1p_mask+', '+ch0_2o_1p_func
	printf, lun1, ch0_2o_2p_bcd+', '+ch0_2o_2p_mask+', '+ch0_2o_2p_func
	printf, lun1, ch2_1o_1p_bcd+', '+ch2_1o_1p_mask+', '+ch2_1o_1p_func
	printf, lun1, ch2_1o_2p_bcd+', '+ch2_1o_2p_mask+', '+ch2_1o_2p_func
	printf, lun1, ch2_2o_1p_bcd+', '+ch2_2o_1p_mask+', '+ch2_2o_1p_func
	printf, lun1, ch2_2o_2p_bcd+', '+ch2_2o_2p_mask+', '+ch2_2o_2p_func

	printf, lun2, ch1_1p_bcd+', '+ch1_1p_mask+', '+ch1_1p_func
	printf, lun2, ch1_2p_bcd+', '+ch1_2p_mask+', '+ch1_2p_func
	printf, lun2, ch3_1p_bcd+', '+ch3_1p_mask+', '+ch3_1p_func
	printf, lun2, ch3_2p_bcd+', '+ch3_2p_mask+', '+ch3_2p_func

endfor

close, lun1
close, lun2

end
