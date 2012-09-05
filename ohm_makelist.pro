; Program to read in all known OHMs (both the Arecibo survey and previous detections)
; and create a target list to be read using Leopard

path = '~/Astronomy/Research/Spitzer/OHM/tables/'
known = path+'known_ohm.txt'
arecibo = path+'arecibo_ohm.txt'

readcol, known, $
	name_known, ra_h, ra_m, ra_s, dec_h, dec_m, dec_s, z_hel, ref, v_hel, v_cmb, dl, f60, f100, fir, $
	format='a,i,i,f,i,i,f,f,i,a,a,a,a,a,a', skipline=1, /silent

readcol, arecibo, $
	name_arecibo, ra_h, ra_m, ra_s, dec_h, dec_m, dec_s, z_hel, ref, v_hel, v_cmb, dl, f60, f100, fir, $
	format='a,i,i,f,i,i,f,f,i,a,a,a,a,a,a', skipline=1, /silent

list = [[transpose(name_known)],[transpose(name_arecibo)]]


openw, 2, path+'all_ohmnames.txt'
printf, 2, transpose('IRAS'+list(sort(list)))
close, 2



end
