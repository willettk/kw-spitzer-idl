
;+
; NAME:
;       
;	OHM_SFR
;
; PURPOSE:
;
;	
;
; INPUTS:
;
;	
;
; OUTPUTS:
;
;	
;
; KEYWORDS:
;
;	
;
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Jan 10
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Neon SFR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ohm_neii = getlinelist('neII')
ohm_neii_tags = ohm_neii[0,*]
ohm_neii_flux = ohm_neii[1,*]

ohm_neiii = getlinelist('neIII')
ohm_neiii_tags = ohm_neiii[0,*]
ohm_neiii_flux = ohm_neiii[1,*]

n = n_elements(ohm_neii_tags)

obj = strarr(n)
dl = dblarr(n)

for i=0, n-1 do begin
	targets, ohm_neii_tags[i], r, o, d
	obj[i] = o
	dl[i] = d
endfor

n3 = n_elements(ohm_neiii_tags)

obj3 = strarr(n3)
dl3 = dblarr(n3)

for i=0, n3 - 1 do begin
	targets, ohm_neiii_tags[i], r, o, d
	obj3[i] = o
	dl3[i] = d
endfor

mpc2cm = 3.086d24
neii_lum  = 4d * !dpi * (dl * mpc2cm)^2 * ohm_neii_flux * 1d-21 * 1d7
neiii_lum = 4d * !dpi * (dl3 * mpc2cm)^2 * ohm_neiii_flux * 1d-21 * 1d7

match, obj, obj3, ne2ind, ne3ind

ohmobj = ohmdat('obj')
archobj = archdat('obj')

allobj = [transpose(ohmobj),transpose(archobj)]
allobj = allobj[sort(allobj)]

match, allobj, obj[ne2ind], a, b

emptyarr=replicate('         ',n_elements(allobj))
emptyarr[a] = string((neii_lum[ne2ind]+neiii_lum[ne3ind]) * 4.73d-41,format='(f9.0)')

; Non-masing galaxies

con_neii = getlinelist('neII',/con)
con_neii_tags = con_neii[0,*]
con_neii_flux = con_neii[1,*]

con_neiii = getlinelist('neIII',/con)
con_neiii_tags = con_neiii[0,*]
con_neiii_flux = con_neiii[1,*]

ncon = n_elements(con_neii_tags)

obj = strarr(ncon)
dl = dblarr(ncon)

for i=0, ncon-1 do begin
	targets, con_neii_tags[i], r, o, d
	obj[i] = o
	dl[i] = d
endfor

n3 = n_elements(con_neiii_tags)

obj3 = strarr(n3)
dl3 = dblarr(n3)

for i=0, n3 - 1 do begin
	targets, con_neiii_tags[i], r, o, d
	obj3[i] = o
	dl3[i] = d
endfor

conneii_lum  = 4d * !dpi * (dl * mpc2cm)^2 * con_neii_flux * 1d-21 * 1d7
conneiii_lum = 4d * !dpi * (dl3 * mpc2cm)^2 * con_neiii_flux * 1d-21 * 1d7

match, obj, obj3, ne2ind, ne3ind

conobj = condat('obj')
archobj = archdat('obj')

conallobj = [transpose(conobj)]
conallobj = conallobj[sort(conallobj)]

match, conallobj, obj[ne2ind], a, b

conemptyarr=replicate('         ',n_elements(conallobj))
conemptyarr[a] = string((neii_lum[ne2ind]+neiii_lum[ne3ind]) * 4.73d-41,format='(f9.0)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PAH SFR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ohm_pah62flux = ohmdat('pah62flux',/obj)
arch_pah62flux = archdat('pah62flux',/obj)

all_pah62flux = [transpose(ohm_pah62flux[1,*]),transpose(arch_pah62flux[1,*])]
all_pah62obj = [transpose(ohm_pah62flux[0,*]),transpose(arch_pah62flux[0,*])]

all_pah62flux = all_pah62flux[sort(all_pah62obj)]
all_pah62obj = all_pah62obj[sort(all_pah62obj)]

ohm_pah11flux = ohmdat('pah11flux',/obj)
arch_pah11flux = archdat('pah11flux',/obj)

all_pah11flux = [transpose(ohm_pah11flux[1,*]),transpose(arch_pah11flux[1,*])]
all_pah11obj = [transpose(ohm_pah11flux[0,*]),transpose(arch_pah11flux[0,*])]

all_pah11flux = all_pah11flux[sort(all_pah11obj)]
all_pah11obj = all_pah11obj[sort(all_pah11obj)]

match, all_pah62obj, all_pah11obj, a4, b4

allpahflux = all_pah62flux[a4] + all_pah11flux[b4]
allpahobj = all_pah62obj[a4]

match, allobj, allpahobj, aa, bb

dl_ohm = ohmdat('dl',/obj)
dl_arch = archdat('dl',/obj)

dl_all = [transpose(dl_ohm[1,*]),transpose(dl_arch[1,*])]
dl_all_obj = [transpose(dl_ohm[0,*]),transpose(dl_arch[0,*])]

match, allpahobj, dl_all_obj, a3, b3

dl_matched = float(dl_all[b3])

pah62lum  = 4d * !dpi * (dl_matched * mpc2cm)^2 * all_pah62flux *  1d7

emptyarr2=replicate('         ',n_elements(allobj))
emptyarr2[a3] = string((pah62lum) * 1.18d-41,format='(f9.0)')

; Non-masing galaxies

con_pah62flux = condat('pah62flux',/obj)

all_pah62flux = [transpose(con_pah62flux[1,*])]
all_pah62obj = [transpose(con_pah62flux[0,*])]

all_pah62flux = all_pah62flux[sort(all_pah62obj)]
all_pah62obj = all_pah62obj[sort(all_pah62obj)]

con_pah11flux = condat('pah11flux',/obj)

all_pah11flux = [transpose(con_pah11flux[1,*])]
all_pah11obj = [transpose(con_pah11flux[0,*])]

all_pah11flux = all_pah11flux[sort(all_pah11obj)]
all_pah11obj = all_pah11obj[sort(all_pah11obj)]

match, all_pah62obj, all_pah11obj, a4, b4

allpahflux = all_pah62flux[a4] + all_pah11flux[b4]
conallpahobj = all_pah62obj[a4]

match, allobj, conallpahobj, aa, bb

dl_con = condat('dl',/obj)

dl_all = [transpose(dl_con[1,*])]
dl_all_obj = [transpose(dl_con[0,*])]

match, conallpahobj, dl_all_obj, a3, b3

dl_matched = float(dl_all[b3])

pah62lum  = 4d * !dpi * (dl_matched * mpc2cm)^2 * all_pah62flux *  1d7

conemptyarr2=replicate('         ',n_elements(conallpahobj))
conemptyarr2[a3] = string((pah62lum) * 1.18d-41,format='(f9.0)')

print,[transpose(allobj), transpose(emptyarr), transpose(emptyarr2)]
print,''
print,[transpose(conallpahobj), transpose(conemptyarr), transpose(conemptyarr2)]

end
