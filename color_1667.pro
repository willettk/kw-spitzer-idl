pro color_1667, ps = ps
;+
; NAME:
;       COLOR_1667
;
; PURPOSE:
; 	Plot optical depth as a function of dust temperature, coded by the measured optical
;	depth in the 1667 OH maser line.
;
; INPUTS:
;
; KEYWORDS:
;
; EXAMPLE:
;	IDL> color_1667
;
; REQUIRES:
;
; NOTES:
;
; REVISION HISTORY
;       Written by K. Willett                Dec 2007
;-

; Load dust temperatures

temp1 = ohmdat('dtemp',sz=2)

dtemp = temp1(*,0)
dtemp_err = temp1(*,1)

; Load tau_v

temp2 = ohmdat('sil', sz=2)

tau_v = -18.5 * temp2(*,0)		; Conversion from tau_sil to tau_V is from Roche & Aiken (1984)
tau_v_err = -18.5 * temp2(*,1) 

; Load maser optical depth

tabledir = '~/Astronomy/Research/Spitzer/OHM/tables/'
newohm_tbl = tabledir+'newohm_radio.tbl'
oldohm_tbl = tabledir+'oldohm_radio.tbl'

readcol, newohm_tbl, new_obj, new_cz, new_f1667, new_w1667, new_delnu1667, new_delv1667, new_rh, new_lfir, new_loh, new_f1420, $
	format = 'a,a,f,f,f,f,i', /silent

readcol, oldohm_tbl, old_obj, old_z, old_lfir, old_loh_pred, old_loh, old_f1667, old_ref, old_f1420, $
	format = 'a,f,f,f,f,f,i,f', /silent

; Find the OHMs that have both 1667 and 1420 MHz flux density measurements

names = ohmdat('obj')
names = strmid(names,5,10)
nobj = n_elements(names)

names_tbl = [new_obj,old_obj]
f1667 = [new_f1667,old_f1667]
f1420 = [new_f1420,old_f1420]

index = intarr(nobj)
for i = 0, nobj-1 do begin
	index(i) = where(names(i) eq names_tbl)
	if index(i) eq -1 then print, 'No match found for '+names(i)
endfor

tau_ohm = -alog(f1667(index)/f1420(index)+1d)

; Plot results

plotdir = '~/Astronomy/Research/Spitzer/OHM/plots/'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotdir+'color_1667.ps'
	lthick = 2
	cthick = 2
endif else begin
	lthick = 1
	cthick = 1
endelse

!p.multi=[0,1,1]
plot, dtemp, tau_v, /nodata, $
	xtitle = 'T!Id!N [K]', $
	ytitle = '!7s!3!IV!N', $
	title = '1667 OHM optical depth vs. dust temperature and opacity', $
	charsize = 1.5, $
	xrange = [30,180], /xstyle, $
	yrange = [0,150], /ystyle, $
	/isotropic, $
	background = fsc_color("White"), $
	color=fsc_color("Black"), $
	charthick = cthick, $
	thick = lthick

; Recreate contours from Lockett and Elitzur

yescont = 1
if yescont eq 1 then begin
cont_10_x = [48,51,59,63,80,180,180,80,63,58,54,52,48]
cont_10_y = [150,35,16,8,5,6,7,6,12,30,60,150,150]
polyfill,cont_10_x,cont_10_y,color=fsc_color("Orange")
cont_10_x = [170,180,180,170]
cont_10_y = [150,150,130,150]
polyfill,cont_10_x,cont_10_y,color=fsc_color("Orange")
cont_15_x = [170,180,180,155,122,170]
cont_15_y = [150,130,68,90,150,150]
polyfill,cont_15_x,cont_15_y,color=fsc_color("Yellow")
cont_15_x = [51,52,59,63,80,179,179,80,63,60,56,55,51]+1
cont_15_y = [149,35,16,8,5,6,7,6,12,30,60,149,149]+1
polyfill,cont_15_x,cont_15_y,color=fsc_color("Yellow")

cont_20_x = [54, 54,55,58,59,63,80,178,178,70,62,61, 61, 52]+2
cont_20_y = [148,80,60,35,16,8, 5,   6,  8,20,30,60,148,148]+2
polyfill,cont_20_x,cont_20_y,color=fsc_color("Green")

cont_20_x = [122,155,180,180,153,135,121,93,122]
cont_20_y = [150,90,68,40,51,69,94,150,150]
polyfill,cont_20_x,cont_20_y,color=fsc_color("Green")
cont_25_x = [93,121,135,153,180,180,130,80,74,68,63,62,63,93]
cont_25_y = [150,94,69,51,40,10,9,10,11,15,35,61,150,150]
polyfill,cont_25_x,cont_25_y,color=fsc_color("Pale Green")
cont_30_x = [78,85,102,122,145,180,180,130,100,80,75,67,78]
cont_30_y = [109,104,80,53,32,24,12,11,11,12,13,51,109]
polyfill,cont_30_x,cont_30_y,color=fsc_color("Cyan")
cont_35_x = [84,93,103,117,160,130,97,81,79,84]
cont_35_y = [53,54,48,31,18,13,14,17,30,53]
polyfill,cont_35_x,cont_35_y,color=fsc_color("Blue")
endif

tau_blue = where(tau_ohm lt -3.0)
tau_cyan = where(tau_ohm gt -3.0 and tau_ohm lt -2.5)
tau_pgreen = where(tau_ohm gt -2.5 and tau_ohm lt -2.0)
tau_green = where(tau_ohm gt -2.0 and tau_ohm lt -1.5)
tau_yellow = where(tau_ohm gt -1.5 and tau_ohm lt -1.0)
tau_orange = where(tau_ohm gt -1.0 and tau_ohm lt -0.5)
tau_white = where(tau_ohm gt -0.5)

if tau_blue(0) ne -1 then begin
	oplot, dtemp(tau_blue), tau_v(tau_blue), psym = symcat(16), color=fsc_color("Blue"),symsize=2
	oplot, dtemp(tau_blue), tau_v(tau_blue), psym = symcat(9), color=fsc_color("Black"),symsize=2
endif
if tau_cyan(0) ne -1 then begin
	oplot, dtemp(tau_cyan), tau_v(tau_cyan), psym = symcat(16), color=fsc_color("Cyan"),symsize=2
	oplot, dtemp(tau_cyan), tau_v(tau_cyan), psym = symcat(9), color=fsc_color("Black"),symsize=2
endif
if tau_pgreen(0) ne -1 then begin
	oplot, dtemp(tau_pgreen), tau_v(tau_pgreen), psym = symcat(16), color=fsc_color("Pale Green"),symsize=2
	oplot, dtemp(tau_pgreen), tau_v(tau_pgreen), psym = symcat(9), color=fsc_color("Black"),symsize=2
endif
if tau_green(0) ne -1 then begin
	oplot, dtemp(tau_green), tau_v(tau_green), psym = symcat(16), color=fsc_color("Green"),symsize=2
	oplot, dtemp(tau_green), tau_v(tau_green), psym = symcat(9), color=fsc_color("Black"),symsize=2
endif
if tau_yellow(0) ne -1 then begin
	oplot, dtemp(tau_yellow), tau_v(tau_yellow), psym = symcat(16), color=fsc_color("Yellow"),symsize=2
	oplot, dtemp(tau_yellow), tau_v(tau_yellow), psym = symcat(9), color=fsc_color("Black"),symsize=2
endif
if tau_orange(0) ne -1 then begin
	oplot, dtemp(tau_orange), tau_v(tau_orange), psym = symcat(16), color=fsc_color("Orange"),symsize=2
	oplot, dtemp(tau_orange), tau_v(tau_orange), psym = symcat(9), color=fsc_color("Black"),symsize=2
endif
if tau_white(0) ne -1 then begin
	oplot, dtemp(tau_white), tau_v(tau_white), psym = symcat(16), color=fsc_color("White"),symsize=2
	oplot, dtemp(tau_white), tau_v(tau_white), psym = symcat(9), color=fsc_color("Black"),symsize=2
endif

;nc=8
;loadct,4,ncolors=nc,/silent
;colorbar, ncolors=nc, /vertical, position=[0.93, 0.40, 0.98, 0.90], color=fsc_color("Black"), $
;	ticknames=string(findgen(8)*(-0.5),format='(f4.1)')
y1 = 0.899
y2 = 0.817
x1 = 0.931
x2 = 0.98
;polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=fsc_color("Blue"),/normal

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

stop
end
