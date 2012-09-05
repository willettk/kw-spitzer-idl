;+
; NAME:
;       
;	OHWRITE
;
; PURPOSE:
;
;	Lookup table for data on OH megamasers and control sample objects
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
;	IDL> .r ohwrite
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett             Mar 08
; 	Modified to account for new luminosity distances (assumed H_) = 70.5, WMAP5 cosmology) - Feb 09
; 	Added hyperfine ratios (R_H) from OHMCAT - Mar 09
;	Reformatted mega OHMs in table, added heliocentric optical velocities and errors from Darling papers and literature - Apr 09
;-


;####################
; OHM (Darling sample)
;####################

nohm = n_elements(ohmdat('tag'))
meganames = ohmdat('tag')

; OHM data

;					log OH (H75)	f1667	f1420		cz_oh	czoh_err    zerr_opt cz_opt	cz_opt_err
;
megaarr = [$
['mega001', 'IRAS 01562+2528'        ,'3.25'       ,'6.95 '     , '6.3 '  ,   '49814'  ,   '15'  , '14 '    , '49665' ,      '45 ' ],$
['mega002', 'IRAS 02524+2046'        ,'3.74'       ,'39.82'     , '2.9 '  ,   '54162'  ,   '15'  , '18 '    , '54383' ,      '55 ' ],$
['mega004', 'IRAS 04121+0223'        ,'2.32'       ,'2.52 '     , '3.1 '  ,   '36590'  ,   '14'  , '11 '    , '36702' ,      '34 ' ],$
['mega005', 'IRAS 06487+2208'        ,'2.82'       ,'7.60 '     , '10.8'  ,   '42972'  ,   '14'  , '8  '    , '42972' ,      '26 ' ],$
['mega006', 'IRAS 07163+0817'        ,'2.37'       ,'4.00 '     , '3.5 '  ,   '33150'  ,   '14'  , '22 '    , '33269' ,      '67 ' ],$
['mega007', 'IRAS 07572+0533'        ,'2.74'       ,'2.26 '     , '5.0 '  ,   '56845'  ,   '15'  , '29 '    , '56900' ,      '90 ' ],$
['mega008', 'IRAS 08201+2801'        ,'3.45'       ,'14.67'     , '16.7'  ,   '50325'  ,   '15'  , '8  '    , '50272' ,      '26 ' ],$
['mega009', 'IRAS 08449+2332'        ,'2.59'       ,'2.49 '     , '6.1 '  ,   '45424'  ,   '14'  , '10 '    , '45406' ,      '31 ' ],$
['mega010', 'IRAS 08474+1813'        ,'2.70'       ,'2.20 '     , '4.2 '  ,   '43750'  ,   '14'  , '17 '    , '43591' ,      '53 ' ],$
['mega013', 'IRAS 10035+2740'        ,'2.50'       ,'2.29 '     , '6.3 '  ,   '50065'  ,   '14'  , '100'    , '49826' ,      '300' ],$
['mega014', 'IRAS 10339+1548'        ,'2.65'       ,'6.26 '     , '5.1 '  ,   '58983'  ,   '15'  , '6  '    , '59130' ,      '18 ' ],$
['mega016', 'IRAS 11028+3130'        ,'2.97'       ,'4.27 '     , '5.0 '  ,   '59619'  ,   '15'  , '9  '    , '59540' ,      '27 ' ],$
['mega017', 'IRAS 11180+1623'        ,'2.34'       ,'1.82 '     , '4.2 '  ,   '49783'  ,   '14'  , '23 '    , '49766' ,      '70 ' ],$
['mega018', 'IRAS 11524+1058'        ,'2.98'       ,'3.17 '     , '5.0 '  ,   '53404'  ,   '15'  , '29 '    , '53564' ,      '90 ' ],$
['mega020', 'IRAS 14059+2000'        ,'3.34'       ,'15.20'     , '7.5 '  ,   '37246'  ,   '14'  , '30 '    , '37084' ,      '89 ' ],$
['mega022', 'IRAS 14553+1245'        ,'2.27'       ,'2.93 '     , '3.8 '  ,   '37462'  ,   '14'  , '35 '    , '37449' ,      '133' ],$
['mega023', 'IRAS 16255+2801'        ,'2.57'       ,'7.02 '     , '5.0 '  ,   '40076'  ,   '14'  , '7  '    , '40056' ,      '24 ' ],$
;['mega025', 'IRAS 17540+2935'        ,'1.76'       ,'0.76 '     , '4.0 '  ,   '32522'  ,   '14'  , '34 '    , '32410' ,      '104' ],$
['mega026', 'IRAS 18368+3549'        ,'2.85'       ,'4.58 '     , '21.0'  ,   '34832'  ,   '14'  , '48 '    , '34813' ,      '146' ],$
['mega027', 'IRAS 18588+3517'        ,'2.52'       ,'7.37 '     , '5.9 '  ,   '31686'  ,   '14'  , '9  '    , '31996' ,      '28 ' ],$
['mega028', 'IRAS 20286+1846'        ,'3.41'       ,'15.58'     , '5.0 '  ,   '40471'  ,   '14'  , '10 '    , '40700' ,      '32 ' ],$
['mega029', 'IRAS 21077+3358'        ,'3.26'       ,'5.04 '     , '9.4 '  ,   '52987'  ,   '15'  , '15 '    , '52973' ,      '48 ' ],$
['mega031', 'IRAS 22055+3024'        ,'2.73'       ,'6.35 '     , '6.4 '  ,   '37965'  ,   '14'  , '9  '    , '38071' ,      '27 ' ],$
['mega032', 'IRAS 22116+0437'        ,'2.77'       ,'1.76 '     , '8.4 '  ,   '58180'  ,   '15'  , '16 '    , '58098' ,      '49 ' ],$
;['mega033', 'IRAS 23019+3405'        ,'2.12'       ,'3.58 '     , '7.7 '  ,   '32294'  ,   '14'  , '23 '    , '32360' ,      '69 ' ],$
['mega034', 'IRAS 23028+0725'        ,'3.29'       ,'8.69 '     , '19.5'  ,   '44529'  ,   '14'  , '28 '    , '44784' ,      '86 ' ]]


; Convert H_0 = 75 to H_0 = 70

readcol,'~/Astronomy/Research/Spitzer/dl_h70.txt',delim='&',format='a,i,i', /silent, $
	dl_obj, d75, d70

mega_taglist = ohmdat('tag')
obj_mega = ohmdat('obj')
obj_mega = strmid(obj_mega,1,strlen(obj_mega[0]))
logoh = float(megaarr[2,*])

meganames = megaarr[0,*]

for i = 0, n_elements(meganames) - 1 do begin

	match, meganames[i], mega_taglist, m1, m2
	match, dl_obj, obj_mega[m2], a, b
	oldoh = logoh[i]
	newoh = oldoh + 2d * alog10(double(d70[a]) / double(d75[a]))
	logoh[i] = newoh

endfor

; Hyperfine ratios

mega_rh = strarr(n_elements(meganames))
for i = 0, n_elements(meganames) - 1 do begin

	targets, meganames[i], z, obj, dl
	a = ohmcat(strtrim(obj,2))
	rh_temp = a.rh
	rhm_temp = a.rh_marker
	if rhm_temp eq '' then mega_rh[i] = rh_temp else mega_rh[i] = 0.
	
endfor


; Peak fluxes of 1667 MHz line and 1.4 GHz continuum [mJy]

f1667_ohm = float(megaarr[3,*])
f1420_ohm = float(megaarr[4,*])

; Velocity and error for OHM emission

czoh_ohm = float(megaarr[5,*])
czoh_err_ohm = float(megaarr[6,*])

; Optical velocities

zerr_opt_ohm = 1d-5 * float(megaarr[7,*])

cz_opt_ohm = float(megaarr[8,*])
cz_opt_err_ohm = float(megaarr[9,*])

;####################
; Archived OHMs
;####################

narch = n_elements(archdat('tag'))

; 50 km/s is default error for objects w/o defined uncertainties in cz_oh

; cz_opt_err uses Darling & Giovanelli (2006) optical velocities where available, literature values compiled in
;	Darling and Giovanelli Paper III where not measured with Palomar. -KW, Apr 09

; arch025 is from Dickey et al (1990)

;					log OH		f1667	f1420		cz_oh	czoh_err    zerr_opt cz_opt	cz_opt_err
;
archarr = [$
['arch003', 'IRAS 01418+1651 '        ,'2.65'       ,'240.'      ,'40.6' ,  '8260',      '30',     ''    ,  '8245 ',    '5  '     ],$
['arch004', 'IRAS 03521+0028 '        ,'2.44'       ,'2.77'      ,'6.7'  , '45512',      '14',     '11'  ,  '45585',    '35 '     ],$
['arch005', 'IRAS 04454-4838 '        ,'2.88'       ,'142.'      ,'0.'   , '15920',      '50',     ''    ,  '15682',    '30 '     ],$
['arch007', 'IRAS 09039+0503 '        ,'2.80'       ,'5.17'      ,'6.6'  , '37720',      '14',     '11'  ,  '37516',    '34 '     ],$
;['arch008', 'IRAS 09320+6134 '        ,'1.61'       ,'12.'       ,'170.9',      '',      '',       ''    ,  '11809',    '9  '     ],$
['arch009', 'IRAS 09539+0857 '        ,'3.45'       ,'14.32'     ,'9.5'  , '38455',      '14',     '15'  ,  '38643',    '15 '     ],$
['arch010', 'IRAS 10039-3338 '        ,'2.92'       ,'315.'      ,'24.7' , '10090',      '50',     ''    ,  '10223',    '36 '     ],$
['arch012', 'IRAS 10173+0828 '        ,'2.68'       ,'105.'      ,'10.8' , '14720',      '39',     ''    ,  '14390',    '300'     ],$
['arch013', 'IRAS 10378+1109 '        ,'3.27'       ,'13.'       ,'8.9'  , '40811',      '14',     '9'   ,  '40854',    '29 '     ],$
['arch014', 'IRAS 10485-1447 '        ,'2.91'       ,'0.'        ,'4.4'  ,      '',      '',       ''    ,  '39872',    '70 '     ],$
['arch016', 'IRAS 12018+1941 '        ,'2.87'       ,'3.'        ,'6.5'  , '50350',      '21',     ''    ,  '50559',    '65 '     ],$
['arch017', 'IRAS 12032+1707 '        ,'4.11'       ,'16.27'     ,'28.7' , '64920',      '15',     '16'  ,  '65291',    '48 '     ],$
['arch018', 'IRAS 12112+0305 '        ,'2.96'       ,'45.'       ,'23.8' ,      '',      '',       ''    ,  '21885',    '70 '     ],$
['arch020', 'IRAS 12540+5708 '        ,'2.87'       ,'50.'       ,'309.9', '12650',      '50',     ''    ,  '12642',    '4  '     ],$
['arch023', 'IRAS 13218+0552 '        ,'3.45'       ,'4.'        ,'5.3'  , '61268',      '15',     ''    ,  '61488',    '58 '     ],$
['arch024', 'IRAS 13428+5608 '        ,'2.55'       ,'70.'       ,'145.4', '10838',      '50',     ''    ,  '11326',    '9  '     ],$
['arch025', 'IRAS 13451+1232 '        ,'2.38'       ,'1.7'       ,'5398.', '36650',      '18',     ''    ,  '36575',    '70 '     ],$
['arch026', 'IRAS 14070+0525 '        ,'4.13'       ,'10.'       ,'5.2'  , '79760',      '20',     ''    ,  '79259',    '9  '     ],$
['arch029', 'IRAS 15327+2340 '        ,'2.59'       ,'280.'      ,'326.8', '5371',       '15',     ''    ,  '5426 ',    '7  '     ],$
['arch030', 'IRAS 16300+1558 '        ,'2.81'       ,'3.12'      ,'7.9'  , '72528',      '15',     '38'  ,  '72457',    '116'     ],$
;['arch031', 'IRAS 16399-0937 '        ,'1.69'       ,'25.'       ,'57.9' ,      '',      '',       ''    ,  '8098 ',    '28 '     ],$
['arch032', 'IRAS 17207-0014 '        ,'3.04'       ,'131.'      ,'82.4' , '12730',      '90',     ''    ,  '12834',    '9  '     ],$
['arch033', 'IRAS 20100-4156 '        ,'4.05'       ,'200.'      ,'0.'   , '38700',      '50',     ''    ,  '38848',    '51 '     ],$
;['arch034', 'IRAS 20550+1656 '        ,'2.13'       ,'40.'       ,'43.9' , '10900',      '50',     ''    ,  '10822',    '10 '     ],$
['arch035', 'IRAS 21272+2514 '        ,'3.63'       ,'16.33'     ,'4.4'  , '45032',      '14',     '12'  ,  '45275',    '39 '     ],$
['arch036', 'IRAS 22491-1808 '        ,'2.39'       ,'11.'       ,'5.9'  , '23116',      '50',     ''    ,  '23312',    '22 '     ],$
['arch039', 'IRAS 23233+0946 '        ,'2.72'       ,'3.32'      ,'11.6' , '38240',      '14',     '18'  ,  '38384',    '56 '     ],$
['arch040', 'IRAS 23365+3604 '        ,'2.45'       ,'0.'        ,'28.7' ,      '',      '',       ''    ,  '19331',    '9  '     ],$
['arch045', 'IRAS 01355-1814 '        ,'2.69'       ,'11.7'      ,'0.'   ,      '',      '',       ''    ,  '57410',    '57 '     ],$
['arch048', 'IRAS 16090-0139 '        ,'3.46'       ,'9.29'      ,'20.9' ,      '',      '',       ''    ,  '40029',    '40 '     ]]

; log(L_OH) [L_sun]

oharchnames = archarr[0,*]
oharchobj = transpose(archarr[1,*])
logoh_arch = float(archarr[2,*])

	; Convert H_0 = 75 to H_0 = 70
	
	for i = 0, n_elements(oharchnames) - 1 do begin
	
		match, oharchobj[i], dl_obj, a1, a2
		oldoh = logoh_arch[i]
		newoh = oldoh + 2d * alog10(double(d70[a2]) / double(d75[a2]))
;		print,oharchobj[i],oldoh, newoh
		logoh_arch[i] = newoh
	
	endfor

; Hyperfine ratios

arch_rh = strarr(n_elements(oharchnames))
for i = 0, n_elements(oharchnames) - 1 do begin

	targets, oharchnames[i], z, obj, dl
	a = ohmcat(strtrim(obj,2))
	rh_temp = a.rh
	rhm_temp = a.rh_marker
	if rhm_temp eq '' then arch_rh[i] = rh_temp else arch_rh[i] = 0.
	
endfor

; Peak fluxes of 1667 MHz line and 1.4 GHz continuum [mJy]

f1667_arch = float(archarr[3,*])
f1420_arch = float(archarr[4,*])

; Velocity and error for OHM emission

czoh_arch = float(archarr[5,*])
czoh_err_arch = float(archarr[6,*])

; Optical velocities

zerr_opt_arch = 1d-5 * float(archarr[7,*])

cz_opt_arch = float(archarr[8,*])
cz_opt_err_arch = float(archarr[9,*])

;####################
; Control sample limits
;####################

; f1420 for zero flux is a limit of 5 mJy (no detection in NVSS)

;					log OH limit	f1667_limit	f1420       cz_opt	cz_opt_err	opt vel ref
;
conarr = [$
['control004', 'IRAS 08559+1053'        , '1.69'      ,  '0.35'   ,    '0.  '   ,   '44369'   ,  '70 '  ],$
['control008', 'IRAS 13349+2438'        , '1.69'      ,  '0.66'   ,    '20.0'   ,   '32270'   ,  '81 '  ],$
['control013', 'IRAS 15206+3342'        , '1.72'      ,  '0.54'   ,    '11.2'   ,   '37297'   ,  '63 '  ],$
['control023', 'IRAS 00163-1039'        , '1.22'      ,  '4.2 '   ,    '0.  '   ,   '8140 '   ,  '8  '  ],$	; NED
['control024', 'IRAS 06538+4628'        , '0.86'      ,  '2.6 '   ,    '64.3'   ,   '6401 '   ,  '16 '  ],$	; Huchra et al 1983
['control025', 'IRAS 09437+0317'        , '0.98'      ,  '3.7 '   ,    '0.  '   ,   '6136 '   ,  '95 '  ],$	; de Vaucouleurs et al 1976
['control026', 'IRAS 05083+7936'        , '1.91'      ,  '4.3 '   ,    '41.4'   ,   '16090'   ,  '9  '  ],$	; Downes et al 1993
['control028', 'IRAS 10565+2448'        , '1.63'      ,  '3.9 '   ,    '57.0'   ,   '12921'   ,  '9  '  ],$	; Downes et al 1993
['control033', 'IRAS 20460+1925'        , '2.12'      ,  '0.64'   ,    '18.9'   ,   '54262'   ,  '300'  ],$
['control034', 'IRAS 15001+1433'        , '2.01'      ,  '0.61'   ,    '16.9'   ,   '48790'   ,  '70 '  ],$
['control035', 'IRAS 11119+3257'        , '2.04'      ,  '0.49'   ,    '110.4'  ,   '56661'   ,  '70 '  ],$
['control036', 'IRAS 23498+2423'        , '2.22'      ,  '0.59'   ,    '6.8 '   ,   '63556'   ,  '70 '  ],$
['control037', 'IRAS 01572+0009'        , '2.09'      ,  '0.73'   ,    '26.7'   ,   '48869'   ,  '55 '  ],$
['control039', 'IRAS 23007+0836'        , '0.59'      ,  '2.4 '   ,    '181.0'  ,   '4891 '   ,  '2  '  ],$     ; Keel 1996
['control040', 'IRAS 23394-0353'        , '1.15'      ,  '4.2 '   ,    '0.  '   ,   '6966 '   ,  '19 '  ]]	; Huchra et al 1993

ncon = n_elements(condat('tag'))
contags = condat('tag')		
ohconnames = contags
obj = condat('obj')
conobj = strmid(obj,1,strlen(obj[0]))

logoh_limit_con = float(conarr[2,*])
f1667_limit_con = float(conarr[3,*])
f1420_con = float(conarr[4,*])

cz_opt_con = float(conarr[5,*])
cz_opt_err_con = float(conarr[6,*])

for i = 0, n_elements(contags) - 1 do begin

	match, dl_obj, conobj[i], a, b
	oldoh = logoh_limit_con[i]
	newoh = oldoh + 2d * alog10(double(d70[a]) / double(d75[a]))
	logoh_limit_con[i] = newoh

endfor

; Check if numbers agree for all data arrays

if nohm ne n_elements(logoh) then print,'Wrong number of data points in L_OH array for OHMs'
if nohm ne n_elements(f1667_ohm) then print,'Wrong number of data points in 1667 MHz array for OHMs'
if nohm ne n_elements(f1420_ohm) then print,'Wrong number of data points in 1420 MHz array for OHMs'
if nohm ne n_elements(czoh_ohm) then print,'Wrong number of data points in cz_OH array for OHMs'
if nohm ne n_elements(czoh_err_ohm) then print,'Wrong number of data points in cz_OH_err array for OHMs'
if nohm ne n_elements(zerr_opt_ohm) then print,'Wrong number of data points in zerr_opt array for OHMs'

if narch ne n_elements(logoh_arch) then print,'Wrong number of data points in L_OH array for archived OHMs'
if narch ne n_elements(f1667_arch) then print,'Wrong number of data points in 1667 MHz array for archived OHMs'
if narch ne n_elements(f1420_arch) then print,'Wrong number of data points in 1420 MHz array for archived OHMs'
if narch ne n_elements(czoh_arch) then print,'Wrong number of data points in cz_OH array for archived OHMs'
if narch ne n_elements(czoh_err_arch) then print,'Wrong number of data points in cz_OH_err array for archived OHMs'
if narch ne n_elements(zerr_opt_arch) then print,'Wrong number of data points in zerr_opt array for archived OHMs'

if ncon ne n_elements(logoh_limit_con) then print,'Wrong number of data points in limits on L_OH array for control sample'
if ncon ne n_elements(f1667_limit_con) then print,'Wrong number of data points in 1667 MHz limits array for control sample'
if ncon ne n_elements(f1420_con) then print,'Wrong number of data points in 1420 MHz array for control sample'


; Save data to IDL files

; OHMs (Darling)

save,logoh,f1667_ohm,f1420_ohm,czoh_ohm,czoh_err_ohm,zerr_opt_ohm,meganames,mega_rh, cz_opt_ohm, cz_opt_err_ohm, $
	filename='~/Astronomy/Research/Spitzer/OHM/ohlum.sav'

; Archived OHMs

save,logoh_arch,f1667_arch,f1420_arch,oharchnames,czoh_arch,czoh_err_arch,zerr_opt_arch,arch_rh,cz_opt_arch, cz_opt_err_arch, $
	filename='~/Astronomy/Research/Spitzer/archived/data/idl_sav/ohlum_arch.sav'

; Control sample

save,logoh_limit_con,f1667_limit_con,f1420_con,ohconnames,cz_opt_con, cz_opt_err_con,$
	filename='~/Astronomy/Research/Spitzer/control/data/idl_sav/ohlum_limit_con.sav'

;stop
end
