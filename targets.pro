pro targets, fname, redshift, obj, dl
;+
; NAME:
;       TARGETS
;
; PURPOSE:
; 	Lookup table for the targets in the Spitzer archive
;
; INPUTS:
;
;	FNAME - string detailing filename for target (eg, 'mega001' or 'cso001')
;
; OUTPUTS:
;
;	REDSHIFT - redshift of target
;
;	OBJ - string with NED-readable name of target
;
;	DL - luminosity distance of object (Darling 2002a or NED; difference of h_75 for Darling, h_73 for NED)
;
; KEYWORDS:
;
; EXAMPLE:
;	IDL> targets, 'mega005', redshift, obj, dl
;
; REQUIRES:
;
; REVISION HISTORY
;       Written by K. Willett                May 2007
;	Added ULIRGs - KW, Jun 07
;	Added luminosity distances - KW, Sep 07
; 	Changed this to an actual lookup table, with columns and such - Feb 08
;-

file = '~/Astronomy/Research/Spitzer/ohm_spitz.dat'

readcol, file, format='a,a,a,a,a,a', comment=';', delim='&', $
	taglist, obj_list, redshift_list, dl75_list, dl70_list, /silent

for i = 0, n_elements(taglist)-1 do taglist[i] = strmid(taglist[i],0,strlen(taglist[i])-1)
redshift_list = double(redshift_list)
dl75_list = float(dl75_list)
dl70_list = float(dl70_list)

; Find the appropriate object

match = strmatch(taglist, fname)
match_ind = where(match eq 1, mcount)

if mcount ne 1 then begin
	on_error, 2
	message, 'Did not find matching source in Spitzer list'
endif

redshift = redshift_list[match_ind]
obj = obj_list[match_ind]

; Takes dl70 where available

dl = dl70_list[match_ind] > dl75_list[match_ind]

; Return values for redshift, object designation, and luminosity distance

redshift = redshift[0]
obj = obj[0]
dl = dl[0]

end
