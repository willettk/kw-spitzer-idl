
;+
; NAME:
;       
;	DLERR_CSO
;
; PURPOSE:
;
;	Save the CSO luminosity distance errors from NED to an IDL save file
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
;	IDL> .r dlerr_cso
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Aug 09
;-

; From NED

dlerrarr = [['cso001','0.9'], $
            ['cso002','4.7'], $
            ['cso003','1.4'], $
            ['cso004','6.1'], $
            ['cso005','19.5'], $
            ['cso006','9.7'], $
            ['cso007','7.9'], $
            ['cso008','6.0'], $
            ['cso009','1.1'], $
            ['cso010','8.5']]

save, filename='~/Astronomy/Research/Spitzer/cso/data/idl_sav/dlerr.sav',dlerrarr

end
