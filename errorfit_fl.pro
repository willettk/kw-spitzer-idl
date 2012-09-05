;+
; NAME:
;       
;	ERRORFIT_FL
;
; PURPOSE:
;
;	Function to minimize for fitting IRS spectra to CLUMPY models with data in F_lambda
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
;	Used with CLUMPY.pro
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
;-

function errorfit_fl, p

common func_fl, fmodel_fl, fdata_fl, ferror_fl

n = n_elements(fdata_fl)

return, max(1./n * sqrt( total( ( (p[0] * fmodel_fl - fdata_fl)/ferror_fl)^2) ) )

end
