;+
; NAME:
;       
;	ERRORFIT_LFL
;
; PURPOSE:
;
;	Function to minimize for fitting IRS spectra to CLUMPY models with data in lambda * F_lambda
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
;       Written by K. Willett            Nov 09    
;-

function errorfit_lfl, p

common func_lfl, fmodel_lfl, fdata_lfl, ferror_lfl, function_lambda 

n = n_elements(fdata_lfl)

return, max(1./n * sqrt( total( ( (p[0] * fmodel_lfl - function_lambda * fdata_lfl)/ferror_lfl)^2) ) )

end
