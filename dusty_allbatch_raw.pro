
;+
; NAME:
;       
;	DUSTY_ALLBATCH_RAW
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
;       Written by K. Willett                
;-

ohmtag = ohmdat('tag')
archtag = archdat('tag')
contag = condat('tag')

alltag = [transpose(ohmtag),transpose(archtag),transpose(contag)]

for i=0, n_elements(alltag) - 1 do dusty_compare_batch, alltag[i]

end
