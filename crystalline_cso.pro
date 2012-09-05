
;+
; NAME:
;       
;	CRYSTALLINE_CSO
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

ps = 1

tags = csodat('tag')
tags = tags[0:8]
ncso = n_elements(tags)

	if keyword_set(ps) then begin
		set_plot,'ps'
		device,filename = '~/Astronomy/Research/Spitzer/cso/plots/cryst_sil_cso.ps', /color, /portrait, xs=18, ys=24, xoff=1, yoff=1
		cs = 1
		ls = 2
	endif else begin
		cs = 1
		ls = 1
	endelse
	
;!p.multi = [0,1,ncso]
erase 
multiplot, [1,ncso], $
	mxtitle = 'Wavelength [!7l!3m]', mytitle = 'Residual optical depth', mtitle='Crystalline silicates in CSOs', $
	mxtitsize = 1.2, $
	mytitsize = 1.2, $
	mtitsize = 1.2

for i = 0, ncso - 1 do begin

	crystalline, tags[i]
	multiplot
endfor

multiplot,/reset

	if keyword_set(ps) then begin
		device, /close
		set_plot,'x'
	endif

; Deepest silicate absorbers in the OHM sample

;tags = ['arch005','mega008','arch009','arch010','arch014','mega016','arch017','arch026','arch029','arch030','arch033','arch035','mega032']

tagohm = ohmdat('tag')
tagarch = archdat('tag')

alltag = [transpose(tagohm),transpose(tagarch)]
alltag = alltag[where(alltag ne 'mega034')]

for j = 0, 4 do begin
	
	nobj = 11
	
		if keyword_set(ps) then begin
			set_plot,'ps'
			device,filename = '~/Astronomy/Research/Spitzer/ohm/plots/cryst_sil_ohm'+strtrim(j+1,2)+'.ps', /color, /portrait, xs=18, ys=24, xoff=1, yoff=1
			cs = 1
			ls = 2
		endif else begin
			cs = 1
			ls = 1
		endelse
		
	erase 
	multiplot, [1,nobj], $
		mxtitle = 'Wavelength [!7l!3m]', mytitle = 'Residual optical depth', mtitle='Crystalline silicates in OHMs', $
		mxtitsize = 1.2, $
		mytitsize = 1.2, $
		mtitsize = 1.2
	
	for i = 0, nobj - 1 do begin
		crystalline, alltag[11*j + i]
		multiplot
	endfor
	
	multiplot,/reset
	
		if keyword_set(ps) then begin
			device, /close
			set_plot,'x'
		endif
	
endfor

; Control sample


contag = condat('tag')

	
	ncon = n_elements(contag)
	
		if keyword_set(ps) then begin
			set_plot,'ps'
			device,filename = '~/Astronomy/Research/Spitzer/control/plots/cryst_sil_con.ps', /color, /portrait, xs=18, ys=24, xoff=1, yoff=1
			cs = 1
			ls = 2
		endif else begin
			cs = 1
			ls = 1
		endelse
		
	erase 
	multiplot, [1,ncon], $
		mxtitle = 'Wavelength [!7l!3m]', mytitle = 'Residual optical depth', mtitle='Crystalline silicates in non-masing galaxies', $
		mxtitsize = 1.2, $
		mytitsize = 1.2, $
		mtitsize = 1.2
	
	for i = 0, ncon - 1 do begin
		crystalline, contag[i],/batch
		multiplot
	endfor
	
	multiplot,/reset
	
		if keyword_set(ps) then begin
			device, /close
			set_plot,'x'
		endif
	

end
