restore,'~/Astronomy/Research/Spitzer/OHM/fname.sav'
fnames=[fname,'control004','control008','control013']
nf = n_elements(fnames)
for i = 0, nf - 1 do begin
	fullql_hr,fnames(i),xr=[13.5,16.5],yr=[0,0.5]
	ver,14.365,linestyle=2
	junk = get_kbrd()
endfor

end
