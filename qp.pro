pro qp, fname

; Very quick routine to plot the individual, fully calibrated models (so I can see what pixels to remove in observed frame)

!p.multi=[0,1,1]
tag,fname,dirtag
targets,fname,r,obj
file = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'+fname+'_lh_cal.tbl'

readcol, file, o, w, f, skipline = 1, /silent
noneg = where(f gt 0)
plot,w[noneg],f[noneg], psym = 10, title=fname+' - '+obj;, /ylog

end
