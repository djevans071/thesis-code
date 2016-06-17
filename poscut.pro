function poscut, RF

R_ = reform(RF[0,*])

rmax = max(R_)                         ;relative displacement binning
rmin = min(R_)
rwidth = (rmax - rmin)
bins = 500
x = rwidth*(findgen(bins)/(bins-1.)) + rmin

;h = histogram(R_, nbins=bins, min=rmin, max=rmax)
    ;plot, x,h
;w = where(h gt 0.)
;x = x[w]

;k = kde(R_,x)
;oplot, x, k, line=3
;stop
return, x
end