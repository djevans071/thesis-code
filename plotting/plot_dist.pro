;@plot_dist
;
rf = read_gdf('../../2014.12.10 - PMMA1 in CHB/trial 03/RF1rm.gdf')

R_ = reform(RF[0,*])

rmax = max(R_)                         ;relative displacement binning
rmin = min(R_)
rwidth = (rmax - rmin)
bins = 500
x = rwidth*(findgen(bins)/(bins-1.)) + rmin
    binwidth = rwidth/(bins-1.)


bins2 = 50
x2 = rwidth*(findgen(bins2)/(bins2-1.)) + rmin 
    binwidth2 = rwidth/(bins2-1.)
h = histogram(R_, nbins=bins, min=rmin, max=rmax)
h2 = histogram(R_, nbins=bins2, min=rmin, max=rmax) 
plot, x,h/total(h*binwidth)
oplot, x2, h2/total(h2*binwidth2), line=3

k = kde(R_,x)
k2 = kde(R_,x2)

end