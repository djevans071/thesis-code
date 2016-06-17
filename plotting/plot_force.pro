;@plot_vel
;
home = '../../2016.05.09 - 3 PMMA3 in CXB/' & $
  forces = file_search(home + '*/forces.gdf', count=n) & $
  
  a = 1.5 & $
  astring = strtrim(string(a,format='(f5.2)'),2) & $
  eps = 5.989 & $
  
  rf1 = read_gdf(forces[0]) & $
  
  r1 = reform(rf1[0,*]) & $
  f1 = reform(rf1[1,*]) & $
  
  marg = {xtitle:'r [ $\mu$m]', font_size:16, xminor:1, yminor:2, $
  xrange:[4,22]} & $
  
  ;initial plot  
  p1 = plot(r1, f1, 'k1', ytitle='F [pN]',yrange=[-0.05,0.8],$
    title='Forces for a='+astring + ' $\mu m$', $
    aspect_ratio=12, $
    _extra=marg,xminor=1) & $
      (p1.axes)[0].minor=1 & (p1.axes)[2].minor=1 & $
      (p1.axes)[1].minor=1 & (p1.axes)[3].minor=1 & $

  zero = plot(findgen(30), fltarr(30), 'k--2', overplot=1) & $      
  
for i=0,n-1 do begin $
  
  rfi = read_gdf(forces[i]) & $
  
  ri = reform(rfi[0,*]) & $
  fi = reform(rfi[1,*]) & $
  fei = reform(rfi[2,*]) & $
  
  pi_err = polygon([ri, reverse(ri)],[fi+fei,reverse(fi-fei)], $
    /data, /fill_background, fill_color=[255,127+20*i,80], linestyle=6, $
    fill_transparency=50) & $
    
endfor & $

for i=0,n-1 do begin $

  rfi = read_gdf(forces[i]) & $

  ri = reform(rfi[0,*]) & $
  fi = reform(rfi[1,*]) & $
  fei = reform(rfi[2,*]) & $

  ;force fitting
  yi = force_fit(rfi, a, eps=eps) & $

  ;force plot
  pi = plot(ri,fi, thick=3, color=[0,30*i,0+70*i], overplot=1) & $
  ;force fit plot
  pfi = plot(ri,yi, 'r-2', overplot=1) & $

endfor & $

home2 = '../../2015.03.18 - PMMA3 in CHB/trial 01'
forces2 = file_search(home2 + '*/forces.gdf', count=n) & $

  rf2 = read_gdf(forces2[0]) & $

  r2 = reform(rf2[0,*]) & $
  f2 = reform(rf2[1,*]) & $
  fe2 = reform(rf2[2,*]) & $
  
  p2 = plot(r2, f2, thick=3, overplot=1)
  p2_err = polygon([r2, reverse(r2)],[f2+fe2, reverse(f2-fe2)], $
    /data, /fill_background, fill_color=[255,127+20*i,80], linestyle=6, $
    fill_transparency=50) &

end