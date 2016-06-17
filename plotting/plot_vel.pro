;@plot_vel
;
home = '../../2014.12.10 - PMMA1 in CHB/' & $
  diffusivity = file_search(home + '*/diffusivity.gdf', count=n) & $
    
  a = 0.55 & $
    astring = strtrim(string(a,format='(f5.2)'),2) & $
  eps = 5.989 & $
  
  rd1 = read_gdf(diffusivity[0]) & $

  r1 = reform(rd1[0,*]) & $
  v1 = reform(rd1[3,*]) & $
  
  marg = {xtitle:'r [ $\mu$m]', font_size:16, xminor:1, yminor:2, $
  xrange:[4,22]} & $
  
  ;initial plot
  vp1 = plot(r1, v1, 'k3',ytitle='V [ $\mum$/s]', $
    title='Velocities for a='+astring + ' $\mu m$', $
    aspect_ratio=0.55, $
    _extra=marg, yrange=[-1,16]) & $
      (vp1.axes)[0].minor=1 & (vp1.axes)[2].minor=1 & $
      (vp1.axes)[1].minor=1 & (vp1.axes)[3].minor=1 & $
  zero = plot(findgen(30), fltarr(30), 'k--2', overplot=1) & $
  
for i=0,n-1 do begin $
  
  rdi = read_gdf(diffusivity[i]) & $
  
  ri = reform(rdi[0,*]) & $
  vi = reform(rdi[3,*]) & $ 
  vei = reform(rdi[4,*]) & $
  
  vpi_err = polygon([ri, reverse(ri)],[vi+vei,reverse(vi-vei)], $
    /data, /fill_background, fill_color=[255,215 - 30*i,0], linestyle=3, $
    fill_transparency=60) & $
    
  endfor & $

for i=0,n-1 do begin $
  
  rdi = read_gdf(diffusivity[i]) & $
  
  ri = reform(rdi[0,*]) & $
  vi = reform(rdi[3,*]) & $
  vei = reform(rdi[4,*]) & $
  
  vpi = plot(ri,vi, thick=3, color=[0,30*i,0+70*i], overplot=1, /current) & $
  
endfor & $
  
  end