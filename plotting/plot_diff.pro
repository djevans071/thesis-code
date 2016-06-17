;plot_diff
;

home = '../../2014.12.10 - PMMA1 in CHB/' & $
diffusivity = file_search(home + '*/diffusivity.gdf', count=n) & $
  
  a = 0.55 & $
    astring = strtrim(string(a,format='(f5.2)'),2) & $
  eps = 5.989 & $
  
  rd1 = read_gdf(diffusivity[0]) & $
  
  r1 = reform(rd1[0,*]) & $
  d1 = reform(rd1[1,*]) & $

  marg = {xtitle:'r [ $\mu$m]', font_size:16, xminor:1, yminor:2, $
    xrange:[4,22]} & $

  ;initial plot
  dp1 = plot(r1, d1*100., 'k1', ytitle='D [$\times 10^{-2}$ $\mum^2$/s]', $
    title='Diffusivities for a='+astring + ' $\mu m$', $
    aspect_ratio=1, $
    _extra = marg, xminor=1, yrange=[7,15.9]) & $
      (dp1.axes)[0].minor=1 & (dp1.axes)[2].minor=1 & $
      (dp1.axes)[1].minor=1 & (dp1.axes)[3].minor=1 & $
  
for i=0,n-1 do begin $
  
  rdi = read_gdf(diffusivity[i]) & $
  
  ri = reform(rdi[0,*]) & $
  di = reform(rdi[1,*]) & $
  dei = reform(rdi[2,*]) & $
  
  ;error plot
  dpi_err = polygon([ri, reverse(ri)],[100.*(di+dei),reverse(100.*(di-dei))], $
    /data, /fill_background, fill_color=[219,112+20*i,147+30*i], $
    fill_transparency=70, linestyle=3) & $
endfor & $


for i=0,n-1 do begin $

  rdi = read_gdf(diffusivity[i]) & $

  ri = reform(rdi[0,*]) & $
  di = reform(rdi[1,*]) & $
  dei = reform(rdi[2,*]) & $

  ;diffusion fitting
  zi = diff_fit(rdi, a, motion='rm_p') & $  

  ;diffusivity plot
  dpi = plot(ri,100.*di, 'o1', color=[0,30*i,0+70*i], overplot=1, /current) & $
  ;diffusivity fit plot
  dpfi = plot(ri,100.*zi, 'r-3', /current, overplot=1) & $
endfor & $

end