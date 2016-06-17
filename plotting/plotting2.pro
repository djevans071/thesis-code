;plotting2
;
home = '../../2014.11.04 - PMMA2.3 in CHB/'  
diffusivity = file_search(home + '*/diffusivity.gdf', count=n)  
forces = file_search(home + '*/forces.gdf', count=n)  
  
a = 1.13
astring = strtrim(string(a,format='(f5.2)'),2)  
eps = 5.989  
  
rd1 = read_gdf(diffusivity[0])  
rf1 = read_gdf(forces[0])  
  
r1 = reform(rd1[0,*])
rf1 = reform(rf1[0,*])  
v1 = reform(rd1[3,*])  
d1 = reform(rd1[1,*])  
f1 = reform(rf1[1,*])  
  
marg = {xtitle:'r [ $\mu$m]', font_size:16, xminor:1, yminor:2, $
  xrange:[4,21]}  
  
;--------VELOCITY PLOT-------------------------------------------
;initial plot
vp1 = plot(r1, v1, 'k3',ytitle='V [ $\mum$/s]', $
  title='a='+astring + ' $\mu m$', $
  _extra=marg, yrange=[-1,15.5], ytickinterval=4, $
  position=[0.15,0.65,0.95,0.94])  
    (vp1.axes)[0].hide=1  
    ;(vp1.axes)[0].minor=1 & (vp1.axes)[2].minor=1  
    (vp1.axes)[1].minor=1 & (vp1.axes)[3].minor=1  
    
  zero = plot(findgen(30), fltarr(30), 'k--2', overplot=1)  
    
for i=0,n-1 do begin $ ;plot velocity errors
  rdi = read_gdf(diffusivity[i])  
  
  ri = reform(rdi[0,*])  
  vi = reform(rdi[3,*])  
  vei = reform(rdi[4,*])  
  
  vpi_err = polygon([ri, reverse(ri)],[vi+vei,reverse(vi-vei)], $
    /data, /fill_background, fill_color=[255,215 - 30*i,0], linestyle=2, $
    fill_transparency=60, target=vp1)  
endfor  

for i=0,n-1 do begin $  ;plot velocities
  rdi = read_gdf(diffusivity[i])  
  
  ri = reform(rdi[0,*])  
  vi = reform(rdi[3,*])  
  vei = reform(rdi[4,*])  

  vpi = plot(ri,vi, thick=3, color=[0,30*i,0+70*i], overplot=1)  
endfor  


;--------DIFFUSIVITY PLOT-------------------------------------------
  ;initial plot
  dp1 = plot(r1, d1, 'k1', ytitle='D [$\mum^2$/s]', $
    _extra = marg, xminor=1, yrange=[.215,.325], $
    position=[0.15,0.36,0.95,0.65], ytickinterval=0.02, /current)  
      (dp1.axes)[0].hide=1  
      ;(dp1.axes)[0].minor=1 & (dp1.axes)[2].minor=1  
      (dp1.axes)[1].minor=1 & (dp1.axes)[3].minor=1  
  
for i=0,n-1 do begin $
  rdi = read_gdf(diffusivity[i])  
  
  ri = reform(rdi[0,*])  
  di = reform(rdi[1,*])  
  dei = reform(rdi[2,*])  

  ;error plot
  dpi_err = polygon([ri, reverse(ri)],[di+dei,reverse(di-dei)], $
    /data, /fill_background, fill_color=[219,112+20*i,147+30*i], $
    fill_transparency=70, linestyle=2, target=dp1)  
endfor  

for i=0,n-1 do begin $
  rdi = read_gdf(diffusivity[i])  

  ri = reform(rdi[0,*])  
  di = reform(rdi[1,*])  
  dei = reform(rdi[2,*])  

  ;diffusivity plot
  dpi = plot(ri,di, 'o1', color=[0,30*i,0+70*i], overplot=1)  
  ;diffusivity fit plot
  zi = diff_fit(rdi, a, motion='rm_p')  
  dpfi = plot(ri,zi, 'r-3', overplot=1)  
endfor  


;--------FORCE PLOT-------------------------------------------
;initial plot
p1 = plot(rf1, f1, 'k1', ytitle='F [pN]', $
  yrange=[-0.05,0.29], ytickinterval=0.2, $
  _extra=marg, xminor=1, $
  position=[0.15,0.07,0.95,0.36], /current)  
  
    (p1.axes)[0].minor=1 & (p1.axes)[2].minor=1  
    (p1.axes)[1].minor=1 & (p1.axes)[3].minor=1  

zero = plot(findgen(30), fltarr(30), 'k--2', overplot=1)  

for i=0,n-1 do begin $ ;plot force errors
  rfi = read_gdf(forces[i])  

  ri = reform(rfi[0,*])  
  fi = reform(rfi[1,*])  
  fei = reform(rfi[2,*])  

  pi_err = polygon([ri, reverse(ri)],[fi+fei,reverse(fi-fei)], $
    /data, /fill_background, fill_color=[255,127+20*i,80], linestyle=2, $
    fill_transparency=50, target=p1)  
endfor  

for i=0,n-1 do begin $
  rfi = read_gdf(forces[i])  

  ri = reform(rfi[0,*])  
  fi = reform(rfi[1,*])  
  fei = reform(rfi[2,*])  

  ;force plot
  pi = plot(ri,fi, thick=3, color=[0,30*i,0+70*i], overplot=1)  
  ;force fit plot
  yi = force_fit(rfi, a, eps=eps)  
  pfi = plot(ri,yi, 'r-2', overplot=1)
endfor  

end  