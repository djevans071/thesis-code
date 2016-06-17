;plotting       creates rho_vdf_plot
;
home = '../../2014.11.04 - PMMA2.3 in CHB/' 
  folder1 = 'trial 01' 
  forces = '/forces.gdf' 
  diffusivity = '/diffusivity.gdf' 
  diff0 = '/diff0.gdf'

a = 1.13
eps = 5.989 

rf1 = read_gdf(home + folder1 + forces) 
rd1 = read_gdf(home + folder1 + diffusivity) 
dd0 = read_gdf(home + folder1 + diff0) 
  ;rd1 = rd1[*,20:-20] 

rf1 = rf1[*,0:-40] 
rd1 = rd1[*,0:-40]   

r1 = reform(rf1[0,*]) & 
f1 = reform(rf1[1,*]) & fe1 = reform(rf1[2,*])  
d1 = reform(rd1[1,*]) & de1 = reform(rd1[2,*]) 
v1 = reform(rd1[3,*]) & ve1 = reform(rd1[4,*]) 

dist = reform(rf1[3,*])
  
;----------------DIFFUSION FITTING---------------------------

z1 = diff_fit(rd1, a, motion='rm_p') 
d0 = fltarr(n_elements(r1)-1) + dd0[0] 

;---------------FORCE FITTING------------------------------

y1 = force_fit(rf1, a, eps=eps) 

marg = {xtitle:'r [ $\mu$m]', font_size:16, xminor:1, yminor:2, $
        xrange:[4,21], font_name:'Times'} 

;-----------------PROB DIST PLOT--------------------------

distp1 = plot(r1, dist, 'k3', _extra=marg, yrange=[0.01,.13],$
  position = [0.20,0.863,0.95,0.98], ytitle='$\rho(r)$',ytickinterval=0.04)
  (distp1.axes)[0].showtext=0
  (distp1.axes)[0].minor=1 & (distp1.axes)[2].minor=1
  ;(distp1.axes)[1].major=4 & (distp1.axes)[3].major=4 
;-----------------VELOCITY PLOT---------------------------

vp1 = plot(r1, v1, 'k1',ytitle='$v(r)$ [ $\mum$/s]', $
  _extra=marg,ylog=1, yrange=[0.2,30], /current, $
  position=[0.20,0.613,0.95,0.86]) 
  ;position=[0.15,0.54,0.95,0.98])
    (vp1.axes)[0].showtext=0
    (vp1.axes)[0].minor=1 & (vp1.axes)[2].minor=1 
    (vp1.axes)[1].minor=9 & (vp1.axes)[3].minor=9 
  
vp1_err = polygon([r1, reverse(r1)],[v1+ve1,reverse(v1-ve1)], $
    target=vp1, $
    /data, /fill_background, fill_color='gold', linestyle=3) 
vp1 = plot(r1,v1, 'k3', overplot=1, font_name='Times') 
  
;-----------------DIFFUSION PLOT---------------------------

dp1 = plot(r1, d1, 'k1', ytitle='$D(r)$ [$\mum^2$/s]', $
  _extra = marg, yrange=[.09,0.23], /current,$
  position=[0.20,0.353,0.95,0.61], $
  ;position=[0.15,0.10,0.95,0.54], $
  ytickinterval=0.04) 
    (dp1.axes)[0].showtext=0
    ;(dp1.axes)[0].minor=1 & (dp1.axes)[2].minor=1 
    (dp1.axes)[1].minor=3 & (dp1.axes)[3].minor=3 
  
dp1_err = polygon([r1, reverse(r1)],[d1+de1,reverse(d1-de1)], $
  target=dp1, fill_transparency=60, $
  /data, /fill_background, fill_color='medium violet red', linestyle=3) 
dp1 = plot(r1,d1, 'ko1', overplot=1, font_name='Times') 

dp2 = plot(r1,2.*d0, 'b--3', overplot=1, font_name='Times') 


dpf1 = plot(r1,z1, 'r-3', overplot=1, font_name='Times') 
;dpf2 = plot(r2,z2, 'b-3', overplot=1, /current) 

;------------------FORCE PLOT------------------------------

p1 = plot(r1, 1000*f1, 'k1', ytitle='$F(r)$ [fN]',yrange=[3,900], $
  _extra=marg, ylog=1, /current, $
  position=[0.20,0.08,0.95,0.35]) 
  
    (p1.axes)[0].minor=1 & (p1.axes)[2].minor=1 
    (p1.axes)[1].minor=9 & (p1.axes)[3].minor=9 

p1_err = polygon([r1, reverse(r1)],[1000*(f1+fe1),reverse(1000.*(f1-fe1))], $
  target=p1, fill_transparency=20, $
  /data, /fill_background, fill_color='coral', linestyle=3) 
p1 = plot(r1,1000.*f1, 'k3', overplot=1, font_name='Times') 

pf1 = plot(r1,1000.*y1, 'r-3', /current, overplot=1, font_name='Times') 
  
  end