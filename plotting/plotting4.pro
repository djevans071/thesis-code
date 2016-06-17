;plotting
;
home = '../../2014.11.04 - PMMA2.3 in CHB/'
folder1 = 'trial 02'
forces = '/forces.gdf'
diffusivity = '/diffusivity.gdf'
diff0 = '/diff0.gdf'

a = 1.15
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

;finding the potential energy
  uu = (r1[1]-r1[0])*reverse(total(reverse(f1),/cum))
  uferr = (r1[1]-r1[0])*reverse(total(reverse(f1 + fe1), /cum))
  u_err = uferr - uu
  
;fitting
y1 = force_fit(rf1, a, eps=eps)   ;screened-Coulomb force fit
y2 = force_fit2(rf1, a, eps=eps)   ;unscreened-Coulomb force fit
zp1 = pot_fit(rf1, a, eps=eps)     ;screened-Coulomb potential fit
zp2 = pot_fit2(rf1, a, eps=eps)    ;unscreened-Coulomb potential fit

marg = {xtitle:'r [ $\mu$m]', font_size:16, xminor:1, yminor:2, $
  xrange:[4,20], font_name:'Times'}
  
;-----------------COLOR TABLES----------------------------
;diffusion color table
loadct, 59, rgb_table=table_d
;force color table
loadct, 55, rgb_table=table_f

colors = bytscl(dist)
  
;-----------------PROB DIST PLOT--------------------------

distp1 = plot(r1, dist, 'k3', _extra=marg, yrange=[0.01,.13],$
  position = [0.20,0.863,0.95,0.98], ytitle='$\rho(r)$',ytickinterval=0.04)
  (distp1.axes)[0].showtext=0
(distp1.axes)[0].minor=1 & (distp1.axes)[2].minor=1
;(distp1.axes)[1].major=4 & (distp1.axes)[3].major=4
;-----------------VELOCITY PLOT---------------------------
v1 = v1+0.2

vp1 = plot(r1, v1, 'k1',ytitle='$v(r)$ [ $\mum$/s]', $
  _extra=marg,ylog=1, yrange=[0.2,20], /current, $
  position=[0.20,0.613,0.95,0.86])
;position=[0.15,0.54,0.95,0.98])
(vp1.axes)[0].showtext=0
(vp1.axes)[0].minor=1 & (vp1.axes)[2].minor=1
(vp1.axes)[1].minor=9 & (vp1.axes)[3].minor=9

vp1_err = polygon([r1, reverse(r1)],[v1+ve1,reverse(v1-ve1)], $
  target=vp1, /data, /fill_background, $
  fill_color='gold', linestyle=3)
vp1 = plot(r1,v1, 'k2', overplot=1, font_name='Times')

;-----------------DIFFUSION PLOT---------------------------

dp1 = plot(r1, d1, 'k1', ytitle='$D(r)$ [$\mum^2$/s]', $
  _extra = marg, yrange=[.09,0.23], /current,$
  position=[0.20,0.353,0.95,0.61], $
  ytickinterval=0.04)

(dp1.axes)[0].showtext=0
(dp1.axes)[1].minor=3 & (dp1.axes)[3].minor=3

;error estimate with densities as a color gradient  
 for j=0,n_elements(r1)-2 do begin $  
  xpoly = [ r1[j], r1[j], r1[j+1], r1[j+1] ]
  ypoly = [ (d1-de1)[j], (d1+de1)[j], (d1+de1)[j+1], (d1-de1)[j+1] ]

  poly_test = polygon(xpoly, ypoly, target=dp1, /data, $
    linestyle=6, fill_background=1, $
    fill_color=reform(table_d(colors[j],*)) )
endfor
dp1_err = polygon([r1, reverse(r1)],[d1+de1,reverse(d1-de1)], $
  target=dp1, /data, linestyle=3, fill_transparency=100)
  
dp1 = plot(r1,d1, 'k-2', overplot=1, font_name='Times')
dp2 = plot(r1,2.*d0, 'b--2', overplot=1, font_name='Times')
dpf1 = plot(r1,z1, 'r-2', overplot=1, font_name='Times')
;dpf2 = plot(r2,z2, 'b-3', overplot=1, /current)

;------------------FORCE PLOT------------------------------
f1 = f1*1000 
fe1 = fe1*1000 
y1 =  y1*1000 

p1 = plot(r1, f1, 'k1', ytitle='$F(r)$ [fN]',yrange=[1,900], $
  _extra=marg, ylog=1, /current, $
  position=[0.20,0.08,0.95,0.35])
  
(p1.axes)[0].minor=1 & (p1.axes)[2].minor=1
(p1.axes)[1].minor=9 & (p1.axes)[3].minor=9

;error estimate with densities as a color gradient
for j=0,n_elements(r1)-2 do begin $
    xpoly = [ r1[j], r1[j], r1[j+1], r1[j+1] ]
    ypoly = [(f1-fe1)[j],(f1+fe1)[j],(f1+fe1)[j+1],(f1-fe1)[j+1] ]

  poly_test = polygon(xpoly, ypoly, target=p1, /data, $
    linestyle=6, fill_background=1, $
    fill_color=reform(table_f(colors[j],*)) )
endfor  
p1_err = polygon([r1, reverse(r1)],[(f1+fe1),reverse((f1-fe1))], $
  target=p1, /data, linestyle=3, fill_transparency=100)

p1 = plot(r1,f1, 'k2', overplot=1, font_name='Times')

pf1 = plot(r1,y1, 'r-2', /current, overplot=1, font_name='Times')

end