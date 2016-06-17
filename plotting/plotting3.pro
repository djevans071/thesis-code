;plotting3
;
home = '../../2014.11.04 - PMMA2.3 in CHB/'
folder1 = 'trial 02'
forces = '/forces.gdf'
charges = '/charges.gdf'
diffusivity = '/diffusivity.gdf'
diff0 = '/diff0.gdf'
charges_p = '/p_charges.gdf'

a = 1.13
eps = 5.989
T = 296.
kb = 1.3806488e-5       ;Boltzmann's constant (pN um/K)
kbt = kb * T            ;Boltzmann energy (pN um)

rf1 = read_gdf(home + folder1 + forces)
rd1 = read_gdf(home + folder1 + diffusivity)
dd0 = read_gdf(home + folder1 + diff0)
;rd1 = rd1[*,20:-20]

;rf1 = rf1[*,0:-40]
;rd1 = rd1[*,0:-40]

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
z1 = pot_fit(rf1, a, eps=eps)     ;screened-Coulomb potential fit
z2 = pot_fit2(rf1, a, eps=eps)    ;unscreened-Coulomb potential fit


marg = {xtitle:'r [ $\mu$m]', font_size:16, xminor:1, yminor:2, $
  xrange:[4,21], font_name:'Times'}
  
;-----------------PROB DIST PLOT--------------------------

;distp1 = plot(r1, dist, 'k3', _extra=marg, yrange=[0.01,.09],$
;  position = [0.20,0.863,0.95,0.98], ytitle='$\rho(r)$',ytickinterval=0.04)
;  (distp1.axes)[0].showtext=0
;(distp1.axes)[0].minor=1 & (distp1.axes)[2].minor=1
;(distp1.axes)[1].major=4 & (distp1.axes)[3].major=4

;-----------------FORCE PLOT---------------------------

p1 = plot(r1, 1000*f1, 'k1', ytitle='$F(r)$ [fN]',yrange=[-5,600], $
  _extra=marg, position=[0.20,0.533,0.95,0.98])
 (p1.axes)[0].showtext=0
(p1.axes)[0].minor=1 & (p1.axes)[2].minor=1
(p1.axes)[1].minor=5 & (p1.axes)[3].minor=5

p1_err = polygon([r1, reverse(r1)],[1000*(f1+fe1),reverse(1000.*(f1-fe1))], $
  target=p1, fill_transparency=20, $
  /data, /fill_background, fill_color='coral', linestyle=3)
p1 = plot(r1,1000.*f1, 'k3', overplot=1, font_name='Times')

pf1 = plot(r1,1000.*y1, 'r-3', /current, overplot=1, $
  name='Screened', font_name='Times' )
pf2 = plot(r1,1000.*y2, 'b-3', /current, overplot=1, $
  name='Unscreened', font_name='Times')

leg = legend(target=[pf1,pf2], position=[20,550], $
  /data,font_name='Times', shadow=0)

;-----------------Screened-POTENTIAL PLOT---------------------------

pot_p1 = plot(r1, uu/kbt, 'k1', ytitle='$U(r)/k_BT$',yrange=[-50,450], $
  _extra=marg, /current, $
  position=[0.20,0.08,0.95,0.53])
 ;(pot_p1.axes)[0].showtext=0 
(pot_p1.axes)[0].minor=1 & (pot_p1.axes)[2].minor=1
(pot_p1.axes)[1].minor=5 & (pot_p1.axes)[3].minor=5

pot_p1_err = polygon([r1, reverse(r1)],[(uu+u_err)/kbt, $
  reverse((uu-u_err)/kbt)], $
  target=pot_p1, fill_transparency=20, $
  /data, /fill_background, fill_color='gold', linestyle=3)
pot_p1 = plot(r1, uu/kbt, 'k3', overplot=1, font_name='Times')

pot_pf1 = plot(r1,z1/kbt, 'r-3', /current, overplot=1, font_name='Times')

;------------------Unscreened-POTENTIAL PLOT------------------------------

;pot_p2 = plot(r1, uu/kbt, 'k1', ytitle='$U(r)/k_bT$',yrange=[-50,450], $
;  _extra=marg, /current, $
;  position=[0.20,0.08,0.95,0.53])
  
;(pot_p2.axes)[0].minor=1 & (pot_p2.axes)[2].minor=1
;(pot_p2.axes)[1].minor=5 & (pot_p2.axes)[3].minor=5

;pot_p2_err = polygon([r1, reverse(r1)],[(uu+u_err)/kbt, $
;  reverse((uu-u_err)/kbt)], $
;  target=pot_p2, fill_transparency=20, $
;  /data, /fill_background, fill_color='gold', linestyle=3)
;pot_p2 = plot(r1, uu/kbt, 'k3', overplot=1, font_name='Times')

pot_pf2 = plot(r1,z2/kbt, 'b-3', /current, overplot=1, font_name='Times')

end