;@plot_charge            
;
function raman_1, x, surf

  Zlamb_a = -2.*alog(x) + 2.*alog(-alog(x)) + 4.*alog(2.) - $
    0.5*alog( (surf - 2.*alog(x))/(surf + 2.*alog(x))  )

return, Zlamb_a
end


function raman, x, surf

  Zlamb_a = fltarr(n_elements(x)) + surf
  w = where(raman_1(x,surf) gt surf)
  Zlamb_a[w] = raman_1(x[w],surf)

  return, Zlamb_a
end


s = file_search('../../*/*/charge.gdf', count=n) 
ss = file_search('../../*/*/charge2.gdf', count=n) 


rz = fltarr(8,1) 

;rz = [a, eps, pot,pot_err, Z, Z_err, kappa, kappa_err]
;      0  1    2   3        4  5      6      7
for i=0,n-1 do begin $
  rzi = reform(read_gdf(s[i])) 
  rz = [[rz],[rzi]] 
endfor 


rz = rz[*,1:*] 
eps = [6.32,4.241,4.843] 

rz1 = rz[*,where(rz[1,*] eq 6.32)]
rz1 = rz1[*,sort(rz1[0,*])] 

kb = 1.3806488e-5        ;Boltzmann's constant (pN um/K)
kbt = kb * 293.15          ;Boltzmann energy (pN um)
ep0 = (8.854e-12)*(1.e-24)   ;vacuum permittivity (C^2/(pN*um))

kk = 1./(4.*!pi*eps*ep0)     ;Coulomb constant
e = 1.62e-19                 ;electron charge (C)
lamb = kk*(e^2.)/kbt     ;Bjerrum length (um)

a = reform(rz1[0,*])
  a_err = 0.03*a
kappa = rz1[6,*]
  kap_err = rz1[7,*]

akap = a*kappa
akap_err = akap*sqrt((kap_err/kappa)^2. + (a_err/a)^2.)

;a[12] = 1.6
;a[13] = 1.6

;x = max(a)*findgen(100)/99. 
;kap = mean(rz1[6,*])
;pot = mean(rz1[2,*])
;Z1 = (x + x^2.)*pot/(lamb[0]*kap)
;Z2 = x*pot/(lamb[0]*kap)

;  zmin= x*pot/(lamb[0]*(kap + stddev(rz1[6,*])))
;  zmax= x*pot/(lamb[0]*(kap - stddev(rz1[6,*])))
  
;  z_std = stddev(x*pot/(lamb[0]*rz1[6,*]))

xx = a*(1.+kappa*a) 
;x_err = xx*sqrt( (a_err/a)^2. + (kap_err/kappa)^2. )
xx_err = sqrt( ((1.+2.*kappa*a)*a_err)^2. + (a*a*kap_err)^2. )

nn = n_elements(a)
;akap_err = a_err*kap

pots = rz1[2,*]
  pots_err = rz1[3,*]

Z_eff = a*(1.+xx)*pots/lamb(0)
Zeff_err = Z_eff*sqrt( (pots_err/pots)^2. + (kap_err*a/(1.+xx))^2. + $
            (a_err*(pots+2.*a*pots*kap_err)/(a*pots*(1.+xx) ) )^2. ) 

; surface charge density
sigma = rz1[4,*]/(4.*!pi*a^2.)
zerr = rz1[5,*]
sig_err = sigma*sqrt( (zerr/rz1[4,*])^2. + (2.*a_err/a)^2.)

;----------- Unscreened Charges ----------------
rz2 = fltarr(4,1)

for i=0,n-1 do begin $
  rzi = reform(read_gdf(ss[i]))
  rz2 = [[rz2],[rzi]]
endfor

rz2 = rz2[*,1:*] 
rz2 = rz2[*,where(rz2[1,*] eq 5.84)]
rz2 = rz2[*,sort(rz2[0,*])]

;stop


;----------- chargefit-------------
expr = 'P[0]*X'
start = 550.

;pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0] },1)
;pi[1].limited[0]=1 & pi[1].limits[0]=0.0      ;setting P[1] to be positive

fit = mpfitexpr(expr, a, rz1[4,*], rz1[5,*], start, perror=pe, $
  bestnorm=chi, dof=dof) 

;---linear least squares error (Taylor Ch. 8)
nn = n_elements(a)
y_err = sqrt( total((rz1[4,*] - fit[0]*a)^2.)/(nn-2.) )   
; Taylor eq. (8.15)
delt = nn*total(a^2.) - total(a)^2.
; Taylor eq. (8.12)
fit_err = y_err*sqrt(nn/delt) 
; Taylor eq. (8.17)

pot1 = fit[0]*lamb[0]
  pot1_err = fit_err*lamb[0]
x = 1.1*max(xx)*findgen(100)/99. + 0.01
;aa = 1.5*findgen(100)/99.
Z2 = fit[0]*x
  zmin = float(fit[0]-2.*fit_err)*x
  zmax = float(fit[0]+2.*fit_err)*x
;print, fit[0]*lamb[0]*kap

;stop

;--------------chargefit 2--------------------------------------
expr2 = 'P[0]*X'
start2 = 550.

fit2 = mpfitexpr(expr, rz2[0,*], rz2[2,*], rz2[3,*], start2)

;nn = n_elements(a)
;fit_err = sqrt( total((rz1[4,*] - fit[0]*x)^2.)/(nn-2.) )

;pot1 = fit[0]*lamb[0]
;pot1_err = fit_err*lamb[0]
;x = 1.1*max(xx)*findgen(100)/99.

Z3 = fit2[0]*x
;zmin = float(fit[0]-fit_err)*x
;zmax = float(fit[0]+fit_err)*x

;-------------chargefit 3-------------------------------------
expr3 = 'P[0]*X'
start3 = 500.

fit3 = mpfitexpr(expr3, xx, rz1[4,*], rz1[5,*], start3)
zz3 = fit3[0]*x

;---linear least squares error (Taylor Ch. 8)
;nn = n_elements(xx)
y_err3 = sqrt( total((rz1[4,*] - fit3[0]*xx)^2.)/(nn-2.) )
; Taylor eq. (8.15)
delt3 = nn*total(xx^2.) - total(xx)^2.
; Taylor eq. (8.12)
fit_err3 = y_err3*sqrt(nn/delt3)
; Taylor eq. (8.17)

pot3 = fit3[0]*lamb[0]
pot3_err = fit_err3*lamb[0]
x = 1.1*max(xx)*findgen(100)/99. + 0.01
;aa = 1.5*findgen(100)/99.

zmin3 = float(fit3[0]-2.*fit_err3)*x
zmax3 = float(fit3[0]+2.*fit_err3)*x
;print, fit[0]*lamb[0]*kap


;---------------------------------------------------------------

marg = {font_name:'Times', xminor:4, yminor:1, sym_size:1, $
  sym_fill_color:'orange', font_size:22, sym_filled:1}
marg2 = {font_name:'Times', xminor:3, yminor:3, sym_size:3, $
  sym_color:'dark cyan', font_size:16, sym_filled:0,sym_thick:2}
marg3 = {font_name:'Times', xminor:3, yminor:3, sym_size:1, $
    sym_fill_color:'aquamarine', font_size:16, sym_filled:1}  

;----------------CHARGE PLOT-----------------------------------------

  
;pl_z1 = errorplot(xx, rz1[4,*], a_err, rz1[5,*], 'ko7', $
;  xtitle='$a(1+\kappa a)$ [$\mu$m]', ytitle='$Z*[a(1+\kappa a)]$', $
;  xtickinterval=1., ytickinterval=500, $
;  xminor=3, yminor=3, xrange=[0.02,4.7], aspect_ratio=4.7/1600, $
;  yrange=[0,1600], _extra=marg, $;name='$Z$ from force fits',$
;  errorbar_capsize=0.5)
;p1b = errorplot(a, Z_eff, rz1[0,*]*0.03, Zeff_err, $
;  'ktu4', sym_size=2, overplot=1,  name='$Z_{eff}$ from surface potential',$
;  sym_fill_color='deep sky blue', sym_filled = 1,errorbar_color='blue')  

;p1_err = polygon([x/kap, reverse(x/kap)],[zmin,reverse(zmax)], $
;  target=p1_z1, fill_transparency=60, overplot=1,$
;  /data, /fill_background, fill_color='medium violet red', linestyle=6)

;p1_err = polygon([x, reverse(x)],[zmin,reverse(zmax)], $
;  target=p1_z1, fill_transparency=60, overplot=1,$
;  /data, /fill_background, fill_color='medium violet red', linestyle=6)


;pl_z2 = plot(xx, rz1[4,*], 'ko4',overplot=1, _extra=marg)
;plc = plot(rz1[0,*], Z_eff, 'ktu4',overplot=1, $
;  sym_size=2, sym_fill_color='orange', sym_filled = 1)
;   p1_one = plot([xx[11],xx[11]], [rz1[4,11],rz1[4,11]],'ko5', $
;    _extra=marg2, /current, overplot=1)
   
;pl_z1f = plot(x/kap, Z1, 'k-2', /current, overplot = 1, $
;  name = 'fit for $\epsilon = 5.84$')
;  pl_z1f = plot(x, Z2, 'r-2', /current, overplot = 1, $
;    _extra=marg)

  

;--------------CHARGE PLOT 2-----------------------------------------
;----------------error-plot of Z* vs a--------------------------

p2_z1 = errorplot(a, rz1[4,*], a_err, rz1[5,*], 'ko7', $
  xtitle='$a$ [$\mu$m]', ytitle='$Z*(a)$', $
  xtickinterval=1., ytickinterval=500, $
  xminor=3, yminor=3, xrange=[0.02,2.8], aspect_ratio=2.8/1600, $
  yrange=[0,1600], _extra=marg, $;name='$Z$ from force fits',$
  errorbar_capsize=0.5)

;p2_z2 = errorplot(rz2[0,*], rz2[2,*], a_err, 5.*rz2[3,*], 'ko7', $
;  errorbar_capsize=0.5, overplot=1)

p2_err = polygon([x, reverse(x)],[zmin,reverse(zmax)], $
  target=p1_z1, fill_transparency=65, overplot=1,$
  /data, /fill_background, fill_color='medium violet red', linestyle=6)
;for i=0,nn-1 do begin $
;  p22 = plot(x, x*(1.+rz1[6,i]*x)*rz1[2,i]/lamb[0], '-3', /current, $
;    overplot=1, transparency=70, color='medium violet red')
;endfor

p2_z1f = plot(x, Z2, 'r-2', /current, overplot = 1, $
  _extra=marg)
;p2_z2f = plot(x, Z3, 'b--2', /current, overplot = 1, $
;    _extra=marg, name='Unscreened')  
  
  ;p2_mean = plot(x, x*(1.+mean(rz1[6,*])*x)*mean(rz1[2,*])/lamb[0], $
  ;  'r-3', /current, overplot=1)

;p2_z3 = plot(rz2[0,*], rz2[2,*], 'ko4',overplot=1, _extra=marg3)

p2_z1a = errorplot(a, rz1[4,*], a_err, rz1[5,*], 'ko7', overplot=1, $
  errorbar_capsize=0.5 )
p2_z4 = plot(a, rz1[4,*], 'ko4',overplot=1, _extra=marg)
p2_one = plot([a[10],a[10]], [rz1[4,10],rz1[4,10]],'ko5', $
  _extra=marg2, /current, overplot=1)

;  leg = legend(target=p2_z1f, position=[2.8,350], $
;    /data,font_name='Times', shadow=0, $
;    horizontal_spacing=0.05, font_size=14)

;----------plot of Z*\lambda/a vs a\kappa-----------------------

x = 0.25*findgen(1000)/999. + 0.02

surf = rz1[4,*]*lamb[0]/a
surf_err = surf*sqrt((rz1[5,*]/rz1[4,*])^2. + (a_err/a)^2.) 

p3 = errorplot(akap, surf, akap_err, surf_err,$
  _extra=marg, xrange=[0.02,0.25], yrange=[3,14],$
  xtickinterval=0.05, 'ko7', $
  ytickinterval=2, $
  xtitle = '$\kappa a$', ytitle= '$Z*(a)\lambda/a$')
  
mean_surf = fltarr(n_elements(x)) + mean(surf)

max_Z = -2.*alog(x) + 2.*alog(-alog(x)) + 4.*alog(2.) ;effective valence
  zet = mean(surf)
min_Z = raman(x,zet) 
;from Ramanathan (1988)
;err_max_Z = 4*akap_err/akap
;  p3_zz1 = polygon([x, reverse(x)],$
;    [min_Z,reverse(max_Z)], $
;    target=p3, fill_transparency=65, overplot=1,$
;    /data, /fill_background, fill_color='gold', linestyle=6)
p3_ZZ = plot(x, max_Z, 'k-1', overplot=1)
;p3_ZZ1 = plot(x, min_Z, 'k--1', overplot=1)

p3_zz2 = polygon([x, reverse(x)],$
  [mean_surf - 0.4,reverse(mean_surf + 0.4)], $
  target=p3, fill_transparency=65, overplot=1,$
  /data, /fill_background, fill_color='gold', linestyle=6)
p3_s = plot(x, mean_surf, 'b--1', overplot=1)

p3 = errorplot(akap, surf, akap_err, surf_err,$
  errorbar_capsize=0.3, 'ko7', overplot=1)  
p3_z1 = plot(akap, surf, 'ko10', _extra=marg, overplot=1) 

;----------plot of Z*\lambda/a vs a\kappa where \kappa = mean(kappa)-----

;surf = rz1[4,*]*lamb[0]/a
;surf_err = surf*sqrt((rz1[5,*]/rz1[4,*])^2. + (a_err/a)^2.)

;p3 = errorplot(a*mean(kappa), surf, akap_err, surf_err,$
;  _extra=marg, xrange=[0.02,0.25], yrange=[3,16],$
;  xtickinterval=0.05, 'ko7', $
;  ytickinterval=2, xtitle = '$\kappa a$', ytitle= '$Z*(a)\lambda/a$')

;mean_surf = fltarr(n_elements(x)) + mean(surf)
;p3_s = plot(x, mean_surf, 'b--2', overplot=1)

;max_Z = -2.*alog(x) + 2.*alog(-alog(x)) + 4.*alog(2.) ;effective valence
;from Ramanathan (1988)

;p3_z1 = plot(a*mean(kappa), surf, 'ko10', _extra=marg, overplot=1)


;--------------------------SURF POT PLOT-------------------------------
  
;p3 = errorplot(rz1[0,*], rz1[2,*], rz1[0,*]*0.03, rz1[3,*], _extra=marg, 'ko4',$
;  xtitle='a ($\mu$m)', ytitle='surface potential [$k_B T$]', $
;  sym_size=3, $
;  name='$\kappa$ from force fits', yrange=[3,8])
  
;  p3_one = plot([rz1[0,11],rz1[0,11]], [rz1[2,11],rz1[2,11]],'ko5', $
;  _extra=marg2, /current, overplot=1)  
;p3a = errorplot(a*kap, rz1[2,*], akap_err, rz1[3,*], $
;  'ktu4', sym_size=2, overplot=1, name='$\kappa = 0.13$ $\mu m^{-1}$',$
;    sym_fill_color='orange', sym_filled = 1)
;p3b = plot(x, 0.4/x, 'k-2', overplot=1)    



;p4 = errorplot(a, sigma, a_err, sig_err, _extra=marg, 'ko4',$
;  xtitle='a [$\mu m$]', ytitle='surface charge density [$Z/\mu m^2$]', $
;  name='$\kappa$ from force fits')
  ;p3a = errorplot(a*kap, rz1[2,*], akap_err, rz1[3,*], $
  ;'ktu4', sym_size=2, overplot=1, name='$\kappa = 0.13$ $\mu m^{-1}$',$
  ;sym_fill_color='orange', sym_filled = 1)
;p3b = plot(x, 0.4/x, 'k-2', overplot=1)

;=-----------------plot of Z* vs. a(1 + a\kappa)----------------

;p4 = errorplot(xx, rz1[4,*], xx_err, rz1[5,*], _extra=marg, 'ko4',$
;  xtitle='$a(1+\kappa a)$ [$\mu$m]', ytitle='$Z*(a)$',$
;  name='$Z$ from force fits',xrange=[0.02,3.25],xmajor=3, yrange=[0,1600],$
;   errorbar_capsize=0.5,xtickinterval=1, aspect_ratio=3.323/1600) 

  ;p4_err = polygon([x, reverse(x)],[zmin3,reverse(zmax3)], $
 ;   target=p4, fill_transparency=65, overplot=1,$
;    /data, /fill_background, fill_color='medium violet red', linestyle=6)
 
; p4_z1f = plot(x, zz3, 'r-2', /current, overplot=1)  
 
; p4_z4 = plot(xx, rz1[4,*], 'ko4',overplot=1, _extra=marg)
; p4_one = plot([xx[10],xx[10]], [rz1[4,10],rz1[4,10]],'ko5', $
;   _extra=marg2, /current, overplot=1) 

;p4a = errorplot(a*kap, rz1[4,*], akap_err, rz1[5,*], $
;    'ktu4', sym_size=2, overplot=1,  name='$\kappa = 0.13$ $\mu m^{-1}$',$
;    sym_fill_color='orange', sym_filled = 1)
; p4b = errorplot(xx, Z_eff, akap_err, Zeff_err, $
;    'ktu4', sym_size=2, overplot=1,  name='$Z_{eff}$ from surface potential',$
;    sym_fill_color='orange', sym_filled = 1)

;pl_z2 = errorplot(a, 1./kappa, a_err, kap_err/(kappa^2.), 'ko4', $
;  xtitle='a ($\mu m$)', ytitle='$\kappa^{-1}$ ($\mu m$)', _extra=marg) 

;pl_z3 = errorplot(rz3[0,*], rz3[3,*], rz3[4,*],'btu4', _extra=marg, $
;  name = '$\epsilon = 4.843$', /current, overplot=1) 
;pl_z3f = plot(a, Z3, 'b-2', /current, overplot = 1, $
;  name = 'fit for $\epsilon = 4.843$')   

;leg = legend(target=[pl_z1,p1b], font_size=16, $
;  position=[1.1,200],/data)
;leg1 = legend(target=[p3,p3a], font_size=16, $
  ;position=[0.19,8.8],/data) 
;leg2 = legend(target=[p4,p4b], font_size=16, $
;  position=[0.19,990],/data) 

end

