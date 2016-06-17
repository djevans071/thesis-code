;plotting7      ;produces full spectrum potential plots
;
s = file_search('../../*/*/forces.gdf', count=fn)
ss = file_search('../../*/*/charge.gdf', count=fn) 

rf = fltarr(4, 500, 1)
rz = fltarr(8,1)

for i=0,fn-1 do begin
  rfi = read_gdf(s[i])
  rf  = [[[rf]],[[rfi]]]
  
  rzi = read_gdf(ss[i])
  rz  = [[rz],[rzi]]
endfor
rf = rf[*,*,1:*]
rz = rz[*,1:*]

;home1 = '../../2014.12.10 - PMMA1 in CHB/'
;home2 = '../../2015.03.18 - PMMA3 in CHB/'
;home5 = '../../2015.04.01 - PMMA5 in CHB/'

;folder = 'trial 01'
;folder2 = 'trial 02'

;forces = '/forces.gdf'
;charges = '/charges.gdf'
;diffusivity = '/diffusivity.gdf'
;diff0 = '/diff0.gdf'
;charges_p = '/p_charges.gdf'

a = reform(rz[0,*])
eps = 5.93
T = 293.15
kb = 1.3806488e-5       ;Boltzmann's constant (pN um/K)
kbt = kb * T            ;Boltzmann energy (pN um)
;n = 500

;rf_1 = read_gdf(home1 + folder2 + forces)
;rf_2 = read_gdf(home2 + folder + forces)
;rf_5 = read_gdf(home5 + folder2 + forces)

;rf = [[[rf_1]],[[rf_2]],[[rf_5]]]
;rf = rf[*,0:n-1,*]
m = n_elements(a)
n = n_elements(rf[0,*,0])

z_scr = fltarr(n,m)
z_ucr = fltarr(n,m)

uu_ = fltarr(n,m)
uerr_ = fltarr(n,m)

;stop

for i=0,m-1 do begin
  r1 = reform(rf[0,*,i]) &
  f1 = reform(rf[1,*,i]) & fe1 = reform(rf[2,*,i])
  
  ;---------------FORCE FITTING------------------------------
  
  ;finding the potential energy
  r = reform(rf[0,*,i])
  uu = (r1[1]-r1[0])*reverse(total(reverse(f1),/cum))
  uferr = (r1[1]-r1[0])*reverse(total(reverse(f1 + fe1), /cum))
  u_err = uferr - uu
  
  potential = transpose([[uu],[u_err]])
  ;fitting
  ;y1 = force_fit(rf1, a, eps=eps)   ;screened-Coulomb force fit
  ;y2 = force_fit2(rf1, a, eps=eps)   ;unscreened-Coulomb force fit
  z1 = pot_fit(r, potential, a[i], eps=eps)     ;screened-Coulomb potential fit
  ;z2 = pot_fit2(reform(rf[*,*,i]), a[i], eps=eps)    ;unscreened-Coulomb potential fit
  
  z1_offset = min(z1)
  uu = uu + z1_offset
  potential = transpose([[uu],[u_err]])
  
  z1 = pot_fit(r, potential, a[i], eps=eps)     ;screened-Coulomb potential fit
  z2 = pot_fit2(r, potential, a[i], eps=eps)    ;unscreened-Coulomb potential fit
  
  uu_[*,i] = uu
  uerr_[*,i] = u_err
  
  z_scr[*,i] = z1
  z_ucr[*,i] = z2
  ;stop
endfor

ar = 25./1600.

marg = {xtitle:'r [ $\mu$m]', font_size:16, xminor:1, yminor:4, $
  xrange:[4,25], font_name:'Times'}
  
;-----------------Screened-POTENTIAL PLOT---------------------------
p = plot(rf[0,*,m-1], uu_[*,m-1]/kbt, 'k1', ylog=0,$
  ytitle='$U(r)/k_BT$',yrange=[-10,1600], $
  _extra=marg, aspect_ratio=ar)
  
(p.axes)[0].minor=1 & (p.axes)[2].minor=1
(p.axes)[1].minor=4 & (p.axes)[3].minor=4

for i=0,m-1 do begin

  ;(pot_p1.axes)[0].showtext=0
  r1 = reform(rf[0,*,i])
  uu = uu_[*,i]
  u_err = uerr_[*,i]
  
  p1 = plot(rf[0,*,i], uu_[*,i]/kbt, 'k1', /current, $
    overplot=1)
    
  p_err = polygon([r1, reverse(r1)],[(uu+u_err)/kbt, $
    reverse((uu-u_err)/kbt)], $
    target=p1, fill_transparency=70, $
    /data, /fill_background, fill_color='gold', linestyle=6)
  

  ;pf_ucr_1 = plot(r1,z2/kbt, color=color, $
  ;  '--1', /current, overplot=1, font_name='Times')
    
  ;pot_p1 = plot(r1, uu/kbt, color=color, '-2', /current, $
  ;  overplot=1, font_name='Times')
    
endfor

color_name = strarr(m,2)
pot_p1 = objarr(m)

for i=0,m-1 do begin
  
  if (i ge 0 and i le 2) then begin
    color = 'deep sky blue' 
    name = '1.13 $\mu$m' 
  endif
  
  if (i ge 3 and i le 5) then begin
    color = 'red'
    name = '0.55 $\mu$m'
  endif
  
  if (i ge 6 and i le 9) then begin
    color = 'orange'
    name = '0.65 $\mu$m'
  endif
  
  if (i ge 10 and i le 11) then begin
    color = 'yellow green'
    name = '0.90 $\mu$m'
  endif
  
  if (i ge 12 and i le 13) then begin 
    color = 'medium violet red'
    name = '1.60 $\mu$m'
  endif
  
  if (i ge 14 and i le 15) then begin
    color = 'indigo'
    name = '2.70 $\mu$m'
  endif
  
  r1 = reform(rf[0,*,i])
  uu = uu_[*,i]
  
  z1 = z_scr[*,i]
  ;z2 = z_ucr[*,i]
  
  pf_scr_1 = plot(r1,z1/kbt, $
    'k--1', /current, overplot=1, font_name='Times')
  
  pot_p1[i] = plot(r1, uu/kbt, color=color, '-2', /current, $
    overplot=1, font_name='Times', name=name)
  
endfor

leg = legend(target=reverse([pot_p1[3],pot_p1[6],pot_p1[10],$
   pot_p1[0], pot_p1[12], pot_p1[14]]), $
  position=[26,1500], $
  /data,font_name='Times', shadow=0, $
  horizontal_spacing=0.05, font_size=16)
p.close

;Tuinier fitting--------------
;finding the potential energy
j=14

f_t = reform(rf[1,*,j]) & fe1_t = reform(rf[2,*,j])
r_t = reform(rf[0,*,j])
uu_t = uu_[*,j]
uerr_t = uerr_[*,j]

potential_t = transpose([[uu_t],[uerr_t]])
z1_t = pot_fit_t(r_t, potential_t, a[j], eps=eps)     ;screened-Coulomb potential fit
;z2 = pot_fit2(reform(rf[*,*,i]), a[i], eps=eps)    ;unscreened-Coulomb potential fit
;z1_offset_t = min(z1_t)
;uu_t = uu_t + z1_offset_t

;z1_offset = min(z1)
;uu = uu + z1_offset
;potential = transpose([[uu],[u_err]])

;z1 = pot_fit(r, potential, a[i], eps=eps)     ;screened-Coulomb potential fit

;uu_[*,i] = uu
;uerr_[*,i] = u_err

;z_scr[*,i] = z1




;------------single potential plot



rr0 = reform(rf[0,*,j])

  p0 = plot(rf[0,*,j], uu_[*,j]/kbt, 'k2', ylog=0,$
    ytitle='$U(r)/k_BT$', $
    _extra=marg, xrange=[4,25])
  p_err = polygon([reform(rf[0,*,j]), reverse(reform(rf[0,*,j]))],[(uu_[*,j]+uerr_[*,j])/kbt, $
    reverse((uu_[*,j]-uerr_[*,j])/kbt)], $
    target=p0, fill_transparency=70, $
    /data, /fill_background, fill_color='gold', linestyle=6)
  pf_scr_1 = plot(rf[0,*,j],z_scr[*,j]/kbt, $
    'r2', /current, overplot=1, font_name='Times')
  p0 = plot(rf[0,*,j], uu_[*,j]/kbt, 'k2', ylog=0,$
    /overplot)
   p0a = plot(rf[0,*,j], z1_t/kbt, 'b--2', /overplot)
  
    
  (p0.axes)[2].hide=1
  (p0.axes)[3].hide=1   

  pp0 = plot(rf[0,*,j], rr0*uu_[*,j]/kbt, 'k1', ylog=1,$
    ytitle='$y$', $
    ;aspect_ratio=33/(alog10(10000)-alog10(70)),
    _extra=marg, /CURRENT, POSITION=[.45,.4,.90,.85])
    
  pp0_err = polygon([reform(rf[0,*,j]), reverse(reform(rf[0,*,j]))],[rr0*(uu_[*,j]+uerr_[*,j])/kbt, $
    reverse(rr0*(uu_[*,j]-uerr_[*,j])/kbt)], $
    target=pp0, fill_transparency=70, $
    /data, /fill_background, fill_color='gold', linestyle=6)
  ppf_scr_1 = plot(rf[0,*,j],rr0*z_scr[*,j]/kbt, $
    'r2', /current, overplot=1, font_name='Times')
  ;ppf_ucr_1 = plot(rf[0,*,0],rr0*z_ucr[*,0]/kbt, $
  ;  'r--2', /current, overplot=1, font_name='Times')
  pp0 = plot(rf[0,*,j], rr0*uu_[*,j]/kbt, 'k2',$
    /current,/overplot)
    pp0a = plot(rf[0,*,j], rr0*z1_t/kbt, 'b--2', /current, /overplot)
  
    
end
