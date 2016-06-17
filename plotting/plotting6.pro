;plotting6        ;produces representative potential plots
;
home1 = '../../2014.12.10 - PMMA1 in CHB/'
home2 = '../../2015.03.18 - PMMA3 in CHB/'
home5 = '../../2015.04.01 - PMMA5 in CHB/'

folder1 = 'trial 01'
folder2 = 'trial 02'
folder3 = 'trial 03'
folder4 = 'trial 04'

forces = '/forces.gdf'
;charges = '/charges.gdf'
;diffusivity = '/diffusivity.gdf'
;diff0 = '/diff0.gdf'
;charges_p = '/p_charges.gdf'

a = [0.55,1.6,2.7]
eps = 6.32
T = 296.
kb = 1.3806488e-5       ;Boltzmann's constant (pN um/K)
kbt = kb * T            ;Boltzmann energy (pN um)
n = 450

rf_1 = read_gdf(home1 + folder2 + forces)
rf_2 = read_gdf(home2 + folder1 + forces)
rf_5 = read_gdf(home5 + folder1 + forces)

rf = [[[rf_1]],[[rf_2]],[[rf_5]]]
rf = rf[*,0:n-1,*]


;rd1 = rd1[*,20:-20]

;rf1 = rf1[*,0:-40]
;rd1 = rd1[*,0:-40]


z_scr = fltarr(n,3)
z_ucr = fltarr(n,3)

uu_ = fltarr(n,3)
uerr_ = fltarr(n,3)

for i=0,2 do begin
  r1 = reform(rf[0,*,i]) &
  f1 = reform(rf[1,*,i]) & fe1 = reform(rf[2,*,i])

x = 50.*findgen(n)/(n-1)

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
  z1 = pot_fit(r,potential, a[i], eps=eps, x=0)     ;screened-Coulomb potential fit
  ;z2 = pot_fit2(reform(rf[*,*,i]), a[i], eps=eps)    ;unscreened-Coulomb potential fit

  z1_offset = min(z1)
  uu = uu + z1_offset
  potential = transpose([[uu],[u_err]])
  
    z1 = pot_fit(r, potential, a[i], eps=eps, x=x);screened-Coulomb potential fit
    z2 = pot_fit2(r, potential, a[i], eps=eps)    ;unscreened-Coulomb potential fit

 ; stop

 uu_[*,i] = uu
 uerr_[*,i] = u_err

  z_scr[*,i] = z1
  z_ucr[*,i] = z2
endfor

marg = {xtitle:'r [$\mu$m]', font_size:16, xminor:1, yminor:9, $
  xrange:[2,35], font_name:'Times'}
  
;-----------------Screened-POTENTIAL PLOT---------------------------
p = plot(rf[0,*,2], (r1)*uu_[*,2]/kbt, 'k1', ylog=1,$
    ytitle='$(r[\mu$m$])\times(U(r)/k_BT)$',yrange=[70,20000], $
    aspect_ratio=33/(alog10(20000)-alog10(70)),_extra=marg)

  (p.axes)[0].minor=1 & (p.axes)[2].minor=1
  (p.axes)[1].minor=9 & (p.axes)[3].minor=9

pf_scr = objarr(3)

for i=2,0,-1 do begin $

;(pot_p1.axes)[0].showtext=0
  r1 = reform(rf[0,*,i])
  uu = uu_[*,i]
  u_err = uerr_[*,i]

  p1 = plot(rf[0,*,i], (r1)*uu_[*,i]/kbt, 'k1', /current, $
    overplot=1)

  p_err = polygon([r1, reverse(r1)],[(r1)*(uu+u_err)/kbt, $
    reverse((r1)*(uu-u_err)/kbt)], $
    target=p1, fill_transparency=50, $
    /data, /fill_background, fill_color='gold', linestyle=6)
    
  if (i eq 0) then begin
    color = 'red'
    name = '0.55 $\mu$m'
    endif
    
  if (i eq 1) then begin
    color = 'blue'
    name = '1.60 $\mu$m'
    endif
    
    if (i eq 2) then begin
    color = 'green'
    name = '2.70 $\mu$m'
    endif

  z1 = z_scr[*,i]
  z2 = z_ucr[*,i]

  pf_scr[i] = plot(x,(x)*z1/kbt, color=color, $ 
    '-2', /current, overplot=1, font_name='Times',$
    name = name)
  pf_ucr = plot(r1,(r1)*z2/kbt, color=color, $
    '--2', /current, overplot=1, font_name='Times')    
  
  pot_p1 = plot(r1, (r1)*uu/kbt, 'k3', /current, $
    overplot=1, font_name='Times')
    
endfor

leg = legend(target=reverse(pf_scr),position=[16,8000], $
  /data,font_name='Times', shadow=0, $
  horizontal_spacing=0.05, font_size=17)


end