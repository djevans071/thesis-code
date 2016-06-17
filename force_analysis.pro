function force_analysis, home, folder, a,$
                          fpb, mpp, $
                          num=num

if not keyword_set(num) then num = 2
;@force_analysis

;home = '../2014.11.04 - PMMA2.3 in CHB/'  
;folder = 'trial 02'  
dir = home + folder  

;stop

fps = 1500                  ;frames per second
;fpb = 130                   ;frames per blink (now an input)
;mpp = 0.103                 ;mpp for 60x (now an input)
;mpp = 0.0685               ;mpp for 100x

eps = 5.93                 ;dielectric constant of medium
;a = 1.13                     ;particle radius (now an input)

x0 =  30  &  x1 = 430
y0 = 50  &  y1 = 500

window = [x0,x1,y0,y1]    ;viewing region for pretrack (in pixels)
diam   = 61              ;pixel diameter input for feature in pretrack


;bg = bg(dir)  
;print, 'background saved'  
;pretrack = pmma_pretrack(dir, fpb, diam=diam, window=window, num=num)  
pretrack = read_gdf(dir + '/pretrack.gdf')

print, dir

pairs = combigen(num,2)       ;creates pairs (0, 1), (0, 2), (1, 2)
     ;for num = 3 and (0, 1) for num = 2
n_pairs = n_elements(pairs)/2

rf1 = track_pmma(dir, mpp=mpp, motion='rm', num=num)
rf2 = track_pmma(dir, mpp=mpp, motion='rm', num=num, tau=2)
rf3 = track_pmma(dir, mpp=mpp, motion='rm', num=num, tau=3)
;stop
;  for j=0,n_pairs-1 do begin $
;    pair = reform(pairs[j,*])
;    print, pair
;    rf1 = track_pmma(dir, mpp=mpp, motion='rm', pair=pair, num=num)  
;    rf2 = track_pmma(dir, mpp=mpp, tau=2, motion='rm', pair=pair, num=num)  
;    rf3 = track_pmma(dir, mpp=mpp, tau=3, motion='rm', pair=pair, num=num)


  


;--------------------CALCULATE DIFFUSIVITIES----------------------------

 ; pair_string = 'pair' + strtrim(pair[0],2) + strtrim(pair[1],2)
  diff1_p = diffusivity(rf1[0:1,*],rf2[0:1,*],rf3[0:1,*],fps,mpp)        ;rm parallel
    ;diff1_p = diff1_p[*,15:-15]  
 ; diff1_n = diffusivity(rf1[0:2:2,*],rf2[0:2:2,*],rf3[0:2:2,*], fps, mpp)   ;rm normal
  write_gdf, diff1_p, dir + '/diffusivity.gdf'  
    print, 'diffusivity.gdf recorded'

  ff = forces(diff1_p)  
  ;ff_n = forces(diff1_n)
  ;ff_n = forces(diff1_n)
    ;ff = ff[*,20:-20]  
    ;stop
  write_gdf, ff, dir + '/forces.gdf'  
  ;write_gdf, ff_n, dir + '/forces_n_' + pair_string + '.gdf'  
    print, 'forces.gdf recorded'
    
  ;endfor
  ;stop

  r1_p = diff1_p[0,*]   
    ;r1_n = diff1_n[0,*]  
  ;r2_p = diff2_p[0,*] & r2_n = diff2_n[0,*]      
  d1_p = diff1_p[1,*]  
    ;d1_n = diff1_n[1,*]  
  ;d2_p = diff2_p[1,*] & d2_n = diff2_n[1,*]  
  de1_p = diff1_p[2,*]  
    ;de1_n = diff1_n[2,*]  
  ;de2_p = diff2_p[2,*] & de2_n = diff2_n[2,*]  

  ;------------------VELOCITY PLOT-----------------------------------

  v = diff1_p[3,*]  ;& v_n = diff1_n[3,*] 
  ve = diff1_p[4,*]  ;& ve_n = diff1_n[4,*] 
  ;window, 2 & plot, r1_p, v, ps=circ()  
  ;  oplot, r1_p, v + ve, line=4  
  ;  oplot, r1_p, v - ve, line=4  
  
  
  
  ;-------------------DIFFUSIVITY PLOTS---------------------------------
  
  window, 0 & plot, r1_p, d1_p, ps=circ()  
    ;oplot, r1_p, d1_p, ps=5  
    oplot, r1_p, d1_p + de1_p, line=1  
    oplot, r1_p, d1_p - de1_p, line=1  

  y1_p = diff_fit(diff1_p,a,motion='rm_p', folder=dir)  
    d0 = read_gdf(dir + '/diff0.gdf')  
    d0 = fltarr(n_elements(r1_p)) + d0[0]
    oplot, r1_p, y1_p, line=3, thick=3, color=150
    oplot, r1_p, 2.*d0, color=150, line=4 
   ; oplot, r1_n, y1_n, line=4  

  ;stop

;-------------------FORCE PLOT-----------------------------------------
  ;ff = ff[*,40:-40]
  
  r = reform(ff[0,*])  
  f = reform(ff[1,*])  
    fe = ff[2,*]  
    z = force_fit(ff,a, eps=eps, folder=dir) 
    ;z2 = force_fit2(ff,a, eps=eps, folder=dir)

  
    ;u = pot_fit(ff,a, eps=eps, folder=dir)  
    ;u2 = pot_fit2(ff,a, eps=eps, folder=dir)
    
  window, 1 & plot, r, f, ps=circ()  
    oplot, r, f + fe, line=4  
    oplot, r, f - fe, line=4  
    oplot, r, z, thick=3, color=150


  
;stop   
return, 1 

end