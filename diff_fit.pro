function diff_fit, diff, a, $
                   motion = motion, $
                   folder = folder  
   ;motion keyword: motion = 'rm_p' for RM parallel to axis of separation
   ;                motion = 'rm_n' for RM normal to axis of separation
   ;                motion = 'cm_p' for CM parallel to axis of separation
   ;                motion = 'cm_n' for CM normal to axis of separation

;a = 0.9 & $
diff0 = 0.13  ;single particle diffusion estimate (um^2/sec)
kb = 1.3806488e-5
kbt = kb * 293.15
eta = 1.8e-3

r = reform(diff[0,*])
d = reform(diff[1,*])
error = reform(diff[2,*])

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0] },2)
pi[0].limited[0]=1 & pi[0].limits[0]=0.      ;setting P[0] to be positive
;pi[1].limited[0]=1 & pi[1].limits[0]=0.      ;setting P[1] to be positive
pi[1].fixed=1                                ;setting P[1] to be a

if not keyword_set(motion) then motion = 'rm_p'

;---------------------RELATIVE MOTION DIFFUSION-------------------

if (motion eq 'rm_p') then begin $      ;parallel to axis of separation
  start = [diff0, a]
  diff_expr = '(P[0]*2.)*(1. - (3./2.)*(P[1]/X)' + $
     ' + (P[1]/X)^3. - (15./4.)*(P[1]/X)^4.)'
  fit = mpfitexpr(diff_expr,r,d,error,start, parinfo=pi)       ;diffusion fit for all pts
  diff = (fit[0]*2.)*(1. - (3./2.)*(fit[1]/r) + $
    (fit[1]/r)^3. - (15./4.)*(fit[1]/r)^4.)
endif

if (motion eq 'rm_n') then begin $      ;normal to axis of separation
  start = [diff0, a]
  diff_expr = '(P[0]*2.)*(1. - (3./4.)*(P[1]/X)' + $
    ' - (1./2.)*(P[1]/X)^3.)'
  fit = mpfitexpr(diff_expr,r,d,error,start, parinfo=pi)       ;diffusion fit for all pts
  diff = (fit[0]*2.)*(1. - (3./4.)*(fit[1]/r) - $
    (1./2.)*(fit[1]/r)^3.)
endif  

;---------------------CENTER OF MASS MOTION DIFFUSION--------------

if (motion eq 'cm_p') then begin $      ;parallel to axis of separation
  start = [diff0, a]
  diff_expr = '(P[0]*2.)*(1. + (3./2.)*(P[1]/X)' + $
    ' - (P[1]/X)^3. - (15./4.)*(P[1]/X)^4.)'
  fit = mpfitexpr(diff_expr,r,d,error,start, parinfo=pi)       ;diffusion fit for all pts
  diff = (fit[0]*2.)*(1. + (3./2.)*(fit[1]/r) - $
    (fit[1]/r)^3. - (15./4.)*(fit[1]/r)^4.)
endif

if (motion eq 'cm_n') then begin $      ;normal to axis of separation
  start = [diff0, a]
  diff_expr = '(P[0]*2.)*(1. + (3./4.)*(P[1]/X)' + $
    ' + (1./2.)*(P[1]/X)^3.)'
  fit = mpfitexpr(diff_expr,r,d,error,start, parinfo=pi)       ;diffusion fit for all pts
  diff = (fit[0]*2.)*(1. + (3./4.)*(fit[1]/r) + $
    (1./2.)*(fit[1]/r)^3.)
endif

;write_gdf, fit, folder + '/diff0.gdf'
;stop

print, [fit[0], fit[1]]  
return, diff

end