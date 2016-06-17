function force_fit2, force, a, $
  eps=eps, $
  folder=folder
  
  ;a = 0.9                ;radius of particles
  T = 293.15
  kb = 1.3806488e-5       ;Boltzmann's constant (pN um/K)
  kbt = kb * T            ;Boltzmann energy (pN um)
  ;stop
  ;force = force[*,20:-20]
  
  if n_elements(eps) ne 1 then eps = 5.84
  ;dielectric constant
  ep0 = (8.854e-12)*(1.e-24)  ;vacuum permittivity (C^2/(pN*um))
  kk = 1./(4.*!pi*eps*ep0)    ;Coulomb constant
  e = 1.62e-19                ;electron charge (C)
  lamb = kk*(e^2.)/kbt        ;Bjerrum length (um)
  
  r = reform(force[0,*])
  f = reform(force[1,*])
  error = reform(force[2,*])
      
  start = [100,0.15]
  expr = 'P[0]/X/X + P[1]'
  
  
  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0] },2)
  pi[1].limited[0]=1 & pi[1].limits[0]=0.0      ;setting P[1] to be positive
  ;pi[2].fixed=1                                ;setting P[2] to be a
  
  FIT = MPFITEXPR(expr, r, f, error, start, parinfo=pi, $
    bestnorm=chi, perror=pe, dof=dof)
  y = FIT[0]/r/r + FIT[1]
  
  ;----surface potential-----
  ;pot = sqrt(FIT[0]*lamb/kbt)*exp(-FIT[1]*a)/(FIT[1]*a)
        
  ;----number of charges (assuming constant charge)-----
  ;Z = sqrt(FIT[0]/(lamb*kbt))*exp(-FIT[1]*a)*(1.+FIT[1]*a)/FIT[1]
  Z = sqrt(FIT[0]/(lamb*kbt))
  
  
  a_err=0.1
  T_err=1.
  ;ka = FIT[1]*a
  ;ka1 = 1. + ka
  ;ka2 = 1. + 2*ka
  F_err = pe*sqrt(chi/dof)
  
  ;pot_err = pot*sqrt( (pe[0]/(2.*FIT[0]))^2. + (T_err/(2.*T))^2. $
  ;      + (a_err*ka2/a)^2. + (pe[1]*ka2/FIT[1])^2. )
    
  Z_err = Z*sqrt( (F_err[0]/(2.*FIT[0]))^2. + (T_err/(2.*T))^2. )
   
  ;-----goodness of fit-------
  
  q = n_elements(FIT)
  n = n_elements(r)
  
  AIC = chi + 2.*q/(1.-(q-1.)/n)
  BIC = chi + q*alog(n)
  
  
  print, 'Fit Parameters'
  print, '     [F[0]       y-offset]'
  print, [FIT[0], FIT[1]] 
  print, '     [Z          Z error]'
  print, [Z, Z_err]
  print, '    [AIC      BIC]'
  print, [AIC, BIC]  
  ;stop
  
;  write_gdf, [a, eps, Z, Z_err], $
;    folder + '/charge2.gdf'
  
  return, y
  
end
