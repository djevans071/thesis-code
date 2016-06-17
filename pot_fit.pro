 function pot_fit, r, potential, a, $
                          eps=eps, $
                          folder=folder,$
                          x=x

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

  ;r = reform(force[0,*])
  ;f = reform(force[1,*])
  ;error = reform(force[2,*])

  ;uu = (r[1]-r[0])*reverse(total(reverse(f),/cum))
  ;uferr = (r[1]-r[0])*reverse(total(reverse(f + error), /cum))
  ;u_err = uferr - uu
 
uu = reform(potential[0,*])
u_err = reform(potential[1,*]) 
  
;stop

start = [100,0.1,-0.3]        
expr = 'P[0]*exp(-P[1]*X)/X + P[2]'

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0] },3)
  pi[0].limited[0]=1 & pi[0].limits[0]=0.0      ;setting P[0] to be positive
  pi[1].limited[0]=1 & pi[1].limits[0]=0.0      ;setting P[1] to be positive
  ;pi[2].limited[0]=1 & pi[2].limits[0]=0.0      ;setting P[2] to be positive

FIT = MPFITEXPR(expr, r, uu, u_err, start, parinfo=pi, perror=pe, $
  bestnorm=chi, dof=dof)  

if n_elements(x) gt 1 then y = FIT[0]*exp(-FIT[1]*x)/x; + FIT[2]
;if n_elements(x) eq 1 then y = FIT[0]*exp(-FIT[1]*r)/r; + FIT[2]
;y = FIT[0]*exp(-FIT[1]*r)/r + FIT[2]
y = FIT[0]*exp(-FIT[1]*r)/r

;y = FIT[0]*(FIT[1]/r)^2.

;pot = sqrt(lamb*FIT[0]/kbt)         ;surface potential               
pot = (exp(-FIT[1]*a)/a)*sqrt(lamb*FIT[0]/kbt)

;Z = FIT[1]*pot/lamb            ;number of charges
Z = (1.+FIT[1]*a)*exp(-FIT[1]*a)*sqrt(FIT[0]/(lamb*kbt) )

a_err=0.1
T_err=1.
ka = FIT[1]*a
ka1 = 1. + ka
F_err = pe*sqrt(chi/dof)

;pot_err = 0.5*pot*pe[0]/FIT[0]
;Z_err = Z*sqrt( (pot_err/pot)^2. + (a_err/a)^2. )

pot_err = pot*sqrt( (F_err[0]/(2.*FIT[0]))^2. + (T_err/(2.*T))^2. $
  + (a_err*ka1/a)^2. + (F_err[1]*a)^2. )

Z_err = Z*sqrt( (F_err[0]/(2.*FIT[0]))^2. + (T_err/(2.*T))^2. $
    + (ka*FIT[1]*a_err/ka1)^2. + (ka*a*F_err[1]/ka1)^2. )
    
  ;-----goodness of fit-------
  q = n_elements(FIT)
  n = n_elements(r)
  
  AIC = chi + 2.*q/(1.-(q-1.)/n)
  BIC = chi + q*alog(n)
      

  print, 'Fit Parameters'
  print, '     [U[0]       kappa        kappa error]'
  print, [FIT[0], FIT[1], F_err[1]] 
  print, '    [pot          Z]'
  print, [pot, Z]
  print, '   [pot error       Z error     bestnorm]'
  print, [pot_err, Z_err, chi]
  print, '    [AIC      BIC]'
  print, [AIC, BIC]
;window, 2
;  plot, r, uu, ps=circ()
;    oplot, r, uu+u_err, line=2
;    oplot, r, uu-u_err, line=2
;  oplot, r, y, thick=2, color=150


 ;write_gdf, [a, eps, pot,pot_err, Z, Z_err, FIT[1], F_err[1]], $
 ;   folder + '/p_charge.gdf'
;stop

return, y


end
