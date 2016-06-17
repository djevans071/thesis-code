function pot_fit2, r, potential, a, $
  eps=eps, $
  folder=folder
  
  ;a = 0.9                ;radius of particles
  T = 293.15
  kb = 1.3806488e-5       ;Boltzmann's constant (pN um/K)
  kbt = kb * T            ;Boltzmann energy (pN um)
  ;stop
  ;force = force[*,20:-20]
  
  if n_elements(eps) ne 1 then eps = 6.32
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

start = [100,-0.3]
expr = 'P[0]/X + P[1]'

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0] },2)
pi[0].limited[0]=1 & pi[0].limits[0]=0.0      ;setting P[0] to be positive
;pi[1].limited[0]=1 & pi[1].limits[0]=0.0      ;setting P[1] to be positive
;pi[2].limited[0]=1 & pi[2].limits[0]=0.0      ;setting P[2] to be positive

FIT = MPFITEXPR(expr, r, uu, u_err, start, parinfo=pi, perror=pe, $
  bestnorm=chi, dof=dof)
y = FIT[0]/r + FIT[1]

  ;number of charges
Z = sqrt(FIT[0]/(lamb*kbt))

a_err=0.1
T_err=1.
F_err = pe*sqrt(chi/dof)
  
Z_err = Z*sqrt( (F_err[0]/(2.*FIT[0]))^2. + (T_err/(2.*T))^2.)
  
  ;-----goodness of fit-------
q = n_elements(FIT)
n = n_elements(r)

AIC = chi + 2.*q/(1.-(q-1.)/n)
BIC = chi + q*alog(n)
  
print, 'Fit Parameters'
print, '     [U[0]       y-offset]'
print, [FIT[0], FIT[1]]
print, '     [Z          Z error]'
print, [Z, Z_err, chi]
print, '    [AIC      BIC]'
print, [AIC, BIC]
;window, 2
;plot, r, uu, ps=circ()
;oplot, r, uu+u_err, line=2
;oplot, r, uu-u_err, line=2
;oplot, r, y, thick=2, color=150


;write_gdf, [a, eps, Z, Z_err], $
;  folder + '/p_charge2.gdf'
  
return, y

end
