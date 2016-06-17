;
;--BAYES--
; Calculates the kernel density estimators necessary to find the 
; diffusivity (+ error) with respect to relative displacement
; r between two charged particles in a dielectric fluid.
;
; SYNTAX:
;   rve = bayes(RF,fps)
;   
; ARGUMENTS:
;   RF: Output of TRACK_PMMA. Tracked position and position difference
;       data where
;       RF[0,*]: relative position of particles
;       RF[1,*]: relative difference in position by a given time step
;       RF[2,*]: frame number
;   fps: frames per second
;   
; RETURN VALUE:
;   rve: data to feed to FORCES to find the force profile where
;        rve[0,*]: relative displacement (in microns)
;        rve[1,*]: average velocity (microns/sec)
;        rve[2,*]: error in average velocity (microns/sec)
;        rve[3,*]: variance of velocity (microns^2)
;        rve[4,*]: error in variance of velocity (microns^2) 
;        rve[5,*]: relative displacement probability distribution function
;        
; 2013.10.07 -- D J Evans
; 2014.09.12 -- D J Evans (used kde to bin mean velocity) 
; 2015.04.29 -- D J Evans (simplified and speeded up the calculation of dEta) 

function bayes, RF, x, fps, mpp              ; x from poscut

R_ = reform(RF[0,*])                ;separation between sphere centers [um]                  
v_ = reform(RF[1,*])*fps            ;radial velocity at each measurement [um/s]
;t_ = reform(RF[3,*]) 
    
dV = kde(R_,x, weight=v_, sigma=ds)  ;unnormalized average velocity flux
d  = kde(R_,x, sigma=s)              ;position pdf

v = reform(dV/d)           ;average drift velocity
verr = abs(v)*sqrt((ds/dV)^2. + (s/d)^2. )    ;statistical error in v

dx = 0.14*mpp                ;measurement error in x
de = err_kde(R_,x,dx)        ;error in the position pdf
dve = err_kde(R_,x,dx,weight=v_) ;error in the unnormed mean vel. flux

ve = abs(v)*sqrt( (dve/dV)^2 + (de/d)^2 ) ;measurement error in v
verr = sqrt(verr^2 + ve^2)    ;total error in v


;---------variance of velocity
gam = kde(R_, x, weight=v_^2., sigma=dgam)
dgame = err_kde(R_,x,dx,weight=v_^2.)
dEta = gam/d - v^2. 

;dEtaerr  = sqrt( (det/dEta)^2. + (s/d)^2. )
dEtaerr1 = sqrt( (gam*s/d^2.)^2. + (dgam/d)^2. $
                + (2.*v*verr)^2. )
dEtaerr2 = sqrt( (dgame/d)^2 +  (gam*de/d^2)^2 )
dEtaerr = sqrt(dEtaerr1^2 + dEtaerr2^2)                

Eta  = (1./fps^2.)*dEta
Etaerr = (1./fps^2.)*dEtaerr

;stop

return, [transpose(x), $
         transpose(v), transpose(verr), $
         transpose(Eta), transpose(Etaerr), $
         transpose(d)]
end