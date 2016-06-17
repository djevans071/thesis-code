;
;--FORCES--
; Calculates kernel density estimations of diffusivity and force between
; two charged particles as a function of their relative displacement using
; the method of Sainis, Germain, and Dufresne 2007 (PRL 99, 018303). 
;
; SYNTAX:
;   rfd = forces(rf1, rf2, fps, [cut=value])
;   
; ARGUMENTS:
;   rf1: Output of TRACK_PMMA where tau=1
;   rf2: Output of TRACK_PMMA where tau=2
;   fps: frames per second
; 
; KEYWORD PARAMETER:
;   cut: specifies how much data is randomly extracted for analysis 
;        (e.g. cut=1000 only analyzes 1000 random points from rf1 and rf2)
;        Default cut=0
; 
; RETURN VALUE:
;   rfd: force and diffusivity array where
;       rfd[0,*]: relative position of particles (microns)
;       rfd[1,*]: force (pN)
;       rfd[2,*]: diffusivity (microns^2/sec)
;       rdf[3,*]: relative position PDF
; 
; 2014.05.06 -- D J Evans
; 2014.09.12 -- D J Evans - cleanup


function forces, diff   ;output of diffusivity

  ;fps = 1000  (now an input)
  kb = 1.3806488e-5     ;Boltzmann's constant (pN um/sec K)
  kbt = kb * 293.15

  r   = reform(diff[0,*])     ;separation
  v   = reform(diff[3,*])     ;velocity
  v_e = reform(diff[4,*])     ;error in velocity
  d   = reform(diff[1,*])     ;diffusivity
  d_e = reform(diff[2,*])     ;error in diffusivity
  posdist = reform(diff[5,*])                              
     
    ; stop
                              
  f = kbt * (v/d) / sqrt(3.)                      ;forces
  f_err = abs(f)*sqrt( (v_e/v)^2. + (d_e/d)^2. )  ;error in force
  
  return, [transpose(r), $
           transpose(f),  transpose(f_err), $
           transpose(posdist)]
end   