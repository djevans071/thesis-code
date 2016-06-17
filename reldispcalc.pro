;
;--RELDISPCALC--
; Tracks particles in a given blink and calculates the relative
; displacement and relative displacement difference over a given step
; size.
;
; SYNTAX:
;   RF_i = reldispcalc(posi, mpp, [lags=value])
;
; ARGUMENT:
;   posi: pretracked position data from PMMA_PRETRACK for the ith blink
;     where
;       posi[0,*]: x-position
;       posi[1,*]: y-position
;       posi[2,*]: frame number
;   mpp: microns per pixel (for scaling dimensions)
;
; KEYWORD PARAMETERS:
;   lags: specifies step size in time (default value of tau is 1 frame)
;   motion: specifies the type of motion to be tracked:
;           motion = 'rm' -- relative motion parallel and normal to line
;                            of separation
;           motion = 'cm' -- center-of-mass motion parallel and normal
;                            to line of separation
;   pair: specifies which pair of particles to track. In the form of a 
;         two elements, each with the label of the particle. For example:
;         pair = [0,1] = [1,0] is a pair of two particles where the first
;         one is particle 0 and the second one is particle 1 from track. 
;
; RETURN VALUE:
;   RF_i: array of tracked particle data for the ith blink where
;   RF_i[0,*]: relative displacement in microns
;   RF_i[1,*]: relative displacement difference by tau steps parallel
;              to the line of separation in microns
;   RF_i[2,*]: relative displacement difference by tau steps normal
;              to the line of separation in microns
;   RF_i[3,*]: frame number
;
;2013.07.02 -- D J Evans
;2014.09.12 -- Added MOTION keyword in order to measure all modes of motion
;              (2 types each of relative and center-of-mass motion)
;              - D J Evans

function reldispcalc, posi, $
  mpp=mpp, $
  lags = lags, $
  motion = motion
  ;pair = pair

  frame = 4    ;frame memory for track
  p = 10       ;maxdisp for track
  cull = 3     ;goodenough for track
  ;mpp=0.0685   ;microns per pixel
  thresh = 35  ;threshold for cutting out bad blinks

  if n_elements(lags) ne 1 then lags=1

  if (n_elements(posi[0,*]) le 2.*thresh) then return, [0,0,0,-1]
    ;if (n_elements(posi[0,*]) gt 2.*thresh) then begin $

    trackposi = track(posi, p, dim=2, mem=frame, good=cull,/quiet)  ;tracking data for particles in the ith blink

  label = reform(trackposi[3,*])

  x0y0 = trackposi[0:1, where(label eq pair[0])]*mpp    ;xy positions for particle 0
  ;if max(label) gt 1. then x0y0 = trackposi[0:1, where(label eq 2.)]*mpp
  x1y1 = trackposi[0:1, where(label eq pair[1])]*mpp    ;xy positions for particle 1
  ;if max(label) gt 1. then x1y1 = trackposi[0:1, where(label eq 3.)]*mpp

  rx = reform(x1y1[0,*] - x0y0[0,*])
  ry = reform(x1y1[1,*] - x0y0[1,*])
  cx = reform(x1y1[0,*] + x0y0[0,*])/2.
  cy = reform(x1y1[1,*] + x0y0[1,*])/2.

  ri = [transpose(rx), transpose(ry)]   ;relative position array
  rp = [-transpose(ry), transpose(rx)]  ;normal relative position array
  ci = [transpose(cx), transpose(cy)]   ;cm position array
  cp = [-transpose(cy), transpose(cx)]  ;normal cm position array

  nri = n_elements(rx)

  if (nri le thresh) then return, [0,0,0,-1]

  ri = ri[*,5:*:lags]
  rp = rp[*,5:*:lags]
  ci = ci[*,5:*:lags]
  cp = cp[*,5:*:lags]

  ;---------------RELATIVE POSITION CALCULATIONS----------------------------
  RDif_i = ri[*,1:*] - ri        ;position difference array
  R_i   = (ri[*,1:*] + ri)/2.    ;midpoint relative position array
  Rp_i  = (rp[*,1:*] + rp)/2.    ;midpoint normal relative position array
  n = n_elements(R_i[0,*])
  t = lags*findgen(n)

  modR_i = reform(sqrt(R_i[0,*]^2. + R_i[1,*]^2.))
  ;magnitude of the midpoint relative position vector R_i
  modRp_i = reform(sqrt(Rp_i[0,*]^2. + Rp_i[1,*]^2.))
  ;magnitude of the midpoint normal relative position vector Rp_i

  modRDif_i = total(RDif_i*R_i,1)/modR_i
  ;position difference array projected onto the direction
  ;of the relative position vector
  modRpDif_i = total(RDif_i*Rp_i,1)/modRp_i
  ;position difference array projected onto the direction
  ;of the normal relative position vector

  ;---------------CENTER OF MASS POSITION CALCULATIONS---------------------
  CDif_i = ci[*,1:*] - ci        ;position difference array
  C_i   = (ci[*,1:*] + ci)/2.    ;midpoint cm position array
  Cp_i  = (cp[*,1:*] + cp)/2.    ;midpoint normal cm position array
  n = n_elements(C_i[0,*])

  modC_i = reform(sqrt(C_i[0,*]^2. + C_i[1,*]^2.))
  ;magnitude of the midpoint cm position vector R_i
  modCp_i = reform(sqrt(Cp_i[0,*]^2. + Cp_i[1,*]^2.))
  ;magnitude of the midpoint normal cm position vector Rp_i

  modCDif_i = total(CDif_i*C_i,1)/modC_i
  ;position difference array projected onto the direction
  ;of the cm position vector
  modCpDif_i = total(CDif_i*Cp_i,1)/modCp_i
  ;position difference array projected onto the direction
  ;of the normal cm position vector

  ;window, 8 & plot, modR_i, ps=4
  ;window, 9 & plot, modRDif_i
  ;wait, 0.5
  ;
  ;stop
  rmaster = [[modR_i], [modRDif_i], [modRpDif_i], [t]]
  cmaster = [[modR_i], [modCDif_i], [modCpDif_i], [t]]

if motion eq 'rm' then return, transpose(rmaster)
if motion eq 'cm' then return, transpose(cmaster)
if n_elements(motion) ne 1 then return, transpose(rmaster)

end