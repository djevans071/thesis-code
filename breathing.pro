;
;--BREATHING--
; Tracks particles in a given blink and calculates the radial displacement
; and radial displacement difference of the breathing mode of three 
; particles (the effective radius of a circle that circumscribes the particles).
;
; SYNTAX:
;   RF_i = reldispcalc(posi, mpp, [lags=value])
;
; ARGUMENTS:
;   posi: pretracked position data from PMMA_PRETRACK for the ith blink
;     where
;       posi[0,*]: x-position
;       posi[1,*]: y-position
;       posi[2,*]: frame number
;   mpp = microns per pixen
;   
; KEYWORD PARAMETERS:
;   lags: specifies step size in time (default value of tau is 1 frame)
;   pair: specifies which pair of particles to track. In the form of a
;         two elements, each with the label of the particle. For example:
;         pair = [0,1] = [1,0] is a pair of two particles where the first
;         one is particle 0 and the second one is particle 1 from track.
;
; RETURN VALUE:
;   RF_i: array of tracked particle data for the ith blink where
;   RF_i[0,*]: effective relative displacement in microns
;   RF_i[1,*]: effective relative displacement difference by tau steps
;   RF_i[2,*]: frame number
;
;2016.06.13 -- D J Evans

function breathing, posi, $
  mpp=mpp, $
  lags = lags
  ;motion = motion
  ;pair = pair

  frame = 4    ;frame memory for track
  p = 10       ;maxdisp for track
  cull = 3     ;goodenough for track
  ;mpp=0.0685   ;microns per pixel
  thresh = 35  ;threshold for cutting out bad blinks

  if n_elements(lags) ne 1 then lags=1
  if (n_elements(posi[0,*]) le 2.*thresh) then return, [0,0,-1] 

  trackposi = track(posi, p, dim=2, mem=frame, good=cull,/quiet)  
  ;tracking data for particles in the ith blink

  label = reform(trackposi[3,*])
  
  ;xy coordinates of each particle, where
  ; xy[a, b, c] = b-position of particle c at time a
  xy = []
  for i=0,2 do begin
    xyi = transpose(reform(trackposi[0:1, where(label eq i)])*mpp)
    xy = [[[xy]],[[xyi]]]
  endfor

  ;calculate legs of triangle 012 formed by particles 0,1,2
  ;in each frame
  pairs = combigen(3,2)
  legs = []
  for i=0,2 do begin
    pairi = pairs[i,*]
    legi = reform(sqrt( (xy[*,0,pairi[0]] - xy[*,0,pairi[1]])^2. $
      + (xy[*,1,pairi[0]] - xy[*,1,pairi[1]])^2.))
    legs = [[legs],[legi]]
  endfor
  
  legs = transpose(legs)
  a = legs[0,*]
  b = legs[1,*]
  c = legs[2,*]
  
  ;radius of circle that circumscribes the triangle formed by [xy0,xy1,xy2]
  ;this is called the radial array
  ri = sqrt(3.)*reform(a*b*c / sqrt( (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c) ))
  ;ri = total(legs, 1)/3.
  
  nri = n_elements(ri)
  if (nri le thresh) then return, [0,0,-1]

  ri = ri[5:*:lags] 
  ;----------BREATHING MODE CALCULATIONS----------------------
  rdif = ri[1:*] - ri           ;radial difference array
  r_i = (ri[1:*] + ri)/2.       ;midpoint radial array
  n = n_elements(r_i)
  t = lags*findgen(n)

  master = [[r_i], [rdif], [t]]
;stop

return, transpose(master)

end