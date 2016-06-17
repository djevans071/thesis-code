;
;--TRACK_PMMA--
; Master routine that tracks particles detected by PMMA_PRETRACK. Forms 
; continuous trajectories of the relative positions of particles 
; for each blink (with RELDISPCALC). 
; 
; SYNTAX:
;   RF = track_pmma(folder, [tau=value, motion=stringvalue])
;   
; ARGUMENT:
;   folder: a string specifying the location of image files
;      (to identify the output file of PMMA_PRETRACK)
; 
; KEYWORD PARAMETERS:
;   tau: specifies step size in time (default value of tau is 1 frame)
;   motion: specifies the type of motion to be tracked:
;           motion = 'rm' -- relative motion parallel and normal to line
;                            of separation
;           motion = 'cm' -- center-of-mass motion parallel and normal 
;                            to line of separation
;   
; RETURN VALUE:
;   RF: array of tracked particle data where
;   RF[0,*]: relative displacement (in microns)
;   RF[1,*]: relative displacement difference by tau steps parallel
;              to the line of separation (in microns)
;   RF[2,*]: relative displacement difference by tau steps normal
;              to the line of separation (in microns)
;   RF[3,*]: frame number
; 
; 2013.07.02 -- D J Evans
; 2014.09.12 -- Added MOTION keyword in order to measure all modes of motion
;              (2 types each of relative and center-of-mass motion) 
;              - D J Evans

function track_pmma,  folder, $
                      mpp = mpp, $
                      tau = tau, $
                      motion = motion, $
                      ;pair = pair, $
                      num = num

COMPILE_OPT IDL2

if n_elements(tau) ne 1 then tau=1
if not keyword_set(num) then num = 2

if num eq 2 then RF = fltarr(4,1)
if num gt 2 then RF = fltarr(3,1) 

  posb = read_gdf(folder + '/pretrack.gdf') 
  posb = labeling(posb)
  posb = cleanup(posb,num) 
  
  ;window,4 & plot_hist, posb[0,*] mod 1.       ;check pixel bias
  
  label = reform(posb[3,*]) 
  m = max(label)                    ;number of total data sets from labeling()
  
  i=1
  while (i le m) do begin $            ;loop over blinks
    posi = posb[0:2, where(label eq i)] 
    n = n_elements(posi[0,*]) 
    ;min_maxx = [min(posi[0,*]*mpp),max(posi[0,*]*mpp)]
    ;min_maxy = [min(posi[1,*]*mpp),max(posi[1,*]*mpp)]
    ;x_center = mean(min_maxx)
    ;y_center = mean(min_maxy)
    ;print, i
        ;xy0_i = xy_track(posi, mpp=mpp, lags = tau, part = 0)
        ;xy1_i = xy_track(posi, mpp=mpp, lags = tau, part = 1)
        ;xy2_i = xy_track(posi, mpp=mpp, lags = tau, part = 2)
        
        ;plot, xy0_i[0,*]-x_center, xy0_i[2,*]-y_center, xrange=[min_maxx[0]-x_center,min_maxx[1]-x_center], $
        ;  yrange=[min_maxy[0]-y_center,min_maxy[1]-y_center], /iso, /yno
        ;wait, 1
        ;oplot, xy1_i[0,*]-x_center, xy1_i[2,*]-y_center
        ;wait, 1
        ;oplot, xy2_i[0,*]-x_center, xy2_i[2,*]-y_center
        
        
        if (num eq 2) then RF_i = reldispcalc(posi,mpp=mpp, lags=tau, motion=motion)
        if (num gt 2) then RF_i = breathing(posi, mpp=mpp, lags=tau)
      ;finding relative displacements
      ;stop
      
      if (RF_i[-1,0] lt 0) then i++ else begin $ ;if (RF_i[3,0] ge 0) then begin $
         RF = [[RF],[RF_i]]
         i++
         endelse

  endwhile 
  
RF = RF[*,1:*] 

write_gdf, RF, folder + '/RF' + strtrim(tau,2) + motion + '.gdf' 
print, 'trajectories recorded for tau = ' + strtrim(tau,2)
 
return, RF
end