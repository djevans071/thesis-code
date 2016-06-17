;
;--CLEANUP--
; Removes bad data from LABELING. Bad data meaning frames where FEATURE
; returns (x,y) = (0,0) for at least one particle where there should be
; NUM particles in total.
;
; SYNTAX:
;   newdat = cleanup(data, num)
; 
; ARGUMENT:
;   data: output from LABELING, where
;   data[0,*]: x-position
;   data[1,*]: y-position
;   data[2,*]: frame number (time)
;   data[3,*]: blink number
; 
; RETURN VALUE:
;   newdat: cleaned up dataset with same columns as input
;   
; PROCEDURE:
;   Within each blink with bad data, identify frames with bad data. 
;   Find maximum frame number where bad data occurs and keep 
;   all data after it. Return cleaned set. 
;
;2014.05.20 - D J Evans
;2016.05.11 - Added NUM keyword to take into account the number of
;             particles per frame


function cleanup, data,$
                  num        
;cleanup of bad frames from LABELING
;;(where feature returns x=0, y=0)

newdat = fltarr(3,1)

bl_col = data[3,*]            ;blink column
n = max(bl_col)               ;number of blinks

i=1
while (i le n) do begin $      ;loop over blinks
  ;print, i
  bl_dati = data[*,where(bl_col eq i)]  ;ith blink data
  xi = reform(bl_dati[0,*])  ;x-values for the ith blink
  ti = bl_dati[2,*]        ;t-values for the ith blink
  
  ti = ti - min(ti)
  bl_dati[2,*] = ti
  
  mt = max(ti)
  ;print, mt
    
  tt = transpose([[findgen(mt+1)], [fltarr(mt+1)]])
  ;zero = [0,0]
  ;stop
  
  j=1                        ;loop over frames within ith blink
  while (j le mt) do begin $
    ;print, j
    tt[1,j] = n_elements(where(ti eq j))
    j++
  endwhile
  
  t0 = tt[0,where(tt[1,*] ne num)]      ;t-values where x=0
  m = max(t0)                       ;max t-value where x=0
  ;if (m gt 70.) then m=0.
  if (m eq mt) then begin $
    bl_dati=bl_dati[*,0:-2]
    m = t0[-2]
    endif
  ;m = max(bl_dati[2,where(xi eq 0.)])      
  
  
    blk = [0,0,-1]
  bl_trunci = bl_dati[0:2,where(ti gt m)]
    bl_trunci = [[bl_trunci], [blk]]
  newdat = [[newdat],[bl_trunci]]
  i++
    
endwhile
newdat = newdat[*,1:*]

;-----relabeling
dat = labeling(newdat)

return, dat
end