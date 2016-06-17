;
;--LABELING--
;
; Places a label on each data set consisting of points continuously 
; gridded in time, in a masterset with discontinuities in the time 
; data (i.e. values of -1 in data[2,*]). 
; Each new data set after a disconinuity will be given a new label. 
;
; SYNTAX:
;   truncdat = labeling(data)
;   
; ARGUMENT:
;   data: Input data must be in the form [xp, yp, time] where xp is the 
;     x-position column of the data and so on. 
;     Usually the output of PMMA_PRETRACK
; 
; RETURN VALUE:
;   truncdat: Returns the truncated data set (minus discontinuities) 
;             plus the labeled data sets for use in the tracking phase
;   truncdat[0,*]: x-position
;   truncdat[1,*]: y-position
;   truncdat[2,*]: frame number (time)
;   truncdat[3,*]: blink number
;
;2012.04.02 - D J Evans

function labeling, data

truncdat = data[*,where(data[2,*] gt 0)]     ;truncated data set (all -1 values removed)

w = where(data[2,*] gt 0)        ;index array where time is gridded continuously
dw = w - w[1:*]                  ;difference in offset index arrays
                                  ;(provides a list of trigger points)

index = where(dw lt -1)          ;set of trigger points in w
index = [0,index,n_elements(dw)]

nsets = n_elements(index) - 1        ;number of data sets

indexcol = 1

;tic
i=1
while (i le nsets) do begin
  t = w[index[i-1]+1 : index[i]]     ;time data in between trigger points
  n = n_elements(t) 
  ind = i # replicate(1,n)     ;making index column for the ith data set
  indexcol = [indexcol,reform(ind)] 
  i++
  ;stop
endwhile
;toc

truncdat = [truncdat,transpose(indexcol)] 

return, truncdat

end 

