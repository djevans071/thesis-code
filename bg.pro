;
;--BG--
; Creates an average background image from a series of images
; of the same size. A simplified form of DHMBACKGROUND
; 
; SYNTAX:
;   background = bg(folder)
; 
; ARGUMENT:
;   folder: a string specifying the location of the image files
;      to be processed.
;   
; RETURN VALUE:
;   res: the average background image
;     (the function also writes an output file in GDF format)
; 
; 2013.06.13 - D J Evans

function bg, folder

s=file_search('./' + folder + '/bg*/*.png', count = n) 

sample = (read_image(s[0]))
w = n_elements(sample[*,0])     ;width of background image
h = n_elements(sample[0,*])     ;height of background image
npts = w * h                    ;number of total pixels 

b = fltarr(npts, 256)           
ndx = lindgen(npts)             
;c = fltarr(w,h)

i=0
while i lt n do begin
  if i mod 100 eq 0 then print, i
    a = (read_image(s[i]))
    b[ndx, a]++
    i++
endwhile

b = total(b, 2, /cumulative) - n/2.
m = min(b, med, /absolute, dim = 2)
med = float(med)/npts > 1
res = reform(med,w,h)

write_gdf, res, './' + folder + '/bg.gdf'
return, res

end