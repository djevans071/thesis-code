;
;--PMMA_PRETRACK--
; Detects features in a series of image files taken from a blinking
; optical trapping setup, with each blink a subset of frames where
; one continuous motion occurs (when the traps are turned off). 
; Creates an array of raw position information used to find 
; trajectories of particles in a given blink.
; 
; SYNTAX:
;   p = pmma_pretrack(folder, fpb, [diam=value, motion=4-value])
;   
; ARGUMENTS:
;   folder: a string specifying the location of image files
;      to be processed
;   fpb: number of frames per blink
;
; KEYWORD PARAMETERS:
;   diam: specifies the diameter of the particles for FEATURE to detect
;         (in pixels)
;   window: a four element array specifying the desired image window to
;           detect features. 
;             Form: window = [x0,x1,y0,y1],
;             where x0:x1 is the x-range and y0:y1 is the y-range.
;             
; RETURN VALUE:
;   d: pretrack data containing position data for all particles detected
;     with a delimiter row to differentiate blinks
;   d[0,*]: x-position
;   d[1,*]: y-position
;   d[2,*]: frame number (or time)
;      (the function also writes a GDF file containing this data) 
;      
; PROCEDURE:
;   Search folder for all image files needed. For each image, cut background 
;   and detect features using FEATURE routine. Concatenate data table 
;   with new feature data for each image. After fpb frames, create 
;   delimiter row [-1,-1,-1] in the data variable to specify the end 
;   of a blink (LABELING will sort out and label blinks). 
;   Return masterdata and write GDF. 
;
; 2014.06.24 -- D J Evans
; 2014.10.28 -- Added DIAM and WINDOW keywords - D J Evans
; 2016.05.09 -- Added NUM keyword for handling more than two particle detection - DJE


function pmma_pretrack, folder, $   ;folder should name the folder the files are located
                        fpb, $      ;e.g. folder = '50px' corresponds to the path /50px
                        diam=diam, $;image size input for feature (diameter in pixels)
                        window=window, $ ;desired viewing region in the form of a
                        num=num   ;number of particles to be detected
                        ;4-element array, [x0,x1,y0,y1]


s=file_search('./' + folder + '/data/*/*.png', count = n) 
bg = read_gdf('./' + folder + '/bg.gdf')

;fpb = 80   ;frames per blink (now an input)

x0 =  10  &  x1 = 450
y0 = 50  &  y1 = 500
if n_elements(window) ne 4 then window=[x0,x1,y0,y1]
if n_elements(diam) ne 1 then diam=41


t = fltarr(n)

d = [-1,-1,-1]
blk = [0,0,-1]
win = window

;loop for finding features
tic
i=0.
print, ["blink no.", "frame no."]
while (i lt n) do begin $                   
  
  a = (float(read_image(s[i])/bg))[win[0]:win[1],win[2]:win[3]]   ;read images and cut background
  dd = byte(diam)
    
  b = 255-bytscl(a)

  if i mod fpb eq 0 then begin $
    d = [[d],[blk]] 
    print, [i/fpb, i] 
    endif

  f = ctfeature(a, /gradient, smoothing=2, pickn=num)
  npts = n_elements(f[0,*])
  time = i # replicate(1,npts)        ;create time column
      x = f[0,*] & y = f[1,*]
      r = reform(sqrt(x^2. + y^2.))
      ;stop
;----if ctfeature only returns one feature, count
;    the frame as zero (it'll be cut out with cleanup)      
  if npts lt num then begin $  
    f=[0,0,time[0]]
    f0 = fover2d(bytscl(a), f, /circle, radius=dd/2)
        
    d = [[d],[f]]
    i++
    endif else begin $
  
    ;----if ctfeature finds features too close together,
    ;    count the frame as zero
      ;if (r[1]-r[0]) lt dd then begin $
      ;  ;f = [0,0,time[0]] 
      ;  f0 = fover2d(bytscl(a), f, /circle, radius=dd/2)
      ;  stop
      ;  d = [[d],[f]]
      ;  i++
      ;  endif else begin $    
    
      f0 = fover2d(bytscl(a), f, /circle, radius=dd/2)  
    
      thisd = f[0:1,*]
      thisd = [thisd, time]
      
      d = [[d],[thisd]]
      i++
      ;stop 
   endelse

endwhile
toc

write_gdf, d, './' + folder + '/pretrack.gdf'    
            ;saving pretrack data
print, 'GDF files written'

return, d

end
  
  
