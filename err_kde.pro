function err_kde, x, y, dx, $
  weight = weight, $
  scale = scale
    
  COMPILE_OPT IDL2, HIDDEN
  
  nx = float(n_elements(x))       ; number of data points
  ny = n_elements(y)              ; number of samples
  
  if n_elements(weight) ne nx then $
    weight = replicate(1., nx)
    
  ; optimal smoothing parameter
  ; Silverman Eqs. (3.30) and (3.31)
  sx = stddev(x)                  ; standard deviation
  rx = iqr(x)                     ; interquartile range
  if rx lt 1e-10 then $           ; more than 3/4 data have same value
    h = 0.9 * sx / nx^0.2 $
  else $
    h = 0.9 * (sx < rx/1.34) / nx^0.2
    
  if arg_present(scale) then scale = h
  
  ; density estimate
  ; Silverman Eq. (2.15) and Table 3.1
  t = x/h
  s = y/h
  res = fltarr(ny)                ; result
  variance = fltarr(ny)           ; variance in result
  bias = fltarr(ny)
  mse = fltarr(ny)                ; asymptotic mean-squared error
  sigma = fltarr(ny)              ; statistical error of density
  
  
  ; Epanechnikov
    norm = (3./4./sqrt(5.))
    norm1 = 1./(h*nx)
    for j = 0L, ny-1L do begin
      z = 2.*(t - s[j])/5.
      w = where(abs(z) lt 1., ngood)
      if ngood gt 0L then begin
        ker = norm * (1. - z[w])
        val = weight[w]*ker
        res[j] = dx*norm1*sqrt(total(val^2.))
        
        sigma[j] = total(val^2)
        variance[j] = total((val - res[j])^2)/nx
        bias[j] = norm/5. * total(weight[w])
        mse[j] = norm * res[j]^2 / total(ker) + bias[j]^2
        
      endif
    endfor
  ;stop
  sigma = sqrt(sigma)
  
  return, res
end