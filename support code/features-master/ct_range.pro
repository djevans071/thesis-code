;+
; NAME:
;    ct_range
;
; PURPOSE:
;    Identify range of pixels in an image that are transformed to a particular
;    feature by circletransform
;
; CATEGORY:
;    Image analysis, feature detection
;
; CALLING SEQUENCE:
;    b = ct_range(a, p)
;
; INPUTS:
;    a: [nx,ny] gray-scale image data
;    p: [2, npts] coordinates of target pixels
;
; KEYWORD PARAMETERS:
;    noise: estimate for additive pixel noise.
;        Default: noise estimated by MAD().
;
;    deinterlace: if set to an odd number, then only perform
;        transform on odd field of an interlaced image.
;        If set to an even number, transform even field.
;        Default: Not set or set to 0: transform entire frame.
;
; OUTPUTS:
;    rad: [npts] array of range estimates for each input point.
;
; PROCEDURE:
;    Compute the gradient of the image.  The local gradient at each
;    pixel defines a line along which the center of a circle may
;    lie.  Pixel "hits" if the line passes within range of the
;    specified feature.
;
; REFERENCE:
; F. C. Cheong, B. Sun, R. Dreyfus, J. Amato-Grill, K. Xiao, L. Dixon
; & D. G. Grier,
; Flow visualization and flow cytometry with holographic video
; microscopy, Optics Express 17, 13071-13079 (2009)
;
; MODIFICATION HISTORY:
; 01/16/2013 Written by David G. Grier, New York University.
; 01/22/2013 DGG Use CLUSTER() to restrict hits to
;    nearest-neighborhood.  Removed RANGE keyword.
; 01/24/2013 DGG Correct test for deinterlace = 0.
; 01/25/2013 DGG Fix hit test for single-target case.
; 02/09/2013 DGG Use SAVGOL2D() to compute derivative kernel.
;    Displace by half a pixel to center.
; 02/15/2013 DGG Weight radius calculation by magnitude of gradient.
;
; Copyright (c) 2013 David G. Grier
;-

function ct_range, a_, p, $
                   noise = noise, $
                   deinterlace = deinterlace

COMPILE_OPT IDL2

umsg = 'USAGE: rad = ct_range(a, p)'

if n_params() ne 2 then begin
   message, umsg, /inf
   return, -1
endif

if ~isa(a_, /number, /array) then begin
   message, umsg, /inf
   return, -1
endif
if size(a_, /n_dimensions) ne 2 then begin
   message, umsg, /inf
   message, 'A must be a two-dimensional numeric array', /inf
   return, -1
endif
sz = size(a_, /dimensions)
nx = sz[0]
ny = sz[1]

if ~isa(p, /number, /array) then begin
   message, umsg, /inf
   message, 'P must be an array of target coordinates', /inf
   return, -1
endif
sz = size(p)
if sz[0] gt 2 or sz[1] lt 2 then begin
   message, umsg, /inf
   message, 'P should be a [2,npts] array of coordinates', /inf
   return, -1
endif
npts = (sz[0] eq 2) ? sz[2] : 1

dodeinterlace = isa(deinterlace, /scalar, /number) ? deinterlace gt 0 : 0
if dodeinterlace then begin
   n0 = deinterlace mod 2
   a = float(a_[*, n0:*:2])
endif else $
   a = float(a_)

if ~isa(noise, /scalar, /number) then $
   noise = mad(a)

rad = intarr(npts)   ; the answer

dx = savgol2d(7, 3, dx = 1)
dadx = convol(a, dx, /edge_truncate)
dady = convol(a, transpose(dx), /edge_truncate)
if dodeinterlace then dady /= 2.
grada = sqrt(dadx^2 + dady^2)           ; magnitude of the gradient
dgrada = noise * sqrt(2. * total(dx^2)) ; error in gradient estimate due to noise

w = where(grada gt 2.*dgrada, ngood)    ; select points with small angular uncertainty
if ngood le 0 then $
   return, rad - 1

grada = grada[w]
dgrada = dgrada[w]/grada
costheta = dadx[w]/grada
sintheta = dady[w]/grada

xy = array_indices(a, w)
if dodeinterlace then xy[1,*] = 2.*xy[1,*] + n0
xy += 0.5                       ; to center over pixels

id = (npts gt 1) ? cluster(xy, p[0:1, *]) : intarr(ngood)

for n = 0, npts-1 do begin
   qx = xy[0,*] - p[0, n]
   qy = xy[1,*] - p[1, n]
   rsq = qx^2 + qy^2
   delta = abs(qx * sintheta - qy * costheta) ; distance from voting line to target
   ddelta = dgrada * abs((qx * costheta + qy * sintheta)) < 2.5 ; uncertainty
   ww = where((delta le ddelta) and (id eq n), nhits) ; must be closer than uncertainty
   rad[n] = (nhits gt 0) ? sqrt(2.*total(grada[ww]*rsq[ww])/total(grada[ww])) : -1
endfor

return, rad
end
