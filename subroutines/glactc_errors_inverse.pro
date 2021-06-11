Pro glactc_errors_inverse, gl, gb, egl, egb, ra, dec, era, edec, FAST=fast
  ;+
  ; NAME:
  ;       GLACTC_ERRORS
  ;
  ; PURPOSE:
  ;       This procedure is similar to glactc.pro (j=1 direction) from the astrolib, but it takes measurement errors into account. Also, it does not need to
  ;       precess coordinates (which significantly fastens cases with many objects).
  ;
  ; CALLING SEQUENCE:
  ;       GLACTC_ERRORS, ra, dec, era, edec, gl, gb[, egl, egb, /FAST]
  ;
  ; INPUTS:
  ;       RA = Scalar or 1D array indicating right ascension (in decimal degrees).
  ;       DEC = Scalar or 1D array indicating declination (in decimal degrees).
  ;       ERA = Scalar or 1D array indicating the error on the right ascension (in decimal degrees).
  ;       EDEC = Scalar or 1D array indicating the error on the declination (in decimal degrees).
  ;
  ; OPTIONAL INPUT KEYWORD:
  ;       /FAST - If this keyword is set, then no error propagation is computed.
  ;
  ; OUTPUTS:
  ;       GL, GB = Galactic longitude and latitude, respectively (in degrees).
  ;       EGL, EGB = Errors on Galactic longitude and latitude, respectively (in degrees).

  ; NOTES :
  ;       This routine is based on equations from glactc.pro (astrolib), except they were
  ;       modified to avoid the need of precession, and analytical propagation of errors
  ;       were carried and added by the authors.
  ;
  ; RESTRICTIONS:
  ;       (1) If you give 1D arrays as inputs, each of them must have the same size.
  ;
  ; PROCEDURES USED:
  ;       NONE
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October 10th, 2012.
  ;-

  ;

  radeg = 1.8d2/!dpi
  
;  ;Si on ne veut pas calculer les erreurs
;  if keyword_set(fast) then begin
;    gl = double(!values.f_nan)
;    gb = gl
;  endif

  ;Positions J2000.0 du pôle Nord Galactique (b=90 degres). Carrol et Ostlie.
  rapol = 192.8595d0 ;total([12d0,51d0,26.28d0]/[1d0,60d0,3600d0],/double)*15d0
  decpol = 27.12825d0 ;total([27d0,7d0,41.7d0]/[1d0,60d0,3600d0],/double)
  ;Latitude galactique du pôle Nord céleste (delta=90 degres J2000.0). Carrol et Ostlie.
  ;! Erreur dans le Carrol et Ostlie ?? C'est écrit 123,55,55.2 dans le livre (p900), ce qui est incoherent avec glactc.pro
  lnord = 122.932d0 ;total([122d0,55d0,55.2d0]/[1d0,60d0,3600d0],/double)
  
  ;Raccourcis
  sdp = sin(decpol/radeg)
  cdp = cos(decpol/radeg)
  
  sgb = sin(gb/radeg)
  cgb = cos(gb/radeg)
  ;dlon --> lnord
  sdec = sgb*sdp + cgb*cdp*cos((lnord-gl)/radeg)
  dec = radeg * asin(sdec)
  cdec = sqrt(1.0d0-sdec^2)
  sinf = cgb * sin((lnord-gl)/radeg) / cdec
  cosf = (sgb-sdp*sdec) / (cdp*cdec)
  ra = rapol + radeg*atan(sinf,cosf)
  gt36 = where(ra gt 360.0, Ngt36)
  if Ngt36 ge 1 then ra[gt36] = ra[gt36] - 360.0d0
  
End