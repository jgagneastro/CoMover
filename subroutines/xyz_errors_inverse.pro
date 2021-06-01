Pro xyz_errors_inverse, X, Y, Z, EX, EY, EZ, ra, dec, dist, era, edec, edist, FAST=fast, GL=gl, GB=gb
  ;+
  ; NAME:
  ;       XYZ_ERRORS_INVERSE
  ;
  ; PURPOSE:
  ;       Computes astrometry (ra, dec) and distance (in parsec) and their errors from heliocentric galactic coordinates (X,Y,Z) and their errors (EX,EY,EZ)
  ;       X points towards the galactic center, Y towards the local direction of rotation in the plane of the galaxy, and Z towards the North Galactic Pole
  ;       (so that the XYZ system is right-handed).
  ;
  ; CALLING SEQUENCE:
  ;       XYZ_ERRORS, X, Y, Z, EX, EY, EZ, ra, dec, dist[, era, edec, edist, /FAST] 
  ;
  ; OUTPUTS:
  ;       RA = Scalar or 1D array indicating right ascension (in decimal degrees).
  ;       DEC = Scalar or 1D array indicating declination (in decimal degrees).
  ;       DIST = Scalar or 1D array indicating distance (in parsec).
  ;       ERA = Scalar or 1D array indicating the error on the right ascension (in decimal degrees).
  ;       EDEC = Scalar or 1D array indicating the error on the declination (in decimal degrees).
  ;       EDIST = Scalar or 1D array indicating the error on distance (in parsec).
  ;
  ; OPTIONAL INPUT KEYWORD:
  ;       /FAST - If this keyword is set, then no error propagation is computed.
  ;
  ; INPUTS:
  ;       X, Y, Z = Heliocentric galactic position. See "Purpose" for more details".
  ;       EX, EY, EZ = Errors on heliocentric galactic position. See "Purpose" for more details".
  ;
  ; RESTRICTIONS:
  ;       (1) If you give 1D arrays as inputs, each of them must have the same size.
  ;
  ; PROCEDURES USED:
  ;       GLACTC_ERRORS
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, April 7th, 2013.
  ;-
  
  forward_function glactc_errors, glactc
  
  ;On calcule la distance
  dist = sqrt(X^2 + Y^2 + Z^2)
  
  ;On Calcule les coordonnees galactiques
  gl = atan(Y,X)/!dtor
  gb = 90d0 - atan(sqrt(X^2+Y^2),Z)/!dtor
  
  ;On transforme les coordonnees galactiques en (RA,DEC)
  ;glactc, ra, dec, 2000., gl, gb, 2, /DEGREE
  glactc_errors_inverse, gl, gb, !NULL, !NULL, ra, dec, /FAST
  
End