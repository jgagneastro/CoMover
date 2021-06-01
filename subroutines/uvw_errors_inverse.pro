Pro uvw_errors_inverse, u, v, w, ra, dec, eu, ev, ew, era, edec, vrad, dpmra, dpmdec, evrad, edpmra, edpmdec, LSR=lsr, DIST=dist
  ;+
  ;Pas termine. Transforme UVW et ra, dec vers (pmra,pmdec)_reduits,vrad
  
  ;Inverse Rotation Matrix for the coordinate change (U towards Gal center)
  AI = [[-0.054875560d0, -0.87343709, -0.48383502], $
    [0.49410943d0, -0.44482963, 0.74698224],$
    [-0.86766615d0, -0.19807637, 0.45598378]]

  ;1 A.U/yr in km/s
  k = 4.743717361d0;. from Wolfram Alpha.
  
  ;Shortcuts
  nra = n_elements(ra)
  cosd = cos(dec*!dtor)
  sind = sin(dec*!dtor)
  cosa = cos(ra*!dtor)
  sina = sin(ra*!dtor)
  T1 =  AI[0,0]*cosa*cosd + AI[1,0]*sina*cosd + AI[2,0]*sind
  T2 = -AI[0,0]*sina      + AI[1,0]*cosa
  T3 = -AI[0,0]*cosa*sind - AI[1,0]*sina*sind + AI[2,0]*cosd
  T4 =  AI[0,1]*cosa*cosd + AI[1,1]*sina*cosd + AI[2,1]*sind
  T5 = -AI[0,1]*sina      + AI[1,1]*cosa
  T6 = -AI[0,1]*cosa*sind - AI[1,1]*sina*sind + AI[2,1]*cosd
  T7 =  AI[0,2]*cosa*cosd + AI[1,2]*sina*cosd + AI[2,2]*sind
  T8 = -AI[0,2]*sina      + AI[1,2]*cosa
  T9 = -AI[0,2]*cosa*sind - AI[1,2]*sina*sind + AI[2,2]*cosd
  
  ;Compute observables
  vrad   = T1*U + T4*V + T7*W
  dpmra  = T2*U + T5*V + T8*W
  dpmdec = T3*U + T6*V + T9*W
  
  if keyword_set(dist) then begin
    ;Compute real proper motions if distances are given
    dpmra = dpmra / (k * dist) * 1d3
    dpmdec = dpmdec / (k * dist) * 1d3
  endif
  
  return
End