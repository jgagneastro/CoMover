Function bsigma_par, name
  par = banyan_sigma(ra=1,dec=1,pmra=1,pmdec=1,epmra=1,epmdec=1,/return_parameters)
  if keyword_set(name) then begin
    gg = where(strtrim(par.NAME,2) eq strupcase(strtrim(name,2)), ngg)
    if ngg eq 0 then return, -1
    par = par[gg]
  endif
  return, par
End