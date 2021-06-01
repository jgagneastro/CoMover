Function sample_banyan_models, hyp, nMonte=nMonte
  if ~keyword_set(nMonte) then nMonte = 1d6
  par = bsigma_par(hyp)
  npar = n_elements(par)
  sample = dblarr(nMonte, 6L,npar)+!values.d_nan
  for i=0L, npar-1L do $
    sample[*,*,i] = mrandomn(seed, par[i].COVARIANCE_MATRIX, nMonte) + make_array(nMonte,value=1d0,/double)#par[i].CENTER_VEC
  
  sample_concat = reform(transpose(sample,[0,2,1]),nMonte*npar,6L)
  return, sample_concat
;  sample2 = dblarr(nMonte*npar,6L)+!values.d_nan
;  for i=0L, 5L do sample2[*,i] = reform(sample[*,i,*])
End