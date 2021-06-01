Pro common_pm_banyan, input_data, SIZE_XYZ_MODEL=size_xyz_model, SIZE_UVW_MODEL=size_uvw_model, BINARY_FRACTION=binary_fraction, MASS_HOST=mass_host, out
  ;This takes a pair of stars and calculates the probability that they are comoving with a BANYAN-style comparison of their XYZ and UVWs
  ;The mass of the host allows to estimate a correction to inflate the star's UVW model to account for possible orbital motion

; SIZE_XYZ_MODEL controls the characteristic width of the spatial XYZ part of the host-size model in parsecs (by default 0.1 pc).
; SIZE_UVW_MODEL controls the characteristic width of the kinematic UVW part of the host-size model in km/s (by default 1 km/s), to which the contribution to possible orbital motion will be added in quadrature.
; MASS_HOST is used to input a mass estimate for the host star in Solar mass units (by default 1.0 solar mass). This will be used to estimate the effect of the orbital motion of the companion on how discrepant the UVW motion is allowed to be. 
; OUT is an IDL structure that contains the outputs of this code

  ;Input data must be an IDL structure with the following tags:
  ;HOST_RA: right ascension of the host star (degrees)
  ;HOST_DEC: declination of the host star (degrees)
  ;HOST_PMRA: proper motion in right ascension of the host star (includes the cos(dec) term; units of milli arc seconds per year)
  ;HOST_PMDEC: proper motion in declination of the host star (units of milli arc seconds per year)
  ;HOST_EPMRA: measurement error of the proper motion in right ascension of the host star (units of milli arc seconds per year)
  ;HOST_EPMDEC: measurement error of the proper motion in declination of the host star (units of milli arc seconds per year)
  ;HOST_RV: heliocentric radial velocity of the host star (units of km/s; facultative keyword)
  ;HOST_ERV: measurement error of the heliocentric radial velocity of the host star (units of km/s; facultative keyword) 
  ;HOST_PLX: trigonometric parallax of the host star (units of milli arcseconds)
  ;HOST_EPLX: measurement error of the trigonometric parallax of the host star (units of milli arcseconds)
  
  ;COMPANION_RA: right ascension of the potential companion (degrees)
  ;COMPANION_DEC: declination of the potential companion (degrees)
  ;COMPANION_PMRA: proper motion in right ascension of the potential companion (includes the cos(dec) term; units of milli arc seconds per year)
  ;COMPANION_PMDEC: proper motion in declination of the potential companion (units of milli arc seconds per year)
  ;COMPANION_EPMRA: measurement error of the proper motion in right ascension of the potential companion (units of milli arc seconds per year)
  ;COMPANION_EPMDEC: measurement error of the proper motion in declination of the potential companion (units of milli arc seconds per year)
  ;COMPANION_RV: heliocentric radial velocity of the potential companion (units of km/s; facultative keyword)
  ;COMPANION_ERV: measurement error of the heliocentric radial velocity of the potential companion (units of km/s; facultative keyword) 
  ;COMPANION_PLX: trigonometric parallax of the potential companion (units of milli arcseconds; facultative keyword)
  ;COMPANION_EPLX: measurement error of the trigonometric parallax of the potential companion (units of milli arcseconds; facultative keyword)
  
  ;Here is an example input for this code (it should return a 100% co-moving probability
;  
;  nan = !values.d_nan
;  snan = 'NaN'
;  input_data = {host_name:snan, host_ra:nan, host_dec:nan, host_pmra:nan, host_pmdec:nan, host_epmra:nan, host_epmdec:nan, host_rv:nan, host_erv:nan, $
;    host_plx:nan, host_eplx:nan, companion_name:snan, companion_ra:nan, companion_dec:nan, $
;    companion_pmra:nan, companion_epmra:nan, companion_pmdec:nan, companion_epmdec:nan, $
;    companion_rv:nan, companion_erv:nan, companion_plx:nan, companion_eplx:nan}
;
;  input_data.HOST_RA = 190.88751d0
;  input_data.HOST_DEC = 60.014355d0
;  input_data.HOST_PMRA = -126.40d0
;  input_data.HOST_PMDEC = -64.14d0
;  input_data.HOST_EPMRA = 0.01d0
;  input_data.HOST_EPMDEC = 0.02d0
;  input_data.HOST_RV = -10.2
;  input_data.HOST_ERV = 0.2
;  input_data.HOST_PLX = 22.24d0
;  input_data.HOST_EPLX = 0.01d0
;  
;  input_data.COMPANION_RA = 190.8838648d0
;  input_data.COMPANION_DEC = 60.0239567d0
;  input_data.COMPANION_PMRA = -133d0
;  input_data.COMPANION_EPMRA = 8d0
;  input_data.COMPANION_PMDEC = -55d0
;  input_data.COMPANION_EPMDEC = 8d0
;  input_data.COMPANION_PLX = 1d3/44d0
;  input_data.COMPANION_EPLX = 1d3/44d0^2*4d0

;  common_pm_banyan, input_data, out
;  help, out

; Those should be the outputs of the sample code above:
;  FIELD P = 0.0%
;  Comoving probability = 100.0%
;  LN_P comover = 0.00
;  LN_P field = -9.91
;  XYZ Separation = 1.0 pc
;  UVW Separation = 3.0 km/s
;  % Program caused arithmetic error: Floating underflow
;  ** Structure <12051a08>, 10 tags, length=696, data length=696, refs=1:
;  NAME            STRING    'NaN'
;  ALL             STRUCT    -> <Anonymous> Array[1]
;  METRICS         STRUCT    -> <Anonymous> Array[1]
;  STAR_A          STRUCT    -> <Anonymous> Array[1]
;  FIELD           STRUCT    -> <Anonymous> Array[1]
;  BESTYA_STR      STRUCT    -> <Anonymous> Array[1]
;  YA_PROB         DOUBLE          0.99995056
;  LIST_PROB_YAS   STRING    'NaN'
;  BEST_HYP        STRING    'NaN'
;  BEST_YA         STRING    'NaN'
  
  ;List all the subroutines
  forward_function sample_banyan_models,xyz_errors_inverse,uvw_errors_inverse,plothist,supersmooth,xyz_errors,uvw_errors,bsigma_par,pm_min_error,alog_sum,banyan_sigma,nan_str

  ;This is the size of the UVW bubble to be tested around each "host star" model hypothesis. Set this to the very minimum kinematic accuracy you expect to obtain
  ;  accounting systematics. It may be unwise to set this smaller than 0.5 km/s unless gravitational reddening is carefully accounted for especially when the RV of the
  ;  host star is known.
  if ~keyword_set(size_uvw_model) then $
    size_uvw_model = 1.0;km/s due to intrinsic errors

  ;This related to the maximum allowed separation. A true distribution following https://arxiv.org/pdf/1912.04150.pdf would be preferred, assuming correct normalization with "binary_fraction"
  ;The recent paper by Nelson et al suggests that pairs could exist all the way to 10^7 AU = 40 pc but are those just moving groups or something https://arxiv.org/pdf/2104.12883.pdf
  if ~keyword_set(size_xyz_model) then $
    size_xyz_model = 0.1;pc

  ;This is the fraction of stars expected to have a companion (assuming the right respective spectral types)
  if ~keyword_set(binary_fraction) then $
    binary_fraction = 1.0

  ;Issue an error when no input data are given
  if ~keyword_set(input_data) then $
    message, ' You must input some observables - please see the header of this code for more information'
  ;Issue an error when no host star parallax is available (for now)
  if ~finite(input_data.HOST_PLX) then $
    message, ' This code currently requires a parallax measurement for the host star'
  if ~finite(input_data.COMPANION_PLX) then $
    message, ' If you do not use even a photometric estimate for the companion parallax, expect very small co-mover probabilities because there are a lot of background stars that may have motions consistent with the host star when you completely ignore photometry', /continue

  ;Assume 1 Msun if mass was not given (which is better than no assumption at all)
  if ~keyword_set(mass_host) then begin
    mass_host = 1.0
    message, ' Assuming the host is 1.0 Msun for the orbital velocity component of the host-star UVW model width', /continue
  endif
  
  ;Determine the effect of the host mass and projected orbital separation on the potential orbital velocity assuming a circular orbit with Kepler's Third Law 
  if keyword_set(mass_host) then begin
    ;Determine projected separation and use that to estimate the worst-case-scenario (smallest) physical separation between host and companion
    separation_as = sqrt((input_data.HOST_RA-input_data.COMPANION_RA)^2*cos(input_data.HOST_DEC*!dtor)^2+(input_data.HOST_DEC-input_data.COMPANION_DEC)^2)*3600
    host_dist_pc = 1d3/input_data.HOST_PLX
    separation_au = separation_as*host_dist_pc
    orbital_motion_kms = sqrt(!const.G*(mass_host*!const.M_SUN)/(separation_au*!const.AU))*1d-3
    size_uvw_model = sqrt(size_uvw_model^2 + orbital_motion_kms^2)
  endif
  
  ;Number of Monte Carlo steps to account for host-star measurement errors
  nMonte_host_errors = 1d5

  ;No-RV case for the host: marginalize over all plausible field RVs at the host star's sky position
  if ~finite(input_data.HOST_RV) then begin
  
    ;Create a Monte Carlo realization of field stars using the BANYAN Sigma models to determine the RV distribution of field stars in the direction of the host star
    field_sample = sample_banyan_models('FIELD', nMonte=1d7)
    xyz_errors_inverse, field_sample[*,0], field_sample[*,1], field_sample[*,2], 0., 0., 0., ra, dec, dist, /fast
    radius_starA = 10d0;degree
    dA = sqrt((input_data.HOST_RA-ra)^2*cos(input_data.HOST_DEC*!dtor)^2+(input_data.HOST_DEC-dec)^2)

    ;Take the 1d3 closest stars that are within 50 degrees on the sky
    gg_field0 = where(dA lt 50., ngg0)
    gg_field = gg_field0[(sort(dA[gg_field0]))[0:(1d3<(ngg0-1L))]]
    
    uvw_errors_inverse, field_sample[gg_field,3], field_sample[gg_field,4], field_sample[gg_field,5], ra[gg_field], dec[gg_field], $
      0, 0, 0, 0, 0, vrad, dpmra, dpmdec, dist=dist[gg_field]
    rv_bin = 1.0;km/s
    plothist, vrad, bin=rv_bin, xhistrv, yhistrv, /noplot
    yhistrv = double(yhistrv)
    
    ;Smooth the RV distribution
    yhistrv = supersmooth(yhistrv,8,/trim)
    bad = where(~finite(yhistrv), nbad)
    if nbad ne 0L then yhistrv[bad] = 0.0
  endif

  ;With-RV case for the host: the RV prior is set to the RV measurement
  if finite(input_data.HOST_RV) then begin
    rv_bin = 0.1;km/s
    rv_noise = randomn(seed, nMonte_host_errors)*input_data.HOST_ERV + input_data.HOST_RV
    plothist, rv_noise, bin=rv_bin, xhistrv, yhistrv, /noplot
    yhistrv = double(yhistrv)
    
  endif

  ;Renormalize the field star density of the Besançon model to that of Kirkpatrick et al. (2012)
  norm_field = 314159.27d0/8.8d-6
  
  ;Create a series of spatial-kinematic models for the host star
  nhyp = n_elements(xhistrv)
  plx_noise = randomn(seed, nMonte_host_errors)*input_data.HOST_EPLX + input_data.HOST_PLX
  pmra_noise = randomn(seed, nMonte_host_errors)*input_data.HOST_EPMRA + input_data.HOST_PMRA
  pmdec_noise = randomn(seed, nMonte_host_errors)*input_data.HOST_EPMDEC + input_data.HOST_PMDEC
  ra_noise = replicate(input_data.HOST_RA, nMonte_host_errors)
  dec_noise = replicate(input_data.HOST_DEC, nMonte_host_errors)
  rv_noise = (make_array(nMonte_host_errors,value=1d0,/double)#xhistrv) + ((randomn(seed, nMonte_host_errors)*rv_bin)#make_array(nhyp,value=1d0,/double))

  xyz_errors, ra_noise, dec_noise, 1d3/plx_noise, 0d0, 0d0, 0d0, X, Y, Z, /FAST
  uvw_errors, (ra_noise#make_array(nhyp,value=1d0,/double)), (dec_noise#make_array(nhyp,value=1d0,/double)), $
    (pmra_noise#make_array(nhyp,value=1d0,/double)), (pmdec_noise#make_array(nhyp,value=1d0,/double)), $
    rv_noise, ((1d3/plx_noise)#make_array(nhyp,value=1d0,/double)), 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, $
    U, V, W, /FAST

  par_obj = bsigma_par('HYA')
  par_field = bsigma_par('FIELD')
  cov_obj = diag_matrix([replicate(size_xyz_model^2,3),replicate(size_uvw_model^2,3)])
  prec_obj = abs(la_invert(cov_obj))
  par_obj.NAME = 'STAR_A'
  par_obj.ln_prior = !values.d_nan
  par_obj.ln_alpha_k = !values.d_nan
  par_obj = replicate(par_obj,nhyp)

  ;Fit a UVW ellipse to each RV case
  covar0 = diag_matrix([replicate(size_xyz_model^2,3),replicate(size_uvw_model^2,3)])
  for i=0L, nhyp-1L do begin
    data = transpose([[X],[Y],[Z],[U[*,i]],[V[*,i]],[W[*,i]]])

    model_center = mean(data,dim=2)
    model_covar = correlate(data,/covariance)

    model_covar += covar0

    par_obj[i].CENTER_VEC = model_center
    par_obj[i].COVARIANCE_MATRIX = model_covar
    par_obj[i].COVARIANCE_DETERM = la_determ(model_covar)
    prec = la_invert(model_covar)
    par_obj[i].PRECISION_MATRIX = prec
    par_obj[i].PRECISION_DETERM = la_determ(prec)
  endfor

  par_obj.COEFFICIENT = double(yhistrv)/total(double(yhistrv))
  par_obj.ln_nobj = alog(1d0)

  par = [par_obj,par_field]
  ln_priors = {STAR_A:alog(binary_fraction),FIELD:alog(norm_field)}
  
  ;Create a monte carlo for the proper motion of the companion
  nMonte_pm = 1d4
  pm_min_error = 2.;mas/yr
  plx_min_error = 0.1;mas
  pmra_noise = randomn(seed, nMonte_pm)*input_data.COMPANION_EPMRA + input_data.COMPANION_PMRA
  pmdec_noise = randomn(seed, nMonte_pm)*input_data.COMPANION_EPMDEC + input_data.COMPANION_PMDEC
  ra_noise = replicate(input_data.COMPANION_RA, nMonte_pm)
  dec_noise = replicate(input_data.COMPANION_DEC, nMonte_pm)
  epmra_noise = replicate(pm_min_error,nMonte_pm)
  epmdec_noise = replicate(pm_min_error,nMonte_pm)
  plx_noise = replicate(input_data.COMPANION_PLX, nMonte_pm)
  eplx_noise = replicate(input_data.COMPANION_EPLX, nMonte_pm)
  out_mc_lnp = banyan_sigma(CUSTOM_MODELS=par,ra=ra_noise,dec=dec_noise,pmra=pmra_noise,epmra=epmra_noise,pmdec=pmdec_noise,epmdec=epmdec_noise,plx=plx_noise,eplx=eplx_noise,/unit_priors,ln_priors=ln_priors,/lnp_only,/no_normalization)
  
  ;Marginalize first and then normalize
  out_mc_marginalized = alog_sum_2d_fixed(out_mc_lnp,dim=1)
  out_mc_norm = out_mc_marginalized - alog_sum(out_mc_marginalized)
  out_p = exp(out_mc_norm)
  out_mc = banyan_sigma(CUSTOM_MODELS=par,ra=ra_noise,dec=dec_noise,pmra=pmra_noise,epmra=epmra_noise,pmdec=pmdec_noise,epmdec=epmdec_noise,plx=plx_noise,eplx=eplx_noise,/unit_priors,ln_priors=ln_priors)
  
  out = banyan_sigma(CUSTOM_MODELS=par,ra=ra_noise[0],dec=dec_noise[0],pmra=pmra_noise[0],epmra=epmra_noise[0],pmdec=pmdec_noise[0],epmdec=epmdec_noise[0],plx=plx_noise[0],eplx=eplx_noise[0],/unit_priors,ln_priors=ln_priors)
  nan_str, out
  weights = 1d0
  out.YA_PROB = out_p[0L]
  out.FIELD.PROB = out_p[1L]
  out.ALL.STAR_A = out_p[0L]
  out.ALL.FIELD = out_p[1L]
  out.STAR_A.PROB = out_p[0L]
  out.FIELD.PROB = out_p[1L]
  out.STAR_A.RV_OPT = mean(out_mc.STAR_A.RV_OPT, /nan)
  out.STAR_A.ERV_OPT = mean(out_mc.STAR_A.ERV_OPT, /nan)
  out.STAR_A.UVW_SEP = mean(out_mc.STAR_A.UVW_SEP, /nan)
  out.STAR_A.XYZ_SEP = mean(out_mc.STAR_A.XYZ_SEP, /nan)

  print, 'FIELD P = '+rtrim(out.FIELD.PROB*1d2,1)+'%'

  print, 'Comoving probability = '+rtrim(out.STAR_A.PROB*1d2,1)+'%'
  print, 'LN_P comover = '+rtrim(alog(out.ALL.STAR_A),2)
  print, 'LN_P field = '+rtrim(alog(out.ALL.FIELD),2)

  print, 'XYZ Separation = '+rtrim(out.STAR_A.XYZ_SEP,1)+' pc'
  print, 'UVW Separation = '+rtrim(out.STAR_A.UVW_SEP,1)+' km/s'
  return

End