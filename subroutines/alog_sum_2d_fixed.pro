Function alog_sum_2d_fixed, x, DIM=dim
  ;This function was written by J. Gagne to sum a 2D array X along a given axis in log space.
  ;X must already be a natural log quantity.
  ;The dimension index works like that of total.pro.

  forward_function alog_sum

  if (size(x))[0] eq 1L then return, alog_sum(x)

  ;Find a normalization factor that avoids numerical explosions
  norm_factor = max(x,dim=dim,/nan)
  norm_factor_2D = norm_factor#make_array(((size(x)))[dim],value=1d0,/double)
  if dim eq 1L then norm_factor_2d = transpose(norm_factor_2D)
  
  ;Calculate the log of the sum
  return, alog(total(exp(x-norm_factor_2D),/nan,dim))+norm_factor

End