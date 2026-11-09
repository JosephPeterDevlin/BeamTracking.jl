faddeeva(z::TPS) = wf(z)


# Faddeeva implementation from Yi-Kai Kan (BNL)
@inline faddeeva(z) = faddeeva(z, Val(48))


function __faddeeva_coeff(n, type::Type=Float64)
  N = n
  M = 2*N
  L = sqrt(N/sqrt(type(2)))
  theta_arr = [k*type(pi)/M for k in (-M+1):(M-1)]
  t_arr = L*tan.(theta_arr/2)
  f_arr = [(L^2 + t^2)*exp(-t^2) for t in t_arr]
  f_arr = vcat(zero(type), f_arr)
  a_arr = real.(fft(fftshift(f_arr)))/(2*M)
  return a_arr[2:(N+1)], L
end


@generated function faddeeva(z::Complex{T}, ::Val{N}) where {T,N}
  a_arr, L = __faddeeva_coeff(N, T)
  a_arr = reverse(a_arr)
  ex = :($(a_arr[1]))
  for i in 2:N
    ex = :(muladd($ex, Z, $(a_arr[i])))
  end
  return quote
    c1 = imag(z) >= 0 ? 0 : 1
    c2 = 1 - 2*c1
    z = c2*z
    lmiz_inv = 1/($L - im*z)
    Z = ($L + im*z)*lmiz_inv
    return c1*2*exp(-z^2) + c2*(2*$ex*lmiz_inv + $(1/sqrt(T(pi))))*lmiz_inv
  end
end


function bb_derivatives(x, y, sigma_x, sigma_y)
  if sigma_x ≈ sigma_y
    return bb_derivatives_round(x, y, sigma_x)
  elseif sigma_x > sigma_y
    return bb_derivatives_elliptical(x, y, sigma_x, sigma_y)
  else
    fy, fx, fyy, fxx = bb_derivatives_elliptical(y, x, sigma_y, sigma_x)
    return fx, fy, fxx, fyy
  end
end


expm1(x::TPS) = exp(x) - 1


# (exp(-x) - 1)/x
@generated function expm1_neg_over(x::T) where {T}
  TS = ForwardDiff.valtype(T)
  if TS == Float16 || TS == Float32
    threshold = (120*eps(TS))^(1//4)
  else
    threshold = 0.0004
  end
  return :(vifelse(abs(x) > $threshold, expm1(-x)/x, -1 + x/2))
end
# For TPSA, must sum series (TBD)


function bb_derivatives_round(x, y, sigma)
  x2 = x*x
  y2 = y*y
  r2 = x2 + y2
  coeff = 2*expm1_neg_over(r2/(2*sigma*sigma))
  return fx, fy, fxx, fyy
end


function bb_derivatives_elliptical(x, y, sigma_x, sigma_y)
  sigma_x2 = sigma_x*sigma_x
  sigma_y2 = sigma_y*sigma_y
  sigma_diff = sigma_x2 - sigma_y2
  sq = sqrt(2*sigma_diff)
  dens = exp(-x*x/(2*sigma_x2) - y*y/(2*sigma_y2))

  arg1 = (x + 1im*y)/sq
  arg2 = (x*sigma_y/sigma_x + 1im*y*sigma_x/sigma_y)/sq
  w1 = faddeeva(arg1)
  w2 = faddeeva(arg2)

  g = -sqrt(2*pi/sigma_diff)*(w1 - dens*w2)
  fx = imag(g) 
  fy = real(g)

  helper = x*fx + y*fy 
  fxx = -(helper + 2*(1 - sigma_y/sigma_x*dens))/sigma_diff
  fyy =  (helper + 2*(1 - sigma_x/sigma_y*dens))/sigma_diff
  return fx, fy, fxx, fyy
end


bb_derivatives_elliptical(x::SIMD.Vec, y::SIMD.Vec, sigma_x, sigma_y) = vmap((x, y) -> bb_derivatives_elliptical(x, y, sigma_x, sigma_y), x, y)