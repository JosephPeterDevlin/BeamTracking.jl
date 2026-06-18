@makekernel fastgtpsa=true function deposit_long!(i, coords::Coords, tilde_m, tmin, tmax, dt)
  v = coords.v
  alive = coords.state[i] == STATE_ALIVE
  rel_p = 1 + v[i,PZI]
  beta = rel_p / sqrt(rel_p*rel_p + tilde_m*tilde_m)
  t = -v[i,ZI]/(beta*C_LIGHT)
  k = floor((t-tmin)/dt)
  tk = tmin + k*dt
  k = Int(k)
  w = (t-tk)/dt
  weight_factor = 1
  if !isnothing(coords.weight)
    if coords.weight isa Number
      weight_factor = coords.weight
    else
      weight_factor = coords.weight[i]
    end
  end
  coords.longitudinal_density[k]   += vifelse(alive, weight_factor*(1-w), 0)
  coords.longitudinal_density[k+1] += vifelse(alive, weight_factor*w,     0)
end