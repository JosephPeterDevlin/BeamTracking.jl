@makekernel function fastgtpsa=true deposit_long!(i, coords::Coords, tilde_m, tmin, tmax, dt)
  v = coords.v
  alive = coords.state[i] == STATE_ALIVE
  rel_p = 1 + v[i,PZI]
  beta = rel_p / sqrt(rel_p*rel_p + tilde_m*tilde_m)
  t = -v[i,ZI]/(beta*C_LIGHT)
  k = floor(Int, (t-tmin)/dt)
  tk = tmin + k*dt
  w = (t-tk)/dt
  coords.longitudinal_density[k]   += vifelse(alive, 1-w, 0)
  coords.longitudinal_density[k+1] += vifelse(alive, w,   0)
end