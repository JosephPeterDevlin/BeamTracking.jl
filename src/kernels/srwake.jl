@makekernel fastgtpsa=true function srwake_long!(i, coords::Coords, tilde_m, tmin, dt, voltage)
  v = coords.v
  alive = coords.state[i] == STATE_ALIVE

  rel_p = 1 + v[i,PZI]
  beta = rel_p / sqrt(rel_p*rel_p + tilde_m*tilde_m)
  t = -v[i,ZI]/(beta*C_LIGHT)
  k = floor((t-tmin)/dt)
  tk = tmin + k*dt
  k = Int(k)
  w = (t-tk)/dt
  
  new_pz = v[i,PZI] + beta*((1-w)*voltage[k] + w*voltage[k+1])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end