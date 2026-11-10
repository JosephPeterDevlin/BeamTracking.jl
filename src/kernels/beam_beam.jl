@makekernel fastgtpsa=true function bb_step!(i, coords::Coords, s, first, last, order, kick_coeff, beta_x, beta_y, sigma_x, sigma_y, sigma_z, L)
  u = 2*s - 1 # transforming (0,1) to (-1,1)
  L = 2*L

  # spin and radiation

  if first
    first_bb_drift!(i, coords, order, sigma_z, L/2)
  else
    bb_drift!(i, coords, u, order, sigma_z, L/2)
  end

  bb_kick!(i, coords, u, order, kick_coeff, beta_x, beta_y, sigma_x, sigma_y, sigma_z, L)

  if last
    last_bb_drift!(i, coords, order, sigma_z, L/2)
  else
    bb_drift!(i, coords, u + L/2, order, sigma_z, L/2)
  end

  # spin and radiation
end


@makekernel fastgtpsa=true function first_bb_drift!(i, coords::Coords, order, sigma_z, L)
  v = coords.v

  L_eff = v[i,ZI] - 2*order*sigma_z*erfcinv(L)
  rel_p = 1 + v[i,PZI]
  pt2 = v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
  ps2 = rel_p*rel_p - pt2
  good_momenta = (ps2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps2_1 = one(ps2)
  ps = sqrt(vifelse(good_momenta, ps2, ps2_1))

  denom = rel_p + ps

  new_x = v[i,XI] + v[i,PXI]/denom*L_eff
  new_y = v[i,YI] + v[i,PYI]/denom*L_eff
  new_z = v[i,ZI] - pt2/(denom*denom)*L_eff
  new_pz = (v[i,PZI] + ps)/2
  v[i,XI]  = vifelse(alive, new_x, v[i,XI])
  v[i,YI]  = vifelse(alive, new_y, v[i,YI])
  v[i,ZI]  = vifelse(alive, new_z, v[i,ZI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end


@makekernel fastgtpsa=true function last_bb_drift!(i, coords::Coords, order, sigma_z, L)
  v = coords.v

  L_eff = v[i,ZI] + 2*order*sigma_z*erfcinv(L)
  rel_p_unsafe = 1 + 2*v[i,PZI]
  good_momenta = (rel_p_unsafe > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  rel_p_1 = one(rel_p_unsafe)
  rel_p = vifelse(good_momenta, rel_p_unsafe, rel_p_1)

  pt2 = v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
  ps_unsafe = v[i,PZI] + 1/2 - pt2/(2*rel_p)
  good_momenta = (ps_unsafe > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps_1 = one(ps_unsafe)
  ps = vifelse(good_momenta, ps_unsafe, ps_1)

  denom = 2*ps

  new_x = v[i,XI] - v[i,PXI]/denom*L_eff
  new_y = v[i,YI] - v[i,PYI]/denom*L_eff
  new_z = v[i,ZI] + pt2/(rel_p*denom)*L_eff
  new_pz = v[i,PZI] - 1/2 + pt2/(2*rel_p)
  v[i,XI]  = vifelse(alive, new_x, v[i,XI])
  v[i,YI]  = vifelse(alive, new_y, v[i,YI])
  v[i,ZI]  = vifelse(alive, new_z, v[i,ZI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end


@makekernel fastgtpsa=true function bb_drift!(i, coords::Coords, u, order, sigma_z, L)
  v = coords.v

  L_eff = order*sigma_z*(erfinv(u + L) - erfinv(u))
  pt2 = v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
  rel_p_unsafe = 1 + 2*v[i,PZI]
  good_momenta = (rel_p_unsafe > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  rel_p_1 = one(rel_p_unsafe)
  rel_p = vifelse(good_momenta, rel_p_unsafe, rel_p_1)

  new_x = v[i,XI] + 2*v[i,PXI]/rel_p*L_eff
  new_y = v[i,YI] + 2*v[i,PYI]/rel_p*L_eff
  new_z = v[i,ZI] - 2*pt2/(rel_p*rel_p)*L_eff
  v[i,XI] = vifelse(alive, new_x, v[i,XI])
  v[i,YI] = vifelse(alive, new_y, v[i,YI])
  v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
end


@makekernel fastgtpsa=true function bb_kick!(i, coords::Coords, u, order, kick_coeff, beta_x, beta_y, sigma_x, sigma_y, sigma_z, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  s_mid = v[i,ZI]/2 + order*sigma_z*erfinv(u + L/2)
  s_mid2 = s_mid*s_mid
  density_weight = (erf(sqrt(2)*order*erfinv(min(u + L, one(u + L)))) - erf(sqrt(2)*order*erfinv(u)))/4 # min prevents rounding overflow

  beta_x2 = beta_x*beta_x
  beta_y2 = beta_y*beta_y

  sigma_x_cur = sigma_x*sqrt(1 + s_mid2/beta_x2)
  sigma_y_cur = sigma_y*sqrt(1 + s_mid2/beta_y2)
  sigma_x_prime_times_sigma_x = sigma_x*sigma_x*s_mid/beta_x2
  sigma_y_prime_times_sigma_y = sigma_y*sigma_y*s_mid/beta_y2

  phi_x, phi_y, phi_xx, phi_yy = kick_coeff .* bb_derivatives(v[i,XI], v[i,YI], sigma_x_cur, sigma_y_cur)
  phi_z = sigma_x_prime_times_sigma_x*phi_xx + sigma_y_prime_times_sigma_y*phi_yy

  new_px = v[i,PXI] - 2*density_weight*phi_x
  new_py = v[i,PYI] - 2*density_weight*phi_y
  new_pz = v[i,PZI] -   density_weight*phi_z

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end