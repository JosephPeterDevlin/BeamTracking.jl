# =========== HELPER FUNCTIONS ============= #
@inline function get_thin_multipoles(bunch, bm, L)
  brho_0 = bunch.Brho_ref
  mm = sort(collect(keys(bm)))
  kn = [get_thin_strength(bm[i], L, brho_0)*cos(-i*bm[i].tilt) for i in mm]
  sn = [get_thin_strength(bm[i], L, brho_0)*sin(-i*bm[i].tilt) for i in mm]
  return mm, kn, sn
end

@inline function get_thick_multipoles(bunch, bm, L)
  brho_0 = bunch.Brho_ref
  mm = sort(collect(keys(bm)))
  kn = [get_thick_strength(bm[i], L, brho_0)*cos(-i*bm[i].tilt) for i in mm]
  sn = [get_thick_strength(bm[i], L, brho_0)*sin(-i*bm[i].tilt) for i in mm]
  return mm, kn, sn
end

@inline function get_thick_multipoles_no_dipole(bunch, bm, L)
  brho_0 = bunch.Brho_ref
  mm = sort(setdiff(collect(keys(bm))))[2:end]
  kn = [get_thick_strength(bm[i], L, brho_0)*cos(-i*bm[i].tilt) for i in mm]
  sn = [get_thick_strength(bm[i], L, brho_0)*sin(-i*bm[i].tilt) for i in mm]
  return mm, kn, sn
end

@inline function integration_launcher!(ker, params, order, n_steps, L)
  ds_step = L/n_steps
  if order == 2
    return KernelCall(IntegrationTracking.order_two_integrator!, (ker, params, ds_step, n_steps, L))
  elseif order == 4
    return KernelCall(IntegrationTracking.order_four_integrator!, (ker, params, ds_step, n_steps, L))
  elseif order == 6
    return KernelCall(IntegrationTracking.order_six_integrator!, (ker, params, ds_step, n_steps, L))
  elseif order == 8
    return KernelCall(IntegrationTracking.order_eight_integrator!, (ker, params, ds_step, n_steps, L))
  else
    error("Symplectic integration only supports orders 2, 4, 6, and 8")
  end
end


# =========== STRAIGHT ELEMENTS ============= #
# === Thin elements === #
@inline function thin_pure_bdipole(tm::Integration, bunch, bm1)
  brho_0 = bunch.Brho_ref
  mm = @SArray [bm1.order]
  kn = @SArray [get_thin_strength(bm1, 0.0, brho_0)*cos(-bm1.order*bm1.tilt)]
  sn = @SArray [get_thin_strength(bm1, 0.0, brho_0)*sin(-bm1.order*bm1.tilt)]
  return KernelCall(ExactTracking.multipole_kick!, (mm, kn, sn))
end

@inline function thin_bdipole(tm::Integration, bunch, bm)
  return KernelCall(ExactTracking.multipole_kick!, get_thin_multipoles(bunch,bm,0.0))
end

@inline thin_pure_bquadrupole(tm::Integration, bunch, bm1) = thin_pure_bdipole(tm, bunch, bm1)

@inline thin_bquadrupole(tm::Integration, bunch, bm1) = thin_bdipole(tm, bunch, bm1)

@inline thin_pure_bmultipole(tm::Integration, bunch, bm1) = thin_pure_bdipole(tm, bunch, bm1)

@inline thin_bmultipole(tm::Integration, bunch, bm1) = thin_bdipole(tm, bunch, bm1)


# === Thick elements === #
@inline drift(tm::Integration, bunch, L) = drift(Exact(), bunch, L)

@inline function thick_pure_bsolenoid(tm::Integration, bunch, bm0, L)
  if isnothing(bunch.q)
    return thick_pure_bsolenoid(Exact(), bunch, bm0, L)
  else
    tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
    Ks = get_thick_strength(bm0, L, bunch.Brho_ref)
    params = (Ks, beta_0, gamsqr_0, tilde_m)
    return integration_launcher!(ExactTracking.solenoid_spin!, params, tm.order, tm.n_steps, L)
  end
end

@inline function thick_bsolenoid(tm::Integration, bunch, bm, L)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  Ks = get_thick_strength(bm[0], L, bunch.Brho_ref)
  mm, kn, sn = get_thick_multipoles(bunch, bm, L)
  params = (Ks, beta_0, gamsqr_0, tilde_m, mm, kn, sn)
  return integration_launcher!(ExactTracking.sks_multipole!, params, tm.order, tm.n_steps, L)
end

@inline function thick_pure_bdipole(tm::Integration, bunch, bm1, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  b0 = brho_0 * get_thick_strength(bm1, L, brho_0)
  return KernelCall(ExactTracking.exact_sbend!, (beta_0, brho_0, 0.0, b0, 0.0, 0.0, L))
end

@inline function thick_bdipole(tm::Integration, bunch, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  b0 = brho_0 * get_thick_strength(bm[1], L, brho_0)
  mm, kn, sn = get_thick_multipoles_no_dipole(bunch, bm, L)
  params = (beta_0, brho_0, 0.0, b0, 0.0, 0.0, mm, kn, sn)
  return integration_launcher!(ExactTracking.bkb_multipole!, params, tm.order, tm.n_steps, L)
end

@inline function thick_pure_bquadrupole(tm::Integration, bunch, bm2, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  K1 = get_thick_strength(bm2, L, brho_0)
  params = (beta_0, gamsqr_0, tilde_m, K1)
  return integration_launcher!(ExactTracking.mkm_quadrupole!, params, tm.order, tm.n_steps, L)
end

@inline function thick_bquadrupole(tm::Integration, bunch, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  mm, kn, sn = get_thick_multipoles(bunch, bm, L)
  params = (beta_0, gamsqr_0, tilde_m, mm, kn, sn)
  return integration_launcher!(ExactTracking.dkd_multipole!, params, tm.order, tm.n_steps, L)
end

@inline function thick_pure_bmultipole(tm::Integration, bunch, bm1, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  mm = @SArray [bm1.order]
  kn = @SArray [get_thick_strength(bm1, L, bunch.Brho_ref)*cos(-bm1.order*bm1.tilt)]
  sn = @SArray [get_thick_strength(bm1, L, bunch.Brho_ref)*sin(-bm1.order*bm1.tilt)]
  params = (beta_0, gamsqr_0, tilde_m, mm, kn, sn)
  return integration_launcher!(ExactTracking.dkd_multipole!, params, tm.order, tm.n_steps, L)
end

@inline thick_bmultipole(tm::Integration, bunch, bm, L) = thick_bquadrupole(tm, bunch, bm, L)


# =========== BENDING ELEMENTS ============= #
@inline function thick_bend_no_field(tm::Integration, bunch, bendparams, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  return KernelCall(ExactTracking.exact_sbend!, (beta_0, brho_0, hc, 0.0, e1, e2, L))
end

@inline function thick_bend_pure_bdipole(tm::Integration, bunch, bendparams, bm1, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  b0 = brho_0 * get_thick_strength(bm1, L, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  return KernelCall(ExactTracking.exact_sbend!, (beta_0, brho_0, hc, b0, e1, e2, L))
end

@inline function thick_bend_bdipole(tm::Integration, bunch, bendparams, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  b0 = brho_0 * get_thick_strength(bm.bdict[1], L, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  mm, kn, sn = get_thick_multipoles_no_dipole(bunch, bm, L)
  params = (beta_0, brho_0, hc, b0, e1, e2, mm, kn, sn)
  return integration_launcher!(ExactTracking.bkb_multipole!, params, tm.order, tm.n_steps, L)
end

@inline function thick_bend_pure_bquadrupole(tm::Integration, bunch, bendparams, bm2, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  K1 = get_thick_strength(bm2, L, brho_0)
  mm = @SArray [bm2.order]
  kn = @SArray [K1 * cos(-bm2.order*bm2.tilt)]
  sn = @SArray [K1 * sin(-bm2.order*bm2.tilt)]
  params = (beta_0, brho_0, hc, 0.0, e1, e2, mm, kn, sn)
  return integration_launcher!(ExactTracking.bkb_multipole!, params, tm.order, tm.n_steps, L)
end

@inline function thick_bend_bquadrupole(tm::Integration, bunch, bendparams, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  mm, kn, sn = get_thick_multipoles(bunch, bm, L)
  params = (beta_0, brho_0, hc, 0.0, e1, e2, mm, kn, sn)
  return integration_launcher!(ExactTracking.bkb_multipole!, params, tm.order, tm.n_steps, L)
end

@inline thick_bend_pure_bmultipole(tm::Integration, bunch, bendparams, bm1, L) = thick_bend_pure_bquadrupole(tm, bunch, bendparams, bm1, L)

@inline thick_bend_bmultipole(tm::Integration, bunch, bendparams, bdict, L) = thick_bend_bquadrupole(tm, bunch, bendparams, bdict, L)
