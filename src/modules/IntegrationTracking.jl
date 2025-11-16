#=

Tracking using symplectic integration with splits.

=#

macro def_integrator_struct(name)
  quote
    struct $(esc(name))
      order::Int
      num_steps::Int 
      ds_step::Float64
      radiation_damping_on::Bool
      radiation_fluctuations_on::Bool
  
      function $(esc(name))(; order::Int=4, num_steps::Int=-1, ds_step::Float64=-1.0, radiation_damping_on::Bool=false, radiation_fluctuations_on::Bool=false)
        _order = order
        _num_steps = num_steps
        _ds_step = ds_step
        if _order ∉ (2, 4, 6, 8)
          error("Symplectic integration only supports orders 2, 4, 6, and 8")
        elseif _num_steps == -1 && _ds_step == -1.0
          _num_steps = 1
        elseif _num_steps > 0 && _ds_step > 0
          error("Only one of num_steps or ds_step should be specified")
        elseif _num_steps < 1 && _ds_step <= 0
          error("Invalid step size")
        elseif _num_steps > 0
          _ds_step = -1.0
        elseif _ds_step > 0
          _num_steps = -1
        end
        return new(_order, _num_steps, _ds_step, radiation_damping_on, radiation_fluctuations_on)
      end
    end
  end
end

@def_integrator_struct(SplitIntegration)
@def_integrator_struct(MatrixKick)
@def_integrator_struct(BendKick)
@def_integrator_struct(SolenoidKick)
@def_integrator_struct(DriftKick)

module IntegrationTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions, ..SIMDMathFunctions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, STATE_ALIVE, STATE_LOST, @makekernel, Coords, vifelse, C_LIGHT, bessel01_RF

#
# ===============  I N T E G R A T O R S  ===============
#

@makekernel fastgtpsa=true function order_two_integrator!(i, coords::Coords, ker, params, ds_step, num_steps, f1, f2, L)
  ExactTracking.linear_bend_fringe!(i, coords, f1)
  for _ in 1:num_steps
    ker(i, coords, update_t0(ker, params, ds)..., ds_step)
  end
  ExactTracking.linear_bend_fringe!(i, coords, f2)
end


@makekernel fastgtpsa=true function order_four_integrator!(i, coords::Coords, ker, params, ds_step, num_steps, f1, f2, L)
  w0 = -1.7024143839193153215916254339390434324741363525390625*ds_step
  w1 =  1.3512071919596577718181151794851757586002349853515625*ds_step
  ds = 0.0
  ExactTracking.linear_bend_fringe!(i, coords, f1)
  for _ in 1:num_steps
    ker(i, coords, update_t0(ker, params, ds)..., w1)
    ds += w1
    ker(i, coords, update_t0(ker, params, ds)..., w0)
    ds += w0
    ker(i, coords, update_t0(ker, params, ds)..., w1)
    ds += w1
  end
  ExactTracking.linear_bend_fringe!(i, coords, f2)
end


@makekernel fastgtpsa=true function order_six_integrator!(i, coords::Coords, ker, params, ds_step, num_steps, f1, f2, L)
  w0 =  1.315186320683911169737712043570355*ds_step
  w1 = -1.17767998417887100694641568096432*ds_step
  w2 =  0.235573213359358133684793182978535*ds_step
  w3 =  0.784513610477557263819497633866351*ds_step
  ExactTracking.linear_bend_fringe!(i, coords, f1)
  for _ in 1:num_steps
    ker(i, coords, update_t0(ker, params, ds)..., w3)
    ker(i, coords, update_t0(ker, params, ds)..., w2)
    ker(i, coords, update_t0(ker, params, ds)..., w1)
    ker(i, coords, update_t0(ker, params, ds)..., w0)
    ker(i, coords, update_t0(ker, params, ds)..., w1)
    ker(i, coords, update_t0(ker, params, ds)..., w2)
    ker(i, coords, update_t0(ker, params, ds)..., w3)
  end
  ExactTracking.linear_bend_fringe!(i, coords, f2)
end


@makekernel fastgtpsa=true function order_eight_integrator!(i, coords::Coords, ker, params, ds_step, num_steps, f1, f2, L)
  w0 =  1.7084530707869978*ds_step
  w1 =  0.102799849391985*ds_step
  w2 = -1.96061023297549*ds_step
  w3 =  1.93813913762276*ds_step
  w4 = -0.158240635368243*ds_step
  w5 = -1.44485223686048*ds_step
  w6 =  0.253693336566229*ds_step
  w7 =  0.914844246229740*ds_step
  ExactTracking.linear_bend_fringe!(i, coords, f1)
  for _ in 1:num_steps
    ker(i, coords, update_t0(ker, params, ds)..., w7)
    ker(i, coords, update_t0(ker, params, ds)..., w6)
    ker(i, coords, update_t0(ker, params, ds)..., w5)
    ker(i, coords, update_t0(ker, params, ds)..., w4)
    ker(i, coords, update_t0(ker, params, ds)..., w3)
    ker(i, coords, update_t0(ker, params, ds)..., w2)
    ker(i, coords, update_t0(ker, params, ds)..., w1)
    ker(i, coords, update_t0(ker, params, ds)..., w0)
    ker(i, coords, update_t0(ker, params, ds)..., w1) 
    ker(i, coords, update_t0(ker, params, ds)..., w2)
    ker(i, coords, update_t0(ker, params, ds)..., w3)
    ker(i, coords, update_t0(ker, params, ds)..., w4)
    ker(i, coords, update_t0(ker, params, ds)..., w5)
    ker(i, coords, update_t0(ker, params, ds)..., w6)
    ker(i, coords, update_t0(ker, params, ds)..., w7)
  end
  ExactTracking.linear_bend_fringe!(i, coords, f2)
end



function update_t0(ker, params, ds)
  @FastGTPSA begin @inbounds begin
    if ker == cavity!
      t0 = params[13] + ds/(params[5]*C_LIGHT)
      new_params = (params[1], params[2], params[3], params[4], params[5], 
      params[6], params[7], params[8], params[9], params[10], params[11], 
      params[12], t0, params[14], params[15], params[16])
    else
      new_params = params
    end
  end end
  return new_params
end



#
# ===============  Q U A D R U P O L E  ===============
#
"""
mkm_quadrupole!()

This integrator uses Matrix-Kick-Matrix to implement a quadrupole
integrator accurate though second-order in the integration step-size. The vectors
kn and ks contain the normal and skew multipole strengths.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
ks: vector of skew multipole strengths scaled by Bρ0
L: element length
"""
@makekernel fastgtpsa=true function mkm_quadrupole!(i, coords::Coords, q, mc2, radiation_damping, radiation_fluctuations, beta_0, gamsqr_0, tilde_m, a, w, w_inv, k1, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  rel_p = 1 + coords.v[i,PZI]
  px = coords.v[i,PXI]
  py = coords.v[i,PYI]
  P_s2 = rel_p*rel_p - px*px - py*py
  good_momenta = (P_s2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])

  E0 = mc2/tilde_m/beta_0 # could probably exclude beta_0 because ultrarelativistic radiation

  #println(kn)

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, 0, tilde_m, mm, kn, ks, L / 2)
  end

  if radiation_damping
    deterministic_radiation!(   i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end
  #if radiation_fluctuations
  #  stochastic_radiation!(      i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  #end

  ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 2)
  quadrupole_kick!(             i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  BeamTracking.track_rotation!( i, coords, w, 0)
  quadrupole_matrix!(           i, coords, k1, L)
  BeamTracking.track_rotation!( i, coords, w_inv, 0)
  quadrupole_kick!(             i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 2)

  #if radiation_fluctuations
  #  stochastic_radiation!(      i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  #end
  if radiation_damping
    deterministic_radiation!(   i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, 0, tilde_m, mm, kn, ks, L / 2)
  end
end 


"""
quadrupole_matrix!()

Track "matrix part" of quadrupole.

Arguments
—————————
k1:  g / Bρ0 = g / (p0 / q)
         where g and Bρ0 respectively denote the quadrupole gradient
         and (signed) reference magnetic R_ref.
s: element length
"""
@makekernel fastgtpsa=true function quadrupole_matrix!(i, coords::Coords, k1, s)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  focus = k1 >= 0  # horizontally focusing if positive

  rel_p = 1 + v[i,PZI]
  xp = v[i,PXI] / rel_p  # x'
  yp = v[i,PYI] / rel_p  # y'
  sqrtks = sqrt(abs(k1 / rel_p)) * s  # |κ|s

  cosine = cos(sqrtks)
  coshine = cosh(sqrtks)
  sinecu = sincu(sqrtks)
  shinecu = sinhcu(sqrtks)
  cx = vifelse(focus, cosine, coshine)
  cy = vifelse(focus, coshine, cosine)
  sx = vifelse(focus, sinecu, shinecu)
  sy = vifelse(focus, shinecu, sinecu)

  new_px = v[i,PXI] * cx - k1 * s * v[i,XI] * sx
  new_py = v[i,PYI] * cy + k1 * s * v[i,YI] * sy
  new_z = v[i,ZI]  - (s / 4) * (  xp*xp * (1 + sx * cx)
                                    + yp*yp * (1 + sy * cy)
                                    + k1 / (1 + v[i,PZI])
                                        * ( v[i,XI]*v[i,XI] * (1 - sx * cx)
                                          - v[i,YI]*v[i,YI] * (1 - sy * cy) )
                                  ) + sign(k1) * ( v[i,XI] * xp * (sqrtks * sx)* 
                                  (sqrtks * sx) - v[i,YI] * yp * (sqrtks * sy)*
                                  (sqrtks * sy) ) / 2
  new_x = v[i,XI] * cx + xp * s * sx
  new_y = v[i,YI] * cy + yp * s * sy
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,ZI]  = vifelse(alive, new_z, v[i,ZI])
  v[i,XI]  = vifelse(alive, new_x, v[i,XI])
  v[i,YI]  = vifelse(alive, new_y, v[i,YI])
end 


"""
quadrupole_kick!()

Track "remaining part" of quadrupole —— a position kick.

### Note re implementation:
A common factor that appears in the expressions for `zf.x` and `zf.y`
originally included a factor with the generic form ``1 - \\sqrt{1 - A}``,
which suffers a loss of precision when ``|A| \\ll 1``. To combat that
problem, we rewrite it in the form ``A / (1 + \\sqrt{1-A})``---more
complicated, yes, but far more accurate.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
s: element length
"""
@makekernel fastgtpsa=true function quadrupole_kick!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, s)
  v = coords.v

  P      = 1 + v[i,PZI]             # [scaled] total momentum, P/P0 = 1 + δ
  PtSqr  = v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]  # (transverse momentum)^2, P⟂^2 = (Px^2 + Py^2) / P0^2
  Ps2    = P*P - PtSqr        
  good_momenta = (Ps2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  Ps2_1 = one(Ps2)
  Ps = sqrt(vifelse(good_momenta, Ps2, Ps2_1)) # longitudinal momentum,   Ps   = √[(1 + δ)^2 - P⟂^2]
  alive = (coords.state[i] == STATE_ALIVE)

  new_x = v[i,XI] + s * v[i,PXI] * PtSqr / (P * Ps * (P + Ps))
  new_y = v[i,YI] + s * v[i,PYI] * PtSqr / (P * Ps * (P + Ps))
  new_z = v[i,ZI] - s * (P * (PtSqr - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                                / ( beta_0 * sqrt(P*P + tilde_m*tilde_m) * Ps
                                    * (beta_0 * sqrt(P*P + tilde_m*tilde_m) + Ps)
                                  )
                            - PtSqr / (2 * P*P))
  v[i,XI] = vifelse(alive, new_x, v[i,XI])
  v[i,YI] = vifelse(alive, new_y, v[i,YI])
  v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
end # function quadrupole_kick!()


#
# ===============  B E N D  ===============
#

"""
bkb_multipole!()

This integrator uses Bend-Kick-Bend to track a beam through
a curved, finite-length multipole magnet. The method is
accurate through second order in the step size. The vectors
kn and ks contain the normal and skew multipole strengths,
starting with the quadrupole component.

Arguments
—————————
- 'tilde_m'  -- mc2/p0c
- 'beta_0'   -- p0c/E0
- 'theta'    -- 'g' * 'L'
- 'g'        -- curvature
- 'w'        -- rotation matrix into curvature/field plane
- 'w_inv'    -- rotation matrix out of curvature/field plane
- 'k0'       -- dipole strength
- 'mm'       -- order of multipoles
- 'kn'       -- normal multipole strengths 
- 'ks'       -- skew multipole strengths 
- 'L'        -- length
"""
@makekernel fastgtpsa=true function bkb_multipole!(i, coords::Coords, q, mc2, radiation_damping, radiation_fluctuations, tilde_m, beta_0, a, g, w, w_inv, k0, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  E0 = mc2/tilde_m/beta_0 # could probably exclude beta_0 because ultrarelativistic radiation

  BeamTracking.track_rotation!( i, coords, w, 0)

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, g, tilde_m, mm, kn, ks, L / 2)
  end

  if radiation_damping
    deterministic_radiation!(   i, coords, q, mc2, E0, g, mm, kn, ks, L / 2)
  end
  #if radiation_fluctuations
  #  stochastic_radiation!(      i, coords, q, mc2, E0, g, mm, kn, ks, L / 2)
  #end

  ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 1)
  ExactTracking.exact_bend!(    i, coords, g*L, g, k0, tilde_m, beta_0, L)
  ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 1)

  #if radiation_fluctuations
  #  stochastic_radiation!(      i, coords, q, mc2, E0, g, mm, kn, ks, L / 2)
  #end
  if radiation_damping
    deterministic_radiation!(   i, coords, q, mc2, E0, g, mm, kn, ks, L / 2)
  end

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, g, tilde_m, mm, kn, ks, L / 2)
  end

  BeamTracking.track_rotation!( i, coords, w_inv, 0)
end 


#
# ===============  S O L E N O I D  ===============
#

"""
sks_multipole!()

This integrator uses Solenoid-Kick-Solenoid to track a beam through
a straight, finite-length multipole magnet with a solenoid field. The method is
accurate through second order in the step size. The vectors
kn and ks contain the normal and skew multipole strengths,
starting with the dipole component.

Arguments
—————————
Ksol: solenoid strength
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
sn: vector of skew multipole strengths scaled by Bρ0
L:  element length
"""
@makekernel fastgtpsa=true function sks_multipole!(i, coords::Coords, q, mc2, radiation_damping, radiation_fluctuations, beta_0, gamsqr_0, tilde_m, a, Ksol, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  E0 = mc2/tilde_m/beta_0 # could probably exclude beta_0 because ultrarelativistic radiation

  ExactTracking.exact_solenoid!(  i, coords, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)

  if radiation_damping
    deterministic_radiation!(     i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end
  #if radiation_fluctuations
  #  stochastic_radiation!(        i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  #end

  if isnothing(coords.q)
    ExactTracking.multipole_kick!(i, coords, mm, knl * 2, ksl * 2, -1)
  else
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, -1)
    rotate_spin!(                 i, coords, a, 0, tilde_m, mm, kn, ks, L)
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, -1)
  end

  #if radiation_fluctuations
  #  stochastic_radiation!(        i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  #end
  if radiation_damping
    deterministic_radiation!(     i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end

  ExactTracking.exact_solenoid!(  i, coords, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)
end 


#
# ===============  D R I F T  ===============
#
"""
dkd_multipole!()

This integrator uses Drift-Kick-Drift to track a beam through
a straight, finite-length multipole magnet. The method is
accurate through second order in the step size. The vectors
kn and ks contain the normal and skew multipole strengths,
starting with the dipole component. (For example, kn[3] denotes
the normal sextupole strength in Tesla/m^2.) The argument ns
denotes the number of slices.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
ks: vector of skew multipole strengths scaled by Bρ0
L:  element length
"""
@makekernel fastgtpsa=true function dkd_multipole!(i, coords::Coords, q, mc2, radiation_damping, radiation_fluctuations, beta_0, gamsqr_0, tilde_m, a, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  E0 = mc2/tilde_m/beta_0 # could probably exclude beta_0 because ultrarelativistic radiation

  ExactTracking.exact_drift!(     i, coords, beta_0, gamsqr_0, tilde_m, L / 2)

  if radiation_damping
    deterministic_radiation!(     i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end
  #if radiation_fluctuations
  #  stochastic_radiation!(        i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  #end

  if isnothing(coords.q)
    ExactTracking.multipole_kick!(i, coords, mm, knl * 2, ksl * 2, -1)
  else
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, -1)
    rotate_spin!(                 i, coords, a, 0, tilde_m, mm, kn, ks, L)
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, -1)
  end

  #if radiation_fluctuations
  #  stochastic_radiation!(        i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  #end
  if radiation_damping
    deterministic_radiation!(     i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end

  ExactTracking.exact_drift!(     i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
end


#
# ===============  S P I N  ===============
#
"""
This function computes the integrated spin-precession vector using the multipole 
coefficients kn and ks indexed by mm, i.e., knl[i] is the normal 
coefficient of order mm[i].
"""
function omega_multipole(i, coords::Coords, a, g, tilde_m, mm, kn, ks, L)
  @FastGTPSA begin @inbounds begin
    v = coords.v

    rel_p = 1 + v[i,PZI]
    beta_gamma = rel_p / tilde_m
    gamma = sqrt(1 + beta_gamma*beta_gamma)
    beta = beta_gamma / gamma

    if mm[1] == 0
      ax = -v[i,YI] * kn[1] / 2
      ay =  v[i,XI] * kn[1] / 2
    else
      ax = zero(v[i,XI])
      ay = ax
    end

    bx, by = ExactTracking.normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
    bz_0 = zero(kn[1])
    bz = mm[1] == 0 ? kn[1] : bz_0
    b_vec = (bx, by, bz)
    e_vec = (bz_0, bz_0, bz_0)

    omega = omega_field(i, coords, a, g, beta, gamma, ax, ay, e_vec, b_vec, L)
  end end

  return omega
end


"""
This function computes the integrated spin-precession vector using the fields.
"""
function omega_field(i, coords::Coords, a, g, beta, gamma, ax, ay, e_vec, b_vec, L)
  @FastGTPSA begin @inbounds begin
    v = coords.v
    px = v[i,PXI] - ax
    py = v[i,PYI] - ay
    rel_p = 1 + v[i,PZI]

    pl2 = rel_p*rel_p - px*px - py*py
    pl2_0 = zero(pl2)
    good_momenta = (pl2 > pl2_0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    pl2_1 = one(pl2)
    pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

    coeff = -(1 + g*v[i,XI])/pl
    coeff1 = coeff * (1 + a*gamma)
    coeff2 = coeff * (1 + a)
    coeff3 = -coeff * beta / C_LIGHT * gamma * (a + 1/(1+gamma))/rel_p

    betax = px / rel_p
    betay = py / rel_p
    betaz = pl / rel_p

    dot = b_vec[1]*betax + b_vec[2]*betay + b_vec[3]*betaz

    b_para_x = dot * betax
    b_para_y = dot * betay
    b_para_z = dot * betaz

    b_perp_x = (b_vec[1] - b_para_x) * coeff1
    b_perp_y = (b_vec[2] - b_para_y) * coeff1
    b_perp_z = (b_vec[3] - b_para_z) * coeff1

    b_para_x = b_para_x * coeff2
    b_para_y = b_para_y * coeff2
    b_para_z = b_para_z * coeff2

    e_part_x = (py*e_vec[3] - pl*e_vec[2]) * coeff3
    e_part_y = (pl*e_vec[1] - px*e_vec[3]) * coeff3
    e_part_z = (px*e_vec[2] - py*e_vec[1]) * coeff3

    ox = (b_perp_x + b_para_x + e_part_x) * L        
    oy = (b_perp_y + b_para_y + e_part_y + g) * L
    oz = (b_perp_z + b_para_z + e_part_z) * L

    omega = (ox, oy, oz)
  end end
  return omega
end


"""
This function rotates particle i's quaternion according to the multipoles present.
"""
@makekernel fastgtpsa=true function rotate_spin!(i, coords::Coords, a, g, tilde_m, mm, kn, ks, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_multipole(i, coords, a, g, tilde_m, mm, kn, ks, L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end


@makekernel fastgtpsa=true function integrate_with_spin_thin!(i, coords::Coords, ker, params, a, g, tilde_m, mm, knl, ksl)
  rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
  ker(i, coords, update_t0(ker, params, ds)...)
  rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
end


#
# ===============  R F  ===============
#
@makekernel fastgtpsa=true function cavity!(i, coords::Coords, q, mc2, radiation_damping, radiation_fluctuations, beta_0, gamsqr_0, tilde_m, E_ref, p0c, a, omega, E0_over_Rref, t0, mm, kn, ks, L)
  cavity_fringe!(i, coords, beta_0, tilde_m, omega, E0_over_Rref, t0)
  multipoles = (length(mm) > 0)
  sol = (multipoles && mm[1] == 0)
  if sol
    ExactTracking.exact_solenoid!(i, coords, kn[1], beta_0, gamsqr_0, tilde_m, L / 2)
  else
    ExactTracking.exact_drift!(   i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  end
  t0 = t0 + (L/2)/(beta_0*C_LIGHT)

  if multipoles
    if radiation_damping
      deterministic_radiation!(   i, coords, q, mc2, E_ref, 0, mm, kn, ks, L / 2)
    end
    #if radiation_fluctuations
    #  stochastic_radiation!(      i, coords, q, mc2, E_ref, 0, mm, kn, ks, L / 2)
    #end
    ExactTracking.multipole_kick!(i, coords, mm, kn * L / 2, ks * L / 2, -1)
  end

  if isnothing(coords.q)
    cavity_kick!(                 i, coords, beta_0, tilde_m, E_ref, p0c, omega, E0_over_Rref, t0, L)
  else
    cavity_kick!(                 i, coords, beta_0, tilde_m, E_ref, p0c, omega, E0_over_Rref, t0, L / 2)
    rotate_spin_cavity!(          i, coords, a, tilde_m, omega, E0_over_Rref, t0, mm, kn, ks, L)
    cavity_kick!(                 i, coords, beta_0, tilde_m, E_ref, p0c, omega, E0_over_Rref, t0, L / 2)
  end

  if multipoles
    ExactTracking.multipole_kick!(i, coords, mm, kn * L / 2, ks * L / 2, -1)
    #if radiation_fluctuations
    #  stochastic_radiation!(      i, coords, q, mc2, E_ref, 0, mm, kn, ks, L / 2)
    #end
    if radiation_damping
      deterministic_radiation!(   i, coords, q, mc2, E_ref, 0, mm, kn, ks, L / 2)
    end
  end

  if sol
    ExactTracking.exact_solenoid!(i, coords, kn[1], beta_0, gamsqr_0, tilde_m, L / 2)
  else
    ExactTracking.exact_drift!(   i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  end
  t0 = t0 + (L/2)/(beta_0*C_LIGHT)
  cavity_fringe!(i, coords, beta_0, tilde_m, omega, -E0_over_Rref, t0)
end


@makekernel fastgtpsa=true function bmad_to_mad!(i, coords::Coords, beta_0, tilde_m, E_ref, p0c)
  v = coords.v

  rel_p = 1 + v[i,PZI]
  beta_gamma = rel_p/tilde_m
  gamma = sqrt(1 + beta_gamma*beta_gamma)
  beta = beta_gamma/gamma
  tau = v[i,ZI]/beta

  gamma_0_inv = tilde_m*beta_0
  E = E_ref*gamma*gamma_0_inv

  v[i,ZI]  = tau
  v[i,PZI] = E/p0c - 1/beta_0
end


@makekernel fastgtpsa=true function mad_to_bmad!(i, coords::Coords, beta_0, tilde_m, E_ref, p0c)
  v = coords.v

  E = E_ref + p0c*v[i,PZI]
  gamma_0_inv = tilde_m*beta_0
  gamma = E/E_ref/gamma_0_inv
  beta = sqrt(1-1/(gamma*gamma))
  z = v[i,ZI]*beta
  
  pc = beta*E

  v[i,ZI]  =  z
  v[i,PZI] = (pc-p0c)/p0c
end


@makekernel fastgtpsa=true function cavity_kick!(i, coords::Coords, beta_0, tilde_m, E_ref, p0c, omega, E0_over_Rref, t0, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  bmad_to_mad!(i, coords, beta_0, tilde_m, E_ref, p0c)
  r2 = v[i,XI]*v[i,XI] + v[i,YI]*v[i,YI]
  b01 = 2.404825557695773 # first zero of J0
  d = C_LIGHT*b01/omega
  arg = (b01*b01)/(d*d)*r2
  b0, b1 = bessel01_RF(arg)
  b1 = b1 * b01/d

  t = t0 - v[i,ZI]/C_LIGHT

  px_0 = v[i,PXI]
  py_0 = v[i,PYI]
  pz_0 = v[i,PZI]

  phi_particle = omega*t
  s, c = sincos(phi_particle)

  coeff = L*E0_over_Rref*b01/(omega*d)*b1*c

  new_px = px_0 - coeff*v[i,XI]
  new_py = py_0 - coeff*v[i,YI]
  new_pz = pz_0 + L*E0_over_Rref*b0*s/C_LIGHT

  v[i,PXI] = vifelse(alive, new_px, px_0)
  v[i,PYI] = vifelse(alive, new_py, py_0)
  v[i,PZI] = vifelse(alive, new_pz, pz_0)

  mad_to_bmad!(i, coords, beta_0, tilde_m, E_ref, p0c)
end


@makekernel fastgtpsa=true function cavity_fringe!(i, coords::Coords, beta_0, tilde_m, omega, E0_over_Rref, t0)
  v = coords.v 
  alive = (coords.state[i] == STATE_ALIVE)

  r2 = v[i,XI]*v[i,XI] + v[i,YI]*v[i,YI]
  b01 = 2.404825557695773 # first zero of J0
  d = C_LIGHT*b01/omega
  arg = (b01*b01)/(d*d)*r2
  b0, b1 = bessel01_RF(arg)
  #b1 = b1 * b01/d I don't want that here
  beta_gamma = (1 + v[i,PZI])/tilde_m
  gamma = sqrt(1 + beta_gamma*beta_gamma)
  beta = beta_gamma/gamma
  vel = beta*C_LIGHT
  t = t0 - v[i,ZI]/vel

  px_0 = v[i,PXI]
  py_0 = v[i,PYI]

  phi_particle = omega*t
  s, c = sincos(phi_particle)

  coeff = E0_over_Rref*b1*s/vel

  new_px = px_0 - coeff*v[i,XI]
  new_py = py_0 - coeff*v[i,YI]

  v[i,PXI] = vifelse(alive, new_px, px_0)
  v[i,PYI] = vifelse(alive, new_py, py_0)
end


function omega_cavity(i, coords::Coords, a, tilde_m, omega, E0_over_Rref, t0, mm, kn, ks, L)
  @FastGTPSA begin @inbounds begin
    v = coords.v
    alive = (coords.state[i] == STATE_ALIVE)
    r2 = v[i,XI]*v[i,XI] + v[i,YI]*v[i,YI]
    b01 = 2.404825557695773 # first zero of J0
    d = C_LIGHT*b01/omega
    arg = (b01*b01)/(d*d)*r2
    b0, b1 = bessel01_RF(arg)
    b1 = b1 * b01/d
    beta_gamma = (1 + v[i,PZI])/tilde_m
    gamma = sqrt(1 + beta_gamma*beta_gamma)
    beta = beta_gamma/gamma
    vel = beta*C_LIGHT
    t = t0 - v[i,ZI]/vel

    phi_particle = omega*t
    s, c = sincos(phi_particle)

    ez = E0_over_Rref*sin(phi_particle)
    ex = zero(ez)
    ey = ex
    e_vec = (ex, ey, ez)

    coeff = E0_over_Rref/C_LIGHT*b1*c

    bx = -coeff*v[i,YI]
    by = coeff*v[i,XI]
    bz = ex
    b_vec = (bx, by, bz)

    if length(mm) > 0 && mm[1] == 0
      ax = -v[i,YI] * kn[1] / 2
      ay =  v[i,XI] * kn[1] / 2
    else
      ax = ex
      ay = ex
    end

    ox, oy, oz = omega_field(i, coords, a, 0, beta, gamma, ax, ay, e_vec, b_vec, L)
    if length(mm) > 0
      ox1, oy1, oz1 = omega_multipole(i, coords, a, 0, tilde_m, mm, kn, ks, L)
      omega = (ox + ox1, oy + oy1, oz + oz1)
    else
      omega = (ox, oy, oz)
    end
  end end
  return omega
end


"""
This function rotates particle i's quaternion in a cavity.
"""
@makekernel fastgtpsa=true function rotate_spin_cavity!(i, coords::Coords, a, tilde_m, omega, E0_over_Rref, t0, mm, kn, ks, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_cavity(i, coords, a, tilde_m, omega, E0_over_Rref, t0, mm, kn, ks, L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end


#
# ===============  R A D I A T I O N  ===============
#
@makekernel fastgtpsa=true function canonical_to_prime!(i, coords::Coords, g, ax, ay)
  v = coords.v

  rel_p = 1 + v[i,PZI]
  px = v[i,PXI] - ax
  py = v[i,PYI] - ay

  pl2 = rel_p*rel_p - px*px - py*py
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  h = (1 + g*v[i,XI])/pl

  new_px = h*px 
  new_py = h*py

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
end


@makekernel fastgtpsa=true function prime_to_canonical!(i, coords::Coords, g, ax, ay)
  v = coords.v

  h = 1 + g*v[i,XI]
  rel_p = 1 + v[i,PZI]

  pl2 = h*h + v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  new_px = rel_p*v[i,PXI]/pl + ax
  new_py = rel_p*v[i,PYI]/pl + ay

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
end


@makekernel fastgtpsa=true function deterministic_radiation!(i, coords::Coords, q, mc2, E0, g, mm, kn, ks, L)
  v = coords.v

  if mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
  else
    ax = zero(v[i,XI])
    ay = ax
  end

  canonical_to_prime!(i, coords, g, ax, ay)

  h = 1 + g*v[i,XI]
  rel_p = 1 + v[i,PZI]

  pl2 = h*h + v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  bx, by = ExactTracking.normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
  bz_0 = zero(kn[1])
  bz = mm[1] == 0 ? kn[1] : bz_0

  betax = v[i,PXI] / pl
  betay = v[i,PYI] / pl
  betaz = h / pl

  dot = bx*betax + by*betay + bz*betaz

  b_perp_x = bx - dot*betax
  b_perp_y = by - dot*betay
  b_perp_z = bz - dot*betaz

  b_perp_2 = b_perp_x*b_perp_x + b_perp_y*b_perp_y + b_perp_z*b_perp_z

  coeff = E_CHARGE/(4*pi*EPS_0) # 8.9875517862e9 * 1.602176634e-19

  K = -pl * coeff * 2/3 * (q*q)/(mc2*mc2*mc2*mc2) * (E0*E0*E0) * b_perp_2 * L

  new_pz = (v[i,PZI] + rel_p*K)/(1 - rel_p*K)
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])

  prime_to_canonical!(i, coords, g, ax, ay)
end

#=
@makekernel fastgtpsa=true function stochastic_radiation!(i, coords::Coords, q, mc2, E0, g, mm, kn, ks, L)
  v = coords.v

  if mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
  else
    ax = zero(v[i,XI])
    ay = ax
  end

  h = 1 + g*v[i,XI]
  rel_p = 1 + v[i,PZI]
  gamma = rel_p * E0 / mc2
  px = v[i,PXI] - ax
  py = v[i,PYI] - ay

  pl2 = rel_p*rel_p - px*px - py*py
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  bx, by = ExactTracking.normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
  bz_0 = zero(kn[1])
  bz = mm[1] == 0 ? kn[1] : bz_0

  betax = v[i,PXI] / pl
  betay = v[i,PYI] / pl
  betaz = h / pl

  dot = bx*betax + by*betay + bz*betaz

  b_perp_x = bx - dot*betax
  b_perp_y = by - dot*betay
  b_perp_z = bz - dot*betaz

  b_perp_2 = b_perp_x*b_perp_x + b_perp_y*b_perp_y + b_perp_z*b_perp_z
  b_perp = sqrt(b_perp_2)

  dt_ds = h * rel_p / pl

  coeff = 55/(24*sqrt(3))*1/(4*pi*EPS_0)*H_BAR*C_LIGHT

  mc27 = mc2*mc2*mc2*mc2*mc2*mc2*mc2
  E05 = E0*E0*E0*E0*E0
  rel_p4 = rel_p*rel_p*rel_p*rel_p
  b_perp_3 = b_perp_2*b_perp
  q2 = q*q

  sigma2 = dt_ds * coeff * q2/mc27 * E05 * rel_p4 * b_perp_3 * abs(L)

  dpz   = randn() * sqrt(sigma2)
  theta = randn() / gamma
  s, c  = sincos(theta)

  b_perp_hat_x = b_perp_x / b_perp
  b_perp_hat_y = b_perp_y / b_perp

  new_px = v[i,PXI] + dpz * (c*betax + s*b_perp_hat_x)
  new_py = v[i,PYI] + dpz * (c*betay + s*b_perp_hat_y)
  new_pz = v[i,PZI] + dpz

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end
=#

end