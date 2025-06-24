#=

Tracking using symplectic integration.

=#

macro def_integrator_struct(name)
  quote
    struct $(esc(name))
      order::Int
      num_steps::Int 
      ds_step::Float64
  
      function $(esc(name))(; order::Int=2, num_steps::Int=-1, ds_step::Float64=-1.0)
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
        return new(_order, _num_steps, _ds_step)
      end
    end
  end
end

@def_integrator_struct(Standard)
@def_integrator_struct(MKM)
@def_integrator_struct(BKB)
@def_integrator_struct(SKS)
@def_integrator_struct(DKD)

module IntegrationTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, @makekernel, BunchView

#
# ===============  I N T E G R A T O R S  ===============
#

@makekernel function order_two_integrator!(i, b::BunchView, ker, params, ds_step, num_steps, L)
  for _ in 1:num_steps
    ker(i, b, params..., ds_step)
  end
end


@makekernel function order_four_integrator!(i, b::BunchView, ker, params, ds_step, num_steps, L)
  w0 = 1.3512071919596577718181151794851757586002349853515625*ds_step
  w1 = -1.7024143839193153215916254339390434324741363525390625*ds_step
  for _ in 1:num_steps
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
  end
end


@makekernel function order_six_integrator!(i, b, ker, params, ds_step, num_steps, L)
  w0 = 1.3151863206857402*ds_step
  w1 = -1.17767998417887*ds_step
  w2 = 0.235573213359*ds_step
  w3 = 0.784513610477*ds_step
  for _ in 1:num_steps
    ker(i, b, params..., w3)
    ker(i, b, params..., w2)
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
    ker(i, b, params..., w2)
    ker(i, b, params..., w3)
  end
end


@makekernel function order_eight_integrator!(i, b, ker, params, ds_step, num_steps, L)
  w0 = -1.7808286265894515*ds_step
  w1 = -1.61582374150097*ds_step
  w2 = -2.44699182370524*ds_step
  w3 = -0.716989419708120e-2*ds_step
  w4 = 2.44002732616735*ds_step
  w5 = 0.157739928123617*ds_step
  w6 = 1.82020630970714*ds_step
  w7 = 1.04242620869991*ds_step
  for _ in 1:num_steps
    ker(i, b, params..., w7)
    ker(i, b, params..., w6)
    ker(i, b, params..., w5)
    ker(i, b, params..., w4)
    ker(i, b, params..., w3)
    ker(i, b, params..., w2)
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
    ker(i, b, params..., w2)
    ker(i, b, params..., w3)
    ker(i, b, params..., w4)
    ker(i, b, params..., w5)
    ker(i, b, params..., w6)
    ker(i, b, params..., w7)
  end
end


#
# ===============  Q U A D R U P O L E  ===============
#
"""
mkm_quadrupole!()

This integrator uses Matrix-Kick-Matrix to implement a quadrupole
integrator accurate though second-order in the integration step-size. The vectors
kn and ks contain the normal and skew multipole strengths,
excluding the quadrupole component.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
k1:   g / Bρ0 = g / (p0 / q)
          where g and Bρ0 respectively denote the quadrupole gradient
          and (signed) reference magnetic rigidity.
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
sn: vector of skew multipole strengths scaled by Bρ0
L: element length
"""
@makekernel fastgtpsa=true function mkm_quadrupole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, k1, mm, kn, sn, L)
  quadrupole_matrix!(i, b, k1, L / 2)
  if length(mm) > 0
    quadrupole_kick!(  i, b, beta_0, gamsqr_0, tilde_m, L/2)
    ExactTracking.multipole_kick!(i, b, mm, kn * L, sn * L)
    quadrupole_kick!(  i, b, beta_0, gamsqr_0, tilde_m, L/2)
  else
    quadrupole_kick!(  i, b, beta_0, gamsqr_0, tilde_m, L)
  end
  quadrupole_matrix!(i, b, k1, L / 2)
end 


"""
quadrupole_matrix!()

Track "matrix part" of quadrupole.

Arguments
—————————
k1:  g / Bρ0 = g / (p0 / q)
         where g and Bρ0 respectively denote the quadrupole gradient
         and (signed) reference magnetic rigidity.
s: element length
"""
@makekernel fastgtpsa=true function quadrupole_matrix!(i, b::BunchView, k1, s)
  v = b.v

  sgn = sign(k1)
  focus = k1 >= 0  # horizontally focusing if positive

  xp = v[i,PXI] / (1 + v[i,PZI])  # x'
  yp = v[i,PYI] / (1 + v[i,PZI])  # y'
  sqrtks = sqrt(abs(k1 / (1 + v[i,PZI]))) * s  # |κ|s
  cx = focus ? cos(sqrtks) : cosh(sqrtks)
  cy = focus ? cosh(sqrtks) : cos(sqrtks)
  sx = focus ? sincu(sqrtks) : sinhcu(sqrtks)
  sy = focus ? sinhcu(sqrtks) : sincu(sqrtks)

  v[i,PXI] = v[i,PXI] * cx - k1 * s * v[i,XI] * sx
  v[i,PYI] = v[i,PYI] * cy + k1 * s * v[i,YI] * sy
  v[i,ZI]  = (v[i,ZI] - (s / 4) * (  xp^2 * (1 + sx * cx)
                                    + yp^2 * (1 + sy * cy)
                                    + k1 / (1 + v[i,PZI])
                                        * ( v[i,XI]^2 * (1 - sx * cx)
                                          - v[i,YI]^2 * (1 - sy * cy) )
                                  )
                      + sgn * ( v[i,XI] * xp * (sqrtks * sx)^2
                              - v[i,YI] * yp * (sqrtks * sy)^2 ) / 2
              )
  v[i,XI]  = v[i,XI] * cx + xp * s * sx
  v[i,YI]  = v[i,YI] * cy + yp * s * sy
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
@makekernel fastgtpsa=true function quadrupole_kick!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, s)
  v = b.v
  rP0 = 1 + v[i,PZI]              # reduced total momentum,  P/P0 = 1 + δ
  sqrPt = v[i,PXI]^2 + v[i,PYI]^2   # (transverse momentum)^2, P⟂^2 = (Px^2 + Py^2) / P0^2
  Ps = sqrt(rP0^2 - sqrPt)          # longitudinal momentum,   Ps = √[(1 + δ)^2 - P⟂^2]
  v[i,XI] = v[i,XI] + s * v[i,PXI] / rP0 * sqrPt / (Ps * (rP0 + Ps))
  v[i,YI] = v[i,YI] + s * v[i,PYI] / rP0 * sqrPt / (Ps * (rP0 + Ps))
  v[i,ZI] = v[i,ZI] - s * ( rP0
                              * (sqrPt - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                                / ( beta_0 * sqrt(rP0^2 + tilde_m^2) * Ps
                                    * (beta_0 * sqrt(rP0^2 + tilde_m^2) + Ps)
                                  )
                            - sqrPt / (2 * rP0^2)
                          )
end


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
beta_0: β_0 = (βγ)_0 / √(γ_0^2)
brho_0: Bρ_0,  reference magnetic rigidity
hc: coordinate frame curvature
b0: dipole field strength
e1: entrance face angle (+ve angle <=> toward rbend)
e2: exit face angle (+ve angle <=> toward rbend)
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
sn: vector of skew multipole strengths scaled by Bρ0
L:  element arc length
"""
@makekernel fastgtpsa=true function bkb_multipole!(i, b::BunchView, beta_0, brho_0, hc, b0, e1, e2, mm, kn, sn, L)
  ExactTracking.exact_sbend!(   i, b, beta_0, brho_0, hc, b0, e1, e2, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L, sn * L)
  ExactTracking.exact_sbend!(   i, b, beta_0, brho_0, hc, b0, e1, e2, L / 2)
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
ks: solenoid strength
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
sn: vector of skew multipole strengths scaled by Bρ0
L:  element length
"""
@makekernel fastgtpsa=true function sks_multipole!(i, b::BunchView, ks, beta_0, gamsqr_0, tilde_m, mm, kn, sn, L)
  ExactTracking.exact_solenoid!(i, b, ks, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L, sn * L)
  ExactTracking.exact_solenoid!(i, b, ks, beta_0, gamsqr_0, tilde_m, L / 2)
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
starting with the dipole component. (For example, b[3] denotes
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
@makekernel fastgtpsa=true function dkd_multipole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, mm, kn, ks, L)
  ExactTracking.exact_drift!(   i, b, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L, ks * L)
  ExactTracking.exact_drift!(   i, b, beta_0, gamsqr_0, tilde_m, L / 2)
end


#
# ===============  S P I N  ===============
#
@inline function binom(m::Integer, x, y)
  """
  This function computes the real and imaginary parts of
  (x + i y)^m. One can use these components to compute,
  among other things, the multipole kick induced by two-
  dimensional multipole magnets.
  """
  if m == 0
    return (1.0, 0.0)
  end
  ar = x
  ai = y
  mm = m
  while mm > 1
    mm -= 1
    t  = x * ar - y * ai
    ai = y * ar + x * ai
    ar = t
  end
  return ar, ai
end 


@inline function field(ks, mm, kn, sn, x, y)
  """
  This function computes the magnetic field (Bx, By, Bz)/Brho at position (x,y)
  using the multipole coefficients kn and sn indexed by mm, i.e.,
  kn[i] is the normal coefficient of order mm[i]. Additionally, ks is the
  solenoid strength.
  """
  factorials = @SArray [
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5040,
    40320,
    362880,
    3628800,
    39916800,
    479001600,
    6227020800,
    87178291200,
    1307674368000,
    20922789888000,
    355687428096000,
    6402373705728000,
    121645100408832000,
    2432902008176640000,
    51090942171709440000
  ]
  Bx, By = 0, 0
  for i in 1:length(mm)
    ar, ai = binom(mm[i]-1, x, y)
    Bx += (kn[i] * ai + sn[i] * ar)/factorials[mm[i]]
    By += (kn[i] * ar - sn[i] * ai)/factorials[mm[i]]
  end
  return Bx, By, ks
end 


@inline function omega(i, b::BunchView, a, g, beta_0, gamma, ks, mm, kn, sn)
  """
  This function computes the spin-precession vector using the multipole 
  coefficients kn and sn indexed by mm, i.e., kn[i] is the normal 
  coefficient of order mm[i].
  """
  v = b.v

  rel_p = 1 + v[i,PZI]
  pl = sqrt(rel_p^2-v[i,PXI]^2-v[i,PYI]^2)
  Bx, By, Bz = field(ks, mm, kn, sn, v[i,XI], v[i,YI])
  Bdotbeta = 1/rel_p*(Bx*v[i,PXI] + By*v[i,PYI] + Bz*pl) 

  chi = 1 + a*gamma
  zeta = a*gamma^2/(1+gamma)*Bdotbeta/rel_p
  scale = -(1+g*v[i,XI])/pl
  
  ox = scale*(chi*Bx - zeta*v[i,PXI])
  oy = scale*(chi*By - zeta*v[i,PYI]) + g 
  oz = scale*(chi*Bz - zeta*pl)
  return @SArray [ox, oy, oz]
end


@inline function expq(v)
  """
  This function computes exp(i v⋅σ) as a quaternion, where σ is the 
  vector of Pauli matrices.
  """
  n = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
  c = cos(n)
  s = sincu(n)
  v2 = s .* v
  return @SArray [-c, v2[1], v2[2], v2[3]]
end


@inline function rotate_spin!(i, b::BunchView, a, g, beta_0, gamma, ks, mm, kn, sn, L)
  """
  This function rotates b.q according to the multipoles present.
  """
  q2 = b.q
  q1 = expq(-L/2 .* omega(i, b, a, g, beta_0, gamma, ks, mm, kn, sn))
  a1, b1, c1, d1 = q1[1], q1[2], q1[3], q1[4]
  a2, b2, c2, d2 = q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ]
  q2[i,Q0] = a1*a2 - b1*b2 - c1*c2 - d1*d2
  q2[i,QX] = a1*b2 + b1*a2 + c1*d2 - d1*c2
  q2[i,QY] = a1*c2 - b1*d2 + c1*a2 + d1*b2
  q2[i,QZ] = a1*d2 + b1*c2 - c1*b2 + d1*a2
end

@makekernel fastgtpsa=true function integrate_with_spin!(i, b::BunchView, ker, params, beta_0, gamsqr_0, tilde_m, g, ks, mm, kn, sn, L)
  a = 0.00115965218128                                # b.species.a
  betagamma = beta_0*sqrt(gamsqr_0) * (1 + b.v[i,PZI])
  gamma = sqrt(1+betagamma^2)
  if L != 0
    rotate_spin!(i, b, a, g, beta_0, gamma, ks, mm, kn, sn, L / 2)
    ker(i, b, params..., L)
    rotate_spin!(i, b, a, g, beta_0, gamma, ks, mm, kn, sn, L / 2)
  else
    rotate_spin!(i, b, a, g, beta_0, gamma, ks, mm, kn, sn, L)
    ker(i, b, params..., L)
  end
end

end