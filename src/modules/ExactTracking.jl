#=

Exact tracking methods

=#
# Define the Exact tracking method, and number of columns in the work matrix
# (equal to number of temporaries needed for a single particle)
struct Exact end

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..ReferenceFrameRotations, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView
using ..BeamTracking: C_LIGHT
const TRACKING_METHOD = Exact

# Update the reference energy of the canonical coordinates
# BUG: z and pz are not updated correctly
#=
@makekernel fastgtpsa=true function update_P0!(i, b, Brho_initial, Brho_final)
  @inbounds begin
    @FastGTPSA! v[i,PXI] = v[i,PXI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PYI] = v[i,PYI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PZI] = v[i,PZI] * Brho_initial / Brho_final
  end
  return v
end
=#

#
# ===============  E X A C T   D R I F T  ===============
#
"""
exact_drift!()

In the computation of z_final, we use the fact that
    1/√a - 1/√b == (b - a)/(√a √b (√a + √b))
to avoid the potential for severe cancellation when
a and b both have the form 1 + ε for different small
values of ε.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
L: element length
"""
@makekernel fastgtpsa=true function exact_drift!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, L)
  v = b.v
  P_s = sqrt((1 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2)) 
  v[i,XI]   = v[i,XI] + v[i,PXI] * L / P_s
  v[i,YI]   = v[i,YI] + v[i,PYI] * L / P_s
  # high-precision computation of z_final
  # vf.z = vi.z - (1 + δ) * L * (1 / Ps - 1 / (β0 * sqrt((1 + δ)^2 + tilde_m^2)))
  v[i,ZI]   = v[i,ZI] - ( (1 + v[i,PZI]) * L
                * ((v[i,PXI]^2 + v[i,PYI]^2) - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                / ( beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) * P_s
                    * (beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) + P_s)
                  )
              )
end # function exact_drift!()


"""
    multipole_kick!(i, b, ms, knl, ksl)

Track a beam of particles through a thin-lens multipole
having integrated normal and skew strengths listed in the
coefficient vectors knl and ksl respectively. The vector ms
lists the orders of the corresponding entries in knl and ksl.

The algorithm used in this function takes advantage of the
complex representation of the vector potential Az,
  - ``-Re{ sum_m (b_m + i a_m) (x + i y)^m / m }``,
and uses a Horner-like scheme (see Shachinger and Talman
[SSC-52]) to compute the transverse kicks induced by a pure
multipole magnet. This method supposedly has good numerical
properties, though I've not seen a proof of that claim.

DTA: Ordering matters!
DTA: Add thin dipole kick.

### Arguments
 - ms:  vector of m values for multipole coefficients
 - knl: vector of normal integrated multipole strengths
 - ksl: vector of skew integrated multipole strengths


     NB: Here the j-th component of knl (ksl) denotes the
       normal (skew) component of the multipole strength of
       order mm[j] (after scaling by the reference Bρ).
       For example, if mm[j] = 3, then knl[j] denotes the
       normal integrated sextupole strength scaled by Bρo.
"""
#@inline function multipole_kick!(i, b, mm, knl, ksl)
#  @inbounds begin #@FastGTPSA! begin
#  jm = length(mm)
#  m = mm[jm]
#  if m == 1 return v end
#  work[i,1] = ar = knl[jm] * v[i,XI] - ksl[jm] * v[i,YI]
#  work[i,2] = ai = knl[jm] * v[i,YI] + ksl[jm] * v[i,XI]
#  jm -= 1
#  while m > 2
#    m -= 1
#    work[i,3] = work[i,1] * v[i,XI] - work[i,2] * v[i,YI]
#    work[i,2] = work[i,1] * v[i,YI] + work[i,2] * v[i,XI]
#    work[i,1] = work[i,3]
#    if m == mm[jm]
#      work[i,1] += knl[jm] * v[i,XI] - ksl[jm] * v[i,YI]
#      work[i,2] += knl[jm] * v[i,YI] + ksl[jm] * v[i,XI]
#      jm -= 1
#    end
#  end
#  v[i,PXI] -= work[i,1]
#  v[i,PYI] += work[i,2]
#  end #end
#  return v
#end # function multipole_kick!()
#
@makekernel fastgtpsa=true function multipole_kick!(i, b::BunchView, ms, knl, ksl)
  v = b.v
  jm = length(ms)
  m  = ms[jm]
  knl_tot = knl[jm]
  ksl_tot = ksl[jm]
  jm -= 1
  while 2 <= m
    m -= 1
    tmp_knl_tot = (knl_tot * v[i,XI] - ksl_tot * v[i,YI]) / m
    ksl_tot     = (knl_tot * v[i,YI] + ksl_tot * v[i,XI]) / m
    knl_tot = tmp_knl_tot
    if 0 < jm && m == ms[jm]
      knl_tot += knl[jm]
      ksl_tot += ksl[jm]
      jm -= 1
    end
  end
  v[i,PXI] -= knl_tot
  v[i,PYI] += ksl_tot
end # function multipole_kick!()


#function binom(m::Integer, x, y)
#  """
#  This function computes the real and imaginary parts of
#  (x + i y)^m. One can use these components to compute,
#  among other things, the multipole kick induced by two-
#  dimensional multipole magnets.
#  """
#  if m == 0
#    return [ 1, 0.0 ]
#  end
#  ar = x
#  ai = y
#  mm = m
#  while mm > 1
#    mm -= 1
#    t  = x * ar - y * ai
#    ai = y * ar + x * ai
#    ar = t
#  end
#  return [ ar, ai ]
#end # function binom()


#
# ===============  E X A C T   S E C T O R   B E N D  ===============
#
"""
This function implements exact symplectic tracking through a
sector bend, derived using the Hamiltonian (25.9) given in the
Bmad manual. As a consequence of using that Hamiltonian, the
reference value of βγ must be that of a particle with the
design energy.  Should we wish to change that, we shall need
to carry both reference and design values.

Arguments
—————————
beta_0: β_0 = (βγ)_0 / √(γ_0^2)
brho_0: Bρ_0,  reference magnetic rigidity
hc: coordinate frame curvature
b0: magnet field strength
e1: entrance face angle (+ve angle <=> toward rbend)
e2: exit face angle (+ve angle <=> toward rbend)
Lr: element arc length
"""
@makekernel fastgtpsa=true function exact_sbend!(i, b::BunchView, beta_0, brho_0, hc, b0, e1, e2, Lr)
  v = b.v

  rho = brho_0 / b0
  ang = hc * Lr
  c1 = cos(ang)
  s1 = sin(ang)

  P_s     = sqrt((1 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2))  # P_s
  P_alpha = sqrt((1 + v[i,PZI])^2 - v[i,PYI]^2)                 # P_α
  s1phx = (1 + hc * v[i,XI]) / (hc * rho)                     # scaled (1 + h x)
  Pxpph = P_s - s1phx                                 # Px'/h
  ang_eff = ang + asin(v[i,PXI] / P_alpha) - asin((v[i,PXI] * c1 + Pxpph * s1) / P_alpha)  # α + φ1 - φ2
  # high-precision computation of x-final
  v[i,XI] = (v[i,XI] * c1 - Lr * sin(ang / 2) * sincu(ang / 2)
             + rho * (v[i,PXI] + ((v[i,PXI]^2 + (P_s + Pxpph) * s1phx) * s1 - 2v[i,PXI] * Pxpph * c1)
                             / (sqrt(P_alpha^2 - (v[i,PXI] * c1 + Pxpph * s1)^2) + P_s * c1)) * s1)
  v[i,PXI] = v[i,PXI] * c1 + Pxpph * s1
  v[i,YI] = v[i,YI] + rho * v[i,PYI] * ang_eff
  # high-precision computation of z-final
  v[i,ZI] = (v[i,ZI] - rho * (1 + v[i,PZI]) * ang_eff
               + (1 + v[i,PZI]) * Lr / (beta_0 * sqrt(1 / beta_0^2 + (2 + v[i,PZI]) * v[i,PZI])))
end # function exact_sbend!()


@makekernel fastgtpsa=true function exact_solenoid!(i, b::BunchView, ks, beta_0, gamsqr_0, tilde_m, L)
  v = b.v
  # Recurring variables 
  rel_p = 1 + v[i,PZI]    
  pr = sqrt(rel_p^2 - (v[i,PXI] + v[i,YI] * ks / 2)^2 - (v[i,PYI] - v[i,XI] * ks / 2)^2)
  s = sin(ks * L / pr)   
  cp = 1 + cos(ks * L / pr)                
  cm = 2 - cp
  # Temporaries
  x_0 = v[i,XI] 
  px_0 = v[i,PXI] 
  y_0 = v[i,YI]  
  # Update
  v[i,ZI]  -= rel_p * L *
                ((v[i,PXI] + v[i,YI] * ks / 2)^2 + (v[i,PYI] - v[i,XI] * ks / 2)^2 - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0) /
                ( beta_0 * sqrt(rel_p^2 + tilde_m^2) * pr * (beta_0 * sqrt(rel_p^2 + tilde_m^2) + pr) )
  v[i,XI] = cp * x_0 / 2 + s * (px_0 / ks + y_0 / 2) + cm * v[i,PYI] / ks
  v[i,PXI] = s * (v[i,PYI] / 2 - ks * x_0 / 4) + cp * px_0 / 2 - ks * cm * y_0 / 4
  v[i,YI] = s * (v[i,PYI] / ks - x_0 / 2) + cp * y_0 / 2 - cm * px_0 / ks
  v[i,PYI]  = ks * cm * x_0 / 4 - s * (px_0 / 2 + ks * y_0 / 4) + cp * v[i,PYI] / 2
end


@makekernel fastgtpsa=true function patch!(i, b::BunchView, tilde_m, dt, dx, dy, dz, winv::Union{StaticMatrix{3,3},Nothing}, L)
  # Temporary momentum [1+δp, ps_0]
  v = b.v
  rel_p = 1 + v[i,PZI]                                 
  ps_0 = sqrt(rel_p^2 - v[i,PXI]^2 - v[i,PYI]^2)  
  # Only apply rotations if needed
  if isnothing(winv)
    # No rotation case
    v[i,XI] -= dx
    v[i,YI] -= dy

    # Apply t_offset
    v[i,ZI] += rel_p/sqrt(rel_p^2+tilde_m^2)*C_LIGHT*dt

    # Drift to face
    v[i,XI]   += v[i,PXI] * dz / ps_0
    v[i,YI]   += v[i,PYI] * dz / ps_0
    v[i,ZI]   -=  dz * rel_p / ps_0 - L*rel_p*sqrt((1+tilde_m^2)/(rel_p^2+tilde_m^2))
  else
    # Translate position vector [x, y]
    x_0 = v[i,XI] - dx                                # x_0
    y_0 = v[i,YI] - dy                                # y_0

    # Temporary momentum vector [px, py]
    px_0 = v[i,PXI]                                    # px_0
    py_0 = v[i,PYI]                                    # py_0

    # Transform position vector [x - dx, y - dy, -dz]
    v[i,XI]   = winv[1,1]*x_0 + winv[1,2]*y_0 - winv[1,3]*dz
    v[i,YI]   = winv[2,1]*x_0 + winv[2,2]*y_0 - winv[2,3]*dz
    s_f = winv[3,1]*x_0 + winv[3,2]*y_0 - winv[3,3]*dz  # s_f

    # Transform momentum vector [px, py, ps]
    v[i,PXI]  = winv[1,1]*px_0 + winv[1,2]*py_0 + winv[1,3]*ps_0
    v[i,PYI]  = winv[2,1]*px_0 + winv[2,2]*py_0 + winv[2,3]*ps_0
    ps_f = winv[3,1]*px_0 + winv[3,2]*py_0 + winv[3,3]*ps_0 # ps_f

    # Apply t_offset
    v[i,ZI] += rel_p/sqrt(rel_p^2+tilde_m^2)*C_LIGHT*dt

    # Drift to face
    v[i,XI] -= s_f * v[i,PXI] / ps_f
    v[i,YI] -= s_f * v[i,PYI] / ps_f
    v[i,ZI] += s_f * rel_p / ps_f + L*rel_p*sqrt((1+tilde_m^2)/(rel_p^2+tilde_m^2))
  end
end


# Utility functions ============================================================

# Rotation matrix
"""
  w_matrix(x_rot, y_rot, z_rot)

Constructs a rotation matrix based on the given Bryan-Tait angles.

Bmad/SciBmad follows the MAD convention of applying z, x, y rotations in that order.
Furthermore, in ReferenceFrameRotations, the rotation angles are defined as negative
of the SciBmad rotation angles `x_rot`, `y_rot`, and `z_rot`.

The inverse matrix reverses the order of operations and their signs.


Arguments:
- `x_rot::Number`: Rotation angle around the x-axis.
- `y_rot::Number`: Rotation angle around the y-axis.
- `z_rot::Number`: Rotation angle around the z-axis.

Returns:
- `DCM{Float64}`: ReferenceFrameRotations.DCM (direct cosine matrix), rotation matrix.
"""
function w_matrix(x_rot, y_rot, z_rot)
  return ReferenceFrameRotations.angle_to_rot(-z_rot, -x_rot, -y_rot, :ZXY)
end

# Inverse rotation matrix
function w_inv_matrix(x_rot, y_rot, z_rot)
  return ReferenceFrameRotations.angle_to_rot(y_rot, x_rot, z_rot, :YXZ)
end

function drift_params(species::Species, Brho)
  beta_gamma_0 = BeamTracking.calc_beta_gammma(species, Brho)
  tilde_m = 1/beta_gamma_0
  gamsqr_0 = @FastGTPSA 1+beta_gamma_0^2
  beta_0 = @FastGTPSA beta_gamma_0/sqrt(gamsqr_0)
  return tilde_m, gamsqr_0, beta_0
end


end

