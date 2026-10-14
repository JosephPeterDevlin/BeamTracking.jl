module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays, KernelAbstractions, AtomicAndPhysicalConstants, LinearAlgebra, Statistics, FFTW
using Beamlines: isactive, deval, unsafe_getparams, isnullspecies
using BeamTracking: R_to_E, R_to_beta_gamma, R_to_gamma, R_to_pc, R_to_v,
                    beta_gamma_to_v, E_to_R, E_to_v,
                    @makekernel, Coords, make_kernel_call, KernelCall, KernelChain, push, TimeDependentParam, RefState, 
<<<<<<< HEAD
                    launch!, AbstractYoshida, rot_quaternion, inv_rot_quaternion, atan2, 
                    get_N_particle, mean_and_cov, ibs_integrals, remake, push_transforms_in, push_transforms_out,
                    ZI
=======
                    launch!, AbstractSymplectic, rot_quaternion, inv_rot_quaternion, atan2, exact_wedge!,
                    get_N_particle, mean_and_cov, ibs_integrals, remake, push_transforms_in, push_transforms_out
>>>>>>> 65e343f9afda8a208fdfaeb487d615b0ca5a408e
                    
import BeamTracking: track!

function track!(
  bunch::Bunch, 
  ele::LineElement;
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  rf_on::Bool=true,
  _p_over_q_ref=nothing,
  context=(haskey(getfield(ele, :pdict), BeamlineParams) ? ele.beamline.context : Beamlines.NULL_CONTEXT),
  batch_start=1,
  kwargs...
)
  if isnothing(_p_over_q_ref)
    p_over_q_ref = ele.p_over_q_ref
  else
    p_over_q_ref = _p_over_q_ref
  end
  coords = bunch.coords
<<<<<<< HEAD
  @noinline _track!(coords, bunch, ele, context, p_over_q_ref, ele.tracking_method, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on; kwargs...)
=======
  @noinline _track!(coords, bunch, ele, context, p_over_q_ref, ele.tracking_method, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on, batch_start; kwargs...)
>>>>>>> 65e343f9afda8a208fdfaeb487d615b0ca5a408e
  return bunch
end

function track!(
  bunch::Bunch, 
  bl::Beamline; 
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  rf_on::Bool=true,
<<<<<<< HEAD
=======
  batch_start=1,
>>>>>>> 65e343f9afda8a208fdfaeb487d615b0ca5a408e
  kwargs...
)
  if length(bl.line) == 0
    return bunch
  end
  __, p_over_q_ref = check_bl_bunch!(bunch, bl)
  context = bl.context
  
  for ele in bl.line
<<<<<<< HEAD
    track!(bunch, ele; context, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on, kwargs...)
=======
    track!(bunch, ele; context, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on, batch_start, kwargs...)
>>>>>>> 65e343f9afda8a208fdfaeb487d615b0ca5a408e
  end

  return bunch
end

include("utils_bl.jl")
include("unpack_bl.jl")
include("scibmadstandard_bl.jl")
include("exact_bl.jl")
include("symplectic_bl.jl")
include("sagan_cavity_bl.jl")
include("general_bl.jl")
include("rungekutta.jl")

end