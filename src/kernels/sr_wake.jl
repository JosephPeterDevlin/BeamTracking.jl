function bin!(i, coords::Coords, tilde_m, tmin, tmax, dt, ::Val{ix}) where {ix}
  @inbounds begin
    v = coords.v
    alive = coords.state[i] == STATE_ALIVE

    rel_p = 1 + v[i,PZI]
    beta = rel_p/sqrt(rel_p*rel_p + tilde_m*tilde_m)
    t = -v[i,ZI]/(beta*C_LIGHT)
    k = floor((t-tmin)/dt)
    tk = tmin + k*dt
    k = Int(k)
    w = (t-tk)/dt

    if ix == 0
      weight_factor = 1
    elseif ix == 1
      weight_factor = v[i,XI]
    elseif ix == 2
      weight_factor = v[i,YI]
    else
      error("Help! I shouldn't be here.")
    end

    if !isnothing(coords.weight)
      if coords.weight isa Number
        weight_factor *= coords.weight
      else
        weight_factor *= coords.weight[i]
      end
    end

    if ix == 0
      arr = coords.longitudinal_density
    elseif ix == 1
      arr = coords.dipole_density_x
    elseif ix == 2
      arr = coords.dipole_density_y
    else
      error("Help! I shouldn't be here.")
    end

    arr[k]   += vifelse(alive, weight_factor*(1-w), 0)
    arr[k+1] += vifelse(alive, weight_factor*w,     0)
  end
  return nothing
end


@inline bin_z!(i, coords::Coords, tilde_m, tmin, tmax, dt) = bin!(i, coords, tilde_m, tmin, tmax, dt, Val{0}())
@inline bin_x!(i, coords::Coords, tilde_m, tmin, tmax, dt) = bin!(i, coords, tilde_m, tmin, tmax, dt, Val{1}())
@inline bin_y!(i, coords::Coords, tilde_m, tmin, tmax, dt) = bin!(i, coords, tilde_m, tmin, tmax, dt, Val{2}())


@makekernel fastgtpsa=true function srwake_long!(i, coords::Coords, tilde_m, tmin, dt, voltage)
  v = coords.v
  alive = coords.state[i] == STATE_ALIVE

  rel_p = 1 + v[i,PZI]
  beta = rel_p/sqrt(rel_p*rel_p + tilde_m*tilde_m)
  t = -v[i,ZI]/(beta*C_LIGHT)
  k = floor((t-tmin)/dt)
  tk = tmin + k*dt
  k = Int(k)
  w = (t-tk)/dt
  
  new_pz = v[i,PZI] + beta*((1-w)*voltage[k] + w*voltage[k+1])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end