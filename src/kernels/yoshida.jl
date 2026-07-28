#
# ===============  I N T E G R A T O R S  ===============
#

@inline function order_two_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, L) where {fringe_in,fringe_out}
  @inbounds begin
    if !isnothing(edge_params) && fringe_in
      a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
      linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e1, 1)
    end
    s = 0
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    for step in 1:(n_steps-1)
      ker(i, coords, s, params..., ds_step)
      s += ds_step
      if !isnothing(photon_params)
        stochastic_radiation!(i, coords, s, photon_params..., ds_step)
      end
      dt_ref = compute_dt_ref(s, ker, params)
      execute_callbacks(i, coords, s, dt_ref)
    end
    ker(i, coords, s, params..., ds_step)
    s += ds_step
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    if !isnothing(edge_params) && fringe_out
      a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
      linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e2, -1)
    end
  end
  return nothing
end


@inline function order_four_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, L) where {fringe_in,fringe_out}
  @inbounds begin
    w0 = -1.7024143839193153215916254339390434324741363525390625*ds_step
    w1 =  1.3512071919596577718181151794851757586002349853515625*ds_step

    g1 =  0.414490771794375711944979912004782818257808685302734375*ds_step
    g0 = -0.6579630871775028477799196480191312730312347412109375*ds_step
    if !isnothing(edge_params) && fringe_in
      a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
      linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e1, 1)
    end
    s = 0
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    for step in 1:(n_steps-1)
      #=
      ker(i, coords, s, params..., w1)
      s += w1
      ker(i, coords, s, params..., w0)
      s += w0
      ker(i, coords, s, params..., w1)
      s += w1
      =#
      ker(i, coords, s, params..., g1)
      s += g1
      ker(i, coords, s, params..., g1)
      s += g1
      ker(i, coords, s, params..., g0)
      s += g0
      ker(i, coords, s, params..., g1)
      s += g1
      ker(i, coords, s, params..., g1)
      s += g1
      
      if !isnothing(photon_params)
        stochastic_radiation!(i, coords, s, photon_params..., ds_step)
      end
      dt_ref = compute_dt_ref(s, ker, params)
      execute_callbacks(i, coords, s, dt_ref)
    end
    #=
    ker(i, coords, s, params..., w1)
    s += w1
    ker(i, coords, s, params..., w0)
    s += w0
    ker(i, coords, s, params..., w1)
    s += w1
    =#
    ker(i, coords, s, params..., g1)
    s += g1
    ker(i, coords, s, params..., g1)
    s += g1
    ker(i, coords, s, params..., g0)
    s += g0
    ker(i, coords, s, params..., g1)
    s += g1
    ker(i, coords, s, params..., g1)
    s += g1
    
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    if !isnothing(edge_params) && fringe_out
      a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
      linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e2, -1)
    end
  end
  return nothing
end


@inline function order_six_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, L) where {fringe_in,fringe_out}
  @inbounds begin
    w0 =  1.315186320683911169737712043570355*ds_step
    w1 = -1.17767998417887100694641568096432*ds_step
    w2 =  0.235573213359358133684793182978535*ds_step
    w3 =  0.784513610477557263819497633866351*ds_step

    g1 =  0.13861930854051695245808013042625*ds_step
    g2 =  0.13346562851074760407046858832209*ds_step
    g3 =  0.13070531011449225190542755785015*ds_step
    g4 =  0.12961893756907034772505366537091*ds_step
    g5 = -0.35000324893920896516170830911323*ds_step
    g6 =  0.11805530653002387170273438954049*ds_step
    g7 =  0.39907751534871587459988795520665*ds_step
    if !isnothing(edge_params)  && fringe_in
      a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
      linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e1, 1)
    end
    s = zero(w0)
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    for step in 1:(n_steps-1)
      #=
      ker(i, coords, s, params..., w3)
      s += w3
      ker(i, coords, s, params..., w2)
      s += w2
      ker(i, coords, s, params..., w1)
      s += w1
      ker(i, coords, s, params..., w0)
      s += w0
      ker(i, coords, s, params..., w1)
      s += w1
      ker(i, coords, s, params..., w2)
      s += w2
      ker(i, coords, s, params..., w3)
      s += w3
      =#
      ker(i, coords, s, params..., g1)
      s += g1
      ker(i, coords, s, params..., g2)
      s += g2
      ker(i, coords, s, params..., g3)
      s += g3
      ker(i, coords, s, params..., g4)
      s += g4
      ker(i, coords, s, params..., g5)
      s += g5
      ker(i, coords, s, params..., g6)
      s += g6
      ker(i, coords, s, params..., g7)
      s += g7
      ker(i, coords, s, params..., g6)
      s += g6
      ker(i, coords, s, params..., g5)
      s += g5
      ker(i, coords, s, params..., g4)
      s += g4
      ker(i, coords, s, params..., g3)
      s += g3
      ker(i, coords, s, params..., g2)
      s += g2
      ker(i, coords, s, params..., g1)
      s += g1
      
      if !isnothing(photon_params)
        stochastic_radiation!(i, coords, s, photon_params..., ds_step)
      end
      dt_ref = compute_dt_ref(s, ker, params)
      execute_callbacks(i, coords, s, dt_ref)
    end
    #=
    ker(i, coords, s, params..., w3)
    s += w3
    ker(i, coords, s, params..., w2)
    s += w2
    ker(i, coords, s, params..., w1)
    s += w1
    ker(i, coords, s, params..., w0)
    s += w0
    ker(i, coords, s, params..., w1)
    s += w1
    ker(i, coords, s, params..., w2)
    s += w2
    ker(i, coords, s, params..., w3)
    s += w3
    =#
    ker(i, coords, s, params..., g1)
    s += g1
    ker(i, coords, s, params..., g2)
    s += g2
    ker(i, coords, s, params..., g3)
    s += g3
    ker(i, coords, s, params..., g4)
    s += g4
    ker(i, coords, s, params..., g5)
    s += g5
    ker(i, coords, s, params..., g6)
    s += g6
    ker(i, coords, s, params..., g7)
    s += g7
    ker(i, coords, s, params..., g6)
    s += g6
    ker(i, coords, s, params..., g5)
    s += g5
    ker(i, coords, s, params..., g4)
    s += g4
    ker(i, coords, s, params..., g3)
    s += g3
    ker(i, coords, s, params..., g2)
    s += g2
    ker(i, coords, s, params..., g1)
    s += g1
    
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    if !isnothing(edge_params) && fringe_out
      a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
      linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e2, -1)
    end
  end
  return nothing
end


@inline function order_eight_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, L) where {fringe_in,fringe_out}
  @inbounds begin
    w0 =  1.7084530707869978*ds_step
    w1 =  0.102799849391985*ds_step
    w2 = -1.96061023297549*ds_step
    w3 =  1.93813913762276*ds_step
    w4 = -0.158240635368243*ds_step
    w5 = -1.44485223686048*ds_step
    w6 =  0.253693336566229*ds_step
    w7 =  0.914844246229740*ds_step

    g1 =  0.10647728984550031823931967854896*ds_step
    g2 =  0.10837408645835726397433410591546*ds_step
    g3 =  0.35337821052654342419534541324080*ds_step
    g4 = -0.23341414023165082198780281128319*ds_step
    g5 = -0.24445266791528841269462171413216*ds_step
    g6 =  0.11317848435755633314700952515599*ds_step
    g7 =  0.11892905625000350062692972283951*ds_step
    g8 =  0.12603912321825988140305670268365*ds_step
    g9 =  0.12581718736176041804392391641587*ds_step
    g10 = 0.11699135019217642180722881433533*ds_step
    g11 = -0.38263596012643665350944670744040*ds_step
    if !isnothing(edge_params) && fringe_in
      a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
      linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e1, 1)
    end
    s = 0
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    for step in 1:(n_steps-1)
      #=
      ker(i, coords, s, params..., w7)
      s += w7
      ker(i, coords, s, params..., w6)
      s += w6
      ker(i, coords, s, params..., w5)
      s += w5
      ker(i, coords, s, params..., w4)
      s += w4
      ker(i, coords, s, params..., w3)
      s += w3
      ker(i, coords, s, params..., w2)
      s += w2
      ker(i, coords, s, params..., w1)
      s += w1
      ker(i, coords, s, params..., w0)
      s += w0
      ker(i, coords, s, params..., w1) 
      s += w1
      ker(i, coords, s, params..., w2)
      s += w2
      ker(i, coords, s, params..., w3)
      s += w3
      ker(i, coords, s, params..., w4)
      s += w4
      ker(i, coords, s, params..., w5)
      s += w5
      ker(i, coords, s, params..., w6)
      s += w6
      ker(i, coords, s, params..., w7)
      s += w7
      =#
      ker(i, coords, s, params..., g1)
      s += g1
      ker(i, coords, s, params..., g2)
      s += g2
      ker(i, coords, s, params..., g3)
      s += g3
      ker(i, coords, s, params..., g4)
      s += g4
      ker(i, coords, s, params..., g5)
      s += g5
      ker(i, coords, s, params..., g6)
      s += g6
      ker(i, coords, s, params..., g7)
      s += g7
      ker(i, coords, s, params..., g8)
      s += g8
      ker(i, coords, s, params..., g9)
      s += g9
      ker(i, coords, s, params..., g10)
      s += g10
      ker(i, coords, s, params..., g11)
      s += g11
      ker(i, coords, s, params..., g10)
      s += g10
      ker(i, coords, s, params..., g9)
      s += g9
      ker(i, coords, s, params..., g8)
      s += g8
      ker(i, coords, s, params..., g7)
      s += g7
      ker(i, coords, s, params..., g6)
      s += g6
      ker(i, coords, s, params..., g5)
      s += g5
      ker(i, coords, s, params..., g4)
      s += g4
      ker(i, coords, s, params..., g3)
      s += g3
      ker(i, coords, s, params..., g2)
      s += g2
      ker(i, coords, s, params..., g1)
      s += g1
      
      if !isnothing(photon_params)
        stochastic_radiation!(i, coords, s, photon_params..., ds_step)
      end
      dt_ref = compute_dt_ref(s, ker, params)
      execute_callbacks(i, coords, s, dt_ref)
    end
    #=
    ker(i, coords, s, params..., w7)
    s += w7
    ker(i, coords, s, params..., w6)
    s += w6
    ker(i, coords, s, params..., w5)
    s += w5
    ker(i, coords, s, params..., w4)
    s += w4
    ker(i, coords, s, params..., w3)
    s += w3
    ker(i, coords, s, params..., w2)
    s += w2
    ker(i, coords, s, params..., w1)
    s += w1
    ker(i, coords, s, params..., w0)
    s += w0
    ker(i, coords, s, params..., w1) 
    s += w1
    ker(i, coords, s, params..., w2)
    s += w2
    ker(i, coords, s, params..., w3)
    s += w3
    ker(i, coords, s, params..., w4)
    s += w4
    ker(i, coords, s, params..., w5)
    s += w5
    ker(i, coords, s, params..., w6)
    s += w6
    ker(i, coords, s, params..., w7)
    s += w7
    =#
    ker(i, coords, s, params..., g1)
    s += g1
    ker(i, coords, s, params..., g2)
    s += g2
    ker(i, coords, s, params..., g3)
    s += g3
    ker(i, coords, s, params..., g4)
    s += g4
    ker(i, coords, s, params..., g5)
    s += g5
    ker(i, coords, s, params..., g6)
    s += g6
    ker(i, coords, s, params..., g7)
    s += g7
    ker(i, coords, s, params..., g8)
    s += g8
    ker(i, coords, s, params..., g9)
    s += g9
    ker(i, coords, s, params..., g10)
    s += g10
    ker(i, coords, s, params..., g11)
    s += g11
    ker(i, coords, s, params..., g10)
    s += g10
    ker(i, coords, s, params..., g9)
    s += g9
    ker(i, coords, s, params..., g8)
    s += g8
    ker(i, coords, s, params..., g7)
    s += g7
    ker(i, coords, s, params..., g6)
    s += g6
    ker(i, coords, s, params..., g5)
    s += g5
    ker(i, coords, s, params..., g4)
    s += g4
    ker(i, coords, s, params..., g3)
    s += g3
    ker(i, coords, s, params..., g2)
    s += g2
    ker(i, coords, s, params..., g1)
    s += g1
    
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    if !isnothing(edge_params)  && fringe_out
      a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
      linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e2, -1)
    end
  end
  return nothing
end