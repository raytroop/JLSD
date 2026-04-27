module BlkRX
using UnPack, DSP, Random, Interpolations
include("../util/Util_JLSD.jl"); using .Util_JLSD
include("./BlkElasticBuffer.jl"); using .BlkElasticBuffer

export clkgen_pi_itp_top!
export sample_itp_top!, sample_phi_top!, slicers_top!, sample_filter_top!
export cdr_top!, adpt_top!
# Re-export elastic-buffer helpers so callers (TB.jl, Widget.jl) only need
# `using .BlkRX` to get everything they need for the dynamic scheduling loop.
export eb_write!, eb_can_read_times, eb_interp


"""
    clkgen_pi_itp_top!(clkgen, eb; pi_code)

Generate one sub-block of RX sampling instants in absolute TX-grid
sample-time units.  The base position is the persistent RX cursor
`eb.t_rx`; the inter-symbol spacing is `param.osr_rx` (so a non-zero
`freq_offset_ppm` naturally compresses or stretches the per-symbol
stride).

CDR phase offset `Φ0`, static phase skew `Φskew`, and random jitter
`Φrj` are all defined in **RX UI** and converted to TX-grid samples
using `osr_rx` (not the nominal `osr`).  This keeps the physical
meaning consistent in the presence of a frequency offset.

Notes:
  * `eb.t_rx` is **not** advanced here — the dynamic scheduling loop
    (in TB.jl / Widget.jl) is responsible for advancing it by
    `subblk_size * osr_rx` only after readiness is confirmed.
  * `clkgen.Φo` history is **not** appended here either — the caller
    appends it only after `eb_can_read_times` confirms the candidate
    times can be served.  This guarantees that random jitter is *not*
    re-sampled between readiness check and actual read.
"""
function clkgen_pi_itp_top!(clkgen, eb; pi_code)
    @unpack tui, osr_rx, subblk_size = clkgen.param
    @unpack nphases, rj, skews = clkgen
    @unpack pi_code_prev, pi_wrap_ui, pi_wrap_ui_Δcode = clkgen
    @unpack pi_nonlin_lut, pi_ui_cover, pi_codes_per_ui = clkgen

    Δpi_code = pi_code - pi_code_prev
    if abs(Δpi_code) > pi_wrap_ui_Δcode
        pi_wrap_ui -= sign(Δpi_code)*pi_ui_cover
    end

    # Φ0 lives in RX-UI units; scaling by nominal `osr` would inject a
    # constant timing error proportional to the frequency offset.
    Φ0    = osr_rx*(pi_wrap_ui + (pi_code + pi_nonlin_lut[pi_code+1])/pi_codes_per_ui)
    Φskew = kron(ones(Int(subblk_size/nphases)), skews/tui*osr_rx)
    Φrj   = (rj/tui)*osr_rx*randn(subblk_size)

    t_rx = eb.t_rx
    @inbounds for j = 0:subblk_size-1
        clkgen.Φo_subblk[j+1] = Φ0 + t_rx + j*osr_rx + Φskew[j+1] + Φrj[j+1]
    end

    clkgen.pi_code_prev = pi_code
    clkgen.pi_wrap_ui   = pi_wrap_ui
    # NOTE: append!(clkgen.Φo, ...) is intentionally NOT done here.
    # See the docstring above.
end



"""
    sample_filter_top!(splr, Vi)

Apply the sampler's bandwidth-limiting impulse response to the channel
waveform `Vi`, updating `splr.Vo` in-place.  Inter-block memory is
preserved through `splr.Vo_mem`.  This is the new entry point used by
the elastic-buffer pipeline: the filter output is then pushed into the
buffer via `eb_write!`.
"""
function sample_filter_top!(splr, Vi)
    @unpack dt = splr.param
    @unpack ir, Vo_conv, Vo_mem = splr
    u_conv!(Vo_conv, Vi, ir, Vi_mem=Vo_mem, gain=dt)
    return nothing
end


"""
    sample_itp_top!(splr, Vi)

Legacy single-block path: apply the RX bandwidth filter then build a
block-local linear interpolant `splr.itp_Vext`.  Retained for backward
compatibility with code paths that still use `sample_phi_top!(splr, Φi)`
without an elastic buffer; not used in the new dynamic-scheduling flow.
"""
function sample_itp_top!(splr, Vi)
    @unpack ir, Vo, Vo_mem = splr
    @unpack V_prev_nui, Vext, tt_Vext = splr

    sample_filter_top!(splr, Vi)

    Vext[eachindex(V_prev_nui)] .= V_prev_nui
    Vext[lastindex(V_prev_nui)+1:end] .= Vo
    splr.itp_Vext = linear_interpolation(tt_Vext, Vext)
end


"""
    sample_phi_top!(splr, eb, Φi)

Sample the elastic buffer `eb` at the absolute TX-grid float times in
`Φi`, writing into `splr.So_subblk` and appending to `splr.So`.  The
caller must have already verified `eb_can_read_times(eb, Φi)`.
"""
function sample_phi_top!(splr, eb, Φi)
    @inbounds for j = eachindex(Φi)
        splr.So_subblk[j] = eb_interp(eb, Φi[j])
    end
    append!(splr.So, splr.So_subblk)
    return nothing
end

function sample_phi_top!(splr, Φi)
    @unpack itp_Vext = splr
    splr.So_subblk .= itp_Vext.(Φi)
    append!(splr.So, splr.So_subblk)
    return nothing
end

function slicers_top!(slc, Si; ref_code)
    @unpack nphases, noise_rms, dac_min, dac_lsb = slc
    @unpack ofsts, N_per_phi = slc

    ref_lvl = [(dac_min .+ dac_lsb * ref_code[n]) for n in 1:nphases]

    for n = eachindex(Si)
        phi_idx = (n-1)%nphases + 1
        nslc = N_per_phi[phi_idx]
        if nslc != 0
            slc.So[n] .=  (ref_lvl[phi_idx]
                            + ofsts[phi_idx]
                            + (noise_rms * randn(nslc))
                            ) .< Si[n]

        end
    end
end


function cdr_top!(cdr, Sd, Se)
    @unpack Neslc_per_phi, Sd_prev = cdr
    @unpack eslc_nvec, filt_patterns, kp, ki = cdr
    @unpack pd_accum, ki_accum, pd_gain, pi_res = cdr

    pi_bnd = 2^pi_res
    Sd_val = [Sd_prev; [sum(dvec) for dvec in Sd]]

    for n = findall(eslc_nvec.!=0)
        if Sd_val[n:n+2] in filt_patterns
            vote = sign(Se[n][1].-0.5)*sign(Sd_val[n]-Sd_val[n+2])
            ki_accum += ki*vote
            pd_accum += pd_gain*(kp*vote + ki_accum)
        end
    end

    cdr.pd_accum = (pd_accum < 0) ? pi_bnd + pd_accum : (pd_accum >= pi_bnd) ? pd_accum - pi_bnd : pd_accum
    cdr.pi_code = Int(floor(cdr.pd_accum))


    cdr.Sd_prev = Sd_val[end]
    cdr.ki_accum = ki_accum

end

function adpt_top!(adpt, Sd, Se)
    @unpack Neslc_per_phi, Sd_prev = adpt
    @unpack eslc_nvec, eslc_filt_patterns, eslc_ref_max, mu_eslc = adpt

    Sd_val = [Sd_prev; [sum(dvec) for dvec in Sd]]

    ref_accum = adpt.eslc_ref_accum

    for n = findall(eslc_nvec.!=0)
        ref_accum +=  (Sd_val[n:n+2] in eslc_filt_patterns) ?
                        mu_eslc*sign(Se[n][1].-0.5) : 0
    end
    adpt.eslc_ref_accum =   ref_accum < 0 ? 0 :
                            ref_accum > eslc_ref_max ? eslc_ref_max :
                            ref_accum
    adpt.eslc_ref_code = floor(adpt.eslc_ref_accum)
    adpt.eslc_ref_vec = [adpt.eslc_ref_code*ones(Int,n) for n in Neslc_per_phi]

    adpt.Sd_prev = Sd_val[end]
end

end
