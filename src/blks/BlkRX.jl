module BlkRX
using UnPack, DSP, Random, Interpolations
include("../util/Util_JLSD.jl"); using .Util_JLSD

export clkgen_pi_itp_top!
export sample_itp_top!, sample_phi_top!, sample_filter_top!, slicers_top!
export cdr_top!, adpt_top!
export rxw_extend!, rxw_covers, rxw_interp, rxw_phase_margin
export clkgen_commit_subblk!


# ── RxWindow helpers ─────────────────────────────────────────────────────────

"""
    rxw_extend!(rxw, waveform; phase_margin = 0.0)

Append `waveform` (one TX block of filtered samples on the TX grid) to
the RX window, advancing the write frontier `t_max`.

1. Reclaims storage already consumed by RX, retaining `phase_margin` samples
   below `t_rx` so negative PI/skew/RJ excursions are still serviceable.
2. Hard-errors on overrun — silent overwrite would corrupt the continuous
   analog signal with a step discontinuity and invalidate downstream
   CDR / BER results.

Underrun (RX cursor outran TX writer) is not handled here; the inner-loop
`rxw_covers` check breaks the scheduler benignly when that happens.
"""
function rxw_extend!(rxw, waveform::AbstractVector{Float64}; phase_margin::Real = 0.0)
    reclaim_at = floor(Int, rxw.t_rx - max(Float64(phase_margin), 0.0))
    rxw.t_min = max(rxw.t_min, reclaim_at)

    n = length(waveform)
    if (rxw.t_max + n) - rxw.t_min > rxw.capacity
        error("""
        RxWindow overrun: TX about to overwrite analog samples that RX has
        not yet consumed.  This would inject a step discontinuity into the
        sampled waveform and invalidate downstream BER / CDR results.

        t_min = $(rxw.t_min)   t_max = $(rxw.t_max)   t_rx = $(rxw.t_rx)
        next write of $n samples would push t_max to $(rxw.t_max + n)
        capacity = $(rxw.capacity)

        Resolve by reducing |freq_offset_ppm|, ensuring the RX scheduler
        is greedy, or constructing RxWindow with a larger capacity.
        """)
    end

    cap = rxw.capacity
    t_start = rxw.t_max
    @inbounds for i = 0:n-1
        rxw.buf[(t_start + i) % cap + 1] = waveform[i + 1]
    end
    rxw.t_max += n
    return nothing
end

"""
    rxw_phase_margin(clkgen; rj_sigma = 8.0) -> Float64

Conservative low-side timing margin, in TX-grid samples, retained by
`rxw_extend!` below `rxw.t_rx`.  It covers a PI wrap/phase excursion,
deterministic skew, and a finite Gaussian RJ envelope.
"""
function rxw_phase_margin(clkgen; rj_sigma::Real = 8.0)
    @unpack tui, osr, osr_rx = clkgen.param
    @unpack pi_ui_cover, skews, rj = clkgen

    pi_margin = pi_ui_cover * osr_rx
    skew_margin = isempty(skews) ? 0.0 : maximum(abs.(skews)) / tui * osr
    rj_margin = max(Float64(rj_sigma), 0.0) * abs(rj) / tui * osr
    return pi_margin + skew_margin + rj_margin
end

"""
    rxw_covers(rxw, Φi) -> Bool

Return `true` when the RX window contains enough data for linear
interpolation at every time in `Φi`.  Has no side effects.

Checks
- `floor(min(Φi))     >= t_min`
- `floor(max(Φi)) + 1 <  t_max`         (the `+1` is the k+1 interp tap)

`Φi` already encodes Φ0 (PI/CDR phase), Φskew (clock skew), and Φrj
(random jitter), so the check uses the *actual* upcoming sample times
rather than nominal-spacing predictions.
"""
function rxw_covers(rxw, Φi)
    t_need_min = floor(Int, minimum(Φi))
    t_need_max = floor(Int, maximum(Φi)) + 1
    return t_need_min >= rxw.t_min && t_need_max < rxw.t_max
end

"""
    rxw_interp(rxw, t) -> Float64

Linear interpolation of the RX window at Float64 TX-grid time `t`.
Caller must ensure `floor(t) >= t_min` and `floor(t) + 1 < t_max`
(use `rxw_covers` to verify), because this kernel always reads both
the `k0` and `k0 + 1` taps.
"""
@inline function rxw_interp(rxw, t::Float64)
    k0   = floor(Int, t)
    frac = t - k0
    cap  = rxw.capacity
    v0   = rxw.buf[k0       % cap + 1]
    v1   = rxw.buf[(k0 + 1) % cap + 1]
    return muladd(frac, v1 - v0, v0)
end


# ── RX block functions ───────────────────────────────────────────────────────


"""
    clkgen_pi_itp_top!(clkgen, rxw; pi_code)

Populate `clkgen.Φo_subblk` with the absolute TX-grid sampling times for
the next RX sub-block, anchored at the RX read cursor `rxw.t_rx`.

- Inter-sample spacing is `osr_rx` (TX-grid samples per RX UI).
- CDR/PI phase offset `Φ0` is scaled by `osr_rx` because `pi_ui_cover` is
  expressed in RX-clock UI.
- Skew and RJ are absolute timing deviations, converted via `osr/t_ui`.

This function prepares a pending candidate only.  It intentionally does not
advance `rxw.t_rx` or update `clkgen.pi_code_prev`; those are committed by
`clkgen_commit_subblk!` only after `rxw_covers` succeeds and the sub-block is
sampled.  If the candidate is not ready yet, the same `Φo_subblk` is retried
after the next TX block instead of regenerating RJ.
"""
function clkgen_pi_itp_top!(clkgen, rxw; pi_code)
    clkgen.Φo_subblk_valid && return nothing

    @unpack tui, osr, osr_rx, subblk_size = clkgen.param
    @unpack nphases, rj, skews = clkgen
    @unpack pi_code_prev, pi_wrap_ui_Δcode = clkgen
    @unpack pi_nonlin_lut, pi_ui_cover, pi_codes_per_ui = clkgen

    Δpi_code = pi_code - pi_code_prev
    if abs(Δpi_code) > pi_wrap_ui_Δcode
        # When pi_code wraps (e.g. 255→0 for ppm<0), absorb the pi_ui_cover
        # discontinuity into the candidate cursor so the absolute sampling
        # phase remains continuous.  This replaces the legacy pi_wrap_ui term
        # in Φ0 — otherwise, for non-zero ppm the CDR's cumulative drift
        # compensation would double-count with the dynamic scheduler.
        t_rx = rxw.t_rx - sign(Δpi_code) * pi_ui_cover * osr_rx
    else
        t_rx = rxw.t_rx
    end

    # Φ0 is now bounded in [0, pi_ui_cover · osr_rx) — a fine within-PI phase
    # correction, not a cumulative drift term.
    Φ0    = osr_rx*(pi_code + pi_nonlin_lut[pi_code+1])/pi_codes_per_ui
    Φskew = kron(ones(Int(subblk_size/nphases)), skews/tui*osr)
    Φrj   = rj/tui*osr*randn(subblk_size)

    @inbounds for j = 0:subblk_size-1
        clkgen.Φo_subblk[j+1] = Φ0 + t_rx + j*osr_rx + Φskew[j+1] + Φrj[j+1]
    end

    clkgen.t_rx_subblk = t_rx
    clkgen.t_rx_next = t_rx + subblk_size * osr_rx
    clkgen.pi_code_subblk = pi_code
    clkgen.Φo_subblk_valid = true
    return nothing
end

"""
    clkgen_commit_subblk!(clkgen, rxw)

Commit the pending RX sample-time candidate after the sub-block has been
sampled.  This is the only place that advances `rxw.t_rx` and
`clkgen.pi_code_prev` for frequency-offset scheduling.
"""
function clkgen_commit_subblk!(clkgen, rxw)
    clkgen.Φo_subblk_valid || error("clkgen_commit_subblk!: no pending RX sub-block")

    rxw.t_rx = clkgen.t_rx_next
    clkgen.pi_code_prev = clkgen.pi_code_subblk
    clkgen.Φo_subblk_valid = false
    return nothing
end


"""
    sample_filter_top!(splr, Vi)

Convolve the channel waveform `Vi` with the sampler bandwidth-limit
impulse response, writing the result into `splr.Vo` (block-stitched via
`splr.Vo_mem`).  This is the analog-domain pre-sampling filter; the
result is then handed to the RxWindow.
"""
function sample_filter_top!(splr, Vi)
    @unpack dt = splr.param
    @unpack ir, Vo_conv, Vo_mem = splr
    u_conv!(Vo_conv, Vi, ir, Vi_mem=Vo_mem, gain=dt)
    return nothing
end


"""
    sample_itp_top!(splr, Vi)

Legacy path: filter the channel waveform and build a block-local
interpolation object.  Retained for compatibility; the freq-offset path
uses `sample_filter_top!` + `rxw_extend!` + `rxw_interp` instead.
"""
function sample_itp_top!(splr, Vi)
    @unpack osr, dt, blk_size_osr = splr.param
    @unpack ir, Vo_conv, Vo, Vo_mem = splr
    @unpack prev_nui, V_prev_nui, Vext, tt_Vext = splr

    sample_filter_top!(splr, Vi)

    Vext[eachindex(V_prev_nui)] .= V_prev_nui
    Vext[lastindex(V_prev_nui)+1:end] .= Vo
    splr.itp_Vext = linear_interpolation(tt_Vext, Vext)
end


"""
    sample_phi_top!(splr, rxw, Φi)

Sample the RX window `rxw` at the absolute TX-grid times in `Φi` using
linear interpolation.  Result is written to `splr.So_subblk` and appended
to `splr.So`.

Caller must have already verified `rxw_covers(rxw, Φi)`.
"""
function sample_phi_top!(splr, rxw, Φi)
    @inbounds for j = eachindex(Φi)
        splr.So_subblk[j] = rxw_interp(rxw, Φi[j])
    end
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
