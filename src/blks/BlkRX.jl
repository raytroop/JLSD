module BlkRX
using UnPack, DSP, Random, Interpolations
include("../util/Util_JLSD.jl"); using .Util_JLSD

export clkgen_pi_itp_top!
export sample_itp_top!, sample_phi_top!, slicers_top!, sample_filter_top!
export cdr_top!, adpt_top!
export eb_write!, eb_can_read_subblk, eb_can_read_times


# ── Elastic Buffer helpers ───────────────────────────────────────────────────

"""
    eb_write!(eb, waveform)

Append `waveform` (a block of filtered samples on the TX grid) to the
elastic buffer, advancing the write frontier `t_tx_max`.

Before writing:
- Trims `t_tx_min` to `floor(t_rx)` to reclaim samples already consumed by RX.
- Detects **underflow**: if `t_rx > t_tx_max` the RX cursor has genuinely
  outrun the TX writer; `eb.underflow_cnt` is incremented.

After writing:
- Detects **overflow**: if the stored window exceeds `eb.capacity`, the
  oldest samples are dropped and `eb.overflow_cnt` is incremented.
"""
function eb_write!(eb, waveform::AbstractVector{Float64})
    # Reclaim samples already consumed by RX
    eb.t_tx_min = max(eb.t_tx_min, floor(Int, eb.t_rx))

    # Genuine underflow: RX cursor has outrun the TX write frontier
    if eb.t_rx > eb.t_tx_max
        eb.underflow_cnt += 1
    end

    n        = length(waveform)
    capacity = eb.capacity
    t_start  = eb.t_tx_max

    for i = 0:n-1
        eb.buf[(t_start + i) % capacity + 1] = waveform[i + 1]
    end
    eb.t_tx_max += n

    # Overflow: stored window exceeds capacity, drop oldest samples
    if eb.t_tx_max - eb.t_tx_min > capacity
        eb.overflow_cnt += 1
        eb.t_tx_min = eb.t_tx_max - capacity
    end
    return nothing
end

"""
    eb_can_read_subblk(eb; margin=1) -> Bool

Return `true` when the elastic buffer contains enough data to support
the next RX sub-block with linear interpolation.

Checks that:
- `floor(t_rx) >= t_tx_min`:            oldest needed sample is available
- `floor(last sample) + 1 < t_tx_max`:  the k+1 term for interpolation is available

This function has **no side effects** (no counters modified).  Underflow and
overflow are counted exclusively in `eb_write!` at TX-block boundaries.
"""
function eb_can_read_subblk(eb; margin::Int = 1)
    osr_rx      = eb.param.osr_rx
    subblk_size = eb.param.subblk_size
    # Last sample time needed (plus one for the k+1 linear-interpolation term)
    t_need_max  = eb.t_rx + (subblk_size - 1) * osr_rx + margin
    # Oldest sample time needed
    t_need_min  = floor(Int, eb.t_rx)

    return t_need_min >= eb.t_tx_min && t_need_max < eb.t_tx_max
end

"""
    eb_can_read_times(eb, Φi) -> Bool

Return `true` when the elastic buffer contains enough data to support linear
interpolation at all times in `Φi`.

Unlike `eb_can_read_subblk`, this function uses the **actual upcoming sample
times** (including PI/CDR phase offset Φ0, static skew Φskew, and random
jitter Φrj) rather than only nominal RX cursor + spacing.  The check is:

- `floor(minimum(Φi)) >= t_tx_min`:   oldest needed sample is available
- `floor(maximum(Φi)) + 1 < t_tx_max`: the k+1 term for interpolation is present

This function has **no side effects**.
"""
function eb_can_read_times(eb, Φi)
    t_need_min = floor(Int, minimum(Φi))
    t_need_max = floor(Int, maximum(Φi)) + 1
    return t_need_min >= eb.t_tx_min && t_need_max < eb.t_tx_max
end

"""
    eb_interp(eb, t) -> Float64

Linear interpolation of the elastic buffer waveform at float TX-grid time
`t`.  Caller must ensure `t` is within `[eb.t_tx_min, eb.t_tx_max - 1]`.
"""
function eb_interp(eb, t::Float64)
    k0       = floor(Int, t)
    frac     = t - k0
    capacity = eb.capacity
    v0 = eb.buf[k0       % capacity + 1]
    v1 = eb.buf[(k0 + 1) % capacity + 1]
    return muladd(frac, v1 - v0, v0)
end

# ── RX block functions ────────────────────────────────────────────────────────


"""
    clkgen_pi_itp_top!(clkgen, eb; pi_code)

Generate RX sampling instants for one sub-block in absolute TX-grid
sample-time units, using the RX cursor `eb.t_rx` as the base position and
`param.osr_rx` as the inter-symbol spacing.

`pi_ui_cover = 4` means the PI covers 4 **RX-clock** UI units, so the CDR
phase offset `Φ0` is scaled by `osr_rx` (TX-grid samples per RX UI), not
the nominal `osr`.  Clock skew (`Φskew`) and random jitter (`Φrj`) remain
expressed in absolute TX-grid sample units (converted via `osr`, which is
the TX-grid sample rate).

After returning, `eb.t_rx` is **not** advanced here; the caller (the
block-iteration loop) is responsible for doing
    eb.t_rx += subblk_size * osr_rx
so that cursor management is centralised in one place.

The caller must also check `eb_can_read_times(eb, clkgen.Φo_subblk)` on the
generated times **before** sampling, and append to `clkgen.Φo` only after
readiness is confirmed, so that `Φo` history contains only committed
sub-blocks.
"""
function clkgen_pi_itp_top!(clkgen, eb; pi_code)
    @unpack tui, osr, osr_rx, subblk_size = clkgen.param
    @unpack nphases, rj, skews = clkgen
    @unpack pi_code_prev, pi_wrap_ui, pi_wrap_ui_Δcode = clkgen
    @unpack pi_nonlin_lut, pi_ui_cover, pi_codes_per_ui = clkgen

    Δpi_code = pi_code-pi_code_prev
    if abs(Δpi_code) > pi_wrap_ui_Δcode
        pi_wrap_ui -= sign(Δpi_code)*pi_ui_cover
    end

    # CDR phase offset in TX-grid sample units.
    # pi_ui_cover = 4 is in RX-clock UI units, so scale by osr_rx (TX-grid
    # samples per RX UI), not by nominal osr.
    Φ0    = osr_rx*(pi_wrap_ui + (pi_code + pi_nonlin_lut[pi_code+1])/pi_codes_per_ui)
    # Skew and jitter are absolute timing deviations: convert seconds → TX-grid
    # samples using osr (samples per nominal UI) and tui (nominal UI period).
    Φskew = kron(ones(Int(subblk_size/nphases)), skews/tui*osr)
    Φrj   = rj/tui*osr*randn(subblk_size)

    # Nominal sampling instants in absolute TX-grid time, stepping by osr_rx
    t_rx = eb.t_rx
    for j = 0:subblk_size-1
        clkgen.Φo_subblk[j+1] = Φ0 + t_rx + j*osr_rx + Φskew[j+1] + Φrj[j+1]
    end

    clkgen.pi_code_prev = pi_code
    clkgen.pi_wrap_ui   = pi_wrap_ui
    # NOTE: append!(clkgen.Φo, ...) is intentionally NOT done here.
    # The caller must append to Φo only after eb_can_read_times confirms the
    # buffer is ready, so that Φo history contains only committed sub-blocks.

end



"""
    sample_filter_top!(splr, Vi)

Apply the sampler bandwidth-limiting impulse response to the channel
waveform `Vi`, updating `splr.Vo` in-place (convolution with inter-block
memory preserved via `splr.Vo_mem`).
"""
function sample_filter_top!(splr, Vi)
    @unpack dt = splr.param
    @unpack ir, Vo_conv, Vo_mem = splr
    u_conv!(Vo_conv, Vi, ir, Vi_mem=Vo_mem, gain=dt)
    return nothing
end

"""
    sample_itp_top!(splr, Vi)

Apply RX bandwidth filter and build a block-local interpolation object
(`splr.itp_Vext`) that can be used by `sample_phi_top!` (backward-
compatible path, not used when the elastic buffer is active).
"""
function sample_itp_top!(splr, Vi)
    @unpack osr,dt, blk_size_osr = splr.param
    @unpack ir, Vo_conv, Vo, Vo_mem = splr
    @unpack prev_nui, V_prev_nui, Vext, tt_Vext = splr

    sample_filter_top!(splr, Vi)

    Vext[eachindex(V_prev_nui)] .= V_prev_nui
    Vext[lastindex(V_prev_nui)+1:end] .= Vo
    splr.itp_Vext = linear_interpolation(tt_Vext, Vext)
end

"""
    sample_phi_top!(splr, eb, Φi)

Sample the elastic buffer `eb` at the absolute TX-grid times in `Φi`
using linear interpolation, storing results in `splr.So_subblk` and
appending to `splr.So`.
"""
function sample_phi_top!(splr, eb, Φi)
    for j = eachindex(Φi)
        splr.So_subblk[j] = eb_interp(eb, Φi[j])
    end
    append!(splr.So, splr.So_subblk)
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