module BlkElasticBuffer

export eb_write!, eb_can_read_times, eb_interp


"""
    eb_write!(eb, waveform)

Append `waveform` (one TX block of RX-filtered samples on the TX grid)
to the elastic buffer, advancing the write frontier `t_tx_max`.

Before writing:
  * `t_tx_min` is bumped forward to `floor(t_rx)` so RX-consumed samples
    can be reclaimed (this is the natural time to do it because TX is
    about to overwrite their ring slots).
  * If `t_rx > t_tx_max` the RX cursor has genuinely outrun the writer
    — a real underflow — `underflow_cnt` is incremented.

After writing:
  * If the live window `t_tx_max - t_tx_min` exceeds `capacity` the ring
    has wrapped past unconsumed data; `t_tx_min` is snapped forward to
    `t_tx_max - capacity` and `overflow_cnt` is incremented.

These per-block side effects are the *only* place the counters are
touched — readiness checks remain side-effect free.
"""
function eb_write!(eb, waveform::AbstractVector{Float64})
    # Reclaim slots already consumed by RX
    eb.t_tx_min = max(eb.t_tx_min, floor(Int, eb.t_rx))

    # Genuine underflow: RX cursor outran TX writer
    if eb.t_rx > eb.t_tx_max
        eb.underflow_cnt += 1
    end

    n        = length(waveform)
    capacity = eb.capacity
    t_start  = eb.t_tx_max
    @inbounds for i = 0:n-1
        eb.buf[(t_start + i) % capacity + 1] = waveform[i + 1]
    end
    eb.t_tx_max += n

    # Overflow: live window exceeds capacity, drop oldest samples
    if eb.t_tx_max - eb.t_tx_min > capacity
        eb.overflow_cnt += 1
        eb.t_tx_min = eb.t_tx_max - capacity
    end
    return nothing
end


"""
    eb_can_read_times(eb, Φi) -> Bool

Return `true` iff every floating-point candidate sample time in `Φi`
falls strictly inside the live ring window so that linear interpolation
(which needs both `buf[k]` and `buf[k+1]`) is safe.

This check uses the *actual* upcoming sample times — including PI/CDR
phase offset Φ0, static skews Φskew and random jitter Φrj — rather than
nominal `t_rx + n·osr_rx` only.  This guarantees that no jitter or skew
can push a sample outside the window between the readiness check and
the actual read.

  * `floor(min(Φi))     >= t_tx_min`  — oldest needed sample present
  * `floor(max(Φi)) + 1 <  t_tx_max`  — k+1 interpolation term present

`< t_tx_max` (strict) plus the explicit `+1` margin protect against
both off-by-one at the write frontier and the linear-interp k+1 access.
"""
function eb_can_read_times(eb, Φi)
    t_need_min = floor(Int, minimum(Φi))
    t_need_max = floor(Int, maximum(Φi)) + 1
    return t_need_min >= eb.t_tx_min && t_need_max < eb.t_tx_max
end


"""
    eb_interp(eb, t) -> Float64

Linear interpolation of the buffered waveform at the float TX-grid time
`t`.  Caller must ensure `t` lies within `[t_tx_min, t_tx_max - 1]`
(use `eb_can_read_times` to guarantee this).
"""
@inline function eb_interp(eb, t::Float64)
    k0       = floor(Int, t)
    α        = t - k0
    capacity = eb.capacity
    # eb_can_read_times guarantees k0 and k0+1 are present in the live ring window.
    @inbounds v0 = eb.buf[ k0       % capacity + 1]
    @inbounds v1 = eb.buf[(k0 + 1)  % capacity + 1]
    return muladd(α, v1 - v0, v0)
end


end
