module BlkElasticBuf
using UnPack, DataStructures

export elastic_buf_write!, elastic_buf_read!

# Conversion factor from parts-per-million to a dimensionless fraction.
const PPM_TO_FRACTION = 1.0e-6

"""
    elastic_buf_write!(ebuf, Vi)

Push all OSR-rate samples in `Vi` into the elastic buffer's internal FIFO `ebuf.buf`.
Call once per block immediately after channel processing and before `elastic_buf_read!`.

If the FIFO is full (capacity reached) the oldest sample is silently overwritten by
the `CircularBuffer`, which models a hard-overflow condition.
"""
function elastic_buf_write!(ebuf, Vi)
    for v in Vi
        push!(ebuf.buf, v)
    end
    return nothing
end

"""
    elastic_buf_read!(ebuf)

Read exactly `blk_size_osr` samples from the elastic buffer into `ebuf.Vo`,
inserting (duplicating) or removing (dropping) individual OSR-rate samples as
needed to model the TX/RX clock frequency offset specified by `ebuf.ppm`.

### Algorithm
An accumulator `ebuf.accum` is incremented by `ppm × PPM_TO_FRACTION` for every
output sample requested.

- **Drop** (positive ppm, TX faster than RX):
  When `accum ≥ 1.0`, one extra sample is popped from the FIFO without being
  forwarded to `Vo`, and `accum` is decremented by 1.  This models the elastic
  buffer discarding a TX sample because the TX has produced one more sample than
  the RX has consumed.

- **Duplicate** (negative ppm, RX faster than TX):
  When `accum ≤ -1.0`, the most-recently read sample is written to `Vo` again
  without advancing the FIFO, and `accum` is incremented by 1.  This models the
  elastic buffer repeating a sample because the RX is demanding one more sample
  than the TX has produced.

- **Pass-through** (ppm = 0 or accumulator within (−1, 1)):
  One sample is popped from the FIFO and copied to `Vo` normally.

On FIFO underflow (empty buffer) the last successfully read value is repeated,
preventing simulation crashes during warm-up.
"""
function elastic_buf_read!(ebuf)
    @unpack ppm, Vo, buf = ebuf
    ppm_per_sample = ppm * PPM_TO_FRACTION
    n_out = length(Vo)
    last_val = isempty(buf) ? 0.0 : first(buf)
    idx = 1

    while idx <= n_out
        ebuf.accum += ppm_per_sample

        # TX faster than RX: drop one extra FIFO sample.
        if ebuf.accum >= 1.0 && !isempty(buf)
            popfirst!(buf)
            ebuf.accum -= 1.0
            ebuf.n_dropped += 1
        end

        # RX faster than TX: duplicate the last read sample.
        if ebuf.accum <= -1.0
            Vo[idx] = last_val
            ebuf.accum += 1.0
            ebuf.n_duplicated += 1
            idx += 1
            continue
        end

        # Normal read.
        if !isempty(buf)
            last_val = popfirst!(buf)
        end
        Vo[idx] = last_val
        idx += 1
    end

    return nothing
end

end
