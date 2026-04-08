module BlkElasticBuffer

export eb_write!, eb_read!, eb_status

function eb_write!(eb, data::AbstractVector)
    # Write oversampled samples into the ring buffer (TX side).
    # If the buffer would overflow, drop oldest samples and increment overflow_cnt.
    for s in data
        if eb.occupancy >= eb.capacity
            # Overflow: advance read pointer to drop oldest sample
            eb.rd_idx = mod(eb.rd_idx, eb.capacity) + 1
            eb.occupancy -= 1
            eb.overflow_cnt += 1
        end
        eb.wr_idx = mod(eb.wr_idx, eb.capacity) + 1
        eb.buffer[eb.wr_idx] = s
        eb.occupancy += 1
    end
end

function eb_read!(eb, n_samples::Int)
    # Read n_samples from the ring buffer (RX side).
    # If not enough samples, pad by repeating the last valid sample and increment underflow_cnt.
    out = Vector{Float64}(undef, n_samples)
    last_val = 0.0
    for i in 1:n_samples
        if eb.occupancy > 0
            eb.rd_idx = mod(eb.rd_idx, eb.capacity) + 1
            last_val = eb.buffer[eb.rd_idx]
            eb.occupancy -= 1
            out[i] = last_val
        else
            eb.underflow_cnt += 1
            out[i] = last_val  # repeat last sample
        end
    end
    return out
end

function eb_status(eb)
    return (fill_level=eb.occupancy, overflow_cnt=eb.overflow_cnt, underflow_cnt=eb.underflow_cnt)
end

end
