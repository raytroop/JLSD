module BlkElasticBuffer
using DataStructures

export eb_write!, eb_read!, eb_status

function eb_write!(eb, data::AbstractVector)
    for sample in data
        if length(eb.fifo) >= eb.depth
            eb.overflow_cnt += 1
            popfirst!(eb.fifo)
        end
        push!(eb.fifo, sample)
    end
    eb.fill_level = length(eb.fifo)
    return nothing
end

function eb_read!(eb, n_samples::Int)
    out = Vector{Float64}(undef, n_samples)
    last_val = 0.0  # pad with zero on underflow before any sample is read
    for i in 1:n_samples
        if isempty(eb.fifo)
            eb.underflow_cnt += 1
            out[i] = last_val  # repeat last sample (or 0.0 if none yet)
        else
            val = popfirst!(eb.fifo)
            last_val = val
            out[i] = val
        end
    end
    eb.fill_level = length(eb.fifo)
    return out
end

function eb_status(eb)
    return (fill_level=eb.fill_level,
            overflow_cnt=eb.overflow_cnt,
            underflow_cnt=eb.underflow_cnt)
end

end
