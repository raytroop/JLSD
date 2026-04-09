module BlkElasticBuffer
using DataStructures

export eb_write!, eb_read!, eb_occupancy

function eb_write!(ebuf, data::Vector)
    for s in data
        if isfull(ebuf.buf)
            ebuf.n_overflow += 1
            popfirst!(ebuf.buf)  # drop oldest sample
        end
        push!(ebuf.buf, s)
    end
end

function eb_read!(ebuf, n_samples::Int)
    out = Vector{Float64}(undef, n_samples)
    available = length(ebuf.buf)
    n_read = min(n_samples, available)
    for i in 1:n_read
        out[i] = popfirst!(ebuf.buf)
    end
    if n_read < n_samples
        ebuf.n_underflow += 1
        last_val = n_read > 0 ? out[n_read] : 0.0
        for i in n_read+1:n_samples
            out[i] = last_val
        end
    end
    return out
end

function eb_occupancy(ebuf)
    return length(ebuf.buf)
end

end
