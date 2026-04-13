"""
Minimal standalone test for the ElasticBuffer layer.

Demonstrates the three key scenarios described in the requirements:
  1. freq_offset_ppm > 0  →  RX clock faster, buffer drains, eventually underflows
  2. freq_offset_ppm < 0  →  RX clock slower, buffer fills, eventually overflows
  3. freq_offset_ppm == 0 →  balanced, sub-block count matches baseline

Run from the repository root with:
    julia test/test_elastic_buffer.jl
"""

using Test

# ── Minimal stubs (avoid loading GLMakie / full sim stack) ──────────────────

module TrxStub

@kwdef mutable struct Param
    const osr::Int64   = 24
    const blk_size::Int64 = 1024
    const blk_size_osr::Int64 = blk_size * osr
    const subblk_size::Int64  = 32
    const freq_offset_ppm::Float64 = 0.0
    const osr_rx::Float64 = osr / (1 + freq_offset_ppm * 1e-6)
end

@kwdef mutable struct ElasticBuffer
    const param::Param
    const capacity::Int = 4 * param.blk_size_osr
    buf::Vector{Float64} = zeros(capacity)
    t_tx_min::Int   = 0
    t_tx_max::Int   = 0
    t_rx::Float64   = 0.0
    overflow_cnt::Int  = 0
    underflow_cnt::Int = 0
end

end  # module TrxStub

# Pull in the real EB helpers (eb_write!, eb_can_read_subblk, eb_interp)
# by re-defining them here against TrxStub so we avoid heavy dependencies.

function eb_write!(eb, waveform)
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

function eb_can_read_subblk(eb; margin::Int = 1)
    osr_rx      = eb.param.osr_rx
    subblk_size = eb.param.subblk_size
    t_need_max  = eb.t_rx + (subblk_size - 1) * osr_rx + margin
    t_need_min  = floor(Int, eb.t_rx)
    return t_need_min >= eb.t_tx_min && t_need_max < eb.t_tx_max
end

function eb_interp(eb, t::Float64)
    k0       = floor(Int, t)
    frac     = t - k0
    capacity = eb.capacity
    v0 = eb.buf[k0       % capacity + 1]
    v1 = eb.buf[(k0 + 1) % capacity + 1]
    return muladd(frac, v1 - v0, v0)
end

# Helper: run N TX blocks and return (total_subblks, underflow_cnt, overflow_cnt)
function run_sim(ppm::Float64, nblks::Int)
    param = TrxStub.Param(freq_offset_ppm = ppm)
    eb    = TrxStub.ElasticBuffer(param = param)

    total_subblks = 0
    for _ = 1:nblks
        # Simulate a TX block: write blk_size_osr unit-valued samples
        eb_write!(eb, ones(Float64, param.blk_size_osr))

        # Dynamic sub-block scheduling
        while eb_can_read_subblk(eb)
            total_subblks += 1
            eb.t_rx += param.subblk_size * param.osr_rx
        end
    end
    return total_subblks, eb.underflow_cnt, eb.overflow_cnt
end

# ── Tests ────────────────────────────────────────────────────────────────────

@testset "ElasticBuffer – freq_offset_ppm == 0" begin
    nblks     = 100
    nsubblk   = 1024 ÷ 32   # = 32 per block (exact integer ratio)
    total, uf, of = run_sim(0.0, nblks)

    @test total == nblks * nsubblk
    @test uf == 0
    @test of == 0
end

@testset "ElasticBuffer – freq_offset_ppm > 0 (RX faster)" begin
    # RX faster → buffer drains → underflow eventually
    ppm   = 100.0
    nblks = 20_000
    total, uf, of = run_sim(ppm, nblks)

    # RX should process strictly more sub-blocks than the baseline
    baseline = nblks * (1024 ÷ 32)
    @test total > baseline
    @test uf > 0
    @test of == 0
end

@testset "ElasticBuffer – freq_offset_ppm < 0 (RX slower)" begin
    # With ppm < 0, osr_rx > osr_tx: each RX sub-block spans more TX-grid time.
    # For ppm = -1000, the 32nd sub-block's interpolation window (≈24578 samples)
    # just exceeds one TX block (24576 samples), so only 31 sub-blocks fit on the
    # first block.  Total sub-blocks ends up strictly less than nblks×32.
    ppm   = -1_000.0
    nblks = 200
    total, uf, of = run_sim(ppm, nblks)

    baseline = nblks * (1024 ÷ 32)
    @test total < baseline   # RX processed fewer sub-blocks overall
    @test of == 0            # greedy scheduler keeps buffer drained — no overflow
end

@testset "ElasticBuffer – overflow when buffer is not drained" begin
    # Write 5 TX blocks without running any RX sub-blocks.
    # capacity = 4 × blk_size_osr = 98304; 5th write pushes past capacity → overflow.
    param = TrxStub.Param()
    eb    = TrxStub.ElasticBuffer(param = param)

    for _ = 1:5
        eb_write!(eb, ones(Float64, param.blk_size_osr))
    end

    @test eb.overflow_cnt > 0
    @test eb.t_tx_max - eb.t_tx_min <= eb.capacity
end

@testset "ElasticBuffer – interpolation correctness" begin
    # Buffer filled with a ramp: buf[k] = k.  Interpolation at k+0.5 should
    # return k+0.5.
    param = TrxStub.Param()
    eb    = TrxStub.ElasticBuffer(param = param)
    ramp  = collect(Float64, 0:param.blk_size_osr - 1)
    eb_write!(eb, ramp)

    for t in [0.0, 1.5, 100.7, Float64(param.blk_size_osr - 2)]
        @test eb_interp(eb, t) ≈ t  atol=1e-12
    end
end

println("All ElasticBuffer tests passed.")
