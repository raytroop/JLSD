"""
Standalone test for the ElasticBuffer layer.

Re-implements the lightweight Param/ElasticBuffer types so the test does
not pull in GLMakie/Makie/DSP, and exercises the same algorithms as
`src/blks/BlkElasticBuffer.jl`.

Covers the boundary conditions called out in the design:
  1. freq_offset_ppm ==  0  → integer ratio, exact baseline sub-block count
  2. freq_offset_ppm  >  0  → RX faster, buffer drains, eventually underflow
  3. freq_offset_ppm  <  0  → RX slower, buffer fills (can overflow if
     never drained)
  4. eb_can_read_times: actual candidate Φi (with Φ0/skew/jitter) is the
     thing that gates readiness, not nominal cursor + spacing
  5. Linear interpolation correctness on a ramp
  6. Ring wrap-around correctness across many capacity-spans of writes

Run with:
    julia test/test_elastic_buffer.jl
"""

using Test, Random


# ── Standalone copies of Param + ElasticBuffer (no external deps) ───────────

Base.@kwdef mutable struct TParam
    osr::Int             = 24
    blk_size::Int        = 1024
    blk_size_osr::Int    = blk_size * osr
    subblk_size::Int     = 32
    freq_offset_ppm::Float64 = 0.0
    osr_rx::Float64      = osr / (1 + freq_offset_ppm * 1e-6)
end

Base.@kwdef mutable struct TEB
    param::TParam
    capacity::Int        = 4 * param.blk_size_osr
    buf::Vector{Float64} = zeros(capacity)
    t_tx_min::Int        = 0
    t_tx_max::Int        = 0
    t_rx::Float64        = 0.0
    overflow_cnt::Int    = 0
    underflow_cnt::Int   = 0
end


# ── Algorithms (mirror src/blks/BlkElasticBuffer.jl) ────────────────────────

function eb_write!(eb, waveform::AbstractVector{Float64})
    eb.t_tx_min = max(eb.t_tx_min, floor(Int, eb.t_rx))

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

    if eb.t_tx_max - eb.t_tx_min > capacity
        eb.overflow_cnt += 1
        eb.t_tx_min = eb.t_tx_max - capacity
    end
    return nothing
end

function eb_can_read_times(eb, Φi)
    t_need_min = floor(Int, minimum(Φi))
    t_need_max = floor(Int, maximum(Φi)) + 1
    return t_need_min >= eb.t_tx_min && t_need_max < eb.t_tx_max
end

function eb_interp(eb, t::Float64)
    k0       = floor(Int, t)
    α        = t - k0
    capacity = eb.capacity
    v0 = eb.buf[ k0       % capacity + 1]
    v1 = eb.buf[(k0 + 1)  % capacity + 1]
    return muladd(α, v1 - v0, v0)
end


# ── Helpers ─────────────────────────────────────────────────────────────────

function make_phi_nominal(eb, Φ0_offset::Float64 = 0.0)
    osr_rx      = eb.param.osr_rx
    subblk_size = eb.param.subblk_size
    return [Φ0_offset + eb.t_rx + j * osr_rx for j in 0:subblk_size-1]
end

function run_sim(ppm::Float64, nblks::Int)
    param = TParam(freq_offset_ppm = ppm)
    eb    = TEB(param = param)

    total_subblks = 0
    for _ = 1:nblks
        eb_write!(eb, ones(Float64, param.blk_size_osr))
        while true
            Φo = make_phi_nominal(eb)
            eb_can_read_times(eb, Φo) || break
            total_subblks += 1
            eb.t_rx += param.subblk_size * param.osr_rx
        end
    end
    return total_subblks, eb.underflow_cnt, eb.overflow_cnt
end


# ── Tests ───────────────────────────────────────────────────────────────────

@testset "ElasticBuffer — zero ppm matches baseline exactly" begin
    nblks   = 100
    nsubblk = 1024 ÷ 32
    total, uf, of = run_sim(0.0, nblks)

    @test total == nblks * nsubblk
    @test uf == 0
    @test of == 0
end

@testset "Variable-length RX symbol buffers resize per delivered block" begin
    bits_per_sym = 1
    si = UInt8[]
    si_bits = Bool[]
    ref_bits = Bool[]

    function resize_rx_buffers!(si, si_bits, ref_bits, received)
        append!(si, received)
        nsym = length(si)
        if nsym == 0
            empty!(si)
            return length(si_bits), length(ref_bits), length(si)
        end
        nbits = bits_per_sym * nsym
        resize!(si_bits, nbits)
        resize!(ref_bits, nbits)
        empty!(si)
        return length(si_bits), length(ref_bits), length(si)
    end

    @test resize_rx_buffers!(si, si_bits, ref_bits, UInt8[]) == (0, 0, 0)
    @test resize_rx_buffers!(si, si_bits, ref_bits, UInt8[0, 1, 1]) == (3, 3, 0)
    @test resize_rx_buffers!(si, si_bits, ref_bits, UInt8[1]) == (1, 1, 0)
    @test resize_rx_buffers!(si, si_bits, ref_bits, UInt8.(ones(1024))) == (1024, 1024, 0)
end

@testset "ElasticBuffer — positive ppm (RX faster) eventually underflows" begin
    freq_offset_ppm = 100.0
    nblks = 20_000
    total, uf, of = run_sim(freq_offset_ppm, nblks)

    baseline = nblks * (1024 ÷ 32)
    @test total > baseline
    @test uf > 0
    @test of == 0
end

@testset "ElasticBuffer — negative ppm (RX slower) drops sub-blocks" begin
    # With ppm = -1000, osr_rx ≈ 24.024 > osr.  Each sub-block now spans
    # subblk_size * osr_rx = 32*24.024 ≈ 768.77 TX-grid samples plus the
    # +1 interpolation margin.  Exactly 32 sub-blocks would need
    # 32*768.77 ≈ 24600.7 samples — slightly more than one TX block
    # (24576 samples).  Greedy scheduler therefore commits 31 sub-blocks
    # in some blocks and 32 in others, so total < nblks*32.
    freq_offset_ppm = -1_000.0
    nblks = 200
    total, uf, of = run_sim(freq_offset_ppm, nblks)

    baseline = nblks * (1024 ÷ 32)
    @test total < baseline
    @test of == 0
end

@testset "ElasticBuffer — overflow when buffer is never drained" begin
    param = TParam()
    eb    = TEB(param = param)

    for _ = 1:5
        eb_write!(eb, ones(Float64, param.blk_size_osr))
    end

    @test eb.overflow_cnt > 0
    @test eb.t_tx_max - eb.t_tx_min <= eb.capacity
end

@testset "ElasticBuffer — interpolation accuracy on a ramp" begin
    param = TParam()
    eb    = TEB(param = param)
    test_ramp = collect(Float64, 0:param.blk_size_osr - 1)
    eb_write!(eb, test_ramp)

    for t in (0.0, 1.5, 100.7, Float64(param.blk_size_osr - 2))
        @test eb_interp(eb, t) ≈ t  atol=1e-12
    end
end

@testset "ElasticBuffer — ring wrap stays consistent over many capacities" begin
    # Drive in unique data each cycle and advance t_rx between cycles so
    # the ring really wraps; check that the most recent block still
    # round-trips through eb_interp at integer slots.
    param = TParam()
    eb    = TEB(param = param)
    blk_n = param.blk_size_osr

    for cycle in 1:6
        wave = Float64.((cycle * 1000) .+ collect(0:blk_n-1))
        eb_write!(eb, wave)
        eb.t_rx += blk_n
    end

    cycle = 6
    base  = eb.t_tx_max - blk_n
    for k in (0, 1, 7, 1234, blk_n - 2)
        t   = Float64(base + k)
        exp = Float64(cycle * 1000 + k)
        @test eb_interp(eb, t) ≈ exp atol=1e-12
    end
end


# ── eb_can_read_times: actual-time-based readiness ─────────────────────────

@testset "eb_can_read_times — nominal aligns with greedy schedule" begin
    param = TParam()
    eb    = TEB(param = param)
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    nominal_steps = 0
    while true
        Φo = make_phi_nominal(eb)
        eb_can_read_times(eb, Φo) || break
        nominal_steps += 1
        eb.t_rx += param.subblk_size * param.osr_rx
    end
    @test nominal_steps == 1024 ÷ 32
end

@testset "eb_can_read_times — sampler pre-history preserves Φ0 headroom" begin
    param = TParam()
    initial_delay = 8 * param.osr
    eb = TEB(param = param, t_tx_max = initial_delay)
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    # The default CDR code in the full model starts near the center of a
    # 4-UI PI range, i.e. about +2 UI.  Initializing t_tx_max with the
    # sampler's legacy pre-history keeps all 32 zero-ppm sub-blocks
    # schedulable in the first TX block despite that phase offset.
    Φ0 = 2 * param.osr_rx
    steps = 0
    while true
        Φo = make_phi_nominal(eb, Φ0)
        eb_can_read_times(eb, Φo) || break
        steps += 1
        eb.t_rx += param.subblk_size * param.osr_rx
    end
    @test steps == param.blk_size ÷ param.subblk_size
end

@testset "eb_can_read_times — large positive Φ0 must block the read" begin
    param = TParam()
    eb    = TEB(param = param)
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    @test eb_can_read_times(eb, make_phi_nominal(eb)) == true
    big_Φ0 = Float64(param.blk_size_osr)
    @test eb_can_read_times(eb, make_phi_nominal(eb, big_Φ0)) == false
end

@testset "eb_can_read_times — large negative Φ0 must block the read" begin
    param = TParam()
    eb    = TEB(param = param)
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    n_nominal = floor(Int, param.blk_size_osr / (param.subblk_size * param.osr_rx))
    eb.t_rx   += (n_nominal - 1) * param.subblk_size * param.osr_rx
    @test eb_can_read_times(eb, make_phi_nominal(eb)) == true

    big_neg_Φ0 = -(eb.t_rx + 1.0)
    Φo_early   = make_phi_nominal(eb, big_neg_Φ0)
    @test minimum(Φo_early) < eb.t_tx_min
    @test eb_can_read_times(eb, Φo_early) == false
end

@testset "eb_can_read_times — jitter past t_tx_max blocks the read" begin
    param = TParam()
    eb    = TEB(param = param)
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    eb.t_rx = Float64(eb.t_tx_max) - (param.subblk_size - 1) * param.osr_rx - 2.0

    Φo = make_phi_nominal(eb)
    @test eb_can_read_times(eb, Φo) == true

    Φo2 = copy(Φo)
    Φo2[end] += 3.0
    @test eb_can_read_times(eb, Φo2) == false
end

@testset "eb_can_read_times — moderate jitter inside buffer is always OK" begin
    param = TParam()
    eb    = TEB(param = param)
    eb_write!(eb, ones(Float64, param.blk_size_osr))
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    eb.t_rx = Float64(param.blk_size_osr)

    Random.seed!(42)
    jitter = 0.5 .* (2 .* rand(param.subblk_size) .- 1)
    Φo     = make_phi_nominal(eb) .+ jitter

    @test minimum(Φo) >= eb.t_tx_min
    @test floor(Int, maximum(Φo)) + 1 < eb.t_tx_max
    @test eb_can_read_times(eb, Φo) == true
end


println("All ElasticBuffer tests passed.")
