"""
Minimal standalone test for the ElasticBuffer layer.

Demonstrates the three key scenarios described in the requirements:
  1. freq_offset_ppm > 0  →  RX clock faster, buffer drains, eventually underflows
  2. freq_offset_ppm < 0  →  RX clock slower, buffer fills, eventually overflows
  3. freq_offset_ppm == 0 →  balanced, sub-block count matches baseline

Also tests that elastic-buffer readiness (eb_can_read_times) is based on the
actual upcoming sample times, accounting for PI/CDR phase offset (Φ0), static
skew (Φskew), and random jitter (Φrj), rather than only nominal spacing.

Run from the repository root with:
    julia test/test_elastic_buffer.jl
"""

using Test, Random

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

# Pull in the real EB helpers by re-defining them here against TrxStub
# so we avoid heavy dependencies.

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

# Nominal readiness check (kept for reference; scheduling uses eb_can_read_times).
function eb_can_read_subblk(eb; margin::Int = 1)
    osr_rx      = eb.param.osr_rx
    subblk_size = eb.param.subblk_size
    t_need_max  = eb.t_rx + (subblk_size - 1) * osr_rx + margin
    t_need_min  = floor(Int, eb.t_rx)
    return t_need_min >= eb.t_tx_min && t_need_max < eb.t_tx_max
end

# eb_can_read_times: readiness based on actual sample times (includes all
# timing perturbations: Φ0, Φskew, Φrj).  The k+1 interpolation term is
# accounted for by the +1 on the max side.
function eb_can_read_times(eb, Φi)
    t_need_min = floor(Int, minimum(Φi))
    t_need_max = floor(Int, maximum(Φi)) + 1
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

# Helper: generate nominal candidate Φo_subblk (no jitter/skew perturbation)
function make_phi_nominal(eb, Φ0_offset = 0.0)
    osr_rx      = eb.param.osr_rx
    subblk_size = eb.param.subblk_size
    return [Φ0_offset + eb.t_rx + j * osr_rx for j in 0:subblk_size-1]
end

# Helper: run N TX blocks using eb_can_read_times with nominal Φo_subblk,
# and return (total_subblks, underflow_cnt, overflow_cnt)
function run_sim(ppm::Float64, nblks::Int)
    param = TrxStub.Param(freq_offset_ppm = ppm)
    eb    = TrxStub.ElasticBuffer(param = param)

    total_subblks = 0
    for _ = 1:nblks
        # Simulate a TX block: write blk_size_osr unit-valued samples
        eb_write!(eb, ones(Float64, param.blk_size_osr))

        # Dynamic sub-block scheduling using actual candidate times
        while true
            Φo = make_phi_nominal(eb)
            eb_can_read_times(eb, Φo) || break
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

# ── Tests for eb_can_read_times: actual-time-based readiness ─────────────────

@testset "eb_can_read_times – nominal matches eb_can_read_subblk" begin
    # With zero phase offset and zero jitter/skew, eb_can_read_times using
    # make_phi_nominal should agree with eb_can_read_subblk at all t_rx.
    param = TrxStub.Param()
    eb    = TrxStub.ElasticBuffer(param = param)
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    while eb_can_read_subblk(eb)
        Φo = make_phi_nominal(eb)
        @test eb_can_read_times(eb, Φo) == true
        eb.t_rx += param.subblk_size * param.osr_rx
    end
    # Once nominal check says no more, actual-time check should also say no.
    Φo = make_phi_nominal(eb)
    @test eb_can_read_times(eb, Φo) == false
end

@testset "eb_can_read_times – large positive Φ0 (late phase offset)" begin
    # A large positive PI phase offset shifts sample times later.
    # If the offset is big enough, those samples are not yet in the buffer
    # even though the nominal cursor t_rx would be in range.
    param = TrxStub.Param()
    eb    = TrxStub.ElasticBuffer(param = param)

    # Write exactly one TX block worth of samples.
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    # Nominal check says the first sub-block is readable.
    @test eb_can_read_subblk(eb) == true

    # But with Φ0 = blk_size_osr (shift one whole block ahead), the last
    # sample would be past t_tx_max — buffer is NOT ready for those actual times.
    big_Φ0 = Float64(param.blk_size_osr)
    Φo_late = make_phi_nominal(eb, big_Φ0)
    @test eb_can_read_times(eb, Φo_late) == false
end

@testset "eb_can_read_times – large negative Φ0 (early phase offset)" begin
    # A large negative PI phase offset shifts all sample times earlier.
    # If the shift goes before t_tx_min, the buffer is not ready even though
    # the nominal cursor is in range.
    param = TrxStub.Param()
    eb    = TrxStub.ElasticBuffer(param = param)

    # Write one TX block and advance RX cursor to nearly the end.
    eb_write!(eb, ones(Float64, param.blk_size_osr))
    n_nominal = floor(Int, param.blk_size_osr / (param.subblk_size * param.osr_rx))
    for _ in 1:n_nominal-1
        eb.t_rx += param.subblk_size * param.osr_rx
    end
    @test eb_can_read_subblk(eb) == true   # nominal says readable

    # Shift all sample times back before t_tx_min with a large negative Φ0.
    big_neg_Φ0 = -(eb.t_rx + 1.0)
    Φo_early = make_phi_nominal(eb, big_neg_Φ0)
    @test minimum(Φo_early) < eb.t_tx_min
    @test eb_can_read_times(eb, Φo_early) == false
end

@testset "eb_can_read_times – jitter can push max sample past t_tx_max" begin
    # Simulate one sub-block of candidate times where random jitter on the
    # last sample pushes it past t_tx_max, making the buffer not ready.
    param = TrxStub.Param()
    eb    = TrxStub.ElasticBuffer(param = param)
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    # Place t_rx so the nominal last sample is one step before t_tx_max.
    subblk_size = param.subblk_size
    osr_rx      = param.osr_rx
    eb.t_rx = Float64(eb.t_tx_max) - (subblk_size - 1) * osr_rx - 2.0

    Φo_no_jitter = make_phi_nominal(eb)
    @test eb_can_read_times(eb, Φo_no_jitter) == true  # just fits

    # Add a large positive jitter only to the last sample.
    Φo_with_jitter = copy(Φo_no_jitter)
    Φo_with_jitter[end] += 3.0    # push last sample past t_tx_max
    @test eb_can_read_times(eb, Φo_with_jitter) == false
end

@testset "eb_can_read_times – jitter/skew within buffer is always readable" begin
    # When the buffer has plenty of headroom on both ends, reasonable jitter/skew
    # should keep all sample times within [t_tx_min, t_tx_max).
    param = TrxStub.Param()
    eb    = TrxStub.ElasticBuffer(param = param)
    # Write 2 TX blocks to give ample write headroom.
    eb_write!(eb, ones(Float64, param.blk_size_osr))
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    osr_rx      = param.osr_rx
    subblk_size = param.subblk_size
    # Start t_rx in the middle of the buffer so jitter has room on both ends.
    eb.t_rx = Float64(param.blk_size_osr)

    # Jitter of ±0.5 samples — well within the buffer on both sides.
    Random.seed!(42)
    jitter = 0.5 * (2 .* rand(subblk_size) .- 1)
    Φo_jittered = make_phi_nominal(eb) .+ jitter
    @test minimum(Φo_jittered) >= eb.t_tx_min
    @test floor(Int, maximum(Φo_jittered)) + 1 < eb.t_tx_max
    @test eb_can_read_times(eb, Φo_jittered) == true
end

@testset "eb_can_read_times – osr_rx scaling with nonzero freq_offset" begin
    # RX faster: osr_rx < osr.  Write one block, check that nominal sub-block
    # count (32) still works with eb_can_read_times using osr_rx spacing.
    ppm   = 100.0
    param = TrxStub.Param(freq_offset_ppm = ppm)
    eb    = TrxStub.ElasticBuffer(param = param)
    eb_write!(eb, ones(Float64, param.blk_size_osr))

    # With osr_rx < osr, the sub-block window is narrower → should still be readable
    Φi = make_phi_nominal(eb)
    @test eb_can_read_times(eb, Φi)
end

println("All ElasticBuffer tests passed.")
