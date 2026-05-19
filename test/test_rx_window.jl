"""
Standalone tests for the RxWindow layer.

Demonstrates:
  1. freq_offset_ppm > 0  →  RX clock faster, RX processes slightly more
                             sub-blocks than nsubblk per TX block.
  2. freq_offset_ppm < 0  →  RX clock slower, RX processes slightly fewer.
  3. freq_offset_ppm == 0 →  exact integer count.
  4. Greedy scheduler keeps buffer occupancy bounded; overrun is fatal.
  5. rxw_covers uses the actual upcoming sample times (including Φ0,
     Φskew, Φrj), not nominal spacing.

Run from the repo root:
    julia test/test_rx_window.jl
"""

using Test, Random

# ── Minimal stubs (avoid loading GLMakie / full sim stack) ───────────────────

module TrxStub

@kwdef mutable struct Param
    const osr::Int64           = 24
    const blk_size::Int64      = 1024
    const blk_size_osr::Int64  = blk_size * osr
    const subblk_size::Int64   = 32
    freq_offset_ppm::Float64   = 0.0
    osr_rx::Float64            = osr / (1 + freq_offset_ppm * 1e-6)
end

@kwdef mutable struct RxWindow
    const param::Param
    const capacity::Int = 4 * param.blk_size_osr
    buf::Vector{Float64} = zeros(capacity)
    t_min::Int     = 0
    t_max::Int     = 0
    t_rx::Float64  = 0.0
end

end  # module TrxStub

# Re-define rxw helpers against TrxStub so we don't import the full module
# graph (which would pull in GLMakie etc.).

function rxw_extend!(rxw, waveform)
    rxw.t_min = max(rxw.t_min, floor(Int, rxw.t_rx))

    n = length(waveform)
    if (rxw.t_max + n) - rxw.t_min > rxw.capacity
        error("RxWindow overrun: t_min=$(rxw.t_min) t_max=$(rxw.t_max) " *
              "t_rx=$(rxw.t_rx) n=$n capacity=$(rxw.capacity)")
    end

    cap = rxw.capacity
    t_start = rxw.t_max
    for i = 0:n-1
        rxw.buf[(t_start + i) % cap + 1] = waveform[i + 1]
    end
    rxw.t_max += n
    return nothing
end

function rxw_covers(rxw, Φi)
    t_need_min = floor(Int, minimum(Φi))
    t_need_max = floor(Int, maximum(Φi)) + 1
    return t_need_min >= rxw.t_min && t_need_max < rxw.t_max
end

function rxw_interp(rxw, t::Float64)
    k0   = floor(Int, t)
    frac = t - k0
    cap  = rxw.capacity
    v0   = rxw.buf[k0       % cap + 1]
    v1   = rxw.buf[(k0 + 1) % cap + 1]
    return muladd(frac, v1 - v0, v0)
end

# Nominal sample times: t_rx, t_rx + osr_rx, ..., t_rx + (subblk_size-1)*osr_rx
function make_phi_nominal(rxw, Φ0_offset = 0.0)
    osr_rx      = rxw.param.osr_rx
    subblk_size = rxw.param.subblk_size
    return [Φ0_offset + rxw.t_rx + j * osr_rx for j in 0:subblk_size-1]
end

# Run nblks TX-block iterations with the greedy scheduler.
# Returns (total_subblks_completed, did_overrun_fire)
function run_sim(ppm::Float64, nblks::Int)
    param = TrxStub.Param(freq_offset_ppm = ppm)
    rxw   = TrxStub.RxWindow(param = param)

    total_subblks = 0
    for _ = 1:nblks
        rxw_extend!(rxw, ones(Float64, param.blk_size_osr))
        while true
            Φo = make_phi_nominal(rxw)
            rxw_covers(rxw, Φo) || break
            total_subblks += 1
            rxw.t_rx += param.subblk_size * param.osr_rx
        end
    end
    return total_subblks
end

# ── Tests ────────────────────────────────────────────────────────────────────

@testset "RxWindow — freq_offset_ppm = 0" begin
    nblks   = 100
    nsubblk = 1024 ÷ 32
    @test run_sim(0.0, nblks) == nblks * nsubblk
end

@testset "RxWindow — freq_offset_ppm > 0 (RX faster)" begin
    # Greedy scheduler should fit slightly more sub-blocks than baseline.
    ppm   = 100.0
    nblks = 20_000
    baseline = nblks * (1024 ÷ 32)
    total = run_sim(ppm, nblks)
    @test total > baseline
end

@testset "RxWindow — freq_offset_ppm < 0 (RX slower)" begin
    # Greedy scheduler should fit slightly fewer sub-blocks than baseline.
    ppm   = -1_000.0
    nblks = 200
    baseline = nblks * (1024 ÷ 32)
    @test run_sim(ppm, nblks) < baseline
end

@testset "RxWindow — overrun is fatal (no silent drop)" begin
    # Make capacity tight (one TX block) so a second un-drained write must
    # overrun.  The default 4*blk_size_osr is by design impossible to overrun
    # under any reasonable use, so we shrink it here to exercise the guard.
    param = TrxStub.Param()
    rxw   = TrxStub.RxWindow(param = param,
                             capacity = param.blk_size_osr,
                             buf = zeros(param.blk_size_osr))
    rxw_extend!(rxw, ones(Float64, param.blk_size_osr))
    @test_throws ErrorException rxw_extend!(rxw, ones(Float64, param.blk_size_osr))
end

@testset "RxWindow — capacity 4*blk_size_osr never overruns at max |ppm|" begin
    # Slider range is ±500 ppm; verify the default capacity survives both ends.
    for ppm in (500.0, -500.0)
        @test run_sim(ppm, 200) > 0  # completed without ErrorException
    end
end

@testset "RxWindow — linear interpolation correctness" begin
    param = TrxStub.Param()
    rxw   = TrxStub.RxWindow(param = param)
    ramp  = collect(Float64, 0:param.blk_size_osr - 1)
    rxw_extend!(rxw, ramp)

    for t in [0.0, 1.5, 100.7, Float64(param.blk_size_osr - 2)]
        @test rxw_interp(rxw, t) ≈ t  atol = 1e-12
    end
end

@testset "rxw_covers — large positive Φ0 fails readiness" begin
    param = TrxStub.Param()
    rxw   = TrxStub.RxWindow(param = param)
    rxw_extend!(rxw, ones(Float64, param.blk_size_osr))

    # Nominal sub-block is ready.
    @test rxw_covers(rxw, make_phi_nominal(rxw))

    # Shift sub-block one whole block into the future → last sample is past t_max.
    Φo_late = make_phi_nominal(rxw, Float64(param.blk_size_osr))
    @test rxw_covers(rxw, Φo_late) == false
end

@testset "rxw_covers — large negative Φ0 fails readiness" begin
    param = TrxStub.Param()
    rxw   = TrxStub.RxWindow(param = param)
    rxw_extend!(rxw, ones(Float64, param.blk_size_osr))

    # Walk t_rx forward then push it back via negative Φ0.
    n_steps = floor(Int, param.blk_size_osr / (param.subblk_size * param.osr_rx)) - 1
    for _ in 1:n_steps
        rxw.t_rx += param.subblk_size * param.osr_rx
    end

    Φo_early = make_phi_nominal(rxw, -(rxw.t_rx + 1.0))
    @test minimum(Φo_early) < rxw.t_min
    @test rxw_covers(rxw, Φo_early) == false
end

@testset "rxw_covers — jitter that pushes max past t_max fails" begin
    param = TrxStub.Param()
    rxw   = TrxStub.RxWindow(param = param)
    rxw_extend!(rxw, ones(Float64, param.blk_size_osr))

    subblk_size = param.subblk_size
    osr_rx      = param.osr_rx
    rxw.t_rx = Float64(rxw.t_max) - (subblk_size - 1) * osr_rx - 2.0

    Φo_clean = make_phi_nominal(rxw)
    @test rxw_covers(rxw, Φo_clean) == true

    Φo_jittered = copy(Φo_clean)
    Φo_jittered[end] += 3.0
    @test rxw_covers(rxw, Φo_jittered) == false
end

# ── Production BlkRX state-machine regressions ───────────────────────────────

include(joinpath(@__DIR__, "..", "src", "blks", "BlkRX.jl"))

module ProdRxStub

@kwdef mutable struct Param
    const osr::Int64           = 24
    const blk_size::Int64      = 1024
    const blk_size_osr::Int64  = blk_size * osr
    const subblk_size::Int64   = 32
    const tui::Float64         = 1 / 56e9
    freq_offset_ppm::Float64   = 0.0
    osr_rx::Float64            = osr / (1 + freq_offset_ppm * 1e-6)
end

@kwdef mutable struct RxWindow
    const param::Param
    const capacity::Int = 4 * param.blk_size_osr
    buf::Vector{Float64} = zeros(capacity)
    t_min::Int     = 0
    t_max::Int     = 0
    t_rx::Float64  = 0.0
end

@kwdef mutable struct Clkgen
    const param::Param
    nphases::Int             = 4
    rj::Float64              = 0.0
    skews::Vector{Float64}   = zeros(nphases)
    pi_res::Int              = 8
    pi_max_code::Int         = 2^pi_res - 1
    pi_ui_cover::Int         = 4
    pi_codes_per_ui::Float64 = 2^pi_res / pi_ui_cover
    pi_nonlin_lut::Vector{Float64} = zeros(2^pi_res)
    pi_code_prev::Int        = 0
    pi_wrap_ui_Δcode::Int    = pi_max_code - 10
    Φo_subblk::Vector{Float64} = zeros(param.subblk_size)
    Φo_subblk_valid::Bool    = false
    t_rx_subblk::Float64     = 0.0
    t_rx_next::Float64       = 0.0
    pi_code_subblk::Int      = 0
end

end

@testset "production clkgen candidate is not committed on underrun" begin
    param  = ProdRxStub.Param()
    clkgen = ProdRxStub.Clkgen(param = param)
    rxw    = ProdRxStub.RxWindow(param = param)

    clkgen.pi_code_prev = 250
    BlkRX.clkgen_pi_itp_top!(clkgen, rxw; pi_code = 0)

    @test BlkRX.rxw_covers(rxw, clkgen.Φo_subblk) == false
    @test rxw.t_rx == 0.0
    @test clkgen.pi_code_prev == 250
    @test clkgen.Φo_subblk_valid

    Φo_pending = copy(clkgen.Φo_subblk)
    BlkRX.clkgen_pi_itp_top!(clkgen, rxw; pi_code = 0)
    @test clkgen.Φo_subblk == Φo_pending

    BlkRX.rxw_extend!(rxw, ones(Float64, param.blk_size_osr))
    @test BlkRX.rxw_covers(rxw, clkgen.Φo_subblk)

    t_rx_next = clkgen.t_rx_next
    BlkRX.clkgen_commit_subblk!(clkgen, rxw)
    @test rxw.t_rx == t_rx_next
    @test clkgen.pi_code_prev == 0
    @test clkgen.Φo_subblk_valid == false
end

@testset "production rxw_extend keeps low-side phase margin" begin
    param = ProdRxStub.Param()
    rxw   = ProdRxStub.RxWindow(param = param)

    BlkRX.rxw_extend!(rxw, ones(Float64, param.blk_size_osr))
    rxw.t_rx = 10.0
    BlkRX.rxw_extend!(rxw, ones(Float64, param.blk_size_osr);
                      phase_margin = 1.0)

    Φi = [9.75, 10.0]
    @test rxw.t_min == 9
    @test BlkRX.rxw_covers(rxw, Φi)
end

println("All RxWindow tests passed.")
