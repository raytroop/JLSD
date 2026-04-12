using Test

stub_root = mktempdir()

mkpath(joinpath(stub_root, "UnPack", "src"))
write(joinpath(stub_root, "UnPack", "src", "UnPack.jl"), raw"""
module UnPack
export @unpack
macro unpack(ex)
    lhs, rhs = ex.args
    vars = lhs isa Expr && lhs.head == :tuple ? lhs.args : [lhs]
    Expr(:block, [:( $(esc(var)) = $(Expr(:., esc(rhs), QuoteNode(var))) ) for var in vars]...)
end
end
""")

mkpath(joinpath(stub_root, "Distributions", "src"))
write(joinpath(stub_root, "Distributions", "src", "Distributions.jl"), """
module Distributions
end
""")

pushfirst!(LOAD_PATH, stub_root)

include(joinpath(@__DIR__, "..", "src", "blks", "BlkBIST.jl"))
using .BlkBIST

mutable struct StubBist
    polynomial::Vector{UInt8}
    inv::Bool
    chk_seed::Vector{Bool}
    ref_bits::Vector{Bool}
    Si_bits::Vector{Bool}
    chk_lock_status::Bool
    chk_lock_cnt::Int
    chk_lock_cnt_threshold::Int
    ber_err_cnt::Int
    ber_bit_cnt::Int
end

@testset "ber_check_prbs! counts trailing bits after lock" begin
    polynomial = UInt8[6, 7]
    seed = ones(Bool, 7)
    nbits = 10
    threshold = 3
    initial_lock_cnt = 1

    si_bits, _ = BlkBIST.bist_prbs_gen(poly=polynomial, inv=false,
                                       Nsym=nbits, seed=copy(seed))

    bist = StubBist(polynomial, false, copy(seed), falses(nbits), copy(si_bits),
                    false, initial_lock_cnt, threshold, 0, 0)

    BlkBIST.ber_check_prbs!(bist)

    @test bist.chk_lock_status
    @test bist.ber_err_cnt == 0
    @test bist.ber_bit_cnt == nbits - (threshold - initial_lock_cnt)
end

# --------------------------------------------------------------------------
# Full ber_checker_top! stub – needs param, Si, and extra fields
# --------------------------------------------------------------------------
mutable struct StubParam
    cur_blk::Int
    pam::Int
    bits_per_sym::Int
end

mutable struct FullStubBist
    param::StubParam
    polynomial::Vector{UInt8}
    inv::Bool
    gen_gray_map::Vector{UInt8}
    gen_en_precode::Bool
    gen_precode_prev_sym::Int
    chk_precode_prev_sym::Int
    chk_start_blk::Int
    chk_seed::Vector{Bool}
    chk_lock_status::Bool
    chk_lock_cnt::Int
    chk_lock_cnt_threshold::Int
    ber_err_cnt::Int
    ber_bit_cnt::Int
    ref_bits::Vector{Bool}
    So_bits::Vector{Bool}
    So::Vector{Float64}
    Si::Vector{UInt8}
    Si_bits::Vector{Bool}
end

function make_full_bist(;polynomial, order, threshold, chk_start_blk=1,
                         blk_size=1024, pam=2, bits_per_sym=1)
    param = StubParam(0, pam, bits_per_sym)
    FullStubBist(
        param,
        polynomial, false,
        UInt8[],        # gen_gray_map (empty → no gray)
        false, 0, 0,    # precode off
        chk_start_blk,
        zeros(Bool, order),
        false, 0, threshold,
        0, 0,
        falses(bits_per_sym * blk_size),
        falses(bits_per_sym * blk_size),
        zeros(blk_size),
        UInt8[],
        falses(bits_per_sym * blk_size),
    )
end

# --------------------------------------------------------------------------
# Variable-length input tests (frequency offset simulation)
# --------------------------------------------------------------------------

@testset "ber_check_prbs! – locked path with variable-length chunks" begin
    # Once locked, feed chunks of different sizes from a continuous PRBS.
    # All bits should compare correctly → zero BER.
    polynomial = UInt8[6, 7]
    seed = ones(Bool, 7)

    total_bits = 5000
    all_bits, _ = BlkBIST.bist_prbs_gen(poly=polynomial, inv=false,
                                        Nsym=total_bits, seed=copy(seed))

    # Start already locked with a seed that matches the TX seed
    bist = StubBist(polynomial, false, copy(seed),
                    falses(0), falses(0),
                    true, 0, 128, 0, 0)   # chk_lock_status = true

    # Simulate variable chunk sizes (as if freq offset causes ±6 syms / block)
    chunk_sizes = [1024, 1018, 1030, 1020, 908]
    @assert sum(chunk_sizes) == total_bits

    offset = 0
    for cs in chunk_sizes
        chunk = all_bits[offset+1:offset+cs]
        resize!(bist.Si_bits, cs)
        resize!(bist.ref_bits, cs)
        bist.Si_bits .= chunk
        BlkBIST.ber_check_prbs!(bist)
        offset += cs
    end

    @test bist.ber_err_cnt == 0
    @test bist.ber_bit_cnt == total_bits
end

@testset "ber_checker_top! – variable nsym across blocks (freq offset)" begin
    # End-to-end test through ber_checker_top! with PAM2.
    # Generates a continuous PRBS, splits into variable-size symbol blocks,
    # and checks zero BER after lock.
    polynomial = UInt8[6, 7]; order = 7
    tx_seed = ones(Bool, order)

    total_syms = 5000
    all_bits, _ = BlkBIST.bist_prbs_gen(poly=polynomial, inv=false,
                                        Nsym=total_syms, seed=copy(tx_seed))
    # For PAM2 with bits_per_sym=1, symbol value == bit value
    all_syms = UInt8.(all_bits)

    bist = make_full_bist(polynomial=polynomial, order=order, threshold=128,
                          chk_start_blk=1, blk_size=1024)

    # Variable chunk sizes simulating freq offset
    chunk_sizes = [1024, 1018, 1030, 1020, 908]
    @assert sum(chunk_sizes) == total_syms

    offset = 0
    for (blk, cs) in enumerate(chunk_sizes)
        bist.param.cur_blk = blk
        bist.Si = all_syms[offset+1:offset+cs]
        BlkBIST.ber_checker_top!(bist)
        offset += cs
    end

    @test bist.chk_lock_status
    @test bist.ber_err_cnt == 0
    @test bist.ber_bit_cnt > 0
    @test bist.ber_bit_cnt < total_syms   # some bits used for lock acquisition
end

@testset "ber_checker_top! – nsym == 0 block is harmless" begin
    polynomial = UInt8[6, 7]; order = 7
    bist = make_full_bist(polynomial=polynomial, order=order, threshold=128,
                          chk_start_blk=1)
    bist.param.cur_blk = 5
    bist.Si = UInt8[]   # empty – no RX symbols this block
    BlkBIST.ber_checker_top!(bist)
    @test bist.ber_bit_cnt == 0
    @test bist.ber_err_cnt == 0
end

@testset "ber_checker_top! – blocks before chk_start_blk are skipped" begin
    polynomial = UInt8[6, 7]; order = 7
    tx_seed = ones(Bool, order)
    bits, _ = BlkBIST.bist_prbs_gen(poly=polynomial, inv=false,
                                    Nsym=1024, seed=copy(tx_seed))

    bist = make_full_bist(polynomial=polynomial, order=order, threshold=128,
                          chk_start_blk=100)
    bist.param.cur_blk = 50   # before start
    bist.Si = UInt8.(bits)
    BlkBIST.ber_checker_top!(bist)
    @test bist.ber_bit_cnt == 0
    @test !bist.chk_lock_status
    @test isempty(bist.Si)     # Si should still be emptied
end
