"""
Regression test for ber_check_prbs! lock-instant accounting.

The pre-fix code threw away the bits between the lock instant and the end
of the block — `bist_prbs_gen` was called only to advance the seed, and
nothing was added to ber_err_cnt / ber_bit_cnt.

This test feeds a PRBS sequence, configures lock to fire mid-block, and
checks that the trailing bits ARE counted in ber_bit_cnt.

Run from the repo root:
    julia test/test_ber_checker.jl
"""

using Test

# Stub a minimal load path so BlkBIST can be included without pulling in
# the full GLMakie stack.
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
    nbits     = 10
    threshold = 3
    initial_lock_cnt = 1

    si_bits, _ = BlkBIST.bist_prbs_gen(poly = polynomial, inv = false,
                                       Nsym = nbits, seed = copy(seed))

    bist = StubBist(polynomial, false, copy(seed), falses(nbits), copy(si_bits),
                    false, initial_lock_cnt, threshold, 0, 0)

    BlkBIST.ber_check_prbs!(bist)

    @test bist.chk_lock_status
    @test bist.ber_err_cnt == 0
    @test bist.ber_bit_cnt == nbits - (threshold - initial_lock_cnt)
end

println("BER checker test passed.")
