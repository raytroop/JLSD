module BlkBIST
using UnPack, Random, Distributions

export pam_gen_top!, ber_checker_top!


function bist_prbs_gen(;poly, inv, Nsym, seed)
    seq = Vector{Bool}(undef,Nsym)
    for n = 1:Nsym
        seq[n] = inv
        for p in poly
            seq[n] ⊻= seed[p]
        end
        seed .= [seq[n]; seed[1:end-1]]
    end
    return seq, seed
end

function bist_prbs_gen!(seq; poly, inv, Nsym, seed)
    for n = 1:Nsym
        seq[n] = inv
        for p in poly
            seq[n] ⊻= seed[p]
        end
        seed .= [seq[n]; seed[1:end-1]]
    end
    return nothing
end


function pam_gen_top!(bist)
    @unpack pam, bits_per_sym, blk_size = bist.param
    @unpack polynomial, inv, gen_seed = bist
    @unpack gen_gray_map, gen_en_precode, gen_precode_prev_sym, So_bits, So = bist

    bist_prbs_gen!(So_bits, poly=polynomial, inv=inv,
                    Nsym=bits_per_sym*blk_size, seed=gen_seed)


    fill!(So, zero(Float64))
    for n = 1:bits_per_sym
        @. So = So + 2^(bits_per_sym-n)*So_bits[n:bits_per_sym:end]
    end

    #gray encoding
    if ~isempty(gen_gray_map)
        for n in 1:blk_size
            So[n] = gen_gray_map[So[n] + 1]
        end
    end

    if gen_en_precode
        for n = 1:blk_size
            So[n] = mod(So[n]-gen_precode_prev_sym , pam)
            gen_precode_prev_sym = So[n]
        end
        bist.gen_precode_prev_sym = gen_precode_prev_sym #need to write back
    end
    @. So = 2/(pam-1)*So - 1

    return nothing
end

function ber_checker_top!(bist)
    @unpack cur_blk, pam, bits_per_sym = bist.param
    @unpack gen_gray_map, chk_precode_prev_sym, chk_start_blk, Si = bist

    nsym = length(Si)

    # Skip until the warm-up window is over, or if the elastic-buffer
    # scheduler delivered no symbols this block (possible during a
    # transient stall when freq_offset_ppm is large).
    if cur_blk < chk_start_blk || nsym == 0
        empty!(Si)
        return
    end

    if bist.gen_en_precode
        bist.chk_precode_prev_sym = Si[end]
        Si .= mod.([chk_precode_prev_sym; Si[1:end-1]] .+ Si , pam)
    end

    if ~isempty(gen_gray_map)
        # Broadcasted getindex applies the lookup table to the variable-length block.
        Si .= getindex.(Ref(gen_gray_map), Si .+ 1)
    end

    # Resize working buffers to match the actual number of symbols this
    # block — under elastic-buffer scheduling the block size varies.
    nbits = bits_per_sym * nsym
    resize!(bist.Si_bits, nbits)
    resize!(bist.ref_bits, nbits)

    bist.Si_bits .= vec(stack(int2bits.(Si, bits_per_sym)))

    ber_check_prbs!(bist)

    # Drain so the next block starts fresh
    empty!(Si)
    return nothing
end

function ber_check_prbs!(bist)
    @unpack polynomial, inv, chk_seed, ref_bits, Si_bits = bist
    nbits_rcvd = lastindex(Si_bits)

    # err_loc = rand(Uniform(0,1.0), nbits_rcvd).< 1e-4;
    # Si_bits .= Si_bits .⊻ err_loc

    if bist.chk_lock_status
        bist_prbs_gen!(ref_bits, poly=polynomial, inv=inv,
                        Nsym=nbits_rcvd,seed=chk_seed)


        bist.ber_err_cnt += sum(Si_bits .⊻ ref_bits)
        bist.ber_bit_cnt += nbits_rcvd
    else
        for n = 1:nbits_rcvd
            brcv = Si_bits[n]
            btst = inv
			for p in polynomial
            	btst ⊻= chk_seed[p]
			end

            #need consecutive non-error for lock. reset when error happens
            bist.chk_lock_cnt = (btst == brcv) ? bist.chk_lock_cnt+1 : 0

            chk_seed .= [brcv; chk_seed[1:end-1]]

            if bist.chk_lock_cnt == bist.chk_lock_cnt_threshold
                bist.chk_lock_status = true
                println("prbs locked")
                # After lock, also score the bits remaining in this block
                # against the locked seed so they aren't lost (important
                # under elastic-buffer scheduling where block sizes vary).
                remaining = nbits_rcvd - n
                if remaining > 0
                    ref_bits_rem = @view ref_bits[1:remaining]
                    bist_prbs_gen!(ref_bits_rem, poly=polynomial, inv=inv,
                                   Nsym=remaining, seed=chk_seed)
                    bist.ber_err_cnt += sum((@view Si_bits[n+1:nbits_rcvd]) .⊻ ref_bits_rem)
                    bist.ber_bit_cnt += remaining
                end
                break
            end
        end
    end
end

function int2bits(num, nbit)
    return [Bool((num>>k)%2) for k in nbit-1:-1:0]
end



end
