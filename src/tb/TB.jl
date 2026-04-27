using GLMakie, Makie
using UnPack, Random, Interpolations, ColorSchemes, Colors
include("../structs/TrxStruct.jl"); #using .TrxStruct
includet("../util/Util_JLSD.jl"); using .Util_JLSD
includet("../blks/BlkBIST.jl"); using .BlkBIST
includet("../blks/BlkTX.jl"); using .BlkTX
includet("../blks/BlkCH.jl"); using .BlkCH
includet("../blks/BlkRX.jl"); using .BlkRX
includet("../blks/WvfmGen.jl"); using .WvfmGen



function init_trx()
    #global param for simulation
    param = TrxStruct.Param(
                data_rate = 56e9,
                pam = 2,
                osr = 24,
                blk_size = 2^10,
                subblk_size = 32,
                nsym_total = Int(1e6),
                freq_offset_ppm = 0.0)
    Random.seed!(param.rand_seed)

    #bist param
    bist = TrxStruct.Bist(
                param = param,
                polynomial = TrxStruct.PRBS31)


    #TX driver param
    drv = TrxStruct.Drv(
                param = param,
                ir = u_gen_ir_rc(param.dt, param.fbaud, 20*param.tui),
                fir = [1., -0.2],
                swing = 0.8,
                jitter_en = true,
                dcd = 0.03,
                rj_s = 300e-15,
                sj_amp_ui = 0.1,
                sj_freq = 2e5)

    #AWGN ch param
    ch = TrxStruct.Ch(
                param = param,
                ir_ch = u_fr_to_imp("./channel_data/TF_data/channel_4inch.mat",
                        param.tui, param.osr, npre = 20, npost= 79),
                ir_pad = u_gen_ir_rc(param.dt, param.fbaud, 20*param.tui),
                noise_en = true,
                noise_dbm_hz = -150 )

    #clkgen param
    clkgen = TrxStruct.Clkgen(
                param = param,
                nphases = 4,
                rj = .3e-12,
                skews = [0e-12, 0e-12, 0e-12, 0e-12])

    #sampler param
    splr = TrxStruct.Splr(
                param = param,
                ir = u_gen_ir_rc(param.dt, param.fbaud, 20*param.tui))

    #slicer param
    dslc = TrxStruct.Slicers(
                param = param,
                N_per_phi = ones(UInt8, clkgen.nphases),
                noise_rms = 1.5e-3,
                ofst_std = 7e-3,
                dac_min = -0.1,
                dac_max = 0.1)

    eslc = TrxStruct.Slicers(
                param = param,
                N_per_phi = UInt8.([1,0,0,0]),
                noise_rms = 1.5e-3,
                ofst_std = 7e-3,
                dac_min = 0,
                dac_max = 0.5)

    cdr = TrxStruct.Cdr(
                param = param,
                Neslc_per_phi = eslc.N_per_phi,
                kp = 1/2^6,
                ki = 1/2^14,
                pi_res = clkgen.pi_res)

    adpt = TrxStruct.Adpt(
                param = param,
                Neslc_per_phi = eslc.N_per_phi,
                mu_eslc = 1/64)

    #waveform plotting param
    wvfm = TrxStruct.Wvfm(
                param = param,
                en_plot = true,
                nrow = 3,
                ncol = 2)
    wvfm.eye1.clk_skews = clkgen.skews
    wvfm.eye1.clk_rj = clkgen.rj
    wvfm.eye1.noise_rms = dslc.noise_rms
    cm_turbo = colorschemes[:turbo]
    cm_turbo.colors[1] = RGB(1.0,1.0,1.0)
    wvfm.eye1.colormap = cm_turbo

    # Elastic buffer: decouples TX block production from RX sub-block
    # consumption so a non-zero freq_offset_ppm (TX↔RX clock-rate
    # mismatch) can be modelled.  At 0 ppm behaviour is identical (up
    # to floating-point round-off) to the previous fixed-stride loop.
    eb = TrxStruct.ElasticBuffer(
                param = param,
                t_tx_max = Int(splr.prev_nui/2*param.osr))

    init_plot(wvfm)

    println("init done")

    return (;param, bist, drv, ch, clkgen, splr, dslc, eslc, cdr, adpt, eb, wvfm)
end

function sim_subblk(trx, blk_idx)
    @unpack param, bist, drv, ch, clkgen, splr, eb = trx
    @unpack dslc, eslc, cdr, adpt, wvfm = trx

    param.cur_subblk = blk_idx

    # NOTE: caller has already generated Φo_subblk and confirmed readiness.

    sample_phi_top!(splr, eb, clkgen.Φo_subblk)

    slicers_top!(dslc, splr.So_subblk, ref_code=[[128],[128],[128],[128]])
    slicers_top!(eslc, splr.So_subblk, ref_code=adpt.eslc_ref_vec)

    cdr_top!(cdr, dslc.So, eslc.So)

    adpt_top!(adpt, dslc.So, eslc.So)

    append!(bist.Si, [sum(dvec) for dvec in dslc.So])

    push!(wvfm.buffer12, adpt.eslc_ref_code)
    push!(wvfm.buffer22, cdr.pi_code)
end

function sim_blk(trx, blk_idx)
    @unpack param, bist, drv, ch, clkgen, splr, eb = trx
    @unpack dslc, eslc, cdr, adpt, wvfm = trx

    param.cur_blk = blk_idx

    pam_gen_top!(bist)

    dac_drv_top!(drv, bist.So)

    # append!(drv.buffer_debug, mod.(u_find_0x(drv.Vo), param.osr) ./ param.osr)

    ch_top!(ch, drv.Vo)

    # Apply the RX bandwidth filter and push the filtered TX block into
    # the elastic buffer.  RX consumption is then driven by the dynamic
    # scheduling loop below: at +ppm we may run >nsubblk sub-blocks per
    # TX block; at -ppm we may run fewer (or zero on a given block) and
    # catch up on subsequent blocks.
    sample_filter_top!(splr, ch.Vo)
    eb_write!(eb, splr.Vo)

    subblk_count = 0
    while true
        # Generate candidate Φo_subblk including all timing perturbations
        # (Φ0, Φskew, Φrj).  Then check whether those exact times can
        # currently be served by the buffer.  Critically, we use the
        # same Φo_subblk for both the readiness check and the actual
        # sample so jitter is *not* re-randomised between the two.
        clkgen_pi_itp_top!(clkgen, eb, pi_code=cdr.pi_code)
        eb_can_read_times(eb, clkgen.Φo_subblk) || break

        # Commit: append history, run the sub-block, advance RX cursor.
        append!(clkgen.Φo, clkgen.Φo_subblk)
        sim_subblk(trx, subblk_count + 1)
        eb.t_rx += param.subblk_size * param.osr_rx
        subblk_count += 1
    end

    ber_checker_top!(bist)

    #record waveform here
    append!(wvfm.buffer11, drv.Vo)
    append!(wvfm.buffer21, ch.Vo)
    append!(wvfm.buffer31, splr.So)
    append!(wvfm.eye1.buffer, splr.Vo)
    wvfm.eslc_ref_ob.val = adpt.eslc_ref_code*eslc.dac_lsb+eslc.dac_min
    wvfm.eye1.x_ofst = 0#round(-(cdr.pi_code/clkgen.pi_codes_per_ui)*wvfm.eye1.x_npts_ui[]+wvfm.eye1.x_npts[]/2)
    w_plot_test(wvfm, cond=(param.cur_blk%wvfm.plot_every_nblk==0))
end


