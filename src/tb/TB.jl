using GLMakie, Makie
using UnPack, Random, Interpolations, ColorSchemes, Colors
include("../structs/TrxStruct.jl"); #using .TrxStruct
includet("../util/Util_JLSD.jl"); using .Util_JLSD
includet("../blks/BlkBIST.jl"); using .BlkBIST
includet("../blks/BlkTX.jl"); using .BlkTX
includet("../blks/BlkCH.jl"); using .BlkCH
includet("../blks/BlkRX.jl"); using .BlkRX
includet("../blks/BlkElasticBuffer.jl"); using .BlkElasticBuffer
includet("../blks/WvfmGen.jl"); using .WvfmGen



function init_trx()
    #global param for simulation
    param = TrxStruct.Param(
                data_rate = 56e9,
                pam = 2,
                osr = 24,
                blk_size = 2^10,
                subblk_size = 32,
                nsym_total = Int(1e6))
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

    init_plot(wvfm)

    println("init done")

    ebuf = TrxStruct.ElasticBuffer(
        param = param,
        ppm = param.freq_offset_ppm)

    return (;param, bist, drv, ch, clkgen, splr, dslc, eslc, cdr, adpt, wvfm, ebuf)
end

function sim_subblk(trx, blk_idx)
    @unpack param, bist, drv, ch, clkgen, splr = trx
    @unpack dslc, eslc, cdr, adpt, wvfm, ebuf = trx

    param.cur_subblk = blk_idx

    clkgen_pi_itp_top!(clkgen, pi_code=cdr.pi_code)

    sample_phi_top!(splr, clkgen.Φo_subblk)

    slicers_top!(dslc, splr.So_subblk, ref_code=[[128],[128],[128],[128]])
    slicers_top!(eslc, splr.So_subblk, ref_code=adpt.eslc_ref_vec)

    cdr_top!(cdr, dslc.So, eslc.So)

    adpt_top!(adpt, dslc.So, eslc.So)

    append!(bist.Si, [sum(dvec) for dvec in dslc.So])

    push!(wvfm.buffer12, adpt.eslc_ref_code)
    push!(wvfm.buffer22, cdr.pi_code)
end

function sim_blk(trx, blk_idx)
    @unpack param, bist, drv, ch, clkgen, splr = trx
    @unpack dslc, eslc, cdr, adpt, wvfm, ebuf = trx

    param.cur_blk = blk_idx

    pam_gen_top!(bist)

    dac_drv_top!(drv, bist.So)

    ch_top!(ch, drv.Vo)

    eb_write!(ebuf, ch.Vo)

    # Read blk_size_osr samples from elastic buffer to build interpolator
    blk_data = eb_read!(ebuf, param.blk_size_osr)

    t_start = (blk_idx - 1) * param.blk_size_osr
    sample_itp_top!(splr, blk_data, t_start=t_start)

    # Variable sub-block loop driven by RX rate
    freq_offset_ppm = param.freq_offset_ppm
    osr_rx = param.osr / (1 + freq_offset_ppm * 1e-6)
    subblk_size_osr_min = floor(Int, param.subblk_size * osr_rx * 0.5)

    subblk_idx = 0
    clkgen.Φo_start = t_start
    while true
        subblk_idx += 1
        next_start = clkgen.Φo_start + param.subblk_size * osr_rx
        if next_start > t_start + param.blk_size_osr + subblk_size_osr_min
            break
        end
        if subblk_idx > param.nsubblk + 2
            break  # safety guard
        end
        sim_subblk(trx, subblk_idx)
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




