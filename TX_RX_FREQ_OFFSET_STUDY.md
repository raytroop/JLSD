# TX/RX Frequency Offset Study Notes

This note explains what was added on top of Kevin Zheng's original JLSD
SerDes framework to model TX/RX clock frequency offset.  It is written as a
reading guide: what the original model assumed, what Claude Code added, what
ChatGPT corrected, and how the pieces fit together.

## 1. Original Model Assumption

The original JLSD flow was block aligned:

```text
TX/BIST -> driver -> channel -> sampler filter -> sample one fixed RX block
```

The important assumption was that TX and RX advanced in lockstep.  A TX block
with `blk_size` symbols produced exactly `blk_size` RX sample decisions.  The
old code could sample inside the current block with interpolation, but it did
not represent a receiver clock that slowly drifts relative to the transmitter
clock.

That works for `freq_offset_ppm = 0`, but it cannot model this SerDes reality:

```text
f_rx != f_tx
```

With a real TX/RX clock offset, the receiver eventually samples slightly more
or slightly fewer symbols than the transmitter produces over the same wall
time.

## 2. Core Frequency Offset Model

The added model uses the TX waveform grid as the reference time axis.

```text
f_rx = f_tx * (1 + ppm * 1e-6)
osr_rx = osr / (1 + ppm * 1e-6)
```

Sign convention:

```text
ppm > 0  -> RX clock faster -> osr_rx < osr
ppm < 0  -> RX clock slower -> osr_rx > osr
```

`osr` is the number of TX-grid samples per nominal UI.  `osr_rx` is the RX
sample spacing measured in TX-grid samples.

The RX sample times for one sub-block are:

```text
phi[j] = phi0 + t_rx_candidate + j * osr_rx + phi_skew[j] + phi_rj[j]
```

where `j = 0:(subblk_size - 1)`.

## 3. Why RxWindow Exists

Once RX samples are not block aligned, a sample for the next RX sub-block may
need waveform data from:

- the previous TX block,
- the current TX block,
- or a future TX block not written yet.

So the sampler cannot use only a block-local interpolation object.  The added
`RxWindow` is an analog waveform ring buffer:

```text
sample_filter_top!(splr, ch.Vo)
rxw_extend!(rxw, splr.Vo)
```

It stores filtered analog samples on the TX grid.  The RX sampler reads the
window at floating-point RX sample times with linear interpolation.

Main fields:

```text
rxw.buf      ring storage
rxw.t_min    oldest valid TX-grid sample index
rxw.t_max    exclusive write frontier
rxw.t_rx     committed RX read cursor
```

Important distinction:

```text
RxWindow is not a digital elastic buffer.
It stores oversampled analog waveform, not recovered symbols.
```

An overrun is therefore fatal.  Silently dropping analog samples would create a
step discontinuity in the waveform and invalidate CDR/BER results.

## 4. Dynamic Scheduler

The TX side still writes one TX block per outer loop.  After each write, the RX
side drains as many sub-blocks as are ready:

```julia
sample_filter_top!(splr, ch.Vo)
rxw_extend!(rxw, splr.Vo; phase_margin = rxw_phase_margin(clkgen))

while true
    clkgen_pi_itp_top!(clkgen, rxw; pi_code = cdr.pi_code)
    rxw_covers(rxw, clkgen.phi_subblk) || break
    append!(clkgen.phi_history, clkgen.phi_subblk)
    sim_subblk(...)
    clkgen_commit_subblk!(clkgen, rxw)
end
```

In the real code the variables use Greek names:

```text
clkgen.phi_subblk  -> clkgen.Phi-o_subblk in code
clkgen.phi_history -> clkgen.Phi-o in code
```

The scheduler is greedy: if the next RX sub-block is available, process it
immediately.  This keeps the window occupancy bounded for both positive and
negative ppm.

## 5. Claude Code's Main Additions

The Claude-generated commit is:

```text
42015f1 implement frequency offset by Claude Opus 4.7
```

Main additions from that work:

1. `Param.freq_offset_ppm` and `Param.osr_rx`

   File: `src/structs/TrxStruct.jl`

   These define the RX clock rate relative to the TX grid.

2. `RxWindow`

   File: `src/structs/TrxStruct.jl`

   A ring buffer for filtered analog waveform across TX block boundaries.

3. RxWindow helpers

   File: `src/blks/BlkRX.jl`

   ```julia
   rxw_extend!
   rxw_covers
   rxw_interp
   sample_filter_top!
   sample_phi_top!
   ```

4. Dynamic RX sub-block scheduler

   Files:

   ```text
   src/tb/TB.jl
   src/tb/Widget.jl
   ```

   The old fixed `nsubblk` loop was replaced by a readiness-driven loop.

5. BER checker changes

   File: `src/blks/BlkBIST.jl`

   The RX can now produce a variable number of symbols per TX block, so
   `Bist.Si` became a growable `Vector{UInt8}` and `ber_checker_top!` checks
   `length(Si)` instead of assuming `blk_size`.

6. Documentation and tests

   Files:

   ```text
   README.md
   freq_offset_formulas.md
   test/test_rx_window.jl
   test/test_ber_checker.jl
   ```

## 6. Problems Found In The Generated Version

The generated design was mostly correct conceptually, but there were subtle
state-machine issues.

### Problem A: Candidate Generation Mutated Committed State Too Early

The scheduler called `clkgen_pi_itp_top!` to generate the next sample-time
candidate, then called `rxw_covers` to check whether that candidate was ready.

The problem was that `clkgen_pi_itp_top!` updated committed state before the
readiness check:

```text
rxw.t_rx
clkgen.pi_code_prev
```

If `rxw_covers` returned false, the sub-block was not sampled, but the clock
state had already moved.  That is incorrect.

Correct rule:

```text
Do not commit RX cursor or PI state until the sub-block is actually sampled.
```

### Problem B: Rejected Candidate Could Regenerate RJ

If a candidate was not ready, the next TX block should retry the same future
sampling event.  The generated code could regenerate the candidate, including
new random jitter.

That means an event that did not happen yet could get a different random jitter
realization simply because the waveform window was not ready.  That is a
simulation artifact.

Correct rule:

```text
If a candidate is not ready, keep it pending and retry it unchanged.
```

### Problem C: RxWindow Reclaimed Too Aggressively

The first implementation reclaimed data below:

```text
floor(t_rx)
```

But a valid next sample can be earlier than `t_rx` because of negative PI,
skew, or random jitter.  Reclaiming exactly at `t_rx` can produce unnecessary
`rxw_covers == false` failures.

Correct rule:

```text
Keep a low-side phase margin below t_rx.
```

## 7. ChatGPT Fixes

The current working tree contains fixes for the above issues.

### 7.1 Candidate/Commit State Machine

File: `src/blks/BlkRX.jl`

`clkgen_pi_itp_top!` now only prepares a pending candidate:

```julia
clkgen.t_rx_subblk = t_rx_candidate
clkgen.t_rx_next = t_rx_candidate + subblk_size * osr_rx
clkgen.pi_code_subblk = pi_code
clkgen.Phi-o_subblk_valid = true
```

It does not commit `rxw.t_rx` or `clkgen.pi_code_prev`.

The commit happens only here:

```julia
clkgen_commit_subblk!(clkgen, rxw)
```

That function runs only after:

```text
rxw_covers == true
sample_phi_top! has sampled the waveform
slicers/CDR/adaptation have processed the sub-block
```

### 7.2 Pending Candidate Retried Unchanged

At the top of `clkgen_pi_itp_top!`:

```julia
clkgen.Phi-o_subblk_valid && return nothing
```

So if the previous candidate was not ready, the next scheduler iteration keeps
the same `Phi-o_subblk` instead of drawing new RJ.

### 7.3 Low-Side Phase Margin

File: `src/blks/BlkRX.jl`

`rxw_extend!` now accepts:

```julia
phase_margin
```

The production scheduler passes:

```julia
rxw_phase_margin(clkgen)
```

That margin covers:

```text
PI phase range
deterministic clock skew
finite RJ envelope
```

The reclaim rule is now:

```text
t_min = max(t_min, floor(t_rx - phase_margin))
```

instead of:

```text
t_min = max(t_min, floor(t_rx))
```

## 8. Important Code Reading Order

Read in this order.

### Step 1: Parameters and State

File: `src/structs/TrxStruct.jl`

Look at:

```julia
Param.freq_offset_ppm
Param.osr_rx
Clkgen.Phi-o_subblk_valid
Clkgen.t_rx_subblk
Clkgen.t_rx_next
Clkgen.pi_code_subblk
RxWindow
```

This tells you what persistent state exists.

### Step 2: RxWindow Helpers

File: `src/blks/BlkRX.jl`

Read:

```julia
rxw_extend!
rxw_phase_margin
rxw_covers
rxw_interp
```

These functions define the waveform window contract.

Key readiness rule:

```text
floor(min(phi)) >= t_min
floor(max(phi)) + 1 < t_max
```

The `+1` is required because linear interpolation reads both `k0` and
`k0 + 1`.

### Step 3: RX Clock Candidate Generation

File: `src/blks/BlkRX.jl`

Read:

```julia
clkgen_pi_itp_top!
clkgen_commit_subblk!
```

This is the most important part of the fix.

### Step 4: Main Scheduler

File: `src/tb/TB.jl`

Read:

```julia
sim_blk
sim_subblk
```

This shows the order:

```text
write waveform
prepare candidate
check readiness
sample
run slicers/CDR/adaptation
commit candidate
repeat
```

### Step 5: Widget Scheduler

File: `src/tb/Widget.jl`

Read:

```julia
step_sim_blk
step_sim_subblk
```

It mirrors the batch scheduler.

### Step 6: Variable-Length BER

File: `src/blks/BlkBIST.jl`

Read:

```julia
ber_checker_top!
ber_check_prbs!
```

The receiver no longer produces exactly `blk_size` decisions each TX block,
so the BER checker must handle variable `length(Si)`.

## 9. Formula Sheet Map

Read `freq_offset_formulas.md` alongside the code.

Useful sections:

```text
Section 1:  ppm -> osr_rx relation
Section 3:  full sample-time model
Section 4:  rxw_covers readiness check
Section 5:  linear interpolation
Section 8:  why greedy scheduling bounds drift
Section 9:  capacity derivation
Section 13: candidate/commit state machine
Section 14: reclaim margin
Section 16: CDR ki sizing
```

The most important conceptual point:

```text
The average frequency drift is modeled by osr_rx-spaced t_rx movement.
The PI code models bounded fine phase correction.
Do not also accumulate pi_wrap_ui inside phi0, or drift is double-counted.
```

## 10. Why pi_wrap_ui Was Removed From phi0

In the original block-aligned model, `pi_wrap_ui` helped represent accumulated
phase wrapping.

In the frequency-offset model, the long-term clock-rate mismatch is already
represented by:

```text
t_rx advances by osr_rx-spaced RX samples
```

If `pi_wrap_ui` is also included in `phi0`, the same long-term drift is counted
twice.

The current model handles PI wrap like this:

```text
1. Compute candidate t_rx.
2. If pi_code wraps, shift candidate t_rx by pi_ui_cover * osr_rx.
3. Keep phi0 bounded inside one PI coverage range.
4. Commit only after sampling.
```

## 11. Tests Added For Confidence

File: `test/test_rx_window.jl`

The early tests use a small local stub to validate the scheduler math without
loading GLMakie.

The later production regression tests import the real `BlkRX.jl` and check:

```text
candidate is not committed on underrun
pending candidate is retried unchanged
rxw_extend keeps low-side phase margin
```

File: `test/test_ber_checker.jl`

Checks BER lock-instant accounting so trailing bits after PRBS lock are counted
instead of being thrown away.

## 12. Verification Performed

Commands run after the fixes:

```bash
julia --project=. test/test_rx_window.jl
julia --project=. test/test_ber_checker.jl
```

Both passed.

A full default simulation was also run:

```text
nblk = 977
freq_offset_ppm = 100
checked bits = 898,913
BER = 0.0
```

## 13. What To Watch If You Extend This

1. Do not mutate committed RX clock state before `rxw_covers` succeeds.

2. Do not regenerate random jitter for a pending candidate.

3. If you add new timing perturbations, include them in:

   ```text
   rxw_covers candidate sample times
   rxw_phase_margin
   capacity derivation
   documentation
   ```

4. If you change `subblk_size`, `pi_ui_cover`, skew, or RJ envelope, revisit
   `RxWindow.capacity`.  The default `4 * blk_size_osr` is intentionally
   conservative for the current configuration.

5. If you change CDR gains, rerun full BER simulations.  Frequency-offset
   support requires enough CDR integral range to track the ppm.

## 14. One-Sentence Summary

The implementation models SerDes TX/RX frequency offset by keeping the analog
TX waveform on a fixed TX sample grid, sampling it through `RxWindow` at an RX
grid spaced by `osr_rx`, and using a candidate/commit scheduler so RX clock
state advances only when the corresponding analog samples have actually been
processed.
