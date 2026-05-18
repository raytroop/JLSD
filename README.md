# JLSD - Julia SerDes
This repository hosts Pluto notebooks and relevant source codes for a personal project - using Julia to build fast and lightweight SerDes models and simulation framework

The goals for this project (still in its infant stage) are
- Document and share my personal Julia journey so far
- Using SerDes simulation as an example to demonstrate the pros/cons of Julia
- Begin a open-source and expandable SerDes simulation framework that encourages academia and industry adoption to evaluate more sophositicated architectures and algorithms
- Expose my fellow analog/mixed-signal designers to Julia (because not everyone can get a MATLAB license). It's much easier to design circuits when one can play with the specifications instead of taking it from others at face value

## Simulation framework
The modeling framework is based on custom structs and block simulations. The code is not heavily documentated yet, but should be self-explanatory. Go through the Pluto notebooks to understand the key concepts in the models.

## Notebooks
In the Pluto Notebooks directory, you will find the .jl files for the notebooks to be viewed and played with on your local machine. .html and .pdf files are also included in the directory.

## Standalone widget
Currently, there is a demo widget powered by Makie (shown below. Simply run `Main_UI.jl` and start playing with a basic SerDes model's parameters. The model consists of a relatively detailed transmitter, a low-loss channel, and a simple receiver w/ baud-rate CDR. Note that the widget might be continuously updated to include more (common) features. Use this as an example to extend to your own models.
![widget_ui](https://github.com/kevjzheng/JLSD/blob/main/img/widget_ui.png)

---

## TX/RX frequency-offset modelling (`RxWindow`)

This fork adds support for an arbitrary frequency offset between the TX
and RX clocks (`freq_offset_ppm`).  The original framework assumed
`blk_size` is an integer multiple of `subblk_size` and that TX and RX
advance in lockstep, foreclosing any rate mismatch.  The note in
`JLSD_pt2_framework.jl` calls this out as a planned extension; the
implementation is now in this fork.

For the underlying math (sample-time model, readiness check, capacity
derivation, etc.) see [freq_offset_formulas.md](./freq_offset_formulas.md).

### Sign convention

$$f_{rx} = f_{tx}\,(1 + \text{ppm}\cdot 10^{-6}),\qquad \text{osr}\_\text{rx} = \frac{\text{osr}}{1 + \text{ppm}\cdot 10^{-6}}$$

- `freq_offset_ppm > 0` → RX clock is faster than TX
- `freq_offset_ppm < 0` → RX clock is slower than TX
- `osr_rx` is the inter-RX-sample spacing on the TX simulation grid

### Architecture: from block-aligned to RX-clock-aligned

The TX side and the channel still produce one block of waveform per
outer iteration on the TX grid.  An **`RxWindow`** sits between the
sampler bandwidth filter and the RX sampler, holding the analog
waveform across TX block boundaries so the RX can interpolate at
arbitrarily-spaced RX-clock instants.

```
   TX-grid (fixed osr)                                RX-grid (osr_rx)
   ──────────────────                                 ────────────────
   bist → drv → ch → splr (filter) ──► ┌───────────┐ ──► sampler.interp
                                       │ RxWindow  │     ↓
                                       │  ring buf │     slicers → CDR
                                       │  t_min,   │
                                       │  t_max,   │     t_rx advances by
                                       │  t_rx     │     subblk_size·osr_rx
                                       └───────────┘     per sub-block
```

The outer block loop is replaced by a **dynamic sub-block scheduler**:
each TX block, the scheduler drains as many sub-blocks as fit, then
hands control back to TX.  Under `freq_offset_ppm` ≠ 0 the per-block count
oscillates around `nsubblk·(1+ppm·10⁻⁶)` (e.g. 31 / 32 / 33 for
`subblk_size=32`).

### `RxWindow` — what it is and what it isn't

`RxWindow` is **not** a digital elastic buffer.  Real-link elastic
buffers sit after deserialization, hold recovered symbols, and are
drained by SKP ordered sets (PCIe SKP, USB SOF, etc.).  This struct sits
at the **analog sampling layer**, holds **oversampled waveform**, and
its only role is to let the RX sampler interpolate at the actual
RX-clock instants which drift relative to the TX simulation grid.

- `buf::Vector{Float64}` — ring storage, indexed `(t mod C) + 1`
- `t_min::Int`  — oldest sample still available
- `t_max::Int`  — exclusive write frontier (TX-grid samples written so far)
- `t_rx::Float64` — RX read cursor, persists across TX block iterations
- `capacity::Int` — `4 · blk_size_osr` by default (see formula sheet §10)

### Underrun is benign, overrun is fatal

These are **asymmetric** for an analog window — unlike a digital FIFO:

| Condition | What it means | Handling |
|---|---|---|
| **Underrun** | `rxw_covers(rxw, Φi)` returns `false` — the next sub-block would read past `t_max` (or below `t_min`). | The inner scheduler loop simply `break`s.  The sub-block is processed when the next TX block arrives.  **No data is lost**, no counter is needed. |
| **Overrun** | About-to-be-written samples would overwrite ring positions that the RX has not yet consumed. | `rxw_extend!` raises a **hard `error()`** with diagnostics.  Silent overwrite would inject a step discontinuity into the analog signal and invalidate downstream CDR/BER results. |

The original silent-drop-and-count design was changed because for an
analog signal an overrun is not a metric to track — it's a hard
invariant violation.

### Why `pi_wrap_ui` is dropped from `Φ0`

In the legacy block-aligned scheduler, `Φ0` carried the full CDR phase
correction:

$$\Phi_0^{\text{legacy}} = \text{osr}\cdot\bigl(\text{pi}\_\text{wrap}\_\text{ui} + \text{pi}\_\text{code}/\text{pi}\_\text{codes}\_\text{per}\_\text{ui}\bigr)$$

`pi_wrap_ui` accumulated the CDR's cumulative ppm-drift compensation
across PI wraps.

Under the new $t_{rx}$-based scheduler, the ppm drift is **already
tracked by $t_{rx}$ advancing at `osr_rx` per RX sample**.  Including
`pi_wrap_ui` in $\Phi_0$ would double-count this drift; for non-zero
ppm, $\Phi_0$ would grow without bound and eventually throttle the
scheduler into overrun.

The fix: when `pi_code` wraps (`|Δpi_code| > pi_wrap_ui_Δcode`), absorb
the `pi_ui_cover` discontinuity into `t_rx` instead of `pi_wrap_ui`:

$$t_{rx} \leftarrow t_{rx} - \text{sign}(\Delta\text{pi}\_\text{code})\cdot\text{pi}\_\text{ui}\_\text{cover}\cdot\text{osr}\_\text{rx}$$

$$\Phi_0 = \text{osr}\_\text{rx}\cdot\frac{\text{pi}\_\text{code}+\text{pi}\_\text{nonlin}\_\text{lut}[\text{pi}\_\text{code}+1]}{\text{pi}\_\text{codes}\_\text{per}\_\text{ui}}$$

This keeps absolute sample positions continuous while bounding $\Phi_0$
in $[0,\;\text{pi}\_\text{ui}\_\text{cover}\cdot\text{osr}\_\text{rx})$.  The CDR's fine
within-PI correction still works; the cumulative drift is tracked
exactly once.

### CDR loop bandwidth and ppm tracking

A consequence of the `pi_wrap_ui → t_rx` absorption: under the new
scheduler, `pi_code` is the **only** mechanism the CDR has for fine
drift tracking (the legacy unbounded `pi_wrap_ui` accumulator is gone).
The integrator gain `ki` therefore directly bounds the maximum
trackable ppm.

The integrator-equilibrium drift rate is approximately

$$\text{drift rate}\;[\text{pi}\_\text{codes / sub-block}] \;\approx\; \text{ki}\cdot N_{\text{votes/sub-block}}\cdot \text{vote bias}$$

Required drift for ppm tracking:

$$\text{required}\;=\;\frac{\text{ppm}\cdot 10^{-6}\cdot \text{osr}\cdot \text{subblk}\_\text{size}}{\text{osr}\_\text{rx}}\cdot\text{pi}\_\text{codes}\_\text{per}\_\text{ui}\;\approx\;\text{ppm}\cdot 10^{-6}\cdot 32\cdot 64$$

For `subblk_size=32`, `pi_codes_per_ui=64`: required ≈ `2.05 × ppm × 10⁻³` `pi_codes`/sub-block.

Practical tracking range (empirical, default kp=1/2⁶, full
`nblk = 977` run):

| `ki` | Trackable \|ppm\| | BER @ ppm = 100 |
|---|---|---|
| `1/2¹⁴` (original) | $\le 30$ | $0.49$ (no lock) |
| `1/2¹²` (**default in this fork**) | $\le 100$ | $\sim 10^{-3}$ |
| `1/2¹⁰` | $\le 300$ | — |

The original `ki = 1/2¹⁴` was tuned for the ppm-free simulator.  This
fork raises `init_trx()`'s default to `ki = 1/2¹²` so the default
`freq_offset_ppm = 100` simulation gives BER ≈ 0 out of the box.
Users running the Widget at higher \|ppm\| can crank `ki` further via
direct edit (no slider yet).

### BER checker — variable-length input + lock-instant accounting

Two BER-checker changes were needed:

1. **Variable RX symbol count per TX block.**  With `freq_offset_ppm ≠ 0`
   each TX block produces 31 / 32 / 33 RX symbols (etc.) instead of a
   fixed `blk_size`.  `Bist.Si` is now a growable `Vector{UInt8}` (was a
   fixed-capacity `CircularBuffer` that silently dropped tail symbols).
   `ber_checker_top!` reads `nsym = length(Si)` and resizes working
   buffers each call; `empty!(Si)` at the end so the next block starts
   fresh.

2. **Trailing-bit accounting at the lock instant.**  The pre-fix
   `ber_check_prbs!` advanced the seed past the just-locked bit but
   then `break`ed without counting the trailing bits in `ber_err_cnt` /
   `ber_bit_cnt`.  Fixed to compare the remaining `nbits_rcvd - n` bits
   against a freshly-generated reference and update both counters
   before breaking.

### Files changed

| File | Change |
|---|---|
| `src/structs/TrxStruct.jl` | `freq_offset_ppm`, `osr_rx` (mutable) in `Param`; new `RxWindow` struct; `Bist.Si` → `Vector{UInt8}` |
| `src/blks/BlkRX.jl` | `rxw_extend!`, `rxw_covers`, `rxw_interp`, `sample_filter_top!`; `clkgen_pi_itp_top!` rewired (wrap absorbed into `t_rx`); `sample_phi_top!` reads from `rxw` |
| `src/blks/BlkBIST.jl` | Variable-`nsym` handling in `ber_checker_top!`; trailing-bit accounting in `ber_check_prbs!` |
| `src/tb/TB.jl` | `freq_offset_ppm = 100.0` default; `cdr.ki` bumped `1/2^14 → 1/2^12` so the CDR can track that ppm; dynamic sub-block scheduler in `sim_blk`; constructs `rxw` |
| `src/tb/Widget.jl` | "RX freq offset" slider (`-500:10:500` ppm); slider callback keeps `osr_rx` consistent; matching dynamic scheduler in `step_sim_blk` |
| `src/Main_UI.jl` | Constructs `rxw`, adds to NamedTuple |
| `test/test_rx_window.jl` | New; 9 standalone test sets (no GLMakie dep) |
| `test/test_ber_checker.jl` | New; lock-instant accounting regression |
| `freq_offset_formulas.md` | New; full derivation of sample-time model, readiness check, capacity bound |

### Validation

End-to-end BER over the full `nblk = 977` TX-block sim (defaults except
where noted):

| Scenario | `kp` | `ki` | ppm | BER |
|---|---|---|---|---|
| **`init_trx()` default in this fork** | $1/2^{6}$ | $1/2^{12}$ | $100$ | $\sim 10^{-3}$ |
| Same gains | $1/2^{6}$ | $1/2^{12}$ | $\le \|100\|$ | $\sim 10^{-3}$ |
| Cranked CDR | $1/2^{4}$ | $1/2^{10}$ | $\le \|300\|$ | $\sim 10^{-3}$ |
| Cranked CDR | $1/2^{4}$ | $1/2^{10}$ | $\pm 500$ | $\sim 0.1$ (CDR-loop bandwidth limit) |
| Legacy `ki` (ppm-free regime) | $1/2^{6}$ | $1/2^{14}$ | $0$ | $0$ |
| Legacy `ki` (would be wrong) | $1/2^{6}$ | $1/2^{14}$ | $100$ | $0.49$ (CDR can't track) |

No `RxWindow` overrun fires anywhere in the slider range ±500 ppm — the
hard-error guard exists for code regressions, not normal operation.

All 19 unit-test assertions in `test/` pass.

### Running

```bash
# Unit tests (no GLMakie dependency)
julia --project=. test/test_rx_window.jl
julia --project=. test/test_ber_checker.jl

# End-to-end widget (drag the "RX freq offset" slider)
julia --project=. src/Main_UI.jl
```
