# Frequency Offset & RxWindow — Formula Sheet

Reference notation for the TX/RX frequency-offset design (analog waveform
window + dynamic sub-block scheduling).  All time variables are measured in
**TX-grid sample-index units** unless noted.

---

## 1. Frequency offset relation

$$f_{rx} = f_{tx}\,(1 + \text{ppm}\cdot 10^{-6})$$

$$\text{osr\_rx} = \frac{\text{osr}}{1 + \text{ppm}\cdot 10^{-6}}$$

$$\text{osr\_rx}_{\max} = \frac{\text{osr}}{1 - |\text{ppm}|_{\max}\cdot 10^{-6}} \approx \text{osr}\,(1 + |\text{ppm}|_{\max}\cdot 10^{-6})$$

Sign convention: $\text{ppm} > 0$ means RX clock is faster than TX
($\text{osr\_rx} < \text{osr}$).

**Greedy scheduler note.** Under the dynamic scheduler used here, the
number of sub-blocks consumed per outer iteration

$$k_n = \left\lfloor \frac{(t_{\max}-t_{rx})_{\text{after write}}}{\text{subblk\_size}\cdot\text{osr\_rx}} \right\rfloor$$

is **not constant** — it adapts every iteration. The §8 proof then bounds
$t_{\max}-t_{rx}$ by one sub-block **for both signs of ppm**; neither
direction accumulates and neither produces a sustained risk.

The "ppm $\to$ failure" mapping below applies only to a hypothetical
**static** schedule with $k_n \equiv \text{nsubblk}$. The §15 hard-error
overrun guard exists to catch that regression, not to express a current
risk.

| Static schedule | per-block RX consumption | per-block TX production | failure mode |
|---|---|---|---|
| $\text{ppm} > 0$ | $\dfrac{\text{blk\_size\_osr}}{1+\text{ppm}\cdot 10^{-6}} < \text{blk\_size\_osr}$ | $\text{blk\_size\_osr}$ | **overrun** (writer outpaces reader) |
| $\text{ppm} < 0$ | $> \text{blk\_size\_osr}$ | $\text{blk\_size\_osr}$ | **underrun** (reader outpaces writer) |

---

## 2. Sub-block nominal sample times

Ignoring $\Phi_0$, $\Phi_{\text{skew}}$, $\Phi_{rj}$:

$$\Phi_i[j] = t_{rx} + (j-1)\cdot\text{osr\_rx}, \qquad j = 1, 2, \ldots, \text{subblk\_size}$$

$$\min_j \Phi_i[j] = t_{rx}$$

$$\max_j \Phi_i[j] = t_{rx} + (\text{subblk\_size} - 1)\cdot\text{osr\_rx}$$

A sub-block contains $\text{subblk\_size}$ **samples** spanning
$\text{subblk\_size}-1$ intervals — hence the $(j-1)$ multiplier.

---

## 3. Full sample-time model (with perturbations)

$$\Phi_i[j] = \underbrace{\Phi_0}_{\text{CDR/PI}} + \underbrace{t_{rx,\text{cand}} + (j-1)\cdot\text{osr\_rx}}_{\text{nominal}} + \underbrace{\Phi_{\text{skew}}[j]}_{\text{clock skew}} + \underbrace{\Phi_{rj}[j]}_{\text{random jitter}}$$

$$\Phi_0 = \text{osr\_rx}\cdot\frac{\text{pi\_code} + \text{pi\_nonlin\_lut}[\text{pi\_code}+1]}{\text{pi\_codes\_per\_ui}}$$

The legacy `pi_wrap_ui` term is **omitted** here.  In the legacy
block-aligned framework, `pi_wrap_ui` absorbed the CDR's cumulative
phase drift; including it in $\Phi_0$ under the new $t_{rx}$-based
scheduler would double-count the ppm drift (already tracked by
$\text{osr\_rx}$-spaced $t_{rx}$ advancement) and inflate $\Phi_0$ until
the greedy scheduler is throttled into overrun.  Instead, the scheduler
first prepares a **candidate** cursor.  When the CDR's $\text{pi\_code}$
wraps ($|\Delta\text{pi\_code}| > \text{pi\_wrap\_ui\_}\Delta\text{code}$),
the discontinuity is absorbed into that candidate cursor:

$$t_{rx,\text{cand}} \;=\; t_{rx} \;-\; \text{sign}(\Delta\text{pi\_code})\cdot\text{pi\_ui\_cover}\cdot\text{osr\_rx}$$

If no wrap occurs, $t_{rx,\text{cand}} = t_{rx}$.

This keeps absolute sample positions continuous while leaving $\Phi_0$
bounded in $[0,\;\text{pi\_ui\_cover}\cdot\text{osr\_rx})$.  The candidate
is committed only after `rxw_covers` succeeds and the sub-block is sampled;
otherwise it remains pending and is retried unchanged after the next TX
write.  This avoids committing PI-wrap state or regenerating RJ for a
sub-block that has not actually occurred.

$$\Phi_{\text{skew}}[j] = \frac{\text{skews}\bigl[((j-1)\bmod n_{\text{phases}})+1\bigr]}{t_{ui}}\cdot\text{osr}$$

$$\Phi_{rj}[j] = \frac{r_j}{t_{ui}}\cdot\text{osr}\cdot \mathcal{N}(0,1)$$

$\Phi_0$ is scaled by $\text{osr\_rx}$ because the PI covers
$\text{pi\_ui\_cover}$ **RX-clock UI**.  Skew and jitter remain in absolute
TX-grid sample units (converted via $\text{osr}/t_{ui}$).

---

## 4. RxWindow readiness check (`rxw_covers`)

$$t_{\text{need,min}} = \lfloor \min_j \Phi_i[j] \rfloor$$

$$t_{\text{need,max}} = \lfloor \max_j \Phi_i[j] \rfloor + 1$$

$$\text{covers} \iff (t_{\text{need,min}} \ge t_{\min}) \wedge (t_{\text{need,max}} < t_{\max})$$

The $+1$ on $t_{\text{need,max}}$ accounts for the $k{+}1$ tap consumed by
linear interpolation.

---

## 5. Linear interpolation kernel (`rxw_interp`)

$$k_0 = \lfloor t \rfloor, \qquad \alpha = t - k_0$$

$$\text{rxw\_interp}(t) = (1-\alpha)\cdot\text{buf}[k_0 \bmod C + 1] + \alpha\cdot\text{buf}[(k_0{+}1) \bmod C + 1]$$

where $C = \text{capacity}$.

---

## 6. Nominal inner-loop exit condition

$$\text{available} = t_{\max} - t_{rx}$$

Ignoring PI/skew/RJ and assuming no pending candidate shift, the greedy
scheduler continues while `rxw_covers` is true:

$$t_{\max} > \lfloor t_{rx} + (\text{subblk\_size} - 1)\cdot\text{osr\_rx} \rfloor + 1$$

Approximating $\lfloor x\rfloor \approx x$, the loop exits when

$$\text{available} \lesssim (\text{subblk\_size} - 1)\cdot\text{osr\_rx} + 1$$

Difference vs. the looser bound $\text{subblk\_size}\cdot\text{osr\_rx}$:

$$\text{subblk\_size}\cdot\text{osr\_rx} - \bigl[(\text{subblk\_size} - 1)\cdot\text{osr\_rx} + 1\bigr] = \text{osr\_rx} - 1$$

Using the loose bound would leave nearly one RX UI unused per iteration —
slight under-utilisation, not a correctness bug.

---

## 7. Per-block sub-block count

$$k = \left\lfloor \frac{\text{blk\_size\_osr}}{\text{subblk\_size}\cdot\text{osr\_rx}} \right\rfloor \approx \text{nsubblk}\,(1 + \text{ppm}\cdot 10^{-6})$$

where $\text{nsubblk} = \dfrac{\text{blk\_size}}{\text{subblk\_size}}$.

For $\text{ppm}>0$: $k$ tends to be larger than nominal $\text{nsubblk}$.
For $\text{ppm}<0$: $k$ tends to be smaller.

---

## 8. Steady-state occupancy analysis

After $n$ outer iterations, define

$$m = n\cdot\text{nsubblk}\cdot(1 + \text{ppm}\cdot 10^{-6}), \qquad r = m - \lfloor m\rfloor \in [0, 1)$$

$$t_{\max} = n\cdot\text{blk\_size\_osr}$$

$$t_{rx} = \lfloor m\rfloor\cdot\text{subblk\_size}\cdot\text{osr\_rx} = (m - r)\cdot\text{subblk\_size}\cdot\text{osr\_rx}$$

Substituting:

$$
\begin{aligned}
t_{\max} - t_{rx}
&= n\cdot\text{blk\_size\_osr} - (m - r)\cdot\text{subblk\_size}\cdot\text{osr\_rx} \\
&= n\cdot\text{blk\_size\_osr} - n\cdot\text{nsubblk}(1+\text{ppm}\cdot 10^{-6})\cdot\text{subblk\_size}\cdot\text{osr\_rx} + r\cdot\text{subblk\_size}\cdot\text{osr\_rx} \\
&= n\cdot\text{blk\_size\_osr} - n\cdot\text{blk\_size\_osr} + r\cdot\text{subblk\_size}\cdot\text{osr\_rx} \\
&= r\cdot\text{subblk\_size}\cdot\text{osr\_rx}
\end{aligned}
$$

Used the identity
$\text{nsubblk}(1+\text{ppm}\cdot 10^{-6})\cdot\text{subblk\_size}\cdot\text{osr\_rx} = \text{blk\_size\_osr}$.

Therefore the steady-state pre-write occupancy is bounded:

$$t_{\max} - t_{rx} < \text{subblk\_size}\cdot\text{osr\_rx}, \qquad \forall\,n$$

The fractional residue $r$ oscillates within $[0,1)$ but **does not
accumulate** — the greedy scheduler self-stabilises.

---

## 9. Maximum occupancy — derivation of `capacity`

The held window is $[t_{\min},\,t_{\max})$.  Its width just after a write is
bounded by the sum of three terms:

$$\text{occupancy}_{\max} \;=\; \underbrace{\text{blk\_size\_osr}}_{\text{(A) the write just done}} \;+\; \underbrace{(t_{\max}-t_{rx})_{\text{pre-write}}}_{\text{(B) post-drain residue}} \;+\; \underbrace{(t_{rx}-t_{\min})_{\text{reclaim slack}}}_{\text{(C) negative-phase reserve}}$$

### Term-by-term

**(A) Write — exact.**
$$(A) \;=\; \text{blk\_size\_osr}$$

**(B) Residue left by the inner-loop exit — bounded by one "fat" sub-block.**

The inner loop exits when the next sub-block fails `rxw_covers`, i.e. when

$$t_{\max} \;\le\; \lfloor t_{rx} + (\text{subblk\_size}-1)\cdot\text{osr\_rx} + \Phi_0^{\max} + \Phi_{\text{skew}}^{\max} + \Phi_{rj}^{\max}\rfloor + 1$$

so

$$(B) \;\le\; (\text{subblk\_size}-1)\cdot\text{osr\_rx} \;+\; \Phi_0^{\max} \;+\; \Phi_{\text{skew}}^{\max} \;+\; \Phi_{rj}^{\max} \;+\; 1$$

with the per-perturbation envelopes:

$$\Phi_0^{\max} \;=\; \text{pi\_ui\_cover}\cdot\text{osr\_rx}_{\max}$$

$$\Phi_{\text{skew}}^{\max} \;=\; \frac{\max|\text{skews}|}{t_{ui}}\cdot\text{osr}$$

$$\Phi_{rj}^{\max} \;=\; n_\sigma\cdot\frac{r_j}{t_{ui}}\cdot\text{osr}$$

Gaussian RJ is unbounded, so $n_\sigma$ is an engineering envelope.  The
implementation's reclaim helper uses $n_\sigma=8$ by default.

**(C) Negative-phase reserve — the data BELOW $t_{rx}$ that the next
sub-block could still need.**

If reclaim used $t_{\min}\!\leftarrow\!\max(t_{\min},\,\lfloor t_{rx}\rfloor)$,
this term would be **0**, and a sub-block with sufficiently negative
phase could fail `rxw_covers` even with infinite capacity.  The production
scheduler therefore keeps the full low-side phase envelope:

$$t_{\min} \;\le\; \lfloor t_{rx} - (\Phi_0^{\min}{+}\Phi_{\text{skew}}^{\min}{+}\Phi_{rj}^{\min})\rfloor$$

so

$$(C) \;\le\; \Phi_0^{\max} + \Phi_{\text{skew}}^{\max} + \Phi_{rj}^{\max}$$

(magnitudes; same envelopes as (B)).

### Combined bound — `capacity_min`

Summing (A)+(B)+(C):

$$\boxed{\;\text{capacity}_{\min} \;=\; \text{blk\_size\_osr} \;+\; (\text{subblk\_size}-1)\cdot\text{osr\_rx}_{\max} \;+\; 2\,(\Phi_0^{\max} + \Phi_{\text{skew}}^{\max} + \Phi_{rj}^{\max}) \;+\; 1\;}$$

For the default `Param` values ($\text{osr}{=}24$, $\text{subblk\_size}{=}32$,
$\text{pi\_ui\_cover}{=}4$, $|\text{ppm}|_{\max}{=}500$, 3 ps skew, 0.3 ps RJ,
$56\,\text{Gb/s}$):

| Term | Value (TX-grid samples) |
|---|---|
| $\text{blk\_size\_osr}$ | $24576$ |
| $(\text{subblk\_size}-1)\cdot\text{osr\_rx}_{\max}$ | $\approx 744$ |
| $2\,\Phi_0^{\max}$ | $\approx 192$ |
| $2\,\Phi_{\text{skew}}^{\max}$ | $\approx 8$ |
| $2\,\Phi_{rj}^{\max}$ | $\approx 6.5$ |
| **`capacity_min`** | $\approx 25\,528$ |

---

## 10. Recommended `capacity`

Three options, in increasing safety margin:

**(i) Tight derived bound** — use the boxed `capacity_min` from §9.
Adaptive to all `Param` values; recompute on construction.

**(ii) Simplified safe default for the current parameter regime**:

$$\text{capacity} \;=\; 2\cdot\text{blk\_size\_osr}$$

Two TX blocks of headroom.  For the default and slider-range parameters it is
about $1.9\times$ the tight bound and is sufficient when the CDR is converged.
For unusually large `subblk_size`, `pi_ui_cover`, skew, or RJ envelopes, use
the §9 bound directly instead of assuming `2*blk_size_osr` is universal.

**(iii) Recommended in practice** — robust against slider extremes:

$$\boxed{\;\text{capacity} \;=\; 4\cdot\text{blk\_size\_osr}\;}$$

The greedy-scheduler proof assumes the CDR is converged.  At slider
extremes ($|\text{ppm}|\sim 500$ with the default CDR gains), the CDR may
fail to track and $\Phi_0$ may swing widely, transiently throttling the
scheduler and growing the residue.  $4\cdot\text{blk\_size\_osr}$ absorbs
this margin and is the value `RxWindow` constructs by default.  Memory
cost: $\approx 768\,\text{kB}$ of Float64 — negligible.

### What the implemented default absorbs

$4\cdot\text{blk\_size\_osr} - \text{capacity}_{\min} \approx 73\,000$ samples
of headroom — covers any plausible future increase in `subblk_size`,
$\text{pi\_ui\_cover}$, skew/RJ, transient deviations from the greedy-scheduler
invariant (e.g. an outer loop that batches two writes before draining), and
CDR-mistracking-induced $\Phi_0$ swings at slider extremes.

---

## 11. Worst-case drift bound (steady state, nominal)

Without $\Phi_0$/skew/RJ, the §8 analysis gives:

$$\bigl|\,t_{\max} - t_{rx}\,\bigr|_{\text{nominal}} \;<\; \text{subblk\_size}\cdot\text{osr\_rx}$$

Independent of $n$ and $\text{nsym\_total}$ under greedy scheduling.  With
perturbations, replace the right-hand side with the §9-(B) envelope.

---

## 12. Ring-buffer indexing

$$\text{buf}[\,t \bmod C + 1\,], \qquad C = \text{capacity}$$

(1-based indexing, Julia convention.)

---

## 13. Candidate/commit state machine (per accepted sub-block)

Prepare a candidate sample-time vector using §3:

$$t_{rx,\text{next}} = t_{rx,\text{cand}} + \text{subblk\_size}\cdot\text{osr\_rx}$$

If `rxw_covers(rxw, Φi)` is false, do not modify committed state:

$$t_{rx}\;\text{unchanged},\qquad \text{pi\_code\_prev}\;\text{unchanged},\qquad \Phi_i\;\text{kept pending}$$

If `rxw_covers` returns true and the sub-block sampling has actually
executed, commit:

$$t_{rx} \leftarrow t_{rx,\text{next}}, \qquad \text{pi\_code\_prev} \leftarrow \text{pi\_code}_{\text{used}}$$

---

## 14. Reclaim before write (`rxw_extend!`)

The scheduler calls `rxw_extend!` with a phase margin from
`rxw_phase_margin(clkgen)` so data below $t_{rx}$ is retained for negative
PI/skew/RJ excursions:

$$t_{\min} \;\leftarrow\; \max\bigl(t_{\min},\;\lfloor t_{rx} - (\Phi_0^{\max}+\Phi_{\text{skew}}^{\max}+\Phi_{rj}^{\max})\rfloor\bigr)$$

This is the §9-(C) reserve.  The low-level `rxw_extend!` API still accepts a
zero margin for tests or custom schedulers, but the production `sim_blk` and
Widget paths pass the phase-envelope margin.

---

## 15. Overrun guard (hard error)

$$\text{abort if:}\quad (t_{\max} + n) - t_{\min} > \text{capacity}$$

where $n$ is the length of the waveform being written.

Rationale: silent overwrite would corrupt the continuous-time analog
signal RX still needs to interpolate, injecting a step discontinuity that
invalidates downstream CDR/BER measurements.  Underrun, by contrast, is
benign — the scheduler loop simply breaks and resumes on the next block.

---

## 16. CDR loop bandwidth — sizing `ki` for a given ppm

Once `pi_wrap_ui` is dropped from $\Phi_0$ (see §3), `pi_code` is the
only fine-tracking mechanism the CDR has.  The integrator gain `ki`
therefore directly bounds the maximum trackable ppm.

### 16.1 Input drift seen by the CDR

Per RX sample the ppm-induced phase drift (relative to TX symbols) is
$\text{ppm}\cdot 10^{-6}$ UI.  In `pi_code` units:

$$\boxed{\;\Delta_{\text{drift}} \;=\; \text{subblk\_size}\cdot\text{pi\_codes\_per\_ui}\cdot\text{ppm}\cdot 10^{-6} \quad [\text{pi\_codes/sub-block}]\;}$$

For the default config ($\text{subblk\_size}=32$, $\text{pi\_codes\_per\_ui}=64$):

$$\Delta_{\text{drift}} \;=\; 2.048\times 10^{-3}\cdot\text{ppm}$$

For $\text{ppm}=100$: $\Delta_{\text{drift}} = 0.205$ pi_codes/sub-block.

### 16.2 Steady-state `ki_accum`

Per CDR call (per vote), $\text{pd\_accum}$ increments by
$\text{kp}\cdot\text{vote} + \text{ki\_accum}$.  In equilibrium
$\text{vote\_avg}\approx 0$ (samples mid-eye on average), so the
integrator must supply the drift entirely:

$$\boxed{\;\text{ki\_accum}_{\infty} \;=\; \frac{\Delta_{\text{drift}}}{N_v}\;}$$

where $N_v$ is the average **biased** votes per sub-block.  For this
simulator:
- $\text{eslc.N\_per\_phi} = [1,0,0,0]$ ⇒ 8 candidate edge positions per sub-block
- $\text{filt\_patterns} = \{[0,1,1],\,[1,1,0]\}$ ⇒ ~25 % of patterns fire on PRBS data
- $N_v \approx 8\times 0.25 = 2$ votes/sub-block

So $\text{ki\_accum}_\infty(\text{ppm}=100) = 0.205/2 \approx 0.103$ —
matches the trace.

### 16.3 Transient overshoot (the binding constraint)

Reaching $\text{ki\_accum}_\infty$ takes time.  During the transient,
the phase error $e$ grows.  If $|e|$ exceeds half a UI ($\sim 32$
pi_codes), the BB-PD votes on the wrong symbol boundary and the loop
locks to a bogus point — the empirical "BER = 0.5" failure mode.

Starting from $\text{ki\_accum}(0)=0$ with saturated votes ($|\text{vote}|=1$)
during pull-in:

$$\text{ki\_accum}(n) \;\approx\; n\cdot\text{ki}\cdot N_v$$

Tracking rate at sub-block $n$ is $\text{ki\_accum}(n)\cdot N_v = n\cdot\text{ki}\cdot N_v^2$, so

$$\frac{de}{dn} \;=\; \Delta_{\text{drift}} - n\cdot\text{ki}\cdot N_v^2$$

$$e(n) \;=\; n\,\Delta_{\text{drift}} \;-\; \tfrac{1}{2}\,n^2\,\text{ki}\,N_v^2$$

Maximum at $de/dn = 0 \Rightarrow n^{\star} = \Delta_{\text{drift}}/(\text{ki}\,N_v^2)$:

$$\boxed{\;e_{\text{peak}} \;=\; \frac{\Delta_{\text{drift}}^2}{2\,\text{ki}\,N_v^2}\;}$$

### 16.4 Lock-loss threshold ⇒ `ki_min`

Require $e_{\text{peak}} < e_{\text{lock}}$ where $e_{\text{lock}} \approx 32$ pi_codes ($=$ half UI):

$$\frac{\Delta_{\text{drift}}^2}{2\,\text{ki}\,N_v^2} \;<\; e_{\text{lock}}$$

$$\boxed{\;\text{ki}_{\min} \;=\; \frac{\Delta_{\text{drift}}^2}{2\,e_{\text{lock}}\,N_v^2} \;=\; \frac{(2.048\times 10^{-3}\cdot\text{ppm})^2}{64\,N_v^2}\;}$$

For $N_v = 2$:

$$\text{ki}_{\min}(\text{ppm}) \;\approx\; 1.64\times 10^{-8}\cdot\text{ppm}^2$$

**`ki` scales quadratically with maximum trackable ppm** — same scaling
as the hold-in range of a classical 2nd-order PLL.

### 16.5 Formula-vs-empirical table

Full $\text{nblk}=977$ sim, default $\text{kp} = 1/2^6$:

| $\text{ki}$ | decimal | $\text{ppm}_{\max}=\sqrt{\text{ki}/1.64{\times}10^{-8}}$ | Empirical break-ppm | Match |
|---|---|---|---|---|
| $1/2^{14}$ | $6.1\times 10^{-5}$ | $\approx 61$ | $\sim 30$–$50$ | $\checkmark$ (within $2\times$) |
| $1/2^{12}$ | $2.44\times 10^{-4}$ | $\approx 122$ | $\sim 100$–$150$ | $\checkmark$ |
| $1/2^{10}$ | $9.77\times 10^{-4}$ | $\approx 244$ | $\sim 300$ | $\checkmark$ |
| $1/2^{8}$  | $3.9\times 10^{-3}$  | $\approx 488$ | $\sim 500$ | $\checkmark$ |

Factor-of-$\sim 2$ spread is explained by: linearization of the BB-PD
(real BB-PD has soft transitions from RJ), the assumption of saturated
votes throughout pull-in, and the binary lock / lose-lock model.

### 16.6 Practical sizing rule

With a safety factor $S \in [2,\,4]$ to absorb the modeling slack:

$$\boxed{\;\text{ki} \;\approx\; S \cdot 1.64\times 10^{-8}\cdot\text{ppm}^2 \;\approx\; 5\times 10^{-8}\cdot\text{ppm}^2\;}$$

The fork's `init_trx()` uses $\text{ki} = 1/2^{12} \approx 2.4\times 10^{-4}$,
which by this rule covers $|\text{ppm}|$ up to roughly
$\sqrt{2.4\times 10^{-4}/(5\times 10^{-8})}\approx 70$ — and empirically
covers the full default $\text{ppm}=100$ run with 0 errors over 898,913
checked bits in the current code.

### 16.7 Role of `kp`

`kp` does **not** affect tracking, only **damping**.  In standard PLL
form the damping factor is

$$\zeta \;\approx\; \tfrac{1}{2}\cdot\frac{\text{kp}}{\sqrt{\text{ki}/K_{pd}}}$$

where $K_{pd}$ is the linearized BB-PD gain (~ $1/\sigma_{rj}$).
Increasing `kp` shrinks overshoot but increases sensitivity to jitter;
typical SerDes CDRs aim for $\zeta \in [0.7,\,1.0]$.  The fork keeps
$\text{kp}=1/2^6$ unchanged.

---

## Symbol reference

| Symbol | Meaning |
|---|---|
| $\text{osr}$ | TX-grid samples per nominal UI |
| $\text{osr\_rx}$ | TX-grid samples per RX-clock UI |
| $\text{ppm}$ | RX/TX frequency offset, parts per million |
| $\text{blk\_size}$ | symbols per TX block |
| $\text{blk\_size\_osr}$ | $\text{blk\_size}\cdot\text{osr}$ — TX-grid samples per block |
| $\text{subblk\_size}$ | symbols per RX sub-block |
| $\text{nsubblk}$ | $\text{blk\_size}/\text{subblk\_size}$ |
| $t_{ui}$ | nominal UI period (seconds) |
| $t_{rx}$ | RX read cursor (TX-grid sample index, Float64) |
| $t_{rx,\text{cand}}$ | pending RX cursor used to build the next candidate sub-block |
| $t_{rx,\text{next}}$ | value committed to $t_{rx}$ after an accepted sub-block |
| $t_{\min}$, $t_{\max}$ | RxWindow time frontiers (TX-grid sample index, Int) |
| $C$ | RxWindow capacity (samples) |
| $\Phi_0$ | CDR / PI phase offset (TX-grid samples) |
| $\Phi_{\text{skew}}$ | per-phase clock skew (TX-grid samples) |
| $\Phi_{rj}$ | per-sample random jitter (TX-grid samples) |
| $\text{pi\_ui\_cover}$ | PI coverage in RX-clock UI |
| $n_{\text{phases}}$ | number of clock phases (e.g. 4) |
| $\text{pi\_codes\_per\_ui}$ | PI resolution: pi_codes per RX UI (e.g. 64) |
| $\text{kp}$, $\text{ki}$ | CDR proportional / integral gains |
| $\text{ki\_accum}$, $\text{pd\_accum}$ | CDR integrator state, total phase accumulator |
| $N_v$ | average biased votes per sub-block (≈ 2 for this config) |
| $\Delta_{\text{drift}}$ | ppm-induced phase drift, pi_codes/sub-block |
| $e_{\text{lock}}$ | lock-loss threshold (~ 32 pi_codes = half UI) |
