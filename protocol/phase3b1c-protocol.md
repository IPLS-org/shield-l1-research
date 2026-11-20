# PROJECT SHIELD L1
## Phase 3B.1C - Pre-registered numerical convergence and reproducibility protocol

**Status:** AUTHORISED DESIGN, NOT EXECUTED. Version 0.2, 20 November 2026 (UTC).  
**Scope:** Pinned PHARE 1D periodic cold-ion / massless-electron dispersion benchmark only.  
**Provenance:** PHAREHUB/PHARE at `7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef` and IPLS-org/shield-l1-research at `8963f401f86f229f74be36389800516c7fdbb007`; successful 3B.1B GitHub run `44618068083`.  
**Authority:** Does not grant permission to claim Project Shield physical validity or launch obstacle studies.

### 1. Observations verified before design

The 3B.1B ZIP SHA-256 is `432cc90cb376cfa69e55b09e2e62780ae4e3dd7c249da89cfb8afba48a7eaa7d`. It contains 24 files, including six EM_B.h5, three plots, the CTest log and a 48-row frequency CSV. All nine recorded diagnostic hashes validated. The archived job reports 3/3 CTests succeeded, in 36.09 seconds total test wall time; the overall Actions workflow was around 23 minutes mainly because of building/installing dependencies. The runner manifest says `nproc=1`, 15 GiB RAM, and 89 GiB free disk; do not assume four simultaneously available CPU cores from hosted-runner marketing specifications.

An independent FFT-based re-extraction of all 48 wave frequencies from the raw HDF5 fields (not the log's printed rows) agrees with the reported rounded frequencies to < 5e-7 in units of Omega_ci. Max absolute relative error against *exact-wavenumber* cold hybrid theory is 5.6624%; against the PHARE **empirical** adjusted reference it is 0.3471%. Amplitude end/start ratios range from 0.7024 to 1.0980. These are **baseline 3B.1B values**, not new convergence outcomes. Largest RMS residual about a linear phase fit is 0.191 radians (high_k, interpolation order 2, AIC, mode 4), which warrants trend diagnostics beyond simple least-squares mean frequency.

The original `low_k` and `high_k` cases do **not** constitute a grid convergence pair: their physical lengths, Fourier modes, cells, dt, ppc, and final times differ. The upstream script fixes wave-mode phases using `np.random.default_rng(0)`, but passes `seed=None` to the particle initialiser. Thus the 3B.1B archive does not document a controlled particle seed.

### 2. Physical problem held constant

Use **high-k** dispersion only for the first convergence screen, with fixed periodic physical length `L = 6.4 d_i`, the same 4 integer Fourier modes `m = 1,2,3,4` on both AIC and whistler branches (`k_m=2*pi*m/L`), `B0=1`, ion number density `n=1`, cold ions `v_th=1e-3 v_A`, massless isothermal electrons at `T_e=0`, transverse modal magnetic amplitude `delta B/B0=0.01`, zero resistivity and hyper-resistivity, simulation end time `10 Omega_ci^-1`, constant recording interval `0.01 Omega_ci^-1` (1001 snapshots including t=0). Wave phases stay fixed. Physical parameters, numerical solver revision, boundary conditions, and analysis window remain unchanged throughout.

**Important correction, v0.2 pre-execution:** In the upstream test, the
initial velocity eigenvectors are computed using grid-dependent `k_grid` and
an empirically adjusted Hall parameter. They therefore change as grid
resolution or interpolation order changes. This is appropriate for an
upstream discrete-eigenmode test but weakens the fixed-physics premise of
spatial convergence. Before any new execution, the study driver has been
adjusted to use **continuum cold-hybrid mode frequencies and exact physical
wavenumbers** when constructing initial magnetic and ion bulk-velocity fields.
The physical functions `B(x)` and `u(x)` and fixed wave phases are identical
for all variants; the particles are still discretely sampled with varying
seeds, ppc, interpolation and grid geometry. This modification touches only
local loaded Python helpers, **not PHARE upstream source on disk or the
compiled solver**. The empirical formula remains available for an explicitly
labelled post-hoc comparison, never for initialisation or acceptance.
A finite-grid run may excite small non-eigenmode components as a result;
assess phase drift and modal amplitudes rather than filtering them out.
A mock-driver test does not prove the compiled solver accepts these settings.


### 3. Smallest robust screening matrix: 18 simulation configurations

Anchor `A`: 64 cells, dx=0.1 d_i, dt=0.0005 Omega_ci^-1, ppc=400, interpolation order 2, particle seed 314159. Each variant changes **one** numerical factor relative to its anchor, except the explicit particle-noise replication groups.

| Family | Configurations | Added runs | Question |
|---|---|---:|---|
| Spatial grid | `N=64,96,128` at seeds `314159,271828,161803` (order=2, ppc=400, dt=0.0005) | 9 incl. anchor | Does measured frequency approach fixed continuum dispersion as dx decreases, accounting for seed variability? |
| Time step | At N=64, seed=314159, dt `0.001,0.0005,0.00025`, all else anchor | 2 | Do time-step effects contract while fixed-grid effects remain? |
| Particle count | At N=64, seed=314159, ppc `100,200,400`, all else anchor | 2 | What changes in sampling noise and fit drift? |
| Interpolation order | At N=64, seed=314159, interpolation orders `1,2,3`, all else anchor | 2 | How strong is order-dependent numerical bias? |
| Additional low-ppc seeds | At N=64, ppc=100, order=2, seeds `271828,161803` | 2 | Compare variability at ppc=100 versus ppc=400 with three seeds each |
| Exact repeat | Repeat anchor with the same fixed seed and settings in a fresh diagnostic directory | 1 | Does identical input reproduce identical field time series? |
| **Total** | | **18** | **144 raw-frequency fits (8 per simulation)** |

At N=64, 96 and 128 the shortest loaded wavelength has respectively **16, 24, 32 cells per wavelength**. The physical domain and each physical `k_m` remain unchanged. All time-step runs retain the fixed **physical** diagnostic interval of 0.01 via `dump_every=10,20,40` respectively. No adaptive step changes or empirical parameter tuning are permitted once results are seen.

Three independent particle seeds per grid level are necessary: a one-seed grid curve could mistake particle noise for systematic discretisation error. Three seeds remain a screening sample, not a precision confidence study. Retain all eight modes: do not silently select favourable branches or modes after measuring results.

### 4. Required execution provenance and resources

Use the existing GitHub Actions `ubuntu-24.04` image and install toolchain/MPI/parallel HDF5 as in 3B.1B. Compile the pinned solver **once per workflow invocation** and execute all 18 configurations sequentially; do not matrix-fan-out 18 separate jobs, each rebuilding PHARE. Use `--parallel 2` for the build and one MPI rank / one thread for controlled tests. Log the real `nproc`, `getconf _NPROCESSORS_ONLN`, `lscpu`, cgroup CPU quota, peak RSS, free memory and disk before execution. Abort before compilation if MemAvailable is below 4 GiB or filesystem free space is below 10 GiB; these are conservative operational thresholds, not measured capacity guarantees. With approximately 15 GiB RAM and 89 GiB free disk observed, 18 one-dimensional cases are expected to fit when run sequentially; this remains a resource **estimate**, not a benchmarked guarantee. Target a workflow timeout of 75-90 minutes for the initial uncached build and tests.

A cross-run binary cache is optional after the first successful 3B.1C workflow; keys must include the exact PHARE commit, Ubuntu runner image, compiler/ABI, CMake configuration and Python versions, and restored artifacts must be validated before use. A **single build inside one job** is the reliable initial approach and already avoids repeated compilation for each case.

Preserve: original research commit, full upstream source SHA, archived benchmark ZIP SHA, exact script hashes, toolchain and Python package manifest, seed, numerical parameters and physical invariants per case, full raw EM_B.h5, unrounded phase-fitted frequencies, phase-residual diagnostics, endpoints and time-series mode amplitudes, all plots, consolidated CSV/JSON, wall times, and checksums. Archive a complete ZIP in the existing Project Shield Audit Reports folder after successful or incomplete execution, clearly identifying failures and missing files.

### 5. Pre-declared analysis metrics

For each `order, N, dt, ppc, seed, branch, m` compute by *independent* FFT over field spatial samples, then linear fit to `unwrap(arg(B_q(t)))`: fitted signed angular frequency (report absolute by branch), residual RMS in radians, initial/end/minimum mode amplitude, and first/second-half fitted frequency drift. Save the complete time-series to allow alternative fitting.

References calculated with exact `k=2*pi*m/L`:

- Continuum hybrid cold-plasma: `omega_aic = k*(sqrt(1+k^2/4)-k/2)` and `omega_whistler = k*(sqrt(1+k^2/4)+k/2)`.
- Staggered-grid derivative: `k_grid=(2/dx)*sin(k*dx/2)` is derived from the ideal centred difference.
- PHARE empirical scheme reference: `di_eff=1-(p+1)*dx^2/12`, inserted into the branch formula with `k_grid`. **Display as an empirical comparison only**. Do not calibrate to 3B.1C observations and do not use its agreement alone as proof of convergence.
- Hall splitting: `omega_whistler - omega_aic`, with continuum target `k^2`. For each branch/mode report signed and absolute relative cold-theory error and residual against empirical scheme separately.

For the 3-seed grids, calculate mean and standard deviation by `(N,branch,m)`, compare errors with seed scatter, fit a documented `omega(dx)=omega_0+C dx^2` extrapolation as a hypothesis rather than a guaranteed convergence order, and report uncertainty and apparent order only when identifiable above noise. For time-step tests compare adjacent finite-dt frequencies at **fixed N**. For ppc=100 versus 400 use three-seed variance by mode; a single ppc=200 run is a directional diagnostic only. Interpolation order effects at fixed dx are not a demonstration of order-independent continuum convergence.

### 6. Pre-declared PASS / FAIL / INCONCLUSIVE gates

**G0. Data completeness and safety.** All 18 runs complete with finite arrays, correct manifest invariants, 1001 ordered physical snapshots, eight identified Fourier modes per run and 144 fitted mode records. Failure to execute is **INCONCLUSIVE (infrastructure)** unless a reproducible numerical instability is established, in which case **FAIL (that configuration)**. Preserve all partial evidence.

**G1. Numerical health of high-resolution/high-ppc tests.** For all ppc=400 modes, no non-finite values, loss of continuous phase tracking, or endpoint amplitude ratio outside `(0.3,3.0)`. RMS phase-fit residual above 0.25 radians or a first/second-half frequency discrepancy above 0.5% forces **INCONCLUSIVE for that mode** pending diagnostics, even if its whole-series average meets a frequency tolerance. Tests deliberately using ppc=100 may exhibit greater noise; record rather than hide any degraded outcomes.

**G2. Exact-repeat reproducibility.** On the same runner and software build the duplicate seeded anchor shall agree to within `1e-8 B0` at all corresponding recorded magnetic field samples and `1e-6` relative in fitted angular frequency. **PASS** if both thresholds hold; **FAIL** if thresholds are exceeded without an identified non-deterministic numerical mechanism. HDF5 file byte hashes alone are not sufficient because file metadata can differ. Different seeds must be demonstrated by non-identical ion-particle realisations or documented model RNG behaviour, not inferred merely from run labels.

**G3. Spatial convergence.** Use 3 seeds at each N=64,96,128. For every one of the eight modes report the extrapolated zero-spacing `omega_0` and a 95% bootstrap uncertainty interval based on resampling the three seeds within each resolution (at least 2000 reproducibly seeded bootstrap replicates). **PASS** only if every interval is wholly inside `+/-1%` of the *uncorrected continuum* frequency **and** finite-spacing errors show no statistically resolved reversal of convergence. **FAIL** if an uncertainty interval lies wholly outside `+/-1%` or the grid trend is reproducibly divergent; otherwise **INCONCLUSIVE**. State when the three-level fit is noise dominated, ill conditioned or not a credible second-order model. The empirical scheme correction is not an acceptance target for G3.

**G4. Time-step stability/convergence.** At fixed N=64, each of eight modes shall differ by no more than 0.25% in fitted frequency between `dt=0.0005` and `0.00025`. A resolvable reduction in error as dt halves is supporting evidence; when differences are below measurement noise, report plateau rather than claim a convergence order. A mode exceeding 0.25% in the initial screen is **INCONCLUSIVE pending an independently pre-registered repeat**; classify **FAIL** only if the excess is confirmed by that repeat. **INCONCLUSIVE** also applies if phase-fit instability or insufficient diagnostic accuracy prevents a clean comparison; otherwise **PASS**.

**G5. Particle noise and interpolation.** Verify that independently seeded particle initialisations produce distinguishable raw fields at the anchor settings; if they do not, classify noise inference INCONCLUSIVE. Report standard deviations across seeds separately at ppc=100 and 400, and the aggregate ratio. A decrease in aggregate seed-to-seed spread at ppc=400 is the expected direction; **PASS** as a *screening noise check* if aggregate spread at ppc=400 is <=80% of the ppc=100 value and no ppc=400 mode violates G1. Do not infer an exact inverse-square-root sampling law from three seeds. A poorly resolved/noise-floor ratio is **INCONCLUSIVE**, not an automatic physical failure. Show p=1,2,3 frequency differences against both references; no order-independent convergence claim at fixed spatial resolution.

**Overall gate:** numerical solver qualification for the *tested 1D problem* is **PASS** only if G0-G5 pass or an explicitly bounded deficiency is closed by a separately pre-registered extension. A scientific result can be **INCONCLUSIVE** without being a software failure. A PASS does not validate Project Shield and does not automatically authorise magnetic obstacle or plasma-inflation research. An independent-code comparison and broader 2D/3D, boundary-condition and conservation tests remain necessary.

### 7. Limits, follow-up and scientific gate

This is a linear, small-amplitude, quasi-cold, 1D periodic wave benchmark with an electron-fluid closure; it does not probe finite-coil topology, injected plasma equilibrium, current closure, nonlinear inflation, IMF coupling, reconnection, open boundaries, transport to Earth or net geoeffectiveness. The source's empirical particle-shape formula may fit the original benchmark population by construction. It is not a separately validated constitutive law. The 3-seed screen has limited precision; expand to at least 5-10 seeds per condition if bootstrap intervals overlap G3 limits or noise masks ordering. Prespecify any expansion before running it.

**Decision at protocol review:** 3B.1B archived raw-field analysis PASS; **3B.1C numerical convergence is NOT YET EXECUTED and therefore INCONCLUSIVE**. Proceed only to a controlled experimental workflow after static verification of driver and analysis code; do not initiate magnetic obstacle or inflation simulations at this gate.
