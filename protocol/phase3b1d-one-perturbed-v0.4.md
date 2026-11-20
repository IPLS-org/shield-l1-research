# Project Shield L1 | Phase 3B.1D: One deliberately perturbed ideal-hybrid conservation screen

**Protocol revision:** v0.4 REVIEW DRAFT, 20 November 2026 UTC. **No PHARE execution requested or performed by this package.** The operator's current decision is to move towards one experiment, not to initiate a broad campaign. Freeze this document, its JSON matrix and hash in an approved Git commit before manually dispatching anything.

## 1. Decision and scientific question

We have independently verified the successful Phase 3B.1D unperturbed 50-step PHARE smoke run, GitHub Actions [44774581522](https://github.com/IPLS-org/shield-l1-research/actions/runs/44774581522), and archived its raw evidence and a separate HDF5 reanalysis in the Project Shield Audit Reports folder. That run established the availability, precision, schema, time alignment and physical-patch ownership of the required diagnostic data. It did **not** challenge an imposed magnetic perturbation or establish conservation convergence.

The next smallest falsifiable question is: **In the pinned PHARE source-free periodic cold-electron hybrid model, can we introduce a non-zero spatially varying transverse magnetic field and observe a bounded drift in global ion number, ion momentum and ion-plus-magnetic energy measured against the perturbation scale?** A single result is a *screen*, not a discretisation error order, formal conservation theorem or shielding result.

## 2. Exactly one physics case and preregistered geometry

- One proton-like population, charge=mass=1 in PHARE normalisation; a single-level periodic 1D box of length `L=6.4 d_i`; 64 cells; `dx=0.1 d_i`; interpolation order 2; 400 particles/cell; fixed seed 314159.
- Cold ions with a finite numerical thermal width `v_th=0.001 v_A` and **zero initial bulk velocity**. Quasineutral isothermal massless electrons, `T_e=0`, `η=0`, hyper-resistivity=0, no source, AMR, obstacle, open boundaries or inflow.
- Background `B_x=1`; **circular transverse perturbation** at Fourier mode `m=1`: `B_y(x,0)=0.01*cos(2πx/L)`, `B_z(x,0)=0.01*sin(2πx/L)`. This is *not* claimed to be an exact AIC eigenmode. It is deliberately out of force balance and expected to produce a transient and dynamic energy exchange.
- `dt=0.001 Ω_ci^-1`, `t_final=1 Ω_ci^-1`, **1,000 solver steps**, 11 joint field/ion-particle snapshots at `t=0,0.1,...,1.0`. No interpolation, inferred intermediate time or change of numeric controls after results.
- The **same PHARE upstream SHA** as the qualified smoke case: `7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef`; diagnostics built using `-DPHARE_DIAG_DOUBLES=1`.

The 50-step unperturbed smoke run at 100 particles/cell took 0.69 seconds in the solver and about 21 minutes overall including build. The proposed 1,000-step, 400-particle/cell test is computationally larger. Its runtime cannot be inferred by linear extrapolation with confidence. Cap its subprocess and preserve failures.

## 3. Continuum balance and numerically measurable invariants

For ideal periodic massless-electron hybrid kinetics, `J=(∇×B)/μ0`, `u_e=(J_i−J)/ρ_i`, `E=−u_e×B` when `T_e=η=ν=0`. The kinetic ion number integrates to constant for a source-free periodic domain. The volume-integrated `dP_i/dt=∫(J×B)dV=0` by periodic magnetic stress balance (assuming ∇·B=0). With Faraday and ion work, `d(K_i+U_B)/dt=∫(J_i−J)·E dV=0` in the cold ideal continuum. These laws **do not guarantee exact discrete invariance** of PHARE's particle-grid solver. In a one-dimensional test, only `B_x` constancy is probed; neither multidimensional divergence nor current closure are qualified.

Use the source-confirmed `w_p=n_cell/ppc` weight convention. Global physical 1D integrals are `N=dx*Σw_p`, `P=dx*Σm*w_p*v_p`, `K=(dx/2)*Σm*w_p*|v_p|²`. Magnetic energy is `UB=(dx/2)[Σ_unique_primal Bx²+Σ_dual By²+Σ_dual Bz²]`. All sums exclude ghost cells and duplicate periodic primal nodes. E is recorded but *not* the hidden pusher-averaged E, so an instantaneous discrete work-budget equality is expressly **out of scope**.

## 4. Independent observation contract

Exactly these five HDF5 files at all 11 common times: `EM_B.h5`, `EM_E.h5`, `ions_charge_density.h5`, `ions_charge_flux.h5`, `ions_pop_protons_domain.h5`. Require proper source/runtime JSON provenance, single level 0 and physical patch coverage 0..63, 64-bit fields/velocities/weights, positive weights, finite data, integer particle cell indices, population mass/charge=1, 25,600 owned macro-particles, and matching root geometry/periodicity/centering. Output exceeding 16 MiB in a single HDF5 file or 64 MiB total is outside the registered run's budget. **Any missing/uncertain schema is INCONCLUSIVE**, not a physical conservation failure.

## 5. Frozen candidate decision thresholds (prospective engineering screens)

These limits are acceptance *policy*, not predictions of PHARE error. Revisions must be made **before** the first real perturbed data are observed. No post hoc adjustment to make a run pass. The JSON matrix is the authoritative machine-readable contract.

| Gate | Registered rule | Outcome if violated with valid measurement |
| --- | --- | --- |
| G0 | Full valid source, data, time, ownership, float64, completeness and resource contract | INCONCLUSIVE for missing data, INVALID_CASE for wrong physics |
| G1 | At `t=0`, observed Fourier first-harmonic complex transverse amplitude within ±2% of 0.01; correct helicity | FAIL: prescribed initial condition not realised |
| G2 | Corrected total ion number initially 6.4 within 1e-10; maximum relative drift ≤ `max(1e-12, 8 ε64*25,600)`; deposited density integrated within 1e-8 | FAIL: trusted particle-number inconsistency |
| G3 | Longitudinal Bx within absolute 1e-12 of 1 across all snapshots; periodic endpoint mismatch ≤1e-12 | FAIL: tested 1D Bx consistency |
| G4 | `max ||P(t)−P(0)||` ≤ 5% of `S_P=L*0.01=0.064`, so absolute limit **0.0032** | FAIL: nontrivial momentum-quality screen |
| G5 | `max |(K(t)−K0)+(UBx−UBx0)+(UBy−UBy0)+(UBz−UBz0)|` ≤ 5% of `S_E=L(0.01)²/2=0.00032`, absolute limit **0.000016** | FAIL: perturbation-scale energy-quality screen |
| G6 | Report ion kinetic and transverse magnetic energy changes and field response | Exploratory; no PASS rule |
| G7 | Charge continuity and discrete work using time-centred pusher fields | NOT TESTED; missing internal averaged fields |

For G5, add individual component changes using compensated summation **instead of taking a difference of two background-dominated totals near 3.2**. Report the absolute magnitude and ratio to `S_E`, not just the insignificant fraction of `U≈3.2`.

A simultaneous G0-G5 PASS is **BOUNDED_SINGLE_CASE_SCREEN_PASS** only. Failure of a reliable quantitative screen is **BOUNDED_SINGLE_CASE_SCREEN_FAIL**, and must be traced to its cause before claiming a solver flaw. Neither outcome demonstrates a mesh/time-step convergence order, seed-to-seed reproducibility, plasma inflation, protective deflection or Project Shield feasibility.

## 6. Resource and evidence controls

One manually triggered `workflow_dispatch` job on `ubuntu-24.04`, 45-minute total timeout, configure/build cap 30 minutes per subprocess, solver cap 300 seconds, minimum available memory 4 GiB and free disk 10 GiB, each HDF5 file capped at 16 MiB, aggregate HDF5 <=64 MiB, exactly one 1,000-step case. Using the permitted full-length pinned SHA for `actions/checkout` and `actions/upload-artifact` is mandatory under IPLS-org organisation policy. No unapproved costly cloud resources or repeated dispatch.

Preserve original configuration JSON and its SHA256, PHARE and research commits, runner/toolchain and package versions, HDF5 datasets, all logs, timing/memory, independent CSV/JSON results, absolute per-file hashes, all evidence even if a gate fails, and a 14-day GitHub Actions artifact. Download and permanently archive the raw bundle in Audit Reports **before GitHub retention expires**. Independently recompute results from raw HDF5 in a separate analysis environment, with no use of the workflow report as ground truth. Record PASS, FAIL or INCONCLUSIVE, never a general conservation-qualified PASS from one case.

## 7. No-run boundary and next decision

The files in this ZIP are **execution-ready proposals**. They have no automatic schedule and were not committed to the repository or run by the assistant. The GitHub integration previously returned 403 for repository writes, so the user may need to add the four root-relative files and manually trigger the reviewed workflow. The user's wording 'move towards' is treated as authorisation to *prepare and review*, not an instruction to dispatch the PHARE job silently. Once reviewed, separately confirm the one-case operator approval and run once. The next decision is whether the resulting one-case quality screen supports authorising a second resolution or noise-control case. No extrapolation to an L1 magnetic obstacle.

### Primary evidence

- https://github.com/IPLS-org/shield-l1-research/actions/runs/44774581522
- https://github.com/PHAREHUB/PHARE/tree/7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef
- https://drive.google.com/file/d/1FpYpDyRoBZJ49hIFTKgKWvipNA2dbK46/view
- https://drive.google.com/file/d/1Hxu1YKwthfPM2rTVPS2tmCWi3gY3TGR-/view
- https://docs.google.com/document/d/1EYcoq7_DyzqUFuiY955OxEPhcXvHP-uiAsT-R12L-8s/edit
