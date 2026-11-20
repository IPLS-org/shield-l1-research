# Project Shield L1 | One perturbed conservation screen, pre-registered v0.4

**Status: prepared for source review, not yet executed. This ZIP contains no new PHARE numerical results.** The completed preceding unperturbed diagnostic qualification is GitHub run 44774581522.

## What this package provides

- Exactly one proposed 1D source-free, zero-electron-temperature PHARE hybrid case, same pinned source as the successful diagnostic smoke. Initial fields `B=(1, 0.01cos(2πx/6.4), 0.01sin(2πx/6.4))` with zero initial ion bulk flow. Grid 64, particles per cell 400, time step 0.001, final time 1.0, 1,000 steps and 11 common HDF5 observation times.
- An executable run script and a **separately implemented**, read-only HDF5 analysis script. The latter calculates weighted N, P and K, Yee-centred B energy, checks source/diagnostic contracts, and compares absolute excursion against *predeclared* scales. It does **not** import PHARE or trust simulation-reported metrics.
- A pinned GitHub Actions workflow, manually dispatched only, with an explicit `confirm_one_case=YES` input, double-precision diagnostics, resource limits and failure-evidence preservation.
- A complete JSON parameter/acceptance contract and readable review protocol.
- Offline tests constructing **synthetic PHARE-shaped HDF5** (not PHARE output). All scientific claims require a future authorised real run and independent readback.

## Place exactly these FOUR files into the GitHub repository at the root-relative paths shown

1. `.github/workflows/phare-perturbed-conservation-3b1d.yml`
2. `scripts/p3b1d_perturbed_run.py`
3. `scripts/p3b1d_perturbed_analyse.py`
4. `protocol/phase3b1d-one-perturbed-v0.4.json`

Do not replace or remove any successful 3B.1C or Phase 3B.1D smoke workflow. The Markdown protocol, tests and README are included for review and reproducibility but are not required to run the workflow. GitHub organisation action policy requires **full commit SHAs**: this workflow pins the two third-party actions accordingly.

## Approval and execution, deliberately separate

1. Review the protocol and its numerical acceptance limits. Commit the four files to the repository and note the final commit SHA. This is pre-registration. If any change to scientific limits is needed, make it **before** an experimental run and document the revised commit.
2. Do **not** dispatch until explicit human approval has been recorded for this **one** case. The workflow starts with input `confirm_one_case=NO` and aborts unless manually changed to `YES`.
3. In GitHub Actions choose **Project Shield L1 - Phase 3B.1D One Perturbed Conservation Screen**, branch `main` and `confirm_one_case=YES` once. Do not start an additional run in response to unknown performance or an early failure without reviewing the evidence.
4. On completion, retrieve the workflow artifact, inspect the raw HDF5 and independent-assessment/assessment.json, independently recompute metrics, archive full raw evidence in Project Shield Audit Reports before 14-day expiry, and update programme continuity in a new version.

## Strict scientific scope

The algorithm's result `BOUNDED_SINGLE_CASE_SCREEN_PASS` (if achieved) means **only** that the observed single case met the prospectively set numerical-quality screen. It is not a grid/time-step convergence study, formal proof of conservation, full PHARE qualification, physical confirmation of plasma inflation, or Project Shield L1 shielding validation.

Resource controls: one workflow job with 45-minute global timeout, one 1,000-step solver invocation capped at 300 seconds, 30-minute build subprocess cap, available memory minimum 4 GiB, free disk minimum 10 GiB, max 16 MiB per simulation output file and 64 MiB aggregate HDF5, 14 days evidence retention. GitHub runner limits and billing depend on account and platform state.

## Synthetic-test command (local only)

```bash
python -m pytest tests/test_p3b1d_perturbed_offline.py -q
```

These tests generate data in temporary directories. They never import the PHARE solver, and their success is no evidence of actual conservation behaviour.

Prepared 20 November 2026 UTC. Pinned PHARE SHA `7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef`.
