# Project Shield L1 | Independent scientific feasibility research

**Current scientific status (21 November 2026):** Phase **3B.2B-05R** has completed a bounded analytical and source-code remediation of magnetic background-field and conservation-accounting questions. **The protective mechanism is unproven. No multidimensional MPI-AMRVAC Shield simulation has been built or run.** The next step is a reviewed implementation design and a separately authorised, prospectively specified small solver-qualification tranche, **not** an obstacle or Earth-protection experiment.

Project Shield L1 investigates whether an artificial magnetic/plasma obstacle near Earth-Sun L1 could mitigate harmful solar-wind interactions with Earth's magnetosphere. This is a research hypothesis, **not** a demonstrated shield, deployable design, operational protection system or validated engineering specification.

## What has actually been established

| Programme stage | Verified scope | What it does **not** establish |
|---|---|---|
| Phase 3B.1A to 3B.1D | Bounded PHARE hybrid solver benchmarks, numerical reproducibility, controlled sensitivity and conservation diagnostics for specified 1D cases | No 3D source-driven obstacle or geoeffectiveness |
| Phase 3B.2B-01 to 03 | Ideal finite circular-current mathematical reference: vacuum field, energy, force and Maxwell-stress consistency | No physical winding, current feed, structure, cryogenic design or plasma inflation |
| Phase 3B.2B-04 | Standalone finite-coil exterior static field mapping: decreasing sampled field errors under grid refinement and discrete curl/divergence consistency | No compiled/run MPI-AMRVAC model or source/plasma boundary qualification |
| Phase 3B.2B-05 | Prospective MPI-AMRVAC v3.2 source/CT/boundary qualification protocol, source audit, failure modes and proposed resource ceiling | No executed MHD benchmark, approved resource spend or empirical solver pass |
| Phase 3B.2B-05R | Analytical compatible-face-flux design, manufactured exact-arithmetic counterexample to point-sampled B0 divergence, and source-level split-field energy/momentum accounting | No installed code patch, solver execution, boundary/AMR validation or protective effect |

**Research integrity:** A successful workflow means a computation ran, not that the physical Shield mechanism works. Distinguish analytical identity, static numerical screen, qualified executable, converged plasma result, and physical validation. The hypotheses about solar-wind-driven magnetic reinforcement (**H-3B2-01**) and expansion/compression/recovery stability (**H-3B2-02**) remain **untested**.

## Current principal blockers

1. **Background B0 and constrained transport:** In the pinned MPI-AMRVAC source, `src/physics/mod_b0.t` populates stationary cell/face background fields by sampling `usr_set_B0`, while `src/mod_constrained_transport.t` creates an evolved staggered field from edge circulation. Controlling `div_h b` does not automatically control `div_h(B0+b)`. A compatible shared-edge, face-flux B0 representation is the preferred *design*, not yet a tested solver feature.
2. **Conservation:** For stationary `B0`, the solver's split energy contains the perturbation magnetic term; full physical energy also contains `|B0|²/2 + B0·b` (in code-normalised units). Boundary energy flux, momentum, discrete source-current/force residuals and numerical corrections must be independently closed.
3. **Source and boundaries:** Real coil current closure, energisation and protection are not engineered. `usr_internal_bc` is a state-overwrite hook, not a validated conservative conductor interface. Large reflecting exclusions may produce false-positive deflection.
4. **Readiness:** The reproducible build, total-field/AMR/BC diagnostics, independently justified quantitative acceptance limits and explicit operator approval remain outstanding.

**No execution authorisation:** Do not run a coupled solar-wind/artificial-obstacle calculation, build a new solver configuration, dispatch GitHub Actions or spend computational resources merely because this README, a protocol or a source-level note exists. The prospective T00-T04 resource ceiling in the protocol is **not** authorisation. All solver tests require a separate, frozen and explicitly approved plan.

## Scientific evidence and programme continuity

- [Programme Continuity folder, latest versioned `memory-v*.md`](https://drive.google.com/drive/folders/1gRs5z67uW4xNc4lyfXGL0Na2A17vqLXx). The detailed continuity record, not this landing page, is authoritative for gate decisions and source provenance.
- [Authoritative Audit Reports folder](https://drive.google.com/drive/folders/1hAXfckULIkFKdpWjhXpZuf_WK824IJN1). Two similarly named folders exist; this exact folder is authoritative.
- [Phase 3B.2B-04 numerical evidence ZIP](https://drive.google.com/file/d/1k1J2WKlFw4CawmjfxeB7Pj0BrfBgEhsk/view) and [scientific assessment](https://drive.google.com/file/d/1nrbJptABh72cPhQpUuZnX0LuBqWorw2i/view).
- [Phase 3B.2B-05 prospective qualification protocol and audit](https://drive.google.com/file/d/1M8dVNez7y_TrAxE0aZHt84uHNH34nztJ/view).
- [Phase 3B.2B-05R source/energy remediation and exact-arithmetic evidence](https://drive.google.com/file/d/1P_E5ocEffJD-29ScRzmqKOmmBLtXILn7/view).
- [Pinned candidate MPI-AMRVAC v3.2 source commit](https://github.com/amrvac/amrvac/tree/c45d84650ed3379eb22dcfa6fb927a45734e5d8c) (upstream third-party solver; not included or qualified here).

## This repository's executable research history

The repository currently contains the earlier PHARE hybrid-solver workflows and test scripts, including:

- `.github/workflows/phare-perturbed-conservation-3b1d.yml`
- `scripts/p3b1d_perturbed_run.py`, `scripts/p3b1d_perturbed_analyse.py`
- `protocol/phase3b1d-one-perturbed-v0.4.json`, plus controlled time-step, grid and seed sensitivity workflows.
- `tests/test_p3b1d_perturbed_offline.py`, which checks synthetic HDF5-shaped data only, **not** a PHARE solver result.

**Historical provenance note:** The former README was a 20 November 2026 *package-specific* Phase 3B.1D pre-registration instruction sheet. That content described a pending one-case workflow and is superseded **as a programme landing page**. The original 3B.1D source files, exact parameter/acceptance contracts and workflows remain in Git history and in the named repository paths. Do not dispatch or repeat them without a fresh, bounded approval and evidence review.

## Participation and standards

Independent contributions are welcome, particularly source-level verification, plasma physics, conservative numerical flux methods, falsification controls, and reproducibility audits. Contributions should preserve original raw evidence and exact solver commit/parameters, pre-register scientific tests, maintain negative controls, distinguish observed results from hypotheses, record explicit human approval, and cite the original audit evidence rather than asserting Shield feasibility from plots.

*Maintained as a research status index. Prepared 21 November 2026 UTC, pending GitHub publication.*
