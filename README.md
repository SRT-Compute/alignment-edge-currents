# Alignment edge currents in turbulence

Code, analysis plans and machine-written results for the manuscript

> J. L. Jensen, *Kinematic orientation and dynamical redistribution of
> vorticity–strain alignment currents in turbulence* (submitted).

The study measures the pairwise currents between the vorticity alignment weights
c_i = (omega_hat . e_i)^2 that are produced by the restricted-Euler terms of the
strain-eigenframe dynamics,

    J[i<-j] = 2 c_i c_j (lambda_i - lambda_j) (1 - chi_ij),   chi_ij = |omega|^2 / (4 (lambda_i - lambda_j)^2),

in forced isotropic turbulence at Re_lambda = 611 and 1300 and in turbulent
channel flow at Re_tau = 1000 (Johns Hopkins Turbulence Databases), and compares
them with a variance-matched Gaussian reference.

## Contents

```
code/
  instrument.py      finite differences, eigenframe, alignment weights, edge currents
  datasets.py        frozen cutout origins, snapshots, channel blocks, file layout
  fetch_data.py      download from JHTDB (requires a JHTDB token)
  gates.py           quality gates
  confirmatory.py    preregistered sign and organisation test (2000 points per cutout)
  dense.py           per-cutout statistics on all interior points, reference ensembles,
                     resolution variants, filtering, (Q, R) quadrants, enstrophy bins
  stats.py           pooled ratios, block bootstrap, sign-flip test, Holm correction
  windows.py         window-level results and verdicts
  derived.py         ratios to the Gaussian reference, edge asymmetry, net drifts N_i,
                     cycle component K, stretching-efficiency rate
  channel.py         channel-flow analysis (stages 1 and 2)
  re_check.py        check of the current formula against direct restricted-Euler evolution
  figures.py         figures of the manuscript
  run_all.py         runs everything
prereg/              English renderings of the locked analysis plans, with checksums
results/             machine-written result files and figures produced by run_all.py
CHECKSUMS.sha256     SHA-256 of every file in this repository
```

## Reproducing the results

Requirements: Python 3.11 or later with numpy, scipy, h5py, matplotlib and the
JHTDB client `givernylocal` (`pip install -r requirements.txt`).

```
cd code
JHTDB_TOKEN=<your token> python fetch_data.py --data-dir ../data
python run_all.py --data-dir ../data --out-dir ../results
```

The download is about 1.5 GB (60 velocity cutouts of 128^3 and 120 channel
blocks of 2000 points). The analysis takes about one hour on a desktop computer.
All random draws are seeded; `run_all.py` reproduces the files in `results/`.

## Where the numbers of the manuscript come from

| Manuscript | File | Entry |
|---|---|---|
| Table 2 (confirmatory test) | `results/confirmatory.json` | `stage1`, `holdout` |
| Table 3 (normalised currents) | `results/edge_balance.json` | `norm`, `O` |
| Table 4 (ratios, asymmetry) | `results/edge_balance.json` | `ratio_to_gaussian`, `asymmetry`, `asymmetry_gaussian` |
| Table 5 (net drifts, cycle component) | `results/edge_balance.json` | `N`, `N_gaussian`, `K`, `K_gaussian`, `K_ratio` |
| significance, heterogeneity | `results/dense_results.json` | `significance` |
| references, resolution, (Q, R), tails, domain | `results/dense_results.json` | `references`, `resolution`, `topology`, `tails`, `domain` |
| stretching-efficiency rate | `results/stretching_efficiency.json` | `rate`, `rate_gaussian`, `rate_edge_*` |
| channel, stage 1 | `results/dense_results.json` | `channel_stage1` |
| channel, stage 2 (confirmation) | `results/channel_stage2.json` | |
| gates | `results/gates.json` | |
| appendix A check | `results/restricted_euler_check.json` | |
| Figures 1-8, S1 | `results/figures/` | |

## Analysis plans

The analyses were specified in advance in four documents that were hash-locked
before the corresponding data were downloaded. English renderings and the
checksums of the originals are in `prereg/`; see `prereg/NOTE_ON_CHECKSUMS.md`.
Quantities that were derived after the measurements (edge asymmetry, net drifts,
cycle component, stretching-efficiency rate) are labelled as such in the
manuscript.

## Data

All flow data are from the Johns Hopkins Turbulence Databases
(https://turbulence.idies.jhu.edu): isotropic4096, isotropic8192,
isotropic1024coarse and channel. Please cite the databases as recommended there.

## License

Code: MIT. Documents and result files: CC BY 4.0. See `LICENSE`.

Copyright (c) 2026 Jesper Lyng Jensen, SRT Compute ApS.
