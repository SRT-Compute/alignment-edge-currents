# Analysis plan 5 — analyses added in revision

English rendering of the locked document (3 November 2026; checksum in
`NOTE_ON_CHECKSUMS.md`). The plan and the analysis script were hash-locked
together before the analyses were run. The rules below were fixed from the
geometry of the cutouts alone; no outcome data had been seen. The analyses are
post hoc revision analyses and are reported regardless of outcome.

Scope: no new DNS download, no Lagrangian campaign and no analysis near R = 0
are part of this plan.

## Data and code

- The repository code is used unchanged (`instrument`, `dense`, `stats`,
  `derived`, `datasets`).
- Per-cutout cache written by `run_all.py` (46 cutouts that passed the gates).
  Before the analysis, the normalised mean currents of the cache are checked
  against the published values (`dense_results.json`) to within 1e-9; otherwise
  the run stops.
- The cutout files are used only in B2.

## A. Reviewer 4, comment 3: spatial dependence between cutouts

**A1. Geometry**

- For all pairs within the same snapshot:
  - periodic centre distance;
  - edge gap, i.e. the distance between the cutout boundaries (0 when they touch).
- Sets:
  - Re 611: 24 cutouts from one snapshot;
  - Re 1300, snapshot 0: 10 cutouts (stage 1);
  - Re 1300, snapshot 3: 12 cutouts (holdout).
- Pairs across snapshots are not included. Snapshots 0 and 3 are 1.12 time
  units apart, about 1.4 large-eddy turnover times.
- Distances are given in code units and in integral lengths L (JHTDB L1):
  snapshot 0: 1.2438; snapshot 3: 1.3192; for isotropic4096 no value is listed
  in the JHTDB documentation, and the value 1.2983 of the 8192^3 run at
  Re_lambda 613 (snapshot 5, eta within 0.2 %) is used as an approximation.

**A2. Dependence test**

- Variable: primarily the block mean m_b per edge (database units); secondarily
  m_b normalised by the block's S_rms. The variable is standardised to z_b within
  the set.
- Statistics:
  - rho_near = mean of z_i z_j over pairs with edge gap < d* = 1.3;
  - combined statistic = mean of rho_near over the three edges;
  - secondary: Spearman correlation between edge gap and z_i z_j.
- Null distribution: 10 000 permutations of the block positions (seed 20261105).
- Primary test: combined rho_near (raw variable) per set, one-sided (> 0), Holm
  over the three sets.
- Verdict: `DEPENDENCE_DETECTED` if a Holm-adjusted p < 0.05, otherwise
  `NO_DEPENDENCE_DETECTED`.
- With 10–24 blocks the power of the test is limited. A non-significant outcome
  is not evidence of independence, and the manuscript states it that way.

**A3. Robustness under geometric rules**

The rules were fixed from the geometry alone; the numbers of units were computed
before this plan was locked. Numbers of units are given in the order
pooled_611 / holdout_611 / pooled_1300 / holdout_1300.

| Rule | Definition | Units |
|---|---|---|
| `cutout` | cutout as unit, as in the manuscript (reference) | 24 / 12 / 22 / 12 |
| `group13` | single linkage within a snapshot: edge gap < 1.3 gives the same group | 10 / 6 / 17 (8+9) / 9 |
| `thin13` | greedy removal of the cutout with most neighbours < 1.3 (ties: smallest neighbour distance, then highest position) until no pair is closer than 1.3 | 14 / 7 / 17 / 9 |
| `thin26` | as `thin13` with 2.6 ≈ 2L (secondary, descriptive) | 5 / 4 / 11 / 6 |

- Grouping at 2.6 merges all cutouts of a set into one group and is not used.
- For each rule and window a unit is the sum of its block statistics (super-block).
- The estimators are those of the manuscript (`derived._quantities`): normalised
  currents, ratios to the Gaussian reference, edge asymmetry, N_i, K and the K ratio.
- Uncertainty: 95 % bootstrap over units (10 000 replicates, seed 20261104).
- Holdout windows: one-sided sign-flip test with units as units, Holm over six
  tests. The smallest attainable p is reported.

Criteria per pooled window:

- (a) the three edges have the signs of the manuscript, and the 95 % interval excludes zero;
- (b) N_2 > 0, and the interval excludes zero;
- (c) all three ratios to the Gaussian reference have an interval excluding 1;
- (d) the asymmetry J_12/J_23 has an interval above 1.

Verdict: `ROBUST_TO_SPATIAL_GROUPING` if (a)–(d) hold under `group13` and
`thin13` at both Reynolds numbers; otherwise `NOT_ROBUST_TO_SPATIAL_GROUPING`.

Consequence fixed in advance: with 6–7 holdout units the sign-flip test cannot
reach the level of the manuscript (Holm p < 0.01); the smallest p is 1/64–1/128
before Holm. The significance statement is therefore reformulated, and intervals
with units are reported side by side with those of the manuscript.

## B1. Reviewer 4, comment 4: N_2 under trimming

- Definition: n_2(x) = −J_12(x) + J_23(x), the pointwise restricted-Euler rate
  of change of c_2; N_2 = <n_2>.
- Data:
  - DNS: the subsamples of the cache, 50 000 points per block, the same as in Figure 6;
  - Gaussian: per block 50 000 new states with the block's <S:S> and <|omega|^2>
    (`gaussian_reference`, seed 20261103 + position) and the same gap selection (g > 0.01).
- Primary variant T0: remove the largest |n_2|.
  - The threshold is global over the window's pooled subsample, as in Figure 6.
  - Levels: 0 (untrimmed), 0.1 % and 1 %; 5 % and 10 % are descriptive.
  - Rationale: N_2 is the mean of n_2. n_2 is odd under A -> −A, so trimming on
    |n_2| preserves the Gaussian symmetry, and the trimmed Gaussian N_2 is zero
    in expectation.
- Secondary variants:
  - V1: edge-wise trimming as in Figure 6, combined into N_2;
  - V2: remove the largest |J_23| from both edges (the tail that reverses 2<-3).
- Normalisation by the window's S_rms.
- Uncertainty: 95 % block bootstrap with 2 000 replicates, the thresholds being
  recomputed in each replicate (seed 20261104); for pooled_611 and pooled_1300.
  The holdout windows get point estimates only.
- Verdict: `N2_ROBUST_UNDER_TRIMMING` if the DNS 95 % interval for T0 lies above
  zero at 0, 0.1 % and 1 % in both pooled windows; otherwise
  `N2_NOT_ROBUST_UNDER_TRIMMING`, and the main conclusion is weakened.
- Gaussian control: the interval is expected to cover zero; the outcome is reported.

## B2. Reviewer 4, comment 4: characterisation of the 2<-3 tail

- Thresholds: tau_q are the manuscript's trimming thresholds for |J_23| (global
  quantile over the window's subsample), q = 0.1 % and 1 %.
- All points with g > 0.01 are recomputed from the cutout files with the
  repository instrument. The tail is the set of points with |J_23| > tau_q.
- Descriptive quantities, tail against all points:
  - share of the sum of J_23;
  - fraction negative;
  - fraction with chi_23 > 1;
  - fraction with |omega|^2/<|omega|^2> > 4;
  - fraction with g_23 = Delta_23/|S| < 0.05;
  - medians of |omega|^2/<|omega|^2>, g_23, Delta_23/S_rms, c_2 c_3 and chi_23;
  - shares of the (Q, R) quadrants in points and in the sum of J_23.
- Decomposition of the size of the rotation term |omega|^2 c_2 c_3/(2 Delta_23),
  as the difference of mean logarithms (tail − all):
  Delta<ln> = E_enstrophy + E_alignment + E_gap + E_block, with
  E_enstrophy = Delta<ln(|omega|^2/<|omega|^2>_b)>,
  E_alignment = Delta<ln c_2 c_3>,
  E_gap = −Delta<ln(Delta_23/S_rms,b)>,
  E_block = Delta<ln(<|omega|^2>_b/S_rms,b)>.
  The identity is exact; c_2 c_3 is clipped at 1e-12 and the clipped fraction is reported.
- Classification per window and level: `ENSTROPHY-DOMINATED` if
  E_enstrophy >= 2 E_gap and E_enstrophy > 0; `GAP-DOMINATED` if
  E_gap >= 2 E_enstrophy and E_gap > 0; otherwise `MIXED`.
- Overall verdict: the common class, otherwise `MIXED_ACROSS_WINDOWS`.

## Reviewer 4, comment 2: confirmation of the method (no new run)

- The "every second grid point" test (`dense.py`) decimates the velocity field by
  a factor of 2 in each direction without prior filtering.
- The gradients are then computed on the coarse grid: sixth-order differences,
  spacing 2 dx, margin 4, 56^3 points per cutout.
- The test is therefore a recomputation on a coarser grid, but without an
  anti-aliasing filter.
- The manuscript states the outcome as stability under the numerical variations
  examined, not as exclusion of subgrid effects.
