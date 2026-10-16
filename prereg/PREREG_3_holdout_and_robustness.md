# Analysis plan 3 — holdout replication, significance, references and domain

English rendering of the locked document (12 October 2026; checksum in
`NOTE_ON_CHECKSUMS.md`). Builds on analysis plan 2 (outcome REPLICATED). Locked
before any new data were downloaded. The analysis ran unattended from start to
finish; any deviation would have been reported as a deviation.

## 1. Frozen object

As in analysis plans 1-2: J[i<-j] = 2 c_i c_j (lambda_i - lambda_j)(1 - chi_ij),
predicted signs (-, +, -), O_ij >= 0.05, selection g > g_min with primary
g_min = 0.01, sixth-order differences, component order matched to axis order.

## 2. Data (downloaded after locking)

Seed 20261012 + N; origins rejected if within 128 points of a stage-1 or earlier
holdout cutout in all three directions.

holdout, Re_lambda 611 (isotropic4096, snapshot 1), 12 x 128^3:
(3275,3767,1618) (1764,3092,2503) ( 473,2180,1145) (2361,3506,2849)
(2293,3865,2074) ( 954,  30, 714) (3771,2531, 837) (2519,  43,2483)
(1305,1707, 262) (1646, 319, 663) (1251,3026,2596) (1860,1710,1899)

holdout, Re_lambda 1300 (isotropic8192, stored snapshot 4, 1-based; stage 1 used 1), 12 x 128^3:
(4449,4284,2942) (2478,3851,1593) (4982, 674,2432) (3613,1944,1588)
(6615,4099,1150) (2958, 767,6665) (6447,6200,6973) (4622,3833,4708)
( 861,4707,6583) (6176,7580,5400) ( 525,5502,5098) (4543, 357,5963)

extension, Re_lambda 433 (isotropic1024coarse, dx/eta 2.14), 12 x 128^3, time
index 1 (c00-c03), 2001 (c04-c07), 4001 (c08-c11):
( 698, 679, 259) (  64, 683, 391) ( 859, 439,  35) ( 246, 343, 324)
( 657, 813, 293) ( 349, 585,  27) ( 780, 213, 230) ( 881, 504, 406)
(  90, 536,  52) ( 457, 127, 614) ( 603, 161, 870) ( 864, 707, 591)
Same database as the quarantined cohort (other times and places); this window
is a boundary/extension window, not confirmatory.

channel (Re_tau 1000), velocity gradients from the database (fd6noint) at
points. Five bands y+ = 15, 30, 100, 300, 1000 (each +/- 10 %; centreline
|y| <= 0.1). 12 blocks per band, 2000 points per block in a 0.4 x 0.2 x-z patch
around the centre, lower wall (y = -1 + y+/Re_tau, Re_tau = 1000). Blocks
(t, x, z):
(1.0, 7.0611, 3.853) (3.0, 14.6972, 7.0311) (5.0, 3.4486, 2.5809)
(7.0, 1.2899, 3.0021) (9.0, 13.9097, 5.8508) (11.0, 23.0263, 1.5872)
(13.0, 6.4042, 6.7532) (15.0, 8.0422, 4.3598) (17.0, 17.5063, 3.1113)
(19.0, 13.3004, 2.8835) (21.0, 14.0915, 9.0557) (23.0, 14.6861, 0.8295)
Points in a block: default_rng(20261012 + 2000 + 100 * band + block).

The stage-1 cutouts of analysis plan 2 enter the pooled windows
(pooled_611 = stage 1 + holdout, pooled_1300 likewise) but not the holdout verdict
V1 or the primary significance V2.

## 3. Gates

Isotropic cutouts: G1-G3 of analysis plan 2; more than 2 of 12 failing makes the
window GATE_FAILURE. Channel: G2 only (G1 and G3 do not apply in wall-bounded
flow), limit 0.05, calibrated before locking on six test blocks with other seeds
and times than those above (G2 = 0.004-0.011 with fd6noint on the non-uniform
wall-normal grid); the limit is meant to catch a wrong reading (order one), not
discretisation error; the test blocks were discarded. More than 2 of 12 failing
makes the band GATE_FAILURE.
G0 (instrument): O(3) invariance on 10 000 real points with g > g_min (random
rotation and reflection of S and omega): max |dJ| / <|J|> < 1e-10 per edge.

## 4. Modules and verdicts

Block = cutout (isotropic) or block (channel). Uncertainty = block bootstrap
(10 000 resamples) of the ratio sum J / sum n. "Cell PASS" means everywhere:
predicted sign, 95 % block interval excluding zero, and O_ij >= 0.05.

V1 HOLDOUT (confirmatory). The procedure of analysis plan 2 verbatim (2000
points, seed 20260811, margin 32, g_min 0.01, unweighted mean) on the two holdout
windows alone. REPLICATED if all three signs and all three O_ij hold in both;
NOT_REPLICATED otherwise.

Dense analysis (V2-V7): all interior points (margin 8) with g > g_min.

V2 SIGNIFICANCE. Per edge and window (holdout 611, holdout 1300, extension 433):
one-sided block sign-flip test (exact when 2^n <= 1e6, otherwise 200 000 Monte
Carlo draws) of the pooled mean in the predicted direction. Family = 3 edges x 3
windows, Holm, alpha 0.01. SIGNIFICANT per cell if the Holm-adjusted p < 0.01.
Secondary (not deciding): binomial block sign test, t test on block means, 99 %
interval, pooled windows, effective sample size from the spatial correlation
length of J, heterogeneity (Cochran Q, I^2) of the normalised mean J/S_rms over
the three Reynolds numbers.

V3 REFERENCES. Three reference ensembles per block, same selection:
  G   isotropic Gaussian gradient (S and omega independent, <|S|^2> and
      <|omega|^2> matched per block): fully kinematic reference.
  xi  real S and |omega|, direction xi random isotropic: removes alignment.
  om  real S and xi, |omega|^2 permuted within the block: removes the coupling
      between enstrophy and the rest of the state.
Per edge and window: paired block difference (real - reference), bootstrap 95 %
interval. TURBULENCE-SPECIFIC if the interval for (real - G) excludes zero,
otherwise KINEMATIC. xi and om are reported as attribution and do not decide V3.
Interpretation fixed in advance: if G has the same sign but V3 is
TURBULENCE-SPECIFIC, the paper must state that the orientation is kinematic and
the magnitude/organisation turbulent. If V3 is KINEMATIC, the paper may not call
the current of that edge a turbulence phenomenon.

V4 RESOLUTION. Per window: (a) fourth- and eighth-order differences instead of
sixth; (b) every second point (2dx) with sixth order. ROBUST if all three edges
are cell PASS in all variants in both pooled windows; otherwise
RESOLUTION-SENSITIVE (edge and variant named). Relative change of the mean is
reported.

V5 (Q, R). Quadrants by the signs of Q = 1/2(|Omega|^2 - |S|^2) and R = -det A.
Per edge, quadrant and window: mean, interval, fraction of points, share of the
total. QR-ROBUST per edge if cell PASS in all four quadrants (with >= 500 points)
in both pooled windows; otherwise QR-PARTIAL with the quadrants named.

V6 TAILS. Per edge: mean after removing the 0.1 / 1 / 5 / 10 % largest |J| (from
50 000 points per block), median, and median of block means. TAIL-CARRIED if the
sign of the mean reverses when <= 1 % is removed, or if the median has the
opposite sign of the mean; otherwise BULK.

V7 DOMAIN (no extrapolation). Axes and levels:
  g_min                    0.001  0.003  0.01*  0.03  0.1  0.2
  filter Delta/eta         0* (unfiltered) 2 4 8 16 32 64 (Gaussian, sigma =
                           Delta/sqrt 12, interior margin ceil(3 sigma/dx) + 4;
                           fewer than 16^3 interior points = not evaluated)
  |omega|^2/<|omega|^2>    [0, 0.5) [0.5, 1) [1, 2) [2, 4) [4, inf)
  Re_lambda                433 611 1300
  dx/eta                   native and 2dx (from V4)
  y+ (channel)             15 30 100 300 1000
(* = reference level). Per edge: valid domain = the largest contiguous run of
levels around the reference where the cell passes. Enstrophy bins and y+ are
reported as the set of passing levels. The paper may only claim the result
within the valid domain.

V8 CHANNEL. Per band: cell PASS per edge with the full gradient (primary) and
with the band-mean gradient removed (secondary). HOLDS / FAILS per band. The
channel is a boundary test of generality; FAILS is a boundary finding, not a
falsification of the isotropic claim, and may not be rescued.

AUDIT. Static check of the analysis code for known defect classes (a superseded
Biot-Savart kernel, hash-based seeding, default grid spacing). CLEAN or FLAGGED.

## 5. Prohibitions

1. No change of parameters, gates, levels, seeds or verdicts after locking.
2. No rerun with other settings to change a verdict. A technical failure
   (network, memory) may be resumed with unchanged settings; it is logged.
3. A NOT_REPLICATED in V1 is not rescued by the extension window, the channel or
   pooling.
4. No post-hoc exclusion beyond the gates of section 3.
5. No interpretive framework beyond fluid mechanics; the pressure Hessian is not
   analysed.
6. The channel and the extension window may not be relabelled as "another
   regime" to excuse a FAIL; they are boundary measurements.
7. Disclosure of the self-test: before locking, the code was run on 4 + 4 stage-1
   cutouts (already seen data) to check that all modules run. Those outcomes are
   a code test only and are not cited. Because the pooled windows contain the
   stage-1 cutouts, V3-V7 harden a known result; the clean confirmatory part is
   V1 and V2 on the holdout windows, which nobody had seen.

## 6. Order (unattended)

1. This document and the analysis code are hash-locked. 2. The code verifies
both hashes, otherwise stop. 3. Audit and G0. 4. Download (retries up to 12 h).
5. Gates. 6. V1. 7. Dense analysis (V2-V7). 8. Channel (V8). 9. Statistics and
verdicts. 10. Report, figures, hash-logged result record.
