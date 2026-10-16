# Analysis plan 2 — clean disjoint replication of the instantaneous current

English rendering of the locked document (11 August 2026; checksum in
`NOTE_ON_CHECKSUMS.md`). Bound by analysis plan 1 with g_min = 0.01. Locked
before any new data were downloaded or inspected.

## 0. Why this plan exists

The hypothesis of analysis plan 1 was first evaluated on two time windows of
isotropic1024coarse with a gap threshold of 0.05, before the selection rule had
been fixed. The analyst has seen that outcome. It is quarantined and may never
be cited. Changing to g_min = 0.01 does not rescue that cohort (the difference
affects less than 1 % of the points). This plan is therefore the confirmatory
test, on data not yet downloaded, with everything frozen in advance.

Provenance chain, to be stated in any report:
exploratory Lagrangian measurement -> quarantined observation -> this replication.

## 1. Data

Two independent DNS runs:

| window | dataset | Re_lambda | N | nu | eta | dx/eta |
|---|---|---|---|---|---|---|
| A | isotropic4096 | 611 | 4096 | 1.732e-4 | 1.3844e-3 | 1.11 |
| B | isotropic8192 | 1300 | 8192 | 4.385e-5 | 5.00e-4 | 1.53 |

The quarantined cohort (isotropic1024coarse, Re_lambda 433, dx/eta 2.14) is a
different realisation at a different Reynolds number and less well resolved.

Cutouts: 12 per window, 128^3, velocity only, stored snapshot 1. isotropic4096
has only one snapshot, so the uncertainty is a spatial block bootstrap in both
windows.

Origins, generated deterministically from default_rng(20260811 + N),
integers(1, N - 128 + 1, 3), written here before download:

isotropic4096:
(3776,3851,3102) (3074,1355,3727) (3786, 525, 479) ( 663,1010,1237)
( 622,3449,3340) (3236,3936,1127) (2389,3127,  46) (3329,1716, 217)
( 592,3440,1150) ( 318,2403, 545) (1734,2478, 105) (2171,3869,3121)

isotropic8192:
(2322,3143, 867) (4932,4272,6794) (6774, 466,7833) (5017,2465,6944)
( 355,4403,3939) (4432,6810,3255) (2465,4028,4256) (7036,3870,1643)
(2080,1810,6660) (6587,7647,5460) (3893,3784,1030) (5106,3254,7070)

Sampling points: 2000 per cutout, default_rng(20260811).integers(32, 96,
size=(4000, 3))[:2000], identical to the quarantined cohort, so that only the
data differ. 24 000 points per window.

## 2. Frozen object and convention

As in analysis plan 1. g_min = 0.01, unchanged.

## 3. Frozen predictions

Signs (-, +, -) on (1<-2, 1<-3, 2<-3); O_ij >= 0.05 on every edge. Primary
aggregate: unweighted mean over eligible points; median and occupation-weighted
mean are secondary. Uncertainty: block bootstrap over the 12 cutouts, 2000
resamples; the confidence interval is reported but the outcome is decided by
signs and O_ij.

REPLICATED: all three signs and all three O_ij hold in both windows.
NOT_REPLICATED: otherwise. No intermediate category.

## 4. Orientation gate (hard stop)

Background: in an earlier, unrelated analysis spatial axes were permuted
without permuting the velocity components, which corrupts the field. The edge
current is fully O(3)-invariant (verified before locking on 20 000 synthetic
tensors: proper rotation 1.7e-14, reflection 8.6e-15, inconsistent permutation
0.92), so the gate must detect corruption, not handedness.

Per cutout, all three must hold:

    G1  <|S|^2> / <|omega|^2>   in [0.45, 0.55]
    G2  <|div u|> / <|grad u|>  <= 1.0e-2
    G3  longitudinal skewness   in [-0.70, -0.30]

A failing cutout is excluded and reported. If more than 2 of 12 fail in a
window, the window stops with outcome GATE_FAILURE.

## 5. Interpretation limits, written before the outcome

The Reynolds numbers differ from the quarantined cohort (611 and 1300 against
433), which makes the test harder. If REPLICATED, the orientation is robust over
a factor of three in Re_lambda and two independent DNS. If NOT_REPLICATED, the
outcome is NOT_REPLICATED; the Re dependence of chi may be reported as an
observation but may not be used to relabel the outcome. The sign of the
stretching part is algebra; what is tested is whether the turbulent
distribution of chi gives the same edge orientation.

## 6. Prohibitions

1. The quarantined outcome may not be opened, consulted or compared before this
   test has been run and judged.
2. No change of g_min, sampling points, origins, seed, convention, signs, O limit
   or gate limits after locking.
3. No new features, coordinates, weightings or forecast models.
4. No claims outside fluid mechanics.
5. If the test fails, no third data set is sought.
6. No post-hoc exclusion of cutouts beyond the three gates.
7. New points in existing fields may be used as additional internal robustness
   after this test, never as the test the main claim rests on.

## 7. Order

1. This document is hash-locked. 2. Download of the 24 cutouts, no analysis
during download. 3. Orientation gate, reported before the currents are computed.
4. Currents and verdict, one pass, no variants.
