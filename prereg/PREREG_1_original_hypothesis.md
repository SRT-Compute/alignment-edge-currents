# Analysis plan 1 — original hypothesis on the edge currents

English rendering of the part of the original analysis plan that concerns this
study (locked 11 August 2026; checksum in `NOTE_ON_CHECKSUMS.md`). The plan was
locked before any outcome on the Johns Hopkins Turbulence Databases (JHTDB) was
opened.

## 0. Declared sequence

Before the plan was written, the mean edge currents had been computed on a
Lagrangian data set (TURB-Lagr). That computation is an exploratory outcome
measurement, not a calibration. It was used only to formulate the hypothesis
below as a replication hypothesis, and it may never be cited as preregistered or
confirmatory.

Exploratory values (TURB-Lagr, for reference only):

| edge | mean current | stretching part alone |
|---|---|---|
| 1<-2 | -6.37 | +2.31 |
| 1<-3 | +1.75 | +2.61 |
| 2<-3 | -0.59 | +2.05 |

Net per weight: out of 1 (-4.62), into 2 (+5.78), out of 3 (-1.16).
Median chi on the edges 1-2, 1-3, 2-3: 1.67, 0.18, 0.45.

## 1. Convention (frozen)

J[i<-j] is the current into alignment weight i from weight j;
dc_i/dt = sum_{j != i} J[i<-j]. A positive J[1<-2] means net transfer 2 -> 1,
a negative value net transfer 1 -> 2. Eigenvalues are ordered
lambda_1 > lambda_2 > lambda_3.

## 2. Hypothesis: replication of the instantaneous current

Object:

    J[i<-j] = 2 c_i c_j (lambda_i - lambda_j)(1 - chi_ij),
    chi_ij  = |omega|^2 / (4 (lambda_i - lambda_j)^2).

No time derivative, no pressure, no forcing.

Frozen predictions, every window:

    (a) <J[1<-2]> < 0      net 1 -> 2
    (b) <J[1<-3]> > 0      net 3 -> 1
    (c) <J[2<-3]> < 0      net 2 -> 3

that is, the oriented pattern 1 -> 2 -> 3 -> 1.

Organisation requirement (preregistered minimal flux fraction):
O_ij = |<J_ij>| / <|J_ij|> >= 0.05 on every edge. Reason: the net current must
not be a vanishing remainder of large cancelling point values. The number was
chosen without having been measured.

Aggregation (one primary): unweighted mean over all eligible points. Secondary
(not deciding): occupation-weighted mean and median.

Uncertainty: block bootstrap, not point-wise N.

Outcome: REPLICATED if all three signs and all three O requirements hold in
both windows; otherwise NOT_REPLICATED.

Note on strength (in the original): the hypothesis replicates an algebraically
composed quantity. The direction of the stretching part is algebra and chi is a
closed form. What is tested is whether the actual turbulent distribution of chi
gives the same edge orientation in another flow.

## 3. Bookkeeping expectation for the remaining terms (not a hypothesis)

In global statistical stationarity the remaining terms of the alignment dynamics
(pressure Hessian, viscosity, forcing) must on average supply minus the net
drift of J, weight by weight. Expected signs of their contribution: positive
into weight 1, negative into weight 2, positive into weight 3 (equivalently, a
net drift of J out of 1, into 2, out of 3). The plan classifies this as a
bookkeeping check that follows from stationarity, not as a finding.

## 4. Out of scope

The plan contains further hypotheses on the attribution of the compensating
terms and on the correspondence with a coarse-grained transition current. They
are outside the scope of the present study and are not reported (see
`NOTE_ON_CHECKSUMS.md`).

## 5. Data and selection (as applied to this hypothesis)

The original plan applied the hypothesis to two time windows of
isotropic1024coarse. That application was run before the gap selection rule had
been fixed and has been quarantined; it is not used anywhere (see analysis
plan 2). The selection rule fixed afterwards in the instrument qualification is
g = min over edges of (lambda_i - lambda_j)/|S|_F > g_min with g_min = 0.01.

## 6. Prohibitions

1. No new coordinates, features or forecast models.
2. No claims outside fluid mechanics.
3. The exploratory TURB-Lagr result may not be cited as confirmatory.
4. No change of g_min, sampling points, convention or thresholds after locking.
5. If the hypothesis fails, no third window is sought to try again.
