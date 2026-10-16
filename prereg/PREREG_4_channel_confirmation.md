# Analysis plan 4 — preregistered confirmation of the channel observation

English rendering of the locked document (12 October 2026; checksum in
`NOTE_ON_CHECKSUMS.md`). Builds on analysis plan 3 and its outcome. Locked before
data were downloaded.

## 1. Background (descriptive observation, not preregistered in plan 3)

Normalised by the band's local S_rms = sqrt(<S:S>) (full gradient), the channel
values of <J>/S_rms at y+ >= 100 were at the isotropic level and were reduced at
y+ = 15-30. This plan tests the observation on new blocks.

## 2. Data

JHTDB channel, same method as plan 3 (fd6noint gradients at points, 2000 points
per block, 0.4 x 0.2 x-z patch, lower wall, Re_tau = 1000). Five bands y+ = 15,
30, 100, 300, 1000 (+/- 10 %; centreline |y| <= 0.1). 12 new blocks per band:
time t_k = 2.0 + 2.0 k (k = 0..11; plan 3 used odd times), x and z centres from
default_rng(20261013 + 1000) as in plan 3. Points from
default_rng(20261013 + 2000 + 100 * band + block).

## 3. Reference (frozen from plan 3, mean of the two pooled isotropic windows, J/S_rms)

    r = (1<-2, 1<-3, 2<-3) = (-0.1742, +0.0474, -0.0154)

## 4. Gate

G2 = <|tr A|>/<|A|> <= 0.05 per block; more than 2 of 12 failing makes the band
GATE_FAILURE.

## 5. Measure

Per band and edge: J/S_rms = (sum J / n) / sqrt(sum S:S / n) over blocks with
g > 0.01 (S_rms over all points); 95 % interval by block bootstrap (10 000).

## 6. Criteria

K1 OUTER AGREEMENT: for y+ = 100, 300, 1000 and the edges 1<-2 and 1<-3 the 95 %
   interval must overlap [0.75 r, 1.25 r]. All six -> K1 PASS.
K2 WALL REDUCTION: at y+ = 15 the upper bound of the |J/S_rms| interval must be
   below 0.5 |r| for both 1<-2 and 1<-3 -> K2 PASS.
VERDICT: CONFIRMED if K1 and K2 pass; PARTIAL if only one; NOT_CONFIRMED otherwise.
2<-3 and y+ = 30 are reported descriptively (2<-3 is tail-carried, plan 3 V6).

## 7. Prohibitions

No change after locking; no rerun; a NOT_CONFIRMED means that the channel
normalisation may only be reported as a descriptive observation.
