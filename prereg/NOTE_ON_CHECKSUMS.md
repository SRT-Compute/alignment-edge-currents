# Note on checksums and on what this package is

## What these documents are

The documents in this folder are **English renderings** of analysis plans that
were written in Danish and hash-locked, in an append-only hash log, before the
corresponding data were downloaded. They are not themselves the locked
artefacts. The locked originals exist unchanged; they also contain references to
unrelated internal projects. They are available confidentially to the editor
and reviewers on request, and the checksums below allow them to be verified.

## Checksums of the locked originals

| Locked original | Rendered here as | SHA-256 |
|---|---|---|
| analysis plan, hypotheses on the edge currents | `PREREG_1_original_hypothesis.md` | `66d4abc703ba376d2bdaeecc518898c48c674179f303b10d47c7d6129b12ed34` |
| instrument qualification, main document | (selection rule quoted in plan 1, section 5) | `1ceec5a4c2f1438f1c3c2db39bcb909df7a11559a12f555e272f065876d4a451` |
| instrument qualification, addendum | (not rendered; concerns the other terms) | `1ce57cdc29fc4823f367b16a6423a669fefcd4a141e6ae411fc52e8bcd252df3` |
| clean disjoint replication | `PREREG_2_replication.md` | `d9b97962b6f532a39e6d6a80cc30a08da3e8e6eaca65b59f614cbe8ab73008ce` |
| holdout, significance, references, domain | `PREREG_3_holdout_and_robustness.md` | `a49548b16da783321cf65194674b8ab43a83bf989c3715a509fcbb2196945a8a` |
| channel confirmation | `PREREG_4_channel_confirmation.md` | `4934731f91a410290569873f30ee9c922a3107adf717a7210cf3512050b872ea` |
| analyses added in revision | `PREREG_5_revision_analyses.md` | `c74069634c93c065133340eef809304f06511287ab59188e0d9187cebe79fb6d` |
| analysis script of the revision analyses | `code/revision_analyses.py` (see below) | `1e48f35b4fa86524266781c7255420b256419f507acf9297341b0805dff8a7e4` |

The analysis code that produced the original results read the checksum of its
own analysis plan at run time and refused to run if it did not match; the
checksums were also written into the result records.

The revision plan and its analysis script were hash-locked together before the
analyses were run. `code/revision_analyses.py` is the locked script with comments
and messages in English and with the input and output locations given as
arguments instead of fixed local paths; the computations, constants, seeds and
order of random draws are unchanged. It reproduces the locked run's results
exactly; `results/revision_analyses.json` was written by it.

## What was changed in preparing the renderings

1. **Language.** Danish to English.
2. **Nomenclature.** The substitutions listed in `NOMENCLATURE.md`, applied
   consistently. Every substitution is a renaming; none is a redefinition.
3. **Scope.** Sections of the first plan that define hypotheses this study does
   not report are summarised rather than rendered in full. Nothing that gates a
   reported result has been shortened.
4. **Internal references.** File names, project labels and references to
   unrelated internal projects are omitted; documents are identified by checksum.

Not changed: any data specification, origin, seed, threshold, gate, null model,
bootstrap procedure, verdict rule, prohibition, or frozen sign.

## Outcomes recorded for hypotheses outside this study

So that no reader supposes a favourable subset was selected: the first plan also
contains hypotheses on the attribution of the compensating terms and on the
correspondence with a coarse-grained transition current. The attribution
hypothesis was stopped at its preregistered gate. None of these outcomes bears
on any claim of the study, and none is used as evidence.
