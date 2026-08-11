# 0002 — Default Tanimoto cutoff is 0.40, not 0.85

**Status:** Accepted · **Supersedes:** the `tanimoto >= 0.85` value in the original blueprint

## Problem

The blueprint specified a similarity cutoff of `0.85`. That number is intuitive —
"85% similar sounds like a strong match" — and it is wrong by roughly a factor of
two for this fingerprint. Shipped as-is it would produce an application that
returns empty results for almost every query while appearing to work correctly.

Tanimoto values are **not comparable across fingerprint types**. Path-based and
MACCS keys produce much higher coefficients than circular fingerprints for the
same pair of molecules. Morgan/ECFP4 is sparse and specific, so its coefficients
run low: 0.85 on ECFP4 means "differs by about a methyl group."

## Evidence

Measured with this repo's own pipeline (`Morgan r=2, 2048-bit`, after
standardisation), on pairs whose pharmacological relationship is not in dispute:

| Pair | Tanimoto | Relationship |
|---|---|---|
| imatinib ↔ nilotinib | **0.517** | Both BCR-ABL inhibitors, heavily overlapping kinase off-target profiles |
| aspirin ↔ salicylic acid | **0.448** | Salicylic acid *is* aspirin's active metabolite |
| ibuprofen ↔ naproxen | **0.421** | Same NSAID class, same primary COX target |
| aspirin ↔ paracetamol | 0.222 | Both analgesics, mechanistically unrelated |
| caffeine ↔ anything | ≤ 0.098 | Unrelated |

Pairs surviving each cutoff, out of 28 tested:

| Cutoff | Pairs surviving |
|---|---|
| 0.85 | **0** |
| 0.70 | **0** |
| 0.55 | **0** |
| 0.40 | 3 |
| 0.30 | 3 |

Imatinib/nilotinib is the textbook polypharmacology example — the exact case this
platform exists to surface. At the blueprint's 0.85 it is invisible.

Note the clean separation: at 0.40 the three surviving pairs are precisely the
three genuine pharmacological relatives, and the highest false pair sits at 0.226.
The gap between 0.42 and 0.23 is where the signal lives.

## Decision

- Default cutoff: **0.40** (`DEFAULT_TANIMOTO_CUTOFF`).
- UI slider range: **0.20 – 1.00**, defaulting to 0.40, with the 0.35–0.55 band
  marked as the useful region.
- The API returns the coefficient on every edge so callers can re-filter without
  re-querying.

## Consequences

- Recall improves enormously; precision drops. This is the correct trade for a
  *hypothesis-generation* tool — a missed off-target is a silent failure, a
  spurious one is visible and dismissible.
- Because similarity is computed at query time ([0001](0001-query-time-similarity.md)),
  the threshold is a live parameter. Users who want 0.85 can drag the slider;
  nothing needs re-ingesting.

## Caveat that belongs in the UI

Chemical similarity implies *possible* shared binding, not confirmed binding, and
absence of a reported interaction is not evidence of absence — ChEMBL is heavily
biased toward well-studied target families (kinases, GPCRs). Output ranks
hypotheses for follow-up; it does not predict safety.
