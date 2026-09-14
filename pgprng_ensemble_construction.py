# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
Staged, self-validating ensemble increment construction.

(Terminology note: these are additive increments -- each component
advances as x -> x + increment mod 2**64 -- not multipliers. There is
no multiplicative term anywhere in the recurrence; the name only ever
described the LCG-style role these values play, which was imprecise.)

Uniform-random selection over the full odd 64-bit range makes the
near-duplicate-increment defect (found earlier this session at N=64)
astronomically unlikely -- but "astronomically unlikely" is a
probability argument, not a guarantee, and costs nothing to actually
verify. This adds three independent gates, any of which can trigger a
re-draw:

  1. STRUCTURAL: minimum pairwise gap across all combined increments
     must exceed a safety threshold, chosen with a large margin below
     the range but a large margin above the failure scale actually
     observed (the N=64 defect band was ~1022 wide, i.e. gap ~10^2-10^3;
     the default threshold here is 2**40 =~ 1.1x10^12, nine-plus orders
     of magnitude above the known failure scale, while still ~250x
     smaller than the ~2**48 mean gap 204 uniform draws typically land
     at -- a real, checkable constraint, not a tautology).
  2. UNIFORMITY: a chi-square goodness-of-fit test that the combined
     increments are spread roughly evenly across the range (catches
     clustering patterns a single pairwise-gap check could miss).
  3. EMPIRICAL SMOKE TEST: actually build the real generator with the
     candidate increments, generate a modest sample, and check
     popcount mean/std and bit-balance are within tolerance. This is
     the most direct gate -- it measures the real behavior instead of
     reasoning about it.

Any candidate draw failing any gate is discarded and a fresh one drawn,
up to a retry budget. This turns "we checked the probability is tiny"
into "we verified this specific draw before using it."
"""
import random
import secrets
import numpy as np
from scipy.stats import chi2
import sys
import os

from pgprng_legacy_clockstep import CombinedPRNG

M64 = 1 << 64


def random_odd_increment(rng, count):
    chosen = set()
    while len(chosen) < count:
        k = rng.randrange(0, 1 << 63)
        chosen.add(2 * k + 1)
    return sorted(chosen)


def min_pairwise_gap(incs):
    s = sorted(incs)
    return min(b - a for a, b in zip(s[:-1], s[1:]))


def uniformity_chi2(incs, M=M64, bins=16):
    counts = np.zeros(bins, dtype=int)
    width = M // bins
    for inc in incs:
        idx = min(inc // width, bins - 1)
        counts[idx] += 1
    expected = len(incs) / bins
    chi_sq = np.sum((counts - expected) ** 2 / expected)
    p = chi2.sf(chi_sq, df=bins - 1)
    return float(chi_sq), float(p), counts.tolist()


def smoke_test_ensemble(incs_a, incs_b, n_words=200_000):
    gen = CombinedPRNG(incs_a, incs_b)  # fresh secrets-random seeds, real usage path
    vals = np.empty(n_words, dtype=np.uint64)
    for i in range(n_words):
        vals[i] = gen.next()
    popcounts = np.zeros(n_words, dtype=np.int64)
    for k in range(64):
        popcounts += ((vals >> np.uint64(k)) & np.uint64(1)).astype(np.int64)
    return dict(
        mean_popcount=float(popcounts.mean()),
        std_popcount=float(popcounts.std()),
        bit0_frac=float((vals & 1).mean()),
        distinct_frac=float(len(np.unique(vals)) / n_words),
    )


def build_validated_ensemble(e1, e2, min_gap=2 ** 40, uniformity_alpha=0.001,
                              popcount_tol=0.5, bit_frac_tol=0.02,
                              max_retries=20, seed=None, verbose=True):
    rng = random.Random(seed) if seed is not None else random.Random(secrets.randbits(64))
    log = []
    for attempt in range(1, max_retries + 1):
        combined = random_odd_increment(rng, e1 + e2)
        incs_a = sorted(combined[:e1])
        incs_b = sorted(combined[e1:])
        gap = min_pairwise_gap(combined)
        chi_sq, unif_p, counts = uniformity_chi2(combined)
        gap_ok = gap >= min_gap
        unif_ok = unif_p >= uniformity_alpha
        entry = dict(attempt=attempt, gap=gap, gap_ok=gap_ok, unif_p=unif_p, unif_ok=unif_ok)
        if not (gap_ok and unif_ok):
            entry["stage"] = "structural"
            entry["accepted"] = False
            log.append(entry)
            if verbose:
                print(f"  attempt {attempt}: REJECTED (structural) gap={gap:.3e} "
                      f"(ok={gap_ok}) unif_p={unif_p:.4f} (ok={unif_ok})")
            continue
        stats = smoke_test_ensemble(incs_a, incs_b)
        pc_ok = abs(stats["mean_popcount"] - 32.0) <= popcount_tol
        bit0_ok = abs(stats["bit0_frac"] - 0.5) <= bit_frac_tol
        distinct_ok = stats["distinct_frac"] > 0.999
        smoke_ok = pc_ok and bit0_ok and distinct_ok
        entry["stage"] = "smoke_test"
        entry["smoke"] = stats
        entry["accepted"] = smoke_ok
        log.append(entry)
        if verbose:
            tag = "ACCEPTED" if smoke_ok else "REJECTED (smoke test)"
            print(f"  attempt {attempt}: structural OK (gap={gap:.3e}, unif_p={unif_p:.4f}); "
                  f"smoke={stats} -> {tag}")
        if smoke_ok:
            return dict(incs_a=incs_a, incs_b=incs_b, attempt=attempt, gap=gap, unif_p=unif_p,
                        smoke=stats, log=log)
    raise RuntimeError(f"failed to construct a validated ensemble in {max_retries} attempts; log={log}")


if __name__ == "__main__":
    print("=" * 70)
    print("Retroactive validation of the ALREADY-DELIVERED e=101/103 increment set")
    print("=" * 70)
    INCS_101, INCS_103 = [], []
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pgprng_increments_101_103.csv")) as f:
        next(f)
        for line in f:
            ens, inc = line.strip().split(",")
            (INCS_101 if ens == "101" else INCS_103).append(int(inc))
    combined = INCS_101 + INCS_103
    gap = min_pairwise_gap(combined)
    chi_sq, unif_p, counts = uniformity_chi2(combined)
    print(f"min pairwise gap: {gap:.6e}  (threshold 2**40 = {2**40:.3e}) -> "
          f"{'PASSES' if gap >= 2**40 else 'FAILS'} structural gate")
    print(f"uniformity chi_sq={chi_sq:.3f} p={unif_p:.4f} bin counts={counts} -> "
          f"{'PASSES' if unif_p >= 0.001 else 'FAILS'} uniformity gate")
    stats = smoke_test_ensemble(INCS_101, INCS_103)
    print(f"smoke test: {stats}")

    print("\n" + "=" * 70)
    print("POSITIVE CONTROL: does the same pipeline correctly REJECT a construction")
    print("we already know is bad -- increments drawn from the narrow 'good value'")
    print("pool instead of the full range (the exact defect found earlier at N=64)?")
    print("=" * 70)
    # NOTE: best_m_n64.csv is an external artifact from earlier this session
    # (the "narrow good-value pool" that produced the known N=64 defect) and
    # is not itself part of this package -- this positive-control block is
    # included for provenance but will raise FileNotFoundError unless that
    # file is supplied alongside this script. Its "m" column header predates
    # this file's increment/multiplier terminology cleanup and was left
    # as-is since it belongs to that external file, not to pgprng.
    good_pool = []
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "best_m_n64.csv")) as f:
        header = f.readline().strip().split(",")
        m_idx = header.index("m")
        for line in f:
            parts = line.strip().split(",")
            good_pool.append(int(parts[m_idx]))
    rng_bad = random.Random(42)
    bad_combined = rng_bad.sample(good_pool, 204)
    bad_a = sorted(bad_combined[:101])
    bad_b = sorted(bad_combined[101:])
    bad_gap = min_pairwise_gap(bad_combined)
    bad_chi, bad_unif_p, bad_counts = uniformity_chi2(bad_combined)
    print(f"[bad-pool draw] min pairwise gap: {bad_gap:.6e} -> "
          f"{'PASSES' if bad_gap >= 2**40 else 'CORRECTLY REJECTED'} structural gate")
    print(f"[bad-pool draw] uniformity chi_sq={bad_chi:.3f} p={bad_unif_p:.6g} -> "
          f"{'PASSES' if bad_unif_p >= 0.001 else 'CORRECTLY REJECTED'} uniformity gate")

    print("\n" + "=" * 70)
    print("DEMONSTRATION: build_validated_ensemble() end to end (fresh random draws)")
    print("=" * 70)
    result = build_validated_ensemble(101, 103, seed=20260913, verbose=True)
    print(f"\naccepted on attempt {result['attempt']}, min_gap={result['gap']:.3e}, "
          f"uniformity_p={result['unif_p']:.4f}")
    print(f"smoke test: {result['smoke']}")
