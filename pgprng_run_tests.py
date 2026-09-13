# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
Thorough statistical test run: applies the NIST-SP-800-22-style battery
(nist_battery.py) to BOTH generator designs at e1=101, e2=103, using the
ACTUAL delivered code (combined_prng_101_103.py) and the ACTUAL delivered
increment sets (pgprng_increments_101_103.csv -- additive increments, not
multipliers; see pgprng_common.py) -- not a re-derived toy version.

Methodology (matches real NIST STS practice): run each test across many
INDEPENDENT sequences (each freshly, genuinely randomly seeded via
secrets.randbits, exactly as real usage would seed it -- no fixed demo
seeds here), then for each test report:
  - proportion of sequences with p >= 0.01 (should fall within the
    standard NIST acceptable band around the sample size)
  - uniformity of the p-value distribution itself (chi-square over a
    10-bin histogram of p-values; report gammaincc(9/2, chi/2))
"""
import sys
import os
import time
import json
import numpy as np
from scipy.special import gammaincc

from pgprng_nist_battery import (
    words_to_bits, run_cheap_battery, test_binary_matrix_rank, test_linear_complexity,
)
from pgprng_legacy_clockstep import CombinedPRNG, SelectionCombinedPRNG

# ---------------------------------------------------------------
# load the actual delivered increment sets
# ---------------------------------------------------------------
INCS_101, INCS_103 = [], []
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pgprng_increments_101_103.csv")) as f:
    next(f)
    for line in f:
        ens, inc = line.strip().split(",")
        (INCS_101 if ens == "101" else INCS_103).append(int(inc))
assert len(INCS_101) == 101 and len(INCS_103) == 103
assert not (set(INCS_101) & set(INCS_103))

N_SEQ = 500
L_BITS = 1_000_000
L_WORDS = L_BITS // 64  # 15625

TEST_NAMES = ["monobit", "block_freq", "runs", "longest_run", "dft",
              "cusum_fwd", "cusum_rev", "approx_entropy", "serial_1", "serial_2",
              "rank", "linear_complexity"]


def gen_bits(gen, n_words):
    vals = np.empty(n_words, dtype=np.uint64)
    for i in range(n_words):
        vals[i] = gen.next()
    return words_to_bits(vals)


def full_battery(bits):
    out = run_cheap_battery(bits)
    p_rank, _ = test_binary_matrix_rank(bits)
    out["rank"] = p_rank
    p_lc, _ = test_linear_complexity(bits, M=500)
    out["linear_complexity"] = p_lc
    return out


def run_design(name, factory):
    print(f"\n{'='*70}\n{name}: {N_SEQ} independent sequences x {L_BITS:,} bits, "
          f"genuinely random-seeded (secrets.randbits)\n{'='*70}", flush=True)
    all_pvals = {t: [] for t in TEST_NAMES}
    t_start = time.time()
    for seq_i in range(N_SEQ):
        gen = factory()  # fresh secrets-random seeds each time -- real usage seeding
        bits = gen_bits(gen, L_WORDS)
        res = full_battery(bits)
        for t in TEST_NAMES:
            all_pvals[t].append(res[t])
        if (seq_i + 1) % 25 == 0:
            elapsed = time.time() - t_start
            print(f"  sequence {seq_i+1}/{N_SEQ}  ({elapsed:.1f}s elapsed)", flush=True)
    total_time = time.time() - t_start
    print(f"  done in {total_time:.1f}s", flush=True)
    return all_pvals


def summarize(name, all_pvals):
    print(f"\n--- {name}: summary over {N_SEQ} sequences ---")
    print(f"{'test':<18} {'pass_prop':>10} {'min_ok':>8} {'uniformity_p':>13} {'verdict':>10}")
    # NIST acceptable proportion band for alpha=0.01, m=N_SEQ
    alpha = 0.01
    p_hat = 1 - alpha
    band = 3 * np.sqrt(p_hat * (1 - p_hat) / N_SEQ)
    min_ok_frac = p_hat - band
    min_ok_count = min_ok_frac * N_SEQ
    summary = {}
    for t in TEST_NAMES:
        pv = np.array(all_pvals[t])
        n_pass = int(np.sum(pv >= 0.01))
        prop = n_pass / N_SEQ
        # uniformity: chi-square over 10 bins
        counts, _ = np.histogram(pv, bins=10, range=(0.0, 1.0))
        expected = N_SEQ / 10
        chi_sq = np.sum((counts - expected) ** 2 / expected)
        unif_p = gammaincc(9 / 2, chi_sq / 2)
        verdict = "OK" if (n_pass >= min_ok_count and unif_p >= 0.0001) else "FLAG"
        print(f"{t:<18} {prop:>10.3f} {min_ok_count:>8.1f} {unif_p:>13.5f} {verdict:>10}")
        summary[t] = dict(n_pass=n_pass, prop=prop, unif_p=float(unif_p), verdict=verdict,
                           mean_p=float(pv.mean()), pvals=pv.tolist())
    return summary, min_ok_count


if __name__ == "__main__":
    results = {}
    results["clock_stepping"] = run_design(
        "CLOCK-STEPPING (CombinedPRNG)", lambda: CombinedPRNG(INCS_101, INCS_103)
    )
    results["prng_selection"] = run_design(
        "PRNG-SELECTION (SelectionCombinedPRNG)", lambda: SelectionCombinedPRNG(INCS_101, INCS_103)
    )

    summaries = {}
    for name, pvals in results.items():
        summaries[name], min_ok = summarize(name, pvals)

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pgprng_nist_battery_results.json"), "w") as f:
        json.dump({"N_SEQ": N_SEQ, "L_BITS": L_BITS, "min_ok_count": min_ok,
                   "results": summaries}, f, indent=1)
    print("\nSaved pgprng_nist_battery_results.json (next to this script)")

    print(f"\n{'='*70}\nSIDE-BY-SIDE COMPARISON (clock-stepping vs PRNG-selection)\n{'='*70}")
    print(f"{'test':<18} {'clock pass%':>12} {'select pass%':>13} {'clock unif_p':>13} {'select unif_p':>14}")
    for t in TEST_NAMES:
        cs = summaries["clock_stepping"][t]
        ss = summaries["prng_selection"][t]
        print(f"{t:<18} {cs['prop']*100:>11.1f}% {ss['prop']*100:>12.1f}% "
              f"{cs['unif_p']:>13.5f} {ss['unif_p']:>14.5f}")
