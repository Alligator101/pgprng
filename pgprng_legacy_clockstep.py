# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
Same two designs (CombinedPRNG = clock-stepping, SelectionCombinedPRNG
= mutual PRNG-selection), enlarged from e=11/e=13 to e=101/e=103 --
twin primes near 100, guaranteeing (since coprime) the full benefit of
multiplying the two component PERIODS together for the clock-stepping
design, rather than a gcd-reduced one (see this session's period
analysis). Note: "multiplying periods together" refers to how the two
ensembles' periods combine -- it has nothing to do with the per-component
recurrence itself, which is purely additive (see pgprng_common.py: these
are increments, not multipliers).

No changes needed to the class internals -- both classes are already
parameterized by e via len(incs). "Enlarging the ensembles" is just:
draw 101+103=204 distinct random odd increments (uniformly over the full
64-bit odd range, as settled on earlier), split them, and instantiate.

Estimated periods (see this session's toy-scale verification):
  - CombinedPRNG (clock-stepping):        EXACT  = 2**64 * 101 * 103
  - SelectionCombinedPRNG (PRNG-select):   ORDER OF MAGNITUDE ~ 2**128
    (extrapolated from a toy-scale pattern -- cycle length == M**2,
    confirmed exactly across 10 diverse (N, e1, e2) configurations
    with zero exceptions -- not proven at full scale)
"""
import random

from pgprng_common import MASK64, mix64, validate_increments, validate_disjoint, prepare_seeds


class RotatingEnsemble:
    def __init__(self, incs, seeds=None):
        incs = validate_increments(incs, "incs")
        seeds = prepare_seeds(seeds, incs, "seeds")
        self.incs = incs
        self.e = len(incs)
        self.states = list(seeds)
        self.cursor = 0

    def next(self):
        i = self.cursor
        self.states[i] = (self.states[i] + self.incs[i]) & MASK64
        out = self.states[i]
        self.cursor = (self.cursor + 1) % self.e
        return out


class CombinedPRNG:
    """Clock-stepping (fixed round-robin) design. Exact period:
    2**64 * e1 * e2 when e1, e2 are coprime (proven, and confirmed
    exactly in every toy-scale test this session)."""

    def __init__(self, incs_11, incs_13, seeds_11=None, seeds_13=None):
        validate_disjoint(incs_11, incs_13, "incs_11", "incs_13")
        self.ens11 = RotatingEnsemble(incs_11, seeds_11)
        self.ens13 = RotatingEnsemble(incs_13, seeds_13)

    def next(self):
        out11 = self.ens11.next()
        out13 = self.ens13.next()
        return mix64(out11 ^ out13)

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()


class SelectionCombinedPRNG:
    """Mutual PRNG-selection design. Period not exactly calculable;
    order-of-magnitude estimate ~2**128 from the toy-scale M**2
    pattern found this session (see module docstring)."""

    def __init__(self, incs_11, incs_13, seeds_11=None, seeds_13=None):
        incs_11 = validate_increments(incs_11, "incs_11")
        incs_13 = validate_increments(incs_13, "incs_13")
        validate_disjoint(incs_11, incs_13, "incs_11", "incs_13")

        self.incs_11 = incs_11
        self.incs_13 = incs_13
        self.e11 = len(incs_11)
        self.e13 = len(incs_13)
        self.states_11 = prepare_seeds(seeds_11, incs_11, "seeds_11")
        self.states_13 = prepare_seeds(seeds_13, incs_13, "seeds_13")

        self.last_out11 = 0
        for s in self.states_11:
            self.last_out11 ^= s
        self.last_out13 = 0
        for s in self.states_13:
            self.last_out13 ^= s

        self.select_counts_11 = [0] * self.e11
        self.select_counts_13 = [0] * self.e13

    def next(self):
        idx13 = self.last_out11 % self.e13
        self.states_13[idx13] = (self.states_13[idx13] + self.incs_13[idx13]) & MASK64
        new_out13 = self.states_13[idx13]
        self.select_counts_13[idx13] += 1

        idx11 = self.last_out13 % self.e11
        self.states_11[idx11] = (self.states_11[idx11] + self.incs_11[idx11]) & MASK64
        new_out11 = self.states_11[idx11]
        self.select_counts_11[idx11] += 1

        self.last_out11 = new_out11
        self.last_out13 = new_out13

        return mix64(new_out11 ^ new_out13)

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()


def random_odd_increment(rng, count):
    chosen = set()
    while len(chosen) < count:
        k = rng.randrange(0, 1 << 63)
        chosen.add(2 * k + 1)
    return sorted(chosen)


SELECTION_SEED = 20260913  # same documented seed as before, now drawing 204 not 24
E1, E2 = 101, 103

if __name__ == "__main__":
    picker = random.Random(SELECTION_SEED)
    combined_incs = random_odd_increment(picker, E1 + E2)
    INCS_101 = sorted(combined_incs[:E1])
    INCS_103 = sorted(combined_incs[E1:])
    assert not (set(INCS_101) & set(INCS_103))
    assert len(INCS_101) == E1 and len(INCS_103) == E2

    for name, incs in (("101-side", INCS_101), ("103-side", INCS_103)):
        residues = set(inc % 8 for inc in incs)
        print(f"{name}: {len(incs)} increments, residues mod 8 present: {sorted(residues)}")

    import os
    _here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(_here, "pgprng_increments_101_103.csv"), "w") as f:
        f.write("ensemble,increment\n")
        for inc in INCS_101:
            f.write(f"101,{inc}\n")
        for inc in INCS_103:
            f.write(f"103,{inc}\n")
    print("wrote pgprng_increments_101_103.csv (next to this script)")

    # =================================================================
    # Statistical comparison, same battery as the e=11/13 run, same
    # increments and seeds for both designs so the comparison isolates
    # the stepping-rule effect alone.
    # =================================================================
    import numpy as np

    FIXED_SEEDS_101 = [i * 13 + 3 for i in range(E1)]  # reproducible test seeds ONLY
    FIXED_SEEDS_103 = [i * 13 + 3 for i in range(E2)]  # -- never use for real output
    N_SAMPLES = 5_000_000

    def run_battery(gen, n=N_SAMPLES):
        vals = np.empty(n, dtype=np.uint64)
        for i in range(n):
            vals[i] = gen.next()

        popcounts = np.zeros(n, dtype=np.int64)
        for k in range(64):
            popcounts += ((vals >> np.uint64(k)) & np.uint64(1)).astype(np.int64)

        bit0 = (vals & 1).astype(np.int64)
        bit1 = ((vals >> 1) & 1).astype(np.int64)
        bit2 = ((vals >> 2) & 1).astype(np.int64)

        vf = vals.astype(np.float64)
        v_next = np.roll(vf, -1)
        corr1 = np.corrcoef(vf[:-1], v_next[:-1])[0, 1]

        distinct = len(np.unique(vals))

        return {
            "popcount_mean": popcounts.mean(),
            "popcount_std": popcounts.std(),
            "bit0_frac": bit0.mean(),
            "bit1_frac": bit1.mean(),
            "bit2_frac": bit2.mean(),
            "lag1_corr": corr1,
            "distinct_frac": distinct / n,
        }

    clock_gen = CombinedPRNG(INCS_101, INCS_103, seeds_11=list(FIXED_SEEDS_101), seeds_13=list(FIXED_SEEDS_103))
    select_gen = SelectionCombinedPRNG(INCS_101, INCS_103, seeds_11=list(FIXED_SEEDS_101), seeds_13=list(FIXED_SEEDS_103))

    clock_stats = run_battery(clock_gen)
    select_stats = run_battery(select_gen)

    print(f"\n{'='*70}\nCOMPARISON over {N_SAMPLES:,} samples each, e1={E1} e2={E2}, identical increments and seeds\n{'='*70}")
    print(f"{'metric':<18} {'clock-step':>14} {'PRNG-select':>14} {'ideal':>10}")
    ideals = {"popcount_mean": 32.0, "popcount_std": 4.0, "bit0_frac": 0.5,
              "bit1_frac": 0.5, "bit2_frac": 0.5, "lag1_corr": 0.0, "distinct_frac": 1.0}
    for key in clock_stats:
        print(f"{key:<18} {clock_stats[key]:>14.5f} {select_stats[key]:>14.5f} {ideals[key]:>10.4f}")

    counts11 = np.array(select_gen.select_counts_11)
    counts13 = np.array(select_gen.select_counts_13)
    exp11 = N_SAMPLES / select_gen.e11
    exp13 = N_SAMPLES / select_gen.e13
    print(f"\nselection fairness (PRNG-selection design only):")
    print(f"  101-side: expected {exp11:.0f} picks/component, "
          f"observed min={counts11.min()} max={counts11.max()} std={counts11.std():.1f}")
    print(f"  103-side: expected {exp13:.0f} picks/component, "
          f"observed min={counts13.min()} max={counts13.max()} std={counts13.std():.1f}")
