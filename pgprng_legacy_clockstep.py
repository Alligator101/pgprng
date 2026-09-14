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
parameterized by ensemble size via len(ensemble). "Enlarging the
ensembles" is just: draw 101+103=204 distinct random odd increments
(uniformly over the full 64-bit odd range, as settled on earlier),
split them, and instantiate.

Object model: each ensemble is a list of State objects (see
pgprng_common.py) -- one per component, bundling its fixed increment
with its current running state, rather than the two parallel lists
(one of increments, one of states) used before this pass. Ensemble
size is always read live via len(ensemble); nothing caches it
separately, since ensemble sizes are meant to be mutable rather than
baked into the code the way the old e=11/e=13-derived attribute names
(incs_11/incs_13, ens11/ens13) implied a fixed scale that had already
stopped being true once this project grew to 101/103 components.

Estimated periods (see this session's toy-scale verification):
  - CombinedPRNG (clock-stepping):        EXACT  = 2**64 * 101 * 103
  - SelectionCombinedPRNG (PRNG-select):   ORDER OF MAGNITUDE ~ 2**128
    (extrapolated from a toy-scale pattern -- cycle length == M**2,
    confirmed exactly across 10 diverse (N, e1, e2) configurations
    with zero exceptions -- not proven at full scale)
"""
import random

from pgprng_common import mix64, validate_disjoint, build_ensemble


class RotatingEnsemble:
    """One ensemble, stepped in a fixed round-robin ("clock-stepping")
    order. `self.ensemble` is a list[State] (see pgprng_common.py) --
    a list of pointers to State objects, with no separate size cached;
    ensemble size is always len(self.ensemble)."""

    def __init__(self, incs, seeds=None):
        self.ensemble = build_ensemble(incs, seeds, "incs", "seeds")
        self.cursor = 0

    def next(self):
        state_obj = self.ensemble[self.cursor]
        out = state_obj.step()
        self.cursor = (self.cursor + 1) % len(self.ensemble)
        return out


class CombinedPRNG:
    """Clock-stepping (fixed round-robin) design. Exact period:
    2**64 * e1 * e2 when e1, e2 are coprime (proven, and confirmed
    exactly in every toy-scale test this session)."""

    def __init__(self, incs_1, incs_2, seeds_1=None, seeds_2=None):
        validate_disjoint(incs_1, incs_2, "incs_1", "incs_2")
        self.ens1 = RotatingEnsemble(incs_1, seeds_1)
        self.ens2 = RotatingEnsemble(incs_2, seeds_2)

    def next(self):
        out1 = self.ens1.next()
        out2 = self.ens2.next()
        return mix64(out1 ^ out2)

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()


class SelectionCombinedPRNG:
    """Mutual PRNG-selection design. Period not exactly calculable;
    order-of-magnitude estimate ~2**128 from the toy-scale M**2
    pattern found this session (see module docstring). ensemble_1/
    ensemble_2 are each a list[State] (see pgprng_common.py) -- a list
    of pointers to State objects, not two parallel lists."""

    def __init__(self, incs_1, incs_2, seeds_1=None, seeds_2=None):
        validate_disjoint(incs_1, incs_2, "incs_1", "incs_2")
        self.ensemble_1 = build_ensemble(incs_1, seeds_1, "incs_1", "seeds_1")
        self.ensemble_2 = build_ensemble(incs_2, seeds_2, "incs_2", "seeds_2")

        self.last_out1 = 0
        for s in self.ensemble_1:
            self.last_out1 ^= s.state
        self.last_out2 = 0
        for s in self.ensemble_2:
            self.last_out2 ^= s.state

        self.select_counts_1 = [0] * len(self.ensemble_1)
        self.select_counts_2 = [0] * len(self.ensemble_2)

    def next(self):
        idx2 = self.last_out1 % len(self.ensemble_2)
        new_out2 = self.ensemble_2[idx2].step()
        self.select_counts_2[idx2] += 1

        idx1 = self.last_out2 % len(self.ensemble_1)
        new_out1 = self.ensemble_1[idx1].step()
        self.select_counts_1[idx1] += 1

        self.last_out1 = new_out1
        self.last_out2 = new_out2

        return mix64(new_out1 ^ new_out2)

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

    clock_gen = CombinedPRNG(INCS_101, INCS_103, seeds_1=list(FIXED_SEEDS_101), seeds_2=list(FIXED_SEEDS_103))
    select_gen = SelectionCombinedPRNG(INCS_101, INCS_103, seeds_1=list(FIXED_SEEDS_101), seeds_2=list(FIXED_SEEDS_103))

    clock_stats = run_battery(clock_gen)
    select_stats = run_battery(select_gen)

    print(f"\n{'='*70}\nCOMPARISON over {N_SAMPLES:,} samples each, e1={E1} e2={E2}, identical increments and seeds\n{'='*70}")
    print(f"{'metric':<18} {'clock-step':>14} {'PRNG-select':>14} {'ideal':>10}")
    ideals = {"popcount_mean": 32.0, "popcount_std": 4.0, "bit0_frac": 0.5,
              "bit1_frac": 0.5, "bit2_frac": 0.5, "lag1_corr": 0.0, "distinct_frac": 1.0}
    for key in clock_stats:
        print(f"{key:<18} {clock_stats[key]:>14.5f} {select_stats[key]:>14.5f} {ideals[key]:>10.4f}")

    counts1 = np.array(select_gen.select_counts_1)
    counts2 = np.array(select_gen.select_counts_2)
    exp1 = N_SAMPLES / len(select_gen.ensemble_1)
    exp2 = N_SAMPLES / len(select_gen.ensemble_2)
    print(f"\nselection fairness (PRNG-selection design only):")
    print(f"  101-side: expected {exp1:.0f} picks/component, "
          f"observed min={counts1.min()} max={counts1.max()} std={counts1.std():.1f}")
    print(f"  103-side: expected {exp2:.0f} picks/component, "
          f"observed min={counts2.min()} max={counts2.max()} std={counts2.std():.1f}")
