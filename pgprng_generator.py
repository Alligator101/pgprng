# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
FilteredSelectionPRNG -- the "filter before combine" redesign of the
PRNG-selection generator, incorporating both changes just requested:

  1. mix64 is applied to EACH ensemble's raw output separately, and the
     two MIXED values are XORed to form the final output -- not (as in
     the original SelectionCombinedPRNG, kept unchanged in
     pgprng_legacy_clockstep.py for future reference/testing) XORing the
     raw values first and mixing once.

     This matters for a reason demonstrated directly in this session:
     mix64(A ^ B) is trivially INVERTIBLE by anyone who knows the
     algorithm -- mix64 (a MurmurHash3-style finalizer) is a published
     bijection with an efficient closed-form inverse (undo each
     xorshift, which is self-inverse at shift=33; undo each multiply
     via the modular inverse of the constant mod 2**64). A
     Kerckhoffs-aware adversary just inverts it and recovers A^B
     exactly, for free -- confirmed empirically against the actual
     delivered CombinedPRNG output stream. mix64(A) ^ mix64(B) has no
     such single-step inversion: there is no one value that was ever
     someone's mix64 output for the adversary to invert their way back
     through.

  2. The mutual-selection indices (which element of the OTHER ensemble
     advances next) are now computed from the MIXED values, not the
     raw ones -- so the index sequence an attacker would need to
     predict is itself hidden behind mix64, not exposed as a raw
     intermediate quantity the way it was in SelectionCombinedPRNG
     (there, last_out11/last_out13 were raw, pre-mix64 values).

Not a claim of proven cryptographic security -- see the caveats
discussed with this change: the two ensembles still evolve via a known,
simple additive recurrence, and this has not been subjected to a real
cryptanalytic attempt. It removes the specific trivial break that
applies to the original design; it does not constitute a security
proof.

The increment sets are the SAME already-validated 101+103 values in
pgprng_increments_101_103.csv (already passed the structural
gap/uniformity gates and the empirical smoke test in
pgprng_ensemble_construction.py) -- reused here rather than redrawn,
since nothing about this change affects increment selection. (These are
additive increments, not multipliers -- see the note in pgprng_common.py.
There is no multiplicative term anywhere in either ensemble's recurrence;
the only real multiplications in this file are inside mix64 itself.)

Object model: ensemble_1 and ensemble_2 are each a list of State
objects (pgprng_common.py) -- one per component, bundling that
component's fixed increment with its current running state. There is
no separate size attribute cached anywhere; ensemble size is always
len(ensemble_1)/len(ensemble_2), read live, since ensemble sizes are
meant to be something the code can change, not something baked into
attribute names the way the old incs_11/incs_13 naming implied a fixed
e=11/e=13 scale that had already stopped being true once this project
grew to 101/103 components.

mix64/mix64_inverse and the input-validation/ensemble-construction
helpers used below live in pgprng_common.py (shared with
pgprng_legacy_clockstep.py) -- they were duplicated inline in both
files until a refactor pass consolidated them; that move changes no
behavior (verified by comparing generator output on a fixed seed
before and after).
"""
from pgprng_common import mix64, mix64_inverse, validate_disjoint, build_ensemble


class FilteredSelectionPRNG:
    """Mutual PRNG-selection with per-ensemble mix64 applied BEFORE
    combination, and selection indices driven by the mixed values.
    ensemble_1/ensemble_2 are each a list[State] (see pgprng_common.py) --
    a list of pointers to State objects, not two parallel lists."""

    def __init__(self, incs_1, incs_2, seeds_1=None, seeds_2=None):
        validate_disjoint(incs_1, incs_2, "incs_1", "incs_2")
        self.ensemble_1 = build_ensemble(incs_1, seeds_1, "incs_1", "seeds_1")
        self.ensemble_2 = build_ensemble(incs_2, seeds_2, "incs_2", "seeds_2")

        # Bootstrap using mixed (not raw) seed values, for consistency
        # with "filter each value before combining it with anything."
        self.last_mixed1 = 0
        for s in self.ensemble_1:
            self.last_mixed1 ^= mix64(s.state)
        self.last_mixed2 = 0
        for s in self.ensemble_2:
            self.last_mixed2 ^= mix64(s.state)

        self.select_counts_1 = [0] * len(self.ensemble_1)
        self.select_counts_2 = [0] * len(self.ensemble_2)

    def next(self):
        idx2 = self.last_mixed1 % len(self.ensemble_2)
        new_mixed2 = mix64(self.ensemble_2[idx2].step())
        self.select_counts_2[idx2] += 1

        idx1 = self.last_mixed2 % len(self.ensemble_1)
        new_mixed1 = mix64(self.ensemble_1[idx1].step())
        self.select_counts_1[idx1] += 1

        self.last_mixed1 = new_mixed1
        self.last_mixed2 = new_mixed2

        return new_mixed1 ^ new_mixed2

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()


if __name__ == "__main__":
    import os
    _here = os.path.dirname(os.path.abspath(__file__))
    INCS_101, INCS_103 = [], []
    with open(os.path.join(_here, "pgprng_increments_101_103.csv")) as f:
        next(f)
        for line in f:
            ens, inc = line.strip().split(",")
            (INCS_101 if ens == "101" else INCS_103).append(int(inc))

    print("=" * 70)
    print("Sanity check 1: mix64_inverse really is the exact inverse of mix64")
    print("=" * 70)
    import random
    rng = random.Random(1)
    ok = all(mix64_inverse(mix64(x)) == x for x in (rng.randrange(0, 1 << 64) for _ in range(100_000)))
    print("mix64_inverse(mix64(x)) == x for 100,000 random x:", ok)

    print("\n" + "=" * 70)
    print("Sanity check 2: inverting the OUTPUT of the new design gives NOTHING")
    print("meaningful (unlike the old design, where inverting output gave the")
    print("exact raw out11^out13 value)")
    print("=" * 70)
    gen = FilteredSelectionPRNG(INCS_101, INCS_103)
    for _ in range(3):
        out = gen.next()
        garbage = mix64_inverse(out)
        # compare against the actual internal mixed values at this step --
        # garbage should match NEITHER of them
        print(f"  output=0x{out:016X}  mix64_inverse(output)=0x{garbage:016X}  "
              f"(this is not last_mixed1=0x{gen.last_mixed1:016X} "
              f"nor last_mixed2=0x{gen.last_mixed2:016X})")

    print("\n" + "=" * 70)
    print("Quick statistical smoke test (200,000 words)")
    print("=" * 70)
    import numpy as np
    gen2 = FilteredSelectionPRNG(INCS_101, INCS_103)
    n = 200_000
    vals = np.empty(n, dtype=np.uint64)
    for i in range(n):
        vals[i] = gen2.next()
    popcounts = np.zeros(n, dtype=np.int64)
    for k in range(64):
        popcounts += ((vals >> np.uint64(k)) & np.uint64(1)).astype(np.int64)
    print(f"mean_popcount={popcounts.mean():.5f} (ideal 32.0)  "
          f"std_popcount={popcounts.std():.5f} (ideal ~4.0)  "
          f"bit0_frac={(vals & 1).mean():.5f} (ideal 0.5)  "
          f"distinct_frac={len(np.unique(vals))/n:.6f} (ideal 1.0)")

    counts1 = np.array(gen2.select_counts_1)
    counts2 = np.array(gen2.select_counts_2)
    print(f"selection fairness: 101-side expected {n/len(gen2.ensemble_1):.0f} min={counts1.min()} max={counts1.max()} std={counts1.std():.1f}")
    print(f"selection fairness: 103-side expected {n/len(gen2.ensemble_2):.0f} min={counts2.min()} max={counts2.max()} std={counts2.std():.1f}")
