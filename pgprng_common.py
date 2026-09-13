# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
Shared building blocks for the pgprng family (pgprng_generator.py,
pgprng_legacy_clockstep.py): the mix64 finalizer and its inverse, and
the input-validation helpers that all three generator classes
(RotatingEnsemble, CombinedPRNG, SelectionCombinedPRNG,
FilteredSelectionPRNG) need in essentially identical form.

Pulled out during a refactor pass whose only goal was removing
duplication -- this file introduces no new behavior. Every function
here previously existed, near-verbatim, inline in two or three places;
moving it here does not change what any generator computes (verified
by comparing generator output on a fixed seed before and after this
refactor).
"""
import secrets

MASK64 = (1 << 64) - 1

C1 = 0xFF51AFD7ED558CCD
C2 = 0xC4CEB9FE1A85EC53
INV_C1 = pow(C1, -1, 1 << 64)
INV_C2 = pow(C2, -1, 1 << 64)


def mix64(value: int) -> int:
    """MurmurHash3-style fmix64 finalizer. A bijection on 64-bit words."""
    z = value & MASK64
    z ^= z >> 33
    z = (z * C1) & MASK64
    z ^= z >> 33
    z = (z * C2) & MASK64
    z ^= z >> 33
    return z


def mix64_inverse(value: int) -> int:
    """Exact inverse of mix64 -- included so the fact that mix64 is
    freely invertible by anyone who knows the algorithm isn't hidden.
    Confirmed exact over 100,000 random values (see pgprng_generator.py's
    __main__ block)."""
    z = value & MASK64
    z ^= z >> 33
    z = (z * INV_C2) & MASK64
    z ^= z >> 33
    z = (z * INV_C1) & MASK64
    z ^= z >> 33
    return z


def validate_increments(incs, name="incs"):
    """Validate a list of per-component increments (the additive step
    each component advances by each time it's selected -- NOT
    multipliers; there is no multiplicative term anywhere in the
    recurrence, only x -> x + increment mod 2**64): at least one, all
    distinct, each an odd integer in [1, 2**64-1]. Returns list(incs).
    `name` is used only to make the error message identify which
    argument was bad."""
    incs = list(incs)
    if len(incs) < 1:
        raise ValueError(f"{name} must contain at least one increment")
    if len(set(incs)) != len(incs):
        raise ValueError(f"all values in {name} must be distinct")
    for inc in incs:
        if not (1 <= inc <= MASK64) or inc % 2 == 0:
            raise ValueError(f"each value in {name} must be an odd integer in [1, 2**64-1], got {inc!r}")
    return incs


def validate_disjoint(incs_a, incs_b, name_a="incs_11", name_b="incs_13"):
    """Raise if the two increment lists share any value."""
    if set(incs_a) & set(incs_b):
        raise ValueError(f"{name_a} and {name_b} must be disjoint")


def prepare_seeds(seeds, incs, name="seeds"):
    """Return a validated seed list matching len(incs): draws fresh
    secrets.randbits(64) seeds if seeds is None, otherwise checks the
    given seeds are the right length and each in [0, 2**64-1]."""
    if seeds is None:
        return [secrets.randbits(64) for _ in incs]
    seeds = list(seeds)
    if len(seeds) != len(incs):
        raise ValueError(f"{name} must have length {len(incs)}, got {len(seeds)}")
    for s in seeds:
        if not (0 <= s <= MASK64):
            raise ValueError(f"each value in {name} must be in [0, 2**64-1], got {s!r}")
    return seeds
