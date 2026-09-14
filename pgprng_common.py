# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
Shared building blocks for the pgprng family (pgprng_generator.py,
pgprng_legacy_clockstep.py): the mix64 finalizer and its inverse, the
State object each ensemble component is built from, and the
input-validation/ensemble-construction helpers all three generator
classes (RotatingEnsemble, CombinedPRNG, SelectionCombinedPRNG,
FilteredSelectionPRNG) need in essentially identical form.

Object model: each ensemble component's fixed increment and current
running state used to be tracked in two separate parallel lists (one
list of increments, one list of states, kept in sync by matching
index) under names like incs_11/states_11 -- names that also baked in
a specific ensemble size (11, from this project's original toy scale)
that stopped being true once it grew to 101/103 components. Increment
and state are now bundled into a single State object per component,
and each ensemble (ensemble_1/ensemble_2 in the generator classes
below) is just one list of pointers to State objects -- there is no
second list to keep in sync, and no separate size ever cached; an
ensemble's size is always len(ensemble), read live, since ensemble
sizes are meant to be something the code can change.

Pulled out during a refactor pass; the mix64/validation logic itself
introduces no new behavior (verified by comparing generator output on
a fixed seed before and after). The State-object restructuring is new
as of this pass -- verified the same way, against the prior two-list
representation.
"""
from dataclasses import dataclass
import secrets

__version__ = "1.1"  # single source of truth for the package version;
                      # see CHANGELOG.md for what changed each release.

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


@dataclass
class State:
    """One ensemble component, as a single object rather than two
    parallel-list entries: its fixed `increment` (set once at
    construction, never changes) and its current running `state`
    (mutates every time this component is selected). ensemble_1/
    ensemble_2 in the generator classes below are each just a list of
    these -- a list of pointers to State objects, nothing more
    structural than that, so ensemble size is simply len(ensemble)."""
    increment: int
    state: int

    def step(self) -> int:
        """Advance state by increment (mod 2**64) and return the new
        state -- the one place the additive recurrence itself lives."""
        self.state = (self.state + self.increment) & MASK64
        return self.state


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


def validate_disjoint(incs_a, incs_b, name_a="incs_1", name_b="incs_2"):
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


def build_ensemble(incs, seeds=None, incs_name="incs", seeds_name="seeds"):
    """Validate increments and seeds, then bundle them into a list of
    State objects -- one ensemble, represented as a single list of
    pointers to State objects rather than two parallel lists. This is
    the only place an ensemble gets constructed anywhere in this
    codebase; every generator class builds ensemble_1/ensemble_2 (or,
    for RotatingEnsemble, its own single `ensemble`) by calling this."""
    incs = validate_increments(incs, incs_name)
    seeds = prepare_seeds(seeds, incs, seeds_name)
    return [State(increment=inc, state=seed) for inc, seed in zip(incs, seeds)]
