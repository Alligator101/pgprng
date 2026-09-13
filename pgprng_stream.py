# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
Streams raw binary output from FilteredSelectionPRNG (the current best
candidate: mix64 per ensemble before XOR, selection driven by mixed
values) to stdout as fast as Python allows, for piping into PractRand's
RNG_test_stdin.

Buffers BATCH words at a time in a numpy uint64 array before writing,
rather than writing one value at a time, since per-call write overhead
would otherwise dominate.

Uses the actual delivered, already-validated increment set
(pgprng_increments_101_103.csv -- additive increments, not multipliers;
see the terminology note in pgprng_common.py) and genuinely random
(secrets.randbits) seeding -- the real production configuration, not a
fixed test seed.
"""
import sys
import os
import numpy as np

from pgprng_generator import FilteredSelectionPRNG

INCS_101, INCS_103 = [], []
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pgprng_increments_101_103.csv")) as f:
    next(f)
    for line in f:
        ens, inc = line.strip().split(",")
        (INCS_101 if ens == "101" else INCS_103).append(int(inc))

BATCH = 1_000_000

if __name__ == "__main__":
    gen = FilteredSelectionPRNG(INCS_101, INCS_103)
    out = sys.stdout.buffer
    buf = np.empty(BATCH, dtype=np.uint64)
    try:
        while True:
            for i in range(BATCH):
                buf[i] = gen.next()
            out.write(buf.tobytes())
    except (BrokenPipeError, KeyboardInterrupt):
        pass
