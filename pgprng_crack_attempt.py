# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
Attempt to crack toy-scale versions of both PRNG designs from output
alone (Kerckhoffs: algorithm fully known, only seeds/increments
secret -- ma/mb below are additive increments, not multipliers; there
is no multiplicative term in either toy recurrence, only A -> A + ma
mod 2**N), using a real attack strategy implemented from scratch --
no SAT/SMT solver was available today (PyPI itself returned 403 on
every package tried, including plain "requests"; z3-solver, pysat,
sympy, fpylll, gmpy2 are all absent and none could be installed this
session -- this is disclosed here, not glossed over).

SCOPE, stated plainly: this tests the ONE-COMPONENT-PER-ENSEMBLE case
(e1=e2=1), which isolates exactly the question this session has been
arguing about -- does "mix-then-combine" resist a real attack better
than "combine-then-mix" -- but is NOT a full attack on the real
204-component design (that also depends on not knowing which of 204
components is active each step, which this simplification removes by
having only one component to begin with). The real system's actual
protection against the attack below comes from a different source:
204 unknowns per side instead of 1, discussed at the end.

Attack strategy: brute-force ALL (A_0, m_a) pairs for side A
(2**(2N) candidates). For each candidate, side A's entire trajectory
is fully determined, so side B's trajectory can be ALGEBRAICALLY
DERIVED from the observed output (via XOR and, for the mixed design,
one mix_N inversion) with no additional guessing -- then checked for
consistency (does the derived B sequence really follow B_0, B_0+m_b,
B_0+2*m_b, ...?). This is a real meet-in-the-middle-style technique,
just implemented directly instead of via a generic solver.
"""
import time


def mix_N(value, N, C, shift):
    mask = (1 << N) - 1
    z = value & mask
    z ^= z >> shift
    z = (z * C) & mask
    z ^= z >> shift
    return z


def mix_N_inverse(value, N, C, shift):
    mask = (1 << N) - 1
    inv_C = pow(C, -1, 1 << N)
    z = value & mask
    z ^= z >> shift  # self-inverse when shift > N/2
    z = (z * inv_C) & mask
    z ^= z >> shift
    return z


def make_mix_params(N):
    # shift chosen > N/2 so the xorshift step is self-inverse in one
    # application (same property mix64 relies on at shift=33, N=64)
    shift = N // 2 + 1
    # any odd constant works; reuse mix64's low bits, truncated to N,
    # forced odd
    C = (0xFF51AFD7ED558CCD & ((1 << N) - 1)) | 1
    return C, shift


def gen_old_design(N, A0, ma, B0, mb, n_samples, C, shift):
    """combine-then-mix: y_n = mix(A_n ^ B_n)"""
    mask = (1 << N) - 1
    A, B = A0, B0
    out = []
    for _ in range(n_samples):
        A = (A + ma) & mask
        B = (B + mb) & mask
        out.append(mix_N(A ^ B, N, C, shift))
    return out


def gen_new_design(N, A0, ma, B0, mb, n_samples, C, shift):
    """mix-then-combine: y_n = mix(A_n) ^ mix(B_n)"""
    mask = (1 << N) - 1
    A, B = A0, B0
    out = []
    for _ in range(n_samples):
        A = (A + ma) & mask
        B = (B + mb) & mask
        out.append(mix_N(A, N, C, shift) ^ mix_N(B, N, C, shift))
    return out


def crack(design, N, observed, C, shift):
    """Brute force (A0, ma); derive B's trajectory algebraically;
    check consistency. Returns (A0, ma, B0, mb) or None."""
    mask = (1 << N) - 1
    n = len(observed)
    for A0 in range(1 << N):
        for ma in range(1 << N):
            A = A0
            traj_A = []
            for _ in range(n):
                A = (A + ma) & mask
                traj_A.append(A)
            if design == "old":
                # raw_n = A_n ^ B_n = mix_inverse(y_n)  -> B_n = raw_n ^ A_n
                B_seq = [mix_N_inverse(y, N, C, shift) ^ a for y, a in zip(observed, traj_A)]
            else:
                # mix(B_n) = y_n ^ mix(A_n) -> B_n = mix_inverse(that)
                B_seq = [mix_N_inverse(y ^ mix_N(a, N, C, shift), N, C, shift)
                         for y, a in zip(observed, traj_A)]
            # check B_seq is a consistent additive sequence
            mb_candidate = (B_seq[1] - B_seq[0]) & mask if n >= 2 else None
            if mb_candidate is None:
                continue
            ok = True
            for i in range(1, n):
                if (B_seq[i - 1] + mb_candidate) & mask != B_seq[i]:
                    ok = False
                    break
            if ok:
                B0 = (B_seq[0] - mb_candidate) & mask
                return (A0, ma, B0, mb_candidate)
    return None


if __name__ == "__main__":
    import random
    rng = random.Random(2026)

    print("=" * 70)
    print("Cracking attempt: brute-force-one-side-derive-the-other attack")
    print("(no SAT/SMT solver available this session -- see module docstring)")
    print("=" * 70)

    for N in (8, 10, 12, 14):
        C, shift = make_mix_params(N)
        mask = (1 << N) - 1
        A0 = rng.randrange(0, 1 << N) | 1
        ma = rng.randrange(1, 1 << N) | 1
        B0 = rng.randrange(0, 1 << N) | 1
        mb = rng.randrange(1, 1 << N) | 1
        n_samples = 6

        for design, gen in (("old", gen_old_design), ("new", gen_new_design)):
            observed = gen(N, A0, ma, B0, mb, n_samples, C, shift)
            t0 = time.time()
            result = crack(design, N, observed, C, shift)
            elapsed = time.time() - t0
            true_vals = (A0, ma, B0, mb)
            got_it = (result == true_vals) or (result == (B0, mb, A0, ma))  # symmetric under A/B swap
            print(f"N={N:2d} design={design:<3s}  search_space=2^{2*N}={1<<(2*N):>12,}  "
                  f"time={elapsed:7.3f}s  recovered_exact_secret={got_it}")
