# ---------------------------------------------------------------------
# Developed collaboratively by Aly Graham and Claude (Anthropic).
# Session work, September 2026.
# ---------------------------------------------------------------------
"""
A core subset of the NIST SP 800-22 statistical test suite for random
and pseudorandom number generators, implemented directly in
Python/numpy (no dependency acquisition needed -- PractRand and every
plausible PyPI equivalent were unreachable/nonexistent from this
sandbox; see session notes).

Tests implemented (bit-level, each returns a p-value in [0,1]; the
standard NIST pass/fail threshold is p >= 0.01):

  1. Frequency (Monobit)
  2. Frequency within a Block
  3. Runs
  4. Longest Run of Ones in a Block
  5. Binary Matrix Rank
  6. Discrete Fourier Transform (Spectral)
  7. Cumulative Sums (forward and reverse)
  8. Approximate Entropy (m=5)
  9. Serial (m=5)
  10. Linear Complexity (Berlekamp-Massey, M=500)

Not implemented (out of scope for this pass -- each needs either a
much larger fixed reference table, non-overlapping template
dictionaries, or millions of extra samples for adequate power):
Non-overlapping Template Matching, Overlapping Template Matching,
Maurer's Universal Statistical test, Random Excursions (+ Variant).
These are noted as a scope limitation in the final report, not
silently skipped.

Methodology follows the NIST STS convention for tests 1-4,6-9: run
each test across many independent sequences and report (a) the
*proportion* of sequences with p >= 0.01 (should fall inside the
standard binomial confidence band around the sample count) and (b)
the *uniformity* of the p-value distribution itself, via a chi-square
goodness-of-fit over 10 equal bins (a healthy generator's p-values
should be uniform on [0,1], not just individually above 0.01).
Tests 5 and 10 (rank, linear complexity) are inherently single-long-
sequence, many-sub-block tests and are reported as one chi-square
p-value from one long stream per generator.
"""
import numpy as np
from scipy.special import gammaincc, erfc


# ---------------------------------------------------------------
# bit unpacking helper
# ---------------------------------------------------------------
def words_to_bits(words):
    """words: uint64 numpy array -> uint8 array of bits, MSB first per word."""
    words = np.asarray(words, dtype=np.uint64)
    bits = np.unpackbits(words[:, None].view(np.uint8), axis=1)
    # words are little-endian on this platform; unpackbits on the raw
    # bytes gives bit order within each byte MSB-first but byte order
    # little-endian-first. Reorder bytes big-endian, keep bit order.
    bits = bits.reshape(-1, 8, 8)  # (n, byte_idx_LE, bit_in_byte_MSBfirst)
    bits = bits[:, ::-1, :]        # byte order -> big-endian (MSB word first)
    return bits.reshape(-1)        # flat bit array, MSB of word0 first


# ---------------------------------------------------------------
# 1. Monobit
# ---------------------------------------------------------------
def test_monobit(bits):
    n = len(bits)
    s = np.sum(2 * bits.astype(np.int64) - 1)
    s_obs = abs(s) / np.sqrt(n)
    p = erfc(s_obs / np.sqrt(2))
    return p


# ---------------------------------------------------------------
# 2. Block frequency
# ---------------------------------------------------------------
def test_block_frequency(bits, M=128):
    n = len(bits)
    N = n // M
    if N < 1:
        return np.nan
    blocks = bits[: N * M].reshape(N, M).astype(np.float64)
    pi = blocks.mean(axis=1)
    chi_sq = 4 * M * np.sum((pi - 0.5) ** 2)
    p = gammaincc(N / 2, chi_sq / 2)
    return p


# ---------------------------------------------------------------
# 3. Runs
# ---------------------------------------------------------------
def test_runs(bits):
    n = len(bits)
    pi = bits.mean()
    tau = 2 / np.sqrt(n)
    if abs(pi - 0.5) >= tau:
        return 0.0  # precondition fails -> frequency test itself would fail
    v = 1 + np.sum(bits[1:] != bits[:-1])
    num = abs(v - 2 * n * pi * (1 - pi))
    den = 2 * np.sqrt(2 * n) * pi * (1 - pi)
    p = erfc(num / den)
    return p


# ---------------------------------------------------------------
# 4. Longest run of ones in a block
# ---------------------------------------------------------------
## Each entry: (K = degrees of freedom = len(cuts), cuts, pi) with
## len(pi) == K+1 categories: <=cuts[0], cuts[1], ..., cuts[K-2], >=cuts[K-1]+1
_LRB_TABLE = {
    8: (3, [1, 2, 3], [0.2148, 0.3672, 0.2305, 0.1875]),
    128: (5, [4, 5, 6, 7, 8], [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124]),
    10000: (6, [10, 11, 12, 13, 14, 15],
            [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]),
}


def _longest_run_of_ones(row):
    # row: 1D array of 0/1
    if not row.any():
        return 0
    idx = np.flatnonzero(np.diff(np.concatenate(([0], row, [0]))))
    idx = idx.reshape(-1, 2)
    return int((idx[:, 1] - idx[:, 0]).max())


def test_longest_run(bits):
    n = len(bits)
    if n >= 750000:
        M = 10000
    elif n >= 6272:
        M = 128
    else:
        M = 8
    K, cuts, pi = _LRB_TABLE[M]
    N = n // M
    if N < 1:
        return np.nan
    blocks = bits[: N * M].reshape(N, M)
    longest = np.array([_longest_run_of_ones(blocks[i]) for i in range(N)])
    # K+1 categories: <=cuts[0], cuts[1], ..., cuts[K-2], >=cuts[K-1]+1
    ncat = K + 1
    counts = np.zeros(ncat, dtype=np.int64)
    counts[0] = np.sum(longest <= cuts[0])
    for i in range(1, ncat - 1):
        counts[i] = np.sum(longest == cuts[i])
    counts[-1] = np.sum(longest >= cuts[-1] + 1)
    pi = np.array(pi)
    assert len(pi) == ncat
    chi_sq = np.sum((counts - N * pi) ** 2 / (N * pi))
    p = gammaincc(K / 2, chi_sq / 2)
    return p


# ---------------------------------------------------------------
# 5. Binary matrix rank (GF(2) rank of 32x32 sub-matrices)
# ---------------------------------------------------------------
def _gf2_rank(mat):
    mat = mat.copy().astype(np.uint8)
    rows, cols = mat.shape
    rank = 0
    for col in range(cols):
        pivot = None
        for r in range(rank, rows):
            if mat[r, col]:
                pivot = r
                break
        if pivot is None:
            continue
        mat[[rank, pivot]] = mat[[pivot, rank]]
        for r in range(rows):
            if r != rank and mat[r, col]:
                mat[r, :] ^= mat[rank, :]
        rank += 1
        if rank == rows:
            break
    return rank


def test_binary_matrix_rank(bits, Q=32):
    n = len(bits)
    block_bits = Q * Q
    N = n // block_bits
    if N < 1:
        return np.nan, None
    ranks = np.empty(N, dtype=np.int64)
    for i in range(N):
        chunk = bits[i * block_bits:(i + 1) * block_bits].reshape(Q, Q)
        ranks[i] = _gf2_rank(chunk)
    full = np.sum(ranks == Q)
    minus1 = np.sum(ranks == Q - 1)
    rest = N - full - minus1
    p_full, p_minus1, p_rest = 0.2888, 0.5776, 0.1336
    chi_sq = ((full - p_full * N) ** 2 / (p_full * N)
              + (minus1 - p_minus1 * N) ** 2 / (p_minus1 * N)
              + (rest - p_rest * N) ** 2 / (p_rest * N))
    p = gammaincc(1, chi_sq / 2)
    return p, dict(N=N, full=int(full), minus1=int(minus1), rest=int(rest))


# ---------------------------------------------------------------
# 6. Discrete Fourier Transform (spectral)
# ---------------------------------------------------------------
def test_dft(bits):
    n = len(bits)
    x = 2 * bits.astype(np.float64) - 1
    fft_vals = np.fft.fft(x)
    m = np.abs(fft_vals[: n // 2])
    T = np.sqrt(np.log(1 / 0.05) * n)
    N0 = 0.95 * n / 2
    N1 = np.sum(m < T)
    d = (N1 - N0) / np.sqrt(n * 0.95 * 0.05 / 4)
    p = erfc(abs(d) / np.sqrt(2))
    return p


# ---------------------------------------------------------------
# 7. Cumulative sums (forward + reverse)
# ---------------------------------------------------------------
def _norm_cdf(x):
    from scipy.special import ndtr
    return ndtr(x)


def test_cusum(bits, mode="forward"):
    n = len(bits)
    x = 2 * bits.astype(np.float64) - 1
    if mode == "reverse":
        x = x[::-1]
    s = np.cumsum(x)
    z = np.max(np.abs(s))
    zi = z / np.sqrt(n)

    def term(k_start, k_end, sign_num):
        ks = np.arange(k_start, k_end + 1)
        return ks

    total = 0.0
    k_lo = int(np.floor((-n / z + 1) / 4)) if z > 0 else 0
    k_hi = int(np.floor((n / z - 1) / 4)) if z > 0 else 0
    for k in range(k_lo, k_hi + 1):
        total += _norm_cdf(((4 * k + 1) * z) / np.sqrt(n)) - _norm_cdf(((4 * k - 1) * z) / np.sqrt(n))
    k_lo2 = int(np.floor((-n / z - 3) / 4)) if z > 0 else 0
    k_hi2 = int(np.floor((n / z - 1) / 4)) if z > 0 else 0
    total2 = 0.0
    for k in range(k_lo2, k_hi2 + 1):
        total2 += _norm_cdf(((4 * k + 3) * z) / np.sqrt(n)) - _norm_cdf(((4 * k + 1) * z) / np.sqrt(n))
    p = 1.0 - total + total2
    return float(np.clip(p, 0.0, 1.0))


# ---------------------------------------------------------------
# 8/9. Approximate entropy + Serial (share the overlapping-pattern machinery)
# ---------------------------------------------------------------
def _psi_sq(bits, m):
    n = len(bits)
    if m == 0:
        return 0.0
    padded = np.concatenate([bits, bits[: m - 1]]) if m > 1 else bits
    counts = np.zeros(2 ** m, dtype=np.int64)
    # build overlapping m-bit patterns as integers
    idx = np.zeros(n, dtype=np.int64)
    for j in range(m):
        idx = (idx << 1) | padded[j: j + n].astype(np.int64)
    np.add.at(counts, idx, 1)
    # standard formula: psi_m^2 = (2^m/n) * sum(counts_i^2) - n
    psi = (2 ** m / n) * np.sum(counts.astype(np.float64) ** 2) - n
    return psi


def test_serial(bits, m=5):
    n = len(bits)
    psi_m = _psi_sq(bits, m)
    psi_m1 = _psi_sq(bits, m - 1)
    psi_m2 = _psi_sq(bits, m - 2) if m >= 2 else 0.0
    d1 = psi_m - psi_m1
    d2 = psi_m - 2 * psi_m1 + psi_m2
    p1 = gammaincc(2 ** (m - 2), d1 / 2) if m >= 2 else np.nan
    p2 = gammaincc(2 ** (m - 3), d2 / 2) if m >= 3 else np.nan
    return p1, p2


def test_approx_entropy(bits, m=5):
    n = len(bits)

    def phi(mm):
        if mm == 0:
            return 0.0
        padded = np.concatenate([bits, bits[: mm - 1]]) if mm > 1 else bits
        counts = np.zeros(2 ** mm, dtype=np.int64)
        idx = np.zeros(n, dtype=np.int64)
        for j in range(mm):
            idx = (idx << 1) | padded[j: j + n].astype(np.int64)
        np.add.at(counts, idx, 1)
        c = counts.astype(np.float64) / n
        c = c[c > 0]
        return np.sum(c * np.log(c))

    phi_m = phi(m)
    phi_m1 = phi(m + 1)
    apen = phi_m - phi_m1
    chi_sq = 2 * n * (np.log(2) - apen)
    p = gammaincc(2 ** (m - 1), chi_sq / 2)
    return p


# ---------------------------------------------------------------
# 10. Linear complexity (Berlekamp-Massey over GF(2))
# ---------------------------------------------------------------
def _berlekamp_massey(bits):
    """Standard Berlekamp-Massey over GF(2), implemented with Python
    big integers as bit-packed polynomials (word-level XOR/shift/
    popcount via CPython's arbitrary-precision int, int.bit_count())
    instead of a per-bit inner loop -- an order of magnitude faster
    than the naive array version at M~500, needed to cover enough
    blocks of the long test streams in reasonable time. Verified
    against a brute-force minimal-LFSR search over thousands of short
    random sequences (zero mismatches) before and after this
    bit-packed rewrite."""
    n = len(bits)
    s_list = bits.tolist() if hasattr(bits, "tolist") else list(bits)
    c = 1  # bit i = coefficient c_i; c_0 = 1 fixed
    b = 1
    L = 0
    m = -1
    H = 0  # bit k (after processing index N-1) = s[N-1-k]; history strictly before current index
    for N in range(n):
        sN = s_list[N]
        parity = ((c >> 1) & H).bit_count() & 1
        d = sN ^ parity
        if d:
            t = c
            shift = N - m
            c ^= (b << shift)
            if 2 * L <= N:
                L = N + 1 - L
                m = N
                b = t
        H = (H << 1) | sN
    return L


def test_linear_complexity(bits, M=500):
    n = len(bits)
    N = n // M
    if N < 200:
        return np.nan, None
    Ls = np.empty(N, dtype=np.int64)
    for i in range(N):
        Ls[i] = _berlekamp_massey(bits[i * M:(i + 1) * M])
    mu = M / 2 + (9 + (-1) ** (M + 1)) / 36 - (3 + (-1) ** (M + 1)) / (2 ** M * 1.0)
    sign = 1 if M % 2 == 0 else -1  # (-1)**M
    T = sign * (Ls.astype(np.float64) - mu) + 2 / 9
    cuts = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
    counts = np.zeros(7, dtype=np.int64)
    counts[0] = np.sum(T <= cuts[0])
    for i in range(1, 6):
        counts[i] = np.sum((T > cuts[i - 1]) & (T <= cuts[i]))
    counts[6] = np.sum(T > cuts[-1])
    pi = np.array([0.010417, 0.031250, 0.125000, 0.500000, 0.250000, 0.062500, 0.020833])
    chi_sq = np.sum((counts - N * pi) ** 2 / (N * pi))
    p = gammaincc(3.0, chi_sq / 2)
    return p, dict(N=N, counts=counts.tolist())


# ---------------------------------------------------------------
# driver: run the "cheap" battery on one bit sequence
# ---------------------------------------------------------------
def run_cheap_battery(bits):
    out = {}
    out["monobit"] = float(test_monobit(bits))
    out["block_freq"] = float(test_block_frequency(bits, M=10000 if len(bits) >= 1000000 else 128))
    out["runs"] = float(test_runs(bits))
    out["longest_run"] = float(test_longest_run(bits))
    out["dft"] = float(test_dft(bits))
    out["cusum_fwd"] = float(test_cusum(bits, "forward"))
    out["cusum_rev"] = float(test_cusum(bits, "reverse"))
    out["approx_entropy"] = float(test_approx_entropy(bits, m=5))
    p1, p2 = test_serial(bits, m=5)
    out["serial_1"] = float(p1)
    out["serial_2"] = float(p2)
    return out
