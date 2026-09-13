# Thorough statistical testing and publication-worthiness analysis
### Coupled additive-ensemble PRNG, e = 101 / e = 103, clock-stepping vs. PRNG-selection

Author: Aly Graham, with Claude (Anthropic) — analysis and code assistance

---

## 1. What was tested

Both generator designs, exactly as delivered in `combined_prng_101_103.py`, using the actual saved increment sets in `prng_m_101_103.csv` (101 + 103 = 204 distinct odd 64-bit increments, drawn uniformly over the full odd range -- these are additive increments, not multipliers; see the terminology note in `pgprng_common.py`):

- **Clock-stepping** (`CombinedPRNG`): fixed round-robin cursor on each ensemble; exact period `2^64 × 101 × 103 ≈ 2^77.3`.
- **PRNG-selection** (`SelectionCombinedPRNG`): each ensemble's previous output selects which element of the *other* ensemble advances; empirically-estimated period `~2^128` (order of magnitude, from the toy-scale `M²` pattern found earlier).

No PractRand binary or NIST STS reference package could be reached from this sandbox (SourceForge is blocked by organization policy, the guessed GitHub mirror didn't exist, and no matching package exists on PyPI under any name tried). Rather than leave the request unanswered, I implemented a core subset of the NIST SP 800-22 test battery directly in Python/NumPy, validated it against known-good and known-bad control sequences and a brute-force minimal-LFSR search (for the linear complexity test) before trusting it, and ran it at real statistical scale.

**Tests implemented** (each yields a p-value; standard threshold p ≥ 0.01 to "pass" a single run):
Monobit (frequency), block frequency, runs, longest run of ones in a block, discrete Fourier transform (spectral), cumulative sums (forward and reverse), approximate entropy (m=5), serial (m=5), binary matrix rank, and linear complexity (Berlekamp–Massey, block size 500).

**Not implemented** (scope limitation, stated plainly rather than glossed over): non-overlapping and overlapping template matching, Maurer's universal statistical test, and the random-excursions tests. These need either large fixed reference tables/dictionaries or many more samples than was practical here. Their omission means "passed every test in this battery," not "passed every test that exists."

**Methodology** (matches real NIST STS practice, not a single lucky run): each design was run as **500 independent sequences of 1,000,000 bits each**, every sequence freshly and genuinely seeded via `secrets.randbits(64)` — the actual production seeding path, not a fixed demo seed. For each test, two things were checked across the 500 sequences:
- **Proportion passing**: how many of the 500 p-values were ≥ 0.01 (expected ~495; the standard acceptable band down to ~488/500).
- **Uniformity of the p-values themselves**: a healthy generator's p-values should be uniformly distributed on [0,1], not just individually above 0.01 — checked via a chi-square goodness-of-fit over 10 bins.

Total data generated and tested: 1,000,000,000 bits (125,000,000 words) per design, 2 billion bits overall.

## 2. Results

| test | clock pass% | select pass% | clock uniformity p | select uniformity p |
|---|---:|---:|---:|---:|
| monobit | 99.0% | 98.6% | 0.929 | 0.809 |
| block frequency | 98.8% | 99.2% | 0.874 | 0.276 |
| runs | 99.6% | 98.6% | 0.542 | 0.370 |
| longest run | 99.6% | 99.0% | 0.898 | 0.237 |
| DFT / spectral | 97.8% | 98.6% | 0.430 | 0.092 |
| cumulative sums (fwd) | 99.6% | 98.2% | 0.703 | 0.313 |
| cumulative sums (rev) | 99.2% | 98.2% | 0.989 | 0.467 |
| approximate entropy | 98.8% | 99.0% | 0.963 | 0.901 |
| serial (∇) | 99.0% | 99.2% | 0.161 | 0.716 |
| serial (∇²) | 99.4% | 99.8% | 0.341 | 0.579 |
| binary matrix rank | 99.2% | 99.2% | 0.014 | 0.959 |
| linear complexity | 99.2% | 98.6% | 0.880 | 0.820 |

Required pass count for 500 sequences at α = 0.01: **≥ 488.3/500 (97.7%)**. Required uniformity p-value: **≥ 0.0001**.

**Every cell in this table clears both thresholds, for both designs, on all 12 tests.** Neither generator failed anything in this battery.

## 3. Is there a real difference between the two designs?

Honestly: no detectable one, at this data volume. The pass rates for the two designs are within a percentage point or two of each other on every test — differences fully consistent with sampling noise at n=500 (the binomial standard deviation on a 99%-expected proportion at n=500 is about ±0.4 percentage points, so the 0.2–1.4 point gaps seen here are unremarkable). The uniformity p-values bounce around with no consistent pattern favoring either design — clock-stepping's rank-test uniformity (0.014) is the single lowest value in the whole table, but it's still nearly 150× above the 0.0001 failure threshold, and with 24 independent uniformity checks run, seeing one fall as low as 0.014 by chance is itself unremarkable (about a 1-in-3 chance that at least one of 24 uniform p-values lands below 0.05). It is not a flag against clock-stepping, and it does not make PRNG-selection "better" — it's noise.

This matches what you predicted before I ran anything: at the level of statistical testing available here, there's no short (or even fairly long) test that distinguishes them. Both are far larger in period and far better-mixed (via `mix64`) than anything this battery, or most practical applications, would ever be able to tell apart. If a real quality difference between clock-stepping and PRNG-selection exists, it did not show up in two billion tested bits, and finding it — if it's there at all — would likely require either the missing template/universal/excursions tests, a purpose-built battery like PractRand run to gigabytes-per-second scale, or a targeted cryptanalytic attack rather than a general statistical battery. I want to be direct that I have not shown PRNG-selection is the better design, even though that was your hope going in; I've shown the two are statistically indistinguishable within what I was able to test.

## 4. Is this worth a publication?

Short answer: **not as a new PRNG design**, but there is a narrower, honest way to think about what's genuinely yours in this work.

**What this session actually discovered, in rough order of how original it is:**

- **Exact bit-level period theorems for XOR-ensembles of additive generators.** For any odd ensemble size e, bit k of the combined output has period exactly 2^(k+1), independent of which odd increments are used — proven via a GF(2) linear-recurrence argument and verified exhaustively. For even e, bit 0 freezes and higher-bit periods depend on residue-parity conditions mod small powers of 2. This is a clean, small, correct piece of applied algebra. It may already be folklore in the PRNG-design community (additive/Fibonacci-style generators' low-bit weaknesses are well known), but the exact closed-form period-per-bit statement and its clean proof are a legitimate, citable technical result — the kind of thing that could be a short note or an appendix, not a standalone paper.
- **N-invariant closed-form statistics for the additive generator's balance and runs behavior**, verified to match simulation exactly at N=16, 24, and 32, enabling O(1)-per-candidate "best increment" ranking at N=64 without simulation. This is a nice, practical piece of engineering, but it's a straightforward consequence of the additive generator's well-understood combinatorics (n_down = m is a textbook fact), not new mathematics.
- **The exact period formula for clock-stepping**, `2^64 × e1 × e2` for coprime e1, e2, is mathematically identical in structure to L'Ecuyer's classical **combined multiple recursive generators** ("Combined Multiple Recursive Random Number Generators," *Operations Research* 44(5), 1996; "Good Parameters and Implementations for Combined Multiple Recursive Random Number Generators," *Operations Research* 47(1), 1999) and the long-standing **combined linear congruential generator** technique — multiplying periods together via coprime component periods is exactly what that literature already does, decades ago, with careful spectral-test-driven parameter selection this session did not attempt. This piece is a rediscovery, not new.
- **The mutual output-driven "PRNG-selection" coupling** — one generator's output deciding how another one steps — is the core idea behind an entire family of **clock-controlled stream cipher generators** going back to the 1980s: the alternating step generator, the shrinking generator, stop-and-go generators, and the self-shrinking generator, along with a directly-named line of work on **coupled/combined congruential generators (CLCG)** in the RNG literature (e.g., "Pseudorandom Bit Generation Using Coupled Congruential Generators," IEEE, and follow-on VLSI-architecture papers). Applying that clocking idea specifically to modular *additive* ensembles, using the full 64-bit output mod ensemble-size as the selector (rather than a single control bit as classical clock-control does), is a genuine variation in the details, but the architectural idea itself is well-established prior art, not a new concept.
- **The empirical `M²` cycle-length pattern for PRNG-selection**, confirmed exactly across 9 of 10 attempted toy-scale configurations (one timed out unresolved), independent of e1 and e2, is the single most interesting unexplained finding from this whole project. It is *not* proven, only demonstrated at tiny scale (N ≤ 8), and I have not searched the CLCG/clock-controlled-generator literature for whether this exact period behavior is already characterized there (it plausibly is, given how mature that literature is — clock-controlled LFSR period behavior has been analyzed in depth since the 1980s). If, after a real literature search and an attempt at a rigorous proof, this pattern turned out to be both correct and not already known for this exact coupling structure, *that specific result* — not the PRNG as a whole — would be the one piece of this project with a plausible path to a short technical note.
- **The N=64 near-duplicate-increment defect** (drawing ensemble increments from a narrow "individually good" pool causes near-identical component generators and catastrophic bit-bias) is a genuinely useful practical lesson, but it's a special case of a well-known general principle in combined-generator design: components must be selected for independence from each other, not just individual quality — again consistent with, not contradicting, the existing CMRG literature's careful multi-component parameter selection.

**Why I don't think this clears the bar for a peer-reviewed publication as presented:** every architectural building block — combining generators for period multiplication, output-driven mutual clocking, additive congruential components, a nonlinear finalizer — is individually well-established, in some cases for 30–40 years. The statistical testing here, while thorough by the standards of a personal project, doesn't show either design outperforming (or even matching in rigor of testing) well-optimized modern generators like PCG64, xoshiro256**, or counter-based designs like Philox, all of which are faster, have been tested to far greater scale (PractRand at terabytes, TestU01's BigCrush), and have published, peer-reviewed spectral-test analyses backing their parameter choices. A publication needs to add something the field doesn't have; right now this project mostly re-derives, in a self-contained and pedagogically clean way, things the field already has.

**What I'd actually recommend, if you want to take this further:**
1. Do a proper literature search (not just what I could find here) specifically on clock-controlled/coupled congruential generator period theory, to see whether the `M²` pattern is already known.
2. If it isn't, try to prove it rather than extrapolate it — the toy-scale evidence is suggestive but thin (10 configurations, one unresolved).
3. If you want a citable artifact from this work regardless, the exact bit-period theorem for XOR ensembles (section 4, bullet 1) is the cleanest, most self-contained, most clearly-yours result here, and would make a solid short note or blog-style writeup even without the rest.

## 5. Files delivered
- `nist_battery.py` — the statistical test implementations (validated against controls before use)
- `run_thorough_tests.py` — the test driver used to produce the results above
- `nist_battery_results.json` — full per-sequence p-values for every test, both designs (500 × 12 × 2 = 12,000 p-values)
