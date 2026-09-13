pgprng -- "Pretty Good PRNG": combined/coupled generator experiments

Author: Aly Graham, with Claude (Anthropic) -- analysis, code assistance,
and testing.

This repository is a personal research project exploring
combined/coupled pseudorandom number generator (PRNG) constructions
built from many small additive congruential components, combined via XOR
and a nonlinear finalizer (mix64, a MurmurHash3-style bit-mixing
function). It grew out of curiosity about whether mutually-clocked,
output-coupled generators -- an old idea in the combined/coupled-LCG and
clock-controlled stream cipher literature -- could be assembled from very
simple parts and still hold up under real statistical and cryptanalytic
scrutiny.

Straight answer up front, so it doesn't get lost below: this is not a
claim of a new, publication-worthy PRNG or a proven-secure cipher. Every
individual architectural idea used here -- combining generators for
period multiplication, output-driven mutual clocking, additive
congruential components, a nonlinear finalizer -- is well-established, in
some cases for 30-40 years (see Prior art, below). What's here is a
clean, tested, honestly-documented implementation and empirical study of
that idea, including the mistakes found and fixed along the way.

The two designs

-   CombinedPRNG (pgprng_legacy_clockstep.py) -- fixed round-robin
    ("clock-stepping") combination of two ensembles of additive
    generators, XORed together and passed through mix64. Its period is
    exactly provable (see Period / cycle structure).
-   FilteredSelectionPRNG (pgprng_generator.py) -- the current,
    recommended design. Each ensemble's raw output is passed through
    mix64 before combination (not after), and which element of the other
    ensemble advances next is chosen by the mixed output, not the raw
    value. This closes a real, previously-present flaw in the
    clock-stepping/earlier SelectionCombinedPRNG design: mix64(A ^ B) is
    freely, exactly invertible by anyone who knows the algorithm (mix64
    is a published, efficiently-invertible bijection), so an attacker
    recovers the raw A ^ B value for free, with zero search.
    mix64(A) ^ mix64(B) has no such single-step inversion.

Statistical testing performed

-   A hand-implemented NIST SP 800-22-style battery (monobit, block
    frequency, runs, longest-run, binary matrix rank, DFT/spectral,
    cumulative sums, approximate entropy, serial, linear complexity via
    Berlekamp-Massey) run across 500 x 1,000,000-bit sequences -- see
    pgprng_nist_battery.py, pgprng_run_tests.py, final_report.md for the
    full write-up and results.
-   A real run of PractRand 0.86 (the purpose-built PRNG statistical
    test suite, not a hand-rolled substitute) against
    FilteredSelectionPRNG's actual output stream, using genuine
    secrets-seeded state and the real 101+103-increment production
    configuration: 8 GiB (2^33 bytes) streamed clean, zero anomalies
    across all 12 standard subtests, over roughly 26 minutes of
    continuous testing. This is real evidence of no obvious statistical
    defect at the gigabyte scale -- it is not a certification, and
    PractRand's expanded test modes and multi-terabyte runs (the
    standard for a serious PRNG candidate) have not been done.

Cryptographic analysis (preliminary -- expect further scrutiny)

Kerckhoffs's-principle assumption throughout: the algorithm is fully
public; only the seeds and increments are secret.

-   The flaw already mentioned above (mix64(A^B) being trivially
    invertible in the older/clock-stepping combination order) was found
    and demonstrated empirically this session, and is the direct
    motivation for FilteredSelectionPRNG's mix-before-combine ordering.
-   A hand-built brute-force / algebraic attack
    (pgprng_crack_attempt.py) was run at toy scale only -- 8-14 bit
    words, one component per ensemble side (not the real 101/103) --
    because no SAT/SMT/lattice solver was available in the development
    environment that day (a package-registry outage blocked installing
    z3, sympy, fpylll, gmpy2, etc.; this is disclosed plainly rather
    than glossed over). At that toy scale:
    -   The attack fully recovered the true secret (up to the expected
        A/B-label symmetry) against the new mix-then-combine design in
        every trial, in well under two minutes even at the largest toy
        size tested. This is expected, not alarming: with only one
        unknown component per side, exhaustive search is trivially
        tractable, and mix64 was never meant to resist brute force once
        every candidate is already being tried -- only to remove the
        free, zero-search shortcut.
    -   Against the older combine-then-mix design, the same search never
        recovered the literal original secret, but reliably found a
        different parameter set that reproduces the exact same output
        forever (verified out to 200,000 samples -- more than a full
        period at that toy word size). This is not evidence the old
        design is harder to crack; it's consistent with the
        free-inversion flaw above -- the attacker already has the
        internal driving value for free and only needs any
        operationally-equivalent parameter set to extrapolate forward,
        which this search finds easily.
    -   This does not constitute an attack on, or a security claim
        about, the real 204-component/64-bit-word configuration. The
        real system's search space (roughly 2^128 per side) plus an
        unknown, hidden per-step component-selection index puts it far
        outside the reach of this or any brute-force/meet-in-the-middle
        technique tested here. Nobody should read "toy version cracked"
        as "real generator cracked."
-   No independent, professional cryptanalytic review has been done.
    This code has not been evaluated by anyone besides its author and an
    AI assistant. Anyone attempting a serious attack against the real
    configuration is doing exactly what this project has not yet been
    able to do, and any findings would be genuinely useful and welcome.

Period / cycle structure -- OPEN QUESTION, read before relying on this

This is the single most important unresolved question about this
project, and it is not fully solved.

-   CombinedPRNG (clock-stepping): solved exactly. Its period is
    provably 2^64 x e1 x e2 when e1 and e2 are coprime (standard
    combined/coupled-generator period theory -- see Prior art). With the
    shipped 101/103 configuration that's 2^64 x 101 x 103.
-   FilteredSelectionPRNG (mutual selection): not solved in closed form
    at real scale. A toy-scale cycle-detection study was run this
    session (exact, exhaustive cycle-tracing over every possible seed at
    the smallest sizes; large random samples at larger sizes) with these
    findings:
    -   At every (word size, ensemble sizes) configuration tested, every
        seed tested converged to the same single cycle length -- no
        evidence was found, at this scale, of multiple different periods
        coexisting for one fixed configuration.
    -   That cycle length's dependence on ensemble size is not yet fully
        characterized. It measured 2 x 2^N (small) when both ensembles
        had exactly 2 components, and jumped to (2^N)^2 once either
        ensemble had 3 or more components -- consistently, across every
        configuration where this was checked, independent of further
        increases in ensemble size up to 5 components per side.
    -   This is a toy-scale empirical finding, not a proof, and it has
        not been extended to real 64-bit words with 101/103 components.
        If the same qualitative pattern held at real scale (single
        global period, growing roughly with the square of the
        per-component state space once ensembles exceed a couple of
        elements), it would suggest an astronomically large period -- but
        this is extrapolation, not demonstrated fact. Resolving this
        properly (ideally with an actual proof, not just larger toy
        runs) is the most valuable open technical task left in this
        project, and exactly the kind of thing worth opening a GitHub
        issue or PR about.

Is the code parameterized, or hard-coded to 101/103?

Fully parameterized -- 101 and 103 are just the choice made for this
project's own testing, not a requirement of the code. Looking directly
at FilteredSelectionPRNG.__init__ (and CombinedPRNG's), the ensemble
sizes are simply len(incs_11) and len(incs_13) -- whatever length lists
of increments you pass in. (These are additive increments -- each
component advances as x -> x + increment mod 2^64 -- not multipliers;
there is no multiplicative term anywhere in this recurrence. See the
terminology note in pgprng_common.py.) The only constraints actually
enforced in code are: each increment must be an odd integer in
[1, 2^64-1] (oddness matters -- see below), all increments within and
across both ensembles must be distinct, and each ensemble must have at
least one element. There is no primality check anywhere in either class.
A developer can use any ensemble sizes they like, prime or not.

Why 101/103 (twin primes) were used here anyway: for CombinedPRNG, the
period formula 2^64 x e1 x e2 only gives the full benefit of multiplying
the two ensembles' periods together when e1 and e2 are coprime -- any
twin-prime pair (or simply any two coprime integers) achieves that;
primality itself isn't required, coprimality is what matters, and any
two primes are automatically coprime, which is the convenient reason
primes were chosen. (This "multiplying periods together" is about how
the two ensembles' periods combine -- it has nothing to do with the
per-component recurrence itself, which is purely additive.) For
FilteredSelectionPRNG, coprimality of the ensemble sizes was carried
over as a reasonable default rather than a proven requirement -- given
the period question above is still open for this design, there's no
proof yet that coprimality helps or matters here at all.

Increment values must each be odd: this ensures each individual additive
component (x -> x + increment mod 2^64) has full period 2^64 on its own
(an odd increment is coprime to 2^64, so the sequence visits all 2^64
residues before repeating; an even increment would only reach a subset).

Performance

Efficiency was not ignored, but it was not a priority during development
-- correctness, testability, and honest documentation were. The generator
classes are plain Python: no numpy vectorization, no compiled extension.
Measured throughput (pgprng_stream.py) is roughly 500,000 64-bit words
per second (~4.4 MB/s). That's fine for the statistical and
cryptanalytic testing this project actually did, but it is not
competitive with a production PRNG library, and it will likely be
obvious to anyone reading the code that speed was never the goal. If
real throughput is ever needed -- large-scale Monte Carlo work, for
instance -- vectorizing the inner loop with numpy or rewriting it as a
small C/Cython extension is the natural next step; the underlying
arithmetic (modular addition, mix64) is simple enough that this would be
a contained piece of work, not a redesign.

Known limitations / open questions

-   The FilteredSelectionPRNG period at real (64-bit, 101/103-component)
    scale is unproven -- see above. This is the main open item.
-   PractRand was run to 8 GiB, not the multi-terabyte scale used to
    seriously stress-test production PRNGs.
-   No SAT/SMT/lattice-based cryptanalysis has been attempted against
    the real configuration (only a hand-built brute-force technique at
    toy scale -- see Cryptographic analysis).
-   No independent/professional review of any kind has been done.
-   The exact-period theorem for CombinedPRNG's XOR-ensemble
    construction has not been checked against the CLCG /
    clock-controlled generator literature for prior characterization --
    it's plausible this is already known there.

Prior art

Combined/coupled congruential generators for period multiplication:
L'Ecuyer's work on combined generators (Operations Research, 1996, 1999)
and coupled congruential generator constructions. Output-driven mutual
clocking has a close historical analog in clock-controlled stream cipher
generators (alternating step, shrinking, self-shrinking, and stop-and-go
generators), including known attack techniques against that class (e.g.,
Golic's edit-distance attack) -- cited here as reason for caution, not as
a claim those specific attacks apply to this exact construction.

File manifest

  -----------------------------------------------------------------------
  File                                What it is
  ----------------------------------- -----------------------------------
  pgprng_common.py                    Shared mix64/mix64_inverse and
                                      input-validation helpers used by
                                      both generators (pulled out during
                                      a refactor pass to remove
                                      duplication -- verified
                                      behavior-preserving against
                                      fixed-seed output)

  pgprng_generator.py                 The recommended generator,
                                      FilteredSelectionPRNG
                                      (mix-then-combine,
                                      selection-from-mixed-values)

  pgprng_legacy_clockstep.py          The simpler clock-stepping
                                      generator, CombinedPRNG, plus basic
                                      increment generation

  pgprng_ensemble_construction.py     Rigorous, self-validating increment
                                      construction
                                      (gap/uniformity/smoke-test gated)

  pgprng_increments_101_103.csv       The actual 101+103 increment set
                                      used in all testing here

  pgprng_stream.py                    Streams raw generator output for
                                      piping into external test tools

  pgprng_nist_battery.py,             The hand-built NIST SP 800-22-style
  pgprng_run_tests.py                 statistical battery

  pgprng_crack_attempt.py             The toy-scale cryptanalytic
                                      brute-force attack described above

  RNG_test_stdin.cpp                  Modified PractRand test harness
                                      (stdin adapter) -- see
                                      SESSION_CHANGES_README.txt if
                                      bundled

  final_report.md                     Full statistical results and honest
                                      publication-worthiness assessment

  LICENSE                             MIT License

  pgprng_README.docx,                 Word and plain-text copies of this
  pgprng_README.txt                   same README, for reading locally
                                      without a Markdown renderer
  -----------------------------------------------------------------------

License

MIT License -- see LICENSE. Anyone may use, modify, and redistribute this
code, including commercially, as long as the copyright notice is kept.
