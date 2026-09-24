# Release checkpoint — 2026-09-23

This is the current package checklist. It supersedes earlier release task lists
in the handoff and audit notes; their technical evidence remains relevant.

## Current snapshots

- Paper: presentation cleanup `57b3be1` on `codex/release-readiness`;
  subsequent checkpoint edits change only this checklist. The PDF is 71 pages.
- libOTe: `334ece9` on `codex/dmpf-package`, pushed to `osu-crypto/libOTe`.
  The implementation worktree is clean at this checkpoint.
- [Current cross-platform CI](https://github.com/osu-crypto/libOTe/actions/runs/35950961030)
  targets `334ece9a828055e96490fe422b47cfb3ac91111c` and is running.
  Dispatch is not a passing result. The preceding run on `7f2b4e5` failed in
  Ubuntu's KOS-Dot test and skipped consumer/install checks. Its macOS and
  fallback jobs passed; Windows was still running at this check.

## Completed fixes and decisions

- [x] **Sparse-DPF activity and correction leakage:** activity masking and
  selected-correction low-bit repair are implemented (`8bcd440`), with
  regression coverage for inactive levels, encoding, and both DMPF backends.
- [x] **AES rekeying:** internal-node and cached-leaf schedules rekey after
  1,024 parent positions or leaf positions, respectively (`8bcd440`). This is
  hardening, not a proof of the concrete evaluator.
- [x] **Fresh Ring-LPN public masks:** resampled for every expansion (`b2e9b7b`).
- [x] **Ring-LPN support filter:** implemented for Goldilocks, ring degree
  2^20, four polynomials, and weight 16 (`c2fa1bc`). Each party resamples its
  whole support tuple until it occupies at least 61 residues modulo 128,
  counted separately by polynomial. Both DMPF backends and OLE/triple modes
  use the filter; accepted supports persist across reuse.
- [x] **OLE composition note:** the conditional simulator argument is recorded
  in `OLE_COMPOSITION_NOTE.md`, outside the PDF as requested. It closes in the
  stated ideal-subprotocol model, not unconditionally for concrete AES.
- [x] **Exact Reverse-Cuckoo permutations:** serial privately controlled
  Waksman passes use independently sampled uniform permutations (`68c1a6b`).
- [x] **Independent Reverse-Cuckoo evaluator seeds:** fresh contributions for
  every `(set, partition)` replace batch-shared evaluators (`68c1a6b`).
- [x] **Zero coefficients:** retain uniform field sampling, including zero.
  This approved choice is reflected in the paper; no rejection fix is pending.
- [x] **PPRF output-size mismatch:** expose only the unpadded PPRF prefix during
  expansion, then restore the encoder buffer size (`597a9cc`). Silent OT/VOLE
  tests now cover the 1,024-output regression.
- [x] **Circuit-test file artifacts:** cryptoTools `d870e28`, pinned by libOTe
  `7f2b4e5`, replaces binary/JSON test files with memory streams. Three binary
  round trips passed from an empty, read-only directory without creating files.
  JSON support was disabled locally. The old generated `.bin` was verified
  against the test circuit and removed from the paper working directory.
- [x] **GCC ASan coroutine stack overflow mitigation:** disable only ASan stack
  instrumentation for GCC 13–16 ASan builds (`597a9cc`), including downstream
  C++ template consumers. Heap/global checks remain enabled; Clang retains
  full ASan. Full cross-platform validation remains pending below.
- [x] **KOS-Dot scratch alignment:** `334ece9` uses 32-byte-aligned arrays in
  the non-inlined, non-coroutine transpose helper. The registered check test
  now covers scalar transpose agreement, empty/partial/full chunks, and
  oversized-input rejection. No allocation or protocol arithmetic changed.
- [x] **Package hygiene:** source, analysis, and regression files are tracked;
  build products are excluded; Python dependencies are pinned. Upgrade notes
  document peer, OT-schedule, cached-state, and serialized-key boundaries.
- [x] **Paper presentation pass:** no draft markers in included TeX, stable
  references, no compilation warnings or overfull boxes, and rendered-page
  inspection completed. Numerical tables and plots are unchanged. These edits
  are included in the separate paper checkpoint commit.

## Validation already obtained

- Full Clang 18 ASan suite with the buffer fix: **333 passed, 3 skipped**, with
  no sanitizer or leak report. This precedes the GCC-only workaround, which
  does not change Clang compilation.
- GCC 13.3 and 15.2 standalone coroutine probes and exported-target consumers
  pass with the workaround. Intentional heap overflow remains detected.
  Compiler/version gating and actual libOTe compile/export options were
  checked. This is not a completed full GCC-ASan suite.
- All five focused KOS-Dot tests pass in each of three builds: plain GCC 13,
  GCC 13 ASan with stack checks disabled, and Clang 18 ASan. The test and
  sender/receiver translation units use those configurations. This runner
  links existing non-ASan dependency archives, not a fully instrumented
  dependency build. Full validation remains assigned to CI.
- Exact-permutation/per-list-seed changes passed the 24-suite focused WSL/GCC
  runner at defaults and domain 4097. Four Ring-LPN suites passed with two
  trials; seven exact Python hash-conditioning checks passed.
- Earlier paper-analysis validation passed 216 tests with pinned dependencies.
  No new end-to-end benchmarks were run for this release pass.

## Remaining release work

1. [ ] **Review current CI:** require Ubuntu/GCC-ASan unit tests, source-tree
   and installed-consumer checks, macOS, Windows, and the fallback job. Fix
   code failures; report infrastructure failures separately. The old failing
   run did not reach Ubuntu's consumer/install checks.
2. [ ] **Refresh implementation validation notes:** after CI completes, record
   its result and tested commit in libOTe's `DMPF_RELEASE_NOTES.md`. Do not mark
   it passed from local probes alone.
3. [x] **Commit paper presentation cleanup separately:** saved in `57b3be1`.
   Unrelated local/generated files are excluded;
   historical-performance qualifications are retained.
4. [ ] **Freeze package references:** record final paper/code commits and
   present remaining failures or decisions to Peter. No release tag, upload,
   or merge into a shared libOTe branch is authorized by this checklist.

The KOS-Dot defect was introduced in `d68d3c7` (2026-08-23): its arrays
guaranteed only 16-byte alignment, whereas AVX transpose requires 32.
A focused probe reproduced it without coroutines. Full ASan happened to
align the arrays; disabling stack instrumentation exposed the existing bug.
The production repair is applied and locally tested, with full CI pending.

The earlier `AnyField_F2Ole_Test` crash has a different cause: GCC's ASan stack
epilogue prevents bounded-stack coroutine transfer. A standard-C++ reproducer
with no libOTe/macoro/coproto dependencies again overflowed under GCC 13.3 and
15.2 ASan, but passed 100,000 awaits without ASan and under Clang 18 ASan.
Disassembly shows indirect calls followed by sanitizer cleanup; disabling
ASan stack instrumentation restores tail jumps. The mitigated CI run passed
`AnyField_F2Ole_Test` before encountering the separate KOS-Dot alignment bug.
Diagnostic probes and logs are retained under libOTe's ignored `out/` directory,
including `crash-root-cause/` and `std_task_*`.

## Accepted limitations and deferred work

- **Other Ring-LPN profiles remain unvalidated.** Do not extrapolate the
  support threshold or security estimates beyond the analyzed profile.
- **Concrete AES security remains conjectural.** The inherited evaluator and
  rekeying schedules are not proved consequences of the abstract DPF theorem.
  Peter accepted documenting this boundary rather than requiring a new proof
  for this research release.
- **Reverse Cuckoo remains application-specific and leaky.** Waterfall is the
  generic construction; the specialized path requires the stated assumption.
- **Performance tables remain historical.** Peter accepted existing numbers
  with explicit qualifications. Refreshing setup, repeated-expansion, and
  communication measurements is deferred, not a current release blocker.
  Future benchmarks must run sequentially, never concurrently.
- **GCC workaround follow-up:** retest GCC 17 before extending the range.
  Clang full-ASan coverage remains important while GCC stack checks are disabled.

Recommended next task: review current CI, especially the previously skipped
consumer/install checks, then finalize validation notes and the paper commit.
