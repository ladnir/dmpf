# Release checkpoint — 2026-09-23

## Saved snapshots

- Paper and analysis snapshot: `ea0ccff` in this repository.
- libOTe implementation: `8bcd440` on `codex/dmpf-package` in
  `osu-crypto/libOTe`. This includes the sparse-DPF correction fixes, activity
  masking, and both AES rekeying schedules.
- libOTe compatibility/release notes: `e46c72d`; only documentation and cache
  exclusions changed after the implementation snapshot.

The implementation snapshot is running through the existing Linux/ASan,
macOS, and Windows CI workflow:
https://github.com/osu-crypto/libOTe/actions/runs/35874069035
Do not infer success from dispatch; check the completed jobs.

## Completed release cleanup

- Committed the previously untracked paper sections, analysis programs,
  recorded experiment results, implementation files, and regression tests.
- Excluded local build products, virtual environments, and Python caches.
- Added pinned Python test dependencies and paper/test CI.
- Added DMPF peer, OT-schedule, cached-state, and serialized-key compatibility
  notes in libOTe. These notes distinguish the generic Waterfall construction
  from the application-specific Reverse-Cuckoo use case.

Local validation: 216 paper-analysis tests pass with Python 3.12 and the pinned
dependencies. Seven libOTe hash-conditioning tests pass. The focused 23-test
WSL/GCC runner passed with defaults and with `-domain 4097` during the release
audit. The existing 72-page PDF has stable references and no overfull boxes;
the current cleanup does not change its TeX content. Remote CI is the remaining
cross-platform check. No new performance benchmarks were run.

## Remaining decisions and work

1. **Ring-LPN support rejection implemented for the analyzed profile.**
   Setup now filters Goldilocks supports at ring degree 2^20, four polynomials,
   and weight 16. It requires at least 61 occupied residues modulo 128,
   counted separately by polynomial, and resamples the whole tuple on
   rejection. Other profiles retain their unfiltered, unvalidated sampler.
   No threshold or security estimate is extrapolated to them.
2. **Complete the OLE composition argument.** The DMPF functionality permits
   the corrupt output share to be selected. Show the target application
   correlation distribution and simulation, rather than inferring them
   from DMPF correctness/privacy alone. Do not strengthen the functionality
   or change the protocol without discussion.
3. **Retain or revise the concrete AES assumption.** The inherited evaluator
   and rekeying schedules are explicitly conjectured instantiations, not
   consequences already proved by the abstract DPF theorem. Current release
   notes preserve that scope; this cleanup does not add a proof requirement
   beyond the paper's stated research claims.
4. **Refresh end-to-end measurements after the sampler is settled.** Existing
   tables describe an earlier revision. Record both constructions, setup,
   repeated expansion, and communication at a fixed commit. Run benchmarks
   sequentially, never concurrently.
5. **Review completed CI and freeze the release.** Fix any code failures,
   distinguish infrastructure failures, and pin the final paper/code commits.
   No release tag or merge into a shared libOTe branch has been made.

Recommended next substantive task: the OLE composition argument. Await
Peter's direction on the remaining gaps.
