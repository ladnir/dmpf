# Distributed multi-point functions

Research sources for *Fully Distributed Multi-Point Functions for PCGs and
Beyond*. Waterfall is the generic construction; Reverse Cuckoo is the
Ring-LPN-specific construction with an additional leakage assumption.

## Build the paper

Install a TeX distribution with `latexmk`, then run from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The output is `main.pdf`. Generated PDFs, auxiliary files, and local rendering
directories are not committed. Figures embedded in the TeX sources and
`paper.bib` remain part of the source distribution.

## Run the analysis tests

Use Python 3.12 and an isolated environment:

```sh
python -m venv .venv
# Linux/macOS:
.venv/bin/python -m pip install -r requirements-test.txt
.venv/bin/python -m pytest -q experiments
# Windows: use .venv/Scripts/python.exe instead.
```

The tests check finite distribution identities, cost calculations, and
analysis code. They do not establish cryptographic security. The CSVs in
`experiments/` preserve prior runs; unit tests do not rerun the production
attack sweeps or benchmarks. Their methodology and limitations are recorded
in `experiments/rev_cuckoo_lpn_attack_study.md` and the other experiment notes.

CI runs these tests and builds the full paper from tracked sources. The build
fails on unresolved references, duplicate labels/destinations, or overfull
boxes. Its PDF artifact is a build check, not a substitute for visual review.

## Release scope

The current paper states its ideal-subprotocol and concrete-generator
assumptions separately. The application-level OLE composition argument is
not complete. The implementation filters supports at the analyzed Goldilocks
point; other parameter sets remain unvalidated. Historical measurements do not describe all
current hardening changes. See `RELEASE_STATUS.md` for the release checkpoint.
