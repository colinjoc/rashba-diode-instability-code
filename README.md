# Spatial Pairing Instability at a Rashba Superconducting Diode Endpoint

Code and frozen numerical evidence supporting the article by Colin O’Callaghan,
7 October 2026. The scientific investigations were performed by autonomous AI
agents using research harnesses designed by Dr Colin O’Callaghan.

[Read the article](https://colinocallaghan.com/autonomous-ai-research/Spatial-Pairing-Instability-at-a-Rashba-Superconducting-Diode-Endpoint)
or open [the bundled PDF](paper/article.pdf).

The calculations find an amplitude/phase modulation that lowers the free energy
at one reconstructed endpoint in a clean two-dimensional Rashba model. The
conclusion is conditional on the model, endpoint uncertainty, and empirical
binary64 numerical allowances. These allowances are not statistical confidence
intervals or rigorous directed-rounding bounds. The calculation does not identify
the replacement phase or establish achievable device current asymmetry.

## Quick start

Use Python 3.12 on Linux (including WSL). Create a virtual environment:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/verify_package.py
.venv/bin/python scripts/check_results.py
.venv/bin/python generate_figures.py
.venv/bin/python scripts/build_endpoint.py
.venv/bin/python -m unittest discover -s tests -v
```

The endpoint build needs `g++` (or set `CXX` to another C++ compiler).

`generate_figures.py` is the original plotting program, unchanged. It reads
`scientific_data.json` and writes seven PDF/PNG figure pairs to `figures/`.
Published figures are preserved in `paper/figures/`. Figure regeneration uses
saved measurements and does not rerun the physics. PDF metadata and rendering
may vary with the platform; numerical inputs are checksum-bound.

`scripts/check_results.py` checks the saved primary/independent agreement,
negative linear margins, original quarter-error criteria, nonlinear finite
amplitude margins, and angular-origin comparisons. It prints the actual values
and fails if an acceptance check fails. It is an audit of saved evidence.

## Reproduce calculations

### Independent linear reassembly

```bash
.venv/bin/python scripts/reproduce.py linear-reassembly --output outputs/linear-reassembly
```

This runs the exact accepted call-271 calculation: it reassembles saved call-264
finite-grid results with the corrected infinite-frequency tail and call-256
source allowance. It does **not** generate fresh finite-grid measurements.

For fresh finite-grid linear calculations with the historical call-264 protocol:

```bash
.venv/bin/python scripts/reproduce.py fresh-linear --output outputs/fresh-linear
```

The call-264 runner still reuses the bundled call-256 source budget. Apply the
call-271 tail corrections when comparing with the article; call-264 alone is
not the final published qualification. Its original source versions, reference
kernel, and protocols are preserved together in `code/linear/fresh264/`.

### Nonlinear baseline and shifted angular origin

```bash
.venv/bin/python scripts/reproduce.py nonlinear-baseline --output outputs/nonlinear-baseline
.venv/bin/python scripts/reproduce.py nonlinear-shifted --output outputs/nonlinear-shifted
```

These execute the complete original numerical schedules. They can take hours;
the recorded shifted run took 14,348 wall seconds (about four hours) with a
process lifetime peak RSS of about 1.5 GiB. This is historical measurement,
not a resource guarantee. The reruns preserve the original reuse of completed
pilot/control data and the original accepted linear anchor. Each mode reads its
own immutable input package; running the baseline does not silently replace
the baseline bundled for the shifted calculation.

The standalone wrapper uses a local wall deadline (24 hours by default;
`--wall-seconds` overrides it), requires a new output directory, and never
accesses the original research workspace. It does not recreate the research
harness's admission, held-out reservation, or controller authority. A rerun of
these already inspected protocols is numerical reproduction, not fresh held-out
scientific confirmation. Original native worker/guard adapters are preserved for
provenance; their engine-dependent entry points are not standalone launchers.

### Endpoint and original analytic kernel

The original source is in `code/endpoint/`, with saved evidence in
`data/endpoint/`. Build the compensated summation library with a C++ compiler:

```bash
.venv/bin/python scripts/build_endpoint.py
```

This builds `sums.so` from `sums.cpp`; no precompiled scientific binaries are
distributed. `uniform.py` supplies the uniform free-energy derivatives and
`primary.py` supplies the original analytic Gaussian response. `endpoint_stage.py`,
`recover.py`, `topology.py`, `screen.py`, and the selected-control scripts retain
their original campaign-relative file conventions. They are archival stage
scripts rather than a single end-to-end launcher. Inspect their frozen protocol
in `provenance/claims/000008.json` before using them. Smoke tests import and
evaluate the kernel after compilation.

The initial call-8 controls were inconclusive; `certificate_audit.json` preserves
that failure. Later source/error qualifications resolved the article's accepted
criteria. The historical failed or partial records must not be read as final
successful confirmations. `selected_checks_executed_reconstruction.py` is
preserved alongside `selected_checks.py` to retain the recorded executed-source
reconstruction.

## Contents and provenance

| Path | Purpose |
| --- | --- |
| `code/endpoint/` | Endpoint reconstruction, analytic primary kernel, compensated sums, historical stage scripts |
| `code/primary_reference.py` | Source-compatible analytic reference implementation |
| `code/linear/source256/` | Continuous source/endpoint allowance calculation |
| `code/linear/fresh264/` | Independent finite-grid response and source controls |
| `code/linear/accepted271/` | Accepted corrected-tail reassembly |
| `code/nonlinear/pilot293/` | Completed scoped pilot and continuous coverage controls |
| `code/nonlinear/partial299/` | Historical partial run, preserved as partial |
| `code/nonlinear/baseline308/` | Completed baseline angular calculation |
| `code/nonlinear/shifted330/` | Completed shifted angular confirmation calculation |
| `data/evidence/` | Published supplement, frozen endpoint/witness, protocols, compact review/source reports |
| `data/runs/` | Original saved qualification results |
| `scientific_data.json` | Exact figure data assembled for the manuscript |
| `paper/` | Published PDF, editable TeX, bibliography, original figures and article README |
| `provenance/sources.json` | Mapping from every copied file to its original path and SHA-256 |
| `provenance/claims/` | Original frozen scientific claims and execution input bindings |
| `MANIFEST.json` | SHA-256 and byte size for every distributed payload |

Copied scientific code, input files, and evidence retain their original bytes.
The root README, dependency list, citation metadata, verification scripts,
standalone reproduction wrapper, and smoke tests are release additions. The
full autonomous harness, private credentials, operational logs, Python caches,
precompiled libraries, and unrelated research are outside this release.
Historical provenance paths identify the source archive; reproduction uses
paths within this repository. Archival manuscript metadata retains its original
URL spelling; the article link above is the current working URL.

## Citation and licensing

Please cite Colin O’Callaghan, *Spatial Pairing Instability at a Rashba
Superconducting Diode Endpoint* (7 October 2026), and the article URL above.
`CITATION.cff` provides GitHub-readable citation metadata. Contact:
colin@colinocallaghan.com.

The manuscript's existing licence is Creative Commons Attribution 4.0
International; see `paper/README.txt`. No additional software licence has been
selected for this GitHub preparation package. Do not assume the manuscript's
licence automatically licenses the code.

## Validation of this release

See `VALIDATION.json` for the actual checks performed. Full multi-hour numerical
calculations are not rerun merely to package the code. The smoke tests and saved
evidence checks do not independently certify the scientific conclusion.
