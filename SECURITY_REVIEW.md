# Security review

Reviewed on 7 October 2026, before public GitHub publication.

## Result

No embedded credentials, malicious network behavior, shell injection, unsafe
deserialization, or known dependency advisories were found in the reviewed
release. Four filesystem/integrity hardening gaps in the release scripts were
fixed. These are local research programs, not an Internet-facing service.

This is a bounded security review, not a guarantee that software is free of all
vulnerabilities. Scientific accuracy is assessed separately.

## Scope and evidence

- Python AST inspection of 61 files (48 distinct source versions), manual review
  of dynamic loading, subprocess calls, native code, and file operations.
- Credential-pattern scanning of release text, extracted PDF text, and all
  committed Git-history blobs available at review time. No values matching
  private-key, GitHub, AWS, Slack, Stripe, Google API-key, or assigned-secret
  patterns were found. Public author contact details and historical provenance
  paths are intentional metadata.
- No network-library imports, direct `eval`/`exec`, shell commands, pickle
  deserialization, or unsafe YAML loading were found in released Python code.
  The two subprocess calls use argument lists with no shell: local C++ compilation
  and the numerical reproduction test.
- The C++ compensated sum consumes a contiguous float64 buffer whose length is
  supplied by the Python wrapper. It has no file or network operations. The
  archived CUDA probe only reports local device availability.
- The manuscript PDF has no JavaScript or embedded attachments according to
  Poppler's `pdfinfo` and `pdfdetach` checks.
- All 14 direct/transitive runtime package versions in the validated environment
  were queried against the OSV vulnerability API. No known advisories were
  returned. The exact query and response are in `security/dependencies-osv.json`.
  Dependencies are pinned to that reviewed set in `requirements.txt`.
- SHA-256 checks bind every copied source/evidence file to its original bytes;
  distribution checks cover the repository and ZIP payloads.

Machine-readable static findings are in `security/static-analysis.json`.

## Fixes made

1. **Manifest and provenance path confinement.** The original release verifier
   checked manifest paths lexically and did not apply the same checks to
   provenance entries. It now rejects absolute paths, traversal, noncanonical
   paths, symlinked files/parents, duplicate entries, and unlisted Python files
   in source/script/test directories. Provenance entries must also match the
   verified release manifest.
2. **Integrity before execution.** Reproduction and compilation now verify the
   release before importing calculation code or invoking the compiler. The
   saved-result audit also verifies its inputs first.
3. **Compiled output safety.** The compiler previously wrote directly to
   `sums.so`, which could follow a pre-existing symlink. The build script now
   refuses symlink outputs, compiles in a private temporary directory, and
   atomically replaces the destination after successful compilation.
4. **Output directory confinement.** Reproduction outputs within the repository
   must now live under `outputs/`, protecting release scripts and other payload
   directories as well as the original frozen scientific data.

Regression tests exercise traversal, absolute paths, symlink files/parents,
tampered content, unlisted executable sources, and unsafe provenance paths.
Original scientific source files were not modified by these fixes.

## Trust boundaries and limits

Checksums detect changes relative to a trusted manifest; they are not digital
signatures. An attacker who replaces both the source and the manifest can bypass
this kind of integrity check. Obtain the repository/archive from its trusted
publisher and use a trusted Python environment and compiler (`CXX` is an
explicit local compiler selection).

Original campaign scripts and native harness adapters are archival executable
code. Some write relative to their working directory, and the adapters require
the original controller environment. Use the documented release entry points;
the archive is not a sandbox for arbitrary untrusted protocols or data. A local
attacker able to modify the repository during execution is outside this review's
confinement guarantees.

The full multi-hour solver schedules and native dependency libraries were not
subjected to fuzzing, formal verification, or a binary audit. OSV reports known
advisories; an empty response does not prove absence of undisclosed defects or
verify package artifact provenance. Installation still uses the configured
Python package index and its normal TLS verification.
