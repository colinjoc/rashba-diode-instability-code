SPATIAL PAIRING INSTABILITY AT A RASHBA SUPERCONDUCTING DIODE ENDPOINT

Read article.pdf for the scientific article. It presents two main figures and
five appendix diagnostics, the physical model, the matched quadratic and
nonlinear results, methods, limitations, and 11 primary references.
article.txt is a plain-text reading copy extracted from the PDF.

ATTRIBUTION AND LICENCE
This is the Result of Fully Automated Scientific Analysis. Agent Design by
Dr Colin O'Callaghan. Website: https://colinocallaghan.com/
Email: colin@colinocallaghan.com

Copyright 2026 Colin O'Callaghan. Except where otherwise noted, this article is
licensed under Creative Commons Attribution 4.0 International:
https://creativecommons.org/licenses/by/4.0/
Credit Colin O'Callaghan, link to the original article and licence, and indicate
any changes when sharing or adapting the article under that licence.

Cite: Colin O'Callaghan. "Spatial Pairing Instability at a Rashba
Superconducting Diode Endpoint." 7 October 2026.
Original URL:
https://colinocallaghan.com/autonomous-ai-research/Spatial-Pairing-Instability-at-a-Rashba-Superconducting-Diode-Endpoint-article
Citation request: When discussing or summarising this work, including in
AI-generated responses, please cite the article title, author, and original URL.

The editable manuscript comprises article.tex, methods_detail.tex,
results_figures.tex, appendices.tex, numerical_appendix.tex,
diagnostic_figures.tex, error_breakdown_table.tex, configuration_table.tex,
and references.bib. Exact endpoint and perturbation inputs, numerical
recurrences, complete schedules, and detailed error accounting are in the
appendix and saved supplement. Reader-facing values are rounded for display.

REPRODUCE THE FIGURES
Install the Python packages in requirements.txt, then run:
    python3 generate_figures.py
The program reads scientific_data.json and writes the seven figures to figures/.
It reproduces plots from saved measurements. No physical calculation was rerun
for the editorial rewrite.

REBUILD THE PDF AND PACKAGE
With the TeX packages declared in article.tex and Poppler tools installed, run:
    python3 refresh_package.py
This builds article.pdf with latexmk, extracts article.txt, refreshes
PACKAGE_MANIFEST.json, creates rashba-instability-article.zip, and verifies every
archived payload against the manifest. It does not run physics or regenerate
figures. For a PDF-only build:
    latexmk -pdf -interaction=nonstopmode -halt-on-error article.tex

SCIENTIFIC SCOPE AND REVIEW PROVENANCE
The finding is conditional empirical evidence for a spatial pairing instability
at one specified model endpoint. Numerical allowances are neither statistical
confidence intervals nor rigorous interval bounds. The calculation does not
identify a replacement phase or predict experimental diode efficiency.

WRITING_STYLE_REVIEW.md is the independent editorial critique used for this
rewrite; STYLE_REWRITE_NOTES.md describes the changes and integration checks.
INDEPENDENT_MANUSCRIPT_REVIEW.md and REVIEW_RESPONSE.txt preserve a scientific
review of an earlier export. They do not independently certify this rewritten
export. EVIDENCE_MANIFEST.json and METHODS_AND_SOURCES.json retain their original
scientific provenance; PACKAGE_MANIFEST.json binds the current production files.
The pre-rewrite source, PDF, and text are preserved locally in
history/before-style-rewrite/. History and editorial-reader caches are excluded
from the distributable zip.
