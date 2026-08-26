# ALOD paper 鈥?sectioned build

This directory is generated from
`../ALOD_paper/helmholtz_lod_exact_solution_certified_amsart_final.tex`.
The original source file is read only and is never modified.

## Layout

- `main.tex`: full-paper entry point.
- `section-master.tex`: internal `\includeonly` entry point used by partial builds.
- `preamble.tex`: packages, theorem declarations, macros, and metadata.
- `frontmatter.tex`: abstract and title material.
- `sections/*.tex`: one file per top-level section.
- `backmatter/references.tex`: the original inline bibliography.
- `section-*.tex`: lightweight entry points for partial builds.
- `build.ps1`: reproducible full/partial build command.
- `generate.ps1`: regenerates this layout from the unchanged monolithic source.

## Build

Run the commands from this directory in PowerShell:

```powershell
# Full paper
.\build.ps1 main

# A single section (01 through 07)
.\build.ps1 04

# Full paper and every section PDF
.\build.ps1 all
```

Output is written to `build/`.  A section build first refreshes the internal
`section-master.tex`, because LaTeX's `\includeonly` mechanism uses its
auxiliary files to preserve section/equation/theorem numbering, citations, and
cross-section references.  Each section PDF contains only that section.
`main.tex` uses `\input`, so the full paper retains the monolithic source's page
flow instead of forcing every section onto a new page.

To refresh the split after intentionally editing the monolithic source, run:

```powershell
.\generate.ps1
```
