param(
    [string]$Source = (Join-Path $PSScriptRoot '..\ALOD_paper\helmholtz_lod_exact_solution_certified_amsart_final.tex')
)

$ErrorActionPreference = 'Stop'
$projectRoot = $PSScriptRoot
$sourcePath = [System.IO.Path]::GetFullPath($Source)
$utf8NoBom = [System.Text.UTF8Encoding]::new($false)

if (-not (Test-Path -LiteralPath $sourcePath -PathType Leaf)) {
    throw "Source file not found: $sourcePath"
}

$sourceText = [System.IO.File]::ReadAllText($sourcePath)
$sourceHashBefore = (Get-FileHash -LiteralPath $sourcePath -Algorithm SHA256).Hash

function Write-Utf8File {
    param(
        [Parameter(Mandatory)] [string]$Path,
        [Parameter(Mandatory)] [string]$Content
    )

    $parent = Split-Path -Parent $Path
    if ($parent) {
        [System.IO.Directory]::CreateDirectory($parent) | Out-Null
    }
    [System.IO.File]::WriteAllText($Path, $Content, $utf8NoBom)
}

function Normalize-Fragment {
    param([Parameter(Mandatory)] [string]$Text)
    return $Text.Trim("`r", "`n") + "`r`n"
}

$documentClassMatch = [regex]::Match($sourceText, '(?m)^\\documentclass[^\r\n]*(?:\r?\n)')
$beginDocumentMatch = [regex]::Match($sourceText, '(?m)^\\begin\{document\}\s*(?:\r?\n)')
$endDocumentMatch = [regex]::Match($sourceText, '(?m)^\\end\{document\}\s*$')
$bibliographyMatch = [regex]::Match($sourceText, '(?m)^\\begin\{thebibliography\}')
$sectionMatches = [regex]::Matches($sourceText, '(?m)^\\section\{')

if (-not $documentClassMatch.Success -or -not $beginDocumentMatch.Success -or
    -not $endDocumentMatch.Success -or -not $bibliographyMatch.Success) {
    throw 'The source does not contain the expected document markers.'
}
if ($sectionMatches.Count -ne 7) {
    throw "Expected 7 top-level sections, found $($sectionMatches.Count)."
}

$sections = @(
    @{ Number = '01'; Slug = 'introduction' },
    @{ Number = '02'; Slug = 'model-and-two-level-setting' },
    @{ Number = '03'; Slug = 'lod-discretization-and-diagnostics' },
    @{ Number = '04'; Slug = 'exact-solution-certification' },
    @{ Number = '05'; Slug = 'adaptive-algorithm-and-cost' },
    @{ Number = '06'; Slug = 'numerical-experiments' },
    @{ Number = '07'; Slug = 'conclusions' }
)

$preambleStart = $documentClassMatch.Index + $documentClassMatch.Length
$preambleLength = $beginDocumentMatch.Index - $preambleStart
$preamble = Normalize-Fragment $sourceText.Substring($preambleStart, $preambleLength)

$frontmatterStart = $beginDocumentMatch.Index + $beginDocumentMatch.Length
$frontmatterLength = $sectionMatches[0].Index - $frontmatterStart
$frontmatter = Normalize-Fragment $sourceText.Substring($frontmatterStart, $frontmatterLength)

Write-Utf8File (Join-Path $projectRoot 'preamble.tex') $preamble
Write-Utf8File (Join-Path $projectRoot 'frontmatter.tex') $frontmatter

for ($index = 0; $index -lt $sections.Count; $index++) {
    $start = $sectionMatches[$index].Index
    if ($index -lt $sections.Count - 1) {
        $end = $sectionMatches[$index + 1].Index
    } else {
        $end = $bibliographyMatch.Index
    }
    $fragment = Normalize-Fragment $sourceText.Substring($start, $end - $start)
    $fileName = "$($sections[$index].Number)-$($sections[$index].Slug).tex"
    Write-Utf8File (Join-Path $projectRoot "sections\$fileName") $fragment
}

$bibliographyLength = $endDocumentMatch.Index - $bibliographyMatch.Index
$bibliography = Normalize-Fragment $sourceText.Substring($bibliographyMatch.Index, $bibliographyLength)
Write-Utf8File (Join-Path $projectRoot 'backmatter\references.tex') $bibliography

$inputLines = foreach ($section in $sections) {
    "\input{sections/$($section.Number)-$($section.Slug)}"
}
$includeLines = foreach ($section in $sections) {
    "\include{sections/$($section.Number)-$($section.Slug)}"
}

$main = @"
\documentclass[11pt]{amsart}

\input{preamble}

\begin{document}

\input{frontmatter}
$($inputLines -join "`r`n")
\input{backmatter/references}

\end{document}
"@
Write-Utf8File (Join-Path $projectRoot 'main.tex') (Normalize-Fragment $main)

$sectionMaster = @"
\documentclass[11pt]{amsart}

% Partial PDFs cannot contain hyperlink destinations from omitted sections.
% Draft mode keeps reference text while suppressing those dead-link warnings.
\ifdefined\ALODIncludeOnly
  \PassOptionsToPackage{draft}{hyperref}
\fi

\input{preamble}

% Section drivers define this macro before loading section-master.tex.
\ifdefined\ALODIncludeOnly
  \expandafter\includeonly\expandafter{\ALODIncludeOnly}
\fi

\begin{document}

\include{frontmatter}
$($includeLines -join "`r`n")
\include{backmatter/references}

\end{document}
"@
Write-Utf8File (Join-Path $projectRoot 'section-master.tex') (Normalize-Fragment $sectionMaster)

foreach ($section in $sections) {
    $baseName = "section-$($section.Number)-$($section.Slug)"
    $includeName = "sections/$($section.Number)-$($section.Slug)"
    $driver = @"
% Compile through build.ps1 so that the full-document auxiliary files exist.
\def\ALODIncludeOnly{$includeName}
\input{section-master.tex}
"@
    Write-Utf8File (Join-Path $projectRoot "$baseName.tex") (Normalize-Fragment $driver)
}

$buildScript = @'
param(
    [ValidateSet('main', 'all', '01', '02', '03', '04', '05', '06', '07')]
    [string]$Target = 'main'
)

$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

$drivers = [ordered]@{
    '01' = 'section-01-introduction.tex'
    '02' = 'section-02-model-and-two-level-setting.tex'
    '03' = 'section-03-lod-discretization-and-diagnostics.tex'
    '04' = 'section-04-exact-solution-certification.tex'
    '05' = 'section-05-adaptive-algorithm-and-cost.tex'
    '06' = 'section-06-numerical-experiments.tex'
    '07' = 'section-07-conclusions.tex'
}

if (-not (Get-Command latexmk -ErrorAction SilentlyContinue)) {
    throw 'latexmk was not found on PATH.'
}

[System.IO.Directory]::CreateDirectory((Join-Path $PSScriptRoot 'build')) | Out-Null
[System.IO.Directory]::CreateDirectory((Join-Path $PSScriptRoot 'build\sections')) | Out-Null
[System.IO.Directory]::CreateDirectory((Join-Path $PSScriptRoot 'build\backmatter')) | Out-Null

function Invoke-Latexmk {
    param([Parameter(Mandatory)] [string]$EntryPoint)

    & latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build $EntryPoint
    if ($LASTEXITCODE -ne 0) {
        throw "latexmk failed for $EntryPoint (exit code $LASTEXITCODE)."
    }
}

Invoke-Latexmk 'main.tex'

if ($Target -eq 'all') {
    # This internal full build supplies the per-include checkpoints,
    # bibliography labels, and cross-section labels used by partial builds.
    Invoke-Latexmk 'section-master.tex'
    foreach ($driver in $drivers.Values) {
        Invoke-Latexmk $driver
    }
} elseif ($Target -ne 'main') {
    Invoke-Latexmk 'section-master.tex'
    Invoke-Latexmk $drivers[$Target]
}
'@
Write-Utf8File (Join-Path $projectRoot 'build.ps1') (Normalize-Fragment $buildScript)

$readme = @'
# ALOD paper — sectioned build

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
'@
Write-Utf8File (Join-Path $projectRoot 'README.md') (Normalize-Fragment $readme)

$sourceHashAfter = (Get-FileHash -LiteralPath $sourcePath -Algorithm SHA256).Hash
if ($sourceHashAfter -ne $sourceHashBefore) {
    throw 'The source hash changed while generating the sectioned project.'
}

Write-Host "Generated sectioned project at $projectRoot"
Write-Host "Source SHA256 (unchanged): $sourceHashAfter"
