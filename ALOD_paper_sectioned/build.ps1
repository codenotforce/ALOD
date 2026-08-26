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
