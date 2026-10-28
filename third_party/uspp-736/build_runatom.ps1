<#
Builds bin\runatom.exe and bin\reform.exe from the unmodified USPP 7.3.6 source
(Source\*.f for runatom, Utils\reform.f for reform, the converter that writes the
formatted pseudopotential CPMD reads).

Needs MSYS2 UCRT64 with gcc-fortran and binutils (see PROVENANCE.md).
The idate patch is applied to a temporary copy; Source\ is never changed.
The result is statically linked: it may import only KERNEL32.dll and the
Windows 10+ UCRT (api-ms-win-crt-*.dll), or the build fails.

Usage (from venv\src):
    powershell -ExecutionPolicy Bypass -File third_party\uspp-736\build_runatom.ps1 [-Msys C:\msys64]
#>
param([string]$Msys = "C:\msys64")
$ErrorActionPreference = "Stop"

$root = $PSScriptRoot
$bin = Join-Path $Msys "ucrt64\bin"
$gfortran = Join-Path $bin "gfortran.exe"
$objdump = Join-Path $bin "objdump.exe"
if (-not (Test-Path $gfortran)) { throw "gfortran not found at ${gfortran}: install mingw-w64-ucrt-x86_64-gcc-fortran" }
if (-not (Test-Path $objdump)) { throw "objdump not found at ${objdump}: install mingw-w64-ucrt-x86_64-binutils" }
$env:PATH = "$bin;$env:PATH"
$flags = @("-O1", "-std=legacy", "-fallow-argument-mismatch", "-static")

$work = Join-Path ([IO.Path]::GetTempPath()) ("uspp-build-" + [guid]::NewGuid())
New-Item -ItemType Directory $work | Out-Null
try {
    Copy-Item (Join-Path $root "Source\*.f") $work
    $main = Join-Path $work "runatom.f"
    $text = [IO.File]::ReadAllText($main)
    $pattern = "(?m)^      external idate(\r?)$"
    if ([regex]::Matches($text, $pattern).Count -ne 1) { throw "the idate patch no longer matches runatom.f" }
    [IO.File]::WriteAllText($main, [regex]::Replace($text, $pattern, 'c     external idate$1'))

    Copy-Item (Join-Path $root "Utils\reform.f") $work      # unmodified: it needs no patch

    function Assert-SelfContained([string]$exe) {
        $dlls = & $objdump -p $exe | Select-String "DLL Name:" | ForEach-Object { ($_ -split ":", 2)[1].Trim() }
        $foreign = @($dlls | Where-Object { $_ -ne "KERNEL32.dll" -and $_ -notlike "api-ms-win-crt-*" })
        if ($foreign.Count -gt 0) { throw "$exe is not self-contained; it imports: $($foreign -join ', ')" }
        return ($dlls -join ', ')
    }

    Push-Location $work
    try {
        $sources = Get-ChildItem -Filter *.f | Where-Object Name -ne "reform.f" | ForEach-Object Name
        & $gfortran @flags -o runatom.exe @sources
        if ($LASTEXITCODE -ne 0) { throw "gfortran failed on runatom with exit code $LASTEXITCODE" }
        & $gfortran @flags -o reform.exe reform.f
        if ($LASTEXITCODE -ne 0) { throw "gfortran failed on reform with exit code $LASTEXITCODE" }
        $imports = Assert-SelfContained "runatom.exe"
        $reformImports = Assert-SelfContained "reform.exe"
    } finally { Pop-Location }

    # Smoke run: the package's own example must give a pseudopotential and its CPMD formatted copy.
    $smoke = Join-Path $work "smoke"
    Copy-Item (Join-Path $root "Work\000-Ex\000-Ex-ca-sp-test") $smoke -Recurse
    Push-Location $smoke
    $savedPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"   # Windows PowerShell 5.1 turns native stderr output into errors; the file checks below decide success
    try {
        & (Join-Path $work "runatom.exe") ti_ae_s2d1.adat ti_ae_s2d1.out ti_ae_s2d1.ae ti_ae_s2d1.atwf ti_ae_s2d1.logd dummy 2>$null
        & (Join-Path $work "runatom.exe") ti_ps.adat ti_ps.out ti_ae_s2d1.ae ti_ps.atwf ti_ps.logd ti_ps.uspp 2>$null
        if (-not (Test-Path ti_ps.uspp)) { throw "smoke run produced no ti_ps.uspp; see $smoke" }
        # stdin comes from a file: piping from Windows PowerShell 5.1 makes reform's list-directed read fail
        [IO.File]::WriteAllText((Join-Path $smoke "reform.in"), "1`nti_ps.uspp`nti_ps_formatted.uspp`n")
        cmd /c "`"$(Join-Path $work 'reform.exe')`" < reform.in" 2>&1 | Out-Null
        if (-not (Test-Path ti_ps_formatted.uspp)) { throw "reform.exe produced no ti_ps_formatted.uspp; see $smoke" }
    } finally { $ErrorActionPreference = $savedPreference; Pop-Location }

    $out = Join-Path $root "bin"
    New-Item -ItemType Directory $out -Force | Out-Null
    $lines = @()
    foreach ($name in "runatom.exe", "reform.exe") {
        Copy-Item (Join-Path $work $name) $out -Force
        $hash = (Get-FileHash (Join-Path $out $name) -Algorithm SHA256).Hash.ToLower()
        [IO.File]::WriteAllText((Join-Path $out "$name.sha256"), "$hash  $name`n")
        $lines += "$name sha256: $hash"
    }
    $version = (& $gfortran --version | Select-Object -First 1)
    [IO.File]::WriteAllText((Join-Path $out "BUILDINFO.txt"),
        "compiler: $version`nflags: $($flags -join ' ')`n" +
        "runatom.exe: Source/*.f + patches/0001-gfortran-intrinsic-idate.patch; imports: $imports`n" +
        "reform.exe: Utils/reform.f (unmodified); imports: $reformImports`n" + ($lines -join "`n") + "`n")
    Write-Host ($lines -join "`n")
} finally {
    Remove-Item $work -Recurse -Force -ErrorAction SilentlyContinue
}
