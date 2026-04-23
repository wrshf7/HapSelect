$ErrorActionPreference = "Stop"

if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
  throw "winget is required for this script. Install App Installer from Microsoft Store or install Rtools and Git manually."
}

Write-Host "Installing Windows system dependencies."

winget install --exact --id RProject.Rtools --accept-package-agreements --accept-source-agreements
winget install --exact --id Git.Git --accept-package-agreements --accept-source-agreements

$RtoolsBinPaths = @(
  "C:\rtools45\ucrt64\bin",
  "C:\rtools45\usr\bin",
  "C:\rtools44\ucrt64\bin",
  "C:\rtools44\usr\bin"
)

$UserPath = [Environment]::GetEnvironmentVariable("Path", "User")
$PathEntries = if ([string]::IsNullOrWhiteSpace($UserPath)) { @() } else { $UserPath -split ";" }

foreach ($PathEntry in $RtoolsBinPaths) {
  if ((Test-Path $PathEntry) -and -not ($PathEntries -contains $PathEntry)) {
    $PathEntries += $PathEntry
  }
}

if (-not ($PathEntries -contains "C:\Program Files\Git\cmd") -and (Test-Path "C:\Program Files\Git\cmd")) {
  $PathEntries += "C:\Program Files\Git\cmd"
}

$NewUserPath = ($PathEntries | Select-Object -Unique) -join ";"
[Environment]::SetEnvironmentVariable("Path", $NewUserPath, "User")

Write-Host "Installed Windows system dependencies."
Write-Host "If R is not installed yet, install it separately before running devtools::install_deps()."
