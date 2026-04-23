$ErrorActionPreference = "Stop"

# Resolve the package root from this script location.
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$PackageDir = (Resolve-Path (Join-Path $ScriptDir "..\..\..")).Path
$PackageDirForR = $PackageDir -replace "\\", "/"

# Add a standard Windows R bin path to PATH if R is installed but Rscript is not yet visible.
# Sometimes the R installer does not add Rscript to the PATH, which is required for installing the R dependencies.
$RRoot = "C:\Program Files\R"
if (Test-Path $RRoot) {
  $RBinCandidates = Get-ChildItem $RRoot -Directory | ForEach-Object {
    @(
      (Join-Path $_.FullName "bin\x64"),
      (Join-Path $_.FullName "bin")
    )
  }

  $UserPath = [Environment]::GetEnvironmentVariable("Path", "User")
  $PathEntries = if ([string]::IsNullOrWhiteSpace($UserPath)) { @() } else { $UserPath -split ";" }

  foreach ($RBinPath in $RBinCandidates) {
    if ((Test-Path (Join-Path $RBinPath "Rscript.exe")) -and -not ($PathEntries -contains $RBinPath)) {
      $PathEntries = @($RBinPath) + $PathEntries
      break
    }
  }

  $NewUserPath = ($PathEntries | Select-Object -Unique) -join ";"
  [Environment]::SetEnvironmentVariable("Path", $NewUserPath, "User")
  $env:Path = "$NewUserPath;$env:Path"
}

# Check that Rscript is available.
if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) {
  throw "Rscript is required but was not found. Please ensure R is installed before running this script."
}

# Check that winget is available for dependency installs.
if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
  throw "winget is required for this script. Install App Installer from Microsoft Store or install Rtools and Git manually."
}

# Install Windows development dependencies.
winget install --exact --id RProject.Rtools --accept-package-agreements --accept-source-agreements
winget install --exact --id Git.Git --accept-package-agreements --accept-source-agreements

# Add common Rtools and Git locations to the user PATH.
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

# Download and install PLINK into the user's bin directory.
$InstallDir = Join-Path $HOME "bin"
$ArchivePath = Join-Path $env:TEMP "plink_win64_20250819.zip"
$ExtractDir = Join-Path $env:TEMP "plink_win64_20250819"
New-Item -ItemType Directory -Force -Path $InstallDir | Out-Null
Remove-Item $ExtractDir -Recurse -Force -ErrorAction SilentlyContinue
Invoke-WebRequest -Uri "https://s3.amazonaws.com/plink1-assets/plink_win64_20250819.zip" -OutFile $ArchivePath
Expand-Archive -Path $ArchivePath -DestinationPath $ExtractDir -Force
Copy-Item `
  (Get-ChildItem $ExtractDir -Recurse -Filter "plink.exe" | Select-Object -First 1).FullName `
  (Join-Path $InstallDir "plink.exe") `
  -Force
Remove-Item $ArchivePath -Force
Remove-Item $ExtractDir -Recurse -Force
if (-not ($PathEntries -contains $InstallDir)) {
  $PathEntries = @($InstallDir) + $PathEntries
}
$NewUserPath = ($PathEntries | Select-Object -Unique) -join ";"
[Environment]::SetEnvironmentVariable("Path", $NewUserPath, "User")
$env:Path = "$InstallDir;$env:Path"
& (Join-Path $InstallDir "plink.exe") --version

# Install the R dependencies required by HapSelect.
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org'); remotes::install_deps('$PackageDirForR', dependencies = TRUE)"

# Install the genomicSimulation R dependency from GitHub.
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org'); remotes::install_github('vllrs/genomicSimulation')"
