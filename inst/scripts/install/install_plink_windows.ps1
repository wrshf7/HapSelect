$ErrorActionPreference = "Stop"

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

$UserPath = [Environment]::GetEnvironmentVariable("Path", "User")
if (-not (($UserPath -split ";") -contains $InstallDir)) {
  $NewUserPath = if ([string]::IsNullOrWhiteSpace($UserPath)) { $InstallDir } else { "$InstallDir;$UserPath" }
  [Environment]::SetEnvironmentVariable("Path", $NewUserPath, "User")
}

$env:Path = "$InstallDir;$env:Path"
& (Join-Path $InstallDir "plink.exe") --version
Write-Host "Installed Plink in $InstallDir. "
