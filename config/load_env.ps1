param([string]$EnvFile = "$PSScriptRoot/.env")

if (!(Test-Path $EnvFile)) { throw "Config not found: $EnvFile" }

Get-Content $EnvFile | ForEach-Object {
  if ($_ -match '^\s*#') { return }
  if ($_ -match '^\s*$') { return }
  if ($_ -match '^\s*([^=]+)\s*=\s*(.*)$') {
    $name  = $matches[1].Trim()
    $value = $matches[2].Trim()

    [System.Environment]::SetEnvironmentVariable($name, $value, 'Process')
  }
}

$WSL_PROJECT_ROOT = (wsl wslpath -a "$env:WIN_PROJECT_ROOT").Trim()
[System.Environment]::SetEnvironmentVariable('WSL_PROJECT_ROOT', $WSL_PROJECT_ROOT, 'Process')

$WSL_REACTIONS    = "$WSL_PROJECT_ROOT/$env:REACTIONS_DIR"
$WSL_REACTION_DIR = "$WSL_REACTIONS/$env:REACTION"
$WSL_INPUTS       = "$WSL_REACTION_DIR/$env:INPUTS_SUBDIR"
$WSL_RUNS         = "$WSL_REACTION_DIR/$env:RUNS_SUBDIR"

$set = @{
  WSL_REACTIONS    = $WSL_REACTIONS
  WSL_REACTION_DIR = $WSL_REACTION_DIR
  WSL_INPUTS       = $WSL_INPUTS
  WSL_RUNS         = $WSL_RUNS
}
$set.GetEnumerator() | ForEach-Object {
  [System.Environment]::SetEnvironmentVariable($_.Key, $_.Value, 'Process')
}