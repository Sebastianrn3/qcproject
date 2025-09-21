@echo off
setlocal EnableExtensions EnableDelayedExpansion

REM
for /f "usebackq tokens=1,* delims== eol=#" %%A in ("%~dp0.env") do (
  if not "%%~A"=="" set "%%~A=%%~B"
)

REM
if "%WIN_PROJECT_ROOT%"=="" (
  echo [ERROR] WIN_PROJECT_ROOT is not set in config/.env
  exit /b 1
)

REM
for /f "usebackq delims=" %%P in (`wsl wslpath -a "%WIN_PROJECT_ROOT%"`) do set "WSL_PROJECT_ROOT=%%P"

REM
set "WSL_REACTIONS=%WSL_PROJECT_ROOT%/%REACTIONS_DIR%"
set "WSL_REACTION_DIR=%WSL_REACTIONS%/%REACTION%"
set "WSL_INPUTS=%WSL_REACTION_DIR%/%INPUTS_SUBDIR%"
set "WSL_RUNS=%WSL_REACTION_DIR%/%RUNS_SUBDIR%"

endlocal & (
  set "WIN_PROJECT_ROOT=%WIN_PROJECT_ROOT%"
  set "WSL_DISTRO=%WSL_DISTRO%"
  set "WSL_NWCHEM_BIN=%WSL_NWCHEM_BIN%"
  set "WSL_MOPAC7_BIN=%WSL_MOPAC7_BIN%"
  set "WSL_SCRATCH=%WSL_SCRATCH%"
  set "WSL_PROJECT_ROOT=%WSL_PROJECT_ROOT%"
  set "WSL_REACTIONS=%WSL_REACTIONS%"
  set "WSL_REACTION_DIR=%WSL_REACTION_DIR%"
  set "WSL_INPUTS=%WSL_INPUTS%"
  set "WSL_RUNS=%WSL_RUNS%"
  set "REACTION=%REACTION%"
  set "ENV_GAS=%ENV_GAS%"
)