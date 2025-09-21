@echo off
set WSL_DIR=/mnt/c/Users/ddizy/Desktop/qcproject/tests
wsl -e bash -lc "cd '%WSL_DIR%' && nwchem f_energy.nw > f_energy.out"
echo Done. Press any key...
pause >nul