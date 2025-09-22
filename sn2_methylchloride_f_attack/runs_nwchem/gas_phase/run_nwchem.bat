@echo off
setlocal

rem
set "WSL_INPUTS=/mnt/c/Users/ddizy/Desktop/qcproject/reactions/sn2_methylchloride_f_attack/inputs_nwchem/gas_phase"
set "WSL_RUNS=/mnt/c/Users/ddizy/Desktop/qcproject/reactions/sn2_methylchloride_f_attack/runs_nwchem/gas_phase"

rem
wsl -e bash -lc "cd '%WSL_RUNS%' && nwchem '%WSL_INPUTS%/f_energy.nw' > f_energy.out" &&
rem
wsl -e bash -lc "cd '%WSL_RUNS%' && nwchem '%WSL_INPUTS%/cl_energy.nw' > cl_energy.out" &&
rem
wsl -e bash -lc "cd '%WSL_RUNS%' && nwchem '%WSL_INPUTS%/ch3cl_opt_freq.nw' > ch3cl_opt_freq.out" &&
rem
wsl -e bash -lc "cd '%WSL_RUNS%' && nwchem '%WSL_INPUTS%/ch3f_opt_freq.nw' > ch3f_opt_freq.out" &&

rem
wsl -e bash -lc "rm -rf '%WSL_RUNS%/scratch' '%WSL_RUNS%/perm'"

echo All done. Press any key to exit
pause >nul