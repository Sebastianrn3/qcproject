from pathlib import Path
import numpy as np
import subprocess
from scipy.optimize import  minimize

#obligatory parameters
MULTI = "SINGLET"
CHARGE = -1
NAME_XYZ = "ClCH3F_TS_preopt.xyz"

#directories
XYZ_PATH = Path("XYZs")
MOPAC_PATH =  Path.cwd().parents[2]/"mopac_runs"
MOPAC_EXE_PATH = Path.cwd().parents[2]/"mopac/bin/mopac.exe"

KCALMOL_TO_EV = 0.0433641153087705
EV_TO_HARTREE = 1 / 27.211386245988
HARTREE_PER_KCALMOL = 1 / 627.5094740631
ANGSTROM_PER_BOHR = 0.529177210903

jobname = NAME_XYZ[:-4]

def main():
    atoms, geometry = import_xyz()
    create_mopac_input(atoms, geometry)
    run_mopac_exe()
    return parse_mopac_output()

def import_xyz():
    with open(XYZ_PATH/NAME_XYZ) as f:
        xyz = f.read()
    xyz_array = xyz.split()
    block = np.array(xyz_array[1:]).reshape(int(xyz_array[0]), 4)
    atoms = block[:, 0]
    geometry = block[:, 1:].astype(float)
    return atoms, geometry

def create_mopac_input(atoms, geometry):
    input_text=f"""AM1 1SCF XYZ GRADIENTS AUX DCART GEO-OK {MULTI} CHARGE={CHARGE}
    {NAME_XYZ}\n\n"""

    for n in range(len(atoms)):
        x = geometry[n]
        input_text += f"{atoms[n]} {x[0]} {x[1]} {x[2]}\n"
    with open(MOPAC_PATH/f"{jobname}.mop", "w") as f:
        f.write(input_text)

def run_mopac_exe():
    mop_input = MOPAC_PATH / f"{jobname}.mop"
    subprocess.run(
        [str(MOPAC_EXE_PATH), mop_input.name],
        cwd=str(MOPAC_PATH),
        text=True,
        check=False
    )

def parse_mopac_output():
    with open(MOPAC_PATH / f"{jobname}.aux") as f:
        output = f.read()

    read_from = output.find("HEAT_OF_FORMATION:KCAL/MOL=")
    read_to = output[read_from:].find(f"\n")+read_from
    E_kcal_mol = float(output[read_from:read_to].replace("D", "E").split("=")[1])
    E_Eh = E_kcal_mol * HARTREE_PER_KCALMOL

    read_from = output.find("GRADIENTS:KCAL/MOL/ANGSTROM")
    read_to = output[read_from:].find("OVERLAP_MATRIX") + read_from
    gradients_kcal_per_A = np.array(output[read_from:read_to].split()[1:], dtype=float)
    gradients_Eh_per_bohr = gradients_kcal_per_A * HARTREE_PER_KCALMOL * ANGSTROM_PER_BOHR

    print(f"Heat of Formation={E_Eh} Eh")
    print(f"Gradients (Eh/bohr)={gradients_Eh_per_bohr}")

    return E_Eh, gradients_Eh_per_bohr

def sci_minimize():
    def rewrite_xyz(x_bohr):
        x_ang = x_bohr * ANGSTROM_PER_BOHR
        with open(XYZ_PATH/ NAME_XYZ) as f:
            lines = [ln.strip() for ln in f.readlines()]
        n = int(lines[0])
        new_xyz = f"{n}\n\n"
        for i in range(n):
            symbol = lines[2 + i].split()
            symbol[1], symbol[2], symbol[3] = x_ang[i*3:i*3+3]
            new_xyz += " ".join(str(x) for x in symbol) + f"\n"
        with open(XYZ_PATH / NAME_XYZ, "w") as f:
            f.write(new_xyz)

    def run_scf(x_bohr):
        rewrite_xyz(x_bohr)
        E_Eh, G_Eh_per_bohr = main()
        return E_Eh, G_Eh_per_bohr

    def optimize(x_start):
        return  minimize(
            fun = run_scf,
            x0 = x_start,
            jac = True
        )

    x0_ang = import_xyz()[1].ravel()
    x0_bohr = x0_ang / ANGSTROM_PER_BOHR

    res = optimize(x0_bohr)
    rewrite_xyz(res.x)
    print(res)

sci_minimize()
