from pathlib import Path
import numpy as np
import subprocess, os
from scipy.optimize import  minimize

#obligatory parameters
MULTI = "SINGLET"
CHARGE = 0
NAME_XYZ = "CH4.xyz"

#directories
XYZ_PATH = Path("XYZs")
MOPAC7_PATH =  Path.cwd().parents[2]/"mopac7"

N=0

def main():
    atoms, geometry = import_xyz()
    create_mopac_input(atoms, geometry)
    clear_mopac_outputs()
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
    input_text=f"""AM1 1SCF XYZ GRADIENTS DEBUG DCART GEO-OK {MULTI} CHARGE={CHARGE}
    {NAME_XYZ}\n\n"""

    for n in range(len(atoms)):
        x = geometry[n]
        input_text += f"{atoms[n]} {x[0]} 0 {x[1]} 0 {x[2]} 0\n"

    with open(MOPAC7_PATH/"FOR005", "w") as f:
        f.write(input_text)

    with open(MOPAC7_PATH/f"FOR005-{N}", "w") as f:
        f.write(input_text)

def clear_mopac_outputs():
    global N
    src = Path(MOPAC7_PATH/"FOR006")  # файл в другой папке
    dst = src.with_name(f"FOR006-{N}")  # то же место, новое имя
    src.rename(dst)
    N += 1
    for file in ("FOR009", "FOR011", "FOR012", "FOR016", "SHUTDOWN"):
        if os.path.exists(MOPAC7_PATH / file):
            os.remove(MOPAC7_PATH / file)

def run_mopac_exe():
    subprocess.run([str(MOPAC7_PATH / "mopac7.exe")], cwd=str(MOPAC7_PATH),text=True, check=False)

def parse_mopac_output():
    with open(MOPAC7_PATH / "FOR006") as f:
        output = f.read()

    read_from = output.find("TOTAL ENERGY")
    read_to = output[read_from:].find("EV")+read_from
    total_energy = float(output[read_from:read_to].split()[-1])

    read_from = output.find("CARTESIAN COORDINATE DERIVATIVES")
    gradients = np.array(output[read_from:].split())[8:]
    index=np.where(gradients=="-------------------------------------------------------------------------------")[0][0]
    gradients_2dim = gradients[:index].reshape(-1, 5)
    gradients_1dim = gradients_2dim[:, -3:].astype(float).ravel()

    return total_energy, gradients_1dim

def sci_minimize():
    # constants
    EV_TO_HARTREE = 1 / 27.211386245988
    HARTREE_PER_KCALMOL = 1 / 627.5094740631
    ANGSTROM_PER_BOHR = 0.529177210903
    global calclog
    calclog = ""

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
        global calclog
        rewrite_xyz(x_bohr)
        E_eV, G_kcal_per_A = main()
        E_Eh = E_eV * EV_TO_HARTREE
        G_Eh_per_bohr = G_kcal_per_A * HARTREE_PER_KCALMOL * ANGSTROM_PER_BOHR
        calclog += f"Tot energy {E_Eh} Eh\n grads, Eh/bohr:\n{G_Eh_per_bohr.reshape(-1, 3)}\n XYZ, bohr: \n{x_bohr.reshape(-1, 3)}\n**********\n"
        return E_Eh, G_Eh_per_bohr

    def lbfgsb(x_start):
        return  minimize(
            fun = run_scf,
            x0 = x_start,
            jac = True,
            method="L-BFGS-B",
            options={'iprint': 2}
        )

    x0_ang = import_xyz()[1].ravel()
    x0_bohr = x0_ang / ANGSTROM_PER_BOHR
    calclog += f"x0, bohr:\n{x0_bohr.reshape(-1, 3)}\n**********\n"

    res = lbfgsb(x0_bohr)
    print(res)
    #print(calclog)

sci_minimize()
