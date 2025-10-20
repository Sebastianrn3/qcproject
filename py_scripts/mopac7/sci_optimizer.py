import numpy as np
import pandas as pd
import subprocess, os, time
from datetime import datetime
from scipy.optimize import  minimize

atoms_list = [""]

MOPAC7_PATH = "C:/Users/ddizy/Desktop/mopac7/"
RESULT_PATH = "C:/Users/ddizy/Desktop/qcproject/py_scripts/mopac7/mopac7_results/"
NAME_XYZ = "TS_FCH3Cl.xyz"
NAME_XYZ_OPTIMIZED = NAME_XYZ.replace(".xyz", "_optimized.xyz")


multiplicity =  ["INVALID_MULTIPLICITY", "SINGLET", "DOUBLET", "TRIPLET", "QUARTET", "QUINTET", "SEXTET"]
PROJECT_NAME=NAME_XYZ
COMMENT="DFT_Energy_Gradient"
CHARGE=-1
MULTI=multiplicity[1]

#**********************************
def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    os.makedirs(RESULT_PATH, exist_ok=True)

    atoms, geometry = import_xyz()
    create_for005(atoms, geometry, timestamp)
    clear_forxxx()
    run_mopac_exe()
    total_energy, gradients = parse_for006()
    create_report(total_energy, gradients, timestamp)
    return total_energy, gradients

#**********************************
def import_xyz():
    with open(MOPAC7_PATH + NAME_XYZ) as f:
        xyz = f.read()
    atoms = int(xyz.split()[0])
    geometry = np.array(xyz.split()[1:]).reshape(atoms, 4)
    return atoms, geometry

def create_for005(atoms, geometry, timestamp):
    input_text=f"""AM1 1SCF XYZ GRADIENTS DEBUG DCART GEO-OK {MULTI} CHARGE={CHARGE}
    {PROJECT_NAME} {timestamp}
    {COMMENT}\n"""

    for n in range(atoms):
        x = geometry[n]
        input_text += f"{x[0]} {x[1]} 0 {x[2]} 0 {x[3]} 0\n"

    with open(MOPAC7_PATH + "FOR005", "w") as f:
        f.write(input_text)

def wait_for_for006(timeout=120, poll=1):
    t0 = time.time()
    while not os.path.exists(MOPAC7_PATH + "FOR006"):
        if time.time() - t0 > timeout:
            raise TimeoutError("FOR 006 was not created")
        time.sleep(poll)

def clear_forxxx():
    for file in ("FOR006", "FOR009", "FOR011", "FOR012", "FOR016", "SHUTDOWN"):
        if os.path.exists(MOPAC7_PATH + file):
            os.remove(MOPAC7_PATH + file)

def run_mopac_exe():
    subprocess.run([str(MOPAC7_PATH + "mopac7.exe")], cwd=str(MOPAC7_PATH),text=True, check=False)

def parse_for006():
    wait_for_for006()

    with open(MOPAC7_PATH + "FOR006") as f:
        text = f.read()

    readFrom = text.find("TOTAL ENERGY")
    readTo = text[readFrom:].find("EV")+readFrom
    total_energy_raw = text[readFrom:readTo]
    total_energy = float(total_energy_raw.split()[-1])

    readFrom = text.find("CARTESIAN COORDINATE DERIVATIVES")
    gradients = np.array(text[readFrom:].split())[8:]
    index=np.where(gradients=="-------------------------------------------------------------------------------")[0][0]
    gradients = gradients[:index].reshape(-1, 5)

    return total_energy, gradients

def create_report(total_energy, gradients, timestamp):
    with open(MOPAC7_PATH + NAME_XYZ) as f:
        xyz = f.read()
    with open(MOPAC7_PATH + "FOR005") as f:
        for006 = f.read()
    result =f"""{PROJECT_NAME}, {timestamp}, {COMMENT} \n
    TOTAL ENERGY = {total_energy} EV\n
    CARTESIAN COORDINATE DERIVATIVES, EV/A
    {pd.DataFrame(gradients)}\n
    ***FOR006 input***
    {for006}\n
    ***  XYZ input ***
    {xyz}"""
    with open(RESULT_PATH+"result"+timestamp+".txt", "w") as f:
        f.write(result)


def sci_optimize():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")

    def get_flat_geometry():
        _, geometry = import_xyz()
        flat_geometry = geometry.reshape(-1, 4)[:, 1:].reshape(-1).astype(float)
        return flat_geometry

    def rewrite_xyz(new_geometry_flat):
        with open(MOPAC7_PATH + NAME_XYZ) as f:
            lines = [ln.strip() for ln in f.readlines()]
        n = int(lines[0])
        new_xyz = f"{n}\n \n"
        for i in range(n):
            symbol = lines[2 + i].split()
            symbol[1] = new_geometry_flat[i * 3]
            symbol[2] = new_geometry_flat[i * 3 + 1]
            symbol[3] = new_geometry_flat[i * 3 + 2]
            new_xyz += " ".join(str(x) for x in symbol) + f"\n"
        with open(MOPAC7_PATH + NAME_XYZ, "w") as f:
            f.write(new_xyz)

    def get_gradient_and_energy_from_flat_geometry(x):
        rewrite_xyz(x)
        total_energy, gradient = main()
        gradient = gradient[:, 2:]

        return float(total_energy), gradient.reshape(-1).astype(float)

    end_result = minimize(
        fun = get_gradient_and_energy_from_flat_geometry,
        x0=get_flat_geometry(),
        jac = True,
        method="L-BFGS-B",
        options = dict(
            maxiter=300,
            gtol=1e-5,
            ftol=0.0,
            maxls=50
        ))
    print(end_result)


sci_optimize()