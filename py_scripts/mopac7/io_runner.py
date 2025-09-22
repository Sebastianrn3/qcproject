import numpy as np
import pandas as pd
import subprocess, os
from datetime import datetime

KCAL_TO_KJ = 4.184
EV_TO_KJMOL = 96.4853321233
EV_TO_KCALMOL = 23.060548


MOPAC7_PATH = "C:/Users/ddizy/Desktop/mopac7/"
RESULT_PATH = "C:/Users/ddizy/Desktop/qcproject/py_scripts/mopac7/mopac7_results/"
NAME_XYZ = "ch3cl_for_grad_test.xyz"


multiplicity =  ["INVALID_MULTIPLICITY", "SINGLET", "DOUBLET", "TRIPLET", "QUARTET", "QUINTET", "SEXTET"]
PROJECT_NAME="CH3Cl_gas"
COMMENT="DFT_Energy_Gradient"
CHARGE=0
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
    print("DONE")

#**********************************
def import_xyz():
    with open(MOPAC7_PATH + NAME_XYZ) as f:
        xyz = f.read()
    atoms = int(xyz.split()[0])
    geometry = np.array(xyz.split()[1:]).reshape(atoms, 4)
    return atoms, geometry

def create_for005(atoms, geometry, timestamp):
    input_text=f"""AM1 1SCF XYZ GRADIENTS DEBUG DCART {MULTI} CHARGE={CHARGE}
    {PROJECT_NAME} {timestamp}
    {COMMENT}\n"""

    for n in range(atoms):
        x = geometry[n]
        input_text += f"{x[0]} {x[1]} 0 {x[2]} 0 {x[3]} 0\n"

    with open(MOPAC7_PATH + "FOR005", "w") as f:
        f.write(input_text)

def clear_forxxx():
    for file in ("FOR006", "FOR009", "FOR011", "FOR012", "FOR016", "SHUTDOWN"):
        if os.path.exists(MOPAC7_PATH + file):
            os.remove(MOPAC7_PATH + file)

def run_mopac_exe():
    subprocess.run([str(MOPAC7_PATH + "mopac7.exe")], cwd=str(MOPAC7_PATH),text=True, check=False)

def parse_for006():
    with open(MOPAC7_PATH + "FOR006") as f:
        text = f.read()

    readFrom = text.find("TOTAL ENERGY")
    readTo = text[readFrom:].find("EV")+readFrom
    total_energy_raw = text[readFrom:readTo]
    total_energy = float(total_energy_raw.split()[-1])*EV_TO_KJMOL

    readFrom = text.find("CARTESIAN COORDINATE DERIVATIVES")
    gradients = np.array(text[readFrom:].split())[8:]
    index=np.where(gradients=="-------------------------------------------------------------------------------")[0][0]
    gradients = gradients[:index].reshape(-1, 5)

    for row in gradients:
        row[2] = float(row[2]) * KCAL_TO_KJ
        row[3] = float(row[3]) * KCAL_TO_KJ
        row[4] = float(row[4]) * KCAL_TO_KJ

    return total_energy, gradients



def create_report(total_energy, gradients, timestamp):
    with open(MOPAC7_PATH + NAME_XYZ) as f:
        xyz = f.read()
    with open(MOPAC7_PATH + "FOR005") as f:
        for006 = f.read()
    result =f"""{PROJECT_NAME}, {timestamp}, {COMMENT} \n
    TOTAL ENERGY = {total_energy} KJ/MOL\n
    CARTESIAN COORDINATE DERIVATIVES, kJ/mol/angstrem
    {pd.DataFrame(gradients)}\n
    ***FOR006 input***
    {for006}\n
    ***  XYZ input ***
    {xyz}"""
    with open(RESULT_PATH+"result"+timestamp+".txt", "w") as f:
        f.write(result)

main()