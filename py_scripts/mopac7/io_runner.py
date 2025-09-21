import numpy as np
import subprocess, os
from datetime import datetime

EV_TO_KJMOL = 96.4853321233
EV_TO_KCALMOL = 23.060548

project_name="TEST"
MOPAC7_PATH = "C:/Users/ddizy/Desktop/mopac7/"
RESULT_PATH = "C:/Users/ddizy/Desktop/qcproject/tests/mopac7_results/"
NAME_XYZ = "molecule.xyz"

multiplicity =  ["INVALID_MULTIPLICITY", "SINGLET", "DOUBLET", "TRIPLET", "QUARTET", "QUINTET", "SEXTET"]
start_time = datetime.now().strftime("%Y%m%d_%H%M%S_%f")

COMMENT="TEST"
CHARGE=-1
MULTI=multiplicity[2]
SPIN=0.5
#**********************************
#importuojama iš molecule.xyz
xyz=""
with open(MOPAC7_PATH + NAME_XYZ) as f:
  xyz = f.read()

molecule = xyz.split()
atoms = int(xyz.split()[0])
geometry = np.array(molecule[1:]).reshape(atoms, 4)

#kuriamas FOR005
input_text=f"""AM1 1SCF XYZ GRADIENTS {MULTI} MS={SPIN} CHARGE={CHARGE}
{project_name} {start_time}
{COMMENT}"""

for n in range(atoms):
    x = geometry[n]
    input_text += f"{x[0]} {x[1]} 0 {x[2]} 0 {x[3]} 0\n"

with open(MOPAC7_PATH + "FOR005", "w") as f:
  f.write(input_text)

#pasalinamas FOR006 ir kt pries atliekant nauja skaiciavima
if os.path.exists(MOPAC7_PATH + "FOR006"):
    os.remove(MOPAC7_PATH + "FOR006")
if os.path.exists(MOPAC7_PATH + "FOR011"):
    os.remove(MOPAC7_PATH + "FOR011")
if os.path.exists(MOPAC7_PATH + "FOR012"):
    os.remove(MOPAC7_PATH + "FOR012")
if os.path.exists(MOPAC7_PATH + "FOR016"):
    os.remove(MOPAC7_PATH + "FOR016")

#leidziamas exe
exe = MOPAC7_PATH + "mopac7.exe"
subprocess.run([str(exe)], cwd=str(MOPAC7_PATH),text=True, check=False)

#parsinam
text=""
with open(MOPAC7_PATH + "FOR006") as f:
  text = f.read()

readFrom = text.find("TOTAL ENERGY")
readTo = text[readFrom:].find("EV")+readFrom

total_energy = float(text[readFrom:readTo].split()[-1])*EV_TO_KCALMOL

readFrom = text.find("FINAL  POINT  AND  DERIVATIVES")
readTo   = text.rfind("KCAL/ANGSTROM")

result =f"""TOTAL ENERGY = {total_energy} KCAL/MOL\n\n {text[readFrom:readTo+15]}
***FOR005 input***
{input_text}
***  XYZ input ***
{xyz}"""

with open(RESULT_PATH + "result"+start_time+".txt", "w") as f:
  f.write(result)
