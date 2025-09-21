import numpy as np
import subprocess, os

EV_TO_KJMOL = 96.4853321233
EV_TO_KCALMOL = 23.060548
#**********************************
#importuojama iš .xyz - atomu sk ir geometriju masyva
xyz=""
with open("C:/Users/User/Desktop/mopac7/molecule.xyz") as f:
  xyz = f.read()
molecule = xyz.split()
atoms = int(xyz.split()[0])
geometry = np.array(molecule[1:]).reshape(atoms, 4)

#kuriamas input_text - FOR005
input_text="AM1 1SCF XYZ GRADIENTS\nDarbas\nStartuoja\n"

for n in range(atoms):
    x = geometry[n]
    new_line = x[0] + " " + x[1] + " 0 " + x[2] + " 0 " + x[3] + " 0\n"
    input_text += new_line
with open("C:/Users/User/Desktop/mopac7/FOR005", "w") as f:
  f.write(input_text)

#pasalinamas FOR006 ir kt pries atliekant nauja skaiciavima
if os.path.exists("C:/Users/User/Desktop/mopac7/FOR006"):
    os.remove("C:/Users/User/Desktop/mopac7/FOR006")
if os.path.exists("C:/Users/User/Desktop/mopac7/FOR011"):
    os.remove("C:/Users/User/Desktop/mopac7/FOR011")
    os.remove("C:/Users/User/Desktop/mopac7/FOR012")
    os.remove("C:/Users/User/Desktop/mopac7/FOR016")

#leidziamas exe
exe = "C:/Users/User/Desktop/mopac7/mopac7.exe"
subprocess.run([str(exe)], cwd=str("C:/Users/User/Desktop/mopac7"),text=True, check=False)

#parsinam
text=""
with open("C:/Users/User/Desktop/mopac7/FOR006") as f:
  text = f.read()

readFrom = text.find("TOTAL ENERGY")
readTo = text[readFrom:].find("EV")+readFrom
total_energy = float(text[readFrom:readTo].split()[-1])*EV_TO_KJMOL

readFrom = text.find("FINAL  POINT  AND  DERIVATIVES")
readTo   = text.rfind("KCAL/ANGSTROM")
result = "TOTAL ENERGY = "+ str(total_energy)+" KCAL/MOL\n\n" + text[readFrom:readTo+15]

with open("C:/Users/User/Desktop/mopac7/result.txt", "w") as f:
  f.write(result)
