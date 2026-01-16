from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

#PM7 method and MS=? keyword

REACTANT_XYZ = Path("XYZinterpol\\h2o_start.xyz")
PRODUCT_XYZ = Path("XYZinterpol\\h2o_end.xyz")

n_images = 3
fixed_atoms = [24,27,62,59,45,54,52,10,8,17,1,4,38,36]
total_charge = 0
unpaired_electrons = 10



def import_xyz(path: Path):
    with open(path) as f:
        xyz = f.read().strip()

    xyz_lines = xyz.split("\n")

    n_atoms = int(xyz_lines[0])
    atoms = []
    geometry = []

    for line in range(2, len(xyz_lines)):
        c = xyz_lines[line].split()
        atoms.append(c[0])
        geometry.extend(c[1:])

    return np.array(atoms), np.array(geometry, dtype=float)


def interpolate_linearly(xyz_start, xyz_end, n_img=n_images):
    images = []

    for k in range(1, n_img + 1):
        t = k / (n_img + 1)
        geometry = (1 - t) * xyz_start + t * xyz_end
        images.append(geometry)
    return np.array(images)





_, start = import_xyz(REACTANT_XYZ)
_, end = import_xyz(PRODUCT_XYZ)
r = interpolate_linearly(start, end)

print(r)