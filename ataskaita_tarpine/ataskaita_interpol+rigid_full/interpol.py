from pathlib import Path
import numpy as np
import subprocess
from scipy.optimize import  minimize
from math import sqrt

N_IMAGES = 3 # n of intermediate images
CHARGE = 0
UNPAIRED_ELECTRONS = 10

FIXED_ATOMS_1BASED = [24,27,62,59,45,54,52,10,8,17,1,4,38,36]
FIXED_ATOMS_0BASED = [i - 1 for i in FIXED_ATOMS_1BASED]

SHIFT_LIMIT_ANG = 0.05

KCALMOL_TO_EV = 0.0433641153087705
EV_TO_HARTREE = 1 / 27.211386245988
HARTREE_PER_KCALMOL = 1 / 627.5094740631
ANGSTROM_PER_BOHR = 0.529177210903

#directories
XYZ_PATH = Path("XYZs")
MOPAC_PATH =  Path("runs")
MOPAC_EXE_PATH = Path.cwd().parents[0]/"mopac/bin/mopac.exe"
SAVES_PATH = Path("runs")

REACTANT_XYZ = Path("XYZs\\newpath_2inp_R.xyz")
PRODUCT_XYZ = Path("XYZs\\newpath_2inp_P.xyz")
NAME_XYZ = "test.xyz"

jobname = NAME_XYZ[:-4]

#--- high-level workflows ---

def prepare_initial_images():
    r_atoms, r_xyz_flat = import_xyz(REACTANT_XYZ)
    p_atoms, p_xyz_flat = import_xyz(PRODUCT_XYZ)
    assert np.array_equal(r_atoms, p_atoms)
    if check_fixed_atom_shifts(r_xyz_flat, p_xyz_flat, FIXED_ATOMS_1BASED, SHIFT_LIMIT_ANG):
        print(f"Too large >{SHIFT_LIMIT_ANG} shifts found among reactant and product fixed atoms geometries.")
        assert False
    n_atoms = len(r_atoms)
    free_atoms_mask = create_free_mask(n_atoms, FIXED_ATOMS_1BASED)
    free_xyz_mask = np.repeat(free_atoms_mask, 3)
    images = interpolate_linearly(r_xyz_flat, p_xyz_flat, N_IMAGES)

    return r_atoms, r_xyz_flat, p_xyz_flat, free_atoms_mask, free_xyz_mask, images

def create_dataset_of_images_from_xyz_files_via_1scf():
    atoms, r_xyz_flat, p_xyz_flat, free_atoms_mask, free_xyz_mask, images = prepare_initial_images()

    parameters_of_images = perform_all_image_1scf(images, atoms)
    save_images_parameters(
        path=SAVES_PATH,
        mask=free_atoms_mask,
        params = parameters_of_images,
        name="images_raw"
    )

def create_dataset_of_images_from_xyz_files_via_optimizer(): #also relaxes
    atoms, r_xyz_flat, p_xyz_flat, free_atoms_mask, free_coords_mask, images = prepare_initial_images()

    parameters_of_images = {
        "atoms": atoms,
        "coords_ang": images,
        "energies_Eh": [],
        "grads_Eh_per_bohr": [],
    }

    for image_id, image in enumerate(images):
        x0_ang = np.asarray(image,float).reshape(-1, 3)

        #1st SCFs results goes to dataset
        E0, G0 = main_mopac(x0_ang, atoms)
        parameters_of_images["energies_Eh"].append(E0)
        parameters_of_images["grads_Eh_per_bohr"].append(G0)

        #full optimization without rotors
        x_opt_ang, res = sci_minimize(x0_ang, atoms, free_coords_mask)

        x_opt_flat = x_opt_ang.ravel()
        diff_R = x_opt_flat[free_coords_mask] - r_xyz_flat[free_coords_mask]
        diff_P = x_opt_flat[free_coords_mask] - p_xyz_flat[free_coords_mask]
        dist_R = np.sqrt(np.sum(diff_R**2))
        dist_P = np.sqrt(np.sum(diff_P**2))

        if dist_R < dist_P:
            direction = "reactant"
        else:
            direction = "product"

        print(
            f"Image {image_id}: optimization converged towards {direction} "
            f"(dist to R = {dist_R:.3f} Å, dist to P = {dist_P:.3f} Å)"
        )

    save_images_parameters(SAVES_PATH, free_atoms_mask, parameters_of_images, name="images_raw_opt")

def import_xyz(path: Path):
    with open(path) as f:
        xyz = f.read().strip()
    xyz_lines = xyz.split("\n")
    atoms, geometry = [], []
    for line in range(2, len(xyz_lines)):
        c = xyz_lines[line].split()
        atoms.append(c[0])
        geometry.extend(c[1:])
    return np.array(atoms), np.array(geometry, dtype=float)

def interpolate_linearly(xyz_start, xyz_end, n_mid_img=10):
    images = [xyz_start.reshape(-1, 3)]
    for k in range(1, n_mid_img + 1):
        t = k / (n_mid_img + 1)
        geometry = (1 - t) * xyz_start + t * xyz_end
        images.append(geometry.reshape(-1, 3))
    images.append(xyz_end.reshape(-1, 3))
    return np.array(images)

def perform_all_image_1scf(images_geom, atoms):
    parameters_of_images = {
        "atoms": atoms,
        "coords_ang": images_geom,
        "energies_Eh": [],
        "grads_Eh_per_bohr": [],
    }
    for image in images_geom:
        E, G = main_mopac(image.reshape(-1, 3), atoms)
        parameters_of_images["energies_Eh"].append(E)
        parameters_of_images["grads_Eh_per_bohr"].append(G)
    return parameters_of_images

def save_images_parameters(path, mask, params, name="images"):
    path = Path(path)/name
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        atoms=np.asarray(params["atoms"]),
        coords_ang=np.asarray(params["coords_ang"], dtype=float),
        energies_Eh=np.asarray(params["energies_Eh"], dtype=float),
        grads_Eh_per_bohr=np.asarray(params["grads_Eh_per_bohr"], dtype=float),
        mask=np.asarray(mask, dtype=bool),
        units=np.array(["coords: angstrom", "energy: hartree", "grad: hartree/bohr"]),
    )
    print(f"{N_IMAGES}+2 images saved into {name}.npz")

def load_images_parameters(path):
    d = np.load(path)
    return {
        "atoms": d["atoms"],
        "coords_ang": d["coords_ang"],
        "energies_Eh": d["energies_Eh"],
        "grads_Eh_per_bohr": d["grads_Eh_per_bohr"],
        "mask": d["mask"],
        "units": d["units"],
    }
# --- MOPAC runs ---
def main_mopac(geometry, atoms):
    create_mopac_input(atoms, geometry)
    run_mopac_exe()
    return parse_mopac_output()

def create_mopac_input(atoms, geometry):
    input_text=f"""PM7 UHF 1SCF XYZ GRADIENTS AUX DCART GEO-OK MS={UNPAIRED_ELECTRONS/2} CHARGE={CHARGE}
    {NAME_XYZ}\n\n"""
    for n in range(len(atoms)):
        x = geometry[n]
        input_text += f"{atoms[n]} {x[0]} {x[1]} {x[2]}\n"
    with open(MOPAC_PATH/f"{jobname}.mop", "w") as f:
        f.write(input_text)

def run_mopac_exe():
    mopac_input = MOPAC_PATH / f"{jobname}.mop"
    subprocess.run(
        [str(MOPAC_EXE_PATH), mopac_input.name],
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

    return E_Eh, gradients_Eh_per_bohr.reshape(-1, 3)

# --- optimization ---

def sci_minimize(x0_ang, atoms, coords_mask): #with fixed atomos
    n_atoms = len(atoms)
    x0_bohr = (x0_ang / ANGSTROM_PER_BOHR).ravel()
    x0_bohr_free = x0_bohr[coords_mask]
    print(x0_bohr_free)

    def run_1scf(x_bohr_free, atoms, x0_bohr, coords_mask):

        x_bohr = x0_bohr.copy()
        x_bohr[coords_mask] = x_bohr_free

        geom_ang =  (x_bohr * ANGSTROM_PER_BOHR).reshape(-1, 3)
        E_Eh, G_Eh_per_bohr = main_mopac(geom_ang, atoms)
        grad_full = G_Eh_per_bohr.ravel()
        grad_free = grad_full[coords_mask]

        return E_Eh, grad_free

    res = minimize(
            fun = run_1scf,
            x0 = x0_bohr_free,
            method="L-BFGS-B",
            jac = True,
            args=(atoms, x0_bohr, coords_mask),
        )

    x_bohr_opt = x0_bohr.copy()
    x_bohr_opt[coords_mask] = res.x

    x_ang_opt = (x_bohr_opt * ANGSTROM_PER_BOHR).reshape(-1, 3)
    print(res)

    return x_ang_opt, res

# --- helper functions ---

def check_fixed_atom_shifts(r_xyz, p_xyz, fixed_list, shift_limit_ang):
    fixed0 = [(i-1) for i in fixed_list]
    for atom in fixed0:
        shift = sqrt(
            (r_xyz[atom*3] - p_xyz[atom*3])**2 +
            (r_xyz[atom*3+1] - p_xyz[atom*3+1])**2 +
            (r_xyz[atom*3+2] - p_xyz[atom*3+2])**2
        )
        if shift > shift_limit_ang:
            return True
    return False

def create_free_mask(n_atoms: int, fixed_atoms):
    fixed0 = [(i - 1) for i in fixed_atoms]
    mask = np.ones(n_atoms, dtype=bool)
    for fix in fixed0:
        mask[fix] = False
    return mask


#create_dataset_of_images_from_xyz_files_via_1scf() #sugeneruoja images, paskaičiuoja 1scf kiekvienam ir įrašo
#create_dataset_of_images_from_xyz_files_via_optimizer() #relaksina kiekvieną image (užtrunka)
#print(load_images_parameters("runs\\images_raw.npz"))  #arba images_raw_opt.npz
#-----------

