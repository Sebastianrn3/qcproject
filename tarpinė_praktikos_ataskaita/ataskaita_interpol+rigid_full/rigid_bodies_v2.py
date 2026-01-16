#params
FIXED_ATOMS_1BASED = [1,4,8,10,17,24,27,36,38,45,52,54,59,62]
ANCHOR_IDS_1BASED = [1, 8, 17, 24, 36, 45, 52, 59, 73]
RIGID_GROUPS_1BASED = [
    list(range(1, 8)), # acetate_1 = [1:7]
    list(range(8, 17)),  # asn_residue_1 = [8:16]
    list(range(17, 24)), # acetate_2 = [17:23]
    list(range(24, 36)),  # his_residue_1 = [24:35]
    list(range(36, 45)),  # asn_residue_2 = [36:44]
    list(range(45, 52)), # acetate_3 = [45:51]
    list(range(52, 59)), # acetate_4 = [52:58]
    list(range(59, 71)), # his_residue_2 = [59:70]
    list(range(73, 76)), # water = [73:75]
]

FIXED_ATOMS_0BASED = [i - 1 for i in FIXED_ATOMS_1BASED]
RIGID_GROUPS_0BASED = [[i - 1 for i in g] for g in RIGID_GROUPS_1BASED]
ANCHOR_IDS_0BASED = [i - 1 for i in ANCHOR_IDS_1BASED]


# ---rigid-rotor-related ---

def prepare_rigid_rotors(x0_bohr_1d, rigid_groups, anchor_ids):
    assert len(rigid_groups) == len(anchor_ids)

    x0_bohr = np.asarray(x0_bohr_1d).reshape(-1, 3)
    groups_data = []

    for atoms, anchor_id in zip(rigid_groups, anchor_ids):
        anchor0 = x0_bohr[anchor_id].copy()
        rel0 = x0_bohr[atoms] - anchor0

        groups_data.append({
            "atoms": atoms,
            "anchor_id": anchor_id,
            "anchor0": anchor0,
            "rel0": rel0,
        })
    return groups_data


def build_geom_from_rotors(x, x0_bohr_1d, rotor_data, free_atoms):
    x0_bohr = np.asarray(x0_bohr_1d).reshape(-1, 3)
    full_geom= x0_bohr.copy()

    offset = 0
    for group in rotor_data: #rigids
        w = np.asarray(x[offset:offset+3])
        offset += 3
        R = rodrigues_from_vector(w)

        xyz_group = group["anchor0"] + group["rel0"] @ R.T
        full_geom[group["atoms"]] = xyz_group

    xyz_free_atoms = np.asarray(x[offset:]).reshape(-1, 3)
    full_geom[free_atoms] = xyz_free_atoms

    return full_geom.ravel()

def rodrigues_from_vector(w):
    w = np.asarray(w, float)
    teta = np.linalg.norm(w) #teta = |w|
    if teta < 1e-12:
        return np.eye(3)
    k = w / teta
    kx, ky, kz = k
    K = np.array([
        [0, -kz,  ky],
        [kz,  0, -kx],
        [-ky, kx,  0],
    ]) #Kv = k*v
    I = np.eye(3)
    R = I + np.sin(teta) * K + (1 - np.cos(teta)) * (K @ K)
    return R

def pack_gradients_multi_rotor(
        G_full,
        full_bohr,
        rotor_data,
        free_atoms,
):
    G_full = np.asarray(G_full).reshape(-1, 3)
    coords = np.asarray(full_bohr).reshape(-1, 3)

    G_parts = []

    for group in rotor_data: #rigids
        group_atoms = group["atoms"]
        anchor_id = group["anchor_id"]

        anchor0 = coords[anchor_id]
        r_rel = coords[group_atoms] - anchor0

        torque_g = np.sum(np.cross(r_rel, G_full[group_atoms]), axis=0)
        G_parts.append(torque_g)

    G_parts.append(G_full[free_atoms].ravel())

    return np.concatenate(G_parts)

# functions for testing

def create_dataset_of_images_from_xyz_files_via_optimizer():
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
        #x_opt_ang, res = sci_minimize(x0_ang, atoms, free_coords_mask)

        x_opt_ang, res = sci_minimize_multi_rotors(
            x0_ang,
            atoms,
            fixed_atoms=FIXED_ATOMS_0BASED,
            anchor_ids=ANCHOR_IDS_0BASED,
            rigid_groups=RIGID_GROUPS_0BASED,
        )

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

    save_images_parameters(SAVES_PATH, free_atoms_mask, parameters_of_images, name="images_raw_opt_rigid")

def sci_minimize_multi_rotors( #with fixed and rigid bodies
        x0_ang,
        atoms,
        fixed_atoms=None,
        anchor_ids=None,
        rigid_groups=None,
):
    fixed_atoms = [] if fixed_atoms is None else fixed_atoms
    anchor_ids = [] if anchor_ids is None else anchor_ids
    rigid_groups = [] if rigid_groups is None else rigid_groups
    n_atoms = len(atoms)
    x0_bohr = (x0_ang/ANGSTROM_PER_BOHR).ravel()
    fixed_atoms = np.asarray(fixed_atoms)
    anchor_ids = np.asarray(anchor_ids)

    #sorting rotor's, non-anchor fixed and free atoms:
    mask_rotor = np.zeros(n_atoms, dtype=bool)
    for group in rigid_groups:
        group = np.asarray(group, dtype=int)
        mask_rotor[group] = True

    mask_fixed = np.zeros(n_atoms, dtype=bool)
    if fixed_atoms.size:
        mask_fixed[fixed_atoms] = True

    if anchor_ids.size:
        mask_fixed[anchor_ids] = False

    free_atoms = np.arange(n_atoms, dtype=int)[~mask_rotor & ~mask_fixed]

    rotors_data = prepare_rigid_rotors(x0_bohr, rigid_groups, anchor_ids)

    x0_free_coords = x0_bohr.reshape(-1, 3)[free_atoms].ravel()
    x0 = np.concatenate([
        np.zeros(3 * len(rotors_data)),
        x0_free_coords
    ])

    def run_1scf(x):
        full_bohr = build_geom_from_rotors(x, x0_bohr, rotors_data, free_atoms)
        geom_ang = (full_bohr * ANGSTROM_PER_BOHR).reshape(-1, 3)

        E_Eh, G_Eh_per_bohr = main_mopac(geom_ang, atoms)

        G_with_packed_rotors = pack_gradients_multi_rotor(
            G_Eh_per_bohr,
            full_bohr,
            rotors_data,
            free_atoms
        )
        return E_Eh, G_with_packed_rotors

    res = minimize(
        fun=run_1scf,
        x0=x0,
        method="L-BFGS-B",
        jac=True,
    )

    full_bohr_opt = build_geom_from_rotors(res.x, x0_bohr,rotors_data, free_atoms)
    full_ang_opt = (full_bohr_opt * ANGSTROM_PER_BOHR).reshape(-1, 3)

    return full_ang_opt, res

create_dataset_of_images_from_xyz_files_via_optimizer()
print(load_images_parameters("runs\\images_raw_opt_rigid.npz"))