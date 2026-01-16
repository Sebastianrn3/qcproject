from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

#PM7 method and MS=? keyword

REACTANT_XYZ = Path("XYZinterpol\\h2o_start.xyz")
PRODUCT_XYZ = Path("XYZinterpol\\h2o_end.xyz")

n_images = 10
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
    return images


def write_xyz(path: Path, atoms, geometry, comment="generated"):
    n = len(atoms)
    text = f"{n}\n{comment}\n"
    for line in n:
        text += f"{atoms[line] + geometry[line*3 : line*3 + 2]}"


def interpolate_spline(geom_list):
    """
    Возвращает функцию geom_at(s), где s в [0, 1].
    """
    geom_list = [np.asarray(g) for g in geom_list]
    k = len(geom_list)               # сколько у нас дискретных точек
    n_atoms = geom_list[0].shape[0]

    # параметр по точкам: 0, 1, ..., k-1
    t = np.arange(k, dtype=float)

    # складываем в один большой массив (k, n_atoms, 3)
    G = np.stack(geom_list, axis=0)

    # расплющим атомы и координаты, чтобы сплайн строить одномерно
    G_flat = G.reshape(k, -1)  # (k, n_atoms*3)

    # строим сплайн
    spl = CubicSpline(t, G_flat, axis=0)

    def geom_at(s):
        """
        s: число в [0, 1]
        возвращает геометрию shape=(n_atoms, 3)
        """
        s_scaled = s * (k - 1)         # переводим в [0, k-1]
        g_flat = spl(s_scaled)         # (n_atoms*3,)
        return g_flat.reshape(n_atoms, 3)

    return geom_at


# =============================
# CONSTRAINTS (замороженные атомы)
# =============================
def split_free_frozen(coords, frozen_indices):
    """
    coords: np.array (n_atoms, 3)
    frozen_indices: список атомов, которые нельзя двигать
    Возвращает:
      free_vec  -- 1D массив координат только свободных атомов
      frozen_coords -- сами координаты замороженных атомов (та же форма, что и coords)
    """
    n_atoms = coords.shape[0]
    mask = np.ones(n_atoms, dtype=bool)
    mask[frozen_indices] = False   # False там, где заморожено

    free_coords = coords[mask]     # (n_free, 3)
    free_vec = free_coords.ravel() # (3*n_free,)
    return free_vec, coords.copy()  # второе — полный снимок, чтобы собрать обратно


def assemble_free_and_fixed(free_vec, frozen_snapshot, frozen_indices):
    """
    Обратно собирает полную геометрию:
      free_vec: 1D массив координат только свободных атомов
      frozen_snapshot: (n_atoms, 3) - полный снимок, из которого мы возьмём замороженные атомы
      frozen_indices: какие именно атомы заморожены
    Возвращает: coords_full shape=(n_atoms, 3)
    """
    coords_full = frozen_snapshot.copy()
    n_atoms = coords_full.shape[0]

    mask = np.ones(n_atoms, dtype=bool)
    mask[frozen_indices] = False

    # куда вставлять свободные
    coords_full[mask] = free_vec.reshape(-1, 3)
    return coords_full


# =============================
# ЗАГЛУШКА ДЛЯ MOPAC
# (сюда ты вставишь свои create_mopac_input/run_mopac/parse_mopac_output)
# =============================
def run_mopac_on_geom(atoms, coords):
    """
    Здесь у тебя должен быть твой рабочий пайплайн:
    - записать .mop
    - запустить MOPAC
    - распарсить .aux
    Ниже — просто имитация, чтобы скрипт был цельным.
    """
    # TODO: заменить на твои функции
    E = 0.0
    G = np.zeros_like(coords)
    return E, G


# =============================
# MAIN
# =============================
def main():
    # 1. читаем две крайности
    R_atoms, R_xyz = read_xyz(REACTANT_XYZ)
    P_atoms, P_xyz = read_xyz(PRODUCT_XYZ)

    # 2. проверяем, что порядок атомов совпадает
    assert (R_atoms == P_atoms).all(), "Reactant and product must have the same atom order"

    # 3. строим линейные промежуточные
    intermediates = linear_interpolate_geoms(R_xyz, P_xyz, N_INTERMEDIATE)

    # 4. собираем полный список геометрий: R, ...intermediates..., P
    all_geoms = [R_xyz] + intermediates + [P_xyz]

    # 5. строим сплайн поверх этих точек (он нам понадобится позже)
    spline_geom = build_spline_from_geoms(all_geoms)

    # 6. показываем, как можно сохранить все структуры в папку
    out_dir = Path("generated_path")
    out_dir.mkdir(exist_ok=True)
    for i, geom in enumerate(all_geoms):
        write_xyz(out_dir / f"image_{i:02d}.xyz", R_atoms, geom,
                  comment=f"image {i}")

    # 7. для каждой структуры — считаем энергию и градиент (пока заглушка)
    energies = []
    for i, geom in enumerate(all_geoms):
        E, G = run_mopac_on_geom(R_atoms, geom)
        energies.append(E)
        print(f"Image {i}: E = {E:.6f} (mock)")

        # 8. пример: как вырезать свободные атомы из этой геометрии
        free_vec, frozen_snapshot = split_free_frozen(geom, FROZEN_ATOMS)
        # ... здесь бы мы вызвали scipy.optimize.minimize(f, free_vec, args=(frozen_snapshot, FROZEN_ATOMS)) ...
        # и внутри f мы бы снова собрали полную геометрию:
        # full = assemble_free_frozen(free_vec, frozen_snapshot, FROZEN_ATOMS)
        # и уже её отдали бы в MOPAC

    # 9. пример: как взять точку со сплайна между нашими дискретными
    mid_geom = spline_geom(0.37)
    write_xyz(out_dir / "image_spline_037.xyz", R_atoms, mid_geom,
              comment="spline at s=0.37")

    print("Done.")


if __name__ == "__main__":
    main()

