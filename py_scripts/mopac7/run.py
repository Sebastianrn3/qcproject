import os
import subprocess
import numpy as np
from datetime import datetime
from scipy.optimize import minimize

# ===== Константы единиц =====
KCAL_TO_KJ = 4.184
EV_TO_KJMOL = 96.4853321233
EV_TO_KCALMOL = 23.060548

# ===== Пути и настройки =====
MOPAC7_PATH = r"C:\Users\ddizy\Desktop\mopac7"
RESULT_PATH = r"C:\Users\ddizy\Desktop\qcproject\py_scripts\mopac7\mopac7_results"
NAME_XYZ    = "ch3cl_for_grad_test.xyz"

PROJECT_NAME = "CH3Cl_gas"
COMMENT      = "DFT_Energy_Gradient"
CHARGE       = 0
MULTI_NAME   = "SINGLET"

def ensure_dirs():
    os.makedirs(RESULT_PATH, exist_ok=True)

def read_xyz():
    with open(r"C:\Users\ddizy\Desktop\mopac7\ch3cl_for_grad_test.xyz", "r") as f:
        lines = f.read()
    atoms = ['C', 'H', 'H', 'H', 'Cl']
    arr = np.array(lines.split()[1:])

    coords = arr.reshape(int(lines[0]), 4)[:, 1:].astype(float)
    print(atoms, coords)

    return atoms, np.array(coords, float), (lines[1] if len(lines) > 1 else "")

def mopac_mult_keywords(mult_name: str):
    """Ключевые слова спина/мультиплетности для MOPAC7."""
    mult_name = mult_name.upper()
    if mult_name in ("", "SINGLET"):
        return []                      # синглет по умолчанию
    # Для открытых оболочек нужен UHF
    allowed = {"DOUBLET","TRIPLET","QUARTET","QUINTET","SEXTET"}
    if mult_name not in allowed:
        raise ValueError(f"Неизвестная мультиплетность: {mult_name}")
    return ["UHF", mult_name]

def build_for005(method, atoms, coords, charge, mult_name, title, comment, want_grad=True):
    """
    Формируем FOR005 в XYZ-формате MOPAC7.
    Флаги '1' после координат — разрешаем изменять координаты (нужно и для градиента).
    """
    kw = [method, "1SCF", "XYZ", "DCART"]   # DCART печатает декартовы градиенты
    if want_grad:
        kw.append("GRADIENTS")
    if charge != 0:
        kw.append(f"CHARGE={charge}")
    kw += mopac_mult_keywords(mult_name)

    head = " ".join(kw) + f"\n{title}\n{comment}\n"
    lines = []
    for el, (x, y, z) in zip(atoms, coords):
        lines.append(f"{el:<2}  {x:.6f} 1  {y:.6f} 1  {z:.6f} 1")
    return head + "\n".join(lines) + "\n"

def clear_forxxx(folder):
    for name in ("FOR006","FOR009","FOR011","FOR012","FOR016","SHUTDOWN"):
        p = os.path.join(folder, name)
        if os.path.exists(p):
            try: os.remove(p)
            except: pass

def run_mopac_exe(folder):
    exe = os.path.join(folder, "mopac7.exe")
    proc = subprocess.run([exe], cwd=folder, text=True, capture_output=True, check=False)
    outp = os.path.join(folder, "FOR006")
    if not os.path.exists(outp):
        raise RuntimeError("MOPAC7 не создал FOR006. Консоль:\n" + (proc.stdout or ""))
    with open(outp, "r", encoding="utf-8", errors="replace") as f:
        return f.read()

def parse_total_energy_kJ(for006_text: str) -> float:
    """
    Ищем TOTAL ENERGY в EV или KCAL/MOL и конвертируем в kJ/mol.
    """
    i = for006_text.find("TOTAL ENERGY")
    if i != -1:
        tail = for006_text[i:i+200]
        j_ev = tail.find("EV")
        if j_ev != -1:
            val_ev = float(tail[:j_ev].split()[-1])
            return val_ev * EV_TO_KJMOL
        j_kcal = tail.find("KCAL/MOL")
        if j_kcal != -1:
            val_kcal = float(tail[:j_kcal].split()[-1])
            return val_kcal * KCAL_TO_KJ
    # fallback: Heat of formation
    k = for006_text.find("TOTAL ENERGY")
    if k != -1:
        tail = for006_text[k:k+120]
        end = tail.find("KCAL/MOL")
        if end != -1:
            val_kcal = float(tail[:end].split()[-1])
            return val_kcal * KCAL_TO_KJ
    raise ValueError("Энергия не найдена в FOR006.")

def parse_gradients_kJ_per_A(for006_text: str, n_atoms: int) -> np.ndarray:
    """
    Парсим блок 'CARTESIAN COORDINATE DERIVATIVES' (kcal/mol/Å) → переводим в kJ/mol/Å.
    Возвращаем массив (N,3).
    """
    key = "CARTESIAN COORDINATE DERIVATIVES"
    a = for006_text.find(key)
    if a == -1:
        raise ValueError("Градиенты не найдены (нет блока CARTESIAN COORDINATE DERIVATIVES).")
    # Отрежем блок до следующей линии дефисов (----)
    tail = for006_text[a:].splitlines()
    # Найдём первую подчеркивающую линию после заголовков
    rows = []
    # пропустим первые ~5–8 строк, затем читаем до линии дефисов
    start = 0
    for idx, ln in enumerate(tail[:15]):
        if set(ln.strip()) == {"-"}:
            start = idx + 1
            break
    for ln in tail[start:]:
        if set(ln.strip()) == {"-"} or not ln.strip():
            break
        toks = ln.split()
        # ожидаем: idx, sym, gx, gy, gz  (или idx, sym, x, y, z, gx, gy, gz — варьирует)
        # безопасно возьмём последние 3 токена как градиент:
        try:
            gx, gy, gz = map(float, toks[-3:])
            rows.append([gx * KCAL_TO_KJ, gy * KCAL_TO_KJ, gz * KCAL_TO_KJ])
        except:
            pass
    if len(rows) < n_atoms:
        # иногда печатается другой блок "GRADIENTS" — можно попробовать его как запасной
        raise ValueError(f"Считано {len(rows)} строк градиента, ожидал {n_atoms}.")
    return np.array(rows[:n_atoms], float)

def write_report(total_E_kJ, G_kJ_per_A, for005_text, xyz_text, timestamp):
    out_path = os.path.join(RESULT_PATH, f"result_{timestamp}.txt")
    lines = []
    lines.append(f"{PROJECT_NAME}, {timestamp}, {COMMENT}")
    lines.append(f"TOTAL ENERGY = {total_E_kJ:.6f} kJ/mol")
    lines.append("")
    lines.append("CARTESIAN COORDINATE DERIVATIVES (kJ/mol/Å):")
    for i, (gx,gy,gz) in enumerate(G_kJ_per_A, 1):
        lines.append(f"{i:3d}  {gx: .6f}  {gy: .6f}  {gz: .6f}")
    lines.append("\n*** FOR005 input ***")
    lines.append(for005_text)
    lines.append("\n*** XYZ input ***")
    lines.append(xyz_text)
    with open(out_path, "w", encoding="utf-8", errors="replace") as f:
        f.write("\n".join(lines))
    print("Report:", out_path)

# ===================== ВЫЗОВ ОДНОГО РАСЧЁТА =====================

def singlepoint_energy_grad(atoms, coords_A, method="AM1", charge=0, mult_name="SINGLET"):
    """
    Собирает FOR005, запускает MOPAC7, возвращает (E_kJ/mol, G(N,3) kJ/mol/Å, for005_text, xyz_text).
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    title = f"{PROJECT_NAME} {timestamp}"

    # Текст FOR005
    for005 = build_for005(method, atoms, coords_A, charge, mult_name, title, COMMENT, want_grad=True)

    # Запись и запуск
    clear_forxxx(MOPAC7_PATH)
    with open(os.path.join(MOPAC7_PATH, "FOR005"), "w", encoding="utf-8") as f:
        f.write(for005)
    out_text = run_mopac_exe(MOPAC7_PATH)

    # Парсинг
    E_kJ  = parse_total_energy_kJ(out_text)
    G_kJ  = parse_gradients_kJ_per_A(out_text, n_atoms=len(atoms))

    # Восстановим XYZ как текст (для отчёта)
    xyz_lines = [str(len(atoms)), COMMENT]
    xyz_lines += [f"{el} {x:.6f} {y:.6f} {z:.6f}" for el,(x,y,z) in zip(atoms, coords_A)]
    xyz_text = "\n".join(xyz_lines)
    return E_kJ, G_kJ, for005, xyz_text, timestamp

# ===================== ОПТИМИЗАЦИЯ SCiPy =====================

def make_objective(atoms, charge, mult_name, method="AM1"):
    """
    Возвращает пару функций (fun, jac) для SciPy.
    Чтобы не считать два раза на одной точке, используем примитивный кэш.
    """
    N = len(atoms)
    cache = {"x": None, "E": None, "G": None}

    def evaluate(x):
        # простейший кэш — если то же x, возвращаем без пересчёта
        if cache["x"] is not None and np.allclose(x, cache["x"]):
            return cache["E"], cache["G"]
        xyz = x.reshape(N, 3)
        E, G, _, _, _ = singlepoint_energy_grad(atoms, xyz, method=method, charge=charge, mult_name=mult_name)
        cache["x"], cache["E"], cache["G"] = x.copy(), E, G.reshape(-1)
        return cache["E"], cache["G"]

    def fun(x):
        E, _ = evaluate(x)
        return E

    def jac(x):
        _, g = evaluate(x)
        return g

    return fun, jac

def optimize_geometry(xyz_atoms, xyz_coords, charge=0, mult_name="SINGLET", method="AM1",
                      gtol=1e-3, maxiter=200):
    fun, jac = make_objective(xyz_atoms, charge, mult_name, method=method)
    x0 = np.asarray(xyz_coords, float).reshape(-1)
    res = minimize(fun=fun, x0=x0, jac=jac, method="L-BFGS-B",
                   options=dict(gtol=gtol, maxiter=maxiter))
    return res, res.x.reshape(len(xyz_atoms), 3)

# ===================== MAIN =====================

def main():
    ensure_dirs()
    atoms, coords, _ = read_xyz()

    # --- одиночный расчёт энергии и градиента ---
    E_kJ, G_kJ, for005, xyz_text, ts = singlepoint_energy_grad(
        atoms, coords, method="AM1", charge=CHARGE, mult_name=MULTI_NAME
    )
    write_report(E_kJ, G_kJ, for005, xyz_text, ts)

    # --- оптимизация геометрии через SciPy ---
    print("Запускаю SciPy.optimize (L-BFGS-B)...")
    res, opt_coords = optimize_geometry(
        atoms, coords, charge=CHARGE, mult_name=MULTI_NAME, method="AM1",
        gtol=1e-3, maxiter=100
    )
    print("Converged:", res.success, res.message)
    print("Final E (kJ/mol):", res.fun)
    print("Optimized geometry (Å):")
    for el, (x,y,z) in zip(atoms, opt_coords):
        print(f"{el:2}  {x: .6f}  {y: .6f}  {z: .6f}")

if __name__ == "__main__":
    main()

def make_objective(atoms):
    N = 5
    def fun(x):
        xyz = x.reshape(N, 3)
        E, _ = energy, gradients = parse_for006()
        return E
    def jac(x):
        xyz = x.reshape(N, 3)
        _, G = parse_for006()
        return G.reshape(-1)
    return fun, jac

def optimize_with_scipy(xyz_coords,
                        charge=0,
                        mult="SINGLET",
                        xyz_atoms=['C', 'H', 'H', 'H', 'Cl'],
                        method="AM1",
                        gtol=1e-3, maxiter=200):

    fun, jac = make_objective(xyz_atoms)

    x0 = np.asarray(xyz_coords, float).reshape(-1)

    res = minimize(fun=fun, x0=x0, jac=jac, method="L-BFGS-B",
                   options=dict(gtol=gtol, maxiter=maxiter))
    return res, res.x.reshape(len(xyz_atoms), 3)



ats = optimize_with_scipy()
print(ats)

