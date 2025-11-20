import os
import sys
import re
import subprocess
import numpy as np
import shutil

# ==============================================================================
# 1. БАЗА ДАННЫХ ГАЗОВ
# ==============================================================================
GAS_DB = {
    "CO2": {
        "atoms": [{"eps": 0.05365, "sig": 2.80}, {"eps": 0.1569, "sig": 3.05}],
        "rigid": True,
        "content": """# CO2 TraPPE
3 atoms
2 bonds
1 angles
Coords
1 0.0 0.0 0.0
2 1.16 0.0 0.0
3 -1.16 0.0 0.0
Types
1 1 # C
2 2 # O
3 2 # O
Charges
1 0.70
2 -0.35
3 -0.35
Bonds
1 1 1 2
2 1 1 3
Angles
1 1 2 1 3
"""
    },
    "N2": {
        "atoms": [{"eps": 0.0, "sig": 0.0}, {"eps": 0.0715, "sig": 3.31}],
        "rigid": True,
        "content": """# N2 TraPPE (3-site)
3 atoms
2 bonds
1 angles
Coords
1 0.0 0.0 0.0
2 0.0 0.0 0.55
3 0.0 0.0 -0.55
Types
1 1 # COM
2 2 # N
3 2 # N
Charges
1 0.964
2 -0.482
3 -0.482
Bonds
1 1 1 2
2 1 1 3
Angles
1 1 2 1 3
"""
    },
    "CH4": {
        "atoms": [{"eps": 0.294, "sig": 3.73}],
        "rigid": False,
        "content": """# CH4 TraPPE-UA
1 atoms
0 bonds
0 angles
Coords
1 0.0 0.0 0.0
Types
1 1 # CH4
Charges
1 0.0
"""
    },
    "H2S": {
        "atoms": [{"eps": 0.208, "sig": 3.72}, {"eps": 0.0, "sig": 0.0}],
        "rigid": True,
        "content": """# H2S TraPPE
4 atoms
3 bonds
2 angles
Coords
1 0.0 0.0 0.0
2 0.937 0.966 0.0
3 0.937 -0.966 0.0
4 0.4 0.0 0.0
Types
1 1 # S
2 2 # H
3 2 # H
4 3 # M
Charges
1 0.0
2 0.25
3 0.25
4 -0.5
Bonds
1 1 1 2
2 1 1 3
3 1 1 4
Angles
1 2 1 3
2 2 1 4
"""
    },
    "SO2": {
        "atoms": [{"eps": 0.168, "sig": 3.58}, {"eps": 0.078, "sig": 3.05}],
        "rigid": True,
        "content": """# SO2 TraPPE
3 atoms
2 bonds
1 angles
Coords
1 0.0 0.0 0.0
2 1.43 0.0 0.0
3 -0.71 1.24 0.0
Types
1 1 # S
2 2 # O
3 2 # O
Charges
1 1.14
2 -0.57
3 -0.57
Bonds
1 1 1 2
2 1 1 3
Angles
1 2 1 3
"""
    }
}

# ==============================================================================
# 2. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# ==============================================================================
def parse_charges_from_cif(cif_path):
    charges = []
    try:
        with open(cif_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError: return []

    loop_started = False
    charge_idx = -1
    headers = []
    data_start_idx = 0
    
    for i, line in enumerate(lines):
        clean = line.strip()
        if clean.startswith("loop_"):
            loop_started = True
            headers = []
            continue
        if clean.startswith("_atom_site"):
            if loop_started: headers.append(clean)
        if loop_started and clean.startswith("_atom_site"): pass
        elif loop_started and not clean.startswith("_atom_site") and not clean.startswith("#") and len(clean) > 0:
            for h_i, h in enumerate(headers):
                if "_atom_site_charge" in h:
                    charge_idx = h_i; break
            data_start_idx = i; break
    
    if charge_idx == -1: return []
    for i in range(data_start_idx, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith("loop_") or line.startswith("#"): break
        parts = line.split()
        if len(parts) > charge_idx:
            try: charges.append(float(parts[charge_idx]))
            except ValueError: pass
    return charges

# ==============================================================================
# 3. МЕНЕДЖЕР СИМУЛЯЦИИ
# ==============================================================================
class SimulationManager:
    def __init__(self, cif_file, gas_name, T=298.0, ff="UFF", cutoff=12.5):
        self.cif = cif_file
        self.gas = gas_name
        self.T = T
        self.ff = ff
        self.cutoff = cutoff
        
        if gas_name not in GAS_DB:
            raise ValueError(f"Газ {gas_name} не найден.")
        self.props = GAS_DB[gas_name]
        base = os.path.splitext(os.path.basename(cif_file))[0]
        self.data_raw = f"data.{base}"
        self.data_final = f"data.{base}_final"
        self.gas_file = f"{gas_name.lower()}.mol"
        
        # Стили по умолчанию (будут перезаписаны автодетектором)
        self.styles = {
            "bond": "harmonic",
            "angle": "cosine/periodic", # Часто используется в UFF
            "dihedral": "harmonic",
            "improper": "fourier"       # Часто используется в UFF
        }

    def step_1_topology(self):
        print(f"\n[1/5] Топология MOF (lammps-interface)...")
        cmd = [
            "lammps-interface", self.cif,
            "--force_field", self.ff,
            "--cutoff", str(self.cutoff),
            "--replication", "2x2x2",
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                if "--output" in res.stderr:
                    cmd = cmd[:-2]
                    subprocess.run(cmd, check=True)
                    default = f"data.{os.path.basename(self.cif)}"
                    if os.path.exists(default): os.rename(default, self.data_raw)
                else:
                    print(f"ОШИБКА: {res.stderr}"); sys.exit(1)
        except FileNotFoundError:
            print("ОШИБКА: lammps-interface не найден."); sys.exit(1)

    def step_2_charges(self):
        print(f"[2/5] Расчет зарядов MOF (PyEQEQ)...")
        cmd = ["eqeq", self.cif]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            base = os.path.splitext(self.cif)[0]
            cands = [f for f in os.listdir('.') if f.startswith(base) and "eqeq" in f and f.endswith(".cif")]
            
            if not cands:
                print("   [WARN] Заряды не найдены. Исп. 0.0")
                shutil.copy(self.data_raw, self.data_final); return
            
            charges = parse_charges_from_cif(cands[0])
            if not charges:
                print("   [WARN] Ошибка чтения зарядов. Исп. 0.0")
                shutil.copy(self.data_raw, self.data_final); return

            with open(self.data_raw, 'r') as f: lines = f.readlines()
            num_atoms = 0
            for line in lines:
                if "atoms" in line: num_atoms = int(line.split()[0]); break
            
            factor = int(num_atoms / len(charges))
            full_charges = charges * factor
            
            with open(self.data_final, 'w') as f:
                in_atoms = False
                idx = 0
                for line in lines:
                    if line.startswith("Atoms"):
                        in_atoms = True; f.write(line); continue
                    if in_atoms and len(line.strip()) == 0: f.write(line); continue
                    if in_atoms:
                        if line[0].isalpha(): in_atoms=False; f.write(line); continue
                        parts = line.split()
                        if len(parts) >= 6:
                            if idx < len(full_charges): parts[3] = f"{full_charges[idx]:.5f}"; idx += 1
                            f.write(" ".join(parts) + "\n")
                        else: f.write(line)
                    else: f.write(line)
        except Exception as e:
            print(f"   [WARN] {e}. Исп 0.0"); shutil.copy(self.data_raw, self.data_final)

    # === ИСПРАВЛЕННЫЙ АВТОДЕТЕКТОР СТИЛЕЙ ===
    def step_3_detect_styles(self):
        print(f"[3/5] Автоопределение стилей силового поля...")
        with open(self.data_final, 'r') as f:
            lines = f.readlines()
            
        def get_param_count(section_name):
            """Считает количество параметров (чисел) в первой строке секции Coeffs."""
            found_section = False
            for i, line in enumerate(lines):
                if section_name in line:
                    found_section = True
                    continue
                
                if found_section:
                    clean = line.strip()
                    if not clean: continue # Пустая строка
                    if clean[0].isalpha(): return 0 # Началась новая секция, данные кончились
                    
                    # Это строка с данными: ID param1 param2 ...
                    # Удаляем комментарии, если есть
                    clean = clean.split('#')[0].strip()
                    parts = clean.split()
                    # Первое число - это ID типа, остальные - параметры
                    return len(parts) - 1 
            return 0
            
        n_bond = get_param_count("Bond Coeffs")
        n_angle = get_param_count("Angle Coeffs")
        n_dihedral = get_param_count("Dihedral Coeffs")
        n_improper = get_param_count("Improper Coeffs")
        
        print(f"   Найдено параметров (колонок): Bond={n_bond}, Angle={n_angle}, Dihed={n_dihedral}, Improp={n_improper}")

        # --- ЛОГИКА ВЫБОРА СТИЛЕЙ ---
        
        # Bond
        if n_bond == 2: self.styles["bond"] = "harmonic"
        elif n_bond == 4: self.styles["bond"] = "morse"
        
        # Angle
        if n_angle == 2: self.styles["angle"] = "cosine/squared"
        elif n_angle == 3: self.styles["angle"] = "cosine/periodic"
        
        # Dihedral
        if n_dihedral == 3: self.styles["dihedral"] = "harmonic"
        elif n_dihedral == 1 or n_dihedral == 4: self.styles["dihedral"] = "fourier"
        
        # Improper (Здесь была ошибка)
        if n_improper == 2: 
            # 2 параметра может быть umbrella или harmonic. Umbrella чаще для inversion.
            self.styles["improper"] = "umbrella" 
        elif n_improper == 3: 
            # 3 параметра - это точно Fourier (K, d, n)
            self.styles["improper"] = "fourier"
        elif n_improper == 0:
            # Если импроперов нет, стиль не важен, но лучше оставить безопасный
            self.styles["improper"] = "umbrella" 

        print(f"   Принятые стили: {self.styles}")

    def step_4_gas(self):
        print(f"[4/5] Подготовка газа {self.gas}...")
        with open(self.gas_file, "w") as f: f.write(self.props["content"])

    def step_5_script(self, pressures):
        print(f"[5/5] Генерация скрипта LAMMPS...")
        with open(self.data_final) as f:
            txt = f.read()
            mof_types = int(re.search(r'(\d+)\s+atom types', txt).group(1))
        
        atoms = self.props["atoms"]
        start_gas = mof_types + 1
        
        script = [
            "log log.adsorption append",
            "units real",
            "atom_style full",
            "boundary p p p",
            "",
            f"pair_style lj/cut/coul/long {self.cutoff}",
            "kspace_style ewald 1.0e-4",
            "pair_modify mix arithmetic",
            "",
            # Стили из автодетектора
            f"bond_style      {self.styles['bond']}",
            f"angle_style     {self.styles['angle']}",
            f"dihedral_style  {self.styles['dihedral']}",
            f"improper_style  {self.styles['improper']}",
            "",
            f"read_data {self.data_final} extra/atom/types {len(atoms)} extra/bond/types 2 extra/angle/types 2",
            f"molecule gas_mol {self.gas_file}",
            ""
        ]
        
        for i, a in enumerate(atoms):
            tid = start_gas + i
            script.append(f"pair_coeff {tid} {tid} {a['eps']} {a['sig']}")
            
        script += [
            "",
            f"group gas type {start_gas}:{start_gas + len(atoms) - 1}",
            f"variable T equal {self.T}",
            "variable p_step index " + " ".join(map(str, pressures)),
            "label loop",
            "    variable P_atm equal v_p_step*0.986923",
            "    reset_timestep 0"
        ]
        
        kw = f"mol gas_mol pressure ${{P_atm}} fugacity_coeff 1.0"
        if self.props["rigid"]:
            kw += " rigid/small molecule"
            fix = f"fix myint gas rigid/nvt/small molecule temp ${{T}} ${{T}} 100.0"
        else:
            fix = f"fix myint gas nvt temp ${{T}} ${{T}} 100.0"
            
        seed = np.random.randint(10000, 99999)
        script.append(f"    fix mygcmc gas gcmc 10 500 500 {start_gas} {seed} ${{T}} -1.0 1.0 {kw}")
        script.append(f"    {fix}")
        
        div = 1
        if "N2" in self.gas or "CO2" in self.gas or "SO2" in self.gas: div = 3
        if "H2S" in self.gas: div = 4
        
        script += [
            f"    variable n_mols equal count(gas)/{div}",
            "    thermo_style custom step temp press v_P_atm v_n_mols",
            "    thermo 1000",
            "    run 5000",
            "    fix avg all ave/time 10 100 1000 v_n_mols file adsorption.txt mode scalar append",
            "    run 5000",
            "    unfix mygcmc",
            "    unfix myint",
            "    unfix avg",
            "    next p_step",
            "jump SELF loop"
        ]
        
        with open("in.adsorption", "w") as f:
            f.write("\n".join(script))

if __name__ == "__main__":
    cif_file = "MOF5.cif"
    gas = "CO2" 
    T = 298.0
    pressures = [0.1, 0.5, 1.0, 5.0, 10.0, 20.0]

    if os.path.exists(cif_file):
        if os.path.exists("adsorption.txt"): os.remove("adsorption.txt")
        
        sim = SimulationManager(cif_file, gas, T)
        sim.step_1_topology()
        sim.step_2_charges()
        sim.step_3_detect_styles() # Этот шаг теперь работает корректно
        sim.step_4_gas()
        sim.step_5_script(pressures)
        
        print("\nГОТОВО! Запустите: lmp_mpi -in in.adsorption")
    else:
        print(f"Файл {cif_file} не найден.")