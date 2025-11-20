import os
import sys
import re
import subprocess
import numpy as np
import shutil

# ==============================================================================
# 1. БАЗА ДАННЫХ ГАЗОВ (TraPPE Force Field)
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
    """Парсит заряды из CIF файла, созданного PyEQEQ."""
    charges = []
    try:
        with open(cif_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        return []

    loop_started = False
    charge_idx = -1
    headers = []
    data_start_idx = 0
    
    # Поиск колонки _atom_site_charge
    for i, line in enumerate(lines):
        clean = line.strip()
        if clean.startswith("loop_"):
            loop_started = True
            headers = []
            continue
        if loop_started and clean.startswith("_atom_site"):
            headers.append(clean)
        elif loop_started and not clean.startswith("_atom_site") and not clean.startswith("#") and len(clean) > 0:
            # Конец заголовков
            for h_i, h in enumerate(headers):
                if "_atom_site_charge" in h:
                    charge_idx = h_i
                    break
            data_start_idx = i
            break
    
    if charge_idx == -1:
        return []

    for i in range(data_start_idx, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith("loop_") or line.startswith("#"): break
        parts = line.split()
        if len(parts) > charge_idx:
            try:
                charges.append(float(parts[charge_idx]))
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
            raise ValueError(f"Газ {gas_name} не найден. Доступны: {list(GAS_DB.keys())}")
        
        self.props = GAS_DB[gas_name]
        base = os.path.splitext(os.path.basename(cif_file))[0]
        self.data_raw = f"data.{base}"
        self.data_final = f"data.{base}_final"
        self.gas_file = f"{gas_name.lower()}.mol"

    def step_1_topology(self):
        print(f"\n[1/4] Топология MOF (lammps-interface)...")
        # Генерируем без зарядов (они будут 0), чтобы не вызывать ошибки зависимостей
        cmd = [
            "lammps-interface", self.cif,
            "--force_field", self.ff,
            "--cutoff", str(self.cutoff),
            "--replication", "2x2x2",
            # "--outputcif", self.data_raw
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                # Fallback: пробуем без флага output для старых версий
                if "--output" in res.stderr:
                    cmd = cmd[:-2]
                    subprocess.run(cmd, check=True)
                    default = f"data.{os.path.basename(self.cif)}"
                    if os.path.exists(default): os.rename(default, self.data_raw)
                else:
                    print(f"ОШИБКА lammps-interface:\n{res.stderr}")
                    sys.exit(1)
        except FileNotFoundError:
            print("ОШИБКА: lammps-interface не найден. (pip install lammps-interface)")
            sys.exit(1)

    def step_2_charges(self):
        print(f"[2/4] Расчет зарядов MOF (PyEQEQ)...")
        
        # --- ИСПРАВЛЕНИЕ: Запуск через subprocess вместо pyeqeq.run() ---
        # PyEQEQ устанавливается как консольная команда 'pyeqeq'
        cmd = ["eqeq", self.cif]
        
        try:
            # Запускаем расчет
            subprocess.run(cmd, check=True, capture_output=True)
            
            # PyEQEQ создает файл с суффиксом _eqeq.cif или _1.2_eqeq.cif
            # Нам нужно найти этот файл
            base_name = os.path.splitext(self.cif)[0]
            candidates = [f for f in os.listdir('.') if f.startswith(base_name) and "eqeq" in f and f.endswith(".cif")]
            
            if not candidates:
                print("   [WARN] Файл с зарядами не найден. Будут использованы заряды 0.0")
                shutil.copy(self.data_raw, self.data_final)
                return
                
            # Берем самый свежий или первый попавшийся
            eqeq_cif = candidates[0]
            print(f"   Используем файл зарядов: {eqeq_cif}")

            # Парсинг
            charges = parse_charges_from_cif(eqeq_cif)
            
            if len(charges) == 0:
                print("   [WARN] Не удалось прочитать заряды из CIF. Используем 0.0")
                shutil.copy(self.data_raw, self.data_final)
                return

            print(f"   Успешно считано зарядов: {len(charges)}")
            
            # Внедрение в data файл
            with open(self.data_raw, 'r') as f: lines = f.readlines()
            
            # Считаем атомы
            num_atoms = 0
            for line in lines:
                if "atoms" in line and len(line.split()) == 2: 
                    num_atoms = int(line.split()[0])
                    break
            
            # Репликация (суперъячейка)
            factor = int(num_atoms / len(charges))
            if factor * len(charges) != num_atoms:
                print("   [WARN] Несовпадение числа атомов. Возможна ошибка топологии.")
            
            full_charges = charges * factor
            
            # Запись финального файла
            with open(self.data_final, 'w') as f:
                in_atoms = False
                idx = 0
                for line in lines:
                    if line.startswith("Atoms"):
                        in_atoms = True
                        f.write(line); continue
                    
                    if in_atoms and len(line.strip()) == 0:
                        f.write(line); continue
                        
                    if in_atoms:
                        if line[0].isalpha(): in_atoms = False; f.write(line); continue
                        parts = line.split()
                        if len(parts) >= 6:
                            # atomID molID type CHARGE x y z
                            if idx < len(full_charges):
                                parts[3] = f"{full_charges[idx]:.5f}"
                                idx += 1
                            f.write(" ".join(parts) + "\n")
                        else:
                            f.write(line)
                    else:
                        f.write(line)
            
            # Очистка временного файла eqeq (по желанию)
            # os.remove(eqeq_cif)

        except subprocess.CalledProcessError:
            print("   [WARN] Ошибка при запуске PyEQEQ. Проверьте установку (pip install pyeqeq).")
            shutil.copy(self.data_raw, self.data_final)
        except FileNotFoundError:
             print("   [WARN] Команда 'pyeqeq' не найдена.")
             shutil.copy(self.data_raw, self.data_final)

    def step_3_gas(self):
        print(f"[3/4] Подготовка газа {self.gas}...")
        with open(self.gas_file, "w") as f:
            f.write(self.props["content"])

    def step_4_script(self, pressures):
        print(f"[4/4] Генерация скрипта LAMMPS...")
        
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
            f"read_data {self.data_final} extra/atom/types {len(atoms)} extra/bond/types 2 extra/angle/types 2",
            f"molecule gas_mol {self.gas_file}",
            "",
            f"pair_style lj/cut/coul/long {self.cutoff}",
            "kspace_style ewald 1.0e-4",
            "pair_modify mix arithmetic",
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
    # 1. Имя CIF файла
    cif_file = "MOF5.cif"
    
    # 2. Газ (CO2, N2, CH4, H2S, SO2)
    gas = "CO2" 
    
    # 3. Параметры
    T = 298.0
    pressures = [0.1, 0.5, 1.0, 5.0, 10.0, 20.0]

    if os.path.exists(cif_file):
        if os.path.exists("adsorption.txt"): os.remove("adsorption.txt")
        
        sim = SimulationManager(cif_file, gas, T)
        sim.step_1_topology()
        sim.step_2_charges()
        sim.step_3_gas()
        sim.step_4_script(pressures)
        
        print("\nГОТОВО! Запустите: lmp_mpi -in in.adsorption")
    else:
        print(f"Файл {cif_file} не найден.")