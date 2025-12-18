import os
import subprocess
import shutil
import sys
import glob
import csv

# ================= НАСТРОЙКИ ЭКСПЕРИМЕНТА =================
MOF_FILE = "MOF5.cif"           # Имя файла в папке MOF/
GAS_NAME = "CO2"                # Имя газа
TEMPERATURES = [273, 298, 323]  # Кельвины
PRESSURES = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0] # Бары
# ==========================================================

# Пути
BASE_DIR = os.getcwd()
MOF_PATH = os.path.join(BASE_DIR, "MOF", MOF_FILE)
RASPA_DIR = os.environ.get("RASPA_DIR")
SIMULATE_BIN = os.path.join(RASPA_DIR, "bin", "simulate") if RASPA_DIR else None

def check_setup():
    if not RASPA_DIR:
        print("ОШИБКА: RASPA_DIR не установлена.")
        sys.exit(1)
    if not os.path.exists(MOF_PATH):
        print(f"ОШИБКА: Файл {MOF_PATH} не найден.")
        sys.exit(1)

def parse_output(directory, temp, press):
    """Читает результаты из файла .data"""
    out_dir = os.path.join(directory, "Output", "System_0")
    files = glob.glob(os.path.join(out_dir, "*.data"))
    
    loading_abs = 0.0
    loading_exc = 0.0
    heat = 0.0
    
    if files:
        with open(files[0], 'r') as f:
            for line in f:
                if "Average loading absolute [mol/kg framework]" in line:
                    try: loading_abs = float(line.split()[5])
                    except: pass
                elif "Average loading excess [mol/kg framework]" in line:
                    try: loading_exc = float(line.split()[5])
                    except: pass
                elif "Heat of desorption" in line and "[K]" in line:
                    try:
                        # Ищем число перед [K]
                        parts = line.split()
                        for i, p in enumerate(parts):
                            if "K" in p or "[K]" in p:
                                heat = float(parts[i-1])
                                break
                    except: pass

    return {"T": temp, "P": press, "Abs_mol_kg": loading_abs, "Exc_mol_kg": loading_exc, "Heat_K": heat}

def run_isotherms():
    print("\n--- ЭТАП: Расчет изотерм адсорбции ---")
    
    results = [] 

    # 1. Готовим структуру и газ один раз (через наш конвертер)
    # Это создаст simulation_structure.cif и CO2.def в корне
    print(f"-> Подготовка CIF и газа {GAS_NAME}...")
    subprocess.call([sys.executable, "convert_to_raspa.py", MOF_PATH, GAS_NAME])
    
    # 2. Читаем шаблон
    with open("simulation.input.template", 'r') as f:
        template_content = f.read()

    # 3. Цикл по точкам
    for temp in TEMPERATURES:
        for press_bar in PRESSURES:
            press_pa = press_bar * 100000.0
            print(f"   -> Запуск: T={temp}K, P={press_bar} bar...")

            # Создаем папку для точки
            point_dir = os.path.join(BASE_DIR, "Results", f"{temp}K", f"{press_bar}bar")
            if not os.path.exists(point_dir): os.makedirs(point_dir)

            # Копируем рабочие файлы в папку точки
            shutil.copy("simulation_structure.cif", point_dir)
            shutil.copy(f"{GAS_NAME}.def", point_dir)
            # Копируем всё из ForceField
            for f in glob.glob("ForceField/*.def"): shutil.copy(f, point_dir)

            # Настраиваем Input (подстановка значений)
            movie_flag = "yes" if press_bar == PRESSURES[0] else "no"

            sim_content = template_content \
                .replace("TEMP_PLACEHOLDER", str(temp)) \
                .replace("PRES_PLACEHOLDER", str(press_pa)) \
                .replace("REPLACE_ME", GAS_NAME) \
                .replace("MOVIE_FLAG", movie_flag)

            with open(os.path.join(point_dir, "simulation.input"), 'w') as f:
                f.write(sim_content)

            # Запуск RASPA (без вывода в консоль)
            subprocess.run([SIMULATE_BIN, "simulation.input"], cwd=point_dir, stdout=subprocess.DEVNULL)

            # Сбор результатов
            res = parse_output(point_dir, temp, press_bar)
            results.append(res)

    return results

def save_csv(results):
    filename = "final_results.csv"
    print(f"\n--- Сохранение результатов в {filename} ---")
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["T", "P", "Abs_mol_kg", "Exc_mol_kg", "Heat_K"])
        writer.writeheader()
        writer.writerows(results)
    print("Готово!")

if __name__ == "__main__":
    check_setup()
    # Запускаем сразу изотермы (без Void Fraction)
    data = run_isotherms()
    save_csv(data)