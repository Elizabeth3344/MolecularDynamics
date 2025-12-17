import os
import glob

# Настройки: какие папки объединять
# GarciaPerez2006 = UFF для MOF
# TraPPE = Точный CO2
SOURCES = ["GarciaPerez2006", "TraPPE"] 

def parse_def_file(filepath, file_type):
    """Читает файл .def и возвращает список строк с данными (без заголовков)"""
    data_lines = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    start_reading = False
    for line in lines:
        # Ищем начало данных
        if file_type == "pseudo_atoms":
            # Данные начинаются после строки заголовка колонок (#type print...)
            if line.strip().startswith("#type"):
                start_reading = True
                continue
        elif file_type == "mixing_rules":
            # Данные начинаются после заголовка (# type interaction)
            if line.strip().startswith("# type") or line.strip().startswith("#type"):
                start_reading = True
                continue
        
        # Если читаем данные и строка не комментарий и не mixing rule type
        if start_reading:
            if "Lorentz-Berthelot" in line: continue
            if line.strip() and not line.startswith("#"):
                data_lines.append(line)
    return data_lines

def build():
    raspa_dir = os.environ.get("RASPA_DIR")
    if not raspa_dir:
        print("ОШИБКА: Не найдена переменная окружения RASPA_DIR")
        return

    ff_path = os.path.join(raspa_dir, "share", "raspa", "forcefield")
    output_dir = "ForceField"
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    print(f"--- Сборка ForceField из: {SOURCES} ---")

    # 1. СБОРКА PSEUDO_ATOMS.DEF
    all_atoms = []
    # Сначала добавляем Гелий (его нет в папках)
    # Формат: type print as chem oxidation mass charge polarization B-factor radii connectivity anisotropic anisotropic-type tinker-type
    all_atoms.append("He         yes     He    He    0           4.0026    0.0      0.0          1.0      1.0    0            0           relative         0\n")
    
    for source in SOURCES:
        fpath = os.path.join(ff_path, source, "pseudo_atoms.def")
        if os.path.exists(fpath):
            print(f"Reading atoms from {source}...")
            all_atoms.extend(parse_def_file(fpath, "pseudo_atoms"))
        else:
            print(f"WARNING: {source} not found in {ff_path}")

    # Запись pseudo_atoms.def
    with open(os.path.join(output_dir, "pseudo_atoms.def"), 'w') as f:
        f.write(f"# number of pseudo atoms\n{len(all_atoms)}\n")
        f.write("#type      print   as    chem  oxidation   mass      charge   polarization B-factor radii  connectivity anisotropic anisotropic-type   tinker-type\n")
        for line in all_atoms:
            f.write(line)

    # 2. СБОРКА MIXING_RULES.DEF
    all_rules = []
    # Добавляем Гелий
    all_rules.append("He         lennard-jones    10.9      2.64\n")

    for source in SOURCES:
        fpath = os.path.join(ff_path, source, "force_field_mixing_rules.def")
        if os.path.exists(fpath):
            print(f"Reading rules from {source}...")
            all_rules.extend(parse_def_file(fpath, "mixing_rules"))

    # Запись force_field_mixing_rules.def
    with open(os.path.join(output_dir, "force_field_mixing_rules.def"), 'w') as f:
        f.write("# general interaction definitions\n1\n0\n")
        f.write(f"# number of interaction definitions\n{len(all_rules)}\n")
        f.write("# type     interaction  parameters\n")
        for line in all_rules:
            f.write(line)
        f.write("# mixing rules\nLorentz-Berthelot\n")

    # 3. КОПИРУЕМ force_field.def (Обычно берем из первого источника)
    base_ff = os.path.join(ff_path, SOURCES[0], "force_field.def")
    if os.path.exists(base_ff):
        # Копируем и переименовываем в force_field.def (на всякий случай)
        with open(base_ff, 'r') as f_in, open(os.path.join(output_dir, "force_field.def"), 'w') as f_out:
            f_out.write(f_in.read())
    else:
        # Создаем заглушку, если файла нет (он часто не нужен, но RASPA может просить)
        with open(os.path.join(output_dir, "force_field.def"), 'w') as f:
            f.write("# rules\nTruncated\nshifted\n")

    print("--- Успешно! Файлы сохранены в папку ForceField/ ---")

if __name__ == "__main__":
    build()