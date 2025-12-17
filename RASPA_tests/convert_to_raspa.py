import sys
import os
import shutil
import glob

# Пробуем подключить PyEQEq
try:
    from pyeqeq import run_eqeq
    HAS_PYEQEQ = True
except ImportError:
    HAS_PYEQEQ = False
    print("ВНИМАНИЕ: Библиотека 'pyeqeq' не установлена. Заряды будут 0.0.")

def handle_gas_file(gas_name, raspa_dir):
    """
    1. Проверяет наличие газа в локальной папке ForceField.
    2. Если нет - ищет в системе RASPA и сохраняет в ForceField.
    3. Копирует из ForceField в корень для работы.
    """
    ff_dir = "ForceField"
    # Создаем хранилище, если нет
    if not os.path.exists(ff_dir):
        os.makedirs(ff_dir)

    target_def = f"{gas_name}.def"
    local_storage_path = os.path.join(ff_dir, target_def)
    
    # --- ШАГ 1: Поиск файла ---
    source_to_copy = None

    if os.path.exists(local_storage_path):
        # Если файл уже есть у нас в папке ForceField
        print(f"-> Газ найден локально в архиве: {local_storage_path}")
        source_to_copy = local_storage_path
    else:
        # Если нет, ищем в системной библиотеке RASPA
        print(f"-> Газ не найден локально, ищем в системе RASPA...")
        search_path = os.path.join(raspa_dir, "share", "raspa", "molecules")
        found = glob.glob(os.path.join(search_path, "**", target_def), recursive=True)
        
        if not found:
            print(f"ОШИБКА: Файл {target_def} не найден ни в ForceField, ни в RASPA.")
            sys.exit(1)
        
        # Сохраняем найденное в наш архив на будущее
        print(f"-> Найдено в системе: {found[0]}")
        shutil.copy(found[0], local_storage_path)
        print(f"-> Сохранено в архив: {local_storage_path}")
        source_to_copy = local_storage_path

    # --- ШАГ 2: Копирование в рабочую зону ---
    # RASPA видит файлы только в корне
    shutil.copy(source_to_copy, target_def)
    print(f"-> Газ готов к работе: {target_def}")

def clean_and_charge_cif(input_file, output_file):
    """Читает CIF, чистит, считает заряды, сохраняет."""
    print(f"Обработка структуры: {input_file}")
    
    with open(input_file, 'r') as f:
        lines = f.readlines()

    cell_params = {}
    atoms = []
    
    cell_keys = ['_cell_length_a', '_cell_length_b', '_cell_length_c', 
                 '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma']

    in_loop = False
    headers = []
    
    # --- Парсинг ---
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'): continue

        for key in cell_keys:
            if line.startswith(key):
                cell_params[key] = line.split()[1].split('(')[0]

        if line.startswith('loop_'):
            in_loop = True; headers = []; continue

        if in_loop:
            if line.startswith('_'): headers.append(line)
            else:
                try:
                    parts = line.split()
                    idx_label = -1; idx_x = -1; idx_y = -1; idx_z = -1; idx_type = -1

                    for i, h in enumerate(headers):
                        if 'label' in h: idx_label = i
                        if 'fract_x' in h: idx_x = i
                        if 'fract_y' in h: idx_y = i
                        if 'fract_z' in h: idx_z = i
                        if 'type_symbol' in h: idx_type = i

                    if idx_x != -1 and idx_y != -1 and idx_z != -1:
                        label = parts[idx_label] if idx_label != -1 else "X"
                        # Чистим тип (Zn1 -> Zn)
                        if idx_type != -1: raw_type = parts[idx_type]
                        else: raw_type = label
                        clean_type = ''.join([c for c in raw_type if c.isalpha()])

                        atoms.append({
                            'label': label,
                            'type': clean_type,
                            'x': parts[idx_x].split('(')[0],
                            'y': parts[idx_y].split('(')[0],
                            'z': parts[idx_z].split('(')[0]
                        })
                except: in_loop = False

    # --- Временный файл ---
    temp_cif = "temp_clean.cif"
    with open(temp_cif, 'w') as f:
        f.write("data_Clean\n_symmetry_space_group_name_H-M 'P 1'\n_symmetry_Int_Tables_number 1\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n")
        for k in cell_keys:
            val = cell_params.get(k, "90.0" if "angle" in k else "10.0")
            f.write(f"{k:<30} {val}\n")
        f.write("\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_charge\n")
        for a in atoms:
            f.write(f"{a['label']:<10} {a['type']:<6} {a['x']:<12} {a['y']:<12} {a['z']:<12} 0.0\n")

    # --- Заряды (PyEQEq) ---
    if HAS_PYEQEQ:
        print("Запуск PyEQEq...")
        try:
            run_eqeq(temp_cif, output_file, method='eqeq')
            print("Заряды успешно посчитаны.")
        except Exception as e:
            print(f"Ошибка PyEQEq: {e}. Оставляем 0.0.")
            shutil.copy(temp_cif, output_file)
    else:
        print("PyEQEq не найден. Оставляем заряды 0.0.")
        shutil.copy(temp_cif, output_file)
    
    if os.path.exists(temp_cif): os.remove(temp_cif)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 convert_to_raspa.py <cif_path> <gas_name>")
        sys.exit(1)
        
    convert_to_raspa_path = sys.argv[1]
    gas_name = sys.argv[2]
    raspa_dir = os.environ.get("RASPA_DIR")
    
    clean_and_charge_cif(convert_to_raspa_path, "simulation_structure.cif")
    handle_gas_file(gas_name, raspa_dir)