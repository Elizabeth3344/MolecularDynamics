import sys
import os

def clean_cif(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    cell_params = {}
    atoms = []
    
    # Ключевые слова для параметров ячейки
    cell_keys = [
        '_cell_length_a', '_cell_length_b', '_cell_length_c',
        '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'
    ]

    # Парсинг файла
    in_loop = False
    headers = []
    
    print(f"Reading {input_file}...")

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # Считываем параметры ячейки
        for key in cell_keys:
            if line.startswith(key):
                parts = line.split()
                if len(parts) >= 2:
                    # Убираем скобки погрешности, если есть (например 10.123(5))
                    val = parts[1].split('(')[0]
                    cell_params[key] = val

        # Обработка лупов (таблиц с атомами)
        if line.startswith('loop_'):
            in_loop = True
            headers = []
            continue

        if in_loop:
            if line.startswith('_'):
                headers.append(line)
            else:
                # Это строка с данными (координатами атомов)
                # Нам нужно найти индексы колонок
                try:
                    # Попытка найти нужные колонки
                    # VESTA обычно дает: label, occ, x, y, z, ..., type
                    idx_label = -1
                    idx_x = -1
                    idx_y = -1
                    idx_z = -1
                    idx_type = -1

                    for i, h in enumerate(headers):
                        if '_atom_site_label' in h: idx_label = i
                        if '_atom_site_fract_x' in h: idx_x = i
                        if '_atom_site_fract_y' in h: idx_y = i
                        if '_atom_site_fract_z' in h: idx_z = i
                        if '_atom_site_type_symbol' in h: idx_type = i

                    parts = line.split()
                    
                    # Если нашли все координаты
                    if idx_x != -1 and idx_y != -1 and idx_z != -1:
                        label = parts[idx_label] if idx_label != -1 else "X"
                        
                        # Определяем тип атома (элемент)
                        if idx_type != -1:
                            atom_type = parts[idx_type]
                        else:
                            # Пытаемся угадать из лейбла (убираем цифры)
                            atom_type = ''.join([i for i in label if not i.isdigit()])

                        # Очистка типа от цифр и +/-, если они попали туда (для RASPA forcefield)
                        # Например Zn1 -> Zn
                        clean_type = ''.join([i for i in atom_type if i.isalpha()])

                        x = parts[idx_x].split('(')[0]
                        y = parts[idx_y].split('(')[0]
                        z = parts[idx_z].split('(')[0]

                        atoms.append({
                            'label': label,
                            'type': clean_type,
                            'x': x, 'y': y, 'z': z,
                            'charge': '0.0' # Дефолтный заряд
                        })
                except Exception as e:
                    # Если строка не парсится, возможно это конец лупа или мусор
                    in_loop = False
                    pass

    # Запись чистого файла для RASPA
    print(f"Writing {output_file} ({len(atoms)} atoms)...")
    
    with open(output_file, 'w') as f:
        f.write("data_RASPA_Simulation\n\n")
        f.write("_audit_creation_method 'Python Converter for RASPA'\n")
        
        # Записываем ячейку
        for key in cell_keys:
            if key in cell_params:
                f.write(f"{key:<30} {cell_params[key]}\n")
            else:
                # Дефолты если не найдено
                val = '90.0' if 'angle' in key else '10.0'
                f.write(f"{key:<30} {val}\n")

        # Принудительная P1 симметрия
        f.write("\n_symmetry_space_group_name_H-M  'P 1'\n")
        f.write("_symmetry_Int_Tables_number     1\n")
        f.write("\nloop_\n_symmetry_equiv_pos_as_xyz\n  'x,y,z'\n\n")

        # Атомы
        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")
        f.write("_atom_site_charge\n")

        for atom in atoms:
            f.write(f"{atom['label']:<10} {atom['type']:<6} {atom['x']:<12} {atom['y']:<12} {atom['z']:<12} {atom['charge']}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python convert_to_raspa.py <input_file.cif> <output_filename_without_extension>")
        sys.exit(1)
    
    input_cif = sys.argv[1]
    output_name = sys.argv[2]
    
    clean_cif(input_cif, output_name + ".cif")