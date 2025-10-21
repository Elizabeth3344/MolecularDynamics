
"""
Скрипт для добавления зарядов в LAMMPS data файл MOF-5
Два метода: простая замена и EQeq
"""

import numpy as np
import re
from collections import defaultdict

def read_lammps_data(filename):
    """
    Читает LAMMPS data файл и парсит всю необходимую информацию
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    data = {
        'header': [],
        'atoms': [],
        'bonds': [],
        'angles': [],
        'dihedrals': [],
        'impropers': [],
        'masses': {},
        'coeffs': defaultdict(list),
        'box': {}
    }
    
    # Находим секции
    section = None
    atom_style = None
    
    for i, line in enumerate(lines):
        line_strip = line.strip()
        
        # Пропускаем пустые строки и комментарии в начале
        if not line_strip or (line_strip.startswith('#') and section is None):
            data['header'].append(line)
            continue
            
        # Определяем секции (проверяем оба регистра)
        if 'Atoms' in line or ('atoms' in line_strip and 'atom types' not in line_strip):
            section = 'Atoms'
            continue
        elif 'Bonds' in line or ('bonds' in line_strip and 'bond types' not in line_strip):
            section = 'Bonds'
            continue
        elif 'Angles' in line or ('angles' in line_strip and 'angle types' not in line_strip):
            section = 'Angles'
            continue
        elif 'Dihedrals' in line or ('dihedrals' in line_strip and 'dihedral types' not in line_strip):
            section = 'Dihedrals'
            continue
        elif 'Impropers' in line or ('impropers' in line_strip and 'improper types' not in line_strip):
            section = 'Impropers'
            continue
        elif 'Masses' in line:
            section = 'Masses'
            continue
        elif 'Coeffs' in line:
            section = line.strip()
            continue
        elif line_strip in ['', ' ']:
            section = None
            continue
            
        # Парсим данные в зависимости от секции
        if section == 'Atoms':
            parts = line.split()
            if len(parts) >= 7:  # full atom style
                atom_data = {
                    'id': int(parts[0]),
                    'mol_id': int(parts[1]),
                    'type': int(parts[2]),
                    'charge': float(parts[3]),
                    'x': float(parts[4]),
                    'y': float(parts[5]),
                    'z': float(parts[6])
                }
                data['atoms'].append(atom_data)
                
        elif section == 'Masses':
            parts = line.split()
            if len(parts) >= 2:
                type_id = int(parts[0])
                mass = float(parts[1])
                # Сохраняем комментарий с названием типа атома
                comment = ' '.join(parts[2:]) if len(parts) > 2 else ''
                data['masses'][type_id] = {'mass': mass, 'comment': comment}
    
    return data

def assign_charges_simple(data):
    """
    Простое присвоение зарядов на основе известных значений для MOF-5
    """
    # Известные заряды для MOF-5 (из литературы)
    charge_map = {
        1: 1.411,    # Zn3+2
        2: -1.493,   # O_2 (кислород в Zn4O)
        3: -0.659,   # O_R (кислород в карбоксилате)
        4: 0.0,      # C_R (будем корректировать)
        5: 0.150     # H_
    }
    
    # Присваиваем базовые заряды
    for atom in data['atoms']:
        atom['charge'] = charge_map.get(atom['type'], 0.0)
    
    # Теперь нужно различить углероды в карбоксильных группах и в бензольных кольцах
    # Для этого анализируем координационное окружение
    
    # Создаем словарь позиций атомов
    positions = {}
    for atom in data['atoms']:
        positions[atom['id']] = np.array([atom['x'], atom['y'], atom['z']])
    
    # Для каждого углерода проверяем, связан ли он с кислородами
    for atom in data['atoms']:
        if atom['type'] == 4:  # C_R
            atom_pos = positions[atom['id']]
            
            # Считаем кислороды в радиусе 1.5 Å
            oxygen_count = 0
            for other in data['atoms']:
                if other['type'] in [2, 3]:  # O_2 или O_R
                    dist = np.linalg.norm(positions[other['id']] - atom_pos)
                    if dist < 1.5:
                        oxygen_count += 1
            
            # Если углерод связан с 2 кислородами - это карбоксильная группа
            if oxygen_count >= 2:
                atom['charge'] = 0.868  # C в карбоксилате
            else:
                atom['charge'] = -0.186  # C в бензольном кольце
    
    # Проверяем электронейтральность
    total_charge = sum(atom['charge'] for atom in data['atoms'])
    print(f"Суммарный заряд после присвоения: {total_charge:.6f}")
    
    # Небольшая корректировка для точной нейтральности
    if abs(total_charge) > 0.001:
        correction = -total_charge / len(data['atoms'])
        for atom in data['atoms']:
            atom['charge'] += correction
        
        total_charge = sum(atom['charge'] for atom in data['atoms'])
        print(f"Суммарный заряд после корректировки: {total_charge:.6f}")
    
    return data

def assign_charges_eqeq(data, chi_values=None, hardness_values=None):
    """
    Присвоение зарядов методом EQeq (Electronegativity Equalization)
    Упрощенная версия для демонстрации принципа
    """
    
    # Параметры электроотрицательности (χ) и химической жесткости (J) для элементов
    # Значения из статьи Wilmer et al., J. Phys. Chem. Lett. 2012
    if chi_values is None:
        chi_values = {
            'Zn': 4.45,
            'O': 8.50,
            'C': 5.34,
            'H': 4.53
        }
    
    if hardness_values is None:
        hardness_values = {
            'Zn': 6.35,
            'O': 10.83,
            'C': 5.34,
            'H': 6.90
        }
    
    # Определяем соответствие типов атомов элементам
    type_to_element = {
        1: 'Zn',
        2: 'O',
        3: 'O',
        4: 'C',
        5: 'H'
    }
    
    n_atoms = len(data['atoms'])
    
    # Создаем матрицы для системы линейных уравнений
    # A·q = b, где q - вектор зарядов
    A = np.zeros((n_atoms + 1, n_atoms + 1))
    b = np.zeros(n_atoms + 1)
    
    # Координаты атомов
    coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in data['atoms']])
    
    # Заполняем матрицу A
    for i in range(n_atoms):
        element_i = type_to_element[data['atoms'][i]['type']]
        
        # Диагональные элементы: 2·J_i
        A[i, i] = 2 * hardness_values[element_i]
        
        # Недиагональные элементы: кулоновское взаимодействие
        for j in range(n_atoms):
            if i != j:
                r_ij = np.linalg.norm(coords[i] - coords[j])
                # Используем экранированное кулоновское взаимодействие
                # с параметром экранирования λ = 0.1
                lambda_screen = 0.1
                A[i, j] = 14.4 * np.exp(-lambda_screen * r_ij) / r_ij
        
        # Последний столбец - условие электронейтральности
        A[i, n_atoms] = 1.0
        A[n_atoms, i] = 1.0
        
        # Правая часть: -χ_i
        b[i] = -chi_values[element_i]
    
    # Последнее уравнение: сумма зарядов = 0
    A[n_atoms, n_atoms] = 0.0
    b[n_atoms] = 0.0
    
    # Решаем систему линейных уравнений
    try:
        solution = np.linalg.solve(A, b)
        charges = solution[:n_atoms]
        
        # Присваиваем заряды атомам
        for i, atom in enumerate(data['atoms']):
            atom['charge'] = charges[i]
        
        # Проверка
        total_charge = sum(atom['charge'] for atom in data['atoms'])
        print(f"EQeq: Суммарный заряд = {total_charge:.6f}")
        
        # Статистика по типам атомов
        charge_stats = defaultdict(list)
        for atom in data['atoms']:
            charge_stats[atom['type']].append(atom['charge'])
        
        print("\nСредние заряды по типам атомов (EQeq):")
        for type_id, charges in charge_stats.items():
            element = type_to_element[type_id]
            avg_charge = np.mean(charges)
            std_charge = np.std(charges)
            print(f"  Тип {type_id} ({element}): {avg_charge:+.4f} ± {std_charge:.4f}")
            
    except np.linalg.LinAlgError:
        print("Ошибка при решении системы уравнений EQeq")
        print("Используем запасной метод - итерационный")
        
        # Простой итерационный метод
        charges = np.zeros(n_atoms)
        for iteration in range(100):
            new_charges = np.zeros(n_atoms)
            
            for i in range(n_atoms):
                element_i = type_to_element[data['atoms'][i]['type']]
                
                # Вычисляем эффективную электроотрицательность
                chi_eff = chi_values[element_i]
                for j in range(n_atoms):
                    if i != j:
                        r_ij = np.linalg.norm(coords[i] - coords[j])
                        chi_eff += 14.4 * charges[j] / r_ij
                
                # Обновляем заряд
                new_charges[i] = -chi_eff / (2 * hardness_values[element_i])
            
            # Нормализуем для нейтральности
            new_charges -= np.mean(new_charges)
            
            # Проверяем сходимость
            if np.max(np.abs(new_charges - charges)) < 1e-6:
                break
                
            charges = new_charges
        
        # Присваиваем заряды
        for i, atom in enumerate(data['atoms']):
            atom['charge'] = charges[i]
    
    return data

def write_lammps_data(data, filename, original_filename):
    """
    Записывает обновленный LAMMPS data файл с новыми зарядами
    """
    # Читаем оригинальный файл для сохранения форматирования
    with open(original_filename, 'r') as f:
        original_lines = f.readlines()
    
    # Находим секцию Atoms в оригинальном файле
    atoms_start = None
    for i, line in enumerate(original_lines):
        if 'Atoms' in line and 'atoms' not in line:
            atoms_start = i + 2  # Пропускаем заголовок и пустую строку
            break
    
    # Создаем новый файл
    with open(filename, 'w') as f:
        # Копируем все до секции Atoms
        atoms_written = False
        i = 0
        
        while i < len(original_lines):
            line = original_lines[i]
            
            # Если дошли до атомов
            if atoms_start and i >= atoms_start and not atoms_written:
                # Проверяем, что это строка с атомом
                parts = line.split()
                if len(parts) >= 7 and parts[0].isdigit():
                    # Записываем все атомы с новыми зарядами
                    for atom in data['atoms']:
                        f.write(f"     {atom['id']:3d}      {atom['mol_id']:3d}        "
                               f"{atom['type']:1d}     {atom['charge']:8.5f}    "
                               f"{atom['x']:8.5f}   {atom['y']:8.5f}   {atom['z']:8.5f}\n")
                    
                    # Пропускаем оригинальные строки с атомами
                    while i < len(original_lines):
                        parts = original_lines[i].split()
                        if len(parts) < 7 or not parts[0].isdigit():
                            break
                        i += 1
                    
                    atoms_written = True
                    continue
            
            # Записываем строку как есть
            f.write(line)
            i += 1

# Основная функция
def main():
    """
    Основная функция для обработки MOF-5 data файла
    """
    input_file = "/home/chernysheva/MolecularDynamics/Tasks/MOF-5/last_vers/data.lammpsInterfaceMOF5"
    
    print("=== Добавление зарядов в MOF-5 ===\n")
    
    # Читаем данные
    print("Чтение файла...")
    data = read_lammps_data(input_file)
    print(f"Прочитано {len(data['atoms'])} атомов")
    
    # Метод 1: Простые заряды из литературы
    print("\n--- Метод 1: Заряды из литературы ---")
    data_simple = assign_charges_simple(data.copy())
    write_lammps_data(data_simple, "/home/chernysheva/MolecularDynamics/Tasks/MOF-5/last_vers/data.MOF5_charged_simple", input_file)
    print("Файл сохранен: data.MOF5_charged_simple")
    
    # Статистика
    charge_stats = defaultdict(list)
    for atom in data_simple['atoms']:
        charge_stats[atom['type']].append(atom['charge'])
    
    print("\nСредние заряды по типам атомов (простой метод):")
    for type_id in sorted(charge_stats.keys()):
        charges = charge_stats[type_id]
        avg_charge = np.mean(charges)
        element = data['masses'][type_id]['comment']
        print(f"  Тип {type_id} {element}: {avg_charge:+.4f} ({len(charges)} атомов)")
    
    # Метод 2: EQeq
    print("\n--- Метод 2: Расчет зарядов методом EQeq ---")
    data_eqeq = assign_charges_eqeq(read_lammps_data(input_file))
    write_lammps_data(data_eqeq, "/home/chernysheva/MolecularDynamics/Tasks/MOF-5/last_vers/data.MOF5_charged_eqeq", input_file)
    print("Файл сохранен: data.MOF5_charged_eqeq")
    
    print("\n=== Готово! ===")
    print("Созданы два файла с зарядами:")
    print("1. data.MOF5_charged_simple - заряды из литературы")  
    print("2. data.MOF5_charged_eqeq - заряды рассчитаны методом EQeq")

if __name__ == "__main__":
    main()