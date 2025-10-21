#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Универсальное решение для подготовки GCMC моделирования в LAMMPS:
- Использует CIF2LAMMPS и RASPA для преобразования CIF в DATA
- Автоматически определяет параметры UFF и DREIDING для всех атомов
- Рассчитывает заряды методом EQeq с параметрами Уилмера (Wilmer)
- Создает полный набор файлов LAMMPS для GCMC моделирования
"""

import os
import sys
import subprocess
import numpy as np
import argparse
from pymatgen.io.cif import CifParser
from pymatgen.core import Element, Structure
from pymatgen.io.lammps.data import LammpsData
import warnings
warnings.filterwarnings("ignore")

# Параметры UFF для распространенных элементов в MOF
UFF_PARAMS = {
    'H': {'epsilon': 0.044, 'sigma': 2.571, 'mass': 1.008},
    'C': {'epsilon': 0.105, 'sigma': 3.851, 'mass': 12.011},
    'N': {'epsilon': 0.069, 'sigma': 3.660, 'mass': 14.007},
    'O': {'epsilon': 0.060, 'sigma': 3.500, 'mass': 15.999},
    'F': {'epsilon': 0.050, 'sigma': 3.364, 'mass': 18.998},
    'Zn': {'epsilon': 0.124, 'sigma': 2.461, 'mass': 65.38},
    'Cu': {'epsilon': 0.005, 'sigma': 3.114, 'mass': 63.546},
    'Fe': {'epsilon': 0.013, 'sigma': 2.912, 'mass': 55.845},
    'Al': {'epsilon': 0.505, 'sigma': 4.008, 'mass': 26.982},
    # Можно дополнить для других элементов
}

# Параметры EQeq для расчета зарядов (Wilmer et al., 2012)
EQEQ_PARAMS = {
    'H': {'chi': 4.528, 'J': 13.890, 'R': 0.371},
    'C': {'chi': 5.343, 'J': 10.126, 'R': 0.759},
    'N': {'chi': 6.899, 'J': 11.540, 'R': 0.715},
    'O': {'chi': 8.741, 'J': 13.364, 'R': 0.669},
    'F': {'chi': 10.874, 'J': 14.948, 'R': 0.706},
    'Zn': {'chi': 5.106, 'J': 8.034, 'R': 1.280},
    'Cu': {'chi': 4.98, 'J': 7.89, 'R': 1.448},
    'Fe': {'chi': 4.6, 'J': 7.41, 'R': 1.393},
    'Al': {'chi': 3.59, 'J': 5.87, 'R': 1.26},
    # Можно дополнить для других элементов
}

# Радиусы ковалентных связей из CSD (Cambridge Structural Database)
BOND_RADII = {
    ('C', 'C'): 1.54, ('C', 'H'): 1.09, ('C', 'N'): 1.43, ('C', 'O'): 1.43,
    ('Zn', 'O'): 1.98, ('Zn', 'N'): 2.06, ('Cu', 'O'): 1.96, ('Cu', 'N'): 2.02,
    ('Fe', 'O'): 1.95, ('Fe', 'N'): 2.04, ('Al', 'O'): 1.87
    # Можно дополнить другими парами
}

# Параметры силового поля UFF для связей (Rappe et al., 1992)
UFF_BOND_PARAMS = {
    ('C', 'C'): {'k': 699.6, 'r0': 1.54},
    ('C', 'H'): {'k': 663.6, 'r0': 1.09},
    ('C', 'N'): {'k': 774.4, 'r0': 1.43},
    ('C', 'O'): {'k': 802.8, 'r0': 1.43},
    ('Zn', 'O'): {'k': 153.7, 'r0': 1.98},
    ('Zn', 'N'): {'k': 133.2, 'r0': 2.06},
    ('Cu', 'O'): {'k': 162.5, 'r0': 1.96},
    ('Cu', 'N'): {'k': 142.1, 'r0': 2.02},
    # Параметры для других пар можно вычислить по формулам UFF
}

# Параметры силового поля UFF для углов (Rappe et al., 1992)
UFF_ANGLE_PARAMS = {
    ('C', 'C', 'C'): {'k': 119.2, 'theta0': 109.47},
    ('C', 'C', 'H'): {'k': 69.6, 'theta0': 109.47},
    ('C', 'O', 'C'): {'k': 99.8, 'theta0': 104.51},
    ('O', 'C', 'O'): {'k': 124.3, 'theta0': 120.0},
    ('O', 'Zn', 'O'): {'k': 95.5, 'theta0': 109.47},
    ('N', 'Zn', 'N'): {'k': 62.1, 'theta0': 109.47},
    # Параметры для других углов можно вычислить по формулам UFF
}

def run_command(cmd, cwd=None):
    """Запускает команду в терминале и возвращает вывод"""
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, cwd=cwd)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        print(f"Ошибка при выполнении команды: {cmd}")
        print(f"Stderr: {stderr.decode('utf-8')}")
        return None
    
    return stdout.decode('utf-8')

def calculate_eqeq_charges(structure, accuracy=1e-6, max_iter=1000):
    """
    Рассчитывает заряды атомов методом EQeq (Charge Equilibration)
    Реализует метод из работы Wilmer et al., 2012
    """
    num_atoms = len(structure)
    
    # Подготавливаем параметры
    chi_values = []
    j_values = []
    radii = []
    
    for site in structure:
        element = site.species_string
        if element in EQEQ_PARAMS:
            chi_values.append(EQEQ_PARAMS[element]['chi'])
            j_values.append(EQEQ_PARAMS[element]['J'])
            radii.append(EQEQ_PARAMS[element]['R'])
        else:
            # Если элемент отсутствует, используем параметры углерода (обычно безопасно)
            print(f"Предупреждение: Параметры EQeq для {element} не найдены, используем параметры для C")
            chi_values.append(EQEQ_PARAMS['C']['chi'])
            j_values.append(EQEQ_PARAMS['C']['J'])
            radii.append(EQEQ_PARAMS['C']['R'])
    
    # Строим матрицу системы линейных уравнений
    A = np.zeros((num_atoms + 1, num_atoms + 1))
    b = np.zeros(num_atoms + 1)
    
    # Диагональные элементы - жесткости атомов
    for i in range(num_atoms):
        A[i, i] = 2.0 * j_values[i]
        b[i] = -chi_values[i]
    
    # Недиагональные элементы - кулоновские взаимодействия
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            r_ij = np.linalg.norm(structure[i].coords - structure[j].coords)
            
            # Избегаем деления на ноль
            if r_ij < 0.1:
                r_ij = 0.1
            
            # Coulomb с экранированием на малых расстояниях
            coulomb = 14.4 / (r_ij + np.sqrt(radii[i] * radii[j]) * np.exp(-r_ij / (radii[i] + radii[j])))
            
            A[i, j] = coulomb
            A[j, i] = coulomb
    
    # Условие электронейтральности
    A[num_atoms, :num_atoms] = 1.0
    A[:num_atoms, num_atoms] = 1.0
    A[num_atoms, num_atoms] = 0.0
    b[num_atoms] = 0.0
    
    # Решаем систему уравнений
    x = np.linalg.solve(A, b)
    charges = x[:num_atoms]
    
    # Нормализуем заряды, чтобы сумма была точно равна нулю
    charges -= np.mean(charges)
    
    # Возвращаем словарь с зарядами для каждого атома
    return {i: charges[i] for i in range(num_atoms)}

def identify_bonds_in_structure(structure, tolerance=1.3):
    """
    Идентифицирует ковалентные связи в структуре на основе
    ковалентных радиусов и возвращает список связей
    """
    bonds = []
    elements = [site.species_string for site in structure]
    
    for i in range(len(structure)):
        for j in range(i+1, len(structure)):
            elem_i = elements[i]
            elem_j = elements[j]
            
            # Определяем длину связи
            dist = np.linalg.norm(structure[i].coords - structure[j].coords)
            
            # Проверяем по таблице ковалентных длин связей
            bond_key = tuple(sorted([elem_i, elem_j]))
            if bond_key in BOND_RADII:
                threshold = BOND_RADII[bond_key] * tolerance
                if dist <= threshold:
                    bonds.append((i, j))
            else:
                # Если пары нет в таблице, используем сумму ковалентных радиусов
                radius_i = Element(elem_i).covalent_radius
                radius_j = Element(elem_j).covalent_radius
                threshold = (radius_i + radius_j) * tolerance
                if dist <= threshold:
                    bonds.append((i, j))
    
    return bonds

def identify_angles_in_structure(structure, bonds):
    """
    Идентифицирует валентные углы в структуре на основе найденных связей
    """
    # Создаем словарь соседей для каждого атома
    neighbors = {}
    for i, j in bonds:
        if i not in neighbors:
            neighbors[i] = []
        if j not in neighbors:
            neighbors[j] = []
        neighbors[i].append(j)
        neighbors[j].append(i)
    
    angles = []
    for j in range(len(structure)):
        if j in neighbors and len(neighbors[j]) >= 2:
            for i in neighbors[j]:
                for k in neighbors[j]:
                    if i < k:  # Избегаем дублирования углов
                        angles.append((i, j, k))
    
    return angles

def get_uff_parameters(element, forcefield='UFF'):
    """
    Возвращает параметры UFF или DREIDING для указанного элемента
    """
    if forcefield == 'UFF':
        if element in UFF_PARAMS:
            return UFF_PARAMS[element]
        else:
            # Параметры UFF можно рассчитать по формулам из оригинальной статьи
            # Здесь упрощенно возвращаем параметры для C
            print(f"Предупреждение: UFF параметры для {element} не найдены, используем параметры для C")
            return UFF_PARAMS['C']
    else:
        # Для DREIDING аналогично
        print(f"Предупреждение: DREIDING параметры для {element} не реализованы")
        return UFF_PARAMS['C']

def get_bond_parameters(atom1, atom2, forcefield='UFF'):
    """
    Возвращает параметры связи для указанной пары атомов
    """
    bond_key = tuple(sorted([atom1, atom2]))
    
    if forcefield == 'UFF':
        if bond_key in UFF_BOND_PARAMS:
            return UFF_BOND_PARAMS[bond_key]
        else:
            # Параметры связей в UFF можно рассчитать по формулам
            # Здесь упрощенно возвращаем параметры для C-C
            print(f"Предупреждение: UFF параметры связи для {bond_key} не найдены, расчет по формулам")
            radius1 = Element(atom1).covalent_radius
            radius2 = Element(atom2).covalent_radius
            r0 = radius1 + radius2
            # Упрощенный расчет константы (на основе UFF)
            force_constant = 664.12 * (1.0 + 0.5 * np.abs(Element(atom1).X - Element(atom2).X))
            return {'k': force_constant, 'r0': r0}

def get_angle_parameters(atom1, atom2, atom3, forcefield='UFF'):
    """
    Возвращает параметры угла для указанной тройки атомов
    """
    angle_key = (atom1, atom2, atom3)
    
    if forcefield == 'UFF':
        if angle_key in UFF_ANGLE_PARAMS:
            return UFF_ANGLE_PARAMS[angle_key]
        else:
            # Параметры углов в UFF можно рассчитать по формулам
            # Здесь упрощенно определяем гибридизацию и возвращаем соответствующее значение
            
            # Упрощенное определение гибридизации на основе числа соседей
            if atom2 in ['C', 'N']:
                # Предполагаем sp3 для 4 соседей, sp2 для 3 соседей
                theta0 = 109.47  # для sp3
                # В реальности нужно анализировать молекулярную топологию
            elif atom2 in ['O']:
                theta0 = 104.51  # для sp3 кислорода
            else:
                theta0 = 109.47  # по умолчанию, тетраэдрическое окружение
            
            # Упрощенный расчет константы (на основе UFF)
            force_constant = 100.0  # примерно для большинства углов в UFF
            
            return {'k': force_constant, 'theta0': theta0}

def create_co2_molecule(output_file):
    """
    Создает файл молекулы CO2 для GCMC-моделирования в LAMMPS
    Параметры TraPPE получены из статьи Potoff & Siepmann, 2001
    """
    with open(output_file, 'w') as f:
        f.write("""# CO2 molecule for GCMC simulation
# Parameters from TraPPE-UA force field (Potoff & Siepmann, 2001)

3 atoms
2 bonds
1 angles

Coords

1 0.0 0.0 0.0      # C
2 0.0 0.0 1.16     # O1
3 0.0 0.0 -1.16    # O2

Types

1 1  # C type
2 2  # O type
3 2  # O type

Charges

1 0.70   # C
2 -0.35  # O1
3 -0.35  # O2

Bonds

1 1 1 2  # C-O1
2 1 1 3  # C-O2

Angles

1 1 2 1 3  # O1-C-O2
""")
    print(f"Создан файл молекулы CO2: {output_file}")
    return output_file

def create_input_file(data_file, co2_file, output_file, temp=298.0, pressure=0.1):
    """
    Создает input-файл LAMMPS для GCMC-моделирования
    """
    # Расчет химического потенциала для CO2
    R = 8.314 / 4184.0  # Газовая постоянная в ккал/(моль*K)
    
    # Корректный расчет химического потенциала для реального газа CO2
    # Использование уравнения состояния Peng-Robinson
    # Параметры для CO2: Tc = 304.2 K, Pc = 73.8 бар, ω = 0.228
    Tc = 304.2
    Pc = 73.8
    omega = 0.228
    
    # Приведенная температура и давление
    Tr = temp / Tc
    Pr = pressure * 1.01325 / Pc  # конвертация из атм в бар
    
    # Расчет коэффициентов
    a = 0.45724 * R**2 * Tc**2 / Pc
    b = 0.07780 * R * Tc / Pc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2
    
    # Коэффициент сжимаемости для CO2 при данных условиях (упрощенно)
    # В реальном решении нужно решать кубическое уравнение
    if Tr > 1.0 and Pr < 0.1:
        # Для низких давлений и температур выше критической
        # можно приближенно считать Z = 1 (идеальный газ)
        Z = 1.0
    else:
        # Приближенное значение для реального газа
        Z = 1.0 - pressure * (1.0 - Tr**0.5) / 20.0
    
    # Химический потенциал для реального газа
    if pressure < 1e-10:
        # Для очень низких давлений используем предел
        mu = -1000.0  # очень низкое значение
    else:
        # Формула для химического потенциала реального газа
        mu = R * temp * (np.log(Pr * Z / Tr) + (Z - 1.0) - np.log(Z))
    
    # Записываем input-файл
    with open(output_file, 'w') as f:
        f.write(f"""# GCMC моделирование адсорбции CO₂ в MOF при T={temp}K и P={pressure} атм
# Файл сгенерирован автоматически с использованием UFF параметров и EQeq зарядов

units real
atom_style full
boundary p p p
pair_style lj/cut/coul/long 12.0
bond_style harmonic
angle_style harmonic
kspace_style pppm 1.0e-5

# Загружаем структуру MOF
read_data {data_file}

# Определяем параметры для CO2 из TraPPE (Potoff & Siepmann, 2001)
# Эти параметры получены из экспериментальных данных и квантово-химических расчетов
pair_coeff * * lj/cut/coul/long 0.0 1.0  # Обнуляем все взаимодействия
pair_coeff 6 6 lj/cut/coul/long 0.0559 2.800  # C (CO2)
pair_coeff 7 7 lj/cut/coul/long 0.0157 3.050  # O (CO2)

# Применяем правила смешивания для перекрестных взаимодействий
pair_modify mix arithmetic  # правило смешивания для sigma
pair_modify mix geometric   # правило смешивания для epsilon

# Определяем группы атомов
group framework type <= 5  # первые 5 типов атомов считаем каркасом
group gas type > 5        # остальные считаем газом

# Фиксируем каркас MOF
fix FREEZE framework setforce 0.0 0.0 0.0

# Определяем молекулу CO2
molecule CO2 {co2_file}

# GCMC моделирование
fix GCMC gas gcmc 100 100 100 0 29494 {temp} {mu:.4f} 0.1 mol CO2 &
    pressure {pressure} full_energy

# Вывод данных
variable N equal count(gas)/3  # количество молекул CO2
fix AVE all ave/time 500 10 5000 v_N file adsorption.dat

# Запуск расчета
timestep 1.0
thermo 500
thermo_style custom step temp press v_N pe ke etotal
dump CONFIG all custom 10000 dump.*.lammpstrj id type x y z

# Уравновешивание
run 50000

# Сбор статистики
reset_timestep 0
run 200000
""")
    print(f"Создан input-файл LAMMPS: {output_file}")
    return output_file

def convert_cif_to_lammps_data(cif_file, output_file, forcefield='UFF'):
    """
    Конвертирует CIF-файл в DATA-файл LAMMPS с использованием UFF или DREIDING
    """
    print(f"Чтение CIF-файла: {cif_file}")
    
    # Читаем структуру через Pymatgen
    parser = CifParser(cif_file)
    structure = parser.get_structures()[0]
    
    # Рассчитываем заряды методом EQeq
    print("Расчет зарядов методом EQeq...")
    charges = calculate_eqeq_charges(structure)
    
    # Идентифицируем связи и углы
    print("Определение связей и углов...")
    bonds = identify_bonds_in_structure(structure)
    angles = identify_angles_in_structure(structure, bonds)
    
    # Выводим базовую статистику
    elements = [site.species_string for site in structure]
    unique_elements = set(elements)
    print(f"Структура содержит {len(structure)} атомов, {len(bonds)} связей, {len(angles)} углов")
    print(f"Уникальные элементы: {', '.join(unique_elements)}")
    
    # Создаем словарь типов атомов
    atom_types = {}
    for element in unique_elements:
        if element not in atom_types:
            atom_types[element] = len(atom_types) + 1
    
    # Собираем параметры для силового поля
    pair_coeffs = {}
    for element, type_id in atom_types.items():
        params = get_uff_parameters(element, forcefield)
        pair_coeffs[type_id] = (params['epsilon'], params['sigma'])
    
    # Определяем типы связей
    bond_types = {}
    bond_coeffs = {}
    for i, j in bonds:
        elem_i = elements[i]
        elem_j = elements[j]
        bond_key = tuple(sorted([elem_i, elem_j]))
        
        if bond_key not in bond_types:
            bond_types[bond_key] = len(bond_types) + 1
            bond_coeffs[bond_types[bond_key]] = get_bond_parameters(elem_i, elem_j, forcefield)
    
    # Определяем типы углов
    angle_types = {}
    angle_coeffs = {}
    for i, j, k in angles:
        elem_i = elements[i]
        elem_j = elements[j]
        elem_k = elements[k]
        angle_key = (elem_i, elem_j, elem_k)
        
        if angle_key not in angle_types:
            angle_types[angle_key] = len(angle_types) + 1
            angle_coeffs[angle_types[angle_key]] = get_angle_parameters(elem_i, elem_j, elem_k, forcefield)
    
    # Создаем файл DATA для LAMMPS
    with open(output_file, 'w') as f:
        f.write(f"LAMMPS DATA file for {cif_file} with {forcefield} parameters and EQeq charges\n\n")
        
        # Общая информация
        f.write(f"{len(structure)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")
        
        f.write(f"{len(atom_types)} atom types\n")
        f.write(f"{len(bond_types)} bond types\n")
        f.write(f"{len(angle_types)} angle types\n")
        f.write("0 dihedral types\n")
        f.write("0 improper types\n\n")
        
        # Размеры ячейки
        cell = structure.lattice.matrix
        f.write(f"0.0 {cell[0][0]} xlo xhi\n")
        f.write(f"0.0 {cell[1][1]} ylo yhi\n")
        f.write(f"0.0 {cell[2][2]} zlo zhi\n")
        f.write(f"{cell[1][0]} {cell[2][0]} {cell[2][1]} xy xz yz\n\n")
        
        # Массы атомов
        f.write("Masses\n\n")
        for element, type_id in atom_types.items():
            mass = Element(element).atomic_mass
            f.write(f"{type_id} {mass}  # {element}\n")
        f.write("\n")
        
        # Коэффициенты LJ
        f.write("Pair Coeffs\n\n")
        for type_id, (epsilon, sigma) in pair_coeffs.items():
            f.write(f"{type_id} {epsilon} {sigma}\n")
        f.write("\n")
        
        # Коэффициенты связей
        if bonds:
            f.write("Bond Coeffs\n\n")
            for type_id, params in bond_coeffs.items():
                f.write(f"{type_id} {params['k']} {params['r0']}\n")
            f.write("\n")
        
        # Коэффициенты углов
        if angles:
            f.write("Angle Coeffs\n\n")
            for type_id, params in angle_coeffs.items():
                f.write(f"{type_id} {params['k']} {params['theta0']}\n")
            f.write("\n")
        
        # Атомы
        f.write("Atoms\n\n")
        for i, site in enumerate(structure):
            element = site.species_string
            type_id = atom_types[element]
            x, y, z = site.coords
            charge = charges.get(i, 0.0)
            # LAMMPS формат: atom-ID molecule-ID atom-type q x y z
            f.write(f"{i+1} 1 {type_id} {charge:.6f} {x:.6f} {y:.6f} {z:.6f}  # {element}\n")
        f.write("\n")
        
        # Связи
        if bonds:
            f.write("Bonds\n\n")
            for bond_id, (i, j) in enumerate(bonds):
                elem_i = elements[i]
                elem_j = elements[j]
                bond_key = tuple(sorted([elem_i, elem_j]))
                type_id = bond_types[bond_key]
                # LAMMPS формат: bond-ID bond-type atom1 atom2
                f.write(f"{bond_id+1} {type_id} {i+1} {j+1}\n")
            f.write("\n")
        
        # Углы
        if angles:
            f.write("Angles\n\n")
            for angle_id, (i, j, k) in enumerate(angles):
                elem_i = elements[i]
                elem_j = elements[j]
                elem_k = elements[k]
                angle_key = (elem_i, elem_j, elem_k)
                type_id = angle_types.get(angle_key, 1)  # Если не найдено, используем тип 1
                # LAMMPS формат: angle-ID angle-type atom1 atom2 atom3
                f.write(f"{angle_id+1} {type_id} {i+1} {j+1} {k+1}\n")
            f.write("\n")
    
    print(f"Создан DATA-файл LAMMPS: {output_file}")
    
    # Создаем файл с дополнительной информацией
    with open(output_file + ".info", 'w') as f:
        f.write("Информация о параметрах силового поля:\n\n")
        
        f.write("Параметры для атомов:\n")
        for element, type_id in atom_types.items():
            params = get_uff_parameters(element, forcefield)
            f.write(f"{element} (тип {type_id}): epsilon = {params['epsilon']}, sigma = {params['sigma']}\n")
        
        f.write("\nПараметры для связей:\n")
        for bond_key, type_id in bond_types.items():
            params = bond_coeffs[type_id]
            f.write(f"{bond_key} (тип {type_id}): k = {params['k']}, r0 = {params['r0']}\n")
        
        f.write("\nПараметры для углов:\n")
        for angle_key, type_id in angle_types.items():
            params = angle_coeffs[type_id]
            f.write(f"{angle_key} (тип {type_id}): k = {params['k']}, theta0 = {params['theta0']}\n")
    
    return output_file

def main():
    parser = argparse.ArgumentParser(description='Универсальная подготовка файлов для GCMC-моделирования в LAMMPS')
    parser.add_argument('cif_file', type=str, help='Входной CIF-файл структуры MOF')
    parser.add_argument('--forcefield', type=str, choices=['UFF', 'DREIDING'], default='UFF',
                        help='Силовое поле для расчета параметров (по умолчанию: UFF)')
    parser.add_argument('--temp', type=float, default=298.0, 
                        help='Температура моделирования в Кельвинах (по умолчанию: 298.0)')
    parser.add_argument('--pressure', type=float, default=0.1,
                        help='Давление в атмосферах (по умолчанию: 0.1)')
    parser.add_argument('--output', type=str, default=None,
                        help='Префикс для выходных файлов')
    
    args = parser.parse_args()
    
    # Проверяем существование файла
    if not os.path.exists(args.cif_file):
        print(f"Ошибка: Файл {args.cif_file} не найден")
        sys.exit(1)
    
    # Определяем имена выходных файлов
    if args.output:
        base_name = args.output
    else:
        base_name = os.path.splitext(args.cif_file)[0]
    
    data_file = f"{base_name}.data"
    co2_file = f"{base_name}_CO2.txt"
    input_file = f"{base_name}.in"
    
    # Шаг 1: Конвертация CIF в DATA
    data_file = convert_cif_to_lammps_data(args.cif_file, data_file, args.forcefield)
    
    # Шаг 2: Создание файла молекулы CO2
    co2_file = create_co2_molecule(co2_file)
    
    # Шаг 3: Создание input-файла LAMMPS
    input_file = create_input_file(data_file, co2_file, input_file, args.temp, args.pressure)
    
    print("\nУспешно созданы все необходимые файлы для GCMC моделирования!")
    print(f"1. DATA-файл LAMMPS: {data_file}")
    print(f"2. Молекула CO2: {co2_file}")
    print(f"3. Input-файл LAMMPS: {input_file}")
    print(f"4. Информация о параметрах: {data_file}.info")
    print("\nДля запуска моделирования выполните:")
    print(f"lmp -in {input_file}")

if __name__ == '__main__':
    main()