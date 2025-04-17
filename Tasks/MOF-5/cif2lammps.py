from ase.io import read, write
mof5 = read("MOF5.cif")
write("MOF5.data", mof5, format="lammps-data")

# from pymatgen.io.cif import CifParser
# from pymatgen.io.lammps.data import LammpsData

# # Загрузка CIF файла
# parser = CifParser('MOF-5.cif')
# structure = parser.get_structures()[0]

# # Создание объекта LammpsData
# lammps_data = LammpsData.from_structure(structure, atom_style='full')

# # Запись в файл
# lammps_data.write_file('MOF-5.data')
