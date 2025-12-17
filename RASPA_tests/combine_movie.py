import os
import sys
import glob

def combine_raspa_output():
    # Путь к папке с фильмами
    movie_dir = os.path.join("Movies", "System_0")
    
    # 1. Ищем файл каркаса (Framework) .pdb
    framework_files = glob.glob(os.path.join(movie_dir, "Framework_final.pdb"))
    if not framework_files:
        print("Ошибка: Не найден .pdb файл каркаса в папке Movies/System_0")
        print("Убедитесь, что в simulation.input стоит 'Movies yes'")
        return
    framework_file = framework_files[0]

    # 2. Ищем файл видео газа (Movie) .pdb
    movie_files = glob.glob(os.path.join(movie_dir, "Movie_simulation_structure_*.pdb"))
    if not movie_files:
        print("Ошибка: Не найден .pdb файл газа в папке Movies/System_0")
        return
    movie_file = movie_files[0] # Берем первый найденный компонент

    output_file = "FULL_MOVIE.pdb"

    print(f"Каркас: {os.path.basename(framework_file)}")
    print(f"Газ:    {os.path.basename(movie_file)}")
    print("Объединение...")

    # Читаем атомы каркаса
    framework_atoms = []
    with open(framework_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                framework_atoms.append(line)

    # Читаем видео газа и внедряем каркас в каждый кадр
    with open(movie_file, 'r') as f_in, open(output_file, 'w') as f_out:
        frame_count = 0
        for line in f_in:
            if line.startswith("MODEL"):
                f_out.write(line) # Пишем начало кадра
                # Вставляем статичный каркас в этот кадр
                for atom in framework_atoms:
                    f_out.write(atom)
                frame_count += 1
            else:
                f_out.write(line)
    
    print(f"Готово! Создан файл: {output_file}")
    print(f"Обработано кадров: {frame_count}")
    print("Перетащите FULL_MOVIE.pdb в Ovito.")

if __name__ == "__main__":
    combine_raspa_output()