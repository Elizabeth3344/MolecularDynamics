import os
import sys
import glob

def make_cloud():
    base_results = "Results"
    if not os.path.exists(base_results):
        print("Ошибка: Папка Results не найдена.")
        return

    print("--- Генерация облака точек (Binding Sites) из Видео ---")

    # 1. Ищем папку с результатами (где есть Movies)
    target_dir = None
    movie_dir = None
    
    # Бегаем по папкам, ищем ту, где лежит Movie файл (он должен быть не пустой)
    for root, dirs, files in os.walk(base_results):
        if "Movies" in dirs:
            check_path = os.path.join(root, "Movies", "System_0")
            # Ищем файл с газом (Movie...pdb)
            movies = glob.glob(os.path.join(check_path, "Movie_simulation_*.pdb"))
            if movies:
                target_dir = root
                movie_dir = check_path
                movie_file = movies[0]
                break
    
    if not target_dir:
        print("ОШИБКА: Не найдено папок с записанными видео (Movies).")
        print("Убедитесь, что для одной из точек (0.1 bar) стояло 'Movies yes'.")
        return

    # Ищем файл каркаса там же
    fw_files = glob.glob(os.path.join(movie_dir, "Framework_*_final.pdb"))
    if not fw_files:
        print("ОШИБКА: Не найден файл каркаса (Framework...pdb).")
        return
    framework_file = fw_files[0]

    output_file = "BINDING_SITES_CLOUD.pdb"
    print(f"Обработка папки: {target_dir}")
    print(f"Видео газа: {os.path.basename(movie_file)}")
    print("Создание облака...")

    atom_count = 0
    
    with open(output_file, 'w') as f_out:
        # 1. Записываем каркас (он статичен)
        with open(framework_file, 'r') as f_in:
            for line in f_in:
                if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("CRYST"):
                    f_out.write(line)

        # 2. Записываем ВСЕ кадры газа как одно большое облако
        with open(movie_file, 'r') as f_in:
            for line in f_in:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Меняем название атома, чтобы в Ovito можно было покрасить отдельно
                    # (В PDB колонки фиксированы, меняем символы 17-20 на 'GAS ')
                    # Но проще просто записать как есть, Ovito разберется по типу
                    
                    # Чтобы газ отличался от каркаса, мы можем хакнуть название остатка (Residue Name)
                    # Обычно там "MOL". Мы заменим на "GAS".
                    new_line = line.replace("MOL", "GAS")
                    f_out.write(new_line)
                    atom_count += 1

    print(f"\nГОТОВО! Файл создан: {output_file}")
    print(f"В облаке: {atom_count} точек.")
    print("--- Инструкция для Ovito ---")
    print("1. Откройте BINDING_SITES_CLOUD.pdb")
    print("2. В Pipeline Browser -> Particles -> Particle Types:")
    print("   - Найдите типы атомов газа (обычно C_co2, O_co2).")
    print("   - Сделайте их ЯРКИМИ, увеличьте Radius и добавьте прозрачность (Transparency).")
    print("   - Атомы каркаса сделайте мелкими (0.3).")
    print("3. Вы увидите 'густые' зоны скопления молекул.")

if __name__ == "__main__":
    make_cloud()