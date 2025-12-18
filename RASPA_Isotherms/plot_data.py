import pandas as pd
import matplotlib.pyplot as plt
import sys

# Газовая постоянная для перевода K -> kJ/mol
R = 8.314462618

def analyze_results():
    csv_file = "final_results.csv"
    
    try:
        data = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"ОШИБКА: Файл {csv_file} не найден. Сначала запустите main_experiment.py")
        sys.exit(1)

    # 1. Пересчет теплоты (Heat of Adsorption)
    # RASPA выдает <U_gh> - <U_h> в Кельвинах.
    # Q_st (kJ/mol) = - (Heat_K * R) / 1000 + RT (примерно)
    # Но обычно для оценки берут просто конвертацию энергии.
    # Если в CSV записано положительное число (Heat of Desorption), то это и есть Qst.
    
    # Добавляем колонку в кДж/моль
    if "Heat_K" in data.columns:
        data["Qst_kJ"] = data["Heat_K"] * R / 1000.0
    
    print("\n--- Результаты (Теплота адсорбции) ---")
    # Выводим среднюю теплоту для каждой температуры (обычно она падает с давлением)
    for temp in data["T"].unique():
        subset = data[data["T"] == temp]
        avg_q = subset["Qst_kJ"].mean()
        print(f"T = {temp} K: Средняя Qst ≈ {avg_q:.2f} кДж/моль")

    # 2. Построение Изотерм
    plt.figure(figsize=(10, 6))
    
    # Получаем список температур
    temps = sorted(data["T"].unique())
    
    colors = ['b', 'r', 'g', 'k', 'm'] # Цвета линий
    
    for i, temp in enumerate(temps):
        subset = data[data["T"] == temp].sort_values("P")
        
        # P - давление (бар)
        # Loading - выберем Абсолютную (можно поменять на Exc_mol_kg)
        plt.plot(
            subset["P"], 
            subset["Abs_mol_kg"], 
            marker='o', 
            linestyle='-', 
            linewidth=2,
            label=f"{temp} K",
            color=colors[i % len(colors)]
        )

    plt.title("Изотермы адсорбции", fontsize=14)
    plt.xlabel("Давление (бар)", fontsize=12)
    plt.ylabel("Абсолютная адсорбция (моль/кг)", fontsize=12)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.legend(title="Температура")
    
    # Сохранение
    plt.savefig("isotherms.png", dpi=300)
    print(f"\n--- График сохранен как 'isotherms.png' ---")
    plt.show()

if __name__ == "__main__":
    # Проверка на наличие библиотек
    try:
        import pandas
        import matplotlib
    except ImportError:
        print("Для работы нужны pandas и matplotlib.")
        print("pip install pandas matplotlib")
        sys.exit(1)
        
    analyze_results()