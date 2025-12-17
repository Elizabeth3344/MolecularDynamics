#! /bin/bash

# Проверка
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Использование: ./run.sh <имя_файла_в_папке_MOF> <газ>"
    echo "Пример: ./run.sh MOF5.cif CO2"
    exit 1
fi

MOF_NAME=$1
GAS_NAME=$2
MOF_PATH="MOF/$MOF_NAME"

if [ ! -f "$MOF_PATH" ]; then
    echo "Ошибка: Файл '$MOF_NAME' не найден в папке 'MOF/'"
    exit 1
fi

# Проверяем, существует ли шаблон
if [ ! -f "simulation.input.template" ]; then
    echo "Ошибка: Не найден файл-шаблон 'simulation.input.template'!"
    exit 1
fi

# Пути
export RASPA_DIR=${HOME}/RASPA/simulations/
export DYLD_LIBRARY_PATH=${RASPA_DIR}/lib
export LD_LIBRARY_PATH=${RASPA_DIR}/lib

# 1. КОПИРОВАНИЕ ДЕФОВ
cp ForceField/*.def . 2>/dev/null
if [ ! -f "force_field_mixing_rules.def" ]; then
    echo "Ошибка: В папке ForceField/ нет .def файлов!"
    exit 1
fi

# 2. ПОДГОТОВКА (Python)
python3 convert_to_raspa.py "$MOF_PATH" "$GAS_NAME"

# 3. СОЗДАНИЕ INPUT ФАЙЛА ИЗ ШАБЛОНА
# Берем чистый шаблон -> создаем рабочий файл simulation.input
cp simulation.input.template simulation.input

# Меняем REPLACE_ME на имя газа ТОЛЬКО в рабочем файле
# Проверяем ОС
if [[ "$OSTYPE" == "darwin"* ]]; then
    # Это macOS
    sed -i '' "s/REPLACE_ME/$GAS_NAME/g" simulation.input
else
    # Это Linux
    sed -i "s/REPLACE_ME/$GAS_NAME/g" simulation.input
fi

# 4. ЗАПУСК
echo ">>> Запуск симуляции ($GAS_NAME)..."
$RASPA_DIR/bin/simulate simulation.input

# 5. УБОРКА
rm *.def
echo ">>> Готово."