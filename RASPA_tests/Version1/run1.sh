#! /bin/sh -f

# Проверяем, подали ли имя файла
if [ -z "$1" ]; then
    echo "Usage: ./run.sh <path_to_cif_file>"
    exit 1
fi

INPUT_CIF=$1

# 1. Запускаем питон-конвертер
# Он берет ваш INPUT_CIF и создает красивый 'simulation_structure.cif'
echo "Converting $INPUT_CIF to RASPA format..."
python3 convert_to_raspa.py "$INPUT_CIF" simulation_structure

# Проверяем, создался ли файл
if [ ! -f "simulation_structure.cif" ]; then
    echo "Error: conversion failed."
    exit 1
fi

# 2. Обязательно: меняем имя фреймворка в simulation.input на то, которое мы создали
# (simulation_structure), если оно вдруг другое.
# Но проще просто в simulation.input один раз написать:
# FrameworkName simulation_structure

# 3. Запускаем RASPA
echo "Running simulation..."

# На всякий случай экспортируем пути (хотя при Local они менее важны, но бинарник нужен)
export RASPA_DIR=${HOME}/RASPA/simulations/
export DYLD_LIBRARY_PATH=${RASPA_DIR}/lib
export LD_LIBRARY_PATH=${RASPA_DIR}/lib

$RASPA_DIR/bin/simulate simulation.input