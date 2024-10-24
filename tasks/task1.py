import sys
from lammps import lammps

# Создаем экземпляр LAMMPS с выводом в консоль
lmp = lammps(cmdargs=["-echo", "both"])

# Команды LAMMPS для задания кубической ячейки ГЦК меди и растяжения
lmp.commands_string("""
# 1. Инициализация
units metal
dimension 3
boundary p p p
atom_style atomic

# 2. Создание атомов и задания начальной структуры меди (ГЦК)
lattice fcc 3.615
region box block 0 10 0 10 0 10
create_box 1 box
create_atoms 1 box

# 3. Задание взаимодействий через потенциал внедренного атома
pair_style eam
pair_coeff * * Cu_u3.eam

# 4. Определение вычислений и условий
mass 1 63.546
velocity all create 300.0 12345 mom yes rot yes

# 5. Управление шагами симуляции и вывод результатов
thermo 100
thermo_style custom step temp pe ke etotal press
dump 1 all atom 100 dump.cube_stretch2
dump_modify 1 scale yes

# 6. Фиксация границ и растяжение куба по оси z
fix 1 all nvt temp 300.0 300.0 100.0
fix 2 all deform 1 z erate 0.001 units box

# 7. Запуск симуляции
run 5000
""")

# Завершение работы
lmp.close()


if __name__ == '__main__':
    sys.exit()
