from lammps import lammps


def relax_box(type_lattice, parameter_lattice, size, potential_file, step, mass_atom, temp):
    # Создаем экземпляр LAMMPS с выводом в консоль
    lmp = lammps(cmdargs=["-echo", "both"])

    # Команды LAMMPS для задания кубической ячейки ГЦК меди и растяжения
    lmp.commands_string(f"""
    # 1. Инициализация
    units metal
    dimension 3
    boundary p p p
    atom_style atomic
    
    # 2. Создание атомов и задания начальной структуры меди (ГЦК)
    lattice {type_lattice} {parameter_lattice}
    region box block 0 {size} 0 {size} 0 {size}
    create_box 1 box
    create_atoms 1 box
    
    # 3. Задание взаимодействий через потенциал внедренного атома
    pair_style eam
    pair_coeff * * {potential_file}
    
    # 4. Определение вычислений и условий
    mass 1 {mass_atom}
    velocity all create {temp} 12345 mom yes rot yes
    
    fix 1 all nvt temp {temp} {temp} 0.05 
    thermo 1000  
    thermo_style custom step temp etotal pe press vol  
    timestep 0.001
    run {step}  
    write_restart relaxed_box.restart """)

    # Завершение работы
    lmp.close()

