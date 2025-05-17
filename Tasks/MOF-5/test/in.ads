units           real
atom_style      full
boundary        p p p

pair_style      lj/cut 12.500
bond_style      harmonic
angle_style     hybrid fourier cosine/periodic
dihedral_style  harmonic
improper_style  fourier

pair_modify     tail yes mix arithmetic
special_bonds   lj/coul 0.0 0.0 1.0
dielectric      1.0
box tilt        large
read_data       data.MOF5

# ******* molecules ******
molecule co2_mol CO2.txt 

group           mof     type 1 2 3 4 5
group           co2     type 6 7

compute         co2T co2 temp
compute_modify  co2T dynamic/dof yes
compute         mofT mof temp

variable N_co2  equal count(co2)/3

# timestep 0.1

dump dco2 co2 atom 100 dump_co2.lammpstrj
dump dall all atom 100 dump_all.lammpstrj

fix             NPT all npt temp 323.0 323.0 100.0 iso 10.0 10.0 1000.0
thermo          10
thermo_style    custom step temp press density
run             10000

fix GCMC co2 gcmc 10 100 0 0 1263 323.0 -5 0.5 mol co2_mol pressure 10.0 &
    full_energy tfac_insert 1.6 group co2 overlap_cutoff 1.5 
    
fix AVE all ave/time 10 100 1000 v_N_co2
thermo          10
thermo_style    custom step temp press density atoms v_N_co2
run            100000