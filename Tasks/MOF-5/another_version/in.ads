######################   SYSTEM CREATION   #################################
units           real
atom_style      full
boundary        p p p

pair_style      lj/cut/coul/long 12
bond_style      harmonic
angle_style     hybrid cosine/squared harmonic
dihedral_style  charmm
improper_style  umbrella

special_bonds   dreiding
dielectric      1.0
box tilt        large
read_data       /home/chernysheva/MolecularDynamics/Tasks/MOF-5/another_version/thermo.data
# read_data data.MOF5

neigh_modify    every 1 delay 0 check yes

kspace_style    pppm 1e-6
pair_modify     tail yes mix arithmetic

#******* molecules ******
molecule co2_mol /home/chernysheva/MolecularDynamics/Tasks/MOF-5/another_version/CO2.txt

############################   SIMULATION   #################################
group           mof     type 1 2 3 4 5
group           co2     type 6 7

compute         co2T co2 temp
compute_modify  co2T dynamic/dof yes
compute         mofT mof temp

variable N_co2  equal count(co2)/3

timestep 0.1

dump dco2 co2 atom 100 dump_co2.lammpstrj
dump dall all atom 100 dump_all.lammpstrj

fix NVT all nvt temp 273.0 273.0 100.0
# fix below works with the MPI
# because this not allow translations (only вставки и удаления)
# fix GCMC co2 gcmc 10 100 *0*  0 1263 298 -16.2 0.5 mol co2_mol full_energy tfac_insert 1.6 group co2 overlap_cutoff 1.5 

fix GCMC co2 gcmc 10 100 0  0 1263 273.0 -16.2 0.5 mol co2_mol pressure 10.0 full_energy tfac_insert 1.6 group co2 overlap_cutoff 1.5
fix AVE all ave/time 10 100 1000 v_N_co2

thermo          10
thermo_style    custom step temp press density atoms v_N_co2 f_GCMC[1] f_GCMC[2] 

run            100000