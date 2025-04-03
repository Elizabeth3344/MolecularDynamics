All modeling projects are done in folder "TIP4P-Lammps".

1) Create a file named `in.<file_name>` with your Lammps script. If you mention another auxiliary code files in your main one, adress them as `../../TIP4P-Lammps/in.<file_name>`. `../../` helps you to get out of lmp, which runs your script.
2) To run your script in terminal change your current folder to "lammps/build" and go `./lmp -in ../../TIP4P-Lammps/in.<file_name>`. You may find "lmp" in the "build" folder. You can do `lmp` from only this location in "lammps".

There may be some problems during the first launches:
1) "ERROR: Unrecognized fix style 'rigid/small' is part of the RIGID package which is not enabled in this LAMMPS binary".
This error means that you need to install the rigid package for your lammps. To install the package from "lammps/build" go `cmake -D PKG_RIGID=yes ../cmake` and then `cmake --build .`.
2) Rigid may also need mpi-cxx package for parallel computing. The error "No package 'mpi-cxx' found" may appear during the trying to go `cmake ../cmake`. To solve the problem go `sudo apt install mpich` from your any current path. Now mpi is located in "lammps/src". You may check it's presents by command `ls | grep lmp_mpi` in "src".
3) You may not to pay your attention on the messages such as "Could NOT find ClangFormat". This tool is needed for formatting and doesn't affect your codework.

Running with mpi:
1) Only once from home go to hidden file .bashrc, with the command `vim .bashrc` or `nano .bashrc` open it for additing (in `vim` use `insert`) and write at the end of the file: `export LD_LIBRARY_PATH=/home/lammps/src:$LD_LIBRARY_PATH` and `export PATH=/home/lammps/src:$PATH`. It is necessary for being able to use mpi from any of your location.
2) Now to run the lammps script go in terminal from "TIP4P-Lammps" folder: `mpirun -np 4 -in in.<file_name>`.

Use ovito:
`Downloads/ovito-basic-3.10.6-x86_64/bin/ovito`.
