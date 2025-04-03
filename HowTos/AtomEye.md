AtomEye is used for visualisation your Lammps models.<br />

Installation AtomEye:<br />
1) Follow the link: http://li.mit.edu/A/Graphics/A/#download and download the file under the name "i686 Linux" for Linux. You may find downloaded file in "Downloads".
2) To grant launch rights, write: ``sudo chmod +x ./Downloads/A.i686``.
3) Now you need to install the X terminal for some input / output as well. For this go ``sudo apt update`` and ``sudo apt install xterm``.
4) You may download this file of graphene nanotube for check: [test.cfg](Test/test.cfg). To open this file with Atomeye go: ``./Downloads/A.i686 ./test.cfg``. Do not forget to specify the correct path to this file on your computer.
5) Now you can use Atomeye only from this path "./Downloads/A.i686". To be able to use it from any location go: ``sudo mv ./Downloads/A.i686 /usr/bin/``. This will move Atomeye into the standard paths.
6) This item is optional. You may change the name of the command which helps you to run AtomEye from "A.i686" to "Atomeye": ``sudo mv /usr/bin/A.i686 /usr/bin/Atomeye``.
7) Now use AtomEye with command: ``Atomeye file_location/file_name.cfg``.

To run AtomEye from every place:<br />
``Atomeye file_direction/file_name.cfg`` or ``.pdb``.
