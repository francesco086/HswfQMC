# HswfQMC: Hydrogen Shadow Wave Function Quantum Monte Carlo


##DOWNLOAD:

Simply execute the command

  `git clone https://github.com/francesco086/HswfQMC`


##INSTALL:

Follow these simple steps:

* Set the PATH. To do so, execute the command (from the root folder of HswfQMC)

  `sh pilot-HswfQMC.sh setpath`
  
  `source ~/.bashrc` (or `source ~/.bash_profile` if you are using a Mac)

  After running this command you can use the script `pilot-HswfQMC` (in the following we will refer to it simply as `pilot`) as an executable from anywhere. To give it a try and see which options you have, go to your home directory (`cd`) and run
  
  `pilots-HswfQMC help`

* Install the Lapack library. We reccomend to use `pilot` in order to have the latest, fully optimized version of the Lapack library. Run

  `pilot-HswfQMC install_lapack``

  and follow the instructions. If you prefer using the system library or your own version, you can skip this step.
  
* Set the Makefile. The Makefile is located in source/Makefile. We reccomend to set it with `pilot` by running
  
 `pilots-HswfQMC setmakefile``

  If you prefer to set it by hand you can simply edit it by yourself.

* Compile the HswfQMC code. This can be done with `pilot` by running:

  `pilot-HswfQMC build``

  If you prefer to do it manually, just run make inside the source/ folder.


##RUN:

`cd`to the folder where you want to run your simulation and execute

  `pilot-HswfQMC.sh setdir`

Then set the input files "dati_*.d" and finally start the simulation by

  `mpirun HswfQMC_exe`  (add the option `-np N` to use multiple cores)

