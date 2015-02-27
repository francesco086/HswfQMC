#compile for debugging (useful for finding arrays overflow)
DEBUG=false
#compile without optimization or debugging flags, for quick compiling
QUICK_COMPILE=false
#standard executable name
EXEC2=HswfQMC-deb

COMPUTER:=$(shell hostname)
USERNAME:=$(shell whoami)
IDENTIFIER:=$(USERNAME)@$(COMPUTER)

#set the flags for the compiler 
FFLAGS=-O3
#standard: append _exe to the executable name
EXEC1=$(EXEC2)_exe

#quick compile compiling mode
ifeq ($(QUICK_COMPILE),true)
	FFLAGS=-O0
	EXEC1=$(EXEC2)_qc
endif

#debug compiling mode
ifeq ($(DEBUG),true)
  	FFLAGS=-O0 -g -fbounds-check
	EXEC1=$(EXEC2)_debug
endif



#set the libraries and machine-depending variables

#stored settings for users generated with pilotHswfQMC setmakefile
include makefile.users_settings

#by hand
ifeq ($(IDENTIFIER),username@hostname)    #insert your username and hostname
	EXEC=             #Name of the executable. Standard is: EXEC=$(EXEC1)
	LIBS=             #necessary library, BLAS AND LAPACK. Therefore usually: LIBS=-llapack -lblas
	LDFLAGS=          #path to the libraries, if needed. With system libraries is empty.
	FC=               #fortran compiler. Normally: FC=gfortran
endif


#OLD SETTINGS
#ifeq ($(COMPUTER),macbook)
#	EXEC=$(EXEC1)_macbook
#	# elenco librerie dinamiche (quelle che aggiungi tu a mano per l'eseguibile in questione
#	# che non centrano con quelle di sistema statiche)
#   LIBS=-lrefblas -llapack
#	# questo è per includere il percorso locale dove hai le tue lapack o quantaltro compilato
#	# correttamente; se c'è una versione di sistema delle lapack non serve esplicitarlo
#   LDFLAGS=-L/Users/kenzo/lib/lapack-3.5.0/
#   FC=mpif90-openmpi-gcc49
#endif
#ifeq ($(COMPUTER),login01.mogon)
#	EXEC=$(EXEC1)
#	FC=mpif90
#    	LIBS=-lmylapack -lmyblas
#    	LDFLAGS=-L/home/calcavec/lapack-3.4.2/
#endif


# elenco degli oggetti da creare (c'è sempre l'oggetto del programma principale e gli oggetti
# dei vari moduli)
OBJECTS=main.o module_dft.o module_dati_fisici.o module_dati_mc.o module_dnfH.o module_stati_eccitati.o module_generic_tools.o module_lattice.o module_VMC.o module_estimatori.o module_momenta.o module_calcola_accettazione.o module_funzione_onda.o module_walkers.o module_variational_opt.o module_variational_calculations.o


# crea tutti gli oggetti e li linka correttamente con le librerie; scrivere "make"
all: $(OBJECTS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

# rimuove tutti gli oggetti e moduli creati lasciando solo i sorgenti; "make clean"
.PHONY : clean
clean:
	rm -f *\.mod *\.o *~ HswfQMC*

# crea singolarmente i vari oggetti (se per caso vuoi crearne solo alcuni a mano); "make nome"
# attenzione a inserire tutti i moduli da cui ogni modulo può dipendere
# la sintassi vuole il tab nella seconda riga
main.o: module_dft.o module_stati_eccitati.o module_walkers.o module_dnfH.o module_funzione_onda.o module_calcola_accettazione.o module_estimatori.o module_generic_tools.o module_lattice.o module_dati_mc.o module_dati_fisici.o module_VMC.o module_variational_opt.o main.f90 
	$(FC) $(FFLAGS) -c main.f90 

module_variational_calculations.o: module_walkers.o module_dati_mc.o module_dati_fisici.o module_funzione_onda.o module_momenta.o module_estimatori.o module_calcola_accettazione.o module_estimatori.o module_variational_calculations.f90
	$(FC) $(FFLAGS) -c module_variational_calculations.f90 

module_VMC.o: module_stati_eccitati.o module_walkers.o module_calcola_accettazione.o module_estimatori.o module_generic_tools.o module_dati_mc.o module_dati_fisici.o module_momenta.o module_variational_calculations.o module_VMC.f90
	$(FC) $(FFLAGS) -c module_VMC.f90 

module_estimatori.o: module_calcola_accettazione.o module_momenta.o module_funzione_onda.o module_dati_mc.o module_dati_fisici.o module_walkers.o module_generic_tools.o module_estimatori.f90
	$(FC) $(FFLAGS) -c module_estimatori.f90

module_calcola_accettazione.o: module_walkers.o module_dati_mc.o module_dati_fisici.o module_funzione_onda.o module_momenta.o module_calcola_accettazione.f90
	$(FC) $(FFLAGS) -c module_calcola_accettazione.f90

module_walkers.o: module_generic_tools.o module_dati_fisici.o module_dati_mc.o module_walkers.f90
	$(FC) $(FFLAGS) -c module_walkers.f90 

module_lattice.o: module_lattice.f90
	$(FC) $(FFLAGS) -c module_lattice.f90

module_generic_tools.o: module_generic_tools.f90
	$(FC) $(FFLAGS) -c module_generic_tools.f90 

module_dati_mc.o: module_dati_fisici.o module_generic_tools.o module_dati_mc.f90
	$(FC) $(FFLAGS) -c module_dati_mc.f90

module_dati_fisici.o: module_lattice.o module_generic_tools.o module_dati_fisici.f90
	$(FC) $(FFLAGS) -c module_dati_fisici.f90 

module_dnfH.o: module_dati_fisici.o module_walkers.o module_dnfH.f90
	$(FC) $(FFLAGS) -c module_dnfH.f90 

module_funzione_onda.o: module_dft.o module_momenta.o module_generic_tools.o module_dati_mc.o module_dati_fisici.o module_dnfH.o module_walkers.o module_funzione_onda.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -c module_funzione_onda.f90 $(LIBS)

module_momenta.o: module_momenta.f90
	$(FC) $(FFLAGS) -c module_momenta.f90 

module_variational_opt.o: module_walkers.o module_VMC.o module_funzione_onda.o module_dati_mc.o module_dati_fisici.o module_variational_calculations.o module_generic_tools.o module_variational_opt.f90
	$(FC) $(FFLAGS) -c module_variational_opt.f90 

module_stati_eccitati.o: module_dati_mc.o module_dati_fisici.o module_funzione_onda.o module_walkers.o module_stati_eccitati.f90
	$(FC) $(FFLAGS) -c module_stati_eccitati.f90 

module_dft.o: module_dft.f90
	$(FC) $(FFLAGS) -c module_dft.f90
