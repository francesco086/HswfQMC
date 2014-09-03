DEBUG=false
QUICK_COMPILE=false
EXEC2=HswfQMC
COMPUTER=$(COMPUTER_MACHINE)

ifeq ($(DEBUG),true)
  	FFLAGS=-O0 -g -fbounds-check
	EXEC1=$(EXEC2)_debug
else
	FFLAGS=-O3
	EXEC1=$(EXEC2)
endif
ifeq ($(QUICK_COMPILE),true)
	FFLAGS=-O0
	EXEC1=$(EXEC2)_qc
endif

ifeq ($(COMPUTER),macbook)
	EXEC=$(EXEC1)_macbook
	# elenco librerie dinamiche (quelle che aggiungi tu a mano per l'eseguibile in questione
	# che non centrano con quelle di sistema statiche)
   LIBS=-lrefblas -llapack
	# questo è per includere il percorso locale dove hai le tue lapack o quantaltro compilato
	# correttamente; se c'è una versione di sistema delle lapack non serve esplicitarlo
   LDFLAGS=-L/Users/kenzo/lib/lapack-3.5.0/
   FC=mpif90-openmpi-gcc49
endif
ifeq ($(COMPUTER),hall)
        EXEC=$(EXEC1)_hall
        LIBS=-lblas -llapack
        FC=mpif90
endif
ifeq ($(COMPUTER),langevin)
	EXEC=$(EXEC1)_langevin
	LIBS=-lblas -llapack
	FC=mpif90
endif
ifeq ($(COMPUTER),roothan)
	EXEC=$(EXEC1)_roothan
	LIBS=-lblas -llapack
	FC=mpif90
endif
ifeq ($(COMPUTER),lc2)
	EXEC=$(EXEC1)_lc2
	LIBS=-lblas -llapack
	FC=mpif90
endif
ifeq ($(COMPUTER),mogon)
	EXEC=$(EXEC1)_mogon
	FC=mpif90
    LIBS=-lmylapack -lmyblas
    LDFLAGS=-L/home/calcavec/lapack-3.4.2/
endif
ifeq ($(COMPUTER),Majorana)
	EXEC=$(EXEC1)_majorana
	FC=mpif90
	LIBS=-lblas -llapack
endif

# elenco degli oggetti da creare (c'è sempre l'oggetto del programma principale e gli oggetti
# dei vari moduli)
OBJECTS=main.o module_dft.o module_dati.o module_hartree.o module_stati_eccitati.o module_generic_tools.o module_VMC.o module_estimatori.o module_momenta.o module_calcola_accettazione.o module_funzione_onda.o module_walkers.o module_variational_opt.o module_variational_calculations.o


# crea tutti gli oggetti e li linka correttamente con le librerie; scrivere "make"
all: $(OBJECTS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o ./$(EXEC) $(OBJECTS) $(LIBS); cp HswfQMC* ../

# rimuove tutti gli oggetti e moduli creati lasciando solo i sorgenti; "make clean"
clean:
	rm -f *\.mod *\.o *~ HswfQMC*


# crea singolarmente i vari oggetti (se per caso vuoi crearne solo alcuni a mano); "make nome"
# attenzione a inserire tutti i moduli da cui ogni modulo può dipendere
# la sintassi vuole il tab nella seconda riga
main.o: module_dft.o module_stati_eccitati.o module_walkers.o module_hartree.o module_funzione_onda.o module_calcola_accettazione.o module_estimatori.o module_generic_tools.o module_dati.o module_VMC.o module_variational_opt.o main.f90 
	$(FC) $(FFLAGS) -c main.f90 

module_variational_calculations.o: module_walkers.o module_dati.o module_funzione_onda.o module_momenta.o module_estimatori.o module_calcola_accettazione.o module_estimatori.o module_variational_calculations.f90
	$(FC) $(FFLAGS) -c module_variational_calculations.f90 

module_VMC.o: module_stati_eccitati.o module_walkers.o module_calcola_accettazione.o module_estimatori.o module_generic_tools.o module_dati.o module_momenta.o module_variational_calculations.o module_VMC.f90
	$(FC) $(FFLAGS) -c module_VMC.f90 

module_estimatori.o: module_calcola_accettazione.o module_momenta.o module_funzione_onda.o module_dati.o module_walkers.o module_generic_tools.o module_estimatori.f90
	$(FC) $(FFLAGS) -c module_estimatori.f90 

module_calcola_accettazione.o: module_walkers.o module_dati.o module_funzione_onda.o module_momenta.o module_calcola_accettazione.f90
	$(FC) $(FFLAGS) -c module_calcola_accettazione.f90

module_walkers.o: module_generic_tools.o module_dati.o module_walkers.f90
	$(FC) $(FFLAGS) -c module_walkers.f90 

module_generic_tools.o: module_generic_tools.f90
	$(FC) $(FFLAGS) -c module_generic_tools.f90 

module_dati.o: module_generic_tools.o module_dati.f90
	$(FC) $(FFLAGS) -c module_dati.f90 

module_hartree.o: module_dati.o module_walkers.o module_hartree.f90
	$(FC) $(FFLAGS) -c module_hartree.f90 

module_funzione_onda.o: module_dft.o module_momenta.o module_generic_tools.o module_dati.o module_hartree.o module_walkers.o module_funzione_onda.f90
	$(FC) $(FFLAGS) -c module_funzione_onda.f90 

module_momenta.o: module_momenta.f90
	$(FC) $(FFLAGS) -c module_momenta.f90 

module_variational_opt.o: module_walkers.o module_VMC.o module_funzione_onda.o module_dati.o module_variational_calculations.o module_generic_tools.o module_variational_opt.f90
	$(FC) $(FFLAGS) -c module_variational_opt.f90 

module_stati_eccitati.o: module_dati.o module_funzione_onda.o module_walkers.o module_stati_eccitati.f90
	$(FC) $(FFLAGS) -c module_stati_eccitati.f90 

module_dft.o: module_dft.f90
	$(FC) $(FFLAGS) -c module_dft.f90 


# se invece vuoi fare tutto a mano basta che fai.:
# gfortran (-flags) -c mod1.f90
# gfortran (-flags) -c mod2.f90
# gfortran (-flags) -c mod3.f90
# gfortran (-flags) -c main.f90
# gfortran (-flags) -o exe mod1.o mod2.o mod3.o main.o
