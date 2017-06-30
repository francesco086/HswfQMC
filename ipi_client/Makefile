FC=gfortran
FFLAGS=-O2 -ffree-line-length-none
EXEC=ipi_client.x

OBJECTS=ipi_client.o sockets.o molshift.o

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS)

.PHONY : clean
clean:
	rm -f *\.mod *\.o ipi_client.x

ipi_client.o: sockets.o molshift.o ipi_client.f90
	$(FC) $(FFLAGS) -c ipi_client.f90 

sockets.o: sockets.c
	$(FC) $(FFLAGS) -c sockets.c 

molshift.o: molshift.f90
	$(FC) $(FFLAGS) -c molshift.f90
