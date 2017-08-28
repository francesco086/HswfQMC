PROGRAM main
	USE dati_mc
	USE VMC
	USE variational_opt
	IMPLICIT NONE
	LOGICAL :: success, flag_file
	CHARACTER(LEN=2) :: id
	CHARACTER(LEN=7) :: task_to_perform
	INTEGER :: i, n1_char, n2_char
	REAL (KIND=8) :: accuracy_energy
	REAL :: time1, time2, delta_t, t_sum

   OPEN(UNIT=8, FILE='state.out', STATUS='REPLACE', ACTION='WRITE')
   WRITE(8,*) 'exec'
   CLOSE(UNIT=8)

	!initialize MPI (in module_dati_mc.f90)
	CALL inizializza_MPI()
	
	IF (flag_output) THEN
		INQUIRE(FILE='output.d',EXIST=flag_file)
		IF ((flag_file).AND.(mpi_myrank==0)) THEN
			OPEN (UNIT=7, FILE='output.d', STATUS='OLD', POSITION='APPEND')
			IF (mpi_myrank==0) THEN
				PRINT * , 'Esisteva giá un file output.d, proseguo.'
				IF (flag_output) WRITE (7, *) 'Esisteva giá un file output.d, proseguo.'
			END IF
		ELSE IF (mpi_myrank==0) THEN
			OPEN (UNIT=7, FILE='output.d', STATUS='NEW')
		END IF
	END IF	
	IF (mpi_myrank==0) CALL stampa_logo_HswfQMC()
	
	CALL inizializza_dati_mc()
	IF (flag_mpi) THEN
		IF (flag_random_file) THEN
			CALL mpi_random_seed_da_file(random_seed_path)
		ELSE
			CALL mpi_random_seed(mpi_myrank)
		END IF
	END IF
	IF (mpi_myrank/=0) flag_output=.FALSE. 
	CALL chiudi_dati_mc
	
	CALL inizializza_dati_mc()
	task_to_perform=what_to_do
	accuracy_energy=accuracy_energy_opt
	CALL chiudi_dati_mc()

	SELECT CASE(task_to_perform)
	CASE('simpcal')
		CALL inizializza_VMC('eva01')
		CALL valuta_step_mc(success)
		CALL sampling(.FALSE.)
		CALL CPU_TIME(time1)
		CALL sampling(.TRUE.)
		CALL CPU_TIME(time2)
		CALL trascrivi_dati()
		CALL stampa_risultati()
		CALL chiudi_VMC()
		delta_t=time2-time1
		CALL MPI_REDUCE(delta_t,t_sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
		IF (mpi_myrank==0) THEN
			PRINT *, ""
			PRINT '(A39,F10.1,A6)' , '< < < TEMPO EFFETTIVO CAMPIONAMENTO = ', t_sum, ' > > >'
			IF (flag_output) WRITE (7, '(A39,F10.1,A6)') &
			  '< < < TEMPO EFFETTIVO CAMPIONAMENTO = ', t_sum, ' > > >'
			PRINT '(1X,A19,E12.5,A6)' , '< < < EFFICIENZA = ', 1./(t_sum*variance_efficiency), ' > > >'
			IF (flag_output) WRITE (7, '(1X,A19,E12.5,A6)') &
			  '< < < EFFICIENZA = ', 1./(t_sum*variance_efficiency), ' > > >'
		END IF
   CASE('q-e_gen')
      CALL inizializza_VMC('QEgen')
      CALL chiudi_VMC()
	CASE('congrad')
		CALL minimizza_energia(accuracy_energy,task_to_perform)
	CASE('axisopt')
		CALL minimizza_energia(accuracy_energy,task_to_perform)
	CASE('stocrec')
		CALL minimizza_energia(accuracy_energy,task_to_perform)
	CASE('stoc_ns')
		CALL minimizza_energia(accuracy_energy,task_to_perform)
	CASE('stoc_av')
		CALL minimizza_energia(accuracy_energy,task_to_perform)
	CASE('excstat')
		CALL inizializza_VMC('eva01')
		CALL inizializza_stati_eccitati()
		CALL load_stato_eccitato(1)
		CALL chiudi_stati_eccitati()
		CALL valuta_step_mc(success)
		CALL sampling(.FALSE.)
		CALL sampling(.TRUE.)
		CALL trascrivi_dati()
		CALL stampa_risultati()
		CALL chiudi_VMC()
	END SELECT
	
	IF (mpi_myrank==0 .AND. flag_output) CLOSE (7)
	
	CALL termina_MPI()

    OPEN(UNIT=8, FILE='state.out', STATUS='REPLACE', ACTION='WRITE')
    WRITE(8,*) 'done'
    CLOSE(UNIT=8)

END PROGRAM main
