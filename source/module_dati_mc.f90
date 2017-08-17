MODULE dati_mc
	IMPLICIT NONE
	INCLUDE 'mpif.h'
	CHARACTER(LEN=100) :: random_seed_path, path_dati_funzione_onda
	LOGICAL, SAVE :: flag_output, flag_email_notification
	LOGICAL, PROTECTED, SAVE :: iniz_dati_mc=.FALSE., iniz_MPI=.FALSE.
	LOGICAL, PROTECTED, SAVE :: flag_E_kin, flag_E_pot, flag_E_tot, flag_gr, flag_posizioni, flag_continua
	LOGICAL, PROTECTED, SAVE :: flag_TABC, flag_elettroni, flag_protoni, flag_shadow, flag_mpi, flag_normalizza_pos
	LOGICAL, PROTECTED, SAVE :: flag_disk, flag_random_file, flag_write_parallel, flag_somme_ewald, trimer_steps
   LOGICAL, PROTECTED, SAVE :: shadow_constr_domain
	LOGICAL, PROTECTED, SAVE :: stampa_dati_funzione_onda
	CHARACTER(LEN=4), SAVE :: howtomove, propmove
	CHARACTER(LEN=7), SAVE :: what_to_do
	INTEGER, PROTECTED, SAVE :: N_1ppt, N_TABC, num_k_ewald, N_AV, mpi_ierr, mpi_myrank, mpi_nprocs
	INTEGER (KIND=8), PROTECTED, SAVE :: N_mc, N_blank, N_hist, N_mc_relax_TABC
	INTEGER, SAVE :: cont_N_mc_increase, quick_error, acceptance_rate
	REAL (KIND=8), PROTECTED, SAVE :: step_e, step_p, step_se
	REAL (KIND=8), PROTECTED, SAVE :: alpha_ewald, accuracy_energy_opt
	LOGICAL, PROTECTED, SAVE :: opt_c_eff_dnfH, opt_A_Jee, opt_F_Jee, opt_A_Jep, opt_F_Jep, opt_Jse, opt_Kse, opt_Jsesp
	LOGICAL, PROTECTED, SAVE :: opt_rp, opt_SDse, opt_SDe, opt_orbital, opt_dynamic_backflow, opt_L
	REAL (KIND=8), SAVE :: time_VMC_start
   CHARACTER (LEN=9) :: SR_kind
   LOGICAL, PROTECTED, SAVE :: SR_adaptative_beta, SR_lambda, SR_lambda_Rp
   LOGICAL, PROTECTED, SAVE :: SR_change_bound
   REAL(KIND=8), SAVE :: SR_beta, SR_beta_Rp, SR_max_change, SR_min_change, SR_max_SVD_MIN, SR_maxdeltaPsi
   REAL(KIND=8), SAVE :: lambda_init, min_lambda, max_lambda, lambda_Rp_init, min_lambda_Rp, max_lambda_Rp
   INTEGER, PROTECTED, SAVE :: SR_num_max_WO_MIN
	
	CONTAINS
	
	SUBROUTINE inizializza_dati_mc()
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		NAMELIST /dati_mc/ N_mc, N_blank, N_1ppt, flag_TABC, N_TABC, N_mc_relax_TABC, step_e, step_se, step_p, &
 		  acceptance_rate, flag_continua, howtomove, propmove, trimer_steps, shadow_constr_domain, flag_elettroni, &
		  flag_protoni, flag_shadow, flag_E_tot, flag_E_kin, flag_E_pot, flag_somme_ewald, alpha_ewald, &
		  num_k_ewald, flag_gr, N_hist, flag_posizioni, flag_normalizza_pos, N_AV, flag_mpi, what_to_do, &
		  stampa_dati_funzione_onda, path_dati_funzione_onda, accuracy_energy_opt, &
		  flag_disk, flag_output, &
		  quick_error, flag_random_file, random_seed_path
		
		NAMELIST /dati_ottimizzazione/ opt_SDe, opt_orbital, opt_dynamic_backflow, opt_A_Jee, opt_F_Jee, opt_A_Jep, &
         opt_F_Jep, opt_Jse, opt_Kse, opt_Jsesp, opt_SDse, opt_c_eff_dnfH, opt_rp, opt_L

      NAMELIST /dati_SR/ SR_kind, SR_num_max_WO_MIN, SR_beta, SR_beta_Rp, SR_maxdeltaPsi, SR_max_SVD_MIN, &
         SR_change_bound, SR_min_change, SR_max_change, SR_adaptative_beta, &
         SR_lambda, lambda_init, min_lambda, max_lambda, &
         SR_lambda_Rp, lambda_Rp_init, min_lambda_Rp, max_lambda_Rp
		
		CALL CPU_TIME(time_VMC_start)
		
		OPEN (2, FILE='dati_mc.d',STATUS='OLD')
		READ (2, NML=dati_mc)
		CLOSE (2)
		
		OPEN (2, FILE='dati_ottimizzazione.d',STATUS='OLD')
		READ (2, NML=dati_ottimizzazione)
		CLOSE (2)

		OPEN (2, FILE='dati_SR.d',STATUS='OLD')
		READ (2, NML=dati_SR)
		CLOSE (2)
		
		IF (flag_continua .AND. (.NOT. flag_disk)) STOP 'Non puoi continuare se non hai scritto su disco i dati &
		  [ module_dati.f90 > inizializza_dati_mc ]'
		IF (alpha_ewald==-1.d0) alpha_ewald=5.d0/MIN(L(1),L(2),L(3))
		IF ((iniz_MPI).AND.(N_mc>0)) THEN
			N_mc=N_mc/mpi_nprocs
		END IF
		IF (N_mc<0) N_mc=-N_mc 
		IF (N_AV<0) N_AV=MAX(-N_mc/N_AV,1)
		IF (N_1ppt<0) THEN
			N_1ppt=CEILING(-(REAL(N_1ppt)/100.)*REAL(N_part))
			IF (flag_shadow) N_1ppt=N_1ppt*3
		END IF
		IF (N_TABC<0) THEN
			IF (-N_TABC>N_mc) THEN
				N_TABC=1
			ELSE
				N_TABC=-CEILING(REAL(N_mc)/REAL(N_TABC))
			END IF
		END IF

		iniz_dati_mc=.TRUE.
		
	END SUBROUTINE inizializza_dati_mc
!-----------------------------------------------------------------------

	SUBROUTINE inizializza_MPI()
		IMPLICIT NONE
		
		CALL inizializza_dati_mc()
		
		IF (flag_mpi) THEN
			IF (.NOT. iniz_dati_mc) STOP 'Non puoi inizializzare MPI senza prima leggere i dati iniziali &
			  [ module_dati.f90 > inizializza_MPI ]'
			CALL MPI_INIT(mpi_ierr)
			IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore nell inizializzazione di MPI &
			  [ module_dati.f90 > inizializza_MPI ]'
			CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nprocs, mpi_ierr)
			IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_COMM_SIZE &
			  [ module_dati.f90 > inizializza_MPI ]'
			CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_myrank, mpi_ierr)
			IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_COMM_RANK &
			  [ module_dati.f90 > inizializza_MPI ]'
			iniz_MPI=.TRUE.
		ELSE
			mpi_nprocs=1
			mpi_myrank=0
		END IF
		
	END SUBROUTINE inizializza_MPI
!-----------------------------------------------------------------------

	SUBROUTINE change_step(rapp_step,quale)
      USE dati_fisici
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: rapp_step
		CHARACTER(LEN=3) :: quale
		SELECT CASE (quale)
		CASE ('e__')
			step_e=MIN(step_e+step_e*rapp_step,MAXVAL(L))
		CASE ('p__')
			step_p=MIN(step_p+step_p*rapp_step,MAXVAL(L))
		CASE ('se_')
			step_se=MIN(step_se+step_se*rapp_step,MAXVAL(L))
		END SELECT
	END SUBROUTINE change_step
!-----------------------------------------------------------------------

	SUBROUTINE replace_step(nuovo_step,quale)
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: nuovo_step
		CHARACTER(LEN=3) :: quale
		SELECT CASE (quale)
		CASE ('e__')
			step_e=nuovo_step
		CASE ('p__')
			step_p=nuovo_step
		CASE ('se_')
			step_se=nuovo_step
		END SELECT
	END SUBROUTINE replace_step
!-----------------------------------------------------------------------

	SUBROUTINE change_flag_continua(new_flag)
		IMPLICIT NONE
		LOGICAL, INTENT(IN) :: new_flag
		flag_continua=new_flag
	END SUBROUTINE change_flag_continua	
!-----------------------------------------------------------------------

	SUBROUTINE change_flag_TABC(new_flag)
		IMPLICIT NONE
		LOGICAL, INTENT(IN) :: new_flag
		flag_TABC=new_flag
	END SUBROUTINE change_flag_TABC
!-----------------------------------------------------------------------

	SUBROUTINE cambia_N_mc(N_mc_new)
		IMPLICIT NONE
		INTEGER (KIND=8) :: N_mc_new
		N_mc=N_mc_new 
	END SUBROUTINE cambia_N_mc
!-----------------------------------------------------------------------

	SUBROUTINE termina_MPI()
		IMPLICIT NONE
		IF (flag_mpi) THEN
			IF (.NOT. iniz_MPI) STOP 'Non puoi inizializzare MPI senza prima averlo inizializzato &
			  [ module_dati.f90 > termina_MPI ]'
			CALL MPI_FINALIZE(mpi_ierr)
			IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore nel finalizzare MPI &
			  [ module_dati.f90 > termina_MPI ]'
		END IF
	END SUBROUTINE termina_MPI
!-----------------------------------------------------------------------

	SUBROUTINE chiudi_dati_mc()
		IMPLICIT NONE
		iniz_dati_mc=.FALSE.
	END SUBROUTINE chiudi_dati_mc
	
END MODULE dati_mc


