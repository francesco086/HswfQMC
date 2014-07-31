MODULE dati_fisici
	USE lattice
	IMPLICIT NONE
	LOGICAL, PROTECTED, SAVE :: iniz_dati_fisici=.FALSE., flag_molecular, flag_2D=.FALSE.
	CHARACTER(LEN=5), PROTECTED, SAVE :: crystal_cell
	INTEGER, PROTECTED, SAVE :: N_part, H_N_part, N_cell_side
	REAL (KIND=8), PARAMETER, PRIVATE :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
	REAL (KIND=8), PARAMETER :: MASS_e=0.0005485899094d0, MASS_p=1.00727646677d0   !in uma
	REAL (KIND=8), PROTECTED, SAVE :: r_s, L(1:3), H_L(1:3), L_cov_bond
	REAL (KIND=8), PROTECTED, SAVE :: hbar, K_coulomb, strecthing_cov_bond
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: r_crystal(:,:)
	
	CONTAINS
	
	SUBROUTINE inizializza_dati_fisici(file_reticolo_opt)
		USE generic_tools
		IMPLICIT NONE
		CHARACTER(LEN=*), OPTIONAL :: file_reticolo_opt     !per ottimizzare le posizioni protoniche con SR
		INTEGER :: i, j, i_seed
		INTEGER, ALLOCATABLE :: seed(:), seed_provv(:)
		CHARACTER(LEN=100) :: file_reticolo
		REAL (KIND=8) :: vect(1:3), sigma_w, eta_w(1:1,1:1), L_w, dist(0:3)
		REAL (KIND=8), ALLOCATABLE :: app(:,:)
		
		NAMELIST /dati_fisici/ r_s, crystal_cell, file_reticolo, flag_molecular, &
		  strecthing_cov_bond, N_cell_side
		OPEN (2, FILE='dati_fisici.d',STATUS='OLD')
		READ (2,NML=dati_fisici)
		CLOSE (2)
		
		hbar=DSQRT(10.9748d-4)   !definisco la costante di plank riscalata per avere l'energia in Rydberg
		K_coulomb=2.0d0    !definisco la costante di Coulomb riscalata per avere l'energia in Rydberg
		
		L_cov_bond=strecthing_cov_bond*0.74d0/0.53d0    !in bohr, 0.53 é il fattore di conversione fra bohr e Angstrom
		
		IF ( crystal_cell=='bcc__' ) THEN
			N_part=2*N_cell_side**3
		ELSE IF ( crystal_cell=='fcc__' ) THEN
			N_part=4*N_cell_side**3
		ELSE IF ( (crystal_cell=='hcp__') .AND. (.NOT. flag_molecular) ) THEN
			N_part=8*N_cell_side**3
		ELSE IF ( (crystal_cell=='hcp__') .AND. (flag_molecular) ) THEN
			N_part=16*N_cell_side**3
		ELSE IF ( (crystal_cell=='hcp_w') .AND. (flag_molecular) ) THEN
			N_part=16*N_cell_side**3
		ELSE IF ( (crystal_cell=='mhcpo') .AND. (.NOT. flag_molecular) ) THEN
			STOP 'mhcpo deve essere associato ad una fase molecolare'
		ELSE IF ( (crystal_cell=='mhcpo') .AND. (flag_molecular) ) THEN
			N_part=16*N_cell_side**3
		ELSE IF ( crystal_cell=='sc___' ) THEN
			N_part=N_cell_side**3
		ELSE IF ( crystal_cell=='mol__' ) THEN
			N_part=2
		ELSE IF ( crystal_cell=='dat__' ) THEN
			N_part=N_cell_side
		ELSE IF ( crystal_cell=='datex' ) THEN
			N_part=N_cell_side
		ELSE IF ( crystal_cell=='grp__' ) THEN
			N_part=4*(N_cell_side**2)
			flag_2D=.TRUE.
		ELSE
			STOP "DATI_FISICI: scegli un reticolo accettabile"
		END IF
		H_N_part=N_part/2
		
		IF (PRESENT(file_reticolo_opt)) THEN
			crystal_cell='datex'
			file_reticolo=file_reticolo_opt
		END IF
		
		L(1:3)=r_s*(4.d0*PI*N_part/3.d0)**(1.d0/3.d0)            !in bohr units
		IF ((crystal_cell=='hcp__').OR.(crystal_cell=='mhcpo').OR.(crystal_cell=='hcp_w')) THEN
			L(1)=L(1)*(2.d0**(1.d0/6.d0))
			L(2)=L(2)*(((3.d0**3.d0)/(2.d0**5.d0))**(1.d0/6.d0))
			L(3)=L(3)*((2.d0**(2.d0/3.d0))/(3**(1.d0/2.d0)))   *(1.58d0/DSQRT(8.d0/3.d0))
		END IF
		H_L=0.5d0*L
		IF ( crystal_cell=='grp__' ) THEN
			L(1:2)=r_s*DSQRT(PI*N_part)
			L(1)=L(1)*3.d0/DSQRT(DSQRT(27.d0))
			L(2)=L(2)*DSQRT(3.d0)/DSQRT(DSQRT(27.d0))
			L(3)=L(1)
		END IF
		
		ALLOCATE(r_crystal(1:3,1:N_part),app(1:3,1:N_part))
		IF ( crystal_cell=='bcc__' ) THEN
			!posiziono i protoni in modo che la prima metà e la seconda metà siano ben equidistribuiti nel box di simulazione (conviene per via dello spin)
			app=bodycenteredcubic(N_cell_side)
			j=1
			DO i = 1, N_part, 2
				r_crystal(1:3,j)=app(1:3,i)
				r_crystal(1:3,j+H_N_part)=app(1:3,i+1)
				j=j+1
			END DO
			DO j = 1, N_part, 1
				r_crystal(1:3,j)=(r_crystal(1:3,j)-0.5d0)*L(1)
			END DO
		ELSE IF ( crystal_cell=='fcc__' ) THEN
			app=facecenteredcubic(N_cell_side)
			j=1
			DO i = 1, N_part, 2
				r_crystal(1:3,j)=app(1:3,i)
				r_crystal(1:3,j+H_N_part)=app(1:3,i+1)
				j=j+1
			END DO
			DO j = 1, N_part, 1
				r_crystal(1:3,j)=(r_crystal(1:3,j)-0.5d0)*L(1)
			END DO
		ELSE IF ( (crystal_cell=='hcp__') .AND. (.NOT. flag_molecular) ) THEN
			app=hexagonalclosepacking(N_cell_side)
			DO i = 1, N_part, 1
				app(1,i)=app(1,i)*0.5d0
				app(2,i)=app(2,i)/DSQRT(3.d0)
				app(3,i)=app(3,i)/(2.d0*DSQRT(2.d0/3.d0))
			END DO
			j=1
			DO i = 1, N_part, 2
				r_crystal(1:3,j)=app(1:3,i)
				r_crystal(1:3,j+H_N_part)=app(1:3,i+1)
				j=j+1
			END DO
			DO j = 1, N_part, 1
				DO i = 1, 3, 1
					r_crystal(i,j)=(r_crystal(i,j)-0.5d0)*L(i)
				END DO
			END DO
		ELSE IF ( (crystal_cell=='hcp__') .AND. (flag_molecular) ) THEN
			app(1:3,1:H_N_part)=hexagonalclosepacking(N_cell_side)
			DO i = 1, N_part, 1
				app(1,i)=app(1,i)*0.5d0
				app(2,i)=app(2,i)/DSQRT(3.d0)
				app(3,i)=app(3,i)/(2.d0*DSQRT(2.d0/3.d0))
			END DO
			DO j = 1, N_part, 1
				DO i = 1, 3, 1
					app(i,j)=(app(i,j)-0.5d0)*L(i)
				END DO
			END DO
			DO i = 1, H_N_part, 1
				IF (MOD(i-1,8)<4) THEN
					vect=0.5d0*L_cov_bond*(/DCOS(PI/6.d0)*DCOS(PI/3.d0),DSIN(PI/6.d0)*DCOS(PI/3.d0),-DSIN(PI/3.d0)/)   !vecchi calcoli
				ELSE
					vect=0.5d0*L_cov_bond*(/0.d0,DCOS(PI/3.d0),DSIN(PI/3.d0)/)
				END IF
				r_crystal(1:3,i)=app(1:3,i)-vect(1:3)
				r_crystal(1:3,i+H_N_part)=app(1:3,i)+vect(1:3)
			END DO
		ELSE IF ( (crystal_cell=='hcp_w') .AND. (flag_molecular) ) THEN
			
			CALL RANDOM_SEED(size = i_seed)
			ALLOCATE(seed(i_seed),seed_provv(i_seed))
			CALL RANDOM_SEED(GET=seed)
			seed_provv = 7 + 199 * (/ (i - 1, i = 1, i_seed) /)
			CALL RANDOM_SEED(PUT = seed_provv)
			
			app(1:3,1:H_N_part)=hexagonalclosepacking(N_cell_side)
			DO i = 1, N_part, 1
				app(1,i)=app(1,i)*0.5d0
				app(2,i)=app(2,i)/DSQRT(3.d0)
				app(3,i)=app(3,i)/(2.d0*DSQRT(2.d0/3.d0))
			END DO
			DO j = 1, N_part, 1
				DO i = 1, 3, 1
					app(i,j)=(app(i,j)-0.5d0)*L(i)
				END DO
			END DO
			DO i = 1, H_N_part, 1
				sigma_w=DSQRT(hbar/(MASS_p*2.d0*0.166569))    !2*0.166569 é la w per il legame H2
				CALL gaussian_sample(eta_w,1,1,sigma_w)
				eta_w(1,1)=eta_w(1,1)+1.d0
				L_w=eta_w(1,1)
				IF (MOD(i-1,8)<4) THEN
					vect=0.5d0*L_w*L_cov_bond*(/DCOS(PI/6.d0)*DCOS(PI/3.d0),DSIN(PI/6.d0)*DCOS(PI/3.d0),-DSIN(PI/3.d0)/)   !vecchi calcoli
				ELSE
					vect=0.5d0*L_w*L_cov_bond*(/0.d0,DCOS(PI/3.d0),DSIN(PI/3.d0)/)
				END IF
				r_crystal(1:3,i)=app(1:3,i)-vect(1:3)
				r_crystal(1:3,i+H_N_part)=app(1:3,i)+vect(1:3)
			END DO
			
			CALL RANDOM_SEED(PUT = seed)
		ELSE IF ( (crystal_cell=='mhcpo') .AND. (flag_molecular) ) THEN
			app(1:3,1:H_N_part)=hexagonalclosepacking(N_cell_side)
			DO i = 1, N_part, 1
				app(1,i)=app(1,i)*0.5d0
				app(2,i)=app(2,i)/DSQRT(3.d0)
				app(3,i)=app(3,i)/(2.d0*DSQRT(2.d0/3.d0))
			END DO
			DO j = 1, N_part, 1
				DO i = 1, 3, 1
					app(i,j)=(app(i,j)-0.5d0)*L(i)
				END DO
			END DO
			DO i = 1, H_N_part, 1
				IF (MOD(i-1,8)<4) THEN
					vect=0.5d0*L_cov_bond*(/DCOS(-PI/6.d0)*DCOS(-PI/6.d0),DSIN(-PI/6.d0),DSIN(-PI/6.d0)/)
				ELSE
					vect=0.5d0*L_cov_bond*(/DCOS(PI/6.d0)*DCOS(-PI/6.d0),DSIN(-PI/6.d0),DSIN(PI/6.d0)/)
				END IF
				r_crystal(1:3,i)=app(1:3,i)-vect(1:3)
				r_crystal(1:3,i+H_N_part)=app(1:3,i)+vect(1:3)
			END DO
		ELSE IF ( crystal_cell=='sc___' ) THEN
			app=simplecubic(N_cell_side)
			j=1
			DO i = 1, N_part, 2
				r_crystal(1:3,j)=app(1:3,i)
				r_crystal(1:3,j+H_N_part)=app(1:3,i+1)
				j=j+1
			END DO
			DO j = 1, N_part, 1
				DO i = 1, 3, 1
					r_crystal(i,j)=(r_crystal(i,j)-0.5d0)*L(i)
				END DO
			END DO
		ELSE IF ( crystal_cell=='mol__' ) THEN
			r_crystal(1:3,1)=(/L_cov_bond*0.5d0,0.d0,0.d0/)
			r_crystal(1:3,2)=(/-L_cov_bond*0.5d0,0.d0,0.d0/)
		ELSE IF ( crystal_cell=='dat__' ) THEN
			OPEN (UNIT=2, FILE=TRIM(file_reticolo), STATUS='OLD')
			DO i = 1, N_part, 1
				READ (2, *), r_crystal(1:3,i)
				r_crystal(1:3,i)=r_crystal(1:3,i)+(/ 0.5d0 , 0.5d0 , 0.5d0 /)
			END DO
			CLOSE (2)
			DO j = 1, N_part, 1
				DO i = 1, 3, 1
					r_crystal(i,j)=(r_crystal(i,j)-0.5d0)*L(i)
				END DO
			END DO
		ELSE IF ( crystal_cell=='datex' ) THEN
			OPEN (UNIT=2, FILE=TRIM(file_reticolo), STATUS='OLD')
			READ (2, *), L(1:3)
			H_L=0.5d0*L
			DO i = 1, N_part, 1
				READ (2, *), r_crystal(1:3,i)
			END DO
			CLOSE (2)
		ELSE IF ( crystal_cell=='grp__' ) THEN
			r_crystal=graphene_layer(N_cell_side)
			DO i = 1, N_part, 1
				DO j = 1, 2, 1
					r_crystal(j,i)=r_crystal(j,i)*L(j)
				END DO
			END DO
			IF ( flag_molecular ) THEN
				DO i = 1, H_N_part, 1
					dist(1:3)=r_crystal(1:3,i)-r_crystal(1:3,i+H_N_part)
					dist(0)=DSQRT(DOT_PRODUCT(dist(1:3),dist(1:3)))
					dist(1:3)=dist(1:3)/dist(0)
					IF ( dist(0)>L_cov_bond ) THEN
						r_crystal(1:3,i)=r_crystal(1:3,i)-dist(1:3)*(dist(0)-L_cov_bond)*0.5d0
						r_crystal(1:3,i+H_N_part)=r_crystal(1:3,i+H_N_part)+dist(1:3)*(dist(0)-L_cov_bond)*0.5d0
					ELSE
						STOP 'distanze giá minori che per la fase molecolare &
						  [ module_funzione_onda.f90 > inizializza_funzione_onda ]'
					END IF
				END DO
			END IF
		ELSE
			STOP "scegli un reticolo accettabile"
		END IF
		
		CALL applica_pbc(r_crystal,N_part,L)
		
		IF (MOD(N_part,2)/=0) STOP 'Stai lavorando con un numero di particelle non pari!!! &
		  [ module_funzione_onda.f90 > inizializza_funzione_onda ]'
		
		DEALLOCATE(app)
		
		iniz_dati_fisici=.TRUE.
	END SUBROUTINE inizializza_dati_fisici
!-----------------------------------------------------------------------
	SUBROUTINE chiudi_dati_fisici()
		IMPLICIT NONE
		IF (.NOT. iniz_dati_fisici) STOP 'Prima di chiudere avresti dovuto inizializzare &
		  [ module_dati.f90 > chiudi_dati_fisici ]'
		DEALLOCATE(r_crystal)
		iniz_dati_fisici=.FALSE.
	END SUBROUTINE chiudi_dati_fisici
	
END MODULE dati_fisici


MODULE dati_simulazione_mc
	IMPLICIT NONE
	INCLUDE 'mpif.h'
	CHARACTER(LEN=100) :: random_seed_path, path_dati_funzione_onda
	LOGICAL, SAVE :: flag_output
	LOGICAL, PROTECTED, SAVE :: iniz_dati_simulazione_mc=.FALSE., iniz_MPI=.FALSE.
	LOGICAL, PROTECTED, SAVE :: flag_E_kin, flag_E_pot, flag_E_tot, flag_gr, flag_posizioni, flag_continua
	LOGICAL, PROTECTED, SAVE :: flag_TABC, flag_elettroni, flag_protoni, flag_shadow, flag_mpi, flag_normalizza_pos
	LOGICAL, PROTECTED, SAVE :: flag_disk, flag_random_file, flag_write_parallel, flag_somme_ewald, trimer_steps
	LOGICAL, PROTECTED, SAVE :: stampa_dati_funzione_onda, flag_email_notification
	CHARACTER(LEN=4), SAVE :: howtomove, propmove
	CHARACTER(LEN=7), SAVE :: what_to_do
	INTEGER, PROTECTED, SAVE :: N_1ppt, N_TABC, num_k_ewald, N_AV, mpi_ierr, mpi_myrank, mpi_nprocs
	INTEGER (KIND=8), PROTECTED, SAVE :: N_mc, N_blank, N_hist, N_mc_relax_TABC
	INTEGER, SAVE :: cont_N_mc_increase, quick_error, acceptance_rate
	REAL (KIND=8), PROTECTED, SAVE :: step_e, step_p, step_se
	REAL (KIND=8), PROTECTED, SAVE :: alpha_ewald, accuracy_energy_opt
	LOGICAL, PROTECTED, SAVE :: opt_c_eff_hartree, opt_A_Jee, opt_F_Jee, opt_A_Jep, opt_F_Jep, opt_Jse, opt_Kse, opt_Jsesp
	LOGICAL, PROTECTED, SAVE :: opt_rp, opt_SDse, opt_SDe
	REAL (KIND=8), SAVE :: time_VMC_start
	
	CONTAINS
	
	SUBROUTINE inizializza_dati_simulazione_mc()
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		NAMELIST /dati_simulazione_mc/ N_mc, N_blank, N_1ppt, flag_TABC, N_TABC, N_mc_relax_TABC, step_e, step_se, step_p, &
		  acceptance_rate, flag_continua, howtomove, propmove, trimer_steps, flag_elettroni, &
		  flag_protoni, flag_shadow, flag_E_tot, flag_E_kin, flag_E_pot, flag_somme_ewald, alpha_ewald, &
		  num_k_ewald, flag_gr, N_hist, flag_posizioni, flag_normalizza_pos, N_AV, flag_mpi, what_to_do, &
		  stampa_dati_funzione_onda, path_dati_funzione_onda, accuracy_energy_opt, &
		  flag_disk, flag_output, &
		  quick_error, flag_random_file, random_seed_path, flag_email_notification
		NAMELIST /dati_ottimizzazione/ opt_SDe, opt_A_Jee, opt_F_Jee, opt_A_Jep, opt_F_Jep, opt_Jse, opt_Kse, opt_Jsesp, &
		  opt_SDse, opt_c_eff_hartree, opt_rp
		
		CALL CPU_TIME(time_VMC_start)
		
		OPEN (2, FILE='dati_mc.d',STATUS='OLD')
		READ (2, NML=dati_simulazione_mc)
		CLOSE (2)
		
		OPEN (2, FILE='dati_ottimizzazione.d',STATUS='OLD')
		READ (2, NML=dati_ottimizzazione)
		CLOSE (2)
		
		IF (flag_continua .AND. (.NOT. flag_disk)) STOP 'Non puoi continuare se non hai scritto su disco i dati &
		  [ module_dati.f90 > inizializza_dati_simulazione_mc ]'
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
		
		iniz_dati_simulazione_mc=.TRUE.
		
	END SUBROUTINE inizializza_dati_simulazione_mc
!-----------------------------------------------------------------------

	SUBROUTINE inizializza_MPI()
		IMPLICIT NONE
		
		CALL inizializza_dati_simulazione_mc()
		
		IF (flag_mpi) THEN
			IF (.NOT. iniz_dati_simulazione_mc) STOP 'Non puoi inizializzare MPI senza prima leggere i dati iniziali &
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
		IF ((mpi_nprocs>1400).AND.(flag_random_file)) THEN
			STOP 'Il file random_seed.d contiene dati per al massimo 1400 processori &
					  [ module_dati.f90 > inizializza_MPI ]'
		END IF
	END SUBROUTINE inizializza_MPI
!-----------------------------------------------------------------------

	SUBROUTINE change_step(rapp_step,quale)
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: rapp_step
		CHARACTER(LEN=3) :: quale
		SELECT CASE (quale)
		CASE ('e__')
			step_e=step_e+step_e*rapp_step
		CASE ('p__')
			step_p=step_p+step_p*rapp_step
		CASE ('se_')
			step_se=step_se+step_se*rapp_step
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

	SUBROUTINE chiudi_dati_simulazione_mc()
		IMPLICIT NONE
		iniz_dati_simulazione_mc=.FALSE.
	END SUBROUTINE chiudi_dati_simulazione_mc
	
END MODULE dati_simulazione_mc


