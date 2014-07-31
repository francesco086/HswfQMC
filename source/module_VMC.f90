MODULE VMC
	USE stati_eccitati
	USE calcola_accettazione
	USE dati_fisici
	USE dati_simulazione_mc
	USE estimatori
	USE generic_tools
	USE walkers
	IMPLICIT NONE
	LOGICAL, PRIVATE :: verbose_mode=.TRUE.
	LOGICAL, PRIVATE, SAVE :: iniz_VMC=.FALSE.
	CHARACTER(LEN=6), PROTECTED, SAVE :: codice_simulazione
	REAL (KIND=8), ALLOCATABLE, PRIVATE, SAVE :: geeuu(:), geeuu_new(:), geeud(:), geeud_new(:), gep(:), gep_new(:)
	REAL (KIND=8), ALLOCATABLE, PRIVATE, SAVE :: gppuu(:), gppuu_new(:), gppud(:), gppud_new(:)
	REAL (KIND=8), ALLOCATABLE, PRIVATE, SAVE :: gseuu1(:), gseuu1_new(:), gseud1(:), gseud1_new(:), gese1(:), gese1_new(:)
	REAL (KIND=8), ALLOCATABLE, PRIVATE, SAVE :: gseuu2(:), gseuu2_new(:), gseud2(:), gseud2_new(:), gese2(:), gese2_new(:)
	REAL (KIND=8), ALLOCATABLE, PRIVATE, SAVE :: gsesp1(:), gsesp1_new(:), gsesp2(:), gsesp2_new(:)
	REAL (KIND=8), SAVE :: variance_efficiency
	
	CONTAINS
	
	SUBROUTINE inizializza_VMC(codice, CONTINUA)
		IMPLICIT NONE
		INTEGER :: i, j, cont_i1, cont_i2
		LOGICAL :: flag_file_dat, flag_file_pos, flag_scrivi
		LOGICAL, OPTIONAL :: CONTINUA
		CHARACTER(LEN=*) :: codice
		CHARACTER(LEN=4) :: istring
		CHARACTER(LEN=8) :: istring1, istring2
		INTEGER :: status(MPI_STATUS_SIZE)
		INTEGER :: dum_myrank, dum_nprocs
		INTEGER, ALLOCATABLE :: seed(:)
		REAL (KIND=8) :: step_e_, step_se_, step_p_
		
		!PRINT * , '000: ', mpi_myrank
		
		codice_simulazione='-'//codice
		IF (mpi_myrank==0) THEN
			PRINT * , ' vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ---> ', codice
			IF (flag_output) THEN
				WRITE (7, *), ' vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ---> ', codice
				CLOSE (7)
				OPEN (UNIT=7, FILE='output.d', STATUS='OLD', POSITION='APPEND')
			END IF
		END IF
		CALL inizializza_dati_fisici()
		CALL inizializza_dati_simulazione_mc()
		IF (PRESENT(CONTINUA)) THEN
			IF (CONTINUA .AND. (.NOT. flag_disk)) CALL cambia_N_mc(N_mc*(2**(cont_N_mc_increase)))
			CALL change_flag_continua(CONTINUA)
		END IF
		IF (flag_mpi) THEN
			IF (mpi_myrank/=0) verbose_mode=.FALSE.
		END IF
		IF (verbose_mode) THEN
			PRINT * , 'VMC: Ho inizializzato i dati fisici e i dati della simulazione MC.'
			PRINT * , 'N_part=', N_part, '       r_s=', r_s
			IF (flag_output) THEN
				WRITE (7, *) , 'VMC: Ho inizializzato i dati fisici e i dati della simulazione MC.'
				WRITE (7, *) , 'N_part=', N_part, '       r_s=', r_s
			END IF
			PRINT * , 'L=', L, '   [bohr]'
			IF (flag_output) WRITE (7, *), 'L=', L, '   [bohr]'
			PRINT * , 'Tipo di reticolo: ', crystal_cell, '       Molecolare? ', flag_molecular
			IF (flag_output) WRITE (7, *), 'Tipo di reticolo: ', crystal_cell, '       Molecolare? ', flag_molecular
		END IF
		
		!PRINT * , '111: ', mpi_myrank
		
		!verifico che si possa continuare dai dati precedenti
		IF (flag_continua) THEN
			WRITE (istring, '(I4.4)'), mpi_myrank
			IF (flag_elettroni) THEN
				INQUIRE(FILE='posizioni/fine'//codice_simulazione//'_e_'//istring//'.pos',EXIST=flag_file_pos)
				IF (.NOT. flag_file_pos) CALL change_flag_continua(.FALSE.)
			END IF
			IF (flag_protoni) THEN
				INQUIRE(FILE='posizioni/fine'//codice_simulazione//'_p_'//istring//'.pos',EXIST=flag_file_pos)
				IF (.NOT. flag_file_pos) CALL change_flag_continua(.FALSE.)
			END IF
			IF (flag_shadow) THEN
				INQUIRE(FILE='posizioni/fine'//codice_simulazione//'_se1_'//istring//'.pos',EXIST=flag_file_pos)
				IF (.NOT. flag_file_pos) CALL change_flag_continua(.FALSE.)
				INQUIRE(FILE='posizioni/fine'//codice_simulazione//'_se2_'//istring//'.pos',EXIST=flag_file_pos)
				IF (.NOT. flag_file_pos) CALL change_flag_continua(.FALSE.)
			END IF
			INQUIRE(FILE='posizioni/seed_and_step'//codice_simulazione//'.dat',EXIST=flag_file_dat)
			IF (.NOT. flag_file_dat) CALL change_flag_continua(.FALSE.)
			
			IF (.NOT. flag_continua) THEN
				IF (verbose_mode) THEN
					PRINT * , 'VMC: Non é possibilie continuare dai dati precedenti'
					IF (flag_output) WRITE (7, *), 'VMC: Non é possibilie continuare dai dati precedenti'
				END IF
			END IF
		END IF
		IF (flag_mpi) CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
		
		!carico i seed e gli step dal file seed_and_step
		IF (flag_continua) THEN
			CALL RANDOM_SEED(size = j)
			ALLOCATE(seed(j))
			OPEN (UNIT=77, FILE='posizioni/seed_and_step'//codice_simulazione//'.dat', STATUS='UNKNOWN')
			READ (77, *), dum_nprocs
			IF (mpi_myrank==0) THEN
				flag_scrivi=.TRUE.
			ELSE
				flag_scrivi=.FALSE.
			END IF
			DO j = 0, MIN(dum_nprocs, mpi_nprocs)-1, 1
				IF (mpi_myrank==j .AND. flag_scrivi) THEN
					DO i = 0, mpi_myrank-1, 1
						READ (77, *)
					END DO
					READ (77, *), dum_myrank, seed, step_e_, step_se_, step_p_
					IF (dum_myrank/=mpi_myrank) STOP 'Errore nel leggere da seed_and_step.dat &
					  [ module_VMC.f90 > inizializza_VMC ]'
					IF (mpi_myrank<mpi_nprocs-1) THEN
						CALL MPI_SEND(flag_scrivi,1,MPI_LOGICAL,mpi_myrank+1,10+j,MPI_COMM_WORLD,mpi_ierr)
						IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_SEND &
						  [ module_VMC.f90 > inizializza_VMC ]'
					END IF
				END IF
				IF (mpi_myrank==j+1 .AND. j<mpi_nprocs) THEN
					CALL MPI_RECV(flag_scrivi,1,MPI_LOGICAL,mpi_myrank-1,10+j,MPI_COMM_WORLD,status,mpi_ierr)
					IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_RECV &
					  [ module_VMC.f90 > inizializza_VMC ]'
				END IF
			END DO
			DO j = dum_nprocs, mpi_nprocs-1, 1
				IF (mpi_myrank==j-1) THEN
				 	CALL MPI_SEND(step_e_,1,MPI_DOUBLE_PRECISION,mpi_myrank+1,10+j,MPI_COMM_WORLD,mpi_ierr)
					CALL MPI_SEND(step_se_,1,MPI_DOUBLE_PRECISION,mpi_myrank+1,10+j,MPI_COMM_WORLD,mpi_ierr)
					CALL MPI_SEND(step_p_,1,MPI_DOUBLE_PRECISION,mpi_myrank+1,10+j,MPI_COMM_WORLD,mpi_ierr)
				END IF
				IF (mpi_myrank==j) THEN
					CALL mpi_random_seed(mpi_myrank)
					CALL MPI_RECV(step_e_,1,MPI_DOUBLE_PRECISION,mpi_myrank-1,10+j,MPI_COMM_WORLD,status,mpi_ierr)
					CALL MPI_RECV(step_se_,1,MPI_DOUBLE_PRECISION,mpi_myrank-1,10+j,MPI_COMM_WORLD,status,mpi_ierr)
					CALL MPI_RECV(step_p_,1,MPI_DOUBLE_PRECISION,mpi_myrank-1,10+j,MPI_COMM_WORLD,status,mpi_ierr)
				END IF
			END DO
			CLOSE (77)
			CALL change_step((step_e_-step_e)/step_e,'e__')
			CALL change_step((step_se_-step_se)/step_se,'se_')
			CALL change_step((step_p_-step_p)/step_p,'p__')
			IF (verbose_mode .AND. flag_elettroni) THEN
				PRINT * , 'VMC: da file, step_e=', step_e
				IF (flag_output) WRITE (7, *), 'VMC: da file, step_e=', step_e
			END IF
			IF (verbose_mode .AND. flag_shadow) THEN
				PRINT * , 'VMC: da file, step_se=', step_se
				IF (flag_output) WRITE (7, *), 'VMC: da file, step_se=', step_se
			END IF
			IF (verbose_mode .AND. flag_protoni) THEN
				PRINT * , 'VMC: da file, step_p=', step_p
				IF (flag_output) WRITE (7, *), 'VMC: da file, step_p=', step_p
			END IF
			CALL RANDOM_SEED(PUT=seed)
			DEALLOCATE(seed)
		END IF
		
		CALL inizializza_walkers(codice_simulazione)
		IF (verbose_mode) THEN
			PRINT * , 'VMC: Ho inizializzato i walkers'
			IF (flag_output) WRITE (7, *), 'VMC: Ho inizializzato i walkers'
		END IF
		IF (flag_posizioni .AND. flag_disk) CALL salva_posizione_walkers('posizioni/inizio'//codice_simulazione)
		
		CALL inizializza_funzione_onda()
		IF (verbose_mode) THEN
			PRINT * , 'VMC: Ho inizializzato la funzione d onda'
			IF (flag_output) WRITE (7, *), 'VMC: Ho inizializzato la funzione d onda'
		END IF
		IF (verbose_mode .AND. flag_TABC) THEN
			cont_i1=0
			DO WHILE ( N_TABC/(10**cont_i1)>0 )
				cont_i1=cont_i1+1
			END DO
			cont_i2=0
			DO WHILE ( N_mc_relax_TABC/(10**cont_i2)>0 )
				cont_i2=cont_i2+1
			END DO
			WRITE (istring1, '(I8.8)'), cont_i1
			WRITE (istring2, '(I8.8)'), cont_i2
			PRINT '(A31,1X,I'//istring1//'.1,5X,A22,1X,I'//istring2//'.1)', &
			  'VMC: Utilizzo TABC con N_TABC=', N_TABC, 'rilassando di N_relax=', N_mc_relax_TABC
			IF (flag_output) WRITE (7, '(A31,1X,I'//istring1//'.1,5X,A22,1X,I'//istring2//'.1)'), &
			  'VMC: Utilizzo TABC con N_TABC=', N_TABC, 'rilassando di N_relax=', N_mc_relax_TABC
		END IF
		IF ( verbose_mode .AND. flag_traccia_coppie_mol_ss ) THEN
			cont_i1=0
			DO WHILE ( N_ritraccia_coppie/(10**cont_i1)>0 )
				cont_i1=cont_i1+1
			END DO
			cont_i2=0
			DO WHILE ( N_mc_relax_traccia_coppie/(10**cont_i2)>0 )
				cont_i2=cont_i2+1
			END DO
			WRITE (istring1, '(I8.8)'), cont_i1
			WRITE (istring2, '(I8.8)'), cont_i2
			IF ( N_mc_relax_traccia_coppie>0 ) THEN
				PRINT '(1X,A52,1X,I'//istring1//'.1,5X,A22,1X,I'//istring2//'.1)', &
				  'VMC: Ritraccio le coppie molecolari ss ogni mc_step=', N_ritraccia_coppie, &
				  'rilassando di N_relax=', N_mc_relax_traccia_coppie
				IF (flag_output) WRITE (7, '(1X,A52,1X,I'//istring1//'.1,5X,A22,1X,I'//istring2//'.1)'), &
				  'VMC: Ritraccio le coppie molecolari ss ogni mc_step=', N_ritraccia_coppie, &
				  'rilassando di N_relax=', N_mc_relax_traccia_coppie
			ELSE
				PRINT * , 'VMC: Ritraccio le coppie molecolari ss ogni mc_step=', N_ritraccia_coppie, &
				  ' senza rilassare'
				IF (flag_output) WRITE (7, *), 'VMC: Ritraccio le coppie molecolari ss ogni mc_step=', N_ritraccia_coppie, &
				  ' senza rilassare'
			END IF
		END IF
		
		CALL inizializza_estimatori()
		
		IF (flag_gr) THEN
			ALLOCATE(geeuu(1:N_hist),geeuu_new(1:N_hist),geeud(1:N_hist),geeud_new(1:N_hist))
			ALLOCATE(gep(1:N_hist),gep_new(1:N_hist))
			ALLOCATE(gppuu(1:N_hist),gppuu_new(1:N_hist),gppud(1:N_hist),gppud_new(1:N_hist))
			geeuu=0.d0
			geeuu_new=0.d0
			geeud=0.d0
			geeud_new=0.d0
			gep=0.d0
			gep_new=0.d0
			gppuu=0.d0
			gppuu_new=0.d0
			gppud=0.d0
			gppud_new=0.d0
			IF (flag_shadow) THEN
				ALLOCATE(gseuu1(1:N_hist), gseuu1_new(1:N_hist), gseud1(1:N_hist), gseud1_new(1:N_hist), gese1(1:N_hist), gese1_new(1:N_hist))
				ALLOCATE(gseuu2(1:N_hist), gseuu2_new(1:N_hist), gseud2(1:N_hist), gseud2_new(1:N_hist), gese2(1:N_hist), gese2_new(1:N_hist))
				ALLOCATE(gsesp1(1:N_hist),gsesp1_new(1:N_hist),gsesp2(1:N_hist),gsesp2_new(1:N_hist))
				gseuu1=0.d0
				gseuu1_new=0.d0
				gseud1=0.d0
				gseud1_new=0.d0
				gese1=0.d0
				gese1_new=0.d0
				gsesp1=0.d0
				gsesp1_new=0.d0
				gseuu2=0.d0
				gseuu2_new=0.d0
				gseud2=0.d0
				gseud2_new=0.d0
				gese2=0.d0
				gese2_new=0.d0
				gsesp2=0.d0
				gsesp2_new=0.d0
			END IF
		END IF
		IF (verbose_mode) THEN
			PRINT * , 'VMC: Ho allocato gli array degli estimatori'
			IF (flag_output) WRITE (7, *), 'VMC: Ho allocato gli array degli estimatori'
		END IF
		
		iniz_VMC=.TRUE.
	END SUBROUTINE inizializza_VMC
!-----------------------------------------------------------------------
	SUBROUTINE cambia_verbosity_VMC(flag)
		IMPLICIT NONE
		LOGICAL, INTENT(IN) :: flag
		verbose_mode=flag
	END SUBROUTINE cambia_verbosity_VMC
!-----------------------------------------------------------------------
	SUBROUTINE valuta_step_mc(success)
		USE funzione_onda
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: STEP_MIN=-0.0000001d0
		INTEGER, PARAMETER :: num_ok=10, num_min=1
		LOGICAL :: accettazione, eseguito_correttamente_e, eseguito_correttamente_se, eseguito_correttamente_ese
		LOGICAL:: mpi_success
		CHARACTER (LEN=3) :: tipo_mossa
		INTEGER :: num_1ppt, ieta
		INTEGER (KIND=8) :: i, acc_e, acc_se, acc_ese, cont_e, cont_se, cont_ese
		INTEGER (KIND=8) :: contacc_e, contacc_se, contacc_ese
		REAL (KIND=8) :: eta, step_partenza_p, step_partenza_e, step_partenza_se
		REAL (KIND=8) :: step_e_new, step_se_new
		REAL (KIND=8) :: AV_step_e, AV_step_se
		LOGICAL, INTENT(OUT) :: success
		
		IF (.NOT. iniz_VMC) STOP 'Prima di valutare lo step ideale devi inizializzare il VMC &
		  [ module_VMC.f90 > valuta_step_mc ]'
		
		!CALL sampling(.FALSE.,500_8)
		
		IF (mpi_myrank==0) THEN
			CALL stampa_parametri_variazionali()
		END IF
		
		cont: IF (.NOT. flag_continua) THEN
		
		eseguito_correttamente_e=.NOT. flag_elettroni
		eseguito_correttamente_se= .NOT. flag_shadow
		eseguito_correttamente_ese= .NOT. flag_shadow
						
		SELECT CASE (howtomove)
		CASE ('allp')
			IF (flag_elettroni) THEN
				step_partenza_e=step_e
				step_partenza_se=step_se
				acc_e=0
				acc_se=0
				cont_e=1
				cont_se=1
				contacc_e=0
				contacc_se=0
				loop_e: DO i = 1, MAX(N_blank,num_min), 1
					IF (flag_shadow) THEN
						CALL RANDOM_NUMBER(eta)
					ELSE
						eta=0.1d0
					END IF
					IF (eta<0.5d0) THEN
						tipo_mossa='e__'
						cont_e=cont_e+1
					ELSE
						tipo_mossa='se_'
						cont_se=cont_se+1
					END IF
					CALL proponi_mossa(tipo_mossa,-1)
					CALL valuta_accettazione(i,tipo_mossa,-1,accettazione)
					IF (accettazione) THEN
						CALL metti_rnew_in_rold(tipo_mossa,-1)
						IF (eta<0.5d0) THEN
							acc_e=acc_e+1
						ELSE
							acc_se=acc_se+1
						END IF
					END IF
					IF (MOD(cont_e,101)==0) THEN
						IF (acc_e>acceptance_rate+10) THEN
							contacc_e=0
							CALL change_step((REAL(acc_e-(acceptance_rate+10),8)/REAL(100-acceptance_rate-10,8)),'e__')
						ELSE IF (acc_e<acceptance_rate-10) THEN
							contacc_e=0
							CALL change_step((REAL(acc_e-(acceptance_rate-10),8)/REAL(acceptance_rate-10,8)),'e__')
						ELSE
							contacc_e=contacc_e+1
							IF (contacc_e==num_ok) THEN
								eseguito_correttamente_e=.TRUE.
								IF (eseguito_correttamente_se) THEN
									EXIT loop_e
								END IF
							END IF
						END IF
						acc_e=0
						cont_e=1
					END IF
					IF (MOD(cont_se,101)==0) THEN
						IF (acc_se>acceptance_rate+10) THEN
							contacc_se=0
							CALL change_step((REAL(acc_se-(acceptance_rate+10),8)/REAL(100-acceptance_rate-10,8)),'se_')
						ELSE IF (acc_se<acceptance_rate-10) THEN
							contacc_se=0
							CALL change_step((REAL(acc_se-(acceptance_rate-10),8)/REAL(acceptance_rate-10,8)),'se_')
						ELSE
							contacc_se=contacc_se+1
							IF (contacc_se==num_ok) THEN
								eseguito_correttamente_se=.TRUE.
								IF (eseguito_correttamente_e) THEN
									EXIT loop_e
								END IF
							END IF
						END IF
						acc_se=0
						cont_se=1
					END IF			
				END DO loop_e
				
				IF (flag_shadow) THEN
					acc_ese=0
					cont_ese=1
					contacc_ese=0
					loop_ese: DO i = 1, MAX(N_blank,num_min), 1
						CALL proponi_mossa('e__',-1)
						CALL proponi_mossa('se_',-1)
						cont_ese=cont_ese+1
						CALL valuta_accettazione(i,'all',-1,accettazione)
						IF (accettazione) THEN
							CALL metti_rnew_in_rold('all',-1)
							acc_ese=acc_ese+1
						END IF
						IF (MOD(cont_ese,101)==0) THEN
							IF (acc_ese>acceptance_rate+10) THEN
								contacc_ese=0
								CALL change_step((REAL(acc_ese-(acceptance_rate+10),8)/REAL(100-acceptance_rate-10,8)),'se_')
								CALL change_step((REAL(acc_ese-(acceptance_rate+10),8)/REAL(100-acceptance_rate-10,8)),'e__')
							ELSE IF (acc_ese<acceptance_rate-10) THEN
								contacc_se=0
								CALL change_step((REAL(acc_ese-(acceptance_rate-10),8)/REAL(acceptance_rate-10,8)),'se_')
								CALL change_step((REAL(acc_ese-(acceptance_rate-10),8)/REAL(acceptance_rate-10,8)),'e__')
							ELSE
								contacc_ese=contacc_ese+1
								IF (contacc_ese==num_ok) THEN
									eseguito_correttamente_ese=.TRUE.
									EXIT loop_ese
								END IF
							END IF
							acc_ese=0
							cont_ese=1
						END IF				
					END DO loop_ese
				END IF
			END IF
		CASE ('1ppt')
			IF (flag_elettroni) THEN
				step_partenza_e=step_e
				step_partenza_se=step_se
				acc_e=0
				acc_se=0
				cont_e=1
				cont_se=1
				contacc_e=0
				contacc_se=0
				loop_1ppt: DO i = 1, MAX(N_blank,num_min)*N_1ppt, 1
					CALL RANDOM_NUMBER(eta)
					IF (flag_shadow) THEN
						ieta=CEILING(eta*3.d0*N_part)
					ELSE
						ieta=CEILING(eta*N_part)
					END IF
					IF (ieta<=N_part) THEN
						IF ( trimer_steps ) THEN
							tipo_mossa='tre'
						ELSE
							tipo_mossa='e__'
						END IF
						num_1ppt=ieta
						cont_e=cont_e+1
					ELSE IF (ieta<=2*N_part) THEN
						tipo_mossa='se1'
						num_1ppt=ieta-N_part
						cont_se=cont_se+1
					ELSE
						tipo_mossa='se2'
						num_1ppt=ieta-2*N_part
						cont_se=cont_se+1
					END IF
					CALL proponi_mossa(tipo_mossa,num_1ppt)
					CALL valuta_accettazione(i,tipo_mossa,num_1ppt,accettazione)
					IF (accettazione) THEN
						CALL metti_rnew_in_rold(tipo_mossa,num_1ppt)
						IF (ieta<=N_part) THEN
							acc_e=acc_e+1
						ELSE
							acc_se=acc_se+1
						END IF
					ELSE
						CALL riporta_indietro_walker(tipo_mossa,num_1ppt)
					END IF
					IF (MOD(cont_e,101)==0) THEN
						IF (acc_e>acceptance_rate+10) THEN
							contacc_e=0
							CALL change_step((REAL(acc_e-(acceptance_rate+10),8)/REAL(100-acceptance_rate-10,8)),'e__')
						ELSE IF (acc_e<acceptance_rate-10) THEN
							contacc_e=0
							CALL change_step((REAL(acc_e-(acceptance_rate-10),8)/REAL(acceptance_rate-10,8)),'e__')
						ELSE
							contacc_e=contacc_e+1
							IF (contacc_e==num_ok) THEN
								eseguito_correttamente_e=.TRUE.
								IF (eseguito_correttamente_se) THEN
									EXIT loop_1ppt
								END IF
							END IF
						END IF
						acc_e=0
						cont_e=1
					END IF
					IF (MOD(cont_se,101)==0) THEN
						IF (acc_se>acceptance_rate+10) THEN
							contacc_se=0
							CALL change_step((REAL(acc_se-(acceptance_rate+10),8)/REAL(100-acceptance_rate-10,8)),'se_')
						ELSE IF (acc_se<acceptance_rate-10) THEN
							contacc_se=0
							CALL change_step((REAL(acc_se-(acceptance_rate-10),8)/REAL(acceptance_rate-10,8)),'se_')
						ELSE
							contacc_se=contacc_se+1
							IF (contacc_se==num_ok) THEN
								eseguito_correttamente_se=.TRUE.
								IF (eseguito_correttamente_e) THEN
									EXIT loop_1ppt
								END IF
							END IF
						END IF
						acc_se=0
						cont_se=1
					END IF			
				END DO loop_1ppt
			END IF
		END SELECT
		!step_e_new=MIN(step_e,H_L)
		!CALL replace_step(step_e_new,'e__')
		!IF (flag_shadow) THEN
		!	step_se_new=MIN(step_se,H_L)
		!	CALL replace_step(step_se_new,'se_')
		!END IF
		IF ( trimer_steps ) THEN
			IF (verbose_mode .AND. eseguito_correttamente_e .AND. flag_elettroni) THEN
				PRINT * , 'VMC: Ho calibrato lo step_tre con successo: ', step_partenza_e, ' ---> ', step_e
				IF (flag_output) WRITE (7, *), 'VMC: Ho calibrato lo step_e con successo: ', step_partenza_e, ' ---> ', step_e
			ELSE IF (verbose_mode .AND. (.NOT. eseguito_correttamente_e)) THEN
				PRINT * , 'VMC: Ho calibrato lo step_tre (circa): ', step_partenza_e, ' ---> ', step_e
				IF (flag_output) WRITE (7, *), 'VMC: Ho calibrato lo step_e (circa): ', step_partenza_e, ' ---> ', step_e
			END IF
		ELSE
			IF (verbose_mode .AND. eseguito_correttamente_e .AND. flag_elettroni) THEN
				PRINT * , 'VMC: Ho calibrato lo step_e con successo: ', step_partenza_e, ' ---> ', step_e
				IF (flag_output) WRITE (7, *), 'VMC: Ho calibrato lo step_e con successo: ', step_partenza_e, ' ---> ', step_e
			ELSE IF (verbose_mode .AND. (.NOT. eseguito_correttamente_e)) THEN
				PRINT * , 'VMC: Ho calibrato lo step_e (circa): ', step_partenza_e, ' ---> ', step_e
				IF (flag_output) WRITE (7, *), 'VMC: Ho calibrato lo step_e (circa): ', step_partenza_e, ' ---> ', step_e
			END IF
		END IF
		IF (verbose_mode .AND. eseguito_correttamente_se .AND. flag_shadow) THEN
			PRINT * , 'VMC: Ho calibrato lo step_se con successo: ', step_partenza_se, ' ---> ', step_se
			IF (flag_output) WRITE (7, *), 'VMC: Ho calibrato lo step_se con successo: ', step_partenza_se, ' ---> ', step_se
		ELSE IF (verbose_mode .AND. (.NOT. eseguito_correttamente_se) .AND. flag_shadow) THEN
			PRINT * , 'VMC: Ho calibrato lo step_se (circa): ', step_partenza_se, ' ---> ', step_se
			IF (flag_output) WRITE (7, *), 'VMC: Ho calibrato lo step_se (circa): ', step_partenza_se, ' ---> ', step_se
		END IF
		IF (verbose_mode .AND. eseguito_correttamente_ese .AND. flag_shadow) THEN
			PRINT * , 'VMC: Ho calibrato bene entrambi gli step'
			IF (flag_output) WRITE (7, *), 'VMC: Ho calibrato bene entrambi gli step'
		ELSE IF (verbose_mode .AND. (.NOT. eseguito_correttamente_se) .AND. flag_shadow) THEN
			PRINT * , 'VMC: I due step non sono stati calibrati bene'
			IF (flag_output) WRITE (7, *), 'VMC: I due step non sono stati calibrati bene'
		END IF
		
		END IF cont
		
		success=.TRUE.
		IF (step_e<STEP_MIN) THEN
			success=.FALSE.
		ELSE
			IF ((flag_shadow).AND.(step_se<STEP_MIN)) success=.FALSE.
		END IF
		
		CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
		CALL MPI_REDUCE(success, mpi_success, 1, MPI_LOGICAL, MPI_LAND, 0, MPI_COMM_WORLD, mpi_ierr)
		IF (mpi_myrank==0) THEN
			success=mpi_success
		END IF
		CALL MPI_BCAST(success,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_ierr)
		CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
				
		CALL MPI_REDUCE(step_e,AV_step_e,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
		CALL MPI_BCAST(AV_step_e,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
		CALL replace_step(AV_step_e/REAL(mpi_nprocs),'e__')
		CALL MPI_REDUCE(step_se,AV_step_se,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
		CALL MPI_BCAST(AV_step_se,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
		CALL replace_step(AV_step_se/REAL(mpi_nprocs),'se_')
							
	END SUBROUTINE valuta_step_mc
!-----------------------------------------------------------------------
		
	SUBROUTINE sampling(val_estimatori,NUMBER_OF_SAMPLINGS)
		USE momenta
		IMPLICIT NONE
		LOGICAL, INTENT(IN) :: val_estimatori
		INTEGER (KIND=8), OPTIONAL :: NUMBER_OF_SAMPLINGS
		LOGICAL :: accettazione
		CHARACTER(LEN=3) :: tipo_mossa
		INTEGER :: ieta, num_1ppt 
		INTEGER (KIND=8) :: num_sampling, i, j, acc, i1
		INTEGER (KIND=8) :: cont_e(1:2), cont_se1(1:2), cont_se2(1:2)
		REAL (KIND=8) :: eta
				
		IF (.NOT. iniz_VMC) STOP 'Prima di eseguire il sampling devi inizializzare il VMC &
		  [ module_VMC.f90 > mosse_a_vuoto ]'
		
		cont_e=0
		cont_se1=0
		cont_se2=0
				
		IF (val_estimatori) THEN
			acc=1
			CALL calcola_nuovi_estimatori(acc)
			num_sampling=N_mc
		ELSE
			num_sampling=N_blank
		END IF
				
		IF (PRESENT(NUMBER_OF_SAMPLINGS)) num_sampling=NUMBER_OF_SAMPLINGS
		
		skip: IF (((.NOT. flag_continua) .AND. (.NOT. val_estimatori)) .OR. (val_estimatori) ) THEN
		
			IF ((mpi_myrank==0).AND.(val_estimatori).AND.(verbose_mode)) &
			  PRINT * , ' PROGRESS:   |0%       |20%      |40%      |60%      |80%      |100%'
			IF ((mpi_myrank==0).AND.(val_estimatori).AND.(verbose_mode)) &
			  WRITE (*, FMT='(A15)', ADVANCE='NO') , '              |'
					
			SELECT CASE (howtomove)
			CASE ('allp')
				DO i = 2, num_sampling, 1
					IF (flag_TABC .AND. (MOD(i,N_TABC)==0)) THEN     !applica TABC
						IF (SDe_kind=='pw_') CALL applica_twist(H_N_part, L)
						IF ((SDe_kind=='prf').OR.(SDe_kind=='fre')) CALL applica_twist_hartree()
						IF (SDe_kind=='lda') CALL applica_twist_lda()
						DO i1 = 1, N_mc_relax_TABC, 1
							IF (flag_elettroni) THEN
								CALL proponi_mossa('e__',-1)
								IF (flag_shadow) THEN
									CALL proponi_mossa('se_',-1)
									tipo_mossa='all'
								ELSE
									tipo_mossa='e__'
								END IF
							END IF
							CALL valuta_accettazione(i,tipo_mossa,-1,accettazione)
							IF (accettazione) THEN
								CALL metti_rnew_in_rold(tipo_mossa,-1)
							END IF
						END DO
					END IF
					IF ( flag_traccia_coppie .AND. (MOD(i,N_ritraccia_coppie)==0) .AND. (N_ritraccia_coppie>0)) THEN      !ritrova le coppie molecolari shadow-shadow
						CALL ritraccia_coppie_molecolari_ss()
						DO i1 = 1, N_mc_relax_traccia_coppie, 1
							IF (flag_elettroni) THEN
								CALL proponi_mossa('e__',-1)
								IF (flag_shadow) THEN
									CALL proponi_mossa('se_',-1)
									tipo_mossa='all'
								ELSE
									tipo_mossa='e__'
								END IF
							END IF
							CALL valuta_accettazione(i,tipo_mossa,-1,accettazione)
							IF (accettazione) THEN
								CALL metti_rnew_in_rold(tipo_mossa,-1)
							END IF
						END DO
					END IF
					IF (flag_elettroni) THEN
						CALL proponi_mossa('e__',-1)
						IF (flag_shadow) THEN
							CALL proponi_mossa('se_',-1)
							tipo_mossa='all'
						ELSE
							tipo_mossa='e__'
						END IF
					END IF
					CALL valuta_accettazione(i,tipo_mossa,-1,accettazione)
					IF (accettazione) THEN
						CALL metti_rnew_in_rold(tipo_mossa,-1)
						IF (val_estimatori) THEN
							CALL calcola_nuovi_estimatori(i)
							acc=acc+1
						END IF
					ELSE
						IF (val_estimatori) CALL conferma_estimatori(i)
					END IF
					IF (MOD(i,N_mc/50)==0 .AND. (mpi_myrank==0).AND.(val_estimatori).AND.(verbose_mode) .AND. (i<N_mc)) THEN
						WRITE (*, FMT='(A1)', ADVANCE='NO'), '='
					END IF
				END DO
				IF (val_estimatori) THEN
					IF ((mpi_myrank==0).AND.(val_estimatori).AND.(verbose_mode)) WRITE (*, FMT='(A1)', ADVANCE='YES'), '|'
					IF (verbose_mode) THEN
						PRINT * , 'VMC: Ho eseguito un campionamento MC con accettazione ', 100.*REAL(acc,4)/REAL(num_sampling,4), '%'
						IF (flag_output) WRITE (7, *), 'VMC: Ho eseguito un campionamento MC con accettazione ', &
						                   100.*REAL(acc,4)/REAL(num_sampling,4), '%'
					END IF
				ELSE
					IF (verbose_mode) THEN
						PRINT * , 'VMC: Ho eseguito ', num_sampling ,' mosse a vuoto con metodo allp'
						IF (flag_output) WRITE (7, *), 'VMC: Ho eseguito ', num_sampling ,' mosse a vuoto con metodo allp'
					END IF
				END IF
			CASE ('1ppt')
				DO i = 2, num_sampling, 1
					IF ((MOD(i,N_TABC)==0) .AND. flag_TABC) THEN
						IF (SDe_kind=='pw_') CALL applica_twist(H_N_part, L)
						IF ((SDe_kind=='prf').OR.(SDe_kind=='fre')) CALL applica_twist_hartree()
						IF (SDe_kind=='lda') CALL applica_twist_lda()
						DO i1 = 1, N_mc_relax_TABC, 1
							DO j = 1, N_1ppt, 1
								IF (flag_elettroni) THEN
									CALL RANDOM_NUMBER(eta)
									IF (flag_shadow) THEN
										ieta=CEILING(eta*3.d0*N_part)
									ELSE
										ieta=CEILING(eta*N_part)
									END IF
								END IF
								IF (ieta<=N_part) THEN
									IF ( trimer_steps ) THEN
										tipo_mossa='tre'
									ELSE
										tipo_mossa='e__'
									END IF
									num_1ppt=ieta
								ELSE IF (ieta<=2*N_part) THEN
									tipo_mossa='se1'
									num_1ppt=ieta-N_part
								ELSE
									tipo_mossa='se2'
									num_1ppt=ieta-2*N_part
								END IF
								CALL proponi_mossa(tipo_mossa,num_1ppt)
								CALL valuta_accettazione(i,tipo_mossa,num_1ppt,accettazione)
								IF (accettazione) THEN
									CALL metti_rnew_in_rold(tipo_mossa,num_1ppt)
								ELSE
									CALL riporta_indietro_walker(tipo_mossa,num_1ppt)
								END IF
							END DO
						END DO
					END IF
					IF ( flag_traccia_coppie .AND. (N_ritraccia_coppie>0) .AND. (MOD(i,N_ritraccia_coppie)==0)) THEN      !ritrova le coppie molecolari shadow-shadow
						CALL ritraccia_coppie_molecolari_ss()
						DO i1 = 1, N_mc_relax_traccia_coppie, 1
							DO j = 1, N_1ppt, 1
								CALL RANDOM_NUMBER(eta)
								ieta=CEILING(eta*2.d0*N_part)
								IF (ieta<=N_part) THEN
									tipo_mossa='se1'
									num_1ppt=ieta
								ELSE IF (ieta<=2*N_part) THEN
									tipo_mossa='se2'
									num_1ppt=ieta-N_part
								END IF
								CALL proponi_mossa(tipo_mossa,num_1ppt)
								CALL valuta_accettazione(i,tipo_mossa,num_1ppt,accettazione)
								IF (accettazione) THEN
									CALL metti_rnew_in_rold(tipo_mossa,num_1ppt)
								ELSE
									CALL riporta_indietro_walker(tipo_mossa,num_1ppt)
								END IF
							END DO
						END DO
					END IF
					DO j = 1, N_1ppt, 1
						IF (flag_elettroni) THEN
							CALL RANDOM_NUMBER(eta)
							IF (flag_shadow) THEN
								ieta=CEILING(eta*3.d0*N_part)
							ELSE
								ieta=CEILING(eta*N_part)
							END IF
						END IF
						IF (ieta<=N_part) THEN
							IF ( trimer_steps ) THEN
								tipo_mossa='tre'
							ELSE
								tipo_mossa='e__'
							END IF
							num_1ppt=ieta
							cont_e(1)=cont_e(1)+1
						ELSE IF (ieta<=2*N_part) THEN
							tipo_mossa='se1'
							num_1ppt=ieta-N_part
							cont_se1(1)=cont_se1(1)+1
						ELSE
							tipo_mossa='se2'
							num_1ppt=ieta-2*N_part
							cont_se2(1)=cont_se2(1)+1
						END IF
						CALL proponi_mossa(tipo_mossa,num_1ppt)
						CALL valuta_accettazione(i,tipo_mossa,num_1ppt,accettazione)
						IF (accettazione) THEN
							CALL metti_rnew_in_rold(tipo_mossa,num_1ppt)
							IF (val_estimatori) acc=acc+1
							IF (ieta<=1*N_part) THEN
								cont_e(2)=cont_e(2)+1
							ELSE IF (ieta<=2*N_part) THEN
								cont_se1(2)=cont_se1(2)+1
							ELSE
								cont_se2(2)=cont_se2(2)+1
							END IF
						ELSE
							CALL riporta_indietro_walker(tipo_mossa,num_1ppt)
						END IF
					END DO
					IF (val_estimatori) CALL calcola_nuovi_estimatori(i)
					IF (MOD(i,N_mc/50)==0 .AND. ((mpi_myrank==0).AND.(val_estimatori).AND.(verbose_mode)) .AND. (i<N_mc)) THEN
						WRITE (*, FMT='(A1)', ADVANCE='NO'), '='
					END IF
				END DO
				IF (verbose_mode .AND. (val_estimatori)) THEN
					IF ((mpi_myrank==0).AND.(val_estimatori)) WRITE (*, FMT='(A1)', ADVANCE='YES'), '|'
					PRINT * , 'VMC: Ho eseguito un campionamento MC 1ppt su ', mpi_nprocs,'processori per ',N_mc, &
					  'mosse ognuno, con accettazione ', &
				 	  100.*REAL(acc,4)/REAL(N_mc*N_1ppt,4), '%'
					IF (flag_output) WRITE (7, *), 'VMC: Ho eseguito un campionamento MC 1ppt su ', mpi_nprocs, &
					  'processori per ',N_mc, 'mosse ognuno, con accettazione ', 100.*REAL(acc,4)/REAL(N_mc*N_1ppt,4), '%'
				ELSE IF (verbose_mode .AND. (.NOT. val_estimatori)) THEN
					PRINT * , 'VMC: Ho eseguito ', N_blank ,' per ', N_1ppt, ' mosse a vuoto con metodo ', howtomove
					IF (flag_output) WRITE (7, *), 'VMC: Ho eseguito ', N_blank ,' per ', N_1ppt, &
					  ' mosse a vuoto con metodo ', howtomove
				END IF
				IF (verbose_mode .AND. val_estimatori) THEN
					PRINT * , 'VMC: accettazione e: ', 100.*REAL(cont_e(2))/REAL(cont_e(1))
					IF (flag_output) WRITE (7, *), 'VMC: accettazione e: ', 100.*REAL(cont_e(2))/REAL(cont_e(1))
					IF (verbose_mode .AND. flag_shadow) THEN
						PRINT * , 'VMC: accettazione se1: ', 100.*REAL(cont_se1(2))/REAL(cont_se1(1))
						PRINT * , 'VMC: accettazione se2: ', 100.*REAL(cont_se2(2))/REAL(cont_se2(1))
						IF (flag_output) THEN
							WRITE (7, *), 'VMC: accettazione se1: ', 100.*REAL(cont_se1(2))/REAL(cont_se1(1))
							WRITE (7, *), 'VMC: accettazione se2: ', 100.*REAL(cont_se2(2))/REAL(cont_se2(1))
						END IF
					END IF
				END IF
				!!!PRINT * , mpi_myrank, '   - - -   VMC: accettazione e: ', 100.*REAL(cont_e(2))/REAL(cont_e(1))
			END SELECT
			IF (val_estimatori .AND. flag_posizioni .AND. flag_disk) CALL salva_posizione_walkers('posizioni/fine'//codice_simulazione)
		
		END IF skip
		
	END SUBROUTINE sampling
!-----------------------------------------------------------------------

	SUBROUTINE calcola_nuovi_estimatori(i)
		USE variational_calculations
		IMPLICIT NONE
		INTEGER :: m, n
		INTEGER (KIND=8), INTENT(IN) :: i
		INTEGER :: j
		IF (.NOT. iniz_VMC) STOP 'Prima di valutare gli estimatori devi inizializzare il VMC &
		  [ module_VMC.f90 > calcola_nuovi_estimatori ]'
		
		CALL valuta_estimatori(i)
				
		IF (flag_gr) THEN
			CALL g_ii_correlation(rij_ee_old,N_part,L,r_s,N_hist,geeuu_new,geeud_new)
			geeuu=geeuu+geeuu_new
			geeud=geeud+geeud_new
			CALL g_ij_correlation(rij_ep_old,N_part,L,r_s,N_hist,gep_new)
			gep=gep+gep_new
			CALL g_ii_correlation(rij_pp_old,N_part,L,r_s,N_hist,gppuu_new,gppud_new)
			gppuu=gppuu+gppuu_new
			gppud=gppud+gppud_new
			IF (flag_shadow) THEN
				CALL g_ii_correlation(rij_se1_old,N_part,L,r_s,N_hist,gseuu1_new,gseud1_new)
				gseuu1=gseuu1+gseuu1_new
				gseud1=gseud1+gseud1_new
				CALL g_ij_correlation(rij_ese1_old,N_part,L,r_s,N_hist,gese1_new)
				gese1=gese1+gese1_new
				CALL g_ij_correlation(rij_sesp1_old,N_part,L,r_s,N_hist,gsesp1_new)
				gsesp1=gsesp1+gsesp1_new
				CALL g_ii_correlation(rij_se2_old,N_part,L,r_s,N_hist,gseuu2_new,gseud2_new)
				gseuu2=gseuu2+gseuu2_new
				gseud2=gseud2+gseud2_new
				CALL g_ij_correlation(rij_ese2_old,N_part,L,r_s,N_hist,gese2_new)
				gese2=gese2+gese2_new
				CALL g_ij_correlation(rij_sesp2_old,N_part,L,r_s,N_hist,gsesp2_new)
				gsesp2=gsesp2+gsesp2_new
			END IF
		END IF
		
		IF (flag_gradiente) THEN
			CALL calcola_estimatori_per_gradiente(i)
		END IF
		
		IF (flag_derivate_var) THEN
			CALL calcola_estimatori_per_derivate_variazionali(i)
		END IF
		
	END SUBROUTINE calcola_nuovi_estimatori
!-----------------------------------------------------------------------	

	SUBROUTINE conferma_estimatori(i)
		USE variational_calculations
		IMPLICIT NONE
		INTEGER (KIND=8), INTENT(IN) :: i
		IF (.NOT. iniz_VMC) STOP 'Prima di valutare gli estimatori devi inizializzare il VMC &
		  [ module_VMC.f90 > conferma_estimatori ]'
		
		CALL mantieni_stessi_estimatori(i)
		
		IF (flag_gr) THEN
			geeuu=geeuu+geeuu_new
			geeud=geeud+geeud_new
			gep=gep+gep_new
			gppuu=gppuu+gppuu_new
			gppud=gppud+gppud_new
			IF (flag_shadow) THEN
				gseuu1=gseuu1+gseuu1_new
				gseud1=gseud1+gseud1_new
				gese1=gese1+gese1_new
				gsesp1=gsesp1+gsesp1_new
				gseuu2=gseuu2+gseuu2_new
				gseud2=gseud2+gseud2_new
				gese2=gese2+gese2_new
				gsesp2=gsesp2+gsesp2_new
			END IF
		END IF
		
		IF (flag_gradiente) THEN
			CALL conferma_estimatori_per_gradiente(i)
		END IF
		
		IF (flag_derivate_var) THEN
			CALL conferma_estimatori_per_derivate_variazionali(i)
		END IF
		
	END SUBROUTINE conferma_estimatori
!-----------------------------------------------------------------------

	SUBROUTINE trascrivi_dati()
		USE variational_calculations
		IMPLICIT NONE
		LOGICAL :: flag_scrivi, flag_continua_rank
		INTEGER :: status(MPI_STATUS_SIZE)
		INTEGER :: j
		INTEGER (KIND=8) :: i
		REAL (KIND=8) :: x(1:N_hist), gr_sum(1:N_hist)
				
		IF (.NOT. iniz_VMC) STOP 'Prima di salvare i dati devi inizializzare il VMC &
		  [ module_VMC.f90 > conferma_estimatori ]'
				
		IF (mpi_myrank==0) THEN
			flag_scrivi=.TRUE.
			flag_continua_rank=flag_continua
		ELSE
			flag_scrivi=.FALSE.
			flag_continua_rank=.TRUE.
		END IF
		
		IF (flag_disk) THEN
			DO j = 0, mpi_nprocs-1, 1
				IF (mpi_myrank==j .AND. flag_scrivi) THEN
					IF (mpi_myrank<mpi_nprocs-1) THEN
						CALL MPI_SEND(flag_scrivi,1,MPI_LOGICAL,mpi_myrank+1,10+j,MPI_COMM_WORLD,mpi_ierr)
						IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_SEND &
						  [ module_VMC.f90 > trascrivi_dati ]'
					END IF
				END IF
				IF (mpi_myrank==j+1 .AND. j<mpi_nprocs) THEN
					CALL MPI_RECV(flag_scrivi,1,MPI_LOGICAL,mpi_myrank-1,10+j,MPI_COMM_WORLD,status,mpi_ierr)
					IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_RECV &
					  [ module_VMC.f90 > trascrivi_dati ]'
				END IF
			END DO
			CALL trascrivi_estimatori(codice_simulazione)
			IF (mpi_myrank==0) THEN
				PRINT * , 'VMC: Ho trascritto i dati'
				IF (flag_output) WRITE (7, *), 'VMC: Ho trascritto i dati'
			END IF
			IF (flag_gradiente) CALL trascrivi_dati_per_gradiente()
			IF (flag_derivate_var) CALL trascrivi_dati_per_derivate_variazionali()
		END IF
		
		IF ( flag_disk .AND. flag_gr) THEN
			j=N_hist
			CALL MPI_REDUCE(geeuu(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
			IF (mpi_myrank==0) geeuu=gr_sum/REAL(N_mc*mpi_nprocs,8)
			CALL MPI_REDUCE(geeud(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
			IF (mpi_myrank==0) geeud=gr_sum/REAL(N_mc*mpi_nprocs,8)
			CALL MPI_REDUCE(gep(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
			IF (mpi_myrank==0) gep=gr_sum/REAL(N_mc*mpi_nprocs,8)
			CALL MPI_REDUCE(gppuu(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
			IF (mpi_myrank==0) gppuu=gr_sum/REAL(N_mc*mpi_nprocs,8)
			CALL MPI_REDUCE(gppud(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
			IF (mpi_myrank==0) gppud=gr_sum/REAL(N_mc*mpi_nprocs,8)
			IF ( flag_shadow ) THEN
				CALL MPI_REDUCE(gseuu1(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
				IF (mpi_myrank==0) gseuu1=gr_sum/REAL(N_mc*mpi_nprocs,8)
				CALL MPI_REDUCE(gseud1(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
				IF (mpi_myrank==0) gseud1=gr_sum/REAL(N_mc*mpi_nprocs,8)
				CALL MPI_REDUCE(gese1(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
				IF (mpi_myrank==0) gese1=0.5d0*gr_sum/REAL(N_mc*mpi_nprocs,8)
				CALL MPI_REDUCE(gsesp1(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
				IF (mpi_myrank==0) gsesp1=0.5d0*gr_sum/REAL(N_mc*mpi_nprocs,8)
				CALL MPI_REDUCE(gseuu2(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
				IF (mpi_myrank==0) gseuu2=gr_sum/REAL(N_mc*mpi_nprocs,8)
				CALL MPI_REDUCE(gseud2(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
				IF (mpi_myrank==0) gseud2=gr_sum/REAL(N_mc*mpi_nprocs,8)
				CALL MPI_REDUCE(gese2(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
				IF (mpi_myrank==0) gese2=0.5d0*gr_sum/REAL(N_mc*mpi_nprocs,8)
				CALL MPI_REDUCE(gsesp2(1:N_hist),gr_sum(1:N_hist),j,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
				IF (mpi_myrank==0) gsesp2=0.5d0*gr_sum/REAL(N_mc*mpi_nprocs,8)
			END IF
			
			IF ( mpi_myrank==0 ) THEN
				x(1:N_hist-1)=(/ ( REAL(i,8)*0.5d0*MINVAL(L)/REAL(N_hist,8) , i=2,N_hist ) /)
				CALL salva_vettore_dati_dat(x,geeuu(2:N_hist),N_hist-1,'estimatori/geeuu'//codice_simulazione)
				CALL salva_vettore_dati_dat(x,geeud(2:N_hist),N_hist-1,'estimatori/geeud'//codice_simulazione)
				CALL salva_vettore_dati_dat(x,geeuu(2:N_hist)+geeud(2:N_hist),N_hist-1,'estimatori/gee'//codice_simulazione)
				CALL salva_vettore_dati_dat(x,gep(2:N_hist),N_hist-1,'estimatori/gep'//codice_simulazione)
				CALL salva_vettore_dati_dat(x,gppuu(2:N_hist),N_hist-1,'estimatori/gppuu'//codice_simulazione)
				CALL salva_vettore_dati_dat(x,gppud(2:N_hist),N_hist-1,'estimatori/gppud'//codice_simulazione)
				CALL salva_vettore_dati_dat(x,gppuu(2:N_hist)+gppud(2:N_hist),N_hist-1,'estimatori/gpp'//codice_simulazione)
				IF (flag_shadow) THEN
					CALL salva_vettore_dati_dat(x,gseuu1(2:N_hist),N_hist-1,'estimatori/gseuu1'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gseud1(2:N_hist),N_hist-1,'estimatori/gseud1'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gseuu1(2:N_hist)+gseud1(2:N_hist),N_hist-1,'estimatori/gse1'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gese1(2:N_hist),N_hist-1,'estimatori/gese1'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gsesp1(2:N_hist),N_hist-1,'estimatori/gsesp1'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gseuu2(2:N_hist),N_hist-1,'estimatori/gseuu2'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gseud2(2:N_hist),N_hist-1,'estimatori/gseud2'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gseuu2(2:N_hist)+gseud2(2:N_hist),N_hist-1,'estimatori/gse2'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gese2(2:N_hist),N_hist-1,'estimatori/gese2'//codice_simulazione)
					CALL salva_vettore_dati_dat(x,gsesp2(2:N_hist),N_hist-1,'estimatori/gsesp2'//codice_simulazione)
				END IF
			END IF
		END IF
				
	END SUBROUTINE trascrivi_dati
!-----------------------------------------------------------------------

	SUBROUTINE stampa_risultati()
		IMPLICIT NONE
		REAL (KIND=8) :: media, errore, estimatore(1:2)
		
		IF (quick_error==0) THEN
			IF (flag_disk) THEN
				!IF (SDse_kind=='no_') THEN
				!	IF ((flag_mpi .AND. mpi_myrank==0) .OR. ( .NOT. flag_mpi)) THEN
				!	IF (flag_E_kin) THEN
				!		CALL calcola_estimatori('estimatori/E_kin'//codice_simulazione,media,errore)
				!		PRINT * , '   --->   E_kin=', media, '+-', errore
				!		IF (flag_output) WRITE (7, *), '   --->   E_kin=', media, '+-', errore
				!		CALL calcola_estimatori('estimatori/E_JF'//codice_simulazione,media,errore)
				!		PRINT * , '   --->   E_JF= ', media, '+-', errore
				!		IF (flag_output) WRITE (7, *), '   --->   E_JF= ', media, '+-', errore
				!	END IF
				!	IF (flag_E_pot) THEN
				!		CALL calcola_estimatori('estimatori/E_pot'//codice_simulazione,media,errore)
				!		PRINT * , '   --->   E_pot=', media, '+-', errore
				!		IF (flag_output) WRITE (7, *), '   --->   E_pot=', media, '+-', errore
				!	END IF
				!	IF (flag_E_tot) THEN
				!		CALL calcola_estimatori('estimatori/E_tot'//codice_simulazione,media,errore)
				!		PRINT * , '   --->   E_tot=', media, '+-', errore
				!		IF (flag_output) WRITE (7, *), '   --->   E_tot=', media, '+-', errore
				!	END IF
				!	END IF
				!ELSE
					IF ((flag_mpi .AND. mpi_myrank==0) .OR. ( .NOT. flag_mpi)) THEN
					IF (flag_E_kin) THEN
						CALL calcola_estimatori('estimatori/E_kin'//codice_simulazione,media,errore,'estimatori/w'//codice_simulazione)
						PRINT * , '   --->   E_kin=', media, '+-', errore
						IF (flag_output) WRITE (7, *), '   --->   E_kin=', media, '+-', errore
						CALL calcola_estimatori('estimatori/E_JF'//codice_simulazione,media,errore,'estimatori/w'//codice_simulazione)
						PRINT * , '   --->   E_JF= ', media, '+-', errore
						IF (flag_output) WRITE (7, *), '   --->   E_JF= ', media, '+-', errore
					END IF
					IF (flag_E_pot) THEN
						CALL calcola_estimatori('estimatori/E_pot'//codice_simulazione,media,errore,'estimatori/w'//codice_simulazione)
						PRINT * , '   --->   E_pot=', media, '+-', errore
						IF (flag_output) WRITE (7, *), '   --->   E_pot=', media, '+-', errore
					END IF
					IF (flag_E_tot) THEN
						CALL calcola_estimatori('estimatori/E_tot'//codice_simulazione,media,errore,'estimatori/w'//codice_simulazione)
						PRINT * , '   --->   E_tot=', media, '+-', errore
						IF (flag_output) WRITE (7, *), '   --->   E_tot=', media, '+-', errore
						variance_efficiency=errore*errore
					END IF
					END IF
				!END IF
			ELSE
				!IF (SDse_kind=='no_') THEN
				!	IF (flag_E_kin) THEN
				!		CALL calcola_estimatore_da_RAM(estimatore,E_kin,N_mc)
				!		IF (mpi_myrank==0) THEN
				!			PRINT * , '   --->   E_kin=', estimatore(1), '+-', estimatore(2)
				!			IF (flag_output) WRITE (7, *), '   --->   E_kin=', estimatore(1), '+-', estimatore(2)
				!		END IF
				!		CALL calcola_estimatore_da_RAM(estimatore,E_JF,N_mc)
				!		IF (mpi_myrank==0) THEN
				!			PRINT * , '   --->   E_JF= ', estimatore(1), '+-', estimatore(2)
				!			IF (flag_output) WRITE (7, *), '   --->   E_JF= ', estimatore(1), '+-', estimatore(2)
				!		END IF
				!	END IF
				!	IF (flag_E_pot) THEN
				!		CALL calcola_estimatore_da_RAM(estimatore,E_pot,N_mc)
				!		IF (mpi_myrank==0) PRINT * , '   --->   E_pot=', estimatore(1), '+-', estimatore(2)
				!		IF (flag_output) WRITE (7, *), '   --->   E_pot=', estimatore(1), '+-', estimatore(2)
				!	END IF
				!	IF (flag_E_tot) THEN
				!		CALL calcola_estimatore_da_RAM(estimatore,E_tot,N_mc)
				!		IF (mpi_myrank==0) THEN
				!			PRINT * , '   --->   E_tot=', estimatore(1), '+-', estimatore(2)
				!			IF (flag_output) WRITE (7, *), '   --->   E_tot=', estimatore(1), '+-', estimatore(2)
				!		END IF
				!	END IF
				!ELSE
					IF (flag_E_kin) THEN
						CALL calcola_estimatore_da_RAM(estimatore,E_kin,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , '   --->   E_kin=', estimatore(1), '+-', estimatore(2)
							IF (flag_output) WRITE (7, *), '   --->   E_kin=', estimatore(1), '+-', estimatore(2)
						END IF
						CALL calcola_estimatore_da_RAM(estimatore,E_JF,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , '   --->   E_JF= ', estimatore(1), '+-', estimatore(2)
							IF (flag_output) WRITE (7, *), '   --->   E_JF= ', estimatore(1), '+-', estimatore(2)
						END IF
					END IF
					IF (flag_E_pot) THEN
						CALL calcola_estimatore_da_RAM(estimatore,E_pot,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , '   --->   E_pot=', estimatore(1), '+-', estimatore(2)
							IF (flag_output) WRITE (7, *), '   --->   E_pot=', estimatore(1), '+-', estimatore(2)
						END IF
					END IF
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM(estimatore,E_tot,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , '   --->   E_tot=', estimatore(1), '+-', estimatore(2)
							IF (flag_output) WRITE (7, *), '   --->   E_tot=', estimatore(1), '+-', estimatore(2)
						END IF
						variance_efficiency=estimatore(2)*estimatore(2)
					END IF
				!END IF
			END IF
		ELSE IF (quick_error>3) THEN
			IF (flag_disk) THEN
				!IF (SDse_kind=='no_') THEN
				!	IF ((flag_mpi .AND. mpi_myrank==0) .OR. ( .NOT. flag_mpi)) THEN
				!	IF (flag_E_kin) THEN
				!		CALL calcola_estimatori_quick(quick_error,'estimatori/E_kin'//codice_simulazione,media,errore)
				!		PRINT * , '   --->   E_kin=', media, '+-', errore
				!		IF (flag_output) WRITE (7, *), '   --->   E_kin=', media, '+-', errore
				!		CALL calcola_estimatori_quick(quick_error,'estimatori/E_JF'//codice_simulazione,media,errore)
				!		PRINT * , '   --->   E_JF= ', media, '+-', errore
				!		IF (flag_output) WRITE (7, *), '   --->   E_JF= ', media, '+-', errore
				!	END IF
				!	IF (flag_E_pot) THEN
				!		CALL calcola_estimatori_quick(quick_error,'estimatori/E_pot'//codice_simulazione,media,errore)
				!		PRINT * , '   --->   E_pot=', media, '+-', errore
				!		IF (flag_output) WRITE (7, *), '   --->   E_pot=', media, '+-', errore
				!	END IF
				!	IF (flag_E_tot) THEN
				!		CALL calcola_estimatori_quick(quick_error,'estimatori/E_tot'//codice_simulazione,media,errore)
				!		PRINT * , '   --->   E_tot=', media, '+-', errore
				!		IF (flag_output) WRITE (7, *), '   --->   E_tot=', media, '+-', errore
				!	END IF
				!	END IF
				!ELSE
					IF ((flag_mpi .AND. mpi_myrank==0) .OR. ( .NOT. flag_mpi)) THEN
					IF (flag_E_kin) THEN
						CALL calcola_estimatori_quick(quick_error,'estimatori/E_kin'//codice_simulazione,media,errore, &
						  'estimatori/w'//codice_simulazione)
						PRINT * , '   --->   E_kin=', media, '+-', errore
						IF (flag_output) WRITE (7, *), '   --->   E_kin=', media, '+-', errore
						CALL calcola_estimatori_quick(quick_error,'estimatori/E_JF'//codice_simulazione,media,errore, &
						  'estimatori/w'//codice_simulazione)
						PRINT * , '   --->   E_JF= ', media, '+-', errore
						IF (flag_output) WRITE (7, *), '   --->   E_JF= ', media, '+-', errore
					END IF
					IF (flag_E_pot) THEN
						CALL calcola_estimatori_quick(quick_error,'estimatori/E_pot'//codice_simulazione,media,errore, &
						  'estimatori/w'//codice_simulazione)
						PRINT * , '   --->   E_pot=', media, '+-', errore
						IF (flag_output) WRITE (7, *), '   --->   E_pot=', media, '+-', errore
					END IF
					IF (flag_E_tot) THEN
						CALL calcola_estimatori_quick(quick_error,'estimatori/E_tot'//codice_simulazione,media,errore, &
						  'estimatori/w'//codice_simulazione)
						PRINT * , '   --->   E_tot=', media, '+-', errore
						IF (flag_output) WRITE (7, *), '   --->   E_tot=', media, '+-', errore
						variance_efficiency=errore*errore
					END IF
					END IF
				!END IF
			ELSE
				!IF (SDse_kind=='no_') THEN
				!	IF (flag_E_kin) THEN
				!		CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,E_kin,N_mc)
				!		IF (mpi_myrank==0) THEN
				!			PRINT * , '   --->   E_kin=', estimatore(1), '+-', estimatore(2)
				!			IF (flag_output) WRITE (7, *), '   --->   E_kin=', estimatore(1), '+-', estimatore(2)
				!		END IF
				!		CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,E_JF,N_mc)
				!		IF (mpi_myrank==0) THEN
				!			PRINT * , '   --->   E_JF= ', estimatore(1), '+-', estimatore(2)
				!			IF (flag_output) WRITE (7, *), '   --->   E_JF= ', estimatore(1), '+-', estimatore(2)
				!		END IF
				!	END IF
				!	IF (flag_E_pot) THEN
				!		CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,E_pot,N_mc)
				!		IF (mpi_myrank==0) PRINT * , '   --->   E_pot=', estimatore(1), '+-', estimatore(2)
				!		IF (flag_output) WRITE (7, *), '   --->   E_pot=', estimatore(1), '+-', estimatore(2)
				!	END IF
				!	IF (flag_E_tot) THEN
				!		CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,E_tot,N_mc)
				!		IF (mpi_myrank==0) THEN
				!			PRINT * , '   --->   E_tot=', estimatore(1), '+-', estimatore(2)
				!			IF (flag_output) WRITE (7, *), '   --->   E_tot=', estimatore(1), '+-', estimatore(2)
				!		END IF
				!	END IF
				!ELSE
					IF (flag_E_kin) THEN
						CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,E_kin,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , '   --->   E_kin=', estimatore(1), '+-', estimatore(2)
							IF (flag_output) WRITE (7, *), '   --->   E_kin=', estimatore(1), '+-', estimatore(2)
						END IF
						CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,E_JF,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , '   --->   E_JF= ', estimatore(1), '+-', estimatore(2)
							IF (flag_output) WRITE (7, *), '   --->   E_JF= ', estimatore(1), '+-', estimatore(2)
						END IF
					END IF
					IF (flag_E_pot) THEN
						CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,E_pot,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , '   --->   E_pot=', estimatore(1), '+-', estimatore(2)
							IF (flag_output) WRITE (7, *), '   --->   E_pot=', estimatore(1), '+-', estimatore(2)
						END IF
					END IF
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,E_tot,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , '   --->   E_tot=', estimatore(1), '+-', estimatore(2)
							IF (flag_output) WRITE (7, *), '   --->   E_tot=', estimatore(1), '+-', estimatore(2)
						END IF
						variance_efficiency=estimatore(2)*estimatore(2)
					END IF
				!END IF
			END IF
		ELSE
			STOP 'Il parametro quick_error non é stato scelto opportunamente'
		END IF
		
	END SUBROUTINE stampa_risultati
!-----------------------------------------------------------------------
	SUBROUTINE restituisci_energia(E_tot_out)
		IMPLICIT NONE
		REAL (KIND=8) :: media, errore, estimatore(1:2)
		REAL (KIND=8), INTENT(OUT) :: E_tot_out(1:2)
		
		IF (quick_error==0) THEN
			IF (flag_disk) THEN
				IF (SDse_kind=='no_') THEN
					IF (flag_E_tot) THEN
						IF (mpi_myrank==0) THEN
							CALL calcola_estimatori('estimatori/E_tot'//codice_simulazione,media,errore)
							IF (mpi_myrank==0) THEN
								PRINT * , 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
								IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
							END IF
							E_tot_out(1:2)=(/media,errore/)
						END IF
						CALL MPI_BCAST(E_tot_out(1:2),2,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				ELSE
					IF (flag_E_tot) THEN
						IF (mpi_myrank==0) THEN
							CALL calcola_estimatori('estimatori/E_tot'//codice_simulazione,media,errore,'estimatori/w'//codice_simulazione)
							IF (mpi_myrank==0) THEN
								PRINT * , 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
								IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
							END IF
							E_tot_out(1:2)=(/media,errore/)
						END IF
						CALL MPI_BCAST(E_tot_out(1:2),2,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				END IF
			ELSE
				IF (SDse_kind=='no_') THEN
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM(E_tot_out,E_tot,N_mc)
						IF (mpi_myrank==0) THEN
							PRINT * , 'Ho calcolato l energia --->   E_tot=', E_tot_out(1), '+-', E_tot_out(2)
							IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', E_tot_out(1), '+-', E_tot_out(2)
						END IF
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				ELSE
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM(E_tot_out,E_tot,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , 'Ho calcolato l energia --->   E_tot=', E_tot_out(1), '+-', E_tot_out(2)
							IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', E_tot_out(1), '+-', E_tot_out(2)
						END IF
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				END IF
			END IF
		ELSE IF (quick_error>3) THEN
			IF (flag_disk) THEN
				IF (SDse_kind=='no_') THEN
					IF (flag_E_tot) THEN
						IF (mpi_myrank==0) THEN
							CALL calcola_estimatori_quick(quick_error,'estimatori/E_tot'//codice_simulazione,media,errore)
							IF (mpi_myrank==0) THEN
								PRINT * , 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
								IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
							END IF
							E_tot_out(1:2)=(/media,errore/)
						END IF
						CALL MPI_BCAST(E_tot_out(1:2),2,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				ELSE
					IF (flag_E_tot) THEN
						IF (mpi_myrank==0) THEN
							CALL calcola_estimatori_quick(quick_error,'estimatori/E_tot'//codice_simulazione,media,errore, &
							  'estimatori/w'//codice_simulazione)
							IF (mpi_myrank==0) THEN
								PRINT * , 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
								IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
							END IF
							E_tot_out(1:2)=(/media,errore/)
						END IF
						CALL MPI_BCAST(E_tot_out(1:2),2,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				END IF
			ELSE
				IF (SDse_kind=='no_') THEN
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM_quick(quick_error,E_tot_out,E_tot,N_mc)
						IF (mpi_myrank==0) THEN
							PRINT * , 'Ho calcolato l energia --->   E_tot=', E_tot_out(1), '+-', E_tot_out(2)
							IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', E_tot_out(1), '+-', E_tot_out(2)
						END IF
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				ELSE
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM_quick(quick_error,E_tot_out,E_tot,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , 'Ho calcolato l energia --->   E_tot=', E_tot_out(1), '+-', E_tot_out(2)
							IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', E_tot_out(1), '+-', E_tot_out(2)
						END IF
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				END IF
			END IF
		ELSE
			STOP 'Il parametro quick_error non é stato scelto opportunamente'
		END IF
			
	END SUBROUTINE restituisci_energia
!-----------------------------------------------------------------------
	SUBROUTINE restituisci_energie(E_tot_out)        !restituisce le differenze di energie fra E e E'
		USE variational_calculations
		IMPLICIT NONE
		CHARACTER(LEN=4) :: codice_numerico
		INTEGER :: i_mc
		REAL (KIND=8) :: media, errore, dummy1(1:N_mc), dummy2(1:N_mc), dummy3(1:N_mc), dummy4(1:N_mc)
		REAL (KIND=8), INTENT(OUT) :: E_tot_out(1:2,0:num_wf)
		INTEGER :: i
		
		IF (quick_error==0) THEN
			IF (flag_disk) THEN
				IF (SDse_kind=='no_') THEN
					IF (flag_E_tot) THEN
						IF (mpi_myrank==0) THEN
							CALL calcola_estimatori('estimatori/E_tot'//codice_simulazione,media,errore)
							IF (mpi_myrank==0) THEN
								PRINT * , 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
								IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
							END IF
							E_tot_out(1:2,0)=(/media,errore/)
							DO i = 1, num_wf, 1
								WRITE (codice_numerico, '(I4.4)'), i
								CALL calcola_estimatori_differenza('estimatori/gradiente/E_tot'//codice_numerico,'estimatori/E_tot'//codice_simulazione, &
								  media,errore,'estimatori/gradiente/w'//codice_numerico)
								E_tot_out(1:2,i)=(/media,errore/)
							END DO
						END IF
						CALL MPI_BCAST(E_tot_out(1:2,0:num_wf),2*(num_wf+1),MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				ELSE
					IF (flag_E_tot) THEN
						IF (mpi_myrank==0) THEN
							CALL calcola_estimatori('estimatori/E_tot'//codice_simulazione,media,errore,'estimatori/w'//codice_simulazione)
							IF (mpi_myrank==0) THEN
								PRINT * , 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
								IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
							END IF
							E_tot_out(1:2,0)=(/media,errore/)
							DO i = 1, num_wf, 1
								WRITE (codice_numerico, '(I4.4)'), i
								CALL calcola_estimatori_differenza('estimatori/gradiente/E_tot'//codice_numerico,'estimatori/E_tot'//codice_simulazione, &
								  media,errore,'estimatori/gradiente/w'//codice_numerico,'estimatori/w'//codice_simulazione)
								E_tot_out(1:2,i)=(/media,errore/)
							END DO
						END IF
						CALL MPI_BCAST(E_tot_out(1:2,0:num_wf),2*(num_wf+1),MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				END IF
			ELSE
				IF (SDse_kind=='no_') THEN
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM(E_tot_out(1:2,0),E_tot,N_mc)
						IF (mpi_myrank==0) THEN
							PRINT * , 'Ho calcolato l energia --->   E_tot=', E_tot_out(1,0), '+-', E_tot_out(2,0)
							IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', E_tot_out(1,0), '+-', E_tot_out(2,0)
						END IF
						DO i_mc = 1, N_mc, 1
							dummy1(i_mc)=E_tot(i_mc)
							dummy2(i_mc)=1.d0
						END DO
						DO i = 1, num_wf, 1
							DO i_mc = 1, N_mc, 1
								dummy3(i_mc)=vwf_data(i)%w(i_mc)*vwf_data(i)%E_tot(i_mc)
								dummy4(i_mc)=vwf_data(i)%w(i_mc)
							END DO
							CALL calcola_estimatore_differenza_da_RAM(E_tot_out(1:2,i),dummy1,dummy2,dummy3,dummy4,N_mc)
						END DO
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				ELSE
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM(E_tot_out(1:2,0),E_tot,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , 'Ho calcolato l energia --->   E_tot=', E_tot_out(1,0), '+-', E_tot_out(2,0)
							IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', E_tot_out(1,0), '+-', E_tot_out(2,0)
						END IF
						DO i_mc = 1, N_mc, 1
							dummy1(i_mc)=w(i_mc)*E_tot(i_mc)
							dummy2(i_mc)=w(i_mc)
						END DO
						DO i = 1, num_wf, 1
							DO i_mc = 1, N_mc, 1
								dummy3(i_mc)=w(i_mc)*vwf_data(i)%w(i_mc)*vwf_data(i)%E_tot(i_mc)
								dummy4(i_mc)=w(i_mc)*vwf_data(i)%w(i_mc)
							END DO
							CALL calcola_estimatore_differenza_da_RAM(E_tot_out(1:2,i),dummy1,dummy2,dummy3,dummy4,N_mc)
						END DO
						ELSE
							STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
							  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				END IF
			END IF
		ELSE IF (quick_error>3) THEN
			IF (flag_disk) THEN
				IF (SDse_kind=='no_') THEN
					IF (flag_E_tot) THEN
						IF (mpi_myrank==0) THEN
							CALL calcola_estimatori_quick(quick_error,'estimatori/E_tot'//codice_simulazione,media,errore)
							IF (mpi_myrank==0) THEN
								PRINT * , 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
								IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
							END IF
							E_tot_out(1:2,0)=(/media,errore/)
							DO i = 1, num_wf, 1
								WRITE (codice_numerico, '(I4.4)'), i
								CALL calcola_estimatori_differenza_quick(quick_error,'estimatori/gradiente/E_tot'//codice_numerico, &
								  'estimatori/E_tot'//codice_simulazione,media,errore,'estimatori/gradiente/w'//codice_numerico)
								E_tot_out(1:2,i)=(/media,errore/)
							END DO
						END IF
						CALL MPI_BCAST(E_tot_out(1:2,0:num_wf),2*(num_wf+1),MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				ELSE
					IF (flag_E_tot) THEN
						IF (mpi_myrank==0) THEN
							CALL calcola_estimatori_quick(quick_error,'estimatori/E_tot'//codice_simulazione,media,errore, &
							  'estimatori/w'//codice_simulazione)
							IF (mpi_myrank==0) THEN
								PRINT * , 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
								IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', media, '+-', errore
							END IF
							E_tot_out(1:2,0)=(/media,errore/)
							DO i = 1, num_wf, 1
								WRITE (codice_numerico, '(I4.4)'), i
								CALL calcola_estimatori_differenza_quick(quick_error,'estimatori/gradiente/E_tot'//codice_numerico, &
								  'estimatori/E_tot'//codice_simulazione, media,errore,'estimatori/gradiente/w'//codice_numerico, &
								  'estimatori/w'//codice_simulazione)
								E_tot_out(1:2,i)=(/media,errore/)
							END DO
						END IF
						CALL MPI_BCAST(E_tot_out(1:2,0:num_wf),2*(num_wf+1),MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				END IF
			ELSE
				IF (SDse_kind=='no_') THEN
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM_quick(quick_error,E_tot_out(1:2,0),E_tot,N_mc)
						IF (mpi_myrank==0) THEN
							PRINT * , 'Ho calcolato l energia --->   E_tot=', E_tot_out(1,0), '+-', E_tot_out(2,0)
							IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', E_tot_out(1,0), '+-', E_tot_out(2,0)
						END IF
						DO i_mc = 1, N_mc, 1
							dummy1(i_mc)=E_tot(i_mc)
							dummy2(i_mc)=1.d0
						END DO
						DO i = 1, num_wf, 1
							DO i_mc = 1, N_mc, 1
								dummy3(i_mc)=vwf_data(i)%w(i_mc)*vwf_data(i)%E_tot(i_mc)
								dummy4(i_mc)=vwf_data(i)%w(i_mc)
							END DO
							CALL calcola_estimatore_differenza_da_RAM_quick(quick_error,E_tot_out(1:2,i),dummy1,dummy2,dummy3,dummy4,N_mc)
						END DO
					ELSE
						STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
						  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				ELSE
					IF (flag_E_tot) THEN
						CALL calcola_estimatore_da_RAM_quick(quick_error,E_tot_out(1:2,0),E_tot,N_mc,w)
						IF (mpi_myrank==0) THEN
							PRINT * , 'Ho calcolato l energia --->   E_tot=', E_tot_out(1,0), '+-', E_tot_out(2,0)
							IF (flag_output) WRITE (7, *), 'Ho calcolato l energia --->   E_tot=', E_tot_out(1,0), '+-', E_tot_out(2,0)
						END IF
						DO i_mc = 1, N_mc, 1
							dummy1(i_mc)=w(i_mc)*E_tot(i_mc)
							dummy2(i_mc)=w(i_mc)
						END DO
						DO i = 1, num_wf, 1
							DO i_mc = 1, N_mc, 1
								dummy3(i_mc)=w(i_mc)*vwf_data(i)%w(i_mc)*vwf_data(i)%E_tot(i_mc)
								dummy4(i_mc)=w(i_mc)*vwf_data(i)%w(i_mc)
							END DO
							CALL calcola_estimatore_differenza_da_RAM_quick(quick_error,E_tot_out(1:2,i),dummy1,dummy2,dummy3,dummy4,N_mc)
						END DO
						ELSE
							STOP 'Hai richiesto il calcolo dell energia ma non hai settato bene le flag &
							  [ module_VMC.f90 > restituisci_risultati ]'
					END IF
				END IF
			END IF
		ELSE
			STOP 'Il parametro quick_error non é stato scelto opportunamente'
		END IF
			
	END SUBROUTINE restituisci_energie
!-----------------------------------------------------------------------
	SUBROUTINE chiudi_VMC()
		IMPLICIT NONE
		LOGICAL :: flag_scrivi
		INTEGER :: status(MPI_STATUS_SIZE)
		INTEGER :: i, j, k
		INTEGER, ALLOCATABLE :: seed(:)
		REAL (KIND=8) :: dum_step1, dum_step2, dum_step3
		
		IF (.NOT. iniz_VMC) STOP 'Prima di chiudere devi inizializzare il VMC &
		  [ module_VMC.f90 > chiudi_VMC ]'
		
		IF (flag_disk) THEN
			CALL RANDOM_SEED(size = j)
			ALLOCATE(seed(j))
			CALL RANDOM_SEED(GET=seed)
			OPEN (UNIT=77, FILE='posizioni/seed_and_step'//codice_simulazione//'.dat', STATUS='UNKNOWN')
			IF (mpi_myrank==0) THEN
				flag_scrivi=.TRUE.
				WRITE (77, *), mpi_nprocs
				WRITE (77, *), mpi_myrank, seed, step_e, step_se, step_p
			ELSE
				flag_scrivi=.FALSE.
			END IF
					
			DO i = 1, mpi_nprocs-1, 1
				IF (mpi_myrank==i) THEN
					CALL MPI_SEND(seed,j,MPI_INTEGER,0,10+i,MPI_COMM_WORLD,mpi_ierr)
					CALL MPI_SEND(step_e,1,MPI_DOUBLE_PRECISION,0,10+i+100,MPI_COMM_WORLD,mpi_ierr)
					CALL MPI_SEND(step_se,1,MPI_DOUBLE_PRECISION,0,10+i+200,MPI_COMM_WORLD,mpi_ierr)
					CALL MPI_SEND(step_p,1,MPI_DOUBLE_PRECISION,0,10+i+300,MPI_COMM_WORLD,mpi_ierr)
				ELSE IF (mpi_myrank==0) THEN
					CALL MPI_RECV(seed,j,MPI_INTEGER,i,10+i,MPI_COMM_WORLD,status,mpi_ierr)
					CALL MPI_RECV(dum_step1,1,MPI_DOUBLE_PRECISION,i,10+i+100,MPI_COMM_WORLD,status,mpi_ierr)
					CALL MPI_RECV(dum_step2,1,MPI_DOUBLE_PRECISION,i,10+i+200,MPI_COMM_WORLD,status,mpi_ierr)
					CALL MPI_RECV(dum_step3,1,MPI_DOUBLE_PRECISION,i,10+i+300,MPI_COMM_WORLD,status,mpi_ierr)
					WRITE (77, *), i, seed, dum_step1, dum_step2, dum_step3
				END IF
			END DO
					
			DEALLOCATE(seed)
			CLOSE (77)
		END IF
		
		CALL chiudi_estimatori()
		
		IF (flag_gr) THEN
			DEALLOCATE(geeuu,geeuu_new,geeud,geeud_new)
			DEALLOCATE(gep,gep_new)
			DEALLOCATE(gppuu,gppuu_new,gppud,gppud_new)
			IF (flag_shadow) THEN
				DEALLOCATE(gseuu1,gseuu1_new,gseud1,gseud1_new,gese1,gese1_new)
				DEALLOCATE(gseuu2,gseuu2_new,gseud2,gseud2_new,gese2,gese2_new)
				DEALLOCATE(gsesp1,gsesp1_new,gsesp2,gsesp2_new)
			END IF
		END IF
				
		CALL chiudi_funzione_onda()
		CALL chiudi_walkers()
		CALL chiudi_dati_fisici()
		CALL chiudi_dati_simulazione_mc()
		
		IF (verbose_mode) THEN
			PRINT * , 'VMC: Ho chiuso il VMC'
			IF (flag_output) WRITE (7, *), 'VMC: Ho chiuso il VMC'
		END IF
				
		IF (mpi_myrank==0) THEN
			PRINT * , '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
			IF (flag_output) THEN
				WRITE (7, *), '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
				CLOSE (7)
				OPEN (UNIT=7, FILE='output.d', STATUS='OLD', POSITION='APPEND')
			END IF
		END IF
		
		iniz_VMC=.FALSE.
	END SUBROUTINE chiudi_VMC
	
END MODULE VMC
