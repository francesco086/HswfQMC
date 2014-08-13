MODULE variational_calculations
	USE funzione_onda
	USE momenta
	USE walkers
	USE dati_fisici
	USE dati_simulazione_mc
	IMPLICIT NONE
	TYPE variational_wave_function
		COMPLEX (KIND=8), ALLOCATABLE :: SDe_up(:,:), ISDe_up(:,:)
		INTEGER, ALLOCATABLE :: pvte_up(:)
		COMPLEX (KIND=8), ALLOCATABLE :: SDe_dw(:,:), ISDe_dw(:,:)
		INTEGER, ALLOCATABLE :: pvte_dw(:)
		COMPLEX (KIND=8) ::  detSDe_up, detSDe_dw
		REAL (KIND=8), ALLOCATABLE :: u_ee(:,:), u_ep(:,:)
		REAL (KIND=8) :: Uee, Uep
		COMPLEX (KIND=8), ALLOCATABLE :: SDse1_up(:,:), ISDse1_up(:,:)
		INTEGER, ALLOCATABLE :: pvtse1_up(:)
		COMPLEX (KIND=8), ALLOCATABLE :: SDse1_dw(:,:), ISDse1_dw(:,:)
		INTEGER, ALLOCATABLE :: pvtse1_dw(:)
		COMPLEX (KIND=8) ::  detSDse1_up, detSDse1_dw
		COMPLEX (KIND=8), ALLOCATABLE :: SDse2_up(:,:), ISDse2_up(:,:)
		INTEGER, ALLOCATABLE :: pvtse2_up(:)
		COMPLEX (KIND=8), ALLOCATABLE :: SDse2_dw(:,:), ISDse2_dw(:,:)
		INTEGER, ALLOCATABLE :: pvtse2_dw(:)
		COMPLEX (KIND=8) ::  detSDse2_up, detSDse2_dw
		REAL (KIND=8), ALLOCATABLE :: u_se1(:,:), u_se2(:,:)
		REAL (KIND=8) :: Use1, Use2
		REAL (KIND=8), ALLOCATABLE :: u_ese1(:), u_ese2(:)
		REAL (KIND=8) :: Uese1, Uese2
		REAL (KIND=8), ALLOCATABLE :: GDse1_up(:,:), IGDse1_up(:,:)
		INTEGER, ALLOCATABLE :: pvtgdse1_up(:)
		REAL (KIND=8), ALLOCATABLE :: GDse1_dw(:,:), IGDse1_dw(:,:)
		INTEGER, ALLOCATABLE :: pvtgdse1_dw(:)
		REAL (KIND=8) ::  detGDse1_up, detGDse1_dw
		REAL (KIND=8), ALLOCATABLE :: GDse2_up(:,:), IGDse2_up(:,:)
		INTEGER, ALLOCATABLE :: pvtgdse2_up(:)
		REAL (KIND=8), ALLOCATABLE :: GDse2_dw(:,:), IGDse2_dw(:,:)
		INTEGER, ALLOCATABLE :: pvtgdse2_dw(:)
		REAL (KIND=8) ::  detGDse2_up, detGDse2_dw
		REAL (KIND=8), ALLOCATABLE :: u_sesp1(:,:), u_sesp2(:,:)
		REAL (KIND=8) :: Usesp1, Usesp2
		REAL (KIND=8), ALLOCATABLE :: GDsp1_up(:,:), IGDsp1_up(:,:)
		INTEGER, ALLOCATABLE :: pvtgdsp1_up(:)
		REAL (KIND=8), ALLOCATABLE :: GDsp1_dw(:,:), IGDsp1_dw(:,:)
		INTEGER, ALLOCATABLE :: pvtgdsp1_dw(:)
		REAL (KIND=8) ::  detGDsp1_up, detGDsp1_dw
		REAL (KIND=8), ALLOCATABLE :: GDsp2_up(:,:), IGDsp2_up(:,:)
		INTEGER, ALLOCATABLE :: pvtgdsp2_up(:)
		REAL (KIND=8), ALLOCATABLE :: GDsp2_dw(:,:), IGDsp2_dw(:,:)
		INTEGER, ALLOCATABLE :: pvtgdsp2_dw(:)
		REAL (KIND=8) ::  detGDsp2_up, detGDsp2_dw
		REAL (KIND=8), ALLOCATABLE :: u_POT_se1(:), u_POT_se2(:)
		REAL (KIND=8), ALLOCATABLE :: E_tot(:), w(:)
	END TYPE variational_wave_function
	TYPE (variational_wave_function), ALLOCATABLE, SAVE :: vwf_data(:)
	INTEGER :: num_wf, num_par
	LOGICAL, SAVE, PUBLIC :: flag_gradiente=.FALSE., flag_derivate_var=.FALSE.
	LOGICAL, PRIVATE, SAVE :: iniz_variational_calculations=.FALSE., verbose_mode=.FALSE.
	REAL (KIND=8), ALLOCATABLE :: parametri_var(:,:)
	REAL (KIND=8), ALLOCATABLE :: O(:,:)
	
	CONTAINS
	
	SUBROUTINE inizializza_variational_calculations(num,par_pt,num_pt)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: num_pt, num                  !num_pt=numero punti da considerare per il gradiente (generalmente = num), num=numero di parametri variazionali
		REAL (KIND=8) :: par_pt(1:num,0:num_pt)
		INTEGER :: i
		
		num_par=num
		IF (what_to_do=='congrad') THEN
			flag_gradiente=.TRUE.
			num_wf=num_pt
			ALLOCATE(parametri_var(1:num_par,0:num_wf))
			parametri_var=par_pt
			ALLOCATE(vwf_data(1:num_wf))
			DO i = 1, num_wf, 1
				IF ((SDe_kind/='no_').AND.(SDe_kind/='pw_').AND.(SDe_kind/='lda')) THEN
					ALLOCATE(vwf_data(i)%SDe_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%SDe_up=0.d0
					ALLOCATE(vwf_data(i)%ISDe_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%ISDe_up=0.d0
					ALLOCATE(vwf_data(i)%pvte_up(1:H_N_part))
					vwf_data(i)%pvte_up=0.d0
					ALLOCATE(vwf_data(i)%SDe_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%SDe_dw=0.d0
					ALLOCATE(vwf_data(i)%ISDe_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%ISDe_dw=0.d0
					ALLOCATE(vwf_data(i)%pvte_dw(1:H_N_part))
					vwf_data(i)%pvte_dw=0.d0
				END IF
				IF (Jee_kind/='no_') THEN
					ALLOCATE(vwf_data(i)%u_ee(1:N_part,1:N_part))
					vwf_data(i)%u_ee=0.d0
				END IF
				IF (Jep_kind/='no_') THEN
					ALLOCATE(vwf_data(i)%u_ep(1:N_part,1:N_part))
					vwf_data(i)%u_ep=0.d0
				END IF
				IF ((SDse_kind/='no_').AND.(SDse_kind/='pw_').AND.(SDse_kind/='lda')) THEN
					ALLOCATE(vwf_data(i)%SDse1_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%SDse1_up=0.d0
					ALLOCATE(vwf_data(i)%ISDse1_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%ISDse1_up=0.d0
					ALLOCATE(vwf_data(i)%pvtse1_up(1:H_N_part))
					vwf_data(i)%pvtse1_up=0.d0
					ALLOCATE(vwf_data(i)%SDse1_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%SDse1_dw=0.d0
					ALLOCATE(vwf_data(i)%ISDse1_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%ISDse1_dw=0.d0
					ALLOCATE(vwf_data(i)%pvtse1_dw(1:H_N_part))
					vwf_data(i)%pvtse1_dw=0.d0
					ALLOCATE(vwf_data(i)%SDse2_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%SDse2_up=0.d0
					ALLOCATE(vwf_data(i)%ISDse2_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%ISDse2_up=0.d0
					ALLOCATE(vwf_data(i)%pvtse2_up(1:H_N_part))
					vwf_data(i)%pvtse2_up=0.d0
					ALLOCATE(vwf_data(i)%SDse2_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%SDse2_dw=0.d0
					ALLOCATE(vwf_data(i)%ISDse2_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%ISDse2_dw=0.d0
					ALLOCATE(vwf_data(i)%pvtse2_dw(1:H_N_part))
					vwf_data(i)%pvtse2_dw=0.d0
				END IF
				IF (Jse_kind/='no_') THEN
					IF ((Jse_kind/='no_').AND.(Jse_kind/='pot')) THEN
						ALLOCATE(vwf_data(i)%u_POT_se1(1:H_N_part))
						vwf_data(i)%u_POT_se1=0.d0
						ALLOCATE(vwf_data(i)%u_POT_se2(1:H_N_part))
						vwf_data(i)%u_POT_se2=0.d0
					END IF
				END IF
				IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
					ALLOCATE(vwf_data(i)%u_ese1(1:N_part))
					vwf_data(i)%u_ese1=0.d0
					ALLOCATE(vwf_data(i)%u_ese2(1:N_part))
					vwf_data(i)%u_ese2=0.d0
				ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
					ALLOCATE(vwf_data(i)%GDse1_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%GDse1_up=0.d0
					ALLOCATE(vwf_data(i)%IGDse1_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%IGDse1_up=0.d0
					ALLOCATE(vwf_data(i)%pvtgdse1_up(1:H_N_part))
					vwf_data(i)%pvtgdse1_up=0.d0
					ALLOCATE(vwf_data(i)%GDse1_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%GDse1_dw=0.d0
					ALLOCATE(vwf_data(i)%IGDse1_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%IGDse1_dw=0.d0
					ALLOCATE(vwf_data(i)%pvtgdse1_dw(1:H_N_part))
					vwf_data(i)%pvtgdse1_dw=0.d0
					ALLOCATE(vwf_data(i)%GDse2_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%GDse2_up=0.d0
					ALLOCATE(vwf_data(i)%IGDse2_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%IGDse2_up=0.d0
					ALLOCATE(vwf_data(i)%pvtgdse2_up(1:H_N_part))
					vwf_data(i)%pvtgdse2_up=0.d0
					ALLOCATE(vwf_data(i)%GDse2_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%GDse2_dw=0.d0
					ALLOCATE(vwf_data(i)%IGDse2_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%IGDse2_dw=0.d0
					ALLOCATE(vwf_data(i)%pvtgdse2_dw(1:H_N_part))
					vwf_data(i)%pvtgdse2_dw=0.d0
				END IF
				IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
					ALLOCATE(vwf_data(i)%u_sesp1(1:N_part,1:N_part))
					vwf_data(i)%u_sesp1=0.d0
					ALLOCATE(vwf_data(i)%u_sesp2(1:N_part,1:N_part))
					vwf_data(i)%u_sesp2=0.d0
				ELSE IF (Jsesp_kind=='gsd') THEN
					ALLOCATE(vwf_data(i)%GDsp1_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%GDsp1_up=0.d0
					ALLOCATE(vwf_data(i)%IGDsp1_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%IGDsp1_up=0.d0
					ALLOCATE(vwf_data(i)%pvtgdsp1_up(1:H_N_part))
					vwf_data(i)%pvtgdsp1_up=0.d0
					ALLOCATE(vwf_data(i)%GDsp1_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%GDsp1_dw=0.d0
					ALLOCATE(vwf_data(i)%IGDsp1_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%IGDsp1_dw=0.d0
					ALLOCATE(vwf_data(i)%pvtgdsp1_dw(1:H_N_part))
					vwf_data(i)%pvtgdsp1_dw=0.d0
					ALLOCATE(vwf_data(i)%GDsp2_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%GDsp2_up=0.d0
					ALLOCATE(vwf_data(i)%IGDsp2_up(1:H_N_part,1:H_N_part))
					vwf_data(i)%IGDsp2_up=0.d0
					ALLOCATE(vwf_data(i)%pvtgdsp2_up(1:H_N_part))
					vwf_data(i)%pvtgdsp2_up=0.d0
					ALLOCATE(vwf_data(i)%GDsp2_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%GDsp2_dw=0.d0
					ALLOCATE(vwf_data(i)%IGDsp2_dw(1:H_N_part,1:H_N_part))
					vwf_data(i)%IGDsp2_dw=0.d0
					ALLOCATE(vwf_data(i)%pvtgdsp2_dw(1:H_N_part))
					vwf_data(i)%pvtgdsp2_dw=0.d0
				END IF
				
				ALLOCATE(vwf_data(i)%E_tot(1:N_mc),vwf_data(i)%w(1:N_mc))
			END DO
		ELSE IF ((what_to_do=='stocrec').OR.(what_to_do=='stoc_ns').OR.(what_to_do=='stoc_av').OR.(what_to_do=='pure_sr')) THEN
			flag_derivate_var=.TRUE.
			ALLOCATE(O(1:num_par,1:N_mc))
		END IF
		
		iniz_variational_calculations=.TRUE.
		
	END SUBROUTINE inizializza_variational_calculations
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatori_per_gradiente(i_mc)
		IMPLICIT NONE
		INTEGER (KIND=8), INTENT(IN) :: i_mc
		INTEGER :: i
		IF (.NOT.iniz_variational_calculations) STOP 'Prima di calcolare gli estimatori devi inizializzare &
		  [ module_variational_calculations.f90 > calcola_estimatori_per_gradiente ]'
		
		DO i = 1, num_wf, 1
			CALL calcola_termini_funzione_onda_var(i,i_mc)
			CALL calcola_estimatori_var(i,i_mc)
		END DO
		
	END SUBROUTINE calcola_estimatori_per_gradiente
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatori_per_derivate_variazionali(i_mc)
		IMPLICIT NONE
		INTEGER (KIND=8), INTENT(IN) :: i_mc
		IF (.NOT.iniz_variational_calculations) STOP 'Prima di calcolare gli estimatori devi inizializzare &
		  [ module_variational_calculations.f90 > calcola_estimatori_per_gradiente ]'
		
		CALL calcola_termini_derivate_variazionali(i_mc)
		
	END SUBROUTINE calcola_estimatori_per_derivate_variazionali
!-----------------------------------------------------------------------
	SUBROUTINE conferma_estimatori_per_gradiente(i_mc)
		IMPLICIT NONE
		INTEGER (KIND=8), INTENT(IN) :: i_mc
		INTEGER :: i
		IF (.NOT.iniz_variational_calculations) STOP 'Prima di calcolare gli estimatori devi inizializzare &
		  [ module_variational_calculations.f90 > conferma_estimatori_per_gradiente ]'
		
		DO i = 1, num_wf, 1
			vwf_data(i)%w(i_mc)=vwf_data(i)%w(i_mc-1)
			vwf_data(i)%E_tot(i_mc)=vwf_data(i)%E_tot(i_mc-1)
		END DO
		
	END SUBROUTINE conferma_estimatori_per_gradiente
!-----------------------------------------------------------------------
	SUBROUTINE conferma_estimatori_per_derivate_variazionali(i_mc)
		IMPLICIT NONE
		INTEGER (KIND=8), INTENT(IN) :: i_mc
		IF (.NOT.iniz_variational_calculations) STOP 'Prima di calcolare gli estimatori devi inizializzare &
		  [ module_variational_calculations.f90 > conferma_estimatori_per_gradiente ]'
		
		O(1:num_par,i_mc)=O(1:num_par,i_mc-1)
		
	END SUBROUTINE conferma_estimatori_per_derivate_variazionali
!-----------------------------------------------------------------------
	SUBROUTINE trascrivi_dati_per_gradiente()
		USE estimatori
		IMPLICIT NONE
		LOGICAL :: flag_scrivi, flag_continua_rank
		CHARACTER(LEN=4) :: codice_numerico
		INTEGER (KIND=8) :: i_mc
		INTEGER :: i, j
		INTEGER :: status(MPI_STATUS_SIZE)
		
		IF (flag_disk) THEN 
		IF (mpi_myrank==0) THEN
			flag_scrivi=.TRUE.
			flag_continua_rank=flag_continua
		ELSE
			flag_scrivi=.FALSE.
			flag_continua_rank=.TRUE.
		END IF
		
		DO j = 0, mpi_nprocs-1, 1
			IF (mpi_myrank==j .AND. flag_scrivi) THEN
				IF (SDse_kind=='no_') THEN
					DO i = 1, num_wf, 1
						DO i_mc = 1, N_mc, 1
							vwf_data(i)%E_tot(i_mc)=vwf_data(i)%w(i_mc)*vwf_data(i)%E_tot(i_mc)
						END DO
					END DO
				ELSE
					DO i = 1, num_wf, 1
						DO i_mc = 1, N_mc, 1
							vwf_data(i)%E_tot(i_mc)=vwf_data(i)%w(i_mc)*vwf_data(i)%E_tot(i_mc)*w(i_mc)
							vwf_data(i)%w(i_mc)=vwf_data(i)%w(i_mc)*w(i_mc)
						END DO
					END DO
				END IF
				DO i = 1, num_wf, 1
					WRITE (codice_numerico, '(I4.4)'), i
					CALL salva_vettore_dati_bin(vwf_data(i)%E_tot,N_mc,N_AV,flag_continua_rank,'estimatori/gradiente/E_tot'//codice_numerico)
					IF (mpi_myrank==0) THEN
						CALL salva_vettore_dati_dat_ridotto(vwf_data(i)%E_tot,N_mc,10000_8,'estimatori/gradiente/E_tot-sample'//codice_numerico)
					END IF
					CALL salva_vettore_dati_bin(vwf_data(i)%w,N_mc,N_AV,flag_continua_rank,'estimatori/gradiente/w'//codice_numerico)
					IF (mpi_myrank==0) THEN
						CALL salva_vettore_dati_dat_ridotto(vwf_data(i)%w,N_mc,10000_8,'estimatori/gradiente/w-sample'//codice_numerico)
					END IF
				END DO
				IF (mpi_myrank<mpi_nprocs-1) THEN
					CALL MPI_SEND(flag_scrivi,1,MPI_LOGICAL,mpi_myrank+1,10+j,MPI_COMM_WORLD,mpi_ierr)
					IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_SEND &
					  [ module_variational_calculations.f90 > trascrivi_dati_per_gradiente ]'
				END IF
			END IF
			IF (mpi_myrank==j+1 .AND. j<mpi_nprocs) THEN
				CALL MPI_RECV(flag_scrivi,1,MPI_LOGICAL,mpi_myrank-1,10+j,MPI_COMM_WORLD,status,mpi_ierr)
				IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_RECV &
				  [ module_variational_calculations.f90 > trascrivi_dati_per_gradiente ]'
			END IF
		END DO
		
		IF (flag_mpi) CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
		END IF
		
	END SUBROUTINE trascrivi_dati_per_gradiente
!-----------------------------------------------------------------------
	SUBROUTINE trascrivi_dati_per_derivate_variazionali()
		USE estimatori
		IMPLICIT NONE
		LOGICAL :: flag_scrivi, flag_continua_rank
		CHARACTER(LEN=4) :: codice_numerico, codice_numerico2
		INTEGER (KIND=8) :: i_mc
		INTEGER :: i, j, i_cpu
		INTEGER :: status(MPI_STATUS_SIZE)
		REAL (KIND=8) :: vec_save(1:N_mc)
		
		IF (mpi_myrank==0) THEN
			flag_scrivi=.TRUE.
			flag_continua_rank=flag_continua
		ELSE
			flag_scrivi=.FALSE.
			flag_continua_rank=.TRUE.
		END IF
		
		DO i_cpu = 0, mpi_nprocs-1, 1
			IF (mpi_myrank==i_cpu .AND. flag_scrivi) THEN
				DO j = 1, num_par, 1
					DO i = j, num_par, 1
						IF (SDse_kind=='no_') THEN
							DO i_mc = 1, N_mc, 1
								vec_save(i_mc)=O(i,i_mc)*O(j,i_mc)
							END DO
						ELSE
							DO i_mc = 1, N_mc, 1
								vec_save(i_mc)=O(i,i_mc)*O(j,i_mc)*w(i_mc)
							END DO
						END IF
						WRITE (codice_numerico, '(I4.4)'), i
						WRITE (codice_numerico2, '(I4.4)'), j
						CALL salva_vettore_dati_bin(vec_save(1:N_mc),N_mc,N_AV,flag_continua_rank, &
						  'estimatori/gradiente/O_'//codice_numerico2//'-'//codice_numerico)
					END DO
				END DO
				DO i = 1, num_par, 1
					IF (SDse_kind=='no_') THEN
						DO i_mc = 1, N_mc, 1
							vec_save(i_mc)=O(i,i_mc)*E_tot(i_mc)
						END DO
					ELSE
						DO i_mc = 1, N_mc, 1
							vec_save(i_mc)=O(i,i_mc)*E_tot(i_mc)*w(i_mc)
						END DO
					END IF
					WRITE (codice_numerico, '(I4.4)'), i
					CALL salva_vettore_dati_bin(vec_save(1:N_mc),N_mc,N_AV,flag_continua_rank,'estimatori/gradiente/EO_'//codice_numerico)
					IF (SDse_kind=='no_') THEN
						WRITE (codice_numerico, '(I4.4)'), i
						CALL salva_vettore_dati_bin(O(i,1:N_mc),N_mc,N_AV,flag_continua_rank,'estimatori/gradiente/O_'//codice_numerico)
					ELSE
						DO i_mc = 1, N_mc, 1
							vec_save(i_mc)=O(i,i_mc)*w(i_mc)
						END DO
						WRITE (codice_numerico, '(I4.4)'), i
						CALL salva_vettore_dati_bin(vec_save(1:N_mc),N_mc,N_AV,flag_continua_rank,'estimatori/gradiente/O_'//codice_numerico)
					END IF
				END DO
				IF (mpi_myrank<mpi_nprocs-1) THEN
					CALL MPI_SEND(flag_scrivi,1,MPI_LOGICAL,mpi_myrank+1,10+i_cpu,MPI_COMM_WORLD,mpi_ierr)
					IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_SEND &
					  [ module_variational_calculations.f90 > trascrivi_dati_per_derivate_variazionali ]'
				END IF
			END IF
			IF (mpi_myrank==i_cpu+1 .AND. i_cpu<mpi_nprocs) THEN
				CALL MPI_RECV(flag_scrivi,1,MPI_LOGICAL,mpi_myrank-1,10+i_cpu,MPI_COMM_WORLD,status,mpi_ierr)
				IF (mpi_ierr/=MPI_SUCCESS) STOP 'Errore in MPI_RECV &
				  [ module_variational_calculations.f90 > trascrivi_dati_per_derivate_variazionali ]'
			END IF
		END DO
		
		IF (flag_mpi) CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
	END SUBROUTINE trascrivi_dati_per_derivate_variazionali
!-----------------------------------------------------------------------
	SUBROUTINE calcola_termini_funzione_onda_var(i_wf,i_mc)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: i_wf
		INTEGER (KIND=8), INTENT(IN) :: i_mc
		INTEGER :: info
		COMPLEX (KIND=8) :: work(1:3*H_N_part)
		
		IF (.NOT.iniz_variational_calculations) STOP 'Prima di calcolare i termini delle funzioni d onda devi inizializzare &
		  [ module_variational_calculations.f90 > calcola_termini_funzioni_onda_var ]'
		
		CALL setta_parametri(parametri_var(1:num_par,i_wf),num_par)
		
		SELECT CASE (SDe_kind)
		CASE ('pw_') 
			vwf_data(i_wf)%detSDe_up=1.d0
			vwf_data(i_wf)%detSDe_dw=1.d0
		CASE ('lda') 
			vwf_data(i_wf)%detSDe_up=1.d0
			vwf_data(i_wf)%detSDe_dw=1.d0
		CASE ('no_') 
			vwf_data(i_wf)%detSDe_up=1.d0
			vwf_data(i_wf)%detSDe_dw=1.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di SDe_kind accettabile &
			  [ module_variational_calculations.f90 > calcola_termini_funzione_onda_var ]'
		END SELECT
        
		SELECT CASE (Jee_kind)
		CASE ('yuk')
			CALL valuta_Uee_YUK(-1,rij_ee_old,N_part,vwf_data(i_wf)%u_ee,vwf_data(i_wf)%Uee)
		CASE ('yup') 
			CALL valuta_Uee_YUK(-1,rijpc_ee_old,N_part,vwf_data(i_wf)%u_ee,vwf_data(i_wf)%Uee)
		CASE ('no_') 
			vwf_data(i_wf)%Uee=0.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Jee_kind accettabile &
			  [ module_variational_calculations.f90 > calcola_termini_funzione_onda_var ]'
		END SELECT
        
		SELECT CASE (Jep_kind)
		CASE ('yuk')
			CALL valuta_Uep_YUK(-1,0,rij_ep_old,N_part,vwf_data(i_wf)%u_ep,vwf_data(i_wf)%Uep)
		CASE ('yup') 
			CALL valuta_Uep_YUK(-1,0,rijpc_ep_old,N_part,vwf_data(i_wf)%u_ep,vwf_data(i_wf)%Uep)
		CASE ('no_') 
			vwf_data(i_wf)%Uep=0.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Jep_kind accettabile &
			  [ module_variational_calculations.f90 > calcola_termini_funzione_onda_var ]'
		END SELECT
        
		SELECT CASE (SDse_kind)
		CASE ('pw_') 
			vwf_data(i_wf)%detSDse1_up=1.d0
			vwf_data(i_wf)%detSDse1_dw=1.d0
			vwf_data(i_wf)%detSDse2_up=1.d0
			vwf_data(i_wf)%detSDse2_dw=1.d0
		CASE ('lda')
			vwf_data(i_wf)%detSDse1_up=1.d0
			vwf_data(i_wf)%detSDse1_dw=1.d0
			vwf_data(i_wf)%detSDse2_up=1.d0
			vwf_data(i_wf)%detSDse2_dw=1.d0
		CASE ('no_') 
			vwf_data(i_wf)%detSDse1_up=1.d0
			vwf_data(i_wf)%detSDse1_dw=1.d0
			vwf_data(i_wf)%detSDse2_up=1.d0
			vwf_data(i_wf)%detSDse2_dw=1.d0
		CASE ('gss') 
			vwf_data(i_wf)%detSDse1_up=1.d0
			vwf_data(i_wf)%detSDse1_dw=1.d0
			vwf_data(i_wf)%detSDse2_up=1.d0
			vwf_data(i_wf)%detSDse2_dw=1.d0
		CASE ('atm') 
			vwf_data(i_wf)%detSDse1_up=1.d0
			vwf_data(i_wf)%detSDse1_dw=1.d0
			vwf_data(i_wf)%detSDse2_up=1.d0
			vwf_data(i_wf)%detSDse2_dw=1.d0
		CASE ('atp') 
			vwf_data(i_wf)%detSDse1_up=1.d0
			vwf_data(i_wf)%detSDse1_dw=1.d0
			vwf_data(i_wf)%detSDse2_up=1.d0
			vwf_data(i_wf)%detSDse2_dw=1.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di SDse_kind accettabile &
			  [ module_variational_calculations.f90 > calcola_termini_funzione_onda_var ]'
		END SELECT
        
		SELECT CASE (Jse_kind)
		CASE ('pot')
			CALL valuta_Use_POT(-1,dist_mol_ss1_old,vwf_data(i_wf)%u_POT_se1,vwf_data(i_wf)%Use1)
			CALL valuta_Use_POT(-1,dist_mol_ss2_old,vwf_data(i_wf)%u_POT_se2,vwf_data(i_wf)%Use2)
		CASE ('yuk')
			CALL valuta_Usese_YUK(-1,rij_se1_old,N_part,vwf_data(i_wf)%u_se1,vwf_data(i_wf)%Use1)
			CALL valuta_Usese_YUK(-1,rij_se2_old,N_part,vwf_data(i_wf)%u_se2,vwf_data(i_wf)%Use2)
		CASE ('yup')
			CALL valuta_Usese_YUK(-1,rijpc_se1_old,N_part,vwf_data(i_wf)%u_se1,vwf_data(i_wf)%Use1)
			CALL valuta_Usese_YUK(-1,rijpc_se2_old,N_part,vwf_data(i_wf)%u_se2,vwf_data(i_wf)%Use2)
		CASE ('no_')
			vwf_data(i_wf)%Use1=0.d0
			vwf_data(i_wf)%Use2=0.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Jse_kind accettabile &
			  [ module_variational_calculations.f90 > calcola_termini_funzione_onda_var ]'
		END SELECT
        
		SELECT CASE (Kse_kind)
		CASE ('gss')
			CALL valuta_KERNse(-1,rij_ese1_old,N_part,vwf_data(i_wf)%u_ese1,vwf_data(i_wf)%Uese1)
			CALL valuta_KERNse(-1,rij_ese2_old,N_part,vwf_data(i_wf)%u_ese2,vwf_data(i_wf)%Uese2)
			vwf_data(i_wf)%detGDse1_up=1.d0
			vwf_data(i_wf)%detGDse1_dw=1.d0
			vwf_data(i_wf)%detGDse2_up=1.d0
			vwf_data(i_wf)%detGDse2_dw=1.d0
		CASE ('gsc')
			CALL valuta_KERNse(-1,rij_ese1_old,N_part,vwf_data(i_wf)%u_ese1,vwf_data(i_wf)%Uese1)
			CALL valuta_KERNse(-1,rij_ese2_old,N_part,vwf_data(i_wf)%u_ese2,vwf_data(i_wf)%Uese2)
			vwf_data(i_wf)%detGDse1_up=1.d0
			vwf_data(i_wf)%detGDse1_dw=1.d0
			vwf_data(i_wf)%detGDse2_up=1.d0
			vwf_data(i_wf)%detGDse2_dw=1.d0
		CASE ('gsp')
			CALL valuta_KERNse(-1,rijpc_ese1_old,N_part,vwf_data(i_wf)%u_ese1,vwf_data(i_wf)%Uese1)
			CALL valuta_KERNse(-1,rijpc_ese2_old,N_part,vwf_data(i_wf)%u_ese2,vwf_data(i_wf)%Uese2)
			vwf_data(i_wf)%detGDse1_up=1.d0
			vwf_data(i_wf)%detGDse1_dw=1.d0
			vwf_data(i_wf)%detGDse2_up=1.d0
			vwf_data(i_wf)%detGDse2_dw=1.d0
		CASE ('gsd')
			CALL valuta_GDse(-1,0,rij_ese1_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  vwf_data(i_wf)%GDse1_up,vwf_data(i_wf)%detGDse1_up,vwf_data(i_wf)%IGDse1_up, &
			  vwf_data(i_wf)%pvtgdse1_up,vwf_data(i_wf)%GDse1_up,vwf_data(i_wf)%detGDse1_up)
			CALL valuta_GDse(-1,0,rij_ese1_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  vwf_data(i_wf)%GDse1_dw,vwf_data(i_wf)%detGDse1_dw,vwf_data(i_wf)%IGDse1_dw, &
			  vwf_data(i_wf)%pvtgdse1_dw,vwf_data(i_wf)%GDse1_dw,vwf_data(i_wf)%detGDse1_dw)
			CALL valuta_GDse(-1,0,rij_ese2_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  vwf_data(i_wf)%GDse2_up,vwf_data(i_wf)%detGDse2_up,vwf_data(i_wf)%IGDse2_up, &
			  vwf_data(i_wf)%pvtgdse2_up,vwf_data(i_wf)%GDse2_up,vwf_data(i_wf)%detGDse2_up)
			CALL valuta_GDse(-1,0,rij_ese2_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  vwf_data(i_wf)%GDse2_dw,vwf_data(i_wf)%detGDse2_dw,vwf_data(i_wf)%IGDse2_dw, &
			  vwf_data(i_wf)%pvtgdse2_dw,vwf_data(i_wf)%GDse2_dw,vwf_data(i_wf)%detGDse2_dw)
			vwf_data(i_wf)%Uese1=0.d0
			vwf_data(i_wf)%Uese2=0.d0
		 CASE ('gdc')
		 	CALL valuta_GDse(-1,0,rij_ese1_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
		 	  vwf_data(i_wf)%GDse1_up,vwf_data(i_wf)%detGDse1_up,vwf_data(i_wf)%IGDse1_up, &
		 	  vwf_data(i_wf)%pvtgdse1_up,vwf_data(i_wf)%GDse1_up,vwf_data(i_wf)%detGDse1_up)
		 	CALL valuta_GDse(-1,0,rij_ese1_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
		 	  vwf_data(i_wf)%GDse1_dw,vwf_data(i_wf)%detGDse1_dw,vwf_data(i_wf)%IGDse1_dw, &
		 	  vwf_data(i_wf)%pvtgdse1_dw,vwf_data(i_wf)%GDse1_dw,vwf_data(i_wf)%detGDse1_dw)
		 	CALL valuta_GDse(-1,0,rij_ese2_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
		 	  vwf_data(i_wf)%GDse2_up,vwf_data(i_wf)%detGDse2_up,vwf_data(i_wf)%IGDse2_up, &
		 	  vwf_data(i_wf)%pvtgdse2_up,vwf_data(i_wf)%GDse2_up,vwf_data(i_wf)%detGDse2_up)
		 	CALL valuta_GDse(-1,0,rij_ese2_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
		 	  vwf_data(i_wf)%GDse2_dw,vwf_data(i_wf)%detGDse2_dw,vwf_data(i_wf)%IGDse2_dw, &
		 	  vwf_data(i_wf)%pvtgdse2_dw,vwf_data(i_wf)%GDse2_dw,vwf_data(i_wf)%detGDse2_dw)
		 	vwf_data(i_wf)%Uese1=0.d0
		 	vwf_data(i_wf)%Uese2=0.d0
		CASE ('gdp')
		 	CALL valuta_GDse(-1,0,rij_ese1_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
		 	  vwf_data(i_wf)%GDse1_up,vwf_data(i_wf)%detGDse1_up,vwf_data(i_wf)%IGDse1_up, &
		 	  vwf_data(i_wf)%pvtgdse1_up,vwf_data(i_wf)%GDse1_up,vwf_data(i_wf)%detGDse1_up)
		 	CALL valuta_GDse(-1,0,rij_ese1_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
		 	  vwf_data(i_wf)%GDse1_dw,vwf_data(i_wf)%detGDse1_dw,vwf_data(i_wf)%IGDse1_dw, &
		 	  vwf_data(i_wf)%pvtgdse1_dw,vwf_data(i_wf)%GDse1_dw,vwf_data(i_wf)%detGDse1_dw)
		 	CALL valuta_GDse(-1,0,rij_ese2_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
		 	  vwf_data(i_wf)%GDse2_up,vwf_data(i_wf)%detGDse2_up,vwf_data(i_wf)%IGDse2_up, &
		 	  vwf_data(i_wf)%pvtgdse2_up,vwf_data(i_wf)%GDse2_up,vwf_data(i_wf)%detGDse2_up)
		 	CALL valuta_GDse(-1,0,rij_ese2_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
		 	  vwf_data(i_wf)%GDse2_dw,vwf_data(i_wf)%detGDse2_dw,vwf_data(i_wf)%IGDse2_dw, &
		 	  vwf_data(i_wf)%pvtgdse2_dw,vwf_data(i_wf)%GDse2_dw,vwf_data(i_wf)%detGDse2_dw)
		 	vwf_data(i_wf)%Uese1=0.d0
		 	vwf_data(i_wf)%Uese2=0.d0
		CASE ('no_')
			vwf_data(i_wf)%Uese1=0.d0
			vwf_data(i_wf)%Uese2=0.d0
			vwf_data(i_wf)%detGDse1_up=1.d0
			vwf_data(i_wf)%detGDse1_dw=1.d0
			vwf_data(i_wf)%detGDse2_up=1.d0
			vwf_data(i_wf)%detGDse2_dw=1.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Kse_kind accettabile &
			  [ module_variational_calculations.f90 > calcola_termini_funzione_onda_var ]'
		END SELECT
        
		SELECT CASE (Jsesp_kind)
		CASE ('pot')
			CALL valuta_Usesp_POT(-1,0,rij_sesp1_old,N_part,vwf_data(i_wf)%u_sesp1,vwf_data(i_wf)%Usesp1)
			CALL valuta_Usesp_POT(-1,0,rij_sesp2_old,N_part,vwf_data(i_wf)%u_sesp2,vwf_data(i_wf)%Usesp2)
			vwf_data(i_wf)%detGDsp1_up=1.d0
			vwf_data(i_wf)%detGDsp1_dw=1.d0
			vwf_data(i_wf)%detGDsp2_up=1.d0
			vwf_data(i_wf)%detGDsp2_dw=1.d0
		CASE ('yuk')
			CALL valuta_Usesp_YUK(-1,0,rij_sesp1_old,N_part,vwf_data(i_wf)%u_sesp1,vwf_data(i_wf)%Usesp1)
			CALL valuta_Usesp_YUK(-1,0,rij_sesp2_old,N_part,vwf_data(i_wf)%u_sesp2,vwf_data(i_wf)%Usesp2)
			vwf_data(i_wf)%detGDsp1_up=1.d0
			vwf_data(i_wf)%detGDsp1_dw=1.d0
			vwf_data(i_wf)%detGDsp2_up=1.d0
			vwf_data(i_wf)%detGDsp2_dw=1.d0
		CASE ('yup')
			CALL valuta_Usesp_YUK(-1,0,rijpc_sesp1_old,N_part,vwf_data(i_wf)%u_sesp1,vwf_data(i_wf)%Usesp1)
			CALL valuta_Usesp_YUK(-1,0,rijpc_sesp2_old,N_part,vwf_data(i_wf)%u_sesp2,vwf_data(i_wf)%Usesp2)
			vwf_data(i_wf)%detGDsp1_up=1.d0
			vwf_data(i_wf)%detGDsp1_dw=1.d0
			vwf_data(i_wf)%detGDsp2_up=1.d0
			vwf_data(i_wf)%detGDsp2_dw=1.d0
		CASE ('gss')
			CALL valuta_Usesp_GSS(-1,0,rij_sesp1_old,N_part,vwf_data(i_wf)%u_sesp1,vwf_data(i_wf)%Usesp1)
			CALL valuta_Usesp_GSS(-1,0,rij_sesp2_old,N_part,vwf_data(i_wf)%u_sesp2,vwf_data(i_wf)%Usesp2)
			vwf_data(i_wf)%detGDsp1_up=1.d0
			vwf_data(i_wf)%detGDsp1_dw=1.d0
			vwf_data(i_wf)%detGDsp2_up=1.d0
			vwf_data(i_wf)%detGDsp2_dw=1.d0
		CASE ('gsd')
			CALL valuta_GDsp(-1,0,rij_sesp1_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  vwf_data(i_wf)%GDsp1_up,vwf_data(i_wf)%detGDsp1_up,vwf_data(i_wf)%IGDsp1_up, &
			  vwf_data(i_wf)%pvtgdsp1_up,vwf_data(i_wf)%GDsp1_up,vwf_data(i_wf)%detGDsp1_up)
			CALL valuta_GDsp(-1,0,rij_sesp1_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  vwf_data(i_wf)%GDsp1_dw,vwf_data(i_wf)%detGDsp1_dw,vwf_data(i_wf)%IGDsp1_dw, &
			  vwf_data(i_wf)%pvtgdsp1_dw,vwf_data(i_wf)%GDsp1_dw,vwf_data(i_wf)%detGDsp1_dw)
			CALL valuta_GDsp(-1,0,rij_sesp2_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  vwf_data(i_wf)%GDsp2_up,vwf_data(i_wf)%detGDsp2_up,vwf_data(i_wf)%IGDsp2_up, &
			  vwf_data(i_wf)%pvtgdsp2_up,vwf_data(i_wf)%GDsp2_up,vwf_data(i_wf)%detGDsp2_up)
			CALL valuta_GDsp(-1,0,rij_sesp2_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  vwf_data(i_wf)%GDsp2_dw,vwf_data(i_wf)%detGDsp2_dw,vwf_data(i_wf)%IGDsp2_dw, &
			  vwf_data(i_wf)%pvtgdsp2_dw,vwf_data(i_wf)%GDsp2_dw,vwf_data(i_wf)%detGDsp2_dw)
			vwf_data(i_wf)%Usesp1=0.d0
			vwf_data(i_wf)%Usesp2=0.d0
		CASE ('no_') 
			vwf_data(i_wf)%Usesp1=0.d0
			vwf_data(i_wf)%Usesp2=0.d0
			vwf_data(i_wf)%detGDsp1_up=1.d0
			vwf_data(i_wf)%detGDsp1_dw=1.d0
			vwf_data(i_wf)%detGDsp2_up=1.d0
			vwf_data(i_wf)%detGDsp2_dw=1.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Jsesp_kind accettabile &
			  [ module_variational_calculations.f90 > calcola_termini_funzione_onda_var ]'
		END SELECT
        
		IF (verbose_mode) PRINT *, 'calcola_accettazione: Ho inizializzato la funzione d onda'
		
		CALL setta_parametri(parametri_var(1:num_par,0),num_par)
				
	END SUBROUTINE calcola_termini_funzione_onda_var
!-----------------------------------------------------------------------

	SUBROUTINE calcola_termini_derivate_variazionali(i_mc)
		USE estimatori
		IMPLICIT NONE
		INTEGER (KIND=8), INTENT(IN) :: i_mc
		INTEGER :: cont, i
		
		cont=1
		
		SELECT CASE (Jee_kind)
		CASE ('yuk')
			IF (opt_A_Jee) THEN
				IF (opt_F_Jee) THEN
					IF (split_Aee) THEN
						IF (split_Fee) THEN
							CALL derivata_Jee_YUK(O(cont:cont+3,i_mc))
							cont=cont+4
						ELSE
							CALL derivata_Jee_YUK(O(cont:cont+2,i_mc))
							cont=cont+3
						END IF
					ELSE
						IF (split_Fee) THEN
							CALL derivata_Jee_YUK(O(cont:cont+2,i_mc))
							cont=cont+3
						ELSE
							CALL derivata_Jee_YUK(O(cont:cont+1,i_mc))
							cont=cont+2
						END IF
					END IF
				ELSE
					IF (split_Aee) THEN
						CALL derivata_Jee_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					ELSE
						CALL derivata_Jee_YUK(O(cont:cont,i_mc))
						cont=cont+1
					END IF
				END IF
			ELSE
				IF (opt_F_Jee) THEN
					IF (split_Fee) THEN
						CALL derivata_Jee_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					ELSE
						CALL derivata_Jee_YUK(O(cont:cont,i_mc))
						cont=cont+1
					END IF
				END IF
			END IF
		CASE ('yup')
			IF (opt_A_Jee) THEN
				IF (opt_F_Jee) THEN
					IF (split_Aee) THEN
						IF (split_Fee) THEN
							CALL derivata_Jee_YUK(O(cont:cont+3,i_mc))
							cont=cont+4
						ELSE
							CALL derivata_Jee_YUK(O(cont:cont+2,i_mc))
							cont=cont+3
						END IF
					ELSE
						IF (split_Fee) THEN
							CALL derivata_Jee_YUK(O(cont:cont+2,i_mc))
							cont=cont+3
						ELSE
							CALL derivata_Jee_YUK(O(cont:cont+1,i_mc))
							cont=cont+2
						END IF
					END IF
				ELSE
					IF (split_Aee) THEN
						CALL derivata_Jee_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					ELSE
						CALL derivata_Jee_YUK(O(cont:cont,i_mc))
						cont=cont+1
					END IF
				END IF
			ELSE
				IF (opt_F_Jee) THEN
					IF (split_Fee) THEN
						CALL derivata_Jee_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					ELSE
						CALL derivata_Jee_YUK(O(cont:cont,i_mc))
						cont=cont+1
					END IF
				END IF
			END IF
		END SELECT
		
		SELECT CASE (Jep_kind)
		CASE ('yuk')
			IF (opt_A_Jep) THEN
				IF (opt_F_Jep) THEN
					IF (split_Aep) THEN
						IF (split_Fep) THEN
							CALL derivata_Jep_YUK(O(cont:cont+3,i_mc))
							cont=cont+4
						ELSE
							CALL derivata_Jep_YUK(O(cont:cont+2,i_mc))
							cont=cont+3
						END IF
					ELSE
						IF (split_Fep) THEN
							CALL derivata_Jep_YUK(O(cont:cont+2,i_mc))
							cont=cont+3
						ELSE
							CALL derivata_Jep_YUK(O(cont:cont+1,i_mc))
							cont=cont+2
						END IF
					END IF
				ELSE
					IF (split_Aep) THEN
						CALL derivata_Jep_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					ELSE
						CALL derivata_Jep_YUK(O(cont:cont,i_mc))
						cont=cont+1
					END IF
				END IF
			ELSE
				IF (opt_F_Jep) THEN
					IF (split_Fep) THEN
						CALL derivata_Jep_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					ELSE
						CALL derivata_Jep_YUK(O(cont:cont,i_mc))
						cont=cont+1
					END IF
				END IF
			END IF
		CASE ('yup')
			IF (opt_A_Jep) THEN
				IF (opt_F_Jep) THEN
					IF (split_Aep) THEN
						IF (split_Fep) THEN
							CALL derivata_Jep_YUK(O(cont:cont+3,i_mc))
							cont=cont+4
						ELSE
							CALL derivata_Jep_YUK(O(cont:cont+2,i_mc))
							cont=cont+3
						END IF
					ELSE
						IF (split_Fep) THEN
							CALL derivata_Jep_YUK(O(cont:cont+2,i_mc))
							cont=cont+3
						ELSE
							CALL derivata_Jep_YUK(O(cont:cont+1,i_mc))
							cont=cont+2
						END IF
					END IF
				ELSE
					IF (split_Aep) THEN
						CALL derivata_Jep_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					ELSE
						CALL derivata_Jep_YUK(O(cont:cont,i_mc))
						cont=cont+1
					END IF
				END IF
			ELSE
				IF (opt_F_Jep) THEN
					IF (split_Fep) THEN
						CALL derivata_Jep_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					ELSE
						CALL derivata_Jep_YUK(O(cont:cont,i_mc))
						cont=cont+1
					END IF
				END IF
			END IF
		CASE ('atm')
			CALL derivata_Jep_ATM(O(cont:cont,i_mc))
			cont=cont+1
		CASE ('atp')
			CALL derivata_Jep_ATM(O(cont:cont,i_mc))
			cont=cont+1
		END SELECT
		
		IF (opt_Kse) THEN
			SELECT CASE (Kse_kind)
			CASE ('gss')
				CALL derivata_KERNese(O(cont,i_mc))
				cont=cont+1
			CASE ('gsc')
				CALL derivata_KERNese(O(cont,i_mc))
				cont=cont+1
			CASE ('gsp')
				CALL derivata_KERNese(O(cont,i_mc))
				cont=cont+1
			CASE ('gsd')
				CALL derivata_KERNese(O(cont,i_mc))
				cont=cont+1
			CASE ('gdc')
				CALL derivata_KERNese(O(cont,i_mc))
				cont=cont+1
			CASE ('gdp')
				CALL derivata_KERNese(O(cont,i_mc))
				cont=cont+1
			CASE ('atm')
				CALL derivata_atmKERNese(O(cont,i_mc))
				cont=cont+1
			END SELECT
		END IF
		
		IF (opt_Jse) THEN
			SELECT CASE (Jse_kind)
			CASE ('pot')
				CALL derivata_Jsese_POT(O(cont:cont+1,i_mc))
				cont=cont+2
			CASE ('bou')
				STOP 'Jsese bou non ancora implementato per lo SR &
				[ module_variational_calculations.f90 > calcola_estimatori_per_derivate_variazionali ]'
			CASE ('ppb')
				STOP 'Jsese ppb non ancora implementato per lo SR &
				[ module_variational_calculations.f90 > calcola_estimatori_per_derivate_variazionali ]'
			CASE ('yuk')
				IF (split_Asese) THEN
					IF (split_Fsese) THEN
						CALL derivata_Jsese_YUK(O(cont:cont+3,i_mc))
						cont=cont+4
					ELSE
						CALL derivata_Jsese_YUK(O(cont:cont+2,i_mc))
						cont=cont+3
					END IF
				ELSE
					IF (split_Fsese) THEN
						CALL derivata_Jsese_YUK(O(cont:cont+2,i_mc))
						cont=cont+3
					ELSE
						CALL derivata_Jsese_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					END IF
				END IF
			CASE ('yup')
				IF (split_Asese) THEN
					IF (split_Fsese) THEN
						CALL derivata_Jsese_YUK(O(cont:cont+3,i_mc))
						cont=cont+4
					ELSE
						CALL derivata_Jsese_YUK(O(cont:cont+2,i_mc))
						cont=cont+3
					END IF
				ELSE
					IF (split_Fsese) THEN
						CALL derivata_Jsese_YUK(O(cont:cont+2,i_mc))
						cont=cont+3
					ELSE
						CALL derivata_Jsese_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					END IF
				END IF
			END SELECT
		END IF
		
		IF (opt_Jsesp) THEN
			SELECT CASE (Jsesp_kind)
			CASE ('pot')
				STOP 'Jsesp pot non ancora implementato per lo SR &
				[ module_variational_calculations.f90 > calcola_estimatori_per_derivate_variazionali ]'
			CASE ('yuk')
				IF (split_Asesp) THEN
					IF (split_Fsesp) THEN
						CALL derivata_Jsesp_YUK(O(cont:cont+3,i_mc))
						cont=cont+4
					ELSE
						CALL derivata_Jsesp_YUK(O(cont:cont+2,i_mc))
						cont=cont+3
					END IF
				ELSE
					IF (split_Fsesp) THEN
						CALL derivata_Jsesp_YUK(O(cont:cont+2,i_mc))
						cont=cont+3
					ELSE
						CALL derivata_Jsesp_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					END IF
				END IF
			CASE ('yup')
				IF (split_Asesp) THEN
					IF (split_Fsesp) THEN
						CALL derivata_Jsesp_YUK(O(cont:cont+3,i_mc))
						cont=cont+4
					ELSE
						CALL derivata_Jsesp_YUK(O(cont:cont+2,i_mc))
						cont=cont+3
					END IF
				ELSE
					IF (split_Fsesp) THEN
						CALL derivata_Jsesp_YUK(O(cont:cont+2,i_mc))
						cont=cont+3
					ELSE
						CALL derivata_Jsesp_YUK(O(cont:cont+1,i_mc))
						cont=cont+2
					END IF
				END IF
			CASE ('gss')
				CALL derivata_KERNsesp(O(cont,i_mc))
				cont=cont+1
			CASE ('gsd')
				CALL derivata_KERNsesp(O(cont,i_mc))
				cont=cont+1
			END SELECT
		END IF
		
		IF ( opt_SDe ) THEN
			SELECT CASE (SDe_kind)
			CASE ('prf')
				IF (opt_c_eff_hartree) THEN
					CALL derivata_SDe_HAR(O(cont:cont+N_pw_lda*N_pw_lda-1,i_mc))
					cont=cont+N_pw_lda*N_pw_lda
				END IF
			CASE ('fre')
				IF (opt_c_eff_hartree) THEN
					CALL derivata_SDe_HAR(O(cont:cont+N_pw_lda*N_pw_lda-1,i_mc))
					cont=cont+N_pw_lda*N_pw_lda
				END IF
			CASE ('atm')
				CALL derivata_SDe_atm(O(cont,i_mc))
				cont=cont+1
			CASE ('bat')
				CALL derivata_SDe_bat(O(cont,i_mc))
				cont=cont+1
			CASE ('atp')
				CALL derivata_SDe_atp(O(cont,i_mc))
				cont=cont+1
			END SELECT
		END IF
		
		IF ( SDse_kind/='no_' ) THEN
			IF ( opt_SDse ) THEN
				SELECT CASE (SDse_kind)
				CASE ('atm')
					CALL derivata_SDse_atm(O(cont,i_mc))
					cont=cont+1
				CASE ('atp')
					CALL derivata_SDse_atp(O(cont,i_mc))
					cont=cont+1
				END SELECT
			END IF
		END IF
		
		IF ( opt_rp ) THEN
			CALL derivata_psi_Rp(O(cont:cont+3*N_part-1,i_mc))
			!DO i = 1, N_part, 1
			!	PRINT * , i, O(i*3-2:i*3,i_mc)
			!END DO
			cont=cont+3*N_part
		END IF
		
	END SUBROUTINE calcola_termini_derivate_variazionali
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatori_var(i_wf,i_mc)
		USE estimatori
		USE calcola_accettazione
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: i_wf
		INTEGER (KIND=8), INTENT(IN) :: i_mc
		REAL (KIND=8) :: dE_kin, dE_JF, dE_pot
		REAL (KIND=8) :: Uep_, u_ep_(1:N_part,1:N_part)
		
		IF (.NOT.iniz_variational_calculations) STOP 'Prima di calcolare gli estimatori devi inizializzare &
		  [ module_variational_calculations.f90 > calcola_estimatori_var ]'
		
		vwf_data(i_wf)%w(i_mc)=REAL( DEXP(-(vwf_data(i_wf)%Uee-Uee_old)-(vwf_data(i_wf)%Uep-Uep_old)- &
		  (vwf_data(i_wf)%Use1-Use1_old)-(vwf_data(i_wf)%Use2-Use2_old)- &
		  (vwf_data(i_wf)%Uese1-Uese1_old)-(vwf_data(i_wf)%Uese2-Uese2_old)- &
		  (vwf_data(i_wf)%Usesp1-Usesp1_old)-(vwf_data(i_wf)%Usesp2-Usesp2_old)) ,8)* &
		  (vwf_data(i_wf)%detGDse1_up/detGDse1_up_old)*(vwf_data(i_wf)%detGDse1_dw/detGDse1_dw_old)* &
		  (vwf_data(i_wf)%detGDse2_up/detGDse2_up_old)*(vwf_data(i_wf)%detGDse2_dw/detGDse2_dw_old)* &
		  (vwf_data(i_wf)%detGDsp1_up/detGDsp1_up_old)*(vwf_data(i_wf)%detGDse1_dw/detGDsp1_dw_old)* &
		  (vwf_data(i_wf)%detGDsp2_up/detGDsp2_up_old)*(vwf_data(i_wf)%detGDse2_dw/detGDsp2_dw_old)
		
		CALL setta_parametri(parametri_var(1:num_par,i_wf),num_par)
		CALL energia_cinetica(dE_kin,dE_JF)
		CALL energia_potenziale(dE_pot)
		vwf_data(i_wf)%E_tot(i_mc)=dE_kin+dE_pot
		CALL setta_parametri(parametri_var(1:num_par,0),num_par)
				
	END SUBROUTINE calcola_estimatori_var
!-----------------------------------------------------------------------
	SUBROUTINE chiudi_variational_calculations()
		IMPLICIT NONE
		INTEGER :: i
		
		IF (.NOT.iniz_variational_calculations) STOP 'Prima di chiudere il modulo devi inizializzarlo &
		  [ module_variational_calculations.f90 > calcola_estimatori_var ]'
		
		IF (what_to_do=='congrad') THEN
			DO i = 1, num_wf, 1
				IF ((SDe_kind/='no_').AND.(SDe_kind/='pw_').AND.(SDe_kind/='lda')) THEN
					DEALLOCATE(vwf_data(i)%SDe_up)
					DEALLOCATE(vwf_data(i)%ISDe_up)
					DEALLOCATE(vwf_data(i)%pvte_up)
					DEALLOCATE(vwf_data(i)%SDe_dw)
					DEALLOCATE(vwf_data(i)%ISDe_dw)
					DEALLOCATE(vwf_data(i)%pvte_dw)
				END IF
				IF (Jee_kind/='no_') THEN
					DEALLOCATE(vwf_data(i)%u_ee)
				END IF
				IF (Jep_kind/='no_') THEN
					DEALLOCATE(vwf_data(i)%u_ep)
				END IF
				IF ((SDse_kind/='no_').AND.(SDse_kind/='pw_').AND.(SDse_kind/='lda')) THEN
					DEALLOCATE(vwf_data(i)%SDse1_up)
					DEALLOCATE(vwf_data(i)%ISDse1_up)
					DEALLOCATE(vwf_data(i)%pvtse1_up)
					DEALLOCATE(vwf_data(i)%SDse1_dw)
					DEALLOCATE(vwf_data(i)%ISDse1_dw)
					DEALLOCATE(vwf_data(i)%pvtse1_dw)
					DEALLOCATE(vwf_data(i)%SDse2_up)
					DEALLOCATE(vwf_data(i)%ISDse2_up)
					DEALLOCATE(vwf_data(i)%pvtse2_up)
					DEALLOCATE(vwf_data(i)%SDse2_dw)
					DEALLOCATE(vwf_data(i)%ISDse2_dw)
					DEALLOCATE(vwf_data(i)%pvtse2_dw)
				END IF
				IF (Jse_kind/='no_') THEN
					IF (Jse_kind/='pot') THEN
						DEALLOCATE(vwf_data(i)%u_se1)
						DEALLOCATE(vwf_data(i)%u_se2)
					ELSE IF (Jse_kind=='pot') THEN
						DEALLOCATE(vwf_data(i)%u_POT_se1)
						DEALLOCATE(vwf_data(i)%u_POT_se2)
					END IF
				END IF
				IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
					DEALLOCATE(vwf_data(i)%u_ese1)
					DEALLOCATE(vwf_data(i)%u_ese2)
				ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
					DEALLOCATE(vwf_data(i)%GDse1_up)
					DEALLOCATE(vwf_data(i)%IGDse1_up)
					DEALLOCATE(vwf_data(i)%pvtgdse1_up)
					DEALLOCATE(vwf_data(i)%GDse1_dw)
					DEALLOCATE(vwf_data(i)%IGDse1_dw)
					DEALLOCATE(vwf_data(i)%pvtgdse1_dw)
					DEALLOCATE(vwf_data(i)%GDse2_up)
					DEALLOCATE(vwf_data(i)%IGDse2_up)
					DEALLOCATE(vwf_data(i)%pvtgdse2_up)
					DEALLOCATE(vwf_data(i)%GDse2_dw)
					DEALLOCATE(vwf_data(i)%IGDse2_dw)
					DEALLOCATE(vwf_data(i)%pvtgdse2_dw)
				END IF
				IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
					DEALLOCATE(vwf_data(i)%u_sesp1)
					DEALLOCATE(vwf_data(i)%u_sesp2)
				ELSE IF (Jsesp_kind=='gsd') THEN
					DEALLOCATE(vwf_data(i)%GDsp1_up)
					DEALLOCATE(vwf_data(i)%IGDsp1_up)
					DEALLOCATE(vwf_data(i)%pvtgdsp1_up)
					DEALLOCATE(vwf_data(i)%GDsp1_dw)
					DEALLOCATE(vwf_data(i)%IGDsp1_dw)
					DEALLOCATE(vwf_data(i)%pvtgdsp1_dw)
					DEALLOCATE(vwf_data(i)%GDsp2_up)
					DEALLOCATE(vwf_data(i)%IGDsp2_up)
					DEALLOCATE(vwf_data(i)%pvtgdsp2_up)
					DEALLOCATE(vwf_data(i)%GDsp2_dw)
					DEALLOCATE(vwf_data(i)%IGDsp2_dw)
					DEALLOCATE(vwf_data(i)%pvtgdsp2_dw)
				END IF
			END DO
			DEALLOCATE(vwf_data)
			DEALLOCATE(parametri_var)
			iniz_variational_calculations=.FALSE.
			flag_gradiente=.FALSE.
		ELSE IF ((what_to_do=='stocrec').OR.(what_to_do=='stoc_ns').OR.(what_to_do=='stoc_av').OR.(what_to_do=='pure_sr')) THEN
			flag_derivate_var=.FALSE.
			DEALLOCATE(O)
		END IF
		
	END SUBROUTINE chiudi_variational_calculations
	
END MODULE variational_calculations
