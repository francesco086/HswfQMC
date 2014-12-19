MODULE estimatori
	USE dati_fisici
	USE dati_mc
	USE momenta
	USE funzione_onda
	USE walkers
	USE calcola_accettazione
	USE generic_tools
	IMPLICIT NONE
	LOGICAL, PRIVATE, SAVE :: iniz_estimatori=.FALSE., verbose_mode=.FALSE.
	REAL (KIND=8), PARAMETER, PRIVATE :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: E_tot(:), E_kin(:), E_JF(:), E_pot(:), w(:)
	
	CONTAINS
	
	SUBROUTINE inizializza_estimatori()
		IMPLICIT NONE
		
		IF (flag_E_tot) THEN
			ALLOCATE(E_tot(1:N_mc))
			E_tot=0.d0
		END IF
		IF (flag_E_kin) THEN
			ALLOCATE(E_kin(1:N_mc), E_JF(1:N_mc))
			E_kin=0.d0
			E_JF=0.d0
		END IF
		IF (flag_E_pot) THEN
			ALLOCATE(E_pot(1:N_mc))
			E_pot=0.d0
		END IF
		ALLOCATE(w(1:N_mc))
		w=0.d0
		
		iniz_estimatori=.TRUE.
	END SUBROUTINE inizializza_estimatori
!-----------------------------------------------------------------------
	SUBROUTINE valuta_estimatori(i)
		IMPLICIT NONE
		INTEGER (KIND=8), INTENT(IN) :: i
		REAL (KIND=8) :: app1, app2
		IF (.NOT. iniz_estimatori) STOP 'Prima di valutare gli estimatori devi inizializzarli &
		  [ module_estimatori.f90 > valuta_estimatori ]'
		
		IF (flag_E_kin) CALL energia_cinetica(E_kin(i),E_JF(i))
		IF (flag_E_pot) CALL energia_potenziale(E_pot(i))
		IF (flag_E_tot) THEN
			IF (flag_E_kin) THEN
				E_tot(i)=E_kin(i)
			ELSE
				CALL energia_cinetica(app1,app2)
				E_tot(i)=app1
			END IF
			IF (flag_E_pot) THEN
				E_tot(i)=E_tot(i)+E_pot(i)
			ELSE
				CALL energia_potenziale(app1)
				E_tot(i)=E_tot(i)+app1
			END IF
		END IF
		CALL peso_funzione_onda(w(i))
				
	END SUBROUTINE valuta_estimatori
!-----------------------------------------------------------------------
	SUBROUTINE mantieni_stessi_estimatori(i)
		IMPLICIT NONE
		INTEGER (KIND=8), INTENT(IN) :: i
		IF (.NOT. iniz_estimatori) STOP 'Prima di confermare gli estimatori devi inizializzarli &
		  [ module_estimatori.f90 > valuta_estimatori ]'
		
		IF (flag_E_kin) THEN
			E_kin(i)=E_kin(i-1)
			E_JF(i)=E_JF(i-1)
		END IF
		IF (flag_E_pot) E_pot(i)=E_pot(i-1)
		IF (flag_E_tot) E_tot(i)=E_tot(i-1)
		w(i)=w(i-1)
		
	END SUBROUTINE mantieni_stessi_estimatori
!-----------------------------------------------------------------------
	SUBROUTINE trascrivi_estimatori(codice)
		IMPLICIT NONE
		LOGICAL :: flag_scrivi, flag_continua_rank
		CHARACTER(LEN=4) :: istring
		CHARACTER(LEN=*) :: codice
		INTEGER :: status(MPI_STATUS_SIZE)
		INTEGER :: j
		INTEGER (KIND=8) :: i
		REAL (KIND=8) :: vec_save(1:N_mc)
		
		WRITE (istring, '(I4.4)'), mpi_myrank
		
		IF ((.NOT. flag_disk) .AND. (mpi_myrank==0)) STOP 'Stai trascrivendo gli estimatori quando non dovrebbe succedere &
		  [ module_VMC.f90 > conferma_estimatori ]'
		
		IF (.NOT. iniz_estimatori) STOP 'Prima di salvare gli estimatori devi inizializzarli &
		  [ module_VMC.f90 > conferma_estimatori ]'
		
		IF (mpi_myrank==0) THEN
			flag_scrivi=.TRUE.
			flag_continua_rank=flag_continua
		ELSE
			flag_scrivi=.FALSE.
			flag_continua_rank=.TRUE.
		END IF
		
		DO j = 0, mpi_nprocs-1, 1
			IF (mpi_myrank==j .AND. flag_scrivi) THEN
				IF (verbose_mode) PRINT *, 'ESTIMATORI: processo ', mpi_myrank, ' trascrive i suoi dati'
				IF (flag_E_kin) THEN
					!IF (SDse_kind=='no_') THEN
					!	CALL salva_vettore_dati_bin(E_kin,N_mc,N_AV,flag_continua_rank,'estimatori/E_kin'//codice)
					!	IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(E_kin,N_mc,10000_8,'estimatori/E_kin-sample'//istring//codice)
					!	CALL salva_vettore_dati_bin(E_JF,N_mc,N_AV,flag_continua_rank,'estimatori/E_JF'//codice)
					!	IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(E_JF,N_mc,10000_8,'estimatori/E_JF-sample'//istring//codice)
					!ELSE
						DO i = 1, N_mc, 1
							vec_save(i)=E_kin(i)*w(i)
						END DO
						CALL salva_vettore_dati_bin(vec_save,N_mc,N_AV,flag_continua_rank,'estimatori/E_kin'//codice)
						IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(vec_save,N_mc,10000_8,'estimatori/E_kin-sample'//istring//codice)
						DO i = 1, N_mc, 1
							vec_save(i)=E_JF(i)*w(i)
						END DO
						CALL salva_vettore_dati_bin(vec_save,N_mc,N_AV,flag_continua_rank,'estimatori/E_JF'//codice)
						IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(vec_save,N_mc,10000_8,'estimatori/E_JF-sample'//istring//codice)
					!END IF
				END IF
				IF (flag_E_pot) THEN
					!IF (SDse_kind=='no_') THEN
					!	CALL salva_vettore_dati_bin(E_pot,N_mc,N_AV,flag_continua_rank,'estimatori/E_pot'//codice)
					!	IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(E_pot,N_mc,10000_8,'estimatori/E_pot-sample'//istring//codice)
					!ELSE
						DO i = 1, N_mc, 1
							vec_save(i)=E_pot(i)*w(i)
						END DO
						CALL salva_vettore_dati_bin(vec_save,N_mc,N_AV,flag_continua_rank,'estimatori/E_pot'//codice)
						IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(vec_save,N_mc,10000_8,'estimatori/E_pot-sample'//istring//codice)
					!END IF
				END IF
				IF (flag_E_tot) THEN
					!IF (SDse_kind=='no_') THEN
					!	CALL salva_vettore_dati_bin(E_tot,N_mc,N_AV,flag_continua_rank,'estimatori/E_tot'//codice)
					!	IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(E_tot,N_mc,10000_8,'estimatori/E_tot-sample'//istring//codice)
					!ELSE
						DO i = 1, N_mc, 1
							vec_save(i)=E_tot(i)*w(i)
						END DO
						CALL salva_vettore_dati_bin(vec_save,N_mc,N_AV,flag_continua_rank,'estimatori/E_tot'//codice)
						IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(vec_save,N_mc,10000_8,'estimatori/E_tot-sample'//istring//codice)
					!END IF
				END IF
				!IF (SDse_kind/='no_') THEN
					CALL salva_vettore_dati_bin(w,N_mc,N_AV,flag_continua_rank,'estimatori/w'//codice)
					IF (mpi_myrank==0) CALL salva_vettore_dati_dat_ridotto(w,N_mc,10000_8,'estimatori/w-sample'//istring//codice)
				!END IF
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
		
		!PRINT * , 'finito!', mpi_myrank
		
		IF (flag_mpi) CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
		
	END SUBROUTINE trascrivi_estimatori
!-----------------------------------------------------------------------

	SUBROUTINE energia_cinetica(E_kin,E_JF)
		IMPLICIT NONE
		INTEGER :: i, j, i3, idw, jdw, ik, info, i_SD, j_SD, alpha, beta, gamm, il, iadd, spin, t, m, n
      REAL(KIND=8) :: q(0:3,H_N_part),rq(0:3,H_N_part,H_N_part),SDe(H_N_part,H_N_part)
      REAL(KIND=8) :: sigm(N_part,N_part),sigm1(N_part,N_part),sigm2(N_part,N_part)
      REAL(KIND=8) :: Gvrq(3,3,H_N_part,H_N_part,H_N_part)
      REAL(KIND=8) :: Grq(3,H_N_part,H_N_part,H_N_part,2), Gphi(3,H_N_part,H_N_part,H_N_part,2)
      REAL(KIND=8) :: Lrq(3,H_N_part,H_N_part,H_N_part,2), Lphi(3,H_N_part,H_N_part,H_N_part,2)
      REAL(KIND=8) :: Dd(3,N_part,N_part), D2d(3,N_part,N_part)
		REAL (KIND=8) :: frf1(0:3), frf2(1:3), frf3, frf5, frf6, mix_prod, norm
		REAL (KIND=8) :: frfs1(0:3), frfs2(0:3), frf2s1(1:3), frf2s2(1:3), frf3s1(1:3), frf3s2(1:3)
		REAL (KIND=8) :: uee1, uee2, uep1, uep2
		REAL (KIND=8) :: gee(1:3,1:N_part), lee
		REAL (KIND=8) :: gep(1:3,1:N_part), lep
		REAL (KIND=8) :: gkernse1(1:3,1:N_part), gkernse2(1:3,1:N_part), lkernse
		REAL (KIND=8) :: work_R(1:3*H_N_part)
		REAL (KIND=8) :: minL, L38
		COMPLEX (KIND=8) :: lsdee, gsdee_up(1:3,1:H_N_part), gsdee_dw(1:3,1:H_N_part)
		COMPLEX (KIND=8) :: work_C(1:3*H_N_part)
		COMPLEX (KIND=8) :: der1_up(1:3), der2_up, der1_dw(1:3), der2_dw, frfc1, frfc2, frfc3(1:3), frfc4(1:3)
		REAL (KIND=8), INTENT(OUT) :: E_kin, E_JF
		
		E_kin=0.d0
		E_JF=0.d0
		
		elettroni: IF (flag_elettroni) THEN
		
		IF (howtomove=='allp') THEN
			IF (SDe_kind/='no_') THEN
				CALL ZGETRI( H_N_part, ISDe_up_old, H_N_part, pvte_up_old, work_C, 3*H_N_part, info )
				IF (info/=0) THEN
					PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA UP. INFO=', info
					STOP
				END IF
				CALL ZGETRI( H_N_part, ISDe_dw_old, H_N_part, pvte_dw_old, work_C, 3*H_N_part, info )
				IF (info/=0) THEN
					PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA DOWN. INFO=', info
					STOP
				END IF
			END IF
			IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc')) THEN	
				CALL DGETRI( H_N_part, IGDse1_up_old, H_N_part, pvtgdse1_up_old, work_R, 3*H_N_part, info )
				IF (info/=0) THEN
					PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GD1 UP. INFO=', info
					STOP
				END IF
				CALL DGETRI( H_N_part, IGDse1_dw_old, H_N_part, pvtgdse1_dw_old, work_R, 3*H_N_part, info )
				IF (info/=0) THEN
					PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GD1 DOWN. INFO=', info
					STOP
				END IF
				CALL DGETRI( H_N_part, IGDse2_up_old, H_N_part, pvtgdse2_up_old, work_R, 3*H_N_part, info )
				IF (info/=0) THEN
					PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GD2 UP. INFO=', info
					STOP
				END IF
				CALL DGETRI( H_N_part, IGDse2_dw_old, H_N_part, pvtgdse2_dw_old, work_R, 3*H_N_part, info )
				IF (info/=0) THEN
					PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GD2 DOWN. INFO=', info
					STOP
				END IF
			END IF
		END IF
		
		!calcolo la parte del Determinante di Slater
		lsdee=0.d0
		gsdee_up=0.d0
		gsdee_dw=0.d0
		SELECT CASE (SDe_kind)
		CASE ('pw_')
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					lsdee=lsdee-k_pw(0,i)*SDe_up_old(j,i)*ISDe_up_old(i,j)           !posso mettere k(0)
					gsdee_up(1:3,j)=gsdee_up(1:3,j)+(0.d0,1.d0)*k_pw(1:3,i)*SDe_up_old(j,i)*ISDe_up_old(i,j)
					lsdee=lsdee-k_pw(0,i)*SDe_dw_old(j,i)*ISDe_dw_old(i,j)
					gsdee_dw(1:3,j)=gsdee_dw(1:3,j)+(0.d0,1.d0)*k_pw(1:3,i)*SDe_dw_old(j,i)*ISDe_dw_old(i,j)
				END DO
			END DO
		CASE ('lda')
			IF ( flag_simm_lda ) THEN
				DO j = 1, H_N_part, 1
					DO i = 1, H_N_part, 1
						der1_up=0.d0
						der1_dw=0.d0
						der2_up=0.d0
						der2_dw=0.d0
						!DO ik = 1, N_pw_lda, 1
						!	frfc1=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_pw_lda(1:3,ik),re_old(1:3,j)))
						!	frfc2=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_pw_lda(1:3,ik),re_old(1:3,j+H_N_part)))
						!	der1_up(1:3)=der1_up(1:3)+(0.d0,1.d0)*k_pw_lda(1:3,ik)*fattori_pw_lda(ik,i)*frfc1 - &
						!	                          (0.d0,1.d0)*k_pw_lda(1:3,ik)*DCONJG(fattori_pw_lda(ik,i)*frfc1)
						!	der1_dw(1:3)=der1_dw(1:3)+(0.d0,1.d0)*k_pw_lda(1:3,ik)*fattori_pw_lda(ik,i)*frfc2 - &
						!	                          (0.d0,1.d0)*k_pw_lda(1:3,ik)*DCONJG(fattori_pw_lda(ik,i)*frfc2)
						!	der2_up=der2_up-k_pw_lda(0,ik)*(fattori_pw_lda(ik,i)*frfc1+DCONJG(fattori_pw_lda(ik,i)*frfc1))
						!	der2_dw=der2_dw-k_pw_lda(0,ik)*(fattori_pw_lda(ik,i)*frfc2+DCONJG(fattori_pw_lda(ik,i)*frfc2))
						!END DO
						DO ik = 1, num_fpw_lda(i), 1
							frfc1=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,i),re_old(1:3,j)))
							frfc2=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,i),re_old(1:3,j+H_N_part)))
							der1_up(1:3)=der1_up(1:3)+(0.d0,1.d0)*k_fpw_lda(1:3,ik,i)*fattori_fpw_lda(ik,i)*frfc1 - &
							                          (0.d0,1.d0)*k_fpw_lda(1:3,ik,i)*DCONJG(fattori_fpw_lda(ik,i)*frfc1)
							der1_dw(1:3)=der1_dw(1:3)+(0.d0,1.d0)*k_fpw_lda(1:3,ik,i)*fattori_fpw_lda(ik,i)*frfc2 - &
							                          (0.d0,1.d0)*k_fpw_lda(1:3,ik,i)*DCONJG(fattori_fpw_lda(ik,i)*frfc2)
							der2_up=der2_up-k_fpw_lda(0,ik,i)*(fattori_fpw_lda(ik,i)*frfc1+DCONJG(fattori_fpw_lda(ik,i)*frfc1))
							der2_dw=der2_dw-k_fpw_lda(0,ik,i)*(fattori_fpw_lda(ik,i)*frfc2+DCONJG(fattori_fpw_lda(ik,i)*frfc2))
						END DO
						gsdee_up(1:3,j)=gsdee_up(1:3,j)+der1_up*ISDe_up_old(i,j)
						gsdee_dw(1:3,j)=gsdee_dw(1:3,j)+der1_dw*ISDe_dw_old(i,j)
						lsdee=lsdee+(der2_up*ISDe_up_old(i,j)+der2_dw*ISDe_dw_old(i,j))
					END DO
				END DO
			ELSE
				DO j = 1, H_N_part, 1
					DO i = 1, H_N_part, 1
						der1_up=0.d0
						der1_dw=0.d0
						der2_up=0.d0
						der2_dw=0.d0
						DO ik = 1, num_fpw_lda(i), 1
							frfc1=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,i),re_old(1:3,j)))
							frfc2=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,i),re_old(1:3,j+H_N_part)))
							der1_up(1:3)=der1_up(1:3)+(0.d0,1.d0)*k_fpw_lda(1:3,ik,i)*fattori_fpw_lda(ik,i)*frfc1
							der1_dw(1:3)=der1_dw(1:3)+(0.d0,1.d0)*k_fpw_lda(1:3,ik,i)*fattori_fpw_lda(ik,i)*frfc2
							der2_up=der2_up-k_fpw_lda(0,ik,i)*(fattori_fpw_lda(ik,i)*frfc1)
							der2_dw=der2_dw-k_fpw_lda(0,ik,i)*(fattori_fpw_lda(ik,i)*frfc2)
						END DO
						gsdee_up(1:3,j)=gsdee_up(1:3,j)+der1_up*ISDe_up_old(i,j)
						gsdee_dw(1:3,j)=gsdee_dw(1:3,j)+der1_dw*ISDe_dw_old(i,j)
						lsdee=lsdee+(der2_up*ISDe_up_old(i,j)+der2_dw*ISDe_dw_old(i,j))
					END DO
				END DO
			END IF
		CASE ('prf')
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					der1_up=0.d0
					der1_dw=0.d0
					der2_up=0.d0
					der2_dw=0.d0
					DO ik = 1, num_pw_dnfH(i), 1
						frfc1=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_pw_dnfH(1:3,ik,i),re_old(1:3,j)))
						frfc2=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_pw_dnfH(1:3,ik,i),re_old(1:3,j+H_N_part)))
						der1_up(1:3)=der1_up(1:3)+(0.d0,1.d0)*k_pw_dnfH(1:3,ik,i)* &
						  fattori_pw_dnfH(ik,i)*frfc1
						der1_dw(1:3)=der1_dw(1:3)+(0.d0,1.d0)*k_pw_dnfH(1:3,ik,i)* &
						  fattori_pw_dnfH(ik,i)*frfc2
						der2_up=der2_up-k_pw_dnfH(0,ik,i)*(fattori_pw_dnfH(ik,i)*frfc1)
						der2_dw=der2_dw-k_pw_dnfH(0,ik,i)*(fattori_pw_dnfH(ik,i)*frfc2)
					END DO
					gsdee_up(1:3,j)=gsdee_up(1:3,j)+der1_up*ISDe_up_old(i,j)
					gsdee_dw(1:3,j)=gsdee_dw(1:3,j)+der1_dw*ISDe_dw_old(i,j)
					lsdee=lsdee+(der2_up*ISDe_up_old(i,j)+der2_dw*ISDe_dw_old(i,j))
				END DO
			END DO
		CASE ('fre')
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					der1_up=0.d0
					der1_dw=0.d0
					der2_up=0.d0
					der2_dw=0.d0
					DO ik = 1, num_pw_dnfH(i), 1
						frfc1=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_pw_dnfH(1:3,ik,i),re_old(1:3,j)))
						frfc2=CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_pw_dnfH(1:3,ik,i),re_old(1:3,j+H_N_part)))
						der1_up(1:3)=der1_up(1:3)+(0.d0,1.d0)*k_pw_dnfH(1:3,ik,i)* &
						  fattori_pw_dnfH(ik,i)*frfc1
						der1_dw(1:3)=der1_dw(1:3)+(0.d0,1.d0)*k_pw_dnfH(1:3,ik,i)* &
						  fattori_pw_dnfH(ik,i)*frfc2
						der2_up=der2_up-k_pw_dnfH(0,ik,i)*(fattori_pw_dnfH(ik,i)*frfc1)
						der2_dw=der2_dw-k_pw_dnfH(0,ik,i)*(fattori_pw_dnfH(ik,i)*frfc2)
					END DO
					gsdee_up(1:3,j)=gsdee_up(1:3,j)+der1_up*ISDe_up_old(i,j)
					gsdee_dw(1:3,j)=gsdee_dw(1:3,j)+der1_dw*ISDe_dw_old(i,j)
					lsdee=lsdee+(der2_up*ISDe_up_old(i,j)+der2_dw*ISDe_dw_old(i,j))
				END DO
			END DO
		CASE ('atm')
			DO i = 1, H_N_part, 1
				i_SD=i+H_N_part
				DO j = 1, H_N_part, 1
					j_SD=j+H_N_part
					der1_up(1:3)=-rij_ep_old(1:3,j,i)*C_atm/rij_ep_old(0,j,i)
					der1_dw(1:3)=-rij_ep_old(1:3,j_SD,i_SD)*C_atm/rij_ep_old(0,j_SD,i_SD)
					der2_up=-2.d0*C_atm/rij_ep_old(0,j,i)+C_atm*C_atm
					der2_dw=-2.d0*C_atm/rij_ep_old(0,j_SD,i_SD)+C_atm*C_atm
					gsdee_up(1:3,j)=gsdee_up(1:3,j)+SDe_up_old(j,i)*der1_up*ISDe_up_old(i,j)
					gsdee_dw(1:3,j)=gsdee_dw(1:3,j)+SDe_dw_old(j,i)*der1_dw*ISDe_dw_old(i,j)
					lsdee=lsdee+(SDe_up_old(j,i)*der2_up*ISDe_up_old(i,j)+SDe_dw_old(j,i)*der2_dw*ISDe_dw_old(i,j))
				END DO
			END DO
		CASE ('bat')
			DO i = 1, H_N_part, 1
				i_SD=i+H_N_part
				DO j = 1, H_N_part, 1
					j_SD=j+H_N_part
					der1_up(1:3)=-(rij_ep_old(1:3,j,i)*C_atm/rij_ep_old(0,j,i))*DEXP(-C_atm*rij_ep_old(0,j,i))  & 
					             -(rij_ep_old(1:3,j,i_SD)*C_atm/rij_ep_old(0,j,i_SD))*DEXP(-C_atm*rij_ep_old(0,j,i_SD))
					der1_dw(1:3)=-(rij_ep_old(1:3,j_SD,i_SD)*C_atm/rij_ep_old(0,j_SD,i_SD))*DEXP(-C_atm*rij_ep_old(0,j_SD,i_SD))  &
					             -(rij_ep_old(1:3,j_SD,i)*C_atm/rij_ep_old(0,j_SD,i))*DEXP(-C_atm*rij_ep_old(0,j_SD,i))
					der2_up=(-2.d0*C_atm/rij_ep_old(0,j,i)+C_atm*C_atm)*DEXP(-C_atm*rij_ep_old(0,j,i))  &
					             +(-2.d0*C_atm/rij_ep_old(0,j,i_SD)+C_atm*C_atm)*DEXP(-C_atm*rij_ep_old(0,j,i_SD))
					der2_dw=(-2.d0*C_atm/rij_ep_old(0,j_SD,i_SD)+C_atm*C_atm)*DEXP(-C_atm*rij_ep_old(0,j_SD,i_SD))  &
					             +(-2.d0*C_atm/rij_ep_old(0,j_SD,i)+C_atm*C_atm)*DEXP(-C_atm*rij_ep_old(0,j_SD,i))
					gsdee_up(1:3,j)=gsdee_up(1:3,j)+der1_up*ISDe_up_old(i,j)
					gsdee_dw(1:3,j)=gsdee_dw(1:3,j)+der1_dw*ISDe_dw_old(i,j)
					lsdee=lsdee+(der2_up*ISDe_up_old(i,j)+der2_dw*ISDe_dw_old(i,j))
				END DO
			END DO

		CASE ('1sb')

         norm=1.d0

         DO j = 1, N_part, 1
         DO i = 1, N_part, 1
            frf5=A_POT_se*(rij_ep_old(0,i,j)-D_POT_se)
            frf3=DEXP(frf5)
            sigm(i,j)=1.d0/(1.d0+frf3)
            sigm1(i,j)=-frf3/((1.d0+frf3)*(1.d0+frf3))
            sigm2(i,j)=(2.d0*frf3*frf3-frf3*(1.d0+frf3))/((1.d0+frf3)**3)
            Dd(1:3,i,j)=A_POT_se*rij_ep_old(1:3,i,j)/rij_ep_old(0,i,j)
            !D2d(1:3,i,j)=A_POT_se*(1.d0-rij_ep_old(1:3,i,j)*rij_ep_old(1:3,i,j)/(rij_ep_old(0,i,j)*rij_ep_old(0,i,j)))/&
            !   rij_ep_old(0,i,j)
         END DO
         END DO
         
         DO iadd = 0, H_N_part, H_N_part
            
            IF (iadd==0) THEN
               spin=1
               SDe=SDe_up_old
            ELSE IF (iadd==H_N_part) THEN
               spin=2
               SDe=SDe_dw_old
            ELSE
               STOP "ORRORE SPIN"
            END IF
            

            DO i = 1, H_N_part, 1
               q(1:3,i)=0.d0
               DO ik = 1, N_part, 1
                  q(1:3,i)=q(1:3,i)+rp_old(1:3,ik)*sigm(i+iadd,ik)
               END DO
            END DO
            CALL valuta_distanza_ij(re_old(1:3,1+iadd:H_N_part+iadd),q(1:3,1:H_N_part),&
               H_N_part,L(1:3),rq(0:3,1:H_N_part,1:H_N_part))

            DO il = 1, H_N_part, 1
            DO j = 1, H_N_part, 1
            DO i = 1, H_N_part, 1
            DO alpha = 1, 3, 1
               Grq(alpha,i,j,il,spin)=0.d0
               IF (il==i) Grq(alpha,i,j,il,spin)=Grq(alpha,i,j,il,spin)+rq(alpha,i,j)
               IF (j==il) THEN
                  DO ik = 1, N_part, 1
                  DO beta = 1, 3, 1
                     Grq(alpha,i,j,il,spin)=Grq(alpha,i,j,il,spin)-&
                        rq(beta,i,j)*rp_old(beta,ik)*sigm1(j+iadd,ik)*Dd(alpha,j+iadd,ik)
                  END DO
                  END DO
               END IF
               Grq(alpha,i,j,il,spin)=Grq(alpha,i,j,il,spin)/rq(0,i,j)
               Gphi(alpha,i,j,il,spin)=-C_atm*SDe(i,j)*Grq(alpha,i,j,il,spin)
            END DO
            END DO
            END DO
            END DO

            DO il = 1, H_N_part, 1
            DO j = 1, H_N_part, 1
            DO i = 1, H_N_part, 1
            DO alpha = 1, 3, 1
            DO beta = 1, 3, 1
               Gvrq(beta,alpha,i,j,il)=0.d0
               IF ((i==il).AND.(alpha==beta)) Gvrq(beta,alpha,i,j,il)=Gvrq(beta,alpha,i,j,il)+1.d0
               IF (il==j) THEN
                  DO ik = 1, N_part, 1
                     Gvrq(beta,alpha,i,j,il)=Gvrq(beta,alpha,i,j,il)-rp_old(beta,ik)*sigm1(j+iadd,ik)*Dd(alpha,j+iadd,ik)
                  END DO
               END IF
            END DO
            END DO
            END DO
            END DO
            END DO

            DO il = 1, H_N_part, 1
            DO j = 1, H_N_part, 1
            DO i = 1, H_N_part, 1
            DO alpha = 1, 3, 1
               Lrq(alpha,i,j,il,spin)=0.d0
               !1
               IF ((il==i).OR.(il==j)) Lrq(alpha,i,j,il,spin)=Lrq(alpha,i,j,il,spin)-&
                  Grq(alpha,i,j,il,spin)*Grq(alpha,i,j,il,spin)/rq(0,i,j)
               IF (i==il) THEN
                  !2
                  Lrq(alpha,i,j,il,spin)=Lrq(alpha,i,j,il,spin)+Gvrq(alpha,alpha,i,j,il)/rq(0,i,j)
               END IF
               DO ik = 1, N_part, 1
                  DO beta = 1, 3, 1
                     !3
                     Lrq(alpha,i,j,il,spin)=Lrq(alpha,i,j,il,spin)-Gvrq(beta,alpha,i,j,il)*rp_old(beta,ik)*&
                        sigm1(j+iadd,ik)*Dd(alpha,j+iadd,ik)/rq(0,i,j)
                  END DO
               END DO
               IF (j==il) THEN
                  DO ik = 1, N_part, 1
                     DO beta = 1, 3, 1
                        !4
                        Lrq(alpha,i,j,il,spin)=Lrq(alpha,i,j,il,spin)-rq(beta,i,j)*rp_old(beta,ik)*sigm2(j+iadd,ik)*&
                           Dd(alpha,j+iadd,ik)*Dd(alpha,j+iadd,ik)/rq(0,i,j)
                        !5
                        Lrq(alpha,i,j,il,spin)=Lrq(alpha,i,j,il,spin)-rq(beta,i,j)*rp_old(beta,ik)*sigm1(j+iadd,ik)*&
                           A_POT_se*(1.d0-(rij_ep_old(alpha,j+iadd,ik)/rij_ep_old(0,j+iadd,ik))**2)/&
                           (rq(0,i,j)*rij_ep_old(0,j+iadd,ik))
                     END DO
                  END DO
               END IF
               Lphi(alpha,i,j,il,spin)=C_atm*SDe(i,j)*(C_atm*Grq(alpha,i,j,il,spin)*Grq(alpha,i,j,il,spin)-Lrq(alpha,i,j,il,spin))
            END DO
            END DO
            END DO
            END DO

         END DO

         DO il = 1, H_N_part, 1
            DO j = 1, H_N_part, 1
            DO i = 1, H_N_part, 1
               gsdee_up(1:3,il)=gsdee_up(1:3,il)+Gphi(1:3,i,j,il,1)*ISDe_up_old(j,i)
            END DO
            END DO
         END DO
         DO il = 1, H_N_part, 1
            DO j = 1, H_N_part, 1
            DO i = 1, H_N_part, 1
               gsdee_dw(1:3,il)=gsdee_dw(1:3,il)+Gphi(1:3,i,j,il,2)*ISDe_dw_old(j,i)
            END DO
            END DO
         END DO


         DO il = 1, H_N_part, 1
            DO j = 1, H_N_part, 1
            DO i = 1, H_N_part, 1
            DO alpha = 1, 3, 1
               lsdee=lsdee+Lphi(alpha,i,j,il,1)*ISDe_up_old(j,i)
            END DO
            END DO
            END DO
            DO j = 1, H_N_part, 1
            DO i = 1, H_N_part, 1
            DO alpha = 1, 3, 1
               lsdee=lsdee+Lphi(alpha,i,j,il,2)*ISDe_dw_old(j,i)
            END DO
            END DO
            END DO
         END DO

         !!!!CHECK GRADIENTE
         !!!il=1
         !!!PRINT *, "CHECK GRADIENT "
         !!!PRINT *, "Analytical: ", REAL(gsdee_up(1:3,1))
         !!!PRINT *, " -> ", Grq(1:3,1,1,il,1)
         !!!PRINT *, 

         !!!frf6=0.00001d0

         !!!q(1:3,1)=re_old(1:3,il)
         !!!DO ik = 1, N_part, 1
         !!!   q(1:3,1)=q(1:3,1)-rp_old(1:3,ik)/(1.d0+DEXP(A_POT_se*(rij_ep_old(0,il,ik)-D_POT_se)))
         !!!END DO
         !!!q(1:3,1)=q(1:3,1)-L(1:3)*DNINT(q(1:3,1)/L(1:3))
         !!!q(0,1)=DSQRT(DOT_PRODUCT(q(1:3,1),q(1:3,1)))
			!!!frf3=q(0,1)
			!!!!frf3=norm*DEXP(-C_atm*q(0,1))
			!!!frf3=q(0,1)
			!!!
         !!!re_new(1:3,il)=re_old(1:3,il)+(/ frf6, 0.d0, 0.d0 /)
         !!!CALL valuta_distanza_ij(re_new,rp_old,N_part,L,rij_ep_new)
         !!!q(1:3,1)=re_new(1:3,il)
         !!!DO ik = 1, N_part, 1
         !!!   q(1:3,1)=q(1:3,1)-rp_new(1:3,ik)/(1.d0+DEXP(A_POT_se*(rij_ep_new(0,il,ik)-D_POT_se)))
         !!!END DO
         !!!q(1:3,1)=q(1:3,1)-L(1:3)*DNINT(q(1:3,1)/L(1:3))
         !!!q(0,1)=DSQRT(DOT_PRODUCT(q(1:3,1),q(1:3,1)))
			!!!!frf2(1)=norm*DEXP(-C_atm*q(0,1))
			!!!frf2(1)=q(0,1)

         !!!re_new(1:3,il)=re_old(1:3,il)+(/ 0.d0, frf6, 0.d0 /)
         !!!CALL valuta_distanza_ij(re_new,rp_old,N_part,L,rij_ep_new)
         !!!q(1:3,1)=re_new(1:3,il)
         !!!DO ik = 1, N_part, 1
         !!!   q(1:3,1)=q(1:3,1)-rp_new(1:3,ik)/(1.d0+DEXP(A_POT_se*(rij_ep_new(0,il,ik)-D_POT_se)))
         !!!END DO
         !!!q(1:3,1)=q(1:3,1)-L(1:3)*DNINT(q(1:3,1)/L(1:3))
         !!!q(0,1)=DSQRT(DOT_PRODUCT(q(1:3,1),q(1:3,1)))
			!!!!frf2(2)=norm*DEXP(-C_atm*q(0,1))
			!!!frf2(2)=q(0,1)

         !!!re_new(1:3,il)=re_old(1:3,il)+(/ 0.d0, 0.d0, frf6 /)
         !!!CALL valuta_distanza_ij(re_new,rp_old,N_part,L,rij_ep_new)
         !!!q(1:3,1)=re_new(1:3,il)
         !!!DO ik = 1, N_part, 1
         !!!   q(1:3,1)=q(1:3,1)-rp_new(1:3,ik)/(1.d0+DEXP(A_POT_se*(rij_ep_new(0,il,ik)-D_POT_se)))
         !!!END DO
         !!!q(1:3,1)=q(1:3,1)-L(1:3)*DNINT(q(1:3,1)/L(1:3))
         !!!q(0,1)=DSQRT(DOT_PRODUCT(q(1:3,1),q(1:3,1)))
			!!!!frf2(3)=norm*DEXP(-C_atm*q(0,1))
			!!!frf2(3)=q(0,1)

         !!!PRINT *, "Numerical:  ", ((DEXP(-C_atm*frf2(1:3))-DEXP(-C_atm*frf3))/frf6)/DEXP(-C_atm*frf3)
         !!!PRINT *, " -> ", (frf2(1:3)-frf3)/frf6
         !!!PRINT *, 
         !!!PRINT *, 
         !!!!FINE CHECK GRADIENTE

         !!!!CHECK LAPLACIANO
         !!!PRINT *, "CHECK LAPLACIAN"

         !!!il=1
         !!!PRINT *, "Analytical: ", Lphi(1:3,1,1,il,1)/REAL(SDe_up_old(1,1))
         !!!PRINT *, " -> ", Lrq(1:3,1,1,il,1)
         !!!PRINT *, 

         !!!re_new(1:3,il)=re_old(1:3,il)-(/ frf6, 0.d0, 0.d0 /)
         !!!CALL valuta_distanza_ij(re_new,rp_old,N_part,L,rij_ep_new)
         !!!q(1:3,1)=re_new(1:3,il)
         !!!DO ik = 1, N_part, 1
         !!!   q(1:3,1)=q(1:3,1)-rp_new(1:3,ik)/(1.d0+DEXP(A_POT_se*(rij_ep_new(0,il,ik)-D_POT_se)))
         !!!END DO
         !!!q(1:3,1)=q(1:3,1)-L(1:3)*DNINT(q(1:3,1)/L(1:3))
         !!!q(0,1)=DSQRT(DOT_PRODUCT(q(1:3,1),q(1:3,1)))
			!!!!frf1(1)=norm*DEXP(-C_atm*q(0,1))
			!!!frf1(1)=q(0,1)

         !!!re_new(1:3,il)=re_old(1:3,il)-(/ 0.d0, frf6, 0.d0 /)
         !!!CALL valuta_distanza_ij(re_new,rp_old,N_part,L,rij_ep_new)
         !!!q(1:3,1)=re_new(1:3,il)
         !!!DO ik = 1, N_part, 1
         !!!   q(1:3,1)=q(1:3,1)-rp_new(1:3,ik)/(1.d0+DEXP(A_POT_se*(rij_ep_new(0,il,ik)-D_POT_se)))
         !!!END DO
         !!!q(1:3,1)=q(1:3,1)-L(1:3)*DNINT(q(1:3,1)/L(1:3))
         !!!q(0,1)=DSQRT(DOT_PRODUCT(q(1:3,1),q(1:3,1)))
			!!!!frf1(2)=norm*DEXP(-C_atm*q(0,1))
			!!!frf1(2)=q(0,1)

         !!!re_new(1:3,il)=re_old(1:3,il)-(/ 0.d0, 0.d0, frf6 /)
         !!!CALL valuta_distanza_ij(re_new,rp_old,N_part,L,rij_ep_new)
         !!!q(1:3,1)=re_new(1:3,il)
         !!!DO ik = 1, N_part, 1
         !!!   q(1:3,1)=q(1:3,1)-rp_new(1:3,ik)/(1.d0+DEXP(A_POT_se*(rij_ep_new(0,il,ik)-D_POT_se)))
         !!!END DO
         !!!q(1:3,1)=q(1:3,1)-L(1:3)*DNINT(q(1:3,1)/L(1:3))
         !!!q(0,1)=DSQRT(DOT_PRODUCT(q(1:3,1),q(1:3,1)))
			!!!!frf1(3)=norm*DEXP(-C_atm*q(0,1))
			!!!frf1(3)=q(0,1)

         !!!PRINT *, "Numerical:  ", (DEXP(-C_atm*frf2(1:3))-2.d0*DEXP(-C_atm*frf3)+DEXP(-C_atm*frf1(1:3)))/&
         !!!   (frf6*frf6*DEXP(-C_atm*frf3))
         !!!PRINT *, " -> ", (frf2(1:3)-2.d0*frf3+frf1(1:3))/(frf6*frf6)
         !!!PRINT *, 
         !!!!FINE CHECK LAPLACIANO
         !!!STOP


		CASE ('atp')
			frfs1(1:3)=PI/L(1:3)
			DO i = 1, H_N_part, 1
				i_SD=i+H_N_part
				DO j = 1, H_N_part, 1
					j_SD=j+H_N_part
					
					der1_up(1:3)=-DCOS(frfs1(1:3)*rij_ep_old(1:3,j,i))* &
					  rijpc_ep_old(1:3,j,i)*C_atm/rijpc_ep_old(0,j,i)
					der1_dw(1:3)=-DCOS(frfs1(1:3)*rij_ep_old(1:3,j_SD,i_SD))* &
					  rijpc_ep_old(1:3,j_SD,i_SD)*C_atm/rijpc_ep_old(0,j_SD,i_SD)
					der2_up=DOT_PRODUCT(der1_up(1:3),der1_up(1:3))
					der2_dw=DOT_PRODUCT(der1_dw(1:3),der1_dw(1:3))
					
					DO i3= 1, 3, 1
						frf1(i3)=((DCOS(frfs1(i3)*rij_ep_old(i3,j,i)))**2) * &
						  (1.d0-(rijpc_ep_old(i3,j,i)/rijpc_ep_old(0,j,i))**2)/rijpc_ep_old(0,j,i)
						frf1(i3)=frf1(i3)-frfs1(i3)*DSIN(frfs1(i3)*rij_ep_old(i3,j,i))* &
						  rijpc_ep_old(i3,j,i)/rijpc_ep_old(0,j,i)
						frf2(i3)=((DCOS(frfs1(i3)*rij_ep_old(i3,j_SD,i_SD)))**2) * &
						  (1.d0-(rijpc_ep_old(i3,j_SD,i_SD)/rijpc_ep_old(0,j_SD,i_SD))**2)/rijpc_ep_old(0,j_SD,i_SD)
						frf2(i3)=frf1(i3)-frfs1(i3)*DSIN(frfs1(i3)*rij_ep_old(i3,j_SD,i_SD))* &
						  rijpc_ep_old(i3,j_SD,i_SD)/rijpc_ep_old(0,j_SD,i_SD)
						der2_up=der2_up-C_atm*frf1(i3)
						der2_dw=der2_dw-C_atm*frf2(i3)
					END DO
					
					gsdee_up(1:3,j)=gsdee_up(1:3,j)+SDe_up_old(j,i)*der1_up*ISDe_up_old(i,j)
					gsdee_dw(1:3,j)=gsdee_dw(1:3,j)+SDe_dw_old(j,i)*der1_dw*ISDe_dw_old(i,j)
					lsdee=lsdee+(SDe_up_old(j,i)*der2_up*ISDe_up_old(i,j)+SDe_dw_old(j,i)*der2_dw*ISDe_dw_old(i,j))
					
				END DO
			END DO
		END SELECT
		
		!calcolo la parte del Jastrow ee
		lee=0.d0
		gee=0.d0
		SELECT CASE (Jee_kind)
		CASE ('yuk')
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
					IF (i/=j) THEN
						frf6=1.d0/rij_ee_old(0,i,j)
						IF (split_Aee) THEN
							IF (((i<=H_N_part) .AND. (j<=H_N_part)) .OR. ((i>H_N_part) .AND. (j>H_N_part))) THEN
								frf1(0)=Aee_yuk/(rij_ee_old(0,i,j)*rij_ee_old(0,i,j))
							ELSE
								frf1(0)=Aee_ud_yuk/(rij_ee_old(0,i,j)*rij_ee_old(0,i,j))
							END IF
						ELSE
							frf1(0)=Aee_yuk/(rij_ee_old(0,i,j)*rij_ee_old(0,i,j))
						END IF
						IF (split_Fee) THEN
							IF (((i<=H_N_part) .AND. (j<=H_N_part)) .OR. ((i>H_N_part) .AND. (j>H_N_part))) THEN
								frf5=Fee_yuk*rij_ee_old(0,i,j)
							ELSE
								frf5=Fee_ud_yuk*rij_ee_old(0,i,j)
							END IF
						ELSE
							frf5=Fee_yuk*rij_ee_old(0,i,j)
						END IF
						frf3=frf1(0)*DEXP(-frf5)
						uee1=frf3*(frf5+1.d0)-frf1(0)
						uee2=2.d0*frf1(0)*frf6-frf3*frf6*(frf5*frf5+2.d0*frf5+2.d0)
						frf1(0)=uee1*frf6
						gee(1:3,j)=gee(1:3,j)-0.5d0*rij_ee_old(1:3,i,j)*frf1(0)
						IF (i>j) THEN
							lee=lee-uee2-2.d0*frf1(0)
						END IF
					END IF
				END DO
				lee=lee+DOT_PRODUCT(gee(1:3,j),gee(1:3,j))
			END DO
		CASE ('yup')
			frf1(1:3)=PI/L(1:3)
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
					IF (i/=j) THEN
						frf2(1:3)=(/ rijpc_ee_old(1,i,j)*DCOS(frf1(1)*rij_ee_old(1,i,j)), &
						  rijpc_ee_old(2,i,j)*DCOS(frf1(2)*rij_ee_old(2,i,j)), &
						  rijpc_ee_old(3,i,j)*DCOS(frf1(3)*rij_ee_old(3,i,j)) /)/rijpc_ee_old(0,i,j)
						IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
							frf3=DEXP(-Fee_yuk*rijpc_ee_old(0,i,j))
							uee1=Aee_yuk*(-(1.d0-frf3)/rijpc_ee_old(0,i,j)+Fee_yuk*frf3)/rijpc_ee_old(0,i,j)
							uee2=Aee_yuk*(2.d0*(1.d0-frf3)/(rijpc_ee_old(0,i,j)**2)-2.d0*Fee_yuk* &
							  frf3/rijpc_ee_old(0,i,j)-Fee_yuk*Fee_yuk*frf3)/rijpc_ee_old(0,i,j)
						ELSE
							IF (split_Fee) THEN
								IF (split_Aee) THEN
									frf3=DEXP(-Fee_ud_yuk*rijpc_ee_old(0,i,j))
									uee1=Aee_ud_yuk*(-(1.d0-frf3)/rijpc_ee_old(0,i,j)+Fee_ud_yuk*frf3)/rijpc_ee_old(0,i,j)
									uee2=Aee_ud_yuk*(2.d0*(1.d0-frf3)/(rijpc_ee_old(0,i,j)**2)-2.d0*Fee_ud_yuk* &
									  frf3/rijpc_ee_old(0,i,j)-Fee_ud_yuk*Fee_ud_yuk*frf3)/rijpc_ee_old(0,i,j)
								ELSE
									frf3=DEXP(-Fee_ud_yuk*rijpc_ee_old(0,i,j))
									uee1=Aee_yuk*(-(1.d0-frf3)/rijpc_ee_old(0,i,j)+Fee_ud_yuk*frf3)/rijpc_ee_old(0,i,j)
									uee2=Aee_yuk*(2.d0*(1.d0-frf3)/(rijpc_ee_old(0,i,j)**2)-2.d0*Fee_ud_yuk* &
									  frf3/rijpc_ee_old(0,i,j)-Fee_ud_yuk*Fee_ud_yuk*frf3)/rijpc_ee_old(0,i,j)
								END IF
							ELSE
								IF (split_Aee) THEN
									frf3=DEXP(-Fee_yuk*rijpc_ee_old(0,i,j))
									uee1=Aee_ud_yuk*(-(1.d0-frf3)/rijpc_ee_old(0,i,j)+Fee_yuk*frf3)/rijpc_ee_old(0,i,j)
									uee2=Aee_ud_yuk*(2.d0*(1.d0-frf3)/(rijpc_ee_old(0,i,j)**2)-2.d0*Fee_yuk* &
									  frf3/rijpc_ee_old(0,i,j)-Fee_yuk*Fee_yuk*frf3)/rijpc_ee_old(0,i,j)
								ELSE
									frf3=DEXP(-Fee_yuk*rijpc_ee_old(0,i,j))
									uee1=Aee_yuk*(-(1.d0-frf3)/rijpc_ee_old(0,i,j)+Fee_yuk*frf3)/rijpc_ee_old(0,i,j)
									uee2=Aee_yuk*(2.d0*(1.d0-frf3)/(rijpc_ee_old(0,i,j)**2)-2.d0*Fee_yuk* &
									  frf3/rijpc_ee_old(0,i,j)-Fee_yuk*Fee_yuk*frf3)/rijpc_ee_old(0,i,j)
								END IF
							END IF
						END IF
						gee(1:3,j)=gee(1:3,j)-0.5d0*uee1*frf2(1:3)
						IF (i>j) THEN
							lee=lee-( uee2*DOT_PRODUCT(frf2(1:3),frf2(1:3))+uee1*( -((rijpc_ee_old(1,i,j)**2-rijpc_ee_old(0,i,j)**2)* &
							  DCOS(frf1(1)*rij_ee_old(1,i,j))**2+(rijpc_ee_old(2,i,j)**2-rijpc_ee_old(0,i,j)**2)*DCOS(frf1(2)*rij_ee_old(2,i,j))**2+&
							  (rijpc_ee_old(3,i,j)**2-rijpc_ee_old(0,i,j)**2)*DCOS(frf1(3)*rij_ee_old(3,i,j))**2)/(rijpc_ee_old(0,i,j)**3) -&
							  (frf1(1)*rijpc_ee_old(1,i,j)*DSIN(frf1(1)*rij_ee_old(1,i,j))+frf1(2)*rijpc_ee_old(2,i,j)*DSIN(frf1(2)*rij_ee_old(2,i,j))+&
							  frf1(3)*rijpc_ee_old(3,i,j)*DSIN(frf1(3)*rij_ee_old(3,i,j)))/rijpc_ee_old(0,i,j) ) )
						END IF
					END IF
				END DO
				lee=lee+DOT_PRODUCT(gee(1:3,j),gee(1:3,j))
			END DO
		CASE ('spl')
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
               IF (i/=j) THEN
                  IF (split_Aee.OR.split_Fee) THEN
                     IF (((i<=H_N_part).AND.(j<=H_N_part)).OR.((i>H_N_part).AND.(j>H_N_part))) THEN
                        CALL MSPL_compute(SPL=Jsplee, DERIV=1, R=rij_ee_old(0,i,j), VAL=frf3)
                        CALL MSPL_compute(SPL=Jsplee, DERIV=2, R=rij_ee_old(0,i,j), VAL=frf5)
                     ELSE
                        CALL MSPL_compute(SPL=Jsplee_ud, DERIV=1, R=rij_ee_old(0,i,j), VAL=frf3)
                        CALL MSPL_compute(SPL=Jsplee_ud, DERIV=2, R=rij_ee_old(0,i,j), VAL=frf5)
                     END IF
                  ELSE
                     CALL MSPL_compute(SPL=Jsplee, DERIV=1, R=rij_ee_old(0,i,j), VAL=frf3)
                     CALL MSPL_compute(SPL=Jsplee, DERIV=2, R=rij_ee_old(0,i,j), VAL=frf5)
                  END IF
                  frf1(1:3)=rij_ee_old(1:3,i,j)/rij_ee_old(0,i,j)
                  gee(1:3,j)=gee(1:3,j)-0.5d0*frf3*frf1(1:3)

                  lee=lee-0.5d0*( 2.d0*frf3/rij_ee_old(0,i,j) + frf5 )
					END IF
				END DO
				lee=lee+DOT_PRODUCT(gee(1:3,j),gee(1:3,j))

            !!!!CHECK GRADIENTE (ricorda di rimettere PROTECTED nei walkers)
            !!!frf6=0.0001d0
            !!!frf3=0.d0
            !!!frf2=0.d0
            !!!frf1=0.d0

            !!!DO i = 1, N_part, 1
            !!!IF (i/=j) THEN
            !!!   CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij_ee_old(0,i,j), VAL=frf3, RESET=.FALSE.)
            !!!   re_new=re_old

            !!!   re_new(1:3,j)=re_old(1:3,j)+(/ frf6, 0.d0, 0.d0/)
            !!!   CALL valuta_distanza_ii(re_new,N_part,L,rij_ee_new)
            !!!   CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij_ee_new(0,i,j), VAL=frf2(1), RESET=.FALSE.)
            !!!   re_new(1:3,j)=re_old(1:3,j)

            !!!   re_new(1:3,j)=re_old(1:3,j)+(/ 0.d0, frf6, 0.d0/)
            !!!   CALL valuta_distanza_ii(re_new,N_part,L,rij_ee_new)
            !!!   CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij_ee_new(0,i,j), VAL=frf2(2), RESET=.FALSE.)
            !!!   re_new(1:3,j)=re_old(1:3,j)

            !!!   re_new(1:3,j)=re_old(1:3,j)+(/ 0.d0, 0.d0, frf6/)
            !!!   CALL valuta_distanza_ii(re_new,N_part,L,rij_ee_new)
            !!!   CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij_ee_new(0,i,j), VAL=frf2(3), RESET=.FALSE.)
            !!!   re_new(1:3,j)=re_old(1:3,j)
            !!!END IF
            !!!END DO

            !!!PRINT *, "Check gradiente"
            !!!PRINT *, "Analitico : ", gee(1:3,j)
            !!!PRINT *, "Numerico : ", & 
            !!!(DEXP(-0.5d0*frf2(1:3))-DEXP(-0.5d0*frf3))/(DEXP(-0.5d0*frf3)*frf6)
            !!!PRINT *, 
            !!!!FINE CHECK GRADIENTE
            !!!
            !!!!CHECK LAPLACIANO
            !!!DO i = 1, N_part, 1
            !!!IF (i/=j) THEN
            !!!   re_new(1:3,j)=re_old(1:3,j)-(/ frf6, 0.d0, 0.d0/)
            !!!   CALL valuta_distanza_ii(re_new,N_part,L,rij_ee_new)
            !!!   CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij_ee_new(0,i,j), VAL=frf1(1), RESET=.FALSE.)
            !!!   re_new(1:3,j)=re_old(1:3,j)

            !!!   re_new(1:3,j)=re_old(1:3,j)-(/ 0.d0, frf6, 0.d0/)
            !!!   CALL valuta_distanza_ii(re_new,N_part,L,rij_ee_new)
            !!!   CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij_ee_new(0,i,j), VAL=frf1(2), RESET=.FALSE.)
            !!!   re_new(1:3,j)=re_old(1:3,j)

            !!!   re_new(1:3,j)=re_old(1:3,j)-(/ 0.d0, 0.d0, frf6/)
            !!!   CALL valuta_distanza_ii(re_new,N_part,L,rij_ee_new)
            !!!   CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij_ee_new(0,i,j), VAL=frf1(3), RESET=.FALSE.)
            !!!   re_new(1:3,j)=re_old(1:3,j)
            !!!END IF
            !!!END DO

            !!!PRINT *, "Check laplaciano"
            !!!PRINT *, "Analitico : ", lee
            !!!PRINT *, "Numerico : ", &
            !!!   SUM(DEXP(-0.5d0*frf2(1:3))-2.d0*DEXP(-0.5d0*frf3)+DEXP(-0.5d0*frf1(1:3)))/(DEXP(-0.5d0*frf3)*frf6*frf6)

            !!!!FINE CHECK LAPLACIANO
            !!!STOP

			END DO
		CASE ('spp')
         frf2(1:3)=PI/L(1:3)
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
               IF (i/=j) THEN
                  IF (split_Aee.OR.split_Fee) THEN
                     IF (((i<=H_N_part).AND.(j<=H_N_part)).OR.((i>H_N_part).AND.(j>H_N_part))) THEN
                        CALL MSPL_compute(SPL=Jsplee, DERIV=1, R=rijpc_ee_old(0,i,j), VAL=frf3)
                        CALL MSPL_compute(SPL=Jsplee, DERIV=2, R=rijpc_ee_old(0,i,j), VAL=frf5)
                     ELSE
                        CALL MSPL_compute(SPL=Jsplee_ud, DERIV=1, R=rijpc_ee_old(0,i,j), VAL=frf3)
                        CALL MSPL_compute(SPL=Jsplee_ud, DERIV=2, R=rijpc_ee_old(0,i,j), VAL=frf5)
                     END IF
                  ELSE
                     CALL MSPL_compute(SPL=Jsplee, DERIV=1, R=rijpc_ee_old(0,i,j), VAL=frf3)
                     CALL MSPL_compute(SPL=Jsplee, DERIV=2, R=rijpc_ee_old(0,i,j), VAL=frf5)
                  END IF
                  frf1(1:3)=DCOS(frf2(1:3)*rij_ee_old(1:3,i,j))*rijpc_ee_old(1:3,i,j)/rijpc_ee_old(0,i,j)
                  gee(1:3,j)=gee(1:3,j)-0.5d0*frf3*frf1(1:3)

                  lee=lee-0.5d0*( frf3*(   &
                     DOT_PRODUCT(DCOS(rij_ee_old(1:3,i,j)*frf2(1:3)),DCOS(rij_ee_old(1:3,i,j)*frf2(1:3)))/rijpc_ee_old(0,i,j) &
                        - DOT_PRODUCT(rijpc_ee_old(1:3,i,j)*DCOS(rij_ee_old(1:3,i,j)*frf2(1:3)), &
                           rijpc_ee_old(1:3,i,j)*DCOS(rij_ee_old(1:3,i,j)*frf2(1:3)))/(rijpc_ee_old(0,i,j)**3) &
                        - DOT_PRODUCT(rijpc_ee_old(1:3,i,j)*frf2(1:3),rijpc_ee_old(1:3,i,j)*frf2(1:3))/rijpc_ee_old(0,i,j) &
                     ) + frf5*DOT_PRODUCT(frf1(1:3),frf1(1:3)) )
					END IF
				END DO
				lee=lee+DOT_PRODUCT(gee(1:3,j),gee(1:3,j))
         END DO
		END SELECT
		
		!calcolo la parte del Jastrow ep
		lep=0.d0
		gep=0.d0
		SELECT CASE (Jep_kind)
		CASE ('yuk')
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
					IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
						uep1=Aep_yuk*( -(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_yuk* &
						  DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
						uep2=Aep_yuk*( 2.d0*(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/(rij_ep_old(0,i,j)**2)+ &
						  (-2.d0*Fep_yuk/rij_ep_old(0,i,j)-Fep_yuk*Fep_yuk)*DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
					ELSE
						IF (split_Fee) THEN
							IF (split_Aee) THEN
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_ud_yuk* &
								  DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
								uep2=Aep_ud_yuk*( 2.d0*(1.d0-DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)))/(rij_ep_old(0,i,j)**2)+ &
								  (-2.d0*Fep_ud_yuk/rij_ep_old(0,i,j)-Fep_ud_yuk*Fep_ud_yuk)*DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
							ELSE
								uep1=Aep_yuk*( -(1.d0-DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_ud_yuk* &
								  DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
								uep2=Aep_yuk*( 2.d0*(1.d0-DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)))/(rij_ep_old(0,i,j)**2)+ &
								  (-2.d0*Fep_ud_yuk/rij_ep_old(0,i,j)-Fep_ud_yuk*Fep_ud_yuk)*DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
							END IF
						ELSE
							IF (split_Aee) THEN
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_yuk* &
								  DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
								uep2=Aep_ud_yuk*( 2.d0*(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/(rij_ep_old(0,i,j)**2)+ &
								  (-2.d0*Fep_yuk/rij_ep_old(0,i,j)-Fep_yuk*Fep_yuk)*DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
							ELSE
								uep1=Aep_yuk*( -(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_yuk* &
								  DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
								uep2=Aep_yuk*( 2.d0*(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/(rij_ep_old(0,i,j)**2)+ &
								  (-2.d0*Fep_yuk/rij_ep_old(0,i,j)-Fep_yuk*Fep_yuk)*DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
							END IF
						END IF
					END IF
					gep(1:3,i)=gep(1:3,i)-0.5d0*uep1*rij_ep_old(1:3,i,j)/rij_ep_old(0,i,j)
					lep=lep-0.5d0*uep2-uep1/rij_ep_old(0,i,j)
				END DO
			END DO
			DO j = 1, N_part, 1
				lep=lep+DOT_PRODUCT(gep(1:3,j),gep(1:3,j))
			END DO
		CASE ('yup')
			frf1(1:3)=PI/L(1:3)
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
					frf2(1:3)=(/ rijpc_ep_old(1,i,j)*DCOS(frf1(1)*rij_ep_old(1,i,j)), &
					  rijpc_ep_old(2,i,j)*DCOS(frf1(2)*rij_ep_old(2,i,j)), &
					  rijpc_ep_old(3,i,j)*DCOS(frf1(3)*rij_ep_old(3,i,j)) /)/rijpc_ep_old(0,i,j)
					
					IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
						uep1=Aep_yuk*( -(1.d0-DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)))/rijpc_ep_old(0,i,j)+Fep_yuk* &
						  DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
						uep2=Aep_yuk*( 2.d0*(1.d0-DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)))/(rijpc_ep_old(0,i,j)**2)+ &
						  (-2.d0*Fep_yuk/rijpc_ep_old(0,i,j)-Fep_yuk*Fep_yuk)*DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
					ELSE
						IF (split_Fee) THEN
							IF (split_Aee) THEN
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)))/rijpc_ep_old(0,i,j)+Fep_ud_yuk* &
								  DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
								uep2=Aep_ud_yuk*( 2.d0*(1.d0-DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)))/(rijpc_ep_old(0,i,j)**2)+ &
								  (-2.d0*Fep_ud_yuk/rijpc_ep_old(0,i,j)-Fep_ud_yuk*Fep_ud_yuk)*DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
							ELSE
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)))/rijpc_ep_old(0,i,j)+Fep_ud_yuk* &
								  DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
								uep2=Aep_ud_yuk*( 2.d0*(1.d0-DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)))/(rijpc_ep_old(0,i,j)**2)+ &
								  (-2.d0*Fep_ud_yuk/rijpc_ep_old(0,i,j)-Fep_ud_yuk*Fep_ud_yuk)*DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
							END IF
						ELSE
							IF (split_Aee) THEN
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)))/rijpc_ep_old(0,i,j)+Fep_yuk* &
								  DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
								uep2=Aep_ud_yuk*( 2.d0*(1.d0-DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)))/(rijpc_ep_old(0,i,j)**2)+ &
								  (-2.d0*Fep_yuk/rijpc_ep_old(0,i,j)-Fep_yuk*Fep_yuk)*DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
							ELSE
								frf3=DEXP(-Fep_yuk*rijpc_ep_old(0,i,j))
								uep1=Aep_yuk*( -(1.d0-frf3)/rijpc_ep_old(0,i,j)+Fep_yuk*frf3)/rijpc_ep_old(0,i,j)
								uep2=Aep_yuk*( 2.d0*(1.d0-frf3)/(rijpc_ep_old(0,i,j)**2)+ &
								  (-2.d0*Fep_yuk/rijpc_ep_old(0,i,j)-Fep_yuk*Fep_yuk)*frf3 )/rijpc_ep_old(0,i,j)
							END IF
						END IF
					END IF
				
					gep(1:3,i)=gep(1:3,i)-0.5d0*uep1*frf2(1:3)
					lep=lep-0.5d0*( uep2*DOT_PRODUCT(frf2(1:3),frf2(1:3))+uep1*( -((rijpc_ep_old(1,i,j)**2-rijpc_ep_old(0,i,j)**2)* &
					  DCOS(frf1(1)*rij_ep_old(1,i,j))**2+(rijpc_ep_old(2,i,j)**2-rijpc_ep_old(0,i,j)**2)*DCOS(frf1(2)*rij_ep_old(2,i,j))**2+&
					  (rijpc_ep_old(3,i,j)**2-rijpc_ep_old(0,i,j)**2)*DCOS(frf1(3)*rij_ep_old(3,i,j))**2)/(rijpc_ep_old(0,i,j)**3) -&
					  (frf1(1)*rijpc_ep_old(1,i,j)*DSIN(frf1(1)*rij_ep_old(1,i,j))+frf1(2)*rijpc_ep_old(2,i,j)*DSIN(frf1(2)*rij_ep_old(2,i,j))+&
					  frf1(3)*rijpc_ep_old(3,i,j)*DSIN(frf1(3)*rij_ep_old(3,i,j)))/rijpc_ep_old(0,i,j) ) )
				END DO
			END DO
			DO j = 1, N_part, 1
				lep=lep+DOT_PRODUCT(gep(1:3,j),gep(1:3,j))
			END DO
		CASE ('atm')
			DO i = 1, N_part, 1
				uep1=-0.5d0*Fep_yuk/rij_ep_old(0,i,i)
				uep2=2.d0*uep1
				gep(1:3,i)=gep(1:3,i)+uep1*rij_ep_old(1:3,i,i)
				lep=lep+uep2
			END DO
			DO j = 1, N_part, 1
				lep=lep+DOT_PRODUCT(gep(1:3,j),gep(1:3,j))
			END DO
		CASE ('atp')
			DO i = 1, N_part, 1
				uep1=-0.5d0*Fep_yuk/rijpc_ep_old(0,i,i)
				uep2=2.d0*uep1
				gep(1:3,i)=gep(1:3,i)+uep1*rijpc_ep_old(1:3,i,i)
				lep=lep+uep2
			END DO
			DO j = 1, N_part, 1
				lep=lep+DOT_PRODUCT(gep(1:3,j),gep(1:3,j))
			END DO
		END SELECT
		
		!calcolo la parte del kernel per le shadow-e
		gkernse1=0.d0
		gkernse2=0.d0
		lkernse=0.d0
		SELECT CASE (Kse_kind)
		CASE ('gss')
			DO i = 1, N_part, 1
				gkernse1(1:3,i)=-2.d0*C_kern_e*rij_ese1_old(1:3,i,i)
				gkernse2(1:3,i)=-2.d0*C_kern_e*rij_ese2_old(1:3,i,i)
				lkernse=lkernse+0.5d0*(DOT_PRODUCT(gkernse1(1:3,i),gkernse1(1:3,i))+DOT_PRODUCT(gkernse2(1:3,i),gkernse2(1:3,i)))
			END DO
			lkernse=lkernse-6.d0*C_kern_e*REAL(N_part,8)
			!ALTERNATIVA PIÙ LEGGIBILE
			!!!DO i = 1, N_part, 1
			!!!	
			!!!	frf2s1(1:3)=rij_ese1_old(1:3,i,i)/rij_ese1_old(0,i,i)       !dr/dx
			!!!	frf2s2(1:3)=rij_ese2_old(1:3,i,i)/rij_ese2_old(0,i,i)       !dr/dx
			!!!	frf3s1(1:3)=(1.d0-(rij_ese1_old(1:3,i,i)/rij_ese1_old(0,i,i))**2)/rij_ese1_old(0,i,i)   !d/dx dr/dx
			!!!	frf3s2(1:3)=(1.d0-(rij_ese2_old(1:3,i,i)/rij_ese2_old(0,i,i))**2)/rij_ese2_old(0,i,i)   !d/dx dr/dx
			!!!	
			!!!	gkernse1(1:3,i)=frf2s1(1:3)*(-2.d0*C_kern_e*rij_ese1_old(0,i,i))
			!!!	gkernse2(1:3,i)=frf2s2(1:3)*(-2.d0*C_kern_e*rij_ese2_old(0,i,i))
			!!!	
			!!!	DO i3 = 1, 3, 1
			!!!		lkernse=lkernse+frf3s1(i3)*(-2.d0*C_kern_e*rij_ese1_old(0,i,i))      !ok
			!!!		lkernse=lkernse+frf3s2(i3)*(-2.d0*C_kern_e*rij_ese2_old(0,i,i))      !ok
			!!!		lkernse=lkernse+frf2s1(i3)**2 * (-2.d0*C_kern_e)
			!!!		lkernse=lkernse+frf2s2(i3)**2 * (-2.d0*C_kern_e)
			!!!	END DO
			!!!	
			!!!	lkernse=lkernse+(DOT_PRODUCT(gkernse1(1:3,i),gkernse1(1:3,i))+DOT_PRODUCT(gkernse2(1:3,i),gkernse2(1:3,i)))
			!!!END DO
			!!!lkernse=lkernse*0.5d0
		CASE ('gsc')
			minL=MINVAL(L)
			DO i = 1, N_part, 1
				IF ( rij_ese1_old(0,i,i)<l1_kern_e ) THEN
					gkernse1(1:3,i)=-2.d0*C_kern_e*rij_ese1_old(1:3,i,i)
					lkernse=lkernse-3.d0*C_kern_e
				ELSE IF (rij_ese1_old(0,i,i)<minL*0.499) THEN
					gkernse1(1:3,i)=alpha1_kern_e*rij_ese1_old(1:3,i,i)/( rij_ese1_old(0,i,i)*((rij_ese1_old(0,i,i)-0.5d0*minL)**2) )
					lkernse=lkernse+(alpha1_kern_e/((rij_ese1_old(0,i,i)-0.5d0*minL)**2))* &
					  (1.d0/rij_ese1_old(0,i,i)-1.d0/(rij_ese1_old(0,i,i)-0.5d0*minL))
				ELSE
					gkernse1(1:3,i)=0.d0
				END IF
				IF ( rij_ese2_old(0,i,i)<l1_kern_e ) THEN
					gkernse2(1:3,i)=-2.d0*C_kern_e*rij_ese2_old(1:3,i,i)
					lkernse=lkernse-3.d0*C_kern_e
				ELSE IF (rij_ese2_old(0,i,i)<minL*0.499) THEN
					gkernse2(1:3,i)=alpha1_kern_e*rij_ese2_old(1:3,i,i)/( rij_ese2_old(0,i,i)*((rij_ese2_old(0,i,i)-0.5d0*minL)**2) )
					lkernse=lkernse+(alpha1_kern_e/((rij_ese2_old(0,i,i)-0.5d0*minL)**2))* &
					  (1.d0/rij_ese2_old(0,i,i)-1.d0/(rij_ese2_old(0,i,i)-0.5d0*minL))
				ELSE
					gkernse2(1:3,i)=0.d0
				END IF
				lkernse=lkernse+0.5d0*(DOT_PRODUCT(gkernse1(1:3,i),gkernse1(1:3,i))+DOT_PRODUCT(gkernse2(1:3,i),gkernse2(1:3,i)))
			END DO
		CASE ('gsp')
			frf2(1:3)=PI/L(1:3)
			frf6=-2.d0*C_kern_e
			DO i = 1, N_part, 1
				frfs1(1:3)=DCOS(rij_ese1_old(1:3,i,i)*frf2(1:3))
				frfs2(1:3)=DCOS(rij_ese2_old(1:3,i,i)*frf2(1:3))
				frf2s1(1:3)=frfs1(1:3)*rij_ese1_old(1:3,i,i)/rij_ese1_old(0,i,i)       !dr/dx
				frf2s2(1:3)=frfs2(1:3)*rij_ese2_old(1:3,i,i)/rij_ese2_old(0,i,i)       !dr/dx
				frf3s1(1:3)=-rijpc_ese1_old(1:3,i,i)*(frf2(1:3)**2)*rijpc_ese1_old(1:3,i,i)/rijpc_ese1_old(0,i,i) + &     !d/dx dr'/dx
				  (frfs1(1:3)**2) * (1.d0 - (rijpc_ese1_old(1:3,i,i)/rijpc_ese1_old(0,i,i))**2 )/rijpc_ese1_old(0,i,i)
				frf3s2(1:3)=-rijpc_ese2_old(1:3,i,i)*(frf2(1:3)**2)*rijpc_ese2_old(1:3,i,i)/rijpc_ese2_old(0,i,i) + &     !d/dx dr'/dx
				  (frfs2(1:3)**2) * (1.d0 - (rijpc_ese2_old(1:3,i,i)/rijpc_ese2_old(0,i,i))**2 )/rijpc_ese2_old(0,i,i)
				
				frfs1(1:3)=-2.d0*C_kern_e*rijpc_ese1_old(0,i,i)
				frfs2(1:3)=-2.d0*C_kern_e*rijpc_ese2_old(0,i,i)
				
				gkernse1(1:3,i)=frf2s1(1:3)*frfs1(1:3)
				gkernse2(1:3,i)=frf2s2(1:3)*frfs2(1:3)
				lkernse=lkernse+(DOT_PRODUCT(gkernse1(1:3,i),gkernse1(1:3,i))+DOT_PRODUCT(gkernse2(1:3,i),gkernse2(1:3,i)))
				DO i3 = 1, 3, 1
					lkernse=lkernse+frf3s1(i3)*frfs1(i3)      !ok
					lkernse=lkernse+frf3s2(i3)*frfs2(i3)      !ok
					lkernse=lkernse+(frf2s1(i3)*frf2s1(i3))*frf6
					lkernse=lkernse+(frf2s2(i3)*frf2s2(i3))*frf6
				END DO
			END DO
			lkernse=lkernse*0.5d0
		CASE ('gsd')
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					lkernse=lkernse+((2.d0*C_kern_e*rij_ese1_old(0,j,i))**2 - 6.d0*C_kern_e)*GDse1_up_old(j,i)*IGDse1_up_old(i,j)
					gkernse1(1:3,j)=gkernse1(1:3,j)-2.d0*C_kern_e*rij_ese1_old(1:3,j,i)*GDse1_up_old(j,i)*IGDse1_up_old(i,j)
				END DO
			END DO
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					lkernse=lkernse+((2.d0*C_kern_e*rij_ese1_old(0,j+H_N_part,i+H_N_part))**2 - 6.d0*C_kern_e)*GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)
					gkernse1(1:3,j+H_N_part)=gkernse1(1:3,j+H_N_part)-2.d0*C_kern_e*rij_ese1_old(1:3,j+H_N_part,i+H_N_part)* &
					  GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)
				END DO
			END DO
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					lkernse=lkernse+((2.d0*C_kern_e*rij_ese2_old(0,j,i))**2 - 6.d0*C_kern_e)*GDse2_up_old(j,i)*IGDse2_up_old(i,j)
					gkernse2(1:3,j)=gkernse2(1:3,j)-2.d0*C_kern_e*rij_ese2_old(1:3,j,i)*GDse2_up_old(j,i)*IGDse2_up_old(i,j)
				END DO
			END DO
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					lkernse=lkernse+((2.d0*C_kern_e*rij_ese2_old(0,j+H_N_part,i+H_N_part))**2 - 6.d0*C_kern_e)*GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)
					gkernse2(1:3,j+H_N_part)=gkernse2(1:3,j+H_N_part)-2.d0*C_kern_e*rij_ese2_old(1:3,j+H_N_part,i+H_N_part)* &
					  GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)
				END DO
			END DO
			lkernse=lkernse*0.5d0
		CASE ('gdc')
			minL=MINVAL(L)
			DO j = 1, H_N_part, 1
				jdw=j+H_N_part
				DO i = 1, H_N_part, 1
					idw=i+H_N_part					
					IF ( rij_ese1_old(0,j,i)<l1_kern_e ) THEN
						frf1(1:3)=-2.d0*C_kern_e*rij_ese1_old(1:3,j,i)
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse1(1:3,j)=gkernse1(1:3,j)+GDse1_up_old(j,i)*IGDse1_up_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse1_up_old(j,i)*IGDse1_up_old(i,j)*frf1(0)
						lkernse=lkernse+GDse1_up_old(j,i)*IGDse1_up_old(i,j)*(-6.d0*C_kern_e)
					ELSE IF (rij_ese1_old(0,j,i)<l2_kern_e) THEN
						frf1(1:3)=alpha1_kern_e*rij_ese1_old(1:3,j,i)/( rij_ese1_old(0,j,i)*((rij_ese1_old(0,j,i)-0.5d0*minL)**2) )
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse1(1:3,j)=gkernse1(1:3,j)+GDse1_up_old(j,i)*IGDse1_up_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse1_up_old(j,i)*IGDse1_up_old(i,j)*frf1(0)
						lkernse=lkernse+GDse1_up_old(j,i)*IGDse1_up_old(i,j)* &
						  (alpha1_kern_e/((rij_ese1_old(0,j,i)-0.5d0*minL)**2))*(2.d0/rij_ese1_old(0,j,i)- &
						  2.d0/(rij_ese1_old(0,j,i)-0.5d0*minL))
					ELSE
						frf1(1:3)=-n_kern_e*beta1_kern_e*rij_ese1_old(1:3,j,i)*(rij_ese1_old(0,j,i)**(n_kern_e-2))
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse1(1:3,j)=gkernse1(1:3,j)+GDse1_up_old(j,i)*IGDse1_up_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse1_up_old(j,i)*IGDse1_up_old(i,j)*frf1(0)
						lkernse=lkernse+GDse1_up_old(j,i)*IGDse1_up_old(i,j)* &
						  n_kern_e*(n_kern_e+1)*beta1_kern_e*(rij_ese1_old(0,j,i)**(n_kern_e-2))						
					END IF
					IF ( rij_ese1_old(0,jdw,idw)<l1_kern_e ) THEN
						frf1(1:3)=-2.d0*C_kern_e*rij_ese1_old(1:3,jdw,idw)
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse1(1:3,jdw)=gkernse1(1:3,jdw)+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)*frf1(0)
						lkernse=lkernse+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)*(-6.d0*C_kern_e)
					ELSE IF (rij_ese1_old(0,jdw,idw)<l2_kern_e) THEN
						frf1(1:3)=alpha1_kern_e*rij_ese1_old(1:3,jdw,idw)/( rij_ese1_old(0,jdw,idw)*((rij_ese1_old(0,jdw,idw)-0.5d0*minL)**2) )
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse1(1:3,jdw)=gkernse1(1:3,jdw)+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)*frf1(0)
						lkernse=lkernse+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)* &
						  (alpha1_kern_e/((rij_ese1_old(0,jdw,idw)-0.5d0*minL)**2))*(2.d0/rij_ese1_old(0,jdw,idw)- &
						  2.d0/(rij_ese1_old(0,jdw,idw)-0.5d0*minL))
					ELSE
						frf1(1:3)=-n_kern_e*beta1_kern_e*rij_ese1_old(1:3,jdw,idw)*(rij_ese1_old(0,jdw,idw)**(n_kern_e-2))
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse1(1:3,jdw)=gkernse1(1:3,jdw)+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)*frf1(0)
						lkernse=lkernse+GDse1_dw_old(j,i)*IGDse1_dw_old(i,j)* &
						  n_kern_e*(n_kern_e+1)*beta1_kern_e*(rij_ese1_old(0,jdw,idw)**(n_kern_e-2))
					END IF
					IF ( rij_ese2_old(0,j,i)<l1_kern_e ) THEN
						frf1(1:3)=-2.d0*C_kern_e*rij_ese2_old(1:3,j,i)
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse2(1:3,j)=gkernse2(1:3,j)+GDse2_up_old(j,i)*IGDse2_up_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse2_up_old(j,i)*IGDse2_up_old(i,j)*frf1(0)
						lkernse=lkernse+GDse2_up_old(j,i)*IGDse2_up_old(i,j)*(-6.d0*C_kern_e)
					ELSE IF (rij_ese2_old(0,j,i)<l2_kern_e) THEN
						frf1(1:3)=alpha1_kern_e*rij_ese2_old(1:3,j,i)/( rij_ese2_old(0,j,i)*((rij_ese2_old(0,j,i)-0.5d0*minL)**2) )
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse2(1:3,j)=gkernse2(1:3,j)+GDse2_up_old(j,i)*IGDse2_up_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse2_up_old(j,i)*IGDse2_up_old(i,j)*frf1(0)
						lkernse=lkernse+GDse2_up_old(j,i)*IGDse2_up_old(i,j)* &
						  (alpha1_kern_e/((rij_ese2_old(0,j,i)-0.5d0*minL)**2))*(2.d0/rij_ese2_old(0,j,i)- &
						  2.d0/(rij_ese2_old(0,j,i)-0.5d0*minL))
					ELSE
						frf1(1:3)=-n_kern_e*beta1_kern_e*rij_ese2_old(1:3,j,i)*(rij_ese2_old(0,j,i)**(n_kern_e-2))
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse2(1:3,j)=gkernse2(1:3,j)+GDse2_up_old(j,i)*IGDse2_up_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse2_up_old(j,i)*IGDse2_up_old(i,j)*frf1(0)
						lkernse=lkernse+GDse2_up_old(j,i)*IGDse2_up_old(i,j)* &
						  n_kern_e*(n_kern_e+1)*beta1_kern_e*(rij_ese2_old(0,j,i)**(n_kern_e-2))
					END IF
					IF ( rij_ese2_old(0,jdw,idw)<l1_kern_e ) THEN
						frf1(1:3)=-2.d0*C_kern_e*rij_ese2_old(1:3,jdw,idw)
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse2(1:3,jdw)=gkernse2(1:3,jdw)+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)*frf1(0)
						lkernse=lkernse+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)*(-6.d0*C_kern_e)
					ELSE IF (rij_ese2_old(0,jdw,idw)<l2_kern_e) THEN
						frf1(1:3)=alpha1_kern_e*rij_ese2_old(1:3,jdw,idw)/( rij_ese2_old(0,jdw,idw)*((rij_ese2_old(0,jdw,idw)-0.5d0*minL)**2) )
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse2(1:3,jdw)=gkernse2(1:3,jdw)+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)*frf1(0)
						lkernse=lkernse+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)* &
						  (alpha1_kern_e/((rij_ese2_old(0,jdw,idw)-0.5d0*minL)**2))*(2.d0/rij_ese2_old(0,jdw,idw)- &
						  2.d0/(rij_ese2_old(0,jdw,idw)-0.5d0*minL))
					ELSE
						frf1(1:3)=-n_kern_e*beta1_kern_e*rij_ese2_old(1:3,jdw,idw)*(rij_ese2_old(0,jdw,idw)**(n_kern_e-2))
						frf1(0)=DOT_PRODUCT(frf1(1:3),frf1(1:3))
						gkernse2(1:3,jdw)=gkernse2(1:3,jdw)+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)*frf1(1:3)
						lkernse=lkernse+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)*frf1(0)
						lkernse=lkernse+GDse2_dw_old(j,i)*IGDse2_dw_old(i,j)* &
						  n_kern_e*(n_kern_e+1)*beta1_kern_e*(rij_ese2_old(0,jdw,idw)**(n_kern_e-2))
					END IF
				END DO
			END DO
			lkernse=lkernse*0.5d0
		CASE ('gdp')
			frf2(1:3)=PI/L(1:3)
			frf6=-2.d0*C_kern_e
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1				
					frf5=GDse1_up_old(j,i)*IGDse1_up_old(i,j)
					frfs1(1:3)=DCOS(rij_ese1_old(1:3,j,i)*frf2(1:3))
					frf2s1(1:3)=frfs1(1:3)*rij_ese1_old(1:3,j,i)/rij_ese1_old(0,j,i)        !dr/dx
					frf3s1(1:3)=-rijpc_ese1_old(1:3,j,i)*(frf2(1:3)**2)*rijpc_ese1_old(1:3,j,i)/rijpc_ese1_old(0,j,i) + &     !d/dx dr'/dx
					  (frfs1(1:3)**2) * (1.d0 - (rijpc_ese1_old(1:3,j,i)/rijpc_ese1_old(0,j,i))**2 )/rijpc_ese1_old(0,j,i)
					frfs1(1:3)=frf6*rijpc_ese1_old(0,j,i)
					
					DO i3 = 1, 3, 1
						gkernse1(i3,j)=gkernse1(i3,j)+frf2s1(i3)*frfs1(i3)*frf5
						lkernse=lkernse+frf3s1(i3)*frfs1(i3)*frf5      !ok
						lkernse=lkernse+(frf2s1(i3)*frf2s1(i3))*frf6*frf5
					END DO
					lkernse=lkernse+(DOT_PRODUCT(gkernse1(1:3,j),gkernse1(1:3,j)))*frf5
				END DO
			END DO
			DO j = H_N_part+1, N_part, 1
				j_SD=j-H_N_part
				DO i = H_N_part+1, N_part, 1
					i_SD=i-H_N_part
					frf5=GDse1_dw_old(j_SD,i_SD)*IGDse1_dw_old(i_SD,j_SD)
					frfs1(1:3)=DCOS(rij_ese1_old(1:3,j,i)*frf2(1:3))
					frf2s1(1:3)=frfs1(1:3)*rij_ese1_old(1:3,j,i)/rij_ese1_old(0,j,i)
					frf3s1(1:3)=-rijpc_ese1_old(1:3,j,i)*(frf2(1:3)**2)*rijpc_ese1_old(1:3,j,i)/rijpc_ese1_old(0,j,i) + &     !d/dx dr'/dx
					  (frfs1(1:3)**2) * (1.d0 - (rijpc_ese1_old(1:3,j,i)/rijpc_ese1_old(0,j,i))**2 )/rijpc_ese1_old(0,j,i)
					frfs1(1:3)=frf6*rijpc_ese1_old(0,j,i)
					
					DO i3 = 1, 3, 1
						gkernse1(i3,j)=gkernse1(i3,j)+frf2s1(i3)*frfs1(i3)*frf5
						lkernse=lkernse+frf3s1(i3)*frfs1(i3)*frf5      !ok
						lkernse=lkernse+(frf2s1(i3)*frf2s1(i3))*frf6*frf5
					END DO
					lkernse=lkernse+(DOT_PRODUCT(gkernse1(1:3,j),gkernse1(1:3,j)))*frf5
				END DO
			END DO
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					frf5=GDse2_up_old(j,i)*IGDse2_up_old(i,j)
					frfs2(1:3)=DCOS(rij_ese2_old(1:3,j,i)*frf2(1:3))
					frf2s2(1:3)=frfs2(1:3)*rij_ese2_old(1:3,j,i)/rij_ese2_old(0,j,i)
					frf3s2(1:3)=-rijpc_ese2_old(1:3,j,i)*(frf2(1:3)**2)*rijpc_ese2_old(1:3,j,i)/rijpc_ese2_old(0,j,i) + &     !d/dx dr'/dx
					  (frfs2(1:3)**2) * (1.d0 - (rijpc_ese2_old(1:3,j,i)/rijpc_ese2_old(0,j,i))**2 )/rijpc_ese2_old(0,j,i)
					frfs2(1:3)=frf6*rijpc_ese2_old(0,j,i)
					
					DO i3 = 1, 3, 1
						gkernse2(i3,j)=gkernse2(i3,j)+frf2s2(i3)*frfs2(i3)*frf5
						lkernse=lkernse+frf3s2(i3)*frfs2(i3)*frf5      !ok
						lkernse=lkernse+(frf2s2(i3)*frf2s2(i3))*frf6*frf5
					END DO
					lkernse=lkernse+(DOT_PRODUCT(gkernse2(1:3,j),gkernse2(1:3,j)))*frf5
				END DO
			END DO
			DO j = H_N_part+1, N_part, 1
				j_SD=j-H_N_part
				DO i = H_N_part+1, N_part, 1
					i_SD=i-H_N_part
					frf5=GDse2_dw_old(j_SD,i_SD)*IGDse2_dw_old(i_SD,j_SD)
					frfs2(1:3)=DCOS(rij_ese2_old(1:3,j,i)*frf2(1:3))
					frf2s2(1:3)=frfs2(1:3)*rij_ese2_old(1:3,j,i)/rij_ese2_old(0,j,i)
					frf3s2(1:3)=-rijpc_ese2_old(1:3,j,i)*(frf2(1:3)**2)*rijpc_ese2_old(1:3,j,i)/rijpc_ese2_old(0,j,i) + &     !d/dx dr'/dx
					  (frfs2(1:3)**2) * (1.d0 - (rijpc_ese2_old(1:3,j,i)/rijpc_ese2_old(0,j,i))**2 )/rijpc_ese2_old(0,j,i)
					frfs2(1:3)=frf6*rijpc_ese2_old(0,j,i)
					
					DO i3 = 1, 3, 1
						gkernse2(i3,j)=gkernse2(i3,j)+frf2s2(i3)*frfs2(i3)*frf5
						lkernse=lkernse+frf3s2(i3)*frfs2(i3)*frf5      !ok
						lkernse=lkernse+(frf2s2(i3)*frf2s2(i3))*frf6*frf5
					END DO
					lkernse=lkernse+(DOT_PRODUCT(gkernse2(1:3,j),gkernse2(1:3,j)))*frf5
				END DO
			END DO
			lkernse=lkernse*0.5d0
		CASE ('atm')
			DO i = 1, N_part, 1
				gkernse1(1:3,i)=-C_kern_e*rij_ese1_old(1:3,i,i)/rij_ese1_old(0,i,i)
				gkernse2(1:3,i)=-C_kern_e*rij_ese2_old(1:3,i,i)/rij_ese2_old(0,i,i)
				lkernse=lkernse-C_kern_e*(1.d0/rij_ese1_old(0,i,i)+1.d0/rij_ese2_old(0,i,i))
			END DO
			lkernse=lkernse+REAL(N_part,8)*C_kern_e*C_kern_e
		CASE ('atc')
			L38=3.d0*MINVAL(L)/8.d0
			DO i = 1, N_part, 1
				IF ( rij_ese1_old(0,i,i)<=L38 ) THEN
					gkernse1(1:3,i)=-C_kern_e*rij_ese1_old(1:3,i,i)/rij_ese1_old(0,i,i)
					lkernse=lkernse-C_kern_e*1.d0/rij_ese1_old(0,i,i)
				ELSE
					gkernse1(1:3,i)=-alpha1_kern_e*10.d0*((rij_ese1_old(0,i,i)-L38)**9)*rij_ese1_old(1:3,i,i)/rij_ese1_old(0,i,i)
					lkernse=lkernse+alpha1_kern_e*10.d0*((rij_ese1_old(0,i,i)-L38)**8)* &
					   ( (rij_ese1_old(0,i,i)-L38)/rij_ese1_old(0,i,i) + 4.5d0 )
				END IF
				IF ( rij_ese2_old(0,i,i)<=L38 ) THEN
					gkernse2(1:3,i)=-C_kern_e*rij_ese2_old(1:3,i,i)/rij_ese2_old(0,i,i)
					lkernse=lkernse-C_kern_e*1.d0/rij_ese2_old(0,i,i)
				ELSE
					gkernse2(1:3,i)=-alpha1_kern_e*10.d0*((rij_ese2_old(0,i,i)-L38)**9)*rij_ese2_old(1:3,i,i)/rij_ese2_old(0,i,i)
					lkernse=lkernse+alpha1_kern_e*10.d0*((rij_ese2_old(0,i,i)-L38)**8)* &
					   ( (rij_ese2_old(0,i,i)-L38)/rij_ese2_old(0,i,i) + 4.5d0  )
				END IF
				lkernse=lkernse+0.5d0*(DOT_PRODUCT(gkernse1(1:3,i),gkernse1(1:3,i))+DOT_PRODUCT(gkernse2(1:3,i),gkernse2(1:3,i)))
			END DO
		END SELECT
		
		!calcolo i prodotti misti
		mix_prod=0.d0
		DO j = 1, H_N_part, 1
			mix_prod=mix_prod+DOT_PRODUCT(gsdee_up(1:3,j),gee(1:3,j)+gep(1:3,j)+0.5d0*(gkernse1(1:3,j)+gkernse2(1:3,j)))
			mix_prod=mix_prod+DOT_PRODUCT(gee(1:3,j),gep(1:3,j)+0.5d0*(gkernse1(1:3,j)+gkernse2(1:3,j)))
			mix_prod=mix_prod+DOT_PRODUCT(gep(1:3,j),0.5d0*(gkernse1(1:3,j)+gkernse2(1:3,j)))
		END DO
		DO j = H_N_part+1, N_part, 1
			mix_prod=mix_prod+DOT_PRODUCT(gsdee_dw(1:3,j-H_N_part),gee(1:3,j)+gep(1:3,j)+0.5d0*(gkernse1(1:3,j)+gkernse2(1:3,j)))
			mix_prod=mix_prod+DOT_PRODUCT(gee(1:3,j),gep(1:3,j)+0.5d0*(gkernse1(1:3,j)+gkernse2(1:3,j)))
			mix_prod=mix_prod+DOT_PRODUCT(gep(1:3,j),0.5d0*(gkernse1(1:3,j)+gkernse2(1:3,j)))
		END DO
		
		!Scrivo il contributo all'energia cinetica dato dagli elettroni
		!Energia cinetica di Pandariphande-Bethe		
		E_kin=E_kin-(lsdee+lee+lep+lkernse+2.d0*mix_prod)*hbar*hbar/(2.d0*MASS_e*N_part)
		
		!Energia cinetica di Jackson-Feenberg
		frf1(0)=0.d0
		DO j = 1, H_N_part, 1
			frfc3(1:3)=gsdee_up(1:3,j)+gee(1:3,j)+gep(1:3,j)
			frf1(0)=frf1(0)+DOT_PRODUCT(frfc3(1:3)+gkernse1(1:3,j),frfc3(1:3)+gkernse2(1:3,j))
		END DO
		DO j = H_N_part+1, N_part, 1
			frfc3(1:3)=gsdee_dw(1:3,j-H_N_part)+gee(1:3,j)+gep(1:3,j)
			frf1(0)=frf1(0)+DOT_PRODUCT(frfc3(1:3)+gkernse1(1:3,j),frfc3(1:3)+gkernse2(1:3,j))
		END DO
		
		E_JF=E_JF+frf1(0)*hbar*hbar/(2.d0*MASS_e*N_part)
		END IF elettroni
				
	END SUBROUTINE energia_cinetica
!-----------------------------------------------------------------------
	SUBROUTINE energia_potenziale(E_pot)
		USE generic_tools
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
		INTEGER :: i, j, k
		REAL (KIND=8) :: frf1, frf3, delta_R(1:3,1:num_k_ewald)
		REAL (KIND=8) :: rho_k(1:num_k_ewald)
		REAL (KIND=8), INTENT(OUT) :: E_pot
		REAL (KIND=8) :: r_prova(0:3)
		
		E_pot=0.d0
		
		IF (flag_somme_ewald) THEN
			!Calcolo il potenziale di Coulomb fra protoni ed elettroni con le somme di Ewald
			rho_k=0.d0
			DO k = 2, num_k_ewald, 1
				frf1=0.d0
				frf3=0.d0
				DO i = 1, N_part, 1
					frf1=frf1-DSIN(DOT_PRODUCT(k_fermi(1:3,k),re_old(1:3,i)))+ &
					  DSIN(DOT_PRODUCT(k_fermi(1:3,k),rp_old(1:3,i)))
					frf3=frf3-DCOS(DOT_PRODUCT(k_fermi(1:3,k),re_old(1:3,i)))+ &
					  DCOS(DOT_PRODUCT(k_fermi(1:3,k),rp_old(1:3,i)))
				END DO
				rho_k(k)=frf1*frf1+frf3*frf3
			END DO
			DO k = 2, num_k_ewald, 1
				E_pot=E_pot+rho_k(k)*DEXP(-k_fermi(0,k)*0.25d0/(alpha_ewald*alpha_ewald))/k_fermi(0,k)
			END DO
			E_pot=E_pot*2.d0*PI/(L(1)*L(2)*L(3))
			DO j = 1, N_part-1, 1
				DO i = j+1, N_part, 1
					E_pot=E_pot+DERFC(alpha_ewald*rij_ee_old(0,i,j))/rij_ee_old(0,i,j)
				END DO
			END DO
			DO j = 1, N_part-1, 1
				DO i = j+1, N_part, 1
					E_pot=E_pot+DERFC(alpha_ewald*rij_pp_old(0,i,j))/rij_pp_old(0,i,j)
				END DO
			END DO
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
					E_pot=E_pot-DERFC(alpha_ewald*rij_ep_old(0,i,j))/rij_ep_old(0,i,j)
				END DO
			END DO
			E_pot=E_pot-2.d0*N_part*alpha_ewald/DSQRT(PI)
			
		ELSE
			DO i = 1, N_part-1, 1
				DO j = i+1, N_part, 1
					E_pot=E_pot+1.d0/rij_ee_old(0,i,j)
					E_pot=E_pot+1.d0/rij_pp_old(0,i,j)
				END DO
			END DO
			DO i = 1, N_part, 1
				DO j = 1, N_part, 1
					E_pot=E_pot-1.d0/rij_ep_old(0,i,j)
				END DO
			END DO
			
			DO k = 2, num_k_ewald, 1
				DO i = 1, 3, 1
					delta_R(i,k)=k_fermi(i,k)*L(i)*MINVAL(L)/(2.d0*PI)
				END DO
				DO i = 1, N_part, 1
					DO j = 1, N_part, 1
						E_pot=E_pot+1.d0/DSQRT(DOT_PRODUCT(re_old(1:3,i)-re_old(1:3,j)+delta_R(1:3,k), &
						  re_old(1:3,i)-re_old(1:3,j)+delta_R(1:3,k)))
						E_pot=E_pot+1.d0/DSQRT(DOT_PRODUCT(rp_old(1:3,i)-rp_old(1:3,j)+delta_R(1:3,k), &
						  rp_old(1:3,i)-rp_old(1:3,j)+delta_R(1:3,k)))
						E_pot=E_pot-1.d0/DSQRT(DOT_PRODUCT(re_old(1:3,i)-rp_old(1:3,j)+delta_R(1:3,k), &
						  re_old(1:3,i)-rp_old(1:3,j)+delta_R(1:3,k)))
						E_pot=E_pot-1.d0/DSQRT(DOT_PRODUCT(rp_old(1:3,i)-re_old(1:3,j)+delta_R(1:3,k), &
						  rp_old(1:3,i)-re_old(1:3,j)+delta_R(1:3,k)))
					END DO
				END DO
			END DO
			
			!!!!CON PC
			!!!DO i = 1, N_part-1, 1
			!!!	DO j = i+1, N_part, 1
			!!!		E_pot=E_pot+1.d0/rijpc_ee_old(0,i,j)
			!!!		r_prova(1:3)=L(1:3)*DSIN(PI*rij_pp_old(1:3,i,j)/L(1:3))/PI
			!!!		r_prova(0)=DSQRT(DOT_PRODUCT(r_prova(1:3),r_prova(1:3)))
			!!!		E_pot=E_pot+1.d0/r_prova(0)
			!!!	END DO
			!!!END DO
			!!!DO i = 1, N_part, 1
			!!!	DO j = 1, N_part, 1
			!!!		E_pot=E_pot-1.d0/rijpc_ep_old(0,i,j)
			!!!	END DO
			!!!END DO
			
		END IF
		
		E_pot=E_pot*K_coulomb/N_part
		
	END SUBROUTINE energia_potenziale
!-----------------------------------------------------------------------

	SUBROUTINE peso_funzione_onda(w)
		IMPLICIT NONE
		REAL (KIND=8), INTENT(OUT) :: w
		w=1.d0
		IF (SDse_kind=='pw2') THEN
		  !  w=w*DSIGN(1.d0,REAL(((detSDse1_up_new**2))*(DCONJG((detSDse2_up_new**2)))* &
		  !((detSDse1_dw_new**2))*(DCONJG((detSDse2_dw_new**2))),8))
		  !  w=w*DSIGN(1.d0,REAL(((detSDse1_up_new**2)),8)*REAL((DCONJG((detSDse2_up_new**2))),8)* &
		  !REAL(((detSDse1_dw_new**2)),8)*REAL((DCONJG((detSDse2_dw_new**2))),8))
		ELSE IF (SDse_kind/='no_') THEN
			w=w*DSIGN(1.d0,REAL((detSDse1_up_new)*(DCONJG(detSDse2_up_new))* &
			  (detSDse1_dw_new)*(DCONJG(detSDse2_dw_new)),8))
		END IF
		IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc')) THEN
		 	w=w*DSIGN(1.d0,detGDse1_up_new*detGDse2_up_new*detGDse1_dw_new*detGDse2_dw_new)
		END IF
		IF (Jsesp_kind=='gsd') THEN
			w=w*DSIGN(1.d0,detGDsp1_up_new*detGDsp2_up_new*detGDsp1_dw_new*detGDsp2_dw_new)
		END IF
		
	END SUBROUTINE peso_funzione_onda
!-----------------------------------------------------------------------

	SUBROUTINE derivata_SDe_atm(O)
		USE walkers
		IMPLICIT NONE
		INTEGER :: i, j
		REAL (KIND=8) :: O      !O(1) - Gsesp
		O=0.d0
		
		DO j = 1, H_N_part, 1
			DO i = 1, H_N_part, 1
				O=O-rij_ep_old(0,i,j)*SDe_up_old(i,j)*ISDe_up_old(j,i)
			END DO
		END DO
		
		DO j = 1, H_N_part, 1
			DO i = 1, H_N_part, 1
				O=O-rij_ep_old(0,i+H_N_part,j+H_N_part)*SDe_dw_old(i,j)*ISDe_dw_old(j,i)
			END DO
		END DO
		
	END SUBROUTINE derivata_SDe_atm
	
!-----------------------------------------------------------------------

	SUBROUTINE derivata_SDe_bat(O)
		USE walkers
		IMPLICIT NONE
		INTEGER :: i, j
		REAL (KIND=8) :: O      !O(1) - Gsesp
		O=0.d0
	
		DO j = 1, H_N_part, 1
			DO i = 1, H_N_part, 1
				O=O-(rij_ep_old(0,i,j)*DEXP(-C_atm*rij_ep_old(0,i,j)) +  &
				  rij_ep_old(0,i,j+H_N_part)*DEXP(-C_atm*rij_ep_old(0,i,j+H_N_part)))*ISDe_up_old(j,i)
			END DO
		END DO
	
		DO j = H_N_part+1, N_part, 1
			DO i = H_N_part+1, N_part, 1
				O=O-(rij_ep_old(0,i,j)*DEXP(-C_atm*rij_ep_old(0,i,j)) + &
				  rij_ep_old(0,i,j-H_N_part)*DEXP(-C_atm*rij_ep_old(0,i,j-H_N_part)) )*ISDe_dw_old(j-H_N_part,i-H_N_part)
			END DO
		END DO
	
	END SUBROUTINE derivata_SDe_bat
!-----------------------------------------------------------------------

	SUBROUTINE derivata_SDe_1sb(O)
		USE walkers
		IMPLICIT NONE
		INTEGER :: i, j, ip, iadd, beta
		REAL (KIND=8) :: O(1:3)
		REAL(KIND=8) :: phi1(1:H_N_part,1:H_N_part) !per derivata C_atm
		REAL(KIND=8) :: phi2(1:H_N_part,1:H_N_part) !per derivata A_POT_se
		REAL(KIND=8) :: phi3(1:H_N_part,1:H_N_part) !per derivata D_POT_se
		REAL(KIND=8) :: norm, q(0:3), sigm, deriv
		
		norm=1.d0
		O=0.d0
		
		DO iadd = 0, H_N_part, H_N_part
		   DO j = 1, H_N_part, 1
		      DO i = 1, H_N_part, 1
	             q(1:3)=re_old(1:3,i+iadd)
	             DO ip = 1, N_part, 1
	                sigm=1.d0/(1.d0+DEXP(A_POT_se*(rij_ep_old(0,j+iadd,ip)-D_POT_se)))
	                q(1:3)=q(1:3)-rp_old(1:3,ip)*sigm
	             END DO
	             q(1:3)=q(1:3)-L(1:3)*DNINT(q(1:3)/L(1:3))
	             q(0)=DSQRT(DOT_PRODUCT(q(1:3),q(1:3)))
	             !derivata di C_atm
			     deriv=-q(0)
		         phi1(i,j)=deriv*norm*DEXP(-C_atm*q(0))
		         !derivata di A_POT_se
		         deriv=0.d0
	             DO ip = 1, N_part, 1
					sigm=DEXP(A_POT_se*(rij_ep_old(0,j+iadd,ip)-D_POT_se))
	                DO beta = 1, 3, 1
	                   deriv=deriv+q(beta)*rp_old(beta,ip)* &
	                      (- sigm * (rij_ep_old(0,j+iadd,ip)-D_POT_se) / ((1.d0 + sigm)**2) )
	                END DO
	             END DO
	             deriv=deriv*C_atm/q(0)
		         phi2(i,j)=deriv*norm*DEXP(-C_atm*q(0))
		         !derivata di A_POT_se
		         deriv=0.d0
	             DO ip = 1, N_part, 1
					sigm=DEXP(A_POT_se*(rij_ep_old(0,j+iadd,ip)-D_POT_se))
	                DO beta = 1, 3, 1
	                   deriv=deriv+q(beta)*rp_old(beta,ip)* &
	                      ( sigm * A_POT_se / ((1.d0 + sigm)**2) )
	                END DO
	             END DO
	             deriv=deriv*C_atm/q(0)
		         phi3(i,j)=deriv*norm*DEXP(-C_atm*q(0))
		      END DO
		   END DO
		   
		   IF (iadd==0) THEN
		      DO j = 1, H_N_part, 1
		         DO i = 1, H_N_part, 1
		            O(1)=O(1)+phi1(i,j)*ISDe_up_old(j,i)
		            O(2)=O(2)+phi2(i,j)*ISDe_up_old(j,i)
		            O(3)=O(3)+phi3(i,j)*ISDe_up_old(j,i)
		         END DO
		      END DO
		   ELSE IF (iadd==H_N_part) THEN
		      DO j = 1, H_N_part, 1
		         DO i = 1, H_N_part, 1
		            O(1)=O(1)+phi1(i,j)*ISDe_dw_old(j,i)
		            O(2)=O(2)+phi2(i,j)*ISDe_dw_old(j,i)
		            O(3)=O(3)+phi3(i,j)*ISDe_dw_old(j,i)
		         END DO
		      END DO
		   END IF
		END DO
				
	END SUBROUTINE derivata_SDe_1sb
	
!-----------------------------------------------------------------------

	SUBROUTINE derivata_SDe_atp(O)
		USE walkers
		IMPLICIT NONE
		INTEGER :: i, j, i_sd, j_sd
		REAL (KIND=8) :: O      !O(1) - Gsesp
		INTEGER :: pvt(1:H_N_part), info, perm
		REAL (KIND=8) :: SD(1:H_N_part,1:H_N_part), ISD(1:H_N_part,1:H_N_part), detSD
		O=0.d0
		
		!up
		DO j = 1, H_N_part, 1
			DO i = 1, H_N_part, 1
				SD(i,j)=-rijpc_ep_old(0,i,j)*DEXP(-C_atm*rijpc_ep_old(0,i,j))
				ISD(i,j)=SD(i,j)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU SDe_atp_up'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDe_up_old,8)
		
		!dw
		DO j = H_N_part+1, N_part, 1
			j_sd=j-H_N_part
			DO i = H_N_part+1, N_part, 1
				i_sd=i-H_N_part
				SD(i_sd,j_sd)=-rijpc_ep_old(0,i,j)*DEXP(-C_atm*rijpc_ep_old(0,i,j))
				ISD(i_sd,j_sd)=SD(i_sd,j_sd)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU SDe_atp_dw'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDe_dw_old,8)
		
	END SUBROUTINE derivata_SDe_atp
!-----------------------------------------------------------------------

	SUBROUTINE derivata_SDe_HAR(O)
		IMPLICIT NONE
		INTEGER :: j, i1, i2, i3, alpha, beta, cont
		REAL (KIND=8), DIMENSION(:) :: O(:)
		REAL (KIND=8) :: diff
		COMPLEX (KIND=8) :: zzz
		
		cont=0
		DO beta = 1, N_pw_lda, 1
			DO alpha = 1, beta-1, 1
			cont=cont+2
				DO i1 = 1, H_N_part, 1
					DO i2 = 1, H_N_part, 1
						DO i3 = 1, H_N_part, 1
						diff=(autoenergie_dnfH(i3)-autoenergie_dnfH(i1))
						IF (DABS(diff)>0.001) THEN
							zzz=ISDe_up_old(i1,i2)*SDe_up_old(i2,i3)* &
							  DCONJG(fattori_pw_lda(alpha,i3))*fattori_pw_lda(beta,i1)/diff
							O(cont-1)=O(cont-1)-REAL(zzz,8)
							O(cont)=O(cont)-DIMAG(zzz)
						END IF
						END DO
					END DO
				END DO
			END DO
		END DO
		DO beta = 1, N_pw_lda, 1
			alpha=beta
			cont=cont+1
				DO i1 = 1, H_N_part, 1
					DO i2 = 1, H_N_part, 1
						DO i3 = 1, H_N_part, 1
						diff=(autoenergie_dnfH(i3)-autoenergie_dnfH(i1))
						IF (DABS(diff)>0.001) THEN
							zzz=ISDe_up_old(i1,i2)*SDe_up_old(i2,i3)* &
							  DCONJG(fattori_pw_lda(alpha,i3))*fattori_pw_lda(beta,i1)/diff
							O(cont-1)=O(cont-1)-REAL(zzz,8)
							O(cont)=O(cont)-DIMAG(zzz)
						END IF
						END DO
					END DO
				END DO
		END DO
				
	END SUBROUTINE derivata_SDe_HAR
!-----------------------------------------------------------------------

	SUBROUTINE derivata_Jee_YUK(O)
		IMPLICIT NONE
		INTEGER :: i, j
		REAL (KIND=8) :: frf1
		REAL (KIND=8), DIMENSION(:) :: O(:)
		!O(1) - Aee  ;  O(2) - Fee_uu  ;  O(3) - Fee_ud
		!O(1) - Aee  ;  O(2) - Aee_ud  ;  O(3) - Fee_uu  ;  O(4) - Fee_ud
		
		!split_Aee			
		!	split_Fee		O(1) - Aep  ;  O(2) - Aep_ud  ;  O(3) - Fep_uu  ;  O(4) - Fep_ud
		!	not split_Fee	O(1) - Aep  ;  O(2) - Aep_ud  ;  O(3) - Fep_uu
		!not split_Aee		
		!	split_Fee		O(1) - Aep  ;  O(2) - Fep_uu  ;  O(3) - Fep_ud
		!	not split_Fee	O(1) - Aep  ;  O(2) - Fep_uu
		
		O=0.d0
		IF (opt_A_Jee) THEN
			IF (opt_F_Jee) THEN
				SELECT CASE(Jee_kind)
				CASE('yuk')
					IF (split_Aee) THEN
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
										O(3)=O(3)-Aee_yuk*frf1
										O(1)=O(1)-(1.d0-frf1)/rij_ee_old(0,i,j)
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
										O(4)=O(4)-Aee_ud_yuk*frf1
										O(2)=O(2)-(1.d0-frf1)/rij_ee_old(0,i,j)
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
										O(3)=O(3)-Aee_yuk*frf1
										O(1)=O(1)-(1.d0-frf1)/rij_ee_old(0,i,j)
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
										O(3)=O(3)-Aee_ud_yuk*frf1
										O(2)=O(2)-(1.d0-frf1)/rij_ee_old(0,i,j)
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
										O(2)=O(2)-Aee_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
										O(3)=O(3)-Aee_ud_yuk*frf1
									END IF
									O(1)=O(1)-(1.d0-frf1)/rij_ee_old(0,i,j)
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
									O(2)=O(2)-Aee_yuk*frf1
									O(1)=O(1)-(1.d0-frf1)/rij_ee_old(0,i,j)
								END DO
							END DO
						END IF
					END IF
				CASE('yup')
					IF (split_Aee) THEN
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
										O(3)=O(3)-Aee_yuk*frf1
										O(1)=O(1)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
										O(4)=O(4)-Aee_ud_yuk*frf1
										O(2)=O(2)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
										O(3)=O(3)-Aee_yuk*frf1
										O(1)=O(1)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
										O(3)=O(3)-Aee_ud_yuk*frf1
										O(2)=O(2)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
										O(2)=O(2)-Aee_yuk*frf1
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
										O(3)=O(3)-Aee_ud_yuk*frf1
									END IF
									O(1)=O(1)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
									O(2)=O(2)-Aee_yuk*frf1
									O(1)=O(1)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
								END DO
							END DO
						END IF
					END IF
				END SELECT
			ELSE
				SELECT CASE(Jee_kind)
				CASE('yuk')
					IF (split_Aee) THEN
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-(1.d0-frf1)/rij_ee_old(0,i,j)
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
										O(2)=O(2)-(1.d0-frf1)/rij_ee_old(0,i,j)
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-(1.d0-frf1)/rij_ee_old(0,i,j)
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
										O(2)=O(2)-(1.d0-frf1)/rij_ee_old(0,i,j)
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
									END IF
									O(1)=O(1)-(1.d0-frf1)/rij_ee_old(0,i,j)
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
									O(1)=O(1)-(1.d0-frf1)/rij_ee_old(0,i,j)
								END DO
							END DO
						END IF
					END IF
				CASE('yup')
					IF (split_Aee) THEN
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
										O(2)=O(2)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
										O(2)=O(2)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
									END IF
									O(1)=O(1)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
									O(1)=O(1)-(1.d0-frf1)/rijpc_ee_old(0,i,j)
								END DO
							END DO
						END IF
					END IF
				END SELECT
			END IF
		ELSE
			IF (opt_F_Jee) THEN
				SELECT CASE(Jee_kind)
				CASE('yuk')
					IF (split_Aee) THEN
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-Aee_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
										O(2)=O(2)-Aee_ud_yuk*frf1
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-Aee_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
										O(1)=O(1)-Aee_ud_yuk*frf1
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-Aee_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ee_old(0,i,j)*Fee_ud_yuk)
										O(2)=O(2)-Aee_ud_yuk*frf1
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									frf1=DEXP(-rij_ee_old(0,i,j)*Fee_yuk)
									O(1)=O(1)-Aee_yuk*frf1
								END DO
							END DO
						END IF
					END IF
				CASE('yup')
					IF (split_Aee) THEN
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-Aee_yuk*frf1
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
										O(2)=O(2)-Aee_ud_yuk*frf1
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-Aee_yuk*frf1
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
										O(1)=O(1)-Aee_ud_yuk*frf1
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fee) THEN
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
										O(1)=O(1)-Aee_yuk*frf1
									ELSE
										frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_ud_yuk)
										O(2)=O(2)-Aee_ud_yuk*frf1
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part-1, 1
								DO i = j+1, N_part, 1
									frf1=DEXP(-rijpc_ee_old(0,i,j)*Fee_yuk)
									O(1)=O(1)-Aee_yuk*frf1
								END DO
							END DO
						END IF
					END IF
				END SELECT
			END IF
		END IF
		
		O=O*0.5d0
		
	END SUBROUTINE derivata_Jee_YUK
!-----------------------------------------------------------------------
   SUBROUTINE  derivata_Jee_SPL(O)
      IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: O(:)
       INTEGER :: i, j
       REAL(KIND=8) :: td(0:Jsplee%m,0:Jsplee%Nknots), td_ud(0:Jsplee_ud%m,0:Jsplee_ud%Nknots)

       SELECT CASE(Jee_kind)
       CASE('spl')
          IF (split_Aee.OR.split_Fee) THEN
             td=0.d0
             td_ud=0.d0
             DO j = 1, N_part-1, 1
             DO i = j+1, N_part, 1
               IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
                  CALL MSPL_t_deriv(SPL=Jsplee,R=rij_ee_old(0,i,j),&
                     T_DERIV=td(0:Jsplee%m,0:Jsplee%Nknots),RESET=.FALSE.)
               ELSE
                  CALL MSPL_t_deriv(SPL=Jsplee_ud,R=rij_ee_old(0,i,j),&
                     T_DERIV=td_ud(0:Jsplee_ud%m,0:Jsplee_ud%Nknots),RESET=.FALSE.)
               END IF
             END DO
             END DO
             O(1:(Jsplee%m+1)*(Jsplee%Nknots+1))=&
                -0.5d0*RESHAPE(td(0:Jsplee%m,0:Jsplee%Nknots),(/(Jsplee%m+1)*(Jsplee%Nknots+1)/))
             O((Jsplee%m+1)*(Jsplee%Nknots+1)+1:(Jsplee%m+1)*(Jsplee%Nknots+1)+(Jsplee_ud%m+1)*(Jsplee_ud%Nknots+1))=&
                -0.5d0*RESHAPE(td(0:Jsplee_ud%m,0:Jsplee_ud%Nknots),(/(Jsplee_ud%m+1)*(Jsplee_ud%Nknots+1)/))
          ELSE
             td=0.d0
             DO j = 1, N_part-1, 1
             DO i = j+1, N_part, 1
                CALL MSPL_t_deriv(SPL=Jsplee,R=rij_ee_old(0,i,j),&
                   T_DERIV=td(0:Jsplee%m,0:Jsplee%Nknots),RESET=.FALSE.)
             END DO
             END DO
             O(1:(Jsplee%m+1)*(Jsplee%Nknots+1))=&
                -0.5d0*RESHAPE(td(0:Jsplee%m,0:Jsplee%Nknots),(/(Jsplee%m+1)*(Jsplee%Nknots+1)/))
          END IF
       CASE('spp')
          IF (split_Aee.OR.split_Fee) THEN
             td=0.d0
             td_ud=0.d0
             DO j = 1, N_part-1, 1
             DO i = j+1, N_part, 1
               IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
                  CALL MSPL_t_deriv(SPL=Jsplee,R=rijpc_ee_old(0,i,j),&
                     T_DERIV=td(0:Jsplee%m,0:Jsplee%Nknots),RESET=.FALSE.)
               ELSE
                  CALL MSPL_t_deriv(SPL=Jsplee_ud,R=rijpc_ee_old(0,i,j),&
                     T_DERIV=td_ud(0:Jsplee_ud%m,0:Jsplee_ud%Nknots),RESET=.FALSE.)
               END IF
             END DO
             END DO
             O(1:(Jsplee%m+1)*(Jsplee%Nknots+1))=&
                -0.5d0*RESHAPE(td(0:Jsplee%m,0:Jsplee%Nknots),(/(Jsplee%m+1)*(Jsplee%Nknots+1)/))
             O((Jsplee%m+1)*(Jsplee%Nknots+1)+1:(Jsplee%m+1)*(Jsplee%Nknots+1)+(Jsplee_ud%m+1)*(Jsplee_ud%Nknots+1))=&
                -0.5d0*RESHAPE(td(0:Jsplee_ud%m,0:Jsplee_ud%Nknots),(/(Jsplee_ud%m+1)*(Jsplee_ud%Nknots+1)/))
          ELSE
             td=0.d0
             DO j = 1, N_part-1, 1
             DO i = j+1, N_part, 1
                CALL MSPL_t_deriv(SPL=Jsplee,R=rijpc_ee_old(0,i,j),&
                   T_DERIV=td(0:Jsplee%m,0:Jsplee%Nknots),RESET=.FALSE.)
             END DO
             END DO
             O(1:(Jsplee%m+1)*(Jsplee%Nknots+1))=&
                -0.5d0*RESHAPE(td(0:Jsplee%m,0:Jsplee%Nknots),(/(Jsplee%m+1)*(Jsplee%Nknots+1)/))
          END IF
       END SELECT
      
   END SUBROUTINE  derivata_Jee_SPL
!-----------------------------------------------------------------------
	SUBROUTINE derivata_Jep_YUK(O)
		IMPLICIT NONE
		INTEGER :: i, j
		REAL (KIND=8) :: frf1
		REAL (KIND=8), DIMENSION(:) :: O(:)
		!split_Aep			
		!	split_Fep		O(1) - Aep  ;  O(2) - Aep_ud  ;  O(3) - Fep_uu  ;  O(4) - Fep_ud
		!	not split_Fep	O(1) - Aep  ;  O(2) - Aep_ud  ;  O(3) - Fep_uu
		!not split_Aep		
		!	split_Fep		O(1) - Aep  ;  O(2) - Fep_uu  ;  O(3) - Fep_ud
		!	not split_Fep	O(1) - Aep  ;  O(2) - Fep_uu  
		O=0.d0
		
		IF (opt_A_Jep) THEN
			IF (opt_F_Jep) THEN
				SELECT CASE(Jep_kind)
				CASE('yuk')
					IF (split_Aep) THEN
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(3)=O(3)-Aep_yuk*frf1
										O(1)=O(1)-(1.d0-frf1)/rij_ep_old(0,i,j)
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_ud_yuk)
										O(4)=O(4)-Aep_ud_yuk*frf1
										O(2)=O(2)-(1.d0-frf1)/rij_ep_old(0,i,j)
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(3)=O(3)-Aep_yuk*frf1
										O(1)=O(1)-(1.d0-frf1)/rij_ep_old(0,i,j)
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(3)=O(3)-Aep_ud_yuk*frf1
										O(2)=O(2)-(1.d0-frf1)/rij_ep_old(0,i,j)
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(2)=O(2)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_ud_yuk)
										O(3)=O(3)-Aep_yuk*frf1
									END IF
									O(1)=O(1)-(1.d0-frf1)/rij_ep_old(0,i,j)
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
									O(2)=O(2)-Aep_yuk*frf1
									O(1)=O(1)-(1.d0-frf1)/rij_ep_old(0,i,j)
								END DO
							END DO
						END IF
					END IF
				CASE('yup')
					IF (split_Aep) THEN
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(3)=O(3)-Aep_yuk*frf1
										O(1)=O(1)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_ud_yuk)
										O(4)=O(4)-Aep_ud_yuk*frf1
										O(2)=O(2)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(3)=O(3)-Aep_yuk*frf1
										O(1)=O(1)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(3)=O(3)-Aep_ud_yuk*frf1
										O(2)=O(2)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(2)=O(2)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_ud_yuk)
										O(3)=O(3)-Aep_yuk*frf1
									END IF
									O(1)=O(1)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
									O(2)=O(2)-Aep_yuk*frf1
									O(1)=O(1)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
								END DO
							END DO
						END IF
					END IF
				END SELECT
			ELSE
				SELECT CASE(Jep_kind)
				CASE('yuk')
					IF (split_Aep) THEN
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-(1.d0-frf1)/rij_ep_old(0,i,j)
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_ud_yuk)
										O(2)=O(2)-(1.d0-frf1)/rij_ep_old(0,i,j)
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-(1.d0-frf1)/rij_ep_old(0,i,j)
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(2)=O(2)-(1.d0-frf1)/rij_ep_old(0,i,j)
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_ud_yuk)
									END IF
									O(1)=O(1)-(1.d0-frf1)/rij_ep_old(0,i,j)
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
									END IF
									O(1)=O(1)-(1.d0-frf1)/rij_ep_old(0,i,j)
								END DO
							END DO
						END IF
					END IF
				CASE('yup')
					IF (split_Aep) THEN
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_ud_yuk)
										O(2)=O(2)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(2)=O(2)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_ud_yuk)
									END IF
									O(1)=O(1)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
									O(1)=O(1)-(1.d0-frf1)/rijpc_ep_old(0,i,j)
								END DO
							END DO
						END IF
					END IF
				END SELECT
			END IF
		ELSE
			IF (opt_F_Jep) THEN
				SELECT CASE(Jep_kind)
				CASE('yuk')
					IF (split_Aep) THEN
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_ud_yuk)
										O(2)=O(2)-Aep_ud_yuk*frf1
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_ud_yuk*frf1
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_ud_yuk)
										O(2)=O(2)-Aep_yuk*frf1
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rij_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_yuk*frf1
									END IF
								END DO
							END DO
						END IF
					END IF
				CASE('yup')
					IF (split_Aep) THEN
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_ud_yuk)
										O(2)=O(2)-Aep_ud_yuk*frf1
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_ud_yuk*frf1
									END IF
								END DO
							END DO
						END IF
					ELSE
						IF (split_Fep) THEN
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
										O(1)=O(1)-Aep_yuk*frf1
									ELSE
										frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_ud_yuk)
										O(2)=O(2)-Aep_yuk*frf1
									END IF
								END DO
							END DO
						ELSE
							DO j = 1, N_part, 1
								DO i = 1, N_part, 1
									frf1=DEXP(-rijpc_ep_old(0,i,j)*Fep_yuk)
									O(1)=O(1)-Aep_yuk*frf1
								END DO
							END DO
						END IF
					END IF
				END SELECT
			END IF
		END IF
		
		O=O*0.5d0
		
	END SUBROUTINE derivata_Jep_YUK
!-----------------------------------------------------------------------
SUBROUTINE derivata_Jep_ATM(O)
	IMPLICIT NONE
	INTEGER :: i, j
	REAL (KIND=8) :: frf1
	REAL (KIND=8), DIMENSION(:) :: O

	O=0.d0
	
	SELECT CASE(Jep_kind)
	CASE('atm')
		DO i = 1, N_part, 1
			O=O-rij_ep_old(0,i,i)
		END DO
	CASE('atp')
		DO i = 1, N_part, 1
			O=O-rijpc_ep_old(0,i,i)
		END DO
	END SELECT
	
	O=O*0.5d0
	
END SUBROUTINE derivata_Jep_ATM
!-----------------------------------------------------------------------

	SUBROUTINE derivata_Jsese_YUK(O)
		IMPLICIT NONE
		INTEGER :: i, j
		REAL (KIND=8) :: frf1
		REAL (KIND=8), DIMENSION(:) :: O(:)
		!flag_Asese
		!	flag_Fsese			O(1) - Asese_yuk  ;  O(2) - Asese_ud_yuk  ;  O(3) - Fsese_yuk  ;  O(4) - Fsese_ud_yuk
		!	not flag_Fsese		O(1) - Asese_yuk  ;  O(2) - Asese_ud_yuk  ;  O(3) - Fsese_yuk
		!not flag_Asese
		!	flag_Fsese			O(1) - Asese_yuk  ;  O(2) - Fsese_yuk  ;  O(3) - Fsese_ud_yuk
		!	not flag_Fsese		O(1) - Asese_yuk  ;  O(2) - Fsese_yuk
		O=0.d0
		SELECT CASE(Jse_kind)
		CASE('yuk')
			IF (split_Asese) THEN
				IF (split_Fsese) THEN
					DO j = 1, N_part-1, 1
						DO i = j+1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_se1_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rij_se1_old(0,i,j)
							ELSE
								frf1=DEXP(-rij_se1_old(0,i,j)*Fsese_ud_yuk)
								O(4)=O(4)-Asese_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rij_se1_old(0,i,j)
							END IF
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_se2_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rij_se2_old(0,i,j)
							ELSE
								frf1=DEXP(-rij_se2_old(0,i,j)*Fsese_ud_yuk)
								O(4)=O(4)-Asese_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rij_se2_old(0,i,j)
							END IF
						END DO
					END DO
				ELSE
					DO j = 1, N_part-1, 1
						DO i = j+1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_se1_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rij_se1_old(0,i,j)
							ELSE
								frf1=DEXP(-rij_se1_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rij_se1_old(0,i,j)
							END IF
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_se2_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rij_se2_old(0,i,j)
							ELSE
								frf1=DEXP(-rij_se2_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rij_se2_old(0,i,j)
							END IF
						END DO
					END DO
				END IF
			ELSE
				IF (split_Fsese) THEN
					DO j = 1, N_part-1, 1
						DO i = j+1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_se1_old(0,i,j)*Fsese_yuk)
								O(2)=O(2)-Asese_yuk*frf1
							ELSE
								frf1=DEXP(-rij_se1_old(0,i,j)*Fsese_ud_yuk)
								O(3)=O(3)-Asese_yuk*frf1
							END IF
							O(1)=O(1)-(1.d0-frf1)/rij_se1_old(0,i,j)
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_se2_old(0,i,j)*Fsese_yuk)
								O(2)=O(2)-Asese_yuk*frf1
							ELSE
								frf1=DEXP(-rij_se2_old(0,i,j)*Fsese_ud_yuk)
								O(3)=O(3)-Asese_yuk*frf1
							END IF
							O(1)=O(1)-(1.d0-frf1)/rij_se2_old(0,i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N_part-1, 1
						DO i = j+1, N_part, 1
							frf1=DEXP(-rij_se1_old(0,i,j)*Fsese_yuk)
							O(2)=O(2)-Asese_yuk*frf1
							O(1)=O(1)-(1.d0-frf1)/rij_se1_old(0,i,j)
							frf1=DEXP(-rij_se2_old(0,i,j)*Fsese_yuk)
							O(2)=O(2)-Asese_yuk*frf1
							O(1)=O(1)-(1.d0-frf1)/rij_se2_old(0,i,j)
						END DO
					END DO
				END IF
			END IF
		CASE('yup')
			IF (split_Asese) THEN
				IF (split_Fsese) THEN
					DO j = 1, N_part-1, 1
						DO i = j+1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_se1_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rijpc_se1_old(0,i,j)
							ELSE
								frf1=DEXP(-rijpc_se1_old(0,i,j)*Fsese_ud_yuk)
								O(4)=O(4)-Asese_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rijpc_se1_old(0,i,j)
							END IF
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_se2_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rijpc_se2_old(0,i,j)
							ELSE
								frf1=DEXP(-rijpc_se2_old(0,i,j)*Fsese_ud_yuk)
								O(4)=O(4)-Asese_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rijpc_se2_old(0,i,j)
							END IF
						END DO
					END DO
				ELSE
					DO j = 1, N_part-1, 1
						DO i = j+1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_se1_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rijpc_se1_old(0,i,j)
							ELSE
								frf1=DEXP(-rijpc_se1_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rijpc_se1_old(0,i,j)
							END IF
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_se2_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rijpc_se2_old(0,i,j)
							ELSE
								frf1=DEXP(-rijpc_se2_old(0,i,j)*Fsese_yuk)
								O(3)=O(3)-Asese_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rijpc_se2_old(0,i,j)
							END IF
						END DO
					END DO
				END IF
			ELSE
				IF (split_Fsese) THEN
					DO j = 1, N_part-1, 1
						DO i = j+1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_se1_old(0,i,j)*Fsese_yuk)
								O(2)=O(2)-Asese_yuk*frf1
							ELSE
								frf1=DEXP(-rijpc_se1_old(0,i,j)*Fsese_ud_yuk)
								O(3)=O(3)-Asese_yuk*frf1
							END IF
							O(1)=O(1)-(1.d0-frf1)/rijpc_se1_old(0,i,j)
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_se2_old(0,i,j)*Fsese_yuk)
								O(2)=O(2)-Asese_yuk*frf1
							ELSE
								frf1=DEXP(-rijpc_se2_old(0,i,j)*Fsese_ud_yuk)
								O(3)=O(3)-Asese_yuk*frf1
							END IF
							O(1)=O(1)-(1.d0-frf1)/rijpc_se2_old(0,i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N_part-1, 1
						DO i = j+1, N_part, 1
							frf1=DEXP(-rijpc_se1_old(0,i,j)*Fsese_yuk)
							O(2)=O(2)-Asese_yuk*frf1
							O(1)=O(1)-(1.d0-frf1)/rijpc_se1_old(0,i,j)
							frf1=DEXP(-rijpc_se2_old(0,i,j)*Fsese_yuk)
							O(2)=O(2)-Asese_yuk*frf1
							O(1)=O(1)-(1.d0-frf1)/rijpc_se2_old(0,i,j)
						END DO
					END DO
				END IF
			END IF
		END SELECT
		
		O=O*0.25d0
		
	END SUBROUTINE derivata_Jsese_YUK
!-----------------------------------------------------------------------

	SUBROUTINE derivata_Jsese_POT(O)
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: A1=4.64355, B1=3.00909, A2=0.902629, B2=0.924515
		INTEGER :: i, j
		REAL (KIND=8) :: dist1, dist2, O(1:2)      !O(1) - A_POT_se, O(2) - D_POT_se
		O=0.d0
		DO i = 1, H_N_part, 1
			O(1)=O(1)-(A1*DEXP(-D_POT_se*B1*dist_mol_ss1_old(i))-A2*DEXP(-D_POT_se*B2*dist_mol_ss1_old(i)))
			O(2)=O(2)-A_POT_se*(-A1*B1*dist_mol_ss1_old(i)*DEXP(-D_POT_se*B1*dist_mol_ss1_old(i))+ &
			  A2*B2*dist_mol_ss1_old(i)*DEXP(-D_POT_se*B2*dist_mol_ss1_old(i)))
			O(1)=O(1)-(A1*DEXP(-D_POT_se*B1*dist_mol_ss2_old(i))-A2*DEXP(-D_POT_se*B2*dist_mol_ss2_old(i)))
			O(2)=O(2)-A_POT_se*(-A1*B1*dist_mol_ss2_old(i)*DEXP(-D_POT_se*B1*dist_mol_ss2_old(i))+ &
			  A2*B2*dist_mol_ss2_old(i)*DEXP(-D_POT_se*B2*dist_mol_ss2_old(i)))
		END DO
		O=O*0.5d0
	END SUBROUTINE derivata_Jsese_POT
!-----------------------------------------------------------------------

	SUBROUTINE derivata_KERNese(O)
		IMPLICIT NONE
		INTEGER :: i, j, i_sd, j_sd
		REAL (KIND=8) :: minL
		INTEGER :: pvt(1:H_N_part), info, perm
		REAL (KIND=8) :: SD(1:H_N_part,1:H_N_part), ISD(1:H_N_part,1:H_N_part), detSD
		REAL (KIND=8) :: O      !O(1) - C_kern_e
		O=0.d0
		
		minL=MINVAL(L)
		IF ( Kse_kind=='gss' ) THEN
			DO i = 1, N_part, 1
				O=O-rij_ese1_old(0,i,i)*rij_ese1_old(0,i,i)
			END DO
			DO i = 1, N_part, 1
				O=O-rij_ese2_old(0,i,i)*rij_ese2_old(0,i,i)
			END DO
		ELSE IF ( Kse_kind=='gsc' ) THEN
			DO i = 1, N_part, 1
				IF ( rij_ese1_old(0,i,i)<l1_kern_e ) THEN
					O=O-rij_ese1_old(0,i,i)*rij_ese1_old(0,i,i)
				ELSE IF ( rij_ese1_old(0,i,i)<minL*0.499 ) THEN
					!O=O-(alpha0_kern_e+alpha1_kern_e/(rij_ese1_old(0,i,i)-0.5d0*minL))/C_kern_e
					O = O + minL*minL/12.d0 + ((minL/3.d0)**3)/(rij_ese1_old(0,i,i)-0.5d0*minL)
				ELSE
					!O=O-(alpha0_kern_e+alpha1_kern_e/(-0.001d0*minL))/C_kern_e
					O = O + minL*minL/12.d0 + ((minL/3.d0)**3)/(-0.001d0*minL)
				END IF
			END DO
			DO i = 1, N_part, 1
				IF ( rij_ese2_old(0,i,i)<l1_kern_e ) THEN
					O=O-rij_ese2_old(0,i,i)*rij_ese2_old(0,i,i)
				ELSE IF ( rij_ese2_old(0,i,i)<minL*0.499 ) THEN
					!O=O-(alpha0_kern_e+alpha1_kern_e/(rij_ese2_old(0,i,i)-0.5d0*minL))/C_kern_e
					O = O + minL*minL/12.d0 + ((minL/3.d0)**3)/(rij_ese2_old(0,i,i)-0.5d0*minL)
				ELSE
					!O=O-(alpha0_kern_e+alpha1_kern_e/(-0.001d0*minL))/C_kern_e
					O = O + minL*minL/12.d0 + ((minL/3.d0)**3)/(-0.001d0*minL)
				END IF
			END DO
		ELSE IF ( Kse_kind=='gsp' ) THEN
			DO i = 1, N_part, 1
				O=O-rijpc_ese1_old(0,i,i)*rijpc_ese1_old(0,i,i)
			END DO
			DO i = 1, N_part, 1
				O=O-rijpc_ese2_old(0,i,i)*rijpc_ese2_old(0,i,i)
			END DO
		ELSE IF ( Kse_kind=='gdp') THEN
			!s1 up
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					SD(i,j)=-rijpc_ese1_old(0,i,j)*rijpc_ese1_old(0,i,j)* &
					  DEXP(-C_kern_e*rijpc_ese1_old(0,i,j)*rijpc_ese1_old(0,i,j))
					ISD(i,j)=SD(i,j)
				END DO
			END DO
			!Calcolo il determinante di SD_new
			CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
			IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
			perm=0
			detSD=1.d0
			DO  i = 1, H_N_part, 1
				IF (pvt(i) /= i) perm=perm+1
			END DO
			IF (MOD(perm,2) == 1 ) detSD=-detSD
			DO  i = 1, H_N_part, 1
				detSD=detSD*ISD(i,i)
			END DO
			O=O+detSD/detGDse1_up_old

			!s1 dw
			DO j = H_N_part+1, N_part, 1
				j_sd=j-H_N_part
				DO i = H_N_part+1, N_part, 1
					i_sd=i-H_N_part
					SD(i_sd,j_sd)=-rijpc_ese1_old(0,i,j)*rijpc_ese1_old(0,i,j)* &
					  DEXP(-C_kern_e*rijpc_ese1_old(0,i,j)*rijpc_ese1_old(0,i,j))
					ISD(i_sd,j_sd)=SD(i_sd,j_sd)
				END DO
			END DO
			!Calcolo il determinante di SD_new
			CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
			IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
			perm=0
			detSD=1.d0
			DO  i = 1, H_N_part, 1
				IF (pvt(i) /= i) perm=perm+1
			END DO
			IF (MOD(perm,2) == 1 ) detSD=-detSD
			DO  i = 1, H_N_part, 1
				detSD=detSD*ISD(i,i)
			END DO
			O=O+detSD/detGDse1_dw_old

			!s2 up
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					SD(i,j)=-rijpc_ese2_old(0,i,j)*rijpc_ese2_old(0,i,j)* &
					  DEXP(-C_kern_e*rijpc_ese2_old(0,i,j)*rijpc_ese2_old(0,i,j))
					ISD(i,j)=SD(i,j)
				END DO
			END DO
			!Calcolo il determinante di SD_new
			CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
			IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
			perm=0
			detSD=1.d0
			DO  i = 1, H_N_part, 1
				IF (pvt(i) /= i) perm=perm+1
			END DO
			IF (MOD(perm,2) == 1 ) detSD=-detSD
			DO  i = 1, H_N_part, 1
				detSD=detSD*ISD(i,i)
			END DO
			O=O+detSD/detGDse2_up_old

			!s2 dw
			DO j = H_N_part+1, N_part, 1
				j_sd=j-H_N_part
				DO i = H_N_part+1, N_part, 1
					i_sd=i-H_N_part
					SD(i_sd,j_sd)=-rijpc_ese2_old(0,i,j)*rijpc_ese2_old(0,i,j)* &
					  DEXP(-C_kern_e*rijpc_ese2_old(0,i,j)*rijpc_ese2_old(0,i,j))
					ISD(i_sd,j_sd)=SD(i_sd,j_sd)
				END DO
			END DO
			!Calcolo il determinante di SD_new
			CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
			IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
			perm=0
			detSD=1.d0
			DO  i = 1, H_N_part, 1
				IF (pvt(i) /= i) perm=perm+1
			END DO
			IF (MOD(perm,2) == 1 ) detSD=-detSD
			DO  i = 1, H_N_part, 1
				detSD=detSD*ISD(i,i)
			END DO
			O=O+detSD/detGDse2_dw_old
		END IF
				
		O=O*0.5d0
	END SUBROUTINE derivata_KERNese
!-----------------------------------------------------------------------

	SUBROUTINE derivata_atmKERNese(O)
		IMPLICIT NONE
		INTEGER :: i
		REAL (KIND=8) :: L38
		REAL (KIND=8) :: O      !O(1) - C_kern_e
		O=0.d0
		
		IF ( Kse_kind=='atm' ) THEN
			DO i = 1, N_part, 1
				O=O-rij_ese1_old(0,i,i)
			END DO
			DO i = 1, N_part, 1
				O=O-rij_ese2_old(0,i,i)
			END DO
		ELSE IF (Kse_kind=='atc') THEN
			L38=3.d0*MINVAL(L)/8.d0
			DO i = 1, N_part, 1
				IF ( rij_ese1_old(0,i,i)<L38 ) THEN
					O=O-rij_ese1_old(0,i,i)
				ELSE
					O=O-(alpha0_kern_e+alpha1_kern_e*((rij_ese1_old(0,i,i)-L38)**10))/C_kern_e
				END IF
			END DO
			DO i = 1, N_part, 1
				IF ( rij_ese2_old(0,i,i)<L38 ) THEN
					O=O-rij_ese2_old(0,i,i)
				ELSE
					O=O-(alpha0_kern_e+alpha1_kern_e*((rij_ese2_old(0,i,i)-L38)**10))/C_kern_e
				END IF
			END DO
		END IF
		
		O=O*0.5d0
	END SUBROUTINE derivata_atmKERNese
!-----------------------------------------------------------------------

	SUBROUTINE derivata_Jsesp_YUK(O)
		IMPLICIT NONE
		INTEGER :: i, j
		REAL (KIND=8) :: frf1
		REAL (KIND=8), DIMENSION(:) :: O(:)
		!flag_Asesp
		!	flag_Fsesp			O(1) - Asesp_yuk  ;  O(2) - Asesp_ud_yuk  ;  O(3) - Fsesp_yuk  ;  O(4) - Fsesp_ud_yuk
		!	not flag_Fsesp		O(1) - Asesp_yuk  ;  O(2) - Asesp_ud_yuk  ;  O(3) - Fsesp_yuk
		!not flag_Asesp
		!	flag_Fsesp			O(1) - Asesp_yuk  ;  O(2) - Fsesp_yuk  ;  O(3) - Fsesp_ud_yuk
		!	not flag_Fsesp		O(1) - Asesp_yuk  ;  O(2) - Fsesp_yuk
		O=0.d0
		SELECT CASE(Jsesp_kind)
		CASE('yuk')
			IF (split_Asesp) THEN
				IF (split_Fsesp) THEN
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_sesp1_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rij_sesp1_old(0,i,j)
							ELSE
								frf1=DEXP(-rij_sesp1_old(0,i,j)*Fsesp_ud_yuk)
								O(4)=O(4)-Asesp_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rij_sesp1_old(0,i,j)
							END IF
						END DO
					END DO
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_sesp2_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rij_sesp2_old(0,i,j)
							ELSE
								frf1=DEXP(-rij_sesp2_old(0,i,j)*Fsesp_ud_yuk)
								O(4)=O(4)-Asesp_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rij_sesp2_old(0,i,j)
							END IF
						END DO
					END DO
				ELSE
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_sesp1_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rij_sesp1_old(0,i,j)
							ELSE
								frf1=DEXP(-rij_sesp1_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rij_sesp1_old(0,i,j)
							END IF
						END DO
					END DO
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_sesp2_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rij_sesp2_old(0,i,j)
							ELSE
								frf1=DEXP(-rij_sesp2_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rij_sesp2_old(0,i,j)
							END IF
						END DO
					END DO
				END IF
			ELSE
				IF (split_Fsesp) THEN
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_sesp1_old(0,i,j)*Fsesp_yuk)
								O(2)=O(2)-Asesp_yuk*frf1
							ELSE
								frf1=DEXP(-rij_sesp1_old(0,i,j)*Fsesp_ud_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
							END IF
							O(1)=O(1)-(1.d0-frf1)/rij_sesp1_old(0,i,j)
						END DO
					END DO
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rij_sesp2_old(0,i,j)*Fsesp_yuk)
								O(2)=O(2)-Asesp_yuk*frf1
							ELSE
								frf1=DEXP(-rij_sesp2_old(0,i,j)*Fsesp_ud_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
							END IF
							O(1)=O(1)-(1.d0-frf1)/rij_sesp2_old(0,i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							frf1=DEXP(-rij_sesp1_old(0,i,j)*Fsesp_yuk)
							O(2)=O(2)-Asesp_yuk*frf1
							O(1)=O(1)-(1.d0-frf1)/rij_sesp1_old(0,i,j)
						END DO
					END DO
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							frf1=DEXP(-rij_sesp2_old(0,i,j)*Fsesp_yuk)
							O(2)=O(2)-Asesp_yuk*frf1
							O(1)=O(1)-(1.d0-frf1)/rij_sesp2_old(0,i,j)
						END DO
					END DO
				END IF
			END IF
		CASE('yup')
			IF (split_Asesp) THEN
				IF (split_Fsesp) THEN
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_sesp1_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rijpc_sesp1_old(0,i,j)
							ELSE
								frf1=DEXP(-rijpc_sesp1_old(0,i,j)*Fsesp_ud_yuk)
								O(4)=O(4)-Asesp_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rijpc_sesp1_old(0,i,j)
							END IF
						END DO
					END DO
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_sesp2_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rijpc_sesp2_old(0,i,j)
							ELSE
								frf1=DEXP(-rijpc_sesp2_old(0,i,j)*Fsesp_ud_yuk)
								O(4)=O(4)-Asesp_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rijpc_sesp2_old(0,i,j)
							END IF
						END DO
					END DO
				ELSE
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_sesp1_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rijpc_sesp1_old(0,i,j)
							ELSE
								frf1=DEXP(-rijpc_sesp1_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rijpc_sesp1_old(0,i,j)
							END IF
						END DO
					END DO
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_sesp2_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
								O(1)=O(1)-(1.d0-frf1)/rijpc_sesp2_old(0,i,j)
							ELSE
								frf1=DEXP(-rijpc_sesp2_old(0,i,j)*Fsesp_yuk)
								O(3)=O(3)-Asesp_ud_yuk*frf1
								O(2)=O(2)-(1.d0-frf1)/rijpc_sesp2_old(0,i,j)
							END IF
						END DO
					END DO
				END IF
			ELSE
				IF (split_Fsesp) THEN
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_sesp1_old(0,i,j)*Fsesp_yuk)
								O(2)=O(2)-Asesp_yuk*frf1
							ELSE
								frf1=DEXP(-rijpc_sesp1_old(0,i,j)*Fsesp_ud_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
							END IF
							O(1)=O(1)-(1.d0-frf1)/rijpc_sesp1_old(0,i,j)
						END DO
					END DO
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
								frf1=DEXP(-rijpc_sesp2_old(0,i,j)*Fsesp_yuk)
								O(2)=O(2)-Asesp_yuk*frf1
							ELSE
								frf1=DEXP(-rijpc_sesp2_old(0,i,j)*Fsesp_ud_yuk)
								O(3)=O(3)-Asesp_yuk*frf1
							END IF
							O(1)=O(1)-(1.d0-frf1)/rijpc_sesp2_old(0,i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							frf1=DEXP(-rijpc_sesp1_old(0,i,j)*Fsesp_yuk)
							O(2)=O(2)-Asesp_yuk*frf1
							O(1)=O(1)-(1.d0-frf1)/rijpc_sesp1_old(0,i,j)
						END DO
					END DO
					DO j = 1, N_part, 1
						DO i = 1, N_part, 1
							frf1=DEXP(-rijpc_sesp2_old(0,i,j)*Fsesp_yuk)
							O(2)=O(2)-Asesp_yuk*frf1
							O(1)=O(1)-(1.d0-frf1)/rijpc_sesp2_old(0,i,j)
						END DO
					END DO
				END IF
			END IF
		END SELECT
		
		O=O*0.25d0
		
	END SUBROUTINE derivata_Jsesp_YUK
!-----------------------------------------------------------------------

	SUBROUTINE derivata_KERNsesp(O)
		IMPLICIT NONE
		INTEGER :: i
		REAL (KIND=8) :: O      !O(1) - Gsesp
		O=0.d0
		DO i = 1, N_part, 1
			O=O-rij_sesp1_old(0,i,i)*rij_sesp1_old(0,i,i)
		END DO
		DO i = 1, N_part, 1
			O=O-rij_sesp2_old(0,i,i)*rij_sesp2_old(0,i,i)
		END DO
		O=O*0.5d0
	END SUBROUTINE derivata_KERNsesp
!-----------------------------------------------------------------------

	SUBROUTINE derivata_SDse_atm(O)
		USE walkers
		IMPLICIT NONE
		INTEGER :: i, j, i_sd, j_sd
		REAL (KIND=8) :: O      !O(1) - Gsesp
		INTEGER :: pvt(1:H_N_part), info, perm
		REAL (KIND=8) :: SD(1:H_N_part,1:H_N_part), ISD(1:H_N_part,1:H_N_part), detSD
		O=0.d0
		
		!s1 up
		DO j = 1, H_N_part, 1
			DO i = 1, H_N_part, 1
				SD(i,j)=-rij_sesp1_old(0,i,j)*DEXP(-C_atm*rij_sesp1_old(0,i,j))
				ISD(i,j)=SD(i,j)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDse1_up_old,8)
		
		!s1 dw
		DO j = H_N_part+1, N_part, 1
			j_sd=j-H_N_part
			DO i = H_N_part+1, N_part, 1
				i_sd=i-H_N_part
				SD(i_sd,j_sd)=-rij_sesp1_old(0,i,j)*DEXP(-C_atm*rij_sesp1_old(0,i,j))
				ISD(i_sd,j_sd)=SD(i_sd,j_sd)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDse1_dw_old,8)
		
		!s2 up
		DO j = 1, H_N_part, 1
			DO i = 1, H_N_part, 1
				SD(i,j)=-rij_sesp2_old(0,i,j)*DEXP(-C_atm*rij_sesp2_old(0,i,j))
				ISD(i,j)=SD(i,j)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDse2_up_old,8)
		
		!s2 dw
		DO j = H_N_part+1, N_part, 1
			j_sd=j-H_N_part
			DO i = H_N_part+1, N_part, 1
				i_sd=i-H_N_part
				SD(i_sd,j_sd)=-rij_sesp2_old(0,i,j)*DEXP(-C_atm*rij_sesp2_old(0,i,j))
				ISD(i_sd,j_sd)=SD(i_sd,j_sd)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDse2_dw_old,8)
		
		O=O*0.5d0
	END SUBROUTINE derivata_SDse_atm
!-----------------------------------------------------------------------

	SUBROUTINE derivata_SDse_atp(O)
		USE walkers
		IMPLICIT NONE
		INTEGER :: i, j, i_sd, j_sd
		REAL (KIND=8) :: O      !O(1) - Gsesp
		INTEGER :: pvt(1:H_N_part), info, perm
		REAL (KIND=8) :: SD(1:H_N_part,1:H_N_part), ISD(1:H_N_part,1:H_N_part), detSD
		O=0.d0
		
		!s1 up
		DO j = 1, H_N_part, 1
			DO i = 1, H_N_part, 1
				SD(i,j)=-rijpc_sesp1_old(0,i,j)*DEXP(-C_atm*rijpc_sesp1_old(0,i,j))
				ISD(i,j)=SD(i,j)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDse1_up_old,8)
		
		!s1 dw
		DO j = H_N_part+1, N_part, 1
			j_sd=j-H_N_part
			DO i = H_N_part+1, N_part, 1
				i_sd=i-H_N_part
				SD(i_sd,j_sd)=-rijpc_sesp1_old(0,i,j)*DEXP(-C_atm*rijpc_sesp1_old(0,i,j))
				ISD(i_sd,j_sd)=SD(i_sd,j_sd)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDse1_dw_old,8)
		
		!s2 up
		DO j = 1, H_N_part, 1
			DO i = 1, H_N_part, 1
				SD(i,j)=-rijpc_sesp2_old(0,i,j)*DEXP(-C_atm*rijpc_sesp2_old(0,i,j))
				ISD(i,j)=SD(i,j)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDse2_up_old,8)
		
		!s2 dw
		DO j = H_N_part+1, N_part, 1
			j_sd=j-H_N_part
			DO i = H_N_part+1, N_part, 1
				i_sd=i-H_N_part
				SD(i_sd,j_sd)=-rijpc_sesp2_old(0,i,j)*DEXP(-C_atm*rijpc_sesp2_old(0,i,j))
				ISD(i_sd,j_sd)=SD(i_sd,j_sd)
			END DO
		END DO
		!Calcolo il determinante di SD_new
		CALL DGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		perm=0
		detSD=1.d0
		DO  i = 1, H_N_part, 1
			IF (pvt(i) /= i) perm=perm+1
		END DO
		IF (MOD(perm,2) == 1 ) detSD=-detSD
		DO  i = 1, H_N_part, 1
			detSD=detSD*ISD(i,i)
		END DO
		O=O+detSD/REAL(detSDse2_dw_old,8)
		
		O=O*0.5d0
	END SUBROUTINE derivata_SDse_atp
!-----------------------------------------------------------------------

	SUBROUTINE derivata_psi_Rp(O)
		IMPLICIT NONE
		REAL (KIND=8), DIMENSION(:) :: O(:)
		INTEGER :: i, j
		REAL (KIND=8) :: uep1, frf1(1:3), frf2(1:3), frf3
		
		SELECT CASE (Jep_kind)
		CASE ('yuk')
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
					IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
						uep1=Aep_yuk*( -(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_yuk* &
						  DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
					ELSE
						IF (split_Fee) THEN
							IF (split_Aee) THEN
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_ud_yuk* &
								  DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
							ELSE
								uep1=Aep_yuk*( -(1.d0-DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_ud_yuk* &
								  DEXP(-Fep_ud_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
							END IF
						ELSE
							IF (split_Aee) THEN
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_yuk* &
								  DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
							ELSE
								uep1=Aep_yuk*( -(1.d0-DEXP(-Fep_yuk*rij_ep_old(0,i,j)))/rij_ep_old(0,i,j)+Fep_yuk* &
								  DEXP(-Fep_yuk*rij_ep_old(0,i,j)) )/rij_ep_old(0,i,j)
							END IF
						END IF
					END IF
					O(3*(j-1)+1:3*j)=O(3*(j-1)+1:3*j)-0.5d0*uep1*rij_ep_old(1:3,i,j)/rij_ep_old(0,i,j)
				END DO
			END DO
		CASE ('yup')
			frf1(1:3)=PI/L(1:3)
			DO j = 1, N_part, 1
				DO i = 1, N_part, 1
					frf2(1:3)=(/ rijpc_ep_old(1,i,j)*DCOS(frf1(1)*rij_ep_old(1,i,j)), &
					  rijpc_ep_old(2,i,j)*DCOS(frf1(2)*rij_ep_old(2,i,j)), &
					  rijpc_ep_old(3,i,j)*DCOS(frf1(3)*rij_ep_old(3,i,j)) /)/rijpc_ep_old(0,i,j)
					
					IF ( (i<=H_N_part .AND. j<=H_N_part) .OR. (i>H_N_part .AND. j>H_N_part) ) THEN
						uep1=Aep_yuk*( -(1.d0-DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)))/rijpc_ep_old(0,i,j)+Fep_yuk* &
						  DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
					ELSE
						IF (split_Fee) THEN
							IF (split_Aee) THEN
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)))/rijpc_ep_old(0,i,j)+Fep_ud_yuk* &
								  DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
							ELSE
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)))/rijpc_ep_old(0,i,j)+Fep_ud_yuk* &
								  DEXP(-Fep_ud_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
							END IF
						ELSE
							IF (split_Aee) THEN
								uep1=Aep_ud_yuk*( -(1.d0-DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)))/rijpc_ep_old(0,i,j)+Fep_yuk* &
								  DEXP(-Fep_yuk*rijpc_ep_old(0,i,j)) )/rijpc_ep_old(0,i,j)
							ELSE
								frf3=DEXP(-Fep_yuk*rijpc_ep_old(0,i,j))
								uep1=Aep_yuk*( -(1.d0-frf3)/rijpc_ep_old(0,i,j)+Fep_yuk*frf3)/rijpc_ep_old(0,i,j)
							END IF
						END IF
					END IF
					O(3*(j-1)+1:3*j)=O(3*(j-1)+1:3*j)-0.5d0*uep1*frf2(1:3)
				END DO
			END DO
		END SELECT
	END SUBROUTINE derivata_psi_Rp
!-----------------------------------------------------------------------

	SUBROUTINE chiudi_estimatori()
		IMPLICIT NONE
		IF (flag_E_tot) DEALLOCATE(E_tot)
		IF (flag_E_kin) DEALLOCATE(E_kin, E_JF)
		IF (flag_E_pot) DEALLOCATE(E_pot)
		DEALLOCATE(w)
		iniz_estimatori=.FALSE.
	END SUBROUTINE chiudi_estimatori
	
END MODULE estimatori
