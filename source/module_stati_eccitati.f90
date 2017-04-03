MODULE stati_eccitati
	USE dati_mc
	USE dati_fisici
	USE funzione_onda
	IMPLICIT NONE
	CHARACTER(LEN=31), PROTECTED, SAVE :: file_exc_stat='orbitals/stati_eccitati.xxxx'
	LOGICAL, PROTECTED, SAVE :: iniz_stati_Eccitati=.FALSE.
	LOGICAL, PROTECTED, SAVE, ALLOCATABLE :: iniz_excwf(:)
	TYPE exc_wave_function
		COMPLEX (KIND=8), ALLOCATABLE :: SDe_up_new(:,:), SDe_up_old(:,:), ISDe_up_new(:,:), ISDe_up_old(:,:)
		INTEGER, ALLOCATABLE :: pvte_up_new(:), pvte_up_old(:)
		COMPLEX (KIND=8), ALLOCATABLE :: SDe_dw_new(:,:), SDe_dw_old(:,:), ISDe_dw_new(:,:), ISDe_dw_old(:,:)
		INTEGER, ALLOCATABLE :: pvte_dw_new(:), pvte_dw_old(:)
		COMPLEX (KIND=8) ::  detSDe_up_new, detSDe_up_old, detSDe_dw_new, detSDe_dw_old
		INTEGER, ALLOCATABLE :: num_pw_dnfH(:)
		REAL (KIND=8), ALLOCATABLE :: k_pw_dnfH(:,:,:), fattori_pw_dnfH(:,:)
	END TYPE exc_wave_function
	TYPE (exc_wave_function), ALLOCATABLE, SAVE, PROTECTED :: excwf(:)
	INTEGER, SAVE :: N_exc, N_exc_stat, N_orb
	INTEGER, SAVE, ALLOCATABLE :: i_exc_stat(:,:)
	REAL (KIND=8) :: E0, E_ctf
	REAL (KIND=8), SAVE, ALLOCATABLE :: energie_exc_stat(:)
	!N_exc contiene il numero di stati eccitati da considerare, N_exc_stat quelli consideraibili con il cutoff
	
	CONTAINS
	
	SUBROUTINE inizializza_stati_eccitati()
		IMPLICIT NONE
		CHARACTER(LEN=4) :: istring
		REAL (KIND=8) :: cutoff    !fattore da moltiplicare a E0 per ottenere l'energia massima da considerare
		NAMELIST /dati_conduttivita/ cutoff, N_exc
		
		OPEN (2, FILE='dati_conduttivita.d',STATUS='OLD')
		READ (2, NML=dati_conduttivita)
		CLOSE (2)
		
		IF (.NOT. iniz_funzione_onda) STOP 'non puoi trovare gli stati eccitati se non hai inizializzato la funzione d onda [ module_stati_eccitati.f90 > inizializza_stati_eccitati ]'
		
		WRITE (istring, '(I4.4)') mpi_myrank
		file_exc_stat(28:31)=istring
		
		E0=SUM(autoenergie_dnfH(1:H_N_part))
		E_ctf=E0+cutoff*DABS(E0)
		N_orb=N_pw_lda
		CALL trova_stati_eccitati()
		
		ALLOCATE(excwf(1:N_exc),iniz_excwf(1:N_exc))
		iniz_excwf=.FALSE.
				
		iniz_stati_eccitati=.TRUE.
		
	END SUBROUTINE inizializza_stati_eccitati
!-----------------------------------------------------------------------
	SUBROUTINE trova_stati_eccitati()
		IMPLICIT NONE
		INTEGER :: i1, i2, num
		INTEGER :: dummy_i(1:H_N_part)
		INTEGER, ALLOCATABLE :: big_i_exc_stat(:,:)
		REAL (KIND=8) :: dummy_r, min_E
		REAL (KIND=8), ALLOCATABLE :: big_energie_exc_stat(:)
		
		N_exc_stat=0
		i1=0
		OPEN (UNIT=99, FILE=file_exc_stat, STATUS='UNKNOWN')
		DO WHILE (num>0)
			CALL trova_stati_eccitati_Np_Nh(i1,num)
			i1=i1+1
			N_exc_stat=N_exc_stat+num
		END DO
		N_exc_stat=N_exc_stat-1
		CLOSE (99)
		IF (mpi_myrank==0) THEN
			PRINT * , 'EXC_STAT: Numero stati eccitati: ', N_exc_stat
			IF (flag_output) THEN
				WRITE (7, *) 'EXC_STAT: Numero stati eccitati: ', N_exc_stat
			END IF
		END IF
		
		ALLOCATE(big_energie_exc_stat(0:N_exc_stat),big_i_exc_stat(1:H_N_part,0:N_exc_stat))
		OPEN (UNIT=99, FILE=file_exc_stat, STATUS='OLD')
		DO i1 = 0, N_exc_stat, 1
			READ (99,*) big_energie_exc_stat(i1), big_i_exc_stat(1:H_N_part,i1)
		END DO
		CLOSE(99)
		CALL SYSTEM('rm -f '//file_exc_stat)
		
		IF (N_exc==-1.d0) N_exc=N_exc_stat
		N_exc=MIN(N_exc_stat,N_exc)
		N_exc=N_exc
		DO i1 = 0, N_exc, 1
			min_E=big_energie_exc_stat(i1)
			DO i2 = i1, N_exc_stat, 1
				IF (big_energie_exc_stat(i2)<min_E) THEN
					min_E=big_energie_exc_stat(i2)
					dummy_i(1:H_N_part)=big_i_exc_stat(1:H_N_part,i2)
					dummy_r=big_energie_exc_stat(i2)
					big_i_exc_stat(1:H_N_part,i2)=big_i_exc_stat(1:H_N_part,i1)
					big_energie_exc_stat(i2)=big_energie_exc_stat(i1)
					big_i_exc_stat(1:H_N_part,i1)=dummy_i(1:H_N_part)
					big_energie_exc_stat(i1)=dummy_r
				END IF
			END DO
		END DO
		
		ALLOCATE(energie_exc_stat(0:N_exc),i_exc_stat(1:H_N_part,0:N_exc))
		energie_exc_stat(0:N_exc)=big_energie_exc_stat(0:N_exc)
		i_exc_stat(1:H_N_part,0:N_exc)=big_i_exc_stat(1:H_N_part,0:N_exc)
		DEALLOCATE(big_energie_exc_stat,big_i_exc_stat)
				
	END SUBROUTINE trova_stati_eccitati
!-----------------------------------------------------------------------
	SUBROUTINE trova_stati_eccitati_Np_Nh(Np_Nh,num)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: Np_Nh
		LOGICAL :: flag_loop_h, flag_loop_p, flag_ext
		INTEGER :: stati_occ(1:H_N_part), ii(1:Np_Nh), i1, i2, i3
		REAL (KIND=8) :: E_stato
		INTEGER, INTENT(OUT) :: num
		
		stati_occ(1:H_N_part-Np_Nh)=(/(i1,i1=1,H_N_part-Np_Nh)/)
		stati_occ(H_N_part-Np_Nh+1:H_N_part)=(/(i1,i1=H_N_part+1,H_N_part+Np_Nh)/)
		
		num=0
		
		flag_loop_h=.TRUE.
		DO WHILE (flag_loop_h)
			flag_loop_p=.TRUE.
			DO WHILE (flag_loop_p)
				E_stato=SUM((/(autoenergie_dnfH(stati_occ(i3)),i3=1,H_N_part)/))
				IF (E_stato<=E_ctf) THEN
					num=num+1
					!PRINT * , E_stato, stati_occ
					WRITE (99, *) E_stato, stati_occ
				END IF
				IF (Np_Nh>0) THEN
					flag_ext=.FALSE.
					IF (stati_occ(H_N_part)+1>N_orb) THEN
						flag_ext=.TRUE.
					ELSE
						stati_occ(H_N_part)=stati_occ(H_N_part)+1
					END IF
					i1=H_N_part
					DO WHILE (flag_ext)
						!i1 é l'indice della p che é fuori range
						IF (i1-1<=H_N_part-Np_Nh) flag_loop_p=.FALSE.
						IF (flag_loop_p) THEN
							IF (stati_occ(i1-1)+1<=N_orb-(H_N_part-i1+1)) THEN
								stati_occ(i1-1)=stati_occ(i1-1)+1
								DO i2 = i1, H_N_part, 1
									stati_occ(i2)=stati_occ(i2-1)+1
								END DO
								IF (stati_occ(H_N_part)<=N_orb) THEN
									flag_ext=.FALSE.
								ELSE
									i1=i1-1
								END IF
							ELSE
								flag_ext=.FALSE.
								flag_loop_p=.FALSE.
							END IF
						ELSE
							flag_ext=.FALSE.
						END IF
					END DO
				ELSE
					flag_loop_p=.FALSE.
				END IF
			END DO
								
			flag_ext=.FALSE.
			IF (stati_occ(H_N_part-Np_Nh)+1>H_N_part) THEN
				flag_ext=.TRUE.
			ELSE
				stati_occ(H_N_part-Np_Nh)=stati_occ(H_N_part-Np_Nh)+1
			END IF
			i1=H_N_part-Np_Nh
			DO WHILE (flag_ext)
				IF (i1-1<=0) flag_loop_h=.FALSE.
				IF (flag_loop_h) THEN
					IF (stati_occ(i1-1)+1<=H_N_part-(H_N_part-Np_Nh-i1+1)) THEN
						stati_occ(i1-1)=stati_occ(i1-1)+1
						DO i2 = i1, H_N_part-Np_Nh, 1
							stati_occ(i2)=stati_occ(i2-1)+1
						END DO
						IF (stati_occ(H_N_part-Np_Nh)<=H_N_part) THEN
							flag_ext=.FALSE.
						ELSE
							i1=i1-1
						END IF
					ELSE
						i1=i1-1
					END IF
				ELSE
					flag_ext=.FALSE.
				END IF 
			END DO
						
			stati_occ(H_N_part-Np_Nh+1:H_N_part)=(/(i1,i1=H_N_part+1,H_N_part+Np_Nh)/)
		END DO
		
	END SUBROUTINE trova_stati_eccitati_Np_Nh
!-----------------------------------------------------------------------
	SUBROUTINE load_stato_eccitato(n_stat)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: n_stat
		INTEGER :: i, j, j1
		
		IF (.NOT. iniz_stati_eccitati) STOP 'non puoi caricare uno stato eccitato se prima non hai inizializzato il modulo [ module_stati_eccitati.f90 > load_stato_eccitato ]'
		IF (n_stat>N_exc) STOP 'stai cercando di caricare uno stato di energia troppo alta [ module_stati_eccitati.f90 > load_stato_eccitato ]'
		
		IF (mpi_myrank==0) THEN
			PRINT * , 'EXC_STAT: Calcolo dello stato eccitato ', n_stat ,' con energia ',  &
			  energie_exc_stat(n_stat), ', usando gli orbitali: ' 
			PRINT * , i_exc_stat(1:H_N_part,n_stat)
			IF (flag_output) THEN
				WRITE (7, *) 'EXC_STAT: Calcolo dello stato eccitato con energia ',  &
				  energie_exc_stat(n_stat), ', usando gli orbitali: '
				WRITE (7, *) i_exc_stat(1:H_N_part,n_stat)
			END IF
		END IF
		
		DO j = 1, H_N_part, 1
			num_pw_dnfH(j)=num_pw_orbit(i_exc_stat(j,n_stat))
		END DO
		DO j = 1, H_N_part, 1
			j1=i_exc_stat(j,n_stat)
			DO i = 1, num_pw_dnfH(j), 1
				k_pw_dnfH(0:3,i,j)=k_pw_lda(0:3,indice_pw_dnfH(i,j1))
			END DO
		END DO
		DO j = 1, H_N_part, 1
			j1=i_exc_stat(j,n_stat)
			DO i = 1, num_pw_dnfH(j), 1
				fattori_pw_dnfH(i,j)=fattori_pw_lda(indice_pw_dnfH(i,j1),j1)
			END DO
		END DO
		
	END SUBROUTINE load_stato_eccitato
!-----------------------------------------------------------------------
	SUBROUTINE inizializza_funzione_onda_eccitata(iexc)
		USE walkers
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: iexc
		INTEGER :: i1, i2, i3, max_num_pw
		
		IF (.NOT. iniz_stati_eccitati) STOP 'devi prima inizializzare [ module_stati_eccitati.f90 > inizializza_funzioni_onda_eccitate ]'
				
		ALLOCATE(excwf(iexc)%SDe_up_new(1:H_N_part,1:H_N_part),excwf(iexc)%SDe_up_old(1:H_N_part,1:H_N_part))
		excwf(iexc)%SDe_up_new=0.d0
		excwf(iexc)%SDe_up_old=0.d0
		ALLOCATE(excwf(iexc)%ISDe_up_new(1:H_N_part,1:H_N_part),excwf(iexc)%ISDe_up_old(1:H_N_part,1:H_N_part))
		excwf(iexc)%ISDe_up_new=0.d0
		excwf(iexc)%ISDe_up_old=0.d0
		ALLOCATE(excwf(iexc)%pvte_up_new(1:H_N_part),excwf(iexc)%pvte_up_old(1:H_N_part))
		excwf(iexc)%pvte_up_new=0.d0
		excwf(iexc)%pvte_up_old=0.d0
		ALLOCATE(excwf(iexc)%SDe_dw_new(1:H_N_part,1:H_N_part),excwf(iexc)%SDe_dw_old(1:H_N_part,1:H_N_part))
		excwf(iexc)%SDe_dw_new=0.d0
		excwf(iexc)%SDe_dw_old=0.d0
		ALLOCATE(excwf(iexc)%ISDe_dw_new(1:H_N_part,1:H_N_part),excwf(iexc)%ISDe_dw_old(1:H_N_part,1:H_N_part))
		excwf(iexc)%ISDe_dw_new=0.d0
		excwf(iexc)%ISDe_dw_old=0.d0
		ALLOCATE(excwf(iexc)%pvte_dw_new(1:H_N_part),excwf(iexc)%pvte_dw_old(1:H_N_part))
		excwf(iexc)%pvte_dw_new=0.d0
		excwf(iexc)%pvte_dw_old=0.d0
		
		IF ((SDe_kind=='prf').OR.(SDe_kind=='fre')) THEN
			ALLOCATE(excwf(iexc)%num_pw_dnfH(H_N_part))
			DO i1 = 1, H_N_part, 1
				excwf(iexc)%num_pw_dnfH(i1)=num_pw_orbit(i_exc_stat(i1,iexc))
			END DO
			max_num_pw=MAXVAL(excwf(iexc)%num_pw_dnfH(1:H_N_part))
			
			ALLOCATE(excwf(iexc)%k_pw_dnfH(0:3,max_num_pw,H_N_part))
			DO i3 = 1, H_N_part, 1
				i2=i_exc_stat(i3,iexc)
				DO i1 = 1, excwf(iexc)%num_pw_dnfH(i3), 1
					excwf(iexc)%k_pw_dnfH(0:3,i1,i3)=k_pw_lda(0:3,indice_pw_dnfH(i1,i2))
				END DO
			END DO
			
			ALLOCATE(excwf(iexc)%fattori_pw_dnfH(max_num_pw,H_N_part))
			DO i3 = 1, H_N_part, 1
				i2=i_exc_stat(i3,iexc)
				DO i1 = 1, excwf(iexc)%num_pw_dnfH(i3), 1
					excwf(iexc)%fattori_pw_dnfH(i1,i3)=fattori_pw_lda(indice_pw_dnfH(i1,i2),i2)
				END DO
			END DO
		END IF
		
		CALL prima_valutazione_funzione_onda_eccitata(iexc)
		iniz_excwf(iexc)=.TRUE.
		
	END SUBROUTINE inizializza_funzione_onda_eccitata
!-----------------------------------------------------------------------
	SUBROUTINE prima_valutazione_funzione_onda_eccitata(iexc)
		USE walkers
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: iexc
		INTEGER :: info
		COMPLEX (KIND=8) :: work(1:3*H_N_part)
		
		SELECT CASE (SDe_kind)
		CASE ('pw_') 
        
		CASE ('lda') 
        
		CASE ('prf')
			CALL valuta_SD_har(-1,re_old(1:3,1:H_N_part),H_N_part, &
			  excwf(iexc)%SDe_up_old,excwf(iexc)%detSDe_up_old,excwf(iexc)%ISDe_up_old, &
			  excwf(iexc)%pvte_up_old,excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new)
			CALL valuta_SD_har(-1,re_old(1:3,H_N_part+1:N_part),H_N_part, &
			  excwf(iexc)%SDe_dw_old,excwf(iexc)%detSDe_dw_old,excwf(iexc)%ISDe_dw_old, &
			  excwf(iexc)%pvte_dw_old,excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new)
			CALL ZGETRI( H_N_part, excwf(iexc)%ISDe_up_old, H_N_part, excwf(iexc)%pvte_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'EXC: ERRORE NEL TROVARE LA MATRICE INVERSA R UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, excwf(iexc)%ISDe_dw_old, H_N_part, excwf(iexc)%pvte_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'EXC: ERRORE NEL TROVARE LA MATRICE INVERSA R DOWN. INFO=', info
				STOP
			END IF
		CASE ('fre')
			CALL valuta_SD_har(-1,re_old(1:3,1:H_N_part),H_N_part, &
			  excwf(iexc)%SDe_up_old,excwf(iexc)%detSDe_up_old,excwf(iexc)%ISDe_up_old, &
			  excwf(iexc)%pvte_up_old,excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new)
			CALL valuta_SD_har(-1,re_old(1:3,H_N_part+1:N_part),H_N_part, &
			  excwf(iexc)%SDe_dw_old,excwf(iexc)%detSDe_dw_old,excwf(iexc)%ISDe_dw_old, &
			  excwf(iexc)%pvte_dw_old,excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new)
			CALL ZGETRI( H_N_part, excwf(iexc)%ISDe_up_old, H_N_part, excwf(iexc)%pvte_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'EXC: ERRORE NEL TROVARE LA MATRICE INVERSA R UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, excwf(iexc)%ISDe_dw_old, H_N_part, excwf(iexc)%pvte_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'EXC: ERRORE NEL TROVARE LA MATRICE INVERSA R DOWN. INFO=', info
				STOP
			END IF
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di SDe_kind accettabile [ module_stati_eccitati.f90 > prima_valutazione_funzioni_onda_eccitate ]'
		END SELECT
		IF (SDe_kind/='no_') THEN
			excwf(iexc)%SDe_up_new=excwf(iexc)%SDe_up_old
			excwf(iexc)%ISDe_up_new=excwf(iexc)%ISDe_up_old
			excwf(iexc)%SDe_dw_new=excwf(iexc)%SDe_dw_old
			excwf(iexc)%ISDe_dw_new=excwf(iexc)%ISDe_dw_old
		END IF
		excwf(iexc)%detSDe_up_new=excwf(iexc)%detSDe_up_old
		excwf(iexc)%detSDe_dw_new=excwf(iexc)%detSDe_dw_old
		
	END SUBROUTINE prima_valutazione_funzione_onda_eccitata
!-----------------------------------------------------------------------
	SUBROUTINE ricalcola_funzione_onda_eccitata(iexc,tipo,num)
		USE walkers
		IMPLICIT NONE
		CHARACTER(LEN=3) :: tipo
		INTEGER, INTENT(IN) :: iexc, num
		
		SELECT CASE (howtomove)
		CASE ('allp')
			!SDe
			IF ((tipo=='all') .OR. (tipo=='e__')) THEN
				SELECT CASE (SDe_kind)
				CASE ('pw_')   !ATTENZIONE, pw da sistemare!
					
				CASE ('lda')
					CALL valuta_SD_lda(num,re_new(1:3,1:H_N_part),H_N_part, &
					  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
					  excwf(iexc)%pvte_up_new,excwf(iexc)%SDe_up_old,excwf(iexc)%detSDe_up_old)
					CALL valuta_SD_lda(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
					  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
					  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
				CASE ('prf')
					CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
					  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
					  excwf(iexc)%pvte_up_new,excwf(iexc)%SDe_up_old,excwf(iexc)%detSDe_up_old)
					CALL valuta_SD_har(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
					  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
					  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
				CASE ('fre')
					CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
					  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
					  excwf(iexc)%pvte_up_new,excwf(iexc)%SDe_up_old,excwf(iexc)%detSDe_up_old)
					CALL valuta_SD_har(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
					  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
					  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
				CASE ('no_')
					excwf(iexc)%detSDe_up_new=1.d0
					excwf(iexc)%detSDe_dw_new=1.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di SDe_kind accettabile [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
		CASE ('1ppt')
			!SDe
			IF (tipo=='e__') THEN
				SELECT CASE (SDe_kind)
				CASE ('pw_')   !attenzione, ps da sistemare
					
				CASE ('lda')
					IF (num==-1) THEN
						CALL valuta_SD_lda(num,re_new(1:3,1:H_N_part),H_N_part, &
						  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
						  excwf(iexc)%pvte_up_new,excwf(iexc)%SDe_up_old,excwf(iexc)%detSDe_up_old)
						CALL valuta_SD_lda(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
						  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_lda(num,re_new(1:3,1:H_N_part),H_N_part, &
						  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
						  excwf(iexc)%pvte_up_new,excwf(iexc)%ISDe_up_old,excwf(iexc)%detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_lda(num-H_N_part,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
						  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
					END IF
				CASE ('prf')
					IF (num==-1) THEN
						CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
						  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
						  excwf(iexc)%pvte_up_new,excwf(iexc)%SDe_up_old,excwf(iexc)%detSDe_up_old)
						CALL valuta_SD_har(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
						  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
						  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
						  excwf(iexc)%pvte_up_new,excwf(iexc)%ISDe_up_old,excwf(iexc)%detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_har(num-H_N_part,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
						  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
					END IF
				CASE ('fre')
					IF (num==-1) THEN
						CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
						  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
						  excwf(iexc)%pvte_up_new,excwf(iexc)%SDe_up_old,excwf(iexc)%detSDe_up_old)
						CALL valuta_SD_har(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
						  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
						  excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new,excwf(iexc)%ISDe_up_new, &
						  excwf(iexc)%pvte_up_new,excwf(iexc)%ISDe_up_old,excwf(iexc)%detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_har(num-H_N_part,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new,excwf(iexc)%ISDe_dw_new, &
						  excwf(iexc)%pvte_dw_new,excwf(iexc)%ISDe_dw_old,excwf(iexc)%detSDe_dw_old)
					END IF
				CASE ('no_')
					excwf(iexc)%detSDe_up_new=1.d0
					excwf(iexc)%detSDe_dw_new=1.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di SDe_kind accettabile [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
		END SELECT
		
	END SUBROUTINE ricalcola_funzione_onda_eccitata
!-----------------------------------------------------------------------
	SUBROUTINE roll_back_exc(iexc,tipo,num)
		IMPLICIT NONE
		CHARACTER(LEN=3) :: tipo
		INTEGER, INTENT(IN) :: iexc, num
		
		IF (num==-1) THEN
			!SDe
			IF ((tipo=='all') .OR. (tipo=='e__')) THEN
				IF (SDe_kind/='no_') THEN
					excwf(iexc)%SDe_up_new=excwf(iexc)%SDe_up_old
					excwf(iexc)%ISDe_up_new=excwf(iexc)%ISDe_up_old
					excwf(iexc)%pvte_up_new=excwf(iexc)%pvte_up_old
					excwf(iexc)%SDe_dw_new=excwf(iexc)%SDe_dw_old
					excwf(iexc)%ISDe_dw_new=excwf(iexc)%ISDe_dw_old
					excwf(iexc)%pvte_dw_new=excwf(iexc)%pvte_dw_old
				END IF
				excwf(iexc)%detSDe_up_new=excwf(iexc)%detSDe_up_old
				excwf(iexc)%detSDe_dw_new=excwf(iexc)%detSDe_dw_old
			END IF
		ELSE IF ((num>0) .AND. (num<=N_part)) THEN
			SELECT CASE (tipo)
			CASE ('e__')
				IF (SDe_kind/='no_') THEN
					IF (num<=H_N_part) THEN
						excwf(iexc)%SDe_up_new(num,1:H_N_part)=excwf(iexc)%SDe_up_old(num,1:H_N_part)
						excwf(iexc)%detSDe_up_new=excwf(iexc)%detSDe_up_old
					ELSE
						excwf(iexc)%SDe_dw_new(num-H_N_part,1:H_N_part)=excwf(iexc)%SDe_dw_old(num-H_N_part,1:H_N_part)
						excwf(iexc)%detSDe_dw_new=excwf(iexc)%detSDe_dw_old
					END IF
				END IF
			END SELECT
		END IF
		
	END SUBROUTINE roll_back_exc
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_funzione_onda_eccitata(iexc,tipo,num)
		USE generic_tools
		IMPLICIT NONE
		CHARACTER(LEN=3) :: tipo
		INTEGER, INTENT(IN) :: iexc, num
		
		IF ((tipo=='all') .OR. (tipo=='e__')) THEN
			IF (SDe_kind/='no_') THEN
				IF (num==-1) THEN
					excwf(iexc)%SDe_up_old=excwf(iexc)%SDe_up_new
					excwf(iexc)%ISDe_up_old=excwf(iexc)%ISDe_up_new
					excwf(iexc)%pvte_up_old=excwf(iexc)%pvte_up_new
					excwf(iexc)%SDe_dw_old=excwf(iexc)%SDe_dw_new
					excwf(iexc)%ISDe_dw_old=excwf(iexc)%ISDe_dw_new
					excwf(iexc)%pvte_dw_old=excwf(iexc)%pvte_dw_new
				ELSE IF ((num>0) .AND. (num<=N_part)) THEN
					IF (num<=H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num,excwf(iexc)%ISDe_up_old, &
						  excwf(iexc)%detSDe_up_old,excwf(iexc)%SDe_up_new,excwf(iexc)%detSDe_up_new, &
						  excwf(iexc)%ISDe_up_new)
						excwf(iexc)%SDe_up_old(num,1:H_N_part)=excwf(iexc)%SDe_up_new(num,1:H_N_part)
						excwf(iexc)%ISDe_up_old=excwf(iexc)%ISDe_up_new
					ELSE IF (num>H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num-H_N_part,excwf(iexc)%ISDe_dw_old, &
						  excwf(iexc)%detSDe_dw_old,excwf(iexc)%SDe_dw_new,excwf(iexc)%detSDe_dw_new, &
						  excwf(iexc)%ISDe_dw_new)
						excwf(iexc)%SDe_dw_old(num-H_N_part,1:H_N_part)=excwf(iexc)%SDe_dw_new(num-H_N_part,1:H_N_part)
						excwf(iexc)%ISDe_dw_old=excwf(iexc)%ISDe_dw_new
					END IF
				END IF
			END IF
			excwf(iexc)%detSDe_up_old=excwf(iexc)%detSDe_up_new
			excwf(iexc)%detSDe_dw_old=excwf(iexc)%detSDe_dw_new
		END IF
		
	END SUBROUTINE aggiorna_funzione_onda_eccitata
!-----------------------------------------------------------------------
	SUBROUTINE valuta_excSD_har(iexc,num,r,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, iexc
		REAL (KIND=8), INTENT(IN) :: r(1:3,1:N)
		COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
		INTEGER :: i, j, ik, info, perm
		INTEGER, INTENT(OUT) :: pvt(1:N)
		COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato[ module_funzione_onda.f90 > valuta_SD_har ]'
		
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			 DO j = 1, N, 1
			 	DO i = 1, N, 1
			 		SD(i,j)=0.d0
			 		DO ik = 1, excwf(iexc)%num_pw_dnfH(j), 1
			 			SD(i,j)=SD(i,j)+excwf(iexc)%fattori_pw_dnfH(ik,j)*CDEXP((0.d0,1.d0)* &
			 			  DOT_PRODUCT(excwf(iexc)%k_pw_dnfH(1:3,ik,j),r(1:3,i)))
			 		END DO
			 		ISD(i,j)=SD(i,j)
			 	END DO
			 END DO
			!Calcolo il determinante di SD_new
			CALL ZGETRF( N, N, ISD, N, pvt, info )
			IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
			perm=0
			detSD=(1.d0,0.d0)
			DO  i = 1, N, 1
				IF (pvt(i) /= i) perm=perm+1
			END DO
			IF (MOD(perm,2) == 1 ) detSD=-detSD
			DO  i = 1, N, 1
				detSD=detSD*ISD(i,i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			DO j = 1, N, 1
				SD(num,j)=0.d0
				DO ik = 1, excwf(iexc)%num_pw_dnfH(j), 1
					SD(num,j)=SD(num,j)+excwf(iexc)%fattori_pw_dnfH(ik,j)*CDEXP((0.d0,1.d0)* &
					  DOT_PRODUCT(excwf(iexc)%k_pw_dnfH(1:3,ik,j),r(1:3,num)))
				END DO
			END DO
			CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
		ELSE
			STOP 'num non accettabile [ module_funzione_onda.f90 > valuta_SD_pw ]'
		END IF
		
	END SUBROUTINE valuta_excSD_har
!-----------------------------------------------------------------------
	SUBROUTINE chiudi_stati_eccitati()
		IMPLICIT NONE
		
		IF (.NOT. iniz_stati_eccitati) STOP 'prima di chiudere devi aver aperto [ module_stati_eccitati.f90 > chiudi_stati_eccitati ]'
		
		DEALLOCATE(i_exc_stat,energie_exc_stat)
		
		iniz_stati_Eccitati=.FALSE.
		
	END SUBROUTINE chiudi_stati_eccitati

END MODULE stati_eccitati

