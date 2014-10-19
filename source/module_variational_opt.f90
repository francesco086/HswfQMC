MODULE variational_opt
	USE dati_mc
	IMPLICIT NONE
	LOGICAL, PARAMETER :: verbose_mode=.TRUE.
	INTEGER, SAVE :: num_par_var, num_coord_Rp
	REAL (KIND=8) :: E_min(1:2), liv_precisione
	REAL (KIND=8), ALLOCATABLE, SAVE :: parametri_var(:)
	REAL (KIND=8), ALLOCATABLE :: s_kl(:,:), f_k(:)
	
	CONTAINS

	SUBROUTINE inizializza_variational_opt()
		USE dati_fisici
		USE funzione_onda
		USE walkers
		IMPLICIT NONE
		INTEGER :: cont, i, j, ih
		
		cont=1
		num_par_var=0
		num_coord_Rp=0
		
		CALL inizializza_dati_fisici()
		CALL inizializza_dati_mc()
		CALL inizializza_walkers('iniz_opt')
		CALL inizializza_dati_funzione_onda()
		
		IF (Jee_kind/='no_') THEN
			SELECT CASE (Jee_kind)
			CASE ('yuk')
				IF (opt_A_Jee) THEN
					IF (opt_F_Jee) THEN
						IF (split_Aee) num_par_var=num_par_var+1
						IF (split_Fee) num_par_var=num_par_var+1
						num_par_var=num_par_var+2
					ELSE
						IF (split_Aee) num_par_var=num_par_var+1
						num_par_var=num_par_var+1
					END IF
				ELSE
					IF (opt_F_Jee) THEN
						IF (split_Fee) num_par_var=num_par_var+1
						num_par_var=num_par_var+1
					END IF
				END IF
			CASE ('yup')
				IF (opt_A_Jee) THEN
					IF (opt_F_Jee) THEN
						IF (split_Aee) num_par_var=num_par_var+1
						IF (split_Fee) num_par_var=num_par_var+1
						num_par_var=num_par_var+2
					ELSE
						IF (split_Aee) num_par_var=num_par_var+1
						num_par_var=num_par_var+1
					END IF
				ELSE
					IF (opt_F_Jee) THEN
						IF (split_Fee) num_par_var=num_par_var+1
						num_par_var=num_par_var+1
					END IF
				END IF
			END SELECT
		END IF
		IF (Jep_kind/='no_') THEN
			SELECT CASE (Jep_kind)
			CASE ('yuk')
				IF (opt_A_Jep) THEN
					IF (opt_F_Jep) THEN
						IF (split_Aep) num_par_var=num_par_var+1
						IF (split_Fep) num_par_var=num_par_var+1
						num_par_var=num_par_var+2
					ELSE
						IF (split_Aep) num_par_var=num_par_var+1
						num_par_var=num_par_var+1
					END IF
				ELSE
					IF (opt_F_Jep) THEN
						IF (split_Fep) num_par_var=num_par_var+1
						num_par_var=num_par_var+1
					END IF
				END IF
			CASE ('yup')
				IF (opt_A_Jep) THEN
					IF (opt_F_Jep) THEN
						IF (split_Aep) num_par_var=num_par_var+1
						IF (split_Fep) num_par_var=num_par_var+1
						num_par_var=num_par_var+2
					ELSE
						IF (split_Aep) num_par_var=num_par_var+1
						num_par_var=num_par_var+1
					END IF
				ELSE
					IF (opt_F_Jep) THEN
						IF (split_Fep) num_par_var=num_par_var+1
						num_par_var=num_par_var+1
					END IF
				END IF
			CASE ('atm')
				num_par_var=num_par_var+1
			CASE ('atp')
				num_par_var=num_par_var+1
			END SELECT
		END IF
		IF (opt_Kse) THEN
			IF (Kse_kind/='no_') THEN
				SELECT CASE (Kse_kind)
				CASE ('gss')
					num_par_var=num_par_var+1
				CASE ('gsc')
					num_par_var=num_par_var+1
				CASE ('gsp')
					num_par_var=num_par_var+1
				CASE ('gsd')
					num_par_var=num_par_var+1
				CASE ('gdc')
					num_par_var=num_par_var+1
				CASE ('gdp')
					num_par_var=num_par_var+1
				CASE ('atm')
					num_par_var=num_par_var+1
				END SELECT
			END IF
		END IF
		IF (opt_Jse) THEN
			IF (Jse_kind/='no_') THEN
				SELECT CASE (Jse_kind)
				CASE ('pot')
					num_par_var=num_par_var+2
				CASE ('bou')
					num_par_var=num_par_var+2
				CASE ('ppb')
					num_par_var=num_par_var+3
				CASE ('yuk')
					IF (split_Asese) num_par_var=num_par_var+1
					IF (split_Fsese) num_par_var=num_par_var+1
					num_par_var=num_par_var+2
				CASE ('yup')
					IF (split_Asese) num_par_var=num_par_var+1
					IF (split_Fsese) num_par_var=num_par_var+1
					num_par_var=num_par_var+2
				END SELECT
			END IF
		END IF
		IF (opt_Jsesp) THEN
			IF (Jsesp_kind/='no_') THEN
				SELECT CASE (Jsesp_kind)
				CASE ('pot')
					num_par_var=num_par_var+1
				CASE ('yuk')
					IF (split_Asesp) num_par_var=num_par_var+1
					IF (split_Fsesp) num_par_var=num_par_var+1
					num_par_var=num_par_var+2
				CASE ('yup')
					IF (split_Asesp) num_par_var=num_par_var+1
					IF (split_Fsesp) num_par_var=num_par_var+1
					num_par_var=num_par_var+2
				CASE ('gss')
					num_par_var=num_par_var+1
				CASE ('gsd')
					num_par_var=num_par_var+1
				END SELECT
			END IF
		END IF
		IF (SDe_kind/='no_') THEN
			IF ( opt_SDe ) THEN
				SELECT CASE (SDe_kind)
				CASE ('prf')
					IF (opt_c_eff_dnfH) THEN
						ih=0
						DO j = 1, N_pw_lda, 1   !colonna
							DO i = 1, j-1, 1   !riga
								ih=ih+2
							END DO
						END DO
						DO j = 1, N_pw_lda, 1   !colonna
							ih=ih+1
						END DO
						num_par_var=num_par_var+ih
					END IF
				CASE ('fre')
					IF (opt_c_eff_dnfH) THEN
						ih=0
						DO j = 1, N_pw_lda, 1   !colonna
							DO i = 1, j-1, 1   !riga
								ih=ih+2
							END DO
						END DO
						DO j = 1, N_pw_lda, 1   !colonna
							ih=ih+1
						END DO
						num_par_var=num_par_var+ih
					END IF
				CASE ('atm')
					num_par_var=num_par_var+1
				CASE ('bat')
					num_par_var=num_par_var+1
				CASE ('atp')
					num_par_var=num_par_var+1
				END SELECT
			END IF
		END IF
		IF ( SDse_kind/='no_' ) THEN
			IF ( opt_SDse ) THEN
				SELECT CASE (SDse_kind)
				CASE ('atm')
					num_par_var=num_par_var+1
				CASE ('atp')
					num_par_var=num_par_var+1
				END SELECT
			END IF
		END IF
		
		IF ( opt_rp ) THEN
			num_par_var=num_par_var+3*N_part
			num_coord_Rp=3*N_part
		END IF
		
		!!!!!!!
		ALLOCATE(parametri_var(1:num_par_var))
		!!!!!!!
		
		IF (Jee_kind/='no_') THEN
			SELECT CASE (Jee_kind)
			CASE ('yuk')
				IF (opt_A_Jee) THEN
					IF (opt_F_Jee) THEN
						parametri_var(cont)=Aee_yuk
						cont=cont+1
						IF (split_Aee) THEN
							parametri_var(cont)=Aee_ud_yuk
							cont=cont+1
						END IF
						parametri_var(cont)=Fee_yuk
						cont=cont+1
						IF (split_Fee) THEN
							parametri_var(cont)=Fee_ud_yuk
							cont=cont+1
						END IF
					ELSE
						parametri_var(cont)=Aee_yuk
						cont=cont+1
						IF (split_Aee) THEN
							parametri_var(cont)=Aee_ud_yuk
							cont=cont+1
						END IF
					END IF
				ELSE
					IF (opt_F_Jee) THEN
						parametri_var(cont)=Fee_yuk
						cont=cont+1
						IF (split_Fee) THEN
							parametri_var(cont)=Fee_ud_yuk
							cont=cont+1
						END IF
					END IF
				END IF
			CASE ('yup')
				IF (opt_A_Jee) THEN
					IF (opt_F_Jee) THEN
						parametri_var(cont)=Aee_yuk
						cont=cont+1
						IF (split_Aee) THEN
							parametri_var(cont)=Aee_ud_yuk
							cont=cont+1
						END IF
						parametri_var(cont)=Fee_yuk
						cont=cont+1
						IF (split_Fee) THEN
							parametri_var(cont)=Fee_ud_yuk
							cont=cont+1
						END IF
					ELSE
						parametri_var(cont)=Aee_yuk
						cont=cont+1
						IF (split_Aee) THEN
							parametri_var(cont)=Aee_ud_yuk
							cont=cont+1
						END IF
					END IF
				ELSE
					IF (opt_F_Jee) THEN
						parametri_var(cont)=Fee_yuk
						cont=cont+1
						IF (split_Fee) THEN
							parametri_var(cont)=Fee_ud_yuk
							cont=cont+1
						END IF
					END IF
				END IF
			END SELECT
		END IF
		IF (Jep_kind/='no_') THEN
			SELECT CASE (Jep_kind)
			CASE ('yuk')
				IF (opt_A_Jep) THEN
					IF (opt_F_Jep) THEN
						parametri_var(cont)=Aep_yuk
						cont=cont+1
						IF (split_Aep) THEN
							parametri_var(cont)=Aep_ud_yuk
							cont=cont+1
						END IF
						parametri_var(cont)=Fep_yuk
						cont=cont+1
						IF (split_Fep) THEN
							parametri_var(cont)=Fep_ud_yuk
							cont=cont+1
						END IF
					ELSE
						parametri_var(cont)=Aep_yuk
						cont=cont+1
						IF (split_Aep) THEN
							parametri_var(cont)=Aep_ud_yuk
							cont=cont+1
						END IF
					END IF
				ELSE
					IF (opt_F_Jep) THEN
						parametri_var(cont)=Fep_yuk
						cont=cont+1
						IF (split_Fep) THEN
							parametri_var(cont)=Fep_ud_yuk
							cont=cont+1
						END IF
					END IF
				END IF
			CASE ('yup')
				IF (opt_A_Jep) THEN
					IF (opt_F_Jep) THEN
						parametri_var(cont)=Aep_yuk
						cont=cont+1
						IF (split_Aep) THEN
							parametri_var(cont)=Aep_ud_yuk
							cont=cont+1
						END IF
						parametri_var(cont)=Fep_yuk
						cont=cont+1
						IF (split_Fep) THEN
							parametri_var(cont)=Fep_ud_yuk
							cont=cont+1
						END IF
					ELSE
						parametri_var(cont)=Aep_yuk
						cont=cont+1
						IF (split_Aep) THEN
							parametri_var(cont)=Aep_ud_yuk
							cont=cont+1
						END IF
					END IF
				ELSE
					IF (opt_F_Jep) THEN
						parametri_var(cont)=Fep_yuk
						cont=cont+1
						IF (split_Fep) THEN
							parametri_var(cont)=Fep_ud_yuk
							cont=cont+1
						END IF
					END IF
				END IF
			CASE ('atm')
				parametri_var(cont)=Fep_yuk
				cont=cont+1
			CASE ('atp')
				parametri_var(cont)=Fep_yuk
				cont=cont+1
			END SELECT
		END IF
		IF (opt_Kse) THEN
			IF (Kse_kind/='no_') THEN
				SELECT CASE (Kse_kind)
				CASE ('gss')
					parametri_var(cont)=C_kern_e
					cont=cont+1
				CASE ('gsc')
					parametri_var(cont)=C_kern_e
					cont=cont+1
				CASE ('gsp')
					parametri_var(cont)=C_kern_e
					cont=cont+1
				CASE ('gsd')
					parametri_var(cont)=C_kern_e
					cont=cont+1
				CASE ('gdc')
					parametri_var(cont)=C_kern_e
					cont=cont+1
				CASE ('gdp')
					parametri_var(cont)=C_kern_e
					cont=cont+1
				CASE ('atm')
					parametri_var(cont)=C_kern_e
					cont=cont+1
				END SELECT
			END IF
		END IF
		IF (opt_Jse) THEN
			IF (Jse_kind/='no_') THEN
				SELECT CASE (Jse_kind)
				CASE ('pot')
					parametri_var(cont)=A_POT_se
					parametri_var(cont+1)=D_POT_se
					cont=cont+2
				CASE ('bou')
					parametri_var(cont)=B_se
					parametri_var(cont+1)=D_se
					cont=cont+2
				CASE ('ppb')
					parametri_var(cont)=c_se
					parametri_var(cont+1)=B_se
					parametri_var(cont+2)=D_se
					cont=cont+3
				CASE ('yuk')
					parametri_var(cont)=Asese_yuk
					cont=cont+1
					IF (split_Asese) THEN
						parametri_var(cont)=Asese_ud_yuk
						cont=cont+1
					END IF
					parametri_var(cont)=Fsese_yuk
					cont=cont+1
					IF (split_Fsese) THEN
						parametri_var(cont)=Fsese_ud_yuk
						cont=cont+1
					END IF
				CASE ('yup')
					parametri_var(cont)=Asese_yuk
					cont=cont+1
					IF (split_Asese) THEN
						parametri_var(cont)=Asese_ud_yuk
						cont=cont+1
					END IF
					parametri_var(cont)=Fsese_yuk
					cont=cont+1
					IF (split_Fsese) THEN
						parametri_var(cont)=Fsese_ud_yuk
						cont=cont+1
					END IF
				END SELECT
			END IF
		END IF
		IF (opt_Jsesp) THEN
			IF (Jsesp_kind/='no_') THEN
				SELECT CASE (Jsesp_kind)
				CASE ('pot')
					parametri_var(cont)=c_sesp
					cont=cont+1
				CASE ('yuk')
					parametri_var(cont)=Asesp_yuk
					cont=cont+1
					IF (split_Asesp) THEN
						parametri_var(cont)=Asesp_ud_yuk
						cont=cont+1
					END IF
					parametri_var(cont)=Fsesp_yuk
					cont=cont+1
					IF (split_Fsesp) THEN
						parametri_var(cont)=Fsesp_ud_yuk
						cont=cont+1
					END IF
				CASE ('yup')
					parametri_var(cont)=Asesp_yuk
					cont=cont+1
					IF (split_Asesp) THEN
						parametri_var(cont)=Asesp_ud_yuk
						cont=cont+1
					END IF
					parametri_var(cont)=Fsesp_yuk
					cont=cont+1
					IF (split_Fsesp) THEN
						parametri_var(cont)=Fsesp_ud_yuk
						cont=cont+1
					END IF
				CASE ('gss')
					parametri_var(cont)=Gsesp
					cont=cont+1
				CASE ('gsd')
					parametri_var(cont)=Gsesp
					cont=cont+1
				END SELECT
			END IF
		END IF
		IF (SDe_kind/='no_') THEN
			IF ( opt_SDe ) THEN
				SELECT CASE (SDe_kind)
				CASE ('prf')
					IF (opt_c_eff_dnfH) THEN
						ih=0
						DO j = 1, N_pw_lda, 1   !colonna
							DO i = 1, j-1, 1   !riga
								ih=ih+2
								parametri_var(ih-1)=REAL(c_eff_dnfH(i,j),8)
								parametri_var(ih)=DIMAG(c_eff_dnfH(i,j))
							END DO
						END DO
						DO j = 1, N_pw_lda, 1   !colonna
							i=j
							ih=ih+1
							parametri_var(ih)=REAL(c_eff_dnfH(i,j),8)
						END DO
					END IF
				CASE ('fre')
					IF (opt_c_eff_dnfH) THEN
						ih=0
						DO j = 1, N_pw_lda, 1   !colonna
							DO i = 1, j-1, 1   !riga
								ih=ih+2
								parametri_var(ih-1)=REAL(c_eff_dnfH(i,j),8)
								parametri_var(ih)=DIMAG(c_eff_dnfH(i,j))
							END DO
						END DO
						DO j = 1, N_pw_lda, 1   !colonna
							i=j
							ih=ih+1
							parametri_var(ih)=REAL(c_eff_dnfH(i,j),8)
						END DO
					END IF
				CASE ('atm')
					parametri_var(cont)=C_atm
					cont=cont+1
				CASE ('bat')
					parametri_var(cont)=C_atm
					cont=cont+1
				CASE ('atp')
					parametri_var(cont)=C_atm
					cont=cont+1
				END SELECT
			END IF
		END IF
		IF ( SDse_kind/='no_' ) THEN
			IF ( opt_SDse ) THEN
				SELECT CASE (SDse_kind)
				CASE ('atm')
					parametri_var(cont)=C_atm
					cont=cont+1
				CASE ('atp')
					parametri_var(cont)=C_atm
					cont=cont+1
				END SELECT
			END IF
		END IF
				
		IF ( opt_rp ) THEN
			CALL inizializza_dati_mc()
			CALL inizializza_walkers('opt_Rp')
			DO i = 1, N_part, 1
				parametri_var(cont:cont+2)=rp_old(1:3,i)
				cont=cont+3
			END DO
			CALL chiudi_walkers()
			CALL chiudi_dati_mc()
		END IF
						
		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: il numero di parametri da ottimizzare é ', num_par_var
			IF (flag_output) WRITE (7, *), 'VAR_OPT: il numero di parametri da ottimizzare é ', num_par_var
		END IF
		
		CALL chiudi_dati_funzione_onda()
		CALL chiudi_walkers()
		CALL chiudi_dati_mc()
		CALL chiudi_dati_fisici()
		
	END SUBROUTINE inizializza_variational_opt
	
	SUBROUTINE minimizza_energia(target_precisione,metodo)
		USE dati_fisici
		USE dati_mc
		USE funzione_onda
		USE walkers
		IMPLICIT NONE
		CHARACTER(LEN=4) :: stringa
		CHARACTER(LEN=7) :: metodo
		REAL (KIND=8), INTENT(IN) :: target_precisione
		
		liv_precisione=target_precisione
		CALL inizializza_variational_opt()
		SELECT CASE(metodo)
		CASE('congrad')
			CALL conj_grad(parametri_var,num_par_var,E_min)
		CASE('axisopt')
			CALL orthogonal_basis_search(parametri_var,num_par_var,E_min)
		CASE('stocrec')
			CALL stochastic_reconfiguration(parametri_var,num_par_var,.TRUE.,.FALSE.)
		CASE('pure_sr')
			CALL pure_stochastic_reconfiguration(parametri_var,num_par_var,.TRUE.)
		CASE('stoc_ns')   !SR senza fermarsi
			CALL stochastic_reconfiguration(parametri_var,num_par_var,.FALSE.,.FALSE.)
		CASE('stoc_av')   !SR senza fermarsi
			CALL stochastic_reconfiguration(parametri_var,num_par_var,.FALSE.,.TRUE.)
		END SELECT
		IF (mpi_myrank==0) THEN
			WRITE (stringa, '(I4.4)'), num_par_var
			PRINT * , '    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # '
			PRINT * , 'VAR_OPT: Ottimizzazione terminata.'
			PRINT * , 'Optimal parameters: '
			PRINT '('//stringa//'(F9.3,1X))' , parametri_var
			PRINT '(11x,A9,F9.6,A6,F9.6)' , 'Energia: ', E_min(1), '  +-  ', E_min(2)
			PRINT * , '    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # '
			IF (flag_output) THEN
				WRITE (7, *) , '    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # '
				WRITE (7, *) , 'VAR_OPT: Ottimizzazione terminata.'
				WRITE (7, *) , 'Optimal parameters: '
				WRITE (7, '('//stringa//'(F9.3,1X))') , parametri_var
				WRITE (7, '(11X,A9,F9.6,A6,F9.6)') , 'Energia: ', E_min(1), '  +-  ', E_min(2)
				WRITE (7, *) , '    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # '
			END IF
		END IF
		IF (stampa_dati_funzione_onda) THEN
			CALL inizializza_dati_fisici()
			CALL inizializza_dati_mc()
			CALL inizializza_walkers('opt_Rp')
			CALL inizializza_dati_funzione_onda()
			CALL setta_parametri(parametri_var,num_par_var)
			CALL stampa_file_dati_funzione_onda('ottimizzazione/OPT_wf.d')
			CALL chiudi_dati_funzione_onda()
			CALL chiudi_walkers()
			CALL chiudi_dati_mc()
			CALL chiudi_dati_fisici()
		END IF
		CALL chiudi_variational_opt()
	END SUBROUTINE minimizza_energia

	SUBROUTINE calcola_energia(id,par,num,energia,accettabile)
		USE VMC
		IMPLICIT NONE
		LOGICAL, INTENT(OUT) :: accettabile
		CHARACTER(LEN=2) :: id
		INTEGER :: num
		REAL (KIND=8) :: par(1:num), energia(1:2)
		
		cont_N_mc_increase=0
				
		CALL inizializza_VMC('opt'//id,CONTINUA=.FALSE.)
		CALL setta_parametri(par,num)
		CALL prima_valutazione_funzione_onda()
		CALL valuta_step_mc(accettabile)
		IF (accettabile) THEN
			CALL sampling(.FALSE.)
			CALL sampling(.TRUE.)
			CALL trascrivi_dati()
			CALL restituisci_energia(energia)
			CALL chiudi_VMC()
		ELSE
			CALL chiudi_VMC()
		END IF
		
	END SUBROUTINE calcola_energia

	SUBROUTINE aumenta_precisione_energia(id,par,num,energia)
		USE VMC
		IMPLICIT NONE
		LOGICAL :: accettabile
		CHARACTER(LEN=2) :: id
		INTEGER :: num
		REAL (KIND=8) :: par(1:num), energia(1:2)
		
		IF (flag_disk) THEN
			IF ((mpi_myrank==0) .AND. (verbose_mode)) THEN
				PRINT * , 'VAR_OPT: aumento l accuratezza della simulazione ', id
				IF (flag_output) WRITE (7, *) 'VAR_OPT: aumento l accuratezza della simulazione ', id
			END IF
			CALL inizializza_VMC('opt'//id,CONTINUA=.TRUE.)
			CALL setta_parametri(par,num)
			CALL prima_valutazione_funzione_onda()
			CALL valuta_step_mc(accettabile)
			CALL sampling(.FALSE.)
			CALL sampling(.TRUE.)
			CALL trascrivi_dati()
			CALL restituisci_energia(energia)
			CALL chiudi_VMC()
		ELSE
			cont_N_mc_increase=cont_N_mc_increase+1
			IF (cont_N_mc_increase>5) STOP 'N_mc é diventato troppo grande &
			  [ module_variational_opt.f90 > aumenta_precisione_energie ]'
			IF ((mpi_myrank==0) .AND. (verbose_mode)) THEN
				PRINT * , 'VAR_OPT: aumento l accuratezza della simulazione ', id
				IF (flag_output) WRITE (7, *), 'VAR_OPT: aumento l accuratezza della simulazione ', id
			END IF
			CALL inizializza_VMC('opt'//id,CONTINUA=.TRUE.)
			CALL setta_parametri(par,num)
			CALL prima_valutazione_funzione_onda()
			CALL valuta_step_mc(accettabile)
			CALL sampling(.FALSE.)
			CALL sampling(.TRUE.)
			CALL trascrivi_dati()
			CALL restituisci_energia(energia)
			CALL chiudi_VMC()
		END IF

	END SUBROUTINE aumenta_precisione_energia

	SUBROUTINE calcola_energie(id,par_pt,num_pt,num,energie,accettabile)
		USE VMC
		USE variational_calculations
		IMPLICIT NONE
		LOGICAL, INTENT(OUT) :: accettabile
		CHARACTER(LEN=2) :: id
		INTEGER :: num_pt, num                  !num_pt=numero punti da considerare per il gradiente (generalmente = num)
		REAL (KIND=8) :: par_pt(1:num,0:num_pt)
		REAL (KIND=8), INTENT(OUT) :: energie(1:2,0:num_pt)
		
		cont_N_mc_increase=0
		
		CALL inizializza_VMC('opt'//id,CONTINUA=.FALSE.)
		CALL setta_parametri(par_pt(1:num,0),num)
		CALL prima_valutazione_funzione_onda()
		CALL valuta_step_mc(accettabile)
		IF (accettabile) THEN
			CALL sampling(.FALSE.)
			CALL inizializza_variational_calculations(num,par_pt,num_pt)
			CALL sampling(.TRUE.)                 !devo aggiungere qualcosa qui per far calcolare anche i punti per il gradiente
			CALL trascrivi_dati()
			CALL restituisci_energie(energie)     !devo modificarlo
			CALL chiudi_VMC()
			CALL chiudi_variational_calculations()
		ELSE
			CALL chiudi_VMC()
		END IF
	END SUBROUTINE calcola_energie

	SUBROUTINE aumenta_precisione_energie(id,par_pt,num_pt,num,energie,accettabile)
		USE VMC
		USE variational_calculations
		IMPLICIT NONE
		LOGICAL, INTENT(OUT) :: accettabile
		CHARACTER(LEN=2) :: id
		INTEGER :: num_pt, num                  !num_pt=numero punti da considerare per il gradiente (generalmente = num)
		REAL (KIND=8) :: par_pt(1:num,0:num_pt)
		REAL (KIND=8), INTENT(OUT) :: energie(1:2,0:num_pt)
		
		IF (flag_disk) THEN
			CALL inizializza_VMC('opt'//id,CONTINUA=.TRUE.)
			CALL setta_parametri(par_pt(1:num,0),num)
			CALL prima_valutazione_funzione_onda()
			CALL valuta_step_mc(accettabile)
			CALL sampling(.FALSE.)
			CALL inizializza_variational_calculations(num,par_pt,num_pt)
			CALL sampling(.TRUE.)                 !devo aggiungere qualcosa qui per far calcolare anche i punti per il gradiente
			CALL trascrivi_dati()
			CALL restituisci_energie(energie)     !devo modificarlo
			CALL chiudi_VMC()
			CALL chiudi_variational_calculations()
		ELSE
			cont_N_mc_increase=cont_N_mc_increase+1
			IF (cont_N_mc_increase>5) STOP 'N_mc é diventato troppo grande &
			  [ module_variational_opt.f90 > aumenta_precisione_energie ]'
			CALL inizializza_VMC('opt'//id,CONTINUA=.TRUE.)
			CALL setta_parametri(par_pt(1:num,0),num)
			CALL prima_valutazione_funzione_onda()
			CALL valuta_step_mc(accettabile)
			CALL sampling(.FALSE.)
			CALL inizializza_variational_calculations(num,par_pt,num_pt)
			CALL sampling(.TRUE.)                 !devo aggiungere qualcosa qui per far calcolare anche i punti per il gradiente
			CALL trascrivi_dati()
			CALL restituisci_energie(energie)     !devo modificarlo
			CALL chiudi_VMC()
			CALL chiudi_variational_calculations()
		END IF
	END SUBROUTINE aumenta_precisione_energie

	SUBROUTINE SR_campiona_wf_attuale(id,par,num,energia,accettabile)
		USE VMC
		USE variational_calculations
		IMPLICIT NONE
		LOGICAL, INTENT(OUT) :: accettabile
		CHARACTER(LEN=2) :: id
		INTEGER :: num, i
		REAL (KIND=8) :: par(1:num), energia(1:2), par_pt(1:num,1)
		
		CALL inizializza_VMC('opt'//id,CONTINUA=.FALSE.)

		CALL setta_parametri(par,num)
		CALL prima_valutazione_funzione_onda()
		CALL valuta_step_mc(accettabile)
		IF (accettabile) THEN
			CALL sampling(.FALSE.)
			CALL inizializza_variational_calculations(num,par_pt,1)
			CALL sampling(.TRUE.)
			CALL trascrivi_dati()
			CALL restituisci_energia(energia)
			IF ((flag_disk .AND. mpi_myrank==0) .OR. (.NOT. flag_disk)) THEN
				CALL estrai_s_kl_e_f_k(energia)
			END IF
			CALL chiudi_VMC()
			CALL chiudi_variational_calculations()
		ELSE
			CALL chiudi_VMC()
		END IF
	END SUBROUTINE SR_campiona_wf_attuale
!-----------------------------------------------------------------------

	SUBROUTINE calcola_g_Rp(Rp,guu_Rp,gud_Rp)
		USE dati_fisici
		USE generic_tools
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: Rp(1:3*N_part)
		INTEGER :: i1, i2
		REAL (KIND=8) :: r_lattice(1:3,1:N_part), rij_lattice(0:3,1:N_part,1:N_part)
		REAL (KIND=8), INTENT(OUT) :: guu_Rp(1:N_hist), gud_Rp(1:N_hist)
		
		!trasporto l'array Rp in coordinate 3D nel vettore r_lattice
		DO i1 = 1, N_part, 1
			r_lattice(1:3,i1)=Rp(1+(i1-1)*3:i1*3)
		END DO
		
		!calcolo le distanze fra le coordinate r_lattice (considenrando le boundary conditions)
		CALL valuta_distanza_ii(r_lattice,N_part,L,rij_lattice)
		
		!trovo la g(r) corrispondente
		CALL g_ii_correlation(rij_lattice,N_part,L,r_s,N_hist,guu_Rp,gud_Rp)
		
	END SUBROUTINE calcola_g_Rp
!-----------------------------------------------------------------------

	SUBROUTINE orthogonal_basis_search(p0,N,minimo)
		USE generic_tools
		IMPLICIT NONE
		INTEGER, PARAMETER :: N_max_rot=5
		INTEGER (KIND=4), INTENT(IN) :: N
		LOGICAL :: flag_loop
		INTEGER :: i, contatore
		REAL (KIND=8) :: p0(1:N), minimo(1:2), minimo_old(1:2)
		REAL (KIND=8) :: basis(1:N,1:N), teta(1:N,1:N)
		
		basis=0.d0
		DO i = 1, N, 1
			basis(i,i)=1.d0
		END DO
		
		CALL genera_matrice_rotazione_eulero(N,teta)
		DO i = 1, N, 1
			IF (mpi_myrank==0) CALL ruota_vettore_Nd_eulero(basis(1:N,i),N,teta)
			CALL MPI_BCAST(basis(1:N,i),N,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
		END DO
		CALL axis_search(basis,p0,N,minimo)
		minimo_old=1000.d0
		flag_loop=.TRUE.
		contatore=0
		DO WHILE (flag_loop)
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: ruoto la base per la ricerca del minimo'
				IF (flag_output) WRITE (7, *), 'VAR_OPT: ruoto la base per la ricerca del minimo'
			END IF
			CALL genera_matrice_rotazione_eulero(N,teta)
			DO i = 1, N, 1
				IF (mpi_myrank==0) CALL ruota_vettore_Nd_eulero(basis(1:N,i),N,teta)
				CALL MPI_BCAST(basis(1:N,i),N,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
			END DO
			minimo_old=minimo
			CALL axis_search(basis,p0,N,minimo)
			IF (DABS(minimo_old(1)-minimo(1))>0.5d0*liv_precisione) THEN
				contatore=0
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: Ruotando la base ho trovato un nuovo minimo'
					IF (flag_output) WRITE (7, *), 'VAR_OPT: Ruotando la base ho trovato un nuovo minimo'
				END IF
			ELSE
				contatore=contatore+1
			END IF
			IF (contatore>N_max_rot) flag_loop=.FALSE.
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: Ruotare la base é stato inutile'
				IF (flag_output) WRITE (7, *), 'VAR_OPT: Ruotare la base é stato inutile'
			END IF
		END DO
		
	END SUBROUTINE orthogonal_basis_search

	SUBROUTINE axis_search(basis,p0,N,minimo)
		USE dati_fisici
		USE funzione_onda
		USE VMC
		IMPLICIT NONE
		INTEGER (KIND=4), INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: basis(1:N,1:N)
		LOGICAL :: accettabile
		CHARACTER(LEN=2) :: id, ida, idb, idc
		INTEGER (KIND=4) :: i
		REAL (KIND=8) :: a, fa(1:2), b, fb(1:2), c, fc(1:2), fold(1:2), fnew(1:2)
		REAL (KIND=8) :: p0(1:N), minimo(1:2), p_old(1:N)
		REAL (KIND=8) :: first_step, direzione(1:N), cambiamento
		
		!cambiamento=1000.d0
		!DO WHILE (cambiamento>liv_precisione*0.25d0)
			DO i = 1, N, 1
				IF (i==1) THEN
					ida='00'
					CALL controlla_punto(p0)
					CALL calcola_energia(ida,p0,N,fa(1:2),accettabile)
					fold=fa
				ELSE
				 	ida=idb
					fa=fb
				END IF
				a=0.d0
				direzione=basis(1:N,i)
				first_step=0.1d0*DSQRT(REAL(N,8))      !per la subroutine bracket
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: cerco il minimo lungo l asse ', i
					PRINT "(8X,20(2X,F5.2))" , direzione
					IF (flag_output) THEN
						WRITE (7, *), 'VAR_OPT: cerco il minimo lungo l asse ', i
						WRITE (7, "(8X,20(2X,F5.2))"), direzione
					END IF
				END IF
				CALL bracket_min_multidim_v2(N,p0,direzione,first_step,a,fa,ida,b,fb,idb,c,fc,idc)
				CALL parabgold_search_multidim(N,p0,direzione,a,fa,ida,b,fb,idb,c,fc,idc)
				!CALL SYSTEM ('export ID_TO_ERASE=opt00; find estimatori/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \; find posizioni/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \')
				p0=p0+b*direzione
				IF (stampa_dati_funzione_onda) THEN
					CALL inizializza_dati_fisici()
					CALL inizializza_dati_funzione_onda()
					CALL setta_parametri(p0,N)
					CALL stampa_file_dati_funzione_onda('ottimizzazione/wf_opt_axis.d')
					CALL chiudi_dati_fisici()
					CALL chiudi_dati_funzione_onda()
				END IF
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: Seguendo l asse ', i, ' ho trovato come nuovo minimo ', &
					  fb(1), ' +- ', fb(2)
					PRINT * , '         usando il fattore b=', b
					IF (flag_output) THEN
						WRITE (7, *), 'VAR_OPT: Seguendo l asse ', i, ' ho trovato come nuovo minimo ', &
						  fb(1), ' +- ', fb(2)
						WRITE (7, *), '         usando il fattore b=', b
					END IF
				END IF
			END DO
			fnew=fb
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: Alla fine di una routine su tutti gli assi ho trovato come nuovo minimo ', &
				  fnew(1), ' +- ', fnew(2)
				PRINT * , '                                                                     partendo da ', &
				  fold(1), ' +- ', fold(2)
				IF (flag_output) THEN
					WRITE (7, *) , 'VAR_OPT: Alla fine di una routine su tutti gli assi ho trovato come nuovo minimo ', &
					  fnew(1), ' +- ', fnew(2)
					WRITE (7, *) , '                                                                     partendo da ', &
					  fold(1), ' +- ', fold(2)
				END IF
			END IF
			cambiamento=fold(1)-fnew(1)
		!END DO
		
		minimo(1:2)=fb(1:2)
		
	END SUBROUTINE axis_search

	SUBROUTINE conj_grad(p0,N,minimo)
		USE dati_fisici
		USE funzione_onda
		IMPLICIT NONE
		LOGICAL :: gradiente_nullo, accettabile, migliora_gradiente, flag_file
		CHARACTER(LEN=2) :: id, ida, idb, idc
		REAL (KIND=8), PARAMETER :: step_min=0.01d0
		INTEGER (KIND=4), INTENT(IN) :: N
		REAL (KIND=8) :: p0(1:N)
		INTEGER (KIND=4) :: n1_char, n2_char
		INTEGER (KIND=8) :: i, contatore, cont_grad
		REAL (KIND=8) :: p(1:N), step, f0(1:2,0:N), grad(1:N), h(1:N), g(1:N), g_old(1:N), gamma
		REAL (KIND=8) :: pts(1:N,0:N)
		REAL (KIND=8) :: a, b, c, minimo(1:2)
		REAL (KIND=8) :: fa(1:2), fb(1:2), fc(1:2)
		REAL (KIND=8) :: direzione(1:N), dir(1:N), partenza(1:N), first_step
		
		contatore=1
		p=p0
		step=0.005d0          !per calcolare il gradiente
		id='00'
		
		!PRINT * , "p0=", p0
		
		CALL controlla_punto(p0)
		
		!PRINT * , "p*=", p0
		
		pts(1:N,0)=p0(1:N)
		DO i = 1, N, 1
			pts(1:N,i)=p0(1:N)
			pts(i,i)=pts(i,i)+step
		END DO
		id='00'
		CALL calcola_energie(id,pts(1:N,0:N),N,N,f0(1:2,0:N),accettabile)		
		IF (.NOT. accettabile) STOP 'I parametri variazionali iniziali non sono accettabili &
		  [ module_variational_opt.f90 > conj_grad ]'
		DO i = 1, N, 1
			grad(i)=f0(1,i)/step
		END DO
		
		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Ho calcolato il primo gradiente. La direzione da seguire é '
			IF (flag_output) WRITE (7, *), 'VAR_OPT: Ho calcolato il primo gradiente. La direzione da seguire é '
			DO i = 1, N, 1
				PRINT * , i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
				IF (flag_output) WRITE (7, *), i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
			END DO
		END IF
		
		migliora_gradiente=.FALSE.
		DO i = 1, N, 1
			IF (DABS(f0(1,i))<2.d0*f0(2,i)) migliora_gradiente=.TRUE.
		END DO
		cont_grad=1
		DO WHILE (migliora_gradiente)
			CALL aumenta_precisione_energie(id,pts(1:N,0:N),N,N,f0(1:2,0:N),accettabile)
			DO i = 1, N, 1
				grad(i)=f0(1,i)/step
			END DO
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: Ho migliorato la precisione del gradiente. La direzione da seguire é '
				IF (flag_output) WRITE (7, *), 'VAR_OPT: Ho migliorato la precisione del gradiente. La direzione da seguire é '
				DO i = 1, N, 1
					PRINT * , i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
					WRITE (7, *), i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
				END DO
			END IF
			migliora_gradiente=.FALSE.
			DO i = 1, N, 1
				IF (DABS(f0(1,i))<2.d0*f0(2,i)) migliora_gradiente=.TRUE.
			END DO
			cont_grad=cont_grad+1
			IF (cont_grad>2) migliora_gradiente=.FALSE. 
		END DO
		
		g=-grad
		h=g
		
		partenza=p0
		direzione=h
		
		cerca_minimo: IF (.NOT. gradiente_nullo) THEN
		
		a=0.d0
		fa=f0(1:2,0)
		first_step=0.1d0*N      !per la subroutine bracket
		CALL bracket_min_multidim(N,partenza,direzione,first_step,a,fa,ida,b,fb,idb,c,fc,idc)
		CALL parabgold_search_multidim(N,partenza,direzione,a,fa,ida,b,fb,idb,c,fc,idc)
		IF (mpi_myrank==0) THEN 
			INQUIRE(FILE='posizioni/inizio-opt00_e_0000.pos',EXIST=flag_file)
			IF (flag_file) CALL SYSTEM ("export ID_TO_ERASE=opt00; &
			  find estimatori/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \; &
			  find posizioni/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \")
		END IF
		p0=partenza+b*direzione

		DO WHILE ((b*MAXVAL(DABS(direzione))>step) .AND. (.NOT. gradiente_nullo))	
			CALL controlla_punto(p0)
			pts(1:N,0)=p0(1:N)
			DO i = 1, N, 1
				pts(1:N,i)=p0(1:N)
				pts(i,i)=pts(i,i)+step
			END DO
			id='00'
			IF (stampa_dati_funzione_onda) THEN
				CALL inizializza_dati_fisici()
				CALL inizializza_dati_funzione_onda()
				CALL setta_parametri(p0,N)
				CALL stampa_file_dati_funzione_onda('ottimizzazione/wf_conjg.d')
				CALL chiudi_dati_fisici()
				CALL chiudi_dati_funzione_onda()
			END IF
			CALL calcola_energie(id,pts(1:N,0:N),N,N,f0(1:2,0:N),accettabile)
			DO i = 1, N, 1
				grad(i)=f0(1,i)/step
			END DO
			
			migliora_gradiente=.FALSE.
			DO i = 1, N, 1
				IF (DABS(f0(1,i))<2.d0*f0(2,i)) migliora_gradiente=.TRUE.
			END DO
			cont_grad=1
			DO WHILE (migliora_gradiente)
				CALL aumenta_precisione_energie(id,pts(1:N,0:N),N,N,f0(1:2,0:N),accettabile)
				DO i = 1, N, 1
					grad(i)=f0(1,i)/step
				END DO
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: Ho migliorato la precisione del gradiente. La direzione da seguire é '
					IF (flag_output) WRITE (7, *), 'VAR_OPT: Ho migliorato la precisione del gradiente. La direzione da seguire é ' 
					DO i = 1, N, 1
						PRINT * , i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
						IF (flag_output) WRITE (7, *), i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
					END DO
				END IF
				migliora_gradiente=.FALSE.
				DO i = 1, N, 1
					IF (DABS(f0(1,i))<2.d0*f0(2,i)) migliora_gradiente=.TRUE.
				END DO
				cont_grad=cont_grad+1
				IF (cont_grad>5) migliora_gradiente=.FALSE.
			END DO
			
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: Ho calcolato il gradiente. La direzione da seguire é '
				IF (flag_output) WRITE (7, *), 'VAR_OPT: Ho calcolato il gradiente. La direzione da seguire é '
				DO i = 1, N, 1
					PRINT * , i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
					IF (flag_output) WRITE (7, *), i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
				END DO
			END IF
			g_old=g
			g=-grad
			gamma=DOT_PRODUCT((g-g_old),g)/DOT_PRODUCT(g,g)
			h=g+gamma*h
			partenza=p0
			direzione=h
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: ho calcolato la nuova direzione'
				PRINT * , direzione
				IF (flag_output) THEN
					WRITE (7, *) , 'VAR_OPT: ho calcolato la nuova direzione'
					WRITE (7, *) , direzione
				END IF
			END IF
			IF (.NOT. gradiente_nullo) THEN
				a=0.d0
				fa=f0(1:2,0)
				first_step=0.1d0*N      !per la subroutine bracket
				CALL bracket_min_multidim(N,partenza,direzione,first_step,a,fa,ida,b,fb,idb,c,fc,idc)
				CALL parabgold_search_multidim(N,partenza,direzione,a,fa,ida,b,fb,idb,c,fc,idc)
				p0=partenza+b*direzione
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: il nuovo punto per il minimo provvisorio é '
					PRINT * , p0
					IF (flag_output) THEN
						WRITE (7, *) , 'VAR_OPT: il nuovo punto per il minimo provvisorio é '
						WRITE (7, *) , p0
					END IF
				END IF
				contatore=contatore+1
			END IF
		END DO
		minimo=fb(1:2)
		
		END IF cerca_minimo
		
		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Conjugate Gradient - sono stati necessari ', contatore, ' passi'
			IF (flag_output) WRITE (7, *), 'VAR_OPT: Conjugate Gradient - sono stati necessari ', contatore, ' passi'
		END IF

	END SUBROUTINE conj_grad

	SUBROUTINE stochastic_reconfiguration(p0,N,auto_stop,AV_accu)
		USE dati_fisici
		USE funzione_onda
		USE walkers
		USE generic_tools
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: MIN_VAR_NEC=0.025, MIN_LAMBDA2=0.5, MIN_VAR_ACC=0.001, MAX_VAR_ACC=0.25, MAX_VAR_PROT=0.1
		LOGICAL, INTENT(IN) :: auto_stop, AV_accu
		INTEGER (KIND=4), INTENT(IN) :: N
		LOGICAL :: flag_loop, loop_stop(1:6), accettabile, flag_file, flag_trovato_minimo, flag_freeze, flag_AV_wf
		CHARACTER(LEN=2) :: id
		CHARACTER(LEN=4) :: stringa, istring
		INTEGER :: i, j, info, contatore, IO, i_orbit, AV_cont
		REAL (KIND=8) :: lambda, eps, lambda2, lambda3, boundd, dummy(1:4), dummy2(0:N)
		REAL (KIND=8) :: lambda_Rp, lambda2_Rp, flambda2_Rp
		REAL (KIND=8) :: p0(1:N), dp_iniz(1:N), AV_P(1:N)
		REAL (KIND=8) :: dp(1:N), Is_kl(1:N,1:N), pvt(1:N), cont_step
		REAL (KIND=8) :: work(1:3*N), D_tr(1:N), E_tr(1:N-1), TAU_tr(1:N-1)
		REAL (KIND=8) :: energia(1:2), vec_app(1:N), vec_app_old(1:N), flambda2
		REAL (KIND=8) :: old_step(1:5), old_p0(1:N,1:5), old_E(1:2,1:5), derivE(1:4)
		REAL (KIND=8) :: old_step_Rp(1:5)
		REAL (KIND=8) :: dist_poss, dist_perc, dist_poss_5, dist_perc_5, av_errore, derE, nextdE
		REAL (KIND=8) :: old_Emin(1:2), p0_minimo(1:N), Pdist_min(1:5)
		REAL (KIND=8), ALLOCATABLE :: gRpuu(:), gRpuu_new(:), gRpud(:), gRpud_new(:), x_gr(:)
		lambda2=1.d0
		lambda2_Rp=1.d0
		id='SR'
		i_orbit=1
		AV_cont=0
		AV_P=0.d0
		contatore=0
		
		INQUIRE(FILE='ottimizzazione/SR_energies.dat',EXIST=flag_file)
		IF ((flag_file).AND.(mpi_myrank==0)) THEN
			OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='OLD')
			IO=0
			DO WHILE (IO>=0)
				READ (8, *, IOSTAT = IO), dummy(1:3)
				IF (IO>=0) dummy(4)=dummy(1) 
			END DO
			cont_step=dummy(4)
			CLOSE(8)
			OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='OLD', POSITION='APPEND')
			IF (mpi_myrank==0) PRINT * , 'VAR_OPT: esisteva giá un file, parto da ', cont_step, '.'
		ELSE IF (mpi_myrank==0) THEN
			OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='NEW')
			cont_step=0.d0
			IF (mpi_myrank==0) PRINT * , 'VAR_OPT: parto da zero.'
		END IF
				
		OPEN (UNIT=9, FILE='ottimizzazione/SR_var_parameters.dat', STATUS='UNKNOWN', POSITION='APPEND')
		IF ( num_coord_Rp>0 ) THEN   !stampo la configurazione Rp_0 nel file 
			IF (stampa_dati_funzione_onda) THEN
				CALL inizializza_dati_fisici()
				CALL inizializza_dati_mc()
				CALL inizializza_walkers('opt_Rp')
				CALL inizializza_dati_funzione_onda()
				IF (mpi_myrank==0) WRITE (istring, '(I4.4)'), contatore
				IF (mpi_myrank==0) CALL stampa_file_Rp('reticolo/SR_Rp-'//istring//'.d')
				CALL chiudi_dati_funzione_onda()
				CALL chiudi_walkers()
				CALL chiudi_dati_mc()
				CALL chiudi_dati_fisici()
			END IF
			IF ( ( AV_accu ).AND.(mpi_myrank==0) ) THEN
				ALLOCATE(gRpuu(1:N_hist),gRpuu_new(1:N_hist),gRpud(1:N_hist),gRpud_new(1:N_hist),x_gr(1:N_hist))
				x_gr(1:N_hist)=(/ ( REAL(i,8)*0.5d0*MINVAL(L)/REAL(N_hist,8) , i=1,N_hist ) /)
				gRpuu=0.d0
				gRpuu_new=0.d0
				gRpud=0.d0
				gRpud_new=0.d0
				CALL calcola_g_Rp(p0(N-num_coord_Rp+1:N),gRpuu,gRpud)
			END IF
		END IF
		
		ALLOCATE(s_kl(1:N,1:N), f_k(1:N))
		
		IF (N-num_coord_Rp>0) lambda2=lambda2*DSQRT((DOT_PRODUCT(p0(1:N-num_coord_Rp),p0(1:N-num_coord_Rp)))/REAL(N-num_coord_Rp,8))
		IF (num_coord_Rp>0) lambda2_Rp=lambda2_Rp*DSQRT(num_coord_Rp*MINVAL(L)*0.01d0)
		old_step=0.d0
		old_p0(1:N,1)=p0(1:N)
		old_p0(1:N,2)=p0(1:N)
		old_p0(1:N,3)=p0(1:N)
		old_p0(1:N,4)=p0(1:N)
		old_p0(1:N,5)=p0(1:N)
		old_E=0.d0
		Pdist_min=0.d0
		flag_loop=.TRUE.
		loop_stop=.FALSE.
		IF ( AV_accu ) THEN
			flag_AV_wf=.TRUE.
		ELSE
			flag_AV_wf=.FALSE.
		END IF
				
		DO WHILE (flag_loop)
			contatore=contatore+1
			CALL SR_campiona_wf_attuale(id,p0,N,energia,accettabile)
			
			IF (mpi_myrank==0) WRITE (8, *), cont_step, energia(1:2)
			IF (mpi_myrank==0) CLOSE (8)
			IF (mpi_myrank==0) OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='OLD', POSITION='APPEND')
			IF (mpi_myrank==0) WRITE (9, *), p0
			IF (mpi_myrank==0) CLOSE (9)
			IF (mpi_myrank==0) OPEN (UNIT=9, FILE='ottimizzazione/SR_var_parameters.dat', STATUS='OLD', POSITION='APPEND')
			IF ((.NOT. accettabile).AND.(mpi_myrank==0)) STOP 'Campionamento non accettabile &
			  [ module_variational_opt.f90 > stochastic_reconfiguration ]'
			DO i = 5, 2, -1
				old_E(1:2,i)=old_E(1:2,i-1)
			END DO
			old_E(1:2,1)=energia(1:2)
			
			IF ((energia(1)<old_Emin(1)).OR.(contatore==1)) THEN
				old_Emin(1:2)=energia(1:2)
				p0_minimo(1:N)=p0(1:N)
				flag_trovato_minimo=.TRUE.
				IF (mpi_myrank==0) PRINT * , 'VAR_OPT: Trovato un nuovo minimo assoluto.'
				IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: Trovato un nuovo minimo assoluto.'
				IF (stampa_dati_funzione_onda) THEN
					CALL inizializza_dati_fisici()
					CALL inizializza_dati_mc()
					CALL inizializza_walkers('opt_Rp')
					CALL inizializza_dati_funzione_onda()
					CALL setta_parametri(p0,N)
					IF (mpi_myrank==0) CALL stampa_file_dati_funzione_onda('ottimizzazione/OPT_wf.d')
					IF ((opt_Rp).AND.(mpi_myrank==0)) CALL stampa_file_Rp('reticolo/OPT_Rp.d')
					CALL chiudi_dati_funzione_onda()
					CALL chiudi_walkers()
					CALL chiudi_dati_mc()
					CALL chiudi_dati_fisici()
				END IF
				CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
			ELSE
				flag_trovato_minimo=.FALSE.
				flag_AV_wf=.TRUE.
			END IF
			
			IF ((flag_disk .AND. mpi_myrank==0) .OR. (.NOT. flag_disk)) THEN
				Is_kl=s_kl
				CALL DGETRF( N, N, Is_kl, N, pvt, info )
				CALL DGETRI( N, Is_kl, N, pvt, work, 3*N, info )
				!!! per evitare ill-conditions
				CALL DSYTRD( 'U', N, s_kl, N, D_tr, E_tr, TAU_tr, work, 3*N, info )
				CALL DSTERF( N, D_tr, E_tr, info )
				eps=0.d0
				DO i = 1, N, 1
					eps=MAX(D_tr(i),eps)
				END DO
				eps=eps*0.001
				DO i = 1, N, 1
					s_kl(i,i)=s_kl(i,i)+eps
				END DO	
						
				dp=0.d0
				DO j = 1, N, 1
					DO i = 1, N, 1
						dp(j)=dp(j)+Is_kl(j,i)*f_k(i)
					END DO
				END DO
				vec_app=dp/DSQRT(DOT_PRODUCT(dp,dp))
				
				IF ( contatore==1 ) THEN
					dp_iniz=dp
				END IF
				
				IF ((contatore>2).AND.(N-num_coord_Rp>0)) THEN       !parametri variazionali wf
					flambda2=(1.d0+0.1d0*DOT_PRODUCT(vec_app(1:N-num_coord_Rp),vec_app_old(1:N-num_coord_Rp)))
					IF (contatore<5) flambda2=(1.d0+0.50d0*DOT_PRODUCT(vec_app(1:N-num_coord_Rp),vec_app_old(1:N-num_coord_Rp)))
					lambda2=lambda2*flambda2
					IF (flambda2>1.d0) THEN
						loop_stop(1)=.FALSE.
					ELSE
						loop_stop(1)=.TRUE.
					END IF
					IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [PS] cambio lambda2 di', flambda2
					IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: [PS] cambio lambda2 di', flambda2
				END IF
				IF ((contatore>2).AND.(num_coord_Rp>0)) THEN       !parametri Rp
					flambda2_Rp=(1.d0+0.1d0*DOT_PRODUCT(vec_app(N-num_coord_Rp+1:N),vec_app_old(N-num_coord_Rp+1:N)))
					IF (contatore<5) flambda2_Rp=(1.d0+0.50d0*DOT_PRODUCT(vec_app(N-num_coord_Rp+1:N),vec_app_old(N-num_coord_Rp+1:N)))
					lambda2_Rp=lambda2_Rp*flambda2_Rp
					IF (flambda2_Rp>1.d0) THEN
						loop_stop(1)=.FALSE.
					ELSE
						loop_stop(1)=.TRUE.
					END IF
					IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [PS] cambio lambda2_Rp di', flambda2_Rp
					IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: [PS] cambio lambda2_Rp di', flambda2_Rp
				END IF
				vec_app_old=vec_app
				
				IF ((contatore>3).AND.(N-num_coord_Rp>0)) THEN          !parametri variazionali wf
					dist_poss=SUM(old_step(1:3))
					dist_perc=DSQRT(DOT_PRODUCT(old_p0(1:N-num_coord_Rp,1)-old_p0(1:N-num_coord_Rp,3), &
					                            old_p0(1:N-num_coord_Rp,1)-old_p0(1:N-num_coord_Rp,3)))
					flambda2=(0.85d0+0.3d0*dist_perc/dist_poss)
					lambda2=lambda2*flambda2
					IF (flambda2>1.d0) THEN
						loop_stop(2)=.FALSE.
					ELSE
						loop_stop(2)=.TRUE.
					END IF
					IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DIST3] cambio lambda2 di', flambda2
					IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: [DIST3] cambio lambda2 di', flambda2
				END IF
				IF ((contatore>3).AND.(num_coord_Rp>0)) THEN       !parametri Rp
					dist_poss=SUM(old_step_Rp(1:3))
					dist_perc=DSQRT(DOT_PRODUCT(old_p0(N-num_coord_Rp+1:N,1)-old_p0(N-num_coord_Rp+1:N,3), &
					                            old_p0(N-num_coord_Rp+1:N,1)-old_p0(N-num_coord_Rp+1:N,3)))
					flambda2_Rp=(0.85d0+0.3d0*dist_perc/dist_poss)
					lambda2_Rp=lambda2_Rp*flambda2_Rp
					IF (flambda2_Rp>1.d0) THEN
						loop_stop(2)=.FALSE.
					ELSE
						loop_stop(2)=.TRUE.
					END IF
					IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DIST3] cambio lambda2_Rp di', flambda2_Rp
					IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: [DIST3] cambio lambda2_Rp di', flambda2_Rp
				END IF
				
				IF ((contatore>5).AND.(N-num_coord_Rp>0)) THEN          !parametri variazionali wf
					dist_poss_5=SUM(old_step(1:5))
					dist_perc_5=DSQRT(DOT_PRODUCT(old_p0(1:N-num_coord_Rp,1)-old_p0(1:N-num_coord_Rp,5), &
					                              old_p0(1:N-num_coord_Rp,1)-old_p0(1:N-num_coord_Rp,5)))
					flambda2=(0.75d0+0.50d0*dist_perc_5/dist_poss_5)
					lambda2=lambda2*flambda2
					IF (flambda2>1.d0) THEN
						loop_stop(3)=.FALSE.
					ELSE
						loop_stop(3)=.TRUE.
					END IF
					!IF ( mpi_myrank==0 ) THEN
					!	PRINT * , 'dist_poss_5 = ', dist_poss_5
					!	PRINT * , 'old_step(1:5) = ', old_step(1:5)
					!	PRINT * , 'dist_perc_5 = ', dist_perc_5
					!END IF
					IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DIST5] cambio lambda2 di', flambda2
					IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: [DIST5] cambio lambda2 di', flambda2
				END IF
				IF ((contatore>5).AND.(num_coord_Rp>0)) THEN       !parametri Rp
					dist_poss_5=SUM(old_step_Rp(1:5))
					dist_perc_5=DSQRT(DOT_PRODUCT(old_p0(N-num_coord_Rp+1:N,1)-old_p0(N-num_coord_Rp+1:N,5), &
					                              old_p0(N-num_coord_Rp+1:N,1)-old_p0(N-num_coord_Rp+1:N,5)))
					flambda2_Rp=(0.75d0+0.50d0*dist_perc_5/dist_poss_5)
					lambda2_Rp=lambda2_Rp*flambda2_Rp
					IF (flambda2_Rp>1.d0) THEN
						loop_stop(3)=.FALSE.
					ELSE
						loop_stop(3)=.TRUE.
					END IF
					IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DIST5] cambio lambda2_Rp di', flambda2_Rp
					IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: [DIST5] cambio lambda2_Rp di', flambda2_Rp
				END IF
				
				IF (contatore>5) THEN
					DO i = 1, 4, 1
						derivE(i)=((old_E(1,i)-old_E(1,i+1)))/old_step(i)
						dummy(1)=derivE(i)-DSIGN(1.d0,derivE(i))*(old_E(2,i)+old_E(2,i+1))/old_step(i)
						IF (derivE(i)>0) THEN
							derivE(i)=MAX(dummy(1),0.d0)
						ELSE
							derivE(i)=MIN(dummy(1),0.d0)
						END IF
						derivE(i)=REAL(5-i,8)*derivE(i)
					END DO
				END IF
				IF (contatore>5) THEN
					av_errore=SUM(old_E(2,1:3))/3.d0
					derE=SUM(derivE(1:2))/7.d0
					nextdE=DABS(DABS(derE)-av_errore)
					flambda2=MIN(1.10d0,MAX(0.900d0,nextdE/av_errore))
					IF (flambda2>1.d0) THEN
						loop_stop(4)=.FALSE.
					ELSE
						loop_stop(4)=.TRUE.
					END IF
				END IF
				IF (contatore>5) THEN
					av_errore=SUM(old_E(2,1:5))/5.d0
					derE=SUM(derivE(1:4))/10.d0
					nextdE=DABS(DABS(derE)-av_errore)
					flambda2=MIN(1.150d0,MAX(0.85d0,nextdE/av_errore))
					IF (flambda2>1.d0) THEN
						loop_stop(5)=.FALSE.
					ELSE
						loop_stop(5)=.TRUE.
					END IF
				END IF
				IF (contatore>5) THEN
					av_errore=SUM(old_E(2,1:5))/5.d0
					derE=0.5d0*(old_E(1,1)+old_E(1,2)-old_E(1,5)-old_E(1,4)) / &
					  DSQRT(DOT_PRODUCT(old_p0(1:N,1)-old_p0(1:N,5),old_p0(1:N,1)-old_p0(1:N,5)))
					nextdE=2.d0*DABS(derE*DSQRT(DOT_PRODUCT(old_p0(1:N,1)-old_p0(1:N,5),old_p0(1:N,1)-old_p0(1:N,5)))*0.2d0)
					flambda2=MIN(1.5d0,MAX(0.5d0,nextdE/av_errore))
					IF (flambda2>1.d0) THEN
						loop_stop(6)=.FALSE.
					ELSE
						loop_stop(6)=.TRUE.
					END IF
				END IF
				IF (contatore>5) THEN
					IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DER] Possibile fermarsi? ', loop_stop(4:6)
					IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: [DER] Possibile fermarsi? ', loop_stop(4:6)
				END IF
				
				IF ( N-num_coord_Rp>0 ) THEN
					lambda=0.25d0
					lambda=lambda/DSQRT(DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp)))
				END IF
				IF ( num_coord_Rp>0 ) THEN
					lambda_Rp=0.25d0
					lambda_Rp=lambda_Rp/DSQRT(DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N)))
				END IF
								
				!controllo se posso uscire dal loop
				flag_loop=.TRUE.
				IF (loop_stop(4).AND.(loop_stop(5)).AND.(loop_stop(6))) flag_loop=.FALSE.
				DO i = 1, N, 1        !evita cambi di parametri troppo bruschi
					IF ( i<=N-num_coord_Rp ) THEN
						vec_app(i)=dp(i)*lambda*lambda2/p0(i)       
						IF (DABS(vec_app(i))>MAX_VAR_ACC) THEN
							IF (mpi_myrank==0) PRINT '(A46,I3,A12,F9.3,A3,F9.3)', &
							  'VAR_OPT: lambda2 troppo grande per parametro ', i, '. Riduco da ', lambda2, ' a ', &
							  DABS(p0(i)*MAX_VAR_ACC/(dp(i)*lambda))
							IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, '(A46,I3,A12,F9.3,A3,F9.3)'), &
							  'VAR_OPT: lambda2 troppo grande per parametro ', i, '. Riduco da ', lambda2, ' a ', &
							  DABS(p0(i)*MAX_VAR_ACC/(dp(i)*lambda))
							lambda2=DABS(p0(i)*MAX_VAR_ACC/(dp(i)*lambda))
						END IF
						IF ((DABS(vec_app(i))>MIN_VAR_NEC).AND.(i<=N-num_coord_Rp)) flag_loop=.TRUE.
					ELSE          !evita cambi di Rp troppo bruschi
						vec_app(i)=dp(i)*lambda_Rp*lambda2_Rp/MINVAL(L)   
						IF (DABS(vec_app(i))>MAX_VAR_PROT) THEN
							IF (mpi_myrank==0) PRINT '(A49,I3,A12,F9.3,A3,F9.3)', &
							  'VAR_OPT: lambda2_Rp troppo grande per parametro ', i, '. Riduco da ', lambda2_Rp, ' a ', &
							  DABS(MINVAL(L)*MAX_VAR_PROT/(dp(i)*lambda_Rp))
							IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, '(A49,I3,A12,F9.3,A3,F9.3)'), &
							  'VAR_OPT: lambda2_Rp troppo grande per parametro ', i, '. Riduco da ', lambda2_Rp, ' a ', &
							  DABS(MINVAL(L)*MAX_VAR_PROT/(dp(i)*lambda_Rp))
							lambda2_Rp=DABS(MINVAL(L)*MAX_VAR_PROT/(dp(i)*lambda_Rp))
						END IF
						IF ((DABS(vec_app(i))>MIN_VAR_NEC).AND.(i<=N-num_coord_Rp)) flag_loop=.TRUE.
					END IF
				END DO
				IF (flag_loop) THEN       !nel caso le derivate non facciano fermare ma i parametri siano quasi congelati
					flag_freeze=.TRUE.
					DO i = 1, N-num_coord_Rp, 1
						IF (DABS(vec_app(i))>MIN_VAR_ACC) flag_freeze=.FALSE.
					END DO
					flag_loop=(.NOT. flag_freeze)
				END IF
				IF (contatore<=5) flag_loop=.TRUE.
				
				IF (.NOT. auto_stop) flag_loop=.TRUE.
												
				!determino e stampo lo spostamento variazionale
				IF (flag_loop) THEN
					IF ( N-num_coord_Rp>0 ) THEN
						cont_step=cont_step+lambda*lambda2*DSQRT(DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp)))
						old_step(5)=old_step(4)
						old_step(4)=old_step(3)
						old_step(3)=old_step(2)
						old_step(2)=old_step(1)
						old_step(1)=lambda*lambda2*DSQRT(DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp)))
						IF (mpi_myrank==0) PRINT * , 'VAR_OPT: fattore dp: ', &
						  DSQRT(DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp))/ &
						  DOT_PRODUCT(dp_iniz(1:N-num_coord_Rp),dp_iniz(1:N-num_coord_Rp)))
						dp(1:N-num_coord_Rp)=dp(1:N-num_coord_Rp)*lambda*lambda2					
					END IF
					IF ( num_coord_Rp>0 ) THEN
						cont_step=cont_step+lambda_Rp*lambda2_Rp*DSQRT(DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N)))
						old_step_Rp(5)=old_step_Rp(4)
						old_step_Rp(4)=old_step_Rp(3)
						old_step_Rp(3)=old_step_Rp(2)
						old_step_Rp(2)=old_step_Rp(1)
						old_step_Rp(1)=lambda_Rp*lambda2_Rp*DSQRT(DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N)))
						IF (mpi_myrank==0) PRINT * , 'VAR_OPT: fattore dp_Rp: ', &
						  DSQRT(DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N)) / &
						  DOT_PRODUCT(dp_iniz(N-num_coord_Rp+1:N),dp_iniz(N-num_coord_Rp+1:N)))
						dp(N-num_coord_Rp+1:N)=dp(N-num_coord_Rp+1:N)*lambda_Rp*lambda2_Rp
					END IF
					
					IF ((mpi_myrank==0).AND.(N-num_coord_Rp>0)) THEN
						PRINT '(A54,F9.3,A2)', ' VAR_OPT: seguendo SR, cambio i parametri:   (lambda2=',lambda2,') '!'  (lambda3=', lambda3,')'
						WRITE (stringa, '(I4.4)'), N-num_coord_Rp
						PRINT '(1X,'//stringa//'(F9.3,1X))' , p0(1:N-num_coord_Rp)
						PRINT '(1X,'//stringa//'(F9.3,1X))' , p0(1:N-num_coord_Rp)+dp(1:N-num_coord_Rp)
						IF (flag_output) THEN
							WRITE (7, '(A54,F9.3,A2)'), &
							  ' VAR_OPT: seguendo SR, cambio i parametri:   (lambda2=',lambda2,') '!'  (lambda3=', lambda3,')'
							WRITE (7, '(1X,'//stringa//'(F9.3,1X))'), p0(1:N-num_coord_Rp)
							WRITE (7, '(1X,'//stringa//'(F9.3,1X))'), p0(1:N-num_coord_Rp)+dp(1:N-num_coord_Rp)
						END IF
					END IF
					IF ((mpi_myrank==0).AND.(num_coord_Rp>0)) THEN
						PRINT '(A54,F9.3,A2)', ' VAR_OPT: seguendo SR, muovo i protoni:   (lambda2_Rp=',lambda2_Rp,') '!'  (lambda3=', lambda3,')'
						!WRITE (stringa, '(I4.4)'), num_coord_Rp
						!DO i = 1, N_part, 1
						!	PRINT '(1X,A8,I3.3,10X,3(F9.3,5X))' , 'Protone ', i, dp(N-num_coord_Rp+i*3-2:N-num_coord_Rp+i*3)
						!END DO
						IF (flag_output) THEN
							WRITE (7, '(A54,F9.3,A2)'), &
							  ' VAR_OPT: seguendo SR, muovo i protoni:   (lambda2_Rp=',lambda2_Rp,') '!'  (lambda3=', lambda3,')'
							!DO i = 1, N_part, 1
							!	WRITE (7, '(1X,A8,I3.3,10X,3(F9.3,5X))'), 'Protone ', i, dp(N-num_coord_Rp+i*3-2:N-num_coord_Rp+i*3)
							!END DO
						END IF
					END IF
				END IF
			END IF
			CALL MPI_BCAST(dp(1:N),N,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
			CALL MPI_BCAST(flag_loop,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_ierr)	
					
			IF (flag_loop) THEN
				p0=p0+dp
				IF ( flag_AV_wf ) THEN
					AV_P=AV_P+p0
					AV_cont=AV_cont+1
				ELSE
					AV_P=p0
					AV_cont=1
				END IF
			END IF
			old_p0(1:N,5)=old_p0(1:N,4)
			old_p0(1:N,4)=old_p0(1:N,3)
			old_p0(1:N,3)=old_p0(1:N,2)
			old_p0(1:N,2)=old_p0(1:N,1)
			old_p0(1:N,1)=p0(1:N)
			
			!stampa gli ultimi dati variazionali SR
			IF (stampa_dati_funzione_onda) THEN
				CALL inizializza_dati_fisici()
				CALL inizializza_dati_mc()
				CALL inizializza_walkers('opt_Rp')
				CALL inizializza_dati_funzione_onda()
				CALL setta_parametri(p0,N)
				IF (mpi_myrank==0) CALL stampa_file_dati_funzione_onda('ottimizzazione/SR_wf.d')
				IF ((mpi_myrank==0).AND.(opt_Rp)) WRITE (istring, '(I4.4)'), contatore
				IF ((mpi_myrank==0).AND.(opt_Rp)) CALL stampa_file_Rp('reticolo/SR_Rp-'//istring//'.d')
				CALL chiudi_dati_fisici()
				CALL chiudi_dati_mc()
				CALL chiudi_walkers()
				CALL chiudi_dati_funzione_onda()
			END IF
			
			!stampa i dati variazionali AVERAGE
			IF ((stampa_dati_funzione_onda)) THEN
				CALL inizializza_dati_fisici()
				CALL inizializza_dati_mc()
				CALL inizializza_walkers('opt_Rp')
				CALL inizializza_dati_funzione_onda()
				AV_P=AV_P/REAL(AV_cont,8)
				CALL setta_parametri(AV_P,N)
				AV_P=AV_P*REAL(AV_cont,8)
				IF (mpi_myrank==0) CALL stampa_file_dati_funzione_onda('ottimizzazione/AV_wf.d')
				IF (mpi_myrank==0) PRINT *, 'VAR_OPT: salvato file AV_wf.d, considerando gli ultimi ', AV_cont, 'punti.'
				IF ((flag_output).AND.(mpi_myrank==0)) THEN
					WRITE (7, *), ' VAR_OPT: salvato file AV_wf.d, considerando gli ultimi ', AV_cont, 'punti.'
				END IF
				IF ( opt_Rp ) THEN
					IF (mpi_myrank==0) CALL stampa_file_Rp('reticolo/AV_Rp.d')
				END IF
				CALL chiudi_dati_fisici()
				CALL chiudi_dati_mc()
				CALL chiudi_walkers()
				CALL chiudi_dati_funzione_onda()
			END IF
			
			!accumulo la g(r) di Rp e la trascrivo
			IF (( AV_accu ).AND.(mpi_myrank==0)) THEN
				CALL calcola_g_Rp(p0(N-num_coord_Rp+1:N),gRpuu_new,gRpud_new)
				gRpuu=gRpuu+gRpuu_new
				gRpud=gRpud+gRpud_new
				CALL salva_vettore_dati_dat(x_gr,gRpuu/REAL(contatore+1,8),N_hist,'reticolo/AV-guu_Rp')
				CALL salva_vettore_dati_dat(x_gr,gRpud/REAL(contatore+1,8),N_hist,'reticolo/AV-gud_Rp')
				CALL salva_vettore_dati_dat(x_gr,(gRpuu+gRpud)/REAL(contatore+1,8),N_hist,'reticolo/AV-g_Rp')
			END IF
			
			CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
		END DO
		
		!salvo la g(r) dei Rp
		IF ( ( AV_accu ).AND.(mpi_myrank==0) ) THEN
			DEALLOCATE(gRpuu,gRpuu_new,gRpud,gRpud_new,x_gr)
		END IF
		
		E_min(1:2)=old_Emin(1:2)
		parametri_var=p0_minimo
		IF (mpi_myrank==0) CLOSE (8)
		DEALLOCATE(s_kl, f_k)

	END SUBROUTINE stochastic_reconfiguration

	SUBROUTINE pure_stochastic_reconfiguration(p0,N,auto_stop)
		USE VMC
		USE dati_fisici
		USE funzione_onda
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: MIN_VAR_NEC=0.025, MIN_LAMBDA2=0.5, MIN_VAR_ACC=0.001, MAX_VAR_ACC=0.25
		LOGICAL, INTENT(IN) :: auto_stop
		INTEGER (KIND=4), INTENT(IN) :: N
		LOGICAL :: flag_loop, loop_stop(1:6), accettabile, flag_file, flag_trovato_minimo, flag_freeze
		CHARACTER(LEN=2) :: id, ida, idb, idc
		CHARACTER(LEN=4) :: stringa
		INTEGER :: i, j, info, contatore, IO, i_orbit, AV_cont
		REAL (KIND=8) :: lambda, eps, boundd, dummy(1:4), dummy2(0:N)
		REAL (KIND=8) :: p0(1:N), dp_iniz(1:N), AV_P(1:N)
		REAL (KIND=8) :: dp(1:N), Is_kl(1:N,1:N), pvt(1:N), cont_step
		REAL (KIND=8) :: work(1:3*N), D_tr(1:N), E_tr(1:N-1), TAU_tr(1:N-1)
		REAL (KIND=8) :: energia(1:2), vec_app(1:N), vec_app_old(1:N), flambda2
		REAL (KIND=8) :: old_step(1:5), old_p0(1:N,1:5), old_E(1:2,1:5), derivE(1:4)
		REAL (KIND=8) :: dist_poss, dist_perc, dist_poss_5, dist_perc_5, av_errore, derE, nextdE
		REAL (KIND=8) :: old_Emin(1:2), p0_minimo(1:N), Pdist_min(1:5)
		REAL (KIND=8) :: a, Ea(1:2), b, Eb(1:2), c, Ec(1:2), first_step
		id='SR'
				
		INQUIRE(FILE='ottimizzazione/SR_energies.dat',EXIST=flag_file)
		IF ((flag_file).AND.(mpi_myrank==0)) THEN
			OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='OLD')
			IO=0
			DO WHILE (IO>=0)
				READ (8, *, IOSTAT = IO), dummy(1:3)
				IF (IO>=0) dummy(4)=dummy(1) 
			END DO
			cont_step=dummy(4)
			CLOSE(8)
			OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='OLD', POSITION='APPEND')
			IF (mpi_myrank==0) PRINT * , 'VAR_OPT: esisteva giá un file, parto da ', cont_step, '.'
		ELSE IF (mpi_myrank==0) THEN
			OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='NEW')
			cont_step=0.d0
			IF (mpi_myrank==0) PRINT * , 'VAR_OPT: parto da zero.'
		END IF
		
		OPEN (UNIT=9, FILE='ottimizzazione/SR_var_parameters.dat', STATUS='UNKNOWN', POSITION='APPEND')
		
		ALLOCATE(s_kl(1:N,1:N), f_k(1:N))
		contatore=0
		lambda=-666.d0
		old_step=0.d0
		old_p0(1:N,1)=p0(1:N)
		old_p0(1:N,2)=p0(1:N)
		old_p0(1:N,3)=p0(1:N)
		old_p0(1:N,4)=p0(1:N)
		old_p0(1:N,5)=p0(1:N)
		old_E=0.d0
		Pdist_min=0.d0
		flag_loop=.TRUE.
		loop_stop=.FALSE.
		DO WHILE (flag_loop)
			contatore=contatore+1
			IF (stampa_dati_funzione_onda) THEN
				CALL inizializza_dati_fisici()
				CALL inizializza_dati_funzione_onda()
				CALL setta_parametri(p0,N)
				CALL stampa_file_dati_funzione_onda('ottimizzazione/SR_wf.d')
				CALL chiudi_dati_fisici()
				CALL chiudi_dati_funzione_onda()
			END IF
			CALL SR_campiona_wf_attuale(id,p0,N,energia,accettabile)
			IF (mpi_myrank==0) WRITE (8, *), cont_step, energia(1:2)
			IF (mpi_myrank==0) CLOSE (8)
			IF (mpi_myrank==0) OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='OLD', POSITION='APPEND')
			IF (mpi_myrank==0) WRITE (9, *), p0
			IF (mpi_myrank==0) CLOSE (9)
			IF (mpi_myrank==0) OPEN (UNIT=9, FILE='ottimizzazione/SR_var_parameters.dat', STATUS='OLD', POSITION='APPEND')
			IF ((.NOT. accettabile).AND.(mpi_myrank==0)) STOP 'Campionamento non accettabile &
			  [ module_variational_opt.f90 > stochastic_reconfiguration ]'
			DO i = 5, 2, -1
				old_E(1:2,i)=old_E(1:2,i-1)
			END DO
			old_E(1:2,1)=energia(1:2)
			IF ((energia(1)<old_Emin(1)).OR.(contatore==1)) THEN
				old_Emin(1:2)=energia(1:2)
				p0_minimo(1:N)=p0(1:N)
				flag_trovato_minimo=.TRUE.
				IF (mpi_myrank==0) PRINT * , 'VAR_OPT: Trovato un nuovo minimo assoluto.'
				IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: Trovato un nuovo minimo assoluto.'
				IF (stampa_dati_funzione_onda) THEN
					CALL inizializza_dati_fisici()
					CALL inizializza_dati_funzione_onda()
					CALL setta_parametri(p0,N)
					CALL stampa_file_dati_funzione_onda('ottimizzazione/OPT_wf.d')
					CALL chiudi_dati_fisici()
					CALL chiudi_dati_funzione_onda()
				END IF
				!contatore=1
			ELSE
				flag_trovato_minimo=.FALSE.
			END IF
			
			Is_kl=s_kl
			CALL DGETRF( N, N, Is_kl, N, pvt, info )
			CALL DGETRI( N, Is_kl, N, pvt, work, 3*N, info )
			CALL DSYTRD( 'U', N, s_kl, N, D_tr, E_tr, TAU_tr, work, 3*N, info )
			CALL DSTERF( N, D_tr, E_tr, info )
			eps=0.d0
			DO i = 1, N, 1
				eps=MAX(D_tr(i),eps)
			END DO
			eps=eps*0.001
			DO i = 1, N, 1
				s_kl(i,i)=s_kl(i,i)+eps
			END DO			
			dp=0.d0
			DO j = 1, N, 1
				DO i = 1, N, 1
					dp(j)=dp(j)+Is_kl(j,i)*f_k(i)
				END DO
			END DO
			!vec_app=dp/DSQRT(DOT_PRODUCT(dp,dp))
			
			IF ( contatore==1 ) THEN
				dp_iniz=dp
			END IF
			
			IF (contatore>5) THEN
				DO i = 1, 4, 1
					derivE(i)=((old_E(1,i)-old_E(1,i+1)))/old_step(i)
					dummy(1)=derivE(i)-DSIGN(1.d0,derivE(i))*(old_E(2,i)+old_E(2,i+1))/old_step(i)
					IF (derivE(i)>0) THEN
						derivE(i)=MAX(dummy(1),0.d0)
					ELSE
						derivE(i)=MIN(dummy(1),0.d0)
					END IF
					derivE(i)=REAL(5-i,8)*derivE(i)
				END DO
			END IF
			IF (contatore>5) THEN
				av_errore=SUM(old_E(2,1:3))/3.d0
				derE=SUM(derivE(1:2))/7.d0
				nextdE=DABS(DABS(derE)-av_errore)
				flambda2=MIN(1.10d0,MAX(0.900d0,nextdE/av_errore))
				IF (flambda2>1.d0) THEN
					loop_stop(4)=.FALSE.
				ELSE
					loop_stop(4)=.TRUE.
				END IF
			END IF
			IF (contatore>5) THEN
				av_errore=SUM(old_E(2,1:5))/5.d0
				derE=SUM(derivE(1:4))/10.d0
				nextdE=DABS(DABS(derE)-av_errore)
				flambda2=MIN(1.150d0,MAX(0.85d0,nextdE/av_errore))
				IF (flambda2>1.d0) THEN
					loop_stop(5)=.FALSE.
				ELSE
					loop_stop(5)=.TRUE.
				END IF
			END IF
			IF (contatore>5) THEN
				av_errore=SUM(old_E(2,1:5))/5.d0
				derE=0.5d0*(old_E(1,1)+old_E(1,2)-old_E(1,5)-old_E(1,4)) / &
				  DSQRT(DOT_PRODUCT(old_p0(1:N,1)-old_p0(1:N,5),old_p0(1:N,1)-old_p0(1:N,5)))
				nextdE=2.d0*DABS(derE*DSQRT(DOT_PRODUCT(old_p0(1:N,1)-old_p0(1:N,5),old_p0(1:N,1)-old_p0(1:N,5)))*0.2d0)
				flambda2=MIN(1.5d0,MAX(0.5d0,nextdE/av_errore))
				IF (flambda2>1.d0) THEN
					loop_stop(6)=.FALSE.
				ELSE
					loop_stop(6)=.TRUE.
				END IF
			END IF
			IF (contatore>5) THEN
				IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DER] Possibile fermarsi? ', loop_stop(4:6)
				IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: [DER] Possibile fermarsi? ', loop_stop(4:6)
			END IF
			
			IF (lambda==-666.d0) THEN
				CALL cambia_verbosity_VMC(.FALSE.)
				first_step=0.1d0
				CALL bracket_min_multidim(N,p0,dp,first_step,a,Ea,ida,b,Eb,idb,c,Ec,idc)
				CALL parabgold_search_multidim(N,p0,dp,a,Ea,ida,b,Eb,idb,c,Ec,idc)
				lambda=b
				IF ( mpi_myrank==0 ) PRINT * , 'VAR_OPT: LAMBDA= ', lambda
				IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *), 'VAR_OPT: LAMBDA= ', lambda
				CALL cambia_verbosity_VMC(.TRUE.)
			END IF
							
			!controllo se posso uscire dal loop
			flag_loop=.TRUE.
			IF (loop_stop(4).AND.(loop_stop(5)).AND.(loop_stop(6))) flag_loop=.FALSE.
			DO i = 1, N, 1
				vec_app(i)=dp(i)*lambda/p0(i)        !!!!DA RIPRISTINARE!!!!!!
				IF (DABS(vec_app(i))>MIN_VAR_NEC) flag_loop=.TRUE.
			END DO
			IF (flag_loop) THEN       !nel caso le derivate non facciano fermare ma i parametri siano quasi congelati
				flag_freeze=.TRUE.
				DO i = 1, N, 1
					IF (DABS(vec_app(i))>MIN_VAR_ACC) flag_freeze=.FALSE.
				END DO
				flag_loop=(.NOT. flag_freeze)
			END IF
			IF (contatore<=5) flag_loop=.TRUE.
			
			IF (.NOT. auto_stop) flag_loop=.TRUE. 
			
			!determino e stampo lo spostamento variazionale
			IF (flag_loop) THEN
				cont_step=cont_step+lambda*DSQRT(DOT_PRODUCT(dp,dp))
				old_step(5)=old_step(4)
				old_step(4)=old_step(3)
				old_step(3)=old_step(2)
				old_step(2)=old_step(1)
				old_step(1)=DSQRT(DOT_PRODUCT(dp,dp))
				IF (mpi_myrank==0) PRINT * , 'VAR_OPT: fattore dp: ', DSQRT(DOT_PRODUCT(dp,dp)/DOT_PRODUCT(dp_iniz,dp_iniz))
				
				IF (mpi_myrank==0) THEN
					!OPEN(UNIT=44, FILE='ottimizzazione/SR_lambda.d', STATUS='UNKNOWN', POSITION='APPEND')
					!WRITE(UNIT=44, FMT=*) lambda, DSQRT(DOT_PRODUCT(dp,dp)), lambda*DSQRT(DOT_PRODUCT(dp,dp))
					!CLOSE(44)
					PRINT '(A54,F9.3,A13,F9.3,A1)', ' VAR_OPT: seguendo SR, cambio i parametri:'
					WRITE (stringa, '(I4.4)'), N
					PRINT '(1X,'//stringa//'(F9.3,1X))' , p0
					PRINT '(1X,'//stringa//'(F9.3,1X))' , p0+lambda*dp
					IF (flag_output) THEN
						WRITE (7, '(A54,F9.3,A13,F9.3,A1)'), &
						  ' VAR_OPT: seguendo SR, cambio i parametri:'
						WRITE (7, '(1X,'//stringa//'(F9.3,1X))'), p0
						WRITE (7, '(1X,'//stringa//'(F9.3,1X))'), p0+lambda*dp
					END IF
				END IF
			END IF
			!CALL MPI_BCAST(dp(1:N),N,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
			!CALL MPI_BCAST(flag_loop,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_ierr)
			IF (flag_loop) THEN
				p0=p0+lambda*dp
				IF ((loop_stop(4)).AND.(loop_stop(5)).AND.(loop_stop(6))) THEN
					AV_cont=AV_cont+1
					AV_P=AV_P+p0
				ELSE
					AV_cont=0
					AV_P=0.d0
				END IF
			END IF
			old_p0(1:N,5)=old_p0(1:N,4)
			old_p0(1:N,4)=old_p0(1:N,3)
			old_p0(1:N,3)=old_p0(1:N,2)
			old_p0(1:N,2)=old_p0(1:N,1)
			old_p0(1:N,1)=p0(1:N)
		END DO
		
		IF (stampa_dati_funzione_onda) THEN
			CALL inizializza_dati_fisici()
			CALL inizializza_dati_funzione_onda()
			AV_P=AV_P/REAL(AV_cont,8)
			CALL setta_parametri(AV_P,N)
			CALL stampa_file_dati_funzione_onda('ottimizzazione/AV_wf.d')
			CALL chiudi_dati_fisici()
			CALL chiudi_dati_funzione_onda()
		END IF
		
		E_min(1:2)=old_Emin(1:2)
		parametri_var=p0_minimo
		IF (mpi_myrank==0) CLOSE (8)
		DEALLOCATE(s_kl, f_k)

	END SUBROUTINE pure_stochastic_reconfiguration

	SUBROUTINE bracket_min_multidim(N,r0,dir,first_step,a,fa,ida,b,fb,idb,c,fc,idc)
		IMPLICIT NONE
		LOGICAL :: flag_too_accurate, accettabile
		CHARACTER(LEN=2) :: ida, idb, idc, id
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: r0(1:N)
		REAL (KIND=8) :: first_step, dir(1:N), ra(1:N), rb(1:N), rc(1:N)
		REAL (KIND=8) :: a, b, c, xm, derivata
		REAL (KIND=8) :: fa(1:2), fb(1:2), fc(1:2), fxm(1:2)
		INTEGER (KIND=4) :: contatore, num_val_en
						
		num_val_en=1
		ida='00'
		ra=r0
								
		first_step=first_step/DSQRT(DOT_PRODUCT(dir,dir))
		idb='b'//CHAR(MOD(num_val_en,10)+48)
		b=a+first_step
		rb=r0+b*dir
		CALL controlla_punto(rb)
		IF (mpi_myrank==0) PRINT * , 'Ho controllato il punto'
		CALL calcola_energia(idb,rb,N,fb(1:2),accettabile)
		IF (mpi_myrank==0) PRINT * , 'Ho calcolato la prima energia'
		DO WHILE (.NOT. accettabile)
			first_step=0.5d0*first_step
			b=a+first_step
			rb=r0+b*dir
			CALL controlla_punto(rb)
			CALL calcola_energia(idb,rb,N,fb(1:2),accettabile)
		END DO
		DO WHILE ((DABS(fa(1)-fb(1))<0.25d0*liv_precisione))
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
				IF (flag_output) WRITE (7, *), 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
			END IF
			b=b+first_step
			rb=r0+b*dir
			CALL controlla_punto(rb)
			CALL calcola_energia(idb,rb,N,fb(1:2),accettabile)
			IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 1 &
			  [ module_variational_opt.f90 > bracket_min_multidim ]'
		END DO
		flag_too_accurate=.FALSE.
		IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
		DO WHILE ((fa(2)+fb(2)>0.25d0*DABS(fa(1)-fb(1))) .AND. (.NOT. flag_too_accurate))
			IF (fa(2)>fb(2)) THEN
				CALL controlla_punto(ra)
				CALL aumenta_precisione_energia(ida,ra,N,fa(1:2))
			ELSE
				CALL controlla_punto(rb)
				CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
			END IF
			IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
		END DO
		derivata=(fb(1)-fa(1))/(b-a)
		
		IF (derivata>0.d0) THEN
			c=b
			fc=fb
			idc=idb
			rc=rb
			b=a
			fb=fa
			idb=ida
			rb=ra
			a=b-first_step
			ra=r0+a*dir
			num_val_en=num_val_en+1
			ida='b'//CHAR(MOD(num_val_en,10)+48)
			DO WHILE ((ida==idb) .OR. (ida==idc))
				num_val_en=num_val_en+1
				ida='b'//CHAR(MOD(num_val_en,10)+48)
			END DO
			CALL controlla_punto(ra)
			CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
			DO WHILE (.NOT. accettabile)
				a=(b+a)*0.5d0
				ra=r0+a*dir
				CALL controlla_punto(ra)
				CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
			END DO 
			DO WHILE (DABS(fa(1)-fb(1))<0.5d0*liv_precisione)
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
					IF (flag_output) WRITE (7, *), 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
				END IF
				a=a-first_step
				ra=r0+a*dir
				CALL controlla_punto(ra)
				CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
				IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 2 &
				  [ module_variational_opt.f90 > bracket_min_multidim ]'
			END DO
			flag_too_accurate=.FALSE.
			IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
			DO WHILE ((fa(2)+fb(2)>0.25d0*DABS(fa(1)-fb(1))) .AND. (.NOT. flag_too_accurate))
				IF (fa(2)>fb(2)) THEN
					CALL controlla_punto(ra)
					CALL aumenta_precisione_energia(ida,ra,N,fa(1:2))
				ELSE
					CALL controlla_punto(rb)
					CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
				END IF
				IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
			END DO
			contatore=1
			DO WHILE (fa(1)<fb(1))
				contatore=contatore+1
				c=b
				fc=fb
				idc=idb
				rc=rb
				b=a
				fb=fa
				idb=ida
				rb=ra
				a=a-(2.d0**contatore)*DABS(first_step)
				ra=r0+a*dir
				num_val_en=num_val_en+1
				ida='b'//CHAR(MOD(num_val_en,10)+48)
				DO WHILE ((ida==idb) .OR. (ida==idc))
					num_val_en=num_val_en+1
					ida='b'//CHAR(MOD(num_val_en,10)+48)
				END DO
				CALL controlla_punto(ra)
				CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
				DO WHILE (.NOT. accettabile)
					a=(b+a)*0.5d0
					ra=r0+a*dir
					IF (a<b-(2.d0**(contatore-1))*DABS(first_step)) THEN
						CALL controlla_punto(ra)
						CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
					END IF
				END DO
				flag_too_accurate=.FALSE.
				IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
				DO WHILE ((fa(2)+fb(2)>0.25d0*DABS(fa(1)-fb(1))) .AND. (.NOT. flag_too_accurate))
					IF (fa(2)>fb(2)) THEN
						CALL controlla_punto(ra)
						CALL aumenta_precisione_energia(ida,ra,N,fa(1:2))
					ELSE
						CALL controlla_punto(rb)
						CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
					END IF
					IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
				END DO
			END DO
		ELSE
			c=b+first_step
			rc=r0+c*dir
			num_val_en=num_val_en+1
			idc='b'//CHAR(MOD(num_val_en,10)+48)
			DO WHILE ((idc==ida) .OR. (idc==idb))
				num_val_en=num_val_en+1
				idc='b'//CHAR(MOD(num_val_en,10)+48)
			END DO
			CALL controlla_punto(rc)
			CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
			DO WHILE (.NOT. accettabile)
				c=(b+c)*0.5d0
				rc=r0+c*dir
				CALL controlla_punto(rc)
				CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
			END DO
			DO WHILE (DABS(fb(1)-fc(1))<0.25d0*liv_precisione)
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
					IF (flag_output) WRITE (7, *), 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
				END IF
				!first_step=2.d0*first_step
				c=c+first_step
				rc=r0+c*dir
				CALL controlla_punto(rc)
				CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
				IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 3 &
				  [ module_variational_opt.f90 > bracket_min_multidim ]'
			END DO
			flag_too_accurate=.FALSE.
			IF (DABS(fb(1)-fc(1))<liv_precisione) flag_too_accurate=.TRUE.
			DO WHILE ((fb(2)+fc(2)>0.25d0*DABS(fb(1)-fc(1))) .AND. (.NOT. flag_too_accurate))
				IF (fb(2)>fc(2)) THEN
					CALL controlla_punto(rb)
					CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
				ELSE
					CALL controlla_punto(rc)
					CALL aumenta_precisione_energia(idc,rc,N,fc(1:2))
				END IF
				IF (DABS(fb(1)-fc(1))<liv_precisione) flag_too_accurate=.TRUE.
			END DO
			contatore=1
			DO WHILE (fc(1)<fb(1))
				contatore=contatore+1
				a=b
				fa=fb
				ida=idb
				ra=rb
				b=c
				fb=fc
				idb=idc
				rb=rc
				c=c+(2.d0**contatore)*first_step
				rc=r0+c*dir
				num_val_en=num_val_en+1
				idc='b'//CHAR(MOD(num_val_en,10)+48)
				DO WHILE ((idc==ida) .OR. (idc==idb))
					num_val_en=num_val_en+1
					idc='b'//CHAR(MOD(num_val_en,10)+48)
				END DO
				CALL controlla_punto(rc)
				CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
				IF (.NOT. accettabile) THEN
					c=(b+c)*0.5d0
					rc=r0+c*dir
					IF (c<b+(2.d0**(contatore-1))*DABS(first_step)) THEN
						CALL controlla_punto(rc)
						CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
					END IF
				END IF
				flag_too_accurate=.FALSE.
				IF (DABS(fb(1)-fc(1))<liv_precisione) flag_too_accurate=.TRUE.
				DO WHILE ((fb(2)+fc(2)>0.25d0*DABS(fb(1)-fc(1))) .AND. (.NOT. flag_too_accurate))
					IF (fb(2)>fc(2)) THEN
						CALL controlla_punto(rb)
						CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
					ELSE
						CALL controlla_punto(rc)
						CALL aumenta_precisione_energia(idc,rc,N,fc(1:2))
					END IF
					IF (DABS(fb(1)-fc(1))<liv_precisione) flag_too_accurate=.TRUE.
				END DO
			END DO
		END IF
		IF (mpi_myrank==0) THEN
			DO num_val_en = 1, 10, 1                    !elimino i dati ormai inutili
				id='b'//CHAR(MOD(num_val_en,10)+48)
				IF ((id/=ida) .AND. (id/=idb) .AND. (id/=idc)) CALL SYSTEM ("export ID_TO_ERASE=opt'//id//'; &
				  find estimatori/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \; &
				  find posizioni/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \")
			END DO
		END IF

		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Ho eseguito il bracket con ', INT(contatore), ' passi.',' I tre punti sono '
			PRINT * , '         punti:   ', a, b ,c
			PRINT * , '         energie: ', fa(1), fb(1), fc(1)
			IF (flag_output) THEN
				WRITE (7, *), 'VAR_OPT: Ho eseguito il bracket con ', INT(contatore), ' passi.',' I tre punti sono '
				WRITE (7, *), '         punti:   ', a, b ,c
				WRITE (7, *), '         energie: ', fa(1), fb(1), fc(1)
			END IF
		END IF
				
	END SUBROUTINE bracket_min_multidim

	SUBROUTINE bracket_min_multidim_v2(N,r0,dir,first_step,a,fa,ida,b,fb,idb,c,fc,idc)
		IMPLICIT NONE
		LOGICAL :: flag_too_accurate, accettabile, flag_border
		CHARACTER(LEN=2) :: ida, idb, idc, id
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: r0(1:N)
		REAL (KIND=8) :: first_step, dir(1:N), ra(1:N), rb(1:N), rc(1:N)
		REAL (KIND=8) :: a, b, c, xm, derivata, b_app
		REAL (KIND=8) :: fa(1:2), fb(1:2), fc(1:2), fxm(1:2)
		INTEGER (KIND=4) :: contatore, num_val_en
				
		num_val_en=1
		!ida='00'
		ra=r0
		first_step=first_step/DSQRT(DOT_PRODUCT(dir,dir))
			
		idb='b'//CHAR(MOD(num_val_en,10)+48)
		DO WHILE ((ida==idb))
			num_val_en=num_val_en+1
			ida='b'//CHAR(MOD(num_val_en,10)+48)
		END DO
		b=a+first_step
		rb=r0+b*dir
		b_app=b
		CALL controlla_punto_e_parametro(r0,rb,b,dir)
		flag_border=.FALSE.
		IF (b_app/=b) flag_border=.TRUE.
		rb=r0+b*dir
		CALL calcola_energia(idb,rb,N,fb(1:2),accettabile)
		DO WHILE ((.NOT. accettabile).AND.(.NOT. flag_border))
			first_step=0.5d0*first_step
			b=a+first_step
			rb=r0+b*dir
			b_app=b
			CALL controlla_punto_e_parametro(r0,rb,b,dir)
			flag_border=.FALSE.
			IF (b_app/=b) flag_border=.TRUE.
			rb=r0+b*dir
			CALL calcola_energia(idb,rb,N,fb(1:2),accettabile)
		END DO
		DO WHILE ((DABS(fa(1)-fb(1))<0.25d0*liv_precisione).AND.(.NOT. flag_border))
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
				IF (flag_output) WRITE (7, *), 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
			END IF
			b=b+first_step
			rb=r0+b*dir
			b_app=b
			CALL controlla_punto_e_parametro(r0,rb,b,dir)
			flag_border=.FALSE.
			IF (b_app/=b) flag_border=.TRUE.
			rb=r0+b*dir
			CALL calcola_energia(idb,rb,N,fb(1:2),accettabile)
			IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 1 &
			  [ module_variational_opt.f90 > bracket_min_multidim ]'
		END DO
		flag_too_accurate=.FALSE.
		IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
		DO WHILE ((fa(2)+fb(2)>0.25d0*DABS(fa(1)-fb(1))) .AND. (.NOT. flag_too_accurate) .AND. (.NOT. flag_border))
			IF (fa(2)>fb(2)) THEN
				b_app=a
				CALL controlla_punto_e_parametro(r0,ra,a,dir)
				flag_border=.FALSE.
				IF (b_app/=a) flag_border=.TRUE.
				ra=r0+a*dir
				CALL aumenta_precisione_energia(ida,ra,N,fa(1:2))
			ELSE
				b_app=b
				CALL controlla_punto_e_parametro(r0,rb,b,dir)
				flag_border=.FALSE.
				IF (b_app/=b) flag_border=.TRUE.
				rb=r0+b*dir
				CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
			END IF
			IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
		END DO
		derivata=(fb(1)-fa(1))/(b-a)
		
		IF (derivata>0.d0) THEN
			c=b
			fc=fb
			idc=idb
			rc=rb
			b=a
			fb=fa
			idb=ida
			rb=ra
			a=b-first_step
			ra=r0+a*dir
			num_val_en=num_val_en+1
			ida='b'//CHAR(MOD(num_val_en,10)+48)
			DO WHILE ((ida==idb) .OR. (ida==idc))
				num_val_en=num_val_en+1
				ida='b'//CHAR(MOD(num_val_en,10)+48)
			END DO
			CALL controlla_punto_e_parametro(r0,ra,a,dir)
			ra=r0+a*dir
			CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
			DO WHILE ((.NOT. accettabile) .AND. (.NOT. flag_border))
				a=(b+a)*0.5d0
				ra=r0+a*dir
				b_app=a
				CALL controlla_punto_e_parametro(r0,ra,a,dir)
				flag_border=.FALSE.
				IF (b_app/=a) flag_border=.TRUE.
				ra=r0+a*dir
				CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
			END DO 
			DO WHILE ((DABS(fa(1)-fb(1))<0.5d0*liv_precisione) .AND. (.NOT. flag_border))
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
					IF (flag_output) WRITE (7, *), 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
				END IF
				!first_step=2.d0*first_step
				a=a-first_step
				ra=r0+a*dir
				b_app=a
				CALL controlla_punto_e_parametro(r0,ra,a,dir)
				flag_border=.FALSE.
				IF (b_app/=a) flag_border=.TRUE.
				ra=r0+a*dir
				CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
				IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 2 &
				  [ module_variational_opt.f90 > bracket_min_multidim ]'
			END DO
			flag_too_accurate=.FALSE.
			IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
			DO WHILE ((fa(2)+fb(2)>0.25d0*DABS(fa(1)-fb(1))) .AND. (.NOT. flag_too_accurate) .AND. (.NOT. flag_border))
				IF (fa(2)>fb(2)) THEN
					b_app=a
					CALL controlla_punto_e_parametro(r0,ra,a,dir)
					flag_border=.FALSE.
					IF (b_app/=a) flag_border=.TRUE.
					ra=r0+a*dir
					CALL aumenta_precisione_energia(ida,ra,N,fa(1:2))
				ELSE
					b_app=b
					CALL controlla_punto_e_parametro(r0,rb,b,dir)
					flag_border=.FALSE.
					IF (b_app/=b) flag_border=.TRUE.
					rb=r0+b*dir
					CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
				END IF
				IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
			END DO
			contatore=1
			DO WHILE ((fa(1)<fb(1)) .AND. (.NOT. flag_border))
				contatore=contatore+1
				c=b
				fc=fb
				idc=idb
				rc=rb
				b=a
				fb=fa
				idb=ida
				rb=ra
				a=a-(2.d0**contatore)*DABS(first_step)
				ra=r0+a*dir
				num_val_en=num_val_en+1
				ida='b'//CHAR(MOD(num_val_en,10)+48)
				DO WHILE ((ida==idb) .OR. (ida==idc))
					num_val_en=num_val_en+1
					ida='b'//CHAR(MOD(num_val_en,10)+48)
				END DO
				b_app=a
				CALL controlla_punto_e_parametro(r0,ra,a,dir)
				flag_border=.FALSE.
				IF (b_app/=a) flag_border=.TRUE.
				ra=r0+a*dir
				CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
				DO WHILE ((.NOT. accettabile) .AND. (.NOT. flag_border))
					a=(b+a)*0.5d0
					ra=r0+a*dir
					IF (a<b-(2.d0**(contatore-1))*DABS(first_step)) THEN
						b_app=a
						CALL controlla_punto_e_parametro(r0,ra,a,dir)
						flag_border=.FALSE.
						IF (b_app/=a) flag_border=.TRUE.
						ra=r0+a*dir
						CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
					END IF
				END DO
				flag_too_accurate=.FALSE.
				IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
				DO WHILE ((fa(2)+fb(2)>0.25d0*DABS(fa(1)-fb(1))) .AND. (.NOT. flag_too_accurate) .AND. (.NOT. flag_border))
					IF (fa(2)>fb(2)) THEN
						b_app=a
						CALL controlla_punto_e_parametro(r0,ra,a,dir)
						flag_border=.FALSE.
						IF (b_app/=a) flag_border=.TRUE.
						ra=r0+a*dir
						CALL aumenta_precisione_energia(ida,ra,N,fa(1:2))
					ELSE
						b_app=b
						CALL controlla_punto_e_parametro(r0,rb,b,dir)
						flag_border=.FALSE.
						IF (b_app/=b) flag_border=.TRUE.
						rb=r0+b*dir
						CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
					END IF
					IF (DABS(fa(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
				END DO
			END DO
		ELSE
			c=b+first_step
			rc=r0+c*dir
			num_val_en=num_val_en+1
			idc='b'//CHAR(MOD(num_val_en,10)+48)
			DO WHILE ((idc==ida) .OR. (idc==idb))
				num_val_en=num_val_en+1
				idc='b'//CHAR(MOD(num_val_en,10)+48)
			END DO
			b_app=c
			CALL controlla_punto_e_parametro(r0,rc,c,dir)
			flag_border=.FALSE.
			IF (b_app/=c) flag_border=.TRUE.
			rc=r0+c*dir
			CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
			DO WHILE ((.NOT. accettabile).AND.(.NOT. flag_border))
				c=(b+c)*0.5d0
				rc=r0+c*dir
				b_app=c
				CALL controlla_punto_e_parametro(r0,rc,c,dir)
				flag_border=.FALSE.
				IF (b_app/=c) flag_border=.TRUE.
				rc=r0+c*dir
				CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
			END DO
			DO WHILE ((DABS(fb(1)-fc(1))<0.25d0*liv_precisione).AND.(.NOT. flag_border))
				IF (mpi_myrank==0) THEN
					PRINT * , 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
					IF (flag_output) WRITE (7, *), 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
				END IF
				!first_step=2.d0*first_step
				c=c+first_step
				rc=r0+c*dir
				b_app=c
				CALL controlla_punto_e_parametro(r0,rc,c,dir)
				flag_border=.FALSE.
				IF (b_app/=c) flag_border=.TRUE.
				rc=r0+c*dir
				CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
				IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 3 &
				  [ module_variational_opt.f90 > bracket_min_multidim ]'
			END DO
			flag_too_accurate=.FALSE.
			IF (DABS(fb(1)-fc(1))<liv_precisione) flag_too_accurate=.TRUE.
			DO WHILE ((fb(2)+fc(2)>0.25d0*DABS(fb(1)-fc(1))) .AND. (.NOT. flag_too_accurate) .AND. (.NOT. flag_border))
				IF (fb(2)>fc(2)) THEN
					b_app=b
					CALL controlla_punto_e_parametro(r0,rb,b,dir)
					flag_border=.FALSE.
					IF (b_app/=b) flag_border=.TRUE.
					rb=r0+b*dir
					CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
				ELSE
					b_app=c
					CALL controlla_punto_e_parametro(r0,rc,c,dir)
					flag_border=.FALSE.
					IF (b_app/=c) flag_border=.TRUE.
					rc=r0+c*dir
					CALL aumenta_precisione_energia(idc,rc,N,fc(1:2))
				END IF
				IF (DABS(fb(1)-fc(1))<liv_precisione) flag_too_accurate=.TRUE.
			END DO
			contatore=1
			DO WHILE ((fc(1)<fb(1)).AND.(.NOT. flag_border))
				contatore=contatore+1
				a=b
				fa=fb
				ida=idb
				ra=rb
				b=c
				fb=fc
				idb=idc
				rb=rc
				c=c+(2.d0**contatore)*first_step
				rc=r0+c*dir
				num_val_en=num_val_en+1
				idc='b'//CHAR(MOD(num_val_en,10)+48)
				DO WHILE ((idc==ida) .OR. (idc==idb))
					num_val_en=num_val_en+1
					idc='b'//CHAR(MOD(num_val_en,10)+48)
				END DO
				b_app=c
				CALL controlla_punto_e_parametro(r0,rc,c,dir)
				flag_border=.FALSE.
				IF (b_app/=c) flag_border=.TRUE.
				rc=r0+c*dir
				CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
				IF (.NOT. accettabile) THEN
					c=(b+c)*0.5d0
					rc=r0+c*dir
					IF (c<b+(2.d0**(contatore-1))*DABS(first_step)) THEN
						b_app=c
						CALL controlla_punto_e_parametro(r0,rc,c,dir)
						flag_border=.FALSE.
						IF (b_app/=c) flag_border=.TRUE.
						rc=r0+c*dir
						CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
					END IF
				END IF
				flag_too_accurate=.FALSE.
				IF (DABS(fb(1)-fc(1))<liv_precisione) flag_too_accurate=.TRUE.
				DO WHILE ((fb(2)+fc(2)>0.25d0*DABS(fb(1)-fc(1))) .AND. (.NOT. flag_too_accurate) .AND. (.NOT. flag_border))
					IF (fb(2)>fc(2)) THEN
						b_app=b
						CALL controlla_punto_e_parametro(r0,rb,b,dir)
						flag_border=.FALSE.
						IF (b_app/=b) flag_border=.TRUE.
						rb=r0+b*dir
						CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
					ELSE
						b_app=c
						CALL controlla_punto_e_parametro(r0,rc,c,dir)
						flag_border=.FALSE.
						IF (b_app/=c) flag_border=.TRUE.
						rc=r0+c*dir
						CALL aumenta_precisione_energia(idc,rc,N,fc(1:2))
					END IF
					IF (DABS(fb(1)-fc(1))<liv_precisione) flag_too_accurate=.TRUE.
				END DO
			END DO
		END IF
		IF (mpi_myrank==0) THEN
			DO num_val_en = 1, 10, 1                    !elimino i dati ormai inutili
				id='b'//CHAR(MOD(num_val_en,10)+48)
				IF ((id/=ida) .AND. (id/=idb) .AND. (id/=idc)) CALL SYSTEM ("export ID_TO_ERASE=opt'//id//'; &
				  find estimatori/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \; &
				  find posizioni/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \")
			END DO
		END IF

		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Ho eseguito il bracket con ', INT(contatore), ' passi.',' I tre punti sono '
			PRINT * , '         punti:   ', a, b ,c
			PRINT * , '         energie: ', fa(1), fb(1), fc(1)
			IF (flag_output) THEN
				WRITE (7, *) , 'VAR_OPT: Ho eseguito il bracket con ', INT(contatore), ' passi.',' I tre punti sono '
				WRITE (7, *) , '         punti:   ', a, b ,c
				WRITE (7, *) , '         energie: ', fa(1), fb(1), fc(1)
			END IF
		END IF
		
	END SUBROUTINE bracket_min_multidim_v2

	SUBROUTINE parabgold_search_multidim(N,r0,dir,a,fa,ida,b,fb,idb,c,fc,idc)
		IMPLICIT NONE
		LOGICAL :: flag_too_accurate, accettabile
		CHARACTER(LEN=2) :: ida, idb, idc, idx, id
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8) :: a, b, c, dir(1:N), r0(1:N), rx(1:N), ra(1:N), rb(1:N), rc(1:N)
		INTEGER :: i, contatore, sx, dx, ierr
		REAL (KIND=8) :: fa(1:2), fb(1:2), fc(1:2), fx(1:2), delta, x, eps, p_old, f_old(1:2)
				
		contatore=0
		sx=0
		dx=0
		ra=r0+a*dir
		rb=r0+b*dir
		rc=r0+c*dir

		IF (fb(1)>fa(1)) THEN
			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
			STOP 'OPT: parabgold - a e b sbagliati'
		END IF
		IF (fb(1)>fc(1)) THEN
			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
			STOP 'OPT: parabgold - a e b sbagliati'
		END IF

		eps=MAX(fa(1),fc(1))-fb(1)

		DO WHILE (eps>liv_precisione)
			IF ((fc(1)<fb(1)) .OR. (fa(1)<fb(1))) THEN
				STOP 'OPT: parabgold - errore'
			END IF
			IF ((c<b) .OR. (a>b)) THEN
				STOP 'OPT: parabgold - errore'
			END IF
			contatore=contatore+1
			x=b-0.5d0*((b-a)*(b-a)*(fb(1)-fc(1))-(b-c)*(b-c)*(fb(1)-fa(1)))/((b-a)*(fb(1)-fc(1))-(b-c)*(fb(1)-fa(1)))
			IF (((x<b).AND.(sx>1)).OR.((x>b).AND.(dx>1))) THEN
				IF (b-a>c-b) THEN
					delta=(b-a)*0.38197d0*0.5d0
					x=b-delta
				ELSE
					delta=(c-b)*0.38197d0*0.5d0
					x=b+delta
				END IF
			END IF
			rx=r0+x*dir
			i=contatore
			idx='p'//CHAR(MOD(i,10)+48)
			DO WHILE ((idx==ida) .OR. (idx==idb) .OR. (idx==idc))
				i=i+1
				idx='p'//CHAR(MOD(i,10)+48)
			END DO
			CALL controlla_punto_e_parametro(r0,rx,x,dir)
			rx=r0+x*dir
			CALL calcola_energia(idx,rx,N,fx(1:2),accettabile)
			flag_too_accurate=.FALSE.
			IF (DABS(fx(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
			DO WHILE ((fx(2)+fb(2)>0.25d0*DABS(fx(1)-fb(1))) .AND. (.NOT. flag_too_accurate))
				IF (fx(2)>fb(2)) THEN
					CALL aumenta_precisione_energia(idx,rx,N,fx(1:2))
				ELSE
					CALL aumenta_precisione_energia(idb,rb,N,fb(1:2))
				END IF
				IF (DABS(fx(1)-fb(1))<liv_precisione) flag_too_accurate=.TRUE.
			END DO
			IF (x<b) THEN
				sx=sx+1
				dx=0
				IF (fx(1)>fb(1)) THEN
					a=x
					fa=fx
					ida=idx
					ra=rx
				ELSE
					c=b
					fc=fb
					idc=idb
					rc=rb
					b=x
					fb=fx
					idb=idx
					rb=rx
				END IF
			ELSE
				sx=0
				dx=dx+1
				IF (fx(1)>fb(1)) THEN
					c=x
					fc=fx
					idc=idx
					rc=rx
				ELSE
					a=b
					fa=fb
					ida=idb
					ra=rb
					b=x
					fb=fx
					idb=idx
					rb=rx
				END IF
			END IF
			eps=MAX(fa(1),fc(1))-fb(1)
			IF ((mpi_myrank==0) .AND. (contatore<10)) THEN
				PRINT * , 'VAR_OPT: Parabgold Search (a,b,c) - ', contatore
				PRINT * , '         ', a, b, c
				IF (flag_output) THEN
					WRITE (7, *), 'VAR_OPT: Parabgold Search (a,b,c) - ', contatore
					WRITE (7, *), '         ', a, b, c
				END IF
			END IF
		END DO
		IF (mpi_myrank==0) THEN
			DO i = 1, 10, 1                    !elimino i dati ormai inutili
				id='p'//CHAR(MOD(i,10)+48)
				IF (id/=idb) CALL SYSTEM ("export ID_TO_ERASE=opt'//id//'; &
				  find estimatori/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \; &
				  find posizioni/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \")
			END DO
			IF (idb/='b') CALL SYSTEM ("export ID_TO_ERASE=optb; &
			  find estimatori/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \; &
			  find posizioni/ -name '*'$ID_TO_ERASE'*' -exec rm -f '{}' \")
		END IF
		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Parabgold Search - sono stati necessari ', contatore, ' passi'
			IF (flag_output) WRITE (7, *), 'VAR_OPT: Parabgold Search - sono stati necessari ', contatore, ' passi'
		END IF

	END SUBROUTINE parabgold_search_multidim

	SUBROUTINE controlla_punto(nuovi_parametri)
		USE funzione_onda
		IMPLICIT NONE
		REAL (KIND=8), INTENT(INOUT) :: nuovi_parametri(1:num_par_var)
		INTEGER :: cont
		
		cont=1
		
		IF (SDe_kind/='no_') THEN
			SELECT CASE (SDe_kind)
			CASE ('pw_')
				
			END SELECT
		END IF
		IF (Jee_kind/='no_') THEN
			SELECT CASE (Jee_kind)
			CASE ('yuk')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Aep) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Fep) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
			CASE ('yup')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Aep) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Fep) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
			END SELECT
		END IF
		IF (Jep_kind/='no_') THEN
			SELECT CASE (Jep_kind)
			CASE ('yuk')
				IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Aep) THEN
					IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Fep) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
			CASE ('yup')
				IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Aep) THEN
					IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Fep) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
			END SELECT
		END IF
		IF (Kse_kind/='no_') THEN
			SELECT CASE (Kse_kind)
			CASE ('gss')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			CASE ('gsc')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			CASE ('gsp')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			CASE ('gsd')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			CASE ('gdc')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			CASE ('gdp')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			CASE ('atm')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			END SELECT
		END IF
		IF (Jse_kind/='no_') THEN
			SELECT CASE (Jse_kind)
			CASE ('pot')
				IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
				IF (nuovi_parametri(cont+1)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+2
			CASE ('bou')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				IF (nuovi_parametri(cont+1)<0.d0) nuovi_parametri(cont+1)=0.d0
				cont=cont+2
			CASE ('ppb')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				IF (nuovi_parametri(cont+1)<0.d0) nuovi_parametri(cont+1)=0.d0
				IF (nuovi_parametri(cont+2)<0.d0) nuovi_parametri(cont+2)=0.d0
				cont=cont+3
			CASE ('yuk')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Asese) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Fsese) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
			CASE ('yup')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Asese) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Fsese) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
			END SELECT
		END IF
		IF (Jsesp_kind/='no_') THEN
			SELECT CASE (Jsesp_kind)
			CASE ('pot')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			CASE ('yuk')
				IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Asese) THEN
					IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Fsese) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
			CASE ('yup')
				IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Asese) THEN
					IF (nuovi_parametri(cont)>0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
				IF (split_Fsese) THEN
					IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
					cont=cont+1
				END IF
			CASE ('gss')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			CASE ('gsd')
				IF (nuovi_parametri(cont)<0.d0) nuovi_parametri(cont)=0.d0
				cont=cont+1
			END SELECT
		END IF
		
	END SUBROUTINE controlla_punto

	SUBROUTINE controlla_punto_e_parametro(P_old,P_new,a,dir)
		USE funzione_onda
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: P_old(1:num_par_var), P_new(1:num_par_var)
		REAL (KIND=8), INTENT(IN) :: dir(1:num_par_var)
		REAL (KIND=8), INTENT(INOUT) :: a
		INTEGER :: cont
		
		cont=1
		
		IF (SDe_kind/='no_') THEN
			SELECT CASE (SDe_kind)
			CASE ('pw_')
				
			END SELECT
		END IF
		IF (Jee_kind/='no_') THEN
			SELECT CASE (Jee_kind)
			CASE ('yuk')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Aee_yuk fuori range'
				cont=cont+1
				IF (split_Aee) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Aee_ud_yuk fuori range'
					cont=cont+1
				END IF
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fee_yuk fuori range'
				cont=cont+1
				IF (split_Fee) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fee_ud_yuk fuori range'
					cont=cont+1
				END IF
			CASE ('yup')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Aee_yuk fuori range'
				cont=cont+1
				IF (split_Aee) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Aee_ud_yuk fuori range'
					cont=cont+1
				END IF
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fee_yuk fuori range'
				cont=cont+1
				IF (split_Fee) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fee_ud_yuk fuori range'
					cont=cont+1
				END IF
			END SELECT
		END IF
		IF (Jep_kind/='no_') THEN
			SELECT CASE (Jep_kind)
			CASE ('yuk')
				IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Aep_yuk fuori range'
				cont=cont+1
				IF (split_Aep) THEN
					IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Aep_ud_yuk fuori range'
					cont=cont+1
				END IF
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fep_yuk fuori range'
				cont=cont+1
				IF (split_Fep) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fep_ud_yuk fuori range'
					cont=cont+1
				END IF
			CASE ('yup')
				IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Aep_yuk fuori range'
				cont=cont+1
				IF (split_Aep) THEN
					IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Aep_ud_yuk fuori range'
					cont=cont+1
				END IF
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fep_yuk fuori range'
				cont=cont+1
				IF (split_Fep) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fep_ud_yuk fuori range'
					cont=cont+1
				END IF
			END SELECT
		END IF
		IF (Kse_kind/='no_') THEN
			SELECT CASE (Kse_kind)
			CASE ('gss')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro C_kern_e fuori range'
				cont=cont+1
			CASE ('gsc')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro C_kern_e fuori range'
				cont=cont+1
			CASE ('gsd')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro C_kern_e fuori range'
				cont=cont+1
			CASE ('gdc')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro C_kern_e fuori range'
				cont=cont+1
			CASE ('gdp')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro C_kern_e fuori range'
				cont=cont+1
			CASE ('atm')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro C_kern_e fuori range'
				cont=cont+1
			END SELECT
		END IF
		IF (Jse_kind/='no_') THEN
			SELECT CASE (Jse_kind)
			CASE ('pot')
				IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro A_POT_se fuori range'
				cont=cont+1
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro D_POT_se fuori range'
				cont=cont+1
			!CASE ('bou')
			!	IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
			!	cont=cont+1
			!	IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
			!	cont=cont+1
			!CASE ('ppb')
			!	IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
			!	cont=cont+1
			!	IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
			!	cont=cont+1
			!	IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
			!	cont=cont+1
			CASE ('yuk')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Asese_yuk fuori range'
				cont=cont+1
				IF (split_Asese) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Asese_ud_yuk fuori range'
					cont=cont+1
				END IF
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fsese_yuk fuori range'
				cont=cont+1
				IF (split_Fsese) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fsese_ud_yuk fuori range'
					cont=cont+1
				END IF
			CASE ('yup')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Asese_yuk fuori range'
				cont=cont+1
				IF (split_Asese) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Asese_ud_yuk fuori range'
					cont=cont+1
				END IF
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fsese_yuk fuori range'
				cont=cont+1
				IF (split_Fsese) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fsese_ud_yuk fuori range'
					cont=cont+1
				END IF
			END SELECT
		END IF
		IF (Jsesp_kind/='no_') THEN
			SELECT CASE (Jsesp_kind)
			!CASE ('pot')
			!	IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
			!	cont=cont+1
			CASE ('yuk')
				IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Asesp_yuk fuori range'
				cont=cont+1
				IF (split_Asesp) THEN
					IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Asesp_ud_yuk fuori range'
					cont=cont+1
				END IF
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fsesp_yuk fuori range'
				cont=cont+1
				IF (split_Fsesp) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fsesp_ud_yuk fuori range'
					cont=cont+1
				END IF
			CASE ('yup')
				IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Asesp_yuk fuori range'
				cont=cont+1
				IF (split_Asesp) THEN
					IF (P_new(cont)>0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)>0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Asesp_ud_yuk fuori range'
					cont=cont+1
				END IF
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fsesp_yuk fuori range'
				cont=cont+1
				IF (split_Fsesp) THEN
					IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
					IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Fsesp_ud_yuk fuori range'
					cont=cont+1
				END IF
			CASE ('gss')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Gsesp fuori range'
				cont=cont+1
			CASE ('gsd')
				IF (P_new(cont)<0.d0) a=MIN(DABS(a),DABS(-P_old(cont)/dir(cont)))*DSIGN(1.d0,a)
				IF ((P_new(cont)<0.d0).AND.(mpi_myrank==0)) PRINT * , 'VAR_OPT: parametro Gsesp fuori range'
				cont=cont+1
			END SELECT
		END IF
		
	END SUBROUTINE controlla_punto_e_parametro

	SUBROUTINE estrai_s_kl_e_f_k(E)
		USE VMC
		USE variational_calculations
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: E(1:2)
		CHARACTER(LEN=4) :: codice_numerico1, codice_numerico2
		INTEGER :: i, j
		INTEGER (KIND=8) :: i_mc
		REAL (KIND=8) :: Oi(1:num_par_var), OiOj(1:num_par_var,1:num_par_var), OiE(1:num_par_var)
		REAL (KIND=8) :: media, errore, estimatore(1:2), dummy1(1:N_mc)
		
		s_kl=0.d0
		f_k=0.d0
		
		IF (quick_error==0) THEN
			IF (flag_disk) THEN
				!IF (SDse_kind=='no_') THEN
				!	IF (mpi_myrank==0) THEN
				!		DO j = 1, num_par_var, 1
				!			WRITE (codice_numerico1, '(I4.4)'), j
				!			CALL calcola_estimatori('estimatori/gradiente/O_'//codice_numerico1,media,errore)
				!			Oi(j)=media
				!			CALL calcola_estimatori('estimatori/gradiente/EO_'//codice_numerico1,media,errore)
				!			OiE(j)=media
				!			DO i = j, num_par_var, 1
				!				WRITE (codice_numerico2, '(I4.4)'), i
				!				CALL calcola_estimatori('estimatori/gradiente/O_'//codice_numerico1//'-'//codice_numerico2,media,errore)
				!				OiOj(i,j)=media
				!				OiOj(j,i)=OiOj(i,j)
				!			END DO
				!		END DO
				!	END IF
				!ELSE
					IF (mpi_myrank==0) THEN
						DO j = 1, num_par_var, 1
							WRITE (codice_numerico1, '(I4.4)'), j
							CALL calcola_estimatori('estimatori/gradiente/O_'//codice_numerico1,media,errore, &
							  'estimatori/w'//codice_simulazione)
							Oi(j)=media
							CALL calcola_estimatori('estimatori/gradiente/EO_'//codice_numerico1,media,errore, &
							  'estimatori/w'//codice_simulazione)
							OiE(j)=media
							DO i = j, num_par_var, 1
								WRITE (codice_numerico2, '(I4.4)'), i
								CALL calcola_estimatori('estimatori/gradiente/O_'//codice_numerico1//'-'//codice_numerico2,media,errore, &
								  'estimatori/w'//codice_simulazione)
								OiOj(i,j)=media
							END DO
						END DO
					END IF
				!END IF
			ELSE
				!IF (SDse_kind=='no_') THEN
				!	DO j = 1, num_par_var, 1
				!		CALL calcola_estimatore_da_RAM(estimatore,O(j,1:N_mc),N_mc)
				!		Oi(j)=estimatore(1)
				!		DO i_mc = 1, N_mc, 1
				!			dummy1(i_mc)=E_tot(i_mc)*O(j,i_mc)
				!		END DO
				!		CALL calcola_estimatore_da_RAM(estimatore,dummy1(1:N_mc),N_mc)
				!		OiE(j)=estimatore(1)
				!		DO i = j, num_par_var, 1
				!			DO i_mc = 1, N_mc, 1
				!				dummy1(i_mc)=O(i,i_mc)*O(j,i_mc)
				!			END DO
				!			CALL calcola_estimatore_da_RAM(estimatore,dummy1(1:N_mc),N_mc)
				!			OiOj(i,j)=estimatore(1)
				!			OiOj(j,i)=estimatore(1)
				!		END DO
				!	END DO
				!ELSE
					DO j = 1, num_par_var, 1
						CALL calcola_estimatore_da_RAM(estimatore,O(j,1:N_mc),N_mc,w)
						Oi(j)=estimatore(1)
						DO i_mc = 1, N_mc, 1
							dummy1(i_mc)=E_tot(i_mc)*O(j,i_mc)
						END DO
						CALL calcola_estimatore_da_RAM(estimatore,dummy1(1:N_mc),N_mc,w)
						OiE(j)=estimatore(1)
						DO i = j, num_par_var, 1
							DO i_mc = 1, N_mc, 1
								dummy1(i_mc)=O(i,i_mc)*O(j,i_mc)
							END DO
							CALL calcola_estimatore_da_RAM(estimatore,dummy1(1:N_mc),N_mc,w)
							OiOj(i,j)=estimatore(1)
							OiOj(j,i)=estimatore(1)
						END DO
					END DO
				!END IF
			END IF

			DO j = 1, num_par_var, 1
				DO i = j, num_par_var, 1
					s_kl(i,j)=OiOj(i,j)-Oi(i)*Oi(j)
					s_kl(j,i)=s_kl(i,j)
				END DO
			END DO
			DO i = 1, num_par_var, 1
				f_k(i)=Oi(i)*E(1)-OiE(i)
			END DO
			
		ELSE IF (quick_error>3) THEN
			IF (flag_disk) THEN
				!IF (SDse_kind=='no_') THEN
				!	IF (mpi_myrank==0) THEN
				!		DO j = 1, num_par_var, 1
				!			WRITE (codice_numerico1, '(I4.4)'), j
				!			CALL calcola_estimatori_quick(quick_error,'estimatori/gradiente/O_'//codice_numerico1,media,errore)
				!			Oi(j)=media
				!			CALL calcola_estimatori_quick(quick_error,'estimatori/gradiente/EO_'//codice_numerico1,media,errore)
				!			OiE(j)=media
				!			DO i = j, num_par_var, 1
				!				WRITE (codice_numerico2, '(I4.4)'), i
				!				CALL calcola_estimatori_quick(quick_error,'estimatori/gradiente/O_'//codice_numerico1//'-'//codice_numerico2,media,errore)
				!				OiOj(i,j)=media
				!				OiOj(j,i)=OiOj(i,j)
				!			END DO
				!		END DO
				!	END IF
				!ELSE
					IF (mpi_myrank==0) THEN
						DO j = 1, num_par_var, 1
							WRITE (codice_numerico1, '(I4.4)'), j
							CALL calcola_estimatori_quick(quick_error,'estimatori/gradiente/O_'//codice_numerico1,media,errore, &
							  'estimatori/w'//codice_simulazione)
							Oi(j)=media
							CALL calcola_estimatori_quick(quick_error,'estimatori/gradiente/EO_'//codice_numerico1,media,errore, &
							  'estimatori/w'//codice_simulazione)
							OiE(j)=media
							DO i = j, num_par_var, 1
								WRITE (codice_numerico2, '(I4.4)'), i
								CALL calcola_estimatori_quick(quick_error,'estimatori/gradiente/O_'//codice_numerico1//'-'//codice_numerico2,media,errore, &
								  'estimatori/w'//codice_simulazione)
								OiOj(i,j)=media
							END DO
						END DO
					END IF
				!END IF
			ELSE
				!IF (SDse_kind=='no_') THEN
				!	DO j = 1, num_par_var, 1
				!		CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,O(j,1:N_mc),N_mc)
				!		Oi(j)=estimatore(1)
				!		DO i_mc = 1, N_mc, 1
				!			dummy1(i_mc)=E_tot(i_mc)*O(j,i_mc)
				!		END DO
				!		CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,dummy1(1:N_mc),N_mc)
				!		OiE(j)=estimatore(1)
				!		DO i = j, num_par_var, 1
				!			DO i_mc = 1, N_mc, 1
				!				dummy1(i_mc)=O(i,i_mc)*O(j,i_mc)
				!			END DO
				!			CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,dummy1(1:N_mc),N_mc)
				!			OiOj(i,j)=estimatore(1)
				!			OiOj(j,i)=estimatore(1)
				!		END DO
				!	END DO
				!ELSE
					DO j = 1, num_par_var, 1
						CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,O(j,1:N_mc),N_mc,w)
						Oi(j)=estimatore(1)
						DO i_mc = 1, N_mc, 1
							dummy1(i_mc)=E_tot(i_mc)*O(j,i_mc)
						END DO
						CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,dummy1(1:N_mc),N_mc,w)
						OiE(j)=estimatore(1)
						DO i = j, num_par_var, 1
							DO i_mc = 1, N_mc, 1
								dummy1(i_mc)=O(i,i_mc)*O(j,i_mc)
							END DO
							CALL calcola_estimatore_da_RAM_quick(quick_error,estimatore,dummy1(1:N_mc),N_mc,w)
							OiOj(i,j)=estimatore(1)
							OiOj(j,i)=estimatore(1)
						END DO
					END DO
				!END IF
			END IF

			DO j = 1, num_par_var, 1
				DO i = j, num_par_var, 1
					s_kl(i,j)=OiOj(i,j)-Oi(i)*Oi(j)
					s_kl(j,i)=s_kl(i,j)
				END DO
			END DO
			DO i = 1, num_par_var, 1
				f_k(i)=Oi(i)*E(1)-OiE(i)
			END DO
			
		ELSE
			STOP 'Il parametro quick_error non é stato scelto opportunamente'
		END IF
				
	END SUBROUTINE estrai_s_kl_e_f_k

	SUBROUTINE chiudi_variational_opt()
		IMPLICIT NONE
		DEALLOCATE(parametri_var)
	END SUBROUTINE chiudi_variational_opt
	
END MODULE variational_opt
