MODULE variational_opt
	USE dati_mc
	IMPLICIT NONE
	LOGICAL, PARAMETER :: verbose_mode=.TRUE.
   REAL(KIND=8), PRIVATE, SAVE :: SVD_MIN=1.d-15
	INTEGER, SAVE :: num_par_var, num_coord_Rp, num_coord_L, num_par_var_eff
	REAL (KIND=8) :: E_min(1:2), liv_precisione
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: parametri_var(:)
   REAL(KIND=8), PROTECTED, SAVE :: H
   REAL(KIND=8), ALLOCATABLE, PROTECTED, SAVE :: Oi(:), HOi(:), OiOj(:,:)
   REAL(KIND=8), ALLOCATABLE, PROTECTED, SAVE :: Hi(:), OiHj(:,:)
   REAL(KIND=8), ALLOCATABLE, PROTECTED, SAVE :: OiHOj(:,:)
   !REAL(KIND=8), PROTECTED, SAVE :: sigmaH
   !REAL(KIND=8), ALLOCATABLE, PROTECTED, SAVE :: sigmaHOi(:)
   REAL(KIND=8), ALLOCATABLE, PROTECTED, SAVE :: Oi_eff(:), HOi_eff(:), OiOj_eff(:,:)
   REAL(KIND=8), ALLOCATABLE, PROTECTED, SAVE :: OiHOj_eff(:,:)
   !REAL(KIND=8), ALLOCATABLE, PROTECTED, SAVE :: sigmaHOi_eff(:)
	!!!REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: c_knm(:,:,:)
   REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: s_kn(:,:), f_k(:)  !parametri per SR
   REAL(KIND=8), PROTECTED, SAVE :: f_SR_beta    !dynamical factor for SR_beta
   LOGICAL, ALLOCATABLE :: used_par(:), M_used_par(:,:)	

	CONTAINS

	SUBROUTINE inizializza_variational_opt()
		USE dati_fisici
		!USE funzione_onda
		USE walkers
      		USE variational_calculations
		IMPLICIT NONE
		INTEGER :: cont, i, j, ih
		
		cont=1
		num_par_var=0
		num_coord_Rp=0
      num_coord_L=0
		
		CALL inizializza_dati_fisici()
		CALL inizializza_dati_mc()
		CALL inizializza_walkers('iniz_opt')
		CALL inizializza_dati_funzione_onda()
		
      IF (opt_L) THEN
         num_par_var=num_par_var+3
         num_coord_L=3
      END IF
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
         CASE ('spl','spp')
            IF (opt_A_Jee.OR.opt_F_Jee) THEN
               IF (split_Aee.OR.split_Fee) THEN
                  num_par_var=num_par_var+(Jsplee%Nknots+1)*(Jsplee%m+1)
                  num_par_var=num_par_var+(Jsplee_ud%Nknots+1)*(Jsplee_ud%m+1)
               ELSE
                  num_par_var=num_par_var+(Jsplee%Nknots+1)*(Jsplee%m+1)
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
         CASE ('spl','spp')
            IF (opt_A_Jep.OR.opt_F_Jep) THEN
               IF (split_Aep.OR.split_Fep) THEN
                  num_par_var=num_par_var+(Jsplep%Nknots+1)*(Jsplep%m+1)
                  num_par_var=num_par_var+(Jsplep_ud%Nknots+1)*(Jsplep_ud%m+1)
               ELSE
                  num_par_var=num_par_var+(Jsplep%Nknots+1)*(Jsplep%m+1)
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
				CASE ('atm','atp')
					num_par_var=num_par_var+1
				CASE ('bat','bap')
					num_par_var=num_par_var+1
				CASE ('hl_')
					num_par_var=num_par_var+1
				CASE ('1sb')
               IF (opt_orbital) num_par_var=num_par_var+1
               IF (opt_dynamic_backflow) num_par_var=num_par_var+2 
				CASE ('spb')
               IF (opt_orbital) num_par_var=num_par_var+(Bsplep%m+1)*(Bsplep%nknots+1)+2
               IF (opt_dynamic_backflow) num_par_var=num_par_var+2 
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
      ELSE
         num_coord_Rp=0
		END IF
		
		!!!!!!!
		ALLOCATE(parametri_var(1:num_par_var))
      ALLOCATE(used_par(1:num_par_var),M_used_par(1:num_par_var,1:num_par_var))
      ALLOCATE(Oi(1:num_par_var),HOi(1:num_par_var),OiOj(1:num_par_var,1:num_par_var))
      ALLOCATE(OiHOj(1:num_par_var,1:num_par_var))
      ALLOCATE(Hi(1:num_coord_L+num_coord_Rp),OiHj(1:num_coord_L+num_coord_Rp,1:num_coord_L+num_coord_Rp))
		!!!!!!!
		

      IF ( opt_L ) THEN
         parametri_var(cont:cont+2)=L(1:3)
         cont=cont+3
      END IF
						
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
         CASE ('spl','spp')
            IF (opt_A_Jee.OR.opt_F_Jee) THEN
               IF (split_Aee.OR.split_Fee) THEN
                  parametri_var(cont:cont+(Jsplee%Nknots+1)*(Jsplee%m+1)-1)=&
                     RESHAPE(Jsplee%t(0:Jsplee%m,0:Jsplee%Nknots),(/(Jsplee%Nknots+1)*(Jsplee%m+1)/))
                  cont=cont+(Jsplee%Nknots+1)*(Jsplee%m+1)
                  parametri_var(cont:cont+(Jsplee_ud%nknots+1)*(Jsplee_ud%m+1)-1)=&
                     RESHAPE(Jsplee_ud%t(0:Jsplee_ud%m,0:Jsplee_ud%Nknots),(/(Jsplee%Nknots+1)*(Jsplee%m+1)/))
                  cont=cont+(Jsplee_ud%Nknots+1)*(Jsplee_ud%m+1)
               ELSE
                  parametri_var(cont:cont+(Jsplee%Nknots+1)*(Jsplee%m+1)-1)=&
                     RESHAPE(Jsplee%t(0:Jsplee%m,0:Jsplee%Nknots),(/(Jsplee%Nknots+1)*(Jsplee%m+1)/))
                  cont=cont+(Jsplee%Nknots+1)*(Jsplee%m+1)
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
         CASE ('spl','spp')
            IF (opt_A_Jep.OR.opt_F_Jep) THEN
               IF (split_Aep.OR.split_Fep) THEN
                  parametri_var(cont:cont+(Jsplep%Nknots+1)*(Jsplep%m+1)-1)=&
                     RESHAPE(Jsplep%t(0:Jsplep%m,0:Jsplep%Nknots),(/(Jsplep%Nknots+1)*(Jsplep%m+1)/))
                  cont=cont+(Jsplep%Nknots+1)*(Jsplep%m+1)
                  parametri_var(cont:cont+(Jsplep_ud%nknots+1)*(Jsplep_ud%m+1)-1)=&
                     RESHAPE(Jsplep_ud%t(0:Jsplep_ud%m,0:Jsplep_ud%Nknots),(/(Jsplep%Nknots+1)*(Jsplep%m+1)/))
                  cont=cont+(Jsplep_ud%Nknots+1)*(Jsplep_ud%m+1)
               ELSE
                  parametri_var(cont:cont+(Jsplep%Nknots+1)*(Jsplep%m+1)-1)=&
                     RESHAPE(Jsplep%t(0:Jsplep%m,0:Jsplep%Nknots),(/(Jsplep%Nknots+1)*(Jsplep%m+1)/))
                  cont=cont+(Jsplep%Nknots+1)*(Jsplep%m+1)
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
				CASE ('atm','atp')
					parametri_var(cont)=C_atm
					cont=cont+1
				CASE ('bat','bap')
					parametri_var(cont)=C_atm
					cont=cont+1
				CASE ('hl_')
					parametri_var(cont)=C_atm
					cont=cont+1
				CASE ('1sb')
               IF (opt_orbital) THEN
					   parametri_var(cont)=C_atm
					   cont=cont+1
               END IF
               IF (opt_dynamic_backflow) THEN
					   parametri_var(cont)=A_POT_se
					   cont=cont+1
					   parametri_var(cont)=D_POT_se
					   cont=cont+1
               END IF
				CASE ('spb')
               IF (opt_orbital) THEN
					   parametri_var(cont:cont+(Bsplep%m+1)*(Bsplep%nknots+1)-1)=&
                     RESHAPE(Bsplep%t(0:Bsplep%m,0:Bsplep%nknots),(/(Bsplep%m+1)*(Bsplep%nknots+1)/))
					   cont=cont+(Bsplep%m+1)*(Bsplep%nknots+1)
               END IF
               IF (opt_dynamic_backflow) THEN
					   parametri_var(cont)=A_POT_se
					   cont=cont+1
					   parametri_var(cont)=D_POT_se
					   cont=cont+1
               END IF
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
			DO i = 1, N_part, 1
				parametri_var(cont:cont+2)=rp_old(1:3,i)
				cont=cont+3
			END DO
		END IF

		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: il numero di parametri da ottimizzare é ', num_par_var
			IF (flag_output) WRITE (7, *) 'VAR_OPT: il numero di parametri da ottimizzare é ', num_par_var
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
		CASE('stoc_ns')   !SR senza fermarsi
			CALL stochastic_reconfiguration(parametri_var,num_par_var,.FALSE.,.FALSE.)
		CASE('stoc_av')   !SR senza fermarsi accumulando AV_wf.d
			CALL stochastic_reconfiguration(parametri_var,num_par_var,.FALSE.,.TRUE.)
		END SELECT
		IF (mpi_myrank==0) THEN
			WRITE (stringa, '(I4.4)') num_par_var
			PRINT * , '    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # '
			PRINT * , 'VAR_OPT: Ottimizzazione terminata.'
			!PRINT * , 'Optimal parameters: '
			!PRINT '('//stringa//'(F9.3,1X))' , parametri_var
			PRINT '(11x,A9,F9.6,A6,F9.6)' , 'Energia: ', E_min(1), '  +-  ', E_min(2)
			PRINT * , '    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # '
			IF (flag_output) THEN
				WRITE (7, *) '    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # '
				WRITE (7, *) 'VAR_OPT: Ottimizzazione terminata.'
				!WRITE (7, *) 'Optimal parameters: '
				!WRITE (7, '('//stringa//'(F9.3,1X))') parametri_var
				WRITE (7, '(11X,A9,F9.6,A6,F9.6)') 'Energia: ', E_min(1), '  +-  ', E_min(2)
				WRITE (7, *) '    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # '
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
			IF (cont_N_mc_increase>5) STOP 'N_mc é diventato troppo grande[ module_variational_opt.f90 > aumenta_precisione_energie ]'
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
			IF (cont_N_mc_increase>5) STOP 'N_mc é diventato troppo grand [ module_variational_opt.f90 > aumenta_precisione_energie ]'
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
		CHARACTER(LEN=*) :: id
		INTEGER :: num, i
		REAL (KIND=8) :: par(1:num), energia(1:2), par_pt(1:num,1)
      REAL(KIND=8) :: foo1, foo2
		
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
		CALL inizializza_VMC('opt'//id,CONTINUA=.FALSE.)

		CALL setta_parametri(par,num)
      IF (mpi_myrank==0) THEN
         PRINT *, "Cambiato L = ", L, "   [bohr]"
         IF (flag_output) WRITE(UNIT=7, FMT=*) "Cambiato L = ", L, "   [bohr]"
      END IF
		CALL prima_valutazione_funzione_onda()
		CALL valuta_step_mc(accettabile)
		IF (accettabile) THEN
         CALL sampling(.FALSE.)
         CALL inizializza_variational_calculations(num,par_pt,1)
         CALL sampling(.TRUE.)
         CALL trascrivi_dati()
         CALL restituisci_energia(energia)
         CALL estrai_medie_SR()
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
				IF (flag_output) WRITE (7, *) 'VAR_OPT: ruoto la base per la ricerca del minimo'
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
					IF (flag_output) WRITE (7, *) 'VAR_OPT: Ruotando la base ho trovato un nuovo minimo'
				END IF
			ELSE
				contatore=contatore+1
			END IF
			IF (contatore>N_max_rot) flag_loop=.FALSE.
			IF (mpi_myrank==0) THEN
				PRINT * , 'VAR_OPT: Ruotare la base é stato inutile'
				IF (flag_output) WRITE (7, *) 'VAR_OPT: Ruotare la base é stato inutile'
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
						WRITE (7, *) 'VAR_OPT: cerco il minimo lungo l asse ', i
						WRITE (7, "(8X,20(2X,F5.2))") direzione
					END IF
				END IF
				CALL bracket_min_multidim_v2(N,p0,direzione,first_step,a,fa,ida,b,fb,idb,c,fc,idc)
				CALL parabgold_search_multidim(N,p0,direzione,a,fa,ida,b,fb,idb,c,fc,idc)
				CALL SYSTEM ('export ID_TO_ERASE=opt00; find estimatori/ -name *$ID_TO_ERASE* -exec rm -f {} \\; find posizioni/ -name *$ID_TO_ERASE* -exec rm -f {} \\')
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
						WRITE (7, *) 'VAR_OPT: Seguendo l asse ', i, ' ho trovato come nuovo minimo ', &
						  fb(1), ' +- ', fb(2)
						WRITE (7, *) '         usando il fattore b=', b
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
					WRITE (7, *) 'VAR_OPT: Alla fine di una routine su tutti gli assi ho trovato come nuovo minimo ', &
					  fnew(1), ' +- ', fnew(2)
					WRITE (7, *) '                                                                     partendo da ', &
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
		IF (.NOT. accettabile) STOP 'I parametri variazionali iniziali non sono accettabili [ module_variational_opt.f90 > conj_grad ]'
		DO i = 1, N, 1
			grad(i)=f0(1,i)/step
		END DO
		
		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Ho calcolato il primo gradiente. La direzione da seguire é '
			IF (flag_output) WRITE (7, *) 'VAR_OPT: Ho calcolato il primo gradiente. La direzione da seguire é '
			DO i = 1, N, 1
				PRINT * , i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
				IF (flag_output) WRITE (7, *) i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
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
				IF (flag_output) WRITE (7, *) 'VAR_OPT: Ho migliorato la precisione del gradiente. La direzione da seguire é '
				DO i = 1, N, 1
					PRINT * , i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
					WRITE (7, *) i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
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
			IF (flag_file) CALL SYSTEM ('export ID_TO_ERASE=opt00; find estimatori/ -name *$ID_TO_ERASE* -exec rm -f {} \\; find posizioni/ -name *$ID_TO_ERASE* -exec rm -f {} \\')
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
					IF (flag_output) WRITE (7, *) 'VAR_OPT: Ho migliorato la precisione del gradiente. La direzione da seguire é ' 
					DO i = 1, N, 1
						PRINT * , i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
						IF (flag_output) WRITE (7, *) i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
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
				IF (flag_output) WRITE (7, *) 'VAR_OPT: Ho calcolato il gradiente. La direzione da seguire é '
				DO i = 1, N, 1
					PRINT * , i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
					IF (flag_output) WRITE (7, *) i, ' - ', f0(1,i)/step, ' +- ', f0(2,i)/step
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
					WRITE (7, *) 'VAR_OPT: ho calcolato la nuova direzione'
					WRITE (7, *) direzione
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
						WRITE (7, *) 'VAR_OPT: il nuovo punto per il minimo provvisorio é '
						WRITE (7, *) p0
					END IF
				END IF
				contatore=contatore+1
			END IF
		END DO
		minimo=fb(1:2)
		
		END IF cerca_minimo
		
		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Conjugate Gradient - sono stati necessari ', contatore, ' passi'
			IF (flag_output) WRITE (7, *) 'VAR_OPT: Conjugate Gradient - sono stati necessari ', contatore, ' passi'
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
      INTEGER :: N_eff
		LOGICAL :: flag_loop, loop_stop(1:6), accettabile, flag_file, flag_trovato_minimo, flag_freeze
		CHARACTER(LEN=6) :: id
		CHARACTER(LEN=4) :: stringa, istring
		INTEGER :: i, j, k, info, contatore, IO, i_orbit, AV_cont, WO_MIN_cont
      INTEGER :: lwork
		REAL (KIND=8) :: lambda, lambda2, flambda, eps, boundd, foo, dummy(1:4), dummy2(0:N)
		REAL (KIND=8) :: lambda_Rp, lambda2_Rp, flambda_Rp
		REAL (KIND=8) :: p0(1:N), dp_iniz(1:N), AV_P(1:N)
		REAL (KIND=8) :: dp(1:N), dp0(1:N), dp_old(1:N), Is_kn(1:N,1:N), pvt(1:N), cont_step
      REAL(KIND=8) :: var_change, var_change_min, Hmin, Hnext, dpnext(1:N), opt_f_SR_beta, static_f_SR_beta
      REAL(KIND=8) :: SRopt_count
      REAL(KIND=8) :: opt_SVD_MIN, static_SVD_MIN
      REAL(KIND=8) :: lambda3, opt_lambda3, static_lambda3
      REAL(KIND=8) :: SRtarget, SRtargetmin
      REAL(KIND=8) :: array_SRtargetmin(0:mpi_nprocs-1)
      REAL(KIND=8), ALLOCATABLE :: work(:)
		REAL (KIND=8) :: D_tr(1:N), E_tr(1:N-1), TAU_tr(1:N-1)
		REAL (KIND=8) :: energia(1:2), minE_AV(1:2), vec_app(1:N), vec_app_old(1:N)
		REAL (KIND=8) :: old_step(1:5), old_p0(1:N,1:5), old_E(1:2,1:5), derivE(1:4)
		REAL (KIND=8) :: old_step_Rp(1:5)
		REAL (KIND=8) :: dist_poss, dist_perc, dist_poss_5, dist_perc_5, av_errore, derE, nextdE
		REAL (KIND=8) :: old_Emin(1:2), p0_minimo(1:N), Pdist_min(1:5)
		REAL (KIND=8), ALLOCATABLE :: gRpuu(:), gRpuu_new(:), gRpud(:), gRpud_new(:), x_gr(:)
      REAL(KIND=8) :: Usvd(1:N,1:N), VTsvd(1:N,1:N), Ssvd(1:N)
      REAL(KIND=8) :: a, b, c, fa(1:2), fb(1:2), fc(1:2)

      lambda=lambda_init
      lambda2=1.d0
      lambda3=1.d0
      lambda_Rp=lambda_Rp_init
      lambda2_Rp=1.d0
      f_SR_beta=1.d0
		i_orbit=1
		AV_cont=0
      WO_MIN_cont=0
		AV_P=0.d0
		contatore=0
      WRITE (istring, '(I4.4)') contatore
      id='SR'//istring
      minE_AV=(/ HUGE(0.d0), 0.d0 /)
		
      lwork=10*N
      ALLOCATE(work(1:lwork))

		INQUIRE(FILE='ottimizzazione/SR_energies.dat',EXIST=flag_file)
		IF ((flag_file).AND.(mpi_myrank==0)) THEN
			OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='OLD')
			IO=0
			DO WHILE (IO>=0)
				READ (8, *, IOSTAT = IO) dummy(1:3)
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
            CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
				CALL inizializza_dati_fisici()
				CALL inizializza_dati_mc()
				CALL inizializza_walkers('opt_Rp')
				CALL inizializza_dati_funzione_onda()
				IF (mpi_myrank==0) WRITE (istring, '(I4.4)') contatore
				IF (mpi_myrank==0) CALL stampa_file_Rp('reticolo/SR_Rp-'//istring//'.d')
            IF (((Jee_kind=='spl').OR.(Jee_kind=='spp')).AND.(opt_A_Jee.OR.opt_F_Jee)) THEN
               IF (mpi_myrank==0) THEN
                  CALL MSPL_print_on_file(SPL=Jsplee, DERIV=0, FILENAME='ottimizzazione/splines/Jsplee'//istring//'.d',&
                     NPOINTS=INT(N_hist,4) ) 
                  IF (split_Aee.OR.split_Fee) CALL MSPL_print_on_file(SPL=Jsplee_ud, DERIV=0,&
                     FILENAME='ottimizzazione/splines/Jsplee_ud.'//istring, NPOINTS=INT(N_hist,4) )
               END IF
            END IF
            IF (((Jep_kind=='spl').OR.(Jep_kind=='spp')).AND.(opt_A_Jep.OR.opt_F_Jep)) THEN
               IF (mpi_myrank==0) THEN
                  CALL MSPL_print_on_file(SPL=Jsplep, DERIV=0, FILENAME='ottimizzazione/splines/Jsplep'//istring//'.d',&
                     NPOINTS=INT(N_hist,4) ) 
                  IF (split_Aep.OR.split_Fep) CALL MSPL_print_on_file(SPL=Jsplep_ud, DERIV=0,&
                     FILENAME='ottimizzazione/splines/Jsplep_ud.'//istring, NPOINTS=INT(N_hist,4) )
               END IF
            END IF
            IF (((SDe_kind=='spb')).AND.(opt_SDe)) THEN
               IF (mpi_myrank==0) THEN
                  CALL MSPL_print_on_file(SPL=Bsplep, DERIV=0, FILENAME='ottimizzazione/splines/Bsplep'//istring//'.d',&
                     NPOINTS=INT(N_hist,4) ) 
               END IF
            END IF
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

      !stampo le spline iniziali
		IF (stampa_dati_funzione_onda) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
			CALL inizializza_dati_fisici()
			CALL inizializza_dati_mc()
			CALL inizializza_walkers('opt_Rp')
			CALL inizializza_dati_funzione_onda()
			IF (mpi_myrank==0) WRITE (istring, '(I4.4)') contatore
         IF (((Jee_kind=='spl').OR.(Jee_kind=='spp')).AND.(opt_A_Jee.OR.opt_F_Jee)) THEN
            IF (mpi_myrank==0) THEN
               CALL MSPL_print_on_file(SPL=Jsplee, DERIV=0, FILENAME='ottimizzazione/splines/Jsplee'//istring//'.d',&
                  NPOINTS=INT(N_hist,4) ) 
               IF (split_Aee.OR.split_Fee) CALL MSPL_print_on_file(SPL=Jsplee_ud, DERIV=0,&
                  FILENAME='ottimizzazione/splines/Jsplee_ud.'//istring, NPOINTS=INT(N_hist,4) )
            END IF
         END IF
         IF (((Jep_kind=='spl').OR.(Jep_kind=='spp')).AND.(opt_A_Jep.OR.opt_F_Jep)) THEN
            IF (mpi_myrank==0) THEN
               CALL MSPL_print_on_file(SPL=Jsplep, DERIV=0, FILENAME='ottimizzazione/splines/Jsplep'//istring//'.d',&
                  NPOINTS=INT(N_hist,4) ) 
               IF (split_Aep.OR.split_Fep) CALL MSPL_print_on_file(SPL=Jsplep_ud, DERIV=0,&
                  FILENAME='ottimizzazione/splines/Jsplep_ud.'//istring, NPOINTS=INT(N_hist,4) )
            END IF
         END IF
         IF (((SDe_kind=='spb')).AND.(opt_SDe)) THEN
            IF (mpi_myrank==0) THEN
               CALL MSPL_print_on_file(SPL=Bsplep, DERIV=0, FILENAME='ottimizzazione/splines/Bsplep'//istring//'.d',&
                  NPOINTS=INT(N_hist,4) ) 
            END IF
         END IF
			CALL chiudi_dati_funzione_onda()
			CALL chiudi_walkers()
			CALL chiudi_dati_mc()
			CALL chiudi_dati_fisici()
		END IF
		
		!IF (N-num_coord_Rp>0) lambda2=lambda2*DSQRT((DOT_PRODUCT(p0(1:N-num_coord_Rp),p0(1:N-num_coord_Rp)))/REAL(N-num_coord_Rp,8))
		!IF (num_coord_Rp>0) lambda2_Rp=lambda2_Rp*DSQRT(num_coord_Rp*MINVAL(L)*0.01d0)
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

         !CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr) 
         !CALL SLEEP(1)
         !PRINT *, mpi_myrank, p0
         !CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr) 
         !CALL SLEEP(1)
         !IF (mpi_myrank==0) PRINT *, "> > > DISTANZA pp = ", DSQRT(DOT_PRODUCT(p0(1:3)-p0(4:6),p0(1:3)-p0(4:6)))
         
			contatore=contatore+1
         WRITE (istring, '(I4.4)') contatore
         id='SR'//istring
			CALL SR_campiona_wf_attuale(id,p0,N,energia,accettabile)
         IF (mpi_myrank==0) WRITE (8, *) cont_step, energia(1:2)
			IF (mpi_myrank==0) CLOSE (8)
			IF (mpi_myrank==0) OPEN (UNIT=8, FILE='ottimizzazione/SR_energies.dat', STATUS='OLD', POSITION='APPEND')
			IF (mpi_myrank==0) WRITE (9, *) p0
			IF (mpi_myrank==0) CLOSE (9)
			IF (mpi_myrank==0) OPEN (UNIT=9, FILE='ottimizzazione/SR_var_parameters.dat', STATUS='OLD', POSITION='APPEND')
			IF ((.NOT. accettabile).AND.(mpi_myrank==0)) STOP 'Campionamento non accettabile [ module_variational_opt.f90 > stochastic_reconfiguration ]'

         !salvo le ultime 5 energie
			DO i = 5, 2, -1
				old_E(1:2,i)=old_E(1:2,i-1)
			END DO
			old_E(1:2,1)=energia(1:2)
			
			IF ((energia(1)<old_Emin(1)).OR.(contatore==1)) THEN
				old_Emin(1:2)=energia(1:2)
				p0_minimo(1:N)=p0(1:N)
				flag_trovato_minimo=.TRUE.
				IF (mpi_myrank==0) PRINT * , 'VAR_OPT: Trovato un nuovo minimo assoluto.'
				IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *) 'VAR_OPT: Trovato un nuovo minimo assoluto.'
				IF (stampa_dati_funzione_onda) THEN
               CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
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
			END IF

         !Stampo il numero effettivo di parametri variazionali
         IF (num_par_var_eff/=num_par_var) THEN
            IF (mpi_myrank==0) PRINT * , 'VAR_OPT: Numero parametri variazionali effettivi: ', num_par_var_eff
            IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *) &
               'VAR_OPT: Numero parametri variazionali effettivi: ', num_par_var_eff
         END IF

         !regolo il beta
         IF (SR_adaptative_beta.AND.(AV_cont>1)) THEN
            f_SR_beta=f_SR_beta*(0.9d0+MIN(DBLE(AV_cont),5.d0)*0.1d0)
            IF (mpi_myrank==0) PRINT '(1X,A33,I2.2,A1)', "VAR_OPT: [BETA] aumento beta del ",MIN(AV_cont-1,5)*10,"%"
            IF ((flag_output).AND.(mpi_myrank==0)) THEN
               WRITE (7,'(1X,A33,I2.2,A1)') "VAR_OPT: [BETA] aumento beta del ",MIN(AV_cont-1,5)*5,"%"
            END IF
         END IF

         !Calcolo il dp
         SELECT CASE(SR_kind)
         !CASE("Hgrad_min")
         !   CALL gradH_trova_dp(num_par_var_eff,num_par_var_eff-num_coord_Rp,num_coord_Rp,dp)
         !CASE("Lagr_molt")
         !   CALL moltiplicatori_lagrange_trova_dp(num_par_var_eff,num_par_var_eff-num_coord_Rp,num_coord_Rp,dp)
         !CASE("LinMolLag")
         !   CALL lin_molt_lagr(num_par_var_eff,num_par_var_eff-num_coord_Rp,num_coord_Rp,dp)
         CASE("LagrDynam")
            CALL lagrangian_dynamic_forces(num_par_var_eff,num_coord_L,num_coord_Rp,dp)
            IF (mpi_myrank==0 .and. opt_rp) THEN
                WRITE (istring, '(I4.4)') contatore - 1
                CALL stampa_file_vettores3D('reticolo/LagrDyn_Frp-'//istring//'.d', &
                RESHAPE(dp(num_par_var-num_coord_rp+1 : num_par_var),(/ 3, num_coord_rp/3 /)), num_coord_rp/3)
                WRITE (istring, '(I4.4)') contatore
            END IF
         !CASE("fast_SR__")
         !   !Determino il beta ottimale usando il principio variazionale stimando la E successiva tramite un'espansione
         !   !IF (mpi_myrank==0) OPEN(UNIT=666,FILE="ottimizzazione/SR_beta."//istring,STATUS='UNKNOWN',POSITION='ASIS')
         !   SRopt_count=0.d0
         !   SRtargetmin=HUGE(0.d0)
         !   static_f_SR_beta=f_SR_beta
         !   static_SVD_MIN=SVD_MIN
         !   static_lambda3=lambda3
         !   var_change_min=HUGE(0.d0)
         !   DO i = 0, 100, 1
         !      IF (mpi_myrank==MOD(i,mpi_nprocs)) THEN
         !         f_SR_beta = static_f_SR_beta * ( 1.1d0**( ((-1.d0)**i) * (i/2) ) )
         !         CALL calcola_skn_fk()
         !         DO j = 0, 100, 1
         !            SVD_MIN = MIN(SR_max_SVD_MIN , static_SVD_MIN * ( 1.1d0**( ((-1.d0)**j) * (j/2) ) ) )
         !            CALL plain_SR_trova_dp(dpnext)
         !            DO k = 0, 100, 1
         !               lambda3 = MAX(static_lambda3 * ( 1.1d0**( ((-1.d0)**k) * (k/2) ) ) , 1.d0 )
         !               dpnext=lambda3*dpnext
         !               var_change=DSQRT(DOT_PRODUCT(dpnext,dpnext)/DOT_PRODUCT(p0,p0))*100.d0
         !               !var_change=DSQRT(MAXVAL(dpnext*dpnext/(p0*p0)))*100.d0
         !               IF ((i==0).AND.(j==0).AND.(k==0)) dp0=dp
         !               IF ( ( (SR_change_bound.AND.(var_change<SR_max_change)) .OR. &
         !                (.NOT.SR_change_bound) ) .AND. (DABS(DOT_PRODUCT(dpnext,Oi))<SR_maxdeltaPsi) ) THEN
         !                  Hnext=( H+DOT_PRODUCT(dpnext,HOi) )/( 1.d0+DOT_PRODUCT(dpnext,Oi) )
         !                  SRtarget=Hnext
         !                  !WRITE(UNIT=666, FMT=*) var_change,f_SR_beta,SVD_MIN,lambda,Hnext,SRtarget
         !                  IF (SRtarget<SRtargetmin) THEN
         !                     SRtargetmin=SRtarget
         !                     dp=dpnext
         !                     opt_f_SR_beta=f_SR_beta
         !                     opt_SVD_MIN=SVD_MIN
         !                     opt_lambda3=lambda3
         !                     SRopt_count=1.d0
         !                     var_change_min=var_change
         !                  ELSE IF ((SRtarget==SRtargetmin).AND.(var_change<var_change_min)) THEN
         !                     dp=dp+dpnext
         !                     opt_f_SR_beta=opt_f_SR_beta+f_SR_beta
         !                     opt_SVD_MIN=opt_SVD_MIN+SVD_MIN
         !                     opt_lambda3=opt_lambda3+lambda3
         !                     SRopt_count=SRopt_count+1.d0

         !                     !SRtargetmin=SRtarget
         !                     !dp=dpnext
         !                     !opt_f_SR_beta=f_SR_beta
         !                     !opt_SVD_MIN=SVD_MIN
         !                     !opt_lambda3=lambda3
         !                     !SRopt_count=1.d0
         !                     !var_change_min=var_change
         !                  END IF
         !               END IF
         !            END DO
         !         END DO
         !      END IF
         !   END DO

         !   !check if at least one SR configuration has been accepted
         !   !PRINT *, mpi_myrank, SRopt_count
         !   CALL MPI_reduce(SRopt_count,foo,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
         !   !IF (mpi_myrank==0) PRINT *, "foo=",foo
         !   CALL MPI_BCAST(foo,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)

         !   IF (foo>0.5d0) THEN

         !      dp=dp/SRopt_count
         !      opt_f_SR_beta=opt_f_SR_beta/SRopt_count
         !      opt_SVD_MIN=opt_SVD_MIN/SRopt_count
         !      opt_lambda3=opt_lambda3/SRopt_count

         !      !PRINT *, mpi_myrank, SRtargetmin
         !      !Determino quale processore ha trovato il target minimo
         !      CALL MPI_GATHER(SRtargetmin,1,MPI_REAL8,array_SRtargetmin(0:mpi_nprocs-1),1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
         !      IF (mpi_myrank==0) THEN
         !         i=MINLOC(array_SRtargetmin(0:mpi_nprocs-1),1)-1
         !         !PRINT *, array_SRtargetmin
         !         !PRINT *, "MINIMO ", MINVAL(array_SRtargetmin(0:mpi_nprocs-1))
         !         !PRINT *, "MINIMO IN ", i
         !         !IF ((i<0).OR.(i>mpi_nprocs-1)) PRINT *, "ERRORE i fuori dal range"
         !         !CLOSE(666)
         !      END IF
         !      CALL MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)

         !      !aggiorno i valori con quelli ottimali per il processore con il target piu' basso
         !      IF (mpi_myrank==i) THEN
         !         f_SR_beta=opt_f_SR_beta
         !         SVD_MIN=opt_SVD_MIN
         !         lambda3=opt_lambda3
         !      END IF

         !      !distribuisco i valori ottimali a tutti i processori
         !      CALL MPI_BCAST(dp,num_par_var,MPI_REAL8,i,MPI_COMM_WORLD,mpi_ierr)
         !      CALL MPI_BCAST(f_SR_beta,1,MPI_REAL8,i,MPI_COMM_WORLD,mpi_ierr)
         !      CALL MPI_BCAST(SVD_MIN,1,MPI_REAL8,i,MPI_COMM_WORLD,mpi_ierr)
         !      CALL MPI_BCAST(lambda3,1,MPI_REAL8,i,MPI_COMM_WORLD,mpi_ierr)
         !      CALL MPI_BCAST(SRopt_count,1,MPI_REAL8,i,MPI_COMM_WORLD,mpi_ierr)

         !      !stampo le informazioni
         !      IF (mpi_myrank==0) THEN
         !         PRINT *, "VAR_OPT: [fSR] usato fSR per determinare i parametri SR ottimali. Numero configurazioni:", &
         !          INT2(SRopt_count)
         !         IF (flag_output) WRITE(UNIT=7, FMT=*) &
         !            "VAR_OPT: [fSR] usato fSR per determinare i parametri SR ottimali. Numero configurazioni:", &
         !            INT2(SRopt_count)
         !         PRINT *, "VAR_OPT: [BETA] regolato beta tramite metodo di proiezione: ", REAL(f_SR_beta*SR_beta,4)
         !         IF (flag_output) WRITE(UNIT=7, FMT=*) &
         !            "VAR_OPT: [BETA] regolato beta tramite metodo di proiezione: ", REAL(f_SR_beta*SR_beta,4)
         !         PRINT *, "VAR_OPT: [SVD_MIN] regolato SVD_MIN tramite metodo di proiezione: ", REAL(SVD_MIN,4)
         !         IF (flag_output) WRITE(UNIT=7, FMT=*) &
         !            "VAR_OPT: [SVD_MIN] regolato SVD_MIN tramite metodo di proiezione: ", REAL(SVD_MIN,4)
         !         PRINT *, "VAR_OPT: [lambda3] regolato lambda3 tramite metodo di proiezione: ", REAL(lambda3,4)
         !         IF (flag_output) WRITE(UNIT=7, FMT=*) &
         !            "VAR_OPT: [lambda3] regolato lambda3 tramite metodo di proiezione: ", REAL(lambda3,4)
         !      END IF
         !   ELSE
         !      IF (mpi_myrank==0) PRINT *, "VAR_OPT: [fSR] nessuna configurazione ottimale trovata,",&
         !        " quindi uso i parametri dati in input"
         !      IF (flag_output) WRITE(UNIT=7, FMT=*) &
         !         "VAR_OPT: [fSR] nessuna configurazione ottimale trovata, quindi uso i parametri dati in input"
         !      IF (mpi_myrank==0) dp=dp
         !      CALL MPI_BCAST(dp,num_par_var,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
         !   END IF

         CASE("plain_SR_")

         !   IF (mpi_myrank==0) THEN
         !      CALL calcola_skn_fk()
         !      CALL plain_SR_trova_dp(dp)
         !   END IF
         !   CALL MPI_BCAST(dp,num_par_var,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)

         CASE DEFAULT
            IF (mpi_myrank==0) PRINT *, "Errore: SR_kind non accettabile"
            IF (mpi_myrank==0) PRINT *, "[module_variational_opt.f90 > stochastic_reconfiguration]"
            IF (flag_output) WRITE(UNIT=7, FMT=*) "Errore: SR_kind non accettabile"
            IF (flag_output) WRITE(UNIT=7, FMT=*) "[module_variational_opt.f90 > stochastic_reconfiguration]"
            CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
            STOP
         END SELECT

         vec_app=dp/DSQRT(DOT_PRODUCT(dp,dp))

         CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)

         !salva il dp iniziale e il primo dp_old
			IF ( contatore==1 ) THEN
				dp_iniz=dp
            dp_old=dp
         END IF
			
         !parte relativa al lambda

         !regolo lambda a seconda di come sono stati gli ultimi spostamenti variazionali
         IF (SR_lambda) THEN

			   IF ((SR_lambda).AND.(contatore>2).AND.(N-num_coord_Rp>0)) THEN       !parametri variazionali wf
			   	flambda=(1.d0+0.1d0*DOT_PRODUCT(vec_app(1:N-num_coord_Rp),vec_app_old(1:N-num_coord_Rp)))
			   	IF (contatore<5) flambda=(1.d0+0.50d0*DOT_PRODUCT(vec_app(1:N-num_coord_Rp),vec_app_old(1:N-num_coord_Rp)))
			   	!flambda=(1.d0+0.5d0*DOT_PRODUCT(vec_app(1:N-num_coord_Rp),vec_app_old(1:N-num_coord_Rp)))
               IF (SR_lambda) lambda=lambda*flambda
			   	IF (flambda>1.d0) THEN
			   		loop_stop(1)=.FALSE.
			   	ELSE
			   		loop_stop(1)=.TRUE.
			   	END IF
			   	IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [PS] cambio lambda di', flambda
			   	IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *) 'VAR_OPT: [PS] cambio lambda di', flambda
			   END IF
			   IF ((SR_lambda_Rp).AND.(contatore>2).AND.(num_coord_Rp>0)) THEN       !parametri Rp
			   	flambda_Rp=(1.d0+0.1d0*DOT_PRODUCT(vec_app(N-num_coord_Rp+1:N),vec_app_old(N-num_coord_Rp+1:N)))
			   	IF (contatore<5) flambda_Rp=(1.d0+0.50d0*DOT_PRODUCT(vec_app(N-num_coord_Rp+1:N),vec_app_old(N-num_coord_Rp+1:N)))
               IF (SR_lambda_Rp) lambda_Rp=lambda_Rp*flambda_Rp
			   	IF (flambda_Rp>1.d0) THEN
			   		loop_stop(1)=.FALSE.
			   	ELSE
			   		loop_stop(1)=.TRUE.
			   	END IF
			   	IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [PS] cambio lambda_Rp di', flambda_Rp
			   	IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *) 'VAR_OPT: [PS] cambio lambda_Rp di', flambda_Rp
			   END IF
			   vec_app_old=vec_app
			   
			   IF ((SR_lambda).AND.(contatore>3).AND.(N-num_coord_Rp>0)) THEN          !parametri variazionali wf
			   	dist_poss=SUM(old_step(1:3))
			   	dist_perc=DSQRT(DOT_PRODUCT(old_p0(1:N-num_coord_Rp,1)-old_p0(1:N-num_coord_Rp,3), &
			   	                            old_p0(1:N-num_coord_Rp,1)-old_p0(1:N-num_coord_Rp,3)))
			   	flambda=(0.85d0+0.3d0*dist_perc/dist_poss)
			   	!flambda=(0.5d0+1.d0*dist_perc/dist_poss)
			   	IF (SR_lambda) lambda=lambda*flambda
			   	IF (flambda>1.d0) THEN
			   		loop_stop(2)=.FALSE.
			   	ELSE
			   		loop_stop(2)=.TRUE.
			   	END IF
			   	IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DIST3] cambio lambda di', flambda
			   	IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *) 'VAR_OPT: [DIST3] cambio lambda di', flambda
			   END IF
			   IF ((SR_lambda_Rp).AND.(contatore>3).AND.(num_coord_Rp>0)) THEN       !parametri Rp
			   	dist_poss=SUM(old_step_Rp(1:3))
			   	dist_perc=DSQRT(DOT_PRODUCT(old_p0(N-num_coord_Rp+1:N,1)-old_p0(N-num_coord_Rp+1:N,3), &
			   	                            old_p0(N-num_coord_Rp+1:N,1)-old_p0(N-num_coord_Rp+1:N,3)))
			   	flambda_Rp=(0.85d0+0.3d0*dist_perc/dist_poss)
               IF (SR_lambda_Rp) lambda_Rp=lambda_Rp*flambda_Rp
			   	IF (flambda_Rp>1.d0) THEN
			   		loop_stop(2)=.FALSE.
			   	ELSE
			   		loop_stop(2)=.TRUE.
			   	END IF
			   	IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DIST3] cambio lambda_Rp di', flambda_Rp
			   	IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *) 'VAR_OPT: [DIST3] cambio lambda_Rp di', flambda_Rp
			   END IF
			   
			   IF ((SR_lambda).AND.(contatore>5).AND.(N-num_coord_Rp>0)) THEN          !parametri variazionali wf
			   	dist_poss_5=SUM(old_step(1:5))
			   	dist_perc_5=DSQRT(DOT_PRODUCT(old_p0(1:N-num_coord_Rp,1)-old_p0(1:N-num_coord_Rp,5), &
			   	                              old_p0(1:N-num_coord_Rp,1)-old_p0(1:N-num_coord_Rp,5)))
               flambda=(0.75d0+0.5d0*dist_perc_5/dist_poss_5)
			   	IF (SR_lambda) lambda=lambda*flambda
			   	IF (flambda>1.d0) THEN
			   		loop_stop(3)=.FALSE.
			   	ELSE
			   		loop_stop(3)=.TRUE.
			   	END IF
			   	IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DIST5] cambio lambda di', flambda
			   	IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *) 'VAR_OPT: [DIST5] cambio lambda di', flambda
			   END IF
			   IF ((SR_lambda_Rp).AND.(contatore>5).AND.(num_coord_Rp>0)) THEN       !parametri Rp
			   	dist_poss_5=SUM(old_step_Rp(1:5))
			   	dist_perc_5=DSQRT(DOT_PRODUCT(old_p0(N-num_coord_Rp+1:N,1)-old_p0(N-num_coord_Rp+1:N,5), &
			   	                              old_p0(N-num_coord_Rp+1:N,1)-old_p0(N-num_coord_Rp+1:N,5)))
			   	flambda_Rp=(0.75d0+0.50d0*dist_perc_5/dist_poss_5)
               IF (SR_lambda_Rp) lambda_Rp=lambda_Rp*flambda_Rp
			   	IF (flambda_Rp>1.d0) THEN
			   		loop_stop(3)=.FALSE.
			   	ELSE
			   		loop_stop(3)=.TRUE.
			   	END IF
			   	IF (mpi_myrank==0) PRINT * , 'VAR_OPT: [DIST5] cambio lambda_Rp di', flambda_Rp
			   	IF ((flag_output).AND.(mpi_myrank==0)) WRITE (7, *) 'VAR_OPT: [DIST5] cambio lambda_Rp di', flambda_Rp
			   END IF

            !limito lambda fra min_lambda e max_lambda
            lambda=MIN(MAX(min_lambda,lambda),max_lambda)
            !limito lambda_Rp fra min_lambda_Rp e max_lambda_Rp
            lambda_Rp=MIN(MAX(min_lambda_Rp,lambda_Rp),max_lambda_Rp)
         END IF

         !Controllo manuale sul cambio in percentuale dei parametri
         lambda2=1.d0
         lambda2_Rp=1.d0
         IF (SR_change_bound) THEN
            !setto il lambda in modo che il cambio dei parametri variazionali non sia maggiore di SR_max_change%
            IF ( lambda*lambda2 > DSQRT(((0.01d0*SR_max_change)**2)*&
                 DOT_PRODUCT(p0(1:N-num_coord_Rp),p0(1:N-num_coord_Rp))/&
                 DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp))) ) THEN
               lambda2=(0.01d0*SR_max_change/lambda)*&
                  DSQRT( DOT_PRODUCT(p0(1:N-num_coord_Rp),p0(1:N-num_coord_Rp))/&
                         DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp)) )
               IF (mpi_myrank==0) THEN
                  PRINT *, &
                     "VAR_OPT: [change_bound] Cambio dei parametri sarebbe troppo grande, lo riduco di ", lambda2
                  IF (flag_output) WRITE(UNIT=7, FMT=*) &
                     "VAR_OPT: [change_bound] Cambio dei parametri sarebbe troppo grande, lo riduco di ", lambda2
               END IF
            END IF
            !setto il lambda in modo che il cambio dei parametri variazionali non sia minore di SR_min_change%
            IF ( lambda*lambda2 < DSQRT(((0.01d0*SR_min_change)**2)*&
                 DOT_PRODUCT(p0(1:N-num_coord_Rp),p0(1:N-num_coord_Rp))/&
                 DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp))) ) THEN
               lambda2=(0.01d0*SR_min_change/lambda)*&
                  DSQRT( DOT_PRODUCT(p0(1:N-num_coord_Rp),p0(1:N-num_coord_Rp))/&
                         DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp)) )
               IF (mpi_myrank==0) THEN
                  PRINT *, &
                     "VAR_OPT: [change_bound] Cambio dei parametri sarebbe troppo piccolo, lo aumento di ", lambda2
                  IF (flag_output) WRITE(UNIT=7, FMT=*) &
                     "VAR_OPT: [change_bound] Cambio dei parametri sarebbe troppo piccolo, lo aumento di ", lambda2
               END IF
            END IF
         END IF
         IF (SR_change_bound_Rp) THEN
            !setto il lambda_Rp in modo che il cambio delle posizioni protoniche non sia maggiore di SR_max_change%
            IF ( lambda_Rp*lambda2_Rp > DSQRT(((0.01d0*SR_max_change)**2)*&
                 DOT_PRODUCT(p0(N-num_coord_Rp+1:N),p0(N-num_coord_Rp+1:N))/&
                 DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N))) ) THEN
               lambda2_Rp=(0.01d0*SR_max_change/lambda_Rp)*&
                  DSQRT( DOT_PRODUCT(p0(N-num_coord_Rp+1:N),p0(N-num_coord_Rp+1:N))/&
                  DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N)) )
               IF (mpi_myrank==0) THEN
                  PRINT *, &
                     "VAR_OPT: [change_bound] Cambio delle posizioni protoniche troppo grande, lo riduco di ", lambda2_Rp
                  IF (flag_output) WRITE(UNIT=7, FMT=*) &
                     "VAR_OPT: [change_bound] Cambio delle posizioni protoniche troppo grande, lo riduco di ", lambda2_Rp
               END IF
            END IF
            !setto il lambda in modo che il cambio delle posizioni protoniche non sia minore di SR_min_change%
            IF ( lambda_Rp*lambda2_Rp < DSQRT(((0.01d0*SR_min_change)**2)*&
                 DOT_PRODUCT(p0(N-num_coord_Rp+1:N),p0(N-num_coord_Rp+1:N))/&
                 DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N))) ) THEN
               lambda2_Rp=(0.01d0*SR_min_change/lambda_Rp)*&
                  DSQRT( DOT_PRODUCT(p0(N-num_coord_Rp+1:N),p0(N-num_coord_Rp+1:N))/&
                  DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N)) )
               IF (mpi_myrank==0) THEN
                  PRINT *, &
                     "VAR_OPT: [change_bound] Cambio delle posizioni protoniche troppo piccolo, lo aumento di ", lambda2_Rp
                  IF (flag_output) WRITE(UNIT=7, FMT=*) &
                     "VAR_OPT: [change_bound] Cambio delle posizioni protoniche troppo piccolo, lo aumento di ", lambda2_Rp
               END IF
            END IF
         END IF
			
			!moltiplico lo spostamento variazionale per le lambda [, lo stampo *eliminato*] e salvo gli old_step e cont_step
         IF (N-num_coord_Rp>0) THEN
            dp(1:N-num_coord_Rp)=dp(1:N-num_coord_Rp)*lambda*lambda2
			   old_step(5)=old_step(4)
			   old_step(4)=old_step(3)
			   old_step(3)=old_step(2)
			   old_step(2)=old_step(1)
			   old_step(1)=DSQRT(DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp)))
			   cont_step=cont_step+old_step(1)
         END IF
         IF (num_coord_Rp>0) THEN
            dp(N-num_coord_Rp+1:N)=dp(N-num_coord_Rp+1:N)*lambda_Rp*lambda2_Rp
			   old_step_Rp(5)=old_step_Rp(4)
			   old_step_Rp(4)=old_step_Rp(3)
			   old_step_Rp(3)=old_step_Rp(2)
			   old_step_Rp(2)=old_step_Rp(1)
			   old_step_Rp(1)=DSQRT(DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N)))
			   cont_step=cont_step+old_step_Rp(1)
         END IF

         !Stampo alcune informazioni riguardo SR
         IF (mpi_myrank==0) THEN
            PRINT '(5X,A92)', "- - - - - - - - - - - - - - - - - - - - - -"//&
              " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
            PRINT '(5X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X, A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1)', &
              "|","beta","|","SVD_MIN","|","lambda3","|","lambda2","|","lambda","|","lambda2_Rp","|","lambda_Rp","|"
            PRINT '(5X,A1,1X,E10.4,1X,A1,1X,E10.4,1X,A1,1X,E10.4,1X,A1, 1X,E10.4,1X,A1,1X,E10.4,1X,A1,1X,E10.4,1X,A1,1X,E10.4,1X,A1)', &
               "|",f_SR_beta*SR_beta,"|",SVD_MIN,"|",lambda3,"|",lambda2,"|",lambda,"|",lambda2_Rp,"|",lambda_Rp,"|"
            PRINT '(5X,A92)', "- - - - - - - - - - - - - - - - - - - - - -"//&
               " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
            WRITE(7, '(5X,A92)') "- - - - - - - - - - - - - - - - - - - - - - -"//&
               " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
            WRITE(7, '(5X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1, 1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1)') &
              "|","beta","|","SVD_MIN","|","lambda3","|","lambda2","|","lambda","|","lambda2_Rp","|","lambda_Rp","|"
            WRITE(7, '(5X,A1,1X,E10.4,1X,A1,1X,E10.4,1X,A1,1X,E10.4,1X,A1, 1X,E10.4,1X,A1,1X,E10.4,1X,A1,1X,E10.4,1X,A1,1X,E10.4,1X,A1)') &
               "|",f_SR_beta*SR_beta,"|",SVD_MIN,"|",lambda3,"|",lambda2,"|",lambda,"|",lambda2_Rp,"|",lambda_Rp,"|"
            WRITE(7, '(5X,A92)') "- - - - - - - - - - - - - - - - - - - - - - -"//&
               " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
         END IF

         !calcolo il cambiamento percentuale dei parametri
         IF (mpi_myrank==0) THEN
            PRINT '(14X,A38,F5.1,A9)', "> > >   VAR_OPT: Cambio dei parametri ",&
               DSQRT(DOT_PRODUCT(dp,dp)/DOT_PRODUCT(p0,p0))*100.d0,"%   < < <"
            IF (flag_output) WRITE(UNIT=7, FMT='(14X,A38,F5.1,A9)') "> > >   VAR_OPT: Cambio dei parametri ",&
               DSQRT(DOT_PRODUCT(dp,dp)/DOT_PRODUCT(p0,p0))*100.d0,"%   < < <"
            IF (N-num_coord_Rp>0) THEN
               PRINT '(14X,A38,F5.1,A9)', "> > >          Parametri variazionali ",&
                  DSQRT(DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp))/&
                  DOT_PRODUCT(p0(1:N-num_coord_Rp),p0(1:N-num_coord_Rp)))*100.d0,"%   < < <"
               IF (flag_output) WRITE(UNIT=7, FMT='(14X,A38,F5.1,A9)') &
                  "> > >          Parametri variazionali ",&
                  DSQRT(DOT_PRODUCT(dp(1:N-num_coord_Rp),dp(1:N-num_coord_Rp))/&
                  DOT_PRODUCT(p0(1:N-num_coord_Rp),p0(1:N-num_coord_Rp)))*100.d0,"%   < < <"
            END IF
            IF (num_coord_Rp>0) THEN
               PRINT '(14X,A38,F5.1,A9)', "> > >            Posizioni protoniche ",&
                  DSQRT(DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N))/&
                  DOT_PRODUCT(p0(N-num_coord_Rp+1:N),p0(N-num_coord_Rp+1:N)))*100.d0,"%   < < <"
               IF (flag_output) WRITE(UNIT=7, FMT='(14X,A38,F5.1,A9)') &
                  "> > >            Posizioni protoniche ",&
                  DSQRT(DOT_PRODUCT(dp(N-num_coord_Rp+1:N),dp(N-num_coord_Rp+1:N))/&
                  DOT_PRODUCT(p0(N-num_coord_Rp+1:N),p0(N-num_coord_Rp+1:N)))*100.d0,"%   < < <"
            END IF
         END IF
		   
         !controllo da quanto sono nello stesso minimo e in caso fermo il loop
         IF (energia(1)+2.d0*energia(2)<minE_AV(1)-2.d0*minE_AV(2)) THEN
            !nuovo minimo
            WO_MIN_cont=0
         ELSE
            WO_MIN_cont=WO_MIN_cont+1
         END IF
         IF (WO_MIN_cont>SR_num_max_WO_MIN .AND. SR_num_max_WO_MIN>=0) flag_loop=.FALSE.
         IF (mpi_myrank==0) THEN
            IF (WO_MIN_cont<=SR_num_max_WO_MIN) THEN
               PRINT *, "VAR_OPT: [LOOP] WO_MIN_count = ", INT2(WO_MIN_cont)
               WRITE(UNIT=7, FMT=*) "VAR_OPT: [LOOP] WO_MIN_count = ", INT2(WO_MIN_cont)
            END IF
            IF (WO_MIN_cont>SR_num_max_WO_MIN) THEN
               PRINT *, "VAR_OPT: [LOOP] Raggiunto il limite di SR-steps senza nuovo minimo. FINE SR. "
               WRITE(UNIT=7, FMT=*) "VAR_OPT: [LOOP] Raggiunto il limite di SR-steps senza nuovo minimo. FINE SR. "
            END IF
         END IF

         !se tutto e' ok ed e' necessario un altro passo SR
			IF (flag_loop) THEN

            !Salvo il volume
            IF (opt_L) THEN
               foo=L(1)*L(2)*L(3)
               IF (flag_2D) foo=L(1)*L(2)
            END IF

            !aggiorno i parametri variazionali
            p0=p0+dp

            IF ( opt_rp ) THEN
                call costrizione_rp(p0(num_par_var-num_coord_rp+1 : num_par_var))
            END IF

            !Se L e' stato cambiato, faccio in modo che il volume rimanga invariato
            IF (opt_L) THEN
               IF (flag_2D) THEN
                  foo=DSQRT(foo/(p0(1)*p0(2)))
                  p0(1:2)=p0(1:2)*foo
               ELSE
                  foo=(foo/(p0(1)*p0(2)*p0(3)))**(1./3.)
                  p0(1:3)=p0(1:3)*foo
               END IF
            END IF
            
            !parte che gestisce AVERAGE
				IF ( AV_accu ) THEN
               !se AV_accu=T (stoc_av) allora accumulo tutto indiscriminatamente
					AV_P=AV_P+p0
					AV_cont=AV_cont+1
				ELSE
               !altrimenti, distinguo

               !se l'energia e` piu bassa di 2 sigma del minimo precedente, resetto i valori e ricomincio ad accumulare
               IF (energia(1)+2.d0*energia(2)<minE_AV(1)-2.d0*minE_AV(2)) THEN
                  minE_AV=energia
                  AV_P=p0
                  AV_cont=1
                  WO_MIN_cont=0
                  IF (mpi_myrank==0) PRINT *, "VAR_OPT: [AV] Nuovo minimo assoluto, resetto l'AVErage."
                  IF (flag_output.AND.(mpi_myrank==0)) WRITE(UNIT=7, FMT=*) &
                     "VAR_OPT: [AV] Nuovo minimo assoluto, resetto l'AVErage."

               !altrimenti, se l'energia e' rimasta la stessa, accumulo i valori in AV_P
               ELSE IF (energia(1)-2.d0*energia(2)<minE_AV(1)+2.d0*minE_AV(2)) THEN
                  AV_P=AV_P+p0
                  AV_cont=AV_cont+1
               !altrimenti, se l'energia e' troppo alta, non accumulo in AV_P
               ELSE
                  IF (mpi_myrank==0) THEN
                     PRINT *,"VAR_OPT: [AV] Energia troppo alta, non la considero per l'AVerage"
                     WRITE (7, *) "VAR_OPT: [AV] Energia troppo alta, non considero la configurazione per l'AVerage"
                  END IF
               END IF
				END IF
            !salvo i parametri old_p0
			   old_p0(1:N,5)=old_p0(1:N,4)
			   old_p0(1:N,4)=old_p0(1:N,3)
			   old_p0(1:N,3)=old_p0(1:N,2)
			   old_p0(1:N,2)=old_p0(1:N,1)
			   old_p0(1:N,1)=p0(1:N)

            CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)

			   !salvo su file gli ultimi dati variazionali SR
			   IF (stampa_dati_funzione_onda) THEN
               CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
			   	CALL inizializza_dati_fisici()
			   	CALL inizializza_dati_mc()
			   	CALL inizializza_walkers('opt_Rp')
			   	CALL inizializza_dati_funzione_onda()
			   	CALL setta_parametri(p0,N)
			   	IF (mpi_myrank==0) CALL stampa_file_dati_funzione_onda('ottimizzazione/SR_wf.d')
			   	IF (mpi_myrank==0) WRITE (istring, '(I4.4)') contatore
			   	IF ((mpi_myrank==0).AND.(opt_Rp)) CALL stampa_file_Rp('reticolo/SR_Rp-'//istring//'.d')
               IF (((Jee_kind=='spl').OR.(Jee_kind=='spp')).AND.(opt_A_Jee.OR.opt_F_Jee)) THEN
                  IF (mpi_myrank==0) THEN
                     CALL MSPL_print_on_file(SPL=Jsplee, DERIV=0, FILENAME='ottimizzazione/splines/Jsplee'//istring//'.d',&
                        NPOINTS=INT(N_hist,4) ) 
                     IF (split_Aee.OR.split_Fee) CALL MSPL_print_on_file(SPL=Jsplee_ud, DERIV=0,&
                        FILENAME='ottimizzazione/splines/Jsplee_ud.'//istring, NPOINTS=INT(N_hist,4) )
                  END IF
               END IF
               IF (((Jep_kind=='spl').OR.(Jep_kind=='spp')).AND.(opt_A_Jep.OR.opt_F_Jep)) THEN
                  IF (mpi_myrank==0) THEN
                     CALL MSPL_print_on_file(SPL=Jsplep, DERIV=0, FILENAME='ottimizzazione/splines/Jsplep'//istring//'.d',&
                        NPOINTS=INT(N_hist,4) ) 
                     IF (split_Aep.OR.split_Fep) CALL MSPL_print_on_file(SPL=Jsplep_ud, DERIV=0,&
                        FILENAME='ottimizzazione/splines/Jsplep_ud.'//istring, NPOINTS=INT(N_hist,4) )
                  END IF
               END IF
               IF (((SDe_kind=='spb')).AND.(opt_SDe)) THEN
                  IF (mpi_myrank==0) THEN
                     CALL MSPL_print_on_file(SPL=Bsplep, DERIV=0, FILENAME='ottimizzazione/splines/Bsplep'//istring//'.d',&
                        NPOINTS=INT(N_hist,4) ) 
                  END IF
               END IF
			   	CALL chiudi_dati_fisici()
			   	CALL chiudi_dati_mc()
			   	CALL chiudi_walkers()
			   	CALL chiudi_dati_funzione_onda()
               CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
			   END IF

			   !stampa i dati variazionali AVERAGE
			   IF (stampa_dati_funzione_onda) THEN
               CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
			   	CALL inizializza_dati_fisici()
			   	CALL inizializza_dati_mc()
			   	CALL inizializza_walkers('opt_Rp')
			   	CALL inizializza_dati_funzione_onda()
			   	AV_P=AV_P/REAL(AV_cont,8)
			   	CALL setta_parametri(AV_P,N)
			   	AV_P=AV_P*REAL(AV_cont,8)
			   	IF (mpi_myrank==0) CALL stampa_file_dati_funzione_onda('ottimizzazione/AV_wf.d')
			   	IF (mpi_myrank==0) PRINT *, 'VAR_OPT: [AV] salvato file AV_wf.d, considerando gli ultimi ', INT2(AV_cont), 'punti.'
               IF (mpi_myrank==0) PRINT *, '          --  minE_AV=',minE_AV(1), "+-", minE_AV(2)
			   	IF ((flag_output).AND.(mpi_myrank==0)) THEN
			   		WRITE (7, *) 'VAR_OPT: [AV] salvato file AV_wf.d, considerando gli ultimi ', INT2(AV_cont), 'punti.'
                  			WRITE (7, *) '          --  minE_AV=',minE_AV(1), "+-", minE_AV(2)
			   	END IF
			   	IF ( opt_Rp ) THEN
			   		IF (mpi_myrank==0) CALL stampa_file_Rp('reticolo/AV_Rp.d')
			   	END IF
			   	CALL chiudi_dati_fisici()
			   	CALL chiudi_dati_mc()
			   	CALL chiudi_walkers()
			   	CALL chiudi_dati_funzione_onda()
               CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
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

         END IF !(se flag_loop==T)

			IF (SR_num_max_WO_MIN<0) flag_loop=.FALSE.			

			CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
		END DO
		
		!salvo la g(r) dei Rp
		IF ( ( AV_accu ).AND.(mpi_myrank==0) ) THEN
			DEALLOCATE(gRpuu,gRpuu_new,gRpud,gRpud_new,x_gr)
		END IF
		
		E_min(1:2)=old_Emin(1:2)
		parametri_var=p0_minimo
		IF (mpi_myrank==0) CLOSE (8)
      DEALLOCATE(work)

	END SUBROUTINE stochastic_reconfiguration

   !!forma quadratica della funzione per lo SR, di cui bisogna trovare lo zero per trovare il dp ottimale
   !SUBROUTINE SR_forma_quadratica(dp,f,k)
   !   IMPLICIT NONE
   !   REAL(KIND=8), INTENT(IN) :: dp(1:num_par_var)
   !   REAL(KIND=8), INTENT(OUT) :: f(1:num_par_var)
   !   INTEGER, OPTIONAL :: k
   !   INTEGER :: i, j, kk

   !   IF (PRESENT(k)) THEN
   !      f(k)=f_k(k)
   !      DO i = 1, num_par_var, 1
   !         f(k)=f(k)+dp(i)*s_kn(i,k)
   !      END DO
   !      DO j = 1, num_par_var, 1
   !      DO i = 1, num_par_var, 1
   !         f(k)=f(k)+dp(i)*dp(j)*c_knm(i,j,k)
   !      END DO
   !      END DO
   !   ELSE
   !      DO kk = 1, num_par_var, 1
   !         f(kk)=f_k(kk)
   !         DO i = 1, num_par_var, 1
   !            f(kk)=f(kk)+dp(i)*s_kn(i,kk)
   !         END DO
   !         DO j = 1, num_par_var, 1
   !         DO i = 1, num_par_var, 1
   !            f(kk)=f(kk)+dp(i)*dp(j)*c_knm(i,j,kk)
   !         END DO
   !         END DO
   !      END DO
   !   END IF
   !   
   !END SUBROUTINE SR_forma_quadratica


   !!massimo fra i termini della forma quadratica
   !SUBROUTINE max_SR_forma_quadratica(dp,f_max)
   !   IMPLICIT NONE
   !   REAL(KIND=8), INTENT(IN) :: dp(1:num_par_var)
   !   REAL(KIND=8), INTENT(OUT) :: f_max
   !   REAL(KIND=8) :: f(1:num_par_var)

   !   CALL SR_forma_quadratica(dp,f)
   !   f_max=MAXVAL(f)
   !   
   !END SUBROUTINE max_SR_forma_quadratica

   !!subroutine che prova a minimizzare il massimo della forma quadratica dello SR
   !SUBROUTINE SR_minimizza_forma_quadratica(dp)
   !   IMPLICIT NONE
   !   INTEGER, PARAMETER :: MAX_NO_UPDATE=100
   !   REAL(KIND=8), INTENT(INOUT) :: dp(1:num_par_var)
   !   REAL(KIND=8) :: dp2(1:num_par_var), dp0(1:num_par_var)
   !   REAL(KIND=8) :: q_k(1:num_par_var), grad(1:num_par_var), q_max, q_max2
   !   REAL(KIND=8) :: lambda
   !   INTEGER :: i1
   !   INTEGER :: k_max, contatore, no_update
   !   
   !   dp0=dp

   !   !trova i valori q_k delle forme quadratiche usando i dp ottenuti usando la forma lineare dello SR
   !   CALL SR_forma_quadratica(dp,q_k)

   !   !trova l'indice k_max per cui q_k e' massimo e il suo valore
   !   k_max=MAXLOC(q_k,1)
   !   q_max=q_k(k_max)

   !   !determina il gradiente da seguire
   !   DO i1 = 1, num_par_var, 1
   !      grad(i1)=DOT_PRODUCT(dp,c_knm(k_max,:,i1))+s_kn(k_max,i1)
   !   END DO

   !   lambda=1.d0
   !   contatore=0
   !   no_update=0
   !   DO WHILE (no_update<MAX_NO_UPDATE)

   !      contatore=contatore+1
   !      !PRINT *, contatore, q_max
   !      
   !      !cambio dp seguendo il gradiente
   !      dp2=dp-grad*lambda*DSIGN(1.d0,q_max)

   !      !calcolo il nuovo valore della forma quadratica
   !      CALL max_SR_forma_quadratica(dp2,q_max2)

   !      !se ho ottenuto un miglioramento aggiorno dp, altrimenti diminuisco lambda e ripeto
   !      IF (DABS(q_max2)>=DABS(q_max)) THEN
   !         lambda=lambda*0.75d0
   !         no_update=no_update+1
   !      ELSE
   !         no_update=0
   !         lambda=lambda*1.5d0
   !         dp=dp2
   !         CALL SR_forma_quadratica(dp,q_k)
   !         k_max=MAXLOC(q_k,1)
   !         q_max=q_k(k_max)
   !         DO i1 = 1, num_par_var, 1
   !            grad(i1)=DOT_PRODUCT(dp,c_knm(k_max,:,i1))+s_kn(k_max,i1)
   !         END DO
   !      END IF
   !      
   !   END DO

   !   !PRINT *, 
   !   !PRINT *, "dp0="
   !   !PRINT *, dp0
   !   !PRINT *, 
   !   !PRINT *, "dp="
   !   !PRINT *, dp

   !   !PRINT *, 
   !   !k_max=MAXLOC(DABS(dp-dp0),1)
   !   !PRINT *, dp0(k_max), dp(k_max)
   !   !
   !   !STOP

   !   IF (mpi_myrank==0) THEN
   !      PRINT *, "VAR_OPT: usato lo sviluppo al secondo ordine della funzione onda. Cambio in dp=",&
   !         DSQRT(DOT_PRODUCT(dp-dp0,dp-dp0))/DSQRT(DOT_PRODUCT(dp0,dp0))
	!		IF (flag_output) THEN
	!		   WRITE (7, *) "VAR_OPT: usato lo sviluppo al secondo ordine della funzione onda. Cambio in dp=",&
   !            DSQRT(DOT_PRODUCT(dp-dp0,dp-dp0))/DSQRT(DOT_PRODUCT(dp0,dp0))
   !      END IF
   !   END IF
   !   
   !END SUBROUTINE SR_minimizza_forma_quadratica


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
				IF (flag_output) WRITE (7, *) 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
			END IF
			b=b+first_step
			rb=r0+b*dir
			CALL controlla_punto(rb)
			CALL calcola_energia(idb,rb,N,fb(1:2),accettabile)
			IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 1[ module_variational_opt.f90 > bracket_min_multidim ]'
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
					IF (flag_output) WRITE (7, *) 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
				END IF
				a=a-first_step
				ra=r0+a*dir
				CALL controlla_punto(ra)
				CALL calcola_energia(ida,ra,N,fa(1:2),accettabile)
				IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 2[ module_variational_opt.f90 > bracket_min_multidim ]'
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
					IF (flag_output) WRITE (7, *) 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
				END IF
				!first_step=2.d0*first_step
				c=c+first_step
				rc=r0+c*dir
				CALL controlla_punto(rc)
				CALL calcola_energia(idc,rc,N,fc(1:2),accettabile)
				IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 3 [ module_variational_opt.f90 > bracket_min_multidim ]'
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
				IF ((id/=ida) .AND. (id/=idb) .AND. (id/=idc)) CALL SYSTEM ('export ID_TO_ERASE=opt'//id//'; find estimatori/ -name *$ID_TO_ERASE* -exec rm -f {} \\; find posizioni/ -name *$ID_TO_ERASE* -exec rm -f {} \\')
			END DO
		END IF

		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Ho eseguito il bracket con ', INT(contatore), ' passi.',' I tre punti sono '
			PRINT * , '         punti:   ', a, b ,c
			PRINT * , '         energie: ', fa(1), fb(1), fc(1)
			IF (flag_output) THEN
				WRITE (7, *) 'VAR_OPT: Ho eseguito il bracket con ', INT(contatore), ' passi.',' I tre punti sono '
				WRITE (7, *) '         punti:   ', a, b ,c
				WRITE (7, *) '         energie: ', fa(1), fb(1), fc(1)
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
				IF (flag_output) WRITE (7, *) 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
			END IF
			b=b+first_step
			rb=r0+b*dir
			b_app=b
			CALL controlla_punto_e_parametro(r0,rb,b,dir)
			flag_border=.FALSE.
			IF (b_app/=b) flag_border=.TRUE.
			rb=r0+b*dir
			CALL calcola_energia(idb,rb,N,fb(1:2),accettabile)
			IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 1[ module_variational_opt.f90 > bracket_min_multidim ]'
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
					IF (flag_output) WRITE (7, *) 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
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
				IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 2[ module_variational_opt.f90 > bracket_min_multidim ]'
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
					IF (flag_output) WRITE (7, *) 'VAR_OPT: bracket - first step era tropo piccolo, l ho raddoppiato'
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
				IF (.NOT. accettabile) STOP 'Problema accettabilitá-risoluzione 3 [ module_variational_opt.f90 > bracket_min_multidim ]'
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
				IF ((id/=ida) .AND. (id/=idb) .AND. (id/=idc)) CALL SYSTEM ('export ID_TO_ERASE=opt'//id//'; find estimatori/ -name *$ID_TO_ERASE* -exec rm -f {} \\; find posizioni/ -name *$ID_TO_ERASE* -exec rm -f {} \\')
			END DO
		END IF

		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Ho eseguito il bracket con ', INT(contatore), ' passi.',' I tre punti sono '
			PRINT * , '         punti:   ', a, b ,c
			PRINT * , '         energie: ', fa(1), fb(1), fc(1)
			IF (flag_output) THEN
				WRITE (7, *) 'VAR_OPT: Ho eseguito il bracket con ', INT(contatore), ' passi.',' I tre punti sono '
				WRITE (7, *) '         punti:   ', a, b ,c
				WRITE (7, *) '         energie: ', fa(1), fb(1), fc(1)
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
					WRITE (7, *) 'VAR_OPT: Parabgold Search (a,b,c) - ', contatore
					WRITE (7, *) '         ', a, b, c
				END IF
			END IF
		END DO
		IF (mpi_myrank==0) THEN
			DO i = 1, 10, 1                    !elimino i dati ormai inutili
				id='p'//CHAR(MOD(i,10)+48)
				IF (id/=idb) CALL SYSTEM ('export ID_TO_ERASE=opt'//id//'; find estimatori/ -name *$ID_TO_ERASE* -exec rm -f {} \\; find posizioni/ -name *$ID_TO_ERASE* -exec rm -f {} \\')
			END DO
			IF (idb/='b') CALL SYSTEM ('export ID_TO_ERASE=optb; find estimatori/ -name *$ID_TO_ERASE* -exec rm -f {} \\; find posizioni/ -name *$ID_TO_ERASE* -exec rm -f {} \\')
		END IF
		IF (mpi_myrank==0) THEN
			PRINT * , 'VAR_OPT: Parabgold Search - sono stati necessari ', contatore, ' passi'
			IF (flag_output) WRITE (7, *) 'VAR_OPT: Parabgold Search - sono stati necessari ', contatore, ' passi'
		END IF

	END SUBROUTINE parabgold_search_multidim

    SUBROUTINE costrizione_rp(nuovi_rp)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(INOUT) :: nuovi_rp(1:num_coord_rp)

        REAL (KIND=8)                :: rp_rshaped(3,num_coord_rp/3)
        REAL (KIND=8), ALLOCATABLE   :: h2coms(:,:), h2vecs(:,:)
        INTEGER                      :: itp, nump, numph

        nump = num_coord_rp/3
        numph = nump/2
        rp_rshaped = RESHAPE(nuovi_rp, (/3,nump/))

        SELECT CASE (costri_rp)
        CASE ('ring__')
            rp_rshaped(3,:) = 0.d0  !prevent movement in z dir
            DO itp=1,nump
                rp_rshaped(1:2,itp) = costri_rp_param * rp_rshaped(1:2, itp) / &
                sqrt(dot_product(rp_rshaped(1:2, itp), rp_rshaped(1:2, itp)))
            ENDDO

        CASE ('h2ring')
            ALLOCATE(h2coms(3, num_coord_rp/6), h2vecs(3, num_coord_rp/6))
            h2coms = 0.5d0 * (rp_rshaped(:,1:numph) + rp_rshaped(:,numph+1:nump))
            h2vecs = rp_rshaped(:, 1:numph) - h2coms
            h2coms(3,:) = 0.d0  !prevent movement in z dir
            DO itp=1,numph
                h2coms(1:2,itp) = costri_rp_param * h2coms(1:2, itp) / &
                sqrt(dot_product(h2coms(1:2, itp), h2coms(1:2, itp)))
                rp_rshaped(:,itp) = h2coms(:, itp) + h2vecs(:, itp)
                rp_rshaped(:,itp+numph) = h2coms(:, itp) - h2vecs(:, itp)
            ENDDO

        CASE ('hprism')
            ALLOCATE(h2coms(3, num_coord_rp/6), h2vecs(3, num_coord_rp/6))
            h2coms = 0.5d0 * (rp_rshaped(:,1:numph) + rp_rshaped(:,numph+1:nump))
            h2vecs = rp_rshaped(:, 1:numph) - h2coms
            h2vecs(1:2,:) = 0.d0 !prevent movement in xy dir for relative hvecs
            h2coms(3,:) = 0.d0  !prevent movement in z dir for com
            DO itp=1,numph
                h2coms(1:2,itp) = costri_rp_param * h2coms(1:2, itp) / &
                sqrt(dot_product(h2coms(1:2, itp), h2coms(1:2, itp)))
                rp_rshaped(:,itp) = h2coms(:, itp) + h2vecs(:, itp)
                rp_rshaped(:,itp+numph) = h2coms(:, itp) - h2vecs(:, itp)
            ENDDO

        END SELECT

        nuovi_rp = RESHAPE(rp_rshaped, (/3*nump/))

    END SUBROUTINE

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

		IF ( opt_rp ) THEN
		    call costrizione_rp(nuovi_parametri(num_par_var-num_coord_rp+1 : num_par_var))
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

	SUBROUTINE estrai_medie_SR()
		USE VMC
		USE variational_calculations
		IMPLICIT NONE
      LOGICAL, PARAMETER :: NORMALIZE=.TRUE.
		INTEGER :: i, j, k, i1, j1
      REAL(KIND=8) :: acca_i(1:num_par_var), eta
		REAL (KIND=8) :: dummy1(1:N_mc), quoz
      REAL(KIND=8), ALLOCATABLE :: O_fast(:,:)
      REAL(KIND=8), ALLOCATABLE :: Hi_fast(:,:)
      LOGICAL :: mask(1:num_par_var)
      REAL(KIND=8) :: a, b, c

      !!!Definisco il quoziente per calcolare le medie
      quoz=1.d0/REAL(N_mc*mpi_nprocs,8)

      !!!Inizializzo le medie
      H=0.d0
      Oi=0.d0
      HOi=0.d0
      OiOj=0.d0
      Hi=0.d0
      OiHj=0.d0
      OiHOj=0.d0

      !!!Metto insieme i flag_O di tutti i processori
      CALL MPI_REDUCE(flag_O,mask,num_par_var,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,mpi_ierr)
      IF (mpi_myrank==0) THEN
         !assicuro che gli spostamenti per i parametri associati ai Rp non siano mascherati
         mask(num_par_var-num_coord_Rp+1:num_par_var)=.TRUE.
         !aggiorno flag_O dai flag raccolti
         flag_O=mask
      END IF
      !!!Li ridistribuisco a tutti i processori
      CALL MPI_BCAST(flag_O,num_par_var,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_ierr)
      !genero la matrice dei parametri usati
      used_par=flag_O
      DO j = 1, num_par_var, 1
      DO i = 1, num_par_var, 1
         M_used_par(i,j)=used_par(i).AND.used_par(j)
      END DO
      END DO

      !!!Determino il numero di parametri variazionali effettivi
      num_par_var_eff=COUNT(used_par(1:num_par_var))
      
      !!!Calcolo H
      dummy1(1)=SUM(E_tot(1:N_mc)*w(1:N_mc))
      CALL MPI_REDUCE(dummy1(1),H,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
      IF (NORMALIZE) H=H*quoz
      CALL MPI_BCAST(H,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)

      IF (OAV_ON_THE_FLY) THEN
         DO i = 1, num_par_var, 1
            IF (flag_O(i)) THEN
               !!!Raccolgo i termini Oi
               CALL MPI_REDUCE(Oi_av(i),Oi(i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
               IF (NORMALIZE) Oi(i)=Oi(i)*quoz
               !!!Raccolgo i termini OiH
               CALL MPI_REDUCE(OiH_av(i),HOi(i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
               IF (NORMALIZE) HOi(i)=HOi(i)*quoz
               !!!Raccolgo i termini OiOj
               DO j = i, num_par_var, 1
                  IF (flag_O(j)) THEN
                     CALL MPI_REDUCE(OiOj_av(j,i),OiOj(j,i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
                     IF (NORMALIZE) OiOj(j,i)=OiOj(j,i)*quoz
                     OiOj(i,j)=OiOj(j,i)
                     CALL MPI_REDUCE(OiHOj_av(j,i),OiHOj(j,i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
                     IF (NORMALIZE) OiHOj(j,i)=OiHOj(j,i)*quoz
                     OiHOj(i,j)=OiHOj(j,i)
                  END IF
               END DO
         END IF
         END DO
         IF (opt_L) THEN
            DO i = 1, 3, 1
               !!!Raccolgo i termini Hi
               CALL MPI_REDUCE(HgradL_av(i),Hi(i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
               IF (NORMALIZE) Hi(i)=Hi(i)*quoz
               DO j = 1, 3, 1
                  IF (flag_O(j)) THEN
                     !!!Raccolgo i termini OiHj
                     CALL MPI_REDUCE(OiHgradL_av(j,i),OiHj(j,i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
                     IF (NORMALIZE) OiHj(j,i)=OiHj(j,i)*quoz
                  END IF
               END DO
            END DO
         END IF
         IF (opt_Rp) THEN
            DO i = 1, num_coord_Rp, 1
               !!!Raccolgo i termini Hi
               CALL MPI_REDUCE(Hgradp_av(i),Hi(i+num_coord_L),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
               IF (NORMALIZE) Hi(i+num_coord_L)=Hi(i+num_coord_L)*quoz
               DO j = 1, num_coord_Rp, 1
                  IF (flag_O(j)) THEN
                     !!!Raccolgo i termini OiHj
                     CALL MPI_REDUCE(OiHgradp_av(j,i),OiHj(j+num_coord_L,i+num_coord_L),&
                        1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
                     IF (NORMALIZE) OiHj(j+num_coord_L,i+num_coord_L)=OiHj(j+num_coord_L,i+num_coord_L)*quoz
                  END IF
               END DO
            END DO
         END IF
      ELSE
         ALLOCATE(O_fast(1:N_mc,1:num_par_var))
         O_fast(1:N_mc,1:num_par_var)=TRANSPOSE(O(1:num_par_var,1:N_mc))
         DO i = 1, num_par_var, 1
            IF (flag_O(i)) THEN
               !!!Calcolo i termini Oi
               dummy1(1)=SUM(O_fast(1:N_mc,i)*w(1:N_mc))
               CALL MPI_REDUCE(dummy1(1),Oi(i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
               IF (NORMALIZE) Oi(i)=Oi(i)*quoz
               !!!Calcolo i termini HOi
               dummy1(1)=SUM(E_tot(1:N_mc)*O_fast(1:N_mc,i)*w(1:N_mc))
               CALL MPI_REDUCE(dummy1(1),HOi(i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
               IF (NORMALIZE) HOi(i)=HOi(i)*quoz
               !!!Calcolo i termini OiOj
               DO j = i, num_par_var, 1
                  IF (flag_O(j)) THEN
                     dummy1(1)=SUM(O_fast(1:N_mc,i)*O_fast(1:N_mc,j)*w(1:N_mc))
                     CALL MPI_REDUCE(dummy1(1),OiOj(j,i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
                     IF (NORMALIZE) OiOj(j,i)=OiOj(j,i)*quoz
                     OiOj(i,j)=OiOj(j,i)
                     dummy1(1)=SUM(O_fast(1:N_mc,i)*E_tot(1:N_mc)*O_fast(1:N_mc,j)*w(1:N_mc))
                     CALL MPI_REDUCE(dummy1(1),OiHOj(j,i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
                     IF (NORMALIZE) OiHOj(j,i)=OiHOj(j,i)*quoz
                     OiHOj(i,j)=OiHOj(j,i)
                  END IF
               END DO
            END IF
         END DO
         IF (opt_L) THEN
            ALLOCATE(Hi_fast(1:N_mc,1:3))
            Hi_fast(1:N_mc,1:3)=TRANSPOSE(HgradL(1:3,1:N_mc))
            DO i = 1, 3, 1
               !!!Calcolo i termini Hi
               dummy1(1)=SUM(Hi_fast(1:N_mc,i)*w(1:N_mc))
               CALL MPI_REDUCE(dummy1(1),Hi(i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
               IF (NORMALIZE) Hi(i)=Hi(i)*quoz
               DO j = 1, 3, 1
                  !!!Calcolo i termini OiHj
                  dummy1(1)=SUM(O_fast(1:N_mc,j)*Hi_fast(1:N_mc,i)*w(1:N_mc))
                  CALL MPI_REDUCE(dummy1(1),OiHj(j,i),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
                  IF (NORMALIZE) OiHj(j,i)=OiHj(j,i)*quoz
               END DO
            END DO
            DEALLOCATE(Hi_fast)
         END IF
         IF (opt_Rp) THEN
            ALLOCATE(Hi_fast(1:N_mc,1:num_coord_Rp))
            Hi_fast(1:N_mc,1:num_coord_Rp)=TRANSPOSE(Hgradp(1:num_coord_Rp,1:N_mc))
            DO i = 1, num_coord_Rp, 1
               !!!Calcolo i termini Hi
               dummy1(1)=SUM(Hi_fast(1:N_mc,i)*w(1:N_mc))
               CALL MPI_REDUCE(dummy1(1),Hi(i+num_coord_L),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
               IF (NORMALIZE) Hi(i+num_coord_L)=Hi(i+num_coord_L)*quoz
               DO j = 1, num_coord_Rp, 1
                  !!!Calcolo i termini OiHj
                  dummy1(1)=SUM(O_fast(1:N_mc,j)*Hi_fast(1:N_mc,i)*w(1:N_mc))
                  CALL MPI_REDUCE(dummy1(1),OiHj(j+num_coord_L,i+num_coord_L),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
                  IF (NORMALIZE) OiHj(j+num_coord_L,i+num_coord_L)=OiHj(j+num_coord_L,i+num_coord_L)*quoz
               END DO
            END DO
            DEALLOCATE(Hi_fast)
         END IF
         DEALLOCATE(O_fast)
      END IF

      !Distribuisco H, Oi, HOi, OiOj a tutti i processori
      CALL MPI_BCAST(Oi,num_par_var,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
      CALL MPI_BCAST(HOi,num_par_var,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
      CALL MPI_BCAST(OiOj,num_par_var*num_par_var,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
      CALL MPI_BCAST(OiHOj,num_par_var*num_par_var,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
      IF (opt_Rp) CALL MPI_BCAST(Hi,num_coord_Rp,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
      IF (opt_Rp) CALL MPI_BCAST(OiHj,num_coord_Rp*num_coord_Rp,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)

      !Alloco s_kn e f_k
      IF (ALLOCATED(s_kn)) DEALLOCATE(s_kn)
      ALLOCATE(s_kn(1:num_par_var_eff,1:num_par_var_eff))
      IF (ALLOCATED(f_k)) DEALLOCATE(f_k)
      ALLOCATE(f_k(1:num_par_var_eff))

      !Alloco e inizializzo gli Oi, HOi e OiOj effettivi
      !Oi
      IF (ALLOCATED(Oi_eff)) DEALLOCATE(Oi_eff)
      ALLOCATE(Oi_eff(1:num_par_var_eff))
      Oi_eff=PACK(Oi,used_par)
      !HOi
      IF (ALLOCATED(HOi_eff)) DEALLOCATE(HOi_Eff)
      ALLOCATE(HOi_eff(1:num_par_var_eff))
      HOi_eff=PACK(HOi,used_par)
      !OiOj
      IF (ALLOCATED(OiOj_eff)) DEALLOCATE(OiOj_eff)
      ALLOCATE(OiOj_eff(1:num_par_var_eff,1:num_par_var_eff))
      OiOj_eff=RESHAPE( PACK(OiOj,M_used_par) , (/ num_par_var_eff,num_par_var_eff /)  )
      !OiHOj
      IF (ALLOCATED(OiHOj_eff)) DEALLOCATE(OiHOj_eff)
      ALLOCATE(OiHOj_eff(1:num_par_var_eff,1:num_par_var_eff))
      OiHOj_eff=RESHAPE( PACK(OiHOj,M_used_par) , (/ num_par_var_eff,num_par_var_eff /)  )

      !IF (mpi_myrank==0) PRINT *, " > > > Oi: ", Oi
      !IF (mpi_myrank==0) PRINT *, " > > > HOi: ", HOi

   END SUBROUTINE estrai_medie_SR


   SUBROUTINE calcola_skn_fk()
      IMPLICIT NONE
      INTEGER :: i, j
      REAL(KIND=8) :: eta
      REAL(KIND=8) :: Hk(1:num_par_var_eff)
      
      !!!Costruisco i termini Hk e eta
      eta=1.d0-f_SR_beta*SR_beta*H
      DO i = 1, num_par_var_eff, 1
         Hk(i)=Oi_eff(i)-f_SR_beta*SR_beta*HOi_eff(i)
      END DO
      !!!Costruisco s_kn e f_k
      DO i = 1, num_par_var_eff, 1
         f_k(i)=Oi_eff(i)*eta-Hk(i)
         !f_k(i)=-f_SR_beta*SR_beta*(Oi_eff(i)*H-HOi_eff(i))  !ORIGINAL SR
         DO j = 1, num_par_var_eff, 1
            s_kn(j,i)=OiOj_eff(j,i)*eta-Oi_eff(i)*Hk(j)
            !s_kn(j,i)=OiOj_eff(j,i)-Oi_eff(j)*Oi_eff(i)  !ORIGINAL SR
         END DO
      END DO

   END SUBROUTINE calcola_skn_fk
   

   !SUBROUTINE plain_SR_trova_dp(dp)
   !   IMPLICIT NONE
   !   INTEGER, PARAMETER :: LWORK=10
   !   REAL(KIND=8), INTENT(OUT) :: dp(1:num_par_var)
   !   REAL(KIND=8) :: dp_eff(1:num_par_var_eff)
   !   INTEGER :: i, j, N
   !   INTEGER :: info
   !   REAL(KIND=8) :: Is_kn(1:num_par_var_eff,1:num_par_var_eff)
   !   REAL(KIND=8) :: work(1:LWORK*num_par_var_eff)
   !   REAL(KIND=8) :: Usvd(1:num_par_var_eff,1:num_par_var_eff)
   !   REAL(KIND=8) :: VTsvd(1:num_par_var_eff,1:num_par_var_eff)
   !   REAL(KIND=8) :: Ssvd(1:num_par_var_eff)

   !   STOP "MUST BE CHECKED AFTER INTRODUCING OPTIMIZATION OF L"

   !   N=num_par_var_eff

   !   Is_kn=s_kn
   !   CALL DGESVD('A','A',N,N,Is_kn,N,Ssvd,Usvd,N,VTsvd,N,work,LWORK*N,info)
   !   IF (info /= 0) THEN
   !      PRINT *, "Error in SVD [ stochastic_reconfiguration > module_variational_opt.f90 ]"
	!	 PRINT*, "Info=", info
   !      STOP
   !   END IF
   !   Is_kn=0.d0
   !   DO i = 1, N, 1
   !      IF (Ssvd(i)>Ssvd(1)*SVD_MIN) THEN
   !         Is_kn(i,i)=1.d0/Ssvd(i)
   !      ELSE
   !         Is_kn(i,i)=0.d0
   !      END IF
   !   END DO
   !   Is_kn=MATMUL(Usvd,MATMUL(Is_kn,VTsvd))

	!	dp_eff=0.d0
	!	DO j = 1, N, 1
	!		DO i = 1, N, 1
	!			dp_eff(j)=dp_eff(j)-Is_kn(j,i)*f_k(i)
	!		END DO
	!	END DO

   !   dp=0.d0
   !   dp=UNPACK(dp_eff,used_par,dp)
   !   
   !END SUBROUTINE plain_SR_trova_dp

   !use le equazioni derivanti dalla lagrangiana scritta da Markus
   SUBROUTINE lagrangian_dynamic_forces(N,NL,Np,dp)
      USE dati_fisici
      IMPLICIT NONE
      INTEGER, PARAMETER :: LWORK=10
      INTEGER, INTENT(IN) :: N, NL, Np !number of effective variational parameters, electronic, and protonic
      REAL(KIND=8), INTENT(OUT) :: dp(1:num_par_var)
      REAL(KIND=8) :: dp_eff(1:N)
      INTEGER :: Ne
      INTEGER :: i, j
      INTEGER :: info
      REAL(KIND=8), ALLOCATABLE :: M(:,:), v(:)
      REAL(KIND=8), ALLOCATABLE :: IM(:,:)
      REAL(KIND=8), ALLOCATABLE :: work(:)
      REAL(KIND=8), ALLOCATABLE :: Usvd(:,:)
      REAL(KIND=8), ALLOCATABLE :: VTsvd(:,:)
      REAL(KIND=8), ALLOCATABLE :: Ssvd(:)

      dp_eff=0.d0

      Ne=N-NL-Np
      ALLOCATE( M(1:Ne,1:Ne), v(1:Ne))
      ALLOCATE(IM(1:Ne,1:Ne))
      ALLOCATE(work(1:LWORK*(Ne)))
      ALLOCATE(Usvd(1:Ne,1:Ne))
      ALLOCATE(VTsvd(1:Ne,1:Ne))
      ALLOCATE(Ssvd(1:Ne))

      IF (Ne>0) THEN
         !costruisco la matrice per la parte elettronica
         M(1:Ne,1:Ne)=OiOj_eff(NL+1:NL+Ne,NL+1:NL+Ne)
         v(1:Ne)=HOi_eff(NL+1:NL+Ne)-H*Oi_eff(NL+1:NL+Ne)

         !inverto M
         IM=M
         CALL DGESVD('A','A',Ne,Ne,IM,Ne,Ssvd,Usvd,Ne,VTsvd,Ne,work,LWORK*Ne,info)
         IF (info /= 0) THEN
            PRINT *, "Error in SVD [ lagrangian_dynamic_forces > module_variational_opt.f90 ]"
		    PRINT*, "Info=", info
            STOP
         END IF
         IM=0.d0
         DO i = 1, Ne, 1
            IF (Ssvd(i)>Ssvd(1)*SVD_MIN) THEN
               IM(i,i)=1.d0/Ssvd(i)
            ELSE
               IM(i,i)=0.d0
            END IF
         END DO
         IM=MATMUL(Usvd,MATMUL(IM,VTsvd))

         !assegno i nuovi dp elettronici
         dp_eff(NL+1:NL+Ne)=MATMUL(IM(1:Ne,1:Ne),-v(1:Ne))
      END IF

      IF (NL>0) THEN
         !assegno i nuovi dp protonici
         dp_eff(1:3)=-Hi(1:3)-2.d0*HOi_eff(1:3)+2.d0*H*Oi_eff(1:3)
         !if a 2D structure is used, then the z component is set to zero
         IF (flag_2D) dp_eff(3)=0.d0
      END IF

      IF (Np>0) THEN
         !assegno i nuovi dp protonici
         dp_eff(NL+Ne+1:N)=-Hi(NL+1:NL+Np)-2.d0*HOi_eff(NL+Ne+1:N)+2.d0*H*Oi_eff(NL+Ne+1:N)
         !if a 2D structure is used, then the z component is set to zero
         IF (flag_2D) THEN
            DO i = 1, N_part, 1
               dp_eff(Ne+3+(i-1)*3)=0.d0
            END DO
         END IF
      END IF

      !!riporto i dp_eff al dp globale
      dp=0.d0
      dp=UNPACK(dp_eff,used_par,dp)

      DEALLOCATE(M, v)
      DEALLOCATE(IM)
      DEALLOCATE(work)
      DEALLOCATE(Usvd)
      DEALLOCATE(VTsvd)
      DEALLOCATE(Ssvd)
      
   END SUBROUTINE lagrangian_dynamic_forces

   !!risolve le quazioni risultati dalla minimizzazione di H'-H
   !SUBROUTINE gradH_trova_dp(N,Ne,Np,dp)
   !   IMPLICIT NONE
   !   INTEGER, PARAMETER :: LWORK=10
   !   INTEGER, INTENT(IN) :: N, Ne, Np !number of effective variational parameters, electronic, and protonic
   !   REAL(KIND=8), INTENT(OUT) :: dp(1:num_par_var)
   !   REAL(KIND=8) :: dp_eff(1:N)
   !   INTEGER :: i, j, k
   !   INTEGER :: info
   !   REAL(KIND=8) :: Hi_gen(1:N), OiHj_gen(1:N,1:N)
   !   REAL(KIND=8) :: M(1:N,1:N), v(1:N)
   !   REAL(KIND=8) :: IM(1:N,1:N)
   !   REAL(KIND=8) :: work(1:LWORK*(N))
   !   REAL(KIND=8) :: Usvd(1:N,1:N)
   !   REAL(KIND=8) :: VTsvd(1:N,1:N)
   !   REAL(KIND=8) :: Ssvd(1:N)

   !   STOP "MUST BE CHECKED AFTER INTRODUCING OPTIMIZATION OF L"

   !   !costruisco il vettore generalizzato Hi
   !   Hi_gen=0.d0
   !   IF (Np>0) Hi_gen(Ne+1:N)=Hi(1:Np)

   !   !costruisco la matrice generalizzata HiOj
   !   OiHj_gen=0.d0
   !   IF (Np>0) OiHj_gen(1:N,Ne+1:N)=OiHj_eff(1:N,1:Np)

   !   !costruisco la matrice M
   !   M=0.d0
   !   DO j = 1, N, 1
   !   DO i = 1, N, 1
   !      !M(i,j)=4.d0*HOi_eff(j)*Oi_eff(i)-4.d0*HOi_eff(i)*Oi_eff(j)-2.d0*H*OiOj_eff(i,j)
   !      M(i,j) = 2.d0*Oi_eff(i)*(Hi_gen(j)+2.d0*HOi_eff(j)-2.d0*H*Oi_eff(j)) &
   !               - 2.d0*Oi_eff(j)*(Hi_gen(i)+2.d0*HOi_eff(i)-2.d0*H*Oi_eff(i)) &
   !               - 2.d0*H*OiOj_eff(i,j) & !manca il termine OiHOj
   !               + 2.d0*(OiHj_gen(i,j)+OiHj_gen(j,i))
   !   END DO
   !   END DO

   !   !costruisco il vettore v
   !   v=0.d0
   !   v(1:N) = 2.d0*(HOi_eff(1:N)-H*Oi_eff(1:N)) &
   !            + Hi_gen(1:N)

   !   !inverto M
   !   IM=M
   !   CALL DGESVD('A','A',N,N,IM,N,Ssvd,Usvd,N,VTsvd,N,work,LWORK*N,info)
   !   IF (info /= 0) THEN
   !      PRINT *, "Error in SVD [ gradH_trova_dp > module_variational_opt.f90 ]"
	!	 PRINT*, "Info=", info
   !      STOP
   !   END IF
   !   IM=0.d0
   !   DO i = 1, N, 1
   !      IF (Ssvd(i)>Ssvd(1)*SVD_MIN) THEN
   !         IM(i,i)=1.d0/Ssvd(i)
   !      ELSE
   !         IM(i,i)=0.d0
   !      END IF
   !   END DO
   !   IM=MATMUL(Usvd,MATMUL(IM,VTsvd))

   !   !assegno i nuovi dp elettronici
   !   dp_eff(1:N)=MATMUL(IM(1:N,1:N),-v(1:N))

   !   !riporto i dp_eff al dp globale
   !   dp=0.d0
   !   dp=UNPACK(dp_eff,used_par,dp)
   !   
   !END SUBROUTINE gradH_trova_dp

   !!risolve l'equazione dei moltiplicatori di lagrange linearizzata
   !SUBROUTINE lin_molt_lagr(N,Ne,Np,dp)
   !   IMPLICIT NONE
   !   INTEGER, PARAMETER :: LWORK=10
   !   INTEGER, INTENT(IN) :: N, Ne, Np !number of effective variational parameters, electronic, and protonic
   !   REAL(KIND=8), INTENT(OUT) :: dp(1:num_par_var)
   !   REAL(KIND=8) :: dp_eff(1:N)
   !   INTEGER :: i, j, k
   !   INTEGER :: info
   !   REAL(KIND=8) :: M(0:N,0:N), v(0:N)
   !   REAL(KIND=8) :: IM(1:N+1,0:N+1)
   !   REAL(KIND=8) :: work(1:LWORK*(N+1))
   !   REAL(KIND=8) :: Usvd(0:N,0:N)
   !   REAL(KIND=8) :: VTsvd(0:N,0:N)
   !   REAL(KIND=8) :: Ssvd(0:N)
   !   INTEGER :: pvt(1:N+1)
   !   REAL(KIND=8) :: soluzione(0:N)

   !   STOP "MUST BE CHECKED AFTER INTRODUCING OPTIMIZATION OF L"

   !   !costruisco la matrice M
   !   M(0,0)=0.d0
   !   M(0,1:N)=Oi_eff(1:N)
   !   M(1:N,0)=Oi_eff(1:N)
   !   DO j = 1, N, 1
   !   DO i = 1, N, 1
   !      M(i,j)=(HOi_eff(j)-H*Oi_eff(j))*Oi_eff(i)-(HOi_eff(i)-H*Oi_eff(i))*Oi_eff(j)
   !   END DO
   !   END DO

   !   !CALL RANDOM_NUMBER(M)
   !   !DO i = 0, N, 1
   !   !   M(i,i)=M(i,i)+5.d0
   !   !END DO
   !   

   !   !costruisco il vettore v
   !   v(0)=0.d0
   !   v(1:N)=HOi_eff(1:N)-H*Oi_eff(1:N)

   !   !!!inverto M
   !   IM(0:N,0:N)=M(0:N,0:N)

   !   !!Uso il metodo singular value decomposition
   !   !CALL DGESVD('A','A',N+1,N+1,IM(0:N,0:N),N+1,Ssvd(0:N),Usvd(0:N,0:N),N+1,&
   !   !   VTsvd(0:N,0:N),N+1,work(1:LWORK*(N+1)),LWORK*(N+1),info)
   !   !IF (info /= 0) THEN
   !   !   PRINT *, "Error in SVD [ gradH_trova_dp > module_variational_opt.f90 ]"
	!	! PRINT*, "Info=", info
   !   !   STOP
   !   !END IF
   !   !IM=0.d0
   !   !DO i = 0, N, 1
   !   !   IF (Ssvd(i)>Ssvd(0)*SVD_MIN) THEN
   !   !      IM(i,i)=1.d0/Ssvd(i)
   !   !   ELSE
   !   !      IM(i,i)=0.d0
   !   !   END IF
   !   !END DO
   !   !IM(0:N,0:N)=MATMUL(Usvd(0:N,0:N),MATMUL(IM(0:N,0:N),VTsvd(0:N,0:N)))
   !   
   !   !Diagonalizzazione diretta
   !   pvt=0
   !   work=0.d0
   !   CALL DGETRF( N+1, N+1, IM(0:N,0:N), N+1, pvt, info )
   !   IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
   !   CALL DGETRI( N+1, IM(0:n,0:N), N+1, pvt, work, LWORK*(N+1), info )
   !   IF (info/=0) STOP 'ERRORE NELLA INVERSIONE'

   !   DO i = 0, N, 1
   !      PRINT *, MATMUL(M(i,0:N),IM(0:N,0:N))
   !   END DO

   !   !assegno i nuovi dp elettronici
   !   soluzione(0:N)=MATMUL(-v(0:N),IM(0:N,0:N))
   !   dp_eff(1:N)=soluzione(1:N)


   !   IF (mpi_myrank==0) THEN
   !      PRINT *, "lambda = ", soluzione(0)
   !      PRINT *, "dp_eff = ", dp_eff
   !      PRINT *, MATMUL(soluzione(0:N),M(0:N,0:N))+v(0:N)
   !      STOP
   !   END IF

   !   !riporto i dp_eff al dp globale
   !   dp=0.d0
   !   dp=UNPACK(dp_eff,used_par,dp)
   !   
   !END SUBROUTINE lin_molt_lagr

   !SUBROUTINE moltiplicatori_lagrange_trova_dp(N,Ne,Np,dp)
   !   IMPLICIT NONE
   !   INTEGER, PARAMETER :: LWORK=10
   !   INTEGER, INTENT(IN) :: N, Ne, Np !number of effective variational parameters, electronic, and protonic
   !   REAL(KIND=8), INTENT(OUT) :: dp(1:num_par_var)
   !   INTEGER :: i1
   !   REAL(KIND=8) :: dp_eff(1:N), lambda
   !   REAL(KIND=8) :: p(0:N)
   !   REAL(KIND=8) :: EqMoltLagr(0:N), MoltLagrTarget, GradEqMoltLagr(0:N,0:N)
   !   REAL(KIND=8) :: norm_EqMoltLagr, weight_EqMoltLagr(0:N), sign_EqMoltLagr(0:N)
   !   LOGICAL :: flag_loop1, flag_loop2
   !   INTEGER(KIND=8) :: cont_loop1, cont_loop2
   !   REAL(KIND=8) :: dir(0:N), dt1, dt2, dt3, f1, f2, f3, sensibility
   !   REAL(KIND=8) :: EqMoltLagr_new(0:N), MoltLagrTarget_new, p_new(0:N)
   !   INTEGER :: imax

   !   STOP "MUST BE CHECKED AFTER INTRODUCING OPTIMIZATION OF L"

   !   sensibility=0.0000001d0

   !   !trovo dei dp ragionevoli
   !   CALL gradH_trova_dp(N,Ne,Np,dp)
   !   dp_eff=PACK(dp,used_par)
   !   !CALL RANDOM_NUMBER(dp_eff)

   !   !inizializzo lambda con un valore ragionevole
   !   lambda=0.d0
   !   CALL RANDOM_NUMBER(lambda)

   !   !inizializzo p
   !   p(0)=lambda
   !   p(1:N)=dp_eff(1:N)
   !   CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
   !   PRINT *, mpi_myrank, "p(0:N) = ", p(0:N)
   !   CALL eq_molt_lagr(N,Ne,Np,p(0:N),EqMoltLagr(0:N))
   !   MoltLagrTarget=MAXVAL(DABS(EqMoltLagr(0:N)))
   !   CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
   !   PRINT *, mpi_myrank, "Target = ", MoltLagrTarget

   !   flag_loop1=.TRUE.
   !   cont_loop1=0
   !   DO WHILE (flag_loop1)
   !      cont_loop1=cont_loop1+1

   !      !calcolo i moliplicatori di Lagrange
   !      CALL eq_molt_lagr(N,Ne,Np,p(0:N),EqMoltLagr(0:N))
   !      norm_EqMoltLagr=SUM(DABS(EqMoltLagr(0:N)))
   !      weight_EqMoltLagr(0:N)=DABS(EqMoltLagr(0:N))/norm_EqMoltLagr
   !      sign_EqMoltLagr(0:N)=-DSIGN(1.d0,EqMoltLagr(0:N))
   !      MoltLagrTarget=MAXVAL(DABS(EqMoltLagr(0:N)))
   !      imax=MAXLOC(DABS(EqMoltLagr(0:N)),1)-1
   !      !PRINT *, "MoltLagr = ", EqMoltLagr
   !      !PRINT *, "Target = ", MoltLagrTarget

   !      !!!calcolo il gradiente dei moltiplicatori di lagrange
   !      CALL grad_eq_molt_lagr(N,Ne,Np,p(0:N),GradEqMoltLagr(0:N,0:N))
   !      !faccio combinazione lineare di tutti i gradienti
   !      IF (mpi_myrank==0) THEN
   !         dir=0.d0
   !         DO i1 = 0, N, 1 !loop sulle equazioni
   !            dir(0:N)=dir(0:N)+GradEqMoltLagr(i1,0:N)*weight_EqMoltLagr(i1)*sign_EqMoltLagr(i1)
   !         END DO
   !      !uso semplicemente il gradiente dell'equazione piu' diversa da zero
   !      ELSE
   !         dir=GradEqMoltLagr(imax,0:N)*weight_EqMoltLagr(imax)*sign_EqMoltLagr(imax)
   !      END IF

   !      !trovo di quanto spostare nella direzione data
   !      dt1=0.d0
   !      dt2=0.1d0
   !      dt3=10.d0
   !      p_new(0:N)=p(0:N)+dt1*dir(0:N)
   !      CALL eq_molt_lagr(N,Ne,Np,p_new(0:N),EqMoltLagr_new(0:N))
   !      f1=MAXVAL(DABS(EqMoltLagr_new(0:N)))
   !      p_new(0:N)=p(0:N)+dt2*dir(0:N)
   !      CALL eq_molt_lagr(N,Ne,Np,p_new(0:N),EqMoltLagr_new(0:N))
   !      f2=MAXVAL(DABS(EqMoltLagr_new(0:N)))
   !      p_new(0:N)=p(0:N)+dt3*dir(0:N)
   !      CALL eq_molt_lagr(N,Ne,Np,p_new(0:N),EqMoltLagr_new(0:N))
   !      f3=MAXVAL(DABS(EqMoltLagr_new(0:N)))
   !      flag_loop2=.TRUE.
   !      cont_loop2=0
   !      DO WHILE (flag_loop2)
   !         cont_loop2=cont_loop2+1
   !         
   !         IF ((f3>f2).AND.(f2>f1)) THEN
   !            !funzione sta crescendo
   !            dt3=dt2
   !            f3=f2
   !            dt2=0.5d0*dt2
   !            p_new(0:N)=p(0:N)+dt2*dir(0:N)
   !            CALL eq_molt_lagr(N,Ne,Np,p_new(0:N),EqMoltLagr_new(0:N))
   !            f2=MAXVAL(DABS(EqMoltLagr_new(0:N)))
   !         END IF

   !         IF ((f2<f1).AND.(f3<f2)) THEN
   !            !funzione sta decrescendo
   !            dt2=dt3
   !            f2=f3
   !            dt1=dt2
   !            f1=f2
   !            dt3=dt3*2.d0
   !            p_new(0:N)=p(0:N)+dt3*dir(0:N)
   !            CALL eq_molt_lagr(N,Ne,Np,p_new(0:N),EqMoltLagr_new(0:N))
   !            f3=MAXVAL(DABS(EqMoltLagr_new(0:N)))
   !         END IF

   !         IF ((f2>f1).AND.(f2>f3)) THEN
   !            !funzione ha un massimo locale
   !            dt3=dt2
   !            f3=f2
   !            dt2=0.5d0*dt2
   !            p_new(0:N)=p(0:N)+dt2*dir(0:N)
   !            CALL eq_molt_lagr(N,Ne,Np,p_new(0:N),EqMoltLagr_new(0:N))
   !            f2=MAXVAL(DABS(EqMoltLagr_new(0:N)))
   !         END IF

   !         IF ((f2<f1).AND.(f2<f3)) THEN
   !            !funzione ha un minimo locale
   !            IF (f1<f3) THEN
   !               !lato sx minore
   !               dt3=dt2+0.5d0*(dt3-dt2)
   !               p_new(0:N)=p(0:N)+dt3*dir(0:N)
   !               CALL eq_molt_lagr(N,Ne,Np,p_new(0:N),EqMoltLagr_new(0:N))
   !               f3=MAXVAL(DABS(EqMoltLagr_new(0:N)))
   !            ELSE
   !               !lato dx minore
   !               dt1=dt2-0.5d0*(dt2-dt1)
   !               p_new(0:N)=p(0:N)+dt1*dir(0:N)
   !               CALL eq_molt_lagr(N,Ne,Np,p_new(0:N),EqMoltLagr_new(0:N))
   !               f3=MAXVAL(DABS(EqMoltLagr_new(0:N)))
   !            END IF
   !         END IF

   !         IF ((MIN(dt3-dt2,dt2-dt1)<sensibility).AND.((f2<=f1).AND.(f2<=f3))) flag_loop2=.FALSE.
   !         IF (cont_loop2>10000) THEN
   !            !non sono riuscito a trovare il dt
   !            !PRINT *, "Non riesco a trovare il dt corretto"
   !            !PRINT *, "f1, f2, f3 = ", f1, f2 ,f3
   !            !PRINT *, "dt1, dt2, dt3 = ", dt1, dt2, dt3
   !            flag_loop2=.FALSE.
   !            flag_loop1=.FALSE. 
   !         END IF

   !      END DO

   !      !PRINT *, "dt = ", dt2

   !      !aggiorno i parametri
   !      p(0:N)=p(0:N)+dt2*dir(0:N)

   !      IF (cont_loop1>10000) flag_loop1=.FALSE.

   !   END DO

   !   CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
   !   PRINT *, mpi_myrank, "p(0:N) = ", p(0:N)
   !   CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
   !   PRINT *, mpi_myrank, "Target = ", MoltLagrTarget
   !   
   !   dp_eff(1:N)=p(1:N)

   !   !riporto i dp_eff al dp globale
   !   dp=0.d0
   !   dp=UNPACK(dp_eff,used_par,dp)
   !   
   !END SUBROUTINE moltiplicatori_lagrange_trova_dp

   !SUBROUTINE eq_molt_lagr(N,Ne,Np,p,EqMoltLagr)
   !   IMPLICIT NONE
   !   INTEGER, INTENT(IN) :: N, Ne, Np
   !   REAL(KIND=8), INTENT(IN) :: p(0:N)
   !   REAL(KIND=8), INTENT(OUT) :: EqMoltLagr(0:N)
   !   INTEGER :: i1, i2, i3
   !   REAL(KIND=8) :: rho, drho, dHOi(1:N), vi(1:N), mij(1:N,1:N), tijk(1:N,1:N,1:N)
   !   REAL(KIND=8) :: f1(1:N), f2, tc, vtc(1:N)

   !   STOP "MUST BE CHECKED AFTER INTRODUCING OPTIMIZATION OF L"

   !   !precalculations
   !   rho=DOT_PRODUCT(p(1:N),Oi_eff(1:N))
   !   !
   !   !drho=2.d0*rho+DOT_PRODUCT(p(1:N),MATMUL(OiOj_eff(1:N,1:N),p(1:N)))
   !   drho=2.d0*rho
   !   !
   !   rho=1.d0+rho
   !   dHOi(1:N)=HOi_eff(1:N)-H*Oi_eff(1:N)
   !   vi(1:N)=dHOi(1:N)+2.d0*p(0)*Oi_eff(1:N)
   !   DO i1 = 1, N, 1
   !      mij(1:N,i1)=p(0)*( OiOj_eff(1:N,i1)+2.d0*Oi_eff(1:N)*Oi_eff(i1) )
   !      tijk(1:N,1:N,i1)=p(0)*OiOj_eff(1:N,1:N)*Oi_eff(i1)
   !   END DO
   !   tc=0.d0
   !   DO i1 = 1, N, 1
   !   DO i2 = 1, N, 1
   !   DO i3 = 1, N, 1
   !      tc=tc+p(i1)*p(i2)*p(i3)*tijk(i3,i2,i1)
   !   END DO
   !   END DO
   !   END DO
   !   vtc=0.d0
   !   DO i1 = 1, N, 1
   !   DO i2 = 1, N, 1
   !      vtc(1:N)=vtc(1:N)+p(i1)*p(i2)*(2.d0*tijk(1:N,i2,i1)+tijk(i2,i1,1:N))
   !   END DO
   !   END DO
   !   !
   !   f1(1:N) = vi(1:N) + 2.d0*MATMUL(p(1:N),mij(1:N,1:N)) + vtc(1:N)
   !   !f1(1:N) = vi(1:N)! + 2.d0*MATMUL(p(1:N),mij(1:N,1:N)) + vtc(1:N)
   !   !
   !   !
   !   f2 = DOT_PRODUCT(p(1:N),vi(1:N)) + DOT_PRODUCT(p(1:N),MATMUL(p(1:N),mij(1:N,1:N))) + tc
   !   !f2 = DOT_PRODUCT(p(1:N),vi(1:N)) ! + DOT_PRODUCT(p(1:N),MATMUL(p(1:N),mij(1:N,1:N))) + tc
   !   !
   !   !calcola il gradiente della lagrangiana vincolata
   !   !
   !   EqMoltLagr(0)=drho/rho
   !   !EqMoltLagr(0)=drho
   !   !
   !   !
   !   EqMoltLagr(1:N)= ( f1(1:N) - Oi_eff(1:N)*f2/rho )/rho
   !   !EqMoltLagr(1:N)= f1(1:N) - Oi_eff(1:N)*f2/rho
   !   !

   !END SUBROUTINE eq_molt_lagr

   !SUBROUTINE grad_eq_molt_lagr(N,Ne,Np,p,dir)
   !   IMPLICIT NONE
   !   INTEGER, INTENT(IN) :: N, Ne, Np
   !   REAL(KIND=8), INTENT(IN) :: p(0:N)
   !   REAL(KIND=8), INTENT(OUT) :: dir(0:N,0:N)  !primo indice: equazioni, secondo indice: componenti p della direzione
   !   INTEGER :: i1
   !   REAL(KIND=8), PARAMETER :: dx=0.000001d0
   !   REAL(KIND=8) :: Idx, pnew(0:N), f0(0:N), f(0:N,0:N)

   !   CALL eq_molt_lagr(N,Ne,Np,p(0:N),f0(0:N))
   !   pnew(0:N)=p(0:N)
   !   DO i1 = 0, N, 1
   !      pnew(i1)=p(i1)+dx
   !      CALL eq_molt_lagr(N,Ne,Np,pnew(0:N),f(0:N,i1))
   !      pnew(i1)=p(i1)
   !   END DO

   !   Idx=1.d0/dx
   !   DO i1 = 0, N, 1
   !      dir(0:N,i1)=(f(0:N,i1)-f0(0:N))*Idx
   !   END DO
   !   
   !END SUBROUTINE grad_eq_molt_lagr

   !SUBROUTINE lagrangian2_gradient(N,Ne,Np,l,p,GradLagr)
   !   IMPLICIT NONE
   !   REAL(KIND=8), INTENT(OUT) :: GradLagr(0:N) !componente zero per lambda
   !   INTEGER, INTENT(IN) :: N, Ne, Np
   !   REAL(KIND=8), INTENT(IN) :: l, p(1:N) !l=lambda, p=dp
   !   REAL(KIND=8) :: dPsi, Psi, dHO(1:N), dHOO(1:N)

   !   dPsi=2.d0*DOT_PRODUCT(p(1:N),Oi_eff(1:N))+DOT_PRODUCT(p(1:N),MATMUL(OiOj_eff(1:N,1:N),p(1:N)))
   !   Psi=1.d0+dPsi
   !   dHO(1:N)=HOi_eff(1:N)-H*Oi_eff(1:N)
   !   dHOO(1:N)=H*MATMUL(p(1:N),OiOj_eff(1:N,1:N))

   !   GradLagr(0)=dPsi
   !   GradLagr(1:N) = (2.d0*dHO(1:N)-2.d0*dHOO(1:N))/(Psi) &
   !                   - 2.d0*(Oi_eff(1:N)+MATMUL(p(1:N),OiOj_eff(1:N,1:N))) &
   !                     *(2.d0*DOT_PRODUCT(p(1:N),dHO(1:N))-DOT_PRODUCT(p(1:N),dHOO(1:N)))/(Psi*Psi) &
   !                   + l*2.d0*(Oi_eff(1:N)+MATMUL(p(1:N),OiOj_eff(1:N,1:N)))
   !
   !END SUBROUTINE lagrangian2_gradient

   !SUBROUTINE lagrangian_gradient(N,Ne,Np,l,p,GradLagr)
   !   IMPLICIT NONE
   !   REAL(KIND=8), INTENT(OUT) :: GradLagr(0:N) !componente zero per lambda
   !   INTEGER, INTENT(IN) :: N, Ne, Np
   !   REAL(KIND=8), INTENT(IN) :: l, p(1:N) !l=lambda, p=dp
   !   INTEGER :: i1, i2
   !   REAL(KIND=8) :: dPsi, Psi, HdPsi, H_dPsi, Oi_Oj(1:N,1:N)

   !   dPsi=DOT_PRODUCT(p(1:N),Oi_eff(1:N))
   !   Psi=1.d0+dPsi
   !   HdPsi=DOT_PRODUCT(p(1:N),HOi_eff(1:N))
   !   H_dPsi=H*dPsi
   !   DO i2 = 1, N, 1
   !   DO i1 = 1, N, 1
   !      Oi_Oj(i1,i2)=Oi_eff(i1)*Oi_eff(i2)
   !   END DO
   !   END DO

   !   GradLagr(0)=( dPsi+DOT_PRODUCT(p(1:N),MATMUL(Oi_Oj(1:N,1:N),p(1:N))) )&
   !                /dPsi
   !   GradLagr(1:N)=  ( (HOi_eff(1:N)-H*Oi_eff(1:N)+l*Oi_eff(1:N) &
   !                       +2.d0*l*MATMUL(p(1:N),Oi_Oj(1:N,1:N))) &
   !                     /(Psi) ) &
   !                 - ( Oi_eff(1:N)*(HdPsi-H_dPsi+l*dPsi &
   !                      +l*DOT_PRODUCT(p(1:N),MATMUL(Oi_Oj(1:N,1:N),p(1:N)))) &
   !                     /(Psi*Psi) ) 
   !
   !END SUBROUTINE lagrangian_gradient
   !
   !SUBROUTINE trova_dp_Rp_con_forze(dp)
   !   IMPLICIT NONE
   !   REAL(KIND=8) :: dp(1:num_coord_Rp)

   !   dp(1:num_coord_Rp)=-Hi(1:num_coord_Rp)
   !   
   !END SUBROUTINE trova_dp_Rp_con_forze

	SUBROUTINE chiudi_variational_opt()
      USE variational_calculations
		IMPLICIT NONE
		DEALLOCATE(parametri_var)
      DEALLOCATE(used_par,M_used_par)
      DEALLOCATE(Oi,HOi,OiOj)
      DEALLOCATE(OiHOj)
      DEALLOCATE(Hi,OiHj)

      IF (ALLOCATED(s_kn)) DEALLOCATE(s_kn)
      IF (ALLOCATED(f_k)) DEALLOCATE(f_k)
      IF (ALLOCATED(Oi_eff)) DEALLOCATE(Oi_eff)
      IF (ALLOCATED(HOi_eff)) DEALLOCATE(HOi_Eff)
      IF (ALLOCATED(OiOj_eff)) DEALLOCATE(OiOj_eff)
	END SUBROUTINE chiudi_variational_opt
	
END MODULE variational_opt
