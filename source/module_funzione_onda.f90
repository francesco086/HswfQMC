MODULE funzione_onda
   USE markuspline
	IMPLICIT NONE
	REAL, PARAMETER :: CUT_LDA=0.00005  !CUT_LDA=0.000001
	LOGICAL, SAVE :: iniz_funzione_onda=.FALSE.
	INTEGER, SAVE :: num_chiamata_twist_lda
	COMPLEX (KIND=8), SAVE, ALLOCATABLE :: c_eff_dnfH(:,:)
	REAL (KIND=8), SAVE :: Aee_yuk, Aee_ud_yuk, Fee_yuk, Fee_ud_yuk                      !per gli pseudopotenziali di Yukawa
   TYPE(MSPLINE) :: foo_spl
   TYPE(MSPLINE) :: Jsplee, Jsplee_ud, Jsplep, Jsplep_ud
   INTEGER :: m_Jsplee, nknots_Jsplee, m_Jsplep, nknots_Jsplep
   LOGICAL :: cutoff_Jsplee, cutoff_Jsplep
   TYPE(MSPLINE) :: Bsplep
   INTEGER :: m_Bsplep, nknots_Bsplep
   LOGICAL :: cutoff_Bsplep
	REAL (KIND=8), SAVE :: Aep_yuk, Aep_ud_yuk, Fep_yuk, Fep_ud_yuk
	REAL (KIND=8), SAVE :: C_kern_e                                             !per il kernel
	REAL (KIND=8), SAVE :: alpha0_kern_e, alpha1_kern_e, beta0_kern_e, beta1_kern_e
	REAL (KIND=8), SAVE :: l1_kern_e, l2_kern_e, eps_kern_e
	INTEGER, SAVE :: n_kern_e
	LOGICAL, PROTECTED, SAVE :: flag_traccia_coppie
	INTEGER , SAVE :: N_ritraccia_coppie, N_mc_relax_traccia_coppie
	REAL (KIND=8), SAVE :: A_POT_se, D_POT_se                                   !per il Jastrow V_eff per le Shadow
	REAL (KIND=8), SAVE :: c_se													!per il caso ppb
	REAL (KIND=8), SAVE :: Asese_yuk, Asese_ud_yuk, Fsese_yuk, Fsese_ud_yuk  !per il Jastrow di Yukawa per le Shadow
	REAL (KIND=8), SAVE :: B_se, D_se								!per il Jastrow per le Shadow nel caso si voglia includere il bounding fra coppie
	REAL (KIND=8), SAVE :: c_sesp                               			    !per il Jastrow tra shadow elettroniche e protoniche
	REAL (KIND=8), SAVE :: Asesp_yuk, Asesp_ud_yuk, Fsesp_yuk, Fsesp_ud_yuk     !per il Jastrow di Yukawa per le Shadow
	REAL (KIND=8), SAVE :: Gsesp												!per il "Jastrow" gaussiano. In questo caso S_p==R_p
	REAL (KIND=8), SAVE :: Gswf, C_atm
	LOGICAL, PARAMETER, PRIVATE :: verbose_mode=.FALSE.
	LOGICAL :: split_Aee, split_Aep, split_Asese, split_Asesp, split_Fee, split_Fep, split_Fsese, split_Fsesp
	LOGICAL :: flag_usa_coeff_dnfH
	INTEGER, ALLOCATABLE :: coppie(:,:,:), k_pw_int_lda(:,:)
	INTEGER :: N_pw_lda, num_K_points
	CHARACTER(LEN=3), PROTECTED, SAVE :: SDe_kind, Jee_kind, Jep_kind, Jpp_kind, SDse_kind, Jse_kind, Kse_kind, Jsesp_kind
	CHARACTER(LEN=120) :: lda_path
	REAL (KIND=8), ALLOCATABLE :: fattori_orb_lda(:), k_pw_lda(:,:), twist_lda(:), pesi_K_points(:)
	COMPLEX (KIND=8), ALLOCATABLE :: fattori_pw_lda(:,:)						!phi_i(x)= fatt_orb_i * ( sum_k [fatt_pw_i_k*EXP(ikx)] )
	REAL (KIND=8), PROTECTED, SAVE :: kf_coeff_dnfH
	INTEGER, ALLOCATABLE, SAVE :: num_pw_orbit(:), num_pw_dnfH(:), indice_pw_dnfH(:,:)
	REAL (KIND=8), ALLOCATABLE :: k_pw_dnfH(:,:,:), fattori_pw_dnfH(:,:), autoenergie_dnfH(:)
	LOGICAL, SAVE :: flag_simm_lda
	INTEGER, ALLOCATABLE, SAVE :: num_fpw_lda(:)
	REAL (KIND=8), ALLOCATABLE :: k_fpw_lda(:,:,:)
	COMPLEX (KIND=8), ALLOCATABLE :: fattori_fpw_lda(:,:)
	
	CONTAINS
	
	SUBROUTINE inizializza_dati_funzione_onda()
		USE dati_mc
		USE dati_fisici
		USE generic_tools
		USE momenta
		USE dnfH
		USE walkers
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
		LOGICAL :: flag_file
		CHARACTER(LEN=4) :: istring
		INTEGER :: M, i, j, max_num_pw, i1, i_twist, ios
		REAL (KIND=8) :: dummy1, SL
		REAL (KIND=8), ALLOCATABLE :: dummy(:,:,:)
		NAMELIST /dati_funzione_onda/ SDe_kind, Jee_kind, Jep_kind, Jpp_kind, SDse_kind, Jse_kind, &
		  Kse_kind, Jsesp_kind, split_Aee, split_Aep, split_Asese, split_Asesp, split_Fee, split_Fep, split_Fsese, &
        split_Fsesp, m_Bsplep, nknots_Bsplep, cutoff_Bsplep, &
        m_Jsplee, nknots_Jsplee, cutoff_Jsplee, m_Jsplep, nknots_Jsplep, cutoff_Jsplep, &
		  Aee_yuk, Aee_ud_yuk, Fee_yuk, Fee_ud_yuk, Aep_yuk, Aep_ud_yuk, Fep_yuk, Fep_ud_yuk, &
		  C_kern_e, Asese_yuk, Asese_ud_yuk, Fsese_yuk, Fsese_ud_yuk, Asesp_yuk, Asesp_ud_yuk, Fsesp_yuk, &
		  Fsesp_ud_yuk, Gswf, C_atm, N_ritraccia_coppie, N_mc_relax_traccia_coppie, A_POT_se, D_POT_se, Gsesp, c_se, &
		  B_se, D_se, c_sesp, lda_path, kf_coeff_dnfH, flag_usa_coeff_dnfH
		
		IF (iniz_funzione_onda) STOP 'funzione_onda é giá inizializzato &
		  [ module_funzione_onda.f90 > inizializza_dati_funzione_onda ]'

      CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
				
		OPEN (2, FILE=path_dati_funzione_onda,STATUS='OLD')
		READ (2,NML=dati_funzione_onda)
		CLOSE (2)

		IF (.NOT. iniz_walkers) STOP 'Non puoi inizializzare la funzione d onda prima di aver inizializzato i walkers &
		  [ module_funzione_onda.f90 > inizializza_dati_funzione_onda ]'
		IF ( Jse_kind=='pot' ) THEN
			flag_traccia_coppie=.TRUE.
			CALL attiva_traccia_coppie_mol_ss()
		ELSE
			flag_traccia_coppie=.FALSE.
		END IF
												
		IF (.NOT. iniz_dati_fisici) STOP 'Non puoi inizializzare la funzione d onda prima di aver letto i dati fisici &
		  [ module_funzione_onda.f90 > inizializza_dati_funzione_onda ]'
				
		IF ((SDe_kind=='lda').OR.(SDse_kind=='lda')) THEN
			
         !!DEPRECATED, it does not work to call Quantum Espresso internally
			!IF ( TRIM(lda_path)=="genera_on_the_fly" ) THEN
			!	CALL genera_orbitali_lda()
			!END IF
			
			IF (flag_TABC) THEN
            lda_path=TRIM(lda_path)//"/OUT.save"
				flag_simm_lda=.FALSE.
				IF ( mpi_myrank==0 ) THEN
					INQUIRE(FILE=TRIM(lda_path)//'/.numero_cartelle_K',EXIST=flag_file)
					IF ( flag_file ) CALL SYSTEM ('rm -f '//TRIM(lda_path)//'/.numero_cartelle_K')
					CALL SYSTEM ('ls -d  '//TRIM(lda_path)//'/K* | wc -w > '//TRIM(lda_path)//'/.numero_cartelle_K')
					OPEN(UNIT=66, FILE=TRIM(lda_path)//'/.numero_cartelle_K', STATUS='OLD', IOSTAT=ios)
					IF ( ios /= 0 ) STOP "Errore ad aprire il file per determinare il numero di cartelle K &
					  [inizializza_dati_funzione_onda.f90]"
					READ (UNIT=66, FMT=*, IOSTAT=ios) num_K_points
					CLOSE (66)
				END IF
				!num_K_points=6
				!PRINT * , 'iniz_wf: ', mpi_myrank, num_K_points
				CALL MPI_BCAST(num_K_points,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)
				!PRINT * , 'iniz_wf fine. ', mpi_myrank, num_K_points
				!CALL determina_numero_cartelle_K(TRIM(lda_path),mpi_myrank,num_K_points)
				ALLOCATE(pesi_K_points(1:num_K_points))
				CALL leggi_pesi_K(TRIM(lda_path),num_K_points,pesi_K_points(1:num_K_points))
				dummy1=SUM(pesi_K_points)
				pesi_K_points=pesi_K_points/dummy1
				DO i1 = 2, num_K_points, 1
					pesi_K_points(i1)=pesi_K_points(i1)+pesi_K_points(i1-1)
				END DO
				CALL RANDOM_NUMBER(dummy1)
				i_twist=1
				DO WHILE ( dummy1>pesi_K_points(i_twist+1) )
					i_twist=i_twist+1
				END DO
				IF (i_twist>num_K_points) STOP 'Errore: i_twist é maggiore di num_K_points &
				  [ module_funzione_onda.f90 > inizializza_dati_funzione_onda ]'
				WRITE (istring, '(I4.4)'), i_twist
				CALL leggi_N_pw(TRIM(lda_path)//'/K0'//istring//'/gkvectors.xml',N_pw_lda)
				ALLOCATE(k_pw_lda(0:3,1:N_pw_lda),fattori_orb_lda(1:H_N_part),fattori_pw_lda(1:N_pw_lda,1:H_N_part))
				ALLOCATE(k_pw_int_lda(1:3,1:N_pw_lda),twist_lda(1:3))
				CALL leggi_evc_xml(H_N_part,N_pw_lda,TRIM(lda_path)//'/K0'//istring//'/evc.xml',fattori_pw_lda)       !o!
				CALL leggi_gkvectors_xml(N_pw_lda,TRIM(lda_path)//'/K0'//istring//'/gkvectors.xml',k_pw_int_lda(1:3,1:N_pw_lda),twist_lda(1:3))
				DO i = 1, 3, 1
					twist_lda(i)=twist_lda(i)/r_s
				END DO
				DO j = 1, N_pw_lda, 1
					DO i = 1, 3, 1
						k_pw_lda(i,j)=(2.d0*PI/L(i))*REAL(k_pw_int_lda(i,j),8)+twist_lda(i)
					END DO
				END DO
			ELSE
            lda_path=TRIM(lda_path)//"/OUT.save/K00001"
				CALL leggi_N_pw(TRIM(lda_path)//'/gkvectors.xml',N_pw_lda)
				ALLOCATE(k_pw_lda(0:3,1:N_pw_lda),fattori_orb_lda(1:H_N_part),fattori_pw_lda(1:N_pw_lda,1:H_N_part))
				ALLOCATE(k_pw_int_lda(1:3,1:N_pw_lda),twist_lda(1:3))
				CALL leggi_evc_xml(H_N_part,N_pw_lda,TRIM(lda_path)//'/evc.xml',fattori_pw_lda)
				CALL leggi_gkvectors_xml(N_pw_lda,TRIM(lda_path)//'/gkvectors.xml',k_pw_int_lda(1:3,1:N_pw_lda),twist_lda(1:3))
				DO i = 1, 3, 1
					twist_lda(i)=twist_lda(i)/r_s
				END DO
				DO j = 1, N_pw_lda, 1
					DO i = 1, 3, 1
						k_pw_lda(i,j)=(2.d0*PI/L(i))*REAL(k_pw_int_lda(i,j),8)+twist_lda(i)
					END DO
				END DO
				IF (k_pw_int_lda(1,2)==-1) THEN
					flag_simm_lda=.FALSE.
				ELSE
					flag_simm_lda=.TRUE.
				END IF
			END IF
			
			DEALLOCATE(k_pw_int_lda)
			DO i = 1, N_pw_lda, 1
				k_pw_lda(0,i)=DOT_PRODUCT(k_pw_lda(1:3,i),k_pw_lda(1:3,i))
			END DO
			
			!trovo i coefficienti per velocizzare il metodo
			max_num_pw=0
			DO j = 1, H_N_part, 1
				i1=0
				DO i = 1, N_pw_lda, 1
					IF (REAL(fattori_pw_lda(i,j)*DCONJG(fattori_pw_lda(i,j)),4)<CUT_LDA) THEN
						
					ELSE
						i1=i1+1
					END IF
				END DO
				max_num_pw=MAX(max_num_pw,i1)
			END DO
			ALLOCATE(num_fpw_lda(1:H_N_part),k_fpw_lda(0:3,1:max_num_pw,1:H_N_part),fattori_fpw_lda(1:max_num_pw,1:H_N_part))
			num_fpw_lda=0
			DO j = 1, H_N_part, 1
				i1=0
				DO i = 1, N_pw_lda, 1
					IF (REAL(fattori_pw_lda(i,j)*DCONJG(fattori_pw_lda(i,j)),4)<CUT_LDA) THEN
						
					ELSE
						i1=i1+1
						num_fpw_lda(j)=num_fpw_lda(j)+1
						k_fpw_lda(0:3,num_fpw_lda(j),j)=k_pw_lda(0:3,i)
						fattori_fpw_lda(num_fpw_lda(j),j)=fattori_pw_lda(i,j)
					END IF
				END DO
			END DO
			
			!CALL normalizza_coefficienti_lda()   !se attivato, attivarlo anche in applica_twist_lda
			
		ELSE IF ((SDe_kind=='prf').OR.(SDe_kind=='fre')) THEN
			CALL inizializza_dnfH(kf_coeff_dnfH)     !inizializzo dnfH, trovano N_pw_lda
			IF (flag_TABC) CALL twist_k_dnfH(kf_coeff_dnfH)
			IF (flag_usa_coeff_dnfH) THEN                        !se voglio usare i coefficienti nella matrice dnfH, li carico
				!ALLOCATE(c_eff_dnfH(1:N_pw_lda,1:N_pw_lda))
				!INQUIRE(FILE=lda_path,EXIST=flag_file)
				!IF (flag_file) THEN
				!	OPEN (UNIT=18, FILE=lda_path, STATUS='OLD')
				!	READ (18,*), c_eff_dnfH
				!	CLOSE(18)
				!ELSE
				!	!IF (mpi_myrank==0) PRINT * , 'funzione_onda: non essendo presente un file, setto c_eff random.'
				!	!IF (flag_output) WRITE (7, *), 'funzione_onda: non essendo presente un file, setto c_eff random.'
				!	ALLOCATE(dummy(1:2,1:N_pw_lda,1:N_pw_lda))
				!	CALL RANDOM_NUMBER(dummy)
				!	dummy=dummy-0.5d0
				!	DO j = 1, N_pw_lda, 1
				!		DO i = 1, N_pw_lda, 1
				!			c_eff_dnfH(i,j)=dummy(1,i,j)+(0.d0,1.d0)*dummy(2,i,j)
				!		END DO
				!	END DO
				!	DEALLOCATE(dummy)
				!	c_eff_dnfH=1.d0*c_eff_dnfH*(hbar*hbar*k_Hartree(0,H_N_part)/(2.d0*MASS_e))
				!END IF
				!CALL costruisci_matrice_Hartree(SDe_kind,c_eff_dnfH)
				STOP 'coefficienti dnfH non ancora implementati'
			ELSE
				CALL costruisci_matrice_Hartree(SDe_kind)
			END IF
			CALL trova_soluzioni_dnfH()
			N_pw_lda=N_M_Hartree
			ALLOCATE(k_pw_lda(0:3,1:N_pw_lda),num_pw_orbit(1:N_pw_lda),fattori_pw_lda(1:N_pw_lda,1:N_pw_lda))
			ALLOCATE(autoenergie_dnfH(1:N_pw_lda),num_pw_dnfH(1:H_N_part))
			k_pw_lda(0:3,1:N_pw_lda)=k_Hartree
			autoenergie_dnfH=autovalori_Hartree
			num_pw_orbit=n_fpw_Hartee
			num_pw_dnfH(1:H_N_part)=num_pw_orbit(1:H_N_part)
			max_num_pw=0
			DO j = 1, N_pw_lda, 1
				IF (num_pw_orbit(j)>max_num_pw) max_num_pw=num_pw_orbit(j)
			END DO
			ALLOCATE(indice_pw_dnfH(1:max_num_pw,1:N_pw_lda))
			indice_pw_dnfH(1:max_num_pw,1:N_pw_lda)=i_fpw_Hartree(1:max_num_pw,1:N_pw_lda)
			fattori_pw_lda=M_Hartree
			ALLOCATE(k_pw_dnfH(0:3,1:max_num_pw,1:N_pw_lda))
			DO j = 1, N_pw_lda, 1
				DO i = 1, num_pw_orbit(j), 1
					k_pw_dnfH(0:3,i,j)=k_pw_lda(0:3,indice_pw_dnfH(i,j))
				END DO
			END DO
			ALLOCATE(fattori_pw_dnfH(1:max_num_pw,1:N_pw_lda))
			DO j = 1, N_pw_lda, 1
				DO i = 1, num_pw_orbit(j), 1
					fattori_pw_dnfH(i,j)=fattori_pw_lda(indice_pw_dnfH(i,j),j)
				END DO
			END DO
			CALL chiudi_dnfH()

      ELSE IF (SDe_kind=='spb') THEN
         !costruisce spline per backflow
         SL=DSQRT(DOT_PRODUCT(H_L,H_L))
         CALL MSPL_new(M=m_Bsplep , NKNOTS=nknots_Bsplep , LA=0.d0 , LB=SL , SPL=Bsplep, &
            CUTOFF=cutoff_Bsplep )
         INQUIRE(FILE=TRIM(path_dati_funzione_onda)//"-Bsplep",EXIST=flag_file)
         IF (flag_file) THEN
            CALL MSPL_load(SPL=foo_spl,FILENAME=TRIM(path_dati_funzione_onda)//"-Bsplep")
            CALL MSPL_carbon_copy(ORIGINAL_SPL=foo_spl,CC_SPL=Bsplep)
            CALL MSPL_deallocate(SPL=foo_spl)
         ELSE
            CALL MSPL_fit_function(SPL=Bsplep,F=Bep1s)
         END IF
         !CALL MSPL_print_on_file(SPL=Bsplep,DERIV=0,FILENAME="prova.d",NPOINTS=1000)
         !STOP
		END IF

      !Costruisce spline per il Jee
      IF ((Jee_kind=='spl').OR.(Jee_kind=='spp')) THEN
         SELECT CASE(Jee_kind)
         CASE('spl')
            SL=DSQRT(DOT_PRODUCT(H_L,H_L))
            !IF (crystal_cell=='mol__') SL=MIN(SL,DSQRT(3.d0*20.d0))
         CASE('spp')
            SL=DSQRT(DOT_PRODUCT(L,L))/PI
         END SELECT

         CALL MSPL_new(M=m_Jsplee , NKNOTS=nknots_Jsplee , LA=0.d0 , LB=SL , SPL=Jsplee, &
            CUTOFF=cutoff_Jsplee )
         INQUIRE(FILE=TRIM(path_dati_funzione_onda)//"-Jsplee",EXIST=flag_file)
         IF (flag_file) THEN
            CALL MSPL_load(SPL=foo_spl,FILENAME=TRIM(path_dati_funzione_onda)//"-Jsplee")
            CALL MSPL_carbon_copy(ORIGINAL_SPL=foo_spl,CC_SPL=Jsplee)
            CALL MSPL_deallocate(SPL=foo_spl)
         ELSE
            !IF (mpi_myrank==0) PRINT *, "Dati spline Jee non presenti. Fitto il Jastrow Yukawa ",&
            !   TRIM(TRIM(path_dati_funzione_onda))//"-Jsplee"
            CALL MSPL_fit_function(SPL=Jsplee,F=Jeeyuk)
         END IF

         IF (split_Aee.OR.split_Fee) THEN
            CALL MSPL_new(M=m_Jsplee , NKNOTS=nknots_Jsplee , LA=0.d0 , LB=SL , SPL=Jsplee_ud, &
               CUTOFF=cutoff_Jsplee ) 
            INQUIRE(FILE=TRIM(path_dati_funzione_onda)//"-Jsplee_ud",EXIST=flag_file)
            IF (flag_file) THEN
               CALL MSPL_load(SPL=foo_spl,FILENAME=TRIM(path_dati_funzione_onda)//"-Jsplee_ud")
               CALL MSPL_carbon_copy(ORIGINAL_SPL=foo_spl,CC_SPL=Jsplee_ud)
               CALL MSPL_deallocate(SPL=foo_spl)
            ELSE
               !IF (mpi_myrank==0) PRINT *, "Dati spline Jee_ud non presenti. Fitto il Jastrow Yukawa ",&
               !   TRIM(TRIM(path_dati_funzione_onda))//"-Jsplee_ud"
               CALL MSPL_fit_function(SPL=Jsplee_ud,F=Jeeyuk_ud)
            END IF
         END IF
      END IF

      !Costruisce spline per il Jep
      IF ((Jep_kind=='spl').OR.(Jep_kind=='spp')) THEN
         SELECT CASE(Jep_kind)
         CASE('spl')
            SL=DSQRT(DOT_PRODUCT(H_L,H_L))
            !IF (crystal_cell=='mol__') SL=MIN(SL,DSQRT(3.d0*20.d0))
         CASE('spp')
            SL=DSQRT(DOT_PRODUCT(L,L))/PI
         END SELECT

         CALL MSPL_new(M=m_Jsplep , NKNOTS=nknots_Jsplep , LA=0.d0 , LB=SL , SPL=Jsplep, &
            CUTOFF=cutoff_Jsplep )
         INQUIRE(FILE=TRIM(path_dati_funzione_onda)//"-Jsplep",EXIST=flag_file)
         IF (flag_file) THEN
            CALL MSPL_load(SPL=foo_spl,FILENAME=TRIM(path_dati_funzione_onda)//"-Jsplep")
            CALL MSPL_carbon_copy(ORIGINAL_SPL=foo_spl,CC_SPL=Jsplep)
            CALL MSPL_deallocate(SPL=foo_spl)
         ELSE
            !IF (mpi_myrank==0) PRINT *, "Dati spline Jep non presenti. Fitto il Jastrow Yukawa ",&
            !   TRIM(TRIM(path_dati_funzione_onda))//"-Jsplep"
            CALL MSPL_fit_function(SPL=Jsplep,F=Jepyuk)
         END IF

         IF (split_Aee.OR.split_Fee) THEN
            CALL MSPL_new(M=m_Jsplep , NKNOTS=nknots_Jsplep , LA=0.d0 , LB=SL , SPL=Jsplep_ud, &
               CUTOFF=cutoff_Jsplep ) 
            INQUIRE(FILE=TRIM(path_dati_funzione_onda)//"-Jsplep_ud",EXIST=flag_file)
            IF (flag_file) THEN
               CALL MSPL_load(SPL=foo_spl,FILENAME=TRIM(path_dati_funzione_onda)//"-Jsplep_ud")
               CALL MSPL_carbon_copy(ORIGINAL_SPL=foo_spl,CC_SPL=Jsplep_ud)
               CALL MSPL_deallocate(SPL=foo_spl)
            ELSE
               !IF (mpi_myrank==0) PRINT *, "Dati spline Jep_ud non presenti. Fitto il Jastrow Yukawa ",&
               !   TRIM(TRIM(path_dati_funzione_onda))//"-Jsplep_ud"
               CALL MSPL_fit_function(SPL=Jsplep_ud,F=Jepyuk_ud)
            END IF
         END IF
      END IF

		IF ((Jsesp_kind=='bou') .OR. (Jsesp_kind=='ppb')) THEN   
			M=fattoriale_doppio(N_part-1)
			ALLOCATE(coppie(1:2,1:N_part/2,1:M))
			CALL crea_tutte_le_coppie(N_part,M,coppie)
		END IF
		
		IF (( Kse_kind=='gsc' ).OR.( Kse_kind=='gdc' )) THEN
			l1_kern_e=MINVAL(L)/6.d0
			alpha1_kern_e=C_kern_e*((MINVAL(L)/6.d0-MINVAL(L)/2.d0)**3)
			alpha0_kern_e=C_kern_e*( MINVAL(L)*MINVAL(L)/36.d0 - (MINVAL(L)/6.d0-MINVAL(L)/2.d0)**2 )
			n_kern_e=12     !minimo 2
			eps_kern_e=MINVAL(L)/REAL(n_kern_e+1,8)
			l2_kern_e=MINVAL(L)*0.5d0-eps_kern_e
			IF (l2_kern_e<l1_kern_e) STOP 'l2_kern_e minore di l1_kern_e &
			  [ module_funzione_onda.f90 > inizializza_dati_funzione_onda ]'
			beta1_kern_e=-2.d0*alpha1_kern_e/((eps_kern_e**3)*n_kern_e*(n_kern_e-1)*(l2_kern_e**(n_kern_e-2)))
			beta0_kern_e=alpha0_kern_e-alpha1_kern_e/eps_kern_e-beta1_kern_e*(l2_kern_e**n_kern_e)
			!IF ( mpi_myrank==0 ) THEN
			!	PRINT * , "alpha1_kern_e=", alpha1_kern_e
			!	PRINT * , "alpha0_kern_e=", alpha0_kern_e
			!	PRINT * , "beta1_kern_e=", beta1_kern_e
			!	PRINT * , "beta0_kern_e=", beta0_kern_e
			!END IF
		END IF
		IF ( Kse_kind=='atc' ) THEN
			alpha1_kern_e=C_kern_e*(1000.d0*MINVAL(L)*0.5d0-3.d0*MINVAL(L)/8.d0)/((MINVAL(L)*0.5d0-3.d0*MINVAL(L)/8.d0)**10)
			alpha0_kern_e=C_kern_e*3.d0*MINVAL(L)/8.d0
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: Ho letto le informazioni riguardo la funzione d onda da usare da &
		  file e ho inizializzato la funzione d onda'
		
		iniz_funzione_onda=.TRUE.
				
	END SUBROUTINE inizializza_dati_funzione_onda

!-----------------------------------------------------------------------
   FUNCTION Jeeyuk(i,x)
      USE dati_fisici
      IMPLICIT NONE
      REAL(KIND=8) :: Jeeyuk
      INTEGER, INTENT(IN) :: i
      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8) :: y
      IF (x<1.d-5) THEN
         y=1.d-5
      ELSE 
         y=x
      END IF
      SELECT CASE(i)
      CASE(0)
         Jeeyuk = Aee_yuk*(1.d0-DEXP(-Fee_yuk*y))/y
      CASE(1)
         Jeeyuk = Aee_yuk*(1.d0 - DEXP(Fee_yuk*y) + Fee_yuk*y)/(DEXP(Fee_yuk*y)*(y*y)) 
      CASE(2)
         Jeeyuk = (Aee_yuk*(-2.d0 + 2.d0* DEXP(Fee_yuk* y) - 2.d0* Fee_yuk* y - (Fee_yuk*y)**2))/(DEXP(Fee_yuk*y)*(y**3))
      CASE DEFAULT
		   STOP 'La spline richiede una derivata di ordine superiore al 2, che non e stata implementata &
		      [ module_funzione_onda.f90 > Jeeyuk ]'
      END SELECT
   END FUNCTION Jeeyuk
!-----------------------------------------------------------------------
   FUNCTION Jeeyuk_ud(i,x)
      USE dati_fisici
      IMPLICIT NONE
      REAL(KIND=8) :: Jeeyuk_ud
      INTEGER, INTENT(IN) :: i
      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8) :: y
      IF (x<1.d-5) THEN
         y=1.d-5
      ELSE 
         y=x
      END IF
      SELECT CASE(i)
      CASE(0)
         Jeeyuk_ud = Aee_yuk*(1.d0-DEXP(-Fee_yuk*y))/y
      CASE(1)
         Jeeyuk_ud = Aee_yuk*(1.d0 - DEXP(Fee_yuk*y) + Fee_yuk*y)/(DEXP(Fee_yuk*y)*(y*y)) 
      CASE(2)
         Jeeyuk_ud = (Aee_yuk*(-2.d0 + 2.d0* DEXP(Fee_yuk* y) - 2.d0* Fee_yuk* y - (Fee_yuk*y)**2))/(DEXP(Fee_yuk*y)*(y**3))
      CASE DEFAULT
		   STOP 'La spline richiede una derivata di ordine superiore al 2, che non e stata implementata &
		      [ module_funzione_onda.f90 > Jeeyuk ]'
      END SELECT
   END FUNCTION Jeeyuk_ud
!-----------------------------------------------------------------------
   FUNCTION Jepyuk(i,x)
      USE dati_fisici
      IMPLICIT NONE
      REAL(KIND=8) :: Jepyuk
      INTEGER, INTENT(IN) :: i
      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8) :: y
      IF (x<1.d-5) THEN
         y=1.d-5
      ELSE 
         y=x
      END IF
      SELECT CASE(i)
      CASE(0)
         Jepyuk = Aep_yuk*(1.d0-DEXP(-Fep_yuk*y))/y
      CASE(1)
         Jepyuk = Aep_yuk*(1.d0 - DEXP(Fep_yuk*y) + Fep_yuk*y)/(DEXP(Fep_yuk*y)*(y*y)) 
      CASE(2)
         Jepyuk = (Aep_yuk*(-2.d0 + 2.d0* DEXP(Fep_yuk* y) - 2.d0* Fep_yuk* y - (Fep_yuk*y)**2))/(DEXP(Fep_yuk*y)*(y**3))
      CASE DEFAULT
		   STOP 'La spline richiede una derivata di ordine superiore al 2, che non e stata implementata &
		      [ module_funzione_onda.f90 > Jeeyuk ]'
      END SELECT
   END FUNCTION Jepyuk
!-----------------------------------------------------------------------
   FUNCTION Jepyuk_ud(i,x)
      USE dati_fisici
      IMPLICIT NONE
      REAL(KIND=8) :: Jepyuk_ud
      INTEGER, INTENT(IN) :: i
      REAL(KIND=8), INTENT(IN) :: x
      REAL(KIND=8) :: y
      IF (x<1.d-5) THEN
         y=1.d-5
      ELSE 
         y=x
      END IF
      SELECT CASE(i)
      CASE(0)
         Jepyuk_ud = Aep_yuk*(1.d0-DEXP(-Fep_yuk*y))/y
      CASE(1)
         Jepyuk_ud = Aep_yuk*(1.d0 - DEXP(Fep_yuk*y) + Fep_yuk*y)/(DEXP(Fep_yuk*y)*(y*y)) 
      CASE(2)
         Jepyuk_ud = (Aep_yuk*(-2.d0 + 2.d0* DEXP(Fep_yuk* y) - 2.d0* Fep_yuk* y - (Fep_yuk*y)**2))/(DEXP(Fep_yuk*y)*(y**3))
      CASE DEFAULT
		   STOP 'La spline richiede una derivata di ordine superiore al 2, che non e stata implementata &
		      [ module_funzione_onda.f90 > Jeeyuk ]'
      END SELECT
   END FUNCTION Jepyuk_ud
!-----------------------------------------------------------------------
   FUNCTION Bep1s(i,x)
      IMPLICIT NONE
      REAL(KIND=8) :: Bep1s
      INTEGER, INTENT(IN) :: i
      REAL(KIND=8), INTENT(IN) :: x

      SELECT CASE(i)
      CASE(0)
         Bep1s=C_atm*x
      CASE(1)
         Bep1s=C_atm
      CASE DEFAULT
         Bep1s=0.d0
      END SELECT
   
   END FUNCTION Bep1s
!-----------------------------------------------------------------------

	SUBROUTINE stampa_file_dati_funzione_onda(nome_file)
		IMPLICIT NONE
		CHARACTER(LEN=*) :: nome_file
		NAMELIST /dati_funzione_onda/ SDe_kind, Jee_kind, Jep_kind, Jpp_kind, SDse_kind, Jse_kind, &
		  Kse_kind, Jsesp_kind, split_Aee, split_Aep, split_Asese, split_Asesp, split_Fee, split_Fep, split_Fsese, &
		  split_Fsesp, m_Bsplep, nknots_Bsplep, cutoff_Bsplep, m_Jsplee, nknots_Jsplee, &
        cutoff_Jsplee, m_Jsplep, nknots_Jsplep, cutoff_Jsplep, &
        Aee_yuk, Aee_ud_yuk, Fee_yuk, Fee_ud_yuk, Aep_yuk, Aep_ud_yuk, Fep_yuk, Fep_ud_yuk, &
		  C_kern_e, Asese_yuk, Asese_ud_yuk, Fsese_yuk, Fsese_ud_yuk, Asesp_yuk, Asesp_ud_yuk, Fsesp_yuk, &
		  Fsesp_ud_yuk, Gswf, C_atm, N_ritraccia_coppie, N_mc_relax_traccia_coppie, A_POT_se, D_POT_se, Gsesp, c_se, &
		  B_se, D_se, c_sesp, lda_path, kf_coeff_dnfH, flag_usa_coeff_dnfH
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > stampa_file_dati_funzione_onda ]'
		
		OPEN (2, FILE=nome_file,STATUS='UNKNOWN')
		IF (.NOT. split_Aee) Aee_ud_yuk=Aee_yuk
		IF (.NOT. split_Aep) Aep_ud_yuk=Aep_yuk
		IF (.NOT. split_Asese) Asese_ud_yuk=Asese_yuk
		IF (.NOT. split_Asesp) Asesp_ud_yuk=Asesp_yuk
		IF (.NOT. split_Fee) Fee_ud_yuk=Fee_yuk
		IF (.NOT. split_Fep) Fep_ud_yuk=Fep_yuk
		IF (.NOT. split_Fsese) Fsese_ud_yuk=Fsese_yuk
		IF (.NOT. split_Fsesp) Fsesp_ud_yuk=Fsesp_yuk
		kf_coeff_dnfH=DABS(kf_coeff_dnfH)
		WRITE (2,NML=dati_funzione_onda)
		CLOSE (2)
		
		IF (((SDe_kind=='prf').OR.(SDe_kind=='fre')).AND.(flag_usa_coeff_dnfH)) THEN
			OPEN (2, FILE=nome_file//'.c_eff_dnfH',STATUS='UNKNOWN')
			WRITE (2, *), c_eff_dnfH
			CLOSE(2)
		END IF

      IF (SDe_kind=='spb') THEN
         CALL MSPL_store(SPL=Bsplep,FILENAME=nome_file//"-Bsplep")
      END IF

      IF ((Jee_kind=='spl').OR.(Jee_kind=='spp')) THEN
         CALL MSPL_store(SPL=Jsplee,FILENAME=nome_file//"-Jsplee")
         IF (split_Aee.OR.split_Fee) CALL MSPL_store(SPL=Jsplee_ud,FILENAME=nome_file//"-Jsplee_ud")
      END IF

      IF ((Jep_kind=='spl').OR.(Jep_kind=='spp')) THEN
         CALL MSPL_store(SPL=Jsplep,FILENAME=nome_file//"-Jsplep")
         IF (split_Aep.OR.split_Fep) CALL MSPL_store(SPL=Jsplep_ud,FILENAME=nome_file//"-Jsplep_ud")
      END IF
		
	END SUBROUTINE stampa_file_dati_funzione_onda
!-----------------------------------------------------------------------

	SUBROUTINE setta_parametri(nuovi_parametri,numero_parametri)
		USE dati_mc
		USE walkers
		USE dnfH
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: numero_parametri
		REAL (KIND=8), INTENT(INOUT) :: nuovi_parametri(1:numero_parametri)
		INTEGER :: cont, i, j, max_num_pw
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > setta_parametri ]'
		
		cont=1
		
		IF (Jee_kind/='no_') THEN
			SELECT CASE (Jee_kind)
			CASE ('yuk')
				IF (opt_A_Jee) THEN
					IF (opt_F_Jee) THEN
						Aee_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Aee) THEN
							Aee_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
						Fee_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Fee) THEN
							Fee_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					ELSE
						Aee_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Aee) THEN
							Aee_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					END IF
				ELSE
					IF (opt_F_Jee) THEN
						Fee_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Fee) THEN
							Fee_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					END IF
				END IF
			CASE ('yup')
				IF (opt_A_Jee) THEN
					IF (opt_F_Jee) THEN
						Aee_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Aee) THEN
							Aee_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
						Fee_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Fee) THEN
							Fee_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					ELSE
						Aee_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Aee) THEN
							Aee_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					END IF
				ELSE
					IF (opt_F_Jee) THEN
						Fee_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Fee) THEN
							Fee_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					END IF
				END IF
         CASE ('spl','spp')
            IF (opt_A_Jee.OR.opt_F_Jee) THEN
               IF (split_Aee.OR.split_Fee) THEN
                  Jsplee%t(0:Jsplee%m,0:Jsplee%Nknots)=&
                     RESHAPE(nuovi_parametri(cont:cont+(Jsplee%m+1)*(Jsplee%Nknots+1)-1),&
                     (/Jsplee%m+1,Jsplee%Nknots+1/))
                  cont=cont+(Jsplee%m+1)*(Jsplee%Nknots+1)
                  Jsplee_ud%t(0:Jsplee_ud%m,0:Jsplee_ud%Nknots)=&
                     RESHAPE(nuovi_parametri(cont:cont+(Jsplee_ud%m+1)*(Jsplee_ud%Nknots+1)-1),&
                     (/Jsplee_ud%m+1,Jsplee_ud%Nknots+1/))
                  cont=cont+(Jsplee_ud%m+1)*(Jsplee_ud%Nknots+1)
               ELSE
                  Jsplee%t(0:Jsplee%m,0:Jsplee%Nknots)=&
                     RESHAPE(nuovi_parametri(cont:cont+(Jsplee%m+1)*(Jsplee%Nknots+1)-1),&
                     (/Jsplee%m+1,Jsplee%Nknots+1/))
                  cont=cont+(Jsplee%m+1)*(Jsplee%Nknots+1)
               END IF
            END IF
			END SELECT
		END IF
		IF (Jep_kind/='no_') THEN
			SELECT CASE (Jep_kind)
			CASE ('yuk')
				IF (opt_A_Jep) THEN
					IF (opt_F_Jep) THEN
						Aep_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Aep) THEN
							Aep_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
						Fep_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Fep) THEN
							Fep_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					ELSE
						Aep_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Aep) THEN
							Aep_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					END IF
				ELSE
					IF (opt_F_Jep) THEN
						Fep_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Fep) THEN
							Fep_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					END IF
				END IF
			CASE ('yup')
				IF (opt_A_Jep) THEN
					IF (opt_F_Jep) THEN
						Aep_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Aep) THEN
							Aep_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
						Fep_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Fep) THEN
							Fep_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					ELSE
						Aep_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Aep) THEN
							Aep_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					END IF
				ELSE
					IF (opt_F_Jep) THEN
						Fep_yuk=nuovi_parametri(cont)
						cont=cont+1
						IF (split_Fep) THEN
							Fep_ud_yuk=nuovi_parametri(cont)
							cont=cont+1
						END IF
					END IF
				END IF
         CASE ('spl','spp')
            IF (opt_A_Jep.OR.opt_F_Jep) THEN
               IF (split_Aep.OR.split_Fep) THEN
                  Jsplep%t(0:Jsplep%m,0:Jsplep%Nknots)=&
                     RESHAPE(nuovi_parametri(cont:cont+(Jsplep%m+1)*(Jsplep%Nknots+1)-1),&
                     (/Jsplep%m+1,Jsplep%Nknots+1/))
                  cont=cont+(Jsplep%m+1)*(Jsplep%Nknots+1)
                  Jsplep_ud%t(0:Jsplep_ud%m,0:Jsplep_ud%Nknots)=&
                     RESHAPE(nuovi_parametri(cont:cont+(Jsplep_ud%m+1)*(Jsplep_ud%Nknots+1)-1),&
                     (/Jsplep_ud%m+1,Jsplep_ud%Nknots+1/))
                  cont=cont+(Jsplep_ud%m+1)*(Jsplep_ud%Nknots+1)
               ELSE
                  Jsplep%t(0:Jsplep%m,0:Jsplep%Nknots)=&
                     RESHAPE(nuovi_parametri(cont:cont+(Jsplep%m+1)*(Jsplep%Nknots+1)-1),&
                     (/Jsplep%m+1,Jsplep%Nknots+1/))
                  cont=cont+(Jsplep%m+1)*(Jsplep%Nknots+1)
               END IF
            END IF
			CASE ('atm')
				Fep_yuk=nuovi_parametri(cont)
				cont=cont+1
			CASE ('atp')
				Fep_yuk=nuovi_parametri(cont)
				cont=cont+1
			END SELECT
		END IF

		IF (Kse_kind/='no_') THEN
			IF (opt_Kse) THEN
				SELECT CASE (Kse_kind)
				CASE ('gss')
					C_kern_e=nuovi_parametri(cont)
					cont=cont+1
				CASE ('gsc')
					C_kern_e=nuovi_parametri(cont)
					cont=cont+1
					alpha1_kern_e=C_kern_e*((MINVAL(L)/6.d0-MINVAL(L)/2.d0)**3)
					alpha0_kern_e=C_kern_e*( MINVAL(L)*MINVAL(L)/36.d0 - (MINVAL(L)/6.d0-MINVAL(L)/2.d0)**2 )
				CASE ('gsp')
					C_kern_e=nuovi_parametri(cont)
					cont=cont+1
				CASE ('gsd')
					C_kern_e=nuovi_parametri(cont)
					cont=cont+1
				CASE ('gdc')
					C_kern_e=nuovi_parametri(cont)
					cont=cont+1
				CASE ('gdp')
					C_kern_e=nuovi_parametri(cont)
					cont=cont+1
				CASE ('atm')
					C_kern_e=nuovi_parametri(cont)
					cont=cont+1
				CASE ('atc')
					C_kern_e=nuovi_parametri(cont)
					cont=cont+1
				END SELECT
			END IF
		END IF

		IF (Jse_kind/='no_') THEN
			IF (opt_Jse) THEN
				SELECT CASE (Jse_kind)
				CASE ('pot')
					!PRINT * , 'QUA CI SONO'
					A_POT_se=nuovi_parametri(cont)
					!PRINT * , 'A_POT_se=', A_POT_se
					D_POT_se=nuovi_parametri(cont+1)
					!PRINT * , 'D_POT_se=', A_POT_se
					cont=cont+2
				CASE ('bou')
					B_se=nuovi_parametri(cont)
					D_se=nuovi_parametri(cont+1)
					cont=cont+2
				CASE ('ppb')
					c_se=nuovi_parametri(cont)
					B_se=nuovi_parametri(cont+1)
					D_se=nuovi_parametri(cont+2)
					cont=cont+3
				CASE ('yuk')
					Asese_yuk=nuovi_parametri(cont)
					cont=cont+1
					IF (split_Asese) THEN
						Asese_ud_yuk=nuovi_parametri(cont)
						cont=cont+1
					END IF
					Fsese_yuk=nuovi_parametri(cont)
					cont=cont+1
					IF (split_Fsese) THEN
						Fsese_ud_yuk=nuovi_parametri(cont)
						cont=cont+1
					END IF
				CASE ('yup')
					Asese_yuk=nuovi_parametri(cont)
					cont=cont+1
					IF (split_Asese) THEN
						Asese_ud_yuk=nuovi_parametri(cont)
						cont=cont+1
					END IF
					Fsese_yuk=nuovi_parametri(cont)
					cont=cont+1
					IF (split_Fsese) THEN
						Fsese_ud_yuk=nuovi_parametri(cont)
						cont=cont+1
					END IF
				END SELECT
			END IF
		END IF
		IF (Jsesp_kind/='no_') THEN
			IF (opt_Jsesp) THEN
				SELECT CASE (Jsesp_kind)
				CASE ('pot')
					c_sesp=nuovi_parametri(cont)
					cont=cont+1
				CASE ('yuk')
					Asesp_yuk=nuovi_parametri(cont)
					cont=cont+1
					IF (split_Asesp) THEN
						Asesp_ud_yuk=nuovi_parametri(cont)
						cont=cont+1
					END IF
					Fsesp_yuk=nuovi_parametri(cont)
					cont=cont+1
					IF (split_Fsesp) THEN
						Fsesp_ud_yuk=nuovi_parametri(cont)
						cont=cont+1
					END IF
				CASE ('yup')
					Asesp_yuk=nuovi_parametri(cont)
					cont=cont+1
					IF (split_Asesp) THEN
						Asesp_ud_yuk=nuovi_parametri(cont)
						cont=cont+1
					END IF
					Fsesp_yuk=nuovi_parametri(cont)
					cont=cont+1
					IF (split_Fsesp) THEN
						Fsesp_ud_yuk=nuovi_parametri(cont)
						cont=cont+1
					END IF
				CASE ('gss')
					Gsesp=nuovi_parametri(cont)
					cont=cont+1
				CASE ('gsd')
					Gsesp=nuovi_parametri(cont)
					cont=cont+1
				END SELECT
			END IF
		END IF

		IF (SDe_kind/='no_') THEN
			IF ((SDe_kind=='prf').OR.(SDe_kind=='fre')) THEN
				IF (opt_c_eff_dnfH) THEN
					cont=0
					DO j = 1, N_pw_lda, 1   !colonna
						DO i = 1, j-1, 1   !riga
							cont=cont+2
							c_eff_dnfH(i,j)=nuovi_parametri(cont-1)+(0.d0,1.d0)*nuovi_parametri(cont)
						END DO
					END DO
					DO j = 1, N_pw_lda, 1   !colonna
						i=j
						cont=cont+1
						c_eff_dnfH(i,j)=nuovi_parametri(cont)
					END DO
            
					CALL inizializza_dnfH(kf_coeff_dnfH)     !inizializzo dnfH, trovano N_pw_lda
					CALL costruisci_matrice_Hartree(SDe_kind,c_eff_dnfH)
					CALL trova_soluzioni_dnfH()
					N_pw_lda=N_M_Hartree
					autoenergie_dnfH=autovalori_Hartree
					num_pw_orbit=n_fpw_Hartee
					num_pw_dnfH(1:H_N_part)=num_pw_orbit(1:H_N_part)
					max_num_pw=0
					DO j = 1, N_pw_lda, 1
						IF (num_pw_orbit(j)>max_num_pw) max_num_pw=num_pw_orbit(j)
					END DO
					DEALLOCATE(indice_pw_dnfH)
					ALLOCATE(indice_pw_dnfH(1:max_num_pw,1:N_pw_lda))
					indice_pw_dnfH(1:max_num_pw,1:N_pw_lda)=i_fpw_Hartree(1:max_num_pw,1:N_pw_lda)
					fattori_pw_lda=M_Hartree
					DEALLOCATE(k_pw_dnfH)
					ALLOCATE(k_pw_dnfH(0:3,1:max_num_pw,1:N_pw_lda))
					DO j = 1, N_pw_lda, 1
						DO i = 1, num_pw_orbit(j), 1
							k_pw_dnfH(0:3,i,j)=k_pw_lda(0:3,indice_pw_dnfH(i,j))
						END DO
					END DO
					DEALLOCATE(fattori_pw_dnfH)
					ALLOCATE(fattori_pw_dnfH(1:max_num_pw,1:N_pw_lda))
					DO j = 1, N_pw_lda, 1
						DO i = 1, num_pw_orbit(j), 1
							fattori_pw_dnfH(i,j)=fattori_pw_lda(indice_pw_dnfH(i,j),j)
						END DO
					END DO
					CALL chiudi_dnfH()
				END IF
			ELSE IF ((SDe_kind=='atm').OR.(SDe_kind=='atp').OR.(SDe_kind=='bat').OR.(SDe_kind=='hl_')) THEN
				IF ( opt_SDe ) THEN
					C_atm=nuovi_parametri(cont)
					cont=cont+1
				END IF
			ELSE IF (SDe_kind=='1sb') THEN
				IF ( opt_SDe ) THEN
               IF (opt_orbital) THEN
					   C_atm=nuovi_parametri(cont)
					   cont=cont+1
               END IF
               IF (opt_dynamic_backflow) THEN
					   A_POT_se=nuovi_parametri(cont)
					   cont=cont+1
					   D_POT_se=nuovi_parametri(cont)
					   cont=cont+1
               END IF
				END IF
			ELSE IF (SDe_kind=='spb') THEN
				IF ( opt_SDe ) THEN
               IF (opt_orbital) THEN
					   Bsplep%t(0:Bsplep%m,0:Bsplep%nknots)=&
                     RESHAPE(nuovi_parametri(cont:cont+(Bsplep%m+1)*(Bsplep%nknots+1)-1),(/Bsplep%m+1,Bsplep%nknots+1/))
					   cont=cont+(Bsplep%m+1)*(Bsplep%nknots+1)
               END IF
               IF (opt_dynamic_backflow) THEN
					   A_POT_se=nuovi_parametri(cont)
					   cont=cont+1
					   D_POT_se=nuovi_parametri(cont)
					   cont=cont+1
               END IF
				END IF
			END IF
		END IF
		IF ( SDse_kind/='no_' ) THEN
			IF ( opt_SDse ) THEN
				SELECT CASE(SDse_kind)
				CASE ('atm')
					C_atm=nuovi_parametri(cont)
					cont=cont+1
				CASE ('atp')
					C_atm=nuovi_parametri(cont)
					cont=cont+1
				CASE DEFAULT
					STOP 'Non puoi minimizzare SDse se non contiene parametri variazionali &
					  [ module_funzione_onda.f90 > setta_parametri ]'
				END SELECT
			END IF
		END IF
		
		IF ( opt_Rp ) THEN
			IF (.NOT. iniz_walkers) STOP 'Prima di settare Rp devi aver inizializzato i walkers &
			  [ module_funzione_onda.f90 > setta_parametri ]'
			CALL setta_Rcrystal(nuovi_parametri(cont:cont+3*N_part-1))
		END IF
		
	END SUBROUTINE setta_parametri
!-----------------------------------------------------------------------

	SUBROUTINE stampa_parametri_variazionali()
		USE dati_mc
		USE dati_fisici
		IMPLICIT NONE
		INTEGER :: media_num_pw
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > stampa_parametri_variazionali ]'
				
		IF (SDe_kind/='no_') THEN
			SELECT CASE (SDe_kind)
			CASE ('pw_')
				PRINT * , '     SDe: ', SDe_kind
				IF (flag_output) WRITE (7, *), '     SDe: ', SDe_kind
			CASE ('lda')
				!PRINT * , '     SDe: ', SDe_kind, '          Numero onde piane:', N_pw_lda
				!IF (flag_output) WRITE (7, *), '     SDe: ', SDe_kind, '          Numero onde piane:', N_pw_lda
				PRINT * , '     SDe: ', SDe_kind, '          Numero onde piane:', REAL(SUM(num_fpw_lda))/REAL(H_N_part)
				IF (flag_output) WRITE (7, *), '     SDe: ', SDe_kind, '          Numero onde piane:', REAL(SUM(num_fpw_lda))/REAL(H_N_part)
			CASE ('prf')
				media_num_pw=SUM(num_pw_orbit(1:N_pw_lda))/N_pw_lda
				PRINT * , '     SDe: ', SDe_kind, '          Numero onde piane:', N_pw_lda, &
				  '     (', media_num_pw, ' per orbitale)      ctf=', REAL(kf_coeff_dnfH,4)
				IF (flag_output) WRITE (7, *), '     SDe: ', SDe_kind, '          Numero onde piane:', N_pw_lda, &
				  '     (', media_num_pw, ' per orbitale)      ctf=', REAL(kf_coeff_dnfH,4)
			CASE ('fre')
				media_num_pw=SUM(num_pw_orbit(1:N_pw_lda))/N_pw_lda
				PRINT * , '     SDe: ', SDe_kind, '          Numero onde piane:', N_pw_lda, &
				  '     (', media_num_pw, ' per orbitale)      ctf=', REAL(kf_coeff_dnfH,4)
				IF (flag_output) WRITE (7, *), '     SDe: ', SDe_kind, '          Numero onde piane:', N_pw_lda, &
				  '     (', media_num_pw, ' per orbitale)      ctf=', REAL(kf_coeff_dnfH,4)
			CASE ('atm')
				PRINT '(6X,A5,A3,A11,F9.3)' , 'SDe: ', SDe_kind,'  -  C_atm=', C_atm
				IF (flag_output) WRITE (7, '(6X,A5,A3,A11,F9.3)'), &
				  'SDe: ', SDe_kind,'  -  C_atm=', C_atm
			CASE ('atp')
				PRINT '(6X,A5,A3,A11,F9.3)' , 'SDe: ', SDe_kind,'  -  C_atm=', C_atm
				IF (flag_output) WRITE (7, '(6X,A5,A3,A11,F9.3)'), &
				  'SDe: ', SDe_kind,'  -  C_atm=', C_atm
	  		CASE ('bat','hl_')
	  			PRINT '(6X,A5,A3,A11,F9.3)' , 'SDe: ', SDe_kind,'  -  C_atm=', C_atm
	  			IF (flag_output) WRITE (7, '(6X,A5,A3,A11,F9.3)'), &
	  			  'SDe: ', SDe_kind,'  -  C_atm=', C_atm
	  		CASE ('1sb')
	  			PRINT '(6X,A5,A3,A11,F9.3,2(5X,A9,F9.3))' , 'SDe: ', SDe_kind,'  -  C_atm=', C_atm, &
               'A_POT_se=', A_POT_se, 'D_POT_se=', D_POT_se
	  			IF (flag_output) WRITE (7, '(6X,A5,A3,A11,F9.3,2(5X,A9,F9.3))'), &
	  			  'SDe: ', SDe_kind,'  -  C_atm=', C_atm, 'A_POT_se=', A_POT_se, 'D_POT_se=', D_POT_se
         CASE ('spb')
	  			PRINT '(6X,A5,A3,A9,I1,5X,A7,I4,5X,A7,L1,2(5X,A9,F9.3))' , 'SDe: ', SDe_kind,'   -   m=', m_Bsplep, &
               'nknots=', nknots_Bsplep, 'cutoff=', cutoff_Bsplep, 'A_POT_se=', A_POT_se, 'D_POT_se=', D_POT_se
	  			IF (flag_output) WRITE (7, '(6X,A5,A3,A7,I1,5X,A7,I4,5X,A7,L1,2(5X,A9,F9.3))'), &
	  			  'SDe: ', SDe_kind,'  -  m=', m_Bsplep, 'nknots=', nknots_Bsplep, 'cutoff=', cutoff_Bsplep,&
              'A_POT_se=', A_POT_se, 'D_POT_se=', D_POT_se
			END SELECT
		END IF

		IF (Jee_kind/='no_') THEN
			SELECT CASE (Jee_kind)
			CASE ('yuk')
				IF (split_Aee) THEN
					IF (split_Fee) THEN
						PRINT '(6X,A5,A3,A15,F9.3,5X,A11,F9.3,5X,A8,F9.3,5X,A11,F9.3)' , &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, 'Aee_ud_yuk=', Aee_ud_yuk, 'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,F9.3,5X,A11,F9.3,5X,A8,F9.3,5X,A11,F9.3)'), &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk,'Aee_ud_yuk=', Aee_ud_yuk, 'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A15,F9.3,5X,A11,F9.3,5X,A8,F9.3)' , &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, 'Aee_ud_yuk=', Aee_ud_yuk, 'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,F9.3,5X,A11,F9.3,5X,A8,F9.3)'), 'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, &
						                  'Aee_ud_yuk=', Aee_ud_yuk, 'Fee_yuk=', Fee_yuk
					END IF
				ELSE
					IF (split_Fee) THEN
						PRINT '(6X,A5,A3,A15,F9.3,5X,A8,F9.3,5X,A11,F9.3)' , &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, 'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,F9.3,5X,A8,F9.3,5X,A11,F9.3)'), 'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, &
						                  'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A15,F9.3,5X,A8,F9.3)' , &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, 'Fee_yuk=', Fee_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,F9.3,5X,A8,F9.3)'), 'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, &
						                  'Fee_yuk=', Fee_yuk
					END IF
				END IF
			CASE ('yup')
				IF (split_Aee) THEN
					IF (split_Fee) THEN
						PRINT '(6X,A5,A3,A15,F9.3,5X,A11,F9.3,5X,A8,F9.3,5X,A11,F9.3)' , &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, 'Aee_ud_yuk=', Aee_ud_yuk, 'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,F9.3,5X,A11,F9.3,5X,A8,F9.3,5X,A11,F9.3)'), &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk,'Aee_ud_yuk=', Aee_ud_yuk, 'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A15,F9.3,5X,A11,F9.3,5X,A8,F9.3)' , &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, 'Aee_ud_yuk=', Aee_ud_yuk, 'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,F9.3,5X,A11,F9.3,5X,A8,F9.3)'), 'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, &
						                  'Aee_ud_yuk=', Aee_ud_yuk, 'Fee_yuk=', Fee_yuk
					END IF
				ELSE
					IF (split_Fee) THEN
						PRINT '(6X,A5,A3,A15,F9.3,5X,A8,F9.3,5X,A11,F9.3)' , &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, 'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,F9.3,5X,A8,F9.3,5X,A11,F9.3)'), 'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, &
						                  'Fee_yuk=', Fee_yuk, 'Fee_ud_yuk=', Fee_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A15,F9.3,5X,A8,F9.3)' , &
						  'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, 'Fee_yuk=', Fee_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,F9.3,5X,A8,F9.3)'), 'Jee: ',Jee_kind,'   -   Aee_yuk=', Aee_yuk, &
						                  'Fee_yuk=', Fee_yuk
					END IF
				END IF
         CASE ('spl','spp')
            PRINT '(6X,A5,A3,A9,I1,5X,A7,I4,5X,A7,L1,5X,A11,L1)', &
               'Jee: ', Jee_kind, '   -   m=', m_Jsplee,'nknots=',nknots_Jsplee,'cutoff=',cutoff_Jsplee,'spin-split=',(split_Aee&
               .OR.split_Fee)
            IF (flag_output)  WRITE(UNIT=7, FMT='(6X,A5,A3,A9,I1,5X,A7,I4,5X,A7,L1,5X,A11,L1)'), & 
               'Jee: ', Jee_kind, '   -   m=', m_Jsplee,'nknots=',nknots_Jsplee,'cutoff=',cutoff_Jsplee,'spin-split=',(split_Aee&
               .OR.split_Fee)
			END SELECT
		END IF
		
		IF (Jep_kind/='no_') THEN
			SELECT CASE (Jep_kind)
			CASE ('yuk')
				IF (split_Aep) THEN
					IF (split_Fep) THEN
						PRINT '(6X,A5,A3,A15,1F9.3,5X,A11,1F9.3,5X,A8,1F9.3,5X,A11,1F9.3)', &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Aep_ud_yuk=', Aep_ud_yuk, 'Fep_yuk=', Fep_yuk, 'Fep_ud_yuk=', Fep_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3,5X,A11,1F9.3,5X,A8,1F9.3,5X,A11,1F9.3)'), &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Aep_ud_yuk=', Aep_ud_yuk, 'Fep_yuk=', Fep_yuk, 'Fep_ud_yuk=', Fep_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A15,1F9.3,5X,A11,1F9.3,5X,A8,1F9.3)', &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Aep_ud_yuk=', Aep_ud_yuk, 'Fep_yuk=', Fep_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3,5X,A11,1F9.3,5X,A8,1F9.3)'), &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Aep_ud_yuk=', Aep_ud_yuk, 'Fep_yuk=', Fep_yuk
					END IF
				ELSE
					IF (split_Fep) THEN
						PRINT '(6X,A5,A3,A15,1F9.3,5X,A8,1F9.3,5X,A11,1F9.3)', &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Fep_yuk=', Fep_yuk, 'Fep_ud_yuk=', Fep_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3,5X,A8,1F9.3,5X,A11,1F9.3)'), &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Fep_yuk=', Fep_yuk, 'Fep_ud_yuk=', Fep_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A15,1F9.3,5X,A8,1F9.3)', &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Fep_yuk=', Fep_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3,5X,A8,1F9.3)'), &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Fep_yuk=', Fep_yuk
					END IF
				END IF
			CASE ('yup')
				IF (split_Aep) THEN
					IF (split_Fep) THEN
						PRINT '(6X,A5,A3,A15,1F9.3,5X,A11,1F9.3,5X,A8,1F9.3,5X,A11,1F9.3)', &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Aep_ud_yuk=', Aep_ud_yuk, 'Fep_yuk=', Fep_yuk, 'Fep_ud_yuk=', Fep_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3,5X,A11,1F9.3,5X,A8,1F9.3,5X,A11,1F9.3)'), &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Aep_ud_yuk=', Aep_ud_yuk, 'Fep_yuk=', Fep_yuk, 'Fep_ud_yuk=', Fep_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A15,1F9.3,5X,A11,1F9.3,5X,A8,1F9.3)', &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Aep_ud_yuk=', Aep_ud_yuk, 'Fep_yuk=', Fep_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3,5X,A11,1F9.3,5X,A8,1F9.3)'), &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Aep_ud_yuk=', Aep_ud_yuk, 'Fep_yuk=', Fep_yuk
					END IF
				ELSE
					IF (split_Fep) THEN
						PRINT '(6X,A5,A3,A15,1F9.3,5X,A8,1F9.3,5X,A11,1F9.3)', &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Fep_yuk=', Fep_yuk, 'Fep_ud_yuk=', Fep_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3,5X,A8,1F9.3,5X,A11,1F9.3)'), &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Fep_yuk=', Fep_yuk, 'Fep_ud_yuk=', Fep_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A15,1F9.3,5X,A8,1F9.3)', &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Fep_yuk=', Fep_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3,5X,A8,1F9.3)'), &
						  'Jep: ',Jep_kind,'   -   Aep_yuk=', Aep_yuk, 'Fep_yuk=', Fep_yuk
					END IF
				END IF
         CASE ('spl','spp')
            PRINT '(6X,A5,A3,A9,I1,5X,A7,I4,5X,A7,L1,5X,A11,L1)', &
               'Jep: ', Jep_kind, '   -   m=', m_Jsplep,'nknots=',nknots_Jsplep,'cutoff=',cutoff_Jsplep,&
               'spin-split=',(split_Aep.OR.split_Fep)
            IF (flag_output)  WRITE(UNIT=7, FMT='(6X,A5,A3,A9,I1,5X,A7,I4,5X,A7,L1,5X,A11,L1)'), & 
               'Jep: ', Jep_kind, '   -   m=', m_Jsplep,'nknots=',nknots_Jsplep,'cutoff=',cutoff_Jsplep,&
               'spin-split=',(split_Aep.OR.split_Fep)
			CASE ('atm')
				PRINT '(6X,A5,A3,A15,1F9.3)', &
				  'Jep: ',Jep_kind,'   -   Fep_yuk=', Fep_yuk
				IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3)'), &
				  'Jep: ',Jep_kind,'   -   Fep_yuk=', Fep_yuk
			CASE ('atp')
				PRINT '(6X,A5,A3,A15,1F9.3)', &
				  'Jep: ',Jep_kind,'   -   Fep_yuk=', Fep_yuk
				IF (flag_output) WRITE (7, '(6X,A5,A3,A15,1F9.3)'), &
				  'Jep: ',Jep_kind,'   -   Fep_yuk=', Fep_yuk
			END SELECT
		END IF
		
		IF (Kse_kind/='no_') THEN
			SELECT CASE (Kse_kind)
			CASE ('gss')
				PRINT '(6X,A5,A3,A16,F9.3)' , &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3)'), &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
			CASE ('gsc')
				PRINT '(6X,A5,A3,A16,F9.3)' , &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3)'), &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
			CASE ('gsp')
				PRINT '(6X,A5,A3,A16,F9.3)' , &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3)'), &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
			CASE ('gsd')
				PRINT '(6X,A5,A3,A16,F9.3)' , &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3)'), &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
			CASE ('gdc')
				PRINT '(6X,A5,A3,A16,F9.3)' , &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3)'), &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
			CASE ('gdp')
				PRINT '(6X,A5,A3,A16,F9.3)' , &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3)'), &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
			CASE ('atm')
				PRINT '(6X,A5,A3,A16,F9.3)' , &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3)'), &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
			CASE ('atc')
				PRINT '(6X,A5,A3,A16,F9.3)' , &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3)'), &
				  'Kse: ',Kse_kind,'   -   C_kern_e=', C_kern_e
			END SELECT
		END IF
		
		IF (Jse_kind/='no_') THEN
			SELECT CASE (Jse_kind)
			CASE ('pot')
				PRINT '(6X,A5,A3,A16,F9.3,5X,A9,F9.3)' , &
				  'Jse: ',Jse_kind,'   -   A_POT_se=', A_POT_se, 'D_POT_se=', D_POT_se
				IF (flag_output) WRITE (7, '(6X,A5,A3,A16,F9.3,5X,A9,F9.3)'), &
				  'Jse: ',Jse_kind,'   -   A_POT_se=', A_POT_se, 'D_POT_se=', D_POT_se
			CASE ('bou')
				PRINT '(6X,A5,A3,A12,F9.3,A5,F9.3)' , &
				  'Jse: ',Jse_kind,'   -   B_se=', B_se, 'D_se=', D_se
				IF (flag_output) WRITE (7, '(6X,A5,A3,A12,F9.3,A5,F9.3)'), &
				  'Jse: ',Jse_kind,'   -   B_se=', B_se, 'D_se=', D_se
			CASE ('ppb')
				PRINT '(6X,A5,A3,A12,F9.3)' , &
				  'Jse: ',Jse_kind,'   -   c_se=', c_se
				PRINT '(6X,A5,A3,A12,F9.3,5X,A5,F9.3)' , &
				  'Jse: ',Jse_kind,'   -   B_se=', B_se, 'D_se=', D_se
				IF (flag_output) THEN
					WRITE (7, '(6X,A5,A3,A12,F9.3)'), &
					  'Jse: ',Jse_kind,'   -   c_se=', c_se
					WRITE (7, '(6X,A5,A3,A12,F9.3,5X,A5,F9.3)'), &
					  'Jse: ',Jse_kind,'   -   B_se=', B_se, 'D_se=', D_se
				END IF
			CASE ('yuk')
				IF (split_Asese) THEN
					IF (split_Fsese) THEN
						PRINT '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3,5X,A13,F9.3)' , &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=',Asese_yuk,'Asese_ud_yuk=',Asese_ud_yuk,&
						  'Fsese_yuk=',Fsese_yuk,'Fsese_ud_yuk=',Fsese_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3,5X,A13,F9.3)'), &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=',Asese_yuk,'Asese_ud_yuk=',Asese_ud_yuk,&
						  'Fsese_yuk=',Fsese_yuk,'Fsese_ud_yuk=',Fsese_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3)' , &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=',Asese_yuk,'Asese_ud_yuk=',Asese_ud_yuk,'Fsese_yuk=',Fsese_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3)'), &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=',Asese_yuk,'Asese_ud_yuk=',Asese_ud_yuk,'Fsese_yuk=',Fsese_yuk
					END IF
				ELSE
					IF (split_Fsese) THEN
						PRINT '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3)' , &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=', Asese_yuk, 'Fsese_yuk=', Fsese_yuk, 'Fsese_ud_yuk=', Fsese_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3)'), &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=', Asese_yuk, 'Fsese_yuk=', Fsese_yuk, 'Fsese_ud_yuk=', Fsese_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A17,F9.3,5X,A13,F9.3)' , &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=', Asese_yuk, 'Fsese_yuk=', Fsese_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A17,F9.3,5X,A13,F9.3)'), &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=', Asese_yuk, 'Fsese_yuk=', Fsese_yuk
					END IF
				END IF
			CASE ('yup')
				IF (split_Asese) THEN
					IF (split_Fsese) THEN
						PRINT '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3,5X,A13,F9.3)' , &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=',Asese_yuk,'Asese_ud_yuk=',Asese_ud_yuk,&
						  'Fsese_yuk=',Fsese_yuk,'Fsese_ud_yuk=',Fsese_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3,5X,A13,F9.3)'), &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=',Asese_yuk,'Asese_ud_yuk=',Asese_ud_yuk,&
						  'Fsese_yuk=',Fsese_yuk,'Fsese_ud_yuk=',Fsese_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3)' , &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=',Asese_yuk,'Asese_ud_yuk=',Asese_ud_yuk,'Fsese_yuk=',Fsese_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3)'), &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=',Asese_yuk,'Asese_ud_yuk=',Asese_ud_yuk,'Fsese_yuk=',Fsese_yuk
					END IF
				ELSE
					IF (split_Fsese) THEN
						PRINT '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3)' , &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=', Asese_yuk, 'Fsese_yuk=', Fsese_yuk, 'Fsese_ud_yuk=', Fsese_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A17,F9.3,5X,A13,F9.3,5X,A13,F9.3)'), &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=', Asese_yuk, 'Fsese_yuk=', Fsese_yuk, 'Fsese_ud_yuk=', Fsese_ud_yuk
					ELSE
						PRINT '(6X,A5,A3,A17,F9.3,5X,A13,F9.3)' , &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=', Asese_yuk, 'Fsese_yuk=', Fsese_yuk
						IF (flag_output) WRITE (7, '(6X,A5,A3,A17,F9.3,5X,A13,F9.3)'), &
						  'Jse: ',Jse_kind,'   -   Asese_yuk=', Asese_yuk, 'Fsese_yuk=', Fsese_yuk
					END IF
				END IF
			END SELECT
		END IF
		
		IF (Jsesp_kind/='no_') THEN
			SELECT CASE (Jsesp_kind)
			CASE ('pot')
				PRINT '(6X,A7,A3,A14,F9.3)' , &
				  'Jsesp: ',Jsesp_kind,'   -   c_sesp=', c_sesp
				IF (flag_output) WRITE (7, '(6X,A7,A3,A14,F9.3)'), &
				  'Jsesp: ',Jsesp_kind,'   -   c_sesp=', c_sesp
			CASE ('yuk')
				IF (split_Asesp) THEN
					IF (split_Fsesp) THEN
						PRINT '(6X,A7,A3,A17,F9.3,5X,A13,F9.3,5X,A10,F9.3,5X,A13,F9.3)' , &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Asesp_ud_yuk=', Asesp_ud_yuk,&
						  'Fsesp_yuk=', Fsesp_yuk, 'Fsesp_ud_yuk=', Fsesp_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A7,A3,A17,F9.3,5X,A13,F9.3,5X,A10,F9.3,5X,A13,F9.3)'), &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Asesp_ud_yuk=', Asesp_ud_yuk,&
						  'Fsesp_yuk=', Fsesp_yuk, 'Fsesp_ud_yuk=', Fsesp_ud_yuk
					ELSE
						PRINT '(6X,A7,A3,A17,F9.3,5X,A13,F9.3,5X,A10,F9.3)' , &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Asesp_ud_yuk=', Asesp_ud_yuk, 'Fsesp_yuk=', Fsesp_yuk
						IF (flag_output) WRITE (7, '(6X,A7,A3,A17,F9.3,5X,A13,F9.3,5X,A10,F9.3)'), &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Asesp_ud_yuk=', Asesp_ud_yuk, 'Fsesp_yuk=', Fsesp_yuk
					END IF
				ELSE
					IF (split_Fsesp) THEN
						PRINT '(6X,A7,A3,A17,F9.3,5X,A10,F9.3,5X,A13,F9.3)' , &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Fsesp_yuk=', Fsesp_yuk, 'Fsesp_ud_yuk=', Fsesp_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A7,A3,A17,F9.3,5X,A10,F9.3)'), &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Fsesp_yuk=', Fsesp_yuk, 'Fsesp_ud_yuk=', Fsesp_ud_yuk
					ELSE
						PRINT '(6X,A7,A3,A17,F9.3,5X,A10,F9.3)' , &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Fsesp_yuk=', Fsesp_yuk
						IF (flag_output) WRITE (7, '(6X,A7,A3,A17,F9.3,5X,A10,F9.3,5X,A13,F9.3)'), &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Fsesp_yuk=', Fsesp_yuk
					END IF
				END IF
			CASE ('yup')
				IF (split_Asesp) THEN
					IF (split_Fsesp) THEN
						PRINT '(6X,A7,A3,A17,F9.3,5X,A13,F9.3,5X,A10,F9.3,5X,A13,F9.3)' , &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Asesp_ud_yuk=', Asesp_ud_yuk,&
						  'Fsesp_yuk=', Fsesp_yuk, 'Fsesp_ud_yuk=', Fsesp_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A7,A3,A17,F9.3,5X,A13,F9.3,5X,A10,F9.3,5X,A13,F9.3)'), &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Asesp_ud_yuk=', Asesp_ud_yuk,&
						  'Fsesp_yuk=', Fsesp_yuk, 'Fsesp_ud_yuk=', Fsesp_ud_yuk
					ELSE
						PRINT '(6X,A7,A3,A17,F9.3,5X,A13,F9.3,5X,A10,F9.3)' , &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Asesp_ud_yuk=', Asesp_ud_yuk, 'Fsesp_yuk=', Fsesp_yuk
						IF (flag_output) WRITE (7, '(6X,A7,A3,A17,F9.3,5X,A13,F9.3,5X,A10,F9.3)'), &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Asesp_ud_yuk=', Asesp_ud_yuk, 'Fsesp_yuk=', Fsesp_yuk
					END IF
				ELSE
					IF (split_Fsesp) THEN
						PRINT '(6X,A7,A3,A17,F9.3,5X,A10,F9.3,5X,A13,F9.3)' , &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Fsesp_yuk=', Fsesp_yuk, 'Fsesp_ud_yuk=', Fsesp_ud_yuk
						IF (flag_output) WRITE (7, '(6X,A7,A3,A17,F9.3,5X,A10,F9.3)'), &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Fsesp_yuk=', Fsesp_yuk, 'Fsesp_ud_yuk=', Fsesp_ud_yuk
					ELSE
						PRINT '(6X,A7,A3,A17,F9.3,5X,A10,F9.3)' , &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Fsesp_yuk=', Fsesp_yuk
						IF (flag_output) WRITE (7, '(6X,A7,A3,A17,F9.3,5X,A10,F9.3,5X,A13,F9.3)'), &
						  'Jsesp: ',Jsesp_kind,'   -   Asesp_yuk=', Asesp_yuk, 'Fsesp_yuk=', Fsesp_yuk
					END IF
				END IF
			CASE ('gss')
				PRINT '(6X,A7,A3,A13,F9.3)' , &
				  'Jsesp: ',Jsesp_kind,'   -   Gsesp=', Gsesp
				IF (flag_output) WRITE (7, '(6X,A7,A3,A13,F9.3)'), &
				  'Jsesp: ',Jsesp_kind,'   -   Gsesp=', Gsesp
			CASE ('gsd')
				PRINT '(6X,A7,A3,A13,F9.3)' , &
				  'Jsesp: ',Jsesp_kind,'   -   Gsesp=', Gsesp
				IF (flag_output) WRITE (7, '(6X,A7,A3,A13,F9.3)'), &
				  'Jsesp: ',Jsesp_kind,'   -   Gsesp=', Gsesp
			END SELECT
		END IF
		
		IF (SDse_kind=='gem') THEN
			PRINT '(6X,A6,A3,A10,F9.3)' , 'SDse: ', SDse_kind,'  -  Gswf=', Gswf
			IF (flag_output) WRITE (7, '(6X,A6,A3,A10,F9.3)'), &
			  'SDse: ', SDse_kind,'  -  Gswf=', Gswf
		ELSE IF (SDse_kind=='gss') THEN
			PRINT '(6X,A6,A3,A10,F9.3)' , 'SDse: ', SDse_kind,'  -  Gswf=', Gswf
			IF (flag_output) WRITE (7, '(6X,A6,A3,A10,F9.3)'), &
			  'SDse: ', SDse_kind,'  -  Gswf=', Gswf
		ELSE IF (SDse_kind=='gsp') THEN
			PRINT '(6X,A6,A3,A10,F9.3)' , 'SDse: ', SDse_kind,'  -  Gswf=', Gswf
			IF (flag_output) WRITE (7, '(6X,A6,A3,A10,F9.3)'), &
			  'SDse: ', SDse_kind,'  -  Gswf=', Gswf
		ELSE IF ((SDse_kind=='atm').OR.(SDse_kind=='atp')) THEN
			PRINT '(6X,A6,A3,A11,F9.3)' , 'SDse: ', SDse_kind,'  -  C_atm=', C_atm
			IF (flag_output) WRITE (7, '(6X,A6,A3,A11,F9.3)'), &
			  'SDse: ', SDse_kind,'  -  C_atm=', C_atm
		ELSE IF (SDse_kind/='no_') THEN
			PRINT '(6X,A6,A3)' , 'SDse: ', SDse_kind
			IF (flag_output) WRITE (7, '(6X,A6,A3)'), 'SDse: ', SDse_kind
		END IF
	END SUBROUTINE stampa_parametri_variazionali
!-----------------------------------------------------------------------

	SUBROUTINE genera_orbitali_lda()
		USE dft
		USE dati_fisici
		USE dati_mc
		IMPLICIT NONE
		
		IF ( mpi_myrank==0 ) THEN
			WRITE (*, FMT='(A42)', ADVANCE='NO'), 'Genero orbitali DFT con Quantum-Espresso '
			CALL inizializza_dft()
			WRITE (*, FMT='(A1)', ADVANCE='NO'), '.'
			CALL genera_scf_in(L,N_part,r_crystal)
			WRITE (*, FMT='(A1)', ADVANCE='NO'), '.'
			CALL run_quantum_espresso()
			DO WHILE ( dft_lock )
				CALL SLEEP(1)
			END DO
			WRITE (*, FMT='(A1)', ADVANCE='NO'), '.'
			CALL chiudi_dft()
			lda_path=TRIM(dft_work_dir)//"OUT.save/"
			IF (.NOT. flag_TABC) CALL change_flag_TABC(.TRUE.)
			WRITE (*, FMT='(A1)', ADVANCE='YES'), '!'
		END IF
		CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
		
		lda_path=TRIM(dft_work_dir)//"OUT.save"
		IF (.NOT. flag_TABC) CALL change_flag_TABC(.TRUE.)
		
	END SUBROUTINE genera_orbitali_lda
!-----------------------------------------------------------------------

	SUBROUTINE normalizza_coefficienti_lda() !necessaria per non far divergere a infinito (NaN) nel caso HCP con 432 particelle
		USE dati_mc
		USE dati_fisici
		USE walkers
		IMPLICIT NONE
		INTEGER :: i, j, ik, info
		INTEGER :: pvt(1:H_N_part)
		REAL (KIND=8) :: norm_fact_up, norm_fact_dw, shared_norm_fact_up, shared_norm_fact_dw
		COMPLEX (KIND=8) :: SD(1:H_N_part,1:H_N_part), ISD(1:H_N_part,1:H_N_part)
		
		IF (.NOT. iniz_walkers) STOP 'Non puoi normalizzare i coefficienti lda prima di aver inizializzato i walkers &
		  [ module_funzione_onda.f90 > normalizza_coefficienti_lda ]'				
		IF (.NOT. iniz_dati_fisici) STOP 'Non puoi normalizzare i coefficienti lda prima di aver letto i dati fisici &
		  [ module_funzione_onda.f90 > normalizza_coefficienti_lda ]'
		
		IF ( flag_simm_lda ) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					SD(i,j)=0.d0
					DO ik = 1, num_fpw_lda(j), 1
						SD(i,j)=SD(i,j)+fattori_fpw_lda(ik,j)*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),re_new(1:3,i))) + &
						                DCONJG(fattori_fpw_lda(ik,j))*CDEXP(-(0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),re_new(1:3,i)))
					END DO
					ISD(i,j)=SD(i,j)
				END DO
			END DO
		ELSE
			!Calcolo i termini matriciali di SD_new
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					SD(i,j)=0.d0
					DO ik = 1, num_fpw_lda(j), 1
						SD(i,j)=SD(i,j)+fattori_fpw_lda(ik,j)*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),re_new(1:3,i)))
					END DO
					ISD(i,j)=SD(i,j)
				END DO
			END DO
		END IF
		!Calcolo il determinante di SD_new
		CALL ZGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		norm_fact_up=0.d0
		DO i = 1, H_N_part, 1
			!norm_fact_up=norm_fact_up+DSQRT(REAL(ISD(i,i)*DCONJG(ISD(i,i)),8))
			!norm_fact_up=norm_fact_up+DABS(REAL(ISD(i,i),8))
			norm_fact_up=norm_fact_up+(DABS(REAL(ISD(i,i),8))+DABS(AIMAG(ISD(i,i))))*0.5d0
		END DO
		norm_fact_up=norm_fact_up/REAL(H_N_part,8)
		shared_norm_fact_up=0.d0
		CALL MPI_REDUCE(norm_fact_up,shared_norm_fact_up,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
		CALL MPI_BCAST(shared_norm_fact_up, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_ierr)
		shared_norm_fact_up=shared_norm_fact_up/REAL(mpi_nprocs,8)
		
		IF ( flag_simm_lda ) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					SD(i,j)=0.d0
					DO ik = 1, num_fpw_lda(j), 1
						SD(i,j)=SD(i,j)+fattori_fpw_lda(ik,j)*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),re_new(1:3,i+H_N_part))) + &
						                DCONJG(fattori_fpw_lda(ik,j))*CDEXP(-(0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),re_new(1:3,i+H_N_part)))
					END DO
					ISD(i,j)=SD(i,j)
				END DO
			END DO
		ELSE
			!Calcolo i termini matriciali di SD_new
			DO j = 1, H_N_part, 1
				DO i = 1, H_N_part, 1
					SD(i,j)=0.d0
					DO ik = 1, num_fpw_lda(j), 1
						SD(i,j)=SD(i,j)+fattori_fpw_lda(ik,j)*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),re_new(1:3,i+H_N_part)))
					END DO
					ISD(i,j)=SD(i,j)
				END DO
			END DO
		END IF
		!Calcolo il determinante di SD_new
		CALL ZGETRF( H_N_part, H_N_part, ISD, H_N_part, pvt, info )
		IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
		norm_fact_dw=0.d0
		DO i = 1, H_N_part, 1                                                                 !gli altri due sono formule alternative, a cui ho preferito l'ultima. Gli SD riportati sono quelli ottenibili con HCP-432, ctf=10
			!norm_fact_dw=norm_fact_dw+DSQRT(REAL(ISD(i,i)*DCONJG(ISD(i,i)),8))               !DET= (-2.60030859368052841E-030, 3.52638097693409618E-031)
			!norm_fact_dw=norm_fact_dw+DABS(REAL(ISD(i,i),8))                                 !DET= ( -293469390.51675814     ,  39798540.778806597     )
			norm_fact_dw=norm_fact_dw+(DABS(REAL(ISD(i,i),8))+DABS(AIMAG(ISD(i,i))))*0.5d0   !DET= ( -342130869.60716617     ,  46397715.760995239     )
		END DO
		norm_fact_dw=norm_fact_dw/REAL(H_N_part,8)
		shared_norm_fact_dw=0.d0
		CALL MPI_REDUCE(norm_fact_dw,shared_norm_fact_dw,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
		CALL MPI_BCAST(shared_norm_fact_dw, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_ierr)
		shared_norm_fact_dw=shared_norm_fact_dw/REAL(mpi_nprocs,8)
				
		fattori_fpw_lda=fattori_fpw_lda/(0.5d0*(shared_norm_fact_dw+shared_norm_fact_dw))
		
	END SUBROUTINE normalizza_coefficienti_lda
!-----------------------------------------------------------------------

	SUBROUTINE applica_twist_dnfH()
		USE dnfH
		IMPLICIT NONE
		INTEGER :: i, j, max_num_pw
				
		CALL inizializza_dnfH(kf_coeff_dnfH)     !inizializzo dnfH, trovano N_pw_lda
		CALL twist_k_dnfH(kf_coeff_dnfH)
		IF (flag_usa_coeff_dnfH) THEN                        !se voglio usare i coefficienti nella matrice dnfH, li carico
			CALL costruisci_matrice_Hartree(SDe_kind,c_eff_dnfH)
		ELSE
			CALL costruisci_matrice_Hartree(SDe_kind)
		END IF
		CALL trova_soluzioni_dnfH()
		N_pw_lda=N_M_Hartree
		DEALLOCATE(k_pw_lda,num_pw_orbit,fattori_pw_lda,autoenergie_dnfH,num_pw_dnfH)
		ALLOCATE(k_pw_lda(0:3,1:N_pw_lda),num_pw_orbit(1:N_pw_lda),fattori_pw_lda(1:N_pw_lda,1:N_pw_lda))
		ALLOCATE(autoenergie_dnfH(1:N_pw_lda),num_pw_dnfH(1:H_N_part))
		k_pw_lda(0:3,1:N_pw_lda)=k_Hartree
		autoenergie_dnfH=autovalori_Hartree
		num_pw_orbit=n_fpw_Hartee
		num_pw_dnfH(1:H_N_part)=num_pw_orbit(1:H_N_part)
		max_num_pw=0
		DO j = 1, N_pw_lda, 1
			IF (num_pw_orbit(j)>max_num_pw) max_num_pw=num_pw_orbit(j)
		END DO
		DEALLOCATE(indice_pw_dnfH)
		ALLOCATE(indice_pw_dnfH(1:max_num_pw,1:N_pw_lda))
		indice_pw_dnfH(1:max_num_pw,1:N_pw_lda)=i_fpw_Hartree(1:max_num_pw,1:N_pw_lda)
		fattori_pw_lda=M_Hartree
		DEALLOCATE(k_pw_dnfH)
		ALLOCATE(k_pw_dnfH(0:3,1:max_num_pw,1:N_pw_lda))
		DO j = 1, N_pw_lda, 1
			DO i = 1, num_pw_orbit(j), 1
				k_pw_dnfH(0:3,i,j)=k_pw_lda(0:3,indice_pw_dnfH(i,j))
			END DO
		END DO
		DEALLOCATE(fattori_pw_dnfH)
		ALLOCATE(fattori_pw_dnfH(1:max_num_pw,1:N_pw_lda))
		DO j = 1, N_pw_lda, 1
			DO i = 1, num_pw_orbit(j), 1
				fattori_pw_dnfH(i,j)=fattori_pw_lda(indice_pw_dnfH(i,j),j)
			END DO
		END DO
		CALL chiudi_dnfH()
	END SUBROUTINE applica_twist_dnfH
!-----------------------------------------------------------------------

	SUBROUTINE applica_twist_lda()
		USE dati_mc
		USE dati_fisici
		USE generic_tools
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
		INTEGER :: i, j, i1, max_num_pw, i_twist
		CHARACTER(LEN=4) :: istring
		REAL (KIND=8) :: dummy1
		
		CALL RANDOM_NUMBER(dummy1)
		i_twist=1
		DO WHILE ( dummy1>pesi_K_points(i_twist+1) )
			i_twist=i_twist+1
		END DO
		IF (i_twist>num_K_points) STOP 'Errore: i_twist é maggiore di num_K_points &
		  [ module_funzione_onda.f90 > inizializza_dati_funzione_onda ]'
		WRITE (istring, '(I4.4)'), i_twist
		
		!num_chiamata_twist_lda=num_chiamata_twist_lda+1
		!WRITE (istring, '(I4.4)'), MOD(mpi_nprocs*num_chiamata_twist_lda+mpi_myrank,260)+1
		
		CALL leggi_N_pw(TRIM(lda_path)//'/K0'//istring//'/gkvectors.xml',N_pw_lda)
		DEALLOCATE(k_pw_lda,fattori_orb_lda,fattori_pw_lda)
		DEALLOCATE(twist_lda)
		ALLOCATE(k_pw_lda(0:3,1:N_pw_lda),fattori_orb_lda(1:H_N_part),fattori_pw_lda(1:N_pw_lda,1:H_N_part))
		ALLOCATE(k_pw_int_lda(1:3,1:N_pw_lda),twist_lda(1:3))
		CALL leggi_evc_xml(H_N_part,N_pw_lda,TRIM(lda_path)//'/K0'//istring//'/evc.xml',fattori_pw_lda)       !o!
		CALL leggi_gkvectors_xml(N_pw_lda,TRIM(lda_path)//'/K0'//istring//'/gkvectors.xml',k_pw_int_lda(1:3,1:N_pw_lda),twist_lda(1:3))
		DO i = 1, 3, 1
			twist_lda(i)=twist_lda(i)/r_s
		END DO
		DO j = 1, N_pw_lda, 1
			DO i = 1, 3, 1
				k_pw_lda(i,j)=(2.d0*PI/L(i))*REAL(k_pw_int_lda(i,j),8)+twist_lda(i)
			END DO
		END DO
		DEALLOCATE(k_pw_int_lda)
		DO i = 1, N_pw_lda, 1
			k_pw_lda(0,i)=DOT_PRODUCT(k_pw_lda(1:3,i),k_pw_lda(1:3,i))
		END DO
		
		!trovo i coefficienti per velocizzare il metodo
		max_num_pw=0
		DO j = 1, H_N_part, 1
			i1=0
			DO i = 1, N_pw_lda, 1
				IF (REAL(fattori_pw_lda(i,j)*DCONJG(fattori_pw_lda(i,j)),4)<CUT_LDA) THEN
					
				ELSE
					i1=i1+1
				END IF
			END DO
			max_num_pw=MAX(max_num_pw,i1)
		END DO
		DEALLOCATE(num_fpw_lda,k_fpw_lda,fattori_fpw_lda)
		ALLOCATE(num_fpw_lda(1:H_N_part),k_fpw_lda(0:3,1:max_num_pw,1:H_N_part),fattori_fpw_lda(1:max_num_pw,1:H_N_part))
		num_fpw_lda=0
		DO j = 1, H_N_part, 1
			i1=0
			DO i = 1, N_pw_lda, 1
				IF (REAL(fattori_pw_lda(i,j)*DCONJG(fattori_pw_lda(i,j)),4)<CUT_LDA) THEN
					
				ELSE
					i1=i1+1
					num_fpw_lda(j)=num_fpw_lda(j)+1
					k_fpw_lda(0:3,num_fpw_lda(j),j)=k_pw_lda(0:3,i)
					fattori_fpw_lda(num_fpw_lda(j),j)=fattori_pw_lda(i,j)
				END IF
			END DO
		END DO
		!CALL normalizza_coefficienti_lda()
		
	END SUBROUTINE applica_twist_lda
!-----------------------------------------------------------------------

	SUBROUTINE valuta_SD_pw(num,r,N,k,SD,detSD,ISD,pvt,ISD_old,detSD_old)
		USE generic_tools
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: r(1:3,1:N), k(1:3,1:N)
		COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
		INTEGER :: i, j, info, perm
		INTEGER, INTENT(OUT) :: pvt(1:N)
		REAL (KIND=8) :: norm
		COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_SD_pw ]'
		
		norm=1.d0/SQRT(REAL(N))
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, N, 1
				DO i = 1, N, 1
					SD(i,j)=norm*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k(1:3,j),r(1:3,i)))
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
			DO i = 1, N, 1
				SD(num,i)=norm*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k(1:3,i),r(1:3,num)))
			END DO
			CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_SD_pw ]'
		END IF
						
		IF (verbose_mode) PRINT * , 'funzione_onda: detSD(pw)=', detSD
	END SUBROUTINE valuta_SD_pw
!-----------------------------------------------------------------------

	SUBROUTINE valuta_SD_gauss(num,rij,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
		USE generic_tools
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(1:N,1:N)
		COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
		INTEGER :: i, j, info, perm
		INTEGER, INTENT(OUT) :: pvt(1:N)
		REAL (KIND=8) :: norm
		COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_SD_pw ]'
		
		norm=2.d0**(1./REAL(N))
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, N, 1
				DO i = 1, N, 1
					SD(i,j)=norm*DEXP(-Gswf*rij(i,j)*rij(i,j))
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
			DO i = 1, N, 1
				SD(num,i)=norm*norm*DEXP(-Gswf*rij(num,i)*rij(num,i))
			END DO
			CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_SD_pw ]'
		END IF
				
		IF (verbose_mode) PRINT * , 'funzione_onda: detSD(gss)=', detSD
				
		!PRINT * , 'funzione_onda: detSD(gss)=', detSD_old, detSD
		!
		!STOP
		
	END SUBROUTINE valuta_SD_gauss
!-----------------------------------------------------------------------

	SUBROUTINE valuta_SD_atm(num,rij,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
		USE generic_tools
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(1:N,1:N)
		COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
		INTEGER :: i, j, info, perm
		INTEGER, INTENT(OUT) :: pvt(1:N)
		REAL (KIND=8) :: norm
		COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_SD_atm ]'
		
		norm=1.d0 !/DSQRT(PI)
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, N, 1
				DO i = 1, N, 1
					SD(i,j)=(1.d0,0.d0)*norm*DEXP(-C_atm*rij(i,j))
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
			DO i = 1, N, 1
				SD(num,i)=(1.d0,0.d0)*norm*DEXP(-C_atm*rij(num,i))
			END DO
			CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_SD_atm ]'
		END IF
				
		IF (verbose_mode) PRINT * , 'funzione_onda: detSD(gss)=', detSD
				
	END SUBROUTINE valuta_SD_atm
	!-----------------------------------------------------------------------

		SUBROUTINE valuta_SD_bat(num,updw,rij,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
			USE generic_tools
			IMPLICIT NONE
			CHARACTER (LEN=2) :: updw
			INTEGER, INTENT(IN) :: N, num
			REAL (KIND=8), INTENT(IN) :: rij(1:N+N,1:N+N)
			COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
			INTEGER :: i, j, info, perm
			INTEGER, INTENT(OUT) :: pvt(1:N)
			REAL (KIND=8) :: norm
			COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
		
			IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
			  [ module_funzione_onda.f90 > valuta_SD_bat ]'
		
			norm=1.d0
			
			SELECT CASE(updw)
			CASE('up')
				IF (num==-1) THEN
					!Calcolo i termini matriciali di SD_new
					DO j = 1, N, 1
						DO i = 1, N, 1
							SD(i,j)=(1.d0,0.d0)*norm*( DEXP(-C_atm*rij(i,j)) + DEXP(-C_atm*rij(i,j+N)) )
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
					DO i = 1, N, 1
						SD(num,i)=(1.d0,0.d0)*norm*( DEXP(-C_atm*rij(num,i)) + DEXP(-C_atm*rij(num,i+N)) )
					END DO
					CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
				ELSE
					STOP 'num non accettabile &
					  [ module_funzione_onda.f90 > valuta_SD_bat ]'
				END IF
			CASE('dw')
				IF (num==-1) THEN
					!Calcolo i termini matriciali di SD_new
					DO j = 1, N, 1
						DO i = 1, N, 1
							SD(i,j)=(1.d0,0.d0)*norm*( DEXP(-C_atm*rij(i+N,j+N)) + DEXP(-C_atm*rij(i+N,j)) )
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
					DO i = 1, N, 1
						SD(num,i)=(1.d0,0.d0)*norm*( DEXP(-C_atm*rij(num+N,i+N)) + DEXP(-C_atm*rij(num+N,i)) )
					END DO
					CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
				ELSE
					STOP 'num non accettabile &
					  [ module_funzione_onda.f90 > valuta_SD_bat ]'
				END IF
				CASE DEFAULT
					STOP 'serve specificare up o dw [ module_funzione_onda.f90 > valuta_SD_bat ]'
			END SELECT
				
			IF (verbose_mode) PRINT * , 'funzione_onda: detSD(gss)=', detSD
				
		END SUBROUTINE valuta_SD_bat	
	
	!-----------------------------------------------------------------------

		SUBROUTINE valuta_SD_HL(rij,SD,detSD,ISD)
			USE generic_tools
			IMPLICIT NONE
			REAL (KIND=8), INTENT(IN) :: rij(1:2,1:2)
			COMPLEX (KIND=8) :: SD(1:1,1:1), ISD(1:1,1:1), detSD
		
			IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
			  [ module_funzione_onda.f90 > valuta_SD_HL ]'
		
         SD(1,1)=(1.d0,0.d0)*(DEXP(-C_atm*(rij(1,1)+rij(2,2)))+DEXP(-C_atm*(rij(1,2)+rij(2,1))))
         detSD=SD(1,1)
         ISD(1,1)=(1.d0,0.d0)/detSD

			IF (verbose_mode) PRINT * , 'funzione_onda: detSD(gss)=', detSD
				
		END SUBROUTINE valuta_SD_HL	
	
	!-----------------------------------------------------------------------

		SUBROUTINE valuta_SD_1s_backflow(num,updw,L,re,rp,rij,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
			USE generic_tools
			IMPLICIT NONE
			REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
         CHARACTER (LEN=2) :: updw
			INTEGER, INTENT(IN) :: N, num   !num e' compreso fra 1 e 2*N
         REAL(KIND=8), INTENT(IN) :: L(1:3)
			REAL (KIND=8), INTENT(IN) :: re(1:3,1:N+N), rp(1:3,1:N+N), rij(1:N+N,1:N+N)
			COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
			INTEGER :: i, j, ip, info, perm, iadd
			INTEGER, INTENT(OUT) :: pvt(1:N)
			REAL (KIND=8) :: norm
			COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
         COMPLEX (KIND=8) :: ISD_now(1:N,1:N), detSD_now
         REAL(KIND=8) :: q(0:3), dist(0:3)
		
			IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
			  [ module_funzione_onda.f90 > valuta_SD_1s_backflow ]'
		
			norm=1.d0 !/DSQRT(PI)
			
         IF (updw=='up') iadd=0
         IF (updw=='dw') iadd=N

			IF (num==-1) THEN
				!Calcolo i termini matriciali di SD_new
				DO j = 1, N, 1
					DO i = 1, N, 1
                  q(1:3)=re(1:3,i+iadd)
                  DO ip = 1, N+N, 1
                  	q(1:3)=q(1:3)-rp(1:3,ip)/(1.d0+DEXP(A_POT_se*(rij(j+iadd,ip)-D_POT_se)))
                  END DO
                  q(1:3)=q(1:3)-L(1:3)*DNINT(q(1:3)/L(1:3))
                  q(0)=DSQRT(DOT_PRODUCT(q(1:3),q(1:3)))
						SD(i,j)=(1.d0,0.d0)*norm*DEXP(-C_atm*q(0))
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
			ELSE IF ((num>0) .AND. (num<=N+N)) THEN
				DO i = 1+iadd, N+iadd, 1
               q(1:3)=re(1:3,num+iadd)
               DO ip = 1, N+N, 1
                  q(1:3)=q(1:3)-rp(1:3,ip)/(1.d0+DEXP(A_POT_se*(rij(i,ip)-D_POT_se)))
               END DO
               q(1:3)=q(1:3)-L(1:3)*DNINT(q(1:3)/L(1:3))
               q(0)=DSQRT(DOT_PRODUCT(q(1:3),q(1:3)))
					SD(num,i-iadd)=(1.d0,0.d0)*norm*DEXP(-C_atm*q(0))
				END DO
				CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD_now)  !aggiorna det
            CALL aggiorna_matrice_inversa_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD_now,ISD_now) !aggiorna ISD
				DO i = 1+iadd, N+iadd, 1
               q(1:3)=re(1:3,i)
               DO ip = 1, N+N, 1
                  q(1:3)=q(1:3)-rp(1:3,ip)/(1.d0+DEXP(A_POT_se*(rij(num+iadd,ip)-D_POT_se)))
               END DO
               q(1:3)=q(1:3)-L(1:3)*DNINT(q(1:3)/L(1:3))
               q(0)=DSQRT(DOT_PRODUCT(q(1:3),q(1:3)))
					SD(i-iadd,num)=(1.d0,0.d0)*norm*DEXP(-C_atm*q(0))
				END DO
            CALL aggiorna_determinante_C_col_1ppt(N,num,ISD_now,detSD_now,SD,detSD)
			ELSE
				STOP 'num non accettabile &
				  [ module_funzione_onda.f90 > valuta_SD_1s_backflow ]'
			END IF
				
			IF (verbose_mode) PRINT * , 'funzione_onda: detSD(gss)=', detSD
				
		END SUBROUTINE valuta_SD_1s_backflow	
	
!-----------------------------------------------------------------------

	SUBROUTINE valuta_SD_spl_backflow(num,updw,L,re,rp,rij,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
		USE generic_tools
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
      CHARACTER (LEN=2) :: updw
		INTEGER, INTENT(IN) :: N, num   !num e' compreso fra 1 e 2*N
      REAL(KIND=8), INTENT(IN) :: L(1:3)
		REAL (KIND=8), INTENT(IN) :: re(1:3,1:N+N), rp(1:3,1:N+N), rij(1:N+N,1:N+N)
		COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
		INTEGER :: i, j, ip, info, perm, iadd
		INTEGER, INTENT(OUT) :: pvt(1:N)
		REAL (KIND=8) :: norm
		COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
      COMPLEX (KIND=8) :: ISD_now(1:N,1:N), detSD_now
      REAL(KIND=8) :: q(0:3), dist(0:3), cspl
	
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_SD_1s_backflow ]'
	
		norm=1.d0 !/DSQRT(PI)
		
      IF (updw=='up') iadd=0
      IF (updw=='dw') iadd=N

		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, N, 1
				DO i = 1, N, 1
               q(1:3)=re(1:3,i+iadd)
               DO ip = 1, N+N, 1
               	q(1:3)=q(1:3)-rp(1:3,ip)/(1.d0+DEXP(A_POT_se*(rij(j+iadd,ip)-D_POT_se)))
               END DO
               q(1:3)=q(1:3)-L(1:3)*DNINT(q(1:3)/L(1:3))
               q(0)=DSQRT(DOT_PRODUCT(q(1:3),q(1:3)))
               CALL MSPL_compute(SPL=Bsplep, DERIV=0, R=q(0), VAL=cspl)
					SD(i,j)=(1.d0,0.d0)*norm*DEXP(-cspl)
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
		ELSE IF ((num>0) .AND. (num<=N+N)) THEN
			DO i = 1+iadd, N+iadd, 1
            q(1:3)=re(1:3,num+iadd)
            DO ip = 1, N+N, 1
               q(1:3)=q(1:3)-rp(1:3,ip)/(1.d0+DEXP(A_POT_se*(rij(i,ip)-D_POT_se)))
            END DO
            q(1:3)=q(1:3)-L(1:3)*DNINT(q(1:3)/L(1:3))
            q(0)=DSQRT(DOT_PRODUCT(q(1:3),q(1:3)))
            CALL MSPL_compute(SPL=Bsplep, DERIV=0, R=q(0), VAL=cspl)
				SD(num,i-iadd)=(1.d0,0.d0)*norm*DEXP(-cspl)
			END DO
			CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD_now)  !aggiorna det
         CALL aggiorna_matrice_inversa_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD_now,ISD_now) !aggiorna ISD
			DO i = 1+iadd, N+iadd, 1
            q(1:3)=re(1:3,i)
            DO ip = 1, N+N, 1
               q(1:3)=q(1:3)-rp(1:3,ip)/(1.d0+DEXP(A_POT_se*(rij(num+iadd,ip)-D_POT_se)))
            END DO
            q(1:3)=q(1:3)-L(1:3)*DNINT(q(1:3)/L(1:3))
            q(0)=DSQRT(DOT_PRODUCT(q(1:3),q(1:3)))
            CALL MSPL_compute(SPL=Bsplep, DERIV=0, R=q(0), VAL=cspl)
				SD(i-iadd,num)=(1.d0,0.d0)*norm*DEXP(-cspl)
			END DO
         CALL aggiorna_determinante_C_col_1ppt(N,num,ISD_now,detSD_now,SD,detSD)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_SD_1s_backflow ]'
		END IF
			
		IF (verbose_mode) PRINT * , 'funzione_onda: detSD(gss)=', detSD
			
	END SUBROUTINE valuta_SD_spl_backflow	

!-----------------------------------------------------------------------

	SUBROUTINE valuta_SD_gem(num,rij,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
		USE generic_tools
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(1:2*N,1:2*N)
		COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
		INTEGER :: i, j, info, perm
		INTEGER, INTENT(OUT) :: pvt(1:N)
		REAL (KIND=8) :: norm
		COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_SD_pw ]'
		
		norm=2.d0**(1./REAL(N))
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, N, 1
				DO i = 1, N, 1
					SD(i,j)=norm*DEXP(-Gswf*rij(i,j+N)*rij(i,j+N))
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
		ELSE IF ((num>0) .AND. (num<=2*N)) THEN
			IF ( num<=N ) THEN
				DO j = 1, N, 1
					SD(num,j)=norm*DEXP(-Gswf*rij(num,j+N)*rij(num,j+N))
				END DO
				CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
			ELSE
				DO i = 1, N, 1
					SD(i,num-N)=norm*DEXP(-Gswf*rij(i,num)*rij(i,num))
				END DO
				CALL aggiorna_determinante_C_col_1ppt(N,num-N,ISD_old,detSD_old,SD,detSD)
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_SD_gem ]'
		END IF
				
		IF (verbose_mode) PRINT * , 'funzione_onda: detSD(gem)=', detSD
				
	END SUBROUTINE valuta_SD_gem
!-----------------------------------------------------------------------

	SUBROUTINE valuta_SD_lda(num,r,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		LOGICAL, SAVE :: flag_normalizza_coeff_lda=.FALSE.
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: r(1:3,1:N)
		COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
		INTEGER :: i, j, ik, info, perm
		INTEGER, INTENT(OUT) :: pvt(1:N)
		COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD, norm_fact
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_SD_lda ]'
		
		IF ( flag_simm_lda ) THEN
			IF (num==-1) THEN
				!Calcolo i termini matriciali di SD_new
				DO j = 1, N, 1
					DO i = 1, N, 1
						SD(i,j)=0.d0
						DO ik = 1, num_fpw_lda(j), 1
							SD(i,j)=SD(i,j)+fattori_fpw_lda(ik,j)*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),r(1:3,i))) + &
							                DCONJG(fattori_fpw_lda(ik,j))*CDEXP(-(0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),r(1:3,i)))
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
					DO ik = 1, num_fpw_lda(j), 1
						SD(num,j)=SD(num,j)+fattori_fpw_lda(ik,j)*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),r(1:3,num))) + &
						                    DCONJG(fattori_fpw_lda(ik,j))*CDEXP(-(0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),r(1:3,num)))
					END DO
				END DO
				CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
			ELSE
				STOP 'num non accettabile &
				  [ module_funzione_onda.f90 > valuta_SD_pw ]'
			END IF
		ELSE
			IF (num==-1) THEN
				!Calcolo i termini matriciali di SD_new
				DO j = 1, N, 1
					DO i = 1, N, 1
						SD(i,j)=0.d0
						DO ik = 1, num_fpw_lda(j), 1
							SD(i,j)=SD(i,j)+fattori_fpw_lda(ik,j)*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),r(1:3,i)))
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
					DO ik = 1, num_fpw_lda(j), 1
						SD(num,j)=SD(num,j)+fattori_fpw_lda(ik,j)*CDEXP((0.d0,1.d0)*DOT_PRODUCT(k_fpw_lda(1:3,ik,j),r(1:3,num)))
					END DO
				END DO
				CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
			ELSE
				STOP 'num non accettabile &
				  [ module_funzione_onda.f90 > valuta_SD_pw ]'
			END IF
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: detSD(pw)=', detSD
	END SUBROUTINE valuta_SD_lda
!-----------------------------------------------------------------------
	SUBROUTINE valuta_SD_har(num,r,N,SD,detSD,ISD,pvt,ISD_old,detSD_old)
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: r(1:3,1:N)
		COMPLEX (KIND=8), INTENT(IN) :: ISD_old(1:N,1:N), detSD_old
		INTEGER :: i, j, ik, info, perm
		INTEGER, INTENT(OUT) :: pvt(1:N)
		COMPLEX (KIND=8) :: SD(1:N,1:N), ISD(1:N,1:N), detSD
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_SD_har ]'
		
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			 DO j = 1, N, 1
			 	DO i = 1, N, 1
			 		SD(i,j)=0.d0
			 		DO ik = 1, num_pw_dnfH(j), 1
			 			SD(i,j)=SD(i,j)+fattori_pw_dnfH(ik,j)*CDEXP((0.d0,1.d0)* &
			 			  DOT_PRODUCT(k_pw_dnfH(1:3,ik,j),r(1:3,i)))
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
				DO ik = 1, num_pw_dnfH(j), 1
					SD(num,j)=SD(num,j)+fattori_pw_dnfH(ik,j)*CDEXP((0.d0,1.d0)* &
					  DOT_PRODUCT(k_pw_dnfH(1:3,ik,j),r(1:3,num)))
				END DO
			END DO
			CALL aggiorna_determinante_C_1ppt(N,num,ISD_old,detSD_old,SD,detSD)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_SD_pw ]'
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: detSD(pw)=', detSD
	END SUBROUTINE valuta_SD_har
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Uee_YUK(num,rij,N,u_ee,Uee)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		INTEGER :: i, j, H_N
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		REAL (KIND=8) :: u_ee(1:N,1:N), Uee
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Uee_YUK ]'
		
		H_N=N/2
		IF (num==-1) THEN
			Uee=0.d0
			IF (split_Aee) THEN
				IF (split_Fee) THEN
					DO j = 1, N-1, 1
						DO i = j+1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_ee(i,j)=Aee_yuk*(1.d0-DEXP(-rij(0,i,j)*Fee_yuk))/rij(0,i,j)
							ELSE
								u_ee(i,j)=Aee_ud_yuk*(1.d0-DEXP(-rij(0,i,j)*Fee_ud_yuk))/rij(0,i,j)
							END IF
							Uee=Uee+u_ee(i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N-1, 1
						DO i = j+1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_ee(i,j)=Aee_yuk*(1.d0-DEXP(-rij(0,i,j)*Fee_yuk))/rij(0,i,j)
							ELSE
								u_ee(i,j)=Aee_ud_yuk*(1.d0-DEXP(-rij(0,i,j)*Fee_yuk))/rij(0,i,j)
							END IF
							Uee=Uee+u_ee(i,j)
						END DO
					END DO
				END IF
			ELSE
				IF (split_Fee) THEN
					DO j = 1, N-1, 1
						DO i = j+1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_ee(i,j)=Aee_yuk*(1.d0-DEXP(-rij(0,i,j)*Fee_yuk))/rij(0,i,j)
							ELSE
								u_ee(i,j)=Aee_yuk*(1.d0-DEXP(-rij(0,i,j)*Fee_ud_yuk))/rij(0,i,j)
							END IF
							Uee=Uee+u_ee(i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N-1, 1
						DO i = j+1, N, 1
							u_ee(i,j)=Aee_yuk*(1.d0-DEXP(-rij(0,i,j)*Fee_yuk))/rij(0,i,j)
							Uee=Uee+u_ee(i,j)
						END DO
					END DO
				END IF
			END IF
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (split_Aee) THEN
				IF (split_Fee) THEN
					DO i = num+1, N, 1
						Uee=Uee-u_ee(i,num)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_ee(i,num)=Aee_yuk*(1.d0-DEXP(-rij(0,i,num)*Fee_yuk))/rij(0,i,num)
						ELSE
							u_ee(i,num)=Aee_ud_yuk*(1.d0-DEXP(-rij(0,i,num)*Fee_ud_yuk))/rij(0,i,num)
						END IF
						Uee=Uee+u_ee(i,num)
					END DO
					DO i = 1, num-1, 1
						Uee=Uee-u_ee(num,i)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_ee(num,i)=Aee_yuk*(1.d0-DEXP(-rij(0,num,i)*Fee_yuk))/rij(0,num,i)
						ELSE
							u_ee(num,i)=Aee_ud_yuk*(1.d0-DEXP(-rij(0,num,i)*Fee_ud_yuk))/rij(0,num,i)
						END IF
						Uee=Uee+u_ee(num,i)
					END DO
				ELSE
					DO i = num+1, N, 1
						Uee=Uee-u_ee(i,num)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_ee(i,num)=Aee_yuk*(1.d0-DEXP(-rij(0,i,num)*Fee_yuk))/rij(0,i,num)
						ELSE
							u_ee(i,num)=Aee_ud_yuk*(1.d0-DEXP(-rij(0,i,num)*Fee_yuk))/rij(0,i,num)
						END IF
						Uee=Uee+u_ee(i,num)
					END DO
					DO i = 1, num-1, 1
						Uee=Uee-u_ee(num,i)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_ee(num,i)=Aee_yuk*(1.d0-DEXP(-rij(0,num,i)*Fee_yuk))/rij(0,num,i)
						ELSE
							u_ee(num,i)=Aee_ud_yuk*(1.d0-DEXP(-rij(0,num,i)*Fee_yuk))/rij(0,num,i)
						END IF
						Uee=Uee+u_ee(num,i)
					END DO
				END IF
			ELSE
				IF (split_Fee) THEN
					DO i = num+1, N, 1
						Uee=Uee-u_ee(i,num)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_ee(i,num)=Aee_yuk*(1.d0-DEXP(-rij(0,i,num)*Fee_yuk))/rij(0,i,num)
						ELSE
							u_ee(i,num)=Aee_yuk*(1.d0-DEXP(-rij(0,i,num)*Fee_ud_yuk))/rij(0,i,num)
						END IF
						Uee=Uee+u_ee(i,num)
					END DO
					DO i = 1, num-1, 1
						Uee=Uee-u_ee(num,i)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_ee(num,i)=Aee_yuk*(1.d0-DEXP(-rij(0,num,i)*Fee_yuk))/rij(0,num,i)
						ELSE
							u_ee(num,i)=Aee_yuk*(1.d0-DEXP(-rij(0,num,i)*Fee_ud_yuk))/rij(0,num,i)
						END IF
						Uee=Uee+u_ee(num,i)
					END DO
				ELSE
					DO i = num+1, N, 1
						Uee=Uee-u_ee(i,num)
						u_ee(i,num)=Aee_yuk*(1.d0-DEXP(-rij(0,i,num)*Fee_yuk))/rij(0,i,num)
						Uee=Uee+u_ee(i,num)
					END DO
					DO i = 1, num-1, 1
						Uee=Uee-u_ee(num,i)
						u_ee(num,i)=Aee_yuk*(1.d0-DEXP(-rij(0,num,i)*Fee_yuk))/rij(0,num,i)
						Uee=Uee+u_ee(num,i)
					END DO
				END IF
			END IF
			
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Uee_YUK ]'
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: Uee(yuk)=', Uee
	END SUBROUTINE valuta_Uee_YUK
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Uee_SPL(num,rij,N,u_ee,Uee)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		INTEGER :: i, j, H_N
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		REAL (KIND=8) :: u_ee(1:N,1:N), Uee
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Uee_YUK ]'
		
		H_N=N/2
		IF (num==-1) THEN
			Uee=0.d0
			IF (split_Aee.OR.split_Fee) THEN
				DO j = 1, N-1, 1
					DO i = j+1, N, 1
						IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
                     CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij(0,i,j), VAL=u_ee(i,j))
						ELSE
                     CALL MSPL_compute(SPL=Jsplee_ud, DERIV=0, R=rij(0,i,j), VAL=u_ee(i,j))
						END IF
						Uee=Uee+u_ee(i,j)
					END DO
				END DO
			ELSE
				DO j = 1, N-1, 1
					DO i = j+1, N, 1
                  CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij(0,i,j), VAL=u_ee(i,j))
						Uee=Uee+u_ee(i,j)
					END DO
				END DO
			END IF
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (split_Aee.OR.split_Fee) THEN
				DO i = num+1, N, 1
					Uee=Uee-u_ee(i,num)
					IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
                  CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij(0,i,num), VAL=u_ee(i,num))
					ELSE
                  CALL MSPL_compute(SPL=Jsplee_ud, DERIV=0, R=rij(0,i,num), VAL=u_ee(i,num))
					END IF
					Uee=Uee+u_ee(i,num)
				END DO
				DO i = 1, num-1, 1
					Uee=Uee-u_ee(num,i)
					IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
                  CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij(0,num,i), VAL=u_ee(num,i))
					ELSE
                  CALL MSPL_compute(SPL=Jsplee_ud, DERIV=0, R=rij(0,num,i), VAL=u_ee(num,i))
					END IF
					Uee=Uee+u_ee(num,i)
				END DO
			ELSE
				DO i = num+1, N, 1
					Uee=Uee-u_ee(i,num)
               CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij(0,i,num), VAL=u_ee(i,num))
					Uee=Uee+u_ee(i,num)
				END DO
				DO i = 1, num-1, 1
					Uee=Uee-u_ee(num,i)
               CALL MSPL_compute(SPL=Jsplee, DERIV=0, R=rij(0,num,i), VAL=u_ee(num,i))
					Uee=Uee+u_ee(num,i)
				END DO
			END IF
			
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Uee_YUK ]'
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: Uee(yuk)=', Uee
	END SUBROUTINE valuta_Uee_SPL
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Uep_YUK(num,eop,rij,N,u_ep,Uep)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, eop    !eop indica se si é mosso l'elettrone (=1) o il protone (=2), nel caso num>0
		INTEGER :: i, j, H_N
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		REAL (KIND=8) :: u_ep(1:N,1:N), Uep
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
		
		H_N=N/2
		IF (num==-1) THEN
			Uep=0.d0
			IF (split_Aep) THEN
				IF (split_Fep) THEN
					DO j = 1, N, 1
						DO i = 1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_ep(i,j)=Aep_yuk*(1.d0-DEXP(-rij(0,i,j)*Fep_yuk))/rij(0,i,j)
							ELSE
								u_ep(i,j)=Aep_ud_yuk*(1.d0-DEXP(-rij(0,i,j)*Fep_ud_yuk))/rij(0,i,j)
							END IF
							Uep=Uep+u_ep(i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N, 1
						DO i = 1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_ep(i,j)=Aep_yuk*(1.d0-DEXP(-rij(0,i,j)*Fep_yuk))/rij(0,i,j)
							ELSE
								u_ep(i,j)=Aep_ud_yuk*(1.d0-DEXP(-rij(0,i,j)*Fep_yuk))/rij(0,i,j)
							END IF
							Uep=Uep+u_ep(i,j)
						END DO
					END DO
				END IF
			ELSE
				IF (split_Fep) THEN
					DO j = 1, N, 1
						DO i = 1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_ep(i,j)=Aep_yuk*(1.d0-DEXP(-rij(0,i,j)*Fep_yuk))/rij(0,i,j)
							ELSE
								u_ep(i,j)=Aep_yuk*(1.d0-DEXP(-rij(0,i,j)*Fep_ud_yuk))/rij(0,i,j)
							END IF
							Uep=Uep+u_ep(i,j)
						END DO
					END DO
				ELSE
					!PRINT * , 'qui1'
					!PRINT * , 'Uep=', Uep, "Aep_yuk=", Aep_yuk, "Fep_yuk=", Fep_yuk
					DO j = 1, N, 1
						DO i = 1, N, 1
							u_ep(i,j)=Aep_yuk*(1.d0-DEXP(-rij(0,i,j)*Fep_yuk))/rij(0,i,j)
							Uep=Uep+u_ep(i,j)
							!PRINT * , i, j, rij(0,i,j), u_ep(i,j)
						END DO
					END DO
					!PRINT * , 'Uep=', Uep
				END IF
			END IF
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (split_Aep) THEN
				IF (split_Fep) THEN
					IF (eop==1) THEN
						DO i = 1, N, 1
							Uep=Uep-u_ep(num,i)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_ep(num,i)=Aep_yuk*(1.d0-DEXP(-rij(0,num,i)*Fep_yuk))/rij(0,num,i)
							ELSE
								u_ep(num,i)=Aep_ud_yuk*(1.d0-DEXP(-rij(0,num,i)*Fep_ud_yuk))/rij(0,num,i)
							END IF
							Uep=Uep+u_ep(num,i)
						END DO
					ELSE IF (eop==2) THEN
						DO i = 1, N, 1
							Uep=Uep-u_ep(i,num)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_ep(i,num)=Aep_yuk*(1.d0-DEXP(-rij(0,i,num)*Fep_yuk))/rij(0,i,num)
							ELSE
								u_ep(i,num)=Aep_ud_yuk*(1.d0-DEXP(-rij(0,i,num)*Fep_ud_yuk))/rij(0,i,num)
							END IF
							Uep=Uep+u_ep(i,num)
						END DO
					ELSE
						STOP 'eop non accettabile &
						  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
					END IF
				ELSE
					IF (eop==1) THEN
						DO i = 1, N, 1
							Uep=Uep-u_ep(num,i)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_ep(num,i)=Aep_yuk*(1.d0-DEXP(-rij(0,num,i)*Fep_yuk))/rij(0,num,i)
							ELSE
								u_ep(num,i)=Aep_ud_yuk*(1.d0-DEXP(-rij(0,num,i)*Fep_yuk))/rij(0,num,i)
							END IF
							Uep=Uep+u_ep(num,i)
						END DO
					ELSE IF (eop==2) THEN
						DO i = 1, N, 1
							Uep=Uep-u_ep(i,num)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_ep(i,num)=Aep_yuk*(1.d0-DEXP(-rij(0,i,num)*Fep_yuk))/rij(0,i,num)
							ELSE
								u_ep(i,num)=Aep_ud_yuk*(1.d0-DEXP(-rij(0,i,num)*Fep_yuk))/rij(0,i,num)
							END IF
							Uep=Uep+u_ep(i,num)
						END DO
					ELSE
						STOP 'eop non accettabile &
						  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
					END IF
				END IF
			ELSE
				IF (split_Fep) THEN
					IF (eop==1) THEN
						DO i = 1, N, 1
							Uep=Uep-u_ep(num,i)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_ep(num,i)=Aep_yuk*(1.d0-DEXP(-rij(0,num,i)*Fep_yuk))/rij(0,num,i)
							ELSE
								u_ep(num,i)=Aep_yuk*(1.d0-DEXP(-rij(0,num,i)*Fep_ud_yuk))/rij(0,num,i)
							END IF
							Uep=Uep+u_ep(num,i)
						END DO
					ELSE IF (eop==2) THEN
						DO i = 1, N, 1
							Uep=Uep-u_ep(i,num)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_ep(i,num)=Aep_yuk*(1.d0-DEXP(-rij(0,i,num)*Fep_yuk))/rij(0,i,num)
							ELSE
								u_ep(i,num)=Aep_yuk*(1.d0-DEXP(-rij(0,i,num)*Fep_ud_yuk))/rij(0,i,num)
							END IF
							Uep=Uep+u_ep(i,num)
						END DO
					ELSE
						STOP 'eop non accettabile &
						  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
					END IF
				ELSE
					IF (eop==1) THEN
						DO i = 1, N, 1
							Uep=Uep-u_ep(num,i)
							u_ep(num,i)=Aep_yuk*(1.d0-DEXP(-rij(0,num,i)*Fep_yuk))/rij(0,num,i)
							Uep=Uep+u_ep(num,i)
						END DO
					ELSE IF (eop==2) THEN
						DO i = 1, N, 1
							Uep=Uep-u_ep(i,num)
							u_ep(i,num)=Aep_yuk*(1.d0-DEXP(-rij(0,i,num)*Fep_yuk))/rij(0,i,num)
							Uep=Uep+u_ep(i,num)
						END DO
					ELSE
						STOP 'eop non accettabile &
						  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
					END IF
				END IF
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: Uep(yuk)=', Uep
	END SUBROUTINE valuta_Uep_YUK
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Uep_SPL(num,rij,N,u_ep,Uep)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		INTEGER :: i, j, H_N
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		REAL (KIND=8) :: u_ep(1:N,1:N), Uep
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'

      !PRINT *, MAXVAL(rij(0,:,:))
      !IF (MAXVAL(rij(0,:,:))>100.) STOP 
	
		H_N=N/2
		IF (num==-1) THEN
			Uep=0.d0
			IF (split_Aep.OR.split_Fep) THEN
				DO j = 1, N, 1
					DO i = 1, N, 1
						IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
                     CALL MSPL_compute(SPL=Jsplep, DERIV=0, R=rij(0,i,j), VAL=u_ep(i,j))
						ELSE
                     CALL MSPL_compute(SPL=Jsplep_ud, DERIV=0, R=rij(0,i,j), VAL=u_ep(i,j))
						END IF
						Uep=Uep+u_ep(i,j)
					END DO
				END DO
			ELSE
				DO j = 1, N, 1
					DO i = 1, N, 1
                  CALL MSPL_compute(SPL=Jsplep, DERIV=0, R=rij(0,i,j), VAL=u_ep(i,j))
						Uep=Uep+u_ep(i,j)
					END DO
				END DO
			END IF
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (split_Aep.OR.split_Fep) THEN
				DO i = 1, N, 1
					Uep=Uep-u_ep(num,i)
					IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
                  CALL MSPL_compute(SPL=Jsplep, DERIV=0, R=rij(0,num,i), VAL=u_ep(num,i))
					ELSE
                  CALL MSPL_compute(SPL=Jsplep_ud, DERIV=0, R=rij(0,num,i), VAL=u_ep(num,i))
					END IF
					Uep=Uep+u_ep(num,i)
				END DO
			ELSE
				DO i = 1, N, 1
					Uep=Uep-u_ep(num,i)
               CALL MSPL_compute(SPL=Jsplep, DERIV=0, R=rij(0,num,i), VAL=u_ep(num,i))
					Uep=Uep+u_ep(num,i)
				END DO
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: Uep(yuk)=', Uep
	END SUBROUTINE valuta_Uep_SPL
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Uep_ATM(num,eop,rij,N,u_ep,Uep)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, eop    !eop indica se si é mosso l'elettrone (=1) o il protone (=2), nel caso num>0
		INTEGER :: i, j, H_N
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		REAL (KIND=8) :: u_ep(1:N,1:N), Uep

		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'

		H_N=N/2
		IF (num==-1) THEN
			Uep=0.d0
			DO i = 1, N, 1
				u_ep(i,i)=Fep_yuk*rij(0,i,i)
				Uep=Uep+u_ep(i,i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (eop==1) THEN
				Uep=Uep-u_ep(num,num)
				u_ep(num,num)=Fep_yuk*rij(0,num,num)
				Uep=Uep+u_ep(num,num)
			ELSE IF (eop==2) THEN
				Uep=Uep-u_ep(num,num)
				u_ep(num,num)=Fep_yuk*rij(0,num,num)
				Uep=Uep+u_ep(num,num)
			ELSE
				STOP 'eop non accettabile &
				  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Uep_YUK ]'
		END IF

		IF (verbose_mode) PRINT * , 'funzione_onda: Uep(atm)=', Uep
	END SUBROUTINE valuta_Uep_ATM	
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Use_POT(i_num,dist_mol_ss,u_ss,Uss)
		USE dati_fisici
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: A1=4.64355, B1=3.00909, A2=0.902629, B2=0.924515
		INTEGER, INTENT(IN) :: i_num
		REAL (KIND=8), INTENT(IN) :: dist_mol_ss(1:H_N_part)
		INTEGER :: i, j
		REAL (KIND=8) :: u_ss(1:H_N_part), Uss
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Use_POT ]'
		
		IF (i_num==-1) THEN
			Uss=0.d0
			DO i = 1, H_N_part, 1
				u_ss(i)=A_POT_se*(A1*DEXP(-D_POT_se*B1*dist_mol_ss(i))-A2*DEXP(-D_POT_se*B2*dist_mol_ss(i)))
				Uss=Uss+u_ss(i)
			END DO
		ELSE IF ((i_num>0) .AND. (i_num<=H_N_part)) THEN
			Uss=Uss-u_ss(i_num)
			u_ss(i_num)=A_POT_se*(A1*DEXP(-D_POT_se*B1*dist_mol_ss(i_num))-A2*DEXP(-D_POT_se*B2*dist_mol_ss(i_num)))
			Uss=Uss+u_ss(i_num)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Use_POT ]'
		END IF
		IF (verbose_mode) PRINT * , 'funzione_onda: Use(pot)=', Uss
				
	END SUBROUTINE valuta_Use_POT
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Usese_YUK(num,sij,N,u_sese,Usese)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		INTEGER :: i, j, H_N
		REAL (KIND=8), INTENT(IN) :: sij(0:3,1:N,1:N)
		REAL (KIND=8) :: u_sese(1:N,1:N), Usese
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Usese_YUK ]'
		
		H_N=N/2
		IF (num==-1) THEN
			Usese=0.d0
			IF (split_Asese) THEN
				IF (split_Fsese) THEN
					DO j = 1, N-1, 1
						DO i = j+1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_sese(i,j)=Asese_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsese_yuk))/sij(0,i,j)
							ELSE
								u_sese(i,j)=Asese_ud_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsese_ud_yuk))/sij(0,i,j)
							END IF
							Usese=Usese+u_sese(i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N-1, 1
						DO i = j+1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_sese(i,j)=Asese_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsese_yuk))/sij(0,i,j)
							ELSE
								u_sese(i,j)=Asese_ud_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsese_yuk))/sij(0,i,j)
							END IF
							Usese=Usese+u_sese(i,j)
						END DO
					END DO
				END IF
			ELSE
				IF (split_Fsese) THEN
					DO j = 1, N-1, 1
						DO i = j+1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_sese(i,j)=Asese_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsese_yuk))/sij(0,i,j)
							ELSE
								u_sese(i,j)=Asese_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsese_ud_yuk))/sij(0,i,j)
							END IF
							Usese=Usese+u_sese(i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N-1, 1
						DO i = j+1, N, 1
							u_sese(i,j)=Asese_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsese_yuk))/sij(0,i,j)
							Usese=Usese+u_sese(i,j)
						END DO
					END DO
				END IF
			END IF
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (split_Asese) THEN
				IF (split_Fsese) THEN
					DO i = num+1, N, 1
						Usese=Usese-u_sese(i,num)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_sese(i,num)=Asese_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsese_yuk))/sij(0,i,num)
						ELSE
							u_sese(i,num)=Asese_ud_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsese_ud_yuk))/sij(0,i,num)
						END IF
						Usese=Usese+u_sese(i,num)
					END DO
					DO i = 1, num-1, 1
						Usese=Usese-u_sese(num,i)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_sese(num,i)=Asese_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsese_yuk))/sij(0,num,i)
						ELSE
							u_sese(num,i)=Asese_ud_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsese_ud_yuk))/sij(0,num,i)
						END IF
						Usese=Usese+u_sese(num,i)
					END DO
				ELSE
					DO i = num+1, N, 1
						Usese=Usese-u_sese(i,num)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_sese(i,num)=Asese_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsese_yuk))/sij(0,i,num)
						ELSE
							u_sese(i,num)=Asese_ud_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsese_yuk))/sij(0,i,num)
						END IF
						Usese=Usese+u_sese(i,num)
					END DO
					DO i = 1, num-1, 1
						Usese=Usese-u_sese(num,i)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_sese(num,i)=Asese_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsese_yuk))/sij(0,num,i)
						ELSE
							u_sese(num,i)=Asese_ud_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsese_yuk))/sij(0,num,i)
						END IF
						Usese=Usese+u_sese(num,i)
					END DO
				END IF
			ELSE
				IF (split_Fsese) THEN
					DO i = num+1, N, 1
						Usese=Usese-u_sese(i,num)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_sese(i,num)=Asese_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsese_yuk))/sij(0,i,num)
						ELSE
							u_sese(i,num)=Asese_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsese_ud_yuk))/sij(0,i,num)
						END IF
						Usese=Usese+u_sese(i,num)
					END DO
					DO i = 1, num-1, 1
						Usese=Usese-u_sese(num,i)
						IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
							u_sese(num,i)=Asese_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsese_yuk))/sij(0,num,i)
						ELSE
							u_sese(num,i)=Asese_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsese_ud_yuk))/sij(0,num,i)
						END IF
						Usese=Usese+u_sese(num,i)
					END DO
				ELSE
					DO i = num+1, N, 1
						Usese=Usese-u_sese(i,num)
						u_sese(i,num)=Asese_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsese_yuk))/sij(0,i,num)
						Usese=Usese+u_sese(i,num)
					END DO
					DO i = 1, num-1, 1
						Usese=Usese-u_sese(num,i)
						u_sese(num,i)=Asese_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsese_yuk))/sij(0,num,i)
						Usese=Usese+u_sese(num,i)
					END DO
				END IF
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Usese_YUK ]'
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: Usese(yuk)=', Usese
	END SUBROUTINE valuta_Usese_YUK
!-----------------------------------------------------------------------
	SUBROUTINE seleziona_coppie_migliori(rij,N,partner)
		USE generic_tools
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		INTEGER, SAVE :: N_save=0.d0, coppie_old(1:2)
		INTEGER :: M, i, j, min
		REAL (KIND=8), ALLOCATABLE :: distanza_coppie(:)
		INTEGER, INTENT(OUT) :: partner(1:N)
		
		M=fattoriale_doppio(N-1)
		ALLOCATE(distanza_coppie(1:M))
		distanza_coppie=0.d0
		min=1
		DO i = 1, M, 1
			j=1
			DO WHILE ((j<=N).AND.(distanza_coppie(i)<=distanza_coppie(min)))
				distanza_coppie(i)=distanza_coppie(i)+rij(0,coppie(1,j,i),coppie(2,j,i))
				j=j+1
			END DO
			IF (i>1) THEN
				IF (distanza_coppie(i)<distanza_coppie(min)) min=i
			END IF
		END DO
		DO i = 1, N, 1
			DO j = 1, N/2, 1
				IF (coppie(1,j,min)==i) partner(i)=coppie(2,j,min)
				IF (coppie(2,j,min)==i) partner(i)=coppie(1,j,min)
			END DO
		END DO
		DEALLOCATE(distanza_coppie)
	END SUBROUTINE seleziona_coppie_migliori
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Use_BOUND(num,rij,N,partner,u_ss,Uss)
		USE dati_fisici
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, partner(1:N)
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		INTEGER :: i, j
		REAL (KIND=8) :: u_ss(1:N,1:N), Uss
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Use_BOUND ]'
		
		IF (num==-1) THEN
			Uss=0.d0
			DO i = 1, N, 1
				u_ss(i,partner(i))=B_se*(rij(0,i,partner(i))-D_se)*(rij(0,i,partner(i))-D_se)
				u_ss(partner(i),i)=u_ss(i,partner(i))
				Uss=Uss+u_ss(i,partner(i))
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			Uss=Uss-u_ss(num,partner(num))
			u_ss(num,partner(num))=B_se*(rij(0,num,partner(num))-D_se)*(rij(0,num,partner(num))-D_se)
			u_ss(partner(num),num)=u_ss(num,partner(num))
			Uss=Uss+u_ss(num,partner(num))
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Use_POT ]'
		END IF
		IF (verbose_mode) PRINT * , 'funzione_onda: Use(pot)=', Uss
	END SUBROUTINE valuta_Use_BOUND
!-----------------------------------------------------------------------

	SUBROUTINE valuta_KERNse(num,rij,N,k_se,Kse)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		INTEGER :: i
		REAL (KIND=8) :: k_se(1:N), Kse
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_KERNse ]'
		
		IF (num==-1) THEN
			Kse=0.d0
			DO i = 1, N, 1
				k_se(i)=C_kern_e*rij(0,i,i)*rij(0,i,i)
				Kse=Kse+k_se(i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			Kse=Kse-k_se(num)
			k_se(num)=C_kern_e*rij(0,num,num)*rij(0,num,num)
			Kse=Kse+k_se(num)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_KERNse ]'
		END IF
		IF (verbose_mode) PRINT * , 'funzione_onda: Kse(pot)=', Kse
	END SUBROUTINE valuta_KERNse
!-----------------------------------------------------------------------

	SUBROUTINE valuta_KERNse_ctf(num,rij,N,k_se,Kse)
		USE dati_fisici
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		INTEGER :: i
		REAL (KIND=8) :: k_se(1:N), Kse, minL

		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_KERNse ]'
		
		minL=MINVAL(L)

		IF (num==-1) THEN
			Kse=0.d0
			DO i = 1, N, 1
				IF ( rij(0,i,i)<l1_kern_e ) THEN
					k_se(i)=C_kern_e*rij(0,i,i)*rij(0,i,i)
				ELSE IF ( rij(0,i,i)<minL*0.499 ) THEN     !epsilon=0.001*L
					k_se(i)=alpha0_kern_e+alpha1_kern_e/(rij(0,i,i)-0.5d0*minL)
				ELSE
					k_se(i)=alpha0_kern_e+alpha1_kern_e/(-0.001d0*minL)
				END IF
				Kse=Kse+k_se(i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			Kse=Kse-k_se(num)
			IF ( rij(0,num,num)<l1_kern_e ) THEN
				k_se(num)=C_kern_e*rij(0,num,num)*rij(0,num,num)
			ELSE IF ( rij(0,num,num)<minL*0.499 ) THEN
				k_se(num)=alpha0_kern_e+alpha1_kern_e/(rij(0,num,num)-0.5d0*minL)
			ELSE
				k_se(num)=alpha0_kern_e+alpha1_kern_e/(-0.001d0*minL)
			END IF
			Kse=Kse+k_se(num)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_KERNse ]'
		END IF
				
		IF (verbose_mode) PRINT * , 'funzione_onda: Kse(pot)=', Kse
	END SUBROUTINE valuta_KERNse_ctf	
!-----------------------------------------------------------------------

	SUBROUTINE valuta_KERNse_periodic(num,rij,N,k_se,Kse)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		INTEGER :: i
		REAL (KIND=8) :: k_se(1:N), Kse
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_KERNse ]'
		
		IF (num==-1) THEN
			Kse=0.d0
			DO i = 1, N, 1
				k_se(i)=C_kern_e*rij(0,i,i)*rij(0,i,i)
				Kse=Kse+k_se(i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			Kse=Kse-k_se(num)
			k_se(num)=C_kern_e*rij(0,num,num)*rij(0,num,num)
			Kse=Kse+k_se(num)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_KERNse ]'
		END IF
		IF (verbose_mode) PRINT * , 'funzione_onda: Kse(pot)=', Kse
	END SUBROUTINE valuta_KERNse_periodic	
!-----------------------------------------------------------------------

	SUBROUTINE valuta_atmKERNse(num,rij,N,k_se,Kse)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		INTEGER :: i
		REAL (KIND=8) :: k_se(1:N), Kse
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_atmKERNse ]'
		
		IF (num==-1) THEN
			Kse=0.d0
			DO i = 1, N, 1
				k_se(i)=C_kern_e*rij(0,i,i)
				Kse=Kse+k_se(i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			Kse=Kse-k_se(num)
			k_se(num)=C_kern_e*rij(0,num,num)
			Kse=Kse+k_se(num)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_KERNse ]'
		END IF
		IF (verbose_mode) PRINT * , 'funzione_onda: Kse(pot)=', Kse
	END SUBROUTINE valuta_atmKERNse
!-----------------------------------------------------------------------
	SUBROUTINE valuta_atmKERNse_ctf(num,rij,N,k_se,Kse)
		USE dati_fisici
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		INTEGER :: i
		REAL (KIND=8) :: L38, k_se(1:N), Kse

		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_atmKERNse ]'
		
		L38=3.d0*MINVAL(L)/8.d0
		IF (num==-1) THEN
			Kse=0.d0
			DO i = 1, N, 1
				IF ( rij(0,i,i)<=L38 ) THEN
					k_se(i)=C_kern_e*rij(0,i,i)
				ELSE
					k_se(i)=alpha0_kern_e+alpha1_kern_e*((rij(0,i,i)-L38)**10)
				END IF
				Kse=Kse+k_se(i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			Kse=Kse-k_se(num)
			IF ( rij(0,num,num)<=L38 ) THEN
				k_se(num)=C_kern_e*rij(0,num,num)
			ELSE
				k_se(num)=alpha0_kern_e+alpha1_kern_e*((rij(0,num,num)-L38)**10)
			END IF
			Kse=Kse+k_se(num)
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_KERNse ]'
		END IF
		IF (verbose_mode) PRINT * , 'funzione_onda: Kse(pot)=', Kse
	END SUBROUTINE valuta_atmKERNse_ctf
!-----------------------------------------------------------------------
	SUBROUTINE valuta_GDse(num,eos,rij,N,GD,detGD,IGD,pvt,IGD_old,detGD_old)
		USE generic_tools
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, eos      !N=N_part/2, eop=1 se real, 2 se shadow
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		REAL (KIND=8), INTENT(IN) :: IGD_old(1:N,1:N), detGD_old
		INTEGER :: i, j, index, info, perm
		REAL (KIND=8) :: GD(1:N,1:N), IGD(1:N,1:N), detGD
		INTEGER, INTENT(OUT) :: pvt(1:N)
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_GDse ]'
		
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, N, 1
				DO i = 1, N, 1
					GD(i,j)=DEXP(-C_kern_e*rij(0,i,j)*rij(0,i,j))
					IGD(i,j)=GD(i,j)
				END DO
			END DO
			!Calcolo il determinante di SD_new
			CALL DGETRF( N, N, IGD, N, pvt, info )
			IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
			perm=0
			detGD=1.d0
			DO  i = 1, N, 1
				IF (pvt(i) /= i) perm=perm+1
			END DO
			IF (MOD(perm,2) == 1 ) detGD=-detGD
			DO  i = 1, N, 1
				detGD=detGD*IGD(i,i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (eos==1) THEN
				DO j = 1, N, 1
					GD(num,j)=DEXP(-C_kern_e*rij(0,num,j)*rij(0,num,j))
				END DO
				CALL aggiorna_determinante_R_1ppt(N,num,IGD_old,detGD_old,GD,detGD)
			ELSE IF (eos==2) THEN
				DO j = 1, N, 1
					GD(j,num)=DEXP(-C_kern_e*rij(0,j,num)*rij(0,j,num))
				END DO
				CALL aggiorna_determinante_R_col_1ppt(N,num,IGD_old,detGD_old,GD,detGD)
			ELSE
				STOP 'eop non accettabile &
				  [ module_funzione_onda.f90 > valuta_GDse ]'
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_GDse ]'
		END IF
		
	END SUBROUTINE valuta_GDse
!-----------------------------------------------------------------------
	SUBROUTINE valuta_GDse_ctf(num,eos,rij,N,GD,detGD,IGD,pvt,IGD_old,detGD_old)
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, eos      !N=N_part/2, eop=1 se real, 2 se shadow
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		REAL (KIND=8), INTENT(IN) :: IGD_old(1:N,1:N), detGD_old
		INTEGER :: i, j, index, info, perm
		REAL (KIND=8) :: GD(1:N,1:N), IGD(1:N,1:N), detGD, minL, exp_term
		INTEGER, INTENT(OUT) :: pvt(1:N)

		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_GDse ]'
		
		minL=MINVAL(L)
		
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, N, 1
				DO i = 1, N, 1
					IF ( rij(0,i,j)<l1_kern_e ) THEN
						exp_term=C_kern_e*rij(0,i,j)*rij(0,i,j)
					ELSE IF ( rij(0,i,j)<l2_kern_e ) THEN
						exp_term=alpha0_kern_e+alpha1_kern_e/(rij(0,i,j)-0.5d0*minL)
					ELSE
						exp_term=beta0_kern_e+beta1_kern_e*(rij(0,i,j)**n_kern_e)
					END IF
					GD(i,j)=DEXP(-exp_term)
					IGD(i,j)=GD(i,j)
				END DO
			END DO
			!Calcolo il determinante di SD_new
			CALL DGETRF( N, N, IGD, N, pvt, info )
			IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
			perm=0
			detGD=1.d0
			DO  i = 1, N, 1
				IF (pvt(i) /= i) perm=perm+1
			END DO
			IF (MOD(perm,2) == 1 ) detGD=-detGD
			DO  i = 1, N, 1
				detGD=detGD*IGD(i,i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (eos==1) THEN
				DO j = 1, N, 1
					IF ( rij(0,num,j)<l1_kern_e ) THEN
						exp_term=C_kern_e*rij(0,num,j)*rij(0,num,j)
					ELSE IF ( rij(0,num,j)<l2_kern_e ) THEN
						exp_term=alpha0_kern_e+alpha1_kern_e/(rij(0,num,j)-0.5d0*minL)
					ELSE
						exp_term=beta0_kern_e+beta1_kern_e*(rij(0,num,j)**n_kern_e)
					END IF
					GD(num,j)=DEXP(-exp_term)
				END DO
				CALL aggiorna_determinante_R_1ppt(N,num,IGD_old,detGD_old,GD,detGD)
			ELSE IF (eos==2) THEN
				DO j = 1, N, 1
					IF ( rij(0,j,num)<l1_kern_e ) THEN
						exp_term=C_kern_e*rij(0,j,num)*rij(0,j,num)
					ELSE IF ( rij(0,j,num)<l2_kern_e ) THEN
						exp_term=alpha0_kern_e+alpha1_kern_e/(rij(0,j,num)-0.5d0*minL)
					ELSE
						exp_term=beta0_kern_e+beta1_kern_e*(rij(0,j,num)**n_kern_e)
					END IF
					GD(j,num)=DEXP(-exp_term)
				END DO
				CALL aggiorna_determinante_R_col_1ppt(N,num,IGD_old,detGD_old,GD,detGD)
			ELSE
				STOP 'eop non accettabile &
				  [ module_funzione_onda.f90 > valuta_GDse ]'
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_GDse ]'
		END IF

	END SUBROUTINE valuta_GDse_ctf
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Usesp_POT(num,eop,rij,N,u_ss,Uss)
		USE dati_fisici
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: dist_lim=0.1d0
		INTEGER, INTENT(IN) :: N, num, eop
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		INTEGER :: i, j
		REAL (KIND=8) :: u_ss(1:N,1:N), Uss
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Usesp_POT ]'
		
		IF (num==-1) THEN
			Uss=0.d0
			DO j = 1, N, 1
				DO i = 1, N, 1
					IF (rij(0,i,j)>dist_lim) THEN
						u_ss(i,j)=c_sesp/rij(0,i,j)
					ELSE
						u_ss(i,j)=-rij(0,i,j)/(dist_lim*dist_lim)+2.d0*c_sesp/dist_lim
					END IF
					u_ss(j,i)=u_ss(i,j)
					Uss=Uss+u_ss(i,j)
				END DO
			END DO
		ELSE IF ((num>0) .AND. (num<=N) .AND. (eop==1)) THEN
			DO i = 1, N, 1
				Uss=Uss-u_ss(num,i)
				IF (rij(0,num,i)>dist_lim) THEN
					u_ss(num,i)=c_sesp/rij(0,num,i)
				ELSE
					u_ss(num,i)=-rij(0,num,i)/(dist_lim*dist_lim)+2.d0*c_sesp/dist_lim
				END IF
				u_ss(i,num)=u_ss(num,i)
				Uss=Uss+u_ss(num,i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N) .AND. (eop==2)) THEN
			DO i = 1, N, 1
				Uss=Uss-u_ss(i,num)
				IF (rij(0,i,num)>dist_lim) THEN
					u_ss(i,num)=c_sesp/rij(0,i,num)
				ELSE
					u_ss(i,num)=-rij(0,i,num)/(dist_lim*dist_lim)+2.d0*c_sesp/dist_lim
				END IF
				u_ss(num,i)=u_ss(i,num)
				Uss=Uss+u_ss(i,num)
			END DO
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Usesp_POT ]'
		END IF
		IF (verbose_mode) PRINT * , 'funzione_onda: Usesp(pot)=', Uss
	END SUBROUTINE valuta_Usesp_POT
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Usesp_GSS(num,eop,sij,N,u_sesp,Usesp)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, eop    !eop indica se si é mosso la shadow-e (=1) o la shadow-p (=2), nel caso num>0
		INTEGER :: i, j
		REAL (KIND=8), INTENT(IN) :: sij(0:3,1:N,1:N)
		REAL (KIND=8) :: u_sesp(1:N,1:N), Usesp
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Usesp_GSS ]'
		
		IF (num==-1) THEN
			Usesp=0.d0
			!*!DO j = 1, N, 1            *sommando su tutte le possibili combinazioni le cose non funzionano
			!*!	DO i = 1, N, 1
			!*!		u_sesp(i,j)=Gsesp*sij(0,i,j)*sij(0,i,j)
			!*!		Usesp=Usesp+u_sesp(i,j)
			!*!	END DO
			!*!END DO
			DO i = 1, N, 1
				u_sesp(i,i)=Gsesp*sij(0,i,i)*sij(0,i,i)
				Usesp=Usesp+u_sesp(i,i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (eop==1) THEN
				Usesp=Usesp-u_sesp(num,num)
				u_sesp(num,num)=Gsesp*sij(0,num,num)*sij(0,num,num)
				Usesp=Usesp+u_sesp(num,num)
			ELSE IF (eop==2) THEN
				Usesp=Usesp-u_sesp(num,num)
				u_sesp(num,num)=Gsesp*sij(0,num,num)*sij(0,num,num)
				Usesp=Usesp+u_sesp(num,num)
			ELSE
				STOP 'eop non accettabile &
				  [ module_funzione_onda.f90 > valuta_Usesp_GSS ]'
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Usesp_GSS ]'
		END IF

		!IF (verbose_mode) PRINT * , 'funzione_onda: Usesp(gss)=', Usesp
		!PRINT * , 'funzione_onda: Usesp(gss)=', Usesp
	END SUBROUTINE valuta_Usesp_GSS
!-----------------------------------------------------------------------
	SUBROUTINE valuta_GDsp(num,eos,rij,N,GD,detGD,IGD,pvt,IGD_old,detGD_old)
		USE generic_tools
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: fatt_Gsesp=1.d0
		INTEGER, INTENT(IN) :: N, num, eos      !N=N_part/2, eop=1 se real, 2 se shadow
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N)
		REAL (KIND=8), INTENT(IN) :: IGD_old(1:N,1:N), detGD_old
		INTEGER :: i, j, index, info, perm
		REAL (KIND=8) :: GD(1:N,1:N), IGD(1:N,1:N), detGD
		INTEGER, INTENT(OUT) :: pvt(1:N)
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_GDsp ]'
		
		IF (num==-1) THEN
			!Calcolo i termini matriciali di SD_new
			DO j = 1, N, 1
				DO i = 1, N, 1
					GD(i,j)=DEXP(-fatt_Gsesp*Gsesp*rij(0,i,j)*rij(0,i,j))
					IGD(i,j)=GD(i,j)
				END DO
			END DO
			CALL DGETRF( N, N, IGD, N, pvt, info )
			IF (info/=0) STOP 'ERRORE NELLA DECOMPOSIZIONE LU'
			perm=0
			detGD=1.d0
			DO  i = 1, N, 1
				IF (pvt(i) /= i) perm=perm+1
			END DO
			IF (MOD(perm,2) == 1 ) detGD=-detGD
			DO  i = 1, N, 1
				detGD=detGD*IGD(i,i)
			END DO
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (eos==1) THEN
				DO j = 1, N, 1
					GD(num,j)=DEXP(-fatt_Gsesp*Gsesp*rij(0,num,j)*rij(0,num,j))
				END DO
				CALL aggiorna_determinante_R_1ppt(N,num,IGD_old,detGD_old,GD,detGD)
			ELSE IF (eos==2) THEN
				DO j = 1, N, 1
					GD(j,num)=DEXP(-Gsesp*rij(0,j,num)*rij(0,j,num))
				END DO
				CALL aggiorna_determinante_R_col_1ppt(N,num,IGD_old,detGD_old,GD,detGD)
			ELSE
				STOP 'eop non accettabile &
				  [ module_funzione_onda.f90 > valuta_GDsp ]'
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_GDsp ]'
		END IF
		
	END SUBROUTINE valuta_GDsp
!-----------------------------------------------------------------------
	SUBROUTINE valuta_Usesp_YUK(num,eop,sij,N,u_sesp,Usesp)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, eop    !eop indica se si é mosso la shadow-e (=1) o la shadow-p (=2), nel caso num>0
		INTEGER :: i, j, H_N
		REAL (KIND=8), INTENT(IN) :: sij(0:3,1:N,1:N)
		REAL (KIND=8) :: u_sesp(1:N,1:N), Usesp
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > valuta_Usesp_YUK ]'
		
		H_N=N/2
		IF (num==-1) THEN
			Usesp=0.d0
			IF (split_Asesp) THEN
				IF (split_Fsesp) THEN
					DO j = 1, N, 1
						DO i = 1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_sesp(i,j)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsesp_yuk))/sij(0,i,j)
							ELSE
								u_sesp(i,j)=Asesp_ud_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsesp_ud_yuk))/sij(0,i,j)
							END IF
							Usesp=Usesp+u_sesp(i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N, 1
						DO i = 1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_sesp(i,j)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsesp_yuk))/sij(0,i,j)
							ELSE
								u_sesp(i,j)=Asesp_ud_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsesp_yuk))/sij(0,i,j)
							END IF
							Usesp=Usesp+u_sesp(i,j)
						END DO
					END DO
				END IF
			ELSE
				IF (split_Fsesp) THEN
					DO j = 1, N, 1
						DO i = 1, N, 1
							IF ( (i<=H_N .AND. j<=H_N) .OR. (i>H_N .AND. j>H_N) ) THEN
								u_sesp(i,j)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsesp_yuk))/sij(0,i,j)
							ELSE
								u_sesp(i,j)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsesp_ud_yuk))/sij(0,i,j)
							END IF
							Usesp=Usesp+u_sesp(i,j)
						END DO
					END DO
				ELSE
					DO j = 1, N, 1
						DO i = 1, N, 1
							u_sesp(i,j)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,j)*Fsesp_yuk))/sij(0,i,j)
							Usesp=Usesp+u_sesp(i,j)
						END DO
					END DO
				END IF
			END IF
		ELSE IF ((num>0) .AND. (num<=N)) THEN
			IF (split_Asesp) THEN
				IF (split_Fsesp) THEN
					IF (eop==1) THEN
						DO i = 1, N, 1
							Usesp=Usesp-u_sesp(num,i)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_sesp(num,i)=Asesp_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsesp_yuk))/sij(0,num,i)
							ELSE
								u_sesp(num,i)=Asesp_ud_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsesp_ud_yuk))/sij(0,num,i)
							END IF
							Usesp=Usesp+u_sesp(num,i)
						END DO
					ELSE IF (eop==2) THEN
						DO i = 1, N, 1
							Usesp=Usesp-u_sesp(i,num)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_sesp(i,num)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsesp_yuk))/sij(0,i,num)
							ELSE
								u_sesp(i,num)=Asesp_ud_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsesp_ud_yuk))/sij(0,i,num)
							END IF
							Usesp=Usesp+u_sesp(i,num)
						END DO
					ELSE
						STOP 'eop non accettabile &
						  [ module_funzione_onda.f90 > valuta_Usesp_YUK ]'
					END IF
				ELSE
					IF (eop==1) THEN
						DO i = 1, N, 1
							Usesp=Usesp-u_sesp(num,i)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_sesp(num,i)=Asesp_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsesp_yuk))/sij(0,num,i)
							ELSE
								u_sesp(num,i)=Asesp_ud_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsesp_yuk))/sij(0,num,i)
							END IF
							Usesp=Usesp+u_sesp(num,i)
						END DO
					ELSE IF (eop==2) THEN
						DO i = 1, N, 1
							Usesp=Usesp-u_sesp(i,num)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_sesp(i,num)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsesp_yuk))/sij(0,i,num)
							ELSE
								u_sesp(i,num)=Asesp_ud_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsesp_yuk))/sij(0,i,num)
							END IF
							Usesp=Usesp+u_sesp(i,num)
						END DO
					ELSE
						STOP 'eop non accettabile &
						  [ module_funzione_onda.f90 > valuta_Usesp_YUK ]'
					END IF
				END IF
			ELSE
				IF (split_Fsesp) THEN
					IF (eop==1) THEN
						DO i = 1, N, 1
							Usesp=Usesp-u_sesp(num,i)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_sesp(num,i)=Asesp_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsesp_yuk))/sij(0,num,i)
							ELSE
								u_sesp(num,i)=Asesp_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsesp_ud_yuk))/sij(0,num,i)
							END IF
							Usesp=Usesp+u_sesp(num,i)
						END DO
					ELSE IF (eop==2) THEN
						DO i = 1, N, 1
							Usesp=Usesp-u_sesp(i,num)
							IF ( (i<=H_N .AND. num<=H_N) .OR. (i>H_N .AND. num>H_N) ) THEN
								u_sesp(i,num)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsesp_yuk))/sij(0,i,num)
							ELSE
								u_sesp(i,num)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsesp_ud_yuk))/sij(0,i,num)
							END IF
							Usesp=Usesp+u_sesp(i,num)
						END DO
					ELSE
						STOP 'eop non accettabile &
						  [ module_funzione_onda.f90 > valuta_Usesp_YUK ]'
					END IF
				ELSE
					IF (eop==1) THEN
						DO i = 1, N, 1
							Usesp=Usesp-u_sesp(num,i)
							u_sesp(num,i)=Asesp_yuk*(1.d0-DEXP(-sij(0,num,i)*Fsesp_yuk))/sij(0,num,i)
							Usesp=Usesp+u_sesp(num,i)
						END DO
					ELSE IF (eop==2) THEN
						DO i = 1, N, 1
							Usesp=Usesp-u_sesp(i,num)
							u_sesp(i,num)=Asesp_yuk*(1.d0-DEXP(-sij(0,i,num)*Fsesp_yuk))/sij(0,i,num)
							Usesp=Usesp+u_sesp(i,num)
						END DO
					ELSE
						STOP 'eop non accettabile &
						  [ module_funzione_onda.f90 > valuta_Usesp_YUK ]'
					END IF
				END IF
			END IF
		ELSE
			STOP 'num non accettabile &
			  [ module_funzione_onda.f90 > valuta_Usesp_YUK ]'
		END IF
		
		IF (verbose_mode) PRINT * , 'funzione_onda: Usesp(yuk)=', Usesp
	END SUBROUTINE valuta_Usesp_YUK
!-----------------------------------------------------------------------

	SUBROUTINE chiudi_dati_funzione_onda()
		USE dnfH
		USE dati_mc
		IMPLICIT NONE
		
		IF (.NOT. iniz_funzione_onda) STOP 'funzione_onda non é inizializzato &
		  [ module_funzione_onda.f90 > chiudi_dati_funzione_onda ]'
		
		IF ((SDe_kind=='lda').OR.(SDse_kind=='lda')) THEN
			IF ( flag_TABC ) THEN
				DEALLOCATE(pesi_K_points)
			END IF
			DEALLOCATE(k_pw_lda,fattori_orb_lda,fattori_pw_lda,twist_lda)
			DEALLOCATE(num_fpw_lda,k_fpw_lda,fattori_fpw_lda)
			num_chiamata_twist_lda=0
		ELSE IF ((SDe_kind=='prf').OR.(SDe_kind=='fre')) THEN
			DEALLOCATE(k_pw_lda,fattori_pw_lda)
			DEALLOCATE(num_pw_orbit,indice_pw_dnfH)
			DEALLOCATE(k_pw_dnfH,fattori_pw_dnfH,num_pw_dnfH)
			IF (flag_usa_coeff_dnfH) THEN                        !se voglio usare i coefficienti nella matrice dnfH, li carico
				DEALLOCATE(c_eff_dnfH)
			END IF
			DEALLOCATE(autoenergie_dnfH)
      ELSE IF (SDe_kind=='spb') THEN
         CALL MSPL_deallocate(SPL=Bsplep)
		END IF
      IF ((Jee_kind=='spl').OR.(Jee_kind=='spp')) THEN
         CALL MSPL_deallocate(SPL=Jsplee)
         IF (split_Aee.OR.split_Fee) CALL MSPL_deallocate(SPL=Jsplee_ud) 
      END IF
      IF ((Jep_kind=='spl').OR.(Jep_kind=='spp')) THEN
         CALL MSPL_deallocate(SPL=Jsplep)
         IF (split_Aep.OR.split_Fep) CALL MSPL_deallocate(SPL=Jsplep_ud) 
      END IF
		IF ((Jsesp_kind=='bou') .OR. (Jsesp_kind=='ppb')) THEN
			DEALLOCATE(coppie)
		END IF
		
		iniz_funzione_onda=.FALSE.
		
	END SUBROUTINE chiudi_dati_funzione_onda
	
END MODULE funzione_onda
