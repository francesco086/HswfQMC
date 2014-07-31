MODULE calcola_accettazione
	IMPLICIT NONE
	LOGICAL, PARAMETER, PRIVATE :: verbose_mode=.FALSE.
	LOGICAL, SAVE, PRIVATE :: iniz_calcola_accettazione=.FALSE.
	COMPLEX (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: SDe_up_new(:,:), SDe_up_old(:,:), ISDe_up_new(:,:), ISDe_up_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvte_up_new(:), pvte_up_old(:)
	COMPLEX (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: SDe_dw_new(:,:), SDe_dw_old(:,:), ISDe_dw_new(:,:), ISDe_dw_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvte_dw_new(:), pvte_dw_old(:)
	COMPLEX (KIND=8), SAVE, PROTECTED ::  detSDe_up_new, detSDe_up_old, detSDe_dw_new, detSDe_dw_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: u_ee_new(:,:), u_ee_old(:,:), u_ep_new(:,:), u_ep_old(:,:)
	REAL (KIND=8), SAVE, PROTECTED :: Uee_new, Uee_old, Uep_new, Uep_old
	COMPLEX (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: SDse1_up_new(:,:), SDse1_up_old(:,:), ISDse1_up_new(:,:), ISDse1_up_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtse1_up_new(:), pvtse1_up_old(:)
	COMPLEX (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: SDse1_dw_new(:,:), SDse1_dw_old(:,:), ISDse1_dw_new(:,:), ISDse1_dw_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtse1_dw_new(:), pvtse1_dw_old(:)
	COMPLEX (KIND=8), SAVE, PROTECTED ::  detSDse1_up_new, detSDse1_up_old, detSDse1_dw_new, detSDse1_dw_old
	COMPLEX (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: SDse2_up_new(:,:), SDse2_up_old(:,:), ISDse2_up_new(:,:), ISDse2_up_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtse2_up_new(:), pvtse2_up_old(:)
	COMPLEX (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: SDse2_dw_new(:,:), SDse2_dw_old(:,:), ISDse2_dw_new(:,:), ISDse2_dw_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtse2_dw_new(:), pvtse2_dw_old(:)
	COMPLEX (KIND=8), SAVE, PROTECTED ::  detSDse2_up_new, detSDse2_up_old, detSDse2_dw_new, detSDse2_dw_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: u_se1_new(:,:), u_se1_old(:,:), u_se2_new(:,:), u_se2_old(:,:)
	REAL (KIND=8), SAVE, PROTECTED :: Use1_new, Use1_old, Use2_new, Use2_old
	INTEGER, ALLOCATABLE :: partner_se1(:), partner_se2(:)
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: b_se1_new(:,:), b_se1_old(:,:), b_se2_new(:,:), b_se2_old(:,:)
	REAL (KIND=8), SAVE, PROTECTED :: Bse1_new, Bse1_old, Bse2_new, Bse2_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: u_ese1_new(:), u_ese1_old(:), u_ese2_new(:), u_ese2_old(:)
	REAL (KIND=8), SAVE, PROTECTED :: Uese1_new, Uese1_old, Uese2_new, Uese2_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: GDse1_up_new(:,:), GDse1_up_old(:,:), IGDse1_up_new(:,:), IGDse1_up_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtgdse1_up_new(:), pvtgdse1_up_old(:)
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: GDse1_dw_new(:,:), GDse1_dw_old(:,:), IGDse1_dw_new(:,:), IGDse1_dw_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtgdse1_dw_new(:), pvtgdse1_dw_old(:)
	REAL (KIND=8), SAVE, PROTECTED ::  detGDse1_up_new, detGDse1_up_old, detGDse1_dw_new, detGDse1_dw_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: GDse2_up_new(:,:), GDse2_up_old(:,:), IGDse2_up_new(:,:), IGDse2_up_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtgdse2_up_new(:), pvtgdse2_up_old(:)
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: GDse2_dw_new(:,:), GDse2_dw_old(:,:), IGDse2_dw_new(:,:), IGDse2_dw_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtgdse2_dw_new(:), pvtgdse2_dw_old(:)
	REAL (KIND=8), SAVE, PROTECTED ::  detGDse2_up_new, detGDse2_up_old, detGDse2_dw_new, detGDse2_dw_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: u_sesp1_new(:,:), u_sesp1_old(:,:), u_sesp2_new(:,:), u_sesp2_old(:,:)
	REAL (KIND=8), SAVE, PROTECTED :: Usesp1_new, Usesp1_old, Usesp2_new, Usesp2_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: GDsp1_up_new(:,:), GDsp1_up_old(:,:), IGDsp1_up_new(:,:), IGDsp1_up_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtgdsp1_up_new(:), pvtgdsp1_up_old(:)
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: GDsp1_dw_new(:,:), GDsp1_dw_old(:,:), IGDsp1_dw_new(:,:), IGDsp1_dw_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtgdsp1_dw_new(:), pvtgdsp1_dw_old(:)
	REAL (KIND=8), SAVE, PROTECTED ::  detGDsp1_up_new, detGDsp1_up_old, detGDsp1_dw_new, detGDsp1_dw_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: GDsp2_up_new(:,:), GDsp2_up_old(:,:), IGDsp2_up_new(:,:), IGDsp2_up_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtgdsp2_up_new(:), pvtgdsp2_up_old(:)
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: GDsp2_dw_new(:,:), GDsp2_dw_old(:,:), IGDsp2_dw_new(:,:), IGDsp2_dw_old(:,:)
	INTEGER, ALLOCATABLE, SAVE, PROTECTED :: pvtgdsp2_dw_new(:), pvtgdsp2_dw_old(:)
	REAL (KIND=8), SAVE, PROTECTED ::  detGDsp2_up_new, detGDsp2_up_old, detGDsp2_dw_new, detGDsp2_dw_old
	REAL (KIND=8), ALLOCATABLE, SAVE, PROTECTED :: u_POT_se1_new(:), u_POT_se1_old(:), u_POT_se2_new(:), u_POT_se2_old(:)
	
	CONTAINS
	
	SUBROUTINE inizializza_funzione_onda()
		USE walkers
		USE dati_fisici
		USE dati_simulazione_mc
		USE momenta
		USE funzione_onda
		IMPLICIT NONE
		INTEGER :: info, i
		COMPLEX (KIND=8) :: work(1:3*H_N_part)
				
		IF (.NOT. iniz_dati_fisici) STOP 'Prima devi leggere i dati fisici &
		  [ module_calcola_accettazione.f90 > inizializza_funzione_onda ]'
		IF (.NOT. iniz_walkers) STOP 'Prima devi inizializzare i walkers &
		  [ module_calcola_accettazione.f90 > inizializza_funzione_onda ]'
		
		CALL inizializza_momenta(H_N_part, num_k_ewald, L, mpi_myrank)
		CALL inizializza_dati_funzione_onda()
				
		IF (flag_TABC) THEN
			CALL applica_twist(H_N_part, L)
			IF ((SDe_kind=='prf').OR.(SDe_kind=='fre')) CALL applica_twist_hartree()
		END IF
						
		IF (SDe_kind/='no_') THEN
			ALLOCATE(SDe_up_new(1:H_N_part,1:H_N_part), SDe_up_old(1:H_N_part,1:H_N_part))
			SDe_up_new=0.d0
			SDe_up_old=0.d0
			ALLOCATE(ISDe_up_new(1:H_N_part,1:H_N_part), ISDe_up_old(1:H_N_part,1:H_N_part))
			ISDe_up_new=0.d0
			ISDe_up_old=0.d0
			ALLOCATE(pvte_up_new(1:H_N_part), pvte_up_old(1:H_N_part))
			pvte_up_new=0.d0
			pvte_up_old=0.d0
			ALLOCATE(SDe_dw_new(1:H_N_part,1:H_N_part), SDe_dw_old(1:H_N_part,1:H_N_part))
			SDe_dw_new=0.d0
			SDe_dw_old=0.d0
			ALLOCATE(ISDe_dw_new(1:H_N_part,1:H_N_part), ISDe_dw_old(1:H_N_part,1:H_N_part))
			ISDe_dw_new=0.d0
			ISDe_dw_old=0.d0
			ALLOCATE(pvte_dw_new(1:H_N_part), pvte_dw_old(1:H_N_part))
			pvte_dw_new=0.d0
			pvte_dw_old=0.d0
		END IF
		IF (Jee_kind/='no_') THEN
			ALLOCATE(u_ee_new(1:N_part,1:N_part), u_ee_old(1:N_part,1:N_part))
			u_ee_new=0.d0
			u_ee_old=0.d0
		END IF
		IF (Jep_kind/='no_') THEN
			ALLOCATE(u_ep_new(1:N_part,1:N_part), u_ep_old(1:N_part,1:N_part))
			u_ep_new=0.d0
			u_ep_old=0.d0
		END IF
		IF (SDse_kind/='no_') THEN
			ALLOCATE(SDse1_up_new(1:H_N_part,1:H_N_part), SDse1_up_old(1:H_N_part,1:H_N_part))
			SDse1_up_new=0.d0
			SDse1_up_old=0.d0
			ALLOCATE(ISDse1_up_new(1:H_N_part,1:H_N_part), ISDse1_up_old(1:H_N_part,1:H_N_part))
			ISDse1_up_new=0.d0
			ISDse1_up_old=0.d0
			ALLOCATE(pvtse1_up_new(1:H_N_part), pvtse1_up_old(1:H_N_part))
			pvtse1_up_new=0.d0
			pvtse1_up_old=0.d0
			ALLOCATE(SDse1_dw_new(1:H_N_part,1:H_N_part), SDse1_dw_old(1:H_N_part,1:H_N_part))
			SDse1_dw_new=0.d0
			SDse1_dw_old=0.d0
			ALLOCATE(ISDse1_dw_new(1:H_N_part,1:H_N_part), ISDse1_dw_old(1:H_N_part,1:H_N_part))
			ISDse1_dw_new=0.d0
			ISDse1_dw_old=0.d0
			ALLOCATE(pvtse1_dw_new(1:H_N_part), pvtse1_dw_old(1:H_N_part))
			pvtse1_dw_new=0.d0
			pvtse1_dw_old=0.d0
			ALLOCATE(SDse2_up_new(1:H_N_part,1:H_N_part), SDse2_up_old(1:H_N_part,1:H_N_part))
			SDse2_up_new=0.d0
			SDse2_up_old=0.d0
			ALLOCATE(ISDse2_up_new(1:H_N_part,1:H_N_part), ISDse2_up_old(1:H_N_part,1:H_N_part))
			ISDse2_up_new=0.d0
			ISDse2_up_old=0.d0
			ALLOCATE(pvtse2_up_new(1:H_N_part), pvtse2_up_old(1:H_N_part))
			pvtse2_up_new=0.d0
			pvtse2_up_old=0.d0
			ALLOCATE(SDse2_dw_new(1:H_N_part,1:H_N_part), SDse2_dw_old(1:H_N_part,1:H_N_part))
			SDse2_dw_new=0.d0
			SDse2_dw_old=0.d0
			ALLOCATE(ISDse2_dw_new(1:H_N_part,1:H_N_part), ISDse2_dw_old(1:H_N_part,1:H_N_part))
			ISDse2_dw_new=0.d0
			ISDse2_dw_old=0.d0
			ALLOCATE(pvtse2_dw_new(1:H_N_part), pvtse2_dw_old(1:H_N_part))
			pvtse2_dw_new=0.d0
			pvtse2_dw_old=0.d0
		END IF
		IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
			ALLOCATE(u_ese1_new(1:N_part), u_ese1_old(1:N_part))
			u_ese1_new=0.d0
			u_ese1_old=0.d0
			ALLOCATE(u_ese2_new(1:N_part), u_ese2_old(1:N_part))
			u_ese2_new=0.d0
			u_ese2_old=0.d0
		ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
			ALLOCATE(GDse1_up_new(1:H_N_part,1:H_N_part), GDse1_up_old(1:H_N_part,1:H_N_part))
			GDse1_up_new=0.d0
			GDse1_up_old=0.d0
			ALLOCATE(IGDse1_up_new(1:H_N_part,1:H_N_part), IGDse1_up_old(1:H_N_part,1:H_N_part))
			IGDse1_up_new=0.d0
			IGDse1_up_old=0.d0
			ALLOCATE(pvtgdse1_up_new(1:H_N_part), pvtgdse1_up_old(1:H_N_part))
			pvtgdse1_up_new=0.d0
			pvtgdse1_up_old=0.d0
			ALLOCATE(GDse1_dw_new(1:H_N_part,1:H_N_part), GDse1_dw_old(1:H_N_part,1:H_N_part))
			GDse1_dw_new=0.d0
			GDse1_dw_old=0.d0
			ALLOCATE(IGDse1_dw_new(1:H_N_part,1:H_N_part), IGDse1_dw_old(1:H_N_part,1:H_N_part))
			IGDse1_dw_new=0.d0
			IGDse1_dw_old=0.d0
			ALLOCATE(pvtgdse1_dw_new(1:H_N_part), pvtgdse1_dw_old(1:H_N_part))
			pvtgdse1_dw_new=0.d0
			pvtgdse1_dw_old=0.d0
			ALLOCATE(GDse2_up_new(1:H_N_part,1:H_N_part), GDse2_up_old(1:H_N_part,1:H_N_part))
			GDse2_up_new=0.d0
			GDse2_up_old=0.d0
			ALLOCATE(IGDse2_up_new(1:H_N_part,1:H_N_part), IGDse2_up_old(1:H_N_part,1:H_N_part))
			IGDse2_up_new=0.d0
			IGDse2_up_old=0.d0
			ALLOCATE(pvtgdse2_up_new(1:H_N_part), pvtgdse2_up_old(1:H_N_part))
			pvtgdse2_up_new=0.d0
			pvtgdse2_up_old=0.d0
			ALLOCATE(GDse2_dw_new(1:H_N_part,1:H_N_part), GDse2_dw_old(1:H_N_part,1:H_N_part))
			GDse2_dw_new=0.d0
			GDse2_dw_old=0.d0
			ALLOCATE(IGDse2_dw_new(1:H_N_part,1:H_N_part), IGDse2_dw_old(1:H_N_part,1:H_N_part))
			IGDse2_dw_new=0.d0
			IGDse2_dw_old=0.d0
			ALLOCATE(pvtgdse2_dw_new(1:H_N_part), pvtgdse2_dw_old(1:H_N_part))
			pvtgdse2_dw_new=0.d0
			pvtgdse2_dw_old=0.d0
		ELSE IF ((Kse_kind=='atm').OR.(Kse_kind=='atc')) THEN
			ALLOCATE(u_ese1_new(1:N_part), u_ese1_old(1:N_part))
			u_ese1_new=0.d0
			u_ese1_old=0.d0
			ALLOCATE(u_ese2_new(1:N_part), u_ese2_old(1:N_part))
			u_ese2_new=0.d0
			u_ese2_old=0.d0
		END IF
		IF (Jse_kind/='no_') THEN
			IF ((Jse_kind/='no_') .AND. (Jse_kind/='bou') .AND. (Jse_kind/='pot')) THEN
				ALLOCATE(u_se1_new(1:N_part,1:N_part), u_se1_old(1:N_part,1:N_part))
				u_se1_new=0.d0
				u_se1_old=0.d0
				ALLOCATE(u_se2_new(1:N_part,1:N_part), u_se2_old(1:N_part,1:N_part))
				u_se2_new=0.d0
				u_se2_old=0.d0
			END IF
			IF ((Jse_kind=='bou') .OR. (Jse_kind=='ppb')) THEN
				STOP 'Jse_kin=bou e ppt ancora da implementare &
				  [ module_calcola_accettazione.f90 > inizializza_funzione_onda ]'
				!ALLOCATE(partner_se1(1:N_part),partner_se2(1:N_part))
				!ALLOCATE(b_se1_new(1:N_part,1:N_part), b_se1_old(1:N_part,1:N_part))
				!b_se1_new=0.d0
				!b_se1_old=0.d0
				!ALLOCATE(b_se2_new(1:N_part,1:N_part), b_se2_old(1:N_part,1:N_part))
				!b_se2_new=0.d0
				!b_se2_old=0.d0
			END IF
			IF ( Jse_kind=='pot' ) THEN
				ALLOCATE(u_POT_se1_new(1:H_N_part), u_POT_se1_old(1:H_N_part))
				u_POT_se1_new=0.d0
				u_POT_se1_old=0.d0
				ALLOCATE(u_POT_se2_new(1:H_N_part), u_POT_se2_old(1:H_N_part))
				u_POT_se2_new=0.d0
				u_POT_se2_old=0.d0
			END IF
			Use1_old=0.d0
			Use1_new=0.d0
			Use2_old=0.d0
			Use2_new=0.d0
		END IF
		IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
			ALLOCATE(u_sesp1_new(1:N_part,1:N_part), u_sesp1_old(1:N_part,1:N_part))
			u_sesp1_new=0.d0
			u_sesp1_old=0.d0
			ALLOCATE(u_sesp2_new(1:N_part,1:N_part), u_sesp2_old(1:N_part,1:N_part))
			u_sesp2_new=0.d0
			u_sesp2_old=0.d0
		ELSE IF (Jsesp_kind=='gsd') THEN
			ALLOCATE(GDsp1_up_new(1:H_N_part,1:H_N_part), GDsp1_up_old(1:H_N_part,1:H_N_part))
			GDsp1_up_new=0.d0
			GDsp1_up_old=0.d0
			ALLOCATE(IGDsp1_up_new(1:H_N_part,1:H_N_part), IGDsp1_up_old(1:H_N_part,1:H_N_part))
			IGDsp1_up_new=0.d0
			IGDsp1_up_old=0.d0
			ALLOCATE(pvtgdsp1_up_new(1:H_N_part), pvtgdsp1_up_old(1:H_N_part))
			pvtgdsp1_up_new=0.d0
			pvtgdsp1_up_old=0.d0
			ALLOCATE(GDsp1_dw_new(1:H_N_part,1:H_N_part), GDsp1_dw_old(1:H_N_part,1:H_N_part))
			GDsp1_dw_new=0.d0
			GDsp1_dw_old=0.d0
			ALLOCATE(IGDsp1_dw_new(1:H_N_part,1:H_N_part), IGDsp1_dw_old(1:H_N_part,1:H_N_part))
			IGDsp1_dw_new=0.d0
			IGDsp1_dw_old=0.d0
			ALLOCATE(pvtgdsp1_dw_new(1:H_N_part), pvtgdsp1_dw_old(1:H_N_part))
			pvtgdsp1_dw_new=0.d0
			pvtgdsp1_dw_old=0.d0
			ALLOCATE(GDsp2_up_new(1:H_N_part,1:H_N_part), GDsp2_up_old(1:H_N_part,1:H_N_part))
			GDsp2_up_new=0.d0
			GDsp2_up_old=0.d0
			ALLOCATE(IGDsp2_up_new(1:H_N_part,1:H_N_part), IGDsp2_up_old(1:H_N_part,1:H_N_part))
			IGDsp2_up_new=0.d0
			IGDsp2_up_old=0.d0
			ALLOCATE(pvtgdsp2_up_new(1:H_N_part), pvtgdsp2_up_old(1:H_N_part))
			pvtgdsp2_up_new=0.d0
			pvtgdsp2_up_old=0.d0
			ALLOCATE(GDsp2_dw_new(1:H_N_part,1:H_N_part), GDsp2_dw_old(1:H_N_part,1:H_N_part))
			GDsp2_dw_new=0.d0
			GDsp2_dw_old=0.d0
			ALLOCATE(IGDsp2_dw_new(1:H_N_part,1:H_N_part), IGDsp2_dw_old(1:H_N_part,1:H_N_part))
			IGDsp2_dw_new=0.d0
			IGDsp2_dw_old=0.d0
			ALLOCATE(pvtgdsp2_dw_new(1:H_N_part), pvtgdsp2_dw_old(1:H_N_part))
			pvtgdsp2_dw_new=0.d0
			pvtgdsp2_dw_old=0.d0
		END IF
		
		iniz_calcola_accettazione=.TRUE.
		CALL prima_valutazione_funzione_onda()
		!!!CALL aggiorna_funzione_onda('all',-1)
				
		IF (verbose_mode) PRINT *, 'calcola_accettazione: Ho inizializzato la funzione d onda'
		
		!IF (verbose_mode) THEN
		!	IF (SDe_kind/='no_') THEN
		!		PRINT * , 'detSDe_up_new -> ', REAL(detSDe_up_old*DCONJG(detSDe_up_old))
		!		PRINT * , 'detSDe_dw_new -> ', REAL(detSDe_dw_old*DCONJG(detSDe_dw_old))
		!	END IF
		!	IF (Jee_kind/='no_') PRINT * , 'Uee -> ', Uee_old
		!	IF (Jep_kind/='no_') PRINT * , 'Uep -> ', Uep_old
		!	IF (SDse_kind/='no_') THEN
		!		PRINT * , 'detSDse1_up_new -> ', REAL(detSDse1_up_old*DCONJG(detSDse1_up_old))
		!		PRINT * , 'detSDse1_dw_new -> ', REAL(detSDse1_dw_old*DCONJG(detSDse1_dw_old))
		!		PRINT * , 'detSDse2_up_new -> ', REAL(detSDse2_up_old*DCONJG(detSDse2_up_old))
		!		PRINT * , 'detSDse2_dw_new -> ', REAL(detSDse2_dw_old*DCONJG(detSDse2_dw_old))
		!	END IF
		!	IF (Jse_kind/='no_') PRINT * , 'Use1 -> ', Use1_old
		!	IF (Jse_kind/='no_') PRINT * , 'Use2 -> ', Use2_old
		!	IF ((Kse_kind/='no_').AND.((Kse_kind/='gsd').OR.(Kse_kind/='gdc'))) PRINT * , 'Kse1 -> ', Uese1_old
		!	IF ((Kse_kind/='no_').AND.((Kse_kind/='gsd').OR.(Kse_kind/='gdc'))) PRINT * , 'Kse2 -> ', Uese2_old
		!	IF ((Kse_kind/='gsd').OR.(Kse_kind/='gdc')) PRINT * , 'Kse1 -> ', detGDse1_up_old*detGDse1_dw_old
		!	IF ((Kse_kind/='gsd').OR.(Kse_kind/='gdc')) PRINT * , 'Kse2 -> ', detGDse2_up_old*detGDse2_dw_old
		!	IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) PRINT * , 'Usesp1 -> ', Usesp1_old
		!	IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) PRINT * , 'Usesp2 -> ', Usesp2_old
		!	IF (Jsesp_kind=='gsd') PRINT * , 'Ksp1 -> ', detGDsp1_up_old*detGDsp1_dw_old
		!	IF (Jsesp_kind=='gsd') PRINT * , 'Ksp2 -> ', detGDsp2_up_old*detGDsp2_dw_old
		!END IF
												
	END SUBROUTINE inizializza_funzione_onda
!-----------------------------------------------------------------------

	SUBROUTINE prima_valutazione_funzione_onda()
		USE walkers
		USE dati_fisici
		USE dati_simulazione_mc
		USE momenta
		USE funzione_onda
		IMPLICIT NONE
		INTEGER :: info, i, j
		COMPLEX (KIND=8) :: work(1:3*H_N_part)
				
		IF (.NOT. iniz_calcola_accettazione) STOP 'Prima devi inizializzare la funzione d onda &
		  [ module_calcola_accettazione.f90 > reinizializza_funzion_onda ]'
		IF (.NOT. iniz_dati_fisici) STOP 'Prima devi leggere i dati fisici &
		  [ module_calcola_accettazione.f90 > reinizializza_funzion_onda ]'
		IF (.NOT. iniz_momenta) STOP 'Prima devi trovare i momenti di fermi &
		  [ module_calcola_accettazione.f90 > reinizializza_funzion_onda ]'
		IF (.NOT. iniz_walkers) STOP 'Prima devi inizializzare i walkers &
		  [ module_calcola_accettazione.f90 > reinizializza_funzion_onda ]'
		
		SELECT CASE (SDe_kind)
		CASE ('pw_') 
			CALL valuta_SD_pw(-1,re_old(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDe_up_old,detSDe_up_old,ISDe_up_old,pvte_up_old,SDe_up_new,detSDe_up_new)
			CALL valuta_SD_pw(-1,re_old(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDe_dw_old,detSDe_dw_old,ISDe_dw_old,pvte_dw_old,SDe_dw_new,detSDe_dw_new)
			CALL ZGETRI( H_N_part, ISDe_up_old, H_N_part, pvte_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA R UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDe_dw_old, H_N_part, pvte_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA R DOWN. INFO=', info
				STOP
			END IF
		CASE ('lda') 
			CALL valuta_SD_lda(-1,re_old(1:3,1:H_N_part),H_N_part, &
			  SDe_up_old,detSDe_up_old,ISDe_up_old,pvte_up_old,SDe_up_new,detSDe_up_new)
			CALL valuta_SD_lda(-1,re_old(1:3,H_N_part+1:N_part),H_N_part, &
			  SDe_dw_old,detSDe_dw_old,ISDe_dw_old,pvte_dw_old,SDe_dw_new,detSDe_dw_new)
			CALL ZGETRI( H_N_part, ISDe_up_old, H_N_part, pvte_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA R UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDe_dw_old, H_N_part, pvte_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA R DOWN. INFO=', info
				STOP
			END IF
		CASE ('prf')
			CALL valuta_SD_har(-1,re_old(1:3,1:H_N_part),H_N_part, &
			  SDe_up_old,detSDe_up_old,ISDe_up_old,pvte_up_old,SDe_up_new,detSDe_up_new)
			CALL valuta_SD_har(-1,re_old(1:3,H_N_part+1:N_part),H_N_part, &
			  SDe_dw_old,detSDe_dw_old,ISDe_dw_old,pvte_dw_old,SDe_dw_new,detSDe_dw_new)
			CALL ZGETRI( H_N_part, ISDe_up_old, H_N_part, pvte_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA R UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDe_dw_old, H_N_part, pvte_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA R DOWN. INFO=', info
				STOP
			END IF
		CASE ('fre')
			CALL valuta_SD_har(-1,re_old(1:3,1:H_N_part),H_N_part, &
			  SDe_up_old,detSDe_up_old,ISDe_up_old,pvte_up_old,SDe_up_new,detSDe_up_new)
			CALL valuta_SD_har(-1,re_old(1:3,H_N_part+1:N_part),H_N_part, &
			  SDe_dw_old,detSDe_dw_old,ISDe_dw_old,pvte_dw_old,SDe_dw_new,detSDe_dw_new)
			CALL ZGETRI( H_N_part, ISDe_up_old, H_N_part, pvte_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA R UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDe_dw_old, H_N_part, pvte_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA R DOWN. INFO=', info
				STOP
			END IF
		CASE ('atm')
			CALL valuta_SD_atm(-1,rij_ep_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDe_up_old,detSDe_up_old,ISDe_up_old,pvte_up_old,SDe_up_new,detSDe_up_new)
			CALL valuta_SD_atm(-1,rij_ep_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDe_dw_old,detSDe_dw_old,ISDe_dw_old,pvte_dw_old,SDe_dw_new,detSDe_dw_new)
			CALL ZGETRI( H_N_part, ISDe_up_old, H_N_part, pvte_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA E UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDe_dw_old, H_N_part, pvte_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA E DOWN. INFO=', info
				STOP
			END IF
		CASE ('atp')
			CALL attiva_pc()
			CALL valuta_SD_atm(-1,rijpc_ep_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDe_up_old,detSDe_up_old,ISDe_up_old,pvte_up_old,SDe_up_new,detSDe_up_new)
			CALL valuta_SD_atm(-1,rijpc_ep_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDe_dw_old,detSDe_dw_old,ISDe_dw_old,pvte_dw_old,SDe_dw_new,detSDe_dw_new)
			CALL ZGETRI( H_N_part, ISDe_up_old, H_N_part, pvte_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA E UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDe_dw_old, H_N_part, pvte_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA E DOWN. INFO=', info
				STOP
			END IF
		CASE ('no_') 
			detSDe_up_old=1.d0
			detSDe_dw_old=1.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di SDe_kind accettabile &
			  [ module_calcola_accettazione.f90 > prima_valutazione_funzione_onda ]'
		END SELECT
		IF (SDe_kind/='no_') THEN
			SDe_up_new=SDe_up_old
			ISDe_up_new=ISDe_up_old
			SDe_dw_new=SDe_dw_old
			ISDe_dw_new=ISDe_dw_old
		END IF
		detSDe_up_new=detSDe_up_old
		detSDe_dw_new=detSDe_dw_old
		
		SELECT CASE (Jee_kind)
		CASE ('yuk')
			CALL valuta_Uee_YUK(-1,rij_ee_old,N_part,u_ee_old,Uee_old)
		CASE ('yup') 
			CALL attiva_pc()
			CALL valuta_Uee_YUK(-1,rijpc_ee_old,N_part,u_ee_old,Uee_old)
		CASE ('no_') 
			Uee_old=0.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Jee_kind accettabile &
			  [ module_calcola_accettazione.f90 > prima_valutazione_funzione_onda ]'
		END SELECT
		IF (Jee_kind/='no_') u_ee_new=u_ee_old
		Uee_new=Uee_old
		
		SELECT CASE (Jep_kind)
		CASE ('yuk')
			CALL valuta_Uep_YUK(-1,0,rij_ep_old,N_part,u_ep_old,Uep_old)
		CASE ('yup') 
			CALL attiva_pc()
			CALL valuta_Uep_YUK(-1,0,rijpc_ep_old,N_part,u_ep_old,Uep_old)
		CASE ('atm')
			CALL valuta_Uep_ATM(-1,0,rij_ep_old,N_part,u_ep_old,Uep_old)
		CASE ('atp') 
			CALL attiva_pc()
			CALL valuta_Uep_ATM(-1,0,rijpc_ep_old,N_part,u_ep_old,Uep_old)
		CASE ('no_') 
			Uep_old=0.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Jep_kind accettabile &
			  [ module_calcola_accettazione.f90 > prima_valutazione_funzione_onda ]'
		END SELECT
		IF (Jep_kind/='no_') u_ep_new=u_ep_old
		Uep_new=Uep_old
		
		SELECT CASE (SDse_kind)
		CASE ('pw_') 
			CALL valuta_SD_pw(-1,se1_old(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDse1_up_old,detSDse1_up_old,ISDse1_up_old,pvtse1_up_old,SDse1_up_new,detSDse1_up_new)
			CALL valuta_SD_pw(-1,se1_old(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDse1_dw_old,detSDse1_dw_old,ISDse1_dw_old,pvtse1_dw_old,SDse1_dw_new,detSDse1_dw_new)
			CALL ZGETRI( H_N_part, ISDse1_up_old, H_N_part, pvtse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse1_dw_old, H_N_part, pvtse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_SD_pw(-1,se2_old(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDse2_up_old,detSDse2_up_old,ISDse2_up_old,pvtse2_up_old,SDse2_up_new,detSDse2_up_new)
			CALL valuta_SD_pw(-1,se2_old(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDse2_dw_old,detSDse2_dw_old,ISDse2_dw_old,pvtse2_dw_old,SDse2_dw_new,detSDse2_dw_new)
			CALL ZGETRI( H_N_part, ISDse2_up_old, H_N_part, pvtse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse2_dw_old, H_N_part, pvtse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 DOWN. INFO=', info
				STOP
			END IF
		CASE ('pw2') 
			CALL valuta_SD_pw(-1,se1_old(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDse1_up_old,detSDse1_up_old,ISDse1_up_old,pvtse1_up_old,SDse1_up_new,detSDse1_up_new)
			CALL valuta_SD_pw(-1,se1_old(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDse1_dw_old,detSDse1_dw_old,ISDse1_dw_old,pvtse1_dw_old,SDse1_dw_new,detSDse1_dw_new)
			CALL ZGETRI( H_N_part, ISDse1_up_old, H_N_part, pvtse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse1_dw_old, H_N_part, pvtse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_SD_pw(-1,se2_old(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDse2_up_old,detSDse2_up_old,ISDse2_up_old,pvtse2_up_old,SDse2_up_new,detSDse2_up_new)
			CALL valuta_SD_pw(-1,se2_old(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
			  SDse2_dw_old,detSDse2_dw_old,ISDse2_dw_old,pvtse2_dw_old,SDse2_dw_new,detSDse2_dw_new)
			CALL ZGETRI( H_N_part, ISDse2_up_old, H_N_part, pvtse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse2_dw_old, H_N_part, pvtse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 DOWN. INFO=', info
				STOP
			END IF
		CASE ('lda')
			CALL valuta_SD_lda(-1,se1_old(1:3,1:H_N_part),H_N_part, &
			  SDse1_up_old,detSDse1_up_old,ISDse1_up_old,pvtse1_up_old,SDse1_up_new,detSDse1_up_new)
			CALL valuta_SD_lda(-1,se1_old(1:3,H_N_part+1:N_part),H_N_part, &
			  SDse1_dw_old,detSDse1_dw_old,ISDse1_dw_old,pvtse1_dw_old,SDse1_dw_new,detSDse1_dw_new)
			CALL ZGETRI( H_N_part, ISDse1_up_old, H_N_part, pvtse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse1_dw_old, H_N_part, pvtse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_SD_lda(-1,se2_old(1:3,1:H_N_part),H_N_part, &
			  SDse2_up_old,detSDse2_up_old,ISDse2_up_old,pvtse2_up_old,SDse2_up_new,detSDse2_up_new)
			CALL valuta_SD_lda(-1,se2_old(1:3,H_N_part+1:N_part),H_N_part, &
			  SDse2_dw_old,detSDse2_dw_old,ISDse2_dw_old,pvtse2_dw_old,SDse2_dw_new,detSDse2_dw_new)
			CALL ZGETRI( H_N_part, ISDse2_up_old, H_N_part, pvtse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse2_dw_old, H_N_part, pvtse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 DOWN. INFO=', info
				STOP
			END IF
		CASE ('gem')
			CALL valuta_SD_gem(-1,rij_se1_old(0,1:N_part,1:N_part),H_N_part, &
			  SDse1_up_old,detSDse1_up_old,ISDse1_up_old,pvtse1_up_old,SDse1_up_new,detSDse1_up_new)
			CALL ZGETRI( H_N_part, ISDse1_up_old, H_N_part, pvtse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 UP. INFO=', info
				STOP
			END IF
			CALL valuta_SD_gem(-1,rij_se2_old(0,1:N_part,1:N_part),H_N_part, &
			  SDse2_up_old,detSDse2_up_old,ISDse2_up_old,pvtse2_up_old,SDse2_up_new,detSDse2_up_new)
			CALL ZGETRI( H_N_part, ISDse2_up_old, H_N_part, pvtse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 UP. INFO=', info
				STOP
			END IF
			detSDse1_dw_old=1.d0
			detSDse2_dw_old=1.d0
		CASE ('gss')
			CALL valuta_SD_gauss(-1,rij_sesp1_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDse1_up_old,detSDse1_up_old,ISDse1_up_old,pvtse1_up_old,SDse1_up_new,detSDse1_up_new)
			CALL valuta_SD_gauss(-1,rij_sesp1_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDse1_dw_old,detSDse1_dw_old,ISDse1_dw_old,pvtse1_dw_old,SDse1_dw_new,detSDse1_dw_new)
			CALL ZGETRI( H_N_part, ISDse1_up_old, H_N_part, pvtse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse1_dw_old, H_N_part, pvtse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_SD_gauss(-1,rij_sesp2_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDse2_up_old,detSDse2_up_old,ISDse2_up_old,pvtse2_up_old,SDse2_up_new,detSDse2_up_new)
			CALL valuta_SD_gauss(-1,rij_sesp2_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDse2_dw_old,detSDse2_dw_old,ISDse2_dw_old,pvtse2_dw_old,SDse2_dw_new,detSDse2_dw_new)
			CALL ZGETRI( H_N_part, ISDse2_up_old, H_N_part, pvtse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse2_dw_old, H_N_part, pvtse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 DOWN. INFO=', info
				STOP
			END IF
		CASE ('gsp')
			CALL attiva_pc()
			CALL valuta_SD_gauss(-1,rijpc_sesp1_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDse1_up_old,detSDse1_up_old,ISDse1_up_old,pvtse1_up_old,SDse1_up_new,detSDse1_up_new)
			CALL valuta_SD_gauss(-1,rijpc_sesp1_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDse1_dw_old,detSDse1_dw_old,ISDse1_dw_old,pvtse1_dw_old,SDse1_dw_new,detSDse1_dw_new)
			CALL ZGETRI( H_N_part, ISDse1_up_old, H_N_part, pvtse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse1_dw_old, H_N_part, pvtse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_SD_gauss(-1,rijpc_sesp2_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDse2_up_old,detSDse2_up_old,ISDse2_up_old,pvtse2_up_old,SDse2_up_new,detSDse2_up_new)
			CALL valuta_SD_gauss(-1,rijpc_sesp2_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDse2_dw_old,detSDse2_dw_old,ISDse2_dw_old,pvtse2_dw_old,SDse2_dw_new,detSDse2_dw_new)
			CALL ZGETRI( H_N_part, ISDse2_up_old, H_N_part, pvtse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse2_dw_old, H_N_part, pvtse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 DOWN. INFO=', info
				STOP
			END IF
		CASE ('atm')     !dovrei usare rij_sep ma per il momento coincide con rij_sesp, quindi meglio fare cosí
			CALL valuta_SD_atm(-1,rij_sesp1_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDse1_up_old,detSDse1_up_old,ISDse1_up_old,pvtse1_up_old,SDse1_up_new,detSDse1_up_new)
			CALL valuta_SD_atm(-1,rij_sesp1_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDse1_dw_old,detSDse1_dw_old,ISDse1_dw_old,pvtse1_dw_old,SDse1_dw_new,detSDse1_dw_new)
			CALL ZGETRI( H_N_part, ISDse1_up_old, H_N_part, pvtse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse1_dw_old, H_N_part, pvtse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_SD_atm(-1,rij_sesp2_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDse2_up_old,detSDse2_up_old,ISDse2_up_old,pvtse2_up_old,SDse2_up_new,detSDse2_up_new)
			CALL valuta_SD_atm(-1,rij_sesp2_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDse2_dw_old,detSDse2_dw_old,ISDse2_dw_old,pvtse2_dw_old,SDse2_dw_new,detSDse2_dw_new)
			CALL ZGETRI( H_N_part, ISDse2_up_old, H_N_part, pvtse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse2_dw_old, H_N_part, pvtse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 DOWN. INFO=', info
				STOP
			END IF
		CASE ('atp')     !dovrei usare rij_sep ma per il momento coincide con rij_sesp, quindi meglio fare cosí
			CALL attiva_pc()
			CALL valuta_SD_atm(-1,rijpc_sesp1_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDse1_up_old,detSDse1_up_old,ISDse1_up_old,pvtse1_up_old,SDse1_up_new,detSDse1_up_new)
			CALL valuta_SD_atm(-1,rijpc_sesp1_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDse1_dw_old,detSDse1_dw_old,ISDse1_dw_old,pvtse1_dw_old,SDse1_dw_new,detSDse1_dw_new)
			CALL ZGETRI( H_N_part, ISDse1_up_old, H_N_part, pvtse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse1_dw_old, H_N_part, pvtse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_SD_atm(-1,rijpc_sesp2_old(0,1:H_N_part,1:H_N_part),H_N_part, &
			  SDse2_up_old,detSDse2_up_old,ISDse2_up_old,pvtse2_up_old,SDse2_up_new,detSDse2_up_new)
			CALL valuta_SD_atm(-1,rijpc_sesp2_old(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  SDse2_dw_old,detSDse2_dw_old,ISDse2_dw_old,pvtse2_dw_old,SDse2_dw_new,detSDse2_dw_new)
			CALL ZGETRI( H_N_part, ISDse2_up_old, H_N_part, pvtse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 UP. INFO=', info
				STOP
			END IF
			CALL ZGETRI( H_N_part, ISDse2_dw_old, H_N_part, pvtse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA SE2 DOWN. INFO=', info
				STOP
			END IF
		CASE ('no_') 
			detSDse1_up_old=1.d0
			detSDse1_dw_old=1.d0
			detSDse2_up_old=1.d0
			detSDse2_dw_old=1.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di SDse_kind accettabile &
			  [ module_calcola_accettazione.f90 > prima_valutazione_funzione_onda ]'
		END SELECT
		IF (SDse_kind/='no_') THEN
			SDse1_up_new=SDse1_up_old
			ISDse1_up_new=ISDse1_up_old
			SDse1_dw_new=SDse1_dw_old
			ISDse1_dw_new=ISDse1_dw_old
			SDse2_up_new=SDse2_up_old
			ISDse2_up_new=ISDse2_up_old
			SDse2_dw_new=SDse2_dw_old
			ISDse2_dw_new=ISDse2_dw_old
		END IF
		detSDse1_up_new=detSDse1_up_old
		detSDse1_dw_new=detSDse1_dw_old
		detSDse2_up_new=detSDse2_up_old
		detSDse2_dw_new=detSDse2_dw_old
		
		Use1_old=0.d0
		Use2_old=0.d0
		Bse1_old=0.d0
		Bse2_old=0.d0
		SELECT CASE (Jse_kind)
		CASE ('no_')
		CASE ('pot')
			CALL valuta_Use_POT(-1,dist_mol_ss1_old,u_POT_se1_old,Use1_old)
			CALL valuta_Use_POT(-1,dist_mol_ss2_old,u_POT_se2_old,Use2_old)
			u_POT_se1_new=u_POT_se1_old
			u_POT_se2_new=u_POT_se2_old
		CASE ('bou')
			!CALL seleziona_coppie_migliori(rij_se1_old,N_part,partner_se1)
			!CALL valuta_Use_BOUND(-1,rij_se1_old,N_part,partner_se1,b_se1_old,Bse1_old)
			!CALL seleziona_coppie_migliori(rij_se2_old,N_part,partner_se2)
			!CALL valuta_Use_BOUND(-1,rij_se2_old,N_part,partner_se2,b_se2_old,Bse2_old)
			!b_se1_new=b_se1_old
			!b_se2_new=b_se2_old
		CASE ('ppb')
			!CALL valuta_Use_POT(-1,rij_se1_old,N_part,u_se1_old,Use1_old)
			!CALL valuta_Use_POT(-1,rij_se2_old,N_part,u_se2_old,Use2_old)
			!u_se1_new=u_se1_old
			!u_se2_new=u_se2_old
			!CALL seleziona_coppie_migliori(rij_se1_old,N_part,partner_se1)
			!CALL valuta_Use_BOUND(-1,rij_se1_old,N_part,partner_se1,b_se1_old,Bse1_old)
			!CALL seleziona_coppie_migliori(rij_se2_old,N_part,partner_se2)
			!CALL valuta_Use_BOUND(-1,rij_se2_old,N_part,partner_se2,b_se2_old,Bse2_old)
			!b_se1_new=b_se1_old
			!b_se2_new=b_se2_old
		CASE ('yuk')
			CALL valuta_Usese_YUK(-1,rij_se1_old,N_part,u_se1_old,Use1_old)
			CALL valuta_Usese_YUK(-1,rij_se2_old,N_part,u_se2_old,Use2_old)
			u_se1_new=u_se1_old
			u_se2_new=u_se2_old
		CASE ('yup') 
			CALL attiva_pc()
			CALL valuta_Usese_YUK(-1,rijpc_se1_old,N_part,u_se1_old,Use1_old)
			CALL valuta_Usese_YUK(-1,rijpc_se2_old,N_part,u_se2_old,Use2_old)
			u_se1_new=u_se1_old
			u_se2_new=u_se2_old
		CASE DEFAULT
			PRINT * , ' ### ', Jse_kind
			STOP 'Non hai selezionato un valore di Jse_kind accettabile &
			  [ module_calcola_accettazione.f90 > prima_valutazione_funzione_onda ]'
		END SELECT
		Use1_new=Use1_old
		Use2_new=Use2_old
		Bse1_new=Bse1_old
		Bse2_new=Bse2_old
				
		SELECT CASE (Kse_kind)
		CASE ('gss')
			CALL valuta_KERNse(-1,rij_ese1_old,N_part,u_ese1_old,Uese1_old)
			CALL valuta_KERNse(-1,rij_ese2_old,N_part,u_ese2_old,Uese2_old)
			detGDse1_up_old=1.d0
			detGDse1_dw_old=1.d0
			detGDse1_up_new=1.d0
			detGDse1_dw_new=1.d0
			detGDse2_up_old=1.d0
			detGDse2_dw_old=1.d0
			detGDse2_up_new=1.d0
			detGDse2_dw_new=1.d0
		CASE ('gsc')
			CALL valuta_KERNse_ctf(-1,rij_ese1_old,N_part,u_ese1_old,Uese1_old)
			CALL valuta_KERNse_ctf(-1,rij_ese2_old,N_part,u_ese2_old,Uese2_old)
			detGDse1_up_old=1.d0
			detGDse1_dw_old=1.d0
			detGDse1_up_new=1.d0
			detGDse1_dw_new=1.d0
			detGDse2_up_old=1.d0
			detGDse2_dw_old=1.d0
			detGDse2_up_new=1.d0
			detGDse2_dw_new=1.d0
		CASE ('gsp')
			CALL attiva_pc()
			CALL valuta_KERNse(-1,rijpc_ese1_old,N_part,u_ese1_old,Uese1_old)
			CALL valuta_KERNse(-1,rijpc_ese2_old,N_part,u_ese2_old,Uese2_old)
			detGDse1_up_old=1.d0
			detGDse1_dw_old=1.d0
			detGDse1_up_new=1.d0
			detGDse1_dw_new=1.d0
			detGDse2_up_old=1.d0
			detGDse2_dw_old=1.d0
			detGDse2_up_new=1.d0
			detGDse2_dw_new=1.d0
		CASE ('gsd')
			CALL valuta_GDse(-1,0,rij_ese1_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  GDse1_up_old,detGDse1_up_old,IGDse1_up_old,pvtgdse1_up_old,GDse1_up_old,detGDse1_up_old)
			CALL valuta_GDse(-1,0,rij_ese1_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  GDse1_dw_old,detGDse1_dw_old,IGDse1_dw_old,pvtgdse1_dw_old,GDse1_dw_old,detGDse1_dw_old)
			CALL DGETRI( H_N_part, IGDse1_up_old, H_N_part, pvtgdse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE1 UP. INFO=', info
				STOP
			END IF
			CALL DGETRI( H_N_part, IGDse1_dw_old, H_N_part, pvtgdse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_GDse(-1,0,rij_ese2_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  GDse2_up_old,detGDse2_up_old,IGDse2_up_old,pvtgdse2_up_old,GDse2_up_old,detGDse2_up_old)
			CALL valuta_GDse(-1,0,rij_ese2_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  GDse2_dw_old,detGDse2_dw_old,IGDse2_dw_old,pvtgdse2_dw_old,GDse2_dw_old,detGDse2_dw_old)
			CALL DGETRI( H_N_part, IGDse2_up_old, H_N_part, pvtgdse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE2 UP. INFO=', info
				STOP
			END IF
			CALL DGETRI( H_N_part, IGDse2_dw_old, H_N_part, pvtgdse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE2 DOWN. INFO=', info
				STOP
			END IF
			GDse1_up_new=GDse1_up_old
			IGDse1_up_new=IGDse1_up_old
			GDse1_dw_new=GDse1_dw_old
			IGDse1_dw_new=IGDse1_dw_old
			pvtgdse1_up_new=pvtgdse1_up_old
			pvtgdse1_dw_new=pvtgdse1_dw_old
			detGDse1_up_new=detGDse1_up_old
			detGDse1_dw_new=detGDse1_dw_old
			GDse2_up_new=GDse2_up_old
			IGDse2_up_new=IGDse2_up_old
			GDse2_dw_new=GDse2_dw_old
			IGDse2_dw_new=IGDse2_dw_old
			pvtgdse2_up_new=pvtgdse2_up_old
			pvtgdse2_dw_new=pvtgdse2_dw_old
			detGDse2_up_new=detGDse2_up_old
			detGDse2_dw_new=detGDse2_dw_old
			Uese1_old=0.d0
			Uese2_old=0.d0
			Uese1_new=0.d0
			Uese2_new=0.d0
		CASE ('gdc')
			CALL valuta_GDse_ctf(-1,0,rij_ese1_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  GDse1_up_old,detGDse1_up_old,IGDse1_up_old,pvtgdse1_up_old,GDse1_up_old,detGDse1_up_old)
			CALL valuta_GDse_ctf(-1,0,rij_ese1_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  GDse1_dw_old,detGDse1_dw_old,IGDse1_dw_old,pvtgdse1_dw_old,GDse1_dw_old,detGDse1_dw_old)
			CALL DGETRI( H_N_part, IGDse1_up_old, H_N_part, pvtgdse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE1 UP. INFO=', info
				STOP
			END IF
			CALL DGETRI( H_N_part, IGDse1_dw_old, H_N_part, pvtgdse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_GDse_ctf(-1,0,rij_ese2_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  GDse2_up_old,detGDse2_up_old,IGDse2_up_old,pvtgdse2_up_old,GDse2_up_old,detGDse2_up_old)
			CALL valuta_GDse_ctf(-1,0,rij_ese2_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  GDse2_dw_old,detGDse2_dw_old,IGDse2_dw_old,pvtgdse2_dw_old,GDse2_dw_old,detGDse2_dw_old)
			CALL DGETRI( H_N_part, IGDse2_up_old, H_N_part, pvtgdse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE2 UP. INFO=', info
				STOP
			END IF
			CALL DGETRI( H_N_part, IGDse2_dw_old, H_N_part, pvtgdse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE2 DOWN. INFO=', info
				STOP
			END IF
			GDse1_up_new=GDse1_up_old
			IGDse1_up_new=IGDse1_up_old
			GDse1_dw_new=GDse1_dw_old
			IGDse1_dw_new=IGDse1_dw_old
			pvtgdse1_up_new=pvtgdse1_up_old
			pvtgdse1_dw_new=pvtgdse1_dw_old
			detGDse1_up_new=detGDse1_up_old
			detGDse1_dw_new=detGDse1_dw_old
			GDse2_up_new=GDse2_up_old
			IGDse2_up_new=IGDse2_up_old
			GDse2_dw_new=GDse2_dw_old
			IGDse2_dw_new=IGDse2_dw_old
			pvtgdse2_up_new=pvtgdse2_up_old
			pvtgdse2_dw_new=pvtgdse2_dw_old
			detGDse2_up_new=detGDse2_up_old
			detGDse2_dw_new=detGDse2_dw_old
			Uese1_old=0.d0
			Uese2_old=0.d0
			Uese1_new=0.d0
			Uese2_new=0.d0
		CASE ('gdp')
			CALL attiva_pc()
			CALL valuta_GDse(-1,0,rijpc_ese1_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  GDse1_up_old,detGDse1_up_old,IGDse1_up_old,pvtgdse1_up_old,GDse1_up_old,detGDse1_up_old)
			CALL valuta_GDse(-1,0,rijpc_ese1_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  GDse1_dw_old,detGDse1_dw_old,IGDse1_dw_old,pvtgdse1_dw_old,GDse1_dw_old,detGDse1_dw_old)
			CALL DGETRI( H_N_part, IGDse1_up_old, H_N_part, pvtgdse1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE1 UP. INFO=', info
				STOP
			END IF
			CALL DGETRI( H_N_part, IGDse1_dw_old, H_N_part, pvtgdse1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_GDse(-1,0,rijpc_ese2_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  GDse2_up_old,detGDse2_up_old,IGDse2_up_old,pvtgdse2_up_old,GDse2_up_old,detGDse2_up_old)
			CALL valuta_GDse(-1,0,rijpc_ese2_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  GDse2_dw_old,detGDse2_dw_old,IGDse2_dw_old,pvtgdse2_dw_old,GDse2_dw_old,detGDse2_dw_old)
			CALL DGETRI( H_N_part, IGDse2_up_old, H_N_part, pvtgdse2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE2 UP. INFO=', info
				STOP
			END IF
			CALL DGETRI( H_N_part, IGDse2_dw_old, H_N_part, pvtgdse2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE2 DOWN. INFO=', info
				STOP
			END IF
			GDse1_up_new=GDse1_up_old
			IGDse1_up_new=IGDse1_up_old
			GDse1_dw_new=GDse1_dw_old
			IGDse1_dw_new=IGDse1_dw_old
			pvtgdse1_up_new=pvtgdse1_up_old
			pvtgdse1_dw_new=pvtgdse1_dw_old
			detGDse1_up_new=detGDse1_up_old
			detGDse1_dw_new=detGDse1_dw_old
			GDse2_up_new=GDse2_up_old
			IGDse2_up_new=IGDse2_up_old
			GDse2_dw_new=GDse2_dw_old
			IGDse2_dw_new=IGDse2_dw_old
			pvtgdse2_up_new=pvtgdse2_up_old
			pvtgdse2_dw_new=pvtgdse2_dw_old
			detGDse2_up_new=detGDse2_up_old
			detGDse2_dw_new=detGDse2_dw_old
			Uese1_old=0.d0
			Uese2_old=0.d0
		CASE ('atm')
			CALL valuta_atmKERNse(-1,rij_ese1_old,N_part,u_ese1_old,Uese1_old)
			CALL valuta_atmKERNse(-1,rij_ese2_old,N_part,u_ese2_old,Uese2_old)
			detGDse1_up_old=1.d0
			detGDse1_dw_old=1.d0
			detGDse2_up_old=1.d0
			detGDse2_dw_old=1.d0
		CASE ('atc')
			CALL valuta_atmKERNse_ctf(-1,rij_ese1_old,N_part,u_ese1_old,Uese1_old)
			CALL valuta_atmKERNse_ctf(-1,rij_ese2_old,N_part,u_ese2_old,Uese2_old)
			detGDse1_up_old=1.d0
			detGDse1_dw_old=1.d0
			detGDse2_up_old=1.d0
			detGDse2_dw_old=1.d0
		CASE ('no_')
			Uese1_old=0.d0
			Uese2_old=0.d0
			detGDse1_up_old=1.d0
			detGDse1_dw_old=1.d0
			detGDse2_up_old=1.d0
			detGDse2_dw_old=1.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Kse_kind accettabile &
			  [ module_calcola_accettazione.f90 > prima_valutazione_funzione_onda ]'
		END SELECT
		Uese1_new=Uese1_old
		Uese2_new=Uese2_old
		detGDse1_up_new=detGDse1_up_old
		detGDse1_dw_new=detGDse1_dw_old
		detGDse2_up_new=detGDse2_up_old
		detGDse2_dw_new=detGDse2_dw_old
		
		SELECT CASE (Jsesp_kind)
		CASE ('pot')
			CALL valuta_Usesp_POT(-1,0,rij_sesp1_old,N_part,u_sesp1_old,Usesp1_old)
			CALL valuta_Usesp_POT(-1,0,rij_sesp2_old,N_part,u_sesp2_old,Usesp2_old)
			detGDsp1_up_old=1.d0
			detGDsp1_dw_old=1.d0
			detGDsp1_up_new=1.d0
			detGDsp1_dw_new=1.d0
			detGDsp2_up_old=1.d0
			detGDsp2_dw_old=1.d0
			detGDsp2_up_new=1.d0
			detGDsp2_dw_new=1.d0
		CASE ('yuk')
			CALL valuta_Usesp_YUK(-1,0,rij_sesp1_old,N_part,u_sesp1_old,Usesp1_old)
			CALL valuta_Usesp_YUK(-1,0,rij_sesp2_old,N_part,u_sesp2_old,Usesp2_old)
			detGDsp1_up_old=1.d0
			detGDsp1_dw_old=1.d0
			detGDsp1_up_new=1.d0
			detGDsp1_dw_new=1.d0
			detGDsp2_up_old=1.d0
			detGDsp2_dw_old=1.d0
			detGDsp2_up_new=1.d0
			detGDsp2_dw_new=1.d0
		CASE ('yup') 
			CALL attiva_pc()
			CALL valuta_Usesp_YUK(-1,0,rijpc_sesp1_old,N_part,u_sesp1_old,Usesp1_old)
			CALL valuta_Usesp_YUK(-1,0,rijpc_sesp2_old,N_part,u_sesp2_old,Usesp2_old)
			detGDsp1_up_old=1.d0
			detGDsp1_dw_old=1.d0
			detGDsp1_up_new=1.d0
			detGDsp1_dw_new=1.d0
			detGDsp2_up_old=1.d0
			detGDsp2_dw_old=1.d0
			detGDsp2_up_new=1.d0
			detGDsp2_dw_new=1.d0
		CASE ('gss')
			CALL valuta_Usesp_GSS(-1,0,rij_sesp1_old,N_part,u_sesp1_old,Usesp1_old)
			CALL valuta_Usesp_GSS(-1,0,rij_sesp2_old,N_part,u_sesp2_old,Usesp2_old)
			detGDsp1_up_old=1.d0
			detGDsp1_dw_old=1.d0
			detGDsp1_up_new=1.d0
			detGDsp1_dw_new=1.d0
			detGDsp2_up_old=1.d0
			detGDsp2_dw_old=1.d0
			detGDsp2_up_new=1.d0
			detGDsp2_dw_new=1.d0
		CASE ('gsd')
			CALL valuta_GDsp(-1,0,rij_sesp1_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  GDsp1_up_old,detGDsp1_up_old,IGDsp1_up_old,pvtgdsp1_up_old,GDsp1_up_old,detGDsp1_up_old)
			CALL valuta_GDsp(-1,0,rij_sesp1_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  GDsp1_dw_old,detGDsp1_dw_old,IGDsp1_dw_old,pvtgdsp1_dw_old,GDsp1_dw_old,detGDsp1_dw_old)
			CALL DGETRI( H_N_part, IGDsp1_up_old, H_N_part, pvtgdsp1_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSP1 UP. INFO=', info
				STOP
			END IF
			CALL DGETRI( H_N_part, IGDsp1_dw_old, H_N_part, pvtgdsp1_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSP1 DOWN. INFO=', info
				STOP
			END IF
			CALL valuta_GDsp(-1,0,rij_sesp2_old(0:3,1:H_N_part,1:H_N_part),H_N_part, &
			  GDsp2_up_old,detGDsp2_up_old,IGDsp2_up_old,pvtgdsp2_up_old,GDsp2_up_old,detGDsp2_up_old)
			CALL valuta_GDsp(-1,0,rij_sesp2_old(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
			  GDsp2_dw_old,detGDsp2_dw_old,IGDsp2_dw_old,pvtgdsp2_dw_old,GDsp2_dw_old,detGDsp2_dw_old)
			CALL DGETRI( H_N_part, IGDsp2_up_old, H_N_part, pvtgdsp2_up_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE2 UP. INFO=', info
				STOP
			END IF
			CALL DGETRI( H_N_part, IGDsp2_dw_old, H_N_part, pvtgdsp2_dw_old, work, 3*H_N_part, info )
			IF (info/=0) THEN
				PRINT *, 'ERRORE NEL TROVARE LA MATRICE INVERSA GDSE2 DOWN. INFO=', info
				STOP
			END IF
			GDsp1_up_new=GDsp1_up_old
			IGDsp1_up_new=IGDsp1_up_old
			GDsp1_dw_new=GDsp1_dw_old
			IGDsp1_dw_new=IGDsp1_dw_old
			pvtgdsp1_up_new=pvtgdsp1_up_old
			pvtgdsp1_dw_new=pvtgdsp1_dw_old
			detGDsp1_up_new=detGDsp1_up_old
			detGDsp1_dw_new=detGDsp1_dw_old
			GDsp2_up_new=GDsp2_up_old
			IGDsp2_up_new=IGDsp2_up_old
			GDsp2_dw_new=GDsp2_dw_old
			IGDsp2_dw_new=IGDsp2_dw_old
			pvtgdsp2_up_new=pvtgdsp2_up_old
			pvtgdsp2_dw_new=pvtgdsp2_dw_old
			detGDsp2_up_new=detGDsp2_up_old
			detGDsp2_dw_new=detGDsp2_dw_old
			Usesp1_old=0.d0
			Usesp2_old=0.d0
			Usesp1_new=0.d0
			Usesp2_new=0.d0
		CASE ('no_') 
			Usesp1_old=0.d0
			Usesp2_old=0.d0
			detGDsp1_up_old=1.d0
			detGDsp1_dw_old=1.d0
			detGDsp1_up_new=1.d0
			detGDsp1_dw_new=1.d0
			detGDsp2_up_old=1.d0
			detGDsp2_dw_old=1.d0
			detGDsp2_up_new=1.d0
			detGDsp2_dw_new=1.d0
		CASE DEFAULT
			STOP 'Non hai selezionato un valore di Jsesp_kind accettabile &
			  [ module_calcola_accettazione.f90 > prima_valutazione_funzione_onda ]'
		END SELECT
		IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
			u_sesp1_new=u_sesp1_old
			u_sesp2_new=u_sesp2_old
		END IF
		Usesp1_new=Usesp1_old
		Usesp2_new=Usesp2_old
				
		IF (verbose_mode) PRINT *, 'calcola_accettazione: Ho reinizializzato la funzione d onda'
		
		!IF (verbose_mode) THEN
		!	IF (SDe_kind/='no_') THEN
		!		PRINT * , 'detSDe_up_new -> ', REAL(detSDe_up_old*DCONJG(detSDe_up_old))
		!		PRINT * , 'detSDe_dw_new -> ', REAL(detSDe_dw_old*DCONJG(detSDe_dw_old))
		!	END IF
		!	IF (Jee_kind/='no_') PRINT * , 'Uee -> ', Uee_old
		!	IF (Jep_kind/='no_') PRINT * , 'Uep -> ', Uep_old
		!	IF (SDse_kind/='no_') THEN
		!		PRINT * , 'detSDse1_up_new -> ', REAL(detSDse1_up_old*DCONJG(detSDse1_up_old))
		!		PRINT * , 'detSDse1_dw_new -> ', REAL(detSDse1_dw_old*DCONJG(detSDse1_dw_old))
		!		PRINT * , 'detSDse2_up_new -> ', REAL(detSDse2_up_old*DCONJG(detSDse2_up_old))
		!		PRINT * , 'detSDse2_dw_new -> ', REAL(detSDse2_dw_old*DCONJG(detSDse2_dw_old))
		!	END IF
		!	IF (Jse_kind/='no_') PRINT * , 'Use1 -> ', Use1_old
		!	IF (Jse_kind/='no_') PRINT * , 'Use2 -> ', Use2_old
		!	IF (Kse_kind/='no_') PRINT * , 'Kse1 -> ', Uese1_old
		!	IF (Kse_kind/='no_') PRINT * , 'Kse2 -> ', Uese2_old
		!	IF (Jsesp_kind/='no_') PRINT * , 'Usesp1 -> ', Usesp1_old
		!	IF (Jsesp_kind/='no_') PRINT * , 'Usesp2 -> ', Usesp2_old
		!END IF
		
	END SUBROUTINE prima_valutazione_funzione_onda
!-----------------------------------------------------------------------

	SUBROUTINE valuta_accettazione(num_mc,tipo,num,accettazione)
		USE dati_fisici
		USE dati_simulazione_mc
		USE funzione_onda
		USE momenta
		USE walkers
		IMPLICIT NONE
		LOGICAL :: flag_coppie
		CHARACTER(LEN=3) :: tipo
		INTEGER :: num
		INTEGER (KIND=8), PARAMETER :: n_coppie=1000
		INTEGER (KIND=8) :: num_mc, contatore_interno=0
		INTEGER (KIND=8), SAVE :: num_mc_old=0
		REAL (KIND=8) :: prob_acc, random
		LOGICAL, INTENT(OUT) :: accettazione
		INTEGER :: i, info, M_pvt(1:H_N_part), perm
		REAL (KIND=8) :: M(1:H_N_part,1:H_N_part), work(1:3*H_N_part), detM
				
		contatore_interno=contatore_interno+1
		
		flag_coppie=.FALSE.
		IF (num_mc/=num_mc_old) THEN
			IF (MOD(num_mc,n_coppie)==0) flag_coppie=.TRUE.
		END IF
		num_mc_old=num_mc
		
		IF (.NOT. iniz_calcola_accettazione) STOP 'Prima devi inizializzare la funzione d onda &
		  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
		IF (.NOT. iniz_dati_fisici) STOP 'Prima devi leggere i dati fisici &
		  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
		IF (.NOT. iniz_momenta) STOP 'Prima hai bisogno di trovare i momenti di fermi &
		  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
		IF (.NOT. iniz_walkers) STOP 'Prima devi inizializzare i walkers &
		  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				
		SELECT CASE (howtomove)
		CASE ('allp')
			!SDe
			IF ((tipo=='all') .OR. (tipo=='e__')) THEN
				SELECT CASE (SDe_kind)
				CASE ('pw_')
					CALL valuta_SD_pw(num,re_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,SDe_up_old,detSDe_up_old)
					CALL valuta_SD_pw(num,re_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
				CASE ('lda')
					CALL valuta_SD_lda(num,re_new(1:3,1:H_N_part),H_N_part, &
					  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,SDe_up_old,detSDe_up_old)
					CALL valuta_SD_lda(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
					  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
				CASE ('prf')
					CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
					  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,SDe_up_old,detSDe_up_old)
					CALL valuta_SD_har(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
					  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
				CASE ('fre')
					CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
					  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,SDe_up_old,detSDe_up_old)
					CALL valuta_SD_har(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
					  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
				CASE ('atm')
					CALL valuta_SD_atm(num,rij_ep_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,SDe_up_old,detSDe_up_old)
					CALL valuta_SD_atm(num,rij_ep_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
				CASE ('atp')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_p_')
					CALL valuta_SD_atm(num,rijpc_ep_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,SDe_up_old,detSDe_up_old)
					CALL valuta_SD_atm(num,rijpc_ep_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
				CASE ('no_')
					detSDe_up_new=1.d0
					detSDe_dw_new=1.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di SDe_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
			!Jee
			IF ((tipo=='all') .OR. (tipo=='e__')) THEN
				SELECT CASE (Jee_kind)
				CASE ('yuk')
					CALL valuta_Uee_YUK(num,rij_ee_new,N_part,u_ee_new,Uee_new)
				CASE ('yup')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_e_')
					CALL valuta_Uee_YUK(num,rijpc_ee_new,N_part,u_ee_new,Uee_new)
				CASE ('no_')
					Uee_old=0.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di Jee_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
			!Jep
			IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='p__')) THEN
				SELECT CASE (Jep_kind)
				CASE ('yuk')
					IF ((tipo=='all') .OR. (tipo=='e__')) THEN
						CALL valuta_Uep_YUK(num,1,rij_ep_new,N_part,u_ep_new,Uep_new)
					ELSE IF (tipo=='p__') THEN
						CALL valuta_Uep_YUK(num,2,rij_ep_new,N_part,u_ep_new,Uep_new)
					END IF
				CASE ('yup')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_p_')
					IF ((tipo=='all') .OR. (tipo=='e__')) THEN
						CALL valuta_Uep_YUK(num,1,rijpc_ep_new,N_part,u_ep_new,Uep_new)
					ELSE IF (tipo=='p__') THEN
						CALL valuta_Uep_YUK(num,2,rijpc_ep_new,N_part,u_ep_new,Uep_new)
					END IF
				CASE ('atm')
					IF ((tipo=='all') .OR. (tipo=='e__')) THEN
						CALL valuta_Uep_ATM(num,1,rij_ep_new,N_part,u_ep_new,Uep_new)
					ELSE IF (tipo=='p__') THEN
						CALL valuta_Uep_ATM(num,2,rij_ep_new,N_part,u_ep_new,Uep_new)
					END IF
				CASE ('atp')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_p_')
					IF ((tipo=='all') .OR. (tipo=='e__')) THEN
						CALL valuta_Uep_ATM(num,1,rijpc_ep_new,N_part,u_ep_new,Uep_new)
					ELSE IF (tipo=='p__') THEN
						CALL valuta_Uep_ATM(num,2,rijpc_ep_new,N_part,u_ep_new,Uep_new)
					END IF
				CASE ('no_') 
					Uep_new=0.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di Jep_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
			!SDse
			IF ((tipo=='all') .OR. (tipo=='se_')) THEN
				SELECT CASE (SDse_kind)
				CASE ('pw_')
					CALL valuta_SD_pw(num,se1_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,SDse1_up_old,detSDse1_up_old)
					CALL valuta_SD_pw(num,se1_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
					CALL valuta_SD_pw(num,se2_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
					CALL valuta_SD_pw(num,se2_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
				CASE ('pw2')
					CALL valuta_SD_pw(num,se1_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,SDse1_up_old,detSDse1_up_old)
					CALL valuta_SD_pw(num,se1_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
					CALL valuta_SD_pw(num,se2_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
					CALL valuta_SD_pw(num,se2_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
					  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
				CASE ('lda')
					CALL valuta_SD_lda(num,se1_new(1:3,1:H_N_part),H_N_part, &
					  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,SDse1_up_old,detSDse1_up_old)
					CALL valuta_SD_lda(num,se1_new(1:3,H_N_part+1:N_part),H_N_part, &
					  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
					CALL valuta_SD_lda(num,se2_new(1:3,1:H_N_part),H_N_part, &
					  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
					CALL valuta_SD_lda(num,se2_new(1:3,H_N_part+1:N_part),H_N_part, &
					  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
				CASE ('gem')
					CALL valuta_SD_gem(num,rij_se1_new(0,1:N_part,1:N_part),H_N_part, &
					  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,SDse1_up_old,detSDse1_up_old)
					CALL valuta_SD_gem(num,rij_se2_new(0,1:N_part,1:N_part),H_N_part, &
					  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
				CASE ('gss')
					CALL valuta_SD_gauss(num,rij_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,SDse1_up_old,detSDse1_up_old)
					CALL valuta_SD_gauss(num,rij_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
					CALL valuta_SD_gauss(num,rij_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
					CALL valuta_SD_gauss(num,rij_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
				CASE ('gsp')
					CALL calcola_nuove_distanze_pc(tipo,num,'sesp')
					CALL valuta_SD_gauss(num,rijpc_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,SDse1_up_old,detSDse1_up_old)
					CALL valuta_SD_gauss(num,rijpc_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
					CALL valuta_SD_gauss(num,rijpc_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
					CALL valuta_SD_gauss(num,rijpc_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
				CASE ('atm')
					CALL valuta_SD_atm(num,rij_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,SDse1_up_old,detSDse1_up_old)
					CALL valuta_SD_atm(num,rij_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
					CALL valuta_SD_atm(num,rij_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
					CALL valuta_SD_atm(num,rij_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
				CASE ('atp')
					CALL calcola_nuove_distanze_pc(tipo,num,'sesp')
					CALL valuta_SD_atm(num,rijpc_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,SDse1_up_old,detSDse1_up_old)
					CALL valuta_SD_atm(num,rijpc_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
					CALL valuta_SD_atm(num,rijpc_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
					  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
					CALL valuta_SD_atm(num,rijpc_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
				CASE ('no_')
					detSDse1_up_new=1.d0
					detSDse1_dw_new=1.d0
					detSDse2_up_new=1.d0
					detSDse2_dw_new=1.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di SDse_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
			!Jse
			IF ((tipo=='all') .OR. (tipo=='se_') .OR. ((tipo=='e__').AND.(flag_shadow))) THEN
				SELECT CASE (Jse_kind)
				CASE ('no_')
				CASE ('pot')
					CALL valuta_Use_POT(index_mol_num,dist_mol_ss1_new,u_POT_se1_new,Use1_new)
					CALL valuta_Use_POT(index_mol_num,dist_mol_ss2_new,u_POT_se2_new,Use2_new)
				CASE ('bou')
					!IF (flag_coppie) CALL seleziona_coppie_migliori(rij_ee_new,N_part,partner_se1)
					!CALL valuta_Use_BOUND(-1,rij_se1_new,N_part,partner_se1,b_se1_new,Bse1_new)
					!IF (flag_coppie) CALL seleziona_coppie_migliori(rij_ee_new,N_part,partner_se2)
					!CALL valuta_Use_BOUND(-1,rij_se2_new,N_part,partner_se2,b_se2_new,Bse2_new)
				CASE ('ppb')
					!CALL valuta_Use_POT(-1,rij_se1_new,N_part,u_se1_new,Use1_new)
					!CALL valuta_Use_POT(-1,rij_se2_new,N_part,u_se2_new,Use2_new)
					!IF (flag_coppie) CALL seleziona_coppie_migliori(rij_ee_new,N_part,partner_se1)
					!CALL valuta_Use_BOUND(-1,rij_se1_new,N_part,partner_se1,b_se1_new,Bse1_new)
					!IF (flag_coppie) CALL seleziona_coppie_migliori(rij_ee_new,N_part,partner_se2)
					!CALL valuta_Use_BOUND(-1,rij_se2_new,N_part,partner_se2,b_se2_new,Bse2_new)
				CASE ('yuk')
					CALL valuta_Usese_YUK(num,rij_se1_new,N_part,u_se1_new,Use1_new)
					CALL valuta_Usese_YUK(num,rij_se2_new,N_part,u_se2_new,Use2_new)
				CASE ('yup')
					CALL calcola_nuove_distanze_pc(tipo,num,'sese')
					CALL valuta_Usese_YUK(num,rijpc_se1_new,N_part,u_se1_new,Use1_new)
					CALL valuta_Usese_YUK(num,rijpc_se2_new,N_part,u_se2_new,Use2_new)
				CASE DEFAULT
					PRINT * , ' ### ', Jse_kind
					STOP 'Non hai selezionato un valore di Jse_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
			!Kse
			IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se_')) THEN
				SELECT CASE (Kse_kind)
				CASE ('gss')
					CALL valuta_KERNse(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
					CALL valuta_KERNse(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
				CASE ('gsc')
					CALL valuta_KERNse_ctf(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
					CALL valuta_KERNse_ctf(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
				CASE ('gsp')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_se')
					CALL valuta_KERNse(num,rijpc_ese1_new,N_part,u_ese1_new,Uese1_new)
					CALL valuta_KERNse(num,rijpc_ese2_new,N_part,u_ese2_new,Uese2_new)
				CASE ('gsd')
					CALL valuta_GDse(-1,0,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
					  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,GDse1_up_old,detGDse1_up_old)
					CALL valuta_GDse(-1,0,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,GDse1_dw_old,detGDse1_dw_old)
					CALL valuta_GDse(-1,0,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
					  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,GDse2_up_old,detGDse2_up_old)
					CALL valuta_GDse(-1,0,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,GDse2_dw_old,detGDse2_dw_old)
				CASE ('gdc')
					CALL valuta_GDse_ctf(-1,0,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
					  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,GDse1_up_old,detGDse1_up_old)
					CALL valuta_GDse_ctf(-1,0,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,GDse1_dw_old,detGDse1_dw_old)
					CALL valuta_GDse_ctf(-1,0,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
					  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,GDse2_up_old,detGDse2_up_old)
					CALL valuta_GDse_ctf(-1,0,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,GDse2_dw_old,detGDse2_dw_old)
				CASE ('gdp')
					CALL valuta_GDse(-1,0,rijpc_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
					  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,GDse1_up_old,detGDse1_up_old)
					CALL valuta_GDse(-1,0,rijpc_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,GDse1_dw_old,detGDse1_dw_old)
					CALL valuta_GDse(-1,0,rijpc_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
					  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,GDse2_up_old,detGDse2_up_old)
					CALL valuta_GDse(-1,0,rijpc_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,GDse2_dw_old,detGDse2_dw_old)
				CASE ('atm')
					CALL valuta_atmKERNse(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
					CALL valuta_atmKERNse(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
				CASE ('atc')
					CALL valuta_atmKERNse_ctf(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
					CALL valuta_atmKERNse_ctf(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
				CASE ('no_')
					Uese1_new=0.d0
					Uese2_new=0.d0
					detGDse1_up_new=1.d0
					detGDse1_dw_new=1.d0
					detGDse2_up_new=1.d0
					detGDse2_dw_new=1.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di Kse_kind accettabile &
					  [ module_calcola_accettazione.f90 > inizializza_funzione_onda ]'
				END SELECT
			END IF
			!Jsesp
			IF ((tipo=='all') .OR. (tipo=='se_')) THEN
				SELECT CASE (Jsesp_kind)
				CASE ('pot')
					CALL valuta_Usesp_POT(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
					CALL valuta_Usesp_POT(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
				CASE ('yuk')
					CALL valuta_Usesp_YUK(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
					CALL valuta_Usesp_YUK(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
				CASE ('yup')
					CALL calcola_nuove_distanze_pc(tipo,num,'sesp')
					CALL valuta_Usesp_YUK(num,1,rijpc_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
					CALL valuta_Usesp_YUK(num,1,rijpc_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
				CASE ('gss')
					CALL valuta_Usesp_GSS(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
					CALL valuta_Usesp_GSS(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
				CASE ('gsd')
					CALL valuta_GDsp(-1,0,rij_sesp1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
					  GDsp1_up_new,detGDsp1_up_new,IGDsp1_up_new,pvtgdsp1_up_new,GDsp1_up_old,detGDsp1_up_old)
					CALL valuta_GDsp(-1,0,rij_sesp1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  GDsp1_dw_new,detGDsp1_dw_new,IGDsp1_dw_new,pvtgdsp1_dw_new,GDsp1_dw_old,detGDsp1_dw_old)
					CALL valuta_GDsp(-1,0,rij_sesp2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
					  GDsp2_up_new,detGDsp2_up_new,IGDsp2_up_new,pvtgdsp2_up_new,GDsp2_up_old,detGDsp2_up_old)
					CALL valuta_GDsp(-1,0,rij_sesp2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
					  GDsp2_dw_new,detGDsp2_dw_new,IGDsp2_dw_new,pvtgdsp2_dw_new,GDsp2_dw_old,detGDsp2_dw_old)
				CASE ('no_') 
					Usesp1_old=0.d0
					Usesp2_old=0.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di Jsesp_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
		CASE ('1ppt')
			!SDe
			IF ((tipo=='e__').OR.(tipo=='tre')) THEN
				SELECT CASE (SDe_kind)
				CASE ('pw_')
					IF (num==-1) THEN
						CALL valuta_SD_pw(num,re_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
						CALL valuta_SD_pw(num,re_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_pw(num,re_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_pw(num-H_N_part,re_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					END IF
				CASE ('lda')
					IF (num==-1) THEN
						CALL valuta_SD_lda(num,re_new(1:3,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
						CALL valuta_SD_lda(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_lda(num,re_new(1:3,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_lda(num-H_N_part,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					END IF
				CASE ('prf')
					IF (num==-1) THEN
						CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,SDe_up_old,detSDe_up_old)
						CALL valuta_SD_har(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_har(num-H_N_part,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					END IF
				CASE ('fre')
					IF (num==-1) THEN
						CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,SDe_up_old,detSDe_up_old)
						CALL valuta_SD_har(num,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_har(num,re_new(1:3,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_har(num-H_N_part,re_new(1:3,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					END IF
				CASE ('atm')
					IF (num==-1) THEN
						CALL valuta_SD_atm(num,rij_ep_new(0,1:H_N_part,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
						CALL valuta_SD_atm(num,rij_ep_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_atm(num,rij_ep_new(0,1:H_N_part,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_atm(num-H_N_part,rij_ep_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					END IF
				CASE ('atp')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_p_')
					IF (num==-1) THEN
						CALL valuta_SD_atm(num,rijpc_ep_new(0,1:H_N_part,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
						CALL valuta_SD_atm(num,rijpc_ep_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL valuta_SD_atm(num,rijpc_ep_new(0,1:H_N_part,1:H_N_part),H_N_part, &
						  SDe_up_new,detSDe_up_new,ISDe_up_new,pvte_up_new,ISDe_up_old,detSDe_up_old)
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL valuta_SD_atm(num-H_N_part,rijpc_ep_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
						  SDe_dw_new,detSDe_dw_new,ISDe_dw_new,pvte_dw_new,ISDe_dw_old,detSDe_dw_old)
					END IF
				CASE ('no_')
					detSDe_up_new=1.d0
					detSDe_dw_new=1.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di SDe_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
			!Jee
			IF ((tipo=='e__').OR.(tipo=='tre')) THEN
				SELECT CASE (Jee_kind)
				CASE ('yuk')
					CALL valuta_Uee_YUK(num,rij_ee_new,N_part,u_ee_new,Uee_new)
				CASE ('yup')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_e_')
					CALL valuta_Uee_YUK(num,rijpc_ee_new,N_part,u_ee_new,Uee_new)
				CASE ('no_')
					Uee_old=0.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di Jee_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
			!Jep
			IF ((tipo=='e__') .OR. (tipo=='tre') .OR. (tipo=='p__')) THEN
				SELECT CASE (Jep_kind)
				CASE ('yuk')
					IF ((tipo=='e__') .OR. (tipo=='tre')) THEN
						CALL valuta_Uep_YUK(num,1,rij_ep_new,N_part,u_ep_new,Uep_new)
					ELSE IF (tipo=='p__') THEN
						CALL valuta_Uep_YUK(num,2,rij_ep_new,N_part,u_ep_new,Uep_new)
					END IF
				CASE ('yup')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_p_')
					IF ((tipo=='e__') .OR. (tipo=='tre')) THEN
						CALL valuta_Uep_YUK(num,1,rijpc_ep_new,N_part,u_ep_new,Uep_new)
					ELSE IF (tipo=='p__') THEN
						CALL valuta_Uep_YUK(num,2,rijpc_ep_new,N_part,u_ep_new,Uep_new)
					END IF
				CASE ('atm')
					IF ((tipo=='e__') .OR. (tipo=='tre')) THEN
						CALL valuta_Uep_ATM(num,1,rij_ep_new,N_part,u_ep_new,Uep_new)
					ELSE IF (tipo=='p__') THEN
						CALL valuta_Uep_ATM(num,2,rij_ep_new,N_part,u_ep_new,Uep_new)
					END IF
				CASE ('atp')
					CALL calcola_nuove_distanze_pc(tipo,num,'e_p_')
					IF ((tipo=='e__') .OR. (tipo=='tre')) THEN
						CALL valuta_Uep_ATM(num,1,rijpc_ep_new,N_part,u_ep_new,Uep_new)
					ELSE IF (tipo=='p__') THEN
						CALL valuta_Uep_ATM(num,2,rijpc_ep_new,N_part,u_ep_new,Uep_new)
					END IF
				CASE ('no_') 
					Uep_new=0.d0
				CASE DEFAULT
					STOP 'Non hai selezionato un valore di Jep_kind accettabile &
					  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
				END SELECT
			END IF
			IF (flag_shadow) THEN
				!SDse
				IF ((tipo=='se1') .OR. (tipo=='tre')) THEN
					SELECT CASE (SDse_kind)
					CASE ('pw_')
						IF (num==-1) THEN
							CALL valuta_SD_pw(num,se1_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
							CALL valuta_SD_pw(num,se1_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_pw(num,se1_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_pw(num-H_N_part,se1_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						END IF
					CASE ('pw2')
						IF (num==-1) THEN
							CALL valuta_SD_pw(num,se1_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
							CALL valuta_SD_pw(num,se1_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_pw(num,se1_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_pw(num-H_N_part,se1_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						END IF
					CASE ('lda')
						IF (num==-1) THEN
							CALL valuta_SD_lda(num,se1_new(1:3,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
							CALL valuta_SD_lda(num,se1_new(1:3,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_lda(num,se1_new(1:3,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_lda(num-H_N_part,se1_new(1:3,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						END IF
					CASE ('gem')
						IF (num==-1) THEN
							CALL valuta_SD_gem(num,rij_se1_new(0,1:N_part,1:N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						ELSE IF ((num>0) .AND. (num<=N_part)) THEN
							CALL valuta_SD_gem(num,rij_se1_new(0,1:N_part,1:N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						END IF
					CASE ('gss')
						IF (num==-1) THEN
							CALL valuta_SD_gauss(num,rij_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
							CALL valuta_SD_gauss(num,rij_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_gauss(num,rij_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_gauss(num-H_N_part,rij_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						END IF
					CASE ('gsp')
						CALL calcola_nuove_distanze_pc(tipo,num,'sesp')
						IF (num==-1) THEN
							CALL valuta_SD_gauss(num,rijpc_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
							CALL valuta_SD_gauss(num,rijpc_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_gauss(num,rijpc_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_gauss(num-H_N_part,rijpc_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						END IF
					CASE ('atm')
						IF (num==-1) THEN
							CALL valuta_SD_atm(num,rij_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
							CALL valuta_SD_atm(num,rij_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_atm(num,rij_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_atm(num-H_N_part,rij_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						END IF
					CASE ('atp')
						CALL calcola_nuove_distanze_pc(tipo,num,'sesp')
						IF (num==-1) THEN
							CALL valuta_SD_atm(num,rijpc_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
							CALL valuta_SD_atm(num,rijpc_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_atm(num,rijpc_sesp1_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse1_up_new,detSDse1_up_new,ISDse1_up_new,pvtse1_up_new,ISDse1_up_old,detSDse1_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_atm(num-H_N_part,rijpc_sesp1_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new,pvtse1_dw_new,ISDse1_dw_old,detSDse1_dw_old)
						END IF
					CASE ('no_')
						detSDse1_up_new=1.d0
						detSDse1_dw_new=1.d0
					CASE DEFAULT
						STOP 'Non hai selezionato un valore di SDse1_kind accettabile &
						  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
					END SELECT
				END IF
				IF ((tipo=='se2') .OR. (tipo=='tre')) THEN
					SELECT CASE (SDse_kind)
					CASE ('pw_')
						IF (num==-1) THEN
							CALL valuta_SD_pw(num,se2_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
							CALL valuta_SD_pw(num,se2_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_pw(num,se2_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_pw(num-H_N_part,se2_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						END IF
					CASE ('pw2')
						IF (num==-1) THEN
							CALL valuta_SD_pw(num,se2_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
							CALL valuta_SD_pw(num,se2_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_pw(num,se2_new(1:3,1:H_N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_pw(num-H_N_part,se2_new(1:3,H_N_part+1:N_part),H_N_part,k_pw(1:3,1:H_N_part), &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						END IF
					CASE ('lda')
						IF (num==-1) THEN
							CALL valuta_SD_lda(num,se2_new(1:3,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
							CALL valuta_SD_lda(num,se2_new(1:3,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_lda(num,se2_new(1:3,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_lda(num-H_N_part,se2_new(1:3,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						END IF
					CASE ('gem')
						IF (num==-1) THEN
							CALL valuta_SD_gem(num,rij_se2_new(0,1:N_part,1:N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,SDse2_up_old,detSDse2_up_old)
						ELSE IF ((num>0) .AND. (num<=N_part)) THEN
							CALL valuta_SD_gem(num,rij_se2_new(0,1:N_part,1:N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
						END IF
					CASE ('gss')
						IF (num==-1) THEN
							CALL valuta_SD_gauss(num,rij_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
							CALL valuta_SD_gauss(num,rij_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_gauss(num,rij_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_gauss(num-H_N_part,rij_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						END IF
					CASE ('gsp')
						CALL calcola_nuove_distanze_pc(tipo,num,'sesp')
						IF (num==-1) THEN
							CALL valuta_SD_gauss(num,rijpc_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
							CALL valuta_SD_gauss(num,rijpc_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_gauss(num,rijpc_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_gauss(num-H_N_part,rijpc_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						END IF
					CASE ('atm')
						IF (num==-1) THEN
							CALL valuta_SD_atm(num,rij_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
							CALL valuta_SD_atm(num,rij_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_atm(num,rij_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_atm(num-H_N_part,rij_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						END IF
					CASE ('atp')
						CALL calcola_nuove_distanze_pc(tipo,num,'sesp')
						IF (num==-1) THEN
							CALL valuta_SD_atm(num,rijpc_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
							CALL valuta_SD_atm(num,rijpc_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
							CALL valuta_SD_atm(num,rijpc_sesp2_new(0,1:H_N_part,1:H_N_part),H_N_part, &
							  SDse2_up_new,detSDse2_up_new,ISDse2_up_new,pvtse2_up_new,ISDse2_up_old,detSDse2_up_old)
						ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
							CALL valuta_SD_atm(num-H_N_part,rijpc_sesp2_new(0,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
							  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new,pvtse2_dw_new,ISDse2_dw_old,detSDse2_dw_old)
						END IF
					CASE ('no_')
						detSDse2_up_new=1.d0
						detSDse2_dw_new=1.d0
					CASE DEFAULT
						STOP 'Non hai selezionato un valore di SDse2_kind accettabile &
						  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
					END SELECT
				END IF
				!Jse
				IF ((tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
					SELECT CASE (Jse_kind)
					CASE ('no_')
					CASE ('pot')
						IF (tipo=='se1') THEN
							CALL valuta_Use_POT(index_mol_num,dist_mol_ss1_new,u_POT_se1_new,Use1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_Use_POT(index_mol_num,dist_mol_ss2_new,u_POT_se2_new,Use2_new)
						ELSE IF (tipo=='tre') THEN
							CALL valuta_Use_POT(index_mol_num,dist_mol_ss1_new,u_POT_se1_new,Use1_new)
							CALL valuta_Use_POT(index_mol_num,dist_mol_ss2_new,u_POT_se2_new,Use2_new)
						END IF
					CASE ('bou')
						!IF (tipo=='se1') THEN
						!	IF (flag_coppie) CALL seleziona_coppie_migliori(rij_ee_new,N_part,partner_se1)
						!	CALL valuta_Use_BOUND(num,rij_se1_new,N_part,partner_se1,b_se1_new,Bse1_new)
						!ELSE IF (tipo=='se2') THEN
						!	IF (flag_coppie) CALL seleziona_coppie_migliori(rij_ee_new,N_part,partner_se2)
						!	CALL valuta_Use_BOUND(num,rij_se2_new,N_part,partner_se2,b_se2_new,Bse2_new)
						!END IF
					CASE ('ppb')
						!IF (tipo=='se1') THEN
						!	CALL valuta_Use_POT(num,rij_se1_new,u_se1_new,Use1_new)
						!	IF (flag_coppie) CALL seleziona_coppie_migliori(rij_ee_new,N_part,partner_se1)
						!	CALL valuta_Use_BOUND(num,rij_se1_new,N_part,partner_se1,b_se1_new,Bse1_new)
						!ELSE IF (tipo=='se2') THEN
						!	CALL valuta_Use_POT(num,rij_se2_new,u_se2_new,Use2_new)
						!	IF (flag_coppie) CALL seleziona_coppie_migliori(rij_ee_new,N_part,partner_se2)
						!	CALL valuta_Use_BOUND(num,rij_se2_new,N_part,partner_se2,b_se2_new,Bse2_new)
						!END IF
					CASE ('yuk')
						IF (tipo=='se1') THEN
							CALL valuta_Usese_YUK(num,rij_se1_new,N_part,u_se1_new,Use1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_Usese_YUK(num,rij_se2_new,N_part,u_se2_new,Use2_new)
						ELSE IF (tipo=='tre') THEN
							CALL valuta_Usese_YUK(num,rij_se1_new,N_part,u_se1_new,Use1_new)
							CALL valuta_Usese_YUK(num,rij_se2_new,N_part,u_se2_new,Use2_new)
						END IF
					CASE ('yup')
						CALL calcola_nuove_distanze_pc(tipo,num,'sese')
						IF (tipo=='se1') THEN
							CALL valuta_Usese_YUK(num,rijpc_se1_new,N_part,u_se1_new,Use1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_Usese_YUK(num,rijpc_se2_new,N_part,u_se2_new,Use2_new)
						ELSE IF (tipo=='tre') THEN
							CALL valuta_Usese_YUK(num,rijpc_se1_new,N_part,u_se1_new,Use1_new)
							CALL valuta_Usese_YUK(num,rijpc_se2_new,N_part,u_se2_new,Use2_new)
						END IF
					CASE DEFAULT
						PRINT * , ' ### ', Jse_kind
						STOP 'Non hai selezionato un valore di Jse_kind accettabile &
						  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
					END SELECT
				END IF
				!Kse
				IF ( (tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
					SELECT CASE (Kse_kind)
					CASE ('gss')
						IF ((tipo=='e__') .OR. (tipo=='all')) THEN
							CALL valuta_KERNse(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
							CALL valuta_KERNse(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
						ELSE IF (tipo=='se1') THEN
							CALL valuta_KERNse(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_KERNse(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
						ELSE IF (tipo=='tre') THEN
							Uese1_new=Uese1_old
							Uese2_new=Uese2_old
						END IF
					CASE ('gsc')
						IF ((tipo=='e__') .OR. (tipo=='all')) THEN
							CALL valuta_KERNse_ctf(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
							CALL valuta_KERNse_ctf(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
						ELSE IF (tipo=='se1') THEN
							CALL valuta_KERNse_ctf(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_KERNse_ctf(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
						ELSE IF (tipo=='tre') THEN
							Uese1_new=Uese1_old
							Uese2_new=Uese2_old
						END IF
					CASE ('gsp')
						CALL calcola_nuove_distanze_pc(tipo,num,'e_se')
						IF ((tipo=='e__') .OR. (tipo=='all')) THEN
							CALL valuta_KERNse(num,rijpc_ese1_new,N_part,u_ese1_new,Uese1_new)
							CALL valuta_KERNse(num,rijpc_ese2_new,N_part,u_ese2_new,Uese2_new)
						ELSE IF (tipo=='se1') THEN
							CALL valuta_KERNse(num,rijpc_ese1_new,N_part,u_ese1_new,Uese1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_KERNse(num,rijpc_ese2_new,N_part,u_ese2_new,Uese2_new)
						ELSE IF (tipo=='tre') THEN
							Uese1_new=Uese1_old
							Uese2_new=Uese2_old
						END IF
					CASE ('gsd')
						IF ((tipo=='e__') .OR. (tipo=='all')) THEN
							IF (num==-1) THEN
								CALL valuta_GDse(-1,1,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse(-1,1,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
								CALL valuta_GDse(-1,1,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
								CALL valuta_GDse(-1,1,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse(num,1,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse(num,1,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse(num-H_N_part,1,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
								CALL valuta_GDse(num-H_N_part,1,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							END IF
						ELSE IF (tipo=='se1') THEN
							IF (num==-1) THEN
								CALL valuta_GDse(-1,2,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse(-1,2,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse(num,2,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse(num-H_N_part,2,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
							END IF
						ELSE IF (tipo=='se2') THEN
							IF (num==-1) THEN
								CALL valuta_GDse(-1,2,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
								CALL valuta_GDse(-1,2,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse(num,2,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse(num-H_N_part,2,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							END IF
						ELSE IF (tipo=='tre') THEN
							!TODO
						END IF
					CASE ('gdc')
						IF ((tipo=='e__') .OR. (tipo=='all')) THEN
							IF (num==-1) THEN
								CALL valuta_GDse_ctf(-1,1,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse_ctf(-1,1,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
								CALL valuta_GDse_ctf(-1,1,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
								CALL valuta_GDse_ctf(-1,1,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse_ctf(num,1,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse_ctf(num,1,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse_ctf(num-H_N_part,1,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
								CALL valuta_GDse_ctf(num-H_N_part,1,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							END IF
						ELSE IF (tipo=='se1') THEN
							IF (num==-1) THEN
								CALL valuta_GDse_ctf(-1,2,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse_ctf(-1,2,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse_ctf(num,2,rij_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse_ctf(num-H_N_part,2,rij_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
							END IF
						ELSE IF (tipo=='se2') THEN
							IF (num==-1) THEN
								CALL valuta_GDse_ctf(-1,2,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
								CALL valuta_GDse_ctf(-1,2,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse_ctf(num,2,rij_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse_ctf(num-H_N_part,2,rij_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							END IF
						ELSE IF (tipo=='tre') THEN
							!TODO
						END IF
					CASE ('gdp')
						CALL calcola_nuove_distanze_pc(tipo,num,'e_se')
						IF ((tipo=='e__') .OR. (tipo=='all')) THEN
							IF (num==-1) THEN
								CALL valuta_GDse(-1,1,rijpc_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse(-1,1,rijpc_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
								CALL valuta_GDse(-1,1,rijpc_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
								CALL valuta_GDse(-1,1,rijpc_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse(num,1,rijpc_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse(num,1,rijpc_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse(num-H_N_part,1,rijpc_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
								CALL valuta_GDse(num-H_N_part,1,rijpc_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							END IF
						ELSE IF (tipo=='se1') THEN
							IF (num==-1) THEN
								CALL valuta_GDse(-1,2,rijpc_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
								CALL valuta_GDse(-1,2,rijpc_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse(num,2,rijpc_ese1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse1_up_new,detGDse1_up_new,IGDse1_up_new,pvtgdse1_up_new,IGDse1_up_old,detGDse1_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse(num-H_N_part,2,rijpc_ese1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new,pvtgdse1_dw_new,IGDse1_dw_old,detGDse1_dw_old)
							END IF
						ELSE IF (tipo=='se2') THEN
							IF (num==-1) THEN
								CALL valuta_GDse(-1,2,rijpc_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
								CALL valuta_GDse(-1,2,rijpc_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDse(num,2,rijpc_ese2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDse2_up_new,detGDse2_up_new,IGDse2_up_new,pvtgdse2_up_new,IGDse2_up_old,detGDse2_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDse(num-H_N_part,2,rijpc_ese2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new,pvtgdse2_dw_new,IGDse2_dw_old,detGDse2_dw_old)
							END IF
						ELSE IF (tipo=='tre') THEN
							!TODO
						END IF
					CASE ('atm')   !OUTDATED
						IF ((tipo=='e__') .OR. (tipo=='all')) THEN
							CALL valuta_atmKERNse(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
							CALL valuta_atmKERNse(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
						ELSE IF (tipo=='se1') THEN
							CALL valuta_atmKERNse(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_atmKERNse(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
						END IF
					CASE ('atc')
						IF ((tipo=='e__') .OR. (tipo=='all')) THEN
							CALL valuta_atmKERNse_ctf(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
							CALL valuta_atmKERNse_ctf(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
						ELSE IF (tipo=='se1') THEN
							CALL valuta_atmKERNse_ctf(num,rij_ese1_new,N_part,u_ese1_new,Uese1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_atmKERNse_ctf(num,rij_ese2_new,N_part,u_ese2_new,Uese2_new)
						END IF
					CASE ('no_')
						Uese1_new=0.d0
						Uese2_new=0.d0
						detGDse1_up_new=1.d0
						detGDse1_dw_new=1.d0
						detGDse2_up_new=1.d0
						detGDse2_dw_new=1.d0
					CASE DEFAULT
						STOP 'Non hai selezionato un valore di Kse_kind accettabile &
						  [ module_calcola_accettazione.f90 > inizializza_funzione_onda ]'
					END SELECT
				END IF
				!Jsesp
				IF ((tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
					SELECT CASE (Jsesp_kind)
					CASE ('pot')
						IF (tipo=='se1') THEN
							CALL valuta_Usesp_POT(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_Usesp_POT(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
						ELSE IF (tipo=='tre') THEN
							CALL valuta_Usesp_POT(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
							CALL valuta_Usesp_POT(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
						END IF
					CASE ('yuk')
						IF (tipo=='se1') THEN
							CALL valuta_Usesp_YUK(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_Usesp_YUK(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
						ELSE IF (tipo=='tre') THEN
							CALL valuta_Usesp_YUK(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
							CALL valuta_Usesp_YUK(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
						END IF
					CASE ('yup')
						CALL calcola_nuove_distanze_pc(tipo,num,'sesp')
						IF (tipo=='se1') THEN
							CALL valuta_Usesp_YUK(num,1,rijpc_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_Usesp_YUK(num,1,rijpc_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
						ELSE IF (tipo=='tre') THEN
							CALL valuta_Usesp_YUK(num,1,rijpc_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
							CALL valuta_Usesp_YUK(num,1,rijpc_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
						END IF
					CASE ('gss')
						IF (tipo=='se1') THEN
							CALL valuta_Usesp_GSS(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
						ELSE IF (tipo=='se2') THEN
							CALL valuta_Usesp_GSS(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
						ELSE IF (tipo=='tre') THEN
							CALL valuta_Usesp_GSS(num,1,rij_sesp1_new,N_part,u_sesp1_new,Usesp1_new)
							CALL valuta_Usesp_GSS(num,1,rij_sesp2_new,N_part,u_sesp2_new,Usesp2_new)
						END IF
					CASE ('gsd')   !OUTDATED
						IF ((tipo=='se_') .OR. (tipo=='all')) THEN
							IF (num==-1) THEN
								CALL valuta_GDsp(-1,0,rij_sesp1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDsp1_up_new,detGDsp1_up_new,IGDsp1_up_new,pvtgdsp1_up_new,IGDsp1_up_old,detGDsp1_up_old)
								CALL valuta_GDsp(-1,0,rij_sesp1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDsp1_dw_new,detGDsp1_dw_new,IGDsp1_dw_new,pvtgdsp1_dw_new,IGDsp1_dw_old,detGDsp1_dw_old)
								CALL valuta_GDsp(-1,0,rij_sesp2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDsp2_up_new,detGDsp2_up_new,IGDsp2_up_new,pvtgdsp2_up_new,IGDsp2_up_old,detGDsp2_up_old)
								CALL valuta_GDsp(-1,0,rij_sesp2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDsp2_dw_new,detGDsp2_dw_new,IGDsp2_dw_new,pvtgdsp2_dw_new,IGDsp2_dw_old,detGDsp2_dw_old)
							END IF
						ELSE IF (tipo=='se1') THEN
							IF (num==-1) THEN
								CALL valuta_GDsp(-1,0,rij_sesp1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDsp1_up_new,detGDsp1_up_new,IGDsp1_up_new,pvtgdsp1_up_new,IGDsp1_up_old,detGDsp1_up_old)
								CALL valuta_GDsp(-1,0,rij_sesp1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDsp1_dw_new,detGDsp1_dw_new,IGDsp1_dw_new,pvtgdsp1_dw_new,IGDsp1_dw_old,detGDsp1_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDsp(num,1,rij_sesp1_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDsp1_up_new,detGDsp1_up_new,IGDsp1_up_new,pvtgdsp1_up_new,IGDsp1_up_old,detGDsp1_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDsp(num-H_N_part,1,rij_sesp1_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDsp1_dw_new,detGDsp1_dw_new,IGDsp1_dw_new,pvtgdsp1_dw_new,IGDsp1_dw_old,detGDsp1_dw_old)
							END IF
						ELSE IF (tipo=='se2') THEN
							IF (num==-1) THEN
								CALL valuta_GDsp(-1,1,rij_sesp2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDsp2_up_new,detGDsp2_up_new,IGDsp2_up_new,pvtgdsp2_up_new,IGDsp2_up_old,detGDsp2_up_old)
								CALL valuta_GDsp(-1,1,rij_sesp2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDsp2_dw_new,detGDsp2_dw_new,IGDsp2_dw_new,pvtgdsp2_dw_new,IGDsp2_dw_old,detGDsp2_dw_old)
							ELSE IF ((num>0).AND.(num<=H_N_part)) THEN
								CALL valuta_GDsp(num,1,rij_sesp2_new(0:3,1:H_N_part,1:H_N_part),H_N_part, &
								  GDsp2_up_new,detGDsp2_up_new,IGDsp2_up_new,pvtgdsp2_up_new,IGDsp2_up_old,detGDsp2_up_old)
							ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
								CALL valuta_GDsp(num-H_N_part,1,rij_sesp2_new(0:3,H_N_part+1:N_part,H_N_part+1:N_part),H_N_part, &
								  GDsp2_dw_new,detGDsp2_dw_new,IGDsp2_dw_new,pvtgdsp2_dw_new,IGDsp2_dw_old,detGDsp2_dw_old)
							END IF
						END IF
					CASE ('no_') 
						Usesp1_old=0.d0
						Usesp2_old=0.d0
					CASE DEFAULT
						STOP 'Non hai selezionato un valore di Jsesp_kind accettabile &
						  [ module_calcola_accettazione.f90 > valuta_accettazione ]'
					END SELECT
				END IF
			END IF
		END SELECT
		
		IF (SDse_kind/='pw2') THEN
			prob_acc=REAL( DEXP(-(Uee_new-Uee_old)-(Uep_new-Uep_old)-0.5*(Use1_new-Use1_old)-0.5*(Use2_new-Use2_old)&
			  -(Bse1_new-Bse1_old)-(Bse2_new-Bse2_old)-(Uese1_new-Uese1_old)-(Uese2_new-Uese2_old)-&
			  0.5*(Usesp1_new-Usesp1_old)-0.5*(Usesp2_new-Usesp2_old)) * &
			  ( (detSDe_up_new*DCONJG(detSDe_up_new)/(detSDe_up_old*DCONJG(detSDe_up_old)))* &
			  (detSDe_dw_new*DCONJG(detSDe_dw_new)/(detSDe_dw_old*DCONJG(detSDe_dw_old))) )* &
			  DABS( REAL((detSDse1_up_new/detSDse1_up_old)*(DCONJG(detSDse2_up_new)/DCONJG(detSDse2_up_old))* &
			  (detSDse1_dw_new/detSDse1_dw_old)*(DCONJG(detSDse2_dw_new)/DCONJG(detSDse2_dw_old)),8) ) * &
			  DABS( ((detGDse1_up_new/detGDse1_up_old)*((detGDse2_up_new)/(detGDse2_up_old))* &
			  (detGDse1_dw_new/detGDse1_dw_old)*((detGDse2_dw_new)/(detGDse2_dw_old))) )* &
			  DABS( ((detGDsp1_up_new/detGDsp1_up_old)*((detGDsp2_up_new)/(detGDsp2_up_old))* &
			  (detGDsp1_dw_new/detGDsp1_dw_old)*((detGDsp2_dw_new)/(detGDsp2_dw_old))) ) ,8)
		ELSE
			prob_acc=REAL( DEXP(-(Uee_new-Uee_old)-(Uep_new-Uep_old)-0.5*(Use1_new-Use1_old)-0.5*(Use2_new-Use2_old)&
			  -(Bse1_new-Bse1_old)-(Bse2_new-Bse2_old)-(Uese1_new-Uese1_old)-(Uese2_new-Uese2_old)-&
			  0.5*(Usesp1_new-Usesp1_old)-0.5*(Usesp2_new-Usesp2_old)) * &
			  ( (detSDe_up_new*DCONJG(detSDe_up_new)/(detSDe_up_old*DCONJG(detSDe_up_old)))* &
			  (detSDe_dw_new*DCONJG(detSDe_dw_new)/(detSDe_dw_old*DCONJG(detSDe_dw_old))) )* &
			  DABS( REAL(((detSDse1_up_new**2)/(detSDse1_up_old**2))*(DCONJG((detSDse2_up_new**2))/DCONJG((detSDse2_up_old**2)))* &
			  ((detSDse1_dw_new**2)/(detSDse1_dw_old**2))*(DCONJG((detSDse2_dw_new**2))/DCONJG((detSDse2_dw_old**2))),8) ) * &
			  DABS( ((detGDse1_up_new/detGDse1_up_old)*((detGDse2_up_new)/(detGDse2_up_old))* &
			  (detGDse1_dw_new/detGDse1_dw_old)*((detGDse2_dw_new)/(detGDse2_dw_old))) )* &
			  DABS( ((detGDsp1_up_new/detGDsp1_up_old)*((detGDsp2_up_new)/(detGDsp2_up_old))* &
			  (detGDsp1_dw_new/detGDsp1_dw_old)*((detGDsp2_dw_new)/(detGDsp2_dw_old))) ) ,8)
		END IF
		
		CALL RANDOM_NUMBER(random)
				
		!prob_acc=REAL( DEXP(-(Uee_new-Uee_old)-(Uep_new-Uep_old)-(Use1_new-Use1_old)-(Use2_new-Use2_old)-(Bse1_new-Bse1_old)- &
		!  (Bse2_new-Bse2_old)-(Uese1_new-Uese1_old)-(Uese2_new-Uese2_old)-(Usesp1_new-Usesp1_old)-(Usesp2_new-Usesp2_old)) * &
		!  ( (detSDe_up_new*DCONJG(detSDe_up_new)/(detSDe_up_old*DCONJG(detSDe_up_old)))* &
		!  (detSDe_dw_new*DCONJG(detSDe_dw_new)/(detSDe_dw_old*DCONJG(detSDe_dw_old))) )* &
		!  DABS( REAL((detSDse1_up_new/detSDse1_up_old)*(DCONJG(detSDse2_up_new)/DCONJG(detSDse2_up_old))* &
		!  (detSDse1_dw_new/detSDse1_dw_old)*(DCONJG(detSDse2_dw_new)/DCONJG(detSDse2_dw_old)),8) )  &
		!   ,8)
		
		!IF (verbose_mode) PRINT * , '-------------------------------------------'
		!IF (verbose_mode) PRINT * , 'prob_acc -> ', prob_acc, '   :::'
		!IF (verbose_mode) THEN
		!IF (mpi_myrank==0) THEN
		!IF (MOD(contatore_interno,100)==0) THEN 
			!PRINT * , '-------------------------------------------'
			!PRINT * , 'tipo: ', tipo, '           num= ', num
			!PRINT * , 'prob_acc -> ', prob_acc, '   :::'
			!IF (Jee_kind/='no_') PRINT * , 'Uee -> ', Uee_old, Uee_new, ' : ', DEXP(-Uee_new+Uee_old)
			!IF (Jep_kind/='no_') PRINT * , 'Uep -> ', Uep_old, Uep_new
			!IF (SDe_kind/='no_') THEN
				!PRINT * , 'detSDe_up_new -> ', REAL(detSDe_up_old*DCONJG(detSDe_up_old)), &
				!  REAL(detSDe_up_new*DCONJG(detSDe_up_new))
				!PRINT * , 'detSDe_dw_new -> ', REAL(detSDe_dw_old*DCONJG(detSDe_dw_old)), &
				!  REAL(detSDe_dw_new*DCONJG(detSDe_dw_new))
				!PRINT * , 'detSDe_up_new -> ', REAL(detSDe_up_old), REAL(detSDe_up_new)
				!PRINT * , 'detSDe_dw_new -> ', REAL(detSDe_dw_old), REAL(detSDe_dw_new)
			!END IF
			!IF (Jse_kind/='no_') PRINT * , 'Use1 -> ', Use1_old, Use1_new
			!IF (Jse_kind/='no_') PRINT * , 'Use2 -> ', Use2_old, Use2_new
			!IF (Jse_kind/='no_') PRINT * , 'Bse1 -> ', Bse1_old, Bse1_new
			!IF (Jse_kind/='no_') PRINT * , 'Bse2 -> ', Bse2_old, Bse2_new
			!IF (Kse_kind/='no_') PRINT * , 'Kse1 -> ', Uese1_old, Uese1_new
			!IF (Kse_kind/='no_') PRINT * , 'Kse2 -> ', Uese2_old, Uese2_new
			!IF (Jsesp_kind/='no_') PRINT * , 'Usesp1 -> ', Usesp1_old, Usesp1_new
			!IF (Jsesp_kind/='no_') PRINT * , 'Usesp2 -> ', Usesp2_old, Usesp2_new
			!PRINT * , tipo, num
			!IF (SDse_kind/='no_') THEN
			!	PRINT * , 'detSDse1_up_new -> ', REAL(detSDse1_up_old*DCONJG(detSDse1_up_old)), &
			!	  REAL(detSDse1_up_new*DCONJG(detSDse1_up_new))
			!	PRINT * , 'detSDse1_dw_new -> ', REAL(detSDse1_dw_old*DCONJG(detSDse1_dw_old)), &
			!	  REAL(detSDse1_dw_new*DCONJG(detSDse1_dw_new))
			!	PRINT * , 'detSDse2_up_new -> ', REAL(detSDse2_up_old*DCONJG(detSDse2_up_old)), &
			!	  REAL(detSDse2_up_new*DCONJG(detSDse2_up_new))
			!	PRINT * , 'detSDse2_dw_new -> ', REAL(detSDse2_dw_old*DCONJG(detSDse2_dw_old)), &
			!	  REAL(detSDse2_dw_new*DCONJG(detSDse2_dw_new))
			!	PRINT * , 'detSDse1_up_new -> ', REAL(detSDse1_up_old), REAL(detSDse1_up_new)
			!	PRINT * , 'detSDse1_dw_new -> ', REAL(detSDse1_dw_old), REAL(detSDse1_dw_new)
			!	PRINT * , 'detSDse2_up_new -> ', REAL(detSDse2_up_old), REAL(detSDse2_up_new)
			!	PRINT * , 'detSDse2_dw_new -> ', REAL(detSDse2_dw_old), REAL(detSDse2_dw_new)
			!END IF
			!IF (prob_acc>random) THEN
			!	PRINT * , 'ACCETTATO'
			!ELSE
			!	PRINT * , 'RIFIUTATO'
			!END IF	
			!PRINT * , '-------------------------------------------'
		!END IF
		!END IF
		!END IF
		
		IF (prob_acc>random) THEN
			accettazione=.TRUE.
			!PRINT * , 'PROB_ACC=', prob_acc, '     random=', random, '     ACCETTATO'
			!prob_acc=REAL( DEXP(-(Uee_new-Uee_old)-(Uep_new-Uep_old)-(Use1_new-Use1_old)-(Use2_new-Use2_old)-(Bse1_new-Bse1_old)- &
			!  (Bse2_new-Bse2_old)-(Uese1_new-Uese1_old)-(Uese2_new-Uese2_old)-(Usesp1_new-Usesp1_old)-(Usesp2_new-Usesp2_old)) * &
			!  ( (detSDe_up_new*DCONJG(detSDe_up_new)/(detSDe_up_old*DCONJG(detSDe_up_old)))* &
			!  (detSDe_dw_new*DCONJG(detSDe_dw_new)/(detSDe_dw_old*DCONJG(detSDe_dw_old))) )* &
			!  DABS( REAL((detSDse1_up_new/detSDse1_up_old)*(DCONJG(detSDse2_up_new)/DCONJG(detSDse2_up_old))* &
			!  (detSDse1_dw_new/detSDse1_dw_old)*(DCONJG(detSDse2_dw_new)/DCONJG(detSDse2_dw_old)),8) )  &
			!   ,8)
			!PRINT * , 'PROB_ACC_SENZA_GD=', prob_acc
			CALL aggiorna_funzione_onda(tipo,num)
		ELSE
			accettazione=.FALSE.
			!PRINT * , 'PROB_ACC=', prob_acc, '     random=', random, '     RIFIUTATO'
			!prob_acc=REAL( DEXP(-(Uee_new-Uee_old)-(Uep_new-Uep_old)-(Use1_new-Use1_old)-(Use2_new-Use2_old)-(Bse1_new-Bse1_old)- &
			!  (Bse2_new-Bse2_old)-(Uese1_new-Uese1_old)-(Uese2_new-Uese2_old)-(Usesp1_new-Usesp1_old)-(Usesp2_new-Usesp2_old)) * &
			!  ( (detSDe_up_new*DCONJG(detSDe_up_new)/(detSDe_up_old*DCONJG(detSDe_up_old)))* &
			!  (detSDe_dw_new*DCONJG(detSDe_dw_new)/(detSDe_dw_old*DCONJG(detSDe_dw_old))) )* &
			!  DABS( REAL((detSDse1_up_new/detSDse1_up_old)*(DCONJG(detSDse2_up_new)/DCONJG(detSDse2_up_old))* &
			!  (detSDse1_dw_new/detSDse1_dw_old)*(DCONJG(detSDse2_dw_new)/DCONJG(detSDse2_dw_old)),8) )  &
			!   ,8)
			!PRINT * , 'PROB_ACC_SENZA_GD=', prob_acc
			CALL roll_back(tipo,num)
		END IF
				
		!IF (verbose_mode) PRINT *, 'calcola_accettazione: Ho calcolato l accettazione, ', accettazione
				
	END SUBROUTINE valuta_accettazione
!-----------------------------------------------------------------------

	SUBROUTINE roll_back(tipo,num)
		USE dati_fisici
		USE funzione_onda
		USE walkers
		IMPLICIT NONE
		CHARACTER(LEN=3) :: tipo
		INTEGER, INTENT(IN) :: num
				
		IF (num==-1) THEN
			!SDe
			IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='tre')) THEN
				IF (SDe_kind/='no_') THEN
					SDe_up_new=SDe_up_old
					ISDe_up_new=ISDe_up_old
					pvte_up_new=pvte_up_old
					SDe_dw_new=SDe_dw_old
					ISDe_dw_new=ISDe_dw_old
					pvte_dw_new=pvte_dw_old
				END IF
				detSDe_up_new=detSDe_up_old
				detSDe_dw_new=detSDe_dw_old
			END IF
			
			!Jee
			IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='tre')) THEN
				IF (Jee_kind/='no_') THEN	
					u_ee_new=u_ee_old
				END IF
				Uee_new=Uee_old
			END IF
			
			!Jep
			IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='p__') .OR. (tipo=='tre')) THEN
				IF (Jep_kind/='no_') THEN
					u_ep_new=u_ep_old
				END IF
				Uep_new=Uep_old
			END IF
			
			!SDse
			IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
				IF (SDse_kind/='no_') THEN
					IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se1')) THEN
						SDse1_up_new=SDse1_up_old
						ISDse1_up_new=ISDse1_up_old
						pvtse1_up_new=pvtse1_up_old
						SDse1_dw_new=SDse1_dw_old
						ISDse1_dw_new=ISDse1_dw_old
						pvtse1_dw_new=pvtse1_dw_old
					ELSE IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se2')) THEN
						SDse2_up_new=SDse2_up_old
						ISDse2_up_new=ISDse2_up_old
						pvtse2_up_new=pvtse2_up_old
						SDse2_dw_new=SDse2_dw_old
						ISDse2_dw_new=ISDse2_dw_old
						pvtse2_dw_new=pvtse2_dw_old
					END IF
					
				END IF
				detSDse1_up_new=detSDse1_up_old
				detSDse1_dw_new=detSDse1_dw_old
				detSDse2_up_new=detSDse2_up_old
				detSDse2_dw_new=detSDse2_dw_old
			END IF
			
			!Jse
			IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
				IF ((Jse_kind/='no_') .AND. (Jse_kind/='bou').AND.(Jse_kind/='pot')) THEN
					IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='tre')) THEN
						u_se1_new=u_se1_old
						u_se2_new=u_se2_old
					ELSE IF (tipo=='se1') THEN
						u_se1_new=u_se1_old
					ELSE IF (tipo=='se2') THEN
						u_se2_new=u_se2_old
					END IF
				ELSE IF ((Jse_kind=='bou') .OR. (Jse_kind=='ppb')) THEN
					!IF ((tipo=='all') .OR. (tipo=='se_')) THEN
					!	b_se1_new=b_se1_old
					!	b_se2_new=b_se2_old
					!ELSE IF (tipo=='se1') THEN
					!	b_se1_new=b_se1_old
					!ELSE IF (tipo=='se2') THEN
					!	b_se2_new=b_se2_old
					!END IF
				ELSE IF (Jse_kind=='pot') THEN
					IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='tre')) THEN
						u_POT_se1_new=u_POT_se1_old
						u_POT_se2_new=u_POT_se2_old
					ELSE IF (tipo=='se1') THEN
						u_POT_se1_new=u_POT_se1_old
					ELSE IF (tipo=='se2') THEN
						u_POT_se2_new=u_POT_se2_old
					END IF
				END IF
				Use1_new=Use1_old
				Use2_new=Use2_old
				Bse1_new=Bse1_old
				Bse2_new=Bse2_old
			END IF
			
			!Kse
			IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se_') .OR. (tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
				IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
					IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se_') .OR. (tipo=='tre')) THEN
						u_ese1_new=u_ese1_old
						u_ese2_new=u_ese2_old
					ELSE IF (tipo=='se1') THEN
						u_ese1_new=u_ese1_old
					ELSE IF (tipo=='se2') THEN
						u_ese2_new=u_ese2_old
					END IF
					Uese1_new=Uese1_old
					Uese2_new=Uese2_old
				ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
					IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se_')) THEN
						GDse1_up_new=GDse1_up_old
						IGDse1_up_new=IGDse1_up_old
						pvtgdse1_up_new=pvtgdse1_up_old
						GDse1_dw_new=GDse1_dw_old
						IGDse1_dw_new=IGDse1_dw_old
						pvtgdse1_dw_new=pvtgdse1_dw_old
						GDse2_up_new=GDse2_up_old
						IGDse2_up_new=IGDse2_up_old
						pvtgdse2_up_new=pvtgdse2_up_old
						GDse2_dw_new=GDse2_dw_old
						IGDse2_dw_new=IGDse2_dw_old
						pvtgdse2_dw_new=pvtgdse2_dw_old
					ELSE IF (tipo=='se1') THEN
						GDse1_up_new=GDse1_up_old
						IGDse1_up_new=IGDse1_up_old
						pvtgdse1_up_new=pvtgdse1_up_old
						GDse1_dw_new=GDse1_dw_old
						IGDse1_dw_new=IGDse1_dw_old
						pvtgdse1_dw_new=pvtgdse1_dw_old
					ELSE IF (tipo=='se2') THEN
						GDse2_up_new=GDse2_up_old
						IGDse2_up_new=IGDse2_up_old
						pvtgdse2_up_new=pvtgdse2_up_old
						GDse2_dw_new=GDse2_dw_old
						IGDse2_dw_new=IGDse2_dw_old
						pvtgdse2_dw_new=pvtgdse2_dw_old
					END IF
					detGDse1_up_new=detGDse1_up_old
					detGDse1_dw_new=detGDse1_dw_old
					detGDse2_up_new=detGDse2_up_old
					detGDse2_dw_new=detGDse2_dw_old
				ELSE IF ((Kse_kind=='atm').OR.(Kse_kind=='atc')) THEN
					IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se_')) THEN
						u_ese1_new=u_ese1_old
						u_ese2_new=u_ese2_old
					ELSE IF (tipo=='se1') THEN
						u_ese1_new=u_ese1_old
					ELSE IF (tipo=='se2') THEN
						u_ese2_new=u_ese2_old
					END IF
					Uese1_new=Uese1_old
					Uese2_new=Uese2_old
				END IF
			END IF
			
			!Jsesp
			IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
				IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
					IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='tre')) THEN
						u_sesp1_new=u_sesp1_old
						u_sesp2_new=u_sesp2_old
					ELSE IF (tipo=='se1') THEN
						u_sesp1_new=u_sesp1_old
					ELSE IF (tipo=='se2') THEN
						u_sesp2_new=usesp2_old
					END IF
				ELSE IF (Jsesp_kind=='gsd') THEN    !OUTDATED
					IF ((tipo=='all') .OR. (tipo=='se_')) THEN
						GDsp1_up_new=GDsp1_up_old
						IGDsp1_up_new=IGDsp1_up_old
						pvtgdsp1_up_new=pvtgdsp1_up_old
						GDsp1_dw_new=GDsp1_dw_old
						IGDsp1_dw_new=IGDsp1_dw_old
						pvtgdsp1_dw_new=pvtgdsp1_dw_old
						GDsp2_up_new=GDsp2_up_old
						IGDsp2_up_new=IGDsp2_up_old
						pvtgdsp2_up_new=pvtgdsp2_up_old
						GDsp2_dw_new=GDsp2_dw_old
						IGDsp2_dw_new=IGDsp2_dw_old
						pvtgdsp2_dw_new=pvtgdsp2_dw_old
					ELSE IF (tipo=='se1') THEN
						GDsp1_up_new=GDsp1_up_old
						IGDsp1_up_new=IGDsp1_up_old
						pvtgdsp1_up_new=pvtgdsp1_up_old
						GDsp1_dw_new=GDsp1_dw_old
						IGDsp1_dw_new=IGDsp1_dw_old
						pvtgdsp1_dw_new=pvtgdsp1_dw_old
					ELSE IF (tipo=='se2') THEN
						GDsp2_up_new=GDsp2_up_old
						IGDsp2_up_new=IGDsp2_up_old
						pvtgdsp2_up_new=pvtgdsp2_up_old
						GDsp2_dw_new=GDsp2_dw_old
						IGDsp2_dw_new=IGDsp2_dw_old
						pvtgdsp2_dw_new=pvtgdsp2_dw_old
					END IF
					detGDsp1_up_new=detGDsp1_up_old
					detGDsp1_dw_new=detGDsp1_dw_old
					detGDsp2_up_new=detGDsp2_up_old
					detGDsp2_dw_new=detGDsp2_dw_old
				END IF
				Usesp1_new=Usesp1_old
				Usesp2_new=Usesp2_old
			END IF
		ELSE IF ((num>0) .AND. (num<=N_part)) THEN
			SELECT CASE (tipo)
			CASE ('e__')
				IF (SDe_kind/='no_') THEN
					IF (num<=H_N_part) THEN
						SDe_up_new(num,1:H_N_part)=SDe_up_old(num,1:H_N_part)
						detSDe_up_new=detSDe_up_old
					ELSE
						SDe_dw_new(num-H_N_part,1:H_N_part)=SDe_dw_old(num-H_N_part,1:H_N_part)
						detSDe_dw_new=detSDe_dw_old
					END IF
				END IF
				IF (Jee_kind/='no_') THEN
					u_ee_new(num,1:num-1)=u_ee_old(num,1:num-1)
					u_ee_new(num+1:N_part,num)=u_ee_old(num+1:N_part,num)
					Uee_new=Uee_old
				END IF
				IF (Jep_kind/='no_') THEN
					u_ep_new(num,1:N_part)=u_ep_old(num,1:N_part)
					Uep_new=Uep_old
				END IF
				IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
					u_ese1_new(num)=u_ese1_old(num)
					Uese1_new=Uese1_old
					u_ese2_new(num)=u_ese2_old(num)
					Uese2_new=Uese2_old
				ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
					IF (num<=H_N_part) THEN
						GDse1_up_new(num,1:H_N_part)=GDse1_up_old(num,1:H_N_part)
						detGDse1_up_new=detGDse1_up_old
						GDse2_up_new(num,1:H_N_part)=GDse2_up_old(num,1:H_N_part)
						detGDse2_up_new=detGDse2_up_old
					ELSE
						GDse1_dw_new(num-H_N_part,1:H_N_part)=GDse1_dw_old(num-H_N_part,1:H_N_part)
						detGDse1_dw_new=detGDse1_dw_old
						GDse2_dw_new(num-H_N_part,1:H_N_part)=GDse2_dw_old(num-H_N_part,1:H_N_part)
						detGDse2_dw_new=detGDse2_dw_old
					END IF
				ELSE IF ((Kse_kind=='atm').OR.(Kse_kind=='atc')) THEN
					u_ese1_new(num)=u_ese1_old(num)
					Uese1_new=Uese1_old
					u_ese2_new(num)=u_ese2_old(num)
					Uese2_new=Uese2_old
				END IF
			CASE ('se1')
				IF (SDse_kind/='no_') THEN
					IF (num<=H_N_part) THEN
						SDse1_up_new(num,1:H_N_part)=SDse1_up_old(num,1:H_N_part)
						detSDse1_up_new=detSDse1_up_old
					ELSE
						SDse1_dw_new(num-H_N_part,1:H_N_part)=SDse1_dw_old(num-H_N_part,1:H_N_part)
						detSDse1_dw_new=detSDse1_dw_old
					END IF
				END IF
				IF ((Jse_kind/='no_') .AND. (Jse_kind/='bou') .AND. (Jse_kind/='pot')) THEN
					u_se1_new(num,1:N_part)=u_se1_old(num,1:N_part)
					u_se1_new(1:N_part,num)=u_se1_old(1:N_part,num)
					Use1_new=Use1_old
				END IF
				IF ((Jse_kind=='bou').OR.(Jse_kind=='ppb')) THEN
					!b_se1_new(num,1:N_part)=b_se1_old(num,1:N_part)
					!b_se1_new(1:N_part,num)=b_se1_old(1:N_part,num)
					!Bse1_new=Bse1_old
				END IF
				IF (Jse_kind=='pot') THEN
					u_POT_se1_new(index_mol_num)=u_POT_se1_old(index_mol_num)
					u_POT_se1_new(index_mol_num)=u_POT_se1_old(index_mol_num)
					Use1_new=Use1_old
				END IF
				IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
					u_ese1_new(num)=u_ese1_old(num)
					Uese1_new=Uese1_old
				ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
					IF (num<=H_N_part) THEN
						GDse1_up_new(1:H_N_part,num)=GDse1_up_old(1:H_N_part,num)
						detGDse1_up_new=detGDse1_up_old
					ELSE
						GDse1_dw_new(1:H_N_part,num-H_N_part)=GDse1_dw_old(1:H_N_part,num-H_N_part)
						detGDse1_dw_new=detGDse1_dw_old
					END IF
				ELSE IF ((Kse_kind=='atm').OR.(Kse_kind=='atc')) THEN
					u_ese1_new(num)=u_ese1_old(num)
					Uese1_new=Uese1_old
				END IF
				IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
					u_sesp1_new(num,1:N_part)=u_sesp1_old(num,1:N_part)
					Usesp1_new=Usesp1_old
				ELSE IF (Jsesp_kind=='gsd') THEN
					IF (num<=H_N_part) THEN
						GDsp1_up_new(num,1:H_N_part)=GDsp1_up_old(num,1:H_N_part)
						detGDsp1_up_new=detGDsp1_up_old
					ELSE
						GDsp1_dw_new(num-H_N_part,1:H_N_part)=GDsp1_dw_old(num-H_N_part,1:H_N_part)
						detGDsp1_dw_new=detGDsp1_dw_old
					END IF
				END IF
			CASE ('se2')
				IF (SDse_kind/='no_') THEN
					IF (num<=H_N_part) THEN
						SDse2_up_new(num,1:H_N_part)=SDse2_up_old(num,1:H_N_part)
						detSDse2_up_new=detSDse2_up_old
					ELSE
						SDse2_dw_new(num-H_N_part,1:H_N_part)=SDse2_dw_old(num-H_N_part,1:H_N_part)
						detSDse2_dw_new=detSDse2_dw_old
					END IF
				END IF
				IF ((Jse_kind/='no_') .AND. (Jse_kind/='bou') .AND. (Jse_kind/='pot')) THEN
					u_se2_new(num,1:N_part)=u_se2_old(num,1:N_part)
					u_se2_new(1:N_part,num)=u_se2_old(1:N_part,num)
					Use2_new=Use2_old
				END IF
				IF ((Jse_kind=='bou').OR.(Jse_kind=='ppb')) THEN
					!b_se2_new(num,1:N_part)=b_se2_old(num,1:N_part)
					!b_se2_new(1:N_part,num)=b_se2_old(1:N_part,num)
					!Bse2_new=Bse2_old
				END IF
				IF (Jse_kind=='pot') THEN
					u_POT_se2_new(index_mol_num)=u_POT_se2_old(index_mol_num)
					u_POT_se2_new(index_mol_num)=u_POT_se2_old(index_mol_num)
					Use2_new=Use2_old
				END IF
				IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
					u_ese2_new(num)=u_ese2_old(num)
					Uese2_new=Uese2_old
				ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
					IF (num<=H_N_part) THEN
						GDse2_up_new(1:H_N_part,num)=GDse2_up_old(1:H_N_part,num)
						detGDse2_up_new=detGDse2_up_old
					ELSE
						GDse2_dw_new(1:H_N_part,num-H_N_part)=GDse2_dw_old(1:H_N_part,num-H_N_part)
						detGDse2_dw_new=detGDse2_dw_old
					END IF
				ELSE IF ((Kse_kind=='atm').OR.(Kse_kind=='atc')) THEN
					u_ese2_new(num)=u_ese2_old(num)
					Uese2_new=Uese2_old
				END IF
				IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
					u_sesp2_new(num,1:N_part)=u_sesp2_old(num,1:N_part)
					Usesp2_new=Usesp2_old
				ELSE IF (Jsesp_kind=='gsd') THEN
					IF (num<=H_N_part) THEN
						GDsp2_up_new(num,1:H_N_part)=GDsp2_up_old(num,1:H_N_part)
						detGDsp2_up_new=detGDsp2_up_old
					ELSE
						GDsp2_dw_new(num-H_N_part,1:H_N_part)=GDsp2_dw_old(num-H_N_part,1:H_N_part)
						detGDsp2_dw_new=detGDsp2_dw_old
					END IF
				END IF
			CASE ('tre')
				IF (SDe_kind/='no_') THEN
					IF (num<=H_N_part) THEN
						SDe_up_new(num,1:H_N_part)=SDe_up_old(num,1:H_N_part)
						detSDe_up_new=detSDe_up_old
					ELSE
						SDe_dw_new(num-H_N_part,1:H_N_part)=SDe_dw_old(num-H_N_part,1:H_N_part)
						detSDe_dw_new=detSDe_dw_old
					END IF
				END IF
				IF (Jee_kind/='no_') THEN
					u_ee_new(num,1:num-1)=u_ee_old(num,1:num-1)
					u_ee_new(num+1:N_part,num)=u_ee_old(num+1:N_part,num)
					Uee_new=Uee_old
				END IF
				IF (Jep_kind/='no_') THEN
					u_ep_new(num,1:N_part)=u_ep_old(num,1:N_part)
					Uep_new=Uep_old
				END IF
				IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
					Uese1_new=Uese1_old
					Uese2_new=Uese2_old
				ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
					!TO DO
				ELSE IF ((Kse_kind=='atm').OR.(Kse_kind=='atc')) THEN
					!TO DO
				END IF
				IF (SDse_kind/='no_') THEN
					IF (num<=H_N_part) THEN
						SDse1_up_new(num,1:H_N_part)=SDse1_up_old(num,1:H_N_part)
						detSDse1_up_new=detSDse1_up_old
					ELSE
						SDse1_dw_new(num-H_N_part,1:H_N_part)=SDse1_dw_old(num-H_N_part,1:H_N_part)
						detSDse1_dw_new=detSDse1_dw_old
					END IF
				END IF
				IF ((Jse_kind/='no_') .AND. (Jse_kind/='bou') .AND. (Jse_kind/='pot')) THEN
					u_se1_new(num,1:N_part)=u_se1_old(num,1:N_part)
					u_se1_new(1:N_part,num)=u_se1_old(1:N_part,num)
					Use1_new=Use1_old
				END IF
				IF ((Jse_kind=='bou').OR.(Jse_kind=='ppb')) THEN
					!b_se1_new(num,1:N_part)=b_se1_old(num,1:N_part)
					!b_se1_new(1:N_part,num)=b_se1_old(1:N_part,num)
					!Bse1_new=Bse1_old
				END IF
				IF (Jse_kind=='pot') THEN
					u_POT_se1_new(index_mol_num)=u_POT_se1_old(index_mol_num)
					u_POT_se1_new(index_mol_num)=u_POT_se1_old(index_mol_num)
					Use1_new=Use1_old
				END IF
				IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
					u_sesp1_new(num,1:N_part)=u_sesp1_old(num,1:N_part)
					Usesp1_new=Usesp1_old
				ELSE IF (Jsesp_kind=='gsd') THEN
					IF (num<=H_N_part) THEN
						GDsp1_up_new(num,1:H_N_part)=GDsp1_up_old(num,1:H_N_part)
						detGDsp1_up_new=detGDsp1_up_old
					ELSE
						GDsp1_dw_new(num-H_N_part,1:H_N_part)=GDsp1_dw_old(num-H_N_part,1:H_N_part)
						detGDsp1_dw_new=detGDsp1_dw_old
					END IF
				END IF
				IF (SDse_kind/='no_') THEN
					IF (num<=H_N_part) THEN
						SDse2_up_new(num,1:H_N_part)=SDse2_up_old(num,1:H_N_part)
						detSDse2_up_new=detSDse2_up_old
					ELSE
						SDse2_dw_new(num-H_N_part,1:H_N_part)=SDse2_dw_old(num-H_N_part,1:H_N_part)
						detSDse2_dw_new=detSDse2_dw_old
					END IF
				END IF
				IF ((Jse_kind/='no_') .AND. (Jse_kind/='bou') .AND. (Jse_kind/='pot')) THEN
					u_se2_new(num,1:N_part)=u_se2_old(num,1:N_part)
					u_se2_new(1:N_part,num)=u_se2_old(1:N_part,num)
					Use2_new=Use2_old
				END IF
				IF ((Jse_kind=='bou').OR.(Jse_kind=='ppb')) THEN
					!b_se2_new(num,1:N_part)=b_se2_old(num,1:N_part)
					!b_se2_new(1:N_part,num)=b_se2_old(1:N_part,num)
					!Bse2_new=Bse2_old
				END IF
				IF (Jse_kind=='pot') THEN
					u_POT_se2_new(index_mol_num)=u_POT_se2_old(index_mol_num)
					u_POT_se2_new(index_mol_num)=u_POT_se2_old(index_mol_num)
					Use2_new=Use2_old
				END IF
				IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
					u_sesp2_new(num,1:N_part)=u_sesp2_old(num,1:N_part)
					Usesp2_new=Usesp2_old
				ELSE IF (Jsesp_kind=='gsd') THEN
					IF (num<=H_N_part) THEN
						GDsp2_up_new(num,1:H_N_part)=GDsp2_up_old(num,1:H_N_part)
						detGDsp2_up_new=detGDsp2_up_old
					ELSE
						GDsp2_dw_new(num-H_N_part,1:H_N_part)=GDsp2_dw_old(num-H_N_part,1:H_N_part)
						detGDsp2_dw_new=detGDsp2_dw_old
					END IF
				END IF
			CASE ('p__')
				IF (Jep_kind/='no_') THEN
					u_ep_new(1:N_part,num)=u_ep_old(1:N_part,num)
					Uep_new=Uep_old
				END IF
			END SELECT
		END IF
				
	END SUBROUTINE roll_back
!-----------------------------------------------------------------------

	SUBROUTINE aggiorna_funzione_onda(tipo,num)
		USE generic_tools
		USE dati_simulazione_mc
		USE dati_fisici
		USE funzione_onda
		USE walkers
		IMPLICIT NONE
		CHARACTER(LEN=3) :: tipo
		INTEGER, INTENT(IN) :: num
		INTEGER :: i, info, M_pvt(1:H_N_part), perm
		REAL (KIND=8) :: M(1:H_N_part,1:H_N_part), work(1:3*H_N_part), detM
		
		IF (.NOT. iniz_calcola_accettazione) STOP 'Prima di aggiornare la funzione d onda devi inizializzarla &
		  [ module_calcola_accettazione.f90 > aggiorna_funzione_onda ]'
				
		IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='tre')) THEN
			IF (SDe_kind/='no_') THEN
				IF (num==-1) THEN
					SDe_up_old=SDe_up_new
					ISDe_up_old=ISDe_up_new
					pvte_up_old=pvte_up_new
					SDe_dw_old=SDe_dw_new
					ISDe_dw_old=ISDe_dw_new
					pvte_dw_old=pvte_dw_new
				ELSE IF ((num>0) .AND. (num<=N_part)) THEN
					IF (num<=H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num,ISDe_up_old,detSDe_up_old,SDe_up_new,detSDe_up_new,ISDe_up_new)
						SDe_up_old(num,1:H_N_part)=SDe_up_new(num,1:H_N_part)
						ISDe_up_old=ISDe_up_new
					ELSE IF (num>H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num-H_N_part,ISDe_dw_old,detSDe_dw_old,SDe_dw_new,detSDe_dw_new,ISDe_dw_new)
						SDe_dw_old(num-H_N_part,1:H_N_part)=SDe_dw_new(num-H_N_part,1:H_N_part)
						ISDe_dw_old=ISDe_dw_new
					END IF
				END IF
			END IF
			detSDe_up_old=detSDe_up_new
			detSDe_dw_old=detSDe_dw_new
			
			IF (Jee_kind/='no_') THEN	
				IF (num==-1) THEN
					u_ee_old=u_ee_new
				ELSE IF ((num>0) .AND. (num<=N_part)) THEN
					u_ee_old(num,1:num-1)=u_ee_new(num,1:num-1)
					u_ee_old(num+1:N_part,num)=u_ee_new(num+1:N_part,num)
				END IF
			END IF
			Uee_old=Uee_new
		END IF
		
		IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='p__') .OR. (tipo=='tre')) THEN
			IF (Jep_kind/='no_') THEN
				IF (num==-1) THEN
					u_ep_old=u_ep_new
				ELSE IF ((num>0) .AND. (num<=N_part)) THEN
					IF ((tipo=='e__') .OR. (tipo=='tre')) THEN
						u_ep_old(num,1:N_part)=u_ep_new(num,1:N_part)
					ELSE IF (tipo=='p__') THEN
						u_ep_old(1:N_part,num)=u_ep_new(1:N_part,num)
					ELSE
						STOP 'Tipo non accettabile --code1-- &
						  [ module_calcola_accettazione.f90 > aggiorna_funzione_onda ]'
					END IF
				END IF
			END IF
			Uep_old=Uep_new
		END IF
		
		!SDse
		IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se1') .OR. (tipo=='tre')) THEN
			IF ( SDse_kind=='gem' ) THEN
				IF (num==-1) THEN
					SDse1_up_old=SDse1_up_new
					ISDse1_up_old=ISDse1_up_new
					pvtse1_up_old=pvtse1_up_new
				ELSE IF ((num>0) .AND. (num<=N_part)) THEN
					IF (num<=H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num,ISDse1_up_old,detSDse1_up_old, &
						  SDse1_up_new,detSDse1_up_new,ISDse1_up_new)
						SDse1_up_old(num,1:H_N_part)=SDse1_up_new(num,1:H_N_part)
						ISDse1_up_old=ISDse1_up_new
					ELSE IF (num>H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_col_1ppt(H_N_part,num-H_N_part,ISDse1_up_old,detSDse1_up_old, &
						  SDse1_up_new,detSDse1_up_new,ISDse1_up_new)
						SDse1_up_old(1:H_N_part,num-H_N_part)=SDse1_up_new(1:H_N_part,num-H_N_part)
						ISDse1_up_old=ISDse1_up_new
					END IF
				END IF
			ELSE IF (SDse_kind/='no_') THEN
				IF (num==-1) THEN
					SDse1_up_old=SDse1_up_new
					ISDse1_up_old=ISDse1_up_new
					pvtse1_up_old=pvtse1_up_new
					SDse1_dw_old=SDse1_dw_new
					ISDse1_dw_old=ISDse1_dw_new
					pvtse1_dw_old=pvtse1_dw_new
				ELSE IF ((num>0) .AND. (num<=N_part)) THEN
					IF (num<=H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num,ISDse1_up_old,detSDse1_up_old, &
						  SDse1_up_new,detSDse1_up_new,ISDse1_up_new)
						SDse1_up_old(num,1:H_N_part)=SDse1_up_new(num,1:H_N_part)
						ISDse1_up_old=ISDse1_up_new
					ELSE IF (num>H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num-H_N_part,ISDse1_dw_old,detSDse1_dw_old, &
						  SDse1_dw_new,detSDse1_dw_new,ISDse1_dw_new)
						SDse1_dw_old(num-H_N_part,1:H_N_part)=SDse1_dw_new(num-H_N_part,1:H_N_part)
						ISDse1_dw_old=ISDse1_dw_new
					END IF
				END IF
			END IF
			detSDse1_up_old=detSDse1_up_new
			detSDse1_dw_old=detSDse1_dw_new
		END IF
		IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
			IF (SDse_kind=='gem') THEN
				IF (num==-1) THEN
					SDse2_up_old=SDse2_up_new
					ISDse2_up_old=ISDse2_up_new
					pvtse2_up_old=pvtse2_up_new
				ELSE IF ((num>0) .AND. (num<=N_part)) THEN
					IF (num<=H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num,ISDse2_up_old,detSDse2_up_old, &
						  SDse2_up_new,detSDse2_up_new,ISDse2_up_new)
						SDse2_up_old(num,1:H_N_part)=SDse2_up_new(num,1:H_N_part)
						ISDse2_up_old=ISDse2_up_new
					ELSE IF (num>H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_col_1ppt(H_N_part,num-H_N_part,ISDse2_up_old,detSDse2_up_old, &
						  SDse2_up_new,detSDse2_up_new,ISDse2_up_new)
						SDse2_up_old(1:H_N_part,num-H_N_part)=SDse2_up_new(1:H_N_part,num-H_N_part)
						ISDse2_up_old=ISDse2_up_new
					END IF
				END IF
			ELSE IF (SDse_kind/='no_') THEN
				IF (num==-1) THEN
					SDse2_up_old=SDse2_up_new
					ISDse2_up_old=ISDse2_up_new
					pvtse2_up_old=pvtse2_up_new
					SDse2_dw_old=SDse2_dw_new
					ISDse2_dw_old=ISDse2_dw_new
					pvtse2_dw_old=pvtse2_dw_new
				ELSE IF ((num>0) .AND. (num<=N_part)) THEN
					IF (num<=H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num,ISDse2_up_old,detSDse2_up_old, &
						  SDse2_up_new,detSDse2_up_new,ISDse2_up_new)
						SDse2_up_old(num,1:H_N_part)=SDse2_up_new(num,1:H_N_part)
						ISDse2_up_old=ISDse2_up_new
					ELSE IF (num>H_N_part) THEN
						CALL aggiorna_matrice_inversa_C_1ppt(H_N_part,num-H_N_part,ISDse2_dw_old,detSDse2_dw_old, &
						  SDse2_dw_new,detSDse2_dw_new,ISDse2_dw_new)
						SDse2_dw_old(num-H_N_part,1:H_N_part)=SDse2_dw_new(num-H_N_part,1:H_N_part)
						ISDse2_dw_old=ISDse2_dw_new
					END IF
				END IF
			END IF
			detSDse2_up_old=detSDse2_up_new
			detSDse2_dw_old=detSDse2_dw_new
		END IF
		
		!Jsese
		IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
			IF ((Jse_kind/='no_') .AND. (Jse_kind/='bou') .AND. (Jse_kind/='pot')) THEN
				IF ((tipo=='all') .OR. (tipo=='se_')) THEN
					IF (num==-1) THEN
						u_se1_old=u_se1_new
						u_se2_old=u_se2_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_se1_old(num,1:N_part)=u_se1_new(num,1:N_part)
						u_se1_old(1:N_part,num)=u_se1_new(1:N_part,num)
						u_se2_old(num,1:N_part)=u_se2_new(num,1:N_part)
						u_se2_old(1:N_part,num)=u_se2_new(1:N_part,num)
					END IF
				ELSE IF (tipo=='se1') THEN
					IF (num==-1) THEN
						u_se1_old=u_se1_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_se1_old(num,1:N_part)=u_se1_new(num,1:N_part)
						u_se1_old(1:N_part,num)=u_se1_new(1:N_part,num)
					END IF
				ELSE IF (tipo=='se2') THEN
					IF (num==-1) THEN
						u_se2_old=u_se2_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_se2_old(num,1:N_part)=u_se2_new(num,1:N_part)
						u_se2_old(1:N_part,num)=u_se2_new(1:N_part,num)
					END IF
				ELSE IF (tipo=='tre') THEN
					IF ((num>0) .AND. (num<=N_part)) THEN
						u_se1_old(num,1:N_part)=u_se1_new(num,1:N_part)
						u_se1_old(1:N_part,num)=u_se1_new(1:N_part,num)
						u_se2_old(num,1:N_part)=u_se2_new(num,1:N_part)
						u_se2_old(1:N_part,num)=u_se2_new(1:N_part,num)
					END IF
				END IF
			END IF
			IF (Jse_kind=='pot') THEN
				IF ((tipo=='all') .OR. (tipo=='se_')) THEN
					IF (index_mol_num==-1) THEN
						u_POT_se1_old=u_POT_se1_new
						u_POT_se2_old=u_POT_se2_new
					ELSE IF ((index_mol_num>0) .AND. (index_mol_num<=H_N_part)) THEN
						u_POT_se1_old(index_mol_num)=u_POT_se1_new(index_mol_num)
						u_POT_se1_old(index_mol_num)=u_POT_se1_new(index_mol_num)
						u_POT_se2_old(index_mol_num)=u_POT_se2_new(index_mol_num)
						u_POT_se2_old(index_mol_num)=u_POT_se2_new(index_mol_num)
					END IF
				ELSE IF (tipo=='se1') THEN
					IF (index_mol_num==-1) THEN
						u_POT_se1_old=u_POT_se1_new
					ELSE IF ((index_mol_num>0) .AND. (index_mol_num<=H_N_part)) THEN
						u_POT_se1_old(index_mol_num)=u_POT_se1_new(index_mol_num)
						u_POT_se1_old(index_mol_num)=u_POT_se1_new(index_mol_num)
					END IF
				ELSE IF (tipo=='se2') THEN
					IF (index_mol_num==-1) THEN
						u_POT_se2_old=u_POT_se2_new
					ELSE IF ((index_mol_num>0) .AND. (index_mol_num<=H_N_part)) THEN
						u_POT_se2_old(index_mol_num)=u_POT_se2_new(index_mol_num)
						u_POT_se2_old(index_mol_num)=u_POT_se2_new(index_mol_num)
					END IF
				ELSE IF (tipo=='tre') THEN
					IF (index_mol_num==-1) THEN
						u_POT_se1_old=u_POT_se1_new
						u_POT_se2_old=u_POT_se2_new
					ELSE IF ((index_mol_num>0) .AND. (index_mol_num<=H_N_part)) THEN
						u_POT_se1_old(index_mol_num)=u_POT_se1_new(index_mol_num)
						u_POT_se1_old(index_mol_num)=u_POT_se1_new(index_mol_num)
						u_POT_se2_old(index_mol_num)=u_POT_se2_new(index_mol_num)
						u_POT_se2_old(index_mol_num)=u_POT_se2_new(index_mol_num)
					END IF
				END IF
			END IF
			Use1_old=Use1_new
			Use2_old=Use2_new
			IF ((Jse_kind=='bou').OR.(Jse_kind=='ppb')) THEN
				!IF ((tipo=='all') .OR. (tipo=='se_')) THEN
				!	IF (num==-1) THEN
				!		b_se1_old=b_se1_new
				!		b_se2_old=b_se2_new
				!	ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				!		b_se1_old(num,1:N_part)=b_se1_new(num,1:N_part)
				!		b_se1_old(1:N_part,num)=b_se1_new(1:N_part,num)
				!		b_se2_old(num,1:N_part)=b_se2_new(num,1:N_part)
				!		b_se2_old(1:N_part,num)=b_se2_new(1:N_part,num)
				!	END IF
				!ELSE IF (tipo=='se1') THEN
				!	IF (num==-1) THEN
				!		b_se1_old=b_se1_new
				!	ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				!		b_se1_old(num,1:N_part)=b_se1_new(num,1:N_part)
				!		b_se1_old(1:N_part,num)=b_se1_new(1:N_part,num)
				!	END IF
				!ELSE IF (tipo=='se2') THEN
				!	IF (num==-1) THEN
				!		b_se2_old=b_se2_new
				!	ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				!		b_se2_old(num,1:N_part)=b_se2_new(num,1:N_part)
				!		b_se2_old(1:N_part,num)=b_se2_new(1:N_part,num)
				!	END IF
				!END IF
			END IF
			Bse1_old=Bse1_new
			Bse2_old=Bse2_new
		END IF
		
		!Kse
		IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se_') .OR. (tipo=='se1') .OR. (tipo=='se2') .OR. (tipo=='tre')) THEN
			IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
				IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se_')) THEN
					IF (num==-1) THEN
						u_ese1_old=u_ese1_new
						u_ese2_old=u_ese2_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_ese1_old(num)=u_ese1_new(num)
						u_ese2_old(num)=u_ese2_new(num)	
					END IF
				ELSE IF (tipo=='se1') THEN
					!PRINT * , 'aggiorno se1'
					IF (num==-1) THEN
						u_ese1_old=u_ese1_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_ese1_old(num)=u_ese1_new(num)
					END IF
				ELSE IF (tipo=='se2') THEN
					!PRINT * , 'aggiorno se2'
					IF (num==-1) THEN
						u_ese2_old=u_ese2_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_ese2_old(num)=u_ese2_new(num)
					END IF
				ELSE IF (tipo=='tre') THEN
					!PRINT * , 'aggiorno tre'
					IF ((num>0) .AND. (num<=N_part)) THEN
						u_ese1_old(num)=u_ese1_new(num)
						u_ese2_old(num)=u_ese2_new(num)
					END IF
				END IF
				Uese1_old=Uese1_new
				Uese2_old=Uese2_new
			ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
				IF (num==-1) THEN
					IF ((tipo=='all').OR.(tipo=='e__')) THEN
						GDse1_up_old=GDse1_up_new
						IGDse1_up_old=IGDse1_up_new
						pvtgdse1_up_old=pvtgdse1_up_new
						GDse1_dw_old=GDse1_dw_new
						IGDse1_dw_old=IGDse1_dw_new
						pvtgdse1_dw_old=pvtgdse1_dw_new
						GDse2_up_old=GDse2_up_new
						IGDse2_up_old=IGDse2_up_new
						pvtgdse2_up_old=pvtgdse2_up_new
						GDse2_dw_old=GDse2_dw_new
						IGDse2_dw_old=IGDse2_dw_new
						pvtgdse2_dw_old=pvtgdse2_dw_new
					ELSE IF (tipo=='se1') THEN
						GDse1_up_old=GDse1_up_new
						IGDse1_up_old=IGDse1_up_new
						pvtgdse1_up_old=pvtgdse1_up_new
						GDse1_dw_old=GDse1_dw_new
						IGDse1_dw_old=IGDse1_dw_new
						pvtgdse1_dw_old=pvtgdse1_dw_new
					ELSE IF (tipo=='se2') THEN
						GDse2_up_old=GDse2_up_new
						IGDse2_up_old=IGDse2_up_new
						pvtgdse2_up_old=pvtgdse2_up_new
						GDse2_dw_old=GDse2_dw_new
						IGDse2_dw_old=IGDse2_dw_new
						pvtgdse2_dw_old=pvtgdse2_dw_new
					END IF
				ELSE IF (num<=H_N_part) THEN
					IF (tipo=='e__') THEN
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num,IGDse1_up_old,detGDse1_up_old,GDse1_up_new,detGDse1_up_new,IGDse1_up_new)
						GDse1_up_old(num,1:H_N_part)=GDse1_up_new(num,1:H_N_part)
						IGDse1_up_old=IGDse1_up_new
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num,IGDse2_up_old,detGDse2_up_old,GDse2_up_new,detGDse2_up_new,IGDse2_up_new)
						GDse2_up_old(num,1:H_N_part)=GDse2_up_new(num,1:H_N_part)
						IGDse2_up_old=IGDse2_up_new
					ELSE IF (tipo=='se1') THEN
						CALL aggiorna_matrice_inversa_R_col_1ppt(H_N_part,num,IGDse1_up_old,detGDse1_up_old,GDse1_up_new, &
						  detGDse1_up_new,IGDse1_up_new)
						GDse1_up_old(1:H_N_part,num)=GDse1_up_new(1:H_N_part,num)
						IGDse1_up_old=IGDse1_up_new
					ELSE IF (tipo=='se2') THEN	
						CALL aggiorna_matrice_inversa_R_col_1ppt(H_N_part,num,IGDse2_up_old,detGDse2_up_old,GDse2_up_new, &
						  detGDse2_up_new,IGDse2_up_new)
						GDse2_up_old(1:H_N_part,num)=GDse2_up_new(1:H_N_part,num)
						IGDse2_up_old=IGDse2_up_new
					ELSE IF (tipo=='tre') THEN
						!TODO
					END IF
				ELSE IF (num>H_N_part) THEN
					IF (tipo=='e__') THEN
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num-H_N_part,IGDse1_dw_old,detGDse1_dw_old, &
						  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new)
						GDse1_dw_old(num-H_N_part,1:H_N_part)=GDse1_dw_new(num-H_N_part,1:H_N_part)
						IGDse1_dw_old=IGDse1_dw_new		
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num-H_N_part,IGDse2_dw_old,detGDse2_dw_old, &
						  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new)
						GDse2_dw_old(num-H_N_part,1:H_N_part)=GDse2_dw_new(num-H_N_part,1:H_N_part)
						IGDse2_dw_old=IGDse2_dw_new
					ELSE IF (tipo=='se1') THEN
						CALL aggiorna_matrice_inversa_R_col_1ppt(H_N_part,num-H_N_part,IGDse1_dw_old,detGDse1_dw_old, &
						  GDse1_dw_new,detGDse1_dw_new,IGDse1_dw_new)
						GDse1_dw_old(1:H_N_part,num-H_N_part)=GDse1_dw_new(1:H_N_part,num-H_N_part)
						IGDse1_dw_old=IGDse1_dw_new
					ELSE IF (tipo=='se2') THEN	
						CALL aggiorna_matrice_inversa_R_col_1ppt(H_N_part,num-H_N_part,IGDse2_dw_old,detGDse2_dw_old, &
						  GDse2_dw_new,detGDse2_dw_new,IGDse2_dw_new)
						GDse2_dw_old(1:H_N_part,num-H_N_part)=GDse2_dw_new(1:H_N_part,num-H_N_part)
						IGDse2_dw_old=IGDse2_dw_new
					ELSE IF (tipo=='tre') THEN
						!TODO
					END IF
				END IF
				detGDse1_up_old=detGDse1_up_new
				detGDse1_dw_old=detGDse1_dw_new
				detGDse2_up_old=detGDse2_up_new
				detGDse2_dw_old=detGDse2_dw_new
			ELSE IF ((Kse_kind=='atm').OR.(Kse_kind=='atc')) THEN
				IF ((tipo=='all') .OR. (tipo=='e__') .OR. (tipo=='se_')) THEN
					IF (num==-1) THEN
						u_ese1_old=u_ese1_new
						u_ese2_old=u_ese2_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_ese1_old(num)=u_ese1_new(num)
						u_ese2_old(num)=u_ese2_new(num)
					END IF
				ELSE IF (tipo=='se1') THEN
					IF (num==-1) THEN
						u_ese1_old=u_ese1_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_ese1_old(num)=u_ese1_new(num)
					END IF
				ELSE IF (tipo=='se2') THEN
					IF (num==-1) THEN
						u_ese2_old=u_ese2_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_ese2_old(num)=u_ese2_new(num)
					END IF
				ELSE IF (tipo=='tre') THEN
					IF ((num>0) .AND. (num<=N_part)) THEN
						u_ese1_old(num)=u_ese1_new(num)
						u_ese2_old(num)=u_ese2_new(num)
					END IF
				END IF
				Uese1_old=Uese1_new
				Uese2_old=Uese2_new
			END IF
		END IF
		
		!Jsesp
		IF ((tipo=='all') .OR. (tipo=='se_') .OR. (tipo=='se1') .OR. (tipo=='se2')) THEN
			IF ((Jsesp_kind/='no_') .AND. (Jsesp_kind/='gsd')) THEN
				IF ((tipo=='all') .OR. (tipo=='se_')) THEN
					IF (num==-1) THEN
						u_sesp1_old=u_sesp1_new
						u_sesp2_old=u_sesp2_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_sesp1_old(num,1:N_part)=u_sesp1_new(num,1:N_part)
						u_sesp2_old(num,1:N_part)=u_sesp2_new(num,1:N_part)
					END IF
				ELSE IF (tipo=='se1') THEN
					IF (num==-1) THEN
						u_sesp1_old=u_sesp1_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_sesp1_old(num,1:N_part)=u_sesp1_new(num,1:N_part)
					END IF
				ELSE IF (tipo=='se2') THEN
					IF (num==-1) THEN
						u_sesp2_old=u_sesp2_new
					ELSE IF ((num>0) .AND. (num<=N_part)) THEN
						u_sesp2_old(num,1:N_part)=u_sesp2_new(num,1:N_part)
					END IF
				ELSE IF (tipo=='tre') THEN
					IF ((num>0) .AND. (num<=N_part)) THEN
						u_sesp1_old(num,1:N_part)=u_sesp1_new(num,1:N_part)
						u_sesp2_old(num,1:N_part)=u_sesp2_new(num,1:N_part)
					END IF
				END IF
			ELSE IF (Jsesp_kind=='gsd') THEN !OUTDATED
				IF ((tipo=='all') .OR. (tipo=='se_')) THEN
					IF (num==-1) THEN
						GDsp1_up_old=GDsp1_up_new
						IGDsp1_up_old=IGDsp1_up_new
						pvtgdsp1_up_old=pvtgdsp1_up_new
						GDsp1_dw_old=GDsp1_dw_new
						IGDsp1_dw_old=IGDsp1_dw_new
						pvtgdsp1_dw_old=pvtgdsp1_dw_new
						GDsp2_up_old=GDsp2_up_new
						IGDsp2_up_old=IGDsp2_up_new
						pvtgdsp2_up_old=pvtgdsp2_up_new
						GDsp2_dw_old=GDsp2_dw_new
						IGDsp2_dw_old=IGDsp2_dw_new
						pvtgdsp2_dw_old=pvtgdsp2_dw_new
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num,IGDsp1_up_old,detGDsp1_up_old, &
						  GDsp1_up_new,detGDsp1_up_new,IGDsp1_up_new)
						GDsp1_up_old(num,1:H_N_part)=GDsp1_up_new(num,1:H_N_part)
						IGDsp1_up_old=IGDsp1_up_new
						pvtgdsp1_up_old=pvtgdsp1_up_new
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num,IGDsp2_up_old,detGDsp2_up_old, &
						  GDsp2_up_new,detGDsp2_up_new,IGDsp2_up_new)
						GDsp2_up_old(num,1:H_N_part)=GDsp2_up_new(num,1:H_N_part)
						IGDsp2_up_old=IGDsp2_up_new
						pvtgdsp2_up_old=pvtgdsp2_up_new
					ELSE IF ((num>H_N_part).AND.(num<=N_part)) THEN
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num-H_N_part,IGDsp1_dw_old,detGDsp1_dw_old, &
						  GDsp1_dw_new,detGDsp1_dw_new,IGDsp1_dw_new)
						GDsp1_dw_old(num,1:H_N_part)=GDsp1_dw_new(num,1:H_N_part)
						IGDsp1_dw_old=IGDsp1_dw_new
						pvtgdsp1_dw_old=pvtgdsp1_dw_new
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num-H_N_part,IGDsp2_dw_old,detGDsp2_dw_old, &
						  GDsp2_dw_new,detGDsp2_dw_new,IGDsp2_dw_new)
						GDsp2_dw_old(num,1:H_N_part)=GDsp2_dw_new(num,1:H_N_part)
						IGDsp2_dw_old=IGDsp2_dw_new
						pvtgdsp2_dw_old=pvtgdsp2_dw_new
					END IF
				ELSE IF (tipo=='se1') THEN
					IF (num==-1) THEN
						GDsp1_up_old=GDsp1_up_new
						IGDsp1_up_old=IGDsp1_up_new
						pvtgdsp1_up_old=pvtgdsp1_up_new
						GDsp1_dw_old=GDsp1_dw_new
						IGDsp1_dw_old=IGDsp1_dw_new
						pvtgdsp1_dw_old=pvtgdsp1_dw_new
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num,IGDsp1_up_old,detGDsp1_up_old, &
						  GDsp1_up_new,detGDsp1_up_new,IGDsp1_up_new)
						GDsp1_up_old(num,1:H_N_part)=GDsp1_up_new(num,1:H_N_part)
						IGDsp1_up_old=IGDsp1_up_new
						pvtgdsp1_up_old=pvtgdsp1_up_new
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num-H_N_part,IGDsp1_dw_old,detGDsp1_dw_old, &
						  GDsp1_dw_new,detGDsp1_dw_new,IGDsp1_dw_new)
						GDsp1_dw_old(num,1:H_N_part)=GDsp1_dw_new(num,1:H_N_part)
						IGDsp1_dw_old=IGDsp1_dw_new
						pvtgdsp1_dw_old=pvtgdsp1_dw_new
					END IF
				ELSE IF (tipo=='se2') THEN
					IF (num==-1) THEN
						GDsp2_up_old=GDsp2_up_new
						IGDsp2_up_old=IGDsp2_up_new
						pvtgdsp2_up_old=pvtgdsp2_up_new
						GDsp2_dw_old=GDsp2_dw_new
						IGDsp2_dw_old=IGDsp2_dw_new
						pvtgdsp2_dw_old=pvtgdsp2_dw_new
					ELSE IF ((num>0) .AND. (num<=H_N_part)) THEN
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num,IGDsp2_up_old,detGDsp2_up_old, &
						  GDsp2_up_new,detGDsp2_up_new,IGDsp2_up_new)
						GDsp2_up_old(num,1:H_N_part)=GDsp2_up_new(num,1:H_N_part)
						IGDsp2_up_old=IGDsp2_up_new
						pvtgdsp2_up_old=pvtgdsp2_up_new
					ELSE IF ((num>H_N_part) .AND. (num<=N_part)) THEN
						CALL aggiorna_matrice_inversa_R_1ppt(H_N_part,num-H_N_part,IGDsp2_dw_old,detGDsp2_dw_old, &
						  GDsp2_dw_new,detGDsp2_dw_new,IGDsp2_dw_new)
						GDsp2_dw_old(num,1:H_N_part)=GDsp2_dw_new(num,1:H_N_part)
						IGDsp2_dw_old=IGDsp2_dw_new
						pvtgdsp2_dw_old=pvtgdsp2_dw_new
					END IF
				END IF
				detGDsp1_up_old=detGDsp1_up_new
				detGDsp1_dw_old=detGDsp1_dw_new
				detGDsp2_up_old=detGDsp2_up_new
				detGDsp2_dw_old=detGDsp2_dw_new
			END IF
			Usesp1_old=Usesp1_new
			Usesp2_old=Usesp2_new
		END IF
		
		IF (verbose_mode) PRINT *, 'calcola_accettazione: Ho aggiornato la funzione d onda'
				
	END SUBROUTINE aggiorna_funzione_onda
!-----------------------------------------------------------------------

	SUBROUTINE chiudi_funzione_onda()
		USE walkers
		USE funzione_onda
		USE momenta
		IMPLICIT NONE
		
		IF (.NOT. iniz_calcola_accettazione) STOP 'Prima di chiudere la funzione d onda devi inizializzarla &
		  [ module_calcola_accettazione.f90 > chiudi_funzione_onda ]'
		IF (SDe_kind/='no_') THEN
			DEALLOCATE(SDe_up_new, SDe_up_old)
			DEALLOCATE(ISDe_up_new, ISDe_up_old)
			DEALLOCATE(pvte_up_new, pvte_up_old)
			DEALLOCATE(SDe_dw_new, SDe_dw_old)
			DEALLOCATE(ISDe_dw_new, ISDe_dw_old)
			DEALLOCATE(pvte_dw_new, pvte_dw_old)
		END IF
		IF (Jee_kind/='no_') THEN
			DEALLOCATE(u_ee_new, u_ee_old)
		END IF
		IF (Jep_kind/='no_') THEN
			DEALLOCATE(u_ep_new, u_ep_old)
		END IF
		IF (SDse_kind/='no_') THEN
			DEALLOCATE(SDse1_up_new, SDse1_up_old)
			DEALLOCATE(ISDse1_up_new, ISDse1_up_old)
			DEALLOCATE(pvtse1_up_new, pvtse1_up_old)
			DEALLOCATE(SDse1_dw_new, SDse1_dw_old)
			DEALLOCATE(ISDse1_dw_new, ISDse1_dw_old)
			DEALLOCATE(pvtse1_dw_new, pvtse1_dw_old)
			DEALLOCATE(SDse2_up_new, SDse2_up_old)
			DEALLOCATE(ISDse2_up_new, ISDse2_up_old)
			DEALLOCATE(pvtse2_up_new, pvtse2_up_old)
			DEALLOCATE(SDse2_dw_new, SDse2_dw_old)
			DEALLOCATE(ISDse2_dw_new, ISDse2_dw_old)
			DEALLOCATE(pvtse2_dw_new, pvtse2_dw_old)
		END IF
		IF (Jse_kind/='no_') THEN
			IF (((Jse_kind/='no_') .AND. (Jse_kind/='bou') .AND. (Jse_kind/='pot'))) THEN
				DEALLOCATE(u_se1_new, u_se1_old)
				DEALLOCATE(u_se2_new, u_se2_old)
			END IF
			IF ((Jse_kind=='bou') .OR. (Jse_kind=='ppb')) THEN
				DEALLOCATE(partner_se1,partner_se2)
				DEALLOCATE(b_se1_new, b_se1_old)
				DEALLOCATE(b_se2_new, b_se2_old)
			END IF
			IF ( Jse_kind=='pot' ) THEN
				DEALLOCATE(u_POT_se1_new, u_POT_se1_old)
				DEALLOCATE(u_POT_se2_new, u_POT_se2_old)
			END IF
		END IF
		IF ((Kse_kind=='gss').OR.(Kse_kind=='gsc').OR.(Kse_kind=='gsp')) THEN
			DEALLOCATE(u_ese1_new, u_ese1_old,u_ese2_new,u_ese2_old)
		ELSE IF ((Kse_kind=='gsd').OR.(Kse_kind=='gdc').OR.(Kse_kind=='gdp')) THEN
			DEALLOCATE(GDse1_up_new, GDse1_up_old)
			DEALLOCATE(IGDse1_up_new, IGDse1_up_old)
			DEALLOCATE(pvtgdse1_up_new, pvtgdse1_up_old)
			DEALLOCATE(GDse1_dw_new, GDse1_dw_old)
			DEALLOCATE(IGDse1_dw_new, IGDse1_dw_old)
			DEALLOCATE(pvtgdse1_dw_new, pvtgdse1_dw_old)
			DEALLOCATE(GDse2_up_new, GDse2_up_old)
			DEALLOCATE(IGDse2_up_new, IGDse2_up_old)
			DEALLOCATE(pvtgdse2_up_new, pvtgdse2_up_old)
			DEALLOCATE(GDse2_dw_new, GDse2_dw_old)
			DEALLOCATE(IGDse2_dw_new, IGDse2_dw_old)
			DEALLOCATE(pvtgdse2_dw_new, pvtgdse2_dw_old)
		ELSE IF ((Kse_kind=='atm').OR.(Kse_kind=='atc')) THEN
			DEALLOCATE(u_ese1_new, u_ese1_old,u_ese2_new,u_ese2_old)
		END IF
		IF ((Jsesp_kind/='no_').AND.(Jsesp_kind/='gsd')) THEN
			DEALLOCATE(u_sesp1_new,u_sesp1_old,u_sesp2_new,u_sesp2_old)
		ELSE IF (Jsesp_kind=='gsd') THEN
			DEALLOCATE(GDsp1_up_new, GDsp1_up_old)
			DEALLOCATE(IGDsp1_up_new, IGDsp1_up_old)
			DEALLOCATE(pvtgdsp1_up_new, pvtgdsp1_up_old)
			DEALLOCATE(GDsp1_dw_new, GDsp1_dw_old)
			DEALLOCATE(IGDsp1_dw_new, IGDsp1_dw_old)
			DEALLOCATE(pvtgdsp1_dw_new, pvtgdsp1_dw_old)
			DEALLOCATE(GDsp2_up_new, GDsp2_up_old)
			DEALLOCATE(IGDsp2_up_new, IGDsp2_up_old)
			DEALLOCATE(pvtgdsp2_up_new, pvtgdsp2_up_old)
			DEALLOCATE(GDsp2_dw_new, GDsp2_dw_old)
			DEALLOCATE(IGDsp2_dw_new, IGDsp2_dw_old)
			DEALLOCATE(pvtgdsp2_dw_new, pvtgdsp2_dw_old)
		END IF
		CALL chiudi_momenta()
		IF (verbose_mode) PRINT *, 'calcola_accettazione: chiudo la funzione d onda'
		
		CALL chiudi_dati_funzione_onda()
				
		iniz_calcola_accettazione=.FALSE.
	END SUBROUTINE chiudi_funzione_onda
	
END MODULE calcola_accettazione