MODULE walkers
	IMPLICIT NONE
	LOGICAL, PROTECTED, SAVE :: iniz_walkers=.FALSE., iniz_pc=.FALSE.
	LOGICAL :: flag_traccia_coppie_mol_ss
	REAL (KIND=8), ALLOCATABLE, SAVE :: re_new(:,:), re_old(:,:), rp_new(:,:), rp_old(:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: se1_new(:,:), se1_old(:,:), se2_new(:,:), se2_old(:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: sp1_new(:,:), sp1_old(:,:), sp2_new(:,:), sp2_old(:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_ee_new(:,:,:), rij_ee_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_ee_new(:,:,:), rijpc_ee_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_se1_new(:,:,:), rij_se1_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_se1_new(:,:,:), rijpc_se1_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_se2_new(:,:,:), rij_se2_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_se2_new(:,:,:), rijpc_se2_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_ep_new(:,:,:), rij_ep_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_ep_new(:,:,:), rijpc_ep_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_ese1_new(:,:,:), rij_ese1_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_ese2_new(:,:,:), rij_ese2_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_ese1_new(:,:,:), rijpc_ese1_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_ese2_new(:,:,:), rijpc_ese2_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_sesp1_new(:,:,:), rij_sesp1_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_sesp2_new(:,:,:), rij_sesp2_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_sesp1_new(:,:,:), rijpc_sesp1_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_sesp2_new(:,:,:), rijpc_sesp2_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rij_pp_new(:,:,:), rij_pp_old(:,:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: rijpc_pp_new(:,:,:), rijpc_pp_old(:,:,:)
	INTEGER, ALLOCATABLE, PROTECTED, SAVE :: index_mol_ss1(:,:), index_mol_ss2(:,:)
	INTEGER, PROTECTED, SAVE :: index_mol_num
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: dist_mol_ss1_new(:), dist_mol_ss2_new(:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: dist_mol_ss1_old(:), dist_mol_ss2_old(:)
	
	CONTAINS
	
	SUBROUTINE inizializza_walkers(codice)
		USE dati_fisici
		USE dati_mc
		USE generic_tools
		IMPLICIT NONE
		LOGICAL :: flag_file_pos
		CHARACTER(LEN=*) :: codice
		CHARACTER(LEN=4) :: istring
		INTEGER :: i, j
		REAL (KIND=8) :: net_dist
		REAL (KIND=8), ALLOCATABLE :: rij_crystal(:,:,:)
		
		WRITE (istring, '(I4.4)'), mpi_myrank
		IF (.NOT. iniz_dati_fisici) STOP 'Non hai inizializzato dati_fisici &
		  [ module_walkers.f90 > inizializza_walkers ]'
		
		ALLOCATE(re_new(1:3,1:N_part),re_old(1:3,1:N_part))
		ALLOCATE(se1_new(1:3,1:N_part),se1_old(1:3,1:N_part),se2_new(1:3,1:N_part),se2_old(1:3,1:N_part))
		ALLOCATE(rij_ee_new(0:3,1:N_part,1:N_part),rij_ee_old(0:3,1:N_part,1:N_part))
		ALLOCATE(rij_se1_new(0:3,1:N_part,1:N_part),rij_se1_old(0:3,1:N_part,1:N_part))
		ALLOCATE(rij_se2_new(0:3,1:N_part,1:N_part),rij_se2_old(0:3,1:N_part,1:N_part))
		IF (flag_continua) THEN
			INQUIRE(FILE='posizioni/fine'//codice//'_e_'//istring//'.pos',EXIST=flag_file_pos)
			IF (flag_file_pos) THEN
				CALL leggi_posizioni_da_file(re_old,N_part,'posizioni/fine'//codice//'_e_'//istring)
				IF (flag_normalizza_pos) THEN
					re_old(1,1:N_part)=(re_old(1,1:N_part))*L(1)
					re_old(2,1:N_part)=(re_old(2,1:N_part))*L(2)
					re_old(3,1:N_part)=(re_old(3,1:N_part))*L(3)
					CALL applica_pbc(re_old,N_part,L)
				END IF
				IF (flag_shadow) THEN
					CALL leggi_posizioni_da_file(se1_old,N_part,'posizioni/fine'//codice//'_se1_'//istring)
					IF (flag_normalizza_pos) THEN
						se1_old(1,1:N_part)=(se1_old(1,1:N_part))*L(1)
						se1_old(2,1:N_part)=(se1_old(2,1:N_part))*L(2)
						se1_old(3,1:N_part)=(se1_old(3,1:N_part))*L(3)
						CALL applica_pbc(se1_old,N_part,L)
					END IF
					CALL leggi_posizioni_da_file(se2_old,N_part,'posizioni/fine'//codice//'_se2_'//istring)
					IF (flag_normalizza_pos) THEN
						se2_old(1,1:N_part)=(se2_old(1,1:N_part))*L(1)
						se2_old(2,1:N_part)=(se2_old(2,1:N_part))*L(2)
						se2_old(3,1:N_part)=(se2_old(3,1:N_part))*L(3)
						CALL applica_pbc(se2_old,N_part,L)
					END IF
				END IF
			ELSE
				STOP 'Non puoi continuare dalle posizioni precedenti, manca il file fine_e.pos &
				  [ module_walkers.f90 > inizializza_walkers ]'
			END IF
		ELSE
			!!! OPT1
			!ALLOCATE(rij_crystal(0:3,1:N_part,1:N_part))
			!CALL valuta_distanza_ii(r_crystal,N_part,L,rij_crystal)
			!net_dist=0.d0
			!DO j = 1, N_part-1, 1
			!	DO i = j+1, N_part, 1
			!		net_dist=MIN(net_dist,rij_crystal(i,j))
			!	END DO
			!END DO
			!!! OPT2
			!net_dist=MINVAL(L)/(N_part**(1.d0/3.d0))
			!!! OPT3, uso il raggio di Bohr
			net_dist=1.0d0
			IF (crystal_cell=='mol__') net_dist=0.74d0
			CALL RANDOM_NUMBER(re_new)
			DO j = 1, N_part, 1
				DO i = 1, 3, 1
					DO WHILE (DABS(re_new(i,j)-0.5d0)<0.001d0)
						CALL RANDOM_NUMBER(re_new(i,j))
					END DO
				END DO
			END DO
			re_old=r_crystal+(re_new-0.5d0)*net_dist
			CALL applica_pbc(re_old,N_part,L)
			IF (flag_shadow) THEN
				CALL RANDOM_NUMBER(se1_old)
				se1_old=(se1_old-0.5d0)*0.01d0*net_dist
				se1_old=re_old+se1_old
				CALL RANDOM_NUMBER(se2_old)
				se2_old=(se2_old-0.5d0)*0.01d0*net_dist
				se2_old=re_old+se2_old
			END IF
		END IF
		re_new=re_old
		CALL valuta_distanza_ii(re_old,N_part,L,rij_ee_old)
		rij_ee_new=rij_ee_old
		IF (flag_shadow) THEN
			se1_new=se1_old
			se2_new=se2_old
			CALL valuta_distanza_ii(se1_old,N_part,L,rij_se1_old)
			CALL valuta_distanza_ii(se2_old,N_part,L,rij_se2_old)
			rij_se1_new=rij_se1_old
			rij_se2_new=rij_se2_old
		END IF
		
		ALLOCATE(rp_new(1:3,1:N_part),rp_old(1:3,1:N_part))
		ALLOCATE(sp1_new(1:3,1:N_part),sp1_old(1:3,1:N_part),sp2_new(1:3,1:N_part),sp2_old(1:3,1:N_part))
		ALLOCATE(rij_pp_new(0:3,1:N_part,1:N_part),rij_pp_old(0:3,1:N_part,1:N_part))
		IF (flag_continua) THEN
			INQUIRE(FILE='posizioni/fine'//codice//'_p_'//istring//'.pos',EXIST=flag_file_pos)
			IF (flag_file_pos) THEN
				CALL leggi_posizioni_da_file(rp_old,N_part,'posizioni/fine'//codice//'_p_'//istring)
				IF (flag_normalizza_pos) THEN
					rp_old(1,1:N_part)=(rp_old(1,1:N_part))*L(1)
					rp_old(2,1:N_part)=(rp_old(2,1:N_part))*L(2)
					rp_old(3,1:N_part)=(rp_old(3,1:N_part))*L(3)
					CALL applica_pbc(rp_old,N_part,L)
				END IF
			ELSE
				STOP 'Non puoi continuare dalle posizioni precedenti, manca il file fine_p.pos &
				  [ module_walkers.f90 > inizializza_walkers ]'
			END IF
		ELSE
			rp_old=r_crystal
		END IF
		rp_new=rp_old

		CALL valuta_distanza_ii(rp_old,N_part,L,rij_pp_old)
		rij_pp_new=rij_pp_old
		IF (flag_shadow) THEN
			sp1_old=rp_old
			sp2_old=rp_old
			sp1_new=sp1_old
			sp2_new=sp2_old
		END IF
		
		ALLOCATE(rij_ep_new(0:3,1:N_part,1:N_part),rij_ep_old(0:3,1:N_part,1:N_part))
		ALLOCATE(rij_ese1_new(0:3,1:N_part,1:N_part),rij_ese1_old(0:3,1:N_part,1:N_part))
		ALLOCATE(rij_ese2_new(0:3,1:N_part,1:N_part),rij_ese2_old(0:3,1:N_part,1:N_part))
		ALLOCATE(rij_sesp1_new(0:3,1:N_part,1:N_part),rij_sesp1_old(0:3,1:N_part,1:N_part))
		ALLOCATE(rij_sesp2_new(0:3,1:N_part,1:N_part),rij_sesp2_old(0:3,1:N_part,1:N_part))
		CALL valuta_distanza_ij(re_old,rp_old,N_part,L,rij_ep_old)
		CALL valuta_distanza_ij(re_old,se1_old,N_part,L,rij_ese1_old)
		CALL valuta_distanza_ij(re_old,se2_old,N_part,L,rij_ese2_old)
		CALL valuta_distanza_ij(se1_old,sp1_old,N_part,L,rij_sesp1_old)
		CALL valuta_distanza_ij(se2_old,sp2_old,N_part,L,rij_sesp2_old)
		rij_ep_new=rij_ep_old
		rij_ese1_new=rij_ese1_old
		rij_ese2_new=rij_ese2_old
		rij_sesp1_new=rij_sesp1_old
		rij_sesp2_new=rij_sesp2_old
		
		iniz_walkers=.TRUE.
	END SUBROUTINE inizializza_walkers
!-----------------------------------------------------------------------

	SUBROUTINE attiva_traccia_coppie_mol_ss()
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		
		flag_traccia_coppie_mol_ss=.TRUE.
		ALLOCATE(index_mol_ss1(1:H_N_part,1:2),dist_mol_ss1_new(1:H_N_part),dist_mol_ss1_old(1:H_N_part))
		ALLOCATE(index_mol_ss2(1:H_N_part,1:2),dist_mol_ss2_new(1:H_N_part),dist_mol_ss2_old(1:H_N_part))
		CALL coppie_minima_distanza_stoc(N_part,rij_se1_new(0,1:N_part,1:N_part),index_mol_ss1,dist_mol_ss1_new)
		dist_mol_ss1_old=dist_mol_ss1_new
		CALL coppie_minima_distanza_stoc(N_part,rij_se2_new(0,1:N_part,1:N_part),index_mol_ss2,dist_mol_ss2_new)
		dist_mol_ss2_old=dist_mol_ss2_new
		
	END SUBROUTINE attiva_traccia_coppie_mol_ss
!-----------------------------------------------------------------------

	SUBROUTINE attiva_pc()
		USE dati_fisici
		USE dati_mc
		USE generic_tools
		IMPLICIT NONE
		IF (.NOT. iniz_dati_fisici) STOP 'Non hai inizializzato dati_fisici &
		  [ module_walkers.f90 > attiva_pc ]'
		IF (.NOT. iniz_pc) THEN
			IF (flag_elettroni) THEN
				ALLOCATE(rijpc_ee_new(0:3,1:N_part,1:N_part),rijpc_ee_old(0:3,1:N_part,1:N_part))
				CALL valuta_distanza_pc_ii(rij_ee_old,N_part,L,rijpc_ee_old)
				rijpc_ee_new=rijpc_ee_old
				ALLOCATE(rijpc_ep_new(0:3,1:N_part,1:N_part),rijpc_ep_old(0:3,1:N_part,1:N_part))
				CALL valuta_distanza_pc_ij(rij_ep_old,N_part,L,rijpc_ep_old)
				rijpc_ep_new=rijpc_ep_old
				IF (flag_shadow) THEN
					ALLOCATE(rijpc_se1_new(0:3,1:N_part,1:N_part),rijpc_se1_old(0:3,1:N_part,1:N_part))
					ALLOCATE(rijpc_se2_new(0:3,1:N_part,1:N_part),rijpc_se2_old(0:3,1:N_part,1:N_part))
					CALL valuta_distanza_pc_ii(rij_se1_old,N_part,L,rijpc_se1_old)
					rijpc_se1_new=rijpc_se1_old
					CALL valuta_distanza_pc_ii(rij_se2_old,N_part,L,rijpc_se2_old)
					rijpc_se2_new=rijpc_se2_old
					ALLOCATE(rijpc_sesp1_new(0:3,1:N_part,1:N_part),rijpc_sesp1_old(0:3,1:N_part,1:N_part))
					ALLOCATE(rijpc_sesp2_new(0:3,1:N_part,1:N_part),rijpc_sesp2_old(0:3,1:N_part,1:N_part))
					CALL valuta_distanza_pc_ij(rij_sesp1_old,N_part,L,rijpc_sesp1_old)
					rijpc_sesp1_new=rijpc_sesp1_old
					CALL valuta_distanza_pc_ij(rij_sesp2_old,N_part,L,rijpc_sesp2_old)
					rijpc_sesp2_new=rijpc_sesp2_old
					ALLOCATE(rijpc_ese1_new(0:3,1:N_part,1:N_part),rijpc_ese1_old(0:3,1:N_part,1:N_part))
					ALLOCATE(rijpc_ese2_new(0:3,1:N_part,1:N_part),rijpc_ese2_old(0:3,1:N_part,1:N_part))
					CALL valuta_distanza_pc_ij(rij_ese1_old,N_part,L,rijpc_ese1_old)
					rijpc_ese1_new=rijpc_ese1_old
					CALL valuta_distanza_pc_ij(rij_ese2_old,N_part,L,rijpc_ese2_old)
					rijpc_ese2_new=rijpc_ese2_old
				END IF
			END IF
			IF (flag_protoni) THEN
				ALLOCATE(rijpc_pp_new(0:3,1:N_part,1:N_part),rijpc_pp_old(0:3,1:N_part,1:N_part))
				CALL valuta_distanza_pc_ii(rij_pp_old,N_part,L,rijpc_pp_old)
				rijpc_pp_new=rijpc_pp_old
			END IF
		END IF
		iniz_pc=.TRUE.
	END SUBROUTINE attiva_pc
!-----------------------------------------------------------------------

	SUBROUTINE proponi_mossa(tipo,num)
		USE dati_mc
		USE dati_fisici
		USE generic_tools
		IMPLICIT NONE
		CHARACTER(LEN=3), INTENT(IN) :: tipo   !'e__', 'p__', 'se_', 'se1', 'se2'
		INTEGER, INTENT(IN) :: num
		IF (.NOT. iniz_dati_fisici) STOP 'Non hai inizializzato dati_fisici &
		  [ module_walkers.f90 > proponi_mossa ]'
		IF (.NOT. iniz_dati_mc) STOP 'Non hai inizializzato dati_mc &
		  [ module_walkers.f90 > proponi_mossa ]'
		IF (.NOT. iniz_walkers) STOP 'Prima di muovere i walkers devi inizializzarli &
		  [ module_walkers.f90 > proponi_mossa ]'
						
		SELECT CASE (tipo)
		CASE ('e__')
			IF ((howtomove=='allp') .AND. (num==-1)) THEN
				SELECT CASE (propmove)
				CASE ('flat')
					CALL RANDOM_NUMBER(re_new)
					re_new=(re_new-0.5d0)*step_e
				CASE ('gaus')
					CALL gaussian_sample(re_new,3,N_part,step_e)
				CASE DEFAULT
					STOP 'Il propmove non é accettabile'
				END SELECT
				IF (flag_shadow) THEN
					se1_new=se1_old+re_new
					CALL applica_pbc(se1_new,N_part,L)
					se2_new=se2_old+re_new
					CALL applica_pbc(se2_new,N_part,L)
					CALL calcola_nuove_distanze(tipo,num,'sese')
					CALL calcola_nuove_distanze(tipo,num,'sesp')
				END IF
				re_new=re_old+re_new
				CALL applica_pbc(re_new,N_part,L)
			ELSE IF ((howtomove=='1ppt') .AND. (num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (propmove)
				CASE ('flat')
					CALL RANDOM_NUMBER(re_new(1:3,num))
					re_new(1:3,num)=(re_new(1:3,num)-0.5d0)*step_e
				CASE ('gaus')
					CALL gaussian_sample(re_new(1:3,num),3,1,step_e)
				CASE DEFAULT
					STOP 'Il propmove non é accettabile'
				END SELECT
				re_new(1:3,num)=re_old(1:3,num)+re_new(1:3,num)
				CALL applica_pbc(re_new,N_part,L)
			END IF
			CALL calcola_nuove_distanze(tipo,num,'e_e_')
			CALL calcola_nuove_distanze(tipo,num,'e_p_')
			IF (flag_shadow) CALL calcola_nuove_distanze(tipo,num,'e_se')
		CASE ('se_')
			IF ((howtomove=='allp') .AND. (num==-1)) THEN
				SELECT CASE (propmove)
				CASE ('flat')
					CALL RANDOM_NUMBER(se1_new)
					se1_new=(se1_new-0.5d0)*step_se
					CALL RANDOM_NUMBER(se2_new)
					se2_new=(se2_new-0.5d0)*step_se
				CASE ('gaus')
					CALL gaussian_sample(se1_new,3,N_part,step_se)
					CALL gaussian_sample(se2_new,3,N_part,step_se)
				CASE DEFAULT
					STOP 'Il propmove non é accettabile'
				END SELECT
				se1_new=se1_old+se1_new
				CALL applica_pbc(se1_new,N_part,L)
				se2_new=se2_old+se2_new
				CALL applica_pbc(se2_new,N_part,L)
			END IF
			CALL calcola_nuove_distanze(tipo,num,'sese')
			CALL calcola_nuove_distanze(tipo,num,'sesp')
			CALL calcola_nuove_distanze(tipo,num,'e_se')
		CASE ('se1')
			IF ((howtomove=='1ppt') .AND. (num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (propmove)
				CASE ('flat')
					CALL RANDOM_NUMBER(se1_new(1:3,num))
					se1_new(1:3,num)=(se1_new(1:3,num)-0.5d0)*step_se
				CASE ('gaus')
					CALL gaussian_sample(se1_new(1:3,num),3,1,step_se)
				CASE DEFAULT
					STOP 'Il propmove non é accettabile'
				END SELECT
				se1_new(1:3,num)=se1_old(1:3,num)+se1_new(1:3,num)
				CALL applica_pbc(se1_new,N_part,L)
			END IF
			CALL calcola_nuove_distanze(tipo,num,'sese')
			CALL calcola_nuove_distanze(tipo,num,'sesp')
			CALL calcola_nuove_distanze(tipo,num,'e_se')
		CASE ('se2')
			IF ((howtomove=='1ppt') .AND. (num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (propmove)
				CASE ('flat')
					CALL RANDOM_NUMBER(se2_new(1:3,num))
					se2_new(1:3,num)=(se2_new(1:3,num)-0.5d0)*step_se
				CASE ('gaus')
					CALL gaussian_sample(se2_new(1:3,num),3,1,step_se)
				CASE DEFAULT
					STOP 'Il propmove non é accettabile'
				END SELECT
				se2_new(1:3,num)=se2_old(1:3,num)+se2_new(1:3,num)
				CALL applica_pbc(se2_new,N_part,L)
			END IF
			CALL calcola_nuove_distanze(tipo,num,'sese')
			CALL calcola_nuove_distanze(tipo,num,'sesp')
			CALL calcola_nuove_distanze(tipo,num,'e_se')
		CASE ('tre')
			IF ((howtomove=='1ppt') .AND. (num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (propmove)
				CASE ('flat')
					CALL RANDOM_NUMBER(re_new(1:3,num))
					re_new(1:3,num)=(re_new(1:3,num)-0.5d0)*step_e
				CASE ('gaus')
					CALL gaussian_sample(re_new(1:3,num),3,1,step_e)
				CASE DEFAULT
					STOP 'Il propmove non é accettabile'
				END SELECT
				se1_new(1:3,num)=se1_old(1:3,num)+re_new(1:3,num)
				CALL applica_pbc(se1_new,N_part,L)
				se2_new(1:3,num)=se2_old(1:3,num)+re_new(1:3,num)
				CALL applica_pbc(se2_new,N_part,L)
				re_new(1:3,num)=re_old(1:3,num)+re_new(1:3,num)
				CALL applica_pbc(re_new,N_part,L)
			END IF
			CALL calcola_nuove_distanze(tipo,num,'e_e_')
			CALL calcola_nuove_distanze(tipo,num,'e_p_')
			CALL calcola_nuove_distanze(tipo,num,'e_se')
			CALL calcola_nuove_distanze(tipo,num,'sese')
			CALL calcola_nuove_distanze(tipo,num,'sesp')
		END SELECT
				
	END SUBROUTINE proponi_mossa
!-----------------------------------------------------------------------

	SUBROUTINE calcola_nuove_distanze(tipo,num,tra_chi)
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		CHARACTER(LEN=3), INTENT(IN) :: tipo              !tipo di particella mossa 'e__', 'p__', 'se1', 'se2'
		CHARACTER(LEN=4), INTENT(IN) :: tra_chi           !distanza fra quali tipi di particelle 'e_e_', 'e_p_', 'e_se', 'sese', 'sesp'
		INTEGER, INTENT(IN) :: num
		INTEGER :: i, j
		
		IF (.NOT. iniz_dati_fisici) STOP 'Non hai inizializzato dati_fisici &
		  [ module_walkers.f90 > calcola_nuove_distanze ]'
		IF (.NOT. iniz_walkers) STOP 'Prima di calcolare le distanze devi aver inizializzato i walkers &
		  [ module_walkers.f90 > calcola_nuove_distanze ]'
				
		SELECT CASE (tra_chi)
		CASE ('e_e_')
			IF (num==-1) THEN
				CALL valuta_distanza_ii(re_new,N_part,L,rij_ee_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				CALL valuta_distanza_ii_1ppt(num,re_new,N_part,L,rij_ee_new)
			END IF
		CASE ('p_p_')
			IF (num==-1) THEN
				CALL valuta_distanza_ii(rp_new,N_part,L,rij_pp_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				CALL valuta_distanza_ii_1ppt(num,rp_new,N_part,L,rij_pp_new)
			END IF
		CASE ('e_p_')
			IF (num==-1) THEN
				CALL valuta_distanza_ij(re_new,rp_new,N_part,L,rij_ep_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (tipo)
				CASE ('e__')
					CALL valuta_distanza_ij_1ppt(num,1,re_new,rp_new,N_part,L,rij_ep_new)
				CASE ('p__')
					CALL valuta_distanza_ij_1ppt(num,2,re_new,rp_new,N_part,L,rij_ep_new)
				CASE ('tre')
					CALL valuta_distanza_ij_1ppt(num,1,re_new,rp_new,N_part,L,rij_ep_new)
				END SELECT
			END IF
		CASE ('sese')
			IF (num==-1) THEN
				CALL valuta_distanza_ii(se1_new,N_part,L,rij_se1_new)
				CALL valuta_distanza_ii(se2_new,N_part,L,rij_se2_new)
				IF ( flag_traccia_coppie_mol_ss ) THEN
					CALL coppie_minima_distanza_stoc(N_part,rij_se1_new(0,1:N_part,1:N_part),index_mol_ss1,dist_mol_ss1_new)
					CALL coppie_minima_distanza_stoc(N_part,rij_se2_new(0,1:N_part,1:N_part),index_mol_ss2,dist_mol_ss2_new)
					index_mol_num=-1
				END IF
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (tipo)
				CASE ('se1')
					CALL valuta_distanza_ii_1ppt(num,se1_new,N_part,L,rij_se1_new)
					IF (flag_traccia_coppie_mol_ss) THEN
					 	CALL aggiorna_coppie_minima_distanza_stoc_1ppt(num,N_part,&
					         rij_se1_new(0,1:N_part,1:N_part),index_mol_ss1,dist_mol_ss1_new,index_mol_num)
					END IF
				CASE ('se2')
					CALL valuta_distanza_ii_1ppt(num,se2_new,N_part,L,rij_se2_new)
					IF (flag_traccia_coppie_mol_ss) THEN
					 	CALL aggiorna_coppie_minima_distanza_stoc_1ppt(num,N_part,&
							rij_se2_new(0,1:N_part,1:N_part),index_mol_ss2,dist_mol_ss2_new,index_mol_num)
					END IF
				CASE ('tre')
					!se1
					CALL valuta_distanza_ii_1ppt(num,se1_new,N_part,L,rij_se1_new)
					IF (flag_traccia_coppie_mol_ss) THEN
					 	CALL aggiorna_coppie_minima_distanza_stoc_1ppt(num,N_part,&
					         rij_se1_new(0,1:N_part,1:N_part),index_mol_ss1,dist_mol_ss1_new,index_mol_num)
					END IF
					!se2
					CALL valuta_distanza_ii_1ppt(num,se2_new,N_part,L,rij_se2_new)
					IF (flag_traccia_coppie_mol_ss) THEN
					 	CALL aggiorna_coppie_minima_distanza_stoc_1ppt(num,N_part,&
							rij_se2_new(0,1:N_part,1:N_part),index_mol_ss2,dist_mol_ss2_new,index_mol_num)
					END IF
				END SELECT
			END IF
		CASE ('e_se')
			IF (num==-1) THEN
				CALL valuta_distanza_ij(re_new,se1_new,N_part,L,rij_ese1_new)
				CALL valuta_distanza_ij(re_new,se2_new,N_part,L,rij_ese2_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (tipo)
				CASE ('e__')
					CALL valuta_distanza_ij_1ppt(num,1,re_new,se1_new,N_part,L,rij_ese1_new)
					CALL valuta_distanza_ij_1ppt(num,1,re_new,se2_new,N_part,L,rij_ese2_new)
				CASE ('se1')
					CALL valuta_distanza_ij_1ppt(num,2,re_new,se1_new,N_part,L,rij_ese1_new)
				CASE ('se2')
					CALL valuta_distanza_ij_1ppt(num,2,re_new,se2_new,N_part,L,rij_ese2_new)
				CASE ('tre')
					!re
					CALL valuta_distanza_ij_1ppt(num,1,re_new,se1_new,N_part,L,rij_ese1_new)
					CALL valuta_distanza_ij_1ppt(num,1,re_new,se2_new,N_part,L,rij_ese2_new)
					!se1
					CALL valuta_distanza_ij_1ppt(num,2,re_new,se1_new,N_part,L,rij_ese1_new)
					!se2
					CALL valuta_distanza_ij_1ppt(num,2,re_new,se2_new,N_part,L,rij_ese2_new)
				END SELECT
			END IF
		CASE ('sesp')
			IF (num==-1) THEN
				CALL valuta_distanza_ij(se1_new,sp1_new,N_part,L,rij_sesp1_new)
				CALL valuta_distanza_ij(se2_new,sp2_new,N_part,L,rij_sesp2_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (tipo)
				CASE ('se1')
					CALL valuta_distanza_ij_1ppt(num,1,se1_new,sp1_new,N_part,L,rij_sesp1_new)
				CASE ('se2')
					CALL valuta_distanza_ij_1ppt(num,1,se2_new,sp2_new,N_part,L,rij_sesp2_new)
				CASE ('tre')
					!se1
					CALL valuta_distanza_ij_1ppt(num,1,se1_new,sp1_new,N_part,L,rij_sesp1_new)
					!se2
					CALL valuta_distanza_ij_1ppt(num,1,se2_new,sp2_new,N_part,L,rij_sesp2_new)
				END SELECT
			END IF			
		END SELECT
						
	END SUBROUTINE calcola_nuove_distanze
!-----------------------------------------------------------------------

	SUBROUTINE ritraccia_coppie_molecolari_ss()
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		
		CALL coppie_updw_minima_distanza_stoc(N_part,rij_se1_new(0,1:N_part,1:N_part),index_mol_ss1,dist_mol_ss1_new)
		CALL coppie_updw_minima_distanza_stoc(N_part,rij_se2_new(0,1:N_part,1:N_part),index_mol_ss2,dist_mol_ss2_new)
		
	END SUBROUTINE ritraccia_coppie_molecolari_ss
!-----------------------------------------------------------------------

	SUBROUTINE calcola_nuove_distanze_pc(tipo,num,tra_chi)
		USE generic_tools
		USE dati_fisici
		IMPLICIT NONE
		CHARACTER(LEN=3), INTENT(IN) :: tipo
		CHARACTER(LEN=4), INTENT(IN) :: tra_chi
		INTEGER, INTENT(IN) :: num
		IF (.NOT. iniz_dati_fisici) STOP 'Non hai inizializzato dati_fisici &
		  [ module_walkers.f90 > calcola_nuove_distanze_pc ]'
		IF (.NOT. iniz_walkers) STOP 'Prima di calcolare le distanze devi aver inizializzato i walkers &
		  [ module_walkers.f90 > calcola_nuove_distanze_pc ]'
		IF (.NOT. iniz_pc) STOP 'Prima di calcolare le distanze PC devi inizializzarle &
		  [ module_walkers.f90 > calcola_nuove_distanze_pc ]'
		
		SELECT CASE (tra_chi)
		CASE ('e_e_')
			IF (num==-1) THEN
				CALL valuta_distanza_pc_ii(rij_ee_new,N_part,L,rijpc_ee_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				CALL valuta_distanza_pc_ii_1ppt(num,rij_ee_new,N_part,L,rijpc_ee_new)
			END IF
			!!!WRITE(UNIT=166, FMT=*) rij_ee_new(1,1,2), rijpc_ee_new(1,1,2)
		CASE ('sese')
			IF (num==-1) THEN
				CALL valuta_distanza_pc_ii(rij_se1_new,N_part,L,rijpc_se1_new)
				CALL valuta_distanza_pc_ii(rij_se2_new,N_part,L,rijpc_se2_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (tipo)
				CASE ('se1')
					CALL valuta_distanza_pc_ii_1ppt(num,rij_se1_new,N_part,L,rijpc_se1_new)
				CASE ('se2')
					CALL valuta_distanza_pc_ii_1ppt(num,rij_se2_new,N_part,L,rijpc_se2_new)
				CASE ('tre')
					CALL valuta_distanza_pc_ii_1ppt(num,rij_se1_new,N_part,L,rijpc_se1_new)
					CALL valuta_distanza_pc_ii_1ppt(num,rij_se2_new,N_part,L,rijpc_se2_new)
				END SELECT
			END IF
		CASE ('p_p_')
			IF (num==-1) THEN
				CALL valuta_distanza_pc_ii(rij_pp_new,N_part,L,rijpc_pp_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				CALL valuta_distanza_pc_ii_1ppt(num,rij_pp_new,N_part,L,rijpc_pp_new)
			END IF
		CASE ('e_p_')
			IF (num==-1) THEN
				CALL valuta_distanza_pc_ij(rij_ep_new,N_part,L,rijpc_ep_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (tipo)
				CASE ('e__')
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_ep_new,N_part,L,rijpc_ep_new)
				CASE ('p__')
					CALL valuta_distanza_pc_ij_1ppt(num,2,rij_ep_new,N_part,L,rijpc_ep_new)
				CASE ('tre')
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_ep_new,N_part,L,rijpc_ep_new)
				END SELECT
			END IF
			!!!WRITE(UNIT=666, FMT=*) rij_ep_new(1,1,2), rijpc_ep_new(1,1,2)
		CASE ('sesp')
			IF (num==-1) THEN
				CALL valuta_distanza_pc_ij(rij_sesp1_new,N_part,L,rijpc_sesp1_new)
				CALL valuta_distanza_pc_ij(rij_sesp2_new,N_part,L,rijpc_sesp2_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (tipo)
				CASE ('se1')
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_sesp1_new,N_part,L,rijpc_sesp1_new)
				CASE ('se2')
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_sesp2_new,N_part,L,rijpc_sesp2_new)
				CASE ('tre')
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_sesp1_new,N_part,L,rijpc_sesp1_new)
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_sesp2_new,N_part,L,rijpc_sesp2_new)
				END SELECT
			END IF
			!!!WRITE(UNIT=466, FMT=*) rij_sesp1_new(1,1,2), rijpc_sesp1_new(1,1,2)
			!!!WRITE(UNIT=566, FMT=*) rij_sesp2_new(1,1,2), rijpc_sesp2_new(1,1,2)
		CASE ('e_se')
			IF (num==-1) THEN
				CALL valuta_distanza_pc_ij(rij_ese1_new,N_part,L,rijpc_ese1_new)
				CALL valuta_distanza_pc_ij(rij_ese2_new,N_part,L,rijpc_ese2_new)
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				SELECT CASE (tipo)
				CASE ('e__')
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_ese1_new,N_part,L,rijpc_ese1_new)
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_ese2_new,N_part,L,rijpc_ese2_new)
				CASE ('se1')
					CALL valuta_distanza_pc_ij_1ppt(num,2,rij_ese1_new,N_part,L,rijpc_ese1_new)
				CASE ('se2')
					CALL valuta_distanza_pc_ij_1ppt(num,2,rij_ese2_new,N_part,L,rijpc_ese2_new)
				CASE ('tre')
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_ese1_new,N_part,L,rijpc_ese1_new)
					CALL valuta_distanza_pc_ij_1ppt(num,1,rij_ese2_new,N_part,L,rijpc_ese2_new)
					CALL valuta_distanza_pc_ij_1ppt(num,2,rij_ese1_new,N_part,L,rijpc_ese1_new)
					CALL valuta_distanza_pc_ij_1ppt(num,2,rij_ese2_new,N_part,L,rijpc_ese2_new)
				END SELECT
			END IF
		CASE DEFAULT
			STOP 'Hai richiesto un calcolo pc che non é stato considerato &
			  [ module_walkers.f90 > calcola_nuove_distanze_pc ]'
		END SELECT
	END SUBROUTINE calcola_nuove_distanze_pc
!-----------------------------------------------------------------------

	SUBROUTINE metti_rnew_in_rold(tipo,num)
		USE dati_mc
		USE dati_fisici
		IMPLICIT NONE
		CHARACTER(LEN=3), INTENT(IN) :: tipo
		INTEGER, INTENT(IN) :: num
		IF (.NOT. iniz_walkers) STOP 'Prima di farlo devi aver inizializzato i walkers &
		  [ module_walkers.f90 > metti_new_in_old ]'
		
		SELECT CASE (tipo)
		CASE ('e__')
			IF (num==-1) THEN
				re_old=re_new
				rij_ee_old=rij_ee_new
				rij_ep_old=rij_ep_new
				IF (iniz_pc) THEN
					rijpc_ee_old=rijpc_ee_new
					rijpc_ep_old=rijpc_ep_new
				END IF
				IF (flag_shadow) THEN
					rij_ese1_old=rij_ese1_new
					rij_ese2_old=rij_ese2_new
				END IF
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				re_old(1:3,num)=re_new(1:3,num)
				rij_ee_old(0:3,1:N_part,num)=rij_ee_new(0:3,1:N_part,num)
				rij_ee_old(0:3,num,1:N_part)=rij_ee_new(0:3,num,1:N_part)
				rij_ep_old(0:3,num,1:N_part)=rij_ep_new(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_ee_old(0:3,1:N_part,num)=rijpc_ee_new(0:3,1:N_part,num)
					rijpc_ee_old(0:3,num,1:N_part)=rijpc_ee_new(0:3,num,1:N_part)
					rijpc_ep_old(0:3,num,1:N_part)=rijpc_ep_new(0:3,num,1:N_part)
				END IF
				IF (flag_shadow) THEN
					rij_ese1_old(0:3,num,1:N_part)=rij_ese1_new(0:3,num,1:N_part)
					rij_ese2_old(0:3,num,1:N_part)=rij_ese2_new(0:3,num,1:N_part)
					IF ( iniz_pc ) THEN
						rijpc_ese1_old(0:3,num,1:N_part)=rijpc_ese1_new(0:3,num,1:N_part)
						rijpc_ese2_old(0:3,num,1:N_part)=rijpc_ese2_new(0:3,num,1:N_part)
					END IF
				END IF
			END IF
		CASE ('se_')
			IF (num==-1) THEN
				se1_old=se1_new
				rij_se1_old=rij_se1_new
				se2_old=se2_new
				rij_se2_old=rij_se2_new
				IF (flag_traccia_coppie_mol_ss) THEN
					dist_mol_ss1_old=dist_mol_ss1_new
					dist_mol_ss2_old=dist_mol_ss2_new
				END IF
				IF (iniz_pc) THEN
					rijpc_se1_old=rijpc_se1_new
					rijpc_se2_old=rijpc_se2_new
				END IF
				rij_ese1_old=rij_ese1_new
				rij_ese2_old=rij_ese2_new
				IF ( iniz_pc ) THEN
					rijpc_ese1_old=rijpc_ese1_new
					rijpc_ese2_old=rijpc_ese2_new
				END IF
				rij_sesp1_old=rij_sesp1_new
				rij_sesp2_old=rij_sesp2_new
				IF (iniz_pc) THEN
					rijpc_sesp1_old(0:3,num,1:N_part)=rijpc_sesp1_new(0:3,num,1:N_part)
					rijpc_sesp2_old(0:3,num,1:N_part)=rijpc_sesp2_new(0:3,num,1:N_part)
				END IF
			END IF
		CASE ('se1')
			IF ((num>0) .AND. (num<=N_part)) THEN
				se1_old(1:3,num)=se1_new(1:3,num)
				rij_se1_old(0:3,1:N_part,num)=rij_se1_new(0:3,1:N_part,num)
				rij_se1_old(0:3,num,1:N_part)=rij_se1_new(0:3,num,1:N_part)
				IF (flag_traccia_coppie_mol_ss) dist_mol_ss1_old(index_mol_num)=dist_mol_ss1_new(index_mol_num)
				IF (iniz_pc) THEN
					rijpc_se1_old(0:3,1:N_part,num)=rijpc_se1_new(0:3,1:N_part,num)
					rijpc_se1_old(0:3,num,1:N_part)=rijpc_se1_new(0:3,num,1:N_part)
				END IF
				rij_ese1_old(0:3,1:N_part,num)=rij_ese1_new(0:3,1:N_part,num)
				IF ( iniz_pc ) THEN
					rijpc_ese1_old(0:3,1:N_part,num)=rijpc_ese1_new(0:3,1:N_part,num)
				END IF
				rij_sesp1_old(0:3,num,1:N_part)=rij_sesp1_new(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_sesp1_old(0:3,num,1:N_part)=rijpc_sesp1_new(0:3,num,1:N_part)
				END IF
			END IF
		CASE ('se2')
			IF ((num>0) .AND. (num<=N_part)) THEN
				se2_old(1:3,num)=se2_new(1:3,num)
				rij_se2_old(0:3,1:N_part,num)=rij_se2_new(0:3,1:N_part,num)
				rij_se2_old(0:3,num,1:N_part)=rij_se2_new(0:3,num,1:N_part)
				IF (flag_traccia_coppie_mol_ss) dist_mol_ss2_old(index_mol_num)=dist_mol_ss2_new(index_mol_num)
				IF (iniz_pc) THEN
					rijpc_se2_old(0:3,1:N_part,num)=rijpc_se2_new(0:3,1:N_part,num)
					rijpc_se2_old(0:3,num,1:N_part)=rijpc_se2_new(0:3,num,1:N_part)
				END IF
				rij_ese2_old(0:3,1:N_part,num)=rij_ese2_new(0:3,1:N_part,num)
				IF ( iniz_pc ) THEN
					rijpc_ese2_old(0:3,1:N_part,num)=rijpc_ese2_new(0:3,1:N_part,num)
				END IF
				rij_sesp2_old(0:3,num,1:N_part)=rij_sesp2_new(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_sesp2_old(0:3,num,1:N_part)=rijpc_sesp2_new(0:3,num,1:N_part)
				END IF
			END IF
		CASE ('tre')
			!re
			IF (num==-1) THEN
				re_old=re_new
				rij_ee_old=rij_ee_new
				rij_ep_old=rij_ep_new
				IF (iniz_pc) THEN
					rijpc_ee_old=rijpc_ee_new
					rijpc_ep_old=rijpc_ep_new
				END IF
				IF (flag_shadow) THEN
					rij_ese1_old=rij_ese1_new
					rij_ese2_old=rij_ese2_new
				END IF
			ELSE IF ((num>0) .AND. (num<=N_part)) THEN
				re_old(1:3,num)=re_new(1:3,num)
				rij_ee_old(0:3,1:N_part,num)=rij_ee_new(0:3,1:N_part,num)
				rij_ee_old(0:3,num,1:N_part)=rij_ee_new(0:3,num,1:N_part)
				rij_ep_old(0:3,num,1:N_part)=rij_ep_new(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_ee_old(0:3,1:N_part,num)=rijpc_ee_new(0:3,1:N_part,num)
					rijpc_ee_old(0:3,num,1:N_part)=rijpc_ee_new(0:3,num,1:N_part)
					rijpc_ep_old(0:3,num,1:N_part)=rijpc_ep_new(0:3,num,1:N_part)
				END IF
				IF (flag_shadow) THEN
					rij_ese1_old(0:3,num,1:N_part)=rij_ese1_new(0:3,num,1:N_part)
					rij_ese2_old(0:3,num,1:N_part)=rij_ese2_new(0:3,num,1:N_part)
					IF ( iniz_pc ) THEN
						rijpc_ese1_old(0:3,num,1:N_part)=rijpc_ese1_new(0:3,num,1:N_part)
						rijpc_ese2_old(0:3,num,1:N_part)=rijpc_ese2_new(0:3,num,1:N_part)
					END IF
				END IF
			END IF
			!se1
			IF ((num>0) .AND. (num<=N_part)) THEN
				se1_old(1:3,num)=se1_new(1:3,num)
				rij_se1_old(0:3,1:N_part,num)=rij_se1_new(0:3,1:N_part,num)
				rij_se1_old(0:3,num,1:N_part)=rij_se1_new(0:3,num,1:N_part)
				IF (flag_traccia_coppie_mol_ss) dist_mol_ss1_old(index_mol_num)=dist_mol_ss1_new(index_mol_num)
				IF (iniz_pc) THEN
					rijpc_se1_old(0:3,1:N_part,num)=rijpc_se1_new(0:3,1:N_part,num)
					rijpc_se1_old(0:3,num,1:N_part)=rijpc_se1_new(0:3,num,1:N_part)
				END IF
				rij_ese1_old(0:3,1:N_part,num)=rij_ese1_new(0:3,1:N_part,num)
				IF ( iniz_pc ) THEN
					rijpc_ese1_old(0:3,1:N_part,num)=rijpc_ese1_new(0:3,1:N_part,num)
				END IF
				rij_sesp1_old(0:3,num,1:N_part)=rij_sesp1_new(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_sesp1_old(0:3,num,1:N_part)=rijpc_sesp1_new(0:3,num,1:N_part)
				END IF
			END IF
			!se2
			IF ((num>0) .AND. (num<=N_part)) THEN
				se2_old(1:3,num)=se2_new(1:3,num)
				rij_se2_old(0:3,1:N_part,num)=rij_se2_new(0:3,1:N_part,num)
				rij_se2_old(0:3,num,1:N_part)=rij_se2_new(0:3,num,1:N_part)
				IF (flag_traccia_coppie_mol_ss) dist_mol_ss2_old(index_mol_num)=dist_mol_ss2_new(index_mol_num)
				IF (iniz_pc) THEN
					rijpc_se2_old(0:3,1:N_part,num)=rijpc_se2_new(0:3,1:N_part,num)
					rijpc_se2_old(0:3,num,1:N_part)=rijpc_se2_new(0:3,num,1:N_part)
				END IF
				rij_ese2_old(0:3,1:N_part,num)=rij_ese2_new(0:3,1:N_part,num)
				IF ( iniz_pc ) THEN
					rijpc_ese2_old(0:3,1:N_part,num)=rijpc_ese2_new(0:3,1:N_part,num)
				END IF
				rij_sesp2_old(0:3,num,1:N_part)=rij_sesp2_new(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_sesp2_old(0:3,num,1:N_part)=rijpc_sesp2_new(0:3,num,1:N_part)
				END IF
			END IF
		CASE ('all')
			IF (num==-1) THEN
				IF (flag_elettroni) THEN
					re_old=re_new
					rij_ee_old=rij_ee_new
					rij_ep_old=rij_ep_new
					IF (iniz_pc) THEN
						rijpc_ep_old=rijpc_ep_new
						rijpc_ee_old=rijpc_ee_new
					END IF
					IF (flag_shadow) THEN
						se1_old=se1_new
						rij_se1_old=rij_se1_new
						se2_old=se2_new
						rij_se2_old=rij_se2_new
						IF (flag_traccia_coppie_mol_ss) THEN
							dist_mol_ss1_old=dist_mol_ss1_new
							dist_mol_ss2_old=dist_mol_ss2_new
						END IF
						IF (iniz_pc) THEN
							rijpc_se1_old=rijpc_se1_new
							rijpc_se2_old=rijpc_se2_new
						END IF
						rij_ese1_old=rij_ese1_new
						rij_ese2_old=rij_ese2_new
						IF ( iniz_pc ) THEN
							rijpc_ese1_old=rijpc_ese1_new
							rijpc_ese2_old=rijpc_ese2_new
						END IF
						rij_sesp1_old=rij_sesp1_new
						rij_sesp2_old=rij_sesp2_new
						IF (iniz_pc) THEN
							rijpc_sesp1_old=rijpc_sesp1_new
							rijpc_sesp2_old=rijpc_sesp2_new
						END IF
					END IF
				END IF
			END IF
		END SELECT
				
	END SUBROUTINE metti_rnew_in_rold
!-----------------------------------------------------------------------

	SUBROUTINE riporta_indietro_walker(tipo,num)
		USE dati_fisici
		USE dati_mc
		IMPLICIT NONE
		CHARACTER(LEN=3), INTENT(IN) :: tipo
		INTEGER, INTENT(IN) :: num
				
		IF (num>0 .AND. num<=N_part) THEN
			SELECT CASE (tipo)
			CASE ('e__')
				re_new(1:3,num)=re_old(1:3,num)
				rij_ee_new(0:3,1:N_part,num)=rij_ee_old(0:3,1:N_part,num)
				rij_ee_new(0:3,num,1:N_part)=rij_ee_old(0:3,num,1:N_part)
				rij_ep_new(0:3,num,1:N_part)=rij_ep_old(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_ee_new(0:3,1:N_part,num)=rijpc_ee_old(0:3,1:N_part,num)
					rijpc_ee_new(0:3,num,1:N_part)=rijpc_ee_old(0:3,num,1:N_part)
					rijpc_ep_new(0:3,num,1:N_part)=rijpc_ep_old(0:3,num,1:N_part)
				END IF
				IF (flag_shadow) THEN
					rij_ese1_new(0:3,num,1:N_part)=rij_ese1_old(0:3,num,1:N_part)
					rij_ese2_new(0:3,num,1:N_part)=rij_ese2_old(0:3,num,1:N_part)
					IF ( iniz_pc ) THEN
						rijpc_ese1_new(0:3,num,1:N_part)=rijpc_ese1_old(0:3,num,1:N_part)
						rijpc_ese2_new(0:3,num,1:N_part)=rijpc_ese2_old(0:3,num,1:N_part)
					END IF
				END IF
			CASE ('se1')
				se1_new(1:3,num)=se1_old(1:3,num)
				rij_se1_new(0:3,1:N_part,num)=rij_se1_old(0:3,1:N_part,num)
				rij_se1_new(0:3,num,1:N_part)=rij_se1_old(0:3,num,1:N_part)
				IF (flag_traccia_coppie_mol_ss) dist_mol_ss1_new(index_mol_num)=dist_mol_ss1_old(index_mol_num)
				IF (iniz_pc) THEN
					rijpc_se1_new(0:3,1:N_part,num)=rijpc_se1_old(0:3,1:N_part,num)
					rijpc_se1_new(0:3,num,1:N_part)=rijpc_se1_old(0:3,num,1:N_part)
				END IF
				rij_ese1_new(0:3,1:N_part,num)=rij_ese1_old(0:3,1:N_part,num)
				IF ( iniz_pc ) THEN
					rijpc_ese1_new(0:3,1:N_part,num)=rijpc_ese1_old(0:3,1:N_part,num)
				END IF
				rij_sesp1_new(0:3,num,1:N_part)=rij_sesp1_old(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_sesp1_new(0:3,num,1:N_part)=rijpc_sesp1_old(0:3,num,1:N_part)
				END IF
			CASE ('se2')
				se2_new(1:3,num)=se2_old(1:3,num)
				rij_se2_new(0:3,1:N_part,num)=rij_se2_old(0:3,1:N_part,num)
				rij_se2_new(0:3,num,1:N_part)=rij_se2_old(0:3,num,1:N_part)
				IF (flag_traccia_coppie_mol_ss) dist_mol_ss2_new(index_mol_num)=dist_mol_ss2_old(index_mol_num)
				IF (iniz_pc) THEN
					rijpc_se2_new(0:3,1:N_part,num)=rijpc_se2_old(0:3,1:N_part,num)
					rijpc_se2_new(0:3,num,1:N_part)=rijpc_se2_old(0:3,num,1:N_part)
				END IF
				rij_ese2_new(0:3,1:N_part,num)=rij_ese2_old(0:3,1:N_part,num)
				IF ( iniz_pc ) THEN
					rijpc_ese2_new(0:3,1:N_part,num)=rijpc_ese2_old(0:3,1:N_part,num)
				END IF
				rij_sesp2_new(0:3,num,1:N_part)=rij_sesp2_old(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_sesp2_new(0:3,num,1:N_part)=rijpc_sesp2_old(0:3,num,1:N_part)
				END IF
			CASE ('tre')
				!re
				re_new(1:3,num)=re_old(1:3,num)
				rij_ee_new(0:3,1:N_part,num)=rij_ee_old(0:3,1:N_part,num)
				rij_ee_new(0:3,num,1:N_part)=rij_ee_old(0:3,num,1:N_part)
				rij_ep_new(0:3,num,1:N_part)=rij_ep_old(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_ee_new(0:3,1:N_part,num)=rijpc_ee_old(0:3,1:N_part,num)
					rijpc_ee_new(0:3,num,1:N_part)=rijpc_ee_old(0:3,num,1:N_part)
					rijpc_ep_new(0:3,num,1:N_part)=rijpc_ep_old(0:3,num,1:N_part)
				END IF
				IF (flag_shadow) THEN
					rij_ese1_new(0:3,num,1:N_part)=rij_ese1_old(0:3,num,1:N_part)
					rij_ese2_new(0:3,num,1:N_part)=rij_ese2_old(0:3,num,1:N_part)
					IF ( iniz_pc ) THEN
						rijpc_ese1_new(0:3,num,1:N_part)=rijpc_ese1_old(0:3,num,1:N_part)
						rijpc_ese2_new(0:3,num,1:N_part)=rijpc_ese2_old(0:3,num,1:N_part)
					END IF
				END IF
				!se1
				se1_new(1:3,num)=se1_old(1:3,num)
				rij_se1_new(0:3,1:N_part,num)=rij_se1_old(0:3,1:N_part,num)
				rij_se1_new(0:3,num,1:N_part)=rij_se1_old(0:3,num,1:N_part)
				IF (flag_traccia_coppie_mol_ss) dist_mol_ss1_new(index_mol_num)=dist_mol_ss1_old(index_mol_num)
				IF (iniz_pc) THEN
					rijpc_se1_new(0:3,1:N_part,num)=rijpc_se1_old(0:3,1:N_part,num)
					rijpc_se1_new(0:3,num,1:N_part)=rijpc_se1_old(0:3,num,1:N_part)
				END IF
				rij_ese1_new(0:3,1:N_part,num)=rij_ese1_old(0:3,1:N_part,num)
				IF ( iniz_pc ) THEN
					rijpc_ese1_new(0:3,1:N_part,num)=rijpc_ese1_old(0:3,1:N_part,num)
				END IF
				rij_sesp1_new(0:3,num,1:N_part)=rij_sesp1_old(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_sesp1_new(0:3,num,1:N_part)=rijpc_sesp1_old(0:3,num,1:N_part)
				END IF
				!se2
				se2_new(1:3,num)=se2_old(1:3,num)
				rij_se2_new(0:3,1:N_part,num)=rij_se2_old(0:3,1:N_part,num)
				rij_se2_new(0:3,num,1:N_part)=rij_se2_old(0:3,num,1:N_part)
				IF (flag_traccia_coppie_mol_ss) dist_mol_ss2_new(index_mol_num)=dist_mol_ss2_old(index_mol_num)
				IF (iniz_pc) THEN
					rijpc_se2_new(0:3,1:N_part,num)=rijpc_se2_old(0:3,1:N_part,num)
					rijpc_se2_new(0:3,num,1:N_part)=rijpc_se2_old(0:3,num,1:N_part)
				END IF
				rij_ese2_new(0:3,1:N_part,num)=rij_ese2_old(0:3,1:N_part,num)
				IF ( iniz_pc ) THEN
					rijpc_ese2_new(0:3,1:N_part,num)=rijpc_ese2_old(0:3,1:N_part,num)
				END IF
				rij_sesp2_new(0:3,num,1:N_part)=rij_sesp2_old(0:3,num,1:N_part)
				IF (iniz_pc) THEN
					rijpc_sesp2_new(0:3,num,1:N_part)=rijpc_sesp2_old(0:3,num,1:N_part)
				END IF
			END SELECT
		END IF
				
	END SUBROUTINE riporta_indietro_walker
!-----------------------------------------------------------------------

	SUBROUTINE salva_posizione_walkers(nome_file)
		USE dati_fisici
		USE dati_mc
		USE generic_tools
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: nome_file
		CHARACTER(LEN=4) :: istring
		INTEGER :: i, j
		REAL (KIND=8), ALLOCATABLE :: x(:,:)
				
		IF (flag_normalizza_pos) ALLOCATE(x(1:3,1:N_part))
		WRITE (istring, '(I4.4)'), mpi_myrank
		
		IF (flag_normalizza_pos) THEN
			x(1,1:N_part)=re_old(1,1:N_part)/L(1)
			x(2,1:N_part)=re_old(2,1:N_part)/L(2)
			x(3,1:N_part)=re_old(3,1:N_part)/L(3)
			DO i = 1, N_part, 1
				DO j = 1, 3, 1
					IF (x(j,i)<0.d0) x(j,i)=x(j,i)+1.d0
				END DO
			END DO
			CALL salva_posizioni_su_file(x,N_part,nome_file//'_e_'//istring)
		ELSE
			CALL salva_posizioni_su_file(re_old,N_part,nome_file//'_e_'//istring)
		END IF
		IF (flag_shadow) THEN
			IF (flag_normalizza_pos) THEN
				x(1,1:N_part)=se1_old(1,1:N_part)/L(1)
				x(2,1:N_part)=se1_old(2,1:N_part)/L(2)
				x(3,1:N_part)=se1_old(3,1:N_part)/L(3)
				DO i = 1, N_part, 1
					DO j = 1, 3, 1
						IF (x(j,i)<0.d0) x(j,i)=x(j,i)+1.d0
					END DO
				END DO
				CALL salva_posizioni_su_file(x,N_part,nome_file//'_se1_'//istring)
			ELSE
				CALL salva_posizioni_su_file(se1_old,N_part,nome_file//'_se1_'//istring)
			END IF
			IF (flag_normalizza_pos) THEN
				x(1,1:N_part)=se2_old(1,1:N_part)/L(1)
				x(2,1:N_part)=se2_old(2,1:N_part)/L(2)
				x(3,1:N_part)=se2_old(3,1:N_part)/L(3)
				DO i = 1, N_part, 1
					DO j = 1, 3, 1
						IF (x(j,i)<0.d0) x(j,i)=x(j,i)+1.d0
					END DO
				END DO
				CALL salva_posizioni_su_file(x,N_part,nome_file//'_se2_'//istring)
			ELSE
				CALL salva_posizioni_su_file(se2_old,N_part,nome_file//'_se2_'//istring)
			END IF
			IF (flag_normalizza_pos) THEN
				x(1,1:N_part)=sp1_old(1,1:N_part)/L(1)
				x(2,1:N_part)=sp1_old(2,1:N_part)/L(2)
				x(3,1:N_part)=sp1_old(3,1:N_part)/L(3)
				DO i = 1, N_part, 1
					DO j = 1, 3, 1
						IF (x(j,i)<0.d0) x(j,i)=x(j,i)+1.d0
					END DO
				END DO
				CALL salva_posizioni_su_file(x,N_part,nome_file//'_sp1_'//istring)
			ELSE
				CALL salva_posizioni_su_file(sp1_old,N_part,nome_file//'_sp1_'//istring)
			END IF
			IF (flag_normalizza_pos) THEN
				x(1,1:N_part)=sp2_old(1,1:N_part)/L(1)
				x(2,1:N_part)=sp2_old(2,1:N_part)/L(2)
				x(3,1:N_part)=sp2_old(3,1:N_part)/L(3)
				DO i = 1, N_part, 1
					DO j = 1, 3, 1
						IF (x(j,i)<0.d0) x(j,i)=x(j,i)+1.d0
					END DO
				END DO
				CALL salva_posizioni_su_file(x,N_part,nome_file//'_sp2_'//istring)
			ELSE
				CALL salva_posizioni_su_file(sp2_old,N_part,nome_file//'_sp2_'//istring)
			END IF
		END IF
		IF (flag_normalizza_pos) THEN
			x(1,1:N_part)=rp_old(1,1:N_part)/L(1)
			x(2,1:N_part)=rp_old(2,1:N_part)/L(2)
			x(3,1:N_part)=rp_old(3,1:N_part)/L(3)
			DO i = 1, N_part, 1
				DO j = 1, 3, 1
					IF (x(j,i)<0.d0) x(j,i)=x(j,i)+1.d0
				END DO
			END DO
			CALL salva_posizioni_su_file(x,N_part,nome_file//'_p_'//istring)
		ELSE
			CALL salva_posizioni_su_file(rp_old,N_part,nome_file//'_p_'//istring)
		END IF
		
		IF (flag_normalizza_pos) DEALLOCATE(x)
				
	END SUBROUTINE salva_posizione_walkers
!-----------------------------------------------------------------------
	
	SUBROUTINE setta_Rcrystal(vett)
		USE generic_tools
		USE dati_fisici
		USE dati_mc
		IMPLICIT NONE
		REAL (KIND=8), DIMENSION(:) :: vett(:)
		REAL (KIND=8) :: r_crystal_new(1:3,1:N_part)
		REAL (KIND=8), ALLOCATABLE :: rij_crystal(:,:,:)
		INTEGER :: i, j, cont
		REAL (KIND=8) :: net_dist
				
		cont=1
		DO i = 1, N_part, 1
			r_crystal_new(1:3,i)=vett(cont:cont+2)
			cont=cont+3
		END DO
		CALL applica_pbc(r_crystal_new,N_part,L)
		
		!!! OPT1
		!ALLOCATE(rij_crystal(0:3,1:N_part,1:N_part))
		!CALL valuta_distanza_ii(r_crystal,N_part,L,rij_crystal)
		!net_dist=0.d0
		!DO j = 1, N_part-1, 1
		!	DO i = j+1, N_part, 1
		!		net_dist=MIN(net_dist,rij_crystal(i,j))
		!	END DO
		!END DO
		!!! OPT2
		!net_dist=MINVAL(L)/(N_part**(1.d0/3.d0))
		!!! OPT3, uso il raggio di Bohr
		net_dist=1.d0
		
		!MINVAL(L)/(N_part**(1.d0/3.d0))
		CALL RANDOM_NUMBER(re_new)
		DO j = 1, N_part, 1
			DO i = 1, 3, 1
				DO WHILE (DABS(re_new(i,j)-0.5d0)<0.001d0)
					CALL RANDOM_NUMBER(re_new(i,j))
				END DO
			END DO
		END DO
		re_old=r_crystal_new+(re_new-0.5d0)*net_dist
		CALL applica_pbc(re_old,N_part,L)
		IF (flag_shadow) THEN
			CALL RANDOM_NUMBER(se1_old)
			se1_old=(se1_old-0.5d0)*0.01d0*net_dist
			se1_old=re_old+se1_old
			CALL RANDOM_NUMBER(se2_old)
			se2_old=(se2_old-0.5d0)*0.01d0*net_dist
			se2_old=re_old+se2_old
		END IF
		re_new=re_old
		CALL valuta_distanza_ii(re_old,N_part,L,rij_ee_old)
		rij_ee_new=rij_ee_old
		IF (flag_shadow) THEN
			se1_new=se1_old
			se2_new=se2_old
			CALL valuta_distanza_ii(se1_old,N_part,L,rij_se1_old)
			CALL valuta_distanza_ii(se2_old,N_part,L,rij_se2_old)
			rij_se1_new=rij_se1_old
			rij_se2_new=rij_se2_old
		END IF
		
		rp_old=r_crystal_new
		rp_new=rp_old

		CALL valuta_distanza_ii(rp_old,N_part,L,rij_pp_old)
		rij_pp_new=rij_pp_old
		
		CALL valuta_distanza_ij(re_old,rp_old,N_part,L,rij_ep_old)
		rij_ep_new=rij_ep_old
		
		IF ( flag_shadow ) THEN
			sp1_old=rp_old
			sp2_old=rp_old
			sp1_new=sp1_old
			sp2_new=sp2_old
			CALL valuta_distanza_ij(re_old,se1_old,N_part,L,rij_ese1_old)
			CALL valuta_distanza_ij(re_old,se2_old,N_part,L,rij_ese2_old)
			CALL valuta_distanza_ij(se1_old,sp1_old,N_part,L,rij_sesp1_old)
			CALL valuta_distanza_ij(se2_old,sp2_old,N_part,L,rij_sesp2_old)
			rij_ese1_new=rij_ese1_old
			rij_ese2_new=rij_ese2_old
			rij_sesp1_new=rij_sesp1_old
			rij_sesp2_new=rij_sesp2_old
		END IF
		
		IF (iniz_pc) THEN
			IF (flag_elettroni) THEN
				CALL valuta_distanza_pc_ii(rij_ee_old,N_part,L,rijpc_ee_old)
				rijpc_ee_new=rijpc_ee_old
				IF (flag_shadow) THEN
					CALL valuta_distanza_pc_ii(rij_se1_old,N_part,L,rijpc_se1_old)
					rijpc_se1_new=rijpc_se1_old
					CALL valuta_distanza_pc_ii(rij_se2_old,N_part,L,rijpc_se2_old)
					rijpc_se2_new=rijpc_se2_old
					CALL valuta_distanza_pc_ij(rij_sesp1_old,N_part,L,rijpc_sesp1_old)
					rijpc_sesp1_new=rijpc_sesp1_old
					CALL valuta_distanza_pc_ij(rij_sesp2_old,N_part,L,rijpc_sesp2_old)
					rijpc_sesp2_new=rijpc_sesp2_old
				END IF
			END IF
			IF (flag_protoni) THEN
				CALL valuta_distanza_pc_ii(rij_pp_old,N_part,L,rijpc_pp_old)
				rijpc_pp_new=rijpc_pp_old
			END IF
			IF (flag_protoni .OR. flag_elettroni) THEN
				CALL valuta_distanza_pc_ij(rij_ep_old,N_part,L,rijpc_ep_old)
				rijpc_ep_new=rijpc_ep_old
			END IF
		END IF
						
	END SUBROUTINE setta_Rcrystal
!-----------------------------------------------------------------------

	SUBROUTINE stampa_file_Rp(nome_file)
		USE dati_fisici
		IMPLICIT NONE
		CHARACTER(LEN=*) :: nome_file
		INTEGER :: i, ios
		
		IF (.NOT. iniz_walkers) STOP 'Prima avresti dovuto inizializzare &
		  [ module_walkers.f90 > stampa_file_Rp ]'
		
		OPEN(UNIT=37, FILE=nome_file, STATUS='UNKNOWN', IOSTAT=ios)
		IF ( ios /= 0 ) STOP 'Errore nell aprire il file per salvare Rp &
		  [ module_walkers.f90 > stampa_file_Rp ]'
		WRITE(UNIT=37, FMT=*), L(1:3)
		DO i = 1, N_part, 1
			WRITE(UNIT=37, FMT=*), rp_old(1:3,i)
		END DO
		CLOSE (37)
		
	END SUBROUTINE stampa_file_Rp
!-----------------------------------------------------------------------

	SUBROUTINE chiudi_walkers()
		USE dati_mc
		USE dati_fisici
		IMPLICIT NONE
		IF (.NOT. iniz_walkers) STOP 'Prima di chiudere avresti dovuto inizializzare &
		  [ module_walkers.f90 > chiudi_walkers ]'
		
		DEALLOCATE(re_new,re_old,se1_new,se1_old,se2_new,se2_old)
		DEALLOCATE(rij_ee_new,rij_ee_old,rij_se1_new,rij_se1_old,rij_se2_new,rij_se2_old)
		DEALLOCATE(rij_ep_new,rij_ep_old,rij_ese1_new,rij_ese1_old,rij_ese2_new,rij_ese2_old)
		DEALLOCATE(rp_new,rp_old,sp1_new,sp1_old,sp2_new,sp2_old)
		DEALLOCATE(rij_pp_new,rij_pp_old)
		DEALLOCATE(rij_sesp1_new,rij_sesp1_old,rij_sesp2_new,rij_sesp2_old)
		IF (iniz_pc) THEN
			IF (flag_elettroni) THEN
				DEALLOCATE(rijpc_ee_new,rijpc_ee_old)
				IF (flag_shadow) THEN
					DEALLOCATE(rijpc_se1_new,rijpc_se1_old,rijpc_se2_new,rijpc_se2_old)
					DEALLOCATE(rijpc_sesp1_new,rijpc_sesp1_old,rijpc_sesp2_new,rijpc_sesp2_old)
					DEALLOCATE(rijpc_ese1_new,rijpc_ese1_old,rijpc_ese2_new,rijpc_ese2_old)
				END IF
			END IF
			IF (flag_protoni) THEN
				DEALLOCATE(rijpc_pp_new,rijpc_pp_old)
			END IF
			IF (flag_protoni .OR. flag_elettroni) THEN
				DEALLOCATE(rijpc_ep_new,rijpc_ep_old)
			END IF
			iniz_pc=.FALSE.
		END IF
		IF ( flag_shadow .AND. flag_traccia_coppie_mol_ss ) THEN
			DEALLOCATE(index_mol_ss1,dist_mol_ss1_new,dist_mol_ss1_old)
			DEALLOCATE(index_mol_ss2,dist_mol_ss2_new,dist_mol_ss2_old)
		END IF
		
		iniz_walkers=.FALSE.
	END SUBROUTINE chiudi_walkers
	
END MODULE walkers
