MODULE dati_fisici
	USE lattice
	IMPLICIT NONE
	LOGICAL, PROTECTED, SAVE :: iniz_dati_fisici=.FALSE., flag_molecular=.FALSE., flag_2D=.FALSE., flag_tilted=.FALSE.
	CHARACTER(LEN=5), PROTECTED, SAVE :: crystal_cell
	INTEGER, PROTECTED, SAVE :: N_part, H_N_part, N_cell_side
	REAL (KIND=8), PARAMETER, PRIVATE :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
	REAL (KIND=8), PARAMETER :: MASS_e=0.0005485899094d0, MASS_p=1.00727646677d0   !in uma
	REAL (KIND=8), PROTECTED, SAVE :: r_s, L(1:3), H_L(1:3), L_cov_bond, L_mat(1:3,1:3), L_mati(1:3,1:3)
	REAL (KIND=8), PROTECTED, SAVE :: hbar, K_coulomb, strecthing_cov_bond
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: r_crystal(:,:)
	
	CONTAINS
	
	SUBROUTINE inizializza_dati_fisici(file_reticolo_opt)
		USE generic_tools
		IMPLICIT NONE
		CHARACTER(LEN=*), OPTIONAL :: file_reticolo_opt     !per ottimizzare le posizioni protoniche con SR
		INTEGER :: i, j, i1, i_seed
		INTEGER, ALLOCATABLE :: seed(:), seed_provv(:)
		CHARACTER(LEN=100) :: file_reticolo
		REAL (KIND=8) :: vect(1:3), sigma_w, eta_w(1:1,1:1), L_w, dist(0:3)
		REAL (KIND=8), ALLOCATABLE :: app(:,:)
		
		NAMELIST /dati_fisici/ r_s, crystal_cell, flag_2D, flag_tilted, file_reticolo, flag_molecular, &
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
		ELSE IF ( crystal_cell=='hring' ) THEN
			N_part=6
		ELSE IF ( crystal_cell=='grp__' ) THEN
			N_part=4*(N_cell_side**2)
      flag_2D=.TRUE.
    ELSE IF ( crystal_cell=='quadr' ) THEN
      N_part=N_cell_side**2
      flag_2D=.TRUE.
    ELSE IF ( crystal_cell=='trian' ) THEN
      N_part=2*(N_cell_side**2)
      flag_2D=.TRUE.
    ELSE IF ( crystal_cell=='tilt_') THEN
      N_part=N_cell_side
      flag_tilted=.TRUE.
    ELSE
			STOP "DATI_FISICI: scegli un reticolo accettabile"
		END IF
		H_N_part=N_part/2

		IF (PRESENT(file_reticolo_opt)) THEN
			crystal_cell='datex'
			file_reticolo=file_reticolo_opt
		END IF

      !3D structures
		L(1:3)=r_s*(4.d0*PI*N_part/3.d0)**(1.d0/3.d0)            !in bohr units
		IF ((crystal_cell=='hcp__').OR.(crystal_cell=='mhcpo').OR.(crystal_cell=='hcp_w')) THEN
			L(1)=L(1)*(2.d0**(1.d0/6.d0))
			L(2)=L(2)*(((3.d0**3.d0)/(2.d0**5.d0))**(1.d0/6.d0))
			L(3)=L(3)*((2.d0**(2.d0/3.d0))/(3**(1.d0/2.d0)))   *(1.58d0/DSQRT(8.d0/3.d0))
		END IF

      ! 2D structures
		IF ( crystal_cell=='grp__' ) THEN
			L(1:2)=r_s*DSQRT(PI*N_part)
			L(1)=L(1)*3.d0/DSQRT(DSQRT(27.d0))
			L(2)=L(2)*DSQRT(3.d0)/DSQRT(DSQRT(27.d0))
         L(3)=100.d0 	!L(3)=L(1)
		END IF
      IF ( crystal_cell=='quadr' ) THEN
         L(1:2)=DSQRT(PI*REAL(N_part,8))*r_s
         L(3)=100.d0 !very large number, so that the layer can be considered isolated
      END IF
      IF ( crystal_cell=='trian' ) THEN
         L(1)=DSQRT(PI*REAL(N_part,8)/DSQRT(2.d0))*r_s
         L(2)=DSQRT(2.d0)*L(1)
         L(3)=100.d0 !very large number, so that the layer can be considered isolated
      END IF
		H_L=0.5d0*L

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

		ELSE IF ( crystal_cell=='hring' ) THEN

      ! Spin is alternating around the hydrogen ring (up: 1,2,3, dw: 4,5,6)
			r_crystal(1:3,1)=(/ L_cov_bond		,	0.d0						,	0.d0/)
			r_crystal(1:3,4)=(/ L_cov_bond*0.5d0,	sqrt(3.d0)*L_cov_bond*0.5d0	,	0.d0/)
			r_crystal(1:3,2)=(/-L_cov_bond*0.5d0,	sqrt(3.d0)*L_cov_bond*0.5d0	,	0.d0/)
			r_crystal(1:3,5)=(/-L_cov_bond		,	0.d0						,	0.d0/)
			r_crystal(1:3,3)=(/-L_cov_bond*0.5d0,	-sqrt(3.d0)*L_cov_bond*0.5d0,	0.d0/)
			r_crystal(1:3,6)=(/ L_cov_bond*0.5d0,	-sqrt(3.d0)*L_cov_bond*0.5d0,	0.d0/)

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
				READ (2, *) r_crystal(1:3,i)
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
			READ (2, *) L(1:3)
			H_L=0.5d0*L
			DO i = 1, N_part, 1
				READ (2, *) r_crystal(1:3,i)
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
						STOP 'distanze giá minori che per la fase molecolar [ module_funzione_onda.f90 > inizializza_funzione_onda ]'
					END IF
				END DO
			END IF
      ELSE IF ( crystal_cell=='quadr' ) THEN
         IF (.NOT. flag_molecular) THEN
            i1=0
            DO j = 1, N_cell_side, 1
            DO i = 1, N_cell_side, 1
               i1=i1+1
               IF (MOD(i1,2)==1) THEN
                  r_crystal(1:3,i1/2+1)=(/ REAL(i,8), REAL(j,8), 0.d0 /)*L(1:3)/REAL(N_cell_side,8)
               ELSE
                  r_crystal(1:3,i1/2+H_N_part)=(/ REAL(i,8), REAL(j,8), 0.d0 /)*L(1:3)/REAL(N_cell_side,8)
               END IF
            END DO
            END DO
         ELSE
            STOP "crystal_cell='quad' con flag_molecular non ancora implementato"
         END IF
         CALL applica_pbc(r_crystal,N_part,L, L_mat,L_mati,flag_tilted)
      ELSE IF ( crystal_cell=='trian' ) THEN
         IF (.NOT. flag_molecular) THEN
            i1=1
            DO j = 1, N_cell_side, 1
            DO i = 1, N_cell_side, 1
                  r_crystal(1:3,i1)=(/ REAL(i,8), REAL(j,8), 0.d0 /)*L(1:3)/REAL(N_cell_side,8)
                  r_crystal(1:3,i1+H_N_part)=(  (/ REAL(i,8), REAL(j,8), 0.d0 /)*L(1:3) &
                                                + (/ 0.5d0*L(1), 0.5d0*DSQRT(2.d0)*L(1), 0.d0 /)  ) &
                                             /REAL(N_cell_side,8)
                  i1=i1+1
            END DO
            END DO
         ELSE
            STOP "crystal_cell='trian' con flag_molecular non ancora implementato"
         END IF
         CALL applica_pbc(r_crystal,N_part,L, L_mat,L_mati,flag_tilted)
      ELSE IF ( crystal_cell=='tilt_' ) THEN
         OPEN (UNIT=2, FILE=TRIM(file_reticolo), STATUS='OLD')
         READ (2, *) L(1:3)
         H_L=0.5d0*L
         DO i = 1, N_part, 1
            READ (2, *) r_crystal(1:3,i)
         END DO
         CLOSE (2)
         L_mat(:,:) = 0.d0
         L_mat(1,1) = L(1)
         L_mat(2,2) = L(2)
         L_mat(3,3) = L(3)
         L_mati(1,1) = 1.d0/L(1)
         L_mati(2,2) = 1.d0/L(2)
         L_mati(3,3) = 1.d0/L(3)
		ELSE
			STOP "scegli un reticolo accettabile"
		END IF
		
  CALL applica_pbc(r_crystal,N_part,L, L_mat,L_mati,flag_tilted)
		
		IF (MOD(N_part,2)/=0) STOP 'Stai lavorando con un numero di particelle non pari!!! [ module_funzione_onda.f90 > inizializza_funzione_onda ]'
		
		DEALLOCATE(app)
		
		iniz_dati_fisici=.TRUE.
	END SUBROUTINE inizializza_dati_fisici
!-----------------------------------------------------------------------
   SUBROUTINE setta_L(Lnew)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: Lnew(1:3)
      
      L=Lnew
      H_L=L*0.5d0

   END SUBROUTINE setta_L
!-----------------------------------------------------------------------
	SUBROUTINE chiudi_dati_fisici()
		IMPLICIT NONE
		IF (.NOT. iniz_dati_fisici) STOP 'Prima di chiudere avresti dovuto inizializzare [ module_dati.f90 > chiudi_dati_fisici ]'
		DEALLOCATE(r_crystal)
		iniz_dati_fisici=.FALSE.
	END SUBROUTINE chiudi_dati_fisici
	
END MODULE dati_fisici



