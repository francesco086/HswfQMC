MODULE momenta
	IMPLICIT NONE
	LOGICAL, PROTECTED, SAVE :: iniz_momenta=.FALSE.
	LOGICAL, PROTECTED, SAVE :: fqt_called
	LOGICAL, PROTECTED :: save_momenta=.FALSE.
	INTEGER, PRIVATE, SAVE :: magic_number(1:86)=(/ 1,7,19,27,33,57,81,93,123,147,171,179,203,251,257,305,341,365,389,437, &
	  461,485,515,587,619,691,739,751,799,847,895,925,949,1021,1045,1141,1189,1213,1237,1309,1357,1365,1419, &
	  1503,1551,1575,1647,1743,1791,1839,1863,1935,2007,2103,2109,2205,2301,2325,2373,2469,2517,2553,2601, &
	  2721,2777,2801,2897,2945,2969,3071,3119,3191,3239,3287,3407,3431,3575,3695,3743,3791,3887,3911,3959, &
	  4067,4139,4169 /)
	INTEGER, PRIVATE, SAVE :: n_shell(1:86)=(/0,1,2,3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20,21,22,24,25,26,27,29,30,32,33, &
	  34,35,36,37,38,40,41,42,43,44,45,46,48,49,50,51,52,53,54,56,57,58,59,61,62,64,65,66,67,68,69,70,72,73, &
	  74,75,76,77,78,80,81,82,83,84,85,86,88,89,90,91,93,94,96,97,98,99,100 /)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: k_pw(:,:), k_fermi(:,:)
	REAL (KIND=8), ALLOCATABLE, PROTECTED, SAVE :: k_fermi_bigger(:,:)
	
	CONTAINS
	
	SUBROUTINE inizializza_momenta(H_N, num_k_ewald, L, mpi_myrank)
		USE generic_tools
		IMPLICIT NONE
		LOGICAL :: flag_file
		INTEGER, INTENT(IN) :: H_N, num_k_ewald, mpi_myrank
		INTEGER :: i
		REAL (KIND=8), INTENT(IN) :: L(1:3)
		
		ALLOCATE(k_pw(0:3,1:H_N),k_fermi(0:3,1:num_k_ewald))
		
		CALL fermi_quantization(H_N, L, k_pw)
		CALL fermi_quantization(num_k_ewald, L, k_fermi)
		
		IF (save_momenta) THEN
			IF (mpi_myrank==0) THEN
				CALL salva_posizioni_su_file(k_pw(1:3,1:H_N),H_N,'posizioni/k_pw')
			END IF
			INQUIRE(FILE='posizioni/k_pw_TABC.pos',EXIST=flag_file)
			IF (flag_file) CALL SYSTEM ('rm posizioni/k_pw_TABC.pos')
			CALL SYSTEM ('> posizioni/k_pw'//CHAR(mpi_myrank+65)//'.pos')
			CALL SYSTEM ('> posizioni/k_pw_TABC'//CHAR(mpi_myrank+65)//'.pos')
			CALL SYSTEM ('> posizioni/app'//CHAR(mpi_myrank+65)//'.pos')
		END IF
		
		fqt_called=.FALSE.
		iniz_momenta=.TRUE.
	END SUBROUTINE inizializza_momenta
!-----------------------------------------------------------------------
	SUBROUTINE applica_twist(H_N, L)
		USE generic_tools
		USE dati_mc
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: H_N
		REAL (KIND=8), INTENT(IN) :: L(1:3)
		IF (.NOT. iniz_momenta) STOP 'Prima di applicare un twist devi inizializzare le k &
		  [ module_momenta.f90 > applica_twist ]'
		CALL fermi_quantization_twist(H_N, L, k_pw)
		IF (save_momenta) THEN
			CALL salva_posizioni_su_file(k_pw(1:3,1:H_N),H_N,'posizioni/k_pw'//CHAR(mpi_myrank+65))
			CALL SYSTEM ('cat posizioni/k_pw'//CHAR(mpi_myrank+65)//'.pos posizioni/k_pw_TABC'//CHAR(mpi_myrank+65)//&
			   '.pos > posizioni/app'//CHAR(mpi_myrank+65)//'.pos')
			CALL SYSTEM ('cat posizioni/app'//CHAR(mpi_myrank+65)//'.pos > posizioni/k_pw_TABC'//CHAR(mpi_myrank+65)//'.pos')
		END IF
	END SUBROUTINE applica_twist
!-----------------------------------------------------------------------
	!calcola i primi N momenti di fermi quantizzati in una scatola di lato L
	SUBROUTINE fermi_quantization(N_part, L, k)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N_part
		REAL (KIND=8), INTENT(IN) :: L(1:3)
		INTEGER :: N, N_liv, res, p1, p2, p3, step_acc, part_num
		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
		REAL (KIND=8) :: k_min(1:3)
		REAL (KIND=8), INTENT(OUT) :: k(0:3, 1:N_part)

		!inizializzo il contatore sul numero di particelle
		part_num=0
		!determino il k minimo dato dalla quantizzazione
		k_min(1:3)=2.D0*PI/L(1:3)
		!k_min=1.D0
		!fisso il numero massimo del livello che voglio raggiungere (NB deve essere abbastanza grande)
		N_liv=10000

		loop1: DO  N = 0, N_liv, 1
			p1=0
			DO WHILE ( p1*p1 < N )
				p1=p1+1
				IF ((p1+1)*(p1+1)>N) EXIT
			END DO
			DO p1=p1, 0, -1
				res=N-p1*p1
				p2=0
				DO WHILE ( p2*p2 < res )
					p2=p2+1
					IF ((p2+1)*(p2+1)>res) EXIT
				END DO
				DO p2=p2, 0, -1
					res=N-p1*p1-p2*p2
					p3=0
					DO WHILE ( p3*p3 < res )
						p3=p3+1
						IF ((p3+1)*(p3+1)>res) EXIT
					END DO
					DO p3=p3, 0, -1
						res=N-p1*p1-p2*p2-p3*p3
						IF (res==0) THEN
							step_acc=1
						ELSE
							step_acc=0
						END IF
						IF (step_acc==1) THEN !qui dentro ho la combinazione giusta, solo attenzione ai segni
							part_num=part_num+1
							IF ( part_num>N_part ) THEN
								EXIT loop1
							END IF
							k(1:3,part_num)=(/ k_min(1)*p1,k_min(2)*p2,k_min(3)*p3 /)
							IF (p3/=0) THEN
								IF ( p2/=0 ) THEN
									IF ( p1/=0 ) THEN
										part_num=part_num+1
										IF ( part_num>N_part ) THEN
											EXIT loop1
										END IF
										k(1:3,part_num)=(/ -k_min(1)*p1,-k_min(2)*p2,-k_min(3)*p3 /)
									END IF
									part_num=part_num+1
									IF ( part_num>N_part ) THEN
										EXIT loop1
									END IF
									k(1:3,part_num)=(/ k_min(1)*p1,-k_min(2)*p2,-k_min(3)*p3 /)
								END IF
								IF (p1/=0) THEN
									part_num=part_num+1
									IF ( part_num>N_part ) THEN
										EXIT loop1
									END IF
									k(1:3,part_num)=(/ -k_min(1)*p1,k_min(2)*p2,-k_min(3)*p3 /)
								END IF
								part_num=part_num+1
								IF ( part_num>N_part ) THEN
									EXIT loop1
								END IF
								k(1:3,part_num)=(/ k_min(1)*p1,k_min(2)*p2,-k_min(3)*p3 /)
							END IF
							IF (p2/=0) THEN
								IF ( p1/=0 ) THEN
									part_num=part_num+1
									IF ( part_num>N_part ) THEN
										EXIT loop1
									END IF
									k(1:3,part_num)=(/ -k_min(1)*p1,-k_min(2)*p2,k_min(3)*p3 /)
								END IF
								part_num=part_num+1
								IF ( part_num>N_part ) THEN
									EXIT loop1
								END IF
								k(1:3,part_num)=(/ k_min(1)*p1,-k_min(2)*p2,k_min(3)*p3 /)
							END IF
							IF (p1/=0) THEN
								part_num=part_num+1
								IF ( part_num>N_part ) THEN
									EXIT loop1
								END IF
								k(1:3,part_num)=(/ -k_min(1)*p1,k_min(2)*p2,k_min(3)*p3 /)
							END IF
						END IF
					END DO
				END DO
			END DO
		END DO loop1
		DO part_num = 1, N_part, 1
			k(0,part_num)=DOT_PRODUCT(k(1:3,part_num),k(1:3,part_num))
		END DO
	END SUBROUTINE fermi_quantization
!-----------------------------------------------------------------------
	!calcola i primi N momenti di fermi in una scatola di lato L, applicando un twist
	SUBROUTINE fermi_quantization_twist(N_part,L,k_fermi_twist)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N_part
		REAL (KIND=8), INTENT(IN) :: L(1:3)
		INTEGER, SAVE :: index, inc2, index_bigger, N_part_bigger, i
		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
		REAL :: random
		REAL (KIND=8) :: vect(0:3), twist(1:3)
		REAL (KIND=8), ALLOCATABLE :: k(:,:)
		REAL (KIND=8), INTENT(OUT) :: k_fermi_twist(0:3,1:N_part)
		
		IF ( .NOT. fqt_called ) THEN
			fqt_called=.TRUE.
			index=1
			DO WHILE (magic_number(index)<N_part)
				index=index+1
				IF (index>86) STOP 'N_part troppo grande [ module_momenta.f90 > fermi_quantization_twist ]'
			END DO
			inc2=FLOOR(0.75+SQRT(3.)*SQRT(REAL(n_shell(index))))
			index_bigger=index
			DO WHILE (n_shell(index_bigger)<n_shell(index)+inc2)
				index_bigger=index_bigger+1
			END DO
			N_part_bigger=magic_number(index_bigger)
			ALLOCATE(k_fermi_bigger(0:3,1:N_part_bigger))
			CALL fermi_quantization(N_part_bigger,L,k_fermi_bigger)
		END IF

		ALLOCATE(k(0:3,1:N_part_bigger))
		CALL RANDOM_NUMBER(twist)
		twist(1:3)=(twist(1:3)-0.5d0)*2.d0*PI/L(1:3)
		!traslo k
		DO i = 1, N_part_bigger, 1
			k(1:3,i)=k_fermi_bigger(1:3,i)+twist(1:3)
			k(0,i)=DOT_PRODUCT(k(1:3,i),k(1:3,i))
		END DO
		!disordino in modo random i k oltre N_part con norma uguale
		i=N_part+1
		DO WHILE (i<N_part_bigger)
			IF (k(0,i)==k(0,i+1)) THEN
				CALL RANDOM_NUMBER(random)
				IF (random>0.5) THEN
					vect(0:3)=k(0:3,i)
					k(0:3,i)=k(0:3,i+1)
					k(0:3,i+1)=vect(0:3)
					IF (i>N_part+1) i=i-2
				END IF
			END IF
			i=i+1
		END DO
		!ordino in modo crescente le k
		i=1
		DO WHILE (i<N_part_bigger)
			IF (k(0,i)>k(0,i+1)) THEN
				vect(0:3)=k(0:3,i)
				k(0:3,i)=k(0:3,i+1)
				k(0:3,i+1)=vect(0:3)
				IF (i>1) i=i-2
			END IF
			i=i+1
		END DO
		k_fermi_twist(0:3,1:N_part)=k(0:3,1:N_part)
		DEALLOCATE(k)
	END SUBROUTINE fermi_quantization_twist
!-----------------------------------------------------------------------
	SUBROUTINE chiudi_momenta()
		USE dati_mc
		IMPLICIT NONE
		IF (.NOT. iniz_momenta) STOP 'Prima di chiudere devi inizializzare i momenti &
		  [ module_momenta.f90 > chiudi_momenta ]'
		
		IF (save_momenta) THEN
			IF (mpi_myrank==0) CALL SYSTEM ('cat posizioni/k_pw_TABC*.pos > posizioni/k_pw_TABC.pos')
			CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
			CALL SYSTEM ('rm posizioni/k_pw'//CHAR(mpi_myrank+65)//'.pos')
			CALL SYSTEM ('rm posizioni/k_pw_TABC'//CHAR(mpi_myrank+65)//'.pos')
			CALL SYSTEM ('rm posizioni/app'//CHAR(mpi_myrank+65)//'.pos')
		END IF
		
		DEALLOCATE(k_pw,k_fermi)
		IF (fqt_called) DEALLOCATE(k_fermi_bigger)
		iniz_momenta=.FALSE.
		fqt_called=.FALSE.
	END SUBROUTINE chiudi_momenta

END MODULE momenta
