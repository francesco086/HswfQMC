MODULE dnfH
	USE dati_fisici
   USE fermi_k
	IMPLICIT NONE
	REAL (KIND=8), PARAMETER, PRIVATE :: cutoff_quick_pw=0.0001d0
	LOGICAL, PRIVATE, SAVE :: flag_inizializza=.FALSE.
	INTEGER, PRIVATE, SAVE :: magic_number(1:86)=(/ 1,7,19,27,33,57,81,93,123,147,171,179,203,251,257,305,341,365,389,437, &
	  461,485,515,587,619,691,739,751,799,847,895,925,949,1021,1045,1141,1189,1213,1237,1309,1357,1365,1419, &
	  1503,1551,1575,1647,1743,1791,1839,1863,1935,2007,2103,2109,2205,2301,2325,2373,2469,2517,2553,2601, &
	  2721,2777,2801,2897,2945,2969,3071,3119,3191,3239,3287,3407,3431,3575,3695,3743,3791,3887,3911,3959, &
	  4067,4139,4169 /)
	INTEGER, PRIVATE, SAVE :: n_shell(1:86)=(/0,1,2,3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20,21,22,24,25,26,27,29,30,32,33, &
	  34,35,36,37,38,40,41,42,43,44,45,46,48,49,50,51,52,53,54,56,57,58,59,61,62,64,65,66,67,68,69,70,72,73, &
	  74,75,76,77,78,80,81,82,83,84,85,86,88,89,90,91,93,94,96,97,98,99,100 /)
	INTEGER, SAVE :: N_M_Hartree
	INTEGER, ALLOCATABLE, SAVE :: n_fpw_Hartee(:), i_fpw_Hartree(:,:)
	REAL (KIND=8) :: k_cutoff
	REAL (KIND=8), ALLOCATABLE, SAVE :: k_Hartree(:,:), autovalori_Hartree(:)
   TYPE(KWaVe) :: kwv_Hartree
	COMPLEX (KIND=8), ALLOCATABLE, SAVE :: M_Hartree(:,:)
		
	CONTAINS
	
	SUBROUTINE inizializza_dnfH(k2_f_factor)
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: k2_f_factor
		INTEGER :: i1, i2, i
		
		IF (.NOT. iniz_dati_fisici) STOP 'Devi prima inizializzare i dati fisici &
		  [ module_dnfH.f90 > inizializza_dnfH ]' 
				
		i1=2   !indice per il numero di particelle esistenti
		DO WHILE (magic_number(i1)<H_N_part)
			i1=i1+1
		END DO
		i2=1   !i2 indice per il numero di onde piane da considerare, per prendere anche gli stati eccitati
		DO WHILE (n_shell(i2)<k2_f_factor*n_shell(i1))
			i2=i2+1
		END DO
		
		N_M_Hartree=magic_number(i2)
		ALLOCATE(M_Hartree(1:N_M_Hartree,1:N_M_Hartree),k_Hartree(0:3,1:N_M_Hartree),n_fpw_Hartee(1:N_M_Hartree))
		ALLOCATE(i_fpw_Hartree(1:N_M_Hartree,1:N_M_Hartree),autovalori_Hartree(1:N_M_Hartree))

      !!!Made obsolete from the KWaVe class
      !!!CALL fermi_quantization(N_M_Hartree,k_Hartree)
      CALL kwv_Hartree%initializeKWaVe(N_M_Hartree,3)
      IF (flag_2D) THEN
         CALL kwv_Hartree%buildBoxK((/ L(1), L(2), MIN(L(1),L(2)) /))
      ELSE 
         CALL kwv_Hartree%buildBoxK(L)
      END IF
      k_Hartree=kwv_Hartree%k
		
		k_cutoff=k_Hartree(0,H_N_part)*k2_f_factor
				
		flag_inizializza=.TRUE.
				
	END SUBROUTINE inizializza_dnfH
!-----------------------------------------------------------------------
	SUBROUTINE twist_k_dnfH(k2_f_factor)
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: k2_f_factor
		INTEGER :: i1, i2, new_N_M_Hartree
		REAL (KIND=8), ALLOCATABLE :: k_dummy(:,:)
		
		
		!!Se si vuole tenere costante il cutoff
		
		!i1=2       !indice per il numero di particelle esistenti
		!DO WHILE (magic_number(i1)<H_N_part)
		!	i1=i1+1
		!END DO
		!i2=1       !i2 indice per il numero di onde piane da considerare, per prendere anche i primi stati eccitati che rientrano nel cutoff
		!DO WHILE (n_shell(i2)<k2_f_factor*n_shell(i1))
		!	i2=i2+1
		!END DO
		!i2=i2+1      ! in questo modo ne prendo un po' di piú
		!N_M_Hartree=magic_number(i2)                  !cambio il valore di N_M_Hartree
		!DEALLOCATE(k_Hartree)
		!ALLOCATE(k_Hartree(0:3,1:N_M_Hartree))
		!CALL fermi_quantization_twist(N_M_Hartree,k_Hartree)
		!DEALLOCATE(M_Hartree,n_fpw_Hartee)
		!DEALLOCATE(i_fpw_Hartree,autovalori_Hartree)
		!ALLOCATE(M_Hartree(1:N_M_Hartree,1:N_M_Hartree),n_fpw_Hartee(1:N_M_Hartree))
		!ALLOCATE(i_fpw_Hartree(1:N_M_Hartree,1:N_M_Hartree),autovalori_Hartree(1:N_M_Hartree))
		
		
		!!Se si vuole tenere costante il numero di k considerati
		
      !!!Made obsolete from the KWaVe class
		!!!CALL fermi_quantization_twist(N_M_Hartree,k_Hartree)
      CALL kwv_Hartree%twistK()
      k_Hartree=kwv_Hartree%k
				
	END SUBROUTINE twist_k_dnfH
!-----------------------------------------------------------------------
	SUBROUTINE costruisci_matrice_Hartree(hamiltoniana,c_eff)
		USE walkers
		IMPLICIT NONE
		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
		CHARACTER(LEN=3), INTENT(IN) :: hamiltoniana
		COMPLEX (KIND=8), OPTIONAL :: c_eff(1:N_M_Hartree,1:N_M_Hartree)
		INTEGER :: i, j, i_p!, contatore_print=0, ios
		REAL (KIND=8) :: k_app(1:3), eta(1:N_M_Hartree,1:N_M_Hartree)
		COMPLEX (KIND=8) :: proton_fact, zzz
		
		IF (.NOT. flag_inizializza) STOP 'Devi prima inizializzare &
		  [ module_dnfH.f90 > trova_soluzioni_dnfH ]'
				
		SELECT CASE(hamiltoniana)
		CASE('fre')
			M_Hartree=(0.d0,0.d0)
			DO i = 1, N_M_Hartree, 1
				M_Hartree(i,i)=hbar*hbar*k_Hartree(0,i)/(2.d0*MASS_e)
			END DO
		CASE('prf')    !protonic field
			M_Hartree=(0.d0,0.d0)
			DO j = 1, N_M_Hartree, 1
				DO i = 1, N_M_Hartree, 1
					IF (i/=j) THEN
						k_app(1:3)=k_Hartree(1:3,i)-k_Hartree(1:3,j)
						M_Hartree(i,j)=-4.d0*PI*K_coulomb/(DOT_PRODUCT(k_app(1:3),k_app(1:3)))
						proton_fact=(0.d0,0.d0)
						DO i_p = 1, N_part, 1
							proton_fact=proton_fact+CDEXP(-(0.d0,1.d0)*DOT_PRODUCT(k_app(1:3),rp_old(1:3,i_p)))
						END DO
						M_Hartree(i,j)=M_Hartree(i,j)*proton_fact/(L(1)*L(2)*L(3))
					ELSE
						M_Hartree(i,i)=hbar*hbar*k_Hartree(0,i)/(2.d0*MASS_e)
					END IF
				END DO
			END DO
		END SELECT
				
		IF (PRESENT(c_eff)) THEN
			DO j = 1, N_M_Hartree, 1
				DO i = 1, j-1, 1
					c_eff(j,i)=DCONJG(c_eff(i,j))
				END DO
			END DO
			DO i = 1, N_M_Hartree, 1
				zzz=(0.d0,0.d0)
				zzz=REAL(c_eff(i,i),8)
				c_eff(i,i)=zzz
			END DO
			M_Hartree=M_Hartree+c_eff
		END IF
		
		!contatore_print=contatore_print+1
		!OPEN(UNIT=20+contatore_print, STATUS='UNKNOWN', IOSTAT=ios)
		!IF ( ios /= 0 ) STOP "Error opening file name"
		!DO j = 1, 3, 1
		!	WRITE(UNIT=20+contatore_print, FMT=*), M_Hartree(1:3,j)
		!END DO
		!CLOSE (20+contatore_print)
						
	END SUBROUTINE costruisci_matrice_Hartree
!-----------------------------------------------------------------------
	SUBROUTINE trova_soluzioni_dnfH()
		IMPLICIT NONE
		INTEGER :: info, i, j, cont
		REAL (KIND=8) :: rwork(1:3*N_M_Hartree-1)
		COMPLEX (KIND=8) :: work(1:3*N_M_Hartree-1)
		
		IF (.NOT. flag_inizializza) STOP 'Devi prima inizializzare &
		  [ module_dnfH.f90 > trova_soluzioni_dnfH ]'
				
		CALL ZHEEV('V','U',N_M_Hartree,M_Hartree,N_M_Hartree,autovalori_Hartree,work,3*N_M_Hartree-1,rwork,info)
		IF (info/=0) THEN
			PRINT * , 'ERRORE ', info
			STOP 'Errore nel calcolo degli autovalori. &
			  [module_dnfH.f90 > trova_soluzioni_dnfH ]'
		END IF
						
		DO i = 1, N_M_Hartree, 1
			cont=0
			DO j = 1, N_M_Hartree, 1
				IF (DSQRT(REAL(M_Hartree(j,i)*DCONJG(M_Hartree(j,i)),8))>cutoff_quick_pw) THEN
					cont=cont+1
					i_fpw_Hartree(cont,i)=j
				END IF
			END DO
			n_fpw_Hartee(i)=cont                      !numero di pw effettive
		END DO
		        		 
	END SUBROUTINE trova_soluzioni_dnfH
!!-----------------------------------------------------------------------
!	SUBROUTINE fermi_quantization(N_fq,k)
!		IMPLICIT NONE
!		INTEGER, INTENT(IN) :: N_fq
!		INTEGER :: N, N_liv, res, p1, p2, p3, step_acc, part_num
!		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
!		REAL (KIND=8) :: k_min(1:3)
!		REAL (KIND=8), INTENT(OUT) :: k(0:3, 1:N_fq)
!
!		!inizializzo il contatore sul numero di particelle
!		part_num=0
!		!determino il k minimo dato dalla quantizzazione
!		k_min(1:3)=2.D0*PI/L(1:3)
!		!k_min=1.D0
!		!fisso il numero massimo del livello che voglio raggiungere (NB deve essere abbastanza grande)
!		N_liv=10000
!
!		loop1: DO  N = 0, N_liv, 1
!			p1=0
!			DO WHILE ( p1*p1 < N )
!				p1=p1+1
!				IF ((p1+1)*(p1+1)>N) EXIT
!			END DO
!			DO p1=p1, 0, -1
!				res=N-p1*p1
!				p2=0
!				DO WHILE ( p2*p2 < res )
!					p2=p2+1
!					IF ((p2+1)*(p2+1)>res) EXIT
!				END DO
!				DO p2=p2, 0, -1
!					res=N-p1*p1-p2*p2
!					p3=0
!					DO WHILE ( p3*p3 < res )
!						p3=p3+1
!						IF ((p3+1)*(p3+1)>res) EXIT
!					END DO
!					DO p3=p3, 0, -1
!						res=N-p1*p1-p2*p2-p3*p3
!						IF (res==0) THEN
!							step_acc=1
!						ELSE
!							step_acc=0
!						END IF
!						IF (step_acc==1) THEN !qui dentro ho la combinazione giusta, solo attenzione ai segni
!							part_num=part_num+1
!							IF ( part_num>N_fq ) THEN
!								EXIT loop1
!							END IF
!							k(1:3,part_num)=(/ k_min(1)*p1,k_min(2)*p2,k_min(3)*p3 /)
!							IF (p3/=0) THEN
!								IF ( p2/=0 ) THEN
!									IF ( p1/=0 ) THEN
!										part_num=part_num+1
!										IF ( part_num>N_fq ) THEN
!											EXIT loop1
!										END IF
!										k(1:3,part_num)=(/ -k_min(1)*p1,-k_min(2)*p2,-k_min(3)*p3 /)
!									END IF
!									part_num=part_num+1
!									IF ( part_num>N_fq ) THEN
!										EXIT loop1
!									END IF
!									k(1:3,part_num)=(/ k_min(1)*p1,-k_min(2)*p2,-k_min(3)*p3 /)
!								END IF
!								IF (p1/=0) THEN
!									part_num=part_num+1
!									IF ( part_num>N_fq ) THEN
!										EXIT loop1
!									END IF
!									k(1:3,part_num)=(/ -k_min(1)*p1,k_min(2)*p2,-k_min(3)*p3 /)
!								END IF
!								part_num=part_num+1
!								IF ( part_num>N_fq ) THEN
!									EXIT loop1
!								END IF
!								k(1:3,part_num)=(/ k_min(1)*p1,k_min(2)*p2,-k_min(3)*p3 /)
!							END IF
!							IF (p2/=0) THEN
!								IF ( p1/=0 ) THEN
!									part_num=part_num+1
!									IF ( part_num>N_fq ) THEN
!										EXIT loop1
!									END IF
!									k(1:3,part_num)=(/ -k_min(1)*p1,-k_min(2)*p2,k_min(3)*p3 /)
!								END IF
!								part_num=part_num+1
!								IF ( part_num>N_fq ) THEN
!									EXIT loop1
!								END IF
!								k(1:3,part_num)=(/ k_min(1)*p1,-k_min(2)*p2,k_min(3)*p3 /)
!							END IF
!							IF (p1/=0) THEN
!								part_num=part_num+1
!								IF ( part_num>N_fq ) THEN
!									EXIT loop1
!								END IF
!								k(1:3,part_num)=(/ -k_min(1)*p1,k_min(2)*p2,k_min(3)*p3 /)
!							END IF
!						END IF
!					END DO
!				END DO
!			END DO
!		END DO loop1
!		DO part_num = 1, N_fq, 1
!			k(0,part_num)=DOT_PRODUCT(k(1:3,part_num),k(1:3,part_num))
!		END DO
!	END SUBROUTINE fermi_quantization
!!!-----------------------------------------------------------------------
!	SUBROUTINE fermi_quantization_twist(N,k_fermi_twist)
!		IMPLICIT NONE
!		INTEGER, INTENT(IN) :: N
!		INTEGER, SAVE :: index, inc2, index_bigger, N_bigger, i
!		REAL (KIND=8), PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
!		REAL :: random
!		REAL (KIND=8) :: vect(0:3), twist(1:3)
!		REAL (KIND=8), ALLOCATABLE :: k(:,:), k_fermi_bigger(:,:)
!		REAL (KIND=8), INTENT(OUT) :: k_fermi_twist(0:3,1:N)
!		
!		index=1
!		DO WHILE (magic_number(index)<N)
!			index=index+1
!			IF (index>86) STOP 'N troppo grande [ module_momenta.f90 > fermi_quantization_twist ]'
!		END DO
!		inc2=FLOOR(0.75+SQRT(3.)*SQRT(REAL(n_shell(index))))
!		index_bigger=index
!		DO WHILE (n_shell(index_bigger)<n_shell(index)+inc2)
!			index_bigger=index_bigger+1
!		END DO
!		N_bigger=magic_number(index_bigger)
!		ALLOCATE(k_fermi_bigger(0:3,1:N_bigger))
!		CALL fermi_quantization(N_bigger, k_fermi_bigger)
!
!		ALLOCATE(k(0:3,1:N_bigger))
!		CALL RANDOM_NUMBER(twist)
!		twist(1:3)=(twist(1:3)-0.5d0)*2.d0*PI/L(1:3)
!		!traslo k
!		DO i = 1, N_bigger, 1
!			k(1:3,i)=k_fermi_bigger(1:3,i)+twist(1:3)
!			k(0,i)=DOT_PRODUCT(k(1:3,i),k(1:3,i))
!		END DO
!		!disordino in modo random i k oltre N con norma uguale
!		i=N+1
!		DO WHILE (i<N_bigger)
!			IF (k(0,i)==k(0,i+1)) THEN
!				CALL RANDOM_NUMBER(random)
!				IF (random>0.5) THEN
!					vect(0:3)=k(0:3,i)
!					k(0:3,i)=k(0:3,i+1)
!					k(0:3,i+1)=vect(0:3)
!					IF (i>N+1) i=i-2
!				END IF
!			END IF
!			i=i+1
!		END DO
!		!ordino in modo crescente le k
!		i=1
!		DO WHILE (i<N_bigger)
!			IF (k(0,i)>k(0,i+1)) THEN
!				vect(0:3)=k(0:3,i)
!				k(0:3,i)=k(0:3,i+1)
!				k(0:3,i+1)=vect(0:3)
!				IF (i>1) i=i-2
!			END IF
!			i=i+1
!		END DO
!		k_fermi_twist(0:3,1:N)=k(0:3,1:N)
!		DEALLOCATE(k,k_fermi_bigger)
!	END SUBROUTINE fermi_quantization_twist
!!-----------------------------------------------------------------------
	SUBROUTINE chiudi_dnfH()
		IMPLICIT NONE
		
		DEALLOCATE(M_Hartree,k_Hartree,i_fpw_Hartree,autovalori_Hartree,n_fpw_Hartee)
      CALL kwv_Hartree%deallocateKWaVe()
		flag_inizializza=.FALSE.
		
	END SUBROUTINE chiudi_dnfH
END MODULE dnfH
