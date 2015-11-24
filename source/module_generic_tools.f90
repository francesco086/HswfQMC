MODULE generic_tools
	IMPLICIT NONE
	REAL (KIND=8), PRIVATE, PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0
	
	CONTAINS
	
	SUBROUTINE genera_matrice_rotazione_eulero(N,teta)
		IMPLICIT NONE
		INTEGER (KIND=4), INTENT(IN) :: N
		INTEGER :: i, j
		REAL (KIND=8), INTENT(OUT) :: teta(1:N,1:N)
		
		DO j = 1, N-1, 1
			DO i = j+1, N, 1
				CALL RANDOM_NUMBER(teta(i,j))
				teta(i,j)=teta(i,j)*2.d0*PI
			END DO
		END DO
		
	END SUBROUTINE genera_matrice_rotazione_eulero
!-----------------------------------------------------------------------
	SUBROUTINE ruota_vettore_Nd_eulero(vec,N,teta)
		IMPLICIT NONE
		INTEGER (KIND=4), INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: teta(1:N,1:N)
		INTEGER :: i, j
		REAL (KIND=8) :: vec(1:N), vec_app(1:2)

		DO j = 1, N-1, 1
			DO i = j+1, N, 1
				vec_app(1:2)=vec((/i,j/))
				CALL ruota_vettore_2d(vec_app(1:2),teta(i,j))
				vec((/i,j/))=vec_app(1:2)
			END DO
		END DO

	END SUBROUTINE ruota_vettore_Nd_eulero
!-----------------------------------------------------------------------
	SUBROUTINE ruota_vettore_2d(vec,teta)
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: teta
		REAL (KIND=8) :: vec(1:2), vec_old(1:2)

		vec_old=vec
		vec(1)=vec_old(1)*DCOS(teta)-vec_old(2)*DSIN(teta)
		vec(2)=vec_old(1)*DSIN(teta)+vec_old(2)*DCOS(teta)

	END SUBROUTINE ruota_vettore_2d
!-----------------------------------------------------------------------
	RECURSIVE FUNCTION fattoriale_doppio(n) RESULT(fatt)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: n
		INTEGER :: fatt
		IF ((n<=2) .AND. (n>0)) THEN
			fatt=n
		ELSE IF (n>2) THEN
			fatt=n*fattoriale_doppio(n-2)
		ELSE
			STOP 'Non puoi fare il fattoriale di un numero minore o uguale a [ module_generic_tools.f90 > fattoriale_doppio ]'
		END IF
	END FUNCTION fattoriale_doppio
!-----------------------------------------------------------------------
	SUBROUTINE stampa_logo_HswfQMC()
		IMPLICIT NONE
		
		PRINT * , '                   ____        ______                                                          '
		PRINT * , 'HHHHH     HHHHH   |  __|      |   ___|   QQQQQQQQQQQQQ                                         '
		PRINT * , 'HHHHH     HHHHH    \_  \      |   ___|  QQQQQQQQQQQQQQQ                           CCCCCCCCCCCCC'
		PRINT * , 'HHHHH     HHHHH    __| |      |  |      QQQQQ      QQQQQ      MMMMM     MMMMM    CCCCCCCCCCCCCC'
		PRINT * , 'HHHHH     HHHHH   |___/       |__|      QQQQQ      QQQQQ      MMMMMM   MMMMMM   CCCCCC'
		PRINT * , 'HHHHHHHHHHHHHHH     ____    ____        QQQQQ      QQQQQ      MMMMMMM MMMMMMM   CCCCCC'
		PRINT * , 'HHHHHHHHHHHHHHH    \     \/     /       QQQQQQQQQQQQQQQQQQ    MMMMMMMMMMMMMMM   CCCCCC'
		PRINT * , 'HHHHH     HHHHH     \          /         QQQQQQQQQQQQQ QQQQ   MMMMM      MMMM    CCCCCCCCCCCCCC'
		PRINT * , 'HHHHH     HHHHH      \   /\   /                               MMMMM      MMMM     CCCCCCCCCCCCC'
		PRINT * , 'HHHHH     HHHHH       \_/  \_/                                MMMMM      MMMM   '
		PRINT * , ''
		
	END SUBROUTINE stampa_logo_HswfQMC
!-----------------------------------------------------------------------
	!dato il numero di particelle N (che deve essere pari) salva tutte M=(N-1)!! le possibilitá di creare coppie nel vettore coppie(1:2,1:N/2,1:M)
	SUBROUTINE crea_tutte_le_coppie(N,M,coppie)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N
		INTEGER :: M, num_sub, i              !num_sub é il numero di possibilitá che ci sono fissata la prima coppia
		INTEGER :: coppie(1:2,1:N/2,1:M)
		IF (MOD(N,2)/=0) STOP 'N non é pari [ module_generic_tools.f90 > crea_tutte_le_coppie ]'
		num_sub=fattoriale_doppio(N-3)
		DO i = 2, N, 1
			coppie(1,1,(i-2)*num_sub+1:(i-1)*num_sub)=1
			coppie(2,1,(i-2)*num_sub+1:(i-1)*num_sub)=i
			!PRINT * , 'A: ', 1, i, ' da riga ', (i-2)*num_sub+1, ' a ', (i-1)*num_sub
			CALL crea_coppie_con_esclusioni(N,num_sub,coppie(1:2,2:N/2,(i-2)*num_sub+1:(i-1)*num_sub),2,(/1,i/))
		END DO
	END SUBROUTINE crea_tutte_le_coppie
!-----------------------------------------------------------------------
	RECURSIVE SUBROUTINE crea_coppie_con_esclusioni(N,num_sub,coppie,num_escl,v_escl)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num_sub, num_escl
		LOGICAL :: flag
		INTEGER :: num_sub_new
		INTEGER :: coppie(1:2,1:(N-num_escl)/2,1:num_sub), v_escl(1:num_escl), v_escl_new(1:num_escl+2)
		INTEGER :: vett_new(1:N-num_escl), i, j, cont
		cont=1
		DO i = 1, N, 1
			flag=.TRUE.
			DO j = 1, num_escl, 1
				IF (i==v_escl(j)) THEN
					flag=.FALSE.
				END IF
			END DO
			IF (flag) THEN
				vett_new(cont)=i
				cont=cont+1
			END IF
		END DO		
		IF (N-num_escl-3>0) THEN
			num_sub_new=fattoriale_doppio(N-num_escl-3)
			v_escl_new(1:num_escl)=v_escl(1:num_escl)
			DO i = 2, N-num_escl, 1
				coppie(1,1,(i-2)*num_sub_new+1:(i-1)*num_sub_new)=vett_new(1)
				coppie(2,1,(i-2)*num_sub_new+1:(i-1)*num_sub_new)=vett_new(i)
				v_escl_new(num_escl+1)=vett_new(1)
				v_escl_new(num_escl+2)=vett_new(i)
				!PRINT * , 'B:', vett_new(1), vett_new(i), ' da riga ', (i-2)*num_sub_new+1, ' a ', (i-1)*num_sub_new
				CALL crea_coppie_con_esclusioni(N,num_sub_new,coppie(1:2,2:(N-num_escl)/2,(i-2)*num_sub_new+1:(i-1)*num_sub_new),&
				   num_escl+2,v_escl_new(1:num_escl+2))
			END DO
		ELSE
			coppie(1:2,1,1)=vett_new(1:2)
			!PRINT * , 'C:', coppie(1,1,1), coppie(2,1,1), ' da riga 1 a 2'
		END IF

	END SUBROUTINE crea_coppie_con_esclusioni
!-----------------------------------------------------------------------
	!trova il path alla posizione attuale
	SUBROUTINE trova_pwd(path)
		IMPLICIT NONE
		CHARACTER(LEN=256), INTENT(OUT) :: path
		!Nella funzione richiamante si puó poi usare TRIM(path) in modo da eliminare gli spazi bianchi
		CALL SYSTEM ('pwd > .mypwd')
		OPEN (11, file='.mypwd', status='OLD')
		READ (11, '(A)') path
		CLOSE (11)
	END SUBROUTINE trova_pwd
!!-----------------------------------------------------------------------
!	!trova il path alla posizione attuale
!	SUBROUTINE determina_numero_cartelle_K(path,myrank,num_K_points)
!		IMPLICIT NONE
!		INTEGER :: ios
!		LOGICAL :: flag_file_creato, flag_lock, flag_read
!		CHARACTER(LEN=*) :: path
!		CHARACTER(LEN=4) :: string_myrank
!		INTEGER, INTENT(IN) :: myrank
!		INTEGER, INTENT(OUT) :: num_K_points
!		!Nella funzione richiamante si puó poi usare TRIM(path) in modo da eliminare gli spazi bianchi
!		IF (myrank>9999) STOP 'Stai usando troppi processori - determina_numero_cartelle_K [module_generic_tools.f90]'
!		WRITE (string_myrank, '(I4.4)'), myrank
!		INQUIRE(FILE=path//'/.lock_qmc',EXIST=flag_lock)
!		DO WHILE (flag_lock )
!			CALL SLEEP(1)
!			INQUIRE(FILE=path//'/.lock_qmc',EXIST=flag_lock)
!		END DO
!		CALL SYSTEM('touch '//path//'/.lock_qmc')
!		CALL SYSTEM ('rm -f '//path//'/.numero_cartelle_K'//string_myrank)
!		CALL SYSTEM ('ls -d  '//path//'/K* | wc -w > '//path//'/.numero_cartelle_K'//string_myrank)
!		INQUIRE(FILE=path//'/.numero_cartelle_K'//string_myrank,EXIST=flag_file_creato)
!		DO WHILE ( .NOT. flag_file_creato )
!			CALL SLEEP(1)
!			INQUIRE(FILE=path//'/.numero_cartelle_K'//string_myrank,EXIST=flag_file_creato)
!		END DO
!		OPEN(UNIT=66, FILE=path//'/.numero_cartelle_K'//string_myrank, STATUS='OLD', IOSTAT=ios)
!		IF ( ios /= 0 ) STOP "Error opening file in determina_numero_cartelle_K [module_generic_tools.f90]"
!		flag_read=.FALSE.
!		DO WHILE ( .NOT. flag_read )
!			READ (UNIT=66, FMT=*, IOSTAT=ios) num_K_points
!			IF ( ios==0 ) THEN
!				flag_read=.TRUE.
!				CLOSE (66)
!			ELSE
!				flag_read=.FALSE.
!				CLOSE (66)
!				CALL SYSTEM ('rm -f '//path//'/.numero_cartelle_K'//string_myrank)
!				CALL SYSTEM ('ls -d  '//path//'/K* | wc -w > '//path//'/.numero_cartelle_K'//string_myrank)
!				INQUIRE(FILE=path//'/.numero_cartelle_K'//string_myrank,EXIST=flag_file_creato)
!				DO WHILE ( .NOT. flag_file_creato )
!					CALL SLEEP(1)
!					INQUIRE(FILE=path//'/.numero_cartelle_K'//string_myrank,EXIST=flag_file_creato)
!				END DO
!				OPEN(UNIT=66, FILE=path//'/.numero_cartelle_K'//string_myrank, STATUS='OLD', IOSTAT=ios)
!				IF ( ios /= 0 ) STOP "Error opening file in determina_numero_cartelle_K [module_generic_tools.f90]"
!			END IF
!		END DO
!		CALL SYSTEM ('rm -f '//path//'/.numero_cartelle_K'//string_myrank)
!		CALL SYSTEM('rm -f '//path//'/.lock_qmc')
!	END SUBROUTINE determina_numero_cartelle_K
!-----------------------------------------------------------------------
	SUBROUTINE leggi_pesi_K(path,num_K,pesi)
		IMPLICIT NONE
		CHARACTER(LEN=*) :: path
		INTEGER, INTENT(IN) :: num_K
		LOGICAL :: flag_find
		CHARACTER(LEN=27) :: string_find
		CHARACTER(LEN=150) :: string_weight
		INTEGER :: ios, cont
		INTEGER (KIND=8) :: i1, i2, i3
		REAL :: weight
		REAL (KIND=8) :: pesi(1:num_K)
		
		OPEN(UNIT=66, FILE=path//'/data-file.xml', STATUS='OLD', IOSTAT=ios)
		IF ( ios /= 0 ) STOP "Error opening file data-file.xml leggi_pesi_K [module_generic_tools.f90]"
		flag_find=.FALSE.
		cont=0
		DO WHILE ( .NOT. flag_find )
			cont=cont+1
			READ (66, *), string_find
			IF ( string_find=='<MONKHORST_PACK_OFFSET' ) THEN
				flag_find=.TRUE.
			END IF
		END DO
		DO i1 = 1, num_K, 1
			READ (66, '(A150)'), string_weight(1:150)
			flag_find=.FALSE.
			i2=0
			DO WHILE ( (.NOT. flag_find) .AND. (i2<130) )
				i2=i2+1
				IF ( string_weight(i2:i2+6)=='WEIGHT=' ) THEN
					i2=i2+8
					i3=1
					DO WHILE ( string_weight(i2+i3+2:i2+i3+3)/='/>' )
						i3=i3+1
					END DO
					READ (string_weight(i2:i2+i3), *), weight
					pesi(i1)=weight
				END IF
			END DO
		END DO
		CLOSE (66)

	END SUBROUTINE leggi_pesi_K
!-----------------------------------------------------------------------
	SUBROUTINE leggi_eigenval_xml(N,file,fattori)
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: file
		INTEGER, INTENT(IN) :: N
		LOGICAL :: flag_file
		INTEGER :: i
		REAL (KIND=8), INTENT(OUT) :: fattori(1:N)

		INQUIRE(FILE=file,EXIST=flag_file)
		IF (flag_file) THEN
			OPEN (UNIT=89, FILE=file, STATUS='OLD')
			DO i = 1, 9, 1
				READ (89,*), 
			END DO
			DO i = 1, N, 1
				READ (89,*), fattori(i)
			END DO
			CLOSE (89)
		ELSE
			STOP 'Manca il file [leggi_eigenval_xml]'
		END IF

	END SUBROUTINE leggi_eigenval_xml
!-----------------------------------------------------------------------
	SUBROUTINE leggi_N_pw(file,N_pw)
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: file
		LOGICAL :: flag_file
		CHARACTER(LEN=22) :: stringa
		INTEGER :: i, j, num_ch, cont
		INTEGER, INTENT(OUT) :: N_pw

		INQUIRE(FILE=file,EXIST=flag_file)
		IF ( flag_file ) THEN
			!OPEN (UNIT=89, FILE=file, STATUS='OLD')   !leggeva da evc.xml
			!DO i = 1, 6, 1
			!	READ (89,*), 
			!END DO
			!READ (89,89), stringa;		89 FORMAT(2X,20A)
			!cont=0
			!num_ch=50
			!DO WHILE ((num_ch<58).AND.(num_ch>47))
			!	num_ch=IACHAR(stringa(12+cont:12+cont))
			!	cont=cont+1
			!END DO
			!cont=cont-1
			!N_pw=0
			!DO i = 0, cont-1, 1
			!	num_ch=IACHAR(stringa(12+i:12+i))
			!	N_pw=N_pw+(num_ch-48)*(10**(cont-i-1))             !N_pw é il numero di pw usate per ogni orbitale
			!END DO
			!CLOSE (89)
			OPEN (UNIT=89, FILE=file, STATUS='OLD')         !legge da gkvectors.xml
			DO i = 1, 7, 1
				READ (89,*), 
			END DO
			READ (89, *), N_pw
			!PRINT * , 'N_pw=', N_pw
			CLOSE (89)
		ELSE
			PRINT * , file
			STOP 'Manca il file [leggi_N_pw]'
		END IF
	END SUBROUTINE leggi_N_pw
!-----------------------------------------------------------------------
	SUBROUTINE leggi_evc_xml(N,N_pw,file,fattori)
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: file
		INTEGER, INTENT(IN) :: N, N_pw
		LOGICAL :: flag_file
		INTEGER :: i, j
		REAL (KIND=8) :: dato(1:2)
		COMPLEX (KIND=8) :: fattori(1:N_pw,1:N)
				
		INQUIRE(FILE=file,EXIST=flag_file)
		IF ( flag_file ) THEN
			OPEN (UNIT=89, FILE=file, STATUS='OLD')
			DO i = 1, 6, 1
				READ (89,*), 
			END DO
			DO i = 1, N, 1
				READ(89,*),
				READ(89,*),
				!PRINT * , 'i=', i
				DO j = 1, N_pw, 1
					!PRINT * , 'j=', j
					READ (89,*), dato(1:2)
					fattori(j,i)=(1.d0,0.d0)*dato(1)+(0.d0,1.d0)*dato(2)
				END DO
			END DO
			CLOSE (89)
		ELSE
			STOP 'Manca il file [leggi_evc_xml]'
		END IF

	END SUBROUTINE leggi_evc_xml
!-----------------------------------------------------------------------
	SUBROUTINE leggi_gkvectors_xml(N_pw,file,k,trasl_k)
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: file
		INTEGER, INTENT(IN) :: N_pw
		LOGICAL :: flag_file
		INTEGER :: i, j
		REAL (KIND=8) :: dato(1:2)
		INTEGER :: k(1:3,1:N_pw)
		REAL (KIND=8) :: trasl_k(1:3)
		
		INQUIRE(FILE=file,EXIST=flag_file)
		IF ( flag_file ) THEN
			OPEN (UNIT=89, FILE=file, STATUS='OLD')
			DO i = 1, 16, 1
				READ (89,*), 
			END DO
			READ (89, *), trasl_k(1)
			READ (89, *), trasl_k(2)
			READ (89, *), trasl_k(3)
			DO i = 1, N_pw+4, 1
				READ (89,*), 
			END DO
			DO i = 1, N_pw, 1
				READ (89,*), k(1:3,i)
			END DO
			CLOSE (89)
		ELSE
			STOP 'Manca il file [leggi_gkvectors_xml]'
		END IF

	END SUBROUTINE leggi_gkvectors_xml
!-----------------------------------------------------------------------
	!dato il vettore con il numero di particelle e le tre dimensioni della scatola di simulazione, usa le PBC per eventualmente rimettere il vettore dentro la scatola
	SUBROUTINE applica_pbc(v,N,L)
		IMPLICIT NONE
		INTEGER :: N, i, j
		REAL (KIND=8) :: L(1:3), v(1:3,1:N)
		DO j = 1, N, 1
			DO i = 1, 3, 1
				v(i,j)=v(i,j)-L(i)*DNINT(v(i,j)/L(i))
			END DO
		END DO
	END SUBROUTINE applica_pbc
!-----------------------------------------------------------------------
	!calcola la matrice delle distanze fra un set di particelle
	SUBROUTINE valuta_distanza_ii(r,N,L,rij)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: r(1:3,1:N), L(1:3)
		LOGICAL :: flag_rpa
		INTEGER :: i, j
		REAL (KIND=8) :: frf1
		REAL (KIND=8), INTENT(OUT) :: rij(0:3,1:N,1:N)
		DO j = 1, N-1, 1
			DO i = j+1, N, 1
				rij(1:3,i,j)=r(1:3,j)-r(1:3,i)
				rij(1:3,i,j)=rij(1:3,i,j)-L(1:3)*DNINT(rij(1:3,i,j)/L(1:3))
				rij(0,i,j)=DSQRT(DOT_PRODUCT(rij(1:3,i,j),rij(1:3,i,j)))
				rij(0,j,i)=rij(0,i,j)
				rij(1:3,j,i)=-rij(1:3,i,j)
			END DO
			rij(0:3,j,j)=0.d0
		END DO
	END SUBROUTINE valuta_distanza_ii
!-----------------------------------------------------------------------
	!calcola la matrice delle distanze fra un set di particelle, quando ne é stata spostata solo una
	SUBROUTINE valuta_distanza_ii_1ppt(num,r,N,L,rij)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: r(1:3,1:N), L(1:3)
		LOGICAL :: flag_rpa
		INTEGER :: i
		REAL (KIND=8), INTENT(OUT) :: rij(0:3,1:N,1:N)
		DO i = 1, num-1, 1
			rij(1:3,i,num)=r(1:3,num)-r(1:3,i)
			rij(1:3,i,num)=rij(1:3,i,num)-L(1:3)*DNINT(rij(1:3,i,num)/L(1:3))
			rij(0,i,num)=DSQRT(DOT_PRODUCT(rij(1:3,i,num),rij(1:3,i,num)))
			rij(0,num,i)=rij(0,i,num)
			rij(1:3,num,i)=-rij(1:3,i,num)
		END DO
		rij(0:3,num,num)=0.d0
		DO i = num+1, N, 1
			rij(1:3,i,num)=r(1:3,num)-r(1:3,i)
			rij(1:3,i,num)=rij(1:3,i,num)-L(1:3)*DNINT(rij(1:3,i,num)/L(1:3))
			rij(0,i,num)=DSQRT(DOT_PRODUCT(rij(1:3,i,num),rij(1:3,i,num)))
			rij(0,num,i)=rij(0,i,num)
			rij(1:3,num,i)=-rij(1:3,i,num)
		END DO
	END SUBROUTINE valuta_distanza_ii_1ppt
!-----------------------------------------------------------------------
	!calcola il corrispondente periodico della matrice delle distanze fra un set di particelle
	SUBROUTINE valuta_distanza_pc_ii(rij,N,L,rij_pc)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N), L(1:3)
		INTEGER :: i, j
		REAL (KIND=8) :: frf1(1:3), frf2(1:3)
		REAL (KIND=8), INTENT(OUT) :: rij_pc(0:3,1:N,1:N)
		frf1(1:3)=PI/L(1:3)
		frf2(1:3)=L(1:3)/PI
		DO j = 1, N, 1
			DO i = j+1, N, 1
				rij_pc(1:3,i,j)=frf2(1:3)*DSIN(frf1(1:3)*rij(1:3,i,j))
				rij_pc(0,i,j)=DSQRT(DOT_PRODUCT(rij_pc(1:3,i,j),rij_pc(1:3,i,j)))
				rij_pc(0,j,i)=rij_pc(0,i,j)
				rij_pc(1:3,j,i)=-rij_pc(1:3,i,j)
			END DO
		END DO
	END SUBROUTINE valuta_distanza_pc_ii
!-----------------------------------------------------------------------
	!calcola il corrispondente periodico della matrice delle distanze fra un set di particelle
	SUBROUTINE valuta_distanza_pc_ii_1ppt(num,rij,N,L,rij_pc)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N), L(1:3)
		INTEGER :: i
		REAL (KIND=8) :: frf1(1:3), frf2(1:3)
		REAL (KIND=8), INTENT(OUT) :: rij_pc(0:3,1:N,1:N)
		frf1=PI/L(1:3)
		frf2=L(1:3)/PI
		DO i = 1, num-1, 1
			rij_pc(1:3,i,num)=frf2(1:3)*DSIN(frf1(1:3)*rij(1:3,i,num))
			rij_pc(0,i,num)=DSQRT(DOT_PRODUCT(rij_pc(1:3,i,num),rij_pc(1:3,i,num)))
			rij_pc(0,num,i)=rij_pc(0,i,num)
			rij_pc(1:3,num,i)=-rij_pc(1:3,i,num)
		END DO
		rij_pc(0:3,num,num)=0.d0
		DO i = num+1, N, 1
			rij_pc(1:3,i,num)=frf2(1:3)*DSIN(frf1(1:3)*rij(1:3,i,num))
			rij_pc(0,i,num)=DSQRT(DOT_PRODUCT(rij_pc(1:3,i,num),rij_pc(1:3,i,num)))
			rij_pc(0,num,i)=rij_pc(0,i,num)
			rij_pc(1:3,num,i)=-rij_pc(1:3,i,num)
		END DO
	END SUBROUTINE valuta_distanza_pc_ii_1ppt
!-----------------------------------------------------------------------
	!calcola la matrice delle distanze fra due set di particelle
	SUBROUTINE valuta_distanza_ij(re,rp,N,L,rij)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: re(1:3,1:N), rp(1:3,1:N), L(1:3)
		LOGICAL :: flag_rpa
		INTEGER :: i, j
		REAL (KIND=8), INTENT(OUT) :: rij(0:3,1:N,1:N)
		DO j = 1, N, 1
			DO i = 1, N, 1
				rij(1:3,i,j)=re(1:3,i)-rp(1:3,j)
				rij(1:3,i,j)=rij(1:3,i,j)-L(1:3)*DNINT(rij(1:3,i,j)/L(1:3))
				rij(0,i,j)=DSQRT(DOT_PRODUCT(rij(1:3,i,j),rij(1:3,i,j)))
			END DO
		END DO
	END SUBROUTINE valuta_distanza_ij
!-----------------------------------------------------------------------
	!calcola la matrice delle distanze fra due set di particelle
	SUBROUTINE valuta_distanza_ij_1ppt(num,eop,re,rp,N,L,rij)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, eop
		REAL (KIND=8), INTENT(IN) :: re(1:3,1:N), rp(1:3,1:N), L(1:3)
		LOGICAL :: flag_rpa
		INTEGER :: i
		REAL (KIND=8), INTENT(OUT) :: rij(0:3,1:N,1:N)
		IF (eop==1) THEN
			DO i = 1, N, 1
				rij(1:3,num,i)=re(1:3,num)-rp(1:3,i)
				rij(1:3,num,i)=rij(1:3,num,i)-L(1:3)*DNINT(rij(1:3,num,i)/L(1:3))
				rij(0,num,i)=DSQRT(DOT_PRODUCT(rij(1:3,num,i),rij(1:3,num,i)))
			END DO
		ELSE IF (eop==2) THEN
			DO i = 1, N, 1
				rij(1:3,i,num)=re(1:3,i)-rp(1:3,num)
				rij(1:3,i,num)=rij(1:3,i,num)-L(1:3)*DNINT(rij(1:3,i,num)/L(1:3))
				rij(0,i,num)=DSQRT(DOT_PRODUCT(rij(1:3,i,num),rij(1:3,i,num)))
			END DO
		END IF
	END SUBROUTINE valuta_distanza_ij_1ppt
!-----------------------------------------------------------------------
	!calcola la matrice delle distanze fra due set di particelle
	SUBROUTINE valuta_distanza_diagonale_ij(re,rp,N,L,rij)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: re(1:3,1:N), rp(1:3,1:N), L(1:3)
		LOGICAL :: flag_rpa
		INTEGER :: i, j
		REAL (KIND=8), INTENT(OUT) :: rij(0:3,1:N)
		DO i = 1, N, 1
			rij(1:3,i)=re(1:3,i)-rp(1:3,i)
			rij(1:3,i)=rij(1:3,i)-L(1:3)*DNINT(rij(1:3,i)/L(1:3))
			rij(0,i)=DSQRT(DOT_PRODUCT(rij(1:3,i),rij(1:3,i)))
		END DO
	END SUBROUTINE valuta_distanza_diagonale_ij
!-----------------------------------------------------------------------
	!calcola il corrispondente periodico della matrice delle distanze fra due set di particelle
	SUBROUTINE valuta_distanza_pc_ij(rij,N,L,rij_pc)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N), L(1:3)
		INTEGER :: i, j
		REAL (KIND=8) :: frf1(1:3), frf2(1:3)
		REAL (KIND=8), INTENT(OUT) :: rij_pc(0:3,1:N,1:N)
		frf1(1:3)=PI/L(1:3)
		frf2(1:3)=L(1:3)/PI
		DO j = 1, N, 1
			DO i = 1, N, 1
				rij_pc(1:3,i,j)=frf2(1:3)*DSIN(frf1(1:3)*rij(1:3,i,j))
				rij_pc(0,i,j)=DSQRT(DOT_PRODUCT(rij_pc(1:3,i,j),rij_pc(1:3,i,j)))
			END DO
		END DO
	END SUBROUTINE valuta_distanza_pc_ij
!-----------------------------------------------------------------------
	!calcola il corrispondente periodico della matrice delle distanze fra due set di particelle
	SUBROUTINE valuta_distanza_pc_ij_1ppt(num,eop,rij,N,L,rij_pc)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, num, eop
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N), L(1:3)
		INTEGER :: i
		REAL (KIND=8) :: frf1(1:3), frf2(1:3)
		REAL (KIND=8), INTENT(OUT) :: rij_pc(0:3,1:N,1:N)
		frf1(1:3)=PI/L(1:3)
		frf2(1:3)=L(1:3)/PI
		IF (eop==1) THEN
			DO i = 1, N, 1
				rij_pc(1:3,num,i)=frf2(1:3)*DSIN(frf1(1:3)*rij(1:3,num,i))
				rij_pc(0,num,i)=DSQRT(DOT_PRODUCT(rij_pc(1:3,num,i),rij_pc(1:3,num,i)))
			END DO
		ELSE IF (eop==2) THEN
			DO i = 1, N, 1
				rij_pc(1:3,i,num)=frf2(1:3)*DSIN(frf1(1:3)*rij(1:3,i,num))
				rij_pc(0,i,num)=DSQRT(DOT_PRODUCT(rij_pc(1:3,i,num),rij_pc(1:3,i,num)))
			END DO
		END IF
	END SUBROUTINE valuta_distanza_pc_ij_1ppt
!-----------------------------------------------------------------------
	!trova le coppie che minimizzano la distanza procedendo in ordine stocastico

	SUBROUTINE coppie_minima_distanza_stoc(N,r_ij,i_coppie,r_coppie)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N        !numero di particelle, DEVE ESSERE UN NUMERO PARI
		REAL (KIND=8), INTENT(IN) :: r_ij(1:N,1:N)      !distanze fra le particelle
		INTEGER :: i1, i_eta, num_free_p, index_free(1:N), index_min, cont
		REAL :: eta, R_N
		REAL (KIND=8) :: dist_min
		INTEGER, INTENT(OUT) :: i_coppie(1:N/2,2)   !riporta gli indici delle coppie
		REAL (KIND=8), INTENT(OUT) :: r_coppie(1:N/2)     !distanza delle coppie piú vicine

		IF ( MOD(N,2)==1 ) THEN
			STOP "Errore, N deve essere un numero pari &
				[ module_generic_tools.f90 > crea_coppie_minima_distanza ]"
		END IF
		R_N=REAL(N)
		num_free_p=N
		DO i1 = 1, N, 1
			index_free(i1)=i1
		END DO

		cont=1
		DO WHILE ( num_free_p>0 )
			i_eta=0
			DO WHILE ( (i_eta<1).OR.(i_eta>N) )
				CALL RANDOM_NUMBER(eta)
				i_eta=CEILING(eta*R_N)
				IF ( index_free(i_eta)==0 ) i_eta=0
			END DO
			index_free(i_eta)=0

			dist_min=-1.d0
			DO i1 = 1, N, 1
				IF ( index_free(i1)>0 ) THEN
					IF ( (r_ij(i1,i_eta)<dist_min).OR.(dist_min<0.d0) ) THEN
						dist_min=r_ij(i1,i_eta)
						index_min=i1
					END IF
				END IF
			END DO
			index_free(index_min)=0

			i_coppie(cont,1)=i_eta
			i_coppie(cont,2)=index_min
			r_coppie(cont)=dist_min

			num_free_p=num_free_p-2
			cont=cont+1
		END DO

	END SUBROUTINE coppie_minima_distanza_stoc
!-----------------------------------------------------------------------
	!trova le coppie che minimizzano la distanza procedendo in ordine stocastico, e considerando solo coppie spin-up/spin-down

	SUBROUTINE coppie_updw_minima_distanza_stoc(N,r_ij,i_coppie,r_coppie)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N        !numero di particelle, DEVE ESSERE UN NUMERO PARI
		REAL (KIND=8), INTENT(IN) :: r_ij(1:N,1:N)      !distanze fra le particelle
		INTEGER :: i1, i_eta, num_free_p, index_free(1:N), index_min, cont, H_N
		REAL :: eta, R_N
		REAL (KIND=8) :: dist_min
		INTEGER, INTENT(OUT) :: i_coppie(1:N/2,2)   !riporta gli indici delle coppie
		REAL (KIND=8), INTENT(OUT) :: r_coppie(1:N/2)     !distanza delle coppie piú vicine

		IF ( MOD(N,2)==1 ) THEN
			STOP "Errore, N deve essere un numero pari &
				[ module_generic_tools.f90 > crea_coppie_minima_distanza ]"
		END IF
		R_N=REAL(N)
		H_N=N/2
		num_free_p=N
		DO i1 = 1, N, 1
			index_free(i1)=i1
		END DO

		cont=1
		DO WHILE ( num_free_p>0 )
			i_eta=0
			DO WHILE ( (i_eta<1).OR.(i_eta>N) )
				CALL RANDOM_NUMBER(eta)
				i_eta=CEILING(eta*R_N)
				IF ( index_free(i_eta)==0 ) i_eta=0
			END DO
			index_free(i_eta)=0

			dist_min=-1.d0
			IF ( i_eta>H_N ) THEN
				DO i1 = H_N+1, N, 1
					IF ( index_free(i1)>0 ) THEN
						IF ( (r_ij(i1,i_eta)<dist_min).OR.(dist_min<0.d0) ) THEN
							dist_min=r_ij(i1,i_eta)
							index_min=i1
						END IF
					END IF
				END DO
			ELSE
				DO i1 = 1, H_N, 1
					IF ( index_free(i1)>0 ) THEN
						IF ( (r_ij(i1,i_eta)<dist_min).OR.(dist_min<0.d0) ) THEN
							dist_min=r_ij(i1,i_eta)
							index_min=i1
						END IF
					END IF
				END DO
			END IF
			index_free(index_min)=0

			i_coppie(cont,1)=i_eta
			i_coppie(cont,2)=index_min
			r_coppie(cont)=dist_min

			num_free_p=num_free_p-2
			cont=cont+1
		END DO

	END SUBROUTINE coppie_updw_minima_distanza_stoc
!-----------------------------------------------------------------------

	SUBROUTINE aggiorna_coppie_minima_distanza_stoc_1ppt(num,N,r_ij,i_coppie,r_coppie,i_num)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N        !numero di particelle (DEVE ESSERE UN NUMERO PARI)
		INTEGER, INTENT(IN) :: num		!indice della particella mossa
		INTEGER :: i_num        !indice della coppia corrispondente a num
		INTEGER, INTENT(IN) :: i_coppie(1:N/2,1:2)   !riporta gli indici delle coppie
		REAL (KIND=8), INTENT(IN) :: r_ij(1:N,1:N)      !distanze fra le particelle
		INTEGER :: i1
		REAL (KIND=8), INTENT(OUT) :: r_coppie(1:N/2)     !distanza delle coppie piú vicine

		IF ( MOD(N,2)==1 ) THEN
			STOP "Errore, N deve essere un numero pari &
				[ module_generic_tools.f90 > crea_coppie_minima_distanza ]"
		END IF
				
		DO i1 = 1, N/2, 1
			IF ( (i_coppie(i1,1)==num).OR.(i_coppie(i1,2)==num) ) THEN
				r_coppie(i1)=r_ij(i_coppie(i1,1),i_coppie(i1,2))
				i_num=i1
			END IF
		END DO

	END SUBROUTINE aggiorna_coppie_minima_distanza_stoc_1ppt
!-----------------------------------------------------------------------
	!calcola la g(r) per un set di particelle
	SUBROUTINE g_ii_correlation(rij,N,L,r_s,N_hist,guu,gud)
		INTEGER, INTENT(IN) :: N
		INTEGER (KIND=8), INTENT(IN) :: N_hist
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N), L(1:3), r_s
		INTEGER :: i, j, k, i_hist
		REAL (KIND=8) :: r_hist, dL, minL, rs3
		REAL (KIND=8), INTENT(OUT) :: guu(0:N_hist), gud(0:N_hist)
		
		INTEGER :: contatore
    	
		minL=MINVAL(L)
		rs3=r_s*r_s*r_s
		dL=0.5d0*minL/N_hist
		guu=0.d0
		gud=0.d0
		        
		DO  j = 1, N/2, 1
			DO  i = j+1, N/2, 1
				IF ( rij(0,i,j)<0.5d0*minL ) THEN
					i_hist=NINT(N_hist*rij(0,i,j)/(0.5d0*minL))
					r_hist=i_hist*dL
					guu(i_hist)=guu(i_hist)+rs3*2.d0/((r_hist)**3-(r_hist-dL)**3)
				END IF
			END DO
		END DO
        
		DO  j = N/2+1, N, 1
			DO  i = j+1, N, 1
				IF ( rij(0,i,j)<0.5d0*minL ) THEN
					i_hist=NINT(N_hist*rij(0,i,j)/(0.5D0*minL))
					r_hist=i_hist*dL
					guu(i_hist)=guu(i_hist)+rs3*2.d0/((r_hist)**3-(r_hist-dL)**3)
				END IF
			END DO
		END DO
		    
		DO  j = 1, N/2, 1
			DO  i = N/2+1, N, 1
				IF ( rij(0,i,j)<0.5d0*minL ) THEN
					i_hist=NINT(N_hist*rij(0,i,j)/(0.5d0*minL))
					r_hist=i_hist*dL
					gud(i_hist)=gud(i_hist)+rs3*2.d0/((r_hist)**3-(r_hist-dL)**3)
				END IF
			END DO
		END DO
        
		guu=guu/REAL(N,8)
		gud=gud/REAL(N,8)
		    
	END SUBROUTINE g_ii_correlation
!-----------------------------------------------------------------------
	!calcola la g(r) tra due set di particelle
	SUBROUTINE g_ij_correlation(rij,N,L,r_s,N_hist,gep)
		INTEGER, INTENT(IN) :: N
		INTEGER (KIND=8), INTENT(IN) :: N_hist
		REAL (KIND=8), INTENT(IN) :: rij(0:3,1:N,1:N), L(1:3), r_s
		INTEGER :: i, j, k, i_hist
		REAL (KIND=8) :: r_hist, dL, minL, rs3
		REAL (KIND=8), INTENT(OUT) :: gep(0:N_hist)
    	
		minL=MINVAL(L)
		rs3=r_s*r_s*r_s
		dL=0.5d0*minL/N_hist
		gep=0.d0
    	
		DO  j = 1, N, 1
			DO  i = 1, N, 1
				IF ( rij(0,i,j)<0.5d0*minL ) THEN
					i_hist=NINT(N_hist*rij(0,i,j)/(0.5d0*minL))
					r_hist=i_hist*dL
					gep(i_hist)=gep(i_hist)+rs3*2.d0/((r_hist)**3-(r_hist-dL)**3)
				END IF
			END DO
		END DO
		    
		gep=gep/REAL(N,8)
    
	END SUBROUTINE g_ij_correlation
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_determinante_C_1ppt(N,riga,IM_old,det_old,M_new,det_new)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, riga       !riga=riga che é stata cambiata
		COMPLEX (KIND=8), INTENT(IN) :: IM_old(1:N,1:N), M_new(1:N,1:N), det_old       !IM_old=matrice inversa old, M_new=nuova matrice
		INTEGER :: i
		COMPLEX (KIND=8), INTENT(OUT) :: det_new
		det_new=0.d0
		DO i = 1, N, 1
			det_new=det_new+IM_old(i,riga)*M_new(riga,i)
		END DO
		det_new=det_new*det_old
	END SUBROUTINE aggiorna_determinante_C_1ppt
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_determinante_C_col_1ppt(N,col,IM_old,det_old,M_new,det_new)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, col       !col=colonna che é stata cambiata
		COMPLEX (KIND=8), INTENT(IN) :: IM_old(1:N,1:N), M_new(1:N,1:N), det_old       !IM_old=matrice inversa old, M_new=nuova matrice
		INTEGER :: i
		COMPLEX (KIND=8), INTENT(OUT) :: det_new
		IF ( col>N ) STOP 'col>N      aggiorna_determinante_C_col_1ppt[module_generic_tools.f90]'
		det_new=0.d0
		DO i = 1, N, 1
			det_new=det_new+IM_old(col,i)*M_new(i,col)
		END DO
		det_new=det_new*det_old
	END SUBROUTINE aggiorna_determinante_C_col_1ppt
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_determinante_R_1ppt(N,riga,IM_old,det_old,M_new,det_new)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, riga       !riga=riga che é stata cambiata
		REAL (KIND=8), INTENT(IN) :: IM_old(1:N,1:N), M_new(1:N,1:N), det_old       !IM_old=matrice inversa old, M_new=nuova matrice
		INTEGER :: i
		REAL (KIND=8), INTENT(OUT) :: det_new
		det_new=0.d0
		DO i = 1, N, 1
			det_new=det_new+IM_old(i,riga)*M_new(riga,i)
		END DO
		det_new=det_new*det_old
	END SUBROUTINE aggiorna_determinante_R_1ppt
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_determinante_R_col_1ppt(N,col,IM_old,det_old,M_new,det_new)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, col       !col=colonna che é stata cambiata
		REAL (KIND=8), INTENT(IN) :: IM_old(1:N,1:N), M_new(1:N,1:N), det_old       !IM_old=matrice inversa old, M_new=nuova matrice
		INTEGER :: i
		REAL (KIND=8), INTENT(OUT) :: det_new
		det_new=0.d0
		DO i = 1, N, 1
			det_new=det_new+IM_old(col,i)*M_new(i,col)
		END DO
		det_new=det_new*det_old
	END SUBROUTINE aggiorna_determinante_R_col_1ppt
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_matrice_inversa_C_1ppt(N,riga,IM_old,det_old,M_new,det_new,IM_new)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, riga       !riga=riga che é stata cambiata
		COMPLEX (KIND=8), INTENT(IN) :: IM_old(1:N,1:N), M_new(1:N,1:N), det_old, det_new       !IM_old=matrice inversa old, M_new=nuova matrice
		COMPLEX (KIND=8) :: quoz, summ
		INTEGER :: i, j, l
		COMPLEX (KIND=8), INTENT(OUT) :: IM_new(1:N,1:N)
		
		quoz=det_old/det_new
		
		DO j = 1, riga-1, 1
			summ=0.d0
			DO i = 1, N, 1
				summ=summ+IM_old(i,j)*M_new(riga,i)
			END DO
			DO i = 1, N, 1
				IM_new(i,j)=IM_old(i,j)-IM_old(i,riga)*quoz*summ
			END DO
		END DO
		
		DO i = 1, N, 1
			IM_new(i,riga)=IM_old(i,riga)*quoz
		END DO
		
		DO j = riga+1, N, 1
			summ=0.d0
			DO i = 1, N, 1
				summ=summ+IM_old(i,j)*M_new(riga,i)
			END DO
			DO i = 1, N, 1
				IM_new(i,j)=IM_old(i,j)-IM_old(i,riga)*quoz*summ
			END DO
		END DO
		
	END SUBROUTINE aggiorna_matrice_inversa_C_1ppt
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_matrice_inversa_C_col_1ppt(N,col,IM_old,det_old,M_new,det_new,IM_new)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, col       !col=colonna che é stata cambiata
		COMPLEX (KIND=8), INTENT(IN) :: IM_old(1:N,1:N), M_new(1:N,1:N), det_old, det_new       !IM_old=matrice inversa old, M_new=nuova matrice
		COMPLEX (KIND=8) :: quoz, summ
		INTEGER :: i, j, l
		COMPLEX (KIND=8), INTENT(OUT) :: IM_new(1:N,1:N)

		quoz=det_old/det_new

		DO i = 1, col-1, 1
			summ=0.d0
			DO j = 1, N, 1
				summ=summ+IM_old(i,j)*M_new(j,col)
			END DO
			DO j = 1, N, 1
				IM_new(i,j)=IM_old(i,j)-IM_old(col,j)*quoz*summ
			END DO
		END DO

		DO j = 1, N, 1
			IM_new(col,j)=IM_old(col,j)*quoz
		END DO

		DO i = col+1, N, 1
			summ=0.d0
			DO j = 1, N, 1
				summ=summ+IM_old(i,j)*M_new(j,col)
			END DO
			DO j = 1, N, 1
				IM_new(i,j)=IM_old(i,j)-IM_old(col,j)*quoz*summ
			END DO
		END DO

	END SUBROUTINE aggiorna_matrice_inversa_C_col_1ppt
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_matrice_inversa_R_1ppt(N,riga,IM_old,det_old,M_new,det_new,IM_new)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, riga       !riga=riga che é stata cambiata
		REAL (KIND=8), INTENT(IN) :: IM_old(1:N,1:N), M_new(1:N,1:N), det_old, det_new       !IM_old=matrice inversa old, M_new=nuova matrice
		REAL (KIND=8) :: quoz, summ
		INTEGER :: i, j, l
		REAL (KIND=8), INTENT(OUT) :: IM_new(1:N,1:N)
		
		quoz=det_old/det_new
		
		DO j = 1, riga-1, 1
			summ=0.d0
			DO i = 1, N, 1
				summ=summ+IM_old(i,j)*M_new(riga,i)
			END DO
			DO i = 1, N, 1
				IM_new(i,j)=IM_old(i,j)-IM_old(i,riga)*quoz*summ
			END DO
		END DO
		
		DO i = 1, N, 1
			IM_new(i,riga)=IM_old(i,riga)*quoz
		END DO
		
		DO j = riga+1, N, 1
			summ=0.d0
			DO i = 1, N, 1
				summ=summ+IM_old(i,j)*M_new(riga,i)
			END DO
			DO i = 1, N, 1
				IM_new(i,j)=IM_old(i,j)-IM_old(i,riga)*quoz*summ
			END DO
		END DO

	END SUBROUTINE aggiorna_matrice_inversa_R_1ppt
!-----------------------------------------------------------------------
	SUBROUTINE aggiorna_matrice_inversa_R_col_1ppt(N,col,IM_old,det_old,M_new,det_new,IM_new)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N, col       !col=colonna che é stata cambiata
		REAL (KIND=8), INTENT(IN) :: IM_old(1:N,1:N), M_new(1:N,1:N), det_old, det_new       !IM_old=matrice inversa old, M_new=nuova matrice
		REAL (KIND=8) :: quoz, summ
		INTEGER :: i, j, l
		REAL (KIND=8), INTENT(OUT) :: IM_new(1:N,1:N)
		
		quoz=det_old/det_new
		
		DO i = 1, col-1, 1
			summ=0.d0
			DO j = 1, N, 1
				summ=summ+IM_old(i,j)*M_new(j,col)
			END DO
			DO j = 1, N, 1
				IM_new(i,j)=IM_old(i,j)-IM_old(col,j)*quoz*summ
			END DO
		END DO
		
		DO j = 1, N, 1
			IM_new(col,j)=IM_old(col,j)*quoz
		END DO
		
		DO i = col+1, N, 1
			summ=0.d0
			DO j = 1, N, 1
				summ=summ+IM_old(i,j)*M_new(j,col)
			END DO
			DO j = 1, N, 1
				IM_new(i,j)=IM_old(i,j)-IM_old(col,j)*quoz*summ
			END DO
		END DO

	END SUBROUTINE aggiorna_matrice_inversa_R_col_1ppt
!-----------------------------------------------------------------------
	SUBROUTINE salva_vettore_dati_dat(x,y,N,nome_file)
		IMPLICIT NONE
		CHARACTER (LEN=*), INTENT(IN) :: nome_file
		INTEGER (KIND=8), INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN) :: x(1:N), y(1:N)
		LOGICAL, PARAMETER :: flag_debug=.FALSE.
		INTEGER (KIND=8) :: i
		
		OPEN (UNIT=77, FILE=nome_file//'.dat', STATUS='UNKNOWN')
		DO i = 1, N, 1
			WRITE (77, *), x(i), y(i)
		END DO
		CLOSE(77)
		
	END SUBROUTINE salva_vettore_dati_dat
!-----------------------------------------------------------------------
	SUBROUTINE salva_vettore_dati_dat_ridotto(x,N,N_dati_nel_file,nome_file)
		IMPLICIT NONE
		CHARACTER (LEN=*), INTENT(IN) :: nome_file
		INTEGER (KIND=8), INTENT(IN) :: N, N_dati_nel_file
		REAL (KIND=8), INTENT(IN) :: x(1:N)
		LOGICAL, PARAMETER :: flag_debug=.FALSE.
		INTEGER (KIND=8) :: i
		
		OPEN (UNIT=77, FILE=nome_file//'.dat', STATUS='UNKNOWN')
		IF (N>N_dati_nel_file) THEN
			DO i = 1, N, N/N_dati_nel_file
				WRITE (77, *), x(i)
			END DO
		ELSE
			DO i = 1, N, 1
				WRITE (77, *), x(i)
			END DO
		END IF
		CLOSE(77)
		
	END SUBROUTINE salva_vettore_dati_dat_ridotto
!-----------------------------------------------------------------------
	SUBROUTINE salva_vettore_dati_bin(F,N,N_AV,flag_continua,nome_file)
		IMPLICIT NONE
		LOGICAL, INTENT(IN) :: flag_continua
		CHARACTER (LEN=*), INTENT(IN) :: nome_file
		INTEGER, INTENT(IN) :: N_AV
		INTEGER (KIND=8) :: N
		REAL (KIND=8), INTENT(IN) :: F(1:N)
		LOGICAL, PARAMETER :: flag_debug=.FALSE.
		LOGICAL :: flag_file_dat, flag_file_w, f1, f2
		INTEGER, PARAMETER :: N_file=512  !numero di dati di tipo double che ci saranno in ogni record dei file accumulatori
		INTEGER :: N_AV_old
		INTEGER (KIND=8) :: i, j, N_old, N_resto
		REAL (KIND=8) :: F_AV(1:N/N_AV), F_AV_prima_riga(1:N_file), R_N_AV
		
		N=N-MOD(N,N_AV)
		R_N_AV=REAL(N_AV,8)
		INQUIRE(FILE=nome_file//'.bin',EXIST=flag_file_dat)
		IF (flag_continua) THEN
			IF (.NOT. flag_file_dat) PRINT *, 'Hai detto di voler continuare a salvare dai dati precedenti ma non ci sono. &
			  Trascrivo ugualmente i dati. [salva_vettore_dati_bin]: ', nome_file
		END IF		
		
		IF ((.NOT. flag_continua) .AND. (.NOT. flag_file_dat)) THEN
			OPEN (77, FILE=nome_file//'.bin', STATUS='NEW',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		ELSE IF ((.NOT. flag_continua) .AND. (flag_file_dat)) THEN
			OPEN (77, FILE=nome_file//'.bin', STATUS='REPLACE',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		ELSE IF ((flag_continua) .AND. (.NOT. flag_file_dat)) THEN
			OPEN (77, FILE=nome_file//'.bin', STATUS='NEW',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		ELSE IF ((flag_continua) .AND. (flag_file_dat)) THEN
			OPEN (77, FILE=nome_file//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		END IF

		F_AV=0.d0
		DO i = 1, N/N_AV, 1
			DO j = 1, N_AV, 1
				F_AV(i)=F_AV(i)+F((i-1)*N_AV+j)
			END DO
			F_AV(i)=F_AV(i)/R_N_AV
		END DO

		IF ( (.NOT. flag_continua) .OR. (.NOT. flag_file_dat) ) THEN
			!scrivo un nuovo file dati
			IF (N_AV<1) STOP 'Non puoi salvare i dati se non setti un N_AV adeguato &
			  [ module_generic_tools.f90 > salva_vettore_dati_bin ]'
			IF (N<1) STOP 'Non hai dati da salvare &
			  [ module_generic_tools.f90 > salva_vettore_dati_bin ]'

			WRITE(77, REC=1), N, N_AV

			DO i = 1, (N/N_AV)/N_file+1, 1
				IF ((i-1)*N_file+1 > N/N_AV) EXIT
				WRITE(77, REC=i+1), F_AV((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
				IF (flag_debug) PRINT * , 'ho scritto i dati sulla riga ', i, ' da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
			END DO
		ELSE
			!continuo da un file dati vecchio
			READ(77, REC=1), N_old, N_AV_old
			IF (flag_debug) PRINT * , 'Leggo: N_old = ', N_old

			IF (N_AV/=N_AV_old) STOP 'Stai provando ad aggiunger dati ma con un N_AV diverso &
			  [ generic_tools.f90 > salva_vettore_dati_bin ]'

			WRITE (77, REC=1), N_old+N, N_AV
			IF (flag_debug) PRINT * , 'Trascrivo: N = ', N_old+N

			IF (MOD(N_old/N_AV_old,N_file)==0) THEN
				!nel caso il file vecchio terminasse senza una linea lasciata a metá
				DO i = (N_old/N_AV_old)/N_file+1, (N_old/N_AV_old)/N_file+(N/N_AV)/N_file+1, 1
					IF ((i-1)*N_file+1 > N/N_AV+N_old/N_AV_old) EXIT
					WRITE(77, REC=i+1), F_AV((i-1)*N_file-N_old/N_AV_old+1:MIN(i*N_file-N_old/N_AV_old,N/N_AV))
					IF (flag_debug) PRINT * , 'ho scritto i dati sulla riga ', i, ' da ', (i-1)*N_file-N_old/N_AV_old+1, &
					  ' a ', MIN(i*N_file-N_old/N_AV_old,N/N_AV)
				END DO
			ELSE
				!nel caso il file vecchio avesse una linea lasciata a metá
				N_resto=N_old/N_AV_old-((N_old/N_AV_old)/N_file)*N_file
				IF (flag_debug) PRINT * , 'N_resto= ', N_resto
				READ(77, REC=(N_old/N_AV)/N_file+2), F_AV_prima_riga(1:N_resto)
				IF (flag_debug) PRINT * , 'ho letto dalla riga ', (N_old/N_AV)/N_file+1, ' i vecchi dati da 1 a ', N_resto
				F_AV_prima_riga(N_resto+1:MIN(N_file,N/N_AV+N_resto))=F_AV(1:MIN(N_file-N_resto,N/N_AV))
				IF (flag_debug) PRINT * , 'ai primi ', N_resto, ' dati, aggiungo i dati da 1 a ', MIN(N_file-N_resto,N/N_AV)
				WRITE(77, REC=(N_old/N_AV)/N_file+2), F_AV_prima_riga(1:MIN(N_file,N/N_AV+N_resto))
				IF (flag_debug) PRINT * , 'ho ri-scritto sulla riga ', (N_old/N_AV)/N_file+1, ' da 1 a ', MIN(N_file,N/N_AV+N_resto)
				DO i = (N_old/N_AV_old)/N_file+2, (N_old/N_AV_old)/N_file+(N/N_AV)/N_file+2, 1
					IF ((i-1)*N_file+1 > N/N_AV+N_old/N_AV_old) EXIT
					WRITE(77, REC=i+1), F_AV((i-1)*N_file-N_old/N_AV_old+1:MIN(i*N_file-N_old/N_AV_old,N/N_AV))
					IF (flag_debug) PRINT * , 'ho scritto i dati sulla riga ', i, ' da ', (i-1)*N_file-N_old/N_AV_old+1, &
					  ' a ', MIN(i*N_file-N_old/N_AV_old,N/N_AV)
				END DO
			END IF

		END IF

		CLOSE(77,STATUS='KEEP')
	END SUBROUTINE salva_vettore_dati_bin
!-----------------------------------------------------------------------
	!comprime un file dati .bin, N_AV é il numero di dati che verranno compattati in uno
	SUBROUTINE comprimi_dati_bin(N_AV,nome_file)
		IMPLICIT NONE
		CHARACTER (LEN=*), INTENT(IN) :: nome_file
		INTEGER, INTENT(IN) :: N_AV
		LOGICAL, PARAMETER :: flag_debug=.FALSE.
		INTEGER, PARAMETER :: N_file=512  !numero di dati di tipo double che ci saranno in ogni record dei file accumulatori
		LOGICAL :: flag_file_dat
		INTEGER :: N_AV_old
		INTEGER (KIND=8) :: i, j, N_old
		REAL (KIND=8) :: R_N_AV
		REAL (KIND=8), ALLOCATABLE :: dati(:), dati_new(:)

		R_N_AV=REAL(N_AV,8)

		INQUIRE(FILE=nome_file//'.bin',EXIST=flag_file_dat)
		IF (.NOT. flag_file_dat) STOP 'Manca il file da comprimere &
		  [comprimi_dati_bin] '
		OPEN (77, FILE=nome_file//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)

		READ(77, REC=1), N_old, N_AV_old
		ALLOCATE(dati(1:N_old/N_AV_old), dati_new(1:(N_old/N_AV_old)/N_AV))
		DO i = 1, (N_old/N_AV_old)/N_file+1, 1
			IF ((i-1)*N_file+1>N_old/N_AV_old) EXIT
			IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N_old/N_AV_old)
			READ (77, REC=i+1), dati((i-1)*N_file+1:MIN(i*N_file,N_old/N_AV_old))
		END DO

		CLOSE (77, STATUS='KEEP')

		DO i = 1, (N_old/N_AV_old)/N_AV, 1
			dati_new(i)=dati((i-1)*N_AV+1)
			IF (flag_debug) PRINT * , 'inizio: ', i, ' parte da elemento ', (i-1)*N_AV+1, dati((i-1)*N_AV+1)
			DO j = (i-1)*N_AV+2, i*N_AV, 1
				dati_new(i)=dati_new(i)+dati(j)
				IF (flag_debug) PRINT * , 'continua: a ', i, ' si aggiunge elemento ', j, dati(j)
			END DO
			dati_new(i)=dati_new(i)/R_N_AV
		END DO

		OPEN (77, FILE=nome_file//'.bin', STATUS='REPLACE',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		WRITE(77, REC=1), N_old-MOD(N_old,N_AV), N_AV*N_AV_old
		DO i = 1, (N_old/N_AV*N_AV_old)/N_file+1, 1
			IF ((i-1)*N_file+1 > N_old/N_AV) EXIT
			WRITE(77, REC=i+1), dati_new((i-1)*N_file+1:MIN(i*N_file,N_old/(N_AV*N_AV_old)))
			IF (flag_debug) PRINT * , 'ho scritto i dati sulla riga ', i, ' da ', (i-1)*N_file+1, &
			  ' a ', MIN(i*N_file,N_old/(N_AV*N_AV_old))
		END DO
		CLOSE(77,STATUS='KEEP')
		DEALLOCATE(dati, dati_new)
	END SUBROUTINE comprimi_dati_bin
!-----------------------------------------------------------------------
	SUBROUTINE trasforma_bin_in_dat(nome_file,N_rid)
		IMPLICIT NONE
		CHARACTER (LEN=*), INTENT(IN) :: nome_file
		INTEGER, INTENT(IN) :: N_rid
		LOGICAL, PARAMETER :: flag_debug=.FALSE.
		INTEGER, PARAMETER :: N_file=512  !numero di dati di tipo double che ci saranno in ogni record dei file accumulatori
		INTEGER :: N_AV
		INTEGER (KIND=8) :: N, i, j
		REAL (KIND=8) :: frf
		REAL (KIND=8), ALLOCATABLE :: dati(:)

		OPEN (77, FILE=nome_file//'.bin', STATUS='UNKNOWN',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		OPEN (UNIT=78, FILE=nome_file//'.dat', STATUS='UNKNOWN')
		READ (77, REC=1), N, N_AV
		IF (flag_debug) PRINT * , 'Il file da trsformare ha N=',N,'    N_AV=', N_AV
		ALLOCATE(dati(1:N/N_AV))
		DO i = 1, (N/N_AV)/N_file+1, 1
			IF ((i-1)*N_file+1>N/N_AV) EXIT
			IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
			READ (77, REC=i+1), dati((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
		END DO
		DO i = 1, (N/N_AV)/N_rid, 1
			frf=0.d0
			DO j = 1, N_rid, 1
				frf=frf+dati((i-1)*N_rid+j)/REAL(N_rid,8)
			END DO
			WRITE (78, *), frf
		END DO
		CLOSE(77,STATUS='KEEP')
		CLOSE(78)
		DEALLOCATE(dati)
	END SUBROUTINE trasforma_bin_in_dat
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatori(nome_file,media,errore,nome_file_pesi)
		IMPLICIT NONE
		CHARACTER (LEN=*), INTENT(IN) :: nome_file           !file con il numeratore
		CHARACTER (LEN=*), OPTIONAL :: nome_file_pesi        !file con il denominatore
		INTEGER, PARAMETER :: num_block_max=50, num_block_min=6, num_media_max=4
		INTEGER :: n_media
		INTEGER :: i, j, i_min
		REAL (KIND=8) :: F(1:8), Q(1:8)
		REAL (KIND=8) :: average, error, delta, delta_(num_block_min:num_block_max)
		REAL (KIND=8) :: error_(num_block_min:num_block_max), average_(num_block_min:num_block_max)
		REAL (KIND=8) :: media, errore
		
		OPEN (UNIT=40, FILE=nome_file//'_details.d', STATUS='UNKNOWN')
		WRITE (40, *), '                                  - - - ESTIMATION FOR ',nome_file,' - - -'
		WRITE (40, *), 
		IF (PRESENT(nome_file_pesi)) THEN
			i=-1
			CALL calcola_estimatore_con_blocking_da_file_bin(nome_file_pesi,i,average,error)
			WRITE (40, *), 'Denominator: ', average ,' +- ', error*error, &
			  ' . Ratio: ', average/error
			WRITE (40, *),
		END IF
		DO i = num_block_min, num_block_max, 1
			IF (PRESENT(nome_file_pesi)) THEN
				CALL calcola_estimatore_con_blocking_da_file_bin(nome_file,i,average_(i),error_(i),nome_file_pesi)
			ELSE
				CALL calcola_estimatore_con_blocking_da_file_bin(nome_file,i,average_(i),error_(i))
			END IF
			WRITE (40, *), i, ' blocks: ', average_(i), ' +- ', error_(i)
		END DO
		WRITE (40, *), 
		IF (PRESENT(nome_file_pesi)) THEN 
			CALL estrai_valori_con_blocking_da_file_bin(nome_file,8,F,nome_file_pesi,Q)
			WRITE (40, *), 'NUM with 8 blocks: ', REAL(F,4)
			WRITE (40, *), 'DEN with 8 blocks: ', REAL(Q,4)
			WRITE (40, *), 'N/D with 8 blocks: ', REAL((/(F(i)/Q(i),i=1,8)/),4)
			WRITE (40, *), 
		END IF
		i_min=8
		DO i = 9, 16, 1
			IF (error_(i)>error_(i_min)) i_min=i
		END DO
		WRITE (40, *), '### CON METODO KALOS-PEDERIVA: ', i_min, ' BLOCCHI ###'
		WRITE (40, *), nome_file, ': ', average_(i_min), '+-', error_(i_min)
		delta_=0.d0
		i_min=num_block_min+num_media_max
		DO n_media = 1, num_media_max, 1
			DO i = num_block_min+num_media_max, num_block_max-num_media_max, 1
				IF (n_media==1) delta=(-0.5d0*error_(i-1)+0.5d0*error_(i+1))
				IF (n_media==2) delta=((1.d0/12.d0)*error_(i-2)-(2.d0/3.d0)*error_(i-1)+ &
				  (2.d0/3.d0)*error_(i+1)-(1.d0/12.d0)*error_(i+2))
				IF (n_media==3) delta=(-(1.d0/60.d0)*error_(i-3)+(3.d0/20.d0)*error_(i-2)-0.75d0*error_(i-1)+ &
				  0.75d0*error_(i+1)-(3.d0/20.d0)*error_(i+2)+(1.d0/60.d0)*error_(i+3))
				IF (n_media==4) delta=((1.d0/280.d0)*error_(i-4)-(4.d0/105.d0)*error_(i-3)+0.2d0*error_(i-2)- &
				  0.8d0*error_(i-1)+0.8d0*error_(i+1)-0.2d0*error_(i+2)+(4.d0/105.d0)*error_(i+3)-(1.d0/280.d0)*error_(i+4))
				delta_(i)=delta_(i)+delta
				IF (n_media==num_media_max) THEN
					IF (DABS(delta_(i))<DABS(delta_(i_min))) i_min=i
				END IF
			END DO
		END DO
		WRITE (40, *), '### CON METODO PLATEAU: ', i_min, ' BLOCCHI ###'
		average=0.2d0*(average_(i_min-2)+average_(i_min-1)+average_(i_min)+average_(i_min+1)+average_(i_min+2))
		error=0.2d0*(error_(i_min-2)+error_(i_min-1)+error_(i_min)+error_(i_min+1)+error_(i_min+2))
		WRITE (40, *), nome_file, ': ', average, ' +- ', error
		WRITE (40,*), 'errore overstimato:', MAXVAL(error_(i_min-2:i_min+2))
		WRITE (40, *), 
		WRITE (40, *), 
		media=average
		errore=error
		CLOSE (40)
	END SUBROUTINE calcola_estimatori
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatori_quick(num_blocchi,nome_file,media,errore,nome_file_pesi)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: num_blocchi
		CHARACTER (LEN=*), INTENT(IN) :: nome_file           !file con il numeratore
		CHARACTER (LEN=*), OPTIONAL :: nome_file_pesi        !file con il denominatore
		REAL (KIND=8) :: media, errore
		
		IF (PRESENT(nome_file_pesi)) THEN
			CALL calcola_estimatore_con_blocking_da_file_bin(nome_file,num_blocchi,media,errore,nome_file_pesi)
		ELSE
			CALL calcola_estimatore_con_blocking_da_file_bin(nome_file,num_blocchi,media,errore)
		END IF
		
	END SUBROUTINE calcola_estimatori_quick
!-----------------------------------------------------------------------
	SUBROUTINE estrai_valori_con_blocking_da_file_bin(nome_file,num_block,F,nome_file_pesi,Q)
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: nome_file
		CHARACTER(LEN=*), OPTIONAL :: nome_file_pesi
		INTEGER :: num_block        !se num_block=-1 allora num_block=N/N_AV (cioé ogni blocco é formato da un dato)
		LOGICAL, PARAMETER :: flag_debug=.FALSE.
		INTEGER, PARAMETER :: N_file=512  !numero di dati di tipo double che ci saranno in ogni record dei file accumulatori
		LOGICAL :: flag_file_dat, flag_file_dat_pesi
		INTEGER :: N_AV
		INTEGER (KIND=8) :: N, i, j, index
		REAL (KIND=8) :: dato_block, peso_block
		REAL (KIND=8), ALLOCATABLE :: numeratore(:), denominatore(:)
		REAL (KIND=8) :: F(1:num_block)
		REAL (KIND=8), OPTIONAL :: Q(1:num_block)     !numeratore e (eventualmente) denominatore per ogni blocco

		INQUIRE(FILE=nome_file//'.bin',EXIST=flag_file_dat)
		IF (PRESENT(nome_file_pesi)) INQUIRE(FILE=nome_file_pesi//'.bin',EXIST=flag_file_dat_pesi)
		flag_file_dat=flag_file_dat.AND.flag_file_dat_pesi
		IF (.NOT. flag_file_dat) THEN
			PRINT *, nome_file_pesi
		  	STOP 'Manca il file da cui prendere i dati per calcolare il valore di aspettazione &
		  	 [estrai_valori_con_blocking_da_file_bin] '
		END IF

		OPEN (77, FILE=nome_file//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		READ (77, REC=1), N, N_AV
		IF (PRESENT(nome_file_pesi)) OPEN (78,FILE=nome_file_pesi//'.bin',STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED',RECL=N_file*8)

		IF (num_block==-1) num_block=N/N_AV
		IF (flag_debug) PRINT * , 'Il file da cui calcolare la media ha N=', N,'    N_AV=', N_AV
		
		ALLOCATE(numeratore(1:N/N_AV))
		DO i = 1, (N/N_AV)/N_file+1, 1
			IF ((i-1)*N_file+1>N/N_AV) EXIT
			IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
			READ (77, REC=i+1), numeratore((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
		END DO
		CLOSE(77,STATUS='KEEP')
		IF (PRESENT(nome_file_pesi)) THEN
			ALLOCATE(denominatore(1:N/N_AV))
			DO i = 1, (N/N_AV)/N_file+1, 1
				IF ((i-1)*N_file+1>N/N_AV) EXIT
				IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
				READ (78, REC=i+1), denominatore((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
			END DO
			CLOSE(78,STATUS='KEEP')
		END IF

		DO i = 1, num_block, 1
			dato_block=0.d0
			index=(i-1)*((N/N_AV)/num_block)
			DO j = 1, (N/N_AV)/num_block, 1
				dato_block=dato_block+numeratore(index+j)
			END DO
			F(i)=dato_block
			IF (flag_debug) PRINT * , 'Ho calcolato il dato per il blocco ', i, ' con i dati da ', &
			  (i-1)*((N/N_AV)/num_block)+1, ' a ', i*((N/N_AV)/num_block)
		END DO
		IF (PRESENT(nome_file_pesi)) THEN
			DO i = 1, num_block, 1
				peso_block=0.d0
				index=(i-1)*((N/N_AV)/num_block)
				DO j = 1, (N/N_AV)/num_block, 1
					peso_block=peso_block+denominatore(index+j)
				END DO
				Q(i)=peso_block
				IF (flag_debug) PRINT * , 'Ho calcolato il dato per il blocco ', i, ' con i dati da ', &
				  (i-1)*((N/N_AV)/num_block)+1, ' a ', i*((N/N_AV)/num_block)
			END DO
		END IF

		DEALLOCATE(numeratore)
		IF (PRESENT(nome_file_pesi)) DEALLOCATE(denominatore)

	END SUBROUTINE estrai_valori_con_blocking_da_file_bin
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_con_blocking_da_file_bin(nome_file,num_block,media,errore,nome_file_pesi)
		IMPLICIT NONE
		CHARACTER (LEN=*), INTENT(IN) :: nome_file           !file con il numeratore
		CHARACTER (LEN=*), OPTIONAL :: nome_file_pesi        !file con il denominatore
		INTEGER :: num_block        !se num_block=-1 allora num_block=N/N_AV (cioé ogni blocco é formato da un dato)
		LOGICAL, PARAMETER :: flag_debug=.FALSE.
		INTEGER, PARAMETER :: N_file=512  !numero di dati di tipo double che ci saranno in ogni record dei file accumulatori
		LOGICAL :: flag_file_dat, flag_file_dat_pesi
		INTEGER :: N_AV
		INTEGER (KIND=8) :: N, i, j, index
		REAL (KIND=8) :: dato_block, peso_block
		REAL (KIND=8), ALLOCATABLE :: numeratore(:), denominatore(:)
		REAL (KIND=8), INTENT(OUT) :: media, errore
		
		IF (PRESENT(nome_file_pesi)) THEN
			INQUIRE(FILE=nome_file//'.bin',EXIST=flag_file_dat)
			INQUIRE(FILE=nome_file_pesi//'.bin',EXIST=flag_file_dat_pesi)
		ELSE
			INQUIRE(FILE=nome_file//'.bin',EXIST=flag_file_dat)
			flag_file_dat_pesi=.TRUE.
		END IF
		!PRINT * , '*******', nome_file
		IF ((.NOT. flag_file_dat).OR.((.NOT. flag_file_dat_pesi) .AND. (PRESENT(nome_file_pesi)))) THEN
                        IF (.NOT. flag_file_dat) PRINT *, nome_file//'.bin'
			IF ((.NOT. flag_file_dat_pesi) .AND. (PRESENT(nome_file_pesi))) PRINT *, nome_file_pesi//'.bin'
                        STOP 'Manca il file da cui prendere i dati per calcolare il valore di aspettazione &
                         [estrai_valori_con_blocking_da_file_bin] '
                END IF
	
	
		IF (PRESENT(nome_file_pesi)) THEN
			OPEN (77, FILE=nome_file//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
			READ (77, REC=1), N, N_AV
			OPEN (78, FILE=nome_file_pesi//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		ELSE
			OPEN (77, FILE=nome_file//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
			READ (77, REC=1), N, N_AV
		END IF
	
		IF (num_block==-1) num_block=N/N_AV
	
		IF (flag_debug) PRINT * , 'Il file da cui calcolare la media ha N=', N,'    N_AV=', N_AV
		IF (PRESENT(nome_file_pesi)) THEN
			ALLOCATE(numeratore(1:N/N_AV),denominatore(1:N/N_AV))
			DO i = 1, (N/N_AV)/N_file+1, 1
				IF ((i-1)*N_file+1>N/N_AV) EXIT
				IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
				READ (77, REC=i+1), numeratore((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
				READ (78, REC=i+1), denominatore((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
			END DO
			CLOSE(77,STATUS='KEEP')
			CLOSE(78,STATUS='KEEP')
		ELSE
			ALLOCATE(numeratore(1:N/N_AV))
			DO i = 1, (N/N_AV)/N_file+1, 1
				IF ((i-1)*N_file+1>N/N_AV) EXIT
				IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
				READ (77, REC=i+1), numeratore((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
			END DO
			CLOSE(77,STATUS='KEEP')
		END IF
	
		media=0.d0
		errore=0.d0
		IF (PRESENT(nome_file_pesi)) THEN
			DO i = 1, num_block, 1
				dato_block=0.d0
				peso_block=0.d0
				index=(i-1)*((N/N_AV)/num_block)
				DO j = 1, (N/N_AV)/num_block, 1
					dato_block=dato_block+numeratore(index+j)
					peso_block=peso_block+denominatore(index+j)
				END DO
				IF (flag_debug) PRINT * , 'Ho calcolato il dato per il blocco ', i, ' con i dati da ', &
				  (i-1)*((N/N_AV)/num_block)+1, ' a ', i*((N/N_AV)/num_block)
				dato_block=dato_block/peso_block
				media=media+dato_block
				errore=errore+dato_block*dato_block
			END DO
			media=media/REAL(num_block,8)
			errore=DSQRT((errore/REAL(num_block,8)-media*media)/(REAL(num_block,8)-1.d0))
		ELSE
			DO i = 1, num_block, 1
				dato_block=0.d0
				index=(i-1)*((N/N_AV)/num_block)
				DO j = 1, (N/N_AV)/num_block, 1
					dato_block=dato_block+numeratore(index+j)
				END DO
				IF (flag_debug) PRINT * , 'Ho calcolato il dato per il blocco ', i, ' con i dati da ', &
				  (i-1)*((N/N_AV)/num_block)+1, ' a ', i*((N/N_AV)/num_block)
				dato_block=dato_block/(REAL((N/N_AV)/num_block,8))
				media=media+dato_block
				errore=errore+dato_block*dato_block
			END DO
			media=media/REAL(num_block,8)
			errore=DSQRT((errore/REAL(num_block,8)-media*media)/(REAL(num_block,8)-1.d0))
		END IF
	
		DEALLOCATE(numeratore)
		IF (PRESENT(nome_file_pesi)) DEALLOCATE(denominatore)
	
	END SUBROUTINE calcola_estimatore_con_blocking_da_file_bin
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_da_RAM(stima,dati,N_cpu,dati_pesi)
		IMPLICIT NONE
		LOGICAL, PARAMETER :: flag_report=.FALSE.
		INTEGER, PARAMETER :: num_block_max=50, num_block_min=6, num_media_max=4
		INTEGER (KIND=8), INTENT(IN) :: N_cpu
		INTEGER :: n_media
		INTEGER :: i
		INTEGER (KIND=8) :: j, i_min
		REAL (KIND=8) :: F(1:8), Q(1:8)
		REAL (KIND=8), INTENT(IN) :: dati(1:N_cpu)    !dati da cui calcolare l'estimatore
		REAL (KIND=8), OPTIONAL :: dati_pesi(1:N_cpu)
		REAL (KIND=8) :: average, error, delta, delta_(num_block_min:num_block_max)
		REAL (KIND=8) :: error_(num_block_min:num_block_max), average_(num_block_min:num_block_max)
		REAL (KIND=8) :: media, errore
		REAL (KIND=8), INTENT(OUT) :: stima(1:2)
		
		IF (flag_report) OPEN (UNIT=40, FILE='report.d', STATUS='UNKNOWN')
		IF (flag_report) WRITE (40, *), '                                  - - - ESTIMATION - - -'
		IF (flag_report) WRITE (40, *), 
		DO i = num_block_min, num_block_max, 1
			IF (PRESENT(dati_pesi)) THEN
				CALL calcola_estimatore_blocco_con_pesi_RAM(dati,dati_pesi,N_cpu,i,average_(i),error_(i))
			ELSE
				CALL calcola_estimatore_blocco_RAM(dati,N_cpu,i,average_(i),error_(i))
			END IF
			IF (flag_report) WRITE (40, *), i, ' blocks: ', average_(i), ' +- ', error_(i)
		END DO
		IF (flag_report) WRITE (40, *), 
		IF (flag_report) i=8
		IF (flag_report) CALL ottieni_dati_blocco_RAM(dati,N_cpu,i,F)
		IF (flag_report) WRITE (40, *), 'NUM with 8 blocks: ', REAL(F,4)
		IF (PRESENT(dati_pesi) .AND. flag_report) THEN
			CALL ottieni_dati_blocco_RAM(dati_pesi,N_cpu,i,Q)
			WRITE (40, *), 'DEN with 8 blocks: ', REAL(Q,4)
			WRITE (40, *), 'N/D with 8 blocks: ', REAL((/(F(i)/Q(i),i=1,8)/),4)
		END IF
		IF (flag_report) WRITE (40, *), 
		i_min=8
		DO i = 9, 16, 1
			IF (error_(i)>error_(i_min)) i_min=i
		END DO
		IF (flag_report) WRITE (40, *), '### CON METODO KALOS-PEDERIVA: ', i_min, ' BLOCCHI ###'
		IF (flag_report) WRITE (40, *), 'ESTIMATOR: ', average_(i_min), '+-', error_(i_min)
		delta_=0.d0
		i_min=num_block_min+num_media_max
		DO n_media = 1, num_media_max, 1
			DO i = num_block_min+num_media_max, num_block_max-num_media_max, 1
				IF (n_media==1) delta=(-0.5d0*error_(i-1)+0.5d0*error_(i+1))
				IF (n_media==2) delta=((1.d0/12.d0)*error_(i-2)-(2.d0/3.d0)*error_(i-1)+ &
				  (2.d0/3.d0)*error_(i+1)-(1.d0/12.d0)*error_(i+2))
				IF (n_media==3) delta=(-(1.d0/60.d0)*error_(i-3)+(3.d0/20.d0)*error_(i-2)-0.75d0*error_(i-1)+ &
				  0.75d0*error_(i+1)-(3.d0/20.d0)*error_(i+2)+(1.d0/60.d0)*error_(i+3))
				IF (n_media==4) delta=((1.d0/280.d0)*error_(i-4)-(4.d0/105.d0)*error_(i-3)+0.2d0*error_(i-2)- &
				  0.8d0*error_(i-1)+0.8d0*error_(i+1)-0.2d0*error_(i+2)+(4.d0/105.d0)*error_(i+3)-(1.d0/280.d0)*error_(i+4))
				delta_(i)=delta_(i)+delta
				IF (n_media==num_media_max) THEN
					IF (DABS(delta_(i))<DABS(delta_(i_min))) i_min=i
				END IF
			END DO
		END DO
		IF (flag_report) WRITE (40, *), '### CON METODO PLATEAU: ', i_min, ' BLOCCHI ###'
		average=0.2d0*(average_(i_min-2)+average_(i_min-1)+average_(i_min)+average_(i_min+1)+average_(i_min+2))
		error=0.2d0*(error_(i_min-2)+error_(i_min-1)+error_(i_min)+error_(i_min+1)+error_(i_min+2))
		IF (flag_report) WRITE (40, *), 'ESTIMATOR: ', average, ' +- ', error
		IF (flag_report) WRITE (40,*), 'errore overstimato:', MAXVAL(error_(i_min-2:i_min+2))
		IF (flag_report) WRITE (40, *), 
		IF (flag_report) WRITE (40, *), 
		media=average
		errore=error
		IF (flag_report) CLOSE (40)
		stima(1)=media
		stima(2)=error
	END SUBROUTINE calcola_estimatore_da_RAM
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_da_RAM_quick(num_blocchi,stima,dati,N_cpu,dati_pesi)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: num_blocchi
		INTEGER (KIND=8), INTENT(IN) :: N_cpu
		REAL (KIND=8), INTENT(IN) :: dati(1:N_cpu)    !dati da cui calcolare l'estimatore
		REAL (KIND=8), OPTIONAL :: dati_pesi(1:N_cpu)
		REAL (KIND=8), INTENT(OUT) :: stima(1:2)
		
		IF (PRESENT(dati_pesi)) THEN
			CALL calcola_estimatore_blocco_con_pesi_RAM(dati,dati_pesi,N_cpu,num_blocchi,stima(1),stima(2))
		ELSE
			CALL calcola_estimatore_blocco_RAM(dati,N_cpu,num_blocchi,stima(1),stima(2))
		END IF
		
	END SUBROUTINE calcola_estimatore_da_RAM_quick
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_blocco_RAM(dati,N_cpu,num_blocks,media,errore)
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER (KIND=8), INTENT(IN) :: N_cpu     !numero di dati in ogni CPU
		INTEGER, INTENT(IN) :: num_blocks    !numero di blocchi
		REAL (KIND=8), INTENT(IN) :: dati(1:N_cpu)    !dati da cui calcolare l'estimatore
		INTEGER :: mpi_ierr, mpi_myrank, mpi_nprocs
		INTEGER (KIND=8) :: N_tot, N_block            !numero totale dati, numero di dati per ogni blocco
		INTEGER (KIND=8) :: i, j, index
		REAL (KIND=8) :: data_block(1:num_blocks), tot_data_block(1:num_blocks)
		REAL (KIND=8), INTENT(OUT) :: media, errore
		
		CALL ottieni_dati_blocco_RAM(dati,N_cpu,num_blocks,tot_data_block)
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nprocs, mpi_ierr)
		N_block=N_cpu*mpi_nprocs/num_blocks
		media=0.d0
		errore=0.d0
		DO i = 1, num_blocks, 1
			media=media+tot_data_block(i)
			errore=errore+tot_data_block(i)**2
		END DO
		media=media/REAL(num_blocks,8)
		errore=errore/REAL(num_blocks,8)
		errore=DSQRT((errore-media*media)/(REAL(num_blocks,8)-1.d0))
		media=media/REAL(N_block,8)
		errore=errore/REAL(N_block,8)
		
	END SUBROUTINE calcola_estimatore_blocco_RAM
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_blocco_con_pesi_RAM(dati,dati_peso,N_cpu,num_blocks,media,errore)
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER (KIND=8), INTENT(IN) :: N_cpu     !numero di dati in ogni CPU
		INTEGER, INTENT(IN) :: num_blocks    !numero di blocchi
		REAL (KIND=8), INTENT(IN) :: dati(1:N_cpu), dati_peso(1:N_cpu)    !dati da cui calcolare l'estimatore
		INTEGER :: mpi_ierr, mpi_myrank, mpi_nprocs
		INTEGER (KIND=8) :: N_tot, N_block            !numero totale dati, numero di dati per ogni blocco
		INTEGER (KIND=8), ALLOCATABLE :: vec_cpu(:), vec_blocks(:)      !vettori che conterranno gli indici da cui partono i dati delle cpu e dei blocchi
		INTEGER (KIND=8) :: i, j, index
		REAL (KIND=8) :: app(1:N_cpu), data_block(1:num_blocks), num(1:num_blocks), den(1:num_blocks)
		REAL (KIND=8), INTENT(OUT) :: media, errore
		
		DO i = 1, N_cpu, 1
			app(i)=dati(i)*dati_peso(i)
		END DO
		CALL ottieni_dati_blocco_RAM(app,N_cpu,num_blocks,num)
		CALL ottieni_dati_blocco_RAM(dati_peso,N_cpu,num_blocks,den)
		DO i = 1, num_blocks, 1
			data_block(i)=num(i)/den(i)
		END DO
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nprocs, mpi_ierr)
		N_block=N_cpu*mpi_nprocs/num_blocks
		media=0.d0
		errore=0.d0
		DO i = 1, num_blocks, 1
			media=media+data_block(i)
			errore=errore+data_block(i)**2
		END DO
		media=media/REAL(num_blocks,8)
		errore=errore/REAL(num_blocks,8)
		errore=DSQRT((errore-media*media)/(REAL(num_blocks,8)-1.d0))
		media=media!/REAL(N_block,8)
		errore=errore!/REAL(N_block,8)

	END SUBROUTINE calcola_estimatore_blocco_con_pesi_RAM
!-----------------------------------------------------------------------
	SUBROUTINE ottieni_dati_blocco_RAM(dati,N_cpu,num_blocks,dati_blocchi)
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER (KIND=8), INTENT(IN) :: N_cpu     !numero di dati in ogni CPU
		INTEGER, INTENT(IN) :: num_blocks    !numero di blocchi
		REAL (KIND=8), INTENT(IN) :: dati(1:N_cpu)    !dati da cui calcolare l'estimatore
		INTEGER :: mpi_ierr, mpi_myrank, mpi_nprocs
		INTEGER (KIND=8) :: N_tot, N_block            !numero totale dati, numero di dati per ogni blocco
		INTEGER (KIND=8), ALLOCATABLE :: vec_cpu(:), vec_blocks(:)      !vettori che conterranno gli indici da cui partono i dati delle cpu e dei blocchi
		INTEGER (KIND=8) :: i, j, index
		REAL (KIND=8) :: data_block(1:num_blocks)
		REAL (KIND=8), INTENT(OUT) :: dati_blocchi(1:num_blocks)
		
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nprocs, mpi_ierr)
		CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_myrank, mpi_ierr)
		ALLOCATE(vec_cpu(1:mpi_nprocs+1), vec_blocks(1:num_blocks+1))

		N_tot=mpi_nprocs*N_cpu
		N_block=N_tot/num_blocks
		vec_cpu=(/(i,i=1,N_tot,N_cpu)/)
		vec_blocks=(/(i,i=1,N_tot,N_block)/)
		data_block=0.d0
		dati_blocchi=0.d0
		DO j = 1, num_blocks, 1
			DO i = 1, N_block, 1
				index=(j-1)*N_block+i
				IF ((index>=vec_cpu(mpi_myrank+1)).AND.(index<vec_cpu(mpi_myrank+1)+N_cpu)) THEN
					data_block(j)=data_block(j)+dati(index-vec_cpu(mpi_myrank+1)+1)
				END IF
			END DO
			CALL MPI_REDUCE(data_block(j),dati_blocchi(j),1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
			CALL MPI_BCAST(dati_blocchi(j),1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_ierr)
		END DO
		
	END SUBROUTINE ottieni_dati_blocco_RAM
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatori_differenza(nome_file1,nome_file2,media,errore,nome_file_pesi1,nome_file_pesi2)
	!calcola lo stimatore di dato1-dato2
		IMPLICIT NONE
		CHARACTER (LEN=*), INTENT(IN) :: nome_file1, nome_file2, nome_file_pesi1
		CHARACTER (LEN=*), OPTIONAL :: nome_file_pesi2               !file con il denominatore
		INTEGER, PARAMETER :: num_block_max=50, num_block_min=6, num_media_max=4
		INTEGER :: n_media
		INTEGER :: i, j, i_min
		REAL (KIND=8) :: F(1:8), Q(1:8)
		REAL (KIND=8) :: average, error, delta, delta_(num_block_min:num_block_max)
		REAL (KIND=8) :: error_(num_block_min:num_block_max), average_(num_block_min:num_block_max)
		REAL (KIND=8) :: media, errore
		
		OPEN (UNIT=40, FILE=nome_file1//'_details.d', STATUS='UNKNOWN')
		WRITE (40, *), '                                  - - - ESTIMATION FOR DIFFERENCE BETWEEN ', &
		  nome_file1, ' AND ', nome_file2, ' - - -'
		WRITE (40, *), 
		DO i = num_block_min, num_block_max, 1
			IF (PRESENT(nome_file_pesi2)) THEN
				CALL calcola_estimatore_differenza_con_blocking_da_file_bin(nome_file1,nome_file2,i, &
				  average_(i),error_(i),nome_file_pesi1,nome_file_pesi2)
			ELSE
				CALL calcola_estimatore_differenza_con_blocking_da_file_bin(nome_file1,nome_file2,i, &
				  average_(i),error_(i),nome_file_pesi1)
			END IF
			WRITE (40, *), i, ' blocks: ', average_(i), ' +- ', error_(i)
		END DO
		WRITE (40, *), 
		i_min=8
		DO i = 9, 16, 1
			IF (error_(i)>error_(i_min)) i_min=i
		END DO
		WRITE (40, *), '### CON METODO KALOS-PEDERIVA: ', i_min, ' BLOCCHI ###'
		WRITE (40, *), nome_file1, ': ', average_(i_min), '+-', error_(i_min)
		delta_=0.d0
		i_min=num_block_min+num_media_max
		DO n_media = 1, num_media_max, 1
			DO i = num_block_min+num_media_max, num_block_max-num_media_max, 1
				IF (n_media==1) delta=(-0.5d0*error_(i-1)+0.5d0*error_(i+1))
				IF (n_media==2) delta=((1.d0/12.d0)*error_(i-2)-(2.d0/3.d0)*error_(i-1)+ &
				  (2.d0/3.d0)*error_(i+1)-(1.d0/12.d0)*error_(i+2))
				IF (n_media==3) delta=(-(1.d0/60.d0)*error_(i-3)+(3.d0/20.d0)*error_(i-2)-0.75d0*error_(i-1)+ &
				  0.75d0*error_(i+1)-(3.d0/20.d0)*error_(i+2)+(1.d0/60.d0)*error_(i+3))
				IF (n_media==4) delta=((1.d0/280.d0)*error_(i-4)-(4.d0/105.d0)*error_(i-3)+0.2d0*error_(i-2)- &
				  0.8d0*error_(i-1)+0.8d0*error_(i+1)-0.2d0*error_(i+2)+(4.d0/105.d0)*error_(i+3)-(1.d0/280.d0)*error_(i+4))
				delta_(i)=delta_(i)+delta
				IF (n_media==num_media_max) THEN
					IF (DABS(delta_(i))<DABS(delta_(i_min))) i_min=i
				END IF
			END DO
		END DO
		WRITE (40, *), '### CON METODO PLATEAU: ', i_min, ' BLOCCHI ###'
		average=0.2d0*(average_(i_min-2)+average_(i_min-1)+average_(i_min)+average_(i_min+1)+average_(i_min+2))
		error=0.2d0*(error_(i_min-2)+error_(i_min-1)+error_(i_min)+error_(i_min+1)+error_(i_min+2))
		WRITE (40, *), nome_file1, ': ', average, ' +- ', error
		WRITE (40,*), 'errore overstimato:', MAXVAL(error_(i_min-2:i_min+2))
		WRITE (40, *), 
		WRITE (40, *), 
		media=average
		errore=error
	END SUBROUTINE calcola_estimatori_differenza
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatori_differenza_quick(num_blocchi,nome_file1,nome_file2,media,errore,nome_file_pesi1,nome_file_pesi2)
	!calcola lo stimatore di dato1-dato2
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: num_blocchi
		CHARACTER (LEN=*), INTENT(IN) :: nome_file1, nome_file2, nome_file_pesi1
		CHARACTER (LEN=*), OPTIONAL :: nome_file_pesi2               !file con il denominatore
		REAL (KIND=8) :: media, errore
		
		IF (PRESENT(nome_file_pesi2)) THEN
			CALL calcola_estimatore_differenza_con_blocking_da_file_bin(nome_file1,nome_file2,num_blocchi, &
			  media,errore,nome_file_pesi1,nome_file_pesi2)
		ELSE
			CALL calcola_estimatore_differenza_con_blocking_da_file_bin(nome_file1,nome_file2,num_blocchi, &
			  media,errore,nome_file_pesi1)
		END IF
		
	END SUBROUTINE calcola_estimatori_differenza_quick
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_differenza_con_blocking_da_file_bin(nome_file1,nome_file2,num_block, &
	  media,errore,nome_file_pesi1,nome_file_pesi2)
		IMPLICIT NONE
		CHARACTER (LEN=*), INTENT(IN) :: nome_file1, nome_file2, nome_file_pesi1           !file con il numeratore
		CHARACTER (LEN=*), OPTIONAL :: nome_file_pesi2        !file con il denominatore
		INTEGER :: num_block        !se num_block=-1 allora num_block=N/N_AV (cioé ogni blocco é formato da un dato)
		LOGICAL, PARAMETER :: flag_debug=.FALSE.
		INTEGER, PARAMETER :: N_file=512  !numero di dati di tipo double che ci saranno in ogni record dei file accumulatori
		LOGICAL :: flag_file_dat, flag_file_dat_pesi
		INTEGER :: N_AV
		INTEGER (KIND=8) :: N, i, j, index
		REAL (KIND=8) :: dato_block1, peso_block1, dato_block2, peso_block2
		REAL (KIND=8) :: dato_block
		REAL (KIND=8), ALLOCATABLE :: numeratore1(:), denominatore1(:)
		REAL (KIND=8), ALLOCATABLE :: numeratore2(:), denominatore2(:)
		REAL (KIND=8), INTENT(OUT) :: media, errore
		
		INQUIRE(FILE=nome_file1//'.bin',EXIST=flag_file_dat)
		INQUIRE(FILE=nome_file_pesi1//'.bin',EXIST=flag_file_dat_pesi)
		IF ((.NOT. flag_file_dat).OR.(.NOT. flag_file_dat_pesi)) &
		  STOP 'Manca il file da cui prendere i dati per calcolare il valore di aspettazione &
		  [calcola_estimatore_con_blocking_da_file_bin] '
		IF (PRESENT(nome_file_pesi2)) THEN
			INQUIRE(FILE=nome_file2//'.bin',EXIST=flag_file_dat)
			INQUIRE(FILE=nome_file_pesi2//'.bin',EXIST=flag_file_dat_pesi)
		ELSE
			INQUIRE(FILE=nome_file2//'.bin',EXIST=flag_file_dat)
			flag_file_dat_pesi=.TRUE.
		END IF
		IF ((.NOT. flag_file_dat).OR.(.NOT. flag_file_dat_pesi)) &
		  STOP 'Manca il file da cui prendere i dati per calcolare il valore di aspettazione &
		  [calcola_estimatore_con_blocking_da_file_bin] '
	
		
		OPEN (75, FILE=nome_file1//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		READ (75, REC=1), N, N_AV
		OPEN (76, FILE=nome_file_pesi1//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		OPEN (77, FILE=nome_file2//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		IF (PRESENT(nome_file_pesi2)) THEN
			OPEN (78, FILE=nome_file_pesi2//'.bin', STATUS='OLD',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=N_file*8)
		END IF
	
		IF (num_block==-1) num_block=N/N_AV
	
		IF (flag_debug) PRINT * , 'Il file da cui calcolare la media ha N=', N,'    N_AV=', N_AV
		
		ALLOCATE(numeratore1(1:N/N_AV),denominatore1(1:N/N_AV))
		DO i = 1, (N/N_AV)/N_file+1, 1
			IF ((i-1)*N_file+1>N/N_AV) EXIT
			IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
			READ (75, REC=i+1), numeratore1((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
			READ (76, REC=i+1), denominatore1((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
		END DO
		CLOSE(75,STATUS='KEEP')
		CLOSE(76,STATUS='KEEP')
		
		ALLOCATE(numeratore2(1:N/N_AV),denominatore2(1:N/N_AV))
		DO i = 1, (N/N_AV)/N_file+1, 1
			IF ((i-1)*N_file+1>N/N_AV) EXIT
			IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
			READ (77, REC=i+1), numeratore2((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
		END DO
		CLOSE(77,STATUS='KEEP')
		IF (PRESENT(nome_file_pesi2)) THEN
			DO i = 1, (N/N_AV)/N_file+1, 1
				IF ((i-1)*N_file+1>N/N_AV) EXIT
				IF (flag_debug) PRINT * , 'Leggo dalla riga ', i, '  i dati da ', (i-1)*N_file+1, ' a ', MIN(i*N_file,N/N_AV)
				READ (78, REC=i+1), denominatore2((i-1)*N_file+1:MIN(i*N_file,N/N_AV))
			END DO
			CLOSE(78,STATUS='KEEP')
		ELSE
			denominatore2=N_AV
		END IF
	
		media=0.d0
		errore=0.d0
		DO i = 1, num_block, 1
			dato_block1=0.d0
			peso_block1=0.d0
			dato_block2=0.d0
			peso_block2=0.d0
			index=(i-1)*((N/N_AV)/num_block)
			DO j = 1, (N/N_AV)/num_block, 1
				dato_block1=dato_block1+numeratore1(index+j)
				peso_block1=peso_block1+denominatore1(index+j)
				dato_block2=dato_block2+numeratore2(index+j)
				peso_block2=peso_block2+denominatore2(index+j)
			END DO
			IF (flag_debug) PRINT * , 'Ho calcolato il dato per il blocco ', i, ' con i dati da ', &
			  (i-1)*((N/N_AV)/num_block)+1, ' a ', i*((N/N_AV)/num_block)
			dato_block=(dato_block1*peso_block2-dato_block2*peso_block1)/(peso_block1*peso_block2)
			media=media+dato_block
			errore=errore+dato_block*dato_block
		END DO
		media=media/REAL(num_block,8)
		errore=DSQRT((errore/REAL(num_block,8)-media*media)/(REAL(num_block,8)-1.d0))
	
		DEALLOCATE(numeratore1,numeratore2,denominatore1)
		IF (PRESENT(nome_file_pesi2)) DEALLOCATE(denominatore2)
	
	END SUBROUTINE calcola_estimatore_differenza_con_blocking_da_file_bin
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_differenza_da_RAM(stima,dati1,dati2,dati3,dati4,N_cpu)
		IMPLICIT NONE
		LOGICAL, PARAMETER :: flag_report=.FALSE.
		INTEGER, PARAMETER :: num_block_max=50, num_block_min=6, num_media_max=4
		INTEGER (KIND=8), INTENT(IN) :: N_cpu
		INTEGER :: n_media
		INTEGER :: i
		INTEGER (KIND=8) :: j, i_min
		REAL (KIND=8) :: F(1:8), Q(1:8)
		REAL (KIND=8), INTENT(IN) :: dati1(1:N_cpu),dati2(1:N_cpu),dati3(1:N_cpu),dati4(1:N_cpu)    !dati1 e dati2 riferiti a E, dati3 e dati4 a E'
		REAL (KIND=8) :: average, error, delta, delta_(num_block_min:num_block_max)
		REAL (KIND=8) :: error_(num_block_min:num_block_max), average_(num_block_min:num_block_max)
		REAL (KIND=8) :: media, errore
		REAL (KIND=8), INTENT(OUT) :: stima(1:2)
		
		IF (flag_report) OPEN (UNIT=40, FILE='report.d', STATUS='UNKNOWN')
		IF (flag_report) WRITE (40, *), '                                  - - - ESTIMATION - - -'
		IF (flag_report) WRITE (40, *), 
		DO i = num_block_min, num_block_max, 1
			CALL calcola_estimatore_differenza_blocco_RAM(dati1,dati2,dati3,dati4,N_cpu,i,average_(i),error_(i))
			IF (flag_report) WRITE (40, *), i, ' blocks: ', average_(i), ' +- ', error_(i)
		END DO
		IF (flag_report) WRITE (40, *), 
		i_min=8
		DO i = 9, 16, 1
			IF (error_(i)>error_(i_min)) i_min=i
		END DO
		IF (flag_report) WRITE (40, *), '### CON METODO KALOS-PEDERIVA: ', i_min, ' BLOCCHI ###'
		IF (flag_report) WRITE (40, *), 'ESTIMATOR: ', average_(i_min), '+-', error_(i_min)
		delta_=0.d0
		i_min=num_block_min+num_media_max
		DO n_media = 1, num_media_max, 1
			DO i = num_block_min+num_media_max, num_block_max-num_media_max, 1
				IF (n_media==1) delta=(-0.5d0*error_(i-1)+0.5d0*error_(i+1))
				IF (n_media==2) delta=((1.d0/12.d0)*error_(i-2)-(2.d0/3.d0)*error_(i-1)+ &
				  (2.d0/3.d0)*error_(i+1)-(1.d0/12.d0)*error_(i+2))
				IF (n_media==3) delta=(-(1.d0/60.d0)*error_(i-3)+(3.d0/20.d0)*error_(i-2)-0.75d0*error_(i-1)+ &
				  0.75d0*error_(i+1)-(3.d0/20.d0)*error_(i+2)+(1.d0/60.d0)*error_(i+3))
				IF (n_media==4) delta=((1.d0/280.d0)*error_(i-4)-(4.d0/105.d0)*error_(i-3)+0.2d0*error_(i-2)- &
				  0.8d0*error_(i-1)+0.8d0*error_(i+1)-0.2d0*error_(i+2)+(4.d0/105.d0)*error_(i+3)-(1.d0/280.d0)*error_(i+4))
				delta_(i)=delta_(i)+delta
				IF (n_media==num_media_max) THEN
					IF (DABS(delta_(i))<DABS(delta_(i_min))) i_min=i
				END IF
			END DO
		END DO
		IF (flag_report) WRITE (40, *), '### CON METODO PLATEAU: ', i_min, ' BLOCCHI ###'
		average=0.2d0*(average_(i_min-2)+average_(i_min-1)+average_(i_min)+average_(i_min+1)+average_(i_min+2))
		error=0.2d0*(error_(i_min-2)+error_(i_min-1)+error_(i_min)+error_(i_min+1)+error_(i_min+2))
		IF (flag_report) WRITE (40, *), 'ESTIMATOR: ', average, ' +- ', error
		IF (flag_report) WRITE (40,*), 'errore overstimato:', MAXVAL(error_(i_min-2:i_min+2))
		IF (flag_report) WRITE (40, *), 
		IF (flag_report) WRITE (40, *), 
		media=average
		errore=error
		IF (flag_report) CLOSE (40)
		stima(1)=media
		stima(2)=error
	END SUBROUTINE calcola_estimatore_differenza_da_RAM
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_differenza_da_RAM_quick(num_blocchi,stima,dati1,dati2,dati3,dati4,N_cpu)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: num_blocchi
		LOGICAL, PARAMETER :: flag_report=.FALSE.
		INTEGER, PARAMETER :: num_block_max=50, num_block_min=6, num_media_max=4
		INTEGER (KIND=8), INTENT(IN) :: N_cpu
		REAL (KIND=8), INTENT(IN) :: dati1(1:N_cpu),dati2(1:N_cpu),dati3(1:N_cpu),dati4(1:N_cpu)    !dati1 e dati2 riferiti a E, dati3 e dati4 a E'
		REAL (KIND=8), INTENT(OUT) :: stima(1:2)
		
		CALL calcola_estimatore_differenza_blocco_RAM(dati1,dati2,dati3,dati4,N_cpu,num_blocchi,stima(1),stima(2))
		
	END SUBROUTINE calcola_estimatore_differenza_da_RAM_quick
!-----------------------------------------------------------------------
	SUBROUTINE calcola_estimatore_differenza_blocco_RAM(dati1,dati2,dati3,dati4,N_cpu,num_blocks,media,errore)
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER (KIND=8), INTENT(IN) :: N_cpu     !numero di dati in ogni CPU
		INTEGER, INTENT(IN) :: num_blocks    !numero di blocchi
		REAL (KIND=8), INTENT(IN) :: dati1(1:N_cpu),dati2(1:N_cpu),dati3(1:N_cpu),dati4(1:N_cpu)     !dati da cui calcolare l'estimatore
		INTEGER :: mpi_ierr, mpi_myrank, mpi_nprocs
		INTEGER (KIND=8) :: N_tot, N_block            !numero totale dati, numero di dati per ogni blocco
		INTEGER (KIND=8) :: i, j, index
		REAL (KIND=8) :: tot_data_block1(1:num_blocks),tot_data_block2(1:num_blocks)
		REAL (KIND=8) :: tot_data_block3(1:num_blocks),tot_data_block4(1:num_blocks)
		REAL (KIND=8) :: dummy1
		REAL (KIND=8), INTENT(OUT) :: media, errore

		CALL ottieni_dati_blocco_RAM(dati1,N_cpu,num_blocks,tot_data_block1)
		CALL ottieni_dati_blocco_RAM(dati2,N_cpu,num_blocks,tot_data_block2)
		CALL ottieni_dati_blocco_RAM(dati3,N_cpu,num_blocks,tot_data_block3)
		CALL ottieni_dati_blocco_RAM(dati4,N_cpu,num_blocks,tot_data_block4)		
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nprocs, mpi_ierr)
		CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_myrank, mpi_ierr)
		N_block=N_cpu*mpi_nprocs/num_blocks
		media=0.d0
		errore=0.d0
		DO i = 1, num_blocks, 1
			dummy1=(tot_data_block1(i)*tot_data_block4(i)-tot_data_block3(i)*tot_data_block2(i)) / &
			  (tot_data_block2(i)*tot_data_block4(i))
			media=media+dummy1
			errore=errore+(dummy1)**2
		END DO
		media=media/REAL(num_blocks,8)
		errore=errore/REAL(num_blocks,8)
		errore=DSQRT((errore-media*media)/(REAL(num_blocks,8)-1.d0))
		
	END SUBROUTINE calcola_estimatore_differenza_blocco_RAM
!-----------------------------------------------------------------------
	SUBROUTINE salva_posizioni_su_file(r,N,nome_file)
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: nome_file
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(IN):: r(1:3,1:N)
		INTEGER :: i
		
		OPEN (UNIT=77, FILE=nome_file//'.pos', STATUS='UNKNOWN')
		DO i = 1, N, 1
			WRITE (77, *), r(1:3,i)
		END DO
		CLOSE(77)
		
	END SUBROUTINE salva_posizioni_su_file
!-----------------------------------------------------------------------
	SUBROUTINE leggi_posizioni_da_file(r,N,nome_file)
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: nome_file
		INTEGER, INTENT(IN) :: N
		REAL (KIND=8), INTENT(OUT):: r(1:3,1:N)
		INTEGER :: i
    
		OPEN (UNIT=77, FILE=nome_file//'.pos', STATUS='UNKNOWN')
		DO i = 1, N, 1
			READ (77, *), r(1:3,i)
		END DO
		CLOSE(77)
    
	END SUBROUTINE leggi_posizioni_da_file
!-----------------------------------------------------------------------
 !campionamento da exp(-x^2/(2*sigma^2)) di una matrice M per N
	SUBROUTINE gaussian_sample(x,M,N,sigma)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: M, N
		REAL (KIND=8), INTENT(IN) :: sigma
		INTEGER :: i, j
		REAL (KIND=8) :: xi1, xi2
		REAL (KIND=8), INTENT(OUT) :: x(M,N)

		DO j = 1, N, 1
			DO i = 1, M, 1
				CALL RANDOM_NUMBER(xi1)
				CALL RANDOM_NUMBER(xi2)
				x(i,j)=sigma*DSQRT(-2.d0*DLOG(1.d0-xi1))*DCOS(2.d0*PI*xi2)
			END DO
		END DO
	END SUBROUTINE gaussian_sample
!-----------------------------------------------------------------------
	SUBROUTINE mpi_random_seed(rank)
		INTEGER :: i, n, clock, rank
		INTEGER, DIMENSION(:), ALLOCATABLE :: seed

		CALL RANDOM_SEED(size = n)
		ALLOCATE(seed(n))

		CALL SYSTEM_CLOCK(COUNT=clock)
		
		seed = rank * 271 + 199 * (/ (i - 1, i = 1, n) /)
		CALL RANDOM_SEED(PUT = seed)

		DEALLOCATE(seed)
	END SUBROUTINE mpi_random_seed
!-----------------------------------------------------------------------
	SUBROUTINE mpi_random_seed_da_file(nome_file)
		INCLUDE 'mpif.h'
		INTEGER, PARAMETER :: num_max_colonne=100
		LOGICAL :: flag
		CHARACTER(LEN=*) :: nome_file
		INTEGER :: mpi_ierr, mpi_myrank, mpi_nprocs
		INTEGER :: n_seed
		INTEGER, DIMENSION(:), ALLOCATABLE :: seed
		INTEGER (KIND=8) :: numero_righe, numero_colonne, num_tot, i, index
		INTEGER, ALLOCATABLE :: random_numbers(:)
		REAL (KIND=8) :: dummy(1:num_max_colonne)
		
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nprocs, mpi_ierr)
		CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_myrank, mpi_ierr)
		
		CALL RANDOM_SEED(size = n_seed)
		ALLOCATE(seed(n_seed))

		IF (mpi_myrank==0) THEN
			OPEN (UNIT=37, FILE=nome_file, STATUS='OLD')
			READ (37, *), numero_righe, numero_colonne
			num_tot=numero_righe*numero_colonne
			IF (num_tot<n_seed*mpi_nprocs) THEN
            PRINT *, "Non ci sono abbastanza numeri random nel file fornito"
            PRINT *, "numero_righe = ", numero_righe
            PRINT *, "numero_colonne = ", numero_colonne
            PRINT *, "n_seed = ", n_seed
            PRINT *, "mpi_nprocs = ", mpi_nprocs
            STOP
         END IF 
			ALLOCATE(random_numbers(1:num_tot))
			DO i = 1, numero_righe, 1
				index=(i-1)*numero_colonne+1
				READ (37, *), random_numbers(index:index+numero_colonne-1)
			END DO
			CLOSE (37)
		END IF
		
		CALL MPI_SCATTER(random_numbers,n_seed,MPI_INTEGER, &
		  seed,n_seed,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)
		
		IF (mpi_myrank==0) DEALLOCATE(random_numbers)
		
		CALL RANDOM_SEED(PUT = seed)

		DEALLOCATE(seed)
	END SUBROUTINE mpi_random_seed_da_file
!-----------------------------------------------------------------------
	function alngam ( xvalue, ifault )

	!*****************************************************************************80
	!
	!! ALNGAM computes the logarithm of the gamma function.
	!
	!  Modified:
	!
	!    13 January 2008
	!
	!  Author:
	!
	!    Allan Macleod
	!    FORTRAN90 version by John Burkardt
	!
	!  Reference:
	!
	!    Allan Macleod,
	!    Algorithm AS 245,
	!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
	!    Applied Statistics,
	!    Volume 38, Number 2, 1989, pages 397-402.
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
	!
	!    Output, integer ( kind = 4 ) IFAULT, error flag.
	!    0, no error occurred.
	!    1, XVALUE is less than or equal to 0.
	!    2, XVALUE is too big.
	!
	!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
	!
	  implicit none

	  real    ( kind = 8 ) alngam
	  real    ( kind = 8 ), parameter :: alr2pi = 0.918938533204673D+00
	  integer ( kind = 4 ) ifault
	  real    ( kind = 8 ), dimension ( 9 ) :: r1 = (/ &
	    -2.66685511495D+00, &
	    -24.4387534237D+00, &
	    -21.9698958928D+00, &
	     11.1667541262D+00, &
	     3.13060547623D+00, &
	     0.607771387771D+00, &
	     11.9400905721D+00, &
	     31.4690115749D+00, &
	     15.2346874070D+00 /)
	  real    ( kind = 8 ), dimension ( 9 ) :: r2 = (/ &
	    -78.3359299449D+00, &
	    -142.046296688D+00, &
	     137.519416416D+00, &
	     78.6994924154D+00, &
	     4.16438922228D+00, &
	     47.0668766060D+00, &
	     313.399215894D+00, &
	     263.505074721D+00, &
	     43.3400022514D+00 /)
	  real    ( kind = 8 ), dimension ( 9 ) :: r3 = (/ &
	    -2.12159572323D+05, &
	     2.30661510616D+05, &
	     2.74647644705D+04, &
	    -4.02621119975D+04, &
	    -2.29660729780D+03, &
	    -1.16328495004D+05, &
	    -1.46025937511D+05, &
	    -2.42357409629D+04, &
	    -5.70691009324D+02 /)
	  real    ( kind = 8 ), dimension ( 5 ) :: r4 = (/ &
	     0.279195317918525D+00, &
	     0.4917317610505968D+00, &
	     0.0692910599291889D+00, &
	     3.350343815022304D+00, &
	     6.012459259764103D+00 /)
	  real    ( kind = 8 ) x
	  real    ( kind = 8 ) x1
	  real    ( kind = 8 ) x2
	  real    ( kind = 8 ), parameter :: xlge = 5.10D+05
	  real    ( kind = 8 ), parameter :: xlgst = 1.0D+30
	  real    ( kind = 8 ) xvalue
	  real    ( kind = 8 ) y

	  x = xvalue
	  alngam = 0.0D+00
	!
	!  Check the input.
	!
	  if ( xlgst <= x ) then
	    ifault = 2
	    return
	  end if

	  if ( x <= 0.0D+00 ) then
	    ifault = 1
	    return
	  end if

	  ifault = 0
	!
	!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
	!
	  if ( x < 1.5D+00 ) then

	    if ( x < 0.5D+00 ) then

	      alngam = - log ( x )
	      y = x + 1.0D+00
	!
	!  Test whether X < machine epsilon.
	!
	      if ( y == 1.0D+00 ) then
	        return
	      end if

	    else

	      alngam = 0.0D+00
	      y = x
	      x = ( x - 0.5D+00 ) - 0.5D+00

	    end if

	    alngam = alngam + x * (((( &
	        r1(5)   * y &
	      + r1(4) ) * y &
	      + r1(3) ) * y &
	      + r1(2) ) * y &
	      + r1(1) ) / (((( &
	                  y &
	      + r1(9) ) * y &
	      + r1(8) ) * y &
	      + r1(7) ) * y &
	      + r1(6) )

	    return

	  end if
	!
	!  Calculation for 1.5 <= X < 4.0.
	!
	  if ( x < 4.0D+00 ) then

	    y = ( x - 1.0D+00 ) - 1.0D+00

	    alngam = y * (((( &
	        r2(5)   * x &
	      + r2(4) ) * x &
	      + r2(3) ) * x &
	      + r2(2) ) * x &
	      + r2(1) ) / (((( &
	                  x &
	      + r2(9) ) * x &
	      + r2(8) ) * x &
	      + r2(7) ) * x &
	      + r2(6) )
	!
	!  Calculation for 4.0 <= X < 12.0.
	!
	  else if ( x < 12.0D+00 ) then

	    alngam = (((( &
	        r3(5)   * x &
	      + r3(4) ) * x &
	      + r3(3) ) * x &
	      + r3(2) ) * x &
	      + r3(1) ) / (((( &
	                  x &
	      + r3(9) ) * x &
	      + r3(8) ) * x &
	      + r3(7) ) * x &
	      + r3(6) )
	!
	!  Calculation for 12.0 <= X.
	!
	  else

	    y = log ( x )
	    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

	    if ( x <= xlge ) then

	      x1 = 1.0D+00 / x
	      x2 = x1 * x1

	      alngam = alngam + x1 * ( ( &
	             r4(3)   * &
	        x2 + r4(2) ) * &
	        x2 + r4(1) ) / ( ( &
	        x2 + r4(5) ) * &
	        x2 + r4(4) )

	    end if

	  end if

	  return
	end function alngam
!-----------------------------------------------------------------------
	function alnorm ( x, upper )

	!*****************************************************************************80
	!
	!! ALNORM computes the cumulative density of the standard normal distribution.
	!
	!  Modified:
	!
	!    13 January 2008
	!
	!  Author:
	!
	!    David Hill
	!    FORTRAN90 version by John Burkardt
	!
	!  Reference:
	!
	!    David Hill,
	!    Algorithm AS 66:
	!    The Normal Integral,
	!    Applied Statistics,
	!    Volume 22, Number 3, 1973, pages 424-427.
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
	!    over which the integration takes place.
	!
	!    Input, logical UPPER, determines whether the upper or lower
	!    interval is to be integrated:
	!    .TRUE.  => integrate from X to + Infinity;
	!    .FALSE. => integrate from - Infinity to X.
	!
	!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
	!    distribution over the desired interval.
	!
	  implicit none

	  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
	  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
	  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
	  real ( kind = 8 ) alnorm
	  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
	  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
	  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
	  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
	  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
	  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
	  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
	  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
	  real ( kind = 8 ), parameter :: con = 1.28D+00
	  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
	  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
	  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
	  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
	  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
	  real ( kind = 8 ), parameter :: ltone = 7.0D+00
	  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
	  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
	  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
	  logical up
	  logical upper
	  real ( kind = 8 ), parameter :: utzero = 18.66D+00
	  real ( kind = 8 ) x
	  real ( kind = 8 ) y
	  real ( kind = 8 ) z

	  up = upper
	  z = x

	  if ( z < 0.0D+00 ) then
	    up = .not. up
	    z = - z
	  end if

	  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

	    if ( up ) then
	      alnorm = 0.0D+00
	    else
	      alnorm = 1.0D+00
	    end if

	    return

	  end if

	  y = 0.5D+00 * z * z

	  if ( z <= con ) then

	    alnorm = 0.5D+00 - z * ( p - q * y &
	      / ( y + a1 + b1 &
	      / ( y + a2 + b2 & 
	      / ( y + a3 ))))

	  else

	    alnorm = r * exp ( - y ) &
	      / ( z + c1 + d1 &
	      / ( z + c2 + d2 &
	      / ( z + c3 + d3 &
	      / ( z + c4 + d4 &
	      / ( z + c5 + d5 &
	      / ( z + c6 ))))))

	  end if

	  if ( .not. up ) then
	    alnorm = 1.0D+00 - alnorm
	  end if

	  return
	end function alnorm
!-----------------------------------------------------------------------
	subroutine gamma_inc_values ( n_data, a, x, fx )

	!*****************************************************************************80
	!
	!! GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
	!
	!  Discussion:
	!
	!    The (normalized) incomplete Gamma function P(A,X) is defined as:
	!
	!      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
	!
	!    With this definition, for all A and X,
	!
	!      0 <= PN(A,X) <= 1
	!
	!    and
	!
	!      PN(A,INFINITY) = 1.0
	!
	!    In Mathematica, the function can be evaluated by:
	!
	!      1 - GammaRegularized[A,X]
	!
	!  Modified:
	!
	!    20 November 2004
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Reference:
	!
	!    Milton Abramowitz, Irene Stegun,
	!    Handbook of Mathematical Functions,
	!    National Bureau of Standards, 1964,
	!    ISBN: 0-486-61272-4,
	!    LC: QA47.A34.
	!
	!    Stephen Wolfram,
	!    The Mathematica Book,
	!    Fourth Edition,
	!    Cambridge University Press, 1999,
	!    ISBN: 0-521-64314-7,
	!    LC: QA76.95.W65.
	!
	!  Parameters:
	!
	!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
	!    before the first call.  On each call, the routine increments N_DATA by 1,
	!    and returns the corresponding data; when there is no more data, the
	!    output value of N_DATA will be 0 again.
	!
	!    Output, real ( kind = 8 ) A, the parameter of the function.
	!
	!    Output, real ( kind = 8 ) X, the argument of the function.
	!
	!    Output, real ( kind = 8 ) FX, the value of the function.
	!
	  implicit none

	  integer ( kind = 4 ), parameter :: n_max = 20

	  real ( kind = 8 ) a
	  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
	    0.10D+00, &
	    0.10D+00, &
	    0.10D+00, &
	    0.50D+00, &
	    0.50D+00, &
	    0.50D+00, &
	    0.10D+01, &
	    0.10D+01, &
	    0.10D+01, &
	    0.11D+01, &
	    0.11D+01, &
	    0.11D+01, &
	    0.20D+01, &
	    0.20D+01, &
	    0.20D+01, &
	    0.60D+01, &
	    0.60D+01, &
	    0.11D+02, &
	    0.26D+02, &
	    0.41D+02 /)
	  real ( kind = 8 ) fx
	  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
	    0.7382350532339351D+00, &
	    0.9083579897300343D+00, &
	    0.9886559833621947D+00, &
	    0.3014646416966613D+00, &
	    0.7793286380801532D+00, &
	    0.9918490284064973D+00, &
	    0.9516258196404043D-01, &
	    0.6321205588285577D+00, &
	    0.9932620530009145D+00, &
	    0.7205974576054322D-01, &
	    0.5891809618706485D+00, &
	    0.9915368159845525D+00, &
	    0.1018582711118352D-01, &
	    0.4421745996289254D+00, &
	    0.9927049442755639D+00, &
	    0.4202103819530612D-01, &
	    0.9796589705830716D+00, &
	    0.9226039842296429D+00, &
	    0.4470785799755852D+00, &
	    0.7444549220718699D+00 /)
	  integer ( kind = 4 ) n_data
	  real ( kind = 8 ) x
	  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
	    0.30D-01, &
	    0.30D+00, &
	    0.15D+01, &
	    0.75D-01, &
	    0.75D+00, &
	    0.35D+01, &
	    0.10D+00, &
	    0.10D+01, &
	    0.50D+01, &
	    0.10D+00, & 
	    0.10D+01, &
	    0.50D+01, &
	    0.15D+00, &
	    0.15D+01, &
	    0.70D+01, &
	    0.25D+01, &
	    0.12D+02, &
	    0.16D+02, &
	    0.25D+02, &
	    0.45D+02 /)

	  if ( n_data < 0 ) then
	    n_data = 0
	  end if

	  n_data = n_data + 1

	  if ( n_max < n_data ) then
	    n_data = 0
	    a = 0.0D+00
	    x = 0.0D+00
	    fx = 0.0D+00
	  else
	    a = a_vec(n_data)
	    x = x_vec(n_data)
	    fx = fx_vec(n_data)
	  end if

	  return
	end subroutine gamma_inc_values
!-----------------------------------------------------------------------
	function gammad ( x, p, ifault )

	!*****************************************************************************80
	!
	!! GAMMAD computes the Incomplete Gamma Integral
	!
	!  Auxiliary functions:
	!
	!    ALOGAM = logarithm of the gamma function, 
	!    ALNORM = algorithm AS66
	!
	!  Modified:
	!
	!    20 January 2008
	!
	!  Author:
	!
	!    B Shea
	!    FORTRAN90 version by John Burkardt
	!
	!  Reference:
	!
	!    B Shea,
	!    Algorithm AS 239:
	!    Chi-squared and Incomplete Gamma Integral,
	!    Applied Statistics,
	!    Volume 37, Number 3, 1988, pages 466-473.
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete 
	!    gamma ratio.  0 <= X, and 0 < P.
	!
	!    Output, integer ( kind = 4 ) IFAULT, error flag.
	!    0, no error.
	!    1, X < 0 or P <= 0.
	!
	!    Output, real ( kind = 8 ) GAMMAD, the value of the incomplete 
	!    Gamma integral.
	!
	  implicit none

	  real    ( kind = 8 ) a
	  real    ( kind = 8 ) an
	  real    ( kind = 8 ) arg
	  real    ( kind = 8 ) b
	  real    ( kind = 8 ) c
	  real    ( kind = 8 ), parameter :: elimit = - 88.0D+00
	  real    ( kind = 8 ) gammad
	  integer ( kind = 4 ) ifault
	  real    ( kind = 8 ), parameter :: oflo = 1.0D+37
	  real    ( kind = 8 ) p
	  real    ( kind = 8 ), parameter :: plimit = 1000.0D+00
	  real    ( kind = 8 ) pn1
	  real    ( kind = 8 ) pn2
	  real    ( kind = 8 ) pn3
	  real    ( kind = 8 ) pn4
	  real    ( kind = 8 ) pn5
	  real    ( kind = 8 ) pn6
	  real    ( kind = 8 ) rn
	  real    ( kind = 8 ), parameter :: tol = 1.0D-14
	  logical upper
	  real    ( kind = 8 ) x
	  real    ( kind = 8 ), parameter :: xbig = 1.0D+08

	  gammad = 0.0D+00
	!
	!  Check the input.
	!
	  if ( x < 0.0D+00 ) then
	    ifault = 1
	    return
	  end if

	  if ( p <= 0.0D+00 ) then
	    ifault = 1
	    return
	  end if

	  ifault = 0

	  if ( x == 0.0D+00 ) then
	    gammad = 0.0D+00
	    return
	  end if
	!
	!  If P is large, use a normal approximation.
	!
	  if ( plimit < p ) then

	    pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) &
	    + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )

	    upper = .false.
	    gammad = alnorm ( pn1, upper )
	    return

	  end if
	!
	!  If X is large set GAMMAD = 1.
	!
	  if ( xbig < x ) then
	    gammad = 1.0D+00
	    return
	  end if
	!
	!  Use Pearson's series expansion.
	!  (Note that P is not large enough to force overflow in ALOGAM).
	!  No need to test IFAULT on exit since P > 0.
	!
	  if ( x <= 1.0D+00 .or. x < p ) then

	    arg = p * log ( x ) - x - alngam ( p + 1.0D+00, ifault )
	    c = 1.0D+00
	    gammad = 1.0D+00
	    a = p

	    do

	      a = a + 1.0D+00
	      c = c * x / a
	      gammad = gammad + c

	      if ( c <= tol ) then
	        exit
	      end if

	    end do

	    arg = arg + log ( gammad )

	    if ( elimit <= arg ) then
	      gammad = exp ( arg )
	    else
	      gammad = 0.0D+00
	    end if
	!
	!  Use a continued fraction expansion.
	!
	  else 

	    arg = p * log ( x ) - x - alngam ( p, ifault )
	    a = 1.0D+00 - p
	    b = a + x + 1.0D+00
	    c = 0.0D+00
	    pn1 = 1.0D+00
	    pn2 = x
	    pn3 = x + 1.0D+00
	    pn4 = x * b
	    gammad = pn3 / pn4

	    do

	      a = a + 1.0D+00
	      b = b + 2.0D+00
	      c = c + 1.0D+00
	      an = a * c
	      pn5 = b * pn3 - an * pn1
	      pn6 = b * pn4 - an * pn2

	      if ( pn6 /= 0.0D+00 ) then

	        rn = pn5 / pn6

	        if ( abs ( gammad - rn ) <= min ( tol, tol * rn ) ) then
	          exit
	        end if

	        gammad = rn

	      end if

	      pn1 = pn3
	      pn2 = pn4
	      pn3 = pn5
	      pn4 = pn6
	!
	!  Re-scale terms in continued fraction if terms are large.
	!
	      if ( oflo <= abs ( pn5 ) ) then
	        pn1 = pn1 / oflo
	        pn2 = pn2 / oflo
	        pn3 = pn3 / oflo
	        pn4 = pn4 / oflo
	      end if

	    end do

	    arg = arg + log ( gammad )

	    if ( elimit <= arg ) then
	      gammad = 1.0D+00 - exp ( arg )
	    else
	      gammad = 1.0D+00
	    end if

	  end if

	  return
	end function gammad
!-----------------------------------------------------------------------
	subroutine timestamp ( )

	!*****************************************************************************80
	!
	!! TIMESTAMP prints the current YMDHMS date as a time stamp.
	!
	!  Example:
	!
	!    31 May 2001   9:45:54.872 AM
	!
	!  Modified:
	!
	!    06 August 2005
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    None
	!
	  implicit none

	  character ( len = 8 ) ampm
	  integer d
	  integer h
	  integer m
	  integer mm
	  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
	    'January  ', 'February ', 'March    ', 'April    ', &
	    'May      ', 'June     ', 'July     ', 'August   ', &
	    'September', 'October  ', 'November ', 'December ' /)
	  integer n
	  integer s
	  integer values(8)
	  integer y

	  call date_and_time ( values = values )

	  y = values(1)
	  m = values(2)
	  d = values(3)
	  h = values(5)
	  n = values(6)
	  s = values(7)
	  mm = values(8)

	  if ( h < 12 ) then
	    ampm = 'AM'
	  else if ( h == 12 ) then
	    if ( n == 0 .and. s == 0 ) then
	      ampm = 'Noon'
	    else
	      ampm = 'PM'
	    end if
	  else
	    h = h - 12
	    if ( h < 12 ) then
	      ampm = 'PM'
	    else if ( h == 12 ) then
	      if ( n == 0 .and. s == 0 ) then
	        ampm = 'Midnight'
	      else
	        ampm = 'AM'
	      end if
	    end if
	  end if

	  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
	    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

	  return
	end subroutine timestamp
!-----------------------------------------------------------------------
	function up_inc_gamma( n, x, ifault )
	!*****************************************************************************80
	!
	!! UP_INC_GAMMA computes the Upper Incomplete Gamma Integral
	!
	!  Discussion:
	!
	!    The upper incomplete Gamma function G(X,P) is defined as:
	!
	!      G(X,P) = Integral ( X <= T <= INFINITY ) T**(N-1) * exp(-T) dT.
	!
	!    and can be computed simply as:
	!
	!      G(X,P) = DGAMMA(P) - DW_INC_GAMMA(X,P)
	!
	!  Auxiliary functions:
	!
	!    ALOGAM = logarithm of the gamma function, 
	!    ALNORM = algorithm AS66
	!    DGAMMA = intrinsic gfortran function
	!    GAMMAD = algorithm AS239
	!    DW_INC_GAMMA = defined down here
	!
	!  Modified:
	!
	!    9 November 2010
	!
	!  Author:
	!
	!    Francesco Calcavecchia
	!
	!  Reference:
	!
	!    http://en.wikipedia.org/wiki/Incomplete_gamma_function
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) X, P, the parameters of the upper incomplete 
	!    gamma function.  0 <= X, and 0 < P.
	!
	!    Output, integer ( kind = 4 ) IFAULT, error flag.
	!    0, no error.
	!    1, X < 0 or P <= 0.
	!
	!    Output, real ( kind = 8 ) GAMMAD, the value of the incomplete 
	!    Gamma integral.
	!
	  IMPLICIT NONE

	  real    ( kind = 8 ) up_inc_gamma
	  integer ( kind = 4 ) ifault
	  real    ( kind = 8 ) n
	  real    ( kind = 8 ) x

	  up_inc_gamma = DGAMMA(n) - DW_INC_GAMMA(n,x,ifault)

	end function up_inc_gamma
!-----------------------------------------------------------------------
	function dw_inc_gamma( n, x, ifault )
	!*****************************************************************************80
	!
	!! DW_INC_GAMMA computes the Lower Incomplete Gamma Integral
	!
	!  Discussion:
	!
	!    The lower incomplete Gamma function G(X,P) is defined as:
	!
	!      G(X,P) = Integral ( 0 <= T <= X ) T**(N-1) * exp(-T) dT.
	!
	!    and can be computed simply as:
	!
	!      G(X,P) = GAMMAD(X,P) * DGAMMA(P)
	!
	!  Auxiliary functions:
	!
	!    ALOGAM = logarithm of the gamma function, 
	!    ALNORM = algorithm AS66
	!    DGAMMA = intrinsic gfortran function
	!    GAMMAD = algorithm AS239
	!
	!  Modified:
	!
	!    9 November 2010
	!
	!  Author:
	!
	!    Francesco Calcavecchia
	!
	!  Reference:
	!
	!    http://en.wikipedia.org/wiki/Incomplete_gamma_function
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) X, P, the parameters of the upper incomplete 
	!    gamma function.  0 <= X, and 0 < P.
	!
	!    Output, integer ( kind = 4 ) IFAULT, error flag.
	!    0, no error.
	!    1, X < 0 or P <= 0.
	!
	!    Output, real ( kind = 8 ) GAMMAD, the value of the incomplete 
	!    Gamma integral.
	!
	  IMPLICIT NONE

	  real    ( kind = 8 ) dw_inc_gamma
	  integer ( kind = 4 ) ifault
	  real    ( kind = 8 ) n
	  real    ( kind = 8 ) x

	  IF (n>0) THEN
	    dw_inc_gamma = GAMMAD(x,n,ifault)*DGAMMA(n)
	  ELSE IF (n==0) THEN
	    dw_inc_gamma = 0.d0
	    ifault=0
	  END IF

	end function dw_inc_gamma
	
END MODULE generic_tools


