!!!!SAMPLE OF A dati_dft.d file
!&dati_dft
!pw_exe="/home/kenzo/espresso-5.0.2/bin/pw.x"
!pseudopot_dir="/home/kenzo/espresso-5.0.2/pseudo"
!iotk_xml_exe="/home/kenzo/espresso-5.0.2/bin/iotk"
!dft_work_dir="DFT/"
!scf_in_file="scf.in"
!dft_parallel=.FALSE.
!dft_mpi_nprocs=4
!/

MODULE dft
	IMPLICIT NONE
	LOGICAL, PROTECTED, SAVE :: dft_parallel, iniz_dft=.FALSE., dft_lock=.FALSE.
	CHARACTER(LEN=250), PROTECTED, SAVE :: original_work_dir, dft_work_dir, scf_in_file, pw_exe, pseudopot_dir, iotk_xml_exe
	INTEGER, PROTECTED, SAVE :: dft_mpi_nprocs
	
	CONTAINS
		
	SUBROUTINE inizializza_dft()
		IMPLICIT NONE
		LOGICAL :: flag_work_dir
		NAMELIST /dati_dft/ pw_exe, pseudopot_dir, iotk_xml_exe, dft_work_dir, scf_in_file, dft_parallel, dft_mpi_nprocs
		
		IF (iniz_dft) STOP 'dati_dft giá inizializzati &
		  [ module_dft.f90 > inizializza_dft ]'
		
		CALL GETCWD(original_work_dir)
		dft_lock=.FALSE.
		
		OPEN (2, FILE='dati_dft.d',STATUS='OLD')
		READ (2, NML=dati_dft)
		CLOSE (2)

		INQUIRE(FILE=dft_work_dir, EXIST=flag_work_dir)
		IF ( .NOT. flag_work_dir ) THEN
			CALL SYSTEM("mkdir "//dft_work_dir)
		END IF
		
		iniz_dft=.TRUE.
		
	END SUBROUTINE inizializza_dft
		
	SUBROUTINE run_quantum_espresso()
		IMPLICIT NONE
		CHARACTER(LEN=4) :: istring
		CHARACTER(LEN=250) :: original_work_dir
		LOGICAL :: flag_sh_file
		INTEGER :: ios
	
		IF (.NOT. iniz_dft) STOP 'Non hai inizializzato dati_dft &
		  [ module_dft.f90 > run_quantum_espresso ]'
		
		dft_lock=.TRUE.
		
		CALL GETCWD(original_work_dir)
		CALL CHDIR(TRIM(original_work_dir)//"/"//TRIM(dft_work_dir))
		IF ( dft_parallel ) THEN
			WRITE (istring, '(I4.4)'), dft_mpi_nprocs
			!PRINT * , "mpirun -np "//istring//" "//TRIM(pw_exe)//" < "//TRIM(scf_in_file)//" > scf.out"
			CALL SYSTEM("mpirun -np "//istring//" "//TRIM(pw_exe)//" < "//TRIM(scf_in_file)//" > scf.out")
		ELSE
			OPEN(UNIT=99, FILE="starter.sh", STATUS='UNKNOWN', IOSTAT=ios)
			IF ( ios /= 0 ) STOP "Error opening file starter.sh"
			WRITE(UNIT=99, FMT=*) TRIM(pw_exe)//" < "//TRIM(scf_in_file)//" > scf.out"
			CLOSE (99)
			CALL SYSTEM("sh starter.sh")
			CALL SLEEP(5)
			!CALL SYSTEM(TRIM(pw_exe)//" < "//TRIM(scf_in_file)//" > scf.out")
		END IF
		INQUIRE(FILE="out_xml.sh",EXIST=flag_sh_file)
		IF ( .NOT. flag_sh_file ) THEN
			CALL genera_file_out_xml_sh()
		END IF
		CALL SYSTEM("sh out_xml.sh")
    
		CALL CHDIR(TRIM(original_work_dir))
		
		dft_lock=.FALSE.
    
	END SUBROUTINE run_quantum_espresso
    
    
	SUBROUTINE genera_file_out_xml_sh()
		IMPLICIT NONE
		INTEGER :: ios
		
		IF (.NOT. iniz_dft) STOP 'Non hai inizializzato dati_dft &
		  [ module_dft.f90 > genera_file_out_xml_sh ]'
    
		OPEN(UNIT=99, FILE='out_xml.sh', STATUS='UNKNOWN', IOSTAT=ios)
		IF ( ios /= 0 ) STOP "Error opening file out_xml.sh"
    
		WRITE(UNIT=99, FMT=*), "cd OUT.save"
		WRITE(UNIT=99, FMT=*), "for CARTELLA in K*/"
		WRITE(UNIT=99, FMT=*), "do"
		WRITE(UNIT=99, FMT=*), "	cd $CARTELLA"
		WRITE(UNIT=99, FMT=*), "	rm -f *.xml"
		WRITE(UNIT=99, FMT=*), "	"//iotk_xml_exe//" convert gkvectors.dat gkvectors.xml >/dev/null  2>&1"
		WRITE(UNIT=99, FMT=*), "	"//iotk_xml_exe//" convert evc.dat evc.xml >/dev/null  2>&1"
		WRITE(UNIT=99, FMT=*), "	rm -f *.dat"
		WRITE(UNIT=99, FMT=*), "	cd .."
		WRITE(UNIT=99, FMT=*), "done"
		WRITE(UNIT=99, FMT=*), "cd .."
    
		CLOSE (99)
    
	END SUBROUTINE genera_file_out_xml_sh
    
    
	SUBROUTINE genera_scf_in(L,N,x)
		IMPLICIT NONE
		REAL (KIND=8), INTENT(IN) :: L(1:3), x(1:3,N)
		INTEGER, INTENT(IN) :: N
		INTEGER :: ios, i1
		
		IF (.NOT. iniz_dft) STOP 'Non hai inizializzato dati_dft &
		  [ module_dft.f90 > genera_scf_in ]'
    
		OPEN(UNIT=99, FILE=TRIM(dft_work_dir)//TRIM(scf_in_file), STATUS='UNKNOWN', IOSTAT=ios)
		IF ( ios /= 0 ) STOP "Error creating file scf_in_file.in"
    
		WRITE(UNIT=99, FMT=*), " &CONTROL                                                                "
		WRITE(UNIT=99, FMT=*), "   calculation =   'scf'   ,                                           "
		WRITE(UNIT=99, FMT=*), "   restart_mode = 'from_scratch' ,                                     "
		WRITE(UNIT=99, FMT=*), "   pseudo_dir = '"//TRIM(pseudopot_dir)//"' ,    "
		WRITE(UNIT=99, FMT=*), "   prefix = 'OUT',                                                     "
		WRITE(UNIT=99, FMT=*), "   wf_collect = .true. ,                                               "
		WRITE(UNIT=99, FMT=*), " /                                                                      "
		WRITE(UNIT=99, FMT=*), " &SYSTEM                                                                "
		WRITE(UNIT=99, FMT=*), "   ibrav = 0,                                                          "
		WRITE(UNIT=99, FMT=*), "   celldm(1) = ", L(1), ","
		WRITE(UNIT=99, FMT=*), "   nat  = ", N," ,  "
		WRITE(UNIT=99, FMT=*), "   ntyp  =  1 ,"
		WRITE(UNIT=99, FMT=*), "   ecutwfc  =   8., "      !"   ecutwfc  =   8., "
		WRITE(UNIT=99, FMT=*), "   occupations  =  'smearing' ,"
		WRITE(UNIT=99, FMT=*), "   smearing  =  'gaussian'"
		WRITE(UNIT=99, FMT=*), "   degauss  =  0.005,"
		WRITE(UNIT=99, FMT=*), " /"
		WRITE(UNIT=99, FMT=*), " &ELECTRONS"
		WRITE(UNIT=99, FMT=*), "   conv_thr  =  1.0d-6  ,"
		WRITE(UNIT=99, FMT=*), "   mixing_beta = 0.5 ,"
		WRITE(UNIT=99, FMT=*), " /"
		WRITE(UNIT=99, FMT=*), "CELL_PARAMETERS (cubic)"
		WRITE(UNIT=99, FMT=*), "   ", L(1)/L(1), 0.d0, 0.d0
		WRITE(UNIT=99, FMT=*), "   ", 0.d0, L(2)/L(1), 0.d0
		WRITE(UNIT=99, FMT=*), "   ", 0.d0, 0.d0, L(3)/L(1)
		WRITE(UNIT=99, FMT=*), "ATOMIC_SPECIES"
		WRITE(UNIT=99, FMT=*), "  H  1.00794 H.coulomb-ae.UPF"
		WRITE(UNIT=99, FMT=*), "ATOMIC_POSITIONS (crystal)     "
		DO i1 = 1, N, 1
			WRITE(UNIT=99, FMT=*), "H   ", x(1:3,i1)
		END DO
		WRITE(UNIT=99, FMT=*), "K_POINTS automatic"
		WRITE(UNIT=99, FMT=*), "5  5  5  0  0  0"     !"5  5  5  0  0  0"
    
		CLOSE(99)
    
	END SUBROUTINE genera_scf_in
	
	
	SUBROUTINE chiudi_dft()
		IMPLICIT NONE
		
		IF (.NOT. iniz_dft) STOP 'Non hai inizializzato dati_dft &
		  [ module_dft.f90 > chiudi_dft ]'
		
		iniz_dft=.FALSE.
		
	END SUBROUTINE chiudi_dft
		
END MODULE dft

