 ! An adaption of driver code from the i-PI MD software
 ! which serves as HswfQMC's interface to i-PI
 !
 ! Copyright (C) 2013, Joshua More and Michele Ceriotti
 ! Copyright (C) 2017, Jan Kessler
 !
 ! Permission is hereby granted, free of charge, to any person obtaining
 ! a copy of this software and associated documentation files (the
 ! "Software"), to deal in the Software without restriction, including
 ! without limitation the rights to use, copy, modify, merge, publish,
 ! distribute, sublicense, and/or sell copies of the Software, and to
 ! permit persons to whom the Software is furnished to do so, subject to
 ! the following conditions:
 !
 ! The above copyright notice and this permission notice shall be included
 ! in all copies or substantial portions of the Software.
 !
 ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 ! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 ! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 ! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 ! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 ! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 ! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 !
 ! -------------------------------------------------------------------------


 PROGRAM IPI_CLIENT

   USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer
   IMPLICIT NONE

   ! SOCKET VARIABLES
   INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
   INTEGER socket, inet, port        ! socket ID & address of the server
   CHARACTER(LEN=1024) :: host

   ! COMMAND LINE PARSING
   CHARACTER(LEN=1024) :: cmdbuffer
   CHARACTER(LEN=3) :: vstyle
   CHARACTER(LEN=10) :: ccmd
   INTEGER verbose
   INTEGER commas(2), par_count_o, par_count_c      ! stores the index of commas in the parameter string
   INTEGER vpars_o(2), vpars_c(3)         ! array to store the parameters of the potential

   ! SOCKET COMMUNICATION BUFFERS
   CHARACTER(LEN=12) :: header
   LOGICAL :: isinit, hasdata
   INTEGER cbuf, rid, istep
   CHARACTER(LEN=2048) :: initbuffer      ! it's unlikely a string this large will ever be passed...
   DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)

   ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
   INTEGER nat
   DOUBLE PRECISION pot, volume, box(3)
   DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:), forces(:,:)
   DOUBLE PRECISION cell_h(3,3), cell_ih(3,3), virial(3,3), mtxbuf(9)

   ! PARAMETERS CONCERNING VMC (NRANKS,NWOMIN)
   INTEGER nranks, nnodes, mpicom, nwomin
   CHARACTER(LEN=1024) ridstr
   LOGICAL domshift, dowfopt
   INTEGER, ALLOCATABLE :: invindarr(:)
   DOUBLE PRECISION :: flimit

   ! ITERATORS
   INTEGER i

   ! intialize parameter defaults
   inet = 1
   verbose = 0
   host = 'localhost'//achar(0)
   port = 54321
   vstyle = 'nop'
   domshift = .FALSE.
   nwomin = 0
   nranks = 1
   nnodes = 1
   mpicom = 0
   flimit = 0.d0
   dowfopt = .FALSE.

   ! initialize control variables
   isinit = .FALSE.
   hasdata = .FALSE.
   ccmd = ''
   par_count_o = 0
   par_count_c = 0

   nat = -1
   rid = 0
   istep = 0

   ! parse the command line parameters
   DO i = 1, COMMAND_ARGUMENT_COUNT()

      CALL GET_COMMAND_ARGUMENT(i, cmdbuffer)

      IF (cmdbuffer == '-u') THEN ! flag for unix socket
         inet = 0
      ELSEIF (cmdbuffer == '-v') THEN ! flag for verbose standard output
         verbose = 1
      ELSEIF (cmdbuffer == '-vv') THEN ! flag for verbose standard output
         verbose = 2
      ELSEIF (cmdbuffer == '-h') THEN ! read the hostname
         ccmd = 'host'
      ELSEIF (cmdbuffer == '-p') THEN ! reads the port number
         ccmd = 'port'
      ELSEIF (cmdbuffer == '-m') THEN ! reads the mode of electronic parameter optimization
         ccmd = 'wfmode'
      ELSEIF (cmdbuffer == '-o') THEN ! reads the options for mode
         ccmd = 'wfopts'
      ELSEIF (cmdbuffer == '-c') THEN ! reads the MPI/Rank configuration
         ccmd = 'mpiconf'
      ELSEIF (cmdbuffer == '-l') THEN ! reads the force amplitude limit
         ccmd = 'flimit'

      ELSE
         IF (ccmd == 'host') THEN
            host = trim(cmdbuffer)//achar(0)

         ELSEIF (ccmd == 'port') THEN
            READ(cmdbuffer,*) port

         ELSEIF (ccmd == 'wfmode') THEN
            IF (trim(cmdbuffer) == 'alt') THEN
               vstyle = 'alt'
               dowfopt = .TRUE.
            ELSEIF (trim(cmdbuffer) == 'sim') THEN
               vstyle = 'sim'
               dowfopt = .TRUE.
            ELSEIF (trim(cmdbuffer) == 'fix') THEN
               vstyle = 'fix'
               dowfopt = .FALSE.
            ELSEIF (trim(cmdbuffer) == 'ffs') THEN
               vstyle = 'ffs'
               dowfopt = .FALSE.
            ELSEIF (trim(cmdbuffer) == 'nop') THEN
               vstyle = 'nop'
               dowfopt = .FALSE.
            ELSE
               WRITE(*,*) ' Unrecognized wavefunction optimization mode ', trim(cmdbuffer)
               WRITE(*,*) ' Use -m [alt|sim|fix|ffs|nop] '
               STOP 'ENDED'
            ENDIF

         ELSEIF (ccmd == 'wfopts') THEN
            par_count_o = 1
            commas(1) = 0
            DO WHILE (index(cmdbuffer(commas(par_count_o)+1:), ',') > 0)
               commas(par_count_o + 1) = index(cmdbuffer(commas(par_count_o)+1:), ',') + commas(par_count_o)
               READ(cmdbuffer(commas(par_count_o)+1:commas(par_count_o + 1)-1),*) vpars_o(par_count_o)
               par_count_o = par_count_o + 1
            ENDDO
            READ(cmdbuffer(commas(par_count_o)+1:),*) vpars_o(par_count_o)

         ELSEIF (ccmd == 'mpiconf') THEN
            par_count_c = 1
            commas(1) = 0
            DO WHILE (index(cmdbuffer(commas(par_count_c)+1:), ',') > 0)
               commas(par_count_c + 1) = index(cmdbuffer(commas(par_count_c)+1:), ',') + commas(par_count_c)
               READ(cmdbuffer(commas(par_count_c)+1:commas(par_count_c + 1)-1),*) vpars_c(par_count_c)
               par_count_c = par_count_c + 1
            ENDDO
            READ(cmdbuffer(commas(par_count_c)+1:),*) vpars_c(par_count_c)

         ELSEIF (ccmd == 'flimit') THEN
            READ(cmdbuffer(:),*) flimit ! force amplitude limit (dirty fix)

         ELSE
            WRITE(*,*) ' Unrecognized command line argument: ', cmdbuffer
            CALL helpmessage
            STOP 'ENDED'
         ENDIF
         ccmd = ''
      ENDIF
   ENDDO

   IF (vstyle == 'nop' .or. vstyle == 'ffs') THEN
      CONTINUE
   ELSEIF (vstyle == 'fix' .or. vstyle == 'sim') THEN
      IF (par_count_o /= 1) THEN
         WRITE(*,*) 'Error: Wrong number of -o parameters provided. Expected: -o ISHIFT'
         STOP 'ENDED'
      ENDIF
      IF (vpars_o(1) > 0) THEN
         domshift = .TRUE.
      ELSE
         domshift = .FALSE.
      END IF
   ELSEIF (vstyle == 'alt') THEN
      IF (par_count_o /= 2) THEN
         WRITE(*,*) 'Error: Wrong number of -o parameters provided. Expected: -o ISHIFT,NWOMIN'
         STOP 'ENDED'
      ENDIF
      IF (vpars_o(1) > 0) THEN
         domshift = .TRUE.
      ELSE
         domshift = .FALSE.
      END IF
      nwomin = vpars_o(2)
   ENDIF

   IF (par_count_c == 1) THEN
      nranks = vpars_c(1)
   ELSEIF (par_count_c == 2) THEN
      nranks = vpars_c(1)
      mpicom = vpars_c(2)
   ELSEIF (par_count_c == 3) THEN
      nranks = vpars_c(1)
      mpicom = vpars_c(2)
      nnodes = vpars_c(3)
   ELSE
      WRITE(*,*) 'Error: Wrong number of -c parameters provided. Expected: -c NRANKS , -c NRANKS,MPICOM or -c NRANKS,MPICOM,NNODES'
      STOP 'ENDED'
   END IF

   IF (dowfopt) THEN
      isinit = .FALSE. ! WF optimization enabled -> we need to know our bead id
   ELSE
      isinit = .TRUE.
   END IF

   !Backup the initial random seed file
   CALL execute_command_line('cp randomseed.d randomseed.bak', WAIT = .true.)

   IF (verbose > 0) THEN
      WRITE(*,*) ' DRIVER - Connecting to host ', trim(host)
      IF (inet > 0) THEN
         WRITE(*,*) ' on port ', port, ' using an internet socket.'
      ELSE
         WRITE(*,*) ' using an UNIX socket.'
      ENDIF
   ENDIF

   ! Calls the interface to the POSIX sockets library to open a communication channel
   CALL open_socket(socket, inet, port, host)

   DO WHILE (.true.) ! Loops forever (or until the wrapper ends!)

      ! Reads from the socket one message header
      CALL readbuffer(socket, header, MSGLEN)
      IF (verbose > 0) WRITE(*,*) ' Message from server: ', trim(header)

      IF (trim(header) == 'STATUS') THEN
         ! The wrapper is inquiring on what we are doing
         IF (.not. isinit) THEN
            CALL writebuffer(socket,'NEEDINIT    ',MSGLEN)  ! Signals that we need initialization data
            IF (verbose > 1) WRITE(*,*) '    !write!=> ', 'NEEDINIT    '
         ELSEIF (hasdata) THEN
            CALL writebuffer(socket,'HAVEDATA    ',MSGLEN)  ! Signals that we are done computing and can return forces
            IF (verbose > 1) WRITE(*,*) '    !write!=> ', 'HAVEDATA    '
         ELSE
            CALL writebuffer(socket,'READY       ',MSGLEN)  ! We are idling and eager to compute something
            IF (verbose > 1) WRITE(*,*) '    !write!=> ', 'READY       '
         ENDIF

      ELSEIF (trim(header) == 'INIT') THEN     ! The driver is kindly providing a string for initialization
         CALL readbuffer(socket, rid)
         IF (verbose > 1) WRITE(*,*) '    !read!=> RID: ', rid
         CALL readbuffer(socket, cbuf)
         IF (verbose > 1) WRITE(*,*) '    !read!=> init_lenght: ', cbuf
         CALL readbuffer(socket, initbuffer, cbuf)
         IF (verbose > 1) WRITE(*,*) '    !read!=> init_string: ', cbuf
         IF (verbose > 0) WRITE(*,*) ' Initializing system from wrapper, using ', trim(initbuffer)

         WRITE(ridstr, *) rid
         ridstr = trim(adjustl(ridstr))

         isinit=.TRUE. 

      ELSEIF (trim(header) == 'POSDATA') THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!
         ! Parses the flow of data from the socket
         CALL readbuffer(socket, mtxbuf, 9)  ! Cell matrix
         IF (verbose > 1) WRITE(*,*) '    !read!=> cell: ', mtxbuf
         cell_h = RESHAPE(mtxbuf, (/3,3/))
         CALL readbuffer(socket, mtxbuf, 9)  ! Inverse of the cell matrix (so we don't have to invert it every time here)
         IF (verbose > 1) WRITE(*,*) '    !read!=> cell-1: ', mtxbuf
         cell_ih = RESHAPE(mtxbuf, (/3,3/))

         ! We assume an upper triangular cell-vector matrix
         box(1) = cell_h(1,1)
         box(2) = cell_h(2,2)
         box(3) = cell_h(3,3)
         volume = box(1)*box(2)*box(3)

         CALL readbuffer(socket, cbuf)       ! The number of atoms in the cell
         IF (verbose > 1) WRITE(*,*) '    !read!=> cbuf: ', cbuf
         IF (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation, so only does this once
            nat = cbuf
            IF (verbose > 0) WRITE(*,*) ' Allocating buffer and data arrays, with ', nat, ' atoms'
            ALLOCATE(msgbuffer(3*nat))
            ALLOCATE(atoms(3,nat))
            ALLOCATE(forces(3,nat))
            IF (domshift) ALLOCATE(invindarr(nat))
            msgbuffer = 0.0d0
            atoms = 0.0d0
            forces = 0.0d0
            invindarr = 0.0d0
         ENDIF

         CALL readbuffer(socket, msgbuffer, nat*3)
         IF (verbose > 1) WRITE(*,*) '    !read!=> positions: ', msgbuffer
         DO i = 1, nat
            atoms(:,i) = msgbuffer(3*(i-1)+1:3*i)
         ENDDO

         pot = 0.0d0
         virial = 0.0d0

         istep = istep + 1
         IF (verbose > 0) THEN
            WRITE(*,*) ' Force calculation step ', istep
            WRITE(*,*) ' Calculating forces for bead ', rid
         END IF

         IF (vstyle /= 'nop') THEN

            IF (domshift) CALL molshift(nat, box, atoms, invindarr)

            IF (dowfopt) THEN ! set bead's wavefunction for calculation
               CALL execute_command_line('cp ../SR_wf.dir/'//ridstr//' wf_now.d', WAIT = .true.)
            END IF

            CALL calc_pot_forces(nat, mpicom, nranks, nnodes, vstyle, flimit, 1.d-2, box, atoms, pot, forces)

            IF (dowfopt) THEN ! backup new wavefunction for next step
               CALL execute_command_line('cp ottimizzazione/SR_wf.d ../SR_wf.dir/'//ridstr, WAIT = .true.)
            END IF

            IF (domshift) THEN
               CALL molshift_back(nat, atoms, invindarr)
               CALL molshift_back(nat, forces, invindarr)
            END IF

         END IF

         IF (verbose > 0) WRITE(*,*) ' Calculated energy is ', pot

         hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper

      ELSEIF (trim(header) == 'GETFORCE') THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

         DO i = 1, nat
            msgbuffer(3*(i-1)+1:3*i) = forces(:,i)
         ENDDO

         CALL writebuffer(socket,'FORCEREADY  ',MSGLEN)
         IF (verbose > 1) WRITE(*,*) '    !write!=> ', 'FORCEREADY  '
         CALL writebuffer(socket,pot)  ! Writing the potential
         IF (verbose > 1) WRITE(*,*) '    !write!=> pot: ', pot
         CALL writebuffer(socket,nat)  ! Writing the number of atoms
         IF (verbose > 1) WRITE(*,*) '    !write!=> nat:', nat
         CALL writebuffer(socket,msgbuffer,3*nat) ! Writing the forces
         IF (verbose > 1) WRITE(*,*) '    !write!=> forces:', msgbuffer
         CALL writebuffer(socket,reshape(virial,(/9/)),9)  ! Writing the virial tensor, NOT divided by the volume
         IF (verbose > 1) WRITE(*,*) '    !write!=> strss: ', reshape(virial,(/9/))


         cbuf = 7 ! Size of the 'extras' string
         CALL writebuffer(socket,cbuf) ! This would write out the 'extras' string, but in this case we only use a dummy string.
         IF (verbose > 1) WRITE(*,*) '    !write!=> extra_lenght: ', cbuf
         CALL writebuffer(socket,'nothing',7)
         IF (verbose > 1) WRITE(*,*) '    !write!=> extra: nothing'

         hasdata = .FALSE.
         IF (dowfopt) isinit=.FALSE. ! we want to make sure that the beadid didn't change next time

      ELSE
         WRITE(*,*) ' Unexpected header ', header
         STOP 'ENDED'
      ENDIF
   ENDDO

   IF (nat > 0) THEN
      DEALLOCATE(atoms, forces, msgbuffer)
      IF (domshift) DEALLOCATE(invindarr)
   END IF


 CONTAINS


   SUBROUTINE calc_pot_forces(nat, mpicom, nranks, nnodes, vstyle, flimit, fddx, box, atoms, pot, forces)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nat, mpicom, nranks, nnodes
     CHARACTER(LEN=3), INTENT(IN) :: vstyle
     DOUBLE PRECISION, INTENT(IN) :: flimit, fddx, box(3), atoms(3, nat)
     DOUBLE PRECISION, INTENT(OUT) :: pot, forces(3, nat)

     DO WHILE (.TRUE.)

        !CALL calc_stocrec_all(nat, mpicom, nranks, nnodes, box, atoms, pot, forces)
        CALL calc_findiff(nat, mpicom, nranks, nnodes, fddx, box, atoms, pot, forces)

        IF (notlimit(nat, forces, flimit)) THEN
           IF (vstyle == 'ffs') THEN
              CALL seedshift('randomseed.d', 'randomseed.new', 20, 21, .TRUE.)
           END IF
           EXIT
        ELSE
           WRITE(*,*) 'Force was above limit. Recalculating with new seed...'
           CALL seedshift('randomseed.d', 'randomseed.new', 20, 21, .TRUE.)
        END IF

     ENDDO

   END SUBROUTINE calc_pot_forces


   SUBROUTINE calc_stocrec_all(nat, mpicom, nranks, nnodes, box, atoms, pot, forces)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nat, mpicom, nranks, nnodes
     DOUBLE PRECISION, INTENT(IN) :: box(3), atoms(3, nat)
     DOUBLE PRECISION, INTENT(OUT) :: pot, forces(3, nat)

     CALL write_atoms(nat, box, atoms)
     CALL exec_hswfqmc(mpicom, nranks, nnodes)
     CALL read_forces_ld(nat, forces)
     CALL read_epot_sr(nat, pot)

   END SUBROUTINE calc_stocrec_all


   SUBROUTINE calc_stocrec_nofrc(nat, mpicom, nranks, nnodes, box, atoms, pot)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nat, mpicom, nranks, nnodes
     DOUBLE PRECISION, INTENT(IN) :: box(3), atoms(3, nat)
     DOUBLE PRECISION, INTENT(OUT) :: pot

     CALL write_atoms(nat, box, atoms)
     CALL exec_hswfqmc(mpicom, nranks, nnodes)
     CALL read_epot_sr(nat, pot)

   END SUBROUTINE calc_stocrec_nofrc


   SUBROUTINE calc_findiff(nat, mpicom, nranks, nnodes, fddx, box, atoms, pot, forces)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nat, mpicom, nranks, nnodes
     DOUBLE PRECISION, INTENT(IN) :: fddx, box(3), atoms(3, nat)
     DOUBLE PRECISION, INTENT(OUT) :: pot, forces(3, nat)

     INTEGER i,j
     DOUBLE PRECISION :: dforces(3, nat), datoms(3, nat), dpot(2)

     datoms = atoms
     DO i=1,nat
        DO j=1,3
           datoms(j,i) = atoms(j,i) + fddx
           CALL calc_stocrec_nofrc(nat, mpicom, nranks, nnodes, box, datoms, dpot(1))

           datoms(j,i) = atoms(j,i) - fddx
           CALL calc_stocrec_nofrc(nat, mpicom, nranks, nnodes, box, datoms, dpot(2))

           forces(j,i) = (dpot(2) - dpot(1)) / (2*fddx)
           datoms(j,i) = atoms(j,i)
        ENDDO
     ENDDO

     dforces = forces

     CALL calc_stocrec_all(nat, mpicom, nranks, nnodes, box, atoms, pot, forces)
     write(*,*) dforces
     write(*,*) forces
     forces = dforces
     !CALL calc_simpcal(nat, mpicom, nranks, nnodes, box, atoms, pot)

   END SUBROUTINE calc_findiff


   SUBROUTINE exec_hswfqmc(mpicom, nranks, nnodes)
     IMPLICIT NONE

     INTEGER, INTENT(IN) ::  mpicom, nranks, nnodes

     INTEGER rpn
     CHARACTER(LEN=1024) nrankstr, rpnstr

     rpn = nranks / nnodes ! ranks per node

     WRITE(nrankstr, *) nranks
     nrankstr = trim(adjustl(nrankstr))

     WRITE(rpnstr, *) rpn
     rpnstr = trim(adjustl(rpnstr))

     IF (mpicom == 1) THEN
        CALL execute_command_line('srun -n ' // nrankstr // ' HswfQMC_exe', WAIT = .true.)
     ELSE IF (mpicom == 2) THEN
        CALL execute_command_line('runjob --np ' // nrankstr &
             //' --exe HswfQMC_exe --ranks-per-node ' // rpnstr, WAIT = .true.)
     ELSE
        CALL execute_command_line('mpirun -np ' // nrankstr // ' HswfQMC_exe', WAIT = .true.)
     END IF

     CALL check_state()

   END SUBROUTINE exec_hswfqmc


   SUBROUTINE write_atoms(nat, box, atoms)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nat
     DOUBLE PRECISION, INTENT(IN) :: box(3), atoms(3, nat)

     INTEGER :: i

     OPEN (UNIT=20, FILE='reticolo/rp_now.d', STATUS='REPLACE', ACTION='WRITE')
     WRITE(20,*) box
     DO i = 1, nat
        WRITE(20,*) atoms(:,i)
     ENDDO
     CLOSE(20)

   END SUBROUTINE write_atoms


   SUBROUTINE read_forces_ld(nat, forces)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nat
     DOUBLE PRECISION, INTENT(OUT) :: forces(3, nat)

     INTEGER :: i

     OPEN (UNIT=20, FILE='reticolo/LagrDyn_Frp-0000.d', ACTION='READ')
     DO i = 1, nat
        READ(20, *) forces(:,i)
     ENDDO
     CLOSE(20)
     forces(:,:) = forces(:,:) * 0.5*nat !HswfQMC output values are per atom, convert to Hartree

   END SUBROUTINE read_forces_ld


   SUBROUTINE read_epot_sr(nat, pot)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nat
     DOUBLE PRECISION, INTENT(OUT) :: pot

     DOUBLE PRECISION enhelp(3)
     INTEGER ios

     ios = 0
     OPEN (UNIT=20, FILE='ottimizzazione/SR_energies.dat', ACTION='READ')
     DO WHILE(ios.eq.0)
        READ(20,*,iostat=ios) enhelp
     ENDDO
     CLOSE(UNIT=20)
     pot = enhelp(2) * 0.5*nat

   END SUBROUTINE read_epot_sr


   LOGICAL FUNCTION notlimit(nat, forces, flimit)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nat
     DOUBLE PRECISION, INTENT(IN) :: forces(3, nat), flimit

     INTEGER :: i

     notlimit = .TRUE.
     IF (flimit > 0) THEN
        DO i = 1, nat
           IF (NORM2(forces(:,i)) > flimit) THEN
              notlimit = .FALSE.
              RETURN
           END IF
        ENDDO
     END IF

   END FUNCTION notlimit


   SUBROUTINE check_state()
     IMPLICIT NONE

     CHARACTER(LEN=4) state

     OPEN(UNIT=20, FILE='state.out', ACTION='READ')
     READ(20,*) state
     IF (state /= 'done') STOP '[ipi_client] An error ocurred during execution of HswfQMC_exe.'
     CLOSE(UNIT=20)

   END SUBROUTINE check_state


   SUBROUTINE helpmessage
     ! Help banner
     WRITE(*,*) ' SYNTAX: ipi_driver.x [-u] [-v] -h hostname -p port -m [alt|sim|fix|ffs|nop] ' &
          //'[-o comma,separated,options] [-c comma,separated,options] '
     WRITE(*,*) ''
     WRITE(*,*) ' The -u flag enables the use of unix sockets instead of internet sockets. '
     WRITE(*,*) ''
     WRITE(*,*) ' The -v flag enables verbose mode. '
     WRITE(*,*) ''
     WRITE(*,*) ' You may want to provide options via -o according to the wavefunction optimization mode:'
     WRITE(*,*) ' Fixed wavefunction mode (-m fix):                     -o ISHIFT '
     WRITE(*,*) ' Simultaneous wavefunction optimization mode (-m sim): -o ISHIFT '
     WRITE(*,*) ' Alternating wavefunction optimization mode (-m alt):  -o ISHIFT,NWOMIN (not implemented yet!)'
     WRITE(*,*) ' In case of the no operation mode (-m nop) and fixed force sampling (-m ffs) there are no options to be set.'
     WRITE(*,*) ''
     WRITE(*,*) ' Option Documentation:'
     WRITE(*,*) ' ISHIFT = Molecular Partner Shifting Flag. 0 -> disabled, 1 -> enabled'
     WRITE(*,*) ' NWOMIN = Number of wavefunction optimizations without new minimum before termination.'
     WRITE(*,*) ''
     WRITE(*,*) ' You can provide options via -c concerning the MPI/Rank configuration:'
     WRITE(*,*) ' -c NRANKS,MPICOM'
     WRITE(*,*) ''
     WRITE(*,*) ' Option Documentation:'
     WRITE(*,*) ' NRANKS = Number of MPI Ranks'
     WRITE(*,*) ' MPICOM = MPI Execution Command: 0 -> mpirun, 1 -> srun, 2 -> runjob'
     WRITE(*,*) ''
   END SUBROUTINE helpmessage

 END PROGRAM IPI_CLIENT
