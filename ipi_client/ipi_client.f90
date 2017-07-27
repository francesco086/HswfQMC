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
   INTEGER ccmd, vstyle
   INTEGER verbose
   INTEGER commas(2), par_count_o, par_count_c      ! stores the index of commas in the parameter string
   INTEGER vpars_o(2), vpars_c(2)         ! array to store the parameters of the potential

   ! SOCKET COMMUNICATION BUFFERS
   CHARACTER(LEN=12) :: header
   LOGICAL :: isinit, hasdata
   INTEGER cbuf, rid, ios, istep
   CHARACTER(LEN=2048) :: initbuffer      ! it's unlikely a string this large will ever be passed...
   DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)

   ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
   INTEGER nat
   DOUBLE PRECISION pot, volume, box(3), enhelp(3)
   DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:), forces(:,:)
   DOUBLE PRECISION cell_h(3,3), cell_ih(3,3), virial(3,3), mtxbuf(9)

   ! PARAMETERS CONCERNING VMC (NRANKS,NWOMIN)
   INTEGER nranks,mpicom,nwomin
   CHARACTER*1024 strnranks, strrid
   LOGICAL domshift, dolimit, notlimit
   INTEGER, ALLOCATABLE :: invindarr(:)
   DOUBLE PRECISION :: flimit

   ! ITERATORS
   INTEGER i, k, l

   ! intialize parameter defaults
   inet = 1
   verbose = 0
   host = "localhost"//achar(0)
   port = 54321
   vstyle = -1
   domshift = .FALSE.
   nwomin = 0
   nranks = 1
   mpicom = 0
   dolimit = .FALSE.

   ! initialize control variables
   isinit = .FALSE.
   hasdata = .FALSE.
   ccmd = 0
   par_count_o = 0
   par_count_c = 0

   nat = -1
   rid = 0
   istep = 0

   ! parse the command line parameters
   DO i = 1, COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(i, cmdbuffer)
      IF (cmdbuffer == "-u") THEN ! flag for unix socket
         inet = 0
      ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
         verbose = 1
      ELSEIF (cmdbuffer == "-vv") THEN ! flag for verbose standard output
         verbose = 2
      ELSEIF (cmdbuffer == "-h") THEN ! read the hostname
         ccmd = 1
      ELSEIF (cmdbuffer == "-p") THEN ! reads the port number
         ccmd = 2
      ELSEIF (cmdbuffer == "-m") THEN ! reads the mode of electronic parameter optimization
         ccmd = 3
      ELSEIF (cmdbuffer == "-o") THEN ! reads the options for mode
         ccmd = 4
      ELSEIF (cmdbuffer == "-c") THEN ! reads the MPI/Rank configuration
         ccmd = 5
      ELSEIF (cmdbuffer == "-l") THEN ! reads the force amplitude limit
         ccmd = 6
      ELSE
         IF (ccmd == 0) THEN
            WRITE(*,*) " Unrecognized command line argument: ", cmdbuffer
            CALL helpmessage
            STOP "ENDED"
         ELSEIF (ccmd == 1) THEN
            host = trim(cmdbuffer)//achar(0)
         ELSEIF (ccmd == 2) THEN
            READ(cmdbuffer,*) port
         ELSEIF (ccmd == 3) THEN
            IF (trim(cmdbuffer) == "alt") THEN
               vstyle = 3
            ELSEIF (trim(cmdbuffer) == "sim") THEN
               vstyle = 2
            ELSEIF (trim(cmdbuffer) == "fix") THEN
               vstyle = 1
            ELSEIF (trim(cmdbuffer) == "ffs") THEN
               vstyle = 0
            ELSEIF (trim(cmdbuffer) == "nop") THEN
               vstyle = -1
            ELSE
               WRITE(*,*) " Unrecognized wavefunction optimization mode ", trim(cmdbuffer)
               WRITE(*,*) " Use -m [alt|sim|fix|ffs|nop] "
               STOP "ENDED"
            ENDIF
         ELSEIF (ccmd == 4) THEN
            par_count_o = 1
            commas(1) = 0
            DO WHILE (index(cmdbuffer(commas(par_count_o)+1:), ',') > 0)
               commas(par_count_o + 1) = index(cmdbuffer(commas(par_count_o)+1:), ',') + commas(par_count_o)
               READ(cmdbuffer(commas(par_count_o)+1:commas(par_count_o + 1)-1),*) vpars_o(par_count_o)
               par_count_o = par_count_o + 1
            ENDDO
            READ(cmdbuffer(commas(par_count_o)+1:),*) vpars_o(par_count_o)
         ELSEIF (ccmd == 5) THEN
            par_count_c = 1
            commas(1) = 0
            DO WHILE (index(cmdbuffer(commas(par_count_c)+1:), ',') > 0)
               commas(par_count_c + 1) = index(cmdbuffer(commas(par_count_c)+1:), ',') + commas(par_count_c)
               READ(cmdbuffer(commas(par_count_c)+1:commas(par_count_c + 1)-1),*) vpars_c(par_count_c)
               par_count_c = par_count_c + 1
            ENDDO
            READ(cmdbuffer(commas(par_count_c)+1:),*) vpars_c(par_count_c)
         ELSEIF (ccmd == 6) THEN
            READ(cmdbuffer(:),*) flimit ! force amplitude limit (dirty fix)
            IF (flimit > 0) THEN
               dolimit = .TRUE.
            END IF
         ENDIF
         ccmd = 0
      ENDIF
   ENDDO

   IF (vstyle < 1) THEN
      CONTINUE
   ELSEIF (vstyle < 3) THEN
      IF (par_count_o /= 1) THEN
         WRITE(*,*) "Error: Wrong number of -o parameters provided. Expected: -o ISHIFT"
         STOP "ENDED"
      ENDIF
      IF (vpars_o(1) > 0) THEN
         domshift = .TRUE.
      ELSE
         domshift = .FALSE.
      END IF
   ELSEIF (vstyle == 3) THEN
      IF (par_count_o /= 2) THEN
         WRITE(*,*) "Error: Wrong number of -o parameters provided. Expected: -o ISHIFT,NWOMIN"
         STOP "ENDED"
      ENDIF
      IF (vpars_o(1) > 0) THEN
         domshift = .TRUE.
      ELSE
         domshift = .FALSE.
      END IF
      nwomin = vpars_o(2)
   ENDIF

   IF (vstyle < 2) THEN
      isinit = .TRUE.
   ELSE
      isinit = .FALSE. ! WF optimization enabled -> we need to know our bead id
   END IF
   IF (par_count_c == 1) THEN
      nranks = vpars_c(1)
   ELSEIF (par_count_c == 2) THEN
      nranks = vpars_c(1)
      mpicom = vpars_c(2)
   ELSE
      WRITE(*,*) "Error: Wrong number of -c parameters provided. Expected: -c NRANKS or -c NRANKS,MPICOM"
      STOP "ENDED"
   END IF

   IF (verbose > 0) THEN
      WRITE(*,*) " DRIVER - Connecting to host ", trim(host)
      IF (inet > 0) THEN
         WRITE(*,*) " on port ", port, " using an internet socket."
      ELSE
         WRITE(*,*) " using an UNIX socket."
      ENDIF
   ENDIF

   ! Calls the interface to the POSIX sockets library to open a communication channel
   CALL open_socket(socket, inet, port, host)

   DO WHILE (.true.) ! Loops forever (or until the wrapper ends!)

      ! Reads from the socket one message header
      CALL readbuffer(socket, header, MSGLEN)
      IF (verbose > 0) WRITE(*,*) " Message from server: ", trim(header)

      IF (trim(header) == "STATUS") THEN
         ! The wrapper is inquiring on what we are doing
         IF (.not. isinit) THEN
            CALL writebuffer(socket,"NEEDINIT    ",MSGLEN)  ! Signals that we need initialization data
            IF (verbose > 1) WRITE(*,*) "    !write!=> ", "NEEDINIT    "
         ELSEIF (hasdata) THEN
            CALL writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
            IF (verbose > 1) WRITE(*,*) "    !write!=> ", "HAVEDATA    "
         ELSE
            CALL writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
            IF (verbose > 1) WRITE(*,*) "    !write!=> ", "READY       "
         ENDIF
      ELSEIF (trim(header) == "INIT") THEN     ! The driver is kindly providing a string for initialization
         CALL readbuffer(socket, rid)
         IF (verbose > 1) WRITE(*,*) "    !read!=> RID: ", rid
         CALL readbuffer(socket, cbuf)
         IF (verbose > 1) WRITE(*,*) "    !read!=> init_lenght: ", cbuf
         CALL readbuffer(socket, initbuffer, cbuf)
         IF (verbose > 1) WRITE(*,*) "    !read!=> init_string: ", cbuf
         IF (verbose > 0) WRITE(*,*) " Initializing system from wrapper, using ", trim(initbuffer)
         isinit=.TRUE. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver
      ELSEIF (trim(header) == "POSDATA") THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!

         ! Parses the flow of data from the socket
         CALL readbuffer(socket, mtxbuf, 9)  ! Cell matrix
         IF (verbose > 1) WRITE(*,*) "    !read!=> cell: ", mtxbuf
         cell_h = RESHAPE(mtxbuf, (/3,3/))
         CALL readbuffer(socket, mtxbuf, 9)  ! Inverse of the cell matrix (so we don't have to invert it every time here)
         IF (verbose > 1) WRITE(*,*) "    !read!=> cell-1: ", mtxbuf
         cell_ih = RESHAPE(mtxbuf, (/3,3/))

         ! We assume an upper triangular cell-vector matrix
         volume = cell_h(1,1)*cell_h(2,2)*cell_h(3,3)
         box(1) = cell_h(1,1)
         box(2) = cell_h(2,2)
         box(3) = cell_h(3,3)

         CALL readbuffer(socket, cbuf)       ! The number of atoms in the cell
         IF (verbose > 1) WRITE(*,*) "    !read!=> cbuf: ", cbuf
         IF (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation, so only does this once
            nat = cbuf
            IF (verbose > 0) WRITE(*,*) " Allocating buffer and data arrays, with ", nat, " atoms"
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
         IF (verbose > 1) WRITE(*,*) "    !read!=> positions: ", msgbuffer
         DO i = 1, nat
            atoms(:,i) = msgbuffer(3*(i-1)+1:3*i)
         ENDDO

         pot = 0.0d0
         virial = 0.0d0

         enhelp = 0
         ios = 0
         istep = istep + 1
         IF (verbose > 0) THEN
            WRITE(*,*) " Force calculation step ", istep
            WRITE(*,*) " Calculating forces for bead ", rid
         END IF

         IF (vstyle > -1) THEN

            IF (domshift) CALL molshift(nat, box, atoms, invindarr)

            OPEN (UNIT=20, FILE='reticolo/rp_now.d', STATUS='REPLACE', ACTION='WRITE')
            WRITE(20,*) box
            DO i = 1, nat
               WRITE(20,*) atoms(:,i)
            ENDDO
            CLOSE(20)

            IF (vstyle > 1) THEN
               WRITE(strrid, *) rid
               strrid = trim(adjustl(strrid))
               CALL execute_command_line('cp ../SR_wf.dir/'//strrid//' wf_now.d', WAIT = .true.)
            END IF

            WRITE (strnranks, *) nranks
            strnranks = trim(adjustl(strnranks))

            DO WHILE (.TRUE.)

               IF (mpicom == 1) THEN
                  CALL execute_command_line('srun -n '//strnranks//' HswfQMC_exe', WAIT = .true.)
               ELSE IF (mpicom == 2) THEN
                  CALL execute_command_line('runjob --np '//strnranks &
                       //' --exe /homea/hpb01/hpb015/HswfQMC/HswfQMC_exe --ranks-per-node 64', WAIT = .true.)
               ELSE
                  CALL execute_command_line('mpirun -np '//strnranks//' HswfQMC_exe', WAIT = .true.)
               END IF

               IF (vstyle > 1) THEN
                  CALL execute_command_line('cp ottimizzazione/SR_wf.d ../SR_wf.dir/'//strrid, WAIT = .true.)
               END IF

               OPEN (UNIT=20, FILE='reticolo/LagrDyn_Frp-0000.d',ACTION='READ')
               DO i = 1, nat
                  READ(20, *) forces(:,i)
               ENDDO
               forces(:,:) = forces(:,:) *0.5*nat !HswfQMC output values are per atom, convert to Hartree
               CLOSE(20)

               notlimit = .TRUE.
               IF (dolimit) THEN
                  DO i = 1, nat
                     IF (NORM2(forces(:,i)) > flimit) THEN
                        notlimit = .FALSE.
                     END IF
                  ENDDO
               END IF

               IF (notlimit) THEN
                  IF (vstyle == 0) THEN
                     CALL execute_command_line('seedfile_shift.py', WAIT = .true.)
                  END IF
                  EXIT
               ELSE
                  WRITE(*,*) "Force was above limit. Recalculating with new seed..."
                  CALL execute_command_line('seedfile_shift.py', WAIT = .true.)
               END IF

            ENDDO

            IF (vstyle > 1) THEN
               CALL execute_command_line('cp ottimizzazione/SR_wf.d ../SR_wf.dir/'//strrid, WAIT = .true.)
               isinit=.FALSE. ! we want to make sure that the beadid didn't change next time
            END IF

            OPEN (UNIT=20, FILE='ottimizzazione/SR_energies.dat',ACTION='READ')
            DO WHILE(ios.eq.0)
               READ(20,*,iostat=ios) enhelp
            ENDDO
            pot = enhelp(2) *0.5*nat
            CLOSE(20)

            IF (domshift) THEN
               CALL molshift_back(nat, atoms, invindarr)
               CALL molshift_back(nat, forces, invindarr)
            END IF

            DO i = 1, nat
               DO k = 1, 3
                  DO l = k, 3
                     virial(k,l) = virial(k,l) + atoms(k,i)*forces(l,i)
                  ENDDO
               ENDDO
            ENDDO
         END IF

         IF (verbose > 0) WRITE(*,*) " Calculated energy is ", pot

         hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper
      ELSEIF (trim(header) == "GETFORCE") THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

         DO i = 1, nat
            msgbuffer(3*(i-1)+1:3*i) = forces(:,i)
         ENDDO

         CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)
         IF (verbose > 1) WRITE(*,*) "    !write!=> ", "FORCEREADY  "
         CALL writebuffer(socket,pot)  ! Writing the potential
         IF (verbose > 1) WRITE(*,*) "    !write!=> pot: ", pot
         CALL writebuffer(socket,nat)  ! Writing the number of atoms
         IF (verbose > 1) WRITE(*,*) "    !write!=> nat:", nat
         CALL writebuffer(socket,msgbuffer,3*nat) ! Writing the forces
         IF (verbose > 1) WRITE(*,*) "    !write!=> forces:", msgbuffer
         CALL writebuffer(socket,reshape(virial,(/9/)),9)  ! Writing the virial tensor, NOT divided by the volume
         IF (verbose > 1) WRITE(*,*) "    !write!=> strss: ", reshape(virial,(/9/))


         cbuf = 7 ! Size of the "extras" string
         CALL writebuffer(socket,cbuf) ! This would write out the "extras" string, but in this case we only use a dummy string.
         IF (verbose > 1) WRITE(*,*) "    !write!=> extra_lenght: ", cbuf
         CALL writebuffer(socket,"nothing",7)
         IF (verbose > 1) WRITE(*,*) "    !write!=> extra: nothing"

         hasdata = .false.
      ELSE
         WRITE(*,*) " Unexpected header ", header
         STOP "ENDED"
      ENDIF
   ENDDO
   IF (nat > 0) THEN
      DEALLOCATE(atoms, forces, msgbuffer)
      IF (domshift) DEALLOCATE(invindarr)
   END IF
 CONTAINS
   SUBROUTINE helpmessage
     ! Help banner
     WRITE(*,*) " SYNTAX: ipi_driver.x [-u] [-v] -h hostname -p port -m [alt|sim|fix|ffs|nop] " &
          //"[-o comma,separated,options] [-c comma,separated,options] "
     WRITE(*,*) ""
     WRITE(*,*) " The -u flag enables the use of unix sockets instead of internet sockets. "
     WRITE(*,*) ""
     WRITE(*,*) " The -v flag enables verbose mode. "
     WRITE(*,*) ""
     WRITE(*,*) " You may want to provide options via -o according to the wavefunction optimization mode:"
     WRITE(*,*) " Fixed wavefunction mode (-m fix):                     -o ISHIFT "
     WRITE(*,*) " Simultaneous wavefunction optimization mode (-m sim): -o ISHIFT "
     WRITE(*,*) " Alternating wavefunction optimization mode (-m alt):  -o ISHIFT,NWOMIN (not implemented yet!)"
     WRITE(*,*) " In case of the no operation mode (-m nop) and fixed force sampling (-m ffs) there are no options to be set."
     WRITE(*,*) ""
     WRITE(*,*) " Option Documentation:"
     WRITE(*,*) " ISHIFT = Molecular Partner Shifting Flag. 0 -> disabled, 1 -> enabled"
     WRITE(*,*) " NWOMIN = Number of wavefunction optimizations without new minimum before termination."
     WRITE(*,*) ""
     WRITE(*,*) " You can provide options via -c concerning the MPI/Rank configuration:"
     WRITE(*,*) " -c NRANKS,MPICOM"
     WRITE(*,*) ""
     WRITE(*,*) " Option Documentation:"
     WRITE(*,*) " NRANKS = Number of MPI Ranks"
     WRITE(*,*) " MPICOM = MPI Execution Command: 0 -> mpirun, 1 -> srun, 2 -> runjob"
     WRITE(*,*) ""
   END SUBROUTINE helpmessage

 END PROGRAM IPI_CLIENT
