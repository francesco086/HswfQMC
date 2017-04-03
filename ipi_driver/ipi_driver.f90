! The main program which runs our driver test case potentials
! 
! Copyright (C) 2013, Joshua More and Michele Ceriotti
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
!
! Currently the potentials implemented are the Lennard-Jones
! potential, the Silvera-Goldman para-hydrogen potential and
! the ideal gas (i.e. no interaction at all)
!
! -------------------------------------------------------------------------
! Modified by Jan Kessler (2015-2017)
! An adaption of the I-Pi driver example to HswfQMC

PROGRAM IPI_DRIVER

    IMPLICIT NONE

    ! SOCKET VARIABLES
    INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
    INTEGER socket, inet, port        ! socket ID & address of the server
    CHARACTER*1024 :: host

    ! COMMAND LINE PARSING
    CHARACTER*1024 :: cmdbuffer, vops
    INTEGER ccmd, vstyle
    LOGICAL verbose
    INTEGER commas(4), par_count      ! stores the index of commas in the parameter string
    INTEGER vpars(4), nargs           ! array to store the parameters of the potential

    ! SOCKET COMMUNICATION BUFFERS
    CHARACTER*12 :: header
    LOGICAL :: isinit=.false., hasdata=.false.
    INTEGER cbuf
    CHARACTER*2048 :: initbuffer      ! it's unlikely a string this large will ever be passed...
    DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)

    ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
    DOUBLE PRECISION              :: box(3), pot, enhelp(3)
    DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:), forces(:,:)
    INTEGER nat, beadid, ios, istep
    DOUBLE PRECISION cell_h(3,3), cell_ih(3,3), virial(3,3)

    ! PARAMETERS CONCERNING VMC (NCPUS,NWOMIN)
    INTEGER ncpus,nwomin
    CHARACTER*1024 strncpu, strbid
    LOGICAL DOSHIFT
    INTEGER, ALLOCATABLE          :: invindarr(:)

    ! ITERATORS
    INTEGER i,k,l

    ! parse the command line parameters
    ! intialize defaults
    ccmd = 0
    inet = 1
    host = "localhost"//achar(0)
    port = 31415
    verbose = .false.
    par_count = 0
    vstyle = -1
    nargs = IARGC()

    DO i = 1, nargs
        CALL GETARG(i, cmdbuffer)
        IF (cmdbuffer == "-u") THEN ! flag for unix socket
            inet = 0
            ccmd = 0
        ELSEIF (cmdbuffer == "-h") THEN ! read the hostname
            ccmd = 1
        ELSEIF (cmdbuffer == "-p") THEN ! reads the port number
            ccmd = 2
        ELSEIF (cmdbuffer == "-m") THEN ! reads the handling of electronic parameter optimization
            ccmd = 3
        ELSEIF (cmdbuffer == "-o") THEN ! reads the parameters
            ccmd = 4
        ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = .true.
        ELSE
            IF (ccmd == 0) THEN
                WRITE(*,*) " Unrecognized command line argument", ccmd
                WRITE(*,*) " SYNTAX: ipi_driver.x [-u] -h hostname -p port -m [alt|sim|fix|nop] " &
                     //"[-o comma,separated,parameters] [-v] "
                WRITE(*,*) ""
                WRITE(*,*) " The -u flag enables the use of unix sockets instead of internet sockets. "
                WRITE(*,*) ""
                WRITE(*,*) " You may have to provide parameters via -o according to the WF optimization mode:"
                WRITE(*,*) " Alternating wavefunction optimization mode (-m alt):  -o NCPU,ISHIFT,NWOMIN "
                WRITE(*,*) " Simultaneous wavefunction optimization mode (-m sim): -o NCPU,ISHIFT "
                WRITE(*,*) " Fixed wavefunction mode (-m fix):                     -o NCPU,ISHIFT "
                WRITE(*,*) " In case of the no operation mode there are no options to be set."
                CALL EXIT(-1)
            ENDIF
            IF (ccmd == 1) THEN
                host = trim(cmdbuffer)//achar(0)
            ELSEIF (ccmd == 2) THEN
                READ(cmdbuffer,*) port
            ELSEIF (ccmd == 3) THEN
                IF (trim(cmdbuffer) == "alt") THEN
                    vstyle = 2
                ELSEIF (trim(cmdbuffer) == "sim") THEN
                    vstyle = 1
                ELSEIF (trim(cmdbuffer) == "fix") THEN
                    vstyle = 0
                ELSEIF (trim(cmdbuffer) == "nop") THEN
                    vstyle = -1
                ELSE
                    WRITE(*,*) " Unrecognized optimization mode ", trim(cmdbuffer)
                    WRITE(*,*) " Use -m [alt|sim|fix|nop] "
                    CALL EXIT(-1)
                ENDIF
            ELSEIF (ccmd == 4) THEN
                par_count = 1
                commas(1) = 0
                DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                    commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                    READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vpars(par_count)
                    par_count = par_count + 1
                ENDDO
                READ(cmdbuffer(commas(par_count)+1:),*) vpars(par_count)
            ENDIF
            ccmd = 0
        ENDIF
    ENDDO
    IF (vstyle == -1) THEN
        isinit = .true.
    ELSEIF (vstyle == 0 .or. vstyle==1) THEN
        IF (par_count /= 2) THEN
            WRITE(*,*) "Error: Wrong number of -o parameters provided. Expected: -o NCPU,ISHIFT"
            CALL EXIT(-1) ! Note that if initialization from the wrapper is implemented this exit should be removed.
        ENDIF
        ncpus = vpars(1)
        IF (vpars(2) > 0) THEN
            doshift = .TRUE.
        ELSE
            doshift = .FALSE.
        END IF
        isinit = .true.
    ELSEIF (vstyle == 2) THEN
        IF (par_count /= 3) THEN
            WRITE(*,*) "Error: Wrong number of -o parameters provided. Expected: -o NCPU,ISHIFT,NWOMIN"
            CALL EXIT(-1) ! Note that if initialization from the wrapper is implemented this exit should be removed.
        ENDIF
        ncpus = vpars(1)
        IF (vpars(2) > 0) THEN
            doshift = .TRUE.
        ELSE
            doshift = .FALSE.
        END IF
        nwomin = vpars(3)
        isinit = .true.
    ENDIF

    IF (verbose) THEN
        WRITE(*,*) " DRIVER - Connecting to host ", trim(host)
        IF (inet > 0) THEN
            WRITE(*,*) " on port ", port, " using an internet socket."
        ELSE
            WRITE(*,*) " using an UNIX socket."
        ENDIF
    ENDIF

    ! Calls the interface to the C sockets to open a communication channel
    CALL open_socket(socket, inet, port, host)

    nat = -1
    beadid = 0
    istep = 0
    DO WHILE (.TRUE.) ! Loops forever (or until the wrapper ends!)

        ! Reads from the socket one message header
        CALL readbuffer(socket, header, MSGLEN)
        IF (verbose) WRITE(*,*) " Message from server: ", trim(header)

        IF (trim(header) == "STATUS") THEN
            ! The wrapper is inquiring on what we are doing
            IF (.NOT. isinit) THEN
                CALL writebuffer(socket,"NEEDINIT    ",MSGLEN)  ! Signals that we need initialization data
            ELSEIF (hasdata) THEN
                CALL writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
            ELSE
                CALL writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
            ENDIF
        ELSEIF (trim(header) == "INIT") THEN     ! The driver is kindly providing a string for initialization
            CALL readbuffer(socket, cbuf, 4)
            CALL readbuffer(socket, initbuffer, cbuf)
            IF (verbose) WRITE(*,*) " Initializing system from wrapper, using ", trim(initbuffer)
            isinit=.true. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver
        ELSEIF (trim(header) == "POSDATA") THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!

            ! Parses the flow of data from the socket
            CALL readbuffer(socket, cell_h,  9*8)  ! Cell matrix
            CALL readbuffer(socket, cell_ih, 9*8)  ! Inverse of the cell matrix (so we don't have to invert it every time here)

            ! The wrapper uses atomic units for everything, and row major storage.
            ! At this stage one should take care that everything is converted in the
            ! units and storage mode used in the driver.
            cell_h = transpose(cell_h)
            cell_ih = transpose(cell_ih)
            ! We assume an upper triangular cell-vector matrix
            box(1) = cell_h(1,1)
            box(2) = cell_h(2,2)
            box(3) = cell_h(3,3)

            CALL readbuffer(socket, cbuf, 4)       ! The number of atoms in the cell
            IF (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation, so only does this once
                nat = cbuf
                IF (verbose) WRITE(*,*) " Allocating buffer and data arrays, with ", nat, " atoms"
                ALLOCATE(msgbuffer(3*nat))
                ALLOCATE(atoms(3,nat))
                ALLOCATE(forces(3,nat))
                IF (doshift) ALLOCATE(invindarr(nat))
            ENDIF

            CALL readbuffer(socket, msgbuffer, nat*3*8) ! Now receive the current atom positions
            DO i = 1, nat
                atoms(:,i) = msgbuffer(3*(i-1)+1:3*i)
            ENDDO

            CALL readbuffer(socket, beadid, 4) ! And finally the replica/bead id
            
            pot = 0
            forces = 0
            virial = 0
            ios = 0
            enhelp = 0

            istep = istep + 1

            IF (verbose) WRITE(*,*) " Force calculation step ", istep
            IF (verbose) WRITE(*,*) " Calculating forces for bead ", beadid
            
            IF (vstyle > -1) THEN 

                IF (doshift) CALL molshift(nat, box, atoms, invindarr)
               
                OPEN (UNIT=20, FILE='reticolo/rp_now.d', STATUS='REPLACE', ACTION='WRITE')
                write(20,*) box
                DO i = 1, nat
                    write(20,*) atoms(:,i)
                ENDDO
                CLOSE(20)

                call execute_command_line('python seedfile_shift.py', WAIT = .true.)

		IF (vstyle > 0) THEN
                    write(strbid, *) beadid
                    strbid = trim(adjustl(strbid))
		    call execute_command_line('cp ../SR_wf.dir/'//strbid//' wf_now.d', WAIT = .true.)
		END IF
                
                write (strncpu, *) ncpus
                strncpu = trim(adjustl(strncpu))
	        !call execute_command_line('mpirun -np '//strncpu//' HswfQMC_exe', WAIT = .true.)
                call execute_command_line('runjob --np '//strncpu &
                //' --exe /homea/hpb01/hpb015/HswfQMC/HswfQMC_exe --ranks-per-node 64', WAIT = .true.)

		IF (vstyle > 0) THEN
		    call execute_command_line('cp ottimizzazione/SR_wf.d ../SR_wf.dir/'//strbid, WAIT = .true.)
		END IF
                
                OPEN (UNIT=20, FILE='reticolo/LagrDyn_Frp-0000.d',ACTION='READ')
                DO i = 1, nat
                   read(20, *) forces(:,i)
                ENDDO
                forces(:,:) = forces(:,:) *0.5*nat !HswfQMC output values are per atom, convert to Hartree
                CLOSE(20)
  	        
                OPEN (UNIT=20, FILE='ottimizzazione/SR_energies.dat',ACTION='READ')
                DO while(ios.eq.0)
                    read(20,*,iostat=ios) enhelp
                ENDDO
                pot = enhelp(2) *0.5*nat
                CLOSE(20)

                IF (doshift) THEN
                    CALL molshift_back(nat, atoms, invindarr)
                    CALL molshift_back(nat, forces, invindarr)
                END IF
                 
                DO i = 1, nat
                    DO k = 1, 3
                        DO l = k, 3
                            ! Only the upper triangular elements calculated.
                            virial(k,l) = virial(k,l) + forces(k,i)*atoms(l,i)
                        ENDDO
                    ENDDO
                ENDDO
            END IF

            IF (verbose) WRITE(*,*) " Calculated energy is ", pot
            hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper

        ELSEIF (trim(header) == "GETFORCE") THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

            ! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            DO i = 1, nat
                msgbuffer(3*(i-1)+1:3*i) = forces(:,i)
            ENDDO
            virial = transpose(virial)

            CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)
            CALL writebuffer(socket,pot,8)  ! Writing the potential
            CALL writebuffer(socket,nat,4)  ! Writing the number of atoms
            CALL writebuffer(socket,msgbuffer,3*nat*8) ! Writing the forces
            CALL writebuffer(socket,virial,9*8)  ! Writing the virial tensor, NOT divided by the volume
            cbuf = 7 ! Size of the "extras" string
            CALL writebuffer(socket,cbuf,4) ! This would write out the "extras" string, but in this case we only use a dummy string.
            CALL writebuffer(socket,"nothing",7)

            hasdata = .false.
        ELSE
            WRITE(*,*) " Unexpected header ", header
            CALL EXIT(-1)
        ENDIF
    ENDDO
    IF (nat > 0) THEN
        DEALLOCATE(atoms, forces, msgbuffer)
        IF (doshift) DEALLOCATE(invindarr)
    END IF

END PROGRAM IPI_DRIVER
