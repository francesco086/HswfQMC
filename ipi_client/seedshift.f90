! Contains a routine to shift a random seed input file for HswfQMC,
! realised by moving the file's first line to the end of the file.
! 
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
!-------------------------------------------------------------------------


SUBROUTINE seedshift(oldfname, newfname, oldunit, newunit, doreplace)
! Creates a copy of file oldfname, where the first line is moved to the end
!
! oldfname:  Old File name
! newfname:  New File name
! oldunit:   Old File I/O unit
! newunit:   New File I/O unit
! doreplace: Replace the old randomseed file?

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: oldfname, newfname
INTEGER, INTENT(IN) :: oldunit, newunit
LOGICAL, INTENT(IN) :: doreplace

CHARACTER(LEN=1024) :: firststr, linestr

OPEN(FILE = oldfname, UNIT = oldunit, STATUS = 'old', ACTION = 'read', ERR=3)
OPEN(FILE = newfname, UNIT = newunit, STATUS = 'replace', ACTION = 'write', ERR=3)

! copy header
READ(UNIT = oldunit, FMT = '(A)', END=2, ERR=3) linestr
WRITE(UNIT = newunit, FMT = '(A)', ERR=3) trim(linestr)

! backup first seed line
READ(UNIT = oldunit, FMT = '(A)', END=2, ERR=3) firststr

DO
    READ(UNIT = oldunit, FMT = '(A)', END=1, ERR=3) linestr
    WRITE(UNIT = newunit, FMT = '(A)', ERR=3) trim(linestr)

    CYCLE
1   EXIT
ENDDO

! write first seed line to end
WRITE(UNIT = newunit, FMT = '(A)', ERR=3) trim(firststr)

CLOSE(UNIT = oldunit, ERR=3)
CLOSE(UNIT = newunit, ERR=3)

IF (doreplace) THEN
    CALL execute_command_line('mv ' // trim(newfname) // ' ' // trim(oldfname), WAIT=.TRUE.)
END IF

RETURN

2   STOP '[seedshift] Error: Random seed file ended prematurely.'

3   STOP '[seedshift] Error: File I/O error.'

END SUBROUTINE
