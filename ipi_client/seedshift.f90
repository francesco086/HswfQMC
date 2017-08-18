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
! oldfname: Old File name
! newfname: New File name
! oldunit: Old File I/O unit
! newunit: New File I/O unit

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: oldfname, newfname
INTEGER, INTENT(IN) :: oldunit, newunit
LOGICAL, INTENT(IN) :: doreplace

CHARACTER(LEN=1000) :: firststr, linestr
INTEGER :: oldst = 0

OPEN(FILE = oldfname, UNIT = oldunit, STATUS = 'old', ACTION = 'read')
OPEN(FILE = newfname, UNIT = newunit, STATUS = 'replace', ACTION = 'write')

! copy header
READ(UNIT = oldunit, FMT = '(A)') linestr
WRITE(UNIT = newunit, FMT = '(A)') trim(linestr)

! backup first seed line
READ(UNIT = oldunit, FMT = '(A)') firststr

DO WHILE (oldst == 0)
    READ(UNIT = oldunit, IOSTAT = oldst, FMT = '(A)') linestr
    IF (oldst == 0) THEN
        WRITE(UNIT = newunit, FMT = '(A)') trim(linestr)
    END IF
ENDDO

! write first seed line to end
WRITE(UNIT = newunit, FMT = '(A)') trim(firststr)

CLOSE(UNIT = oldunit)
CLOSE(UNIT = newunit)

IF (doreplace) THEN
    CALL execute_command_line('mv ' // trim(newfname) // ' ' // trim(oldfname))
END IF

END SUBROUTINE
