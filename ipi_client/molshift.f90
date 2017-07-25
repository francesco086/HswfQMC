! A routine which takes an array of coordinates, finds pairs of nearest
! neighbors and sorts the array in a way that HswfQMC assigns spin up/down
! to pairs of H atoms (H2 molecules) appropriately.
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
! -------------------------------------------------------------------------

SUBROUTINE molshift(N, L, POSARR, INVINDARR)

IMPLICIT NONE

! Passed variables
INTEGER, INTENT(IN)       ::      N
REAL(8), INTENT(IN)       ::      L(3)
REAL(8), INTENT(INOUT)    ::      POSARR(3,N)
INTEGER, INTENT(OUT)      ::      INVINDARR(N)


! Working variables
INTEGER		::	N_PART
REAL(8)	        ::	NEWARR(3,N), DISTARR(N)
INTEGER         ::      INDARR(N)
LOGICAL         ::	ISPAIR(N)


! Helper variables
REAL(8)		::	distv(3), LR(3)
INTEGER		::	i,j,inew,jnew

N_PART = N/2

INVINDARR(:) = -1

NEWARR(:,:) = 0.d0
INDARR(:) = -1
ISPAIR(:) = .FALSE.

LR(:) = 1.d0 / L(:)

inew=0
DO i=1,N

    IF (ISPAIR(i)) THEN
        CYCLE
    END IF
    IF (inew >= N_PART) THEN
        EXIT
    END IF
    inew = inew + 1
    DISTARR(:) = -1.d0

    DO j=1,N

        IF (.NOT.ISPAIR(j) .AND. i/=j) THEN
            distv(:) = POSARR(:,j) - POSARR(:,i)
            distv(:) = distv(:) - L(:) * NINT(distv(:) * LR(:)) 
            DISTARR(j) = SUM(distv(:)**2) 
        END IF 

    ENDDO

    j = MINLOC(DISTARR, DIM=1, MASK=(DISTARR>=0))
    jnew = inew + N_PART 

    INDARR(i) = inew
    INDARR(j) = jnew

    ISPAIR(i) = .TRUE.
    ISPAIR(j) = .TRUE.
ENDDO

DO i=1,N
    NEWARR(:,INDARR(i)) = POSARR(:, i)  
ENDDO

! Write result to output arrays

POSARR(:,:) = NEWARR(:,:)

DO i=1,N
    INVINDARR(INDARR(i)) = i
ENDDO

! Tests -> Errors

DO i=1,N
   IF (INDARR(i) == -1) THEN
       WRITE(*,*) '(molshift) Error: Index Array entry equal -1 detected. Something went wrong!'
   END IF
ENDDO

DO i=1,N
   IF (INVINDARR(i) == -1) THEN
       WRITE(*,*) '(molshift) Error: Inverse Index Array entry equal -1 detected. Something went wrong!'
   END IF
ENDDO

END SUBROUTINE molshift


SUBROUTINE molshift_back(N, FRCARR, INVINDARR)

IMPLICIT NONE

! Passed variables
INTEGER, INTENT(IN)         ::      N
REAL(8), INTENT(INOUT)      ::      FRCARR(3,N)
INTEGER, INTENT(IN)         ::      INVINDARR(N)

! Working variables
REAL(8)         ::      NEWARR(3,N)

! Helper variables
INTEGER         ::      i

DO i=1,N
    NEWARR(:,INVINDARR(i)) = FRCARR(:,i)  
ENDDO

FRCARR(:,:) = NEWARR(:,:)

END SUBROUTINE
