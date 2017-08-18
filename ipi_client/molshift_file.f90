! This program applies the molshift routine to coordinates read from a given
! input file named 'molxyz_toshift' and writes the result to 'molxyz_shifted'
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

PROGRAM molshift_file
IMPLICIT NONE

REAL(8)    ::    L(3)
INTEGER    ::    N, N_PART

REAL(8), ALLOCATABLE    ::    POSARR(:,:), NEWARR(:,:), DISTARR(:)
INTEGER, ALLOCATABLE    ::    INVINDARR(:)

REAL(8)    ::    distv(3), LR(3)
INTEGER    ::    i,j

N = 64
N_PART = N/2

ALLOCATE(POSARR(3,N), NEWARR(3,N), DISTARR(N), INVINDARR(N))

NEWARR(:,:) = 0.d0

OPEN(20,FILE='molxyz_toshift')
READ(20,*) L(:)
DO i=1,N
    READ(20,*) POSARR(:,i) 
ENDDO
CLOSE(20)

LR(:) = 1.d0 / L(:)

NEWARR(:,:) = POSARR(:,:)

CALL molshift(N, L, NEWARR, INVINDARR)

OPEN(21,FILE='molxyz_shifted')
WRITE(21,*) L(:)
DO i=1,N
    WRITE(21,*) NEWARR(:,i)
ENDDO
CLOSE(21)

! Test Output

WRITE(*,*) INVINDARR

WRITE(*,*)

DISTARR = 0.d0
DO i=1,N_PART
j = i + N_PART
distv(:) = POSARR(:,j)-POSARR(:,i)
distv(:) = distv(:) -  L(:) * NINT(distv(:) * LR(:))
DISTARR(i) = SQRT(SUM(distv(:)**2))
ENDDO
WRITE(*,*) DISTARR(1:N_PART)

WRITE(*,*)

DISTARR = 0.d0
DO i=1,N_PART
j = i + N_PART
distv(:) = NEWARR(:,j)-NEWARR(:,i)
distv(:) = distv(:) -  L(:) * NINT(distv(:) * LR(:))
DISTARR(i) = SQRT(SUM(distv(:)**2))
ENDDO
WRITE(*,*) DISTARR(1:N_PART)

DEALLOCATE(POSARR, NEWARR, DISTARR, INVINDARR)

END PROGRAM molshift_file
