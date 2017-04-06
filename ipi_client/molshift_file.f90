PROGRAM molshift_file
IMPLICIT NONE

REAL(8)		::	L(3)
INTEGER		::	N, N_PART

REAL(8), ALLOCATABLE	::	POSARR(:,:), NEWARR(:,:), DISTARR(:)
INTEGER, ALLOCATABLE	::	INVINDARR(:)

REAL(8)		::	distv(3), LR(3)
INTEGER		::	i,j

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
