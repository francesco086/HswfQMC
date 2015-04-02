MODULE fermi_k
	IMPLICIT NONE
   
! - - - USEFUL CONSTANTS - - - !

	REAL (KIND=8), PRIVATE, PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0


! - - - DECLARATION OF THE CLASS KWaVe - - - !

   TYPE KWaVe
	
		REAL(KIND=8), ALLOCATABLE :: k(:,:)                    !k-vectors. They will be of dimension (0:num_dim,1:num_part)
                                                             !The elements (0,i) are reserved for the squared modules k^2 of the k-vector i
		INTEGER :: num_dim                                     !Number of dimensions
		INTEGER :: num_part                                    !Number of particles
      LOGICAL, PRIVATE :: flag_init=.FALSE.                  !Flag indicating whether the KWaVe object has been initialized or not
      REAL(KIND=8), ALLOCATABLE :: L(:)
		REAL(KIND=8), PRIVATE, ALLOCATABLE :: k_big(:,:)                !Extended grid for k, useful for Twist Averaged Boundary Conditions (TABC)
                                                             !It will have dimension (0:num_dim,num_big)
      INTEGER, PRIVATE :: num_big                                     !Number of k-vectors in k_big
	
	CONTAINS

      PROCEDURE :: initializeKWaVe, deallocateKWaVe   !Initializer and destructor
      PROCEDURE :: buildBoxK                          !Build the the wave vectors for a fermi gas in a box with sides L
      PROCEDURE :: putInOrder                         !Reorder the k-vectors k according to their module, from small values to large values
      PROCEDURE :: shuffleK                           !Shuffle the k-vectors k, i.e. randomly exchange them if they have the same module
      PROCEDURE, PRIVATE :: putInOrderKBig                     !Reorder the k-vectors k_big according to their module
      PROCEDURE, PRIVATE :: shuffleKBig                        !Shuffle the k-vectors k_big
      PROCEDURE, PRIVATE :: copyFromKBig                       !Copy the first num_part k-vectors from k_big to k
      PROCEDURE :: twistK                             !Apply Twist Averaged Boundary Conditions

	END TYPE KWaVe
   

! - - - SUBROUTINES USED INTERNALLY  - - - !

   PRIVATE :: countBoxK, provideBoxK


CONTAINS


! - - - INITIALIZER - - - !
	SUBROUTINE initializeKWaVe(k,n,d)
		IMPLICIT NONE
		CLASS(KWaVe), TARGET :: k
		INTEGER, INTENT(IN) :: n, d
      
      !check that the target object k has not been previously initialized
      IF (k%flag_init) THEN
         PRINT *, "Error called by: quantum_momenta > initializeKWaVe()"
         PRINT *, "You are trying to initialize an already initialized object."
         STOP
      END IF
		
      !Set num_dim and num_part
      k%num_dim=d
		k%num_part=n
      !Allocate the k vectors
		ALLOCATE(k%k(0:d,1:n))
		k%k=0.d0
      !Succesfully initialized
      k%flag_init=.TRUE.

	END SUBROUTINE initializeKWaVe
	
	
! - - - METHODS OF THE CLASS KWaVe - - - !

	SUBROUTINE buildBoxK(k,L)
		!Build the the wave vectors for a fermi gas in a box with sides L
		IMPLICIT NONE
		CLASS(KWaVe), TARGET :: k
		REAL(KIND=8), INTENT(IN) :: L(1:k%num_dim)
      INTEGER :: i1, i2
		INTEGER :: numk
      REAL(KIND=8) :: k2, kmin

      !check that the target object k has been initialized
      IF (.NOT. k%flag_init) THEN
         PRINT *, "Error called by: quantum_momenta > buildBoxK()"
         PRINT *, "You are trying to build the k-vectors before initializing it."
         STOP
      END IF

      !store the L vector
      ALLOCATE(k%L(1:k%num_dim))
      k%L=L
		
      !determine the number of k required to close a shell numk, which is larger than num_part
      numk=0
      kmin=2.d0*PI/MAXVAL(L)
      k2=0.d0
      DO WHILE(numk<k%num_part)
         numk=countBoxK(k2,k%num_dim,L)
         k2=k2+kmin
      END DO

      !set the number of k-vectors num_big for k_big, useful for TABC
      k2=k2-kmin+2.d0*PI/MINVAL(L)
      k%num_big=countBoxK(k2,k%num_dim,L)
      
      !allocate k_big, which will be used to buil the k
      IF (ALLOCATED(k%k_big)) DEALLOCATE(k%k_big) 
      ALLOCATE(k%k_big(0:k%num_dim,1:k%num_big))
      
      !set the k_big vectors
      CALL provideBoxK(k2,k%num_dim,L,k%num_big,k%k_big)
      !put in order the k of k_big according to the 0th element
      CALL k%putInOrderKBig()
      !shuffle them for isotropy
      CALL k%shuffleKBig(10)
      !assign the first num_part k of k_big to k
      CALL k%copyFromKBig()

	END SUBROUTINE buildBoxK
	
	
   SUBROUTINE putInOrder(k)
      IMPLICIT NONE
      CLASS(KWaVe), TARGET :: k
      REAL(KIND=8) :: foo(0:k%num_dim)
      INTEGER :: i1
      
      !check that the target object k has been initialized
      IF (.NOT. k%flag_init) THEN
         PRINT *, "Error called by: quantum_momenta > putInOrder()"
         PRINT *, "You first need to initialize the object."
         STOP
      END IF
		
      i1=0
      DO WHILE (i1<k%num_part-1)
         i1=i1+1
         IF (k%k(0,i1)>k%k(0,i1+1)) THEN
            foo(0:k%num_dim)=k%k(0:k%num_dim,i1)
            k%k(0:k%num_dim,i1)=k%k(0:k%num_dim,i1+1)
            k%k(0:k%num_dim,i1+1)=foo(0:k%num_dim)
            i1=MAX(i1-2,0)
         END IF
      END DO

   END SUBROUTINE putInOrder
   
   
   SUBROUTINE putInOrderKBig(k)
      IMPLICIT NONE
      CLASS(KWaVe), TARGET :: k
      REAL(KIND=8) :: foo(0:k%num_dim)
      INTEGER :: i1
      
      !check that k_big has been allocated
      IF (.NOT. ALLOCATED(k%k_big)) THEN
         PRINT *, "Error called by: quantum_momenta > putInOrderKBig()"
         PRINT *, "The k-vectors k_big do not exist. You first need to build them."
         STOP
      END IF
		
      i1=0
      DO WHILE (i1<k%num_big-1)
         i1=i1+1
         IF (k%k_big(0,i1)>k%k_big(0,i1+1)) THEN
            foo(0:k%num_dim)=k%k_big(0:k%num_dim,i1)
            k%k_big(0:k%num_dim,i1)=k%k_big(0:k%num_dim,i1+1)
            k%k_big(0:k%num_dim,i1+1)=foo(0:k%num_dim)
            i1=MAX(i1-2,0)
         END IF
      END DO

   END SUBROUTINE putInOrderKBig


   SUBROUTINE shuffleK(k,n_iterations)
      IMPLICIT NONE
      CLASS(KWaVe), TARGET :: k
      REAL(KIND=8) :: foo(0:k%num_dim)
      INTEGER, INTENT(IN) :: n_iterations
      INTEGER :: i1, i2
      REAL :: eta
      
      !check that the target object k has been initialized
      IF (.NOT. k%flag_init) THEN
         PRINT *, "Error called by: quantum_momenta > shuffleK()"
         PRINT *, "You first need to initialize the object."
         STOP
      END IF
		
      DO i2 = 1, n_iterations, 1
         i1=0
         DO WHILE (i1<k%num_part-1)
            i1=MAX(i1,0)
            i1=i1+1
            IF (k%k(0,i1)==k%k(0,i1+1)) THEN
               CALL RANDOM_NUMBER(eta)
               IF (eta>0.5) THEN
                   foo(0:k%num_dim)=k%k(0:k%num_dim,i1)
                   k%k(0:k%num_dim,i1)=k%k(0:k%num_dim,i1+1)
                   k%k(0:k%num_dim,i1+1)=foo(0:k%num_dim)
               END IF
            END IF
         END DO
      END DO

   END SUBROUTINE shuffleK


   SUBROUTINE shuffleKBig(k,n_iterations)
      IMPLICIT NONE
      CLASS(KWaVe), TARGET :: k
      REAL(KIND=8) :: foo(0:k%num_dim)
      INTEGER, INTENT(IN) :: n_iterations
      INTEGER :: i1, i2
      REAL :: eta
      
      !check that k_bigger has been allocated
      IF (.NOT. ALLOCATED(k%k_big)) THEN
         PRINT *, "Error called by: quantum_momenta > shuffleKBig()"
         PRINT *, "The k-vectors k_big do not exist. You first need to build them."
         STOP
      END IF
		
      DO i2 = 1, n_iterations, 1
         i1=0
         DO WHILE (i1<k%num_big-1)
            i1=MAX(i1,0)
            i1=i1+1
            IF (k%k_big(0,i1)==k%k_big(0,i1+1)) THEN
               CALL RANDOM_NUMBER(eta)
               IF (eta>0.5) THEN
                   foo(0:k%num_dim)=k%k_big(0:k%num_dim,i1)
                   k%k_big(0:k%num_dim,i1)=k%k_big(0:k%num_dim,i1+1)
                   k%k_big(0:k%num_dim,i1+1)=foo(0:k%num_dim)
               END IF
            END IF
         END DO
      END DO

   END SUBROUTINE shuffleKBig

	
   SUBROUTINE copyFromKBig(k)
      IMPLICIT NONE
      CLASS(KWaVe), TARGET :: k
      
      !check that the target object k has been initialized
      IF (.NOT. k%flag_init) THEN
         PRINT *, "Error called by: quantum_momenta > buildBoxK()"
         PRINT *, "You first need to initialize the object."
         STOP
      END IF
      !check that k_big has been correctly allocated
      IF ((.NOT. ALLOCATED(k%k_big)) .OR. (k%num_big < k%num_part)) THEN
         PRINT *, "Error called by: quantum_momenta > buildBoxK()"
         PRINT *, "The k-vectors k_big do not exist or are not enough."
         STOP
      END IF
		
      k%k(0:k%num_dim,1:k%num_part)=k%k_big(0:k%num_dim,1:k%num_part)

   END SUBROUTINE copyFromKBig


   SUBROUTINE twistK(k)
      IMPLICIT NONE
      CLASS(KWaVe), TARGET :: k
      REAL(KIND=8) :: twist_angle(1:k%num_dim)
      REAL(KIND=8) :: backup_k_big(0:k%num_dim,1:k%num_big)
      INTEGER :: i1

      IF ( (.NOT. ALLOCATED(k%k_big)).OR.(.NOT. ALLOCATED(k%L)) ) THEN
         PRINT *, "Error called by: quantum_momenta > twistK()"
         PRINT *, "Could not find the larger grid and/or the L of the ",&
                  "box, probably the method buildBoxK() was not called before"
         STOP
      END IF
      
      !randomly generate the twist angle
      CALL RANDOM_NUMBER(twist_angle)
      twist_angle=(twist_angle-0.5d0)*2.d0*PI/k%L
      !make a copy of the k%k_big K-vector in backup_k_big
      backup_k_big=k%k_big
      !translate the K-vector k%k_big by the twist angle
      DO i1 = 1, k%num_big, 1
         k%k_big(1:k%num_dim,i1)=k%k_big(1:k%num_dim,i1)+twist_angle(1:k%num_dim)
         k%k_big(0,i1)=DOT_PRODUCT(k%k_big(1:k%num_dim,i1),k%k_big(1:k%num_dim,i1))
      END DO
      !shuffle the K-vector k%k_big
      CALL k%shuffleKBig(10)
      !reorder the K-vector k%k_big
      CALL k%putInOrderKBig()
      !Finally copy the first elements of k%k_big into k%k
      CALL k%copyFromKBig()

      k%k_big=backup_k_big

   END SUBROUTINE twistK


! - - - DESTRUCTOR - - - !

	SUBROUTINE deallocateKWaVe(k)
		!DESTRUCTOR
		IMPLICIT NONE
		CLASS(KWaVe), TARGET :: k
		
      !check that the target object k has been initialized
      IF (.NOT. k%flag_init) THEN
         PRINT *, "Error called by: quantum_momenta > buildBoxK()"
         PRINT *, "You cannot deallocate a non-initialized object."
         STOP
      END IF
		
		IF (ALLOCATED(k%k)) DEALLOCATE(k%k)
      IF (ALLOCATED(k%k_big)) DEALLOCATE(k%k_big) 
      k%flag_init=.FALSE.

	END SUBROUTINE deallocateKWaVe


! - - - SUBROUTINE USED INTERNALLY AND NOT VISIBLE FROM OUTSIDE - - - !

	FUNCTION countBoxK(k2,ndim,L)
		!Counts how many k are possible within a certain k^2
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: ndim
      REAL(KIND=8), INTENT(IN) :: k2, L(1:ndim)
		INTEGER :: countBoxK
		LOGICAL :: flag
		INTEGER :: i1, i2
      INTEGER :: nnet(1:ndim), iv(1:ndim)
      REAL(KIND=8) :: knet(1:ndim), v(0:ndim)
      REAL(KIND=8) :: kmin
		
      DO i2 = 1, ndim, 1
         kmin=2.d0*PI/L(i2)
         i1=FLOOR(DSQRT(k2)/kmin)
         nnet(i2)=i1
         knet(i2)=kmin
      END DO

		countBoxK=0
      iv(1:ndim)=nnet(1:ndim)
		v(1:ndim)=nnet(1:ndim)*knet(1:ndim)
		flag=.TRUE.
		DO WHILE ( flag )
			v(0)=DOT_PRODUCT(v(1:ndim),v(1:ndim))
			IF ( v(0)<=k2 ) THEN
				countBoxK=countBoxK+1
			END IF
         iv(ndim)=iv(ndim)-1
			v(ndim)=iv(ndim)*knet(ndim)
			DO i1 = ndim, 1, -1
				IF ( v(i1) < -k2 ) THEN
					IF ( i1==1 ) THEN
					 	flag=.FALSE.
					ELSE
                  iv(i1)=nnet(i1)
						v(i1)=iv(i1)*knet(i1)
                  iv(i1-1)=iv(i1-1)-1
						v(i1-1)=iv(i1-1)*knet(i1-1)
					END IF
				END IF
			END DO
		END DO
		
	END FUNCTION countBoxK


   SUBROUTINE provideBoxK(k2,ndim,L,num_k,k_box)
		!provides the first num_k k-vectors in the box with sizes L
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: ndim, num_k
      REAL(KIND=8), INTENT(IN) :: k2, L(1:ndim)
		LOGICAL :: flag
		INTEGER :: i1, i2
      INTEGER :: nnet(1:ndim), iv(1:ndim)
      REAL(KIND=8) :: knet(1:ndim), v(0:ndim)
      REAL(KIND=8) :: kmin
      REAL(KIND=8), INTENT(OUT) :: k_box(0:ndim,1:num_k)
		
      DO i2 = 1, ndim, 1
         kmin=2.d0*PI/L(i2)
         i1=FLOOR(DSQRT(k2)/kmin)
         nnet(i2)=i1
         knet(i2)=kmin
      END DO

		i2=0
      iv(1:ndim)=nnet(1:ndim)
		v(1:ndim)=nnet(1:ndim)*knet(1:ndim)
		flag=.TRUE.
		DO WHILE ( flag )
			v(0)=DOT_PRODUCT(v(1:ndim),v(1:ndim))
			IF ( v(0)<=k2 ) THEN
				i2=i2+1
            k_box(0:ndim,i2)=v(0:ndim)
            IF (i2>num_k) flag=.FALSE.
			END IF
         iv(ndim)=iv(ndim)-1
			v(ndim)=iv(ndim)*knet(ndim)
			DO i1 = ndim, 1, -1
				IF ( v(i1) < -k2 ) THEN
					IF ( i1==1 ) THEN
					 	flag=.FALSE.
					ELSE
                  iv(i1)=nnet(i1)
						v(i1)=iv(i1)*knet(i1)
                  iv(i1-1)=iv(i1-1)-1
						v(i1-1)=iv(i1-1)*knet(i1-1)
					END IF
				END IF
			END DO
		END DO
		
	END SUBROUTINE provideBoxK
	
	
END MODULE fermi_k
