!$Id: lattice.f90,v 1.3 2006/06/29 11:59:12 schmidt Exp $ modified Calcavecchia 16/6/2014 (added graphene layer)
module lattice
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

contains

   function simplecubic(nside)												! crea un reticolo cubico semplice con nside siti per lato della cella   
!																			! cubica unitaria
! set up a simple cubic lattice with nside unit cells per side in
! the unit cube
!  
   implicit none
   integer(kind=i4) :: nside
   real(kind=r8), dimension(3,nside**3) :: simplecubic
   real(kind=r8), dimension(3,1) :: basis = &
      reshape((/0.25_r8, 0.25_r8, 0.25_r8/),(/3,1/))						! crea la base simple cubic 3x1 (1 punto nella cella elementare)
   simplecubic=basiscubic(nside,basis)
   return
   end function simplecubic

   function bodycenteredcubic(nside)										! crea un reticolo cubico a corpo centrato con nside siti per lato della cella
!                                                                           ! cubica unitaria
! set up a body-centered cubic lattice with nside unit cells per side in
! the unit cube
!
   implicit none
   integer(kind=i4) :: nside
   real(kind=r8), dimension(3,2*nside**3) :: bodycenteredcubic
   real(kind=r8), dimension(3,2) :: basis = &
      reshape((/0.25_r8, 0.25_r8, 0.25_r8, 0.75_r8, 0.75_r8, 0.75_r8/) &	! crea la base body centered cubic 3x2 (2 punti nella cella elementare)
      ,(/3,2/))
   bodycenteredcubic=basiscubic(nside,basis)
   return
   end function bodycenteredcubic
   
   function facecenteredcubic(nside)										! crea un reticolo cubico a facce centrate con nside siti per lato della cella
!                                                                           ! cubica unitaria
! set up a face-centered cubic lattice with nside unit cells per side in
! the unit cube
!
   implicit none
   integer(kind=i4) :: nside
   real(kind=r8), dimension(3,4*nside**3) :: facecenteredcubic
   real(kind=r8), dimension(3,4) :: basis = &
      reshape((/0.25_r8, 0.25_r8, 0.25_r8, 0.75_r8, 0.75_r8, 0.25_r8, &
       0.75_r8, 0.25_r8, 0.75_r8, 0.25_r8, 0.75_r8, 0.75_r8/),(/3,4/))		! crea la base face centered cubic 3x4 (4 punti nella cella elementare)
   facecenteredcubic=basiscubic(nside,basis)
   return
   end function facecenteredcubic

	FUNCTION hexagonalclosepacking(nside)                                               ! crea un reticolo hcp con nside siti per lato della cella. Il diamentro
																						! delle particelle é 1/nside. una cella ha dimensione: 2 x sqrt(3) x 2 sqrt(2/3).
																						! Le particelle saranno 8*(nside**3)
		IMPLICIT NONE
		INTEGER :: nside
		REAL(KIND=8), DIMENSION(3,8*nside**3) :: hexagonalclosepacking
		REAL (KIND=8), DIMENSION(3,8) :: basis
		REAL (KIND=8) :: vbasis(24)
		REAL (KIND=8), PARAMETER :: const(1:3)=(/0.25d0,DSQRT(3.d0)/6.d0,1.d0/DSQRT(6.d0)/)
		vbasis=(/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 1.5d0,0.5d0*DSQRT(3.d0),0.d0, 0.5d0,0.5d0*DSQRT(3.d0),0.d0, &
		         0.5d0,DSQRT(3.d0)/6.d0,DSQRT(2.d0/3.d0), 1.5d0,DSQRT(3.d0)/6.d0,DSQRT(2.d0/3.d0), &
		         1.d0,2.d0*DSQRT(3.d0)/3.d0,DSQRT(2.d0/3.d0), 0.d0,2.d0*DSQRT(3.d0)/3.d0,DSQRT(2.d0/3.d0) /)
		basis=RESHAPE(vbasis,(/3,8/))
		basis(1,1:8)=0.5d0*basis(1,1:8)
		basis(2,1:8)=basis(2,1:8)/DSQRT(3.d0)
		basis(3,1:8)=0.5d0*basis(3,1:8)/(DSQRT(2.d0/3.d0))
		hexagonalclosepacking=basiscubic(nside,basis)
		hexagonalclosepacking(1,1:8*nside**3)=2.d0*hexagonalclosepacking(1,1:8*nside**3)+const(1)
		hexagonalclosepacking(2,1:8*nside**3)=hexagonalclosepacking(2,1:8*nside**3)*DSQRT(3.d0)+const(2)
		hexagonalclosepacking(3,1:8*nside**3)=2.d0*hexagonalclosepacking(3,1:8*nside**3)*(DSQRT(2.d0/3.d0))+const(3)
	END FUNCTION hexagonalclosepacking                                                   !autore: Francesco Calcavecchia

   function basiscubic(nside,basis)											! data una base, crea un reticolo cubico con con nside siti per lato della cella
!                                                                           ! cubica unitaria
! given a basis in the cell, set up a cubic lattice with nside unit cells
! per side in the unit cube
!
   implicit none
   integer(kind=i4) :: nside												! nside: siti per lato della cella cubica unitaria
   real(kind=r8), dimension(:,:) :: basis
   real(kind=r8), dimension(3,nside**3*size(basis)/3) :: basiscubic
   real(kind=r8) :: scale
   integer(kind=i4) :: i,i1,i2,i3,ib,nbasis
   nbasis=size(basis)/3														! nbasis: punti della base
   scale=1.0_r8/nside														! passo reticolare
   i=0
   do i1=0,nside-1															! ciclo lungo x
      do i2=0,nside-1														! ciclo lungo y	
         do i3=0,nside-1													! ciclo lungo z		
            do ib=1,nbasis													! ciclo sui punti della base
               i=i+1
               basiscubic(1,i)=(i1+basis(1,ib))*scale						! reticolo x
               basiscubic(2,i)=(i2+basis(2,ib))*scale						! reticolo y
               basiscubic(3,i)=(i3+basis(3,ib))*scale						! reticolo z
            enddo
         enddo
      enddo
   enddo
   return
   end function basiscubic

   function bestcubic(npart)												! dato il numero di siti npart, costruisce il miglior reticolo cubico, quello che
!																			! minimizza gli spazi vuoti nella struttura reticolare
! try simple, body-centered, and face-centered cubic and return npart
! positions in the unit cube that has the fewest vacancies for npart
! particles
!
   implicit none
   integer(kind=i4) :: npart
   real(kind=r8), dimension(3,npart) :: bestcubic
   integer(kind=i4) :: nsc,nbcc,nfcc,nsidesc,nsidebcc,nsidefcc
   real(kind=r8), dimension(:,:), allocatable :: x
   real(kind=r8) :: en
   en=npart
   nsidesc=nint(en**(1.0_r8/3.0_r8))										! nsidesc:  siti per lato sc	=>	npart^(1/3)
   if (nsidesc**3.lt.npart) nsidesc=nsidesc+1								! 			nsidesc^3 deve essere maggiore o ugaule a npart	
   nsc=nsidesc**3															! nsc: particelle sc
   nsidebcc=nint((en*0.5_r8)**(1.0_r8/3.0_r8))								! nsidebcc: siti per lato bcc	=>	(1/2*npart)^(1/3)
   if (nsidebcc**3*2.lt.npart) nsidebcc=nsidebcc+1                          ! 			2*nsidebcc^3 deve essere maggiore o ugaule a npart
   nbcc=nsidebcc**3*2                                                       ! nbcc: particelle bcc
   nsidefcc=nint((en*0.25_r8)**(1.0_r8/3.0_r8))								! nsidefcc: siti per lato fcc	=>	(1/4*npart)^(1/3)
   if (nsidefcc**3*4.lt.npart) nsidefcc=nsidefcc+1                          ! 			4*nsidefcc^3 deve essere maggiore o ugaule a npart
   nfcc=nsidefcc**3*4                                                       ! nfcc: particelle fcc
   if (nsc.le.nbcc.and.nsc.le.nfcc) then									! numero particelle sc minore 	=>	più vicino a npart
      allocate(x(3,nsc))													!								=>	scelta: simple cubic
      x = simplecubic(nsidesc)												
      else if (nbcc.le.nfcc) then											! numero particelle bcc minore	=>	più vicino a npart
      allocate(x(3,nbcc))													!								=>	scelta: body centered cubic
      x = bodycenteredcubic(nsidebcc)
      else   																! numero particelle fcc minore	=>	più vicino a npart
      allocate(x(3,nfcc))													!								=>	scelta: face centered cubic
      x = facecenteredcubic(nsidefcc)
      endif
   bestcubic(:,:)=x(:,1:npart:1)											! registrazione del cubo migliore fino a npart
   deallocate(x)
   return
   end function bestcubic


	FUNCTION graphene_layer(nside)
		INTEGER, INTENT(IN) :: nside
		REAL (KIND=8) :: graphene_layer(1:3,4*(nside**(2)))
		INTEGER :: nx, ny, cont
		REAL (KIND=8) :: lxy
		
		lxy=1.d0/REAL(nside,8)
		
		cont=1
		!spinup
		DO nx = 0, nside-1, 1
			DO ny = 0, nside-1, 1
				graphene_layer(1:2,cont)=(/0.d0,0.5d0/)*lxy+nx*(/lxy,0.d0/)+ny*(/0.d0,lxy/)
				graphene_layer(3,cont)=0.d0
				cont=cont+1
				graphene_layer(1:2,cont)=(/0.5d0,0.d0/)*lxy+nx*(/lxy,0.d0/)+ny*(/0.d0,lxy/)
				graphene_layer(3,cont)=0.d0
				cont=cont+1
			END DO
		END DO
		!spindown
		DO nx = 0, nside-1, 1
			DO ny = 0, nside-1, 1
				graphene_layer(1:2,cont)=(/1.d0/6.d0,0.d0/)*lxy+nx*(/lxy,0.d0/)+ny*(/0.d0,lxy/)
				graphene_layer(3,cont)=0.d0
				cont=cont+1
				graphene_layer(1:2,cont)=(/2.d0/3.d0,0.5d0/)*lxy+nx*(/lxy,0.d0/)+ny*(/0.d0,lxy/)
				graphene_layer(3,cont)=0.d0
				cont=cont+1
			END DO
		END DO
		
	END FUNCTION graphene_layer

end module lattice
