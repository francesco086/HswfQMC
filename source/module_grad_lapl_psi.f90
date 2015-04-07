!!!NOT FINISHED TO IMPLEMENT YET, THERE 
!!!ARE SOME TRACKS IN module_estimatori.f90


MODULE grad_lapl_psi
   USE funzione_onda
   USE walkers
   USE dati_fisici
   USE dati_mc
   IMPLICIT NONE
   REAL(KIND=8), PROTECTED, SAVE, ALLOCATABLE :: grad_Jee(:,:,:)
   REAL(KIND=8), PROTECTED, SAVE, ALLOCATABLE :: lapl_Jee(:,:,:)
   REAL(KIND=8), PROTECTED, SAVE, ALLOCATABLE :: grad_Jep(:,:,:)
   REAL(KIND=8), PROTECTED, SAVE, ALLOCATABLE :: lapl_Jep(:,:,:)

CONTAINS

   SUBROUTINE inizializza_grad_lapl_psi()
      IMPLICIT NONE
      
      IF (Jee_kind/="no_") THEN
         ALLOCATE(grad_Jee(1:3,1:N_part,0:N_mc))
         ALLOCATE(lapl_Jee(1:3,1:N_part,0:N_mc))
      END IF
      IF (Jep_kind/="no_") THEN
         ALLOCATE(grad_Jep(1:3,1:N_part,0:N_mc))
         ALLOCATE(lapl_Jep(1:3,1:N_part,0:N_mc))
      END IF
      
   END SUBROUTINE inizializza_grad_lapl_psi
!###################################################################### 
   SUBROUTINE calcola_grad_lapl_psi(i_mc)
      IMPLICIT NONE
      INTEGER(KIND=8) :: i_mc

      !Jee
      CALL calcola_grad_lapl_Jee(i_mc)
      !Jep
      CALL calcola_grad_lapl_Jep(i_mc)
      
   END SUBROUTINE calcola_grad_lapl_psi
!###################################################################### 
   SUBROUTINE conferma_grad_lapl_psi(i_mc)
      IMPLICIT NONE
      INTEGER(KIND=8) :: i_mc

      !Jee
      IF (Jee_kind/="no_") THEN
         grad_Jee(1:3,1:N_part,i_mc)=grad_Jee(1:3,1:N_part,i_mc-1)
         lapl_Jee(1:3,1:N_part,i_mc)=lapl_Jee(1:3,1:N_part,i_mc-1)
      END IF
      !Jep
      IF (Jep_kind/="no_") THEN
         grad_Jep(1:3,1:N_part,i_mc)=grad_Jep(1:3,1:N_part,i_mc-1)
         lapl_Jep(1:3,1:N_part,i_mc)=lapl_Jep(1:3,1:N_part,i_mc-1)
      END IF
      
   END SUBROUTINE conferma_grad_lapl_psi
!###################################################################### 
   SUBROUTINE calcola_grad_lapl_Jee(i_mc)
      IMPLICIT NONE
      INTEGER(KIND=8) :: i_mc
      
      SELECT CASE(Jee_kind)
      CASE ("yuk")
         CALL grad_lapl_Jee_yuk(grad_Jee(:,:,i_mc),lapl_Jee(:,:,i_mc))
      END SELECT

      IF (Jee_kind/="no_") THEN
         grad_Jee(1:3,1:N_part,0)=grad_Jee(1:3,1:N_part,i_mc)
         lapl_Jee(1:3,1:N_part,0)=lapl_Jee(1:3,1:N_part,i_mc)
      END IF
      
   END SUBROUTINE calcola_grad_lapl_Jee
!###################################################################### 
   SUBROUTINE calcola_grad_lapl_Jep(i_mc)
      IMPLICIT NONE
      INTEGER(KIND=8) :: i_mc
      
      SELECT CASE(Jep_kind)
      CASE ("yuk")
         CALL grad_lapl_Jep_yuk(grad_Jep(:,:,i_mc),lapl_Jep(:,:,i_mc))
      END SELECT

      IF (Jep_kind/="no_") THEN
         grad_Jep(1:3,1:N_part,0)=grad_Jep(1:3,1:N_part,i_mc)
         lapl_Jep(1:3,1:N_part,0)=lapl_Jep(1:3,1:N_part,i_mc)
      END IF
      
   END SUBROUTINE calcola_grad_lapl_Jep
!###################################################################### 
   SUBROUTINE grad_lapl_generic_J(within_same_particles,u1,u2,rij,grad,lapl)
   !calcola il gradiente e laplaciano di un jastrow generico di forma
   !J(R) = exp( -u(R) )
   !In input bisogna fornire la derivata prima e seconda di u
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: within_same_particles
      REAL(KIND=8), INTENT(IN) :: u1(1:N_part,1:N_part),u2(1:N_part,1:N_part)
      REAL(KIND=8), INTENT(IN) :: rij(0:3,1:N_part,1:N_part)
      REAL(KIND=8), INTENT(OUT) :: grad(1:3,1:N_part),lapl(1:3,1:N_part)
      INTEGER :: i1,i2,li1,ui1,li2,ui2
      REAL(KIND=8) :: du(1:3,1:2),d2u(1:3,1:2)
      REAL(KIND=8) :: frf(1:3), frf0

      IF (within_same_particles) THEN
         li1=1; ui1=N_part-1; li2=N_part; ui2=N_part
      ELSE
         li1=1; ui1=N_part; li2=1; ui2=N_part
      END IF

      !NOTA: tenendo conto che r(i,j)=r(j)-r(i), bisogna tenere conto di
      !un segno meno in frf
      grad=0.d0; lapl=0.d0
      DO i1 = li1, ui1, 1
      DO i2 = MIN(li2,i1+1), ui2, 1
         !derivate prime du(R)/dr_alpha_i
         frf(1:3) = -rij(1:3,i2,i1)*u1(i2,i1)/rij(0,i2,i1)
         du(1:3,2) = frf(1:3)
         du(1:3,1) = -frf(1:3)
         !derivate seconde d2u(R)/dr2_alpha_i
         frf0 = u1(i2,i1)/rij(0,i2,i1)
         frf(1:3) = ((rij(1:3,i2,i1)**2)/rij(0,i2,i1)) &
                      * ( - u1(i2,i1)/(rij(0,i2,i1)**2) &
                          + u2(i2,i1)/rij(0,i2,i1) )
         d2u(1:3,2) = (frf0+frf(1:3))
         d2u(1:3,1) = (frf0+frf(1:3))
         !gradiente di J diviso per J
         grad(1:3,i2) = grad(1:3,i2) - du(1:3,2)
         IF (within_same_particles) grad(1:3,i1) = grad(1:3,i1) - du(1:3,1)
         !laplaciano di J diviso per J
         lapl(1:3,i2) = lapl(1:3,i2) + du(1:3,2)**2 - d2u(1:3,2)
         IF (within_same_particles) lapl(1:3,i1) = lapl(1:3,i1) + du(1:3,1)**2 - d2u(1:3,1)
      END DO
      END DO
      
   END SUBROUTINE grad_lapl_generic_J
!###################################################################### 
   SUBROUTINE grad_lapl_Jee_yuk(grad,lapl)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(OUT) :: grad(1:3,1:N_part), lapl(1:3,1:N_part)
      INTEGER :: i1,i2
      REAL(KIND=8) :: A,F,r
      REAL(KIND=8) :: u1(1:N_part,1:N_part),u2(1:N_part,1:N_part)

      !calcolo le derivate dello pseudopotenziale
      DO i1 = 1, N_part-1, 1
      DO i2 = i1+1, N_part, 1
         !setto i termini A e F
         A=Aee_yuk
         IF ( (split_Aee) .AND. &
              ( ((i2<=H_N_part) .AND. (i1>H_N_part)) .OR. &
                ((i2>H_N_part) .AND. (i1<=H_N_part)) ) ) A=Aee_ud_yuk
         F=Fee_yuk
         IF ( (split_Fee) .AND. &
              ( ((i2<=H_N_part) .AND. (i1>H_N_part)) .OR. &
                ((i2>H_N_part) .AND. (i1<=H_N_part)) ) ) F=Fee_ud_yuk
         !derivata prima e seconda di u
         r=rij_ee_old(0,i2,i1)
         u1(i2,i1) = A*F*DEXP(-F*r)/r &
              - A*(1.d0-DEXP(-F*r))/(r**2)
         u2(i2,i1) = - A*F*F*DEXP(-F*r)/r &
              - 2.d0*A*F*DEXP(-F*r)/(r**2) &
              + 2.d0*A*(1.d0-DEXP(-F*r))/(r**3) 
         !u e' diviso per 2 nel Jastrow
         u1(i2,i1)=u1(i2,i1)*0.5d0; u2(i2,i1)=u2(i2,i1)*0.5d0
      END DO
      END DO

      CALL grad_lapl_generic_J(.TRUE.,u1,u2,rij_ee_old,grad,lapl)

   END SUBROUTINE grad_lapl_Jee_yuk
!###################################################################### 
   SUBROUTINE grad_lapl_Jep_yuk(grad,lapl)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(OUT) :: grad(1:3,1:N_part), lapl(1:3,1:N_part)
      INTEGER :: i1,i2
      REAL(KIND=8) :: A,F,r
      REAL(KIND=8) :: u1(1:N_part,1:N_part),u2(1:N_part,1:N_part)

      !calcolo le derivate dello pseudopotenziale
      DO i1 = 1, N_part, 1
      DO i2 = 1, N_part, 1
         !setto i termini A e F
         A=Aep_yuk
         IF ( (split_Aep) .AND. &
              ( ((i2<=H_N_part) .AND. (i1>H_N_part)) .OR. &
                ((i2>H_N_part) .AND. (i1<=H_N_part)) ) ) A=Aep_ud_yuk
         F=Fep_yuk
         IF ( (split_Fep) .AND. &
              ( ((i2<=H_N_part) .AND. (i1>H_N_part)) .OR. &
                ((i2>H_N_part) .AND. (i1<=H_N_part)) ) ) F=Fep_ud_yuk
         !derivata prima e seconda di u
         r=rij_ep_old(0,i2,i1)
         u1(i2,i1) = A*F*DEXP(-F*r)/r &
              - A*(1.d0-DEXP(-F*r))/(r**2)
         u2(i2,i1) = - A*F*F*DEXP(-F*r)/r &
              - 2.d0*A*F*DEXP(-F*r)/(r**2) &
              + 2.d0*A*(1.d0-DEXP(-F*r))/(r**3) 
         !u e' diviso per 2 nel Jastrow
         u1(i2,i1)=u1(i2,i1)*0.5d0; u2(i2,i1)=u2(i2,i1)*0.5d0
      END DO
      END DO

      CALL grad_lapl_generic_J(.FALSE.,u1,u2,rij_ep_old,grad,lapl)

   END SUBROUTINE grad_lapl_Jep_yuk
!###################################################################### 
   SUBROUTINE chiudi_grad_lapl_psi()
      IMPLICIT NONE
      
      IF (Jee_kind/="no_") THEN
         DEALLOCATE(grad_Jee)
         DEALLOCATE(lapl_Jee)
      END IF
      IF (Jep_kind/="no_") THEN
         DEALLOCATE(grad_Jep)
         DEALLOCATE(lapl_Jep)
      END IF
      
   END SUBROUTINE chiudi_grad_lapl_psi

END MODULE grad_lapl_psi
