cc ---------------------------------------------------------
cc Code realeased under Creative Common Licence
cc The code was used to generate the figures available in the Paper:
cc "Oscillators and relaxation phenomena in Pleistocene climate theory "
cc published in the Philosophical Transactions of the Royal Society
cc Code delivered without Garantee and Support whatsoever
cc Contact michel.crucifix for details 
cc Published on arXiv.org on 21 DEC 2012.
cc Michel Crucifix
cc I keep intellectual ownership of the work
cc ---------------------------------------------------------
cc ---------------------------------------------------------
cc Code realeased under Creative Common Licence
cc The code was used to generate the figures available in the Paper:
cc "Oscillators and relaxation phenomena in Pleistocene climate theory "
cc published in the Philosophical Transactions of the Royal Society
cc Code delivered without Garantee and Support whatsoever
cc Contact michel.crucifix for details 
cc Published on arXiv.org on 21 DEC 2012.
cc Michel Crucifix
cc I keep intellectual ownership of the work
cc ---------------------------------------------------------
!     generates a file of initial conditions for lor.f90
!        explore empirically equilibrium points
!        generate initial conditions in auto. 
!        PAR(2) is the forcing amplitude. 

      PROGRAM ITER

      IMPLICIT NONE
      INTEGER, PARAMETER :: NDIM = 5
      INTEGER, PARAMETER :: NPAR = 3
      DOUBLE PRECISION  U(NDIM),PAR(11)
      DOUBLE PRECISION  T
      DOUBLE PRECISION  F(NDIM)
      DOUBLE PRECISION  DFDU(NDIM,NDIM), DFDP(NDIM,NPAR)
      DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897D0
      DOUBLE PRECISION, PARAMETER :: TPI = 2 * PI
      DOUBLE PRECISION  DELTAT
      INTEGER NINT
      INTEGER ICP(1), IJAC
      INTEGER I,J,K,L,M

      IJAC=0
      CALL STPNT(NDIM,U,PAR,T)
      NINT = 600000
      DELTAT = PAR(11) / NINT
      DO I=1,50
         DO J=1,NINT 
              CALL FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
              U(:) = U(:) + F(:) * DELTAT
         ENDDO
      ENDDO


      DO J=0,NINT 
           CALL FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
           U(:) = U(:) + F(:) * DELTAT
           IF (MOD(J,1000).EQ.0) WRITE(*,'(6e15.7)') J*DELTAT,U(:)
      ENDDO
      END PROGRAM
   



