cc ---------------------------------------------------------
cc Code realeased under Creative Common Licence
cc The code was used to generate the figures available in the Paper:
cc "Oscillators and relaxation phenomena in Pleistocene climate theory "
cc published in the Philosophical Transactions of the Royal Society
cc Code delivered without Garantee and Support whatsoever
cc Contact michel.crucifix for details 
cc Published on arXiv.org on 21 DEC 2012.
cc Michel Crucifix
cc ---------------------------------------------------------
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   vdp : the Van der pol oscillator in the slow-fast regime
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)


      DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897D0
      DOUBLE PRECISION, PARAMETER :: TPI = 2.D0 * 3.1415926535897D0

      DOUBLE PRECISION  PHI, DPHI, ALPHA, BET
      DOUBLE PRECISION, PARAMETER ::   ALPHA1=1.673915E-2
      DOUBLE PRECISION, PARAMETER ::   ALPHA2=9.523810E-3
      DOUBLE PRECISION, PARAMETER ::   ALPHA3=1.E-4
      DOUBLE PRECISION, PARAMETER ::   BETA1=5.118377E-1
      DOUBLE PRECISION, PARAMETER ::   BETA2=6.258680E-3
      DOUBLE PRECISION, PARAMETER ::   BETA3=2.639456E-5
      DOUBLE PRECISION, PARAMETER ::   BETA4=3.628118E-8
      DOUBLE PRECISION, PARAMETER ::   BETA5=5.833333E-3
      DOUBLE PRECISION, PARAMETER ::   GAMMA1=1.851250E-3
      DOUBLE PRECISION, PARAMETER ::   GAMMA2=1.125000E-5
      DOUBLE PRECISION, PARAMETER ::   GAMMA3=2.5E-4
      DOUBLE PRECISION, PARAMETER ::   CMU   =4.0E-3
      DOUBLE PRECISION, PARAMETER ::   KTHETA=4.4444444E-2

      DOUBLE PRECISION I, MU, THETA
      BET   = PAR(1)
      
      I = U(1)
      MU=U(2)
      THETA=U(3) 

      F(1) =  ALPHA1 - ALPHA2 * CMU * MU - ALPHA3*I - ALPHA2*KTHETA*  THETA 
      F(2) =  BETA1 - (BETA2 - BETA3*MU + BETA4*MU*MU)*MU - BETA5  * THETA + BET
      F(3) =  GAMMA1 - GAMMA2 * I - GAMMA3 * THETA

       IF(IJAC.EQ.0)RETURN

!         DFDU(1,1) = - 0.8
!          DFDU(1,2) = -1

!        DFDU(2,1) = +  1
!        DFDU(2,2) = - 3*U(2)*U(2) - 1.2*U(2) + 1.3


      IF(IJAC.EQ.1)RETURN

!         DFDP(1,1) = 0
!         DFDP(2,1) = 1


      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897D0
      DOUBLE PRECISION   PHI
         PAR(1) = 0.02
         U(1) = 14.477431
         U(2) = 326.3605
         U(3) = 6.753516

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
