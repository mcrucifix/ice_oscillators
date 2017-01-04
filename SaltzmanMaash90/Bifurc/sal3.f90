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
      DOUBLE PRECISION, PARAMETER ::   ALPHA1=1.807500E-2*1E3
      DOUBLE PRECISION, PARAMETER ::   ALPHA2=1.275000E-2*1E3
      DOUBLE PRECISION, PARAMETER ::   ALPHA3=1.E-4*1E3
      DOUBLE PRECISION, PARAMETER ::   BETA1=1.35560800*1E3
      DOUBLE PRECISION, PARAMETER ::   BETA2=0.0054688000*1E3
      DOUBLE PRECISION, PARAMETER ::   BETA3=0.0022130000*1E3
      DOUBLE PRECISION, PARAMETER ::   BETA4=0.0002200000*1E3
      DOUBLE PRECISION, PARAMETER ::   BETA5=0.5410550000*1E3
      DOUBLE PRECISION, PARAMETER ::   BETA6=0.0532000000*1E3
      DOUBLE PRECISION, PARAMETER ::   GAMMA1=0.0018360000*1E3
      DOUBLE PRECISION, PARAMETER ::   GAMMA2=0.0000120000*1E3
      DOUBLE PRECISION, PARAMETER ::   GAMMA3=0.0002400000*1E3
      DOUBLE PRECISION, PARAMETER ::   CMU   =4.E-3
      DOUBLE PRECISION, PARAMETER ::   KTHETA=3.333E-2

      DOUBLE PRECISION I, MU, THETA
      BET   = PAR(1)
      
      I = U(1)
      MU=U(2)
      THETA=U(3) 

      F(1) =  ALPHA1 - ALPHA2 * CMU * MU - ALPHA3*I - ALPHA2*KTHETA*  THETA 
      F(2) =  BETA1 - (BETA2 - BETA3*THETA + BETA4*THETA*THETA)*MU - (BETA5 - BETA6 * THETA) * THETA + BET
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
         PAR(1) = 0.06*1E3
         U(1) = 4.40479852
         U(2) = 283.8602
         U(3) = 7.429760

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
