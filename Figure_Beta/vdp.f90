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

      BET   = PAR(1)
      ALPHA = PAR(2)

        F(1) =  +  ALPHA * (U(2) -  PHI(U(1)))
        F(2) =  -  ( U(1) + BET  ) 

       IF(IJAC.EQ.0)RETURN

        DFDU(1,1) = -   DPHI(U(1)) *ALPHA
        DFDU(1,2) = + ALPHA

        DFDU(2,1) = -  1
        DFDU(2,2) = 0.


      IF(IJAC.EQ.1)RETURN

         DFDP(1,1) = 0.
         DFDP(1,2) =F(1) / ALPHA
         DFDP(2,1) =  - 1
         DFDP(2,2) =  0.


      END SUBROUTINE FUNC

      FUNCTION PHI(X)
      DOUBLE PRECISION PHI, X 
      PHI = (X*X*X/3.D0 - X)
      RETURN
      END FUNCTION

      FUNCTION DPHI(X)
      DOUBLE PRECISION DPHI, X 
      DPHI = (X*X - 1)
      RETURN
      END FUNCTION
      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897D0
      DOUBLE PRECISION   PHI
         PAR(1) = 1.5
         PAR(2)=1.
         U(1) = -PAR(1)
         U(2) = PHI(U(1))

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
