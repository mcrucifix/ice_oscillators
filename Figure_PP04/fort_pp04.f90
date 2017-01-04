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

      double precision ALPHA, BETA, GAMMA, OMEGA
      double precision TAUV, TAUC, TAUA, XX, YY, ZZ
      double precision GGAMMA, DDELTA, AA, BB, CC, D
      double precision A, V, C
      double precision  dX, H, FF
      integer          nw,ix,NN
 

      PARAMETER  (TAUV = 1.5)
!      PARAMETER  (TAUC = .5)
!      PARAMETER  (TAUA = 1.2)
      PARAMETER  (XX   = 1.3)
      PARAMETER  (YY   = 0.5)
!      PARAMETER  (ZZ   = 0.8)
      PARAMETER  (ALPHA = 0.15)
      PARAMETER  (BETA  = 0.5)
!      PARAMETER  (GGAMMA = 0.5)
!      PARAMETER  (DDELTA = 0.4)
      PARAMETER  (AA = 0.3)
!      PARAMETER  (BB = 0.7)
      PARAMETER  (CC = 0.00)
!      PARAMETER  (D = 0.27)


        D  = par(1)
        DDELTA  = par(2)
        zz  = par(3)
        bb  = par(4)
        GGAMMA  = par(5)
        TAUA  = par(6)
        TAUC  = par(7)

        FF = aa*U(1) - bb*U(2) + d
        H = (ATAN(-500.*FF)/3.141592+0.5)
!        if (FF > 0) H = 0.


        F(1) = (- xx * U(3)   + zz - U(1))/TAUV
        F(2) = (U(1) - U(2) ) / TAUA
        F(3) = ( - BETA*U(1) - U(3) + GGAMMA * H + DDELTA ) / TAUC

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897D0
      DOUBLE PRECISION   PHI
         PAR(1) = 0.
         PAR(2) = 0.
         PAR(3) = 0.
         PAR(4) = 0.
         PAR(5) = 0.
         PAR(6) = 1.2
         PAR(7) = 0.5
         U(1) = 0.
         U(2) = 0.
         U(3) = 0.

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
