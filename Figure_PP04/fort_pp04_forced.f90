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

      double precision ALPHA, BETA, OMEGA
      double precision TAUV, TAUC, TAUA, XX, YY, ZZ
      double precision GGAMMA, DDELTA, AA, BB, CC, D
      double precision R
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
        OMEGA = 2*pi/4.1 
        FF = aa*U(1) - bb*U(2) + d
        H = (ATAN(-500.*FF)/3.141592+0.5)
!        if (FF > 0) H = 0.

        F(1) = (- xx * U(3) + yy * U(4)  + zz - U(1))/TAUV
        F(2) = (U(1) - U(2) ) / TAUA
        F(3) = ( ALPHA * U(4)  - BETA*U(1) - U(3) + GGAMMA * H + DDELTA ) / TAUC
        R    = U(4)*U(4)+U(5)*U(5)
        F(4) = U(4) * (1.-R)  + OMEGA * U(5) 
        F(5) = U(5) * (1.-R)  - OMEGA * U(4)

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897D0
      DOUBLE PRECISION   D, DDELTA, zz, BB, GGAMMA, Q

         OPEN(10, file='vdp.startpar', status='old')

         READ(10,*) D
         READ(10,*) DDELTA
         READ(10,*) zz
         READ(10,*) BB
         READ(10,*) GGAMMA
         READ(10,*) Q



         PAR(1) = D
         PAR(2) = DDELTA
         PAR(3) = zz
         PAR(4) = BB
         PAR(5) = GGAMMA
         PAR(6) = 1.2
         PAR(7) = 0.5
         
         PAR(11) = 4.1 * Q

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
