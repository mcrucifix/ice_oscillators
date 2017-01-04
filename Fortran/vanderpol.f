! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!! ATTENTION 
!  THE ROUTINE IS ENTIRELY WRITTEN ASSUMING  TIME
!  AND AMPLITUDE SCALINGS SUCH THAT THE
!  DOMINANT TERM OF THE INSOLATION EXPANSION
!  HAS AN AMPLITUDE = 1 AND ANGULAR VELOCITY
!  = 1. 
!  IF ONLY ONE TERM IS RETAINED IN THE EXPANSIOT
!  OF INSOLATION [NA=1] , THIS IS EQUIVALENT
!  TO A FORCING OF THE FORM SIN(T). 
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


      subroutine read_insol(na,ampsin,ampcos,ome)

       implicit none
       integer na
       integer i    ! counter
       double precision ampsin(na),ampcos(na),ome(na)
       double precision tmp1,tmp2,tmp3 ! tmp
       double precision scale
       parameter (scale = 1d4 / 1.5324947854)

      open (10, file='../Data/insol65nber78.dat'   ,status='old')
       read (10,*)  ! skip first line
       do i=1,na
         read (10,*) tmp1,tmp2,tmp3
         ome(i)=tmp1*scale
         ampsin(i)=tmp2
         ampcos(i)=tmp3
       enddo
       close(10)
       end

      subroutine insol(na,t,ampsin,ampcos,ome,x,dxdt)
        implicit none
        integer  na
        integer  i                                     ! counter
        double precision t,ampsin(na),ampcos(na),ome(na) ! input
        double precision x,dxdt                        ! output
        double precision some(na),come(na)             ! tmp

        some = dsin(ome*t)
        come = dcos(ome*t)
!        come = dsqrt(1.d0-some*some)
        x=sum(ampsin*some + ampcos*come)
        dxdt = sum(ome*(ampsin*come-ampcos*some))
      end 

      subroutine propagate_lin (N, state, par, t0, t1, b, ix,  output, 
     c              icalclyap, ds, lyap, na) 
!     MAIN ROUTIN called f  rom R 
!     IN PARAMETERS 
!     -------------
!     N : number of parameters
!     par   (N,4)  : parameters : amplitude and natural omega
!     t0 , t1      : initial and final time
!     ix           : seed
!     b            : noise amplitude
!     na           : number of terms retained in coef. expansion
!     icalclyap    : set 1 if you want it to be calcuated. 

!     INOUT
!     -----
!     state (N,2)  : initial conditions
!     lyap         : greatest lyapunov coefficient. Set eq 1 
!     ds           : direction of the most unstable direction
!                    (if not known : may be obtained by spin-up)

       integer I,ix,l, IMAX 
       integer icalclyap
       integer output
       integer N   ! number of particles
       integer na  ! number of terms in the development of insolation
       double precision t0,t1, lyap(N)
       double precision par(N,4),state(N,2),dx,dY(N,2), dL(N,2)
       double precision ds(N,2), norm(N)
       double precision t,deltat,b(2)
       double precision f(2),dfdt(2)
       double precision ampsin(na),ampcos(na),ome(na)

       !! scaled such as first term has amplitude 1
       double precision insol_scale
       parameter (insol_scale = 1.d0/11.766587361d0)


       double precision x,dxdt           ! temporary variables for exch.
                                         ! with insol
       deltat  = 0.02

       if (output.EQ.1) open(20, file = 'trajectory.dat', 
     c                  status='unknown')

       !! adjusts imax to have an int. numb. of timesteps
       IMAX = max(int(((t1-t0)/deltat)+0.5),1)
       deltat  = (t1-t0)/IMAX     

       !! read insolation data
       call read_insol(na,ampsin,ampcos,ome)

       f = 0
       dfdt = 0
        
       DO I=1,(IMAX)
         t = t0+(I-1)*deltat
         call insol(na,t, ampsin,ampcos, ome, x, dxdt)
         f(1)    = x    * insol_scale
         dfdt(1) = dxdt * insol_scale
         call ddtlin_F(state,N,par,b,f,dfdt,deltat,dY,ix,ds,dl)
         state = state + dY

         if (icalclyap .EQ. 1) THEN
          ds = ds + dL
          norm = (sqrt(ds(:,1)*ds(:,1) +
     c    ds(:,2)*ds(:,2)))
          lyap = lyap + dlog(norm)
          ds(:,1) = ds(:,1) / norm
          ds(:,2) = ds(:,2) / norm
         endif
         if (output.EQ.1)  write(20,*) t,state, f(1)
       ENDDO
       if (output.EQ.1)  close(20)
      end 
 
      function P(x) ! phi(x)
       double precision P,x
        P=x*x*x/3.-x
      end function
 
      function Pp(x) ! phi(x)
       double precision Pp,x
        Pp=x*x-1.
      end function


      subroutine ddtlin_F(state,N,par,b,f,dfdt,deltat,dY,ix,ds,dl)
      implicit none
      integer K,J,l,M,nwork,il
      integer N
      parameter (nwork=1e4)
      double precision par(N,4),state(N,2), ds(N,2),dl(N,2)
      double precision X,Y
      double precision q,r,eps
      double precision deltat,sdeltat
      double precision a(2),b(2), da(2)
      double precision f(2), dfdt(2)
      double precision Wt(N,2),dY(N,2)
      double precision dwork(nwork)
      double precision P
      double precision Pp, DXT, DYT
      double precision ALPHA, BETA, GAMMA, OMEGA
      integer          nw,ix,NN
      sdeltat = dsqrt(deltat)
      NN = 2*N
      call rannw (0.d0,sdeltat, ix, Wt,NN, dwork, nwork)        

      

      DO m=1,N

        X=state(m,1)
        Y=state(m,2)

        DXT = ds(m,1)
        DYT = ds(m,2)
        
        ALPHA = par(m,1)
        BETA  = par(m,2)
        GAMMA = par(m,3)
        OMEGA = par(m,4)

        q=par(m,1)

        a(1)   =   -  (Y + BETA)
        a(2)   =    ALPHA * ( X - P(Y))

        da(1)  =    - DYT
        da(2)  =    ALPHA * (DXT - Pp(Y) * DYT)
        
!       simple Euler integration
!       sums
        DO k=1,2
         dY(m,k) = ((- GAMMA * f(k) + a(k)) / OMEGA )*deltat 
     c              + b(k)*Wt(m,k)
         dL(m,k) = da(k)/ OMEGA * deltat
 
        ENDDO  
      ENDDO

      END

     
