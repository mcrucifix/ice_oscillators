! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!  LORENZ MODEL TO TEST PROCEDURE OF CALCULATING
!  LARGEST LYAPUNOV EXPONENT WITH EULER SCHEME
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


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
       double precision t0,t1, lyap(N)
       double precision par(N,4),state(N,3),dx,dY(N,3), dL(N,3)
       double precision ds(N,3), norm(N)
       double precision t,deltat,b(3)
       double precision f(2),dfdt(2)
       double precision ampsin(na),ampcos(na),ome(na)

       !! scaled such as first term has amplitude 1


       double precision x,dxdt           ! temporary variables for exch.
       deltat  = 0.0002

       if (output.EQ.1) open(20, file = 'trajectory.dat', 
     c                  status='unknown')

       !! adjusts imax to have an int. numb. of timesteps
       IMAX = max(int(((t1-t0)/deltat)+0.5),1)
       deltat  = (t1-t0)/IMAX     

        
       DO I=1,(IMAX)
         t = t0+(I-1)*deltat
         call ddtlin_F(state,N,par,b,f,dfdt,deltat,dY,ix,ds,dl)
         state = state + dY

         if (icalclyap .EQ. 1) THEN
          ds = ds + dL
          if ((MOD(I,100).EQ.0) .OR. (I.EQ.IMAX)) THEN
          norm = (sqrt(ds(:,1)*ds(:,1) +  ds(:,2)*ds(:,2) +
     c                 ds(:,3)*ds(:,3)))
            lyap = lyap + dlog(norm) 
            ds(:,1) = ds(:,1) / norm
            ds(:,2) = ds(:,2) / norm
            ds(:,3) = ds(:,3) / norm
          endif
         endif
         if (output.EQ.1)  write(20,*) t,state, f(1)
       ENDDO
       if (output.EQ.1)  close(20)
      end 
 

      subroutine ddtlin_F(state,N,par,b,f,dfdt,deltat,dY,ix,ds,dl)
      implicit none
      integer K,J,l,M,nwork,il
      integer N
      parameter (nwork=1e4)
      double precision par(N,4),state(N,3), ds(N,3),dl(N,3)
      double precision X,Y,Z
      double precision q,r,eps
      double precision deltat,sdeltat
      double precision a(3),b(3), da(3), f(2), dfdt(2)
      double precision Wt(N,3),dY(N,3)
      double precision dwork(nwork)
      double precision P
      double precision Pp, DXT, DYT, DZT
      double precision SIGMA, RHO, BETA
      integer          nw,ix,NN

      SIGMA = 10.0
      RHO   = 28.0
      BETA  = 8./3.

      DO m=1,N

        X=state(m,1)
        Y=state(m,2)
        Z=state(m,3)

        DXT = ds(m,1)
        DYT = ds(m,2)
        DZT = ds(m,3)
        

        q=par(m,1)

        a(1)   =   SIGMA*(Y - X)
        a(2)   =   X * (RHO - Z) - Y
        a(3)   =   X*Y - BETA * Z

        da(1)  =   SIGMA*(DYT - DXT)
        da(2)  =   DXT*(RHO-Z) - X*DZT - DYT
        da(3)  =   DXT*Y + DYT*X - BETA * DZT
        
!       simple Euler integration
!       sums
        DO k=1,3
         dY(m,k) =  a(k) * deltat 
         dL(m,k) = da(k) * deltat
        ENDDO  
      ENDDO

      END

     
