      subroutine ODE(t,Nx,Ny,Nc,gamma,ro,cells)
      
      implicit none
      double precision t
      integer Nx, Ny,Nc,i,j

      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision gamma0(Nx,Ny),ro0(Nx,Ny)
      double precision gammak1(Nx,Ny),rok1(Nx,Ny)
      double precision gammak2(Nx,Ny),rok2(Nx,Ny)
      double precision gammak3(Nx,Ny),rok3(Nx,Ny)
      double precision gammak4(Nx,Ny),rok4(Nx,Ny)
      double precision gammak5(Nx,Ny),rok5(Nx,Ny)
      double precision g1(Nx,Ny),r1(Nx,Ny)
      double precision vdx(Nx,Ny), vdy(Nx,Ny)
      double precision cells(Nc,6)
      integer grid(Nx,Ny)

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision tau,h,err
      integer iteration, index

      tau=0.d0
      h=dt
      call flow(t,Nx,Ny,vdx,vdy)

 13   call FromCellToGrid(Nx,Ny,Nc,cells,grid)

      do j=1,Ny
       do i=1,Nx
        gamma0(i,j)=gamma(i,j)
        ro0(i,j)=ro(i,j)
       enddo
      enddo
      iteration=0

      call rs(t,Nx,Ny,gamma0,ro0,gammak1,rok1,grid,vdx)
!     Runge-Kutta-Merson Method

 16   do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*gammak1(i,j)/3
        ro(i,j)=ro0(i,j)+h*rok1(i,j)/3
       enddo
      enddo

      call rs(t+h/3,Nx,Ny,gamma,ro,gammak2,rok2,grid,vdx)

      do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+gammak2(i,j))/6
        ro(i,j)=ro0(i,j)+h*(rok1(i,j)+rok2(i,j))/6
       enddo
      enddo

      call rs(t+h/3,Nx,Ny,gamma,ro,gammak3,rok3,grid,vdx)

      do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+3*gammak3(i,j))/8
        ro(i,j)=ro0(i,j)+h*(rok1(i,j)+3*rok3(i,j))/8
       enddo
      enddo

      call rs(t+h/2,Nx,Ny,gamma,ro,gammak4,rok4,grid,vdx)

       do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)-3*gammak3(i,j)
     .   +4*gammak4(i,j))/2
        ro(i,j)=ro0(i,j)+h*(rok1(i,j)-3*rok3(i,j)
     .   +4*rok4(i,j))/2
       enddo
      enddo

      call rs(t+h,Nx,Ny,gamma,ro,gammak5,rok5,grid,vdx)

      do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+4*gammak4(i,j)
     .   +gammak5(i,j))/6
        ro(i,j)=ro0(i,j)+h*(rok1(i,j)+4*rok4(i,j)
     .   +rok5(i,j))/6
       enddo
      enddo

      do j=1,Ny
       do i=1,Nx
        g1(i,j)=gamma(i,j)-h*(gammak1(i,j)+gammak2(i,j)+gammak3(i,j)
     .   +gammak4(i,j)+gammak5(i,j))/5-gamma0(i,j)
        r1(i,j)=ro(i,j)-h*(rok1(i,j)+rok2(i,j)+rok2(i,j)+rok3(i,j)
     .   +rok4(i,j)+rok5(i,j))/5-ro0(i,j)
       enddo
      enddo

      err=0.d0
      index=0
      do j=1,Ny
       do i=1,Nx
       err=max(abs(g1(i,j)),abs(r1(i,j)),err)
      if (gamma(i,j) .lt. 0 .or. ro(i,j) .lt. 0)
     . then
      index=1

      endif
       enddo
      enddo

      if (err .gt. tol .or. index .eq. 1) then
      h=h/2
      iteration=iteration+1

        if (iteration .gt. 2) then
            write(6,*) 't =',t,' index =',index, 'iteration=',iteration
        endif
        if (iteration .gt. 40) then
            write(6,*) 'Emergency Exit'
            call EXIT(0)
        endif

      go to 16
      endif


      t=t+h
      tau=tau+h

!      call UpdateDiscreteCells(h,Nx,Ny,Nc, gamma, ro, cells,gamma0)
      call UpdateWithMemory(h,Nx,Ny,Nc, gamma, ro, cells,gamma0)

      h=dt

      if (tau + h .le. tout+tol*dt) then


       go to 13
      elseif(tau .lt. tout-tol*dt)then
         h = tout - tau

         go to 13
      endif


      return
      end




