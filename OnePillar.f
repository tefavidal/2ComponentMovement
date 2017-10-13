      subroutine Pillars(Nx,Ny,d,grid)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      integer Nx,Ny,i,j

      double precision d,r,r0
      integer grid(Nx,Ny)
      integer x0, y0

!      d=1.0
      r=1.0/2.0*dk1/(dke0*Diffgamma)**0.5
      r=r*r

      x0=ceiling(Nx/2.0)
      y0=ceiling(Ny/2.0)

      do j= 1,Ny
        do i= 1, Nx
            r0=(i-x0)*(i-x0)*dx*dx + (j-y0)*(j-y0)*dy*dy
            if (r0 .le. r) then
                grid(i,j)=-1
            else
                grid(i,j)=0
                endif
        enddo
      enddo

      return
      end
