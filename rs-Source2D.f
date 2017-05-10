      subroutine rs(t,Nx,Ny,gamma,ro,gammaprime,roprime
     . ,grid,vdx)

      implicit none
      double precision t, factor, aux
      integer Nx, Ny, i, j,ID
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision gammaprime(Nx,Ny),roprime(Nx,Ny)
      double precision f1,f2,Phi,Y
      double precision gLaplace(Nx,Ny)
      double precision xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision vdx(Nx,Ny),vdy
      double precision dke(Nx,Ny),dsigma(Nx,Ny)
      integer grid(Nx,Ny), factor2



      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)


!      call Development(t,Nx,Ny,TS,dke,dsigma)

      do j=1,Ny
       do i=1,Nx
!       Extra variables calculation
        if(grid(i,j) .gt. 0.5)then
            factor=grid(i,j)
        else
            factor=0.0
        endif

        factor2=grid(i,j)
        if (factor2 .lt. 0)then
            factor2=0
        else
            factor2=1
        endif

        vdy=0.0
        aux=gamma(i,j)
        dke(i,j)=1.0
        dsigma(i,j)=1.0
        f1=(1.d0+dk*aux)/(1.d0+aux)
        f2=(dL1+dk*dL2*dc*aux)/(1.d0+dc*aux)
        Y=ro(i,j)*aux/(1.d0+aux)
        Phi=(dlambda1+Y*Y)/(dlambda2+Y*Y)

!       Right hand side
        roprime(i,j)=factor2*(-f1*ro(i,j)+f2*(1.d0-ro(i,j)))
        gammaprime(i,j)=factor/depsilon*
     .              (s2*s1*Phi-dke(i,j)*gamma(i,j))
     .              +depsilon*gLaplace(i,j)
     .          -  (vdx(i,j)*xgradeC(i,j)+vdy*ygradeC(i,j))

!        gammaprime(i,j)=1.0/depsilon*
!     .              (factor*s2*s1*Phi-dke(i,j)*gamma(i,j))
!     .              +depsilon*gLaplace(i,j)
!     .          -  (vdx(i,j)*xgradeC(i,j)+vdy*ygradeC(i,j))


       enddo
      enddo

      return

      end
!      ***********************************************************
!     *************************************************
      subroutine functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      implicit none
      integer Nx, Ny,i,j

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision gamma(Nx,Ny)
      double precision gLaplace(Nx,Ny),xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision gammaim2,gammaim1,gammaip1,gammajm1,gammajp1, d
      double precision gLapX,gLapY,thetai,thetaim1,psii,psiim1
      integer grid(Nx,Ny)

       d=1.0
      call Pillars(Nx,Ny,d,grid)

      do j=1,Ny
       do i=1,Nx
!       No-Flux boundary condition
       if(i .eq. 1 .or. grid(i-1,j) .le.-0.5) then
!       if(i .eq. 1) then
        gammaim2=gamma(i+1,j)
        gammaim1=gamma(i+1,j)
        gammaip1=gamma(i+1,j)

!        gammaim2=0
!        gammaim1=0

       elseif(i .eq. 2 .or. grid(i-2,j) .le.-0.5) then
!       elseif(i .eq. 2) then
        gammaim2=-gamma(i,j)+2*gamma(i-1,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i+1,j)

!        gammaim2=0

       elseif(i .eq. Nx .or. grid(i+1,j) .le.-0.5) then
!       elseif(i .eq. Nx) then
        gammaim2=gamma(i-2,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i-1,j)
       else
        gammaim2=gamma(i-2,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i+1,j)

       endif

!      Periodic Boundary
!       if(i .eq. 1) then
!        gammaim2=gamma(Nx-1,j)
!        gammaim1=gamma(Nx,j)
!        gammaip1=gamma(2,j)
!
!       elseif(i .eq. 2) then
!        gammaim2=gamma(Nx,j)
!        gammaim1=gamma(i-1,j)
!        gammaip1=gamma(i+1,j)
!
!
!       elseif(i .eq. Nx) then
!        gammaim2=gamma(i-2,j)
!        gammaim1=gamma(i-1,j)
!        gammaip1=gamma(1,j)
!       else
!        gammaim2=gamma(i-2,j)
!        gammaim1=gamma(i-1,j)
!        gammaip1=gamma(i+1,j)
!
!       endif


!     Non Flux Boundary Cond
       if(Ny .eq. 1) then
        gammajm1=gamma(i,j)
        gammajp1=gamma(i,j)

!       elseif(j .eq. 1) then
       elseif(j .eq. 1 .or. grid(i,j-1) .le. -0.5) then
        gammajm1=gamma(i,j+1)
        gammajp1=gamma(i,j+1)
!       elseif(j .eq. Ny) then
       elseif(j .eq. Ny .or. grid(i,j+1) .le. -0.5) then
        gammajp1=gamma(i,j-1)
        gammajm1=gamma(i,j-1)
       else
        gammajp1=gamma(i,j+1)
        gammajm1=gamma(i,j-1)
       endif

!     Periodic Boundary
!       if(Ny .eq. 1) then
!        gammajm1=gamma(i,j)
!        gammajp1=gamma(i,j)
!       elseif(j .eq. 1) then
!        gammajm1=gamma(i,Ny)
!        gammajp1=gamma(i,j+1)
!       elseif(j .eq. Ny) then
!        gammajp1=gamma(i,1)
!        gammajm1=gamma(i,j-1)
!       else
!        gammajp1=gamma(i,j+1)
!        gammajm1=gamma(i,j-1)
!       endif


        gLapX=(gammaip1+gammaim1-2*gamma(i,j))/(dx**2)
        gLapY=(gammajp1+gammajm1-2*gamma(i,j))/(dy**2)
        gLaplace(i,j)=gLapX+gLapY


!
!        if(gammaip1 .eq. gamma(i,j)) then
!        thetai=1.d-10
!        else
!        thetai=(gamma(i,j)-gammaim1)/(gammaip1-gamma(i,j))+1.d-10
!        endif
!
!        if(gamma(i,j) .eq. gammaim1)then
!        thetaim1=1.d-10
!        else
!        thetaim1=(gammaim1-gammaim2)/(gamma(i,j)-gammaim1)+1.d-10
!        endif
!
!        psii=max(0.0,min(1.0,1.0/3.0+thetai/6.0,thetai))
!        psiim1=max(0.0,min(1.0,1.0/3.0+thetaim1/6.0,thetaim1))
!
!      xgradeC(i,j)=(1.0-psiim1+psii/thetai)*(-gammaim1+gamma(i,j))/(dx)


        ygradeC(i,j)=(gammajp1-gammajm1)/(2*dy)
        xgradeC(i,j)=(gammaip1-gammaim1)/(2*dx)
!        ygradeC(i,j)=(gammajp1-gamma(i,j))/(dy)
       enddo
      enddo

      return
      end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateCells(h,Nx,Ny,Nc, gamma, ro, cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      integer Nx, Ny,Nc, i, j, k
      double precision cells(Nc,6)
      double precision gamma(Nx,Ny), ro(Nx,Ny)
      double precision gLaplace(Nx,Ny),xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision grad, angle, h,speed


      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

        speed=0.05/sqrt(dke0*0.024)
!        speed=0.0;
      do k= 1,Nc
        i=ceiling(cells(k,1)/dx)
        j=ceiling(cells(k,2)/dy)
        grad=sqrt(xgradeC(i,j)**2 +ygradeC(i,j)**2)
        if(ro(i,j) .ge. 0.6 .and. grad .ge. 1.0)then
            angle=atan2(ygradeC(i,j),xgradeC(i,j))
        cells(k,1)=cells(k,1)+h*speed*cos(angle)
        cells(k,2)=cells(k,2)+h*speed*sin(angle)
        endif
      enddo

      return
      end

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine UpdateDiscreteCells(h,Nx,Ny,Nc, gamma, ro, cells,g0)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      integer Nx, Ny,Nc, i, j, k, oldi, oldj
      double precision cells(Nc,6), g0(Nx,Ny)
      double precision gamma(Nx,Ny), ro(Nx,Ny)
      double precision gLaplace(Nx,Ny),xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision grad, angle, h,speed, newx, newy
      integer grid(Nx,Ny)


      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      call FromCellToGrid(Nx,Ny,Nc,cells,grid)


        speed=0.05/sqrt(dke0*0.024)
!        speed=0.0;

      do k= 1,Nc
        i=ceiling(cells(k,1)/dx)
        j=ceiling(cells(k,2)/dy)
        grad=sqrt(xgradeC(i,j)**2 +ygradeC(i,j)**2)
        if(ro(i,j) .ge. 0.6 .and. grad .ge. 1.0)then
!        if(ro(i,j) .ge. 0.4 .and. grad .ge. 1.0)then
!        if((gamma(i,j)-g0(i,j)) .gt. 0.0 .and. grad .ge. 1.0)then
!          write(6,*) 'Atempt move'
            angle=atan2(ygradeC(i,j),xgradeC(i,j))
            newx=cells(k,1)+h*speed*cos(angle)
            newy=cells(k,2)+h*speed*sin(angle)
            oldi=ceiling(cells(k,1)/dx)
            oldj=ceiling(cells(k,2)/dy)
            i=ceiling(newx/dx)
            j=ceiling(newy/dy)
            if (oldi .eq. i .and. oldj .eq. j)then
                cells(k,1)=newx
                cells(k,2)=newy
            else

                    if(i .gt.Nx)then
                        newx=newx-Nx*dx
                        i=ceiling(newx/dx)
                    endif

                    if(i .lt.1)then
                        newx=newx+Nx*dx
                        i=ceiling(newx/dx)
                    endif

                    if(j .gt. Ny)then
                        newy=newy-Ny*dy
                        j=ceiling(newy/dy)
                    endif

                    if(j .lt. 1)then
                        newy=newy+Ny*dy
                        j=ceiling(newy/dy)
                    endif

                    if(grid(i,j) .eq. 0) then
!                  write(6,*) 'Actual move'
                        cells(k,1)=newx
                        cells(k,2)=newy
                        grid(i,j)=grid(i,j)+1
                        grid(oldi,oldj)=grid(oldi,oldj)-1
                    endif
            endif


        endif
      enddo

      return
      end
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine UpdateWithMemory(h,Nx,Ny,Nc, gamma, ro, cells,g0)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      integer Nx, Ny,Nc, i, j, k, oldi, oldj
      double precision cells(Nc,6), g0(Nx,Ny)
      double precision gamma(Nx,Ny), ro(Nx,Ny)
      double precision gLaplace(Nx,Ny),xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision grad, angle, h,speed, newx, newy, Memory
      integer grid(Nx,Ny)


      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      call FromCellToGrid(Nx,Ny,Nc,cells,grid)


!        speed=0.05/sqrt(dke0*Diffgamma)
        speed=0.02/sqrt(dke0*Diffgamma)
!        speed=0.0;
         Memory=2.0*dk1

      do k= 1,Nc
        i=ceiling(cells(k,1)/dx)
        j=ceiling(cells(k,2)/dy)
        grad=sqrt(xgradeC(i,j)**2 +ygradeC(i,j)**2)
!        if((gamma(i,j)-g0(i,j)) .gt. 0.0 .and. grad .ge. 1.0
!     .    .and. cells(k,5) .le. 0.0 .and. gamma(i,j) .gt. 1.0)then
        if(ro(i,j) .ge. 0.6 .and. grad .ge. 1.0
     .    .and. cells(k,5) .le. 0.0)then
            cells(k,5)=ceiling(Memory/h)
            cells(k,6)=atan2(ygradeC(i,j),xgradeC(i,j))
        endif


        if(cells(k,5) .gt. 0.0)then
            cells(k,5)=cells(k,5)-1;
            angle=cells(k,6)
            newx=cells(k,1)+h*speed*cos(angle)
            newy=cells(k,2)+h*speed*sin(angle)
            oldi=ceiling(cells(k,1)/dx)
            oldj=ceiling(cells(k,2)/dy)
            i=ceiling(newx/dx)
            j=ceiling(newy/dy)
            if (oldi .eq. i .and. oldj .eq. j)then
                cells(k,1)=newx
                cells(k,2)=newy
            else

                    if(i .gt.Nx)then
                        newx=newx-Nx*dx
                        i=ceiling(newx/dx)
                    endif

                    if(i .lt.1)then
                        newx=newx+Nx*dx
                        i=ceiling(newx/dx)
                    endif

                    if(j .gt. Ny)then
                        newy=newy-Ny*dy
                        j=ceiling(newy/dy)
                    endif

                    if(j .lt. 1)then
                        newy=newy+Ny*dy
                        j=ceiling(newy/dy)
                    endif

                    if(grid(i,j) .lt. 3 .and. grid(i,j) .ge. 0) then
!                  write(6,*) 'Actual move'
                        cells(k,1)=newx
                        cells(k,2)=newy
                        grid(i,j)=grid(i,j)+1
                        grid(oldi,oldj)=grid(oldi,oldj)-1
                    endif
            endif


        endif
      enddo

      return
      end
      
