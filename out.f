      subroutine out(t,Nx,Ny,Nc,gamma,ro,cells)

      implicit none

      double precision t
      integer Nx, Ny, Nc,i, j

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision cells(Nc,6)
      integer grid (Nx,Ny)

      call FromCellToGrid(Nx,Ny,Nc,cells,grid)

        do j=1,Ny
            do i=1,Nx
            write(10,*) t/dk1,gamma(i,j),ro(i,j),grid(i,j)
            enddo
        enddo

      write(10,*)
      return
      end


!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine outFinal(t,Nx,Ny,Nc,gamma,ro,cells)

      implicit none
      integer Nx,Ny, Nc, i, j
      double precision t

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision cells(Nc,6)



      do i=1,Nx
        do j=1,Ny
            write(42,*) t/dk1,gamma(i,j),ro(i,j)
        enddo
      enddo
      write(42,*)


      do i=1,Nc
        write(43,*) cells(i,1), cells(i,2)
      enddo
      write(43,*)


      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine loadState(t,Nx,Ny,Nc,gamma,ro,cells)

      implicit none

      integer Nx,Ny, Nc, i, j
      double precision t, aux

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision cells(Nc,6)

      do i=1,Nx
        do j=1,Ny
            read(7,*) aux,gamma(i,j),ro(i,j)
        enddo
      enddo
      close(7)
      t=aux*dk1

      do i=1,Nc
        read(8,*) cells(i,1), cells(i,2)
      enddo
      close(8)


      return
      end

