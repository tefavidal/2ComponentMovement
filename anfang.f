      subroutine anfang(t,Nx,Ny,Nc,gamma,ro,cells)
      
      implicit none

      integer Nx,Ny, Nc
      double precision t

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision dh,tEprime,dk2,dki,dkt,dq,depsilono,dlambda
      double precision dtheta,Vmax,tendprime,toutprime,dtprime
      double precision cells(Nc,6)

      open(8,file = 'startdata',status = 'unknown', form = 'formatted')

   
      read(8,*) dc,dh,tEprime
      read(8,*) dk1,dk2
      read(8,*) dki,dke0,dkt
      read(8,*) dL1,dL2,Diffgamma
      read(8,*) dq,depsilono,dlambda
      read(8,*) dsigma0,dtheta,dalpha
      read(8,*) Vmax
      read(8,*) tendprime,toutprime,dtprime
      read(8,*) amplit
      read(8,*) tol,prob
      read(8,*) isf,itstart,pi

      close(8)
      dk=dk2/dk1
      dlambda1=dlambda*dtheta/depsilono
      dlambda2=(1+dalpha*dtheta)/(depsilono*(1+dalpha))
      s1=dq*dsigma0/(dkt+dki)*dalpha/(dalpha+1)
      s2=dkt/(dke0*dh)
      depsilon=dk1/dke0
      depsilonp=dk1/(dki+dkt)


      vd=Vmax/(dke0*Diffgamma)**0.5
      dt=dtprime*dk1
      tend=tendprime*dk1
      tE=tEprime*dk1
      tout=toutprime*dk1

      write(6,*) 'dk=',dk
      write(6,*) 'dlambda1=',dlambda1
      write(6,*) 'dlambda2=',dlambda2
      write(6,*) 's1=',s1
      write(6,*) 's2=',s2
      write(6,*) 's=',s1*s2
      write(6,*) 'depsilon=',depsilon
      write(6,*) 'depsilonp=',depsilonp
      write(6,*) 'vd=',vd,'dt=',dt,'tend=',tend,'tE=',tE

       dx=0.05*dk1/(dke0*Diffgamma)**0.5
       dy=0.05*dk1/(dke0*Diffgamma)**0.5

      write(6,*)'dimensionless dx=', dx
      write(6,*)'dimensionless dy=', dy

!     %%%%% Optional load of initial conditions %%%%%%%%%%

      open(7,file = 'OutputData2D/Initial-State', err = 20,
     .       status = 'old', form = 'formatted')

      open(8,file = 'OutputData2D/Initial-Positions', err = 20,
     .       status = 'old', form = 'formatted')

      call loadState(t,Nx,Ny,Nc,gamma,ro,cells)
      tE=0
      tend=tend+t
      go to 10

!     %%%%% Random positions one cell per grid point %%%%%%%%%%

 20   call initialDiscreteDistribution(Nx,Ny,Nc,cells)
      call ic(t,Nx,Ny,Nc,gamma,ro,cells)


!      call StartingTime(t,Nx,Ny,TS)

 10   return
      
      end
