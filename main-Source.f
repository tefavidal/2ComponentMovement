      program main
      
      implicit none


      INTEGER, PARAMETER :: Nx=200
      INTEGER, PARAMETER :: Ny=200
!         FOUR PILLARS   40 PERCENT
!      INTEGER, PARAMETER :: Nc=15512

!           ONE PILLAR 40 PERCENT
!      INTEGER, PARAMETER :: Nc=15878

!           ONE PILLAR 35 PERCENT
!      INTEGER, PARAMETER :: Nc=13893

!           ONE PILLAR 30 PERCENT
!      INTEGER, PARAMETER :: Nc=11908

!           ONE PILLAR 25 PERCENT
!      INTEGER, PARAMETER :: Nc=9923

!           ZERO PILLAR 20 PERCENT
      INTEGER, PARAMETER :: Nc=8000

!          WALL 40 PERCEN
!      INTEGER, PARAMETER :: Nc=14480

!          WALL 3mm 40 PERCEN
!      INTEGER, PARAMETER :: Nc=11280

!     Nc decided as Nx*Ny*percentageofcells

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision t
      double precision cells(Nc,6)
      integer pgbeg, i, j, counter, factor, k
      integer grid(Nx,Ny)
      character(len=4) ct1
      character(len=17) ct2
      character(len=21) ct3
      real aux

      double precision gamma0(10),ro0(10)

  ! ----- variables for portable seed setting -----
      INTEGER :: i_seed
      INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
      INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----

 ! ----- Set up random seed portably -----
      CALL RANDOM_SEED(size=i_seed)
      ALLOCATE(a_seed(1:i_seed))
      CALL RANDOM_SEED(get=a_seed)
      CALL DATE_AND_TIME(values=dt_seed)
      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)
     .*dt_seed(6)
       write(6,*) 'seed=',a_seed(i_seed)
      CALL RANDOM_SEED(put=a_seed)
      DEALLOCATE(a_seed)
  ! ----- Done setting up random seed ----



      t=0.d0
      counter=0;
      call anfang(t,Nx,Ny,Nc,gamma,ro,cells)
          open(10,file ='OutputData2D/data   0'
     .     ,status = 'unknown',form = 'formatted')
              call out(t,Nx,Ny,Nc,gamma,ro,cells)
          close(10)
      
      ct2='OutputData2D/data'
 5    continue


      call ODE(t,Nx,Ny,Nc,gamma,ro,cells)
      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1
      counter=counter+1
      write(ct1,'(I4)') counter
      ct3 = ct2 // ct1

      write(6,*) ct3

          open(10,file =ct3
     .     ,status = 'unknown',form = 'formatted')
              call out(t,Nx,Ny,Nc,gamma,ro,cells)
          close(10)


      if(mod(counter,50) .eq. 0)then
        write(6,*) 'Random Firing'

!%%%%%%Second Iteration with fixed firing cells
        do k=1,Nc
            if(cells(k,3) .gt. 0.5)then
                i=ceiling(cells(k,1)/dx)
                j=ceiling(cells(k,2)/dy)
                gamma(i,j)=gamma(i,j)+2
            endif
        enddo
      endif

!%%%%%%%%%%%First Iteration any cell might fire

!        call FromCellToGrid(Nx,Ny,Nc,cells,grid)
!        do j=1,Ny
!            do i=1,Nx
!                    if(grid(i,j) .gt. 0.5)then
!                     do k=1,grid(i,j)
!                        call random_number(aux)
!                        if (aux .lt. 0.1) then
!                            gamma(i,j)=gamma(i,j)+2
!                        endif
!                     enddo
!                    endif
!            enddo
!        enddo
!      endif



      if (t+dt .lt. tend) then
         goto 5

      endif

      close(10)

!

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     WRITES FINAL STATE
      open(42,file ='OutputData2D/Final-State'
     . ,status = 'unknown',form = 'formatted')

      open(43,file ='OutputData2D/Final-Positions'
     . ,status = 'unknown',form = 'formatted')

      call outFinal(t,Nx,Ny,Nc,gamma,ro,cells)
      close(42)
      close(43)

      end
      


