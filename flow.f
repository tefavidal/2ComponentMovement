      subroutine flow(t,Nx,Ny,vdx,vdy)

      implicit none
      double precision t
      integer Nx, Ny, i, j
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision vdx(Nx,Ny),vdy(Nx,Ny)



!     %%%%%%%%%%%%%%%%parabolic flow
!      do i=1,Nx
!       do j=1,Ny
!
!      vdx(i,j)=vd*
!     . (1.d0-(j-1-(Ny-1.d0)/2)**2/((Ny-1.d0)/2)**2)
!        vdy(i,j)=0
!       enddo
!      enddo
!     %%%%%%%%%%%%%%%%Realistic Flow
!      if(Ny .eq. 80)then
!        do i=1,Nx
!            vdx(i,1)=0
!            vdx(i,2)=vd*0.7946
!            vdx(i,3)=vd*0.9720
!            vdx(i,4)=vd*0.9976
!            vdx(i,5)=vd*0.9999
!
!            vdx(i,76)=vd*0.9999
!            vdx(i,77)=vd*0.9976
!            vdx(i,78)=vd*0.9720
!            vdx(i,79)=vd*0.7946
!            vdx(i,80)=0
!
!            vdy(i,1)=0
!            vdy(i,2)=0
!            vdy(i,3)=0
!            vdy(i,4)=0
!            vdy(i,5)=0
!
!            vdy(i,76)=0
!            vdy(i,77)=0
!            vdy(i,78)=0
!            vdy(i,79)=0
!            vdy(i,80)=0
!
!            do j=6,75
!                vdx(i,j)=vd
!                vdy(i,j)=0
!            enddo
!        enddo
!      elseif(Ny .eq. 100)then
!        do i=1,Nx
!            vdx(i,1)=0
!            vdx(i,2)=vd*0.7004
!            vdx(i,3)=vd*0.9317
!            vdx(i,4)=vd*0.9895
!            vdx(i,5)=vd*0.9987
!            vdx(i,6)=vd*0.9999
!
!            vdx(i,95)=vd*0.9999
!            vdx(i,96)=vd*0.9987
!            vdx(i,97)=vd*0.9895
!            vdx(i,98)=vd*0.9317
!            vdx(i,99)=vd*0.7004
!            vdx(i,100)=0
!
!            vdy(i,1)=0
!            vdy(i,2)=0
!            vdy(i,3)=0
!            vdy(i,4)=0
!            vdy(i,5)=0
!            vdy(i,6)=0
!
!            vdy(i,95)=0
!            vdy(i,96)=0
!            vdy(i,97)=0
!            vdy(i,98)=0
!            vdy(i,99)=0
!            vdy(i,100)=0
!
!            do j=7,94
!                vdx(i,j)=vd
!                vdy(i,j)=0
!            enddo
!        enddo
!
!
!      elseif(Ny .eq. 120)then
!        do i=1,Nx
!            vdx(i,1)=0
!            vdx(i,2)=vd*0.6350
!            vdx(i,3)=vd*0.8903
!            vdx(i,4)=vd*0.9708
!            vdx(i,5)=vd*0.9943
!            vdx(i,6)=vd*0.9991
!            vdx(i,7)=vd*0.9999
!
!            vdx(i,114)=vd*0.9999
!            vdx(i,115)=vd*0.9991
!            vdx(i,116)=vd*0.9943
!            vdx(i,117)=vd*0.9708
!            vdx(i,118)=vd*0.8903
!            vdx(i,119)=vd*0.6350
!            vdx(i,120)=0
!
!            vdy(i,1)=0
!            vdy(i,2)=0
!            vdy(i,3)=0
!            vdy(i,4)=0
!            vdy(i,5)=0
!            vdy(i,6)=0
!            vdy(i,7)=0
!
!            vdy(i,114)=0
!            vdy(i,115)=0
!            vdy(i,116)=0
!            vdy(i,117)=0
!            vdy(i,118)=0
!            vdy(i,119)=0
!            vdy(i,120)=0
!
!            do j=8,113
!                vdx(i,j)=vd
!                vdy(i,j)=0
!            enddo
!        enddo
!
!
!
!      elseif(Ny .eq. 160)then
!        do i=1,Nx
!            vdx(i,1)=0
!            vdx(i,2)=vd*0.5144
!            vdx(i,3)=vd*0.7946
!            vdx(i,4)=vd*0.9165
!            vdx(i,5)=vd*0.9721
!            vdx(i,6)=vd*0.9910
!            vdx(i,7)=vd*0.9977
!            vdx(i,8)=vd*0.9995
!            vdx(i,9)=vd*0.9999
!
!            vdx(i,152)=vd*0.9999
!            vdx(i,153)=vd*0.9995
!            vdx(i,154)=vd*0.9977
!            vdx(i,155)=vd*0.9910
!            vdx(i,156)=vd*0.9721
!            vdx(i,157)=vd*0.9165
!            vdx(i,158)=vd*0.7946
!            vdx(i,159)=vd*0.5144
!            vdx(i,160)=0
!
!            vdy(i,1)=0
!            vdy(i,2)=0
!            vdy(i,3)=0
!            vdy(i,4)=0
!            vdy(i,5)=0
!            vdy(i,6)=0
!            vdy(i,7)=0
!            vdy(i,8)=0
!            vdy(i,9)=0
!
!            vdy(i,152)=0
!            vdy(i,153)=0
!            vdy(i,154)=0
!            vdy(i,155)=0
!            vdy(i,156)=0
!            vdy(i,157)=0
!            vdy(i,158)=0
!            vdy(i,159)=0
!            vdy(i,160)=0
!
!            do j=10,151
!                vdx(i,j)=vd
!                vdy(i,j)=0
!            enddo
!        enddo
!
!!
!      elseif(Ny .eq. 200)then
!        do i=1,Nx
!            vdx(i,1)=0
!            vdx(i,2)=vd*0.4347
!            vdx(i,3)=vd*0.7004
!            vdx(i,4)=vd*0.8517
!            vdx(i,5)=vd*0.9317
!            vdx(i,6)=vd*0.9708
!            vdx(i,7)=vd*0.9885
!            vdx(i,8)=vd*0.9958
!            vdx(i,9)=vd*0.9986
!            vdx(i,10)=vd*0.9996
!            vdx(i,11)=vd*0.9999
!
!            vdx(i,190)=vd*0.9999
!            vdx(i,191)=vd*0.9996
!            vdx(i,192)=vd*0.9986
!            vdx(i,193)=vd*0.9958
!            vdx(i,194)=vd*0.9885
!            vdx(i,195)=vd*0.9708
!            vdx(i,196)=vd*0.9317
!            vdx(i,197)=vd*0.8517
!            vdx(i,198)=vd*0.7004
!            vdx(i,199)=vd*0.4347
!            vdx(i,200)=0
!
!            vdy(i,1)=0
!            vdy(i,2)=0
!            vdy(i,3)=0
!            vdy(i,4)=0
!            vdy(i,5)=0
!            vdy(i,6)=0
!            vdy(i,7)=0
!            vdy(i,8)=0
!            vdy(i,9)=0
!            vdy(i,10)=0
!            vdy(i,11)=0
!
!            vdy(i,190)=0
!            vdy(i,191)=0
!            vdy(i,192)=0
!            vdy(i,193)=0
!            vdy(i,194)=0
!            vdy(i,195)=0
!            vdy(i,196)=0
!            vdy(i,197)=0
!            vdy(i,198)=0
!            vdy(i,199)=0
!            vdy(i,200)=0
!
!            do j=12,189
!                vdx(i,j)=vd
!                vdy(i,j)=0
!            enddo
!        enddo
!!
!
!
!      else
        do i=1,Nx
            do j=1,Ny
                vdx(i,j)=vd
                vdy(i,j)=0
            enddo
        enddo

!      endif

      return
      end
