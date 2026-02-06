  subroutine idata_TOVstar2

! *********************************
! ***   TOV2 STAR INITIAL DATA   ***
! *********************************

! Include modules.

  use procinfo
  use param
  use arrays
  use derivatives

! Extra variables.

  implicit none
  logical foundsurf                     ! Did we find the star's surface?
  logical found_P_null                  ! Did we find when P is negative?

  integer i,l,iter,j                    ! Counters.
  integer imin                          ! Leftmost grid point.
  integer iaux                          ! Auxiliary quantity.
  integer irad, irad_p                  ! Position of star's surface, and core surface.

  integer :: size_filter
  real(8) sigma,sigma2, conv_sum,conv_sum2,conv_sum3, conv_sum4
  real(8), dimension (:), allocatable :: kernel(:),kernel2


  integer :: xf, half_size
  real(8) :: sum_val, sum_val2, g1, g2, g

  ! Declare variables for writing Delete 
  integer :: unit_number,uni
  character(len=20) :: filename
  real(8) auxxx,auxx2


  real(8) r0,delta                      ! Local radius and grid spacing.
  real(8) A0,rho0_c,P0,rho0,rho_c       ! Initial values of variables.
  real(8) k11,k12,k13,k14               ! Runge-Kutta sources for A.
  real(8) k21,k22,k23,k24               ! Runge-Kutta sources for rho0.
  real(8) P0_rk,rho0_rk,Pi              ! Runge-Kutta values of variables.
  real(8) drp,drho0dr                   ! Functions for sources of differential equations.
  real(8) rm,aux                        ! Auxiliary quantities.
  real(8) half,smallpi                  ! Numbers.
  real(8) A_function,e_function
  real(8) Rho_function
  real(8) rad_P_null
  real(8) TOV_rad2
  real(8) h,h0,vs0


  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rr             ! Radial coordinate.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: A_g            ! Radial metric global array.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: e_g            ! Specifical energy.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rho0_g         ! Rest mass density global array.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: rho_g          ! Density global array.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: P_g            ! Pressure global array.
  real(8), dimension (0:Nl-1,1-ghost:Nrtotal) :: fluid_vsAux    ! Sound velocity auxilary.

! NUMBERS   
  half = 0.5d0
  smallpi = acos(-1.d0)
 
! SURFACE PARAMETERS
! Initialize flag and surface position.

  foundsurf = .false.
  found_P_null = .false.
  irad = 0                ! To save the location i for the star radius
  irad_p = 0.d0           ! To save the location i for the core radius
  
  print *, "Solving initial data for a TOV2 star ..."
  print *


! FIND GRID POINT POSITIONS  

  do l=0,Nl-1
     do i=1-ghost,Nrtotal
        rr(l,i) = (dble(Nmin(rank) + i) - half)*dr(l)
     end do
  end do


! INTEGRATION
! We solve the system of equations using fourth order Runge--Kutta.
! Initialize arrays.
  rho_g = 1.d0
  A_g = 1.d0
  e_g = 1.d0

  rho0_g = 0.d0
  P_g = 0.d0 

  Omega_TOV2 = 1.d-20
  fluid_vs = 1.d-4
  rho0_c = TOV_rho_c/(1.d0 + TOV_e_c)


  vs0= sqrt(abs(2.d0*smallpi*(TOV_P0+TOV_rho_c/3.d0)*(TOV_P0+TOV_rho_c) / (C2_tov)))
  fluid_vs(l,0)=vs0


  if (rank==0) then ! ONLY PROCESSOR 0 SOLVES THE ODE's

!    LOOP OVER GRID LEVELS
!    We solve from fine to coarse grid.

     do l=Nl-1,0,-1

!       Find initial point. Only the finest grid integrates from the origin.

        if (l==Nl-1) then
           imin = 1
        else
           imin = Nrtotal/2
        end if

!       For coarse grids we interpolate the initial point (which is half way 
!       through the grid) from the last points of the higher resolution grid.  
!       Here we use cubic interpolation.
        if (l<Nl-1) then

           A_g(l,imin-1) = (9.d0*(A_g(l+1,Nrtotal-2)+A_g(l+1,Nrtotal-3)) &
                        - (A_g(l+1,Nrtotal-4)+A_g(l+1,Nrtotal-1)))/16.d0

           rho0_g(l,imin-1) = (9.d0*(rho0_g(l+1,Nrtotal-2)+rho0_g(l+1,Nrtotal-3)) &
                            - (rho0_g(l+1,Nrtotal-4)+rho0_g(l+1,Nrtotal-1)))/16.d0

           rho_g(l,imin-1) = (9.d0*(rho_g(l+1,Nrtotal-2)+rho_g(l+1,Nrtotal-3)) &
                            - (rho_g(l+1,Nrtotal-4)+rho_g(l+1,Nrtotal-1)))/16.d0

           P_g(l,imin-1) = (9.d0*(P_g(l+1,Nrtotal-2)+P_g(l+1,Nrtotal-3)) &
                            - (P_g(l+1,Nrtotal-4)+P_g(l+1,Nrtotal-1)))/16.d0

           e_g(l,imin-1) = (9.d0*(e_g(l+1,Nrtotal-2)+e_g(l+1,Nrtotal-3)) &
                            - (e_g(l+1,Nrtotal-4)+e_g(l+1,Nrtotal-1)))/16.d0

           Omega_TOV2(l,imin-1) = (9.d0*(Omega_TOV2(l+1,Nrtotal-2)+Omega_TOV2(l+1,Nrtotal-3)) &
                            - (Omega_TOV2(l+1,Nrtotal-4)+Omega_TOV2(l+1,Nrtotal-1)))/16.d0

           fluid_vs(l,imin-1) = (9.d0*(fluid_vs(l+1,Nrtotal-2)+fluid_vs(l+1,Nrtotal-3)) &
                            - (fluid_vs(l+1,Nrtotal-4)+fluid_vs(l+1,Nrtotal-1)))/16.d0

        end if


!       ************************************
!       ***   FOURTH ORDER RUNGE-KUTTA   ***
!       ************************************

        do i=imin,Nrtotal

!          Grid spacing and values at first point
!          if we start from the origin (finer grid).

           if (i==1) then

!             For the first point we use dr/2.
              delta = half*dr(l)
              r0    = 0.d0 !then r0 is the previus point.

!             Values that evolve.
              rho0 = rho0_c
              P0   = TOV_P0
              fluid_vs(l,i)=vs0
              Omega_TOV2(l,i)=TOV_P0/rho0

!          Grid spacing and values at previous grid point.

           else
              delta = dr(l)
              r0    = rr(l,i-1)
              rho0 = rho0_g(l,i-1)
              P0   = P_g(l,i-1)

              if (P_g(l,i-1)<0) then  ! Here we make sure that
                P_g(l,i-1) = 0.d0     ! the pressure is not
                P_g(l,i) = 0.d0       ! negative at any point.
                P0 = 0.d0
              end if
           end if

!          We only solve the equations inside the star.
!          Section I. Core:
           if (.not.found_P_null) then

!             I) First Runge-Kutta step.
              rm = r0
              P0_rk = P0
              rho0_rk = rho0
           
              k21 = delta*drp(rm,P0_rk) 
              k11 = delta*drho0dr(rm,P0_rk,rho0_rk)

!             II) Second Runge-Kutta step.
              rm      = r0   + half*delta
              P0_rk   = P0   + half*k21
              rho0_rk = rho0 + half*k11

              k22 = delta*drp(rm,P0_rk)
              k12 = delta*drho0dr(rm,P0_rk,rho0_rk)

!             III) Third Runge-Kutta step.
              P0_rk   = P0   + half*k22
              rho0_rk = rho0 + half*k12

              k23 = delta*drp(rm,P0_rk)
              k13 = delta*drho0dr(rm,P0_rk,rho0_rk)
 
!             IV) Fourth Runge-Kutta step.
              rm      = r0 + delta
              P0_rk   = P0   + k23
              rho0_rk = rho0 + k13

              k24 = delta*drp(rm,P0_rk)
              k14 = delta*drho0dr(rm,P0_rk,rho0_rk)

!             Advance eho0 to next grid point.
              Pi          = P0   + (k21 + 2.d0*(k22 + k23) + k24)/6.d0
              rho0_g(l,i) = rho0 + (k11 + 2.d0*(k12 + k13) + k14)/6.d0

!             Filling extra arrays.
              A_g(l,i) = A_function(rm)
              rho_g(l,i) = Rho_function(rm)
              e_g(l,i)   = e_function(rho_g(l,i), rho0_g(l,i))

!             If the P0 is positive update the pressure.
!             RK4 is also cutt off when we reach the star radius.
              if (P0>0.d0 .and. rr(l,i)<=TOV_Rt ) then
                 P_g(l,i) = Pi
                 Omega_TOV2(l,i) = P_g(l,i)/rho0_g(l,i)

!             Otherwise we are at the end of the core, so we set the pressure to zero.

              else
                 P_g(l,i) = 0.d0
                 Omega_TOV2(l,i) = 0.0
!                Set flag to true and save this position.

                 found_P_null = .true.
                 irad_p = i-1

!                Find the position using linear extrapolation from
!                previous two points.

                 rad_P_null = rr(l,irad_p-1) + dr(l)*P_g(l,irad_p-i)/(P_g(l,irad_p) - P_g(l,irad_p-1))

                 write(*,'(A,E19.12)') ' Point at which pressure is zero    = ',rad_P_null

                 ! If we reach the radius of the star, 
                 ! then we go directly to the outer part, section III.
                 if (rr(l,i)>=TOV_Rt) then
                    foundsurf =.true.
                 end if
              end if

!          Section II. Dust
!          It just solves the energy density. The pressure is zer
           else if (.not.foundsurf) then

              P0 = 0.d0
              P_g(l,i) = 0.d0
              Omega_TOV2(l,i)=0.d0

!             I) First Runge-Kutta step.
              rm = r0
              rho0_rk = rho0
              P0_rk = P0
           
              k21 = delta*drp(rm,P0_rk) 
              k11 = delta*drho0dr(rm,P0_rk,rho0_rk)

!             II) Second Runge-Kutta step.
              rm = r0 + half*delta
              rho0_rk = rho0 + half*k11
              P0_rk   = P0   + half*k21

              k22 = delta*drp(rm,P0_rk)
              k12 = delta*drho0dr(rm,P0_rk,rho0_rk)

!             III) Third Runge-Kutta step.
              rho0_rk = rho0 + half*k12
              P0_rk    = P0   + half*k22

              k23 = delta*drp(rm,P0_rk)
              k13 = delta*drho0dr(rm,P0_rk,rho0_rk)

!             IV) Fourth Runge-Kutta step.
              rm = r0 + delta
              P0_rk   = P0   + k23
              rho0_rk = rho0 + k13

              k24 = delta*drp(rm,P0_rk)
              k14 = delta*drho0dr(rm,P0_rk,rho0_rk)

!             Advance variables to next grid point.
              rho0_g(l,i) = rho0 + (k11 + 2.d0*(k12 + k13) + k14)/6.d0

!             Extra arrays
              A_g(l,i) = A_function(rm)
              rho_g(l,i) = Rho_function(rm)
              e_g(l,i)   = e_function(rho_g(l,i), rho0_g(l,i))

!             Have we reached the outside of the star if...
              if (rho_g(l,i) <= fluid_atmos) then

!                Set flag to true and save position
!                of last grid point inside star.

                 foundsurf = .true.

                 irad = i-1                

!                Find star radius using linear extrapolation from
!                previous two points.

                 TOV_rad  = rr(l,irad-1) + dr(l)*rho0_g(l,irad-i)/(rho0_g(l,irad) - rho0_g(l,irad-1))
                 TOV_rad2 = rr(l,i)

                 write(*,'(A,E19.12)') ' Radius of TOV star    = ',TOV_rad

!                Set fluid density equal to zero.

                 rho0_g(l,i) = 0.d0
                 rho_g(l,i)  = 0.d0
                 e_g(l,i)    = 0.d0
                 A_g(l,i)    = A_function(rr(l,i))

               end if

!           Section III. Outside.
            else
                rm = r0 + delta
                P_g(l,i)    = 0.d0
                rho0_g(l,i) = 0.d0

                rho_g(l,i)  = 0.d0
                e_g(l,i)    = 0.d0 
                A_g(l,i)    = A_function(rr(l,i))

            end if

        end do

!       GHOST ZONES   

        do i=1,ghost 
           P_g(l,1-i)    = P_g(l,i)
           A_g(l,1-i)    = A_g(l,i)
           rho0_g(l,1-i) = rho0_g(l,i)
           rho_g(l,1-i)  = rho_g(l,i)
           Omega_TOV2(l,1-i)  = Omega_TOV2(l,i)
           e_g(l,1-i)    = e_g(l,i)
        end do
     end do

!    RESTRICT TO COARSE GRIDS

     do l=Nl-1,1,-1
        do i=1,Nrtotal-ghost,2
           iaux = i/2 + 1
           rm = rr(l-1,iaux)

           A_g(l-1,iaux) = (9.d0*(A_g(l,i)+A_g(l,i+1)) - (A_g(l,i-1)+A_g(l,i+2)))/16.d0
           rho0_g(l-1,iaux) = (9.d0*(rho0_g(l,i)+rho0_g(l,i+1)) - (rho0_g(l,i-1)+rho0_g(l,i+2)))/16.d0
           P_g(l-1,iaux) = (9.d0*(P_g(l,i)+P_g(l,i+1)) - (P_g(l,i-1)+P_g(l,i+2)))/16.d0
           rho_g(l-1,iaux) = (9.d0*(rho_g(l,i)+rho_g(l,i+1)) - (rho_g(l,i-1)+rho_g(l,i+2)))/16.d0
           e_g(l-1,iaux) = (9.d0*(e_g(l,i)+e_g(l,i+1)) - (e_g(l,i-1)+e_g(l,i+2)))/16.d0
           fluid_vs(l-1,iaux) = (9.d0*(fluid_vs(l,i)+fluid_vs(l,i+1)) - (fluid_vs(l,i-1)+fluid_vs(l,i+2)))/16.d0
        end do

!       Fix ghost zones.
        do i=1,ghost
           A_g(l-1,1-i)    = A_g(l-1,i)
           rho0_g(l-1,1-i) = rho0_g(l-1,i)

           P_g(l-1,1-i)    = P_g(l-1,i)
           rho_g(l-1,1-i)  = rho_g(l-1,i)
           fluid_vs(l-1,1-i) = fluid_vs(l-1,i)
           Omega_TOV2(l-1,1-i)  = Omega_TOV2(l-1,i)
           e_g(l-1,1-i)  = e_g(l-1,i)
        end do
     end do

  end if


! DISTRIBUTE SOLUTION AMONG PROCESSORS 

! For parallel runs, when we get here the solution is known only on 
! processor zero for the full size arrays with dimensions Nrtotal.  
! We must now distribute the solution among all other processors.

  if (size==1) then
     A = A_g
     fluid_rho = rho0_g
  else
     call distribute(0,Nl-1,A,A_g)
     call distribute(0,Nl-1,fluid_rho,rho0_g)
  end if


! FIND ALL OTHER FLUID VARIABLES

! Add atmosphere (we just add it everywhere).
  fluid_rho = fluid_rho + fluid_atmos

! Find pressure and internal energy.
  fluid_p = P_g
  fluid_e = e_g    


! Find enthalpy: h = 1 + e + p/rho0.
  fluid_h = 1.d0 + fluid_e + fluid_p/fluid_rho

! Set fluid velocity to zero.
  fluid_v = 0.d0

! Set Sound Velocity to second order
 do l=0,Nl-1
    do i=1-ghost,Nrtotal

        ! At ghosts and i=1 points we set vs=vs_0
        if (i<=1) then
            fluid_vs(l,i) = vs0

        ! At the Dust zone we set vs=0
        else if (rr(l,i)>=rad_P_null) then
            fluid_vs(l,i) = 0.d0

        ! When Nrtotal doesn't find the dust part
        else if (i==Nrtotal .and. P_g(l,i)/= 0.d0 ) then
            fluid_vs(l,i) = sqrt(abs(  ( P_g(l,i-1)   - P_g(l,i) ) & 
                                     / ( rho_g(l,i-1) - rho_g(l,i) )  ))
        else 
            fluid_vs(l,i) = sqrt(abs(  ( P_g(l,i-1)   - P_g(l,i+1) ) & 
                                     / ( rho_g(l,i-1) - rho_g(l,i+1) )  ))
        end if
    end do
  end do 

! Set p_eos and rho_eos
  if (fluid_EOS=="TOV2_EoS") then
     do l=0,Nl-1
        do i=1,Nrtotal
            if (P_g(l,i)>0.d0) then
                p_eos(l,i)   = P_g(l,i)
                rho_eos(l,i) = rho0_g(l,i)
            else
                p_eos(l,i) = 0.d0
                rho_eos(l,i)=0.d0
            end if

            if (i == 1) rho_maxim = rho0_g(l,i)
            if (rho0_g(l,i) > rho_maxim) rho_maxim = rho0_g(l,i)
        end do
     end do
  end if

! Set Lorentz factor to one.
  fluid_W = 1.d0

! Conserved quantities.
  fluid_cD = fluid_rho
  fluid_cE = fluid_rho*fluid_e
  fluid_cS = 0.d0

! Initialize artificial viscosity to zero.
  fluid_q = 0.d0

! GAUSSIAN FILTER          
  if (filterTOV2) then
      print*, "Entering Gaussian Filter(rho & omega)"
      size_filter = ghost*2+1
      print*,"Filter Size=",size_filter

      sigma  = 2.0d0
      half_size = size_filter / 2
      sum_val  = 0.0d0

      allocate(kernel(size_filter))

      do xf = -half_size,half_size
            g1=1.0d0/(sigma*sqrt(2.0d0*smallpi))
            g2=-real(xf**2,8)/(2.0d0*sigma**2)
            g=g1*exp(g2)

            kernel(xf+half_size+1)=g
            sum_val = sum_val+g        
      end do

    ! Nomalize the kernel
      kernel  = kernel  / sum_val
      do l=Nl-1,0,-1
      do i = 1,Nrtotal-ghost
        conv_sum = 0
        conv_sum2= 0
        conv_sum3= 0
        conv_sum4= 0

        do j = 1,size_filter
            conv_sum  = conv_sum  +  fluid_rho(l,i+j-ghost-1)*kernel(j)
            conv_sum2 = conv_sum2  +  fluid_e(l,i+j-ghost-1)*kernel(j)
            conv_sum3 = conv_sum3 + Omega_TOV2(l,i+j-ghost-1)*kernel(j)
            conv_sum4 = conv_sum4 +   fluid_p(l,i+j-ghost-1)*kernel(j)
        end do

        fluid_rho(l,i)  = conv_sum
        fluid_e(l,i)    = conv_sum2
        Omega_TOV2(l,i) = conv_sum3
        fluid_p(l,i)    = conv_sum4
      
      end do
      end do
  end if

! FIND LAPSE 

  do l=0,Nl-1
     call auxiliary(l) !auxiliary.f90 -> fluidprimitive
  end do

  call alphamaximal(0,Nl-1,"robin",1.d0)

  print *
  print *, 'idata_TOVstar2 module done...'
  print *, '-----------------------------'
  print *

  end subroutine idata_TOVstar2


! ******************************
! ***   Auxiliary functions  ***
! ******************************

! Ansatz. rho_value
! It gives the total energy density, rho, at the point rm. 
  function Rho_function(rm)

      use param
      implicit none
      real(8) Rho_function
      real(8) rm
      
      if (rm<TOV_Rt) then
        Rho_function = TOV_rho_c - C2_tov*rm**2 - C4_tov*rm**4
      else
        Rho_function = 0.d0
      end if

      if (Rho_function < 0) then
        print*, "Warning! Funtion Rho_value is giving a negative value!",rm
      end if
  end function Rho_function

! A_function   
! It gives the radial metric funtion A = (1-2m(r)/r)^-1
! At r=0, A=1. If rm > Rt, m = M
  function A_function(rm)

      use param
      implicit none
      real(8) A_function
      real(8) rm
      real(8) smallpi
      real(8) m
      smallpi = acos(-1.d0)
      
      if (rm==0.d0) then
         A_function = 1.d0
      else
         if (rm < TOV_Rt) then
            m = 4.d0*smallpi * rm**3 * ( (TOV_rho_c/3.d0) - (C2_tov * (rm**2)/5.d0) &
                                                          - (C4_tov * (rm**4)/7.d0) )  
         else
            m = 4.d0*smallpi * TOV_Rt**3 * ( (TOV_rho_c/3.d0) - (C2_tov * (TOV_Rt**2)/5.d0) &
                                                              - (C4_tov * (TOV_Rt**4)/7.d0) )  
         end if

         A_function = (1.d0 - 2.0*m/rm)**(-1)
      end if

  end function A_function

! **********************************
! ***     Principal Funtions     ***
! **********************************

! The radial derivative of P comes from the TOV equations.
! Explicit expresion is needed to eliminate the problem at r=0 
! The explicit form comes form the use of m and rho.

  function drp(rm,P)

      use param
      implicit none
      real(8) drp
      real(8) rm,P
      real(8) smallpi
      real(8) ee,a,b,c,d
      smallpi = acos(-1.d0)

      ee = TOV_rho_c/3.d0 - C2_tov*(rm**2)/5.d0 - C4_tov*(rm**4)/7.d0

      a = P + TOV_rho_c - C2_tov*(rm**2) - C4_tov*(rm**4) 
      b = 4.d0*smallpi*rm
      c = ee + P
      d = 1.d0 - 8.d0*smallpi*(rm**2)*ee
      drp = - a * b * c / d

      if (d == 0.d0) then
        print*,"Warning, problems in dp/dr. Dividing by zero"
      end if
  end function drp


! Radial rest energy density derivative  d rho0/dr
  function drho0dr(rm,P,rho0)
 
      use param
      implicit none
      real(8) rm,P,rho0,drho0dr
      real(8) Rho_function

      if (rm<TOV_Rt) then
        drho0dr = rho0*(-2.d0*C2_tov*rm - 4.d0*C4_tov*(rm**3) ) / (P + Rho_function(rm))
      else
        drho0dr = 0.d0
      end if

      if (P==0.d0 .and. Rho_function(rm)==0.d0) then 
      print*, 'Careful. rho and p en d rho0/dr are zero.'
      print*, '          at rm= ', rm
      drho0dr = 1.d0
      end if

  end function drho0dr


! Specific internal energy function e. From rho = rho0(1+e)
  function e_function(rho,rho0)
      implicit none
      real(8) rho,rho0
      real(8) e_function

      if (rho==0.d0 .and. rho0 ==0.d0) then
        e_function = 0.d0
      else if (rho0>rho) then
        e_function = 0.d0
      else
        e_function = (rho/rho0 - 1.d0)
      end if

  end function e_function

