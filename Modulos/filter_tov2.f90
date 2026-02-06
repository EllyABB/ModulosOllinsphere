
subroutine filter_tov2()
    use param
    use arrays

    implicit none

    ! Variables
    integer :: xf, i, j, l,n_vs
    integer :: size_filter, half_size
    integer :: idx
    integer :: start_30,start_10,start_60

    real(8) :: sum_x2,sum_x2y,sum_x4,sum_y
    real(8) :: aa,bb, xx
    real(8) :: sigma, smallpi, sum_val
    real(8) :: g1, g2, g
    real(8) :: conv_sum, conv_sum2, conv_sum3
    real(8) :: conv_sum4, conv_sum5, conv_sum6
    real(8), allocatable :: kernel(:)


    ! 1. Calculate the sound velocity in each step.
    
    if (vsSecondOrderTOV2) then

    !print*, "Entrando a vsSecondOrderTOV2"
    idx=0
    do i=Nrtotal,1-ghost,-1

    ! At ghosts and i=1 points we set vs=vs_0
        if (i<1) then
            fluid_vs(l,i) = fluid_vs(l,1)

        else if (fluid_p(l,i)<=1e-8) then
            fluid_vs(l,i) = 0.d0

        else if (i==Nrtotal) then
            fluid_vs(l,i) = sqrt(abs(  ( fluid_p(l,i-1)   - fluid_p(l,i) ) & 
                                     / (fluid_h(l,i)*( fluid_rho(l,i-1) - fluid_rho(l,i)   ))  ))
        else 
            fluid_vs(l,i) = sqrt(abs(  ( fluid_p(l,i-1)   - fluid_p(l,i+1) ) & 
                                     / (fluid_h(l,i)*( fluid_rho(l,i-1) - fluid_rho(l,i+1) ))  ))
        end if

        if (fluid_p(l,i-1)>1e-15) then
            idx=idx+1
        end if
    end do


    !! Process of smoothing out the speed in the center

    ! Define regions
    start_10 = int(0.05d0 * idx)  ! First 5% to correct
    start_30 = int(0.40d0 * idx)  ! Previous 40% to fit
    start_60 = int(0.55d0 * idx)  ! Previous 55% boundary

    ! Step 1: Fit y = a * x^2 + b using the 30% section
    sum_x2 = 0.0d0
    sum_x4 = 0.0d0
    sum_y = 0.0d0
    sum_x2y = 0.0d0

    do l = 0,Nl - 1
    do i = start_10, start_10+start_30
        xx = real(i, 8) *dr(l)
        sum_x2 = sum_x2 + xx**2
        sum_x4 = sum_x4 + xx**4
        sum_y = sum_y + fluid_vs(l,i)
        sum_x2y = sum_x2y + fluid_vs(l,i) * xx**2
    end do
    end do

    ! Solve for a and b:
    ! a = (sum_x2y - (sum_x2 * sum_y) / n) / (sum_x4 - (sum_x2^2) / n)
    ! b = (sum_y - a * sum_x2) / n
    aa = (sum_x2y - (sum_x2 * sum_y) / (start_30 )) / &
        (sum_x4 - (sum_x2**2) / (start_30 ))


    ! Step 2: Update the first 10% with the quadratic function
    do l = 0,Nl - 1
      bb = fluid_vs(l,start_10) - aa*(dr(l)*real(start_10,8))**2
      do i = 0,start_10
        xx =  real(i, 8) *dr(l)
        fluid_vs(l,i) =  aa * xx**2 + bb
      end do
    end do

    do i=1,ghost
      fluid_vs(l,1-i)  = + fluid_vs(l,i)
    end do
    
    end if



    !! 2.
    !! Filter !!

    if (filterTOV2) then
    ! Constants
    smallpi = acos(-1.0d0)
    size_filter = ghost * 2 + 1
    sigma = 2.0d0
    half_size = size_filter / 2

    allocate(kernel(size_filter))

    ! Construct Gaussian kernel
    sum_val = 0.0d0
    do xf = -half_size, half_size
        g1 = 1.0d0 / (sigma * sqrt(2.0d0 * smallpi))
        g2 = -real(xf**2, 8) / (2.0d0 * sigma**2)
        g = g1 * exp(g2)
        kernel(xf + half_size + 1) = g
        sum_val = sum_val + g
    end do

    ! Normalize kernel
    kernel = kernel / sum_val

    ! Apply Gaussian smoothing to arrays
    do l = Nl - 1, 0, -1
        do i = 1, Nrtotal - ghost
            ! Initialize sums to zero
            conv_sum = 0.0d0
            conv_sum2 = 0.0d0
            conv_sum3 = 0.0d0
            conv_sum4 = 0.0d0
            conv_sum5 = 0.0d0
            conv_sum6 = 0.0d0

            ! Convolution over the kernel window
            do j = 1, size_filter
                !if (i + j - ghost - 1 >= 1 .and. i + j - ghost - 1 <= Nrtotal) then
                    conv_sum  = conv_sum  + fluid_rho(l, i + j - ghost - 1) * kernel(j)
                    conv_sum2 = conv_sum2 + Omega_TOV2(l, i + j - ghost - 1) * kernel(j)
                    conv_sum3 = conv_sum3 + fluid_p(l, i + j - ghost - 1) * kernel(j)
                    conv_sum4 = conv_sum4 + fluid_e(l, i + j - ghost - 1) * kernel(j)
                    conv_sum5 = conv_sum5 + fluid_vs(l, i + j - ghost - 1) * kernel(j)
                    conv_sum6 = conv_sum6 + fluid_v(l, i + j - ghost - 1) * kernel(j)
                !end if
            end do

            ! Update smoothed arrays
            fluid_rho(l, i)  = conv_sum
            Omega_TOV2(l, i) = conv_sum2
            fluid_p(l, i)    = conv_sum3
            fluid_e(l, i)    = conv_sum4
            fluid_vs(l, i)   = conv_sum5
            fluid_v(l, i)    = conv_sum6
        end do
    end do

    ! Clean up
    deallocate(kernel)

    end if






end subroutine filter_tov2