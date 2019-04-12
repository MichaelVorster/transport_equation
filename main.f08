program main

  implicit none

  ! double precision
  integer, parameter :: dp = kind(1.d0)

  ! file containing input data
  character(len=255) :: file_input = 'input'

  ! geometry
  character(len=11) :: geometry                                  ! 'cartesian', 'cylindrical', 'spherical'
  integer :: dimensions                                          ! number of spatial dimensions

  ! time
  integer :: n_time                                              ! number of time steps
  real(dp) :: d_time                                             ! time step size

  ! spatial coordinates
  integer :: n_x1, n_x2                                          ! number of grid points in spatial coordinates - NB grid is defined to have n+1 points
  real(dp) :: x1_start, x1_end                                   ! start and end of numerical domain for first coordinate
  real(dp) :: x2_start, x2_end                                   ! start and end of numerical domain for second coordinate
  real(dp) :: d_x1, d_x2                                         ! step size of spatial coordinates 
  character(len=12) :: x1_inner_boundary_cond                    ! inner boundary condition for first coordinate
  character(len=12) :: x1_outer_boundary_cond                    ! outer boundary condition for first coordinate
  character(len=12) :: x2_inner_boundary_cond                    ! inner boundary condition for second coordinate
  character(len=12) :: x2_outer_boundary_cond                    ! outer boundary condition for second coordinate
  real(dp), dimension(:), allocatable :: x1, x2                  ! arrays containing the spatial grid points
  real(dp), dimension(6) :: v_x1                                 ! vectors incorporating boundary conditions into Thomas algorithm
  real(dp), dimension(6) :: v_x2                                 ! vectors incorporating boundary conditions into Thomas algorithm

  ! momentum
  integer :: n_momentum                                          ! number of grid points in momentum - NB grid is defined to have n+1 points
  real(dp) :: momentum_start, momentum_end                       ! start and end of numerical domain for momentum
  real(dp) :: d_ln_momentum                                      ! logarithmic step size of momentum 
  real(dp), dimension(:), allocatable :: momentum                ! array containing the momentum grid points

  ! particle distribution function
  real(dp), dimension(:,:,:), allocatable :: f                   ! array containing solutions of transport equation

  ! source functions at boundaries
  real(dp), dimension(:,:,:), allocatable :: source              ! array containing sources

  ! coefficients of the partial differentials in the transport equation
  real(dp), dimension(:,:,:), allocatable :: coef_x1x1, coef_x1  ! coefficients of the second and first order partial derivatives, respectively
  real(dp), dimension(:,:,:), allocatable :: coef_x2x2, coef_x2  ! coefficients of the second and first order partial derivatives, respectively
!  real(dp), dimension(:,:,:), allocatable :: coef_mom, coef_f    ! coefficients of the first order partial derivatives w.r.t. momentum and the coefficient of the solution, respectively
  real(dp), dimension(:,:), allocatable :: coef_ad, coef_syn     ! coefficients for adiabatic and synchrotron loss rate (without momentum dependence)

  integer :: i, j, k, time_step
  real(dp) :: plot_fac
  logical :: lexist

  namelist /input/ &
    geometry, &
    dimensions, &
    n_time, &
    d_time, &
    n_x1, &
    n_x2, &
    x1_start, &
    x1_end, &
    x2_start, &
    x2_end, &
    x1_inner_boundary_cond, &
    x1_outer_boundary_cond, &
    x2_inner_boundary_cond, &
    x2_outer_boundary_cond, &
    n_momentum, &
    momentum_start, &
    momentum_end

  ! check existence of file 'input' and read its parameters
  inquire(file=file_input, exist=lexist)
  if (lexist) then
    open(unit=10, file=file_input, form='formatted', status='old', action='read')
    write(*,"('Opening file ',a,'...')", advance='yes') trim(file_input)
  else
    write(*,"('Error: file ',a,' does not exist. Stopping.')") trim(file_input)
    stop
  endif

  write(*,"('Reading namelist /input/...')")
  read(unit=10, nml=input)
    if (dimensions .eq. 1) n_x2 = 0
  write(*, nml=input)
  close(unit=10)

  ! allocate memory for arrays and initialise to zero
  allocate(x1(0:n_x1))
  allocate(x2(0:n_x2))
  allocate(momentum(0:n_momentum))
  allocate(f(0:n_x1, 0:n_x2, 0:n_momentum))
  allocate(coef_x1x1(0:n_x1, 0:n_x2, 0:n_momentum))
  allocate(coef_x1(0:n_x1, 0:n_x2, 0:n_momentum))
  allocate(coef_x2x2(0:n_x1, 0:n_x2, 0:n_momentum))
  allocate(coef_x2(0:n_x1, 0:n_x2, 0:n_momentum))
!  allocate(coef_mom(0:n_x1, 0:n_x2, 0:n_momentum))
!  allocate(coef_f(0:n_x1, 0:n_x2, 0:n_momentum))
  allocate(coef_ad(0:n_x1, 0:n_x2))
  allocate(coef_syn(0:n_x1, 0:n_x2))
  allocate(source(0:n_x1, 0:n_x2, 0:n_momentum))

  ! initialise arrays to zero
  x1 = 0.
  x2 = 0.
  momentum = 0.
  f = 0.  
  coef_x1x1 = 0.
  coef_x1 = 0.
  coef_x2x2 = 0.
  coef_x2 = 0.
!  coef_mom = 0.
!  coef_f = 0.
  coef_ad = 0.
  coef_syn = 0.
  source = 0.

  ! calculate grid step sizes (momentum in logarithmic steps)
  d_x1 = (x1_end - x1_start) / n_x1
  d_x2 = (x2_end - x2_start) / n_x2
    if (dimensions .eq. 1) d_x2 = 1.
  d_ln_momentum = log(momentum_end/momentum_start) / n_momentum

    ! populate grid
  do i = 0, n_x1
    x1(i) = x1_start + d_x1*i
  enddo  
  do j = 0, n_x2 
    x2(j) = x2_start + d_x2*j
  enddo
  do k = 0, n_momentum 
    momentum(k) = momentum_start*exp(d_ln_momentum*k)
  enddo

  print*, 'Step size for first coordinate', d_x1
  print*, 'Step size for second coordinate', d_x2
  print*, 'Logarithmic step size for momentum', d_ln_momentum
  print*, 'Total time', d_time*n_time


  call coefficients( &
    geometry, &
    n_x1, &
    n_x2, &
    n_momentum, &
    x1, &
    x2, &
    momentum, &
    coef_x1x1, &
    coef_x1, &
    coef_x2x2, &
    coef_x2, &
!    coef_mom, &
!    coef_f &
    coef_ad, &
    coef_syn, &
    dimensions &
  )

  call initial_condition(n_x1, n_x2, n_momentum, x1, x2, momentum, f)

  call source_function( &
    x1_inner_boundary_cond, &
    x1_outer_boundary_cond, &
    x2_inner_boundary_cond, &
    x2_outer_boundary_cond, &
    n_x1, &
    n_x2, &
    n_momentum, &
    x1, &
    x2, &
    momentum, &
    source, &
    f &
  )

  call boundary_conditions_thomas( &
    x1_inner_boundary_cond, &
    x1_outer_boundary_cond, &
    v_x1 &
  )
    if (dimensions .eq. 2) then
      call boundary_conditions_thomas( &
        x2_inner_boundary_cond, &
        x2_outer_boundary_cond, &
        v_x2 &
      )
    else
      v_x2 = 0.
    endif


  ! ! start of numerical solution
  do time_step = 0, n_time
!$OMP PARALLEL DEFAULT(SHARED) 

print*, time_step, 'a', f(0,0,n_momentum - 10)
print*, time_step, 'a', source(0,0,n_momentum - 10)

!$OMP DO 
    ! solve along first spatial coordinate
    do k = 0, n_momentum
      do j = 0, 0
        call crank_nicolson( &
          f(:,j,k), &
          n_x1, &
          d_time, &
          d_x1, &
          v_x1, &
          coef_x1x1(:,j,k), &
          coef_x1(:,j,k), &
          source(:,j,k) &
        )
      enddo
    enddo
!$OMP END DO NOWAIT

print*, time_step, 'b',  f(0,0,n_momentum - 10)
print*, time_step, 'b',  source(0,0,n_momentum - 10)

!$OMP DO 
    ! solve momentum equation analytically
    do i = 1, n_x1-1
      do j = 0, 0
        call analytic_momentum_solution( &
          f(i,j,:), &
          d_time, &
          d_ln_momentum, &
          n_momentum, &
          momentum, &
          coef_ad(i,j), &
          coef_syn(i,j) &
        )
      enddo
    enddo
!$OMP END DO NOWAIT

print*, time_step, 'c',  f(0,0,n_momentum - 10)
print*, time_step, 'c',  source(0,0,n_momentum - 10)

!$OMP END PARALLEL
  enddo ! time_step

  open (unit=12,file='results.txt',status='unknown',form='formatted')
    do k = 0, n_momentum
      plot_fac = momentum(k)**4
      write(12, *) momentum(k), plot_fac*f(0,0,k), plot_fac*f(20,0,k), plot_fac*f(30,0,k), plot_fac*f(50,0,k), plot_fac*f(80,0,k)
    enddo
  close(unit=12) 

  ! free memory
  if (allocated(x1)) deallocate(x1)
  if (allocated(x2)) deallocate(x2)
  if (allocated(momentum)) deallocate(momentum)
  if (allocated(f)) deallocate(f)
  if (allocated(coef_x1x1)) deallocate(coef_x1x1)
  if (allocated(coef_x1)) deallocate(coef_x1)
  if (allocated(coef_x2x2)) deallocate(coef_x2x2)
  if (allocated(coef_x2)) deallocate(coef_x2)
!  if (allocated(coef_mom)) deallocate(coef_mom)
!  if (allocated(coef_f)) deallocate(coef_f)
  if (allocated(coef_ad)) deallocate(coef_ad)
  if (allocated(coef_syn)) deallocate(coef_syn)
  if (allocated(source)) deallocate(source)

end program main


! specify the initial condition over the whole computaional grid
subroutine initial_condition(n_x1, n_x2, n_momentum, x1, x2, momentum, f)
  implicit none

  integer, parameter :: dp = kind(1.d0)

  integer, intent(in) :: n_x1, n_x2, n_momentum                          ! number of grid points for spatial and momentum dimensions
  real(dp), dimension(0:n_x1), intent(in) :: x1                          ! array of grid points for first spatial dimension
  real(dp), dimension(0:n_x2), intent(in) :: x2                          ! array of grid points for first spatial dimension
  real(dp), dimension(0:n_momentum), intent(in) :: momentum              ! array containing momentum grid points
  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: f  ! array containing solutions to transport equation (f is the distribution function)

  integer :: i, j, k

  do k = 1, n_momentum-1
    do j = 0, 0
      do i = 1, n_x1-1
        f(i,j,k) = 0.
      enddo  
    enddo
  enddo

  return
end subroutine initial_condition


! specify any source functions, including at the boundaries
subroutine source_function( &
  x1_inner_boundary_cond, &
  x1_outer_boundary_cond, &
  x2_inner_boundary_cond, &
  x2_outer_boundary_cond, &
  n_x1, &
  n_x2, &
  n_momentum, &
  x1, &
  x2, &
  momentum, &
  source, &
  f &
)
  implicit none

  integer, parameter :: dp = kind(1.d0)

  character(len=12), intent(in) :: x1_inner_boundary_cond, x1_outer_boundary_cond  ! boundary conditions of first spatial  coordinate              
  character(len=12), intent(in) :: x2_inner_boundary_cond, x2_outer_boundary_cond  ! boundary conditions of second spatial  coordinate
  integer, intent(in) :: n_x1, n_x2, n_momentum                                    ! number of grid points
  real(dp), dimension(0:n_x1), intent(in) :: x1                                    ! array of grid points for first spatial dimension
  real(dp), dimension(0:n_x2), intent(in) :: x2                                    ! array of grid points for first spatial dimension
  real(dp), dimension(0:n_momentum), intent(in) :: momentum                        ! array of grid points for momentum
  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: source       ! array containing sources
  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: f            ! array containing solutions of transport equation

  real(dp) :: LIS

  integer :: i, j, k

  if (x1_inner_boundary_cond .eq. 'source') then
    do k = 0, n_momentum
      do j = 0, 0
        source(0,j,k) = 1.0e0*momentum(k)**(-4.)
        f(0,j,k) = source(0,j,k)
      enddo
    enddo    
  endif

  if (x1_outer_boundary_cond .eq. 'source') then
    do k = 0, n_momentum
      do j = 0, 0
        source(n_x1,j,k) = LIS(momentum(k))/momentum(k)/momentum(k)
      enddo
    enddo    
  endif

  if (x2_inner_boundary_cond .eq. 'source') then
    do k = 0, n_momentum
      do i = 0, n_x1
        source(i,0,k) = 0.
      enddo
    enddo    
  endif

  if (x2_outer_boundary_cond .eq. 'source') then
    do k = 0, n_momentum
      do i = 0, n_x1
        source(i,n_x2,k) = 0.
      enddo
    enddo    
  endif

  return
end subroutine source_function


! set vectors that are needed to incorporate boundary conditions into the Thomas algorithm
subroutine boundary_conditions_thomas(inner_condition, outer_condition, v_boundary)
  implicit none

  integer, parameter :: dp = kind(1.d0)

  character(len=12), intent(in) :: inner_condition, outer_condition  ! boundary conditions
  real(dp), dimension(6), intent(out) :: v_boundary                  ! array needed to incorporate boundary conditions in Thomas algorithm 

  v_boundary = 0.

  if (inner_condition .eq. 'escape') then
    v_boundary(1:3) = 0.
  elseif (inner_condition .eq. 'reflective') then
    v_boundary(2) = 1.
    print*, 'Current implementation of this boundary condition makes scheme only 1st order accurate'
  elseif (inner_condition .eq. 'source') then
    v_boundary(3) = 1.
  else
    print*, 'Inner boundary condition set incorrectly'
  endif

  if (outer_condition .eq. 'escape') then
    v_boundary(4:6) = 0.
  elseif (outer_condition .eq. 'reflective') then
    v_boundary(5) = 1.
    print*, 'Current implementation of this boundary condition makes scheme only 1st order accurate'
  elseif (outer_condition .eq. 'source') then
    v_boundary(6) = 1.
  else
    print*, 'Outer boundary condition set incorrectly'
  endif

  return
end subroutine boundary_conditions_thomas


! calculate the coefficients of the transport equation
subroutine coefficients( &
  geometry, &
  n_x1, &
  n_x2, &
  n_momentum, &
  x1, &
  x2, &
  momentum, &
  coef_x1x1, &
  coef_x1, &
  coef_x2x2, &
  coef_x2, &
  coef_ad, &
  coef_syn, &
  dimensions &
)
  implicit none

  integer, parameter :: dp = kind(1.d0)

  character(len=11), intent(in) :: geometry                                     ! 'cartesian', 'cylindrical', 'spherical'
  integer, intent(in) :: dimensions                                             ! number of spatial dimensions 
  integer, intent(in) :: n_x1, n_x2, n_momentum                                 ! number of grid points
  real(dp), dimension(0:n_x1), intent(in) :: x1                                 ! array of grid points for first spatial dimension
  real(dp), dimension(0:n_x2), intent(in) :: x2                                 ! array of grid points for first spatial dimension
  real(dp), dimension(0:n_momentum), intent(in) :: momentum                     ! array of grid points for momentum

  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: coef_x1x1 ! coefficient of the second order partial derivative, first spatial coordinate
  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: coef_x1   ! coefficient of the first order partial derivative, first spatial coordinate
  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: coef_x2x2 ! coefficient of the second order partial derivative, second spatial coordinate
  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: coef_x2   ! coefficient of the first order partial derivatives, second spatial coordinate
!  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: coef_mom  ! coefficient of the first order partial derivatives w.r.t. momentum
!  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum), intent(inout) :: coef_f    ! coefficient of the solution
  real(dp), dimension(0:n_x1, 0:n_x2), intent(inout) :: coef_ad                 ! coefficient for adiabatic loss rate (without momentum dependence)
  real(dp), dimension(0:n_x1, 0:n_x2), intent(inout) :: coef_syn                ! coefficient for synchrotron loss rate (without momentum dependence)

  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum) :: kappa_x1x1               ! diffusion coefficient
  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum) :: kappa_x2x2               ! diffusion coefficient
  real(dp), dimension(0:n_x1, 0:n_x2) :: V_x1                                   ! convection velocity  
  real(dp), dimension(0:n_x1, 0:n_x2) :: V_x2                                   ! convection velocity
  real(dp), dimension(0:n_x1, 0:n_x2) :: B_field                                ! magnetic field

  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum) :: dkappa_x1x1_dx1          ! derivative of diffusion coefficient with respect to x1
  real(dp), dimension(0:n_x1, 0:n_x2, 0:n_momentum) :: dkappa_x2x2_dx2          ! derivative of diffusion coefficient with respect to x2
  real(dp), dimension(0:n_x1, 0:n_x2) :: dV_x1_dx1                              ! derivative of velocity componennt with respect to x1   
  real(dp), dimension(0:n_x1, 0:n_x2) :: dV_x2_dx2                              ! derivative of velocity component with respect to x1

  real(dp) :: kappa0 = 5.                                                  ! initial diffusion coefficient 
  real(dp) :: V0 = 1.                                                          ! initial convection velocity
  real(dp) :: alpha = 0.                                                       ! parameter that controls radial velocity profile
  real(dp) :: sync_coef = 0.05*0
  real(dp) :: E0 = 5.11e-7
  real(dp) :: beta
         
  integer :: i, j, k

  do k = 0, n_momentum
    do j = 0, 0
      do i = 0, n_x1
        
        ! x1 = radius, x2 = polar angle 
        if (geometry .eq. 'spherical') then
          ! for heliosphere
          !beta = momentum(k)/sqrt(momentum(k)**2 + E0**2)
          !kappa_x1x1(i,j,k) = kappa0*momentum(k)*beta
          !V_x1(i,j) = V0*(1.0-exp(-13.862*x1(i)))
          !dV_x1_dx1(i,j) = (V0-V_x1(i, j))*13.862

          kappa_x1x1(i,j,k) = kappa0*momentum(k)
          dkappa_x1x1_dx1(i,j,k) = 0.
          kappa_x2x2(i,j,k) = 0.
          dkappa_x2x2_dx2(i,j,k) = 0.
          V_x1(i,j) = V0*(x1(0)/x1(i))**(alpha)
          dV_x1_dx1(i,j) = -alpha*V_x1(i,j)/x1(i)
          V_x2(i,j) = 0.
          dV_x2_dx2(i,j) = 0.
          B_field(i,j) = 1.0

          if (dimensions .eq. 1) then
            kappa_x2x2(i,j,k) = 0.
            dkappa_x2x2_dx2(i,j,k) = 0.
            V_x2(i,j) = 0.
            dV_x2_dx2(i,j) = 0.
          endif  

          coef_x1x1(i,j,k) = kappa_x1x1(i,j,k)
          coef_x2x2(i,j,k) = kappa_x2x2(i,j,k)/(x1(i)**2)
          coef_x1(i,j,k) = &
            2.*kappa_x1x1(i,j,k)/x1(i) &
            + dkappa_x1x1_dx1(i,j,k) &
            - V_x1(i,j)
          coef_x2(i,j,k) = &
            kappa_x2x2(i,j,k)/(x1(i)**2 + tan(x2(j))) &
            + dkappa_x2x2_dx2(i,j,k)/(x1(i)**2) &
            - V_x2(i,j)/x1(i)
          ! coef_mom(i,j,k) = &
          !   1./3.*(2.*V_x1(i,j)/x1(i) + dV_x1_dx1(i,j) + V_x2(i,j)/(x1(i)*tan(x2(j))) + dV_x2_dx2(i,j)/x1(i)) &
          !   + sync_coef*B_field(i,j)**(2.)*momentum(k)
          ! coef_f(i,j,k) = 4.*sync_coef*B_field(i,j)**(2.)*momentum(k)
          coef_ad(i,j) = &
            1./3.*(2.*V_x1(i,j)/x1(i) + dV_x1_dx1(i,j) + V_x2(i,j)/(x1(i)*tan(x2(j))) + dV_x2_dx2(i,j)/x1(i))
          coef_syn(i,j) = sync_coef*B_field(i,j)**(2.)
        endif

      enddo  
    enddo
  enddo

  return
end subroutine coefficients


! solution of a one-dimensional PDE using the Crank_Nicolson method
subroutine crank_nicolson( &
  u, &
  n, &
  d_t, &
  d_x, &
  v_boundary, &
  coef_xx, &
  coef_x, &
  source &
)
  implicit none

  integer, parameter :: dp = kind(1.d0)

  real(dp), intent(inout), dimension(0:n) :: u             ! solution of the transport equation
  integer, intent(in) :: n                                 ! number of grid points along spatial direction
  real(dp), intent(in) :: d_t, d_x                         ! time and spatial step sizes, respectively 
  real(dp), intent(inout), dimension(6) :: v_boundary         ! vector incorporating boundary conditions into Thomas algorithm
  real(dp), intent(in), dimension(0:n) :: coef_xx, coef_x  ! coefficients of the second and first order partial derivatives, respectively 
  real(dp), intent(in), dimension(0:n) :: source           ! array specifying any sources, including at the boundary

  real(dp) coef_i_min_1, coef_i, coef_i_plus_1             ! coefficients in the Crank-Nicolson method, for grid points at i-1, i, i+1, respectively 
  real(dp) RS_known_values                                 ! the known solutions at the previous time step (the right-hand side of the Crank-Nicolson equation)
  real(dp), dimension(0:n) :: Z, d                         ! variables used in Thomas algorithm
  real(dp), dimension(-n:n) :: c                           ! variables used in Thomas algorithm

  integer :: i

  Z = 0.
  d = 0.
  c = 0.

  v_boundary(1) = 1./(1.+1.0*d_x/coef_xx(1))
  v_boundary(3) = d_x/coef_xx(1)/(1.+1.0*d_x/coef_xx(1))


  c(0) = -v_boundary(1)
  c(-1) = -v_boundary(2)
  d(0) = v_boundary(3)*source(0)

  u(0) = v_boundary(1)*u(1) + v_boundary(2)*u(2) + v_boundary(3)*source(0)
  print*, source(0), u(0)
  u(n) = v_boundary(4)*u(n-1) + v_boundary(5)*u(n-2) + v_boundary(6)*source(n)

  do i = 1, n-1
    coef_i_min_1 = 1./(4.*d_x)*coef_x(i) - 1./(2.*d_x*d_x)*coef_xx(i)
    coef_i = 1./(d_t) + 1./(d_x*d_x)*coef_xx(i)
    coef_i_plus_1 = -1./(4.*d_x)*coef_x(i) - 1./(2.*d_x*d_x)*coef_xx(i)

    RS_known_values = &
        - coef_i_min_1*u(i-1) &
        + (1./(d_t) - 1./(d_x*d_x)*coef_xx(i))*u(i) &
        - coef_i_plus_1*u(i+1)

    Z(i) = coef_i - c(i-1)*coef_i_min_1
    c(i) = (coef_i_plus_1 - c(-i)*coef_i_min_1) / Z(i)
    d(i) = (RS_known_values - d(i-1)*coef_i_min_1) / Z(i)
  enddo

  u(n-1) = &
    (d(n-1) - c(n-1)*(v_boundary(6)*source(n) + v_boundary(5)*d(n-2))) &
      / (1. + c(n-1)*(v_boundary(4) - v_boundary(5)*c(n-2)))
  do i = n-2, 1, -1
    u(i) = d(i) - c(i)*u(i+1)
  enddo

  return
end subroutine crank_nicolson


! the energy equation is solved analytically using the method of characteristics
subroutine analytic_momentum_solution( &
  u, &
  d_t, &
  d_ln_momentum, &
  n, &
  momentum, &
  coef_ad, &
  coef_syn &
)
  implicit none

  integer, parameter :: dp = kind(1.d0)

  real(dp), intent(inout), dimension(0:n) :: u  ! the distribution function / solution of the transport equation
  real(dp), intent(in) :: d_t                   ! fractional time step
  real(dp), intent(in) :: d_ln_momentum         ! momentum step size
  integer, intent(in) :: n                      ! number of momentum grid points 
  real(dp), intent(in) :: momentum(0:n)         ! array of momentum values
  real(dp), intent(in) :: coef_ad               ! coefficient of adaibatic enregy loss rate (without momentum dependence)
  real(dp), intent(in) :: coef_syn              ! coefficient of synchrotron loss rate (without momentum dependence)

  real(dp), dimension(0:n) :: u_known           ! the known solution at the previous time
  real(dp) x_0, y_0                             ! variables used in the analytic calculation
  real(dp) momentum_0, u_0                      ! variables used in the analytic calculation
  integer index

  integer :: k


  u_known = u
  u(0) = 0.
  u(n) = 0.

  do k = 1, n-1
    x_0 = coef_ad*momentum(k)
    y_0 = (coef_ad+coef_syn*momentum(k))*exp(-coef_ad*d_t)-coef_syn*momentum(k)
    momentum_0 = x_0 / y_0      
    index = int(log(momentum_0/momentum(0))/d_ln_momentum)

    if ((momentum_0 .le. 0.) .or. (y_0 .lt. 0.) .or. (index .ge. n)) then
      u(k) = 1.0e-30
    else
      u_0 = u_known(index) + (u_known(index+1) - &
        u_known(index))*(log(momentum_0)-log(momentum(index)))/d_ln_momentum

        if (u_0 .le. 0.) then
          u(k) = 1.0e-30
        else
          u(k) = u_0*exp(4.*coef_syn*momentum(k)*d_t)
        endif

    endif
  enddo

  return
end subroutine analytic_momentum_solution


! The local Interstellar cosmic ray spectrum
real*8 function LIS(momentum)
  implicit none

  integer, parameter :: dp = kind(1.d0)

  real(dp) :: beta, fnp, momentum, T
  real(dp) :: E0, AZ

  AZ = 1.0
  E0 = 0.938

  beta = momentum/sqrt(momentum*momentum + E0*E0*AZ*AZ)
  fnp = beta*momentum
  T = momentum/beta/AZ-E0
  LIS = 21.1*exp(log(T)*(-2.8))/(1. + 5.85*exp(log(T)*(-1.22)) + 1.18*exp(log(T)*(-2.54)))

  return
end function LIS
