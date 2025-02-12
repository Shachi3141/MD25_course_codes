module md_simulation
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: N = 216
  real(dp), parameter :: m = 1.0_dp, sig = 1.0_dp, roh = 0.8442_dp
  real(dp), parameter :: epsilon = 1.0_dp, rc = 2.5_dp
  real(dp), parameter :: box_length = 6.348475790538236_dp, k_B  = 1.0_dp, dt = 0.001_dp,  T = 1.0_dp
  real(dp), dimension(N, 3) :: coordinates, velocities, forces, old_coordinates
contains

!--------------------------------------------------------------- POSITION INIT ----------------------------------------------
  subroutine initialize_positions(coordinates)
    implicit none
    integer :: i
    real(dp) :: lattice_spacing
    real(dp), intent(out) :: coordinates(N,3)
    integer :: n_cube, ix, iy, iz
    
    n_cube = nint(N**(1.0_dp/3.0_dp))   ! Ensure integer division
    lattice_spacing = (box_length-sig) / real((n_cube-1), dp)

    open(1,file='initial_posn.txt', status='replace')
    do i = 1, N
      ix = mod(i-1, n_cube)
      iy = mod((i-1) / n_cube, n_cube)
      iz = (i-1) / (n_cube**2)
      
      coordinates(i, 1) = (sig/2.0) + real(ix, dp) * lattice_spacing
      coordinates(i, 2) = (sig/2.0) + real(iy, dp) * lattice_spacing
      coordinates(i, 3) = (sig/2.0) + real(iz, dp) * lattice_spacing
      write(1, *) coordinates(i, :)
    end do
    print*, 'Initial posn saved.'
    close(1)
  end subroutine initialize_positions

!----------------------------------------------------- VELOCITY INITIALIZE -------------------------------------------------------

  subroutine velocity_initialize(velocities)
    implicit none
    real(dp) :: mean, std_dev
    real(dp), intent(out) :: velocities(N,3)
    integer :: i
    mean = 0.0_dp
    std_dev = sqrt(k_B * T / m)
    call random_number(velocities)
    velocities = velocities - 0.5_dp
    velocities = velocities * std_dev
    velocities = velocities - spread(sum(velocities, dim=1), dim=1, ncopies=N) / real(N, dp)
    open(2,file='initial_velocity.txt', status='replace')
    do i = 1,N
      write(2, *) velocities(i, :)
    end do
    print*, 'Initial Velocity saved.'
    close(2)
  end subroutine velocity_initialize
!------------------------------------------------------------ PBC ------------------------------------------------------------
  subroutine apply_pbc(coordinates)
    implicit none
    real(dp), intent(inout) :: coordinates(N,3)
    integer :: i, j
   
    do i = 1, N
      do j = 1, 3
        if (coordinates(i, j) > box_length) then
            coordinates(i, j) = coordinates(i, j) - box_length * floor(coordinates(i, j)/box_length)
        end if
        if (coordinates(i, j) < 0.0_dp) then
            coordinates(i, j) = coordinates(i, j) - box_length * floor(coordinates(i, j)/box_length)
            !print*, ' coor les zero!!'
        end if
      end do
    end do
  end subroutine apply_pbc
  
!----------------------------------------------------------- FORCE CALCULATION -------------------------------------------

  subroutine force_calculation(count, coordinates, forces)
    implicit none
    integer :: i, j
    integer, intent(in) :: count
    real(dp), intent(in) :: coordinates(N,3)
    real(dp), intent(out) ::  forces(N,3)
    real(dp) :: dx, dy, dz, r2, r2i, r6i, ff, total_pot, rc2, r2c, r6c,pot_cut
    
    forces = 0.0_dp
    total_pot = 0.0_dp
    rc2 = rc**2
    r2c = 1/rc2
    r6c = r2c**3
    pot_cut = 4.0*r6c*(r6c-1)
    do i = 1, N-1
      do j = i+1, N
        dx = coordinates(j,1) - coordinates(i,1)
        dy = coordinates(j,2) - coordinates(i,2)
        dz = coordinates(j,3) - coordinates(i,3)
        dx = dx - box_length * nint(dx / box_length)
        dy = dy - box_length * nint(dy / box_length)
        dz = dz - box_length * nint(dz / box_length)
        r2 = dx**2 + dy**2 + dz**2

        if (r2<(sig/10.0)**2) then
            print*, i,'and',j,'particles are too closer, dist = ',sqrt(r2),'coordinates(i,:) =',coordinates(i,:),'coordinates(j,:)=',coordinates(j,:)
            print*, 'count =', count
            stop
        end if
                
        if (r2 < rc**2) then
          r2i = 1.0_dp / r2
          r6i = r2i**3
          ff = 48.0_dp * r2i * r6i * (r6i - 0.5_dp)

          forces(j,1:3) = forces(j,1:3) + ff * (/dx, dy, dz/)
          forces(i,1:3) = forces(i,1:3) - ff * (/dx, dy, dz/)
          
          total_pot = total_pot + 4.0_dp * r6i * (r6i - 1.0_dp) - pot_cut
        end if
      end do
    end do
    
    if (mod(count, 10) == 0) then
    do i= 1,N
        write(3, *) forces(i,:)
    end do
    write(3, *) "  "
    end if
  end subroutine force_calculation
!-------------------------------------------------------------------- INTEGRATE ---------------------------------------------------
  subroutine integrate(count,coordinates,forces)
    implicit none
    integer, intent(in) :: count
    real(dp), intent(inout) :: coordinates(N,3)
    real(dp), intent(in) ::  forces(N,3)
    integer :: i,j
    real(dp) :: vi(N,3), coordinates_new(N,3), vi2tot, kin_tot, poten_tot, etot, temp
    
    coordinates_new = 2.0_dp * coordinates - old_coordinates + dt**2 * forces
    vi = (old_coordinates - coordinates) / (2.0_dp * dt)
    
    vi = vi - spread(spread(sum(vi) / real(N, dp), dim=1, ncopies=N), dim=2, ncopies=3)

    vi2tot =sum(vi**2)
    kin_tot = 0.5_dp * vi2tot / N
    poten_tot = sum(forces**2) / (2.0_dp * N)
    etot = kin_tot + poten_tot
    temp = vi2tot / (3.0_dp * real(N, dp))
    old_coordinates = coordinates
    coordinates = coordinates_new

    
    if (mod(count, 10) == 0) then
      ! write coordinstes
      do i=1,N
          write(10, *) coordinates(i,:)
      end do
      write(10, *) "   "
      
      ! write .gro file 
      write(20, *) 'Appended Frame'
      write(20, *)  N
      do i=1,N
          write(20, '(i6,1x, a,1x, a, i6, 6f22.12 )') 1, "Ar", "Ar", i, coordinates(i,:), vi(i,:)
      end do
      write(20, *)  (/box_length, box_length, box_length/)
      
      ! write energys
      write(30, *) count,kin_tot, poten_tot, etot
      
    end if
  end subroutine integrate

!------------------------------------------------------ RUN SIMULATION ------------------------------------------

  subroutine run_simulation()
    implicit none
    integer :: count, max_count
    real(dp) :: t
    max_count = 2000
    count = 0
    t = 0.0_dp
    call initialize_positions(coordinates)
    print*, "Position initializ done!"
    call velocity_initialize(velocities)
    print*, "Velocity initialize done!"
    old_coordinates = coordinates - velocities * dt

    open(3,file='force.txt', status='replace')
    open(10, file='coordinates.dat', status='replace')
    open(20, file='md.gro', status='replace')
    open(30, file='energy.dat', status='replace')
    do while (count < max_count)
      call force_calculation(count, coordinates, forces)
      call integrate(count,coordinates, forces)
      call apply_pbc(coordinates)
      t = t + dt
      count = count + 1
    end do
    close(3)
    close(10)
    close(20)
    close(30)
  end subroutine run_simulation
end module md_simulation

!--------------------------------------------------------MAIN PROGRAM---------------------------------------------------------

program main
  use md_simulation
  implicit none

  print*, 'dt = ',dt
  call run_simulation()
  call system('gnuplot energy_plot.gnu')
  print*, 'Energy plot done: check energy_plot.png'
  print*, 'Simulation done!!!  # NO ERROR AND OVERLAPES#.'
end program main
