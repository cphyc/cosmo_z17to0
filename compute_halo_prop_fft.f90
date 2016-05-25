program compute_halo_prop
  use io
  use misc
  use types
  use flap, only : command_line_interface
  use compute
  use convolution
  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: param_output_path
  integer            :: param_output_number
  integer, parameter :: NPARTICLE_TO_PROBE_HALO = 50, NCPU_PER_HALO = 10
  real(dp)           :: param_min_m, param_max_m
  integer            :: param_verbosity

  type(command_line_interface)               :: cli
  !-------------------------------------
  ! Brick data
  !-------------------------------------
  real(dp), dimension(:, :), allocatable     :: posDM, LDM
  integer, dimension(:), allocatable         :: idDM
  real(dp), dimension(:), allocatable        :: rvirDM, mDM, mvirDM, TvirDM, hlevel
  type(INFOS_T)                              :: infos
  type(MEMBERS_T), dimension(:), allocatable :: members
  integer                                    :: nbodies, nb_of_halos, nb_of_subhalos, nDM
  real(sp)                                   :: aexp_tmp, age_univ
  !-------------------------------------
  ! Particle data
  !-------------------------------------
  integer                                    :: ndim, nparts
  real(dp), dimension(:, :), allocatable     :: pos, vel
  integer                                    :: nstar
  integer,      dimension(:), allocatable    :: ids, order
  real(dp), dimension(:), allocatable        :: m
  real(sp), dimension(:), allocatable        :: birth_date
  !-------------------------------------
  ! Halo to cpu
  !-------------------------------------
  integer, dimension(:), allocatable         :: cpu_list
  !-------------------------------------
  ! Halo properties
  !-------------------------------------
  real(dp), dimension(:,:,:), allocatable    :: density
  real(dp), dimension(:), allocatable        :: m_in_halo, m_in_box
  real(dp), dimension(:, :), allocatable     :: pos_in_halo,&
       & pos_around_halo, grid
  integer, dimension(:), allocatable         :: ids_in_box, ids_around_halo
  real(dp), dimension(:, :), allocatable     :: pos_in_box
  real(dp)                                   :: mtot
  real(dp), dimension(3)                     :: pos_mean, pos_std
  integer                                    :: halo_i, halo_counter, part_i
  integer                                    :: part_counter, parts_in_region, local_part_counter
  logical, dimension(:), allocatable         :: halo_found_mask
  type(CONV_T)                               :: conv
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                    :: i, counter, cpu_i, x, y, z, dim
  integer                                    :: tmp_int, index, tmp_int2
  character(len=200)                         :: tmp_char
  real(dp), allocatable, dimension(:,:)      :: tmp_dblarr
  real(dp), dimension(3)                     :: X0, X1
  real(dp)                                   :: margin
  integer                                    :: param_nstep, max_nparts, step_i

  !----------------------------------------
  ! FFT variables
  !----------------------------------------
  real(dp), allocatable :: gaussian(:, :, :), conv_dens(:, :, :), edges(:, :)
  real(dp) :: sigma
  integer :: param_nbin


  !-------------------------------------
  ! random
  !-------------------------------------
  call random_seed()

  !-------------------------------------
  ! Read parameters
  !-------------------------------------
  call cli%init(progname='compute_halo_prop_fft', &
       description='Compute the FFT of each halo given and its surrounding')
  call cli%add(switch='--brick', help='Path for brick file that contains the halos', &
       act='store', def='/data52/Horizon-AGN/TREE_DM_celldx2kpc_SC0.9r/tree_bricks782')
  call cli%add(switch='--info-file', help='Information file about simulation', &
       act='store', def='/data52/Horizon-AGN/OUTPUT_DIR/output_00782/info_00782.txt')
  call cli%add(switch='--margin', help='margin', act='store', def='0.15')
  call cli%add(switch='--max-mass', switch_ab='-maxm', help='Maximum mass', &
       act='store',def='1e99')
  call cli%add(switch='--min-mass', switch_ab='-minm', help='Minimum mass', &
       act='store',def='0')
  call cli%add(switch='--nstep', act='store', def='4', &
       help='Number of step in which to compute the data')
  call cli%add(switch='--nbin', act='store', def='32', &
       help='Number of bins to use to estimate density')
  call cli%add(switch='--output', switch_ab='-o', help='Name of the output file', &
       act='store', def='data.out')
  call cli%add(switch='--output-number', help='Number of the output(s)',&
       act='store', def='2', nargs='+')
  call cli%add(switch='--output-path', help='Path of the simulation output',&
       act='store', def='/data52/Horizon-AGN/OUTPUT_DIR')
  call cli%add(switch='--verbose', help='Verbosity', act='store', def='0')

  call cli%get(switch='--min-mass', val=param_min_m)
  call cli%get(switch='--max-mass', val=param_max_m)
  call cli%get(switch='--output-path', val=param_output_path)
  call cli%get(switch='--output-number', val=param_output_number)
  call cli%get(switch='--verbose', val=param_verbosity)
  call cli%get(switch='--margin', val=margin)
  call cli%get(switch='--nstep', val=param_nstep)
  call cli%get(switch='--nbin', val=param_nbin)

  !-------------------------------------
  ! Read brick file
  !-------------------------------------
  print*, ''
  print*, 'Reading info file…'
  call read_info_headers(param_output_path, param_output_number, infos)
  print*, '   …red'

  print*, ''
  call cli%get(switch='--brick', val=tmp_char)
  print*, 'Reading brick file', trim(tmp_char)
  call read_brick_header(tmp_char, infos, nbodies, aexp_tmp, age_univ,&
       nb_of_halos, nb_of_subhalos)
  nDM = nb_of_halos + nb_of_subhalos
  allocate(idDM(nDM))
  allocate(posDM(infos%ndim, nDM))
  allocate(rvirDM(nDM), mDM(nDM), mvirDM(nDM), TvirDM(nDM), hlevel(nDM))
  allocate(LDM(infos%ndim, nDM))
  allocate(members(nDM))
  allocate(cpu_list(infos%ncpu))
  call read_brick_data(nDM, infos, .true., &
       & mDM, posDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LDM, idDM, members)
  print*, '    …red!'
  deallocate(rvirDM, mvirDM, TvirDM, hlevel, LDM)

  call cli%get(switch='--output', val=tmp_char)
  print*, ''
  print*, 'Writing output in ', trim(tmp_char)
  open(unit=10, file=trim(tmp_char))
  open(unit=123, file='missing_parts')
  write(10, '(a9, 20a13)') 'id', 'mass', 'x', 'y', 'z', 'xx', 'xy', 'xz', 'yy', 'yz', 'zz'
  halo_counter = 1

  max_nparts = 2**30


  !----------------------------------------
  ! Allocate data
  !----------------------------------------
  allocate(gaussian(param_nbin, param_nbin, param_nbin), &
       conv_dens(param_nbin, param_nbin, param_nbin), &
       density(param_nbin, param_nbin, param_nbin))
  allocate(edges(infos%ndim, param_nbin+1))
  allocate (halo_found_mask(nDM))
  halo_found_mask = .false.


  !-------------------------------------
  ! Split the box into cells of equal sizes
  !-------------------------------------
  do x = 0, param_nstep-1
     do y = 0, param_nstep-1
        do z = 0, param_nstep-1
           X0 = (/ x*1d0/param_nstep, y*1d0/param_nstep, z*1d0/param_nstep /)
           X1 = (/ (x+1)*1d0/param_nstep, (y+1)*1d0/param_nstep, (z+1)*1d0/param_nstep /) + margin
           !----------------------------------------
           ! Get the list of all cpus
           !----------------------------------------
           call get_cpu_list(X0, X1, infos%levelmax, infos%bound_key, &
                cpu_list, infos%ncpu, infos%ndim)

           ! count the number of cpus
           tmp_int = 0
           do cpu_i = 1, infos%ncpu
              if (cpu_list(cpu_i) > 0) then
                 tmp_int = tmp_int + 1
              end if
           end do

           step_i = x*param_nstep**2 + y*param_nstep + z + 1
           write(*, '(a, i4,    a,   i4,      1x,  a, i5,1x, a, 2(3f7.4,1x))') &
                'step', step_i, '/', param_nstep**3, &
                'ncpu=', tmp_int, ', box', X0, X1
           allocate(m_in_box(max_nparts), pos_in_box(infos%ndim, max_nparts), &
                ids_in_box(max_nparts))
           ids_in_box = 0
           !-------------------------------------
           ! Loop over all cpus listed, read the particle and save them
           !-------------------------------------
           part_counter = 0
           tmp_int2 = 0
           !$OMP PARALLEL DO DEFAULT(none) SHARED(cpu_list, param_output_path) &
           !$OMP SHARED(param_output_number, part_counter, X0, X1, pos_in_box) &
           !$OMP SHARED(m_in_box, ids_in_box, infos, param_verbosity, tmp_int2)&
           !$OMP PRIVATE(nstar, pos, vel, m, ids, birth_date, ndim, nparts, i) &
           !$OMP PRIVATE(tmp_int, local_part_counter)                          &
           !$OMP SCHEDULE(guided, 1)
           do cpu_i = 1, infos%ncpu
              if (cpu_list(cpu_i) == 0) then
                 cycle
              end if

              !$OMP ATOMIC
              tmp_int2 = tmp_int2 + 1

              call read_particle(param_output_path, param_output_number, cpu_list(cpu_i), &
                   nstar, pos, vel, m, ids, birth_date, ndim, nparts)

              if (param_verbosity >= 1) then
                 write(*, '(a,          i5,    a,   i5,      a,    i5, a, i10, a, i10, a)') &
                      'Reading cpu', tmp_int2, '/', infos%ncpu, ' (cpu=', cpu_list(cpu_i), &
                      ', nparts=', nparts, ', loaded=', part_counter,')'
              end if
              !$OMP CRITICAL
              do i = 1, nparts
                 if ( pos(1, i) >= X0(1) .and. pos(1, i) < X1(1) .and. &
                      pos(2, i) >= X0(2) .and. pos(2, i) < X1(2) .and. &
                      pos(3, i) >= X0(3) .and. pos(3, i) < X1(3) .and. &
                      ids(i) > 0 & ! only keep DM particles
                      ) then

                    part_counter = part_counter + 1
                    pos_in_box(:, part_counter) = pos(:, i)
                    m_in_box(part_counter)      = m(i)
                    ids_in_box(part_counter)    = ids(i)
                 end if
              end do
              !$OMP END CRITICAL

           end do
           !$OMP END PARALLEL DO

           allocate(order(part_counter))
           print*, 'Sorting array (len=', part_counter, ')…'
           ! Sort the ids
           call quick_sort(ids_in_box(1:part_counter), order)
           print*, '   … sorted'

           !-------------------------------------
           ! Once the particle read, for all complete halos
           ! compute the properties
           !-------------------------------------
           !$OMP PARALLEL DO DEFAULT(none) &
           !$OMP SHARED(mDM, param_min_m, param_max_m, X0, X1, posDM, halo_found_mask) &
           !$OMP SHARED(infos, param_nbin, ids_in_box, members, part_counter, m_in_box) &
           !$OMP SHARED(nDM, max_nparts, pos_in_box, param_verbosity, order) &
           !$OMP PRIVATE(pos_in_halo, m_in_halo, pos_around_halo, ids_around_halo) &
           !$OMP PRIVATE(mtot, tmp_dblarr, parts_in_region, grid, edges, conv, conv_dens) &
           !$OMP PRIVATE(sigma, density, gaussian, index, pos_mean, pos_std, counter) &
           !$OMP SCHEDULE(guided, 1)
           do halo_i = 1, nDM
              ! filter halos outside mass range and not in box
              if ( mDM(halo_i) > param_min_m .and. mDM(halo_i) < param_max_m .and. &
                   posDM(1, halo_i) >= X0(1) .and. posDM(1, halo_i) < X1(1) .and. &
                   posDM(2, halo_i) >= X0(2) .and. posDM(2, halo_i) < X1(2) .and. &
                   posDM(3, halo_i) >= X0(3) .and. posDM(3, halo_i) < X1(3) .and. &
                   (.not. halo_found_mask(halo_i)) ) then

                 ! allocate data for the halo
                 allocate(&
                      pos_in_halo(infos%ndim, members(halo_i)%parts), &
                      m_in_halo(members(halo_i)%parts),&
                      pos_around_halo(infos%ndim, max_nparts),&
                      ids_around_halo(max_nparts))

                 ! count the number of particles in the region that are in the halo
                 counter = 0
                 do part_i = 1, members(halo_i)%parts
                    index = indexOf(members(halo_i)%ids(part_i), ids_in_box(1:part_counter))
                    if (index > 0) then
                       pos_in_halo(:, part_i) = pos_in_box(:, order(index))
                       m_in_halo(part_i)      = m_in_box(order(index))
                       counter = counter + 1
                    end if
                 end do

                 ! cheers (once) if complete
                 if ((counter == members(halo_i)%parts) .and. & ! cheers if complete …
                      (.not. halo_found_mask(halo_i))) then     ! … only once
                    if (param_verbosity >= 2) then
                       print*, halo_i, 'is complete'
                    end if
                    !$OMP CRITICAL
                    halo_found_mask(halo_i) = .true.
                    !$OMP END CRITICAL

                    ! correct the positions and compute data about it
                    call correct_positions(pos_in_halo)
                    call compute_mean(pos_in_halo, pos_mean)
                    do dim = 1, infos%ndim
                       call stddev(pos_in_halo(dim, :) - pos_mean(dim), pos_std(dim), &
                            members(halo_i)%parts)
                    end do

                    mtot = sum(m_in_halo)
                    allocate(tmp_dblarr(3, 3))

                    !----------------------------------------
                    ! Get particles in neighboorhood
                    !----------------------------------------
                    call filter_region(center=pos_mean, width=pos_std,&
                         idsin=  ids_in_box,      in=  pos_in_box, &
                         idsout= ids_around_halo, out= pos_around_halo, &
                         parts_in_region= parts_in_region, order= order)

                    !-------------------------------------
                    ! Estimate density
                    !-------------------------------------
                    allocate(grid(infos%ndim, 100))
                    call conv_density(data=pos_around_halo(:, :parts_in_region),&
                         nbin=param_nbin, dens=density, edges=edges)

                    !----------------------------------------
                    ! Compute the FFT
                    !----------------------------------------
                    call conv%init_A(density)

                    ! TODO: loop over sigmas
                    ! TODO: convert sizes into sigmas
                    sigma = 1
                    call kernel_gaussian3d(param_nbin, sigma, gaussian)
                    call conv%init_B(gaussian)
                    call conv%execute()

                    ! copy back the convolution data
                    conv_dens = conv%conv
                    call conv%free()

                    deallocate(grid, density)

                 else
                    if (param_verbosity >= 3) then
                       write(*, '(a,i10,a,i10,a,i10,a)') 'halo n°', halo_i, ' is incomplete: ', &
                            counter, '/', members(halo_i)%parts, ' particles'
                    end if
                 end if
                 deallocate(pos_in_halo, m_in_halo)
                 deallocate(pos_around_halo)
              end if
           end do

           deallocate(order)
           deallocate(m_in_box, pos_in_box, ids_in_box)
        end do
     end do
  end do

  close(10)
  close(123)

contains
  subroutine compute_mean(pos, mean)
    real(dp), dimension(:, :), intent(in) :: pos
    real(dp), dimension(3), intent(out)   :: mean

    integer :: i, count
    count = 0
    mean = 0
    do i = 1, size(pos, 2)
       if (sum(pos(:, i)**2) > 0) then
          mean = mean + pos(:, i)
          count = count + 1
       end if
    end do
    if (count > 0) then
       mean = mean / (1d0*count)
    else
       print*, 'W: count is null', pos(:, :10)
    end if
  end subroutine compute_mean

end program compute_halo_prop
