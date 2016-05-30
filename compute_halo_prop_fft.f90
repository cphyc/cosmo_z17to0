program compute_halo_prop
  use io
  use misc
  use types
  use flap, only : command_line_interface
  use compute
  use convolution
  use extrema_mod
  use extrema_types, only : EXT_DATA, CND_CNTRL_TYPE
  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: param_output_path
  integer            :: param_output_number, unit_output
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
  real(sp), dimension(:), allocatable        :: flattened_field
  real(dp), dimension(:), allocatable        :: m_in_halo, m_in_box
  real(dp), dimension(:, :), allocatable     :: pos_in_halo,&
       & pos_around_halo
  integer, dimension(:), allocatable         :: ids_in_box, ids_around_halo
  real(dp), dimension(:, :), allocatable     :: pos_in_box
  real(dp)                                   :: mtot
  real(dp), dimension(3)                     :: pos_mean, pos_std
  integer                                    :: halo_i, halo_counter, part_i
  integer                                    :: part_counter, parts_in_region
  logical, dimension(:), allocatable         :: halo_found_mask
  type(CONV_T)                               :: conv
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                    :: i, counter, cpu_i, x, y, z, dim
  integer                                    :: tmp_int, index, tmp_int2
  character(len=200)                         :: tmp_char
  real(dp), dimension(3)                     :: X0, X1
  integer                                    :: max_nparts, step_i

  !----------------------------------------
  ! FFT variables
  !----------------------------------------
  real(dp), allocatable :: gaussian(:, :, :), conv_dens(:, :, :), edges(:, :)
  real(dp) :: sigma, param_around
  integer  :: param_nbin, param_nsigma, param_sigma_min, param_sigma_max, isigma
  !----------------------------------------
  ! Peaks
  !----------------------------------------
  type(EXT_DATA), allocatable, dimension(:) :: extrema
  type(CND_CNTRL_TYPE)                      :: extrema_ctrl



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
  call cli%add(switch='--max-mass', switch_ab='-maxm', help='Maximum mass', &
       act='store',def='1e99')
  call cli%add(switch='--min-mass', switch_ab='-minm', help='Minimum mass', &
       act='store',def='0')
  call cli%add(switch='--nbin', act='store', def='32', &
       help='Number of bins to use to estimate density')
  call cli%add(switch='--output', switch_ab='-o', help='Name of the output file', &
       act='store', def='data.out')
  call cli%add(switch='--output-number', help='Number of the output(s)',&
       act='store', def='2', nargs='+')
  call cli%add(switch='--output-path', help='Path of the simulation output',&
       act='store', def='/data52/Horizon-AGN/OUTPUT_DIR')
  call cli%add(switch='--verbose', help='Verbosity', act='store', def='0')
  call cli%add(switch='--sigma-min', help='Minimum sigma (ramses units)', act='store', &
       required=.true.)
  call cli%add(switch='--sigma-max', help='Maximum sigma (ramses units)', act='store', &
       required=.true.)
  call cli%add(switch='--nsigma', help='Number of sigma step to do', def='10')
  call cli%add(switch='--around', help='How much space to explore around halos (ramses units)', act='store', required=.true.)

  call cli%get(switch='--min-mass', val=param_min_m)
  call cli%get(switch='--max-mass', val=param_max_m)
  call cli%get(switch='--output-path', val=param_output_path)
  call cli%get(switch='--output-number', val=param_output_number)
  call cli%get(switch='--verbose', val=param_verbosity)
  call cli%get(switch='--nbin', val=param_nbin)
  call cli%get(switch='--sigma-min', val=param_sigma_min)
  call cli%get(switch='--sigma-min', val=param_sigma_max)
  call cli%get(switch='--nsigma', val=param_nsigma)
  call cli%get(switch='--around', val=param_around)

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
  open(newunit=unit_output, file=trim(tmp_char))
  write(unit_output, '(a7, a5, a14, 3a8, 7a14)' ) &
       'halo_i', 'type', 'sigma', 'xpixel', 'ypixel', 'zpixel', &
       'x', 'y', 'z', 'ex', 'ey', 'ez', 'eigval'

  halo_counter = 1

  max_nparts = 2**30


  !----------------------------------------
  ! Allocate data
  !----------------------------------------
  allocate(gaussian(param_nbin, param_nbin, param_nbin), &
       conv_dens(param_nbin, param_nbin, param_nbin), &
       density(param_nbin, param_nbin, param_nbin), &
       flattened_field(param_nbin**3))
  allocate(edges(infos%ndim, param_nbin+1))
  allocate (halo_found_mask(nDM))
  halo_found_mask = .false.


  !-------------------------------------
  ! Split the box into cells of equal sizes
  !-------------------------------------
  X0 = (/0, 0, 0/)
  X1 = (/1, 1, 1/)
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

  write(*, '(a, i5,1x, a, 2(3f7.4,1x))') &
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
  !$OMP PRIVATE(tmp_int)                                              &
  !$OMP FIRSTPRIVATE(unit_output)                                     &
  !$OMP SCHEDULE(guided, 1)
  do cpu_i = 1, 1000!infos%ncpu
     if (cpu_list(cpu_i) == 0) then
        cycle
     end if

     !$OMP ATOMIC
     tmp_int2 = tmp_int2 + 1 ! only a passive counter

     call read_particle(param_output_path, param_output_number, cpu_list(cpu_i), &
          nstar, pos, vel, m, ids, birth_date, ndim, nparts)

     if (param_verbosity >= 1) then
        write(*, '(a,          i5,    a,   i5,      a,    i5, a, i10, a, i10, a)') &
             'Reading cpu', tmp_int2, '/', infos%ncpu, ' (cpu=', cpu_list(cpu_i), &
             ', nparts=', nparts, ', loaded=', part_counter,')'
     end if
     do i = 1, nparts
        if ( pos(1, i) >= X0(1) .and. pos(1, i) < X1(1) .and. &
             pos(2, i) >= X0(2) .and. pos(2, i) < X1(2) .and. &
             pos(3, i) >= X0(3) .and. pos(3, i) < X1(3) .and. &
             ids(i) > 0 & ! only keep DM particles
             ) then

           !$OMP ATOMIC
           part_counter = part_counter + 1
           pos_in_box(:, part_counter) = pos(:, i)
           m_in_box(part_counter)      = m(i)
           ids_in_box(part_counter)    = ids(i)
        end if
     end do
  end do
  !$OMP END PARALLEL DO

  allocate(order(part_counter))
  print*, 'Sorting array (len=', part_counter, ')…'
  ! Sort the ids
  call quick_sort(ids_in_box(1:part_counter), order(1:part_counter))
  print*, '   … sorted'

  !-------------------------------------
  ! Once the particle read, for all complete halos
  ! compute the properties
  !-------------------------------------
  !$OMP PARALLEL DO DEFAULT(none)                                              &
  !$OMP SHARED(mDM, param_min_m, param_max_m, X0, X1, posDM, halo_found_mask)  &
  !$OMP SHARED(infos, param_nbin, ids_in_box, members, part_counter, m_in_box) &
  !$OMP SHARED(nDM, max_nparts, pos_in_box, param_verbosity, order)            &
  !$OMP SHARED(param_nsigma, param_sigma_max, param_sigma_min, unit_output)    &
  !$OMP SHARED(param_around)                                                   &
  !$OMP PRIVATE(pos_in_halo, m_in_halo, pos_around_halo, ids_around_halo)      &
  !$OMP PRIVATE(mtot, parts_in_region, edges, conv, conv_dens)                 &
  !$OMP PRIVATE(sigma, density, gaussian, index, pos_mean, pos_std, counter)   &
  !$OMP PRIVATE(flattened_field, extrema_ctrl, extrema)                        &
  !$OMP SCHEDULE(guided, 1)
  do halo_i = 1, nDM
     ! filter halos outside mass range and not in box
     if ( mDM(halo_i) > param_min_m .and. mDM(halo_i) < param_max_m .and. &
          (.not. halo_found_mask(halo_i)) ) then

        ! allocate data for the halo
        allocate(pos_in_halo(infos%ndim, members(halo_i)%parts))
        allocate(m_in_halo(members(halo_i)%parts))
        allocate(pos_around_halo(infos%ndim, max_nparts))
        allocate(ids_around_halo(max_nparts))

        ! get the particles in the region around the halo
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

           ! correct the positions and get center + stddev
           call correct_positions(pos_in_halo)
           call compute_mean(pos_in_halo, pos_mean)
           ! do dim = 1, infos%ndim
           !    call stddev(pos_in_halo(dim, :) - pos_mean(dim), pos_std(dim), &
           !         members(halo_i)%parts)
           ! end do

           mtot = sum(m_in_halo)

           !----------------------------------------
           ! Get particles in neighboorhood
           !----------------------------------------
           if (param_verbosity >= 4) write(*, *) halo_i, ': filtering'

           call filter_region(center=pos_mean, &
                width=(/param_around, param_around, param_around/), &
                idsin=  ids_in_box(1:part_counter),  in=  pos_in_box(:,1:part_counter), &
                idsout= ids_around_halo, out= pos_around_halo, &
                parts_in_region= parts_in_region, order= order)

           !-------------------------------------
           ! Estimate density
           !-------------------------------------
           if (param_verbosity >= 4) write(*, *) halo_i, ': density estimation'
           call conv_density(data=pos_around_halo(:, :parts_in_region),&
                nbin=param_nbin, dens=density, edges=edges)

           !----------------------------------------
           ! Compute the FFT
           !----------------------------------------
           if (param_verbosity >= 4) write(*, *) halo_i, ': fft of density'
           call conv%init_A(density)

           !----------------------------------------
           ! Iterate over sigma
           !----------------------------------------
           allocate(extrema(param_nbin**3))
           do isigma = 1, param_nsigma
              sigma = (param_sigma_max - param_sigma_min) * (isigma - 1) &
                   / (param_nsigma - 1) + param_sigma_min

              if (param_verbosity >= 4) then
                 write(*, '(i7,a,ES10.3e2,a,i3,a,i3,a)')&
                      halo_i, ': fft of gaussian kernel (sigma=', sigma, ', ',&
                      isigma, '/', param_nsigma, ')'
              end if
              call kernel_gaussian3d(param_nbin, sigma, gaussian)
              call conv%init_B(gaussian)

              if (param_verbosity >= 4) then
                 write(*, '(i7,a,i3,a,i3,a)')&
                      halo_i, ': convolution (', &
                      isigma, '/', param_nsigma, ')'
              end if

              call conv%execute()

              ! compute and get the extrema
              extrema_ctrl%nproc = 1
              extrema_ctrl%justprint = .false.

              flattened_field = reshape(real(conv%conv), (/ size(conv%conv) /))
              call find_extrema(&
                   dt=flattened_field, &
                   nn=(/param_nbin, param_nbin, param_nbin/), &
                   nd=3, &
                   ext=extrema, &
                   cnd_cntrl=extrema_ctrl)

              ! save output
              do i = 1, param_nbin**3
                 if (extrema(i)%typ > 0) then
                    write(unit_output, '(i7, i5, ES14.6e2, 3i8, 7ES14.6e2)' )&
                         halo_i, extrema(i)%typ, sigma, extrema(i)%pix, &
                         extrema(i)%pos, extrema(i)%eig, extrema(i)%val
                 end if
              end do

           end do
           ! free memory
           if (param_verbosity >= 4) write(*, *) halo_i, ': free'
           call conv%free()
           deallocate(extrema)
        else
           if (param_verbosity >= 3) then
              write(*, '(a,i10,a,i10,a,i10,a)') 'halo n°', halo_i, ' is incomplete: ', &
                   counter, '/', members(halo_i)%parts, ' particles'
           end if
        end if
        deallocate(pos_in_halo, m_in_halo)
        deallocate(pos_around_halo)
        deallocate(ids_around_halo)
     end if
  end do

  deallocate(order)
  deallocate(m_in_box, pos_in_box, ids_in_box)

  close(unit_output)

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
