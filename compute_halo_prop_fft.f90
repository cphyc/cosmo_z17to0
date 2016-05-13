program compute_halo_prop
  use io
  use misc
  use types
  use flap, only : command_line_interface
  use compute
  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: halo_to_cpu_file, mergertree_file, param_output_path
  integer :: param_from, param_to, param_output_number
  integer, parameter :: NPARTICLE_TO_PROBE_HALO = 50, NCPU_PER_HALO = 10
  real(dp) :: param_min_m, param_max_m
  integer :: param_verbosity

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
  ! Merger tree data
  !-------------------------------------
  ! integer,  allocatable, dimension(:)      :: mt_nhalos, mt_nsubhalos
  ! real(sp), dimension(:), allocatable      :: mt_aexp, mt_omega_t, mt_age_univ, time
  ! real(sp), dimension(:, :), allocatable   :: mt_pos, mt_vel, initial_pos
  ! integer                                  :: mt_nsteps, nsteps, nhalos
  ! integer                                  :: initial_halo_i, current_halo, parent_halo, istep, prev, father, step, max_nhalo
  ! integer, dimension(:), allocatable       :: tmp_nhalos, halos_z0
  ! integer, dimension(:), allocatable       :: parent, parent_at_step
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
  real(dp), dimension(:,:,:), allocatable    :: I_t
  real(dp), dimension(:), allocatable        :: m_in_halo, m_in_box
  real(dp), dimension(:, :), allocatable     :: pos_in_halo
  integer, dimension(:), allocatable         :: ids_in_box
  real(dp), dimension(:, :), allocatable     :: pos_in_box
  real(dp)                                   :: mtot
  real(dp), dimension(3)                     :: pos_mean, pos_std
  integer                                    :: halo_i, halo_counter, part_i
  integer                                    :: part_counter
  logical, dimension(:), allocatable         :: halo_found_mask
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                    :: i, counter, cpu_i, x, y, z, dim
  integer                                    :: tmp_int, index
  character(len=200)                         :: tmp_char
  real                                       :: factor
  real(dp), allocatable, dimension(:,:)      :: tmp_dblarr
  real(dp), dimension(3)                     :: X0, X1
  real(dp)                                   :: margin
  integer                                    :: param_nstep, max_nparts, step_i

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
  call cli%add(switch='--cpu-from', help='First cpu to use', act='store', def='1')
  call cli%add(switch='--cpu-to', help='Last cpu to use', act='store', def='4096')
  call cli%add(switch='--info-file', help='Information file about simulation', &
       act='store', def='/data52/Horizon-AGN/OUTPUT_DIR/output_00782/info_00782.txt')
  call cli%add(switch='--margin', help='margin', act='store', def='0.15')
  call cli%add(switch='--max-mass', switch_ab='-maxm', help='Maximum mass', &
       act='store',def='1e99')
  call cli%add(switch='--min-mass', switch_ab='-minm', help='Minimum mass', &
       act='store',def='0')
  call cli%add(switch='--nstep', act='store', def='4', &
       help='Number of step in which to compute the data')
  call cli%add(switch='--output', switch_ab='-o', help='Name of the output file', &
       act='store', def='data.out')
  call cli%add(switch='--output-number', help='Number of the output(s)',&
       act='store', def='2', nargs='+')
  call cli%add(switch='--output-path', help='Path of the simulation output',&
       act='store', def='/data52/Horizon-AGN/OUTPUT_DIR')
  call cli%add(switch='--verbose', help='Verbosity', act='store', def='0')

  call cli%get(switch='--cpu-to', val=param_to)
  call cli%get(switch='--cpu-from', val=param_from)
  call cli%get(switch='--min-mass', val=param_min_m)
  call cli%get(switch='--max-mass', val=param_max_m)
  call cli%get(switch='--output-path', val=param_output_path)
  call cli%get(switch='--output-number', val=param_output_number)
  call cli%get(switch='--verbose', val=param_verbosity)
  call cli%get(switch='--margin', val=margin)
  call cli%get(switch='--nstep', val=param_nstep)


  !-------------------------------------
  ! Read brick file
  !-------------------------------------
  print*, ''
  print*, 'Reading info file'
  call read_info_headers(param_output_path, param_output_number, infos)

  call cli%get(switch='--brick', val=tmp_char)
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
  allocate(I_t(3, 3, nDM))
  I_t = 0

  if (param_nstep == 1) then
     max_nparts = 2**30
  else
     max_nparts = 2**31 / (param_nstep**3)
  end if

  allocate (halo_found_mask(nDM))
  halo_found_mask = .false.
  !-------------------------------------
  ! Split the box into cells of equal sizes
  !-------------------------------------
  !$OMP PARALLEL DO DEFAULT(firstprivate) &
  !$OMP SHARED(infos, members, idDM, posDM, halo_found_mask, margin) &
  !$OMP PRIVATE(x, y, z, cpu_list, i)
  do x = 0, param_nstep-1
     do y = 0, param_nstep-1
        do z = 0, param_nstep-1
           X0 = (/ x*1d0/param_nstep, y*1d0/param_nstep, z*1d0/param_nstep /)
           X1 = (/ (x+1)*1d0/param_nstep, (y+1)*1d0/param_nstep, (z+1)*1d0/param_nstep /) + margin

           call get_cpu_list(X0, X1, infos%levelmax, infos%bound_key, &
                cpu_list, infos%ncpu, infos%ndim)

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
           !-------------------------------------
           ! Loop over all cpus listed, read the particle and save them
           !-------------------------------------
           part_counter = 0
           do cpu_i = 1, infos%ncpu
              if (cpu_list(cpu_i) == 0) then
                 cycle
              end if

              call read_particle(param_output_path, param_output_number, cpu_list(cpu_i), &
                   nstar, pos, vel, m, ids, birth_date, ndim, nparts)

              if (param_verbosity >= 1) then
                 write(*, '(a,          i5,    a,   i5,      a,    i5, a, i10, a, i10, a)') &
                      'Reading cpu', cpu_i, '/', tmp_int, ' (cpu=', cpu_list(cpu_i), &
                      ', nparts=', nparts, ', loaded=',part_counter,')'
              end if
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

           end do

           allocate(order(part_counter))
           print*, 'Sorting array…'
           ! Sort the ids
           call quick_sort(ids_in_box(1:part_counter), order)

           !-------------------------------------
           ! Once the particle read, for all complete halos
           ! compute the properties
           !-------------------------------------
           do halo_i = 1, nDM
              if ( mDM(halo_i) > param_min_m .and. mDM(halo_i) < param_max_m .and. &
                   posDM(1, halo_i) >= X0(1) .and. posDM(1, halo_i) < X1(1) .and. &
                   posDM(2, halo_i) >= X0(2) .and. posDM(2, halo_i) < X1(2) .and. &
                   posDM(3, halo_i) >= X0(3) .and. posDM(3, halo_i) < X1(3) .and. &
                   (.not. halo_found_mask(halo_i)) ) then
                 allocate(pos_in_halo(infos%ndim, members(halo_i)%parts), &
                      m_in_halo(members(halo_i)%parts))

                 counter = 0
                 do part_i = 1, members(halo_i)%parts
                    index = indexOf(members(halo_i)%ids(part_i), ids_in_box(1:part_counter))
                    if (index > 0) then
                       pos_in_halo(:, part_i) = pos_in_box(:, order(index))
                       m_in_halo(part_i)      = m_in_box(order(index))
                       counter = counter + 1
                    end if
                 end do
                 if (counter == members(halo_i)%parts) then
                    if (param_verbosity >= 2) then
                       print*, halo_i, 'is complete'
                    end if
                    !$OMP CRITICAL
                    halo_found_mask(halo_i) = .true.
                    !$OMP END CRITICAL
                    call correct_positions(pos_in_halo)
                    call compute_mean(pos_in_halo, pos_mean)
                    do dim = 1, infos%ndim
                       call stddev(pos_in_halo(dim, :) - pos_mean(dim), pos_std(dim), &
                            members(halo_i)%parts)
                    end do

                    mtot = sum(m_in_halo)
                    allocate(tmp_dblarr(3, 3))

                    ! !----------------------------------------
                    ! ! Read particles in neighboorhood
                    ! !----------------------------------------
                    ! call read_region(pos_mean, )
                    !-------------------------------------
                    ! Compute
                    !-------------------------------------
                    call compute_inertia_tensor(m_in_halo, pos_in_halo, tmp_dblarr)
                    I_t(:, :, halo_i) = tmp_dblarr
                    write(10, '(i12, 11ES14.6e2)') idDM(halo_i), mDM(halo_i), &
                         pos_mean(1), pos_mean(2), pos_mean(3), &
                         tmp_dblarr(1, 1), tmp_dblarr(1, 2), tmp_dblarr(1, 3), &
                                           tmp_dblarr(2, 2), tmp_dblarr(2, 3), &
                                                             tmp_dblarr(3, 3)
                    deallocate(tmp_dblarr)
                 else
                    if (param_verbosity >= 3) then
                       write(*, '(a,i10,a,i10,a,i10,a)') 'halo n°', halo_i, ' is incomplete: ', &
                            counter, '/', members(halo_i)%parts, ' particles'
                    end if
                 end if
                 deallocate(pos_in_halo, m_in_halo)
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
    do i = 1, members(halo_i)%parts
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
