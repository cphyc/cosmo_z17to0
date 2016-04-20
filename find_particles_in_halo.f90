program compute_halo_prop
  use io
  use misc
  use flap, only : command_line_interface
  use compute
  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: ramses_output_end, ramses_output_start, gal_list_filename, &
       associations_filename, dm_halo_list_filename, brick_file, info_file_end,&
       info_file_start, outfile, halo_to_cpu_file, mergertree_file, param_output_path
  integer            :: param_from, param_to, param_output_number, param_halo_i
  integer, allocatable :: param_output_number_list(:)
  integer, parameter :: NPARTICLE_TO_PROBE_HALO = 50, NCPU_PER_HALO = 10
  real(kind=8)       :: param_min_m, param_max_m
  integer            :: param_verbosity
  logical            :: stop_flag

  type(command_line_interface) :: cli
  !-------------------------------------
  ! List data
  !-------------------------------------
  real(kind=4), dimension(:,:), allocatable :: data_gal, data_associations, data_halo
  integer                                   :: ngal, nassoc, ndm_halo, gal_cols, assoc_cols, dm_halo_cols
  !-------------------------------------
  ! Brick data
  !-------------------------------------
  real(kind=8), dimension(:, :), allocatable :: posDM, LDM
  integer, dimension(:), allocatable         :: idDM
  real(kind=8), dimension(:), allocatable    :: rvirDM, mDM, mvirDM, TvirDM, hlevel
  type(INFOS_T)                              :: infos
  type(MEMBERS_T), dimension(:), allocatable :: members
  integer                                    :: nbodies, nb_of_halos, nb_of_subhalos, nDM
  real(kind=4)                               :: aexp_tmp, age_univ
  !-------------------------------------
  ! Merger tree data
  !-------------------------------------
  ! integer,  allocatable, dimension(:)        :: mt_nhalos, mt_nsubhalos
  ! real(kind=4), dimension(:), allocatable    :: mt_aexp, mt_omega_t, mt_age_univ, time
  ! real(kind=4), dimension(:, :), allocatable :: mt_pos, mt_vel, initial_pos
  ! integer                                    :: mt_nsteps, nsteps, nhalos
  ! integer                                    :: initial_halo_i, current_halo, parent_halo, istep, prev, father, step, max_nhalo
  ! integer, dimension(:), allocatable         :: tmp_nhalos, halos_z0
  ! integer, dimension(:), allocatable         :: parent, parent_at_step
  !-------------------------------------
  ! Particle data
  !-------------------------------------
  integer                                    :: ndim, nparts
  real(kind=8), dimension(:, :), allocatable :: pos, vel, tmp_pos
  integer                                    :: nstar, halo_found
  integer,      dimension(:), allocatable    :: ids, order
  real(kind=8), dimension(:), allocatable    :: m, tmp_mass
  real(kind=4), dimension(:), allocatable    :: birth_date
  !-------------------------------------
  ! Halo to cpu
  !-------------------------------------
  integer                               :: n_cpu_per_halo
  integer, dimension(:, :), allocatable :: halo_to_cpu
  integer, dimension(:), allocatable    :: cpu_list
  !-------------------------------------
  ! Halo properties
  !-------------------------------------
  real(kind=8), dimension(3, 3)               :: I_t_diag
  real(kind=8), dimension(:,:,:), allocatable :: I_t
  real(kind=8), dimension(:), allocatable     :: m_in_halo, m_in_box
  real(kind=8), dimension(:, :), allocatable  :: pos_in_halo
  integer, dimension(:), allocatable          :: ids_in_box, ids_in_halo
  real(kind=8), dimension(:, :), allocatable  :: pos_in_box
  real(kind=8)                                :: mtot
  real(kind=8), dimension(3)                  :: std, pos_mean
  integer                                     :: halo_i, ntot_halo, halo_counter, part_i
  integer                                     :: part_counter
  integer, dimension(:), allocatable          :: halo_found_mask
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                   :: i, j, cpu, i1, i2, counter, x, y, z
  integer                                   :: tmp_int, unit, tmp_int2, index
  character(len=200)                        :: tmp_char
  real                                      :: tmp_real, factor
  logical                                   :: tmp_bool
  real(kind=8), allocatable, dimension(:)   :: tmp_arr, tmp_arr2
  real(kind=8), allocatable, dimension(:,:) :: tmp_dblarr
  integer, allocatable, dimension(:)        :: tmp_iarr
  logical, dimension(4096)                  :: cpu_read
  real(kind=8), dimension(3)                :: X0, X1, center
  real(kind=8)                              :: margin
  integer                                   :: nstep, max_nparts, step_i

  !-------------------------------------
  ! random
  !-------------------------------------
  call random_seed()

  !-------------------------------------
  ! Read parameters
  !-------------------------------------
  call parse_params(cli)
  call cli%get(switch='--output-path', val=param_output_path)
  call cli%get_varying(switch='--output-number', val=param_output_number_list)
  call cli%get(switch='--verbose', val=param_verbosity)
  call cli%get(switch='--halo-i', val=param_halo_i)

  !-------------------------------------
  ! Read lists
  !-------------------------------------
  call cli%get(switch='--gal-list', val=tmp_char)
  print*, 'Reading file "' // trim(tmp_char) // '"'
  call read_list_header(trim(tmp_char), unit, ngal, gal_cols)
  allocate(data_gal(ngal, gal_cols+1))
  call read_list_data(unit, ngal, gal_cols, data_gal)
  print*, 'Got', ngal, 'galaxies in da pocket'

  call cli%get(switch='--association-list', val=tmp_char)
  print*, ''
  print*, 'Reading file "' // trim(tmp_char) // '"'
  call read_list_header(trim(tmp_char), unit, nassoc, assoc_cols)
  allocate(data_associations(nassoc, assoc_cols))
  call read_list_data(unit, nassoc, assoc_cols, data_associations)
  print*, 'Got', nassoc, 'associations in da pocket'

  call cli%get(switch='--halo-list', val=tmp_char)
  print*, ''
  print*, 'Reading file "' // trim(tmp_char) // '"'
  call read_list_header(trim(tmp_char), unit, ndm_halo, dm_halo_cols)
  allocate(data_halo(ndm_halo, dm_halo_cols))
  call read_list_data(unit, ndm_halo, dm_halo_cols, data_halo)
  print*, 'Got', ndm_halo, 'halos in da pocket'

  !-------------------------------------
  ! Read brick file
  !-------------------------------------
  call cli%get(switch='--info-file', val=tmp_char)
  print*, ''
  print*, 'Reading brick file "' // trim(tmp_char) // '"…'
  call read_info_headers(tmp_char, infos)

  call cli%get(switch='--brick', val=tmp_char)
  call read_brick_header(tmp_char, infos, nbodies, aexp_tmp, age_univ,&
       nb_of_halos, nb_of_subhalos)
  nDM = nb_of_halos + nb_of_subhalos
  allocate(idDM(nDM))
  allocate(posDM(infos%ndim, nDM))
  allocate(rvirDM(nDM))
  allocate(mDM(nDM))
  allocate(mvirDM(nDM))
  allocate(TvirDM(nDM))
  allocate(hlevel(nDM))
  allocate(LDM(infos%ndim, nDM))
  allocate(members(nDM))
  allocate(cpu_list(infos%ncpu))
  call read_brick_data(nDM, infos, .true., &
       & mDM, posDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LDM, idDM, members)
  print*, '    …red!'
  deallocate(rvirDM, mvirDM, TvirDM, hlevel, LDM)

  print*, ''
  call cli%get(switch='--halo-to-cpu', val=tmp_char)
  print*, 'Reading halo to cpu file "' // trim(tmp_char) // '"'
  call read_list_header(trim(tmp_char), unit, nDM, n_cpu_per_halo)
  allocate(halo_to_cpu(nDM, n_cpu_per_halo))
  call read_list_data(unit, nDM, n_cpu_per_halo, halo_to_cpu)
  ntot_halo = 0
  do i = 1, nDM
     if (halo_to_cpu(i, 1) > 0) then
        ntot_halo = ntot_halo + 1
     end if
  end do
  print*, '        Found', ntot_halo, 'halos!'

  cpu_read = .false.

  if (param_halo_i <= 0) then
     print*, 'Halo_i parameter must be positive.'
     stop 0
  end if

  !-------------------------------------
  ! Find the index of the halo given
  !-------------------------------------
  halo_i = indexOf(param_halo_i, idDM)
  if (halo_i <= 0) then
     print*, 'Halo not found. Stopping'
     stop
  end if

  allocate(pos_in_halo(infos%ndim, members(halo_i)%parts), ids_in_halo(members(halo_i)%parts))
  center = posDM(:, halo_i)
  do j = 1, size(param_output_number_list)
     param_output_number = param_output_number_list(j)

     !-------------------------------------
     ! Reading information file
     !-------------------------------------
     call cli%get(switch='--output-path', val=tmp_char)
     write(tmp_char, '(a, a, i0.5,a,i0.5,a)') trim(tmp_char), '/output_', param_output_number, &
          '/info_', param_output_number, '.txt'

     if (param_verbosity >= 2) then
        print*, ''
        print*, 'Reading info file "' // trim(tmp_char) // '".'
     end if
     call read_info_headers(tmp_char, infos)

     !-------------------------------------
     ! Open output file
     !-------------------------------------
     call cli%get(switch='--output', val=tmp_char)
     write(tmp_char, '(a,i0.5)') trim(tmp_char), param_output_number

     if (param_verbosity >= 2) then
        print*, ''
        write(*, '(a,i5,a,a,x,a,3ES14.6e2)') 'Working on output', param_output_number, &
             ', writing output in ', trim(tmp_char), 'center:', center
     end if

     open(unit=10, file=trim(tmp_char))
     write(10, '(a9, 20a13)') 'id', 'x', 'y', 'z'

     !-------------------------------------
     ! Iterate until all the particles are found
     !-------------------------------------
     margin = 0
     counter = 0
     pos_in_halo = 0

     allocate(order(members(halo_i)%parts))
     call quick_sort(members(halo_i)%ids, order)
     deallocate(order)

     do while (counter < members(halo_i)%parts)
        ! Read information of output

        if (param_verbosity >= 3) then
           write(*, '(a,i10,a,i10,a,i10,a)') 'halo n°', halo_i, ' is incomplete: ', &
                counter, '/', members(halo_i)%parts, ' particles'
        end if

        ! Increase margin each time
        margin = margin + 0.01
        X0 = center - margin
        X1 = center + margin

        ! Find the cpu for the box
        call get_cpu_list(X0, X1, infos%levelmax, infos%bound_key, cpu_list, infos%ncpu, infos%ndim)
        n_cpu_per_halo = 0

        ! a bool flag to know if there's any unread cpu
        tmp_bool = .false.

        ! Count the number of cpus that are going to be read
        do cpu = 1, infos%ncpu
           if (cpu_list(cpu) > 0) then
              n_cpu_per_halo = n_cpu_per_halo + 1

              ! if there is any cpu unread in the list, mark as true
              if (cpu_read(cpu) == .false.) then
                 tmp_bool = .true.
              end if
           end if
        end do
        if (n_cpu_per_halo == infos%ncpu .and. tmp_bool == .false.) then
           print*, 'No cpu to read!'
           exit
        end if

        stop_flag = .false.
        ! Read all cpus
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP SHARED(members, X0, X1, stop_flag, counter, infos, halo_i, cpu_read, param_verbosity) &
        !$OMP SHARED(n_cpu_per_halo, param_output_path, param_output_number, cpu_list) &
        !$OMP PRIVATE(order, nparts, tmp_int, ids, m, vel, pos, nstar, birth_date, ndim) &
        !$OMP REDUCTION(+:pos_in_halo) REDUCTION(+:ids_in_halo)&
        !$OMP SCHEDULE(dynamic, 10)
        do cpu = 1, infos%ncpu
           if (stop_flag) then
              cycle
           end if
           if (counter == members(halo_i)%parts) then
              stop_flag = .true.
           else if (counter > members(halo_i)%parts) then
              print*, 'W: something weird happend, found', counter, 'instead of', members(halo_i)%parts
              stop_flag = .true.
           end if

           ! for the cpus not already read
           if (cpu_list(cpu) > 0 .and. (.not. cpu_read(cpu))) then
              if (param_verbosity >= 3) then
                 write(*, '(a,i5,a,i5,a,i5,a,i5,a,i5,a)') 'Reading cpu n°', cpu_list(cpu),&
                      ' (cpu=', cpu, '/', n_cpu_per_halo, &
                      ', nparts=', counter, '/', members(halo_i)%parts,')'
              end if
              ! read the particles
              call read_particle(param_output_path, param_output_number, cpu_list(cpu), &
                   nstar, pos, vel, m, ids, birth_date, ndim, nparts)

              allocate(order(nparts))
              call quick_sort(ids, order)

              ! Store the position of the particles in the halo
              do part_i = 1, nparts
                 tmp_int = indexOf(ids(part_i), members(halo_i)%ids)
                 if (tmp_int > 0) then
                    ids_in_halo(tmp_int)    = ids(part_i)
                    pos_in_halo(:, tmp_int) = pos(:, part_i)
                    !$OMP ATOMIC
                    counter = counter + 1
                 end if
              end do
              deallocate(order)
              cpu_read(cpu_list(cpu)) = .true.
           end if
        end do
        !$OMP END PARALLEL DO
     end do

     call correct_positions(pos_in_halo)

     print*, 'all found :D'
     do i = 1, members(halo_i)%parts
        write(10, '(i12, 3ES14.6e2)') ids_in_halo(i), pos_in_halo(:, i)
     end do
     close(10)

     center = sum(pos_in_halo, 2) / size(pos_in_halo, 2)
  end do


end program compute_halo_prop
