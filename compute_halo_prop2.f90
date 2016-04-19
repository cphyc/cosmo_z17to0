program compute_halo_prop
  use io
  use misc
  use flap, only : command_line_interface
  use compute
  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: ramses_output_end, ramses_output_start, gal_list_filename, associations_filename, dm_halo_list_filename, brick_file, info_file_end, info_file_start, outfile, halo_to_cpu_file, mergertree_file, param_output_path
  integer :: param_from, param_to, param_output_number
  integer, parameter :: NPARTICLE_TO_PROBE_HALO = 50, NCPU_PER_HALO = 10
  real(kind=8) :: param_min_m, param_max_m
  integer :: param_verbosity

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
  integer, dimension(:), allocatable          :: ids_in_box
  real(kind=8), dimension(:, :), allocatable  :: pos_in_box
  real(kind=8)                                :: mtot
  real(kind=8), dimension(3)                  :: std, pos_mean
  integer                                     :: halo_i, ntot_halo, halo_counter, part_i
  integer                                     :: part_counter
  integer, dimension(:), allocatable          :: halo_found_mask
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                   :: i, j, cpu, i1, i2, counter, cpu_i, x, y, z
  integer                                   :: tmp_int, unit, tmp_int2, index
  character(len=200)                        :: tmp_char
  real                                      :: tmp_real, factor
  logical                                   :: tmp_bool
  real(kind=8), allocatable, dimension(:)   :: tmp_arr, tmp_arr2
  real(kind=8), allocatable, dimension(:,:) :: tmp_dblarr
  integer, allocatable, dimension(:)        :: tmp_iarr
  logical, dimension(4096)                  :: cpu_read
  real(kind=8), dimension(3)                :: X0, X1
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
  call cli%add(switch='--margin', help='margin', &
       act='store', def='0.15')
  call cli%add(switch='--nstep', help='Nstep', &
       act='store', def='4')
  call cli%get(switch='--cpu-to', val=param_to)
  call cli%get(switch='--cpu-from', val=param_from)
  call cli%get(switch='--min-mass', val=param_min_m)
  call cli%get(switch='--max-mass', val=param_max_m)
  call cli%get(switch='--output-path', val=param_output_path)
  call cli%get(switch='--output-number', val=param_output_number)
  call cli%get(switch='--verbose', val=param_verbosity)
  call cli%get(switch='--margin', val=margin)
  call cli%get(switch='--nstep', val=nstep)

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

  call cli%get(switch='--output', val=tmp_char)
  print*, ''
  print*, 'Writing output in ', trim(tmp_char)
  open(unit=10, file=trim(tmp_char))
  open(unit=123, file='missing_parts')
  write(10, '(a9, 20a13)') 'id', 'mass', 'x', 'y', 'z', 'xx', 'xy', 'xz', 'yy', 'yz', 'zz'
  halo_counter = 1
  allocate(I_t(3, 3, nDM))
  I_t = 0

  max_nparts = 2**30

  allocate (halo_found_mask(nDM))
  halo_found_mask = .false.
  !-------------------------------------
  ! Split the box into cells of equal sizes
  !-------------------------------------
  !$OMP PARALLEL DO DEFAULT(firstprivate) &
  !$OMP SHARED(infos, members, idDM, posDM, halo_found_mask, margin) &
  !$OMP PRIVATE(x, y, z, cpu_list, i)
  do x = 0, nstep-1
     do y = 0, nstep-1
        do z = 0, nstep-1
           X0 = (/ x*1d0/nstep, y*1d0/nstep, z*1d0/nstep /)
           X1 = (/ (x+1)*1d0/nstep, (y+1)*1d0/nstep, (z+1)*1d0/nstep /) + margin

           call get_cpu_list(X0, X1, infos%levelmax, infos%bound_key, &
                cpu_list, infos%ncpu, infos%ndim)

           tmp_int = 0
           do cpu_i = 1, infos%ncpu
              if (cpu_list(cpu_i) > 0) then
                 tmp_int = tmp_int + 1
              end if
           end do
           step_i = x*nstep**2 + y*nstep + z + 1
           write(*, '(a, i4,    a,   i4,       x,  a, i5, x, a, 2(3f7.4,x))') &
                'step', step_i, '/', nstep**3, &
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

                    mtot = sum(m_in_halo)
                    allocate(tmp_dblarr(3, 3))
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
    real(kind=8), dimension(:, :), intent(in) :: pos
    real(kind=8), dimension(3), intent(out)   :: mean

    integer :: i, j, count
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
