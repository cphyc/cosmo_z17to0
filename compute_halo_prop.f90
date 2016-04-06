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
  real(kind=8), dimension(:), allocatable    :: m, birth_date, tmp_mass
  !-------------------------------------
  ! Halo to cpu
  !-------------------------------------
  integer                               :: n_cpu_per_halo
  integer, dimension(:, :), allocatable :: halo_to_cpu
  !-------------------------------------
  ! Halo properties
  !-------------------------------------
  real(kind=8), dimension(3, 3) :: I_t, I_t_diag
  real(kind=8)                  :: mtot
  integer                       :: halo_i, ntot_halo, halo_counter
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                 :: i, j, cpu, i1, i2, counter
  integer                                 :: tmp_int, unit, tmp_int2
  character(len=200)                      :: tmp_char
  real                                    :: tmp_real
  logical                                 :: tmp_bool

  !-------------------------------------
  ! random
  !-------------------------------------
  call random_seed()

  !-------------------------------------
  ! Read parameters
  !-------------------------------------
  call parse_params(cli)
  call cli%get(switch='--cpu-to', val=param_to)
  call cli%get(switch='--cpu-from', val=param_from)
  call cli%get(switch='--min-mass', val=param_min_m)
  call cli%get(switch='--max-mass', val=param_max_m)
  call cli%get(switch='--output-path', val=param_output_path)
  call cli%get(switch='--output-number', val=param_output_number)

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
  call read_brick_data(nDM, infos, .true., &
       & mDM, posDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LDM, idDM, members)
  print*, '    …red!'
  deallocate(rvirDM, mDM, mvirDM, TvirDM, hlevel, LDM)

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

  print*, ''
  print*, 'Computing inertia tensor (test)'
  call cli%get(switch='--output', val=tmp_char)
  print*, 'Writing in', tmp_char
  open(unit=10, file=trim(tmp_char))
  halo_counter = 1
  !$OMP PARALLEL DO default(firstprivate) shared(halo_to_cpu, members, posDM, infos, halo_counter)
  do halo_i = 1, nDM
     if (halo_to_cpu(halo_i, 1) == 0) then
        cycle
     end if

     write(*, '(a, i7, a, F5.1, a)') 'halo', halo_i, '(', &
          100.*halo_counter/ntot_halo, '%)'
     !$OMP ATOMIC
     halo_counter = halo_counter + 1

     allocate(tmp_pos(infos%ndim, members(halo_i)%parts), tmp_mass(members(halo_i)%parts))
     I_t = 0
     counter = 0
     do j = 1, n_cpu_per_halo
        if (halo_to_cpu(halo_i, j) == 0) then
           exit
        end if
        cpu = halo_to_cpu(halo_i, j)
        !-------------------------------------
        ! Reading cpu
        !-------------------------------------
        call read_particle(param_output_path, param_output_number, cpu, nstar, pos, vel, m,&
             ids, birth_date, ndim, nparts)
        allocate(order(nparts))

        ! Filter out particles not in halo
        call quick_sort(ids, order)

        do i = 1, members(halo_i)%parts
           ! get the position of the halo_id in the ids given
           ! if found, store its velocity into our temporary array
           ! and add it into the total mass
           tmp_int = indexOf(members(halo_i)%ids(i), ids)
           if (tmp_int > 0) then
              tmp_pos(:, i) = pos(:, tmp_int)
              tmp_mass(i) = m(tmp_int)
              counter = counter + 1
           else
              tmp_pos(:, i) = 0
           end if
        end do


        deallocate(order)
     end do
     mtot = sum(tmp_mass)
     ! Iterate over each couples to populate I
     do i1 = 1, ndim
        do i2 = i1, ndim
           tmp_real = sum(tmp_mass*tmp_pos(i1, :)*tmp_pos(i2, :)) / mtot
           I_t(i1, i2) = tmp_real
           I_t(i2, i1) = tmp_real
        end do
     end do
     print*, 'Found', counter, 'particles (expected', members(halo_i)%parts, ')'
     ! print*, I_t
     !!$OMP CRITICAL
     write(10, '(i9, 9ES13.6e2)') idDM(halo_i), I_t
     !!$OMP END CRITICAL
     deallocate(tmp_pos, tmp_mass)
  end do
  !!$OMP END PARALLEL DO
  close(10)

  ! call compute_inertia_tensor(mDM, posDM, I_t, I_t_diag)
  ! do i = 1, 3
  !    print*, I_t(i, :)
  ! end do
  ! print*, ''
  ! do i = 1, 3
  !    print*, I_t_diag(i, :)
  ! end do

  ! tmp_int = 0
  ! do i = 1, nDM
  !    if ( (param_min_m == 0 .or. mDM(i) > param_min_m) .and. &
  !         (param_max_m == 0 .or. mDM(i) < param_max_m)) then
  !       tmp_int = tmp_int + 1
  !    end if
  ! end do

  ! print*, 'Found', tmp_int, 'with m > ', param_min_m
contains

end program compute_halo_prop
