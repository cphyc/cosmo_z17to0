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
  real(kind=8)                                :: mtot
  real(kind=8), dimension(3) :: std, pos_mean
  integer                                     :: halo_i, ntot_halo, halo_counter
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                 :: i, j, cpu, i1, i2, counter
  integer                                 :: tmp_int, unit, tmp_int2
  character(len=200)                      :: tmp_char
  real                                    :: tmp_real, factor
  logical                                 :: tmp_bool
  real(kind=8), allocatable, dimension(:) :: tmp_arr, tmp_arr2
  real(kind=8), allocatable, dimension(:,:) :: tmp_dblarr
  integer, allocatable, dimension(:) :: tmp_iarr
  logical, dimension(4096) :: cpu_read

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
  print*, 'Computing inertia tensor'
  open(unit=10, file=trim(tmp_char))
  open(unit=123, file='missing_parts')
  write(10, '(a12, 13a14)') 'id', 'mass', 'x', 'y', 'z', 'size', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'
  halo_counter = 1
  allocate(I_t(3, 3, nDM))
  I_t = 0
  !$OMP PARALLEL DO default(firstprivate) shared(halo_to_cpu, members, posDM, infos) &
  !$OMP shared(halo_counter, I_t)
  do halo_i = 1, nDM
     cpu_read = .false.
     if (halo_to_cpu(halo_i, 1) == 0) then
        cycle
     end if

     if (mDM(halo_i) < param_min_m .or. mDM(halo_i) > param_max_m) then
        cycle
     end if

     write(*, '(a, i7, a, F5.1, a)') 'halo', halo_i, '(', &
          100.*halo_counter/ntot_halo, '%)'
     !$OMP ATOMIC
     halo_counter = halo_counter + 1

     allocate(tmp_pos(infos%ndim, members(halo_i)%parts), tmp_mass(members(halo_i)%parts))
     counter  = 0
     tmp_pos  = 0
     tmp_mass = 0

     do j = 1, n_cpu_per_halo
        if (halo_to_cpu(halo_i, j) == 0) then
           cycle
        end if
        cpu = halo_to_cpu(halo_i, j)
        call treat_cpu(tmp_pos, tmp_mass, cpu)
        cpu_read(cpu) = .true.
     end do
     mtot = sum(tmp_mass)

     !-------------------------------------
     ! find missing particles, if any (trigger for 1% of missing particles)
     !-------------------------------------
     factor = 0
     call correct_positions(tmp_pos)
     do while (1d0*counter / members(halo_i)%parts < 0.95)
        factor = factor + 1d-4
        call compute_mean(tmp_pos, pos_mean)

        allocate(tmp_arr(infos%ndim), tmp_arr2(infos%ndim))
        tmp_arr = pos_mean - factor
        tmp_arr2 = pos_mean + factor
        call get_cpu_list(tmp_arr, tmp_arr2, infos%levelmax, infos%bound_key, &
             cpu_list, infos%ncpu, infos%ndim)
        ! print*, 'W: missing particles, looking them in '
        ! print*, tmp_arr
        ! print*, tmp_arr2
        ! print*, infos%levelmax, infos%ncpu, infos%ndim
        deallocate(tmp_arr, tmp_arr2)

        do i = 1, infos%ncpu
           if (cpu_list(i) > 0) then
              if (.not. cpu_read(cpu_list(i))) then
                 ! print*, 'Reading further', cpu_list(i)
              end if
           end if
        end do
        do i = 1, infos%ncpu
           if (cpu_list(i) > 0) then
              if (.not. cpu_read(cpu_list(i))) then
                 if (counter == members(halo_i)%parts) then
                    exit
                 end if
                 cpu = cpu_list(i)
                 call treat_cpu(tmp_pos, tmp_mass, cpu)
                 cpu_read(cpu_list(i)) = .true.
              end if
           end if
        end do
        ! print*, 'CPUs'
        ! print*, halo_to_cpu(halo_i, :)
        ! print*, cpu_list(:tmp_int)
        ! write(123, *) idDM(halo_i)
        call correct_positions(tmp_pos)
     end do
     print*, 'Found', counter, 'particles (expected', members(halo_i)%parts, ')'
     !-------------------------------------
     ! compute inertia tensor
     !-------------------------------------
     allocate(tmp_dblarr(3, 3))
     call compute_inertia_tensor(tmp_mass, tmp_pos, tmp_dblarr)
     I_t(:, :, halo_i) = tmp_dblarr
     write(10, '(i12, 20ES14.6e2)') idDM(halo_i), mtot, pos_mean, factor, tmp_dblarr
     deallocate(tmp_dblarr)

     deallocate(tmp_pos, tmp_mass)
  end do
  !$OMP END PARALLEL DO
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
  subroutine treat_cpu(tmp_pos, tmp_mass, cpu)
    integer, intent(in) :: cpu

    real(kind=8), intent(inout), dimension(:, :) :: tmp_pos
    real(kind=8), intent(inout), dimension(:)    :: tmp_mass

    real(kind=8), dimension(:, :), allocatable :: pos, vel
    real(kind=8), dimension(:),    allocatable :: m
    real(kind=4), dimension(:),    allocatable :: birth_date
    integer     , dimension(:),    allocatable :: ids

    integer :: nstar, index
    integer :: i
    !-------------------------------------
    ! Reading cpu
    !-------------------------------------
    ! print*, 'Reading cpu', cpu, '(', counter, '/', members(halo_i)%parts, ')'
    call read_particle(param_output_path, param_output_number, cpu, nstar, pos, vel, m,&
         ids, birth_date, ndim, nparts)
    allocate(order(nparts))

    ! Filter out particles not in halo
    call quick_sort(ids, order)

    do i = 1, members(halo_i)%parts
       ! get the position of the halo_id in the ids given
       ! if found, store its velocity into our temporary array
       ! and add it into the total mass
       index = indexOf(members(halo_i)%ids(i), ids)
       if (index > 0) then
          tmp_pos(:, i) = pos(:, index)
          tmp_mass(i) = m(index)
          counter = counter + 1
       end if
    end do

    deallocate(order)
  end subroutine treat_cpu

end program compute_halo_prop
