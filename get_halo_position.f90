program sort_galaxy
  use io
  use misc
  use flap, only : command_line_interface

  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: ramses_output_end, ramses_output_start, gal_list_filename, &
       associations_filename, dm_halo_list_filename, brick_file, info_file,&
       outfile, halo_to_cpu_file, mergertree_file
  integer            :: param_from, param_to
  integer, parameter :: NCPU_PER_HALO = 10
  integer            :: nparticle_to_probe_halo, frac_particle_to_probe

  type(command_line_interface) :: cli

  !-------------------------------------
  ! List data
  !-------------------------------------
  real(kind=4), dimension(:,:), allocatable :: data_gal, data_associations, data_halo
  integer                                   :: ngal, nassoc, ndm_halo, gal_cols, assoc_cols, dm_halo_cols
  real(kind=8) :: param_min_m, param_max_m
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
  ! integer                                    :: initial_halo_id, current_halo, parent_halo, istep, prev, father, step, max_nhalo
  ! integer, dimension(:), allocatable         :: tmp_nhalos, halos_z0
  ! integer, dimension(:), allocatable         :: parent, parent_at_step
  !-------------------------------------
  ! Particle data
  !-------------------------------------
  integer                                    :: ndim, nparts
  real(kind=8), dimension(:, :), allocatable :: pos, vel
  integer                                    :: nstar, halo_found
  integer,      dimension(:), allocatable    :: ids, order, particles_found
  integer, dimension(:, :), allocatable      :: halo_to_cpu
  real(kind=8), dimension(:), allocatable    :: m, birth_date
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer            :: i, j, cpu
  integer            :: tmp_int, unit, tmp_int2
  character(len=200) :: tmp_char
  real               :: tmp_real
  logical            :: tmp_bool
  integer, dimension(:), allocatable :: tmp_arr
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
  call cli%get(switch='--gal-list', val=gal_list_filename)
  call cli%get(switch='--halo-list', val=dm_halo_list_filename)
  call cli%get(switch='--association-list', val=associations_filename)
  call cli%get(switch='--brick', val=brick_file)
  call cli%get(switch='--info-file', val=info_file)
  call cli%get(switch='--percent-probe-particle', val=frac_particle_to_probe)
  print*, param_min_m

  !-------------------------------------
  ! Read lists
  !-------------------------------------
  print*, 'Reading file "' // trim(gal_list_filename) // '"'
  call read_list_header(trim(gal_list_filename), unit, ngal, gal_cols)
  allocate(data_gal(ngal, gal_cols))
  call read_list_data(unit, ngal, gal_cols, data_gal)
  print*, 'Got', ngal, 'galaxies in da pocket'

  print*, ''
  print*, 'Reading file "' // trim(associations_filename) // '"'
  call read_list_header(trim(associations_filename), unit, nassoc, assoc_cols)
  allocate(data_associations(nassoc, assoc_cols))
  call read_list_data(unit, nassoc, assoc_cols, data_associations)
  print*, 'Got', nassoc, 'associations in da pocket'

  print*, ''
  print*, 'Reading file "' // trim(dm_halo_list_filename) // '"'
  call read_list_header(trim(dm_halo_list_filename), unit, ndm_halo, dm_halo_cols)
  allocate(data_halo(ndm_halo, dm_halo_cols))
  call read_list_data(unit, ndm_halo, dm_halo_cols, data_halo)
  print*, 'Got', ndm_halo, 'halos in da pocket'

  !-------------------------------------
  ! Read brick file
  !-------------------------------------
  print*, ''
  print*, 'Reading brick file…'
  call read_info_headers(info_file, infos)
  call read_brick_header(brick_file, infos, nbodies, aexp_tmp, age_univ,&
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
  allocate(halo_to_cpu(nDM, NCPU_PER_HALO))
  call read_brick_data(nDM, infos, .true., &
       & mDM, posDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LDM, idDM, members)
  print*, '    …red!'

  halo_to_cpu = 0
  tmp_int2 = 0

  if (param_to < 0) then
     param_to = infos%ncpu + param_to + 1
  end if

  !allocate(particles_found(nDM))
  !$OMP PARALLEL DO PRIVATE(halo_found, tmp_char, i, j, tmp_real, tmp_int, ndim, nparts, unit) &
  !$OMP PRIVATE(order, pos, vel, ids, m, birth_date, tmp_arr, nparticle_to_probe_halo) &
  !$OMP SCHEDULE(guided, 5)
  do cpu = param_from, param_to
     write(tmp_char, '(i0.5)') cpu
     tmp_char = "/data52/Horizon-AGN/OUTPUT_DIR/output_00002/part_00002.out" // trim(tmp_char)

     ! call read_particle(param_output_path, param_output_number, cpu, nstar, pos, vel, m,&
     !      ids, birth_date, ndim, nparts)
     call read_particle_header(trim(tmp_char), ndim, nparts, unit)
     allocate(pos(ndim, nparts), vel(ndim, nparts), ids(nparts), m(nparts), birth_date(nparts))
     call read_particle_data(ndim, nparts, unit, nstar, pos, vel, m, ids, birth_date)
     deallocate(pos, vel, m, birth_date)

     allocate(order(nparts))
     call quick_sort(ids, order)
     deallocate(order)

     halo_found = 0
     do i = 1, nDM
        if (mDM(i) < param_min_m .or. mDM(i) > param_max_m) then
           cycle
        end if
        nparticle_to_probe_halo = max(50, &
             floor(frac_particle_to_probe * members(i)%parts / 100.))
        allocate(tmp_arr(nparticle_to_probe_halo))

        ! Pick 10 random particles and see if it's in the CPU
        do j = 1, nparticle_to_probe_halo
           call random_number(tmp_real)
           tmp_int = ceiling(tmp_real*members(i)%parts)
           ! get the position of the random particle in the ids
           tmp_arr(j) = indexOf(members(i)%ids(tmp_int), ids)
           if (tmp_arr(j) > 0) then
              exit
           end if
        end do
        do j = 1, nparticle_to_probe_halo
           if (tmp_arr(j) > 0) then
              call fill(halo_to_cpu(i, :), cpu, halo_found)
              exit
           end if
        end do

        deallocate(tmp_arr)
     end do
     !$OMP ATOMIC
     tmp_int2 = tmp_int2 + 1
     write(*, '(i4,  a,  i4,         a,          i4,  a,x,        i5,x,       a)') &
          tmp_int2, '/', infos%ncpu, ' (cpu n°', cpu ,'), ', halo_found, 'halos found'
     deallocate(ids)
  end do
  !$OMP END PARALLEL DO

  ! !-------------------------------------
  ! ! Find missing particles
  ! !-------------------------------------
  ! do i = 1, nDM
  !    if (particles_found(i) /= members(i)%parts) then
  !       ! Missing particle :(
  !    end if
  ! end do

  !-------------------------------------
  ! Simplify columns (find first empty column)
  !-------------------------------------
  tmp_int = NCPU_PER_HALO
  allocate(tmp_arr(NCPU_PER_HALO))
  tmp_arr = maxval(halo_to_cpu, 1)
  do i = 1, NCPU_PER_HALO
     if (tmp_arr(i) == 0) then
        tmp_int = i-1
        exit
     end if
  end do

  print*, 'Found ', tmp_int, 'empty columns'
  !-------------------------------------
  ! Only save non empty columns
  !-------------------------------------
  write(tmp_char, '(a, i0.5, a, i0.5)') "out", param_from, '-', param_to
  open(10, file=tmp_char)
  write(10, *) nDM, tmp_int
  do i = 1, nDM
     tmp_arr = halo_to_cpu(i, :)
     write(10, '(i10,100i6)') idDM(i), tmp_arr(:tmp_int)
  end do
  close(10)
  deallocate(tmp_arr)

  call write_list(trim(tmp_char) // '.bin', nDM, tmp_int, halo_to_cpu)

contains
  subroutine count_particles(halo_ids, ids, count)
    integer, intent(in), dimension(:) :: halo_ids, ids

    integer, intent(out) :: count

    do i = 1, size(halo_ids)
       if (indexOf(halo_ids(i), ids) > 0) then
          count = count + 1
       end if
    end do
  end subroutine count_particles
  ! subroutine read_params ()
  !   integer            :: i,n
  !   integer            :: iargc
  !   character(len=4)   :: opt
  !   character(len=128) :: arg
  !   n = iargc()

  !   ! if (n < 4) then
  !   !    print*, 'Sort the galaxies between elliptics and spirals'
  !   !    print*, ''
  !   !    print*, 'Usage (all lists should come in binary format):'
  !   !    print*, '\t -fga File to pick galaxies from'
  !   !    print*, '\t -fdm File to pick dark matter halos from'
  !   !    print*, '\t -fas File to get associations between dark matter halos and galaxies'
  !   !    !TODO: write usage
  !   ! end if

  !   !-------------------------------------
  !   ! default values
  !   !-------------------------------------
  !   gal_list_filename = "lists/list_kingal_00782.dat"
  !   dm_halo_list_filename = "lists/list_halo.dat.bin"
  !   associations_filename = "lists/associated_halogal_782.dat.bin"
  !   ramses_output_start = "/data52/Horizon-AGN/OUTPUT_DIR/output_00002"
  !   ramses_output_end = "/data52/Horizon-AGN/OUTPUT_DIR/output_00782"
  !   brick_file = "/data52/Horizon-AGN/TREE_DM_celldx2kpc_SC0.9r/tree_bricks782"
  !   info_file = '/data52/Horizon-AGN/OUTPUT_DIR/output_00782/info_00782.txt'
  !   info_file_start = '/data52/Horizon-AGN/OUTPUT_DIR/output_00002/info_00002.txt'
  !   halo_to_cpu_file = 'lists/halo_to_cpu.00002.raw.dat.bin'
  !   mergertree_file = '/data33/dubois/H-AGN/MergerTree/TreeMaker_HAGN/tree.dat'

  !   param_from = 1
  !   param_to = -1

  !   do i = 1, n, 2
  !      call getarg(i, opt)
  !      if (i == n) then
  !         print*, '("option ",a2," has no argument")', opt
  !         stop
  !      end if

  !      call getarg(i+1, arg)
  !      select case(opt)
  !      case ('-fga')
  !         gal_list_filename = trim(arg)
  !      case ('-fdm')
  !         dm_halo_list_filename = trim(arg)
  !      case('-fas')
  !         associations_filename = trim(arg)
  !      case('-ods') ! Ramses output start
  !         ramses_output_start = trim(arg)
  !      case('-ode') ! Ramses output start
  !         ramses_output_end = trim(arg)
  !      case('-bri')
  !         brick_file = trim(arg)
  !      case('-ifs')
  !         info_file_start = trim(arg)
  !      case('-ife')
  !         info_file = trim(arg)
  !      case ('-out')
  !         outfile = trim(arg)
  !      case ('-fro')
  !         tmp_char = trim(arg)
  !         read(tmp_char, '(i10)') param_from
  !      case ('-to')
  !         tmp_char = trim(arg)
  !         read(tmp_char, '(i10)') param_to
  !      case default
  !         print '("unknown option ",a2," ignored")', opt
  !      end select
  !   end do

  ! end subroutine read_params

end program sort_galaxy
