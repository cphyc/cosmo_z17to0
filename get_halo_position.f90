program sort_galaxy
  use io
  use misc
  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: ramses_output_end, ramses_output_start, gal_list_filename, associations_filename, dm_halo_list_filename, brick_file, info_file_end, info_file_start, outfile, halo_to_cpu_file, mergertree_file
  integer :: param_from, param_to
  integer, parameter :: NPARTICLE_TO_PROBE_HALO = 50, NCPU_PER_HALO = 10
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
  ! integer                                    :: initial_halo_id, current_halo, parent_halo, istep, prev, father, step, max_nhalo
  ! integer, dimension(:), allocatable         :: tmp_nhalos, halos_z0
  ! integer, dimension(:), allocatable         :: parent, parent_at_step
  !-------------------------------------
  ! Particle data
  !-------------------------------------
  integer                                    :: ndim, nparts
  real(kind=8), dimension(:, :), allocatable :: pos, vel
  integer                                    :: nstar, halo_found
  integer,      dimension(:), allocatable    :: ids, order
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
  call read_params()

  !-------------------------------------
  ! Read lists
  !-------------------------------------
  print*, 'Reading file "' // trim(gal_list_filename) // '"'
  call read_list_header(trim(gal_list_filename), ngal, gal_cols)
  allocate(data_gal(ngal, gal_cols))
  call read_list_data(ngal, gal_cols, data_gal)
  print*, 'Got', ngal, 'galaxies in da pocket'

  print*, ''
  print*, 'Reading file "' // trim(associations_filename) // '"'
  call read_list_header(trim(associations_filename), nassoc, assoc_cols)
  allocate(data_associations(nassoc, assoc_cols))
  call read_list_data(nassoc, assoc_cols, data_associations)
  print*, 'Got', nassoc, 'associations in da pocket'

  print*, ''
  print*, 'Reading file "' // trim(dm_halo_list_filename) // '"'
  call read_list_header(trim(dm_halo_list_filename), ndm_halo, dm_halo_cols)
  allocate(data_halo(ndm_halo, dm_halo_cols))
  call read_list_data(ndm_halo, dm_halo_cols, data_halo)
  print*, 'Got', ndm_halo, 'halos in da pocket'

  !-------------------------------------
  ! Read brick file
  !-------------------------------------
  print*, ''
  print*, 'Reading brick file…'
  call read_info_headers(info_file_end, infos)
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

  allocate(tmp_arr(NPARTICLE_TO_PROBE_HALO))
  halo_to_cpu = 0
  tmp_int2 = 0

  if (param_to < 0) then
     param_to = infos%ncpu + param_to + 1
  end if

  !$OMP PARALLEL DO PRIVATE(halo_found, tmp_char, i, j, tmp_real, tmp_int, ndim, nparts, unit) &
  !$OMP PRIVATE(order, pos, vel, ids, m, birth_date, tmp_arr) SCHEDULE(guided, 1)
  do cpu = param_from, param_to
     write(tmp_char, '(i0.5)') cpu
     tmp_char = "/data52/Horizon-AGN/OUTPUT_DIR/output_00002/part_00002.out" // trim(tmp_char)

     call read_particle_header(trim(tmp_char), ndim, nparts, unit)
     allocate(pos(ndim, nparts), vel(ndim, nparts), ids(nparts), m(nparts), birth_date(nparts))
     call read_particle_data(ndim, nparts, unit, nstar, pos, vel, m, ids, birth_date)
     deallocate(pos, vel, m, birth_date)

     allocate(order(nparts))
     call quick_sort(ids, order, nparts)
     deallocate(order)

     halo_found = 0

     do i = 1, nDM
        ! Pick 10 random particles and see if it's in the CPU
        do j = 1, NPARTICLE_TO_PROBE_HALO
           call random_number(tmp_real)
           tmp_int = ceiling(tmp_real*members(i)%parts)
           ! get the position of the random particle in the ids
           tmp_arr(j) = indexOf(members(i)%ids(tmp_int), ids)
           if (tmp_arr(j) > 0) then
              exit
           end if
        end do
        do j = 1, NPARTICLE_TO_PROBE_HALO
           if (tmp_arr(j) > 0) then
              call fill(halo_to_cpu(i, :), cpu, halo_found)
              exit
           end if
        end do
     end do
     !$OMP ATOMIC
     tmp_int2 = tmp_int2 + 1
     write(*, '(i4,  a,  i4,         a,          i4,  a,x,        i5,x,       a)') &
          tmp_int2, '/', infos%ncpu, ' (cpu n°', cpu ,'), ', halo_found, 'halos found'
     deallocate(ids)
  end do
  !$OMP END PARALLEL DO
  deallocate(tmp_arr)

  allocate(tmp_arr(NCPU_PER_HALO))
  write(tmp_char, '(a, i0.5, a, i0.5)') "out", param_from, '-', param_to
  open(10, file=tmp_char)
  write(10, *) nDM, NCPU_PER_HALO
  do i = 1, nDM
     tmp_arr = halo_to_cpu(i, :)
     write(10, *) idDM(i), tmp_arr
  end do
  close(10)
  deallocate(tmp_arr)

  call write_list(trim(tmp_char) // '.bin', nDM, NCPU_PER_HALO, halo_to_cpu)

contains
  subroutine read_params ()
    integer            :: i,n
    integer            :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    n = iargc()

    ! if (n < 4) then
    !    print*, 'Sort the galaxies between elliptics and spirals'
    !    print*, ''
    !    print*, 'Usage (all lists should come in binary format):'
    !    print*, '\t -fga File to pick galaxies from'
    !    print*, '\t -fdm File to pick dark matter halos from'
    !    print*, '\t -fas File to get associations between dark matter halos and galaxies'
    !    !TODO: write usage
    ! end if

    !-------------------------------------
    ! default values
    !-------------------------------------
    gal_list_filename = "lists/list_kingal_00782.dat"
    dm_halo_list_filename = "lists/list_halo.dat.bin"
    associations_filename = "lists/associated_halogal_782.dat.bin"
    ramses_output_start = "/data52/Horizon-AGN/OUTPUT_DIR/output_00002"
    ramses_output_end = "/data52/Horizon-AGN/OUTPUT_DIR/output_00782"
    brick_file = "/data52/Horizon-AGN/TREE_DM_celldx2kpc_SC0.9r/tree_bricks782"
    info_file_end = '/data52/Horizon-AGN/OUTPUT_DIR/output_00782/info_00782.txt'
    info_file_start = '/data52/Horizon-AGN/OUTPUT_DIR/output_00002/info_00002.txt'
    halo_to_cpu_file = 'lists/halo_to_cpu.00002.raw.dat.bin'
    mergertree_file = '/data33/dubois/H-AGN/MergerTree/TreeMaker_HAGN/tree.dat'

    param_from = 1
    param_to = -1

    do i = 1, n, 2
       call getarg(i, opt)
       if (i == n) then
          print*, '("option ",a2," has no argument")', opt
          stop
       end if

       call getarg(i+1, arg)
       select case(opt)
       case ('-fga')
          gal_list_filename = trim(arg)
       case ('-fdm')
          dm_halo_list_filename = trim(arg)
       case('-fas')
          associations_filename = trim(arg)
       case('-ods') ! Ramses output start
          ramses_output_start = trim(arg)
       case('-ode') ! Ramses output start
          ramses_output_end = trim(arg)
       case('-bri')
          brick_file = trim(arg)
       case('-ifs')
          info_file_start = trim(arg)
       case('-ife')
          info_file_end = trim(arg)
       case ('-out')
          outfile = trim(arg)
       case ('-fro')
          tmp_char = trim(arg)
          read(tmp_char, '(i10)') param_from
       case ('-to')
          tmp_char = trim(arg)
          read(tmp_char, '(i10)') param_to
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

  end subroutine read_params

end program sort_galaxy