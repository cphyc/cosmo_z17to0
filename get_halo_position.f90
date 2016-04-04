program sort_galaxy
  use io
  use misc
  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: ramses_output_end, ramses_output_start, gal_list_filename, associations_filename, dm_halo_list_filename, brick_file, info_file_end, info_file_start, outfile, halo_to_cpu_file, mergertree_file
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
  integer, dimension(:, :), allocatable :: halo_to_cpu
  real(kind=8), dimension(:), allocatable    :: m, birth_date
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer            :: i, cpu
  integer            :: tmp_int
  character(len=200) :: tmp_char
  real               :: tmp_real
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
  allocate(data_gal(ngal, gal_cols+1))
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
  ! print*, 'Reading brick file…'
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
  allocate(halo_to_cpu(nDM, 3))
  call read_brick_data(nDM, infos, .true., &
       & mDM, posDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LDM, idDM, members)
  print*, '    …red!'

  halo_to_cpu = 0
  do cpu = 1, infos%ncpu
     write(tmp_char, '(i0.5)') cpu
     tmp_char = "/data52/Horizon-AGN/OUTPUT_DIR/output_00002/part_00002.out" // trim(tmp_char)

     call read_particle_header(trim(tmp_char), ndim, nparts)
     allocate(pos(ndim, nparts), vel(ndim, nparts), ids(nparts), m(nparts), birth_date(nparts),&
          order(nparts))
     call read_particle_data(ndim, nparts, nstar, pos, vel, m, ids, birth_date)

     call quick_sort(ids, order, nparts)

     halo_found = 0
     !!$OMP PARALLEL DO
     do i = 1, nDM

        call random_number(tmp_real)
        tmp_int = floor(tmp_real*members(i)%parts)
        tmp_int = indexOf(members(i)%ids(tmp_int), ids)
        if (halo_to_cpu(i, 3) > 0) then
           cycle
        else if (tmp_int > 0) then
           if (halo_to_cpu(i, 1) > 0) then
              if (halo_to_cpu(i, 2) > 0) then
                 halo_to_cpu(i, 3) = cpu
              else
                 halo_to_cpu(i, 2) = cpu
              end if
           else
              halo_to_cpu(i, 1) = cpu
           end if
           halo_found = halo_found + 1
        end if
     end do
!!$OMP END PARALLEL DO
     write(*, '(a3,x,i4,a1,i4,x,a,i5,x,a)') 'cpu', cpu, '/', infos%ncpu, ', found', halo_found, 'halos'

     deallocate(pos, vel, ids, m, birth_date, order)
  end do

  open(10, file='out')
  do i = 1, nDM
     write(10, *) idDM(i), halo_to_cpu(i, :)
  end do
  close(10)

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
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

  end subroutine read_params

end program sort_galaxy
