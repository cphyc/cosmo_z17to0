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
  integer :: ngal, nassoc, ndm_halo, gal_cols, assoc_cols, dm_halo_cols
  !-------------------------------------
  ! Brick data
  !-------------------------------------
  real(kind=8), dimension(:, :), allocatable :: posDM, LDM
  integer, dimension(:), allocatable :: idDM
  real(kind=8), dimension(:), allocatable :: rvirDM, mDM, mvirDM, TvirDM, hlevel
  type(INFOS_T) :: infos_start
  type(MEMBERS_T), dimension(:), allocatable :: members
  integer :: nbodies, nb_of_halos, nb_of_subhalos, nDM
  real(kind=4) :: aexp_tmp, age_univ
  !-------------------------------------
  ! Merger tree data
  !-------------------------------------
  integer,  allocatable, dimension(:) :: mt_nhalos, mt_nsubhalos
  real(kind=4), dimension(:), allocatable  :: mt_aexp, mt_omega_t, mt_age_univ
  integer(kind=4) :: mt_nsteps
  integer :: initial_halo_id, current_halo, parent_halo, istep, prev, father, step, max_nhalo
  integer, dimension(:), allocatable :: tmp_nhalos, halos_z0
  integer, dimension(:,:), allocatable :: parent, progenitor
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer :: i

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
  ! call read_info_headers(info_file_start, infos_start)
  ! call read_brick_header(brick_file, infos_start, nbodies, aexp_tmp, age_univ,&
  !      nb_of_halos, nb_of_subhalos)
  ! nDM = nb_of_halos + nb_of_subhalos
  ! allocate(idDM(nDM))
  ! allocate(posDM(infos_start%ndim, nDM))
  ! allocate(rvirDM(nDM))
  ! allocate(mDM(nDM))
  ! allocate(mvirDM(nDM))
  ! allocate(TvirDM(nDM))
  ! allocate(hlevel(nDM))
  ! allocate(LDM(infos_start%ndim, nDM))
  ! allocate(members(nDM))
  ! call read_brick_data(nDM, infos_start, .true., &
  !      & mDM, posDM, rvirDM, mvirDM, TvirDM,&
  !      & hlevel, LDM, idDM, members)
  ! print*, '    …red!'

  !-------------------------------------
  ! Read merger tree and reduce it
  !-------------------------------------
  print*, ''
  print*, 'Reading merger tree file "' // trim(mergertree_file) // '"…'
  call read_mergertree_headers_1(mergertree_file, mt_nsteps)
  allocate(mt_nhalos(mt_nsteps), mt_nsubhalos(mt_nsteps), mt_aexp(mt_nsteps),&
       mt_omega_t(mt_nsteps), mt_age_univ(mt_nsteps), tmp_nhalos(mt_nsteps))
  call read_mergertree_headers_2(mt_nhalos, mt_nsubhalos, mt_aexp, mt_omega_t, mt_age_univ,&
       mt_nsteps)
  tmp_nhalos = mt_nhalos + mt_nsubhalos
  max_nhalo = maxval(tmp_nhalos, 1)
  allocate(parent(max_nhalo, mt_nsteps), halos_z0(max_nhalo))
  call read_mergertree_parent_of(max_nhalo, tmp_nhalos, mt_nsteps, halos_z0, parent)
  print*, '      …red!'
  !-------------------------------------
  ! reduce the tree
  !-------------------------------------
  allocate(progenitor(max_nhalo, 2))
  print*, tmp_nhalos(mt_nsteps)
  do i = 1, tmp_nhalos(mt_nsteps)
     prev = i
     father = parent(i, mt_nsteps)
     step = mt_nsteps - 1
     do while (father /= 0 .and. step > 0)
        prev = father
        father = parent(father, step)
        step = step - 1
     end do
     progenitor(i, 1) = prev
     progenitor(i, 2) = step + 1
  end do

  open(unit=10, file='out')
  print*, 'Saving file'
  write(10, '(a, i)') 'nsteps', mt_nsteps
  write(10, '(a, i)') 'nhalos', max_nhalo
  write(10, '(a, i)') 'dimensions', 3
  write(10, '(3a10)') 'halo', 'parent', 'step'
  do i = 1, tmp_nhalos(mt_nsteps)
     write(10, '(3i10)') halos_z0(i), progenitor(i, 1), progenitor(i, 2)
  end do
  close(10)
  open(unit=10, file="out.raw", form="unformatted")
  print*, "Saving raw file"
  write(10) mt_nsteps, max_nhalo, 3
  write(10) halos_z0, progenitor(:, 1), progenitor(:, 2)
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
    ramses_output_start = "/data52/Horizon-AGN/OUTPUT_DIR/output_00032"
    ramses_output_end = "/data52/Horizon-AGN/OUTPUT_DIR/output_00761"
    brick_file = "/data52/Horizon-AGN/TREE_DM_celldx2kpc_SC0.9r/tree_bricks782"
    info_file_end = '/data52/Horizon-AGN/OUTPUT_DIR/output_00761/info_00761.txt'
    info_file_start = '/data52/Horizon-AGN/OUTPUT_DIR/output_00032/info_00032.txt'
    halo_to_cpu_file = 'lists/halo_to_cpu.00032.raw.dat.bin'
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
