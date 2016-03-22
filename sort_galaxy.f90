program sort_galaxy
  use io
  use misc
  implicit none

  character(LEN=128)                         :: gal_list_filename, outfile, fileinfo='none', filename
  character(LEN=128)                         :: dm_halo_list_filename, associations_filename
  character(LEN=128)                         :: ramses_output_end, ramses_output_start, brick_file, info_file_start, info_file_end
  logical                                    :: verbose
  integer                                    :: ngal, gal_cols, i, j, k, index
  integer                                    :: ndm_halo, dm_halo_cols
  integer                                    :: nassoc, assoc_cols
  integer                                    :: ellipticals = 0, spirals = 0, others = 0
  real(kind=4), dimension(:,:), allocatable  :: data_gal, data_halo
  ! contains association flag / dark matter id / mass / gal id / mass
  real(kind=4), dimension(:, :), allocatable :: data_associations
  real(kind=4), dimension(:), allocatable    :: kind ! 0 for elliptical, 1 for spirals, 2 for other
  real(kind=4)                               :: sigma
  real(kind=4)                               :: elliptical_threshold = 1.5d0 ! FIXME
  real(kind=4)                               :: spiral_threshold = 0.8d0 ! FIXME

  real(kind=8), dimension(:), allocatable    :: mDM, rvirDM
  real(kind=8), dimension(:), allocatable    :: mvirDM, TvirDM, hlevel
  real(kind=8), dimension(:, :), allocatable :: LDM, posDM
  integer, dimension(:), allocatable         :: idDM, cpu_list

  real(kind=4):: aexp_tmp, age_univ
  real(kind=8), dimension(:), allocatable :: X0, X1
  integer :: nbodies, nb_of_halos, nb_of_subhalos, nDM

  integer                                    :: gal_id, halo_id, id_in_brick, cpu
  character(len=5) :: cpu_as_str
  integer, dimension(:), allocatable         :: gal_to_halo

  type(MEMBERS_T), dimension(:), allocatable :: members
  type(INFOS_T) :: infos_start, infos_end

  integer :: ndim, nparts, nstars, particles_to_find
  real(kind=8), dimension(:,:), allocatable :: pos, vel, tmp_pos, tmp_vel
  real(kind=8), dimension(:), allocatable :: tmp_m, tmp_birth_date
  integer, dimension(:), allocatable :: tmp_ids

  !-------------------------------------
  ! Read parameters
  !-------------------------------------
  call read_params()

  !-------------------------------------
  ! Read file
  !-------------------------------------
  if (gal_list_filename == 'none') then
     stop
  end if

  print*, 'Reading file "' // trim(gal_list_filename) // '"'
  call read_list_header(trim(gal_list_filename), ngal, gal_cols)
  allocate(data_gal(ngal, gal_cols+1))
  call read_list_data(ngal, gal_cols, data_gal(1:ngal, 1:gal_cols))
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

  allocate(kind(ngal))

  !-------------------------------------
  ! Processing
  !-------------------------------------
  print*, 'Sorting galaxies (ellipticals/spirals) using v/sigma criterion…'
  do i = 1, ngal
     sigma = 1d0/3d0 * sqrt(data_gal(i, 3)**2 + data_gal(i, 4)**2 + data_gal(i, 5)**2)
     if (sigma / data_gal(i, 2) > elliptical_threshold) then
        data_gal(i, gal_cols+1) = 0
        ellipticals = ellipticals + 1
     else if (sigma / data_gal(i, 2) < spiral_threshold) then
        data_gal(i, gal_cols+1) = 0
        spirals = spirals + 1
     else
        others = others + 1
        kind(i) = 2
     end if

  end do

  print*, 'Associating galaxies to halo…'
  allocate(gal_to_halo(int(maxval(data_gal(:, 1))))) ! this is the ids
  gal_to_halo = -1
  do i = 1, nassoc
     gal_id = int(data_associations(i, 4))
     halo_id = int(data_associations(i, 1))
     if (gal_id > 0) then
        gal_to_halo(gal_id) = halo_id
     end if
  end do
  i = -1

  print*, 'Some statistics…'
  ! Sort the galaxies
  print*, int(ellipticals*100./ngal), "% ellipticals"
  print*, int(spirals*100./ngal), "% spirals"
  print*, int(100 - (ellipticals + spirals)*100./ngal), "% others"

  !-------------------------------------
  ! Read brick
  !-------------------------------------
  call read_info_headers(info_file_start, infos_start)
  call read_info_headers(info_file_start, infos_end)
  print*, "Reading brick file " // brick_file
  call read_brick_header(brick_file, infos_start, nbodies, aexp_tmp, age_univ,&
       nb_of_halos, nb_of_subhalos)
  nDM = nb_of_halos + nb_of_subhalos
  allocate(idDM(nDM))
  allocate(posDM(infos_start%ndim, nDM))
  allocate(rvirDM(nDM))
  allocate(mDM(nDM))
  allocate(mvirDM(nDM))
  allocate(TvirDM(nDM))
  allocate(hlevel(nDM))
  allocate(LDM(infos_start%ndim, nDM))
  allocate(members(nDM))
  print*, infos_start%boxlen, infos_start%unit_l/3.085677581e+24,  infos_start%boxlen*infos_start%unit_l/3.085677581e+24
  call read_brick_data(nDM, infos_start, .true., &
       & mDM, posDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LDM, idDM, members)

  print*, '', 'read!'
  ! print*, maxval(posDM(1,:)), maxval(posDM(2,:)), maxval(posDM(3,:))
  ! print*, minval(posDM(1,:)), minval(posDM(2,:)), minval(posDM(3,:))
  !-------------------------------------
  ! Iterate over each galaxy
  !-------------------------------------
  allocate(cpu_list(infos_start%ncpu))
  allocate(X0(infos_start%ndim))
  allocate(X1(infos_start%ndim))

  do gal_id = 1, 1 !ngal
     halo_id = gal_to_halo(gal_id)
     print*, 'Doing galaxy', gal_id, 'halo', halo_id
     if (halo_id <= 0) then
        cycle
     end if

     !-------------------------------------
     ! find particles in halo
     !-------------------------------------
     do id_in_brick = 1, nDM
        if (idDM(id_in_brick) == halo_id) then
           exit
        else
           cycle
        end if
     end do
     X0 = posDM(:, id_in_brick) - rvirDM(id_in_brick)
     X1 = posDM(:, id_in_brick) + rvirDM(id_in_brick)
     print*, X0, X1
     call get_cpu_list(X0, X1, infos_start%levelmax, infos_start%bound_key, &
          cpu_list, infos_start%ncpu, infos_start%ndim)

     !-------------------------------------
     ! Number of particles to read
     !-------------------------------------
     allocate(pos(infos_start%ndim, members(id_in_brick)%parts))
     allocate(vel(infos_start%ndim, members(id_in_brick)%parts))
     pos = 0
     vel = 0
     !-------------------------------------
     ! Reading outputs
     !-------------------------------------
     call quick_sort(members(halo_id)%ids, members(halo_id)%parts)

     print*, 'Looking for'
     print*, members(halo_id)%ids
     do j = 3277,3277!1, infos_start%ncpu
        cpu = 3277!cpu_list(j)
        if (cpu <= 0) then
           exit
        end if

        !-------------------------------------
        ! reading output
        !-------------------------------------
        write(cpu_as_str, '(i0.4)') cpu
        filename = trim(ramses_output_start) // '/part_00032.out0' // trim(cpu_as_str)

        print*, 'Reading', filename

        call read_particle_header(filename, ndim, nparts)
        allocate(tmp_pos(ndim, nparts))
        allocate(tmp_vel(ndim, nparts))
        allocate(tmp_m(nparts))
        allocate(tmp_ids(nparts))
        allocate(tmp_birth_date(nparts))
        call read_particle_data(ndim, nparts, nstars, tmp_pos, tmp_vel, tmp_m,&
             tmp_ids, tmp_birth_date)
        !-------------------------------------
        ! find particles in halo
        !-------------------------------------
        particles_to_find = members(halo_id)%parts - 1
        do k = 1, nparts
           index = indexOf(tmp_ids(k), members(halo_id)%ids)
           if (index > 0) then
              particles_to_find = particles_to_find - 1
              ! Store the data
              pos(:, index) = tmp_pos(:, k)
              vel(:, index) = tmp_vel(:, k)
              if (particles_to_find < 0) then
                 ! No need to read more, we found all particles
                 print*, 'Found all particles!'
                 exit
              end if
           end if
        end do
        open(unit=11, file="cpu3246", form='unformatted')
        write(11) tmp_ids
        close(11)
        open(unit=11, file="part_list", form='unformatted')
        write(11) members(halo_id)%ids
        close(11)
        deallocate(tmp_pos, tmp_vel, tmp_m, tmp_ids, tmp_birth_date)
     end do
     print*, pos
     deallocate(pos, vel)
  end do
  deallocate(X0)
  deallocate(X1)

  !-------------------------------------
  ! Cleanup
  !-------------------------------------
  deallocate(cpu_list)
  deallocate(data_gal)
  deallocate(data_halo)
  deallocate(data_associations)
  deallocate(gal_to_halo)
  deallocate(idDM)
  deallocate(posDM)
  deallocate(rvirDM)
  deallocate(mDM)
  deallocate(mvirDM)
  deallocate(TvirDM)
  deallocate(hlevel)
  deallocate(LDM)
  do i = 1, nDM
     deallocate(members(i)%ids)
  end do
  deallocate(members)

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
    brick_file = "/data40b/Horizon-AGN/TREE_DM_celldx2kpc_Rmax_guess9/tree_bricks761"
    info_file_end = '/data52/Horizon-AGN/OUTPUT_DIR/output_00761/info_00761.txt'
    info_file_start = '/data52/Horizon-AGN/OUTPUT_DIR/output_00032/info_00032.txt'

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
       case ('-fin')
          fileinfo = trim(arg)
       case ('-v')
          read (arg,*) verbose
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

  end subroutine read_params

end program sort_galaxy
