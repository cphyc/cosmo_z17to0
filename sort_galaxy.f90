program sort_galaxy
  use io
  use misc
  implicit none

  character(LEN=128) :: gal_list_filename, repository, outfile, fileinfo='none', halogalfile
  character(LEN=128) :: dm_halo_list_filename, dm_gal_assoc, associations_filename
  logical :: verbose
  integer :: ngal, gal_cols, i
  integer :: ndm_halo, dm_halo_cols
  integer :: nassoc, assoc_cols
  integer :: ellipticals = 0, spirals = 0, others = 0
  real(kind=4), dimension(:,:), allocatable :: data_gal, data_halo
  real(kind=4), dimension(:, :), allocatable :: data_associations ! contains association flag / dark matter id / mass / gal id / mass
  real(kind=4) :: dmhalo_id, lvl, dmhalo_mass, gal_id, gal_mass
  real(kind=4), dimension(:), allocatable :: kind ! 0 for elliptical, 1 for spirals, 2 for other
  real(kind=4) :: sigma
  real(kind=4) :: elliptical_threshold = 1.5d0 ! FIXME
  real(kind=4) :: spiral_threshold = 0.8d0 ! FIXME

  integer :: gal_num, halo_num, unit
  integer, dimension(:), allocatable :: gal_to_halo

  !-------------------------------------
  ! Read parameters
  !-------------------------------------
  call read_params

  !-------------------------------------
  ! Read file
  !-------------------------------------
  if (gal_list_filename == 'none') then
     stop
  end if

  print*, 'Reading file "' // trim(gal_list_filename) // '"'
  call read_list_header(trim(gal_list_filename), ngal, gal_cols, unit)
  allocate(data_gal(ngal, gal_cols))
  call read_list_data(ngal, gal_cols, data_gal, unit)
  print*, 'Got', nassoc, 'galaxies in da pocket'

  print*, ''
  print*, 'Reading file "' // trim(associations_filename) // '"'
  call read_list_header(trim(associations_filename), nassoc, assoc_cols, unit)
  allocate(data_associations(nassoc, assoc_cols))
  call read_list_data(nassoc, assoc_cols, data_associations, unit)
  print*, 'Got', nassoc, 'associations in da pocket'

  print*, ''
  print*, 'Reading file "' // trim(dm_halo_list_filename) // '"'
  call read_list_header(trim(dm_halo_list_filename), ndm_halo, dm_halo_cols, unit)
  allocate(data_halo(ndm_halo, dm_halo_cols))
  call read_list_data(ndm_halo, dm_halo_cols, data_halo, unit)
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
        kind(i) = 2
     end if

  end do

  print*, 'Associating galaxies to halo…'
  allocate(gal_to_halo(int(maxval(data_gal(:, 1)))))
  gal_to_halo = -1
  do i = 1, nassoc
     gal_num = int(data_associations(i, 4))
     halo_num = int(data_associations(i, 1))

     gal_to_halo(gal_num) = halo_num
  end do

  print*, 'Some statistics…'
  ! Sort the galaxies
  print*, int(ellipticals*100./ngal), "% ellipticals"
  print*, int(spirals*100./ngal), "% spirals"
  print*, int(100 - (ellipticals + spirals)*100./ngal), "% others"

  !-------------------------------------
  ! Test
  !-------------------------------------
  ! FIXME: do it for each galaxy!
  do i = 1, ngal
     ! Find halo related to galaxy in association table
     gal_num = int(data_gal(i, 1))
     ! Find particles in halo

     ! Loop back in time
  end do

  !-------------------------------------
  ! Cleanup
  !-------------------------------------
  close(11)
  close(12)
  close(13)
  deallocate(data_gal)
  deallocate(data_halo)
  deallocate(data_associations)
  deallocate(gal_to_halo)

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
       case ('-out')
          outfile = trim(arg)
       case ('-fin')
          fileinfo = trim(arg)
       case ('-halogal')
          halogalfile = trim(arg)
       case ('-v')
          read (arg,*) verbose
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

  end subroutine read_params

end program sort_galaxy
