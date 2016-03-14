program sort_galaxy
  use hashtbl
  use convert
  implicit none
  character(LEN=128) :: gal_list_filename, repository, outfile, fileinfo='none', halogalfile
  character(LEN=128) :: dm_halo_list_filename, dm_gal_assoc, associations_filename
  logical :: verbose
  integer :: ngal, gal_cols, i
  integer :: ndmhalo, dm_halo_cols
  integer :: nassoc, assoc_cols
  integer :: ellipticals = 0, spirals = 0, others = 0
  real(kind=4), dimension(:,:), allocatable :: data_gal, data_halo
  real(kind=4), dimension(:, :), allocatable :: associations ! contains association flag / dark matter id / mass / gal id / mass
  real(kind=4) :: dmhalo_id, lvl, dmhalo_mass, gal_id, gal_mass
  real(kind=4), dimension(:), allocatable :: kind ! 0 for elliptical, 1 for spirals, 2 for other
  real(kind=4) :: sigma
  real(kind=4) :: elliptical_threshold = 1.5d0 ! FIXME
  real(kind=4) :: spiral_threshold = 0.8d0 ! FIXME

  character(LEN=20) :: gal_num, halo_num
  type(hash_tbl_sll) :: gal_to_halo
  character(LEN=:), allocatable :: hashmap_char

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
  open(unit=11, file=gal_list_filename, form='unformatted')
  read(11) ngal, gal_cols
  print*, 'Got', ngal, 'galaxies in da pocket'

  print*, ''
  print*, 'Reading file "' // trim(associations_filename) // '"'
  open(unit=12, file=associations_filename, form='unformatted')
  read(12) nassoc, assoc_cols
  print*, 'Got', nassoc, 'associations in da pocket'

  print*, ''
  print*, 'Reading file "' // trim(dm_halo_list_filename) // '"'
  open(unit=13, file=dm_halo_list_filename, form='unformatted')
  read(13) ndmhalo, dm_halo_cols
  print*, 'Got', ndmhalo, 'halos in da pocket'

  !-------------------------------------
  ! Allocations & datareading
  !-------------------------------------
  allocate(data_gal(ngal, gal_cols+1))
  allocate(data_halo(ndmhalo, dm_halo_cols+1))
  allocate(associations(nassoc, assoc_cols+1))

  allocate(kind(ngal))

  read(11) data_gal(1:ngal, 1:gal_cols)
  read(12) associations(1:nassoc, 1:assoc_cols)
  read(13) data_halo(1:ndmhalo, 1:dm_halo_cols)

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
  call gal_to_halo%init(nassoc)
  do i = 1, nassoc
     if (mod(i, nassoc/25) == 0) then
        print*, int(dble(i) / nassoc * 100), '%'
     end if
     gal_num = itos(int(associations(i, 4)))
     halo_num = itos(int(associations(i, 1)))

     call gal_to_halo%put(gal_num, halo_num)
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
     gal_num = itos(int(data_gal(i, 0)))
     call gal_to_halo%get(gal_num, hashmap_char)
     if (allocated(hashmap_char)) then
        print*, hashmap_char
     end if

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
  deallocate(associations)
  call gal_to_halo%free()

contains
  subroutine read_params ()
    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    n = iargc()

    if (n < 4) then
       print*, 'Sort the galaxies between elliptics and spirals'
       print*, ''
       print*, 'Usage (all lists should come in binary format):'
       print*, '\t -fga File to pick galaxies from'
       print*, '\t -fdm File to pick dark matter halos from'
       print*, '\t -fas File to get associations between dark matter halos and galaxies'
       !TODO: write usage
    end if

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
