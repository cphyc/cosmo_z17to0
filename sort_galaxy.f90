program sort_galaxy
  implicit none
  character(LEN=128) :: gal_list_filename, repository, outfile, fileinfo='none', halogalfile
  character(LEN=128) :: dm_halo_list_filename
  logical :: verbose
  integer :: ngal, gal_cols, i
  integer :: ndmhalo, dm_halo_cols
  integer :: ellipticals = 0, spirals = 0, others = 0
  real(kind=4), dimension(:,:), allocatable :: data
  real(kind=4), dimension(:, :), allocatable :: associations ! contains association flag / dark matter id / mass / gal id / mass
  real(kind=4) :: dmhalo_id, lvl, dmhalo_mass, gal_id, gal_mass
  real(kind=4), dimension(:), allocatable :: kind ! 0 for elliptical, 1 for spirals, 2 for other
  real(kind=4) :: sigma
  real(kind=4) :: elliptical_threshold = 1.5d0 ! FIXME
  real(kind=4) :: spiral_threshold = 0.8d0 ! FIXME

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

  print*, 'Reading file', gal_list_filename
  open(unit=11, file=gal_list_filename, form='unformatted')
  read(11) ngal, gal_cols
  print*, 'Got', ngal, 'galaxies in da pocket'

  print*, 'Reading file', dm_halo_list_filename
  open(unit=12, file=dm_halo_list_filename, form='unformatted')
  read(12) ndmhalo, dm_halo_cols
  !-------------------------------------
  ! Allocations & datareading
  !-------------------------------------
  allocate(data(ngal, gal_cols+1))
  allocate(associations(ndmhalo, dm_halo_cols+1))

  allocate(kind(ngal))

  read(11) data(1:ngal, 1:gal_cols)
  read(12) associations(1:ndmhalo, 1:dm_halo_cols)

  !-------------------------------------
  ! Processing
  !-------------------------------------
  do i = 1, ngal
     sigma = 1d0/3d0 * sqrt(data(i, 3)**2 + data(i, 4)**2 + data(i, 5)**2)
     if (sigma / data(i, 2) > elliptical_threshold) then
        data(i, gal_cols+1) = 0
        ellipticals = ellipticals + 1
     else if (sigma / data(i, 2) < spiral_threshold) then
        data(i, gal_cols+1) = 1
        spirals = spirals + 1
     else
        kind(i) = 2
     end if

  end do

  ! Sort the galaxies
  print*, ellipticals*100./ngal, "% ellipticals"
  print*, spirals*100./ngal, "% spirals"
  print*, 100 - (ellipticals + spirals)*100./ngal, "% others"

  !-------------------------------------
  ! Load halo-galaxy association file halo
  !-------------------------------------
  open(unit=12, file=halogalfile)
  read(12) dmhalo_id, lvl, dmhalo_mass, gal_id, gal_mass
  associations(dmhalo_id, 1) = .true.
  associations(dmhalo_id, 2) = lvl
  associations(dmhalo_id, 3) = dmhalo_mass
  associations(dmhalo_id, 4) = gal_id
  associations(dmhalo_id, 5) = gal_mass

  !-------------------------------------
  ! Cleanup
  !-------------------------------------
  close(11)
  close(12)
  deallocate(data)

contains
  subroutine read_params ()
    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    n = iargc()

    if (n < 4) then
       print*, 'Sort the galaxies between elliptics and spirals'
       !TODO: write usage
    end if

    do i = 1, n, 2
       call getarg(i, opt)
       if (i == n) then
          print*,  '("option ",a2," has no argument")', opt
          stop
       end if

       call getarg(i+1, arg)
       select case(opt)
       case ('-gal_list')
          gal_list_filename = trim(arg)
       case ('-dm_halo_list')
          dm_halo_list_filename = trim(arg)
       case ('-out')
          outfile = trim(arg)
       case ('-fin')
          fileinfo = trim(arg)
       case ('-halogal')
          halogalfile = trim(arg)
       case ('-ver')
          read (arg,*) verbose
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

  end subroutine read_params

end program sort_galaxy
