program sort_galaxy
  implicit none
  character(LEN=128) :: filename, repository, outfile, fileinfo='none', halogalfile
  logical :: verbose
  integer :: ngal, cols, i
  integer :: ellipticals = 0, spirals = 0, others = 0
  real(kind=4), dimension(:,:), allocatable :: data
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
  if (filename == 'none') then
     stop
  end if

  print*, 'Reading file', filename
  open(unit=11, file=filename, form='unformatted')
  read(11) ngal, cols
  print*, 'Got', ngal, 'galaxies in da pocket'
  !-------------------------------------
  ! Allocations & datareading
  !-------------------------------------
  allocate(data(ngal, cols+1))
  allocate(associations(ngal, 5))
  allocate(kind(ngal))

  read(11) data(1:ngal, 1:cols)

  !-------------------------------------
  ! Processing
  !-------------------------------------
  do i = 1, ngal
     sigma = 1d0/3d0 * sqrt(data(i, 3)**2 + data(i, 4)**2 + data(i, 5)**2)
     if (sigma / data(i, 2) > elliptical_threshold) then
        data(i, cols+1) = 0
        ellipticals = ellipticals + 1
     else if (sigma / data(i, 2) < spiral_threshold) then
        data(i, cols+1) = 1
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
  read(12) dmhalo_id(i), lvl(i), dmhalo_mass(i), gal_id(i), gal_mass(i)

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
       case ('-inp')
          filename = trim(arg)
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
