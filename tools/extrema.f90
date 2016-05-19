module extrema
  use extrema_mod
  use extrema_types

  real(dp), allocatable, dimension(:, :) :: peaks, eigvect
  real(dp), allocatable, dimension(:)    :: eigval
  integer, allocatable, dimension(:)     :: peak_type

  integer :: NPEAKS
  
  private :: NPEAKS

contains

  subroutine set_max_npeaks(maxnpeak)
    integer, intent(in) :: maxnpeak

    NPEAKS = maxnpeak
  end subroutine set_max_npeaks
  
  !! Find the maximum of the field, and store them in the common data of the module
  !! args:
  !!   - real(dp) field(N1, N2, N3), the field
  !! computes:
  !!   - real(dp) peaks(3, Npeak) the locations of the peaks
  !!   - real(dp) eigvect(3, Npeak) the eigvector of the hessian of the peak
  !!   - real(dp) eigval(Npeak) the eigen value of the hessian of the peak
  !!   - integer peak_type(Npeak) the type of the peak
  subroutine extrema_get_peaks(field)
    real(kind=8), intent(in)  :: field(:, :, :)

    type(CND_CNTRL_TYPE) :: ctrl
    type(EXT_DATA)       :: ext(NPEAKS)
    real(kind=4) :: flattened_field(size(field))

    integer :: npeak
    ctrl%NPROC = 1
    ctrl%justprint = .false.

    ! flatten the field
    flattened_field = reshape(field, (/size(field, 1) * size(field, 2) * size(field, 3)/))

    ! get the extrema
    call find_extrema(flattened_field, nn=(/NPEAKS, NPEAKS, NPEAKS/), ext=ext, nd=3, cnd_cntrl=ctrl)

    ! store all that in the module vars
    if (allocated(peaks)) deallocate(peaks)
    if (allocated(eigvect)) deallocate(eigvect)
    if (allocated(eigval)) deallocate(eigval)
    if (allocated(peak_type)) deallocate(peak_type)

    npeak = 0
    do i = 1, NPEAKS
       if (ext(i)%typ >= 0) then
          npeak = npeak + 1
       end if
    end do

    allocate(peaks(3, npeak), eigvect(3, npeak), eigval(npeak), peak_type(npeak))
    do i = 1, NPEAKS
       if (ext(i)%typ >= 0) then
          peaks(:, i)   = ext(i)%pos
          eigvect(:, i) = ext(i)%eig
          eigval(i)     = ext(i)%val
          peak_type(i)  = ext(i)%typ
       end if
    end do

  end subroutine extrema_get_peaks

end module extrema
