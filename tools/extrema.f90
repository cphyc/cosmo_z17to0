module extrema
  use convolution
  use types
  use extrema_mod
  use extrema_types, only : EXT_DATA, CND_CNTRL_TYPE
  use extrema_storage

  real(kind=8), allocatable, dimension(:, :) :: mod_peaks, mod_eigvect
  real(kind=8), allocatable, dimension(:)    :: mod_eigval
  integer, allocatable, dimension(:)         :: mod_peak_type, mod_index

  integer :: MOD_NPEAKS, MOD_NBINS(3), MOD_NPROC = 1
  private :: MOD_NPEAKS, MOD_NBINS, MOD_NPROC
contains
  !! Set the values
  subroutine extrema_set(npeaks, nbins, nprocs)
    integer, intent(in), optional :: npeaks, nbins, nprocs

    if (present(npeaks)) MOD_NPEAKS=npeaks
    if (present(nbins))  MOD_NBINS=NBINS
    if (present(nrpcos)) MOD_NPROCS=nprocs
  end subroutine extrema_set

  !! Call find_extrema from extrema/ files
  !! args:
  !!   - real(dp) field(N1, N2, N3) the field
  !! returns:
  !!   - type(EXT_DATA) ext(NPEAK) the peaks
  !!
  !! Note: this routine can be parallelized (but extrema_compute cannot)
  subroutine extrema_compute_ext(field, ext)
    real(kind=8), intent(in)    :: field(:, :, :)
    type(EXT_DATA), intent(out) :: ext(MOD_NPEAKS)

    type(CND_CNTRL_TYPE) :: ctrl
    real(kind=4)         :: flattened_field(size(field))

    ndim = 3

    ctrl%MOD_NPROC = MOD_NPROC
    ctrl%justprint = .false.

    ! flatten the field
    flattened_field = reshape(real(field), (/size(field)/))

    ! get the extrema
    call find_extrema(flattened_field, nn=MOD_NBINS, ext=ext, nd=3, cnd_cntrl=ctrl)

  end subroutine extrema_compute_ext

  !! Find the maximum of the field, and store them in the common data of the module.
  !! it is only a f2py compatible wrapper around extrema_compute_ext
  !! args:
  !!   - real(dp) field(N1, N2, N3), the field
  !! returns:
  !!   - integer ndim, npeak
  !! computes:
  !!   - real(dp) mod_peaks(3, Npeak) the locations of the peaks
  !!   - real(dp) mod_eigvect(3, Npeak) the eigvector of the hessian of the peak
  !!   - real(dp) mod_eigval(Npeak) the eigen value of the hessian of the peak
  !!   - integer mod_peak_type(Npeak) the type of the peak
  subroutine extrema_compute(field, ndim, npeak)
    real(kind=8), intent(in) :: field(:, :, :)

    integer, intent(out) :: ndim, npeak

    type(EXT_DATA)       :: ext(MOD_NPEAKS)

    ! get the extrema
    call extrema_compute_ext(field, ext)

    ! store all that in the module vars
    if (allocated(mod_peaks))     deallocate(mod_peaks)
    if (allocated(mod_eigvect))   deallocate(mod_eigvect)
    if (allocated(mod_eigval))    deallocate(mod_eigval)
    if (allocated(mod_peak_type)) deallocate(mod_peak_type)
    if (allocated(mod_index))     deallocate(mod_index)

    !TODO: remove parallel foobars
    npeak = 0
    do i = 1, MOD_NPEAKS
       if (ext(i)%typ >= 0) then
          npeak = npeak + 1
       end if
    end do

    allocate(mod_peaks(3, npeak), mod_eigvect(3, npeak), &
         mod_eigval(npeak), mod_peak_type(npeak), mod_index(npeak))
    do i = 1, MOD_NPEAKS
       if (ext(i)%typ >= 0) then
          mod_peaks(:, i)   = ext(i)%pos
          mod_eigvect(:, i) = ext(i)%eig
          mod_eigval(i)     = ext(i)%val
          mod_peak_type(i)  = ext(i)%typ
          mod_index(i)      = ext(i)%pix
       end if
    end do

    ndim = 3

  end subroutine extrema_compute

  !! Get the peaks at some stepx
  !! args:
  !!   - integer ndim, npeak
  !! returns:
  !!   - double peakpos(ndim, npeak)    the position of the peak
  !!   - double eigvectors(ndim, npeak) the eigenvector of the hessian
  !!   - double eigvalue(npeak)         the eigenvalue of the hessian
  !!   - integer peaktype(npeak)        the type of the peak
  !!   - integer index(npeak)           the index of the peak (in contiguous array, you should reshape it)
  subroutine extrema_get_i(ndim, npeak, istep, index, peakpos, eigvectors, eigvalues, peaktype)
    integer, intent(in) :: ndim, npeak, istep
    real(kind=8), intent(out), dimension(ndim, npeak) :: peakpos
    real(kind=8), intent(out), dimension(ndim, npeak) :: eigvectors
    real(kind=8), intent(out), dimension(npeak)       :: eigvalues
    integer, intent(out), dimension(npeak)            :: peaktype, index

    type(EXT_DATA), dimension(:), allocatable :: ext

    call extrema_pull(ext, istep)

    do i = 1, npeak
       if (ext(i)%typ >= 0) then
          peakpos(:, i)    = ext(i)%pos
          eigvectors(:, i) = ext(i)%eig
          eigvalues(i)     = ext(i)%val
          peaktype(i)      = ext(i)%typ
          index(i)         = ext(i)%pix
       end if
    end do

  end subroutine extrema_get_i


  !! Get the peaks
  !! args:
  !!   - integer ndim, npeak
  !! returns:
  !!   - double peakpos(ndim, npeak)    the position of the peak
  !!   - double eigvectors(ndim, npeak) the eigenvector of the hessian
  !!   - double eigvalue(npeak)         the eigenvalue of the hessian
  !!   - integer peaktype(npeak)        the type of the peak
  !!   - integer index(npeak)           the index of the peak (in contiguous array, you should reshape it)
  subroutine extrema_get(ndim, npeak, index, peakpos, eigvectors, eigvalues, peaktype)
    integer, intent(in) :: ndim, npeak
    real(kind=8), intent(out), dimension(ndim, npeak) :: peakpos
    real(kind=8), intent(out), dimension(ndim, npeak) :: eigvectors
    real(kind=8), intent(out), dimension(npeak)       :: eigvalues
    integer, intent(out), dimension(npeak)        :: peaktype, index

    peakpos    = mod_peaks
    eigvectors = mod_eigvect
    eigvalues  = mod_eigvals
    peaktype   = mod_peak_type
    index      = mod_index

  end subroutine extrema_get

  subroutine extrema_compute_smooth_space(sigmamin, sigmamax, positions, nbin, nsigmastep, &
       sigmas)
    real(kind=8), intent(in) :: sigmamin, sigmamax
    real(kind=8), intent(in) :: positions(:, :)

    integer, intent(in) :: nbin, nsigmastep

    real(kind=8), intent(out) :: sigmas(nsigmastep)

    real(kind=8), dimension(nbin, nbin, nbin)    :: density, gaussian, smoothed_density
    real(kind=8), dimension(3, nbin+1)           :: edges
    complex(kind=8), dimension(nbin, nbin, nbin) :: fftdensity, fftgaussian, fftconv

    type(EXT_DATA), allocatable, target :: extrema(:, :)

    real(kind=8) :: sigma
    integer  :: ndim, npeak

    MOD_NPEAKS = nbin ** 3
    MOD_NBINS = (/nbin, nbin, nbin/)

    allocate(extrema(MOD_NPEAKS, nsigmastep))

    ! estimate density
    call conv_density(positions, nbin, density, edges)

    ! precompute ffts
    call fft(density, fftdensity)

    ! iterate over sigmas
    do i = 1, nsigmastep
       sigma = (sigmamax-sigmamin)*(i-1)/nsigmastep + sigmamin
       sigmas(i) = sigma

       print*, sigma

       call kernel_gaussian3d(nbin, sigma, gaussian)
       call fft(gaussian, fftgaussian)

       call conv_prod(fftdensity, fftgaussian, fftconv)

       ! get back to real space
       call ifft(fftconv, smoothed_density)

       ! find the peaks
       call extrema_compute_ext(smoothed_density, extrema(:, i))
    end do

    call extrema_push(extrema)

  end subroutine extrema_compute_smooth_space

end module extrema
