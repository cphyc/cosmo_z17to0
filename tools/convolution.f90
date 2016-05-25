module convolution
  use, intrinsic :: iso_c_binding
  use types
  implicit none

  include 'fftw3.f03'
  ! use fftw3

  type :: CONV_T
     real(kind=8), allocatable    :: A(:,:,:), B(:,:,:), conv(:,:,:)
     complex(kind=8), allocatable :: dftA(:,:,:), dftB(:,:,:), dftConv(:,:,:)
   contains
     procedure :: init_A => conv_init_A
     procedure :: init_B => conv_init_B
     procedure :: execute => conv_execute
     procedure :: free => conv_free
  end type CONV_T

  real(kind=8), parameter :: pi = 3.14159265358979d0
  private :: pi

contains
  !> Generate a 3d gaussian kernel
  !! parameters:
  !! s: the number of elements in the kernel
  !! sigma: the sigma of the kernel
  !! kernel: a 3d array containing the kernel
  subroutine kernel_gaussian3d(s, sigma, kernel)
    integer, intent(in) :: s
    real(kind=8), intent(in) :: sigma
    real(kind=8), dimension(s, s, s), intent(out) :: kernel

    integer :: i, j, k
    real(kind=8) :: sizeo2, sq_sum

    sizeo2 = s/2.

    do k = 1, s
       do j = 1, s
          do i = 1, s
             sq_sum = (i - sizeo2)**2 + (j - sizeo2)**2 + (k - sizeo2)**2
             kernel(i,j,k) = 1/(sqrt(2*pi) * sigma)**3 * exp(-sq_sum / (2*sigma**2))
          end do
       end do
    end do

  end subroutine kernel_gaussian3d

  !> Compute the fft of in, putting the result in out
  subroutine fft(in, out)
    real(kind=8), intent(inout), dimension(:,:,:) :: in
    complex(kind=8), intent(out), dimension(size(in, 1), size(in, 2),&
         & size(in, 3)) :: out
    complex(kind=8), dimension(size(in, 1), size(in, 2),&
         & size(in, 3)) :: cplx_in

    type(C_PTR) :: plan

    integer :: L, N, M

    L = size(in, 1)
    M = size(in, 2)
    N = size(in, 3)

    cplx_in = cmplx(in)

    ! create plan
    plan = fftw_plan_dft_3d(N, M, L, cplx_in, out, FFTW_FORWARD, FFTW_ESTIMATE)

    ! execute it
    call fftw_execute_dft(plan, cplx_in, out)

    ! delete plan
    call fftw_destroy_plan(plan)

  end subroutine fft

  subroutine ifft(in, out)
    complex(kind=8), intent(inout), dimension(:,:,:) :: in

    real(kind=8), intent(out), dimension(size(in, 1), size(in, 2), size(in, 3)) :: out

    complex(kind=8), dimension(size(in, 1), size(in, 2), size(in, 3)) :: cplx_out

    type(C_PTR) :: plan

    integer :: L, N, M

    L = size(in, 1)
    M = size(in, 2)
    N = size(in, 3)

    ! create plan
    plan = fftw_plan_dft_3d(N, M, L, in, cplx_out, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! execute it
    call fftw_execute_dft(plan, in, cplx_out)

    out = dble(cplx_out) / (L*M*N)

    ! delete plan
    call fftw_destroy_plan(plan)

  end subroutine ifft

  !> Prepare the convolution A*B by giving A
  subroutine conv_init_A (self, A)
    real(kind=8), intent(in)  :: A(:,:,:)
    class(CONV_T), intent(inout) :: self

    allocate(self%A(size(A, 1), size(A, 2), size(A, 3)))
    allocate(self%dftA(size(A, 1), size(A, 2), size(A, 3)))
    self%A = A

    call fft(self%A, self%dftA)
  end subroutine conv_init_A

  !> Prepare the convolution A*B by giving B
  subroutine conv_init_B (self, B)
    real(kind=8), intent(in)      :: B(:,:,:)
    class(CONV_T), intent(inout) :: self

    allocate(self%B(size(B, 1), size(B, 2), size(B, 3)))
    allocate(self%dftB(size(B, 1), size(B, 2), size(B, 3)))
    self%B = B

    call fft(self%B, self%dftB)
  end subroutine conv_init_B

  !> Execute the convolution
  subroutine conv_execute(self)
    class(CONV_T), intent(inout) :: self

    integer :: L, M, N
    L = size(self%dftA, 1)
    M = size(self%dftA, 2)
    N = size(self%dftA, 3)

    allocate(self%dftConv(L, M, N))
    allocate(self%conv(L, M, N))

    call conv_prod(self%dftA, self%dftB, self%dftConv)
    call ifft(self%dftConv, self%conv)
  end subroutine conv_execute

  !> Compute the convolution product of A, B in fourier space
  subroutine conv_prod(A, B, C)
    complex(kind=8), intent(in)  :: A(:,:,:), B(:,:,:)
    complex(kind=8), intent(out) :: C(size(A, 1), size(A, 2), size(A, 3))

    integer :: i,j,k
    do k = 1, size(A, 3)
       do j = 1, size(A, 2)
          do i = 1, size(A, 1)
             ! from
             ! http://www.fftw.org/faq/section3.html#centerorigin
             ! with special care with the fact that fortran starts
             ! indexing at 1, not 0!
             C(i,j,k) = A(i,j,k) * B(i,j,k) * (-1)**(i+j+k+1)
          end do
       end do
    end do

  end subroutine conv_prod

  !! Compute the 3d histogram of data, using nbin and weights and
  !! store it into bins
  subroutine conv_hist3d(data, nbin, weights, hist, edges)
    real(kind=8), dimension(:,:), intent(in) :: data !! data(ndim, nparts)
    real(kind=8), dimension(size(data, 2)), &
         intent(in), optional            :: weights ! weights(nparts)
    integer, intent(in)                  :: nbin

    real(kind=8), dimension(nbin, nbin, nbin), intent(out)    :: hist
    real(kind=8), dimension(3, nbin+1), intent(out), optional :: edges

    real(kind=8), dimension(3) :: maxis, minis, spans

    integer :: i, j, k, part_i, ndim, nparts

    ! get the dimensions
    ndim = size(data, 1)
    nparts = size(data, 2)
    ! compute the edges
    maxis = maxval(data, 2)
    minis = minval(data, 2)
    spans = maxis - minis

    do i = 1, nbin + 1
       edges(:, i) = minis + (maxis - minis) * (i-1) / nbin
    end do
    print*, minis, maxis

    print*, nparts
    ! project the data onto the edges
    do part_i = 1, nparts
       ! get the position in the grid
       i = floor((data(1, part_i) - minis(1)) * nbin / spans(1)) + 1
       j = floor((data(2, part_i) - minis(2)) * (nbin * 1d0) / spans(2)) + 1
       k = floor((data(3, part_i) - minis(3)) * nbin / spans(3)) + 1

       ! because of rounding errors, the maxima aren't well found,
       ! fix that
       if (data(1, part_i) == maxis(1)) i = nbin
       if (data(2, part_i) == maxis(2)) j = nbin
       if (data(3, part_i) == maxis(3)) k = nbin

       ! print*, i,j,k
       ! TODO: linear approximation of density
       ! increment the histogram
       if ((i <= 0) .or. (i > nbin) .or. (j <= 0) .or. (j > nbin) .or.&
            (k <= 0) .or. (k > nbin)) then
          print*, 'Dafuk ?!'
          print*, (data(:, part_i) - maxis) == 0
          print*, i, j, k
       end if

       if (present(weights)) then
          hist(i, j, k) = hist(i, j, k) + weights(part_i)
       else
          hist(i, j, k) = hist(i, j, k) + 1
       end if
    end do

    print*, nparts, ndim
  end subroutine conv_hist3d

  !! Estimate the density
  subroutine conv_density(data, nbin, dens, edges)
    real(kind=8), dimension(:,:), intent(in) :: data !! data(ndim, nparts)
    integer, intent(in)                      :: nbin

    real(kind=8), dimension(nbin, nbin, nbin), intent(out)   :: dens
    real(kind=8), dimension(3, nbin+1), intent(out), optional:: edges

    real(kind=8), dimension(3) :: maxis, minis, spans

    integer :: i, j, k, part_i, ndim, nparts, i0, j0, k0, imax, jmax, kmax, imin, jmin, kmin, dim
    real(kind=8) :: ri, rj, rk

    ! get the dimensions
    ndim = size(data, 1)
    nparts = size(data, 2)
    ! compute the bins
    maxis = maxval(data, 2)
    minis = minval(data, 2)

    spans = (maxis - minis)

    do i = 1, nbin + 1
       edges(:, i) = minis + spans * (i-1) / nbin
    end do

    ! project the data onto the bins
    do part_i = 1, nparts
       ! get the position in the grid
       ri = (data(1, part_i) - minis(1)) * nbin / spans(1) + 0.5
       rj = (data(2, part_i) - minis(2)) * nbin / spans(2) + 0.5
       rk = (data(3, part_i) - minis(3)) * nbin / spans(3) + 0.5

       ! the max and min are actually only half in leftmost/rightmost bin
       ! so the min has ri =  0.5
       ! and the max has ri = nbin + 0.5

       i0 = floor(ri)
       j0 = floor(rj)
       k0 = floor(rk)

       imin = max(1, i0)
       jmin = max(1, j0)
       kmin = max(1, k0)

       imax = min(nbin, i0 + 1)
       jmax = min(nbin, j0 + 1)
       kmax = min(nbin, k0 + 1)

       do k = kmin, kmax
          do j = jmin, jmax
             do i = imin, imax
                ! add the proportion of the particle in the box
                dens(i, j, k) = dens(i, j, k) + abs((ri - i)*(rj - j)*(rk - k))
             end do
          end do
       end do
    end do

    ! Convert the number density into a mass density (assuming uniform mass)
    ! TODO: don't assume uniform mass
    dens = dens * nbin**3 / product(maxis-minis)

  end subroutine conv_density

  subroutine conv_free(self)
    class(CONV_T), intent(inout) :: self
    if (allocated(self%A)) deallocate(self%A)
    if (allocated(self%B)) deallocate(self%B)
    if (allocated(self%conv)) deallocate(self%conv)

    if (allocated(self%dftA)) deallocate(self%dftA)
    if (allocated(self%dftB)) deallocate(self%dftB)
    if (allocated(self%dftConv)) deallocate(self%dftConv)

  end subroutine conv_free

end module convolution
