module convolution
  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'
  ! use fftw3

  type :: CONV_T(k)
     integer, kind                :: k
     real(kind=k), allocatable    :: A(:,:,:), B(:,:,:), conv(:,:,:)
     complex(kind=k), allocatable :: dftA(:,:,:), dftB(:,:,:), dftConv(:,:,:)
   contains
     procedure :: init_A => conv_init_A
     procedure :: init_B => conv_init_B
     procedure :: execute => conv_execute
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

  !> Compute the fft of int, putting the result in out
  subroutine fft(in, out)
    real(C_DOUBLE), intent(inout), dimension(:,:,:)                       :: in
    complex(C_DOUBLE_COMPLEX), intent(out), dimension(:,:,:), allocatable :: out

    type(C_PTR) :: plan, data

    integer :: L, N, M, flags

    L = size(in, 1)
    M = size(in, 2)
    N = size(in, 3)

    ! from http://www.fftw.org/doc/Multi_002dDimensional-DFTs-of-real-data.html#Multi_002dDimensional-DFTs-of-real-data

    allocate(out(L/2+1, M, N))

    ! create plan
    plan = fftw_plan_dft_r2c_3d(N, M, L, in, out, FFTW_ESTIMATE)

    ! execute it
    call fftw_execute_dft_r2c(plan, in, out)

    ! delete plan
    call fftw_destroy_plan(plan)

  end subroutine fft

  subroutine ifft(in, out)
    complex(C_DOUBLE_COMPLEX), intent(inout), dimension(:,:,:) :: in
    real(C_DOUBLE), intent(out), dimension(:,:,:), allocatable :: out

    type(C_PTR) :: plan, data

    integer :: L, N, M, flags

    L = size(in, 1)
    M = size(in, 2)
    N = size(in, 3)

    ! from http://www.fftw.org/doc/Multi_002dDimensional-DFTs-of-real-data.html#Multi_002dDimensional-DFTs-of-real-data

    allocate(out(L*2, M, N))

    ! create plan
    plan = fftw_plan_dft_c2r_3d(N, M, L, in, out, FFTW_ESTIMATE)

    ! execute it
    call fftw_execute_dft_c2r(plan, in, out)

    ! delete plan
    call fftw_destroy_plan(plan)

  end subroutine ifft

  !> Prepare the convolution A*B by giving A
  subroutine conv_init_A (conv_object, A)
    real(kind=8), intent(in)      :: A(:,:,:)
    class(CONV_T(8)), intent(out) :: conv_object

    allocate(conv_object%A(size(A, 1), size(A, 2), size(A, 3)))
    conv_object%A = A

    call fft(conv_object%A, conv_object%dftA)
  end subroutine conv_init_A

  !> Prepare the convolution A*B by giving B
  subroutine conv_init_B (conv_object, B)
    real(kind=8), intent(in)      :: B(:,:,:)
    class(CONV_T(8)), intent(inout) :: conv_object

    allocate(conv_object%B(size(B, 1), size(B, 2), size(B, 3)))
    conv_object%B = B

    call fft(conv_object%B, conv_object%dftB)
  end subroutine conv_init_B

  !> Execute the convolution
  subroutine conv_execute(conv_object)
    class(CONV_T(8)), intent(inout) :: conv_object

    integer :: L, M, N
    L = size(conv_object%dftA, 1)
    M = size(conv_object%dftA, 2)
    N = size(conv_object%dftA, 3)

    allocate(conv_object%dftConv(L, M, N))

    conv_object%dftConv = conv_object%dftA * conv_object%dftB
    call ifft(conv_object%dftConv, conv_object%conv)
  end subroutine conv_execute

  subroutine conv_hist3d(data, nbin, hist, bins)
    real(kind=8), dimension(:,:), intent(in) :: data ! data(ndim, nparts)
    integer, intent(in)                      :: nbin

    real(kind=8), dimension(nbin, nbin, nbin), intent(out) :: hist
    real(kind=8), dimension(3, nbin), intent(out), optional :: bins

    real(kind=8), dimension(3) :: maxis, minis, spans

    integer :: i, j, k, part_i, ndim, nparts

    ! get the dimensions
    ndim = size(data, 1)
    nparts = size(data, 2)
    ! compute the bins
    maxis = maxval(data, 2)
    minis = minval(data, 2)
    spans = maxis - minis

    do i = 1, nbin
       bins(:, i) = minis + (maxis - minis) * (i-1) / nbin
    end do
    print*, minis, maxis

    print*, nparts
    hist = 0
    ! project the data onto the bins
    do part_i = 1, nparts
       ! get the position in the grid
       i = floor((data(1, part_i) - minis(1)) * nbin / spans(1)) + 1
       j = floor((data(2, part_i) - minis(2)) * (nbin * 1d0) / spans(2)) + 1
       k = floor((data(3, part_i) - minis(3)) * nbin / spans(3)) + 1

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
       hist(i, j, k) = hist(i, j, k) + 1
    end do

    print*, nparts, ndim
  end subroutine conv_hist3d
end module convolution
