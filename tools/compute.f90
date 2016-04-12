module compute
  implicit none

  private

  integer :: opt_lwork = -1

  public :: compute_inertia_tensor, correct_positions

contains
  ! This routine shall be used with positions between 0 and 1. If an halo is close
  ! to the border of the box, it may have positions ~0 and others ~1. This problem
  ! is solved by the routine and the position are changed
  subroutine compute_inertia_tensor(mass, pos, I_t)
    real(kind=8), intent(in), dimension(:)       :: mass
    real(kind=8), intent(inout), dimension(:, :) :: pos

    real(kind=8), intent(out), dimension(3, 3) :: I_t

    real(kind=8), dimension(3) :: I_t_diag, tau
    real(kind=8), dimension(2) :: E
    real(kind=8), dimension(max(1, opt_lwork)) :: work
    real(kind=8) :: mtot, tmp

    integer :: i, j, k, info
    integer :: ndim, nparts, dim, part

    ndim = size(pos, 1)
    nparts = size(pos, 2)

    mtot = sum(mass)

    I_t = 0

    !-------------------------------------
    ! Compute the tensor
    !-------------------------------------
    do i = 1, ndim
       do j = i, ndim
          tmp = sum(mass(:)*pos(i, :)*pos(j, :)) / mtot
          I_t(i, j) = tmp + I_t(i, j)
          I_t(j, i) = tmp + I_t(i, j)
       end do
    end do


    ! call dsytrd('U', 3, I_t, 3, I_t_diag, E, tau, worpart, opt_lwork, info)
    ! if (info == 0) then
    !    opt_lwork = int(work(1))
    ! end if

  end subroutine compute_inertia_tensor

  subroutine correct_positions (pos)
    real(kind=8), intent(inout), dimension(:, :) :: pos

    real(kind=8), dimension(size(pos, 1)) :: span
    integer :: ndim, nparts, part, dim

    ndim = size(pos, 1)
    nparts = size(pos, 2)

    span = maxval(pos, 1) - minval(pos, 1)

    !-------------------------------------
    ! If the particles have a separation > 0.5
    ! then shift 'em
    !-------------------------------------
    do dim = 1, 3
       if (span(dim) > 0.5) then
          print*, 'Correcting on dimension', dim, 'span:', span
          do part = 1, nparts
             if (pos(dim, part) > 0.5) then
                pos(dim, part) = pos(dim, part) - 1d0
             end if
          end do
       end if
    end do
  end subroutine correct_positions

end module compute
