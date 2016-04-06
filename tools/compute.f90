module compute
  implicit none

  private

  integer :: opt_lwork = -1

  public :: compute_inertia_tensor

contains
  ! This routine shall be used with positions between 0 and 1. If an halo is close
  ! to the border of the box, it may have positions ~0 and others ~1. This problem
  ! is solved by the routine
  subroutine compute_inertia_tensor(mass, pos, I_t)
    real(kind=8), intent(in), dimension(:)    :: mass
    real(kind=8), intent(in), dimension(:, :) :: pos

    real(kind=8), intent(out), dimension(3, 3) :: I_t

    real(kind=8), dimension(3) :: I_t_diag
    real(kind=8), dimension(2) :: E

    real(kind=8) :: mtot, tmp
    real(kind=8), dimension(3) :: tau
    real(kind=8), dimension(:, :), allocatable :: corr_pos
    integer :: i, j, info

    real(kind=8), dimension(max(1, opt_lwork)) :: work

    allocate(corr_pos(3, ubound(pos, 2)))
    do i = 1, size(pos, 1)
       call clean_data(pos(i, :), corr_pos(i, :))
    end do

    mtot = sum(mass)

    do i = 1, 3
       do j = i, 3
          tmp = sum(mass*corr_pos(i, :)*corr_pos(j, :)) / mtot
          I_t(i, j) = tmp
          I_t(j, i) = tmp
       end do
    end do

    ! call dsytrd('U', 3, I_t, 3, I_t_diag, E, tau, work, opt_lwork, info)
    ! if (info == 0) then
    !    opt_lwork = int(work(1))
    ! end if
  contains
    subroutine clean_data(data, c_data)
      real(kind=8), intent(in), dimension(:)           :: data
      real(kind=8), intent(out), dimension(size(data)) :: c_data

      c_data = data

    end subroutine clean_data
  end subroutine compute_inertia_tensor
end module compute
