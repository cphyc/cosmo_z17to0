module compute
  implicit none

  private

  integer :: opt_lwork = -1

  public :: compute_inertia_tensor

contains
  subroutine compute_inertia_tensor(mass, pos, I_t, I_t_diag, E)
    real(kind=8), intent(in), dimension(:)    :: mass
    real(kind=8), intent(in), dimension(:, :) :: pos

    real(kind=8), intent(out), dimension(3, 3) :: I_t
    real(kind=8), intent(out), dimension(3) :: I_t_diag
    real(kind=8), intent(out), dimension(2), optional :: E

    real(kind=8) :: mtot
    real(kind=8), dimension(3) :: tau
    integer :: i, j, info

    real(kind=8), dimension(max(1, opt_lwork)) :: work

    mtot = sum(mass)

    do i = 1, 3
       do j = 1, 3
          I_t(i, j) = sum(mass*pos(i, :)*pos(j, :)) / mtot
       end do
    end do

    call dsytrd('U', 3, I_t, 3, I_t_diag, E, tau, work, opt_lwork, info)
    if (info == 0) then
       opt_lwork = int(work(1))
    end if
  end subroutine compute_inertia_tensor
end module compute
