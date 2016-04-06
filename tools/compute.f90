module compute
  implicit none

contains
  subroutine compute_inertia_tensor(mass, pos, I_t, I_t_diag)
    real(kind=8), intent(in), dimension(:)    :: mass
    real(kind=8), intent(in), dimension(:, :) :: pos

    real(kind=8), intent(out), dimension(3, 3) :: I_t, I_t_diag

    real(kind=8) :: mtot
    real(kind=8), dimension(3:3) :: tau
    integer :: i, j
    mtot = sum(mass)

    do i = 1, 3
       do j = 1, 3
          I_t(i, j) = sum(mass*pos(:, i)*pos(:, j)) / mtot
       end do
    end do

    ! Diag the matrix
    I_t_diag = I_t
    call sytrd(I_t_diag, tau)
  end subroutine compute_inertia_tensor
end module compute
