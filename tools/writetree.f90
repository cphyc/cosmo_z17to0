module writetree
  implicit none
  public

  integer :: unit

  private :: unit
contains
  !> Read a merger tree, storing the data about it in the module.
  subroutine write_tree(filename, nsteps, nb_of_halos, nb_of_subhalos, aexp, omega_t, age_univ)

    character(len=*), intent(in) :: filename !> Tree filename

    integer, intent(in) :: nsteps !> Number of steps in the simulation

    integer, dimension(:), intent(in)       :: nb_of_halos    !> halos in the tree
    integer, dimension(:), intent(in)       :: nb_of_subhalos !> subhalos in the tree
    real (kind=4), dimension(:), intent(in) :: aexp     !> expansion factor at each step
    real (kind=4), dimension(:), intent(in) :: age_univ !> age of the universe at each step
    real (kind=4), dimension(:), intent(in) :: omega_t  !> omega_t at at each step

    open(newunit=unit, file=trim(filename), form='unformatted')
    write(unit) nsteps

    write(unit) nb_of_halos, nb_of_subhalos
    write(unit) aexp
    write(unit) omega_t
    write(unit) age_univ

  end subroutine write_tree

  !> Read the next part of the tree. If the end is reached, set stop_now to True.
  !! To get the father / son of the halo, call directly after 'get_father' or 'get_son'
  subroutine write_iter_tree(stop_now, halo_id, BushID, step, level, hosthalo, hostsub, &
       nbsub, nextsub, nb_of_fathers, nb_of_sons, &
       m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz, &
       r, ra, rb, rc, ek, ep, et, spin, &
       id_fathers, m_fathers, id_sons)

    logical, intent(in)         :: stop_now
    integer(kind=4), intent(in) :: halo_id       !> id if the halo
    integer(kind=4), intent(in) :: BushID        !> ???
    integer(kind=4), intent(in) :: step          !> current step
    integer(kind=4), intent(in) :: level         !> level of the halo (1 is halo, 2+ are subhalos)
    integer(kind=4), intent(in) :: hosthalo      !> id of hosting halo
    integer(kind=4), intent(in) :: hostsub       !> ???
    integer(kind=4), intent(in) :: nbsub         !> number of subhalos
    integer(kind=4), intent(in) :: nextsub       !> ???
    integer(kind=4), intent(in) :: nb_of_fathers
    integer(kind=4), intent(in) :: nb_of_sons
    real(kind=4), intent(in)    :: m             !> Mass of the halo
    real(kind=4), intent(in)    :: macc          !> Accreted mass
    real(kind=4), intent(in)    :: px, py, pz    !> position of the halo
    real(kind=4), intent(in)    :: vx, vy, vz    !> velocitiy of the halo
    real(kind=4), intent(in)    :: Lx, Ly, Lz    !> angular momentum of the halo
    real(kind=4), intent(in)    :: r, ra, rb, rc !> ellispoid description
    real(kind=4), intent(in)    :: ek, ep, et    !> Kinetic, potential and thermic energy
    real(kind=4), intent(in)    :: spin          !> spin of the halo
    integer(kind=4), intent(in), dimension(:), optional :: id_fathers, id_sons
    real(kind=4), intent(in), dimension(:), optional    :: m_fathers

    integer :: ntot

    if (stop_now) then
       write(unit) -1
       print*, 'Stopping because halo_id == 1'
       return
    end if
    
    write(unit) halo_id
    write(unit) BushID
    write(unit) step
    write(unit) level, hosthalo, hostsub, nbsub, nextsub
    write(unit) m
    write(unit) macc
    write(unit) px, py, pz
    write(unit) vx, vy, vz
    write(unit) Lx, Ly, Lz
    write(unit) r, ra, rb, rc
    write(unit) ek, ep, et
    write(unit) spin

    write(unit) nb_of_fathers

    if (nb_of_fathers > 0) then
       write(unit) id_fathers
       write(unit) m_fathers
    endif

    write(unit) nb_of_sons
    if(nb_of_sons>0)then
       write(unit) id_sons
    endif

    write(unit)
    write(unit)
    ! if(star) write(unit)
  end subroutine write_iter_tree

  !> Close the tree file
  subroutine close_file()
    logical :: op
    ! TODO: inquire it is opened
    inquire(unit, opened=op)
    if (op) close(unit)
  end subroutine close_file

end module writetree
