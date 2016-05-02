module readTree
  implicit none
  public
  integer :: nsteps

  integer, dimension(:), allocatable       :: nb_of_halos, nb_of_subhalos !> Number of halos and subhalos in the tree
  real (kind=4), dimension(:), allocatable :: aexp, age_univ, omega_t     !> Parameters of each outputs

  ! Internal variables
  integer :: istep = 1
  integer :: j = 1
  integer :: unit

  ! Some temp variable
  integer, allocatable       :: tmp_id_father(:), tmp_id_son(:)
  real(kind=4), allocatable  :: tmp_m_father(:)

contains
  !> Read a merger tree, storing the data about it in the module.
  subroutine read_tree(filename)

    character(len=*), intent(in) :: filename !> Tree filename

    open(newunit=unit, file=trim(filename), status='old', form='unformatted')
    read(unit) nsteps

    ! Deallocate all previously allocated data
    call free()

    allocate(nb_of_halos(nsteps), nb_of_subhalos(nsteps), aexp(nsteps), &
         omega_t(nsteps), age_univ(nsteps))
    read(unit) nb_of_halos, nb_of_subhalos
    read(unit) aexp
    read(unit) omega_t
    read(unit) age_univ

  end subroutine read_tree

  !> Read the next part of the tree. If the end is reached, set stop_now to True.
  !! To get the father / son of the halo, call directly after 'get_father' or 'get_son'
  subroutine iter_tree(stop_now, halo_id, BushID, step, level, hosthalo, hostsub, &
       nbsub, nextsub, nb_of_fathers, nb_of_son, &
       m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz, &
       r, ra, rb, rc, ek, ep, et, spin)

    logical, intent(out), optional         :: stop_now      !> true if end is reached
    integer(kind=4), intent(out), optional :: halo_id       !> id if the halo
    integer(kind=4), intent(out), optional :: BushID        !> ???
    integer(kind=4), intent(out), optional :: step          !> current step
    integer(kind=4), intent(out), optional :: level         !> level of the halo (1 is halo, 2+ are subhalos)
    integer(kind=4), intent(out), optional :: hosthalo      !> id of hosting halo
    integer(kind=4), intent(out), optional :: hostsub       !> ???
    integer(kind=4), intent(out), optional :: nbsub         !> number of subhalos
    integer(kind=4), intent(out), optional :: nextsub       !> ???
    integer(kind=4), intent(out), optional :: nb_of_fathers
    integer(kind=4), intent(out), optional :: nb_of_son
    real(kind=4), intent(out), optional    :: m             !> Mass of the halo
    real(kind=4), intent(out), optional    :: macc          !> Accreted mass
    real(kind=4), intent(out), optional    :: px, py, pz    !> position of the halo
    real(kind=4), intent(out), optional    :: vx, vy, vz    !> velocitiy of the halo
    real(kind=4), intent(out), optional    :: Lx, Ly, Lz    !> angular momentum of the halo
    real(kind=4), intent(out), optional    :: r, ra, rb, rc !> ellispoid description
    real(kind=4), intent(out), optional    :: ek, ep, et    !> Kinetic, potential and thermic energy
    real(kind=4), intent(out), optional    :: spin          !> spin of the halo

    integer :: ntot

    stop_now = .false.
    ! Loop
    j = j + 1
    ntot = nb_of_subhalos(istep) + nb_of_halos(istep)

    if (j > ntot) then
       j = 1
       istep = istep + 1
       if (istep > nsteps) then
          close(unit)
          stop_now = .true.
          return
       end if
    end if

    read(unit) halo_id
    read(unit) BushID
    read(unit) step
    read(unit) level, hosthalo, hostsub, nbsub, nextsub
    read(unit) m
    read(unit) macc
    read(unit) px, py, pz
    read(unit) vx, vy, vz
    read(unit) Lx, Ly, Lz
    read(unit) r, ra, rb, rc
    read(unit) ek, ep, et
    read(unit) spin

    read(unit) nb_of_fathers

    if (nb_of_fathers > 0) then
       if (allocated(tmp_id_father)) deallocate(tmp_id_father)
       if (allocated(tmp_m_father )) deallocate(tmp_m_father)

       allocate(tmp_id_father(nb_of_fathers), tmp_m_father(nb_of_fathers))

       read(unit) tmp_id_father
       read(unit) tmp_m_father
    endif

    read(unit) nb_of_son
    if(nb_of_son>0)then
       if (allocated(tmp_id_son)) then
          deallocate(tmp_id_son)
       end if
       allocate(tmp_id_son(nb_of_son))

       read(unit) tmp_id_son
    endif

    read(unit)
    read(unit)
    ! if(star) read(unit)
  end subroutine iter_tree

  !> Get the father of the last red halo
  subroutine get_father(nb_of_fathers, id_father, m_father)
    integer, intent(in)                                 :: nb_of_fathers !> number of son of the halo
    integer, dimension(nb_of_fathers), intent(out)      :: id_father     !> id of the father
    real(kind=4), dimension(nb_of_fathers), intent(out) :: m_father      !> mass of the father

    id_father = tmp_id_father
    m_father = tmp_m_father
  end subroutine get_father

  !> Get the son of the last red halo
  subroutine get_son(nb_of_son, id_son)
    integer, intent(in)                        :: nb_of_son !> number of son of the halo
    integer, dimension(nb_of_son), intent(out) :: id_son    !> ids of the son of the halo

    id_son = tmp_id_son
  end subroutine get_son

  !> Close the tree file
  subroutine close_file()
    ! TODO: inquire it is opened
    close(unit)
  end subroutine close_file

  !> Free allocated variables
  subroutine free()
    if (allocated(nb_of_halos    )) deallocate(nb_of_halos)
    if (allocated(nb_of_subhalos )) deallocate(nb_of_subhalos)
    if (allocated(aexp           )) deallocate(aexp)
    if (allocated(omega_t        )) deallocate(omega_t)
    if (allocated(age_univ       )) deallocate(age_univ)
  end subroutine free

end module readTree
