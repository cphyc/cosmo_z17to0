module readTree
  implicit none
  public
  integer :: nsteps

  integer, dimension(:), allocatable       :: nb_of_halos, nb_of_subhalos
  real (kind=4), dimension(:), allocatable :: aexp, age_univ, omega_t

  ! type HALO_LIST_T
  !    type(HALO_T), pointer :: halo
  ! end type HALO_LIST_T

  ! type HALO_T
  !    type(HALO_T), pointer          :: parent
  !    type(HALO_LIST_T), allocatable :: childs(:)
  !    integer                        :: mynumber, BushID, mystep, level, hosthalo, hostsub, &
  !         nbsub, nextsub, nb_of_fathers, nb_of_son
  !    integer, allocatable           :: idfather(:), idson(:)
  !    real(kind=4)                   :: m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz, &
  !         r, ra, rb, rc, ek, ep, et, spin
  !    real(kind=4), allocatable      :: m_father(:)
  !    integer, allocatable           :: id_father(:), idson(:)
  ! end type HALO_T

  integer :: istep = 1
  integer :: j = 1
  integer :: unit

  ! Some temp variable
  integer, allocatable       :: tmp_id_father(:), tmp_id_son(:)
  real(kind=4), allocatable  :: tmp_m_father(:)

contains
  subroutine read_tree(filename)

    character(len=*), intent(in) :: filename

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

  subroutine iter_tree(stop_now, mynumber, BushID, mystep, level, hosthalo, hostsub, &
       nbsub, nextsub, nb_of_fathers, nb_of_son, &
       m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz, &
       r, ra, rb, rc, ek, ep, et, spin)
    logical, intent(out)         :: stop_now
    integer(kind=4), intent(out) :: mynumber, BushID, mystep, level, hosthalo, hostsub
    integer(kind=4), intent(out) :: nbsub, nextsub, nb_of_fathers, nb_of_son
    real(kind=4), intent(out)    :: m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz
    real(kind=4), intent(out)    :: r, ra, rb, rc, ek, ep, et, spin

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

    read(unit) mynumber
    read(unit) BushID
    read(unit) mystep
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

  subroutine get_fathers(nb_of_fathers, id_father, m_father)
    integer, intent(in) :: nb_of_fathers

    integer, dimension(nb_of_fathers), intent(out)      :: id_father
    real(kind=4), dimension(nb_of_fathers), intent(out) :: m_father

    !f2py depend(nb_of_fathers) :: id_father, m_father
    id_father = tmp_id_father
    m_father = tmp_m_father
  end subroutine get_fathers

  subroutine get_son(nb_of_son, id_son)
    integer, intent(in) :: nb_of_son

    integer, dimension(nb_of_son), intent(out) :: id_son
    id_son = tmp_id_son
  end subroutine get_son

  subroutine close_file()
    ! TODO: inquire it is opened
    close(unit)
  end subroutine close_file

  subroutine free()
    if (allocated(nb_of_halos    )) deallocate(nb_of_halos)
    if (allocated(nb_of_subhalos )) deallocate(nb_of_subhalos)
    if (allocated(aexp           )) deallocate(aexp)
    if (allocated(omega_t        )) deallocate(omega_t)
    if (allocated(age_univ       )) deallocate(age_univ)
  end subroutine free

end module readTree
