module readtree
  implicit none
  public
  integer :: nsteps !> Number of steps in the simulation

  integer, dimension(:), allocatable       :: nb_of_halos    !> halos in the tree
  integer, dimension(:), allocatable       :: nb_of_subhalos !> subhalos in the tree
  real (kind=4), dimension(:), allocatable :: aexp     !> expansion factor at each step
  real (kind=4), dimension(:), allocatable :: age_univ !> age of the universe at each step
  real (kind=4), dimension(:), allocatable :: omega_t  !> omega_t at at each step

  ! Internal variables
  integer :: istep = 1
  integer :: j = 1
  integer :: unit

  ! Some temp variable
  integer, allocatable       :: tmp_id_fathers(:), tmp_id_sons(:)
  real(kind=4), allocatable  :: tmp_m_fathers(:)

  private :: istep, j, unit, tmp_id_fathers, tmp_id_sons, tmp_m_fathers
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
       nbsub, nextsub, nb_of_fathers, nb_of_sons, &
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
    integer(kind=4), intent(out), optional :: nb_of_sons
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
          call close_file()
          stop_now = .true.
          return
       end if
    end if

    ! Check that the halo_id is non negative (if it is, stop immediately)
    read(unit) halo_id
    if (halo_id == -1) then
       call close_file()
       stop_now = .true.
       return
    end if
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
       if (allocated(tmp_id_fathers)) deallocate(tmp_id_fathers)
       if (allocated(tmp_m_fathers )) deallocate(tmp_m_fathers)

       allocate(tmp_id_fathers(nb_of_fathers), tmp_m_fathers(nb_of_fathers))

       read(unit) tmp_id_fathers
       read(unit) tmp_m_fathers
    endif

    read(unit) nb_of_sons
    if(nb_of_sons>0)then
       if (allocated(tmp_id_sons)) then
          deallocate(tmp_id_sons)
       end if
       allocate(tmp_id_sons(nb_of_sons))

       read(unit) tmp_id_sons
    endif

    read(unit)
    read(unit)
    ! if(star) read(unit)
  end subroutine iter_tree

  !> Get the father of the last red halo
  subroutine get_fathers(nb_of_fathers, id_fathers, m_fathers)
    integer, intent(in)                                 :: nb_of_fathers !> number of son of the halo
    integer, dimension(nb_of_fathers), intent(out)      :: id_fathers     !> id of the father
    real(kind=4), dimension(nb_of_fathers), intent(out) :: m_fathers      !> mass of the father

    id_fathers = tmp_id_fathers
    m_fathers = tmp_m_fathers
  end subroutine get_fathers

  !> Get the son of the last red halo
  subroutine get_sons(nb_of_sons, id_sons)
    integer, intent(in)                         :: nb_of_sons!> number of son of the halo
    integer, dimension(nb_of_sons), intent(out) :: id_sons    !> ids of the son of the halo

    id_sons = tmp_id_sons
  end subroutine get_sons

  !> Close the tree file
  subroutine close_file()
    logical :: op
    ! TODO: inquire it is opened
    inquire(unit, opened=op)
    if (op) close(unit)
  end subroutine close_file

  !> Free allocated variables
  subroutine free()
    if (allocated(nb_of_halos    )) deallocate(nb_of_halos)
    if (allocated(nb_of_subhalos )) deallocate(nb_of_subhalos)
    if (allocated(aexp           )) deallocate(aexp)
    if (allocated(omega_t        )) deallocate(omega_t)
    if (allocated(age_univ       )) deallocate(age_univ)
  end subroutine free

  !> Revert the tree (from end to begin)
  subroutine revert_steps(filetree, outfich, star)
    character(len=*), intent(in) :: outfich, filetree
    logical, intent(in)            :: star

    integer                                    :: unit_in, unit_out
    integer                                    :: ntot, i, j, istep, idone
    integer(kind=4)                            :: nsteps, nbodies, nb_of_fathers, nb_of_sons, ncont
    integer(kind=4)                            :: mynumber, BushID, mystep, level, hosthalo, hostsub, nbsub, nextsub
    real(kind=4)                               ::m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz, r, ra, rb, rc, ek, ep, et, spin
    real(kind=4)                               ::rvir, mvir, tvir, cvel, rho0, rc2
    integer(kind=4), allocatable, dimension(:) :: nb_of_halos, nb_of_subhalos
    integer(kind=4), allocatable, dimension(:) :: idfather, idson
    real(kind=4),    allocatable, dimension(:) :: aexp, age_univ, omega_t, mfather
    integer(kind=4), allocatable, dimension(:) :: nb_of_halos_rev, nb_of_subhalos_rev
    real(kind=4),    allocatable, dimension(:) :: aexp_rev, age_univ_rev, omega_t_rev

    !-----------------------------------------------
    ! Lecture du fichier contenant les galaxies
    !-----------------------------------------------
    write(*, *) filetree

    open(newunit=unit_in, file=filetree, status='old', form='unformatted')
    read(unit_in) nsteps
    allocate(nb_of_halos(1:nsteps), nb_of_subhalos(1:nsteps), aexp(1:nsteps), omega_t(1:nsteps), age_univ(1:nsteps))

    read(unit_in) nb_of_halos, nb_of_subhalos
    read(unit_in) aexp
    read(unit_in) omega_t
    read(unit_in) age_univ
    write(*, *) nsteps, ' steps'
    write(*, *) 'list of timesteps:'
    do i=1, nsteps
       write(*, *) 1d0/aexp(i)-1d0, age_univ(i)
    enddo
    write(*, *) MINVAL(nb_of_halos), MAXVAL(nb_of_halos), ' galaxies'
    write(*, *) MINVAL(nb_of_subhalos), MAXVAL(nb_of_subhalos), ' sub-galaxies'
    close(unit_in)

    allocate(nb_of_halos_rev(1:nsteps), nb_of_subhalos_rev(1:nsteps))
    allocate(aexp_rev(1:nsteps), omega_t_rev(1:nsteps), age_univ_rev(1:nsteps))
    do i=1, nsteps
       nb_of_halos_rev(nsteps-i+1)    = nb_of_halos(i)
       nb_of_subhalos_rev(nsteps-i+1) = nb_of_subhalos(i)
       aexp_rev(nsteps-i+1)     = aexp(i)
       omega_t_rev(nsteps-i+1)  = omega_t(i)
       age_univ_rev(nsteps-i+1) = age_univ(i)
    enddo

    open(newunit=unit_out, file=trim(outfich), form='unformatted')
    write(unit_out) nsteps
    write(unit_out) nb_of_halos_rev, nb_of_subhalos_rev
    write(unit_out) aexp_rev
    write(unit_out) omega_t_rev
    write(unit_out) age_univ_rev
    deallocate(nb_of_halos_rev, nb_of_subhalos_rev, aexp_rev, omega_t_rev, age_univ_rev)

    idone = nsteps

    do while(idone > 0)

       open(newunit=unit_in, file=trim(filetree), status='old', form='unformatted')
       read(unit_in) nsteps
       read(unit_in) nb_of_halos, nb_of_subhalos
       read(unit_in) aexp
       read(unit_in) omega_t
       read(unit_in) age_univ

       do istep=1, idone

          ntot = nb_of_halos(istep) + nb_of_subhalos(istep)

          if (istep == idone) then

             do j = 1, ntot
                read(unit_in) mynumber
                read(unit_in) BushID
                read(unit_in) mystep
                read(unit_in) level, hosthalo, hostsub, nbsub, nextsub
                read(unit_in) m
                read(unit_in) macc
                read(unit_in) px, py, pz
                read(unit_in) vx, vy, vz
                read(unit_in) Lx, Ly, Lz
                read(unit_in) r, ra, rb, rc
                read(unit_in) ek, ep, et
                read(unit_in) spin

                write(unit_out) mynumber
                write(unit_out) BushID
                write(unit_out) mystep
                write(unit_out) level, hosthalo, hostsub, nbsub, nextsub
                write(unit_out) m
                write(unit_out) macc
                write(unit_out) px, py, pz
                write(unit_out) vx, vy, vz
                write(unit_out) Lx, Ly, Lz
                write(unit_out) r, ra, rb, rc
                write(unit_out) ek, ep, et
                write(unit_out) spin

                read(unit_in) nb_of_fathers
                write(unit_out) nb_of_fathers
                if (nb_of_fathers > 0) then
                   allocate(idfather(nb_of_fathers), mfather(nb_of_fathers))
                   read(unit_in)   idfather
                   read(unit_in)   mfather


                   write(unit_out) idfather
                   write(unit_out) mfather
                   deallocate(idfather, mfather)
                endif

                read(unit_in) nb_of_sons
                write(unit_out) nb_of_sons

                if (nb_of_sons > 0) then
                   allocate(idson(nb_of_sons))
                   read(unit_in)   idson

                   write(unit_out) idson
                   deallocate(idson)
                endif

                read(unit_in) rvir, mvir, tvir, cvel
                read(unit_in) rho0, rc2
                if (star) read(unit_in) ncont

                write(unit_out) rvir, mvir, tvir, cvel
                write(unit_out) rho0, rc2
                if (star) write(unit_out) ncont
             enddo

             write(*, *)'wrote ', idone, 1d0/aexp(idone)-1d0
             idone = idone-1

          else if (istep < idone) then

             do j=1, ntot
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)
                read(unit_in)

                read(unit_in) nb_of_fathers
                if (nb_of_fathers > 0) then
                   allocate(idfather(1:nb_of_fathers), mfather(1:nb_of_fathers))
                   read(unit_in) idfather
                   read(unit_in) mfather
                   deallocate(idfather, mfather)
                endif

                read(unit_in) nb_of_sons
                if (nb_of_sons > 0) then
                   allocate(idson(1:nb_of_sons))
                   read(unit_in) idson
                   deallocate(idson)
                endif

                read(unit_in)
                read(unit_in)
                if (star) read(unit_in)
             enddo

          endif

       enddo
       close(unit_in)

    enddo
    close(unit_out)


  end subroutine revert_steps
end module readtree
