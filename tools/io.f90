module io
  implicit none

  integer :: tmp_unit
  logical :: infos_read = .false.

  type :: MEMBERS_T
     integer, dimension(:), allocatable :: ids
     integer :: parts
  end type MEMBERS_T

  type :: INFOS_T
     integer                                   :: ncpu, ndim, levelmin, levelmax
     real(kind = 8)                            :: t, aexp, unit_l, unit_t, boxlen
     real(kind = 8), dimension(:), allocatable :: bound_key
  end type INFOS_T

  private :: infos_read, tmp_unit, read_1, read_2, read_2_dummy
contains

  !! Read information of an output, returning an INFOS_T object
  subroutine read_info_headers(filename, infos)
    character(len=*), intent(in)                           :: filename
    type(INFOS_T), intent(out) :: infos

    logical            :: ok
    integer            :: impi, i
    character(len=128) :: ordering

    inquire(file=filename, exist=ok)
    if (.not. ok) then
       print*, filename // ' not found'
       stop
    end if

    open(unit=10, file=filename, form='formatted', status='old')
    read(10, '("ncpu        =",I11)') infos%ncpu
    read(10, '("ndim        =",I11)') infos%ndim
    read(10, '("levelmin    =",I11)') infos%levelmin
    read(10, '("levelmax    =",I11)') infos%levelmax
    read(10, *)
    read(10, *)
    read(10, *)

    read(10, '("boxlen      =",E23.15)') infos%boxlen
    read(10, '("time        =",E23.15)') infos%t
    read(10, '("aexp        =",E23.15)') infos%aexp
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, '("unit_l      =",E23.15)') infos%unit_l
    read(10, *)
    read(10, '("unit_t      =",E23.15)') infos%unit_t

    read(10, *)
    read(10, '("ordering type=",A80)') ordering
    read(10, *)

    if (TRIM(ordering) == 'hilbert') then
       if (.not. allocated(infos%bound_key)) then
          allocate(infos%bound_key(0:infos%ncpu))
       end if

       do impi = 1, infos%ncpu
          read(10, '(I8,1X,E23.15,1X,E23.15)') i, infos%bound_key(impi-1), infos%bound_key(impi)
       end do
    endif
    close(10)
    infos_read = .true.
  end subroutine read_info_headers

  subroutine assert_infos(status)
    integer, intent(out) :: status
    status = 0
    if (.not. infos_read) then
       print*, 'E: you need first to call read_info(…)'
       status = 1
    end if
  end subroutine assert_infos

  !! Read particle file header, returning the number of particules and the number of dimensions
  subroutine read_particle_header (filename, ndim, nparts)
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: nparts, ndim
    open(newunit = tmp_unit, file=filename, status='old', form='unformatted')
    read(tmp_unit) !ncpu
    read(tmp_unit) ndim
    read(tmp_unit) nparts
  end subroutine read_particle_header

  !! Read particles data, returning the particles information
  subroutine read_particle_data (ndim, nparts, nstar, pos, vel, m, ids, birth_date)
    integer, intent(in)                          :: ndim, nparts
    real(kind=8), dimension(ndim, nparts), intent(out) :: pos, vel
    integer, intent(out)                         :: nstar
    integer,      dimension(nparts), intent(out) :: ids
    real(kind=8), dimension(nparts), intent(out) :: m, birth_date

    real(kind=8), dimension(nparts) :: tmp

    read(tmp_unit) ! ?
    read(tmp_unit) nstar
    read(tmp_unit) ! ?
    read(tmp_unit) ! ?
    read(tmp_unit) ! ?

    read(tmp_unit) tmp
    pos(1,:) = tmp
    read(tmp_unit) tmp
    pos(2,:) = tmp
    read(tmp_unit) tmp
    pos(3,:) = tmp
    read(tmp_unit) tmp
    vel(1,:) = tmp
    read(tmp_unit) tmp
    vel(2,:) = tmp
    read(tmp_unit) tmp
    vel(3,:) = tmp

    read(tmp_unit) m

    read(tmp_unit) ids

    read(tmp_unit) birth_date
    close(tmp_unit)
  end subroutine read_particle_data

  !! Read the header of a brick file
  subroutine read_brick_header(filename, infos, nbodies, aexp, age_univ, nb_of_halos, &
       nb_of_subhalos)
    character(len=*), intent(in) :: filename
    type(INFOS_T), intent(in)    :: infos

    integer, intent(out)         :: nbodies, nb_of_subhalos, nb_of_halos
    real(kind=4), intent(out)    :: aexp
    real(kind=4), intent(out)    :: age_univ

    open(newunit=tmp_unit, file=filename, form='unformatted')

    read(tmp_unit) nbodies
    read(tmp_unit) !massp
    read(tmp_unit) aexp
    read(tmp_unit) !omega_t
    read(tmp_unit) age_univ
    read(tmp_unit) nb_of_halos, nb_of_subhalos

    age_univ = age_univ*infos%unit_t

  end subroutine read_brick_header

  !! Read a brick file, returning the relevant quantities of each halo
  subroutine read_brick_data(nb_of_DM, infos, DM_type, &
       & mDM, posDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LDM, idDM, members)

    integer, intent(in)                            :: nb_of_DM
    logical, intent(in)                            :: DM_type
    type(INFOS_T), intent(in)                      :: infos

    real(kind=8), intent(out), dimension(nb_of_DM) :: mDM, rvirDM
    real(kind=8), intent(out), dimension(nb_of_DM) :: mvirDM, TvirDM, hlevel
    real(kind=8), intent(out), dimension(infos%ndim, nb_of_DM) :: LDM, posDM
    integer, intent(out), dimension(nb_of_DM)      :: idDM
    type(MEMBERS_T), dimension(nb_of_DM), intent(out) :: members

    integer                                        :: nb_of_parts, idh, mylevel, hosthalo
    integer                                        :: hostsub, nbsub, nextsub
    real(kind=4) :: mhalo, rvir, mvir, tvir, cvel
    real(kind=4), dimension(infos%ndim) :: pos, L
    real(kind=8)                  :: Lnorm, csound2


    integer :: i, status

    call assert_infos(status)
    if (status > 0) then
       close(tmp_unit)
       return
    end if

    do i = 1, nb_of_DM
       read(tmp_unit) nb_of_parts
       allocate(members(i)%ids(nb_of_parts))
       members(i)%parts = nb_of_parts
       read(tmp_unit) members(i)%ids
       ! Read properties of each halo
       read(tmp_unit) idh
       read(tmp_unit) !timestep
       read(tmp_unit) mylevel, hosthalo, hostsub, nbsub, nextsub

       read(tmp_unit) mhalo
       read(tmp_unit) pos
       read(tmp_unit) !speed
       read(tmp_unit) L
       read(tmp_unit) !r, a, b, c
       read(tmp_unit) !ek, ep, et
       read(tmp_unit) !spin

       if (.not. DM_type) read(tmp_unit) !sigma stuff
       read(tmp_unit) rvir, mvir, tvir, cvel
       read(tmp_unit)
       if (.not. DM_type) then
          read(tmp_unit) !npoints
          read(tmp_unit) !rdum
          read(tmp_unit) !density
       endif

       ! Convert back to adim units

       hlevel(i) = mylevel
       idDM(i) = idh
       mDM(i) = mhalo*1d11
       posDM(:, i) = pos / (infos%boxlen*infos%unit_l/3.08d24)+0.5d0
       Lnorm = sqrt(L(1)*L(1) + L(2)*L(2) + L(3)*L(3))
       LDM(:, i) = L/Lnorm
       rvirDM(i) = rvir / (infos%boxlen*infos%unit_l/3.08d24)
       if(DM_type) then
          mvirDM(i) = mvir*1d11
       else
          mvirDM(i) = mhalo*1d11
       endif
       csound2 = 6.67d-8*mvirDM(i)*2d33/(rvirDM(i)*3.08d24)
       TvirDM(i) = csound2*1.66d-24/1.666666667/1.38d-16

    end do
    close(tmp_unit)
  end subroutine read_brick_data

  !! Read the headers of the list, returning the shape of the data
  subroutine read_list_header(filename, lines, columns)
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: lines, columns

    open(newunit=tmp_unit, file=filename, form='unformatted')
    read(tmp_unit) lines, columns

  end subroutine read_list_header

  !! Read the data contained in the list, returning it
  subroutine read_list_data(lines, columns, data)
    integer, intent(in)                                  :: lines, columns
    real(kind=4), dimension(lines, columns), intent(out) :: data

    read(tmp_unit) data
    close(tmp_unit)
  end subroutine read_list_data

  !! Read the first part of the mergertree headers, returning the number of steps of the simulation
  subroutine read_mergertree_headers_1 (mergertree_file, nsteps)
    character(len=*), intent(in)                    :: mergertree_file
    integer(kind=4), intent(out) :: nsteps

    open(newunit=tmp_unit, file=trim(mergertree_file), form='unformatted')
    read(tmp_unit) nsteps

  end subroutine read_mergertree_headers_1

  !! Read the second part of the header. You have to call read_mergertree_headers_1 BEFORE
  subroutine read_mergertree_headers_2(nb_of_halos, nb_of_subhalos, aexp, omega_t, age_univ, nsteps)
    integer, intent(in) :: nsteps
    integer, intent(out), dimension(nsteps)      :: nb_of_halos, nb_of_subhalos
    real(kind=4), intent(out), dimension(nsteps) :: aexp, omega_t, age_univ
    read(tmp_unit) nb_of_halos, nb_of_subhalos
    read(tmp_unit) aexp
    read(tmp_unit) omega_t
    read(tmp_unit) age_univ
  end subroutine read_mergertree_headers_2

  !! Read the mergertree and return an array containing each halos parent
  subroutine read_mergertree_parent_of (nhalos, nhalos_at_step, nsteps, halos_z0, parent)
    use misc

    integer, intent(in)  :: nhalos, nsteps
    integer, intent(in), dimension(nsteps) :: nhalos_at_step

    integer, intent(out), dimension(nhalos, nsteps) :: parent
    integer, intent(out), dimension(nhalos_at_step(nsteps)) :: halos_z0

    integer, dimension(nhalos) :: current_children
    integer, dimension(nhalos) :: order

    integer(kind=4), dimension(:), allocatable :: idfather
    real(kind=4), dimension(:), allocatable :: mfather
    integer :: i, step, nb_of_fathers, nb_of_sons, halo_id, imax, child

    ! children_sorted = children

    ! call quick_sort(children_sorted, order, nhalos)
    ! print*, children, children_sorted, order
    halos_z0 = 0
    do step = 1, nsteps
       do i = 1, nhalos_at_step(step)
          call read_1(halo_id, nb_of_fathers)
          if (step == nsteps) then
             halos_z0(i) = halo_id
          end if
          allocate(idfather(nb_of_fathers), mfather(nb_of_fathers))
          call read_2(nb_of_fathers, idfather, mfather)
          call max_index(mfather, imax)

          parent(halo_id, step) = idfather(imax)
          deallocate(idfather, mfather)
       end do
    end do

    close(tmp_unit)

  end subroutine read_mergertree_parent_of

  !! Read a mergertree file and compute the positions of each couple halo, step
  !! parameter: halos, integer array of halos to find
  !! parameter: steps, integer array that specify the step in which to find the halo
  !! parameter: pos, vel, real arrays that gives each halo position and velocity
  subroutine read_mergertree_positions (halos, steps, &
       pos, vel,&
       nhalos, nsteps)
    integer, intent(in) :: nhalos, nsteps
    integer, intent(in), dimension(nhalos) :: halos, steps
    real(kind=8), intent(out), dimension(nhalos) :: pos, vel

    integer :: halo_id, nb_of_father, step
    do step = 1, nsteps
       call read_1(halo_id, nb_of_father)
       allocate(idfather(nb_of_fathers), mfather(nb_of_fathers))
       call read_2(nb_of_fathers, idfather, mfather)
       deallocate(idfather, mfather)
    end do

  end subroutine read_mergertree_positions

  subroutine read_1 (halo_id, nb_of_fathers)
    integer, intent(out) :: halo_id, nb_of_fathers
    read(tmp_unit) halo_id

    read(tmp_unit) ! bushid
    read(tmp_unit) ! mystep
    read(tmp_unit) ! leve, hosthalo, hostsub, nbsub, nextsub
    read(tmp_unit) ! m
    read(tmp_unit) ! macc
    read(tmp_unit) ! px, py, pz
    read(tmp_unit) ! vx, vy, vz
    read(tmp_unit) ! Lx, Ly, Lz
    read(tmp_unit) ! r, ra, rb, rc
    read(tmp_unit) ! ek, ep, et
    read(tmp_unit) ! spin

    read(tmp_unit) nb_of_fathers

  end subroutine read_1

  subroutine read_2_dummy (nb_of_fathers)
    integer, intent(in) :: nb_of_fathers
    integer :: nb_of_sons
    if(nb_of_fathers > 0)then
       read(tmp_unit) !idfather
       read(tmp_unit) !mfather
    endif

    read(tmp_unit) nb_of_sons
    if(nb_of_sons>0)then
       read(tmp_unit) !idson(1:nb_of_sons)
    endif

    read(tmp_unit) ! ??
    read(tmp_unit) ! ??
  end subroutine read_2_dummy

  subroutine read_2 (nb_of_fathers, idfather, mfather)
    integer, intent(in) :: nb_of_fathers
    integer(kind=4), intent(out), dimension(nb_of_fathers) :: idfather
    real(kind=4), intent(out), dimension(nb_of_fathers) :: mfather
    integer :: nb_of_sons

    if(nb_of_fathers > 0)then
       read(tmp_unit) idfather
       read(tmp_unit) mfather
    endif
    read(tmp_unit) nb_of_sons
    if(nb_of_sons>0)then
       ! allocate(idson(1:nb_of_sons))
       read(tmp_unit) !idson(1:nb_of_sons)
       ! deallocate(idson)
    endif

    read(tmp_unit) ! ??
    read(tmp_unit) ! ??

  end subroutine read_2

end module io
