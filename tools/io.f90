module io
  implicit none

  integer :: tmp_unit
  logical :: infos_read = .false.

  type :: MEMBERS_T
     integer, dimension(:), allocatable :: ids
  end type MEMBERS_T

  type :: INFOS_T
     integer                                   :: ncpu, ndim, levelmin, levelmax
     real(kind = 8)                            :: t, aexp, unit_l, unit_t, boxlen
     real(kind = 8), dimension(:), allocatable :: bound_key
  end type INFOS_T

  private :: infos_read, tmp_unit
contains

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

  subroutine read_particle_header (filename, ndim, nparts)
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: nparts, ndim
    open(newunit = tmp_unit, file=filename, status='old', form='unformatted')
    read(tmp_unit) !ncpu
    read(tmp_unit) ndim
    read(tmp_unit) nparts
  end subroutine read_particle_header

  subroutine read_particle_data (ndim, nparts, nstar, pos, vel, m, ids, birth_date)
    integer, intent(in)                          :: ndim, nparts
    real(kind=8), dimension(ndim, nparts), intent(out) :: pos, vel
    integer, intent(out)                         :: nstar
    integer,      dimension(nparts), intent(out) :: ids
    real(kind=8), dimension(nparts), intent(out) :: m, birth_date

    read(tmp_unit) ! ?
    read(tmp_unit) nstar
    read(tmp_unit) ! ?
    read(tmp_unit) ! ?
    read(tmp_unit) ! ?

    read(tmp_unit) pos(1,:)
    read(tmp_unit) pos(2,:)
    read(tmp_unit) pos(3,:)
    read(tmp_unit) vel(1,:)
    read(tmp_unit) vel(2,:)
    read(tmp_unit) vel(3,:)

    read(tmp_unit) m

    read(tmp_unit) ids

    read(tmp_unit) birth_date
    close(tmp_unit)
  end subroutine read_particle_data

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
       pos = pos / (infos%boxlen*infos%unit_l/3.085677581e+24) + 0.5d0
       hlevel(i) = mylevel
       idDM(i) = idh
       mDM(i) = mhalo*1d11
       posDM(:, i) = pos
       Lnorm = sqrt(L(1)*L(1) + L(2)*L(2) + L(3)*L(3))
       LDM(:, i) = L/Lnorm
       rvirDM(i) = rvir / (infos%boxlen*infos%unit_l/3.085677581e+24)
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

  subroutine read_list_header(filename, lines, columns)
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: lines, columns

    open(newunit=tmp_unit, file=filename, form='unformatted')
    read(tmp_unit) lines, columns

  end subroutine read_list_header

  subroutine read_list_data(lines, columns, data)
    integer, intent(in)                                  :: lines, columns
    real(kind=4), dimension(lines, columns), intent(out) :: data

    read(tmp_unit) data
    close(tmp_unit)
  end subroutine read_list_data
end module io
