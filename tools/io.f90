module io
  implicit none
  real(kind = 8), dimension(:), allocatable :: bound_key
  integer                                   :: ncpu, ndim, levelmin, levelmax
  real(kind = 8)                            :: t, aexp, unit_l, unit_t

  logical :: infos_read = .false.

  private :: infos_read
contains
  subroutine read_info_headers(filename)
    character(len=*), intent(in)                           :: filename

    logical                                                :: ok
    integer                                                :: impi, i
    character(len=80)                                      :: ordering

    inquire(file=filename, exist=ok)
    if (.not. ok) then
       print*, filename // ' not found'
       stop
    end if

    open(unit=10, file=filename, form='formatted', status='old')
    read(10, '("ncpu        =",I11)') ncpu
    read(10, '("ndim        =",I11)') ndim
    read(10, '("levelmin    =",I11)') levelmin
    read(10, '("levelmax    =",I11)') levelmax
    read(10, *)
    read(10, *)
    read(10, *)

    read(10, *)
    read(10, '("time        =",E23.15)') t
    read(10, '("aexp        =",E23.15)') aexp
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, '("unit_l      =",E23.15)') unit_l
    read(10, *)
    read(10, '("unit_t      =",E23.15)') unit_t

    read(10, *)
    read(10, '("ordering type=",A80)') ordering
    read(10, *)

    if (TRIM(ordering) == 'hilbert') then
       if (.not. allocated(bound_key)) then
          allocate(bound_key(0:ncpu))
       end if

       do impi = 1, ncpu
          read(10, '(I8,1X,E23.15,1X,E23.15)') i, bound_key(impi-1), bound_key(impi)
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
    open(unit = 20, file=filename, status='old', form='unformatted')
    read(20) !ncpu
    read(20) ndim
    read(20) nparts
  end subroutine read_particle_header

  subroutine read_particle_data (ndim, nparts, nstar, x, y, z, vx, vy, vz, m, ids, birth_date)
    integer, intent(in)                          :: ndim, nparts
    real(kind=8), dimension(nparts), intent(out) :: x, y, z, vx, vy, vz
    integer, intent(out)                         :: nstar
    integer,      dimension(nparts), intent(out) :: ids
    real(kind=8), dimension(nparts), intent(out) :: m, birth_date

    read(20) ! ?
    read(20) nstar
    read(20) ! ?
    read(20) ! ?
    read(20) ! ?

    read(20) x
    read(20) y
    read(20) z
    read(20) vx
    read(20) vy
    read(20) vz

    read(20) m

    read(20) ids

    read(20) birth_date
    close(20)
  end subroutine read_particle_data

  subroutine read_brick_header(filename, nbodies, aexp, age_univ, nb_of_halos, &
       nb_of_subhalos)
    character(len=*), intent(in) :: filename

    integer, intent(out)         :: nbodies, nb_of_subhalos, nb_of_halos
    real(kind=4), intent(out)    :: aexp
    real(kind=4), intent(out)    :: age_univ

    open(20, file=filename, form='unformatted')

    read(20) nbodies
    read(20)
    read(20) aexp
    read(20)
    read(20) age_univ
    read(20) nb_of_halos, nb_of_subhalos

    age_univ = age_univ*unit_t

  end subroutine read_brick_header

  subroutine read_brick_data(nb_of_DM, DM_type, &
       & mDM, xDM, yDM, zDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LxDM, LyDM, LzDM, idDM)

    integer, intent(in)                            :: nb_of_DM
    logical, intent(in)                            :: DM_type
    real(kind=8), intent(out), dimension(nb_of_DM) :: mDM, xDM, yDM, zDM, rvirDM
    real(kind=8), intent(out), dimension(nb_of_DM) :: mvirDM, TvirDM, hlevel
    real(kind=8), intent(out), dimension(nb_of_DM) :: LxDM, LyDM, LzDM
    integer, intent(out), dimension(nb_of_DM)      :: idDM

    integer                                        :: nb_of_parts, idh, mylevel, hosthalo
    integer                                        :: hostsub, nbsub, nextsub
    real(kind=8)                                   :: mhalo, px, py, pz, Lx, Ly, Lz, rvir, mvir
    real(kind=8)                                   :: xx, yy, zz, drr, tvir, cvel, Lnorm, csound2

    integer :: i, status

    call assert_infos(status)
    if (status > 0) then
       close(20)
       return
    end if

    do i = 1, nb_of_DM
       read(20) nb_of_parts
       read(20) ! read members
       ! Read properties of each halo
       read(20) idh
       read(20)
       read(20) mylevel, hosthalo, hostsub, nbsub, nextsub
       read(20) mhalo
       read(20) px, py, pz
       read(20)
       read(20) Lx, Ly, Lz
       read(20)
       read(20)
       read(20)
       print*, 'Here'
       if (.not. DM_type) read(20) !sigma stuff
       read(20) rvir, mvir, tvir, cvel
       read(20)
       if (.not. DM_type) then
          read(20) !npoints
          read(20) !rdum
          read(20) !density
       endif
       close(20)
       xx = px*1d6*3.08d18/unit_l+0.5d0
       yy = py*1d6*3.08d18/unit_l+0.5d0
       zz = pz*1d6*3.08d18/unit_l+0.5d0
       drr = rvir*1d6*3.08d18/unit_l

       hlevel(i) = mylevel
       idDM(i) = idh
       mDM(i) = mhalo*1d11
       xDM(i) = px
       yDM(i) = py
       zDM(i) = pz
       Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)
       LxDM(i) = Lx/Lnorm
       LyDM(i) = Ly/Lnorm
       LzDM(i) = Lz/Lnorm
       rvirDM(i) = rvir
       if(DM_type) then
          mvirDM(i) = mvir*1d11
       else
          mvirDM(i) = mhalo*1d11
       endif
       csound2 = 6.67d-8*mvirDM(i)*2d33/(rvirDM(i)*3.08d24)
       TvirDM(i) = csound2*1.66d-24/1.666666667/1.38d-16

    end do

  end subroutine read_brick_data

  subroutine free ()
    if (allocated(bound_key)) then
       deallocate(bound_key)
    end if
  end subroutine free

  subroutine read_list_header(filename, lines, columns, unit)
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: lines, columns
    integer, intent(out)         :: unit

    open(newunit=unit, file=filename, form='unformatted')
    read(unit) lines, columns

    ! columns = columns + 1 ! firt column: id
    print*, lines, '*', columns

  end subroutine read_list_header

  subroutine read_list_data(lines, columns, data, unit)
    integer, intent(in)                                  :: lines, columns
    integer, intent(out) :: unit
    real(kind=4), dimension(lines, columns), intent(out) :: data

    read(unit) data
    close(unit)
  end subroutine read_list_data
end module io
