module io
  use types
  use misc
  use compute
  implicit none

  integer :: tmp_unit
  logical :: infos_read = .false.

  !> Type holding the members in a halo
  type :: MEMBERS_T
     integer, dimension(:), allocatable :: ids
     integer :: parts
  end type MEMBERS_T

  !> Type containing data about simulation
  type :: INFOS_T
     integer                                   :: ncpu, ndim, levelmin, levelmax
     character(len=128)                        :: ordering
     real(kind = 8)                            :: t, aexp, unit_l, unit_t, boxlen
     real(kind = 8), dimension(:), allocatable :: bound_key
     character(len=128)                        :: basepath
     integer                                   :: output
  end type INFOS_T

  type :: PARTICLE_DATA
     real(kind=8), dimension(:, :), allocatable   :: pos, vel
     real(kind=8), dimension(:), allocatable      :: m
     real(kind=4), dimension(:), allocatable      :: birth_date
     integer, dimension(:), allocatable           :: ids
     integer :: cpu
  end type PARTICLE_DATA

  private :: infos_read, tmp_unit, read_1, read_2, read_2_dummy,&
       & write_list_reals, write_list_ints

  !> Interface to read lists
  interface read_list_data
     subroutine read_list_data_ints(unit, lines, columns, data)
       integer, intent(in)                             :: unit, lines, columns
       integer, dimension(lines, columns), intent(out) :: data
     end subroutine read_list_data_ints

     subroutine read_list_data_reals(unit, lines, columns, data)
       integer, intent(in)                                  :: unit, lines, columns
       real(kind=4), dimension(lines, columns), intent(out) :: data
     end subroutine read_list_data_reals

     subroutine read_list_data_reals8(unit, lines, columns, data)
       integer, intent(in)                                  :: unit, lines, columns
       real(kind=8), dimension(lines, columns), intent(out) :: data
     end subroutine read_list_data_reals8

  end interface read_list_data

  !> Interface to write lists
  interface write_list
     subroutine write_list_reals(filename, lines, columns, data)
       integer, intent(in) :: lines, columns
       character(len=*), intent(in) :: filename

       real(kind=4), intent(in), dimension(lines, columns) :: data
     end subroutine write_list_reals

     subroutine write_list_ints(filename, lines, columns, data)
       integer, intent(in) :: lines, columns
       character(len=*), intent(in) :: filename

       integer, intent(in), dimension(lines, columns) :: data
     end subroutine write_list_ints
  end interface write_list


contains

  !! Read information of an output, returning an INFOS_T object
  subroutine read_info_headers(basepath, output, infos)
    character(len=*), intent(in) :: basepath
    integer, intent(in)          :: output

    class(INFOS_T), intent(out) :: infos

    character(len=256) :: filename

    logical            :: ok
    integer            :: impi, i, unit

    write(filename, '(a,a,i0.5,a,i0.5,a)') trim(basepath), &
         '/output_', output, '/info_',&
         output, '.txt'

    print*, trim(filename)
    inquire(file=trim(filename), exist=ok)
    if (.not. ok) then
       print*, filename // ' not found'
       stop
    end if

    infos%basepath = basepath
    infos%output   = output

    open(newunit=unit, file=trim(filename), form='formatted', status='old')
    read(unit, '("ncpu        =",I11)') infos%ncpu
    read(unit, '("ndim        =",I11)') infos%ndim
    read(unit, '("levelmin    =",I11)') infos%levelmin
    read(unit, '("levelmax    =",I11)') infos%levelmax
    read(unit, *)
    read(unit, *)
    read(unit, *)

    read(unit, '("boxlen      =",E23.15)') infos%boxlen
    read(unit, '("time        =",E23.15)') infos%t
    read(unit, '("aexp        =",E23.15)') infos%aexp
    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, '("unit_l      =",E23.15)') infos%unit_l
    read(unit, *)
    read(unit, '("unit_t      =",E23.15)') infos%unit_t

    read(unit, *)
    read(unit, '("ordering type=",A80)') infos%ordering
    read(unit, *)

    if (TRIM(infos%ordering) == 'hilbert') then
       if (.not. allocated(infos%bound_key)) then
          allocate(infos%bound_key(0:infos%ncpu))
       end if

       do impi = 1, infos%ncpu
          read(unit, '(I8,1X,E23.15,1X,E23.15)') i, infos%bound_key(impi-1), infos%bound_key(impi)
       end do
    endif
    close(unit)
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

  !! Wrapper around particle reader
  subroutine read_particle(basepath, output, cpu, nstar, pos, vel, m, ids, birth_date, ndim, nparts)
    character(len=*), intent(in) :: basepath
    integer, intent(in)          :: output, cpu

    real(kind=8), dimension(:, :), intent(out), allocatable :: pos, vel
    integer, intent(out)                                    :: nstar, ndim, nparts
    integer,      dimension(:), intent(out), allocatable    :: ids
    real(kind=8), dimension(:), intent(out), allocatable    :: m
    real(kind=4), dimension(:), intent(out), allocatable    :: birth_date

    character(len=10)  :: tmp_char1, tmp_char2
    character(len=100) :: path
    integer            :: unit

    write(tmp_char1, '(i0.5)') output
    write(tmp_char2, '(i0.5)') cpu

    path = trim(basepath) // '/output_' // trim(tmp_char1) // '/part_' // trim(tmp_char1) &
         // '.out' // trim(tmp_char2)
    call read_particle_header(path, ndim, nparts, unit)
    allocate(pos(ndim, nparts))
    allocate(vel(ndim, nparts))
    allocate(ids(nparts))
    allocate(m(nparts))
    allocate(birth_date(nparts))
    call read_particle_data(ndim, nparts, unit, nstar, pos, vel, m, ids, birth_date)

  end subroutine read_particle
  !! Read particle file header, returning the number of particules and the number of dimensions
  subroutine read_particle_header (filename, ndim, nparts, tmp_unit)
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: nparts, ndim, tmp_unit
    open(newunit = tmp_unit, file=filename, status='old', form='unformatted')
    read(tmp_unit) !ncpu
    read(tmp_unit) ndim
    read(tmp_unit) nparts
  end subroutine read_particle_header

  !! Read particles data, returning the particles information
  subroutine read_particle_data (ndim, nparts, tmp_unit, nstar, pos, vel, m, ids, birth_date)
    integer, intent(in)                          :: ndim, nparts, tmp_unit
    real(kind=8), dimension(ndim, nparts), intent(out) :: pos, vel
    integer, intent(out)                         :: nstar
    integer,      dimension(nparts), intent(out) :: ids
    real(kind=8), dimension(nparts), intent(out) :: m
    real(kind=4), dimension(nparts), intent(out) :: birth_date

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
    class(INFOS_T), intent(in)    :: infos

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
    class(INFOS_T), intent(in)                     :: infos

    real(kind=8), intent(out), dimension(nb_of_DM)             :: mDM, rvirDM
    real(kind=8), intent(out), dimension(nb_of_DM)             :: mvirDM, TvirDM, hlevel
    real(kind=8), intent(out), dimension(infos%ndim, nb_of_DM) :: LDM, posDM
    integer, intent(out), dimension(nb_of_DM)                  :: idDM
    class(MEMBERS_T), dimension(nb_of_DM), intent(out)         :: members

    integer                             :: nb_of_parts, idh, mylevel, hosthalo
    integer                             :: hostsub, nbsub, nextsub
    real(kind=4)                        :: mhalo, rvir, mvir, tvir, cvel
    real(kind=4), dimension(infos%ndim) :: pos, L
    real(kind=8)                        :: Lnorm, csound2


    integer :: i, status

    call assert_infos(status)
    if (status > 0) then
       close(tmp_unit)
       return
    end if

    do i = 1, nb_of_DM
       ! if (modulo(i, 1000) == 0) then
       !    write(*, *) i
       ! end if
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
       hlevel(i)   = mylevel
       idDM(i)     = idh
       mDM(i)      = mhalo*1d11 ! Msun
       posDM(:, i) = pos / (infos%unit_l/3.08d24)+0.5d0  ! in [0,1]
       Lnorm       = norm2(L)
       LDM(:, i)   = L/Lnorm
       rvirDM(i)   = rvir / (infos%boxlen*infos%unit_l/3.08d24) ! in [0, 1]
       if(DM_type) then
          mvirDM(i) = mvir*1d11 ! Msun
       else
          mvirDM(i) = mhalo*1d11 ! Msun
       endif
       csound2 = 6.67d-8*mvirDM(i)*2d33/(rvirDM(i)*3.08d24)
       TvirDM(i) = csound2*1.66d-24/1.666666667/1.38d-16

    end do
    close(tmp_unit)
  end subroutine read_brick_data

  !! Read the headers of the list, returning the shape of the data
  subroutine read_list_header(filename, unit, lines, columns)
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: unit, lines, columns

    open(newunit=unit, file=filename, form='unformatted')
    read(unit) lines, columns

  end subroutine read_list_header

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

    integer(kind=4), dimension(:), allocatable :: idfather
    real(kind=4), dimension(:), allocatable :: mfather
    integer :: i, step, nb_of_fathers, halo_id, imax

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
       pos, vel, &
       nhalos_at_step, nhalos, nsteps)
    integer, intent(in) :: nhalos, nsteps
    integer, intent(in), dimension(nhalos+1) :: halos, steps
    integer, intent(in), dimension(nsteps) :: nhalos_at_step


    real(kind=4), intent(out), dimension(nhalos, 3) :: pos, vel

    integer :: halo_id, nb_of_fathers, step, i, j
    real(kind=4), dimension(3) :: tmp_pos, tmp_vel

    pos = 0
    vel = 0
    do step = 1, nsteps
       print*, 'Step:', step, nhalos_at_step(step)
       do i = 1, nhalos_at_step(step)
          call read_1_dynamics(halo_id, tmp_pos, tmp_vel, nb_of_fathers)
          call read_2_dummy(nb_of_fathers)
          do j = 1, nhalos
             ! at some point, halos are all 0
             if (halos(j) == 0) then
                exit
             end if
             if (step == steps(j) .and. halo_id == halos(j)) then
                pos(j, :) = tmp_pos
                vel(j, :) = tmp_vel
                exit
             end if
          end do
       end do
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

  subroutine read_1_dynamics (halo_id, pos, vel, nb_of_fathers)
    integer, intent(out) :: halo_id, nb_of_fathers
    real(kind=4), intent(out), dimension(3) :: pos, vel
    read(tmp_unit) halo_id

    read(tmp_unit) ! bushid
    read(tmp_unit) ! mystep
    read(tmp_unit) ! leve, hosthalo, hostsub, nbsub, nextsub
    read(tmp_unit) ! m
    read(tmp_unit) ! macc
    read(tmp_unit) pos
    read(tmp_unit) vel
    read(tmp_unit) ! Lx, Ly, Lz
    read(tmp_unit) ! r, ra, rb, rc
    read(tmp_unit) ! ek, ep, et
    read(tmp_unit) ! spin

    read(tmp_unit) nb_of_fathers

  end subroutine read_1_dynamics

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

  !> Read the particles in region around center with size
  !! arguments:
  !!
  !! center, size : array, center and size of the region
  !! infos : INFOS_T containing the information about an output
  !! callback: a subroutine to be called for each particle found
  subroutine read_region(center, width, infos, data)
    class(INFOS_T), intent(in)                      :: infos
    real(kind=8), dimension(infos%ndim), intent(in) :: center, width

    type(PARTICLE_DATA), allocatable, dimension(:), intent(out) :: data

    integer, dimension(infos%ncpu) :: cpu_list
    integer                                      :: nstar, ndim, nparts
    real(kind=8), dimension(infos%ndim)          :: X0, X1

    type(PARTICLE_DATA)  :: dt
    integer :: counter, cpu, i, j, k, dim, part_i
    logical, dimension(:), allocatable :: mask
    real(8) :: dist

    X0 = center-width
    X1 = center+width

    call get_cpu_list(X0, X1, infos%levelmax, infos%bound_key, &
         & cpu_list, infos%ncpu, infos%ndim)

    ! count cpus
    counter = 0
    do i = 1, infos%ncpu
       if (cpu_list(i) > 0) then
          counter = counter + 1
       end if
    end do

    ! allocate data for the custom array
    allocate(data(counter))

    ! iterate over each cpus
    !$OMP PARALLEL DO DEFAULT(firstprivate) shared(data)
    do i = 1, infos%ncpu
       cpu = cpu_list(i)
       if (cpu_list(i) == 0) then
          cycle
       end if

       call read_particle(infos%basepath, infos%output, cpu, &
            nstar, dt%pos, dt%vel, dt%m, dt%ids, dt%birth_date, ndim, nparts)
       print*, 'cpu ', cpu, nparts, counter

       ! keep the indexes of the particles within bounds
       allocate(mask(nparts))
       mask = .true.
       do part_i = 1, nparts
          do dim = 1, ndim
             ! rule out particles outside region
             ! if the distance is larger than width

             dist = dt%pos(dim, part_i) - center(dim)
             ! correct distances if farther than 0.5 units
             if (dist > 0.5) then
                dist = dist - 1d0
                dt%pos(dim, part_i) = dt%pos(dim, part_i) - 1d0
             else if (dist < -0.5) then
                dist = dist + 1d0
                dt%pos(dim, part_i) = dt%pos(dim, part_i) + 1d0
             else
                dist = abs(dist)
             end if

             if (dist > width(dim)) then
                mask(part_i) = .false.
             end if

          end do
       end do

       ! get the size of the mask
       counter = 0
       do j = 1, nparts
          if (mask(j)) counter = counter + 1
       end do

       ! allocate data
       allocate(data(i)%vel(infos%ndim, counter), &
            data(i)%pos(infos%ndim, counter), &
            data(i)%ids(counter), data(i)%m(counter), &
            data(i)%birth_date(counter))

       !----------------------------------------
       ! Copy the data
       !----------------------------------------
       k = 1
       do j = 1, nparts
          if (mask(j)) then
             data(i)%vel(:, k)     = dt%vel(:, j)
             data(i)%pos(:, k)     = dt%pos(:, j)
             data(i)%ids(k)        = dt%ids(j)
             data(i)%m(k)          = dt%m(j)
             data(i)%birth_date(k) = dt%birth_date(j)

             k = k + 1
          end if
       end do
       data(i)%cpu = cpu

       deallocate(mask)

    end do

  end subroutine read_region
end module io

!-------------------------------------
! Subroutines to interface
!-------------------------------------
!> Read the data contained in the list, returning it
subroutine read_list_data_reals(unit, lines, columns, data)
  integer, intent(in)                                  :: unit, lines, columns
  real(kind=4), dimension(lines, columns), intent(out) :: data

  read(unit) data
  close(unit)
end subroutine read_list_data_reals

subroutine read_list_data_reals8(unit, lines, columns, data)
  integer, intent(in)                                  :: unit, lines, columns
  real(kind=8), dimension(lines, columns), intent(out) :: data

  read(unit) data
  close(unit)
end subroutine read_list_data_reals8

subroutine read_list_data_ints(unit, lines, columns, data)
  integer, intent(in)                             :: unit, lines, columns
  integer, dimension(lines, columns), intent(out) :: data

  read(unit) data
  close(unit)
end subroutine read_list_data_ints

!! Write a list file
subroutine write_list_reals(filename, lines, columns, data)
  integer, intent(in) :: lines, columns
  character(len=*), intent(in) :: filename

  real (kind=4), intent(in), dimension(lines, columns) :: data

  integer :: unit

  open(newunit=unit, file=trim(filename), form='unformatted')

  write(unit) lines, columns
  write(unit) data
  close(unit)
end subroutine write_list_reals

subroutine write_list_ints(filename, lines, columns, data)
  integer, intent(in) :: lines, columns
  character(len=*), intent(in) :: filename

  integer, intent(in), dimension(lines, columns) :: data

  integer :: unit

  open(newunit=unit, file=trim(filename), form='unformatted')

  write(unit) lines, columns
  write(unit) data
  close(unit)
end subroutine write_list_ints
