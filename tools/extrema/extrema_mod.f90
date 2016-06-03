!   this is old cartezian extremum code adapted to F90

!   May  , 2012, Dmitri Pogosyan    -   change to planar input of 3D boxes
!   March, 2012, Dmitri Pogosyan    -   first rewrite of an old code
!
!   Input:    dt(:)        - field map
!             nd           - number of dimensions
!             n1,n2(,n3)   - box dimensions
!
!   Output:   ext%pix      - pixel position of extrema
!             ext%ang(2)   - angular position (more accurate than pix center)
!             ext%hes(2,2) - hessian
!             ext%eig(2)   - eigenvalues
!             ext%type     - min, max, saddle
!

MODULE extrema_mod
  use extrema_types
  use omp_lib

  IMPLICIT NONE
  PRIVATE

  type NEIGH_DATA
     INTEGER(I8B)      :: pix
     REAL(DP)          :: val
     INTEGER(I4B)      :: xyz(3)
  end type NEIGH_DATA

  type EXT_META
     integer(I8B), allocatable :: l_map(:)
     real(DP),     allocatable :: CNA(:, :), AtCNA(:, :)
     integer(I8B) :: JITTER=0
  end type EXT_META


  integer(I4B), allocatable :: mod_ikvadr(:, :)
  public  :: FIND_EXTREMA

CONTAINS


  subroutine FIND_EXTREMA(dt, nn, nd, ext, ctrl)

    real(SP),             intent(in), dimension(0:) :: dt
    integer(I4B),         intent(in)                :: nd, nn(:)
    type(EXT_DATA),       intent(inout), optional   :: ext(:)
    type(CND_CNTRL_TYPE), intent(in), optional      :: ctrl

    integer(I8B)            :: ic, n_ext, n_ext_low, n_ext_up
    integer(I4B)            :: nparam, nneigh
    type(NEIGH_DATA), save, allocatable :: all_neighbour_list(:, :) ! list of the neighours of each points
    type(NEIGH_DATA)        :: neighbour_list(3**nd - 1), neighbour_list_new(3**nd - 1)
    type(EXT_DATA)          :: extc
    type(EXT_META), save    :: extmeta

    logical                 :: ifextremum, ifjustprint
    real(SP)                :: dtc
    real(DP)                :: bfit(nd*(nd+3)/2)
    real(DP)                :: am(nd, nd), vm(nd), xm(nd)
    integer(I8B)            :: NCHUNK, NCHUNKE
    integer(I4B)            :: NPROC, OMP_GET_THREAD_NUM
    integer :: i

    logical :: stop_now

    logical, save :: firstCall = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Data read-in and setup

    NPIX = product(int(nn(1:nd), I8B))

    if (allocated(extmeta%l_map)) deallocate(extmeta%l_map)
    allocate(extmeta%l_map(0:NPIX-1))     ! Allocate index map
    extmeta%l_map = 0


    if ( PRESENT(ext) ) then
       n_ext_low   = lbound(ext, 1)
       n_ext_up    = ubound(ext, 1)
       ext(:)%typ = -1
       ext(:)%pix = 0
       do i = 1, 3
          ext(:)%eig(i) = 0
          ext(:)%pos(i) = 0
       end do
       ext(:)%val = 0

       ifjustprint = .false.
    else
       ! write(0, *) 'array to store extrema is not supplied, just print out'
       ifjustprint  = .true.
    endif

    nneigh = 3**nd - 1
    nparam = nd*(nd+3)/2

    if (firstCall) then
       allocate(all_neighbour_list(nneigh, 0:NPIX-1))
       do ic = 0, NPIX-1
          call preset_neighbours(nd, all_neighbour_list(:, ic), nneigh)
       end do
       ! neighbour_list = all_neighbour_list(:, 1)
       ! call set_basis(nd, neighbour_list(:, 1), nneigh, nparam, extmeta)
    end if

    call preset_neighbours(nd, neighbour_list, nneigh)
    call set_basis(nd, neighbour_list, nneigh, nparam, extmeta)

    if ( PRESENT(ctrl) ) then
       NPROC       = ctrl%nproc
       if (NPROC == -1) NPROC = OMP_GET_MAX_THREADS()
       ifjustprint = ctrl%justprint
    endif
    NCHUNK  = NPIX/NPROC
    NCHUNKE = (n_ext_up - n_ext_low + 1)/NPROC

    n_ext = 0
    stop_now = .false.

    ! Precompute the neighbours of each point
    if (firstCall) then
       ! print*, 'Precomputing neighbours'
       do ic = 0, NPIX-1
          call set_current_neighbours(dt, ic, nn, nd, all_neighbour_list(:, ic), nneigh)
       end do
    end if

    !$OMP PARALLEL DEFAULT(SHARED)                                                 &
    !$OMP private(dtc, ic, bfit, am, vm, xm, ifextremum, extc, neighbour_list_new) &
    !$OMP FIRSTPRIVATE(n_ext, neighbour_list) NUM_THREADS(NPROC)
    !$OMP DO SCHEDULE(DYNAMIC, NCHUNK)
    do ic = 0, NPIX-1
       ! Early break if stop_now flag is true
       if (stop_now) cycle

       ! ! NEW CODE
       ! neighour_list = all_neighbour_list(:, ic)

       ! OLD CODE
       neighbour_list_new = all_neighbour_list(:, ic)
       call set_current_neighbours(dt, ic, nn, nd, neighbour_list, nneigh)

       if (any(neighbour_list_new(:)%pix /= neighbour_list(:)%pix)) then
          write(*,*) 'f**k on pix'
       end if
       if (any(neighbour_list_new(:)%xyz(1) /= neighbour_list(:)%xyz(1))) then
          write(*,*) 'f**k on xyz(1)'
       end if
       if (any(neighbour_list_new(:)%xyz(2) /= neighbour_list(:)%xyz(2))) then
          write(*,*) 'f**k on xyz(2)'
       end if
       if (any(neighbour_list_new(:)%xyz(3) /= neighbour_list(:)%xyz(3))) then
          write(*,*) 'f**k on xyz(3)'
       end if
       if (any(neighbour_list_new(:)%val /= neighbour_list(:)%val)) then
          write(*,*) 'f**k on val'
       end if

       ! fit quadratic to the neightbours
       call quadratic_fit(neighbour_list, nneigh, nparam, bfit, extmeta)
       call setmatrix(bfit, am, vm, nd)
       call findextremum(am, vm, xm, nd)
       ifextremum = checkextremum(xm, ic, nn, nd, neighbour_list, nneigh, extmeta)
       if ( ifextremum ) then               ! do postprocessing
          dtc = dt(ic)
          call setmatrix(bfit, am, vm, nd)
          call extremum_properties(ic, dtc, nn, nd, am, vm, xm, extc)
          call markextremum(ic, neighbour_list, nneigh, extc%typ, extmeta)
          if ( ifjustprint ) then
             ! write(*, '(3(f7.2,1x),4(e10.3,1x),i1)') extc%pos, extc%eig, extc%val, extc%typ
          else
             if ( n_ext >= NCHUNKE ) then
                write(0, *) 'Run out of output storage', n_ext, NCHUNKE, n_ext_low, n_ext_up, 'exiting'
                stop_now = .true.
             endif
             ext(n_ext_low + OMP_GET_THREAD_NUM()*NCHUNKE + n_ext) = extc
             n_ext = n_ext + 1
          endif
       endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(extmeta%l_map)
    firstCall = .false.
    return
  END SUBROUTINE FIND_EXTREMA

  !! Fill the neighbours with the positions of the neighbours (except for central one)
  SUBROUTINE preset_neighbours(nd, neighbour_list, nneigh)
    integer(I4B),     intent(in)    :: nd, nneigh
    TYPE(NEIGH_DATA), intent(out)   :: neighbour_list(:)

    integer(I4B), dimension(nd)     :: nn_neigh
    integer(I8B)                    :: i
    integer(I4B)                    :: nc

    nn_neigh = 3
    nc = 1
    do i = 0, nneigh
       if ( 2*i /= nneigh ) then
          neighbour_list(nc)%xyz(1:nd) = index_to_grid(i, nn_neigh, nd) - 1
          nc = nc+1
       endif
    enddo
    return
  END SUBROUTINE preset_neighbours

  !!
  subroutine set_basis(nd, neighbour_list, nneigh, nparam, extmeta)
    integer(I4B), intent(in)      :: nd, nneigh, nparam
    type(NEIGH_DATA), intent(in)  :: neighbour_list(:)
    type(EXT_META), intent(inout) :: extmeta


    REAL(DP),  dimension(nneigh, nparam)  :: AA
    REAL(DP),  dimension(nneigh, nneigh)  :: CNpp

    integer(I4B)                         :: i, j, ic

    if (allocated(extmeta%CNA)) deallocate(extmeta%CNA)
    if (allocated(extmeta%AtCNA)) deallocate(extmeta%AtCNA)
    allocate(extmeta%CNA(nneigh, nparam))
    allocate(extmeta%AtCNA(nparam, nparam))

    ! Set basis
    do i = 1, nd
       AA(:, i) = neighbour_list(:)%xyz(i)
       AA(:, i+nd) = 0.5_dp*AA(:, i)**2
    enddo

    ic=1
    do j = 2, nd
       do i = 1, j-1
          AA(:, 2*nd+ic) = AA(:, i)*AA(:, j)
          ic = ic+1
       enddo
    enddo

    ! Set weights
    CNpp = 0.d0
    forall(i=1:nneigh) CNpp(i, i) = 1.d0/SUM(neighbour_list(i)%xyz(1:nd)**2)
    !    forall(i=1:nneigh) CNpp(i, i) = 1.d0

    call DSYMM('L', 'L', nneigh, nparam, 1.d0, CNpp, nneigh, AA, nneigh, 0.d0, &
         extmeta%CNA, nneigh)
    call DGEMM('T', 'N', nparam, nparam, nneigh, 1.d0, AA, nneigh, extmeta%CNA, nneigh, 0.d0,&
         extmeta%AtCNA, nparam)
    return
  END SUBROUTINE set_basis

  !> Set the neighbours around a cell
  SUBROUTINE set_current_neighbours(dt, icell, nn, nd, neighbour_list, nneigh)
    real(SP),      intent(in), dimension(0:)      :: dt
    integer(I8B),  intent(in)                     :: icell
    integer(I4B),  intent(in), dimension(:)       :: nn
    integer(I4B),  intent(in)                     :: nd, nneigh
    type(NEIGH_DATA), intent(inout)          :: neighbour_list(:)


    integer(I4B), dimension(nd)            :: ijkc, ijk
    integer(I8B)                           :: ineigh
    integer(I4B)                           :: i

    ijkc = index_to_grid(icell, nn, nd)

    do i = 1, nneigh
       ijk = ikvadr(neighbour_list(i)%xyz + ijkc, nn)
       ineigh = grid_to_index(ijk, nn, nd)
       neighbour_list(i)%pix = ineigh
       neighbour_list(i)%val = dt(ineigh) - dt(icell)
    enddo
    return
  END SUBROUTINE set_current_neighbours

  subroutine precompute_ikvadr(nn, ikvadr, lbnd, ubnd)
    integer(I4B), intent(in) :: nn(:)
    integer, intent(in) :: lbnd, ubnd
    integer(I4B), intent(inout) :: ikvadr(1:size(nn),lbnd:ubnd)

    integer(I4B) :: i, d

    do d = 1, size(nn)
       do i = -nn(d), 2*nn(d)
          ikvadr(d, i) = modulo(i, nn(d))
       end do
    end do
  end subroutine precompute_ikvadr

  FUNCTION ikvadr(ijk, nn)
    integer(I4B), intent(in) :: ijk(:), nn(:)
    integer(I4B)             :: ikvadr(size(nn))
    integer(I4B)             :: ikvadr2(size(nn))

    integer(I4B) :: d

    logical,      save :: firstCall = .true.
    integer(I4B), save, allocatable :: ikvadr_arr(:, :)

    if (firstCall) then
       !$OMP CRITICAL
       if (allocated(ikvadr_arr)) deallocate(ikvadr_arr)
       allocate(ikvadr_arr(1:size(nn), -maxval(nn):2*maxval(nn)))


       firstCall = .false.
       call precompute_ikvadr(nn, ikvadr_arr, lbound(ikvadr_arr, 2), ubound(ikvadr_arr, 2))
       !$OMP END CRITICAL
    end if

    do d = 1, size(nn)
       ikvadr(d) = ikvadr_arr(d, ijk(d))
    end do

    ! where (ijk < 0)
    !    ikvadr2 = ijk + nn
    ! elsewhere (ijk >= nn)
    !    ikvadr2 = ijk - nn
    ! elsewhere
    !    ikvadr2 = ijk
    ! endwhere

    ! ikvadr = ikvadr2
    ! if (any(ikvadr /= ikvadr2)) write(*, *) 'F**k'
    return
  END FUNCTION ikvadr

  FUNCTION fkvadr(xyz, nn)
    real(DP),     intent(in) :: xyz(:)
    integer(I4B), intent(in) :: nn(:)
    real(DP)                 :: fkvadr(size(nn))

    where (xyz < -0.5d0 )
       fkvadr = xyz + nn - 1.d-5
    elsewhere (xyz >= nn-0.5d0)
       fkvadr = xyz - nn + 1.d-5
    elsewhere
       fkvadr = xyz
    endwhere
    return
  END FUNCTION fkvadr

  FUNCTION rkvadr(xyz, nn)
    real(DP),     intent(in) :: xyz(:)
    integer(I4B), intent(in) :: nn(:)
    real(DP)                 :: rkvadr(size(nn))

    integer                  :: i

    do i=1, size(nn)
       if (xyz(i) < -0.5d0 ) then
          rkvadr(i) = xyz(i) + nn(i)
       elseif (xyz(i) >= nn(i)-0.5d0) then
          rkvadr = xyz(i) - nn(i)
       else
          rkvadr(i) = xyz(i)
       endif
    enddo
    return
  END FUNCTION rkvadr

  FUNCTION grid_to_index(ijk, nn, nd)
    integer(I4B), intent(in), dimension(:) :: ijk, nn
    integer(I4B), intent(in)               :: nd
    integer(I8B)                           :: grid_to_index

    integer                           :: i

    integer(I8B), save, allocatable :: prod_dims(:)
    logical, save :: firstCall = .true.

    if (firstCall) then
       allocate(prod_dims(2:size(nn)))
       do i = 2, nd
          prod_dims(i) = product(int(nn(1:i-1)))
       end do
       firstCall = .false.
    end if


    grid_to_index = ijk(1)
    do i = 2, nd
       grid_to_index = grid_to_index + ijk(i)*prod_dims(i)
       ! if (product(int(nn(1:i-1), I8B)) /= prod_dims(i)) write(*, *) 'f**k 2'
    enddo

    return
  END FUNCTION grid_to_index

  FUNCTION index_to_grid(ic, nn, nd)
    integer(I8B), intent(in)                :: ic
    integer(I4B), intent(in), dimension(:)  :: nn
    integer(I4B), intent(in)                :: nd
    integer(I4B),             dimension(nd) :: index_to_grid

    integer(I8B)                            :: icell, ibase
    integer(I4B)                            :: i

    icell=ic
    do i = nd, 2, -1
       ibase=PRODUCT(nn(2:i))
       index_to_grid(i) = icell/ibase
       icell = icell - index_to_grid(i)*ibase
    enddo
    index_to_grid(1) = icell
    return
  END FUNCTION index_to_grid

  subroutine quadratic_fit(neighbour_list, nneigh, nparam, bfit, extmeta)
    INTEGER(I4B),     intent(in)  :: nneigh, nparam
    TYPE(NEIGH_DATA), intent(in)  :: neighbour_list(:)
    real(DP),         intent(out) :: bfit(:)
    type(EXT_META), intent(inout) :: extmeta

    real(dp), allocatable         :: tmp_val(:)

    INTEGER(I4B)                               :: i, INFO
    REAL(DP), DIMENSION(40)                    :: WORK
    INTEGER(I4B), DIMENSION(nparam)            :: IPIV
    REAL(DP), DIMENSION(nparam, nparam)         :: AtCNA_loc

    AtCNA_loc = extmeta%AtCNA

    allocate(tmp_val(size(neighbour_list)))
    tmp_val = neighbour_list(:)%val

    ! bfit := 1*CNA'*val.
    call DGEMV('T', nneigh, nparam, 1._dp, extmeta%CNA, nneigh, tmp_val,&
         1, 0._dp, bfit, 1)
    deallocate(tmp_val)

    ! solves AtCNA_loc * X = bfit, thus getting the location
    call DSYSV('L', nparam, 1, AtCNA_loc, nparam, IPIV, bfit, nparam, WORK, 40, INFO)

    return
  END SUBROUTINE quadratic_fit


  !> Solves A*X = -v, where
  !! :param A the input array
  !! :param v the value
  SUBROUTINE findextremum(a, v, x, nd)
    REAL(DP), intent(inout)  :: a(:, :), v(:)
    REAL(DP), intent(out)    :: x(:)
    INTEGER(I4B), intent(in) :: nd

    INTEGER(I4B)        :: INFO, IPIV(nd)
    integer, save       :: lwork = 20
    real(DP)            :: WORK(lwork)
    ! as it is, matrix 'a' is destroyed and v is overwritten

    x = -v
    call DSYSV( 'L', nd, 1, a, nd, IPIV, x, nd, WORK, lwork, INFO )
    lwork = WORK(1)

    return
  END SUBROUTINE findextremum


  logical function checkextremum(x, icell, nn, nd, neighbour_list, nneigh, extmeta)
    REAL(DP),     intent(in), dimension(:)     :: x
    INTEGER(I8B), intent(in)                   :: icell
    INTEGER(I4B), intent(in)                   :: nn(:), nd, nneigh
    type(NEIGH_DATA), intent(in)               :: neighbour_list(:)
    type(EXT_META), intent(inout) :: extmeta

    INTEGER(I8B)                           :: ic
    INTEGER(I4B)                           :: i

    if (extmeta%l_map(icell) > 0) then
       checkextremum = .false.
       return
    endif

    ! Find extrema coordinates on a grid, and the nearest grid point
    ic = grid_to_index(nint(fkvadr(x+index_to_grid(icell, nn, nd), nn)), nn, nd)

    if ( ic == icell ) then
       checkextremum = .true.
    else
       ! checking for 'jitter'
       checkextremum = .false.
       do i = 1, nneigh
          if ( ic == neighbour_list(i)%pix ) then
             !             if (icell == 11100) write(0, *)'lmap', ic, l_map(ic), checkextremum
             if ( extmeta%l_map(ic) == -icell-1 ) then
                checkextremum = .true.
                !                JITTER = JITTER + 1           ! just counter for info
             else
                extmeta%l_map(icell) = -ic-1
             endif
             exit
          endif
       enddo
    endif
    return

  END FUNCTION checkextremum

  subroutine markextremum(vert, neighbour_list, nneigh, typ, extmeta)
    INTEGER(I8B), intent(in)          :: vert
    INTEGER(I4B), intent(in)          :: nneigh, typ
    TYPE(NEIGH_DATA), intent(in)      :: neighbour_list(:)
    type(EXT_META), intent(inout) :: extmeta

    INTEGER(I4B)                           :: i

    extmeta%l_map(vert) = typ
    forall (i=1:nneigh) extmeta%l_map(neighbour_list(i)%pix) = typ
    return

  END SUBROUTINE markextremum

  SUBROUTINE extremum_properties(vert, dc, nn, nd, a, v, x, ext)
    INTEGER(I8B),    intent(in)    :: vert
    REAL(SP),        intent(in)    :: dc
    INTEGER(I4B),    intent(in)    :: nn(:)
    INTEGER(I4B),    intent(in)    :: nd
    REAL(DP),        intent(in)    :: a(:, :), v(:), x(:)
    TYPE(EXT_DATA),  intent(inout) :: ext

    INTEGER(I4B)      :: WORK(10)
    INTEGER(I4B)      :: INFO

    ext%val       = dc + DOT_PRODUCT(x, 0.5_dp*MATMUL(a, x)+v)
    ext%pos(1:nd) = x + index_to_grid(vert, nn, nd)
    ext%pix       = grid_to_index(nint(fkvadr(ext%pos(1:nd), nn)), nn, nd)

    ! Find eigenvalues
    call DSYEV( 'N', 'L', nd, a, nd, ext%eig, WORK, 10, INFO )

    ! Set type
    ext%typ = nd + 1 - COUNT( ext%eig > 0 )

    return
  END SUBROUTINE extremum_properties

  !> Copies A first nd values into vm, the nd next into am diagonal
  !! and the rest of a into the out-of-diagonal part
  SUBROUTINE setmatrix(a, am, vm, nd)
    real(DP),  intent(in)    :: a(:) !! test
    integer(I4B), intent(in) :: nd
    real(DP),  intent(out)   :: am(nd, nd), vm(nd)

    integer(I4B)             :: i, j, m
    do i=1, nd
       vm(i)=a(i)
       am(i, i)=a(nd+i)
    enddo
    m=2*nd+1
    do j=2, nd
       do i=1, j-1
          am(i, j)=a(m)
          am(j, i)=am(i, j)
          m=m+1
       enddo
    enddo
    return

  END SUBROUTINE setmatrix

END MODULE extrema_mod
