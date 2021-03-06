module misc
  use types
  use omp_lib

  implicit none

  interface minmax
     subroutine minmax_r(array, min, max)
       real(kind=8), intent(in), dimension(:) :: array
       real(kind=8), intent(out) :: min, max
     end subroutine minmax_r

     subroutine minmax_i (array, min, max)
       integer, intent(in), dimension(:) :: array

       integer, intent(out) :: min, max
     end subroutine minmax_i
  end interface minmax

  interface meanval
     subroutine mean_1(array, mean, n)
       integer, intent(in) :: n
       real(kind=8), intent(in), dimension(n) :: array
       real(kind=8), intent(out) :: mean
     end subroutine mean_1
     subroutine mean_2(array, mean, n, m)
       integer, intent(in) :: n, m
       real(kind=8), intent(in), dimension(n, m) :: array
       real(kind=8), intent(out), dimension(n)   :: mean
     end subroutine mean_2
  end interface meanval

  interface stddev
     subroutine stddev_1(array, std, n)
       integer, intent(in) :: n
       real(kind=8), intent(in), dimension(n) :: array
       real(kind=8), intent(out) :: std
     end subroutine stddev_1
     subroutine stddev_2(array, std, n, m)
       integer, intent(in) :: n, m
       real(kind=8), intent(in), dimension(n,m) :: array
       real(kind=8), intent(out), dimension(n)  :: std
     end subroutine stddev_2
  end interface stddev
contains
  subroutine hilbert3D(x, y, z, order, bit_length, npoint)
    implicit none

    integer       , INTENT(IN)                       :: bit_length, npoint
    integer       , INTENT(IN) , dimension(1:npoint) :: x, y, z
    real(kind = 8), INTENT(OUT), dimension(1:npoint) :: order

    logical, dimension(0:3*bit_length-1)           :: i_bit_mask
    logical, dimension(0:1*bit_length-1)           :: x_bit_mask, y_bit_mask, z_bit_mask
    integer, dimension(0:7, 0:1, 0:11)             :: state_diagram
    integer                                        :: i, ip, cstate, nstate, b0, b1, b2
    integer                                        :: sdigit, hdigit

    if (bit_length>bit_size(bit_length)) then
       write(*, *)'Maximum bit length=', bit_size(bit_length)
       write(*, *)'stop in hilbert3d'
       stop
    endif

    state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5, &
         &   0, 1, 3, 2, 7, 6, 4, 5, &
         &   2, 6, 0, 7, 8, 8, 0, 7, &
         &   0, 7, 1, 6, 3, 4, 2, 5, &
         &   0, 9,10, 9, 1, 1,11,11, &
         &   0, 3, 7, 4, 1, 2, 6, 5, &
         &   6, 0, 6,11, 9, 0, 9, 8, &
         &   2, 3, 1, 0, 5, 4, 6, 7, &
         &  11,11, 0, 7, 5, 9, 0, 7, &
         &   4, 3, 5, 2, 7, 0, 6, 1, &
         &   4, 4, 8, 8, 0, 6,10, 6, &
         &   6, 5, 1, 2, 7, 4, 0, 3, &
         &   5, 7, 5, 3, 1, 1,11,11, &
         &   4, 7, 3, 0, 5, 6, 2, 1, &
         &   6, 1, 6,10, 9, 4, 9,10, &
         &   6, 7, 5, 4, 1, 0, 2, 3, &
         &  10, 3, 1, 1,10, 3, 5, 9, &
         &   2, 5, 3, 4, 1, 6, 0, 7, &
         &   4, 4, 8, 8, 2, 7, 2, 3, &
         &   2, 1, 5, 6, 3, 0, 4, 7, &
         &   7, 2,11, 2, 7, 5, 8, 5, &
         &   4, 5, 7, 6, 3, 2, 0, 1, &
         &  10, 3, 2, 6,10, 3, 4, 4, &
         &   6, 1, 7, 0, 5, 2, 4, 3 /), &
         & (/8 ,2, 12 /) )

    do ip = 1, npoint

       ! convert to binary
       do i = 0, bit_length-1
          x_bit_mask(i) = btest(x(ip), i)
          y_bit_mask(i) = btest(y(ip), i)
          z_bit_mask(i) = btest(z(ip), i)
       enddo

       ! interleave bits
       do i = 0, bit_length-1
          i_bit_mask(3*i+2) = x_bit_mask(i)
          i_bit_mask(3*i+1) = y_bit_mask(i)
          i_bit_mask(3*i  ) = z_bit_mask(i)
       end do

       ! build Hilbert ordering using state diagram
       cstate = 0
       do i = bit_length-1, 0,-1
          b2 = 0 ; if (i_bit_mask(3*i+2)) b2 = 1
          b1 = 0 ; if (i_bit_mask(3*i+1)) b1 = 1
          b0 = 0 ; if (i_bit_mask(3*i  )) b0 = 1
          sdigit = b2*4+b1*2+b0
          nstate = state_diagram(sdigit, 0,cstate)
          hdigit = state_diagram(sdigit, 1,cstate)
          i_bit_mask(3*i+2) = btest(hdigit, 2)
          i_bit_mask(3*i+1) = btest(hdigit, 1)
          i_bit_mask(3*i  ) = btest(hdigit, 0)
          cstate = nstate
       enddo

       ! save Hilbert key as double precision real
       order(ip) = 0.
       do i = 0, 3*bit_length-1
          b0 = 0 ; if (i_bit_mask(i)) b0 = 1
          order(ip) = order(ip)+dble(b0)*dble(2)**i
       end do

    end do

  end subroutine hilbert3D

  subroutine quick_sort(list, order)
    ! quick sort routine from:
    ! brainerd, w.s., goldberg, c.h. & adams, j.c. (1990) "programmer's guide to
    ! fortran 90", mcgraw-hill  isbn 0-07-000248-7, pages 149-150.
    ! modified by alan miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    integer, parameter :: i8b = 8

    integer, dimension (:), intent(inout)        :: list
    integer, dimension (size(list)), intent(out) :: order

    integer :: n
    ! local variable
    integer :: i

    n = size(list)
    do i = 1, n
       order(i) = i
    end do

    call quick_sort_1(1, n)

  contains
    recursive subroutine quick_sort_1(left_end, right_end)

      integer, intent(in) :: left_end, right_end

      !     local variables
      integer             :: i, j, itemp
      integer(i8b)        :: reference, temp
      integer, parameter  :: max_simple_sort_size = 6

      if (right_end < left_end + max_simple_sort_size) then
         ! use interchange sort for small lists
         call interchange_sort(left_end, right_end)

      else
         ! use partition ("quick") sort
         reference = list((left_end + right_end)/2)
         i = left_end - 1; j = right_end + 1

         do
            ! scan list from left end until element >= reference is found
            do
               i = i + 1
               if (list(i) >= reference) exit
            end do
            ! scan list from right end until element <= reference is found
            do
               j = j - 1
               if (list(j) <= reference) exit
            end do


            if (i < j) then
               ! swap two out-of-order elements
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            else if (i == j) then
               i = i + 1
               exit
            else
               exit
            end if
         end do

         if (left_end < j) call quick_sort_1(left_end, j)
         if (i < right_end) call quick_sort_1(i, right_end)
      end if

    end subroutine quick_sort_1

    subroutine interchange_sort(left_end, right_end)

      integer, intent(in) :: left_end, right_end

      !     local variables
      integer             :: i, j, itemp
      integer(i8b)        :: temp

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (list(i) > list(j)) then
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort
  end subroutine quick_sort

  recursive subroutine rquick_sort(A)
    real, intent(in out), dimension(:) :: A
    integer :: iq

    if(size(A) > 1) then
       call Partition(A, iq)
       call rquick_sort(A(:iq-1))
       call rquick_sort(A(iq:))
    endif

  contains
    subroutine Partition(A, marker)
      real, intent(in out), dimension(:) :: A
      integer, intent(out) :: marker
      integer :: i, j
      real :: temp
      real :: x      ! pivot point
      x = A(1)
      i= 0
      j= size(A) + 1

      do
         j = j-1
         do
            if (A(j) <= x) exit
            j = j-1
         end do
         i = i+1
         do
            if (A(i) >= x) exit
            i = i+1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
         elseif (i == j) then
            marker = i+1
            return
         else
            marker = i
            return
         endif
      end do

    end subroutine Partition

  end subroutine rquick_sort

  subroutine get_cpu_list(X0, X1, levelmax, bound_key, cpu_list, ncpu, ndim)
    integer, intent(in)                          :: ncpu, ndim, levelmax
    real(kind = 8), dimension(1:ndim), intent(in):: X0, X1
    real(kind = 8), dimension(0:ncpu), intent(in):: bound_key
    integer, dimension(ncpu), intent(out)        :: cpu_list

    logical, dimension(ncpu)                     :: cpu_read
    integer :: ncpu_read

    real(kind = 8), dimension(1:ndim) :: X0_tmp, X1_tmp
    real(kind=8), dimension(3, 2, 2)  :: bounds
    integer :: i, j, k

    cpu_list = 0
    ncpu_read = 0
    bounds = 0
    cpu_read = .false.

    do i = 1, ndim
       if (X0(i) < 0 .and. X1(i) < 0) then
          bounds(i, 1, 1) = X0(i) + 1d0
          bounds(i, 1, 2) = X1(i) + 1d0

          bounds(i, 2, 1) = -1.
          bounds(i, 2, 2) = -1.
       else if (X0(i) < 0) then
          bounds(i, 1, 1) = X0(i) + 1d0
          bounds(i, 1, 2) = 1d0

          bounds(i, 2, 1) = 0d0
          bounds(i, 2, 2) = X1(i)
       else if (X0(i) > 1. .and. X1(i) > 1.) then
          bounds(i, 1, 1) = X0(i) - 1d0
          bounds(i, 1, 2) = X1(i) - 1d0
       else if (X1(i) > 1.) then
          bounds(i, 1, 1) = X0(i)
          bounds(i, 1, 2) = 1d0

          bounds(i, 2, 1) = 0d0
          bounds(i, 2, 2) = X1(i) - 1d0
       else
          bounds(i, 1, 1) = X0(i)
          bounds(i, 1, 2) = X1(i)

          bounds(i, 2, 1) = -1.
          bounds(i, 2, 2) = -1.
       end if
    end do

    do i = 1, 2
       do j = 1, 2
          do k = 1, 2
             X0_tmp(1) = bounds(1, i, 1)
             X1_tmp(1) = bounds(1, i, 2)
             X0_tmp(2) = bounds(2, j, 1)
             X1_tmp(2) = bounds(2, j, 2)
             X0_tmp(3) = bounds(3, k, 1)
             X1_tmp(3) = bounds(3, k, 2)
             ! print*, '------------------------'
             ! print*, X0_tmp, X1_tmp
             call get_it(cpu_list, ncpu_read, cpu_read)
          end do
       end do
    end do

  contains

    subroutine get_it(cpu_list, ncpu_read, cpu_read)
      integer, dimension(ncpu), intent(inout)      :: cpu_list
      logical, dimension(ncpu), intent(inout)      :: cpu_read
      integer, intent(inout)                       :: ncpu_read

      real(kind = 8)                 :: xmin, xmax, ymin, ymax, zmin, zmax
      integer                        :: imin, imax, jmin, jmax, kmin, kmax, lmin
      integer                        :: ilevel, bit_length, maxdom
      real(kind = 8), dimension(1:8) :: bounding_min, bounding_max
      real(kind = 8)                 :: dkey, dmax, deltax
      real(kind = 8), dimension(1:1) :: order_min
      integer, dimension(1:8)        :: idom, jdom, kdom, cpu_min, cpu_max

      integer :: ndom, i, impi, j
      xmin = X0_tmp(1); xmax = X1_tmp(1)
      ymin = X0_tmp(2); ymax = X1_tmp(2)
      zmin = X0_tmp(3); zmax = X1_tmp(3)

      dmax = max(xmax-xmin, ymax-ymin, zmax-zmin)
      do ilevel = 1, levelmax
         deltax = 0.5d0**ilevel
         if (deltax.lt.dmax) exit
      end do
      lmin = ilevel
      bit_length = lmin-1
      maxdom = 2**bit_length
      imin = 0; imax = 0; jmin = 0; jmax = 0; kmin = 0; kmax = 0
      if (bit_length>0) then
         imin = int(xmin*dble(maxdom))
         imax = imin+1
         jmin = int(ymin*dble(maxdom))
         jmax = jmin+1
         kmin = int(zmin*dble(maxdom))
         kmax = kmin+1
      endif

      dkey = (dble(2**(levelmax+1)/dble(maxdom)))**ndim
      ndom = 1
      if (bit_length>0) ndom = 8
      idom(1) = imin; idom(2) = imax
      idom(3) = imin; idom(4) = imax
      idom(5) = imin; idom(6) = imax
      idom(7) = imin; idom(8) = imax
      jdom(1) = jmin; jdom(2) = jmin
      jdom(3) = jmax; jdom(4) = jmax
      jdom(5) = jmin; jdom(6) = jmin
      jdom(7) = jmax; jdom(8) = jmax
      kdom(1) = kmin; kdom(2) = kmin
      kdom(3) = kmin; kdom(4) = kmin
      kdom(5) = kmax; kdom(6) = kmax
      kdom(7) = kmax; kdom(8) = kmax

      do i = 1, ndom
         if (bit_length>0) then
            call hilbert3D(idom(i), jdom(i), kdom(i), order_min, bit_length, 1)
         else
            order_min = 0.0d0
         endif
         bounding_min(i) = (order_min(1))*dkey
         bounding_max(i) = (order_min(1)+1.0D0)*dkey
      end do

      cpu_min = 0; cpu_max = 0
      do impi = 1, ncpu
         do i = 1, ndom
            if (   bound_key(impi-1) <= bounding_min(i).and.&
                 & bound_key(impi  ) > bounding_min(i)) then
               cpu_min(i) = impi
            endif
            if (   bound_key(impi-1) < bounding_max(i).and.&
                 & bound_key(impi  ) >= bounding_max(i)) then
               cpu_max(i) = impi
            endif
         end do
      end do

      do i = 1, ndom
         do j = cpu_min(i), cpu_max(i)
            if (.not. cpu_read(j)) then
               ncpu_read = ncpu_read+1
               cpu_list(ncpu_read) = j
               cpu_read(j) = .true.
            endif
         enddo
      enddo

      ! deallocate(cpu_read)
    end subroutine get_it
  end subroutine get_cpu_list

  function indexOf(element, array)
    ! Please note that array has to be sorted
    integer :: indexOf
    integer, intent(in) :: element
    integer, intent(in), dimension(:) :: array

    integer :: left, right, middle
    left = lbound(array, 1)
    right = ubound(array, 1)
    middle = (left + right) / 2

    indexOf = -1

    if (array(left) > element .or. array(right) < element) then
       return
    end if
    do while (element /= array(middle))
       if (left > right) then
          return
       end if

       if (element > array(middle)) then
          left = middle + 1
       else if (element < array(middle)) then
          right = middle - 1
       end if
       middle = (left + right) / 2

    end do

    indexOf = middle

  end function indexOf

  subroutine max_index (array, imax)
    real(sp), intent(in), dimension(:) :: array
    integer, intent(out) :: imax
    real :: maxval
    integer :: i
    imax = lbound(array, 1)
    maxval = array(imax)
    do i = lbound(array, 1), ubound(array, 1)
       if (array(i) > maxval) then
          imax = i
          maxval = array(i)
       end if
    end do
  end subroutine max_index

  function median (A, n, m)
    integer, intent(in) :: n, m
    real, dimension(n, m), intent(in) :: A

    real, dimension(m) :: A_copy

    integer :: i
    real, dimension(n) :: median

    do i = 1, n
       A_copy = A(i, :)
       call rquick_sort(A_copy)

       median(i) = A_copy(m/ 2)
    end do
  end function median

  ! Return an array containing unique data from array
  subroutine unique (array, u_array)
    integer, intent(in), dimension(:) :: array
    integer, intent(out), dimension(size(array)) :: u_array

    integer, dimension(size(array)) :: order, tmp_array
    integer :: ptr, i

    tmp_array = array
    call quick_sort(tmp_array, order)

    u_array = 0
    u_array(1) = tmp_array(1)

    ptr = 1 ! pointer on unique array
    do i = 1, size(array)
       if (u_array(ptr) /= tmp_array(i)) then
          ptr = ptr + 1
          u_array(ptr) = tmp_array(i)
       end if
    end do

  end subroutine unique

  subroutine parse_params (cli)
    use flap, only : command_line_interface

    type(command_line_interface), intent(out) :: cli

    call cli%init(progname='compute_halo_fft')
    call cli%add(switch='--min-mass', switch_ab='-minm', help='Minimum mass', act='store',&
         def='0')
    call cli%add(switch='--max-mass', switch_ab='-maxm', help='Maximum mass', act='store',&
         def='1e99')
    call cli%add(switch='--gal-list', help='List of galaxies', act='store', &
         def='lists/list_kingal_00782.dat')
    call cli%add(switch='--halo-list', help='List of dark matter halo', act='store', &
         def='lists/list_halo.dat.bin')
    call cli%add(switch='--association-list', help='List of association between galaxy and halo', &
         act='store', def='lists/associated_halogal_782.dat.bin')
    call cli%add(switch='--info-file', help='Information file about simulation', &
         act='store', def='/data52/Horizon-AGN/OUTPUT_DIR/output_00782/info_00782.txt')
    call cli%add(switch='--output-end', help='Path for latest output of simulation', &
         act='store', def='/data52/Horizon-AGN/OUTPUT_DIR/output_00782/')
    call cli%add(switch='--output-start', help='Path for first output of simulation', &
         act='store', def='/data52/Horizon-AGN/OUTPUT_DIR/output_00002/')
    call cli%add(switch='--brick', help='Path for brick file that contains the halos', &
         act='store', def='/data52/Horizon-AGN/TREE_DM_celldx2kpc_SC0.9r/tree_bricks782')
    call cli%add(switch='--cpu-from', help='First cpu to use', &
         act='store', def='1')
    call cli%add(switch='--cpu-to', help='Last cpu to use', &
         act='store', def='4096')
    call cli%add(switch='--halo-to-cpu', help='Halo to cpu list, binary', &
         act='store', def='lists/halo_to_cpu.00002.m<1e12.dat.bin')
    call cli%add(switch='--percent-probe-particle', help='When using random particles, percent of random particle to pick.', &
         act='store', def='5')
    call cli%add(switch='--output-path', help='Path of the output', &
         act='store', def='/data52/Horizon-AGN/OUTPUT_DIR')
    call cli%add(switch='--output-number', help='Number of the output', &
         act='store', def='2', nargs='+')
    call cli%add(switch='--output', switch_ab='-o', help='Name of the output file', &
         act='store', def='data.out')
    call cli%add(switch='--verbose', help='Verbosity', &
         act='store', def='0')
    call cli%add(switch='--halo-i', help='Halo number', &
         act='store', def='-1')
    call cli%add(switch='--halo-to-treat', help='Halo to treat', &
         act='store', nargs='*', def='-1')


  end subroutine parse_params

  ! Fill an array using val. If the value if found in array, or the array
  ! has already been filled do nothing, else add the value and increment the
  ! counter
  subroutine fill(array, val, counter)

    integer, intent(in) :: val
    integer, intent(inout), dimension(:) :: array

    integer, intent(out) :: counter
    integer :: i

    do i = 1, size(array)
       if (array(i) == val) then
          exit
       else if (array(i)  == 0) then
          counter = counter + 1
          array(i) = val
          exit
       end if
    end do

  end subroutine fill

  !> Filter the region around center, with size width
  !! inputs:
  !! -----------
  !!   - center: the center, in ramses units
  !!   - width: the width, in ramses units
  !!   - idsins: the ids of the particles
  !!   - in : the position of the particles
  !!   - order, optional: the order of the ids
  !! outputs:
  !! --------
  !!   - idsout: the ids of the particles found
  !!   - out: the position of the particles
  !!   - parts_in_region: number of particles in the region
  subroutine filter_region(center, width, idsin, in, idsout, out, parts_in_region, &
       order)
    real(dp), dimension(3), intent(in)                         :: center, width
    real(dp), dimension(:, :), intent(in)                      :: in
    integer, dimension(:), intent(in)                          :: idsin
    integer, dimension(:), intent(in), optional                :: order ! if the idsin have been sorted
    real(dp), dimension(size(in, 1), size(in, 2)), intent(out) :: out
    integer, intent(out), optional                             :: parts_in_region
    integer, dimension(:), intent(out)                         :: idsout

    integer :: i, dim, iout, tmp_iout

    real(dp), dimension(3) :: maxis, minis, dist, correction
    logical :: ok

    real(dp) :: start
    integer :: nthreads, thread, from, to, step, nfound_shared, nfound, part_i
    ! integer, allocatable, dimension(:) :: tmp_ids

    iout = 0

    start = OMP_GET_WTIME()
    nthreads = OMP_GET_MAX_THREADS()
    step = ceiling(1d0 * size(idsin, 1) / nthreads)

    nfound_shared = 0

    !$OMP PARALLEL DEFAULT(none)                                       &
    !$OMP PRIVATE(thread, from, to, ok, dist, correction, dim, part_i) &
    !$OMP PRIVATE(nfound)                                              &
    !$OMP FIRSTPRIVATE(nthreads, step, center, width)                  &
    !$OMP SHARED(idsin, nfound_shared, order, in, idsout, out)

    thread = OMP_GET_THREAD_NUM()
    from = step * thread + 1
    to   = min(step * (thread + 1), size(idsin, 1))
    nfound = 0

    !  allocate(tmp_ids(to-from+1))

    print*, width

    do part_i = from, to
       dist = in(:, part_i) - center

       ! correct the distance with boundary conditions (cannot exceed 0.5)
       do dim = 1, 3
          if (dist(dim) > 0.5) then
             dist(dim) = abs(dist(dim) - 1d0)
             correction(dim) = -1d0
          else if (dist(dim) < -0.5) then
             dist(dim) = abs(dist(dim) + 1d0)
             correction(dim) = 1d0
          else
             dist(dim) = abs(dist(dim))
             correction(dim) = 0d0
          end if
       end do

       ok = all(dist < width)

       if (ok) then
          !nfound = nfound + 1
          !tmp_ids(nfound) = part_i ! tmp_ids(i) contains the index

          ! TODO: actually apply the correction...
          !$OMP CRITICAL
          nfound_shared = nfound_shared + 1
          if (present(order)) then
             idsout(nfound_shared) = idsin(order(part_i))
          else
             idsout(nfound_shared) = idsin(part_i)
          end if

          out(:, nfound_shared) = in(:, part_i) + correction(:)
          !$OMP END CRITICAL
       end if
    end do
    !$OMP END PARALLEL

    ! copy the result in the shared array
    ! !$OMP CRITICAL
    ! do i = 1, nfound
    !    nfound_shared = nfound_shared + 1

    !    part_i = tmp_ids(i)
    !    ! write(*, '(i4,*(a,i9))') thread, ': part_i ', part_i, 'i', i

    !    if (present(order)) then
    !       idsout(nfound_shared) = idsin(order(part_i))
    !    else
    !       idsout(nfound_shared) = idsin(part_i)
    !    end if

    !    out(:, nfound_shared) = in(:, part_i)
    ! end do
    ! !$OMP END CRITICAL
    ! deallocate(tmp_ids)

    write(*, *) 'Took', OMP_GET_WTIME() - start

    parts_in_region = nfound_shared
  end subroutine filter_region

   !> Filter the region around center, with size width
  !! inputs:
  !! -----------
  !!   - center: the center, in ramses units
  !!   - width: the width, in ramses units
  !!   - idsins: the ids of the particles
  !!   - in : the position of the particles
  !!   - order, optional: the order of the ids
  !! outputs:
  !! --------
  !!   - idsout: the ids of the particles found
  !!   - out: the position of the particles
  !!   - parts_in_region: number of particles in the region
  subroutine filter_region_no_omp(center, width, idsin, in, idsout, out, parts_in_region, &
       order)
    real(dp), dimension(3), intent(in)                         :: center, width
    real(dp), dimension(:, :), intent(in)                      :: in
    integer, dimension(:), intent(in)                          :: idsin
    integer, dimension(:), intent(in), optional                :: order ! if the idsin have been sorted
    real(dp), dimension(size(in, 1), size(in, 2)), intent(out) :: out
    integer, intent(out), optional                             :: parts_in_region
    integer, dimension(:), intent(out)                         :: idsout

    integer :: i, dim, iout

    real(dp), dimension(3) :: maxis, minis, dist, correction
    logical :: ok

    iout = 0
    ! loop over each positions given in input

    !$nOMP PARALLEL DO default(none)                   &
    !$nOMP SHARED(idsin, idsout, out, order, iout, in) &
    !$nOMP FIRSTPRIVATE(center, width)                 &
    !$nOMP PRIVATE(ok, dist, correction)               &
    !$nOMP SCHEDULE(DYNAMIC, 1000000)
    do i = 1, size(idsin, 1)
       ok = .true.
       dist = in(:, i) - center
       correction = 0

       ! correct the distance with boundary conditions (cannot exceed 0.5)
       do dim = 1, 3
          if (dist(dim) > 0.5) then
             dist(dim) = dist(dim) - 1
             correction(dim) = -1
          else if (dist(dim) < -0.5) then
             dist(dim) = dist(dim) + 1
             correction(dim) = 1
          end if

          dist(dim) = abs(dist(dim))

          if (dist(dim) > width(dim)) then
             ok = .false.
             exit
          end if
       end do

       if (ok) then
          !$OMP CRITICAL

          iout = iout + 1
          if (present(order)) then
             ! get the real id using order array
             idsout(iout) = idsin(order(i))
          else
             idsout(iout) = idsin(i)
          end if

          out(:, iout) = in(:, i) + correction(:)
          !$OMP END CRITICAL
       end if
    end do

    parts_in_region = iout
  end subroutine filter_region_no_omp

  !! Finds the indexes of B in A when both arrays are sorted
  !! The subroutine finds the position of B's leftmost and rightmost values
  !! into A. All the remaining values are within left[B], right[B]
  !! and therefore within A[index(left[B]), index(right[B])]
  !! this is then repeated until all elements of B are exhausted
  subroutine indexOfArr(B, A, indexArr)
    integer, intent(in), dimension(:) :: A, B
    integer, intent(out), dimension(size(B))&
                                      :: indexArr

    integer :: left, right, i, j
    left  = 1
    right = size(A)

    i = 1
    j = size(B)

    do while (i <= j)
       ! Look left one in subset of A
       indexArr(i) = indexOf(B(i), A(left:right))
       if (indexArr(i) > 0) then
          indexArr(i) = indexArr(i) + left - 1
          left = indexArr(i)
       end if

       if (j == i) exit
       i = i + 1

       ! Look right one in subset of A
       indexArr(j) = indexOf(B(j), A(left:right))
       if (indexArr(j) > 0) then
          indexArr(j) = indexArr(j) + left - 1
          right = indexArr(j)
       end if

       j = j - 1
    end do


  end subroutine indexOfArr


end module misc
!-------------------------------------
! To interface
!-------------------------------------
subroutine minmax_r (array, min, max)
  implicit none

  real(kind=8), intent(in), dimension(:) :: array

  real(kind=8), intent(out) :: min, max
  integer :: i

  do i = 1, size(array)
     if (array(i) > min) min = array(i)
     if (array(i) < max) max = array(i)
  end do
end subroutine minmax_r

subroutine minmax_i (array, min, max)
  implicit none

  integer, intent(in), dimension(:) :: array

  integer, intent(out) :: min, max
  integer :: i

  do i = 1, size(array)
     if (array(i) > min) min = array(i)
     if (array(i) < max) max = array(i)
  end do
end subroutine minmax_i


subroutine mean_1 (array, mean, n)
  integer, intent(in) :: n
  real(kind=8), intent(in), dimension(n) :: array
  real(kind=8), intent(out)              :: mean

  mean = sum(array) / size(array)
end subroutine mean_1

subroutine mean_2 (array, mean, n, m)
  integer, intent(in) :: n, m
  real(kind=8), intent(in), dimension(n, m) :: array
  real(kind=8), intent(out), dimension(n)   :: mean
  real(kind=8), dimension(m)   :: tmp_arr
  integer :: i

  do i = 1, n
     tmp_arr = array(i, :)
     call mean_1(tmp_arr, mean(i), m)
  end do
end subroutine mean_2

subroutine stddev_1 (array, std, n)
  integer, intent(in) :: n
  real(kind=8), intent(in), dimension(n) :: array
  real(kind=8), intent(out)              :: std

  real(kind=8), dimension(n) :: tmp_array
  real(kind=8) :: mean

  call mean_1(array, mean, n)
  tmp_array = (array - mean)**2
  call mean_1(tmp_array, std, n)
  std = sqrt(std)
end subroutine stddev_1

subroutine stddev_2 (array, std, n, m)
  integer, intent(in) :: n, m
  real(kind=8), intent(in), dimension(n, m) :: array
  real(kind=8), intent(out), dimension(n)   :: std
  real(kind=8), dimension(m)   :: tmp_arr
  integer :: i

  do i = 1, n
     tmp_arr = array(i, :)
     call stddev_1(tmp_arr, std(i), m)
  end do
end subroutine stddev_2
