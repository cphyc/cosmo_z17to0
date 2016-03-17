module misc
  implicit none
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
       write(*,*)'Maximum bit length=', bit_size(bit_length)
       write(*,*)'stop in hilbert3d'
       stop
    endif

    state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
         &   0, 1, 3, 2, 7, 6, 4, 5,&
         &   2, 6, 0, 7, 8, 8, 0, 7,&
         &   0, 7, 1, 6, 3, 4, 2, 5,&
         &   0, 9,10, 9, 1, 1,11,11,&
         &   0, 3, 7, 4, 1, 2, 6, 5,&
         &   6, 0, 6,11, 9, 0, 9, 8,&
         &   2, 3, 1, 0, 5, 4, 6, 7,&
         &  11,11, 0, 7, 5, 9, 0, 7,&
         &   4, 3, 5, 2, 7, 0, 6, 1,&
         &   4, 4, 8, 8, 0, 6,10, 6,&
         &   6, 5, 1, 2, 7, 4, 0, 3,&
         &   5, 7, 5, 3, 1, 1,11,11,&
         &   4, 7, 3, 0, 5, 6, 2, 1,&
         &   6, 1, 6,10, 9, 4, 9,10,&
         &   6, 7, 5, 4, 1, 0, 2, 3,&
         &  10, 3, 1, 1,10, 3, 5, 9,&
         &   2, 5, 3, 4, 1, 6, 0, 7,&
         &   4, 4, 8, 8, 2, 7, 2, 3,&
         &   2, 1, 5, 6, 3, 0, 4, 7,&
         &   7, 2,11, 2, 7, 5, 8, 5,&
         &   4, 5, 7, 6, 3, 2, 0, 1,&
         &  10, 3, 2, 6,10, 3, 4, 4,&
         &   6, 1, 7, 0, 5, 2, 4, 3 /), &
         & (/8 ,2, 12 /) )

    do ip = 1,npoint

       ! convert to binary
       do i = 0,bit_length-1
          x_bit_mask(i) = btest(x(ip),i)
          y_bit_mask(i) = btest(y(ip),i)
          z_bit_mask(i) = btest(z(ip),i)
       enddo

       ! interleave bits
       do i = 0,bit_length-1
          i_bit_mask(3*i+2) = x_bit_mask(i)
          i_bit_mask(3*i+1) = y_bit_mask(i)
          i_bit_mask(3*i  ) = z_bit_mask(i)
       end do

       ! build Hilbert ordering using state diagram
       cstate = 0
       do i = bit_length-1,0,-1
          b2 = 0 ; if (i_bit_mask(3*i+2))b2 = 1
          b1 = 0 ; if (i_bit_mask(3*i+1))b1 = 1
          b0 = 0 ; if (i_bit_mask(3*i  ))b0 = 1
          sdigit = b2*4+b1*2+b0
          nstate = state_diagram(sdigit,0,cstate)
          hdigit = state_diagram(sdigit,1,cstate)
          i_bit_mask(3*i+2) = btest(hdigit,2)
          i_bit_mask(3*i+1) = btest(hdigit,1)
          i_bit_mask(3*i  ) = btest(hdigit,0)
          cstate = nstate
       enddo

       ! save Hilbert key as double precision real
       order(ip) = 0.
       do i = 0,3*bit_length-1
          b0 = 0 ; if (i_bit_mask(i))b0 = 1
          order(ip) = order(ip)+dble(b0)*dble(2)**i
       end do

    end do

  end subroutine hilbert3D

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

  subroutine get_cpu_list(X0, X1, levelmax, bound_key, cpu_list, ncpu, ndim)
    real(kind = 8), intent(in)                   :: ncpu, ndim, levelmax
    real(kind = 8), dimension(1:ndim), intent(in):: X0, X1
    real(kind = 8), dimension(0:ncpu), intent(in):: bound_key
    integer, dimension(ncpu), intent(out)        :: cpu_list

    real(kind = 8)                               :: xmin, xmax, ymin, ymax, zmin, zmax
    logical, dimension(ncpu)                     :: cpu_read
    integer                                      :: imin, imax, jmin, jmax, kmin, kmax, lmin, ipart
    integer                                      :: ilevel, bit_length, maxdom
    real(kind = 8), dimension(1:8)               :: bounding_min, bounding_max
    real(kind = 8)                               :: dkey, dmax, deltax
    real(kind = 8), dimension(1:1)               :: order_min
    integer, dimension(1:8)                      :: idom, jdom, kdom, cpu_min, cpu_max

    integer :: ndom, i, impi, ncpu_read, j
    xmin = X0(1); xmax = X1(1)
    ymin = X0(2); ymax = X1(2)
    zmin = X0(3); zmax = X1(3)

    dmax = max(xmax-xmin, ymax-ymin, zmax-zmin)
    do ilevel = 1, levelmax
       deltax = 0.5d0**ilevel
       if (deltax.lt.dmax)exit
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
    if (bit_length>0)ndom = 8
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
    ncpu_read = 0

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
  end subroutine get_cpu_list
  
end module cpus

module io
  implicit none
  real(kind = 8), dimension(:), allocatable :: bound_key
  integer                                   :: ncpu, ndim, levelmin, levelmax
contains
  subroutine read_info(filename)
    character(len=*), intent(in)                           :: filename

    logical                                                :: ok
    integer                                                :: impi, i
    real(kind = 8)                                         :: t, aexp, unit_l, unit_t
    character(len=80)                                      :: ordering

    inquire(file=filename, exist=ok)
    if (.not. ok) then
       print*, filename // ' not found'
       stop
    end if

    open(unit=10, file=filename, form='formatted', status='old')
    read(10,'("ncpu        =",I11)') ncpu
    read(10,'("ndim        =",I11)') ndim
    read(10,'("levelmin    =",I11)') levelmin
    read(10,'("levelmax    =",I11)') levelmax
    read(10,*)
    read(10,*)
    read(10,*)

    read(10,*)
    read(10,'("time        =",E23.15)') t
    read(10,'("aexp        =",E23.15)') aexp
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,'("unit_l      =",E23.15)') unit_l
    read(10,*)
    read(10,'("unit_t      =",E23.15)') unit_t

    read(10,*)
    read(10,'("ordering type=",A80)') ordering
    read(10,*)

    if (TRIM(ordering) == 'hilbert') then
       allocate(bound_key(0:ncpu))
       do impi = 1,ncpu
          read(10,'(I8,1X,E23.15,1X,E23.15)') i, bound_key(impi-1), bound_key(impi)
       end do
    endif
    close(10)
  end subroutine read_info

  subroutine free ()
    if (allocated(bound_key)) then
       deallocate(bound_key)
    end if
  end subroutine free

end module io
