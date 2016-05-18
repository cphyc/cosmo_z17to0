!> Get the tree for a given halo
program get_tree
  use readtree
  use writetree
  use misc
  use flap, only : command_line_interface
  implicit none

  !----------------------------------------
  ! cli params
  !----------------------------------------
  type(command_line_interface) :: cli
  integer, allocatable, dimension(:) :: param_halos
  character(len=128) :: param_output
  integer :: param_verbosity

  !----------------------------------------
  ! loops
  !----------------------------------------
  integer :: i, j
  !----------------------------------------
  ! Tree paramters
  !----------------------------------------
  logical          :: stop_now
  integer(kind=4)  :: halo_id, BushID, step, level, hosthalo, hostsub
  integer(kind=4)  :: nbsub, nextsub, nb_of_fathers, nb_of_sons
  real(kind=4)     :: m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz
  real(kind=4)     :: r, ra, rb, rc, ek, ep, et, spin

  integer(kind=4) :: prevStep, index

  real(kind=4), allocatable, dimension(:)    :: m_fathers
  integer(kind=4), allocatable, dimension(:) :: id_fathers, id_sons

  ! The list of halos to find
  type halo_list
     integer(kind=4), dimension(10000) :: list = 0
  end type halo_list
  type(halo_list), dimension(:), allocatable :: halos ! FIXME: be able to have more
  integer(kind=4), dimension(10000) :: order
  !----------------------------------------
  ! Get the parameters
  !---------------------------------------- 
  call cli%init(progname='get_tree_of_halo')
  call cli%add(switch='--halo-to-treat', help='Halo to treat', &
       act='store', nargs='*', def='-1')
  call cli%add(switch='--output', switch_ab='-o', help='Name of the output file', &
       act='store', def='data.out')
  call cli%add(switch='--verbose', help='Verbosity', &
       act='store', def='0')

  call cli%get_varying(switch='--halo-to-treat', val=param_halos)
  call cli%get(switch='--output', val=param_output)
  call cli%get(switch='--verbose', val=param_verbosity)

  !----------------------------------------
  ! Read the reversed tree
  !----------------------------------------
  call read_tree('/home/cadiou/data/H-AGN/tree.reversed.dat')

  ! allocate data
  allocate(halos(nsteps))
  do i = 1, nsteps-1
     halos(i)%list = 0
  end do
  call append(halos(nsteps), param_halos)

  if (param_verbosity >= 1) then
     write(*, *) 'Looking for halos', param_halos
  end if

  ! Opening output file
  if (param_verbosity >= 1) then
     write(*, *) 'Writing output in ', trim(param_output)
  end if
  call write_tree(param_output, nsteps, nb_of_halos, nb_of_subhalos, aexp, omega_t, age_univ)

  ! iter over all halos
  stop_now = .false.
  prevStep = 0
  do while (.not. stop_now)
     ! print*, halo_id, step
     call iter_tree(stop_now, halo_id, BushID, step, level, hosthalo, hostsub, &
          nbsub, nextsub, nb_of_fathers, nb_of_sons, &
          m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz, &
          r, ra, rb, rc, ek, ep, et, spin)

     ! at a new step, reorder the halos to find
     if (step - prevStep /= 0) then
        prevStep = step
        ! sort the halos at the step
        call quick_sort(halos(step)%list, order)

        if (param_verbosity >= 2) then
           write(*, *)    'step :', step
           if (param_verbosity >= 3) then
              write(*, *) 'halos:', halos(step)%list(:sizeOf(halos(step)))
           end if
        end if

     end if

     ! Try to find the current halo in the halo list
     index = indexOf(halo_id, halos(step)%list)
     if (index > 0) then
        allocate(id_fathers(nb_of_fathers), m_fathers(nb_of_fathers))
        allocate(id_sons(nb_of_sons))

        ! add the parents to the next halos to find
        if (nb_of_fathers > 0) then
           if (step-1 >= 1 ) then
              call get_fathers(nb_of_fathers, id_fathers, m_fathers)

              call get_sons(nb_of_sons, id_sons)
              ! write(*, "(1x, a, I8)")    'id      :', halo_id
              ! write(*, "(1x, a, 1000i8)") 'fathers :', id_fathers

              call append(halos(step-1), id_fathers)
              ! write(*, *)
           end if
        end if

        !----------------------------------------
        ! Write the output
        !----------------------------------------
        call write_iter_tree(.false., halo_id, BushID, step, level, hosthalo, hostsub, &
             nbsub, nextsub, nb_of_fathers, nb_of_sons, &
             m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz, &
             r, ra, rb, rc, ek, ep, et, spin, id_fathers, m_fathers, id_sons)

        deallocate(id_fathers, m_fathers, id_sons)
     end if
  end do

  ! Write stop as halo with id == -1
  call write_iter_tree(.true., halo_id, BushID, step, level, hosthalo, hostsub, &
                   nbsub, nextsub, nb_of_fathers, nb_of_sons, &
                   m, macc, px, py, pz, vx, vy, vz, Lx, Ly, Lz, &
                   r, ra, rb, rc, ek, ep, et, spin)

contains
  subroutine append (self, elements)
    type(halo_list), intent(inout) :: self
    integer(kind=4), intent(in), dimension(:)  :: elements
    integer(kind=4), dimension(size(elements)) :: tmp_elements, u_tmp_elements, order_els
    integer(kind=4), dimension(size(self%list)) :: tmp_list, order_list

    integer :: i, from, to, j

    u_tmp_elements = 0
    tmp_elements = elements
    tmp_list     = self%list
    call quick_sort(tmp_elements, order_els)
    call quick_sort(tmp_list, order_list)

    !----------------------------------------
    ! Remove elements already in list
    !----------------------------------------
    j = 0
    do i = 1, size(tmp_elements)
       if (indexOf(tmp_elements(i), tmp_list) <= 0) then
          j = j + 1
          u_tmp_elements(j) = tmp_elements(i)
       else
          ! print*, tmp_elements(i), 'already in list'
       end if
    end do

    ! Reorder array
    call quick_sort(u_tmp_elements, order_els)
    !----------------------------------------
    ! Find first non null index
    !----------------------------------------
    from = 0
    to = size(u_tmp_elements)
    ! Remove 0s from tmp_elements
    do i = 1, to
       if ((u_tmp_elements(i) > 0) .and. from == 0) then
          from = i
       end if
    end do

    ! find first null index in self
    i = sizeOf(self)+1

    ! TODO: check no overflow
    self%list(i:i+to-from) = u_tmp_elements(from:to)
    ! print*, 'WE:', self%list(1:i+to-from)

  end subroutine append

  !> Return the index of the last non null element
  function sizeOf(self)
    type(halo_list), intent(in) :: self
    integer :: sizeOf

    integer :: j

    do j = 1, size(self%list)
       if (self%list(j) == 0) exit
    end do
    sizeOf = j-1
  end function sizeOf

end program get_tree
