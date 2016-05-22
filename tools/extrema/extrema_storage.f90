module extrema_storage
  use extrema_types

  type(EXT_DATA), allocatable :: stored_ext(:, :)

contains

  subroutine extrema_push(newext)
    type(EXT_DATA), intent(in) :: newext(:, :)

    if (allocated(stored_ext)) deallocate(stored_ext)
    allocate(stored_ext(size(newext, 1), size(newext, 2)))

    stored_ext = newext

    print*, 'Extrema has size', size(stored_ext, 1), size(stored_ext, 2)

  end subroutine extrema_push

  subroutine extrema_pull(newext, index)
    type(EXT_DATA), intent(out), allocatable :: newext(:)
    integer, intent(in)         :: index

    if (allocated(newext)) deallocate(newext)
    allocate(newext(size(stored_ext, 2)))

    newext = stored_ext(index, :)

  end subroutine extrema_pull
end module extrema_storage
