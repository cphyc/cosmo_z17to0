program compute_halo_prop
  use io
  use misc
  use flap, only : command_line_interface
  use compute
  ! use hdf5
  implicit none

  !-------------------------------------
  ! Parameters
  !-------------------------------------
  character(len=200) :: param_output_prefix, param_output_path, param_halo_centers
  integer            :: param_output_number, param_halo_i

  real(kind=8), dimension(3) :: center, region_size

  integer, parameter :: NPARTICLE_TO_PROBE_HALO = 50, NCPU_PER_HALO = 10
  real(kind=8)       :: param_min_m, param_max_m
  real, dimension(:), allocatable :: param_center, param_region_size
  integer            :: param_verbosity
  logical            :: stop_flag

  type(command_line_interface) :: cli
  type(PARTICLE_DATA), dimension(:), allocatable :: data
  type(INFOS_T) :: infos

  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                   :: i, j, cpu, i1, i2, counter, x, y, z
  integer                                   :: tmp_int, unit, tmp_int2, index
  character(len=200)                        :: tmp_char
  real                                      :: tmp_real, factor
  logical                                   :: tmp_bool
  real(kind=8), allocatable, dimension(:)   :: tmp_arr, tmp_arr2
  real(kind=8), allocatable, dimension(:,:) :: tmp_dblarr
  integer, allocatable, dimension(:)        :: tmp_iarr
  logical, dimension(4096)                  :: cpu_read

  !-------------------------------------
  ! random
  !-------------------------------------
  call random_seed()

  !-------------------------------------
  ! parameters
  !-------------------------------------
  call cli%init(progname='Get the particles in a region around center')
  call cli%add(switch='--output-path', help='Path to the outputs', &
       act='store', def='/data52/Horizon-AGN/OUTPUT_DIR')
  call cli%add(switch='--output-number', help='Output number', &
       act='store', def='2')
  call cli%add(switch='-x', act='store', required=.true.)
  call cli%add(switch='-y', act='store', required=.true.)
  call cli%add(switch='-z', act='store', required=.true.)
  call cli%add(switch='-dx', act='store', required=.true.)
  call cli%add(switch='-dy', act='store', required=.true.)
  call cli%add(switch='-dz', act='store', required=.true.)

  call cli%add(switch='--verbose', help='Verbosity', &
       act='store', def='-1')
  call cli%add(switch='--output', help='Prefix for output of this program', &
       act='store', required=.true.)


  call cli%get(switch='--output-path',    val=param_output_path)
  call cli%get(switch='--output-number',  val=param_output_number)
  call cli%get(switch='--verbose',        val=param_verbosity)
  call cli%get(switch='--output',         val=param_output_prefix)
  call cli%get(switch='-x', val=center(1))
  call cli%get(switch='-y', val=center(2))
  call cli%get(switch='-z', val=center(3))
  call cli%get(switch='-dx', val=region_size(1))
  call cli%get(switch='-dy', val=region_size(2))
  call cli%get(switch='-dz', val=region_size(3))

  !-------------------------------------
  ! Iterate over each output
  !-------------------------------------
  call read_info_headers(param_output_path, param_output_number, infos)

  call read_region(center, region_size, infos, data)

  write(*,*) 'Output to ', trim(param_output_prefix)
  open(newunit=unit, file=trim(param_output_prefix), form='formatted')
  write(unit, '(a10, a7, 6a14)') 'id', 'cpu', 'x', 'y', 'z'
  do i = 1, size(data)
     do j = 1, size(data(i)%ids)
        write(unit, '(i10, i7, 3ES14.6e2)') data(i)%ids(j), data(i)%cpu, data(i)%pos(:, j)
     end do
  end do

  write(*,*) 'Output to ', trim(param_output_prefix) // '.bin'
  open(newunit=unit, file=trim(param_output_prefix) // '.bin', form='unformatted')
  tmp_int = 1
  do i = 1, size(data)
     tmp_int = tmp_int + size(data(i)%ids)
  end do

  allocate(tmp_dblarr(3, tmp_int), tmp_iarr(tmp_int))

  tmp_int2 = 1
  do i = 1, size(data)
     do j = 1, size(data(i)%ids)
        tmp_iarr(tmp_int) = data(i)%ids(j)
        tmp_dblarr(:, tmp_int) = data(i)%pos(:, j)

        tmp_int2 = tmp_int2 + 1
     end do
  end do
  write(unit) 2, tmp_int2
  write(unit) tmp_iarr
  write(unit) tmp_dblarr

  close(unit)
end program compute_halo_prop
