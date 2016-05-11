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
  integer            :: param_from, param_to, output_number, param_halo_i
  integer, allocatable :: param_output_number_list(:)
  integer, parameter :: NPARTICLE_TO_PROBE_HALO = 50, NCPU_PER_HALO = 10
  real(kind=8)       :: param_min_m, param_max_m
  integer            :: param_verbosity
  logical            :: stop_flag

  type(command_line_interface) :: cli
  type(CPU_INFOS_T) :: infos
  !----------------------------------------
  ! halo centers
  !----------------------------------------
  real(kind=8), allocatable, dimension(:, :) :: centers
  integer :: lines, cols

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
  real(kind=8), dimension(3)                :: X0, X1, center

  !-------------------------------------
  ! random
  !-------------------------------------
  call random_seed()

  !-------------------------------------
  ! parameters
  !-------------------------------------
  call cli%init(progname='Get the particles in a region around')
  call cli%add(switch='--output-path', help='Path of the output', &
       act='store', def='/data52/Horizon-AGN/OUTPUT_DIR')
  call cli%add(switch='--output-number', help='List of the outputs', &
       act='store', def='2', nargs='+')
  call cli%add(switch='--center', help='Center of region', &
       act='store', nargs='3', required=.true.)
  call cli%add(switch='--region', help='Size of the region', &
       act='store', nargs='3', required=.true.)
  call cli%add(switch='--verbose', help='Verbosity', &
       act='store', def='-1')
  call cli%add(switch='--output', help='Prefix for output of this program', &
       act='store', required=.true.)
  call cli%add(switch='--halo-centers', help='List of the centers of the halo at each output')


  call cli%get(switch='--output-path',           val=param_output_path)
  call cli%get_varying(switch='--output-number', val=param_output_number_list)
  call cli%get(switch='--verbose',               val=param_verbosity)
  call cli%get(switch='--output',                val=param_output_prefix)
  call cli%get(switch='--halo-centers',          val=param_halo_centers)

  !----------------------------------------
  ! Reading halo centers
  !----------------------------------------
  call read_list_header(param_halo_centers, unit, lines, cols)
  allocate(centers(lines, cols))
  call read_list_data(unit, lines, cols, centers)

  !-------------------------------------
  ! iterate over each output
  !-------------------------------------
  do i = 1, size(param_output_number_list)
     output_number = param_output_number_list(i)

     ! read the information file of the output
     write(tmp_char, '(a,a,i5,a,i5,a)') trim(param_output_path), '/output_', output_number, &
          '/info_', output_number, '.txt'
     call read_info_headers(tmp_char, infos)


     ! read the 
  end do


end program compute_halo_prop
