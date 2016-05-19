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
  integer, parameter :: NPART = 2**30
  character(len=200) :: ramses_output_end, ramses_output_start, gal_list_filename, &
       associations_filename, dm_halo_list_filename, brick_file, info_file_end,&
       info_file_start, outfile, halo_to_cpu_file, mergertree_file, param_output_path, &
       param_output
  integer            :: param_from, param_to, param_output_number, param_halo_i
  integer, allocatable :: param_output_number_list(:), param_halos(:)
  integer, parameter :: NPARTICLE_TO_PROBE_HALO = 50, NCPU_PER_HALO = 10
  real(kind=8)       :: param_min_m, param_max_m
  integer            :: param_verbosity
  logical            :: stop_flag

  type(command_line_interface) :: cli
  !-------------------------------------
  ! List data
  !-------------------------------------
  real(kind=4), dimension(:,:), allocatable :: data_gal, data_associations, data_halo
  integer                                   :: ngal, nassoc, ndm_halo, gal_cols, assoc_cols, dm_halo_cols
  !-------------------------------------
  ! Brick data
  !-------------------------------------
  real(kind=8), dimension(:, :), allocatable :: posDM, LDM
  integer, dimension(:), allocatable         :: idDM
  real(kind=8), dimension(:), allocatable    :: rvirDM, mDM, mvirDM, TvirDM, hlevel
  type(INFOS_T)                              :: infos
  type(MEMBERS_T), dimension(:), allocatable :: members
  integer                                    :: nbodies, nb_of_halos, nb_of_subhalos, nDM
  real(kind=4)                               :: aexp_tmp, age_univ
  !-------------------------------------
  ! Merger tree data
  !-------------------------------------
  ! integer,  allocatable, dimension(:)        :: mt_nhalos, mt_nsubhalos
  ! real(kind=4), dimension(:), allocatable    :: mt_aexp, mt_omega_t, mt_age_univ, time
  ! real(kind=4), dimension(:, :), allocatable :: mt_pos, mt_vel, initial_pos
  ! integer                                    :: mt_nsteps, nsteps, nhalos
  ! integer                                    :: initial_halo_i, current_halo, parent_halo, istep, prev, father, step, max_nhalo
  ! integer, dimension(:), allocatable         :: tmp_nhalos, halos_z0
  ! integer, dimension(:), allocatable         :: parent, parent_at_step
  !-------------------------------------
  ! Particle data
  !-------------------------------------
  integer                                    :: ndim, nparts
  real(kind=8), dimension(:, :), allocatable :: vel, tmp_pos
  real(kind=8), dimension(3, NPART)          :: pos, pos_around_halo
  real(kind=4), dimension(:), allocatable    :: birth_date
  real(kind=8), dimension(:), allocatable    :: m
  integer                                    :: nstar, halo_found
  integer,      dimension(:), allocatable    :: tmp_ids
  integer,      dimension(NPART)             :: ids, order, ids_around
  integer, dimension(:), allocatable         :: order2
  !-------------------------------------
  ! Halo to cpu
  !-------------------------------------
  integer                               :: n_cpu_per_halo
  integer, dimension(:, :), allocatable :: halo_to_cpu
  integer, dimension(:), allocatable    :: cpu_list
  !-------------------------------------
  ! Halo properties
  !-------------------------------------
  real(kind=8), dimension(3, 3)               :: I_t_diag
  real(kind=8), dimension(:,:,:), allocatable :: I_t
  real(kind=8), dimension(:), allocatable     :: m_in_halo, m_in_box
  real(kind=8), dimension(:, :), allocatable  :: pos_in_halo, prev_pos_in_halo, vel_in_halo
  integer, dimension(:), allocatable          :: ids_in_box, ids_in_halo, prev_ids_in_halo
  real(kind=8), dimension(:, :), allocatable  :: pos_in_box
  real(kind=8)                                :: mtot
  real(kind=8), dimension(3)                  :: pos_span, pos_mean
  integer                                     :: halo_i, ntot_halo, halo_counter, part_i, parts_in_region
  integer                                     :: part_counter
  integer, dimension(:), allocatable          :: halo_found_mask
  !-------------------------------------
  ! Tmp variables
  !-------------------------------------
  integer                                   :: i, j, k, cpu, i1, i2, counter, x, y, z
  integer                                   :: tmp_int, unit,&
       & out_unit, center_output_unit, tmp_int2, index
  character(len=200)                        :: tmp_char
  real                                      :: tmp_real, factor
  logical                                   :: tmp_bool
  real(kind=8), allocatable, dimension(:)   :: tmp_arr, tmp_arr2
  real(kind=8), allocatable, dimension(:,:) :: tmp_dblarr
  integer, allocatable, dimension(:)        :: tmp_iarr
  logical, dimension(4096)                  :: cpu_read
  real(kind=8), dimension(3)                :: X0, X1, center
  real(kind=8), dimension(:, :), allocatable :: centers
  real(kind=8)                              :: margin
  integer                                   :: nstep, max_nparts, step_i

  !-------------------------------------
  ! random
  !-------------------------------------
  call random_seed()

  !-------------------------------------
  ! Read parameters
  !-------------------------------------
  call cli%init(progname='find_particles_in_halo',&
       description='Find particles around many halos')

  call cli%add(switch='--output-path', help='Path of the output', &
       act='store', def='/data52/Horizon-AGN/OUTPUT_DIR')
  call cli%add(switch='--output-number', help='Number of the output', &
       act='store', def='2', nargs='+')
  call cli%add(switch='--verbose', help='Verbosity', &
       act='store', def='0')
  call cli%add(switch='--halos', help='List of halos to find', act='store', &
       nargs='+', def='1')
  call cli%add(switch='--brick', help='Path for brick file that contains the halos', &
       act='store', def='/data52/Horizon-AGN/TREE_DM_celldx2kpc_SC0.9r/tree_bricks782')
  call cli%add(switch='--output', switch_ab='-o', help='Name of the output file', &
       act='store', def='data.out')


  call cli%get(switch='--output-path', val=param_output_path)
  call cli%get_varying(switch='--output-number', val=param_output_number_list)
  call cli%get(switch='--verbose', val=param_verbosity)
  call cli%get_varying(switch='--halos', val=param_halos)
  call cli%get(switch='--output', val=param_output)

  !-------------------------------------
  ! Read brick file
  !-------------------------------------
  print*, ''
  print*, 'Reading info file'
  call read_info_headers(param_output_path, param_output_number_list(1), infos)

  call cli%get(switch='--brick', val=tmp_char)
  print*, ''
  print*, 'Reading brick file "' // trim(tmp_char) // '"…'
  call read_brick_header(tmp_char, infos, nbodies, aexp_tmp, age_univ,&
       nb_of_halos, nb_of_subhalos)
  nDM = nb_of_halos + nb_of_subhalos
  allocate(idDM(nDM))
  allocate(posDM(infos%ndim, nDM))
  allocate(rvirDM(nDM))
  allocate(mDM(nDM))
  allocate(mvirDM(nDM))
  allocate(TvirDM(nDM))
  allocate(hlevel(nDM))
  allocate(LDM(infos%ndim, nDM))
  allocate(members(nDM))
  allocate(cpu_list(infos%ncpu))
  print*, '    …allocated…'
  call read_brick_data(nDM, infos, .true., &
       & mDM, posDM, rvirDM, mvirDM, TvirDM,&
       & hlevel, LDM, idDM, members)
  print*, '    …red!'
  deallocate(rvirDM, mvirDM, TvirDM, hlevel, LDM)


  counter = 1
  !----------------------------------------
  ! Read all the cpus
  !----------------------------------------
  do cpu = 1, infos%ncpu
     if (param_verbosity >= 3) then
        write(*, '(a, i5, a, i5, a, 1x, i10, a, i10, a)') 'cpu', cpu, '/', infos%ncpu,&
             '(parts', counter, '/', NPART, ')'
     end if

     call read_particle(param_output_path, 2, cpu, nstar, &
          tmp_pos, vel, m, tmp_ids, birth_date, ndim, nparts)
     pos(:, counter:counter + nparts) = tmp_pos
     ids(counter:counter + nparts) = tmp_ids
     ! add the particles into the array
     counter = counter + nparts
  end do

  !----------------------------------------
  ! Sort the indexes
  !----------------------------------------
  if (param_verbosity >= 2) then
     write(*, *) 'Sorting'
  end if
  call quick_sort(ids, order)

  !----------------------------------------
  ! Iterate over each given halo
  !----------------------------------------
  do i = 1, size(param_halos)
     halo_i = param_halos(i)

     if (param_verbosity >= 2) then
        write(*, *) 'Getting halo', halo_i
     end if

     ! allocate data for halo
     allocate(pos_in_halo(infos%ndim, members(halo_i)%parts))

     !----------------------------------------
     ! get the members of the halo
     !----------------------------------------
     do j = 1, members(halo_i)%parts
        tmp_int = indexOf(members(halo_i)%ids(j), ids)

        pos_in_halo(:, j) = pos(:, order(tmp_int))
     end do
     !----------------------------------------
     ! Rectify the position and get center
     !----------------------------------------
     call correct_positions(pos_in_halo)
     call meanval(pos_in_halo, pos_mean, infos%ndim, members(halo_i)%parts)
     pos_span = maxval(pos_in_halo, 2) - minval(pos_in_halo, 2)

     !----------------------------------------
     ! get the particles around halo
     !----------------------------------------
     call filter_region(pos_mean, 2*pos_span, &
          ids, pos, &
          ids_around, pos_around_halo, &
          parts_in_region, order)
     call correct_positions(pos_around_halo(:, :parts_in_region))

     if (param_verbosity >= 3) write(*, *) 'Found', parts_in_region, 'around halo'

     !----------------------------------------
     ! write data
     !----------------------------------------
     write(tmp_char, '(a,i0.5)') param_output_path, halo_i
     if (param_verbosity >= 2) write(*, *) 'Writing output to', trim(tmp_char)
     open(newunit=out_unit, file=trim(tmp_char), form='formatted')

     ! headers
     write(out_unit, '(a9, a5, 6a13)') 'id', 'in', 'x', 'y', 'z'

     ! in halo
     do k = 1, members(halo_i)%parts
        write(out_unit, '(i12, L5, 6ES14.6e2)') &
             members(halo_i)%ids(k), .true., pos_in_halo(:, k)
     end do

     ! out of halo (skipping parts in halo)
     allocate(order2(members(halo_i)%parts))
     call quick_sort(members(halo_i)%ids, order2)

     ! for particles around …
     do k = 1, parts_in_region
        ! … that aren't in the halo …
        if (indexOf(ids_around(k), members(halo_i)%ids) == -1) then
           ! … write their data to the file
           write(out_unit, '(i12, L5, 6ES14.6e2)') &
                ids_around(k), .false., pos_around_halo(:, k)
        end if
     end do
     close(out_unit)

     deallocate(order2, pos_in_halo)

  end do

end program compute_halo_prop
