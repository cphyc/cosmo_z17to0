program create_halolist
  !--------------------------------------------------------------------------
  ! This program build up a synthetic list of halos with positions, mass and radius
  ! Version by Y. Dubois
  !--------------------------------------------------------------------------
  implicit none
  integer::i,idh=-1
  integer::nbodies,nb_of_halos,nb_of_subhalos,nb_of_parts,my_number,nbsub
  integer::hosthalo,hostsub,mylevel,nextsub,nDM_tmp,nDM,pmin=0
  integer,allocatable,dimension(:)::members,hlevel,hlevel_tmp,ind_DM,ind_DM_tmp
  real(kind=4)::mhalo,px,py,pz,vx,vy,vz,tvir,cvel,mvir,xx,yy,zz,rvir,aexp,age_univ
  real(kind=4),allocatable,dimension(:)::mDM,xDM,yDM,zDM,rvirDM,mvirDM
  real(kind=4),allocatable,dimension(:)::mDM_tmp,xDM_tmp,yDM_tmp,zDM_tmp,rvirDM_tmp,mvirDM_tmp
  real(kind=4)::unit_l,unit_t,unit_d,time,mmin=0d0,mmax=1d20,xmin,xmax,ymin,ymax,zmin,zmax,rad,dum
  character(LEN=128)::nomfich,repository,outfich,fileDM,filegal='none'
  logical::subhalo=.false.
  integer::npoints
  real(kind=4),allocatable,dimension(:)::rr
  real(kind=4),allocatable,dimension(:,:) :: tmp

  call read_params
  call read_info

  !-----------------------------------------------
  ! Reading file containing halos
  !-----------------------------------------------
  if(filegal.ne.'none')fileDM=filegal
  write(*,*)fileDM
  open(unit=1,file=fileDM,status='old',form='unformatted')
  read(1)nbodies
  read(1)
  read(1)aexp
  read(1)
  read(1)age_univ
  read(1)nb_of_halos,nb_of_subhalos
  nDM_tmp=nb_of_halos+nb_of_subhalos
  write(*,*)nDM_tmp,' halos in total'
  write(*,*)nb_of_halos,' main halos'
  write(*,*)nb_of_subhalos,' sub halos'
  allocate(mDM_tmp(1:nDM_tmp),xDM_tmp(1:nDM_tmp),yDM_tmp(1:nDM_tmp),zDM_tmp(1:nDM_tmp) &
       & ,rvirDM_tmp(1:nDM_tmp),mvirDM_tmp(1:nDM_tmp) &
       & ,hlevel_tmp(1:nDM_tmp),ind_DM_tmp(1:nDM_tmp))
  nDM=0
  do i=1,nDM_tmp
     read(1)nb_of_parts
     allocate(members(1:nb_of_parts))
     read(1)members
     if(i.eq.idh)then
        open(unit=10,file='members.dat',form='unformatted')
        write(10)nb_of_parts
        write(10)members
        close(10)
     endif
     deallocate(members)
     ! Read properties of each halo
     read(1)my_number
     read(1)
     read(1)mylevel,hosthalo,hostsub,nbsub,nextsub
     read(1)mhalo
     read(1)px,py,pz
     read(1)
     read(1)
     read(1)rad,dum,dum,dum
     read(1)
     read(1)
     if(filegal.ne.'none')read(1)
     read(1)rvir,mvir,tvir,cvel
     read(1)
     if(filegal.ne.'none')then
        read(1)
        read(1)
        read(1)
     endif
     xx=px*1d6*3.08d18/unit_l+0.5d0
     yy=py*1d6*3.08d18/unit_l+0.5d0
     zz=pz*1d6*3.08d18/unit_l+0.5d0
     ! Reduce the number of halos to the number some user specified constraints
     if( ((.not.subhalo) .and. mylevel.eq.1) .or. (subhalo))then
     if (mvir*1d11.ge.mmin .and. mvir*1d11.lt.mmax .and.nb_of_parts.gt.pmin)then
        nDM=nDM+1
        ind_DM_tmp(nDM)=my_number
        hlevel_tmp(nDM)=mylevel
        mDM_tmp(nDM)=mhalo*1d11
        xDM_tmp(nDM)=xx
        yDM_tmp(nDM)=yy
        zDM_tmp(nDM)=zz
        if(filegal.ne.'none')then
           rvirDM_tmp(nDM)=rad*3.08d24/unit_l
        else
           rvirDM_tmp(nDM)=rvir*3.08d24/unit_l
        endif
        mvirDM_tmp(nDM)=mvir*1d11
     endif
     endif
  enddo
  write(*,*)nDM,' halos within the boundaries'
  allocate(mDM(1:nDM),xDM(1:nDM),yDM(1:nDM),zDM(1:nDM) &
       & ,rvirDM(1:nDM),mvirDM(1:nDM),hlevel(1:nDM),ind_DM(1:nDM))
  hlevel(1:nDM)=hlevel_tmp(1:nDM)
  ind_DM(1:nDM)=ind_DM_tmp(1:nDM)
  mDM(1:nDM)=mDM_tmp(1:nDM)
  xDM(1:nDM)=xDM_tmp(1:nDM)
  yDM(1:nDM)=yDM_tmp(1:nDM)
  zDM(1:nDM)=zDM_tmp(1:nDM)
  rvirDM(1:nDM)=rvirDM_tmp(1:nDM)
  mvirDM(1:nDM)=mvirDM_tmp(1:nDM)
  deallocate(mDM_tmp,xDM_tmp,yDM_tmp,zDM_tmp,rvirDM_tmp,mvirDM_tmp,hlevel_tmp)

  write(*,*)'==============================='
  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Writing file '//TRIM(nomfich)
  open(unit=20,file=nomfich)
  write(20, *) nDM, 7
  do i=1,nDM
     write(20,'(2I9,1ES10.3,4f12.8)')ind_DM(i),hlevel(i),mvirDM(i),xDM(i),yDM(i),zDM(i),rvirDM(i)
  enddo
  close(20)

  allocate(tmp(nDM,7))
  tmp(1:nDM, 1) = float(ind_DM)
  tmp(1:nDM, 2) = float(hlevel(i))
  tmp(1:nDM, 3) = mvirDM
  tmp(1:nDM, 4) = xDM
  tmp(1:nDM, 5) = yDM
  tmp(1:nDM, 6) = zDM
  tmp(1:nDM, 7) = rvirDM
  nomfich=TRIM(outfich)//'.bin'
  write(*,*)'Writing file '//TRIM(nomfich)
  open(unit=20, file=nomfich, form='unformatted')
  write(20) nDM, 7
  write(20) tmp
  deallocate(tmp)

contains

  subroutine read_info
    implicit none

    integer::ipos
    character(LEN=5)::nchar

    !=======================
    ! Read info file
    !=======================
    ipos=INDEX(repository,'output_')
    nchar=repository(ipos+7:ipos+13)
    nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
    open(unit=10,file=nomfich,form='formatted',status='old')
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,'("time        =",E23.15)')time
    write(*,*)'time=',time
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,'("unit_l      =",E23.15)')unit_l
    read(10,'("unit_d      =",E23.15)')unit_d
    read(10,'("unit_t      =",E23.15)')unit_t
    write(*,*)'unit_l=',unit_l
    write(*,*)'unit_d=',unit_d
    write(*,*)'unit_t=',unit_t
    close(10)
    !=======================
    ! End read info file
    !=======================
  end subroutine read_info

  subroutine read_params

    implicit none

    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok

    n = iargc()
    if (n < 4) then
       print *, 'usage: create_halolist -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-fil filetype] '
       print *, 'ex: create_halolist -inp output_00001 -fDM tree_brick001 -out list_halo.dat'// &
            &   ' -xmi 0.1 -xma 0.7 -mmi 1d11 -mma 1d12'
       stop
    end if

    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case ('-inp')
          repository = trim(arg)
       case ('-out')
          outfich = trim(arg)
       case ('-xmi')
          read (arg,*) xmin
       case ('-xma')
          read (arg,*) xmax
       case ('-ymi')
          read (arg,*) ymin
       case ('-yma')
          read (arg,*) ymax
       case ('-zmi')
          read (arg,*) zmin
       case ('-zma')
          read (arg,*) zmax
       case ('-mmi')
          read (arg,*) mmin
       case ('-mma')
          read (arg,*) mmax
       case ('-sub')
          read (arg,*) subhalo
       case ('-pmi')
          read (arg,*) pmin
       case ('-idh')
          read (arg,*) idh
       case ('-fDM')
          fileDM=trim(arg)
       case ('-fga')
          filegal=trim(arg)
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

    return

  end subroutine read_params


end program create_halolist
