program create_listgal
  implicit none
  integer                                 ::n,i,j,k,nbodies,nb_of_parts,npoints,mynumber
  integer                                 ::nx,ny,nz,ilevel,idim,jdim,kdim,icell
  integer                                 ::nlevelmax,ilevel1,ngrid1
  integer                                 ::nlevelmaxs,nlevel,iout
  integer                                 ::ind,ipos,ngrida,ngridh,ilevela,ilevelh
  integer                                 ::ngridmax,nstep_coarse,icpu,ncpu_read
  integer                                 ::nhx,nhy,ihx,ihy,ivar1,ivar2,mylevel,hosthalo,hostsub,nbsub,nextsub
  real                                    ::gamma,smallr,smallc,gammah
  real                                    ::boxlen,boxlen2
  real                                    ::t,aexp,hexp,t2,aexp2,hexp2
  real                                    ::omega_m,omega_l,omega_k,omega_b
  real                                    ::scale_l,scale_d,scale_t
  real                                    ::omega_m2,omega_l2,omega_k2,omega_b2

  integer                                 ::ngrid,imin,imax,jmin,jmax,kmin,kmax
  integer                                 ::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer                                 ::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer                                 ::nx_full,ny_full,nz_full,lmin,levelmin
  character(LEN=5)                        ::nchar,ncharcpu
  character(LEN=80)                       ::ordering
  character(LEN=128)                      ::nomfich,repository,outfich,filetype='bin',fileDM,filegal
  logical                                 ::ok,ok_part,ok_cell
  real(kind=4)                            ::mhalo,px,py,pz,vx,vy,vz,Lx,Ly,Lz,Lnorm,rad,dum
  real(kind=4)                            ::age_univ,rvir,mvir,tvir,cvel,unit_l,unit_d,unit_t,dr,dr2,radius,r2close,time
  integer                                 ::ngal,nb_of_gal,nb_of_subgal
  integer     ,allocatable,dimension(:)   ::idgal,members
  real(kind=8),allocatable,dimension(:)   ::mgal,xgal,ygal,zgal,Lxgal,Lygal,Lzgal,rgal
  real(kind=4),allocatable,dimension(:,:) ::toto
  real(kind=4)                            ::mstar
  real(kind=4),allocatable,dimension(:)   ::rr,density
  integer,allocatable,dimension(:)        ::hlevel

  call read_params
  call read_info

  !-----------------------------------------------
  ! Lecture du fichier contenant les galaxies
  !-----------------------------------------------
  write(*,*)filegal
  open(unit=1,file=filegal,status='old',form='unformatted')
  read(1)nbodies
  read(1)
  read(1)
  read(1)
  read(1)
  read(1)nb_of_gal,nb_of_subgal
  ngal=nb_of_gal+nb_of_subgal
  write(*,*)ngal,' galaxies'
  allocate(mgal(1:ngal),xgal(1:ngal),ygal(1:ngal),zgal(1:ngal),Lxgal(1:ngal),Lygal(1:ngal),Lzgal(1:ngal),idgal(1:ngal))
  allocate(hlevel(1:ngal),rgal(1:ngal))
  do i=1,ngal
    read(1)nb_of_parts
    allocate(members(1:nb_of_parts))
    read(1)members
    deallocate(members)
    ! Read properties of each galaxies
    read(1)mynumber
    idgal(i)=mynumber
    read(1)
    read(1)mylevel,hosthalo,hostsub,nbsub,nextsub
    read(1)mstar
    mgal(i)=mstar*1d11
    read(1)px,py,pz
    xgal(i)=px
    ygal(i)=py
    zgal(i)=pz
    read(1)
    read(1)Lx,Ly,Lz
    Lnorm=sqrt(Lx*Lx+Ly*Ly+Lz*Lz)
    Lxgal(i)=Lx/Lnorm
    Lygal(i)=Ly/Lnorm
    Lzgal(i)=Lz/Lnorm
    hlevel(i)=mylevel
    read(1)rad,dum,dum,dum
    read(1)
    read(1)
    read(1)
    read(1)
    read(1)
    read(1)npoints
    allocate(rr(1:npoints),density(1:npoints))
    read(1)rr
    read(1)density
    deallocate(rr,density)
    rgal(i)=rad
 enddo
!!$
!!$  ipos=INDEX(repository,'output_')
!!$  nchar=repository(ipos+7:ipos+13)
!!$  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
!!$  open(unit=10,file=nomfich,form='formatted',status='old')
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,'("levelmin    =",I11)')levelmin
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,'("time        =",E23.15)')t
!!$  write(*,*)'time=',t
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,*)
!!$  read(10,'("unit_l      =",E23.15)')unit_l
!!$  read(10,'("unit_d      =",E23.15)')unit_d
!!$  read(10,'("unit_t      =",E23.15)')unit_t
!!$  write(*,*)'unit_l=',unit_l
!!$  write(*,*)'unit_d=',unit_d
!!$  write(*,*)'unit_t=',unit_t
!!$  close(10)

  ! Convert Mpc units into box units
  xgal=xgal*1d6*3.08d18/unit_l+0.5d0
  ygal=ygal*1d6*3.08d18/unit_l+0.5d0
  zgal=zgal*1d6*3.08d18/unit_l+0.5d0
  rgal=rgal*1d6*3.08d18/unit_l

  write(*,*)'==============================='
  write(*,*)'==============================='
  write(*,*)'==============================='
  ! Output file
  nomfich=TRIM(outfich)//'.bin'
  write(*,*)'Writing file '//TRIM(nomfich)
  open(unit=20,file=nomfich,form='unformatted')
  write(20)ngal,8
  allocate(toto(ngal,8))
  toto(1:ngal,1)=float(idgal(1:ngal))
  toto(1:ngal,2)=mgal(1:ngal)
  toto(1:ngal,3)=xgal(1:ngal)
  toto(1:ngal,4)=ygal(1:ngal)
  toto(1:ngal,5)=zgal(1:ngal)
  toto(1:ngal,6)=Lxgal(1:ngal)
  toto(1:ngal,7)=Lygal(1:ngal)
  toto(1:ngal,8)=Lzgal(1:ngal)
  write(20)toto
  close(20)

  write(*,*)'==============================='
  write(*,*)'==============================='
  write(*,*)'==============================='
  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Writing file '//TRIM(nomfich)
  open(unit=20, file=nomfich)
  write(20,'(2i9)') ngal, 7
  do i=1,ngal
     write(20,'(2I9,1ES10.3,4f12.8)')idgal(i),hlevel(i),mgal(i),xgal(i),ygal(i),zgal(i),rgal(i)
  enddo
  close(20)

  deallocate(mgal)
  deallocate(xgal,ygal,zgal)
  deallocate(idgal,Lxgal,Lygal,Lzgal)
  deallocate(toto)

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
       print *, 'usage: create_listgal  -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-fil filetype] '
       print *, 'ex: create_listgal -inp output_00001 -out listgal.dat'// &
            &   ' -fga tree_bricks001'
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
       case ('-fil')
          filetype = trim(arg)
       case ('-fga')
          filegal=trim(arg)
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

    return

  end subroutine read_params

end program create_listgal

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
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

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d
!================================================================
!================================================================
!================================================================
!================================================================
