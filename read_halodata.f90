program read_halodata
  implicit none
  integer::nhalo,io,idh,sub,ihalo,nvarh,ncell,ndm,nstar,iread,nchemo,i,index,nread
  real(kind=8)::mh,xh,yh,zh,rh,xoff,yoff,zoff,dxx,dyy,dzz,dr2,mtot_gas,mtot_dm,mtot_star
  real(kind=8)::aexp,scale_l,scale_d,scale_t,msun
  integer,dimension(:),allocatable::idhalo,idhalo2,subhalo
  real(kind=8),dimension(:),allocatable::mhalo,xhalo,yhalo,zhalo,rvir,rhalo,r2halo
  integer,dimension(:),allocatable::idp,idm,istar
  real(kind=8),dimension(:),allocatable::m,mp,zp,age,dxcell,varread,mdm,mstar,zstar,tstar,chread
  real(kind=8),dimension(:),allocatable::rcell,rdm,rstar
  real(kind=8),dimension(:,:),allocatable::chp,v,xcell,var,xdm,vdm,xstar,vstar,chstar
  real(kind=8),dimension(:),allocatable::dmcell,dmdm,dmstar,dmtot,radius,density
  character(LEN=128)::listhalo='none',nomfich,repository,outdir,outfich
  character(LEN=7)::ncharlong
  integer::ihread=0,nbin=500
  real(kind=8)::xread,yread,zread,dxread,vxread,vyread,vzread,mread,tread,zmread
  real(kind=8)::drbin,omega_m=0.272,rhoc
  logical::chemo=.true.,ok

  call read_params
  ncell=0; ndm=0; nstar=0

  if(listhalo.ne.'none')then
     nhalo=0
     open(unit=1,file=listhalo)
     do
        read(1,*,iostat=io)idh,sub,mh,xh,yh,zh,rh
        if(io > 0)then
           write(*,*)'Something bad happened while reading ',listhalo
           stop
        elseif(io < 0)then
           !write(*,*)'Finished reading ',listhalo
           exit
        else
           nhalo=nhalo+1
        endif
     enddo
     close(1)
     write(*,*)'Found ',nhalo,' halo'

     allocate(idhalo(1:nhalo),idhalo2(1:nhalo),subhalo(1:nhalo),mhalo(1:nhalo))
     allocate(xhalo(1:nhalo),yhalo(1:nhalo),zhalo(1:nhalo),rvir(1:nhalo),rhalo(1:nhalo),r2halo(1:nhalo))
     nhalo=0
     open(unit=1,file=listhalo)
     do
        read(1,*,iostat=io)idh,sub,mh,xh,yh,zh,rh
        if(io > 0)then
           write(*,*)'Something bad happened while reading ',listhalo
           stop
        elseif(io < 0)then
           !write(*,*)'Finished reading ',listhalo
           exit
        else
           nhalo=nhalo+1
           idhalo (nhalo)=idh
        endif
     enddo
     close(1)

  else

     if(ihread.eq.0)then
        write(*,*)'You need to give me one halo to read at least (-idh 1)'
        stop
     else
        nhalo=1
        allocate(idhalo(1:nhalo),idhalo2(1:nhalo),subhalo(1:nhalo),mhalo(1:nhalo))
        allocate(xhalo(1:nhalo),yhalo(1:nhalo),zhalo(1:nhalo),rvir(1:nhalo),rhalo(1:nhalo),r2halo(1:nhalo))
        idhalo(nhalo)=ihread
     endif

  endif

  open(unit=99,file=outfich,form='unformatted')
  write(99)nhalo,nbin

  ! Loop over halos
  do ihalo=1,nhalo
     call title2(idhalo(ihalo),ncharlong)

     !-----------------------------------------------
     ! Read hydro data
     !----------------------------------------------
     nomfich=TRIM(repository)//'/halo_'//TRIM(ncharlong)//'/cells.dat'
     INQUIRE(file=nomfich,exist=ok)
     if(.not.ok)then
        write(*,*)'Skipping file ',TRIM(nomfich)
     else
        write(*,*)'Reading file ',TRIM(nomfich)
        open(unit=10,file=nomfich,status='old',form='unformatted')
        read(10)idhalo2(ihalo),subhalo(ihalo)
        if(idhalo2(ihalo).ne.idhalo(ihalo))then
           write(*,'(A,I10,I10)')'WARNING: Halo id in the file is different from the name file',idhalo2(ihalo),idhalo(ihalo)
           write(*,*)'Something suspicious is happening! I continue nonetheless...'
        endif
        read(10)mhalo(ihalo)
        read(10)xhalo(ihalo)
        read(10)yhalo(ihalo)
        read(10)zhalo(ihalo)
        read(10)rvir(ihalo),rhalo(ihalo)
        read(10)aexp
        read(10)scale_l,scale_d,scale_t
        read(10)nvarh
        read(10)ncell,nread
        write(*,*)'ncell=',ncell
        allocate(xcell(1:ncell,1:3),dxcell(1:ncell),var(1:ncell,1:nvarh),rcell(1:ncell))
        allocate(varread(1:nvarh))
        read(10)xcell(1:ncell,1)
        read(10)xcell(1:ncell,2)
        read(10)xcell(1:ncell,3)
        read(10)dxcell(1:ncell)
        do i=1,nvarh
           read(10)var(1:ncell,i)
        enddo
        close(10)

        do i=1,ncell
           xread=xcell(i,1)
           yread=xcell(i,2)
           zread=xcell(i,3)
           ! Check for periodicity
           dxx=xhalo(ihalo)-xread
           if(dxx> 0.5)then
              dxx=dxx-1d0
           endif
           if(dxx<-0.5)then
              dxx=dxx+1d0
           endif
           dr2=dxx*dxx
           dyy=yhalo(ihalo)-yread
           if(dyy> 0.5)then
              dyy=dyy-1d0
           endif
           if(dyy<-0.5)then
              dyy=dyy+1d0
           endif
           dr2=dyy*dyy+dr2
           dzz=zhalo(ihalo)-zread
           if(dzz> 0.5)then
              dzz=dzz-1d0
           endif
           if(dzz<-0.5)then
              dzz=dzz+1d0
           endif
           dr2=dzz*dzz+dr2
           rcell(i)=sqrt( dr2 )
        enddo
        msun=scale_d*scale_l**3d0/2d33
        write(*,'(A,ES13.6)')'Mgas =',sum(var(1:ncell,1)*dxcell**3d0*msun)
        !-----------------------------------------------
        ! End read hydro data
        !-----------------------------------------------
     endif

     !-----------------------------------------------
     ! Read dm particle data
     !----------------------------------------------
     nomfich=TRIM(repository)//'/halo_'//TRIM(ncharlong)//'/dmpart.dat'
     INQUIRE(file=nomfich,exist=ok)
     if(.not.ok)then
        write(*,*)'Skipping file ',TRIM(nomfich)
     else
        write(*,*)'Reading file ',TRIM(nomfich)
        open(unit=10,file=nomfich,status='old',form='unformatted')
        read(10)idhalo2(ihalo),subhalo(ihalo)
        if(idhalo2(ihalo).ne.idhalo(ihalo))then
           write(*,'(A,I10,I10)')'WARNING: Halo id in the file is different from the name file',idhalo2(ihalo),idhalo(ihalo)
           write(*,*)'Something suspicious is happening! I continue nonetheless...'
        endif
        read(10)mhalo(ihalo)
        read(10)xhalo(ihalo)
        read(10)yhalo(ihalo)
        read(10)zhalo(ihalo)
        read(10)rvir(ihalo),rhalo(ihalo)
        read(10)aexp
        read(10)scale_l,scale_d,scale_t
        read(10)ndm,nread
        write(*,*)'ndm  =',ndm
        allocate(xdm(1:ndm,1:3),vdm(1:ndm,1:3),mdm(1:ndm),idm(1:ndm),rdm(1:ndm))
        read(10)idm(1:ndm)
        read(10)xdm(1:ndm,1)
        read(10)xdm(1:ndm,2)
        read(10)xdm(1:ndm,3)
        read(10)vdm(1:ndm,1)
        read(10)vdm(1:ndm,2)
        read(10)vdm(1:ndm,3)
        read(10)mdm(1:ndm)
        close(10)

        do i=1,ndm
           xread=xdm(i,1)
           yread=xdm(i,2)
           zread=xdm(i,3)
           ! Check for periodicity
           dxx=xhalo(ihalo)-xread
           if(dxx> 0.5)then
              dxx=dxx-1d0
           endif
           if(dxx<-0.5)then
              dxx=dxx+1d0
           endif
           dr2=dxx*dxx
           dyy=yhalo(ihalo)-yread
           if(dyy> 0.5)then
              dyy=dyy-1d0
           endif
           if(dyy<-0.5)then
              dyy=dyy+1d0
           endif
           dr2=dyy*dyy+dr2
           dzz=zhalo(ihalo)-zread
           if(dzz> 0.5)then
              dzz=dzz-1d0
           endif
           if(dzz<-0.5)then
              dzz=dzz+1d0
           endif
           dr2=dzz*dzz+dr2
           rdm(i)=sqrt( dr2 )
        enddo
        msun=scale_d*scale_l**3d0/2d33
        write(*,'(A,ES13.6)')'Mdm  =',sum(mdm(1:ndm)*msun)
        !-----------------------------------------------
        ! End read dm particle data
        !-----------------------------------------------
     endif

     !-----------------------------------------------
     ! Read star particle data
     !-----------------------------------------------
     nomfich=TRIM(repository)//'/halo_'//TRIM(ncharlong)//'/starpart.dat'
     INQUIRE(file=nomfich,exist=ok)
     if(.not.ok)then
        write(*,*)'Skipping file ',TRIM(nomfich)
     else
        write(*,*)'Reading file ',TRIM(nomfich)

        ! First pass to get global info and number of cells
        open(unit=10,file=nomfich,status='old',form='unformatted')
        read(10)idhalo2(ihalo),subhalo(ihalo)
        if(idhalo2(ihalo).ne.idhalo(ihalo))then
           write(*,'(A,I10,I10)')'WARNING: Halo id in the file is different from the name file',idhalo2(ihalo),idhalo(ihalo)
           write(*,*)'Something suspicious is happening! I continue nonetheless...'
        endif
        read(10)mhalo(ihalo)
        read(10)xhalo(ihalo)
        read(10)yhalo(ihalo)
        read(10)zhalo(ihalo)
        read(10)rvir(ihalo),rhalo(ihalo)
        read(10)aexp
        read(10)scale_l,scale_d,scale_t
        read(10)nvarh
        read(10)nstar,nread
        write(*,*)'nstar=',nstar
        nchemo=nvarh-6
        write(*,*)'chemo=',chemo
        allocate(xstar(1:nstar,1:3),vstar(1:nstar,1:3),mstar(1:nstar),istar(1:nstar),tstar(1:nstar),rstar(1:nstar))
        if(nvarh.gt.5)allocate(zstar(1:nstar))
        if(nchemo.gt.0.and.chemo)allocate(chstar(1:nstar,1:nchemo),chread(1:nchemo))
        read(10)istar(1:nstar)
        read(10)xstar(1:nstar,1)
        read(10)xstar(1:nstar,2)
        read(10)xstar(1:nstar,3)
        read(10)vstar(1:nstar,1)
        read(10)vstar(1:nstar,2)
        read(10)vstar(1:nstar,3)
        read(10)mstar(1:nstar)
        read(10)tstar(1:nstar)
        if(nvarh.gt.5)then
           read(10)zstar(1:nstar)
           if(nchemo.gt.0.and.chemo)then
              do i=1,nchemo
                 read(10)chstar(1:nstar,i)
              enddo
           endif
        endif
        close(10)

        do i=1,nstar
           xread=xstar(i,1)
           yread=xstar(i,2)
           zread=xstar(i,3)
           ! Check for periodicity
           dxx=xhalo(ihalo)-xread
           if(dxx> 0.5)then
              dxx=dxx-1d0
           endif
           if(dxx<-0.5)then
              dxx=dxx+1d0
           endif
           dr2=dxx*dxx
           dyy=yhalo(ihalo)-yread
           if(dyy> 0.5)then
              dyy=dyy-1d0
           endif
           if(dyy<-0.5)then
              dyy=dyy+1d0
           endif
           dr2=dyy*dyy+dr2
           dzz=zhalo(ihalo)-zread
           if(dzz> 0.5)then
              dzz=dzz-1d0
           endif
           if(dzz<-0.5)then
              dzz=dzz+1d0
           endif
           dr2=dzz*dzz+dr2
           rstar(i)=sqrt( dr2 )
        enddo
        msun=scale_d*scale_l**3d0/2d33
        write(*,'(A,ES13.6)')'Mstar=',sum(mstar(1:nstar)*msun)
        if(nchemo.gt.0.and.chemo)deallocate(chread)
        !-----------------------------------------------
        ! End read star particle data
        !-----------------------------------------------
     endif

     ! Find r200c
     drbin=rhalo(ihalo)/dble(nbin)
     ! ---- GAS ----
     allocate(dmcell(1:nbin))
     dmcell=0d0
     do i=1,ncell
        index=int(rcell(i)/drbin)+1
        if(index .le. nbin)then
           dmcell(index)=dmcell(index)+var(i,1)*dxcell(i)**3d0
        endif
     enddo
     do i=2,nbin
        dmcell(i)=dmcell(i)+dmcell(i-1)
     enddo
     !write(*,*)dmcell*msun
     ! ---- GAS ----
     ! ---- DM -----
     allocate(dmdm(1:nbin))
     dmdm=0d0
     do i=1,ndm
        index=int(rdm(i)/drbin)+1
        if(index .le. nbin)then
           dmdm(index)=dmdm(index)+mdm(i)
        endif
     enddo
     do i=2,nbin
        dmdm(i)=dmdm(i)+dmdm(i-1)
     enddo
     ! ---- DM -----
     !write(*,*)dmdm*msun/1d13
     ! ---- STAR -----
     allocate(dmstar(1:nbin))
     dmstar=0d0
     do i=1,nstar
        index=int(rstar(i)/drbin)+1
        if(index .le. nbin)then
           dmstar(index)=dmstar(index)+mstar(i)
        endif
     enddo
     do i=2,nbin
        dmstar(i)=dmstar(i)+dmstar(i-1)
     enddo
     !write(*,*)dmstar*msun/1d13
     ! ---- STAR -----
     allocate(dmtot(1:nbin),radius(1:nbin),density(1:nbin))
     dmtot=dmcell+dmdm+dmstar
     do i=1,nbin
        radius (i)=dble(i)*drbin
        density(i)=dmtot(i)/(4d0*acos(-1d0)/3d0*radius(i)**3d0)
     enddo
     !write(*,*)dmtot*msun
     !write(*,*)density*msun/(scale_l/3.08d21)**3d0
     rhoc=scale_d/omega_m !rho_crit at the corresponding redshift

     open(unit=1,file='test.dat')
     do i=1,nbin
        write(1,'(2es14.7)')radius(i)*scale_l/3.08d24,density(i)*scale_d/rhoc
     enddo
     close(1)

     write(99)radius*scale_l/3.08d24,dmcell*msun,dmdm*msun,dmstar*msun

     deallocate(xcell)
     deallocate(dxcell)
     deallocate(var)
     deallocate(rcell)
     deallocate(varread)
     deallocate(xdm,vdm,mdm,idm,rdm)
     deallocate(xstar,vstar,mstar,istar,tstar,rstar)
     if(nvarh.gt.5)deallocate(zstar)
     if(nchemo.gt.0.and.chemo)deallocate(chstar)
     deallocate(dmcell,dmdm,dmstar,dmtot,radius,density)

  enddo

contains

  !=======================================================================
  subroutine read_params

    implicit none

    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok

    n = iargc()
    if (n < 4) then
       print *, 'usage: read_halodata   -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 -lis  halo list'
       print *, 'ex: read_halodata -inp HaloData -lis list_halo_00001.dat -out haloprops.dat'
       print *, 'ex: read_halodata -inp HaloData -idh 1 -out haloprops.dat'
       print *, ' '
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
       case ('-lis')
          listhalo = trim(arg)
       case ('-idh')
          read (arg,*) ihread
       case ('-om')
          read (arg,*) omega_m
       case ('-nbi')
          read (arg,*) nbin
       case ('-che')
          read (arg,*) chemo
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

    return

  end subroutine read_params
  !=======================================================================
  !=======================================================================
  subroutine title2(n,nchar)
    !=======================================================================
    implicit none
    integer::n
    character*7::nchar

    character*1::nchar1
    character*2::nchar2
    character*3::nchar3
    character*4::nchar4
    character*5::nchar5
    character*6::nchar6
    character*7::nchar7

    if(n.ge.1000000)then
       write(nchar7,'(i7)') n
       nchar = nchar7
    elseif(n.ge.100000)then
       write(nchar6,'(i6)') n
       nchar = '0'//nchar6
    elseif(n.ge.10000)then
       write(nchar5,'(i5)') n
       nchar = '00'//nchar5
    elseif(n.ge.1000)then
       write(nchar4,'(i4)') n
       nchar = '000'//nchar4
    elseif(n.ge.100)then
       write(nchar3,'(i3)') n
       nchar = '0000'//nchar3
    elseif(n.ge.10)then
       write(nchar2,'(i2)') n
       nchar = '00000'//nchar2
    else
       write(nchar1,'(i1)') n
       nchar = '000000'//nchar1
    endif

  end subroutine title2
  !================================================================

end program read_halodata
