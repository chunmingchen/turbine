!     ==================================================================
      subroutine out_probe_reading
!     ==================================================================
      use common_area
      use variable_area
      use error_report
      implicit none
!this routine does not support k-partision
!blocks from low_limit to up_limit form the probe ring
!----------------------------------------------------
      integer::num_probe_rings=10
!---------------------------------------------------
      integer,save,allocatable :: my_pring_comm(:)
      type probe_ring
           integer ::low_limit,up_limit,nkbr,num_probes
           integer ::iindex,jindex,bound1,bound2
           integer ::pring_comm
           real*8,pointer ::probe_phases(:)
           real*8,pointer ::xring_i(:)
           real*8,pointer ::xring_iw(:)
           real*8,pointer ::xring_ang(:)
      end type
      type(probe_ring),save,allocatable :: pring(:)

      real*8,allocatable,dimension(:),save::pring_ang
      real*8,allocatable,dimension(:),save::pring_q1
      real*8,allocatable,dimension(:),save::pring_q2
      real*8,allocatable,dimension(:),save::pring_q3
      real*8,allocatable,dimension(:),save::pring_q4
      real*8,allocatable,dimension(:),save::pring_q5
      real*8,allocatable,dimension(:),save::pring_q6
      real*8,allocatable,dimension(:),save::pring_s
      real*8,allocatable,dimension(:),save::r_probe
      real*8,allocatable,dimension(:),save::ru_probe
      real*8,allocatable,dimension(:),save::rvr_probe
      real*8,allocatable,dimension(:),save::rvt_probe
      real*8,allocatable,dimension(:),save::e_probe
      real*8,allocatable,dimension(:),save::p_probe
      real*8,allocatable,dimension(:),save::s_probe
      real*8,allocatable,dimension(:),save::thta

      real*8 dthta,angle,yc,zc,rvcar,rwcar,theta_tmp
      real*8 sinthta,costhta,den,w1,w2,thta1,thta2
      real*8 y1,y2,z1,z2,xp,x1,x2,sector_ang
      real*8 angle1,angle2,thta_shift,thta_target
      real*8,allocatable,save :: buf_tmp_l(:),buf_tmp(:)

      integer,save,allocatable ::my_pring_group(:)
      integer,save,allocatable ::lu_probe_p_hs(:)
      integer,save,allocatable ::lu_probe_s_hs(:)
      integer,save,allocatable ::lu_probe_u_hs(:)
      integer,save,allocatable ::lu_probe_v_hs(:)
      integer,save,allocatable ::lu_probe_w_hs(:)
      integer ierr,ll,i,j,k,l,m,n,ntry,nkbrp1,count
      integer kk,kph1,kph2,iend1,iend2,iend3,i1,i2
      integer bound1,bound2,id_passage,max_nkbr,max_probes
      integer :: data_size
      logical,save :: first_time=.true.
      character :: xch1*13,xch2*20,xch3*3,fname*40

      if(first_time) then
      allocate(pring(num_probe_rings))

! specify probe rings, and probes
!-----------------------------------------
      pring(1)%low_limit=37   !full annulus  starting block # for rotor. 
      pring(1)%up_limit=72    !full annulus  ending block # for rotor
c     pring(1)%low_limit=19   !half annulus
c     pring(1)%up_limit=36    !half annulus
      pring(1)%nkbr=1980      !(nk-1)*nbld(nbr) total # of grid in circumferential direction
      pring(1)%num_probes=8
      pring(1)%iindex=9  !43.7% chord ahead leading edge
!     pring(1)%jindex=70
      pring(1)%jindex=57   !98% span
      pring(1)%bound1=pring(1)%low_limit-1
      pring(1)%bound2=pring(1)%up_limit-1
      ll=pring(1)%num_probes
      allocate(pring(1)%probe_phases(ll))
      pring(1)%probe_phases(1)= 10./360.*pi2
      pring(1)%probe_phases(2)= 70./360.*pi2
      pring(1)%probe_phases(3)=100./360.*pi2
      pring(1)%probe_phases(4)=160./360.*pi2
      pring(1)%probe_phases(5)=190./360.*pi2
      pring(1)%probe_phases(6)=250./360.*pi2
      pring(1)%probe_phases(7)=280./360.*pi2
      pring(1)%probe_phases(8)=340./360.*pi2

      pring(2)%low_limit=37   !full annulus
      pring(2)%up_limit=72    !full annulus
c     pring(2)%low_limit=19   !half annulus
c     pring(2)%up_limit=36    !half annulus
      pring(2)%nkbr=1980      !(nk-1)*nbld(nbr)
      pring(2)%num_probes=8
      pring(2)%iindex=19 !15% chord ahead leading edge
!     pring(2)%jindex=70
      pring(2)%jindex=57   !98% span
      pring(2)%bound1=pring(2)%low_limit-1
      pring(2)%bound2=pring(2)%up_limit-1
      ll=pring(2)%num_probes
      allocate(pring(2)%probe_phases(ll))
      pring(2)%probe_phases(1)= 10./360.*pi2
      pring(2)%probe_phases(2)= 70./360.*pi2
      pring(2)%probe_phases(3)=100./360.*pi2
      pring(2)%probe_phases(4)=160./360.*pi2
      pring(2)%probe_phases(5)=190./360.*pi2
      pring(2)%probe_phases(6)=250./360.*pi2
      pring(2)%probe_phases(7)=280./360.*pi2
      pring(2)%probe_phases(8)=340./360.*pi2

      pring(3)%low_limit=109   !full annulus
      pring(3)%up_limit=144    !full annulus
c     pring(3)%low_limit=55    !half annulus
c     pring(3)%up_limit=72     !half annulus
      pring(3)%nkbr=1980
      pring(3)%num_probes=4
      pring(3)%iindex=41   !10% chord downstream rotor
      pring(3)%jindex=57   !98% span
      pring(3)%bound1=pring(3)%low_limit-1
      pring(3)%bound2=pring(3)%up_limit-1
      ll=pring(3)%num_probes
      allocate(pring(3)%probe_phases(ll))
      pring(3)%probe_phases(1)=  0./360.*pi2
      pring(3)%probe_phases(2)= 90./360.*pi2
      pring(3)%probe_phases(3)=180./360.*pi2
      pring(3)%probe_phases(4)=270./360.*pi2

      do n=4, num_probe_rings
        pring(n)%low_limit=109   !full annulus
        pring(n)%up_limit=144    !full annulus
        pring(n)%nkbr=1980
        pring(n)%num_probes=36
        pring(n)%iindex=(n-4)*10+10   !10% chord downstream rotor
        pring(n)%jindex=57   !98% span
        pring(n)%bound1=pring(n)%low_limit-1
        pring(n)%bound2=pring(n)%up_limit-1
        ll=pring(n)%num_probes
        allocate(pring(n)%probe_phases(ll))
        do ll=0,35
            pring(n)%probe_phases(ll+1)=ll*10./360.*pi2
        enddo
      enddo
!----------------------------------------------
      max_nkbr=0
      max_probes=0
      do n=1,num_probe_rings
         ll=pring(n)%nkbr
         max_nkbr=max(max_nkbr,ll)
         ll=pring(n)%num_probes
         max_probes=max(max_probes,ll)     
      enddo
      nkbrp1=max_nkbr+1
      data_size=max_nkbr*7
      allocate(pring_ang(max_nkbr))
      allocate(pring_q1(max_nkbr))
      allocate(pring_q2(max_nkbr))
      allocate(pring_q3(max_nkbr))
      allocate(pring_q4(max_nkbr))
      allocate(pring_q5(max_nkbr))
      allocate(pring_q6(max_nkbr))
      allocate(pring_s(max_nkbr))
      allocate(buf_tmp_l(data_size))
      allocate(buf_tmp(data_size))  
      allocate(thta(nkbrp1))
      allocate(r_probe(max_probes))
      allocate(ru_probe(max_probes))
      allocate(rvr_probe(max_probes))
      allocate(rvt_probe(max_probes))
      allocate(e_probe(max_probes))
      allocate(p_probe(max_probes))
      allocate(s_probe(max_probes))
!--------------------------------------

c     if(first_time) then

      allocate(my_pring_comm(num_probe_rings))
      allocate(my_pring_group(num_probe_rings))
      allocate(lu_probe_p_hs(num_probe_rings))
      allocate(lu_probe_s_hs(num_probe_rings))
      allocate(lu_probe_u_hs(num_probe_rings))
      allocate(lu_probe_v_hs(num_probe_rings))
      allocate(lu_probe_w_hs(num_probe_rings))

      write(lu_out,*) 'x for probe: ',x(9,57,2)
      do 10 n=1,num_probe_rings
         my_pring_group(n)=MPI_UNDEFINED
         bound1=pring(n)%bound1
         bound2=pring(n)%bound2
         if(my_id.le.bound2 .and. my_id.ge.bound1)then
                my_pring_group(n)=n
         endif 
         call mpi_comm_split(mpi_comm_world,my_pring_group(n),
     &                       my_id,my_pring_comm(n),ierr)
         
! creat output file for each probe ring
         call int_to_char(n,xch2,ll)
         xch3='.hs'

         xch1='probe_p_ring.'
         fname=xch1//xch2(1:ll)//xch3
         if(my_id==bound1)
     &        lu_probe_p_hs(n)=
     &        fopen(fname,'formatted',position='append')

         xch1='probe_s_ring.'
         fname=xch1//xch2(1:ll)//xch3
         if(my_id==bound1)
     &        lu_probe_s_hs(n)=
     &        fopen(fname,'formatted',position='append')

         xch1='probe_u_ring.'
         fname=xch1//xch2(1:ll)//xch3
         if(my_id==bound1)
     &        lu_probe_u_hs(n)=
     &        fopen(fname,'formatted',position='append')

         xch1='probe_v_ring.'
         fname=xch1//xch2(1:ll)//xch3
         if(my_id==bound1)
     &        lu_probe_v_hs(n)=
     &        fopen(fname,'formatted',position='append')

         xch1='probe_w_ring.'
         fname=xch1//xch2(1:ll)//xch3
         if(my_id==bound1)
     &        lu_probe_w_hs(n)=
     &        fopen(fname,'formatted',position='append')



         if(my_pring_group(n)/=n)goto 10 

! at the first time, creat probe ring positioned at x=const by 
!    interpolation along i direction.
        allocate(pring(n)%xring_i(nk))
        allocate(pring(n)%xring_iw(nk))
        allocate(pring(n)%xring_ang(nk))

        pring(n)%xring_ang=0.
        i=pring(n)%iindex
        j=pring(n)%jindex
        xp=x(i,j,1)
          do 20 k=2,nk
             iend1=2
             iend3=ni
             ntry=0
  50         ntry=ntry+1
             if(ntry.gt.20) then
               write(error_message(1),*) 
     &         'Search in I-direction failed for nbr ',nbr
               call fatal_error_trap('out_probe_reading',1)
             endif

          if(iend3.eq.iend1+1) then
             i1=iend1
             i2=iend3
            x1=0.125*(x(i1,j,  k)  +x(i1-1,j,  k)
     .               +x(i1,j-1,k)  +x(i1-1,j-1,k)
     .               +x(i1,j,  k-1)+x(i1-1,j,  k-1)
     .               +x(i1,j-1,k-1)+x(i1-1,j-1,k-1))
            x2=0.125*(x(i2,j,  k)  +x(i2-1,j,  k)
     .               +x(i2,j-1,k)  +x(i2-1,j-1,k)
     .               +x(i2,j,  k-1)+x(i2-1,j,  k-1)
     .               +x(i2,j-1,k-1)+x(i2-1,j-1,k-1))
             goto 60
          endif

            iend2=iend1+(iend3-iend1)/2
            x1=0.125*(x(iend1,j,  k)  +x(iend1-1,j,  k)
     .               +x(iend1,j-1,k)  +x(iend1-1,j-1,k)
     .               +x(iend1,j,  k-1)+x(iend1-1,j,  k-1)
     .               +x(iend1,j-1,k-1)+x(iend1-1,j-1,k-1))
            x2=0.125*(x(iend2,j,  k)  +x(iend2-1,j,  k)
     .               +x(iend2,j-1,k)  +x(iend2-1,j-1,k)
     .               +x(iend2,j,  k-1)+x(iend2-1,j,  k-1)
     .               +x(iend2,j-1,k-1)+x(iend2-1,j-1,k-1))
            if(xp>x1 .and. xp<x2) then
              iend3=iend2
            else
              iend1=iend2
            endif
            goto 50

   60       pring(n)%xring_i(k)=i1
            w1=(xp-x1)/(x2-x1)
            pring(n)%xring_iw(k)=w1

            i=i1
            y1=0.125*(y(i,j,k)+y(i,j,k-1)+y(i-1,j,k-1)+y(i-1,j,k)+
     &         y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k))
            z1=0.125*(z(i,j,k)+z(i,j,k-1)+z(i-1,j,k-1)+z(i-1,j,k)+
     &         z(i,j-1,k)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j-1,k))

            i=i2
            y2=0.125*(y(i,j,k)+y(i,j,k-1)+y(i-1,j,k-1)+y(i-1,j,k)+
     &         y(i,j-1,k)+y(i,j-1,k-1)+y(i-1,j-1,k-1)+y(i-1,j-1,k))
            z2=0.125*(z(i,j,k)+z(i,j,k-1)+z(i-1,j,k-1)+z(i-1,j,k)+
     &         z(i,j-1,k)+z(i,j-1,k-1)+z(i-1,j-1,k-1)+z(i-1,j-1,k))
            
            yc=(y2-y1)*w1+y1
            zc=(z2-z1)*w1+z1
            pring(n)%xring_ang(k)=atan2(yc,zc)
 20     continue  !k=2,nk
 10     continue  !n=1,num_probe_rings
        first_time=.false.
      endif  ! (first_time)

      do 120 n=1,num_probe_rings

      call mpi_barrier(mpi_comm_world,ierr)
      if(my_pring_group(n)/=n)goto 120

            buf_tmp_l=0.
            buf_tmp=0.
            thta=0.
            pring_ang=0.
            pring_q1=0.
            pring_q2=0.
            pring_q3=0.
            pring_q4=0.
            pring_q5=0.
            pring_q6=0.
            pring_s=0.

      nkbrp1=pring(n)%nkbr+1
      dthta=dtdt(nbr)*time
      m=nbld(nbr)/ksym_br(nbr)
      sector_ang=pi2/ksym_br(nbr)
      j=pring(n)%jindex
      bound1=pring(n)%bound1
      bound2=pring(n)%bound2
          do l=1,ksym_br(nbr)  
          do k=2,nk
!do not support k-partision on the next line
!assumed each block has same nk values 
            id_passage=(l-1)*m+(my_id-bound1)
            kk=id_passage*(nk-1)+k-1 
            w1=pring(n)%xring_iw(k)
            w2=1.-w1
            i1=pring(n)%xring_i(k)
            i2=i1+1
            angle=pring(n)%xring_ang(k)
            sinthta=sin(angle)
            costhta=cos(angle)
            rvcar=w1*q(3,i2,j,k)+w2*q(3,i1,j,k)
            rwcar=w1*q(4,i2,j,k)+w2*q(4,i1,j,k)
            theta_tmp=angle+dthta-(l-1)*sector_ang
            theta_tmp=theta_tmp+1000.*pi2
            pring_ang(kk)=mod(theta_tmp,pi2)
            pring_q1(kk)=w1*q(1,i2,j,k)+w2*q(1,i1,j,k)
            pring_q2(kk)=w1*q(2,i2,j,k)+w2*q(2,i1,j,k)
            pring_q3(kk)=rvcar*sinthta+rwcar*costhta
            pring_q4(kk)=rvcar*costhta-rwcar*sinthta
            pring_q5(kk)=w1*q(5,i2,j,k)+w2*q(5,i1,j,k)
            pring_q6(kk)=w1*p(i2,j,k)+w2*p(i1,j,k)
            pring_s(kk) =w1*p(i2,j,k)/q(1,i2,j,k)**gam+
     &                   w2*p(i1,j,k)/q(1,i1,j,k)**gam
          enddo
          enddo
!// pack probe ring data
          count=0
          do k=1,pring(n)%nkbr
             count=count+1 
             buf_tmp_l(count)=pring_ang(k)
             count=count+1
             buf_tmp_l(count)=pring_q1(k)
             count=count+1
             buf_tmp_l(count)=pring_q2(k)
             count=count+1
             buf_tmp_l(count)=pring_q3(k)
             count=count+1
             buf_tmp_l(count)=pring_q4(k)
             count=count+1
             buf_tmp_l(count)=pring_q5(k)
             count=count+1
             buf_tmp_l(count)=pring_q6(k)
             count=count+1
             buf_tmp_l(count)=pring_s(k)
          enddo

          buf_tmp=0.0
       call mpi_reduce(buf_tmp_l,buf_tmp,count,mpi_real8,
     &                mpi_sum,0,my_pring_comm(n),ierr)

          if(my_id /=bound1) goto 120
!// unpack probe ring data
          count=0
          do k=1,pring(n)%nkbr
             count=count+1
             pring_ang(k)=buf_tmp(count)
             count=count+1
             pring_q1(k)=buf_tmp(count)
             count=count+1
             pring_q2(k)=buf_tmp(count)
             count=count+1
             pring_q3(k)=buf_tmp(count)
             count=count+1
             pring_q4(k)=buf_tmp(count)
             count=count+1
             pring_q5(k)=buf_tmp(count)
             count=count+1
             pring_q6(k)=buf_tmp(count)
             count=count+1
             pring_s(k) =buf_tmp(count)
          enddo
!//extract probe readings from probe ring data
       do 200 l=1,pring(n)%num_probes
       thta_shift=pring_ang(1)
       thta_target=pring(n)%probe_phases(l)+
     .             1000.*pi2-thta_shift
       thta_target=mod(thta_target,pi2)
       thta_target=pi2-thta_target

       do k=2,pring(n)%nkbr
          thta(k)=pring_ang(k)+1000*pi2-thta_shift
          thta(k)=mod(thta(k),pi2)
          thta(k)=pi2-thta(k)
       enddo
       thta(1)=0.0
       thta(nkbrp1)=pi2

      iend1=1
      iend3=nkbrp1
      ntry=0
  100 ntry=ntry+1
      if(ntry.gt.20) then
      write(error_message(1),*) 
     &  'Search in I-buffer ring failed for nbr ',nbr
      call fatal_error_trap('out_probe_reading',1)
      endif

      if(iend3.eq.iend1+1) then
        kph1=iend1
        kph2=iend3
        if(kph2.eq.nkbrp1) then
          kph2=1
          thta1=thta(kph1)
          thta2=pi2
        else
          thta1=thta(kph1)
          thta2=thta(kph2)
        endif
        goto 150
      endif

      iend2=iend1+(iend3-iend1)/2
      thta1=thta(iend1)
      thta2=thta(iend2)
      if(thta_target.ge.thta1.and.thta_target.le.thta2) then
        iend3=iend2
      else
        iend1=iend2
      endif
      goto 100

  150 w1 = thta_target-thta1
      w2 = thta2-thta_target
      den = w1 + w2
      w1 = w1/den
      w2 = w2/den
         r_probe(l)=w1*pring_q1(kph2)+w2*pring_q1(kph1)
        ru_probe(l)=w1*pring_q2(kph2)+w2*pring_q2(kph1)
       rvr_probe(l)=w1*pring_q3(kph2)+w2*pring_q3(kph1)
       rvt_probe(l)=w1*pring_q4(kph2)+w2*pring_q4(kph1)
         e_probe(l)=w1*pring_q5(kph2)+w2*pring_q5(kph1)
         p_probe(l)=w1*pring_q6(kph2)+w2*pring_q6(kph1)
         s_probe(l)=w1*pring_s( kph2)+w2*pring_s( kph1)
 200  continue  ! l=1,pring(n)%num_probes

      if(my_id == bound1)then
          write(lu_probe_p_hs(n),400)ncyct,
     &         (p_probe(l)+(l-1),l=1,pring(n)%num_probes)

          write(lu_probe_s_hs(n),400)ncyct,
     &         (s_probe(l)+(l-1),l=1,pring(n)%num_probes)

          write(lu_probe_u_hs(n),400)ncyct,
     &         ( ru_probe(l)/r_probe(l)+(l-1),
     &           l=1,pring(n)%num_probes)

          write(lu_probe_v_hs(n),400)ncyct,
     &         (rvr_probe(l)/r_probe(l)+(l-1),
     &           l=1,pring(n)%num_probes)

          write(lu_probe_w_hs(n),400)ncyct,
     &         (rvt_probe(l)/r_probe(l)+(l-1),
     &           l=1,pring(n)%num_probes)

 400      format(i8,20f8.5)
          call dummy_flush(lu_probe_p_hs(n))
          call dummy_flush(lu_probe_s_hs(n))
          call dummy_flush(lu_probe_u_hs(n))
          call dummy_flush(lu_probe_v_hs(n))
          call dummy_flush(lu_probe_w_hs(n))
      endif

 120  continue !  n=1, num_probe_rings

c     deallocate(pring_q1)
c     deallocate(pring_q2)
c     deallocate(pring_q3)
c     deallocate(pring_q4)
c     deallocate(pring_q5)
c     deallocate(pring_q6)
c     deallocate(pring_6)
c     deallocate(buf_tmp_l)
c     deallocate(buf_tmp)
c     deallocate(thta)

      return
      end
!     ==================================================================
      subroutine int_to_char(num,xch,l1)
!     read inital condition namelist
!     ==================================================================
      integer,intent(in) :: num
      integer:: l1,inum,iord
      character,intent(out):: xch*20

       l1=0
       inum=num
       do iord=1,10
         inum=inum/10
         if(inum.gt.0)l1=l1+1
       enddo
       l1=l1+1
       if(l1.eq.10)write(xch,'(i10)') num
       if(l1.eq.9)write(xch,'(i9)') num
       if(l1.eq.8)write(xch,'(i8)') num
       if(l1.eq.7)write(xch,'(i7)') num
       if(l1.eq.6)write(xch,'(i6)') num
       if(l1.eq.5)write(xch,'(i5)') num
       if(l1.eq.4)write(xch,'(i4)') num
       if(l1.eq.3)write(xch,'(i3)') num
       if(l1.eq.2)write(xch,'(i2)') num
       if(l1.eq.1)write(xch,'(i1)') num
      return
      end subroutine int_to_char

