      subroutine fci_uhf_string(rtdb,
     &    d_v2,k_v2_offset,nos,nvs,noas,nobs,nvas,nvbs,
     &    nstot,
     &    h,v,repulsion)
c
c For fci_uhf  never use the intorb i.e. 
c 2emet /2eorb options in the input.
c Always use spinorbitals form of 4-index transformation
c
c all varaible below have to be defined in the tce_energy.f
c
c please do not freeze the occupied orbitals for now
c
c nstot - total number of spinorbitals nstot=nos+nvs
c nos   - number of occupied spinorbitals
c nvs   - number of virtual spinorbitals
c noas  - number of occupied alpha spinorbitals
c nobs  - number of occupied beta spinorbitals
c nvas  - number of virtual  alhoa spinorbitals
c nvbs  - number of virtual beta spinorbitals
c
c spinorbital convention
c
c | noas | nobs | nvas | nvbs |
c
c ene_orb contains orbital energies
c
c COMMENTS: 
c iext1 & iext2 - redundant
c
c
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_v2    ! GA handle for v2
      integer k_v2_offset
      integer nos,nvs ! # of occupied/virtual spinorbitals
      integer noas,nobs,nvas,nvbs ! # of occupied/virtual alphas betas
      integer nstot
cc      double precision ene_orb(nos+nvs)  ! a copy of dbl_mb(k_sorted)
      integer size
      integer rtdb
      integer l_aux1,k_aux1 ! for local memory allocator: loc. mem. buffer 1
      integer i,j,k,l,m,n,p,q,r,s         ! auxiliary indices
      integer ia,ib,ic,id,ie      !auxiliary indices
      integer iaux
c      integer nact              !number of active virtual orbitals
c h created here from "orbital" matrix horb
      double precision h(nos+nvs,nos+nvs)
c horb valid only for RHF case
      double precision horb((nos+nvs)/2,(nos+nvs)/2)
      double precision v(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
c FCI dimension
      integer qnumber
      integer irrep_g
      integer dim_fci
      integer alpha_str
      integer beta_str
      integer l_str_a,k_str_a
      integer l_str_b,k_str_b
      integer k_matrix,l_matrix
      logical sym_test
c
c
c
      double precision xxx,yyy,zzz
c
      double precision repulsion
c
      logical nodezero
      logical oprint_qa
c
c
      nodezero = (ga_nodeid().eq.0)
      oprint_qa=util_print('ducc_qa', print_high)
c
c
c 
      if(nodezero) then
          write(6,*)'-------- string-based FCI ---------- '
          call util_flush(6)
      endif
c
c
c

      do i=1,(nos+nvs)/2
      do j=1,(nos+nvs)/2
        horb(i,j)=0.0d0
      enddo
      enddo
c
      do i=1,nos+nvs
      do j=1,nos+nvs
        h(i,j)=0.0d0
      enddo
      enddo


c 
c
c########################################################################
c     call kinetic_hcore_1spin(rtdb,horb,h,nos,nvs,noas,nobs,nvas,nvbs)
c########################################################################
c
c     
c forming v matrix
c
      do i=1,nos+nvs
      do j=1,nos+nvs
      do k=1,nos+nvs
      do l=1,nos+nvs
        v(i,j,k,l)=0.0d0 
      enddo
      enddo
      enddo
      enddo
c
c For fci_uhf  never use the intorb i.e. 
c 2emet /2eorb options in the input
c`
      if(nodezero) then
          write(6,*)'step2 '
          call util_flush(6)
      endif
c
c########################################################################
c     call mapping_v2_spin(rtdb,d_v2,k_v2_offset,v,nos,nvs)
c########################################################################
c
c
c
c  you have everything here for FCI: v & h are ready
c
c  ATTENTION:   Works only for the RHF case
c
c
c
c  Matrix dimensions - simplest case Ms symmetry binom(n,k)*binom(n,k) 
c  n=noas+nvas / k=noas
c
c  
        irrep_g = 0
        dim_fci = 0
        alpha_str = 0
        beta_str = 0
        sym_test = .true.
c
      if (.not.MA_PUSH_GET(mt_int,(noas+nvas),'a_string',l_str_a,
     1 k_str_a)) CALL ERRQUIT('a_string',0,MA_ERR)
c
        qnumber = 0
        do m = 0, noas+nvas-1
         qnumber = qnumber + 2**m
        enddo 
c
        do i=1,qnumber   !main loop 
c
         do m=1,noas+nvas
          int_mb(k_str_a+m-1) = 0
         enddo
c         
         k=i
c
         do j = noas+nvas, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str_a+j-1) = 0
            else
             int_mb(k_str_a+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of alpha particle check
         m=0
         do j=1,noas+nvas
c additional symmetry checks possible here
          m = m+int_mb(k_str_a+j-1)
         enddo
         if(m.eq.noas) alpha_str = alpha_str + 1
c
        enddo !main loop
c
      if (.not.MA_POP_STACK(l_str_a))
     1   call errquit("l_str",0,MA_ERR)
c
c
c
      if (.not.MA_PUSH_GET(mt_int,(nobs+nvbs),'b_string',l_str_b,
     1 k_str_b)) CALL ERRQUIT('a_string',0,MA_ERR)
c
        qnumber = 0
        do m = 0, nobs+nvbs-1
         qnumber = qnumber + 2**m
        enddo
c
        do i=1,qnumber   !main loop 
c
         do m=1,nobs+nvbs
          int_mb(k_str_b+m-1) = 0
         enddo
c         
         k=i
c
         do j = nobs+nvbs, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str_b+j-1) = 0
            else
             int_mb(k_str_b+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of alpha particle check
         m=0
         do j=1,noas+nvas
c additional symmetry checks possible here
          m = m+int_mb(k_str_b+j-1)
         enddo
         if(m.eq.nobs) beta_str = beta_str + 1
c
        enddo !main loop
c
      if (.not.MA_POP_STACK(l_str_b))
     1   call errquit("l_str",0,MA_ERR)
c
c
c
c
c
c
c
c
c
c
c
c
c        Determination of dim_fci with symmetry
c        returns dim_fci (with symmetry inclusion)
         call dim_det_uhf(noas,nobs,nvas,nvbs,nos,nvs,nstot,
     1        irrep_g,alpha_str,beta_str,dim_fci)
c
c
c
      if(nodezero) then
        write(6,*)'dim_fci: ',dim_fci
        call util_flush
      endif
c
      if (.not.MA_PUSH_GET(mt_dbl,dim_fci*dim_fci,'fcim',
     1      l_matrix,k_matrix))
     1      call errquit('fcim matrix',1,MA_ERR)
c
c 
c initialization
c
       do i=1,dim_fci*dim_fci
         dbl_mb(k_matrix+i-1) = 0.0d0
       enddo
c
c
c
c
       call fci_matrix_uhf(rtdb,
     &      dim_fci,alpha_str,beta_str,dbl_mb(k_matrix),h,v,
     &      nos,nvs,noas,nobs,nvas,nvbs,nstot,repulsion,irrep_g)
c
c
c
c
      if(.not.MA_POP_STACK(l_matrix))
     &      call errquit('rrst3',4,MA_ERR)
c
c
c
c
c
      return 
      end
c
c
c
c
c ------------- linearized approach RRST(1)
c
c
c
c
c
       subroutine fci_matrix_uhf(rtdb,
     &      dim_fci,alpha_str,beta_str,matrix,h,v,
     &      nos,nvs,noas,nobs,nvas,nvbs,nstot,repulsion,irrep_g)
c
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c auxiliary integers
        integer i,j,k,l,m,n,iaux
        integer p,q,r,s,t,u,w,x,y,z
        integer noccp,nvirt,nh1,nv2
        integer pa,pb,qa,qb
        integer ps,qs,rs,ss
        integer ipl1,ipl2
c
        integer rtdb
        integer nstot
        integer dim_fci
        integer alpha_str
        integer beta_str
        integer nfill(nstot/2)
        integer qq(noas+nobs),so(noas+nobs)
        integer chi1,chi2
        double precision h(nstot,nstot)
        double precision ho(nstot/2,nstot/2)
        double precision hs(nstot,nstot)
        double precision hos(nstot/2,nstot/2)
        double precision v(nstot,nstot,nstot,nstot)
        double precision vc(nstot/2,nstot/2,nstot/2,nstot/2)
        double precision vx(nstot/2,nstot/2,nstot/2,nstot/2)
        double precision matrix(dim_fci,dim_fci)
        double precision trial(dim_fci)
        double precision xxx,yyy,zzz
        double precision repulsion
        double precision my_rhf
c
        integer nos,nvs,noas,nobs,nvas,nvbs
        integer irrep_g
        integer maxm
c
        integer qnumber
        integer k_str_a,l_str_a
        integer k_str_b,l_str_b
        integer half_str_a(alpha_str,noas+nvas)
        integer half_str_b(beta_str,nobs+nvbs)
        integer strings(dim_fci,nstot)
        integer string1(nstot)
        integer string2(nstot)
        integer string(nstot)
        integer string_diff(nstot)
        integer ind_set(noas+nobs)
        integer dim_fci_test
        integer iphase
        integer isum
c
        logical nodezero
        logical sym_test
        logical periodic

c  eigensolver
        double precision wr(dim_fci),wi(dim_fci)
        double precision vl,vr,work(3*dim_fci)
        integer info
c
c
        sym_test=.true.
        nodezero = (ga_nodeid().eq.0)
        periodic=.true.
c       periodic=.false.
c
c
c
*       ***********************
*       one-electron integral
*       ***********************
        open(101,file="one_electron_integrals.dat",status="old")
        read(101,*) noccp,nvirt,nh1

        do i=1,nh1
           read(101,*) p,q,ho(p,q)
           ho(q,p) = ho(p,q)
        end do

        close(101)

        do p=1,noccp+nvirt
        do q=1,noccp+nvirt
           if(p.le.noccp) then
              pa = p
              pb = p + noccp
           end if
           if(p.gt.noccp) then
              pa = p + noccp
              pb = p + noccp + nvirt
           end if
           if(q.le.noccp) then
              qa = q
              qb = q + noccp
           end if
           if(q.gt.noccp) then
              qa = q + noccp
              qb = q + noccp + nvirt
           end if

           h(pa,qa) = ho(p,q)
           h(pa,qb) = 0.0d0
           h(pb,qa) = 0.0d0
           h(pb,qb) = ho(p,q)
           h(qa,pa) = ho(q,p)
           h(qa,pb) = 0.0d0
           h(qb,pa) = 0.0d0
           h(qb,pb) = ho(q,p)
        enddo
        enddo

*       *******************************
*       one-electron screened integral
*       *******************************
        if(periodic) then
           open(105,file="one_electron_screened_integrals.dat",
     >              status="old")
           read(105,*) noccp,nvirt,nh1

           do i=1,nh1
              read(105,*) p,q,hos(p,q)
              hos(q,p) = hos(p,q)
           end do

           close(105)

           do p=1,noccp+nvirt
           do q=1,noccp+nvirt
              if(p.le.noccp) then
                 pa = p
                 pb = p + noccp
              end if
              if(p.gt.noccp) then
                 pa = p + noccp
                 pb = p + noccp + nvirt
              end if
              if(q.le.noccp) then
                 qa = q
                 qb = q + noccp
              end if
              if(q.gt.noccp) then
                 qa = q + noccp
                 qb = q + noccp + nvirt
              end if

              hs(pa,qa) = hos(p,q)
              hs(pa,qb) = 0.0d0
              hs(pb,qa) = 0.0d0
              hs(pb,qb) = hos(p,q)
              hs(qa,pa) = hos(q,p)
              hs(qa,pb) = 0.0d0
              hs(qb,pa) = 0.0d0
              hs(qb,pb) = hos(q,p)
          enddo
          enddo
        end if

*       ***********************
*       two-electron integral
*       ***********************
       open(102,file="two_electron_integrals_coulomb.dat",status="old")
       open(103,file="two_electron_integrals_exchange.dat",status="old")

        read(102,*) noccp,nvirt,nv2
        do i=1,nv2
           read(102,*) p,q,r,s,vc(p,q,r,s)
           vc(p,q,s,r) = vc(p,q,r,s)
           vc(q,p,r,s) = vc(p,q,r,s)
           vc(q,p,s,r) = vc(p,q,r,s)
           vc(r,s,p,q) = vc(p,q,r,s)
           vc(r,s,q,p) = vc(p,q,r,s)
           vc(s,r,p,q) = vc(p,q,r,s)
           vc(s,r,q,p) = vc(p,q,r,s)
        end do
        close(102)

        read(103,*) noccp,nvirt,nv2
        do i=1,nv2
           read(103,*) p,q,r,s,vx(p,q,r,s)
           vx(p,q,s,r) = vx(p,q,r,s)
           vx(q,p,r,s) = vx(p,q,r,s)
           vx(q,p,s,r) = vx(p,q,r,s)
           vx(r,s,p,q) = vx(p,q,r,s)
           vx(r,s,q,p) = vx(p,q,r,s)
           vx(s,r,p,q) = vx(p,q,r,s)
           vx(s,r,q,p) = vx(p,q,r,s)
        end do
        close(103)

*       ***********************
*       ion-ion energy 
*       ***********************
        open(104,file="ion_ion.dat")
        read(104,*) repulsion
        close(104)

        thres1 = noccp
        thres2 = noccp*2
        thres3 = noccp*2+nvirt
        
        my_rhf = 0.0d0
        do i=1,noas+nobs
         my_rhf = my_rhf + h(i,i)
        enddo
        do i=1,noas+nobs
        do j=1,noas+nobs
         my_rhf = my_rhf + (0.50d0)*vx(i,j,i,j)
        enddo
        enddo
c
        my_rhf = my_rhf + repulsion
c
c
        if(nodezero) then 
         write(6,*)'nos: ',nos
         write(6,*)'nvs: ',nvs
         write(6,*)'noas: ',noas
         write(6,*)'nobs: ',nobs
         write(6,*)'nvas: ',nvas
         write(6,*)'nvbs" ',nvbs
         write(6,*)'nstot = ',nstot
         write(6,*)'alpha_str = ',alpha_str
         write(6,*)'beta_str =  ',beta_str
         write(6,*)'dim_fci = ',dim_fci
         write(6,*)'repulsion = ',repulsion
         write(6,*)'my_rhf = ',my_rhf
         call util_flush(6)
        endif
c
        do i=1,noas+nvas
        do j=1,alpha_str
         half_str_a(j,i) = 0
        enddo
        enddo
c
        do i=1,nobs+nvbs
        do j=1,beta_str
         half_str_b(j,i) = 0
        enddo
        enddo
c
        do i=1,nstot
        do j=1,dim_fci
         strings(j,i) = 0
        enddo
        enddo
c
        do i=1,dim_fci
        do j=1,dim_fci
         matrix(j,i)=0.0d0
        enddo
        enddo
c
c
c
        sym_test = .true.
        dim_fci_test = 0
        alpha_str_test = 0
        beta_str_test = 0
c
      if (.not.MA_PUSH_GET(mt_int,(noas+nvas),'a_string',l_str_a,
     1 k_str_a)) CALL ERRQUIT('a_string',0,MA_ERR)
c ----- half_str_a starts here --------------------
        qnumber = 0
        do m = 0, noas+nvas-1
         qnumber = qnumber + 2**m
        enddo
c
        do i=1,qnumber   !main loop 
c
         do m=1,noas+nvas
          int_mb(k_str_a+m-1) = 0
         enddo
c         
         k=i
c
         do j = noas+nvas, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str_a+j-1) = 0
            else
             int_mb(k_str_a+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of alpha particle check
         m=0
         do j=1,noas+nvas
c additional symmetry checks possible here
          m = m+int_mb(k_str_a+j-1)
         enddo
         if(m.eq.noas) then
c additional symmetry checks possible here
           alpha_str_test = alpha_str_test + 1
           do n=1,noas+nvas
            half_str_a(alpha_str_test,n) = 
     1        int_mb(k_str_a+n-1)    
           enddo
         endif
c
       enddo !main loop
c
      if (.not.MA_POP_STACK(l_str_a))
     1   call errquit("l_str_a",0,MA_ERR)
c
      if(alpha_str_test.ne.alpha_str) 
     1   call errquit("dim_fci_test_a",0,MA_ERR)
c
c  ------- done with half_str_a -----------------------
c
      if (.not.MA_PUSH_GET(mt_int,(nobs+nvbs),'a_string',l_str_b,
     1 k_str_b)) CALL ERRQUIT('b_string',0,MA_ERR)
c ----- half_str_b starts here --------------------
        qnumber = 0
        do m = 0, nobs+nvbs-1
         qnumber = qnumber + 2**m
        enddo
c
        do i=1,qnumber   !main loop 
c
         do m=1,nobs+nvbs
          int_mb(k_str_b+m-1) = 0
         enddo
c         
         k=i
c
         do j = nobs+nvbs, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str_b+j-1) = 0
            else
             int_mb(k_str_b+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of beta particle check
         m=0
         do j=1,nobs+nvbs
c additional symmetry checks possible here
          m = m+int_mb(k_str_b+j-1)
         enddo
         if(m.eq.nobs) then
c additional symmetry checks possible here
           beta_str_test = beta_str_test + 1
           do n=1,nobs+nvbs
            half_str_b(beta_str_test,n) =
     1        int_mb(k_str_b+n-1)
           enddo
         endif
c
       enddo !main loop
c
      if (.not.MA_POP_STACK(l_str_b))
     1   call errquit("l_str_b",0,MA_ERR)
c
      if(beta_str_test.ne.beta_str)
     1   call errquit("dim_fci_test_b",0,MA_ERR)
c
c  ------- done with half_str_b -----------------------
c
c
c
c  ------- total number of strings 
c  
      m = 0
      do i=1,alpha_str
      do j=1,beta_str
c 
       do l=1,nstot
        string(l)=0
       enddo
c
       do l=1,noas
        string(l)=half_str_a(i,l)
       enddo
       do l=1,nvas
       string(l+nos)=half_str_a(i,l+noas)
       enddo
       do l=1,nobs
        string(l+noas)=half_str_b(j,l)
       enddo
       do l=1,nvbs
        string(l+nos+nvas)=half_str_b(j,l+nobs)
       enddo
c
         do l=1,noas+nobs
          ind_set(l)=0
         enddo
         k=1
         do l=1,nstot
          if(string(l).eq.1) then
           ind_set(k)=l
           k=k+1
          endif
         enddo
c --- safety check ---
      if(k-1.ne.noas+nobs) then
         write(6,*)'k-1 ',k-1
         write(6,*)'noas+nobs ',noas+nobs
         call util_flush(6)
         call errquit("noas_nobs_test_0",0,MA_ERR)
      endif
c --------------------
c
c determination of irrep_g
         irrep_g = 0
         do k=1,noas+nobs
          irrep_g = ieor(irrep_g,int_mb(k_irs_sorted+k-1))
         enddo
c
         irrep1=0
         do k=1,noas+nobs
          iaux = ind_set(k)
          irrep1 = ieor(irrep1,int_mb(k_irs_sorted+iaux-1))
         enddo
         if(irrep1.eq.irrep_g) then
          m = m + 1
          do k=1,nstot
           strings(m,k)=string(k)
          enddo
         endif
c
      enddo !j
      enddo !i
c
c
c
c ------- done with strings list --------------------
c
c upper half only
c
      do i=1,dim_fci  !main loop1
      do j=i,dim_fci  !main loop2
c
       idiff=0
       do k=1,nstot
        string1(k)=strings(i,k)
        string2(k)=strings(j,k)
        string_diff(k)=string1(k)-string2(k)
        idiff=idiff+abs(string_diff(k))
       enddo
c
       idiff=idiff/2
c ---------------------
       if(idiff.gt.2) then
        matrix(i,j)=0.0d0
        matrix(j,i)=0.0d0
       endif
c ---------------------
       if(idiff.eq.2) then
        p=0
        q=0
        r=0
        s=0
c       v^{pq}_{rs} a_p^+ a_q^+ a_s a_r (p<q & r<s)
        do k=1,nstot !k ----
         if(string_diff(k).eq.1) then
          if(p.eq.0) then
           p=k
          else  
           q=k
          endif
         endif
         if(string_diff(k).eq.-1) then
          if(r.eq.0) then
           r=k
          else
           s=k
          endif
         endif
        enddo !k ---
c
        iphase = 1
c r part 
        isum=0
        do k=1,r-1
         isum=isum+string2(k) 
        enddo 
        string2(r)=0
c s part
        do k=1,s-1
         isum=isum+string2(k) 
        enddo
        string2(s)=0
c q part
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=1
c p part
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
        iphase=(-1)**isum
c      
        if(p.le.noccp) m = p
        if(p.gt.noccp .and. p.le.(noccp*2)) m = p - noccp
        if(p.gt.(noccp*2) .and. p.le.(noccp*2+nvirt)) m = p - noccp
        if(p.gt.(noccp*2+nvirt)) m = p - noccp - nvirt

        if(q.le.noccp) n = q
        if(q.gt.noccp .and. q.le.(noccp*2)) n = q - noccp
        if(q.gt.(noccp*2) .and. q.le.(noccp*2+nvirt)) n = q - noccp
        if(q.gt.(noccp*2+nvirt)) n = q - noccp - nvirt

        if(r.le.noccp) u = r
        if(r.gt.noccp .and. r.le.(noccp*2)) u = r - noccp
        if(r.gt.(noccp*2) .and. r.le.(noccp*2+nvirt)) u = r - noccp
        if(r.gt.(noccp*2+nvirt)) u = r - noccp - nvirt

        if(s.le.noccp) w = s
        if(s.gt.noccp .and. s.le.(noccp*2)) w = s - noccp
        if(s.gt.(noccp*2) .and. s.le.(noccp*2+nvirt)) w = s - noccp
        if(s.gt.(noccp*2+nvirt)) w = s - noccp - nvirt

        if(p.le.thres1 .or. ((p.gt.thres2) .and. (p.le.thres3)))   ps=0
        if(((p.gt.thres1) .and. (p.le.thres2)) .or. (p.gt.thres3)) ps=1
        if(q.le.thres1 .or. ((q.gt.thres2) .and. (q.le.thres3)))   qs=0
        if(((q.gt.thres1) .and. (q.le.thres2)) .or. (q.gt.thres3)) qs=1
        if(r.le.thres1 .or. ((r.gt.thres2) .and. (r.le.thres3)))   rs=0
        if(((r.gt.thres1) .and. (r.le.thres2)) .or. (r.gt.thres3)) rs=1
        if(s.le.thres1 .or. ((s.gt.thres2) .and. (s.le.thres3)))   ss=0
        if(((s.gt.thres1) .and. (s.le.thres2)) .or. (s.gt.thres3)) ss=1

        if(ps.eq.rs .and. qs.eq.ss .and. ps.ne.qs) 
     >     v(p,q,r,s) = vx(m,u,n,w)
        if(ps.eq.ss .and. qs.eq.rs .and. ps.ne.qs) 
     >     v(p,q,r,s) = -vx(m,w,n,u)
        if(ps.eq.rs .and. qs.eq.ss .and. ps.eq.qs)
     >     v(p,q,r,s) = vx(m,u,n,w)-vx(m,w,n,u)

        matrix(i,j)=matrix(i,j)+dble(iphase)*v(p,q,r,s) 
        matrix(j,i)=matrix(j,i)+dble(iphase)*v(p,q,r,s)
c check
        do k=1,nstot
         if((string1(k)-string2(k)).ne.0) 
     1   call errquit("string_1_2_test",0,MA_ERR)
        enddo  
c
       endif !idiff=2
c ------------------------------------------------
       if(idiff.eq.1) then
        p=0
        q=0
        do k=1,nstot !k ----
         if(string_diff(k).eq.1) p=k
         if(string_diff(k).eq.-1) q=k
        enddo
c
         k=1
         do l=1,nstot
          if(string2(l).eq.1) then
           ind_set(k)=l
           k=k+1
          endif
         enddo
c --- safety check ----
      if(k-1.ne.noas+nobs) then
         write(6,*)'k-1 ',k-1
         write(6,*)'noas+nobs ',noas+nobs
         call util_flush(6)
         call errquit("noas_nobs_test_1",0,MA_ERR)
      endif
c ---------------------
c
        iphase = 1
c q part 
        isum=0
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=0
c p part        
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
         iphase=(-1)**isum
         xxx = 0.0d0
c
         do m=1,noas+nobs
          t=ind_set(m)
          if(p.le.noccp) w = p
          if(p.gt.noccp .and. p.le.(noccp*2)) w = p - noccp
          if(p.gt.(noccp*2) .and. p.le.(noccp*2+nvirt)) w = p - noccp
          if(p.gt.(noccp*2+nvirt)) w = p - noccp - nvirt

          if(q.le.noccp) n = q
          if(q.gt.noccp .and. q.le.(noccp*2)) n = q - noccp
          if(q.gt.(noccp*2) .and. q.le.(noccp*2+nvirt)) n = q - noccp
          if(q.gt.(noccp*2+nvirt)) n = q - noccp - nvirt

          if(t.le.noccp) u = t
          if(t.gt.noccp .and. t.le.(noccp*2)) u = t - noccp
          if(t.gt.(noccp*2) .and. t.le.(noccp*2+nvirt)) u = t - noccp
          if(t.gt.(noccp*2+nvirt)) u = t - noccp - nvirt

          if(p.le.thres1 .or. ((p.gt.thres2) .and. (p.le.thres3)))   
     >       ps=0
          if(((p.gt.thres1) .and. (p.le.thres2)) .or. (p.gt.thres3)) 
     >       ps=1
          if(q.le.thres1 .or. ((q.gt.thres2) .and. (q.le.thres3)))   
     >       qs=0
          if(((q.gt.thres1) .and. (q.le.thres2)) .or. (q.gt.thres3)) 
     >       qs=1
          if(t.le.thres1 .or. ((t.gt.thres2) .and. (t.le.thres3)))   
     >       ts=0
          if(((t.gt.thres1) .and. (t.le.thres2)) .or. (t.gt.thres3)) 
     >       ts=1

          if(ps.eq.ts .and. p.ne.t .and. q.ne.t) then
              v(p,t,q,t) = 0.50d0*vc(u,u,n,w) - vx(n,u,w,u)
          else if(ps.ne.ts) then
              v(p,t,q,t) = 0.50d0*vc(u,u,n,w)
          else if(p.eq.t .or. q.eq.t) then
              v(p,t,q,t) = 0.25d0*vc(w,n,w,n) - 0.5*vx(w,n,w,n)
          end if

          xxx = xxx + v(p,t,q,t)*dble(iphase) 
         enddo
c
         xxx = xxx + h(p,q)*dble(iphase)
         matrix(i,j) = matrix(i,j) + xxx
         matrix(j,i) = matrix(j,i) + xxx
c 
         endif !idiff=1
c ------------------------------------------------
       if(idiff.eq.0) then
         k=1
         do l=1,nstot
          if(string1(l).eq.1) then
           ind_set(k)=l
           k=k+1
          endif
         enddo
c -- safety check ----
      if(k-1.ne.noas+nobs) then
         write(6,*)'k-1 ',k-1
         write(6,*)'noas+nobs ',noas+nobs
         call util_flush(6)
         call errquit("noas_nobs_test_2",0,MA_ERR)
      endif

c --------------------
         do k=1,noas+nobs
          p=ind_set(k)
          matrix(i,i) = matrix(i,i) + h(p,p)
c         if(i.eq.2) then
c           write(*,*) "matrix22,h1",p,h(p,p),i,matrix(i,i)
c         end if
          do l=1,noas+nobs
           q=ind_set(l) 

           if(p.le.noccp) m = p
           if(p.gt.noccp .and. p.le.(noccp*2)) m = p - noccp
           if(p.gt.(noccp*2) .and. p.le.(noccp*2+nvirt)) m = p - noccp
           if(p.gt.(noccp*2+nvirt)) m = p - noccp - nvirt

           if(q.le.noccp) n = q
           if(q.gt.noccp .and. q.le.(noccp*2)) n = q - noccp
           if(q.gt.(noccp*2) .and. q.le.(noccp*2+nvirt)) n = q - noccp
           if(q.gt.(noccp*2+nvirt)) n = q - noccp - nvirt
c ---------------------
          if(p.le.thres1 .or. ((p.gt.thres2) .and. (p.le.thres3)))   
     >       ps=0
          if(((p.gt.thres1) .and. (p.le.thres2)) .or. (p.gt.thres3)) 
     >       ps=1
          if(q.le.thres1 .or. ((q.gt.thres2) .and. (q.le.thres3)))   
     >       qs=0
          if(((q.gt.thres1) .and. (q.le.thres2)) .or. (q.gt.thres3)) 
     >       qs=1

           if(ps.eq.qs) then  
               v(p,q,p,q) = 0.50d0*vc(m,m,n,n) - vx(n,m,m,n) 
           else if(ps.ne.qs) then
               v(p,q,p,q) = 0.50d0*vc(m,m,n,n)  
           end if

           v(q,p,q,p) = v(p,q,p,q)
           v(p,q,q,p) = v(p,q,p,q)
           v(q,p,p,q) = v(p,q,p,q)

c          if(i.eq.6) then
c             write(*,*) "matrix22,v2",p,q,v(p,q,p,q),i,matrix(i,i)
c             write(*,*) "matrix22,v2",m,n,vc(m,m,n,n),chi1,vx(m,n,m,n),
c    >                   chi2,vx(m,m,n,n),v(p,q,p,q)
c          end if
           matrix(i,i) = matrix(i,i) + (0.50d0)*v(p,q,p,q)
c          if(i.eq.8) then
c            write(*,*) "v2",p,q,v(p,q,p,q),i,matrix(i,i)
c          end if
          enddo
         enddo
         matrix(i,i) = matrix(i,i) + repulsion
c         if(i.eq.8) then
c          write(*,*) "diagonal ion-ion",repulsion,i,matrix(i,i)
c         end if
       endif !idiff=0
c ------------------------------------------------
c

c     if(i.eq.1) then
c        write(*,*) i,j,idiff
c     end if
      enddo           !main loop2
      enddo           !main loop1
c
c  Moments evaluation
c
        maxm=20
        do m=1,dim_fci
         trial(m) = 0.0d0
        enddo
        do i=1,dim_fci
         n=0
         do j=1,noas+nobs
          n = n + strings(i,j)
         enddo
         if(n.eq.noas+nobs) then
          trial(i)=1.0d0
          go to 30
         endif
        enddo
c
 30     continue
c
cccx        call moments_uhf(matrix,trial,dim_fci,maxm)
c
c diagonalization
c
c     open(303,file="FCI_matrix.dat")
c     do i=1,dim_fci
c     do j=1,dim_fci
c        write(303,*) i,j,matrix(i,j)
c     end do
c     end do
c     close(303)

      call util_dgeev('n','n',dim_fci,matrix,dim_fci,
     1  wr,wi,vl,1,vr,1,
     2  work,3*dim_fci,info)
c
c  *** debug ****
c      if(nodezero) then
c       write(6,*)' wr       wi '
c       do i=1,dim_fci
c        write(6,222) wr(i),wi(i)
c       enddo
c       call util_flush(6)
c      endif 
c 222  format(2f16.6)
c ***************
c
      xxx=wr(1)
c     write(*,*) "Eigenvalues="
      do i=1,dim_fci
c      write(*,*) wr(i)
       if(wr(i).lt.xxx) xxx=wr(i)
      enddo
c
      if(nodezero) then
      write(6,*)'----------'
      write(6,*)'lowest FCI energy: ',xxx
      write(6,*)'----------'
      call util_flush(6)
      endif
c
c
 303  format(f15.6,4x,f15.6)
c
        return 
        end
c
c
c
c
c
c    
c
      subroutine dim_det_uhf(noas,nobs,nvas,nvbs,nos,nvs,nstot,
     1           irrep_g,alpha_str,beta_str,dim_fci)
c
c
c nstot - total number of spinorbitals nstot=nos+nvs
c nos   - number of occupied spinorbitals
c nvs   - number of virtual spinorbitals
c noas  - number of occupied alpha spinorbitals
c nobs  - number of occupied beta spinorbitals
c nvas  - number of virtual  alhoa spinorbitals
c nvbs  - number of virtual beta spinorbitals
c
c spinorbital convention
c
c | noas | nobs | nvas | nvbs |
c
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
        integer nos,nvs,noas,nobs,nvas,nvbs,nstot
        integer i,j,k,l,m,n,iaux
        integer qnumber
        integer alpha_str,alpha_str_test
        integer beta_str,beta_str_test
        integer dim_fci
        integer half_str_a(alpha_str,noas+nvas)
        integer half_str_b(beta_str,nobs+nvbs)
        integer k_str_a,l_str_a
        integer k_str_b,l_str_b
        integer ind_set(noas+nobs)
        integer string(nstot)
        integer irrep1,irrep_g
c
        dim_fci = 0
        alpha_str_test = 0
        beta_str_test = 0
c
      if (.not.MA_PUSH_GET(mt_int,(noas+nvas),'a_string',l_str_a,
     1 k_str_a)) CALL ERRQUIT('a_string',0,MA_ERR)
c ----- half_str_a starts here --------------------
        qnumber = 0
        do m = 0, noas+nvas-1
         qnumber = qnumber + 2**m
        enddo
c
        do i=1,qnumber   !main loop 
c
         do m=1,noas+nvas
          int_mb(k_str_a+m-1) = 0
         enddo
c         
         k=i
c
         do j = noas+nvas, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str_a+j-1) = 0
            else
             int_mb(k_str_a+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of alpha particle check
         m=0
         do j=1,noas+nvas
c additional symmetry checks possible here
          m = m+int_mb(k_str_a+j-1)
         enddo
         if(m.eq.noas) then
c additional symmetry checks possible here
           alpha_str_test = alpha_str_test + 1
           do n=1,noas+nvas
            half_str_a(alpha_str_test,n) =
     1        int_mb(k_str_a+n-1)
           enddo
         endif
c
       enddo !main loop
c
      if (.not.MA_POP_STACK(l_str_a))
     1   call errquit("l_str_a",0,MA_ERR)
c
      if(alpha_str_test.ne.alpha_str)
     1   call errquit("dim_fci_test_adim",0,MA_ERR)
c
c ---------------------
c
c
      if (.not.MA_PUSH_GET(mt_int,(nobs+nvbs),'b_string',l_str_b,
     1 k_str_b)) CALL ERRQUIT('b_string',0,MA_ERR)
c ----- half_str_b starts here --------------------
        qnumber = 0
        do m = 0, nobs+nvbs-1
         qnumber = qnumber + 2**m
        enddo
c
        do i=1,qnumber   !main loop 
c
         do m=1,nobs+nvbs
          int_mb(k_str_b+m-1) = 0
         enddo
c         
         k=i
c
         do j = nobs+nvbs, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str_b+j-1) = 0
            else
             int_mb(k_str_b+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of beta particle check
         m=0
         do j=1,nobs+nvbs
c additional symmetry checks possible here
          m = m+int_mb(k_str_b+j-1)
         enddo
         if(m.eq.nobs) then
c additional symmetry checks possible here
           beta_str_test = beta_str_test + 1
           do n=1,nobs+nvbs
            half_str_b(beta_str_test,n) =
     1        int_mb(k_str_b+n-1)
           enddo
         endif
c
       enddo !main loop
c
      if (.not.MA_POP_STACK(l_str_b))
     1   call errquit("l_str_b",0,MA_ERR)
c
      if(beta_str_test.ne.beta_str)
     1   call errquit("dim_fci_test_bdim",0,MA_ERR)
c
c ------------------- done_kk
c
      dim_fci = 0
      do i=1,alpha_str
      do j=1,beta_str
c 
       do l=1,nstot
        string(l)=0
       enddo
c
       do l=1,noas
        string(l)=half_str_a(i,l)
       enddo
       do l=1,nvas
       string(l+nos)=half_str_a(i,l+noas)
       enddo
       do l=1,nobs
        string(l+noas)=half_str_b(j,l)
       enddo
       do l=1,nvbs
        string(l+nos+nvas)=half_str_b(j,l+nobs)
       enddo
c
         do l=1,noas+nobs
          ind_set(l)=0
         enddo
         k=1
         do l=1,nstot
          if(string(l).eq.1) then
           ind_set(k)=l
           k=k+1
          endif
         enddo
c determination of irrep_g
         irrep_g = 0 
         do k=1,noas+nobs
          irrep_g = ieor(irrep_g,int_mb(k_irs_sorted+k-1))
         enddo
c
         irrep1=0
         do k=1,noas+nobs
          iaux = ind_set(k)  
          irrep1 = ieor(irrep1,int_mb(k_irs_sorted+iaux-1))
         enddo
         if(irrep1.eq.irrep_g) dim_fci = dim_fci+1
      enddo !j
      enddo !i
c
        return
        end
