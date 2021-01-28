module sapt_exch
use types
use tran
use exmisc
use exi
use timing
use sapt_utils

implicit none

contains

subroutine e1exchs2(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: i, j, k, l, ia, jb
integer :: ij,ipr
integer :: ip, iq, ir, is
integer :: iv, iz, iu, it
integer :: iunit
integer :: dimOA,dimOB,NBas
integer :: GdimOA,GdimOB
double precision,allocatable :: S(:,:),Sab(:,:)
double precision,allocatable :: USa(:,:),USb(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:)
double precision,allocatable :: Kb(:,:)
double precision,allocatable :: JJb(:,:)
double precision,allocatable :: Qab(:,:),Qba(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:)
double precision,allocatable :: tmpA(:,:,:,:),tmpB(:,:,:,:), &
                                tmpAB(:,:,:,:)
double precision,allocatable :: work(:,:),RDM2Aval(:,:,:,:), &
                                RDM2Bval(:,:,:,:)
double precision :: tmp,ea,eb,exchs2
double precision :: t1(2),t2a(4),t2b(2),t2c(2),t2d
double precision :: t1f,t2f
double precision :: Tcpu,Twall
double precision,parameter :: Half=0.5d0
double precision,external  :: trace,FRDM2,FRDM2GVB
double precision,allocatable :: work1(:)

! set dimensions
 NBas = A%NBasis 
 dimOA = A%num0+A%num1
 dimOB = B%num0+B%num1
 !  
 !dimOA = A%INAct+A%NAct 
 !dimOB = B%INAct+B%NAct 
 GdimOA = A%num0+A%num1
 GdimOB = B%num0+B%num1

 !print*, 'A-MON'
 !print*, A%num0,A%num1
 !print*, A%INAct,A%NAct

 call clock('START',Tcpu,Twall)

 allocate(S(NBas,NBas),Sab(NBas,NBas),&
          PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          PAbb(NBas,NBas),PBaa(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas))
 allocate(USa(NBas,NBas),USb(NBas,NBas),&
          Qab(NBas,NBas),Qba(NBas,NBas),&
          Kb(NBas,NBas))
 allocate(work(NBas,NBas),tmp1(NBas,NBas),tmp2(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NBas,B%CMO,B%Occ,1d0,PB)

 call get_one_mat('S',S,A%Monomer,NBas)
 call tran2MO(S,A%CMO,B%CMO,Sab,NBas) 

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas) 
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas) 

 allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
          RDM2Bval(dimOB,dimOB,dimOB,dimOB))

 if(Flags%ICASSCF==1) then
    ! CAS
    RDM2Aval = A%RDM2val
    RDM2Bval = B%RDM2val
    ! old:
    !do l=1,dimOA
    !   do k=1,dimOA 
    !      do j=1,dimOA
    !         do i=1,dimOA
    !            RDM2Aval(i,j,k,l) = FRDM2(i,k,j,l,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
    !         enddo
    !      enddo
    !   enddo
    !enddo
    !do l=1,dimOB
    !   do k=1,dimOB 
    !      do j=1,dimOB
    !         do i=1,dimOB
    !            RDM2Bval(i,j,k,l) = FRDM2(i,k,j,l,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)
    !         enddo
    !      enddo
    !   enddo
    !enddo

 elseif(Flags%ICASSCF==0) then

    ! GVB
    RDM2Aval = A%RDM2val
    RDM2Bval = B%RDM2val

 endif

! USa=0; USb=0
! USa,USb in MOAO
! call dgemm('T','N',NBas,NBas,NBas,1d0,A%CMO,NBas,S,NBas,0d0,USa,NBas)
! call dgemm('T','N',NBas,NBas,NBas,1d0,B%CMO,NBas,S,NBas,0d0,USb,NBas)
! USa,USb in AOMO
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,A%CMO,NBas,0d0,USa,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,B%CMO,NBas,0d0,USb,NBas)

! PA(B), PB(A)
 call tran2MO(PA,USb,USb,PAbb,NBas) 
 call tran2MO(PB,USa,USa,PBaa,NBas) 

! Qab=0; Qba=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,USb,NBas,0d0,Qab,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,USa,NBas,0d0,Qba,NBas)

! old (too large)
 print*, 'A: num0, num1',A%num0,A%num1
 print*, 'B: num0, num1',B%num0,B%num1
 print*, 'dimOA,dimOB',dimOA,dimOB
! call tran3MO_Q(NBas,dimOA,A%CMO,Qba,'TWOA3B')
! call tran3MO_Q(NBas,dimOB,B%CMO,Qab,'TWOB3A')

 call tran4_gen(NBas,&
          dimOA,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          dimOA,  Qba(1:NBas,1:(A%num0+A%num1)),&
          dimOA,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          dimOA,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'TWOA3B','AOTWOSORT')
 call tran4_gen(NBas,&
          dimOB,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          dimOB,  Qab(1:NBas,1:(B%num0+B%num1)),&
          dimOB,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          dimOB,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'TWOB3A','AOTWOSORT')

 call make_K(NBas,PB,Kb)

! block
! integer :: ip,iq,ir,is
! integer :: NInte1,NInte2,NOcc
! double precision,allocatable :: TwoMO(:)
! double precision :: ETot
! integer,external :: NAddr3
! double precision,external :: FRDM2
!
! !NOcc=A%NAct+A%INAct
! NOcc=A%NAct+A%INAct
! NInte1 = NBas*(NBas+1)/2
! NInte2 = NInte1*(NInte1+1)/2
!
! allocate(TwoMO(NInte2))
!
! call LoadSaptTwoEl(A%Monomer,TwoMO,NBas,NInte2)
! ETot=0
! do ip=1,NOcc
!    do iq=1,NOcc
!      do ir=1,NOcc
!         do is=1,NOcc
!            ETot=ETot+FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)&
!            *TwoMO(NAddr3(ip,ir,iq,is))
!         enddo
!      enddo
!    enddo
! enddo
! print*, 'Check 2-el part of the energy: ',ETot
!
! deallocate(TwoMO)
!
! end block

! T1a
 t1 = 0
 t1(1) = SAPT%elst
 ! PA(ab).S(bc)
 ! S(ad).PB(dc)
 tmp1=0
 tmp2=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PB,NBas,0d0,tmp2,NBas)
 do j=1,NBas
    do i=1,NBas
       t1(2) = t1(2) + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 t1f = 2d0*t1(1)*t1(2) 
 !write(LOUT,*) 'T1 ',t1f
 if(SAPT%IPrint>=10) write(LOUT,'(/,1x,a,f16.8)') 'ExchS2(T1   ) = ', t1f*1000d0

! T2d
 t2d = -2d0*SAPT%Vnn*t1(2)
 !write(LOUT,*) 'T2d',t2d
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2d  ) = ', t2d*1000d0

! T2c
 t2c=0
 tmp1=0
 tmp2=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,Va,NBas,PB,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp2,NBas)
  do j=1,NBas
    do i=1,NBas
       t2c(1) = t2c(1) + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo
 t2c(1) = -2d0*t2c(1)
 !write(LOUT,*) 'T2c(1)',t2c(1) 
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2c-1) = ', t2c(1)*1000d0

! new
 t2c(2)=0
 do is=1,dimOB
    do iq=1,dimOB
       !t2c(2) = t2c(2) + sum(RDM2Bval(:,:,iq,is)*Vabb(:,:)*PAbb(is,iq))
       t2c(2) = t2c(2) + sum(RDM2Bval(1:dimOB,1:dimOB,iq,is)*Vabb(1:dimOB,1:dimoB)*PAbb(is,iq))
    enddo
 enddo    
 t2c(2) = -2d0*t2c(2)
 !write(LOUT,*) 'T2c(2) ',t2c(2)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2c-2) = ', t2c(2)*1000d0
! 
! T2b
 t2b=0
 tmp1=0
 tmp2=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,Vb,NBas,PA,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,S,NBas,0d0,tmp2,NBas)
  do j=1,NBas
    do i=1,NBas
       t2b(1) = t2b(1) + tmp1(i,j)*tmp2(i,j)
    enddo
  enddo
 t2b(1) = -2d0*t2b(1)
 !print*, 'T2b(1)',t2b(1)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2b-1) = ', t2b(1)*1000d0

!new
 t2b(2)=0
 do is=1,dimOA
    do iq=1,dimOA
       !t2b(2) = t2b(2) + sum(RDM2Aval(:,:,iq,is)*Vbaa(:,:)*PBaa(is,iq))
       t2b(2) = t2b(2) + sum(RDM2Aval(1:dimOA,1:dimOA,iq,is)*Vbaa(1:dimOA,1:dimOA)*PBaa(is,iq))
    enddo
 enddo    
 t2b(2) = -2d0*t2b(2)
 !write(LOUT,*) 'T2b(2) ',t2b(2)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2b-2) = ', t2b(2)*1000d0

! T2a 
 t2a=0
 do jb=1,NBas
    do ia=1,NBas
       t2a(1) = t2a(1) + PA(ia,jb)*Kb(jb,ia)
    enddo
 enddo
 t2a(1) = -2.0d0*t2a(1)
 !write(LOUT,*) 'T2a(1)',t2a(1)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2a-1) = ', t2a(1)*1000d0

 open(newunit=iunit,file='TWOA3B',status='OLD',&
     !access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)
     access='DIRECT',form='UNFORMATTED',recl=8*dimOA**2)

! new
! Qba
 do ir=1,dimOA
    do ip=1,dimOA
       read(iunit,rec=ip+(ir-1)*dimOA) work(1:dimOA,1:dimOA)
           do is=1,dimOA
              do iq=1,dimOA
                     t2a(2) = t2a(2) + work(iq,is)* &
                              RDM2Aval(ip,ir,iq,is)
              enddo
           enddo
    enddo
 enddo

!! old
!! Qba 
! do ir=1,NBas
!    do ip=1,ir
!      read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOA,1:dimOA)
!
!       if(ip==ir) then
!
!         if(ip<=dimOA) then 
!            do is=1,dimOA
!               do iq=1,dimOA
!                    t2a(2) = t2a(2) + work(iq,is)* &
!                             RDM2Aval(ip,ir,iq,is)
!                enddo
!             enddo
!         else
!             do iq=1,dimOA
!                  t2a(2) = t2a(2) + work(iq,iq)* &
!                           2d0*A%Occ(ip)*A%Occ(iq)
!             enddo
!         endif 
!      
!       else
!
!          if(ir<=dimOA) then
!             do is=1,dimOA 
!                do iq=1,dimOA 
!                     t2a(2) = t2a(2) + work(iq,is)* &
!                            (RDM2Aval(ip,ir,iq,is)+RDM2Aval(ir,ip,iq,is))
!                enddo
!            enddo
!          endif
!
!       endif
!
!    enddo
! enddo

 t2a(2) = -2*t2a(2)
 ! write(LOUT,*) 'T2a(2) ',t2a(2)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2a-2) = ', t2a(2)*1000d0
 close(iunit,status='DELETE')

 open(newunit=iunit,file='TWOB3A',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)
     !access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)
! new
! Qab
 do ir=1,dimOB
    do ip=1,dimOB
       read(iunit,rec=ip+(ir-1)*dimOB) work(1:dimOB,1:dimOB)
           do is=1,dimOB
              do iq=1,dimOB
                     t2a(3) = t2a(3) + work(iq,is)* &
                              RDM2Bval(ip,ir,iq,is)
              enddo
           enddo
    enddo
 enddo

!! old
!! Qab
! do ir=1,NBas
!    do ip=1,ir
!      read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOB,1:dimOB)
!
!       if(ip==ir) then
!
!         if(ip<=dimOB) then 
!            do is=1,dimOB
!               do iq=1,dimOB
!                    t2a(3) = t2a(3) + work(iq,is)* &
!                             RDM2Bval(ip,ir,iq,is)
!                enddo
!             enddo
!         else
!             do iq=1,dimOB
!                  t2a(3) = t2a(3) + work(iq,iq)* &
!                           2d0*B%Occ(ip)*B%Occ(iq)
!             enddo
!         endif 
!      
!       else
!
!          if(ir<=dimOB) then
!             do is=1,dimOB 
!                do iq=1,dimOB 
!                     t2a(3) = t2a(3) + work(iq,is)* &
!                              (RDM2Bval(ip,ir,iq,is)+RDM2Bval(ir,ip,iq,is))
!                enddo
!             enddo
!          endif
!
!       endif
!
!    enddo
! enddo

 t2a(3) = -2*t2a(3)
 !write(LOUT,*) 'T2a(3) ',t2a(3)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2a-3) = ', t2a(3)*1000d0
 close(iunit,status='DELETE')

! T2a(4)
 allocate(tmpA(dimOA,dimOA,dimOA,dimOB),tmpB(dimOB,dimOB,dimOA,dimOB),&
          tmpAB(dimOA,dimOA,dimOB,dimOB))
!
! Full NBas check:
!! N^5 
! tmpA = 0
! do iz=1,NBas
!    do ir=1,NBas
!       do iq=1,NBas
!          do ip=1,NBas
!             do is=1,NBas
!                 tmpA(ip,iq,ir,iz) = tmpA(ip,iq,ir,iz) + &
!                                  Sab(is,iz)* &
!                                  FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!! N^5
!tmpB = 0
! do iq=1,NBas
!    do iu=1,NBas
!       do iz=1,NBas
!          do iv=1,NBas
!             do it=1,NBas
!                 tmpB(iv,iz,iu,iq) = tmpB(iv,iz,iu,iq) + &
!                                  Sab(iq,it)* &
!                                  FRDM2(iv,iz,iu,it,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!
!! N^6
! tmpAB=0
!
! do iu=1,NBas
!    do iv=1,NBas
!       do ir=1,NBas
!          do ip=1,NBas
!             do iz=1,NBas
!                do iq=1,NBas
!                   tmpAB(ip,ir,iv,iu) = tmpAB(ip,ir,iv,iu) + &
!                                        tmpA(ip,iq,ir,iz)*tmpB(iv,iz,iu,iq)
!                enddo
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!
! work=0
! open(newunit=iunit,file='TMPMOAB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)
!
! do ir=1,NBas
!    do ip=1,NBas
!       read(iunit,rec=ip+(ir-1)*NBas) work
!      
!       do iv=1,NBas
!          do iu=1,NBas
!             t2a(4)=t2a(4)+work(iv,iu)* &
!                    tmpAB(ip,ir,iv,iu)
!          enddo
!       enddo
!
!     enddo
! enddo
! t2a(4) = -2*t2a(4)
! print*, 't2a(4): ',t2a(4)

 ! dimOA, dimOB 
 ! old
 !tmpA = 0
 !do iz=1,dimOB
 !   do ir=1,dimOA
 !      do iq=1,dimOA
 !         do ip=1,dimOA
 !            do is=1,dimOA
 !                tmpA(ip,iq,ir,iz) = tmpA(ip,iq,ir,iz) + &
 !                                 Sab(is,iz)* &
 !                                 !FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
 !                                 RDM2Aval(ip,ir,iq,is)
 !            enddo
 !         enddo
 !      enddo
 !   enddo
 !enddo
 ! new
 call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,RDM2Aval,dimOA**3,Sab,NBas,0d0,tmpA,dimOA**3)

! old:
! N^5
! tmpB = 0
! do iq=1,dimOA
!    do iu=1,dimOB
!       do iz=1,dimOB
!          do iv=1,dimOB
!             do it=1,dimOB
!                 tmpB(iv,iz,iu,iq) = tmpB(iv,iz,iu,iq) + &
!                                  Sab(iq,it)* &
!                                  !FRDM2(iv,iz,iu,it,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)
!                                  RDM2Bval(iv,iu,iz,it)
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
! new
 do is=1,dimOB
    call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,RDM2Bval(:,:,:,is),dimOB**2,Sab,NBas,0d0,tmpB(:,:,:,is),dimOB**2)
 enddo

! old:
! N^6
! tmpAB=0
! do iu=1,dimOB
!    do iv=1,dimOB
!       do ir=1,dimOA
!          do ip=1,dimOA
!             do iz=1,dimOB
!                do iq=1,dimOA
!                   tmpAB(ip,ir,iv,iu) = tmpAB(ip,ir,iv,iu) + &
!                                        tmpA(ip,ir,iq,iz)*tmpB(iv,iu,iz,iq)
!                enddo
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
! new
 call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,tmpA,dimOA**2,tmpB,dimOB**2,0d0,tmpAB,dimOA**2)

 !do is=1,dimOB 
 !   do iq=1,dimOB 
 !      do ir=1,dimOA 
 !         do ip=1,dimOA 
 !         write(LOUT,'(1x,a,4i2,f12.6)') 'ip,ir,iq,is',ip,ir,iq,is,tmpAB(ip,ir,iq,is)
 !         enddo
 !      enddo
 !   enddo
 !enddo
 !print*, 'tmpAB-exch',norm2(tmpAB)

 work=0
 open(newunit=iunit,file='TMPOOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*GdimOB**2)

! do ir=1,A%INAct+A%NAct
!    do ip=1,A%INAct+A%NAct
 do ir=1,dimOA
    do ip=1,dimOA
!  print*, ip,ir,ip+(ir-1)*GdimOA 
      read(iunit,rec=ip+(ir-1)*GdimOA) work(1:GdimOB,1:GdimOB)

      t2a(4) = t2a(4) + sum(work(1:dimOB,1:dimOB)*tmpAB(ip,ir,1:dimOB,1:dimOB))

    enddo
 enddo
 close(iunit)
 t2a(4) = -2*t2a(4)
 ! write(LOUT,*) 't2a(4): ',t2a(4)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2a-4) = ', t2a(4)*1000d0

 deallocate(tmpAB,tmpB,tmpA)

 exchs2=t1f+sum(t2a)+sum(t2b)+sum(t2c)+t2d
 SAPT%exchs2 = exchs2
 write(LOUT,'(/1x,a,f16.8)') 'ExchS2        = ', exchs2*1000d0 

 deallocate(Vbaa,Vabb,PBaa,PAbb,Vb,Va,PB,PA,Sab,S)
 deallocate(Kb,Qba,Qab,USb,USa) 
 deallocate(tmp2,tmp1,work)
 !deallocate(JJb) 
 deallocate(RDM2Bval,RDM2Aval)

 call clock('E1exch(S2)',Tcpu,Twall)

end subroutine e1exchs2

subroutine hl_2el(Flags,A,B,SAPT)
! calculate Heitler-London energy
! for 2-electron CAS monomers
implicit none

type(FlagsData)   :: Flags
type(SaptData)    :: SAPT
type(SystemBlock) :: A, B

integer :: i,j,ij,k,l,kl,ia,jb
integer :: ip,iq,ir,is
integer :: iunit
integer :: dimOA,dimOB,NBas
double precision :: val,fac,tmp
double precision :: ccs2,nns2,ccaaaa,ccbbbb
double precision :: t12(3),t34(3),t13(6)
double precision :: Dfull,Ds2inv,e1,e2,ehl,ehlint
double precision :: e1s2,e2s2,ehls2,ehls2int
double precision :: deltaMs2
double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:) 
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Ha(:,:),Haa(:,:),&
                                HB(:,:),Hbb(:,:),Hab(:,:)
double precision,allocatable :: Kb(:,:),Jmat(:,:)
double precision,allocatable :: work(:,:),work1(:),ints(:,:)
double precision,allocatable :: RDM2Aval(:,:,:,:), &
                                RDM2Bval(:,:,:,:)
double precision,parameter   :: Half=0.5d0
double precision,external    :: trace
! test
double precision,allocatable :: tmpAB(:,:,:,:),Sab_save(:,:)
double precision,allocatable :: work2(:),work3(:)

! info
 write(LOUT,'(/,1x,a,/)') 'HEITLER-LONDON ENERGY FOR 2-el MONOMERS'

! set dimensions
 NBas = A%NBasis 
 dimOA = A%num0+A%num1
 dimOB = B%num0+B%num1

 allocate(S(NBas,NBas),Sab(NBas,NBas),&
          PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          PAbb(NBas,NBas),PBaa(NBas,NBas))
 allocate(work(NBas,NBas),work1(NBas*NBas),&
          Ha(NBas,NBas),Haa(NBas,NBas),&
          Hb(NBas,NBas),Hbb(NBas,NBas),&
          Hab(Nbas,NBas))
 allocate(Sba(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NBas,B%CMO,B%Occ,1d0,PB)

 call get_one_mat('S',S,A%Monomer,NBas)
 call tran2MO(S,A%CMO,B%CMO,Sab,NBas) 
 call tran2MO(S,B%CMO,A%CMO,Sba,NBas) 

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call get_one_mat('H',Ha,A%Monomer,NBas)
 call get_one_mat('H',Hb,B%Monomer,NBas)

 work = 0
 work = Ha + Vb
 !work = Vb
 call tran2MO(work,A%CMO,A%CMO,Haa,NBas)

 work = 0
 !work = Va
 work = Hb + Va
 call tran2MO(work,B%CMO,B%CMO,Hbb,NBas)

 call tran2MO(work,A%CMO,B%CMO,Hab,NBas)

 allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
          RDM2Bval(dimOB,dimOB,dimOB,dimOB))

 ! CAS
 RDM2Aval = A%RDM2val
 RDM2Bval = B%RDM2val

 ! denominator
 ccs2 = 0
 nns2 = 0
 do iq=1,NBas
    do ip=1,NBas
       nns2 = nns2 + A%Occ(ip)*B%Occ(iq)*Sab(ip,iq)*Sab(ip,iq)
       ccs2 = ccs2 + A%CICoef(ip)*B%CICoef(iq)*Sab(ip,iq)*Sab(ip,iq)
    enddo
 enddo
 Dfull  = 1d0 - 2d0 * nns2 + ccs2**2
 Ds2inv = 1d0 + 2d0 * nns2 

 print*, 'D      ', Dfull
 print*, '1/D(S2)', Ds2inv
 print*, ''

 ! <Psi0|HA|Psi0>
 ! 1-el part
 e1   = 0 
 e1s2 = 0
 val  = 0
 do ip=1,NBas
    val = val + A%Occ(ip)*Haa(ip,ip) + B%Occ(ip)*Hbb(ip,ip)
 enddo
 e1 = 2d0*val

! extra terms in 1-el S2
 e1s2 = 4d0*nns2*val
 !!! TEST - elst
 !e1s2 = e1

 val = 0
 do iq=1,NBas
    do ip=1,NBas
       val = val + A%Occ(ip)*B%Occ(iq)*Hab(ip,iq)*Sab(ip,iq)
    enddo
 enddo
 e1 = e1 - 4d0*val

 val = 0
 do ir=1,NBas
    do ip=1,NBas

       fac = 2d0*A%CICoef(ip)*A%CICoef(ir)*Haa(ip,ir)

       val = 0
       do iq=1,NBas
          val = val + B%Occ(iq)*Sab(ip,iq)*Sab(ir,iq)
       enddo
       e1 = e1 - fac*val

    enddo
 enddo
 !print*, 'E1-3',e1

 val = 0
 do ir=1,NBas
    do ip=1,NBas

       fac = 2d0*B%CICoef(ip)*B%CICoef(ir)*Hbb(ip,ir)

       val = 0
       do iq=1,NBas
          val = val + A%Occ(iq)*Sab(iq,ip)*Sab(iq,ir)
       enddo
       e1 = e1 - fac*val

    enddo
 enddo

 ! S2 approx
 e1s2 = e1s2 + e1

 val = 0
 do iq=1,NBas
    do ip=1,NBas
       val = val + A%CICoef(ip)*B%CICoef(iq)*Hab(ip,iq)*Sab(ip,iq)
    enddo
 enddo
 e1 = e1 + 4d0*ccs2*val 
 print*, '1-el part',e1

 ! two-electron part

 ! test Coulomb and exchange
 allocate(Kb(NBas,NBas),Jmat(NBas,NBas))
 call make_K(NBas,PB,Kb)
 call make_J1(NBas,PB,Jmat,'AOTWOSORT')

 tmp = 0
 do jb=1,NBas
    do ia=1,NBas
       tmp = tmp + PA(ia,jb)*Jmat(jb,ia)
    enddo
 enddo
 tmp = 4.0d0*tmp
 write(LOUT,*) 'Coulomb ',tmp

 tmp = 0
 do jb=1,NBas
    do ia=1,NBas
       tmp = tmp + PA(ia,jb)*Kb(jb,ia)
    enddo
 enddo
 tmp = -2.0d0*tmp
 write(LOUT,*) 'Exchange',tmp
 deallocate(Jmat,Kb)
 ! end test

! Coulomb and exchange
 call tran4_gen(NBas,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'OOOOAABB','AOTWOSORT')
 call tran4_gen(NBas,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'OOOOABAB','AOTWOSORT')

 !call tran4_gen(NBas,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         NBas,A%CMO,&
 !         NBas,B%CMO,&
 !         'FFOOABAB','AOTWOSORT')

! remaining ints
 call tran4_gen(NBas,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'OOOOABBB','AOTWOSORT')
 call tran4_gen(NBas,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'OOOOBAAA','AOTWOSORT')
 call tran4_gen(NBas,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'OOOOAAAA','AOTWOSORT')
 call tran4_gen(NBas,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'OOOOBBBB','AOTWOSORT')

 e2    = 0
 e2s2  = 0

 t12 = 0
 t34 = 0
 t13 = 0

 allocate(ints(NBas,NBas))

! ! test exchange
! !(FF|OO):(AB|AB)
!  open(newunit=iunit,file='FFOOABAB',status='OLD', &
!     access='DIRECT',recl=8*NBas*NBas)
!
! work1 = 0
! ints  = 0
! kl    = 0
! tmp =0 
! do l=1,dimOB
!    do k=1,dimOA
!       kl = kl + 1
!       read(iunit,rec=kl) work1(1:NBas*NBas)
! 
!          do j=1,NBas
!             do i=1,NBas
!                ints(i,j) = work1((j-1)*NBas+i)
!             enddo
!          enddo
!
!          ip = k
!          iq = l
!
!          tmp = tmp + A%Occ(ip)*B%Occ(iq)*ints(ip,iq)
!
!    enddo
! enddo
! print*, 'EXCH-TEST',-2d0*tmp
!
! close(iunit)
! !!!!! end test

 !(OO|OO):(AA|AA)
 open(newunit=iunit,file='OOOOAAAA',status='OLD', &
     access='DIRECT',recl=8*dimOA*dimOA)

 work1  = 0
 ints   = 0
 ccaaaa = 0
 kl     = 0
 do l=1,dimOA
    do k=1,dimOA
       kl = kl + 1
       read(iunit,rec=kl) work1(1:dimOA*dimOA)

       iq = l
       ip = k

       ij = 0
       do j=1,dimOA
          do i=1,dimOA
             ij = ij + 1
             !ints(i,j) = work1((j-1)*dimOA+i)
             ints(i,j) = work1(ij)
          enddo
       enddo

       ccaaaa = ccaaaa + A%CICoef(ip)*A%CICoef(iq)*ints(ip,iq)
       
    enddo
 enddo
 t12(1) = t12(1) + ccaaaa
 close(iunit)

!(OO|OO):(BB|BB)
 open(newunit=iunit,file='OOOOBBBB',status='OLD', &
     access='DIRECT',recl=8*dimOB*dimOB)

 work1  = 0
 ints   = 0
 ccbbbb = 0
 kl     = 0
 do l=1,dimOB
    do k=1,dimOB
       kl = kl + 1
       read(iunit,rec=kl) work1(1:dimOB*dimOB)

       iq = l
       ip = k

       ij = 0
       do j=1,dimOB
          do i=1,dimOB
             ij = ij + 1
             !ints(i,j) = work1((j-1)*dimOB+i)
             ints(i,j) = work1(ij)
          enddo
       enddo

       ccbbbb = ccbbbb + B%CICoef(ip)*B%CICoef(iq)*ints(ip,iq)

    enddo
 enddo
 t34(1) = t34(1) + ccbbbb
 close(iunit)

 !(OO|OO):(AB|BB)
 open(newunit=iunit,file='OOOOABBB',status='OLD', &
     access='DIRECT',recl=8*dimOA*dimOB)

 work1 = 0
 ints  = 0
 kl    = 0
 do l=1,dimOB
    do k=1,dimOB
       kl = kl + 1
       read(iunit,rec=kl) work1(1:dimOA*dimOB)

       ir = l
       iq = k

       fac = B%CICoef(iq)*B%CICoef(ir)

       ij = 0
       do j=1,dimOB
          do i=1,dimOA
             ij = ij + 1
             !ints(i,j) = work1((j-1)*dimOB+i)
             ints(i,j) = work1(ij)
          enddo
       enddo

       do ip=1,dimOA
          t13(3) = t13(3) + fac*A%Occ(ip)*Sab(ip,iq)*ints(ip,ir) 
       enddo

    enddo
 enddo
 close(iunit)
 t13(3) = -2d0*t13(3)
 t34(2) = t13(3)

 !(OO|OO):(BA|AA)
 open(newunit=iunit,file='OOOOBAAA',status='OLD', &
     access='DIRECT',recl=8*dimOA*dimOB)

 work1 = 0
 ints  = 0
 kl    = 0
 do l=1,dimOA
    do k=1,dimOA
       kl = kl + 1
       read(iunit,rec=kl) work1(1:dimOA*dimOB)

       ir = l
       ip = k

       fac = A%CICoef(ip)*A%CICoef(ir)

       ij = 0
       do j=1,dimOA
          do i=1,dimOB
             ij = ij + 1
             !ints(i,j) = work1((j-1)*dimOA+i)
             ints(i,j) = work1(ij)
          enddo
       enddo

       do iq=1,dimOB
          t13(4) = t13(4) + fac*B%Occ(iq)*Sab(ip,iq)*ints(iq,ir)
       enddo

    enddo
 enddo
 close(iunit)
 t13(4) = -2d0*t13(4)
 t12(2) = t13(4)

 ! test!
 open(newunit=iunit,file='TMPOOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)

 allocate(tmpAB(dimOA,dimOA,dimOB,dimOB),Sab_save(Nbas,Nbas))

 tmpAB=0
 do is=1,dimOB
    do iq=1,dimOB
       do ir=1,dimOA
          do ip=1,dimOA
             tmpAB(ip,ir,iq,is) = tmpAB(ip,ir,iq,is) + & 
                                  A%CICoef(ip)*A%CICoef(ir)*Sab(ip,iq) &
                                * B%CICoef(iq)*B%CICoef(is)*Sab(ir,is)
          enddo
       enddo
    enddo
 enddo

 do is=1,dimOB 
    do iq=1,dimOB 
       do ir=1,dimOA 
          do ip=1,dimOA 
          write(LOUT,'(1x,a,4i2,f12.6)') 'ip,ir,iq,is',ip,ir,iq,is,tmpAB(ip,ir,iq,is)
          enddo
       enddo
    enddo
 enddo

 print*, 'tmpAB-hl2e',norm2(tmpAB)

 ints = 0

 do ir=1,dimOA
    do ip=1,dimOA
      read(iunit,rec=ip+(ir-1)*dimOA) ints(1:dimOB,1:dimOB)

      t13(5) = t13(5) + sum(ints(1:dimOB,1:dimOB)*tmpAB(ip,ir,1:dimOB,1:dimOB))
      !t13(5) = t13(5) + sum(tmpAB(ip,ir,1:dimOB,1:dimOB))

    enddo
 enddo
 t13(5) = -2d0*t13(5)
 !print*, 'T2a-4',t13(5) 
 deallocate(tmpAB,Sab_save)

 close(iunit)

 !(OO|OO):(AA|BB)
 open(newunit=iunit,file='OOOOAABB',status='OLD', &
     access='DIRECT',recl=8*dimOA*dimOA)

 work1 = 0
 ints  = 0
 kl    = 0
 val   = 0
 tmp   = 0
 do l=1,dimOB
    do k=1,dimOB

!       kl = kl + 1
!       read(iunit,rec=kl) work1(1:dimOA*dimOA)
       read(iunit,rec=k+(l-1)*dimOB) ints(1:dimOA,1:dimOA)

       is = l
       iq = k

       fac = B%CICoef(iq)*B%CICoef(is)

      ! do ir=1,dimOA
      !    do ip=1,dimOA
      !    !t13(5) = t13(5) + fac*A%CICoef(ip)*A%CICoef(ir)*Sab(ip,iq)*Sab(ir,is)*ints(ip,ir)
      !    t13(5) = t13(5) + ints(ip,ir)
      !    enddo 
      ! enddo 

       if(k==l) then

          ip = k

          do iq=1,dimOA
             t13(1) = t13(1) + A%Occ(iq)*B%Occ(ip)*ints(iq,iq)
          enddo

       endif
    enddo
 enddo
 close(iunit)
 t13(1) = 4d0*t13(1)
 !t13(5) = -2d0*t13(5)

 !(OO|OO):(AB|AB)
 open(newunit=iunit,file='OOOOABAB',status='OLD', &
     access='DIRECT',recl=8*dimOA*dimOB)

 work1 = 0
 ints  = 0
 kl    = 0
 tmp   = 0
 do l=1,dimOB
    do k=1,dimOA

       kl = kl + 1
       read(iunit,rec=kl) work1(1:dimOA*dimOB)

       ip = l
       iq = k

       ij = 0
       do j=1,dimOB
          do i=1,dimOA
             ij = ij + 1
             !ints(i,j) = work1((j-1)*dimOB+i)
             ints(i,j) = work1(ij)
          enddo
       enddo

       ! t13b
       t13(2) = t13(2) + A%Occ(iq)*B%Occ(ip)*ints(iq,ip)
       ! t12c = t34c
       t12(3) = t12(3) + A%CICoef(iq)*B%CICoef(ip)*ints(iq,ip)

       is = l
       ir = k
       
       fac = A%CICoef(ir)*B%CICoef(is)*Sab(ir,is)

       do iq=1,dimOB
          do ip=1,dimOA
             t13(6) = t13(6) + fac*A%CICoef(ip)*B%CICoef(iq)*Sab(ip,iq)*ints(ip,iq)
          enddo
       enddo

    enddo
 enddo
 close(iunit)
 t12(3) = ccs2*t12(3)
 t34(3) = t12(3)

 t13(2) = -2d0*t13(2)
 t13(6) =  2d0*t13(6)

 print*, ''
 print*, 'T12(a)',t12(1)
 print*, 'T12(b)',t12(2)
 print*, 'T12(c)',t12(3)

 print*, 'T34(a)',t34(1)
 print*, 'T34(b)',t34(2)
 print*, 'T34(c)',t34(3)
 print*, 'T13(a)[Coul]',t13(1)
 print*, 'T13(b)[Exch]',t13(2)
 print*, 'T13(c)      ',t13(3)
 print*, 'T13(d)      ',t13(4)
 print*, 'T13(e)      ',t13(5)
 print*, 'T13(f)      ',t13(6)

 !print*, 'T12-sum',sum(t12)
 !print*, 'T34-sum',sum(t34)
 !print*, 'T13-sum',sum(t13)

 e2 = sum(t12) + sum(t34) + sum(t13)
 e2s2 = sum(t12) + sum(t34) - 2d0*t12(3) + sum(t13) - t13(6)
 e2s2 = e2s2 + 2d0*nns2*ccaaaa + 2d0*nns2*ccbbbb + 2d0*nns2*t13(1)

 print*, ''
 print*, 'E1    ',e1
 print*, 'E1(S2)',e1s2

 print*, 'E2    ',e2
 print*, 'E2(S2)',e2s2
 
 eHL = (e1+e2)/Dfull
 ehls2 = e1s2 + e2s2

 print*, 'E(A),E(B)', A%ECASSCF,B%ECASSCF
 print*, 'E(A)+E(B)', A%ECASSCF+B%ECASSCF

 print*, ''
 print*, 'eHL    ', ehl
 print*, 'eHL(S2)', ehls2

 ehlint = ehl - A%ECASSCF - B%ECASSCF + SAPT%Vnn + SAPT%monA%PotNuc + SAPT%monB%PotNuc
 !ehlint = ehl - A%ECASSCF - B%ECASSCF + SAPT%Vnn
 ehls2int = ehls2 - A%ECASSCF - B%ECASSCF + SAPT%Vnn + SAPT%monA%PotNuc + SAPT%monB%PotNuc
 print*, ''
 print*, 'eHLint       ', ehlint*1000d0 
 print*, 'eHL(S2)int   ', ehls2int*1000d0
 print*, 'ELST+EXCH(s2)',(SAPT%elst+SAPT%exchs2)*1000d0

 !print*, 'HL-ELST    ',ehlint-SAPT%elst
 !print*, 'HL(S2)-ELST',ehls2int-SAPT%elst
 !
 !tmp = ehls2int-SAPT%elst-SAPT%exchs2
 !print*, 'test     ',tmp
 tmp = ehlint-SAPT%elst-SAPT%exchs2
 print*, 'delta     ',tmp*1000d0

 print*, 'TEST-exch'
 t12(2) = SAPT%exch_part(2)
 t13(2) = SAPT%exch_part(1)
 t34(2) = SAPT%exch_part(3)
 t13(5) = SAPT%exch_part(4)
 t13(6) = SAPT%exch_part(5)
 e2 = sum(t12) + sum(t34) + sum(t13)
 eHL = (e1+e2)/Dfull
 print*, 'eHL-tst', ehl
 ehlint = ehl - A%ECASSCF - B%ECASSCF + SAPT%Vnn + SAPT%monA%PotNuc + SAPT%monB%PotNuc
 print*, 'eHLint       ', ehlint*1000d0 

 !! TESTY S2
 write(LOUT,'(/,1x,a)') 'TESTY-S2:'
 print*, 'Pb.Ja',t13(1)
 print*, 'nns2 ',nns2
 print*, 'T1   ',2d0*nns2*t13(1)
 print*, 'T2a-1',t13(2)
 print*, 'T2a-2',t12(2)
 print*, 'T2a-3',t34(2)
 print*, 'T2a-4',t13(5)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! test eq 58
 !! MONOMER A
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !write(LOUT,'(/,1x,a)') 'Test Eq. 58' 
 !Haa = 0; Hab = 0
 !call tran2MO(Ha,A%CMO,A%CMO,Haa,NBas)
 !call tran2MO(Ha,B%CMO,A%CMO,Hab,NBas)

 !! test 1-el
 !tmp=0
 !do i=1,dimOA
 !   tmp = tmp + 2d0*A%Occ(i)*Haa(i,i)
 !enddo
 !print*, 'OneEl(A)',tmp

 !work = 0
 !do ir=1,NBas
 !   do iq=1,dimOA
 !   !do iq=1,NBas 

 !     work(iq,ir) = work(iq,ir) +  A%CICoef(iq)*Hab(ir,iq) 

 !     val = 0
 !     do ip=1,dimOA
 !        val = val + A%CICoef(ip)*Sab(ip,ir)*Haa(iq,ip)
 !     enddo
 !     work(iq,ir) = work(iq,ir) + val

 !   enddo
 !enddo

 !!!!(FF|FF):(BA|AA)
 !!call tran4_gen(NBas,&
 !!         NBas,A%CMO,&
 !!         NBas,A%CMO,&
 !!         NBas,B%CMO,&
 !!         NBas,A%CMO,&
 !!         'FFFFBAAA','AOTWOSORT')

 !!open(newunit=iunit,file='FFFFBAAA',status='OLD', &
 !!    access='DIRECT',recl=8*NBas*NBas)

 !!work1 = 0
 !!ints  = 0
 !!kl    = 0
 !!do l=1,NBas 
 !!   do k=1,NBas 
 !!      kl = kl + 1

 !!      if((l.le.dimOA).and.(k.le.dimOA)) then

 !!         read(iunit,rec=kl) work1(1:NBas*NBas)

 !!         ip = l
 !!         iq = k

 !!         do j=1,NBas
 !!            do i=1,NBas
 !!               ints(i,j) = work1((j-1)*NBas+i)
 !!            enddo
 !!         enddo

 !!         do ir=1,NBas
 !!            work(iq,ir) = work(iq,ir) + A%CICoef(ip)*ints(ir,ip)
 !!         enddo

 !!     endif
 !!   enddo
 !!enddo
 !!close(iunit)
 !!call delfile('FFFFBAAA') 
 !!!(FO|OO):(BA|AA)
 !call tran4_gen(NBas,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         NBas,B%CMO,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         'FOOOBAAA','AOTWOSORT')

 !open(newunit=iunit,file='FOOOBAAA',status='OLD', &
 !    access='DIRECT',recl=8*NBas*dimOA)

 !work1 = 0
 !ints  = 0
 !kl    = 0
 !do l=1,dimOA
 !   do k=1,dimOA
 !      kl = kl + 1

 !         read(iunit,rec=kl) work1(1:NBas*dimOA)

 !         ip = l
 !         iq = k

 !         do j=1,NBas
 !            do i=1,NBas
 !               ints(i,j) = work1((j-1)*NBas+i)
 !            enddo
 !         enddo

 !         do ir=1,NBas
 !            work(iq,ir) = work(iq,ir) + A%CICoef(ip)*ints(ir,ip)
 !         enddo

 !   enddo
 !enddo
 !close(iunit)
 !call delfile('FOOOBAAA') 

 !do ir=1,NBas
 !do iq=1,NBas
 !   val = A%ECASSCF*A%CICoef(iq)*Sab(iq,ir)
 !   !print*,'iq,ir',iq,ir,work(iq,ir),val
 !   print*,'iq,ir',iq,ir,work(iq,ir)-val
 !enddo
 !enddo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! test eq 58
 !! MONOMER B
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Hbb = 0; Hab = 0
 call tran2MO(Hb,B%CMO,B%CMO,Hbb,NBas)
 call tran2MO(Hb,A%CMO,B%CMO,Hab,NBas)

 ! test 1-el
 tmp=0
 do i=1,dimOB
    tmp = tmp + 2d0*B%Occ(i)*Hbb(i,i)
 enddo
 print*, 'OneEl(B)',tmp

 !work = 0
 !do ir=1,NBas
 !   do iq=1,dimOB
 !   !do iq=1,NBas 

 !     work(iq,ir) = work(iq,ir) +  B%CICoef(iq)*Hab(iq,ir) 

 !     val = 0
 !     do ip=1,dimOB
 !        val = val + B%CICoef(ip)*Sab(ip,ir)*Hbb(iq,ip)
 !     enddo
 !     work(iq,ir) = work(iq,ir) + val

 !   enddo
 !enddo

!!(FO|OO):(AB|BB)
 !call tran4_gen(NBas,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         NBas,A%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         'FOOOABBB','AOTWOSORT')

 !open(newunit=iunit,file='FOOOABBB',status='OLD', &
 !    access='DIRECT',recl=8*NBas*dimOB)

 !work1 = 0
 !ints  = 0
 !kl    = 0
 !do l=1,dimOB
 !   do k=1,dimOB
 !      kl = kl + 1

 !         read(iunit,rec=kl) work1(1:NBas*dimOB)

 !         ip = l
 !         iq = k

 !         do j=1,NBas
 !            do i=1,NBas
 !               ints(i,j) = work1((j-1)*NBas+i)
 !            enddo
 !         enddo

 !         do ir=1,NBas
 !            work(iq,ir) = work(iq,ir) + B%CICoef(ip)*ints(ir,ip)
 !         enddo

 !   enddo
 !enddo
 !close(iunit)
 !call delfile('FOOOABBB') 

 !do ir=1,NBas
 !do iq=1,NBas
 !   val = B%ECASSCF*B%CICoef(iq)*Sab(iq,ir)
 !   print*,'iq,ir',iq,ir,work(iq,ir),val
 !enddo
 !enddo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! test eq 59
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Haa = 0; Hab = 0
 !call tran2MO(Ha,A%CMO,A%CMO,Haa,NBas)
 !work = 0
 !!(OO|OO):(AA|AA)
 !open(newunit=iunit,file='OOOOAAAA',status='OLD', &
 !    access='DIRECT',recl=8*dimOA*dimOA)

 !work1  = 0
 !ints   = 0
 !kl     = 0
 !do l=1,dimOA
 !   do k=1,dimOA
 !      kl = kl + 1
 !      read(iunit,rec=kl) work1(1:dimOA*dimOA)

 !      ir = l
 !      ip = k

 !      do j=1,dimOA
 !         do i=1,dimOA
 !            ints(i,j) = work1((j-1)*dimOA+i)
 !         enddo
 !      enddo

 !      do iq=1,dimOA
 !         work(iq,ir) = work(iq,ir) + A%CICoef(ip)*ints(iq,ip)
 !      enddo

 !   enddo
 !   ir = l
 !   do iq=1,NBas
 !    work(iq,ir) = work(iq,ir) +  A%CICoef(ir)*Haa(iq,ir) + A%CICoef(iq)*Haa(ir,iq)
 !   enddo

 !enddo
 !close(iunit)

 !do ir=1,Nbas
 !do iq=1,Nbas
 !   print*, 'iq,ir',iq,ir,work(iq,ir)
 !enddo
 !enddo

 !do iq=1,dimOA
 !   print*, 'A-iq',iq,work(iq,iq),A%ECASSCF*A%CICoef(iq)
 !enddo

 !work = 0
 !!(OO|OO):(BB|BB)
 !open(newunit=iunit,file='OOOOBBBB',status='OLD', &
 !    access='DIRECT',recl=8*dimOB*dimOB)

 !work1  = 0
 !ints   = 0
 !kl     = 0
 !do l=1,dimOB
 !   do k=1,dimOB
 !      kl = kl + 1
 !      read(iunit,rec=kl) work1(1:dimOB*dimOB)

 !      ir = l
 !      ip = k

 !      do j=1,dimOB
 !         do i=1,dimOB
 !            ints(i,j) = work1((j-1)*dimOB+i)
 !         enddo
 !      enddo

 !      do iq=1,dimOB
 !         work(iq,ir) = work(iq,ir) + B%CICoef(ip)*ints(iq,ip)
 !      enddo

 !   enddo
 !   ir = l
 !   do iq=1,NBas
 !    work(iq,ir) = work(iq,ir) +  B%CICoef(ir)*Hbb(iq,ir) + B%CICoef(iq)*Hbb(ir,iq)
 !   enddo

 !enddo
 !close(iunit)

 !do iq=1,dimOB
 !   print*, 'B-iq',iq,work(iq,iq),B%ECASSCF*B%CICoef(iq)
 !enddo
 !! end test eq 59

 deallocate(ints)


 deallocate(Hab,Hbb,Hb,Haa,Ha,work1,work)
 deallocate(PBaa,PAbb,Vb,Va,PB,PA,Sab,S)
 deallocate(RDM2Bval,RDM2Aval)

 ! delete ints
 call delfile('OOOOABAB') 
 call delfile('OOOOAABB') 
 call delfile('OOOOABBB') 
 call delfile('OOOOBAAA') 
 call delfile('OOOOAAAA') 
 call delfile('OOOOBBBB') 

end subroutine hl_2el

subroutine e2exind(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT

integer :: NBas,dimOA,dimOB,dimVA,dimVB,nOVA,nOVB
double precision :: nelA,nelB
integer,allocatable :: posA(:,:),posB(:,:)
double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAaa(:,:),PBbb(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:),&
                                Vaba(:,:),Vaab(:,:),Vbab(:,:)
double precision,allocatable :: uA(:),uB(:)
double precision,allocatable :: tindA(:),tindB(:),&
                                tindX(:),VindX(:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: RDM2Aval(:,:,:,:), &
                                RDM2Bval(:,:,:,:)
! unc
type(EBlockData)             :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable :: SBlockA(:),SBlockB(:)
double precision :: termZ_u,termY_u,termX_u
double precision :: e2exi_unc
double precision,allocatable :: uA0(:),uB0(:)
double precision,allocatable :: OmA0(:),OmB0(:)
integer                      :: iblk,nblkA,nblkB
! full test
double precision,allocatable :: AVecX0(:),AVecY0(:), &
                                BVecX0(:),BVecY0(:)

double precision,allocatable :: tmp1(:,:),tmp2(:,:)
double precision,allocatable :: tmpXA(:),tmpYA(:),&
                                tmpXB(:),tmpYB(:)
integer :: i,j,ipq,ip,iq,irs,ir,is
logical :: both
double precision :: termZ,termY,termX 
double precision :: e2exi
double precision :: fact,tmp
! Thresholds
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8 

!avoid tran4 in exch-disp 
 SAPT%noE2exi = .false. 
 both = SAPT%iCpld

! set dimensions
 NBas = A%NBasis 
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

! print thresholds
 if(SAPT%IPrint>5) then 
    write(LOUT,'(/,1x,a)') 'Thresholds in E2exch-ind:'
    write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
    write(LOUT,'(1x,a,2x,e15.4,/)') 'BigE        =', BigE
 endif

 ! unc -- docelowo:
 !if(Flags%ICASSCF==1) then
 !   allocate(OmA0(A%NDimX),OmB0(B%NDimX))

 !   call read_SBlock(SBlockA,SBlockAIV,nblkA,'XY0_A')
 !   call read_SBlock(SBlockB,SBlockBIV,nblkB,'XY0_B')

 !   call unpack_Eig(SBlockA,SBlockAIV,nblkA,OmA0,A%NDimX)
 !   call unpack_Eig(SBlockB,SBlockBIV,nblkB,OmB0,B%NDimX)
 !endif

 ! unc -- full
 if(Flags%ICASSCF==1) then
    allocate(AVecX0(A%NDimX*A%NDimX),OmA0(A%NDimX), &
             AVecY0(A%NDimX*A%NDimX), &
             BVecX0(B%NDimX*B%NDimX),OmB0(B%NDimX), &
             BVecY0(B%NDimX*B%NDimX))

    call unpack_XY0_full(AVecX0,AVecY0,OmA0,A%CICoef,A%IndN,A%NDimX,NBas,'XY0_A')
    call unpack_XY0_full(BVecX0,BVecY0,OmB0,B%CICoef,B%IndN,B%NDimX,NBas,'XY0_B')
 endif

! read EigValA_B
! allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX), &
!          EVecB(B%NDimX,B%NDimX),OmB(B%NDimX))
 allocate(S(NBas,NBas),Sab(NBas,NBas),&
          Sba(NBas,NBas),&
          PA(NBas,NBas),PB(NBas,NBas),&
          PAaa(NBas,NBas),PBbb(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),  &
          Vabb(NBas,NBas),Vbaa(NBas,NBas),&
          Vaba(NBas,NBas),Vaab(NBas,NBas),Vbab(NBas,NBas),&
          WaBB(NBas,NBas),WbAA(NBas,NBas))

 call get_one_mat('S',S,A%Monomer,NBas)
 call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
 call tran2MO(S,B%CMO,A%CMO,Sba,NBas)

 call get_den(NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NBas,B%CMO,B%Occ,1d0,PB)

 allocate(tmp1(NBas,NBas),tmp2(NBas,NBas))

 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,A%CMO,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,B%CMO,NBas,0d0,tmp2,NBas)
 ! PA(B), PB(A)
 call tran2MO(PA,tmp1,tmp1,PAaa,NBas)
 call tran2MO(PB,tmp2,tmp2,PBbb,NBas)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

 call tran2MO(Va,B%CMO,A%CMO,Vaba,NBas)
 call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
 call tran2MO(Vb,A%CMO,B%CMO,Vbab,NBas)

 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)

 allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
          RDM2Bval(dimOB,dimOB,dimOB,dimOB))

 if(Flags%ICASSCF==1) then
    ! CAS
    RDM2Aval = A%RDM2val
    RDM2Bval = B%RDM2val
 elseif(Flags%ICASSCF==0) then
    ! GVB
    RDM2Aval = A%RDM2val
    RDM2Bval = B%RDM2val
 endif

 allocate(posA(NBas,NBas),posB(NBas,NBas))

 if(Flags%ICASSCF==1) then
    posA = 0
    do i=1,A%NDimX
       posA(A%IndN(1,i),A%IndN(2,i)) = A%IndX(i)
    enddo
    posB = 0
    do i=1,B%NDimX
       posB(B%IndN(1,i),B%IndN(2,i)) = B%IndX(i)
    enddo
 elseif(Flags%ICASSCF==0) then
    posA = 0
    do i=1,A%NDimX
       posA(A%IndN(1,i),A%IndN(2,i)) = i
    enddo
    posB = 0
    do i=1,B%NDimX
       posB(B%IndN(1,i),B%IndN(2,i)) = i
    enddo
 endif

 ! termZ
 ! PA(ab).S(bc)
 ! S(ad).PB(dc)
 tmp1 = 0
 tmp2 = 0
 termZ = 0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PB,NBas,0d0,tmp2,NBas)
 do j=1,NBas
    do i=1,NBas
       termZ = termZ + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 if(Flags%ICASSCF==1) termZ_u = 2d0*SAPT%e2ind_unc*termZ
 termZ = 2d0*SAPT%e2ind*termZ

 deallocate(tmp2,tmp1)

!! Y term

! allocate(tindX(A%NDimX))
!
! ! Y(A<--B)
! uA = 0
! tindX = 0
! do i=1,A%NDimX
!    ip = A%IndN(1,i)
!    iq = A%IndN(2,i)
!    ipq = posA(ip,iq)
!
!    tindX(ipq) = tindX(ipq) + (A%Occ(ip)-A%Occ(iq))*WbAA(ip,iq)
!
! enddo
! call dgemv('T',A%NDimX,A%NDimX,1d0,A%EigY-A%EigX,A%NDimX,tindX,1,0d0,uA,1)
! print*, 'uA',norm2(uA)
! !call make_tind2(tindB,B%EigX,B%EigY,Sab,PAaa,B%Occ,B%IndN,posB,dimOA,B%NDimX,NBas)
! !print*, 'tindA-1',norm2(tindA)
! !print*, 'tindB-2',norm2(tindB)
!
! !call make_tind(tindX,A%EigX,A%EigY,Sab,A%Occ,B%Occ,A%IndN,posA,A%NDimX,NBas)
! call make_tind2(tindX,A%EigX,A%EigY,Sab,PBbb,A%Occ,A%IndN,posA,dimOB,A%NDimX,NBas)
! print*, 'tindA',norm2(tindX)
!
! termY=0d0
! do i=1,A%NDimX
!    if(abs(A%Eig(i)).gt.SmallE.and.abs(A%Eig(i)).lt.BigE) then
!       termY = termY + (tindX(i)*uA(i))/A%Eig(i)
!    endif   
! enddo
! print*, 'termY-1',termY
!
! deallocate(tindX)
! allocate(tindX(B%NDimX))
!
! ! Y(A-->B)
! uB = 0
! tindX = 0
! do j=1,B%NDimX
!    ir = B%IndN(1,j)
!    is = B%IndN(2,j)
!    irs = posB(ir,is)
!
!    tindX(irs) = tindX(irs) + (B%Occ(ir)-B%Occ(is))*WaBB(ir,is)
!
! enddo
! call dgemv('T',B%NDimX,B%NDimX,1d0,B%EigY-B%EigX,B%NDimX,tindX,1,0d0,uB,1)
! print*, 'uB',norm2(uB)
!
! !call make_tind(tindX,B%EigX,B%EigY,Sba,B%Occ,A%Occ,B%IndN,posB,B%NDimX,NBas)
! call make_tind2(tindX,A%EigX,A%EigY,Sab,PBbb,A%Occ,A%IndN,posA,dimOB,A%NDimX,NBas)
! print*, 'tindB',norm2(tindX)
!
! do i=1,B%NDimX
!    if(abs(B%Eig(i)).gt.SmallE.and.abs(B%Eig(i)).lt.BigE) then
!       termY = termY + (tindX(i)*uB(i))/B%Eig(i)
!    endif
! enddo
! print*, 'termY-2',termY
!! end split Y

 allocate(tindA(A%NDimX),tindB(B%NDimX))

 ! unc 
 if(Flags%ICASSCF==1) then

    allocate(uA0(A%NDimX),uB0(B%NDimX))
    uA0 = 0
    tindA = 0
    do i=1,A%NDimX
       ip = A%IndN(1,i)
       iq = A%IndN(2,i)
       ipq = posA(ip,iq)

       tindA(ipq) = tindA(ipq) + (A%Occ(ip)-A%Occ(iq))*WbAA(ip,iq)

    enddo
    call dgemv('T',A%NDimX,A%NDimX,1d0,AVecY0-AVecX0,A%NDimX,tindA,1,0d0,uA0,1)
    !print*, 'uA0',norm2(uA0)

    uB0 = 0
    tindB = 0
    do j=1,B%NDimX
       ir = B%IndN(1,j)
       is = B%IndN(2,j)
       irs = posB(ir,is)

       tindB(irs) = tindB(irs) + (B%Occ(ir)-B%Occ(is))*WaBB(ir,is)

    enddo
    call dgemv('T',B%NDimX,B%NDimX,1d0,BVecY0-BVecX0,B%NDimX,tindB,1,0d0,uB0,1)
    !print*, 'uA0',norm2(uB0)

    call make_tind(tindA,AVecX0,AVecY0,Sab,A%Occ,B%Occ,A%IndN,posA,A%NDimX,NBas)
    call make_tind(tindB,BVecX0,BVecY0,Sba,B%Occ,A%Occ,B%IndN,posB,B%NDimX,NBas)

    termY_u = 0d0
    do i=1,A%NDimX
       if(abs(OmA0(i)).gt.SmallE.and.abs(OmA0(i)).lt.BigE) then
          termY_u = termY_u + (tindA(i)*uA0(i))/OmA0(i)
       endif   
    enddo
    do i=1,B%NDimX
       if(abs(OmB0(i)).gt.SmallE.and.abs(OmB0(i)).lt.BigE) then
          termY_u = termY_u + (tindB(i)*uB0(i))/OmB0(i)
       endif
    enddo

    termY_u=-4d0*(SAPT%elst-SAPT%Vnn)*termY_u
    !print*, 'termY_u',termY_u

    ! not working..?
    !allocate(uA0(A%NDimX))
    !   call make_uX_unc(A,WbAA,posA,uA0,A%NDimX,NBas)
    !   print*, 'uA0',norm2(uA0)
    !deallocate(uA0)
 endif

 ! cpld
 if(both) then

    allocate(uA(A%NDimX),uB(B%NDimX))

    uA = 0
    tindA = 0
    do i=1,A%NDimX
       ip = A%IndN(1,i)
       iq = A%IndN(2,i)
       ipq = posA(ip,iq)

       tindA(ipq) = tindA(ipq) + (A%Occ(ip)-A%Occ(iq))*WbAA(ip,iq)

    enddo
    call dgemv('T',A%NDimX,A%NDimX,1d0,A%EigY-A%EigX,A%NDimX,tindA,1,0d0,uA,1)
    !print*, 'uA',norm2(uA)

    uB = 0
    tindB = 0
    do j=1,B%NDimX
       ir = B%IndN(1,j)
       is = B%IndN(2,j)
       irs = posB(ir,is)

       tindB(irs) = tindB(irs) + (B%Occ(ir)-B%Occ(is))*WaBB(ir,is)

    enddo
    call dgemv('T',B%NDimX,B%NDimX,1d0,B%EigY-B%EigX,B%NDimX,tindB,1,0d0,uB,1)
    !print*, 'uB',norm2(uB)

     call make_tind(tindA,A%EigX,A%EigY,Sab,A%Occ,B%Occ,A%IndN,posA,A%NDimX,NBas)
     call make_tind(tindB,B%EigX,B%EigY,Sba,B%Occ,A%Occ,B%IndN,posB,B%NDimX,NBas)
    
    !call make_tind2(tindA,A%EigX,A%EigY,Sab,PBbb,A%Occ,A%IndN,posA,dimOB,A%NDimX,NBas)
    !call make_tind2(tindB,B%EigX,B%EigY,Sba,PAaa,B%Occ,B%IndN,posB,dimOA,B%NDimX,NBas)

    termY=0d0
    do i=1,A%NDimX
       if(abs(A%Eig(i)).gt.SmallE.and.abs(A%Eig(i)).lt.BigE) then
          termY = termY + (tindA(i)*uA(i))/A%Eig(i)
       endif   
    enddo
    do i=1,B%NDimX
       if(abs(B%Eig(i)).gt.SmallE.and.abs(B%Eig(i)).lt.BigE) then
          termY = termY + (tindB(i)*uB(i))/B%Eig(i)
       endif
    enddo

    termY=-4d0*(SAPT%elst-SAPT%Vnn)*termY
    !write(*,*) 'termY',termY

 endif
 deallocate(tindB,tindA)

 ! term A3
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'FOFOAABB','AOTWOSORT')
! test
! call tran4_gen(NBas,&
!          NBas,A%CMO,&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          NBas,B%CMO,&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          'FOFOBBAA','AOTWOSORT')
 ! term A1
 call tran4_gen(NBas,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          NBas,B%CMO,&
          'FFOOABAB','AOTWOSORT')
 !! test-2
 !call tran4_gen(NBas,&
 !         NBas,A%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         NBas,A%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         'FOFOABAB','AOTWOSORT')
 !! test-3
 !call tran4_gen(NBas,&
 !         NBas,B%CMO,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         NBas,A%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         'FOFOABBA','AOTWOSORT')
 ! term A2
 ! A2A(B): XX
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,B%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'FOFOBBBA','AOTWOSORT')
 call tran4_gen(NBas,&
          NBas,A%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'FOFOAAAB','AOTWOSORT')
 !! A2A(B): YY
 call tran4_gen(NBas,&
          NBas,A%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,B%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'FOFOBBAB','AOTWOSORT')
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,A%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'FOFOAABA','AOTWOSORT')

 ! A3 prototype
 allocate(tmpXA(A%NDimX),tmpYA(A%NDimX),&
          tmpXB(B%NDimX),tmpYB(B%NDimX))
 tmpXA = 0
 tmpXB = 0
 tmpYA = 0
 tmpYB = 0
 !call exind_A3_XY(A%NDimX,A%NDimX,tmpXA,tmpYA,RDM2Aval,RDM2Bval,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
 !                 B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)
 !
 !! test B
 !print*, 'test A3-B!'
 !call exind_A3_XY(B%NDimX,B%NDimX,tmpXB,tmpYB,RDM2Bval,RDM2Aval,Sba,nelB,Vbaa,nelA,Vabb,'FOFOBBAA',&
 !                 A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas)

 ! A3 here!!
 call exind_A3_XY_full(A%NDimX,B%NDimX,tmpXA,tmpXB,RDM2Aval,RDM2Bval,Sab, &
                    nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                    B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 tmpYA = -tmpXA
 tmpYB = -tmpXB

! A1
 call exind_A1_AB(A%NDimX,B%NDimX,tmpXA,tmpXB,tmpYA,tmpYB,nelA,nelB,dimOA,dimOB, &
                  A%Occ,B%Occ,PAaa,PBbb,Vaab,Vbab,Sab,A%IndN,B%IndN,posA,posB,NBas)
! A2
 call exind_A2_XX(A%NDimX,B%NDimX,tmpXA,tmpXB,RDM2Bval, &
                  Sab,nelA,Vabb,nelB,Vbab,PAaa,'FOFOBBBA', &
                  B%Occ,A%Occ,B%IndN,A%IndN,posB,posA,  &
                  dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
 call exind_A2_XX(B%NDimX,A%NDimX,tmpXB,tmpXA,RDM2Aval, &
                 Sba,nelB,Vbaa,nelA,Vaba,PBbb,'FOFOAAAB', &
                  A%Occ,B%Occ,A%IndN,B%IndN,posA,posB,  &
                  dimOA,dimOB,A%NDimX,B%NDimX,NBas,.false.)

 !print*, 'XA-tot',norm2(tmpXA)
 !print*, 'XB-tot',norm2(tmpXB)

 call exind_A2_YY(A%NDimX,B%NDimX,tmpYA,tmpYB,RDM2Bval, &
                  Sab,nelA,Vabb,nelB,Vbab,PAaa,'FOFOBBAB', &
                  B%Occ,A%Occ,B%IndN,A%IndN,posB,posA,  &
                  dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)

 !print*, 'YA-tot-1',norm2(tmpYA)
 !print*, 'YB-tot-1',norm2(tmpYB)

 call exind_A2_YY(B%NDimX,A%NDimX,tmpYB,tmpYA,RDM2Aval, &
                  Sba,nelB,Vbaa,nelA,Vaba,PBbb,'FOFOAABA', &
                  A%Occ,B%Occ,A%IndN,B%IndN,posA,posB,  &
                  dimOA,dimOB,A%NDimX,B%NDimX,NBas,.false.)

 !print*, 'YA-tot-2',norm2(tmpYA)
 !print*, 'YB-tot-2',norm2(tmpYB)

 ! unc
 if(Flags%ICASSCF==1) then

    allocate(VindX(A%NDimX))

    call dgemv('T',A%NDimX,A%NDimX,1d0,AVecX0,A%NDimX,tmpXA,1,0d0,VindX,1)
    call dgemv('T',A%NDimX,A%NDimX,1d0,AVecY0,A%NDimX,tmpYA,1,1d0,VindX,1)

    termX_u = 0
    do i=1,A%NDimX
       if(abs(OmA0(i)).gt.SmallE.and.abs(OmA0(i)).lt.BigE) then
          termX_u = termX_u + (VindX(i)*uA0(i))/OmA0(i)
       endif   
    enddo
    !print*, 'termX_unc-1',termX_u

    deallocate(VindX)
    allocate(VindX(B%NDimX))

    call dgemv('T',B%NDimX,B%NDimX,1d0,BVecX0,B%NDimX,tmpXB,1,0d0,VindX,1)
    call dgemv('T',B%NDimX,B%NDimX,1d0,BVecY0,B%NDimX,tmpYB,1,1d0,VindX,1)

    do i=1,B%NDimX
       if(abs(OmB0(i)).gt.SmallE.and.abs(OmB0(i)).lt.BigE) then
          termX_u = termX_u + (VindX(i)*uB0(i))/OmB0(i)
       endif
    enddo
    !print*, 'termX_unc-2',termX_u
    termX_u = -2d0*termX_u

    e2exi_unc = termX_u + termY_u + termZ_u

    SAPT%e2exind_unc = e2exi_unc
    write(LOUT,'(/1x,a,f16.8)') 'E2exch-ind(unc) =', e2exi_unc*1.0d3

    deallocate(VindX)
    deallocate(uB0,uA0)

 endif

 ! cpld
 if(both) then
    allocate(VindX(A%NDimX))

    call dgemv('T',A%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmpXA,1,0d0,VindX,1)
    call dgemv('T',A%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmpYA,1,1d0,VindX,1)

    termX = 0
    do i=1,A%NDimX
       if(abs(A%Eig(i)).gt.SmallE.and.abs(A%Eig(i)).lt.BigE) then
          termX = termX + (VindX(i)*uA(i))/A%Eig(i)
       endif   
    enddo
    !print*, 'termX-1',termX

    deallocate(VindX)
    allocate(VindX(B%NDimX))

    VindX=0
    call dgemv('T',B%NDimX,B%NDimX,1d0,B%EigX,B%NDimX,tmpXB,1,0d0,VindX,1)
    call dgemv('T',B%NDimX,B%NDimX,1d0,B%EigY,B%NDimX,tmpYB,1,1d0,VindX,1)

    do i=1,B%NDimX
       if(abs(B%Eig(i)).gt.SmallE.and.abs(B%Eig(i)).lt.BigE) then
          termX = termX + (VindX(i)*uB(i))/B%Eig(i)
       endif
    enddo
    !print*, 'termX-2',termX
    termX=-2d0*termX

    !print*, 'termX',termX
    if(SAPT%IPrint>2) write(LOUT,'(/1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3
    if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)')  'term Y      = ',  termY*1.0d3
    if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)')  'term X      = ',  termX*1.0d3

    e2exi = termX + termY + termZ

    SAPT%e2exind = e2exi
    write(LOUT,'(/1x,a,f16.8)') 'E2exch-ind  = ', e2exi*1.0d3

    deallocate(VindX)
    deallocate(uB,uA)

 endif

! deallocate(VindB,VindA)
 deallocate(tmpYB,tmpXB)
 deallocate(tmpYA,tmpXA)

 deallocate(posB,posA)
 deallocate(RDM2Bval,RDM2Aval)
 deallocate(PBbb,PAaa)
 deallocate(Vbab,Vaba,Vaab,Vbaa,Vb,Va)
 deallocate(WbAA,WaBB,PB,PA,Sba,Sab,S)

end subroutine e2exind

subroutine e2exdisp(Flags,A,B,SAPT)

use exappr 
use sref_exch 

implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT

integer          :: NBas,NInte1
integer          :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer          :: GdimOA,GdimOB
integer          :: iunit
integer          :: i,j,k,l,kl,ij,ii,jj,pq,rs
integer          :: ip,iq,ir,is,ipq,irs
integer          :: iblk,nblkA,nblkB
logical          :: both,approx,ipropab 
double precision :: nelA,nelB
double precision :: fact,val,termZ,termY,termX
double precision :: termZ_u,termY_u,termX_u
double precision :: e2exd_u,e2exd
integer,allocatable          :: posA(:,:),posB(:,:)
double precision,allocatable :: ints(:,:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),tmp3(:,:),&
                                sij(:,:),  &
                                tmp_u(:,:),&
                                work(:),workSq(:,:)
double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:),&
                                Vaba(:,:),Vaab(:,:),Vbab(:,:)
double precision,allocatable :: RDM2Aval(:,:,:,:), &
                                RDM2Bval(:,:,:,:)
type(EBlockData)             :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable :: SBlockA(:),SBlockB(:)
! uncoupled full test
double precision,allocatable :: AVecX0(:),AVecY0(:), &
                                BVecX0(:),BVecY0(:), &
                                OmA0(:),OmB0(:)
double precision :: Tcpu,Twall
double precision,external  :: trace,FRDM2,FRDM2GVB
! test for Be
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: BigE   = 1.D8 
double precision,parameter :: SmallE = 1.D-3

! set dimensions
 NBas = A%NBasis 
 !dimOA = A%INAct+A%NAct
 dimOA = A%num0+A%num1
 GdimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 !dimOB = B%INAct+B%NAct
 dimOB = B%num0+B%num1
 GdimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

! print thresholds
 if(SAPT%IPrint>5) then 
    write(LOUT,'(/,1x,a)') 'Thresholds in E2exch-disp:'
    write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
    write(LOUT,'(1x,a,2x,e15.4,/)') 'BigE        =', BigE
 endif

! set e2exd_version
 !approx=.true.
 approx = .false.
 both   = SAPT%iCpld

! approximate RDM2
 if(approx) then
    allocate(A%Fmat(NBas,NBas),B%Fmat(NBas,NBas))
    call fill_Fmat(A%Fmat,A%Occ,NBas,1)
    call fill_Fmat(B%Fmat,B%Occ,NBas,1)
 endif

 ! uncoupled - works for CAS only!
 if(Flags%ICASSCF==1) then
    allocate(OmA0(A%NDimX),OmB0(B%NDimX))

    call read_SBlock(SBlockA,SBlockAIV,nblkA,'XY0_A')
    call read_SBlock(SBlockB,SBlockBIV,nblkB,'XY0_B')

    call unpack_Eig(SBlockA,SBlockAIV,nblkA,OmA0,A%NDimX)
    call unpack_Eig(SBlockB,SBlockBIV,nblkB,OmB0,B%NDimX)
 endif
 !! uncoupled-ver0
 !allocate(AVecX0(A%NDimX*A%NDimX),OmA0(A%NDimX), &
 !         AVecY0(A%NDimX*A%NDimX), &
 !         BVecX0(B%NDimX*B%NDimX),OmB0(B%NDimX), &
 !         BVecY0(B%NDimX*B%NDimX))
 !
 !call unpack_XY0_full(AVecX0,AVecY0,OmA0,A%CICoef,A%IndN,A%NDimX,NBas,'XY0_A') 
 !call unpack_XY0_full(BVecX0,BVecY0,OmB0,B%CICoef,B%IndN,B%NDimX,NBas,'XY0_B') 

 allocate(S(NBas,NBas),&
          Sab(NBas,NBas),Sba(NBas,NBas),&
          PA(NBas,NBas),PB(NBas,NBas),  &
          Va(NBas,NBas),Vb(NBas,NBas),  &
          Vabb(NBas,NBas),Vbaa(NBas,NBas),&
          Vaba(NBas,NBas),Vaab(NBas,NBas),Vbab(NBas,NBas))
 allocate(tmp1(NBas,NBas),tmp2(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NBas,B%CMO,B%Occ,1d0,PB)

 call get_one_mat('S',S,A%Monomer,NBas)
 call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
 call tran2MO(S,B%CMO,A%CMO,Sba,NBas)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

 call tran2MO(Va,B%CMO,A%CMO,Vaba,NBas)
 call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
 call tran2MO(Vb,A%CMO,B%CMO,Vbab,NBas)

 allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
          RDM2Bval(dimOB,dimOB,dimOB,dimOB))

 if(Flags%ICASSCF==1) then
    ! CAS
    RDM2Aval = A%RDM2val
    RDM2Bval = B%RDM2val
 elseif(Flags%ICASSCF==0) then
    ! GVB
    RDM2Aval = A%RDM2val
    RDM2Bval = B%RDM2val
 endif

 allocate(posA(NBas,NBas),posB(NBas,NBas))
 if(Flags%ICASSCF==1) then
    posA = 0
    do i=1,A%NDimX
       posA(A%IndN(1,i),A%IndN(2,i)) = A%IndX(i)
    enddo
    posB = 0
    do i=1,B%NDimX
       posB(B%IndN(1,i),B%IndN(2,i)) = B%IndX(i)
    enddo
 elseif(Flags%ICASSCF==0) then
    posA = 0
    do i=1,A%NDimX
       posA(A%IndN(1,i),A%IndN(2,i)) = i
    enddo
    posB = 0
    do i=1,B%NDimX
       posB(B%IndN(1,i),B%IndN(2,i)) = i
    enddo
 endif

 ! termZ
 ! PA(ab).S(bc)
 ! S(ad).PB(dc)
 tmp1 = 0
 tmp2 = 0
 termZ = 0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PB,NBas,0d0,tmp2,NBas)
 do j=1,NBas
    do i=1,NBas
       termZ = termZ + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 termZ_u = termZ
 termZ = 2d0*SAPT%e2disp*termZ
 termZ_u = 2d0*SAPT%e2disp_unc*termZ_u
 !write(LOUT,*) 'termZ ',termZ
 !write(LOUT,'(1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3

 deallocate(tmp2,tmp1)

 allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX),&
          tmp_u(A%NDimX,B%NDimX),sij(A%NDimX,B%NDimX),work(nOVB),workSq(NBas,NBas))
 if(both) allocate(tmp3(A%NDimX,B%NDimX))

 ! term A1
 ! transform J and K
 if(SAPT%noE2exi) then
    call tran4_gen(NBas,&
             A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
             B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
             NBas,A%CMO,&
             NBas,B%CMO,&
             'FFOOABAB','AOTWOSORT')
    endif
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,A%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'FOFOABBA','AOTWOSORT')
 ! term A2
 ! A2A(B): XX
 if(SAPT%noE2exi) then
    call tran4_gen(NBas,&
             NBas,B%CMO,&
             A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
             NBas,B%CMO,&
             B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
             'FOFOBBBA','AOTWOSORT')
    call tran4_gen(NBas,&
             NBas,A%CMO,&
             B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
             NBas,A%CMO,&
             A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
             'FOFOAAAB','AOTWOSORT')
    ! A2A(B): YY
    call tran4_gen(NBas,&
             NBas,A%CMO,&
             B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
             NBas,B%CMO,&
             B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
             'FOFOBBAB','AOTWOSORT')
    call tran4_gen(NBas,&
             NBas,B%CMO,&
             A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
             NBas,A%CMO,&
             A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
             'FOFOAABA','AOTWOSORT')
 endif
 ! term A3
 if(SAPT%noE2exi) then
    call tran4_gen(NBas,&
             NBas,B%CMO,&
             B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
             NBas,A%CMO,&
             A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
             'FOFOAABB','AOTWOSORT')
 endif
 ! XY and YX, A2
 call tran4_gen(NBas,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          NBas,B%CMO,&
          'FFOOABBB','AOTWOSORT')
 call tran4_gen(NBas,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,B%CMO,&
          NBas,A%CMO,&
          'FFOOBAAA','AOTWOSORT')

 deallocate(work)
 allocate(work(NBas**2),ints(NBas,NBas))

 ! A3: XX
 tmp1=0
 if(approx) then

    call app_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                   B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 else

    call inter_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,RDM2Aval,RDM2Bval,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                     B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 endif

 ! TERMS XX, YY
 !(FO|FO):(AB|BA)
 open(newunit=iunit,file='FOFOABBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

 ints = 0
 kl = 0
 do l=1,dimOA
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       ir = k
       iq = l

       do ip=1,NBas
          do is=1,dimOB

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then

                do j=1,dimOB
                   do i=1,NBas
                      ints(i,j) = work((j-1)*NBas+i)
                   enddo
                enddo

                fact = -2d0*(A%Occ(ip)-A%Occ(iq)) * &
                            (B%Occ(ir)-B%Occ(is)) * &
                            ints(ip,is)

                tmp1(ipq,irs) = tmp1(ipq,irs) + fact
 
             endif
          enddo
       enddo
 
    enddo
 enddo
 close(iunit)

 sij = tmp1

 ! A1:XX
! call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'XX')
 call A1_Mod_1el(tmp1,A%XELE,B%XELE,A%Occ,B%Occ,Vaab,Vbab,Sab,&
                 A%IndN,B%IndN,posA,posB,A%NDimX,B%NDimX,NBas,'XX')

 ! A2:XX
 if(approx) then

    call app_A2_XX(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat, &
                   Sab,nelA,Vabb,nelB,Vbab,'FOFOBBBA',&
                   A%Occ,B%IndN,A%IndN,posB,posA,     &
                   dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_XX(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat, & 
                   Sba,nelB,Vbaa,nelA,Vaba,'FOFOAAAB',&
                   B%Occ,A%IndN,B%IndN,posA,posB,     &
                   dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_XX(A%NDimX,B%NDimX,tmp1,RDM2Bval,     &
                     Sab,nelA,Vabb,nelB,Vbab,'FOFOBBBA',&
                     A%Occ,B%IndN,A%IndN,posB,posA,     &
                     dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_XX(A%NDimX,B%NDimX,tmp1,RDM2Aval,     &
                     Sba,nelB,Vbaa,nelA,Vaba,'FOFOAAAB',&
                     B%Occ,A%IndN,B%IndN,posA,posB,     &
                     dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif
 
 !X_A.I.X_B
 if(both) then
    call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
   ! Agnieszka-test
   ! tmp2=0
   ! call mnoz_jeden(NBas,A%NDimX,B%NDimX,A%EigX,tmp1,tmp2,A,B)
    call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,0d0,tmp3,A%NDimX)
 endif

 if(Flags%ICASSCF==1) then
    ! UNC
    tmp_u = 0
    call abpm_tran_gen(tmp1,tmp_u,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                       nblkA,nblkB,A%NDimX,B%NDimX,'XX')
    !
    ! UNC-test-full
    !!X(0)_A.I.X(0)_B
    !call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,AVecX0,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    !call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,BVecX0,B%NDimX,0d0,tmp_u,A%NDimX)
 endif

  tmp1 = sij

 ! A1:YY
 ! call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'YY')
 call A1_Mod_1el(tmp1,A%XELE,B%XELE,A%Occ,B%Occ,Vaab,Vbab,Sab,&
                 A%IndN,B%IndN,posA,posB,A%NDimX,B%NDimX,NBas,'YY')

 ! A2:YY
 if(approx) then

    call app_A2_YY(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_YY(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_YY(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
             A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_YY(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif

 !Y_A.I.Y_B
 if(both) then
    call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
   ! Agnieszka-test
   ! tmp2=0
   ! call mnoz_jeden(NBas,A%NDimX,B%NDimX,A%EigY,tmp1,tmp2,A,B)
    call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
 endif

 if(Flags%ICASSCF==1) then
    ! UNC
    ! Y(0)_A.I.Y(0)_B
    call abpm_tran_gen(tmp1,tmp_u,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                       nblkA,nblkB,A%NDimX,B%NDimX,'YY')

    ! UNCOUPLED-test-full
    ! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,AVecY0,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    ! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,BVecY0,B%NDimX,1d0,tmp_u,A%NDimX)
    !
 endif

! TERMS XY, YX
! A3:XY
 tmp1 = 0
 if(approx) then

    call app_A3_XY(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                   B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 else

   call inter_A3_XY(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,RDM2Aval,RDM2Bval,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                  B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 endif

! (FF|OO):(AB|AB)
 open(newunit=iunit,file='FFOOABAB',status='OLD', &
     access='DIRECT',recl=8*NBas*NBas)

 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,dimOA
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*NBas)

       iq = k
       is = l

       do j=1,NBas
          do i=1,NBas
             ints(i,j) = work((j-1)*NBas+i)
          enddo
       enddo

       do ir=1,NBas
          do ip=1,NBas

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then


                fact = 2d0*(A%Occ(ip)-A%Occ(iq)) * &
                           (B%Occ(ir)-B%Occ(is)) * &
                           ints(ip,ir)

                tmp1(ipq,irs) = tmp1(ipq,irs) + fact
 
             endif
            enddo
         enddo
 
    enddo
 enddo
 close(iunit)

 sij = tmp1

 ! A1:XY
! call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'XY')
 call A1_Mod_1el(tmp1,A%XELE,B%XELE,A%Occ,B%Occ,Vaab,Vbab,Sab,&
                 A%IndN,B%IndN,posA,posB,A%NDimX,B%NDimX,NBas,'XY')

! A2:XY
 if(approx) then

    call app_A2_XY(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_YX(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA','FFOOBAAA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_XY(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_YX(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA','FFOOBAAA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.,B)

 endif

 !X_A.I.Y_B
 if(both) then
    call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    ! Agnieszka-test
    ! tmp2=0
    ! call mnoz_jeden(NBas,A%NDimX,B%NDimX,A%EigX,tmp1,tmp2,A,B)
    call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
 endif

 if(Flags%ICASSCF==1) then
    ! UNC
    !X(0)_A.I.Y(0)_B
    call abpm_tran_gen(tmp1,tmp_u,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                       nblkA,nblkB,A%NDimX,B%NDimX,'XY')

    !! UNCOUPLED-test-full
    !call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,AVecX0,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    !call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,BVecY0,B%NDimX,1d0,tmp_u,A%NDimX)
 endif

 tmp1 = sij

 ! A1:YX
! call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'YX')
 call A1_Mod_1el(tmp1,A%XELE,B%XELE,A%Occ,B%Occ,Vaab,Vbab,Sab, &
                 A%IndN,B%IndN,posA,posB,A%NDimX,B%NDimX,NBas,'YX')

 !A2: YX 
 if(approx) then

    call app_A2_YX(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB','FFOOABBB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_XY(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_YX(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB','FFOOABBB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.,B)
    call inter_A2_XY(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif

 !Y_A.I.X_B
 if(both) then
    call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    !! Agnieszka-test
    ! tmp2=0
    ! call mnoz_jeden(NBas,A%NDimX,B%NDimX,A%EigY,tmp1,tmp2,A,B)
    call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,1d0,tmp3,A%NDimX)
 endif

 if(Flags%ICASSCF==1) then
    ! UNC
    call abpm_tran_gen(tmp1,tmp_u,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                       nblkA,nblkB,A%NDimX,B%NDimX,'YX')

    ! UNC-test-full
    !call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,AVecY0,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    !call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,BVecX0,B%NDimX,1d0,tmp_u,A%NDimX)
 endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print*, 'sij',norm2(sij)
! print*, 'tmp3-exd',norm2(tmp3)

!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!! SINGLE REFERENCE CODE !!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THERE IS SOME BUG IN A2_XY...
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call exdisp_sref(SAPT%elst,SAPT%Vnn,termZ,sij,Sab,Sba,&
!                 Vaab,Vaba,Vbab,Vabb,Vbaa,posA,posB,NBas,A,B)


 if(Flags%ICASSCF==1) then

    ! UNC
    sij = 0
    inquire(file='PROP_AB0',EXIST=ipropab)
    if(ipropab) then
       ! read s_ij
       open(newunit=iunit,file='PROP_AB0',form='UNFORMATTED',&
          access='SEQUENTIAL',status='OLD')
       read(iunit) sij
       close(iunit)

    else
    
       ! make s_ij uncoupled
       !call make_sij_Yunc(sij,tmp1,A%Occ,B%Occ,A%EigY,A%EigX,B%EigY,B%EigX,&
       !                A%num0,B%num0,dimOA,dimOB,nOVB,A%IndN,B%IndN,A%NDimX,B%NDimX,NBas)
       write(LOUT,'(1x,a)') 'ERROR! Zrob make_sij(unc)!!!'
       stop

    endif

    ! term X(unc)
    termX_u = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX
           
          if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
             .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then

             termX_u = termX_u + tmp_u(i,j)*sij(i,j)/(OmA0(i)+OmB0(j))

          endif
       enddo
    enddo
    termX_u = -4d0*termX_u

    call make_tij_Y_unc(tmp_u,tmp1,A%Occ,B%Occ,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                        A%IndN,B%IndN,posA,posB,Sab,Sba,nblkA,nblkB,A%NDimX,B%NDimX,NBas)

    ! term Y(unc)
    termY_u = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX
           
          if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
             .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then

             termY_u = termY_u + sij(i,j)*tmp_u(i,j)/(OmA0(i)+OmB0(j))

          endif
       enddo
    enddo

    termY_u = -8d0*(SAPT%elst-SAPT%Vnn)*termY_u

    if(SAPT%IPrint>2) write(LOUT,'(/1x,a,f16.8)') 'term Z(unc) = ',  termZ_u*1.0d3
    if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)')  'term X(unc) = ',  termX_u*1.0d3
    if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)')  'term Y(unc) = ',  termY_u*1.0d3
    
    e2exd_u = termX_u + termY_u + termZ_u
    SAPT%e2exdisp_unc = e2exd_u
    write(LOUT,'(/1x,a,f11.8)') 'E2exch-disp(unc) = ', e2exd_u*1.0d3

    ! deallocate SBLOCK
    do iblk=1,nblkA
       associate(A => SBlockA(iblk))
         deallocate(A%matY,A%matX,A%vec)
         deallocate(A%pos)
       end associate
    enddo
    do iblk=1,nblkB
       associate(B => SBlockB(iblk))
         deallocate(B%matY,B%matX,B%vec)
         deallocate(B%pos)
       end associate
    enddo
    ! deallocate IV part 
    associate(A => SBlockAIV)
      if(A%n>0) then
         deallocate(A%vec)
         deallocate(A%pos)
      endif
    end associate
    associate(B => SBlockBIV)
      if(B%n>0) then
         deallocate(B%vec)
         deallocate(B%pos)
      endif
    end associate

    deallocate(OmA0,OmB0)
    !deallocate(AVecX0,AVecY0,BVecX0,BVecY0)

 ! end UNC (for CAS only)
 endif

 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(both) then
    ! coupled
    sij = 0
    inquire(file='PROP_AB',EXIST=ipropab)
    if(ipropab) then
       ! read s_ij
       open(newunit=iunit,file='PROP_AB',form='UNFORMATTED',&
          access='SEQUENTIAL',status='OLD')
       read(iunit) sij
       close(iunit)

    else
    
       ! make s_ij
       call make_sij_Y(sij,tmp1,A%Occ,B%Occ,A%EigY,A%EigX,B%EigY,B%EigX,&
                       A%num0,B%num0,dimOA,dimOB,nOVB,A%IndN,B%IndN,A%NDimX,B%NDimX,NBas)

    endif

    ! term X
    termX = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX
           
      !    if(A%Eig(i).gt.SmallE.and.B%Eig(j).gt.SmallE&
      !       .and.A%Eig(i).lt.BigE.and.B%Eig(j).lt.BigE) then

          if(abs(A%Eig(i)).gt.SmallE.and.abs(B%Eig(j)).gt.SmallE&
             .and.abs(A%Eig(i)).lt.BigE.and.abs(B%Eig(j)).lt.BigE) then


             termX = termX + tmp3(i,j)*sij(i,j)/(A%Eig(i)+B%Eig(j))

          endif
       enddo
    enddo
    termX = -4d0*termX

    if(SAPT%IPrint>2) write(LOUT,'(/1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3
    !write(LOUT,*) 'termX ', termX
    if(SAPT%IPrint>2)write(LOUT,'(1x,a,f16.8)') 'term X      = ',  termX*1.0d3

    ! termY
    ! term Y: t_ij 
    !call make_tij_Y(tmp3,tmp2,tmp1,posA,posB,Sab,Sba,A,B,NBas)
    call make_tij_Y(tmp3,tmp2,tmp1,A%Occ,B%Occ,A%EigY,A%EigX,B%EigY,B%EigX,&
                    A%IndN,B%IndN,posA,posB,Sab,Sba,A%NDimX,B%NDimX,NBas)

    ! term Y
    termY = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX
           
         ! if(A%Eig(i).gt.SmallE.and.B%Eig(j).gt.SmallE&
         !    .and.A%Eig(i).lt.BigE.and.B%Eig(j).lt.BigE) then

          if(abs(A%Eig(i)).gt.SmallE.and.abs(B%Eig(j)).gt.SmallE&
             .and.abs(A%Eig(i)).lt.BigE.and.abs(B%Eig(j)).lt.BigE) then

             !termY = termY + 1d0/(A%Eig(i)+B%Eig(j))
             termY = termY + sij(i,j)*tmp3(i,j)/(A%Eig(i)+B%Eig(j))

          endif
       enddo
    enddo

    termY = -8d0*(SAPT%elst-SAPT%Vnn)*termY
    if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)') 'term Y      = ',  termY*1.0d3

    e2exd = termX + termY + termZ
    SAPT%e2exdisp     = e2exd
    write(LOUT,'(/1x,a,f16.8)') 'E2exch-disp = ', e2exd*1.0d3

 ! end coupled
 endif
   
 deallocate(sij)

 deallocate(posB,posA,ints)
 deallocate(Vbab,Vaba,Vaab,Vbaa,Vabb,Vb,Va,PB,PA,Sab,S)
 deallocate(workSq,work,tmp_u,tmp2,tmp1)
 if(both) deallocate(tmp3)

end subroutine e2exdisp

subroutine ModABMin(Occ,SRKer,Wt,OrbGrid,TwoNO,TwoElErf,ABMin,IndN,IndX,NDimX,NGrid,NInte2,NBasis)
!     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid,NInte2
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX)
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(in) :: TwoNO(NInte2),TwoElErf(NInte2)
double precision,intent(inout) :: ABMin(NDimX,NDimX)

double precision :: CICoef(NBasis)
double precision,allocatable :: work(:),batch(:,:),ABKer(:,:)
integer :: i,j,IRow,ICol,ia,ib,iab,ic,id,icd
integer :: offset,batchlen
double precision :: XKer1234,TwoSR,CA,CB,CC,CD
integer,external :: NAddr3,NAddrrK

allocate(work(maxlen),batch(maxlen,NBasis),ABKer(NDimX,NDimX))

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

ABKer = 0

do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   work(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,NDimX
      ia = IndN(1,IRow)
      ib = IndN(2,IRow)
      iab = IndX(IRow)
   
      do ICol=1,NDimX
         ic=IndN(1,ICol)
         id=IndN(2,ICol)
         icd=IndX(ICol)
         if(icd.gt.iab) cycle
    
         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + work(i)* &
            batch(i,ia)*batch(i,ib)*batch(i,ic)*batch(i,id)
         enddo
         
         ABKer(iab,icd) = ABKer(iab,icd) + XKer1234
         ABKer(icd,iab) = ABKer(iab,icd)
      
      enddo
   enddo

enddo

do IRow=1,NDimX
   ia = IndN(1,IRow)
   ib = IndN(2,IRow)
   iab = IndX(IRow)
   CA = CICoef(ia)
   CB = CICoef(ib)

   do ICol=1,NDimX
      ic=IndN(1,ICol)
      id=IndN(2,ICol)
      icd=IndX(ICol)
      CC=CICoef(ic)
      CD=CICoef(id)
      !if(icd.gt.iab) cycle
 
      TwoSR=TwoNO(NAddr3(ia,ib,ic,id))-TwoElErf(NAddr3(ia,ib,ic,id))
      
      ABMin(iab,icd) = ABMin(iab,icd) &
                       +4.0d0*(CA+CB)*(CD+CC)*(ABKer(iab,icd)+TwoSR)
      !ABMin(icd,iab) = ABMin(iab,icd)   

   enddo
enddo

deallocate(ABKer,batch,work)

end subroutine ModABMin

subroutine ModABMin_mithap(Occ,SRKer,Wt,OrbGrid,ABMin,IndN,IndX,NDimX,NGrid,NBasis,&
                           twofile,twoerfile)
! ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX)
character(*),intent(in) :: twofile,twoerfile
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(inout) :: ABMin(NDimX,NDimX)

integer :: offset,batchlen,iunit1,iunit2
integer :: i,j,k,l,kl,ip,iq,ir,is,irs,ipq,igrd
integer :: IRow,ICol
double precision :: XKer1234,TwoSR,Cpq,Crs
integer :: pos(NBasis,NBasis)
double precision :: CICoef(NBasis)
double precision,allocatable :: work1(:),work2(:),WtKer(:)
double precision,allocatable :: batch(:,:),ABKer(:,:)
double precision,allocatable :: ints1(:,:),ints2(:,:)

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints1(NBasis,NBasis),ints2(NBasis,NBasis),&
         WtKer(maxlen),batch(maxlen,NBasis),ABKer(NDimX,NDimX))

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

ABKer = 0

!print*, 'ModABMin_mithap'

do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   WtKer(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,NDimX
      ip = IndN(1,IRow)
      iq = IndN(2,IRow)
      ipq = IndX(IRow)
   
      do ICol=1,NDimX
         ir=IndN(1,ICol)
         is=IndN(2,ICol)
         irs=IndX(ICol)
         if(irs.gt.ipq) cycle
    
         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + WtKer(i)* &
            batch(i,ip)*batch(i,iq)*batch(i,ir)*batch(i,is)
         enddo
         
         ABKer(ipq,irs) = ABKer(ipq,irs) + XKer1234
         ABKer(irs,ipq) = ABKer(ipq,irs)
      
      enddo
   enddo

enddo

open(newunit=iunit1,file=trim(twofile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
open(newunit=iunit2,file=trim(twoerfile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1   
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        read(iunit1,rec=kl) work1(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work1,ints1,NBasis)
        read(iunit2,rec=kl) work2(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work2,ints2,NBasis)

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)
                !if(irs.gt.ipq) cycle

                TwoSR = ints1(i,j)-ints2(i,j)

                ABMIN(ipq,irs) = ABMIN(ipq,irs) & 
                               + 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))
                !ABMIN(irs,ipq) = ABMIN(ipq,irs) 

              endif  
           enddo
        enddo

      endif 
   enddo
enddo 

close(iunit1)
close(iunit2)

deallocate(ABKer,batch,WtKer,ints2,ints1,work2,work1)

end subroutine ModABMin_mithap

subroutine ACEneERPA_FFFF(ECorr,EVec,EVal,Occ,IGem, &
                          IndN,IndX,NOccup,NDimX,NBasis,IntFile)
implicit none

integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: IGem(NBasis),IndN(2,NDimX),IndX(NDimX)
integer,intent(in) :: NOccup
character(*),intent(in) :: IntFile
double precision,intent(out) :: ECorr
double precision,intent(in) :: EVec(NDimX,NDimX),EVal(NDimX)
double precision :: Occ(NBasis)

integer :: i,j,k,l,kl,kk,ip,iq,ir,is,ipq,irs
integer :: iunit,ISkippedEig
integer :: pos(NBasis,NBasis)
logical :: AuxCoeff(3,3,3,3)
double precision :: CICoef(NBasis),Cpq,Crs,SumY,Aux
double precision,allocatable :: work(:),ints(:,:),Skipped(:)
double precision,parameter :: SmallE = 1.d-3,BigE = 1.d8

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxCoeff = .true.
do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = .false.
            endif
         enddo
      enddo
   enddo
enddo

allocate(work(NBasis**2),ints(NBasis,NBasis),Skipped(NDimX))

ISkippedEig = 0
ECorr = 0

! FULL INTS
open(newunit=iunit,file=trim(IntFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work,ints,NBasis)

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)

                !if(.not.(IGem(ir).eq.IGem(is).and.IGem(ip).eq.IGem(iq)&
                !.and.IGem(ir).eq.IGem(ip)).and.ir.gt.is.and.ip.gt.iq) then
                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))) then

                   ISkippedEig = 0
                   SumY = 0 
                   do kk=1,NDimX
                      if(EVal(kk).gt.SmallE.and.EVal(kk).lt.BigE) then
                         SumY = SumY + EVec(ipq,kk)*EVec(irs,kk)
                      else
                         ISkippedEig = ISkippedEig + 1
                         Skipped(ISkippedEig) = EVal(kk)
                      endif
                   enddo
          
                   Aux = 2*Crs*Cpq*SumY
          
                   if(iq.Eq.is.and.ip.Eq.ir) then
                      Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
                   endif  

                   ECorr = ECorr + Aux*ints(j,i)
          
                ! endinf of If(IP.Gt.IR.And.IQ.Gt.IS) 
                endif

              endif
           enddo
        enddo

      endif

   enddo
enddo

close(iunit)

if(ISkippedEig/=0) then
  write(LOUT,'(/,1x,"The number of discarded eigenvalues is",i4)') &
       ISkippedEig
  do i=1,ISkippedEig
     write(LOUT,'(1x,a,i4,f15.8)') 'Skipped',i,Skipped(i)
  enddo
endif

deallocate(Skipped)
deallocate(ints,work)

end subroutine ACEneERPA_FFFF

subroutine make_tij_Y_unc(tmp2,tmp1,AOcc,BOcc,SBlockA,SBlockAIV,SBlockB,SBlockBIV,&
                          AIndN,BIndN,posA,posB,Sab,Sba,nblkA,nblkB,ANDimX,BNDimX,NBas)
implicit none

type(EBlockData)               :: SBlockA(nblkA),SBlockAIV,&
                                  SBlockB(nblkB),SBlockBIV
integer,intent(in)             :: nblkA,nblkB,ANDimX,BNDimX,NBas
integer,intent(in)             :: posA(NBas,NBas),posB(NBas,NBas),&
                                  AIndN(2,ANDimX),BIndN(2,BNDimX)
double precision,intent(in)    :: Sab(NBas,NBas),Sba(NBas,NBas),&
                                  AOcc(NBas),BOcc(NBas)
double precision,intent(inout) :: tmp1(ANDimX,BNDimX),&
                                  tmp2(ANDimX,BNDimX)

integer :: i,j,ir,is,irs,ip,iq,ipq
double precision :: fact

tmp1=0
do j=1,BNDimX

   ir = BIndN(1,j)
   is = BIndN(2,j)
   irs = posB(ir,is)

   do i=1,ANDimX

      ip = AIndN(1,i)
      iq = AIndN(2,i)
      ipq = posA(ip,iq)

      fact = (BOcc(ir)-BOcc(is)) * &
             (AOcc(ip)-AOcc(iq))

      tmp1(ipq,irs) = tmp1(ipq,irs) + Sab(iq,ir)*Sba(is,ip)*fact

   enddo
enddo
tmp2 = 0
!X_A.I.X_B
call abpm_tran_gen(tmp1,tmp2,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                   nblkA,nblkB,ANDimX,BNDimX,'XX')
!Y_A.I.Y_B
call abpm_tran_gen(tmp1,tmp2,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                   nblkA,nblkB,ANDimX,BNDimX,'YY')

tmp1 = 0
do j=1,BNDimX

   ir = BIndN(1,j)
   is = BIndN(2,j)
   irs = posB(ir,is)

   do i=1,ANDimX

      ip = AIndN(1,i)
      iq = AIndN(2,i)
      ipq = posA(ip,iq)

      fact = (BOcc(ir)-BOcc(is)) * &
             (AOcc(ip)-AOcc(iq))

      tmp1(ipq,irs) = tmp1(ipq,irs) - Sab(ip,ir)*Sba(is,iq)*fact

   enddo
enddo
!X_A.I.Y_B
call abpm_tran_gen(tmp1,tmp2,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                   nblkA,nblkB,ANDimX,BNDimX,'XY')
!Y_A.I.X_B
call abpm_tran_gen(tmp1,tmp2,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                   nblkA,nblkB,ANDimX,BNDimX,'YX')

end subroutine make_tij_Y_unc

subroutine make_uX_unc(Mon,WPot,pos,uX,NDimX,NBas)
implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: pos(NBas,NBas)
integer,intent(in) :: NBas,NDimX
double precision,intent(in)  :: WPot(NBas,NBas)
double precision,intent(out) :: uX(NDimX)

type(Y01BlockData),allocatable :: Y01BlockM(:)
integer :: i,j,ip,iq,ipq
double precision :: fact
double precision,allocatable :: OmM0(:),tind(:),test(:)
character(:),allocatable :: filexy0 

 if(Mon%Monomer==1) then
   filexy0 = 'XY0_A'
 elseif(Mon%Monomer==2) then
   filexy0 = 'XY0_B'
 endif

 allocate(Y01BlockM(Mon%NDimX),OmM0(Mon%NDimX))

 call convert_XY0_to_Y01(Mon,Y01BlockM,OmM0,NBas,filexy0)

 uX = 0
 do i=1,Mon%NDimX
    ip = Mon%IndN(1,i)
    iq = Mon%IndN(2,i)
    ipq = pos(ip,iq)

    fact = (Mon%Occ(ip)-Mon%Occ(iq)) * WPot(ip,iq)

    associate(Y => Y01BlockM(ipq))
       uX(Y%l1:Y%l2) = uX(Y%l1:Y%l2) + fact*Y%vec0(1:Y%n)
    end associate

 enddo
 ! here I guess mat-vec needed?
 !call dgemv('T',Mon%NDimX,Mon%NDimX,1d0,test,Mon%NDimX,tindA,1,0d0,uX,1)

 do i=1,Mon%NDimX
    associate(Y => Y01BlockM(i))
      deallocate(Y%vec0)
    end associate
 enddo
 deallocate(OmM0,Y01BlockM)

end subroutine make_uX_unc

subroutine mnoz_jeden(NBasis,ANDimX,BNDimX,EigVecA,IntAB,Amat,A,B)
implicit none

type(SystemBlock) :: A, B
integer,intent(in) :: NBasis,ANDimX,BNDimX
double precision,intent(in) :: EigVecA(ANDimX**2),IntAB(ANDimX,BNDimX)

integer :: i,j,ip,iq,ir,is,pq,rs
double precision :: fact
double precision,intent(inout) :: Amat(ANDimX,BNDimX)
double precision,allocatable :: tmp(:,:)

!allocate(tmp(ANDimX,BNDimX))

! (pq,i).(pq,rs)
do pq=1,A%NDimX
   ip = A%IndN(1,pq)
   iq = A%IndN(2,pq)
   do rs=1,B%NDimX
      ir = B%IndN(1,rs)
      is = B%IndN(2,rs)

      fact=1.0d0
      if(A%IGem(ip)==2.and.A%IGem(iq)==2.and. &
         B%IGem(ir)==2.and.B%IGem(is)==2) fact = 0d0

      do i=1,A%NDimX
         Amat(i,rs) = Amat(i,rs) + fact*EigVecA(pq+(i-1)*A%NDimX)*IntAB(pq,rs)
      enddo

   enddo
enddo

!deallocate(tmp)

end subroutine mnoz_jeden

end module sapt_exch
