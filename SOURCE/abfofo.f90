module abfofo
use types,only : LOUT,EblockData
use tran
implicit none


contains 

subroutine AB_CAS_FOFO(ABPLUS,ABMIN,ETot,URe,Occ,XOne, &
     IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     IntJFile,IntKFile,ACAlpha,AB1)
!
! COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
! 
! RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
! THE FOLLOWING SYMMETRY IS ASSUMED
! RDM2(ij,kl) = RDM2(kl,ij)
! ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
! COULOMB INTEGRALS ARE READ FROM IntJFile IN (FF|OO) FORMAT
! EXCHANGE INTEGRALS ARE READ FROM IntKFile IN (FO|FO) FORMAT
 implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1
character(*) :: IntJFile,IntKFile
double precision,intent(out) :: ABPLUS(NDim,NDim),ABMIN(NDim,NDim)
double precision,intent(out) :: ETot
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
double precision,intent(in)  :: ACAlpha
logical,intent(in) :: AB1

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit1,iunit2,ios
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,HNOCoef,val
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2

ABPLUS = 0
ABMIN  = 0
ETot   = 0

! set dimensions
NOccup = NAct + INActive
Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit1,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit1,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit1)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)


call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
ETot = ETot + 2*val

if(AB1) then
   do j=1,NBasis
      do i=1,NBasis
         if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
         if(IGem(i)==IGem(j)) HNO(i,j) = 0
      enddo
   enddo
else
   do j=1,NBasis
      do i=1,NBasis
         if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
      enddo
   enddo
endif

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

if(AB1) then
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               if((i==j).and.(j==k).and.(k==l)) then
                  AuxCoeff(i,j,k,l) = 0
               else
                  AuxCoeff(i,j,k,l) = ACAlpha
               endif
            enddo
         enddo
      enddo
   enddo
else
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               if((i==j).and.(j==k).and.(k==l)) then
                  AuxCoeff(i,j,k,l) = 1
               else
                  AuxCoeff(i,j,k,l) = ACAlpha
               endif
            enddo
         enddo
      enddo
   enddo
endif

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxI  = 0
AuxIO = 0
WMAT  = 0

if(AB1) then
   HNOCoef = -1
else
   HNOCoef = 1 - ACAlpha
endif

! does not work in Y01CAS?
!call AB_CAS_assemb(ABPLUS,ABMIN,Occ,RDM2val,HNO,pos, &
!           HNOCoef,AuxInd,AuxCoeff, &
!           IGem,IndN,IndX,NDimX,NDim,INActive,NOccup,NBasis, &
!           IntJFile,IntKFile)

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do ll=1,NOccup
   do kk=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      !call triang_to_sq2(work1,ints,NBasis)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      k = kk
      l = ll

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = HNOCoef*Occ(l)
         if(IGem(l)==1) then
            if(IGem(k)==2) then
               do i=INActive+1,NOccup
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
         
         if(IGem(l)==2) then
            if(IGem(k)==1) then
               do i=1,INActive
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
      endif
         
      ! AUXILIARY MATRIX AuxI AND AuxIO
      
      ! exchange 
      if(l<=INActive) then
         do i=1,NBasis
            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif
      if(l<=NOccup) then
         do i=1,NBasis
            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif

      ! AUXILIARY MATRIX WMAT
      if(l<=NOccup) then
         do ir=1,NOccup
            val = 0
            val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                 sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                 sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                 sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                 sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
            WMAT(k,ir) = WMAT(k,ir) + val
         enddo
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      ! Coulomb
      if(k>INActive.and.l<=NOccup) then
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then
            
            do iq=1,NOccup
               do ip=INActive+1,NBasis
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                     
                     val = 0
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                  
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     val = 0
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            enddo

         endif
      endif

      ! exchange
      if(l<=NOccup) then
         if(k>INActive) then
            ip = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then                  
                  do iq=1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = - AuxVal*val*ints(ir,iq)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               endif
            enddo
         endif
      endif

      ! T1 p->q
      if((k>INActive).and.(l<=NOccup)) then
         ip = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T3
      if((l<=NOccup).and.(k>INActive)) then
         iq = l
         ir = k
         
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T4+T6: AuxInd = 2
      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
         ir = k
         iu = l

         ! T4
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
         ! T6
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then                        
                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then  
                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif
      
   enddo
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do ll=1,NOccup
   do kk=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      
      k = kk
      l = ll

      if(k>NOccup.or.l>NOccup) cycle

      ! COMPUTE THE ENERGY FOR CHECKING
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*HNOCoef*Occ(k)
         if(k>INActive) then 
            do j=1,INActive
               do i=1,INActive
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         if(k<=INActive) then
            do j=INActive+1,NOccup
               do i=INActive+1,NOccup
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         do j=NOccup+1,NBasis
            do i=NOccup+1,NBasis
               HNO(i,j) = HNO(i,j) + val*ints(i,j)
            enddo
         enddo
         
      endif
          
      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! Coulomb
      if(k==l) then
         val = 2*Occ(k)
         if(k<=INActive) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
               enddo
            enddo
         endif
         if(k>INActive.and.k<=NOccup) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
               enddo
            enddo
         endif
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      
      ! exchange
      if(l<=NOccup) then
         if(k<=NOccup) then
            iq = k
            is = l
            
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = - AuxVal*val*ints(ip,ir)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               endif
            enddo
               
         endif
      endif

      ! T1 
      if((k<=NOccup).and.(l<=NOccup)) then
         iq = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo
            
      endif
         
      ! T3
      if(IGem(k)==2.and.IGem(l)==2) then
         
         do is=INActive+1,NOccup
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           
                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                RDM2val(iq,is,k,l)*ints(ip,ir)
                           
                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif

      ! T2+T5: AuxInd = 2
      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
         is = k
         iu = l
         
         ! T2
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo

         ! T5
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif
      
   enddo
enddo

close(iunit2)

write(LOUT,'(1x,a,5x,f15.8)') "CASSCF Energy (w/o ENuc)", ETot

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = IndX(ICol)
   if(ir/=is) then
      
      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = IndX(IRow)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif
         
         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val

      enddo

   endif
enddo

if(AB1) then
! symmetrize AB1
! hererxxx
   call sq_symmetrize(ABPLUS,NDimX)
   call sq_symmetrize(ABMIN,NDimX)
!  do ip=1,ndimx
!  do is=1,ip-1
!  ABPLUS(ip,is)=ABPLUS(is,ip)
!  ABMIN(ip,is)=ABMIN(is,ip)
!  enddo
!  enddo
endif

!print*, "AB-my",norm2(ABPLUS),norm2(ABMIN)

!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!HNO=transpose(HNO)
!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!!$
!!!$write(LOUT,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)') ETot
!!!$
!call sq_to_triang2(AuxI,work1,NBasis)
!write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!AuxI=transpose(AuxI)
!call sq_to_triang2(AuxI,work1,NBasis)
!write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
!call sq_to_triang2(AuxIO,work1,NBasis)
!write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!AuxIO=transpose(AuxIO)
!call sq_to_triang2(AuxIO,work1,NBasis)
!write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
!write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)

deallocate(RDM2val)
deallocate(ints,work2,work1)

end subroutine AB_CAS_FOFO  

subroutine ACABMAT0_FOFO(AMAT,BMAT,URe,Occ,XOne, &
                         IndN,IndX,IGem,C, &
                         NAct,INActive,NBasis,NDim,NDimX,NInte1,NGem, &
                         IntFileName,IntJFile,IntKFile,ISAPT,ACAlpha,IFlag)
!     IFlag = 1 - AMAT AND BMAT WILL CONTAIN (A+B)/C+/C+ AND (A-B)/C-/C-, RESPECTIVELY
!             0 - AMAT AND BMAT WILL CONTAIN A ANB B MATRICES, RESPECTIVELY
!
!     ACAlpha - Alpha-connection parameter in AC
!     HNO AND TwoMO are modified to correspond to an alpha-Hamiltonian
!
!     STRAIGHTFORWARD IMPLEMENTATION OF THE AMAT AND BMAT DEFINITIONS (NO SPECIAL CASES CONSIDERED) 
!     
!     COMPUTE THE A+B AND A-B MATRICES IN ERPA WITH APSG APPROXIMATION
!
!     NESTED COMMUTATOR 
!
implicit none

integer,intent(in) :: NAct,INActive,NBasis,NDim,NDimX,NInte1,ISAPT,NGem
character(*) :: IntJFile,IntKFile,IntFileName
double precision,intent(out) :: AMAT(NDimX,NDimX),BMAT(NDimX,NDimX)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1),C(NBasis)
double precision,intent(in)  :: ACAlpha
integer,intent(in) :: IndN(2,NDim),IndX(NDim)
integer,intent(in) :: IGem(NBasis),IFlag

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,iunit1,iunit2,ios
integer :: NOccup
integer :: IGemType
integer :: Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: HNO(NBasis,NBasis),HNOCoef
double precision :: AuxCoeff(NGem,NGem,NGem,NGem),AuxVal,val
double precision :: OccProd(NBasis,NBasis),CProd(NBasis,NBasis)
double precision :: AuxH(NBasis,NBasis,NGem),AuxXC(NBasis,NBasis,NGem)
double precision :: SaveA,SaveB
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
double precision,parameter :: Delta = 1.d-6

if(ISAPT==1) then
   write(6,'(1x,a)') 'Computing response-my'
else
   write(6,'(/,X,"***** COMPUTING AMAT, BMAT IN ACABMAT0 *****",/)')
endif

AMAT = 0
BMAT = 0

! set dimensions
NOccup = NAct + INActive

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

print*, 'HNO-test',norm2(HNO)

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
   enddo
enddo

do l=1,NGem
   do k=1,NGem
      do j=1,NGem
         do i=1,NGem
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = ACAlpha
            endif
         enddo
      enddo
   enddo
enddo

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = i
enddo

CProd = 0
OccProd = 0
do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)==IGem(j)) then
         CProd(i,j) = C(i)*C(j)         
      else
         OccProd(i,j) = Occ(i)*Occ(j)
      endif
   enddo
enddo

AuxH = 0
AuxXC = 0

HNOCoef = 1 - ACAlpha

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do l=1,NOccup
   do k=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      ints(:,NOccup+1:NBasis) = 0

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      if(IGem(k)/=IGem(l)) then
         val = HNOCoef*Occ(l) 
         do i=1,NBasis
            if(IGem(i)==IGem(k)) HNO(i,k) = HNO(i,k) - val*ints(i,l)
         enddo
      endif


      do IGemType=1,NGem-1

         if(IGem(l)==IGemType) then
            do i=1,NBasis
                AuxXC(i,k,IGemType) = AuxXC(i,k,IGemType) + AuxCoeff(IGem(i),IGem(l),IGem(l),IGem(k))*C(l)*ints(i,l)
            enddo
         else
            do i=1,NBasis
                AuxH(i,k,IGemType) = AuxH(i,k,IGemType) - AuxCoeff(IGem(i),IGem(l),IGem(l),IGem(k))*Occ(l)*ints(i,l)
            enddo
         endif

      enddo

      ! INTERGEMINAL PART

      ir = k
      is = l
      irs = pos(ir,is)
      if(irs>0) then
         do iq=1,NOccup
            do ip=1,NBasis
               ipq = pos(ip,iq)
               if(ipq>0) then

                  val = 2*AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                       *(-OccProd(ip,ir)+OccProd(ip,is)+OccProd(iq,ir)-OccProd(iq,is))*ints(ip,iq)   
                  BMAT(ipq,irs) = BMAT(ipq,irs) + val
                  AMAT(ipq,irs) = AMAT(ipq,irs) - val

               endif
            enddo    
         enddo
      endif

      ip = k
      is = l
      do iq=1,min(ip-1,NOccup)
         ipq = pos(ip,iq)
         if(ipq>0) then
            do ir=is+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
 
                  val = -AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                       *(-OccProd(ip,ir)+OccProd(ip,is)+OccProd(iq,ir)-OccProd(iq,is))*ints(ir,iq)
                  BMAT(ipq,irs) = BMAT(ipq,irs) + val
 
               endif
            enddo
         endif
      enddo

      ! INTRAGEMINAL PART

      ir = k
      iq = l
      do ip=iq+1,NBasis
         ipq = pos(ip,iq)
         if(ipq>0) then
            do is=1,min(ir-1,NOccup)
               irs = pos(ir,is)
               if(irs>0) then

                  val = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                       *(CProd(iq,ir)+CProd(ip,is))*ints(ip,is)
                  BMAT(ipq,irs) = BMAT(ipq,irs) + val
                  
               endif
            enddo
         endif
      enddo

      ip = k
      is = l
      do iq=1,min(ip-1,NOccup)
         ipq = pos(ip,iq)
         if(ipq>0) then
            do ir=is+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then

                  val = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                       *(CProd(ip,ir)+CProd(iq,is))*ints(ir,iq)
                  AMAT(ipq,irs) = AMAT(ipq,irs) + val
                  
               endif
            enddo
         endif
      enddo


   enddo
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do l=1,NOccup
   do k=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      if(k==l) then

         val = 2*HNOCoef*Occ(k)
         do j=1,NBasis
            if(IGem(j)/=IGem(k)) then
               do i=1,NBasis
                  if(IGem(i)==IGem(j)) HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            endif   
         enddo

         do IGemType=1,NGem
            if(IGem(k)/=IGemType) then
               val = 2*Occ(k)
               do j=1,NBasis
                  do i=1,NBasis
                     AuxH(i,j,IGemType) = AuxH(i,j,IGemType) & 
                                        + val*AuxCoeff(IGem(i),IGem(j),IGem(k),IGem(k))*ints(i,j)
                  enddo
               enddo
            endif
         enddo

      endif

      ! INTERGEMINAL PART

      iq = k
      is = l
      do ip=iq+1,NBasis
         ipq = pos(ip,iq)
         if(ipq>0) then
            do ir=is+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
 
                  val = -AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                       *(-OccProd(ip,ir)+OccProd(ip,is)+OccProd(iq,ir)-OccProd(iq,is))*ints(ir,ip)
                  AMAT(ipq,irs) = AMAT(ipq,irs) - val
 
               endif
            enddo
         endif
      enddo

      ! INTRAGEMINAL PART

      iq = k
      is = l
      do ir=is+1,NBasis
         irs = pos(ir,is)
         if(irs>0) then
            do ip=iq+1,NBasis
               ipq = pos(ip,iq)
               if(ipq>0) then

                  val = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                       *(CProd(iq,ir)+CProd(ip,is))*ints(ip,ir)
                  BMAT(ipq,irs) = BMAT(ipq,irs) + val
                  
               endif
            enddo
         endif
      enddo

      iq = k
      is = l
      do ir=is+1,NBasis
         irs = pos(ir,is)
         if(irs>0) then
            do ip=iq+1,NBasis
               ipq = pos(ip,iq)
               if(ipq>0) then

                  val = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                       *(CProd(ip,ir)+CProd(iq,is))*ints(ip,ir)
                  AMAT(ipq,irs) = AMAT(ipq,irs) + val
                  
               endif
            enddo
         endif
      enddo

   enddo
enddo

close(iunit2)

do ICol=1,NDimX
   ip = IndN(1,ICol)
   iq = IndN(2,ICol)
   ipq = ICol
 
   if(ip/=iq) then
      do IRow=1,NDimX
         ir = IndN(1,IRow)
         is = IndN(2,IRow)
         irs = IRow

            val = 0
            if(iq==ir) val = val &
                 + (Occ(ir)-Occ(is))*HNO(ip,is) &
                 + Occ(iq)*AuxH(ip,is,IGem(iq)) &
                 - Occ(is)*AuxH(ip,is,IGem(is)) &
                 - C(is)*AuxXC(ip,is,IGem(is))
            if(ip==is) val = val &
                 - (Occ(ir)-Occ(is))*HNO(iq,ir) &
                 + Occ(ip)*AuxH(iq,ir,IGem(ip)) &
                 - Occ(ir)*AuxH(iq,ir,IGem(ir)) &
                 - C(ir)*AuxXC(iq,ir,IGem(ir))
            if(ip==ir) val = val &
                 - C(ip)*AuxXC(iq,is,IGem(ip))
            if(iq==is) val = val &
                 - C(iq)*AuxXC(ip,ir,IGem(iq))
            BMAT(irs,ipq) = BMAT(irs,ipq) + val
         
            val = 0
            if(ip==ir) val = val &
                 + (Occ(ir)-Occ(is))*HNO(iq,is) &
                 + Occ(ip)*AuxH(iq,is,IGem(ip)) &
                 - Occ(is)*AuxH(iq,is,IGem(is)) &
                 - C(is)*AuxXC(iq,is,IGem(is))
            if(iq==is) val = val &
                 - (Occ(ir)-Occ(is))*HNO(ip,ir) &
                 + Occ(iq)*AuxH(ip,ir,IGem(iq)) &
                 - Occ(ir)*AuxH(ip,ir,IGem(ir)) & 
                 - C(ir)*AuxXC(ip,ir,IGem(ir))
            if(iq==ir) val = val &
                 - C(iq)*AuxXC(ip,is,IGem(iq))
            if(ip==is) val = val &
                 - C(ip)*AuxXC(iq,ir,IGem(ip))
            AMAT(irs,ipq) = AMAT(irs,ipq) + val

            if(IFlag/=0) then
!!$               if(abs(abs(C(ip))-abs(C(iq)))<=Delta*abs(C(iq)).or.& 
!!$                    abs(abs(C(ir))-abs(C(is)))<=Delta*abs(C(ir))) then
!!$               
!!$                  AMAT(irs,ipq) = 0
!!$                  BMAT(irs,ipq) = 0
!!$               else
                  
                  SaveA = AMAT(irs,ipq)
                  SaveB = BMAT(irs,ipq)
                  
                  val = (C(ip) + C(iq))*(C(ir) + C(is))
                  if(val==0) then
                     AMAT(irs,ipq) = 0
                  else
                     AMAT(irs,ipq) = (SaveA+SaveB)/val
                  endif
                  val = (C(ip) - C(iq))*(C(ir) - C(is))
                  if(val==0) then
                     BMAT(irs,ipq) = 0
                  else
                     BMAT(irs,ipq) = (SaveA-SaveB)/val
                  endif
               endif
!!$            endif
                        
         enddo
      endif
   enddo

print*, "AB-my",norm2(AMAT),norm2(BMAT)

deallocate(ints,work2,work1)

end subroutine ACABMAT0_FOFO

subroutine Y01CAS_FOFO(Occ,URe,XOne,ABPLUS,ABMIN, &
     propfile0,propfile1, &
     xy0file, &
     IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     NoSt,IntFileName,IntJFile,IntKFile,IFlag0,ETot,ECorr)
!
!     A ROUTINE FOR COMPUTING Y VECTORS AND EIGENVALUES OF ERPA 
!     IN THE 1ST-ORDER APPROXIMATION
!
!     IFlag0 = 1 - compute only 0th-order Y [EigY] and 0th-order omega [Eig] 
!              0 - compute both 0th-order and 1st-order Y [EigY1] and omega [Eig1]
use timing
!
implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1,NoSt
character(*) :: IntFileName,IntJFile,IntKFile
character(*) :: propfile0,propfile1,xy0file
double precision,intent(out) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
!double precision,intent(out) :: EigY(NDimX,NDimX),EigY1(NDimX,NDimX)
!double precision,intent(out) :: Eig(NDimX),Eig1(NDimX)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis),IFlag0
double precision,intent(out),optional :: ETot,ECorr

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: iunit1,iunit2
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,val
double precision :: EnDummy,Aux,Crs,Cpq,EIntra,EAll
double precision :: EMP2
double precision :: work(NBasis,NBasis)
double precision,allocatable :: EigY(:,:),EigY1(:,:)
double precision,allocatable :: Eig(:),Eig1(:)
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
!double precision,allocatable :: EigYt(:,:),EigXt(:,:),Eigt(:)
double precision,allocatable :: ABP(:,:),ABM(:,:),workA(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2
integer :: ii,jj,ipos,jpos,iblk,jblk,nblk
integer :: nAA,tmpAA(NAct*(NAct-1)/2),limAA(2)
integer :: nAI(INActive),tmpAI(NAct,1:INActive),limAI(2,1:INActive)
integer :: nAV(INActive+NAct+1:NBasis),tmpAV(NAct,INActive+NAct+1:NBasis),limAV(2,INActive+NAct+1:NBasis)
integer :: nIV,tmpIV(INActive*(NBasis-NAct-INActive)),limIV(2)
!
type(EblockData),allocatable :: Eblock(:) 
type(EblockData) :: EblockIV
!!
double precision,parameter :: Thresh = 1.D-12
double precision :: tmpEn
!! test timings
double precision :: Tcpu,Twall

! timing 
call clock('START',Tcpu,Twall)

ABPLUS = 0
ABMIN  = 0
if(present(ETot)) ETot = 0

! test EUGENE
 tmpEn = 0

! set dimensions
NOccup = NAct + INActive

Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
tmpEn = tmpEn + 2*val
print*, 'ONE ELECTRON ENERGY:', tmpEn
if(present(ETot)) ETot = ETot + 2*val

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = 0d0
   enddo
enddo

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = 0
            endif
         enddo
      enddo
   enddo
enddo

nAA = 0
nAI = 0
nAV = 0
nIV = 0
tmpAA = 0
tmpAI = 0
tmpAV = 0
tmpIV = 0
limAA = 0
limAI = 0
limAV = 0
limIV = 0

do i=1,NDimX
   ip = IndN(1,i)
   iq = IndN(2,i)
   if(IGem(ip)==2.and.IGem(iq)==2) then
      nAA = nAA + 1
      tmpAA(nAA) = i
    elseif(IGem(ip)==2.and.IGem(iq)==1) then
       nAI(iq) = nAI(iq) + 1
       tmpAI(nAI(iq),iq) = i
    elseif(IGem(ip)==3.and.IGem(iq)==2) then
       nAV(ip) = nAV(ip) + 1
       tmpAV(nAV(ip),ip) = i
    elseif(IGem(ip)==3.and.IGem(iq)==1) then
       nIV = nIV + 1
       tmpIV(nIV) = i
   endif
enddo

pos = 0
ipos = 0
limAA(1) = ipos + 1
do ii=1,nAA
   i = tmpAA(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limAA(2) = ipos
do is=1,INActive
   limAI(1,is) = ipos + 1
   do ii=1,nAI(is)
      i = tmpAI(ii,is)
      ! print*, ii,is,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAI(2,is) = ipos
enddo
!print*, ' '
do ir=NOccup+1,NBasis
   limAV(1,ir) = ipos + 1
   do ii=1,nAV(ir)
      i = tmpAV(ii,ir)
      ! print*, ii,ir,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAV(2,ir) = ipos
enddo
limIV(1) = ipos + 1
do ii=1,nIV
   i = tmpIV(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limIV(2) = ipos

!not working... ???
!call AB_CAS_assemb(ABPLUS,ABMIN,Occ,RDM2val,HNO,pos, &
!           1d0,AuxInd,AuxCoeff, &
!           IGem,IndN,IndX,NDimX,NDim,INActive,NOccup,NBasis, &
!           IntJFile,IntKFile)

AuxI  = 0
AuxIO = 0
WMAT  = 0

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do ll=1,NOccup
   do kk=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      k = kk
      l = ll

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = Occ(l)
         if(IGem(l)==1) then
            if(IGem(k)==2) then
               do i=INActive+1,NOccup
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
         
         if(IGem(l)==2) then
            if(IGem(k)==1) then
               do i=1,INActive
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO
      
      ! exchange 
      if(l<=INActive) then
         do i=1,NBasis
            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif
      if(l<=NOccup) then
         do i=1,NBasis
            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif

! AUXILIARY MATRIX WMAT
         if(l<=NOccup) then
            do ir=1,NOccup
               val = 0
               val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                    sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                    sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                    sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                    sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
               WMAT(k,ir) = WMAT(k,ir) + val
            enddo
         endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      ! Coulomb
      if(k>INActive.and.l<=NOccup) then
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then
            
            do iq=1,NOccup
               do ip=INActive+1,NBasis
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                     
                     val = 0
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                  
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     val = 0
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            enddo

         endif
      endif

      ! exchange
      if(l<=NOccup) then
         if(k>INActive) then
            ip = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then                  
                  do iq=1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = - AuxVal*val*ints(ir,iq)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               endif
            enddo
         endif
      endif

      ! T1 p->q
      if((k>INActive).and.(l<=NOccup)) then
         ip = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T3
      if((l<=NOccup).and.(k>INActive)) then
         iq = l
         ir = k
         
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T4+T6: AuxInd = 2
      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
         ir = k
         iu = l

         ! T4
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
         ! T6
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then                        
                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then  
                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif

   enddo
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do ll=1,NOccup
   do kk=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      
      k = kk
      l = ll

      if(k>NOccup.or.l>NOccup) cycle

      if(present(ETot)) then  
      ! COMPUTE THE ENERGY FOR CHECKING
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))
      endif

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*Occ(k)
         if(k>INActive) then 
            do j=1,INActive
               do i=1,INActive
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         if(k<=INActive) then
            do j=INActive+1,NOccup
               do i=INActive+1,NOccup
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         do j=NOccup+1,NBasis
            do i=NOccup+1,NBasis
               HNO(i,j) = HNO(i,j) + val*ints(i,j)
            enddo
         enddo
         
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! Coulomb
      if(k==l) then
         val = 2*Occ(k)
         if(k<=INActive) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
               enddo
            enddo
         endif
         if(k>INActive.and.k<=NOccup) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
               enddo
            enddo
         endif
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      
      ! exchange
      if(l<=NOccup) then
         if(k<=NOccup) then
            iq = k
            is = l
            
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = - AuxVal*val*ints(ip,ir)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               endif
            enddo
               
         endif
      endif

      ! T1 
      if((k<=NOccup).and.(l<=NOccup)) then
         iq = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo
            
      endif

      ! T3
      if(IGem(k)==2.and.IGem(l)==2) then
         
         do is=INActive+1,NOccup
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           
                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                RDM2val(iq,is,k,l)*ints(ip,ir)
                           
                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif

      ! T2+T5: AuxInd = 2
      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
         is = k
         iu = l
         
         ! T2
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo

         ! T5
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif

   enddo
enddo

close(iunit2)
deallocate(work1)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = pos(ir,is)
      
      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = pos(ip,iq)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val
         
      enddo
enddo

call clock('AB0MAT',Tcpu,Twall)

allocate(EBlock(1+NBasis-NAct))

nblk = 0

! put pack into separate subroutine
!pack AA
if(nAA>0) then

   nblk = nblk + 1
   associate(B => Eblock(nblk))
     B%n  = nAA
     B%l1 = limAA(1)
     B%l2 = limAA(2)
     allocate(B%pos(B%n))
     B%pos(1:B%n) = tmpAA(1:B%n)

     allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
     allocate(ABP(B%n,B%n),ABM(B%n,B%n))

     ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
     ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
     if(NoSt==1) then
        call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
     elseif(NoSt>1) then
        call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
     endif
   
     deallocate(ABM,ABP)
   
   end associate

endif

write(*,'(/,1x,a)') 'Act-Act   block diagonalized!'

!pack AI
do iq=1,INActive
   if(nAI(iq)>0) then

      nblk = nblk + 1
      associate(B => Eblock(nblk))
        B%n  = nAI(iq)
        B%l1 = limAI(1,iq)
        B%l2 = limAI(2,iq)
        allocate(B%pos(B%n))
        B%pos(1:B%n) = tmpAI(1:B%n,iq)

        allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
        allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  
        ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
        ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
        if(NoSt==1) then
           call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        elseif(NoSt>1) then
           call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        endif

        deallocate(ABM,ABP)

      end associate

   endif
enddo

print*, 'Act-InAct block diagonalized!'

!pack AV
do ip=NOccup+1,NBasis
   if(nAV(ip)>0) then

      nblk = nblk + 1
      associate(B => Eblock(nblk))
        B%n  = nAV(ip)
        B%l1 = limAV(1,ip)
        B%l2 = limAV(2,ip)
        allocate(B%pos(B%n))
        B%pos(1:B%n) = tmpAV(1:B%n,ip)

        allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
        allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  
        ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
        ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
        if(NoSt==1) then
           call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        elseif(NoSt>1) then
           call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        endif

        deallocate(ABM,ABP)

      end associate

    endif
enddo

print*, 'Act-Virt  block diagonalized!'

!pack IV
associate(B => EblockIV)

  B%l1 = limIV(1)
  B%l2 = limIV(2)
  B%n  = B%l2-B%l1+1
  allocate(B%pos(B%n))
  B%pos(1:B%n) = tmpIV(1:B%n)

  allocate(B%vec(B%n))
  do i=1,B%n
     ii = B%l1+i-1
     B%vec(i) = ABPLUS(ii,ii)
  enddo

end associate

 call clock('PACK',Tcpu,Twall)

!! dump EBLOCKS: X(0),Y(0)
call dump_Eblock(Eblock,EblockIV,Occ,IndN,nblk,NBasis,NDimX,xy0file)

if(IFlag0==0) then
   
   ! AB(1) PART
   call AB_CAS_FOFO(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,&
                 IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
                 NInte1,IntJFile,IntKFile,1d0,.true.)

  call clock('AB(1)',Tcpu,Twall)
!  ! B = A.X
!  ! C = X^T.B
   allocate(workA(NDimX,NDimX))

   call ABPM_TRAN(ABPLUS,workA,EBlock,EBlockIV,nblk,NDimX,.true.)
   call clock('MULT-1',Tcpu,Twall)
   ABPLUS=workA
   call ABPM_TRAN(ABMIN,workA,EBlock,EBlockIV,nblk,NDimX,.false.)
   ABMIN=workA

   call clock('MULT-2',Tcpu,Twall)

   deallocate(workA)

   allocate(Eig(NDimX),Eig1(NDimX))
 
   ! unpack Eig (1)
   do iblk=1,nblk
      associate(B => Eblock(iblk))
 
   !     do i=1,B%n
   !        ipos = B%pos(i)
   !        EigY(ipos,B%l1:B%l2) = B%matY(i,1:B%n)
   !        if(IFlag0==0) EigY1(ipos,B%l1:B%l2) = B%matX(i,1:B%n)
   !     enddo
        Eig(B%l1:B%l2) = B%vec(1:B%n)

   !     deallocate(B%matY,B%matX,B%vec)
   !     deallocate(B%pos)

      end associate
   enddo

   !unpack Eig (2, IV part)
   associate(B => EblockIV)
 
     do i=1,B%n
        ii = B%l1+i-1
        ipos = B%pos(i)
        !EigY(ipos,ii) = 1d0/sqrt(2d0)
        !if(IFlag0==0) EigY1(ipos,ii) = 1d0/sqrt(2d0)
        Eig(ii) = B%vec(i)
     enddo

   !  deallocate(B%vec)
   !  deallocate(B%pos)

   end associate
   call clock('UNPACK',Tcpu,Twall)

   ! old ways...
   ! allocate(work1(NDimX**2))

   ! work1=ABPLUS.EigX
   ! ABPLUS=work1.EigX
   !call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,EigY1,NDimX,0d0,work1,NDimX)
   !call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY1,NDimX,0d0,ABPLUS,NDimX)

   ! work1=ABMIN.EigY
   ! ABMIN=work1.EigY
   !call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABMIN,NDimX,EigY,NDimX,0d0,work1,NDimX) 
   !call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY,NDimX,0d0,ABMIN,NDimX)

   ! deallocate(work1)

endif

deallocate(RDM2val,work2)
deallocate(ints)

If(IFlag0==1) return

if(.not.present(ECorr)) then
! this part for E2disp
! ------------------------------------------------------------------------------

allocate(EigY(NDimX,NDimX),EigY1(NDimX,NDimX))

EigY1 = 0
EigY  = 0

! unpack (1)
do iblk=1,nblk
   associate(B => Eblock(iblk))

     do i=1,B%n
        ipos = B%pos(i)
        EigY(ipos,B%l1:B%l2) = B%matY(i,1:B%n)
        if(IFlag0==0) EigY1(ipos,B%l1:B%l2) = B%matX(i,1:B%n)
     enddo

     deallocate(B%matY,B%matX,B%vec)
     deallocate(B%pos)

   end associate
enddo

!unpack (2, IV part)
associate(B => EblockIV)

  do i=1,B%n
     ii = B%l1+i-1
     ipos = B%pos(i)
     EigY(ipos,ii) = 1d0/sqrt(2d0)
     if(IFlag0==0) EigY1(ipos,ii) = 1d0/sqrt(2d0)
  enddo

  deallocate(B%vec)
  deallocate(B%pos)

end associate


do i=1,NDimX
   Eig1(i)=ABPLUS(i,i)+ABMIN(i,i)
enddo

EigY1 = 0
!do j=1,NDimX
!   if(Eig(j)/=0d0) then
!      do i=1,NDimX
!         if(Eig(i)/=0d0) then
!            val = (ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
!            if(Abs(Eig(i)-Eig(j))>Thresh) then
!               val = val + (ABPLUS(i,j)+ABMIN(i,j))/(Eig(j)-Eig(i))
!            endif
!            do ii=1,NDimX
!               EigY1(ii,j) = EigY1(ii,j) + val*EigY(ii,i)
!            enddo
!         endif
!      enddo
!   endif
!enddo

allocate(workA(NDimX,NDimX))
do j=1,NDimX
   if(Eig(j)/=0d0) then
      do i=1,NDimX
         if(Eig(i)/=0d0) then
            workA(i,j) = (ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
            if(Abs(Eig(i)-Eig(j))>Thresh) then
               workA(i,j) = workA(i,j) + (ABPLUS(i,j)+ABMIN(i,j))/(Eig(j)-Eig(i))
            endif
         endif
      enddo
   endif
enddo

call dgemm('N','N',NDimX,NDimX,NDimX,1d0,EigY,NDimX,workA,NDimX,0d0,EigY1,NDimX)

 ! dump response to a file!
 open(newunit=iunit,file=propfile0,form='unformatted')
 write(iunit) EigY
 write(iunit) Eig
 close(iunit)
 if(IFlag0==0) then
    open(newunit=iunit,file=propfile1,form='unformatted')
    write(iunit) EigY1
    write(iunit) Eig1
    close(iunit)
 endif

deallocate(workA)
deallocate(EigY1,EigY)

! ------------------------------------------------------------------------------

elseif(present(ECorr)) then
! for AC0Corr
! ------------------------------------------------------------------------------

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo
! original
!EigY1 = 0
!do j=1,NDimX
!   if(Eig(j)/=0d0) then
!      do i=1,NDimX
!         if(Eig(i)/=0d0) then
!            val = 2d0*(ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
!            do ii=1,NDimX
!               EigY1(ii,j) = EigY1(ii,j) + val*EigY(ii,i)
!            enddo
!         endif
!      enddo
!   endif
!enddo

!prepare new AMAT
allocate(workA(NDimX,NDimX))
workA=0
do j=1,NDimX
   if(Eig(j)/=0d0) then
      do i=1,NDimX
         if(Eig(i)/=0d0) workA(i,j) = 2d0*(ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
      enddo
   endif
enddo
call clock('APL-AMIN/EPL',Tcpu,Twall)

! old ways ...
!call dgemm('N','N',NDimX,NDimX,NDimX,1d0,EigY,NDimX,workA,NDimX,0d0,EigY1,NDimX)
!call dgemm('N','T',NDimX,NDimX,NDimX,1d0,EigY,NDimX,EigY1,NDimX,0d0,ABPLUS,NDimX)

! new
call ABPM_BACKTRAN(workA,ABPLUS,EBlock,EBlockIV,nblk,NDimX)

deallocate(workA)

! deallocate blocks
do iblk=1,nblk
   associate(B => Eblock(iblk))

     deallocate(B%matY,B%matX,B%vec)
     deallocate(B%pos)

   end associate
enddo
! (IV part)
associate(B => EblockIV)

  deallocate(B%vec)
  deallocate(B%pos)

end associate

call clock('ABPl-dgemm',Tcpu,Twall)
!print*, 'ABPLUS-MY',norm2(ABPLUS)

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo
!
! energy loop
allocate(work1(NBasis*NBasis),ints(NBasis,NBasis))
EAll = 0
EIntra = 0
open(newunit=iunit,file='FOFO',status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

ints = 0
kl   = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work1(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints(i,j) = work1((j-1)*NBasis+i)
           enddo
        enddo
        ints(:,NOccup+1:NBasis) = 0

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i
                Crs = C(ir)+C(is)
                Cpq = C(ip)+C(iq)

                Aux = Crs*Cpq*ABPLUS(ipq,irs)
                EAll = EAll + Aux*ints(j,i)

                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))==1) EIntra = EIntra + Aux*ints(j,i)

              endif
           enddo
        enddo

      endif
   enddo
enddo

close(iunit)

call clock('ENE-loop',Tcpu,Twall)
!call check_mp2(EMP2,'FOFO',INActive,NOccup,NBasis)

ECorr = EAll-EIntra
!print*, 'EAll,EIntra',EAll,EIntra

deallocate(ints,work1)

endif


!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!HNO=transpose(HNO)
!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
!write(LOUT,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)') ETot
!
!call sq_to_triang2(AuxI,work1,NBasis)
!write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!AuxI=transpose(AuxI)
!call sq_to_triang2(AuxI,work1,NBasis)
!write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
!call sq_to_triang2(AuxIO,work1,NBasis)
!write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!AuxIO=transpose(AuxIO)
!call sq_to_triang2(AuxIO,work1,NBasis)
!write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
!write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)

deallocate(Eig1,Eig)

end subroutine Y01CAS_FOFO

subroutine ABPM_BACKTRAN(AMAT,AOUT,EBlock,EBlockIV,nblk,NDimX)
implicit none

integer,intent(in) :: nblk,NDimX
double precision,intent(in) :: AMAT(NDimX,NDimX)
double precision,intent(inout) :: AOUT(NDimX,NDimX)

type(EBlockData),intent(in) :: EBlock(nblk),EBlockIV

integer :: i,j,ii,jj,ipos,jpos,iblk,jblk
double precision,allocatable :: ABP(:,:),ABM(:,:)
double precision :: fac

!write(*,*) 'ABPM_BACKTRAN-GO, GO!'

fac = 1.d0/sqrt(2.d0)

AOUT=0

do jblk=1,nblk
   associate( jB => Eblock(jblk) )
   do iblk=1,nblk
      associate( iB => Eblock(iblk))

        allocate(ABP(iB%n,jB%n),ABM(iB%n,jB%n))

        ABP = AMAT(iB%l1:iB%l2,jB%l1:jB%l2)

        call dgemm('N','T',iB%n,jB%n,jB%n,1d0,ABP,iB%n,jB%matY,jB%n,0d0,ABM,iB%n)
        call dgemm('N','N',iB%n,jB%n,iB%n,1d0,iB%matY,iB%n,ABM,iB%n,0d0,ABP,iB%n)

        do j=1,jB%n
           jpos = jB%pos(j)
           do i=1,iB%n
              ipos = iB%pos(i)
              AOUT(ipos,jpos) = ABP(i,j)
           enddo
        enddo

        deallocate(ABM,ABP)

      end associate
   enddo
   end associate
enddo

associate(B => EblockIV)

  do iblk=1,nblk
     associate(iB => Eblock(iblk))

       allocate(ABP(iB%n,B%n),ABM(iB%n,B%n))

       ABP = AMAT(iB%l1:iB%l2,B%l1:B%l2)

       call dgemm('N','N',iB%n,B%n,iB%n,fac,iB%matY,iB%n,ABP,iB%n,0d0,ABM,iB%n)

       do j=1,B%n
          jpos = B%pos(j)
          do i=1,iB%n
             ipos = iB%pos(i)
             AOUT(ipos,jpos) = ABM(i,j)
          enddo
       enddo

       deallocate(ABM,ABP)

     end associate
  enddo

  do jblk=1,nblk
     associate(jB => Eblock(jblk))

       allocate(ABP(B%n,jB%n),ABM(B%n,jB%n))

       ABP = AMAT(B%l1:B%l2,jB%l1:jB%l2)

       call dgemm('N','T',B%n,jB%n,jB%n,fac,ABP,B%n,jB%matY,jB%n,0d0,ABM,B%n)

       do j=1,jB%n
          jpos = jB%pos(j)
          do i=1,B%n
             ipos = B%pos(i)
             AOUT(ipos,jpos) = ABM(i,j)
          enddo
       enddo

       deallocate(ABM,ABP)

     end associate
  enddo

  do j=1,B%n
     jj = B%l1+j-1
     jpos = B%pos(j)
     do i=1,B%n
        ii = B%l1+i-1
        ipos = B%pos(i)
        AOUT(ipos,jpos) = AMAT(ii,jj)*0.5d0
     enddo
  enddo

end associate

end subroutine ABPM_BACKTRAN

subroutine ABPM_TRAN(AMAT,AOUT,EBlock,EBlockIV,nblk,NDimX,isPl)
implicit none

integer,intent(in) :: nblk,NDimX
logical,intent(in) :: isPl
double precision,intent(in) :: AMAT(NDimX,NDimX)
double precision,intent(inout) :: AOUT(NDimX,NDimX)

type(EBlockData),intent(in) :: EBlock(nblk),EBlockIV

integer :: i,j,ii,jj,ipos,jpos,iblk,jblk
double precision,allocatable :: ABP(:,:),ABM(:,:)
double precision :: fac

fac = 1.d0/sqrt(2.d0)

AOUT=0

do jblk=1,nblk
   associate( jB => Eblock(jblk) )
   do iblk=1,nblk
      associate( iB => Eblock(iblk))

        allocate(ABP(iB%n,jB%n),ABM(iB%n,jB%n))
        do j=1,jB%n
           jpos = jB%pos(j)
           do i=1,iB%n
              ipos = iB%pos(i)
              ABP(i,j) = AMAT(ipos,jpos)
           enddo
        enddo
        if(isPl) then
           call dgemm('N','N',iB%n,jB%n,jB%n,1d0,ABP,iB%n,jB%matX,jB%n,0d0,ABM,iB%n)
           call dgemm('T','N',iB%n,jB%n,iB%n,1d0,iB%matX,iB%n,ABM,iB%n,0d0,ABP,iB%n)
        else
           call dgemm('N','N',iB%n,jB%n,jB%n,1d0,ABP,iB%n,jB%matY,jB%n,0d0,ABM,iB%n)
           call dgemm('T','N',iB%n,jB%n,iB%n,1d0,iB%matY,iB%n,ABM,iB%n,0d0,ABP,iB%n)
        endif

        AOUT(iB%l1:iB%l2,jB%l1:jB%l2) = ABP
        deallocate(ABM,ABP)

      end associate
   enddo
   end associate
enddo

associate(B => EblockIV)

  if(B%n>0) then
     do iblk=1,nblk
        associate(iB => Eblock(iblk))

          allocate(ABP(iB%n,B%n),ABM(iB%n,B%n))
          do j=1,B%n
             jpos = B%pos(j)
             do i=1,iB%n
                ipos = iB%pos(i)
                ABP(i,j) = AMAT(ipos,jpos)
             enddo
          enddo
          if(isPl) then
             call dgemm('T','N',iB%n,B%n,iB%n,fac,iB%matX,iB%n,ABP,iB%n,0d0,ABM,iB%n) 
          else
             call dgemm('T','N',iB%n,B%n,iB%n,fac,iB%matY,iB%n,ABP,iB%n,0d0,ABM,iB%n) 
          endif

          AOUT(iB%l1:iB%l2,B%l1:B%l2) = ABM
          deallocate(ABM,ABP)

        end associate
     enddo

     do jblk=1,nblk
        associate(jB => Eblock(jblk))

          allocate(ABP(B%n,jB%n),ABM(B%n,jB%n))
          do j=1,jB%n
             jpos = jB%pos(j)
             do i=1,B%n
                ipos = B%pos(i)
                ABP(i,j) = AMAT(ipos,jpos)
             enddo
          enddo
          if(isPl) then
             call dgemm('N','N',B%n,jB%n,jB%n,fac,ABP,B%n,jB%matX,jB%n,0d0,ABM,B%n)
          else
             call dgemm('N','N',B%n,jB%n,jB%n,fac,ABP,B%n,jB%matY,jB%n,0d0,ABM,B%n)
          endif

          AOUT(B%l1:B%l2,jB%l1:jB%l2) = ABM
          deallocate(ABM,ABP)

        end associate
     enddo

     do j=1,B%n
        jj = B%l1+j-1
        jpos = B%pos(j)
        do i=1,B%n
           ii = B%l1+i-1
           ipos = B%pos(i)
           AOUT(ii,jj) = AMAT(ipos,jpos)*0.5d0
        enddo
     enddo
  endif

end associate

end subroutine ABPM_TRAN

subroutine Y01CASLR_FOFO(Occ,URe,XOne,ABPLUS,ABMIN, &
     MultpC,NSymNO, &
     SRKer,Wt,OrbGrid, &
     propfile0,propfile1, & 
     IndN,IndX,IGemIN,NAct,INActive,NGrid,NDimX,NBasis,NDim,NInte1, &
     NoSt,IntFileName,IntJFile,IntKFile,IFlag0,IFunSRKer,ETot,ECorr)
!
!     A ROUTINE FOR COMPUTING Y VECTORS AND EIGENVALUES OF ERPA 
!     IN THE 1ST-ORDER APPROXIMATION
!
!     IFlag0 = 1 - compute only 0th-order Y [EigY] and 0th-order omega [Eig] 
!              0 - compute both 0th-order and 1st-order Y [EigY1] and omega [Eig1]
use timing

implicit none
! here!!!!
integer,intent(in) :: NAct,INActive,NGrid,NDimX,NBasis,NDim,NInte1,NoSt
character(*) :: IntFileName,IntJFile,IntKFile
character(*) :: propfile0,propfile1
double precision,intent(out) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
double precision,intent(in) :: SRKer(NGrid),Wt(NGrid),OrbGrid(NGrid,NBasis)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis),IFlag0,IFunSRKer, &
                      MultpC(15,15),NSymNO(NBasis)
double precision,intent(out),optional :: ETot,ECorr

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: iunit1,iunit2
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
logical :: EigChck(NDimX)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,val
double precision :: EnDummy,Aux,Crs,Cpq,EIntra,EAll
double precision,allocatable :: EigY(:,:),EigY1(:,:)
double precision,allocatable :: Eig(:),Eig1(:)
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
double precision,allocatable :: EigYt(:,:),EigXt(:,:),Eigt(:)
double precision,allocatable :: ABP(:,:),ABM(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2
integer :: ii,jj,ipos
integer :: nAA,tmpAA(NAct*(NAct-1)/2),limAA(2)
integer :: nAI(INActive),tmpAI(NAct,1:INActive),limAI(2,1:INActive)
integer :: nAV(INActive+NAct+1:NBasis),tmpAV(NAct,INActive+NAct+1:NBasis),limAV(2,INActive+NAct+1:NBasis)
integer :: nIV,tmpIV(INActive*(NBasis-NAct-INActive)),limIV(2)
double precision,parameter :: Thresh = 1.D-12
double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

ABPLUS = 0
ABMIN  = 0
ETot = 0

! set dimensions
NOccup = NAct + INActive
Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))
allocate(EigY(NDimX,NDimX),Eig(NDimX))
allocate(EigY1(NDimX,NDimX),Eig1(NDimX))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
ETot = ETot + 2*val

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = 0d0
   enddo
enddo

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = 0
            endif
         enddo
      enddo
   enddo
enddo

nAA = 0
nAI = 0
nAV = 0
nIV = 0
tmpAA = 0
tmpAI = 0
tmpAV = 0
tmpIV = 0
limAA = 0
limAI = 0
limAV = 0
limIV = 0

do i=1,NDimX
   ip = IndN(1,i)
   iq = IndN(2,i)
   if(IGem(ip)==2.and.IGem(iq)==2) then
      nAA = nAA + 1
      tmpAA(nAA) = i
    elseif(IGem(ip)==2.and.IGem(iq)==1) then
       nAI(iq) = nAI(iq) + 1
       tmpAI(nAI(iq),iq) = i
    elseif(IGem(ip)==3.and.IGem(iq)==2) then
       nAV(ip) = nAV(ip) + 1
       tmpAV(nAV(ip),ip) = i
    elseif(IGem(ip)==3.and.IGem(iq)==1) then
       nIV = nIV + 1
       tmpIV(nIV) = i
   endif
enddo

pos = 0
ipos = 0
limAA(1) = ipos + 1
do ii=1,nAA
   i = tmpAA(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limAA(2) = ipos
do is=1,INActive
   limAI(1,is) = ipos + 1
   do ii=1,nAI(is)
      i = tmpAI(ii,is)
      ! print*, ii,is,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAI(2,is) = ipos
enddo
!print*, ' '
do ir=NOccup+1,NBasis
   limAV(1,ir) = ipos + 1
   do ii=1,nAV(ir)
      i = tmpAV(ii,ir)
      ! print*, ii,ir,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAV(2,ir) = ipos
enddo
limIV(1) = ipos + 1
do ii=1,nIV
   i = tmpIV(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limIV(2) = ipos

!not working... ???
!call AB_CAS_assemb(ABPLUS,ABMIN,Occ,RDM2val,HNO,pos, &
!           1d0,AuxInd,AuxCoeff, &
!           IGem,IndN,IndX,NDimX,NDim,INActive,NOccup,NBasis, &
!           IntJFile,IntKFile)

AuxI  = 0
AuxIO = 0
WMAT  = 0

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do ll=1,NOccup
   do kk=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      !call triang_to_sq2(work1,ints,NBasis)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      k = kk
      l = ll

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = Occ(l)
         if(IGem(l)==1) then
            if(IGem(k)==2) then
               do i=INActive+1,NOccup
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
         
         if(IGem(l)==2) then
            if(IGem(k)==1) then
               do i=1,INActive
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO
      
      ! exchange 
      if(l<=INActive) then
         do i=1,NBasis
            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif
      if(l<=NOccup) then
         do i=1,NBasis
            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif

! AUXILIARY MATRIX WMAT
         if(l<=NOccup) then
            do ir=1,NOccup
               val = 0
               val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                    sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                    sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                    sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                    sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
               WMAT(k,ir) = WMAT(k,ir) + val
            enddo
         endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      ! Coulomb
      if(k>INActive.and.l<=NOccup) then
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then
            
            do iq=1,NOccup
               do ip=INActive+1,NBasis
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                     
                     val = 0
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                  
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     val = 0
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            enddo

         endif
      endif

      ! exchange
      if(l<=NOccup) then
         if(k>INActive) then
            ip = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then                  
                  do iq=1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = - AuxVal*val*ints(ir,iq)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               endif
            enddo
         endif
      endif

      ! T1 p->q
      if((k>INActive).and.(l<=NOccup)) then
         ip = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T3
      if((l<=NOccup).and.(k>INActive)) then
         iq = l
         ir = k
         
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T4+T6: AuxInd = 2
      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
         ir = k
         iu = l

         ! T4
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
         ! T6
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then                        
                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then  
                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif

   enddo
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do ll=1,NOccup
   do kk=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      
      k = kk
      l = ll

      if(k>NOccup.or.l>NOccup) cycle

      ! COMPUTE THE ENERGY FOR CHECKING
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*Occ(k)
         if(k>INActive) then 
            do j=1,INActive
               do i=1,INActive
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         if(k<=INActive) then
            do j=INActive+1,NOccup
               do i=INActive+1,NOccup
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         do j=NOccup+1,NBasis
            do i=NOccup+1,NBasis
               HNO(i,j) = HNO(i,j) + val*ints(i,j)
            enddo
         enddo
         
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! Coulomb
      if(k==l) then
         val = 2*Occ(k)
         if(k<=INActive) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
               enddo
            enddo
         endif
         if(k>INActive.and.k<=NOccup) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
               enddo
            enddo
         endif
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      
      ! exchange
      if(l<=NOccup) then
         if(k<=NOccup) then
            iq = k
            is = l
            
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = - AuxVal*val*ints(ip,ir)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               endif
            enddo
               
         endif
      endif

      ! T1 
      if((k<=NOccup).and.(l<=NOccup)) then
         iq = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo
            
      endif

      ! T3
      if(IGem(k)==2.and.IGem(l)==2) then
         
         do is=INActive+1,NOccup
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           
                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                RDM2val(iq,is,k,l)*ints(ip,ir)
                           
                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif

      ! T2+T5: AuxInd = 2
      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
         is = k
         iu = l
         
         ! T2
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo

         ! T5
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif

   enddo
enddo

close(iunit2)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = pos(ir,is)
      
      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = pos(ip,iq)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val
         
      enddo
enddo

call clock('Y01CASLR:ABMAT',Tcpu,Twall)
call clock('START',Tcpu,Twall)

! symmetrize AB(0)
call sq_symmetrize(ABPLUS,NDimX)
call sq_symmetrize(ABMIN,NDimX)

EigY1 = 0
EigY = 0

i = limAA(1)
j = limAA(2)

if(nAA>0) then

   allocate(ABP(nAA,nAA),ABM(nAA,nAA),EigYt(nAA,nAA),EigXt(nAA,nAA),Eigt(nAA))

   ABP = ABPLUS(i:j,i:j)
   ABM = ABMIN(i:j,i:j)
 
   if(IFunSRKer==1) then
      call ModABMin_Act_FOFO(Occ,SRKer,Wt,OrbGrid,ABM,&
             MultpC,NSymNO,tmpAA,&
             IndN,IndX,NDimX,NGrid,NBasis,&
             NOccup,NAct,INActive,nAA,'FOFO','FOFOERF')
   endif

   !print*, 'ACT-ACT-MY',nAA,norm2(ABP),norm2(ABM)
   if(NoSt==1) then
      call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAA)
   elseif(NoSt>1) then
      call ERPAVECYX(EigYt,EigXt,Eigt,ABP,ABM,nAA)
   endif
   
   do ii=1,nAA
      ipos = tmpAA(ii)
      EigY(ipos,i:j) = EigYt(ii,1:nAA)
      if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAA)
   enddo
   Eig(i:j) = Eigt(1:nAA)

   deallocate(Eigt,EigXt,EigYt,ABM,ABP)
endif

do iq=1,INActive
   i = limAI(1,iq)
   j = limAI(2,iq)

   if(nAI(iq)>0) then
      allocate(ABP(nAI(iq),nAI(iq)),ABM(nAI(iq),nAI(iq)),&
           EigYt(nAI(iq),nAI(iq)),EigXt(nAI(iq),nAI(iq)),Eigt(nAI(iq)))

      ABP = ABPLUS(i:j,i:j)
      ABM = ABMIN(i:j,i:j)
      !print*, 'AI-MY',iq,norm2(ABP),norm2(ABM)
      if(NoSt==1) then
         call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAI(iq))
      elseif(NoSt>1) then
         call ERPAVECYX(EigYt,EigXt,Eigt,ABP,ABM,nAI(iq))
      endif

      do ii=1,nAI(iq)
         ipos = tmpAI(ii,iq)
         EigY(ipos,i:j) = EigYt(ii,1:nAI(iq))
         if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAI(iq))
      enddo
      Eig(i:j) = Eigt(1:nAI(iq))
      
      deallocate(Eigt,EigXt,EigYt,ABM,ABP)
   endif 
enddo

do ip=NOccup+1,NBasis
   i = limAV(1,ip)
   j = limAV(2,ip)
   if(nAV(ip)>0) then
      allocate(ABP(nAV(ip),nAV(ip)),ABM(nAV(ip),nAV(ip)),&
           EigYt(nAV(ip),nAV(ip)),EigXt(nAV(ip),nAV(ip)),Eigt(nAV(ip)))

      ABP = ABPLUS(i:j,i:j)
      ABM = ABMIN(i:j,i:j)
      !print*, 'AV-MY',ip,norm2(ABP),norm2(ABM)
      if(NoSt==1) then
         call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAV(ip))
      elseif(NoSt>1) then
         call ERPAVECYX(EigYt,EigXt,Eigt,ABP,ABM,nAV(ip))
      endif

      do ii=1,nAV(ip)
         ipos = tmpAV(ii,ip)
         EigY(ipos,i:j) = EigYt(ii,1:nAV(ip))
         if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAV(ip))
      enddo
      Eig(i:j) = Eigt(1:nAV(ip))

      deallocate(Eigt,EigXt,EigYt,ABM,ABP)
   
   endif
enddo

do i=limIV(1),limIV(2)
   ii = tmpIV(i-limIV(1)+1)
   ip = IndN(1,ii)
   iq = IndN(2,ii)
   ! write(*,*) 'ABiv-my',ip,iq,ABPLUS(i,i)
   Eig(i) = ABPLUS(i,i)
enddo
do ii=1,nIV
   i = tmpiV(ii)
   EigY(i,limIV(1)+ii-1) = 1d0/sqrt(2d0)
   if(IFlag0==0) EigY1(i,limIV(1)+ii-1) = 1d0/sqrt(2d0)
enddo

deallocate(work1)

!print*, 'Eig-MY',norm2(Eig),norm2(EigY),norm2(EigY1)
call clock('Y01CAS:DIAG',Tcpu,Twall)
call clock('START',Tcpu,Twall)

if(IFlag0==1) return
!return

! AB(1) PART
call AB_CAS_FOFO(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,&
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,IntJFile,IntKFile,1d0,.true.)
!print*, 'AB1-MY',norm2(ABPLUS),norm2(ABMIN)

call clock('AB_CAS_FOFO',Tcpu,Twall)

if(IFunSRKer==1) then
   call ModABMin_FOFO(Occ,SRKer,Wt,OrbGrid,ABMIN,&
                 MultpC,NSymNO,&
                 IndN,IndX,NDimX,NGrid,NBasis,&
                 NAct,INActive,'FOFO','FOFOERF',.true.)
   !print*, 'ABM-MY',norm2(ABMIN)
endif
! here!!!! can this be made cheaper?
call clock('START',Tcpu,Twall)
allocate(work1(NDimX**2))
! work1=ABPLUS.EigX
! ABPLUS=work1.EigX
call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,EigY1,NDimX,0d0,work1,NDimX) 
call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY1,NDimX,0d0,ABPLUS,NDimX)

! work1=ABMIN.EigY
! ABMIN=work1.EigY
call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABMIN,NDimX,EigY,NDimX,0d0,work1,NDimX) 
call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY,NDimX,0d0,ABMIN,NDimX)

call clock('Y01CAS:dgemm',Tcpu,Twall)

if(.not.present(ECorr)) then
! this part for E2disp
! ------------------------------------------------------------------------------
do i=1,NDimX
   Eig1(i)=ABPLUS(i,i)+ABMIN(i,i)
enddo

EigY1 = 0
do j=1,NDimX
   if(Eig(j)/=0d0) then
      do i=1,NDimX
         if(Eig(i)/=0d0) then
            val = (ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
            if(Abs(Eig(i)-Eig(j))>Thresh) then
               val = val + (ABPLUS(i,j)+ABMIN(i,j))/(Eig(j)-Eig(i))
            endif
            do ii=1,NDimX
               EigY1(ii,j) = EigY1(ii,j) + val*EigY(ii,i)
            enddo
         endif
      enddo
   endif
enddo

 ! dump response to a file!
 open(newunit=iunit,file=propfile0,form='unformatted')
 write(iunit) EigY
 write(iunit) Eig
 close(iunit)
 if(IFlag0==0) then
 open(newunit=iunit,file=propfile1,form='unformatted')
 write(iunit) EigY1
 write(iunit) Eig1
 close(iunit)
 endif

! ------------------------------------------------------------------------------

elseif(present(ECorr)) then
! for AC0Corr
! ------------------------------------------------------------------------------

call clock('START',Tcpu,Twall)
do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do i=1,NDimX
    EigChck(i) = (Eig(i)/=0d0)
enddo

!print*, 'AB-MY',norm2(ABPLUS),norm2(ABMIN)
EigY1 = 0
do j=1,NDimX
!   if(Eig(j)/=0d0) then
   if(EigChck(j)) then
      do i=1,NDimX
!         if(Eig(i)/=0d0) then
          if(EigChck(i)) then
            val = 2d0*(ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
            call daxpy(NDimX,val,EigY(:,i),1,EigY1(:,j),1)
!            do ii=1,NDimX
!               EigY1(ii,j) = EigY1(ii,j) + val*EigY(ii,i)
!            enddo
         endif
      enddo
   endif
enddo

print*, 'EigY1',norm2(EigY1)

call clock('Y01CAS:EigY1',Tcpu,Twall)
call clock('START',Tcpu,Twall)

!print*, 'FIRST',norm2(EigY1)
! second loop ...
call dgemm('N','T',NDimX,NDimX,NDimX,1d0,EigY,NDimX,EigY1,NDimX,0d0,ABPLUS,NDimX)
!print*, 'ABPLUS-MY',norm2(ABPLUS)

call clock('Y01CAS:dgemm',Tcpu,Twall)
call clock('START',Tcpu,Twall)

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo
!
! energy loop
EAll = 0
EIntra = 0
open(newunit=iunit,file='FOFOERF',status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

kl = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work1(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints(i,j) = work1((j-1)*NBasis+i)
           enddo
        enddo
        ints(:,NOccup+1:NBasis) = 0

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i
                Crs = C(ir)+C(is)
                Cpq = C(ip)+C(iq)

                Aux = Crs*Cpq*ABPLUS(ipq,irs)
                EAll = EAll + Aux*ints(j,i)

                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))==1) EIntra = EIntra + Aux*ints(j,i)

              endif
           enddo
        enddo

      endif
   enddo
enddo

close(iunit)

ECorr = EAll-EIntra

print*, 'EAll,EIntra',EAll,EIntra

endif

!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!HNO=transpose(HNO)
!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
!write(LOUT,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)') ETot
!
!call sq_to_triang2(AuxI,work1,NBasis)
!write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!AuxI=transpose(AuxI)
!call sq_to_triang2(AuxI,work1,NBasis)
!write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
!call sq_to_triang2(AuxIO,work1,NBasis)
!write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!AuxIO=transpose(AuxIO)
!call sq_to_triang2(AuxIO,work1,NBasis)
!write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
!write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)


deallocate(RDM2val)
deallocate(ints,work2,work1)
deallocate(Eig1,EigY1)
deallocate(Eig,EigY)

call clock('Y01CASLR:ENE',Tcpu,Twall)

end subroutine Y01CASLR_FOFO

subroutine ModABMin_Act_FOFO(Occ,SRKer,Wt,OrbGrid,ABMin,&
                         MultpC,NSymNO,posAA,&
                         IndN,IndX,NDimX,NGrid,NBasis,&
                         NOccup,NAct,INActive,nAA,twokfile,twokerf)
! ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO ABMIN
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid,NOccup,NAct,INActive,nAA
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),&
                      MultpC(15,15),NSymNO(NBasis),&
                      posAA(NAct*(NAct-1)/2)
character(*),intent(in) :: twokfile,twokerf
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(inout) :: ABMin(nAA,nAA)

integer :: offset,batchlen,iunit1,iunit2
integer :: i,j,k,l,kl,ip,iq,ir,is,irs,ipq,igrd
integer :: IRow,ICol
integer :: i1i2s,i3i4s,iSym
integer :: pos(NBasis,NBasis)
double precision :: XKer1234,TwoSR,Cpq,Crs
double precision :: CICoef(NBasis)
double precision,allocatable :: work1(:),work2(:),WtKer(:)
double precision,allocatable :: batch(:,:),ABKer(:,:)
double precision,allocatable :: ints1(:,:),ints2(:,:)

pos=0
do i=1,nAA
   j=posAA(i)
   pos(IndN(1,j),IndN(2,j)) = i
enddo

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints1(NBasis,NBasis),ints2(NBasis,NBasis),&
         WtKer(maxlen),batch(maxlen,NBasis),ABKer(nAA,nAA))

ABKer = 0
do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   WtKer(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,nAA
      ipq = posAA(IRow)
      ip = IndN(1,ipq)
      iq = IndN(2,ipq)
   
      do ICol=1,nAA
         irs = posAA(ICol) 
         ir=IndN(1,irs)
         is=IndN(2,irs)

         if(irs.gt.ipq) cycle
 
         i1i2s = MultpC(NSymNO(ip),NSymNO(iq))
         i3i4s = MultpC(NSymNO(ir),NSymNO(is))
         iSym = MultpC(i1i2s,i3i4s)
         if(iSym/=1) cycle

         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + WtKer(i)* &
            batch(i,ip)*batch(i,iq)*batch(i,ir)*batch(i,is)
         enddo
         
         ABKer(IRow,ICol) = ABKer(IRow,ICol) + XKer1234
         ABKer(ICol,IRow) = ABKer(IRow,ICol)
      
      enddo
   enddo

enddo

open(newunit=iunit1,file=trim(twokfile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)
open(newunit=iunit2,file=trim(twokerf),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

kl  = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(k>INActive.and.(l>INActive.and.l<=NOccup).and.pos(l,k)/=0) then
        irs = pos(l,k) 
        ir = l
        is = k
        read(iunit1,rec=kl) work1(1:NBasis*NOccup)
        read(iunit2,rec=kl) work2(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints1(i,j) = work1((j-1)*NBasis+i)
              ints2(i,j) = work2((j-1)*NBasis+i)
           enddo
        enddo
        ints1(:,NOccup+1:NBasis) = 0
        ints2(:,NOccup+1:NBasis) = 0
  
        do j=1,NBasis
           do i=1,j 
              if((j>INActive.and.j<=NOccup).and.i>INActive.and.pos(j,i)/=0) then
                ipq = pos(j,i)
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)
                if(irs.gt.ipq) cycle

                TwoSR = ints1(j,i)-ints2(j,i)

                ABMIN(ipq,irs) = ABMIN(ipq,irs) &
                               + 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))
                ABMIN(irs,ipq) = ABMIN(ipq,irs) 

              endif  
           enddo
        enddo

      endif 
   enddo
enddo 

close(iunit2)
close(iunit1)

deallocate(ABKer,batch,WtKer,ints2,ints1,work2,work1)

end subroutine ModABMin_Act_FOFO

subroutine ModABMin_FOFO(Occ,SRKer,Wt,OrbGrid,ABMin,&
                         MultpC,NSymNO,&
                         IndN,IndX,NDimX,NGrid,NBasis,&
                         NAct,INActive,twokfile,twokerf,&
                         AB1,dfile)
! ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
use timing
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid,NAct,INActive
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),&
                      MultpC(15,15),NSymNO(NBasis)
character(*),intent(in) :: twokfile,twokerf
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
logical,intent(in) :: AB1
double precision,intent(inout) :: ABMin(NDimX,NDimX)
character(*),intent(in),optional :: dfile

integer :: offset,batchlen,iunit1,iunit2
integer :: iunit3,iunit4
integer :: i,j,k,l,kl,ip,iq,ir,is,irs,ipq,igrd
integer :: IRow,ICol
integer :: i1i2s,i3i4s,iSym
double precision :: XKer1234,TwoSR,Cpq,Crs
integer :: pos(NBasis,NBasis)
integer :: NOccup,IGem(NBasis),AuxCoeff(3,3,3,3)
double precision :: CICoef(NBasis)
double precision,allocatable :: work1(:),work2(:),WtKer(:)
double precision,allocatable :: batch(:,:),ABKer(:,:)
double precision,allocatable :: ints1(:,:),ints2(:,:)
logical :: dump
double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

dump = .false.
if(present(dfile)) dump=.true. 
print*, dfile

NOccup = NAct + INActive
! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxCoeff = 1
if(AB1) then
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               if((i==j).and.(j==k).and.(k==l)) then
                  AuxCoeff(i,j,k,l) = 0
               else
                  AuxCoeff(i,j,k,l) = 1
               endif
            enddo
         enddo
      enddo
   enddo
endif

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo


allocate(work1(NBasis**2),work2(NBasis**2),ints1(NBasis,NBasis),ints2(NBasis,NBasis),&
         WtKer(maxlen),batch(maxlen,NBasis),ABKer(NDimX,NDimX))

ABKer = 0

!print*, 'ModABMin_FOFO'
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
 
         i1i2s = MultpC(NSymNO(ip),NSymNO(iq))
         i3i4s = MultpC(NSymNO(ir),NSymNO(is))
         iSym = MultpC(i1i2s,i3i4s)
         if(iSym/=1) cycle

         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + WtKer(i)* &
                       AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))* & 
                       batch(i,ip)*batch(i,iq)*batch(i,ir)*batch(i,is)
         enddo
         
         ABKer(ipq,irs) = ABKer(ipq,irs) + XKer1234
         ABKer(irs,ipq) = ABKer(ipq,irs)
      
      enddo


   enddo

enddo

call clock('ModABMin:Ker',Tcpu,Twall)
call clock('START',Tcpu,Twall)

open(newunit=iunit1,file=trim(twokfile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)
open(newunit=iunit2,file=trim(twokerf),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

if(dump) open(newunit=iunit3,file=trim(dfile),form='unformatted')

kl = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit1,rec=kl) work1(1:NBasis*NOccup)
        read(iunit2,rec=kl) work2(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints1(i,j) = work1((j-1)*NBasis+i)
              ints2(i,j) = work2((j-1)*NBasis+i)
           enddo
        enddo
        ints1(:,NOccup+1:NBasis) = 0
        ints2(:,NOccup+1:NBasis) = 0

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)
                !if(irs.gt.ipq) cycle

                TwoSR = ints1(j,i)-ints2(j,i)

                ABMIN(ipq,irs) = ABMIN(ipq,irs) & 
                               + 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))
                !ABMIN(irs,ipq) = ABMIN(ipq,irs) 
                if(dump) write(iunit3) 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))

              endif  
           enddo
        enddo

      endif 
   enddo
enddo 

if(dump) close(iunit3)
close(iunit2)
close(iunit1)

deallocate(ABKer,batch,WtKer,ints2,ints1,work2,work1)

call clock('ModABMin',Tcpu,Twall)

end subroutine ModABMin_FOFO

subroutine ModABMinSym(Occ,SRKer,Wt,OrbGrid,TwoNO,TwoElErf,ABMin,&
           MultpC,NSymNO,IndN,IndX,NDimX,NGrid,NInte2,NBasis)
!     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid,NInte2
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX), &
                      MultpC(15,15),NSymNO(NBasis)
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(in) :: TwoNO(NInte2),TwoElErf(NInte2)
double precision,intent(inout) :: ABMin(NDimX,NDimX)

double precision :: CICoef(NBasis)
double precision,allocatable :: work(:),batch(:,:),ABKer(:,:)
integer :: i,j,IRow,ICol,ia,ib,iab,ic,id,icd
integer :: i1i2s,i3i4s,iSym
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
    
         i1i2s = MultpC(NSymNO(ia),NSymNO(ib))
         i3i4s = MultpC(NSymNO(ic),NSymNO(id))
         iSym = MultpC(i1i2s,i3i4s)
         if(iSym/=1) cycle

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

end subroutine ModABMinSym

subroutine ACEInteg_FOFO(ECorr,URe,Occ,XOne,UNOAO,&
      ABPLUS,ABMIN,EigVecR,Eig,&
      EGOne,NGOcc,&
      NBasis,NInte1,NDim,NGem,IndAux,ACAlpha,&
      IGemIN,NAct,INActive,IndN,IndX,NDimX,&
      NoSt,ICASSCF,IFlFrag1,IFunSR,IFunSRKer)
!
!  A ROUTINE FOR COMPUTING AC INTEGRAND
!
implicit none
integer,intent(in) :: NGOcc,NBasis,NInte1,NDim,NGem,NDimX
integer,intent(in) :: NAct,INActive,NoSt,ICASSCF,IFlFrag1,IFunSR,IFunSRKer 
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),&
                      IGemIN(NBasis)
double precision,intent(in) :: ACAlpha
double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),&
                               UNOAO(NBasis,NBasis),XONe(NInte1) 
double precision,intent(out) :: ABPLUS(NDim,NDim),ABMIN(NDim,NDim),&
                                EigVecR(NDimX,NDimX),Eig(NDimX)
double precision,intent(inout) :: ECorr,EGOne(NGem)

integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: pos(NBasis,NBasis)
double precision :: ECASSCF,XKer
character(:),allocatable :: twojfile,twokfile
  
 NOccup = NAct + INActive

 if(IFunSR.Ne.0.And.IFunSRKer.Eq.1) then
    twojfile = 'FFOOERF'
    twokfile = 'FOFOERF'
 else
    twojfile = 'FFOO'
    twokfile = 'FOFO'
 endif

 if(IFunSR.Ne.0.And.IFunSRKer.Eq.1) then
    pos = 0
    do ia=1,NDimX
       pos(IndN(1,ia),IndN(2,ia)) = IndX(ia)
    enddo
 endif

 if(ICASSCF==1) then

    call AB_CAS_FOFO(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne, &
                   IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDimX,&
                   NInte1,twojfile,twokfile,ACAlpha,.false.)
    EGOne(1)=ECASSCF

 elseif(ICASSCF==0) then

    write(LOUT,'(1x,a)') 'ERROR! ACABMAT0 NOT AVAILABLE FOR FOFO!'
    stop

 endif

 if(IFlFrag1.Eq.1) then
    write(LOUT,'(1x,a)') 'ERROR! FRAGMENT CALCULATIONS NOT AVAILABLE FOR FOFO!'
    stop
 endif

! ADD A SR KERNEL
 if(IFunSR.Ne.0.And.IFunSRKer.Eq.1) then

    open(newunit=iunit,file="srdump",form='UNFORMATTED')

    kl = 0
    do k=1,NOccup
       do l=1,NBasis
          kl = kl + 1
          if(pos(l,k)/=0) then
            irs = pos(l,k)
            ir = l
            is = k
    
            do j=1,NBasis
               do i=1,j
                  if(pos(j,i)/=0) then
                    ipq = pos(j,i)
                    ip = j
                    iq = i
   
                    read(iunit) XKer
                    if(.not.(IGemIN(ir).Eq.IGemIN(is).And.&
                             IGemIN(is).Eq.IGemIN(ip).And.&
                             IGemIN(ip).Eq.IGemIN(iq))) XKer=XKer*ACAlpha
                   
                    ABMIN(ipq,irs) = ABMIN(ipq,irs) & 
                                   + XKer
   
                  endif  
               enddo
            enddo
    
          endif 
       enddo
    enddo 

    close(iunit)

 endif
!
! FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY

 if(NoSt==1) then
    call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
   else
    call ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
 endif

 call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Occ, &
                     IGemIN,IndN,IndX,INActive+NAct, &
                     NDimX,NBasis,twokfile)

end subroutine ACEInteg_FOFO

subroutine ACEneERPA_FOFO(ECorr,EVec,EVal,Occ,IGem, &
                          IndN,IndX,NOccup,NDimX,NBasis,IntKFile)
implicit none

integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: IGem(NBasis),IndN(2,NDimX),IndX(NDimX)
integer,intent(in) :: NOccup
character(*),intent(in) :: IntKFile
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

open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

kl = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints(i,j) = work((j-1)*NBasis+i)
           enddo
        enddo
        ints(:,NOccup+1:NBasis) = 0

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

end subroutine ACEneERPA_FOFO

subroutine dump_Eblock(Eblock,EblockIV,Occ,IndN,nblk,NBasis,NDimX,xy0file)
implicit none

integer,intent(in)           :: nblk,NBasis,NDimX
integer,intent(in)           :: IndN(2,NDimX)
character(*),intent(in)      :: xy0file
double precision,intent(in)  :: Occ(NBasis)
type(EBlockData),intent(in)  :: EBlock(nblk),EBlockIV
type(EBlockData),allocatable :: SBlock(:)
type(EBlockData)             :: SBlockIV

integer                      :: iunit
integer                      :: i,j,k,ii,ipos,iblk,ip,iq
double precision             :: C(NBasis),fac,valY,valX
double precision,allocatable :: EigY(:,:),EigX(:,:),Eig(:)
double precision,allocatable :: work(:,:)

fac = 1d0/sqrt(2d0)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

allocate(SBlock(nblk))

! repack
do iblk=1,nblk
   associate(B => Eblock(iblk),&
             A => Sblock(iblk))

     A%n  = B%n
     A%l1 = B%l1
     A%l2 = B%l2
     allocate(A%pos(A%n),A%vec(A%n))
     A%pos(1:A%n)=B%pos(1:B%n)
     A%vec(1:A%n)=B%vec(1:B%n)

     allocate(A%matX(A%n,A%n),A%matY(A%n,A%n))
     do i=1,B%n
        ipos = B%pos(i)
        ip = IndN(1,ipos)
        iq = IndN(2,ipos)

        valX = 1d0/(C(ip)+C(iq))
        valY = 1d0/(C(ip)-C(iq))
        A%matX(i,1:B%n) = 0.5d0*(B%matX(i,1:B%n)*valX - B%matY(i,1:B%n)*valY)
        A%matY(i,1:B%n) = 0.5d0*(B%matX(i,1:B%n)*valX + B%matY(i,1:B%n)*valY)

     enddo

   end associate
enddo
! repack IV
associate(B => EblockIV, &
          A => SBlockIV )

  A%n  = B%n
  A%l1 = B%l1
  A%l2 = B%l2
  allocate(A%pos(A%n))
  A%pos(1:A%n) = B%pos(1:B%n)

  allocate(A%vec(A%n),A%matX(A%n,1),A%matY(A%n,1))
  do i=1,A%n
     ipos = B%pos(i)
     ip = IndN(1,ipos)
     iq = IndN(2,ipos)

     valX = 1d0/(C(ip)+C(iq))
     valY = 1d0/(C(ip)-C(iq))

     A%matX(i,1) = 0.5d0*fac*(valX-valY)
     A%matY(i,1) = 0.5d0*fac*(valX+valY)

     A%vec(i) = B%vec(i)
  enddo

end associate

! dump to a file
open(newunit=iunit,file=xy0file,form='unformatted')
write(iunit) nblk
do iblk=1,nblk
   associate(B => SBlock(iblk))
     write(iunit) iblk, B%n, B%l1, B%l2
     write(iunit) B%pos,B%matX,B%matY,B%vec
   end associate
enddo
associate(B => SBlockIV)
  write(iunit) B%n,B%l1,B%l2
  write(iunit) B%pos,B%vec
  write(iunit) B%matX,B%matY
end associate
close(iunit)
! deallocate blocks
do iblk=1,nblk
   associate(B => Sblock(iblk))
     deallocate(B%matY,B%matX,B%vec)
     deallocate(B%pos)
   end associate
enddo
! deallocate IV blocks 
associate(B => SblockIV)
  deallocate(B%vec)
  deallocate(B%pos)
end associate

deallocate(Sblock)

end subroutine dump_Eblock

subroutine check_mp2(EMP2,IntKFIle,INActive,NOccup,NBasis)
! check MP2 energy (only inactive-vritual enter) 
! to compare with molpro use {mp2;core,0}    
implicit none

integer,intent(in) :: INActive,NOccup,NBasis
character(*),intent(in) :: IntKFile
double precision,intent(inout) :: EMP2

integer :: i,j,k,l,a,b
integer :: iunit
double precision :: val,eps(NBasis,NBasis)
double precision,allocatable :: work(:),ints_bj(:,:),ints_bi(:,:)

allocate(work(NBasis**2),ints_bj(NBasis,NBasis),ints_bi(NBasis,NBasis))

eps = 0
open(10,file='fock.dat')
do i=NOccup+1,NBasis
   do j=1,INActive
      read(10,*) k,l,val
      eps(k,l) = val
   enddo
enddo
close(10)

open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

EMP2 = 0
do j=1,INActive
   do b=NOccup+1,NBasis
     
      read(iunit,rec=(b+(j-1)*NBasis)) work(1:NBasis*NOccup)
      ! ints_bj
      do l=1,NOccup
         do k=1,NBasis
            ints_bj(k,l) = work((l-1)*NBasis+k)
         enddo
      enddo
      ints_bj(:,NOccup+1:NBasis) = 0

      do i=1,INActive
         read(iunit,rec=(b+(i-1)*NBasis)) work(1:NBasis*NOccup)
         ! ints_bi
         do l=1,NOccup
            do k=1,NBasis
               ints_bi(k,l) = work((l-1)*NBasis+k)
            enddo
         enddo
         ints_bi(:,NOccup+1:NBasis) = 0

         do a=NOccup+1,NBasis
            EMP2 = EMP2 - ints_bj(a,i)*(2d0*ints_bj(a,i)-ints_bi(a,j)) / &
                          (eps(a,i)+eps(b,j))
         enddo
      enddo

   enddo
enddo

write(LOUT,'(/,1x,a,13x,f16.8,/)') 'EMP2   Energy', EMP2 

deallocate(ints_bi,ints_bj,work)

end subroutine check_mp2

subroutine Y01CASD_FOFO(IH0St,Occ,URe,XOne, &
     propfile0,propfile1, &
     xy0file, &
     UNOAO,IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     NoSt,IntFileName,IntJFile,IntKFile,ETot,ECorr)
!
!     A ROUTINE FOR COMPUTING AC0D CORRELATION ENERGIES AND
!     TRANSITION DIPOLE MOMENTS IN THE 0- AND 1 - ORDER APPROXIMATIONS
!
use timing
!
implicit none
double precision :: UNOAO(NBasis,NBasis),DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1,NoSt
character(*) :: IntFileName,IntJFile,IntKFile
character(*) :: propfile0,propfile1,xy0file
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
double precision,intent(out),optional :: ETot,ECorr

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: iunit1,iunit2
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,val
double precision :: EnDummy,Aux,Crs,Cpq,EIntra,EAll
double precision :: EMP2
double precision :: work(NBasis,NBasis)
double precision,allocatable :: EigY(:,:),EigY1(:,:),EigX(:,:)
integer,allocatable :: iaddr(:)
double precision,allocatable :: Eig(:),Eig1(:)
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
!double precision,allocatable :: EigYt(:,:),EigXt(:,:),Eigt(:)
double precision,allocatable :: ABP(:,:),ABM(:,:),workA(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2
integer :: ii,jj,ipos,jpos,iblk,jblk,nblk
integer :: nAA,tmpAA(NAct*(NAct-1)/2),limAA(2)
integer :: nAI(INActive),tmpAI(NAct,1:INActive),limAI(2,1:INActive)
integer :: nAV(INActive+NAct+1:NBasis),tmpAV(NAct,INActive+NAct+1:NBasis),limAV(2,INActive+NAct+1:NBasis)
integer :: nIV,tmpIV(INActive*(NBasis-NAct-INActive)),limIV(2)
integer :: IH0St,IStERPA,NoStMx,IPr, IPair(NBasis,NBasis),NU,IAB,IA,IB,NUMx
double precision :: GammaS(100,NInte1),XCAS(NBasis,NInte1),YCAS(NBasis,NInte1),EExcit(NInte1)
double precision :: OvMax,OvXMax,OvYMax,SumERY,SumCAY,SumERX,SumCAX,YER,XER,YCA,XCA,SOvY,SOvX,TDIP2,EigOv
double precision,allocatable :: ABPLUS(:,:),ABMIN(:,:)
!
type(EblockData),allocatable :: Eblock(:) 
type(EblockData) :: EblockIV
!!
double precision,parameter :: Thresh = 1.D-12
double precision :: tmpEn
!! test timings
double precision :: Tcpu,Twall

! timing 
call clock('START',Tcpu,Twall)

IStERPA=0
IPair(1:NBasis,1:NBasis)=0
Do II=1,NDimX
  I=IndN(1,II)
  J=IndN(2,II)
  IPair(I,J)=1
  IPair(J,I)=1
EndDo

Do I=1,NBasis
   C(I)=SQRT(Occ(I))
   If(Occ(I).Lt.0.5) C(I)=-C(I)
EndDo

If(IH0St.Lt.NoSt) Then
  Write(6,'(/,X,"IH0St=",I4,"  NoSt=",I4)')IH0St,NoSt
  Stop 'Stop in AC0CDEXCIT: IH0St is smaller than NoSt. AC0 correction would not be reliable '
EndIf

Call ReadDip(DipX,DipY,DipZ,UNOAO,NBasis)

NoStMx=0
Write(6,'(X,"**** SA-CAS FROM MOLPRO ****",/)')
IPr=0
If(IH0St.Eq.NoSt) IPr=1
Call RDM_SACAS(GammaS,XCAS,YCAS,EExcit,C,UNOAO,IPair,DipX,DipY,DipZ,NoSt,NoStMx,NInte1,NBasis,IPr)
Write(6,'(X,"The number of states in SA-CAS: ",I4,/)')NoStMx
If(IH0St.Gt.NoStMx) Stop'Stop in AC0CDEXCIT: IH0St is greater than the number of states in SA-CAS!'
Do NU=1,NoSt-1
  Write(6,'(X,"SA-CAS Deexcitation from ",I4," to",I4,2E15.6)')NoSt,NU,EExcit(NU)
EndDo
Do NU=NoSt+1,NoStMx
  Write(6,'(X,"SA-CAS Excitation   from ",I4," to",I4,2E15.6)') NoSt,NU,EExcit(NU)
EndDo

If(IH0St.Ne.NoSt.And.EExcit(IH0St).Eq.0.0) Then
  Write(6,'(/,X,"Transition energy of interest is zero, quitting")')
  ECorr=0.0
  Return
EndIf

allocate(ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX))
ABPLUS = 0
ABMIN  = 0
if(present(ETot)) ETot = 0

! test EUGENE
 tmpEn = 0

! set dimensions
NOccup = NAct + INActive

Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
tmpEn = tmpEn + 2*val
print*, 'ONE ELECTRON ENERGY:', tmpEn
if(present(ETot)) ETot = ETot + 2*val

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = 0d0
   enddo
enddo

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = 0
            endif
         enddo
      enddo
   enddo
enddo

nAA = 0
nAI = 0
nAV = 0
nIV = 0
tmpAA = 0
tmpAI = 0
tmpAV = 0
tmpIV = 0
limAA = 0
limAI = 0
limAV = 0
limIV = 0

do i=1,NDimX
   ip = IndN(1,i)
   iq = IndN(2,i)
   if(IGem(ip)==2.and.IGem(iq)==2) then
      nAA = nAA + 1
      tmpAA(nAA) = i
    elseif(IGem(ip)==2.and.IGem(iq)==1) then
       nAI(iq) = nAI(iq) + 1
       tmpAI(nAI(iq),iq) = i
    elseif(IGem(ip)==3.and.IGem(iq)==2) then
       nAV(ip) = nAV(ip) + 1
       tmpAV(nAV(ip),ip) = i
    elseif(IGem(ip)==3.and.IGem(iq)==1) then
       nIV = nIV + 1
       tmpIV(nIV) = i
   endif
enddo

pos = 0
ipos = 0
limAA(1) = ipos + 1
do ii=1,nAA
   i = tmpAA(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limAA(2) = ipos
do is=1,INActive
   limAI(1,is) = ipos + 1
   do ii=1,nAI(is)
      i = tmpAI(ii,is)
      ! print*, ii,is,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAI(2,is) = ipos
enddo
!print*, ' '
do ir=NOccup+1,NBasis
   limAV(1,ir) = ipos + 1
   do ii=1,nAV(ir)
      i = tmpAV(ii,ir)
      ! print*, ii,ir,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAV(2,ir) = ipos
enddo
limIV(1) = ipos + 1
do ii=1,nIV
   i = tmpIV(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limIV(2) = ipos

!not working... ???
!call AB_CAS_assemb(ABPLUS,ABMIN,Occ,RDM2val,HNO,pos, &
!           1d0,AuxInd,AuxCoeff, &
!           IGem,IndN,IndX,NDimX,NDim,INActive,NOccup,NBasis, &
!           IntJFile,IntKFile)

AuxI  = 0
AuxIO = 0
WMAT  = 0

Write(6,'(/,X,"*** H0 constructed for SA-CAS state no",I4," ***")') NoSt

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do ll=1,NOccup
   do kk=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      k = kk
      l = ll

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = Occ(l)
         if(IGem(l)==1) then
            if(IGem(k)==2) then
               do i=INActive+1,NOccup
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
         
         if(IGem(l)==2) then
            if(IGem(k)==1) then
               do i=1,INActive
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO
      
      ! exchange 
      if(l<=INActive) then
         do i=1,NBasis
            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif
      if(l<=NOccup) then
         do i=1,NBasis
            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif

! AUXILIARY MATRIX WMAT
         if(l<=NOccup) then
            do ir=1,NOccup
               val = 0
               val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                    sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                    sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                    sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                    sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
               WMAT(k,ir) = WMAT(k,ir) + val
            enddo
         endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      ! Coulomb
      if(k>INActive.and.l<=NOccup) then
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then
            
            do iq=1,NOccup
               do ip=INActive+1,NBasis
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                     
                     val = 0
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                  
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     val = 0
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            enddo

         endif
      endif

      ! exchange
      if(l<=NOccup) then
         if(k>INActive) then
            ip = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then                  
                  do iq=1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = - AuxVal*val*ints(ir,iq)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               endif
            enddo
         endif
      endif

      ! T1 p->q
      if((k>INActive).and.(l<=NOccup)) then
         ip = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T3
      if((l<=NOccup).and.(k>INActive)) then
         iq = l
         ir = k
         
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T4+T6: AuxInd = 2
      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
         ir = k
         iu = l

         ! T4
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
         ! T6
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then                        
                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then  
                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif

   enddo
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do ll=1,NOccup
   do kk=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      
      k = kk
      l = ll

      if(k>NOccup.or.l>NOccup) cycle

      if(present(ETot)) then  
      ! COMPUTE THE ENERGY FOR CHECKING
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))
      endif

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*Occ(k)
         if(k>INActive) then 
            do j=1,INActive
               do i=1,INActive
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         if(k<=INActive) then
            do j=INActive+1,NOccup
               do i=INActive+1,NOccup
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         do j=NOccup+1,NBasis
            do i=NOccup+1,NBasis
               HNO(i,j) = HNO(i,j) + val*ints(i,j)
            enddo
         enddo
         
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! Coulomb
      if(k==l) then
         val = 2*Occ(k)
         if(k<=INActive) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
               enddo
            enddo
         endif
         if(k>INActive.and.k<=NOccup) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
               enddo
            enddo
         endif
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      
      ! exchange
      if(l<=NOccup) then
         if(k<=NOccup) then
            iq = k
            is = l
            
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = - AuxVal*val*ints(ip,ir)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               endif
            enddo
               
         endif
      endif

      ! T1 
      if((k<=NOccup).and.(l<=NOccup)) then
         iq = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo
            
      endif

      ! T3
      if(IGem(k)==2.and.IGem(l)==2) then
         
         do is=INActive+1,NOccup
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           
                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                RDM2val(iq,is,k,l)*ints(ip,ir)
                           
                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif

      ! T2+T5: AuxInd = 2
      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
         is = k
         iu = l
         
         ! T2
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo

         ! T5
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif

   enddo
enddo

close(iunit2)
deallocate(work1)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = pos(ir,is)
      
      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = pos(ip,iq)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val
         
      enddo
enddo

call clock('AB0MAT',Tcpu,Twall)

allocate(EBlock(1+NBasis-NAct))

nblk = 0

! put pack into separate subroutine
!pack AA
if(nAA>0) then

   nblk = nblk + 1
   associate(B => Eblock(nblk))
     B%n  = nAA
     B%l1 = limAA(1)
     B%l2 = limAA(2)
     allocate(B%pos(B%n))
     B%pos(1:B%n) = tmpAA(1:B%n)

     allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
     allocate(ABP(B%n,B%n),ABM(B%n,B%n))

     ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
     ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
     if(NoSt==1) then
        call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
     elseif(NoSt>1) then
        call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
     endif
   
     deallocate(ABM,ABP)
   end associate

endif

write(*,'(/,1x,a)') 'Act-Act   block diagonalized!'

allocate(Eig(nAA),EigY(nAA,nAA),EigX(nAA,nAA),iaddr(nAA))
Eig = 0
EigY = 0
! nAA to liczba aktywnych
if(nAA > 0) then
   associate(B => Eblock(1))
      do i=1,B%n

           If(B%vec(i).lt.0) Then
               Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)') B%vec(i)
               B%vec(i)=0
               B%matY(1:B%n,i)=0
               B%matX(1:B%n,i)=0
           EndIf

           iaddr(i)=B%pos(i)
           EigY(i,B%l1:B%l2) = B%matY(i,1:B%n)
           EigX(i,B%l1:B%l2) = B%matX(i,1:B%n)

     enddo
     Eig(B%l1:B%l2) = B%vec(1:B%n)
   end associate

!   if(NoSt>1) then
!     Call SortEigXY(1,Eig,EigY,EigX,nAA)
!   endif

   Write(6,'(X,"Active ERPA Eigenvalues and Eigenvecs")')
   do i=1,nAA
      Write(6,'(X,I4,E15.6)') i,eig(i)
      do iq=1,nAA
         if(abs(eigy(iq,i))+abs(eigx(iq,i)).gt.1.d-7) Write(6,'(X,"Y_ERPA, X_ERPA",2I3,2E15.6)') & 
         IndN(1,iaddr(iq)),IndN(2,iaddr(iq)),eigy(iq,i),eigx(iq,i)
      enddo
   enddo

endif

If(IH0St.Ne.NoSt) Then

      Write(6,*)
      OvMax=0
      OvXMax=0
      OvYMax=0

      Do NU=1,nAA

      SumERY=0
      SumERX=0
      SumCAY=0
      SumCAX=0

      Do I=1,nAA
      IA=IndN(1,iaddr(i))
      IB=IndN(2,iaddr(i))
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)

      YER=eigy(i,NU)
      YCA=YCAS(IH0St,IAB)
      XER=eigx(i,NU)
      XCA=XCAS(IH0St,IAB)

      SumERX=SumERX+XER*XER
      SumCAX=SumCAX+XCA*XCA
      SumERY=SumERY+YER*YER
      SumCAY=SumCAY+YCA*YCA
      EndDo

      If(Abs(SumERX).Gt.1.D-12) SumERX=1.D0/Sqrt(Abs(SumERX))
      If(Abs(SumCAX).Gt.1.D-12) SumCAX=1.D0/Sqrt(Abs(SumCAX))
      If(Abs(SumERY).Gt.1.D-12) SumERY=1.D0/Sqrt(Abs(SumERY))
      If(Abs(SumCAY).Gt.1.D-12) SumCAY=1.D0/Sqrt(Abs(SumCAY))

      SOvY=0
      SOvX=0

      Do I=1,nAA
      IA=IndN(1,iaddr(i))
      IB=IndN(2,iaddr(i))
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      YER=eigy(i,NU)
      YCA=YCAS(IH0St,IAB)
      XER=eigx(i,NU)
      XCA=XCAS(IH0St,IAB)
      SOvY=SOvY+YER*YCA*SumERY*SumCAY
      SOvX=SOvX+XER*XCA*SumERX*SumCAX
      EndDo

      SOvY=Abs(SOvY)
      SOvX=Abs(SOvX)
      Write(6,'(X,"SA-CAS-Overlap for ERPA vector",I2,3F15.8)') NU,SOvY,SOvX,(SOvY+SovX)/2.D0

      If((SOvY+SovX)/2.D0.Gt.OvMax) Then
      OvMax=(SOvY+SovX)/2.D0
      OvXMax=SOvX
      OvYMax=SOvY
      NUMx=NU
      EndIf

      EndDo 

      If((OvMax.Ge.0.5).Or. &
       ( (OvMax.Lt.0.5).And.(OvXMax.Gt.0.5.Or.OvYMax.Gt.0.5) )) Then

         IStERPA=NUMx
         EigOv=Eig(IStERPA)

         Write(6,'(/,X,"ERPA vector best overlapping with SA-CAS is vector no", &
         I2," of excit energy:",F15.8)') IStERPA,Eig(IStERPA)
         NU=IStERPA

         Do I=1,nAA
         IA=IndN(1,iaddr(i))
         IB=IndN(2,iaddr(i))
         If(Abs(eigy(i,NU))+Abs(eigx(i,NU)).Gt.1.d-7) &
         Write(6,'(X,"Y_ERPA, X_ERPA        ",2I3,2E15.6)')IA,IB,eigy(i,NU),eigx(i,NU)
         EndDo

         Write(6,'(X,"best-matching SA-CAS Excitation from ", &
               I4," to ",I4,23X,F15.8)') NoSt,IH0St,EExcit(IH0St)

         Do I=1,nAA
         IA=IndN(1,iaddr(i))
         IB=IndN(2,iaddr(i))
         IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
         If(Abs(YCAS(IH0St,IAB))+Abs(XCAS(IH0St,IAB)).Gt.1.d-7) &
         Write(6,'(X,"Y_SA-CAS, X_SA-CAS    ",2I3,2E15.6)')IA,IB, &
         YCAS(IH0St,IAB),XCAS(IH0St,IAB)
         EndDo

      Else

         IStERPA=0
         Write(6,'(/,X, "No ERPA vector overlaps with SA-CAS State No",I2)')IH0St

      EndIf


EndIf
deallocate(Eig,EigY,EigX,iaddr)

!pack AI
do iq=1,INActive
   if(nAI(iq)>0) then

      nblk = nblk + 1
      associate(B => Eblock(nblk))
        B%n  = nAI(iq)
        B%l1 = limAI(1,iq)
        B%l2 = limAI(2,iq)
        allocate(B%pos(B%n))
        B%pos(1:B%n) = tmpAI(1:B%n,iq)

        allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
        allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  
        ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
        ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
        if(NoSt==1) then
           call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        elseif(NoSt>1) then
           call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        endif

        deallocate(ABM,ABP)
      end associate

   endif
enddo

print*
print*, 'Act-InAct block diagonalized!'

!pack AV
do ip=NOccup+1,NBasis
   if(nAV(ip)>0) then

      nblk = nblk + 1
      associate(B => Eblock(nblk))
        B%n  = nAV(ip)
        B%l1 = limAV(1,ip)
        B%l2 = limAV(2,ip)
        allocate(B%pos(B%n))
        B%pos(1:B%n) = tmpAV(1:B%n,ip)

        allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
        allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  
        ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
        ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
        if(NoSt==1) then
           call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        elseif(NoSt>1) then
           call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        endif

        deallocate(ABM,ABP)
      end associate

    endif
enddo

print*, 'Act-Virt  block diagonalized!'

!pack IV
associate(B => EblockIV)

  B%l1 = limIV(1)
  B%l2 = limIV(2)
  B%n  = B%l2-B%l1+1
  allocate(B%pos(B%n))
  B%pos(1:B%n) = tmpIV(1:B%n)

  allocate(B%vec(B%n))
  do i=1,B%n
     ii = B%l1+i-1
     B%vec(i) = ABPLUS(ii,ii)
  enddo

end associate

 call clock('PACK',Tcpu,Twall)

!! dump EBLOCKS: X(0),Y(0)
call dump_Eblock(Eblock,EblockIV,Occ,IndN,nblk,NBasis,NDimX,xy0file)

   ! AB(1) PART
   call AB_CAS_FOFO(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,&
                 IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
                 NInte1,IntJFile,IntKFile,1d0,.true.)

  call clock('AB(1)',Tcpu,Twall)
!  ! B = A.X
!  ! C = X^T.B
   allocate(workA(NDimX,NDimX))

   call ABPM_TRAN(ABPLUS,workA,EBlock,EBlockIV,nblk,NDimX,.true.)
   call clock('MULT-1',Tcpu,Twall)
   ABPLUS=workA
   call ABPM_TRAN(ABMIN,workA,EBlock,EBlockIV,nblk,NDimX,.false.)
   ABMIN=workA

   call clock('MULT-2',Tcpu,Twall)

   deallocate(workA)

   allocate(Eig(NDimX),Eig1(NDimX))
 
   ! unpack Eig (1)
   do iblk=1,nblk
      associate(B => Eblock(iblk))
 
        Eig(B%l1:B%l2) = B%vec(1:B%n)

      end associate
   enddo

   !unpack Eig (2, IV part)
   associate(B => EblockIV)
 
     do i=1,B%n
        ii = B%l1+i-1
        ipos = B%pos(i)
        Eig(ii) = B%vec(i)
     enddo

   end associate
   call clock('UNPACK',Tcpu,Twall)

deallocate(RDM2val,work2)
deallocate(ints)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

If(NoSt.Eq.1) Then
! transition dipole moments in 0th and 1st order approximations
! ------------------------------------------------------------------------------

allocate(EigY(NDimX,NDimX),EigY1(NDimX,NDimX))

EigY1 = 0
EigY  = 0
! unpack (1)
do iblk=1,nblk
   associate(B => Eblock(iblk))
     do i=1,B%n
        ipos = B%pos(i)
        EigY(ipos,B%l1:B%l2) = B%matY(i,1:B%n)
     enddo
   end associate
enddo

!unpack (2, IV part)
associate(B => EblockIV)
  do i=1,B%n
     ii = B%l1+i-1
     ipos = B%pos(i)
     EigY(ipos,ii) = 1d0/sqrt(2d0)
  enddo
end associate

do i=1,NDimX
   Eig1(i)=ABPLUS(i,i)+ABMIN(i,i)
enddo

allocate(workA(NDimX,NDimX))
do j=1,NDimX
   if(Eig(j)/=0d0) then
      do i=1,NDimX
         if(Eig(i)/=0d0) then
            workA(i,j) = (ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
            if(Abs(Eig(i)-Eig(j))>Thresh) then
               workA(i,j) = workA(i,j) + (ABPLUS(i,j)+ABMIN(i,j))/(Eig(j)-Eig(i))
            endif
         endif
      enddo
   endif
enddo

call dgemm('N','N',NDimX,NDimX,NDimX,1d0,EigY,NDimX,workA,NDimX,0d0,EigY1,NDimX)

If(IStERPA.Ne.0) Then

      Write(6,'(/,X,"Y(0) corresponding to SA-CAS best-matching vector no",I2)')IStERPA
      Call TrDipMoms(IStERPA,TDIP2,EigY,C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
      Write(6,'(X,I2,"->",I2,"  Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)') NoSt,IH0St,TDIP2

      Write(6,'(X,"Y(0)+Y(1) corresponding to SA-CAS best-matching vector no",I2)')IStERPA

      Do I=1,NdimX
      EigY(I,IStERPA)=EigY(I,IStERPA)+EigY1(I,IStERPA)
      EndDo
      Call TrDipMoms(IStERPA,TDIP2,EigY,C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
      Write(6,'(X,I2,"->",I2,"  Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)') NoSt,IH0St,TDIP2

Else

      If(IH0St.Ne.NoSt) Then
      Write(6,'(X,I2,"->",I2,"  Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)') NoSt,IH0St,0.0
      Write(6,'(X,I2,"->",I2,"  Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)') NoSt,IH0St,0.0
      EndIf

EndIf

deallocate(workA)
deallocate(EigY1,EigY)

! ------------------------------------------------------------------------------
! If(NoSt.Eq.1) 
EndIf
!
! for AC0Corr
! ------------------------------------------------------------------------------

!prepare new AMAT
allocate(workA(NDimX,NDimX))
workA=0

If(IH0St.Ne.NoSt.And.IStERPA.Ne.0) & 
Write(6,'(/,X,"Deexcitation correction is computed for ERPA vector no",I2," Eig=",F15.8)') IStERPA,Eig(IStERPA)

If(IH0St.Ne.NoSt.And.IStERPA.Eq.0) &
 Write(6,'(/," ERPA vector for deexcitation correction could not be determined. The correction will be set to 0.")')

do j=1,NDimX
   if(Eig(j)/=0d0) then
      do i=1,NDimX
         If(IH0St.Eq.NoSt.Or.(I.Eq.IStERPA.Or.J.Eq.IStERPA)) Then
         if(Eig(i)/=0d0) workA(i,j) = 2d0*(ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
         EndIf 
      enddo
   endif
enddo
call clock('APL-AMIN/EPL',Tcpu,Twall)

call ABPM_BACKTRAN(workA,ABPLUS,EBlock,EBlockIV,nblk,NDimX)

deallocate(workA)

! deallocate blocks
do iblk=1,nblk
   associate(B => Eblock(iblk))

     deallocate(B%matY,B%matX,B%vec)
     deallocate(B%pos)

   end associate
enddo
! (IV part)
associate(B => EblockIV)

  deallocate(B%vec)
  deallocate(B%pos)

end associate

call clock('ABPl-dgemm',Tcpu,Twall)
!print*, 'ABPLUS-MY',norm2(ABPLUS)

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo
!
! energy loop
allocate(work1(NBasis*NBasis),ints(NBasis,NBasis))
EAll = 0
EIntra = 0
open(newunit=iunit,file='FOFO',status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

ints = 0
kl   = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work1(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints(i,j) = work1((j-1)*NBasis+i)
           enddo
        enddo
        ints(:,NOccup+1:NBasis) = 0

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i
                Crs = C(ir)+C(is)
                Cpq = C(ip)+C(iq)

                Aux = Crs*Cpq*ABPLUS(ipq,irs)
                EAll = EAll + Aux*ints(j,i)

                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))==1) EIntra = EIntra + Aux*ints(j,i)

              endif
           enddo
        enddo

      endif
   enddo
enddo

close(iunit)

call clock('ENE-loop',Tcpu,Twall)

ECorr = EAll-EIntra
!print*, 'EAll,EIntra',EAll,EIntra

deallocate(ints,work1)
deallocate(Eig1,Eig)
deallocate(ABPLUS,ABMIN)

end subroutine Y01CASD_FOFO

subroutine Y01CASDSYM_FOFO(ICAS,NoStMx,ICORR,EExcit,IStCAS,NSym,NSymNO,MultpC,ECorrSym,  &
     Occ,URe,XOne, &
     propfile0,propfile1, &
     xy0file, &
     UNOAO,IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     NoSt,IntFileName,IntJFile,IntKFile,ETot)
!
!     A ROUTINE FOR COMPUTING AC0D CORRELATION ENERGIES AND
!     TRANSITION DIPOLE MOMENTS IN THE 0- AND 1 - ORDER APPROXIMATIONS
!
use timing
!
implicit none
double precision :: UNOAO(NBasis,NBasis),DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1,NoSt
character(*) :: IntFileName,IntJFile,IntKFile
character(*) :: propfile0,propfile1,xy0file
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
double precision,intent(out),optional :: ETot

integer :: i,j,k,l,ij,kl,kk,ll,klround,ia,ib
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: iunit1,iunit2
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,val
double precision :: EnDummy,Aux,Crs,Cpq,EIntra,EAll
double precision :: EMP2
double precision :: work(NBasis,NBasis)
double precision,allocatable :: EigY(:,:),EigY1(:,:),EigX(:,:)
integer,allocatable :: iaddr(:)
double precision,allocatable :: Eig(:),Eig1(:)
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
!double precision,allocatable :: EigYt(:,:),EigXt(:,:),Eigt(:)
double precision,allocatable :: ABP(:,:),ABM(:,:),workA(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2
integer :: ii,jj,ipos,jpos,iblk,jblk,nblk
integer :: nAA,tmpAA(NAct*(NAct-1)/2),limAA(2)
integer :: nAI(INActive),tmpAI(NAct,1:INActive),limAI(2,1:INActive)
integer :: nAV(INActive+NAct+1:NBasis),tmpAV(NAct,INActive+NAct+1:NBasis),limAV(2,INActive+NAct+1:NBasis)
integer :: nIV,tmpIV(INActive*(NBasis-NAct-INActive)),limIV(2)
integer :: IERPA,NoStMx,IPr,IPair(NBasis,NBasis),NU,IAB,NUMx
!double precision :: GammaS(100,NInte1),XCAS(NBasis,NInte1),YCAS(NBasis,NInte1)
!double precision :: OvMax,OvXMax,OvYMax,SumERY,SumCAY,SumERX,SumCAX,YER,XER,YCA,XCA,SOvY,SOvX,TDIP2,EigOv
double precision,allocatable :: ABPLUS(:,:),ABMIN(:,:)
double precision :: ECorrSym(100),EExcit(100),ECorr,EigMin,Hlp,TDIP2,Eig11(100)
integer :: NSym,NSymNO(NBasis),MultpC(8,8),IStCAS(2,100),ICORR(100),IStERPA(2,100),ISym,IStart,ICount, &
           IndMin,IndHlp,ICAS,IDCORR,NoEig,IAC0,ISt11ERPA(2,100),NoEig11,NegSym(8),NoERPASym(8),IStateInSACAS(100,8)
integer,allocatable :: IZeroNU(:)
logical :: file_exists
!
type(EblockData),allocatable :: Eblock(:) 
type(EblockData) :: EblockIV
!!
double precision,parameter :: Thresh = 1.D-12
double precision :: tmpEn
!! test timings
double precision :: Tcpu,Twall

! timing 
call clock('START',Tcpu,Twall)
Do I=1,NBasis
   C(I)=SQRT(Occ(I))
   If(Occ(I).Lt.0.5) C(I)=-C(I)
EndDo

If(IStCAS(1,ICAS).Eq.1.And.IStCAS(2,ICAS).Eq.1) Call ReadDip(DipX,DipY,DipZ,UNOAO,NBasis)

allocate(ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX))
ABPLUS = 0
ABMIN  = 0
if(present(ETot)) ETot = 0

! test EUGENE
 tmpEn = 0

! set dimensions
NOccup = NAct + INActive

Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
tmpEn = tmpEn + 2*val
print*, 'ONE ELECTRON ENERGY:', tmpEn
if(present(ETot)) ETot = ETot + 2*val

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = 0d0
   enddo
enddo

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = 0
            endif
         enddo
      enddo
   enddo
enddo

nAA = 0
nAI = 0
nAV = 0
nIV = 0
tmpAA = 0
tmpAI = 0
tmpAV = 0
tmpIV = 0
limAA = 0
limAI = 0
limAV = 0
limIV = 0

do i=1,NDimX
   ip = IndN(1,i)
   iq = IndN(2,i)
   if(IGem(ip)==2.and.IGem(iq)==2) then
      nAA = nAA + 1
      tmpAA(nAA) = i
    elseif(IGem(ip)==2.and.IGem(iq)==1) then
       nAI(iq) = nAI(iq) + 1
       tmpAI(nAI(iq),iq) = i
    elseif(IGem(ip)==3.and.IGem(iq)==2) then
       nAV(ip) = nAV(ip) + 1
       tmpAV(nAV(ip),ip) = i
    elseif(IGem(ip)==3.and.IGem(iq)==1) then
       nIV = nIV + 1
       tmpIV(nIV) = i
   endif
enddo

pos = 0
ipos = 0
limAA(1) = ipos + 1
do ii=1,nAA
   i = tmpAA(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limAA(2) = ipos
do is=1,INActive
   limAI(1,is) = ipos + 1
   do ii=1,nAI(is)
      i = tmpAI(ii,is)
      ! print*, ii,is,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAI(2,is) = ipos
enddo
!print*, ' '
do ir=NOccup+1,NBasis
   limAV(1,ir) = ipos + 1
   do ii=1,nAV(ir)
      i = tmpAV(ii,ir)
      ! print*, ii,ir,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAV(2,ir) = ipos
enddo
limIV(1) = ipos + 1
do ii=1,nIV
   i = tmpIV(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limIV(2) = ipos

AuxI  = 0
AuxIO = 0
WMAT  = 0

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do ll=1,NOccup
   do kk=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      k = kk
      l = ll

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = Occ(l)
         if(IGem(l)==1) then
            if(IGem(k)==2) then
               do i=INActive+1,NOccup
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
         
         if(IGem(l)==2) then
            if(IGem(k)==1) then
               do i=1,INActive
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO
      
      ! exchange 
      if(l<=INActive) then
         do i=1,NBasis
            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif
      if(l<=NOccup) then
         do i=1,NBasis
            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif

! AUXILIARY MATRIX WMAT
         if(l<=NOccup) then
            do ir=1,NOccup
               val = 0
               val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                    sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                    sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                    sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                    sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
               WMAT(k,ir) = WMAT(k,ir) + val
            enddo
         endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      ! Coulomb
      if(k>INActive.and.l<=NOccup) then
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then
            
            do iq=1,NOccup
               do ip=INActive+1,NBasis
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                     
                     val = 0
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                  
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     val = 0
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            enddo

         endif
      endif

      ! exchange
      if(l<=NOccup) then
         if(k>INActive) then
            ip = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then                  
                  do iq=1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = - AuxVal*val*ints(ir,iq)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               endif
            enddo
         endif
      endif

      ! T1 p->q
      if((k>INActive).and.(l<=NOccup)) then
         ip = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T3
      if((l<=NOccup).and.(k>INActive)) then
         iq = l
         ir = k
         
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                     
                  endif
               enddo
            endif
         enddo
         
      endif

      ! T4+T6: AuxInd = 2
      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
         ir = k
         iu = l

         ! T4
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
         ! T6
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then                        
                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then  
                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif

   enddo
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do ll=1,NOccup
   do kk=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      
      k = kk
      l = ll

      if(k>NOccup.or.l>NOccup) cycle

      if(present(ETot)) then  
      ! COMPUTE THE ENERGY FOR CHECKING
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))
      endif

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*Occ(k)
         if(k>INActive) then 
            do j=1,INActive
               do i=1,INActive
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         if(k<=INActive) then
            do j=INActive+1,NOccup
               do i=INActive+1,NOccup
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         do j=NOccup+1,NBasis
            do i=NOccup+1,NBasis
               HNO(i,j) = HNO(i,j) + val*ints(i,j)
            enddo
         enddo
         
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! Coulomb
      if(k==l) then
         val = 2*Occ(k)
         if(k<=INActive) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
               enddo
            enddo
         endif
         if(k>INActive.and.k<=NOccup) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
               enddo
            enddo
         endif
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      
      ! exchange
      if(l<=NOccup) then
         if(k<=NOccup) then
            iq = k
            is = l
            
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
                        
                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = - AuxVal*val*ints(ip,ir)
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               endif
            enddo
               
         endif
      endif

      ! T1 
      if((k<=NOccup).and.(l<=NOccup)) then
         iq = k
         is = l
         
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     
                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                     
                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo
            
      endif

      ! T3
      if(IGem(k)==2.and.IGem(l)==2) then
         
         do is=INActive+1,NOccup
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           
                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                RDM2val(iq,is,k,l)*ints(ip,ir)
                           
                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif

      ! T2+T5: AuxInd = 2
      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
         is = k
         iu = l
         
         ! T2
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo

         ! T5
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      endif

   enddo
enddo

close(iunit2)
deallocate(work1)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = pos(ir,is)
      
      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = pos(ip,iq)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val
         
      enddo
enddo

call clock('AB0MAT',Tcpu,Twall)

allocate(EBlock(1+NBasis-NAct))

nblk = 0

! put pack into separate subroutine
!pack AA
if(nAA>0) then

   nblk = nblk + 1
   associate(B => Eblock(nblk))
     B%n  = nAA
     B%l1 = limAA(1)
     B%l2 = limAA(2)
     allocate(B%pos(B%n))
     B%pos(1:B%n) = tmpAA(1:B%n)

     allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
     allocate(ABP(B%n,B%n),ABM(B%n,B%n))

     ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
     ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
     if(NoSt==1) then
        call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
     elseif(NoSt>1) then
        call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
     endif
 
     deallocate(ABM,ABP)
   end associate

endif

write(*,'(/,1x,a)') 'Act-Act   block diagonalized!'

allocate(Eig(nAA),Eig1(nAA),EigY(nAA,nAA),EigX(nAA,nAA),iaddr(nAA))
Eig = 0
EigY = 0
EigX=0
! nAA to liczba aktywnych
if(nAA > 0) then
   associate(B => Eblock(1))
      do i=1,B%n

           If(B%vec(i).lt.0) Then
               Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)') B%vec(i)
               B%vec(i)=0
               B%matY(1:B%n,i)=0
               B%matX(1:B%n,i)=0
           EndIf

           iaddr(i)=B%pos(i)
           EigY(i,B%l1:B%l2) = B%matY(i,1:B%n)
           EigX(i,B%l1:B%l2) = B%matX(i,1:B%n)

     enddo
     Eig(B%l1:B%l2) = B%vec(1:B%n)
   end associate

   Write(6,'(X,"Active ERPA Eigenvalues and Eigenvecs")')
   do i=1,nAA
!      Write(6,'(X,I4,E15.6)') i,eig(i)
      IStERPA(2,i)=0
      do iq=1,nAA
         if(abs(eigy(iq,i))+abs(eigx(iq,i)).gt.0.1) then
               ia=IndN(1,iaddr(iq))
               ib=IndN(2,iaddr(iq))
               ISym=MultpC(NSymNO(ia),NSymNO(ib)) 

               If(IStERPA(2,i).Eq.0) Then
                  IStERPA(2,i)=ISym
               Else
                  If(IStERPA(2,i).Ne.ISym) &
                    Write(6,'("In Y01CASDSYM_FOFO: Symm of act-act excit ",I3, &
                    " cannot be established")')i
               EndIf
          endif
       enddo
    enddo
endif

Do i=1,nAA
Ind(i)=i
EndDo

Do ISym=1,NSym
      Eig1(1:nAA)=Eig(1:nAA)
      IStart=1
      ICount=1
      Do i=1,nAA
         If(IStERPA(2,i).Eq.ISym) Then
            EigMin=Eig1(i)
            IndMin=IStart
            Do NU=IStart,nAA
                 If(IStERPA(2,NU).Eq.ISym.And.Eig1(NU).Lt.EigMin) Then
                 EigMin=Eig1(NU)
                 IndMin=NU
                 EndIf
            EndDo
            Hlp=Eig1(IStart)
            IndHlp=Ind(IStart)
            Eig1(IStart)=EigMin
            Ind(IStart)=Ind(IndMin)
            IStERPA(1,Ind(IndMin))=ICount
            Eig1(IndMin)=Hlp
            Ind(IndMin)=IndHlp
            ICount=ICount+1
         EndIf
         IStart=IStart+1
      EndDo
EndDo

!     Shift labels in each irrep depending on the number of states in SA of the energy
!     lower or equal than that of ICAS

Do I=1,NoStMx
  If(ICORR(I).Eq.0.Or.I.Eq.ICAS) Then
      ISym=IStCAS(2,I)
      Do NU=1,nAA
         If(IStERPA(2,NU).Eq.ISym) IStERPA(1,NU)=IStERPA(1,NU)+1
      EndDo
  EndIf
EndDo
NoEig=nAA
Do NU=1,nAA
   Write(6,'(X,I2,2X,I2,".",I1,E15.6)') NU,IStERPA(1,NU),IStERPA(2,NU),&
                                        Eig(NU)
EndDo

! hererxxx
If (IStCAS(1,ICAS).Eq.1.And.IStCAS(2,ICAS).Eq.1) Then
   Open(10,File='eig_1.1.dat')
   Write(10,*)nAA
   Do NU=1,nAA
      Write(10,*) IStERPA(1,NU),IStERPA(2,NU),Eig(NU)
   EndDo 
   Close(10)
EndIf

If(IStCAS(1,ICAS).Ne.1.Or.IStCAS(2,ICAS).Ne.1) Then
   INQUIRE(FILE="eig_1.1.dat",EXIST=file_exists)
   If(file_exists) Then
      Open(10,File='eig_1.1.dat')
      Read(10,*) NoEig11
      Do NU=1,NoEig11
         Read(10,*) ISt11ERPA(1,NU),ISt11ERPA(2,NU),Eig11(NU) 
         If(ISt11ERPA(1,NU).Eq.IStCAS(1,ICAS).And.ISt11ERPA(2,NU).Eq.IStCAS(2,ICAS)) Then
            AuxVal=Eig11(NU)
         EndIf
      EndDo
      Close(10)
      
      Do NU=1,NoEig11
      Write(6,'(X,"ERPA_1.1 Excit Energy ",I2,".",I1,2F22.12)') &
       ISt11ERPA(1,NU),ISt11ERPA(2,NU),Eig11(NU),Eig11(NU)-AuxVal
      Eig11(NU)=Eig11(NU)-AuxVal
      EndDo

   EndIf 
EndIf

      NegSym(1:8)=0
      Do I=1,NoStMx
      Write(6,'(X,"Excit Energy ",I1,".",I1,2F22.12)') &
      IStCAS(1,I),IStCAS(2,I),EExcit(I),EExcit(I)-EExcit(ICAS)
      If(EExcit(I)-EExcit(ICAS).Lt.0) NegSym(IStCAS(2,I))=NegSym(IStCAS(2,I))+1
      EndDo

deallocate(EigY,EigX,Eig,Eig1,iaddr)

!pack AI
do iq=1,INActive
   if(nAI(iq)>0) then

      nblk = nblk + 1
      associate(B => Eblock(nblk))
        B%n  = nAI(iq)
        B%l1 = limAI(1,iq)
        B%l2 = limAI(2,iq)
        allocate(B%pos(B%n))
        B%pos(1:B%n) = tmpAI(1:B%n,iq)

        allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
        allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  
        ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
        ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
        if(NoSt==1) then
           call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        elseif(NoSt>1) then
           call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        endif

        deallocate(ABM,ABP)
      end associate

   endif
enddo

print*
print*, 'Act-InAct block diagonalized!'

!pack AV
do ip=NOccup+1,NBasis
   if(nAV(ip)>0) then

      nblk = nblk + 1
      associate(B => Eblock(nblk))
        B%n  = nAV(ip)
        B%l1 = limAV(1,ip)
        B%l2 = limAV(2,ip)
        allocate(B%pos(B%n))
        B%pos(1:B%n) = tmpAV(1:B%n,ip)

        allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
        allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  
        ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
        ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
        if(NoSt==1) then
           call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        elseif(NoSt>1) then
           call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        endif

        deallocate(ABM,ABP)
      end associate

    endif
enddo

print*, 'Act-Virt  block diagonalized!'

!pack IV
associate(B => EblockIV)

  B%l1 = limIV(1)
  B%l2 = limIV(2)
  B%n  = B%l2-B%l1+1
  allocate(B%pos(B%n))
  B%pos(1:B%n) = tmpIV(1:B%n)

  allocate(B%vec(B%n))
  do i=1,B%n
     ii = B%l1+i-1
     B%vec(i) = ABPLUS(ii,ii)
  enddo

end associate

 call clock('PACK',Tcpu,Twall)

!! dump EBLOCKS: X(0),Y(0)
call dump_Eblock(Eblock,EblockIV,Occ,IndN,nblk,NBasis,NDimX,xy0file)

   ! AB(1) PART
   call AB_CAS_FOFO(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,&
                 IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
                 NInte1,IntJFile,IntKFile,1d0,.true.)

  call clock('AB(1)',Tcpu,Twall)
!  ! B = A.X
!  ! C = X^T.B
   allocate(workA(NDimX,NDimX))

   call ABPM_TRAN(ABPLUS,workA,EBlock,EBlockIV,nblk,NDimX,.true.)
   call clock('MULT-1',Tcpu,Twall)
   ABPLUS=workA
   call ABPM_TRAN(ABMIN,workA,EBlock,EBlockIV,nblk,NDimX,.false.)
   ABMIN=workA

   call clock('MULT-2',Tcpu,Twall)

   deallocate(workA)

   allocate(Eig(NDimX),Eig1(NDimX))
 
   ! unpack Eig (1)
   do iblk=1,nblk
      associate(B => Eblock(iblk))
 
        Eig(B%l1:B%l2) = B%vec(1:B%n)

      end associate
   enddo

   !unpack Eig (2, IV part)
   associate(B => EblockIV)
 
     do i=1,B%n
        ii = B%l1+i-1
        ipos = B%pos(i)
        Eig(ii) = B%vec(i)
     enddo

   end associate
   call clock('UNPACK',Tcpu,Twall)

deallocate(RDM2val,work2)
deallocate(ints)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

If(IStCAS(1,ICAS).Eq.1.And.IStCAS(2,ICAS).Eq.1) Then
!
! transition dipole moments in 0th and 1st order approximations
! ------------------------------------------------------------------------------

allocate(EigY(NDimX,NDimX),EigY1(NDimX,NDimX))

   EigY1 = 0
   EigY  = 0
! unpack (1)
   do iblk=1,nblk
      associate(B => Eblock(iblk))
        do i=1,B%n
           ipos = B%pos(i)
           EigY(ipos,B%l1:B%l2) = B%matY(i,1:B%n)
        enddo
      end associate
   enddo

!unpack (2, IV part)
   associate(B => EblockIV)
     do i=1,B%n
        ii = B%l1+i-1
        ipos = B%pos(i)
        EigY(ipos,ii) = 1d0/sqrt(2d0)
     enddo
   end associate

   do i=1,NDimX
      Eig1(i)=ABPLUS(i,i)+ABMIN(i,i)
   enddo

   allocate(workA(NDimX,NDimX))
   do j=1,NDimX
      if(Eig(j)/=0d0) then
         do i=1,NDimX
            if(Eig(i)/=0d0) then
               workA(i,j) = (ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
               if(Abs(Eig(i)-Eig(j))>Thresh) then
                  workA(i,j) = workA(i,j) + (ABPLUS(i,j)+ABMIN(i,j))/(Eig(j)-Eig(i))
               endif
            endif
         enddo
      endif
   enddo

   call dgemm('N','N',NDimX,NDimX,NDimX,1d0,EigY,NDimX,workA,NDimX,0d0,EigY1,NDimX)

   Do IDCORR=1,NoStMx
   If (ICORR(IDCORR).Eq.1) Then

      IERPA=0
      Do NU=1,NoEig
         If(IStERPA(1,NU).Eq.IStCAS(1,IDCORR).And.  &
         IStERPA(2,NU).Eq.IStCAS(2,IDCORR)) IERPA=NU
      EndDo

      IAC0=0
      If(IStCAS(1,ICAS).Eq.IStCAS(1,IDCORR).And. &
         IStCAS(2,ICAS).Eq.IStCAS(2,IDCORR)) IAC0=1

      If(IERPA.Ne.0) Then
            Write(6,'(/,X,"Y(0) corresponding to vector no",I2," Sym=",I1,".",I1)') & 
                          IERPA,IStERPA(1,IERPA),IStERPA(2,IERPA)
            Call TrDipMoms(IERPA,TDIP2,EigY,C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
            Write(6,'(X,"1.1->",I1,".",I1,"  Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)')IStERPA(1,IERPA), &
                           IStERPA(2,IERPA),TDIP2

            Write(6,'(/,X,"Y(0)+Y(1) corresponding to vector no",I2, " Sym=",I1,".",I1)') &
                           IERPA,IStERPA(1,IERPA),IStERPA(2,IERPA)
   
            Do I=1,NdimX
            EigY(I,IERPA)=EigY(I,IERPA)+EigY1(I,IERPA)
            EndDo
            Call TrDipMoms(IERPA,TDIP2,EigY,C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
            Write(6,'(X,"1.1->",I1,".", I1,"  Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)')IStERPA(1,IERPA), &
                           IStERPA(2,IERPA),TDIP2
      Else
            If(IAC0.Ne.1) Then
                Write(6,'(X,"1.1->",I1,".",I1,"  Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)')IStERPA(1,IERPA), &
                           IStERPA(2,IERPA),0.0
                Write(6,'(X,"1.1->",I1,".",I1,"  Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)')IStERPA(1,IERPA), &
                           IStERPA(2,IERPA),0.0
            EndIf

      EndIf

   EndIf
   EndDo

   deallocate(workA)
!
deallocate(EigY1,EigY)
! ------------------------------------------------------------------------------
! If(IStCAS(1,ICAS).Eq.1. ....
EndIf
!
! for AC0Corr
! ------------------------------------------------------------------------------

allocate(ABP(NDimX,NDimX))
ABP=0
IStateInSACAS=0
! hererxxx replace eig's with excitation energies from SA-CAS
Do I=1,NoStMx
    AuxVal=EExcit(I)-EExcit(ICAS) 
    If (AuxVal.Gt.0) Then
        Do NU=1,NoEig
             If(IStERPA(1,NU).Eq.IStCAS(1,I).And.IStERPA(2,NU).Eq.IStCAS(2,I)) Then
                     IStateInSACAS(IStERPA(1,NU),IStERPA(2,NU))=1
! hererxxx
!                    Write(6,'(X,"Replacing Eig ",I1,".",I1,F22.12," with SA-CAS ",F22.12)') &
!                    IStERPA(1,NU),IStERPA(2,NU),Eig(NU),AuxVal
!                    Eig(NU)=AuxVal
                    Write(6,'(X,"Ratio of Eig ",I1,".",I1,F22.12," to SA-CAS ",2F22.12)') &
                    IStERPA(1,NU),IStERPA(2,NU),Eig(NU),AuxVal,Eig(NU)/AuxVal
             EndIf  
        EndDo
    EndIf
EndDo

If(IStCAS(1,ICAS).Ne.1.Or.IStCAS(2,ICAS).Ne.1) Then
Do I=1,NoEig11
    AuxVal=Eig11(I)
! hererxxx
!    If (AuxVal.Gt.0) Then
        Do NU=1,NoEig
             If(IStERPA(1,NU).Eq.ISt11ERPA(1,I).And.IStERPA(2,NU).Eq.ISt11ERPA(2,I)) Then

                If(IStateInSACAS(IStERPA(1,NU),IStERPA(2,NU)).Eq.1) Then    

!                    If(Eig(NU).Gt.AuxVal) Then
                     If(Eig(NU).Gt.Abs(AuxVal)) Then
                       Write(6,'(X,"Replacing Eig ",I2,".",I1,F22.12," with ERPA(1.1) ",F22.12)') &
                       IStERPA(1,NU),IStERPA(2,NU),Eig(NU),AuxVal
                       Eig(NU)=AuxVal
                    EndIf

                EndIf

             EndIf
        EndDo
 !   EndIf
EndDo
EndIf

! set to zero egivals and eigvects of the highest value in a given symmetry if deexcitation corrections will be added
NoERPASym(1:8)=0
Do NU=1,NoEig
 If(IStERPA(1,NU).Gt.NoERPASym(IStERPA(2,NU))) NoERPASym(IStERPA(2,NU))=IStERPA(1,NU)
EndDo

allocate(IZeroNU(NDimX))
IZeroNU(1:NDimX)=0
Do I=1,8
   If(NegSym(I).Ne.0) Then
! hererxxx
!      Do J=1,NegSym(I) 
      Do J=1,NegSym(I)-1
         ISym=NoERPASym(I)-J+1
            Do NU=1,NoEig
               If(IStERPA(1,NU).Eq.ISym.And.IStERPA(2,NU).Eq.I) Then
! hererxxx
!                  Write(6,'(X,"Zeroing ERPA EigVal and EigVec in AC0",I2,".",I1,F22.12)') &
!                       IStERPA(1,NU),IStERPA(2,NU),Eig(NU)
!                  IZeroNU(NU)=1
               EndIf
            EndDo
      EndDo
   EndIf
EndDo

Do IDCORR=1,NoStMx
If (ICORR(IDCORR).Eq.1) Then

   IERPA=0
   Do NU=1,NoEig
      If(IStERPA(1,NU).Eq.IStCAS(1,IDCORR).And.  &
      IStERPA(2,NU).Eq.IStCAS(2,IDCORR)) IERPA=NU
   EndDo

   IAC0=0
   If(IStCAS(1,ICAS).Eq.IStCAS(1,IDCORR).And. &
      IStCAS(2,ICAS).Eq.IStCAS(2,IDCORR)) IAC0=1

allocate(workA(NDimX,NDimX))
workA=0

   If(IAC0.Eq.0.And.IERPA.Ne.0) Write(6,'(X, &
      "Deexcitation correction is computed for ERPA vector no", &
      I2," Sym=",I1,".",I1," Eig=",F15.8)')  IERPA, &
      IStERPA(1,IERPA),IStERPA(2,IERPA),Eig(IERPA)
   If(IAC0.Eq.0.And.IERPA.Eq.0) Write(6,'(/, &
    " ERPA vector for deexcitation correction could not &
      be determined. The correction will be set to 0.")')
   If(IAC0.Eq.1) Write(6,'(X,"AC0 correction is computed for SA-CAS state",I2," Sym=", &
      I1,".",I1)')  ICAS, IStCAS(1,ICAS),IStCAS(2,ICAS)

do j=1,NDimX
   if(Eig(j)/=0d0) then
      do i=1,NDimX
! hererxxx
!         If(IAC0.Eq.1.Or.I.Eq.IERPA.Or.J.Eq.IERPA) Then
         If((IAC0.Eq.1.And.IZeroNU(i).Eq.0.And.IZeroNU(j).Eq.0).Or.(IAC0.Eq.0.And.(I.Eq.IERPA.Or.J.Eq.IERPA))) Then
         if(Eig(i)/=0d0) workA(i,j) = 2d0*(ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
         EndIf 
      enddo
   endif
enddo
call clock('APL-AMIN/EPL',Tcpu,Twall)

call ABPM_BACKTRAN(workA,ABP,EBlock,EBlockIV,nblk,NDimX)

deallocate(workA)

call clock('ABPl-dgemm',Tcpu,Twall)
!print*, 'ABPLUS-MY',norm2(ABPLUS)

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo
!
! energy loop
allocate(work1(NBasis*NBasis),ints(NBasis,NBasis))
EAll = 0
EIntra = 0
open(newunit=iunit,file='FOFO',status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

ints = 0
kl   = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work1(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints(i,j) = work1((j-1)*NBasis+i)
           enddo
        enddo
        ints(:,NOccup+1:NBasis) = 0

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i
                Crs = C(ir)+C(is)
                Cpq = C(ip)+C(iq)

                Aux = Crs*Cpq*ABP(ipq,irs)
                EAll = EAll + Aux*ints(j,i)

                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))==1) EIntra = EIntra + Aux*ints(j,i)

              endif
           enddo
        enddo

      endif
   enddo
enddo

close(iunit)

call clock('ENE-loop',Tcpu,Twall)

ECorr = EAll-EIntra
ECorrSym(IDCORR)=ECorr

deallocate(ints,work1)
!
!
EndIf
EndDo
!
! deallocate blocks
do iblk=1,nblk
   associate(B => Eblock(iblk))

     deallocate(B%matY,B%matX,B%vec)
     deallocate(B%pos)

   end associate
enddo
! (IV part)
associate(B => EblockIV)

  deallocate(B%vec)
  deallocate(B%pos)

end associate
!
deallocate(Eig1,Eig)
deallocate(ABPLUS,ABMIN,ABP)

deallocate(IZeroNU)
end subroutine Y01CASDSYM_FOFO

end module
