module sapt_Chol_pol

use tran, only : ABPM_HALFTRAN_GEN_L
use sapt_utils
use timing

implicit none

contains

subroutine e1elst_Chol(A,B,SAPT)
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: i,j,ii,jj
integer :: NBas,NCholesky
double precision,allocatable :: Va(:,:),Vb(:,:)
double precision,allocatable :: Vabb(:,:),Vbaa(:,:)
double precision :: ea,eb,eab,elst
double precision,external  :: ddot

! set dimensions
NBas = A%NBasis
NCholesky = SAPT%NCholesky

allocate(Va(NBas,NBas),Vb(NBas,NBas),&
         Vabb(NBas,NBas),Vbaa(NBas,NBas))

call get_one_mat('V',Va,A%Monomer,NBas)
call get_one_mat('V',Vb,B%Monomer,NBas)

call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

! sum_p n_p v^B_pp
ea = 0
do i=1,A%num0+A%num1
   ea = ea + A%Occ(i)*Vbaa(i,i)
enddo
ea = 2d0*ea
!print*, 'ea',ea

! sum_q n_q v^A_qq
eb = 0
do j=1,B%num0+B%num1
   eb = eb + B%Occ(j)*Vabb(j,j)
enddo
eb = 2d0*eb
!print*, 'eb',eb

! sum_pq n_p n_q v_{pq}^{pq}
eab = 0
do j=1,B%num0+B%num1
   jj = (j-1)*(B%num0+B%num1)+j
   do i=1,A%num0+A%num1
      ii = (i-1)*(A%num0+A%num1)+i
      eab = eab + A%Occ(i)*B%Occ(j)*ddot(NCholesky,A%OO(:,ii),1,B%OO(:,jj),1)
   enddo
enddo

eab = 4d0*eab
!print*, 'eab', eab

elst = ea + eb + eab + SAPT%Vnn

call print_en('V_nn',SAPT%Vnn,.false.)
call print_en('Eelst',elst*1000,.false.)
SAPT%elst = elst

deallocate(Vb,Va,Vbaa,Vabb)

! Visualize (can be moved to a separate subroutine)
if (SAPT%Visual) then

   print*, 'VISUAL not finished for E1elst...'

!   NOccupA = A%num0+A%num1
!   NOccupB = B%num0+B%num1
!
!   allocate(QelA(NOccupA),QelB(NOccupB))
!
!   ! Eq (30) in the note
!   do i=1,NOccupA
!      QelA(i) = 2d0*A%Occ(i)*Vbaa(i,i)
!   enddo
!   ! ... and the remaining pieces from Eq (30)
!
!   ! allocate(SAPT%QelA(NOccupA,SAPT%QelB(NOccupB))
!   ! SAPT%QelA = QelA
!   ! SAPT%QelB = QelB
!
!   deallocate(QelB,QelA)

endif

end subroutine e1elst_Chol

subroutine test_Chol_ints(Flags,A,B,SAPT)
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: iunit
integer :: i,j,ij,ic,id,cd,irec
integer :: nA,nB,nC,nD,nAB,nCD
integer :: NCholesky,NBasis
double precision :: diff,diffOne,tot
double precision,allocatable :: work(:),workAO(:)

NBasis    = A%NBasis
NCholesky = SAPT%NCholesky

nA = NBasis
nB = NBasis
nC = NBasis
nD = NBasis

nAB = nA*nB
nCD = nC*nD

allocate(work(nAB),workAO(nAB))

open(newunit=iunit,file='FFFFAABB',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*NBasis**2)

tot = 0d0
irec = 0
do id=1,nD
   do ic=1,nC
      cd = ic+(id-1)*NBasis
      irec = irec + 1
      !call dgemv('T',NCholesky,nAB,1d0,B%FF,NCholesky,A%FF(1:NCholesky,cd),1,0d0,work,1)
      call dgemv('T',NCholesky,nCD,1d0,A%FF,NCholesky,B%FF(1:NCholesky,cd),1,0d0,work,1)
      read(iunit,rec=irec) workAO(1:nAB)

      ij = 0
      do j=1,nB
      do i=1,nA
         ij = ij + 1
         diffOne = abs(work(ij)) - abs(workAO(ij))
         if(abs(diffOne).gt.1d-7) print*, 'i,j',i,j,diffOne
      enddo
      enddo

      diff = norm2(work) - norm2(workAO)
      tot  = tot + abs(diff)
      write(lout,*) 'cd = ',cd,abs(diff)

   enddo
enddo
print*, 'Total : ', tot
deallocate(work,workAO)
close(iunit)

end subroutine test_Chol_ints

subroutine e2disp_Chol(Flags,A,B,SAPT)
!
! calculate 2nd order dispersion energy
! in coupled and uncoupled approximations
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: NCholesky
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
logical,allocatable          :: condOmA(:),condOmB(:)
double precision,allocatable :: OmA(:), OmB(:), &
                                OmA0(:),OmB0(:)
double precision,allocatable :: EVecA(:), EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),&
                                tmp01(:,:),tmp02(:,:)
double precision,allocatable :: work(:)
double precision :: e2d,fact,tmp
double precision :: e2du,dea,deb
double precision :: inv_omega
! for Be ERPA:
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

! Parameter(SmallE=1.D-3,BigE=1.D8)

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

! print thresholds
 if(SAPT%IPrint>1) then
    write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
    write(LOUT,'(1x,a,t18,a,e15.4)') 'SmallE','=', SmallE
    write(LOUT,'(1x,a,t18,a,e15.4)') 'BigE',  '=', BigE
 endif

! set dimensions
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA  = dimOA*dimVA
 nOVB  = dimOB*dimVB

 NCholesky = SAPT%NCholesky

! read EigValA_B
 allocate(EVecA(A%NDimX*A%NDimX),OmA(A%NDimX),  &
          EVecB(B%NDimX*B%NDimX),OmB(B%NDimX),  &
          OmA0(A%NDimX),OmB0(B%NDimX))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')

 ! uncoupled - works for CAS only
 if(Flags%ICASSCF==1) then
    allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))

    call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
    call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')
 endif

allocate(work(B%NDimX))

allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX),&
        tmp01(A%NDimX,B%NDimX),tmp02(A%NDimX,B%NDimX))

! coupled
do i=1,A%NDimX
   if(OmA(i)<0d0) write(LOUT,*) 'Negative omega A!',i,OmA(i)
enddo
do i=1,B%NDimX
   if(OmB(i)<0d0) write(LOUT,*) 'Negative omega B!',i,OmB(i)
enddo

if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then

 tmp1=0
 tmp01=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    call dgemv('T',NCholesky,B%NDimX,1d0,B%OV,NCholesky,A%OV(:,pq),1,0d0,work,1)

    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(rs)

       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + &
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)
       enddo

       associate(Y => Y01BlockA(pq))
          tmp01(Y%l1:Y%l2,rs) = tmp01(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
       end associate

    enddo
 enddo
 ! coupled
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp1,A%NDimX,EVecB,B%NDimX,0d0,tmp2,A%NDimX)

 ! uncoupled
 tmp02=0
 do rs=1,B%NDimX
    associate(Y => Y01BlockB(rs))
      call dger(A%NDimX,Y%n,1d0,tmp01(:,rs),1,Y%vec0,1,tmp02(:,Y%l1:Y%l2),A%NDimX)
    end associate
 enddo

elseif(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

 tmp1 = 0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    call dgemv('T',NCholesky,B%NDimX,1d0,B%OV,NCholesky,A%OV(:,pq),1,0d0,work,1)

    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(rs)

       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + &
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)
       enddo

    enddo
 enddo

 tmp2=0
 do j=1,B%NDimX
    do i=1,A%NDimX
       do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
       tmp2(i,j) = tmp2(i,j) + &
                    EVecB(rs+(j-1)*B%NDimX)*tmp1(i,rs)
       enddo
    enddo
 enddo

endif ! end GVB select

if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then
   ! uncoupled
    e2du = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX

          if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
             .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then


          inv_omega = 1d0/(OmA0(i)+OmB0(j))
          e2du = e2du + tmp02(i,j)**2*inv_omega

          endif
       enddo
    enddo
    SAPT%e2disp_unc = -16d0*e2du

    e2du = -16d0*e2du*1000d0

    call writeampl(tmp02,'PROP_AB0')

endif

 allocate(condOmA(A%NDimX),condOmB(B%NDimX))
 condOmA = (abs(OmA).gt.SmallE.and.abs(OmA).lt.BigE)
 condOmB = (abs(OmB).gt.SmallE.and.abs(OmB).lt.BigE)

 e2d = 0d0
 do j=1,B%NDimX
    if(condOmB(j)) then
       do i=1,A%NDimX
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
!             .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

             if(condOmA(i)) then
                e2d = e2d + tmp2(i,j)**2/(OmA(i)+OmB(j))
             endif
       enddo
    endif
 enddo
 SAPT%e2disp  = -16d0*e2d

 e2d  = -16d0*e2d*1000d0

 call print_en('E2disp',e2d,.true.)
 call print_en('E2disp(unc)',e2du,.false.)

 ! write amplitude to a file
 call writeampl(tmp2,'PROP_AB')

 !! calucate semicoupled and dexcitations
 !if(SAPT%SemiCoupled) call e2disp_semi(Flags,A,B,SAPT)

 !! calculate extrapolated E2disp
 !if(A%Cubic.or.B%Cubic) call e2disp_cpld(Flags,A,B,SAPT)

 !! calculate Wterms (deexcitations)
 !if(SAPT%Wexcit) call e2inddisp_dexc(Flags,A,B,SAPT)

 deallocate(work)

 if(Flags%ICASSCF==1) then
    ! deallocate Y01Block
    do i=1,A%NDimX
       associate(Y => Y01BlockA(i))
         deallocate(Y%vec0)
       end associate
    enddo
    do i=1,B%NDimX
       associate(Y => Y01BlockB(i))
         deallocate(Y%vec0)
       end associate
    enddo
    deallocate(Y01BlockB,Y01BlockA)
 endif

 deallocate(condOmB,condOmA)
 deallocate(tmp02,tmp01,tmp2,tmp1)
 deallocate(OmB0,OmA0,OmB,EVecB,OmA,EVecA)

end subroutine e2disp_Chol

subroutine e2disp_Cmat(Flags,A,B,SAPT)
!
! THIS IS FOR TESTING ONLY:
! calculate 2nd order dispersion energy
! use C(omega) obained brute-force (i.e.,
! by inversion: C(om)=(APLUS.AMIN + om^2)^-1.APLUS
! in coupled and uncoupled approximations
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
double precision :: fact,val
double precision :: Omega,Pi,e2d
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUSA(:,:),ABPLUSB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: work(:,:),ints(:)

double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

! Parameter(SmallE=1.D-3,BigE=1.D8)

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

! check Cholesky
if(Flags%ICholesky==1) then
   write(LOUT,'(1x,a)') 'Cholesky Cmat not ready yet! Aborting...'
   return
endif

! set dimensions
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1
dimVA = A%num1+A%num2
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

Pi = 4.0d0*atan(1.0)

! get ABPM = ABPLUS.ABMIN
allocate(ABPMA(A%NDimX,A%NDimX),ABPMB(B%NDimX,B%NDimX),&
         ABPLUSA(A%NDimX,A%NDimX),ABPLUSB(B%NDimX,B%NDimX))
allocate(work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)

deallocate(work)
allocate(work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)

deallocate(work)

! frequency integration
NFreq = 18
allocate(XFreq(NFreq),WFreq(NFreq))

call FreqGrid(XFreq,WFreq,NFreq)

allocate(CB(B%NDimX,B%NDimX))
allocate(ints(NBas**2),work(A%NDimX,B%NDimX))

e2d = 0
do ifreq=1,NFreq

   Omega = XFreq(ifreq)

   allocate(CA(A%NDimX,A%NDimX))

   call get_Cmat(CA,A%CICoef,A%IndN,ABPMA,ABPLUSA,Omega,A%NDimX,NBas)
   call get_Cmat(CB,B%CICoef,B%IndN,ABPMB,ABPLUSB,Omega,B%NDimX,NBas)

   open(newunit=iunit,file='TWOMOAB',status='OLD',&
        access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

   ints = 0
   work = 0
   do ipq=1,A%NDimX
      ip = A%IndN(1,ipq)
      iq = A%IndN(2,ipq)
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) ints(1:nOVB)

      do irs=1,B%NDimX
         ir = B%IndN(1,irs)
         is = B%IndN(2,irs)

         fact = ints(is+(ir-B%num0-1)*dimOB)

         do i=1,A%NDimX
            work(i,irs) = work(i,irs) + fact*CA(ipq,i)
         enddo

      enddo
   enddo

   deallocate(CA)
   allocate(CA(B%NDimX,B%NDimX))
   CA = 0
   ints = 0
   do ipq=1,A%NDimX
      ip = A%IndN(1,ipq)
      iq = A%IndN(2,ipq)
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) ints(1:nOVB)

      do irs=1,B%NDimX
         ir = B%IndN(1,irs)
         is = B%IndN(2,irs)

         fact = ints(is+(ir-B%num0-1)*dimOB)
         do j=1,B%NDimX
            CA(irs,j) = CA(irs,j) + fact*work(ipq,j)
         enddo

      enddo

   enddo

   val = 0
   do j=1,B%NDimX
      do i=1,B%NDimX
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo
   e2d = e2d + WFreq(ifreq)*val

   close(iunit)
   deallocate(CA)

enddo ! NFreq

SAPT%e2disp = -8d0/Pi*e2d
e2d = -8d0/Pi*e2d*1d3

call print_en('E2disp',e2d,.true.)

deallocate(ints,work)
deallocate(WFreq,XFreq)
deallocate(CB)
deallocate(ABPMB,ABPMA,ABPLUSB,ABPLUSA)

end subroutine e2disp_Cmat

subroutine e2disp_Cmat_Chol(Flags,A,B,SAPT)
!
! THIS PROCEDURE USES ADAM'S DIAG/PROJECT
! WHICH ASSUME FULL A0(NDIMX,NDIMX) MATRICES
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! with Cholesky vectors

use class_IterStats
use class_IterAlgorithm
use class_IterAlgorithmDIIS
use class_LambdaCalculator
use class_LambdaCalculatorDiag
use class_LambdaCalculatorProjector

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ij,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d
logical          :: project

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: PMatA(:,:),PmatB(:,:)
double precision,allocatable :: AB0A(:,:),AB0B(:,:)
double precision,allocatable :: ABPMA(:),ABPMB(:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:,:),CTildeB(:,:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: work(:,:),work1(:)

type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()
class(IterAlgorithmDIIS), allocatable, target    :: iterAlgo
!class(LambdaCalculatorDiag), allocatable, target :: LambdaCalcA,LambdaCalcB
class(LambdaCalculatorProjector), allocatable, target :: LambdaCalcA,LambdaCalcB

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! use projection
project=.true.
!project=.false.

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! get Pmat
print*, 'to-do: adapt Pmat procedure for SAPT...'

! monomer A
allocate(ABPMA(A%NDimX*A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,A%DChol,NCholesky,0d0,ABPTildeA,A%NDimX)
!print*, 'ABPTildeA',norm2(ABPTildeA)

deallocate(ABPLUSA,work)

! monomer B
allocate(ABPMB(B%NDimX*B%NDimX),ABPTildeB(B%NDimX,NCholesky))
allocate(ABPLUSB(B%NDimX,B%NDimX),work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,B%DChol,NCholesky,0d0,ABPTildeB,B%NDimX)
!print*, 'ABPTildeB',norm2(ABPTildeB)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner

allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

NFreq = 12
write(lout,'(/1x,a,i3)') 'SAPT%NFreq =', NFreq
write(lout,'(/1x,a)',advance='no') 'E2disp(Cmat) with projection: '
write(lout,'(1x,a,i3)') 'Adams diag version'

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(LambdaA(A%NDimX,A%NDimX),LambdaB(B%NDimX,B%NDimX))
allocate(AB0A(A%NDimX,A%NDimX),AB0B(B%NDimX,B%NDimX))

call FreqGrid(XFreq,WFreq,NFreq)

!iterAlgo = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20)
allocate(iterAlgo,SOURCE = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20))

if(project) then
   !LambdaCalcA = LambdaCalculatorProjector(A%NDimX,ABPMA,A%Pmat)
   !LambdaCalcB = LambdaCalculatorProjector(B%NDimX,ABPMB,B%Pmat)
   allocate(LambdaCalcA, SOURCE = LambdaCalculatorProjector(A%NDimX,ABPMA,A%Pmat))
   allocate(LambdaCalcB, SOURCE = LambdaCalculatorProjector(B%NDimX,ABPMB,B%Pmat))
else
 !  LambdaCalcA = LambdaCalculatorDiag(A%NDimX,ABPMA)
 !  LambdaCalcB = LambdaCalculatorDiag(B%NDimX,ABPMB)
endif

iStatsA%maxIterationsLimit = iterAlgo%maxIterations
iStatsB%maxIterationsLimit = iterAlgo%maxIterations

call LambdaCalcA%calculateInitialA()
call LambdaCalcB%calculateInitialA()

do i=1,A%NDimX
   val = ABPMA((i-1)*A%NDimX+i)
   ABPMA((i-1)*A%NDimX+i) = ABPMA((i-1)*A%NDimX+i) - val
enddo
do i=1,B%NDimX
   val = ABPMB((i-1)*B%NDimX+i)
   ABPMB((i-1)*B%NDimX+i) = ABPMB((i-1)*B%NDimX+i) - val
enddo
if(project) then
   allocate(work1(A%NDimX*A%NDimX))
   call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,A%PMat,A%NDimX,ABPMA,A%NDimX,0d0,work1,A%NDimX)
   ABPMA = ABPMA - work1
   deallocate(work1)

   allocate(work1(B%NDimX*B%NDimX))
   call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,B%PMat,B%NDimX,ABPMB,B%NDimX,0d0,work1,B%NDimX)
   ABPMB = ABPMB - work1
   deallocate(work1)
endif

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call LambdaCalcA%calculateLambda(LambdaA, OmI)
   call LambdaCalcB%calculateLambda(LambdaB, OmI)

   if(ifreq==NFreq) call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,LambdaA,A%NDimX,ABPTildeA,A%NDimX,0.0d0,CTildeA,A%NDimX)
   call iterAlgo%iterate(CTildeA, A%NDimX, NCholesky, LambdaA, ABPTildeA, ABPMA, iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,LambdaB,B%NDimX,ABPTildeB,B%NDimX,0.0d0,CTildeB,B%NDimX)
   call iterAlgo%iterate(CTildeB, B%NDimX, NCholesky, LambdaB, ABPTildeB, ABPMB, iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

enddo

write(lout,'(/1x,a)') 'C(omega): Monomer A'
call iStatsA%print()
write(lout,'(1x,a)') 'C(omega): Monomer B'
call iStatsB%print()

e2d = -8d0/Pi*e2d*1d3
!print*, 'E2disp(Cmat)',e2d
call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(ABPMB,ABPMA)
deallocate(AB0B,AB0A)
deallocate(ABPTildeB,ABPTildeA)
deallocate(LambdaB,LambdaA,XFreq,WFreq)
deallocate(CB,CA)
deallocate(CTildeB,CTildeA)

end subroutine e2disp_Cmat_Chol

subroutine e2disp_Cmat_Chol_diag(Flags,A,B,SAPT)
!
! THIS PROCEDURE USES MY DIAG PROCEDURES
! WHICH KEEP AND MULTIPLY ONLY NDIMX ELEMENTS
! PMAT PROJECTION NOT AVAILABLE!
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! with Cholesky vectors

use class_IterStats

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(SaptDIIS)    :: SAPT_DIIS

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ij,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: PMatA(:,:),PmatB(:,:)
double precision,allocatable :: ABPMA(:),ABPMB(:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:,:),CTildeB(:,:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: DiagA(:),DiagB(:)
double precision,allocatable :: LamDiaA(:),LamDiaB(:)
double precision,allocatable :: work(:,:)

type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPMA(A%NDimX*A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,A%DChol,NCholesky,0d0,ABPTildeA,A%NDimX)
print*, 'ABPTildeA',norm2(ABPTildeA)

deallocate(ABPLUSA,work)

! monomer B
allocate(ABPMB(B%NDimX*B%NDimX),ABPTildeB(B%NDimX,NCholesky))
allocate(ABPLUSB(B%NDimX,B%NDimX),work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,B%DChol,NCholesky,0d0,ABPTildeB,B%NDimX)
print*, 'ABPTildeB',norm2(ABPTildeB)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner
! making use of diagonal Lambda matrices

allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

NFreq = 12
write(lout,'(/1x,a,i3)') 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(DiagA(A%NDimX),DiagB(B%NDimX))
allocate(LamDiaA(A%NDimX),LamDiaB(B%NDimX))

call FreqGrid(XFreq,WFreq,NFreq)

iStatsA%maxIterationsLimit = SAPT_DIIS%maxIter
iStatsB%maxIterationsLimit = SAPT_DIIS%maxIter

call calculateInitialA_diag(DiagA,ABPMA,A%NDimX)
call calculateInitialA_diag(DiagB,ABPMB,B%NDimX)

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call calculateLambda_diag(LamDiaA,DiagA,OmI,A%NDimX)
   call calculateLambda_diag(LamDiaB,DiagB,OmI,B%NDimX)

   if(ifreq==NFreq) call MultpDiagMat(LamDiaA,ABPTildeA,0d0,CTildeA,A%NDimX,NCholesky)
   call Cmat_diag_iterDIIS(CTildeA,A%NDimX,NCholesky,LamDiaA,ABPTildeA,ABPMA,iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call MultpDiagMat(LamDiaB,ABPTildeB,0d0,CTildeB,B%NDimX,NCholesky)
   call Cmat_diag_iterDIIS(CTildeB,B%NDimX,NCholesky,LamDiaB,ABPTildeB,ABPMB,iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

enddo

write(lout,'(/1x,a)') 'C(omega): Monomer A'
call iStatsA%print()
write(lout,'(1x,a)') 'C(omega): Monomer B'
call iStatsB%print()

SAPT%e2disp = -8d0/Pi*e2d
e2d = -8d0/Pi*e2d*1d3
call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(ABPMB,ABPMA)
deallocate(ABPTildeB,ABPTildeA)
deallocate(DiagB,DiagA)
deallocate(LamDiaB,LamDiaA)
deallocate(XFreq,WFreq)
deallocate(CB,CA)
deallocate(CTildeB,CTildeA)

end subroutine e2disp_Cmat_Chol_diag

subroutine e2disp_Cmat_Chol_block(Flags,A,B,SAPT)
!
! THIS PROCEDURE SHOULD BE IMPROVED:
! A0 MATRICES ARE KEPT IN DIAGONAL BLOCKS
! AND BLOCK MULTIPLICATION IS USED
! (THE ABPM_HALFTRAN_LR PROCEDURE
! IN sapt_utils.f90);
! IDEALLY, THE BLOCK/DIAGONAL MULTIPLICATIONS
! COMMON FOR THE ENTIRE GAMMCOR SHOULD BE USED
! AND DIAG/BLOCK VERSION OF THE ALGORITHM SHOULD
! BE CHOSEN AUTOMATICALLY, E.G., BASED ON NACT
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! using A0 blocks (not diagonal)
! with Cholesky vectors

use class_IterStats

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n,nblkA,nblkB
integer :: maxIter
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:,:),CTildeB(:,:)
double precision,allocatable :: work(:,:)
double precision,allocatable :: WorkA(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)
type(EBlockData)             :: LambdaIVA,LambdaIVB
type(EBlockData),allocatable :: LambdaA(:),LambdaB(:)

type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPMA(A%NDimX,A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,A%DChol,NCholesky,0d0,ABPTildeA,A%NDimX)

deallocate(ABPLUSA,work)

! monomer B
allocate(ABPMB(B%NDimX,B%NDimX),ABPTildeB(B%NDimX,NCholesky))
allocate(ABPLUSB(B%NDimX,B%NDimX),work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,B%DChol,NCholesky,0d0,ABPTildeB,B%NDimX)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner
NFreq = 12
write(lout,'(/1x,a,i3)') 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! get A2 = ABPM - ABPM0
call calculateInitialA_blk(ABPMA,A0BlkA,A0BlkIVA,nblkA,A%NDimX,'A0BLK_A')
call calculateInitialA_blk(ABPMB,A0BlkB,A0BlkIVB,nblkB,B%NDimX,'A0BLK_B')

maxIter = 20
iStatsA%maxIterationsLimit = maxIter
iStatsB%maxIterationsLimit = maxIter

allocate(LambdaA(nblkA),LambdaB(nblkB))

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call calculateLambda_blk(LambdaA,LambdaIVA,OmI**2,A%NDimX,nblkA,A0BlkA,A0BlkIVA)
   call calculateLambda_blk(LambdaB,LambdaIVB,OmI**2,B%NDimX,nblkB,A0BlkB,A0BlkIVB)

   if(ifreq==NFreq) call ABPM_HALFTRAN_GEN_L(ABPTildeA,CTildeA,0d0,LambdaA,LambdaIVA,nblkA,A%NDimX,NCholesky,'X')
   call Cmat_blk_iterDIIS(CTildeA,A%NDimX,NCholesky,nblkA,LambdaA,LambdaIVA,ABPTildeA,ABPMA,iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call ABPM_HALFTRAN_GEN_L(ABPTildeB,CTildeB,0d0,LambdaB,LambdaIVB,nblkB,B%NDimX,NCholesky,'X')
   call Cmat_blk_iterDIIS(CTildeB,B%NDimX,NCholesky,nblkB,LambdaB,LambdaIVB,ABPTildeB,ABPMB,iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

   call release_ac0block(LambdaA,LambdaIVA,nblkA)
   call release_ac0block(LambdaB,LambdaIVB,nblkB)

enddo

write(lout,'(/1x,a)') 'C(omega): Monomer A'
call iStatsA%print()
write(lout,'(1x,a)') 'C(omega): Monomer B'
call iStatsB%print()

SAPT%e2disp  = -8d0/Pi*e2d

e2d = -8d0/Pi*e2d*1d3
call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(CB,CA)
deallocate(CTildeB,CTildeA)
deallocate(WFreq,XFreq)
!deallocate(DCholB,DCholA)
deallocate(ABPMB,ABPMA)
deallocate(ABPTildeB,ABPTildeA)

end subroutine e2disp_Cmat_Chol_block

subroutine e2disp_Cmat_Chol_proj(Flags,A,B,SAPT)
!
! THIS PROCEDURE USES MY PROJECT PROCEDURES
! BUT THERE IS PERHAPS NO SENSE OF DOING THAT?
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! with Cholesky vectors

use class_IterStats
use class_IterAlgorithm
use class_IterAlgorithmDIIS

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(SaptDIIS)    :: SAPT_DIIS

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ij,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d
logical :: diag, blk

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: PMatA(:,:),PmatB(:,:)
double precision,allocatable :: ABPMA(:),ABPMB(:),&
                                AB0A(:,:),AB0B(:,:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:,:),CTildeB(:,:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: work(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()
class(IterAlgorithmDIIS), allocatable, target    :: iterAlgo

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! diagonal or block
diag = .true.
blk  = .false.
!diag = .false.
!blk  = .true.

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPMA(A%NDimX*A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,A%DChol,NCholesky,0d0,ABPTildeA,A%NDimX)
print*, 'ABPTildeA',norm2(ABPTildeA)

deallocate(ABPLUSA,work)

! monomer B
allocate(ABPMB(B%NDimX*B%NDimX),ABPTildeB(B%NDimX,NCholesky))
allocate(ABPLUSB(B%NDimX,B%NDimX),work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,B%DChol,NCholesky,0d0,ABPTildeB,B%NDimX)
print*, 'ABPTildeB',norm2(ABPTildeB)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner
! making use of diagonal Lambda matrices

allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

NFreq = 12
write(lout,'(/1x,a)',advance='no') 'E2disp(Cmat) with projection: '
if(diag) write(lout,'(1x,a,i3)') 'diagonal version'
if(blk)  write(lout,'(1x,a,i3)') 'block version'

write(lout,'(/1x,a,i3)') 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(AB0A(A%NDimX,A%NDimX),AB0B(B%NDimX,B%NDimX))
allocate(LambdaA(A%NDimX,A%NDimX),LambdaB(B%NDimX,B%NDimX))

call FreqGrid(XFreq,WFreq,NFreq)

!iterAlgo = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20)
allocate(iterAlgo, SOURCE = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20))

iStatsA%maxIterationsLimit = iterAlgo%maxIterations
iStatsB%maxIterationsLimit = iterAlgo%maxIterations

if(diag) then
  call calculateInitialA_diagP(AB0A,ABPMA,A%PMat,A%NDimX)
  call calculateInitialA_diagP(AB0B,ABPMB,B%PMat,B%NDimX)
elseif(blk) then
  call calculateInitialA_blkP(ABPMA,AB0A,A%PMat,A%NDimX,'A0BLK_A')
  call calculateInitialA_blkP(ABPMB,AB0B,B%PMat,B%NDimX,'A0BLK_B')
endif

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call calculateLambda_Pmat(LambdaA,AB0A,OmI,A%NDimX)
   call calculateLambda_Pmat(LambdaB,AB0B,OmI,B%NDimX)

   if(ifreq==NFreq) call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,LambdaA,A%NDimX,ABPTildeA,A%NDimX,0.0d0,CTildeA,A%NDimX)
   call iterAlgo%iterate(CTildeA, A%NDimX, NCholesky, LambdaA, ABPTildeA, ABPMA, iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,LambdaB,B%NDimX,ABPTildeB,B%NDimX,0.0d0,CTildeB,B%NDimX)
   call iterAlgo%iterate(CTildeB, B%NDimX, NCholesky, LambdaB, ABPTildeB, ABPMB, iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

enddo

write(lout,'(/1x,a)') 'C(omega): Monomer A'
call iStatsA%print()
write(lout,'(1x,a)') 'C(omega): Monomer B'
call iStatsB%print()

e2d = -8d0/Pi*e2d*1d3
call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(AB0B,AB0A)
deallocate(ABPMB,ABPMA)
deallocate(ABPTildeB,ABPTildeA)
deallocate(LambdaB,LambdaA)
deallocate(XFreq,WFreq)
deallocate(CB,CA)
deallocate(CTildeB,CTildeA)

end subroutine e2disp_Cmat_Chol_proj

subroutine e2ind_CAlphaTilde_block(Flags,A,B,SAPT)
!
! calculate 2nd order induction energy
! using expansion of C(0) in alpha around alpha=0, up to Max_Cn order
! use A0 blocks (not diagonal)
! with Cholesky vectors
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NCholesky,NBas
integer :: i,j,ip,iq,ipq,ir,is,irs
integer :: nblkA,nblkB
double precision :: e2ba, e2ab
double precision,allocatable :: WaBB(:,:),WbAA(:,:),   &
                                WaChBB(:),WbChAA(:)
double precision,allocatable :: A1A(:,:),A1B(:,:),&
                                A2A(:,:),A2B(:,:),&
                                ABP0TildeA(:,:),ABP0TildeB(:,:), &
                                ABP1TildeA(:,:),ABP1TildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:), &
                                CTildeA(:,:),CTildeB(:,:),   &
                                CA(:,:),CB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

!! set dimensions
!NCholesky = SAPT%NCholesky
!
!! check MCBS/DCBS
!if(A%NBasis.ne.B%NBasis) then
!   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
!   stop
!else
!   NBas = A%NBasis
!endif
!
!allocate(WaBB(NBas,NBas),WbAA(NBas,NBas))
!
!call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)
!call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)
!
!allocate(A1A(A%NDimX,A%NDimX),A2A(A%NDimX,A%NDimX))
!allocate(ABP0TildeA(A%NDimX,NCholesky), &
!         ABP1TildeA(A%NDimX,NCholesky), &
!         DCholA(A%NDimX,NCholesky))
!
!call prepare_resp_Cmat(A,A1A,A2A,ABP0TildeA,ABP1TIldeA,DCholA,A%NDimX,NCholesky,NBas)
!
!allocate(WbChAA(NCholesky))
!do i=1,NCholesky
!   do ipq=1,A%NDimX
!      ip = A%IndN(1,ipq)
!      iq = A%IndN(2,ipq)
!   
!      WbChAA(i) = WbChAA(i) + WbAA(ip,iq)*DCholA(i,ipq)
!
!   enddo
!enddo
!
!allocate(A1B(B%NDimX,B%NDimX),A2B(B%NDimX,B%NDimX))
!allocate(ABP0TildeB(B%NDimX,NCholesky), &
!         ABP1TildeB(B%NDimX,NCholesky), &
!         DCholB(B%NDimX,NCholesky))
!
!call prepare_resp_Cmat(B,A1B,A2B,ABP0TildeB,ABP1TIldeB,DCholB,B%NDimX,NCholesky,NBas)
!
!allocate(CTildeA(A%NDimX,NCholesky),C0TildeA(A%NDimX,NCholesky))
!allocate(CA(A%NDimX,A%NDimX))
!call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
!
!call C_AlphaExpand(CTildeA,C0TildeA,0d0,SAPT%Max_Cn, &
!                   A1A,A2A,ABP0TildeA,ABP1TildeA,    &
!                   A0BlkA,A0BlkIVA,nblkA,NCholesky,A%NDimX)
!
!print*, 'CTildeA',norm2(CTildeA)
!!call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,DCholA,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
!call dgemm('N','N',A%NDimX,A%NDimX,NCholesky,1d0,CTildeA,A%NDimX,DCholA,NCholesky,0d0,CA,A%NDimX)
!
!e2ba = 0
!do ipq=1,A%NDimX
!   ip = A%IndN(1,ipq)
!   iq = A%IndN(2,ipq)
!   do irs=1,A%NDimX
!      ir = A%IndN(1,irs)
!      is = A%IndN(2,irs)
!      e2ba = e2ba + WbAA(ip,iq)*CA(ipq,irs)*WbAA(ir,is)
!   enddo
!enddo
!e2ba = -0.5d0*e2ba
!print*, 'e2ba',e2ba*1000
!
!!call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')
!!call C_AlphaExpand(CTildeB,C0TildeB,OmI,SAPT%Max_Cn,A1B,A2B,ABP0TildeB,ABP1TildeB, &
!!                   A0BlkB,A0BlkIVB,nblkB,NCholesky,B%NDimX)
!
!deallocate(DCholB,DCholA)
!deallocate(WbAA,WaBB)
!deallocate(CA)
!deallocate(C0TildeA,CTildeA)
!deallocate(A2B,A2A,A1B,A1A)
!deallocate(ABP0TildeB,ABP0TildeA,ABP1TildeB,ABP1TIldeA)

end subroutine e2ind_CAlphaTilde_block

subroutine e2disp_CAlphaTilde_block(Flags,A,B,SAPT)
!
! calculate 2nd order dispersion energy
! using expansion of C(w) in alpha around alpha=0, up to Max_Cn order
! https://doi.org/10.1021/acs.jpclett.3c01568
! use A0 blocks (not diagonal)
! with Cholesky vectors
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: nblkA,nblkB
integer :: N,Max_Cn
double precision :: fact,val,val2
double precision :: Cpq,Crs
double precision :: XFactorial,XN1,XN2
double precision :: ACAlpha,OmI,Pi,e2du,e2d
logical :: both

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUS1A(:,:),ABPLUS1B(:,:),&
                                ABMIN1A(:,:),ABMIN1B(:,:),  &
                                A1A(:,:),A1B(:,:),&
                                A2A(:,:),A2B(:,:),&
                                ABP0TildeA(:,:),ABP0TildeB(:,:),&
                                ABP1TildeA(:,:),ABP1TildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:), &
                                C1TildeA(:,:),C1TildeB(:,:), &
                                C2TildeA(:,:),C2TildeB(:,:), &
                                CTildeA(:,:), CTildeB(:,:)
double precision,allocatable :: WorkA(:,:),WorkB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

type(EBlockData)             :: LambdaIVA,LambdaIVB
type(EBlockData),allocatable :: LambdaA(:),LambdaB(:)

! for Visualize
integer          :: dimOA, dimOB
double precision :: e2dv
integer,allocatable          :: posA(:,:),posB(:,:)
double precision,allocatable :: Qmat(:,:),CAB(:,:),DAB(:,:)
! for Visualize with local
character*60 line
integer :: iip,ILoc
double precision,allocatable :: ALOC(:,:),BLOC(:,:)
double precision,allocatable :: ALOAO(:,:),BLOAO(:,:)
double precision,allocatable :: QAux1(:,:,:),QAuxC(:,:,:,:),QAuxD(:,:,:,:)

! test
double precision :: ErrMax
double precision :: Tcpu,Twall

print*, ''
print*, 'Experimental E2disp(CAlpha) procedure...'

! both = coupled + uncoupled
both = SAPT%iCpld

! timing
call clock('START',Tcpu,Twall)

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

! for Visualize
ILoc=0
if (SAPT%Visual) then

   dimOA = A%num0 + A%num1
   allocate(posA(NBas,NBas))
   posA = 0
   do i=1,A%NDimX
      posA(A%IndN(1,i),A%IndN(2,i)) = A%IndX(i)
   enddo

   dimOB = B%num0 + B%num1
   allocate(posB(NBas,NBas))
   posB = 0
   do i=1,B%NDimX
      posB(B%IndN(1,i),B%IndN(2,i)) = B%IndX(i)
   enddo
      allocate(ALOC(A%NBasis,A%NBasis),BLOC(B%NBasis,B%NBasis), &
              ALOAO(A%NBasis,A%NBasis),BLOAO(B%NBasis,B%NBasis))
      ALOAO=A%CMO
      BLOAO=B%CMO 
! for Q_AB in localized orbitals
  if(ILoc==1) then
      allocate(CA(A%NBasis,A%NBasis),CB(B%NBasis,B%NBasis))
      allocate(A1A(A%NBasis,A%NBasis))
      allocate(A1B(A%NBasis,A%NBasis))

      open(10,file="loc_a.dat",status='OLD')
      read(10,'(A10)')line
      read(10,'(A10)')line
      do i=1,A%NBasis
          read(10,*) (CA(j,i),j=1,A%NBasis)
      enddo
      close(10)
      open(10,file="loc_b.dat",status='OLD')
      read(10,'(A10)')line
      read(10,'(A10)')line
      do i=1,B%NBasis
          read(10,*) (CB(j,i),j=1,B%NBasis)
      enddo
      close(10)

      ALOAO=CA
      BLOAO=CB

      call CpyM(A1A,A%CMO,A%NBasis)
      A1A=transpose(A1A)
      call minvr(A1A,1.0d-7,val,i,A%NBasis)
      if (i.ne.0) then
         Write(6,'(/,'' ERROR : transformation from ao to mo basis is singular'')')
         Stop
      endif
      call CpyM(A1B,B%CMO,A%NBasis)
      A1B=transpose(A1B)
      call minvr(A1B,1.0d-7,val,i,B%NBasis)
      if (i.ne.0) then
         Write(6,'(/,'' ERROR : transformation from ao to mo basis is singular'')')
         Stop
      endif
      call dgemm('N','N',A%NBasis,A%NBasis,A%NBasis,1.d0,CA,A%NBasis,A1A,A%NBasis,0.0d0,ALOC,A%NBasis)
      call dgemm('N','N',B%NBasis,B%NBasis,B%NBasis,1.d0,CB,B%NBasis,A1B,B%NBasis,0.0d0,BLOC,B%NBasis)
      deallocate(CA,CB,A1A,A1B)
  endif 

endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPLUS1A(A%NDimX,A%NDimX),ABMIN1A(A%NDimX,A%NDimX))
allocate(A1A(A%NDimX,A%NDimX),A2A(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1A
read(iunit) ABMIN1A

close(iunit)

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkA,A0BlkIVA,A%IndN,A%CICoef,nblkA,NBas,A%NDimX,'XY0_A')

! AB1 = AB1 - A0
call add_blk_right(ABPLUS1A,A0BlkA,A0BlkIVA,-1d0,.false.,nblkA,A%NDimX)
call add_blk_right(ABMIN1A, A0BlkA,A0BlkIVA,-1d0,.true., nblkA,A%NDimX)

print*, 'ABPLUS1A',norm2(ABPLUS1A)
print*, 'ABMIN1A ',norm2(ABMIN1A)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call ABPM_HALFTRAN_GEN_L(ABMIN1A, A1A,0.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,A%NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1A,A1A,1.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,A%NDimX,'X')
print*, 'A1A',norm2(A1A)

!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS1A,A%NDimX,ABMIN1A,A%NDimX,0.0d0,A2A,A%NDimX)
deallocate(ABMIN1A)
print*, 'A2A',norm2(A2A)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeA(A%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(A%DChol),ABP0TildeA,0.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde',norm2(ABP0TildeA)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeA(A%NDimX,NCholesky))
Call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUS1A,A%NDimX,A%DChol,NCholesky,0.0d0,ABP1TildeA,A%NDimX)
print*, 'APLUS1Tilde',norm2(ABP1TildeA)

deallocate(ABPLUS1A)
call release_ac0block(A0BlkA,A0BlkIVA,nblkA)
deallocate(A0BlkA)

! monomer B
allocate(ABPLUS1B(B%NDimX,B%NDimX),ABMIN1B(B%NDimX,B%NDimX))
allocate(A1B(B%NDimX,B%NDimX),A2B(B%NDimX,B%NDimX))

!allocate(ABPMB(B%NDimX,B%NDimX),ABPTildeB(B%NDimX,NCholesky))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1B
read(iunit) ABMIN1B

close(iunit)

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkB,A0BlkIVB,B%IndN,B%CICoef,nblkB,NBas,B%NDimX,'XY0_B')

! AB1 = AB1 - A0
call add_blk_right(ABPLUS1B,A0BlkB,A0BlkIVB,-1d0,.false.,nblkB,B%NDimX)
call add_blk_right(ABMIN1B, A0BlkB,A0BlkIVB,-1d0,.true., nblkB,B%NDimX)

print*, 'ABPLUS1B',norm2(ABPLUS1B)
print*, 'ABMIN1B ',norm2(ABMIN1B)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call ABPM_HALFTRAN_GEN_L(ABMIN1B, A1B,0.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,B%NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1B,A1B,1.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,B%NDimX,'X')
print*, 'A1B',norm2(A1B)

!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS1B,B%NDimX,ABMIN1B,B%NDimX,0.0d0,A2B,B%NDimX)
deallocate(ABMIN1B)
print*, 'A2B',norm2(A2B)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeB(B%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(B%DChol),ABP0TildeB,0.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde-b',norm2(ABP0TildeB)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeB(B%NDimX,NCholesky))
Call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUS1B,B%NDimX,B%DChol,NCholesky,0.0d0,ABP1TildeB,B%NDimX)
print*, 'APLUS1Tilde-b',norm2(ABP1TildeB)

deallocate(ABPLUS1B)
call release_ac0block(A0BlkB,A0BlkIVB,nblkB)
deallocate(A0BlkB)

! get CAlphaTilde in an iterative manner

NFreq = 12
Max_Cn = SAPT%Max_Cn
print*, ''
print*, 'CAlphaTilde: '
print*, 'SAPT%NFreq =', NFreq
print*, 'SAPT%MaxCn =', Max_Cn

allocate(XFreq(NFreq),WFreq(NFreq))

allocate(C0TildeA(A%NDimX,NCholesky),C0TildeB(B%NDimX,NCholesky))
allocate(C1TildeA(A%NDimX,NCholesky),C1TildeB(B%NDimX,NCholesky))
allocate(C2TildeA(A%NDimX,NCholesky),C2TildeB(B%NDimX,NCholesky))
allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(WorkA(A%NDimX,NCholesky),WorkB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

if (SAPT%Visual.and.ILoc.eq.0) allocate(Qmat(dimOA,dimOB),CAB(A%NDimX,B%NDimX),DAB(A%NDimX,B%NDimX))
if (SAPT%Visual.and.ILoc.eq.1) allocate(Qmat(dimOA,dimOB))
if (SAPT%Visual) Qmat = 0d0

call FreqGrid(XFreq,WFreq,NFreq)

! read ABPLUS0.ABMIN0 blocks
call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')

e2d  = 0
e2du = 0
e2dv = 0d0
ErrMax = 0d0
do ifreq=1,NFreq

   OmI = XFreq(ifreq)
   print*,'OmI', OmI

   !if(both) then

      ! coupled
      print*, 'ErrMax', ErrMax
      call C_AlphaExpand(CTildeA,C0TildeA,OmI,Max_Cn,A1A,A2A,ABP0TildeA,ABP1TildeA, &
                         A0BlkA,A0BlkIVA,nblkA,NCholesky,A%NDimX,ErrMax,ifreq)
      call C_AlphaExpand(CTildeB,C0TildeB,OmI,Max_Cn,A1B,A2B,ABP0TildeB,ABP1TildeB, &
                         A0BlkB,A0BlkIVB,nblkB,NCholesky,B%NDimX,ErrMax,ifreq)

      call dgemm('T','T',NCholesky,NCholesky,A%NDimX,1d0,CTildeA,A%NDimX,A%DChol,NCholesky,0d0,CA,NCholesky)
      call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

      val = 0
      do j=1,NCholesky
         do i=1,NCholesky
            val = val + CA(i,j)*CB(i,j)
         enddo
      enddo

      e2d = e2d + WFreq(ifreq)*val

!  for Q_AB in localized orbitals
   if(SAPT%Visual.and.ILoc.eq.1) then

      allocate(CAB(A%NDimX,B%NDimX))
      call dgemm('N','T',A%NDimX,B%NDimX,NCholesky,1d0,CTildeA,A%NDimX,CTildeB,B%NDimX,0d0,CAB,A%NDimX)

      allocate(QAux1(dimOA,NBas,B%NDimX))
      do i=1,dimOA
       do is=1,dimOB
         do ir=1,NBas
            if(posB(ir,is)/=0) then
               irs = posB(ir,is)
               do ip=1,NBas
                 QAux1(i,ip,irs)=0.0
                  do iq=1,dimOA
                     if(posA(ip,iq)/=0) then
                        ipq = posA(ip,iq)
                        QAux1(i,ip,irs)= QAux1(i,ip,irs)+ ALOC(i,iq)*CAB(ipq,irs)
                     endif
                  enddo
               enddo
            endif
         enddo
       enddo
      enddo

      deallocate(CAB)
      allocate(QAuxC(dimOA,NBas,dimOB,NBas))
     
      do j=1,dimOB
       do i=1,dimOA
         do ip=1,NBas
               do ir=1,NBas
                 QAuxC(i,ip,j,ir)=0.0
                  do is=1,dimOB
                     if(posB(ir,is)/=0) then
                        irs = posB(ir,is)
                        QAuxC(i,ip,j,ir)= QAuxC(i,ip,j,ir)+ BLOC(j,is)*QAux1(i,ip,irs)
                     endif
                  enddo
               enddo
         enddo
       enddo
      enddo
    
      allocate(DAB(A%NDimX,B%NDimX))
      call dgemm('T','N',A%NDimX,B%NDimX,NCholesky,1d0,A%DChol,NCholesky,B%DChol,NCholesky,0d0,DAB,A%NDimX)
 
      do i=1,dimOA
       do is=1,dimOB
         do ir=1,NBas
            if(posB(ir,is)/=0) then
               irs = posB(ir,is)
               do ip=1,NBas
                 QAux1(i,ip,irs)=0.0
                  do iq=1,dimOA
                     if(posA(ip,iq)/=0) then
                        ipq = posA(ip,iq)
                        QAux1(i,ip,irs)= QAux1(i,ip,irs)+ ALOC(i,iq)*DAB(ipq,irs)
                     endif
                  enddo
               enddo
            endif
         enddo
       enddo
      enddo 
      deallocate(DAB)

      allocate(QAuxD(dimOA,NBas,dimOB,NBas))

      do j=1,dimOB
       do i=1,dimOA
         do ip=1,NBas
               do ir=1,NBas
                 QAuxD(i,ip,j,ir)=0.0
                  do is=1,dimOB
                     if(posB(ir,is)/=0) then
                        irs = posB(ir,is)
                        QAuxD(i,ip,j,ir)= QAuxD(i,ip,j,ir)+ BLOC(j,is)*QAux1(i,ip,irs)
                     endif
                  enddo
               enddo
         enddo
       enddo
      enddo
      deallocate(QAux1)

      do j=1,dimOB
         do ir=1,NBas
               do i=1,dimOA
                  do ip=1,NBas
                        Qmat(i,j) = Qmat(i,j) + WFreq(ifreq)*QAuxC(i,ip,j,ir)*QAuxD(i,ip,j,ir)
                  enddo
               enddo
         enddo
      enddo
     deallocate(QAuxC,QAuxD)

   endif

   if(SAPT%Visual.and.ILoc.eq.0) then

      call dgemm('T','N',A%NDimX,B%NDimX,NCholesky,1d0,A%DChol,NCholesky,B%DChol,NCholesky,0d0,DAB,A%NDimX)
      call dgemm('N','T',A%NDimX,B%NDimX,NCholesky,1d0,CTildeA,A%NDimX,CTildeB,B%NDimX,0d0,CAB,A%NDimX)
      do is=1,dimOB
         do ir=1,NBas
            if(posB(ir,is)/=0) then
               irs = posB(ir,is)
               do iq=1,dimOA
                  do ip=1,NBas
                     if(posA(ip,iq)/=0) then
                        ipq = posA(ip,iq)
                        Qmat(iq,is) = Qmat(iq,is) + WFreq(ifreq)*CAB(ipq,irs)*DAB(ipq,irs)
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo

   endif

   !endif

   !! uncoupled

   !call C_AlphaExpand_unc(C0TildeA,OmI,A1A,A2A,ABP0TildeA,ABP1TildeA, &
   !                   A0BlkA,A0BlkIVA,nblkA,NCholesky,A%NDimX)
   !call C_AlphaExpand_unc(C0TildeB,OmI,A1B,A2B,ABP0TildeB,ABP1TildeB, &
   !                   A0BlkB,A0BlkIVB,nblkB,NCholesky,B%NDimX)

   !call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,DCholA,NCholesky,C0TildeA,A%NDimX,0d0,CA,NCholesky)
   !call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,DCholB,NCholesky,C0TildeB,B%NDimX,0d0,CB,NCholesky)

   !val2 = 0
   !do j=1,NCholesky
   !   do i=1,NCholesky
   !      val2 = val2+ CA(j,i)*CB(i,j)
   !   enddo
   !enddo
   !e2du = e2du + WFreq(ifreq)*val2

enddo

!SAPT%e2disp_unc = -32d0/Pi*e2du
!e2du = -32d0/Pi*e2du*1d3
!call print_en('E2disp(Alph,unc)',e2du,.false.)
!if(both) then
   SAPT%e2disp  = -32d0/Pi*e2d

   ! skip uncoupled E2disp/E2exdisp with Cholesky...
   SAPT%e2disp_unc = 0d0 !SAPT%e2disp

   e2d = -32d0/Pi*e2d*1d3
   !print*, 'e2d = ',e2d
   call print_en('E2disp(CAlpha)',e2d,.false.)
!endif

! test Visualize
if (SAPT%Visual) then

   Qmat = -32d0/Pi * Qmat
   e2dv = 0
   do ir=1,dimOB
      do ip=1,dimOA
         !print*, 'Qmat(p,r)',ip,ir,Qmat(ip,ir)
         e2dv = e2dv + Qmat(ip,ir)
      enddo
   enddo

   !e2dv = -32d0/Pi*e2dv*1d3
   !call print_en('E2disp(Visual)',e2dv,.false.)
   e2dv = e2dv*1d3
   call print_en('E2disp(Visual)',e2dv,.false.)
   if(abs(e2d-e2dv).gt.1.d-8) write(6,*)'!!! warning: wrong disp_en from Q_AB !!!',e2dv-e2d

   allocate(SAPT%Qmat(NBas,NBas),SAPT%ALOC(NBas,NBas),SAPT%BLOC(NBas,NBas))
   SAPT%Qmat = Qmat
   SAPT%ALOC = transpose(ALOAO)
   SAPT%BLOC = transpose(BLOAO)
   
endif

call clock('E2disp(CAlpha)',Tcpu,Twall)

deallocate(WFreq,XFreq)
deallocate(A2A,A1A)
deallocate(ABP1TildeA,ABP0TildeA)
deallocate(CTildeB,CTildeA)
deallocate(C2TildeB,C2TildeA)
deallocate(C1TildeB,C1TildeA)
deallocate(C0TildeB,C0TildeA)
deallocate(WorkB,WorkA)
deallocate(CB,CA)
!deallocate(B%DCholB,A%DCholA)

end subroutine e2disp_CAlphaTilde_block

subroutine e2disp_CAlphaTilde_unc(Flags,A,B,SAPT)
!
! calculate 2nd order uncoupled dispersion energy
! using expansion of C(w) in alpha around alpha=0, up to Max_Cn order
! use A0 blocks (not diagonal)
! with Cholesky vectors
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBasis,NCholesky
integer :: nblkA,nblkB
integer :: i,j
integer :: ifreq,NFreq

double precision :: Pi
double precision :: OmI,val
double precision :: e2du

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABP0TildeA(:,:),ABP0TildeB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

Pi = 4.0d0*atan(1.0)

! set dimensions
NBasis = A%NBasis
NCholesky = SAPT%NCholesky

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkA,A0BlkIVA,A%IndN,A%CICoef,nblkA,NBasis,A%NDimX,'XY0_A')

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeA(A%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(A%DChol),ABP0TildeA,0.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde-a',norm2(ABP0TildeA)

call release_ac0block(A0BlkA,A0BlkIVA,nblkA)
deallocate(A0BlkA)

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkB,A0BlkIVB,B%IndN,B%CICoef,nblkB,NBasis,B%NDimX,'XY0_B')

allocate(ABP0TildeB(B%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(B%DChol),ABP0TildeB,0.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde-b',norm2(ABP0TildeB)

call release_ac0block(A0BlkB,A0BlkIVB,nblkB)
deallocate(A0BlkB)

! get CAlphaTilde in 0th order
NFreq = 12
print*, ''
print*, 'E2disp, unc from CAlphaTilde(0)'
print*, 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(C0TildeA(A%NDimX,NCholesky),C0TildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! read ABPLUS0.ABMIN0 blocks
call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')

e2du = 0
do ifreq=1,NFreq

   OmI = XFreq(ifreq)

   call C_AlphaExpand_unc(C0TildeA,OmI,ABP0TildeA, &
                      A0BlkA,A0BlkIVA,nblkA,NCholesky,A%NDimX)
   call C_AlphaExpand_unc(C0TildeB,OmI,ABP0TildeB, &
                      A0BlkB,A0BlkIVB,nblkB,NCholesky,B%NDimX)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,C0TildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,C0TildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo
   e2du = e2du + WFreq(ifreq)*val
   !print*,'OmI, E2d,unc=', OmI,-32d0/Pi*e2du*1d3
   write(6, '(1x,"OmI= ",F12.6,", E2d,unc=",F12.6)') OmI,-32d0/Pi*e2du*1d3

enddo

SAPT%e2disp_unc = -32d0/Pi*e2du
e2du = -32d0/Pi*e2du*1d3

!print*, 'E2disp(unc) =',e2du
call print_en('E2disp,unc,C(0)',e2du,.true.)

deallocate(CB,CA)
deallocate(C0TildeB,C0TildeA)
deallocate(ABP0TildeB,ABP0TildeA)
deallocate(WFreq,XFreq)

end subroutine e2disp_CAlphaTilde_unc

subroutine e2disp_CAlphaTilde_full(Flags,A,B,SAPT)
!
! this is CAlphaTilde procedure which uses
! NDimX*NDimX Lambda (no block structure)
! THIS IS HERE FOR TESTING ONLY
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: nblkA,nblkB
integer :: N,Max_Cn
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: XFactorial,XN1,XN2
double precision :: ACAlpha,OmI,Pi,e2d

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUS0A(:,:),ABPLUS0B(:,:),&
                                ABMIN0A(:,:),ABMIN0B(:,:),  &
                                ABPLUS1A(:,:),ABPLUS1B(:,:),&
                                ABMIN1A(:,:),ABMIN1B(:,:),  &
                                A1A(:,:),A1B(:,:),&
                                A2A(:,:),A2B(:,:),&
                                ABP0TildeA(:,:),ABP0TildeB(:,:),&
                                ABP1TildeA(:,:),ABP1TildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:), &
                                C1TildeA(:,:),C1TildeB(:,:), &
                                C2TildeA(:,:),C2TildeB(:,:), &
                                CTildeA(:,:), CTildeB(:,:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: WorkA(:,:),WorkB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

!type(EBlockData)             :: LambdaIVA,LambdaIVB
!type(EBlockData),allocatable :: LambdaA(:),LambdaB(:)


ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPLUS0A(A%NDimX,A%NDimX),ABMIN0A(A%NDimX,A%NDimX))
allocate(ABPLUS1A(A%NDimX,A%NDimX),ABMIN1A(A%NDimX,A%NDimX))
allocate(A1A(A%NDimX,A%NDimX),A2A(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1A
read(iunit) ABMIN1A

close(iunit)

open(newunit=iunit,file='A0MAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS0A
read(iunit) ABMIN0A

close(iunit)

! AB1 = AB1 - A0
ABPLUS1A = ABPLUS1A - ABPLUS0A
ABMIN1A  = ABMIN1A  - ABMIN0A

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS0A,A%NDimX,ABMIN1A,A%NDimX,0.0d0,A1A,A%NDimX)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS1A,A%NDimX,ABMIN0A,A%NDimX,1d0,A1A,A%NDimX)
print*, 'A1A',norm2(A1A)
!
!Calc: A2 = ABP1*ABM1
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS1A,A%NDimX,ABMIN1A,A%NDimX,0.0d0,A2A,A%NDimX)
deallocate(ABMIN1A)
print*, 'A2A',norm2(A2A)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeA(A%NDimX,NCholesky))
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUS0A,A%NDimX,A%DChol,NCholesky,0.0d0,ABP0TildeA,A%NDimX)
print*, 'APLUS0Tilde',norm2(ABP0TildeA)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeA(A%NDimX,NCholesky))
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUS1A,A%NDimX,A%DChol,NCholesky,0.0d0,ABP1TildeA,A%NDimX)
print*, 'APLUS1Tilde',norm2(ABP1TildeA)

deallocate(ABPLUS1A)
!
! monomer B
allocate(ABPLUS0B(B%NDimX,B%NDimX),ABMIN0B(B%NDimX,B%NDimX))
allocate(ABPLUS1B(B%NDimX,B%NDimX),ABMIN1B(B%NDimX,B%NDimX))
allocate(A1B(B%NDimX,B%NDimX),A2B(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1B
read(iunit) ABMIN1B

close(iunit)

open(newunit=iunit,file='A0MAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS0B
read(iunit) ABMIN0B

close(iunit)

! AB1 = AB1 - A0
ABPLUS1B = ABPLUS1B - ABPLUS0B
ABMIN1B  = ABMIN1B  - ABMIN0B

print*, 'ABPLUS1B',norm2(ABPLUS1B)
print*, 'ABMIN1B ',norm2(ABMIN1B)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS0B,B%NDimX,ABMIN1B,B%NDimX,0.0d0,A1B,B%NDimX)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS1B,B%NDimX,ABMIN0B,B%NDimX,1d0,A1B,B%NDimX)
print*, 'A1B',norm2(A1B)
!
!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS1B,B%NDimX,ABMIN1B,B%NDimX,0.0d0,A2B,B%NDimX)
deallocate(ABMIN1B)
print*, 'A2B',norm2(A2B)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeB(B%NDimX,NCholesky))
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUS0B,B%NDimX,B%DChol,NCholesky,0.0d0,ABP0TildeB,B%NDimX)
print*, 'APLUS0Tilde-b',norm2(ABP0TildeB)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeB(B%NDimX,NCholesky))
Call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUS1B,B%NDimX,B%DChol,NCholesky,0.0d0,ABP1TildeB,B%NDimX)
print*, 'APLUS1Tilde-b',norm2(ABP1TildeB)

deallocate(ABPLUS1B)
!
! get CAlphaTilde in an iterative manner

NFreq = 12
Max_Cn = SAPT%Max_Cn
print*, 'SAPT%NFreq =', NFreq
print*, 'SAPT%MaxCn =', Max_Cn

allocate(XFreq(NFreq),WFreq(NFreq))
!
allocate(C0TildeA(A%NDimX,NCholesky),C0TildeB(B%NDimX,NCholesky))
allocate(C1TildeA(A%NDimX,NCholesky),C1TildeB(B%NDimX,NCholesky))
allocate(C2TildeA(A%NDimX,NCholesky),C2TildeB(B%NDimX,NCholesky))
allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(LambdaA(A%NDimX,A%NDimX),LambdaB(B%NDimX,B%NdimX))
allocate(WorkA(A%NDimX,NCholesky),WorkB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! read ABPLUS0.ABMIN0 blocks
call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

!  Calc: LAMBDA=(A0+Om^2)^-1
   Call Inv_AC0Blk(OmI**2,LambdaA,A0BlkA,A0BlkIVA,nblkA,A%NDimX)
   Call Inv_AC0Blk(OmI**2,LambdaB,A0BlkB,A0BlkIVB,nblkB,B%NDimX)

!  Calc: C0Tilde=1/2 LAMBDA.APLUS0Tilde
   call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,0.5d0,LambdaA,A%NDimX,ABP0TildeA,A%NDimX,0.0d0,C0TildeA,A%NDimX)
   call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,0.5d0,LambdaB,B%NDimX,ABP0TildeB,B%NDimX,0.0d0,C0TildeB,B%NDimX)

   !print*, 'C0Tilde-A',OmI,norm2(C0TildeA)
   !print*, 'C0Tilde-B',OmI,norm2(C0TildeB)

!  Calc: C1Tilde=LAMBDA.(1/2 APLUS1Tilde - A1.C0Tilde)
   call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,A1A,A%NDimX,C0TildeA,A%NDimX,0.0d0,CTildeA,A%NDimX)
   CTildeA = 0.5d0*ABP1TildeA - CTildeA
   call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,LambdaA,A%NDimX,CTildeA,A%NDimX,0.0d0,C1TildeA,A%NDimX)

   call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,A1B,B%NDimX,C0TildeB,B%NDimX,0.0d0,CTildeB,B%NDimX)
   CTildeB = 0.5d0*ABP1TildeB - CTildeB
   call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,LambdaB,B%NDimX,CTildeB,B%NDimX,0.0d0,C1TildeB,B%NDimX)

   !print*, 'C1Tilde-A',OmI,norm2(C1TildeA)
   !print*, 'C1Tilde-B',OmI,norm2(C1TildeB)

   ! test uncoupled
   CTildeA = 0
   CTildeB = 0

   CTildeA = C0TildeA
   CTildeB = C0TildeB

   ! test semicoupled
   CTildeA = CTildeA + C1TildeA
   CTildeB = CTildeB + C1TildeB

   XFactorial = 1
   do N=2,Max_Cn

       XFactorial = XFactorial*N
       XN1 = -N
       XN2 = -N*(N-1)

       call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,XN2,A2A,A%NDimX,C0TildeA,A%NDimX,0.0d0,WorkA,A%NDimX)
       call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,XN1,A1A,A%NDimX,C1TildeA,A%NDimX,1.0d0,WorkA,A%NDimX)
       call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.0d0,LambdaA,A%NDimX,WorkA,A%NDimX,0.0d0,C2TildeA,A%NDimX)

       CTildeA  = CTildeA + C2TildeA / XFactorial
       C0TildeA = C1TildeA
       C1TildeA = C2TildeA

       !print*, 'N,COMTildeA',N,norm2(CTildeA)

       call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,XN2,A2B,B%NDimX,C0TildeB,B%NDimX,0.0d0,WorkB,B%NDimX)
       call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,XN1,A1B,B%NDimX,C1TildeB,B%NDimX,1.0d0,WorkB,B%NDimX)
       call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.0d0,LambdaB,B%NDimX,WorkB,B%NDimX,0.0d0,C2TildeB,B%NDimX)

       CTildeB  = CTildeB + C2TildeB / XFactorial
       C0TildeB = C1TildeB
       C1TildeB = C2TildeB

       !print*, 'N,COMTildeB',N,norm2(CTildeB)

   enddo

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

enddo

e2d = -32d0/Pi*e2d*1d3
!print*, '???e2d = ',e2d
call print_en('E2disp(CAlpha)-TEST',e2d,.false.)

deallocate(WFreq,XFreq)
deallocate(A2B,A1B,A2A,A1A)
deallocate(ABP1TildeA,ABP0TildeA,ABP1TildeB,ABP0TildeB)
deallocate(LambdaB,LambdaA)
deallocate(CTildeB,CTildeA)
deallocate(C2TildeB,C2TildeA)
deallocate(C1TildeB,C1TildeA)
deallocate(C0TildeB,C0TildeA)
deallocate(WorkB,WorkA)
deallocate(CB,CA)
!deallocate(DCholB,DCholA)

end subroutine e2disp_CAlphaTilde_full

subroutine get_Cmat(Cmat,CICoef,IndN,ABPM,ABPLUS,Omega,NDimX,NBas)
implicit none

integer,intent(in) :: NDimX,NBas,IndN(2,NBas)
double precision,intent(in)  :: CICoef(NBas),Omega
double precision,intent(in)  :: ABPM(NDimX,NDimX),ABPLUS(NDimX,NDimX)
double precision,intent(out) :: Cmat(NDimX,NDimX)

integer          :: info,lwork
integer          :: i,ipq,irs,ip,iq,ir,is
double precision :: fact_rs,fact_pq
integer,allocatable :: ipiv(:)
double precision,allocatable :: work(:,:),work1(:)

allocate(work(NDimX,NDimX),ipiv(NDimX))
work = ABPM
do i=1,NDimX
   work(i,i) = work(i,i) + Omega**2
enddo

! work=(ABPM+omega^2)-1.ABPLUS
!...

! this is slower
!call dgetrf(NDimX,NDimX,work,NDimX,ipiv,info)
!print*, 'info1',info
!allocate(work1(1))
!lwork = -1
!call dgetri(NDimX,work,NDimX,ipiv,work1,lwork,info)
!lwork = int(work1(1))
!print*, 'lwork:',lwork
!deallocate(work1)
!allocate(work1(lwork))
!call dgetri(NDimX,work,NDimX,ipiv,work1,lwork,info)
!print*, 'info2',info
!!
!! Cmat=work.ABPLUS
!call dgemm('N','N',NDimX,NDimX,NDimX,1d0,work,NDimX,ABPLUS,NDimX,0d0,Cmat,NDimX)
!
! deallocate(work1)
!
Cmat = ABPlus
call dgesv(NDimX,NDimX,work,NDimX,ipiv,Cmat,NDimX,info)
!print*, 'info',info
!print*, 'Cmat',norm2(Cmat)

!
do irs=1,NDimX
   ir = IndN(1,irs)
   is = IndN(2,irs)
   fact_rs = CICoef(ir)+CICoef(is)
   do ipq=1,NDimX
      ip = IndN(1,ipq)
      iq = IndN(2,ipq)
      fact_pq = CICoef(ip)+CICoef(iq)
      Cmat(ipq,irs) = fact_rs*fact_pq*Cmat(ipq,irs)
   enddo
enddo

deallocate(ipiv,work)

end subroutine get_Cmat

end module sapt_Chol_pol
