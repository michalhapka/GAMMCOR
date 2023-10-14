module sapt_utils
use types
use tran
use diis

contains

subroutine get_den(nbas,MO,Occ,Fac,Den)
implicit none

integer,intent(in) :: nbas
double precision, intent(in) :: MO(nbas,nbas)
double precision, intent(in) :: Occ(nbas)
double precision, intent(in) :: Fac
double precision, intent(out) :: Den(nbas,nbas)
integer :: i

Den = 0
do i=1,nbas
   call dger(nbas, nbas, Fac*Occ(i), MO(:, i), 1, MO(:, i), 1, Den, nbas)
enddo

end subroutine get_den

subroutine get_one_mat(var,mat,mono,nbas)
!
! Purpose: get 1-el V, S or H matrix from ONEEL files
!
implicit none

character(1),intent(in)      :: var
integer,intent(in)           :: nbas,mono
double precision,intent(out) :: mat(nbas,nbas)

integer                  :: ione
logical                  :: valid
character(8)             :: label
character(:),allocatable :: onefile

 if(mono==1) then
    onefile = 'ONEEL_A'
 elseif(mono==2) then
    onefile = 'ONEEL_B'
 else
    write(LOUT,'(1x,a)') 'ERROR!!! ONLY 2 MONOMERS ACCEPTED!'
    stop
 endif

 open(newunit=ione,file=onefile,access='sequential',&
      form='unformatted',status='old')

 valid=.false.
 mat=0
 select case(var)
 case('V','v')

    read(ione)
    read(ione) label,mat
    if(label=='POTENTAL') valid=.true.

 case('S','s')

    read(ione) label,mat
    if(label=='OVERLAP ') valid=.true.

 case('H','h')

    read(ione)
    read(ione)
    read(ione) label,mat
    if(label=='ONEHAMIL') valid=.true.

 case default
    write(LOUT,'()')
    write(LOUT,'(1x,a)') 'ERROR IN get_one_max! TYPE '//var//' NOT AVAILABLE!'
    stop
 end select

 if(.not.valid) then
    write(LOUT,'(1x,a)') 'ERROR!!! LABEL MISMATCH IN get_one_mat!'
    stop
 endif

 close(ione)

end subroutine get_one_mat


subroutine read_SBlock(SBlock,SBlockIV,nblk,xy0file)
implicit none

type(EBlockData)             :: SBlockIV
type(EBlockData),allocatable :: SBlock(:)
integer,intent(out)          :: nblk
character(*),intent(in)      :: xy0file

integer                      :: iunit,i,ival,iblk

! read SBLOCK
open(newunit=iunit,file=xy0file,status='OLD',&
     form='UNFORMATTED')
read(iunit) nblk
!print*, 'ival,nblk',ival,nblk
allocate(SBlock(nblk))
do iblk=1,nblk
   associate(A => SBlock(iblk))
     read(iunit) i, A%n, A%l1, A%l2
     allocate(A%pos(A%n),A%matX(A%n,A%n),A%matY(A%n,A%n),A%vec(A%n))
     read(iunit) A%pos,A%matX,A%matY,A%vec
   end associate
enddo
associate(A => SBlockIV)
  read(iunit) A%n,A%l1,A%l2
  allocate(A%pos(A%n),A%vec(A%n))
  allocate(A%matX(A%n,1),A%matY(A%n,1))
  read(iunit) A%pos,A%vec
  read(iunit) A%matX,A%matY
end associate
close(iunit)

end subroutine read_SBlock

subroutine read_ABPM0Block(A0Block,A0BlockIV,nblk,abpm0file)
implicit none

type(EBlockData)             :: A0BlockIV
type(EBlockData),allocatable :: A0Block(:)
integer,intent(out)          :: nblk
character(*),intent(in)      :: abpm0file

integer                      :: iunit,i,ival,iblk

! read ABPM0 = ABPLUS0.ABMIN0
open(newunit=iunit,file=abpm0file,status='OLD',&
     form='UNFORMATTED')
read(iunit) nblk
allocate(A0Block(nblk))
do iblk=1,nblk
   associate(A => A0Block(iblk))
     read(iunit) i, A%n, A%l1, A%l2
     allocate(A%pos(A%n),A%matX(A%n,A%n))
     read(iunit) A%pos,A%matX
   end associate
enddo
associate(A => A0BlockIV)
  read(iunit) A%n,A%l1,A%l2
  allocate(A%pos(A%n),A%vec(A%n))
  read(iunit) A%pos,A%vec
end associate
close(iunit)

end subroutine read_ABPM0Block

subroutine create_Y0tilde(Y01Block,SBlock,SBlockIV,Eig,CICoef,IndN,nblk,NBasis,NDimX)
implicit none

type(EBlockData)   :: Sblock(nblk),SblockIV
type(Y01BlockData) :: Y01Block(NDimX)

integer,intent(in)           :: nblk,NBasis,NDimX
integer,intent(in)           :: IndN(2,NDimX)
double precision,intent(in)  :: CICoef(NBasis)
double precision,intent(out) :: Eig(NDimX)

integer          :: i,ii,ipos,iblk,ip,iq
double precision :: fact

do iblk=1,nblk
   associate(A => SBlock(iblk))

     do i=1,A%n
        ipos = A%pos(i)

        ip   = IndN(1,ipos)
        iq   = IndN(2,ipos)
        fact = CICoef(ip)-CICoef(iq)
        associate(Y => Y01Block(ipos))

          Y%n  = A%n
          Y%l1 = A%l1
          Y%l2 = A%l2
          allocate(Y%vec0(Y%n))
          Y%vec0(1:Y%n) = fact*(A%matY(i,1:Y%n)-A%matX(i,1:Y%n))

        end associate
     enddo
     Eig(A%l1:A%l2) = A%vec(1:A%n)
   end associate
enddo
associate(A => SBlockIV)

  do i=1,A%n
     ii = A%l1+i-1
     ipos = A%pos(i)

     associate(Y => Y01Block(ipos))

       Y%n = 1
       Y%l1 = ii
       Y%l2 = ii
       allocate(Y%vec0(1))
       Y%vec0(1) = 1d0/sqrt(2d0)

     end associate

     Eig(ii) = A%vec(i)
  enddo

end associate

end subroutine create_Y0tilde

subroutine convert_XY0_to_Y01(Mon,Y01Block,Eig0,NBasis,xy0file)
implicit none

type(SystemBlock)              :: Mon
type(EblockData)               :: SBlockIV
type(EblockData),allocatable   :: SBlock(:)
type(Y01BlockData)             :: Y01Block(Mon%NDimX)

integer,intent(in)           :: NBasis
character(*),intent(in)      :: xy0file
double precision,intent(out) :: Eig0(Mon%NDimX)

integer :: iunit
integer :: i,iblk,nblk,ival,ipos

call read_SBlock(SBlock,SBlockIV,nblk,xy0file)

call create_Y0tilde(Y01Block,SBlock,SBlockIV,Eig0,Mon%CICoef,Mon%IndN,nblk,NBasis,Mon%NDimX)

! deallocate SBLOCK
 do iblk=1,nblk
    associate(A => SBlock(iblk))
      deallocate(A%matY,A%matX,A%vec)
      deallocate(A%pos)
    end associate
 enddo
 ! deallocate IV part
 associate(A => SblockIV)
   if(A%n>0) then
      deallocate(A%vec)
      deallocate(A%pos)
   endif
 end associate

end subroutine convert_XY0_to_Y01

subroutine unpack_XY0_full(VecX0,VecY0,Eig0,C,IndN,NDimX,NBasis,xy0file)
implicit none

integer,intent(in)           :: NDimX,NBasis
integer,intent(in)           :: IndN(2,NDimX)
double precision,intent(in)  :: C(NBasis)
character(*),intent(in)      :: xy0file
double precision,intent(out) :: Eig0(NDimX), &
                                VecX0(NDimX,NDimX),VecY0(NDimX,NDimX)

integer                      :: iunit
integer                      :: i,ii,ip,iq,iblk,nblk,ipos
double precision             :: valX,valY
type(EblockData),allocatable :: SBlock(:)
type(EblockData)             :: SBlockIV
double precision,parameter   :: fac = 1d0/sqrt(2d0)

! read SBLOCK
open(newunit=iunit,file=xy0file,status='OLD',&
     form='UNFORMATTED')
read(iunit) nblk
allocate(SBlock(nblk))
!print*, 'ival,nblkA',ival,nblkA
do iblk=1,nblk
   associate(A => SBlock(iblk))
     read(iunit) i, A%n, A%l1, A%l2
     allocate(A%pos(A%n),A%matX(A%n,A%n),A%matY(A%n,A%n),A%vec(A%n))
     read(iunit) A%pos,A%matX,A%matY,A%vec
   end associate
enddo
associate(A => SBlockIV)
  read(iunit) A%n,A%l1,A%l2
  allocate(A%pos(A%n),A%vec(A%n))
  read(iunit) A%pos,A%vec
end associate
close(iunit)

VecX0 = 0
VecY0 = 0
Eig0  = 0

do iblk=1,nblk
   associate(B => Sblock(iblk))

     do i=1,B%n
        ipos = B%pos(i)
        VecY0(ipos,B%l1:B%l2) = B%matY(i,1:B%n)
        VecX0(ipos,B%l1:B%l2) = B%matX(i,1:B%n)
     enddo
     Eig0(B%l1:B%l2) = B%vec(1:B%n)

   end associate
enddo

!unpack (2, IV part)
associate(B => SblockIV)
  if(B%n>0) then
     do i=1,B%n
        ii = B%l1+i-1
        ipos = B%pos(i)
        ip = IndN(1,ipos)
        iq = IndN(2,ipos)

        valX = 1d0/(C(ip)+C(iq))
        valY = 1d0/(C(ip)-C(iq))

        VecX0(ipos,ii) = 0.5d0*fac*(valX-valY)
        VecY0(ipos,ii) = 0.5d0*fac*(valX+valY)
        Eig0(ii) = B%vec(i)
     enddo
  endif

end associate

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

end subroutine unpack_XY0_full

subroutine convert_XY_to_Z(EVecZ,CICoef,IndN,NDimX,NBasis,xyfile)
implicit none

integer,intent(in)           :: NDimX,NBasis
integer,intent(in)           :: IndN(2,NDimX)
character(*),intent(in)      :: xyfile
double precision,intent(in)  :: CICoef(NBasis)
double precision,intent(out) :: EVecZ(NDimX,NDimX)

integer                      :: i,j,ip,iq
double precision             :: fact
double precision,allocatable :: EVecX(:,:),EVecY(:,:)

allocate(EVecX(NDimX,NDimX),EVecY(NDimX,NDimX))

call readEVecXY(EVecX,EVecY,NDimX,xyfile)

do j=1,NDimX
   do i=1,NDimX
      ip = IndN(1,i)
      iq = IndN(2,i)
      fact = CICoef(ip) - CICoef(iq)
      EVecZ(i,j) = fact*(EVecY(i,j)-EVecX(i,j))
   enddo
enddo

deallocate(EVecY,EVecX)

end subroutine convert_XY_to_Z

subroutine writeampl(sij,ampfile)
implicit none

character(*) :: ampfile
double precision :: sij(:,:)
integer :: iunit

 open(newunit=iunit,file=ampfile,form='unformatted')
 write(iunit) sij
 close(iunit)

end subroutine writeampl

subroutine solve_cphf(M,WPot,e2indxy,Flags,NBas)
implicit none

type(SystemBlock) :: M
type(FlagsData)   :: Flags

integer,intent(in)           :: NBas
double precision,intent(in)  :: WPot(NBas,NBas)
double precision,intent(out) :: e2indxy

type(DIISData)                 :: DIISBlock
type(EblockData)               :: ABlockIV
type(EblockData),allocatable   :: ABlock(:)

double precision,allocatable :: WxYY(:,:)
double precision,allocatable :: wVecxYY(:)
double precision,allocatable :: amps(:),vecR(:),delta(:)
double precision,allocatable :: OmM0(:)

integer      :: iunit
integer      :: iter,i,pq,ip,iq
integer      :: iblk,nblk,nmax
logical      :: conv=.FALSE.
character(8) :: nameunc,abfile

double precision :: fact, error
integer,parameter          :: MaxIt = 20
double precision,parameter :: ThrDIIS = 1.d-8

double precision,allocatable :: ABMin(:,:),Work(:)

 if(M%Monomer==1) nameunc='XY0_A'
 if(M%Monomer==2) nameunc='XY0_B'

 if(M%Monomer==1) abfile='ABMAT_A'
 if(M%Monomer==2) abfile='ABMAT_B'

 allocate(ABMin(M%NDimX,M%NDimX))

 open(newunit=iunit,file=abfile,status='OLD',&
      access='SEQUENTIAL',form='UNFORMATTED')

 read(iunit)
 read(iunit) ABMin

 close(iunit)

 call init_DIIS(DIISBlock,M%NDimX,M%NDimX,Flags%DIISN)

 allocate(WxYY(NBas,NBas))
 allocate(wVecxYY(M%NDimX))
 call tran2MO(WPot,M%CMO,M%CMO,WxYY,NBas)
 wVecxYY = 0
 do pq=1,M%NDimX
    ip = M%IndN(1,pq)
    iq = M%IndN(2,pq)
    fact = M%CICoef(ip)+M%CICoef(iq)
    wVecxYY(pq) = fact*WxYY(ip,iq)
 enddo

 call convert_to_ABMIN(ABlock,ABlockIV,M%IndN,M%CICoef,nblk,NBas,M%NDimX,trim(nameunc))

 allocate(vecR(M%NDimX),amps(M%NDimX),delta(M%NDimX))

 ! zero iter
 !call amplitudes_T1_cphf(OmM0,wVecxYY,amps,M%NDimX)
 call amplitudes_T1_cerpa(ABlock,ABlockIV,nblk,wVecxYY,amps,M%NDimX)

 !block
 !print*, 'amp-1',norm2(amps)
 !do pq=1,10
 !   print*, 'pq',pq,amps(pq)
 !enddo
 !end block

 e2indxy = 0
 do pq=1,M%NDimX
    ip = M%IndN(1,pq)
    iq = M%IndN(2,pq)
    fact = M%CICoef(ip)+M%CICoef(iq)
    e2indxy = e2indxy + fact*amps(pq)*WxYY(ip,iq)
 enddo
 e2indxy = 1d0*e2indxy
 print*, 'e2ind(unc)',e2indxy

 write(LOUT,'(1x,a,5x,a)') 'ITER', 'ERROR'
 iter = 0
 do

 vecR = wVecxYY
 !call dgemv('N',M%NDimX,M%NDimX,1.d0,M%PP,M%NDimX,amps,1,1d0,vecR,1)
 call dgemv('N',M%NDimX,M%NDimX,1.d0,ABMin,M%NDimX,amps,1,1d0,vecR,1)

 error = norm2(vecR)
 conv=error.le.ThrDIIS

 write(LOUT,'(1x,i3,f11.6)') iter,error

 if(conv) then
    exit
 elseif(.not.conv.and.iter.le.MaxIt) then
    iter = iter + 1
 elseif(.not.conv.and.iter.gt.MaxIt) then
   write(*,*) 'Error!!! E2ind DIIS not converged!'
   exit
 endif

 ! HF test
 !call amplitudes_T1_cphf(OmM0,vecR,delta,M%NDimX)
 call amplitudes_T1_cerpa(ABlock,ABlockIV,nblk,vecR,delta,M%NDimX)
 amps = amps + delta
 if(iter>Flags%DIISOn) call use_DIIS(DIISBlock,amps,vecR)

 enddo

 !E2ind(X--Y)
 e2indxy = 0
 do pq=1,M%NDimX
    ip = M%IndN(1,pq)
    iq = M%IndN(2,pq)
    fact = M%CICoef(ip)+M%CICoef(iq)
    e2indxy = e2indxy + fact*amps(pq)*WxYY(ip,iq)
 enddo
 e2indxy = 2d0*e2indxy
 print*, 'e2indxy',e2indxy

 call free_DIIS(DIISBlock)

 ! deallocate ABLOCK
 do iblk=1,nblk
    associate(A => ABlock(iblk))
      deallocate(A%matY,A%matX,A%vec)
      deallocate(A%pos)
    end associate
 enddo
 ! deallocate IV part
 associate(A => AblockIV)
   if(A%n>0) then
      deallocate(A%vec)
      deallocate(A%pos)
   endif
 end associate

 deallocate(ABmin)
 deallocate(delta,amps,vecR)
 deallocate(ABlock)
 deallocate(wVecxYY,WxYY)

end subroutine solve_cphf

subroutine amplitudes_T1_cphf(deps,ints,res,NDimX)
implicit none

integer,intent(in) :: NDimX
double precision,intent(in)  :: deps(NDimX),ints(NDimX)
double precision,intent(out) :: res(NDimX)

integer :: pq,ip,iq

res = 0
do pq=1,NDimX
   res(pq) = res(pq) - ints(pq) / deps(pq)
enddo

end subroutine amplitudes_T1_cphf

subroutine amplitudes_T1_cerpa(ABlock,ABlockIV,nblk,ints,res,NDimX)
implicit none

type(EBlockData)   :: ABlock(nblk),AblockIV

integer,intent(in) :: nblk,NDimX
double precision,intent(in)  :: ints(NDimX)
double precision,intent(out) :: res(NDimX)

integer :: iblk,ipos,pq,ip,iq
integer :: i,nmax,info
double precision,allocatable :: tmp(:)

! set nmax
nmax = ABlockIV%n
do iblk=1,nblk
   nmax = max(nmax,ABlock(iblk)%n)
enddo

allocate(tmp(nmax))

res = 0
tmp = 0
do iblk=1,nblk
   associate(A => ABlock(iblk))

     do i=1,A%n
        ipos = A%pos(i)
        tmp(i) = ints(ipos)
     enddo

     call dsytrs('U',A%n,1,A%matY,A%n,A%ipiv,tmp,A%n,info)
!     print*, 'info',info,norm2(A%matY)

     do i=1,A%n
        ipos = A%pos(i)
       ! print*, 'loop',iblk,norm2(tmp)
        res(ipos) = -tmp(i)
     enddo

   end associate
enddo
! block IV
associate(A => ABlockIV)

   do i=1,A%n
     ipos = A%pos(i)
     res(ipos) = -2d0*ints(ipos) / A%vec(i)
  enddo

end associate

deallocate(tmp)

end subroutine amplitudes_T1_cerpa

subroutine Sblock_to_ABMAT(ABlock,ABlockIV,IndN,CICoef,nblk,NBas,NDimX,xy0file)
implicit none

type(EBlockData)             :: ABlockIV
type(EBlockData),allocatable :: Ablock(:)

integer,intent(in)  :: NBas,NDimX
integer,intent(in)  :: IndN(2,NBas)
integer,intent(out) :: nblk
double precision,intent(in) :: CICoef(NBas)
character(*),intent(in)     :: xy0file

type(EBlockData)             :: SblockIV
type(EBlockData),allocatable :: Sblock(:)

integer :: i,iblk,ipos,ip,iq
integer :: j,info
double precision  :: fact,valX,valY
double precision,allocatable :: tmpX(:,:),tmpY(:,:),work(:)
! test
double precision :: tstX,tstY

! why factor sqrt(2)?
fact = sqrt(2d0)

call read_SBlock(SBlock,SBlockIV,nblk,xy0file)

allocate(ABlock(nblk))

tstX = 0
tstY = 0
do iblk=1,nblk
   associate(A => SBlock(iblk),&
             C => ABlock(iblk))

     C%n = A%n
     allocate(C%pos(C%n),C%ipiv(C%n),C%vec(C%n))
     allocate(C%matX(C%n,C%n),C%matY(C%n,C%n))

     C%l1  = A%l1
     C%l2  = A%l2
     C%pos = A%pos
     C%vec = A%vec

     allocate(tmpX(A%n,A%n),tmpY(A%n,A%n),work(A%n))

     do i=1,A%n
        ipos = A%pos(i)
        ip   = IndN(1,ipos)
        iq   = IndN(2,ipos)
        valX = fact*(CICoef(ip)+CICoef(iq))
        valY = fact*(CICoef(ip)-CICoef(iq))
        C%matX(i,1:C%n) = valX*(A%matX(i,1:A%n) + A%matY(i,1:A%n))
        C%matY(i,1:C%n) = valY*(A%matY(i,1:A%n) - A%matX(i,1:A%n))
     enddo

     do i=1,A%n
        tmpX(:,i) = C%matX(:,i) * C%vec(i)
     enddo
     tmpY(1:A%n,1:A%n) = C%matX(1:A%n,1:A%n)

     ! ABMIN in matX
     call dgemm('N','T',A%n,A%n,A%n,1d0,tmpX,A%n,tmpY,A%n,0d0,C%matX,A%n)

     tmpX = 0
     tmpY = 0
     do i=1,A%n
        tmpY(:,i) = C%matY(:,i) * C%vec(i)
     enddo
     tmpX(1:A%n,1:A%n) = C%matY(1:A%n,1:A%n)

     ! ABPLUS in matY
     call dgemm('N','T',A%n,A%n,A%n,1d0,tmpY,A%n,tmpX,A%n,0d0,C%matY,A%n)

     do j=1,A%n
     do i=1,A%n
     tstX = tstX + C%matX(i,j)**2
     tstY = tstY + C%matY(i,j)**2
     enddo
     enddo

     deallocate(work,tmpX,tmpY)

   end associate
enddo

ABlockIV = SBlockIV

! deallocate SBLOCK
do iblk=1,nblk
   associate(A => SBlock(iblk))
     deallocate(A%matY,A%matX,A%vec)
     deallocate(A%pos)
   end associate
enddo
! deallocate IV part
associate(A => SblockIV)
  if(A%n>0) then

    ! test IV
    do i=1,A%n
       tstX = tstX + (1d0*A%vec(i))**2
       tstY = tstY + (1d0*A%vec(i))**2
    enddo

    deallocate(A%vec)
    deallocate(A%pos)

  endif
end associate

!print*, 'something is wrong with this spectral subroutine'
!print*, 'tst-Y(ABPLS)',sqrt(tstY)
!print*, 'tst-X(ABMIN)',sqrt(tstX)

end subroutine SBlock_to_ABMAT

subroutine add_blk_right(AMAT,A0Blk,A0BlkIV,fact,isX,nblk,NDimX)
!
! subtract AMAT = AMAT - ABLOCK
!
implicit none

integer,intent(in) :: nblk,NDimX
double precision,intent(in) :: fact
logical,intent(in) :: isX
double precision,intent(inout) :: AMAT(NDimX,NDimX)

type(EblockData)    :: A0Blk(nblk), A0BlkIV

integer :: i,j
integer :: iblk,ipos,jpos

!print*, 'add: before blocks'
do iblk=1,nblk
   associate(B => A0Blk(iblk))

     do j=1,B%n
        jpos=B%pos(j)
        do i=1,B%n
           ipos=B%pos(i)
           if(isX) AMAT(ipos,jpos) = AMAT(ipos,jpos) + fact * B%matX(i,j)
           if(.not.isX) AMAT(ipos,jpos) = AMAT(ipos,jpos) + fact * B%matY(i,j)
        enddo
     enddo

   end associate
enddo

!IV block
! MH 19.10.2022 change below is a bypass of some weird bug
!               which is to be fixed!  
associate(B => A0BlkIV)
  do i=1,B%n
     !j = B%pos(i)
     !AMAT(j,j) = AMAT(j,j) + fact * B%vec(i)
     AMAT(B%pos(i),B%pos(i)) = AMAT(B%pos(i),B%pos(i)) + fact * B%vec(i)
  enddo
end associate

end subroutine add_blk_right

subroutine test_ABMAT(ABP,ABM,A0Blk,A0BlkIV,nblk,NBasis,NDimX)
! for testing
! fill ABPLUS, ABMIN from blocks

integer,intent(in) :: nblk,NBasis,NDimX

type(EBlockData)             :: A0BlkIV,A0Blk(nblk)

double precision,intent(out) :: ABP(NDimX,NDimX),ABM(NDimX,NDimX)

integer :: i,j,ii,iblk,jpos,ipos

ABP = 0
ABM = 0
do iblk=1,nblk
   associate(B => A0Blk(iblk))

     ! ABPLUS in matY
     ! ABMIN  in matX
     do j=1,B%n
        jpos=B%pos(j)
        do i=1,B%n
           ipos=B%pos(i)
           ABP(ipos,jpos) = B%matY(i,j)
           ABM(ipos,jpos) = B%matX(i,j)
        enddo
     enddo

   end associate
enddo
!IV block
associate(B => A0BlkIV)
  if(B%n>0) then
     do i=1,B%n
       ii = B%pos(i)
       ABP(ii,ii) = 1.0d0* B%vec(i)
       ABM(ii,ii) = 1.0d0* B%vec(i)
     enddo
  endif
end associate

print*, 'ABP',norm2(ABP)
print*, 'ABM',norm2(ABM)

end subroutine test_ABMAT

subroutine convert_to_ABMIN(ABlock,ABlockIV,IndN,CICoef,nblk,NBas,NDimX,xy0file)
implicit none

type(EBlockData)             :: ABlockIV
type(EBlockData),allocatable :: Ablock(:)

integer,intent(in)  :: NBas,NDimX
integer,intent(in)  :: IndN(2,NBas)
integer,intent(out) :: nblk
double precision,intent(in) :: CICoef(NBas)
character(*),intent(in) :: xy0file

type(EBlockData)             :: SblockIV
type(EBlockData),allocatable :: Sblock(:)

integer :: i,iblk,ipos,ip,iq
integer :: j,info
double precision  :: valX,valY
double precision,allocatable :: tmp(:,:),work(:)

call read_SBlock(SBlock,SBlockIV,nblk,xy0file)

allocate(ABlock(nblk))

do iblk=1,nblk
   associate(A => SBlock(iblk),&
             C => ABlock(iblk))

     C%n = A%n
     allocate(C%pos(C%n),C%ipiv(C%n),C%vec(C%n))
     allocate(C%matX(C%n,C%n),C%matY(C%n,C%n))

     C%l1  = A%l1
     C%l2  = A%l2
     C%pos = A%pos
     C%vec = A%vec

     allocate(tmp(A%n,A%n),work(A%n))

     do i=1,A%n
        ipos = A%pos(i)
        ip   = IndN(1,ipos)
        iq   = IndN(2,ipos)
        valX = (CICoef(ip)+CICoef(iq))
        C%matX(i,1:C%n) = valX*(A%matX(i,1:A%n) + A%matY(i,1:A%n))
     enddo

     do i=1,A%n
        tmp(:,i) = C%matX(:,i) * C%vec(i)
     enddo
     call dgemm('N','T',A%n,A%n,A%n,1d0,tmp,A%n,C%matX,A%n,0d0,C%matY,A%n)

     ! CmatY: inverse ABMIN
     !call dgetrf(A%n,A%n,C%matY,A%n,Ipiv,Info)
     !call dgetri(A%n,C%matY,A%n,Ipiv,Work,A%n,Info)

     !! try diagonalization
     !C%matX = C%matY
     !call Diag8(C%matX,A%n,A%n,A%vec,work)

     ! CmatY: factorized ABMIN
     call dsytrf('U',C%n,C%matY,C%n,C%ipiv,Work,C%n,Info)

     deallocate(work,tmp)

   end associate
enddo

ABlockIV = SBlockIV

! deallocate SBLOCK
do iblk=1,nblk
   associate(A => SBlock(iblk))
     deallocate(A%matY,A%matX,A%vec)
     deallocate(A%pos)
   end associate
enddo
! deallocate IV part
associate(A => SblockIV)
  if(A%n>0) then
     deallocate(A%vec)
     deallocate(A%pos)
  endif
end associate

end subroutine convert_to_ABMIN

subroutine readresp(EVec,EVal,NDim,fname)
implicit none

integer :: NDim
double precision :: EVec(NDim,NDim), EVal(NDim)
character(*) :: fname
integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVec
 read(iunit) EVal

 close(iunit)

end subroutine readresp

subroutine readEval(EVal,NDim,fname)
implicit none

integer :: NDim
double precision :: EVal(NDim)
character(*) :: fname
integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVal

 close(iunit)

end subroutine readEval

subroutine readEvalZ(EVal,NDim,fname)
implicit none

integer :: NDim
double precision :: EVal(NDim)
character(*) :: fname
integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit)
 read(iunit) EVal

 close(iunit)

end subroutine readEvalZ

subroutine readEvalXY(EVal,NDim,fname)
implicit none

integer,intent(in)           :: NDim
character(*),intent(in)      :: fname
double precision,intent(out) :: EVal(NDim)

integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit)
 read(iunit)
 read(iunit) EVal

 close(iunit)

end subroutine readEvalXY

subroutine readEvecZ(EVecZ,NDim,fname)
implicit none

integer          :: NDim
character(*)     :: fname
double precision :: EVecZ(NDim,NDim)

integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVecZ

 close(iunit)

end subroutine readEvecZ

subroutine readEvecXY(EVecX,EVecY,NDim,fname)
implicit none

integer,intent(in)           :: NDim
character(*),intent(in)      :: fname
double precision,intent(out) :: EVecX(NDim,NDim),EVecY(NDim,NDim)

integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVecX
 read(iunit) EVecY

 close(iunit)

end subroutine readEvecXY

subroutine unpack_Eig(SBlock,SBlockIV,nblk,Eig,NDimX)
implicit none

type(EBlockData)             :: SBlock(nblk),SBlockIV
integer,intent(in)           :: nblk,NDimX
double precision,intent(out) :: Eig(NDimX)

integer :: i,ii,iblk

do iblk=1,nblk
   associate(B => SBlock(iblk))
     Eig(B%l1:B%l2) = B%vec(1:B%n)
   end associate
enddo
associate(B => SBlockIV)
  do i=1,B%n
     ii = B%l1+i-1
     Eig(ii) = B%vec(i)
  enddo
end associate

end subroutine unpack_Eig

subroutine fill_Fmat(Fmat,Occ,NBas,variant)
implicit none

integer,intent(in)           :: NBas,variant
double precision,intent(in)  :: Occ(NBas)
double precision,intent(out) :: Fmat(NBas,NBas)

integer :: i,j

Fmat = 0d0
select case(variant)
case(0)
   ! HF
   do j=1,NBas
      do i=1,NBas
         Fmat(i,j) = Occ(i)*Occ(j)
      enddo
   enddo
case(1,11)
   ! BB functional
   do j=1,NBas
      do i=1,NBas
         Fmat(i,j) = sqrt(Occ(i)*Occ(j))
      enddo
   enddo
case(3)
   ! Power functional
   write(LOUT,*) 'POWER FUNCITONAL NOT READY YET!'
   stop

end select

end subroutine fill_Fmat

subroutine calc_resp(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
double precision,intent(in) :: EVec(Mon%NDimX,Mon%NDimX), &
                               EVal(Mon%NDimX)
double precision,intent(in) :: Freq
double precision,intent(out) :: Alpha(Mon%NDimX,Mon%NDimX)
double precision :: frac
integer :: pq,rs,t,ip,iq,ir,is
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8

 Alpha = 0
 do t=1,Mon%NDimX

!   if(imag) then
!      v=4d0*omega(p)/(freq**2+omega(p)**2)
!   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
!   end if
    if(EVal(t).gt.SmallE.and.EVal(t).lt.BigE) then

    frac = 8d0*EVal(t)/(EVal(t)**2d0+freq)

    do pq=1,Mon%NDimX
       ip = Mon%IndN(1,pq)
       iq = Mon%IndN(2,pq)
       do rs=1,pq
          ir = Mon%IndN(1,rs)
          is = Mon%IndN(2,rs)

          Alpha(pq,rs) = Alpha(pq,rs) + &
                       (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                       (Mon%CICoef(is)+Mon%CICoef(ir)) * &
                       EVec(pq,t)*EVec(rs,t)*frac
          Alpha(rs,pq) = Alpha(pq,rs)

       enddo
    enddo
    endif
 enddo

end subroutine calc_resp

subroutine calc_resp_apsg(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
!double precision,intent(in) :: EVec(2*(Mon%NDimX+Mon%NDimN),2*(Mon%NDimX+Mon%NDimN)), &
!                               EVal(2*(Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: EVec((Mon%NDimX+Mon%NDimN),(Mon%NDimX+Mon%NDimN)), &
                               EVal((Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: Freq
!double precision,intent(out) :: Alpha(2*(Mon%NDimX+Mon%NDimN),2*(Mon%NDimX+Mon%NDimN))
double precision,intent(out) :: Alpha((Mon%NDimX+Mon%NDimN),(Mon%NDimX+Mon%NDimN))
double precision :: frac
integer :: pq,rs,t,ip,iq,ir,is
integer :: NDimEx
double precision,parameter :: SmallE = 1.d-6

 NDimEx = (Mon%NDimX + Mon%NDimN)

 Alpha = 0
 ! p>q, r>s
 do t=1,NDimEx

!   if(imag) then
!      v=4d0*omega(p)/(freq**2+omega(p)**2)
!   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
!   end if

    if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
       frac = 8d0*EVal(t)/(EVal(t)**2+freq)

       do pq=1,Mon%NDimX
          ip = Mon%IndN(1,pq)
          iq = Mon%IndN(2,pq)
          do rs=1,pq
             ir = Mon%IndN(1,rs)
             is = Mon%IndN(2,rs)

             Alpha(pq,rs) = Alpha(pq,rs) + &
                          (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                          (Mon%CICoef(is)+Mon%CICoef(ir)) * &
                          EVec(pq,t)*EVec(rs,t)*frac
             Alpha(rs,pq) = Alpha(pq,rs)

          enddo
       enddo
    endif
 enddo

! p>q,r=s
!
  do t=1,NDimEx
 !   if(imag) then
 !      v=4d0*omega(p)/(freq**2+omega(p)**2)
 !   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
 !   end if

     if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
        frac = 8d0*EVal(t)/(EVal(t)**2+freq)

        do pq=1,Mon%NDimX
           ip = Mon%IndN(1,pq)
           iq = Mon%IndN(2,pq)
           do ir=1,Mon%NDimN

!              Alpha(pq,ir) = Alpha(pq,ir) + &
!                           (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
!                          Mon%CICoef(ir) * &
!                           EVec(pq,t)*EVec(2*Mon%NDimX+ir,t)*frac
!              Alpha(ir,pq) = Alpha(pq,ir)

              Alpha(pq,Mon%NDimX+ir) = Alpha(pq,Mon%NDimX+ir) + &
                           (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                           Mon%CICoef(ir) * &
                           !EVec(pq,t)*EVec(2*Mon%NDimX+ir,t)*frac
                           EVec(pq,t)*EVec(Mon%NDimX+ir,t)*frac
              Alpha(Mon%NDimX+ir,pq) = Alpha(pq,Mon%NDimX+ir)

           enddo
        enddo

     endif
  enddo

 print*, 'aaa?',2*Mon%NDimX+Mon%NDimN,NDimEx

 ! p=q,r=s
  do t=1,NDimEx

 !   if(imag) then
 !      v=4d0*omega(p)/(freq**2+omega(p)**2)
 !   else
 !      v=4d0*omega(p)/(-freq**2+omega(p)**2)
 !   end if

     if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
        frac = 8d0*EVal(t)/(EVal(t)**2+freq)

        do ip=1,Mon%NDimN
           do ir=1,ip !Mon%NDimN

!              Alpha(ip,ir) = Alpha(ip,ir) + &
!                           Mon%CICoef(ip) * &
!                           Mon%CICoef(ir) * &
!                           EVec(2*Mon%NDimX+ip,t)*EVec(2*Mon%NDimX+ir,t)*frac
!              Alpha(ir,ip) = Alpha(ip,ir)

              Alpha(Mon%NDimX+ip,Mon%NDimX+ir) = Alpha(Mon%NDimX+ip,Mon%NDimX+ir) + &
                           Mon%CICoef(ip) * &
                           Mon%CICoef(ir) * &
                           !EVec(2*Mon%NDimX+ip,t)*EVec(2*Mon%NDimX+ir,t)*frac
                           EVec(Mon%NDimX+ip,t)*EVec(Mon%NDimX+ir,t)*frac
              Alpha(Mon%NDimX+ir,Mon%NDimX+ip) = Alpha(Mon%NDimX+ip,Mon%NDimX+ir)


           enddo
        enddo
     endif
  enddo

print*, 'resp,Alpha',norm2(Alpha)

end subroutine calc_resp_apsg

subroutine calc_resp_apsg2(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
double precision,intent(in) :: EVec((Mon%NDimX+Mon%NDimN)**2), &
                               EVal((Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: Freq
double precision,intent(out) :: Alpha((Mon%NDimX+Mon%NDimN),(Mon%NDimX+Mon%NDimN))
double precision :: frac,fact
integer,allocatable :: MIndEx(:,:)
double precision,allocatable :: VecEx(:)
integer :: pq,rs,t,ip,iq,ir,is
integer :: i
integer :: NDimEx
integer :: coef
double precision,parameter :: SmallE = 1.d-6, BigE = 1.d20

 NDimEx = (Mon%NDimX + Mon%NDimN)
 coef = 1
 ! with PINOVEC
 !coef = 2

 allocate(VecEx(NDimEx**2),MIndEx(2,NDimEx))

 MIndEx = Mon%IndNx

 VecEx = 0
 ! prepare extended vec
 do pq=1,NDimEx
    if(pq<=Mon%NDimX) then

        ip = Mon%IndN(1,pq)
        iq = Mon%IndN(2,pq)

        fact = Mon%CICoef(iq) + Mon%CICoef(ip)
        do i=1,coef*NDimEx
           VecEx((i-1)*coef*NDimEx+pq) = fact * &
                                          EVec((i-1)*coef*NDimEx+pq)
        enddo

    elseif(pq>Mon%NDimX) then

        ir = pq - Mon%NDimX
        fact = Mon%CICoef(ir)
        do i=1,coef*NDimEx
           VecEx((i-1)*coef*NDimEx+pq) = fact * &
                                         EVec((i-1)*coef*NDimEx+pq)
        enddo

    endif
 enddo

 Alpha = 0d0
 do t=1,NDimEx

    if(EVal(t).gt.SmallE.and.EVal(t).lt.BigE) then

       frac = 8d0*EVal(t)/(EVal(t)**2d0+freq)

       do pq=1,NDimEx
          ip = MIndEx(1,pq)
          iq = MIndEx(2,pq)
          do rs=1,NDimEx
          !do rs=1,pq
             ir = MIndEx(1,rs)
             is = MIndEx(2,rs)

             Alpha(pq,rs) = Alpha(pq,rs) + &
                          !EVec(pq,t)*EVec(rs,t)*frac
                          !VecEx((pq-1)*coef*NDimEx+t)*VecEx((rs-1)*coef*NDimEx+t)*frac
                          frac*VecEx((t-1)*coef*NDimEx+pq)*VecEx((t-1)*coef*NDimEx+rs)
             !Alpha(rs,pq) = Alpha(pq,rs)

          enddo
       enddo
    endif
 enddo

print*, 'resp2,Alpha',norm2(Alpha)

deallocate(MIndEx)
deallocate(VecEx)

end subroutine calc_resp_apsg2

subroutine calculateInitialA_diag(Diag,ABPM,NDimX)
! calculate A2 = ABPM - ABPM0
implicit none

integer,intent(in) :: NDimX
double precision,intent(inout) :: Diag(NDimX), ABPM(NDimX,NDimX)

integer :: i

Diag = 0d0
do i=1,NDimX
   Diag(i) = ABPM(i,i)
   ABPM(i,i) = ABPM(i,i) - Diag(i)
enddo

end subroutine calculateInitialA_diag

subroutine calculateInitialA_diagP(AB0,ABPM,Pmat,NDimX)
! calculate A2 = ABPM - ABPM0
implicit none

integer,intent(in) :: NDimX
double precision,intent(in)    :: Pmat(NDimX,NDimX)
double precision,intent(inout) :: AB0(NDimX,NDimX), ABPM(NDimX,NDimX)

integer :: i
double precision,allocatable :: work(:,:)

AB0 = 0d0
do i=1,NDimX
   AB0(i,i)  = ABPM(i,i)
   ABPM(i,i) = ABPM(i,i) - AB0(i,i)
enddo

! PROJECT A2 WITH PMat : work = PMat.A2
allocate(work(NDimX,NDimX))

Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,PMat,NDimX,ABPM,NDimX,0d0,work,NDimX)
AB0 = AB0 + work

print*, 'A2 norm before projection-moje',norm2(ABPM)
ABPM = ABPM - work
print*, 'A2 norm after projection-moje',norm2(ABPM)

deallocate(work)

end subroutine calculateInitialA_diagP

subroutine calculateInitialA_blk(A2,A0Blk,A0BlkIV,nblk,NDimX,abpm0file)
! a) read ABPM0 blocks from file
! b) calculate A2 = ABPM - ABPM0
implicit none

integer,intent(in)  :: NDimX
integer,intent(out) :: nblk
character(*)        :: abpm0file
double precision    :: A2(NDimX,NDimX)
type(EblockData)    :: A0BlkIV
type(EblockData),allocatable :: A0Blk(:)

integer :: i,j,ipos,jpos,iblk,ii
! test
double precision,allocatable :: A0(:,:)

call read_ABPM0Block(A0Blk,A0BlkIV,nblk,abpm0file)

! HERE: replace with call subtr_blk_right(..,true,...)!

!allocate(A0(NDimX,NDimX))
!A0 = 0
do iblk=1,nblk
   associate(B => A0Blk(iblk))

     ! A2 = ABPM - A0PM
     do j=1,B%n
        jpos=B%pos(j)
        do i=1,B%n
           ipos=B%pos(i)
           A2(ipos,jpos) = A2(ipos,jpos) - B%matX(i,j)
           !A0(ipos,jpos) = B%matX(i,j)
        enddo
     enddo

   end associate
enddo
!IV block
associate(B => A0BlkIV)
  do i=1,B%n
     ii = B%pos(i)
     A2(ii,ii) = A2(ii,ii) - B%vec(i)
     !A0(ii,ii) = B%vec(i)
  enddo
end associate

!print*, 'test: A0',norm2(A0)
!deallocate(A0)

end subroutine calculateInitialA_blk

subroutine calculateInitialA_blkP(A2,A0,Pmat,NDimX,abpm0file)
! a) read ABPM0 blocks from file
! b) calculate A2 = ABPM - ABPM0
implicit none

integer,intent(in)  :: NDimX
character(*)        :: abpm0file
double precision    :: Pmat(NDimX,NDimX)
double precision    :: A2(NDimX,NDimX),A0(NDimX,NDimX)

type(EblockData)    :: A0BlkIV
type(EblockData),allocatable :: A0Blk(:)

integer :: i,j,ii,ipos,jpos
integer :: iblk,nblk
double precision,allocatable :: work(:,:)

! A2 = A2 - A0(blk)
call read_ABPM0Block(A0Blk,A0BlkIV,nblk,abpm0file)
call add_blk_right(A2,A0Blk,A0BlkIV,-1d0,.true.,nblk,NDimX)

! PROJECT A2 WITH PMat : A0 = PMat.A2
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,PMat,NDimX,A2,NDimX,0d0,A0,NDimX)

print*, 'A2 norm before projection',norm2(A2)
A2 = A2 - A0
print*, 'A2 norm after projection',norm2(A2)

!A0 = A0 + Pmat.A2
call add_blk_right(A0,A0Blk,A0BlkIV,1d0,.true.,nblk,NDimX)

end subroutine calculateInitialA_blkP

subroutine calculateLambda_Pmat(Lambda,A0,omega,NDimX)
! THIS IS FOR TESTING ONLY!
! calculate Lambda with Pmat
! Lambda = (A0 + omega^2 + A2.Q )-1
!
implicit none

integer,intent(in)  :: NDimX
double precision,intent(in)    :: omega
double precision,intent(in)    :: A0(NDimX,NDimX)
double precision,intent(out)   :: Lambda(NDimX,NDimX)

integer :: i,inf1,inf2
integer,allocatable            :: ipiv(:)
double precision,allocatable   :: work(:)

Lambda = A0
do i=1,NDimX
   Lambda(i,i) = Lambda(i,i) + omega**2
enddo

allocate(ipiv(NDimX),work(NDimX))

call dgetrf(NDimX,NDimX,Lambda,NDimX,ipiv,inf1)
call dgetri(NDimX,Lambda,NDimX,ipiv,work,NDimX,inf2)

deallocate(work,ipiv)

end subroutine calculateLambda_Pmat

subroutine calculateLambda_diag(Lambda,Diag,omega,NDimX)
! calculate Lambda in diagonal form
! Lambda = (A0+omega^2)-1
implicit none

integer,intent(in)  :: NDimX
double precision,intent(in)    :: omega, Diag(NDimX)
double precision,intent(inout) :: Lambda(NDimX)

integer :: i

do i=1,NDimX
   Lambda(i) = 1d0 / (Diag(i) + omega**2)
enddo

end subroutine calculateLambda_diag

subroutine calculateLambda_blk(Lambda,LambdaIV,omega,NDimX,nblk,A0Blk,A0BlkIV)
! calculate Lambda by inversion of blocks ABPM0 blocks
! Lambda = (A0+omega^2)-1
implicit none

integer,intent(in)  :: NDimX,nblk
double precision,intent(in) :: omega
type(EblockData)    :: A0Blk(nblk),A0BlkIV
type(EblockData),allocatable :: Lambda(:)
type(EblockData)             :: LambdaIV

integer :: i,ii,j,ipos,jpos,iblk
integer :: info
integer,allocatable :: ipiv(:)
double precision,allocatable :: work(:)
! for testing full Lambda
!double precision,allocatable :: A0Inv(:,:)

!allocate(A0Inv(NDimX,NDimX))
!A0Inv=0

Lambda = A0Blk
do iblk=1,nblk
   associate(B => Lambda(iblk), A => A0Blk(iblk))

     allocate(ipiv(B%n),work(B%n))

     do i=1,B%n
        B%matX(i,i) = B%matX(i,i) + omega
     enddo

     call dgetrf(B%n,B%n,B%matX,B%n,ipiv,info)
     call dgetri(B%n,B%matX,B%n,ipiv,work,B%n,info)

     !! test Lambda
     !do j=1,B%n
     !   jpos=B%pos(j)
     !   do i=1,B%n
     !      ipos=B%pos(i)
     !      A0Inv(ipos,jpos) = B%matX(i,j)
     !   enddo
     !enddo

     deallocate(work,ipiv)

   end associate
enddo

! IV block
associate(B => A0blkIV, A => LambdaIV)
  A%n = B%n
  A%l1 = B%l1
  A%l2 = B%l2
  allocate(A%pos(A%n),A%vec(A%n))
  A%pos = B%pos
  do i=1,B%n
     ii = B%pos(i)
     A%vec(i) = 1d0 / (B%vec(i) + omega)
     !A0Inv(ii,ii) = 1d0 / (B%vec(i) + omega)
  enddo
end associate

!print*, 'Lambda',omega,norm2(A0Inv)
!deallocate(A0Inv)

end subroutine calculateLambda_blk

subroutine Cmat_blk_iterDIIS(CTilde,NDimX,NCholesky,nblk,Lambda,LambdaIV,ABPTilde,A2,iStats)
use diis
use class_IterStats
implicit none

type(SaptDIIS) :: SAPT_DIIS
type(DIISData) :: DIISBlock

type(EblockData)                :: Lambda(nblk),LambdaIV
class(IterStats), intent(inout) :: iStats
integer,intent(in) :: nblk,NDimX,NCholesky
double precision, intent(in)  :: A2(NDimX,NDimX), ABPTilde(NDimX,NCholesky)
double precision, intent(out) :: CTilde(NDimX,NCholesky)
! for testing full Lambda (passed a work)
!double precision, intent(in)  :: work(NDimX,NDimX)

integer :: N
double precision :: A2CTilde(NDimX,NCholesky), CTilde_prev(NDimX,NCholesky)
double precision :: norm
double precision :: Thresh

!print*, 'SAPTDiis: test'
!print*, 'DIISN',SAPT_DIIS%DIISN
!print*, 'DIISON',SAPT_DIIS%DIISOn
!print*, 'maxIter',SAPT_DIIS%maxIter
!print*, 'Thresh',SAPT_DIIS%Threshold

call init_DIIS(DIISBlock,NDimX*NCholesky,NDimX*NCholesky,SAPT_DIIS%DIISN)

CTilde_prev = 0.0d0
N = 1

do
   call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,A2,NDimX,CTilde,NDimX,0.0d0,A2CTilde,NDimX)
   A2CTilde = ABPTilde - A2CTilde

   ! with full Lambda(=work)
   !call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,work,NDimX,A2CTilde,NDimX,0.0d0,CTilde,NDimX)
   call ABPM_HALFTRAN_GEN_L(A2CTilde,CTilde,0d0,Lambda,LambdaIV,nblk,NDimX,NCholesky,'X')

   if(N > SAPT_DIIS%DIISOn ) then
       call use_DIIS(DIISBlock, CTilde, CTilde - CTilde_prev)
   endif

   norm = norm2(CTilde - CTilde_prev)
   !print '(a, i3, a, e)', "norm (", N, ") = ", norm

   if ((norm < SAPT_DIIS%Thresh) .or. (N .ge. SAPT_DIIS%maxIter .and. SAPT_DIIS%maxIter .ge. 0)) then
       call iStats%addN(N)
       exit
   endif

   CTilde_prev = CTilde
   N = N + 1

enddo

call free_DIIS(DIISBlock)

end subroutine Cmat_blk_iterDIIS

subroutine Cmat_diag_iterDIIS(CTilde,NDimX,NCholesky,Lambda,ABPTilde,A2,iStats)
use diis
use class_IterStats
implicit none

type(SaptDIIS) :: SAPT_DIIS
type(DIISData) :: DIISBlock

class(IterStats), intent(inout) :: iStats
integer,intent(in) :: NDimX,NCholesky
double precision, intent(in)  :: Lambda(NDimX)
double precision, intent(in)  :: A2(NDimX,NDimX), ABPTilde(NDimX,NCholesky)
double precision, intent(out) :: CTilde(NDimX,NCholesky)

integer :: N
double precision :: A2CTilde(NDimX,NCholesky), CTilde_prev(NDimX,NCholesky)
double precision :: norm
double precision :: Thresh

call init_DIIS(DIISBlock,NDimX*NCholesky,NDimX*NCholesky,SAPT_DIIS%DIISN)

CTilde_prev = 0.0d0
N = 1

do
   call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,A2,NDimX,CTilde,NDimX,0.0d0,A2CTilde,NDimX)
   A2CTilde = ABPTilde - A2CTilde

   call MultpDiagMat(Lambda,A2CTilde,0d0,CTilde,NDimX,NCholesky)

   if(N > SAPT_DIIS%DIISOn ) then
       call use_DIIS(DIISBlock, CTilde, CTilde - CTilde_prev)
   endif

   norm = norm2(CTilde - CTilde_prev)
   !print '(a, i3, a, e)', "norm (", N, ") = ", norm

   if ((norm < SAPT_DIIS%Thresh) .or. (N .ge. SAPT_DIIS%maxIter .and. SAPT_DIIS%maxIter .ge. 0)) then
       call iStats%addN(N)
       exit
   endif

   CTilde_prev = CTilde
   N = N + 1

enddo

call free_DIIS(DIISBlock)

end subroutine Cmat_diag_iterDIIS

subroutine C_AlphaExpand(CTilde,C0Tilde,OmI,Max_Cn,A1,A2,ABP0Tilde,ABP1Tilde,&
                         A0Blk,A0BlkIV,nblk,NCholesky,NDimX,ErrMax,ifreq)
!
!  For a given frequency OmI, CTilde(Alpha=1) is computed by expanding
!  around alpha=0 up to the order Max_Cn
!
!  for the original Kasia's version, see, e.g. commit 94fd018ed7efd
!
implicit none

integer,intent(in) :: Max_Cn
integer,intent(in) :: nblk,NCholesky,NDimX
double precision,intent(in)  :: OmI

double precision,intent(in)    :: A1(NDimX,NDimX),A2(NDimX,NDimX)
double precision,intent(in)    :: ABP0Tilde(NDimX,NCholesky),ABP1TIlde(NDimX,NCholesky)
double precision,intent(inout) :: CTilde(NDimX,NCholesky),C0Tilde(NDimX,NCholesky)

integer :: N
double precision :: XFactorial, XN1, XN2
! testing only

type(EBlockData)             :: A0Blk(nblk),A0BlkIV

type(EBlockData),allocatable :: Lambda(:)
type(EBlockData)             :: LambdaIV

double precision, allocatable :: C1Tilde(:,:),C2Tilde(:,:)
double precision, allocatable :: Work(:,:)
! testing ErrMax
integer :: ifreq
double precision :: WFact,XNorm0,XNorm1,PI,FF
double precision :: ErrMax

PI = 4.0*ATAN(1.0)

allocate(C1Tilde(NDimX,NCholesky),C2Tilde(NDimX,NCholesky),Work(NDimX,NCholesky))

! Calc: LAMBDA=(A0+Om^2)^-1
call calculateLambda_blk(Lambda,LambdaIV,OmI**2,NDimX,nblk,A0Blk,A0BlkIV)

! Calc: C0Tilde=1/2 LAMBDA.APLUS0Tilde
call ABPM_HALFTRAN_GEN_L(ABP0Tilde,C0Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NCholesky,'X')
C0Tilde = 0.5d0*C0Tilde
!print*, 'C0Tilde',OmI,norm2(C0Tilde)

! Calc: C1Tilde=LAMBDA.(1/2 APLUS1Tilde - A1.C0Tilde)
call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,A1,NDimX,C0Tilde,NDimX,0.0d0,CTilde,NDimX)
CTilde = 0.5d0*ABP1Tilde - CTilde
call ABPM_HALFTRAN_GEN_L(CTilde,C1Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NCholesky,'X')
!print*, 'C1Tilde',OmI,norm2(C1Tilde)

! test uncoupled
CTilde = C0Tilde
! test semicoupled
CTilde = CTilde + C1Tilde

WFact = 4.D0/PI*OmI
XFactorial = 1d0
do N=2,Max_Cn

   XFactorial = XFactorial*N
   XN1 = -N
   XN2 = -N*(N-1)

   call dgemm('N','N',NDimX,NCholesky,NDimX,XN2,A2,NDimX,C0Tilde,NDimX,0.0d0,Work,NDimX)
   call dgemm('N','N',NDimX,NCholesky,NDimX,XN1,A1,NDimX,C1Tilde,NDimX,1.0d0,Work,NDimX)
   call ABPM_HALFTRAN_GEN_L(Work,C2Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NCholesky,'X')

   FF = WFact / XFactorial / (N+1)
   XNorm1 = norm2(FF*C2Tilde)
   write(6,'(X,"Order (n), |Delta_C|",I3,E14.4)') N, XNorm1
   if(XNorm1.Lt.ErrMax) exit

   CTilde  = CTilde + C2Tilde / XFactorial
   C0Tilde = C1Tilde
   C1Tilde = C2Tilde

   !print*, 'N,COMTilde',N,norm2(CTilde)

enddo

call release_ac0block(Lambda,LambdaIV,nblk)

if(ifreq==1) ErrMax = norm2(FF*C2Tilde)

deallocate(Work,C2Tilde,C1Tilde)

end subroutine C_AlphaExpand

subroutine C_AlphaExpand_unc(C0Tilde,OmI,A1,A2,ABP0Tilde,ABP1Tilde,&
                         A0Blk,A0BlkIV,nblk,NCholesky,NDimX)
!
!  For a given frequency OmI, CTilde(Alpha=1) is computed by expanding
!  around alpha=0, here only in 0-th term, i.e., uncoupled
!
implicit none

integer,intent(in)           :: nblk,NCholesky,NDimX
double precision,intent(in)  :: OmI

double precision,intent(in)    :: A1(NDimX,NDimX),A2(NDimX,NDimX)
double precision,intent(in)    :: ABP0Tilde(NDimX,NCholesky),ABP1TIlde(NDimX,NCholesky)
double precision,intent(inout) :: C0Tilde(NDimX,NCholesky)

integer :: N
double precision :: XFactorial, XN1, XN2

type(EBlockData)             :: A0Blk(nblk),A0BlkIV

type(EBlockData),allocatable :: Lambda(:)
type(EBlockData)             :: LambdaIV

! Calc: LAMBDA=(A0+Om^2)^-1
call calculateLambda_blk(Lambda,LambdaIV,OmI**2,NDimX,nblk,A0Blk,A0BlkIV)

! Calc: C0Tilde=1/2 LAMBDA.APLUS0Tilde
call ABPM_HALFTRAN_GEN_L(ABP0Tilde,C0Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NCholesky,'X')
C0Tilde = 0.5d0*C0Tilde

call release_ac0block(Lambda,LambdaIV,nblk)

end subroutine C_AlphaExpand_unc

subroutine prepare_resp_Cmat(Mon,A1,A2,ABP0Tilde,ABP1TIlde,DChol,NDimX,NCholesky,NBasis)
!
! prepare matrices for construction of CAlpha(omega)
!
implicit none

type(SystemBlock) :: Mon
integer,intent(in) :: NDimX,NCholesky,NBasis

integer :: iunit
integer :: nblk
integer :: i,j,ip,iq,ipq

double precision,intent(out) :: DChol(NCholesky,NDimX)
double precision,intent(out) :: A1(NDimX,NDimX), &
                                A2(NDimX,NDimX), &
                                ABP0Tilde(NDimX,NCholesky), &
                                ABP1Tilde(NDimX,NCholesky)

double precision :: Cpq
character(:),allocatable     :: abfile,xy0file
double precision,allocatable :: ABPLUS1(:,:),ABMIN1(:,:)

type(EBlockData)             :: A0BlkIV
type(EBlockData),allocatable :: A0Blk(:)

! set filenames
if(Mon%Monomer==1) then
   abfile='ABMAT_A'
   xy0file='XY0_A'
elseif(Mon%Monomer==2) then
   abfile='ABMAT_B'
   xy0file='XY0_B'
endif

DChol = 0
do j=1,NDimX
   ip = Mon%IndN(1,j)
   iq = Mon%IndN(2,j)
   ipq = iq + (ip-1)*NBasis
   Cpq = Mon%CICoef(ip) + Mon%CICoef(iq)
   do i=1,NCholesky
      DChol(i,j) = Cpq*Mon%FF(i,ipq)
   enddo
enddo

allocate(ABPLUS1(NDimX,NDimX),ABMIN1(NDimX,NDimX))

open(newunit=iunit,file=abfile,status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1
read(iunit) ABMIN1

close(iunit)

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0Blk,A0BlkIV,Mon%IndN,Mon%CICoef,nblk,NBasis,NDimX,xy0file)

! AB1 = AB1 - A0
call add_blk_right(ABPLUS1,A0Blk,A0BlkIV,-1d0,.false.,nblk,NDimX)
call add_blk_right(ABMIN1, A0Blk,A0BlkIV,-1d0,.true., nblk,NDimX)

print*, 'ABPLUS1',norm2(ABPLUS1)
print*, 'ABMIN1 ',norm2(ABMIN1)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call ABPM_HALFTRAN_GEN_L(ABMIN1, A1,0.0d0,A0Blk,A0BlkIV,nblk,NDimX,NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1,A1,1.0d0,A0Blk,A0BlkIV,nblk,NDimX,NDimX,'X')
print*, 'A1',norm2(A1)

!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,ABMIN1,NDimX,0.0d0,A2,NDimX)
print*, 'A2',norm2(A2)

!Calc: APLUS0Tilde=ABPLUS0.DChol
call ABPM_HALFTRAN_GEN_L(transpose(DChol),ABP0Tilde,0.0d0,A0Blk,A0BlkIV,nblk,NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde',norm2(ABP0Tilde)

!Calc: APLUS1Tilde=ABPLUS1.DChol
Call dgemm('N','T',NDimX,NCholesky,NDimX,1d0,ABPLUS1,NDimX,DChol,NCholesky,0.0d0,ABP1Tilde,NDimX)
print*, 'APLUS1Tilde',norm2(ABP1Tilde)

call release_ac0block(A0Blk,A0BlkIV,nblk)

deallocate(A0Blk)
deallocate(ABPLUS1,ABMIN1)

end subroutine prepare_resp_Cmat

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

subroutine print_en(string,val,nline)
implicit none

character(*),intent(in)     :: string
double precision,intent(in) :: val
logical,intent(in)          :: nline

if(nline) then
   write(LOUT,'(/1x,a,t19,a,f16.8)') string, "=", val
else
   write(LOUT,'(1x,a,t19,a,f16.8)')  string, "=", val
endif

end subroutine print_en

end module sapt_utils
