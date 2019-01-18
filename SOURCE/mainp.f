C     A GENERAL CODE FOR CARRYING OUT LR-DMFT/SR-DFT CALCULATIONS
C
C                 PROGRAM INTERFACED WITH OUTPUT FILES FROM MOLPRO DEVELOP VERSION
C                 ONE- AND TWO-ELECTRON INTEGRALS (IN HF MO) ARE OBTAINED
C                 FROM ATMOL AND THEY ARE USED AS A GUESS WITH U=1
C                 GAUSSIAN BASIS SET USED
C
C
C     K.PERNAL 2018 
C
      Program PRDMFT
C
      use types
      use inputfill
      use systemdef
      use sapt_main
      use timing
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FMultTab
C
      Real*8, Dimension(:), Allocatable :: Occ
      Real*8, Dimension(:), Allocatable :: URe
      Real*8, Dimension(:), Allocatable :: OrbGrid
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: XKin
      Real*8, Dimension(:), Allocatable :: XNuc
      Real*8, Dimension(:), Allocatable :: TwoEl
      Real*8, Dimension(:), Allocatable :: TwoElErf
      Real*8, Dimension(:), Allocatable :: DipX
      Real*8, Dimension(:), Allocatable :: DipY
      Real*8, Dimension(:), Allocatable :: DipZ
      Real*8, Dimension(:), Allocatable :: EpsHF 
      Real*8, Dimension(:), Allocatable :: UMOAO
      Integer, Dimension(:), Allocatable :: NSymMO
C
      type(InputData) :: Input
      type(FlagsData) :: Flags
      type(SystemBlock) :: System
      type(SaptData) :: Sapt
C
      Include 'commons.inc'
C     
C     COMMON BLOCK USED IN OPTCPMFT
C     
      Common/CPMFT/ MFrac,MOcc,NFrac
C     
C     *************************************************************************
C
      Call read_Input(Input)
      Call check_Calc(Input%CalcParams)
      Call fill_Flags(Input,Flags)
      Call create_System(Input,Flags,System,Sapt)
C
      Call free_Input(Input)
C    
C     FILL COMMONS AND CONSTANTS
      XELE = System%XELE
      NELE = System%NELE
      Charge = System%Charge
      NBasis = System%NBasis
      Title = Flags%JobTitle
      ITwoEl = Flags%ITwoEl
      IWarn = 0 
C
C     *************************************************************************
C
C     ICASSCF = 1 - CORRELATION ENE FOR CAS WILL BE COMPUTED
C             = 0 - CORRELATION ENE FOR GVB WILL BE COMPUTED
C
      ICASSCF=Flags%ICASSCF

C     IF IDALTON=1 - READ INTEGRALS AND C COEF, IGem FROM A FILE GENERATED BY DALTON
C
      IDALTON=Flags%IDALTON
C
C     IF IRes=1 - RESTART THE CALCULATIONS FROM A RESTART FILE
C
c herer!!!
      IRes=Flags%IRes
C
C     SET IAO TO ONE IF INTEGRALS IN AO ARE OBTAINED FROM MOLPRO
C
      IAO=Flags%IAO
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
      IFlAC=Flags%IFlAC
C
C     IFlSnd  = 1 - run AC0 (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0
      IFlSnd=Flags%IFlSnd
C
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation 
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
      IFlCore=Flags%IFlCore
C
      IFlFrag1=Flags%IFlFrag
      IFl12=Flags%IFl12
C
C     NoSym = 1 - DO NOT IMPOSE SYMMETRY (RUN molpro WITH 'nosym')
C           = 0 - IMPOSE SYMMETRY (RUN molpro WITHOUT 'nosym')
C
      NoSym=Flags%NoSym
C
C     ************************************************************************* 
C
C     SELECT A LONG-RANGE DMFT FUNCTIONAL
C
C     IFun   = 0 - NO DMFT FUNCTIONAL (PURE DFT CALCULATIONS, EQUIVALENT TO USING mu=0)
C
C     IFun   = 1 - USE THE KUTZELNIGG FUNCTIONAL:
C                  Eee = sqrt(ni nj) <ij|ji> + sqrt(na nb) <ab|ba>
C                      - 2 sqrt(ni na) <ia|ai>
C
C     IFun   = 2 - USE THE BUIJSE-BAERENDS FUNCTIONAL:
C                  Eee = np nq <pq|pq> - sqrt(np nq) <pq|qp>
C
C     IFun   =20 - BB WITH A "+" SIGN FOR THE BONDING-ANTIBONDING PAIR 
C
C     IFun   = 3 - USE THE GOEDECKER-UMRIGAR FUNCTIONAL:
C                  Eee = np nq <pq|pq> - sqrt(np nq) <pq|qp>
C                      + (np-np**2) <pp|pp>
C
C     IFun   = 4 - USE THE BBC1 FUNCTIONAL
C
C     IFun   = 5 - USE THE BBC2 FUNCTIONAL
C
C     IFun   = 6 - USE THE BBC3 FUNCTIONAL
C
C     IFun   = 7 - CORRECTED HARTREE-FOCK FUNCTIONAL
C                  Eee = np nq <pq|pq>
C                      - [(np nq) + sqrt(np(1-np)nq(1-nq)) ]<pq|qp>
C     IFun   = 8 - BB AND HF HYBRID FUNCTIONAL
C                  (1-Cmix)*Eee_BB + Cmix*Eee_HF
C     IFun   = 10 - LR-BB+SR-HF+SR-PBE_CORR
C                  Eee_LRBB + Eee_SRPBE
C     IFun   = 11 - CPMFT method of Scuseria et al. 
C                   (define the active space: NActive, MActive 
C                    - the number of active elecrons and orbitals)
C
C     IFun   = 12 - Hartree-Fock
C
C     IFun   = 13 - APSG
C
C     IGVB   = 1  - APSG WITH ONLY TWO ORBITALS PER GEMINAL
c herer!!!
      IGVB=Flags%IGVB
C
      IFun=Flags%IFun
C
C     *************************************************************************
C      
C     SELECT A SHORT-RANGE DFT FUNCTIONAL 
C
C     IFunSR = 0 - do not include a short-range functional
C     IFunSR = 1 - SR-LSDA, Paziani et al. 
C     IFunSR = 2 - SR-PBE, Goll et al. PCCP 7, (2005) 3917
C     IFunSR = 3 - Gagliardi-Truhlar with PBE
C     IFunSR = 4 - LR-CAS+SR-PBE+LR-AC/AC0 (with full range CAS RDM's)
C     IFunSR = 5 - CASPiDFT
C
      IFunSR=Flags%IFunSR
      IFunSRKer=Flags%IFunSRKer
C
      If(IFunSRKer.Eq.2)Stop'Fatal Error: Kernel for PBE not available!'
      Write(*,'(/,"IFunSR=",I1)')IFunSR
      Write(*,'(/,"IFunSRKer=",I1,/)')IFunSRKer
C
C     *************************************************************************
C
C     NBasis READ FROM SIRIUS.RST
      If(IDALTON.Eq.1) Call basinfo(NBasis,'SIRIUS.RST','DALTON')
C
C     CHECK IF SAPT RUN      
      ISAPT=Flags%ISAPT
C
      If(ISAPT.Eq.1) Call sapt_driver(Flags,Sapt)
C
C     *************************************************************************
C         
C     SELECT ELECTRONIC STATE  
C
      NStates = System%NStates
      InSt(1:2,1:NStates) = System%InSt
C
      If(InSt(2,1).Gt.0) Then
      NoSt = System%InSt(1,1)
C      Print*,'VALUE DECLARED IN INPUT: ',NoSt
C
      ElseIf(IDALTON.Eq.0) Then
      Call read_NoSt_molpro(NoSt,'2RDM')
      ElseIf(IDALTON.Eq.1) Then
      NoSt = 1
      Write(6,'(/,1x,a)') 'WARNING! ASSUMING RMDs CORRESPOND TO
     $                     GROUND STATE FROM DALTON OUTPUT (NoSt=1)!'
C
      EndIf
C
      Write(6,'(1x,"NoSt=",I2)')NoSt
C
C     READ THE INPUT AND PRINT THE INPUT DATA 
C
C     OLD INPUT-READ
C      Call RWInput(Title,ZNucl,Charge,NBasis)
C
C
C     CALCULATE THE DIMENSIONS
      If(IDALTON.Eq.0) then
C        Call CheckNBa(NBasis,Title)
        Call basinfo(NBasis,'AOTWOINT.mol','MOLPRO')
      endif
C
      Call DimSym(NBasis,NInte1,NInte2,MxHVec,MaxXV)
C
      If(IFunSR.Ne.0) Then
C      Call GetNGrid(NGrid,Title)
      Call molprogrid0(NGrid,NBasis)
      Else
      NGrid=1
      EndIf
C
C     GET THE VALUE OF THE SEPARATION PARAMETER OM
C
      Alpha = System%Omega
      If(IFunSR.Ne.0.And.IFunSR.Ne.3) Then
C      Call GetAlpha(Title)
      Call readalphamolpro(Alpha)
      Else
      Alpha=0.D0
      EndIf
C
C     ALLOCATE THE MATRICES
C
      Allocate  (Occ(NBasis))
      Allocate  (URe(NBasis*NBasis))
      Allocate  (XKin(NInte1))
      Allocate  (XNuc(NInte1))
      Allocate  (TwoEl(NInte2))
      Allocate  (UMOAO(NBasis*NBasis))
C
      Occ(1:NBasis)=           0.D0
      URe(1:NBasis*NBasis)=    0.D0
      XKin(1:NInte1)=          0.D0
      XNuc(1:NInte1)=          0.D0
      TwoEl(1:NInte2)=         0.D0
      UMOAO(1:NBasis*NBasis)=  0.D0
C
C     LOAD THE GRID AND THE VALUES OF THE MO ORBITALS ON THE GRID
C
      write(LOUT,'()')
      write(LOUT,'(1x,a)') 'STARTING CALCULATIONS '
      write(LOUT,'(8a10)') ('**********',i=1,8)
C
      Call clock('START',Tcpu,Twall)
C     LOAD THE INTEGRALS
C
      If(IDALTON.Eq.0) Then
C
      Call LdInteg(Title,XKin,XNuc,ENuc,Occ,URe,TwoEl,UMOAO,NInte1,
     $ NBasis,NInte2,NGem)
C
      Else
C
      Call ReadDAL(XKin,XNuc,ENuc,Occ,URe,TwoEl,UMOAO,
     $ NInte1,NBasis,NInte2,NGem,Flags)     
C
      EndIf  
C
C      If(IFunSR.Eq.5) Then
C     $ Call CASPIDFT(ENuc,URe,UMOAO,Occ,XKin,TwoEl,
C     $ NBasis,NInte1,NInte2)
C      Else
      Call DMSCF(Title,URe,Occ,XKin,XNuc,ENuc,UMOAO,
     $ TwoEl,NBasis,NInte1,NInte2,NGem,System)
C      EndIf
C
      If(IWarn.Gt.0) Then
      Write(6,'(/,1x,a,i2,1x,a)') 'CHECK OUTPUT FOR',IWarn,'WARNINGS!'
C      Write(6,'(8a10)') ('**********',i=1,9)
      EndIf
C
      Call free_System(System)
      Call clock(PossibleJobType(Flags%JobType),Tcpu,Twall)
C      Stop
      End
C
