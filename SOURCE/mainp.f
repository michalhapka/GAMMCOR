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
      Call free_Input(Input)
C    
C     FILL SYSTEM COMMONS
      XELE = System%XELE
      NELE = System%NELE
      NBasis = System%NBasis
C     HERE!!!!
      write(LOUT,*) 'XELE-NELE-NBasis'
      write(LOUT,*) XELE, NELE, NBasis 
C
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
C     SET INO TO 1 IF YOU WANT TO USE NATURAL OCCUPATION NUMBERS AND  
C     THE NATURAL ORBITALS FROM THE MOLPRO CALCULATION 
C
C     FOR TD CALCS WHEN THE SYMMETRY IS ON ALWAYS USE INO=0 !!!
C
      INO=Flags%INO
      Write(*,'(/,"INO=",I1,/)')INO
C
c herer!!!
C     NoSym = 1 - DO NOT IMPOSE SYMMETRY (RUN molpro WITH 'nosym')
C           = 0 - IMPOSE SYMMETRY (RUN molpro WITHOUT 'nosym')
C
      NoSym=Flags%NoSym
C
C
C     MxSym = NUMBER OF IRREPS IN A GIVEN GROUP
C
      If(NoSym.Eq.1) Then
C
      MxSym=1
      FMultTab="c1_table.txt"
C
      Else  
C
c herer!!!
C C1
c       MxSym=1
c       FMultTab="c1_table.txt"  
C Cs point group
c       MxSym=2
c       FMultTab="cs_table.txt"
C HeH+,LiH,BH,H2O,CH2O, BeH2
      MxSym=4
      FMultTab="c2v_table.txt"
C H2,Be,N2,Li2, H4
c      MxSym=8
c      FMultTab="d2h_table.txt"
C
      EndIf
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
      Cmix=0.2D0 
C
C     If IFun=6 set which orbitals should be treated as bonding
C
      Do I=1,1000
      IType(I)=1
      EndDo
      If(IFun.Eq.6) Then
c     it is decided 'on the fly' which orbitals are bonding
c     the criterium is based on the value of occ, see Deck GOCC
c      IType(1)=0
c      IType(2)=0
      EndIf
C
C     *************************************************************************
C      
C     SELECT A SHORT-RANGE DFT FUNCTIONAL 
C
C     IFunSR = 0 - do not include a short-range functional
C     IFunSR = 1 - SR-LSDA, Paziani et al. 
C     IFunSR = 2 - SR-PBE, Goll et al. PCCP 7, (2005) 3917
C
      IFunSR=Flags%IFunSR
      IFunSRKer=Flags%IFunSRKer
C
      If(IFun.Eq.10.And.IFunSR.Eq.1) 
     $ Stop 'Fatal Error: conflicting IFun and IFunSR'
      If(IFunSRKer.Eq.2)Stop'Fatal Error: Kernel for PBE not available!'
      Write(*,'(/,"IFunSR=",I1)')IFunSR
      Write(*,'(/,"IFunSRKer=",I1,/)')IFunSRKer
C
C     READ THE INPUT AND PRINT THE INPUT DATA 
C

      Call RWInput(Title,ZNucl,Charge,NBasis)
C
C     *************************************************************************
C
C
C     CHOOSE THE OCCUPATION PATTERN FOR GEMINALS
C
C     IModG = 0 - THE OCCUPATION PATTERN IS SET IN DMSCF AND FIXED. THE PATTERN IS 
C                 CHOSEN AS: GEM FROM 1 TO NGOcc FULLY OCCUPIED (CONSISTING OF TWO EL)
C                 THEN THE ORBITALS ARE ASSIGNED CONSECUTIVELY TO GEMINALS STARTING FROM NGOcc+1
C             1 - OCCUPANCIES OF GEMINALS FROM 1 TO NGOcc IS FIXED TO 1, THE ARAI SPACES OF THE
C                 REST IS OPTIMIZED
C
C     IF NGOcc DIFFERENT FROM ZERO - SOME ORBITALS WILL BE FIXED TO BE FULLY OCCUPIED
C     IF IFreeze = 1 THOSE ORBITALS WILL BE FROZEN, OTHERWISE THEY WILL BE VARIED
C
      IModG=1
C
C     SET THE NUMBER OF BONDS TO ESTABLISH NGOcc
c herer!!!
C
c      NGOcc = NELE
      NGOcc=0
C
C     ILoc = 1 - READ UMOAO MATRIX FROM THE molden FILE WITH LOCALIZED MO'S
C
      ILoc=1
C
C     IFreeze = 1 - WORKS ONLY WITH APSG. FREEZE ORBITALS THAT ARE FULLY OCCUPIED AND 
C                   CONSTITUTE FULLY OCCUPIED GEMINALS (REMEMBER TO MODIFY NGOcc ACCORDINGLY)
c herer!!!
      IFreeze=0
C
C     If ILoc=1 the core orbitals are not put as the first few ones by molrpo. do not use IFreeze=1 (you risk that other than core orbitals will be frozen)
C
      If(ILoc.Eq.1.And.IFreeze.Eq.1.Or.ILoc.Eq.1.And.NGOcc.Gt.0) 
     $ Write(*,*) 'If using IFreeze=1 together with ILoc=1 
     $ or NGOcc>0! make sure that orbitals are ordered so that
     $ core orbitals are the first one'
C
      If (IFreeze.Eq.0) 
     $ Write(*,'(/,X,"IFreeze = ",I1," all orbitals varied",/)')IFreeze
      If (IFreeze.Eq.1) 
     $ Write(*,'(/,X,"IFreeze = ",I1," some orbitals are frozen",/)')
     $ IFreeze
C
c herer!!!
      Small=1.D-8
C
C     SELECT A METHOD FOR ERPA CALCULATIONS
C
C     IAPSG =  1 - ERPA WITH THE APSG 2-RDM (Call APSG_NEST)
C              0 - LINEARLIZED ERPA (Call LinAB_NEST)
C     ISERPA = 1 - EXCITATIONS FROM SERPA (Call SERPAExcit)
C     ISERPA = 0 - EXCITATIONS FROM ERPA/APSG OR LINEARIZED ERPA (RPAExcit)     
C     ISERPA = 2 - EXCITATIONS FROM LINEAR RESPONSE (PINO) WITH APSG (PINOExcit)
C
      IAPSG=1
      ISERPA=0
      If(IAPSG.Eq.0.And.ISERPA.Eq.1)
     $ Stop 'Fatal Error: Conflicting IAPSG and ISERPA!'
C
C     SET QMAX FOR ERPA CALCULATIONS
C
c herer!!!
      QMAX=NBasis
c      QMAX=16
C
C     *************************************************************************
C     SET THE NUMBER OF INACTIVE ORBITALS, ACTIVE ORBITALS, ADN VIRTUAL ORBITALS
C     IN ORBITAL OPTIMIZATION OptAPSGFun ACT-ACT AND VIRT-VIRT BLOCKS WILL BE DISCARDED
C
c herer!!!
C optimize all
      NInActOrb=0
      NActOrb=NBasis
      NVirtOrb=NBasis-NActOrb-NInActOrb
c
c HF-like partitioning
c      NInActOrb=NELE
c      NActOrb=0
c      NVirtOrb=NBasis-NActOrb-NInActOrb
c
c GVB-like
c      NInActOrb=0
c      NActOrb=NELE*2
c      NVirtOrb=NBasis-NActOrb-NInActOrb
c      If(IFreeze.Eq.1) Then
c      NActOrb=NActOrb-NGOcc
c      NVirtOrb=NVirtOrb+NGOcc
c      EndIf
C
C     *************************************************************************
C
C     SET THE NUMBER OF ACTIVE ELECTRONS (NActive) AND ORBITALS (MActive)
C     TO BE USED IN SELECTING INTERACTING PAIRS IN TIME-DEPENDENT CALCULATIONS
C
      NActive=6
      MActive=6
      MFracTD=MActive
      NFrac=NActive/2
      MOccTD=NELE-NFrac
C
C     SET THE NUMBER OF ACTIVE ELECTRONS (NActive) AND ORBITALS (MActive) 
C     TO BE USED IN POWELL OPTIMIZATION OF OCCUPANCIES FOR IFun=11
C
      If(IFun.Eq.11) Then
      IFun=7
c N2
      NActive=6
      MActive=6
c H2
c      NActive=2
c      MActive=2
      Else
      NActive=NELE*2
      MActive=NBasis
      EndIf
C
      If(NActive/2.Gt.NELE) Stop'Fatal Error: Too many active electrons'
      MFrac=MActive
      NFrac=NActive/2
      MOcc=NELE-NFrac
C
C     CHECK IF SAPT RUN      
      ISAPT=Flags%ISAPT
C
C     *************************************************************************
C
C     CALCULATE THE DIMENSIONS
C
      Call DimSym(NBasis,NInte1,NInte2,MxHVec,MaxXV)
      If(IDALTON.Eq.0) Call CheckNBa(NBasis,Title)
      If(IFunSR.Ne.0) Then
      Call GetNGrid(NGrid,Title)
      Else
      NGrid=1
      EndIf
C
C     GET THE VALUE OF THE SEPARATION PARAMETER OM
C
      If(IFun.Ne.0.And.IFunSR.Ne.0) Then
      Call GetAlpha(Title)
      Else
      Alpha=0.D0
      EndIf
C
C     ALLOCATE THE MATRICES
C
      Allocate  (Occ(NBasis))
      Allocate  (URe(NBasis*NBasis))
      Allocate  (OrbGrid(NBasis*NGrid))
      Allocate  (WGrid(NGrid))
      Allocate  (XKin(NInte1))
      Allocate  (XNuc(NInte1))
      Allocate  (TwoEl(NInte2))
      If(IFunSR.Eq.1) Then
      Allocate  (TwoElErf(NInte2))
      Else
      Allocate  (TwoElErf(1))
      EndIf
      Allocate  (DipX(NInte1))
      Allocate  (DipY(NInte1))
      Allocate  (DipZ(NInte1))
      Allocate  (EpsHF(NBasis))
      Allocate  (UMOAO(NBasis*NBasis)) 
      Allocate  (NSymMO(NBasis))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
C
      Occ(1:NBasis)=           0.D0
      URe(1:NBasis*NBasis)=    0.D0
      OrbGrid(1:NBasis*NGrid)= 0.D0 
      WGrid(1:NGrid)=          0.D0
      XKin(1:NInte1)=          0.D0
      XNuc(1:NInte1)=          0.D0
      TwoEl(1:NInte2)=         0.D0
      If(IFunSR.Eq.1) Then
      TwoElErf(1:NInte2)=      0.D0
      Else
      TwoElErf(1)=      0.D0
      EndIf
      DipX(1:NInte1)=          0.D0
      DipY(1:NInte1)=          0.D0
      DipZ(1:NInte1)=          0.D0 
      EpsHF(1:NBasis)=         0.D0
      UMOAO(1:NBasis*NBasis)=  0.D0
      NSymMO(1:NBasis)=        0
      OrbXGrid(1:NBasis*NGrid)= 0.D0
      OrbYGrid(1:NBasis*NGrid)= 0.D0
      OrbZGrid(1:NBasis*NGrid)= 0.D0
C
C     LOAD THE GRID AND THE VALUES OF THE MO ORBITALS ON THE GRID
C
      If(IFunSR.Ne.0) Then
C
      Call GetGrid(OrbGrid,WGrid,NSymMO,NGrid,Title,NBasis)
C
      Else
C
      If(NoSym.Eq.0) Then
C
C     READ SYMMETRIES OF MO FROM A FILE
C
      Open(10,File="sym_mo.dat")
C
      Counter=0
      JEnd=0
      Do I=1,MxSym
      Counter=Counter+1
      Read(10,*) JMax
     
      Do J=1,JMax
      NSymMO(J+JEnd)=Counter
      EndDo
      JEnd=JEnd+JMax
      EndDo
      If(JEnd.Ne.NBasis) Stop 'Fatal Error. Check sym_mo.dat file'
C
      Else 
C
      Do I=1,NBasis
      NSymMO(I)=1
      EndDo
C
      EndIf
C
      EndIf
C
      If(IFunSR.Eq.2)
     $ Call GetGrad(OrbXGrid,OrbYGrid,OrbZGrid,NGrid,Title,NBasis)
C
C     LOAD THE INTEGRALS
C
      If(IDALTON.Eq.0) Then
C
      Call LdInteg(Title,XKin,XNuc,ENuc,Occ,URe,DipX,DipY,DipZ,
     $ TwoEl,TwoElErf,UMOAO,NSymMO,NInte1,NBasis,NInte2)
C
C     RUN HARTREE-FOCK TO PRODUCE A GUESS (epsilons are needed for guess)
C
      If(IRes.Eq.0.And.INO.Eq.0) Then
C
C     HF GUESS
C
      IJ=0
      Do I=1,NBasis
      Do J=1,NBasis
      IJ=IJ+1
      URe(IJ)=0.D0
      If(I.Eq.J) URe(IJ)=1.0D0
      EndDo
      EndDo
C
      Open(10,File='occ_hf.dat')
      Do I=1,NBasis
      Read(10,*) A,B,C,D
      Occ(I)=B/2.D0
      EndDo
      Close(10)
C
      ElseIf(INO.Eq.1) Then
C
C     LOAD NO's AND NON's FROM THE MOLPRO OUTPUT
C
      Call GetNO(URe,Occ,UMOAO,Title,NBasis) 
C
      EndIf
C
c     If(IDALTON.Eq.0)
      Else
C
      Call ReadDAL(XKin,XNuc,ENuc,Occ,URe,TwoEl,UMOAO,
     $ NInte1,NBasis,NInte2,NGem,Flags)     
C
      EndIf  
C
C     RUN DMFT CALCULATION
C
      Call DMSCF(ETot,Title,URe,Occ,EpsHF,XKin,XNuc,ENuc,UMOAO,
     $ DipX,DipY,DipZ,TwoEl,TwoElErf,NBasis,NInte1,NInte2,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NSymMO,NGrid,
     $ IAPSG,ISERPA,QMAX,IModG,NGOcc,Small,NGem,
     $ Flags)
C
      Stop
      End
C
