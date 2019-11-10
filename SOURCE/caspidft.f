*Deck CASPIDFT
      Subroutine CASPIDFT(ENuc,URe,UNOAO,Occ,XOne,TwoNO,
     $ NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
c      Real*8, Dimension(:), Allocatable :: OrbGrid(:,:)
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RhoGrid
      Real*8, Dimension(:), Allocatable :: Sigma
      Real*8, Dimension(:), Allocatable :: Zk,Zk1
      Real*8, Dimension(:), Allocatable :: rhoo,sigmaco,sigmaoo,vrhoc,
     $ vrhoo,vsigmacc,vsigmaco,vsigmaoo
      Real*8, Allocatable :: RDM2Act(:),OrbGrid(:,:)
      Real*8, Allocatable :: RR(:,:)
      Real*8, Allocatable :: OnTop(:)
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),XOne(NInte1),TwoNO(NInte2),
     $ ConCorr(200)
     $ ,EpsC(NBasis,NBasis),EpsCI(NBasis),EpsDiag(NBasis)
      logical fderiv,open
c      double precision rhoo(ngrid)
c      double precision sigmaco(ngrid),sigmaoo(ngrid)
      integer igrad
      character*(30) name
c      double precision vrhoc(ngrid),vrhoo(ngrid)
c      double precision vsigmacc(ngrid),vsigmaco(ngrid),vsigmaoo(ngrid)
C
      FDeriv=.True.
      Open=.False.
c      Alpha=0.5
      Alpha=0.4
      CMix=0.5
      Write(*,'(/,X,"ALPHA CMix",X,2F5.2)')Alpha,CMix
C
      Write(6,'(/,X,"***************** CASPIDFT ***************** ")')
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/,X,"The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate  (rhoo(NGrid))
      Allocate  (sigmaco(NGrid))
      Allocate  (sigmaoo(NGrid))
      Allocate  (vrhoc(NGrid))
      Allocate  (vrhoo(NGrid))
      Allocate  (vsigmacc(NGrid)) 
      Allocate  (vsigmaco(NGrid)) 
      Allocate  (vsigmaoo(NGrid))
C
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NGrid,NBasis))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
      Allocate  (RhoGrid(NGrid))
      Allocate  (Sigma(NGrid))
      Allocate  (Zk(NGrid))
      Allocate  (Zk1(NGrid))
      Allocate  (RR(3,NGrid))
      Allocate  (OnTop(NGrid))
C
C      Call molprogrid1(RR,NGrid)
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind1(I)=INActive+I
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
C
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
      GoTo 10
   40 Continue
      Close(10)

C     COMPUTE THE ENERGY FOR CHECKING
C
      ETot=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      ETot=ETot+Two*Occ(I)*XOne(II)
      EndDo
      Write(6,'(/,1X,''One-Electron Energy'',6X,F15.8)')ETot 
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      ETot=ETot+FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoNO(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(1X,''CASSCF Energy (w/o ENuc)'',X,F15.8)')ETot
      Write(6,'(1X,''Total CASSCF Energy '',5X,F15.8)')ETot+ENuc
C      
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)   
C
      Do I=1,NGrid
C
      vrhoc(I)=Zero
      vsigmacc(i)=Zero 
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
C
      OnTop(I)=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop(I)=OnTop(I)
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo 
C
      EndDo
C
      Call LYP(RhoGrid,Sigma,Zk,NGrid)
c      Call PBECor(RhoGrid,Sigma,Zk,NGrid)
c      Call GGA_SPIN(Zk1,URe,Occ,
c     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C      Call dftfun_ecerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,rhoo,
C     >                   Sigma,sigmaco,sigmaoo,
C     >                   Zk1,vrhoc,vrhoo,
C     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)
C
C
      IVer=1
      Call MCORRECTION(ELSM,Occ,TwoNO,URe,
     $                 OrbGrid,WGrid,
     $                 OrbXGrid,OrbYGrid,OrbZGrid,
     $                 NGrid,NBasis,NInte1,NInte2,IVer)
C
      ConCorr(1:200)=Zero
c      A=0.35D0
      A=0.2D0
      B=A-One
      C=2.6D0
c      c=2.0
      G=1.5D0
      D=(C-One)/(One-G)**2
C
      Write(6,'(/,1X,''Values of A, C, G parameters in CASPIDFT'',
     $ 2X,3F10.3)')A,C,G
C
      EDYN=Zero
      ELYP=Zero
      ESR=Zero
      Sum1=Zero
      Sum2=Zero
      Do I=1,NGrid
C
      If(RhoGrid(I).Ne.Zero) Then
C 
      XX=Two*OnTop(I)/RhoGrid(I)**2
C
      If(XX.Le.One) Then      
      PX=A*XX/(One+B*XX)
      Else
      PX=C*XX**0.25-D*(XX-G)**2
      EndIf
C
      EDYN=EDYN+PX*Zk(I)*WGrid(I)
c      EDYN=EDYN+PX*Zk1(I)*WGrid(I)
      ELYP=ELYP+Zk(I)*WGrid(I)
c      ELYP=ELYP+PX*Zk1(I)*WGrid(I)
c      ESR=ESR+Zk1(I)*WGrid(I)

c      if(abs(Zk(I)*WGrid(I)).gt.1.d-4) 
c       write(*,'(i4,x,6e15.4)')i,RhoGrid(I),xx,PX*Zk(I)*WGrid(I),
c     $ Zk(I)*WGrid(I),EDYN,ELYP
C
      EndIf
C
      EndDo
C     TEST TRDMs
      Call TEST_TRDMs(Occ,UNOAO,NInte1,NBasis)
C
c      Sum=Zero
c      Sum11=Zero
c      Sum07=Zero
c      Do I=1,120
c      if(abs(ConCorr(I)).gt.1.d-5) Write(*,'(i4,x,e17.6)'),I,ConCorr(I)
c      Sum=Sum+ConCorr(I)
c      If(I.Gt.100) Sum11=Sum11+ConCorr(I)
c      If(I.Gt.59.And.I.Lt.91) Sum07=Sum07+ConCorr(I)
c      EndDo
c      Write(*,'("Sum:",x,e17.6)')Sum
c      Write(*,'("Contribution from X>1 region: ",x,e17.6)')Sum11
c      Write(*,'("Contribution from 0.90.0>X>0.60 region: ",x,e17.6)')
c     $ Sum07
C
      Write(6,'(/,1X,''LYP Correlation'',7X,F15.8)')ELYP
      Write(6,'(1X,''CASPIDFT Correlation'',2X,F15.8)')EDYN
C      Write(6,'(1X,''SR  Correlation'',7X,F15.8)')ESR
      Write(6,'(1X,''Total CAS+ LYP Energy'',X,F15.8)')ELYP+ETot+ENuc
      Write(6,'(1X,''Total CASPIDFT Energy'',X,F15.8)')EDYN+ETot+ENuc
      Write(6,'(1X,''Total CASPI(M)DFT Energy'',X,F15.8)')
     $ EDYN+ETot+ENuc+ELSM
C
c      Write(6,'(/,1X,'' 1-P(X) Weight'',2X,F15.8)')Sum1
c      Write(6,'(1X,  ''    W=1 Weight'',2X,F15.8)')Sum2
c      Write(6,'(1X,''Total hybrid'',X,F15.8)')
c     $ CMix*ESR+(One-CMix)*Sum1+ETot+ENuc+Sum2
c      Write(6,'(1X,''Total SR-LR'',X,F15.8)')
c     $ ESR+ETot+ENuc+Sum2
C
      Stop 
C
      Return
      End

*Deck GGA_SPIN
      Subroutine GGA_SPIN(Zk,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     RETURNS CORRELATION ENERGY DENSITY USING SPIN FUNCTIONAL
C
C     XCFUN IS USED !!!
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Allocatable :: RDM2Act(:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ WGrid(NGrid),OrbGrid(NGrid,NBasis),OrbXGrid(NGrid,NBasis),
     $ OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C
      Dimension Zk(NGrid),RhoA(NGrid),RhoB(NGrid),
     & SigmaAA(NGrid),SigmaAB(NGrid),SigmaBB(NGrid)
C
      IOnTop=1
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind1(I)=INActive+I
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
C
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
      GoTo 10
   40 Continue
      Close(10)
C
      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid,Occ,URe,OrbGrid,NGrid,NBasis)
C
      If(RhoGrid.Gt.1.D-12) Then
C
      OnTop=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop=OnTop
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
C
      If(IOnTop.Eq.1) Then
C
      Rho=RhoGrid
      R=Two*OnTop/Rho**2
      XFactor=Zero
      If(R.Lt.One) XFactor=SQRT(One-R)
C
      RhoA(I)=Rho/Two*(One+XFactor)
      RhoB(I)=Rho/Two*(One-XFactor)
C
      RhoXa=RhoX/Two*(One+XFactor)
      RhoXb=RhoX/Two*(One-XFactor)
      RhoYa=RhoY/Two*(One+XFactor)
      RhoYb=RhoY/Two*(One-XFactor)
      RhoZa=RhoZ/Two*(One+XFactor)
      RhoZb=RhoZ/Two*(One-XFactor)
C
      SigmaAA(I)=RhoXa*RhoXa+RhoYa*RhoYa+RhoZa*RhoZa
      SigmaAB(I)=RhoXa*RhoXb+RhoYa*RhoYb+RhoZa*RhoZb
      SigmaBB(I)=RhoXb*RhoXb+RhoYb*RhoYb+RhoZb*RhoZb
C     if ontop.eq.1 
      Else
C full spin polarization
      RhoA(I)=RhoGrid
      RhoB(I)=Zero   
      RhoXa=RhoX
      RhoXb=Zero
      RhoYa=RhoY
      RhoYb=Zero
      RhoZa=RhoZ
      RhoZb=Zero

      Fac=1.999
      RhoA(I)=RhoGrid/2.d0*Fac
      RhoB(I)=RhoGrid/2.D0*(Two-Fac)
      RhoXa=RhoX/2.d0*Fac
      RhoXb=RhoX/2.d0*(Two-Fac)
      RhoYa=RhoY/2.d0*Fac
      RhoYb=RhoY/2.d0*(Two-Fac)
      RhoZa=RhoZ/2.d0*Fac
      RhoZb=RhoZ/2.d0*(Two-Fac)

      SigmaAA(I)=RhoXa*RhoXa+RhoYa*RhoYa+RhoZa*RhoZa
      SigmaAB(I)=RhoXa*RhoXb+RhoYa*RhoYb+RhoZa*RhoZb
      SigmaBB(I)=RhoXb*RhoXb+RhoYb*RhoYb+RhoZb*RhoZb
C
      EndIf
C
      Else
C
      RhoA(I)=Zero
      RhoB(I)=Zero
      SigmaAA(I)=Zero
      SigmaAB(I)=Zero
      SigmaBB(I)=Zero
C
      EndIf
C
      EndDo
C
      Call LYP_SPIN(RhoA,RhoB,SigmaAA,SigmaAB,SigmaBB,Zk,NGrid)
C
      Return
      End

*Deck MCORRECTION
      Subroutine MCORRECTION(ELSM,Occ,TwoNO,URe,OrbGrid,WGrid,
     $                       OrbXGrid,OrbYGrid,OrbZGrid,
     $                       NGrid,NBasis,NInte1,NInte2,IVer)
C
C     "MAP" CAS 2-RDM ON GVB 2-RDM BY COUPLING ORBITALS INTO GEMINALS
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension Occ(NBasis),TwoNO(NInte2),URe(NBasis,NBasis)
      Dimension OrbGrid(NGrid,NBasis),WGrid(NGrid),
     $          OrbXGrid(NGrid,NBasis),
     $          OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension C(NBasis),RDM2(NBasis**2*(NBasis**2+1)/2),Ind1(NBasis)
      Double Precision, Allocatable :: OccGH(:),OrbGem(:,:),
     $                                 Zk(:),Sigma(:),
     $                                 RhoGridG(:),RhoGridH(:)
C
      NRDM2 = NBasis**2*(NBasis**2+1)/2
      RDM2(1:NRDM2)=Zero
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct 
C
C     scaling of 2nd variant
      ACoef = 2.5d0
C
      Do I=1,NAct
      Ind1(I)=INActive+I
      EndDo
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
      RDM2(NAddrRDM(J,L,I,K,NBasis))=Half*X
C
      GoTo 10
   40 Continue
      Close(10)
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      Do IR=1,NBasis
      Do IS=1,NBasis
      IAdd=NAddrRDM(IP,IQ,IR,IS,NBasis)
      Hlp=Zero
      If(IP.Eq.IR.And.IQ.Eq.IS.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $  Hlp=Hlp+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $ Hlp=Hlp-Occ(IP)*Occ(IQ)
      If(Hlp.Ne.Zero) RDM2(IAdd)=Hlp
      EndDo
      EndDo
      EndDo
      EndDo
C
      If(NELE-INActive.Ne.NAct-(NELE-INActive)) Then
      Write(6,*) 'Fatal Error: mapping of CAS
     $ on GVB only defined for CAS(m,m)'
c      Stop
      EndIf
C
      NGemSave=NGem
      Do I=1,NBasis
      IGemSave(I)=IGem(I)
      EndDo
      NGem=INActive+NAct/2
      Do I=1,NBasis
      If(Occ(I).Ne.Zero) IGem(I)=0
      If(Occ(I).Eq.Zero) IGem(I)=NGem+1
      EndDo
      Do I=1,INActive
      IGem(I)=I
      EndDo
C
C     COUPLE ORBITALS
C
      Write(6,'(/,X,"Mappinng of CAS(n,n) 2-RDM on GVB-like 2-RDM")')
C
      Allocate(OrbGem(NELE-INActive,2))
C
      II=1
      Do IP=INActive+1,NELE
      IGem(IP)=IP
C
      XDevMax=Zero
      Do IQ=NELE+1,NOccup
C
      XPQPQ=Abs(RDM2(NAddrRDM(IP,IQ,IP,IQ,NBasis))-
     $ 2.0D0*Occ(IP)*Occ(IQ))
      XPQQP=Abs(RDM2(NAddrRDM(IP,IQ,IQ,IP,NBasis))-
     $ (-Occ(IP)*Occ(IQ)))
      XDev=Half*(XPQPQ+XPQQP)/Occ(IQ)
      If(XDev.Gt.XDevMax) Then
      XDevMax=XDev
      IQMax=IQ
      EndIf
C
      CP=SQRT(Occ(IP))
      If(Occ(IP).Lt.Half) CP=-CP
      CQ=SQRT(Occ(IQ))
      If(Occ(IQ).Lt.Half) CQ=-CQ
      Write(*,*)IP,IQ
      Write(*,*)'PPQQ',RDM2(NAddrRDM(IP,IP,IQ,IQ,NBasis)),CP*CQ
      Write(*,*)'PQPQ',RDM2(NAddrRDM(IP,IQ,IP,IQ,NBasis)),
     $ 2.0D0*Occ(IP)*Occ(IQ)
      Write(*,*)'PQQP',RDM2(NAddrRDM(IP,IQ,IQ,IP,NBasis)),
     $ -Occ(IP)*Occ(IQ)
C
      EndDo
C
      If(IGem(IQMax).Eq.0) Then
      IGem(IQMax)=IP
      Write(6,'(X,"**** Orbital",I2," coupled with ",I2 )')IP,IQMax
      OrbGem(II,1)=IP
      OrbGem(II,2)=IQMax
      II = II + 1
      Else
      Write(6,*)
     $ "Warning: more than 2 orbitals assigned to a geminal no",IP
      EndIf
C
      EndDo
C
C
      Do I=1,NGem
C
      Write(6,'(/,X,"Geminal no",I2," includes")')I
      Sum=Zero
C
      II=0
      Do J=1,NBasis
      If(IGem(J).Eq.I) Then
      Sum=Sum+Occ(J)
      Write(6,'(X,"Orbital No: ",I4)')J
      EndIf
      EndDo
      Write(6,'(X,"Norm: ",F12.6)')Sum
C
      EndDo
C
C      If(IVer.Eq.0) Then
C
C     COMPUTE THE M CORRECTION AS GIVEN IN Eq.(26), van Meer et al. JCP 148, 104102 (2018)
C
      ELSM=Zero
      Do I=1,NOccup
      Do J=1,NOccup
      If(IGem(I).Ne.IGem(J)) Then
      ELSM=ELSM+FM(Occ(I),Occ(J))*TwoNO(NAddr3(I,J,I,J))
      EndIf
      EndDo
      EndDo
C
      Write(6,'(/,1X,''M Correlation Correction '',5X,F15.8)')ELSM
C
C      ElseIf(IVer.Eq.1) Then
C
C     MODIFIED M CORRECTION
c
      Allocate(Zk(NGrid),RhoGridH(NGrid),RhoGridG(NGrid))
      Allocate(OccGH(NBasis),Sigma(NGrid))
C
      ELSM=Zero
      Do IG=1,NELE-INActive
      Do IH=1,IG-1
C
      IQ=OrbGem(IG,1)
      IP=OrbGem(IG,2)
      IS=OrbGem(IH,1)
      IR=OrbGem(IH,2)
C 
      XX=max(Occ(IQ),Occ(IP))*(One-max(Occ(IQ),Occ(IP)))
      YY=max(Occ(IS),Occ(IR))*(One-max(Occ(IS),Occ(IR)))
      Coef=SQRT(XX*YY)
C
      OccGH=0
      OccGH(IP)=Occ(IP)
      OccGH(IQ)=Occ(IQ)
      OccGH(IR)=Occ(IR)
      OccGH(IS)=Occ(IS)
C
      Do I=1,NGrid
      Call DenGrid(I,RhoGridH(I),OccGH,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,OccGH,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,OccGH,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,OccGH,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
      EndDo
C
      Call LYP(RhoGridH,Sigma,Zk,NGrid)
C
      VAL=0
      Do I=1,NGrid
      VAL = VAL + Zk(I)*WGrid(I)
      EndDo
C   
      ELSM = ELSM + Coef*VAL
C
      EndDo
      EndDo
C
      ELSM=ACoef*ELSM
C
      Deallocate(OrbGem,OccGH)
      Deallocate(Sigma)
      Deallocate(Zk,RhoGridH,RhoGridG)
C
C      EndIf
C
      Write(6,'(/,1X,''M Correlation Correction '',5X,F15.8)')ELSM
C
      NGem=NGemSave
      Do I=1,NBasis
      IGem(I)=IGemSave(I)
      EndDo
C
      Return
      End

*Deck FM
      Real*8 Function FM(X,Y)
C
C     RETURNS A VALUE OF THE FM FUNCTION DEFINED IN Eqs.(27),(28) van Meer et al. JCP 148, 104102 (2018)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Gamma=1500.D0
C
      XX=X*(One-X)
      YY=Y*(One-Y)
C
      PX=(One+16.D0/Gamma)*Gamma*XX**2/(One+Gamma*XX**2)
      PY=(One+16.D0/Gamma)*Gamma*YY**2/(One+Gamma*YY**2)
C
      FM=-PX*PY*SQRT(XX*YY)
C
      Return
      End

      Subroutine TEST_TRDMS(Occ,UNOAO,NInte1,NBasis) 
C
      use types
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Big=1.5D0)

c
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ UNOAO(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Double precision,dimension(NBasis,NBasis) :: dipx,dipy,dipz
      Double precision :: GammaAB(NInte1),PC(NBasis)
      Double precision :: trdm(NBasis,NBasis),AUXM(NBasis,NBasis),
     $                    WorkV(NBasis)
C
C     Prepare MO->NO
C
      Call read_1rdm_molpro(GammaAB,InSt(1,1),InSt(2,1),
     $ '2RDM',IWarn,NBasis)
C
      Call CpySym(AUXM,GammaAB,NBasis)
      Call Diag8(AUXM,NBasis,NBasis,PC,WorkV)
      Call SortP(PC,AUXM,NBasis)
C
      Sum=Zero
      NAc=0
      Do I=1,NBasis
      Sum=Sum+PC(I)
      If(PC(I).Gt.Zero) NAc=NAc+1
      EndDo
C
      NInAc=XELE-Sum+1.D-2
      INActive=NInAcCAS
C
C     FULL TRANSFORMATION
      URe=0
      Do I=1,NBasis
      IIAct=I-NInAc
      Do J=1,NBasis
      If(I.Eq.J) URe(I,J)=One
      JJAct=J-NInAc
      If(IIAct.Gt.0.And.IIAct.Le.NAc.And.JJAct.Gt.0.And.JJAct.Le.NAc)
     $ URe(I,J)=AUXM(IIAct,JJAct)
      EndDo
      EndDo
C
      Call read_1trdm_molpro(AUXM,InSt(1,1),InTrSt(1,1),
     $ '2RDM',NBasis)
C
      trdm=0
C
      Do J=1,NAc
      Do I=1,NAc
      trdm(INActive+I,INActive+J) = AUXM(I,J)
      Enddo
      Enddo
C
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,
     $           trdm,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           URe,NBasis,0d0,trdm,NBasis)
C
C     READ DIPOLE MOMENT
      Call read_dip_molpro(dipx,dipy,dipz,'DIP',NBasis)
C
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,UNOAO,NBasis,
     $           dipz,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           UNOAO,NBasis,0d0,dipz,NBasis)
C
      TSDipZ=0
      IJ=0
      Do I=1,NBasis
      Do J=1,NBasis
      TSDipZ = TSDipZ + Two*trdm(J,I)*dipz(I,J)
      Enddo
      Enddo
      Write(6,
     $ '(X,"Transition State DMZ: ",F15.8)')TSDipZ
C 
C     REFERENCE STATE
C
      GSDipZ=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,NBasis
      IJ=IJ+1
      If(I.Eq.J) GSDipZ=GSDipZ+Two*Occ(I)*dipz(I,J)
      Enddo
      Enddo
C
      Write(6,
     $ '(X,"Reference  State DMZ: ",F15.8)')GSDipZ 
C
      End

