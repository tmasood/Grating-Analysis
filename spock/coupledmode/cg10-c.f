C   
C
      PROGRAM CG1
      IMPLICIT REAL*8 (A-H,O-Z)
      WRITE(*,11)
   11 FORMAT(/'                                                    '// 
     +        'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'/
     +        '                                                     '/
     +        ' THE MAIN PROGRAM IS ONLY A CALLING PROGRAM.         '/
     +        '                                                     '/
     +        'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
    1 CALL SELECT (N)
      GOTO (10,20,30,40,50,60,70,80,90), N
   10 CALL INDEXA   ! FIND THE REFRACTIVE INDEX OF DIFFERENT MATERIAL
      GOTO 1
   20 CALL INDEX1   ! TMM FIND THE INDEX,NEAR,FAR FIELD AND CONFINEMENT FACTOR
      GOTO 1
   30 CALL COUPLE   ! FIND COUPLING COEFFICIENT
      GOTO 1
   40 CALL DFB      ! FOR PASSIVE DISTRIBUTED FEEDBACK STRUCTURE
      GOTO 1
   50 CALL DFBSL    ! FOR DFB LASER
      GOTO 1           
C   50 CALL DBRSL
C      GOTO 1
   60 CALL SYNCHRO  ! FOR COUPLER (NO GRATING) 
      GOTO 1
   70 CALL NONSYNCH ! FOR COUPLER (GRATING)
      GOTO 1
   80 CALL BRAGG ! FOR BRAGG REFLECTOR
      GOTO 1
   90 PRINT*, ' THIS PROGRAM STOP HERE !'
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SELECT (N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N                                         
      PRINT*, '                                                      '
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      PRINT*, '                                                      '
      PRINT*, ' FROM THE ITEMS DOWN BELOW,INPUT THE SELECTION!       '
      PRINT*, '                                                      '
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      PRINT*, '                                                      '
      WRITE(*,11)
   11 FORMAT(/' ENTER 1 FOR REFRACTIVE INDEX OF DIFFERENT MATERIAL'/
     +        '       2 FOR EFFECTIVE REFRACTIVE INDEX'/
     +        '       3 FOR COUPLEING COEFFICIENT(DFB)'/
     +        '       4 FOR DFB STRUCTURE (PASSIVE) '/
     +        '         ( CALL SELECTION 3 FIRST)'/
     +        '       5 FOR DFB STRUCTURE (LASER) '/
     +        '         ( CALL SELECTION 3 FIRST)'/
C     +        '       5 FOR DBR LASER STRUCTURE       '/
     +        '       6 FOR COUPLER, POWER TRANSFER (W/O GRATING)'/
     +        '       7 FOR COUPLER, POWER TRANSFER (WITH GRATING)  '/
     +        '       8 FOR BRAGG REFLECTOR (STACK AND VCSEL)'/
     +        '       9 FOR EXIT                      '/
     +        '                                                      '/
     +        '  PLEASE INPUT THE SELECTION !                       ')
      READ(*,*) N
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE INDEX1
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LAMBDA,KO,N(500),W(500),ALPHA(500),K(500),D(500),DX(500),BR,BI
      REAL*8 FF,FX(181),FY(181),T(181),CONF(500)
      COMPLEX*16 ZG,E(500),ZN(500),CI,A(500),B(500),BX,Z1,Z2,ZGRX
      COMPLEX*16 ZGA,ZGB,XDG(500),XBX(500),HH,AA,BB,BXA,H(500)
      INTEGER NN,IC,XI(500),NL 
      CHARACTER IFILE*40
      COMMON DX,E,NL,IC,CI
      COMMON /STATE/ BX,ko
      COMMON /PSIF/HH,AA,BB,CONF 
      COMMON /FAR1/ H,A,B,XNORN
      COMMON /XMAX/ FX
      COMMON /A1/ N
      COMMON /A10/ LAMBDA,BR,BI
	COMMON /A11/ AIN
      EXTERNAL TE,TM,GSP                  
      MIN=5
      PI=3.1415926535897932D0
      PRINT*,'*******************************************************'
      PRINT*,'                                                      '          
      PRINT*,' THIS USE TRANSFER MATRIX METHOD TO FIND THE EFFECTIVE'          
      PRINT*,' REFRACTIVE INDEX. AND COMPLEX ANSWERS ARE CALCULATED.'          
      PRINT*,'                                                      ' 
      PRINT*,' ITTER=2 -- CONVERGENCE CONDITION REACHED FOR F(X)    '
      PRINT*,' ITTER=1 -- CONVERGENCE CONDITION REACHED FOR X       '
      PRINT*,' ITTER=0 -- WRONG ANSWER.                             '
      PRINT*,'                                                      ' 
      PRINT*,'*******************************************************'
      REWIND (1)
      REWIND (3)
      REWIND (4)
 9997 PRINT*,' THE INPUT FILE NAME='
      READ(MIN,'(A40)') IFILE
      OPEN(1,FILE=IFILE,status='unknown')
      OPEN(3,FILE='field.dat',status='unknown')
      OPEN(4,FILE='far.dat',status='unknown')
      OPEN(7,FILE='con-f.dat',status='unknown')
      READ(1,*) NL
      READ(1,*) LAMBDA
      do i=1,nl
      READ(1,*) N(I),ALPHA(I),W(I) 
      WRITE(*,7776) N(I),ALPHA(I),W(I),I
 7776 FORMAT(' N=',D12.6,' ALPHA=',D12.6,' W=',D14.7,' I=',I3)     
      ENDDO
      PRINT*, ' INPUT THE CENTER LAYER IC = '
      READ(*,*) IC 
      TPI=2.D0*PI
      KO=TPI/(LAMBDA) 
      EPO=8.85418782D-12
      UO=4.*PI*1.D-7
      C=3.D8                                     
      CI=DCMPLX(0.D0,1.D0)
      DO I=1,NL
c        K(I)=ALPHA(I)/(2.D0*KO)
	 K(I)=ALPHA(I)/KO
         ZN(I)=DCMPLX(N(I)*KO,K(I)*KO)
         E(I)=ZN(I)*ZN(I) 
         D(I)=W(I)                         ! INPUT THE DISTANCE
      ENDDO
      DX(IC-1)=D(IC)/2.D0                  ! FIND THE CENTER LAYER
      DX(IC)=(-D(IC))/2.D0
      DO I=IC-2,1,-1                       ! FROM LEFT TO CENTER
         DX(I)=DX(I+1)+D(I+1)
      ENDDO
      DO I=IC+1,NL-1,1                     ! FROM RIGHT TO CENTER
         DX(I)=DX(I-1)-D(I)      
      ENDDO 
  109 PRINT*, ' INPUT NN=1 FOR TE MODE, NN=2 FOR TM MODE, NN=3 STOP'
      PRINT*,' NN='       
      READ(*,*) NN
      IF (NN.EQ.1) THEN
         GOTO 309
         ELSEIF (NN.EQ.2) THEN
                GOTO 409
         ELSEIF (NN.EQ.3) THEN
                GOTO 209
         ELSE
      ENDIF
  309 PRINT*, ' INPUT INITIAL GUESS'
      PRINT*, ' Z1R=',' Z1I=',' Z2R=',' Z2I='
      READ(*,*) Z1R,Z1I,Z2R,Z2I
      Z1=DCMPLX(Z1R*KO,Z1I*KO)
      Z2=DCMPLX(Z2R*KO,Z2I*KO)
      AIN=1.D-4
      EPS=1.D-6
      EPSS=1.D-8
      DO I=1,20
         ZGRX=(Z2-Z1)/21
         ZG=Z1+ZGRX*I          
         CALL CROOT (TE,ZG,ITTER,J,EPS,EPSS)
         XDG(I)=ZG                        
         XBX(I)=BX
         XI(I)=ITTER                
      ENDDO                  
      DO I=1,20 
         CALL SEARCH (XDG,XBX,XI)    
      ENDDO
      DO I=20,1,-1        
      WRITE(*,502) XDG(I)/KO,XBX(I),XI(I)
      ENDDO
      GOTO 509 
  409 PRINT*, ' INPUT INITIAL GUESS'
      PRINT*, ' Z1R=',' Z1I=',' Z2R=',' Z2I='
      READ(*,*) Z1R,Z1I,Z2R,Z2I
      Z1=DCMPLX(Z1R*KO,Z1I*KO)
      Z2=DCMPLX(Z2R*KO,Z2I*KO)  
      AIN=1.E-4
      EPS=1.D-6
      EPSS=1.D-8
      DO I=1,20
         ZGRX=(Z2-Z1)/21
         ZG=Z1+ZGRX*I          
         CALL CROOT (TM,ZG,ITTER,J,EPS,EPS)
         XDG(I)=ZG                        
         XBX(I)=BX
         XI(I)=ITTER
      ENDDO
      DO I=1,20
         CALL SEARCH (XDG,XBX,XI)
      ENDDO
      DO I=20,1,-1 
      WRITE(*,502) XDG(I)/KO,XBX(I),XI(I)
      ENDDO
  502 FORMAT(' MODE=',2(1X,(D13.7)),' BX=',2(1X,(D13.7)),' ITTER=',I4)
  509 PRINT*, ' ACCEPT THE ANSEWR INPUT NN=1 OTHER MODE NN=2'
      PRINT*, 'N='
      READ(*,*) NX
      IF (NX.EQ.1) THEN
      GOTO 999
      ELSEIF (NX.EQ.2) THEN
      GOTO 109
      ELSE
      ENDIF
  999 IF (NN.EQ.1) THEN
      PRINT*, ' INPUT THE BETA VALUE, AND BX VALUE'
      PRINT*, ' BR=',' BI=',' BXR=',' BXI='
      READ(*,*) BR,BI,BXR,BXI
      ZGA=DCMPLX(BR*KO,BI*KO)      
      BXA=DCMPLX(BXR,BXI)
      CALL TEFUN (ZGA,BXA)  
      PRINT*, '                                                  '
      DO I=1,NL                                                  
      WRITE(7,503) I,CONF(I)
      WRITE(*,503) I,CONF(I)
      ENDDO
      ELSEIF (NN.EQ.2) THEN
      PRINT*, ' INPUT THE BETA VALUE, AND BX VALUE'
      PRINT*, ' BR=',' BI=',' BXR=',' BXI='
      READ(*,*) BR,BI,BXR,BXI
      ZGB=DCMPLX(BR*KO,BI*KO)        
      BXA=DCMPLX(BXR,BXI)
      CALL TMFUN (ZGB,BXA)   
      PRINT*, '                                                 '
      DO I=1,NL                                                
      WRITE(*,503) I,CONF(I)
      ENDDO                                            
      ELSE
      ENDIF
      PRINT*, ' FOR FAR FIELD INPUT 1, FOR OTHER INPUT 2'
      READ(*,*) IJ
      IF (IJ.EQ.1) THEN
      GOTO 1111
      ELSEIF (IJ.EQ.2) THEN
      GOTO 1112
      ELSE
      ENDIF
 1111 DO I=1,179,1 
         T(I)=(I-90)*1.D0
         THETA=(PI/180)*T(I)
         CALL FFIELD (THETA,BXA,FF)
         FX(I)=FF
         FY(I)=FF  
      ENDDO
         CALL XSEARCH(XFF)
      DO I=1,179,1      
         T(I)=(I-90)*1.D0
         FY(I)=FY(I)/XFF 
         WRITE(4,504) T(I),FY(I)
      ENDDO
      DO I=1,89            
         T(I)=(I-90)*1.d0
         IF (FY(I).GE.0.49D0) THEN
         II=I    
          X2=FY(II)
          X3=FY(II+1)
          X1=FY(II-1)
          Y2=DABS(T(II))
          Y3=DABS(T(II)+1)
          Y1=DABS(T(II)-1)
         GOTO 510
         ELSE
         ENDIF
      ENDDO
  510 CONTINUE 
      CALL CSPLINE (X1,X2,X3,Y1,Y2,Y3,SOL1) 
      DO I=90,1,-1 
         JJ=I+89
         T(JJ)=(JJ)*1.D0
         IF (FY(JJ).GE.0.49D0) THEN 
         KK=JJ
         X2=FY(KK)
         X3=FY(KK+1)
         X1=FY(KK-1)
         Y2=T(KK)-90.D0
         Y3=(T(KK)+1)-90.D0
         Y1=(T(KK)-1)-90.D0
         GOTO 511
         ELSE
         ENDIF                       
      ENDDO
  511 CONTINUE                            
      CALL CSPLINE (X1,X2,X3,Y1,Y2,Y3,SOL2)
      FWHPF=(SOL2+SOL1)                
      WRITE(*,505) FWHPF
  503 FORMAT(2X,'CONFINEMENT FACTOR OF ',I3,' th LAYER=',E14.8)
  504 FORMAT(1X,F8.3,1X,E14.7)                   
  505 FORMAT(2X,'FWHPF=',E12.6/)
  209 CLOSE (1)
      CLOSE (3)
      CLOSE (4)
 1112 PRINT*, ' FOR NEW INPUT DATAS INPUT I=1 FOR OTHER FUNCTIONS I=2'
      PRINT*, ' INPUT = '
      READ(*,*) I
      IF (I.EQ.1) THEN  
         GOTO 9997
      ELSE       
         GOTO 9999
      ENDIF
 9999 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c      SUBROUTINE INDEXA
C
c      IMPLICIT REAL*8 (A-H,O-Z)
c    1 CALL SELECTA(N)
c      GOTO (10,30) N
c   10 CALL SINGLE
c      GOTO 1
c   30 PRINT*, ' THIS PROGRAM STOP HERE!'
c      RETURN
c      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
c      SUBROUTINE SELECTA (N)
c      IMPLICIT REAL*8 (A-H,O-Z)
c      INTEGER N
c      WRITE(*,11)
c   11 FORMAT(/' CALCULATE THE WAVE LENGTH DEPENDENT REFRACTIVE INDEX'/
c     +        ' ENTER 1 FOR SINGLE LAYER REFRACTIVE INDEX'/  
c     +        '       2 FOR STOP')
c      READ(*,*) N
c      RETURN
c      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c      SUBROUTINE SINGLE
      SUBROUTINE INDEXA
      IMPLICIT REAL*8 (A-H,O-Z)
    1 CALL SELECT2 (N)
      GOTO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160) N
   10 CALL sInGaAsP   ! 1.35-0.72*Y+0.12*Y**2
      GOTO 1
   20 CALL sInGaAsP2
      GOTO 1
   30 CALL sAlGaInAs  ! 0.75+1.548*X
      GOTO 1
   40 CALL sAlGaInAs2 ! Adachi's
      GOTO 1
   50 CALL sAlGaPSb   ! 0.885+2.33*X-1.231*X**2
      GOTO 1
   60 CALL sAlGaAsSb  ! 0.829+1.822*X-0.22*X**2
      GOTO 1
   70 CALL sAlInAsSb  ! -0.173+4.014*X-1.418*X**2
      GOTO 1
   80 CALL sGaInPSb   ! 1.326-0.018*X-0.34*X**2  
      GOTO 1
   90 CALL sAlAsSb    ! 1.659+0.523*X-0.002*X**2(1.7+0.53*X)
      GOTO 1              
  100 CALL sAlGaAs    ! 1.424+1.247*X (X < 0.45)
C                     ! 1.424+1.247*X+1.147*(X-0.45)**2
      GOTO 1
  110 CALL sAlGaAs2
      GOTO 1
  120 CALL sAlGaInP   ! 1.9+0.6*X
      GOTO 1
  130 CALL sInGaAs1   !In(1-x)Ga(x)As/InP
      GOTO 1
  140 CALL sInGaAs2   !In(1-x)Ga(x)As/GaAs
      GOTO 1
C 110 CALL sGaInP
C     GOTO 1
  150 CALL sInGaAsN   !In(y)Ga(1-y)As(x)N(1-x)
      GOTO 1
  160 PRINT*, ' THIS PROGRAM STOP HERE!'
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sInGaAsP
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' INPUT X,Y AND PHOTON WAVELENGTH FOR InGaAsP'
      PRINT*, ' X=, Y=, LAMBDA='
      READ(*,*) X,Y,XLAMBDA
      HV=1.24/XLAMBDA
      PRINT*, ' FOR LATTICE MATCHED, INPUT 1,--STRAIN INPUT 2' 
      PRINT*, ' INPUT = ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      A=8.4-3.4*Y
      B=6.6+3.4*Y
      EO=1.35-0.72*Y+0.12*(Y**2)
      EDD=1.466-0.557*Y+0.129*Y**2
      ELSEIF (I.EQ.2) THEN
      A=(1-X)*Y*5.14+(1-X)*(1-Y)*8.4+X*Y*6.3+X*(1-Y)*22.25
      B=(1-X)*Y*10.15+(1-X)*(1-Y)*6.6+X*Y*9.4+X*(1-Y)*0.9
c      EO=(1-X)*Y*0.36+(1-X)*(1-Y)*1.35+X*Y*1.424+X*(1-Y)*2.74
      EO=1.35-0.72*Y+0.12*(Y**2)
      EDO=(1-X)*Y*0.38+(1-X)*(1-Y)*0.11+X*Y*0.34+X*(1-Y)*0.08
c      EDD=1.466-0.557*Y+0.129*Y**2
      EDD=EDO+EO
      ELSE
      ENDIF
      XO=HV/EO
      XSO=HV/EDD
	Q1=MAX(1-XO,0.D0)
	Q2=MAX(1-XSO,0.D0)
	F1=(XO**(-2))*(2-((1+XO)**0.5)-(Q1**0.5))
	F2=(XSO**(-2))*(2-((1+XSO)**0.5)-(Q2**0.5))
C      F1=(XO**(-2))*(2-((1+XO)**0.5)-((1-XO)**0.5))
C      F2=(XSO**(-2))*(2-((1+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/(EDD))**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sInGaAsP2
      IMPLICIT REAL*8 (A-H,O-Z)
      PI=3.14159265D0
    1 PRINT*, ' FOR NON-ADACHI FORMULA.'
      PRINT*, ' 1 FOR BOTH WAVELENGTH AND COMPOSITION REQUIRE'
      PRINT*, ' 2 FOR WAVELENGTH ONLY'
      PRINT*, ' INPUT J FOR SELECTION NUMBER'
      READ(*,*) J
      PRINT *, ' INPUT LAYER PL WAVELENGTH XPL=?'
      READ(*,*) XPL
      EP=1.24/XPL
      PRINT *, ' INPUT TARGET WAVELENGTH, XLAMBDA=?'
      READ(*,*) XLAMBDA
      E=1.24/XLAMBDA
      IF (J.EQ.1) THEN
      PRINT*, ' INPUT X AND Y '
      READ(*,*) X,Y
      EO=0.595*X*X*(1-Y)+1.626*X*Y-1.891*Y+0.524*X+3.391
      ED=(12.36*X-12.71)*Y+7.54*X+28.91
      ETHA=PI*ED/(2*EO**3*(EO**2-EP**2))
      XN=(1+(ED/EO)+(ED*E**2/(EO**3))+(ETHA*E**4/PI)
     +  *DLOG((2*EO**2-EP**2-E**2)/(EP**2-E**2)))
      XN=DSQRT(XN)
      ELSEIF (J.EQ.2) THEN
      E1=2.5048
      E2=0.1638
      A1=13.3510-5.4554*EP+1.2332*EP*EP
      A2=0.7140-0.3606*EP
      ESP=1.D0+(A1/(1-((E/(EP+E1))**2)))+(A2/(1-((E/(EP+E2))**2)))
      XN=DSQRT(ESP)
      ELSE
      ENDIF
      PRINT*, ' INDEX N= ', XN
      PRINT*, ' IF NEW INPUT, I=1, RETURN TO MAIN I=2'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      GOTO 2
      ENDIF
    2 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlGaInAs
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LAMBDA
   10 PRINT*, ' INPUT XX, QY AND PHOTON WAVELENGTH FOR AlGaInAs'
      PRINT*, ' XX= , QY=, LAMBDA='
      READ(*,*) XX,QY,LAMBDA
       XLAMBDA=LAMBDA*1000.D0
       X=XX/0.48
       IF ((X.LT.0.3D0).AND.(LAMBDA.EQ.1.3D0)) THEN
       XN=DSQRT((1.-XX-QY)*14.6+QY*13.2+XX*10.06)-0.17d0
       ELSEIF ((X.LT.0.3D0).AND.(LAMBDA.EQ.1.55D0)) THEN
      XN=DSQRT((1.-XX-QY)*14.6+QY*13.2+XX*10.06)-0.31D0
      GOTO 101
       ELSEIF ((X.LT.0.3D0).AND.(LAMBDA.EQ.1.4D0)) THEN
       XN=DSQRT((1.-XX-QY)*14.6+QY*13.2+XX*10.06)-0.226d0
c      PRINT*, ' THE REFRACTIVE INDEX IS=', XN
       ELSEIF (X.GE.0.3) THEN
       A=9.689-1.012*X
       B=1.590-0.376*X
       C=1102.4-702.0*X+330.4*(X*X)                  
       XN=DSQRT(A+((B*XLAMBDA**2)/(XLAMBDA**2-C**2)))      
C      PRINT*, ' THE REFRACTIVE INDEX IS=', XN
       ELSE
       ENDIF 
cccccccccccccccccccccccccccccccccccccccccccccc
c     HV=1.24D0/LAMBDA
c     EO=0.75D0+1.548D0*XX
c     EDO=XX*0.28+QY*0.34+(1-XX-QY)*0.38
c     XO=HV/EO  
c     EDD=EDO+EO    
c     XSO=HV/EDD                  
c     A=XX*25.30+QY*6.30+(1-XX-QY)*5.14
c     B=XX*(-0.80)+QY*9.40+(1-XX-QY)*10.15
c     F1=(XO**(-2.D0))*(2.D0-((1.D0+XO)**0.5)-((1-XO)**0.5))
c     F2=(XSO**(-2.D0))*(2.D0-((1.D0+XSO)**0.5)-((1-XSO)**0.5))
c     X1=((EO/EDD)**(1.5))/2
c     X2=F1+X1*F2
c     X3=A*X2
c     X4=X3+B
c     XN=DSQRT(X4)
  101 PRINT*, ' THE REFRACTIVE INDEX IS=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I='
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END      
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlGaInAs2
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LAMBDA
   10 PRINT*, ' INPUT XX, QY AND PHOTON WAVELENGTH FOR AlGaInAs'
      PRINT*, ' XX= , QY=, LAMBDA='
      READ(*,*) XX,QY,LAMBDA
      HV=1.24D0/LAMBDA
      EO=0.75D0+1.548D0*XX
      EDO=XX*0.28+QY*0.34+(1-XX-QY)*0.38
      XO=HV/EO
      EDD=EDO+EO
      XSO=HV/EDD
      A=XX*25.30+QY*6.30+(1-XX-QY)*5.14
      B=XX*(-0.80)+QY*9.40+(1-XX-QY)*10.15
	Q1=MAX(1-XO,0.D0)
	Q2=MAX(1-XSO,0.D0)
	F1=(XO**(-2))*(2-((1+XO)**0.5)-(Q1**0.5))
	F2=(XSO**(-2))*(2-((1+XSO)**0.5)-(Q2**0.5))
C      F1=(XO**(-2.D0))*(2.D0-((1.D0+XO)**0.5)-((1-XO)**0.5))
C      F2=(XSO**(-2.D0))*(2.D0-((1.D0+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/EDD)**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' THE REFRACTIVE INDEX IS=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I='
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlGaPSb
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' INPUT X,Y AND PHOTON WAVELENGTH FOR AlGaPSb'
      PRINT*, ' X=, Y=, LAMBDA='
      READ(*,*) X,Y,XLAMBDA
      HV=1.24/XLAMBDA
      A=(1-X)*Y*22.25+(1-X)*(1-Y)*4.05+X*Y*24.10+X*(1-Y)*59.68
      B=(1-X)*Y*0.9+(1-X)*(1-Y)*12.66+X*Y*(-2.0)+X*(1-Y)*(-9.53)
c      EO=(1-X)*Y*2.74+(1-X)*(1-Y)*0.72+X*Y*3.58+X*(1-Y)*2.22
      EO=0.885+2.33*X-1.231*(X**2)
      EDO=(1-X)*Y*0.08+(1-X)*(1-Y)*0.82+X*Y*0.07+X*(1-Y)*0.65
      EDD=EDO+EO
      XO=HV/EO
      XSO=HV/EDD
      F1=(XO**(-2))*(2-((1+XO)**0.5)-((1-XO)**0.5))
      F2=(XSO**(-2))*(2-((1+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/EDD)**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlGaAsSb
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' INPUT X,Y AND PHOTON WAVELENGTH FOR AlGaAsSb'
      PRINT*, ' X=, Y=, LAMBDA='
      READ(*,*) X,Y,XLAMBDA
      HV=1.24/XLAMBDA
      A=(1-X)*Y*6.30+(1-X)*(1-Y)*4.05+X*Y*25.30+X*(1-Y)*59.68
      B=(1-X)*Y*9.4+(1-X)*(1-Y)*12.66+X*Y*(-0.8)+X*(1-Y)*(-9.53)
C      EO=(1-X)*Y*1.42+(1-X)*(1-Y)*0.72+X*Y*2.95+X*(1-Y)*2.22
      EO=0.829+1.822*X-0.22*(X**2)
      EDO=(1-X)*Y*0.34+(1-X)*(1-Y)*0.82+X*Y*0.28+X*(1-Y)*0.65
      EDD=EDO+EO
      XO=HV/EO
      XSO=HV/EDD
	Q1=MAX(1-XO,0.D0)
	Q2=MAX(1-XSO,0.D0)
	F1=(XO**(-2))*(2-((1+XO)**0.5)-(Q1**0.5))
	F2=(XSO**(-2))*(2-((1+XSO)**0.5)-(Q2**0.5))
C      F1=(XO**(-2))*(2-((1+XO)**0.5)-((1-XO)**0.5))
C      F2=(XSO**(-2))*(2-((1+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/EDD)**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlInAsSb
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' INPUT X,Y AND PHOTON WAVELENGTH FOR AlInAsSb'
      PRINT*, ' X=, Y=, LAMBDA='
      READ(*,*) X,Y,XLAMBDA
      HV=1.24/XLAMBDA
      A=(1-X)*Y*5.14+(1-X)*(1-Y)*7.91+X*Y*25.30+X*(1-Y)*59.68
      B=(1-X)*Y*10.15+(1-X)*(1-Y)*13.07+X*Y*(-0.8)+X*(1-Y)*(-9.53)
C      EO=(1-X)*Y*0.36+(1-X)*(1-Y)*0.17+X*Y*2.95+X*(1-Y)*2.22
      EO=-0.173+4.014*X-1.418*(X**2)
      EDO=(1-X)*Y*0.38+(1-X)*(1-Y)*0.81+X*Y*0.28+X*(1-Y)*0.65
      EDD=EO+EDO
      XO=HV/EO
      XSO=HV/EDD
	Q1=MAX(1-XO,0.D0)
	Q2=MAX(1-XSO,0.D0)
	F1=(XO**(-2))*(2-((1+XO)**0.5)-(Q1**0.5))
	F2=(XSO**(-2))*(2-((1+XSO)**0.5)-(Q2**0.5))
C      F1=(XO**(-2))*(2-((1+XO)**0.5)-((1-XO)**0.5))
C      F2=(XSO**(-2))*(2-((1+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/EDD)**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sGaInPSb
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' INPUT X,Y AND PHOTON WAVELENGTH FOR GaInPSb'
      PRINT*, ' X=, Y=, LAMBDA='
      READ(*,*) X,Y,XLAMBDA
      HV=1.24/XLAMBDA
      A=(1-X)*Y*8.4+(1-X)*(1-Y)*7.91+X*Y*22.25+X*(1-Y)*4.05
      B=(1-X)*Y*6.6+(1-X)*(1-Y)*13.07+X*Y*0.9+X*(1-Y)*12.66
C      EO=(1-X)*Y*1.35+(1-X)*(1-Y)*0.17+X*Y*2.74+X*(1-Y)*0.72
      EO=1.326-0.018*X-0.34*(X**2)
      EDO=(1-X)*Y*0.11+(1-X)*(1-Y)*0.81+X*Y*0.08+X*(1-Y)*0.82
      EDD=EDO+EO
      XO=HV/EO
      XSO=HV/EDD
	Q1=MAX(1-XO,0.D0)
	Q2=MAX(1-XSO,0.D0)
	F1=(XO**(-2))*(2-((1+XO)**0.5)-(Q1**0.5))
	F2=(XSO**(-2))*(2-((1+XSO)**0.5)-(Q2**0.5))
C      F1=(XO**(-2))*(2-((1+XO)**0.5)-((1-XO)**0.5))
C      F2=(XSO**(-2))*(2-((1+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/EDD)**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END             
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlAsSb 
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' INPUT X,Y AND PHOTON WAVELENGTH FOR AlAsSb'
      PRINT*, ' X= ?,  LAMBDA='
      READ(*,*) X,XLAMBDA
      HV=1.24/XLAMBDA                
      A=(1-X)*59.68+X*25.30
      B=(1-X)*(-9.53)+X*(-0.80)
      EO=1.7+0.53*X
      EDO=(1-X)*0.65+X*0.28
      EDD=EDO+EO
      XO=HV/EO
      XSO=HV/EDD
	Q1=MAX(1-XO,0.D0)
	Q2=MAX(1-XSO,0.D0)
	F1=(XO**(-2))*(2-((1+XO)**0.5)-(Q1**0.5))
	F2=(XSO**(-2))*(2-((1+XSO)**0.5)-(Q2**0.5))
C      F1=(XO**(-2))*(2-((1+XO)**0.5)-((1-XO)**0.5))
C      F2=(XSO**(-2))*(2-((1+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/EDD)**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlGaAs
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 F1,F2,X2,X3,X4,XN,XO,XSO,CHI,CHISO
   10 PRINT*, ' INPUT X, AND PHOTON WAVELENGTH FOR AlGaAs'
      PRINT*, ' X=, LAMBDA='
      READ(*,*) X,XLAMBDA
      HV=1.24/XLAMBDA
      A=6.3+19.0*X
      B=9.4-10.2*X
c     IF (X.LT.0.45) THEN
c     EO=1.424+1.247*X
c     ELSEIF (X.GT.0.45) THEN
c     EO=1.424+1.247*X+1.147*((X-0.45)**2)
c     EDO=0.34-0.5*X
c     ELSE
c     ENDIF
c     EDD=EDO+EO
      EO=1.425+1.155*X+0.37*X*X
      EDD=1.765+1.115*X+0.37*X*X
      XO=HV/EO
      XSO=HV/EDD
      CHI=DCMPLX(1.D0,0.D0)-XO
      CHISO=DCMPLX(1.D0,0.D0)-XSO
      IF ((DREAL(CHI).LT.0.D0).OR.(DREAL(CHISO).LT.0.D0)) THEN
      PRINT*, ' RANGE ERROR ! COMPLEX SQUARE ROOT!'
      ELSE
      ENDIF
      F1=(XO**(-2))*(2-(CDSQRT(1.d0+XO))-(CDSQRT(CHI)))
      F2=(XSO**(-2))*(2-(CDSQRT(1.d0+XSO))-(CDSQRT(CHISO)))
      X1=((EO/(EDD))**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=CDSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlGaAs2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 ARG1,ARG2,CEOO,CEC,CE,CEO,CEXPO,CEXPOO
      COMPLEX*16 ARG3,ARG4,ARG5,ARG6,ARG7,ARG8
      REAL*8 NEFF,LAMBDA,LAMBDA1
 1000 PRINT*, ' INPUT THE Al FRACTION X, AND WAVELENGTH WVL'
      PRINT*, ' X=?,  WVL=?'
      READ(*,*) X,LAMBDA
      PI=3.141592653589D0
      AP=(X-0.5D0)
      APP=(X-0.5D0)*(X-0.5D0)
      A=4.264374+(4.402701*AP)+(9.160232*APP)
      B1=8.256268+(-0.585350*AP)+(-1.901850*APP)
      B2=0.462866+(0.012206*AP)+(1.047697*APP)
      B11=19.788257+(-3.706511*AP)+(2.990894*APP)
      B22=0.539077+(-0.172307*AP)+(1.031416*APP)
      C=3.078636+(-1.598544*AP)+(-0.742071*APP)
      D=29.803617+(-22.124036*AP)+(-57.863105*APP)
      E1=3.212414+(0.804397*AP)+(0.228317*APP)
      EC=4.724383+(0.024499*AP)+(0.030653*APP)
      EIC=3.160138+(0.138046*AP)+(1.066214*APP)
      EINF=-0.494941+(-0.030802*AP)+(0.486201*APP)
      E1C=6.413109+(0.571119*AP)+(-0.735610*APP)
      T=0.263476+(0.090532*AP)+(-0.254099*APP)
      G=0.147517+(-0.068764*AP)+(0.047345*APP)
      F=1.628226+(-1.682422*AP)+(-2.081273*APP)
      FI=0.507707+(-0.070165*AP)+(-0.122169*APP)
C
      EO=1.425000+(1.155000*X)+(0.370000*(X*X))
      EI=1.734000+(0.574000*X)+(0.055000*(X*X))
      DO=0.340000+(0.0*X)+(0.0*(X*X))
      D1=0.230000+(-0.030000*X)+(0.0*(X*X))
C
      E=(4.135701327D-15*2.99792458D14)/LAMBDA
      EE=E*E
      EEI=1.D0/EE
C
CCCCCC REAL PART OF DIELECTRIC FUNCTION : DIRECT EDGE (EO)
C
      ARG1=DCMPLX(EO+E,0.D0)
      ARG2=DCMPLX(EO-E,0.D0)
      IF (EO.GE.E) THEN
      ARG2=ARG2
      ELSE
      ARG2=DCMPLX(0.D0,0.D0)
      ENDIF
      EP1EO1=A*EEI*DREAL(2.D0*DSQRT(EO)-CDSQRT(ARG1)-CDSQRT(ARG2))
      ARG3=DCMPLX(EO+DO+E,0.D0)
      ARG4=DCMPLX(EO+DO-E,0.D0)
      IF (EO+DO.GE.E) THEN
      ARG4=ARG4
      ELSE
      ARG4=DCMPLX(0.D0,0.D0)
      ENDIF
      EP1EO2=(A*EEI/2.D0)*DREAL(2.D0*DSQRT(EO+DO)-CDSQRT(ARG3)
     +      -CDSQRT(ARG4))    
      EP1EO=EP1EO1+EP1EO2
C
CCCCCC IMAGINARY PART OF DIELECTRIC FUNCTION
C
      ARG5=DCMPLX(E-EO,0.D0)
      ARG6=DCMPLX(E-EO-DO,0.D0)
      IF (E.GE.EO) THEN
      ARG5=ARG5
      ELSE
      ARG5=DCMPLX(0.D0,0.D0)
      ENDIF
      IF (E.GE.(EO+DO)) THEN
      ARG6=ARG6
      ELSE
      ARG6=DCMPLX(0.D0,0.D0)
      ENDIF
      EP2EO=DREAL((2.D0*CDSQRT(ARG5)+CDSQRT(ARG6))*EEI)*A/3.D0  
C
CCCCCC REAL PART OF DIELECTRIC FUNCTION : L=[K=PI/A]
C
      CEOO=E*DCMPLX(1.D0,0.D0)+T*DCMPLX(0.D0,1.D0) !E=hw+iGamma     
      CEC=1.D0/(DCMPLX(1.D0,0.D0)-((CEOO/E1C)*(CEOO/E1C)))
      CE=CEOO*CEOO
      CEO=CE/((E1+D1)*(E1+D1))
      CE=CE/(E1*E1)
      IF (E.GE.E1) THEN
      CEXPO=DEXP(-F*(E-E1))/CE
      ELSE
      CEXPO=1.D0/CE
      ENDIF
      IF (E.GE.(E1+D1)) THEN
      CEXPOO=DEXP(-F*(E-E1-D1))/CEO
      ELSE
      CEXPOO=1.D0/CEO
      ENDIF
      EP1E11=-B1*DREAL(LOG(((1.D0,0.D0)-CE)*CEC)*CEXPO)
      EP1E12=-B2*DREAL(LOG(((1.D0,0.D0)-CEO)*CEC)*CEXPOO)
      EP1E1=EP1E11+EP1E12
C
CCCCCC IMAGINARY PART
C
      ARG7=DCMPLX(E1,0.D0)-CEOO
      ARG8=DCMPLX(E1+D1,0.d0)-CEOO
      EP2E11=PI*MAX(DREAL((B1-B11*CDSQRT(ARG7))*CEXPO),0.D0)
      EP2E12=PI*MAX(DREAL((B2-B22*CDSQRT(ARG8))*CEXPOO),0.D0)
c     EP2E1=MAX(EP2E11+EP2E12,0.D0)
      EP2E1=EP2E11+EP2E12
C
CCCCCC REAL PART
C
      CE=DCMPLX(1.D0,0.D0)-EE/(EC*EC)
      EP2EC=DREAL(C/(CE*CE+EE*((G/EC)*(G/EC))))
      EP1EC=DREAL(CE*EP2EC)
C
CCCCCC IMAGINARY PART
C
      EP2EC=(E*G/EC)*EP2EC
C
CCCCCC IMAGINARY PART
C
      EQ1=E-EI
      IF (E.GE.EIC) THEN
      CEXPOOO=DEXP(-FI*(E-EIC))
      ELSE
      CEXPOOO=1.D0
      ENDIF
      EP2EI=(D*EEI*(EQ1*EQ1))*CEXPOOO
C
CCCCCC HEAVYSIDE FUNCTION
C
      IF (E.LT.EO) THEN
         EP2E1=0.D0
         EP2EC=0.D0
         EP2EI=0.D0
      ENDIF
      IF (E.LT.E1) THEN
         EP2E1=0.D0
      ENDIF
C
      EPSILON1=EINF+EP1EO+EP1E1+EP1EC
      EPSILON2=EP2EO+EP2E1+EP2EC+EP2EI
      ARG1=DCMPLX(EPSILON1,EPSILON2)
      NEFF=DREAL(CDSQRT(ARG1))
      PRINT*, ' N= ', NEFF
      PRINT*, ' FOR NEW CALCULATION, INPUT 1, FOR STOP INPUT 2 '
      PRINT*, ' INPUT =?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1000
      ELSE
      GOTO 2000
      ENDIF
 2000 PRINT*, ' PROGRAM STOP HERE!'
      RETURN
      END
c      SUBROUTINE sAlGaAs2
c      IMPLICIT REAL*8 (A-H,O-Z)
c      COMPLEX*16 ARG1,ARG2,CEOO,CEC,CE,CEO,CEXPO,CEXPOO
c      REAL*8 NEFF
c 1000 PRINT*, ' INPUT THE Al FRACTION X, AND WAVELENGTH WVL'
c      PRINT*, ' X=?,  WVL=?'
c      READ(*,*) P,WVL 
c      PI=3.141592653589D0
c      A=4.264374+(4.402701*(P-0.5))+(9.160232*((P-0.5)*(P-0.5)))
c      B1=8.256268+(-0.585350*(P-0.5))+(-1.901850*((P-0.5)*(P-0.5)))
c      B2=0.462866+(0.012206*(P-0.5))+(1.047697*((P-0.5)*(P-0.5)))
c      B11=19.788257+(-3.706511*(P-0.5))+(2.990894*((P-0.5)*(P-0.5)))
c      B22=0.539077+(-0.172307*(P-0.5))+(1.031416*((P-0.5)*(P-0.5)))
c      C=3.078636+(-1.598544*(P-0.5))+(-0.742071*((P-0.5)*(P-0.5)))
c      D=29.803617+(-22.124036*(P-0.5))+(-57.863105*((P-0.5)*(P-0.5)))
c      E1=3.212414+(0.804397*(P-0.5))+(0.228317*((P-0.5)*(P-0.5)))
c      EC=4.724383+(0.024499*(P-0.5))+(0.030653*((P-0.5)*(P-0.5)))
c      EIC=3.160138+(0.138046*(P-0.5))+(1.066214*((P-0.5)*(P-0.5)))
c      EINF=-0.494941+(-0.030802*(P-0.5))+(0.486201*((P-0.5)*(P-0.5)))
c      E1C=6.413109+(0.571119*(P-0.5))+(-0.735610*((P-0.5)*(P-0.5)))
c      T=0.263476+(0.090532*(P-0.5))+(-0.254099*((P-0.5)*(P-0.5)))
c      G=0.147517+(-0.068764*(P-0.5))+(0.047345*((P-0.5)*(P-0.5)))
c      F=1.628226+(-1.682422*(P-0.5))+(-2.081273*((P-0.5)*(P-0.5)))
c      FI=0.507707+(-0.070165*(P-0.5))+(-0.122169*((P-0.5)*(P-0.5)))
C
c      EO=1.425000+(1.155000*P)+(0.370000*(P*P))
c      EI=1.734000+(0.574000*P)+(0.055000*(P*P))
c      DO=0.340000+(0.0*P)+(0.0*(P*P))
c      D1=0.230000+(-0.030000*P)+(0.0*(P*P))
C
c      E=(4.135701327D-15*2.99792458D14)/WVL
c      EE=E*E
c      EEI=1.D0/EE
C
CCCCCC REAL PART OF DIELECTRIC FUNCTION : DIRECT EDGE (EO)
C
c      ARG1=DCMPLX(1.D0+E/EO,0.D0)
c      ARG2=DCMPLX(1.D0-E/EO,0.D0)
c      EP1EO1=A*EEI*DREAL(DSQRT(EO)*(2.D0-CDSQRT(ARG1)-CDSQRT(ARG2)))
c      ARG1=DCMPLX(1.D0+E/(EO+DO),0.D0)
c      ARG2=DCMPLX(1.D0-E/(EO+DO),0.D0)
c      EP1EO2=A*EEI*DREAL(0.5D0*DSQRT(EO+DO)*(2.D0-CDSQRT(ARG1)
c     +      -CDSQRT(ARG2)))    
c      EP1EO=EP1EO1+EP1EO2
C
CCCCCC IMAGINARY PART OF DIELECTRIC FUNCTION
C
c      ARG1=DCMPLX(E-EO,0.D0)
c      ARG2=DCMPLX(E-EO-DO,0.D0)
c      EP2EO=DREAL((2.D0*CDSQRT(ARG1)+CDSQRT(ARG2))*EEI)*A/3.D0  
C
CCCCCC REAL PART OF DIELECTRIC FUNCTION : L=[K=PI/A]
C
c      CEOO=E*DCMPLX(1.D0,0.D0)+T*DCMPLX(0.D0,1.D0)      
c      CEC=1.D0/(DCMPLX(1.D0,0.D0)-((CEOO/E1C)*(CEOO/E1C)))
c      CE=CEOO*CEOO
c      CEO=CE/((E1+D1)*(E1+D1))
c      CE=CE/(E1*E1)
c      CEXPO=MIN(DEXP(-F*(E-E1)),1.D0)/CE
c      CEXPOO=MIN(DEXP(-F*(E-E1-D1)),1.D0)/CEO
c      EP1E11=-B1*DREAL(LOG(((1.D0,0.D0)-CE)*CEC)*CEXPO)
c      EP1E12=-B2*DREAL(LOG(((1.D0,0.D0)-CEO)*CEC)*CEXPOO)
c      EP1E1=EP1E11+EP1E12
C
CCCCCC IMAGINARY PART
C
c      ARG1=DCMPLX(E1,0.D0)-CEOO
c      ARG2=DCMPLX(E1+D1,0.d0)-CEOO
c      EP2E11=PI*MAX(DREAL((B1-B11*CDSQRT(ARG1))*CEXPO),0.D0)
c      EP2E12=PI*MAX(DREAL((B2-B22*CDSQRT(ARG2))*CEXPOO),0.D0)
c      EP2E1=MAX(EP2E11+EP2E12,0.D0)
C
CCCCCC REAL PART
C
c      CE=DCMPLX(1.D0,0.D0)-EE/(EC*EC)
c      EP2EC=DREAL(C/(CE*CE+EE*((G/EC)*(G/EC))))
c      EP1EC=DREAL(CE*EP2EC)
C
CCCCCC IMAGINARY PART
C
c      EP2EC=(E*G/EC)*EP2EC
C
CCCCCC IMAGINARY PART
C
c      EQ1=E-EI
c      EP2EI=(D*EEI*(EQ1*EQ1))*MIN(DEXP(-FI*(E-EIC)),1.D0)
C
CCCCCC HEAVYSIDE FUNCTION
C
c      IF (E.LT.EO) THEN
c         EP2EO=0.D0
c         EP2E1=0.D0
c         EP2EC=0.D0
c         EP2EI=0.D0
c      ENDIF
c      IF (E.LT.E1) THEN
c         EP2E1=0.D0
c      ENDIF
C
c      EPSILON1=EINF+EP1EO+EP1E1+EP1EC
c      EPSILON2=EP2EO+EP2E1+EP2EC+EP2EI
c      ARG1=DCMPLX(EPSILON1,EPSILON2)
c      NEFF=DREAL(CDSQRT(ARG1))
c      PRINT*, ' N= ', NEFF
c      PRINT*, ' FOR NEW CALCULATION, INPUT 1, FOR STOP INPUT 2 '
c      PRINT*, ' INPUT =?'
c      READ(*,*) I
c      IF (I.EQ.1) THEN
c      GOTO 1000
c      ELSE
c      GOTO 2000
c      ENDIF
c 2000 PRINT*, ' PROGRAM STOP HERE!'
c      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sAlGaInP
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' (AlxGa(1-x))0.5In0.5P--Eg=1.9+0.6*X'
      PRINT*, ' INPUT X AND PHOTON WAVELENGTH FOR AlGaInP'
      PRINT*, ' X= ?,  LAMBDA='
      READ(*,*) X,XLAMBDA
      HV=1.24/XLAMBDA
      EO=3.39+0.62*X
      ED=28.07+1.72*X
      XN=DSQRT((EO*ED/(EO**2-HV**2))+1)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sInGaAs1
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' In(1-x)Ga(x)As/InP'
      PRINT*, ' INPUT X AND PHOTON WAVELENGTH FOR In(1-x)Ga(x)As'
      PRINT*, ' X= ?,  LAMBDA='
      READ(*,*) X,XLAMBDA
      A=(1-X)*5.14+X*6.30
      B=(1-X)*10.15+X*9.40
      HV=1.24/XLAMBDA
      EO=0.324+0.7*X+0.4*(X**2)
      EDO=(1-X)*0.38+X*0.34
      EDD=EDO+EO
      XO=HV/EO
      XSO=HV/EDD
	Q1=MAX(1-XO,0.D0)
	Q2=MAX(1-XSO,0.D0)
	F1=(XO**(-2))*(2-((1+XO)**0.5)-(Q1**0.5))
	F2=(XSO**(-2))*(2-((1+XSO)**0.5)-(Q2**0.5))
C      F1=(XO**(-2))*(2-((1+XO)**0.5)-((1-XO)**0.5))
C      F2=(XSO**(-2))*(2-((1+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/(EDD))**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sInGaAs2
      IMPLICIT REAL*8 (A-H,O-Z)
   10 PRINT*, ' In(1-x)Ga(x)As/GaAs'
      PRINT*, ' INPUT X AND PHOTON WAVELENGTH FOR In(1-x)Ga(x)As'
      PRINT*, ' X= ?,  LAMBDA='
      READ(*,*) X,XLAMBDA
      A=(1-X)*5.14+X*6.30
      B=(1-X)*10.15+X*9.40
      HV=1.24/XLAMBDA
      EO=0.36+0.509*X+0.555*(X**2)
      EDO=(1-X)*0.38+X*0.34
      EDD=EDO+EO
      XO=HV/EO
      XSO=HV/EDD
	Q1=MAX(1-XO,0.D0)
	Q2=MAX(1-XSO,0.D0)
	F1=(XO**(-2))*(2-((1+XO)**0.5)-(Q1**0.5))
	F2=(XSO**(-2))*(2-((1+XSO)**0.5)-(Q2**0.5))
C      F1=(XO**(-2))*(2-((1+XO)**0.5)-((1-XO)**0.5))
C      F2=(XSO**(-2))*(2-((1+XSO)**0.5)-((1-XSO)**0.5))
      X1=((EO/(EDD))**(1.5))/2
      X2=F1+X1*F2
      X3=A*X2
      X4=X3+B
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE sInGaAsN
      IMPLICIT REAL*8 (A-H,O-Z)
   10	PRINT*, ' DUE TO INSUFFICIENT DATA, THE WAVELENGTH DEPENDENT'
	PRINT*, '                                                  '
	PRINT*, ' REFRACTIVE INDEX CAN NOT CALCULATE AT THIS TIME!'
	PRINT*, '                                                 '
      PRINT*, ' In(y)Ga(1-y)As(x)N(1-x)/GaAs WITH SIMPLE INTERPOLATION.'
      PRINT*, '                                                       '
      PRINT*, ' INPUT X and Y FOR In(y)Ga(1-y)As(x)N(1-x)'
      PRINT*, ' X= ?,  Y='
      READ(*,*) X,Y
      X4=(1-Y)*X*12.85+X*Y*15.15+(1-Y)*(1-X)*10+Y*(1-X)*15.3
      XN=DSQRT(X4)
      PRINT*, ' REFRACTIVE INDEX=', XN
      PRINT*, ' FOR NEW INPUT I=1, STOP I=2'
      PRINT*, ' I=' 
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 10
      ELSEIF (I.EQ.2)  THEN
      GOTO 20
      ELSE
      ENDIF
   20 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SELECT2 (N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
      WRITE(*,11)
   11 FORMAT(/' CALCULATE THE WAVE LENGTH DEPENDENT REFRACTIVE INDEX'/
     +        ' ENTER 1 FOR In(1-X)Ga(X)As(Y)P(1-Y), Adachi formula'/
     +        '       2 FOR In(1-X)Ga(X)As(Y)P(1-Y), Non-Adachi '/
     +        '       3 FOR Al(X)Ga(Y)In(1-X-Y)As  Mondry formula'/
     +        '       4 FOR Al(X)Ga(Y)In(1-X-Y)As Adachi formula'/
     +        '       5 FOR Al(X)Ga(1-X)P(Y)Sb(1-Y)'/
     +        '       6 FOR Al(X)Ga(1-X)As(Y)Sb(1-Y)'/
     +        '       7 FOR Al(X)In(1-X)As(Y)Sb(1-Y)'/
     +        '       8 FOR Ga(X)In(1-X)P(Y)Sb(1-Y)'/ 
     +        '       9 FOR AlAs(X)Sb(1-X)         '/
     +        '      10 FOR Al(x)Ga(1-x)As         '/
     +        '      11 FOR Al(x)Ga(1-x)As Jenkins formula '/
     +        '      12 FOR (AlxGa(1-x))0.5P0.5    '/
     +        '      13 FOR In(1-x)Ga(x)As/(matched to InP)        '/
     +        '      14 FOR In(1-x)Ga(x)As/(matched to GaAs)        '/
     +        '      15 FOR In(y)Ga(1-y)As(x)N(1-x)/(matched to GaAs'/
     +        '      16 FOR STOP')
      READ(*,*) N
      RETURN
      END      


C                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE COUPLE 
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(500),B(500),H(500),AAX(500),BBX(500),HHX(500)
      COMPLEX*16 BX,CI,E(500),AUP,BUP,HUP,ABELOW,BBELOW,HBELOW,BETAA
      COMPLEX*16 CXUP,CXBELOW,C1,KAB,BETA
      REAL*8 KO,DX(500),LAMBDA,LAMBDAG,N(500),NUP,NBELOW,DUP
      INTEGER NL,IC,LUP,LBELOW
      COMMON DX,E,NL,IC,CI
      COMMON /A5/ AAX,BBX,HHX
      COMMON /STATE/ BX,KO
      COMMON /RUP/ AUP,BUP,HUP,CXUP
      COMMON /RBELOW/ ABELOW,BBELOW,HBELOW,CXBELOW  
      COMMON /A7/ LUP,LBELOW
      COMMON /A1/ N     
      COMMON /A2/ C1          
      COMMON /PERIOD/ G,GL,M   
      COMMON /A8/ DUP,DMID
      COMMON /A10/ LAMBDA,BR,BI
      COMMON /XCOUPLE/ KAB
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      PRINT*, 'C                                                      C'
      PRINT*, 'C THIS SUBROUTINE FINDS THE COUPLING COEFFICIENT OF    C'
      PRINT*, 'C DFB LASER STRUCTURE.                                 C'
      PRINT*, 'C                                                      C'
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      CALL INDEX1
      PRINT*, ' INPUT WAVELENGTH= (um)?'
      PRINT*, ' LAMBDA=',LAMBDA
      PI=3.1415926535897932D0
      C=3.0D8
      PRINT*, ' INPUT THE EFFECTIVE REFRACTIVE INDEX'
      PRINT*, ' (BETA IS COMPLEX)='
      BETA=DCMPLX(BR,BI)
      PRINT*, ' BETA=',BETA
      BETAA=BETA*KO  
	PRINT*, '                               '
	PRINT*, ' N-SUB IS THE FIRST LAYER'
	PRINT*, '                               '
c      PRINT*, ' INPUT THE LAYER ABOVE AND BELOW THE GRATING' ! N-SUB 1st layer.
c      PRINT*, ' LAYER ABOVE= ? , LAYER BELOW= ?'
c      READ(*,*) LUP,LBELOW
      PRINT*, ' INPUT THE NUMBER OF GRATING LAYER'
	PRINT*, ' LAYER=?'
	READ(*,*) XGRATING 
      PRINT*, ' INPUT THE INDICE WHICH FORMED GRATING LAYER'
      PRINT*, ' N1=?, N2=?'
      READ(*,*) NUP,NBELOW
      DO I=1,NL       
         A(I)=AAX(I)     
         B(I)=BBX(I)    
         H(I)=HHX(I)
      ENDDO             
	LUP=XGRATING
	LBELOW=XGRATING     
      AUP=A(LUP)
      BUP=B(LUP)
c     NUP=N(LUP)
c     NBELOW=N(LBELOW)
      ABELOW=A(LBELOW)
      BBELOW=B(LBELOW) 
      HUP=H(LUP)           
      HBELOW=H(LBELOW)         
 9997 PRINT*, ' INPUT THE TOOTH HEIGHT G= (um), BRAGG CONDITION M=?' 
      READ(*,*) G,M
      PRINT*, ' INPUT RATIO OF (G) FOR UPPER(NU) AND LOWER(NL) LAYER?'
	PRINT*, ' FOR GENERAL CASE THE RATIO FOR BOTH ARE 0.5' 
      PRINT*, ' NU=? , NL=?'
      READ(*,*) XUG,XBG
      G=-G
      GU=XUG*G
      GB=XBG*G
      LAMBDAG=LAMBDA/DREAL(BETA)
      GL=M*LAMBDAG/2.D0
      DUP=DX(LUP)
      DBELOW=DUP
	DMID=DUP-(G/2)
      DXX=DMID+GU
      DYY=DMID-GB
      C1=(KO**2*(NBELOW**2-NUP**2))/(4*PI*BETAA*M)
      PRINT*, ' GRATING PERIOD GL = ?',GL,' um'
      PRINT*, '                                     '
      PRINT*, ' INPUT 1 FOR RECTANGULAR GRATING          ' 
      PRINT*, ' INPUT 2 FOR SINUSOIDAL GRATING           '
      PRINT*, ' INPUT 3 FOR SYMMETRIC TRIANGULAR GRATING ' 
      PRINT*, ' INPUT 4 FOR SAW TOOTH GRATING            '
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
         PRINT*, ' INPUT THE BOUNDARY VARIABLE W1,W2,W3,W4'
         PRINT*, ' W1= ?,W2= ?,W3= ?,W4= ?'
         READ(*,*) W1,W2,W3,W4
         CXUP=(CDEXP(CI*2*PI*M*W2)-CDEXP(CI*2*PI*M*W1))
         CXBELOW=(CDEXP(CI*2*PI*M*W4)-CDEXP(CI*2*PI*M*W3))
         CALL RECT(DXX,DYY,KAB)
       ELSEIF (I.EQ.2) THEN
         CALL SINUS (DXX,DYY,KAB)
       ELSEIF (I.EQ.3) THEN
          CALL TRI (DXX,DYY,KAB)
       ELSEIF (I.EQ.4) THEN
          CALL SAW (DXX,DYY,KAB)
       ELSE
      ENDIF
      PRINT*, ' KAB=',KAB,' 1/um' 
      PRINT*, ' |KAB|=',CDABS(KAB)
      AKAB=CDABS(KAB)
      PRINT*, ' FOR NEW INPUT DATAS INPUT 1 FOR MAIN SCREEN  INPUT 2'
      PRINT*, ' INPUT = '
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 9997
      ELSE
      GOTO 9999
      ENDIF
 9999 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE DFB
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 BPP(200),BPN(200),ZS,AX(200),BX(200),DELTAB
      COMPLEX*16 ZG,BETAB(200),DELTABB,R(200),T(200)
      COMPLEX*16 CI,E(500),BETAA,CDSINH,CDCOSH,CDSINHX,CDCOSHX,ZS1
      COMPLEX*16 KAB,BETA,BETAZ,ZN(500),EX(500)
      REAL*8 KO,DX(500),LAMBDA,LAMBDAG,N(500),ALPHA(500),W(500),K(500)
      REAL*8 KOX(200),D(500)
      INTEGER NL,IC,L
      CHARACTER IFILE*40 
      COMMON DX,E,NL,IC,CI
      COMMON /A10/ LAMBDA,BR,BI
      COMMON /XCOUPLE/ KAB
	COMMON /A11/ AIN
      EXTERNAL TE                            
      PI=3.1415926535897932D0 
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      PRINT*, 'C                                                    C'
      PRINT*, 'C FIND THE W-B, REFLECTION AND TRANSFER PLOT         C'
      PRINT*, 'C                                                    C'
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' 
 9997 PRINT*, ' # OF Lth MODE L= ?'
      READ(*,*) L
      REWIND (1)
      REWIND (20)
      REWIND (23)
      BETA=DCMPLX(BR,BI)
      PRINT*, ' BETA=',BETA
      PRINT*, ' PRINT COUPLEING COEFFICIENT KAB(COMPLEX) = (1/um)?'
      PRINT*, ' KAB=',KAB,' (1/um)'
      AKAB=CDABS(KAB)
      KO=2*PI/LAMBDA
      LAMBDAG=LAMBDA/DREAL(BETA)
      GLL=L*LAMBDAG/2.D0
      PRINT*, 'GLL=',GLL, ' AKAB=',AKAB
      BETAA=BETA*KO       
      DELTAB=BETAA-(L*PI/GLL)
      PRINT*, ' DELTAB=',DELTAB
      ZS1=CDSQRT((KAB*DCONJG(KAB))-DELTAB*DELTAB)
      PRINT*, ' CALCULATE THE W-B PLOT' 
      MIN=5
      PRINT*,' THE INPUT FILE NAME='
      READ(MIN,'(A40)') IFILE
      OPEN(1,FILE=IFILE,status='unknown')
      READ(1,*) NL    
      READ(1,*) LAMBDA
      DO I=1,NL
      READ(1,*) N(I),ALPHA(I),W(I)
      ENDDO              
      OPEN(20,FILE='wb.dat',STATUS='UNKNOWN')
c     OPEN(22,FILE='ampl.dat',STATUS='UNKNOWN')
      OPEN(23,FILE='tranref.dat',STATUS='UNKNOWN')
      PRINT*, ' INPUT THE LENGTH OF GRATING XL= (um)?'
      READ(*,*) XL
c     PRINT*, ' INPUT THE CENTER LAYER IC = '
c     READ(*,*) IC 
      TPI=2.D0*PI
      KO=TPI/(LAMBDA)
      EPO=8.85418782D-12
      UO=4.*PI*1.D-7
      C=3.D8                                     
      CI=DCMPLX(0.D0,1.D0)
      AIN=1.D-4
      EPS=1.D-6
      EPSS=1.D-7
      DO I=1,150 
         XLAMBDA=(LAMBDA-0.01D0)+(I-1)*
     +   (((LAMBDA+0.01D0)-(LAMBDA-0.01D0))/149)
         KOX(I)=2*PI/(XLAMBDA)
         DO J=1,NL
            K(J)=ALPHA(J)/(KOX(I))
            ZN(J)=DCMPLX(N(J)*KOX(I),K(J)*KOX(I))
            E(J)=ZN(J)*ZN(J)
            D(J)=W(J)                         ! INPUT THE DISTANCE
         ENDDO
         DX(IC-1)=D(IC)/2.D0                  ! FIND THE CENTER LAYER
         DX(IC)=(-D(IC))/2.D0
         DO J=IC-2,1,-1                       ! FROM LEFT TO CENTER
            DX(J)=DX(J+1)+D(J+1)
         ENDDO
         DO J=IC+1,NL-1,1                     ! FROM RIGHT TO CENTER
            DX(J)=DX(J-1)-D(J)      
         ENDDO
         ZG=DCMPLX(BR*KOX(I),BI*KOX(I)) 
c        CALL CROOT (TE,ZG,ITTER,J,EPS,EPSS)
c        BETAZ=ZG/KO 
c        BETAB(I)=ZG/KO  
	 BETAZ=ZG
	 BETAB(I)=ZG
c        BPP(I)=(L*PI/GLL)+CDSQRT((BETAZ*(2*PI/XLAMBDA)-
c    +          (L*PI/GLL))**2-AKAB*AKAB)    
	 BPP(I)=(L*PI/GLL)+CDSQRT((BETAZ-(L*PI/GLL))**2-AKAB*AKAB)
c        BPN(I)=(L*PI/GLL)-CDSQRT((BETAZ*(2*PI/XLAMBDA)-
c    +          (L*PI/GLL))**2-AKAB*AKAB) 
	 BPN(I)=(L*PI/GLL)-CDSQRT((BETAZ-(L*PI/GLL))**2-AKAB*AKAB)
c        DELTABB=(BETAB(I)*(2*PI/XLAMBDA))-(L*PI/GLL)  
	 DELTABB=BETAB(I)-(L*PI/GLL)
         ZS=CDSQRT((KAB*DCONJG(KAB))-DELTABB*DELTABB)    
c         XDELTAB=DREAL(DELTABB/AKAB)
         XDELTAB=DREAL(DELTABB)   
         CDSINH=(CDEXP(ZS*XL)-CDEXP(-(ZS*XL)))/2.D0 
         CDCOSH=(CDEXP(ZS*XL)+CDEXP(-(ZS*XL)))/2.D0   
         R(I)=(-DCONJG(KAB)*CDSINH)/(DELTABB*CDSINH
     +   +CI*ZS*CDCOSH)  
         T(I)=(CI*ZS)/(DELTABB*CDSINH+CI*ZS*CDCOSH)    
      WRITE(20,21) (2*PI/XLAMBDA)*GLL,BPP(I),BPN(I)  
      WRITE(23,24) XDELTAB,CDABS(R(I))**2,CDABS(T(I))**2
   21 FORMAT(1X,E10.5,1X,2(1X,E12.6),1X,2(1X,E12.6))
   24 FORMAT(2X,E12.5,2X,E12.6,2X,E12.6)
      ENDDO
      CLOSE (20)
c     PRINT*, ' CALCULATE THE A(Z) AND B(Z) WAVE FUNCTION'
c     XXXL=0.D0
c     DO I=1,150    
c        XXL=XXXL+(I-1)*(XL/149)     
c        CDSINHX=(CDEXP(ZS1*(XXL-XL))-CDEXP(-(ZS1*(XXL-XL))))/2.D0   
c        CDCOSHX=(CDEXP(ZS1*(XXL-XL))+CDEXP(-(ZS1*(XXL-XL))))/2.D0 
c        CDSINH=(CDEXP(ZS1*XL)-CDEXP(-(ZS1*XL)))/2.D0  
c        CDCOSH=(CDEXP(ZS1*XL)+CDEXP(-(ZS1*XL)))/2.D0   
c        AX(I)=(-DELTAB*CDSINHX+CI*ZS1*CDCOSHX)
c    +   /(DELTAB*CDSINH+CI*ZS1*CDCOSH)
c        BX(I)=DCONJG(KAB)*CDSINHX           
c    +   /(DELTAB*CDSINH+CI*ZS1*CDCOSH)     
c        WRITE(22,*) XXL,CDABS(AX(I)),CDABS(BX(I))
c     ENDDO
c     CLOSE (22)
      CLOSE (23)
      CLOSE (1)
      PRINT*, ' FOR NEW INPUT DATAS INPUT 1 FOR MAIN SCREEN  INPUT 2'
      PRINT*, ' INPUT = '
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 9997
      ELSE
      GOTO 9999
      ENDIF
 9999 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine dfbsl
      implicit real*8 (a-h,o-z)
      complex*16 z,w4,r1,r2,ci,kab
      real*8  xy,xx,za,xl,ko,lambdag,lambda,deltabb(30)
      integer icase
      COMMON /A10/ LAMBDA,BR,BI
      COMMON /XCOUPLE/ KAB
      COMMON /INX/ IXX
      common /a6/ r1,r2,xl,xkab
      common /d1/ xx,icase
      external wx
      open(3,file='gain.dat',status='unknown')
      open(4,file='gain2.dat',status='unknown')
      PRINT*, '                                                      '
      PRINT*, ' INPUT # OF Lth MODE OF DFB STRUCTURE ?'
      READ(*,*) L
      PRINT*, '                                                      '
      PRINT*, ' INPUT THE REFELECTIVITY OF BOUNDARIES R1 AND R2 = ? '
      PRINT*, ' R1=?,R2=?,PHASE -->PHI1=?,PHI2=? '
      READ(*,*) XR1,XR2,PHI1,PHI2
      PRINT*, '                                                      '
      PRINT*, ' INPUT THE CAVITY LENGTH XL OF LASER,XL= (um)?'
      READ(*,*) XL
      PRINT*, '                                                      '
      ci=dcmplx(0.d0,1.d0)
      eps=1.e-9
      epss=1.e-11
      r1=xr1*cdexp(ci*phi1)
      r2=xr2*cdexp(ci*phi2)
      AKAB=CDABS(KAB)
      XKAB=DREAL(KAB)*XL
      PI=3.1415926535897932D0
      KO=2*PI/LAMBDA
      BETA=BR
      LAMBDAG=LAMBDA/BETA
      GL=L*LAMBDAG/2.D0
      DO I=1,30
C         XLAMBDA=(LAMBDA-0.05D0)+(I-1)*
C     +   (((LAMBDA+0.05D0)-(LAMBDA-0.05D0))/29)
         XLAMBDA=LAMBDA-0.001*(I-1)
c         KOX(I)=2*PI/(XLAMBDA)
         DELTABB(I)=((BETA*(2*PI/XLAMBDA))-(L*PI/GL))*XL  !DETUNG
      ENDDO
    1 print*, ' FOR R1=R2=0 NO END REFLECTION, CASE1 I=1'
      print*, ' FOR R1,R2,K NOT ZERO, THEN CASE2 I=2'
      print*, ' FOR K=0, NO DISTRIBUTED FEED BACK, CASE3 I=3'
      print*, ' CASE= ?'
      read(*,*) icase
      print*, '                                             '
      print*, ' INPUT ZA (initial guess)'
      READ(*,*) ZA
      ix=1
c  12 xx=DELTABB(IX) 
c  12 xx=(1.+((ix-1)/2.)*1.)*(-1.d0)
   12 xx=-15.d0+((ix-1)/2.d0)*1.d0
      if (xx.lt.0.d0) then
      ixx=1
      elseif (xx.ge.0.d0) then
      ixx=2
      else
      endif
      call xcroot (wx,z,za,w4,i,j,eps,epss)
c     if ((dimag(z).gt.-1.d0).and.(dimag(z).lt.0.d0)) then
      if (i.eq.2) then
      write(3,11) z,i,j,xx
      write(4,16) xx,z
   11 format(2x,'z=',2(d14.7),2x,'i=',i6,2x,' j=',I4,2x,' xx=',f7.3)
   16 format(2x,f7.3,2x,2(d14.7))
      else
      endif
      if (ix.lt.61) then
      ix=ix+1
      za=dreal(z)
      goto 12
      elseif (ix.gt.61) then
      goto 14
      else
      endif
      close (3)
      close (4)
   14 print*, ' FOR I=1, INPUT NEW VALUES, I=2, STOP'
      print*, ' INPUT I=?'
      read(*,*) i
      if (i.eq.1) then
      goto 1
      else
      goto 2
      endif
    2 return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SYNCHRO
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A1(500),B1(500),H1(500),A2(500),B2(500),H2(500)
      COMPLEX*16 BBX(500),HHX(500),E(500),K12,K21,K11,K22,CCOEF3
      COMPLEX*16 CCOEF1,CCOEF2,BX,H(500),A(500),B(500),KAPA1,KAPA2
      COMPLEX*16 CI,KAPA,NA,NB,BETA1,BETA2,BPP(150),BPM(150),ZDELTA
      COMPLEX*16 ZDELTA2,BNA(150),BNB(150),ZN(500),ZG,AAX(500),CCOEF4
      REAL*8 KO,DX(500),DIST,NG1,NG2,NC,KOX(150),LAMBDA,ALPHA(500)
      REAL*8 D(500),W(500),K(500),N(500),KBA,KAB,AP(91),BP(91)
      INTEGER NL,IC
      CHARACTER IFILE*40
      COMMON DX,E,NL,IC,CI
      COMMON /A5/ AAX,BBX,HHX
      COMMON /AC1/ A1,B1,H1
      COMMON /AC2/ A2,B2,H2
      COMMON /AC4/ DIST 
      COMMON /FAR1/ H,A,B,XNORN
      COMMON /STATE/ BX,KO
      COMMON /A10/ LAMBDA,BR,BI
      COMMON /A11/ AIN
      EXTERNAL TE
      OPEN(32,FILE='ptran.dat',status='unknown') 
      PI=3.1415926535897932D0
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      PRINT*, 'C                                                     C'
      PRINT*, 'C THIS SUBROUTINE FIND THE POWER EXCHANGE LENGTH LX   C'
      PRINT*, 'C W-B PLOT, AND POWER TRANSFER OF SYNCHRONOUS CASE.   C'
      PRINT*, 'C                                                     C'
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      PRINT*, '                                                       '
      OPEN(30,FILE='wbs.dat',status='unknown')
 9997 PRINT*, ' FIND THE EFFECTIVE INDEX Na AND COEFF. OF 1st WAVWGUIDE'
      CALL INDEX1
      BR1=BR
      BI1=BI
      DO I=1,NL    
         A1(I)=AAX(I)     
         B1(I)=BBX(I)    
         H1(I)=HHX(I)
      ENDDO          
c     DXA=DX(1)                
      PRINT*, ' FIND THE EFFECTIVE INDEX Nb AND COEFF. OF 2nd WAVWGUIDE'
      CALL INDEX1
      BR2=BR
      BI2=BI
      DO I=1,NL      
         A2(I)=AAX(I)   
         B2(I)=BBX(I)    
         H2(I)=HHX(I)
      ENDDO       
c     DXB=DX(1)   
c     DIST=(DXB-DXA)  
      PRINT*, ' INPUT THE EFFECTIVE REFRACTIVE INDEX OF GUIDE 1 AND 2'
      PRINT*, ' NA= ?, NB= ?'
      READ(*,*) NA,NB                     
      BETA1=NA*KO
      BETA2=NB*KO
      PRINT*, '                                                     '
      PRINT*, ' INPUT THE INDEX OF WG1(NG1),WG2(NG2) AND CORE(NC)'
      READ(*,*) NG1,NG2,NC
      KAPA=(2*KO*3*1.D14*4*PI*1.D-13)/CDSQRT(BETA1*BETA2)
      KAPA1=(2*KO*3*1.D14*4*PI*1.D-13)/CDSQRT(BETA1*BETA1)
      KAPA2=(2*KO*3*1.D14*4*PI*1.D-13)/CDSQRT(BETA2*BETA2)
	print*, ' kapa=',kapa
      CALL COUPL1 (NG1,NG2,NC,CCOEF1,CCOEF2,CCOEF3,CCOEF4)      
      WRITE(*,11) CCOEF1*KAPA,CCOEF2*KAPA
      K12=CCOEF1*KAPA
      K21=CCOEF2*KAPA
   11 FORMAT(1X,'K12=',2(E15.7),1X,' K21=',2(E15.7))
      PL=PI/(2.D0*CDABS(K12))
C      PL=PI/DABS(BETA1-BETA2)
      PRINT*, '                                                     '
      PRINT*, ' THE POWER EXCHANGE LENGTH = ',PL, ' um'
      PRINT*, '                                                     '
      PRINT*, ' FOR W-B PLOT, THAN INPUT 1, OTHERWISE INPUT 2 '
      PRINT*, ' INPUT =?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 7000
      ELSE
      GOTO 8000
      ENDIF
 7000 PRINT*, '                                                     '
      PRINT*, ' CALCULATE THE W DEPENDENT BETA VALUE OF WG1'
      PRINT*, '                                                     '
      MIN=5
      PRINT*,' THE INPUT FILE NAME='
      READ(MIN,'(A40)') IFILE
      OPEN(1,FILE=IFILE,status='unknown')
      READ(1,*) NL
      READ(1,*) LAMBDA
      DO I=1,NL
      READ(1,*) N(I),ALPHA(I),W(I)
      ENDDO
      AIN=1.D-5
      EPS=1.D-8
      EPSS=1.D-9
      DO I=1,150
         XLAMBDA=(LAMBDA-0.01D0)+(I-1)*
     +   (((LAMBDA+0.01D0)-(LAMBDA-0.01D0))/149)
         KOX(I)=2*PI/(XLAMBDA)
         DO J=1,NL
            K(J)=ALPHA(J)/(2.D0*KOX(I))
            ZN(J)=DCMPLX(N(J)*KOX(I),K(J)*KOX(I))
            E(J)=ZN(J)*ZN(J)
            D(J)=W(J)                         ! INPUT THE DISTANCE
         ENDDO
         DX(IC-1)=D(IC)/2.D0                  ! FIND THE CENTER LAYER
         DX(IC)=(-D(IC))/2.D0
         DO J=IC-2,1,-1                       ! FROM LEFT TO CENTER
            DX(J)=DX(J+1)+D(J+1)
         ENDDO
         DO J=IC+1,NL-1,1                     ! FROM RIGHT TO CENTER
            DX(J)=DX(J-1)-D(J)
         ENDDO
         ZG=DCMPLX(BR1*KOX(I),BI1*KOX(I))
         CALL CROOT (TE,ZG,ITTER,J,EPS,EPSS)
         BNA(I)=ZG
	   BR1=DREAL(BNA(I)/KOX(I))
	   BI1=DIMAG(BNA(I)/KOX(I))
      ENDDO
      CLOSE (1)
      PRINT*, '                                            '
      PRINT*, ' CALCULATE THE W DEPENDENT BETA VALUE OF WG2'
      PRINT*, '                                            '
      MIN=5
      PRINT*,' THE INPUT FILE NAME='
      READ(MIN,'(A40)') IFILE
      OPEN(1,FILE=IFILE,status='unknown')
      READ(1,*) NL
      READ(1,*) LAMBDA
      DO I=1,NL
      READ(1,*) N(I),ALPHA(I),W(I)
      ENDDO
      DO I=1,150
         DO J=1,NL
            K(J)=ALPHA(J)/(2.D0*KOX(I))
            ZN(J)=DCMPLX(N(J)*KOX(I),K(J)*KOX(I))
            E(J)=ZN(J)*ZN(J)
            D(J)=W(J)                         ! INPUT THE DISTANCE
         ENDDO
         DX(IC-1)=D(IC)/2.D0                  ! FIND THE CENTER LAYER
         DX(IC)=(-D(IC))/2.D0
         DO J=IC-2,1,-1                       ! FROM LEFT TO CENTER
            DX(J)=DX(J+1)+D(J+1)
         ENDDO
         DO J=IC+1,NL-1,1                     ! FROM RIGHT TO CENTER
            DX(J)=DX(J-1)-D(J)
         ENDDO
         ZG=DCMPLX(BR2*KOX(I),BI2*KOX(I))
         CALL CROOT (TE,ZG,ITTER,J,EPS,EPSS)
         BNB(I)=ZG
	   BR2=DREAL(BNB(I)/KOX(I))
	   BI2=DIMAG(BNB(I)/KOX(I))
      ENDDO
      CLOSE (1)
	! BETA1~BNA(I) 
      DO I=1,150  
         ZDELTA=((BNA(I)-bnb(i))/2.D0)**2
         ZDELTA2=(BNA(I)+bnb(i))/2.D0
         BPP(I)=ZDELTA2+CDSQRT(ZDELTA+K12*K21)
         BPM(I)=ZDELTA2-CDSQRT(ZDELTA+K12*K21) 
         WRITE(30,31) KOX(I),BPP(I),BPM(I),BETA1*KOX(I),BETA2*KOX(I)
   31  FORMAT(2X,E12.6,2(1X,E12.6),2(1X,E12.6),2(1X,E12.6),2(1X,E12.6))
      ENDDO
      CLOSE (30) 
 8000 PRINT*, '                                                       '
      PRINT*, ' CALCULATE THE DATA FOR POWER TRANSFER                 '
      PRINT*, '                                                       '
      BP1=DREAL(NA)*KO
      BP2=DREAL(NB)*KO
      KAB=CDABS(K12)
      KBA=CDABS(K21)
      DELTA=(BP1-BP2)/2.D0
      PHASE=DSQRT(DELTA**2+4*KAB*KBA)
      DO I=1,181                                
         Z=((2*PI/90)*(I-1))/PHASE
         AP(I)=(DELTA**2/(PHASE**2))+(4*KAB*KBA/(PHASE**2))
     +   *(DCOS((PHASE/2)*Z))**2
         BP(I)=(4*KAB*KBA/(PHASE**2))*(DSIN((PHASE/2)*Z))**2
         WRITE(32,33) Z,AP(I),BP(I)
   33    FORMAT(1X,E12.6,1X,E14.7,1X,E14.7)
      ENDDO  
      PRINT*, ' FOR NEW INPUT DATAS INPUT 1 FOR MAIN SCREEN  INPUT 2'
      PRINT*, ' INPUT = '
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 9997
      ELSE
      GOTO 9999
      ENDIF
 9999 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE NONSYNCH
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A1(500),B1(500),H1(500),A2(500),B2(500),H2(500)
      COMPLEX*16 BBX(500),HHX(500),E(500),CCOEF,CI,BETA1,BETA2,AAX(500)
      COMPLEX*16 BX,NA,NB
      REAL*8 DX(500),NG1,NC,GA,KO,KAPA,A(91),B(91),KAB,KBA 
      INTEGER NL,IC,IX
      COMMON DX,E,NL,IC,CI
      COMMON /A5/ AAX,BBX,HHX
      COMMON /AC1/ A1,B1,H1
      COMMON /AC2/ A2,B2,H2
      COMMON /GRATINGA/ GA    
      COMMON /STATE/ BX,KO
      COMMON /A9/ IX
      OPEN(30,FILE='ptrang.dat',status='unknown')
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      PRINT*, 'C                                                     C'
      PRINT*, 'C THIS SUBROUTINE FIND THE COUPLE COEFFICIENT AND     C'
      PRINT*, 'C POWER EXCHANGE LENGTH, AND POWER TRANSFER OF TWO    C'
      PRINT*, 'C COUPLED WAVEGUIDE                                   C'
      PRINT*, 'C (NONSYNCHRONOUS CASE).                              C'
      PRINT*, 'C                                                     C'
      PRINT*, 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      PRINT*, '                                                       '
      PRINT*, '*******************************************************'
      PRINT*, ' TAKE THE SECOND WG AS THE CENTER OF THE LAYER IC      '
      PRINT*, '*******************************************************' 
      PRINT*, '                                                       '  
 9997 PRINT*, ' CALL INDEX1 FIND THE FIRST MODE OF COMPOUND STRUCTURE '  
      CALL INDEX1
      DO I=1,NL    
         A1(I)=AAX(I)     
         B1(I)=BBX(I)    
         H1(I)=HHX(I)
      ENDDO          
      PRINT*, ' CALL INDEX1 FIND THE SECOND MODE OF COMPOUND STRUCTURE'  
      CALL INDEX1
      DO I=1,NL
         A2(I)=AAX(I)
         B2(I)=BBX(I)
         H2(I)=HHX(I)
      ENDDO               
      PRINT*, ' INPUT EFFECTIVE REFRACTIVE INDEX OF COMPOUND MODES'
      PRINT*, ' NA= ? , NB= ?'
      READ(*,*) NA,NB
      BETA1=NA*KO
      BETA2=NB*KO  
      PI=3.1415926535897932D0                      
      W=KO*3.D14
      U=4*PI*1.D-13
      KAPA=(2.*W*U)/(CDSQRT(BETA1*BETA2))
      PRINT*, ' FIND THE COUPLE COEFFICIENT(WITH GRATING) FOR COUPLER'
      PRINT*, ' INPUT THE GRATING WIDTH GA= ? um'
      READ(*,*) GA                                                     
      PRINT*, ' INPUT THE UPPER PORTION OF THE GRATING '
      PRINT*, ' INPUT LAYER= ? '
      READ(*,*) IX
      PRINT*, ' INPUT THE INDEX OF WAVE GUIDE(NG) AND CORE(NC) REGION'
      PRINT*, ' NG1= ?  NC= ?'  
      READ(*,*) NG1,NC 
      CALL COUPLE2 (NG1,NC,CCOEF)   
      CCOEF=(2.D0/PI)*GA*CCOEF*KAPA    
      PRINT*, ' K= ',CCOEF,' |K=|',CDABS(CCOEF),' 1/um'
      PL=PI/(2.D0*(CDABS(CCOEF)))
      PRINT*, '                                                       '
      PRINT*, ' POWER EXCHANGE LENGTH = ',PL*1.D-3,' mm'
      PRINT*, '                                                       '
      PRINT*, ' CALCULATE THE DATA FOR POWER TRANSFER                 '
      PRINT*, '                                                       '
c     BP1=DREAL(NA)*KO
c     BP2=DREAL(NB)*KO
      PRINT*, ' INPUT THE MODES OF SEPARATE SLAB WAVEGUIDES' 
      PRINT*, '                                                       '
      PRINT*, ' BP1=?, BP2=?'
      READ(*,*) BP1,BP2
      KAB=CDABS(CCOEF)
      delta=(BP1-BP2)/2
      PHASE=DSQRT(DELTA*DELTA+KAB*KAB)
      DO I=1,91
         Z=((2*PI/90)*(I-1))/PHASE
c        A(I)=(DELTA**2/(PHASE**2))+(4*KAB*KAB/(PHASE**2))
c    +   *(DCOS((PHASE/2)*Z))**2
	 A(I)=(DCOS(PHASE*Z))**2+((DELTA/PHASE)*DSIN(PHASE*Z))**2
         B(I)=(KAB*KAB/(PHASE**2))*(DSIN(PHASE*Z))**2
         WRITE(30,31) Z,A(I),B(I)
   31    FORMAT(1X,E12.6,1X,E14.7,1X,E14.7)
      ENDDO 
      PRINT*, ' FOR NEW INPUT DATAS INPUT 1 FOR MAIN SCREEN  INPUT 2'
      PRINT*, ' INPUT = '
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 9997
      ELSE
      GOTO 9999
      ENDIF
 9999 RETURN
      END   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
      SUBROUTINE TE (G,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DX(500),ko
      COMPLEX*16 G,W,E(500),H(500),ALT(2,2),BLT(2,2),DLT(2,2),BX,CI
      COMPLEX*16 LT(2,2,500),ELT(2,2),FLT(2,2),GLT(2,2),ZT(2,2),DET
      COMPLEX*16 LLT(2,2,500),RRT(2,2,500),WRT(2,2),WLT(2,2),RT(2,2,500)
      INTEGER I,NL,ic
      COMMON DX,E,NL,IC,CI       
      COMMON /STATE/ BX,ko
      DO I=1,NL,1
         H(I)=CDSQRT(E(I)-G*G)
         IF (I.EQ.1) THEN
            LT(1,1,I)=CDEXP(CI*H(I)*DX(I))
            LT(1,2,I)=DCMPLX(0.D0,0.D0)
            LT(2,1,I)=CI*H(I)*LT(1,1,I)
            LT(2,2,I)=DCMPLX(0.D0,0.D0)
            ELSEIF (I.EQ.NL) THEN
                   RT(1,1,I)=DCMPLX(0.D0,0.D0)
                   RT(1,2,I)=CDEXP((-CI)*H(I)*DX(I-1))
                   RT(2,1,I)=DCMPLX(0.D0,0.D0)
                   RT(2,2,I)=(-CI)*H(I)*RT(1,2,I) 
            ELSE
         ENDIF
      ENDDO
C
      DO I=2,IC
	 LLT(1,1,I)=CDEXP(CI*H(I)*DX(I-1))
	 LLT(1,2,I)=CDEXP(-CI*H(I)*DX(I-1))
	 LLT(2,1,I)=CI*H(I)*LLT(1,1,I)
	 LLT(2,2,I)=-CI*H(I)*LLT(1,2,I)
	 LT(1,1,I)=CDEXP(CI*H(I)*DX(I))
	 LT(1,2,I)=CDEXP(-CI*H(I)*DX(I))
	 LT(2,1,I)=CI*H(I)*LT(1,1,I)
	 LT(2,2,I)=-CI*H(I)*LT(1,2,I) 
      ENDDO
C
      K=1
      DO I=1,2
	 DO J=1,2
	    WLT(I,J)=LT(I,J,K)
         ENDDO
      ENDDO
C
      DO K=2,IC
	 DO I=1,2
	    DO J=1,2
	       ALT(I,J)=LLT(I,J,K)
	       DLT(I,J)=LT(I,J,K)
            ENDDO
         ENDDO
	 CALL ZLINEAR2(ALT,2,2,DET,INFO,1)
	 CALL MULTI (2,ALT,WLT,BLT)
	 IF (K.LT.IC) THEN
	    CALL MULTI(2,DLT,BLT,WLT)
         ELSEIF (K.EQ.IC) THEN
	    GOTO 12
         ELSE
	 ENDIF
      ENDDO
   12 CONTINUE
C
      DO I=NL-1,IC,-1
	 RRT(1,1,I)=CDEXP(CI*H(I)*DX(I))
	 RRT(1,2,I)=CDEXP(-CI*H(I)*DX(I))
	 RRT(2,1,I)=CI*H(I)*RRT(1,1,I)
	 RRT(2,2,I)=-CI*H(I)*RRT(1,2,I)
	 RT(1,1,I)=CDEXP(CI*H(I)*DX(I-1))
	 RT(1,2,I)=CDEXP(-CI*H(I)*DX(I-1))
	 RT(2,1,I)=CI*H(I)*RT(1,1,I)
	 RT(2,2,I)=-CI*H(I)*RT(1,2,I)
      ENDDO
C
      K=NL
      DO I=1,2
	 DO J=1,2
	    WRT(I,J)=RT(I,J,K)
         ENDDO
      ENDDO
C
      DO K=NL-1,IC,-1
	 DO I=1,2
	    DO J=1,2
	       ELT(I,J)=RRT(I,J,K)
	       FLT(I,J)=RT(I,J,K)
            ENDDO
         ENDDO
	 CALL ZLINEAR2(ELT,2,2,DET,INFO,1)
	 CALL MULTI (2,ELT,WRT,GLT)
	 IF (K.GT.IC) THEN
	    CALL MULTI (2,FLT,GLT,WRT)
         ELSEIF (K.EQ.IC) THEN
	    GOTO 14
         ELSE
	 ENDIF
      ENDDO
   14 CONTINUE
C
      DO I=1,2
	 DO J=1,2
          ZT(I,J)=BLT(I,J)-GLT(I,J)
         ENDDO
      ENDDO
      CALL ZLINEAR2 (ZT,2,2,DET,INFO,2)
      BX=-ZT(1,1)/ZT(1,2)
      W=DET
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE TM (G,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DX(500),ko
      COMPLEX*16 G,W,E(500),H(500),ALT(2,2),BLT(2,2),DLT(2,2),BX,CI
      COMPLEX*16 LT(2,2,500),ELT(2,2),FLT(2,2),GLT(2,2),ZT(2,2),DET
      COMPLEX*16 LLT(2,2,500),RRT(2,2,500),WRT(2,2),WLT(2,2),RT(2,2,500)
      INTEGER I,NL,ic
      COMMON DX,E,NL,IC,CI       
      COMMON /STATE/ BX,ko
      DO I=1,NL,1
         H(I)=CDSQRT(E(I)-G*G)
         IF (I.EQ.1) THEN
            LT(1,1,I)=CDEXP(CI*H(I)*DX(I))
            LT(1,2,I)=DCMPLX(0.D0,0.D0)
            LT(2,1,I)=CI*H(I)*LT(1,1,I)/E(I)
            LT(2,2,I)=DCMPLX(0.D0,0.D0)
            ELSEIF (I.EQ.NL) THEN
                   RT(1,1,I)=DCMPLX(0.D0,0.D0)
                   RT(1,2,I)=CDEXP((-CI)*H(I)*DX(I-1))
                   RT(2,1,I)=DCMPLX(0.D0,0.D0)
                   RT(2,2,I)=(-CI)*H(I)*RT(1,2,I)/E(I) 
            ELSE
         ENDIF
      ENDDO
C
      DO I=2,IC
	 LLT(1,1,I)=CDEXP(CI*H(I)*DX(I-1))
	 LLT(1,2,I)=CDEXP(-CI*H(I)*DX(I-1))
	 LLT(2,1,I)=CI*H(I)*LLT(1,1,I)/E(I)
	 LLT(2,2,I)=-CI*H(I)*LLT(1,2,I)/E(I)
	 LT(1,1,I)=CDEXP(CI*H(I)*DX(I))
	 LT(1,2,I)=CDEXP(-CI*H(I)*DX(I))
	 LT(2,1,I)=CI*H(I)*LT(1,1,I)/E(I)
	 LT(2,2,I)=-CI*H(I)*LT(1,2,I)/E(I) 
      ENDDO
C
      K=1
      DO I=1,2
	 DO J=1,2
	    WLT(I,J)=LT(I,J,K)
         ENDDO
      ENDDO
C
      DO K=2,IC
	 DO I=1,2
	    DO J=1,2
	       ALT(I,J)=LLT(I,J,K)
	       DLT(I,J)=LT(I,J,K)
            ENDDO
         ENDDO
	 CALL ZLINEAR2(ALT,2,2,DET,INFO,1)
	 CALL MULTI (2,ALT,WLT,BLT)
	 IF (K.LT.IC) THEN
	    CALL MULTI(2,DLT,BLT,WLT)
         ELSEIF (K.EQ.IC) THEN
	    GOTO 12
         ELSE
	 ENDIF
      ENDDO
   12 CONTINUE
C
      DO I=NL-1,IC,-1
	 RRT(1,1,I)=CDEXP(CI*H(I)*DX(I))
	 RRT(1,2,I)=CDEXP(-CI*H(I)*DX(I))
	 RRT(2,1,I)=CI*H(I)*RRT(1,1,I)/E(I)
	 RRT(2,2,I)=-CI*H(I)*RRT(1,2,I)/E(I)
	 RT(1,1,I)=CDEXP(CI*H(I)*DX(I-1))
	 RT(1,2,I)=CDEXP(-CI*H(I)*DX(I-1))
	 RT(2,1,I)=CI*H(I)*RT(1,1,I)/E(I)
	 RT(2,2,I)=-CI*H(I)*RT(1,2,I)/E(I)
      ENDDO
C
      K=NL
      DO I=1,2
	 DO J=1,2
	    WRT(I,J)=RT(I,J,K)
         ENDDO
      ENDDO
C
      DO K=NL-1,IC,-1
	 DO I=1,2
	    DO J=1,2
	       ELT(I,J)=RRT(I,J,K)
	       FLT(I,J)=RT(I,J,K)
            ENDDO
         ENDDO
	 CALL ZLINEAR2(ELT,2,2,DET,INFO,1)
	 CALL MULTI (2,ELT,WRT,GLT)
	 IF (K.GT.IC) THEN
	    CALL MULTI (2,FLT,GLT,WRT)
         ELSEIF (K.EQ.IC) THEN
	    GOTO 14
         ELSE
	 ENDIF
      ENDDO
   14 CONTINUE
C
      DO I=1,2
	 DO J=1,2
	    ZT(I,J)=BLT(I,J)-GLT(I,J)
         ENDDO
      ENDDO
      CALL ZLINEAR2 (ZT,2,2,DET,INFO,2)
      BX=-ZT(1,1)/ZT(1,2)
      W=DET
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
      SUBROUTINE TEFUN (ZGA,BXA)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DX(500),DY1,DZ1,DY2,DZ2,DY,DYY,DZZ,ko,XC(500),CONF(500)
      REAL*8 N(500)
      COMPLEX*16 ZGA,E(500),H(500),ALT(2,2),FLT(2,2),A(500),B(500)
      COMPLEX*16 LT(2,2,500),GLT(2,2),ELT(2,2),RT(2,2,500),Bx,PHI
      COMPLEX*16 LLT(2,2,500),RRT(2,2,500),BLT(2,2),DLT(2,2),CI,H1,HNL
      COMPLEX*16 AA,BB,HH,BXA,AAX(500),BBX(500),HHX(500)
      COMPLEX*16 WRT(2,2),WLT(2,2)
      INTEGER I,NL,ic
      COMMON DX,E,NL,IC,CI       
      COMMON /STATE/ BX,ko
      COMMON /PSIF/ HH,AA,BB,CONF                                 
      COMMON /FAR/ H1,HNL 
      COMMON /FAR1/ H,A,B,XNORN
      COMMON /A1/ N
      COMMON /A5/ AAX,BBX,HHX                 
      EXTERNAL GSP   
c     open(35,file='coef.dat',status='unknown')
      PI=3.1415926535897932D0
      A(1)=DCMPLX(1.D0,0.D0)
      B(NL)=BXA
      DO I=1,NL,1
         H(I)=CDSQRT(E(I)-ZGA*ZGA)
         IF (I.EQ.1) THEN
            LT(1,1,I)=CDEXP(CI*H(I)*DX(I))
            LT(1,2,I)=DCMPLX(0.D0,0.D0)
            LT(2,1,I)=CI*H(I)*LT(1,1,I)
            LT(2,2,I)=DCMPLX(0.D0,0.D0)
            ELSEIF (I.EQ.NL) THEN
                   RT(1,1,I)=DCMPLX(0.D0,0.D0)
                   RT(1,2,I)=CDEXP((-CI)*H(I)*DX(I-1))
                   RT(2,1,I)=DCMPLX(0.D0,0.D0)
                   RT(2,2,I)=(-CI)*H(I)*RT(1,2,I) 
            ELSE
         ENDIF
      ENDDO
C
      DO I=2,IC
	 LLT(1,1,I)=CDEXP(CI*H(I)*DX(I-1))
	 LLT(1,2,I)=CDEXP(-CI*H(I)*DX(I-1))
	 LLT(2,1,I)=CI*H(I)*LLT(1,1,I)
	 LLT(2,2,I)=-CI*H(I)*LLT(1,2,I)
	 LT(1,1,I)=CDEXP(CI*H(I)*DX(I))
	 LT(1,2,I)=CDEXP(-CI*H(I)*DX(I))
	 LT(2,1,I)=CI*H(I)*LT(1,1,I)
	 LT(2,2,I)=-CI*H(I)*LT(1,2,I) 
      ENDDO
C
      K=1
      DO I=1,2
	 DO J=1,2
	    WLT(I,J)=LT(I,J,K)
         ENDDO
      ENDDO
C
      DO K=2,IC
	 DO I=1,2
	    DO J=1,2
	       ALT(I,J)=LLT(I,J,K)
	       DLT(I,J)=LT(I,J,K)
            ENDDO
         ENDDO
	 CALL ZLINEAR2(ALT,2,2,DET,INFO,1)
	 CALL MULTI (2,ALT,WLT,BLT)
	 A(K)=BLT(1,1)*A(1)+BLT(1,2)*B(NL)
	 B(K)=BLT(2,1)*A(1)+BLT(2,2)*B(NL)
         IF (K.LT.IC) THEN
	    CALL MULTI(2,DLT,BLT,WLT)
         ELSEIF (K.EQ.IC) THEN
            GOTO 12
         ELSE
         ENDIF
      ENDDO
   12 CONTINUE
c      WRITE(32,*) A(K),B(K),K
C
      DO I=NL-1,IC,-1
	 RRT(1,1,I)=CDEXP(CI*H(I)*DX(I))
	 RRT(1,2,I)=CDEXP(-CI*H(I)*DX(I))
	 RRT(2,1,I)=CI*H(I)*RRT(1,1,I)
	 RRT(2,2,I)=-CI*H(I)*RRT(1,2,I)
	 RT(1,1,I)=CDEXP(CI*H(I)*DX(I-1))
	 RT(1,2,I)=CDEXP(-CI*H(I)*DX(I-1))
	 RT(2,1,I)=CI*H(I)*RT(1,1,I)
	 RT(2,2,I)=-CI*H(I)*RT(1,2,I)
      ENDDO
C
      K=NL
      DO I=1,2
	 DO J=1,2
	    WRT(I,J)=RT(I,J,K)
         ENDDO
      ENDDO
C
      DO K=NL-1,IC,-1
	 DO I=1,2
	    DO J=1,2
	       ELT(I,J)=RRT(I,J,K)
	       FLT(I,J)=RT(I,J,K)
            ENDDO
         ENDDO
	 CALL ZLINEAR2(ELT,2,2,DET,INFO,1)
	 CALL MULTI (2,ELT,WRT,GLT)
	 A(K)=GLT(1,1)*A(1)+GLT(1,2)*B(NL)
	 B(K)=GLT(2,1)*A(1)+GLT(2,2)*B(NL)
 	 IF (K.GT.IC) THEN
	    CALL MULTI (2,FLT,GLT,WRT)
         ELSEIF (K.EQ.IC) THEN
            GOTO 14
         ELSE
         ENDIF
      ENDDO
   14 CONTINUE
c      WRITE(32,*) A(K),B(K),K
C
      XNORM=0.D0                               
c     XC(1)=((CDABS(A(1)))**2/(2*CDABS(H(1))))
c    +      *DEXP(-2.*(CDABS(H(1)*DX(1))))
      XC(1)=(A(1)*CDEXP(CI*H(1)*DX(1))
     +      *DCONJG(A(1)*CDEXP(CI*H(1)*DX(1))))
     +       /(2.*CDABS(H(1)))
c     XC(NL)=((CDABS(B(NL)))**2/(2*CDABS(H(NL))))
c    +       *DEXP(-2.*(CDABS(H(NL)*DX(NL-1))))
      XC(NL)=(B(NL)*CDEXP(-CI*H(NL)*DX(NL-1))
     +       *DCONJG(B(NL)*CDEXP(-CI*H(NL)*DX(NL-1))))
     +       /(2.*CDABS(H(NL)))
      H1=H(1)
      HNL=H(NL)
      DO I=2,NL-1
         AA=A(I)
         BB=B(I)     
         HH=H(I)
         xl=dx(i)
         xR=dx(i-1)
         CALL QUANC8(gsp,xl,xr,0.d0,1.d-15,cc,errest,nofun,Xflag)
         XC(I)=CC 
         XNORM=XNORM+CC
      ENDDO
      XNORM=(XNORM+XC(1)+XC(NL))
      XNORN=XNORM                
      DO I=1,NL
         AAX(I)=A(I)/DSQRT(XNORM)
         BBX(I)=B(I)/DSQRT(XNORM)
         HHX(I)=H(I)  
      ENDDO
      DO I=1,NL
      CONF(I)=(XC(I))/XNORN
      ENDDO  
      XNOR=DSQRT(XNORN)
c      WRITE(31,*) XNOR,XC(1),XC(NL)
      I=1
      DY1=DX(I)+0.00001
      DYY1=DX(I)
      DO DZ1=DY1,DYY1,(DYY1-DY1)/10.D0
         PHI=(AAX(1))*CDEXP(CI*H(1)*DZ1)
         XINT=DSQRT(DREAL(PHI)**2+DIMAG(PHI)**2)
         XANG=(DATAN2(DIMAG(PHI),DREAL(PHI)))*180.D0/PI
         WRITE(3,103) DZ1,DREAL(PHI),DIMAG(PHI),XANG,XINT,N(I)
      ENDDO     
C         
      DO I=2,NL-1
         DY=DX(I-1)
         DYY=DX(I)
         DO DZZ=DY,DYY,(DYY-DY)/10.D0
c        PHI=(AAX(I))*CDCOS(H(I)*DZZ)+(BBX(I))*CDSIN(H(I)*DZZ)
	 PHI=(AAX(I)*CDEXP(CI*H(I)*DZZ)+BBX(I)*CDEXP(-CI*H(I)*DZZ))
            XINT=DSQRT(DREAL(PHI)**2+DIMAG(PHI)**2)
            XANG=(DATAN2(DIMAG(PHI),DREAL(PHI)))*180.D0/PI
            WRITE(3,103) DZZ,DREAL(PHI),DIMAG(PHI),XANG,XINT,N(I)
c	    WRITE(31,104) CDCOS(H(I)*DZZ),CDSIN(H(I)*DZZ),DZZ,I
         ENDDO
      ENDDO
C      
      I=NL 
      DY2=DX(I-1)-0.00001
      DYY2=DX(I-1)
      DO DZ2=DYY2,DY2,(DY2-DYY2)/10.D0
         PHI=(BBX(NL))*CDEXP((-CI)*H(NL)*DZ2) 
         XINT=DSQRT(DREAL(PHI)**2+DIMAG(PHI)**2)
         XANG=(DATAN2(DIMAG(PHI),DREAL(PHI)))*180.D0/PI
         WRITE(3,103) DZ2,DREAL(PHI),DIMAG(PHI),XANG,XINT,N(NL)
      ENDDO
  103 FORMAT(1X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8,1X,F8.5)
      RETURN
      END
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
      SUBROUTINE TMFUN (ZGA,BXA)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DX(500),DY1,DZ1,DY2,DZ2,DY,DYY,DZZ,ko,XC(500),CONF(500)
      REAL*8 N(500)
      COMPLEX*16 ZGA,E(500),H(500),ALT(2,2),FLT(2,2),A(500),B(500)
      COMPLEX*16 LT(2,2,500),GLT(2,2),ELT(2,2),RT(2,2,500),Bx,PHI
      COMPLEX*16 LLT(2,2,500),RRT(2,2,500),BLT(2,2),DLT(2,2),CI,H1,HNL
      COMPLEX*16 AA,BB,HH,BXA,AAX(500),BBX(500),HHX(500)
      COMPLEX*16 WRT(2,2),WLT(2,2)
      INTEGER I,NL,ic
      COMMON DX,E,NL,IC,CI       
      COMMON /STATE/ BX,ko
      COMMON /PSIF/ HH,AA,BB,CONF                                 
      COMMON /FAR/ H1,HNL 
      COMMON /FAR1/ H,A,B,XNORN
      COMMON /A1/ N
      COMMON /A5/ AAX,BBX,HHX                 
      EXTERNAL GSP   
c     open(35,file='coef.dat',status='unknown')
      PI=3.1415926535897932D0
      A(1)=DCMPLX(1.D0,0.D0)
      B(NL)=BXA
      DO I=1,NL,1
         H(I)=CDSQRT(E(I)-ZGA*ZGA)
         IF (I.EQ.1) THEN
            LT(1,1,I)=CDEXP(CI*H(I)*DX(I))
            LT(1,2,I)=DCMPLX(0.D0,0.D0)
            LT(2,1,I)=CI*H(I)*LT(1,1,I)/E(I)
            LT(2,2,I)=DCMPLX(0.D0,0.D0)
            ELSEIF (I.EQ.NL) THEN
                   RT(1,1,I)=DCMPLX(0.D0,0.D0)
                   RT(1,2,I)=CDEXP((-CI)*H(I)*DX(I-1))
                   RT(2,1,I)=DCMPLX(0.D0,0.D0)
                   RT(2,2,I)=(-CI)*H(I)*RT(1,2,I)/E(I) 
            ELSE
         ENDIF
      ENDDO
C
      DO I=2,IC
	 LLT(1,1,I)=CDEXP(CI*H(I)*DX(I-1))
	 LLT(1,2,I)=CDEXP(-CI*H(I)*DX(I-1))
	 LLT(2,1,I)=CI*H(I)*LLT(1,1,I)/E(I)
	 LLT(2,2,I)=-CI*H(I)*LLT(1,2,I)/E(I)
	 LT(1,1,I)=CDEXP(CI*H(I)*DX(I))
	 LT(1,2,I)=CDEXP(-CI*H(I)*DX(I))
	 LT(2,1,I)=CI*H(I)*LT(1,1,I)/E(I)
	 LT(2,2,I)=-CI*H(I)*LT(1,2,I)/E(I) 
      ENDDO
C
      K=1
      DO I=1,2
	 DO J=1,2
	    WLT(I,J)=LT(I,J,K)
         ENDDO
      ENDDO
C
      DO K=2,IC
	 DO I=1,2
	    DO J=1,2
	       ALT(I,J)=LLT(I,J,K)
	       DLT(I,J)=LT(I,J,K)
            ENDDO
         ENDDO
	 CALL ZLINEAR2(ALT,2,2,DET,INFO,1)
	 CALL MULTI (2,ALT,WLT,BLT)
	 A(K)=BLT(1,1)*A(1)+BLT(1,2)*B(NL)
	 B(K)=BLT(2,1)*A(1)+BLT(2,2)*B(NL)
         IF (K.LT.IC) THEN
	    CALL MULTI(2,DLT,BLT,WLT)
         ELSEIF (K.EQ.IC) THEN
            GOTO 12
         ELSE
         ENDIF
      ENDDO
   12 CONTINUE
c      WRITE(32,*) A(K),B(K),K
C
      DO I=NL-1,IC,-1
	 RRT(1,1,I)=CDEXP(CI*H(I)*DX(I))
	 RRT(1,2,I)=CDEXP(-CI*H(I)*DX(I))
	 RRT(2,1,I)=CI*H(I)*RRT(1,1,I)/E(I)
	 RRT(2,2,I)=-CI*H(I)*RRT(1,2,I)/E(I)
	 RT(1,1,I)=CDEXP(CI*H(I)*DX(I-1))
	 RT(1,2,I)=CDEXP(-CI*H(I)*DX(I-1))
	 RT(2,1,I)=CI*H(I)*RT(1,1,I)/E(I)
	 RT(2,2,I)=-CI*H(I)*RT(1,2,I)/E(I)
      ENDDO
C
      K=NL
      DO I=1,2
	 DO J=1,2
	    WRT(I,J)=RT(I,J,K)
         ENDDO
      ENDDO
C
      DO K=NL-1,IC,-1
	 DO I=1,2
	    DO J=1,2
	       ELT(I,J)=RRT(I,J,K)
	       FLT(I,J)=RT(I,J,K)
            ENDDO
         ENDDO
	 CALL ZLINEAR2(ELT,2,2,DET,INFO,1)
	 CALL MULTI (2,ELT,WRT,GLT)
	 A(K)=GLT(1,1)*A(1)+GLT(1,2)*B(NL)
	 B(K)=GLT(2,1)*A(1)+GLT(2,2)*B(NL)
 	 IF (K.GT.IC) THEN
	    CALL MULTI (2,FLT,GLT,WRT)
         ELSEIF (K.EQ.IC) THEN
            GOTO 14
         ELSE
         ENDIF
      ENDDO
   14 CONTINUE
c      WRITE(32,*) A(K),B(K),K
C
      XNORM=0.D0                               
c     XC(1)=((CDABS(A(1)))**2/(2*CDABS(H(1))))
c    +      *DEXP(-2.*(CDABS(H(1)*DX(1))))
      XC(1)=(A(1)*CDEXP(CI*H(1)*DX(1))
     +      *DCONJG(A(1)*CDEXP(CI*H(1)*DX(1))))
     +       /(2.*CDABS(H(1)))
c     XC(NL)=((CDABS(B(NL)))**2/(2*CDABS(H(NL))))
c    +       *DEXP(-2.*(CDABS(H(NL)*DX(NL-1))))
      XC(NL)=(B(NL)*CDEXP(-CI*H(NL)*DX(NL-1))
     +       *DCONJG(B(NL)*CDEXP(-CI*H(NL)*DX(NL-1))))
     +       /(2.*CDABS(H(NL)))
      H1=H(1)
      HNL=H(NL)
      DO I=2,NL-1
         AA=A(I)
         BB=B(I)     
         HH=H(I)
         xl=dx(i)
         xR=dx(i-1)
         CALL QUANC8(gsp,xl,xr,0.d0,1.d-15,cc,errest,nofun,Xflag)
         XC(I)=CC 
         XNORM=XNORM+CC
      ENDDO
      XNORM=(XNORM+XC(1)+XC(NL))
      XNORN=XNORM                
      DO I=1,NL
         AAX(I)=A(I)/DSQRT(XNORM)
         BBX(I)=B(I)/DSQRT(XNORM)
         HHX(I)=H(I)  
      ENDDO
      DO I=1,NL
      CONF(I)=(XC(I))/XNORN
      ENDDO  
      XNOR=DSQRT(XNORN)
c      WRITE(31,*) XNOR,XC(1),XC(NL)
      I=1
      DY1=DX(I)+0.0050
      DYY1=DX(I)
      DO DZ1=DY1,DYY1,(DYY1-DY1)/10.D0
         PHI=(AAX(I))*CDEXP(CI*H(I)*DZ1)
         XINT=DSQRT(DREAL(PHI)**2+DIMAG(PHI)**2)
         XANG=(DATAN2(DIMAG(PHI),DREAL(PHI)))*180.D0/PI
         WRITE(3,103) DZ1,DREAL(PHI),DIMAG(PHI),XANG,XINT,N(I)
      ENDDO     
C         
      DO I=2,NL-1
         DY=DX(I-1)
         DYY=DX(I)
         DO DZZ=DY,DYY,(DYY-DY)/10.D0
c        PHI=(AAX(I))*CDCOS(H(I)*DZZ)+(BBX(I))*CDSIN(H(I)*DZZ)
	 PHI=(AAX(I)*CDEXP(CI*H(I)*DZZ)+BBX(I)*CDEXP(-CI*H(I)*DZZ))
            XINT=DSQRT(DREAL(PHI)**2+DIMAG(PHI)**2)
            XANG=(DATAN2(DIMAG(PHI),DREAL(PHI)))*180.D0/PI
            WRITE(3,103) DZZ,DREAL(PHI),DIMAG(PHI),XANG,XINT,N(I)
c	    WRITE(31,104) CDCOS(H(I)*DZZ),CDSIN(H(I)*DZZ),DZZ,I
         ENDDO
      ENDDO
C      
      I=NL 
      DY2=DX(I-1)-0.0050
      DYY2=DX(I-1)
      DO DZ2=DYY2,DY2,(DY2-DYY2)/10.D0
         PHI=(BBX(I))*CDEXP((-CI)*H(I)*DZ2) 
         XINT=DSQRT(DREAL(PHI)**2+DIMAG(PHI)**2)
         XANG=(DATAN2(DIMAG(PHI),DREAL(PHI)))*180.D0/PI
         WRITE(3,103) DZ2,DREAL(PHI),DIMAG(PHI),XANG,XINT,N(I)
      ENDDO
  103 FORMAT(1X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8,1X,F8.5)
      RETURN
      END
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         
      SUBROUTINE SEARCH (XDG,XBX,XI)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER XI(500) 
      LOGICAL FLAG
      COMPLEX*16 XDG(500),XBX(500)
      DO I=19,1,-1
         FLAG=.FALSE.
         DO J=1,I
            IF (DREAL(XDG(J)).GE.DREAL(XDG(J+1))) THEN
               T=XDG(J)         
               TT=XBX(J)
               TTT=XI(J)
               XDG(J)=XDG(J+1)
               XBX(J)=XBX(J+1)
               XI(J)=XI(J+1)
               XDG(J+1)=T
               XBX(J+1)=TT
               XI(J+1)=TTT
               FLAG=.TRUE.
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                     C                                                
      SUBROUTINE XSEARCH (XFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 FX(181)
      LOGICAL FLAG
      COMMON /XMAX/ FX
      DO I=180,1,-1
         FLAG=.FALSE.
         DO J=1,I
            IF (FX(J).GE.FX(J+1)) THEN
               T=FX(J)         
               FX(J)=FX(J+1)
               FX(J+1)=T
               FLAG=.TRUE.
            ENDIF
         ENDDO
      ENDDO
      XFF=FX(181)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
      SUBROUTINE FFIELD (THETA,BXA,FF)                                 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NL,IC
      REAL*8 KO,DX(500),N(500),LAMBDA,XNORN,BR,BI
      COMPLEX*16 A(500),B(500),H(500),PA,XSUM,XHX,XLX
      COMPLEX*16 CI,E(500),BX,BXA,PB,AAX,BBX,HHX,FA,BETA
      COMPLEX*16 AA1(500),A1(500),A2(500)
      COMMON DX,E,NL,IC,CI
      COMMON /STATE/ BX,ko
      COMMON /FAR1/ H,A,B,XNORN       
      COMMON /FF1/ AAX,BBX,HHX,PA
      COMMON /A1/ N
      COMMON /A10/ LAMBDA,BR,BI
      EXTERNAL AAGSP
      XNOR=DSQRT(XNORN) 
      BETA=DCMPLX(BR,BI)
      PI=3.1415926535897932D0
      PA=CI*ko*DSIN(THETA)
      PB=PA*PA
         XSUM=DCMPLX(0.D0,0.D0)  
         XHX=(-1*(A(1))*CDEXP((CI*H(1)+PA)*DX(1))/(CI*H(1)+PA))
         AA1(1)=(2/(1+N(1)))*(N(1)+(BETA/N(1)))
      A1(1)=2*DCOS(THETA)/(DCOS(THETA)+DSQRT(N(1)**2-DSIN(THETA)**2))
	 A2(1)=A1(1)*(DSQRT(N(1)**2-DSIN(THETA)**2)+(BETA/N(1)))
	 XSUM=XSUM+(A2(1)/AA1(1))*XHX
         DO J=2,NL-1
c           FA=DCMPLX(0.D0,0.D0)
            XL=DX(J)
            XH=DX(J-1)
            AAX=A(J)
            BBX=B(J)
            HHX=H(J)
         CALL CQUANC8 (AAgsp,XL,XH,0.d0,1.d-15,FA,errest,nofun,flag)
	 AA1(J)=(2/(1+N(J)))*(N(J)+(BETA/N(J)))
      A1(J)=2*DCOS(THETA)/(DCOS(THETA)+DSQRT(N(J)**2-DSIN(THETA)**2))
      A2(J)=A1(J)*(DSQRT(N(J)**2-DSIN(THETA)**2)+(BETA/N(J)))
         XSUM=XSUM+(A2(J)/AA1(J))*FA
         ENDDO
         XLX=(BXA)*CDEXP((-CI*H(NL)+PA)*DX(NL-1))/(-CI*H(NL)+PA) 
	 AA1(NL)=(2/(1+N(NL)))*(N(NL)+(BETA/N(NL)))
	Z1=N(NL)**2-DSIN(THETA)**2
	 IF (Z1.GT.0.D0) THEN
	    A3=CDSQRT(DCMPLX(Z1,0.D0))
	 ELSE
	    A3=CDSQRT(DCMPLX(0.D0,Z1))
	 ENDIF
      A1(NL)=2*DCOS(THETA)/(DCOS(THETA)+A3) 
	 A2(NL)=A1(NL)*(A3+(BETA/N(NL)))
	 XSUM=XSUM+(A2(NL)/AA1(NL))*XLX
	 FF=(XSUM*DCONJG(XSUM))
      RETURN
      END
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE COUPL1 (NG1,NG2,NC,CCOEF1,CCOEF2,CCOEF3,CCOEF4)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 CCOEF1,E(500),A1(500),B1(500),B2(500),H2(500)
      COMPLEX*16 H1(500),AA1,BB1,HH1,AA2,BB2,HH2,A2(500),CC1,ZD1,ZD2
      COMPLEX*16 C1(500),C2(500),CCOEF2,CC2,BX,CI,NA,NB,CCOEF3,CCOEF4
      COMPLEX*16 CH1,CA1,CBL,CHL,XCCOEF1,XCCOEF2,CC11,CC22
      COMPLEX*16 C11(500),C22(500),XCCOEF3,XCCOEF4  
      REAL*8 DX(500),KO,DIST,LAMBDA,NG1,NG2,NC
      INTEGER NL,IC
      COMMON DX,E,NL,IC,CI
      COMMON /STATE/ BX,KO
      COMMON /AC1/ A1,B1,H1
      COMMON /AC2/ A2,B2,H2
      COMMON /AC3/ AA1,BB1,HH1
      COMMON /AC3A/ AA2,BB2,HH2 
      COMMON /AC4/ DIST
c     COMMON /A21/ NA,NB  
      EXTERNAL GGSP,HGSP,GGSPA,HGSPA    
      W=KO*3.D8                                                   
      U=4*3.1415926535897932D0*1.D-7
      EO=8.854E-12
      CH1=DCONJG(H1(1))
      CA1=DCONJG(A1(1))
      CBL=DCONJG(B1(NL))
      CHL=DCONJG(H1(NL)) 
c     ZD1=-1.*CA1*A2(1)*CDEXP(CI*CH1*DIST)*CDEXP(CI*(H2(1)-CH1)*DX(1))
c    +    /(CI*(H2(1)-CH1))
c     ZD2=CBL*B2(NL)*CDEXP(-1*CI*(CHL*DIST))*CDEXP(CI*(CHL-H2(NL))
c    +    *DX(NL-1))/(CI*(CHL-H2(NL)))
      XCCOEF1=DCMPLX(0.D0,0.D0)
      XCCOEF2=DCMPLX(0.D0,0.D0)
      XCCOEF3=DCMPLX(0.D0,0.D0)
      XCCOEF4=DCMPLX(0.D0,0.D0)
      DO I=2,NL-3
         DL=DX(I)
         DH=DX(I-1)      
         AA1=A1(I)     
         BB1=B1(I)
         HH1=H1(I) 
         AA2=A2(I)   
         BB2=B2(I)
         HH2=H2(I)  
         CALL cquanc8(ggsp,DL,DH,0.d0,1.d-12,cc1,errest,nofun,flag)
         CALL cquanc8(ggspa,DL,DH,0.d0,1.d-12,cc11,errest,nofun,flag) 
         C1(I)=CC1
         C11(I)=CC11
         XCCOEF1=XCCOEF1+C1(I)
         XCCOEF3=XCCOEF3+C11(I)
      ENDDO      
      DO I=4,NL-1
         DL=DX(I)
         DH=DX(I-1)      
         AA1=A1(I)     
         BB1=B1(I)
         HH1=H1(I) 
         AA2=A2(I)   
         BB2=B2(I)
         HH2=H2(I)   
         CALL cquanc8(hgsp,DL,DH,0.d0,1.d-12,cc2,errest,nofun,flag)
         CALL cquanc8(hgspa,DL,DH,0.d0,1.d-12,cc22,errest,nofun,flag)
         C2(I)=CC2
         C22(I)=CC22
         XCCOEF2=XCCOEF2+C2(I)
         XCCOEF4=XCCOEF4+C22(I) 
      ENDDO      
       CCOEF1=(W*EO/4)*(NC**2-NG1**2)*(XCCOEF1)
       CCOEF2=(W*EO/4)*(NC**2-NG2**2)*(XCCOEF2)
       CCOEF3=(W*EO/4)*(NG1**2-NC**2)*(XCCOEF3)
       CCOEF4=(W*EO/4)*(NG2**2-NC**2)*(XCCOEF4)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C
      SUBROUTINE COUPLE2 (NG1,NC,CCOEF)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 CCOEF,E(500),A1(500),B1(500),B2(500),H2(500),CCOEF1
      COMPLEX*16 H1(500),AA1,BB1,HH1,AA2,BB2,HH2,A2(500)
      COMPLEX*16 C1(500),CC1,CC2,C2(500),CCOEF2,CI,BX
      REAL*8 KO,NC,NG1,GA,DX(500),DUP
      INTEGER NL,IC,IX
      COMMON DX,E,NL,IC,CI
      COMMON /STATE/ BX,KO
      COMMON /AC1/ A1,B1,H1
      COMMON /AC2/ A2,B2,H2
      COMMON /AC3/ AA1,BB1,HH1        
      COMMON /AC3A/ AA2,BB2,HH2
      COMMON /GRATINGA/ GA                
      COMMON /A8/ DUP  
      COMMON /A9/ IX
      EXTERNAL XXGSP   
      PRINT*, ' INPUT RATIO FOR GRATING HEIGTH'
      print*, ' INPUT R1 FOR UPPER LEVEL, R2 FOR LOWER LEVEL'
      READ(*,*) RATIO1,RATIO2
      W=KO*(3.D14)  ! (C=3.D14 um)
      EO=8.854E-18  ! ( um unit)
      DH=DX(IX)+RATIO1*GA
      DL=DX(IX)-RATIO2*GA                    
      DUP=DX(IX)             
      CCOEF1=DCMPLX(0.D0,0.D0)
      CCOEF2=DCMPLX(0.D0,0.D0)                    
      DO I=IX,IX
         HH1=H1(I)     
         AA1=A1(I)     
         BB1=B1(I)   
         HH2=H2(I)  
         AA2=A2(I)   
         BB2=B2(I)   
         CALL cquanc8(XXgsp,Dup,Dh,0.d0,1.d-15,cc1,errest,nofun,flag)
         C1(I)=CC1  
         CCOEF1=CCOEF1+C1(I)                                        
      ENDDO       
      DO I=IX+1,IX+1
         HH1=H1(I)     
         AA1=A1(I)     
         BB1=B1(I)   
         HH2=H2(I)  
         AA2=A2(I)   
         BB2=B2(I)   
         CALL cquanc8(XXgsp,Dl,dup,0.d0,1.d-15,cc2,errest,nofun,flag)
         C2(I)=CC2          
         CCOEF2=CCOEF2+C2(I)                                        
      ENDDO
      CCOEF=-(W*EO/(8.D0*GA))*(NG1**2-NC**2)*(CCOEF1+CCOEF2)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
      SUBROUTINE RECT (DXX,DYY,KAB)     
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 KUP,KBELOW,KAB,E(500),CI,C1
      REAL*8 DX(500),DUP         
      INTEGER NL,IC
      COMMON /A2/ C1               
      COMMON /A8/ DUP,DMID
      COMMON DX,E,NL,IC,CI       
      EXTERNAL XGSP,YGSP    
      CALL cquanc8(xgsp,DYY,DMID,0.d0,1.d-15,kup,errest,nofun,flag)
      CALL cquanc8(ygsp,DMID,DXX,0.d0,1.d-15,kbelow,errest,nofun,flag) 
      KAB=-CI*C1*(KUP-KBELOW)
      RETURN
      END               
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SINUS (DXX,DYY,KAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 KAB,CI,KUP,KBELOW,C1,E(500)
      REAL*8 DXX,DYY,DUP,DX(500)                       
      INTEGER M,NL,IC
      COMMON DX,E,NL,IC,CI
      COMMON /A2/ C1 
      COMMON /A8/ DUP,DMID 
      COMMON /period/ G,GL,M
      EXTERNAL AGSP,BGSP
      CALL cquanc8(agsp,DYY,DMID,0.d0,1.d-15,kup,errest,nofun,flag)
      CALL cquanc8(bgsp,DMID,DXX,0.d0,1.d-15,kbelow,errest,nofun,flag)
      KAB=2.D0*C1*(KUP-((-1.D0)**M)*KBELOW)
      RETURN
      END  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE TRI (DXX,DYY,KAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 KAB,CI,KUP,KBELOW,C1,E(500)
      REAL*8 DUP,G,DX(500)
      INTEGER M,NL,IC   
      COMMON DX,E,NL,IC,CI
      COMMON /A2/ C1                               
      COMMON /period/ G,GL,M 
      COMMON /A8/ DUP,DMID
      EXTERNAL CGSP,DGSP
      CALL cquanc8(cgsp,DYY,DMID,0.d0,1.d-15,kup,errest,nofun,flag)
      CALL cquanc8(dgsp,DMID,DXX,0.d0,1.d-15,kbelow,errest,nofun,flag) 
      KAB=-CI*C1*(KUP-KBELOW)
      RETURN
      END  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SAW (DXX,DYY,KAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 KAB,CI,KUP,KBELOW,C1,E(500)
      REAL*8 DX(500),DUP,G
      INTEGER M,NL,IC    
      COMMON DX,E,NL,IC,CI
      COMMON /A2/ C1                               
      COMMON /period/ G,GL,M 
      COMMON /A8/ DUP,DMID
      EXTERNAL EGSP,FGSP       
      CALL cquanc8(egsp,DYY,DMID,0.d0,1.d-15,kup,errest,nofun,flag)  
      CALL cquanc8(fgsp,DMID,DXX,0.d0,1.d-15,kbelow,errest,nofun,flag)
      KAB=-CI*C1*(KUP-KBELOW)
      RETURN
      END
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION AAGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 AAX,BBX,HHX,PA,CI
      COMMON /FF1/ AAX,BBX,HHX,PA
	CI=DCMPLX(0.D0,1.D0)
C      AAGSP=(AAX*CDCOS(HHX*X)+BBX*CDSIN(HHX*X))*CDEXP(PA*X)
      AAGSP=(AAX*CDEXP(CI*HHX*X)+BBX*CDEXP(-CI*HHX*X))*CDEXP(PA*X)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         
      DOUBLE PRECISION FUNCTION GSP(X) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CONF(500)
      COMPLEX*16 HH,AA,BB,CI
      COMMON /PSIF/ HH,AA,BB,CONF
      CI=DCMPLX(0.D0,1.D0)
c     GSP=(CDABS(AA*CDCOS(HH*X)+BB*CDSIN(HH*X)))**2
c     GSP=(AA*CDCOS(HH*X)+BB*CDSIN(HH*X))*
c    +    DCONJG(AA*CDCOS(HH*X)+BB*CDSIN(HH*X))
      GSP=(AA*CDEXP(CI*HH*X)+BB*CDEXP(-CI*HH*X))*
     +    DCONJG(AA*CDEXP(CI*HH*X)+BB*CDEXP(-CI*HH*X)) 
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION XGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16  AUP,BUP,HUP,CXUP,E(500),CI
      REAL*8 X,DX(500),DUP              
      INTEGER LUP,LBELOW,NL,IC 
      COMMON DX,E,NL,IC,CI
      COMMON /RUP/ AUP,BUP,HUP,CXUP  
      COMMON /A7/ LUP,LBELOW
      COMMON /A8/ DUP,DMID
      IF (LUP.EQ.1) THEN
      XGSP=CXUP*(AUP*CDEXP(CI*HUP*X))*(AUP*CDEXP(CI*HUP*X))
      ELSEIF (LUP.EQ.NL) THEN
      XGSP=CXUP*(BUP*CDEXP(-CI*HUP*X))*(BUP*CDEXP(-CI*HUP*X))
      ELSE
C      XGSP=CXUP*(AUP*CDCOS(HUP*X)+BUP*CDSIN(HUP*X))
C     +     *(AUP*CDCOS(HUP*X)+BUP*CDSIN(HUP*X)) 
      XGSP=CXUP*(AUP*CDEXP(CI*HUP*X)+BUP*CDEXP(-CI*HUP*X))
     +	*(AUP*CDEXP(CI*HUP*X)+BUP*CDEXP(-CI*HUP*X)) 
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION YGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)                          
      REAL*8 DX(500),X,DUP
      INTEGER NL,LUP,LBELOW,IC
      COMPLEX*16 ABELOW,BBELOW,HBELOW,CXBELOW,E(500),CI
      COMMON /RBELOW/ ABELOW,BBELOW,HBELOW,CXBELOW
      COMMON DX,E,NL,IC,CI       
      COMMON /A7/ LUP,LBELOW  
      COMMON /A8/ DUP,DMID
      IF (LBELOW.EQ.NL) THEN 
      YGSP=CXBELOW*(BBELOW*CDEXP(-CI*HBELOW*X))
     +     *(BBELOW*CDEXP(-CI*HBELOW*X)) 
      ELSEIF (LBELOW.EQ.1) THEN
      YGSP=CXBELOW*(ABELOW*CDEXP(CI*HBELOW*X))
     +     *(ABELOW*CDEXP(CI*HBELOW*X))
      ELSE
C      YGSP=CXBELOW*((ABELOW*CDCOS(HBELOW*(X))
C     +     +BBELOW*CDSIN(HBELOW*(X)))*
C     +     (ABELOW*CDCOS(HBELOW*(X))
C     +     +BBELOW*CDSIN(HBELOW*(X))))
      YGSP=CXBELOW*(ABELOW*CDEXP(CI*HBELOW*X)
     +     +BBELOW*CDEXP(-CI*HBELOW*X))*
     +	   (ABELOW*CDEXP(CI*HBELOW*X)+
     +     BBELOW*CDEXP(-CI*HBELOW*X))
      ENDIF 
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION AGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 E(500),AUP,BUP,HUP,CXUP,CI 
      REAL*8 DX(500),G,X,GL,DUP,P
      INTEGER M,LUP,LBELOW,NL,IC
      COMMON DX,E,NL,IC,CI    
      COMMON /period/ G,GL,M
      COMMON /RUP/ AUP,BUP,HUP,CXUP 
      COMMON /A7/ LUP,LBELOW  
      COMMON /A8/ DUP,DMID
      P=2*(X-DMID)/G
      IF (P.LT.0.D0) THEN
      P=0.D0
      ELSEIF (P.GT.1.D0) THEN
      P=1.D0
      ELSE
	P=P
      ENDIF
      IF (LUP.EQ.1) THEN
      AGSP=CDSIN(M*DCMPLX(DACOS(P),0.D0))
     +     *(AUP*CDEXP(CI*HUP*(X)))
     +     *(AUP*CDEXP(CI*HUP*(X))) 
      ELSEIF (LUP.EQ.NL) THEN
      AGSP=CDSIN(M*DCMPLX(DACOS(P),0.D0))
     +     *(BUP*CDEXP(-CI*HUP*(X)))
     +     *(BUP*CDEXP(-CI*HUP*(X))) 
      ELSE
      AGSP=CDSIN(M*DCMPLX(DACOS(P),0.D0))
     +      *(AUP*CDEXP(CI*HUP*X)+BUP*CDEXP(-CI*HUP*X))
     +      *(AUP*CDEXP(CI*HUP*X)+BUP*CDEXP(-CI*HUP*X))
C     +     *(AUP*CDCOS(HUP*(X))+BUP*CDSIN(HUP*(X)))
C     +     *(AUP*CDCOS(HUP*(X))+BUP*CDSIN(HUP*(X)))
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION BGSP(X)    
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 E(500),ABELOW,BBELOW,HBELOW,CXBELOW,CI 
      REAL*8 DX(500),G,X,GL,DUP,P
      INTEGER M,NL,LUP,LBELOW,IC     
      COMMON /period/ G,GL,M
      COMMON /RBELOW/ ABELOW,BBELOW,HBELOW,CXBELOW        
      COMMON DX,E,NL,IC,CI
      COMMON /A7/ LUP,LBELOW  
      COMMON /A8/ DUP,DMID
      P=2*(X-DMID)/G
      IF (P.LT.0.D0) THEN
      P=0.D0
      ELSEIF (P.GT.1.D0) THEN
      P=1.D0
      ELSE
	P=P
      ENDIF
      IF (LBELOW.EQ.NL) THEN
      BGSP=CDSIN(M*DCMPLX(DACOS(P),0.D0))
     +  *(BBELOW*CDEXP(-CI*HBELOW*(X)))
     +  *(BBELOW*CDEXP(-CI*HBELOW*(X)))
      ELSEIF (LBELOW.EQ.1) THEN
      BGSP=CDSIN(M*DCMPLX(DACOS(P),0.D0))
     +  *(ABELOW*CDEXP(CI*HBELOW*(X)))
     +  *(ABELOW*CDEXP(CI*HBELOW*(X)))
      ELSE
      BGSP=CDSIN(M*DCMPLX(DACOS(P),0.D0)) 
     +   *(ABELOW*CDEXP(CI*HBELOW*X)
     +   +BBELOW*CDEXP(-CI*HBELOW*X))
     +   *(ABELOW*CDEXP(CI*HBELOW*X)
     +   +BBELOW*CDEXP(-CI*HBELOW*X)) 
C     +  *(ABELOW*CDCOS(HBELOW*(X))
C     +  +BBELOW*CDSIN(HBELOW*(X)))
C     +  *(ABELOW*CDCOS(HBELOW*(X))
C     +  +BBELOW*CDSIN(HBELOW*(X)))
      ENDIF            
      RETURN
      END      
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION CGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 G,X,A1,DX(500),DUP
      COMPLEX*16  AUP,BUP,HUP,CXUP,CI,CX,E(500)
      INTEGER M,LUP,LBELOW,NL,IC
      COMMON /RUP/ AUP,BUP,HUP,CXUP
      COMMON /period/ G,GL,M
      COMMON /A7/ LUP,LBELOW
      COMMON DX,E,NL,IC,CI
      COMMON /A8/ DUP,DMID        
      PI=3.1415926535897932D0
      A1=PI/(2.D0)        
      CX=(CDEXP(CI*M*A1)-CDEXP(-1.D0*CI*M*A1))
     +   *CDEXP(CI*M*(PI/G)*(X-DMID)) 
      IF (LUP.EQ.1) THEN 
      CGSP=CX*(AUP*CDEXP(CI*HUP*(X)))*(AUP*CDEXP(CI*HUP*(X)))
      ELSEIF (LUP.EQ.NL) THEN
      CGSP=CX*(BUP*CDEXP(-CI*HUP*X))*(BUP*CDEXP(-CI*HUP*X))
      ELSE
      CGSP=CX*(AUP*CDEXP(CI*HUP*X)+BUP*CDEXP(-CI*HUP*X))
     +	   *(AUP*CDEXP(CI*HUP*X)+BUP*CDEXP(-CI*HUP*X))
C      CGSP=CX*((AUP*CDCOS(HUP*(X))+BUP*CDSIN(HUP*(X)))*
C     +     (AUP*CDCOS(HUP*(X))+BUP*CDSIN(HUP*(X)))) 
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION DGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)                          
      REAL*8 DX(500),G,X,A1,DUP
      INTEGER NL,M,LUP,LBELOW,IC
      COMPLEX*16 ABELOW,BBELOW,HBELOW,CXBELOW,E(500),CI,CY
      COMMON /RBELOW/ ABELOW,BBELOW,HBELOW,CXBELOW
      COMMON DX,E,NL,IC,CI
      COMMON /period/ G,GL,M  
      COMMON /A7/ LUP,LBELOW   
      COMMON /A8/ DUP,DMID
      PI=3.1415926535897932D0
      A1=PI/2
      CY=(CDEXP(CI*M*3*A1)-CDEXP(CI*M*A1))
     +   *CDEXP(CI*PI*M*((X-DMID)/G)) 
      IF (LBELOW.EQ.1) THEN
      DGSP=CY*((ABELOW*CDEXP(CI*HBELOW*X))*(ABELOW*CDEXP(CI*HBELOW*X)))
      ELSEIF (LBELOW.EQ.NL) THEN   
      DGSP=CY*((BBELOW*CDEXP(-CI*HBELOW*(X)))
     +    *(BBELOW*CDEXP(-CI*HBELOW*(X))))
      ELSE
      DGSP=CY*(ABELOW*CDEXP(CI*HBELOW*X)+BBELOW*CDEXP(-CI*HBELOW*X))
     +       *(ABELOW*CDEXP(CI*HBELOW*X)+BBELOW*CDEXP(-CI*HBELOW*X))
C      DGSP=CY*((ABELOW*CDCOS(HBELOW*(X))+BBELOW*CDSIN(HBELOW*(X)))
C     +         *(ABELOW*CDCOS(HBELOW*(X))+BBELOW*CDSIN(HBELOW*(X))))
      ENDIF 
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION EGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 G,X,AW1,DX(500),DUP
      COMPLEX*16  AUP,BUP,HUP,CXUP,CI,E(500)
      INTEGER M,NL,IC
      COMMON /RUP/ AUP,BUP,HUP,CXUP
      COMMON /period/ G,GL,M
      COMMON /A7/ LUP,LBELOW
      COMMON DX,E,NL,IC,CI
      COMMON /A8/ DUP,DMID        
      PI=3.1415926535897932D0
      AW1=(X-DMID)/G
      CX=CDEXP(CI*M*PI)-CDEXP(CI*2.D0*PI*M*AW1)
      IF (LUP.EQ.1) THEN
      EGSP=CX*(AUP*CDEXP(CI*HUP*X))*(AUP*CDEXP(CI*HUP*X))
      ELSEIF (LUP.EQ.NL) THEN  
      EGSP=CX*(BUP*CDEXP(-CI*HUP*X))*(BUP*CDEXP(-CI*HUP*X))
      ELSE
      EGSP=CX*(AUP*CDEXP(HUP*X)+BUP*CDEXP(-CI*HUP*X))
     +       *(AUP*CDEXP(HUP*X)+BUP*CDEXP(-CI*HUP*X))
C      EGSP=CX*((AUP*CDCOS(HUP*X)+BUP*CDSIN(HUP*X))*
C     +        (AUP*CDCOS(HUP*X)+BUP*CDSIN(HUP*X)))
      ENDIF 
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION FGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)                          
      REAL*8 DX(500),G,X,AW4,DUP
      INTEGER NL,M,LUP,LBELOW,IC
      COMPLEX*16 ABELOW,BBELOW,HBELOW,CXBELOW,E(500),CI,CY
      COMMON /RBELOW/ ABELOW,BBELOW,HBELOW,CXBELOW
      COMMON DX,E,NL,IC,CI
      COMMON /period/ G,GL,M   
      COMMON /A7/ LUP,LBELOW 
      COMMON /A8/ DUP,DMID
      PI=3.1415926535897932D0
      AW4=(X-DMID)/G
      CY=(CDEXP(CI*2.d0*PI*M*AW4)-CDEXP(-CI*PI*M))
      IF (LBELOW.EQ.1) THEN
      FGSP=CY*((ABELOW*CDEXP(CI*HBELOW*X))*(ABELOW*CDEXP(CI*HBELOW*X))) 
      ELSEIF (LBELOW.EQ.NL) THEN 
      FGSP=CY*((BBELOW*CDEXP(-CI*HBELOW*X))
     +    *(BBELOW*CDEXP(-CI*HBELOW*X)))
      ELSE
      FGSP=CY*(ABELOW*CDEXP(CI*HBELOW*X)+BBELOW*CDEXP(-CI*HBELOW*X))
     +       *(ABELOW*CDEXP(CI*HBELOW*X)+BBELOW*CDEXP(-CI*HBELOW*X))
C      FGSP=CY*((ABELOW*CDCOS(HBELOW*X)
C     +     +BBELOW*CDSIN(HBELOW*X))*
C     +     (ABELOW*CDCOS(HBELOW*X)
C     +     +BBELOW*CDSIN(HBELOW*X)))
      ENDIF 
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
      DOUBLE COMPLEX FUNCTION GGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 AA1,AA2,BB1,BB2,HH1,HH2,CI   
      REAL*8 DIST
      COMMON /AC3/ AA1,BB1,HH1
      COMMON /AC3A/ AA2,BB2,HH2                
      COMMON /AC4/ DIST
	CI=DCMPLX(0.D0,1.D0)
C      GGSP=DCONJG(AA1*CDCOS(HH1*(X-DIST))+BB1*CDSIN(HH1*(X-DIST)))
C     +    *(AA2*CDCOS(HH2*X)+BB2*CDSIN(HH2*X))
C      GGSP=DCONJG(AA1*CDCOS(HH1*(X))+BB1*CDSIN(HH1*(X)))
C     +    *(AA2*CDCOS(HH2*X)+BB2*CDSIN(HH2*X))
      GGSP=DCONJG(AA1*CDEXP(CI*HH1*X)+BB1*CDEXP(-CI*HH1*X))
     +     *(AA2*CDEXP(CI*HH2*X)+BB2*CDEXP(-CI*HH2*X))
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   
      DOUBLE COMPLEX FUNCTION GGSPA(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 AA1,BB1,HH1,CI   
      COMMON /AC3/ AA1,BB1,HH1 
	CI=DCMPLX(0.D0,1.D0)
      GGSPA=DCONJG(AA1*CDEXP(CI*HH1*X)+BB1*CDEXP(-CI*HH1*X))
     +     *(AA2*CDEXP(CI*HH2*X)+BB2*CDEXP(-CI*HH2*X))
C      GGSPA=DCONJG(AA1*CDCOS(HH1*(X))+BB1*CDSIN(HH1*(X)))
C     +      *(AA1*CDCOS(HH1*(X))+BB1*CDSIN(HH1*(X)))
      RETURN
      END   
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                      
      DOUBLE COMPLEX FUNCTION HGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 AA1,AA2,BB1,BB2,HH1,HH2,CI 
      REAL*8 DIST
      COMMON /AC3/ AA1,BB1,HH1
      COMMON /AC3A/ AA2,BB2,HH2   
      COMMON /AC4/ DIST
	CI=DCMPLX(0.D0,1.D0)
C      HGSP=DCONJG(AA2*CDCOS(HH2*X)+BB2*CDSIN(HH2*X))
C     +    *(AA1*CDCOS(HH1*(X-DIST))+BB1*CDSIN(HH1*(X-DIST)))
      HGSP=DCONJG(AA2*CDEXP(CI*HH2*X)+BB2*CDEXP(-CI*HH2*X))
     +    *(AA1*CDEXP(CI*HH1*X)+BB1*CDEXP(-CI*HH1*X))
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE COMPLEX FUNCTION HGSPA(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 AA2,BB2,HH2,CI 
      COMMON /AC3A/ AA2,BB2,HH2
	CI=DCMPLX(0.D0,1.D0)   
      HGSPA=DCONJG(AA2*CDEXP(CI*HH2*X)+BB2*CDEXP(-CI*HH2*X))
     +    *(AA2*CDEXP(CI*HH2*X)+BB2*CDEXP(-CI*HH2*X))
      RETURN
      END
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                
      DOUBLE COMPLEX FUNCTION XXGSP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 AA1,AA2,BB1,BB2,HH1,HH2,CI 
      REAL*8 DUP
      COMMON /AC3/ AA1,BB1,HH1  
      COMMON /AC3A/ AA2,BB2,HH2
      COMMON /A8/ DUP
	CI=DCMPLX(0.D0,1.D0) 
      XXGSP=(AA2*CDEXP(CI*HH2*X)+BB2*CDEXP(-CI*HH2*X))
     +    *(AA1*CDEXP(CI*HH1*X)+BB1*CDEXP(-CI*HH1*X))
      RETURN
      END  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                       C
C     THE PROGRAM BRAGG.F IS FOR THE CALCULATION OF REFLICIVITY OF      C
C     BRAGG REFLECTOR IN THE APPLICATION OF VCSEL AND SURFACE EMITTING  C
C     LED.                                                              C
C                                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      SUBROUTINE BRAGG
      IMPLICIT REAL*8 (A-H,O-Z)
    1 CALL SELECTZ (N)
      GOTO (10,20,30,40), N
   10 CALL PERIODIC ! (1/4 LAMBDA) STACK
      GOTO 1
c  20 CALL TEM ! NORMAL INCIDENT
c     GOTO 1
   20 CALL TEZ  ! PERPENDICULAR INCIDENT (S-WAVE)
      GOTO 1
   30 CALL TMZ  ! PARALLEL INCIDENT (P-WAVE)
      GOTO 1
   40 PRINT*, ' THIS PROGRAM STOPS HERE !'
      RETURN 
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SELECTZ (N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
      PRINT*, '                                                 '
      PRINT*, ' PROGRAM PURPOSE -- BRAGG REFLECTOR DESIGN       '
      PRINT*, '                                                 '
      PRINT*, ' FOR NORMAL INCIDENT S-WAVE, P-WAVE WAVE AND THE '
      PRINT*, ' TEM WAVE HAVE THE SAME SOLUTIONS.               '
      PRINT*, '                                                 '
      WRITE(*,11)
   11 FORMAT(/' ENTER 1 FOR REFLECTIVITY OF (1/4 LAMBDA) STACK'/
c    +        '       2 FOR NORMAL INCIDENT'/ 
     +        '       2 FOR PERPENDICULAR INCIDENT (S-WAVE)'/ 
     +        '       3 FOR PARALLEL INCIDENT (P-WAVE)   '/ 
     +        '       4 FOR EXIT     ')
      PRINT*, '                                                 '
      PRINT*, ' THE SELECTION IS --> ?                          '
      PRINT*, '                                                 '
      READ(*,*) N
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PERIODIC
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 N,NO,NS,N1,N2
      PRINT*, '                                                  '
      PRINT*, ' FOR THE REFLECTIVITY OF PERIODIC DIELECTRIC STACK'
      PRINT*, ' (1/4 LAMBDA) MOST IN VCESL APPLICATION.          '
      PRINT*, '                                                  '
      PRINT*, ' INPUT THE REFRACTIVE INDEX OF INCIDENT LAYER NO --> ?'
      READ(*,*) NO
      PRINT*, ' INPUT THE REFRACTIVE INDEX OF SUBSTRATE LAYER NS -->?'
      READ(*,*) NS
      PRINT*, ' INPUT THE REFRACTIVE INDEX OF STACK LAYER 1 AND 2    '
      PRINT*, ' N1 -- > ?, N2 --> ?'
      READ(*,*) N1,N2
 1001 PRINT*, ' INPUT THE NUMBER OF PAIRS OF STACK N -- >?'
      READ(*,*) N
      XNU=(1.D0-(NS/NO)*(N1/N2)**(2.D0*N))
      XDE=(1.D0+(NS/NO)*(N1/N2)**(2.D0*N))
      R=(XNU/XDE)*(XNU/XDE)
      PRINT*, ' REFLECTIVITY R= -->', R
      PRINT*, ' FOR NEW CALCULATION INPUT --> 1, FOR EXIT INPUT --> 2'
      PRINT*, ' INPUT -- >? '
      READ(*,*) I
      IF (I.EQ.1 ) THEN
      GOTO 1001
      ELSEIF (I.EQ.2) THEN
      GOTO 1002
      ELSE
      ENDIF
 1002 PRINT*, ' PROGRAM STOP HERE !'
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE TEM
      IMPLICIT REAL*8 (A-H,O-Z)
      PRINT*, '                                                '
      PRINT*, ' FOR EM WAVE NORMAL INCIDENT INTO THE DIELECTRIC'
      PRINT*, ' MOST IN VCESL APPLICATION.                     '
      PRINT*, '                                                '
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE TEZ
      PARAMETER (IN=100)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 D(2,2,0:IN),P(2,2,0:IN),M(2,2)
      COMPLEX*16 XD(2,2),CNL(0:IN),CKNL(0:IN),DX(2,2),PX(2,2),CI
      COMPLEX*16 A(2,2),B(2,2),BB(2,2,0:IN),PT(2,2),T(2,2)
      COMPLEX*16 XT(2,2),DD(2,2),DT(2,2),BA(2,2),R,TT
      REAL*8 N(0:IN),ALPHA(0:IN),W(0:IN),THETAL(0:IN),LAMBDA,LAMBDA1
      INTEGER NL
      CHARACTER IFILE*40
      PRINT*, '                                                       '
      PRINT*, ' FOR EM WAVE PERPENDICULAR INCIDENT INTO THE DIELECTRIC'
      PRINT*, ' ITS ALSO CALLED TE POLARIZATION OR S-WAVE.           '
      PRINT*, '                                                       '
      PRINT*, ' SET INCIDENT LAYER IS 0th LAYER '
      PRINT*, ' SET FINAL LAYER IS N+1th LAYER, N LAYERS IN BETWEEN'
      PRINT*, '                                                       '
      OPEN(2,FILE='reflect.dat',STATUS='UNKNOWN')
      MIN=5
      PI=3.141592653589D0
      CI=DCMPLX(0.D0,1.D0)
      PRINT*,' THE INPUT FILE NAME='
      READ(MIN,'(A40)') IFILE
      OPEN(1,FILE=IFILE,status='unknown')
      READ(1,*) NL
      READ(1,*) LAMBDA
      DO I=0,NL+1
	 READ(1,*) N(I),ALPHA(I),W(I)
	 WRITE(*,7776) N(I),ALPHA(I),W(I),I
         CNL(I)=DCMPLX(N(I),ALPHA(I))
 7776    FORMAT(' N=',D12.6,' ALPHA=',D12.6,' W=',D14.7,' I=',I3)
      ENDDO
C
      XLL=LAMBDA-0.5
      XLR=LAMBDA+0.5
      PRINT*, '                                                       '
      PRINT*, ' INPUT INCIDENT ANGLE THETAL -->?'
      PRINT*, '                                                       '
      READ(*,*) ANGLE
      THETAL(0)=ANGLE*PI/(180.D0)               ! RADIANS
      DO LAMBDA1=XLL,XLR,0.01
      CKNL(0)=CNL(0)*(2*PI/LAMBDA1)*COS(THETAL(0))
      DO I=1,NL+1
         THETAL(I)=ASIN(SIN(THETAL(I-1)*(N(I-1)/N(I))))
         CKNL(I)=CNL(I)*(2*PI/LAMBDA1)*COS(THETAL(I))
      ENDDO
C
      D(1,1,0)=DCMPLX(1.D0,0.D0)  ! MATRIX D0
      D(1,2,0)=DCMPLX(1.D0,0.D0)
      D(2,1,0)=CNL(0)*COS(THETAL(0))
      D(2,2,0)=-1.D0*CNL(0)*COS(THETAL(0))
C
      DD(1,1)=D(1,1,0)
      DD(1,2)=D(1,2,0)
      DD(2,1)=D(2,1,0)
      DD(2,2)=D(2,2,0)
      CALL ZLINEAR2 (DD,2,2,DET,INFO,01)
      D(1,1,0)=DD(1,1)                      ! MATRIX D0'
      D(1,2,0)=DD(1,2)
      D(2,1,0)=DD(2,1)
      D(2,2,0)=DD(2,2)
C
      DT(1,1)=DCMPLX(1.D0,0.D0) ! MATRIX DS
      DT(1,2)=DCMPLX(1.D0,0.D0)
      DT(2,1)=CNL(NL+1)*COS(THETAL(NL+1))
      DT(2,2)=-1.D0*CNL(NL+1)*COS(THETAL(NL+1))
C
      DO I=1,NL
         D(1,1,I)=DCMPLX(1.D0,0.D0)
         D(1,2,I)=DCMPLX(1.D0,0.D0)
         D(2,1,I)=CNL(I)*COS(THETAL(I))
         D(2,2,I)=-1.D0*CNL(I)*COS(THETAL(I))
	 DO J=1,2
	    DO K=1,2
	       XD(J,K)=D(J,K,I)
	       DX(J,K)=D(J,K,I)
            ENDDO
         ENDDO
         CALL ZLINEAR2 (XD,2,2,DET,INFO,01)
         P(1,1,I)=CDEXP(CI*CKNL(I)*W(I)) ! MATRIX PL
         P(2,2,I)=CDEXP(-1.D0*CI*CKNL(I)*W(I))
         P(1,2,I)=DCMPLX(0.D0,0.D0)
         P(2,1,I)=DCMPLX(0.D0,0.D0)
	 DO J=1,2
	    DO K=1,2
	       PX(J,K)=P(J,K,I)
            ENDDO
         ENDDO
         CALL MULTI (2,PX,XD,A)  ! MATRIX PL*DL'
         CALL MULTI (2,DX,A,B)   ! MATRIX DL*PL*DL'
         BB(1,1,I)=B(1,1)
         BB(1,2,I)=B(1,2)
         BB(2,1,I)=B(2,1)
         BB(2,2,I)=B(2,2)
      ENDDO
      T(1,1)=BB(1,1,NL)
      T(1,2)=BB(1,2,NL)
      T(2,1)=BB(2,1,NL)
      T(2,2)=BB(2,2,NL)
      DO K=NL-1,1,-1
	 DO I=1,2
	    DO J=1,2
	       BA(I,J)=BB(I,J,K)
            ENDDO
         ENDDO
	 CALL MULTI (2,BA,T,PT)
	 DO I=1,2
	    DO J=1,2
	       T(I,J)=PT(I,J)
            ENDDO
         ENDDO
      ENDDO 
      CALL MULTI (2,T,DT,XT)
      M(1,1)=D(1,1,0)*XT(1,1)+D(1,2,0)*XT(2,1)
      M(1,2)=D(1,1,0)*XT(1,2)+D(1,2,0)*XT(2,2)
      M(2,1)=D(2,1,0)*XT(1,1)+D(2,2,0)*XT(2,1)
      M(2,2)=D(2,1,0)*XT(1,2)+D(2,2,0)*XT(2,2) 
      R=(M(2,1)/M(1,1))*DCONJG(M(2,1)/M(1,1))
      TX=(1/(M(1,1)))*DCONJG(1/(M(1,1)))
      TT=(CNL(NL+1)*COS(THETAL(NL+1))/(CNL(0)*COS(THETAL(0))))*TX
C      PRINT*, ' R=',R,' TT=',TT
      WRITE(2,10) LAMBDA1,R
   10 FORMAT(2X,F10.5,2X,2(E16.9))
      ENDDO
      CLOSE (1)
      CLOSE (2)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      SUBROUTINE TMZ
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IN=100)
      COMPLEX*16 D(2,2,0:IN),P(2,2,0:IN),M(2,2)
      COMPLEX*16 XD(2,2),CNL(0:IN),CKNL(0:IN),DX(2,2),PX(2,2),CI
      COMPLEX*16 A(2,2),B(2,2),BB(2,2,0:IN),PT(2,2),T(2,2)
      COMPLEX*16 XT(2,2),DD(2,2),DT(2,2),R,TT
      REAL*8 N(0:IN),ALPHA(0:IN),W(0:IN),THETAL(0:IN),LAMBDA,LAMBDA1
      INTEGER NL
      CHARACTER IFILE*40
      PRINT*, '                                                       '
      PRINT*, ' FOR EM WAVE PARALLEL INCIDENT INTO THE DIELECTRIC     '
      PRINT*, ' ITS ALSO CALLED TM POLARIZATION OR P-WAVE.           '
      PRINT*, '                                                       '
      PRINT*, ' SET INCIDENT LAYER IS 0th LAYER '
      PRINT*, ' SET FINAL LAYER IS N+1th LAYER, N LAYERS IN BETWEEN'
      PRINT*, '                                                       '
      OPEN(2,FILE='reflect.dat',STATUS='UNKNOWN')
      MIN=5
      PI=3.141592653589D0
      CI=DCMPLX(0.D0,1.D0)
      PRINT*,' THE INPUT FILE NAME='
      READ(MIN,'(A40)') IFILE
      OPEN(1,FILE=IFILE,status='unknown')
      READ(1,*) NL
      READ(1,*) LAMBDA
      DO I=0,NL+1
	 READ(1,*) N(I),ALPHA(I),W(I)
	 WRITE(*,7776) N(I),ALPHA(I),W(I),I
         CNL(I)=DCMPLX(N(I),ALPHA(I))
 7776    FORMAT(' N=',D12.6,' ALPHA=',D12.6,' W=',D14.7,' I=',I3)
      ENDDO
C
      XLL=LAMBDA-0.5
      XLR=LAMBDA+0.5
      PRINT*, '                                                       '
      PRINT*, ' INPUT INCIDENT ANGLE THETAL -->?'
      PRINT*, '                                                       '
      READ(*,*) ANGLE
      THETAL(0)=ANGLE*PI/(180.D0)               ! RADIANS
      DO LAMBDA1=XLL,XLR,0.01
      CKNL(0)=CNL(0)*(2*PI/LAMBDA1)*COS(THETAL(0))
      DO I=1,NL+1
         THETAL(I)=ASIN(SIN(THETAL(I-1)*(N(I-1)/N(I))))
         CKNL(I)=CNL(I)*(2*PI/LAMBDA1)*COS(THETAL(I))
      ENDDO
C
      D(1,1,0)=COS(THETAL(0))  ! MATRIX D0
      D(1,2,0)=COS(THETAL(0))
      D(2,1,0)=CNL(0)
      D(2,2,0)=-1.D0*CNL(0)
C
      DD(1,1)=D(1,1,0)
      DD(1,2)=D(1,2,0)
      DD(2,1)=D(2,1,0)
      DD(2,2)=D(2,2,0)
      CALL ZLINEAR2 (DD,2,2,DET,INFO,01)
      D(1,1,0)=DD(1,1)                      ! MATRIX D0'
      D(1,2,0)=DD(1,2)
      D(2,1,0)=DD(2,1)
      D(2,2,0)=DD(2,2)
C
      DT(1,1)=COS(THETAL(NL+1)) ! MATRIX DS
      DT(1,2)=COS(THETAL(NL+1))
      DT(2,1)=CNL(NL+1)
      DT(2,2)=-1.D0*CNL(NL+1)
C
      DO I=1,NL
         D(1,1,I)=COS(THETAL(I))
         D(1,2,I)=COS(THETAL(I))
         D(2,1,I)=CNL(I)
         D(2,2,I)=-1.D0*CNL(I)
         XD(1,1)=D(1,1,I)
         XD(1,2)=D(1,2,I)
         XD(2,1)=D(2,1,I)
         XD(2,2)=D(2,2,I)
         DX(1,1)=D(1,1,I)    ! MATRIX DL
         DX(1,2)=D(1,2,I)
         DX(2,1)=D(2,1,I)
         DX(2,2)=D(2,2,I) 
         CALL ZLINEAR2 (XD,2,2,DET,INFO,01)
         P(1,1,I)=CDEXP(CI*CKNL(I)*W(I)) ! MATRIX PL
         P(2,2,I)=CDEXP(-1.D0*CI*CKNL(I)*W(I))
         P(1,2,I)=DCMPLX(0.D0,0.D0)
         P(2,1,I)=DCMPLX(0.D0,0.D0)
         PX(1,1)=P(1,1,I)
         PX(1,2)=P(1,2,I)
         PX(2,1)=P(2,1,I)
         PX(2,2)=P(2,2,I)
         CALL MULTI (2,PX,XD,A)  ! MATRIX PL*DL'
         CALL MULTI (2,DX,A,B)   ! MATRIX DL*PL*DL'
         BB(1,1,I)=B(1,1)
         BB(1,2,I)=B(1,2)
         BB(2,1,I)=B(2,1)
         BB(2,2,I)=B(2,2)
      ENDDO
      T(1,1)=BB(1,1,NL)
      T(1,2)=BB(1,2,NL)
      T(2,1)=BB(2,1,NL)
      T(2,2)=BB(2,2,NL)
      DO I=NL-1,1,-1
         PT(1,1)=BB(1,1,I)*T(1,1)+BB(1,2,I)*T(2,1)
         PT(1,2)=BB(1,1,I)*T(1,2)+BB(1,2,I)*T(2,2)
         PT(2,1)=BB(2,1,I)*T(1,1)+BB(2,2,I)*T(2,1)
         PT(2,2)=BB(2,1,I)*T(1,2)+BB(2,2,I)*T(2,2)
         T(1,1)=PT(1,1)
         T(2,1)=PT(2,1)
         T(1,2)=PT(1,2)
         T(2,2)=PT(2,2)
      ENDDO 
      CALL MULTI (2,T,DT,XT)
      M(1,1)=D(1,1,0)*XT(1,1)+D(1,2,0)*XT(2,1)
      M(1,2)=D(1,1,0)*XT(1,2)+D(1,2,0)*XT(2,2)
      M(2,1)=D(2,1,0)*XT(1,1)+D(2,2,0)*XT(2,1)
      M(2,2)=D(2,1,0)*XT(1,2)+D(2,2,0)*XT(2,2) 
      R=(M(2,1)/M(1,1))*DCONJG(M(2,1)/M(1,1))
      TX=(1/(M(1,1)))*DCONJG(1/(M(1,1)))
      TT=(CNL(NL+1)*COS(THETAL(NL+1))/(CNL(0)*COS(THETAL(0))))*TX
C      PRINT*, ' R=',R,' TT=',TT
      WRITE(2,10) LAMBDA1,R
   10 FORMAT(2X,F10.5,2X,2(E16.9))
      ENDDO
      CLOSE (1)
      CLOSE (2)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE MULTI (N,A,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(N,N),C(N,N),D(N,N),SUM
      DO 10 I=1,N
         DO 20 J=1,N
            SUM=DCMPLX(0.D0,0.D0)
            DO 30 K=1,N
               SUM=SUM+A(I,K)*C(K,J)
   30 CONTINUE
      D(I,J)=SUM
   20 CONTINUE
   10 CONTINUE
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE ZLINEAR2(A,NA,NB,DET,INFO,II)
      COMPLEX*16 A(NA,NB),DET,XDET(2)
      INTEGER JOB,INFO
c     CALL ZGECO (A,NA,NB,IPVT,RCOND,Z)
      CALL XZGEFA (A,NA,NB,IPVT,INFO)
	IF (II.EQ.1) THEN
      CALL XZGEDI (A,NA,NB,IPVT,XDET,WORK,01)
	ELSEIF (II.EQ.2) THEN
      CALL XZGEDI (A,NA,NB,IPVT,XDET,WORK,10)      ! 01 inverse only
      DET=XDET(1)*10.D0**XDET(2)              ! 11 both can calculate
	ELSE
	ENDIF
      RETURN
      END
c
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine xzgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(n),job
      complex*16 a(lda,n),det(2),work(n)
c
c     zgedi computes the determinant and inverse of a matrix
c     using the factors computed by zgeco or zgefa.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the output from zgeco or zgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from zgeco or zgefa.
c
c        work    complex*16(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     complex*16(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if zgeco has set rcond .gt. 0.0 or zgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,zswap
c     fortran dabs,dcmplx,mod
c
c     internal variables
c
      complex*16 t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = (1.0d0,0.0d0)
         det(2) = (0.0d0,0.0d0)
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (cabs1(det(1)) .eq. 0.0d0) go to 60
   10       if (cabs1(det(1)) .ge. 1.0d0) go to 20
               det(1) = dcmplx(ten,0.0d0)*det(1)
               det(2) = det(2) - (1.0d0,0.0d0)
            go to 10
   20       continue
   30       if (cabs1(det(1)) .lt. ten) go to 40
               det(1) = det(1)/dcmplx(ten,0.0d0)
               det(2) = det(2) + (1.0d0,0.0d0)
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = (1.0d0,0.0d0)/a(k,k)
            t = -a(k,k)
            call zscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = (0.0d0,0.0d0)
               call zaxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = (0.0d0,0.0d0)
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call zaxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call zswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end



      double precision function dzasum(n,zx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double complex zx(n)
      double precision stemp,dcabs1
      integer i,incx,ix,n
c
      dzasum = 0.0d0
      stemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        stemp = stemp + dcabs1(zx(ix))
        ix = ix + incx
   10 continue
      dzasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        stemp = stemp + dcabs1(zx(i))
   30 continue
      dzasum = stemp
      return
      end




      subroutine  zdscal(n,da,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double complex zx(n)
      double precision da
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = dcmplx(da,0.0d0)*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = dcmplx(da,0.0d0)*zx(i)
   30 continue
      return
      end



      double complex function zdotc(n,zx,incx,zy,incy)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c
      double complex zx(n),zy(n),ztemp
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotc = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zy(i)
   30 continue
      zdotc = ztemp
      return
      end



      subroutine zaxpy(n,za,zx,incx,zy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c
      double complex zx(n),zy(n),za
      double precision dcabs1
      if(n.le.0)return
      if (dcabs1(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end




      integer function izamax(n,zx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, 1/15/85.
c     modified 3/93 to return if incx .le. 0.
c
      double complex zx(n)
      double precision smax
      integer i,incx,ix,n
      double precision dcabs1
c
      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = dcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(zx(ix)).le.smax) go to 5
         izamax = i
         smax = dcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dcabs1(zx(1))
      do 30 i = 2,n
         if(dcabs1(zx(i)).le.smax) go to 30
         izamax = i
         smax = dcabs1(zx(i))
   30 continue
      return
      end




      subroutine  zscal(n,za,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double complex za,zx(n)
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end



      subroutine  zswap (n,zx,incx,zy,incy)
c
c     interchanges two vectors.
c     jack dongarra, 3/11/78.
c
      double complex zx(n),zy(n),ztemp
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end



      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
c


      subroutine xzgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(n),info
      complex*16 a(lda,n)
c
c     zgefa factors a complex*16 matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax
c     fortran dabs
c
c     internal variables
c
      complex*16 t
      integer izamax,j,k,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = izamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
      return
      end
*
***********************************************************************
* 
      subroutine wx(z,w)
      implicit real*8 (a-h,o-z)
      complex*16 w,ci,x,xk,xz,z,csinh,ccosh,r1,r2,rp,rm
      real*8 xy,xx,xl
      integer icase
      common /a6/ r1,r2,xl,xkab
      common /d1/ xx,icase
      common /inx/ ixx
      ci=dcmplx(0.d0,1.d0)
      xk=xx-(ci*z) !delta(beta*L)
      xz=cdsqrt(xk**2-xkab**2) !q*L
      if (ixx.eq.1) then
      rp=xkab/(xk-xz) ! rp=xkab/(xk+xz)  xx > 0
      rm=xkab/(xk-xz) ! rm=xkab/(xk+xz)  xx > 0
      elseif (ixx.eq.2) then
      rp=xkab/(xk+xz)
      rm=xkab/(xk+xz)
      else
      endif
      if (icase.eq.1) then
      w=xk+ci*xz*(cdcos(xz)/(cdsin(xz)))
      elseif (icase.eq.2) then
      w=(1-r1*rp)*(1-r2*rm)+(r2-rp)*(rm-r1)
     + *(cdcos(2.*xz)+ci*cdcos(2.*xz))
      else
      w=r1*r2*(cdcos(2.*xz)+ci*cdcos(2.*xz))-1
      endif
      return
      end

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE XCROOT(FCN,Z,ZA,W4,I,J,EPS,EPSS)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,W4,DELX,H1,H2,DISC
      COMPLEX*16 Z1,Z2,Z3,W1,W2,W3,G,A,B,C  
      REAL*8 EPS,EPSS,ZA 
      INTEGER I,M
	COMMON /A11/ AIN
C**********
      M=500
      Z1=ZA-DCMPLX(AIN,0.D0)
      Z2=ZA
      Z3=ZA+DCMPLX(AIN,0.D0)
      CALL FCN(Z1,W1)
      CALL FCN(Z2,W2)
      CALL FCN(Z3,W3)
C     BEGIN ITERATIONS
      DO  J=1,M
          H1=Z2-Z1
          H2=Z3-Z2
          G=H1/H2  
          A=(W3*G-W2*(1.D0+G)+W1)/(H1*(H1+H2))
          B=(W3-W2-A*H2*H2)/H2
          C=W2
          DISC=CDSQRT(B*B-4.D0*A*C)
          IF (CDABS(B-DISC).LT.CDABS(B+DISC)) THEN
	     DISC=DISC
          ELSE
	     DISC=-DISC
          ENDIF
             DELX=-2.D0*C/(B+DISC)
             Z=Z2+DELX
             CALL FCN (Z,W4)
C       IF (I.EQ.0) PRINT 199 ,J,Z,W4
             IF (CDABS(DELX).LE.EPS) THEN
                I=1
C       PRINT 202,J,Z,W4
             RETURN
             ENDIF
             IF (CDABS(W4).LE.EPSS) THEN
                I=2
C       PRINT 203,J,Z,W4
             RETURN
             ENDIF
 	     IF (CDABS(DELX).GE.0.D0) THEN
                Z1=Z2
                W1=W2
                IF (CDABS(DELX).GT.CDABS(H2)) THEN
                   Z2=Z3
                   W2=W3
                   Z3=Z
                   W3=W4
                ELSE
                   Z2=Z
                   W2=W4
                ENDIF
             ELSE
                Z3=Z2
                W3=W2
                IF (CDABS(DELX).GT.CDABS(H1)) THEN
                   Z2=Z1
                   W2=W1
                   Z1=Z
                   W1=W4
                ELSE
                   Z2=Z
                   W2=W4
                ENDIF
             ENDIF
      ENDDO
      I=-1
c      PRINT 200 ,M,Z,W4
C      zr=z
      RETURN
C  199 FORMAT(' AT ITERATION',I3,3X,'Z=',2E12.5,4X,'W4=',2E12.5)   
c  200 FORMAT(/' TOLERANCE NOT MET AFTER',I4,' ITERATIONS X=',
c     +       2E18.8,2X,'W4=',2E18.8)
C  202 FORMAT(/' Z TOLERANCE MET IN ',I4,' ITERATIONS'/ 
C     +     2X,' Z=',2E16.8,2X,'W4=',2E16.8)  
C  203 FORMAT(/' W TOLERANCE MET IN ',I4,' ITERATIONS'/
C     +     2X,' Z=',2E16.8,2X,'W4=',2E16.8)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C            
      SUBROUTINE CSPLINE (X1,X2,X3,Y1,Y2,Y3,SOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,I
      REAL*8 X(500),Y(500),A(500),B(500),C(500),S(500)
      DATA N/3/
      DATA S/500*0.D0/  
      DO I=1,500
      X(I)=0.D0
      Y(I)=0.D0
      ENDDO 
      X(1)=X1
      X(2)=X2
      X(3)=X3
      Y(1)=Y1
      Y(2)=Y2
      Y(3)=Y3
C
      CALL CUBSPL3(X,Y,N,A,B,C,S)
      DO U=0.48D0,0.5201D0,0.001D0 
      TEMP1=SEVAL(N,U,X,Y,A,B,C)
c      IF ((U-0.5000D0).LT.1.E-7) SOL=TEMP1
      IF (U.EQ.0.5000D0) SOL=TEMP1
      ENDDO                           
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CUBSPL3(X,Y,N,A,B,C,S)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(N),Y(N),S(N),A(N),B(N),C(N),SMATRIX(0:100,4)
      REAL*8 DX1,DY1,DX2,DY2,DXN1,DXN2,H(100)
      INTEGER N,NM1,NM2,I,J,FIRST,LAST
      NM2=N-2
      NM1=N-1
      DX1=X(2)-X(1)
      DY1=(Y(2)-Y(1))/DX1*6.D0
      DO I=1,NM2
         DX2=X(I+2)-X(I+1)
         DY2=(Y(I+2)-Y(I+1))/DX2*6.D0
         SMATRIX(I,1)=DX1
         SMATRIX(I,2)=2.D0*(DX1+DX2)
         SMATRIX(I,3)=DX2
         SMATRIX(I,4)=DY2-DY1
         DX1=DX2
         DY1=DY2
      ENDDO
      FIRST=2
      LAST=NM2
   80 DX1=X(2)-X(1)
      DX2=X(3)-X(2)
      SMATRIX(1,2)=(DX1+DX2)*(DX1+2.D0*DX2)/DX2
      SMATRIX(1,3)=(DX2*DX2-DX1*DX1)/DX2
      DXN2=X(NM1)-X(NM2)
      DXN1=X(N)-X(NM1)
      SMATRIX(NM2,1)=(DXN2*DXN2-DXN1*DXN1)/DXN2
      SMATRIX(NM2,2)=(DXN1+DXN2)*(DXN1+2.D0*DXN2)/DXN2
      DO I=FIRST-1,LAST
      ENDDO
      DO I=FIRST,LAST
         SMATRIX(I,1)=SMATRIX(I,1)/SMATRIX(I-1,2)
         SMATRIX(I,2)=SMATRIX(I,2)-SMATRIX(I,1)*SMATRIX(I-1,3)
         SMATRIX(I,4)=SMATRIX(I,4)-SMATRIX(I,1)*SMATRIX(I-1,4)
      ENDDO
      SMATRIX(LAST,4)=SMATRIX(LAST,4)/SMATRIX(LAST,2)
      DO J=LAST-1,FIRST-1,-1
         SMATRIX(J,4)=(SMATRIX(J,4)-SMATRIX(J,3)*SMATRIX(J+1,4))
     +   /SMATRIX(J,2)
      ENDDO
      DO I=FIRST-1,LAST
         S(I+1)=SMATRIX(I,4)
      ENDDO
      S(1)=((DX1+DX2)*S(2)-DX1*S(3))/DX2
      S(N)=((DXN2+DXN1)*S(NM1)-DXN1*S(NM2))/DXN2
      DO I=1,N-1
         H(I)=X(I+1)-X(I)
         A(I)=(S(I+1)-S(I))/(6.D0*H(I))
         B(I)=S(I)/2.D0
         C(I)=((Y(I+1)-Y(I))/H(I))-((2*H(I)*S(I)+H(I)*S(I+1))/6.D0)
      ENDDO
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SEVAL (N,U,X,Y,A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 U,X(N),Y(N),A(N),B(N),C(N)
      INTEGER I,J,K
      REAL*8 DX
      SAVE I
      DATA I/1/
      IF (I.GE.N) THEN
          THEN I=1
      ENDIF
      IF(U.GE.X(N)) THEN
              DX=U-X(N-1)
              SEVAL=Y(N-1)+DX*(C(N-1)+DX*(B(N-1)+DX*A(N-1)))
              RETURN
      ENDIF
      IF(U.GE.X(I)) THEN
              IF(U.LE.X(I+1)) THEN
              DX=U-X(I)
              SEVAL=Y(I)+DX*(C(I)+DX*(B(I)+DX*A(I)))
              RETURN
      ENDIF
      ENDIF
      I=1
      J=N+1
   10 K=(I+J)/2
      IF(U.LT.X(K)) THEN
              J=K
      ELSE
              I=K
      ENDIF
      IF(J.GT.I+1) THEN
              GOTO 10
      ELSE
      DX=U-X(I)
      SEVAL=Y(I)+DX*(C(I)+DX*(B(I)+DX*A(I)))
      RETURN
      ENDIF
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       
      SUBROUTINE CROOT(FCN,Z,I,J,EPS,EPSS)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,W4,DISC,DELX,H1,H2,Z1,Z2,Z3,W1,W2,W3,G,A,B,C
      REAL*8 EPS,EPSS,AIN
C*****
      INTEGER M,I
	COMMON /A11/ AIN   
      M=500
      Z1=Z-DCMPLX(AIN,0.d0)
      Z2=Z
      Z3=Z+DCMPLX(AIN,0.d0)
      CALL FCN(Z1,W1)
      CALL FCN(Z2,W2)
      CALL FCN(Z3,W3)
C     BEGIN ITERATIONS
      DO  J=1,M
          H1=Z2-Z1
          H2=Z3-Z2
          G=H1/H2  
          A=(W3*G-W2*(DCMPLX(1.D0,0.D0)+G)+W1)/(H1*(H1+H2))
          B=(W3-W2-A*H2*H2)/H2
          C=W2        
          DISC=CDSQRT(B*B-(DCMPLX(4.D0,0.D0)*A*C))
          IF (CDABS(B-DISC).LT.CDABS(B+DISC)) THEN
	     DISC=DISC
          ELSE
	     DISC=-DISC
          ENDIF
             DELX=-2.D0*C/(B+DISC)
             Z=Z2+DELX
             CALL FCN (Z,W4)
C      IF (I.EQ.0) PRINT 199 ,J,Z,W4
             IF (CDABS(DELX).LE.EPS) THEN
                I=1
C      PRINT 202,J,Z,W4
             RETURN
             ENDIF
             IF (CDABS(W4).LE.EPSS) THEN
                I=2
C      PRINT 203,J,Z,W4
             RETURN
             ENDIF
c             IF ((Dreal(DELX).GE.0.D0).AND.(DIMAG(DELX).GE.0.D0)) THEN
              IF (Dreal(DELX).GE.0.D0) THEN
                Z1=Z2
                W1=W2
              IF (CDABS(DELX).GT.CDABS(H2)) THEN
c              IF (DREAL(DELX).GT.DREAL(H2)) THEN
                   Z2=Z3
                   W2=W3
                   Z3=Z
                   W3=W4
                ELSE
                   Z2=Z
                   W2=W4
                ENDIF
             ELSE
                Z3=Z2
                W3=W2
                IF (CDABS(DELX).GT.CDABS(H1)) THEN
                   Z2=Z1
                   W2=W1
                   Z1=Z
                   W1=W4
                ELSE
                   Z2=Z
                   W2=W4
                ENDIF
             ENDIF
      ENDDO
      I=-1
      PRINT 200 ,M,Z,W4
      RETURN
C  199 FORMAT(' AT ITERATION',I3,3X,'Z=',2E12.5,4X,'W4=',2E12.5)
  200 FORMAT(/' TOLERANCE NOT MET AFTER',I4,' ITERATIONS X=',
     +       2E18.8,2x,'W4=',2E18.8)
C  202 FORMAT(/' TOLERANCE MET IN ',I4,' ITERATIONS'/ 
C     +     2X,'Z=',2E16.8,2X,'W4=',2E16.8)
C  203 FORMAT(/' W TOLERANCE MET IN ',I4,' ITERATIONS Z=',2E12.5, 
C     +      ' W4=',2E12.5)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine quanc8(fun,a,b,abserr,relerr,result,errest,nofun,flag)
      double precision fun, a, b, abserr, relerr, result, errest, flag
      integer nofun
c
c
      double precision w0,w1,w2,w3,w4,area,x0,f0,stone,step,cor11,temp
      double precision qprev,qnow,qdiff,qleft,esterr,tolerr
      double precision qright(31),f(16),x(16),fsave(8,30),xsave(8,30)
      double precision dabs,dmax1
      integer levmin,levmax,levout,nomax,nofin,lev,nim,i,j
c
c   ***   stage 1 ***   general initialization
c   set constants.
c
      levmin = 1
      levmax = 30
      levout = 6
      nomax = 5000
      nofin = nomax - 8*(levmax-levout+2**(levout+1))
c
c   trouble when nofun reaches nofin
c
      w0 =   3956.0d0 / 14175.0d0
      w1 =  23552.0d0 / 14175.0d0
      w2 =  -3712.0d0 / 14175.0d0
      w3 =  41984.0d0 / 14175.0d0
      w4 = -18160.0d0 / 14175.0d0
c
c   initialize running sums to zero.
c
      flag = 0.0d0
      result = 0.0d0
      cor11  = 0.0d0
      errest = 0.0d0
      area   = 0.0d0
      nofun = 0
      if (a .eq. b) return
c
c   ***   stage 2 ***   initialization for first interval
c
      lev = 0
      nim = 1
      x0 = a
      x(16) = b
      qprev  = 0.0d0
      f0 = fun(x0)
      stone = (b - a) / 16.0d0
      x(8)  =  (x0  + x(16)) / 2.0d0
      x(4)  =  (x0  + x(8))  / 2.0d0
      x(12) =  (x(8)  + x(16)) / 2.0d0
      x(2)  =  (x0  + x(4))  / 2.0d0
      x(6)  =  (x(4)  + x(8))  / 2.0d0
      x(10) =  (x(8)  + x(12)) / 2.0d0
      x(14) =  (x(12) + x(16)) / 2.0d0
      do 25 j = 2, 16, 2
         f(j) = fun(x(j))
   25 continue
      nofun = 9
c
c   ***   stage 3 ***   central calculation
c   requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16.
c   calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area.
c
   30 x(1) = (x0 + x(2)) / 2.0d0
      f(1) = fun(x(1))
      do 35 j = 3, 15, 2
         x(j) = (x(j-1) + x(j+1)) / 2.0d0
         f(j) = fun(x(j))
   35 continue
      nofun = nofun + 8
      step = (x(16) - x0) / 16.0d0
      qleft  =  (w0*(f0 + f(8))  + w1*(f(1)+f(7))  + w2*(f(2)+f(6))
     1  + w3*(f(3)+f(5))  +  w4*f(4)) * step
      qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14))
     1  + w3*(f(11)+f(13)) + w4*f(12)) * step
      qnow = qleft + qright(lev+1)
      qdiff = qnow - qprev
      area = area + qdiff
c
c   ***   stage 4 *** interval convergence test
c
      esterr = dabs(qdiff) / 1023.0d0
      tolerr = dmax1(abserr,relerr*dabs(area)) * (step/stone)
      if (lev .lt. levmin) go to 50
      if (lev .ge. levmax) go to 62
      if (nofun .gt. nofin) go to 60
      if (esterr .le. tolerr) go to 70
c
c   ***   stage 5   ***   no convergence
c   locate next interval.
c
   50 nim = 2*nim
      lev = lev+1
c
c   store right hand elements for future use.
c
      do 52 i = 1, 8
         fsave(i,lev) = f(i+8)
         xsave(i,lev) = x(i+8)
   52 continue
c
c   assemble left hand elements for immediate use.
c
      qprev = qleft
      do 55 i = 1, 8
         j = -i
         f(2*j+18) = f(j+9)
         x(2*j+18) = x(j+9)
   55 continue
      go to 30
c
c   ***   stage 6   ***   trouble section
c   number of function values is about to exceed limit.
c
   60 nofin = 2*nofin
      levmax = levout
      flag = flag + (b - x0) / (b - a)
      go to 70
c
c   current level is levmax.
c
   62 flag = flag + 1.0d0
c
c   ***   stage 7   ***   interval converged
c   add contributions into running sums.
c
   70 result = result + qnow
      errest = errest + esterr
      cor11  = cor11  + qdiff / 1023.0d0
c
c   locate next interval.
c
   72 if (nim .eq. 2*(nim/2)) go to 75
      nim = nim/2
      lev = lev-1
      go to 72
   75 nim = nim + 1
      if (lev .le. 0) go to 80
c
c   assemble elements required for the next interval.
c
      qprev = qright(lev)
      x0 = x(16)
      f0 = f(16)
      do 78 i = 1, 8
         f(2*i) = fsave(i,lev)
         x(2*i) = xsave(i,lev)
   78 continue
      go to 30
c
c   ***   stage 8   ***   finalize and return
c
   80 result = result + cor11
c
c   make sure errest not less than roundoff level.
c
      if (errest .eq. 0.0d0) return
   82 temp = dabs(result) + errest
      if (temp .ne. dabs(result)) return
      errest = 2.0d0*errest
      go to 82
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine cquanc8(fun,a,b,abserr,relerr,result,errest,nofun,flag)
c                        
      implicit real*8 (a-h,o-z)
      REAL*8  a, b, abserr, relerr, flag,errest
      complex*16 fun,result
      integer nofun
c

c
      REAL*8 w0,w1,w2,w3,w4,x0,stone,temp
      REAL*8 esterr,tolerr,step,x(16),xsave(8,30)
      complex*16 qprev,qnow,qdiff,qleft,area,cor11,f0
      complex*16 qright(31),f(16),fsave(8,30)
      integer levmin,levmax,levout,nomax,nofin,lev,nim,i,j
c
c   ***   stage 1 ***   general initialization
c   set constants.
c
      levmin = 1
      levmax = 30
      levout = 6
      nomax = 5000
      nofin = nomax - 8*(levmax-levout+2**(levout+1))
c
c   trouble when nofun reaches nofin
c
      w0 =   3956.0d0 / 14175.0d0
      w1 =  23552.0d0 / 14175.0d0
      w2 =  -3712.0d0 / 14175.0d0
      w3 =  41984.0d0 / 14175.0d0
      w4 = -18160.0d0 / 14175.0d0
c
c   initialize running sums to zero.
c
      flag = 0.0d0
      result = dcmplx(0.0d0,0.d0)
      cor11  = dcmplx(0.0d0,0.d0)
      errest = 0.0d0
      area   = dcmplx(0.0d0,0.d0)
      nofun = 0
      if (a .eq. b) return
c
c   ***   stage 2 ***   initialization for first interval
c
      lev = 0
      nim = 1
      x0 = a
      x(16) = b
      qprev  = dcmplx(0.0d0,0.d0)
      f0 = fun(x0)
      stone = (b - a) / 16.0d0
      x(8)  =  (x0  + x(16)) / 2.0d0
      x(4)  =  (x0  + x(8))  / 2.0d0
      x(12) =  (x(8)  + x(16)) / 2.0d0
      x(2)  =  (x0  + x(4))  / 2.0d0
      x(6)  =  (x(4)  + x(8))  / 2.0d0
      x(10) =  (x(8)  + x(12)) / 2.0d0
      x(14) =  (x(12) + x(16)) / 2.0d0
      do 25 j = 2, 16, 2   
            f(j) = fun(x(j))
   25 continue
      nofun = 9
c
c   ***   stage 3 ***   central calculation
c   requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16.
c   calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area.
c
   30 x(1) = (x0 + x(2)) / 2.0d0
      f(1) = fun(x(1))
      do 35 j = 3, 15, 2            
            x(j) = (x(j-1) + x(j+1)) / 2.0d0    
            f(j) = fun(x(j))
   35 continue
      nofun = nofun + 8
      step = (x(16) - x0) / 16.0d0
      qleft  =  (w0*(f0 + f(8))  + w1*(f(1)+f(7))  + w2*(f(2)+f(6))
     1  + w3*(f(3)+f(5))  +  w4*f(4)) * step
      qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14))
     1  + w3*(f(11)+f(13)) + w4*f(12)) * step
      qnow = qleft + qright(lev+1)
      qdiff = qnow - qprev
      area = area + qdiff
c
c   ***   stage 4 *** interval convergence test
c
      esterr = cdabs(qdiff) / 1023.0d0
      tolerr = dmax1(abserr,relerr*cdabs(area)) * (step/stone)
      if (lev .lt. levmin) go to 50
      if (lev .ge. levmax) go to 62
      if (nofun .gt. nofin) go to 60
      if (esterr .le. tolerr) go to 70
c
c   ***   stage 5   ***   no convergence
c   locate next interval.
c
   50 nim = 2*nim
      lev = lev+1
c
c   store right hand elements for future use.
c
      do 52 i = 1, 8  
            fsave(i,lev) = f(i+8)     
            xsave(i,lev) = x(i+8)
   52 continue
c
c   assemble left hand elements for immediate use.
c
      qprev = qleft
      do 55 i = 1, 8   
            j = -i           
            f(2*j+18) = f(j+9)     
            x(2*j+18) = x(j+9)
   55 continue
      go to 30
c
c   ***   stage 6   ***   trouble section
c   number of function values is about to exceed limit.
c
   60 nofin = 2*nofin
      levmax = levout
      flag = flag + (b - x0) / (b - a)
      go to 70
c
c   current level is levmax.
c
   62 flag = flag + 1.0d0
c
c   ***   stage 7   ***   interval converged
c   add contributions into running sums.
c
   70 result = result + qnow
      errest = errest + esterr
      cor11  = cor11  + qdiff / 1023.0d0
c
c   locate next interval.
c
   72 if (nim .eq. 2*(nim/2)) go to 75
      nim = nim/2
      lev = lev-1
      go to 72
   75 nim = nim + 1
      if (lev .le. 0) go to 80
c
c   assemble elements required for the next interval.
c
      qprev = qright(lev)
      x0 = x(16)
      f0 = f(16)
      do 78 i = 1, 8       
            f(2*i) = fsave(i,lev)     
            x(2*i) = xsave(i,lev)
   78 continue
      go to 30
c
c   ***   stage 8   ***   finalize and return
c
   80 result = result + cor11
c
c   make sure errest not less than roundoff level.
c
      if (errest .eq. 0.0d0) return
   82 temp = cdabs(result) + errest
      if (temp .ne. cdabs(result)) return
      errest = 2.0d0*errest
      go to 82
      end 

