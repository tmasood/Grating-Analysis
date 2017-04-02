CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z) 
      PRINT*, ' **************************************************** '
      PRINT*, '                                                      '    
      PRINT*, ' THE MAIN PROGRAM IS ONLY A CALLING PROGRAM, WHICH    '
      PRINT*, ' CONTAINS THREE MAJOR PARTS(MORE THAN 20 SUBROUTINES) '    
      PRINT*, '                                                      '
      PRINT*, ' FIRST SELECT MATERIAL PARAMETERS.                    '
      PRINT*, '                                                      '
      PRINT*, ' SECOND SELECT ENERGY LEVELS IN BOTH BANDS.           '
      PRINT*, '                                                      '
      PRINT*, ' THIRD FIND THE G(J), G(WAVELENGTH) AND RATE EQUATIONS' 
      PRINT*, '                                                      '
      PRINT*, ' **************************************************** '
      PRINT*, '                                                      '      
      PRINT*, ' MAKE YOUR SELECTION NOW!                             '
      PRINT*, '                                                      ' 
c    1 CALL SEL (N)
c      GOTO (10,20,30) N
c   10 CALL PARAB
c      GOTO 1
c   20 CALL BANDMIX
c      GOTO 1
c   30 PRINT*, ' THIS PROGRAM STOP HERE !'
c      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c      SUBROUTINE SEL (N)
c      IMPLICIT REAL*8 (A-H,O-Z)
c      WRITE(*,11)
c   11 FORMAT(/' ENTER 1 FOR PARABOLIC LASER MODEL'/
c     +        '       2 FOR BAND MIXING LASER MODEL'/
c     +        '       3 FOR EXIT'/)
c      READ(*,*) N
c      RETURN
c      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
c      SUBROUTINE PARAB
c      IMPLICIT REAL*8 (A-H,O-Z)     
    1 CALL SELECT (N)
      GO TO (10,20,30,40,50,60,70) N
   10 CALL PARAMETER     
      GOTO 1
   20 CALL CONDUCT
      GOTO 1
   30 CALL HEAVY
      GOTO 1
   40 CALL LIGHT
      GOTO 1
   50 CALL LASER
      GOTO 1 
   60 CALL RATEEQ
      GOTO 1
   70 PRINT*, ' BACK TO FIRST SELECTION PAGE !'
c      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c      SUBROUTINE BANDMIX
c      IMPLICIT REAL*8 (A-H,O-Z)
c    1 CALL SEL2(N)
c     GOTO (10,20,30,40,50) N
c      GOTO (10,20,50) N
c   10 CALL PARAMETER
c      GOTO 1
c   20 CALL CONDUCT
c      GOTO 1
c  30 CALL VALENCE 
c     GOTO 1
C   40 CALL LIGHT2
C      GOTO 1 
c   50 PRINT*, ' BACK TO FIRST SELECTION PROGRAM PAGE !'
c      RETURN
c      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    
      SUBROUTINE SELECT (N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NA  
      COMMON /PARA1/ NA
      WRITE(*,11)
   11 FORMAT(/' ENTER 1 FOR THE NECESSARY PARAMETERS'/
     +        '       2 FOR THE ENERGY VALUES OF CONDUCTION BAND'/
     +        '       3 FOR THE ENERGY VALUES OF HEAVY HOLE BAND'/
     +        '       4 FOR THE ENERGY VALUES OF LIGHT HOLE BAND'/
     +        '       5 FOR THE LASER G-J AND G(LAMBDA)'/
     +        '       6 FOR RATE EQUATIONS(TWO SECTION MODEL INCLUDED)'/
     +        '       7 FOR EXIT '/)
      READ(*,*) N
      NA=N
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c      SUBROUTINE SEL2 (N)
c      IMPLICIT REAL*8 (A-H,O-Z)
c      INTEGER NB  
c      COMMON /BANDMIX1/ NB
c      WRITE(*,11)
c   11 FORMAT(/' ENTER 1 FOR THE NECESSARY PARAMETERS'/
c     +        '       2 FOR THE ENERGY VALUES OF CONDUCTION BAND'/
C     +        '       3 FOR THE ENERGY VALUES OF HEAVY HOLE BAND'/
C     +        '       4 FOR THE ENERGY VALUES OF LIGHT HOLE BAND'/
C     +        '       5 FOR THE LASER G-J AND G(LAMBDA)'/
C     +        '       6 FOR RATE EQUATIONS(TWO SECTION MODEL INCLUDED)'/
c     +        '       5 FOR EXIT '/)
c      READ(*,*) N
c      NB=N
c      RETURN
c      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PARAMETER
      IMPLICIT REAL*8 (A-H,O-Z)
    1 CALL SELECT2 (NX)
      GOTO (10,20,30,40,50,60,70,80,90,100,110,120) NX
   10 CALL ALGAAS
      GOTO 1
   20 CALL INGAASP
      GOTO 1
   30 CALL INGAAS1
      GOTO 1
   40 CALL ALGAINAS
      GOTO 1
   50 CALL GAINP
      GOTO 1
   60 CALL INGAAS2
      GOTO 1
   70 CALL INGAAS3
      GOTO 1
   80 CALL ALINGAAS
      GOTO 1
   90 CALL INGAAS4
      GOTO 1
  100 CALL INGAALAS
      GOTO 1
  110 CALL INGAAS5
      GOTO 1
  120 PRINT*, ' THIS PROGRAM STOP HERE!, BACK TO MAIN PAGE'
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
      SUBROUTINE SELECT2 (NX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LX,LY,LZ
      INTEGER N,NX,IC      
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
      OPEN(1,FILE='cbandeg.dat',status='unknown')
      OPEN(4,FILE='vbandeg.dat',status='unknown') 
      WRITE(*,11)  
   11 FORMAT(/' ENTER 1 FOR AlGaAs/AlGaAs      '/
     +        '       2 FOR InGaAsP/InGaAs     '/
     +        '       3 FOR InGaAs/InGaAsP/InP '/
     +        '       4 FOR InGaAlAs/InGaAlAs  '/
     +        '       5 FOR GaInP/(AlGa)0.5In0.5P/AlInP '/
     +        '       6 FOR InGaAs/AlGaAs/AlGaAs '/
     +        '       7 FOR InGaAs/InGaAsP/Ga0.51In0.49P(MATCHED GaAs)'/
     +        '       8 FOR AlyInxGa1-x-yAs/AlzGa1-zAs/GaAs '/
     +        '       9 FOR InzGa1-zAs/AlxGayIn1-x-yAs/InP '/
     +        '      10 FOR InGaAlAs/InGaAlAs/AlAsxSb1-x(matched InP)'/
     +        '      11 FOR InzGa1-zAs/AlxGayIn1-x-yAs/AlAsxSb1-x '/
     +        '      12 FOR EXIT, BACK TO MAIN PAGE!      ')
      READ(*,*) NX
      IF (NX.EQ.12) THEN
      GOTO 1000
      ELSE
      ENDIF
C	CALL BANDIN
 1000	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	   
      SUBROUTINE BANDIN
	IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 LX,LY,LZ
      INTEGER N,NX,IC      
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
      OPEN(1,FILE='cbandeg.dat',status='unknown')
      OPEN(4,FILE='vbandeg.dat',status='unknown') 
      PRINT*, 'INPUT THE LAYER # FOR GRIN STRUCTURE(STEP)'
      PRINT*, 'STEP N='
      READ(*,*) N                           
      PRINT*, ' INPUT THE WELL WAVELENGTH (um)'
      READ(*,*) LX
      PRINT*, ' INPUT THE BARRIER WAVELENGTH (um)'
      READ(*,*) LY
      PRINT*, ' INPUT THE CLADDING WAVELENGTH (um)'
      READ(*,*) LZ
      BANDEG=1.24D0/LX
      EGB=1.24D0/LY
      EGC=1.24D0/LZ     
      A=6.0584D0  ! InAs
      B=5.8688D0  ! InP
      C=5.6533D0  ! GaAs
      D=5.4512D0  ! GaP
      E=5.6611D0  ! AlAs                          
      F=5.4512D0  ! AlP
      G=6.136D0   ! AlSb
      H=6.096D0   ! GaSb
      HI=6.479D0  ! InSb
      PRINT*, ' BANDGAP ENERGY OF QUANTUM WELL= ',BANDEG,' eV'
      PRINT*, ' INPUT CLADDING, BARRIER,QUANTUM WELL WIDTH (A)'
      READ(*,*) CLAY,BLAY,QLAY       
      BBLAY=BLAY/(N-1)
 1000 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE BANDOUT
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LAYERX(200),LAYERY(200),DELTAC(200),DELTAV(200),VB(200)
      REAL*8 VH(200)
      INTEGER IC,N
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
      REWIND 1
      REWIND 4
      PRINT*, ' WRITE CONDUCTION BAND PARAMETERS INTO CBANDEG.DAT'
      WRITE(1,10) STRAIN,AZ1*1.D-10   
      WRITE(1,11) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)             
      DO I=N-1,1,-1 
      WRITE(1,11) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I)
      ENDDO
      WRITE(1,11) QLAY,LAYERX(N+1),LAYERY(N+1),VB(IC)
      DO I=1,N-1
      WRITE(1,11) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I)
      ENDDO
      WRITE(1,11) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)
C      ELSEIF (J.EQ.2) THEN
      PRINT*, '                                                  '
      PRINT*, ' WRITE VALENCE BAND PARAMETERS INTO VBANDEG.DAT'
      WRITE(4,10) STRAIN,AZ1*1.D-10
      WRITE(4,11) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)             
      DO I=N-1,1,-1 
      WRITE(4,11) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)
      ENDDO
      WRITE(4,11) QLAY,LAYERX(N+1),LAYERY(N+1),VH(IC)
      DO I=1,N-1
      WRITE(4,11) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)
      ENDDO
      WRITE(4,11) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)
C      ELSE
C      ENDIF
   10 FORMAT(2X,E12.6,2X,E12.6)
   11 FORMAT(2(2X,E15.8),2(2X,F12.7))
      RETURN
      END     
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C             
      SUBROUTINE ALGAAS
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 LAYERX(200),DELTAC(200),VB(200),VH(200)
      REAL*8 DELTAV(200),EG(200),INCRE,LAYERY(200)  
      INTEGER N,IC          
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1	CALL BANDIN
      A1=1.147D0
      B1=0.2147D0
      C1=(1.6562675D0-BANDEG)
      C2=(1.6562675D0-EGB)
      C3=(1.6562675D0-EGC)               
      IC=N+1
      IF ((BANDEG-1.98515D0).LE.0.0000D0) THEN       
      G1=BANDEG-1.424D0                      
      LAYERX(N+1)=G1/1.247D0
      ELSE      
      LAYERX(N+1)=(-B1+DSQRT((B1**2)-4.D0*A1*C1))/(2.D0*A1) 
      ENDIF                   
C      
      IF ((EGB-1.98515D0).LE.0.0000D0) THEN
      LAYERX(1)=(EGB-1.424D0)/(1.247D0)
      ELSE
      LAYERX(1)=(-B1+DSQRT((B1**2)-4.D0*A1*C2))/(2.D0*A1)
      ENDIF  
C                        
      IF ((EGC-1.98515D0).LE.0.0000D0) THEN
      LAYERX(N)=(EGC-1.424D0)/(1.247D0)
      ELSE
      LAYERX(N)=(-B1+DSQRT((B1**2)-4.D0*A1*C3))/(2.D0*A1)
      ENDIF
C
      AZ1=LAYERX(N+1)*E+(1-LAYERX(N+1))*C
C      PRINT*, 'INPUT EX'
C      READ(*,*) EX
   10 FORMAT(2X,E12.6,2X,E12.6)
   11 FORMAT(2(2X,E15.8),2(2X,F12.7)) 
      INCRE=(LAYERX(N)-LAYERX(1))/(N-1)
      DO I=2,N-1
      LAYERX(I)=LAYERX(I-1)+INCRE
      ENDDO
      DO I=1,N
      IF (LAYERX(I).LE.0.45) THEN
      EG(I)=1.424+1.247*LAYERX(I)
      ELSE
      EG(I)=1.424+1.247*LAYERX(I)+1.147*((LAYERX(I)-0.45)**2)
      ENDIF
      DELTAC(I)=(EG(I)-BANDEG)*0.65
      DELTAV(I)=-((EG(I)-BANDEG)*0.35)
      ENDDO            
	STRAIN=0.D0
      PRINT*, '                                                  '
      CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE INGAASP
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VH(200),LAYERX(200),VB(200),DELTAC(200),INCREX,INCREY
      REAL*8 DELTAV(200),LAYERY(200),EG(200),X(200),Y(200)
      INTEGER N,IC
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1	CALL BANDIN
      PRINT*, '                                                      '
      PRINT*, ' For InGaAsP, only compress strain (~1.5%) available'
      PRINT*, '                                                      '
      PRINT*, 'INPUT EX'
      READ(*,*) EX
C      IF (EX.EQ.0.00) THEN
C	  GOTO 31
C      ELSEIF (EX.NE.0.00000) THEN
C	  GOTO 30
C      ELSE
C      ENDIF
C   30 IF ((EX.LT.0.00).AND.(EX.GT.-0.005)) THEN
C      LAYERY(N+1)=((0.72-DSQRT((0.72)**2-4*(0.12)*(1.35-BANDEG)))/0.24)
C     +            -(DABS(EX)*50*0.732)
C      ELSEIF ((EX.LE.-0.005).AND.(EX.GE.-0.015)) THEN
C      LAYERY(N+1)=((0.72-DSQRT((0.72)**2-4*(0.12)*(1.35-BANDEG)))/0.24)
C     +            -(0.183+(DABS(EX)-0.005)*15.6)
   30 IF (EX.LT.0.000) THEN
      XLAYERY=(0.72-DSQRT((0.72)**2-4*(0.12)*(1.35-BANDEG)))/0.24
	DX=EX/(XLAYERY*0.069026)
	LAYERX(N+1)=(0.1894*XLAYERY/(0.4184-0.013*XLAYERY))+DX
	XDX=LAYERX(N+1)
      LAYERY(N+1)=(5.8688*(1-EX)-5.8688+0.4176*XDX)/(0.1896+0.0125*XDX)
      ELSEIF (EX.EQ.0.00) THEN
      LAYERY(N+1)=(0.72-DSQRT((0.72)**2-4*(0.12)*(1.35-BANDEG)))/0.24
	ELSEIF (EX.GT.0.0000) THEN
      XLAYERY=(0.72-DSQRT((0.72)**2-4*(0.12)*(1.35-BANDEG)))/0.24
	DX=EX/(XLAYERY*0.069026)
	LAYERX(N+1)=(0.1894*XLAYERY/(0.4184-0.013*XLAYERY))+DX
	XDX=LAYERX(N+1)
      LAYERY(N+1)=(5.8688*(1-EX)-5.8688+0.4176*XDX)/(0.1896+0.0125*XDX)
	ELSE
      ENDIF
   31 LAYERY(1)=(0.72-DSQRT((0.72)**2-4*(0.12)*(1.35-EGB)))/0.24
      LAYERX(1)=0.1894*LAYERY(1)/(0.4184-0.013*LAYERY(1))
      LAYERY(N)=(0.72-DSQRT((0.72)**2-4*(0.12)*(1.35-EGC)))/0.24
      LAYERX(N)=0.1894*LAYERY(N)/(0.4184-0.013*LAYERY(N))
      IC=N+1
      X1=(1-LAYERX(1))*LAYERY(1)*A
      X2=(1-LAYERX(1))*(1-LAYERY(1))*B
      X3=LAYERX(1)*LAYERY(1)*C
      X4=LAYERX(1)*(1-LAYERY(1))*D
      AZB2=X1+X2+X3+X4
C      LAYERX(N+1)=0.1894*LAYERY(N+1)/(0.4184-0.013*LAYERY(N+1))
	 AZ2=AZB2
C	 IF (EX.EQ.0.00) THEN
C	 GOTO 33
C	 ELSE
C	 GOTO 32
C	 ENDIF
   32 AA=LAYERY(N+1)*A
      BB=(1-LAYERY(N+1))*B
      CC=LAYERY(N+1)*C
      DD=(1-LAYERY(N+1))*D
      LAYERX(N+1)=(AA+BB+AZB2*(EX-1))/(AA+BB-CC-DD)
      XX1=(1-LAYERX(N+1))*LAYERY(N+1)*A
      XX2=(1-LAYERX(N+1))*(1-LAYERY(N+1))*B
      XX3=LAYERX(N+1)*LAYERY(N+1)*C
      XX4=LAYERX(N+1)*(1-LAYERY(N+1))*D
      AZ1=XX1+XX2+XX3+XX4  
   33 INCREX=(LAYERX(N)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N)-LAYERY(1))/(N-1)
      DO I=2,N-1
      LAYERY(I)=LAYERY(I-1)+INCREY
      LAYERX(I)=LAYERX(I-1)+INCREX
      ENDDO
      X(IC)=LAYERX(N+1)
      Y(IC)=LAYERY(N+1)
      C11=(1.-X(IC))*Y(IC)*8.329+(1.-X(IC))*(1.-Y(IC))*10.22
     +    +X(IC)*Y(IC)*11.88+X(IC)*(1.-Y(IC))*14.12
      C12=(1.-X(IC))*Y(IC)*4.526+(1.-X(IC))*(1.-Y(IC))*5.76
     +    +X(IC)*Y(IC)*5.38+X(IC)*(1.-Y(IC))*6.253
      CAA=(C11-C12)/C11
      AX=(1.-X(IC))*Y(IC)*(-6.08)+(1.-X(IC))*(1.-Y(IC))*(-6.31)
     +   +X(IC)*Y(IC)*(-9.77)+X(IC)*(1.-Y(IC))*(-8.83)
      VB(IC)=(2.D0*(2.D0/3.D0)*AX*CAA*EX)
      VH(IC)=(-2.D0*(1.D0/3.D0)*AX*CAA*EX)
      EGN=BANDEG+VB(IC)-VH(IC)
   10 FORMAT(2X,E12.6,2X,E12.6)   
   11 FORMAT(2(2X,E15.8),2(2X,F12.7)) 
   12 FORMAT(2(2X,E14.7),3(2X,F10.5))
      DO I=1,N
      EG(I)=1.35-0.72*LAYERY(I)+0.12*(LAYERY(I)**2)
      DELTAC(I)=(EG(I)-BANDEG)*0.39
      DELTAV(I)=-((EG(I)-BANDEG)*0.61)
      ENDDO
c     EG(IC)=0.012*LAYERY(N+1)**2-0.6228*LAYERY(N+1)+1.35
c     EG(IC)=0.012*LAYERY(N+1)**2-(0.7272-EGB)*LAYERY(N+1)+EGB
      EG(IC)=BANDEG
      PRINT*, '                                                  '
      AZB2=X1+X2+X3+X4
	STRAIN=(AZB2-AZ1)/AZB2
      PRINT*, ' STRAIN=',(AZB2-AZ1)/AZB2
      CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                    
      SUBROUTINE INGAAS1
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VH(200),LAYERX(200),VB(200),DELTAC(200)
      REAL*8 DELTAV(200),LAYERY(200),EG(200),X(200),INCREX,INCREY
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1	CALL BANDIN 
      G2=0.324D0-BANDEG
      LAYERX(N+1)=(-0.7+DSQRT((.7**2)-(4*0.4*G2)))/0.8             
      LAYERY(1)=(0.775-DSQRT((0.775)**2-4*(0.12)*(1.35-EGB)))/0.298
      LAYERX(1)=0.1894*LAYERY(1)/(0.4184-0.013*LAYERY(1))      
      AZ1=6.0584-0.4051*LAYERX(N+1)
      STRAIN=(B-AZ1)/B
      PRINT*, ' STRAIN FOR In1-xGaxAs= ',STRAIN
C      PRINT*, 'INPUT EX'
C      READ(*,*) EX             
      INCREX=(LAYERX(N)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N)-LAYERY(1))/(N-1)
      DO I=2,N-1
      LAYERY(I)=LAYERY(I-1)+INCREY
      LAYERX(I)=LAYERX(I-1)+INCREX
      ENDDO                                  
      IC=(N+1)
      X(IC)=LAYERX(N+1)
      C11=X(IC)*11.88+(1.D0-X(IC))*8.329D0
      C12=X(IC)*5.38D0+(1.D0-X(IC))*4.526D0
      CAA=(C11-C12)/C11
      AA=X(IC)*(-9.77D0)+(1.D0-X(IC))*(-6.D0)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*STRAIN) 
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*STRAIN)
      EGN=BANDEG+VB(IC)-VH(IC)
C      PRINT*, ' NEW BAND GAP=', EGN
      DO I=1,N
      EG(I)=1.35-0.775*LAYERY(I)+0.149*(LAYERY(I)**2)
      DELTAC(I)=(EG(I)-BANDEG)*0.36
      DELTAV(I)=-((EG(I)-BANDEG)*0.64)
      ENDDO   
      PRINT*, '                                                  '
      CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                      
      SUBROUTINE ALGAINAS
      PARAMETER (II=200)
      IMPLICIT REAL*8 (A-H,O-Z)              
      REAL*8 VH(500),LAYERX(500),VB(500),DELTAC(500)
      REAL*8 DELTAV(500),LAYERY(500),X(500),Y(500),INCREX,INCREY
      REAL*8 CT(500),CG(500),CA(500),T(500),TG(500),TC(500)
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
    1	CALL BANDIN         
      PRINT*, 'INPUT I=1 FOR COMPRESS, 2 FOR TENSILE'
      PRINT*, 'I=#'   
      REWIND 1
      REWIND 4 
      IC=(N+1)
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1002
      ELSEIF (I.EQ.2) THEN
      GOTO 1003
      ELSE
      ENDIF
 1002 DO J=1,II
      CT(J)=0.53D0+((0.820D0-0.53D0)/II)*J
      CG(J)=0.75D0-((0.75D0-0.526D0)/II)*J
      CA(J)=1.548D0-((1.548D0-1.516D0)/II)*J
      YY=J*0.01
      WRITE(*,101) CT(J),CG(J),CA(J),YY,J      
  101 FORMAT(2X,'CT',F12.8,2X,'CG',F12.8,2X,'CA',F12.8,'COMPS',F6.3,'%',
     +       2X,I3)
      ENDDO
      PRINT*, 'CALCULATE THE MATCHED BARRIER AND CLADDING COMPOSITION'
      PRINT*, 'X,Z--FOR BARRIER,XC,ZC--FOR CLADDING'
      Z1=(EGB-0.75D0)/1.548D0
      X1=1-0.53D0-Z1
      ZC=(EGC-0.75D0)/1.548D0
      XC=1-0.53D0-ZC
      IF (XC.LT.0.00000D0) THEN
      XC=0.D0
      ELSE
      ENDIF
      WRITE(*,1004) Z1,X1,ZC,XC
 1004 FORMAT(1X,'Z=',F12.7,1X,'X=',F12.7,1X,'ZC=',F12.7,1X,'XC=',F12.7)      
C      PRINT*, 'INPUT STRAIN "+" FOR COMPRESS'
      PRINT*, 'AND CALCULATE THE WELL COMPOSITION XW AND ZW'
      PRINT*, 'J= #'
      READ(*,*) J
      EX=-(J*0.0001)
      ZW=(BANDEG-CG(J))/CA(J)
      XW=1-CT(J)-ZW
      LAYERX(N+1)=XW
      LAYERY(N+1)=ZW
      LAYERX(N)=XC
      LAYERY(N)=ZC
      LAYERX(1)=X1
      LAYERY(1)=Z1      
      X(IC)=LAYERX(N+1)
      BBLAY=BLAY/(N-1)
      INCREX=(LAYERX(N)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N)-LAYERY(1))/(N-1)
      DO I=2,N-1
      LAYERY(I)=LAYERY(I-1)+INCREY
      LAYERX(I)=LAYERX(I-1)+INCREX
      ENDDO            
      C11=(1.-X(IC)-Y(IC))*8.329+Y(IC)*12.02+X(IC)*11.88
      C12=(1.-X(IC)-Y(IC))*4.526+Y(IC)*5.7+X(IC)*5.38
      CAA=(C11-C12)/C11
      AA=(1.-X(IC)-Y(IC))*(-6.08)+Y(IC)*(-8.11)+X(IC)*(-9.77)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*EX)  
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*EX)
      EGN=BANDEG+VB(IC)-VH(IC)
C     PRINT*, ' NEW BAND GAP=', EGN
      DO I=1,N
      EC1=(EGC-BANDEG)*0.72
      EC2=(EGB-BANDEG)*0.72
      DEG=(EC1-EC2)/(N-1)
      DELTAC(I)=EC2+(I-1)*DEG
      EV1=(EGC-BANDEG)*0.28
      EV2=(EGB-BANDEG)*0.28
      DVG=(EV1-EV2)/(N-1)
      DELTAV(I)=-(EV2+(I-1)*DVG)
      ENDDO                   
      AZ3=(1-XW-ZW)*A+ZW*E+XW*C
      PRINT*, 'INPUT STRAIN'                     
      READ(*,*) EX
      PRINT*, '                                                  '
      PRINT*, ' WRITE CONDUCTION BAND PARAMETERS INTO CBANDEG.DAT'
      WRITE(1,10) EX,AZ3*1.D-10
      WRITE(1,12) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)
      DO I=N-1,1,-1
      WRITE(1,12) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I) 
      ENDDO
      WRITE(1,12) QLAY,LAYERX(N+1),LAYERY(N+1),VB(IC)
      DO I=1,N-1
      WRITE(1,12) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I)
      ENDDO
      WRITE(1,12) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)
      PRINT*, '                                                  '
      PRINT*, ' WRITE VALENCE BAND PARAMETERS INTO VBANDEG.DAT'      
      WRITE(4,10) EX,AZ3*1.D-10
      WRITE(4,12) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)
      DO I=N-1,1,-1
      WRITE(4,12) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)  
      ENDDO
      WRITE(4,12) QLAY,LAYERX(N+1),LAYERY(N+1),VH(IC) 
      DO I=1,N-1
      WRITE(4,12) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)   
      ENDDO
      WRITE(4,12) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)
	CALL BANDOUT   
      GOTO 1100 
C
C     FOR TENSILE STRAIN COMPOSITION  !
C 
 1003 DO J=1,II
      T(J)=0.53D0-((0.53D0-0.230D0)/II)*J
      TG(J)=0.75D0+((0.83D0-0.75D0)/II)*J
      TC(J)=1.548D0+((1.588D0-1.548D0)/II)*J
      XX=J*0.01
      WRITE(*,102) T(J),TG(J),TC(J),XX,J
  102 FORMAT(2X,'T',F12.8,2X,'TG',F12.8,2X,'TC',F12.8,'TENSS',F6.3,'%',
     +       2X,I3)
      ENDDO
      PRINT*, 'CALCULATE THE MATCHED BARRIER AND CLADDING COMPOSITION'
      PRINT*, 'X,Z--FOR BARRIER,XC,ZC--FOR CLADDING'
      Z1=(EGB-0.75D0)/1.548D0
      X1=1-0.53D0-Z1
      ZC=(EGC-0.75D0)/1.548D0
      XC=1-0.53D0-ZC
      IF (XC.LT.0.00000D0) THEN
      XC=0.D0
      ELSE
      ENDIF
      WRITE(*,1004) Z1,X1,ZC,XC 
C      PRINT*, 'INPUT STRAIN "-" FOR TENSILE'
      PRINT*, 'AND CALCULATE THE WELL COMPOSITION XW AND ZW'
      PRINT*, 'J= #'
      READ(*,*) J
      EX=J*0.0001
      ZW=(BANDEG-TG(J))/TC(J)
      XW=1-T(J)-ZW
      LAYERX(N+1)=XW
      LAYERY(N+1)=ZW
      LAYERX(N)=XC
      LAYERY(N)=ZC
      LAYERX(1)=X1
      LAYERY(1)=Z1
      X(IC)=LAYERX(N+1)      
      Y(IC)=LAYERY(N+1)
      BBLAY=BLAY/(N-1)
      INCREX=(LAYERX(N)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N)-LAYERY(1))/(N-1)
      DO I=2,N-1
      LAYERY(I)=LAYERY(I-1)+INCREY
      LAYERX(I)=LAYERX(I-1)+INCREX
      ENDDO
      C11=(1.-X(IC)-Y(IC))*8.329+Y(IC)*12.02+X(IC)*11.88
      C12=(1.-X(IC)-Y(IC))*4.526+Y(IC)*5.7+X(IC)*5.38
      CAA=(C11-C12)/C11
      AA=(1.-X(IC)-Y(IC))*(-6.08)+Y(IC)*(-8.11)+X(IC)*(-9.77)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*EX)
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*EX)
      EGN=BANDEG+VB(IC)-VH(IC)
C      PRINT*, ' NEW BAND GAP=', EGN
      DO I=1,N
      EC1=(EGC-BANDEG)*0.72
      EC2=(EGB-BANDEG)*0.72
      DEG=(EC1-EC2)/(N-1)
      DELTAC(I)=EC2+(I-1)*DEG
      EV1=(EGC-BANDEG)*0.28
      EV2=(EGB-BANDEG)*0.28
      DVG=(EV1-EV2)/(N-1)
      DELTAV(I)=-(EV2+(I-1)*DVG)
      ENDDO                   
      AZ3=(1-XW-ZW)*A+ZW*E+XW*C
      PRINT*, 'INPUT STRAIN'                     
      READ(*,*) EX
      PRINT*, '                                                  '
      PRINT*, ' WRITE CONDUCTION BAND PARAMETERS INTO CBANDEG.DAT'
      WRITE(1,10) EX,AZ3*1.D-10
      WRITE(1,12) CLAY,LAYERX(N),LAYERY(N),DELTAC(N) 
      DO I=N-1,1,-1
      WRITE(1,12) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I) 
      ENDDO
      WRITE(1,12) QLAY,LAYERX(N+1),LAYERY(N+1),VB(IC)   
      DO I=1,N-1
      WRITE(1,12) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I)   
      ENDDO
      WRITE(1,12) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)
      PRINT*, '                                                  '
      PRINT*, ' WRITE VALENCE BAND PARAMETERS INTO VBANDEG.DAT'  
      WRITE(4,10) EX,AZ3*1.D-10
      WRITE(4,12) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)  
      DO I=N-1,1,-1
      WRITE(4,12) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I) 
      ENDDO
      WRITE(4,12) QLAY,LAYERX(N+1),LAYERY(N+1),VH(IC)
      DO I=1,N-1
      WRITE(4,12) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)
      ENDDO
      WRITE(4,12) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)  
   10 FORMAT(2X,E12.6,2X,E12.6)   
   12 FORMAT(2(2X,E14.7),3(2X,F10.5))
      CALL BANDOUT           
 1100 PRINT*, ' INPUT 1 FOR NEW CALCULATION'
      PRINT*, '       2 FOR EXIT'
      PRINT*, ' INPUT =?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE GAINP
      IMPLICIT REAL*8 (A-H,O-Z)             
      REAL*8 VH(200),LAYERX(200),VB(200),DELTAC(200),LAYERZ(200)
      REAL*8 DELTAV(200),LAYERY(200),EG(200),X(200),LAYERW(200)
      REAL*8 INCREX,INCREY
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC  
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1	CALL BANDIN 
      IC=(N+1)
c       LAYERX(N)=(EGC-1.351D0)/2.23D0  ! FOR AlInP (CLADDING)
      LAYERX(N)=(EGC-1.91)/0.61
      LAYERY(N)=1.D0-LAYERX(N)  ! FOR AlGaInP (CLADDING)
      LAYERX(1)=(EGB-1.91)/0.61
      LAYERZ(1)=LAYERX(1)*0.5
      LAYERY(1)=1.D0-LAYERX(1)
      LAYERW(1)=LAYERY(1)*0.5
      IF (DABS(BANDEG-1.91D0).LE.0.0001) THEN
          LAYERX(N+1)=0.5D0
      ELSE
          G3=1.351D0-BANDEG        
          LAYERX(N+1)=(-0.643+DSQRT((0.643**2)-4*0.786*G3))/(2*0.786)
      ENDIF
      AZ1=5.8686D0-0.4174*LAYERX(N+1)
      G4=1-LAYERW(1)-LAYERZ(1)
      AZB5=5.6533D0
      STRAIN=(AZB5-AZ1)/AZB5
      PRINT*, ' STRAIN FOR GaInP ',STRAIN           
      INCREX=(LAYERX(N)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N)-LAYERY(1))/(N-1)
      DO I=2,N-1
      LAYERY(I)=LAYERY(I-1)+INCREY
      LAYERX(I)=LAYERX(I-1)+INCREX
      ENDDO        
      X(IC)=LAYERX(N+1)
      C11=10.22+3.9*X(IC)
      C12=5.76+0.49*X(IC)
      CAA=(C11-C12)/C11
      AA=X(IC)*(-8.83)+(1-X(IC))*(-6.31)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*STRAIN)
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*STRAIN)
      EGN=BANDEG+VB(IC)-VH(IC)
      DO I=1,N
      EG(I)=1.91D0+0.61*LAYERX(I)
       DELTAC(I)=(EG(I)-BANDEG)*0.35
       DELTAV(I)=-((EG(I)-BANDEG)*0.65)
      ENDDO
      PRINT*, '                                                  '
      CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                
      SUBROUTINE INGAAS2
      IMPLICIT REAL*8 (A-H,O-Z)               
      REAL*8 VH(200),LAYERX(200),VB(200),DELTAC(200)
      REAL*8 DELTAV(200),EG(200),X(200),INCRE,LAYERY(200)
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1	CALL BANDIN  
      IC=(N+1)
      GX=1.424D0-BANDEG
      A1=1.147D0
      B1=0.2147D0
      C1=(1.6562675D0-BANDEG)
      C2=(1.6562675D0-EGB)
      C3=(1.6562675D0-EGC)
      LAYERX(N+1)=(1.619-DSQRT((1.619**2)-4.0*0.555*GX))/(2.0*0.555) 
      IF ((EGB-1.98515D0).LE.0.0000D0.AND.(EGB-1.424D0).GT.0.005D0) THEN
      LAYERX(1)=(EGB-1.424D0)/(1.247D0)           ! Al
      ELSEIF ((EGB-1.98515D0).LE.0.0000D0.AND.
     +       (EGB-1.424D0).LE.0.005D0) THEN
      EGB=1.424D0
      LAYERX(1)=0.0000D0
      ELSE
      LAYERX(1)=(-B1+DSQRT((B1**2)-4.D0*A1*C2))/(2.D0*A1)  ! Al
      ENDIF
C
      IF ((EGC-1.98515D0).LE.0.0000D0) THEN
      LAYERX(N)=(EGC-1.424D0)/(1.247D0)           ! Al
      ELSE
      LAYERX(N)=(-B1+DSQRT((B1**2)-4.D0*A1*C3))/(2.D0*A1)  ! Al
      ENDIF
      AZ1=5.6533D0+0.405D0*LAYERX(N+1)
      AZB6=5.6533D0+0.0078D0*LAYERX(1)
      STRAIN=(AZB6-AZ1)/AZB6
      PRINT*, ' STRAIN FOR InGaAs/AlGaAs IS ', STRAIN
      INCRE=(LAYERX(N)-LAYERX(1))/(N-1)
      DO I=2,N-1
      LAYERX(I)=LAYERX(I-1)+INCRE   ! Al
      ENDDO                   
      X(IC)=LAYERX(N+1)
      C11=X(IC)*8.33+(1-X(IC))*11.88
      C12=X(IC)*4.53+(1-X(IC))*5.38
      CAA=(C11-C12)/C11
      AA=X(IC)*(-6)+(1-X(IC))*(-9.77)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*STRAIN)
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*STRAIN) 
      EGN=BANDEG+VB(IC)-VH(IC)
C      PRINT*, ' NEW BAND GAP=', EGN
      DO I=1,N
      IF (LAYERX(I).LE.0.45) THEN
      EG(I)=1.424+1.247*LAYERX(I)
      ELSE
      EG(I)=1.424+1.247*LAYERX(I)+1.147*((LAYERX(I)-0.45)**2)
      ENDIF    
      DELTAC(I)=(EG(I)-BANDEG)*0.6
      DELTAV(I)=-((EG(I)-BANDEG)*0.4)
      ENDDO
      PRINT*, '                                                  '
      CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE INGAAS3
      IMPLICIT REAL*8 (A-H,O-Z)             
      REAL*8 VH(200),LAYERX(200),VB(200),DELTAC(200)
      REAL*8 DELTAV(200),LAYERY(200),EG(200),X(200),Y(200)  
      REAL*8 CEL(200),CEX(200),CEY(200),CEG(200),INCREX,INCREY
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC  
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1	CALL BANDIN 
      PRINT*, 'CALCULATE THE InGaAs/InGaAsP/GaInP(MATCHED GaAs)'
      PRINT*, 'FOR GaInP MATCHED TO GaAs LZ=0.6583947'
      G7=1.351-EGC
      GX7=0.36-BANDEG
      IC=N+1
      LAYERX(N+1)=(-0.509+DSQRT((0.509**2)-4*0.555*GX7))/(2*0.555)
      LAYERX(N)=(-0.643+DSQRT((0.643**2)-4*(0.786*G7)))/(2*0.786)
      LAYERY(N)=0.D0
      LAYERY(N+1)=1.D0 
      X(IC)=LAYERX(N+1)
      Y(IC)=LAYERY(N+1)
      EG(N)=1.351+0.643*LAYERX(N)+0.786*(LAYERX(N)**2)
      AZ1=5.6533+0.405*(1-LAYERX(N+1))
      CALL SPLINE (CEL,CEY,CEX)
      DO I=1,100
      WRITE(*,3002) CEL(I),CEY(I),CEX(I),I
      CEG(I)=1.24/CEL(I)
 3002 FORMAT(2X,'CEL=',F10.7,2X,'CEY=',F10.7,2X,'CEX=',F10.7,'  I=',I4)
      ENDDO
      PRINT*, 'INPUT WAVELENGTH FOR BARRIER I = '
      READ(*,*) IX
      BBLAY=BLAY/(N-1)
      LAYERX(1)=CEX(IX)
      LAYERY(1)=CEY(IX)
      EG(1)=CEG(IX)
      CXX=CEX(IX)
      CYY=CEY(IX)
      AZB7=(1-CXX)*CYY*A+(1-CXX)*(1-CYY)*B+CXX*CYY*C+CXX*(1-CYY)*D
      PRINT*, 'LAYERX(1)=',LAYERX(1),' LAYERY(1)=',LAYERY(1)
      STRAIN=(AZB7-AZ1)/AZB7
      PRINT*, ' STRAIN FOR In1-xGaxAs= ',STRAIN
C      PRINT*, 'INPUT EX'
C      READ(*,*) EX
      INCREX=(LAYERX(N)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N)-LAYERY(1))/(N-1)
      INCREG=(EG(N)-EG(1))/(N-1)
      DO I=2,N-1
      LAYERY(I)=LAYERY(I-1)+INCREY
      LAYERX(I)=LAYERX(I-1)+INCREX
      EG(I)=EG(I-1)+INCREG
      ENDDO            
      C11=X(IC)*11.88+(1-X(IC))*8.33
      C12=X(IC)*5.38+(1-X(IC))*4.53
      CAA=(C11-C12)/C11
      AA=X(IC)*(-9.77)+(1-X(IC))*(-6)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*STRAIN)
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*STRAIN)
      EGN=BANDEG+VB(IC)-VH(IC)
C      PRINT*, ' NEW BAND GAP=', EGN
      DO I=1,N
      DELTAC(I)=(EG(I)-BANDEG)*0.68
      DELTAV(I)=-((EG(I)-BANDEG)*0.32)
      ENDDO
      PRINT*, '                                                  '
      CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END       
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      SUBROUTINE ALINGAAS
      IMPLICIT REAL*8 (A-H,O-Z)             
      REAL*8 VH(200),LAYERX(200),VB(200),DELTAC(200)
      REAL*8 DELTAV(200),LAYERY(200),EG(200),X(200),Y(200)  
      REAL*8 INCRE
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC    
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1 CALL BANDIN
      PRINT*, ' CALCULATE AlyInxGa1-x-yAs/AlzGa1-zAs/GaAs'
      PRINT*, ' INPUT Y AND X FOR QUANTUM WELL'
      PRINT*, ' Y= , X= '
      READ(*,*) Y1,X1
      BANDEG=1.424+1.455*Y1+0.191*(Y1*Y1)-1.614*X1+0.55*(X1*X1)
     +       +0.043*X1*Y1
      PRINT*, 'BANDGAP=',BANDEG  
      BBLAY=BLAY/(N-1)
      LAYERX(N+1)=X1
      LAYERY(N+1)=Y1
      IC=N+1
      X(IC)=X1
      Y(IC)=Y1
      A1=1.147D0
      B1=0.2147D0
      C1=(1.6562675D0-BANDEG)
      C2=(1.6562675D0-EGB)
      C3=(1.6562675D0-EGC)
      IF ((EGB-1.98515D0).LE.0.0000D0) THEN
      LAYERX(1)=(EGB-1.424D0)/(1.247D0)
      ELSE
      LAYERX(1)=(-B1+DSQRT((B1**2)-4.D0*A1*C2))/(2.D0*A1)
      ENDIF
C
      IF ((EGC-1.98515D0).LE.0.0000D0) THEN
      LAYERX(N)=(EGC-1.424D0)/(1.247D0)
      ELSE
      LAYERX(N)=(-B1+DSQRT((B1**2)-4.D0*A1*C3))/(2.D0*A1)
      ENDIF
C
      AZ1=Y(IC)*E+X(IC)*A+(1-X(IC)-Y(IC))*C
      AZB8=LAYERX(1)*E+(1-LAYERX(1))*C
      PRINT*, 'WELL LATTICE = ', AZ1 ,' BARRIER LATTICE = ',AZB8
      STRAIN=(AZB8-AZ1)/AZB8
      PRINT*, ' STRAIN = ',STRAIN
C      PRINT*, ' INPUT EX'
C      READ(*,*) EX
      INCRE=(LAYERX(N)-LAYERX(1))/(N-1)
      DO I=2,N-1
	 LAYERX(I)=LAYERX(I-1)+INCRE
      ENDDO
      DO I=1,N
	 IF (LAYERX(I).LE.0.45) THEN
	    EG(I)=1.424+1.247*LAYERX(I)
	    ELSE
	    EG(I)=1.424+1.247*LAYERX(I)+1.147*((LAYERX(I)-0.45)**2)
	 ENDIF                             
	 C11=Y(IC)*12.02+X(IC)*8.33+(1-X(IC)-Y(IC))*11.88
	 C12=Y(IC)*5.7+X(IC)*4.53+(1-X(IC)-Y(IC))*5.38
	 CAA=(C11-C12)/C11
	 AA=Y(IC)*(-8.11)+X(IC)*(-6.08)+(1-X(IC)-Y(IC))*(-9.77)
	 VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*STRAIN)
	 VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*STRAIN)
	 EGN=BANDEG+VB(IC)-VH(IC) 
	 DELTAC(I)=(EG(I)-BANDEG)*0.68
	 DELTAV(I)=-((EG(I)-BANDEG)*0.32)
      ENDDO
       PRINT*, '                                                  '
       CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE INGAAS4
      IMPLICIT REAL*8 (A-H,O-Z)             
      REAL*8 VH(200),LAYERX(200),VB(200),DELTAC(200)
      REAL*8 DELTAV(200),LAYERY(200),X(200),INCREX,INCREY  
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC 
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1	CALL BANDIN                                   
      PRINT*, 'THE In(z)Ga(1-z)As/Al(x)Ga(y)In(1-x-y)As/InP MATERIAL'  
  201 PRINT*, 'CALCULATE THE In(z)Ga(1-z)As QUANTUM WELL--Z'
      Z=(1.5-DSQRT((1.5**2)-(4*0.4*(1.424-BANDEG))))/0.8
      PRINT*, 'BARRIER Al(x)Ga(y)In(1-x-y)As AND'
      PRINT*, 'CLADDING Al(xc)Ga(yc)In(1-xc-yc)As'
      PRINT*, 'X,Y--FOR BARRIER,XC,YC--FOR CLADDING'        
      X1=(EGB-0.75)/1.548  !Al
      Y1=1.D0-X1-0.53      !Ga
      XC=(EGC-0.75)/1.548  !Al
      YC=1.D0-XC-0.53      !Ga
      IF (YC.LT.0.000000D0) THEN
      YC=0.D0
      ELSE
      ENDIF
      WRITE(*,204) X1,Y1,XC,YC
  204 FORMAT(1X,'X=',F12.7,1X,'Y=',F12.7,1X,'XC=',F12.7,1X,'YC=',F12.7)
      PRINT*, 'Z=',Z
      LAYERX(N+1)=Z     !In
      LAYERY(N+1)=0.D0
      LAYERX(N)=XC      !Al
      LAYERY(N)=YC      !Ga
      LAYERX(1)=X1      !Al
      LAYERY(1)=Y1      !Ga
      BBLAY=BLAY/(N-1)
      INCREX=(LAYERX(N)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N)-LAYERY(1))/(N-1)
      DO I=2,N-1
      LAYERX(I)=LAYERX(I-1)+INCREX
      LAYERY(I)=LAYERY(I-1)+INCREY
      ENDDO
      AZ1=Z*A+(1-Z)*C
C      AZB9=(1-X1-Y1)*A+X1*E+Y1*C
      AZB9=5.8688
      EX=(AZB9-AZ1)/AZB9
      X(IC)=LAYERX(N+1)
      C11=X(IC)*8.329+(1-X(IC))*11.88
      C12=X(IC)*4.526+(1-X(IC))*5.38
      CAA=(C11-C12)/C11
      AA=X(IC)*(-6.08)+(1-X(IC))*(-9.77)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*EX)
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*EX)
      EGN=BANDEG+VB(IC)-VH(IC)
C      PRINT*, ' NEW BAND GAP=', EGN
      DO I=1,N
      EC1=(EGC-BANDEG)*0.72
      EC2=(EGB-BANDEG)*0.72
      DEG=(EC1-EC2)/(N-1)
      DELTAC(I)=EC2+(I-1)*DEG
      EV1=(EGC-BANDEG)*0.28
      EV2=(EGB-BANDEG)*0.28
      DVG=(EV1-EV2)/(N-1)
      DELTAV(I)=-(EV2+(I-1)*DVG)
      ENDDO
      STRAIN=(AZB9-AZ1)/AZB9	
      PRINT*,'CHECK STRAIN=',STRAIN
      PRINT*, '                                                  '
      CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE INGAALAS
      PARAMETER (II=200)  
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VH(500),LAYERX(500),VB(500),DELTAC(500)
      REAL*8 DELTAV(500),LAYERY(500),X(500),Y(500),INCREX,INCREY
      REAL*8 CT(500),CG(500),CA(500),T(500),TG(500),TC(500)  
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC 
      COMMON /BANDOUT1/ STRAIN,AZ1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
      IC=N+1
    1 REWIND 1
      REWIND 4    
      EGCC=1.9513   !AlAsSb LATTICE MATCHED TO InP
      PRINT*, 'INPUT I=1 FOR COMPRESS, 2 FOR TENSILE'
      PRINT*, 'I=#'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 2002
      ELSEIF (I.EQ.2) THEN
      GOTO 2003
      ELSE
      ENDIF
 2002 DO J=1,II
      CT(J)=0.53D0+((0.820D0-0.53D0)/II)*J
      CG(J)=0.75D0-((0.75D0-0.526D0)/II)*J
      CA(J)=1.548D0-((1.548D0-1.516D0)/II)*J
      YY=J*0.01
      WRITE(*,101) CT(J),CG(J),CA(J),YY,J 
  101 FORMAT(2X,'CT',F12.8,2X,'CG',F12.8,2X,'CA',F12.8,'COMPS',F6.3,'%',
     +       2X,I3)        
      ENDDO
      PRINT*, 'CALCULATE THE MATCHED BARRIER AND CLADDING COMPOSITION'
      PRINT*, 'X,Z--FOR BARRIER,XC,ZC--FOR CLADDING'
      Z1=(EGB-0.75D0)/1.548D0
      X1=1-0.53D0-Z1
      ZC=(EGC-0.75D0)/1.548D0
      XC=1-0.53D0-ZC
      ZCC=0.D0
      XCC=0.56
      IF (XC.LT.0.00000D0) THEN
      XC=0.D0
      ELSE
      ENDIF
      WRITE(*,1004) Z1,X1,ZC,XC
C      PRINT*, 'INPUT STRAIN "-" FOR COMPRESS'
      PRINT*, 'AND CALCULATE THE WELL COMPOSITION XW AND ZW'
      PRINT*, 'J= #'
      READ(*,*) J
      EX=-(J*0.0001)
      ZW=(BANDEG-CG(J))/CA(J)
      XW=1-CT(J)-ZW
      LAYERX(N+1)=XW
      LAYERY(N+1)=ZW
      LAYERX(N-1)=XC
      LAYERY(N-1)=ZC
      LAYERX(N)=XCC
      LAYERY(N)=ZCC
      LAYERX(1)=X1
      LAYERY(1)=Z1
      X(IC)=LAYERX(N+1)
      Y(IC)=LAYERY(N+1)
C      BBLAY=BLAY/(N-1)
      INCREX=(LAYERX(N-1)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N-1)-LAYERY(1))/(N-1)
      DO I=2,N-1
      LAYERY(I)=LAYERY(I-1)+INCREY
      LAYERX(I)=LAYERX(I-1)+INCREX
      ENDDO            
      C11=(1.-X(IC)-Y(IC))*8.329+Y(IC)*12.02+X(IC)*11.88
      C12=(1.-X(IC)-Y(IC))*4.526+Y(IC)*5.7+X(IC)*5.38
      CAA=(C11-C12)/C11
      AA=(1.-X(IC)-Y(IC))*(-6.08)+Y(IC)*(-8.11)+X(IC)*(-9.77)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*EX)  
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*EX)
      EGN=BANDEG+VB(IC)-VH(IC)
C     PRINT*, ' NEW BAND GAP=', EGN
      DO I=1,N-1
      EC1=(EGC-BANDEG)*0.72
      EC2=(EGB-BANDEG)*0.72
      DEG=(EC1-EC2)/(N-1)
      DELTAC(I)=EC2+(I-1)*DEG
      EV1=(EGC-BANDEG)*0.28
      EV2=(EGB-BANDEG)*0.28
      DVG=(EV1-EV2)/(N-1)
      DELTAV(I)=-(EV2+(I-1)*DVG)
      ENDDO               
      DELTAC(N)=(EGCC-BANDEG)*0.6
      DELTAV(N)=-(EGCC-BANDEG)*0.4    
      AZ3=(1-XW-ZW)*A+ZW*E+XW*C             
      PRINT*, 'INPUT STRAIN'                     
      READ(*,*) EX
      PRINT*, '                                                  '
      PRINT*, ' WRITE CONDUCTION BAND PARAMETERS INTO CBANDEG.DAT'
      WRITE(1,10) EX,AZ3*1.D-10
      WRITE(1,12) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)
      DO I=N-1,1,-1
         WRITE(1,12) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I)
      ENDDO
      WRITE(1,12) QLAY,LAYERX(N+1),LAYERY(N+1),VB(IC)
      DO I=1,N-1
         WRITE(1,12) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I)
      ENDDO
      WRITE(1,12) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)
      PRINT*, '                                                  '
      PRINT*, ' WRITE VALENCE BAND PARAMETERS INTO VBANDEG.DAT'
      WRITE(4,10) EX,AZ3*1.D-10
      WRITE(4,12) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)
      DO I=N-1,1,-1
         WRITE(4,12) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)
      ENDDO
      WRITE(4,12) QLAY,LAYERX(N+1),LAYERY(N+1),VH(IC)
      DO I=1,N-1
         WRITE(4,12) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)
      ENDDO
      WRITE(4,12) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)
      GOTO 4000
 2003 DO J=1,II
         T(J)=0.53D0-((0.53D0-0.230D0)/II)*J
	 TG(J)=0.75D0+((0.83D0-0.75D0)/II)*J
	 TC(J)=1.548D0+((1.588D0-1.548D0)/II)*J
	 XX=J*0.01
	 WRITE(*,102) T(J),TG(J),TC(J),XX,J
      ENDDO
      PRINT*, 'CALCULATE THE MATCHED BARRIER AND CLADDING COMPOSITION'
      PRINT*, 'X,Z--FOR BARRIER,XC,ZC--FOR CLADDING'
      Z1=(EGB-0.75D0)/1.548D0
      X1=1-0.53D0-Z1
      ZC=0.D0
      XC=0.56
      WRITE(*,1004) Z1,X1,ZC,XC
      PRINT*, 'AND CALCULATE THE WELL COMPOSITION XW AND ZW'
      PRINT*, 'J= #'
      READ(*,*) J
      EX=J*0.0001
      ZW=(BANDEG-TG(J))/TC(J)
      XW=1-T(J)-ZW
      LAYERX(N+1)=XW
      LAYERY(N+1)=ZW
      LAYERX(N)=XC
      LAYERY(N)=ZC
      LAYERX(1)=X1
      LAYERY(1)=Z1
      X(IC)=LAYERX(N+1)
      Y(IC)=LAYERY(N+1)
      INCREX=(LAYERX(N)-LAYERX(1))/(N-1)
      INCREY=(LAYERY(N)-LAYERY(1))/(N-1)
      DO I=2,N-1
	 LAYERY(I)=LAYERY(I-1)+INCREY
	 LAYERX(I)=LAYERX(I-1)+INCREX
      ENDDO
      C11=(1.-X(IC)-Y(IC))*8.329+Y(IC)*12.02+X(IC)*11.88
      C12=(1.-X(IC)-Y(IC))*4.526+Y(IC)*5.7+X(IC)*5.38
      CAA=(C11-C12)/C11
      AA=(1.-X(IC)-Y(IC))*(-6.08)+Y(IC)*(-8.11)+X(IC)*(-9.77)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*EX)
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*EX)
      EGN=BANDEG+VB(IC)-VH(IC)
      DO I=1,N-1
	 EC1=(EGC-BANDEG)*0.72
	 EC2=(EGB-BANDEG)*0.72
	 DEG=(EC1-EC2)/(N-1)
	 DELTAC(I)=EC2+(I-1)*DEG
	 EV1=(EGC-BANDEG)*0.28
	 EV2=(EGB-BANDEG)*0.28
	 DVG=(EV1-EV2)/(N-1)
	 DELTAV(I)=-(EV2+(I-1)*DVG)
      ENDDO
      DELTAC(N)=(EGCC-BANDEG)*0.6
      DELTAV(N)=-(EGCC-BANDEG)*0.4
      AZ3=(1-XW-ZW)*A+ZW*E+XW*C
      PRINT*, 'INPUT STRAIN'
      READ(*,*) EX
      PRINT*, '                                                  '
      PRINT*, ' WRITE CONDUCTION BAND PARAMETERS INTO CBANDEG.DAT'
      WRITE(1,10) EX,AZ3*1.D-10
      WRITE(1,12) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)
      DO I=N-1,1,-1
	 WRITE(1,12) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I)
      ENDDO
      WRITE(1,12) QLAY,LAYERX(N+1),LAYERY(N+1),VB(IC)
      DO I=1,N-1
	 WRITE(1,12) BBLAY,LAYERX(I),LAYERY(I),DELTAC(I)
      ENDDO
      WRITE(1,12) CLAY,LAYERX(N),LAYERY(N),DELTAC(N)
      PRINT*, '                                                  '
      PRINT*, ' WRITE VALENCE BAND PARAMETERS INTO VBANDEG.DAT'
      WRITE(4,10) EX,AZ3*1.D-10
      WRITE(4,12) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)
      DO I=N-1,1,-1
         WRITE(4,12) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)
      ENDDO
      WRITE(4,12) QLAY,LAYERX(N+1),LAYERY(N+1),VH(IC)
      DO I=1,N-1
         WRITE(4,12) BBLAY,LAYERX(I),LAYERY(I),DELTAV(I)
      ENDDO
      WRITE(4,12) CLAY,LAYERX(N),LAYERY(N),DELTAV(N)
      GOTO 4000
   10 FORMAT(2X,E12.6,2X,E12.6)     
   11 FORMAT(2(2X,E15.8),2(2X,F12.7)) 
   12 FORMAT(2(2X,E14.7),3(2X,F10.5))
  102 FORMAT(2X,'T',F12.8,2X,'TG',F12.8,2X,'TC',F12.8,'TENSS',F6.3,'%',
     +       2X,I3) 
 1004 FORMAT(1X,'Z=',F12.7,1X,'X=',F12.7,1X,'ZC=',F12.7,1X,'XC=',F12.7)
 4000 PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                          
      SUBROUTINE INGAAS5
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VH(500),LAYERX(500),VB(500),DELTAC(500)
      REAL*8 DELTAV(500),LAYERY(500),X(500),INCREX,INCREY
      INTEGER N,IC                           
      COMMON /ENERGY/ BANDEG,EGB,EGC,N
      COMMON /LATTICE/ A,B,C,D,E,F,G,H,HI
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC 
      common /bandout1/ strain,az1
      COMMON /BANDOUT2/ LAYERX,LAYERY,VB,VH,DELTAC,DELTAV
    1	CALL BANDIN 
      IC=N+1
      EGCC=1.9513
      PRINT*, 'THE In(z)Ga(1-z)As/Al(x)Ga(y)In(1-x-y)As/AlAsSb MATERIAL'
  211 PRINT*, 'CALCULATE THE In(z)Ga(1-z)As QUANTUM WELL--Z'
      Z=(1.5-DSQRT((1.5**2)-(4*0.4*(1.424-BANDEG))))/0.8
      PRINT*, 'BARRIER Al(x)Ga(y)In(1-x-y)As AND'
      PRINT*, 'CLADDING AlAs(xc)Sb(1-xc)'
      PRINT*, 'X,Y--FOR BARRIER,XC,YC--FOR CLADDING'        
      X1=(EGB-0.75)/1.548  !Al
      Y1=1.D0-0.53-X1      !Ga
      XC=(EGC-0.75)/1.548
      YC=1.D0-0.53-XC
      WRITE(*,204) X1,Y1,XC,YC
  204 FORMAT(1X,'X=',F12.7,1X,'Y=',F12.7,1X,'XC=',F12.7,1X,'YC=',F12.7)
      PRINT*, 'Z=',Z
      LAYERX(N+1)=Z    !In
      LAYERY(N+1)=0.D0      
      LAYERX(N-1)=XC      
      LAYERY(N-1)=YC      
      LAYERX(1)=X1
      LAYERY(1)=Y1
      LAYERX(N)=0.56                                     
      LAYERY(N)=0.E0
      INCREX=(LAYERX(N-1)-LAYERX(1))/(N-2)
      INCREY=(LAYERY(N-1)-LAYERY(1))/(N-2)
      DO I=2,N-2
      LAYERX(I)=LAYERX(I-1)+INCREX
      LAYERY(I)=LAYERY(I-1)+INCREY
      ENDDO
      AZ1=Z*A+(1-Z)*C
      AZB9=5.8688 
      EX=(AZB9-AZ1)/AZB9
      X(IC)=LAYERX(N+1)
      C11=X(IC)*8.329+(1-X(IC))*11.88
      C12=X(IC)*4.526+(1-X(IC))*5.38
      CAA=(C11-C12)/C11
      AA=X(IC)*(-6.08)+(1-X(IC))*(-9.77)
      VB(IC)=(2.D0*(2.D0/3.D0)*AA*CAA*EX)
      VH(IC)=(-2.D0*(1.D0/3.D0)*AA*CAA*EX)
      EGN=BANDEG+VB(IC)-VH(IC)
      DO I=1,N-1
      EC1=(EGC-BANDEG)*0.72
      EC2=(EGB-BANDEG)*0.72
      DEG=(EC1-EC2)/(N-2)
      DELTAC(I)=EC2+(I-1)*DEG
      EV1=(EGC-BANDEG)*0.28
      EV2=(EGB-BANDEG)*0.28
      DVG=(EV1-EV2)/(N-2)
      DELTAV(I)=-(EV2+(I-1)*DVG)
      ENDDO
      DELTAC(N)=(EGCC-BANDEG)*0.6
      DELTAV(N)=-(EGCC-BANDEG)*0.4
	STRAIN=(AZB9-AZ1)/AZB9 
      PRINT*,'CHECK STRAIN=',STRAIN
      PRINT*, '                                                  '
      CALL BANDOUT
	PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR EXIT'
      PRINT*, ' I= ?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSE
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CONDUCT
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LW(200),X(200),VB(200),MX(200),MC(200),DX(200)
      REAL*8 A(5),Y(200),WX(200),Z(200),MH(200),ML(200),D1(200)  
      INTEGER NL,IC,IX,IZ 
      COMMON /QUANTUM/ DX,VB,D1,NL       
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC 
      COMMON /MASS/ MC,MH,ML 
      COMMON /IFLAG/ IX
      COMMON /PARA1/ IZ
      COMMON /CENTER/ ICR,NUM
      EXTERNAL QW,ZOUT
      IX=1
      REWIND 1
      OPEN(1,FILE='cbandeg.dat',status='unknown')
	PRINT*, ' INPUT THE NUMBER OF QUANTUM WELLS NUM=?'
      READ(*,*) NUM
      PRINT*, ' INPUT TOTAL LAYERS FOR STRUCTURE--N ODD'
      PRINT*, ' INPUT N='
      READ(*,*) NL
      PRINT*, ' INPUT THE LOWEST POTENTIAL LAYER(1st Q-WELL) IC= ?'
      READ(*,*) IC
      PRINT*, ' INPUT THE SELECTED CENTER LAYER OF STRUCTURE ICR='
      READ(*,*) ICR
      READ(1,*) EX,AZ
      PRINT*, '*******************************************************'
      PRINT*, '  INPUT I=1  FOR AlGaAs '
      PRINT*, '        I=2  FOR InGaAsP '
      PRINT*, '        I=3  FOR In1-xGaxAs/InGaAsP/InP'
      PRINT*, '        I=4  FOR InGaAlAs/InGaAlAs'  
      PRINT*, '        I=5  FOR GaInP/(AlGa)0.5In0.5P/AlInP '
      PRINT*, '        I=6  FOR InGaAs/AlGaAs/AlGaAs '
      PRINT*, '        I=7  FOR InGaAs/InGaAsP/Ga0.51In0.49P(GaAs)'
      PRINT*, '        I=8  FOR AlyInxGa1-x-yAs/AlzGa1-zAs/GaAs '
      PRINT*, '        I=9  FOR InzGa1-zAs/AlxGayIn1-x-yAs/InP '
      PRINT*, '        I=10 FOR InGaAlAs/InGaAlAs/AlAsxSb1-x(InP)'
      PRINT*, '        I=11 FOR InzGa1-zAs/AlxGayIn1-x-yAs/AlAsxSb1-x'
      PRINT*, '  INPUT I= ?'
      PRINT*, '*******************************************************'
      READ(*,*) J
      IF (J.EQ.1) THEN
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      MX(I)=(0.067D0+0.083D0*X(I))
      ENDDO
      ELSEIF (J.EQ.2) THEN
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      MX(I)=(0.08D0-0.039D0*Y(I)) 
      ENDDO
      ELSEIF (J.EQ.3) THEN  
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      MX(I)=(0.08D0-0.039D0*Y(I))
      ENDDO
      DO J=1,NUM
      MX(IC+(J-1)*2)=(0.025D0*(1.D0-X(IC))+0.071D0*X(IC)-0.0163D0*X(IC)
     +       *(1.D0-X(IC)))  
      ENDDO     
      ELSEIF (J.EQ.4) THEN
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      MX(I)=(1.-X(I)-Y(I))*(0.0239)+Y(I)*(0.15)+X(I)*(0.067)
      ENDDO
      DO J=1,NUM
      MX(IC+(J-1)*2)=(1.-X(IC)-Y(IC))*(0.0239)+Y(IC)*(0.15)
     +       +X(IC)*(0.067)  
      ENDDO  
      ELSEIF (J.EQ.5) THEN
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      ENDDO
c      MX(NL)=X(NL)*(0.22D0)+(1.D0-X(NL))*(0.064D0)
c     MX(IC)=X(IC)*(0.15D0)+(1.D0-X(IC))*(0.077D0)
      DO J=1,NUM
c      MX(IC+(J-1)*2)=0.11
      MX(IC+(J-1)*2)=X(IC)*(0.15D0)+(1.D0-X(IC))*(0.077D0)
      ENDDO
      DO I=1,(NL-1)/2
      Z(I)=0.5D0*X(I)
      WX(I)=0.5D0*Y(I)
      MX(I)=(1.D0-Z(I)-WX(I))*(0.077D0)+Z(I)*(0.22D0)
     +      +WX(I)*(0.15D0)
c      MX(I)=(0.11+0.00915*X(I)-0.0024*X(I)**2)
      ENDDO
      DO I=((NL+1)/2)+1,NL
      Z(I)=0.5D0*X(I)
      WX(I)=0.5D0*Y(I)
      MX(I)=(1.D0-Z(I)-WX(I))*(0.077D0)+Z(I)*(0.22D0)
     +      +WX(I)*(0.15D0)
c      MX(I)=(0.11+0.00915*X(I)-0.0024*X(I)**2)
      ENDDO
      ELSEIF (J.EQ.6) THEN  
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      MX(I)=(0.067+0.083*X(I)) !Al
      ENDDO
      DO J=1,NUM
      MX(IC+(J-1)*2)=(0.067-0.04*X(IC))!In
      ENDDO
      ELSEIF (J.EQ.7) THEN
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      ENDDO
      MX(NL)=X(NL)*0.17+(1-X(NL))*0.077  ! GaxIn1-xP
      DO I=1,(NL-1)/2
      MX(I)=((1-X(I))*Y(I)*0.023+(1-X(I))*(1-Y(I))*0.077+X(I)*Y(I)*0.067
     +      +X(I)*(1-Y(I))*0.17)
      ENDDO
      DO J=1,NUM
      MX(IC+(J-1)*2)=(0.067-0.04*(1-X(IC)))
      ENDDO
      DO I=((NL+1)/2)+1,NL
      MX(I)=((1-X(I))*Y(I)*0.023+(1-X(I))*(1-Y(I))*0.077+X(I)*Y(I)*0.067
     +      +X(I)*(1-Y(I))*0.17)
      ENDDO
      ELSEIF (J.EQ.8) THEN
      DO I=1,NL
      READ(1,*) LW(I),Y(I),X(I),VB(I)
      MX(I)=(0.067+0.083*X(I))
      ENDDO
      MX(IC)=Y(IC)*0.15+X(IC)*0.0239+(1-X(IC)-Y(IC))*0.067
      ELSEIF (J.EQ.9) THEN   
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      MX(I)=(1.-X(I)-Y(I))*(0.0239)+X(I)*(0.15)+Y(I)*(0.067)
      ENDDO
      DO J=1,NUM
      MX(IC+(J-1)*2)=(0.025D0*X(IC)+0.071D0*(1-X(IC))-(0.0163D0*X(IC)
     +       *(1.D0-X(IC)))) 
      ENDDO
      ELSEIF (J.EQ.10) THEN
      DO I=1,NL
      READ(1,*) LW(I),Y(I),X(I),VB(I)
      MX(I)=(1.-X(I)-Y(I))*(0.0239)+Y(I)*(0.15)+X(I)*(0.067)
      ENDDO
      DO J=1,NUM
      MX(IC+(J-1)*2)=(1.-X(IC)-Y(IC))*(0.0239)+Y(IC)*(0.15)
     +       +X(IC)*(0.067)
      ENDDO
      ELSEIF (J.EQ.11) THEN
      DO I=1,NL
      READ(1,*) LW(I),X(I),Y(I),VB(I)
      MX(I)=(1.-X(I)-Y(I))*(0.0239)+X(I)*(0.15)+Y(I)*(0.067)
      ENDDO
C      MX(IC)=X(IC)*0.0239+(1-X(IC))*0.067 
      DO J=1,NUM
      MX(IC+(J-1)*2)=(0.025D0*X(IC))+0.071D0*(1-X(IC))-(0.0163D0*X(IC)
     +       *(1.D0-X(IC)))
      ENDDO    
      ELSE
      ENDIF                
      DO I=1,NL    
      DX(I)=LW(I)*DSQRT(MX(I)/3.804451D0)  
      MC(I)=MX(I)   
      VMAX=DMAX1(VMAX,VB(I))
      DO J=1,NUM
      VMIN=DMIN1(VB(IC+(J-1)*2),VB(I))
      A(1)=VMIN+1.D-5
      A(2)=VMAX-1.D-5
      A(3)=0.0001D0*(A(2)-A(1)) 
      A(4)=0.D0
      A(5)=0.D0
      ENDDO
      ENDDO
      D1(ICR-1)=-DX(ICR)/2.D0
      DO I=ICR-2,1,-1
      D1(I)=D1(I+1)-DX(I+1)
      ENDDO   
      D1(ICR)=DX(ICR)/2.D0
      DO I=ICR+1,NL-1,1
      D1(I)=D1(I-1)+DX(I)
      ENDDO
      CALL ZLOC(A,QW,ZOUT)
 9998 PRINT*, ' INPUT THE EIGENVALUE' 
      PRINT*, ' EIGEN VALUE='  
      READ(*,*) XE  
      CALL QW (XE,WELL)
      CALL STATE (XE) 
      PRINT*, ' INPUT NEW EIGENVALUE--> 1, BACK TO MAIN PAGE--> 2'
      PRINT*, ' SELECT=?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 9998
      ELSEIF (I.EQ.2) THEN
      GOTO 9999
      ELSE
      ENDIF
 9999 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE HEAVY
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 LW(200),VB(200),MC(200),DX(200),D1(200),X(200),Y(200)
      REAL*8 A(5),R1(200),R2(200),R3(200),MH(200),ML(200)  
      INTEGER NL,IC,IX,IZ,JJ 
      COMMON /QUANTUM/ DX,VB,D1,NL       
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC 
      COMMON /MASS/ MC,MH,ML
      COMMON /IFLAG/ IX
      COMMON /COEFF/ R1,R2,R3,LW
      COMMON /PARA1/ IZ
      COMMON /STRAIN/ DELTA
      COMMON /CENTER/ ICR,NUM
      COMMON /Z31/ JJ
      EXTERNAL QW,ZOUT
      IX=2
	PRINT*, ' INPUT THE NUMBER OF QUANTUM WELLS NUM=?'
      READ(*,*) NUM
      PRINT*, ' INPUT TOTAL LAYERS FOR STRUCTURE--N ODD'
      PRINT*, ' INPUT N='
      READ(*,*) NL
      PRINT*, ' INPUT THE HIGHEST POTENTIAL(1st Q-WELL) LAYER IC= ?'
      READ(*,*) IC 
      PRINT*, ' INPUT THE SELECTED CENTER OF THE STRUCTURE ICR=?'
      READ(*,*) ICR
      CALL COEF
      REWIND 4
      OPEN(4,FILE='vbandeg.dat',status='unknown')
      READ(4,*) EX,AZ
      DO I=1,NL
         READ(4,*) LW(I),X(I),Y(I),VB(I)
      ENDDO
      VBXX=VB(IC)
c      IF (JJ.NE.5) THEN
         DO I=1,NL
            DX(I)=LW(I)*DSQRT((1.D0/(R1(I)-2.D0*R2(I)))/3.804451D0)  
            MH(I)=1.D0/(R1(I)-2.D0*R2(I))
         ENDDO
c      ELSEIF (JJ.EQ.5) THEN
c         DO I=1,NL
c            MH(I)=0.62+0.05*X(I)
c            DX(I)=LW(I)*DSQRT(MH(I)/3.804451D0)
c         ENDDO
c        DO J=1,NUM
c           MH(IC+(J-1)*2)=0.62
c           DX(IC+(J-1)*2)=LW(IC)*DSQRT(MH(IC)/3.804451D0)
c        ENDDO
c      ELSE
c      ENDIF
      DO I=1,NL
      VMIN=DMIN1(VMIN,VB(I))
      DO J=1,NUM       
      VB(IC+(J-1)*2)=VBXX+DELTA/2.D0
      VMAX=DMAX1(VB(IC+(J-1)*2),VB(I))
      A(1)=VMIN+1.D-5
      A(2)=VMAX-1.D-5
      A(3)=0.0001D0*(A(2)-A(1)) 
      A(4)=0.D0
      A(5)=0.D0
      ENDDO
      ENDDO
      D1(ICR-1)=-DX(ICR)/2.D0
      DO I=ICR-2,1,-1
      D1(I)=D1(I+1)-DX(I+1)
      ENDDO   
      D1(ICR)=DX(ICR)/2.D0
      DO I=ICR+1,NL-1,1
      D1(I)=D1(I-1)+DX(I)
      ENDDO  
      CALL ZLOC(A,QW,ZOUT)
 9998 PRINT*, ' INPUT THE EIGENVALUE' 
      PRINT*, ' EIGEN VALUE='  
      READ(*,*) XE  
      CALL QW (XE,WELL)
      CALL STATE (XE)
      PRINT*, ' INPUT NEW EIGENVALUE--> 1, BACK TO MAIN PAGE--> 2'
      PRINT*, ' SELECT=?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 9998
      ELSEIF (I.EQ.2) THEN
      GOTO 9999
      ELSE
      ENDIF
 9999 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE LIGHT
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LW(200),VB(200),MC(200),DX(200),D1(200),X(200),Y(200)
      REAL*8 A(5),R1(200),R2(200),R3(200),ML(200),MH(200)  
      INTEGER NL,IC,IX,IZ,JJ 
      COMMON /QUANTUM/ DX,VB,D1,NL       
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC 
      COMMON /MASS/ MC,MH,ML
      COMMON /IFLAG/ IX
      COMMON /COEFF/ R1,R2,R3,LW  
      COMMON /PARA1/ IZ
      COMMON /STRAIN/ DELTA
      COMMON /Z31/ JJ
      COMMON /CENTER/ ICR,NUM
      EXTERNAL QW,ZOUT
      IX=2
      PRINT*, ' INPUT THE NUMBER OF QUANTUM WELLS NUM=?'
      READ(*,*) NUM
      PRINT*, ' INPUT TOTAL LAYERS FOR STRUCTURE--N ODD'
      PRINT*, ' INPUT N='
      READ(*,*) NL
      PRINT*, ' INPUT THE HIGHEST POTENTIAL(1st Q-WELL) LAYER IC= ?'
      READ(*,*) IC
      PRINT*, ' INPUT THE SELECTED CENTER OF THE STRUCTURE ICR=?'
      READ(*,*) ICR
      CALL COEF
      REWIND 4
       OPEN(4,FILE='vbandeg.dat',status='unknown')
      READ(4,*) EX,AZ
      DO I=1,NL
         READ(4,*) LW(I),X(I),Y(I),VB(I)
      ENDDO
      VBXX=VB(IC)
c      IF (JJ.NE.5) THEN
         DO I=1,NL
            DX(I)=LW(I)*DSQRT((1.D0/(R1(I)+2.D0*R2(I)))/3.804451D0)  
            ML(I)=1.D0/(R1(I)+2.D0*R2(I))
         ENDDO
c      ELSEIF (JJ.EQ.5) THEN
c         DO I=1,NL
c            ML(I)=0.11+0.03*X(I)
c            DX(I)=LW(I)*DSQRT(ML(I)/3.804451D0)
c            DO J=1,NUM
c               ML(IC+(J-1)*2)=0.11
c               DX(IC+(J-1)*2)=LW(IC)*DSQRT(ML(IC)/3.804451D0)
c            ENDDO
c        ENDDO
c      ELSE
c      ENDIF
      DO I=1,NL
      VMIN=DMIN1(VMIN,VB(I))
      DO J=1,NUM           
C      VB(IC+(J-1)*2)=VBXX-DELTA/2.D0
      VMAX=DMAX1(VB(IC+(J-1)*2),VB(I))
      A(1)=VMIN+1.D-5
      A(2)=VMAX-1.D-5
      A(3)=0.0001D0*(A(2)-A(1)) 
      A(4)=0.D0
      A(5)=0.D0
      ENDDO
      ENDDO
      D1(ICR-1)=-DX(ICR)/2.D0
      DO I=ICR-2,1,-1
      D1(I)=D1(I+1)-DX(I+1)
      ENDDO   
      D1(ICR)=DX(ICR)/2.D0
      DO I=ICR+1,NL-1,1
      D1(I)=D1(I-1)+DX(I)
      ENDDO  
      CALL ZLOC(A,QW,ZOUT)
 9998 PRINT*, ' INPUT THE EIGENVALUE' 
      PRINT*, ' EIGEN VALUE='  
      READ(*,*) XE  
      CALL QW (XE,WELL)
      CALL STATE (XE)
      PRINT*, ' INPUT NEW EIGENVALUE--> 1, BACK TO MAIN PAGE--> 2'
      PRINT*, ' SELECT=?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 9998
      ELSEIF (I.EQ.2) THEN
      GOTO 9999
      ELSE
      ENDIF
 9999 RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE COEF
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 R1(200),R2(200),R3(200),VB(200),DX(200),LW(200),X(200)
      REAL*8 Y(200),Z(200),WX(200),D1(200)
      INTEGER J,IC,NL,NA,JJ
      COMMON /QUANTUM/ DX,VB,D1,NL
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC
      COMMON /COEFF/ R1,R2,R3,LW 
      COMMON /PARA1/ NA
      COMMON /STRAIN/ DELTA
      COMMON /Z31/ JJ
      COMMON /CENTER/ ICR,NUM
      REWIND 4
      OPEN(4,FILE='vbandeg.dat',status='unknown')
      PRINT*, '*******************************************************'
      PRINT*, '  INPUT I=1  FOR AlGaAs '
      PRINT*, '        I=2  FOR InGaAsP ' 
      PRINT*, '        I=3  FOR In(1-x)Ga(x)As/InGaAsP/InP '  
      PRINT*, '        I=4  FOR InGaAlAs/InGaAlAs'  
      PRINT*, '        I=5  FOR GaInP/(AlGa)0.5In0.5P/AlInP '
      PRINT*, '        I=6  FOR InGaAs/AlGaAs/AlGaAs '
      PRINT*, '        I=7  FOR InGaAs/InGaAsP/Ga0.51In0.49P(GaAs)'
      PRINT*, '        I=8  FOR AlyInxGa1-x-yAs/AlzGa1-zAs/GaAs '
      PRINT*, '        I=9  FOR In(z)Ga(1-z)As/AlxGayIn1-x-yAs/InP '
      PRINT*, '        I=10 FOR InGaAlAs/InGaAlAs/AlAsxSb1-x(InP)'
      PRINT*, '        I=11 FOR InzGa1-zAs/AlxGayIn1-x-yAs/AlAsxSb1-x'
      PRINT*, '  INPUT I= ?'
      PRINT*, '*******************************************************'
      READ(*,*) J
      READ(4,*) EX,AZ
      JJ=J
      IF (J.EQ.1) THEN
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      R1(I)=6.85*(1.-X(I))+3.45*X(I)
      R2(I)=2.1*(1.-X(I))+0.68*X(I)
      R3(I)=2.9*(1.-X(I))+1.29*X(I)
      ENDDO
      C11=11.88+0.14*(X(IC))
      C12=5.38+0.32*(X(IC))
      C=(C11-C12)/C11
      CX=(C11+2.D0*C12)/C11
      AA=(-8.11)*X(IC)+(1.-X(IC))*(-9.77)
      BB=-1.7+0.2*X(IC)
      VW=(-2.0*(1.D0/3.D0)*AA*C*EX)
      DELTA=2.D0*BB*CX*EX
      VB(IC)=VW
      ELSEIF (J.EQ.2) THEN
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      ENDDO
      R1(1)=4.95D0
      R2(1)=1.65D0
      R3(1)=2.35D0
      R1(NL)=R1(1)
      R2(NL)=R2(1)
      R3(NL)=R3(1)
      DO I=2,NL-1
      R1(I)=(1.-X(I))*Y(I)*20.4+(1.-X(I))*(1.-Y(I))*4.95+X(I)*Y(I)*6.85
     +      +X(I)*(1.-Y(I))*4.05
      R2(I)=(1.-X(I))*Y(I)*8.3+(1.-X(I))*(1.-Y(I))*1.65+X(I)*Y(I)*2.1
     +      +X(I)*(1.-Y(I))*0.49
      R3(I)=(1.-X(I))*Y(I)*9.1+(1.-X(I))*(1.-Y(I))*2.35+X(I)*Y(I)*2.9
     +      +X(I)*(1.-Y(I))*1.25
      ENDDO
      C11=(1.-X(IC))*Y(IC)*8.329+(1.-X(IC))*(1.-Y(IC))*10.22
     +    +X(IC)*Y(IC)*11.88+X(IC)*(1.-Y(IC))*14.12
      C12=(1.-X(IC))*Y(IC)*4.526+(1.-X(IC))*(1.-Y(IC))*5.76
     +    +X(IC)*Y(IC)*5.38+X(IC)*(1.-Y(IC))*6.253
      CX=(C11+2.D0*C12)/C11
      BB=(1.-X(IC))*Y(IC)*(-1.8)+(1.-X(IC))*(1.-Y(IC))*(-2.0)+X(IC)
     +   *Y(IC)*(-1.7)+X(IC)*(1.-Y(IC))*(-1.5)
      DELTA=2.D0*BB*CX*EX
      ELSEIF (J.EQ.3) THEN                    ! In(1-x)GaxAs -- xGaAs+(1-x)InAs
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      ENDDO
      R1(1)=4.95D0
      R2(1)=1.65D0
      R3(1)=2.35D0
      R1(NL)=R1(1)
      R2(NL)=R2(1)
      R3(NL)=R3(1) 
      DO I=2,NL-1
      R1(I)=(1.-X(I))*Y(I)*20.4+(1.-X(I))*(1.-Y(I))*4.95+X(I)*Y(I)*6.85
     +      +X(I)*(1.-Y(I))*4.05
      R2(I)=(1.-X(I))*Y(I)*8.3+(1.-X(I))*(1.-Y(I))*1.65+X(I)*Y(I)*2.1
     +      +X(I)*(1.-Y(I))*0.49
      R3(I)=(1.-X(I))*Y(I)*9.1+(1.-X(I))*(1.-Y(I))*2.35+X(I)*Y(I)*2.9
     +      +X(I)*(1.-Y(I))*1.25
      ENDDO
      DO J=1,NUM 
      R1(IC+(J-1)*2)=X(IC)*6.85D0+20.4D0*(1.D0-X(IC))
      R2(IC+(J-1)*2)=X(IC)*2.1D0+8.3D0*(1.D0-X(IC))
      R3(IC+(J-1)*2)=X(IC)*2.9D0+9.1D0*(1.D0-X(IC))
      ENDDO
      C11=X(IC)*11.88+(1.D0-X(IC))*8.329D0
      C12=X(IC)*5.38D0+(1.D0-X(IC))*4.526D0
      CX=(C11+2.D0*C12)/C11
      BB=X(IC)*(-1.7D0)+(1.D0-X(IC))*(-1.8D0)
      DELTA=2.D0*BB*CX*EX   
      ELSEIF (J.EQ.4) THEN
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      R1(I)=(1.-X(I)-Y(I))*20.4+Y(I)*3.45+X(I)*6.85
      R2(I)=(1.-X(I)-Y(I))*8.3+Y(I)*0.68+X(I)*2.1
      R3(I)=(1.-X(I)-Y(I))*9.1+Y(I)*1.29+X(I)*2.9
      ENDDO
      C11=(1.-X(IC)-Y(IC))*8.329+Y(IC)*12.02+X(IC)*11.88
      C12=(1.-X(IC)-Y(IC))*4.526+Y(IC)*5.7+X(IC)*5.38
      CX=(C11+2.D0*C12)/C11
      BB=(1.-X(IC)-Y(IC))*(-1.8)+Y(IC)*(-1.5)+X(IC)*(-1.7)
      DELTA=2.D0*BB*CX*EX  
      ELSEIF (J.EQ.5) THEN
      DO I=1,NL 
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      ENDDO
c      Z(NL)=0.5D0*X(NL)        ! FOR AlGaInP CLADDING
c      WX(NL)=0.5D0*Y(NL)
c      R1(NL)=(1-Z(NL)-WX(NL))*6.35D0+Z(NL)*3.47D0+WX(NL)*4.2D0
c      R2(NL)=(1-Z(NL)-WX(NL))*2.08D0+Z(NL)*0.06D0+WX(NL)*0.98D0
c      R3(NL)=(1-Z(NL)-WX(NL))*2.76D0+Z(NL)*1.15D0+WX(NL)*1.66D0
c     DO I=1,(NL-1)/2
c     Z(I)=0.5D0*X(I)
c     WX(I)=0.5D0*Y(I)
c     R1(I)=(1-Z(I)-WX(I))*6.35D0+Z(I)*3.47D0+WX(I)*4.2D0
c     R2(I)=(1-Z(I)-WX(I))*2.08D0+Z(I)*0.06D0+WX(I)*0.98D0
c     R3(I)=(1-Z(I)-WX(I))*2.76D0+Z(I)*1.15D0+WX(I)*1.66D0
c     ENDDO         
c     DO I=((NL+1)/2)+1,NL
c     Z(I)=0.5D0*X(I)
c     WX(I)=0.5D0*Y(I)
c     R1(I)=(1-Z(I)-WX(I))*6.35D0+Z(I)*3.47D0+WX(I)*4.2D0
c     R2(I)=(1-Z(I)-WX(I))*2.08D0+Z(I)*0.06D0+WX(I)*0.98D0
c     R3(I)=(1-Z(I)-WX(I))*2.76D0+Z(I)*1.15D0+WX(I)*1.66D0
c     ENDDO
      DO I=1,NL
      Z(I)=0.5D0*X(I)
      WX(I)=0.5D0*Y(I)
      R1(I)=(1-Z(I)-WX(I))*6.35D0+Z(I)*3.47D0+WX(I)*4.2D0
      R2(I)=(1-Z(I)-WX(I))*2.08D0+Z(I)*0.06D0+WX(I)*0.98D0
      R3(I)=(1-Z(I)-WX(I))*2.76D0+Z(I)*1.15D0+WX(I)*1.66D0
      ENDDO
      DO J=1,NUM
      R1(IC+(J-1)*2)=X(IC)*4.2D0+(1.D0-X(IC))*6.35D0
      R2(IC+(J-1)*2)=X(IC)*0.98D0+(1.D0-X(IC))*2.08D0
      R3(IC+(J-1)*2)=X(IC)*1.66D0+(1.D0-X(IC))*2.76D0
      ENDDO
      C11=10.22+3.9*X(IC)
      C12=5.76+0.49*X(IC)
      CX=(C11+2.D0*C12)/C11
      BB=X(IC)*(-1.8D0)+(1.D0-X(IC))*(-2.0D0)
C      BB=-(2.0-0.2*X(IC))
      DELTA=2.D0*BB*CX*EX
      ELSEIF (J.EQ.6) THEN
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      R1(I)=6.85*(1-X(I))+3.45*X(I)
      R2(I)=2.1*(1-X(I))+0.68*X(I)
      R3(I)=2.9*(1-X(I))+1.29*X(I)
      ENDDO
      DO J=1,NUM
      R1(IC+(J-1)*2)=6.85*(1-X(IC))+19.67*X(IC)
      R2(IC+(J-1)*2)=2.1*(1-X(IC))+8.37*X(IC)
      R3(IC+(J-1)*2)=2.9*(1-X(IC))+9.29*X(IC)
      ENDDO
      C11=X(IC)*8.33+(1-X(IC))*11.88
      C12=X(IC)*4.53+(1-X(IC))*5.38
      CX=(C11+2.D0*C12)/C11
      BB=X(IC)*(-1.8)+(1-X(IC))*(-1.7)
      DELTA=2.D0*BB*CX*EX
      ELSEIF (J.EQ.7) THEN
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      ENDDO
      R1(NL)=X(NL)*4.2+(1-X(NL))*6.35
      R2(NL)=X(NL)*0.98+(1-X(NL))*2.08
      R1(NL)=X(NL)*1.66+(1-X(NL))*2.76
      DO J=1,NUM
      R1(IC+(J-1)*2)=6.85*X(IC)+19.67*(1-X(IC))
      R2(IC+(J-1)*2)=2.1*X(IC)+(1-X(IC))*8.37
      R3(IC+(J-1)*2)=2.9*X(IC)+9.29*(1-X(IC))
      ENDDO
      DO I=1,(NL-1)/2
      R1(I)=(1-X(I))*Y(I)*19.67+(1-X(I))*(1-Y(I))*6.35+X(I)*Y(I)*6.85
     +      +X(I)*(1-Y(I))*4.2
      R2(I)=(1-X(I))*Y(I)*8.37+(1-X(I))*(1-Y(I))*2.08+X(I)*Y(I)*2.1
     +      +X(I)*(1-Y(I))*0.98
      R3(I)=(1-X(I))*Y(I)*9.29+(1-X(I))*(1-Y(I))*2.76+X(I)*Y(I)*2.9
     +      +X(I)*(1-Y(I))*1.66
      ENDDO
      DO I=((NL+1)/2)+1,NL
      R1(I)=(1-X(I))*Y(I)*19.67+(1-X(I))*(1-Y(I))*6.35+X(I)*Y(I)*6.85
     +      +X(I)*(1-Y(I))*4.2
      R2(I)=(1-X(I))*Y(I)*8.37+(1-X(I))*(1-Y(I))*2.08+X(I)*Y(I)*2.1
     +      +X(I)*(1-Y(I))*0.98
      R3(I)=(1-X(I))*Y(I)*9.29+(1-X(I))*(1-Y(I))*2.76+X(I)*Y(I)*2.9
     +      +X(I)*(1-Y(I))*1.66
      ENDDO
      C11=X(IC)*11.88+(1-X(IC))*8.33
      C12=X(IC)*5.38+(1-X(IC))*4.53
      CX=(C11+2.D0*C12)/C11
      BB=X(IC)*(-1.7)+(1-X(IC))*(-1.8)
      DELTA=2.D0*BB*CX*EX
      ELSEIF (J.EQ.8) THEN
      DO I=1,NL
      READ(4,*) LW(I),Y(I),X(I),VB(I)
      R1(I)=6.85*(1.-X(I))+3.45*X(I)
      R2(I)=2.1*(1.-X(I))+0.68*X(I)
      R3(I)=2.9*(1.-X(I))+1.29*X(I)
      ENDDO
      DO J=1,NUM
      R1(IC+(J-1)*2)=Y(IC)*3.45+X(IC)*19.67+(1-X(IC)-Y(IC))*6.85
      R2(IC+(J-1)*2)=Y(IC)*0.68+X(IC)*8.37+(1-X(IC)-Y(IC))*2.1
      R3(IC+(J-1)*2)=Y(IC)*1.29+X(IC)*9.29+(1-X(IC)-Y(IC))*2.9
      ENDDO
      C11=Y(IC)*12.02+X(IC)*8.33+(1-X(IC)-Y(IC))*11.88
      C12=Y(IC)*5.7+X(IC)*4.53+(1-X(IC)-Y(IC))*5.38
      CX=(C11+2.D0*C12)/C11
      BB=Y(IC)*(-1.5)+X(IC)*(-1.8)+(1-X(IC)-Y(IC))*(-1.7)
      DELTA=2.D0*BB*CX*EX
      ELSEIF (J.EQ.9) THEN         ! InxGa(1-x)As -- xInAs+(1-x)GaAs
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      R1(I)=(1.-X(I)-Y(I))*20.4+X(I)*3.45+Y(I)*6.85
      R2(I)=(1.-X(I)-Y(I))*8.3+X(I)*0.68+Y(I)*2.1
      R3(I)=(1.-X(I)-Y(I))*9.1+X(I)*1.29+Y(I)*2.9
      ENDDO
      DO J=1,NUM
      R1(IC+(J-1)*2)=X(IC)*20.4+(1-X(IC))*6.85
      R2(IC+(J-1)*2)=X(IC)*8.3+(1-X(IC))*2.1
      R3(IC+(J-1)*2)=X(IC)*9.1+(1-X(IC))*2.9
      ENDDO
      C11=X(IC)*8.329+(1-X(IC))*11.88
      C12=X(IC)*4.526+(1-X(IC))*5.38
      CX=(C11+2.D0*C12)/C11
      BB=X(IC)*(-1.8)+(1-X(IC))*(-1.7)
      DELTA=2.D0*BB*CX*EX 
      ELSEIF (J.EQ.10) THEN
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      R1(I)=(1.-X(I)-Y(I))*20.4+X(I)*3.45+Y(I)*6.85
      R2(I)=(1.-X(I)-Y(I))*8.3+X(I)*0.68+Y(I)*2.1
      R3(I)=(1.-X(I)-Y(I))*9.1+X(I)*1.29+Y(I)*2.9
      ENDDO
      R1(NL)=X(NL)*3.45+(1-X(NL))*4.93
      R2(NL)=X(NL)*0.68+(1-X(NL))*2.59
C      R3(NL)=
      C11=(1.-X(IC)-Y(IC))*8.329+Y(IC)*12.02+X(IC)*11.88
      C12=(1.-X(IC)-Y(IC))*4.526+Y(IC)*5.7+X(IC)*5.38
      CX=(C11+2.D0*C12)/C11
      BB=(1.-X(IC)-Y(IC))*(-1.8)+Y(IC)*(-1.5)+X(IC)*(-1.7)
      DELTA=2.D0*BB*CX*EX  
      ELSEIF (J.EQ.11) THEN        !InxGa(1-x)As -- xInAs+(1-x)GaAs
      DO I=1,NL
      READ(4,*) LW(I),X(I),Y(I),VB(I)
      R1(I)=(1.-X(I)-Y(I))*20.4+X(I)*3.45+Y(I)*6.85
      R2(I)=(1.-X(I)-Y(I))*8.3+X(I)*0.68+Y(I)*2.1
      R3(I)=(1.-X(I)-Y(I))*9.1+X(I)*1.29+Y(I)*2.9
      ENDDO
      DO J=1,NUM
      R1(IC+(J-1)*2)=X(IC)*20.4+(1-X(IC))*6.85
      R2(IC+(J-1)*2)=X(IC)*8.3+(1-X(IC))*2.1
      R3(IC+(J-1)*2)=X(IC)*9.1+(1-X(IC))*2.9
      ENDDO
      R1(NL)=X(NL)*3.45+(1-X(NL))*4.93
      R2(NL)=X(NL)*0.68+(1-X(NL))*2.59
C      R3(NL)=
      C11=X(IC)*8.329+(1-X(IC))*11.88
      C12=X(IC)*4.526+(1-X(IC))*5.38
      CX=(C11+2.D0*C12)/C11
      BB=X(IC)*(-1.8)+(1-X(IC))*(-1.7)
      DELTA=2.D0*BB*CX*EX 
      ELSE
      ENDIF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
      SUBROUTINE SPLINE (CEL,CEY,CEX)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,I,IEND
      REAL*8 X(100),Y(100),A(100),B(100),C(100),S(100),CEL(100)
      REAL*8 CEY(100),CEX(100)
      DATA N,IEND/11,1/
      DATA X/0.65,0.658,0.668,0.68,0.7,0.72,0.747,0.771,0.8,0.8375,
     +      0.875,89*0.0/
      DATA Y/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,89*0.0/
      DATA S/100*0.0/
      IEND=3
      PRINT '(///)'
      PRINT*, 'OUTPUT FOR IEND= ',IEND
      PRINT '(///)'
C
      CALL CUBSPL2 (X,Y,IEND,N,A,B,C,S)
      PRINT 98
      PRINT 99
      PRINT 101 ,(S(I),I=1,N)
      PRINT 98
      PRINT 200     
      DO I=1,100
      U=0.65+(I-1)*(0.875-0.65)/(100-1)
      TEMP1=SEVAL(N,U,X,Y,A,B,C)
      PRINT 100, U,TEMP1
      CEL(I)=U
      CEY(I)=TEMP1
      CEX(I)=(0.1894*CEY(I)+0.2156)/(0.418-0.013*CEY(I))
      ENDDO      
      PRINT 98
   98 FORMAT(//T2,60('*')/)
   99 FORMAT('  THE SECOND DERIVATIVES S(I) ARE;'/)
  100 FORMAT(3X,F10.4,4X,F10.4)
  101 FORMAT(10F9.4)
  200 FORMAT(I5,'X',T15,'SPLINE(X)'//)
      RETURN
      END
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CUBSPL2(X,Y,IEND,N,A,B,C,S)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(N),Y(N),S(N),A(N),B(N),C(N),SMATRIX(0:20,4)
      REAL*8 DX1,DY1,DX2,DY2,DXN1,DXN2,H(20)
      INTEGER N,IEND,NM1,NM2,I,J,FIRST,LAST
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

C     -------------------------------------------------------------
      IF (IEND.EQ.2) THEN
   50 SMATRIX(1,2)=SMATRIX(1,2)+X(2)-X(1)
      SMATRIX(NM2,2)=SMATRIX(NM2,2)+X(N)-X(NM1)
      ELSE IF (IEND.EQ.3) THEN
   80 DX1=X(2)-X(1)
      DX2=X(3)-X(2)
      SMATRIX(1,2)=(DX1+DX2)*(DX1+2.D0*DX2)/DX2
      SMATRIX(1,3)=(DX2*DX2-DX1*DX1)/DX2
      DXN2=X(NM1)-X(NM2)
      DXN1=X(N)-X(NM1)
      SMATRIX(NM2,1)=(DXN2*DXN2-DXN1*DXN1)/DXN2
      SMATRIX(NM2,2)=(DXN1+DXN2)*(DXN1+2.D0*DXN2)/DXN2
      ENDIF
C
C
      IF(IEND.EQ.4) THEN
      DX1=X(2)-X(1)
      DY1=(Y(2)-Y(1))/DX1
      SMATRIX(0,1)=1.D0
      SMATRIX(0,2)=2.D0*DX1
      SMATRIX(0,3)=DX1
      SMATRIX(0,4)=(DY1-S(1))*6.D0
      DX1=X(N)-X(N-1)
      DY1=(Y(N)-Y(N-1))/DX1
      SMATRIX(NM1,1)=DX1
      SMATRIX(NM1,2)=2.D0*DX1
      SMATRIX(NM1,3)=0.D0
      SMATRIX(NM1,4)=(S(N)-DY1)*6.D0
      FIRST=1
      LAST=N-1
      ENDIF
C
      DO I=FIRST-1,LAST
         PRINT*,  (SMATRIX(I,J),J=1,4)
      ENDDO
C
      DO I=FIRST,LAST
         SMATRIX(I,1)=SMATRIX(I,1)/SMATRIX(I-1,2)
         SMATRIX(I,2)=SMATRIX(I,2)-SMATRIX(I,1)*SMATRIX(I-1,3)
         SMATRIX(I,4)=SMATRIX(I,4)-SMATRIX(I,1)*SMATRIX(I-1,4)
      ENDDO
C
      SMATRIX(LAST,4)=SMATRIX(LAST,4)/SMATRIX(LAST,2)
      DO J=LAST-1,FIRST-1,-1
         SMATRIX(J,4)=(SMATRIX(J,4)-SMATRIX(J,3)*SMATRIX(J+1,4))
     +   /SMATRIX(J,2)
      ENDDO
C
      DO I=FIRST-1,LAST
         S(I+1)=SMATRIX(I,4)
      ENDDO
C
      IF(IEND.EQ.1) THEN
C
      S(1)=0.D0
      S(N)=0.D0
      ELSEIF (IEND.EQ.2) THEN
C
      S(1)=S(2)
      S(N)=S(N-1)
C
      ELSEIF (IEND.EQ.3) THEN
      S(1)=((DX1+DX2)*S(2)-DX1*S(3))/DX2
      S(N)=((DXN2+DXN1)*S(NM1)-DXN1*S(NM2))/DXN2
      ENDIF
C
      PRINT 99
   99 FORMAT(/T2,60('*')//T10,'THE CUBIC POLYNOMIALS,G(X)'
     +,'DEFINED ON THE INTERVALS'//)
      PRINT 101
C
      DO I=1,N-1
         H(I)=X(I+1)-X(I)
         A(I)=(S(I+1)-S(I))/(6.D0*H(I))
         B(I)=S(I)/2.D0
         C(I)=((Y(I+1)-Y(I))/H(I))-((2*H(I)*S(I)+H(I)*S(I+1))/6.D0)
C
         PRINT 102,I,A(I),B(I),C(I),Y(I)
  101    FORMAT(T5,'T',T17,'A',T27,'B',T37,'C',T47,'D'/)
  102    FORMAT(T5,I1,T12,F8.4,T22,F8.4,T32,F8.4,T42,F8.4)
      ENDDO
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION SEVAL (N,U,X,Y,A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 U,X(N),Y(N),A(N),B(N),C(N)
C
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
C
      DX=U-X(I)
      SEVAL=Y(I)+DX*(C(I)+DX*(B(I)+DX*A(I)))
      RETURN
      ENDIF
      END            
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE QW(E,WELL)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 V(200),VB(200),H(100),E,WELL
      REAL*8 D1(200),LT(2,2,100),RT(2,2,100),P(2,2),Q(2,2)
      REAL*8 PT(2,2),QT(2,2),CBX,MM(200)
      REAL*8 MC(200),MH(200),ML(200),DX(200)
      INTEGER NL,IC,IZ
      COMMON /QUANTUM/ DX,VB,D1,NL
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC 
      COMMON /MASS/ MC,MH,ML 
      COMMON /FUN/ CBX
      COMMON /STRAIN/ DELTA
      COMMON /PARA1/ IZ
      COMMON /CENTER/ ICR,NUM
      DO I=1,NL     
      IF (IZ.EQ.2) THEN
      V(I)=VB(I) 
      MM(I)=MC(I)
      ELSEIF (IZ.EQ.3) THEN   
      V(I)=VB(I)
      MM(I)=MH(I)                  
      ELSEIF (IZ.EQ.4) THEN       
      V(I)=VB(I)
      MM(I)=ML(I)
      ELSE
      ENDIF    
      ENDDO             
      IF (IZ.EQ.2) THEN
      H(1)=DSQRT(V(1)-E)
      H(NL)=DSQRT(V(NL)-E)
      ELSE
      H(1)=DSQRT(E-V(1))
      H(NL)=DSQRT(E-V(NL))   
      ENDIF
CCCCC      
      LT(1,1,1)=DEXP(H(1)*D1(1))
      LT(1,2,1)=0.D0
      LT(2,1,1)=H(1)*DEXP(H(1)*D1(1))/MM(1)
      LT(2,2,1)=0.D0
      RT(1,1,NL)=0.D0
      RT(1,2,NL)=DEXP((-1.D0)*H(NL)*D1(NL-1))
      RT(2,1,NL)=0.D0
      RT(2,2,NL)=(-1.D0)*H(NL)*DEXP(-H(NL)*D1(NL-1))/MM(NL)    
CCCCC      
      DO I=2,NL-1         
      IF (IZ.EQ.2) THEN
      H(I)=(E-V(I))    
      ELSE
      H(I)=(V(I)-E)
      ENDIF
      IF (H(I).GT.0.D0) THEN
         H(I)=DSQRT(H(I))
         IF (I.EQ.ICR) THEN
             LT(1,1,I)=DCOS(H(I)*D1(I-1))
             LT(1,2,I)=((-DSIN(H(I)*D1(I-1)))*MM(I))/H(I)
             LT(2,1,I)=DSIN(H(I)*D1(I-1))
             LT(2,2,I)=(DCOS(H(I)*D1(I-1))*MM(I))/H(I)   
CCCCC        
             RT(1,1,I)=DCOS(H(I)*D1(I))
             RT(1,2,I)=((-DSIN(H(I)*D1(I)))*MM(I))/H(I)
             RT(2,1,I)=DSIN(H(I)*D1(I))
             RT(2,2,I)=(DCOS(H(I)*D1(I))*MM(I))/H(I)     
         ELSE
             LT(1,1,I)=DCOS(H(I)*(D1(I)-D1(I-1)))
             LT(1,2,I)=(DSIN(H(I)*(D1(I)-D1(I-1)))*MM(I))/H(I)
             LT(2,1,I)=-H(I)*DSIN(H(I)*(D1(I)-D1(I-1)))/MM(I)
             LT(2,2,I)=DCOS(H(I)*(D1(I)-D1(I-1)))   
CCCCC        
             RT(1,1,I)=DCOS(H(I)*(D1(I-1)-D1(I)))
             RT(1,2,I)=(DSIN(H(I)*(D1(I-1)-D1(I)))*MM(I))/H(I)
             RT(2,1,I)=-H(I)*DSIN(H(I)*(D1(I-1)-D1(I)))/MM(I)
             RT(2,2,I)=DCOS(H(I)*(D1(I-1)-D1(I)))
         ENDIF
      ELSE
         H(I)=DSQRT(-H(I))
         IF (I.EQ.ICR) THEN
             LT(1,1,I)=DCOSH(H(I)*D1(I-1))
             LT(1,2,I)=(-DSINH(H(I)*D1(I-1))*MM(I))/H(I)
             LT(2,1,I)=-DSINH(H(I)*D1(I-1))
             LT(2,2,I)=(DCOSH(H(I)*D1(I-1))*MM(I))/H(I) 
CCCCC        
             RT(1,1,I)=DCOSH(H(I)*D1(I))
             RT(1,2,I)=(-DSINH(H(I)*D1(I))*MM(I))/H(I)
             RT(2,1,I)=-DSINH(H(I)*D1(I))
             RT(2,2,I)=(DCOSH(H(I)*D1(I))*MM(I))/H(I)
         ELSE                           
             LT(1,1,I)=DCOSH(H(I)*(D1(I)-D1(I-1)))
             LT(1,2,I)=(DSINH(H(I)*(D1(I)-D1(I-1)))*MM(I))/H(I)
             LT(2,1,I)=H(I)*DSINH(H(I)*(D1(I)-D1(I-1)))/MM(I)
             LT(2,2,I)=DCOSH(H(I)*(D1(I)-D1(I-1)))   
CCCCC        
             RT(1,1,I)=DCOSH(H(I)*(D1(I-1)-D1(I)))
             RT(1,2,I)=(DSINH(H(I)*(D1(I-1)-D1(I)))*MM(I))/H(I)
             RT(2,1,I)=H(I)*DSINH(H(I)*(D1(I-1)-D1(I)))/MM(I)
             RT(2,2,I)=DCOSH(H(I)*(D1(I-1)-D1(I)))
         ENDIF
      ENDIF
      ENDDO
CCCCC      
      P(1,1)=LT(1,1,1)
      P(1,2)=LT(1,2,1)
      P(2,1)=LT(2,1,1)
      P(2,2)=LT(2,2,1)
CCCCC 
      Q(1,1)=RT(1,1,NL)
      Q(1,2)=RT(1,2,NL)
      Q(2,1)=RT(2,1,NL)
      Q(2,2)=RT(2,2,NL)
CCCCC 
      DO I=2,ICR,1
      PT(1,1)=LT(1,1,I)*P(1,1)+LT(1,2,I)*P(2,1)
      PT(2,1)=LT(2,1,I)*P(1,1)+LT(2,2,I)*P(2,1)
      PT(1,2)=LT(1,1,I)*P(1,2)+LT(1,2,I)*P(2,2)
      PT(2,2)=LT(2,1,I)*P(1,2)+LT(2,2,I)*P(2,2)
      P(1,1)=PT(1,1)
      P(2,1)=PT(2,1)
      P(1,2)=PT(1,2)
      P(2,2)=PT(2,2)
      ENDDO
CCCCC 
      DO I=NL-1,ICR,-1
      QT(1,1)=RT(1,1,I)*Q(1,1)+RT(1,2,I)*Q(2,1)
      QT(2,1)=RT(2,1,I)*Q(1,1)+RT(2,2,I)*Q(2,1)
      QT(1,2)=RT(1,1,I)*Q(1,2)+RT(1,2,I)*Q(2,2)
      QT(2,2)=RT(2,1,I)*Q(1,2)+RT(2,2,I)*Q(2,2)
      Q(1,1)=QT(1,1)
      Q(2,1)=QT(2,1)
      Q(1,2)=QT(1,2)
      Q(2,2)=QT(2,2)
      ENDDO
      CBX=(P(1,1)-P(2,1))/(Q(1,2)-Q(2,2))
      WELL=Q(1,2)*P(2,1)-P(1,1)*Q(2,2)
      RETURN
      END    
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE STATE (E) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 V(200),VB(200),H(100),AX(100),BX(100),CBX
      REAL*8 LT(2,2,100),RT(2,2,100),P(2,2),Q(2,2),DZZ(100),HH2(100)
      REAL*8 PT(2,2),QT(2,2),LLT(2,2,100),RRT(2,2,100),XP(2,2),YP(2,2) 
      REAL*8 MM(200),MC(200),MH(200),ML(200),PHC(1200),CONF(100)
      REAL*8 DZ1(1200),DZ2(1200),DX(200),D1(200),XC(100),HH1(100)
      INTEGER NL,IC,IZ
      CHARACTER OUTPUT*40
      COMMON /QUANTUM/ DX,VB,D1,NL
      COMMON /LAYER/ CLAY,BLAY,QLAY,BBLAY,IC 
      COMMON /MASS/ MC,MH,ML 
      COMMON /FUN/CBX
      COMMON /STRAIN/ DELTA 
      COMMON /PARA1/ IZ
      COMMON /CFM1/ AXX,BXX,HH1X,HH2X
      COMMON /CENTER/ ICR,NUM
      EXTERNAL CFM
      MIN=5 
      MOUT=7
      AX(1)=1.D0
      BX(NL)=CBX     
      DO I=1,NL     
         IF (IZ.EQ.2) THEN
            V(I)=VB(I)
            MM(I)=MC(I)
            ELSEIF (IZ.EQ.3) THEN   
            V(I)=VB(I)
            MM(I)=MH(I)                  
         ELSEIF (IZ.EQ.4) THEN       
            V(I)=VB(I)
            MM(I)=ML(I)
         ELSE
         ENDIF    
      ENDDO                             
      IF (IZ.EQ.2) THEN
         H(1)=DSQRT(V(1)-E)
         H(NL)=DSQRT(V(NL)-E)
      ELSE
         H(1)=DSQRT(E-V(1))
         H(NL)=DSQRT(E-V(NL))
      ENDIF
CCCCC      
      LT(1,1,1)=DEXP(H(1)*D1(1))
      LT(1,2,1)=0.D0
      LT(2,1,1)=H(1)*DEXP(H(1)*D1(1))/MM(1)
      LT(2,2,1)=0.D0
CCCCC      
      RT(1,1,NL)=0.D0
      RT(1,2,NL)=DEXP((-1.D0)*H(NL)*D1(NL-1))
      RT(2,1,NL)=0.D0
      RT(2,2,NL)=(-1.D0)*H(NL)*DEXP(-H(NL)*D1(NL-1))/MM(NL)    
CCCCC      
      DO I=2,NL-1  
         IF (IZ.EQ.2) THEN
            H(I)=(E-V(I))
            HH1(I)=H(I) 
         ELSE
            H(I)=(V(I)-E)
            HH2(I)=H(I)
         ENDIF
         IF (H(I).GT.0.D0) THEN
            H(I)=DSQRT(H(I))
            LT(1,1,I)=DCOS(H(I)*(D1(I)-D1(I-1)))
            LT(1,2,I)=(DSIN(H(I)*(D1(I)-D1(I-1)))*MM(I))/H(I)
            LT(2,1,I)=-H(I)*DSIN(H(I)*(D1(I)-D1(I-1)))/MM(I)
            LT(2,2,I)=DCOS(H(I)*(D1(I)-D1(I-1)))   
CCCCC                             
            LLT(1,1,I)=DCOS(H(I)*D1(I-1))
            LLT(1,2,I)=(-DSIN(H(I)*D1(I-1))*MM(I))/H(I)
            LLT(2,1,I)=DSIN(H(I)*D1(I-1))
            LLT(2,2,I)=(DCOS(H(I)*D1(I-1))*MM(I))/H(I) 
CCCCC                     
            RT(1,1,I)=DCOS(H(I)*(D1(I-1)-D1(I)))
            RT(1,2,I)=(DSIN(H(I)*(D1(I-1)-D1(I)))*MM(I))/H(I)
            RT(2,1,I)=-H(I)*DSIN(H(I)*(D1(I-1)-D1(I)))/MM(I)
            RT(2,2,I)=DCOS(H(I)*(D1(I-1)-D1(I))) 
CCCCC
            RRT(1,1,I)=DCOS(H(I)*D1(I))
            RRT(1,2,I)=(-DSIN(H(I)*D1(I))*MM(I))/H(I)
            RRT(2,1,I)=DSIN(H(I)*D1(I))
            RRT(2,2,I)=(DCOS(H(I)*D1(I))*MM(I))/H(I)              
        ELSE
            H(I)=DSQRT(-H(I))
            LT(1,1,I)=DCOSH(H(I)*(D1(I)-D1(I-1)))
            LT(1,2,I)=(DSINH(H(I)*(D1(I)-D1(I-1)))*MM(I))/H(I)
            LT(2,1,I)=H(I)*DSINH(H(I)*(D1(I)-D1(I-1)))/MM(I)
            LT(2,2,I)=DCOSH(H(I)*(D1(I)-D1(I-1))) 
CCCCC  
            LLT(1,1,I)=DCOSH(H(I)*D1(I-1))
            LLT(1,2,I)=(-DSINH(H(I)*D1(I-1))*MM(I))/H(I)
            LLT(2,1,I)=-DSINH(H(I)*D1(I-1))
            LLT(2,2,I)=(DCOSH(H(I)*D1(I-1))*MM(I))/H(I)    
CCCCC                             
            RRT(1,1,I)=DCOSH(H(I)*D1(I))
            RRT(1,2,I)=(-DSINH(H(I)*D1(I))*MM(I))/H(I)
            RRT(2,1,I)=-DSINH(H(I)*D1(I))
            RRT(2,2,I)=(DCOSH(H(I)*D1(I))*MM(I))/H(I) 
CCCCC              
            RT(1,1,I)=DCOSH(H(I)*(D1(I-1)-D1(I)))
            RT(1,2,I)=(DSINH(H(I)*(D1(I-1)-D1(I)))*MM(I))/H(I)
            RT(2,1,I)=H(I)*DSINH(H(I)*(D1(I-1)-D1(I)))/MM(I)
            RT(2,2,I)=DCOSH(H(I)*(D1(I-1)-D1(I)))
        ENDIF
      ENDDO
CCCCC      
      P(1,1)=LT(1,1,1)
      P(1,2)=LT(1,2,1)
      P(2,1)=LT(2,1,1)
      P(2,2)=LT(2,2,1)
CCCCC 
      Q(1,1)=RT(1,1,NL)
      Q(1,2)=RT(1,2,NL)
      Q(2,1)=RT(2,1,NL)
      Q(2,2)=RT(2,2,NL)
CCCCC  
      PT(1,1)=LLT(1,1,2)*P(1,1)+LLT(1,2,2)*P(2,1)
      PT(1,2)=LLT(1,1,2)*P(1,2)+LLT(1,2,2)*P(2,2)
      PT(2,1)=LLT(2,1,2)*P(1,1)+LLT(2,2,2)*P(2,1)
      PT(2,2)=LLT(2,1,2)*P(1,2)+LLT(2,2,2)*P(2,2)
      AX(2)=PT(1,1)*AX(1)+PT(1,2)*BX(NL)
      BX(2)=PT(2,1)*AX(1)+PT(2,2)*BX(NL) 
CCCCC 
      QT(1,1)=RRT(1,1,NL-1)*Q(1,1)+RRT(1,2,NL-1)*Q(2,1)
      QT(1,2)=RRT(1,1,NL-1)*Q(1,2)+RRT(1,2,NL-1)*Q(2,2)
      QT(2,1)=RRT(2,1,NL-1)*Q(1,1)+RRT(2,2,NL-1)*Q(2,1)
      QT(2,2)=RRT(2,1,NL-1)*Q(1,2)+RRT(2,2,NL-1)*Q(2,2)
      AX(NL-1)=QT(1,1)*AX(1)+QT(1,2)*BX(NL)
      BX(NL-1)=QT(2,1)*AX(1)+QT(2,2)*BX(NL)      
CCCCC           
      DO I=3,ICR
         PT(1,1)=LT(1,1,I-1)*P(1,1)+LT(1,2,I-1)*P(2,1)
         PT(2,1)=LT(2,1,I-1)*P(1,1)+LT(2,2,I-1)*P(2,1)
         PT(1,2)=LT(1,1,I-1)*P(1,2)+LT(1,2,I-1)*P(2,2)
         PT(2,2)=LT(2,1,I-1)*P(1,2)+LT(2,2,I-1)*P(2,2)
         P(1,1)=PT(1,1)
         P(2,1)=PT(2,1)
         P(1,2)=PT(1,2)                                                  
         P(2,2)=PT(2,2)   
         XP(1,1)=LLT(1,1,I)*P(1,1)+LLT(1,2,I)*P(2,1)
         XP(2,1)=LLT(2,1,I)*P(1,1)+LLT(2,2,I)*P(2,1)
         XP(1,2)=LLT(1,1,I)*P(1,2)+LLT(1,2,I)*P(2,2)
         XP(2,2)=LLT(2,1,I)*P(1,2)+LLT(2,2,I)*P(2,2)
         AX(I)=XP(1,1)*AX(1)+XP(1,2)*BX(NL)
         BX(I)=XP(2,1)*AX(1)+XP(2,2)*BX(NL)    
      ENDDO                  
CCCCC 
      DO I=NL-2,ICR,-1
         QT(1,1)=RT(1,1,I+1)*Q(1,1)+RT(1,2,I+1)*Q(2,1)
         QT(2,1)=RT(2,1,I+1)*Q(1,1)+RT(2,2,I+1)*Q(2,1)
         QT(1,2)=RT(1,1,I+1)*Q(1,2)+RT(1,2,I+1)*Q(2,2)
         QT(2,2)=RT(2,1,I+1)*Q(1,2)+RT(2,2,I+1)*Q(2,2)
         Q(1,1)=QT(1,1)
         Q(2,1)=QT(2,1)
         Q(1,2)=QT(1,2)
         Q(2,2)=QT(2,2)   
         YP(1,1)=RRT(1,1,I)*Q(1,1)+RRT(1,2,I)*Q(2,1)
         YP(2,1)=RRT(2,1,I)*Q(1,1)+RRT(2,2,I)*Q(2,1)
         YP(1,2)=RRT(1,1,I)*Q(1,2)+RRT(1,2,I)*Q(2,2)
         YP(2,2)=RRT(2,1,I)*Q(1,2)+RRT(2,2,I)*Q(2,2) 
         AX(I)=YP(1,1)*AX(1)+YP(1,2)*BX(NL)  
         BX(I)=YP(2,1)*AX(1)+YP(2,2)*BX(NL) 
      ENDDO     
CCCCC  
      WELL=Q(1,2)*P(2,1)-P(1,1)*Q(2,2) 
      REWIND MOUT 
      PRINT*, ' INPUT THE NAME OF OUTPUT FILE'
      READ(MIN,'(A40)') OUTPUT
      OPEN(UNIT=MOUT,FILE=OUTPUT, STATUS='UNKNOWN')
      XNORM=0.D0 
      XC(1)=((DABS(AX(1)))**2/(2*DABS(H(1))))
     +      *DEXP(-2.*(DABS(H(1)*D1(1))))
      XC(NL)=((DABS(BX(NL)))**2/(2*DABS(H(NL))))
     +       *DEXP(-2.*(DABS(H(NL)*D1(NL-1))))
      DO I=2,NL-1
         AXX=AX(I)
         BXX=BX(I)
         HH1X=HH1(I)
         HH2X=HH2(I)
         XL=D1(I)
         XR=D1(I-1)
         CALL QUANC8(CFM,xR,xL,0.d0,1.d-15,cc,errest,nofun,Xflag)
         XC(I)=CC 
         XNORM=XNORM+CC
      ENDDO
      XNORM=(XNORM+XC(1)+XC(NL))
      XNORN=DSQRT(DABS(XNORM))
      DO I=1,NL
      CONF(I)=XC(I)/XNORM
      AX(I)=AX(I)/XNORN
      BX(I)=BX(I)/XNORN
      WRITE(*,103) I,CONF(I)
  103 FORMAT(1X,'CONFINEMENT FACTOR OF ',I4,' th LAYER =',E15.8)
      ENDDO 
 4001 DO I=1,NL
         IF (I.EQ.1) THEN
         DY1=D1(1)-10.D0
         N1=ABS((D1(1)-DY1)/((D1(1)-DY1)/5.D0))+1
         DO J=1,N1
            DZ1(J)=DY1+(J-1)*((D1(1)-DY1)/5.D0)
            PHC(J)=AX(1)*DEXP(H(1)*DZ1(J))
            WRITE(MOUT,101) DZ1(J),PHC(J),V(I)
         ENDDO
         ELSEIF (I.EQ.NL) THEN
                DY2=D1(NL-1)+10.D0
                N3=ABS((DY2-D1(NL-1))/((DY2-D1(NL-1))/5.D0))+1
                DO J=1,N3
                   DZ2(J)=D1(NL-1)+(J-1)*((DY2-D1(NL-1))/5.D0)
                   PHC(J)=BX(NL)*DEXP(-H(NL)*DZ2(J))
                   WRITE(MOUT,101) DZ2(J),PHC(J),V(I)
                ENDDO
                PP=1/(2*H(NL))
         ELSE
                DY=D1(I-1)
                DYY=D1(I)
                N2=ABS((DYY-DY)/((DYY-DY)/10.D0))+1
                DO J=1,N2   
                   DZZ(J)=DY+(J-1)*((DYY-DY)/10.D0)
                IF (IZ.EQ.2) THEN
                   IF ((E-V(I)).GT.0.D0) THEN
                PHC(J)=AX(I)*DCOS(H(I)*DZZ(J))+BX(I)*DSIN(H(I)*DZZ(J))
                   WRITE(MOUT,101) DZZ(J),PHC(J),V(I)
                   ELSE                                    
                PHC(J)=AX(I)*DCOSH(H(I)*DZZ(J))+BX(I)*DSINH(H(I)*DZZ(J))
                   WRITE(MOUT,101) DZZ(J),PHC(J),V(I)
                   ENDIF
                ELSEIF (IZ.NE.2) THEN
                   IF ((V(I)-E).GT.0.D0) THEN
                PHC(J)=AX(I)*DCOS(H(I)*DZZ(J))+BX(I)*DSIN(H(I)*DZZ(J))
                   WRITE(MOUT,101) DZZ(J),PHC(J),V(I) 
                   ELSE                                    
                PHC(J)=AX(I)*DCOSH(H(I)*DZZ(J))+BX(I)*DSINH(H(I)*DZZ(J))
                   WRITE(MOUT,101) DZZ(J),PHC(J),V(I)
                   ENDIF
                ELSE
                ENDIF
                ENDDO
         ENDIF
      ENDDO
  101 FORMAT(2X,E20.10,2X,E20.10,2X,E20.10)
      CLOSE(MOUT)        
  102 RETURN
      END 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION CFM(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IZ
      COMMON /PARA1/ IZ
      COMMON /CFM1/ AXX,BXX,HH1X,HH2X
      IF (IZ.EQ.2) THEN
         IF (HH1X.GT.0.D0) THEN
             H11=DSQRT(HH1X)
             CFM=DABS((AXX*DCOS(H11*X)+BXX*DSIN(H11*X))**2)
         ELSEIF (HH1X.LT.0.D0) THEN
             H11=DSQRT(-HH1X)
             CFM=DABS((AXX*DCOSH(H11*X)+BXX*DSINH(H11*X))**2)
         ELSE
         ENDIF
      ELSEIF (IZ.NE.2) THEN
         IF (HH2X.GT.0.D0) THEN
            H22=DSQRT(HH2X)
            CFM=DABS((AXX*DCOS(H22*X)+BXX*DSIN(H22*X))**2)
         ELSEIF (HH2X.LT.0.D0) THEN
            H22=DSQRT(-HH2X)
            CFM=DABS((AXX*DCOSH(H22*X)+BXX*DSINH(H22*X))**2)
         ELSE
         ENDIF
      ELSE
      ENDIF
      RETURN
      END      

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE ZOUT(E,ERR)
      REAL*8 E,ERR
      INTEGER IZ
      COMMON /PARA1/ IZ
      OPEN(20,FILE='energy.dat',status='unknown')
      WRITE(*,500) E,ERR
      IF (IZ.EQ.2) THEN
      WRITE(20,497) E,ERR
      ELSEIF (IZ.EQ.3) THEN
      WRITE(20,498) E,ERR
      ELSEIF (IZ.EQ.4) THEN
      WRITE(20,499) E,ERR
      ELSE
      ENDIF
      WRITE(20,496)
  496 FORMAT(2X,'                                                  ')   
  497 FORMAT(2X,'CONDUCTION BAND ENERGY===>',E20.12,' ERROR= ',E12.7)
  498 FORMAT(2X,'HEAVY HOLE ENERGY===>',E20.12,' ERROR= ',E12.7)
  499 FORMAT(2X,'LIGHT HOLE ENERGY===>',E20.12,' ERROR= ',E12.7)       
  500 FORMAT(2X,'ENERGY EIGENVALUE===>',E20.12,' ERROR= ',E12.7)
      RETURN
      END       
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE VALENCE
      IMPLICIT REAL*8 (A-H,O-Z)
      PRINT*, ' THIS SUBROUTINE USES THE 2X2 PIKUS HAMILTONIAN '
      PRINT*, ' TO FIND BOTH HEAVY AND LIGHT HOLE ENERGY LEVEL '
      PRINT*, ' WITH BAND MIXING EFFECT.                       '
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE LASER
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(KX=500,KK=11)     
      REAL*8 DGN(KX),DGAN(KX),DGAJ(KX),POWF(KX),POWB(KX)
      REAL*8 POW(KX),JX(KX),NC(KX),JUX(KX),DGNX(KX)
      REAL*8 MO,MB,L,LX,N,XEE(KX),FR3(KX)
      REAL*8 AU(KX),JU(KX),JP(KX),MH,ML,CN(KX),DBY(KX),DB2(KX)
      REAL*8 RSP(KX),PP(KX),PQ(KX),QP(KX),DB3(KX)
      REAL*8 JKX(KX),FR2(KX),NC1
      REAL*8 EFCX(KX),EFVX(KX),DN(KX),JP2,JP3
      REAL*8 PFX(KX),ALX(KX),LAM,XEG(KX),FND(KX)
      REAL*8 H1(KX),H2(KX),GZ(KX),LZ,XXLEAK(KX),CPX(KX)
      REAL*8 WM(KX),NCX,X(KX),Y(KX),GAMMA(KX),SDP(KX),RIN(KX)
      REAL*8 PFXX(KX),XNPP(KX),XNC(KX)
      INTEGER LA,FYI,MM,TMODE,NDATA,MZ,NOPT,TEMP000
      CHARACTER OUTPUT*40,IFILE*40,GAINJ*40,GAINJA*40,GAINC*40,RINDEX*40
      CHARACTER DIFGAN*40,ANTIG*40,POWER*40,RELAX*40,XKK*40
      CHARACTER GAINCOM*40,POWER2*40,POWER3*40,YLEAK*40
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /ZM/ MO
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z7/ XEX
      COMMON /Z10/ SFH,SFL
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /Z17/ XX,XY,QY,XZ
      COMMON /Z18/ EC,EV
      COMMON /Z20/ FCX,FPX      
      COMMON /SPLIT/ XMS
      COMMON /Z22/ FF
      COMMON /Z23/ YQ
      COMMON /Z24/ CXZ,CXY,ECC,EVV
      COMMON /Z25/ XXEP
      COMMON /Z26/ NCX,LZ    
      COMMON /SQUARE/ X,Y
      COMMON /XMODE/ MODE,TMODE
      COMMON /WIDTH/ L
      COMMON /DIFFG/ QM,MH,ML
      COMMON /Z27/ MB,EO,N      
      COMMON /Z28/ TP,GR,DI,AAP,DI2,WOSC,XJTH,XJP
      COMMON /NOISE/ A1,A2,WRX,WDX,HIY,PWRF
      COMMON /Z29/ EC2,EH2,EL2
      COMMON /Z32/ PP,PQ,QP,CN,XNC,DN
      COMMON /Z33/ CONFINE
      COMMON /Z34/ EP
      EXTERNAL SSA,SSB,SSC,SSD,SSE,SSF,SSI,SSJ,SSK,SSL,SSM,SSN
      MIN=5         
      MOUT1=11
      MOUT2=12
      MOUT3=13
      MOUT4=14
      MOUT5=15
      MOUT6=16
      MOUT7=17
      MOUT8=18
      MOUT9=19
      MOUT10=20
      MOUT11=21
      MOUT12=22
      MOUT13=23        
      MOUT14=24
      MOUT15=25
      MOUT16=26  
      MOUT17=27   
      MOUT18=28 
      MOUT19=29
      MOUT20=30                   
      MOUT21=31
      MOUT22=32
      MOUT23=33
      MOUT24=34
      MOUT25=35     
      MOUT26=36
      PRINT*, ' THE INPUT FILE NAME='
      READ(MIN,'(A40)') IFILE
      OPEN(1,FILE=IFILE,status='unknown') ! INPUT DATA FILE   
      REWIND 1
      REWIND 11
      REWIND 12
      REWIND 13
      REWIND 14
      REWIND 15
      REWIND 16
      REWIND 17
      REWIND 18
      REWIND 19
      REWIND 20 
      REWIND 21
      REWIND 22
      REWIND 23
      REWIND 25       
      REWIND 26
      REWIND 27
      REWIND 28
      REWIND 29
      REWIND 30
      REWIND 31
      REWIND 32
      REWIND 33
      REWIND 34
      REWIND 35
      REWIND 36
      DO I=1,21
      READ(1,*)
      ENDDO
      READ(1,*) XX,XZ,QY,XY,LX,N,LAM
      READ(1,*) EG,TEMP,EC,EV
      DO I=1,8
      READ(1,*)
      ENDDO
      READ(1,*) EC1,EH1,EL1,EC2,EH2,EL2
      READ(1,*) ALPHA,R1,R2,MM,BETA
      DO I=1,8
      READ(1,*)
      ENDDO
      READ(1,*) CL,CW,ETHA,CA,ES,CONFINE
      READ(1,*) CXZ,CXY,ECC,EVV ! FOR CALCULATE THE LEAKAGE   
      TEMP000=INDEX(IFILE,' ')
      GAINJ=IFILE
      GAINJ(temp000:)='.gjr'
      GAINC=IFILE
      GAINC(temp000:)='.gc'
      RINDEX=IFILE
      RINDEX(temp000:)='.ric'
      GAINJA=IFILE
      GAINJA(temp000:)='.gja'
      DIFGAN=IFILE
      DIFGAN(temp000:)='.dg'
      ANTIG=IFILE
      ANTIG(TEMP000:)='.ant'
      POWER=IFILE
      POWER(TEMP000:)='.pi'
      RELAX=IFILE
      RELAX(TEMP000:)='.rq'
      XKK=IFILE
      XKK(TEMP000:)='.k'
      GAINCOM=IFILE
      GAINCOM(TEMP000:)='.cp'
      POWER2=IFILE
      POWER2(TEMP000:)='.pi2'
      POWER3=IFILE
      POWER3(TEMP000:)='.pi3'
      YLEAK=IFILE
      YLEAK(TEMP000:)='.lek'
      L=LX*1.D-9
      LZ=L
      Q=1.D0
      QX=1.6D-19
      CXL=CL
      K=400
      QQ=401
      T2=0.1D-12
      MO=5.692D-12
      XK=8.617D-5
      XKT=XK*TEMP
      PI=3.1415926535897932D0
      H=4.135D-15
      XH=(H/(2.D0*PI))
      LAM=LAM*1.D-4
CCCC  C,L,EO,LX USE MKS SYSTEM'
      C=2.99D8
      EO=5.5338D7
      PRINT*, ' SELECT MATERIAL=?' 
      PRINT*, ' 1--AlGaAs'
      PRINT*, ' 2--InGaAsP' 
      PRINT*, ' 3--In1-zGazAs/InGaAsP/InP'
      PRINT*, ' 4-- InGaAlAs'
      PRINT*, ' 5--GaInP/AlzGawIn1-z-wP/Al0.5In0.5P'
      PRINT*, ' 6-- InxGa1-xAs/AlxGa1-xAs/AlGaAs '
      PRINT*, ' 7--In1-xGaxAs/InGaAsP/GaxIn1-xP(X=0.51) MATCHED TO GaAs'
      PRINT*, ' 8--AlyInxGa1-x-yAs/AlzGa1-zAs/GaAs '  
      PRINT*, ' 9--InzGa1-zAs/AlxGayIn1-x-yAs/InP'
      PRINT*, ' 10-- InGaAlAs/InGaAlAs/AlAsSb'
      PRINT*, ' 11--InzGa1-zAs/AlxGayIn1-x-yAs/AlAsSb'
      PRINT*, ' INPUT SELECTION'
      READ(*,*) FYI
CCCC  DELTA IS THE SPIN-OFF SPLITTING FOR InGaAsP
      IF (FYI.EQ.1) THEN
         DELTA=XX*0.28+(1-XX)*0.34
         UN=0.03*XX+(1-XX)*0.85
         DES=N**2.D0 
      ELSEIF (FYI.EQ.2) THEN
             DELTA=0.11+0.31*QY-0.09*QY**2
             EPSILON=(1-XX)*QY*14.6+(1-XX)*(1-QY)*12.4+XX*QY*13.2
     +      +XX*(1-QY)*11.1
             DES=EPSILON      
      ELSEIF (FYI.EQ.3) THEN
             DELTA=0.34*XX+(1-XX)*0.38
             DES=N**2.D0
      ELSEIF (FYI.EQ.4) THEN
             DELTA=(1.-XX-QY)*0.38+QY*0.28+XX*0.34
             DES=N**2.D0
      ELSEIF (FYI.EQ.5 ) THEN
c            DELTA=XX*0.08+(1.D0-XX)*0.11
	     DELTA=0.11
             ZZ1=3.39+0.62*XZ
             ZZ2=28.07+1.72*XZ
             DES=N**2.D0
      ELSEIF (FYI.EQ.6 ) THEN
             DELTA=0.38*XX+(1-XX)*0.34
             DES=N**2
      ELSEIF (FYI.EQ.7) THEN
             DELTA=0.34*XX+(1-XX)*0.38
             DES=N**2
      ELSEIF (FYI.EQ.8) THEN
             DELTA=XX*0.28+QY*0.38+(1-XX-QY)*0.34
             DES=N**2
      ELSEIF (FYI.EQ.9) THEN
             DELTA=XX*0.38+(1-XX)*0.34
             DES=N**2
      ELSEIF (FYI.EQ.10) THEN
             DELTA=(1.-XX-QY)*0.38+QY*0.28+XX*0.34
             DES=N**2
      ELSEIF (FYI.EQ.11) THEN
             DELTA=XX*0.38+(1-XX)*0.34
             DES=N**2
      ELSE      
      ENDIF              
      PRINT*, 'INPUT MODE = ? FOR TE--> MODE =1, FOR TM--> MODE =2'
      PRINT*, ' INPUT TE OR TM ?'
      READ(*,*) MODE
      PRINT*, 'IF EL1 BELOW EH1 THEN SELECT 1, OTHERWISE SELECT 2' 
      PRINT*, ' SELECTION=?'
      READ(*,*) TMODE
      XES=(1.6D-19/(2*PI*8.854D-12*DES))*((3/PI)**(1.D0/3.D0))*1.D2
      EP=EG+EC+EV
      VG=(C/N)*1.D2
      WRITE(*,77)
      PRINT*, 'CALCULATE THE EFFECTIVE MASS'
      CALL EMASS (FYI)
      PXX=(1.D0/(2.D0*CL))*DLOG(1.D0/(0.09D0))
      THRGG=ALPHA+PXX
      PX=(1.D0/(2.D0*CL))*DLOG(1.D0/(R1*R2))
      THRG=ALPHA+PX  
      WRITE(*,77)
      PY=PX/THRG
      XS1H=1.D0*XM1CH
      XS1L=1.D0*XM1CL
      PRINT*, ' FOR QUASI-FERMI LEVEL SELECT=1, '
      PRINT*, ' FOR READ EXISTING QUASI-FERMI LEVEL SELECT=2 '
      PRINT*, ' SELECT=?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 7000
      ELSEIF (I.EQ.2) THEN
      GOTO 7002
      ELSE
      ENDIF
 7000 PRINT*, 'CALCULATE THE QUASI FERMIC LEVELS'
c     WRITE(*,77)
  700 CALL XSEFC (K,EFCX,FYI)
  701 CALL XSEFP (K,EFVX,FYI)
      GOTO 7001
 7002 OPEN(UNIT=MOUT12,FILE='cfermil.dat',status='UNKNOWN')
      OPEN(UNIT=MOUT13,FILE='vfermil.dat',status='UNKNOWN')
      DO I=1,K
         READ(MOUT12,40) NC(I),EFCX(I)
         READ(MOUT13,41) NC(I),EFVX(I)
      ENDDO
   40 FORMAT(12X,E20.12,14X,E20.12)
   41 FORMAT(12X,E20.12,16X,E20.12)
 7001 CONTINUE
      CMR=(XMC1*XMV1H)/(XMC1+XMV1H)
      AO=((XH**2*DES*4*PI*EO)/CMR)*1.D2 
      ER=((1/(4*PI*EO))**2*CMR)/(2*DES**2*XH**2) 
      WRITE(*,77)
      OPEN(UNIT=MOUT1,FILE=GAINJ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT2,FILE=GAINC ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT3,FILE=RINDEX ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT4,FILE=GAINJA ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT5,FILE=DIFGAN ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT6,FILE=ANTIG ,STATUS='UNKNOWN') 
      OPEN(UNIT=MOUT17,FILE=GAINCOM ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT26,FILE=YLEAK ,STATUS='UNKNOWN')
  801 DO I=1,K
         EFC=EFCX(I)
         EFV=EFVX(I)
         NC(I)=1.D17+(I-1)*((8.D18-1.D17)/(K-1))
         NCX=NC(I)
C       XEG(I)=-XES*(NC(I)**(1.D0/3.D0))  ! Theory 2 
         CALL DNPEF (L,DNEF,DPEF)       ! Theory 3
C         XEG(I)=-3.2*1.D-8*(NC(I)**(1/3)) ! Theory 4	(FOR AlGaInP)
         XKX=DSQRT((2*PI/EO)*(DNEF+DPEF))
         DCH=1+DSQRT((32*PI*NC(I)*L*1.D2)/(2*XKX**3*AO))
         XEG(I)=-2*ER*AO*XKX*DLOG(DCH)      
         XEP=EP+XEG(I)
         EGX=EG+XEG(I)
         XXEP=EGX+(I-1)*((XEP-EGX)/(K-1))  
CCCCC         
         IF (TMODE.EQ.1) THEN
            H1(I)=EGX+EC1+EH1      ! TE AND FOR EL1 BELOW EH1(TM)
            HO1=H1(I)                 ! EC1 AND EH1 TRANSITION
            H2(I)=EGX+EC1+EL1
            HO2=H2(I)           
            ELSEIF (TMODE.EQ.2) THEN  ! TM AND EL1 ABOVE EH1
                   H1(I)=EGX+EC1+EL1      ! EC1 AND EL1 TRANSITION
                   HO1=H1(I)
                   H2(I)=EGX+EC1+EH1
                   HO2=H2(I)
            ELSE
         ENDIF 
CCCCC         
       MB=(MO**2*EGX*(EGX+DELTA))/(6.D0*XMC1*(EGX+(2.D0/3.D0)*DELTA))
       GGG=Q**4*((DABS(MB))**2)*QX*16*(PI)**2*N
       GXX=(GGG/(EO**2*MO**4*C**4*H**4*XH*N*L))*1.D-6
       GG=((Q**2*(DABS(MB)))/(EO*C*MO**2*N*XH*L))*1.D-2
       GZ(I)=GG
       WS1=((16*(PI**2)*N*Q**2*(DABS(MB)))/(EO*MO**2*C**3*H**4*LX))
     +       *1.D-4
       WS=WS1*QX*LX
       WS3=(16*(PI*PI)*N*Q*Q*(DABS(MB)))/(EO*MO*MO*C**3*H**4*LX*1.D-7)
       WS2=WS3*1.D-4
       call quanc8(ssa,egx,Xep,0.d0,1.d-12,dd,errest,nofun,Xflag)
       call quanc8(ssb,egx,xep,0.d0,1.d-12,ff,errest,nofun,Xflag)  
       call quanc8(ssc,egx,xep,0.d0,1.d-5,yq,errest,nofun,Xflag) 
       call quanc8(sse,egx,xep,0.d0,1.d-12,zq,errest,nofun,Xflag) 
       call quanc8(ssf,egx,xep,0.d0,1.d-5,rp,errest,nofun,Xflag)    
       call quanc8(ssi,egx,Xep,0.d0,1.d-12,pf,errest,nofun,Xflag)  
       call quanc8(ssj,egx,xep,0.d0,1.d-12,pz,errest,nofun,Xflag)
       call quanc8(ssm,egx,xep,0.d0,1.d-20,cp,errest,nofun,Xflag)
C       dd-dG/dN, ff-J(RAD), yq-G(J)W/O L(E), zq-G(J)W/L(E), rp-RSP(E)
C      pf- d(dPhi/dZ)dN, pz-dPhi/dZ, cp-Gain compression  
      CALL LEAKAGE (FYI,XLEAK)
      WRITE(*,628) XLEAK,NC(I)
      WRITE(MOUT26,*) NC(I)*1.E-18,XLEAK
      XXLEAK(I)=XLEAK
  628 FORMAT(2X,'J(LEAKAGE)=',D12.6,' A/cm^2',2X,'N=',D12.6,' 1/cm^3')      
       CALL XXNP (XNP,LAM,N)     
       RSP(I)=RP 
       CPX(I)=CP
       PFXX(I)=PF
       XNPP(I)=XNP
       DGN(I)=DD
       ALX(I)=((PFXX(I))/DGN(I))
       JX(I)=FF
       JKX(I)=FF/NC(I)
       DGAN(I)=ZQ
       DGAJ(I)=ZQ
       XEX=XNP*NC(I)
       PFX(I)=((-PF/(4*PI/LAM))*NC(I)-XEX)*CONFINE
       CB=CA*(1.D0-33.D0*ES)   ! AlGaAs ONLY
       AU(I)=CB*(NC(I)**3.D0)        
       JUX(I)=(JX(I)+QX*L*1.D2*AU(I))+XLEAK
       JU(I)=(JX(I)+QX*L*1.D2*AU(I))*MM+XLEAK
       WRITE(MOUT1,*) JU(I),DGAJ(I)*MM*CONFINE
       WRITE(MOUT2,*) NC(I)*1.d-18,DGAN(I)
       WRITE(MOUT3,*) JU(I),PFX(I)
       WRITE(MOUT4,*) JU(I),DGAJ(I)
       WRITE(MOUT5,*) NC(I)*1.d-18,DGN(I)*1.d16    
       WRITE(MOUT6,*) NC(I)*1.D-18,ALX(I)       
       WRITE(MOUT17,*) NC(I)*1.D-18,CP
      ENDDO
 1003 WRITE(*,77)
      DO I=1,K
         IF (DGAN(I).GE.0.D0) THEN
            JJ=I
            GOTO 888
            ELSE
         ENDIF
      ENDDO
  888 CONTINUE
c     WRITE(*,77) 
      DO I=1,K   
         PP(I)=DGAJ(I)
         PQ(I)=JUX(I)
         QP(I)=DGAJ(I)/JUX(I)                    
         CN(I)=DGAN(I)/NC(I)
         DN(I)=DGAN(I)
         XNC(I)=NC(I)
      ENDDO
      CALL BUBLE (GO,XNGO,XJO,XNO)
      PRINT*, ' G(J) PARAMETERS FROM SINGLE WELL'
      WRITE(*,630) GO,XJO
  630 FORMAT(2X,'Go=',D12.6,' 1/cm',2X,'Jo=',D12.6,' A/cm^2',/)
      PRINT*, ' G(N) PARAMETERS FROM SINGLE WELL'
      WRITE(*,631) XNGO,XNO
  631 FORMAT(2X,'NGo=',D12.6,' 1/cm',2X,'XNo=',D12.6,' 1/cm^3',/)
      XJTR=XJO*DEXP(-1.d0)
      XNTR=XNO*DEXP(-1.d0)
      WRITE(*,301) XJTR,XNTR
  301 FORMAT(2X,'Jtr=',D12.6,' A/cm^2',2X,'NTR=',D12.6,' 1/cm^3',/)
      MK=INT(THRGG/(GO))+1
      PRINT*, ' THE OPTIMUM NUMBER OF QUANTUM WELL Nopt =',MK
      PRINT*, ' INPUT Nopt(CAN BE DIFFERENT FROM ABOVE CALCULATION)=?'  
      READ(*,*) NOPT
      PRINT*, ' NUMBER OF QUANTUM WELL(MAY OR MAY NOT BE Nopt)=?'
      READ(*,*) MZ
      PRINT*, '                                                       '
 1005 WRITE(*,77)                            
      DO I=1,K
         IF ((CONFINE*DGAJ(I)-THRG).GE.0.D0) THEN
            IY=I                                
            HIY=H1(IY)
            GOTO 1004
            ELSE
         ENDIF
      ENDDO                                                             
      PRINT*, ' NOT REACH THRESHOLD CONDITION, TRY ANOTHER POSSIBILITY' 
      GOTO 900     
 1004 XJTH=(JX(IY)+QX*L*1.D2*AU(IY))*NOPT+XXLEAK(IY) 
C      JMK=INT(THRG/GO)+1
      XJJTH=MM*XJO*DEXP((NOPT*1.D0/MM)-1)
c     XJJJTH=(MZ/ETHA)*XJO*DEXP((1.d0*NOPT/MZ)-1)
      XJJJTH=MZ*XJO*DEXP((1.d0*NOPT/MZ)-1)
      WRITE(*,77)
      PRINT*, ' 1ST CHECK USE SINGLE WELL TIMES # OF WELLS'
      WRITE(*,77)
      PRINT*, ' 2ND CHECK USE McIlory METHOD FOR SINGLE WELL'
      PRINT*, ' THEN WITH Nopt'
      WRITE(*,77)
      PRINT*, ' 3RD CHECK FOLLOWS McIlory METHOD'
      WRITE(*,77)
      WRITE(*,75) THRG,NC(IY),IY,XJTH,XJJTH,XJJJTH
   75 FORMAT(2X,'Gth=',F9.4,' 1/cm',' Nth=',D12.6,' 1/cm^3',' IY=',I5/
     +       2X,' 1ST CHECK Jth=',F15.8,' A/cm^2'/
     +       2X,' 2ND CHECK Jth=',F12.5,' A/cm^2'/
     +       2X,' 3RD CHECK Jth=',F12.5,' A/cm^2'/)
      XITH=XJTH*CXL*CW*1.D3
      XIITH=XJJTH*CXL*CW*1.D3
      XIIITH=XJJJTH*CXL*CW*1.D3
      WRITE(*,627) XITH,NOPT,XIITH,XIIITH
  627 FORMAT(2X,' 1ST CHECK Ith=',D12.6,' mA',' NUMBER OF WELLS=',I3/
     +       2X,' 2ND CHECK Ith=',D12.6,' mA'/
     +       2X,' 3RD CHECK Ith=',D12.6,' mA'/)
c      DO I=1,K
c         WRITE(MOUT3,*) JU(I),(-PFX(I)*CONFINE*NOPT)/((2*PI)/LAM)
c      ENDDO
      WRITE(*,77)
 1008 PRINT*, 'CALCULATE THE P-I RELATION'
      OPEN(UNIT=MOUT7,FILE=POWER ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT8,FILE=RELAX ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT24,FILE=POWER2 ,STATUS='UNKNOWN')
      OPEN(UNIT=MOUT25,FILE=POWER3 ,STATUS='UNKNOWN')  
      DEV=EFCX(IY)-EFVX(IY)+EG+XEG(IY) 
      IF (TMODE.EQ.1) THEN
         XNSP=1/(1-DEXP(-(DEV-H1(IY))/XKT))  
         ELSE
         XNSP=1/(1-DEXP(-(DEV-H2(IY))/XKT))
      ENDIF
c      TP=(1.D0/(VG*THRG))
c      ETHAR=RSP(IY)/(RSP(IY)+AU(IY))
      BETA=(QX/(XITH*TP*1.D-3))*(XNSP/(ETHA*ETHAR))
c      PHNUMS=CONFINE*BETA*RSP(IY)/(VG*DABS(THRG-CONFINE*DGAN(IY)))
c       AP=CPX(IY)*DGAN(IY)/((1+CPX(IY)*PHNUMS))
c       AG=DGN(IY)/(1+CPX(IY)*PHNUMS)
c      XB=RSP(IY)/(NC(IY)**2)
c       TAUN=1.D0/(2*XB*NC(IY)+3*CA*NC(IY)**2)
C       PRINT*, ' PHOTON DENSITY= ',PHNUMS,' cm-3'
C       PRINT*, ' NONLINEAR DIFFERENTIAL GAIN(dG/dN)=',AG
C       PRINT*, ' NONLINEAR DIFFERENTIAL GAIN(dG/dNp)=',AP
C       PRINT*, ' DIFFERENTIAL LIFETIME(TAU) =',TAUN,' SEC'
       PRINT*, '                                                       ' 
      DO I=IY,K     
      JP(I)=(XJTH+((XJTH/20.d0)*(I-IY)))
       POW(I)=PY*ETHA*((JP(I)-XJTH)*CXL*CW)*H1(IY)
       IF (POW(I).GE.0.D0) THEN
       FN=((1-R1)*DSQRT(R2))/((DSQRT(R1)+DSQRT(R2))*(1-DSQRT(R1*R2)))
       BA=((1-R2)*DSQRT(R1))/((DSQRT(R1)+DSQRT(R2))*(1-DSQRT(R1*R2)))
       POWF(I)=FN*POW(I)
       POWB(I)=BA*POW(I)
       DGNX(IY)=DGN(IY)*(1./(1+CPX(IY)*PHNUMS))
       FR2(I)=(1./(2.*PI))*DSQRT((VG*CONFINE*DGN(IY)
     +      *(JP(I)-XJTH)*ETHA*CXL*CW)/(QX*L*1.D2*CW*CXL))
       FR3(I)=(1./(2.*PI))*DSQRT((VG*CONFINE*DGNX(IY)
     +      *(JP(I)-XJTH)*ETHA*CXL*CW)/(QX*L*1.D2*CW*CXL))
cc      SOS2(I)=MK*ETHA*(1.D0/(VG*THRG*QX*L*1.D2))*(JP(I)-XJTH)
       WRITE(MOUT7,*) (JP(I))*CXL*CW*1.d3,POWF(I)*1.d3
       WRITE(MOUT8,*) JP(I)*CXL*CW*1.D3,FR2(I)
       ELSE
      ENDIF
      ENDDO   
      DO I=1,KX
         JP2=(XJJTH+((XJJTH/10)*(I-1)))
         JP3=(XJJJTH+((XJJJTH/10)*(I-1)))
         POW2=PY*ETHA*((JP2-XJJTH)*CXL*CW)*H1(IY)
         POW3=PY*ETHA*((JP3-XJJJTH)*CXL*CW)*H1(IY)
         FN=((1-R1)*DSQRT(R2))/((DSQRT(R1)+DSQRT(R2))*(1-DSQRT(R1*R2)))
         IF (POW2.GE.0.D0) THEN
         POWF2=FN*POW2
         ELSE
         POW2=0.D0
         ENDIF
         IF (POW3.GE.0.D0) THEN
         POWF3=FN*POW3
         ELSE
         POW3=0.D0
         ENDIF
         WRITE(MOUT24,*) JP2*CXL*CW*1.d3,POWF2*1.d3
         WRITE(MOUT25,*) JP3*CXL*CW*1.d3,POWF3*1.d3
      ENDDO
      CLOSE (UNIT=1)
       NDATA=K-IY+1
       PRINT*, 'NDATA=',NDATA
      WRITE(*,77)                                  
      PRINT*, ' CALCULATE THE SLOPE: mW/mA Y=A+BX'
      DO I=IY,K
         X(I-(IY-1))=JP(I)*CXL*CW*1.D3
         Y(I-(IY-1))=POWF(I)*1.D3
      ENDDO
c      call fit(ndata,sig,0,a,b,siga,sigb,chi2,q)
      CALL LEST (NDATA,A,B,1,1)
      WRITE(*,5001) A,B
 5001 FORMAT(2X,'CONSTANT A=',F12.7,2X,' SLOPE B=',F12.7,/)  
      WRITE(*,77)
 8885 PRINT*, ' INPUT POWER PO FOR THE LINEWIDTH, PO=0 FOR STOP'
      PRINT*, ' INPUT PO=    mW'
      READ(*,*) PO
      IF (PO.EQ.0.D0) THEN
         GOTO 8887
         ELSE
         GOTO 8886
      ENDIF
 8886 DEV=EFCX(IY)-EFVX(IY)+EG+XEG(IY) 
      IF (TMODE.EQ.1) THEN
         GNSP=THRG/(1-DEXP(-(DEV-H1(IY))/XKT))  
         ELSE
         GNSP=THRG/(1-DEXP(-(DEV-H2(IY))/XKT))
      ENDIF 
      DALX=ALX(IY)      
      DF=(VG**2*H1(IY)*1.6E-19*GNSP*PX*(1+DALX**2))/(8*PI*1.D-3*PO)
      PRINT*, ' THE LINEWIDTH DF AND ALPHA FACTOR AT THRESHOLD ARE'
      WRITE(*,8880) DF*1.d-6,DALX
 8880 FORMAT(2X,'DF=',D18.10,' MHz',2X,' ALPHA FACTOR=',D12.7)
      WRITE(*,77)
      GOTO 8885  
 8887 PRINT*, 'INPUT 1 FOR THE DYNAMIC CALCULATION. 2 FOR SKIP'
      PRINT*, 'INPUT = '
      READ(*,*) I    
      TP=(1.D0/(VG*THRG))
      TSX=(1.D0/(CB*NC(IY)**2+(RSP(IY)/NC(IY))))
      IF (I.EQ.1) THEN
         GOTO 8888
         ELSEIF (I.EQ.2) THEN
                GOTO 8889
                ELSE
      ENDIF
 8888 PRINT*, 'INPUT THE VALUE ABOVE THRESHOLD:-- 0 STOP'
      PRINT*, ' I='
      READ(*,*) I 
      IF (I.EQ.0) THEN
         GOTO 8889
         ELSE
      ENDIF
      WR=FR2(IY+I)*2*PI
      WRX=FR3(IY+I)*2*PI         
      XJP=JP(IY+I)
      PHNUMS=(XJP-XJTH)/((QX/ETHA)*VG*THRG)
      WD=(1.D0/(TP+(1.D0/(WR**2*TSX))))  
      G1=(2*(RSP(IY)/NC(IY))+3*CB*NC(IY)**2)
      GR=4*(PI**2)*TP*(1+(CONFINE*CPX(IY)*DGAN(IY)
     +    /(DGN(IY)*(1+CPX(IY)*PHNUMS))))*(WRX/(2*PI))**2+G1   
      WDX=(WRX**2)/GR                                    
      WOSC=WRX*DSQRT(1-((GR/(2.*WRX))**2))
      DI=(TP*CONFINE*ETHA)/(QX*L*1.D2) 
      DI2=(ETHA)/(wosc*NOPT*QX*L*1.D2)
      PWRF=POWF(IY+I)
      DDALX=ALX(IY+I)
      AAP=(CPX(IY)*DGAN(IY)/(DGN(IY)*(1+CPX(IY)*PHNUMS)))
      A1=2*(BETA*RSP(IY))*FN*PY*(L*1.D2*CXL*CW)*(1/TP)
     +   *(G1**2)+FN*PY*(WRX**4)*((((JP(IY+I)+XJTH)*CXL*CW)
     +   /((JP(IY+I)-XJTH)*CXL*CW))-1)
      A2=2*(BETA*RSP(IY))*FN*PY*(L*1.D2*CXL*CW)*(1/TP)
     +   -2.*FN*PY*(WRX**2)*CONFINE*AAP
C      WRITE(*,8881) PHNUMS,TP,TSX
C      WRITE(*,8882) WD,WR,(JP(IY+I)-XJTH)*1.D3*CW*CXL
C 8881 FORMAT(2X,'PHOTON DENSITY ABOVE THRESHOLD=',D15.10,' 1/cm3',/
C     +       2X,' PHOTON LIFETIME=',D15.10,' SEC',/
C     +       2X,' CARRIER LIFETIME=',D15.10,' SEC',/)
C 8882 FORMAT(2X,' DAMPING ANGULAR FREQ.=',D15.10,'Hz',/
C     +       2X,' RELAX ANGULAR OSCILLATION FREQ.=',D15.10,'Hz',/
C     +       2X,' CURRENT ABOVE THRESHOLD=',F15.8,'mA',/)       
      PRINT*, ' INPUT NAME FOR FREQ. RESPONSE(GAIN COMPRESS INCLUDED)'
      PRINT*, ' DATA OUTPUT -- .GFR'
      READ(MIN,'(A40)') OUTPUT
      OPEN(UNIT=MOUT14,FILE=OUTPUT, STATUS='UNKNOWN')         
      PRINT*, ' INPUT THE NAME FOR FREQUENCY MODULATION'
      PRINT*, ' DATA OUTPUT -- .GFM'
      READ(MIN,'(A40)') OUTPUT
      OPEN(UNIT=MOUT9,FILE=OUTPUT, STATUS='UNKNOWN')                
      PRINT*, 'INPUT THE NAME FOR RIN CHARACTERISTIC(dB/Hz)'
      PRINT*, ' DTAT OUTPUT -- .RIN'
      READ(MIN,'(A40)') OUTPUT
      OPEN(UNIT=MOUT19,FILE=OUTPUT, STATUS='UNKNOWN')
      PRINT*, 'INPUT THE NAME FOR RIN(dB)'
      PRINT*, ' DATA OUTPUT -- . RDB'
      READ(MIN,'(A40)') OUTPUT
      OPEN(UNIT=MOUT20,FILE=OUTPUT, STATUS='UNKNOWN')
C              
      PRINT*, 'INPUT THE NAME FOR FREQUENCY NOISE DENSITY'
      PRINT*, ' DATA OUTPUT -- . FND'
      READ(MIN,'(A40)') OUTPUT
      OPEN(UNIT=MOUT21,FILE=OUTPUT, STATUS='UNKNOWN')
C      
      FM1=(ALX(IY)*CPX(IY)*ETHA*CONFINE)/(QX*4*PI*L*1.D6*CXL*CW
     +     *(1+CPX(IY)*PHNUMS))              
      AP=CPX(IY)*DGAN(IY)/((1+CPX(IY)*PHNUMS)**2)                 
      DO I=1,K   
         WM(I)=5.D9+(I-1)*((100*1.D9-5.D9)/(K-1))
         WWM0=(wm(i))/(2*pi)
         WWM=(5.D9+I*((100*1.D9-5.D9)/(K-1)))/(2*pi)       
         FM2=FM1*(WM(I)/(CONFINE*VG*PHNUMS*AP)) 
       DBY(I)=1.D0/(1-2*((WM(I)/WRX)**2)+(WM(I)/WRX)**4+(WM(I)/WDX)**2)
         DB2(I)=20*DLOG10(DBY(I)**0.5) 
         DB3(I)=10*DLOG10((DBY(I)*(FM1**2+FM2**2))**0.5)
c          DB3(I)=(DBY(I)*(FM1**2+FM2**2))**0.5
         SDP(I)=H1(IY)*PWRF*(((A1+A2*WM(I)**2)/(WRX**4)) 
     +          *(DBY(I))+1)
         RIN(I)=10*DLOG10((2.*SDP(I)/(PWRF**2))*QX)
         call quanc8(ssn,WWM0,WWM,0.d0,1.d-5,RINX,errest,nofun,Xflag)
         RINDB=10*DLOG10((2.*RINX/(PWRF**2))*QX)        
         FND(I)=(1/(2*PI))*(2*(BETA*RSP(IY))*FN*PY
     +          *(L*1.D2*CXL*CW)*(1/TP))*(H1(IY)/(8*PI*PWRF))
     +          *(1+DDALX**2*DBY(I))*1.6D-19
         WRITE(MOUT14,*) (WM(I)*1.D-9)/(2*PI),DB2(I) 
         WRITE(MOUT9,*) (WM(I)*1.D-9)/(2*PI),DB3(I)
         WRITE(MOUT19,*) (WM(I)*1.D-9)/(2*PI),RIN(I)
         WRITE(MOUT20,*) (WM(I)*1.D-9)/(2*PI),RINDB
         WRITE(MOUT21,*) (WM(I)*1.D-9)/(2*PI),FND(I)*1.d-6
      ENDDO   
      GOTO 8888
 8889 CONTINUE  
      DO I=IY,K                                  
         GAMMA(I+1)=4*(PI**2)*TP*(1+(CONFINE*CPX(IY)*DGAN(IY)
     +             /(DGN(IY)*(1+CPX(IY)*PHNUMS))))*FR3(I+1)**2+G1
         X(I-(IY-1))=(FR3(I+1)*1.D-9)**2
         Y(I-(IY-1))=GAMMA(I+1)*1.D-9
c         WRITE(MOUT15,*) X(I-(IY-1)),Y(I-(IY-1))
      ENDDO   
c      call fit(ndata,sig,0,a,b,siga,sigb,chi2,q)
      CALL LEST (NDATA,A,B,1,1)
      FRMAX=2.D0*PI*DSQRT(2.D0)/B          
      WRITE(*,5002) B,FRMAX 
 5002 FORMAT(2X,'K-FACTOR=',F8.5,' nS',2X,'MAXIUM FREQ.=',f8.4,' GHz',/)     
CCCC
  900 WRITE(*,77)
      WRITE(*,3001)
 3001 FORMAT(2X,'INPUT 1 FOR CALCULATE THE GAIN(E) RELATION.'/) 
      WRITE(*,3002) 
 3002 FORMAT(2X,'INPUT 2 FOR CALCULATE THE LINEWIDTH ENHENCEMENT',/ 
     +       2X,'FACTOR AND PHOTON ENERGY RELATION',/) 
      WRITE(*,3004)
 3004 FORMAT(2X,'INPUT 3 FOR EXIT THE PROGRAM',/)
      PRINT*, 'THE INPUT # IS'
      READ(*,*) I
      IF (I.EQ.1) THEN 
      PRINT*, 'INPUT FERMILEVELS IN C-BAND, V-BAND, AND CARRIER DENSITY'
      READ(*,*) FC,FP,NC1
      CALL DNPEF2 (L,FC,FP,DNEF,DPEF)
      XKX=DSQRT((2*PI/EO)*(DNEF+DPEF))
      DCH=1+DSQRT((32*PI*NC1*L*1.D2)/(2*XKX**3*AO))
      XEG1=-2*ER*AO*XKX*DLOG(DCH) 
      GOTO 911            
      ELSEIF (I.EQ.2) THEN
      PRINT*, 'INPUT FERMILEVELS IN C-BAND, V-BAND, AND CARRIER DENSITY'
      READ(*,*) FC,FP,NC1
      CALL DNPEF2 (L,FC,FP,DNEF,DPEF)
      XKX=DSQRT((2*PI/EO)*(DNEF+DPEF))
      DCH=1+DSQRT((32*PI*NC1*L*1.D2)/(2*XKX**3*AO))
      XEG1=-2*ER*AO*XKX*DLOG(DCH)   
      GOTO 912
      ELSEIF (I.EQ.3) THEN
      GOTO 2001
      ELSE
      ENDIF
  911 PRINT*, 'CALCULATE THE CONVOLUTION GAIN(E) COEFFICIENT'
      WRITE(*,77)
      PRINT*, 'INPUT THE NAME FOR THE CONVOLUTION OPTICAL GAIN(LAMBDA)'
      READ(MIN,'(A40)') OUTPUT
C
      OPEN(UNIT=MOUT10,FILE=OUTPUT ,STATUS='UNKNOWN')
C
      WRITE(*,77)
      PRINT*, 'INPUT THE NAME FOR THE CONVOLUTION MODE GAIN(LAMBDA)'
      READ(MIN,'(A40)') OUTPUT
C
      OPEN(UNIT=MOUT11,FILE=OUTPUT ,STATUS='UNKNOWN')
C
      PRINT*, 'INPUT THE NAME FOR THE CONVOLUTION OPTICAL GAIN(E)'
      READ(MIN,'(A40)') OUTPUT
C
      OPEN(UNIT=MOUT22,FILE=OUTPUT ,STATUS='UNKNOWN')
      WRITE(*,77)
      PRINT*, 'INPUT THE NAME FOR THE CONVOLUTION MODE GAIN(E)'
      READ(MIN,'(A40)') OUTPUT
C
      OPEN(UNIT=MOUT23,FILE=OUTPUT ,STATUS='UNKNOWN')
      DO I=1,K
         XEE(I)=EG+(I-1)*((EP-EG)/(K-1))
         XEX=XEE(I)
         IF (TMODE.EQ.1) THEN                   
         HO1=EG+EC1+EH1+XEG1
         HO2=EG+EC1+EL1+XEG1
         ELSEIF (TMODE.EQ.2) THEN
         HO1=EG+EC1+EL1+XEG1
         HO2=EG+EC1+EH1+XEG1
         ENDIF
         XEP=EP+XEG1
         EGX=EG+XEG1              
c        XXEP=0.5+(I-1)*((XEP-0.5)/(K-1)) 
	 XXEP=EGX+(I-1)*((XEP-EGX)/(K-1))
         call quanc8(ssd,egx,xep,0.d0,1.d-5,zy,errest,nofun,Xflag)
         WRITE(MOUT10,*) 1.24/XXEP,ZY
         WRITE(MOUT11,*) 1.24/XXEP,CONFINE*ZY
         WRITE(MOUT22,*) XXEP,ZY
         WRITE(MOUT23,*) XXEP,CONFINE*ZY
      ENDDO
      WRITE(*,77)    
      GOTO 944
  912 PRINT*, 'CALCULATE THE ANTIGUIDING FACTOR AND ENERGY RELATION'
      WRITE(*,77)
      PRINT*, 'INPUT THE NAME FOR THE ANTIGUIDING FACTOR'
      READ(MIN,'(A40)') OUTPUT
C
      OPEN(UNIT=MOUT16,FILE=OUTPUT ,STATUS='UNKNOWN') 
      DO I=1,K
         XEE(I)=EG+(I-1)*((EP-EG)/(K-1))
         XEX=XEE(I)
         IF (TMODE.EQ.1) THEN
         HO1=EG+EC1+EH1+XEG1
         HO2=EG+EC1+EL1+XEG1
         ELSEIF (TMODE.EQ.2) THEN
         HO1=EG+EC1+EL1+XEG1
         HO2=EG+EC1+EH1+XEG1
         ENDIF
         XEP=EP+XEG1
         EGX=EG+XEG1              
         XXEP=0.5+(I-1)*((XEP-0.5)/(K-1))                   
         call quanc8(ssk,HO1,xep,0.d0,1.d-5,af1,errest,nofun,Xflag)
         call quanc8(ssl,HO1,xep,0.d0,1.d-5,af2,errest,nofun,Xflag)  
         ANF=AF1/AF2
         WRITE(MOUT16,*) XXEP,ANF
      ENDDO
      WRITE(*,77)   
C  944 PRINT*, 'INPUT THE CONTROL NUMBER L FOR L=0 STOP' 
  944 PRINT*, 'INPUT 1 FOR REPEAT THE G(E) CALCULATION'
      PRINT*, 'INPUT 2 FOR REPEAT THE ALPHA(E) CALCULATION'
      PRINT*, 'INPUT 3 FOR EXIT' 
      READ(*,*) LA
C      IF (LA.EQ.0) THEN
C         GOTO 2001
      IF (LA.EQ.1) THEN    
                GOTO 900       
         ELSEIF (LA.EQ.2) THEN
                GOTO 900
         ELSEIF (LA.EQ.3) THEN
         GOTO 2001
         ELSE
      ENDIF
C
C      PRINT*, 'INPUT 1 FOR REPEAT G(E) CALCULATION'  
C      PRINT*, 'INPUT 2 FOR REPEAT ANTIGUIDING FACTOR CALCULATION'
C      READ(*,*) I
C      IF (I.EQ.1) THEN
C         GOTO 900                                               
C         ELSEIF (I.EQ.2) THEN
C                 GOTO 900
C         ELSE
C      ENDIF
C      WRITE(*,77)
 2001 CONTINUE   
   77 FORMAT(2X,'**************************************************')
      CLOSE (MOUT10)
      CLOSE (MOUT11)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE BUBLE (GO,XNGO,XJO,XNO)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KX=500,K=400)
      LOGICAL FLAG
      REAL*8 PP(KX),PQ(KX),QP(KX),CN(KX),XNC(KX),DN(KX)
      COMMON /Z32/ PP,PQ,QP,CN,XNC,DN
      COMMON /Z33/ CONFINE
      DO I=K-1,1,-1     
         FLAG=.FALSE.
         DO J=1,I
            IF (QP(J).GE.QP(J+1)) THEN
               T=QP(J)
               TT=PP(J)     
               TTT=PQ(J)
               TTTT=XNC(J) 
               QP(J)=QP(J+1)
               PP(J)=PP(J+1)
               PQ(J)=PQ(J+1)
               XNC(J)=XNC(J+1) 
               QP(J+1)=T
               PP(J+1)=TT   
               PQ(J+1)=TTT  
               XNC(J+1)=TTTT
            ENDIF
        ENDDO 
        IF (.NOT. FLAG) GOTO 629
      ENDDO                         
  629 XJO=TTT
      XNO=TTTT
      GO=TT*CONFINE
      XNGO=TT
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE EMASS (FYI)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MO,LXME,LXMWC,LXMC1,MH,ML
      INTEGER FYI
      COMMON /ZM/ MO
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /CLADDING/ XXMC1,LXMC1,GXMC1,GXMH,GXML
      COMMON /Z17/ XX,XY,QY,XZ
      COMMON /BARRIER/ XME,XMH,XML
      COMMON /SPLIT/ XMS      
      COMMON /Z24/ CXZ,CXY,ECC,EVV
      COMMON /DIFFG/ QM,MH,ML
      IF (FYI.EQ.1) THEN
         XME=(0.067+0.083*XZ)*MO
         XMH=(0.62+0.14*XZ)*MO
         XML=(0.087+0.063*XZ)*MO
         XMWC=(0.067+0.083*XX)*MO
         XMWVH=(0.62+0.14*XX)*MO
         XMWVL=(0.087+0.063*XX)*MO
         GXMH1=(0.62+0.14*CXZ)*MO
         GXML1=(0.087+0.063*CXZ)*MO
         GXMH=(GXMH1*XMH)/(GXMH1+XMH)
         GXML=(GXML1*XML)/(GXML1+XML)
C 
CCC   THE L AND X BAND ALSO CALCULATED, FOR LEAKAGE CURRENT
CCC   WHICH MEANS THE EFFECTIVE MASS OF CLADDING LAYER
C
         LXME=(0.56+0.1*CXZ)*MO  !L BAND CLADDING EFFECTIVE MASS
         LXMWC=(0.56+0.1*XZ)*MO
         LXMC1=(LXMWC*LXME)/(LXMWC+LXME)
         XXME=(0.85-0.14*CXZ)*MO      !X BAND CLADDING EFFECTIVE MASS
         XXMWC=(0.85-0.14*XZ)*MO
         XXMC1=(XXMWC*XXME)/(XXMWC+XXME)
         GXME=(0.067+0.083*CXZ)*MO         !G BAND CLADDING EFFECTIVE MASS
         GXMWC=(0.067+0.083*XZ)*MO
         GXMC1=(GXME*GXMWC)/(GXME+GXMWC)
         QM=(0.067D0+0.083D0*XX)*MO
         MH=(0.62D0+0.14D0*XX)*MO
         ML=(0.087D0+0.063D0*XX)*MO  
         GOTO 11     
CCCCC      
      ELSEIF (FYI.EQ.2) THEN
                XME=(0.08-0.039*XY)*MO
      BH=(1-XZ)*XY*0.6+(1-XZ)*(1-XY)*0.85+XZ*XY*0.62+XZ*(1-XY)*0.79
      BL=(1-XZ)*XY*0.027+(1-XZ)*(1-XY)*0.089+XZ*XY*0.074+XZ*(1-XY)*0.14
                XMH=BH*MO
                XML=BL*MO
CCCC  XMWC,XMWVH,XMWVL FOR THE WELL
CCCC  WH AND WL THE PARAMETERS FOR InGaAsP
CCCC  AWH, AWL, AND AWC THE WELL REGION PARAMETERS FOR THE InGaAs
                XMWC=(0.08D0-0.039D0*QY)*MO
      WH=(1-XX)*QY*0.6+(1-XX)*(1-QY)*0.85+XX*QY*0.62+XX*(1-QY)*0.79
      WL=(1-XX)*QY*0.027+(1-XX)*(1-QY)*0.089+XX*QY*0.074+XX*(1-QY)*0.14
                XMWVH=WH*MO
                XMWVL=WL*MO       
                GXMC=0.077*MO
                GXMC1=(XME*GXMC)/(XME+GXMC)
                GXMH1=0.61*MO
                GXML1=0.12*MO
                GXMH=(GXMH1*XMH)/(GXMH1+XMH)
                GXML=(GXML1*XML)/(GXML1+XML)
      QM=(0.080D0-0.039D0*QY)*MO
      TP=(1.-XX)*QY*0.6+(1.-XX)*(1.-QY)*0.85+XX*QY*0.62+XX*(1-QY)*0.79
      MH=TP*MO
      TR=(1-XX)*QY*0.027+(1-XX)*(1-QY)*0.089+XX*QY*0.074+XX*(1-QY)*0.14
      ML=TR*MO
                GOTO 11
         ELSEIF (FYI.EQ.3) THEN
                XME=(0.08-0.039*XY)*MO
      BH=(1-XZ)*XY*0.6+(1-XZ)*(1-XY)*0.85+XZ*XY*0.62+XZ*(1-XY)*0.79
      BL=(1-XZ)*XY*0.027+(1-XZ)*(1-XY)*0.089+XZ*XY*0.074+XZ*(1-XY)*0.14
                XMH=BH*MO
                XML=BL*MO
                AWC=0.025*(1-XX)+0.071*XX-0.0163*XX*(1-XX)
                AWH=0.62*XX+(1-XX)*0.6
                AWL=0.074*XX+(1-XX)*0.027
                XMWC=AWC*MO
                XMWVH=AWH*MO
                XMWVL=AWL*MO
                GXMH1=0.61*MO
                GXML1=0.12*MO
		GXMC=0.077*MO
		GXMC1=(XME*GXMC)/(XME+GXMC)
                GXMH=(GXMH1*XMH)/(GXMH1+XMH)
                GXML=(GXML1*XML)/(GXML1+XML)
                QM=(0.025D0*(1-XX)+0.071D0*XX-0.0163D0*XX*(1-XX))*MO
                TP=0.62D0*XX+(1.D0-XX)*0.60D0
                TR=0.074D0*XX+(1.D0-XX)*0.027D0
                MH=TP*MO
                ML=TR*MO
                GOTO 11
CCCCC
         ELSEIF (FYI.EQ.4) THEN
                XME=((1.-XY-XZ)*0.0239+XY*0.15+XZ*0.067)*MO
                RB1=(1.-XY-XZ)*20.4+XY*3.45+XZ*6.85
                RB2=(1.-XY-XZ)*8.3+XY*0.68+XZ*2.1
                XMH=(MO)/(RB1-2.D0*RB2)
                XML=(MO)/(RB1+2.D0*RB2)
                XMWC=((1.-XX-QY)*0.0239+QY*0.15+XX*0.067)*MO
                RW1=(1.-QY-XX)*20.4+QY*3.45+XX*6.85
                RW2=(1.-QY-XX)*8.3+QY*0.68+XX*2.1
                XMWVH=(MO)/(RW1-2.D0*RW2)
                XMWVL=(MO)/(RW1+2.D0*RW2)  
                GXMC=(0.48*0.15+0.52*0.023)*MO
                GXMC1=(XME*GXMC)/(XME+GXMC)
                GXMH1=0.61*MO
                GXML1=0.12*MO
                GXMH=(GXMH1*XMH)/(GXMH1+XMH)
                GXML=(GXML1*XML)/(GXML1+XML)
                QM=((1.-XX-QY)*0.0239+QY*0.15+XX*0.067 )*MO
                R1=(1.-QY-XX)*20.4+QY*3.45+XX*6.85
                R2=(1.-QY-XX)*8.3+QY*0.68+XX*2.1
                MH=(MO)/(R1-2.D0*R2)
                ML=(MO)/(R1+2.D0*R2)
                GOTO 11
CCCCC
         ELSEIF (FYI.EQ.5) THEN      
                Z=0.5*XZ	 !Al
                W=0.5*XY	 !Ga
                CZ=0.5*CXZ !Al
                CW=0.5*CXY !Ga
                XME=(0.11+0.00915*XZ-0.0024*(XZ**2))*MO
                RB1=(1.-Z-W)*6.35+Z*3.47+W*4.2
                RB2=(1.-Z-W)*2.08+Z*0.06+W*0.98
C                XMH=MO/(RB1-2.D0*RB2)
C                XML=MO/(RB1+2.D0*RB2)
C                XMH=(0.62+0.05*XZ)*MO
C                XML=(0.11+0.03*XZ)*MO
                XMH=(Z*0.62+0.5*0.45+(1-Z-0.5)*0.54)*MO
                XML=(Z*0.22+0.5*0.12+(1-Z-0.5)*0.16)*MO
C	        XMWC=(0.15*XX+(1-XX)*0.077)*MO
                XMWC=0.11*MO
                RW1=XX*4.2D0+(1.D0-XX)*6.35D0
                RW2=XX*0.98D0+(1.D0-XX)*2.08D0
C                XMWVH=MO/(RW1-2.D0*RW2)
C                XMWVL=MO/(RW1+2.D0*RW2)
                XMWVH=0.62*MO
                XMWVL=0.11*MO
C                XMWVH=(0.54*XX+(1.D0-XX)*0.45)*MO
C	          XMWVL=(0.16*XX+(1.D0-XX)*0.12)*MO
                RC1=(1.-CZ-CW)*6.35+CZ*3.47+CW*4.2
                RC2=(1.-CZ-CW)*2.08+CZ*0.06+CW*0.98
C                GXMH1=MO/(RC1-2.D0*RC2)
C                GXML1=MO/(RC1+2.D0*RC2)
C                GXMH1=(CZ*0.62+0.5*0.45+(1-CZ-0.5)*0.54)*MO
C	           GXML1=(CZ*0.22+0.5*0.12+(1-CZ-0.5)*0.16)*MO
                GXMH1=0.62+0.05*CXZ
                GXML1=0.11+0.03*CXZ
                GXMC=(0.11+0.00915*CXZ-0.0024*(CXZ**2))*MO
                GXMC1=(XME*GXMC)/(XME+GXMC)
                GXMH=(GXMH1*XMH)/(GXMH1+XMH)
                GXML=(GXML1*XML)/(GXML1+XML)
                QM=0.11*MO
                R1=XX*4.2D0+(1.D0-XX)*6.35D0
                R2=XX*0.98D0+(1.D0-XX)*2.08D0
                MH=(MO)/(R1-2.D0*R2)
                ML=(MO)/(R1+2.D0*R2)
                GOTO 11
CCCCCC
         ELSEIF (FYI.EQ.6) THEN
                XME=(0.067+0.083*XY)*MO
                RB1=6.85*(1-XY)+3.45*XY
                RB2=2.1*(1-XY)+0.68*XY
                CRB1=6.85*(1-CXZ)+3.45*CXZ
                CRB2=2.1*(1-CXZ)+0.68*CXZ
                XMH=MO/(RB1-2.D0*RB2)
                XML=MO/(RB1+2.D0*RB2)
                XMWC=(0.067-0.04*XX)*MO
                RW1=6.85*(1-XX)+19.67*XX
                RW2=2.1*(1-XX)+8.37*XX
                XMWVH=MO/(RW1-2.D0*RW2)
                XMWVL=MO/(RW1+2.D0*RW2)
      GXMC1=(XME*((0.067+0.083*CXZ)*MO))/(XME+((0.067*0.083*CXZ)*MO))
                GXMH1=MO/(CRB1-2.D0*CRB2)
                GXMH=(GXMH1*XMH)/(GXMH1+XMH)
                GXML1=MO/(CRB1+2.D0*CRB2)
                GXML=(GXML1*XML)/(GXML1+XML)
                QM=(0.067-0.04*XX)*MO
                R1=6.85*(1.D0-XX)+19.67*XX
                R2=2.1*(1.D0-XX)+8.37*XX
                MH=(MO)/(R1-2.D0*R2)
                ML=(MO)/(R1+2.D0*R2)
                GOTO 11
CCCCC
         ELSEIF (FYI.EQ.7) THEN
                XME=((1-XZ)*XY*0.023+(1-XZ)*(1-XY)*0.077+XZ*XY*0.067
     +          +XZ*(1-XY)*0.17)*MO
                RB1=(1-XZ)*XY*19.67+(1-XZ)*(1-XY)*6.35+XZ*XY*6.85
     +          +XZ*(1-XY)*4.2
                RB2=(1-XZ)*XY*8.37+(1-XZ)*(1-XY)*2.08+XZ*XY*2.1
     +          +XZ*(1-XY)*0.98
                XMH=MO/(RB1-2.D0*RB2)
                XML=MO/(RB1+2.D0*RB2)
                XMWC=(0.067-(1-XX)*0.04)*MO
                RW1=6.85*XX+19.67*(1-XX)
                RW2=2.1*XX+8.37*(1-XX)
                XMWVH=MO/(RW1-2.D0*RW2)
                XMWVL=MO/(RW1+2.D0*RW2)
                QM=(0.067-0.04*(1-XX))*MO
                R1=6.85*XX+(1.D0-XX)*19.67
                R2=2.1*XX+(1.D0-XX)*8.37
                MH=(MO)/(R1-2.D0*R2)
                ML=(MO)/(R1+2.D0*R2)
                GOTO 11
CCCCC
         ELSEIF (FYI.EQ.8) THEN
                XME=(0.067+0.083*XZ)*MO
                XMH=(0.62+0.14*XZ)*MO
                XML=(0.087+0.063*XZ)*MO
                XMWC=((1.-XX-QY)*0.067+XX*0.15+QY*0.0239)*MO
                RW1=XX*3.45+QY*20.40+(1-XX-QY)*6.85
                RW2=XX*0.68+QY*8.37+(1-XX-QY)*2.1
                CRB1=6.85*(1-CXZ)+3.45*CXZ
                CRB2=2.1*(1-CXZ)+0.68*CXZ
                XMWVH=MO/(RW1-2.D0*RW2)
                XMWVL=MO/(RW1+2.D0*RW2)
      GXMC1=(XME*((0.067+0.083*CXZ)*MO))/(XME+((0.067*0.083*CXZ)*MO))
                GXMH1=MO/(CRB1-2.D0*CRB2)
                GXMH=(GXMH1*XMH)/(GXMH1+XMH)
                GXML1=MO/(CRB1+2.D0*CRB2)
                GXML=(GXML1*XML)/(GXML1+XML)
                QM=(XX*0.15+QY*0.0239+(1-XX-QY)*0.067)*MO
                R1=XX*3.45+QY*20.4+(1-XX-QY)*6.85
                R2=XX*0.68+QY*8.3+(1-XX-QY)*2.1
                MH=(MO)/(R1-2.D0*R2)
                ML=(MO)/(R1+2.D0*R2)
                GOTO 11
CCCCC
         ELSEIF (FYI.EQ.9) THEN
                XME=((1.-QY-XY)*0.0239+QY*0.15+XY*0.067)*MO
                RB1=(1.-QY-XY)*20.4+QY*3.45+XY*6.85
                RB2=(1.-QY-XY)*8.3+QY*0.68+XY*2.1
                XMH=(MO)/(RB1-2.D0*RB2)
                XML=(MO)/(RB1+2.D0*RB2)
C                XMWC=(XX*0.0239+(1-XX)*0.067)*MO
                XMWC=(0.025*XX+0.071*(1-XX)-0.0163*(1-XX)*XX)*MO
                RW1=XX*20.4+(1-XX)*6.85
                RW2=XX*8.3+(1-XX)*2.1
                XMWVH=(MO)/(RW1-2.D0*RW2)
                XMWVL=(MO)/(RW1+2.D0*RW2)
                GXME=(CXZ*0.15+(1-CXZ)*0.023)*MO
                GXMC1=(XME*GXME)/(XME+GXME)
                QM=(0.025*XX+0.071*(1-XX)-0.0163*(1-XX)*XX)*MO
                R1=XX*20.4+(1-XX)*6.85
                R2=XX*8.3+(1-XX)*2.1
                MH=(MO)/(R1-2.D0*R2)
                ML=(MO)/(R1+2.D0*R2)
                GOTO 11
         ELSEIF (FYI.EQ.10) THEN
                XME=((1.-XY-XZ)*0.0239+XY*0.15+XZ*0.067)*MO
                RB1=(1.-XY-XZ)*20.4+XY*3.45+XZ*6.85
                RB2=(1.-XY-XZ)*8.3+XY*0.68+XZ*2.1
                XMH=(MO)/(RB1-2.D0*RB2)
                XML=(MO)/(RB1+2.D0*RB2)
                XMWC=((1.-XX-QY)*0.0239+QY*0.15+XX*0.067)*MO
                RW1=(1.-QY-XX)*20.4+QY*3.45+XX*6.85
                RW2=(1.-QY-XX)*8.3+QY*0.68+XX*2.1
                XMWVH=(MO)/(RW1-2.D0*RW2)
                XMWVL=(MO)/(RW1+2.D0*RW2)  
c                GXMC=(0.48*0.15+0.52*0.023)*MO
                GXMC=(0.56*0.15+0.44*0.92)*MO
                GXMC1=(XME*GXMC)/(XME+GXMC)
                GXMH1=0.94*MO
                GXML1=0.14*MO
                GXMH=(GXMH1*XMH)/(GXMH1+XMH)
                GXML=(GXML1*XML)/(GXML1+XML)
                QM=((1.-XX-QY)*0.0239+QY*0.15+XX*0.067 )*MO
                R1=(1.-QY-XX)*20.4+QY*3.45+XX*6.85
                R2=(1.-QY-XX)*8.3+QY*0.68+XX*2.1
                MH=(MO)/(R1-2.D0*R2)
                ML=(MO)/(R1+2.D0*R2)
                GOTO 11
         ELSEIF (FYI.EQ.11) THEN
                XME=((1.-QY-XY)*0.0239+QY*0.15+XY*0.067)*MO
                RB1=(1.-QY-XY)*20.4+QY*3.45+XY*6.85
                RB2=(1.-QY-XY)*8.3+QY*0.68+XY*2.1
                XMH=(MO)/(RB1-2.D0*RB2)
                XML=(MO)/(RB1+2.D0*RB2)
C                XMWC=(XX*0.0239+(1-XX)*0.067)*MO
                XMWC=(0.025*XX+0.071*(1-XX)-0.0163*(1-XX)*XX)*MO
                RW1=XX*20.4+(1-XX)*6.85
                RW2=XX*8.3+(1-XX)*2.1
                XMWVH=(MO)/(RW1-2.D0*RW2)
                XMWVL=(MO)/(RW1+2.D0*RW2)
C                GXME=(CXZ*0.15+(1-CXZ)*0.023)*MO
                GXME=(0.56*0.15+0.44*0.92)*MO
                GXMC1=(XME*GXME)/(XME+GXME)
                GXMH1=0.94*MO
                GXML1=0.14*MO
                GXMH=(GXMH1*XMH)/(GXMH1+XMH)
                GXML=(GXML1*XML)/(GXML1+XML)               
                QM=(0.025*XX+0.071*(1-XX)-0.0163*(1-XX)*XX)*MO
                R1=XX*20.4+(1-XX)*6.85
                R2=XX*8.3+(1-XX)*2.1
                MH=(MO)/(R1-2.D0*R2)
                ML=(MO)/(R1+2.D0*R2)                
                GOTO 11 
      ENDIF   
   11 XMC1=(XMWC*XME)/(XMWC+XME)
      XMV1H=(XMWVH*XMH)/(XMWVH+XMH)
      XMV1L=(XMWVL*XML)/(XMWVL+XML)
      XM1CH=(XMC1*XMV1H)/(XMC1+XMV1H)
      XM1CL=(XMC1*XMV1L)/(XMC1+XMV1L)  
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE LEAKAGE (FYI,NR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NR,LXMC1,GXMC1,GXMH,GXML,L,LZ,NCX 
      INTEGER FYI
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /BARRIER/ XME,XMH,XML
      COMMON /CLADDING/ XXMC1,LXMC1,GXMC1,GXMH,GXML
      COMMON /Z4/ EFC,EFV,FC,FP  
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /Z18/ EC,EV
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /Z22/ FF
      COMMON /Z24/ CXZ,CXY,ECC,EVV
      COMMON /WIDTH/ L
      COMMON /Z26/ NCX,LZ
      Q=1.6D-19         
      IF (FYI.EQ.1) THEN
         UN=CXZ*300+(1-CXZ)*8500 
         UP=CXZ*0.D0+(1-CXZ)*400 
         TB=5.D-9
         TN=2.D-9 
         ELSEIF (FYI.EQ.2) THEN
         TB=5.D-9
         TN=1.3D-8 
                UN=120.d0  
                UP=50.d0
         ELSEIF (FYI.EQ.3) THEN
                UN=120.d0
                UP=50.d0   
		TB=5.D-9
		TN=1.3D-8
         ELSEIF (FYI.EQ.4) THEN
c                UN=CXZ*300.D0+(1-CXZ)*22600.D0
c                UP=CXZ*0.D0+(1-CXZ)*200.D0
                UN=120.D0
                UP=50.D0 
                TB=5.D-9
                TN=2.D-9 
Cc     UN=CXZ*0.03+(1-CXZ)*3.2  
cc     UN=4000.d0
cc     UP=150.d0    
         ELSEIF (FYI.EQ.5) THEN
                UN=100.d0
                up=10.d0  
                TB=5.D-9
                TN=2.D-9  
         ELSEIF (FYI.EQ.6) THEN
		TB=5.E-9
		TN=2.E-9
                UN=CXZ*300+(1-CXZ)*8500.d0 
                UP=CXZ*0.D0+(1-CXZ)*400.d0  
         ELSEIF (FYI.EQ.8) THEN
                UN=CXZ*300.d0+(1-CXZ)*8500.d0
                UP=CXZ*0.D0+(1-CXZ)*400.d0 
		TB=5.E-9
		TN=1.3E-8
         ELSEIF (FYI.EQ.9) THEN
c                UN=CXZ*300.D0+(1-CXZ)*22600.D0
c                UP=CXZ*0.D0+(1-CXZ)*200.D0
         TB=5.D-9
         TN=1.3D-8 
                UN=120.d0  
                UP=50.d0
cc     UN=4000.d0
cc     UP=150.d0 
         ELSEIF (FYI.EQ.10) THEN
c                UN=CXZ*300+(1-CXZ)*200
c                UP=CXZ*0.D0+(1-CXZ)*400
         TB=5.D-9
         TN=1.3D-8 
                UN=120.d0  
                UP=50.d0
         ELSEIF (FYI.EQ.11) THEN
         TB=5.D-9
         TN=1.3D-8 
c                UN=120.d0  
c                UP=50.d0
                UN=CXZ*300+(1-CXZ)*200
                UP=CXZ*0.D0+(1-CXZ)*400
         ELSE
      ENDIF
      DN=UN*XKT  
      DP=UP*XKT
      XLN=DSQRT(DN*TN) 
c      XXLN=(1.D0/XLN)
      XLP=DSQRT(DP*TB)
      H=XH*2*PI
      DP=5.E17
      XJ=FF/(Q*DP*UP)
      XZ=XKT/XJ
      YNN=DEXP((EFC-ECC)/XKT)-(DEXP(2*(EFC-ECC)/XKT)/(2**1.5))
     +    +(DEXP(3*(EFC-ECC)/XKT)/(3**1.5))
      YPP=DEXP((EVV-EFV)/XKT) 
      U=DSQRT((1/(XLN**2))+(1/(4*XZ**2)))*1.E-4
      EU=(DEXP(U)-DEXP(-U))/(DEXP(U)+DEXP(-U))
      TU=DSQRT((1/(XLN**2))+(1/(4*XZ**2)))*EU+(1/(2*XZ))
      NR=((Q*2.D0*(GXMC1*XKT/(2*PI*XH**2))**1.5)*YNN*1.D-6)*DN*TU
      IF (FYI.EQ.2) THEN
      NR=NR*1.D-2
      ELSEIF (FYI.EQ.3) THEN
      NR=NR*1.D-2
      ELSE 
      ENDIF
      RETURN
      END

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSA(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 L,MH,ML
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z12/ EC1 
      COMMON /Z14/ EH1,EL1 
      COMMON /Z25/ XXEP
      COMMON /XMODE/ MODE,TMODE
      COMMON /DIFFG/ QM,MH,ML
      COMMON /WIDTH/ L 
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L                                 
CCCC  XXEA=Ec XXEB=Eh XXEC=El    
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFH=DEXP((XXEB-EFV)/XKT)*(1.D0/(1+DEXP((XXEB-EFV)/XKT)))**2
      XFL=DEXP((XXEC-EFV)/XKT)*(1.D0/(1+DEXP((XXEC-EFV)/XKT)))**2
      XFC=DEXP((XXEA-EFC)/XKT)*(1.D0/(1+DEXP((XXEA-EFC)/XKT)))**2 
      XFD=DEXP((XXED-EFC)/XKT)*(1.D0/(1+DEXP((XXED-EFC)/XKT)))**2  
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)  
      XFH=DEXP((XXEB-EFV)/XKT)*(1.D0/(1+DEXP((XXEB-EFV)/XKT)))**2
      XFL=DEXP((XXEC-EFV)/XKT)*(1.D0/(1+DEXP((XXEC-EFV)/XKT)))**2
      XFC=DEXP((XXEA-EFC)/XKT)*(1.D0/(1+DEXP((XXEA-EFC)/XKT)))**2
      XFD=DEXP((XXED-EFC)/XKT)*(1.D0/(1+DEXP((XXED-EFC)/XKT)))**2  
      ELSE
      ENDIF
c      DFF=(1+DEXP((EFC-EC1)/XKT))/DEXP((EFC-EC1)/XKT)
      DFF=1.D0/((DEXP((EFC-EC1)/XKT)/(1+DEXP((EFC-EC1)/XKT)))
     +   +(DEXP((EFC-EC2)/XKT)/(1+DEXP((EFC-EC2)/XKT))))		
      DN=((PI*XH**2*L)/QM)*1.D6*DFF
C      DPH=(DEXP((EFV-EH1)/XKT)+1)*(DEXP((EFV-EL1)/XKT)+1)
C      DPL=(DEXP((EFV-EH1)/XKT))*(DEXP((EFV-EL1)/XKT)+1)+((ML/MH)
C     +    *(DEXP((EFV-EL1)/XKT))*(DEXP((EFV-EH1)/XKT)+1))       
C      DP=(((PI*XH**2*L)/MH)*1.D6)*(DPH/DPL)
      DPHH=1.D0/((DEXP((EFV-EH1)/XKT)/(1+DEXP((EFV-EH1)/XKT)))
     +   +(DEXP((EFV-EH2)/XKT)/(1+DEXP((EFV-EH2)/XKT))))
      DPLL=1.D0/((DEXP((EFV-EL1)/XKT)/(1+DEXP((EFV-EL1)/XKT)))
     +   +(DEXP((EFV-EL2)/XKT)/(1+DEXP((EFV-EL2)/XKT))))
      DPH=((PI*XH**2*L)/MH)*1.D6*DPHH
	DPL=((PI*XH**2*L)/ML)*1.D6*DPLL
c	DP=DPH+DPL
      IF ((MODE.EQ.1).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/4.D0)*(1.D0+(HO1/X))
      A1L=(1.D0/4.D0)*(5.D0-3.D0*(HO2/X))
      WXS1H=XS1H*A1H*((1/XKT)*(XFC*DN+XFH*DPH))
      WXS1L=XS1L*A1L*((1/XKT)*(XFD*DN+XFL*DPL))
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO1/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO2/X)))
      WXS1H=XS1H*A1H*((1/XKT)*(XFC*DN+XFH*DPH))
      WXS1L=XS1L*A1L*((1/XKT)*(XFD*DN+XFL*DPL))
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.2)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO2/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO1/X)))
      WXS1H=XS1L*A1L*((1/XKT)*(XFC*DN+XFL*DPL))
      WXS1L=XS1H*A1H*((1/XKT)*(XFD*DN+XFH*DPH))           
      ELSE
      ENDIF                      
      SP=XH/T2
      IF (X.LT.HO1) THEN
      XF=0.D0
      ELSE
      XF=WXS1H
      AX=((SP)/((X-X2)**2+SP**2))*(1.d0/pi)
      IF (X.GE.HO2) THEN
      XF=WXS1H+WXS1L  
      AX=((SP)/((X-X3)**2+SP**2))*(1.d0/pi) 
      ELSE
      ENDIF
      ENDIF
      SSA=(GG/X)*XF*AX
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSB(X)           ! SPONTANEOUS EMISSION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP  
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /XMODE/ MODE,TMODE
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L) 
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
CCCC  XXEA=Ec XXEB=Eh XXEC=El 
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))  
      ZFCH=XFC*XFH
      ZFCL=XFD*XFL
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))  
      ZFCH=XFD*XFH
      ZFCL=XFC*XFL
      ELSE
      ENDIF
      IF (TMODE.EQ.1) THEN
      ZXS1H=XS1H*ZFCH
      ZXS1L=XS1L*ZFCL     
      ELSEIF (TMODE.EQ.2) THEN
      ZXS1H=XS1L*ZFCL
      ZXS1L=XS1H*ZFCH
      ELSE
      ENDIF
      IF (X.LT.HO1) THEN
      XF=0.D0
      ELSE
      XF=ZXS1H         
C      SSB=WS*X2*XF
      IF (X.GE.HO2) THEN
      XF=ZXS1H+ZXS1L
C      SSB=WS*X3*XF
      ELSE
      ENDIF
      ENDIF
      SSB=WS*X*XF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSC(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /Z22/ FF 
      COMMON /XMODE/ MODE,TMODE
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
CCCC  XXEA=Ec XXEB=Eh XXEC=El    
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))  
      ZFCH=XFC*XFH
      ZFCL=XFD*XFL
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))  
      ZFCH=XFD*XFH
      ZFCL=XFC*XFL
      ELSE
      ENDIF
      IF (TMODE.EQ.1) THEN
      ZXS1H=XS1H*ZFCH
      ZXS1L=XS1L*ZFCL     
      ELSEIF (TMODE.EQ.2) THEN
      ZXS1H=XS1L*ZFCL
      ZXS1L=XS1H*ZFCH
      ELSE
      ENDIF
      IF (X.LT.HO1) THEN
      XF=0.D0 
      ELSE
      XF=ZXS1H
      IF (X.GE.HO2) THEN
      XF=ZXS1H+ZXS1L
      ELSE
      ENDIF
      ENDIF
      SSC=(GXX/FF)*XF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSE(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /Z23/ YQ 
      COMMON /Z34/ EP
      COMMON /XMODE/ MODE,TMODE
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
CCCC  XXEA=Ec XXEB=Eh XXEC=El
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))  
      ZXFCH=-1.D0+XFC+XFH
      ZXFCL=-1.D0+XFD+XFL
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))  
      ZXFCH=-1.D0+XFD+XFH
      ZXFCL=-1.D0+XFC+XFL
      ELSE
      ENDIF
      IF ((MODE.EQ.1).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/4.D0)*(1.D0+(HO1/X))
      A1L=(1.D0/4.D0)*(5.D0-3.D0*(HO2/X))   
      WXS1H=XS1H*A1H*ZXFCH
      WXS1L=XS1L*A1L*ZXFCL
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO1/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO2/X)))
      WXS1H=XS1H*A1H*ZXFCH
      WXS1L=XS1L*A1L*ZXFCL
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.2)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO2/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO1/X)))
      WXS1H=XS1L*A1L*ZXFCL
      WXS1L=XS1H*A1H*ZXFCH
      ELSE
      ENDIF  
      SP=XH/T2
      IF (X.LT.HO1) THEN
      XF=0.D0
      SSE=0.D0
      ELSE
      XF=WXS1H
      AX=(SP)/((X-X2)**2+(SP)**2)
c     SSE=(GG/EP)*XF*AX*(1.D0/PI)
      IF (X.GE.HO2) THEN
      XF=WXS1H+WXS1L
      AX=(SP)/((X-X3)**2+(SP)**2)
c     SSE=(GG/EP)*XF*AX*(1.D0/PI)
      ELSE
      ENDIF
      ENDIF
      SSE=YQ*XF*AX*(1.D0/PI)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSD(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z7/ XEX
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /Z25/ XXEP                         
      COMMON /XMODE/ MODE,TMODE
      X1=XXEP
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
CCCC  XXEA=Ec XXEB=Eh XXEC=El
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-FC)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-FP)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-FP)/XKT))
      XFD=1.D0/(1+DEXP((XXED-FC)/XKT))  
      ZFCH=XFC*XFH
      ZFCL=XFD*XFL
      ZXFCH=-1.D0+XFC+XFH
      ZXFCL=-1.D0+XFD+XFL
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-FC)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-FP)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-FP)/XKT))
      XFD=1.D0/(1+DEXP((XXED-FC)/XKT))  
      ZFCH=XFD*XFH
      ZFCL=XFC*XFL
      ZXFCH=-1.D0+XFD+XFH
      ZXFCL=-1.D0+XFC+XFL
      ELSE
      ENDIF
      IF ((MODE.EQ.1).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/4.D0)*(1.D0+(HO1/X))
      A1L=(1.D0/4.D0)*(5.D0-3.D0*(HO2/X))   
      WXS1H=XS1H*A1H*ZXFCH
      WXS1L=XS1L*A1L*ZXFCL
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO1/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO2/X)))
      WXS1H=XS1H*A1H*ZXFCH
      WXS1L=XS1L*A1L*ZXFCL
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.2)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO2/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO1/X)))
      WXS1H=XS1L*A1L*ZXFCL
      WXS1L=XS1H*A1H*ZXFCH
      ELSE
      ENDIF     
      SP=XH/T2
      IF (X.LT.HO1) THEN
      XF=0.D0   
      ELSE
      XF=WXS1H
      IF (X.GE.HO2) THEN
      XF=WXS1H+WXS1L             
      ELSE
      ENDIF
      ENDIF
      AX=((SP)/((X-X1)**2+(SP)**2))*(1.D0/PI)
      SSD=(GG/X)*XF*AX
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSF(X)
      IMPLICIT REAL*8 (A-H,O-Z)  
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z2/ WS,GG,GXX,WS2  
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /XMODE/ MODE,TMODE       
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)  
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF    
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))  
      ZFCH=XFC*XFH
      ZFCL=XFD*XFL
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))  
      ZFCH=XFD*XFH
      ZFCL=XFC*XFL
      ELSE
      ENDIF
      IF (TMODE.EQ.1) THEN
      ZXS1H=XS1H*ZFCH
      ZXS1L=XS1L*ZFCL     
      ELSEIF (TMODE.EQ.2) THEN
      ZXS1H=XS1L*ZFCL
      ZXS1L=XS1H*ZFCH
      ELSE
      ENDIF
      IF (X.LT.HO1) THEN
      XF=0.D0
      ELSE
      XF=ZXS1H    
      IF (X.GE.HO2) THEN
      XF=ZXS1H+ZXS1L
      ELSE
      ENDIF
      ENDIF
      SSF=WS2*X*XF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                             
      DOUBLE PRECISION FUNCTION SSI(X)
      IMPLICIT REAL*8 (A-H,O-Z)   
      REAL*8 L,MH,ML
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP  
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /Z25/ XXEP           
      COMMON /XMODE/ MODE,TMODE 
      COMMON /DIFFG/ QM,MH,ML
      COMMON /WIDTH/ L
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF 
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
CCCC  XXEA=Ec XXEB=Eh XXEC=El
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFH=DEXP((XXEB-EFV)/XKT)*(1.D0/(1+DEXP((XXEB-EFV)/XKT)))**2
      XFL=DEXP((XXEC-EFV)/XKT)*(1.D0/(1+DEXP((XXEC-EFV)/XKT)))**2
      XFC=DEXP((XXEA-EFC)/XKT)*(1.D0/(1+DEXP((XXEA-EFC)/XKT)))**2
      XFD=DEXP((XXED-EFC)/XKT)*(1.D0/(1+DEXP((XXED-EFC)/XKT)))**2 
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)   
      XFH=DEXP((XXEB-EFV)/XKT)*(1.D0/(1+DEXP((XXEB-EFV)/XKT)))**2
      XFL=DEXP((XXEC-EFV)/XKT)*(1.D0/(1+DEXP((XXEC-EFV)/XKT)))**2
      XFC=DEXP((XXEA-EFC)/XKT)*(1.D0/(1+DEXP((XXEA-EFC)/XKT)))**2
      XFD=DEXP((XXED-EFC)/XKT)*(1.D0/(1+DEXP((XXED-EFC)/XKT)))**2  
      ELSE
      ENDIF
C      DFF=(1+DEXP((EFC-EC1)/XKT))/DEXP((EFC-EC1)/XKT)
C      DN=((PI*XH**2*L)/QM)*1.D6*DFF
C      DPH=(DEXP((EFV-EH1)/XKT)+1)*(DEXP((EFV-EL1)/XKT)+1)
C      DPL=(DEXP((EFV-EH1)/XKT))*(DEXP((EFV-EL1)/XKT)+1)+((ML/MH)
C     +    *(DEXP((EFV-EL1)/XKT))*(DEXP((EFV-EH1)/XKT)+1))        
C      DP=(((PI*XH**2*L)/MH)*1.D6)*(DPH/DPL)
c      DFF=(1+DEXP((EFC-EC1)/XKT))/DEXP((EFC-EC1)/XKT)
      DFF=1.D0/((DEXP((EFC-EC1)/XKT)/(1+DEXP((EFC-EC1)/XKT)))
     +   +(DEXP((EFC-EC2)/XKT)/(1+DEXP((EFC-EC2)/XKT))))		
      DN=((PI*XH**2*L)/QM)*1.D6*DFF
C      DPH=(DEXP((EFV-EH1)/XKT)+1)*(DEXP((EFV-EL1)/XKT)+1)
C      DPL=(DEXP((EFV-EH1)/XKT))*(DEXP((EFV-EL1)/XKT)+1)+((ML/MH)
C     +    *(DEXP((EFV-EL1)/XKT))*(DEXP((EFV-EH1)/XKT)+1))       
C      DP=(((PI*XH**2*L)/MH)*1.D6)*(DPH/DPL)
      DPHH=1.D0/((DEXP((EFV-EH1)/XKT)/(1+DEXP((EFV-EH1)/XKT)))
     +   +(DEXP((EFV-EH2)/XKT)/(1+DEXP((EFV-EH2)/XKT))))
      DPLL=1.D0/((DEXP((EFV-EL1)/XKT)/(1+DEXP((EFV-EL1)/XKT)))
     +   +(DEXP((EFV-EL2)/XKT)/(1+DEXP((EFV-EL2)/XKT))))
      DPH=((PI*XH**2*L)/MH)*1.D6*DPHH
	DPL=((PI*XH**2*L)/ML)*1.D6*DPLL
C	DP=DPH+DPL
      IF ((MODE.EQ.1).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/4.D0)*(1.D0+(HO1/X))
      A1L=(1.D0/4.D0)*(5.D0-3.D0*(HO2/X))
      WXS1H=XS1H*A1H*((1/XKT)*(XFC*DN+XFH*DPH))
      WXS1L=XS1L*A1L*((1/XKT)*(XFD*DN+XFL*DPL))
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO1/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO2/X)))
      WXS1H=XS1H*A1H*((1/XKT)*(XFC*DN+XFH*DPH))
      WXS1L=XS1L*A1L*((1/XKT)*(XFD*DN+XFL*DPL))
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.2)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO2/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO1/X)))
      WXS1H=XS1L*A1L*((1/XKT)*(XFC*DN+XFL*DPL))
      WXS1L=XS1H*A1H*((1/XKT)*(XFD*DN+XFH*DPH))
      ELSE
      ENDIF                                        
      SP=XH/T2
      IF (X.LT.HO1) THEN
      XF=0.D0
      ELSEIF ((X.GE.HO1).AND.(X.LT.HO2)) THEN
      XF=WXS1H
      AX=((X-X2)/((X-X2)**2+SP**2))*(1.D0/PI)
      ELSEIF (X.GE.HO2) THEN
      XF=WXS1H+WXS1L  
      AX=((X-X3)/((X-X3)**2+SP**2))*(1.D0/PI)
      ELSE
      ENDIF
      SSI=(GG/X)*XF*AX
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSJ(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP 
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /Z25/ XXEP
      COMMON /XMODE/ MODE,TMODE
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
CCCC  XXEA=Ec XXEB=Eh XXEC=El
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))
      ZXFCH=-1.D0+XFC+XFH
      ZXFCL=-1.D0+XFD+XFL
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)
      XFC=1.D0/(1+DEXP((XXEA-EFC)/XKT))
      XFL=1.D0/(1+DEXP((XXEC-EFV)/XKT))
      XFH=1.D0/(1+DEXP((XXEB-EFV)/XKT))
      XFD=1.D0/(1+DEXP((XXED-EFC)/XKT))
      ZXFCH=-1.D0+XFD+XFH
      ZXFCL=-1.D0+XFC+XFL
      ELSE
      ENDIF
      IF ((MODE.EQ.1).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/4.D0)*(1.D0+(HO1/X))
      A1L=(1.D0/4.D0)*(5.D0-3.D0*(HO2/X))
      WXS1H=XS1H*A1H*ZXFCH
      WXS1L=XS1L*A1L*ZXFCL
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO1/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO2/X)))
      WXS1H=XS1H*A1H*ZXFCH
      WXS1L=XS1L*A1L*ZXFCL
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.2)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO2/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO1/X)))
      WXS1H=XS1L*A1L*ZXFCL
      WXS1L=XS1H*A1H*ZXFCH
      ELSE
      ENDIF
      SP=XH/T2
      IF (X.LT.HO1) THEN
      XF=0.D0
      ELSEIF ((X.GE.HO1).AND.(X.LT.HO2)) THEN
      XF=WXS1H
      AX=((X-X2)/((X-X2)**2+SP**2))*(1.D0/PI)
      ELSEIF (X.GE.HO2) THEN
      XF=WXS1H+WXS1L
      AX=((X-X3)/((X-X3)**2+SP**2))*(1.D0/PI) 
      ELSE
      ENDIF
c      IF (XF.LT.0.D0) THEN
c      SSJ=-(GG/X)*XF*AX
c      ELSE
      SSJ=(GG/X)*XF*AX
c      ENDIF
      RETURN
      END  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSK(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 L,MH,ML
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z12/ EC1 
      COMMON /Z14/ EH1,EL1 
      COMMON /Z25/ XXEP
      COMMON /XMODE/ MODE,TMODE
      COMMON /DIFFG/ QM,MH,ML
      COMMON /WIDTH/ L 
      X1=XXEP
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF 
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
CCCC  XXEA=Ec XXEB=Eh XXEC=El
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFH=DEXP((XXEB-FP)/XKT)*(1.D0/(1+DEXP((XXEB-FP)/XKT)))**2
      XFL=DEXP((XXEC-FP)/XKT)*(1.D0/(1+DEXP((XXEC-FP)/XKT)))**2
      XFC=DEXP((XXEA-FC)/XKT)*(1.D0/(1+DEXP((XXEA-FC)/XKT)))**2
      XFD=DEXP((XXED-FC)/XKT)*(1.D0/(1+DEXP((XXED-FC)/XKT)))**2 
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)   
      XFH=DEXP((XXEB-FP)/XKT)*(1.D0/(1+DEXP((XXEB-FP)/XKT)))**2
      XFL=DEXP((XXEC-FP)/XKT)*(1.D0/(1+DEXP((XXEC-FP)/XKT)))**2
      XFC=DEXP((XXEA-FC)/XKT)*(1.D0/(1+DEXP((XXEA-FC)/XKT)))**2
      XFD=DEXP((XXED-FC)/XKT)*(1.D0/(1+DEXP((XXED-FC)/XKT)))**2  
      ELSE
      ENDIF
      DFF=(1+DEXP((FC-EC1)/XKT))/DEXP((FC-EC1)/XKT)
      DN=((PI*XH**2*L)/QM)*1.D6*DFF
      DFH=(1+DEXP((FP-EH1)/XKT))/DEXP((FP-EH1)/XKT)
      DFL=(1+DEXP((FP-EL1)/XKT))/DEXP((FP-EL1)/XKT)
      DP=((PI*XH**2*L)/MH)*1.D6*DFH+((PI*XH**2*L)/ML)*1.D6*DFL
      IF ((MODE.EQ.1).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/4.D0)*(1.D0+(HO1/X))
      A1L=(1.D0/4.D0)*(5.D0-3.D0*(HO2/X))
      WXS1H=XS1H*A1H*((1/XKT)*(XFC*DN+XFH*DP))
      WXS1L=XS1L*A1L*((1/XKT)*(XFD*DN+XFL*DP))
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO1/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO2/X)))
      WXS1H=XS1H*A1H*((1/XKT)*(XFC*DN+XFH*DP))
      WXS1L=XS1L*A1L*((1/XKT)*(XFD*DN+XFL*DP))
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.2)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO2/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO1/X)))
      WXS1H=XS1L*A1L*((1/XKT)*(XFC*DN+XFL*DP))
      WXS1L=XS1H*A1H*((1/XKT)*(XFD*DN+XFH*DP))
      ELSE
      ENDIF
      SP=XH/T2
      IF (X.LT.HO1) THEN
      XF=0.D0
      ELSE
      XF=WXS1H
      IF (X.GE.HO2) THEN
      XF=WXS1H+WXS1L  
      ELSE
      ENDIF
      ENDIF 
      AX=((X-X1)/((X-X1)**2+SP**2))*(1.D0/PI)
      SSK=(GG/X)*XF*AX
      RETURN
      END                      
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      DOUBLE PRECISION FUNCTION SSL(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 L,MH,ML
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z2/ WS,GG,GXX,WS2
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z12/ EC1 
      COMMON /Z14/ EH1,EL1 
      COMMON /Z25/ XXEP
      COMMON /XMODE/ MODE,TMODE
      COMMON /DIFFG/ QM,MH,ML
      COMMON /WIDTH/ L 
      X1=XXEP
      X2=HO1
      X3=HO2
      X4=(XH**2)/(2.D0*XMC1)
      X5=(XH**2)/(2.D0*XMV1H)
      X6=(XH**2)/(2.D0*XMV1L)
      IF (TMODE.EQ.1) THEN
      XK1=(X-X2)/(X4+X5)
      XK2=(X-X3)/(X4+X6)  
      ELSEIF (TMODE.EQ.2) THEN
      XK1=(X-X2)/(X4+X6)
      XK2=(X-X3)/(X4+X5)
      ELSE
      ENDIF 
      XEA=2.D0*XMC1
      XEB=2.D0*XMV1H
      XEC=2.D0*XMV1L
CCCC  XXEA=Ec XXEB=Eh XXEC=El
      IF (TMODE.EQ.1) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-H (C-MASS)
      XXEB=((XH**2)*XK1/XEB)           ! C-H (H-MASS)
      XXEC=((XH**2)*XK2/XEC)           ! C-L (L-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-L (C-MASS)
      XFH=DEXP((XXEB-FP)/XKT)*(1.D0/(1+DEXP((XXEB-FP)/XKT)))**2
      XFL=DEXP((XXEC-FP)/XKT)*(1.D0/(1+DEXP((XXEC-FP)/XKT)))**2
      XFC=DEXP((XXEA-FC)/XKT)*(1.D0/(1+DEXP((XXEA-FC)/XKT)))**2
      XFD=DEXP((XXED-FC)/XKT)*(1.D0/(1+DEXP((XXED-FC)/XKT)))**2 
      ELSEIF (TMODE.EQ.2) THEN
      XXEA=((XH**2)*XK1/XEA)           ! C-L (C-MASS)
      XXEC=((XH**2)*XK1/XEC)           ! C-L (L-MASS)
      XXEB=((XH**2)*XK2/XEB)           ! C-H (H-MASS)
      XXED=((XH**2)*XK2/XEA)           ! C-H (C-MASS)   
      XFH=DEXP((XXEB-FP)/XKT)*(1.D0/(1+DEXP((XXEB-FP)/XKT)))**2
      XFL=DEXP((XXEC-FP)/XKT)*(1.D0/(1+DEXP((XXEC-FP)/XKT)))**2
      XFC=DEXP((XXEA-FC)/XKT)*(1.D0/(1+DEXP((XXEA-FC)/XKT)))**2
      XFD=DEXP((XXED-FC)/XKT)*(1.D0/(1+DEXP((XXED-FC)/XKT)))**2  
      ELSE
      ENDIF
      DFF=(1+DEXP((FC-EC1)/XKT))/DEXP((FC-EC1)/XKT)
      DN=((PI*XH**2*L)/QM)*1.D6*DFF
      DFH=(1+DEXP((FP-EH1)/XKT))/DEXP((FP-EH1)/XKT)
      DFL=(1+DEXP((FP-EL1)/XKT))/DEXP((FP-EL1)/XKT)
      DP=((PI*XH**2*L)/MH)*1.D6*DFH+((PI*XH**2*L)/ML)*1.D6*DFL
      IF ((MODE.EQ.1).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/4.D0)*(1.D0+(HO1/X))
      A1L=(1.D0/4.D0)*(5.D0-3.D0*(HO2/X))
      WXS1H=XS1H*A1H*((1/XKT)*(XFC*DN+XFH*DP))
      WXS1L=XS1L*A1L*((1/XKT)*(XFD*DN+XFL*DP))
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO1/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO2/X)))
      WXS1H=XS1H*A1H*((1/XKT)*(XFC*DN+XFH*DP))
      WXS1L=XS1L*A1L*((1/XKT)*(XFD*DN+XFL*DP))
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.2)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO2/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO1/X)))
      WXS1H=XS1L*A1L*((1/XKT)*(XFC*DN+XFL*DP))
      WXS1L=XS1H*A1H*((1/XKT)*(XFD*DN+XFH*DP))
      ELSE
      ENDIF
      SP=XH/T2
      IF (X.LT.HO1) THEN
      XF=0.D0
      ELSE
      XF=WXS1H
      IF (X.GE.HO2) THEN
      XF=WXS1H+WXS1L  
      ELSE
      ENDIF
      ENDIF 
      AX=((SP)/((X-X1)**2+SP**2))*(1.D0/PI)                   
      SSL=(GG/X)*XF*AX
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SSM(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MODE,TMODE
      REAL*8 MO,MB,N
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZM/ MO
      COMMON /Z27/ MB,EO,N
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /Z1/ XS1H,XS1L
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      COMMON /Z23/ YQ
      COMMON /XMODE/ MODE,TMODE
      TAU=(0.25*1.D-12+0.075*1.D-12)*XH ! (TAUC+TAUV)*XH
      XNG=4.D0
      X2=HO1
      X3=HO2
      IF ((MODE.EQ.1).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/4.D0)*(1.D0+(HO1/X))
      A1L=(1.D0/4.D0)*(5.D0-3.D0*(HO2/X))
      WXS1H=A1H
      WXS1L=A1L
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.1)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO1/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO2/X)))
      WXS1H=A1H
      WXS1L=A1L
      ELSEIF ((MODE.EQ.2).AND.(TMODE.EQ.2)) THEN
      A1H=(3.D0/2.D0)*(1.D0-(HO2/X))
      A1L=(1.D0/2.D0)*(4.D0-3.D0*(1.D0-(HO1/X)))
      WXS1H=A1L
      WXS1L=A1H
      ELSE
      ENDIF                         
      SP=XH/T2
      IF (X.LT.HO1) THEN
      XF=0.D0
      ELSE
      XF=WXS1H
      AX=((SP)/((X-X2)**2+(SP)**2))*(1.D0/PI)
      IF (X.GE.HO2) THEN
      XF=WXS1H+WXS1L
      AX=((SP)/((X-X3)**2+(SP)**2))*(1.D0/PI)
      ELSE
      ENDIF
      ENDIF
      SSM=1.d+6*XF*MB*(1/SP)*(TAU/(MO**2*EO*N*XNG*X))*AX
      RETURN
      END 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                        
      DOUBLE PRECISION FUNCTION SSN(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /NOISE/ A1,A2,WRX,WDX,HIY,PWRF 
      T=2.D0*3.1415926535897932D0
      DDB=1.D0/(1-2*((X/(WRX/T))**2)+(X/(WRX/T))**4+(X/(WDX/T))**2)
      SSN=HIY*PWRF*(((A1+A2*(X/T)**2)/((WRX/T)**4))*DABS(DDB)+1) 
      RETURN                         
      END                         
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                               
      SUBROUTINE DNPEF (L,DNEF,DPEF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 L       
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /Z4/ EFC,EFV,FC,FP
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      DNF=(XMC1/(PI*XH**2*L))*1.D-6
      DNEF=DNF*(DEXP((EFC-EC1)/XKT)/(1+DEXP((EFC-EC1)/XKT)))
      DNPH=(XMV1H/(PI*XH**2*L))*1.D-6
      DNPL=(XMV1L/(PI*XH**2*L))*1.D-6
      DNPH1=DNPH*(DEXP((EFV-EH1)/XKT)/(1+DEXP((EFV-EH1)/XKT)))
      DNPL1=DNPL*(DEXP((EFV-EL1)/XKT)/(1+DEXP((EFV-EL1)/XKT)))
      DPEF=DNPL1+DNPH1
c     DPEF=DNPH1
      RETURN
      END                       
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE DNPEF2 (L,FC,FP,DNEF,DPEF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 L       
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /Z12/ EC1
      COMMON /Z14/ EH1,EL1
      DNF=(XMC1/(PI*XH**2*L))*1.D-6
      DNEF=DNF*(DEXP((FC-EC1)/XKT)/(1+DEXP((FC-EC1)/XKT)))
      DNPH=(XMV1H/(PI*XH**2*L))*1.D-6
      DNPL=(XMV1L/(PI*XH**2*L))*1.D-6
      DNPH1=DNPH*(DEXP((FP-EH1)/XKT)/(1+DEXP((FP-EH1)/XKT)))
      DNPL1=DNPL*(DEXP((FP-EL1)/XKT)/(1+DEXP((FP-EL1)/XKT)))
      DPEF=DNPL1+DNPH1
c     DPEF=DNPH1
      RETURN 
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
C  
      SUBROUTINE XXNP (XNP,LAM,N)   ! Plasma frequency
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LAM,N
      INTEGER MODE,TMODE
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /XMODE/ MODE,TMODE 
      COMMON /Z25/ XXEP
      EO=5.5338D7                       
      CX=2*PI*(3.D10/LAM)
      IF (TMODE.EQ.1) THEN                                  
c        IF (XXEP.LE.HO1) THEN
      XNP=(1.D0/(N*XM1CH*EO*(CX)**2))*1.D6
          ELSEIF (TMODE.EQ.2) THEN
      XNP=(1.D0/(N*XM1CL*EO*(CX)**2))*1.D6
      ELSE  
      XNP=0.D0
      ENDIF                                       
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
      SUBROUTINE XSEFC (K,EFCX,FYI)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KX=500)
      REAL*8 NC(KX),EF(KX),NX(KX),EFCX(KX)
      REAL*8 MO,L,MH,ML
      INTEGER NMAX,FYI
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /XX2/ X1,X2,X3,X4
      COMMON /YY/ E,PP
      COMMON /Z12/ EC1
      COMMON /Z17/ XX,XY,QY,XZ    
      COMMON /WIDTH/ L
      COMMON /DIFFG/ QM,MH,ML
      COMMON /Z29/ EC2,EH2,EL2
      EXTERNAL SSS,SDS
      OPEN(2,FILE='cfermil.dat',status='unknown')
      NMAX=1000
      MO=(0.9109534D-30)/(1.60218D-19)
c     QM=XMC1
      X1=DEXP(-EC1/XKT)
      X2=DEXP(-EC2/XKT)
c	X3=DEXP(-EC3/XKT)
CCCC  QM FOR InGaAsP AND TM FOR InGaAs MATERIAL PARAMETERS
C     L USE CM FOR UNIT
      PSI=((QM*XKT)/(PI*XH**2*L))*1.D-6
CCCC  A,B ARE THE ENDPOINTS
      DO I=1,K
         TOL=1.D-8
         A=-3.0D0
         B=3.0D0
         NC(I)=1.D17+(I-1)*((8.D18-1.D17)/(K-1))
         NX(I)=DEXP(NC(I)/PSI)
         E=NX(I)
  221    CALL NEWBIS (A,B,SSS,SDS,TOL,NMAX,X,IFLG)
         IF(IFLG.LT.0) GOTO 11
           EF(I)=X
           WRITE(2,10) NC(I),EF(I) 
           EFCX(I)=EF(I)
   10      FORMAT(2X,'N-CARRIER=',E20.12,2X,'FERMI-LEVEL=',E20.12)
      ENDDO
   11 WRITE(*,111) I
  111 FORMAT(2x,'FOR EFC SAME SIGN AT THE ENDPOINTS',I5)
      CLOSE (2)
  113 RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION FUNCTION SSS(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /XX2/ X1,X2,X3,X4
      COMMON /YY/ E,PP
      Z1=(1+X1*(DEXP(X/XKT)))
      Z2=(1+X2*(DEXP(X/XKT)))
      SSS=Z1*Z2-E 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION SDS(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /XX2/ X1,X2,X3,X4
      COMMON /YY/ E,PP
      Z1=(1+X1*(DEXP(X/XKT)))
      Z2=(1+X2*(DEXP(X/XKT)))
      ZP1=(X1*(DEXP(X/XKT)))*(1/XKT)
      ZP2=(X2*(DEXP(X/XKT)))*(1/XKT)
      SDS=ZP1*Z2+Z1*ZP2
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE XSEFP (K,EFVX,FYI)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KX=500)
      REAL*8 NP(KX),EP(KX),PX(KX),EFVX(KX)
      REAL*8 MO,L,MH,ML
      INTEGER FYI
      COMMON /ZMASS/ XMC1,XMV1H,XMV1L,XM1CH,XM1CL
      COMMON /ZSSA/ PI,T2,H1,HO2,XH,XKT
      COMMON /XX1/ V1,V2,V3,V4,V5,V6,V7,V8,V9
      COMMON /YY1/ P,D
      COMMON /Z14/ EH1,EL1
      COMMON /Z17/ XX,XY,QY,XZ   
      COMMON /DIFFG/ QM,MH,ML
      COMMON /WIDTH/ L
      COMMON /Z29/ EC2,EH2,EL2
      EXTERNAL PPP,PDP
      OPEN(3,FILE='vfermil.dat',status='unknown')
      NMAX=1000
      MO=(0.9109534D-30)/(1.60218D-19)
c     MH=XMV1H
c     ML=XMV1L
      D=ML/MH
      PSI=(MH*XKT)/(PI*XH**2*L)
      PSI=PSI*1.D-6
      V1=DEXP(-EH1/XKT)
      V2=DEXP(-EL1/XKT)
      V3=DEXP(-EH2/XKT)
      V4=DEXP(-EL2/XKT)
c      V5=DEXP(-EH3/XKT)
c      V6=DEXP(-EL3/XKT)
C	V7=DEXP(-0.103139/XKT)
C	V8=DEXP(-1.57309/XKT)
C	V9=DEXP(-0.208413/XKT)
CCCC  A,B ARE THE ENDPOINTS
      DO I=1,K
         TOL=1.D-8
         A=-3.0D0
         B=3.0D0
         NP(I)=1.D17+(I-1)*((8.D18-1.D17)/(K-1))
c         NP(I)=2.E18
         PX(I)=DEXP(NP(I)/PSI)
         P=PX(I)
         CALL NEWBIS (A,B,PPP,PDP,TOL,NMAX,X,IFLG)
         IF(IFLG.LT.0) GOTO 11
           EP(I)=X
           WRITE(3,10) NP(I),EP(I)
           EFVX(I)=EP(I)
   10      FORMAT(2X,'P-CARRIER=',E20.12,2X,'FERMI-P LEVEL=',E20.12)
      ENDDO
   11 WRITE(*,111) I
  111 FORMAT(2X,'FOR EFV SAME SIGN AT THE ENDPOINTS',I5)
      CLOSE (3)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION PPP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ZSSA/ PI,T2,HO1,HO2,XH,XKT
      COMMON /XX1/ V1,V2,V3,V4,V5,V6,V7,V8,V9
      COMMON /YY1/ P,D
      Z1=(1+V1*(DEXP(X/XKT)))
      Z2=(1+V2*(DEXP(X/XKT)))**D
      Z3=(1+V3*(DEXP(X/XKT)))
      Z4=(1+V4*(DEXP(X/XKT)))**D
      PPP=Z1*Z2*Z3*Z4-P
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION PDP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ZSSA/ PPI,T2,HO1,HO2,XH,XKT
      COMMON /XX1/ V1,V2,V3,V4,V5,V6,V7,V8,V9
      COMMON /YY1/ P,D
      Z1=(1+V1*(DEXP(X/XKT)))
      Z2=(1+V2*(DEXP(X/XKT)))**D
      Z3=(1+V3*(DEXP(X/XKT)))
      Z4=(1+V4*(DEXP(X/XKT)))**D
      ZP1=(V1*(DEXP(X/XKT)))*(1/XKT)
      ZP3=(V3*(DEXP(X/XKT)))*(1/XKT)
      ZP2=D*((1+V2*(DEXP(X/XKT)))**(D-1))*V2*(DEXP(X/XKT))*(1/XKT)
      ZP4=D*((1+V4*(DEXP(X/XKT)))**(D-1))*V4*(DEXP(X/XKT))*(1/XKT)
      PDP=ZP1*Z2*Z3*Z4+Z1*ZP2*Z3*Z4+Z1*Z2*ZP3*Z4+Z1*Z2*Z3*ZP4
      RETURN
      END                        
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          
      SUBROUTINE NEWBIS (A,B,FCN,FDN,TOL,NMAX,X,IFLG)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IFLG,NMAX
      IF(A.LE.B) GOTO 10
      TEMP=A
      A=B
      B=TEMP
   10 FA=FCN(A)
      SFA=SIGN(1.D0,FA)
      FB=FCN(B)
      IFLG=1
      IF(SFA*FB.LE.0.D0) GOTO 3
    2 IFLG=-IFLG
      RETURN
    3 X=B
      FX=FB
      DO 7 K=1,NMAX
         XN=X-FX/FDN(X)
         XTEST=XN
         IF((A.GT.XN).OR.(B.LT.XN)) XTEST=(A+B)/2.D0
           TOLX=TOL*ABS(XTEST)
           DIFF=ABS(XTEST-X)
           IF(XTEST.EQ.XN) GOTO 6
    6        X=XTEST
             IF(DIFF.LE.TOLX) RETURN
               FX=FCN(X)
               IF(SFA*FX.GE.0.) A=X
                 IF(A.EQ.X) GOTO 7
                   B=X
    7 CONTINUE
      WRITE(*,200)
  200 FORMAT(/2X,'MAX. ITERATIONS REACHED'/)
      RETURN
      END 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  
C  
C    __________________________________________________________________  
C   | This is the test program for the self-pulsation problem.         |  
C   |                                                                  |  
C   | Use the subroutine CG (from EISPACK) to find the self-pulsation  |  
C   |                                                                  |  
C   | condition. Also uses the subroutine RK4.f solves the initial     |  
C   |                                                                  |  
C   | value problem.                                                   |  
C   |__________________________________________________________________|  
C                                                                       
      SUBROUTINE RATEEQ
      IMPLICIT REAL*8 (A-H,O-Z)
      WRITE(*,11)
   11 FORMAT(/'                                                     '// 
     +        'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'/
     +        '                                                      '/
     +        '  THE PURPOSES OF THIS PROGRAM ARE FOR THE COUPLED    '/
     +        '  RATE EQUATIONS.  THERE ARE TWO KINDS OF COUPLED     '/
     +        '  RATE EQUATIONS.                                     '/
     +        '                                                      '/
     +        '  THE FIRST ONE ONLY CONSIDER PHOTON AND CARRIER RATE '/
     +        '  EQUATIONS. (GENERAL DYNAMICS ANALYSIS)              '/
     +        '                                                      '/
     +        '  THE SENOND ONE INCLUDES ABSORPTION CARRIER TERM     '/
     +        '  CALLED TWO SECTION MODEL(FOR SELF-PULSATION)        '/
     +        '                                                      '/
     +        'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC') 
    1 CALL SELECTA1(N) 
      GOTO (10,20,30) N
   10 CALL RATE  ! TWO RATE EQUATIONS
      GOTO 1
   20 CALL RATEA ! THREE RATE EQUATIONS
      GOTO 1
   30 PRINT*, ' THIS PROGRAM STOP HERE !'
      RETURN 
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SELECTA1(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
      WRITE(*,11)
   11 FORMAT(/' ENTER 1 FOR GENERAL TWO RATE EQUATIONS   '/
     +        '       2 FOR COUPLED THREE RATE EQUATIONS '/
     +        '       3 FOR EXIT                         '/)
      READ(*,*) N
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE RATE
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NEQN=2)
      REAL*8 T,Y(2),TOUT,RELERR,ABSERR  
      REAL*8 TFINAL,TPRINT,WORK(151) 
      INTEGER IWORK(5),IFLAG
      COMMON /PP1/ D,CS,GACT,BETA1,GOX,XNACT,TP,AA,B,C,xj
      EXTERNAL TRAN
      OPEN(1,FILE='rate.dat',STATUS='unknown')
      OPEN(2,FILE='filex.dat',status='unknown')
    1 REWIND 2
      DO I=1,4
      READ(2,*)
      ENDDO
      READ(2,*) D,GACT,BETA1
      DO I=1,4
      READ(2,*)
      ENDDO
      READ(2,*) GOX,GOXNA,XN
      DO I=1,4
      READ(2,*)
      ENDDO
      READ(2,*) XI,XL,XM1,XM2
      DO I=1,4
      READ(2,*)
      ENDDO
      READ(2,*) AA,B,C
      DO I=1,4
      READ(2,*)
      ENDDO
      READ(2,*) XS,XNACT,XJA
      CS=3.D10/XN
      TP=CS*(XI+(1/(2*XL))*LOG(1/(XM1*XM2)))
      T=0.0  
      Y(1)=XS  
      Y(2)=GOXNA  
      XJ=XJA  
      RELERR=1.e-5 
      ABSERR=0.e0  
      TFINAL=20.e-9  
      TPRINT=0.05e-9  
      IFLAG=1  
      TOUT=T
   10 CALL RKF45 (TRAN,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
      WRITE(6,11) T,Y(1),Y(2)
      WRITE(1,11) T*1.E9,Y(1)/1.e16,Y(2)/1.e18
      GOTO (80,20,30,40,50,60,70,80), IFLAG
   20 TOUT=T+TPRINT
      IF (T.LT.TFINAL) GOTO 10
c     STOP
      GOTO 99
   30 WRITE(6,31) RELERR,ABSERR
      GOTO 10
   40 WRITE(6,41)
      GOTO 10
   50 ABSERR=1.E-9
      WRITE(6,31) RELERR,ABSERR
      GOTO 10
   60 RELERR=10.D0*RELERR
      WRITE(6,31) RELERR,ABSERR
      IFLAG=2
      GOTO 10
   70 WRITE(6,71)
      IFLAG=2
      GOTO 10
   80 WRITE(6,81)
      CLOSE(2)    
c     STOP
C 
   11 FORMAT(3(1x,E19.12))
   31 FORMAT(17H TOLERANCES REST,2E12.3)
   41 FORMAT(11H MANY STEPS)
   71 FORMAT(12H MUCH OUTPUT)
   81 FORMAT(14H IMPROPER CALL)
   99 PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR STOP'
      PRINT*, ' I=?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSEIF (I.EQ.2) THEN
      GOTO 100
      ELSE
      ENDIF
  100 PRINT*, ' PROGRAM STOP HERE, BACK TO MAIN SCREEN!'
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE RATEA   
      IMPLICIT REAL*8 (A-H,O-Z)  
      PARAMETER (ND=3,NEQN=3)  
      COMPLEX*16 A(3,3),EL(ND),EV(ND,ND)  
      REAL*8 T,Y(3),TOUT,RELERR,ABSERR  
      REAL*8 TFINAL,TPRINT,WORK(151),XNACTOA(150),XSOA(150),XXIA(150)    
      REAL*8 XNABSOA(150),XXJA(150),V71A(150)   
      INTEGER IWORK(5),IFLAG  
      COMMON /PARA/ Q,D,CS,TP,GO1,GO2,XN1,XN2,GACT,GABS,TABS,BETA1,  
     +              BETA2,AA,B,C,XXSO,XJ  
      COMMON /PARA2/ XSO,gox,GOY,GOXNA,GOYNA 
      EXTERNAL TRANS,FCN,FDN 
      OPEN(1,FILE='rate.dat',status='unknown')  
      OPEN(2,FILE='data.dat',status='unknown')  
      OPEN(3,FILE='file.dat',status='unknown')  
      OPEN(4,FILE='data1.dat',status='unknown') 
      Q=1.6E-19  
      PRINT*, '                               '  
    1 REWIND 3
      DO I=1,4  
      READ(3,*)  
      ENDDO	      
      READ(3,*) D,XN,TABS  
      CS=3.E10/XN  
      DO I=1,4  
      READ(3,*)  
      ENDDO  
      READ(3,*) GOX,GOY  
      DO I=1,4  
      READ(3,*)  
      ENDDO  
      READ(3,*) GOXNA,GOYNA  
      DO I=1,4  
      READ(3,*)  
      ENDDO  
      READ(3,*) AA,B,C  
      DO I=1,4  
      READ(3,*)  
      ENDDO   
      READ(3,*) BETA1,BETA2  
      DO I=1,4  
      READ(3,*)  
      ENDDO   
      READ(3,*) GACT,GABS   
      DO I=1,4  
      READ(3,*)  
      ENDDO  
      READ(3,*) XL,XW  
      TP=CS*(10+(1/(2*5.E-2))*LOG(1/0.25))  
  100 PRINT*, '                                        ' 
      PRINT*, ' FOR RANGE OF SELF-PULSATION INPUT I=1, '
      PRINT*, ' FOR DATAS OF SELF-PULSATION INPUT I=2  '
	PRINT*, ' INPUT I=?'
	READ(*,*) I
	IF (I.EQ.1) THEN
      PRINT*, ' INPUT THE PHOTON NUMBER XSO=?' 
	READ(*,*) XSO
	GOTO 101
	ELSEIF (I.EQ.2) THEN
	GOTO 200
	ELSE
	ENDIF  
  101 DO I=1,100  
         XSOA(I)=XSO+(I-1)*(XSO/49)  
         XXSO=XSOA(I)  
         AX=2.0E18   
         BX=3.1E18  
         TOL=1.E+3  
         NMAX=150  
         CALL NEWBIS (AX,BX,FCN,FDN,TOL,NMAX,X,IFLG)   
         IF (IFLG.LT.0) GO TO 110  
         WRITE(6,111) AX,BX  
  111    FORMAT(2X,'LEFT END=',2X,E25.20,2X,'RIGHT END=',E25.20/)  
         WRITE(6,112) X  
  112    FORMAT(2x,'TERMINAL ITERATE=',2X,E25.20)  
  110    WRITE(6,114)  
  114    FORMAT(2x,'SAME SIGNS AT THE ENDPOINTS') 
         XNABSOA(I)=X  
	 V71A(I)=(TP-CS*GABS*GOY*(XNABSOA(I)-GOYNA)+CS*GACT*GOX 
     +    *GOXNA)/(CS*GACT*GOX) 
        XNACTOA(I)=V71A(I)   
	XXJA(I)=(Q*D*((AA*XNACTOA(I)**2+B*XNACTOA(I)+C)+CS*GOX 
     +         *(XNACTOA(I)-GOXNA)*xsoa(i))) 
        XXIA(I)=XXJA(I)*XL*XW*1.E-8  
        WRITE(2,15) XXJA(I),XSOA(I),XNACTOA(I),XNABSOA(I),XXIA(I)*1.E3  
	WRITE(4,25) XXJA(I),XSOA(I),XNACTOA(I),XNABSOA(I),XXIA(I)*1.E3 
        WRITE(2,*) '************************************************'  
   15  FORMAT(1X,F10.4,1X,E25.20,1X,E25.20/ 
     +       ,1X,E25.20,1X,E12.6) 
   25  FORMAT(1X,F10.4,1X,E15.8,1X,E15.8,1X,E15.8,1X,F10.4) 
       V1=-(2*AA*XNACTOA(I)+B)-(CS*(GOX)*XSOA(I)) 
       V2=-CS*GOX*(XNACTOA(I)-GOXNA) 
       V3=(-1/TABS)-CS*(GOY)*XSOA(I)  
       V4=-CS*GOY*(XNABSOA(I)-GOYNA)  
       V5=CS*GACT*(GOX)*XSOA(I)+BETA1*(2*AA*XNACTOA(I)+B)  
       V6=CS*GABS*(GOY)*XSOA(I)+BETA2*(1/TABS)  
       V7=CS*(GACT*GOX*(XNACTOA(I)-GOXNA) 
     +     +GABS*GOY*(XNABSOA(I)-GOYNA))-TP       
       A(1,1)=DCMPLX(V1,0.D0)  
       A(1,2)=DCMPLX(0.D0,0.D0)  
       A(1,3)=DCMPLX(V2,0.D0)  
       A(2,1)=DCMPLX(0.D0,0.D0)  
       A(2,2)=DCMPLX(V3,0.D0)  
       A(2,3)=DCMPLX(V4,0.D0)  
       A(3,1)=DCMPLX(V5,0.D0)  
       A(3,2)=DCMPLX(V6,0.D0)  
       A(3,3)=DCMPLX(V7,0.D0)  
       PRINT*, '                                '  
       CALL CGEVV1 (ND,A,EL,EV)  
       DO J=1,3  
          PRINT*, ' EL=',EL(J),' I=',J  
          WRITE(2,16) EL(J),J  
   16     FORMAT(2X,2(E26.20),2X,I3)  
       ENDDO  
      WRITE(2,*) '*********************************************'  
      PRINT*, '                                '  
      WRITE(*,220) XXJA(I),XNACTOA(I) 
      WRITE(*,230) XSOA(I),XNABSOA(I) 
  220 FORMAT(1X,' XXJ=',D15.8,1X,' XNACT=',D15.8) 
  230 FORMAT(1X,' XSO=',D15.8,1X,' XNABS=',D15.8) 
      PRINT*, ' *******************************'    
      ENDDO  
      PRINT*, ' INPUT I=1 FOR NEW SEARCH, I=2 FOR SELF-PULSATION'  
      PRINT*, ' INPUT I=3 STOP'  
      PRINT*, ' I=?'  
      READ(*,*) I  
      IF (I.EQ.1) THEN   
         GOTO 100  
              ELSEIF (I.EQ.2) THEN  
         GOTO 200  
              ELSEIF (I.EQ.3) THEN  
         GOTO 999  
              ELSE  
      ENDIF  
  200 PRINT*, ' INPUT INITIAL VALUE XNACT---ACTIVE CARRIER DENSITY'  
      PRINT*, ' INPUT INITIAL VALUE XS   ---PHOTON DENSITY        '  
      PRINT*, ' INPUT INITIAL VALUE XNABS---ABSORBER CARRIER DENSITY'  
      PRINT*, ' INPUT INITIAL VALUE XXJ  ---CURRENT DENSITY'  
      PRINT*, ' XXJ=?, XS=?, XNACT=?, XNABS=?'  
      READ(*,*) XJA,XS,XNACT,XNABS         
      T=0.E0  
      Y(1)=XS  
      Y(2)=XNACT  
      Y(3)=XNABS  
      XJ=XJA  
      RELERR=1.e-9  
      ABSERR=0.e0  
      TFINAL=50.e-9  
      TPRINT=0.1e-9  
      IFLAG=1  
      TOUT=T  
   10 CALL RKF45 (TRANS,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)  
      WRITE(6,11) T,Y(1),Y(2),Y(3)  
      WRITE(1,11) T*1.E9,Y(1)/1.e15,Y(2)/1.e18,Y(3)/1.e17  
      GOTO (80,20,30,40,50,60,70,80), IFLAG  
   20 TOUT=T+TPRINT  
      IF (T.LT.TFINAL) GOTO 10  
c     STOP  
      GOTO 99
   30 WRITE(6,31) RELERR,ABSERR  
      GOTO 10  
   40 WRITE(6,41)  
      GOTO 10  
   50 ABSERR=1.E-9  
      WRITE(6,31) RELERR,ABSERR  
      GOTO 10  
   60 RELERR=10.D0*RELERR  
      WRITE(6,31) RELERR,ABSERR  
      IFLAG=2  
      GOTO 10  
   70 WRITE(6,71)  
      IFLAG=2  
      GOTO 10  
   80 WRITE(6,81)
      CLOSE(3)
c     STOP  
C  
   11 FORMAT(4(1x,E17.12))  
   31 FORMAT(17H TOLERANCES REST,2E12.3)  
   41 FORMAT(11H MANY STEPS)  
   71 FORMAT(12H MUCH OUTPUT)  
   81 FORMAT(14H IMPROPER CALL)
   99 PRINT*, ' INPUT 1 FOR NEW CALCULATION, 2 FOR STOP'
      PRINT*, ' I=?'
      READ(*,*) I
      IF (I.EQ.1) THEN
      GOTO 1
      ELSEIF (I.EQ.2) THEN
      GOTO 999
      ELSE
      ENDIF
  999 PRINT*, ' THIS PROGRAM STOP HERE'  
      RETURN  
      END  
C  
C
      SUBROUTINE TRAN (T,Y,YP)
      IMPLICIT REAL*8 (A-H,O-Z)  
      REAL*8 T,Y(2),YP(2)
      COMMON /PP1/ D,CS,GACT,BETA1,GOX,XNACT,TP,AA,B,C,XJ
      Q=1.6E-19     
c      YP(1)=(CS*GACT*GOX*(Y(2)-GOXNA)-TP)*Y(1)
c     +      +BETA1*(AA*Y(2)**2+B*Y(2)+C)    
c      YP(2)=(XJ/(Q*D))-(AA*Y(2)**2+B*Y(2)+C)  
c     +      -(CS*GOX*(Y(2)-GOXNA))*Y(1)
c      PRINT*, ' Y1=',Y(1),' Y2=',Y(2)
      YP(1)=(CS*GACT*GOX*(DLOG(Y(2)/XNACT)+1)-TP)*Y(1)
     +      +GACT*BETA1*(AA*Y(2)**3+B*Y(2)**2.D0+C) 
 	YP(2)=(XJ/(Q*D))-(AA*Y(2)**3+B*Y(2)**2.D0+C)
     +      -(CS*GACT*GOX*(DLOG(Y(2)/XNACT)+1))*Y(1)	 
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE TRANS (T,Y,YP)  
      IMPLICIT REAL*8 (A-H,O-Z)  
      REAL*8 T,Y(3),YP(3)   
      COMMON /PARA/ Q,D,CS,TP,GO1,GO2,XN1,XN2,GACT,GABS,TABS,BETA1,  
     +              BETA2,AA,B,C,XXSO,XJ  
      COMMON /PARA2/ XSO,GOX,GOY,GOXNA,GOYNA 
      YP(1)=(CS*(GACT*GOX*(Y(2)-GOXNA)+GABS*GOY  
     +       *(Y(3)-GOYNA))  
     +       -TP)*Y(1)+BETA1*(AA*Y(2)**2+B*Y(2)+C)+  
     +      (BETA2*Y(3)/TABS)  
      YP(2)=(XJ/(Q*D))-(AA*Y(2)**2+B*Y(2)+C)  
     +      -(CS*GOX*(Y(2)-GOXNA))*Y(1)     
      YP(3)=(-Y(3)/TABS)-CS*GOY*(Y(3)-GOYNA)*Y(1)  
      RETURN  
      END  
c  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c  
      REAL*8 FUNCTION FCN(X)  
      IMPLICIT REAL*8 (A-H,O-Z)  
      COMMON /PARA/ Q,D,CS,TP,GO1,GO2,XN1,XN2,GACT,GABS,TABS,BETA1,  
     +              BETA2,AA,B,C,XXSO,XJ  
      COMMON /PARA2/ XSO,GOX,GOY,GOXNA,GOYNA 
      FCN=((-X/TABS)*(1/(CS*GOY*(X-GOYNA))))-XXSO  
      RETURN  
      END  
        
      REAL*8 FUNCTION FDN(X)  
      IMPLICIT REAL*8 (A-H,O-Z)  
      COMMON /PARA/ Q,D,CS,TP,GO1,GO2,XN1,XN2,GACT,GABS,TABS,BETA1,  
     +              BETA2,AA,B,C,XXSO,XJ  
      COMMON /PARA2/ XSO,GOX,GOY,GOXNA,GOYNA 
      FDN=(-1.D0/TABS)*(1/(CS*GOY*(X-GOYNA)))  
     +    +(-1)*(-1*X/TABS)  
     +    *(CS*GOY*(X-GOYNA))**(-2)  
     +    *(CS*GOY)  
      RETURN  
      END  
C  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C  
      subroutine rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)  
c  
c     fehlberg fourth-fifth order runge-kutta method  
c  
c     written by h.a.watts and l.f.shampine  
c                   sandia laboratories  
c                  albuquerque,new mexico  
c  
c    rkf45 is primarily designed to solve non-stiff and mildly stiff  
c    differential equations when derivative evaluations are inexpensive.  
c    rkf45 should generally not be used when the user is demanding  
c    high accuracy.  
c  
c abstract  
c  
c    subroutine  rkf45  integrates a system of neqn first order  
c    ordinary differential equations of the form  
c             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))  
c              where the y(i) are given at t .  
c    typically the subroutine is used to integrate from t to tout but it  
c    can be used as a one-step integrator to advance the solution a  
c    single step in the direction of tout.  on return the parameters in  
c    the call list are set for continuing the integration. the user has  
c    only to call rkf45 again (and perhaps define a new value for tout).  
c    actually, rkf45 is an interfacing routine which calls subroutine  
c    rkfs for the solution.  rkfs in turn calls subroutine  fehl which  
c    computes an approximate solution over one step.  
c  
c    rkf45  uses the runge-kutta-fehlberg (4,5)  method described  
c    in the reference  
c    e.fehlberg , low-order classical runge-kutta formulas with stepsize  
c                 control , nasa tr r-315  
c  
c    the performance of rkf45 is illustrated in the reference  
c    l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary  
c                 differential equations-the state of the art ,  
c                 sandia laboratories report sand75-0182 ,  
c                 to appear in siam review.  
c  
c  
c    the parameters represent-  
c      f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt  
c      neqn -- number of equations to be integrated  
c      y(*) -- solution vector at t  
c      t -- independent variable  
c      tout -- output point at which solution is desired  
c      relerr,abserr -- relative and absolute error tolerances for local  
c            error test. at each step the code requires that  
c                 abs(local error) .le. relerr*abs(y) + abserr  
c            for each component of the local error and solution vectors  
c      iflag -- indicator for status of integration  
c      work(*) -- array to hold information internal to rkf45 which is  
c            necessary for subsequent calls. must be dimensioned  
c            at least  3+6*neqn  
c      iwork(*) -- integer array used to hold information internal to  
c            rkf45 which is necessary for subsequent calls. must be  
c            dimensioned at least  5  
c  
c  
c  first call to rkf45  
c  
c    the user must provide storage in his calling program for the arrays  
c    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,  
c    declare f in an external statement, supply subroutine f(t,y,yp) and  
c    initialize the following parameters-  
c  
c      neqn -- number of equations to be integrated.  (neqn .ge. 1)  
c      y(*) -- vector of initial conditions  
c      t -- starting point of integration , must be a variable  
c      tout -- output point at which solution is desired.  
c            t=tout is allowed on the first call only, in which case  
c            rkf45 returns with iflag=2 if continuation is possible.  
c      relerr,abserr -- relative and absolute local error tolerances  
c            which must be non-negative. relerr must be a variable while  
c            abserr may be a constant. the code should normally not be  
c            used with relative error control smaller than about 1.e-8 .  
c            to avoid limiting precision difficulties the code requires  
c            relerr to be larger than an internally computed relative  
c            error parameter which is machine dependent. in particular,  
c            pure absolute error is not permitted. if a smaller than  
c            allowable value of relerr is attempted, rkf45 increases  
c            relerr appropriately and returns control to the user before  
c            continuing the integration.  
c      iflag -- +1,-1  indicator to initialize the code for each new  
c            problem. normal input is +1. the user should set iflag=-1  
c            only when one-step integrator control is essential. in this  
c            case, rkf45 attempts to advance the solution a single step  
c            in the direction of tout each time it is called. since this  
c            mode of operation results in extra computing overhead, it  
c            should be avoided unless needed.  
c  
c  
c  output from rkf45  
c  
c      y(*) -- solution at t  
c      t -- last point reached in integration.  
c      iflag = 2 -- integration reached tout. indicates successful retur  
c                   and is the normal mode for continuing integration.  
c            =-2 -- a single successful step in the direction of tout  
c                   has been taken. normal mode for continuing  
c                   integration one step at a time.  
c            = 3 -- integration was not completed because relative error  
c                   tolerance was too small. relerr has been increased  
c                   appropriately for continuing.  
c            = 4 -- integration was not completed because more than  
c                   3000 derivative evaluations were needed. this  
c                   is approximately 500 steps.  
c            = 5 -- integration was not completed because solution  
c                   vanished making a pure relative error test  
c                   impossible. must use non-zero abserr to continue.  
c                   using the one-step integration mode for one step  
c                   is a good way to proceed.  
c            = 6 -- integration was not completed because requested  
c                   accuracy could not be achieved using smallest  
c                   allowable stepsize. user must increase the error  
c                   tolerance before continued integration can be  
c                   attempted.  
c            = 7 -- it is likely that rkf45 is inefficient for solving  
c                   this problem. too much output is restricting the  
c                   natural stepsize choice. use the one-step integrator  
c                   mode.  
c            = 8 -- invalid input parameters  
c                   this indicator occurs if any of the following is  
c                   satisfied -   neqn .le. 0  
c                                 t=tout  and  iflag .ne. +1 or -1  
c                                 relerr or abserr .lt. 0.  
c                                 iflag .eq. 0  or  .lt. -2  or  .gt. 8  
c      work(*),iwork(*) -- information which is usually of no interest  
c                   to the user but necessary for subsequent calls.  
c                   work(1),...,work(neqn) contain the first derivatives  
c                   of the solution vector y at t. work(neqn+1) contains  
c                   the stepsize h to be attempted on the next step.  
c                   iwork(1) contains the derivative evaluation counter.  
c  
c  
c  subsequent calls to rkf45  
c  
c    subroutine rkf45 returns with all information needed to continue  
c    the integration. if the integration reached tout, the user need onl  
c    define a new tout and call rkf45 again. in the one-step integrator  
c    mode (iflag=-2) the user must keep in mind that each step taken is  
c    in the direction of the current tout. upon reaching tout (indicated  
c    by changing iflag to 2),the user must then define a new tout and  
c    reset iflag to -2 to continue in the one-step integrator mode.  
c  
c    if the integration was not completed but the user still wants to  
c    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3  
c    the relerr parameter has been adjusted appropriately for continuing  
c    the integration. in the case of iflag=4 the function counter will  
c    be reset to 0 and another 3000 function evaluations are allowed.  
c  
c    however,in the case iflag=5, the user must first alter the error  
c    criterion to use a positive value of abserr before integration can  
c    proceed. if he does not,execution is terminated.  
c  
c    also,in the case iflag=6, it is necessary for the user to reset  
c    iflag to 2 (or -2 when the one-step integration mode is being used)  
c    as well as increasing either abserr,relerr or both before the  
c    integration can be continued. if this is not done, execution will  
c    be terminated. the occurrence of iflag=6 indicates a trouble spot  
c    (solution is changing rapidly,singularity may be present) and it  
c    often is inadvisable to continue.  
c  
c    if iflag=7 is encountered, the user should use the one-step  
c    integration mode with the stepsize determined by the code or  
c    consider switching to the adams codes de/step,intrp. if the user  
c    insists upon continuing the integration with rkf45, he must reset  
c    iflag to 2 before calling rkf45 again. otherwise,execution will be  
c    terminated.  
c  
c    if iflag=8 is obtained, integration can not be continued unless  
c    the invalid input parameters are corrected.  
c  
c    it should be noted that the arrays work,iwork contain information  
c    required for subsequent integration. accordingly, work and iwork  
c    should not be altered.  
c  
c  
      integer neqn,iflag,iwork(5)  
      double precision y(neqn),t,tout,relerr,abserr,work(1)  
c     if compiler checks subscripts, change work(1) to work(3+6*neqn)  
c  
      external f  
c  
      integer k1,k2,k3,k4,k5,k6,k1m  
c  
c  
c     compute indices for the splitting of the work array  
c  
      k1m=neqn+1  
      k1=k1m+1  
      k2=k1+neqn  
      k3=k2+neqn  
      k4=k3+neqn  
      k5=k4+neqn  
      k6=k5+neqn  
c  
c     this interfacing routine merely relieves the user of a long  
c     calling list via the splitting apart of two working storage  
c     arrays. if this is not compatible with the users compiler,  
c     he must use rkfs directly.  
c  
      call rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,work(1),work(k1m),  
     1          work(k1),work(k2),work(k3),work(k4),work(k5),work(k6),  
     2          work(k6+1),iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))  
c  
      return  
      end  
      subroutine rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,yp,h,f1,f2,f3,  
     1                f4,f5,savre,savae,nfe,kop,init,jflag,kflag)  
c  
c     fehlberg fourth-fifth order runge-kutta method  
c  
c  
c     rkfs integrates a system of first order ordinary differential  
c     equations as described in the comments for rkf45 .  
c     the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and  
c     the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used  
c     internally by the code and appear in the call list to eliminate  
c     local retention of variables between calls. accordingly, they  
c     should not be altered. items of possible interest are  
c         yp - derivative of solution vector at t  
c         h  - an appropriate stepsize to be used for the next step  
c         nfe- counter on the number of derivative function evaluations  
c  
c  
      logical hfaild,output  
c  
      integer  neqn,iflag,nfe,kop,init,jflag,kflag  
      double precision  y(neqn),t,tout,relerr,abserr,h,yp(neqn),  
     1  f1(neqn),f2(neqn),f3(neqn),f4(neqn),f5(neqn),savre,  
     2  savae  
c  
      external f  
c  
      double precision  a,ae,dt,ee,eeoet,esttol,et,hmin,remin,rer,s,  
     1  scale,tol,toln,u26,epsp1,eps,ypk  
c  
      integer  k,maxnfe,mflag  
c  
      double precision  dabs,dmax1,dmin1,dsign  
c  
c     remin is the minimum acceptable value of relerr.  attempts  
c     to obtain higher accuracy with this subroutine are usually  
c     very expensive and often unsuccessful.  
c  
      data remin/1.d-12/  
c  
c  
c     the expense is controlled by restricting the number  
c     of function evaluations to be approximately maxnfe.  
c     as set, this corresponds to about 500 steps.  
c  
      data maxnfe/3000/  
c  
c  
c     check input parameters  
c  
c  
      if (neqn .lt. 1) go to 10  
      if ((relerr .lt. 0.0d0)  .or.  (abserr .lt. 0.0d0)) go to 10  
      mflag=iabs(iflag)  
      if ((mflag .eq. 0) .or. (mflag .gt. 8)) go to 10  
      if (mflag .ne. 1) go to 20  
c  
c     first call, compute machine epsilon  
c  
      eps = 1.0d0  
    5 eps = eps/2.0d0  
      epsp1 = eps + 1.0d0  
      if (epsp1 .gt. 1.0d0) go to 5  
      u26 = 26.0d0*eps  
      go to 50  
c  
c     invalid input  
   10 iflag=8  
      return  
c  
c     check continuation possibilities  
c  
   20 if ((t .eq. tout) .and. (kflag .ne. 3)) go to 10  
      if (mflag .ne. 2) go to 25  
c  
c     iflag = +2 or -2  
      if ((kflag .eq. 3) .or. (init .eq. 0)) go to 45  
      if (kflag .eq. 4) go to 40  
      if ((kflag .eq. 5)  .and.  (abserr .eq. 0.0d0)) go to 30  
      if ((kflag .eq. 6)  .and.  (relerr .le. savre)  .and.  
     1    (abserr .le. savae)) go to 30  
      go to 50  
c  
c     iflag = 3,4,5,6,7 or 8  
   25 if (iflag .eq. 3) go to 45  
      if (iflag .eq. 4) go to 40  
      if ((iflag .eq. 5) .and. (abserr .gt. 0.0d0)) go to 45  
c  
c     integration cannot be continued since user did not respond to  
c     the instructions pertaining to iflag=5,6,7 or 8  
   30 stop  
c  
c     reset function evaluation counter  
   40 nfe=0  
      if (mflag .eq. 2) go to 50  
c  
c     reset flag value from previous call  
   45 iflag=jflag  
      if (kflag .eq. 3) mflag=iabs(iflag)  
c  
c     save input iflag and set continuation flag value for subsequent  
c     input checking  
   50 jflag=iflag  
      kflag=0  
c  
c     save relerr and abserr for checking input on subsequent calls  
      savre=relerr  
      savae=abserr  
c  
c     restrict relative error tolerance to be at least as large as  
c     2*eps+remin to avoid limiting precision difficulties arising  
c     from impossible accuracy requests  
c  
      rer=2.0d0*eps+remin  
      if (relerr .ge. rer) go to 55  
c  
c     relative error tolerance too small  
      relerr=rer  
      iflag=3  
      kflag=3  
      return  
c  
   55 dt=tout-t  
c  
      if (mflag .eq. 1) go to 60  
      if (init .eq. 0) go to 65  
      go to 80  
c  
c     initialization --  
c                       set initialization completion indicator,init  
c                       set indicator for too many output points,kop  
c                       evaluate initial derivatives  
c                       set counter for function evaluations,nfe  
c                       estimate starting stepsize  
c  
   60 init=0  
      kop=0  
c  
      a=t  
      call f(a,y,yp)  
      nfe=1  
      if (t .ne. tout) go to 65  
      iflag=2  
      return  
c  
c  
   65 init=1  
      h=dabs(dt)  
      toln=0.  
      do 70 k=1,neqn  
        tol=relerr*dabs(y(k))+abserr  
        if (tol .le. 0.) go to 70  
        toln=tol  
        ypk=dabs(yp(k))  
        if (ypk*h**5 .gt. tol) h=(tol/ypk)**0.2d0  
   70 continue  
      if (toln .le. 0.0d0) h=0.0d0  
      h=dmax1(h,u26*dmax1(dabs(t),dabs(dt)))  
      jflag=isign(2,iflag)  
c  
c  
c     set stepsize for integration in the direction from t to tout  
c  
   80 h=dsign(h,dt)  
c  
c     test to see if rkf45 is being severely impacted by too many  
c     output points  
c  
      if (dabs(h) .ge. 2.0d0*dabs(dt)) kop=kop+1  
      if (kop .ne. 100) go to 85  
c  
c     unnecessary frequency of output  
      kop=0  
      iflag=7  
      return  
c  
   85 if (dabs(dt) .gt. u26*dabs(t)) go to 95  
c  
c     if too close to output point,extrapolate and return  
c  
      do 90 k=1,neqn  
   90   y(k)=y(k)+dt*yp(k)  
      a=tout  
      call f(a,y,yp)  
      nfe=nfe+1  
      go to 300  
c  
c  
c     initialize output point indicator  
c  
   95 output= .false.  
c  
c     to avoid premature underflow in the error tolerance function,  
c     scale the error tolerances  
c  
      scale=2.0d0/relerr  
      ae=scale*abserr  
c  
c  
c     step by step integration  
c  
  100 hfaild= .false.  
c  
c     set smallest allowable stepsize  
c  
      hmin=u26*dabs(t)  
c  
c     adjust stepsize if necessary to hit the output point.  
c     look ahead two steps to avoid drastic changes in the stepsize and  
c     thus lessen the impact of output points on the code.  
c  
      dt=tout-t  
      if (dabs(dt) .ge. 2.0d0*dabs(h)) go to 200  
      if (dabs(dt) .gt. dabs(h)) go to 150  
c  
c     the next successful step will complete the integration to the  
c     output point  
c  
      output= .true.  
      h=dt  
      go to 200  
c  
  150 h=0.5d0*dt  
c  
c  
c  
c     core integrator for taking a single step  
c  
c     the tolerances have been scaled to avoid premature underflow in  
c     computing the error tolerance function et.  
c     to avoid problems with zero crossings,relative error is measured  
c     using the average of the magnitudes of the solution at the  
c     beginning and end of a step.  
c     the error estimate formula has been grouped to control loss of  
c     significance.  
c     to distinguish the various arguments, h is not permitted  
c     to become smaller than 26 units of roundoff in t.  
c     practical limits on the change in the stepsize are enforced to  
c     smooth the stepsize selection process and to avoid excessive  
c     chattering on problems having discontinuities.  
c     to prevent unnecessary failures, the code uses 9/10 the stepsize  
c     it estimates will succeed.  
c     after a step failure, the stepsize is not allowed to increase for  
c     the next attempted step. this makes the code more efficient on  
c     problems having discontinuities and more effective in general  
c     since local extrapolation is being used and extra caution seems  
c     warranted.  
c  
c  
c     test number of derivative function evaluations.  
c     if okay,try to advance the integration from t to t+h  
c  
  200 if (nfe .le. maxnfe) go to 220  
c  
c     too much work  
      iflag=4  
      kflag=4  
      return  
c  
c     advance an approximate solution over one step of length h  
c  
  220 call fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,f1)  
      nfe=nfe+5  
c  
c     compute and test allowable tolerances versus local error estimates  
c     and remove scaling of tolerances. note that relative error is  
c     measured with respect to the average of the magnitudes of the  
c     solution at the beginning and end of the step.  
c  
      eeoet=0.0d0  
      do 250 k=1,neqn  
        et=dabs(y(k))+dabs(f1(k))+ae  
        if (et .gt. 0.0d0) go to 240  
c  
c       inappropriate error tolerance  
        iflag=5  
        return  
c  
  240   ee=dabs((-2090.0d0*yp(k)+(21970.0d0*f3(k)-15048.0d0*f4(k)))+  
     1                        (22528.0d0*f2(k)-27360.0d0*f5(k)))  
  250   eeoet=dmax1(eeoet,ee/et)  
c  
      esttol=dabs(h)*eeoet*scale/752400.0d0  
c  
      if (esttol .le. 1.0d0) go to 260  
c  
c  
c     unsuccessful step  
c                       reduce the stepsize , try again  
c                       the decrease is limited to a factor of 1/10  
c  
      hfaild= .true.  
      output= .false.  
      s=0.1d0  
      if (esttol .lt. 59049.0d0) s=0.9d0/esttol**0.2d0  
      h=s*h  
      if (dabs(h) .gt. hmin) go to 200  
c  
c     requested error unattainable at smallest allowable stepsize  
      iflag=6  
      kflag=6  
      return  
c  
c  
c     successful step  
c                        store solution at t+h  
c                        and evaluate derivatives there  
c  
  260 t=t+h  
      do 270 k=1,neqn  
  270   y(k)=f1(k)  
      a=t  
      call f(a,y,yp)  
      nfe=nfe+1  
c  
c  
c                       choose next stepsize  
c                       the increase is limited to a factor of 5  
c                       if step failure has just occurred, next  
c                          stepsize is not allowed to increase  
c  
      s=5.0d0  
      if (esttol .gt. 1.889568d-4) s=0.9d0/esttol**0.2d0  
      if (hfaild) s=dmin1(s,1.0d0)  
      h=dsign(dmax1(s*dabs(h),hmin),h)  
c  
c     end of core integrator  
c  
c  
c     should we take another step  
c  
      if (output) go to 300  
      if (iflag .gt. 0) go to 100  
c  
c  
c     integration successfully completed  
c  
c     one-step mode  
      iflag=-2  
      return  
c  
c     interval mode  
  300 t=tout  
      iflag=2  
      return  
c  
      end  
      subroutine fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,s)  
c  
c     fehlberg fourth-fifth order runge-kutta method  
c  
c    fehl integrates a system of neqn first order  
c    ordinary differential equations of the form  
c             dy(i)/dt=f(t,y(1),---,y(neqn))  
c    where the initial values y(i) and the initial derivatives  
c    yp(i) are specified at the starting point t. fehl advances  
c    the solution over the fixed step h and returns  
c    the fifth order (sixth order accurate locally) solution  
c    approximation at t+h in array s(i).  
c    f1,---,f5 are arrays of dimension neqn which are needed  
c    for internal storage.  
c    the formulas have been grouped to control loss of significance.  
c    fehl should be called with an h not smaller than 13 units of  
c    roundoff in t so that the various independent arguments can be  
c    distinguished.  
c  
c  
      integer  neqn  
      double precision  y(neqn),t,h,yp(neqn),f1(neqn),f2(neqn),  
     1  f3(neqn),f4(neqn),f5(neqn),s(neqn)  
c  
      double precision  ch  
      integer  k  
c  
      ch=h/4.0d0  
      do 221 k=1,neqn  
  221   f5(k)=y(k)+ch*yp(k)  
      call f(t+ch,f5,f1)  
c  
      ch=3.0d0*h/32.0d0  
      do 222 k=1,neqn  
  222   f5(k)=y(k)+ch*(yp(k)+3.0d0*f1(k))  
      call f(t+3.0d0*h/8.0d0,f5,f2)  
c  
      ch=h/2197.0d0  
      do 223 k=1,neqn  
  223   f5(k)=y(k)+ch*(1932.0d0*yp(k)+(7296.0d0*f2(k)-7200.0d0*f1(k)))  
      call f(t+12.0d0*h/13.0d0,f5,f3)  
c  
      ch=h/4104.0d0  
      do 224 k=1,neqn  
  224   f5(k)=y(k)+ch*((8341.0d0*yp(k)-845.0d0*f3(k))+  
     1                            (29440.0d0*f2(k)-32832.0d0*f1(k)))  
      call f(t+h,f5,f4)  
c  
      ch=h/20520.0d0  
      do 225 k=1,neqn  
  225   f1(k)=y(k)+ch*((-6080.0d0*yp(k)+(9295.0d0*f3(k)-  
     1         5643.0d0*f4(k)))+(41040.0d0*f1(k)-28352.0d0*f2(k)))  
      call f(t+h/2.0d0,f1,f5)  
c  
c     compute approximate solution at t+h  
c  
      ch=h/7618050.0d0  
      do 230 k=1,neqn  
  230   s(k)=y(k)+ch*((902880.0d0*yp(k)+(3855735.0d0*f3(k)-  
     1        1371249.0d0*f4(k)))+(3953664.0d0*f2(k)+  
     2        277020.0d0*f5(k)))  
c  
      return  
      end  
  
C  
C      ________________________________________________________  
  
  
  
        subroutine cgevv1(n, a, el, ev)  
c  
           integer i, j, k, n, mw, ierr  
           parameter(mw=3)  
           real*8 ar(mw,mw), ai(mw,mw), wr(mw), wi(mw), tr, ti  
           real*8 zr(mw,mw), zi(mw,mw), fv1(mw), fv2(mw), fv3(mw)  
           complex*16 a(mw,mw), el(mw), ev(mw,mw)  
c  
           do 100 i=1,n  
              do 100 j=1,n  
                    ar(i,j)=dble(a(i,j))  
                    ai(i,j)=dimag(a(i,j))  
 100       continue  
c  
           call cg(n,n,ar,ai,wr,wi,0,zr,zi,fv1,fv2,fv3,ierr)  
c  
           do 600 i=1,n-1  
              do 500 j=i+1,n  
                    if (wr(i).gt.wr(j)) then  
                          tr=wr(j)  
                          ti=wi(j)  
                          wr(j)=wr(i)  
                          wi(j)=wi(i)  
                          wr(i)=tr  
                          wi(i)=ti  
                       do 400 k=1,n  
                             tr=zr(k,j)  
                             ti=zi(k,j)  
                             zr(k,j)=zr(k,i)  
                             zi(k,j)=zi(k,i)  
                             zr(k,i)=tr  
                             zi(k,i)=ti  
 400                   continue  
                    end if  
 500          continue  
 600       continue  
c  
           do 200 i=1,n  
                 el(i)=dcmplx(wr(i),wi(i))  
              do 200 j=1,n  
                    ev(i,j)=dcmplx(zr(i,j),zi(i,j))  
 200       continue  
c  
        end  
  
  
  
c******************************************************************************  
c    Solving general complex eigen system problem  
c******************************************************************************  
c  
      subroutine cg(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)  
c  
      integer n,nm,is1,is2,ierr,matz  
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),  
     x       fv1(n),fv2(n),fv3(n)  
c  
c     this subroutine calls the recommended sequence of  
c     subroutines from the eigensystem subroutine package (eispack)  
c     to find the eigenvalues and eigenvectors (if desired)  
c     of a complex general matrix.  
c  
c     on input  
c  
c        nm  must be set to the row dimension of the two-dimensional  
c        array parameters as declared in the calling program  
c        dimension statement.  
c  
c        n  is the order of the matrix  a=(ar,ai).  
c  
c        ar  and  ai  contain the real and imaginary parts,  
c        respectively, of the complex general matrix.  
c  
c        matz  is an integer variable set equal to zero if  
c        only eigenvalues are desired.  otherwise it is set to  
c        any non-zero integer for both eigenvalues and eigenvectors.  
c  
c     on output  
c  
c        wr  and  wi  contain the real and imaginary parts,  
c        respectively, of the eigenvalues.  
c  
c        zr  and  zi  contain the real and imaginary parts,  
c        respectively, of the eigenvectors if matz is not zero.  
c  
c        ierr  is an integer output variable set equal to an error  
c           completion code described in the documentation for comqr  
c           and comqr2.  the normal completion code is zero.  
c  
c        fv1, fv2, and  fv3  are temporary storage arrays.  
c  
c     questions and comments should be directed to burton s. garbow,  
c     mathematics and computer science div, argonne national laboratory  
c  
c     this version dated august 1983.  
c  
c     ------------------------------------------------------------------  
c  
      if (n .le. nm) go to 10  
      ierr = 10 * n  
      go to 50  
c  
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)  
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)  
      if (matz .ne. 0) go to 20  
c     .......... find eigenvalues only ..........  
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)  
      go to 50  
c     .......... find both eigenvalues and eigenvectors ..........  
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)  
      if (ierr .ne. 0) go to 50  
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)  
   50 return  
      end  
  
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)  
c  
      integer i,j,k,m,n,ii,nm,igh,low  
      double precision scale(n),zr(nm,m),zi(nm,m)  
      double precision s  
c  
c     this subroutine is a translation of the algol procedure  
c     cbabk2, which is a complex version of balbak,  
c     num. math. 13, 293-304(1969) by parlett and reinsch.  
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).  
c  
c     this subroutine forms the eigenvectors of a complex general  
c     matrix by back transforming those of the corresponding  
c     balanced matrix determined by  cbal.  
c  
c     on input  
c  
c        nm must be set to the row dimension of two-dimensional  
c          array parameters as declared in the calling program  
c          dimension statement.  
c  
c        n is the order of the matrix.  
c  
c        low and igh are integers determined by  cbal.  
c  
c        scale contains information determining the permutations  
c          and scaling factors used by  cbal.  
c  
c        m is the number of eigenvectors to be back transformed.  
c  
c        zr and zi contain the real and imaginary parts,  
c          respectively, of the eigenvectors to be  
c          back transformed in their first m columns.  
c  
c     on output  
c  
c        zr and zi contain the real and imaginary parts,  
c          respectively, of the transformed eigenvectors  
c          in their first m columns.  
c  
c     questions and comments should be directed to burton s. garbow,  
c     mathematics and computer science div, argonne national laboratory  
c  
c     this version dated august 1983.  
c  
c     ------------------------------------------------------------------  
c  
      if (m .eq. 0) go to 200  
      if (igh .eq. low) go to 120  
c  
      do 110 i = low, igh  
         s = scale(i)  
c     .......... left hand eigenvectors are back transformed  
c                if the foregoing statement is replaced by  
c                s=1.0d0/scale(i). ..........  
         do 100 j = 1, m  
            zr(i,j) = zr(i,j) * s  
            zi(i,j) = zi(i,j) * s  
  100    continue  
c  
  110 continue  
c     .......... for i=low-1 step -1 until 1,  
c                igh+1 step 1 until n do -- ..........  
  120 do 140 ii = 1, n  
         i = ii  
         if (i .ge. low .and. i .le. igh) go to 140  
         if (i .lt. low) i = low - ii  
         k = scale(i)  
         if (k .eq. i) go to 140  
c  
         do 130 j = 1, m  
            s = zr(i,j)  
            zr(i,j) = zr(k,j)  
            zr(k,j) = s  
            s = zi(i,j)  
            zi(i,j) = zi(k,j)  
            zi(k,j) = s  
  130    continue  
c  
  140 continue  
c  
  200 return  
      end  
  
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)  
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)  
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)  
c  
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,  
     x        itn,its,low,lp1,enm1,iend,ierr  
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),  
     x       ortr(igh),orti(igh)  
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,  
     x       pythag  
c  
c     this subroutine is a translation of a unitary analogue of the  
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters  
c     and wilkinson.  
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).  
c     the unitary analogue substitutes the qr algorithm of francis  
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.  
c  
c     this subroutine finds the eigenvalues and eigenvectors  
c     of a complex upper hessenberg matrix by the qr  
c     method.  the eigenvectors of a complex general matrix  
c     can also be found if  corth  has been used to reduce  
c     this general matrix to hessenberg form.  
c  
c     on input  
c  
c        nm must be set to the row dimension of two-dimensional  
c          array parameters as declared in the calling program  
c          dimension statement.  
c  
c        n is the order of the matrix.  
c  
c        low and igh are integers determined by the balancing  
c          subroutine  cbal.  if  cbal  has not been used,  
c          set low=1, igh=n.  
c  
c        ortr and orti contain information about the unitary trans-  
c          formations used in the reduction by  corth, if performed.  
c          only elements low through igh are used.  if the eigenvectors  
c          of the hessenberg matrix are desired, set ortr(j) and  
c          orti(j) to 0.0d0 for these elements.  
c  
c        hr and hi contain the real and imaginary parts,  
c          respectively, of the complex upper hessenberg matrix.  
c          their lower triangles below the subdiagonal contain further  
c          information about the transformations which were used in the  
c          reduction by  corth, if performed.  if the eigenvectors of  
c          the hessenberg matrix are desired, these elements may be  
c          arbitrary.  
c  
c     on output  
c  
c        ortr, orti, and the upper hessenberg portions of hr and hi  
c          have been destroyed.  
c  
c        wr and wi contain the real and imaginary parts,  
c          respectively, of the eigenvalues.  if an error  
c          exit is made, the eigenvalues should be correct  
c          for indices ierr+1,...,n.  
c  
c        zr and zi contain the real and imaginary parts,  
c          respectively, of the eigenvectors.  the eigenvectors  
c          are unnormalized.  if an error exit is made, none of  
c          the eigenvectors has been found.  
c  
c        ierr is set to  
c          zero       for normal return,  
c          j          if the limit of 30*n iterations is exhausted  
c                     while the j-th eigenvalue is being sought.  
c  
c     calls cdiv for complex division.  
c     calls csroot for complex square root.  
c     calls pythag for  dsqrt(a*a + b*b) .  
c  
c     questions and comments should be directed to burton s. garbow,  
c     mathematics and computer science div, argonne national laboratory  
c  
c     this version dated october 1989.  
c  
c     ------------------------------------------------------------------  
c  
      ierr = 0  
c     .......... initialize eigenvector matrix ..........  
      do 101 j = 1, n  
c  
         do 100 i = 1, n  
            zr(i,j) = 0.0d0  
            zi(i,j) = 0.0d0  
  100    continue  
         zr(j,j) = 1.0d0  
  101 continue  
c     .......... form the matrix of accumulated transformations  
c                from the information left by corth ..........  
      iend = igh - low - 1  
      if (iend) 180, 150, 105  
c     .......... for i=igh-1 step -1 until low+1 do -- ..........  
  105 do 140 ii = 1, iend  
         i = igh - ii  
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140  
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140  
c     .......... norm below is negative of h formed in corth ..........  
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)  
         ip1 = i + 1  
c  
         do 110 k = ip1, igh  
            ortr(k) = hr(k,i-1)  
            orti(k) = hi(k,i-1)  
  110    continue  
c  
         do 130 j = i, igh  
            sr = 0.0d0  
            si = 0.0d0  
c  
            do 115 k = i, igh  
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)  
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)  
  115       continue  
c  
            sr = sr / norm  
            si = si / norm  
c  
            do 120 k = i, igh  
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)  
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)  
  120       continue  
c  
  130    continue  
c  
  140 continue  
c     .......... create real subdiagonal elements ..........  
  150 l = low + 1  
c  
      do 170 i = l, igh  
         ll = min0(i+1,igh)  
         if (hi(i,i-1) .eq. 0.0d0) go to 170  
         norm = pythag(hr(i,i-1),hi(i,i-1))  
         yr = hr(i,i-1) / norm  
         yi = hi(i,i-1) / norm  
         hr(i,i-1) = norm  
         hi(i,i-1) = 0.0d0  
c  
         do 155 j = i, n  
            si = yr * hi(i,j) - yi * hr(i,j)  
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)  
            hi(i,j) = si  
  155    continue  
c  
         do 160 j = 1, ll  
            si = yr * hi(j,i) + yi * hr(j,i)  
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)  
            hi(j,i) = si  
  160    continue  
c  
         do 165 j = low, igh  
            si = yr * zi(j,i) + yi * zr(j,i)  
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)  
            zi(j,i) = si  
  165    continue  
c  
  170 continue  
c     .......... store roots isolated by cbal ..........  
  180 do 200 i = 1, n  
         if (i .ge. low .and. i .le. igh) go to 200  
         wr(i) = hr(i,i)  
         wi(i) = hi(i,i)  
  200 continue  
c  
      en = igh  
      tr = 0.0d0  
      ti = 0.0d0  
      itn = 30*n  
c     .......... search for next eigenvalue ..........  
  220 if (en .lt. low) go to 680  
      its = 0  
      enm1 = en - 1  
c     .......... look for single small sub-diagonal element  
c                for l=en step -1 until low do -- ..........  
  240 do 260 ll = low, en  
         l = en + low - ll  
         if (l .eq. low) go to 300  
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))  
     x            + dabs(hr(l,l)) + dabs(hi(l,l))  
         tst2 = tst1 + dabs(hr(l,l-1))  
         if (tst2 .eq. tst1) go to 300  
  260 continue  
c     .......... form shift ..........  
  300 if (l .eq. en) go to 660  
      if (itn .eq. 0) go to 1000  
      if (its .eq. 10 .or. its .eq. 20) go to 320  
      sr = hr(en,en)  
      si = hi(en,en)  
      xr = hr(enm1,en) * hr(en,enm1)  
      xi = hi(enm1,en) * hr(en,enm1)  
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340  
      yr = (hr(enm1,enm1) - sr) / 2.0d0  
      yi = (hi(enm1,enm1) - si) / 2.0d0  
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)  
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310  
      zzr = -zzr  
      zzi = -zzi  
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)  
      sr = sr - xr  
      si = si - xi  
      go to 340  
c     .......... form exceptional shift ..........  
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))  
      si = 0.0d0  
c  
  340 do 360 i = low, en  
         hr(i,i) = hr(i,i) - sr  
         hi(i,i) = hi(i,i) - si  
  360 continue  
c  
      tr = tr + sr  
      ti = ti + si  
      its = its + 1  
      itn = itn - 1  
c     .......... reduce to triangle (rows) ..........  
      lp1 = l + 1  
c  
      do 500 i = lp1, en  
         sr = hr(i,i-1)  
         hr(i,i-1) = 0.0d0  
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)  
         xr = hr(i-1,i-1) / norm  
         wr(i-1) = xr  
         xi = hi(i-1,i-1) / norm  
         wi(i-1) = xi  
         hr(i-1,i-1) = norm  
         hi(i-1,i-1) = 0.0d0  
         hi(i,i-1) = sr / norm  
c  
         do 490 j = i, n  
            yr = hr(i-1,j)  
            yi = hi(i-1,j)  
            zzr = hr(i,j)  
            zzi = hi(i,j)  
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr  
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi  
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr  
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi  
  490    continue  
c  
  500 continue  
c  
      si = hi(en,en)  
      if (si .eq. 0.0d0) go to 540  
      norm = pythag(hr(en,en),si)  
      sr = hr(en,en) / norm  
      si = si / norm  
      hr(en,en) = norm  
      hi(en,en) = 0.0d0  
      if (en .eq. n) go to 540  
      ip1 = en + 1  
c  
      do 520 j = ip1, n  
         yr = hr(en,j)  
         yi = hi(en,j)  
         hr(en,j) = sr * yr + si * yi  
         hi(en,j) = sr * yi - si * yr  
  520 continue  
c     .......... inverse operation (columns) ..........  
  540 do 600 j = lp1, en  
         xr = wr(j-1)  
         xi = wi(j-1)  
c  
         do 580 i = 1, j  
            yr = hr(i,j-1)  
            yi = 0.0d0  
            zzr = hr(i,j)  
            zzi = hi(i,j)  
            if (i .eq. j) go to 560  
            yi = hi(i,j-1)  
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi  
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr  
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr  
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi  
  580    continue  
c  
         do 590 i = low, igh  
            yr = zr(i,j-1)  
            yi = zi(i,j-1)  
            zzr = zr(i,j)  
            zzi = zi(i,j)  
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr  
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi  
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr  
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi  
  590    continue  
c  
  600 continue  
c  
      if (si .eq. 0.0d0) go to 240  
c  
      do 630 i = 1, en  
         yr = hr(i,en)  
         yi = hi(i,en)  
         hr(i,en) = sr * yr - si * yi  
         hi(i,en) = sr * yi + si * yr  
  630 continue  
c  
      do 640 i = low, igh  
         yr = zr(i,en)  
         yi = zi(i,en)  
         zr(i,en) = sr * yr - si * yi  
         zi(i,en) = sr * yi + si * yr  
  640 continue  
c  
      go to 240  
c     .......... a root found ..........  
  660 hr(en,en) = hr(en,en) + tr  
      wr(en) = hr(en,en)  
      hi(en,en) = hi(en,en) + ti  
      wi(en) = hi(en,en)  
      en = enm1  
      go to 220  
c     .......... all roots found.  backsubstitute to find  
c                vectors of upper triangular form ..........  
  680 norm = 0.0d0  
c  
      do 720 i = 1, n  
c  
         do 720 j = i, n  
            tr = dabs(hr(i,j)) + dabs(hi(i,j))  
            if (tr .gt. norm) norm = tr  
  720 continue  
c  
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001  
c     .......... for en=n step -1 until 2 do -- ..........  
      do 800 nn = 2, n  
         en = n + 2 - nn  
         xr = wr(en)  
         xi = wi(en)  
         hr(en,en) = 1.0d0  
         hi(en,en) = 0.0d0  
         enm1 = en - 1  
c     .......... for i=en-1 step -1 until 1 do -- ..........  
         do 780 ii = 1, enm1  
            i = en - ii  
            zzr = 0.0d0  
            zzi = 0.0d0  
            ip1 = i + 1  
c  
            do 740 j = ip1, en  
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)  
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)  
  740       continue  
c  
            yr = xr - wr(i)  
            yi = xi - wi(i)  
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765  
               tst1 = norm  
               yr = tst1  
  760          yr = 0.01d0 * yr  
               tst2 = norm + yr  
               if (tst2 .gt. tst1) go to 760  
  765       continue  
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))  
c     .......... overflow control ..........  
            tr = dabs(hr(i,en)) + dabs(hi(i,en))  
            if (tr .eq. 0.0d0) go to 780  
            tst1 = tr  
            tst2 = tst1 + 1.0d0/tst1  
            if (tst2 .gt. tst1) go to 780  
            do 770 j = i, en  
               hr(j,en) = hr(j,en)/tr  
               hi(j,en) = hi(j,en)/tr  
  770       continue  
c  
  780    continue  
c  
  800 continue  
c     .......... end backsubstitution ..........  
c     .......... vectors of isolated roots ..........  
      do  840 i = 1, N  
         if (i .ge. low .and. i .le. igh) go to 840  
c  
         do 820 j = I, n  
            zr(i,j) = hr(i,j)  
            zi(i,j) = hi(i,j)  
  820    continue  
c  
  840 continue  
c     .......... multiply by transformation matrix to give  
c                vectors of original full matrix.  
c                for j=n step -1 until low do -- ..........  
      do 880 jj = low, N  
         j = n + low - jj  
         m = min0(j,igh)  
c  
         do 880 i = low, igh  
            zzr = 0.0d0  
            zzi = 0.0d0  
c  
            do 860 k = low, m  
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)  
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)  
  860       continue  
c  
            zr(i,j) = zzr  
            zi(i,j) = zzi  
  880 continue  
c  
      go to 1001  
c     .......... set error -- all eigenvalues have not  
c                converged after 30*n iterations ..........  
 1000 ierr = en  
 1001 return  
      end  
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)  
c  
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr  
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)  
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,  
     x       pythag  
c  
c     this subroutine is a translation of a unitary analogue of the  
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin  
c     and wilkinson.  
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).  
c     the unitary analogue substitutes the qr algorithm of francis  
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.  
c  
c     this subroutine finds the eigenvalues of a complex  
c     upper hessenberg matrix by the qr method.  
c  
c     on input  
c  
c        nm must be set to the row dimension of two-dimensional  
c          array parameters as declared in the calling program  
c          dimension statement.  
c  
c        n is the order of the matrix.  
c  
c        low and igh are integers determined by the balancing  
c          subroutine  cbal.  if  cbal  has not been used,  
c          set low=1, igh=n.  
c  
c        hr and hi contain the real and imaginary parts,  
c          respectively, of the complex upper hessenberg matrix.  
c          their lower triangles below the subdiagonal contain  
c          information about the unitary transformations used in  
c          the reduction by  corth, if performed.  
c  
c     on output  
c  
c        the upper hessenberg portions of hr and hi have been  
c          destroyed.  therefore, they must be saved before  
c          calling  comqr  if subsequent calculation of  
c          eigenvectors is to be performed.  
c  
c        wr and wi contain the real and imaginary parts,  
c          respectively, of the eigenvalues.  if an error  
c          exit is made, the eigenvalues should be correct  
c          for indices ierr+1,...,n.  
c  
c        ierr is set to  
c          zero       for normal return,  
c          j          if the limit of 30*n iterations is exhausted  
c                     while the j-th eigenvalue is being sought.  
c  
c     calls cdiv for complex division.  
c     calls csroot for complex square root.  
c     calls pythag for  dsqrt(a*a + b*b) .  
c  
c     questions and comments should be directed to burton s. garbow,  
c     mathematics and computer science div, argonne national laboratory  
c  
c     this version dated august 1983.  
c  
c     ------------------------------------------------------------------  
c  
      ierr = 0  
      if (low .eq. igh) go to 180  
c     .......... create real subdiagonal elements ..........  
      l = low + 1  
c  
      do 170 i = l, igh  
         ll = min0(i+1,igh)  
         if (hi(i,i-1) .eq. 0.0d0) go to 170  
         norm = pythag(hr(i,i-1),hi(i,i-1))  
         yr = hr(i,i-1) / norm  
         yi = hi(i,i-1) / norm  
         hr(i,i-1) = norm  
         hi(i,i-1) = 0.0d0  
c  
         do 155 j = i, igh  
            si = yr * hi(i,j) - yi * hr(i,j)  
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)  
            hi(i,j) = si  
  155    continue  
c  
         do 160 j = low, ll  
            si = yr * hi(j,i) + yi * hr(j,i)  
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)  
            hi(j,i) = si  
  160    continue  
c  
  170 continue  
c     .......... store roots isolated by cbal ..........  
  180 do 200 i = 1, n  
         if (i .ge. low .and. i .le. igh) go to 200  
         wr(i) = hr(i,i)  
         wi(i) = hi(i,i)  
  200 continue  
c  
      en = igh  
      tr = 0.0d0  
      ti = 0.0d0  
      itn = 30*n  
c     .......... search for next eigenvalue ..........  
  220 if (en .lt. low) go to 1001  
      its = 0  
      enm1 = en - 1  
c     .......... look for single small sub-diagonal element  
c                for l=en step -1 until low d0 -- ..........  
  240 do 260 ll = low, en  
         l = en + low - ll  
         if (l .eq. low) go to 300  
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))  
     x            + dabs(hr(l,l)) + dabs(hi(l,l))  
         tst2 = tst1 + dabs(hr(l,l-1))  
         if (tst2 .eq. tst1) go to 300  
  260 continue  
c     .......... form shift ..........  
  300 if (l .eq. en) go to 660  
      if (itn .eq. 0) go to 1000  
      if (its .eq. 10 .or. its .eq. 20) go to 320  
      sr = hr(en,en)  
      si = hi(en,en)  
      xr = hr(enm1,en) * hr(en,enm1)  
      xi = hi(enm1,en) * hr(en,enm1)  
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340  
      yr = (hr(enm1,enm1) - sr) / 2.0d0  
      yi = (hi(enm1,enm1) - si) / 2.0d0  
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)  
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310  
      zzr = -zzr  
      zzi = -zzi  
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)  
      sr = sr - xr  
      si = si - xi  
      go to 340  
c     .......... form exceptional shift ..........  
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))  
      si = 0.0d0  
c  
  340 do 360 i = low, en  
         hr(i,i) = hr(i,i) - sr  
         hi(i,i) = hi(i,i) - si  
  360 continue  
c  
      tr = tr + sr  
      ti = ti + si  
      its = its + 1  
      itn = itn - 1  
c     .......... reduce to triangle (rows) ..........  
      lp1 = l + 1  
c  
      do 500 i = lp1, en  
         sr = hr(i,i-1)  
         hr(i,i-1) = 0.0d0  
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)  
         xr = hr(i-1,i-1) / norm  
         wr(i-1) = xr  
         xi = hi(i-1,i-1) / norm  
         wi(i-1) = xi  
         hr(i-1,i-1) = norm  
         hi(i-1,i-1) = 0.0d0  
         hi(i,i-1) = sr / norm  
c  
         do 490 j = i, en  
            yr = hr(i-1,j)  
            yi = hi(i-1,j)  
            zzr = hr(i,j)  
            zzi = hi(i,j)  
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr  
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi  
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr  
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi  
  490    continue  
c  
  500 continue  
c  
      si = hi(en,en)  
      if (si .eq. 0.0d0) go to 540  
      norm = pythag(hr(en,en),si)  
      sr = hr(en,en) / norm  
      si = si / norm  
      hr(en,en) = norm  
      hi(en,en) = 0.0d0  
c     .......... inverse operation (columns) ..........  
  540 do 600 j = lp1, en  
         xr = wr(j-1)  
         xi = wi(j-1)  
c  
         do 580 i = l, j  
            yr = hr(i,j-1)  
            yi = 0.0d0  
            zzr = hr(i,j)  
            zzi = hi(i,j)  
            if (i .eq. j) go to 560  
            yi = hi(i,j-1)  
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi  
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr  
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr  
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi  
  580    continue  
c  
  600 continue  
c  
      if (si .eq. 0.0d0) go to 240  
c  
      do 630 i = l, en  
         yr = hr(i,en)  
         yi = hi(i,en)  
         hr(i,en) = sr * yr - si * yi  
         hi(i,en) = sr * yi + si * yr  
  630 continue  
c  
      go to 240  
c     .......... a root found ..........  
  660 wr(en) = hr(en,en) + tr  
      wi(en) = hi(en,en) + ti  
      en = enm1  
      go to 220  
c     .......... set error -- all eigenvalues have not  
c                converged after 30*n iterations ..........  
 1000 ierr = en  
 1001 return  
      end  
  
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)  
c  
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low  
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)  
      double precision f,g,h,fi,fr,scale,pythag  
c  
c     this subroutine is a translation of a complex analogue of  
c     the algol procedure orthes, num. math. 12, 349-368(1968)  
c     by martin and wilkinson.  
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).  
c  
c     given a complex general matrix, this subroutine  
c     reduces a submatrix situated in rows and columns  
c     low through igh to upper hessenberg form by  
c     unitary similarity transformations.  
c  
c     on input  
c  
c        nm must be set to the row dimension of two-dimensional  
c          array parameters as declared in the calling program  
c          dimension statement.  
c  
c        n is the order of the matrix.  
c  
c        low and igh are integers determined by the balancing  
c          subroutine  cbal.  if  cbal  has not been used,  
c          set low=1, igh=n.  
c  
c        ar and ai contain the real and imaginary parts,  
c          respectively, of the complex input matrix.  
c  
c     on output  
c  
c        ar and ai contain the real and imaginary parts,  
c          respectively, of the hessenberg matrix.  information  
c          about the unitary transformations used in the reduction  
c          is stored in the remaining triangles under the  
c          hessenberg matrix.  
c  
c        ortr and orti contain further information about the  
c          transformations.  only elements low through igh are used.  
c  
c     calls pythag for  dsqrt(a*a + b*b) .  
c  
c     questions and comments should be directed to burton s. garbow,  
c     mathematics and computer science div, argonne national laboratory  
c  
c     this version dated august 1983.  
c  
c     ------------------------------------------------------------------  
c  
      la = igh - 1  
      kp1 = low + 1  
      if (la .lt. kp1) go to 200  
c  
      do 180 m = kp1, la  
         h = 0.0d0  
         ortr(m) = 0.0d0  
         orti(m) = 0.0d0  
         scale = 0.0d0  
c     .......... scale column (algol tol then not needed) ..........  
         do 90 i = m, igh  
   90    scale = scale + dabs(ar(i,m-1)) + dabs(ai(i,m-1))  
c  
         if (scale .eq. 0.0d0) go to 180  
         mp = m + igh  
c     .......... for i=igh step -1 until m do -- ..........  
         do 100 ii = m, igh  
            i = mp - ii  
            ortr(i) = ar(i,m-1) / scale  
            orti(i) = ai(i,m-1) / scale  
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)  
  100    continue  
c  
         g = dsqrt(h)  
         f = pythag(ortr(m),orti(m))  
         if (f .eq. 0.0d0) go to 103  
         h = h + f * g  
         g = g / f  
         ortr(m) = (1.0d0 + g) * ortr(m)  
         orti(m) = (1.0d0 + g) * orti(m)  
         go to 105  
c  
  103    ortr(m) = g  
         ar(m,m-1) = scale  
c     .......... form (i-(u*ut)/h) * a ..........  
  105    do 130 j = m, n  
            fr = 0.0d0  
            fi = 0.0d0  
c     .......... for i=igh step -1 until m do -- ..........  
            do 110 ii = m, igh  
               i = mp - ii  
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)  
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)  
  110       continue  
c  
            fr = fr / h  
            fi = fi / h  
c  
            do 120 i = m, igh  
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)  
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)  
  120       continue  
c  
  130    continue  
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........  
         do 160 i = 1, igh  
            fr = 0.0d0  
            fi = 0.0d0  
c     .......... for j=igh step -1 until m do -- ..........  
            do 140 jj = m, igh  
               j = mp - jj  
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)  
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)  
  140       continue  
c  
            fr = fr / h  
            fi = fi / h  
c  
            do 150 j = m, igh  
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)  
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)  
  150       continue  
c  
  160    continue  
c  
         ortr(m) = scale * ortr(m)  
         orti(m) = scale * orti(m)  
         ar(m,m-1) = -g * ar(m,m-1)  
         ai(m,m-1) = -g * ai(m,m-1)  
  180 continue  
c  
  200 return  
      end  
      subroutine cbal(nm,n,ar,ai,low,igh,scale)  
c  
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc  
      double precision ar(nm,n),ai(nm,n),scale(n)  
      double precision c,f,g,r,s,b2,radix  
      logical noconv  
c  
c     this subroutine is a translation of the algol procedure  
c     cbalance, which is a complex version of balance,  
c     num. math. 13, 293-304(1969) by parlett and reinsch.  
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).  
c  
c     this subroutine balances a complex matrix and isolates  
c     eigenvalues whenever possible.  
c  
c     on input  
c  
c        nm must be set to the row dimension of two-dimensional  
c          array parameters as declared in the calling program  
c          dimension statement.  
c  
c        n is the order of the matrix.  
c  
c        ar and ai contain the real and imaginary parts,  
c          respectively, of the complex matrix to be balanced.  
c  
c     on output  
c  
c        ar and ai contain the real and imaginary parts,  
c          respectively, of the balanced matrix.  
c  
c        low and igh are two integers such that ar(i,j) and ai(i,j)  
c          are equal to zero if  
c           (1) i is greater than j and  
c           (2) j=1,...,low-1 or i=igh+1,...,n.  
c  
c        scale contains information determining the  
c           permutations and scaling factors used.  
c  
c     suppose that the principal submatrix in rows low through igh  
c     has been balanced, that p(j) denotes the index interchanged  
c     with j during the permutation step, and that the elements  
c     of the diagonal matrix used are denoted by d(i,j).  then  
c        scale(j) = p(j),    for j = 1,...,low-1  
c                 = d(j,j)       j = low,...,igh  
c                 = p(j)         j = igh+1,...,n.  
c     the order in which the interchanges are made is n to igh+1,  
c     then 1 to low-1.  
c  
c     note that 1 is returned for igh if igh is zero formally.  
c  
c     the algol procedure exc contained in cbalance appears in  
c     cbal  in line.  (note that the algol roles of identifiers  
c     k,l have been reversed.)  
c  
c     arithmetic is real throughout.  
c  
c     questions and comments should be directed to burton s. garbow,  
c     mathematics and computer science div, argonne national laboratory  
c  
c     this version dated august 1983.  
c  
c     ------------------------------------------------------------------  
c  
      radix = 16.0d0  
c  
      b2 = radix * radix  
      k = 1  
      l = n  
      go to 100  
c     .......... in-line procedure for row and  
c                column exchange ..........  
   20 scale(m) = j  
      if (j .eq. m) go to 50  
c  
      do 30 i = 1, l  
         f = ar(i,j)  
         ar(i,j) = ar(i,m)  
         ar(i,m) = f  
         f = ai(i,j)  
         ai(i,j) = ai(i,m)  
         ai(i,m) = f  
   30 continue  
c  
      do 40 i = k, n  
         f = ar(j,i)  
         ar(j,i) = ar(m,i)  
         ar(m,i) = f  
         f = ai(j,i)  
         ai(j,i) = ai(m,i)  
         ai(m,i) = f  
   40 continue  
c  
   50 go to (80,130), iexc  
c     .......... search for rows isolating an eigenvalue  
c                and push them down ..........  
   80 if (l .eq. 1) go to 280  
      l = l - 1  
c     .......... for j=l step -1 until 1 do -- ..........  
  100 do 120 jj = 1, l  
         j = l + 1 - jj  
c  
         do 110 i = 1, l  
            if (i .eq. j) go to 110  
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120  
  110    continue  
c  
         m = l  
         iexc = 1  
         go to 20  
  120 continue  
c  
      go to 140  
c     .......... search for columns isolating an eigenvalue  
c                and push them left ..........  
  130 k = k + 1  
c  
  140 do 170 j = k, l  
c  
         do 150 i = k, l  
            if (i .eq. j) go to 150  
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170  
  150    continue  
c  
         m = k  
         iexc = 2  
         go to 20  
  170 continue  
c     .......... now balance the submatrix in rows k to l ..........  
      do 180 i = k, l  
  180 scale(i) = 1.0d0  
c     .......... iterative loop for norm reduction ..........  
  190 noconv = .false.  
c  
      do 270 i = k, l  
         c = 0.0d0  
         r = 0.0d0  
c  
         do 200 j = k, l  
            if (j .eq. i) go to 200  
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))  
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))  
  200    continue  
c     .......... guard against zero c or r due to underflow ..........  
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270  
         g = r / radix  
         f = 1.0d0  
         s = c + r  
  210    if (c .ge. g) go to 220  
         f = f * radix  
         c = c * b2  
         go to 210  
  220    g = r * radix  
  230    if (c .lt. g) go to 240  
         f = f / radix  
         c = c / b2  
         go to 230  
c     .......... now balance ..........  
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270  
         g = 1.0d0 / f  
         scale(i) = scale(i) * f  
         noconv = .true.  
c  
         do 250 j = k, n  
            ar(i,j) = ar(i,j) * g  
            ai(i,j) = ai(i,j) * g  
  250    continue  
c  
         do 260 j = 1, l  
            ar(j,i) = ar(j,i) * f  
            ai(j,i) = ai(j,i) * f  
  260    continue  
c  
  270 continue  
c  
      if (noconv) go to 190  
c  
  280 low = k  
      igh = l  
      return  
      end  
      subroutine cdiv(ar,ai,br,bi,cr,ci)  
      double precision ar,ai,br,bi,cr,ci  
c  
c     complex division, (cr,ci) = (ar,ai)/(br,bi)  
c  
      double precision s,ars,ais,brs,bis  
      s = dabs(br) + dabs(bi)  
      ars = ar/s  
      ais = ai/s  
      brs = br/s  
      bis = bi/s  
      s = brs**2 + bis**2  
      cr = (ars*brs + ais*bis)/s  
      ci = (ais*brs - ars*bis)/s  
      return  
      end  
  
      subroutine csroot(xr,xi,yr,yi)  
      double precision xr,xi,yr,yi  
c  
c     (yr,yi) = complex dsqrt(xr,xi)  
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)  
c  
      double precision s,tr,ti,pythag  
      tr = xr  
      ti = xi  
      s = dsqrt(0.5d0*(pythag(tr,ti) + dabs(tr)))  
      if (tr .ge. 0.0d0) yr = s  
      if (ti .lt. 0.0d0) s = -s  
      if (tr .le. 0.0d0) yi = s  
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)  
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)  
      return  
      end  
  
      double precision function pythag(a,b)  
      double precision a,b  
c  
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow  
c  
      double precision p,r,s,t,u  
      p = dmax1(dabs(a),dabs(b))  
      if (p .eq. 0.0d0) go to 20  
      r = (dmin1(dabs(a),dabs(b))/p)**2  
   10 continue  
         t = 4.0d0 + r  
         if (t .eq. 4.0d0) go to 20  
         s = r/t  
         u = 1.0d0 + 2.0d0*s  
         p = u*p  
         r = (s/u)**2 * r  
      go to 10  
   20 pythag = p  
      return  
      end  
c  
c  
c  
  
      subroutine quanc8(fun,a,b,abserr,relerr,result,errest,nofun,flag)
c
      double precision fun, a, b, abserr, relerr, result, errest, flag
      integer nofun
c
c   estimate the integral of fun(x) from a to b
c   to a user provided tolerance.
c   an automatic adaptive routine based on
c   the 8-panel newton-cotes rule.
c
c   input ..
c
c   fun     the name of the integrand function subprogram fun(x).
c   a       the lower limit of integration.
c   b       the upper limit of integration.(b may be less than a.)
c   relerr  a relative error tolerance. (should be non-negative)
c   abserr  an absolute error tolerance. (should be non-negative)
c
c   output ..
c
c   result  an approximation to the integral hopefully satisfying the
c           least stringent of the two error tolerances.
c   errest  an estimate of the magnitude of the actual error.
c   nofun   the number of function values used in calculation of result.
c   flag    a reliability indicator.  if flag is zero, then result
c           probably satisfies the error tolerance.  if flag is
c           xxx.yyy , then  xxx = the number of intervals which have
c           not converged and  0.yyy = the fraction of the interval
c           left to do when the limit on  nofun  was approached.
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	SUBROUTINE LEST (N,A1,B1,MS,MF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(500),Y(500),C(500),A(10,11),XN(500)
	INTEGER N
	COMMON /SQUARE/ X,Y 
      IF (MF.GT.(N-1)) THEN   
         MF=N-1            
         PRINT 200, MF
      ENDIF  
    5 MFP1=MF+1
      MFP2=MF+2
C
      DO 10 I=1,N        
         XN(I)=1.D0
   10 CONTINUE
C
      DO 30 I=1,MFP1   
         A(I,1)=0.D0      
         A(I,MFP2)=0.D0     
         DO 20 J=1,N        
            A(I,1)=A(I,1)+XN(J) 
            A(I,MFP2)=A(I,MFP2)+Y(J)*XN(J)      
            XN(J)=XN(J)*X(J)     
   20    CONTINUE     
   30 CONTINUE       
C   
      DO 50 I=2,MFP1
         A(MFP1,I)=0.D0       
         DO 40 J=1,N       
            A(MFP1,I)=A(MFP1,I)+XN(J)       
            XN(J)=XN(J)*X(J)              
   40    CONTINUE
   50 CONTINUE
C
      DO 70 J=2,MFP1
         DO 60 I=1,MF
            A(I,J)=A(I+1,J-1)
   60    CONTINUE
   70 CONTINUE
C
C      PRINT '(///)'
C      PRINT*, '          THE NORMAL EQUATIONS ARE :'
C      PRINT'(/)'
C      PRINT 201, ((A(I,J),J=1,MFP2),I=1,MFP1) 
C      PRINT '(//)'
C
      CALL LUDCMQ(A,MFP1,10)
C
      MSP1=MS+1
      DO 95 I=MSP1,MFP1
         DO 90 J=1,I
               C(J)=A(J,MFP2)
   90    CONTINUE
         CALL SOLNQ (A,C,I,10)
         IM1=I-1
	   A1=C(1)
	   B1=C(2)
         BETA=0.0
         DO 94 IPT=1,N
               SUM=0.D0
               DO 93 ICOEF=2,I
                     JCOEF=I-ICOEF+2
                     SUM=(SUM+C(JCOEF))*X(IPT)
   93          CONTINUE
               SUM=SUM+C(1)
               BETA=BETA+(Y(IPT)-SUM)**2
   94    CONTINUE
         BETA=BETA/(N-1)
C         PRINT 203, BETA
   95 CONTINUE
C
  200 FORMAT(//' DEGREE OF POLNOMIAL CAN NOT EXCEED N-1.',/
     +      '  REQUESTED MAXIMUM DEGREE TOO LARGE -',
     +      'REDUCED TO ',I3)
C  201 FORMAT(1X,9F8.2)
C  202 FORMAT(/' FOR DEGREE OF ',I2,' COEFFICIENTS ARE'//
C     +      ' ',5X,11F9.3)
C  203 FORMAT(9X,' BETA IS ',F10.5//)
C      STOP
      RETURN
      END
C
      SUBROUTINE LUDCMQ (A,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,NDIM),SUM
      DO 30 I=1,N
         DO 30 J=2,N
            SUM=0.D0
            IF (J.LE.I) THEN
               JM1=J-1
               DO 10 K =1,JM1
                  SUM=SUM+A(I,K)*A(K,J)
   10          CONTINUE
               A(I,J)=A(I,J)-SUM
            ELSE
               IM1=I-1
               IF (IM1.NE.0) THEN
                  DO 20 K=1,IM1
                     SUM=SUM+A(I,K)*A(K,J)
   20             CONTINUE
               ENDIF
   25          IF (DABS(A(I,I)).LT.1.0D-10) THEN
                  PRINT 100,I
                  RETURN
               ELSE
               A(I,J)=(A(I,J)-SUM)/A(I,I)
               ENDIF
            ENDIF
   30 CONTINUE
      RETURN                                   
  100 FORMAT(' REDUCTION NOT COMPLICATED BECAUSE SMALL VALUE',
     +       ' FOUND FOR DIVISOR IN ROW ',I3)
      END
C
      SUBROUTINE SOLNQ (A,B,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,NDIM),B(NDIM)
      B(1)=B(1)/A(1,1)
      DO 20 I=2,N
         IM1=I-1
         SUM=0.D0
         DO 10 K=1,IM1
            SUM=SUM+A(I,K)*B(K)
   10    CONTINUE
         B(I)=(B(I)-SUM)/A(I,I)
   20 CONTINUE
C
      DO 40 J=2,N
         NMJP2=N-J+2
         NMJP1=N-J+1
         SUM=0.D0
         DO 30 K=NMJP2,N
            SUM=SUM+A(NMJP1,K)*B(K)
   30    CONTINUE
         B(NMJP1)=B(NMJP1)-SUM
   40 CONTINUE
      RETURN
      END                                                                         
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                          
      SUBROUTINE ZLOC(ARRAY,APZOO1,APZOO2)
C
C **********************************************************************
C
C     SUBROUTINE ZLOC
C
C     PURPOSE
C  LOCATING THE ZEROS OF A FUNCTION OF
C  ONE VARIABLE IN AN INTERVAL
C
C     USAGE
C  CALL ZLOC(ARRAY,SUB1,SUB2)
C
C     DESCRIPTION
C  ARRAY - ARRAY OF LENGTH 5
C  ARRAY(1) - ENDPOINT OF SEARCH INTERVAL
C  ARRAY(2) - OTHER ENDPOINT
C  ARRAY(3) - TENTATIVE STEP SIZE
C  ARRAY(4) - ZERO COUNTER
C  ARRAY(5) - TERMINATION SWITCH
C  SUB1 - USER SUPPLIED SUBROUTINE SUBR1(X,Y)
C  ZLOC SUPPLIES ARGUMENT X AND USER STORES
C  CORRESPONDING FUNCTION VALUE IN Y
C  SUBR2 - USER SUPPLIED SUBROUTINE SUBR2(Z,E)
C  CURRENT ZERO Z AND ASSOCIATED ABSOLUTE
C  ERROR E IS SUPPLIED TO USER
C
C     REMARKS
C  1.  ARRAY(4) CONTAINS AS A FLOATING NUMBER THE
C  CURRENT NUMBER OF ZEROS FOUND
C  2.  USER MAY END RUN AT ANY TIME PRIOR TO NORMAL
C      TERMINATION BY MAKING ARRAY(5) DIFFERENT FROM 0.
C  3.  SEARCH PROCEEDS FROM ARRAY(1) TO ARRAY(2)
C      REGARDLESS FO THE ORDERING OF THESE ENTRIES
C  4.  SUB1 AND SUB2 MUST BE NAMED IN AN EXTERNAL
C      STATEMENT IN THE CALLING PROGRAM
C  5.  AN ERROR STOP WITH AN APPROPRIATE MESSAGE OCCURS IS
C      A - THE SEARCH INTERVAL BECOMES TOO SMALL
C      B - 30 SUCCESSIVE INTERVAL HALVINGS OCCUR
C      C - THE ENDPOINTS ARE INITIALLY EQUAL
C  6.  A VALID DOUBLE PRECISION ROUTINE RESULTS IF
C      3 REAL DECLARATIONS AT BEGINNING OF PROGRAM
C      ARE CHANGED TO REAL*8 AND THE COMMENT DESIGNATOR
C      C IS DELETED FROM STATEMENT FOLLOWING
C      FIRST EXECUTABLE STATEMENT
C
C     METHOD
C  MUELLER,S METHOD COUPLED WITH ADAPTIVE SEARCH
C  BASED ON QUADRATIC EXTRAPOLATION
C
C **********************************************************************
C      
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ND,NC,J,INST,KRT,IM,IP
      REAL*8 X(4),Y(4),S(4),T(4),F(3),U(3),ARRAY(5),A,AH,ALPHA
     .     ,B,BETA,BH,CH,DEL,DH,DLAMB,E,H,HPR,Q0,Q1,Q2
     .     ,R,RHO,TWO,W,XP,YG,YP,Z
      E=1.D-14
      TWO=2.D0
      ND=0
      NC=0
      ARRAY(4)=0.D0
      ARRAY(5)=0.D0
      A=ARRAY(1)
      B=ARRAY(2)
      H=ARRAY(3)
      IF(H+A.NE.A)GOTO2
    1 H=(B-A)/4.D0
      GOTO3
    2 H=DABS(H)
      IF(H.GT.DABS(B-A)/4.D0)H=DABS(B-A)/4.D0
      IF(B.LT.A)H=-H
    3 IF(H)4,12,4
   12 WRITE(6,13)
   13 FORMAT(' SUBR ZLOC CALLED WITH = ENDPOINTS')
  511 STOP
  507 WRITE(6,508)
  508 FORMAT(' INTERVAL TOO SMALL IN SUBR ZLOC')
      GO TO 511
    4 Z=A
      ALPHA= A
      BETA=B
      IF(B.LT.A)ALPHA=B
      IF(A.GT.B)BETA=A
      J=1
  700 CALL APZOO1((Z),W)
      IF(ARRAY(5).NE.0.D0)GOTO706
      X(J)=Z
      Y(J)=W
      Z=Z+H
      IF(ALPHA.GT.Z)Z=ALPHA
      IF(BETA.LT.Z)Z=BETA
      IF(Z-X(J))5,507,5
    5 H=Z-X(J)
      J=J+1
      IF(J-4)700,700,6
    6 YG=Y(1)+3.*(Y(3)-Y(2))
      DH=DABS(Y(4))
      IF(DABS(Y(3)).GT.DH)DH=DABS(Y(3))
      IF(DABS(Y(2)).GT.DH)DH=DABS(Y(2))
      IF(DABS(Y(1)).GT.DH)DH=DABS(Y(1))
      DH=3.D0*DH
      IF(DABS(Y(4)-YG)-.0333D0*DH)11,11,7
    7 NC=NC+1
      IF(NC-30)9,8,8
    8 WRITE(6,107)
  107 FORMAT(' 30 HALVES IN SUBR ZLOC')
      GO TO 511
    9 H=H/TWO
      X(3)=X(2)
      Y(3)=Y(2)
      Z=A+H
  509 J=2
  701 CALL APZOO1((Z),W)
      IF(ARRAY(5).NE.0.D0)GOTO706
      X(J)=Z
      Y(J)=W
   10 Z=Z+H+H
      J=J+2
      IF(J-4)701,701,5022
 5022 IF(X(1).EQ.X(2))GOTO507
      IF(X(2).EQ.X(3))GOTO507
      IF(X(3).EQ.X(4))GOTO507
      GOTO6
   11 NC=0
      INST=1
      I=1
      Z=X(1)
      W=Y(1)
      KRT=-1
      IF(W)300,250,300
  300 KRT=0
      Z=X(I+1)
      W=Y(I+1)
      IF(W)301,250,301
  301 KRT=1
  302 IF(Y(I).EQ.0.D0)GOTO303
      IF(Y(I)/DABS(Y(I))*Y(I+1).GE.0.D0)GOTO303
      DLAMB=X(I)
      DEL=Y(I)
      RHO=X(I+1)
      R=Y(I+1)
      S(1)=X(1)
      S(2)=X(2)
      S(3)=X(3)
      T(1)=Y(1)
      T(2)=Y(2)
      T(3)=Y(3)
      GOTO200
  303 IF(I.EQ.2)GOTO304
      I=2
      GOTO300
  304 IF(DABS(Y(1))+DABS(Y(2)).EQ.0.D0)GOTO350
      IF(DABS(Y(2))+DABS(Y(3)).EQ.0.D0)GOTO350
      IF(Y(1)/(Y(2)+Y(3)).LT.0.)GOTO350
      IF(DABS(Y(3))+DABS(Y(1)).EQ.0.D0)GOTO350
      U(1)=X(1)
      U(2)=X(2)
      U(3)=X(3)
      F(1)=Y(1)
      F(2)=Y(2)
      F(3)=Y(3)
  310 ND=2
      IF(F(3).NE.0.D0)GOTO311
      ND=1
      IF(U(3).NE.B)GOTO350
      GOTO315
  311 IF(F(1).NE.0.D0)GOTO312
      ND=1
      IF(U(1).NE.A)GOTO350
      GOTO315
  312 IF(DABS(F(1)).GT.DABS(F(2)))GOTO313
      IF(U(1).NE.A)GOTO350
      GOTO314
  313 IF(DABS(F(3)).GT.DABS(F(2)))GOTO314
      IF(U(3).NE.B)GOTO350
  314 IF(F(2).EQ.0.D0)ND=1
      IF(F(1)/F(3).LT.0.D0)GOTO350
  315 HPR=U(2)
      CH=F(2)
 3150 IF(DABS(F(2)).LT.DABS(CH))HPR=U(2)
      IF(DABS(F(2)).LT.DABS(CH))CH=F(2)
      IF((HPR-U(1))*(U(3)-HPR).LT.0.D0)GOTO350
 3151 IM=1
      IP=3
      IF(DABS(F(3)).LE.DABS(F(1)))GOTO3152
      IP=1
      IM=3
 3152 IF(DABS(F(2)).GE.DABS(F(IP)))GOTO3154
      Z=U(IP)-U(2)
      IF(DABS(Z).GE.DABS(U(2)-U(IM)))GOTO3154
      Z=Z+U(IM)
      CALL APZOO1((Z),W)
      IF(ARRAY(5).NE.0.D0)GOTO706
      IF((W-F(2))/(F(1)+F(3)).GT.0.D0)GOTO3153
      F(IP)=F(2)
      U(IP)=U(2)
      U(2)=Z
      F(2)=W
      GOTO3151
 3153 IF((W-F(IP))/(F(1)+F(3)).LT.0.D0)GOTO3155
      U(IM)=Z
      F(IM)=W
      GOTO3151
 3155 U(2)=Z
      F(2)=W
      GOTO3151
 3154 T(1)=(F(1)-F(2))/(U(1)-U(2))
      T(3)=(F(2)-F(3))/(U(2)-U(3))
      S(1)=T(1)/(U(3)-U(1))
      S(3)=T(3)/(U(3)-U(1))
      Z=(T(3)-T(1))/(U(3)-U(1))*TWO
      IF((F(1)+F(3))/DABS(F(1)+F(3))*Z.LE.0.D0)GOTO350
      Z=(S(3)*(U(1)-U(2))+S(1)*(U(2)-U(3)))/Z+U(2)
      IF((Z-U(1))*(U(3)-Z).LE.0.D0)GOTO350
      IF(Z.EQ.HPR)GOTO350
      IF(Z.EQ.U(2))GOTO350
      YG=(S(3)*(Z-U(1))-S(1)*(Z-U(3)))*(Z-U(2))+F(2)
      IF(YG/(F(1)+F(3)).GE..05D0)GOTO350
      CALL APZOO1((Z),W)
      IF(ARRAY(5).NE.0.D0)GOTO706
      IF(W.EQ.0.D0)GOTO250
      IF(W/(F(1)+F(3)).GE..05D0)GOTO350
      IF(W/(F(1)+F(3)).GT..0D0)GOTO319
      T(1)=F(1)
      T(2)=F(2)
      T(3)=F(3)
      S(1)=U(1)
      S(2)=U(2)
      S(3)=U(3)
      IF((Z-S(1))*(S(2)-Z).LT.0.D0)GOTO318
      S(3)=S(2)
      T(3)=T(2)
  316 T(2)=W
      S(2)=Z
      F(2)=W
      U(2)=Z
      IF(T(1).EQ.0.)GOTO317
      DLAMB=S(1)
      DEL=T(1)
      RHO=Z
      R=W
      GOTO200
  317 DLAMB=Z
      DEL=W
      RHO=S(3)
      R=T(3)
      GOTO200
  318 S(1)=S(2)
      T(1)=T(2)
      GOTO316
  319 IF(DABS(YG-W).LE..05D0*(DABS(W)+DABS(YG)))GOTO350
      IM=1
      IF((Z-U(1))*(U(2)-Z).LT.0.D0)IM=3
      IF(DABS(W).GT.DABS(F(IM)))GOTO350
      IF(DABS(W).GT.DABS(F(2 )))GOTO320
      IP=4-IM
      IM=2
      F(IP)=F(2)
      U(IP)=U(2)
  320 F(IM)=W
      U(IM)=Z
      GOTO3150
  325 ND=ND-1
      IF(ND.EQ.0D0)GOTO350
      IF(F(2)/(F(1)+F(3)).GT.0.D0)GOTO326
      S(1)=Z
      T(1)=0.D0
      S(2)=U(2)
      T(2)=F(2)
      S(3)=U(3)
      T(3)=F(3)
      Z=U(2)
      W=F(2)
      GOTO317
  326 U(2)=Z
      F(2)=0.D0
      GOTO310
  330 IF(ND.NE.0)GOTO325
      IF(KRT)300,301,331
  331 IF(I.EQ.2)GOTO350
      I=2
      GOTO300
  350 ND=0
      IF(INST.EQ.1)GOTO353
  351 IF(X(3)-B)100,706,100
  706 IF(H)800,800,801
  800 B=ALPHA
      A=BETA
      RETURN
  801 A=ALPHA
      B=BETA
      RETURN
  353 INST=0
      XP=X(1)
      YP=Y(1)
      X(1)=X(2)
      X(2)=X(3)
      X(3)=X(4)
      Y(1)=Y(2)
      Y(2)=Y(3)
      Y(3)=Y(4)
      GOTO300
  100 W=-YP
      X(4)=X(3)
      Z=0.D0
      J=0
  707 J1=J+1-((J+1)/3)*3+1
      J2=J+2-((J+2)/3)*3+1
      Z=Z+DABS(Y(J+1))
  101 W=W+Y(J+1)*(XP-X(J1))/(X(J+1)-X(J1))*(XP-X(J2))/(X(J+1)-X(J2))
      J=J+1
      IF(J-2)707,707,708
  708 H=(X(3)-X(2))*3./TWO
      IF(W)110,104,110
  110 AH=((X(3)-X(2))+(X(3)-X(1)))
      BH=(X(3)-X(2))*(X(3)-X(1))
      CH=.021517*(X(1)-XP)*(X(2)-XP)*(X(3)-XP)*Z/DABS(W)
      CH=CH*.666666666D0
  102 HPR=((2.D0*H+AH)*H*H+CH)/((3.D0*H+2.D0*AH)*H+BH)
      DH=DABS(HPR)-DABS(H)
      IF(DH)103,104,104
  103 H=HPR
      IF(DH+.1*DABS(H))102,104,104
  104 Z=X(3)+H
      IF(ALPHA.GT.Z)Z=ALPHA
      IF(BETA.LT.Z)Z=BETA
      IF(Z-X(3))510,507,510
  510 IF(Z-X(4))450,507,450
  450 H=Z-X(3)
      CALL APZOO1((Z),W)
      IF(ARRAY(5).NE.0.D0)GOTO706
      X(4)=Z
      Y(4)=W
      J=0
  709 J1=J+1-((J+1)/3)*3+1
      J2=J+2-((J+2)/3)*3+1
  105 W=W-Y(J+1)*(Z-X(J1))/(X(J+1)-X(J1))*(Z-X(J2))/(X(J+1)-X(J2))
      J=J+1
      IF(J-2)709,709,710
  710 DH=DABS(Y(4))
      IF(DABS(Y(3)).GT.DH)DH=DABS(Y(3))
      IF(DABS(Y(2)).GT.DH)DH=DABS(Y(2))
      IF(DABS(Y(1)).GT.DH)DH=DABS(Y(1))
      DH=3.D0*DH
      IF(DABS(W)-.0333D0*DH)109,109,106
  106 NC=NC+1
      IF(NC-30)108,8,8
  108 H=(Z-X(3))/TWO
      GOTO104
  109 NC=0
      GOTO353
  200 Q2=T(3)*(S(2)-S(1))+T(1)*(S(3)-S(2))+T(2)*(S(1)-S(3))
      Q1=(T(2)-T(3))*S(1)**2
      Q1=(T(3)-T(1))*S(2)**2+Q1
      Q1=((T(1)-T(2))*S(3)**2+Q1)/TWO
      Q0=(S(2)-S(1))*S(1)*S(2)*T(3)
      Q0=(S(3)-S(2))*S(2)*S(3)*T(1)+Q0
      Q0=(S(1)-S(3))*S(3)*S(1)*T(2)+Q0
      IF(Q2)210,201,210
  201 IF(Q1)202,222,202
  202 Z=Q0/(Q1*TWO)
      GOTO214
  210 IF(Q1)209,223,209
  223 Z=-(Q0/Q2)
      IF(Z)222,225,224
  224 Z=DSQRT(Z)
      GOTO225
  209 Z=DABS(Q1)-(Q0/DABS(Q1))*Q2
      IF(Z)222,211,211
  211 Z=(DABS(Q1)+DSQRT(DABS(Q1))*DSQRT(Z))/Q2
      IF(Q1.LT.0D0)Z=-Z
  225 IF((Z-DLAMB)*(RHO-Z))212,215,215
  212 IF(Z)213,222,213
  213 Z=Q0/(Q2*Z)
  214 IF((Z-DLAMB)*(RHO-Z))222,222,215
  215 S(4)=Z
      CALL APZOO1((Z),W)
      IF(ARRAY(5).NE.0.D0)GOTO706
      T(4)=W
      IF(W)216,250,216
  216 IF(DABS(Z-S(3))-E*DABS(S(1)))250,250,217
  217 Q0=DABS(W)
      IF(DEL.LT.0D0)Q0=-Q0
      IF(Q0-W)218,221,218
  218 RHO=Z
      R=W
  219 S(1)=S(2)
      S(2)=S(3)
      S(3)=S(4)
      T(1)=T(2)
      T(2)=T(3)
      T(3)=T(4)
      GOTO200
  221 DLAMB=Z
      DEL=W
      GOTO219
  222 Z=(R*DLAMB-DEL*RHO)/(R-DEL)
      IF((Z-DLAMB)*(RHO-Z).LT.0.D0)Z=S(3)
      GOTO215
  250 S(1)=X(1)
      T(1)=Y(1)
      S(2)=X(3)
      T(2)=Y(3)
      IF(Z.EQ.X(2))GOTO252
      IF(Z.EQ.X(1))GOTO251
      S(2)=X(2)
      T(2)=Y(2)
      GOTO252
  251 S(1)=X(2)
      T(1)=Y(2)
  252 W=(S(2)-S(1))/TWO
      S(1)=(S(1)-Z)/W
      S(2)=(S(2)-Z)/W
      T(3)=(DABS(T(1))+DABS(T(2)))/W
      IF(T(3).EQ.0.D0)W=((DABS(S(1))+DABS(S(2)))/TWO+1.D0)/TWO*W
      IF(T(3).EQ.0.D0)GOTO253
      T(1)=T(1)/S(1)/T(3)
      T(2)=T(2)/S(2)/T(3)
      S(3)=DABS(T(2)-T(1))*TWO/E
      T(3)=DABS(T(2)*S(1)-T(1)*S(2))/E
      IF(S(3)+T(3).EQ.0.D0)GOTO253
      S(3)=W/(T(3)+DSQRT(DABS(T(3)**2-S(3))))*(DLOG(E)*.05D0-.0346D0)**2
      IF(DABS(S(3)).LT.DABS(W))W=S(3)
  253 IF(DABS(W).LT.DABS(Z)*E/TWO)W=Z*E/TWO
      ARRAY(4)=ARRAY(4)+1.
      CALL APZOO2((Z),DABS(W))
      IF(ARRAY(5).NE.0.)GOTO706
      GOTO330
      END

