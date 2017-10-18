! A very old set of subroutines from the Harwell Subroutine Library (HSL)
! Free for use in non-commercial applications
!
MODULE liboceqplus
!
! modified 2015 to F90 by Bo Sundman for use in OC
!
!use general_thermodynamic_package
!
! problems as assessment_calfun is defined in matsmin
! I did not manage to pass the subroutine name as argument
  use liboceq
!
CONTAINS

  SUBROUTINE VA05AD(M,N,F,X,DSTEP,DMAX,ACC,MAXFUN,IPRINT,W,&
       IENT,IEXIT)
!  SUBROUTINE VA05AD(M,N,F,X,DSTEP,DMAX,ACC,MAXFUN,IPRINT,W,&
!       IENT,IEXIT,ash)
!  SUBROUTINE VA05AD(M,N,F,X,DSTEP,DMAX,ACC,MAXFUN,IPRINT,W,&
!       IENT,IEXIT,CALFUN)
!  SUBROUTINE VA05AD(M,N,F,X,DSTEP,DMAX,ACC,MAXFUN,IPRINT,W,&
!       IENT,IEXIT,ash,calfun)
!  SUBROUTINE VA05AD(M,N,F,X,DSTEP,DMAX,ACC,MAXFUN,IPRINT,W,&
!       IENT,IEXIT,ash)
!...... 20/05/70 LAST LIBRARY UPDATE
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!   NOTE THAT THE INSTRUCTION CALLING SUBROUTINE 'MB11A',
!   ON LINE NUMBER '138' ,IS NOT STANDARD FORTRAN
    INTEGER, PARAMETER :: LP=6
    DIMENSION F(*),X(*),W(*)
    DIMENSION IEXIT(*)
    integer, parameter :: MSAVE=2500
    DIMENSION WSAVE(MSAVE)
!    EXTERNAL CALFUN
! this contains all necessary information to calculate error for experiments
!    type(gtp_assessmenthead), pointer :: ash
!    type(gtp_assessmenthead) :: ash
!...check that all important local variables are saved here!!!!
    SAVE MPN,NT,NTEST,DTEST,NWI,NWX,NWF,NWC,NWD,NWW,NWV,NWT,NWU,&
         FMIN,DD,DSS,DM,PARM,PAR,PPAR,DPAR,IS,IC,TINC,WSAVE
    if(IENT.ne.0) then
! CONTINUE assessment with same Jacobian when IENT=0
!      IF(IEXIT(4).NE.2) THEN
!         WRITE(LP,*)'Sorry, can continue after TOO MANY ITERATIONS only'
!         RETURN
!      ENDIF
       IENT=2
       IPC=0
       MAXC=-1
! RESTORE THE LAST VALUES OF X AND F
       DO I=1,N
          X(I)=WSAVE(I)
          F(I)=WSAVE(N+I)
       enddo
       WRITE(LP,9998)
9998   FORMAT (4X,'The following output is provided by subroutine',&
           ' VA05A'/4X,'Optimization continuing with same Jacobian'/)
       goto 3
    else
! start of VA05AD, SET VARIOUS PARAMETERS
       MAXC=-1
! MAXC COUNTS THE NUMBER OF CALLS OF CALFUN, initiate to -1
       MPN=M+N
       NT=N+2
!       NTEST=0
       NTEST=NT
! 'NT' AND 'NTEST' CAUSE AN ERROR RETURN IF F(X) DOES NOT DECREASE
       DTEST=FLOAT(N+N)-0.5D0
!     'DTEST' IS USED IN A TEST TO MAINTAIN LINEAR INDEPENDENCE
!     PARTITION THE WORKING SPACE ARRAY W
!     THE FIRST PARTITION HOLDS THE JACOBIAN APPROXIMATION
       NWI=M*N
!     THE NEXT PARTITION HOLDS THE GENERALIZED INVERSE
       NWX=NWI+MPN*N
!     THE NEXT PARTITION HOLDS THE BEST VECTOR X
       NWF=NWX+N
!     THE NEXT PARTITION HOLDS THE BEST VECTOR F
       NWC=NWF+M
!     THE NEXT PARTITION HOLDS THE COUNTS OF THE INDEPENDENT DIRECTIONS
       NWD=NWC+N
!     THE NEXT PARTITION HOLDS THE INDEPENDENT DIRECTIONS
       NWW=NWD+N*N
!     THE REMAINDER OF W IS USED FOR SCRATCH VECTORS
       NWV=NWW+N
       NWT=NWV+M
       NWU=NWT+N
! we could allocate W here but then the continue command is impossible
!       write(*,*)'VA05AD workspace needed ',NWF,NWU
       FMIN=-1.0D0
!     USUALLY 'FMIN' IS THE LEAST CALCULATED VALUE OF F(X)
!     DW USED BEFORE SET?
!      DW=0.0D0
       DD=0.0D0
!     USUALLY 'DD' IS THE SQUARE OF THE CURRENT STEP LENGTH
       DSS=DSTEP*DSTEP
       DM=DMAX*DMAX
       PARM=SQRT(ACC)/DMAX
!     'PARM' IS THE LEAST VALUE OF THE MARQUARDT PARAMETER
       DPAR=10.0D0*DM
!     'DPAR' AND 'NTPAR' ARE USED TO REGULATE THE MARQUARDT PARAMETER
       IS=4
!     'IS' CONTROLS A GO TO STATEMENT FOLLOWING A CALL OF CALFUN
       IC=0
       TINC=1.0D0
!     'TINC' IS USED IN THE CRITERION TO INCREASE THE STEP LENGTH
!     START A NEW PAGE FOR PRINTING
       IF (IPRINT) 1,3,1
1      WRITE(LP,2)
2      FORMAT (4X,'The following output is provided by subroutine VA05A'/)
       IPC=0
       GO TO 3
    endif
! =========================================== return here for each iteration
!     TEST WHETHER THERE HAVE BEEN MAXFUN CALLS OF CALFUN
4   IF (MAXFUN-MAXC) 5,5,3
5   CONTINUE
    IEXIT(4)=2
    IF (IPRINT) 139,140,139
140 IPRINT=2
    GO TO 19
139 WRITE(LP,6)MAXC
6   FORMAT (/1X,'*** ','ERROR RETURN FROM VA05A BECAUSE ',&
         'THERE HAVE BEEN',I5,' CALLS OF CALFUN')
    GO TO 7
!     CALL THE SUBROUTINE CALFUN
3   MAXC=MAXC+1
!============= CONDITION FOR USER BREAK
!    IF(KABORT('TURE').EQ.1)THEN
!       write(LP,*)'user break'
!       IEXIT(4)=3
!       RETURN
!    ENDIF
!    write(*,*)'Calling calfun'
!    CALL CALFUN (M,N,F,X)
!    CALL CALFUN (M,N,F,X,ash)
!    CALL old_assessment_calfun (M,N,F,X,ash)
!    CALL assessment_calfun (M,N,F,X)
!    CALL new_assessment_calfun (M,N,F,X,ash)
!    CALL assessment_calfun (M,N,F,X,ash)
! We use "firstash", in CALFUN a global data structure declared in gtp
    CALL assessment_calfun (M,N,F,X)
!     CALCULATE THE SUM OF SQUARES
    FSQ=0.0D0
    DO I=1,M
       FSQ=FSQ+F(I)*F(I)
    enddo
!    8 CONTINUE
!     TEST FOR ERROR RETURN BECAUSE F(X) DOES NOT DECREASE
    GO TO (9,10,9,10),IS
9   IF (FSQ-FMIN) 11,12,12
12  IF (DD-DSS) 13,13,10
13  NTEST=NTEST-1
    IF (NTEST) 14,14,10
14  CONTINUE
    IEXIT(4)=1
    IF (IPRINT) 15,17,15
17  IPRINT=1
    GO TO 19
15  WRITE(LP,16)
16  FORMAT(/1X,'*** ERROR RETURN FROM VA05A BECAUSE F(X) NO LONGER',&
         ' DECREASES'/5X,'THIS MAY BE DUE TO THE VALUES OF DSTEP',&
         ' AND ACC, OR TO LOSS OF RANK IN THE JACOBIAN MATRIX')
!     PROVIDE PRINTING OF FINAL SOLUTION IF REQUESTED
7   IF (IPRINT) 18,19,18
18  WRITE(LP,20)MAXC
20  FORMAT (/5X,'THE FINAL SOLUTION CALCULATED BY VA05A REQUIRED',&
         I5,' iterations')
! best set of variables
    WRITE(LP,21)(I,W(NWX+I),I=1,N)
21  FORMAT(5(I4,1PE12.4))
! correspoing best values of errors
    WRITE(LP,*)
    WRITE(LP,22)(I,W(NWF+I),I=1,M)
22     FORMAT(5(I4,1PE12.4))
    WRITE(LP,23)FMIN
23  FORMAT (/,15X,'THE SUM OF SQUARES IS',1PE17.8)
!     RESTORE THE BEST VALUES OF X AND F
19  DO I=1,N
       WSAVE(I)=X(I)
       X(I)=W(NWX+I)
    enddo
!  135 CONTINUE
    DO I=1,M
       IF(N+I.LT.MSAVE) WSAVE(N+I)=F(I)
       F(I)=W(NWF+I)
    enddo
!  136 CONTINUE
!...if we have not been able to save X and F then mark that no CONTINUE
    IF(N+M.GT.MSAVE) THEN
       IF(IEXIT(4).EQ.2) WRITE(LP,*)'Too big problem for CONTINUE'
       IEXIT(4)=0
    ENDIF
    RETURN
11  NTEST=NT
!     PROVIDE ORDINARY PRINTING IF REQUESTED
10  IF(IABS(IPRINT)-1) 39,38,40
38  CONTINUE
    IF (MAXC.EQ.0) WRITE(LP,41)MAXC,FSQ
    IF (MAXC.EQ.1) WRITE(LP,1041)MAXC,FSQ
    IF (MAXC.EQ.2) WRITE(LP,1141)MAXC,FSQ
    IF (MAXC.EQ.3) WRITE(LP,1241)MAXC,FSQ
    IF (MAXC.GE.4) WRITE(LP,41)MAXC,FSQ
41  FORMAT (/5X,' AT THE',I5,' TH ITERATION WE HAVE ',&
         'THE SUM OF SQUARES ',1PE17.8)
1041 FORMAT (/5X,' AT THE',I5,' ST ITERATION WE HAVE ',&
          'THE SUM OF SQUARES ',1PE17.8)
1141 FORMAT (/5X,' AT THE',I5,' ND ITERATION WE HAVE ',&
          'THE SUM OF SQUARES ',1PE17.8)
1241 FORMAT (/5X,' AT THE',I5,' RD ITERATION WE HAVE ',&
          'THE SUM OF SQUARES ',1PE17.8)
42  IF(IEXIT(2).EQ.1)WRITE(LP,21)(I,X(I),I=1,N)
    IF (IPRINT) 39,39,142
142 IF(IEXIT(3).EQ.1) THEN
       WRITE(LP,22)(I,F(I),I=1,M)
!...mark functions that changed a lot
       DO I=1,M
          IF(ABS(F(I)-W(NWF+I)).GT.2.0D1) WRITE(LP,7001)I
       enddo
!7000      CONTINUE
7001   FORMAT('Large jump for experiment ',I5)
          
    ENDIF
    GO TO 39
40  IPC=IPC-1
    IF (IPC) 43,43,39
43  WRITE(LP,44)MAXC
44  FORMAT (/5X,'THE BEST ESTIMATE AFTER',I5,' CALLS OF CALFUN IS')
    IPC=IABS(IPRINT)
    IF (FSQ-FMIN) 42,45,45
45  IF (FMIN) 42,46,46
46  WRITE(LP,21)(I,W(NWX+I),I=1,N)
    WRITE(LP,23)FMIN
    IF (IPRINT) 39,39,143
143 WRITE(LP,22) (I,W(NWF+I),I=1,M)
39  GO TO (49,47,47,48),IS
!     STORE THE INITIAL VECTORS X AND F
48  IF (IC) 50,50,51
50  DO I=1,N
       W(NWX+I)=X(I)
    enddo
!   52 CONTINUE
    GO TO 54
!     CALCULATE THE INITIAL JACOBIAN APPROXIMATION
51  K=IC
    DO I=1,M
       W(K)=(F(I)-W(NWF+I))/DSTEP
       K=K+N
    enddo
!   55 CONTINUE
!     TEST WHETHER THE MOST RECENT X IS BEST
    IF (FMIN-FSQ) 56,56,57
56  X(IC)=W(NWX+IC)
    GO TO 58
57  W(NWX+IC)=X(IC)
54  DO  I=1,M
       W(NWF+I)=F(I)
    enddo
!   53 CONTINUE
    FMIN=FSQ
!     SET X FOR THE NEXT CALL OF CALFUN
58  IC=IC+1
    IF (IC-N) 59,59,60
59  X(IC)=W(NWX+IC)+DSTEP
    GO TO 3
!     SET THE DIRECTION MATRIX
60  K=NWD
    DO I=1,N
       DO J=1,N
          K=K+1
          W(K)=0.0D0
       enddo
!   62 CONTINUE
       W(K+I-N)=1.0D0
       W(NWC+I)=1.0D0+FLOAT(N-I)
    enddo
!   61 CONTINUE
!     SET THE MARQUARDT PARAMETER TO ITS LEAST VALUE
24  PAR=PARM
!     COPY THE JACOBIAN AND APPEND THE MARQUARDT MATRIX
25  PPAR=PAR*PAR
    NTPAR=0
63  KK=0
    K=NWI+NWI
!      DO 26 I=1,N
!      DO 141 J=1,M
    DO I=1,N
       DO J=1,M
          KK=KK+1
          W(KK+NWI)=W(KK)
       enddo
!  141 CONTINUE
       DO J=1,N
          K=K+1
          W(K)=0.0D0
       enddo
!   27 CONTINUE
       W(K+I-N)=PAR
    enddo
!   26 CONTINUE
!     CALCULATE THE GENERALIZED INVERSE OF J
    CALL MB11A (N,MPN,W(NWI+1),N,W(NWW+1),ERR)
!     NOTE THAT THE THIRD AND FIFTH ENTRIES OF THIS ARGUMENT LIST
!     STAND FOR ONE-DIMENSIONAL ARRAYS.
!     START THE ITERATION BY TESTING FMIN
64  IF (FMIN-ACC) 7,7,65
!     NEXT PREDICT THE DESCENT AND MARQUARDT MINIMA
65  DS=0.0D0
    DN=0.0D0
    SP=0.0D0
    DO I=1,N
       X(I)=0.0D0
       F(I)=0.0D0
       K=I
       DO J=1,M
          X(I)=X(I)-W(K)*W(NWF+J)
          F(I)=F(I)-W(NWI+K)*W(NWF+J)
          K=K+N
       enddo
!   67 CONTINUE
       DS=DS+X(I)*X(I)
       DN=DN+F(I)*F(I)
       SP=SP+X(I)*F(I)
    enddo
!   66 CONTINUE
!     PREDICT THE REDUCTION IN F(X) DUE TO THE MARQUARDT STEP
!     AND ALSO PREDICT THE LENGTH OF THE STEEPEST DESCENT STEP
    PRED=SP+SP
    DMULT=0.0D0
    K=0
    DO I=1,M
       AP=0.0D0
       AD=0.0D0
       DO J=1,N
          K=K+1
          AP=AP+W(K)*F(J)
          AD=AD+W(K)*X(J)
       enddo
!   69 CONTINUE
       PRED=PRED-AP*AP
       DMULT=DMULT+AD*AD
    enddo
!   68 CONTINUE
!     TEST FOR CONVERGENCE
    IF (DN-DM) 28,28,29
28  AP=SQRT(DN)
    IF (PRED+2.0D0*PPAR*AP*(DMAX-AP)-ACC) 7,7,70
29  IF (PRED+PPAR*(DM-DN)-ACC) 7,7,70
!     TEST WHETHER TO APPLY THE FULL MARQUARDT CORRECTION
70  DMULT=DS/DMULT
    DS=DS*DMULT*DMULT
71  IS=2
    IF (DN-DD) 72,72,73
!     TEST THAT THE MARQUARDT PARAMETER HAS ITS LEAST VALUE
72  IF (PAR-PARM) 30,30,24
30  DD=MAX(DN,DSS)
    DS=0.25D0*DN
    TINC=1.0D0
    IF (DN-DSS) 74,132,132
74  IS=3
    GO TO 103
!     TEST WHETHER TO INCREASE THE MARQUARDT PARAMETER
73  IF (DN-DPAR) 31,31,32
31  NTPAR=0
    GO TO 33
32  IF (NTPAR) 34,34,35
34  NTPAR=1
    PTM=DN
    GO TO 33
35  NTPAR=NTPAR+1
    PTM=MIN (PTM,DN)
    IF (NTPAR-NT) 33,36,36
!     SET THE LARGER VALUE OF THE MARQUARDT PARAMETER
36  PAR=PAR*(PTM/DM)**0.25D0
    IF (6.0D0*DD-DM) 137,25,25
137    AP=SQRT(PRED/DN)
    IF (AP-PAR) 25,25,138
138 PAR=MIN(DBLE(AP),PAR*(DM/(6.0D0*DD))**0.25D0)
    GO TO 25
!     TEST WHETHER TO USE THE STEEPEST DESCENT DIRECTION
33  IF (DS-DD) 75,76,76
!     TEST WHETHER THE INITIAL VALUE OF DD HAS BEEN SET
76  IF (DD) 77,77,78
77  DD=MIN(DM,DS)
    IF (DD-DSS) 79,78,78
79  DD=DSS
    GO TO 71
!     SET THE MULTIPLIER OF THE STEEPEST DESCENT DIRECTION
78  ANMULT=0.0D0
    DMULT=DMULT*SQRT(DD/DS)
    GO TO 80
!     INTERPOLATE BETWEEN THE STEEPEST DESCENT AND MARQUARDT DIRECTIONS
75  SP=SP*DMULT
    ANMULT=(DD-DS)/((SP-DS)+SQRT((SP-DD)**2+(DN-DD)*(DD-DS)))
    DMULT=DMULT*(1.D0-ANMULT)
!     CALCULATE THE CORRECTION TO X, AND ITS ANGLE WITH THE FIRST
!     DIRECTION
80  DN=0.0D0
    SP=0.0D0
    DO I=1,N
       F(I)=DMULT*X(I)+ANMULT*F(I)
       DN=DN+F(I)*F(I)
       SP=SP+F(I)*W(NWD+I)
    enddo
!   81 CONTINUE
    DS=0.25D0*DN
!     TEST WHETHER AN EXTRA STEP IS NEEDED FOR INDEPENDENCE
    IF (W(NWC+1)-DTEST) 132,132,82
82  IF (SP*SP-DS) 83,132,132
!     TAKE THE EXTRA STEP AND UPDATE THE DIRECTION MATRIX
83  DO I=1,N
       X(I)=W(NWX+I)+DSTEP*W(NWD+I)
       W(NWC+I)=W(NWC+I+1)+1.0D0
    enddo
!   84 CONTINUE
    W(NWD)=1.0D0
    DO I=1,N
       K=NWD+I
       SP=W(K)
       DO J=2,N
          W(K)=W(K+N)
          K=K+N
       enddo
!   86 CONTINUE
       W(K)=SP
    enddo
!   85 CONTINUE
    GO TO 4
!     EXPRESS THE NEW DIRECTION IN TERMS OF THOSE OF THE DIRECTION
!     MATRIX, AND UPDATE THE COUNTS IN W(NWC+1) ETC.
132 SP=0.0D0
    K=NWD
!      DO 87 I=1,N
    DO I=1,N
       X(I)=DW
       DW=0.0D0
       DO J=1,N
          K=K+1
          DW=DW+F(J)*W(K)
       enddo
!   88 CONTINUE
       GO TO (89,90),IS
90     W(NWC+I)=W(NWC+I)+1.0D0
       SP=SP+DW*DW
       IF (SP-DS) 87,87,91
91     IS=1
       KK=I
       X(1)=DW
       GO TO 92
89     X(I)=DW
92     W(NWC+I)=W(NWC+I+1)+1.0D0
    enddo
87  CONTINUE
    W(NWD)=1.0D0
!     REORDER THE DIRECTIONS SO THAT KK IS FIRST
    IF (KK-1) 93,93,94
94  KS=NWC+KK*N
    DO I=1,N
       K=KS+I
       SP=W(K)
       DO J=2,KK
          W(K)=W(K-N)
          K=K-N
       enddo
!   96 CONTINUE
       W(K)=SP
    enddo
!   95 CONTINUE
!     GENERATE THE NEW ORTHOGONAL DIRECTION MATRIX
93  DO I=1,N
       W(NWW+I)=0.0D0
    enddo
!   97 CONTINUE
    SP=X(1)*X(1)
    K=NWD
    DO I=2,N
       DS=SQRT(SP*(SP+X(I)*X(I)))
       DW=SP/DS
       DS=X(I)/DS
       SP=SP+X(I)*X(I)
       DO J=1,N
          K=K+1
          W(NWW+J)=W(NWW+J)+X(I-1)*W(K)
          W(K)=DW*W(K+N)-DS*W(NWW+J)
       enddo
!   99 CONTINUE
    enddo
!   98 CONTINUE
    SP=1.0D0/SQRT(DN)
    DO I=1,N
       K=K+1
       W(K)=SP*F(I)
    enddo
!  100 CONTINUE
!     PREDICT THE NEW RIGHT HAND SIDES
    FNP=0.0D0
    K=0
    DO I=1,M
       W(NWW+I)=W(NWF+I)
       DO J=1,N
          K=K+1
          W(NWW+I)=W(NWW+I)+W(K)*F(J)
       enddo
!  102 CONTINUE
       FNP=FNP+W(NWW+I)**2
    enddo
!  101 CONTINUE
!     CALCULATE THE NEXT VECTOR X, AND THEN CALL CALFUN
103 DO I=1,N
       X(I)=W(NWX+I)+F(I)
    enddo
!  104 CONTINUE
    GO TO 4
!     UPDATE THE STEP SIZE
49  DMULT=0.90D0*FMIN+0.1D0*FNP-FSQ
    IF (DMULT) 105,108,108
105 DD=MAX(DBLE(DSS),0.25D0*DD)
    TINC=1.0D0
    IF (FSQ-FMIN) 106,107,107
!     TRY THE TEST TO DECIDE WHETHER TO INCREASE THE STEP LENGTH
108 SP=0.0D0
    SS=0.0D0
    DO I=1,M
       SP=SP+ABS(F(I)*(F(I)-W(NWW+I)))
       SS=SS+(F(I)-W(NWW+I))**2
    enddo
!  109 CONTINUE
    PJ=1.0D0+DMULT/(SP+SQRT(SP*SP+DMULT*SS))
    SP=MIN(4.0D0,DBLE(TINC),DBLE(PJ))
    TINC=PJ/SP
    DD=MIN(DM,SP*DD)
    GO TO 106
!     IF F(X) IMPROVES STORE THE NEW VALUE OF X
47  IF (FSQ-FMIN) 106,110,110
106 FMIN=FSQ
    DO I=1,N
       SP=X(I)
       X(I)=W(NWX+I)
       W(NWX+I)=SP
    enddo
!  111 CONTINUE
    DO I=1,M
       SP=F(I)
       F(I)=W(NWF+I)
       W(NWF+I)=SP
    enddo
!  112 CONTINUE
110 GO TO (107,107,113),IS
113 IS=2
    IF (FMIN-ACC) 7,7,83
!     CALCULATE THE CHANGES IN X AND IN F
107 DS=0.0D0
    DO I=1,N
       X(I)=X(I)-W(NWX+I)
       DS=DS+X(I)*X(I)
    enddo
!  114 CONTINUE
    DO I=1,M
       F(I)=F(I)-W(NWF+I)
    enddo
!  115 CONTINUE
!     CALCULATE THE GENERALIZED INVERSE TIMES THE CHANGE IN X
    K=NWI
    SS=0.0D0
    DO I=1,MPN
       SP=0.0D0
       DO J=1,N
          K=K+1
          SP=SP+W(K)*X(J)
       enddo
!  117 CONTINUE
       W(NWV+I)=SP
       SS=SS+SP*SP
    enddo
!  116 CONTINUE
!     CALCULATE J TIMES THE CHANGE IN F
!     ALSO APPLY PROJECTION TO THE GENERALIZED INVERSE
!      DO 118 I=1,N
    DO I=1,N
       ST=0.0D0
       K=NWI+I
       DO J=1,MPN
          ST=ST+W(K)*W(J+NWV)
          K=K+N
       enddo
!  119 CONTINUE
       ST=ST/SS
       K=NWI+I
       DO J=1,MPN
          W(K)=W(K)-ST*W(J+NWV)
          K=K+N
       enddo
!  120 CONTINUE
       ST=PPAR*X(I)
       K=I
       DO J=1,M
          ST=ST+W(K)*F(J)
          K=K+N
       enddo
!  121 CONTINUE
       W(NWW+I)=ST
    enddo
!  118 CONTINUE
!     REVISE J AND CALCULATE ROW VECTOR FOR CORRECTION TO INVERSE
    IC=0
    K=0
    KK=NWI
    SP=0.0D0
    SPP=0.0D0
    DO I=1,M
       SS=F(I)
       ST=F(I)
       DO J=1,N
          IC=IC+1
          KK=KK+1
          SS=SS-W(IC)*X(J)
          ST=ST-W(KK)*W(NWW+J)
       enddo
!  123 CONTINUE
       SS=SS/DS
       W(NWV+I)=ST
       SP=SP+F(I)*ST
       SPP=SPP+ST*ST
       DO J=1,N
          K=K+1
          W(K)=W(K)+SS*X(J)
       enddo
!  124 CONTINUE
    enddo
!  122 CONTINUE
    DO I=1,N
       ST=PAR*X(I)
       DO J=1,N
          KK=KK+1
          ST=ST-W(KK)*W(NWW+J)
       enddo
!  126 CONTINUE
       W(NWT+I)=ST
       SP=SP+PAR*X(I)*ST
       SPP=SPP+ST*ST
    enddo
!  125 CONTINUE
!     TEST THAT THE SCALAR PRODUCT IS SUFFICIENTLY ACCURATE
    IF (0.01D0*SPP-ABS(SP-SPP)) 63,63,127
!     CALCULATE THE NEW GENERALIZED INVERSE
!  127 DO 128 I=1,N
127 DO I=1,N
       K=NWI+I
       ST=X(I)
       DO J=1,M
          ST=ST-W(K)*F(J)
          K=K+N
       enddo
!  129 CONTINUE
       SS=0.0D0
       DO J=1,N
          SS=SS+W(K)*X(J)
          K=K+N
       enddo
!  130 CONTINUE
       ST=(ST-PAR*SS)/SP
       K=NWI+I
       DO J=1,MPN
          W(K)=W(K)+ST*W(NWV+J)
          K=K+N
       enddo
!  131 CONTINUE
    enddo
!  128 CONTINUE
    GO TO 64
  END SUBROUTINE VA05AD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!......INGA [NDRINGAR I MB11A
  SUBROUTINE MB11A (M,N,A,IA,W,ERR)
!...... 18/05/70 LAST LIBRARY UPDATE
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION A(IA,*),W(*)
!     PARTITION THE WORKING SPACE ARRAY W
!     THE FIRST PARTITION HOLDS THE FIRST COMPONENTS OF THE VECTORS OF
!     THE ELEMENTARY TRANSFORMATIONS
    NRW=M
!     THE SECOND PARTITION RECORDS ROW INTERCHANGES
    NCW=M+M
!     THE THIRD PARTITION RECORDS COLUMN INTERCHANGES
!     SET THE INITIAL RECORDS OF ROW AND COLUMN INTERCHANGES
    DO I=1,M
       W(NRW+I)=0.5D0+FLOAT(I)
    enddo
!    1 CONTINUE
    DO I=1,N
       W(NCW+I)=0.5D0+FLOAT(I)
    enddo
!    2 CONTINUE
!     'KK' COUNTS THE SEPARATE ELEMENTARY TRANSFORMATIONS
    KK=1
!     FIND LARGEST ROW AND MAKE ROW INTERCHANGES
3   RMAX=0.D0
!      DO 4 I=KK,M
    DO I=KK,M
       SUM=0.D0
       DO J=KK,N
          SUM=SUM+A(I,J)**2
       enddo
!    5 CONTINUE
       IF (RMAX-SUM) 6,4,4
6      RMAX=SUM
       IR=I
4      CONTINUE
    enddo
    IF (IR-KK) 7,7,8
8   SUM=W(NRW+KK)
    W(NRW+KK)=W(NRW+IR)
    W(NRW+IR)=SUM
    DO J=1,N
       SUM=A(KK,J)
       A(KK,J)=A(IR,J)
       A(IR,J)=SUM
    enddo
!    9 CONTINUE
!     FIND LARGEST ELEMENT OF PIVOTAL ROW, AND MAKE COLUMN INTERCHANGES
7   RMAX=0.D0
    SUM=0.D0
!      DO 10 J=KK,N
    DO J=KK,N
       SUM=SUM+A(KK,J)**2
       IF (RMAX-ABS(A(KK,J))) 11,10,10
11     RMAX=ABS(A(KK,J))
       IR=J
10     CONTINUE
    enddo
    IF (IR-KK) 12,12,13
13  RMAX=W(NCW+KK)
    W(NCW+KK)=W(NCW+IR)
    W(NCW+IR)=RMAX
    DO  I=1,M
       RMAX=A(I,KK)
       A(I,KK)=A(I,IR)
       A(I,IR)=RMAX
    enddo
!   14 CONTINUE
!     REPLACE THE PIVOTAL ROW BY THE VECTOR OF THE TRANSFORMATION
12  SIGMA=SQRT(SUM)
    BSQ=SQRT(SUM+SIGMA*ABS(A(KK,KK)))
!..ERROR SET IF BSQ=0 ???
    IF(BSQ.EQ.0.0D0)GOTO 9999
    W(KK)=SIGN(SIGMA+ABS(A(KK,KK)),A(KK,KK))/BSQ
    A(KK,KK)=-SIGN(SIGMA,A(KK,KK))
    KP=KK+1
    IF (KP-N) 15,15,16
15  DO J=KP,N
       A(KK,J)=A(KK,J)/BSQ
    enddo
!   17 CONTINUE
!     APPLY THE TRANSFORMATION TO THE REMAINING ROWS OF A
    IF (KP-M) 18,18,16
!18  DO 19 I=KP,M
18  DO I=KP,M
       SUM=W(KK)*A(I,KK)
       DO J=KP,N
          SUM=SUM+A(KK,J)*A(I,J)
       enddo
!   20 CONTINUE
       A(I,KK)=A(I,KK)-SUM*W(KK)
       DO J=KP,N
          A(I,J)=A(I,J)-SUM*A(KK,J)
       enddo
!   21 CONTINUE
    enddo
!   19 CONTINUE
    KK=KP
    GO TO 3
!     AT THIS STAGE THE REDUCTION OF A IS COMPLETE
!     NOW WE BUILD UP THE GENERALIZED INVERSE
!     APPLY THE FIRST ELEMENTARY TRANSFORMATION
16  KK=M
    KP=M+1
    SUM=W(M)/A(M,M)
    IF (N-M) 33,33,34
34  DO J=KP,N
       A(M,J)=-SUM*A(M,J)
    enddo
!   35 CONTINUE
33  A(M,M)=1.D0/A(M,M)-SUM*W(M)
!     NOW APPLY THE OTHER (M-1) TRANSFORMATIONS
36  KP=KK
    KK=KP-1
    IF (KK) 37,37,38
!     FIRST TRANSFORM THE LAST (M-KK) ROWS
!   38 DO 39 I=KP,M
38  DO I=KP,M
       SUM=0.D0
       DO J=KP,N
          SUM=SUM+A(KK,J)*A(I,J)
       enddo
!   40 CONTINUE
       DO J=KP,N
          A(I,J)=A(I,J)-SUM*A(KK,J)
       enddo
!   41 CONTINUE
       W(I)=-SUM*W(KK)
    enddo
!   39 CONTINUE
!     THEN CALCULATE THE NEW ROW IN POSITION KK
    DO J=KP,N
       SUM=-W(KK)*A(KK,J)
       DO I=KP,M
          SUM=SUM-A(I,KK)*A(I,J)
       enddo
!   43 CONTINUE
       A(KK,J)=SUM/A(KK,KK)
    enddo
!   42 CONTINUE
!     AND REVISE THE COLUMN IN POSITION KK
    SUM=1.D0-W(KK)**2
    DO I=KP,M
       SUM=SUM-A(I,KK)*W(I)
       A(I,KK)=W(I)
    enddo
!   44 CONTINUE
    A(KK,KK)=SUM/A(KK,KK)
    GO TO 36
!     RESTORE THE ROW INTERCHANGES
!   37 DO 45 I=1,M
37  DO I=1,M
!   46 IR=IFIX(W(NRW+I))
46     IR=INT(W(NRW+I))
       IF (I-IR) 47,45,45
47     SUM=W(NRW+I)
       W(NRW+I)=W(NRW+IR)
       W(NRW+IR)=SUM
       DO J=1,N
          SUM=A(I,J)
          A(I,J)=A(IR,J)
          A(IR,J)=SUM
       enddo
!   48 CONTINUE
       GO TO 46
45     continue
    enddo
!   45 CONTINUE
!     RESTORE THE COLUMN INTERCHANGES
!      DO 49 J=1,N
    DO J=1,N
!   50 IR=IFIX(W(NCW+J))
50     IR=INT(W(NCW+J))
       IF (J-IR) 51,49,49
51     SUM=W(NCW+J)
       W(NCW+J)=W(NCW+IR)
       W(NCW+IR)=SUM
       DO I=1,M
          SUM=A(I,J)
          A(I,J)=A(I,IR)
          A(I,IR)=SUM
       enddo
!   52 CONTINUE
       GO TO 50
49     continue
    enddo
!   49 CONTINUE
    ERR=1.0D0
9999 RETURN
  END SUBROUTINE MB11A

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE EA06CD(A,VALUE,VECTOR,M,IA,IV,W)
!######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
    DOUBLE PRECISION  A,PP,VALUE,VECTOR,W
    DIMENSION A(IA,M),VALUE(M),VECTOR(IV,M),W(*)
    M1=M+1
    M2=M1+M
    W(1)=A(1,1)
    IF(M-2)60,10,15
10  W(2)=A(2,2)
    W(4)=A(2,1)
    GO TO 60
15  CALL MC04BD(A,W,W(M1),M,IA,W(M2))
60  CALL EA08CD(W,W(M1),VALUE,VECTOR,M,IV,W(M2))
    IF(M.LE.2)RETURN
!      DO 56 L=1,M
!      DO 56 II=3,M
    DO L=1,M
       DO II=3,M
          I=M-II+1
          M3=M1+I
          IF(W(M3))57,56,57
57        PP=0.0D0
          I1=I+1
          DO K=I1,M
             PP=PP+A(I,K)*VECTOR(K,L)
          enddo
!58    PP=PP+A(I,K)*VECTOR(K,L)
          PP=PP/(A(I,I+1)*W(M3))
          DO K=I1,M
             VECTOR(K,L)=VECTOR(K,L)+PP*A(I,K)
          enddo
!59    VECTOR(K,L)=VECTOR(K,L)+PP*A(I,K)
56        continue
       enddo
    enddo
!56    CONTINUE
    RETURN
  END SUBROUTINE EA06CD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE EA08CD(A,B,VALUE,VEC,M,IV,W)
!######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
!  STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
    DOUBLE PRECISION A,A11,A12,A13,A21,A22,A23,A33,A34,B,BB,CC,&
         CO,EPS,ROOT,S,SI,SML,VALUE,VEC,V1,V2,W,XAX,XX
    DIMENSION A(M),B(M),VALUE(M),VEC(1),W(1)
    DATA EPS/2.3D-16/,A34/0.0D0/
!     THIS USES QR ITERATION TO FIND THE EIGENVALUES AND EIGENVECTORS
!  OF THE SYMMETRIC TRIDIAGONAL MATRIX WHOSE DIAGONAL ELEMENTS ARE
!  A(I),I=1,M AND OFF-DIAGONAL ELEMENTS ARE B(I),I=2,M.  THE ARRAY
!  W IS USED FOR WORKSPACE AND MUST HAVE DIMENSION AT LEAST 2*M.
!  WE TREAT VEC AS IF IT HAD DIMENSIONS (IV,M).
    SML=EPS*FLOAT(M)
    CALL EA09CD(A,B,W(M+1),M,W)
!     SET VEC TO THE IDENTITY MATRIX.
    DO I=1,M
       VALUE(I)=A(I)
       W(I)=B(I)
       K=(I-1)*IV+1
       L=K+M-1
       DO J=K,L
!    3 VEC(J)=0.0D0
          VEC(J)=0.0D0
       enddo
       KI=K+I
       VEC(KI-1)=1.D0
!5     VEC(KI-1)=1.D0
    enddo
    K=0
    IF(M.EQ.1)RETURN
!      DO 200 N3=2,M
    DO N3=2,M
       N2=M+2-N3
!     EACH QR ITERATION IS PERFORMED OF ROWS AND COLUMNS N1 TO N2
       MN2=M+N2
       ROOT=W(MN2)
!      DO 190 ITER=1,20
       DO ITER=1,20
          BB=(VALUE(N2)-VALUE(N2-1))*0.5D0
          CC=W(N2)*W(N2)
          A22=VALUE(N2)
          IF(CC.NE.0.0D0)A22=A22+CC/(BB+DSIGN(1.0D0,BB)*DSQRT(BB*BB+CC))
          DO I=1,N2
             MI=M+I
             IF(DABS(ROOT-A22).LE.DABS(W(MI)-A22))GO TO 125
             ROOT=W(MI)
             MN=M+N2
             W(MI)=W(MN)
             W(MN)=ROOT
125          continue
          enddo
!125   CONTINUE
          DO II=2,N2
             N1=2+N2-II
             IF(DABS(W(N1)).LE.(DABS(VALUE(N1-1))+DABS(VALUE(N1)))*SML)GOTO140
          enddo
!130   CONTINUE
          N1=1
140       IF(N2.EQ.N1) GO TO 200
          N2M1=N2-1
          IF(ITER.GE.3)ROOT=A22
          K=K+1
          A22=VALUE(N1)
          A12=A22-ROOT
          A23=W(N1+1)
          A13=A23
          DO180 I=N1,N2M1
          A33=VALUE(I+1)
          IF(I.NE.N2M1)A34=W(I+2)
          S=DSIGN(DSQRT(A12*A12+A13*A13),A12)
          SI=A13/S
          CO=A12/S
          JK=I*IV+1
          J1=JK-IV
          J2=J1+MIN0(M,I+K)-1
!      DO 160 JI=J1,J2
          DO JI=J1,J2
             V1=VEC(JI)
             V2=VEC(JK)
             VEC(JI)=V1*CO+V2*SI
             VEC(JK)=V2*CO-V1*SI
160          JK=JK+1
!160   JK=JK+1
          enddo
          IF(I.NE.N1)  W(I)=S
          A11=CO*A22+SI*A23
          A12=CO*A23+SI*A33
          A13=SI*A34
          A21=CO*A23-SI*A22
          A22=CO*A33-SI*A23
          A23=CO*A34
          VALUE(I)=A11*CO+A12*SI
          A12=-A11*SI+A12*CO
          W(I+1)=A12
180       A22=A22*CO-A21*SI
!  190 VALUE(N2)=A22
190       VALUE(N2)=A22
       enddo
!      WRITE(6,195)
!  195 FORMAT('CYCLE DETECTED IN SUBROUTINE EA08,',&
!           ' TURNING OFF STABILITY CHECK')
!     STOP
!       CALL ST2ERR(1209,'EA08','LOOPING IN STABILITY CHECK')
       RETURN
200    CONTINUE
    enddo
!     RAYLEIGH QUOTIENT!
!    DO 220 J=1,M
    DO J=1,M
       K=(J-1)*IV
       XX=VEC(K+1)**2
       XAX=XX*A(1)
!      DO 210 I=2,M
       DO I=2,M
          KI=K+I
          XX=XX+VEC(KI)**2
210       XAX=XAX+VEC(KI)*(2.0D0*B(I)*VEC(KI-1)+A(I)*VEC(KI))
!  210 XAX=XAX+VEC(KI)*(2.0D0*B(I)*VEC(KI-1)+A(I)*VEC(KI))
       enddo
220    VALUE(J)=XAX/XX
    enddo
    RETURN
  END SUBROUTINE EA08CD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE EA09CD(A,B,VALUE,M,OFF)
!######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
!  STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
    DOUBLE PRECISION A,A11,A12,A13,A21,A22,A23,A33,A34,B,BB,&
         CC,CO,EPS,OFF,ROOT,S,SBB,SI,SML,VALUE
    DIMENSION A(M),B(M),VALUE(M),OFF(M)
    DATA A34/0.0D0/,EPS/2.3D-16/
    SML=EPS*FLOAT(M)
    VALUE(1)=A(1)
    IF(M.EQ.1)RETURN
    DO I=2,M
       VALUE(I)=A(I)
10     OFF(I)=B(I)
    enddo
!10    OFF(I)=B(I)
!     EACH QR ITERATION IS PERFORMED OF ROWS AND COLUMNS N1 TO N2
    MAXIT=10*M
!      DO 90 ITER=1,MAXIT
!      DO 45 N3=2,M
    DO ITER=1,MAXIT
       DO N3=2,M
          N2=M+2-N3
          DO II=2,N2
             N1=2+N2-II
             IF(DABS(OFF(N1)).LE.(DABS(VALUE(N1-1))+DABS(VALUE(N1)))*SML)GOTO40
          enddo
!30    CONTINUE
          N1=1
40        IF(N2.NE.N1) GO TO 50
!   45 CONTINUE
       enddo
       RETURN
!     ROOT  IS THE EIGENVALUE OF THE BOTTOM 2*2 MATRIX THAT IS NEAREST
!     TO THE LAST MATRIX ELEMENT AND IS USED TO ACCELERATE THE
!     CONVERGENCE
50     BB=(VALUE(N2)-VALUE(N2-1))*0.5D0
       CC=OFF(N2)*OFF(N2)
       SBB=1.0D0
       IF(BB.LT.0.0D0)SBB=-1.0D0
       ROOT=VALUE(N2)+CC/(BB+SBB*DSQRT(BB*BB+CC))
       N2M1=N2-1
75     A22=VALUE(N1)
       A12=A22-ROOT
       A23=OFF(N1+1)
       A13=A23
!      DO 80 I=N1,N2M1
       DO I=N1,N2M1
          A33=VALUE(I+1)
          IF(I.NE.N2M1)A34=OFF(I+2)
          S=DSQRT(A12*A12+A13*A13)
          SI=A13/S
          CO=A12/S
          IF(I.NE.N1)OFF(I)=S
          A11=CO*A22+SI*A23
          A12=CO*A23+SI*A33
          A13=SI*A34
          A21=CO*A23-SI*A22
          A22=CO*A33-SI*A23
          A23=CO*A34
          VALUE(I)=A11*CO+A12*SI
          A12=-A11*SI+A12*CO
          OFF(I+1)=A12
!80    A22=A22*CO-A21*SI
80        A22=A22*CO-A21*SI
       enddo
90     VALUE(N2)=A22
!   90 VALUE(N2)=A22
    enddo
!       WRITE(6,100)
!  100 FORMAT('LOOPING DETECTED IN EA09, TURNING OFF STABILITY CHECK')
!     STOP
    RETURN
  END SUBROUTINE EA09CD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE MC04BD(A,ALPHA,BETA,M,IA,Q)
!######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
!     STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
    DOUBLE PRECISION A,ALPHA,BETA,BIGK,H,Q,QJ,PP,PP1
    DIMENSION A(IA,1),ALPHA(1),BETA(1),Q(1)
    ALPHA(1)=A(1,1)
    DO J=2,M
       J1=J-1
       DO I=1,J1
          A(I,J)=A(J,I)
       enddo
!   22 CONTINUE
       ALPHA(J)=A(J,J)
    enddo
!   21 CONTINUE
    M1=M-1
    M2=M-2
!      DO 1 I=1,M2
    DO I=1,M2
       PP=0.0D0
       I1=I+1
       DO J=I1,M
          PP=PP+A(I,J)**2
       enddo
!    2 CONTINUE
       PP1=DSQRT(PP)
       IF(A(I,I+1))3,5,5
5      BETA(I+1)=-PP1
       GO TO 6
3      BETA(I+1)=PP1

6      IF(PP)1,1,17

17     H=PP-BETA(I+1)*A(I,I+1)
       A(I,I+1)=A(I,I+1)-BETA(I+1)
!      DO 7 KI=I1,M
       DO KI=I1,M
          QJ=0.0D0
          DO KJ=I1,KI
             QJ=QJ+A(KJ,KI)*A(I,KJ)
          enddo
!    8 CONTINUE
          IF(KI-M)19,20,20
19        I2=KI+1
          DO KJ=I2,M
             QJ=QJ+A(KI,KJ)*A(I,KJ)
          enddo
!   18 CONTINUE
20        Q(KI)=QJ/H
       enddo
!    7 CONTINUE
       BIGK=0.0D0
       DO KJ=I1,M
          BIGK=BIGK+A(I,KJ)*Q(KJ)
       enddo
!    9 CONTINUE
       BIGK=BIGK/(2.0*H)
       DO KJ=I1,M
          Q(KJ)=Q(KJ)-BIGK*A(I,KJ)
       enddo
!   10 CONTINUE
       DO KI=I1,M
          DO KJ=KI,M
             A(KI,KJ)=A(KI,KJ)-Q(KI)*A(I,KJ)-Q(KJ)*A(I,KI)
          enddo
       enddo
!   12 CONTINUE
!   11 CONTINUE
1      continue
    enddo
!1   CONTINUE
    DO I=2,M
       H=ALPHA(I)
       ALPHA(I)=A(I,I)
       A(I,I)=H
    enddo
!   23 CONTINUE
    BETA(M)=A(M-1,M)
    RETURN
  END SUBROUTINE MC04BD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

end module LIBOCEQPLUS

