!
! general utilities in Fortran 95 a la METLIB upgraded 2015-2019 to
! eliminate most specific F77 features
! 1. All ENTRY removed.  GPARxyz and MACRO routines seems OK
! 2. Problems with getkey developed by John S. Urban has been fixed,
!    It has been renamed getkex in the iso-C interface.
! 3. IMPLICT NONE introduced in the whole module.
! 4. A revised online help system using HTML \hypertarget in user guide
!
! To be done:
! - move constants ZERO, ONE here (from gtp3)
! - use same error code system as in gtp
! - revise the online help system
!
MODULE METLIB
!
! Copyright 1980-2021, Bo Sundman bo.sundman@gmail.com 
! 
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!--------------------------------------------------------------------------
#ifdef lixed
! LINUX: For character by character input allowing emacs editing
  use M_getkey
#endif
!
#ifdef tinyfd
! use tinyfiledialogs to browse for files
  use ftinyopen
#endif
!
!--------------------------------------------------------------------------
!
! Error codes (buperr)
!  1001 Too small stack for sorting integers
!  1002 Too small workspace
!  1003 Pointer outside workspace
!  1004 Free areas not in increasing order
!  1005 Too small or too big free area
!  1006 No space available
!  1007 Free workspace destroyed
!  1008 Attempt to reserve one word or less
!  1009 Released area inside free area
!  1010 Too large character or real arrays in LOADC/STORC or LOADRN/STORRN
!  1030 NO SUCH TYPE OPTION
!  1031 Empty line, expected number
!  1032 PARAMETER VALUE MISSING
!  1033 Decimal point but no digits
!  1034 No digits
!  1035 Positive sign but no digits
!  1036 Negative sign but no digits
!  1037 No sign and no digits
!  1038 NO DIGITS AFTER EXPONENTIAL E
!  1039 Exponent larget then 99
!  1040 NO HELP FOR <COMMAND>
!  1041 NO SUCH QUESTION FOR <COMMAND>
!  1042 TOO LARGE INTEGER VALUE
!  1057 Too long input line
!  1060 illegal bit number
!  1070 Margin error in wrice
!  1083 Too deeply nested macros
!  1100 Not enough space to write number in text
!  1101 Name does not start with letter A-Z
!  1235 too many (
!  1236 too few )
!  1237 error with parenthesis
!  1332 Illegal option
!  1333 Missing option value delimiter
!  1334 Missing option value
!  1350 File system error
!  1360 Missing column number for substitution
!
!----------------------------------------------------------
!
! global variables and constants
!  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  implicit none
! Data structures for putfun below (no change 190802)
!\begin{verbatim}
! Data structures in METLIB
  TYPE putfun_node
! all nodes of function stored as part of a binary tree
! kod is operation kod (0 datanod), links is how many links to this node
     integer kod,links
! this is the sequential order the node is allocated (for debugging)
     integer debug
! each node has a left and a right link.  If the left node is empty the
! right is normally a data node
     TYPE(putfun_node), pointer :: left,right
! A data node can have a  numeric value and/or a link to another function
     double precision value
! this is an identification of external symbols
     integer dataid
  end TYPE putfun_node
!
! BEWARE entering putfuns cannot be made in parallel processing
! but one may evaluate them in different threads
!
! PUTFUNNVAR is associated with external symbols in the LOKV array
  integer, private :: putfunvar
  TYPE PUTFUN_STACK
     type(putfun_node), pointer ::savetop, savebin, saveuni
     type(putfun_stack), pointer :: previous
  end TYPE PUTFUN_STACK
  type(putfun_stack), pointer :: stacktop
! topnod is the current top node
! lastopnod is last binary opkod node
! datanod is last data node
  TYPE(putfun_node), private, pointer :: topnod,datanod,lastopnod
  integer pfnerr,debuginc
!
! end data structures for PUTFUN
!\end{verbatim}
!
!\begin{verbatim}
! data structures for history and help
!  
  integer, parameter :: histlines=100
!
  TYPE CHISTORY
! to save the last 20 lines of commands
     character*80 hline(histlines)
     integer :: hpos=0
  END TYPE CHISTORY
  type(chistory) :: myhistory
!  
    integer, parameter :: maxhelplevel=15
! A help structure used in new on-line help system
! this was designed for both LaTeX and HTML help, now only HTML
    TYPE help_str
       integer :: okinit=0
       character*128 filename
       character*8 type
       integer level
       character*32, dimension(maxhelplevel) :: cpath
    END TYPE help_str
! this record is used to file the appropriate help text
    type(help_str), save :: helprec
! this is useful to add %\section and %\subsection in helpfile
    logical :: helptrace=.FALSE.
!
! using browser and html files for on-line help
  type onlinehelp
! if htmlhelp is TRUE then browser is the path/name of browser
! htmlfile is full path/name of html file
! target is used to find the relevant text the html file
! values of browser and htmlfile set by the main program (and htmlhelp=.TRUE.)
! The value of target is found searching the original LaTeX file!!
! In this file there are \hypertarget{target} which can be searched in the
! html file as <a id="target" />
! Searching the LaTeX file the help system will find a section
! matching the history of commands/questions the user has given
! and the target in the first \hypertarget {target} found within these lines
! will be used for the help displayed in the browser window
     logical :: htmlhelp=.FALSE.
     character*128 browser
     character*128 htmlfile
     character*128 latexfile
     character*64 target
  end type onlinehelp
  type(onlinehelp) :: ochelp
  save ochelp
! end data structures for history and help
!\end{verbatim}
!  
!
! default units for command output, input, error message, list
! and default language
  integer, parameter :: koud=6,kiud=5,keud=6,lutd=6,lerd=6,langd=1
! representation of the numerical value "none"
  double precision, parameter :: RNONE=-1.0D-36,FLTSML=1.0D-36,FLTPRS=1.0D-14
  integer, parameter :: NONE=-2147483647,MAXINT=2147483647,MININT=-2147483646
  character*4, parameter :: CNONE='NONE'
! initiate i/o and error code
  integer :: kou=koud,kiu=kiud,keu=keud,lut=lutd,ler=lerd,lang=langd
  integer :: buperr=0,iox(10)
! LSTCMD is the last command given. Saved by NCOMP, used by help routines
  character, private :: lstcmd*40
! LUN unsed for macros and SAVE files
  integer :: lun=50
! LOGFIL is nonzero if a log file is set
  integer, private :: logfil=0
! global values for history
  CHARACTER, private :: HIST(20)*80
  integer, private :: LHL=0,LHM=0,LHP=0
! terminal charcterististics, koltrm is number of columns, default 80
  integer :: KOLTRM=80
! ECHO on/off
  integer :: JECHO=0
! no idea what is in KFLAGS, has to do with VT200 terminals. Not needed??
!  integer KFLAGS(24)
!
! This is for environment variables used in MACROs
  character, private :: ENVIR(9)*60
!
  character*3, parameter :: MACEXT='OCM'
!
!----------- some added things 190802/BoS
! prevent use of popup windows for open/save file
! logical nopopup
  logical :: NOPENPOPUP=.FALSE.
! character for PATH to macro file in order to open files inside macro
    character macropath(5)*256
! the working directory
    character workingdir*256
! for macros
    integer IUMACLEVL,MACROUNIT(5)
!\begin{verbatim}
! >>>>>>>>>> SYSTEM DEPENDENT <<<<<<<<<<
! nbpw is number if bytes per INTEGER, nwpr number of words per (double) real 
! nbitpw number of bits per word
! USED when WPACK routines store data in integer workspace 
    integer, parameter :: nbpw=4,nwpr=2,nbitpw=32
! >>>>>>>>>> SYSTEM DEPENDENT <<<<<<<<<<
!\end{verbatim}
!
! some constants
  double precision, parameter, private :: ZERO=0.0D0,ONE=1.0D0,TEN=1.0D1
!
! -------------------------------------------------------------------
! GPARxyz routines parameter transfer of integer, real and logical values
!    
!  integer GPARIDEF,GPARITYP
!  double precision GPARRDEF
!  logical GPARWDEF,GPARENTES
!  character GPARCH2*1
!
! added private to avoid problem with modules using metlib
  integer, private :: GPARIDEF,GPARITYP
  double precision, private :: GPARRDEF
  logical, private :: GPARWDEF,GPARENTES
  character, private :: GPARCH2*1
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CONTAINS

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
!      SORTING SUBROUTINES FOR INTEGERS, REALS AND CHARACTERS
!      QUICKSORT ENL KNUTH ALGORTIM Q
!      THE ART OF COMPUTER PROGRAMMING, VOL 3, P 117
!  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine sortrd & Sorting reals
!\begin{verbatim}
  SUBROUTINE SORTRD(ARR,N,IX)
! ...SORTING REAL NUMBERS IN ASCENDING ORDER
! INPUT:
!      ARR   ARRAY TO BE SORTED
!      N     NUMBER OF ELEMENTS TO BE SORTED >1
!      IX    INTEGER ARRAY WITH DIMENSION N
! EXIT:
!      ARR   SORTED ARRAY
!      IX    ARRAY WHERE IX(I) IS THE PREVIOS INDEX OF ARR(I)
    implicit none
    real ARR(*)
    integer n,ix(*)
!\end{verbatim} %+
    integer, parameter :: MSTACK=20
    real part,val
    integer LOW(MSTACK),IGH(MSTACK)
    integer i,is,j,k,m,min,max
!      LOW AND IGH IS USED TO STORE THE LOWER AND HIGHER PARTION BOUNDARIES
    IF(N.LT.1) GOTO 900
!      IX IS ORIGINAL INDEX OF ELEMENT I IN ARR
    DO I=1,N
       IX(I)=I
    enddo
    IF(N.EQ.1) GOTO 900
!      M IS THE PARTION SIZE THAT IS SORTED WITH STRIGHT INSERTION
!      IS POINTS TO FREE STACK, MSTACK IS SIZE OF STACK
    M=1
    IS=1
!******STEP Q1, INITIATING
    IF(N.LE.M) GOTO 900
    MIN=1
    MAX=N
!******STEP Q2, NEW STAGE
!      MIN AND MAX ARE LOWER AND UPPER LIMITS FOR THE PARTION
100 PART=ARR(MIN)
    I=MIN
    J=MAX+1
!******STEP Q3, INCREASE I UNTIL I>J OR ARR(I)>PART
110 I=I+1
    IF(I.GE.J) GOTO 200
    IF(ARR(I).LE.PART) GOTO 110
!******STEP Q4, DECREASE J UNTIL J<I OR ARR(J)<PART
200 J=J-1
    IF(J.LT.I) GOTO 300
    IF(ARR(J).GT.PART) GOTO 200
!******STEP Q6, SWITCH ARR(I) AND ARR(J) AND CONTINUE FROM Q3
    VAL=ARR(I)
    ARR(I)=ARR(J)
    ARR(J)=VAL
    K=IX(I)
    IX(I)=IX(J)
    IX(J)=K
    GOTO 110
!******STEP Q5, I AD J HAVE PASSED EACHOTHER, SWITHCH PART=ARR(MIN) AND ARR(J)
300 ARR(MIN)=ARR(J)
    ARR(J)=PART
    K=IX(MIN)
    IX(MIN)=IX(J)
    IX(J)=K
!******STEP Q7, PUSH THE GREATEST PARTITION ON STACK
    IF(MAX-J.GT.J-MIN) GOTO 350
!      J-MIN GREATEST
    IF(J-MIN.LE.M) GOTO 400
!      PUSH ONLY IF MAX-J>M
    IF(MAX-J.LE.M) GOTO 360
!      BOTH PARTITIONS ARE GREATER THAN M, THE GREATEST IS PUSHED
    IF(IS.GT.MSTACK) GOTO 910
    LOW(IS)=MIN
    IGH(IS)=J-1
    IS=IS+1
!      CONTINUE TO PARTITION THE SMALLEST
310 MIN=J+1
    GOTO 100
!      MAX-J GREATEST
350 IF(MAX-J.LE.M) GOTO 400
!      PUSH ONLY IF J-MIN>M
    IF(J-MIN.LE.M) GOTO 310
!      BOTH PARTITIONS ARE GREATER THAN M, PUSH THE GREATEST
    IF(IS.GT.MSTACK) GOTO 910
    LOW(IS)=J+1
    IGH(IS)=MAX
    IS=IS+1
!      CONTINUE TO PARTITION THE SMALLEST
360 MAX=J-1
    GOTO 100
!******STEP Q8, POP FROM STACK
400 IS=IS-1
    IF(IS.LT.1) GOTO 500
    MIN=LOW(IS)
    MAX=IGH(IS)
    GOTO 100
!******STEP Q9, STRIGHT INSERTION, ONLY NECESSARY IF M>1
500 CONTINUE
900 RETURN
910 buperr=1051
    GOTO 900
  end SUBROUTINE SORTRD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine sortrdd & Sorting doubles
!\begin{verbatim}
  SUBROUTINE SORTRDD(ARR,N,IX)
! ...SORTING DOUBLE PRECISION NUMBERS IN DECENDING ORDER
! INPUT:
!      ARR   ARRAY TO BE SORTED
!      N     NUMBER OF ELEMENTS TO BE SORTED >1
!      IX    INTEGER ARRAY WITH DIMENSION N
! EXIT:
!      ARR   SORTED ARRAY
!      IX    ARRAY WHERE IX(I) IS THE PREVIOS INDEX OF ARR(I)
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    double precision ARR(*)
    integer n,ix(*)
!\end{verbatim} %+
    integer, parameter :: MSTACK=20
    double precision part,val
    integer LOW(MSTACK),IGH(MSTACK)
    integer i,is,j,k,m,min,max
!      LOW AND IGH IS USED TO STORE THE LOWER AND HIGHER PARTION BOUNDARIES
    IF(N.LT.1) GOTO 900
!      IX IS ORIGINAL INDEX OF ELEMENT I IN ARR
    DO I=1,N
       IX(I)=I
    enddo
    IF(N.EQ.1) GOTO 900
!      M IS THE PARTION SIZE THAT IS SORTED WITH STRIGHT INSERTION
!      IS POINTS TO FREE STACK, MSTACK IS SIZE OF STACK
    M=1
    IS=1
!******STEP Q1, INITIATING
    IF(N.LE.M) GOTO 900
    MIN=1
    MAX=N
!******STEP Q2, NEW STAGE
!      MIN AND MAX ARE LOWER AND UPPER LIMITS FOR THE PARTION
100 PART=ARR(MIN)
    I=MIN
    J=MAX+1
!******STEP Q3, INCREASE I UNTIL I>J OR ARR(I)>PART
110 I=I+1
    IF(I.GE.J) GOTO 200
!    IF(ARR(I).LE.PART) GOTO 110
    IF(ARR(I).GE.PART) GOTO 110
!******STEP Q4, DECREASE J UNTIL J<I OR ARR(J)<PART
200 J=J-1
    IF(J.LT.I) GOTO 300
!    IF(ARR(J).GT.PART) GOTO 200
    IF(ARR(J).LT.PART) GOTO 200
!******STEP Q6, SWITCH ARR(I) AND ARR(J) AND CONTINUE FROM Q3
    VAL=ARR(I)
    ARR(I)=ARR(J)
    ARR(J)=VAL
    K=IX(I)
    IX(I)=IX(J)
    IX(J)=K
    GOTO 110
!******STEP Q5, I AD J HAVE PASSED EACHOTHER, SWITHCH PART=ARR(MIN) AND ARR(J)
300 ARR(MIN)=ARR(J)
    ARR(J)=PART
    K=IX(MIN)
    IX(MIN)=IX(J)
    IX(J)=K
!******STEP Q7, PUSH THE GREATEST PARTITION ON STACK
    IF(MAX-J.GT.J-MIN) GOTO 350
!      J-MIN GREATEST
    IF(J-MIN.LE.M) GOTO 400
!      PUSH ONLY IF MAX-J>M
    IF(MAX-J.LE.M) GOTO 360
!      BOTH PARTITIONS ARE GREATER THAN M, THE GREATEST IS PUSHED
    IF(IS.GT.MSTACK) GOTO 910
    LOW(IS)=MIN
    IGH(IS)=J-1
    IS=IS+1
!      CONTINUE TO PARTITION THE SMALLEST
310 MIN=J+1
    GOTO 100
!      MAX-J GREATEST
350 IF(MAX-J.LE.M) GOTO 400
!      PUSH ONLY IF J-MIN>M
    IF(J-MIN.LE.M) GOTO 310
!      BOTH PARTITIONS ARE GREATER THAN M, PUSH THE GREATEST
    IF(IS.GT.MSTACK) GOTO 910
    LOW(IS)=J+1
    IGH(IS)=MAX
    IS=IS+1
!      CONTINUE TO PARTITION THE SMALLEST
360 MAX=J-1
    GOTO 100
!******STEP Q8, POP FROM STACK
400 IS=IS-1
    IF(IS.LT.1) GOTO 500
    MIN=LOW(IS)
    MAX=IGH(IS)
    GOTO 100
!******STEP Q9, STRIGHT INSERTION, ONLY NECESSARY IF M>1
500 CONTINUE
900 RETURN
910 buperr=1051
    GOTO 900
  end SUBROUTINE SORTRDD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine sortin & Sorting integers
!\begin{verbatim}
  SUBROUTINE SORTIN(IARR,N,IX)
! ...SORTING INTEGERS IN ASCENDING ORDER
! INPUT:
!      IARR   ARRAY TO BE SORTED
!      N     NUMBER OF ELEMENTS TO BE SORTED >1
!      IX    INTEGER ARRAY WITH DIMENSION N
! EXIT:
!      IARR   SORTED ARRAY
!      IX    ARRAY WHERE IX(I) IS THE PREVIOS INDEX OF IARR(I)
    implicit none
    integer IARR(*),n,ix(*)
!\end{verbatim} %+
    integer, parameter :: MSTACK=20
    integer ipart,ival
    integer LOW(MSTACK),IGH(MSTACK)
    integer i,is,j,k,m,min,max
!    PARAMETER (MSTACK=20)
!    DIMENSION IARR(*),IX(*),LOW(MSTACK),IGH(MSTACK)
!      LOW AND IGH IS USED TO STORE THE LOWER AND HIGHER PARTISION BOUNDARIES
    IF(N.LT.1) GOTO 900
!      IX IS ORIGINAL INDEX OF ELEMENT I IN IARR
    DO I=1,N
       IX(I)=I
    endDO
    IF(N.EQ.1) GOTO 900
!      M IS THE PARTITION SIZE THAT IS SORTED WITH STRIGHT INSERTION
!      IS POINTS TO FREE STACK, MSTACK IS SIZE OF STACK
    M=1
    IS=1
!******STEP Q1, INITIATING
    IF(N.LE.M) GOTO 900
    MIN=1
    MAX=N
!******STEP Q2, NEW STAGE
!      MIN AND MAX ARE LOWER AND UPPER LIMITS FOR THE PARTISION
100 IPART=IARR(MIN)
    I=MIN
    J=MAX+1
!******STEP Q3, INCREASE I UNTIL I>J OR IARR(I)>IPART
110 I=I+1
    IF(I.GE.J) GOTO 200
    IF(IARR(I).LE.IPART) GOTO 110
!******STEP Q4, DECREASE J UNTIL J<I OR IARR(J)<IPART
200 J=J-1
    IF(J.LT.I) GOTO 300
    IF(IARR(J).GT.IPART) GOTO 200
!******STEP Q6, SWITCH IARR(I) AND IARR(J) AND CONTINUE FROM Q3
    IVAL=IARR(I)
    IARR(I)=IARR(J)
    IARR(J)=IVAL
    K=IX(I)
    IX(I)=IX(J)
    IX(J)=K
    GOTO 110
!*****STEP Q5, I AD J HAVE PASSED EACHOTHER, SWITCH IPART=IARR(MIN) AND IARR(J)
300 IARR(MIN)=IARR(J)
    IARR(J)=IPART
    K=IX(MIN)
    IX(MIN)=IX(J)
    IX(J)=K
!******STEP Q7, PUSH THE GREATEST PARTITION ON STACK
    IF(MAX-J.GT.J-MIN) GOTO 350
!      J-MIN GREATEST
    IF(J-MIN.LE.M) GOTO 400
!      PUSH ONLY IF MAX-J>M
    IF(MAX-J.LE.M) GOTO 360
!      BOTH PARTITIONS ARE GREATER THAN M, THE GREATEST IS PUSHED
    IF(IS.GT.MSTACK) GOTO 910
    LOW(IS)=MIN
    IGH(IS)=J-1
    IS=IS+1
!      CONTINUE TO PARTITION THE SMALLEST
310 MIN=J+1
    GOTO 100
!      MAX-J GREATEST
350 IF(MAX-J.LE.M) GOTO 400
!      PUSH ONLY IF J-MIN>M
    IF(J-MIN.LE.M) GOTO 310
!      BOTH PARTITIONS ARE GREATER THAN M, PUSH THE GREATEST
    IF(IS.GT.MSTACK) GOTO 910
    LOW(IS)=J+1
    IGH(IS)=MAX
    IS=IS+1
!      CONTINUE TO PARTITION THE SMALLEST
360 MAX=J-1
    GOTO 100
!******STEP Q8, POP FROM STACK
400 IS=IS-1
    IF(IS.LT.1) GOTO 500
    MIN=LOW(IS)
    MAX=IGH(IS)
    GOTO 100
!******STEP Q9, STRIGHT INSERTION, ONLY NECESSARY IF M>1
500 CONTINUE
    buperr=0
900 RETURN
910 buperr=1001
    GOTO 900
  END SUBROUTINE SORTIN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine ssort DOES NOT WORK
!\begin{verbatim}
  SUBROUTINE SSORT(CMD,NS,INDEX)
!...SORTING a character array, max 40 characters long
! THIS DOES NOT WORK !!!
    implicit none
    CHARACTER CMD(*)*(*)
    integer ns,index(*)
!\end{verbatim}
    CHARACTER STR*40
    integer l,itop,j,j1,j2,k
!
    L=LEN(CMD(1))
    ITOP=1
    INDEX(ITOP)=1
100 ITOP=ITOP+1
    IF(ITOP.GT.NS) GOTO 900
    STR=CMD(ITOP)
    IF(STR(1:L).GE.CMD(INDEX(ITOP-1))) THEN
       INDEX(ITOP)=ITOP
       GOTO 100
    ENDIF
    J1=1
    J2=ITOP
    J=(J1+J2)/2
200 IF(STR(1:L).LT.CMD(INDEX(J))) THEN
       J2=J
    ELSEIF(J.GT.J1) THEN
       J1=J
    ELSE
       J=J2
       GOTO 300
    ENDIF
    IF(J1.NE.J2) THEN
       K=J
       J=(J1+J2)/2
       IF(K.NE.J) GOTO 200
       J=J2
    ENDIF
!...PLACE FOUND
300 CONTINUE
    MOVE: DO K=ITOP-1,J,-1
       INDEX(K+1)=INDEX(K)
    enddo MOVE
    INDEX(J)=ITOP
    GOTO 100
900 RETURN
  END SUBROUTINE SSORT
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine ssort & Sorting characters
!\begin{verbatim}
  SUBROUTINE SSORT2(CMD,NS,INDX)
!...SORTING a character array, max 40 characters long
! it does not change the position of the texts in CMD but return order in ORDER
    implicit none
    CHARACTER CMD(*)*(*)
    integer ns,indx(*)
!\end{verbatim}
    CHARACTER STR*40
    integer j1,j2,first,previous,next,limit
    integer, allocatable, dimension(:) :: order
!
! links has the the index of the CMD in the increasing order
    allocate(order(ns))
    do j1=1,ns
       order(j1)=-j1
    enddo
!    write(*,'(a,20i4)')'SSORT ',ns,(order(j1),j1=1,ns)
    next=1
    first=1
!    write(*,*)'SSORT first quad ',next,': ',trim(cmd(next))
    all: do j1=2,ns
       previous=-1
       next=first
       limit=0
!       write(*,'(a,i3,a,i3,a,a)')'SSORT loop from ',first,' to ',j1,&
!            ' to find place for ',trim(cmd(j1))
       find: do while(next.le.j1)
          limit=limit+1; if(limit.gt.2*ns) stop 'ininite loop'
          if(next.lt.0) then
! there are no more to compare with, this is the last
             order(previous)=j1
!             write(*,'(a,i3,2x,20i3)')'SSORT insert last at ',&
!                  previous,(order(j2),j2=1,ns)
! do not change sign or order(j1)
             cycle all
          endif
!          write(*,'(a,3i3,1x,a,1x,a,1x,a,1x,a)')'SSORT 1:',previous,j1,next,&
!               ' compare ',trim(cmd(j1)),' and ',trim(cmd(next))
          if(cmd(j1).lt.cmd(next)) then
! insert this after previous, copy link to next to order(j1)
             if(previous.lt.0) then
!                write(*,'(a,a,3i3)')'SSORT 2: insert first before ',&
!                  trim(cmd(next)),previous,j1,next
                order(j1)=first; first=j1
             else
!                write(*,'(a,a,"< ",a," >",a,3i3)')'SSORT 2: insert between ',&
!                     trim(cmd(previous)),trim(cmd(j1)),trim(cmd(next)),&
!                     previous,j1,next
                order(j1)=order(previous);
                order(previous)=j1
             endif
!             write(*,'(a,2i3,2x,20i3)')'SSORT 3:',first,j1,(order(j2),j2=1,ns)
             exit find
          endif
!          write(*,'(a,2i3,2x,20i3)')'SSORT 5:',next,j1,(order(j2),j2=1,ns)
! compare with next
          previous=next
          next=order(next)
       enddo find
!       write(*,'(a,2i3,2x,20i3)')'SSORT 6:',first,0,(order(j2),j2=1,ns)
    enddo all
!    write(*,'(a,2i3,2x,20i3)')'SSORT 7:',first,ns,(order(j2),j2=1,ns)
!
!    next=first
!    limit=1
!    do while(next.gt.0)
!       write(*,*)limit,' ',cmd(next)
!       next=order(next)
!       limit=limit+1
!    enddo
! convert to positions ...
    next=first
    limit=1
    do while(next.gt.0)
       j1=next
       next=order(next)
       indx(j1)=limit
       limit=limit+1
    enddo
!    write(*,'(a,2i3,2x,20i3)')'SSORT 9:',first,ns,(order(j2),j2=1,ns)
!    stop 'ssol'
    return
  end SUBROUTINE SSORT2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
! Routines for manipulation of characters
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable logical function ucletter & Check if character is UPPER case
!\begin{verbatim}
  LOGICAL FUNCTION ucletter(ch1)
! returns TRUE if the character is A to Z
    implicit none
    character ch1*1
!\end{verbatim} %+
    if(lge(ch1,'A') .and. lle(ch1,'Z')) then
       ucletter=.TRUE.
    else
       ucletter=.FALSE.
    endif
  END FUNCTION ucletter
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable character function biglet & Convert one character to UPPER case
!\begin{verbatim}
  CHARACTER FUNCTION BIGLET(CHA)
!...CONVERTS ONE CHARACTER FROM LOWER TO UPPER CASE
    implicit none
    CHARACTER*1 CHA
!\end{verbatim} %+
    CHARACTER*1 CHLAST
    PARAMETER (CHLAST='z')
    IF(CHA.GE.'a' .AND. CHA.LE.CHLAST) THEN
       BIGLET=CHAR(ICHAR(CHA)+ICHAR('A')-ICHAR('a'))
    ELSE
       BIGLET=CHA
    ENDIF
    RETURN
  END FUNCTION BIGLET
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine capson & Convert character to UPPER case
!\begin{verbatim}
  SUBROUTINE capson(text)
! converts lower case ASCII a-z to upper case A-Z, no other changes
    implicit none
    character text*(*)
!\end{verbatim}
    integer, parameter :: lowa=ichar('a'),lowz=ichar('z'),&
         iup=ICHAR('A')-ICHAR('a')
    integer i,ich1
    DO i=1,len(text)
       ich1=ichar(text(i:i))
       IF(ich1.ge.lowa .and. ich1.le.lowz) THEN
          text(i:i)=char(ich1+iup)
       ENDIF
    ENDDO
  END SUBROUTINE capson
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable logical function eolch & TRUE is character is empty after ip
!\begin{verbatim}
  LOGICAL FUNCTION EOLCH(STR,IP)
!...End of Line CHeck, TO SKIP SPACES FROM IP. RETURNS .TRUE. IF ONLY SPACES
!....MODIFIED TO SKIP TAB CHARACTERS ALSO
    implicit none
    CHARACTER STR*(*)
    integer ip
    integer, parameter :: ITAB=9
!\end{verbatim}
!
    EOLCH=.FALSE.
    IF(IP.LE.0) IP=1
100 IF(IP.GT.LEN(STR)) GOTO 110
    IF(STR(IP:IP).NE.' ' .AND. ICHAR(STR(IP:IP)).NE.ITAB) GOTO 900
    IP=IP+1
    GOTO 100
110 EOLCH=.TRUE.
900 RETURN
  END FUNCTION EOLCH

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine getrel & Extrtact real or double
!\begin{verbatim}
  SUBROUTINE GETREL(SVAR,LAST,VALUE)
! extract a real from a character
    implicit none
    character svar*(*)
    integer last
    double precision value
!\end{verbatim} %+
    integer isig
    call getrels(svar,last,value,isig)
    return
  END SUBROUTINE GETREL

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine getrem & Extract real skipping trailing ;
!\begin{verbatim}
  SUBROUTINE GETREM(SVAR,LAST,VAL)
! ...IDENTICAL TO GETREL EXCEPT THAT A TERMINATING COMMA "," IS SKIPPED
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER SVAR*(*)
    integer last
    double precision val
!\end{verbatim} %+
    CALL GETREL(SVAR,LAST,VAL)
    IF(BUPERR.NE.0) RETURN
    IF(SVAR(LAST:LAST).EQ.',') LAST=LAST+1
    RETURN
  END SUBROUTINE GETREM

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine getrels & Extract real
!\begin{verbatim}
  SUBROUTINE GETRELS(SVAR,LAST,VALUE,ISIG)
!...DECODES A REAL NUMBER FROM A TEXT
!      IT MAY BE PRECEEDED BY SPACES AND A + OR -
!      THERE MUST BE AT LEAST ONE NUMBER BEFORE OR AFTER A PERIOD
!      THERE MUST BE AT LEAST ONE NUMBER BEFORE AN "E" OR "D"
!      AFTER AN "E" OR "D" THERE MAY BE A + OR - AND MUST BE ONE OR TWO NUMBERS
! 840310 CHANGE TO ALLOW SPACES AFTER A SIGN I.E. + 2.2 IS ALLOWED
! 860201 EXPONENTIAL D ACCEPTED
! 100910 F95 version
! ISIG is zero if no sign, needed to separte terms inside expressions
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    character svar*(*)
    integer last,isig
    double precision value
!\end{verbatim} %+
! EOLCH is declared as logical function in this module so not needed here
!    LOGICAL EOLCH
    CHARACTER CH*1
    integer i,ierr,jerr
    double precision hel,dec,expo
!    INTEGER GPS,GPN
    CONTINUE
    letter: IF(EOLCH(SVAR,LAST)) THEN
       buperr=1031
       GOTO 9000
    ELSEIF(LAST.LT.LEN(SVAR)-2) THEN
       IF(SVAR(LAST:LAST+3).EQ.'NONE') THEN
          VALUE=RNONE
          buperr=0
          GOTO 9000
       ENDIF
    ENDIF letter
    CH=SVAR(LAST:LAST)
    isig=0
    sign: IF(CH.EQ.'-') THEN
       LAST=LAST+1
       ISIG=-1
       IERR=1036
    ELSE
       ISIG=1
       IF(CH.EQ.'+') THEN
          IERR=1035
          LAST=LAST+1
       ELSE
          IERR=0
       ENDIF
    ENDIF sign
! 840310 NEXT LINE ADDED TO ALLOW FOR SPACES AFTER A SIGN
    IF(EOLCH(SVAR,LAST)) THEN
! if nothing after sign return error
       buperr=IERR; GOTO 9000
    ENDIF
    CH=SVAR(LAST:LAST)
    continue
    nodot: IF(CH.NE.'.') THEN
       JERR=GPN(SVAR,LAST,HEL)
       IF(JERR.NE.0) THEN
! keep error code meaning nagative or positive sign
          if(ierr.eq.0) then
             buperr=jerr
          else
             buperr=ierr
          endif
          GOTO 9000
       ELSE
!...      REMOVE POSSIBLE ERROR CODE SET BY NEGATIVE SIGN
          IERR=0
       ENDIF
    ELSE
!...      MARK THAT THERE WHERE NO DIGITS BEFORE THE DECIMAL POINT
       IF(IERR.EQ.0) IERR=1037
       HEL=ZERO
    ENDIF nodot
!...If the next character is a period then decode the decimal part.
!      If there is no numbers after the period, JERR is nonzero.
!      Then check if IERR is nonzero otherwise there is no numbers at all
!      before or after the period. If so return with error code.
    continue
    decimalpoint: if(last.ge.len(svar)) then
       dec=zero
    elseIF(SVAR(LAST:LAST).EQ.'.') THEN
       LAST=LAST+1
       I=LAST
       JERR=GPN(SVAR,LAST,DEC)
       IF(JERR.EQ.0) THEN
          IERR=0
          I=LAST-I
          DEC=DEC/(TEN**I)
       ELSEIF(IERR.EQ.0) THEN
          DEC=ZERO
       ELSE
!...      NO DIGITS BEFORE OR AFTER A DECIMAL POINT
          buperr=1033
          GOTO 9000
       ENDIF
    ELSE
       DEC=ZERO
    ENDIF decimalpoint
    EXPO=ONE
    exponent: if(last.lt.len(svar)) then
       IF(BIGLET(SVAR(LAST:LAST)).EQ.'E' &
            .OR. BIGLET(SVAR(LAST:LAST)).EQ.'D') THEN
          LAST=LAST+1
          IERR=GPS(SVAR,LAST,EXPO)       
          if(ierr.ne.0) then
             buperr=ierr
             GOTO 9000
          ENDIF
          IF(INT(ABS(EXPO)).GT.99) THEN
             buperr=1039
             GOTO 9000
          ENDIF
          I=INT(EXPO)
          EXPO=TEN**I
       ENDIF
    endif exponent
    VALUE=DBLE(ISIG)*(HEL+DEC)*EXPO
9000 continue
    RETURN
  END SUBROUTINE GETRELS

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable integer function gps & Extract real
!\begin{verbatim}
  INTEGER FUNCTION GPS(SVAR,LAST,VALUE)
!...DECODES A NUMBER WITH OR WITHOUT A SIGN
    implicit none
    DOUBLE PRECISION VALUE
    CHARACTER SVAR*(*)
    integer last
!\end{verbatim} %+
    integer ierr,jerr,isig
    CHARACTER SIG*1
    SIG=SVAR(LAST:LAST)
    IF(SIG.EQ.'-') THEN
       LAST=LAST+1
       ISIG=-1
       IERR=1036
    ELSE
       ISIG=1
       IF(SIG.EQ.'+') THEN
          LAST=LAST+1
          IERR=1035
       ELSE
          IERR=1037
       ENDIF
    ENDIF
    JERR=GPN(SVAR,LAST,VALUE)
    IF(JERR.EQ.0) IERR=0
    GPS=IERR
    VALUE=DBLE(ISIG)*VALUE
    RETURN
  END FUNCTION GPS

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable integer function gpn & Extract real without sign
!\begin{verbatim}
  INTEGER FUNCTION GPN(SVAR,LAST,VALUE)
!...DECODES A NUMBER WITHOUT SIGN
!    DOUBLE PRECISION VALUE
    implicit none
    CHARACTER SVAR*(*)
    integer last
    double precision value
!\end{verbatim}
!    DOUBLE PRECISION, parameter :: ZERO=0.0D0,TEN=1.0D1
    integer l,ierr,n
    L=LEN_TRIM(SVAR)
    VALUE=ZERO
    IERR=1034
    digits: DO LAST=LAST,L
       N=ICHAR(SVAR(LAST:LAST))-ICHAR('0')
       IF(N.LT.0 .OR. N.GT.9) GOTO 800
       IERR=0
       VALUE=TEN*VALUE+DBLE(N)
    enddo digits
800 GPN=IERR
    RETURN
  END FUNCTION GPN
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine getint & Extract integer
!\begin{verbatim}
  SUBROUTINE GETINT(SVAR,LAST,IVAL)
!...DECODES AN INTEGER FROM A TEXT
!      IT MAY BE PRECCEDED BY SPACES AND A + OR -
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER SVAR*(*)
    integer last,ival
!\end{verbatim} %+
    integer ierr
    double precision value
!    
    IF(EOLCH(SVAR,LAST)) THEN
!       CALL ST2ERR(1031,'GETINT','LINE EMPTY')
       buperr=1031
    ELSEIF(SVAR(LAST:MIN(LEN(SVAR),LAST+3)).EQ.'NONE') THEN
       IVAL=NONE
       IERR=0
    ELSE
       IERR=GPS(SVAR,LAST,VALUE)
       IF(IERR.EQ.0) THEN
!          IF(VALUE.GT.FLOAT(MAXINT) .OR. VALUE.LT.FLOAT(MININT)) THEN
          IF(VALUE.GT.DBLE(MAXINT) .OR. VALUE.LT.DBLE(MININT)) THEN
!             CALL ST2ERR(1033,'GETINT','TOO LARGE INTEGER VALUE')
             buperr=1042
             IVAL=0
          ELSE
             IVAL=INT(VALUE)
          ENDIF
       ELSE
!          CALL ST2ERR(IERR,'GETINT','NO DIGIT')
          buperr=ierr
          IVAL=0
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE GETINT

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine getinm & Extract integer and trailing ,
!\begin{verbatim}
  SUBROUTINE GETINM(SVAR,LAST,IVAL)
! ...IDENTICAL TO GETINT EXCEPT THAT A TERMINATING COMMA ",", IS SKIPPED
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER SVAR*(*)
    integer last,ival
!\end{verbatim} %+
    CALL GETINT(SVAR,LAST,IVAL)
    IF(BUPERR.NE.0) RETURN
    IF(SVAR(LAST:LAST).EQ.',') LAST=LAST+1
    RETURN
  END SUBROUTINE GETINM

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine getoct & Extract octal number
!\begin{verbatim}
  SUBROUTINE GETOCT(LINE,IP,IVAL)
!...DECODE AN OCTAL NUMBER
    implicit none
    CHARACTER LINE*(*)
    integer ip,ival
!\end{verbatim} %+
    integer ierr,j
    IERR=0
    IF(EOLCH(LINE,IP)) THEN
!       CALL ST2ERR(1031,'GETOCT','LINE EMPTY')
       buperr=1031
    ELSEIF(LINE(IP:IP+3).EQ.'NONE') THEN
       IVAL=NONE
    ELSE
       IERR=1038
       IVAL=0
100    J=ICHAR(LINE(IP:IP))-ICHAR('0')
       IF(J.GE.0 .AND. J.LE.7) THEN
          IERR=0
          IVAL=8*IVAL+J
       ELSE
          GOTO 800
       ENDIF
       IP=IP+1
       GOTO 100
    ENDIF
!800 IF(IERR.NE.0) CALL ST2ERR(IERR,'GETOCT','NO DIGIT')
800 continue
    RETURN
  END SUBROUTINE GETOCT

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gethex & Extract hexadecimal number
!\begin{verbatim}
  SUBROUTINE GETHEX(LINE,IP,IVAL)
!...DECODE A HEXADECIMAL NUMBER
    implicit none
    CHARACTER LINE*(*)
    integer ip,ival
!\end{verbatim}
    integer bug,ierr,isign,idig,maxdig,j
    CHARACTER CH1*1
!
    IERR=0
    ISIGN=0
    IF(EOLCH(LINE,IP)) THEN
!       CALL ST2ERR(1031,'GETHEX','LINE EMPTY')
       buperr=1031
    ELSEIF(LINE(IP:IP+3).EQ.'NONE') THEN
       IVAL=NONE
    ELSE
       IERR=1038
       IVAL=0
       IDIG=0
       MAXDIG=NBITPW/4
100    CH1=LINE(IP:IP)
       IF(LGE(CH1,'0') .AND. LLE(CH1,'9')) THEN
          J=ICHAR(CH1)-ICHAR('0')
          IERR=0
       ELSEIF(LGE(CH1,'A') .AND. LLE(CH1,'F')) THEN
          J=ICHAR(CH1)-ICHAR('A')+10
          IERR=0
       ELSE
          GOTO 800
       ENDIF
       IDIG=IDIG+1
       IF(IDIG.EQ.1 .AND. J.GE.8) THEN
          ISIGN=1
          J=J-8
       ENDIF
       IVAL=16*IVAL+J
       IP=IP+1
       GOTO 100
    ENDIF
!800 IF(IERR.NE.0) CALL ST2ERR(IERR,'GETHEX','NO DIGIT')
800 continue
!    IF(ISIGN.EQ.1) CALL SETB(1,IVAL)
    bug=ival
! wow, set sign bit of an integer? Assume 32 bits ...
    IF(ISIGN.EQ.1) ival=ibset(ival,31)
    write(*,*)'In metlib4 GETHEX: ',ival,bug
    RETURN
  END SUBROUTINE GETHEX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine getname & Extract a species name
!\begin{verbatim}
  subroutine getname(text,ip,name,mode,ch1)
! reading a species name, this should be incorporated in metlib, 
    implicit none
    character text*(*),name*(*),ch1*1
    integer ip,mode
!\end{verbatim} %+
! Always a letter A-Z as first character
! mode=0 is normal, letters, numbers, "." and "_" allowed ?? should . be allowed
! mode=1 used for species names with "/", "+" and "-" allowed also
    integer jp
    ch1=biglet(text(ip:ip))
    if(ch1.lt.'A' .or. ch1.gt.'Z') then
       write(*,17)ichar(ch1),ch1,text(1:24),ip
17     format('GETNAME error: ',i5,' "',a,'" in "',a,'" at ',i4)
       buperr=1101; goto 1000
    endif
    jp=ip
    do while(ip.lt.len(text))
       ip=ip+1
       ch1=biglet(text(ip:ip))
       if(ch1.ge.'A' .and. ch1.le.'Z') goto 100
       if(ch1.ge.'0' .and. ch1.le.'9') goto 100
       if(ch1.eq.'_' .or. ch1.eq.'.') goto 100
       if(mode.eq.1) then
! special for species names
          if(ch1.eq.'/' .or. ch1.eq.'+' .or. ch1.eq.'-') goto 100
       endif
       goto 200
100    continue
    enddo
200 continue
    name=text(jp:ip-1)
1000 continue
    return
  end subroutine getname

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine getext & Extact a text item
!\begin{verbatim}
  SUBROUTINE GETEXT(SVAR,LAST,JTYP,STRING,CDEF,LENC)
!...SVAR SHALL CONTAIN A TEXT. SCAN STARTS AT POSITION LAST.
!      STRING IS SET TO THE FIRST NONBLANK CHARACTER UP TO THE TERMINATOR.
!      CDEF IS A DEFAULT VAUE IF SVAR IS EMPTY.
!      LENC IS THE LENGTH OF THE TEXT IN STRING
!      JTYP DEFINES THE TERMINATION OF A STRING
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 text terminated by space but if first char is ', " up to next ' or "
!      8 text terminated by space but if first char is (, {, [ or < all text
!             until matching ), }, ] or >. Possibly including more ( ) etc.
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER SVAR*(*),CDEF*(*),STRING*(*)
    integer last,jtyp,lenc
!\end{verbatim}
!    
    CHARACTER CH1*1,CH2*1
    character*1, parameter ::  par(4)=['(','{','[','<']
    character*1, parameter :: ipar(4)=[')','}',']','>']
    integer i,j,k,l1,l2,level,ityp
!    LOGICAL EOLCH,SG2ERR
    IF(JTYP.LE.0) THEN
!       CALL ST2ERR(1030,'GETEXT','NO SUCH TYPE OPTION')
       buperr=1030
       GOTO 900
    ENDIF
    IF(JTYP.LE.6) THEN
       ITYP=JTYP+3
    ELSE
       ITYP=10
       CH2=CHAR(JTYP)
    ENDIF
!...INCREMENT LAST BY ONE TO BYPASS TERMINATOR OF COMMAND OR PREVIOUS
!      ANSWER
    LAST=LAST+1
    IF(LAST.LT.1 .OR. LAST.GT.LEN(SVAR)) LAST=LEN(SVAR)+1
!...SKIP BLANKS STARTING FROM THE POSITION AFTER LAST
!      IF LAST OUTSIDE SVAR THEN ASK QUESTION
    I=LAST
    CONTINUE
    IF(EOLCH(SVAR,I)) GOTO 910
!          STRING=CDEF
!          LENC=LEN(CDEF)
!          LAST=I
!          GOTO 900
!       ENDIF
    CH1=SVAR(I:I)
!...IF FIRST CHARACTER IS "," PUT DEFAULT VALUE IF ANY
    IF(CH1.EQ.',') GOTO 910
! handle ITYP=7 and 8 separately
    if(jtyp.eq.7) then
       if(ch1.eq."'") then
          j=index(svar(i+1:),"'")
          if(j.eq.0) then
! no matching ', return whole string, position after last character
             string=svar(i:len_trim(svar))
             last=len_trim(svar)
             lenc=last-i
             buperr=1032
          else
! return string without ', position after last '
             string=svar(i+1:i+j-1)
             last=i+j+1
             lenc=j-1
          endif
       elseif(ch1.eq.'"') then
          j=index(svar(i+1:),'"')
          if(j.eq.0) then
! no matching ", return whole string, position after last character
             string=svar(i:len_trim(svar))
             last=len_trim(svar)
             lenc=last-i
             buperr=1032
          else
! return string without ", position after last "
             string=svar(i+1:i+j-1)
             last=i+j+1
             lenc=j-1
          endif
       endif
       goto 900
    elseif(jtyp.eq.8) then
! check if first character is ( { or [
       do j=1,4
          if(ch1.eq.par(j)) goto 17
       enddo
       write(*,*)'no open parenthesis ',ch1
! if not ( { [ or < continue with original code
       goto 33
! we must scan svar character by character until matching ipar(j)          
17     continue
       level=1
       k=i
       write(*,*)'jtyp 8, found ',par(j),', in position: ',k
20     k=k+1
       if(k.gt.len(svar)) goto 920
       ch1=svar(k:k)
       if(ch1.eq.par(j)) then
! if we find a new ( { [ or < increase level
          level=level+1
       elseif(ch1.eq.ipar(j)) then
          level=level-1
          if(level.eq.0) then
! we have found matching ) } ] or >
             string=svar(i+1:k-1)
             last=k+1
             lenc=k-i-2
             goto 900
          endif
       endif
       goto 20
    endif
!-------------------------------
! here original code continue
33  continue
!...FETCH THE VALUE FROM SVAR
    LAST=I
    L1=0
    L2=0
    GOTO(40,50,60,70,80,70,100),ITYP-3
40  L1=INDEX(SVAR(LAST:),',')
50  L2=INDEX(SVAR(LAST:),' ')
    GOTO 400
!...
60  L1=INDEX(SVAR(LAST:),'.')
70  L2=INDEX(SVAR(LAST:),';')
!...STRING INCLUDING THE ;
    IF(ITYP.EQ.9 .AND. L2.GT.0) L2=L2+1
    GOTO 400
!...
80  L1=LEN(SVAR)
    GOTO 400
100 L2=INDEX(SVAR(LAST:),CH2)
400 IF(L1.GT.0 .AND. L2.GT.0) THEN
       L1=LAST+MIN(L1,L2)-1
    ELSEIF(L1.LE.0 .AND. L2.LE.0) THEN
       L1=LEN(SVAR)+1
    ELSE
       L1=LAST+MAX(L1,L2)-1
    ENDIF
    IF(L1.GT.LAST) THEN
       STRING=SVAR(LAST:MIN(LEN(SVAR),L1-1))
       LENC=L1-LAST
    ELSE
       STRING=' '
       LENC=0
    ENDIF
    LAST=L1
!
900 RETURN
!...SET DEFAULT VALUE
910 IF(CDEF.NE.CNONE) THEN
       STRING=CDEF
       LENC=LEN(CDEF)
!...SET POSITION IN STRING TO POSITION OF ,
       LAST=I
       GOTO 900
    ENDIF
!...NO ANSWER AND NO DEFAULT VALUE, ERROR RETURN
!920 CALL ST2ERR(1032,'GETEXT','TEXT VALUE MISSING')
920 continue
    buperr=1032
    GOTO 900
  END SUBROUTINE GETEXT

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine wrinum & Write a double left justified
!\begin{verbatim}
  SUBROUTINE WRINUM(STR,IP,NNW,JSIGN,VALUE)
!...EDITS A REAL NUMBER INTO STR WITH LEAST NUMBER OF DIGITS
!      NNW IS MAXIMUM NUMBER OF SIGNIFICANT DIGITS (0<NNW<16)
!      JSIGN >0 INDICATES THAT + SIGN SHOULD BE WRITTEN
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER STR*(*)
    integer ip,nnw,jsign
    double precision value
!\end{verbatim} %+
    CHARACTER CSTR*21,CFRMT*12
!    double precision, parameter :: ZERO=0.0D0,TEN=1.0D1,EPS=1.0D-7
    double precision, parameter :: EPS=1.0D-7
    double precision cc,xx
    integer nw,jj,k,nwd
    CSTR=' '
    NW=NNW
    IF(NW.LE.0) NW=1
    IF(NW.GT.15) NW=15
    IF(IP+NW.GT.LEN(STR)) then
       buperr=1100
       goto 9000
    endif
    IF(VALUE.EQ.ZERO) THEN
       IF(JSIGN.GT.0) THEN
          STR(IP:IP+1)='+0'
          IP=IP+2
       ELSE
          STR(IP:IP)='0'
          IP=IP+1
       ENDIF
       GOTO 9000
    ELSEIF(VALUE.LT.ZERO) THEN
       STR(IP:IP)='-'
       IP=IP+1
    ELSEIF(JSIGN.GT.0) THEN
       STR(IP:IP)='+'
       IP=IP+1
    ENDIF
    CC=ABS(VALUE)
    XX=LOG10(CC+MAX(CC*EPS,EPS))
!    K=INT(XX)
    K=INT(XX)+1
    IF(XX.GT.ZERO) K=K+1
! some problems writing 81000000000 as fixed format ...
! This should be handelled by checking the number of zeroes at the end !!
!    write(*,27)k,nw,cc
!27  format('wrinum: ',2i5,1pe20.12)
    IF(NW.GT.2 .AND. (K.GE.NW .OR. K.LT.-2)) THEN
!...FLOATING FORMAT
       WRITE(CFRMT,100)NW+5,NW-1
100    FORMAT('(1P,E',I2,'.',I2,')')
       WRITE(CSTR,CFRMT)CC
       JJ=NW+1
150    IF(CSTR(JJ:JJ).EQ.'0') THEN
          JJ=JJ-1
          GOTO 150
       ENDIF
       IF(CSTR(JJ:JJ).EQ.'.') JJ=JJ-1
       STR(IP:IP+JJ-1)=CSTR(1:JJ)
       STR(IP+JJ:IP+JJ+3)=CSTR(NW+2:NW+5)
       IP=IP+JJ+4
    ELSE
!...FIXED FORMAT
       NWD=NW-K
       WRITE(CFRMT,200)MAX(NW,NWD)+1,NWD
200    FORMAT('(F',I2,'.',I2,')   ')
       WRITE(CSTR,CFRMT)CC
       JJ=MAX(NW,NWD)+1
250    IF(CSTR(JJ:JJ).EQ.'0') THEN
          JJ=JJ-1
          GOTO 250
       ENDIF
       IF(CSTR(JJ:JJ).EQ.'.') JJ=JJ-1
       if(CSTR(1:1).eq.' ') then
! supress any initial space in CSTR, adjust lenght!!
          CSTR(1:)=CSTR(2:)
          jj=jj-1
       endif
       STR(IP:IP+JJ-1)=CSTR(1:JJ)
       IP=IP+JJ
    ENDIF
9000 continue
    RETURN
  END SUBROUTINE WRINUM

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine wriint & Write an inter left justified
!\begin{verbatim}
  subroutine wriint(text,ipos,int)
! write an integer in text from position ipos (left adjusted)
    implicit none
    character text*(*),number*16
    integer ipos,int,jp
!\end{verbatim} %+
    if(int.lt.0) then
       buperr=1200; text(ipos:ipos)='*'; ipos=ipos+1
    elseif(int.eq.0) then
       text(ipos:ipos)='0'; ipos=ipos+1
    else
       write(number,20)int
20     format(i16)
       jp=1
       if(eolch(number,jp)) then
          buperr=1201; goto 1000
       else
          text(ipos:)=number(jp:)
!          write(*,22)'wriint: ',jp,number(jp:),text(ipos:ipos+16-jp)
!22        format(a,i3,'>',a,'< >',a,'<')
          ipos=ipos+17-jp
!          write(*,30)'wriint: ',ipos,jp,' >'//text(1:ipos+5)//'<'
!30        format(a,2i3,a)
       endif
    endif
1000 continue
    return
  end subroutine wriint

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine wrihex & Write a hexadecimal
!\begin{verbatim}
  SUBROUTINE WRIHEX(STR,IVAL)
!...TO WRITE AN INTEGER AS HEXADECIMAL
!    LOGICAL TESTB
    implicit none
    CHARACTER STR*(*)
    integer ival
!\end{verbatim} %+
    integer j,ip,k
    J=IVAL
    IP=0
10  IP=IP+1
    K=0
    write(*,*)'calling testb from wrihex'
!    IF(TESTB(4*IP-3,J)) K=8
!    IF(TESTB(4*IP-2,J)) K=K+4
!    IF(TESTB(4*IP-1,J)) K=K+2
!    IF(TESTB(4*IP,J)) K=K+1
    IF(btest(4*IP-3,J)) K=8
    IF(btest(4*IP-2,J)) K=K+4
    IF(btest(4*IP-1,J)) K=K+2
    IF(btest(4*IP,J)) K=K+1
    IF(K.GT.9) THEN
       STR(IP:IP)=CHAR(K-10+ICHAR('A'))
    ELSE
       STR(IP:IP)=CHAR(K+ICHAR('0'))
    ENDIF
    IF(IP.LT.LEN(STR)) GOTO 10
    RETURN
  END SUBROUTINE WRIHEX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine wrice & Write a long text
!\begin{verbatim}
  subroutine wrice(lut,margl1,margl2,maxl,str)
! writes str on unit lut with left margin largl1 for first line, margl2 for all
! following lines, max length maxl characters (assuming typewriter font)
    implicit none
    integer lut,margl1,margl2,maxl
    character str*(*)
!\end{verbatim} %+
!    
!    character margx*40
    integer lbreak
    lbreak=0
    call wrice2(lut,margl1,margl2,maxl,lbreak,str)
    continue
    return
  end subroutine wrice

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine wrice2 & Write a long text
!\begin{verbatim}
  subroutine wrice2(lut,margl1,margl2,maxl,lbreak,str)
! writes str on unit lut with left margin largl1 for first line, margl2 for all
! following lines, max length maxl characters (assuming typewriter font)
! lbreak>0 for writing math expression, with stricter linebreak rules
! lbreak<0 for breaking only at space
    implicit none
    character str*(*)
    integer lut,margl1,margl2,maxl,lbreak
!\end{verbatim} %+
!
    character margx*40
    integer lend,nend,lbeg
!
    nend=len_trim(str)
    margx=' '
    lbeg=1
    lend=maxl-margl1
    if(margl1.lt.0.or.margl2.lt.0 .or. maxl.lt.margl1.or.maxl.lt.margl2) then
       buperr=1070; goto 1000
    endif
    if(nend.lt.lend) then
       if(margl1.eq.0) then
          write(lut,10)str(1:nend)
10        format(A)
       else
          write(lut,11)margx(1:margl1),str(1:nend)
11        format(A,A)
       endif
    else
       call cwricend(str,lbeg,lend,lbreak)
       if(margl1.eq.0) then
          write(lut,10)str(1:lend)
       else
          write(lut,11)margx(1:margl1),str(1:lend)
       endif
       do while(lend.lt.nend)
          lbeg=lend+1
          lend=min(lbeg+maxl-margl2-1,nend)
          if(lend.lt.nend) call cwricend(str,lbeg,lend,lbreak)
          if(margl2.eq.0) then
             write(lut,10)str(lbeg:lend)
          else
             write(lut,11)margx(1:margl2),str(lbeg:lend)
          endif
       enddo
    endif
1000 continue
    return
  end subroutine wrice2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine cwicend & Find a possible place for linebreak
!\begin{verbatim}
  subroutine cwricend(str,lbeg,lend,lbreak)
! find a possible place for a newline in str going back from lend
! but not bypassing lbeg.  str is a numerical expression.
! lbreak>0 means stricter rules (mathematical expression)
! lbreak<0 means break only at space
    implicit none
    character str*(*)
    integer lbeg,lend,lbreak
!\end{verbatim}
!
    character ch1*1,ch2*1
    integer ip
! lbreak=0 means
! newline possible at space, ;, +, - (but not sign in exponents like E+02)
    findpos: do ip=lend,lbeg,-1
       ch1=str(ip:ip)
       if(ch1.eq.' ' .or. ch1.eq.';') then
          lend=ip
          goto 1000
       elseif(lbreak.ge.0 .and. (ch1.eq.'+' .or. ch1.eq.'-')) then
          ch2=str(ip-1:ip-1)
!          write(*,*)'cwriceend 3: ',ch2,ch1
          if(ch2.eq.'e' .or. ch2.eq.'E' .or. ch2.eq.'d' .or. ch2.eq.'D' .or. &
               ch2.eq.'(' ) then
             continue
          else
! we cannot find a good breakpoint, break the line here
             lend=ip-1
             goto 1000
          endif
       endif
    enddo findpos
! no position found, just cut at lend
1000 continue
    return
  end subroutine cwricend

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
  
!\addtotable logical function isabbr
!\begin{verbatim}
  INTEGER FUNCTION ISABBR(LONG,SHORT,NC)
! This is for comparing user provided phase names with database phase names
! LONG is a phase name read from database file
! SHORT is an array with abbreviated phase names that should be selected
! Seach array SHORT for any that is an abbreviation of LONG
! Abbreviations between each _ allowed
! NOTE: Any - (minus) in SHORT should have been converted to _
    implicit none
    integer nc
    CHARACTER LONG*(*),SHORT(NC)*(*)
!\end{verbatim} %
    character chs
    integer j1,k1,k2,ltrim,fit,slen
    fit=0
    ltrim=len_trim(long)
!    write(*,*)'M4 ISABBR ***********: ',trim(long),nc,ltrim
! loop to compare LONG with all abbreviations in short
    find: do k1=1,nc
       slen=len_trim(short(k1))
       j1=1
!       write(*,*)'M4 abbr: ',trim(short(k1)),slen
       letter: do k2=1,slen
! this is a loop for all characters in short
          chs=short(k1)(k2:k2)
          uscore: if(chs.eq.'-' .or. chs.eq.'_') then
!             write(*,*)'M4 Found "-" in short, skipping to "-" in long',j1
             long1: do j1=j1,ltrim
                if(long(j1:j1).eq.'-' .or. long(j1:j1).eq.'_') exit long1
             enddo long1
!             write(*,*)'M4 Looking for "-": ',j1,ltrim
             j1=j1+1
! there is no _ or - in long, skip this abbreviation
             if(j1.gt.ltrim) cycle find
! found a - in long, compare letter after - in short and long
             cycle letter
          endif uscore
! accept if next character in short is blank (also if first!)
          if(k2.gt.1 .and. chs.eq.' ') then
             fit=k1; exit find
          endif
! compare letter in short(k1)(k2:k2) with long(j1:j1)
!          write(*,*)'M4 Letter: "',chs,'" and "',long(j1:j1),k2,j1
          if(chs.ne.long(j1:j1)) cycle find
          j1=j1+1
       enddo letter
! accept as all slen letters in SHORT match corresponding leters in LONG
       fit=k1
       exit find
    enddo find
1000 continue
!    if(fit.gt.0) then
!       write(*,*)'Accept abbreviation ',trim(short(fit)),' for ',trim(long),fit
!    endif
    isabbr=fit
    return
  end FUNCTION ISABBR

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
! command interpreters
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
  
!\addtotable integer function ncomp & Top command interpreter
!\begin{verbatim}
  INTEGER FUNCTION NCOMP(SVAR,COMM,NC,NEXT)
! SUBROUTINE NCOMP
    implicit none
    integer nc,next,ient
    CHARACTER SVAR*(*),COMM(NC)*(*)
!\end{verbatim} %+
    IENT=1
    ncomp=ncompx(svar,comm,nc,next,ient)
    return
  end FUNCTION NCOMP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable integer function ncomp2 & Level 1 subcommand
!\begin{verbatim}
  INTEGER FUNCTION NCOMP2(SVAR,COMM,NC,NEXT)
! SUBROUTINE NCOMP2
    implicit none
    integer nc,next,ient
    CHARACTER SVAR*(*),COMM(NC)*(*)
!\end{verbatim} %+
    IENT=2
    ncomp2=ncompx(svar,comm,nc,next,ient)
    return
  end FUNCTION NCOMP2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable integer function ncomp3 & Level 2 subcommand
!\begin{verbatim}
  INTEGER FUNCTION NCOMP3(SVAR,COMM,NC,NEXT)
! SUBROUTINE NCOMP3
    implicit none
    integer nc,next,ient
    CHARACTER SVAR*(*),COMM(NC)*(*)
!\end{verbatim} %+
    IENT=3
    ncomp3=ncompx(svar,comm,nc,next,ient)
    return
  end FUNCTION NCOMP3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable integer function ncompx & Actual command interpreter
!\begin{verbatim}
  INTEGER FUNCTION NCOMPX(SVAR,COMM,NC,NEXT,IENT)
! ...TO DECODE A COMMAND
    implicit none
    CHARACTER SVAR*(*),COMM(NC)*(*)
    integer nc,next,ient
!\end{verbatim}
    character LINE*80
    CHARACTER*1 CH1,CHSEP,CHLAST,CHSYS,CHSEP2,CH2,CHHELP
    CHARACTER*1 CHHIST,CHMAC
    LOGICAL EXAKT,YESLOG
    PARAMETER (CHSEP='_',CHSEP2='-',CHLAST='Z',CHSYS='@',CHHELP='?')
    PARAMETER (CHHIST='!',CHMAC='#')
    integer ls,lc,lika,klik,i,j,last,n,nmatc
!
    if(ient.eq.1) then
       YESLOG=.TRUE.
       last=1
    elseif(ient.eq.2) then
       last=1
    else
       LAST=NEXT+1
    endif
3   LS=LEN(SVAR)
    LC=LEN(COMM(1))
    LIKA=0
    KLIK=0
    IF(EOLCH(SVAR,LAST)) GOTO 300
    IF(LAST.GT.LEN(SVAR)) GOTO 300
!...SPECIAL TREATMENT IF FIRST CHARACTER IS CHSYS OR CHHELP
    CH1=SVAR(LAST:LAST)
    IF(CH1.EQ.CHSYS) GOTO 800
    IF(IENT.EQ.1 .AND. CH1.EQ.CHHELP) THEN
       LAST=LAST+1
       HELP: DO LIKA=1,NC
          IF(COMM(LIKA)(1:5).EQ.'HELP ') THEN
!...            SKIP QUESTION FOR COMMAND IF LINE EMPTY
             IF(SVAR(LAST+1:LAST+1).EQ.' ') SVAR(LAST+1:LAST+2)=',,'
             GOTO 300
          ENDIF
       enddo HELP
       LAST=LAST-1
    ENDIF
    LIKA=0
    IF(IENT.EQ.1 .AND. CH1.EQ.CHHIST) THEN
!... A HISTORY COMMAND. RETURN TO 3 IF A PREVIOUS COMMAND SHALL BE EXEC
       CALL NGHIST(SVAR,LAST)
       IF(LAST.EQ.0) GOTO 3
       NCOMPX=0
       NEXT=0
       GOTO 900
    ENDIF
! FIND LAST CHARACTER IN SVAR THAT IS LEGAL IN A COMMAND
! CONVERT TO CAPITAL LETTERS AT THE SAME TIME
    IF(IENT.EQ.1) LHP=LAST
    LS=LS-LAST+1
    L20: DO N=1,LS
       CH1=BIGLET(SVAR(N+LAST-1:N+LAST-1))
       IF(LGE(CH1,'A') .AND. LLE(CH1,CHLAST)) GOTO 15
       IF(LGE(CH1,'0') .AND. LLE(CH1,'9')) GOTO 15
       IF(CH1.EQ.CHSEP .OR. CH1.EQ.CHSEP2) GOTO 15
       GOTO 50
15     LINE(N:N)=CH1
    enddo L20
!...N UNDEFINED AFTER LOOP?
    N=LS+1
50  LS=N-1
    KLIK=0
    NMATC=0
    IF(N.EQ.1) GOTO 300
    COMPERE: DO N=1,NC
       EXAKT=.TRUE.
       J=0
       LETTER: DO I=1,LS
          CH1=LINE(I:I)
          J=J+1
          IF(CH1.EQ.CHSEP .OR. CH1.EQ.CHSEP2) THEN
!...         PREPARE FOR A "-" JOINING COMMAND AND ARGUMENT
             IF(I.GT.NMATC) THEN
                KLIK=N
                NMATC=I
             ELSEIF(I.EQ.NMATC) THEN
                KLIK=-1
             ENDIF
90           CH2=COMM(N)(J:J)
             IF(CH2.EQ.CHSEP .OR. CH2.EQ.CHSEP2) GOTO 100
             EXAKT=.FALSE.
             J=J+1
             IF(J.GT.LC) GOTO 200
             GOTO 90
          ENDIF
          IF(CH1.NE.COMM(N)(J:J)) GOTO 200
100       CONTINUE
       enddo LETTER
!...A COMMAND THAT CAN FIT, IF EXACTLY EQUAL FINISH
       IF(EXAKT) THEN
          IF(J.EQ.LC) GOTO 500
          IF(COMM(N)(J+1:J+1).EQ.' ') GOTO 500
       ENDIF
!...IF LIKA>0 THE COMMAND IS AMBIGUOUS
       IF(LIKA.GT.0) GOTO 910
       LIKA=N
       LAST=I+LAST-1
200    CONTINUE
    enddo COMPERE
!...ALL COMMANDS COMPERED, IF LIKA=0 THERE WAS NO SUCH COMMAND
300 NEXT=LAST
    IF(LIKA.EQ.0 .AND. KLIK.GT.0) THEN
!...      NO MATCHING COMMAND BUT PART BEFORE A - MATCHES
       LIKA=-(NC+KLIK)
       NEXT=NMATC
    ENDIF
    GOTO 510
500 NEXT=I+LAST-1
    YESLOG=.FALSE.
    LIKA=N
510 CONTINUE
!...RETURN FUNCTION VALUE
    IF(IENT.EQ.1) THEN
       NCOMPX=LIKA
       IF(LIKA.GT.0) THEN
          LSTCMD=COMM(LIKA)
          IF(LOGFIL.GT.0 .AND. YESLOG) WRITE(KOU,517)LSTCMD
517       FORMAT('   ... the command in full is ',A)
          CALL CAPSON(LSTCMD)
       ENDIF
!...SAVE HISTORY, do not save empty lines
       IF(LHP.LE.0)LHP=1
       IF(LEN_TRIM(SVAR(LHP:)).GT.0)THEN
          LHL=LHL+1
          IF(LHL.GT.20) THEN
             LHL=1
             LHM=LHM+20
          ENDIF
          HIST(LHL)=SVAR(LHP:)
       ENDIF
    ELSEIF(IENT.EQ.2) THEN
       NCOMPX=LIKA
    ELSE
       NCOMPX=LIKA
    ENDIF
    GOTO 900
!...A SYSTEM COMMAND OR COMMENT, UPDATE ACCOUNT RECORD BEFORE EXECUTION
!       TWO CHSYS means skip this line (comment)
!       CHSYS followed by CHMAC means macro line, skip it
800 IF(SVAR(LAST+1:LAST+1).EQ.CHSYS) GOTO 810
    IF(SVAR(LAST+1:LAST+1).EQ.CHMAC) GOTO 810
    LINE(1:)=SVAR(LAST+1:)
!       CALL CAPSON(LINE)
!       CALL WPAC2
!      CALL UECOM(LINE)
    write(*,*)'Hit return to continue'
    read(*,808)ch1
808 format(a)
!    CALL COMND(LINE)
810 NEXT=0
    LIKA=0
    GOTO 510
900 RETURN
!...AMBIGUOUS
910 LIKA=-LIKA
    GOTO 510
  END FUNCTION NCOMPX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
! Extracting command arguments from a character
!
! There are two groups, those new finishing with x
! the old without final x (which are listed first)
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparid & Superceeded by gparidx
!\begin{verbatim}
  SUBROUTINE GPARID(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
! ask for integer value with default
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER PROMT*(*),SVAR*(*)
    integer last,ival,idef
    EXTERNAL HELP
!\end{verbatim} %+
    CHARACTER SLIN*512
    integer iflag
! chcek for environment variables
    CALL GQXENV(SVAR)
100 CALL GQARID(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
    RETURN
  END SUBROUTINE GPARID

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gpari
!\begin{verbatim}
  SUBROUTINE GPARI_old(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
! ask for integer value woth no default
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*)
    integer last,ival,idef
    EXTERNAL HELP
!\end{verbatim} %+
    CHARACTER SLIN*80
    integer iflag
! check for environment variables
    CALL GQXENV(SVAR)
100 CALL GQARI(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
    RETURN
  END SUBROUTINE GPARI_OLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparr
!\begin{verbatim}
  SUBROUTINE GPARR_old(PROMT,SVAR,LAST,VAL,RDEF,HELP)
! asks for a double with no default
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER PROMT*(*),SVAR*(*)
    integer last
    double precision val,rdef
    EXTERNAL HELP
!\end{verbatim} %+
    CHARACTER SLIN*80
    integer iflag
! check for environment variables
    CALL GQXENV(SVAR)
100 CALL GQARR(PROMT,SVAR,LAST,VAL,RDEF,HELP)
!    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
    RETURN
  END SUBROUTINE GPARR_OLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparrd
!\begin{verbatim}
  SUBROUTINE GPARRD_old(PROMT,SVAR,LAST,VAL,RDEF,HELP)
! ask for a double with default provided
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER PROMT*(*),SVAR*(*)
    integer last
    EXTERNAL HELP
    double precision val,rdef
!\end{verbatim} %+
    CHARACTER SLIN*80
    integer iflag
! ths checks for environment variables
    CALL GQXENV(SVAR)
100 CALL GQARRD(PROMT,SVAR,LAST,VAL,RDEF,HELP)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
    RETURN
  END SUBROUTINE GPARRD_OLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparc
!\begin{verbatim}
  SUBROUTINE GPARC_old(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
! read a character without default
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CDEF*(*),SVAL*(*)
    integer last,jtyp
    EXTERNAL HELP
!\end{verbatim} %+
    CHARACTER SLIN*80
    integer iflag
! this call handles environment variables
    CALL GQXENV(SVAR)
100 CALL GQARC(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
    IF(BUPERR.NE.0) GOTO 900
    SLIN=SVAL(1:max(1,LEN_TRIM(sval)))
!    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARC_OLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparcd
!\begin{verbatim}
  SUBROUTINE GPARCD_old(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
! read a character with default provided
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CDEF*(*),SVAL*(*)
    integer last,jtyp
    EXTERNAL HELP
!\end{verbatim} %+
!
    CHARACTER SLIN*80
    integer iflag
! this call exchanges environment variables for actual variables
    CALL GQXENV(SVAR)
! this is the real interactive call
100 CALL GQARCD(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
    IF(BUPERR.NE.0) GOTO 900
    SLIN=SVAL
    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN)
!...the next line was missing ...
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARCD_OLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarrd
!\begin{verbatim}
  subroutine GQARRD(PROMT,SVAR,LAST,VAL,RDEF,HELP)
! read real with default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*)
    integer last,ival
    character*1 str,cdef
    double precision val,rdef
    EXTERNAL HELP
!\end{verbatim} %+
    GPARITYP=3
    GPARWDEF=.TRUE.
    GPARRDEF=RDEF
    call gparall(PROMT,SVAR,LAST,IVAL,val,str,cdef,HELP)
    return
  end SUBROUTINE GQARRD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarr
!\begin{verbatim}
  subroutine GQARR(PROMT,SVAR,LAST,VAL,RDEF,HELP)
! read real without default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*)
    integer last,ival
    EXTERNAL HELP
    double precision val,rdef
    character*1 str,cdef
!\end{verbatim} %+
    GPARITYP=3
    GPARWDEF=.FALSE.
    GPARRDEF=RDEF
    call gparall(PROMT,SVAR,LAST,IVAL,val,str,cdef,HELP)
    return
  end SUBROUTINE GQARR

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarid
!\begin{verbatim}
  SUBROUTINE GQARID(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
! previously subroutine GPARID
!...SVAR SHALL CONTAIN A PARAMETER VALUE. IF EMPTY THE PARAMETER IS ASKED FOR
!      USING PROMT AS OUTPUT STRING. IF NO ANSWER THE VALUE IN DEF IS RETURNED
!      INTEGER VALUES. THE DEFAULT VALUE IS DISPLAYED IN THE PROMT WITHIN
!      SLASHES. THE SAME ROUTINES WITHOUT THE FINAL D DOES NOT DISPALY THE
!      DEFAULT VALUE
!      HELP IS A ROUTINE THAT WRITES AN EXPLAINING MESSAGE.
!      LAST IS THE POSITION OF THE TERMINATOR OF THE FORMER PARAMETER OR
!      COMMAND, DECODING STARTS FROM THE POSITION AFTER LAST
    implicit none
    CHARACTER PROMT*(*),SVAR*(*)
    integer last,ival,idef
    character*1 str,cdef
    double precision val
    EXTERNAL HELP
!\end{verbatim} %+
    GPARITYP=1
    GPARWDEF=.TRUE.
    GPARIDEF=IDEF
    call gparall(PROMT,SVAR,LAST,IVAL,val,str,cdef,HELP)
    return
  end SUBROUTINE GQARID

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqari
!\begin{verbatim}
  subroutine GQARI(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
! read integer with no default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*)
    integer last,ival,idef
    character*1 str,cdef
    double precision val
    EXTERNAL HELP
!\end{verbatim} %+
    GPARITYP=1
    GPARWDEF=.FALSE.
    GPARIDEF=IDEF
    call gparall(PROMT,SVAR,LAST,IVAL,val,str,cdef,HELP)
    return
  end SUBROUTINE GQARI

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarcd
!\begin{verbatim}
  subroutine GQARCD(PROMT,SVAR,LAST,JTYP,STR,CDEF,HELP)
! TO READ A STRING VALUE with default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),str*(*),cdef*(*)
    integer last,jtyp
    EXTERNAL HELP
!\end{verbatim} %+
!...SUBROUTINE GQARCD
!      JTYP DEFINES THE TERMINATION OF A STRING
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
    integer ival
    double precision val
    GPARWDEF=.TRUE.
    IF(JTYP.LE.0) THEN
!       CALL ST2ERR(1030,'GPARC','NO SUCH TYPE OPTION')
       buperr=1030
       GOTO 900
    ENDIF
    IF(JTYP.EQ.7) THEN
       GPARENTES=.TRUE.
       GPARITYP=4
    ELSE
       GPARENTES=.FALSE.
       IF(JTYP.LE.6) THEN
! NOTE GPARITYP 1 and 3 used for integer and double precision !!!
          GPARITYP=JTYP+3
       ELSE
          GPARITYP=10
          GPARCH2=CHAR(JTYP)
       ENDIF
    ENDIF
    call gparall(PROMT,SVAR,LAST,IVAL,val,str,cdef,HELP)
    return
900 continue
  end SUBROUTINE GQARCD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparall
!\begin{verbatim}
  SUBROUTINE gparall(PROMT,SVAR,LAST,IVAL,val,string,cdef,HELP)
! previously subroutine GPARID
!...SVAR SHALL CONTAIN A PARAMETER VALUE. IF EMPTY THE PARAMETER IS ASKED FOR
!      USING PROMT AS OUTPUT STRING. IF NO ANSWER THE VALUE IN DEF IS RETURNED
!      INTEGER VALUES. THE DEFAULT VALUE IS DISPLAYED IN THE PROMT WITHIN
!      SLASHES. THE SAME ROUTINES WITHOUT THE FINAL D DOES NOT DISPALY THE
!      DEFAULT VALUE
!      HELP IS A ROUTINE THAT WRITES AN EXPLAINING MESSAGE.
!      LAST IS THE POSITION OF THE TERMINATOR OF THE FORMER PARAMETER OR
!      COMMAND, DECODING STARTS FROM THE POSITION AFTER LAST
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),CH1*1,CDEF*(*),STRING*(*),SSD*30
    CHARACTER PPROMT*132,CH2*1
!    LOGICAL EOLCH,SG2ERR,WDEF,MATP
    LOGICAL WDEF,MATP
    EXTERNAL HELP
!\end{verbatim} %+
    integer last,ival
    double precision val
! local variables
    integer i,ijp,j,jjp,l1,l2,llq,llp,llz,m,kxy,iqq,ityp,idef
    double precision rdef,x
! All routines converge here, update command level and save promt 
! for use by help routines
    CONTINUE
! check if promt already in path, otherwise increase level
    do iqq=2,helprec%level
       kxy=min(len_trim(promt),12)
       if(promt(1:kxy).eq.helprec%cpath(iqq)(1:kxy)) then
          helprec%level=iqq
          goto 991
       endif
    enddo
    if(helprec%level.lt.maxhelplevel) then
       helprec%level=helprec%level+1
       helprec%cpath(helprec%level)=promt
!       write(*,*)'Help level increased to: ',helprec%level
    else
!       write(*,*)'Warning, too many levels in help path'
! This can happen when asking for constitution including the constituent
! just save the last questions
       helprec%cpath(helprec%level)=promt
    endif
991 continue
!-------------------------------
! extract values from calling routines stored in GPARxyz
    wdef=gparwdef
    if(gparityp.eq.1) then
! calling routine wants an integer value
       ityp=1
!       if(wdef) idef=gparidef
! always set the default, it may be used anyway!!
       idef=gparidef
    elseif(gparityp.eq.3) then
! calling routine wants a double precision value
       ityp=3
!       if(wdef) rdef=gparrdef
       rdef=gparrdef
    else
! calling routine wants a string, ityp>3, 
       ityp=gparityp
       matp=gparentes
       ch2=' '
       if(ityp.eq.10) ch2=gparch2
    endif
!-------------------------------------------------------
!...INCREMENT LAST BY ONE TO BYPASS TERMINATOR OF COMMAND OR PREVIOUS ANSWER
    LAST=LAST+1
    IF(LAST.LT.1 .OR. LAST.GT.LEN(SVAR)) LAST=LEN(SVAR)+1
!...SKIP BLANKS STARTING FROM THE POSITION AFTER LAST
!      IF LAST OUTSIDE SVAR THEN ASK QUESTION
    I=LAST
10  CONTINUE
    IF(EOLCH(SVAR,I)) THEN
!...      EMPTY STRING, IF NO BYTES IN PROMT TAKE DEFAULT VALUES
       M=LEN_TRIM(PROMT)
       IF(M.LT.1) GOTO 910
       IF(WDEF) THEN
          PPROMT=PROMT(1:M)//' /'
          JJP=M+3
!...         INSERT DEFAULT VALUE INTO PROMT
          IF(ITYP.EQ.1) THEN
             X=REAL(IDEF)
             CALL WRINUM(PPROMT,JJP,10,0,X)
          ELSEIF(ITYP.EQ.3) THEN
!             IF(RDEF.NE.RNONE) CALL WRINUM(PPROMT,JJP,10,0,RDEF)
! to avoid getting values as 0.0250000004  rather than 0.025 ...
!             IF(RDEF.NE.RNONE) CALL WRINUM(PPROMT,JJP,8,0,RDEF)
             CALL WRINUM(PPROMT,JJP,8,0,RDEF)
          ELSE
             PPROMT(JJP:)=CDEF
             I=LEN_TRIM(CDEF)
             JJP=JJP+I
          ENDIF
          IF(LEN_TRIM(PPROMT).GT.M+2) THEN
             PPROMT(JJP:)='/: '
             JJP=JJP+2
          ELSE
!...            TO AVOID AN EMPTY STRING BETWEEN SLASHES
!             PPROMT(M+1:M+2)=': '
             PPROMT(M+1:)=': '
             JJP=len_trim(ppromt)
          ENDIF
          CALL BOUTXT(KOU,PPROMT(1:JJP))
          IJP=JJP
       ELSE
          PPROMT=PROMT
          JJP=LEN_TRIM(PPROMT)
          CALL BOUTXT(KOU,PPROMT(1:MAX(1,JJP)))
       ENDIF
       SVAR=' '
       CALL BINTXT(KIU,SVAR(1:MIN(LEN(SVAR),130)))
       if(jecho.ne.0) then
! echo ....
          j=len_trim(svar)
          if(j.gt.0) then
             write(kou,77)svar(1:len_trim(svar))
77           format('... echo: ',a)
          endif
       endif
!...WRITE INPUT ON LOG FILE IF ANY
       IF(LOGFIL.GT.0) THEN
!...write on logfile only if question asked !!!
          I=LEN_TRIM(SVAR)
          IF(I.LE.0 .AND. WDEF) THEN
             IF(IJP-3.GE.M+3) THEN
!...WRITE THE DEFAULT ANSWER IF USER INPUT IS EMPTY
                WRITE(LOGFIL,17)PPROMT(M+3:IJP-3)
             ELSE
                WRITE(LOGFIL,17)' '
             ENDIF
          ELSE
             WRITE(LOGFIL,17)SVAR(1:MAX(1,I))
          ENDIF
17        FORMAT(A)
       ENDIF
       I=1
       IF(EOLCH(SVAR,I)) GOTO 910
    ENDIF
!...DECODE THE ANSWER IN SVAR
    CH1=SVAR(I:I)
!...IF FIRST CHARACTER IS "," PUT DEFAULT VALUE IF ANY
    IF(CH1.EQ.',') GOTO 910
!...IF FIRST CHARACTER IS '?' WRITE HELP MESSAGE
    IF(CH1.EQ.'?') THEN
       M=LEN_TRIM(PROMT)
       CALL HELP(PROMT(1:M),SVAR(I:))
       IF(SVAR(I:I+1).EQ.'?!') THEN
!...         THE SPECIAL ROUTINE TOPHLP SHOULD BE USED WHEN HELP IS
!            PROVIDED INSIDE THE CALLING ROUTINE! A ? IS RETURNED
          SVAR(I+1:I+1)=' '
          if (ityp.gt.3) STRING=SVAR(I:)
          GOTO 900
       ENDIF
       I=LEN(SVAR)
       GOTO 10
    ENDIF
!...FETCH THE VALUE FROM SVAR
    LAST=I
    IF(ITYP.EQ.1) THEN
       CALL GETINT(SVAR,LAST,IVAL)
    ELSEIF(ITYP.EQ.3) THEN
       CALL GETREL(SVAR,LAST,VAL)
    ELSEIF(ITYP.LE.10) THEN
!...THE PART HERE IS GLITCHY AS ICE ... handling character input
       L1=0
       L2=0
!  ITYP=      4   5   6   7   8   9  10
       GOTO( 40, 50, 60, 70, 80, 70,100),ITYP-3
! terminate with space or ,
40     L1=INDEX(SVAR(LAST:),',')
! terminate with space
50     L2=INDEX(SVAR(LAST:),' ')
! handle if space or , inside parenthesis should be ignored, like x(fcc,cr)
       IF(MATP) THEN
!...A , OR SPACE INSIDE PARENTHESIS SHOULD BE IGNORED
          LLQ=MIN(L1,L2)
          IF(LLQ.EQ.0) LLQ=MAX(L1,L2)
          LLP=INDEX(SVAR(LAST:),'(')
          IF(LLP.GT.0 .AND. LLP.LT.LLQ) THEN
!...	LLP SHALL BE POSITION OF (, FDMTP UPDATES LLP TO POSITION AFTER )
51           CALL FDMTP(SVAR(LAST:),LLP,SSD)
             IF(BUPERR.NE.0) GOTO 900
             IF(ITYP.EQ.4) L1=INDEX(SVAR(LAST+LLP-1:),',')
             L2=INDEX(SVAR(LAST+LLP-1:),' ')
             IF(L1.GT.0) L1=L1+LLP-1
             IF(L2.GT.0) L2=L2+LLP-1
             LLQ=MIN(L1,L2)
             IF(LLQ.EQ.0) LLQ=MAX(L1,L2)
             LLZ=INDEX(SVAR(LAST+LLP-1:),'(')
             IF(LLZ.GT.0 .AND. LLP+LLZ.LT.LLQ) THEN
!...		   WE HAVE MORE THAN ONE ( BEFORE , OR SPACE
                LLP=LLP+LLZ-1
                GOTO 51
             ENDIF
          ENDIF
       ENDIF
       GOTO 400
! terminale with period
60     L1=INDEX(SVAR(LAST:),'.')
! terminate with semicolon
70     L2=INDEX(SVAR(LAST:),';')
!...STRING INCLUDING THE ;
       IF(ITYP.EQ.9 .AND. L2.GT.0) L2=L2+1
       GOTO 400
!...WHOLE STRING
80     L1=LEN(SVAR)
       GOTO 400
! terminate with provided character >31
100    L2=INDEX(SVAR(LAST:),CH2)
400    IF(L1.GT.0 .AND. L2.GT.0) THEN
          L1=LAST+MIN(L1,L2)-1
       ELSEIF(L1.LE.0 .AND. L2.LE.0) THEN
          L1=LEN(SVAR)+1
       ELSE
!...         BUG FOUND HERE: IF L1 IS LEN(SVAR) AND LAST>1 L1 SET >LEN(SVAR)+1
          L1=MIN(LAST+MAX(L1,L2)-1,LEN(SVAR)+1)
       ENDIF
       IF(L1.GT.LAST) THEN
          STRING=SVAR(LAST:L1-1)
       ELSE
          STRING=' '
       ENDIF
       LAST=L1
    ELSE
!       CALL ST2ERR(1030,'GPAR','NO SUCH TYPE OPTION')
       buperr=1030
    ENDIF
900 RETURN
!...SET DEFAULT VALUE
910 continue
    IF(ITYP.EQ.1) THEN
       IF(IDEF.EQ.NONE) GOTO 920
       IVAL=IDEF
    ELSEIF(ITYP.EQ.3) THEN
       IF(RDEF.EQ.RNONE) GOTO 920
       VAL=RDEF
    ELSE
!       write(*,911)'ML gparqall: ',trim(cdef),trim(cnone),wdef
!911    format(a,a,' - ',a,l2)
       IF(CDEF.NE.CNONE) THEN
          STRING=CDEF
!!       endif
       ELSE
          GOTO 920
       endif
    ENDIF
!...SET POSITION IN STRING TO POSITION OF ,
    LAST=I
    GOTO 900
!...NO ANSWER AND NO DEFAULT VALUE, ERROR RETURN
!920 CALL ST2ERR(1032,'GPAR','PARAMETER VALUE MISSING')
920 buperr=1032
    GOTO 900
  END SUBROUTINE GPARALL

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarc
!\begin{verbatim}
  subroutine GQARC(PROMT,SVAR,LAST,JTYP,STR,CDEF,HELP)
! read a string without default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*)
    integer last,ival,jtyp
    EXTERNAL HELP
    double precision val
    character str*(*),cdef*(*)
!\end{verbatim} %+
    GPARWDEF=.FALSE.
    IF(JTYP.LE.0) THEN
!       CALL ST2ERR(1030,'GPARC','NO SUCH TYPE OPTION')
       buperr=1030
       GOTO 900
    ENDIF
    IF(JTYP.EQ.7) THEN
       GPARENTES=.TRUE.
       GPARITYP=4
    ELSE
       GPARENTES=.FALSE.
       IF(JTYP.LE.6) THEN
          GPARITYP=JTYP+3
       ELSE
          GPARITYP=10
          GPARCH2=CHAR(JTYP)
       ENDIF
    ENDIF
    call gparall(PROMT,SVAR,LAST,IVAL,val,str,cdef,HELP)
900 continue
    return
  end SUBROUTINE GQARC

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparfile
!\begin{verbatim}
    SUBROUTINE GPARFILE(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,TYP,HELP)
! to ask for a file name using command line or external window
! prompt is question
! svar is a character variable which may already contain an answer
! last is position in svar to start searching for an answer
!      JTYP DEFINES THE TERMINATION OF A STRING
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
! sval is the answer either extracted from SVAR or obtained by user input
! cdef is a default answer
! typ  is default file extenion, at present only:
!  1=".TDB", 2=".UNF", 3=".OCM"
! help is a help routine    
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CDEF*(*),SVAL*(*)
    integer last,jtyp
    EXTERNAL HELP
!\end{verbatim}
    CHARACTER SLIN*80
    integer typ,typeahead,kk,iflag
    logical beware
#ifdef tinyfd
! only if we use tinyfiledialogs, check if any character after last+1
    typeahead=last+1
    beware=.FALSE.
! beware set to TRUE if no typeahead (there are non-blanks after positon last+1)
    beware=eolch(svar,typeahead)
!    write(*,*)'M3 gparfile: ',kou,koud,last,eolch(svar,last)
    if(nopenpopup .or. kiu.ne.kiud .or. .not.beware) then
#endif
! If we are not connected to a terminal (reading a macro file) use line input
! Also if there are "type ahead" use the line input
! This call exchanges any macro variables in SVAR for defined macro values
       CALL GQXENV(SVAR)
! If interactive
       if(kiu.eq.kiud .and. beware) write(kou,"(a)") &
            'Beware: you must give the full path unless the file '//&
            'is in working directory!'
100    CALL GQARC(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
       IF(BUPERR.NE.0) GOTO 900
       SLIN=SVAL(1:max(1,LEN_TRIM(sval)))
! This call handles ? @ and other things in SVAR
       CALL GPTCM2(IFLAG,SVAR,LAST,SLIN)
       IF (IFLAG.NE.0) GOTO 100
       if(IUMACLEVL.ge.1) then
          if(sval(1:2).eq.'./') then
! we are running a macro and if SVAL(1:2) is './' replace this with MACROPATH'
             sval=trim(macropath(IUMACLEVL))//sval(3:)
          elseif(sval(1:3).eq.'../') then
! we are running a macro and if SVAL(1:3) is '../' prefix with MACROPATH'
             sval=trim(macropath(IUMACLEVL))//sval
!             write(*,*)'M3 add path: ',trim(sval),IUMACLEVL
!          else
!             write(*,*)'M3 assuming full path or in working directory: '
          endif
       endif
#ifdef tinyfd
    else
! open a popup window to browse directories and files using tinyfiledialogs
! typ<0 means new or old file; 0 old file no filer, 
! typ >0 means old file with filter:
! typ=1 TDB, 2=OCU, 3=OCM, 4=OCD, 5=plt, 6=PDB, 7=DAT
       call getfilename(typ,sval)
       if(sval(1:1).eq.' ') then
          buperr=1020
       elseif(typ.eq.-7) then
! this is for output and file created, if no extension add DAT
          kk=index(sval,'.DAT ')
          if(kk.eq.0) then
             sval(len_trim(sval)+1:)='.DAT'
          endif
       endif
    endif
#endif    
900 RETURN
  END SUBROUTINE GPARFILE

!/!\!/!\!/!\!/!\!/!\!/!\! new X routines /!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
! Extracting command arguments from a character
!
! This is second group with new routines
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparidx & Ask for integer with default
!\begin{verbatim}
  SUBROUTINE GPARIDx(PROMT,SVAR,LAST,IVAL,IDEF,hyper)
! ask for integer value with default
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),hyper*(*)
    integer last,ival,idef
!    EXTERNAL HELP
!\end{verbatim} %+
    CHARACTER SLIN*512
    integer iflag
! chcek for environment variables
    CALL GQXENV(SVAR)
100 CALL GQARIDx(PROMT,SVAR,LAST,IVAL,IDEF,hyper)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
    RETURN
  END SUBROUTINE GPARIDX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparix & Ask for integer no default
!\begin{verbatim}
  SUBROUTINE GPARIx(PROMT,SVAR,LAST,IVAL,IDEF,hyper)
! ask for integer value woth no default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),hyper*(*)
    integer last,ival,idef
!    EXTERNAL HELP
!\end{verbatim} %+
    CHARACTER SLIN*80
    integer iflag
! check for environment variables
    CALL GQXENV(SVAR)
100 CALL GQARIx(PROMT,SVAR,LAST,IVAL,IDEF,hyper)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
    RETURN
  END SUBROUTINE GPARIX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparrx & Ask for double no default
!\begin{verbatim}
  SUBROUTINE GPARRx(PROMT,SVAR,LAST,VAL,RDEF,hyper)
! asks for a double with no default
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),hyper*(*)
    integer last
    double precision val,rdef
!    EXTERNAL HELP
!\end{verbatim} %+
    CHARACTER SLIN*80
    integer iflag
! check for environment variables
    CALL GQXENV(SVAR)
100 CALL GQARRx(PROMT,SVAR,LAST,VAL,RDEF,hyper)
!    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
    RETURN
  END SUBROUTINE GPARRX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparrdx & Ask for double with default
!\begin{verbatim}
  SUBROUTINE GPARRDx(PROMT,SVAR,LAST,VAL,RDEF,hyper)
! ask for a double with default provided
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),hyper*(*)
    integer last
!    EXTERNAL HELP
    double precision val,rdef
!\end{verbatim} %+
    CHARACTER SLIN*80
    integer iflag
! ths checks for environment variables
    CALL GQXENV(SVAR)
100 CALL GQARRDx(PROMT,SVAR,LAST,VAL,RDEF,hyper)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
    RETURN
  END SUBROUTINE GPARRDX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparcdx & Ask for character with default
!\begin{verbatim}
  SUBROUTINE GPARCDx(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,hyper)
! read a character with default provided
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CDEF*(*),SVAL*(*),hyper*(*)
    integer last,jtyp
    EXTERNAL HELP
!\end{verbatim} %+
!
    CHARACTER SLIN*80
    integer iflag
! this call exchanges environment variables for actual variables
    CALL GQXENV(SVAR)
! this is the real interactive call
100 continue
!    CALL GQARCDX(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP,hyper)
    CALL GQARCDX(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,hyper)
    IF(BUPERR.NE.0) GOTO 900
    SLIN=SVAL
    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN)
!...the next line was missing ...
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARCDX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparcx & Ask for character no default
!\begin{verbatim}
  SUBROUTINE GPARCX(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,hyper)
! read a character with default provided and hypertarget
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CDEF*(*),SVAL*(*),hyper*(*)
    integer last,jtyp
!    EXTERNAL HELP now always use Q4HELP
!\end{verbatim} %+
!
    CHARACTER SLIN*80
    integer iflag
! this call exchanges environment variables for actual variables
    CALL GQXENV(SVAR)
! this is the real interactive call ... well no longer ...
!100 CALL GQARCX(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP,hyper)
100 continue
    CALL GQARCX(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,hyper)
    IF(BUPERR.NE.0) GOTO 900
    SLIN=SVAL
    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN)
!...the next line was missing ...
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARCX

!\addtotable subroutine gqaridx & Ask for integer with default
!\begin{verbatim}
  SUBROUTINE GQARIDX(PROMT,SVAR,LAST,IVAL,IDEF,hyper)
!...SVAR SHALL CONTAIN A PARAMETER VALUE. IF EMPTY THE PARAMETER IS ASKED FOR
!      USING PROMT AS OUTPUT STRING. IF NO ANSWER THE VALUE IN DEF IS RETURNED
!      INTEGER VALUES. THE DEFAULT VALUE IS DISPLAYED IN THE PROMT WITHIN
!      SLASHES. THE SAME ROUTINES WITHOUT THE FINAL D DOES NOT DISPALY THE
!      DEFAULT VALUE
!      hyper is a hypertarget for online help
!      LAST IS THE POSITION OF THE TERMINATOR OF THE FORMER PARAMETER OR
!      COMMAND, DECODING STARTS FROM THE POSITION AFTER LAST
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),hyper*(*)
    integer last,ival,idef
    character*1 str,cdef
    double precision val
!\end{verbatim} %+
    GPARITYP=1
    GPARWDEF=.TRUE.
    GPARIDEF=IDEF
    call gparallx(PROMT,SVAR,LAST,IVAL,val,str,cdef,hyper)
    return
  end SUBROUTINE GQARIDX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarix & Ask for integer no default
!\begin{verbatim}
  subroutine GQARIx(PROMT,SVAR,LAST,IVAL,IDEF,hyper)
! read integer with no default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),hyper*(*)
    integer last,ival,idef
!    EXTERNAL HELP
!\end{verbatim} %+
! dummy variables for gparallx
    character*1 str,cdef
    double precision val
    GPARITYP=1
    GPARWDEF=.FALSE.
    GPARIDEF=IDEF
    call gparallx(PROMT,SVAR,LAST,IVAL,val,str,cdef,hyper)
    return
  end SUBROUTINE GQARIX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarrdx & Ask for double with default
!\begin{verbatim}
  subroutine GQARRDx(PROMT,SVAR,LAST,VAL,RDEF,hyper)
! read real with default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),hyper*(*)
    integer last
    double precision val,rdef
    EXTERNAL HELP
!\end{verbatim} %+
    character*1 str,cdef
    integer ival
    GPARITYP=3
    GPARWDEF=.TRUE.
    GPARRDEF=RDEF
    call gparallx(PROMT,SVAR,LAST,IVAL,val,str,cdef,hyper)
    return
  end SUBROUTINE GQARRDX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarrx & Ask for double no default
!\begin{verbatim}
  subroutine GQARRx(PROMT,SVAR,LAST,VAL,RDEF,hyper)
! read real without default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),hyper*(*)
    integer last
    double precision val,rdef
!\end{verbatim} %+
    character*1 str,cdef
    integer ival
    GPARITYP=3
    GPARWDEF=.FALSE.
    GPARRDEF=RDEF
    call gparallx(PROMT,SVAR,LAST,IVAL,val,str,cdef,hyper)
    return
  end SUBROUTINE GQARRX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarcdx & Ask for character with default
!\begin{verbatim}
  subroutine GQARCDX(PROMT,SVAR,LAST,JTYP,STR,CDEF,hyper)
! TO READ A STRING VALUE with default
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),str*(*),cdef*(*),hyper*(*)
    integer last,jtyp
!    EXTERNAL HELP no longer needed
!\end{verbatim} %+
!...SUBROUTINE GQARCDX
!      JTYP DEFINES THE TERMINATION OF A STRING
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
    integer ival
    double precision val
    GPARWDEF=.TRUE.
    IF(JTYP.LE.0) THEN
!       CALL ST2ERR(1030,'GPARC','NO SUCH TYPE OPTION')
       buperr=1030
       GOTO 900
    ENDIF
    IF(JTYP.EQ.7) THEN
       GPARENTES=.TRUE.
       GPARITYP=4
    ELSE
       GPARENTES=.FALSE.
       IF(JTYP.LE.6) THEN
! NOTE GPARITYP 1 and 3 used for integer and double precision !!!
          GPARITYP=JTYP+3
       ELSE
          GPARITYP=10
          GPARCH2=CHAR(JTYP)
       ENDIF
    ENDIF
!    call gparallx(PROMT,SVAR,LAST,IVAL,val,str,cdef,HELP,hyper)
    call gparallx(PROMT,SVAR,LAST,IVAL,val,str,cdef,hyper)
    return
900 continue
  end SUBROUTINE GQARCDX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqarcx & Ask for character no default
!\begin{verbatim}
  subroutine GQARCX(PROMT,SVAR,LAST,JTYP,STR,CDEF,hyper)
! TO READ A STRING VALUE with default user hypertext
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),str*(*),cdef*(*),hyper*(*)
    integer last,jtyp
!    EXTERNAL HELP no longer needed
!\end{verbatim} %+
!...SUBROUTINE GQARCx
!      JTYP DEFINES THE TERMINATION OF A STRING
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
    integer ival
    double precision val
    GPARWDEF=.TRUE.
    IF(JTYP.LE.0) THEN
!       CALL ST2ERR(1030,'GPARC','NO SUCH TYPE OPTION')
       buperr=1030
       GOTO 900
    ENDIF
    IF(JTYP.EQ.7) THEN
       GPARENTES=.TRUE.
       GPARITYP=4
    ELSE
       GPARENTES=.FALSE.
       IF(JTYP.LE.6) THEN
! NOTE GPARITYP 1 and 3 used for integer and double precision !!!
          GPARITYP=JTYP+3
       ELSE
          GPARITYP=10
          GPARCH2=CHAR(JTYP)
       ENDIF
    ENDIF
!    write(*,*)'In GQARCX calling gparallx: ',trim(hyper)
!    call gparallx(PROMT,SVAR,LAST,IVAL,val,str,cdef,HELP,hyper)
    call gparallx(PROMT,SVAR,LAST,IVAL,val,str,cdef,hyper)
    return
900 continue
  end SUBROUTINE GQARCX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparallx & Ask for anything
!\begin{verbatim}
  SUBROUTINE gparallx(PROMT,SVAR,LAST,IVAL,val,string,cdef,hyper)
! this is the focal routine for all variants of GPARxyz
!...SVAR shall contain an answer or command. If EMPTY THE answer IS ASKED FOR
!      USING PROMT AS OUTPUT STRING. IF NO ANSWER THE VALUE IN DEF (default) 
!      is returned if a provided.  The routine can return integer, double or
!      character variables. THE DEFAULT VALUE IS DISPLAYED IN THE PROMT WITHIN
!      SLASHES. if no answer and no defualt an error is returned.
!      HELP is no longer a parameter Q4HELP is always used
!      as hypertarget in a HTML file
!      If hyper contains the character ?TOPHLP and the user has typed a single ?
!      the routine returns with this ? and the calling routine can display
!      a menu.  If the user types two ?? the PROMT is used as hypertarget.
!      LAST IS THE current POSITION IN SVAR, it is incremented by one
!      before looking for an answer (skipping the terminator of any previous
!      input.
! REPEAT:
! when called on top level or from a submenu then hyper='?TOPHLP'
! if user types a single ? only menu listed, with ?? use PROMT as target 
    implicit none
    CHARACTER PROMT*(*),SVAR*(*),CH1*1,CDEF*(*),STRING*(*),hyper*(*)
    integer ival
    double precision val
!    EXTERNAL HELP
!\end{verbatim}
    CHARACTER PPROMT*132,CH2*1,ssd*30
    LOGICAL WDEF,MATP
    integer last,unused
! local variables
!    integer i,ijp,j,jjp,l1,l2,llq,llp,llz,m,kxy,iqq,ityp,idef,kk,nw,kl
    integer i,ijp,j,jjp,l1,l2,llq,llp,llz,m,ityp,idef,kk,nw,kl,qz
    double precision rdef,x
    character hypertarget*(40)
    logical once
! All input routines converge here, update command level and save promt 
! for use by help routines
    once=.TRUE.
    CONTINUE
!    write(*,*)'In gparallx: ',trim(promt),' & ',trim(hyper)
! check if promt already in path, otherwise increase level
    unused=0
! skip old helprec stuff ...
!    do iqq=2,helprec%level
! save first 12 characters in helprec%path (obsolete now?)
!       kxy=min(len_trim(promt),12)
!       if(promt(1:kxy).eq.helprec%cpath(iqq)(1:kxy)) then
!          helprec%level=iqq
!          goto 991
!       endif
!    enddo
!    if(helprec%level.lt.maxhelplevel) then
!       helprec%level=helprec%level+1
!       helprec%cpath(helprec%level)=promt
!       write(*,*)'Help level increased to: ',helprec%level
!    else
!       write(*,*)'Warning, too many levels in help path'
! This can happen when asking for constitution including the constituent
! just save the last questions
!       helprec%cpath(helprec%level)=promt
!    endif
!991 continue
!-------------------------------
! extract values from calling routines stored in GPARxyz
! gparityp is a GLOBAL private variable to transfer type of value to return
! Assuming this is not run in parallel...
    wdef=gparwdef
    if(gparityp.eq.1) then
! calling routine wants an integer value
       ityp=1
!       if(wdef) idef=gparidef
! always set the default, it may be used anyway!!
       idef=gparidef
    elseif(gparityp.eq.3) then
! calling routine wants a double precision value
       ityp=3
!       if(wdef) rdef=gparrdef
       rdef=gparrdef
    else
! calling routine wants a string, ityp>3, 
       ityp=gparityp
       matp=gparentes
       ch2=' '
       if(ityp.eq.10) ch2=gparch2
    endif
!-------------------------------------------------------
!...INCREMENT LAST BY ONE TO BYPASS TERMINATOR OF COMMAND OR PREVIOUS ANSWER
    LAST=LAST+1
    IF(LAST.LT.1 .OR. LAST.GT.LEN(SVAR)) LAST=LEN(SVAR)+1
!...SKIP BLANKS STARTING FROM THE POSITION AFTER LAST
!      IF LAST OUTSIDE SVAR THEN ASK QUESTION
    I=LAST
10  CONTINUE
!    write(*,*)'GPARALLX 10: ',trim(promt)
    IF(EOLCH(SVAR,I)) THEN
!...      EMPTY STRING, IF PROMT empty TAKE DEFAULT VALUES
       M=LEN_TRIM(PROMT)
       IF(M.LT.1) GOTO 910
! avoid promt with double ::
!       if(M.GT.2 .and. promt(M:M).eq.promt(M-1:M-1)) M=M-1
       if(M.GT.2 .and. promt(M:M).eq.':') M=M-1
       ijp=0
       IF(WDEF) THEN
          PPROMT=PROMT(1:M)//' /'
          JJP=M+3
!...         INSERT DEFAULT VALUE INTO PROMT
          IF(ITYP.EQ.1) THEN
             X=REAL(IDEF)
             CALL WRINUM(PPROMT,JJP,10,0,X)
          ELSEIF(ITYP.EQ.3) THEN
!             IF(RDEF.NE.RNONE) CALL WRINUM(PPROMT,JJP,10,0,RDEF)
! to avoid getting values as 0.0250000004  rather than 0.025 ...
!             IF(RDEF.NE.RNONE) CALL WRINUM(PPROMT,JJP,8,0,RDEF)
             CALL WRINUM(PPROMT,JJP,8,0,RDEF)
          ELSE
             PPROMT(JJP:)=CDEF
             I=LEN_TRIM(CDEF)
             JJP=JJP+I
          ENDIF
          IF(LEN_TRIM(PPROMT).GT.M+2) THEN
             PPROMT(JJP:)='/: '
             JJP=JJP+2
          ELSE
!...            TO AVOID AN EMPTY STRING BETWEEN SLASHES
!             PPROMT(M+1:M+2)=': '
             PPROMT(M+1:)=': '
             JJP=len_trim(ppromt)
          ENDIF
          CALL BOUTXT(KOU,PPROMT(1:JJP))
          IJP=JJP
       ELSE
          PPROMT=PROMT
          JJP=LEN_TRIM(PPROMT)
          CALL BOUTXT(KOU,PPROMT(1:MAX(1,JJP)))
       ENDIF
       SVAR=' '
       CALL BINTXT(KIU,SVAR(1:MIN(LEN(SVAR),130)))
! THIS IS USER INPUT
!       write(*,*)'GPARALLX 77: ',trim(svar)
!       if(hyper.eq.'TOPHLP') write(*,*)'gparallx input: ',trim(svar)
       if(jecho.ne.0) then
! echo ....
          j=len_trim(svar)
          if(j.gt.0) then
             write(kou,77)svar(1:len_trim(svar))
77           format('... echo: ',a)
          endif
       endif
!...WRITE INPUT ON LOG FILE IF ANY
       IF(LOGFIL.GT.0) THEN
!...write on logfile only if question asked !!!
          I=LEN_TRIM(SVAR)
          IF(I.LE.0 .AND. WDEF) THEN
             IF(IJP-3.GE.M+3) THEN
!...WRITE THE DEFAULT ANSWER IF USER INPUT IS EMPTY and there is a default
                WRITE(LOGFIL,17)PPROMT(M+3:IJP-3)
             ELSE
                WRITE(LOGFIL,17)' '
             ENDIF
          ELSE
             WRITE(LOGFIL,17)SVAR(1:MAX(1,I))
          ENDIF
17        FORMAT(A)
       ENDIF
       I=1
       IF(EOLCH(SVAR,I)) GOTO 910
    ENDIF
!...DECODE THE ANSWER IN SVAR
!    write(*,*)'GPARALLX 17: ',trim(svar),i
    CH1=SVAR(I:I)
!...IF FIRST CHARACTER IS "," PUT DEFAULT VALUE IF ANY
    IF(CH1.EQ.',') GOTO 910
!...IF FIRST CHARACTER IS '?' WRITE HELP MESSAGE
    IF(CH1.EQ.'?') THEN
       if(hyper.eq.'?TOPHLP') then
          if(SVAR(I+1:I+1).ne.'?') then
! if the user types a single ? then just display the menu, no browserhelp!
             last=i
             string=svar(i:)
! The menu is displayed in calling routine, just return here
!             write(*,*)'GPARALLX quick help exit!',trim(string),i
             goto 900
          else
! user types two ??, generate the hypertarget from prompt
             if(promt(1:2).eq.'--') then
! this is ?? typed at top level
                hypertarget='?All commands'
             else
! when prompting for commands the default must be a character
                write(*,*)'gparallx extract: "',trim(promt),'" and "',&
                     trim(cdef),'"'
! extract part of the promt as hypertarget, 
! for a promt "Amend for phase LIQUID what?" extract "Amend for phase"
! to use as hypertarget.  Use only the 3 first words
                kl=1
                kk=1
                nw=0
                max3: do while(kk.gt.0)
! a subsub command should always have 3 fixed words in the promt!!!
! the last word of a subcommand, "what?", has no " " at the end
                   kk=index(promt(kl:),' ')
                   if(kk.gt.0) then
                      kl=kl+kk
                      nw=nw+1
                   else
                      kl=kl-1
                      exit max3
                   endif
                   if(nw.eq.3) then
                      kl=kl-1
                      exit max3
                   endif
                enddo max3
                if(kl.le.1) kl=len_trim(promt)
                hypertarget='?'//promt(1:kl)
! to handle help when user types two ?? for the promt "amend what? /phase/:"
! then include the default answer "phase" in the hypertarget !!
! A single ? already gives the submenu for "amend phase"
                if(cdef(1:1).ne.' ') then
                   hypertarget(kl+2:)=trim(cdef)
                endif
! convert prompt to lower case except first letter
!                qz=len_trim(hypertarget)
                call lowercase1(hypertarget)
!                write(*,*)'gparallx hypertarget: "',&
!                     trim(hypertarget),'" and "',trim(cdef),'"',nw,kl
             endif
!             write(*,*)'GPARALLX hypertarget: ',trim(hypertarget)
          endif
       else
! normal questions have the hypretarget in hyper to provide help
          hypertarget=hyper
       endif
! if two ?? or a real question (not menu) provide more advanced help
       M=LEN_TRIM(PROMT)
! using q4help the arguments should be hyper and an (so far) unused integer
!       CALL HELP(PROMT(1:M),SVAR(I:))
!       CALL HELP(hypertarget,unused)
       write(*,*)'Calling Q4HELP: "',trim(hypertarget),'"',unused
       CALL q4help(hypertarget,unused)
       if(unused.ne.0) then
! return a ? to calling routine ... why?
          SVAR(I+1:I+1)='?'
          if (ityp.gt.3) STRING=SVAR(I:)
          GOTO 900
       ENDIF
       I=LEN(SVAR)
       if(once) then
          once=.false.; GOTO 10
       else
! A ? more than once as answer, quit with error
          goto 920
       endif
    ENDIF
!...FETCH THE VALUE FROM SVAR
    LAST=I
    IF(ITYP.EQ.1) THEN
       CALL GETINT(SVAR,LAST,IVAL)
    ELSEIF(ITYP.EQ.3) THEN
       CALL GETREL(SVAR,LAST,VAL)
    ELSEIF(ITYP.LE.10) THEN
!...THE PART HERE IS GLITCHY AS ICE ... handling character input
       L1=0
       L2=0
!  ITYP=      4   5   6   7   8   9  10
       GOTO( 40, 50, 60, 70, 80, 70,100),ITYP-3
! terminate with space or ,
40     L1=INDEX(SVAR(LAST:),',')
! terminate with space
50     L2=INDEX(SVAR(LAST:),' ')
! handle if space or , inside parenthesis should be ignored, like x(fcc,cr)
       IF(MATP) THEN
!...A , OR SPACE INSIDE PARENTHESIS SHOULD BE IGNORED
          LLQ=MIN(L1,L2)
          IF(LLQ.EQ.0) LLQ=MAX(L1,L2)
          LLP=INDEX(SVAR(LAST:),'(')
          IF(LLP.GT.0 .AND. LLP.LT.LLQ) THEN
!...	LLP SHALL BE POSITION OF (, FDMTP UPDATES LLP TO POSITION AFTER )
51           CALL FDMTP(SVAR(LAST:),LLP,SSD)
             IF(BUPERR.NE.0) GOTO 900
             IF(ITYP.EQ.4) L1=INDEX(SVAR(LAST+LLP-1:),',')
             L2=INDEX(SVAR(LAST+LLP-1:),' ')
             IF(L1.GT.0) L1=L1+LLP-1
             IF(L2.GT.0) L2=L2+LLP-1
             LLQ=MIN(L1,L2)
             IF(LLQ.EQ.0) LLQ=MAX(L1,L2)
             LLZ=INDEX(SVAR(LAST+LLP-1:),'(')
             IF(LLZ.GT.0 .AND. LLP+LLZ.LT.LLQ) THEN
!...		   WE HAVE MORE THAN ONE ( BEFORE , OR SPACE
                LLP=LLP+LLZ-1
                GOTO 51
             ENDIF
          ENDIF
       ENDIF
       GOTO 400
! input terminates with period
60     L1=INDEX(SVAR(LAST:),'.')
! input terminates with semicolon
70     L2=INDEX(SVAR(LAST:),';')
!...STRING INCLUDING THE ;
       IF(ITYP.EQ.9 .AND. L2.GT.0) L2=L2+1
       GOTO 400
!...WHOLE STRING
80     L1=LEN(SVAR)
       GOTO 400
! input terminates with provided character >31
100    L2=INDEX(SVAR(LAST:),CH2)
400    IF(L1.GT.0 .AND. L2.GT.0) THEN
          L1=LAST+MIN(L1,L2)-1
       ELSEIF(L1.LE.0 .AND. L2.LE.0) THEN
          L1=LEN(SVAR)+1
       ELSE
!...         BUG FOUND HERE: IF L1 IS LEN(SVAR) AND LAST>1 L1 SET >LEN(SVAR)+1
          L1=MIN(LAST+MAX(L1,L2)-1,LEN(SVAR)+1)
       ENDIF
       IF(L1.GT.LAST) THEN
          STRING=SVAR(LAST:L1-1)
       ELSE
          STRING=' '
       ENDIF
       LAST=L1
    ELSE
!       CALL ST2ERR(1030,'GPAR','NO SUCH TYPE OPTION')
       buperr=1030
    ENDIF
900 RETURN
!...SET DEFAULT VALUE
910 continue
    IF(ITYP.EQ.1) THEN
       IF(IDEF.EQ.NONE) GOTO 920
       IVAL=IDEF
    ELSEIF(ITYP.EQ.3) THEN
       IF(RDEF.EQ.RNONE) GOTO 920
       VAL=RDEF
    ELSE
!       write(*,911)'ML gparqall: ',trim(cdef),trim(cnone),wdef
!911    format(a,a,' - ',a,l2)
       IF(CDEF.NE.CNONE) THEN
          STRING=CDEF
!!       endif
       ELSE
          GOTO 920
       endif
    ENDIF
!...SET POSITION IN STRING TO POSITION OF ,
    LAST=I
    GOTO 900
!...NO ANSWER AND NO DEFAULT VALUE, ERROR RETURN
!920 CALL ST2ERR(1032,'GPAR','PARAMETER VALUE MISSING')
920 buperr=1032
    GOTO 900
  END SUBROUTINE GPARALLX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine lowercase1 & convert character to lower case
!\begin{verbatim}
  subroutine lowercase1(text)
    character text*(*)
!\end{verbatim}
    integer ip,jp,kp,ichA,ichZ,chlower,ich1
    ichA=ichar('A')
    ichZ=ichar('Z')
! cha = chA+chlower
    chlower=ichar('a')-ichA
!    write(*,*)'ML text: "',trim(text),'"'
! do not convert first 2 characters or any character not between chA and chZ
    do ip=3,len(text)
       ich1=ichar(text(ip:ip))
       if(ich1.ge.ichA .and. ich1.le.ichZ) then
          text(ip:ip)=char(ich1+chlower)
       endif
    enddo
!    write(*,*)'ML text: "',trim(text),'"'
900 return
  end subroutine lowercase1
    
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gparfilex & Ask for file name
!\begin{verbatim}
    SUBROUTINE GPARFILEx(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,TYP,hyper)
! to ask for a file name using command line or external window
! prompt is question
! svar is a character variable which may already contain an answer
! last is position in svar to start searching for an answer
!      JTYP DEFINES THE TERMINATION OF A STRING (maybe redundant??)
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
! sval is the answer either extracted from SVAR or obtained by user input
! cdef is a default answer
! typ  is default file extenion, at present only:
!  1=".TDB", 2=".UNF", 3=".OCM", 4= , 5= , 6= , 7= , 8=".LOG"
! hyper is a hypertext target for help
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CDEF*(*),SVAL*(*),hyper*(*)
    integer last,jtyp
!    EXTERNAL HELP
!\end{verbatim}
    CHARACTER SLIN*256
    integer typ,typeahead,kk,iflag
    logical beware
    sval=' '
    slin=cdef
!    write(*,10)'M4 in gparfilex: ',typ,trim(cdef),trim(sval),trim(slin)
10  format(a,i3,' "',a,'" "',a,'" "',a,'"')
#ifdef tinyfd
! only if we use tinyfiledialogs, check if any character after last+1
    typeahead=last+1
    beware=.FALSE.
! beware set to TRUE if no typeahead (there are non-blanks after positon last+1)
    beware=eolch(svar,typeahead)
!    write(*,*)'M3 gparfile: ',kou,koud,last,eolch(svar,last)
    if(nopenpopup .or. kiu.ne.kiud .or. .not.beware) then
       continue
#endif
! If we are not connected to a terminal (reading a macro file) use line input
! Also if there are "type ahead" use the line input
! This call exchanges any macro variables in SVAR for defined macro values
       CALL GQXENV(SVAR)
! If interactive
       if(kiu.eq.kiud .and. beware) write(kou,"(a)") &
            'Beware: you must give the full path unless the file '//&
            'is in working directory!'
100    CALL GQARCx(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,hyper)
       IF(BUPERR.NE.0) GOTO 900
       SLIN=SVAL(1:max(1,LEN_TRIM(sval)))
! This call handles ? @ and other things in SVAR
       CALL GPTCM2(IFLAG,SVAR,LAST,SLIN)
       IF (IFLAG.NE.0) GOTO 100
       if(IUMACLEVL.ge.1) then
          if(sval(1:2).eq.'./') then
! we are running a macro and if SVAL(1:2) is './' replace this with MACROPATH'
             sval=trim(macropath(IUMACLEVL))//sval(3:)
          elseif(sval(1:3).eq.'../') then
! we are running a macro and if SVAL(1:3) is '../' prefix with MACROPATH'
             sval=trim(macropath(IUMACLEVL))//sval
!             write(*,*)'M3 add path: ',trim(sval),IUMACLEVL
!          else
!             write(*,*)'M3 assuming full path or in working directory: '
          endif
       endif
#ifdef tinyfd
    else
! open a popup window to browse directories and files using tinyfiledialogs
! typ<0 means new or old file; 0 old file no filer, 
! typ >0 means old file with filter:
! typ=1 TDB, 2=OCU, 3=OCM, 4=OCD, 5=plt, 6=PDB, 7=DAT, 8=LOG
!       write(*,*)'M4 opening popup window',typ
       call getfilename(typ,sval)
!       write(*,*)'M4 From getfilename: "',trim(sval),'"',typ
       if(sval(1:1).eq.' ') then
          buperr=1020
       elseif(typ.eq.-7) then
! this is for output and file created, if no extension add DAT
          kk=index(sval,'.DAT ')
          if(kk.eq.0) then
             sval(len_trim(sval)+1:)='.DAT'
          endif
       elseif(typ.eq.-8) then
! this is for output and file created(?), if no extension add LOG
! Check if last 4 letters are none
          kk=len_trim(sval)
          if(kk.ge.4) then
             if(sval(kk-3:kk).eq.'none') then
                sval='NONE'
                goto 300
             endif
          endif
          kk=index(sval,'.LOG ')
!          write(*,*)'gparfilex: ',trim(sval),kk,trim(cdef)
          if(kk.eq.0) then
             iflag=len_trim(sval)
             sval(iflag+1:)='.LOG'
          endif
300       continue
!          write(*,*)'gparfilex: ',trim(sval),kk
       endif
!       write(*,*)'M4 file: "',trim(sval),'"'
    endif
#endif    
! Can the rather odd ifdef/endif cause problems ???
! if there is a segmentation fault it is after this write statement ... SUCK
!    write(*,*)'M4 exit gparfilex: ',trim(sval)
900 RETURN
  END SUBROUTINE GPARFILEX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
! A new set of on-line help routines using browser
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine init_help & Initiate help and history
!\begin{verbatim}
  subroutine init_help(browser,htmlfile)
! This routine is called from oc_command_monitor to inititate
! the on-line help system. It saves the name of the browser and HTML file
    implicit none
    character*(*) htmlfile,browser
!\end{verbatim} %+
    character noquotes*128
    integer kk
    logical logok
! the latex file no longer used for help
    ochelp%latexfile=' '
! test that file exists
    inquire(file=browser,exist=logok)
!    write(*,*)'m4A: ',trim(browser),logok
    if(.not.logok) then
! This is emergency use of explorer if no Firefox
!    browser='C:\PROGRA~1\INTERN~1\iexplore.exe '
       noquotes=browser
       kk=index(noquotes,'"')
       do while(kk.gt.0)
          noquotes(kk:)=noquotes(kk+1:)
          kk=index(noquotes,'"')
       enddo
       inquire(file=noquotes,exist=logok)
!       write(*,*)'m4C: ',trim(noquotes),logok
    endif
    allok: if(logok) then
       ochelp%browser=browser
       inquire(file=htmlfile,exist=logok)
!       write(*,*)'m4B: ',trim(htmlfile),logok
       if(logok) then
          helprec%okinit=1
          helprec%type='html'
          ochelp%htmlhelp=.TRUE.
          ochelp%htmlfile=htmlfile
          goto 1000
       endif
    endif allok
    helprec%okinit=0
    helprec%type=' '
    ochelp%htmlhelp=.FALSE.
    ochelp%htmlfile=' '
    ochelp%browser=' '
1000 continue
    return
  end subroutine init_help

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine helplevel & Redundant
!\begin{verbatim}
  subroutine helplevel1(line)
! This routine is called from the monitor for the top level command
! It initiates the path to find the correct help text
! In all gparx routines the help level in increased and the question saved
    implicit none
    character*(*) line
!\end{verbatim} %+
    helprec%level=1
    helprec%cpath(1)=line
    return
  end subroutine helplevel1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine q1help & Old help routine 1
!\begin{verbatim}
  subroutine q1help(prompt,line)
! This routine is called from all gparx routines 
! when the user types a ?
! prompt is never used ...
    implicit none
    character*(*) prompt,line
    character hline*80,mtext*12
    integer, parameter :: maxlevel=20
!\end{verbatim} %+
!    character subsec(5)*10,saved(maxlevel)*24
    character subsec(5)*10
    character htmlhelp*256
    integer nsaved(maxlevel)
!    integer izz,jj,kk,kkk,level,nl,l2,np1,np2,nsub,zz
    integer izz,jj,kk,level,nl,np1,np2,nsub,zz
    logical foundall
!
    nsaved=0
    subsec(1)='%\section{'
    subsec(2)='%\subsecti'
    subsec(3)='%\subsubse'
    subsec(4)='%\subsubsu'
    subsec(5)='%\question'
    if(helprec%okinit.eq.0) then
       if(helptrace) write(kou,*)'Sorry no help file'
       goto 1000
    endif
! USEFUL for helptraceging list current search path:
    if(helptrace) then
       do nl=1,helprec%level
          write(*,17)'Search level: ',nl,trim(helprec%cpath(nl))
17        format(a,i3,2x,a)
       enddo
    endif
!
    open(31,file=ochelp%latexfile,status='old',access='sequential')
    nl=0
    level=2
    np1=0
    np2=0
    nsub=1
    foundall=.false.
    ochelp%target=' '
    if(helprec%type.ne.'latex   ') then
       write(*,*)'Sorry only help based on LaTeX implemented'
       goto 900
    endif
! plain LaTeX file. The questions the OC software asks are saved from the
! top level in helprec%cpath(1..level).  This makes it possible to compare
! the these commands with comment lines in the help file to find the relevant 
! helptext.  The comment lines are structured as the LaTeX sections
! %\subsection{question1}, %\subsubse..{questione} etc
! for each match the sublevel is increased and when we find match
! with the last helprec%cpath(helprec%level) we assume the text until
! the next %\sub....  can be provided as help
! If there is an additional HTML help file the text can instead be displayed
! in a browser using \hypertarget{label} from the LaTeX file found after
! the last matching sublevel
! Only first 12 characters in helprec%cpath and %\section{ sublevel are used
! return here when we found match at level
100 continue
    level=level+1
    if(helptrace) write(*,*)'At label 100: ',level,helprec%level,nl
    if(level.gt.helprec%level) then
       foundall=.true.
       if(helptrace) write(*,*)'Foundall 1',nl
       goto 200
    elseif(level.eq.helprec%level .and.&
         helprec%cpath(level)(1:2).eq.'? ') then
! this is when help is asked in a submenue with two ??
! with just one ? the menue is displayed, with ?? the helpfile is used
       foundall=.TRUE.
       if(helptrace) write(*,*)'Foundall 2',nl
       goto 200
    endif
110 continue
! skip cpath levels that contain COMMAND: or WHAT?
! if last level and cpath contain ? we have found all
    if(index(helprec%cpath(level),'COMMAND: ').gt.0 .or. &
         index(helprec%cpath(level),' WHAT? ').gt.0) then
       level=level+1
       if(level.gt.helprec%level) then
          foundall=.TRUE.
          if(helptrace) write(*,*)'Foundall 2',nl
          goto 200
       endif
       goto 110
    endif
    if(helptrace) write(*,*)'Searching for: ',trim(helprec%cpath(level)),level
! return here when last line did not contain any matching subsec
! we can arrive here with np1=0 and foundall==true
! for help at first command level
200 continue
    read(31,210,end=700)hline
210 format(a)
    nl=nl+1
    if(np1.gt.0) then
! np1 is nonzero if we have found a line matching one helprec%cpath
! We will save all hypertarget labels to have some idea what help text
! to provide if we do not find all %cpath
! If we found the helprec%cpath(helprec%level) foundall is set TRUE
! but we continue until we find the following %\section at the same
! or higher sublevel
       kk=index(hline,'\hypertarget{')
       if(kk.gt.0) then
          ochelp%target=hline(kk+13:)
       endif
       if(foundall) then
! terminate at a line with any sublevel
          izz=0
          do kk=1,5
             if(hline(1:10).eq.subsec(6-kk)) izz=1
          enddo
          if(izz.gt.0) then
             np2=nl-1
             goto 700
          endif
          goto 200
       endif
    elseif(foundall) then
! this should give help from user guide for section %\section{All commands}
!       write(*,*)'M3 "All commands"',hline(1:24),nl
       if(hline(1:23).eq.'%\section{All commands}') then
          np1=nl
          np2=nl+20
! next line should be hypertarget
          read(31,210)hline
!          write(*,*)'next line: ',trim(hline),nl
          kk=index(hline,'\hypertarget{')
          if(kk.gt.0) then
             ochelp%target=hline(kk+13:)
             kk=index(ochelp%target,'}')
             ochelp%target(kk:)=' '
             goto 700
          else
! the help file is messed up ...             
             ochelp%target='All commands'
             goto 700
          endif
       endif
       goto 200
!    else
! here we now have np1>0 and use the rest of this routine as usual
    endif
! we are searching for a subsec on the sublevel nsub
! Check if we have a %\section of this sublevel on the line
    kk=index(hline,subsec(nsub))
    section: if(kk.eq.0) then
! if there is none but we already found one sublevel check if we find the same
!       write(*,*)'no subsec: ',nsub
       prevsub: if(nsub.gt.2) then
          kk=index(hline,subsec(nsub-2))
          if(kk.gt.0) then
! we have found a sublevel 2 levels up ... we are out of scope          
             if(helptrace) write(*,*)'Found subsec two levels up!'
             np2=nl
             goto 700
          elseif(nsub.gt.1) then
             kk=index(hline,subsec(nsub-1))
             if(kk.gt.0) then
! we have found a subsec at the same sublevel we already found
! check if we have match with the helprec%cpath, only 12 first characters!
                jj=index(hline,'{')
                if(jj.le.0) then
                   write(*,*)'LaTeX helpfil missing { on line:',nl
                   goto 200
                endif
                mtext=hline(jj+1:)
                kk=index(mtext,'}')
                if(kk.gt.0) mtext(kk:)=' '
                call capson(mtext)
                zz=len_trim(mtext)
                if(helptrace) write(*,300)'same: ',helprec%cpath(level)(1:zz),&
                     ' =?= ',mtext(1:zz),level,nsub,nl
                if(helprec%cpath(level)(1:zz).eq.mtext(1:zz)) then
! we have found match with the next level of user path on same sublevel
                   goto 100
                endif
             endif
          endif
       endif prevsub
! just read another line
       goto 200
    else
! we have found a %\sub... for next level, check if it is %cpath(level)
       jj=index(hline,'{')
       if(jj.le.0) then
          write(*,*)'LaTeX helpfil missing { on line:',nl
          goto 200
       endif
       mtext=hline(jj+1:)
       kk=index(mtext,'}')
       if(kk.gt.0) mtext(kk:)=' '
       call capson(mtext)
       zz=len_trim(mtext)
       if(helptrace) write(*,300)'next: ',helprec%cpath(level)(1:zz),' =?= ',&
            mtext(1:zz),level,nsub,nl
300    format(a,a,a,a,5i5)
       if(helprec%cpath(level)(1:zz).eq.mtext(1:zz)) then
! we have found match with the next level of user path
          if(helptrace) write(*,*)'Match: ',level,nsub,nl
          nsub=nsub+1
          np1=nl
          goto 100
       endif
       goto 200
    endif section
! jump here if we do not search any more
! we should write lines from np1 to np2 from help file or HTML file
700 continue
    if(np1.gt.0) then
       if(np2.le.np1) then
! we found no obvious end of help text
          write(*,*)'Help text range error: ',np1,np2
       endif
! if htmlhelp is true open a browser window and place text at target
       htmlfil: if(ochelp%htmlhelp .and. ochelp%target(1:1).ne.' ') then
! the user has to close the help window to continue ... spawn??
!          write(*,711)np1,np2
!711       format(/' *** You must close the browser window to continue OC',2i5/)
! the \hypertaget should be finished by a }
          kk=index(ochelp%target,'}')-1
          if(kk.le.0) kk=len_trim(ochelp%target)
#ifdef lixhlp
! on linux just ' "file:" as ochelp#htmlfile start with a /
! The & at the end spawns the browser window and furter ? creates new tags !!
          htmlhelp=trim(ochelp%browser)//' "file:'//&
               trim(ochelp%htmlfile)//'#'//ochelp%target(1:kk)//'" &'
#else
! on Windows we need the / after file
! the initial start spawns a new window with the browser, each ? a new browser
          htmlhelp='start '//trim(ochelp%browser)//' "file:/'//&
               trim(ochelp%htmlfile)//'#'//ochelp%target(1:kk)//'"'
#endif
          if(helptrace) write(*,*)'MM: ',trim(htmlhelp)
          call execute_command_line(htmlhelp)
          goto 900
       else
! help in user terminal screen: write a blank line
          write(kou,*)
          write(*,798)np1,np2
798       format(' >>> We should open a help window to display text: ',2i5)
          rewind(31)
          nl=0
800       continue
          read(31,210)hline
          nl=nl+1
          if(nl.ge.np2) then
             goto 900
          elseif(nl.ge.np1) then
             if(hline(1:1).ne.'%') then
! ignore LaTeX comment lines and replace \item with a -
                if(hline(2:5).eq.'item') then
                   write(*,811)trim(hline(6:))
811                format('- ',a)
                else
                   write(*,210)trim(hline)
                endif
             endif
          endif
          goto 800
       endif htmlfil
    else
       write(*,*)'No help found'
    endif
900 continue
    close(31)
!
1000 continue
    return
  end subroutine q1help

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine q2help & Old help routine 2
!\begin{verbatim}
  subroutine q2help(prompt,line)
! This routine is called from submenu
! when the user types a ?
    implicit none
    character*(*) prompt,line
!\end{verbatim} %+
    character helpquest*32
    integer savedlevel,kk,ip
!
    savedlevel=helprec%level-1
! If the ? is followed by a text push that text on the helprec%cpath
    ip=2
! This is to force q2help to work ... otherwise segmentation fault ??!!!    
    if(ip.lt.0) write(*,*)'q2help: ',savedlevel,line(1:20)
    if(.not.eolch(line,ip)) then
!       write(*,*)'q2help: ',helprec%level,ip,helprec%cpath(helprec%level)
       helpquest=line(ip:)
       helpquest=prompt
       call capson(helpquest)
! remove any WHAT? as such levels will be ignored by q1help
       kk=index(helpquest,'WHAT?')
       if(kk.gt.0) then
          helpquest(kk:)='COMMAND '
          if(helptrace) write(*,*)'MM hepquest: ',kk,trim(helpquest)
       endif
! use the saved helprec%level 
       helprec%level=savedlevel
       helprec%cpath(helprec%level)=helpquest
! always upper case ...
       call capson(helprec%cpath(helprec%level))
!       if(helptrace) write(*,11)helprec%level,&
!            (trim(helprec%cpath(i)),i=1,helprec%level)
!11     format('q2help: ',i3,10(', ',a))
    else
! when we are here we have just a ? from user, return to submenu with that
! with two ?? or anything else q1help is called (I hope ...)
       line='?!'
       if(ochelp%htmlhelp) then
          write(*,17)
17        format(/'By typing two ?? you will open the browser')
       endif
       goto 1000
    endif
! this is a dummy line needed to force the MacOS linker to find this routine
!??  if(savedlevel.eq.helprec%level) write(*,*)'Inside q2help: ',trim(prompt)
    if(ip.lt.0) write(*,*)'in q2help calling q1help'
! write help text from help file and then return with ?! to get submenu
    if(helptrace) write(*,*)'q2help calling q1help: ',trim(helpquest)
    call q1help(prompt,line)
    line='?!'
1000 continue
    return
  end subroutine q2help

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine q3help & Old help routine 3
!\begin{verbatim}
  SUBROUTINE Q3HELP(LINE,LAST,COMM,NC)
! used in submeny when user gives "? 'command' " taken as "help 'command'"
!...EXECUTES A HELP COMMAND
    implicit none
    CHARACTER LINE*(*),COMM(NC)*(*)
    integer last
!\end{verbatim} %+
    CHARACTER CMD*40
    integer, parameter :: MC=100
    integer INDX(MC)
    integer nkpl,nc,nlfk,nbk,i,j,k
! To avoid storing "COMMAND" in the helprec%cpath
!    if(helprec%level.gt.2) helprec%level=helprec%level-1 !.. HELP HELP not OK
    if(helprec%level.gt.3) helprec%level=helprec%level-1
!    write(*,*)'q3help: asking for help for command: "',trim(line),'"',last,nc
    CALL GPARC_old('Help for which command? ',LINE,LAST,1,CMD,'*',tophlp)
!    write(*,*)'q3help: command: "',trim(cmd),'"'
    IF(CMD(1:1).EQ.'*' .or. cmd(1:1).eq.'?') THEN
!...LIST ALL COMMANDS IN UNIX ALPHABETICAL ORDER
       NKPL=80/(LEN(COMM(1))+1)
       IF(NKPL*(LEN(COMM(1))+1).GE.80) NKPL=NKPL-1
       IF(NC.LT.MC) THEN
          CALL SSORT(COMM,NC,INDX)
          ALLCOM: DO NBK=1,NC
             IF(COMM(INDX(NBK))(1:1).NE.' ') GOTO 301
          enddo ALLCOM
301       NLFK=(NC+NKPL-NBK)/NKPL
          NBK=NBK-1
          COMLIST: DO I=1,NLFK
             WRITE(KOU,320)(COMM(INDX(NBK+J)),J=I,NC-NBK,NLFK)
          enddo COMLIST
320       FORMAT(10(1X,A))
       ELSE
!...      TOO MANY COMMANDS TO SORT
          NLFK=(NC+NKPL-1)/NKPL
          UNSORTED: DO I=1,NLFK
             WRITE(KOU,320)(COMM(J),J=I,NC,NLFK)
          enddo UNSORTED
       ENDIF
    ELSE
!...HELP <COMMAND>
!      IF UNIQUE LIST DESCRIPTION ON HELP FILE. OTHERWISE
!      ALL COMMANDS THAT MATCHES
       K=NCOMP2(CMD,COMM,NC,I)
       IF(K.GT.0) THEN
! we have to replace HELP by CMD on the stack of commands
! to get the correct help text
          CALL CAPSON(CMD)
          helprec%level=helprec%level-1
          helprec%cpath(helprec%level)=CMD(1:32)
!          write(*,11)helprec%level,(helprec%cpath(i)(1:8),i=1,helprec%level)
!11        format('q3help: ',i3,10(', ',a))
!          write(*,*)helprec%level
!          do ii=1,helprec%level
!             write(*,*)helprec%cpath(ii)
!          enddo
          call q1help(' ',CMD)
       ELSEIF(K.EQ.0 .OR. K.LT.-NC) THEN
          WRITE(KOU,*)'No matching command, use HELP * or ?'
       ELSE
500       WRITE(KOU,*)COMM(-K)
          IF(NC+K.LE.0) GOTO 900
          J=NCOMP2(CMD,COMM(1-K),NC+K,I)
!  ...bugfix for "help s-i" in poly
          IF(j .LT. -(NC+K) ) GOTO 900
          IF(K.EQ.-NC .OR. J.EQ.0) GOTO 900
          K=K-ABS(J)
          GOTO 500
       ENDIF
    ENDIF
900 RETURN
  END SUBROUTINE Q3HELP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine q3helpx & New help routine 3
!\begin{verbatim}
  SUBROUTINE Q3HELPx(LINE,LAST,COMM,NC)
! used in submeny when user gives "? 'command' " taken as "help 'command'"
!...EXECUTES A HELP COMMAND
    implicit none
    CHARACTER LINE*(*),COMM(NC)*(*)
    integer last
!\end{verbatim} %+
    CHARACTER CMD*40
    integer, parameter :: MC=100
    integer INDX(MC)
    integer nkpl,nc,nlfk,nbk,i,j,k
! To avoid storing "COMMAND" in the helprec%cpath
!    if(helprec%level.gt.2) helprec%level=helprec%level-1 !.. HELP HELP not OK
    if(helprec%level.gt.3) helprec%level=helprec%level-1
!    write(*,*)'q3help: asking for help for command: "',trim(line),'"',last,nc
    CALL GPARCx('Help for which command? ',LINE,LAST,1,CMD,'*','?TOPHLP')
!    write(*,*)'q3help: command: "',trim(cmd),'"'
    IF(CMD(1:1).EQ.'*' .or. cmd(1:1).eq.'?') THEN
!...LIST ALL COMMANDS IN UNIX ALPHABETICAL ORDER
       NKPL=80/(LEN(COMM(1))+1)
       IF(NKPL*(LEN(COMM(1))+1).GE.80) NKPL=NKPL-1
       IF(NC.LT.MC) THEN
          CALL SSORT(COMM,NC,INDX)
          ALLCOM: DO NBK=1,NC
             IF(COMM(INDX(NBK))(1:1).NE.' ') GOTO 301
          enddo ALLCOM
301       NLFK=(NC+NKPL-NBK)/NKPL
          NBK=NBK-1
          COMLIST: DO I=1,NLFK
             WRITE(KOU,320)(COMM(INDX(NBK+J)),J=I,NC-NBK,NLFK)
          enddo COMLIST
320       FORMAT(10(1X,A))
       ELSE
!...      TOO MANY COMMANDS TO SORT
          NLFK=(NC+NKPL-1)/NKPL
          UNSORTED: DO I=1,NLFK
             WRITE(KOU,320)(COMM(J),J=I,NC,NLFK)
          enddo UNSORTED
       ENDIF
    ELSE
!...HELP <COMMAND>
!      IF UNIQUE LIST DESCRIPTION ON HELP FILE. OTHERWISE
!      ALL COMMANDS THAT MATCHES
       K=NCOMP2(CMD,COMM,NC,I)
       IF(K.GT.0) THEN
! we have to replace HELP by CMD on the stack of commands
! to get the correct help text
          CMD=COMM(K)
!          CALL CAPSON(CMD)
!          helprec%level=helprec%level-1
!          helprec%cpath(helprec%level)=CMD
!          write(*,11)helprec%level,(helprec%cpath(i)(1:8),i=1,helprec%level)
!11        format('q3help: ',i3,10(', ',a))
!          write(*,*)helprec%level
!          do ii=1,helprec%level
!             write(*,*)helprec%cpath(ii)
!          enddo
          write(*,*)'Calling q4help from q3helpx: ',trim(cmd)
          call q4help(cmd,0)
!          call q1help(' ',CMD)
       ELSEIF(K.EQ.0 .OR. K.LT.-NC) THEN
          WRITE(KOU,*)'No matching command, use HELP * or ?'
       ELSE
500       WRITE(KOU,*)COMM(-K)
          IF(NC+K.LE.0) GOTO 900
          J=NCOMP2(CMD,COMM(1-K),NC+K,I)
!  ...bugfix for "help s-i" in poly
          IF(j .LT. -(NC+K) ) GOTO 900
          IF(K.EQ.-NC .OR. J.EQ.0) GOTO 900
          K=K-ABS(J)
          GOTO 500
       ENDIF
    ENDIF
900 RETURN
  END SUBROUTINE Q3HELPx

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine q4help & New help routine 4
!\begin{verbatim}
  subroutine q4help(hypertarget,extra)
! This routine is adapted to provide help from webrowsers using hypertarget
! when the user types a ? or ??
    implicit none
    integer extra
    character*(*) hypertarget
!\end{verbatim} %+
! this routine independent of the command history.  The GPARX routines
! call this with the hypertarget provided in the call to GPARX routine
! and searches this target in the HTML file
! if extra=0 help has been provided, otherwise the calling routine try to do it
    character htmlhelp*256
    if(helprec%okinit.eq.0) then
       if(helptrace) write(kou,*)'Sorry no help file'
       goto 1000
    endif
    if(helprec%type.ne.'html    ') then
       write(*,*)'Sorry only help based on HTML implemented'
       goto 1000
    endif
!    if(helptrace) then
! helptrace help debugging ...
!       write(*,*)'q4help: ',trim(hypertarget),extra
!    endif
    if(hypertarget(1:1).eq.' ') then
       write(*,*)'Sorry, the software provides no help for this question'
       goto 1000
    endif
!       
! we have tested this file exists when initating help 
!    write(*,*)'Q4HELP: ',trim(hypertarget),extra
!    open(31,file=ochelp%latexfile,status='old',access='sequential')
! if first character in hypertarget is ? remove that
! in OC I try to use a ? in all calls for gparxyz to find hypertargets in 
! the source code.  Seach for "'?" in the source code!
! This ? is not needed or used in the LaTeX file.
    if(hypertarget(1:1).eq.'?') then
       ochelp%target=hypertarget(2:)
    else
       ochelp%target=hypertarget
    endif
! The help system depends on 3 files:
! 1: a plain LaTeX file is the base.  The \hypertarget{target}{text} feature
!    is are used to find a specific help text in the user guide.
! 2: a html file is generated from this LaTeX file with the hypertargets.
!    this allows to locate the text inside the html file and display
!    in a separate windows in a browser.  The user can scroll this 
!    while running the program.
! 3: the same LaTeX can also be used to generate a PDF.  But no one reads 
!    the manual.
! For each command and question the software asks it uses a GPARX routine.
! in the call this subroutine and a hypertarget text is provided.
! When inside this routine the user has typed ? or ?? to get help.
! The browser used depend on compiler options ...
#ifdef winhlp
! on Windows we need the / after file
! the initial start spawns a new window with the browser, each ? a new browser
    htmlhelp='start '//trim(ochelp%browser)//' "file:/'//&
         trim(ochelp%htmlfile)//'#'//ochelp%target//'"'
#else
! on linux or Mac just ' "file:" as ochelp#htmlfile start with a /
! The & at the end spawns the browser window and furter ? creates new tags !!
    htmlhelp=trim(ochelp%browser)//' "file:'//&
         trim(ochelp%htmlfile)//'#'//ochelp%target//'" &'
#endif
    if(helptrace) write(*,*)'QZ: ',trim(htmlhelp)
    call execute_command_line(htmlhelp)
    close(31)
!
1000 continue
    return
  end subroutine q4help

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine nohelp & No help
!\begin{verbatim}
  SUBROUTINE NOHELP(PROMT,LINE)
! no help available
    implicit none
    CHARACTER PROMT*(*),LINE*(*)
!\end{verbatim} %+
    RETURN
  END SUBROUTINE NOHELP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine tophlp & Help from top level
!\begin{verbatim}
  SUBROUTINE TOPHLP(PROMPT,LINE)
! return to calling routine for help, do not save the current command ...
    implicit none
    CHARACTER PROMPT*(*),LINE*(*)
!\end{verbatim} %+
!    helprec%level=helprec%level-1
!    write(*,11)helprec%level,(helprec%cpath(i)(1:8),i=1,helprec%level)
!11  format('tophlp: ',i3,10(', ',a))
    LINE(2:2)='!'
    RETURN
  END SUBROUTINE TOPHLP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable logical function yeschk & Check for Y or y
!\begin{verbatim}
  LOGICAL FUNCTION YESCHK(CH1)
!    returns TRUE if CH1 is Y or y
    CHARACTER CH1*1
!\end{verbatim}
    YESCHK=.FALSE.
    IF(CH1.EQ.'Y' .OR. CH1.EQ.'y') YESCHK=.TRUE.
    RETURN
  END FUNCTION YESCHK

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
! History of commands
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine nghist & Execute history caommand
!\begin{verbatim}
  SUBROUTINE NGHIST(LINE,LAST)
!...EXECUTES A HISTORY COMMAND
!       LAST IS SET TO 0 IF LINE IS SET TO A COMMAND FROM HISTORY LIST
!    CHARACTER HIST*80,LINE*(*),CH1*1
    implicit none
    CHARACTER LINE*(*)
    integer last
!\end{verbatim} %+
    CHARACTER CH1*1
    LOGICAL IED
! LHL, LHM,LHP are prive global variables
    integer LOW,KADD,K,IDIG
!
    CHARACTER*1 CHHIST,CHHELP,CHEDIT
    PARAMETER (CHHIST='!',CHHELP='?',CHEDIT='*')
    LAST=LAST+1
    CH1=LINE(LAST:LAST)
    IF(CH1.EQ.CHHELP) THEN
       WRITE(KOU,100)
100    FORMAT(' History commands are useful when the same command ', &
            ' shall'/' be executed several times. It is also possible to',&
            ' amend or correct'/' a command before it is executed again.'//&
            ' A history command always begins with an !'/&
            ' !?    gives this help'/&
            ' !!    gives a list of the history'/&
            ' !<digit>      executes the command <digit> in the history',&
            ' list'/&
            ' !<text>       executes the most recent command starting with',&
            ' <text>'/&
            ' !*<digit> or !*<text> makes the command available for',&
            ' editing before execution'//&
            ' NOTE the < > around digit and text should not be typed!')
    ELSEIF(CH1.EQ.CHHIST) THEN
!...       A LIST OF THE HISTORY
       IF(LHM.GT.0) THEN
          LOW=LHL
          KADD=-20
       ELSEIF(LHL.EQ.0) THEN
          WRITE(KOU,*)'No history yet!'
          GOTO 900
       ELSE
          LOW=0
          KADD=0
       ENDIF
!...       LOOP
200    LOW=LOW+1
       IF(LOW.GT.20) THEN
          LOW=1
          KADD=KADD+20
       ENDIF
       K=LEN_TRIM(HIST(LOW))
       IF(K.GT.0)WRITE(KOU,210)LHM+LOW+KADD,HIST(LOW)(1:K)
210    FORMAT(I5,'> ',A)
       IF(LOW.NE.LHL) GOTO 200
    ELSE
!...       A COMMAND TO BE EXECUTED OR EDITED SHOULD BE FOUND
       IF(CH1.EQ.CHEDIT) THEN
          LAST=LAST+1
          IED=.TRUE.
       ELSE
          IED=.FALSE.
       ENDIF
       CALL GETINT(LINE,LAST,IDIG)
       IF(BUPERR.NE.0) THEN
          BUPERR=0
          K=LEN_TRIM(LINE)
          LOW=LHL
          IF(K.GT.LAST) THEN
400          CONTINUE
             IF(HIST(LOW)(1:K-LAST+1).NE.LINE(LAST:K)) THEN
                LOW=LOW-1
                IF(LOW.EQ.0) LOW=20
                IF(LOW.NE.LHL) GOTO 400
                WRITE(KOU,*)'No matching command'
                GOTO 900
             ENDIF
          ENDIF
       ELSE
          IF(IDIG.GE.-20 .AND. IDIG.LT.0) IDIG=LHM+LHL+IDIG+1
          IF(IDIG.LE.MAX(0,LHM+LHL-20) .OR. IDIG.GT.LHM+LHL) THEN
             WRITE(KOU,*)'Number outside history'
             GOTO 900
          ENDIF
          LOW=MOD(IDIG,20)
          IF(LOW.EQ.0) LOW=20
       ENDIF
       LINE=HIST(LOW)
       em1: IF(IED) THEN
       ENDIF em1
       LAST=0
    ENDIF
900 RETURN
  END SUBROUTINE NGHIST

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine openlogfile & Opem log file
!\begin{verbatim}
  subroutine openlogfile(name,text,lun)
! opens a logfile for commands, if exits it will overwrite
    implicit none
    character name*(*),text*(*)
    integer lun
!\end{verbatim} %+
    integer kkp,ierr
    if(lun.le.0) then
       if(logfil.gt.0) close(logfil)
       goto 1000
    endif
!    write(*,*)'METLIB: opening logfile: "',trim(name),'"'
    if(len_trim(name).le.0) then
       name='OCLOG.LOG'
       write(*,*)'No logfile name, using default: ',trim(name)
    else
! it seems tinyfiledialogs return working directory ...
       kkp=index(name,'.')
       if(kkp.le.0) then
          kkp=len_trim(name)
          name(kkp+1:)='./OCLOG.LOG'
          write(*,*)'No logfile extention, using: ',trim(name)
       endif
    endif
    open(lun,file=name,access='sequential',status='unknown',&
         err=1100,iostat=ierr)
    write(lun,10)text(1:len_trim(text))
10  format('Logfile title: ',a)
    logfil=lun
1000 continue
    return
! error opening
1100 continue
    buperr=ierr
    goto 1000
  end subroutine openlogfile

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine set_echo & Set/reet echo of commands
!\begin{verbatim}
  subroutine set_echo(ion)
! set echo of command input, does this really work?
    implicit none
    integer ion
!\end{verbatim}
    jecho=ion
    return
  end subroutine set_echo

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!
! >>>> subsection 
!      output of promt for command
!      and input of command including command line editing on Linux
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
  
!\addtotable subroutine boutxt & Write a text noadvance
!\begin{verbatim}
  subroutine boutxt(lut,line)
! writes the text on line on unit lut without CR/LF
    implicit none
    integer lut
    character line*(*)
!\end{verbatim} %+
!    write(*,*)'boutxt; ',lut,line
    write(lut,10,advance='no')line
10  format(a)
    return
  end subroutine boutxt

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine bintxt & Read a text
!\begin{verbatim}
  subroutine bintxt(lin,cline)
! read a command line with or without arguments. On LINUX command line editing
    implicit none
    character cline*(*)
    integer lin
!\end{verbatim} %+
#ifdef lixed
! LINUX: to have command line editing uncomment the line above and comment the 
! line with the call bintxt_nogetkey
    call bintxt_getkey(lin,cline)
#else
! On Windows command line editing is provided by the OS
    call bintxt_nogetkey(lin,cline)
#endif
    return
  end subroutine bintxt

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine bintxt_getkey & Read a text with editing
!\begin{verbatim}
  subroutine bintxt_getkey(lin,cline)
! LINUX subroutine to read a line with history and editing a la emacs
!
    implicit none
    character cline*(*)
    integer lin
!\end{verbatim} %+
!--------------------  
! CONTROL CHARACTERS FROM KEYBOARD
! DEL delete curret character
    integer, parameter :: ctrla=1        ! CTRLA move cursor to first position
    integer, parameter :: backspace2=2   ! CTRLB move cursor one step left
    integer, parameter :: ctrlc=3        ! CTRLC terminate program
    integer, parameter :: ctrld=4        ! CTRLD delete char at cursor
    integer, parameter :: ctrle=5        ! CTRLE move cursor to last position
    integer, parameter :: forward=6      ! CTRLF move cursor one step right
    integer, parameter :: HELP=8         ! CTRLH give coordinates and update
    integer, parameter :: TAB=9          ! CTRLI end of input
    integer, parameter :: ctrlk=11       ! CTRLK delete to end of line
    integer, parameter :: return=13      ! CTRLM end of input
    integer, parameter :: DEL=127        ! DEL delete char left of cursor
    integer, parameter :: mode=17        ! CTRLQ toggle insert/replace
! on MAC same as UP DOWN FORWARD suck
    integer, parameter :: backspace=27   ! CTRL[ previous in history
!--------------------  
! UP previous history line (if any)
! DOWN and LF next history line (if any)
    integer, parameter :: CTRLP=16       ! CTRLP previous in history
!    integer, parameter :: UP=27         ! uparrow previous in history
    integer, parameter :: LF=10          ! CTRLJ next in history
    integer, parameter :: ctrln=14      ! CTRLN next in history
!--------------------
! backspace on a MAC screen
    integer, parameter :: tbackspace=8
!-----------
!
! ip is cursor position (>=1), lastp is last character on line (>=0)
    integer ip,lastp,kiud,jj,kou,size,hlast
    character line*128,ch1*1
!    character getkey
    logical endoftext
! global structure myhistory
!    type(chistory) :: myhistory
    logical insert
!
    kiud=5; kou=6
    size=1
!
!      write(*,*)'Reading input using getkey',lin,kiud
    if(lin.ne.kiud) then
! reading macro from file
       read(lin,10)cline
10     format(a)     
       goto 1000
    endif
!
! input trom terminal with editing
!
    insert=.TRUE.
    ip=1
    lastp=0
    line=' '
    endoftext=.true.
    hlast=myhistory%hpos+1
! read one character at a time without echo
100 continue
#ifdef lixed
! LINUX: read one character at a time without echo and allow editing history
    ch1=getkex()
#endif
110 continue
!    write(*,*)'got from getkey: ',ichar(ch1)
! handle control character
    if(ichar(ch1).ge.32 .and. ichar(ch1).lt.127) then
! printable character, write on screen and store inline     
       if(ip.eq.lastp+1 .or. .not.insert) then
          write(kou,10,advance='no')ch1
          line(ip:ip)=ch1
          if(ip.eq.lastp+1) lastp=lastp+1
          ip=ip+1
       else
! insert a character inside a text
          line(ip+1:)=line(ip:)
          line(ip:ip)=ch1
          lastp=lastp+1
!          write(kou,10,advance='no')tbackspace
          write(kou,10,advance='no')line(ip:lastp)
          do jj=ip,lastp-1
             write(kou,10,advance='no')tbackspace
          enddo
          ip=ip+1
       endif
       goto 100
    endif
!=======================  
!    write(*,*)'control character: ',ichar(ch1)
120 continue
    select case(ichar(ch1))
    case default
! ignore
!       write(*,*)'Ignoring ',ichar(ch1)
       goto 100
!............. OK
    case(ctrla)
! move cursor to first character
       do jj=1,ip-1
          write(kou,10,advance='no')tbackspace
       enddo
       ip=1
#ifdef lixed
!............. NEW handle arrow key on Linix/Mac
    case(backspace) ! try to handle arrow keys sequence of 27, 91, A/B/C/D
       ch1=getkex()
       if(ichar(ch1).ne.91) goto 110
       ch1=getkex()
       if(ch1.eq.'A') then
!          write(*,*)'Arrow up'
          ch1=char(ctrlp)
       elseif(ch1.eq.'B') then
!          write(*,*)'Arrow down'
          ch1=char(ctrln)
       elseif(ch1.eq.'C') then
!          write(*,*)'Arrow forward'
          ch1=char(forward)
       elseif(ch1.eq.'D') then
!          write(*,*)'Arrow backward'
          ch1=char(backspace2)
       else !page up/down which has similar sequences etc ignored
          goto 100
!          write(*,*)'Input messed up ...'
       endif
       goto 120
!...............OK
!    case(backspace,backspace2) ! ctrlb leftarrow (also up/down/right arrow)
    case(backspace2) ! ctrlb leftarrow (also up/down/right arrow)
! move cursor one step back
       if(ip.gt.1) then
          write(kou,10,advance='no')tbackspace
          ip=ip-1
       endif

#else
! this rooutine is never called on Windows
    case(backspace,backspace2) ! ctrlb leftarrow (also up/down/right arrow)
! move cursor one step back
       if(ip.gt.1) then
          write(kou,10,advance='no')tbackspace
          ip=ip-1
       endif
#endif
!............. OK
    case(ctrlc)
! terminate the program
       stop 'User break'
    case(ctrle)
! move cursor after last character
!       if(ip.eq.1) then
!          jj=ip+1
!       else
!          jj=ip
!       endif
       jj=ip
       do jj=jj,lastp
          write(kou,10,advance='no')line(jj:jj)
       enddo
       ip=lastp+1
!............. OK
    case(ctrld)
! delete character at cursor (ctrld)
       if(ip.gt.lastp .or. lastp.eq.0) goto 100
       jj=ip
! remove the character at position jj and write the whole line from jj to end
       line(jj:)=line(jj+1:)
       write(kou,10,advance='no')line(jj:lastp)
       do jj=lastp,ip,-1
          write(kou,10,advance='no')tbackspace
       enddo
       lastp=lastp-1
!............. OK
    case(del)
! delete character to the left of cursor (del), if ip=1 ignore
       if(ip.eq.1) goto 100
       write(kou,10,advance='no')tbackspace
! remove the character at position jj and write the whole line from jj to end
       jj=ip-1
       line(jj:)=line(jj+1:)
       write(kou,10,advance='no')line(jj:lastp)
       ip=ip-1
       lastp=lastp-1
! NOTE lastp can be zero here
! otherwise we should backspace lastp-ip positions
       do jj=lastp+1,ip,-1
          write(kou,10,advance='no')tbackspace
       enddo
!............. OK
    case(ctrlk)
! delete all characters from cursor to end of line
       if(ip.le.lastp) then
          line(ip:)=' '
!          write(kou,10,advance='no')tbackspace
          write(kou,10,advance='no')line(ip:lastp)
          do jj=ip,lastp
             write(kou,10,advance='no')tbackspace
          enddo
          lastp=ip-1
       endif
!.............
    case(help) ! ctrlh
       write(kou,77, advance='no')ip,lastp,line(1:lastp+1)
77     format(/'Current local values are: ',2i4/a)
!       write(kou,10,advance='no')'xyz'
       do jj=lastp,ip,-1
          write(kou,10,advance='no')tbackspace
       enddo
!.............
    case(mode)  ! crtlQ
! change inset mode
       if(insert) then
          insert=.FALSE.
       else
          insert=.TRUE.
       endif
!.............
    case(return,tab)
! save line (if not empty) finish editing and return current line
       cline=line
       if(len_trim(line).eq.0) then
          continue
!            write(*,*)'Not saving empty line'
       elseif(myhistory%hpos.le.0) then
! saving the first line as
          myhistory%hpos=1
          myhistory%hline(myhistory%hpos)=line(1:80)
       elseif(line(1:ip+1).eq.myhistory%hline(myhistory%hpos)(1:ip+1)) then
          continue
!          write(*,*)'Not saving same line'
       else
          if(myhistory%hpos.ge.histlines) then
! history full, the oldest history line deleted
             do jj=2,histlines
                myhistory%hline(jj-1)=myhistory%hline(jj)
             enddo
          else
             myhistory%hpos=myhistory%hpos+1
          endif
          myhistory%hline(myhistory%hpos)=line(1:80)
       endif
! write a CR on screen ... maybe also LF ?? NO!!
       if(ichar(ch1).eq.return) write(kou,*)
       goto 1000
!............. OK
    case(forward)
! move cursor one step right if we are not at lastp
       if(ip.le.lastp) then
          write(kou,10,advance='no')line(ip:ip)
          ip=ip+1
!          if(ip.eq.lastp) endoftext=.true.
!       elseif(.not.endoftext) then
!          write(kou,10,advance='no')line(ip:ip)
!          endoftext=.true.
       endif
!.............
    case(ctrlp)
! copy previous history line to current
! first remove anything on the line (not the question ...)
       if(hlast.gt.1) then
          do jj=1,ip-1
             write(kou,10,advance='no')tbackspace
          enddo
          line=' '
          write(kou,10,advance='no')line(1:lastp)
          do jj=1,lastp
             write(kou,10,advance='no')tbackspace
          enddo
          hlast=hlast-1
          line=myhistory%hline(hlast)
          lastp=len_trim(line)
          ip=lastp+1
          write(kou,10,advance='no')line(1:lastp)
       endif
!.............CTRLJ and CTRLN
    case(lf,ctrln)
! copy next history line to current
       if(hlast.lt.myhistory%hpos) then
          do jj=1,ip-1
             write(kou,10,advance='no')tbackspace
          enddo
          line=' '
          write(kou,10,advance='no')line(1:lastp)
          do jj=1,lastp
             write(kou,10,advance='no')tbackspace
          enddo
          hlast=hlast+1
          line=myhistory%hline(hlast)
          lastp=len_trim(line)
          ip=lastp+1
          write(kou,10,advance='no')line(1:lastp)
       endif
    end select
!-----------------
    goto 100
!=================  
1000 continue
    return
  end subroutine bintxt_getkey

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine bintxt_nogetkey & Read a text
!\begin{verbatim}
  subroutine bintxt_nogetkey(lin,line)
! Reading a command line on Windows with editing provided by the OS
    implicit none
    character line*(*)
    integer lin
!\end{verbatim}
    integer iostatus
    read(lin,10,iostat=iostatus)line
10  format(a)
    if(iostatus.lt.0) then
! reading a macro beyond EOL/EOF ??
       write(*,*)' *** WARNING: MACRO ENDS WITHOUT SET INTERACTIVE!'
       line='set inter '
    endif
    return
  end subroutine bintxt_nogetkey

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection 
!          command line macros
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine macbeg & Start a maro
!\begin{verbatim}
  SUBROUTINE MACBEG(LINE,LAST,OK)
!....subroutine to execute set-interactive allowing nesting of macros
!
! IDEA: addera lablar i macro sa man kan ange MACRO fil LABEL
! och vid stop som @? eller @& man kan interaktivt ange GOTO label
! Ocksa en generisk subrutin som gor att man kan fa fram ett variabelvarde
! call macsymval(package,symbol,ival,rval,cval)
!
    implicit none
    CHARACTER LINE*(*)
    LOGICAL OK
    integer last
!\end{verbatim} %+
!
!    CHARACTER MACFIL*256,FIL*256,CH1*1
    CHARACTER MACFIL*256,FIL*256
    LOGICAL FIRST
    character*3 USEEXT
    character*1 dirsep,backslash
    character dummy*256
    integer ll,kk,ierr
    SAVE FIRST
    DATA FIRST/.TRUE./
    IF(FIRST) THEN
       FIRST=.FALSE.
       IUMACLEVL=0
!       lun=50
    ENDIF
    MACFIL=' '
    useext=macext
!    ipos=index(line(max(1,last):),'.')
!    if (ipos.gt.0) then
!       if (LEN_TRIM(line(ipos:)).gt.1) then
!          useext=' '
!       endif
!    endif
!    write(*,*)'In MACBEG: ',trim(line),last
! added that extension should be OCM by the argument "3"
    CALL GPARFILEx('Macro filename: ',LINE,LAST,1,FIL,MACFIL,3,'?MACRO file')
!    CALL GPARC('Macro filename: ',LINE,LAST,1,FIL,MACFIL,nohelp)
! add default extension if needed
    CALL FXDFLT(FIL,MACEXT)
    IF(BUPERR.NE.0) GOTO 910
!    write(*,*)'open macro: ',lun,IUMACLEVL
!    LUN=50
    OPEN(LUN,FILE=FIL,ACCESS='SEQUENTIAL',STATUS='OLD', &
         FORM='FORMATTED',IOSTAT=IERR,ERR=910)
! we can have macros nested 5 livels deep
    IF(IUMACLEVL.LT.5) THEN
       IUMACLEVL=IUMACLEVL+1
       MACROUNIT(IUMACLEVL)=KIU
    ELSE
!       CALL ST2ERR(1083,'MACBEG','TOO DEEPLY NESTED MACRO FILES')
       buperr=1083
       OK=.FALSE.
       GOTO 900
    ENDIF
! extract the PATH to this macro file, needed to open files inside the macro
    backslash=char(92)
!    write(*,*)'M3 macro file: ',trim(fil),' bacslash: ',backslash
    if(index(fil,backslash).gt.0) then
! this is on Windows
       dirsep=backslash
    else
! this is on UNIX type systems       
       dirsep='/'
    endif
!    write(*,*)'M3 macro file: ',trim(fil),' backslash: ',backslash,IUMACLEVL
    ll=1
    kk=0
    do while(ll.gt.0)
       kk=ll+kk
       dummy=fil(kk:)
       ll=index(dummy,dirsep)
!       ll=index(fil(kk:),dirsep)
    enddo
! we have found the position of the actual filename.  Save the path incl dirsep
    if(kk.gt.1) then
       macropath(IUMACLEVL)=fil(1:kk-1)
    else
       macropath(IUMACLEVL)=' '
    endif
!    write(*,*)'Macro path saved: ',IUMACLEVL,': ',trim(macropath(IUMACLEVL))
!    write(*,*)'Command input set: ',kiu,IUMACLEVL
! this is to suprees "press return to continue" but not implemented ...
    OK=.TRUE.
    KIU=LUN
    LUN=LUN+1
!    write(*,*)'Command input is: ',kiu
900 continue
    return
910 OK=.FALSE.
    write(*,*)'Error ',ierr,' opening macro file: ',trim(fil)
    buperr=1000+ierr
    GOTO 900
  end SUBROUTINE MACBEG

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine macend & End a macro
!\begin{verbatim}
    SUBROUTINE MACEND(LINE,LAST,OK)
! end of macro detected, close file and return to upper level
      implicit none
      CHARACTER LINE*(*)
      LOGICAL OK
      integer last
!\end{verbatim} %+
! set interactive gives back control to calling macro if any
    IF(KIU.NE.KIUD) THEN
       IF(KIU.NE.0) CLOSE(KIU)
!       write(*,*)'end of macro: ',kiu,kiud,IUMACLEVL
       IF(IUMACLEVL.GT.0) THEN
!          write(*,*)'calling macro: ',macrounit(IUMACLEVL)
          KIU=MACROUNIT(IUMACLEVL)
          IUMACLEVL=IUMACLEVL-1
       ELSE
!          write(*,*)'terminal: ',kiud
          KIU=KIUD
       ENDIF
    ENDIF
!    write(*,*)'Output: ',kou,koud
    IF(KOU.NE.KOUD) THEN
       IF(KOU.NE.0) CLOSE(KOU)
       KOU=KOUD
    ENDIF
!...ANYTHING AFTER A SET_INTERACTIVE IS TAKEN AS A MODULE NAME
    OK=.FALSE.
    IF(EOLCH(LINE,LAST)) GOTO 900
    OK=.TRUE.
    LAST=LAST-1
900 continue
!    write(*,*)'Leaving macbeg/macend'
    RETURN
  END SUBROUTINE MACEND

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gptcm1 & Replace macro variables with value 1
!\begin{verbatim}
  SUBROUTINE GPTCM1(IFLAG,SVAR,LAST,SLIN)
!...handling of MACRO directives like @& @? and @# etc
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER SVAR*(*),slin*(*)
    integer iflag,last
!\end{verbatim} %+
!    CHARACTER PP*30,CH1*1
!    CHARACTER ENVIR(9)*60,LABEL*8,LABLIN*60,SYMBOL*60
!    CHARACTER LABEL*8,LABLIN*60,SYMBOL*60
!    LOGICAL SG2ERR,TESTB,EOLCH
!
    IFLAG=0
!...IF NO ERROR RETURN
    IF(.NOT.BUPERR.NE.0) GOTO 900
    IF(LAST.LE.0 .OR. LAST.GE.LEN(SVAR)) GOTO 900
    SLIN=SVAR(LAST:)
!...IF FIRST CHARACTER NOT A @ RETURN
    call gptcm2(iflag,svar,last,slin)
900 continue
    return
  end SUBROUTINE GPTCM1
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gptcm2 & Replace a macro variable with value 2
!\begin{verbatim}
  subroutine GPTCM2(IFLAG,SVAR,LAST,SLIN)
! handling of macro variables
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER SVAR*(*),slin*(*)
    integer iflag,last
!\end{verbatim} %+
!    CHARACTER ENVIR(9)*60,LABEL*8,LABLIN*60,SYMBOL*60
    CHARACTER PP*30,CH1*1
!    CHARACTER LABEL*8,LABLIN*60,SYMBOL*60
    CHARACTER LABLIN*60
!    double precision, parameter :: ZERO=0.0D0,ONE=1.0D0
    integer ienv
!    LOGICAL SG2ERR,TESTB,EOLCH
    IFLAG=0
    IF(SLIN(1:1).NE.'@') GOTO 900
!...A @ MEANS THAT THE FOLLOWING CHARACTER MAY HAVE SPECIAL MEANING
! (many of these not implemented)
!       @$ MEANS A COMMENT LINE
!       @& MEANS PAUSE
!       @?text MEANS QUEARY, value supplied entered to program
!       @#itext MEANS DEFINING MACRO VARIABLE i, queary is "text", 0<i<10
!       @: MEANS LABEL
!       @( MEANS begin multiline comment to be terminated by a line with @)
!       @) MEANS end multiline comment
!       @@ MEANS an OS command
!       ##i MEANS replace with value of macro variable i
    IF(SLIN(2:2).EQ.'$') THEN
!...       A COMMENT. SKIP LINE, RESET ERROR CODE AND PROMPT AGAIN
       LAST=LEN(SVAR)
       BUPERR=0
       IFLAG=1
    ELSEIF(SLIN(2:2).EQ.'&') THEN
!...       A PAUSE REQUESTED. CONTINUE AFTER A RETURN, skipp if iox(8) nonzero
!       write(*,*)'calling testb from gptcm1'
       if(iox(8).eq.0) then
! if iox(8) nonzero do not stop
          WRITE(KOUD,*)'Hit RETURN to continue'
          READ(KIUD,310)CH1
310       FORMAT(A1)
       endif
       BUPERR=0
       IFLAG=1
    ELSEIF(SLIN(2:2).EQ.'?') THEN
!...       A QUEARY, PROMPT AND READ FROM TERMINAL
       PP=SLIN(3:)
       LAST=LEN_TRIM(PP)+1
!..........lh test
       if (kiu.ne.kiud) then
          CALL BOUTXT(KOU,SLIN)
          write(kou,*)
       else
          CALL BOUTXT(KOUD,PP(1:LAST))
       endif
!..........lh test
       CALL BINTXT(KIUD,SVAR)
       BUPERR=0
       LAST=0
       IFLAG=1
    ELSEIF(SLIN(2:2).EQ.'#') THEN
       IENV=ICHAR(SLIN(3:3))-ICHAR('0')
       IF(IENV.GT.9 .OR. IENV.LT.1) GOTO 900
       PP=SLIN(4:)
       LAST=LEN_TRIM(PP)+1
!       write(*,*)'extracting macro variable value 1: ',ienv,slin(1:20)
       CALL BOUTXT(KOUD,PP(1:LAST))
       CALL BINTXT(KIUD,SVAR)
       BUPERR=0
       ENVIR(IENV)=SVAR
       write(*,*)'Macro variable: ',ienv,' value: ',envir(ienv)(1:20)
! to avoid "no such command"
       svar='@$'
    ELSEIF(SLIN(2:2).EQ.':') THEN
       IF(KIU.EQ.KIUD) GOTO 900
       LABLIN=SVAR(3:)
       CALL CAPSON(LABLIN)
    ELSEIF(SLIN(2:2).EQ.'(') THEN
!...begin multiline comment, read until line starting with @)
100    CONTINUE
       CALL BINTXT(KIU,SVAR)
       IF(SVAR(1:2).NE.'@)') GOTO 100
       SVAR(2:2)='@'
    else
! system command
       write(*,*)'system command: ',slin(2:len_trim(slin))
       call execute_command_line(slin(2:))
!       call system(slin(2:))
    ENDIF
900 RETURN
!
!910 WRITE(LER,911)LABEL
!911 FORMAT('No label ',A,' terminating macro')
!    SVAR='set interactive '
!    GOTO 900
  END SUBROUTINE GPTCM2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine gqexv & Replace a macro variable with value
!\begin{verbatim}
  SUBROUTINE GQXENV(SVAR)
!...EXCHANGES REFERENCES TO ENVIRONMENT MACRO VARIABLES TO ACTUAL VALUES
!       REFERENCES ARE FOR EXAMPLE ##4
!    CHARACTER SVAR*(*),ENVIR(*)*(*),CH1*1,LABLIN*60,LABEL*8
    implicit none
    CHARACTER SVAR*(*)
!\end{verbatim}
!    COMMON/TCMACRO/IUL,IUN(5),MACEXT
!    CHARACTER*3 MACEXT
!    CHARACTER CH1*1,LABLIN*60,LABEL*8,HOLDER*200
    CHARACTER CH1*1,HOLDER*200
    integer k,ienv,l
    IF(SVAR(1:2).EQ.'@&') THEN
!       write(*,*)'calling testb from gqxenv'
!       IF(.NOT.TESTB(1,IOX(8))) THEN
          WRITE(KOUD,*)'Hit RETURN to continue'
          READ(KIUD,310)CH1
310       FORMAT(A1)
!       ENDIF
    ENDIF
100 K=INDEX(SVAR,'##')
    IF(K.GT.0) THEN
       IENV=ICHAR(SVAR(K+2:K+2))-ICHAR('0')
       IF(IENV.LT.1 .OR. IENV.GT.9) GOTO 900
       HOLDER=SVAR(K+3:)
!       kk=index(svar,' ')
       write(*,*)'I get: ',k,svar(1:k)
!       SVAR(K:)=ENVIR(IENV)
!       if(kk.gt.0) then
!          svar=svar(kk:k-1)//' '//envir(ienv)
!       else
          SVAR=ENVIR(IENV)
!       endif
       L=LEN_TRIM(SVAR)
       SVAR(L+1:)=HOLDER
       write(*,312)'replaced macro variable: ',k,l,svar(1:20)
312    format(a,2i3,' >',a)
       GOTO 100
    ENDIF
900 RETURN
  END SUBROUTINE GQXENV
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!
!
! >>>> subsection
! PUTFUN can parse a fortran like expression and create a binary tree
! It cannot calculate derivatives.  Used for state variable symbols in OC
! Rather final version of PUTFUN below     
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine putfun & Enter an expression
!\begin{verbatim}
  SUBROUTINE PUTFUN(STRING,L,MAXS,SYMBOL,LOKV,LROT,ALLOWCH,NV)
!...READS AN EXPRESSION FROM STRING POSITION L AND CREATES AN BINARY TREE
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER STRING*(*),SYMBOL(*)*(*)
    integer LOKV(*)
    integer maxs,allowch,nv
!    type(putfun_symlink) :: symlist
    LOGICAL NOTPP
    TYPE(putfun_node), pointer :: LROT
!\end{verbatim}
    CHARACTER CH*1
!    double precision, parameter :: ZERO=0.0D0
    integer l,i,ipn,lab,lq,negmark,iopuni
    double precision val
    TYPE(putfun_node), pointer :: nynod,nonod
! these dummy pointers seem to be redundant ... ???
    TYPE(putfun_node), pointer :: dummy,dummy2
! INPUT:
!      STRING CHARACTER WITH EXPRESSION
!      L     POSITION WHERE THE EXPRESSION STARTS
!      MAXS  NUMBER OF SYMBOLIC VARIABLES THAT MAY BE USED.
!            MAXS=0 NO VARIABLES ALLOWED
!            MAXS>0 INDICATES THAT MAXS ALLOWED NAMES ARE IN SYMBOL.
!            MAXS<0 INDICATES THAT ABS(MAXS) USERDEFINED NAMES ARE ALLOWED.
!            NOTE! SYMBOL AND LOKV MUST HAVE THE DIMENSION ABS(MAXS)
!      SYMBOL CHARACTER ARRAY OF DIMENSION ABS(MAXS)
!      LOKV  ARRAY FOR LOCAL USE OF DIMENSION ABS(MAXS)
! EXIT:
!      L     POSITION AFTER LAST CHARACTER IN THE EXPRESSION
!      LROT  POINTER TO ROOT (MAY BE ZERO IF NO TREE)
!      NV    NUMBER OF SYMBOLIC VARIABLES DEFINED BY USER
!...THE FOLLOWING NODES ARE USED
!      1.    OPKOD,LLINK,MLINK         BIN[R OPERATOR
!      2.    OPKOD,LLINK,0             UNIT[R OPERATOR
!      3.    OPKOD,0,0,VALUE           DATA ELLER VARIABEL
!      The following binary opcodes:
!      1     +
!      2     -
!      3     *
!      4     /
!      5     ** (EXPONENTIERING)
!      The following unary functions:
!      1     - (should not be used)
!      2     SQRT
!      3     EXP
!      4     LOG (Natural log, e-log)
!      5     LOG10 (logaritm bas 10)
!      6     SIN
!      7     COS
!      8     ATAN (Arctangent)
!      9     ERF (Error function)
!     10     IVAN (Ivantsof function)
!      The following data codes
!      -1    Whole number (stored as real in node)
!      0     Real constant stored in node
!      >0    Symbol variable, INT(value) is the index ????
!..INITIERA
    PUTFUNVAR=0
    NOTPP=.TRUE.
    debuginc=0
! IPN associated with the LEVEL array
    IPN=0
    NEGMARK=0
    LAB=1
    nullify(topnod); nullify(datanod); nullify(lastopnod); nullify(nonod)
!...initiate external symbolik links
!    write(*,*)'putfun: ',maxs
    do i=1,abs(maxs)
       LOKV(i)=0
    enddo
! I cannot make it work to have one node for all occurence of one symbol
!    allocate(symlist%symnod(20))
    nullify(stacktop)
    if(eolch(string,l)) then
! expected function, found nothing
       pfnerr=1059; goto 880
    endif
!..IF FIRST CHARACTER ; THEN EXIT with function zero
    IF(STRING(L:L).EQ.';') GOTO 800
    L=L-1
!    write(*,*)'PUTFUN: ',lab,l
    GOTO(100,100,200),LAB
!    if(lab.gt.0) goto 200
! ******* Expecting variable ********* > Expecting next
!      Allowed characters
!      ?                               > Give hekptext
!      - (make negative)               > Expect variable
!      Variabele or constant           > Binary operator
!      Unary operator (includ. '(')    > Expect variable
!      (                               > Expect variable
100 L=L+1
    if(pfnerr.ne.0) goto 900
    IF(L.GT.len(string)) then
       pfnerr=1052; GOTO 900
    endif
    CH=STRING(L:L)
    IF(CH.EQ.' ') GOTO 100
    IF(CH.EQ.'-') THEN
! treat negative sign special as one can have a symbol afterwards
       NEGMARK=-1
110    continue
       L=L+1
       IF(L.GT.len(string)) then
          pfnerr=1052; GOTO 900
       endif
! allow spaces after sign
       CH=STRING(L:L)
       IF(CH.EQ.' ') GOTO 110
    else
       NEGMARK=0
    endif
    IF(CH.EQ.'(') THEN
       CALL NYLP(nonod,IPN,NOTPP)
       GOTO 100
    endif
    LQ=L
    CALL GETREL(STRING,L,VAL)
    if(buperr.ne.0) then
       buperr=0
! not a number, it must be a symbol or unary operator
       L=LQ
!       write(*,*)'PUTFUN buperror: ',l
       CALL NYVAR(STRING,L,IOPUNI,negmark,MAXS,SYMBOL,LOKV,allowch,dummy2)
       IF(pfnerr.ne.0) GOTO 900
!       write(*,*)'After nyvar: ',iopuni,symbol(1)
       IF(IOPUNI.GT.0) THEN
          CALL NYUNI(IOPUNI,negmark,NYNOD,IPN,NOTPP)
          GOTO 100
       ENDIF
!       write(*,*)'PUTFUN buperror: nyvar return iopuni=0, look for operator'
    ELSE
! we have found a symbol
       CALL NYDAT(0,VAL,dummy,negmark)
       if(pfnerr.ne.0) goto 900
       L=L-1
    ENDIF
    LAB=3
! ****** Expecting binary operator **** > Expected next
!      Allowed characters              
!      +,-,*,**,/                       > Expect variable for right tree, LAB=2
!      )                                > Binary operator
!      ;                                > This means end of expression
200 L=L+1
    IF(L.GT.len(string)) GOTO 800
    CH=STRING(L:L)
    IF(CH.EQ.' ') GOTO 200
    IF(CH.EQ.';') GOTO 800
    binop: IF(CH.EQ.'+') THEN
       CALL NYBIN(1,NYNOD,NOTPP)
    ELSEIF(CH.EQ.'-') THEN
       CALL NYBIN(2,NYNOD,NOTPP)
    ELSEIF(CH.EQ.'*') THEN
       IF(STRING(L+1:L+1).EQ.'*') THEN
! exponentiation **
          L=L+1
          CALL NYBIN(5,NYNOD,NOTPP)
       ELSE
          CALL NYBIN(3,NYNOD,NOTPP)
       ENDIF
    ELSEIF(CH.EQ.'/') THEN
       CALL NYBIN(4,NYNOD,NOTPP)
    ELSEIF(CH.EQ.')') THEN
       CALL NYRP(IPN,NOTPP)
       GOTO 200
    ELSE
       write(*,*)'putfun error: "',ch,'" position ',L
       pfnerr=1051; GOTO 900
    ENDIF binop
    NEGMARK=0
    LAB=2
    GOTO 100
! we have evaluated the expression, IPN is parenthesis level
800 L=L+1
!    write(*,*)'PUTFUN: label 800'
    IF(IPN.NE.0) THEN
       pfnerr=1050; GOTO 900
    endif
! expression finished, set lrot
    if(associated(topnod)) then
       if(topnod%kod.ne.0) then
          lrot=>topnod
       else
! topnode has no binary operation, return datanod if any
          if(.not.associated(datanod)) then
             nullify(lrot)
          elseif(datanod%value.eq.zero) then
! single value equal to zero, do not return any node. A symbol would have 1.0
             nullify(lrot)
          else
             lrot=>datanod
          endif
       endif
    else
! there is no topnod
!       write(*,*)'PUTFUN: no topnode'
       if(associated(datanod)) then
          if(datanod%value.eq.zero) then
!             write(*,*)'PUTFUN: datanode with zero'
             nullify(lrot)
          else
!             write(*,*)'PUTFUN: datanode with non-zero'
             lrot=>datanod
          endif
!       else
! no topnode and no datanode, empty function
!          write(*,*)'PUTFUN: no datanode'
       endif
    endif
880 continue
! return number of external variables used
    NV=PUTFUNVAR
900 RETURN
  END SUBROUTINE PUTFUN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine nybin & Found a binary operator + - * /
!\begin{verbatim}
  SUBROUTINE NYBIN(kod,binnod,NOTPP)
!...INSERTS A NEW OPNODE IN THE TREE
    implicit none
    integer kod
    TYPE(putfun_node), pointer :: binnod
    LOGICAL NOTPP
!\end{verbatim} %+
!
    TYPE(putfun_node), pointer :: temp
!    double precision, parameter :: zero=0.0d0
! INPUT:
!      KOD    IS OPERATION CODE: 1 +, 2 -, 3 *, 4 /, 5 **
! EXIT: 
!      LOKBIN    IS NEW BINARY NODE WITH NEW OPERATION
!   1. IF THERE IS NO TOP NODE:
!         INSERT THIS NODE AS TOP NODE
!   2. IF THE NEW NODE IS + OR -:
!         THE PREVIOUS TOP NODE IS SET AS LEFT SUBTREE OF NEW NODE
!         THE NEW NODE IS SET AS TOP NODE
!   3. IF THE NEW NODE IS * OR /:
!         IF THE PREVIOUS TOP NODE IS * OR / OR ** DO AS 2.
!         ELSE THE RIGHT SUBTREE OF THE TOP NODE IS SET AS LEFT SUBTREE
!         OF NEW NODE. NEW NODE IS SET AS RIGHT SUBTREE OF THE TOP NODE
!   4. IF THE NEW NODE IS **:
!         THE RIGHT SUBTREE OF THE TOP NODE IS SET AS LEFT SUBTREE
!         OF NEW NODE. NEW NODE IS SET AS RIGHT SUBTREE OF THE TOP NODE
!
    IF(KOD.LE.0) then
       pfnerr=1058; goto 900
    endif
! one may get error "already allocated??"
    allocate(binnod)
    debuginc=debuginc+1
    binnod%debug=debuginc
    nullify(binnod%left); nullify(binnod%right); binnod%value=zero
    binnod%kod=kod; binnod%links=0
    lastopnod=>binnod
!...arrange binary opnodes according to priorities
    binop: IF(.not.associated(topnod)) THEN
! set this as topnod and link the datanod as left subtree
       topnod=>binnod
       topnod%left=>datanod
    ELSEIF(KOD.LE.2) THEN
!         + OR -
       binnod%left=>topnod
       topnod=>binnod
    ELSEIF(KOD.LE.4) THEN
!         * OR /  one has to consider priorities of operators
       if(topnod%kod.gt.2) then
          binnod%left=>topnod
          topnod=>binnod
       else
          binnod%left=>topnod%right
          topnod%right=>binnod
       endif
       NOTPP=.TRUE.
    ELSEIF(KOD.EQ.5 .AND. NOTPP) THEN
!         ** (TWO ** IN A ROW ILLEGAL)
       if(topnod%kod.gt.2) then
          binnod%left=topnod%right
          topnod%right=binnod
       else
! rearrange
          temp=>topnod%right
          if(associated(temp)) then
             binnod%left=>temp%right
             temp%right=>binnod
          else
             binnod%left=>topnod%right
             topnod%right=>binnod
          endif
       endif
       NOTPP=.FALSE.
    ELSE
       pfnerr=1058
    ENDIF binop
900 RETURN
  END SUBROUTINE NYBIN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine nyuni & Found a unary operator LOG EXP ...
!\begin{verbatim}
  SUBROUTINE NYUNI(KOD,negmark,uninod,IPN,NOTPP)
!   Creates a node with a unary operator
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    TYPE(putfun_node), pointer :: UNINOD
    LOGICAL NOTPP
    integer kod,negmark,ipn
!\end{verbatim} %+
!
    double precision, parameter :: one=1.0D0
    allocate(uninod)
    debuginc=debuginc+1
    uninod%debug=debuginc
    nullify(uninod%left); nullify(uninod%right); uninod%value=zero
    uninod%kod=-kod; uninod%links=0
!    write(*,*)'creating unary node ',kod,debuginc
    if(negmark.lt.0) then
       uninod%value=-one
    else
       uninod%value=one
    endif
! if there is a previous binary operator
    if(associated(lastopnod)) then
!       write(*,*)'linking unary node as right link: ',lastopnod%debug
       lastopnod%right=>uninod
    elseif(.not.associated(topnod)) then
       datanod=>uninod
    else
! this should nover happen
       pfnerr=1064
    endif
    CALL NYLP(uninod,IPN,NOTPP)
    RETURN
  END SUBROUTINE NYUNI

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine nylp & Found a left (
!\begin{verbatim}
  SUBROUTINE NYLP(uninod,IPN,NOTPP)
!...OPENING PARENTHESIS, push links on LEVEL. Also after unary operator
    implicit none
    TYPE(putfun_node), pointer :: uninod
    integer ipn
    LOGICAL NOTPP
!\end{verbatim} %+
    type(putfun_stack), pointer :: temp
!
    IPN=IPN+1
    IF(IPN.GT.20) THEN
       pfnerr=1055; goto 1000
    endif
    if(associated(stacktop)) then
       allocate(temp)
       temp%previous=>stacktop
       stacktop=>temp
    else
       allocate(stacktop)
    endif
    stacktop%savetop=>topnod
    stacktop%savebin=>lastopnod
    stacktop%saveuni=>uninod
!    NSTACK(1,IPN)=>topnod
!    NSTACK(2,IPN)=>lastopnod
! uninod is null if not ( after unary operator
!    NSTACK(3,IPN)=>uninod
! start new expression after (
    NOTPP=.TRUE.
    nullify(topnod)
    nullify(lastopnod)
1000 continue
    return
  end SUBROUTINE NYLP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine nyrp & Found a )
!\begin{verbatim}
  subroutine NYRP(IPN,NOTPP)
!...CLOSING PARENTHESIS
!    implicit double precision (a-h,o-z)
    implicit none
    integer ipn
    LOGICAL NOTPP
!\end{verbatim} %+
    TYPE(putfun_node), pointer :: uninod,subtree
!
    IF(IPN.LE.0) then
       pfnerr=1056; goto 1000
    endif
! save link to expression inside parenthesis
    if(.not.associated(topnod)) then
       subtree=>datanod
    else
       subtree=>topnod
    endif
! POP previous nstack
!    topnod=>NSTACK(1,IPN)
!    lastopnod=>NSTACK(2,IPN)
!    uninod=>nstack(3,IPN)
    topnod=>stacktop%savetop
    lastopnod=>stacktop%savebin
    uninod=>stacktop%saveuni
    stacktop=>stacktop%previous
    IPN=IPN-1
!    write(*,*)'right ): ',topnod%debug,lastopnod%debug,subtree%debug
! I do not understand why this IF is not related to those following
    IF(associated(uninod)) THEN
!       write(*,*)'right ) after unary function: ',uninod%debug
       uninod%left=>subtree
    endif
    if(associated(lastopnod)) then
       if(.not.associated(lastopnod%right)) then
          lastopnod%right=>subtree
       else
          datanod=>subtree
       endif
    elseif(associated(uninod)) then
       datanod=>uninod
    else
!...PARENTHESISED EXPRESSION IS LEFT SUBTREE OF EMPTY BINARY NODE.
       datanod=>subtree
    endif
    NOTPP=.TRUE.
1000 continue
    RETURN
  END SUBROUTINE NYRP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine nyvar & Found a symbol
!\begin{verbatim}
  SUBROUTINE NYVAR(TEXT,L,IOPUNI,negmark,MAXS,SYMBOL,LOKV,allowch,dummy2)
! inserts a symbol in an expression
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER TEXT*(*),SYMBOL(*)*(*)
    integer LOKV(*)
    integer iopuni,negmark,maxs,allowch
    type(putfun_node), pointer :: dummy2
!\end{verbatim} %+
!
    integer, parameter :: NOPER=14
    CHARACTER CH*1,NAME*16
!    double precision, parameter :: ZERO=0.0D0,ONE=1.0D0
    LOGICAL DEL2
    integer l,k,ln,i
!    type(putfun_symlink) :: symlist
    character*6, dimension(noper) :: OPER=&
        ['SQRT  ','EXP   ','LOG   ','LOG10 ','SIN   ','COS   ','ATAN  ',&
         'SIGN  ','ERF   ','IVAN  ','BSUM  ','ABS   ','HS    ','LN    ']
!
    IOPUNI=0
    DEL2=.FALSE.
70  NAME=' '
    K=1
100 CH=BIGLET(TEXT(L:L))
    L=L+1
    IF(K.GT.1 .AND. ((CH.GE.'0' .AND. CH.LE.'9') .or. CH.eq.'_')) THEN
       NAME(K:K)=CH
       K=K+1
    ELSE
! first letter must be A-Z
       IF(CH.GE.'A' .AND. CH.LE.'Z') THEN
          NAME(K:K)=CH
          K=K+1
       ELSEIF(K.GT.1 .AND. allowch.EQ.1) THEN
! allowch=1 means allow & and # in symbol names
          if(ch.eq.'#' .or. ch.eq.'&') then
             name(k:k)=ch
             k=k+1
          else
             goto 200
          endif
       ELSE
          GOTO 200
       ENDIF
    ENDIF
    GOTO 100
!..DEL2 is true if it is the second part of a symbol with a dot
200 IF(DEL2) GOTO 315
!.. if character is ( it must be a unary operator
    IF(CH.EQ.'(') GOTO 300
!..If character is "." it is an external derivative (like H.T)
    IF(CH.EQ.'.') GOTO 311
!
    IF(MAXS.LE.0 .AND. PUTFUNVAR.EQ.0) GOTO 220
!..compare with existing variable symbols
!      If maxs>0 no new symbols allowed, only those in symbol(1..maxs)
    IF(MAXS.GT.0) PUTFUNVAR=MAXS
! exact match required, include final space in name, max 15 characters
    ln=len_trim(name)+1
    if(ln.gt.16) then
       pfnerr=1062; goto 900
    endif
    DO I=1,PUTFUNVAR
       IF(NAME(1:ln).EQ.SYMBOL(I)(1:ln)) GOTO 230
    enddo
!..New symbol, error if MAXS>0, if not add it to SYMBOL and increment PUTFUNVAR
220 IF(MAXS.GT.0) GOTO 910
    IF(PUTFUNVAR.GE.ABS(MAXS)) GOTO 920
    PUTFUNVAR=PUTFUNVAR+1
    SYMBOL(PUTFUNVAR)=NAME
! if this is never set each occurence of the same symbol will have a node
!    LOKV(PUTFUNVAR)=1
! jump here from several places below
224 continue
    I=PUTFUNVAR
225 continue
    CALL NYDAT(I,one,dummy2,negmark)
!    write(*,*)'nyvar assigning symnod: ',dummy2%debug,dummy2%value
!    symlist%symnod(i)=>dummy2
    GOTO 800
!..Known symbol, with index I
!  If LOKV(I)=0 it is a predefined symbol without node and one must created
230 continue
    IF(LOKV(I).EQ.0) GOTO 225
!    CALL SETADS(symnod(I))
!    write(*,*)'nyvar 230: ',i
!    IF(associated(lastopnod)) THEN
!       lastopnod%right=>symlist%symnod(i)
!    ELSE
!       datanod=>symlist%symnod(i)
!    ENDIF
!...Keep track of number of links to this node
!    symlist%symnod(i)%links=symlist%symnod(i)%links-1
!    write(*,*)'nyvar 230: ',i,symlist%symnod(i)%links,symlist%symnod(i)%debug
!    return
!    GOTO 800
!=======================================
!..A unary OPERATOR
300 continue
!    write(*,*)'nyvar found unary operator: ',name
    DO I=1,NOPER
       IF(NAME(1:6).EQ.OPER(I)(1:6)) GOTO 330
    enddo
!..IF USERDEFINED UNARY OPERATORS (OR ARRAYS) ALLOWED (MAXS<0)
311 IF(MAXS.GT.0) GOTO 910
    IF(PUTFUNVAR.GE.ABS(MAXS)) GOTO 920
    PUTFUNVAR=PUTFUNVAR+1
    IF(CH.EQ.'.') THEN
       SYMBOL(PUTFUNVAR)=NAME(1:K-1)
       L=L-1
    ELSE
       SYMBOL(PUTFUNVAR)=NAME(1:K-1)//'('
       L=L-1
! FDMTP extracts text within parenthesis
       CALL FDMTP(TEXT,L,SYMBOL(PUTFUNVAR)(K+1:))
       K=LEN_trim(SYMBOL(PUTFUNVAR))
       SYMBOL(PUTFUNVAR)(K+1:K+1)=')'
       CH=TEXT(L:L)
!...   this line I do not understand
       IF(CH.NE.'.') L=L+1
    ENDIF
    CALL CAPSON(SYMBOL(PUTFUNVAR))
!...check if an external symbol with dot (derivative like H.T)
    IF(CH.NE.'.') GOTO 224
    L=L+1
    DEL2=.TRUE.
    GOTO 70
!
!...Second part of symbolic derivative after a "."
315 K=LEN_TRIM(SYMBOL(PUTFUNVAR))
    SYMBOL(PUTFUNVAR)(K+1:)='.'//NAME
    IF(CH.EQ.'(') THEN
       K=LEN_TRIM(SYMBOL(PUTFUNVAR))
       SYMBOL(PUTFUNVAR)(K+1:K+1)='('
       L=L-1
       CALL FDMTP(TEXT,L,SYMBOL(PUTFUNVAR)(K+2:))
       K=LEN_TRIM(SYMBOL(PUTFUNVAR))
       SYMBOL(PUTFUNVAR)(K+1:K+1)=')'
       L=L+1
    ENDIF
    CALL CAPSON(SYMBOL(PUTFUNVAR))
    GOTO 224
!
!...note unary opcode is one larger than opcode index!
330 continue
    IOPUNI=I+1
    L=L+1
800 L=L-2
900 RETURN
910 continue
    pfnerr=1053; GOTO 900
920 continue
    pfnerr=1054; GOTO 900
  END SUBROUTINE NYVAR

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine nydat & Found a numeric value
!\begin{verbatim}
  SUBROUTINE NYDAT(KOD,VAL,nynod,negmark)
! store a constant or symbol.  The address to the node is returned in lok
! which is used if the symbol is used several times.
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    integer kod,negmark
    TYPE(putfun_node), pointer :: nynod
    double precision val
!\end{verbatim}
!    write(*,*)'nydat 1: ',kod,negmark
    allocate(nynod)
    nynod%kod=kod; nynod%links=0
    nullify(nynod%left); nullify(nynod%right)
    debuginc=debuginc+1
    nynod%debug=debuginc
    if(negmark.lt.0) then
       nynod%value=-val
    else
       nynod%value=val
    endif
!    write(*,*)'nydat 2: ',nynod%kod,debuginc,nynod%value
    if(associated(lastopnod)) then
       if(.not.associated(lastopnod%right)) then
          lastopnod%right=>nynod
       else
! this should never happen
          write(*,*)'PUTFUN never never error 1'
          pfnerr=7777
       endif
!       write(*,*)'nydat 4A: ',lastopnod%kod,lastopnod%value
!       write(*,*)'nydat 4B: ',lastopnod%right%kod,lastopnod%right%value
!       write(*,*)'nydat 4C: ',lastopnod%left%kod,lastopnod%left%value
    else
       datanod=>nynod
!       write(*,*)'nydat 5: ',datanod%kod,datanod%value
    endif
    return
  end SUBROUTINE NYDAT

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable double precision evalf & Evaluate a function
!\begin{verbatim}
  double precision function evalf(LROT,VAR)
!      Calculates the value of an expression MEMORY LEAK 
! ?? I do not know what is the difference with evalf_x ??/BoS 190804
!
! VAR is array with values of symbols that can be referenced
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    double precision VAR(*)
    type(putfun_node), pointer :: lrot
!\end{verbatim} %+
!    character ch1*1
    double precision STACK(20)
!    type(putfun_node), pointer :: llink,current,mlink
    type(putfun_node), pointer :: current,mlink
    TYPE PUTFUN_SAVE
       integer right
       type(putfun_node), pointer :: savecurrent
       type(putfun_save), pointer :: previous
    end TYPE PUTFUN_SAVE
! these pointers are allocated creating memory leaks
    type(putfun_save), pointer :: topsave,temp
!    double precision, parameter :: ZERO=0.0D0
    integer last,lstp,kod
!
!...If LROT<=0 there is no expression, return sero
    STACK(1)=ZERO
    IF(.not.associated(LROT)) THEN
       GOTO 800
    ENDIF
!..INITIATE
    LAST=0
    LSTP=0
    current=>LROT
    nullify(topsave)
!    read(*,72)ch1
!72  format(a)
!71  format(a,5i5,1pe16.6)
!..New node, take is left link if any
100 continue
!    if(associated(current%right)) then
!       write(*,71)'>>>> evalf 100A1: ',current%debug,current%kod,&
!            current%left%debug,current%right%debug,current%links,current%value
!    elseif(associated(current%left)) then
!       write(*,71)'>>>> evalf 100A2: ',current%debug,current%kod,&
!            current%left%debug,0,current%links,current%value
!    else
!       write(*,71)'>>>> evalf 100A3: ',current%debug,current%kod,&
!            0,0,current%links,current%value
!    endif
! ERROR if I first set llink=>current%left and then tested llink if associated
    if(associated(current%left)) then
!       write(*,*)'Taking the left link and pushing current'
       LAST=LAST+1
       if(associated(topsave)) then
          allocate(temp)
          temp%previous=>topsave
          topsave=>temp
       else
          allocate(topsave)
          nullify(topsave%previous)
       endif
       topsave%savecurrent=>current
!       write(*,71)'evalf 100D: ',topsave%savecurrent%debug,current%left%debug
! mark that right link not visited
       topsave%right=1
       current=>current%left
    ELSE
!..If no left link the right link must be a data or unary negation
       KOD=current%kod
       LSTP=LSTP+1
       IF(KOD.GT.0) then
! unary operator, store the operation as a real
          STACK(LSTP)=VAR(KOD)
       else
          stack(lstp)=current%value
       endif
!       write(*,71)'evalf 100X: ',current%debug,kod,lstp,0,0,stack(lstp)
!..When coming here with LAST=0 the expression has been evaluated.
!  If not check if right link of current node has been visited
200    IF(LAST.LE.0) GOTO 800
       current=>topsave%savecurrent
!       write(*,71)'evalf 100YA: ',current%debug,topsave%right,current%kod
       IF(topsave%right.gt.0) THEN
!..Follow the right link
          if(associated(current%right)) then
             MLINK=>current%right
!             write(*,71)'evalf 100YB: ',mlink%debug,mlink%kod
!..Follow the left link of the right link but first mark that the right
! link of current has been visited
             topsave%right=-1
             current=>MLINK
!             write(*,71)'evalf 100Z: ',current%debug,current%kod,topsave%right
             GOTO 100
          ELSE
!..unary operator, in some cases it can have a sign
             CALL EUNARY(current%kod,STACK(LSTP))
             STACK(LSTP)=current%value*STACK(LSTP)
!             write(*,71)'evalf U: ',current%debug,current%kod,&
!                  lstp,0,0,stack(lstp)
          ENDIF
       ELSE
!..Binary operator with both left and right links evaluated
          LSTP=LSTP-1
!          write(*,73)'evalf B: ',current%debug,current%kod,lstp,&
!               stack(lstp),stack(lstp+1)
!73        format(a,3i3,2(1pe14.5))
          CALL EBINRY(current%kod,STACK(LSTP),STACK(LSTP+1))
       ENDIF
       LAST=LAST-1
       topsave=>topsave%previous
       IF(LAST.LT.0) goto 900
       IF(LAST.EQ.0) goto 800
       goto 200
    ENDIF
!    write(*,*)'evalf 799: ',current%debug,current%kod,lstp,current%value
    GOTO 100
!..KLAR
800 EVALF=STACK(1)
900 RETURN
  END FUNCTION EVALF

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable double precision function evalf_x & Evaluate a function
!\begin{verbatim}
  double precision FUNCTION EVALF_X(LROT,VAR)
!      Calculates the value of an expression
! ?? I do not know what is the difference with evalf ??/BoS 190804
!
! VAR is array with values of symbols that can be referenced
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    type(putfun_node), pointer :: lrot
    double precision VAR(*)
!\end{verbatim} %+
    double precision STACK(20)
!    character ch1*1
!    type(putfun_node), pointer :: llink,current,mlink
    type(putfun_node), pointer :: current,mlink
    TYPE PUTFUN_SAVE
       integer right
       type(putfun_node), pointer :: savecurrent
       type(putfun_save), pointer :: previous
    end TYPE PUTFUN_SAVE
! memory leak allocating pointers
!    type(putfun_save), target :: saverec
    type(putfun_save), pointer  :: topsave,temp
!    double precision, parameter :: ZERO=0.0D0
    integer last,lstp,kod
!
!...If LROT<=0 there is no expression, return sero
    IF(.not.associated(LROT)) THEN
       STACK(1)=ZERO
       GOTO 800
    ENDIF
!..INITIATE
    LAST=0
    LSTP=0
    current=>LROT
    nullify(topsave)
!    read(*,72)ch1
!72  format(a)
!71  format(a,5i5,1pe16.6)
!..New node, take is left link if any
100 continue
!    if(associated(current%right)) then
!       write(*,71)'>>>> evalf 100A1: ',current%debug,current%kod,&
!            current%left%debug,current%right%debug,current%links,current%value
!    elseif(associated(current%left)) then
!       write(*,71)'>>>> evalf 100A2: ',current%debug,current%kod,&
!            current%left%debug,0,current%links,current%value
!    else
!       write(*,71)'>>>> evalf 100A3: ',current%debug,current%kod,&
!            0,0,current%links,current%value
!    endif
! ERROR if I first set llink=>current%left and then tested llink if associated
    if(associated(current%left)) then
!       write(*,*)'Taking the left link and pushing current'
       LAST=LAST+1
       if(associated(topsave)) then
          allocate(temp)
          temp%previous=>topsave
          topsave=>temp
       else
          allocate(topsave)
          nullify(topsave%previous)
       endif
       topsave%savecurrent=>current
!       write(*,71)'evalf 100D: ',topsave%savecurrent%debug,current%left%debug
! mark that right link not visited
       topsave%right=1
       current=>current%left
    ELSE
!..If no left link the right link must be a data or unary negation
       KOD=current%kod
       LSTP=LSTP+1
       IF(KOD.GT.0) then
! unary operator, store the operation as a real
          STACK(LSTP)=VAR(KOD)
       else
          stack(lstp)=current%value
       endif
!       write(*,71)'evalf 100X: ',current%debug,kod,lstp,0,0,stack(lstp)
!..When coming here with LAST=0 the expression has been evaluated.
!  If not check if right link of current node has been visited
200    IF(LAST.LE.0) GOTO 800
       current=>topsave%savecurrent
!       write(*,71)'evalf 100YA: ',current%debug,topsave%right,current%kod
       IF(topsave%right.gt.0) THEN
!..Follow the right link
          if(associated(current%right)) then
             MLINK=>current%right
!             write(*,71)'evalf 100YB: ',mlink%debug,mlink%kod
!..Follow the left link of the right link but first mark that the right
! link of current has been visited
             topsave%right=-1
             current=>MLINK
!             write(*,71)'evalf 100Z: ',current%debug,current%kod,topsave%right
             GOTO 100
          ELSE
!..unary operator, in some cases it can have a sign
             CALL EUNARY(current%kod,STACK(LSTP))
             STACK(LSTP)=current%value*STACK(LSTP)
!             write(*,71)'evalf U: ',current%debug,current%kod,&
!                  lstp,0,0,stack(lstp)
          ENDIF
       ELSE
!..Binary operator with both left and right links evaluated
          LSTP=LSTP-1
!          write(*,73)'evalf B: ',current%debug,current%kod,lstp,&
!               stack(lstp),stack(lstp+1)
!73        format(a,3i3,2(1pe14.5))
          CALL EBINRY(current%kod,STACK(LSTP),STACK(LSTP+1))
       ENDIF
       LAST=LAST-1
       topsave=>topsave%previous
       IF(LAST.LT.0) goto 900
       IF(LAST.EQ.0) goto 800
       goto 200
    ENDIF
!    write(*,*)'evalf 799: ',current%debug,current%kod,lstp,current%value
    GOTO 100
!..KLAR
800 EVALF_X=STACK(1)
900 RETURN
  END FUNCTION EVALF_X

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine eunary & Evaluate a unary function
!\begin{verbatim}
  SUBROUTINE EUNARY(KOD,X)
! calculates a unary function such as LOG, EXP etc
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    integer kod
    double precision X
!\end{verbatim} %+
    double precision, parameter :: ONE=1.0D0
!    y=x
    IF(KOD.EQ.-1) X=-X
    IF(KOD.EQ.-2) X=SQRT(X)
    IF(KOD.EQ.-3) X=EXP(X)
    IF(KOD.EQ.-4) X=LOG(X)
    IF(KOD.EQ.-5) X=LOG10(X)
    IF(KOD.EQ.-6) X=SIN(X)
    IF(KOD.EQ.-7) X=COS(X)
    IF(KOD.EQ.-8) X=ATAN(X)
    IF(KOD.EQ.-9) X=SIGN(ONE,X)
    IF(KOD.EQ.-10) X=ERF(X)
    IF(KOD.EQ.-11) X=AIVAN(X)
    IF(KOD.EQ.-12) X=PF_BSUM(X)
    IF(KOD.EQ.-13) X=ABS(X)
    IF(KOD.EQ.-14) X=PF_HS(X)
! this is LN same as LOG, 10th log is LOG10
    IF(KOD.EQ.-15) X=LOG(X)
!    write(*,*)'eunary: ',y,x
    RETURN
  END SUBROUTINE EUNARY

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine ebinary & Evaluate a binary operator
!\begin{verbatim}
  SUBROUTINE EBINRY(KOD,X,Y)
! Calculates the value of a binary node with two data nodes
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    integer kod
    double precision X,Y
!\end{verbatim} %+
    integer nn
!
    IF(KOD.EQ.1) X=X+Y
    IF(KOD.EQ.2) X=X-Y
    IF(KOD.EQ.3) X=X*Y
    IF(KOD.EQ.4) THEN
       IF(Y.ne.zero) then
          X=X/Y
       else
          pfnerr=1063
       endif
    endif
    IF(KOD.EQ.5) THEN
       NN=INT(Y)
       IF(ABS(X).LE.0.1D-30) THEN
          X=0.0D0
       ELSEIF(ABS(DBLE(NN)-Y).LT.1.0D-30) THEN
          X=X**NN
       ELSEIF(X.GT.0.1D-30) THEN
          X=EXP(Y*LOG(X))
       ELSE
          X=0.0D0
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE EBINRY

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable double precision function aivan & Evaluate Ivantsov's function
!\begin{verbatim}
  double precision FUNCTION AIVAN(PECN)
!      CALCULATES THE DIMENSIONLESS SUPERCOOLING OF DIFFUSION BY
!      IVANTSOV'S SOLUTION
!...added by Zikui and also an updated ERF
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    double precision PECN
!      APPROXIMATIVE FORMULA FOR ERROR FUNCTION GIVEN BY:
!      ABRAMOWITZ AND STEGUN: HANDBOOK OF MATHEMATICAL FUNCTIONS,
!      NATIONAL BUREAU OF STANDARDS, 9TH EDITION, 1970
!\end{verbatim} %+
    double precision, parameter :: ONE=1.0D0,TWO=2.0D0,PI=3.141592654D0
    double precision A,C,Q
    integer i
    IF(PECN.LE.8.5D0) THEN
       AIVAN=DSQRT(PI*PECN)*DEXP(PECN)*(ONE-ERF(DSQRT(PECN)))
    ELSE
       A=ONE
       C=ONE
       Q=ONE
       DO I=1,9
          A=A*(TWO*DBLE(I)-ONE)/TWO/PECN
          C=-C
          Q=Q+A*C
       enddo
       AIVAN=Q
    ENDIF
    RETURN
  END FUNCTION AIVAN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable double precision function pf_bsum & Evaluate BSUM
!\begin{verbatim}
  double precision FUNCTION PF_BSUM(FA)
!.. 1993-10-06 20:10:56 /BJ
!
!                                 ( sin(n*pi*f) )^2  
!  Calc. the infinit sum B = sum(-------------------)
!                                    (n*pi)^3        
!
!.. If we truncate the sum at N=200 the relative error is
!   less than ?% for 0.01 < F < 0.99
!
!
    implicit none
!    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    double precision FA
!\end{verbatim} %+
! value of PI not very accurate ....
!    double precision, parameter :: ZERO=0.0D+00,PI=3.14159D0
    double precision, parameter :: PI=3.14159D0
    double precision, parameter :: PI3=PI*PI*PI
!
    integer loopmx
    double precision val,a,b
    integer i
!
    LOOPMX=1000
    VAL=ZERO
    DO I=1,LOOPMX
       A=DBLE(I)
       B=SIN(A*PI*FA)
       VAL=VAL+B*B/(A*A*A*PI3)
    enddo
!
    PF_BSUM=VAL
    RETURN
  END FUNCTION PF_BSUM

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable double precision function pf_hs & Evaluate Heaviside 
!\begin{verbatim}
  double precision FUNCTION PF_HS(X)
!      Calculates Heaviside function
!    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    implicit none
    double precision X
!\end{verbatim} %+
!    double precision, parameter :: ZERO=0.0D+00,ONE=1.0D+00
! BUG!!!!
!    HS=ZERO
    PF_HS=ZERO
    IF (X.GE.ZERO) PF_HS=ONE
    RETURN
  END FUNCTION PF_HS

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable double precision function pf_erf & Evaluate ERF
!\begin{verbatim}
  double precision FUNCTION PF_ERF(X0)
!      CALCULATES ERROR-FUNCTION OF X, USING AN
!      APPROXIMATIVE FORMULA GIVEN BY:
!      ABRAMOWITZ AND STEGUN: HANDBOOK OF MATHEMATICAL FUNCTIONS,
!      NATIONAL BUREAU OF STANDARDS, 9TH EDITION, 1970
    implicit none
    double precision X0
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!\end{verbatim}
    double precision, parameter :: ONE=1.0D0,TWO=2.0D0
    double precision P,A1,A2,A3,A4,A5,PI,S,X,T,Q
    DATA P,A1,A2,A3,A4,A5,PI/.3275911D0,.254829592D0,-.284496736D0, &
         1.421413741D0,-1.453152027D0,1.061405429D0,3.141592654D0/
    S=DSIGN(ONE,X0)
    X=DABS(X0)
    T=ONE/(ONE+P*X)
    Q=T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))
    Q=(ONE-Q*DEXP(-X*X))*S
    PF_ERF=Q
    RETURN
  END FUNCTION PF_ERF

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine wrtfun & Write the function
!\begin{verbatim}
  SUBROUTINE WRTFUN(STRING,IPOS,LROT,SYMBOL)
!      Writes a PUTFUN expression
!
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER STRING*(*),SYMBOL(*)*(*)
    integer ipos
    type(putfun_node), pointer :: lrot,current,lnode,rnode,tnode
!\end{verbatim} %+
    TYPE PUTFUN_SAVE
       integer right
       type(putfun_node), pointer :: savecurrent
       type(putfun_save), pointer :: previous
    end TYPE PUTFUN_SAVE
    type(putfun_save), pointer :: topsave,temp,bug
    integer last,noneg,negmark
    double precision val
!    type(putfun_node), dimension(20) :: link
!    integer doneboth(20)
!    string=' '
!    ipos=1
!   write(*,*)'wrtfun: ',trim(symbol(1))
!...Quick return if no expression
    IF(.not.associated(LROT)) THEN
       CALL CONS(STRING,IPOS,'0')
       GOTO 900
    ENDIF
!..INITIATE
    LAST=0
    current=>LROT
!    write(*,*)'wrtfun 1: ',current%kod
    nullify(topsave)
    IF(current%kod.LT.0) then
!..start with a unary operator
       noneg=1
       if(current%value.lt.zero) noneg=-1
       CALL WRTLPQ(STRING,IPOS,0,0,current%kod,noneg)
!       write(*,76)'wrtfun 99: ',current%kod,current%debug,current%value
!       write(*,77)'wrtfun 100A: ',ipos,string(1:ipos)
    endif
!..new node, visit its left link
100 continue
    lnode=>current%left
!    write(*,76)'wrtfun 99: ',current%debug,current%kod,current%value
    if(associated(topsave)) then
    endif
    bigif: IF(associated(LNODE)) THEN
!..PUSH LINK
       IF(associated(LNODE%left)) &
            CALL WRTLPQ(STRING,IPOS,1,current%kod,lnode%kod,1)
       if(associated(topsave)) then
          allocate(temp)
          temp%previous=>topsave
          topsave=>temp
       else
          allocate(topsave)
          nullify(topsave%previous)
       endif
       topsave%savecurrent=>current
! write all saved links
       bug=>topsave
55     continue
       if(associated(bug%previous)) then
          bug=>bug%previous
          goto 55
       endif
!
! mark that right link not visited
       LAST=LAST+1
       topsave%right=1
       current=>LNODE
!       write(*,77)'wrtfun 100A: ',ipos,string(1:ipos)
!77     format(a,i3,' "',a,'"')
    ELSE !bigif
!..If no left link node must be data
       val=current%value
       negmark=0
       if(last.gt.0) then
! surround a negative value by ( ) if this data in a right link
          if(topsave%right.lt.0) negmark=-1
       endif
       CALL WRTDAQ(STRING,IPOS,current%kod,VAL,SYMBOL,negmark)
!       write(*,77)'wrtfun 100B: ',ipos,string(1:ipos)
       IF(LAST.EQ.0) GOTO 800
!..check if right link has been visited
200    continue
       current=>topsave%savecurrent
       smallif: IF(topsave%right.gt.0) THEN
!..follow right link
          RNODE=>current%right
          hlink: IF(associated(RNODE)) THEN
!..there is a right link, follow its left link
!      mark first the the right link of current has been visited
             topsave%right=-1
!      check if ) is needed
!      then write the operator and possibly a (
             TNODE=>current%left
             IF(associated(tnode%left)) &
                  CALL WRTRPQ(STRING,IPOS,1,current%kod,tnode%kod)
             CALL WRTBIQ(STRING,IPOS,current%kod)
!             write(*,77)'wrtfun 200A: ',ipos,string(1:ipos)
             TNODE=>current%right
             IF(associated(tnode%left)) then
                noneg=1
                if(current%kod.lt.0 .and. &
                     current%value.lt.zero) noneg=-1
                CALL WRTLPQ(STRING,IPOS,2,current%kod,tnode%kod,noneg)
!                write(*,77)'wrtfun 200B: ',ipos,string(1:ipos)
             endif
             current=>RNODE
             GOTO 100
          ELSE
!..unary operator, write ) if necessary
             IF(current%kod.LT.-1) CALL CONS(STRING,IPOS,')')
!             write(*,77)'wrtfun 200C: ',ipos,string(1:ipos)
          ENDIF hlink
       else !smallif
!..binary operator and both links visited, check if ) needed
! IT WAS A DIFFICULT BUG TO FIND WHEN tnode=topsave%savecurrent .....
          tnode=>topsave%savecurrent
          if(associated(tnode%right%left)) &
               call WRTRPQ(STRING,IPOS,2,tnode%kod,tnode%right%kod)
!          write(*,77)'wrtfun 200D: ',ipos,string(1:ipos)
       ENDIF smallif
       LAST=LAST-1
       topsave=>topsave%previous
!       IF(LAST) 900,800,200
!       write(*,*)'wrtfun 798: ',last
       IF(LAST.lt.0) goto 900
       IF(LAST.eq.0) goto 800
       goto 200
    ENDIF bigif
!    write(*,*)'wrtfun 799: ',current%kod,last
    GOTO 100
!..KLAR
800 continue
    string(ipos:ipos)=';'
    ipos=ipos+1
900 RETURN
  END SUBROUTINE WRTFUN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine wrtlpq & Write a keft (
!\begin{verbatim}
  SUBROUTINE WRTLPQ(STRING,IPOS,LINK,KOD,LOD,negmark)
! write a left ( or unary operator followed by (
! the unary operator is in LOD
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER STRING*(*)
    integer ipos,link,kod,lod,negmark
!\end{verbatim} %+
!    IF(LOD) 10,90,20
    IF(LOD.eq.0) goto 90
    if(LOD.gt.0) goto 20
!..unary operator
    IF(LOD.EQ.-1) THEN
       IF(KOD.LE.0) THEN
          CALL CONS(STRING,IPOS,'-')
       ELSE
          IF(LINK.EQ.1) CALL CONS(STRING,IPOS,'-(')
          IF(LINK.EQ.2) CALL CONS(STRING,IPOS,'(-')
       ENDIF
    ELSE
! this is a sign for a unary function
       if(negmark.eq.-1) then
          CALL CONS(STRING,IPOS,'-')
       endif
       IF(LOD.EQ. -2) CALL CONS(STRING,IPOS,'SQRT(')
       IF(LOD.EQ. -3) CALL CONS(STRING,IPOS,'EXP(')
       IF(LOD.EQ. -4) CALL CONS(STRING,IPOS,'LOG(')
       IF(LOD.EQ. -5) CALL CONS(STRING,IPOS,'LOG10(')
       IF(LOD.EQ. -6) CALL CONS(STRING,IPOS,'SIN(')
       IF(LOD.EQ. -7) CALL CONS(STRING,IPOS,'COS(')
       IF(LOD.EQ. -8) CALL CONS(STRING,IPOS,'ATAN(')
       IF(LOD.EQ. -9) CALL CONS(STRING,IPOS,'SIGN(')
       IF(LOD.EQ.-10) CALL CONS(STRING,IPOS,'ERF(')
       IF(LOD.EQ.-11) CALL CONS(STRING,IPOS,'IVAN(')
       IF(LOD.EQ.-12) CALL CONS(STRING,IPOS,'BSUM(')
       IF(LOD.EQ.-13) CALL CONS(STRING,IPOS,'ABS(')
       IF(LOD.EQ.-14) CALL CONS(STRING,IPOS,'HS(')
       IF(LOD.EQ.-15) CALL CONS(STRING,IPOS,'LN(')
    ENDIF
    GOTO 90
!..one must check LOD if left (
20  continue
    IF(KOD.GE.3 .AND. LOD.LT.3) CALL CONS(STRING,IPOS,'(')
    IF(KOD.EQ.5 .AND. LOD.GE.3) CALL CONS(STRING,IPOS,'(')
!..if LINK=2, i.e. a right link, write ( of KOD is - or /
    IF(LINK.EQ.2) THEN
       IF(KOD.EQ.2 .AND. LOD.EQ.1) CALL CONS(STRING,IPOS,'(')
       IF(KOD.EQ.4 .AND. LOD.EQ.3) CALL CONS(STRING,IPOS,'(')
    ENDIF
90  RETURN
  END SUBROUTINE WRTLPQ

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine wrtrpq & Write a right )
!\begin{verbatim}
  SUBROUTINE WRTRPQ(STRING,IPOS,LINK,KOD,LOD)
!  write a right )  but if LOD<-1 do not write (
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    CHARACTER STRING*(*)
    integer ipos,link,kod,lod
!\end{verbatim} %+
!    IF(LOD+1) 90,10,20
    IF(LOD+1.lt.0) goto 90
    if(LOD+1.gt.0) goto 20
!..negation need ) if KOD>0
    IF(KOD.GT.0) CALL CONS(STRING,IPOS,')')
    GOTO 90
20  IF(KOD.GE.3 .AND. LOD.LT.3) CALL CONS(STRING,IPOS,')')
    IF(KOD.EQ.5 .AND. LOD.GE.3) CALL CONS(STRING,IPOS,')')
    IF(LINK.EQ.2) THEN
       IF(KOD.EQ.2 .AND. LOD.EQ.1) CALL CONS(STRING,IPOS,')')
       IF(KOD.EQ.4 .AND. LOD.EQ.3) CALL CONS(STRING,IPOS,')')
    ENDIF
90  RETURN
  END SUBROUTINE WRTRPQ

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine wrtbiq & Write a binary operator 
!\begin{verbatim}
  SUBROUTINE WRTBIQ(STRING,IPOS,KOD)
! write a binary operator
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER STRING*(*)
    integer ipos,kod
!\end{verbatim} %+
!      write a binary operator
!    write(*,*)'wrtbiq 1: ',ipos,kod
    IF(KOD.EQ.1) CALL CONS(STRING,IPOS,'+')
    IF(KOD.EQ.2) CALL CONS(STRING,IPOS,'-')
    IF(KOD.EQ.3) CALL CONS(STRING,IPOS,'*')
    IF(KOD.EQ.4) CALL CONS(STRING,IPOS,'/')
    IF(KOD.EQ.5) CALL CONS(STRING,IPOS,'**')
    RETURN
  END SUBROUTINE WRTBIQ

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine wrtdaq & Write a number
!\begin{verbatim}
  SUBROUTINE WRTDAQ(STRING,IPOS,KOD,VAL,SYMBOL,negmark)
!     write a number, if KOD<0 a whole number
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
!    CHARACTER NAME*8,SYMBOL(*)*(*)
    CHARACTER SYMBOL(*)*(*)
    CHARACTER STRING*(*)
    integer ipos,kod,negmark
    double precision val
!\end{verbatim} %+
!    double precision, PARAMETER :: ZERO=0.0D0
    IF(KOD.EQ.0) THEN
       IF(VAL.GT.ZERO .or. negmark.eq.0) then
          CALL WRINUM(STRING,IPOS,12,0,VAL)
       elseif(VAL.LT.ZERO) then
          CALL CONS(STRING,IPOS,'(')
          CALL WRINUM(STRING,IPOS,12,0,VAL)
          CALL CONS(STRING,IPOS,')')
       endif
    ELSE
!..a name of a variable, the name is in SYMBOL(KOD), skip trailing spaces
! if negated surround by ( )
!       CALL CONS(STRING,IPOS,SYMBOL(KOD))
!       write(*,*)'wrtdaq symbol: ',kod,trim(symbol(kod))
       if(val.lt.zero) then
          CALL CONS(STRING,IPOS,'(')
          CALL CONS(STRING,IPOS,SYMBOL(KOD))
          CALL CONS(STRING,IPOS,')')
       else
          CALL CONS(STRING,IPOS,SYMBOL(KOD))
       endif
    ENDIF
    RETURN
  END SUBROUTINE WRTDAQ

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine cons & Concatinate
!\begin{verbatim}
  SUBROUTINE CONS(STR1,IPOS,STR2)
! used in PUTFUN but should be replaced by //
    implicit none
    CHARACTER STR1*(*),STR2*(*)
    integer ipos
!\end{verbatim}
!      CONS. TWO STRINGS, RESULT IN PARAMETER STR1
!      IPOS = POSITION IN STR1 WHERE STR2 SHOULD BE PUT
!      IPOS IS UPPDATED TO THE FIRST FREE POSITION AT THE END
!      OF STR1. TRAILING SPACES ARE STRIPPED OFF.
!      IF STR2 CONTAINES ONLY SPACES ONE SPACE IS WRITTEN
!      IN TO STR1.
    CHARACTER SPC*(1)
    PARAMETER (SPC=' ')
    integer ilen,k,i
    ILEN=LEN(STR2)
!...FIND THE LENGHT OF STR2
    K=ILEN
    DO I=K,1,-1
       IF(STR2(ILEN:ILEN).EQ.SPC) ILEN=ILEN-1
    enddo
    IF(ILEN.EQ.0)ILEN=1
    STR1(IPOS:IPOS+ILEN-1)=STR2(1:ILEN)
    IPOS=IPOS+ILEN
    RETURN
  END SUBROUTINE CONS

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine exphlp & Provide help
!\begin{verbatim}
!  SUBROUTINE EXPHLP(PROMPT,SVAR)
  SUBROUTINE EXPHLP
! writes help to enter a PUTFUN expression
    implicit none
!    CHARACTER PROMPT*(*),SVAR*(*)
!\end{verbatim} %+
    WRITE(KOU,10)
10  FORMAT(' You are expected to give a formula that shall be', &
         ' evaluated or manipulated.'/ &
         ' The formula shall be written as a FORTRAN statement with the', &
         ' following rules:'/ &
         ' A variable must begin with a letter', &
         ' and a number with a number (not a dot).'/' A real number must', &
         ' have a dot or an exponent (E).'/ &
         ' The operators + , - ,   , / , ** (exponentiation) can be used'/ &
         ' and any level of parenthesis.'/ &
         ' SQRT(X) is the square root'/ &
         ' EXP(X) is the exponential'/ &
         ' LOG(X) or LN(X) is the natural logarithm'/ &
         ' LOG10(X) is the base 10 logarithm'/ &
         ' SIN(X), COS(X), ATAN(X)'/ &
         ' SIGN(X)'/ &
         ' ERF(X) is the error function'/ &
         ' IVAN(X) Ivantsof function'/ &
         ' BSUM(X) is sum(sin(n*pi*f)**2/(n*pi)**3)'/ &
         ' ABS(X) is absolute value'/ &
         ' HS(X) is the Heaviside function'/ &
         ' Notice that these operators must be followed by a (.'// &
         ' The statement must be terminated by a ;'/)
    RETURN
  END SUBROUTINE EXPHLP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine putprp & Asking for a function
!\begin{verbatim}
  SUBROUTINE PUTPRP(NAMN,MAXS,SYMBOL,PROMPT,ILEN)
!...CREATES A PROMPT asking for a putfun expression with formal arguments
    implicit none
    CHARACTER NAMN*(*),PROMPT*(*),SYMBOL(*)*(*)
    integer ilen,maxs
!...write a prompt with name of all variables
!\end{verbatim} %+
    integer i,j
    ILEN=1
    PROMPT=' '
    CALL CONS(PROMPT,ILEN,NAMN)
    IF(MAXS.LE.0) THEN
       CALL CONS(PROMPT,ILEN,'=')
       GOTO 900
    ENDIF
    CALL CONS(PROMPT,ILEN,'(')
    I=1
11  IF(I.GT.MAXS) GOTO 12
    J=LEN_trim(SYMBOL(I))
    IF(LEN(SYMBOL(I)).GT.J) THEN
       J=J+1
       SYMBOL(I)(J:J)=' '
    ENDIF
    CALL CONS(PROMPT,ILEN,SYMBOL(I)(1:J))
    IF(I.NE.MAXS) CALL CONS(PROMPT,ILEN,',')
    I=I+1
    GOTO 11
12  CONTINUE
    CALL CONS(PROMPT,ILEN,')= ')
    ILEN=ILEN-1
900 RETURN
  END SUBROUTINE PUTPRP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

!\addtotable subroutine delfun & Delete a function
!\begin{verbatim}
  SUBROUTINE DELFUN(LROT,IWS)
!   delete a putfun expression :: not converted to structures
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    integer IWS(*)
    integer lrot
!\end{verbatim}
    integer link(20)
    integer nod,last,lnode
!    DIMENSION LINK(20)
!    LOGICAL SG2ERR
    NOD=LROT
    IF(NOD.LE.2 .OR. NOD.GT.IWS(2)) GOTO 800
    IF(NOD.LE.0) GOTO 900
    LAST=0
!..visit left link
100 LNODE=IWS(NOD+1)
    IF(LNODE.LE.0) GOTO 110
    LAST=LAST+1
    LINK(LAST)=NOD
    NOD=LNODE
    GOTO 100
!..data record at bottom
110 IF(IWS(NOD).EQ.0 .OR. IWS(NOD+3).EQ.1) THEN
!       CALL WRELS(NOD,3+NWPR,IWS)
!       IF(SG2ERR(KERR)) GOTO 900
!       call release_pnode(nod)
       if(pfnerr.ne.0) goto 900
   ELSE
       IWS(NOD+3)=IWS(NOD+3)-1
    ENDIF
!..visit right link
200 IF(LAST.LE.0) GOTO 800
    NOD=LINK(LAST)
    IF(NOD.LT.0) THEN
!..right link visited, remove binary operator, CHECK CODE HERE ...
!       CALL WRELS(-NOD,3,IWS)
!       IF(SG2ERR(KERR)) GOTO 900
    ELSE
!..mark right link is now visited
       LINK(LAST)=-NOD
       IF(IWS(NOD+2).LE.0) THEN
!..remove unary operator
!          CALL WRELS(NOD,3,IWS)
!          IF(SG2ERR(KERR)) GOTO 900
       ELSE
!..set node to this and visit its left link
          NOD=IWS(NOD+2)
          GOTO 100
       ENDIF
    ENDIF
    LAST=LAST-1
    IF(LAST.GT.0) GOTO 200
800 LROT=0
900 RETURN
  END SUBROUTINE DELFUN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
!      HPCALC is a screen HP calculator
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine hpcalc & HP calculator
!\begin{verbatim}
  SUBROUTINE HPCALC
!...EMULATES A HP CALCULATOR ON SCREEN
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!\end{verbatim} %+
    integer, parameter :: NOP=40,MAXPRG=100
    double precision stk(4),reg(0:9)
!    DIMENSION STK(4),REG(0:9)
    CHARACTER LINE*80,INPUT*80,OPER(NOP)*10,CH1*1,CH2*2
    CHARACTER PROG(MAXPRG+1)*20
    LOGICAL PROGT,OK,RUN,TRACE
    SAVE STK,REG,PROG
    integer naxop,lprog,kprog,ip,iback,k,next,i,jprog,last
    double precision ss,val
    DATA STK/0.0,0.0,0.0,0.0/
    DATA REG/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
    DATA &
         OPER( 1)/'BACK      '/,OPER( 2)/'HELP      '/, &
         OPER( 3)/'SHOW_STACK'/,OPER( 4)/'EXP       '/, &
         OPER( 5)/'LN        '/,OPER( 6)/'LOG       '/, &
         OPER( 7)/'SIN       '/,OPER( 8)/'COS       '/, &
         OPER( 9)/'TAN       '/,OPER(10)/'[ASIN     '/, &
         OPER(11)/'[ACOS     '/,OPER(12)/'[ATAN     '/, &
         OPER(13)/'SQRT      '/,OPER(14)/'ROT_STACK '/, &
         OPER(15)/'SWITCH_XY '/,OPER(16)/'POWER_2   '/, &
         OPER(17)/'CLX       '/,OPER(18)/'CLSTACK   '/, &
         OPER(19)/'STO_REG   '/,OPER(20)/'RCL_REG   '/ 
    DATA &
         OPER(21)/'CLEAR_REG '/,OPER(22)/'CHSIGN    '/, &
         OPER(23)/'DISPLAYREG'/,OPER(24)/'          '/, &
         OPER(25)/'_STOP     '/,OPER(26)/'_GOTO     '/, &
         OPER(27)/'_IF       '/,OPER(28)/'_PROGRAM  '/, &
         OPER(29)/'ENTER_PUSH'/,OPER(30)/'_LIST     '/, &
         OPER(31)/'_ERASE_PRO'/,OPER(32)/'_STEP     '/, &
         OPER(33)/'_EDIT     '/,OPER(34)/'_BACK     '/, &
         OPER(35)/'_NOOP     '/,OPER(36)/'_TRACE    '/, &
         OPER(37)/'_END      '/,OPER(38)/'_RUN      '/, &
         OPER(39)/'QUIT      '/,OPER(40)/'FIN       '/
!
    WRITE(*,*)'REVERSE POLISH CALCULATOR'
! uninitiated in original?
    jprog=0
    NAXOP=24
    PROGT=.FALSE.
    RUN=.FALSE.
    OK=.FALSE.
    LPROG=0
    PROG(MAXPRG+1)='_END'
    LAST=LEN(LINE)
100 CONTINUE
    IF(PROGT .AND. OK) THEN
!...SAVE PROGRAM STEP
       LPROG=LPROG+1
       IF(LPROG.GT.MAXPRG) THEN
          WRITE(*,*)'TOO MANY PROGRAM STEPS, MAXIMUM IS ',MAXPRG
       ELSE
          PROG(LPROG)=INPUT(1:20)
          WRITE(*,98)LPROG,INPUT(1:LEN_TRIM(INPUT))
98        FORMAT(' STEP',I4,': ',A)
       ENDIF
    ENDIF
    IF(RUN) THEN
       KPROG=KPROG+1
       IF(KPROG.GT.MAXPRG) THEN
          WRITE(*,*)'ILLEGAL STEP'
          RUN=.FALSE.
          GOTO 100
       ENDIF
       INPUT=PROG(KPROG)
       IF(TRACE)WRITE(*,*)KPROG,STK(1),INPUT
    ELSEIF(PROGT) THEN
!       CALL GPARC('HPP>',LINE,LAST,1,INPUT,' ',HPHLP)
       CALL GPARC_old('HPP>',LINE,LAST,1,INPUT,' ',HPHLP)
       CALL CAPSON(INPUT)
       IF(INPUT(1:1).EQ.' ') GOTO 100
       OK=.FALSE.
    ELSE
       WRITE(*,101)STK(1)
101    FORMAT(1PE15.7)
      CALL GPARC_old('HPC>',LINE,LAST,1,INPUT,' ',HPHLP)
       CALL CAPSON(INPUT)
       IF(INPUT(1:1).EQ.' ') GOTO 100
       OK=.FALSE.
    ENDIF
!...MATH OP
    CH1=INPUT(1:1)
    IF(CH1.EQ.'+') THEN
       STK(1)=STK(1)+STK(2)
       GOTO 102
    ELSEIF(CH1.EQ.'-') THEN
       STK(1)=STK(2)-STK(1)
       GOTO 102
    ELSEIF(CH1.EQ.'*') THEN
       STK(1)=STK(1)*STK(2)
       GOTO 102
    ELSEIF(CH1.EQ.'/') THEN
       STK(1)=STK(2)/STK(1)
       GOTO 102
    ELSEIF(CH1.EQ.'^') THEN
       STK(1)=STK(2)**STK(1)
       GOTO 102
    ENDIF
    GOTO 109
!...SHIFT STACK DOWN
102 STK(2)=STK(3)
    STK(3)=STK(4)
    OK=.TRUE.
    GOTO 100
!...NUMBER OR OPCODE
109 IP=0
    buperr=0
!    write(*,*)'number 1: ',input(1:10),ip
    CALL GETREL(INPUT,IP,VAL)
!    write(*,*)'number 2: ',input(1:10),ip,buperr,val
    IF(buperr.eq.0) THEN
!...	   NUMBER, SAVE ON STACK
       STK(4)=STK(3)
       STK(3)=STK(2)
       STK(2)=STK(1)
       STK(1)=VAL
       OK=.TRUE.
!       write(*,*)'pushed: ',ip
       IF(INPUT(IP:IP).NE.' ') THEN
!...to allow input like 3000 3 4 5*/+
          IBACK=LEN_TRIM(INPUT)-IP+1
          LAST=LAST-IBACK-1
       ENDIF
!       write(*,*)'pushed: ',ip
       GOTO 100
    ENDIF
!    CALL RESERR
    K=NCOMP2(INPUT,OPER,NOP,NEXT)
    IF(K.EQ.0) THEN
       WRITE(*,*)'NO SUCH OPCODE'
       RUN=.FALSE.
       GOTO 100
    ELSEIF(K.LT.0) THEN
       WRITE(*,*)'AMBIGUOUS OPCODE'
       RUN=.FALSE.
       GOTO 100
    ELSE
       OK=.TRUE.
       GOTO(110,120,130,140,150,160,170,180,190,200, &
            210,220,230,240,250,260,270,280,290,300, &
            310,320,380,340,350,360,370,330,390,400, &
            410,420,430,440,450,460,470,480,490,500),K
    ENDIF
!...EXIT
110 RETURN
!...HELP
120 INPUT(NEXT:NEXT+1)=',,'
!    CALL NGHELP(INPUT,NEXT,OPER,NAXOP) routine removed
    OK=.FALSE.
    GOTO 100
!...SHOW STACK
130 WRITE(*,131)STK
131 FORMAT(4(1PE15.7))
    GOTO 100
!...EXP
140 STK(1)=EXP(STK(1))
    GOTO 100
!...LN
150 STK(1)=LOG(STK(1))
    GOTO 100
!...LOG10
160 STK(1)=LOG10(STK(1))
    GOTO 100
!...SIN
170 STK(1)=SIN(STK(1))
    GOTO 100
!...COS
180 STK(1)=COS(STK(1))
    GOTO 100
!...TAN
190 STK(1)=SIN(STK(1))/COS(STK(1))
    GOTO 100
!...ASIN
200 CONTINUE
!...ACOS
210 CONTINUE
!...ATAN
220 CONTINUE
    WRITE(*,*)'NOT IMPLEMENTED'
    GOTO 100
!...SQRT
230 STK(1)=SQRT(STK(1))
    GOTO 100
!...ROT
240 CONTINUE
    SS=STK(4)
    STK(4)=STK(3)
    STK(3)=STK(2)
    STK(2)=STK(1)
    STK(1)=SS
    GOTO 100
!...SWITCH_XY
250 CONTINUE
    SS=STK(2)
    STK(2)=STK(1)
    STK(1)=SS
    GOTO 100
!...POWER_2
260 CONTINUE
    STK(1)=STK(1)**2
    GOTO 100
!...CLX
270 STK(1)=0.0
    GOTO 100
!...CLEAR
280 STK(1)=0.0
    STK(2)=0.0
    STK(3)=0.0
    STK(4)=0.0
    GOTO 100
!...STO_REG <N>
290 IF(RUN) THEN
       CALL GETINT(INPUT,NEXT,I)
    ELSE
       CALL GPARI_old('REG#',LINE,LAST,I,-1,NOHELP)
    ENDIF
    IF(I.LT.0 .OR. I.GT.9) THEN
       WRITE(*,*)'REGISTER NUMBER MUST BE 0..9'
    ELSE
       REG(I)=STK(1)
       IF(PROGT) INPUT(NEXT+1:NEXT+1)=CHAR(I+ICHAR('0'))
    ENDIF
    GOTO 100
!...RCL_REG <N>
300 IF(RUN) THEN
       CALL GETINT(INPUT,NEXT,I)
    ELSE
       CALL GPARI_old('REG#',LINE,LAST,I,-1,NOHELP)
    ENDIF
    IF(I.LT.0 .OR. I.GT.9) THEN
       WRITE(*,*)'REGISTER NUMBER MUST BE 0..9'
    ELSE
       STK(1)=REG(I)
       IF(PROGT) INPUT(NEXT+1:NEXT+1)=CHAR(I+ICHAR('0'))
    ENDIF
    GOTO 100
!...CLR_REG
310 DO I=0,9
       REG(I)=0.0
    ENDDO
    GOTO 100
!...CHSIGN
320 STK(1)=-STK(1)
    GOTO 100
!...PROGRAM
330 IF(PROGT)THEN
       WRITE(*,*)'ALREDAY SET'
    ELSE
       PROGT=.TRUE.
       LPROG=0
       RUN=.FALSE.
    ENDIF
    OK=.FALSE.
    GOTO 100
!...
340 CONTINUE
    GOTO 100
!...STOP
350 IF(RUN) THEN
       RUN=.FALSE.
       WRITE(*,*)'PROGRAM STOP AT ',KPROG
    ELSEIF(PROGT) THEN
       PROGT=.FALSE.
       LPROG=LPROG+1
       PROG(LPROG)='STOP'
       WRITE(*,*)'PROGRAM STOP AT ',LPROG
    ENDIF
    GOTO 100
!...GOTO
360 IF(RUN) THEN
       CALL GETINT(INPUT,NEXT,IP)
       IF(IP.EQ.0) THEN
          WRITE(*,*)'PROGRAM STOP AT 0'
          RUN=.FALSE.
       ELSEIF(IP.LE.0 .OR. IP.GT.100 .OR. IP.EQ.KPROG) THEN
          WRITE(*,*)'ILLEGAL GOTO ADDRESS IN STEP ',KPROG
          RUN=.FALSE.
       ELSE
          KPROG=IP-1
       ENDIF
    ELSE
       CALL GPARI_old('STEP ',LINE,LAST,IP,0,NOHELP)
       IF(IP.LT.0 .OR. IP.GT.MAXPRG) THEN
          WRITE(*,*)'ILLEGAL ADDRESS'
          OK=.FALSE.
       ELSE
          WRITE(INPUT(NEXT+1:NEXT+3),365)IP
365       FORMAT(I3)
       ENDIF
    ENDIF
    GOTO 100
!...IF >=<
370 IF(RUN) THEN
       IF(INPUT(NEXT+1:NEXT+1).EQ.'>') THEN
          IF(STK(1).GT.STK(2)) KPROG=KPROG+1
       ELSEIF(INPUT(NEXT+1:NEXT+1).EQ.'=') THEN
          IF(ABS(STK(1)-STK(2)).LT.1.0D-16) KPROG=KPROG+1
       ELSEIF(INPUT(NEXT+1:NEXT+1).EQ.'<') THEN
          IF(STK(1).LT.STK(2)) KPROG=KPROG+1
       ELSE
          WRITE(*,*)'ILLEGAL CONDITION IN STEP',KPROG
          RUN=.FALSE.
       ENDIF
       OK=.TRUE.
    ELSE
       CALL GPARC_old('CONDITION ( > = OR < )',LINE,LAST,1, &
            CH2,'NONE',NOHELP)
       IF(CH2.EQ.'> ' .OR. CH2.EQ.'= ' .OR. CH2.EQ.'< ') THEN
          INPUT(NEXT+1:)=CH2
       ELSE
          WRITE(*,*)'ILLEGAL CONDITION'
          OK=.FALSE.
       ENDIF
    ENDIF
    GOTO 100
!...DISPLAY_REG
380 WRITE(*,381)REG
381 FORMAT(5(1PE15.7)/5(1PE15.7)/)
    GOTO 100
!...ENTER, PUSH STACK
390 STK(4)=STK(3)
    STK(3)=STK(2)
    STK(2)=STK(1)
    GOTO 100
!...LIST
400 WRITE(*,401)(I,PROG(I)(1:LEN_TRIM(PROG(I))),I=1,LPROG)
401 FORMAT(I4,': ',A)
    OK=.FALSE.
    GOTO 100
!...ERASE
410 LPROG=0
    OK=.FALSE.
    GOTO 100
!...STEP
420 WRITE(*,*)'NOT IMPLEMENTED'
    GOTO 100
    KPROG=KPROG+1
    IF(KPROG.GT.LPROG) THEN
       WRITE(*,*)'PROGRAM ENDS AT ',LPROG
    ELSE
       WRITE(*,401)KPROG,PROG(KPROG)
    ENDIF
    OK=.FALSE.
    GOTO 100
!...EDIT
430 IF(JPROG.EQ.0) JPROG=LPROG-2
    CALL GPARI_old('STEP ',LINE,LAST,I,JPROG+1,NOHELP)
    WRITE(*,*)'NO SUCH STEP'
    OK=.FALSE.
    GOTO 100
!...BACK
440 WRITE(*,*)'NOT IMPLEMENTED'
    RETURN
!    GOTO 100
    WRITE(*,401)KPROG,PROG(KPROG)
    IF(KPROG.GT.1) KPROG=KPROG-2
    OK=.FALSE.
!...NOOP
450 GOTO 100
!...TRACE
460 CONTINUE
    IF(TRACE) THEN
       TRACE=.FALSE.
    ELSE
       TRACE=.TRUE.
    ENDIF
    GOTO 100
!...END
470 PROGT=.FALSE.
    LPROG=LPROG+1
    PROG(LPROG)='STOP'
    OK=.FALSE.
    GOTO 100
!...RUN
480 IF(LPROG.EQ.0) THEN
       WRITE(*,*)'NO PROGRAM'
    ELSE
       RUN=.TRUE.
       PROGT=.FALSE.
       KPROG=0
    ENDIF
    OK=.FALSE.
    GOTO 100
!...QUIT
490 CONTINUE
    RETURN
!...FIN
500 CONTINUE
    RETURN
!    GOTO 100
  END SUBROUTINE HPCALC

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine hphelp & Help to HP calculator
!\begin{verbatim}
  SUBROUTINE HPHLP
! writes a help text for using the online HP calculator
    implicit none
!\end{verbatim}
    WRITE(*,10)
10  FORMAT(' This is a revese polish calculator'/&
         ' Input are numbers, + - * / and ^ and OPCODEs.',&
         ' Use HELP to list OPCODEs.'/' Several numbers and operations',&
         ' can be given on one line.'/' The content of the X register',&
         ' is displayed after each operation'//&
         ' Example input: 30000 8 1273 * / chs 1.5 3 ^ + exp 2 *'/&
         ' Computes 2*EXP(1.5**3-30000/(8*1273))'//)
    RETURN
  END SUBROUTINE HPHLP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!
! >>>> subsection
! WPACK can convert from an integer workspace to normal double/character
! used to save data on an unformatted Fortran file
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine winit & Initate integer workspace
!\begin{verbatim}
  SUBROUTINE WINIT(NWT,NWR,IWS)
!...INITIATES A WORKSPACE
! INPUT: NWT IS THE DIMENSION OF THE WORKSPACE
!        NWR IS THE NUMBER OF WORDS TO BE EXCLUDED IN THE BEGINNING
!        IWS IS THE WORKSPACE
! EXIT:  THE FREE LIST IS INITIATED IN IWS
! ERRORS: NWR LESSER THAN ZERO
!         NWT LESSER THAN NWR+100
    implicit none
    integer nwt,nwr,iws(*)
!    DIMENSION IWS(*)
!\end{verbatim} %+
    integer nwres,ifri
    IF(NWR.LT.0) GOTO 910
    IF(NWT.LT.NWR+100) GOTO 920
    NWRES=NWR+3
!...IWS(1) IS PUT TO THE FIRST FREE AREA IN THE WORKSPACE AND IWS(2)
!      TO THE SIZE OF THE WORKSPACE. THE APPLICATION PROGRAM MUST NOT
!      CHANGE THESE LOCATIONS!
    IWS(1)=NWRES
    IWS(2)=NWT
!...PUT ALL WORDS FROM 3 TO NWRES TO ZERO
!      THIS INCLUDES THE FIRST WORD IN THE FREE AREA
    DO IFRI=3,NWRES
       IWS(IFRI)=0
    enddo
!...THE SECOND WORD IN THE FREE AREA IS PUT THE THE NUMBER OF FREE WORDS THIS
!      NUMBER IS NWT-NWR-(TWO WORDS IN THE BEGINNING)-(TWO WORDS IN THE END)
    IWS(NWRES+1)=NWT-NWR-4
900 RETURN
910 buperr=1008
    goto 900
920 buperr=1002
    GOTO 900
  END SUBROUTINE WINIT

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
  
!\addtotable subroutine wold & Read an integer workspace from file
!\begin{verbatim}
  SUBROUTINE WOLD(FIL,NW,IWS)
!...READS A FILE INTO A WORKSPACE. THE FILE MUST HAVE BEEN WRITTEN BY WSAVE
! INPUT: FIL A CHARACTER WITH A LEGAL FILE NAME
!        NW THE DIMENSION OF IWS
!        IWS THE WORKSPACE
! CALLS: WRKCHK TO CHECK THE FREE LIST
! EXIT:  THE CONTENT OF THE FILE IS IN IWS. THE DIMENSION OF IWS IS SET TO
!            NW AND THE LAST FREE AREA IS CORRECTED
    implicit none
    CHARACTER FIL*(*)
    integer nw,iws(*)
!    DIMENSION IWS(*)
!\end{verbatim} %+
    integer ierr,last,j,k
    OPEN(UNIT=LUN,FILE=FIL,ACCESS='SEQUENTIAL',STATUS='OLD',&
         IOSTAT=IERR,ERR=910,FORM='UNFORMATTED')
! note: first integer on file is size of unformatted file
    READ(LUN,END=100,ERR=100)J,(IWS(K),K=1,J)
!...CHECK THE WORKSPACE
    CALL WRKCHK(LAST,NW,IWS)
100 CLOSE(LUN)
    RETURN
910 continue
    GOTO 100
  END SUBROUTINE WOLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine wsave & Save integer worspace to file
!\begin{verbatim}
!  SUBROUTINE WSAVE(FIL,NW,IWS)
  SUBROUTINE WSAVE(FIL,IWS)
!...WRITES A WORKSPACE ON A FILE
! INPUT: FIL IS A CHARACTER WITH A LEGAL FILE NAME
!        NW IS THE DIMENSION OF THE WORKSPACE
!        IWS IS THE WORKSPACE
! CALLS: WRKCHK TO CHECK THE WORKSPACE
! ERROR: IF THE WORKSPACE IS INCORRECT IT CANNOT BE SAVED
    implicit none
!    integer nw,iws(*)
    integer iws(*)
!    DIMENSION IWS(*)
    CHARACTER FIL*(*)
!\end{verbatim}
    integer i,ierr,last
    I=IWS(2)
    CALL WRKCHK(LAST,I,IWS)
!    IF(SG2ERR(IERR)) GOTO 900
    if(buperr.ne.0) goto 900
    OPEN(UNIT=LUN,FILE=FIL,ACCESS='SEQUENTIAL',STATUS='UNKNOWN',&
         IOSTAT=IERR,ERR=910,FORM='UNFORMATTED')
! note: first integer on file is size of unformatted file
    WRITE(LUN,ERR=910)LAST,(IWS(I),I=1,LAST)
800 CLOSE(LUN)
900 RETURN
910 continue
    GOTO 800
  END SUBROUTINE WSAVE

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine wpatch & Patch an integer workspace
!\begin{verbatim}
  SUBROUTINE WPATCH(NW,IWS)
!...ROUTINE TO PATCH A WORKSPACE
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    integer nw,iws(*)
!    DIMENSION IWS(*)
!\end{verbatim} %+
    integer idum,ip,iadr,ival,j
    CHARACTER LINE*80,CHX*(NBPW),CHHEX*(2*NBPW)
    double precision x,z
    IF(IWS(2).NE.NW) THEN
       WRITE(KOU,*)' WORKSPACE DIMENSION INCORRECT, SET TO ',NW
       IWS(2)=NW
    ENDIF
    CALL WRKCHK(IDUM,NW,IWS)
!    IF(SG1ERR(IERR)) THEN
    if(buperr.ne.0) then
       WRITE(KOU,*)' YOU MAY ATTEMPT TO CORRECT THE FREE LIST'
!       CALL RESERR
       buperr=0
    ENDIF
10  WRITE(KOU,*)' ADDRESS: '
    CALL BINTXT(KIU,LINE)
    IP=1
    CALL GETINT(LINE,IP,IADR)
    if(buperr.ne.0) then
       buperr=0
!       IF(LINE(1:1).EQ.'?') CALL WPHLP(IP,LINE)
       IF(LINE(1:1).EQ.'?') CALL WPHLP
       IF(LINE(1:1).EQ.'@') GOTO 900
       WRITE(KOU,20)
20     FORMAT(' TYPE ? FOR HELP'/)
       GOTO 10
    ENDIF
    WRITE(KOU,30)
30  FORMAT(' ADDRESS           INTEGER  CHAR  HEXADEC.  REAL VALUE'/)
100 IF(IADR.LT.1 .OR. IADR.GT.NW) THEN
       WRITE(KOU,*)'OUTSIDE WORKSPACE',1,NW
       GOTO 900
    ENDIF
    CALL LOADR(1,IWS(IADR),X)
    CALL LOADC(NBPW,IWS(IADR),CHX)
    CALL WRIHEX(CHHEX,IWS(IADR))
    DO J=1,NBPW
!...MACHDEP REPLACEMENT OF NON-PRINTABLE ASCII CHARACTERS WITH A PERIOD
       IF(LLT(CHX(J:J),' ') .OR. LGT(CHX(J:J),'~')) CHX(J:J)='.'
    enddo
    WRITE(KOU,110,ERR=911)IADR,IWS(IADR),CHX,CHHEX,X
110 FORMAT('$',I7,5X,I15,2X,A,2X,A,2X,E15.8)
111 CALL GPARR_old('NEW VALUE: ',LINE,IP,Z,RNONE,WPHLP)
!    IF(SG2ERR(IERR)) THEN
!...      NOT A DIGIT: EXIT, STORE AS BYTES, OCTAL OR IGNORE
    if(buperr.ne.0) then
       buperr=0
!       CALL RESERR
       IF(LINE(1:1).EQ.'@') GOTO 900
       IF(LINE(1:2).EQ.'EX'.OR.LINE(1:2).EQ.'ex') GOTO 900
       IF(LINE(1:1).EQ.'"') THEN
          CALL STORC(NBPW,IWS(IADR),LINE(2:))
          IADR=IADR+1
       ELSEIF(LINE(1:1).EQ.'&') THEN
! OCTAL VALUE
          IP=2
          CALL GETOCT(LINE,IP,IVAL)
          if(buperr.ne.0) then
             buperr=0
!          IF(SG2ERR(IERR)) THEN
!             CALL RESERR
             WRITE(KOU,*)'VALUE AFTER & NOT OCTAL'
          ELSE
             IWS(IADR)=IVAL
             IADR=IADR+1
          ENDIF
       ELSEIF(LINE(1:1).EQ.'#') THEN
! HEXADECIMAL VALUE
          IP=2
          CALL GETHEX(LINE,IP,IVAL)
!          IF(SG2ERR(IERR)) THEN
!             CALL RESERR
          if(buperr.ne.0) then
             buperr=0
             WRITE(KOU,*)'VALUE AFTER # NOT HEXADECIMAL'
          ELSE
             IWS(IADR)=IVAL
             IADR=IADR+1
          ENDIF
       ELSEIF(EOLCH(LINE,IP)) THEN
          IADR=IADR+1
       ELSE
          WRITE(KOU,20)
       ENDIF
    ELSE
! DIGIT
       IF(LINE(IP:IP).EQ.'/') THEN
! NEW ADDRESS
          IADR=INT(Z)
       ELSEIF(INDEX(LINE,'.').GT.0) THEN
! REAL VALUE
          CALL STORR(1,IWS(IADR),Z)
          IADR=IADR+NWPR
       ELSE
! INTEGER VALUE
          IWS(IADR)=INT(Z)
          IADR=IADR+1
       ENDIF
    ENDIF
! SKIP REST OF LINE
    IP=LEN(LINE)
    GOTO 100
900 CALL WRKCHK(IDUM,NW,IWS)
    RETURN
911 continue
    GOTO 111
  END SUBROUTINE WPATCH

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine wphlp & Help to patch workspace
!\begin{verbatim}
!  SUBROUTINE WPHLP(ITYP,LINE)
  SUBROUTINE WPHLP
!...HELP ROUTINE FOR WPATCH
    implicit none
!    CHARACTER LINE*(*)
!    integer ityp
!\end{verbatim} %+
    WRITE(KOU,10)
10  FORMAT(' YOU MAY PATCH THE WORKSPACE.'/&
         ' THE VALUE AT THE SPECIFIED ADDRESS IN THE WORKSPACE',&
         ' IS DISPLAYED AS'/' INTEGER, CHAR (NON-PRINTABLE',&
         ' CHAR REPLACED BY .) HEXADECIMAL AND REAL.'//&
         ' THE FOLLOWING INPUT IS LEGAL:'/&
         ' <CR>            VALUE IN NEXT ADDRESS IS DISPLAYED'/&
         ' <NUMBER>        <NUMBER> IS STORED AT THE ADDRESS'/&
         '        A REAL NUMBER MUST INCLUDE A PERIOD (.)'/&
         ' <NUMBER>/       <NUMBER> IS TAKEN AS NEW ADDRESS'/&
         ' &<OCTAL NUMBER> <NUMBER> STORED AS OCTAL'/&
         ' #<HEX NUMBER>   <NUMBER> STORED AS HEXADECIMAL'/&
         ' "<TEXT>         <TEXT> STORED AS BYTES',&
         ' (BYTES FOR ONE WORD ONLY)'/&
         ' @ OR EXIT       EXIT'/&
         ' ?               THIS TEXT'/&
         ' <ANYTHING ELSE> IGNORED'/)
    RETURN
  END SUBROUTINE WPHLP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine wrkchk & Check consistency of workspace
!\begin{verbatim}
  SUBROUTINE WRKCHK(LAST,NW,IWS)
!...CHECKS THE FREE LIST IN A WORKSPACE
! INPUT: NW IS THE DIMENSION
!        IWS IS THE WORKSPACE
! EXIT:  LAST IS PUT TO THE LAST WORD USED IN THE WORKSPACE
! ERRORS: ANY ERROR IN THE FREE LIST (POINTER OUTSIDE WORKSPACE ETC)
    implicit none
    integer last,nw,iws(*)
!    DIMENSION IWS(*)
!\end{verbatim} %+
    integer lok,lfr
    IF(NW.LT.100) GOTO 910
    IWS(2)=NW
!...SERACH THE FREE LIST STARTING IN WORD 1
!      (THERE MUST ALWAYS BE A FREE AREA!)
    LOK=IWS(1)
!...A FREE AREA MUST LIE BETWEEN 2 AND THE DIMENSION
100 IF(LOK.LE.2 .OR. LOK.GE.NW) GOTO 920
    LAST=LOK
    LOK=IWS(LOK)
!...IN THE LAST FREE AREA LOK=0
    IF(LOK.EQ.0) GOTO 200
!...THE FREE AREAS ARE ALWAYS ORDERD INCREASINGLY
    IF(LOK.LT.LAST+2) GOTO 930
    LFR=IWS(LAST+1)
!...A FREE AREA IS AT LEAST TWO WORDS AND NOT PAST THE NEXT AREA
    IF(LFR.LT.2 .OR. LAST+LFR.GT.LOK) GOTO 940
    GOTO 100
!...THE FREE AREA SEEMS CORRECT
200 LFR=LAST+1
    IWS(LFR)=NW-LFR
    LAST=LFR
900 RETURN
910 continue
    buperr=1002
    GOTO 900
920 continue
    buperr=1003
    GOTO 900
930 continue 
    buperr=1004
   GOTO 900
940 continue
    buperr=1005
    GOTO 900
  END SUBROUTINE WRKCHK

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine wlist & List free list in worspace
!\begin{verbatim}
  SUBROUTINE WLIST(IWS)
!...LISTS THE FREE AREAS
    implicit none
    integer iws(*)
!    DIMENSION IWS(*)
!\end{verbatim}
    integer n,nw,nwp
    N=1
    NW=0
    WRITE(KOU,10)IWS(2)
10  FORMAT(/' MAP OF THE FREE SPACE CONTAINING',I12,' WORDS')
100 N=IWS(N)
    IF(N.LE.0) GOTO 200
    NWP=IWS(N+1)
    NW=NW+NWP
    WRITE(KOU,110)N,NWP
110 FORMAT(' FROM ',I12,' ARE ',I12,' WORDS FREE')
    GOTO 100
200 WRITE(KOU,210)NW
210 FORMAT(/' TOTAL NUMBER OF FREE WORDS ARE',I12/)
    RETURN
  END SUBROUTINE WLIST

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine wtrest & Reserve rest of workspace
!\begin{verbatim}
  SUBROUTINE WTREST(NYB,NW,IWS)
!...RESERVES THE LAST PART OF THE WORKSPACE
! INPUT: IWS IS A WORKSPACE
! EXIT:  NYB IS A POINTER TO THE RESERVED PART
!        NW IS THE NUMBER OF RESERVED WORDS
    implicit none
    integer nyb,nw,iws(*)
!    DIMENSION IWS(*)
!\end{verbatim} %+
    integer lok,last
    LOK=1
100 LAST=LOK
    LOK=IWS(LAST)
    IF(LOK.GT.0) GOTO 100
    NW=IWS(LAST+1)
    CALL WTAKE(NYB,NW,IWS)
    RETURN
  END SUBROUTINE WTREST

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine wtake & Reserve a record in workspace
!\begin{verbatim}
  SUBROUTINE WTAKE(NYB,NW,IWS)
!......RESERVS NW WORDS IN THE WORKSPACE
! INPUT: NW IS THE NUMBER OF WORDS TO BE RESERVED
!        IWS IS THE WORKSPACE
! EXIT:  NYB POINTS TO THE FIRST WORD THAT IS RESERVED
! ERROR: TOO SMALL OR TOO LARGE NUMBER OF WORDS TO BE RESERVED
    implicit none
    integer nyb,nw,iws(*)
!    DIMENSION IWS(*)
!...THE FREE LIST START IN THE FIRST WORD
!      IN EACH FREE AREA THE FIRST WORD POINTS TO THE NEXT FREE AREA
!      AND THE SECOND GIVES THE NUMBER OF WORDS IN THIS AREA
!      THE FREE LIST ENDS WITH THE POINTER EQUAL TO ZERO
!\end{verbatim} %+
    integer loka,lokb,next
    IF(NW.LT.2) GOTO 910
    LOKB=1
4   LOKA=IWS(LOKB)
    IF(LOKA.LE.0) GOTO 920
    IF(LOKA.GE.IWS(2)) GOTO 930
! deleted feature !!!  if(X) <0 first label, =0 second label >0 third label
!    IF(IWS(LOKA+1)-NW) 10,20,30
    IF(IWS(LOKA+1)-NW .eq.0) then
       goto 20
    elseif(iws(loka+1)-nw.gt.0) then
       goto 30
    endif
!...TOO SMALL AREA, CONTINUE WITH THE NEXT
10  LOKB=LOKA
    GOTO 4
!...EXACT FIT WITH THE REQUESTED NUMBER OF WORDS
!      IF IWS(LOKA)=0 IT IS THE LAST FREE AREA OF THE WORKSPACE.
!      AS WRELS WILL NOT WORK PROPERLY IF THERE IS NOT A FREE AREA AFTER THE
!      LAST RESERVED A POINTER IS SET. IN WINIT THE SIZE OF THE LAST FREE AREA
!      WAS DECREASED BY TWO TO LEAVE A UNRESERVABLE FREE AREA LAST
20  IF(IWS(LOKA).GT.0) THEN
       IWS(LOKB)=IWS(LOKA)
    ELSE
       GOTO 31
    ENDIF
    GOTO 50
!...LARGER AREA THAN REQUESTED
!      A FREE AREA MUST BE AT LEAST TWO WORDS, IF THIS AREA IS AT LEAST
!      TWO WORDS LARGER THAN THE REQUEST IT IS DIVIDED, OTHERWISE SKIPPED
30  IF(IWS(LOKA+1)-NW-2.LT.0) GOTO 10
31  NEXT=LOKA+NW
    IWS(NEXT)=IWS(LOKA)
    IWS(NEXT+1)=IWS(LOKA+1)-NW
    IWS(LOKB)=NEXT
50  NYB=LOKA
!...THE RESERVED AREA IS ZEROED
    LOKB=NYB+NW-1
    DO LOKA=NYB,LOKB
       IWS(LOKA)=0
    enddo
900 RETURN
910 continue
    buperr=1008
    GOTO 900
920 continue
    buperr=1006
    GOTO 900
930 continue
    buperr=1007
    GOTO 900
  END SUBROUTINE WTAKE

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine wrels & Release a record in workspace
!\begin{verbatim}
  SUBROUTINE WRELS(IDP,NW,IWS)
!......Returns NW words beginning from IDP to the free workspace list
!      The free workspace list is in increasing order
!      IWS(1) points to the first free space
!      IWS(2) gives the total number of words in the workspace
    implicit none
!    DIMENSION IWS(*)
    integer idp,nw,iws(*)
!......Check that the released space is at lest 2 words and that it is
!      inside the workspace (That is between 3 and IWS(2))
!\end{verbatim}
    integer loka,lokb,lokc
    IF(IDP.LT.3.OR.IDP.GE.IWS(2).OR.NW.LT.2.OR.NW.GE.IWS(2)) GOTO 910
    LOKC=IDP
    LOKB=1
100 LOKA=LOKB
    LOKB=IWS(LOKA)
    IF(LOKB.LE.0) GOTO 920
    IF(LOKB.LT.LOKC) GOTO 100
!..LOKA is the address of the nearest free space below LOKC
    IF(LOKA.EQ.1) GOTO 120
!..Check if the two areas can be merged
!    IF(LOKA+IWS(LOKA+1)-LOKC) 120,110,930
    IF(LOKA+IWS(LOKA+1)-LOKC .lt.0) then
       goto 120
    ELSEIF(LOKA+IWS(LOKA+1)-LOKC .gt.0) then
       goto 930
    endif
!..The released space follows directly on LOKA => Merge LOKA and LOKC
    LOKC=LOKA
    IWS(LOKC+1)=IWS(LOKC+1)+NW
    GOTO 130
!..Set the pointer from LOKC to LOKB and from LOKA to LOKC
120 IWS(LOKC)=LOKB
    IWS(LOKA)=LOKC
    IWS(LOKC+1)=NW
!..Check if LOKC now can be merged with LOKB!
! deleted fetaure
!130 IF(LOKC+IWS(LOKC+1)-LOKB) 900,140,940
130 continue
    IF(LOKC+IWS(LOKC+1)-LOKB .lt.0) then
       goto 900
    elseif(LOKC+IWS(LOKC+1)-LOKB.gt.0) then
       goto 940
    endif
!..Merge LOKC and LOKB
    IWS(LOKC)=IWS(LOKB)
    IWS(LOKC+1)=IWS(LOKC+1)+IWS(LOKB+1)
900 RETURN
!...ERRORS
! TOO SMALL OR OUTSIDE WORKSPACE
910 buperr=1008
    GOTO 900
! ABOVE HIGHEST FREE WORKSPACE
920 buperr=1008
    GOTO 900
! FIRST PART ALREADY FREE
930 GOTO 920
! LAST PART ALREADY FREE
940 GOTO 920
  END SUBROUTINE WRELS

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function nwch & Number of words to store a character
!\begin{verbatim}
  INTEGER FUNCTION NWCH(NB)
! number of words to store a character with nb bytes
! nbpw is the number of bytes in a word.  If not even multiple add 1 word
    implicit none
    integer nb
!\end{verbatim} %+
    integer i
    i=nb/nbpw
    if(mod(nb,nbpw).gt.0) then
       i=i+1
    endif
    nwch=i
    return
  end FUNCTION NWCH

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine storc & Store a character in workspace
!\begin{verbatim}
  SUBROUTINE STORC(N,IWS,text)
! Stores a text character in an integer workspace at position N
! The length of the character to store is len(text)
!\end{verbatim} %+
    implicit none
    integer n,iws(*)
    character text*(*)
! maximal size of character, note used also to store functions and bibliography
    integer, parameter :: maxchar=2048,maxequiv=512
! NOTE BELOW DIMENSIONING BLEOW, maxchar=nbpw*maxequiv
!    character (len=:), allocatable :: localtxt
!    integer, allocatable, dimension(:) :: localint
    character*(maxchar) localtxt
! assumed 32 bit integres, 8 bits/character, 4 characters/word =nbpw
    integer llen,j,now
    integer localint(maxequiv)
! equivalence can only be made between local unallocated variables
    equivalence (localtxt,localint)
!    equivalence (localtxt2,localint2)
    if(maxchar.ne.nbpw*maxequiv) then
       write(*,*)'METLIB utility error: maxchar and maxequiv do not match'
       stop
    endif
    llen=len(text)
    if(llen.gt.maxchar) then
       write(*,*)'Cannot store texts larger than ',maxchar,' characters'
       buperr=1010; goto 900
    endif
! due to the equivalence this stores the character bit map into localint2 !!
! the localtext will be padded with spaces after text
    localtxt=text
! number of words to store rounding off?? integers are 4 bytes (32 bits)
    now=nwch(llen)
    do j=1,now
       iws(n+j-1)=localint(j)
    enddo
!    localint2=localint
!    write(*,800)llen,now,text(1:llen),localtxt(1:llen),localtxt2(1:llen)
!800 format('storc: ',2i4,3('"',a),'"')
900 continue
    return
  end SUBROUTINE STORC

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine loadc & Load a character from workspace
!\begin{verbatim}
  SUBROUTINE LOADC(N,IWS,text)
! copies a text from an integer workspace at position N into a character
! The number of characters to copy is len(text)
    implicit none
    integer n,iws(*)
    character text*(*)
!    character (len=:), allocatable :: localtxt
!    integer, allocatable, dimension(:) :: localint
! maximal size of character, note used also to store functions and bibliography
!\end{verbatim} %+
    integer, parameter :: maxchar=2048,maxequiv=512
! NOTE BELOW DIMENSIONING BELOW, maxchar=nbpw*maxequiv
    character*(maxchar) localtxt
! assumed 32 bit integer, 8 bits character, 4 characters/word
    integer llen,j,now
    integer localint(maxequiv)
! equivalence can obly be made between local unallocated variables
    equivalence (localtxt,localint)
    llen=len(text)
    if(llen.gt.maxchar) then
       write(*,*)'Attempt to extract a text larger than ',maxchar
       buperr=1010;; goto 900
    endif
    now=nwch(llen)
    do j=1,now
       localint(j)=iws(n+j-1)
    enddo
!    write(*,800)llen,now,localtxt(1:llen)
!800 format('LOADC: ',2i3,' "',a,'"')
    text=localtxt(1:llen)
900 continue
    return
  end SUBROUTINE LOADC

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine storr & Store a double in workspace
!\begin{verbatim}
  SUBROUTINE STORR(N,IWS,VALUE)
!...STORES A REAL NUMBER IN A WORKSPACE at index N
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!    DIMENSION IWS(*)
    implicit none
    integer iws(*)
    double precision value
    integer n
!\end{verbatim} %+
    INTEGER JWS(2),int(2)
    DOUBLE PRECISION WS,aws
    EQUIVALENCE (WS,JWS),(int,aws)
! move the exact bit pattern from real VALUE to integer IWS(N)
    WS=VALUE
    IWS(N)=JWS(1)
    IWS(N+1)=JWS(2)
!    int=jws
!    write(*,17)value,ws,aws
!17  format('storr: ',3(1pe14.6),/10x,4i14)
    RETURN
  END SUBROUTINE STORR

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine loadr & Load a double from workspace
!\begin{verbatim}
  SUBROUTINE LOADR(N,IWS,VALUE)
!...LOADS A REAL NUMBER FROM A WORKSPACE at index N
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    integer iws(*)
!    DIMENSION IWS(*)
    DOUBLE PRECISION VALUE
    integer N
!\end{verbatim} %+
    DOUBLE PRECISION WS
    INTEGER JWS(2)
    EQUIVALENCE (WS,JWS)
! move the exact bit pattern from integer IWS(N) to real VALUE
    JWS(1)=IWS(N)
    JWS(2)=IWS(N+1)
    VALUE=WS
    RETURN
  END SUBROUTINE LOADR

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine storrn & Store N doubles in workspace
!\begin{verbatim}
  SUBROUTINE STORRN(N,IWS,ARR)
! store N doubles in workspace
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
!    DIMENSION IWS(*),ARR(*)
    integer n,iws(*)
    double precision arr(*)
!\end{verbatim} %+
    integer, parameter :: maxr=256
    double precision dlocal(maxr)
    integer ilocal(maxr*nwpr)
    integer i
    equivalence (dlocal,ilocal)
!    if(n.gt.256) then
    if(n.gt.512) then
       write(*,*)'M4 STORRN cannot handle arrays larger than ',maxr,n
       buperr=1010; goto 900
    endif
    do i=1,n
       dlocal(i)=arr(i)
    enddo
    DO I=1,N*nwpr
       iws(I)=ilocal(I)
    enddo
900 continue
    RETURN
  END SUBROUTINE STORRN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine loadrn & Load N doubles frm workspace
!\begin{verbatim}
  SUBROUTINE LOADRN(N,IWS,ARR)
! load N doubles from workspace
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    double precision ARR(*)
    integer n,iws(*)
    integer, parameter :: maxr=256
    double precision dlocal(maxr)
!\end{verbatim} %+
    integer ilocal(maxr*nwpr)
    integer i
    equivalence (dlocal,ilocal)
    if(n.gt.256) then
       write(*,*)'LOADRN cannot handle arrays larger than ',maxr
       buperr=1010; goto 900
    endif
    do i=1,n*nwpr
       ilocal(i)=iws(i)
    enddo
    DO I=1,N
       ARR(I)=dlocal(I)
    enddo
900 continue
    RETURN
  END SUBROUTINE LOADRN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine storr1 & Store 1 double at current position
!\begin{verbatim}
  SUBROUTINE STORR1(ARR,VAL)
! store a single double in workspace
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    double precision arr,val
!\end{verbatim} %+
    ARR=VAL
    RETURN
  END SUBROUTINE STORR1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine loadr1 & Load 1 double from current position
!\begin{verbatim}
  SUBROUTINE LOADR1(ARR,VAL)
! load a single double from workspace
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    implicit none
    double precision arr,val
!\end{verbatim}
    VAL=ARR
    RETURN
  END SUBROUTINE LOADR1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
! >>>> subsection
!         2D matrix indexing
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable integer function ixsym & Index 2D array stored as upper triangle 
!\begin{verbatim}
  integer function ixsym(ix1,ix2)
! calculates the storage place of value at (i,j) for a symmetrix matrix
! storage order 11, 12, 22, 13, 23, 33, etc
    implicit none
    integer ix1,ix2
!\end{verbatim} %+
    integer, save :: ncall=0, mcall=0
!    integer, allocatable, dimension(:) :: bug
! at a testing ncall=24623, mcall=507127
!    if(ix1.le.ix2) then
!       mcall=mcall+1
!    else
!       write(*,*)'Indices order',ncall,mcall,ix1,ix2
!       ncall=ncall+1
!    endif
    if(ix1.le.0 .or. ix2.le.0) then
       ixsym=0; buperr=1000; goto 1000
    endif
    if(ix1.gt.ix2) then
       ixsym=ix2+ix1*(ix1-1)/2
    else
       ixsym=ix1+ix2*(ix2-1)/2
    endif
1000 continue
    return
  end function ixsym

!\addtotable integer function kxsym & Index 2D array stored as upper triangle
!\begin{verbatim}
  integer function kxsym(ix1,ix2)
! calculates the storage place of value at (i,j) for a symmetrix matrix
! storage order 11, 12, 22, 13, 23, 33, etc
! In OC the calls to ixsym take about 10 % of the CPU time
! I am trying to replace with local indexing but I need a routine
! that calculates the index when both indices are equal or when I know
! the second index is larger
!    if(ix1.le.0 .or. ix2.le.0) then
!       buperr=1000; goto 1000
!    endif
    implicit none
    integer ix1,ix2
!\end{verbatim}
! this if should be removed when all works
    if(ix1.gt.ix2) stop "Illegal call to kxsym"
    kxsym=ix1+ix2*(ix2-1)/2
    return
  end function kxsym

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
!
!  >>>> subsection
!      miscaleneous
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
  
!\addtotable subroutine fxdflt & Add file extension
!\begin{verbatim}
  subroutine fxdflt(file,ext)
! add default file extention, no good as it thinks .. is an externtion
    implicit none
    character file*(*),ext*(*)
!\end{verbatim} %+
    integer kx
    if(len_trim(file).gt.0) then
       kx=index(file,'.')
       if(kx.le.0) then
          kx=len_trim(file)
          file(kx+1:)='.'//ext
       elseif(kx.lt.len(file)) then
          if(file(kx:kx+1).eq.'..') then
             kx=len_trim(file)
             file(kx+1:)='.'//ext
          endif
       endif
    else
       write(*,*)'No file name'
       file=' '
    endif
    return
  end subroutine fxdflt

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine iniio & Initiate I/O variables
!\begin{verbatim}
  subroutine iniio
! initiates i/o variables, they are all global variables
    implicit none
!\end{verbatim} %+
    kou=koud
    kiu=kiud
    ler=lerd
    iox=0
    return
  end subroutine iniio

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine fisepa & Find separator
!\begin{verbatim}
  SUBROUTINE FISEPA(STR,IP0,IP1)
!...FINDS A SEPARATOR AFTER POSITION IP0
!      A separator is:
!      Any character exept A-Z, 0-9 and _
    implicit none
    CHARACTER STR*(*)
    integer IP0,IP1
!\end{verbatim} %+
    CHARACTER CH1*1
    integer l
    L=LEN_TRIM(STR)
    IP1=IP0
100 IP1=IP1+1
    IF(IP1.GT.L) GOTO 900
    CH1=BIGLET(STR(IP1:IP1))
    IF((LGE(CH1,'0') .AND. LLE(CH1,'9')) .OR. &
         (LGE(CH1,'A') .AND. LLE(CH1,'Z')) .OR. CH1.EQ.'_') GOTO 100
!...Return position before separator
    IP1=IP1-1
900 RETURN
  END SUBROUTINE FISEPA

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine fdmtp & Find matching )
!\begin{verbatim}
  SUBROUTINE FDMTP(LINE1,IP,LINE2)
!...FINDS A MATCHING ) AFTER THAT AT IP. IP UPDATED TO POSITION AFTER )
    implicit none
    CHARACTER LINE1*(*),LINE2*(*)
    integer ip
!\end{verbatim} %+
    integer jp,np,kp1,kp2
    IF(IP.LE.0) GOTO 900
    JP=IP+1
!...np is number of inner levels of parenthesis
    NP=0
10  KP1=INDEX(LINE1(JP:),'(')
    KP2=INDEX(LINE1(JP:),')')
    IF(KP1.EQ.0) THEN
       IF(NP.EQ.0) GOTO 100
       NP=NP-1
       IF(KP2.EQ.0) GOTO 910
       JP=JP+KP2
    ELSEIF(KP1.LT.KP2) THEN
!...INNER PAIR OF ()
       JP=JP+KP1
       NP=NP+1
    ELSEIF(KP1.GT.KP2) THEN
       IF(NP.EQ.0) GOTO 100
       NP=NP-1
       IF(KP2.EQ.0) GOTO 910
       JP=JP+KP2
    ELSE
!       STOP 'FDMTP'
       buperr=1237
       goto 900
    ENDIF
    GOTO 10
!...LINE2 SET TO TEXT INSIDE ( ), IP UPDATED TO POSITION BEHIND )
100 IF(KP2.EQ.0) GOTO 920
    LINE2=LINE1(IP+1:JP+KP2-2)
    IP=JP+KP2
900 RETURN
!910 CALL ST2ERR(1235,'FDMTP','TOO MANY (')
910    buperr=1235
    GOTO 900
!920 CALL ST2ERR(1235,'FDMTP','NO MATCHING )')
920 buperr=1236
    GOTO 900
  END SUBROUTINE FDMTP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable integer function kndex & Find substring from current position
!\begin{verbatim}
  INTEGER FUNCTION KNDEX(LINE,IP,SS)
! SUBROUTINE KNDEX
!...SEARCHES FOR STRING SS IN LINE FROM IP
    implicit none
    CHARACTER LINE*(*),SS*(*)
    integer ip
!\end{verbatim} %+
    integer k
    K=INDEX(LINE(IP:),SS)
    IF(K.GT.0) K=IP-1+K
    KNDEX=K
    RETURN
  END FUNCTION KNDEX

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine cpsstr & Remove tabs and multiple spaces
!\begin{verbatim}
  SUBROUTINE CPSSTR(STRING,LC)
!...THIS SUBROUINE COMPRESSES STRING BY REPLACING MULTIPLE SPACES
!	OR TABS WITH A SINGLE SPACE
    implicit none
    CHARACTER STRING*(*)
    integer LC
!\end{verbatim} %+
    integer i,k,l
    CALL UNTAB(STRING)
10  K=INDEX(STRING(1:LC),'  ')
    IF(K.GT.0) THEN
       L=K
       IF(EOLCH(STRING(1:LC),L)) THEN
          LC=K-1
          GOTO 900
       ENDIF
       L=L-K-1
       LC=LC-L
       REMSP: DO I=K+1,LC
          STRING(I:I)=STRING(I+L:I+L)
       enddo REMSP
       STRING(LC+1:)=' '
       GOTO 10
    ENDIF
900 RETURN
  END SUBROUTINE CPSSTR

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine untab & Remove tab characters
!\begin{verbatim}
  SUBROUTINE UNTAB(LINE)
!...REMOVES ALL TABS FROM LINE. INSERTS SPACES UP TO NEXT TAB STOP
!       TAB STOPS GIVEN IN ITABS. TABS AFTER POSITION 80 REPLACED
!       WITH A SPACE
    implicit none
    CHARACTER LINE*(*)
!\end{verbatim}
    CHARACTER CHTAB*1,XLINE*128
    integer ITABS(11)
    DATA ITABS/8,16,24,32,40,48,56,64,72,80,81/
    integer k,i
!    
    CHTAB=CHAR(9)
100 XLINE=LINE
    K=INDEX(XLINE,CHTAB)
    IF(K.GT.0) THEN
       ADDSP: DO I=1,10
          IF(ITABS(I).GE.K) GOTO 120
       enddo ADDSP
!...BEYOND POSITION 80
       I=11
       ITABS(11)=K
120    I=ITABS(I)
       LINE(K:I)=' '
       LINE(I+1:)=XLINE(K+1:)
       GOTO 100
    ENDIF
    RETURN
  END SUBROUTINE UNTAB

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine incnum
!\begin{verbatim}
  SUBROUTINE INCNUM(NUMB)
!...increments the number in character NUMB by 1, if >9 set to 0
! and increment precedent number 
! if first number >9 set all to zero.
    implicit none
    CHARACTER numb*(*)
!\end{verbatim}
    integer clen,ipos,idig,czero
    czero=ichar('0')
    clen=len(numb)
    ipos=clen
    loop: do while(ipos.gt.0)
       idig=ichar(numb(ipos:ipos))-czero
       if(idig.eq.9) then
          numb(ipos:ipos)='0'
          ipos=ipos-1
       else
          numb(ipos:ipos)=char(czero+idig+1)
          exit loop
       endif
    enddo loop
1000 continue
    return
  end SUBROUTINE INCNUM

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

END MODULE METLIB
