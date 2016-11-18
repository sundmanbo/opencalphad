!
! general utilities in Fortran 95 a la METLIB
!
MODULE METLIB
!
! Copyright 1980-2015, Bo Sundman and others, bo.sundman@gmail.com 
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
!  1030  NO SUCH TYPE OPTION
!  1031 Empty line, expected number
!  1032  PARAMETER VALUE MISSING
!  1033 Decimal point but no digits
!  1034 No digits
!  1035 Positive sign but no digits
!  1036 Negative sign but no digits
!  1037 No sign and no digits
!  1038  NO DIGITS AFTER EXPONENTIAL E
!  1039 Exponent larget then 99
!  1040  NO HELP FOR <COMMAND>
!  1041  NO SUCH QUESTION FOR <COMMAND>
!  1042  TOO LARGE INTEGER VALUE
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
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! default units for command output, input, error message, list
! and default language
  PARAMETER (koud=6,kiud=5,keud=6,lutd=6,lerd=6,langd=1)
! representation of the numerical value "none"
  character cnone*4
  PARAMETER (RNONE=-1.0D-36,NONE=-2147483647,CNONE='NONE')
  PARAMETER (MAXINT=2147483647,MININT=-2147483646)
  PARAMETER (FLTSML=1.0D-36,FLTPRS=1.0D-14)
! initiate i/o and error code
  integer :: kou=koud,kiu=kiud,keu=keud,lut=lutd,ler=lerd,lang=langd
  integer :: buperr=0,iox(10)
! LSTCMD is the last command given. Saved by NCOMP, used by help routines
  character, private :: lstcmd*40
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
! all COMMON removed
! This is for environment variables used in MACROs
  character, private :: ENVIR(9)*60
!
!    COMMON/TCMACRO/IUL,IUN(5),MACEXT
    integer, private :: IUL,IUN(5)
    character MACEXT*3
! nbpw is number if bytes per INTEGER
    integer, parameter :: nbpw=4,nwpr=2
    parameter (MACEXT='OCM')
!    parameter (MACEXT='BMM')
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    integer, parameter :: maxhelplevel=15
! A help structure used in new on-line help system
    TYPE help_str
       integer :: okinit=0
       character*128 filename
       character*8 type
       integer level
       character*32, dimension(maxhelplevel) :: cpath
    END TYPE help_str
! this record is used to file the appropriate help text
    type(help_str), save :: helprec
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! Data structures for putfun below

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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CONTAINS

! subroutines and functions
! UCLETTER TRUE if letter is letter is apper case A-Z
! CAPSON   coverts a text to upper case
! SORTRD   sorts an array or double reals in ascending order
! SORTDD   sorts an array or double reals in ascending order
! SORTIN   sorts an array of integers
! GETREL   extract a double number from a text
! GETRELS  extract a double number from a text
! GPS      extract a signed number
! GPN      extract an unsigned number
! WRINUM   writes a double into a text
! WRIINT   writes an integer into a text
! BIGLET   converts a letter to upper case
! EOLCH    TRUE if text is empty
! GETNAME  extracts a name with characters
! IXSYM    index of a symmetric array stored as a single array
! WRICE    writes a text over several lines
! WRICE2   writes a text over several lines
! CWRICEND finds the place to insert newlines in a text
! NGPSBS  Substitutes symbols with values
! NGHELP  Standard HELP command
! NGOPFI  writes an error message
! FILDOC  Routine to read documentation of a command from file
! FILHLP  Routine to read help text for a question from a file
! FILINF  Routine for the INFORMATION command
! NGGOTO  Standard check of a GOTO command
! GPAROP  Input of options
! NXTRCT  Checks if the next nonblank character is one of a given set
! NCOMP   Command interpreter accepting abbreviations
! NGHIST  Executes a history command
! LABBR   Comparares an abbreviation to a full name
! GPARQ   Routine to give information before prompting
! GPARLD  Routine to ask for multi-line input with default value
! GPARL   (extry) Routine to ask for multi-line input
! GQARID  Same as GPARI but default value is given also.
! GQARI   Routine to ask for an integer value from a string.
! GQARRD  Same as GPARR but default value is given also.
! GQARR   Routine to ask for a real     value from a string.
! GQARCD  Same as GPARC but default value is given also.
! GQARC   Routine to ask for a string   value from a string.
! SET_ECHO set echo of commands
! GETEXT  Extracts a text from a character string.
! GETINM  Decodes an integer from a string and skipps a final ,.
! GETREM  Decodes a real     from a string and skipps a final ,.
! GETINT  Decodes an integer from a string.
! GETOCT  Decodes an octal value without sign from a string
! GETHEX  Decodes a hexadecimal value without sign from a string
! NOHELP  Default help routine.
! TOPHLP  Help routine that forces return to calling program on ?
! FISEPA  Finds the first character not a digit, letter or _
! WRIHEX  Edits a hexadecimal number into a string.
! SSORT   SORTS A CHARACTER ARRAY WITHOUT MOVING DATA
! FDMTP   Extracts a sequence surrounded by parenthesis
! YESCHK  Returns TRUE if argument is Y or y
! KNDEX   As INDEX but starts from a given position
! CPSSTR  Compress tabs and multiple spaces in a string to a single space
! UNTAB   Replace all tab characters in a string with spaces
! GPARID calls GQARID after replacing environment variables
! GPARI  calls GQARI after replacing environment variables
! GPARR  calls GQARR after replacing environment variables
! GPARRD calls GQARRD
! GPARC  calls GQARC
! GPARCD  calls GQARCD
! openlogfile opens a log file 
! GPTCM1 handles macro directives
! GQEXENV exchanges macro variables with actual values
! MACBEG starts a macro
! TESTB, setb clrb obsolete, use btest etc
! BOUTXT writes a promt
! BINTXT read a command line
! IONOFF sets ionoff to zero
! COMND  stop
! TCQFFB stop
! FXDFLT add a file extention if none
! GPARFD stop
! INIIO  sets default io routine units
! PUTFUN enters a function as a binary tree
! NYBIN  creates a binary bode in symbol tree
! NYUNI  creates a unary node
! NYLP   handles an opening parenthesis in a function
! NYRP   handles a closing parenthesis in a function
! NYVAR  handles a symbol in a function
! NYDAT  insets a data node
! EVALF  evaluates a function stored as a binary tree
! EUNARY evaluates a unary node
! EBINRY evaluates a binary node
! AIVAN  Evaluates IVANTSSOV's solution of error function
! PF_BSUM calculates a sum
! PF_HS  calculates Heaviside function
! PF_ERF calculates the error function
! WRTFUN writes a function stored as a binary tree
! WRTLPQ writes an opening parenthesis
! WRTRPQ writes a closing parenthesis
! WRTBIQ writes a binary operator
! WRTDAQ writes a numeric value
! DELFUN deletes a function (not implemented??)
! CONS   concatinates two characters (unused??)
! GZRFUN interactive input of a function
! EXPHLP help text for entering function
! PUTPRP creates a prompt for a function including variables
! HPCALC HP calculator
! HPHLP  some help using the HP calculator
! init_help initiates help for OC
! helplevel1 reset helplevel to 1
! q2help  help for submenus
! q1help  help for questions
! q3help  help in submenues
! winit   initates integer workspace
! wold    read a file to workspace
! wsave   write workspace to file
! wpatch  interactive patching a workspace
! wphelp  help for wpatch
! wrkchk  check integrety of workspace
! wlist   list free areas
! wtrest  reserves the rest of a workspace
! wtake   reserves a record in workspace
! wrels   returns a record to workspace
! nwch    number of words needed to store a character
! storc   stores a character in workspace
! loadc   extracts a character from workspace
! storr   stores a double in workspace
! loadr   extracts a double from workspace
! storrn  stores an array of doubles in workspace
! loadrn  extract an array of doubles from workspace
! storr1  stores a double in workspace
! loadr1  extracts a double from workspace
! 
  LOGICAL FUNCTION ucletter(ch1)
! returns TRUE if the character is A to Z
    character ch1*1
    if(lge(ch1,'A') .and. lle(ch1,'Z')) then
       ucletter=.TRUE.
    else
       ucletter=.FALSE.
    endif
  END FUNCTION ucletter
  
  SUBROUTINE capson(text)
! changes lower case ASCII a-z to upper case A-Z, no other changes
    character text*(*)
    parameter (lowa=ichar('a'),lowz=ichar('z'),iup=ICHAR('A')-ICHAR('a'))
    DO i=1,len(text)
       ich1=ichar(text(i:i))
       IF(ich1.ge.lowa .and. ich1.le.lowz) THEN
          text(i:i)=char(ich1+iup)
       ENDIF
    ENDDO
  END SUBROUTINE capson
  
! intrinsic bitroutines
! LOGICAL FUNCTION BTEST(iword,ibit) bit numbered from right
! iword2=IBSET(iword1,ibit) ettar ibit of iword1 and returns in iword2
! iword2=IBCLR(iword1,ibit) nollar ibit of iword1 and returns in iword2
  
!...SORTING SUBROUTINES FOR INTEGERS, REALS AND CHARACTERS
!      QUICKSORT ENL KNUTH ALGORTIM Q
!      THE ART OF COMPUTER PROGRAMMING, VOL 3, P 117
  
  SUBROUTINE SORTRD(ARR,N,IX)
! ...SORTING REAL NUMBERS IN ASCENDING ORDER
! ENTRY:
!      ARR   ARRAY TO BE SORTED
!      N     NUMBER OF ELEMENTS TO BE SORTED >1
!      IX    INTEGER ARRAY WITH DIMENSION N
! EXIT:
!      ARR   SORTED ARRAY
!      IX    ARRAY WHERE IX(I) IS THE PREVIOS INDEX OF ARR(I)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (MSTACK=20)
    DIMENSION ARR(*),IX(*),LOW(MSTACK),IGH(MSTACK)
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
210 VAL=ARR(I)
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

  SUBROUTINE SORTRDD(ARR,N,IX)
! ...SORTING REAL NUMBERS IN DECENDING ORDER
! ENTRY:
!      ARR   ARRAY TO BE SORTED
!      N     NUMBER OF ELEMENTS TO BE SORTED >1
!      IX    INTEGER ARRAY WITH DIMENSION N
! EXIT:
!      ARR   SORTED ARRAY
!      IX    ARRAY WHERE IX(I) IS THE PREVIOS INDEX OF ARR(I)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (MSTACK=20)
    DIMENSION ARR(*),IX(*),LOW(MSTACK),IGH(MSTACK)
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
210 VAL=ARR(I)
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

  SUBROUTINE SORTIN(IARR,N,IX)
! ...SORTING INTEGERS IN ASCENDING ORDER
! ENTRY:
!      IARR   ARRAY TO BE SORTED
!      N     NUMBER OF ELEMENTS TO BE SORTED >1
!      IX    INTEGER ARRAY WITH DIMENSION N
! EXIT:
!      IARR   SORTED ARRAY
!      IX    ARRAY WHERE IX(I) IS THE PREVIOS INDEX OF IARR(I)
    PARAMETER (MSTACK=20)
    DIMENSION IARR(*),IX(*),LOW(MSTACK),IGH(MSTACK)
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
210 IVAL=IARR(I)
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
  
  SUBROUTINE GETREL(SVAR,LAST,VALUE)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    character svar*(*)
    call getrels(svar,last,value,isig)
    return
  END SUBROUTINE GETREL

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
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! EOLCH is declared as logical function in this module so not needed here
!    LOGICAL EOLCH
    PARAMETER (ZERO=0.0D0,ONE=1.0D0,TEN=1.0D1)
    CHARACTER SVAR*(*),CH*1
!    INTEGER GPS,GPN
3   CONTINUE
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
10  continue
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
20  continue
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
30  EXPO=ONE
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
800 VALUE=DBLE(ISIG)*(HEL+DEC)*EXPO
9000 continue
    RETURN
  END SUBROUTINE GETRELS

  INTEGER FUNCTION GPS(SVAR,LAST,VALUE)
!...DECODES A NUMBER WITH OR WITHOUT A SIGN
    DOUBLE PRECISION VALUE
    CHARACTER SVAR*(*),SIG*1
!    INTEGER GPN
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
    VALUE=FLOAT(ISIG)*VALUE
    RETURN
  END FUNCTION GPS

  INTEGER FUNCTION GPN(SVAR,LAST,VALUE)
!...DECODES A NUMBER WITOUT SIGN
    DOUBLE PRECISION ZERO,TEN,VALUE
    PARAMETER (ZERO=0.0D0,TEN=1.0D1)
    CHARACTER SVAR*(*)
    L=LEN_TRIM(SVAR)
    VALUE=ZERO
    IERR=1034
    digits: DO LAST=LAST,L
       N=ICHAR(SVAR(LAST:LAST))-ICHAR('0')
       IF(N.LT.0 .OR. N.GT.9) GOTO 800
       IERR=0
       VALUE=TEN*VALUE+FLOAT(N)
    enddo digits
800 GPN=IERR
    RETURN
  END FUNCTION GPN
  
  SUBROUTINE WRINUM(STR,IP,NNW,JSIGN,VALUE)
!...EDITS A REAL NUMBER INTO STR WITH LEAST NUMBER OF DIGITS
!      NNW IS MAXIMUM NUMBER OF SIGNIFICANT DIGITS (0<NNW<16)
!      JSIGN >0 INDICATES THAT + SIGN SHOULD BE WRITTEN
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER STR*(*),CSTR*21,CFRMT*12
    PARAMETER (ZERO=0.0D0,TEN=1.0D1,EPS=1.0D-7)
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
27  format('wrinum: ',2i5,1pe20.12)
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

  subroutine wriint(text,ipos,int)
! write an integer in text from position ipos (left adjusted)
    implicit none
    character text*(*),number*16
    integer ipos,int,jp
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
22        format(a,i3,'>',a,'< >',a,'<')
          ipos=ipos+17-jp
!          write(*,30)'wriint: ',ipos,jp,' >'//text(1:ipos+5)//'<'
30        format(a,2i3,a)
       endif
    endif
1000 continue
    return
  end subroutine wriint

  CHARACTER FUNCTION BIGLET(CHA)
!...CONVERTS ONE CHARACTER FROM LOWER TO UPPER CASE
    CHARACTER*1 CHA,CHLAST
    PARAMETER (CHLAST='z')
    IF(CHA.GE.'a' .AND. CHA.LE.CHLAST) THEN
       BIGLET=CHAR(ICHAR(CHA)+ICHAR('A')-ICHAR('a'))
    ELSE
       BIGLET=CHA
    ENDIF
900 RETURN
  END FUNCTION BIGLET
  
  LOGICAL FUNCTION EOLCH(STR,IP)
!...TO SKIP SPACES FROM IP. RETURNS .TRUE. IF ONLY SPACES
!....MODIFIED TO SKIP TAB CHARACTERS ALSO
    CHARACTER STR*(*)
    PARAMETER (ITAB=9)
    EOLCH=.FALSE.
    IF(IP.LE.0) IP=1
100 IF(IP.GT.LEN(STR)) GOTO 110
    IF(STR(IP:IP).NE.' ' .AND. ICHAR(STR(IP:IP)).NE.ITAB) GOTO 900
    IP=IP+1
    GOTO 100
110 EOLCH=.TRUE.
900 RETURN
  END FUNCTION EOLCH

  subroutine getname(text,ip,name,mode,ch1)
! this should be incorporated in metlib
    character text*(*),name*(*),ch1*1
! Always a letter A-Z as first character
! mode=0 is normal, letters, numbers, "." and "_" allowed
! mode=1 used for species names with "/", "+" and "-" allowed also
    ch1=biglet(text(ip:ip))
    if(ch1.lt.'A' .or. ch1.gt.'Z') then
       write(*,17)ichar(ch1),ch1,text(1:24),ip
17     format('getname error: ',i5,' "',a,'" in "',a,'" at ',i4)
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

  integer function ixsym(ix1,ix2)
! calculates the storage place of value at (i,j) for a symmetrix matrix
! storage order 11, 12, 22, 13, 23, 33, etc
    if(ix1.le.0 .or. ix2.le.0) then
       bmperr=1000; goto 1000
    endif
    if(ix1.gt.ix2) then
       ixsym=ix2+ix1*(ix1-1)/2
    else
       ixsym=ix1+ix2*(ix2-1)/2
    endif
1000 continue
    return
  end function ixsym

  subroutine wrice(lut,margl1,margl2,maxl,str)
! writes str on unit lut with left margin largl1 for first line, margl2 for all
! following lines, max length maxl characters (assuming typewriter font)
    character str*(*),margx*40
    lbreak=0
    call wrice2(lut,margl1,margl2,maxl,lbreak,str)
1000 continue
    return
  end subroutine wrice

  subroutine wrice2(lut,margl1,margl2,maxl,lbreak,str)
! writes str on unit lut with left margin largl1 for first line, margl2 for all
! following lines, max length maxl characters (assuming typewriter font)
! lbreak>0 for writing math expression, with stricter linebreak rules
! lbreak<0 for breaking only at space
    character str*(*),margx*40
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
20 format('wrice :',10I5)
    return
  end subroutine wrice2

  subroutine cwricend(str,lbeg,lend,lbreak)
! find a possible place for a newline in str going back from lend
! but not bypassing lbeg.  str is a numerical expression.
! lbreak>0 means stricter rules (mathematical expression)
! lbreak<0 means break only at space
    character str*(*),ch1*1,ch2*1
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

!=============================================================
! start on-line help package
!=============================================================
  SUBROUTINE NGPSBS(LINE,CHT,NS,COLUMN,LENC)
!...SUBSTITUTES SYMBOLS WITH VALUES
! ENTRY: LINE   TEXT LINE WITH POSSIBLE PARAMETER REFERENCES
!        CHT    PARAMETER REFERENCE CHARACTER
!        NS     NUMBER OF ENTRIES IN COLUMN (NUMBER OF PARAMETER VALUES)
!        COLUMN CHARACTER ARRAY WITH PARAMETER VALUES
!        LENC   INTEGER ARRAY WITH LENGHT OF VALUE TEXT
! EXIT:  LINE   REFERENCES REPLACED WITH VALUES
!        ALL OTHER ARGUMENTS UNCHANGED
    CHARACTER LINE*(*),CHT*1,COLUMN(*)*(*)
    CHARACTER LLINE*132
    DIMENSION LENC(*)
!...PERFORM SUBSTITUTION
    LLINE=LINE
200 IA=1
    JA=1
300 IB=KNDEX(LLINE,IA,CHT)
    IF(IB.LE.0) GOTO 800
    IC=IB+1
    CALL GETINT(LLINE,IC,KOLUMN)
    IF(buperr.ne.0) then
       buperr=1360; goto 900
    endif
    LINE(JA:)=LLINE(IA:IB-1)//COLUMN(KOLUMN)
    JA=JA+IB-IA+LENC(KOLUMN)
    IA=IC
    GOTO 300
!...COPY REST OF LINE
800 LINE(JA:)=LLINE(IA:)
900 RETURN
  END SUBROUTINE NGPSBS
! ======
  SUBROUTINE NGHELP(LINE,LAST,COMM,NC)
!...EXECUTES A HELP COMMAND
    CHARACTER LINE*(*),COMM(NC)*(*),CMD*40
    PARAMETER (MC=100)
    DIMENSION INDX(MC)
!    EXTERNAL FILHLP
    CALL GPARC('COMMAND: ',LINE,LAST,1,CMD,'*',FILHLP)
    IF(CMD(1:1).EQ.'*') THEN
!...LIST ALL COMMANDS IN UNIX ALPHABETICAL ORDER
       NKPL=80/(LEN(COMM(1))+1)
       IF(NKPL*(LEN(COMM(1))+1).GE.80) NKPL=NKPL-1
100    IF(NC.LT.MC) THEN
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
          CALL CAPSON(CMD)
          CALL FILDOC(CMD)
       ELSEIF(K.EQ.0 .OR. K.LT.-NC) THEN
          WRITE(KOU,*)'NO MATCHING COMMAND, USE HELP *'
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
  END SUBROUTINE NGHELP
! ====================
  SUBROUTINE NGOPFI(LUNDOC)
    write(*,*)'No help file available, sorry'
!    STOP 'NGOPFI'
    RETURN
  END SUBROUTINE NGOPFI
! ====================
  SUBROUTINE FILDOC(STRING)
!...Reads documentation of commands from a file to the terminal
    CHARACTER STRING*(*)
    CHARACTER LINE*80,CH1*1,PAUSE*2
    LOGICAL EXACT
!    EXTERNAL NOHELP
    CALL NGOPFI(LUNDOC)
! FIRST LINE GIVES NUMBER OF COMMANDS
    READ(LUNDOC,99,REC=1,ERR=910,IOSTAT=IERR)LINE
99  FORMAT(A)
    IP=1
    CALL GETINT(LINE,IP,NRC)
    K=LEN_TRIM(LINE)
    PAUSE=LINE(K-1:K)
    K=LEN_TRIM(STRING)
    LOCC=2
    ALLCOM: DO MRC=1,NRC
!...      LOOP TO FIND MATCHING COMMAND LINE
! COMMAND(1:40)  #QUESTIONS  RECORD WITH FIRST INFO  #INFOLINES
       READ(LUNDOC,99,REC=LOCC,ERR=910,IOSTAT=IERR)LINE
       IP=41
       CALL GETINT(LINE,IP,NRQ)
       IF(LABBR(STRING(1:K),LINE(1:40),EXACT)) THEN
!...MATCHING COMMAND, READ THE NEXT TWO NUMBERS
          CALL GETINT(LINE,IP,LOCHLP)
          CALL GETINT(LINE,IP,NRL)
          WRITE(KOU,*)
          RLINE: DO JR=LOCHLP,LOCHLP+NRL-1
             KK=IONOFF('FILDOC')
             IF(KK.EQ.15) GOTO 60
             READ(LUNDOC,99,REC=JR,ERR=910,IOSTAT=IERR)LINE
             IP=LEN_TRIM(LINE)
             WRITE(KOU,40)LINE(1:MAX(1,IP))
40           FORMAT(1X,A)
50           CONTINUE
          enddo RLINE
60        WRITE(KOU,*)
       ENDIF
       LOCC=LOCC+NRQ+1
    enddo ALLCOM
900 RETURN
!910 CALL ST2ERR(IERR,'FILDOC','FILE SYSTEM ERROR')
910 buperr=1350
    GOTO 900
  END SUBROUTINE FILDOC
! ====================
  SUBROUTINE FILHLP(PROMPT,SVAR)
!...TO OBTAIN HELP FROM A FILE WHEN TYPING ?. EXTRA HELP WHEN ??
    CHARACTER PROMPT*(*),SVAR*(*)
    CHARACTER LINE*80,CH1*1,PROMT*40
! NOTE! THE FILE MUST BE WRITTEN WITHOUT PARITY!!!!!!!!!!!!!!!!!!!!!!
! FILE STRUCTURE:
! BLOCKED WITH 80 CHARACTERS PER RECORD.
! FIRST LINE HOLDS THE NUMBER OF COMMANDS
! THE SECOND LINE IS A COMMAND RECORD
! COMMAND IN (1:40), THEN THREE NUMBERS: #QUESTIONS, INFO RECORD #, #INFO LIN
!      THE NEXT COMMAND IS IN RECORD = CURRENT + #QUESTIONS + 1
!      THE QUESTIONS FOLLOWS DIRECTLY AFTER THE COMMAND RECORD
!      ONE QUESTION PER RECORD, INDENTED WITH 3 SPACES
! QUESTION(4:60), THEN 2 OR 4 NUMBERS:
!                 HELP RECORD#, #LINES, EXTRA RECORD#, #LINES
    CALL NGOPFI(LUNDOC)
    READ(LUNDOC,99,REC=1,ERR=910,IOSTAT=IERR)LINE
99  FORMAT(A)
    IP=1
    CALL GETINT(LINE,IP,NRC)
    IF(buperr.ne.0) GOTO 950
    PROMT=PROMPT
    CALL CAPSON(PROMT)
!...FIRST RECORD WITH QUESTION IS NRC+2
    LOCC=2
    FINDCOM: DO I=1,NRC
       READ(LUNDOC,99,REC=LOCC,ERR=910,IOSTAT=IERR)LINE
       IP=41
       CALL GETINT(LINE,IP,NRQ)
       IF(buperr.ne.0) GOTO 950
       CALL CAPSON(LINE)
       IF(LSTCMD.EQ.LINE(1:40)) GOTO 200
       LOCC=LOCC+NRQ+1
    enddo FINDCOM
!...THE COMMAND WAS NOT FOUND
    GOTO 920
!...MATCHING COMMAND FOUND
200 IF(NRQ.EQ.0) GOTO 930
    FINDQ: DO I=LOCC+1,LOCC+NRQ
!...      SEARCH FOR A QUESTION MATCHING PROMT. USE SHORTEST STRING LENGTH!
       READ(LUNDOC,99,REC=I,ERR=910,IOSTAT=IERR)LINE
       J=MIN(LEN_TRIM(LINE(4:40)),LEN_TRIM(PROMT))
       CALL CAPSON(LINE(4:J+3))
       IF(LINE(4:J+3).EQ.PROMT(1:J)) GOTO 300
    enddo FINDQ
    GOTO 930
!...MATCHING QUESTION FOUND
300 IP=51
    CALL GETINT(LINE,IP,LOCHLP)
    CALL GETINT(LINE,IP,NRL)
    if (LEN_TRIM(svar).ge.2) then
       IF(SVAR(1:2).EQ.'??') THEN
          I=LOCHLP
          CALL GETINT(LINE,IP,LOCHLP)
          IF(BUPERR.NE.0) THEN
             BUPERR=0
             LOCHLP=I
          ELSE
             CALL GETINT(LINE,IP,NRL)
          ENDIF
       ENDIF
    endif
    IF(LOCHLP.LE.0 .OR. BUPERR.NE.0) GOTO 940
    WRITE(KOU,*)
    LISTH: DO I=LOCHLP,LOCHLP+NRL-1
!...      LIST ON SCREEN, CANNOT USE GPARC HERE BECAUSE NESTED CALLS!
       KK=IONOFF('FILHLP')
       IF(KK.EQ.15) GOTO 360
       READ(LUNDOC,99,REC=I,ERR=910,IOSTAT=IERR)LINE
!         LINE(79:80)='  '
       J=MAX(1,LEN_TRIM(LINE))
       WRITE(KOU,310)LINE(1:J)
310    FORMAT(1X,A)
    enddo LISTH
360 WRITE(KOU,*)
900 IF(buperr.ne.0) THEN
       IF(LER.GT.0) WRITE(LER,*)'ERROR ',buperr,' READING HELP FILE'
    ENDIF
901 RETURN
910 IF(LER.GT.0) WRITE(LER,*)'FILE SYSTEM ERROR IN FILHLP'
    GOTO 900
920 IF(LER.GT.0) WRITE(LER,*)'NO HELP FOR ',LSTCMD
    GOTO 900
930 IF(LER.GT.0) WRITE(LER,*)'UNKNOWN QUESTION ',PROMPT
    GOTO 900
940 IF(LER.GT.0) WRITE(LER,*)'NO EXPLANATION AVAILABLE, LINE: ',IP
    BUPERR=0
    GOTO 900
950 IF(LER.GT.0) WRITE(LER,*)'Error occurred on line ',IP
    GOTO 900
  END SUBROUTINE FILHLP
! ====================
  SUBROUTINE FILINF(ILINE,LAST)
!...EXECUTES AN INFORMATION COMMAND
    CHARACTER ILINE*(*),LINE*80,SUBJCT*40,SUBDEF*7,CH1*1
!    CHARACTER BIGLET*1
!    EXTERNAL NOHELP
!    EXTERNAL FILHLP
!    LOGICAL SG2ERR,LABBR,EXACT
    LOGICAL EXACT
    CALL NGOPFI(LUNDOC)
    READ(LUNDOC,99,REC=1,ERR=910,IOSTAT=IERR)LINE
99  FORMAT(A)
    IP=1
    CALL GETINT(LINE,IP,NRC)
    IF(BUPERR.NE.0) GOTO 900
!...SCAN THE COMMANDS
    LOCC=2
    ALLCOM: DO I=1,NRC
       READ(LUNDOC,99,REC=LOCC,ERR=910,IOSTAT=IERR)LINE
       IP=41
       CALL GETINT(LINE,IP,NRQ)
       IF(BUPERR.NE.0) GOTO 900
       IF(LSTCMD.EQ.LINE(1:40)) GOTO 200
       LOCC=LOCC+NRQ+1
    enddo ALLCOM
!...THE COMMAND WAS NOT FOUND
    GOTO 920
!...MATCHING COMMAND FOUND
200 IF(NRQ.EQ.0) GOTO 920
    SUBDEF='PURPOSE'
210 CALL GPARCD('WHICH SUBJECT ',ILINE,LAST,1,SUBJCT,SUBDEF,FILHLP)
    IF(SUBJCT(1:1).EQ.' ') GOTO 900
    SUBDEF=' '
    CALL CAPSON(SUBJCT)
!...SCAN THE QUESTIONS
    LOCS=LOCC
    ALLINF: DO NR=1,NRQ
       LOCS=LOCS+1
       READ(LUNDOC,99,REC=LOCS,ERR=910,IOSTAT=IERR)LINE
       IF(LABBR(SUBJCT,LINE(4:60),EXACT)) THEN
          IP=61
          CALL GETINT(LINE,IP,LOCL)
          CALL GETINT(LINE,IP,NRL)
          ALLIN: DO I=1,NRL
             KKB=IONOFF('FILINF')
             IF(KKB.EQ.15) GOTO 440
             READ(LUNDOC,99,REC=LOCL,ERR=910,IOSTAT=IERR)LINE
             IP=MAX(LEN_TRIM(LINE),1)
             WRITE(KOU,310)LINE(1:IP)
310          FORMAT(1X,A)
             LOCL=LOCL+1
          enddo ALLIN
          WRITE(KOU,*)
       ENDIF
    enddo ALLINF
440 LAST=LEN(ILINE)
    GOTO 210
900 RETURN
910 IF(LER.GT.0) WRITE(LER,*)'FILE SYSTEM ERROR ',IERR
    GOTO 900
920 IF(LER.GT.0) WRITE(LER,*)'NO INFORMATION AVAILABLE'
    GOTO 900
  END SUBROUTINE FILINF
!========================================================
! end on-line help package
! =======================================================
! routine to execute a GOTO_MODULE command
!========================================================
  SUBROUTINE NGGOTO(LINE,LAST,NPROG,PROGS,JUMP)
! Compares input with given module names. Used in all module monitors.
! RETURNS JUMP>0 if legal program name given, otherwise JUMP<=0
    CHARACTER LINE*(*),PROG*40,PROGS(NPROG)*(*)
    LOGICAL ONCE
!    EXTERNAL TOPHLP
    ONCE=.TRUE.
100 CALL GPARC('MODULE NAME: ',LINE,LAST,1,PROG,' ',TOPHLP)
    IF(PROG(1:1).NE.' ') THEN
       JUMP=NCOMP2(PROG,PROGS,NPROG,L)
    ELSE
       JUMP=0
    ENDIF
    IF(JUMP.LE.0 .AND. ONCE) THEN
       ONCE=.FALSE.
       WRITE(KOU,110)
110    FORMAT(' NO SUCH MODULE, USE ANY OF THESE:')
       ALLMOD: DO I=1,NPROG
          IF(PROGS(I)(1:1).NE.' ') WRITE(KOU,111)PROGS(I)
111       FORMAT(1X,A)
       enddo ALLMOD
       GOTO 100
    ENDIF
900 RETURN
  END SUBROUTINE NGGOTO
!===================================================
! unused?? routine to read an option 
!===================================================
  SUBROUTINE GPAROP(PROMPT,LINE,LAST,NOP,OPTEXT,OPTVAL,OPTDEF,DV,HELP)
!...OBTAINING OPTIONS
!       THIS SUBROUTINE SHOULD SIMULATE A FACILITY SIMILAR TO THAT
!       USED IN MULTICS AND UNIX. AN EXAMPLE OF OPTIONS THAT COULD
!       BE USEFUL AFTER A LIST COMMAND:
!       Option: ?
!       List of options:
!       Output_to:<file name>           o:line-printer
!       Append_output_to:<file>         a:summary
!       Short                           s
!       Alphabetical_order              al
!       Note that the : is mandatory to separate an option specifier
!       from its value. Many options have no value but are enough in
!       themselvs.
!       Note that a = is equivalent to a : and the = is useful if the
!       value is numeric e.g. temp=1000
!       THE SAME SUBROUTINE SHOULD ALSO BE USEFUL WHEN A COMMAND
!       NEEDS A SPECIFICATION, E.G. IF THERE IS ONLY ONE LIST
!       COMMAND ONE MUST SPECIFY WHAT TO LIST (THIS IS FOR AN OS):
!       List what? ?
!       Files                           f
!       Users                           u
!       System                          s
!       Devices                         d
!       Queues                          q
    CHARACTER PROMPT*(*),LINE*(*)
    CHARACTER OPTEXT(*)*(*),OPTVAL(*)*(*),OPTDEF(*)*(*),CH1(2)*1
! IN OPTEXT THE SPECIFYING TEXT FOR EACH OPTION IS CONTAINED E.G.
!       Files, Output_to:  Max length is 40 characters
! IN OPTVAL THE VALUE OF THE OPTION WILL BE RETURNED
!       Note that numeric values are legal but are stored as texts also
!       Max length is 40 characters
! IN OPTDEF THE DEFAULT VALUE OF THE OPTION IS STORED
!    LOGICAL EOLCH,ONCE,DV
    LOGICAL ONCE,DV
    DATA CH1/'=',':'/
    ONCE=.FALSE.
    DV=.TRUE.
!...SET DEFAULT VALUES ON ALL OPTIONS
    INITVAL: DO I=1,NOP
       OPTVAL(I)=OPTDEF(I)
    enddo INITVAL
!...IF LINE EMPTY WRITE PROMPT ONCE
100 IF(.NOT.EOLCH(LINE,LAST)) GOTO 200
120 IF(PROMPT(1:1).NE.' ') CALL BOUTXT(KOU,PROMPT)
    CALL BINTXT(KIU,LINE)
    LAST=1
    IF(EOLCH(LINE,LAST)) GOTO 900
!...A ? will give help, first only a list of the options.
200 IF(LINE(LAST:LAST).EQ.'?') THEN
!...      Two ?? will give extended help from help routine.
       IF(LINE(LAST+1:LAST+1).EQ.'?') CALL HELP(PROMPT,LINE)
       COMPARAE: DO I=1,NOP
          WRITE(KOU,220)OPTEXT(I),OPTVAL(I)
       enddo COMPARAE
220    FORMAT(1X,A,2X,A)
       WRITE(KOU,*)
       GOTO 120
    ENDIF
300 IF(EOLCH(LINE,LAST)) GOTO 900
!...LAST IS INCREMENTED BY 1 IN NCOMP3
    LAST=LAST-1
    K=NCOMP3(LINE,OPTEXT,NOP,LAST)
    IF(K.LE.0) GOTO 905
    DV=.FALSE.
    I1=INDEX(OPTEXT(K),CH1(1))
    I2=INDEX(OPTEXT(K),CH1(2))
    IF(I1.GT.0 .OR. I2.GT.0) THEN
!...       THIS OPTION REQUIRES A VALUE AFTER A : OR =
       I1=NXTRCT(LINE,LAST,2,CH1)
       IF(I1.EQ.0) GOTO 910
       CALL GETEXT(LINE,LAST,1,OPTVAL(K),OPTDEF(K),LENC)
    ELSE
!...       THIS OPTION REQUIRES NO VALUE BUT A MARK
       OPTVAL(K)='*'
    ENDIF
    GOTO 300
900 RETURN
!905 CALL ST1ERR(1332,'GPAROP','ILLEGAL OPTION')
905 buperr=1332
    GOTO 900
!910 CALL ST1ERR(1333,'GPAROP','MISSING OPTION VALUE DELIMITER')
910 buperr=1333
    GOTO 900
!920 CALL ST1ERR(1334,'GPAROP','MISSING OPTION VALUE')
    buperr=1334
    GOTO 900
  END SUBROUTINE GPAROP
! ====================
  FUNCTION NXTRCT(LINE,LAST,NCH,CHL)
!....SKIPS SPACES AND THEN CHECKS IF THE CHARACTER IS ONE OF THOSE IN CHL
    CHARACTER LINE*(*),CHL(*)*1
!    LOGICAL EOLCH
    IF(EOLCH(LINE,LAST)) THEN
       NXTRCT=0
    ELSE
       COMPARAE: DO I=1,NCH
          IF(LINE(LAST:LAST).EQ.CHL(I)) GOTO 60
       enddo COMPARAE
       I=0
60     NXTRCT=I
    ENDIF
    RETURN
  END FUNCTION NXTRCT
!============================================================
! command interpreter
  !==========================================================
  INTEGER FUNCTION NCOMP(SVAR,COMM,NC,NEXT)
! SUBROUTINE NCOMP
! ...TO DECODE A COMMAND
    CHARACTER SVAR*(*),COMM(NC)*(*),LINE*80
    CHARACTER*1 CH1,CHSEP,CHLAST,CHSYS,CHSEP2,CH2,CHHELP
    CHARACTER*1 CHHIST,CHMAC
!    LOGICAL EXAKT,EOLCH,YESLOG
    LOGICAL EXAKT,YESLOG
!    CHARACTER LSTCMD*40
    PARAMETER (CHSEP='_',CHSEP2='-',CHLAST='Z',CHSYS='@',CHHELP='?')
    PARAMETER (CHHIST='!',CHMAC='#')
    IENT=1
    YESLOG=.TRUE.
! Note the command line should be read by GPARC
    GOTO 2
!...ENTRY SAME AS NCOMP BUT THE COMMAND NOT SAVED
    ENTRY NCOMP2(SVAR,COMM,NC,NEXT)
    IENT=2
! LAST IS SET TO FIRST NONBLANK CHARACTER
2   LAST=1
    GOTO 3
!...ENTRY TO CONTINUE WHEN LAST ALREADY SET
    ENTRY NCOMP3(SVAR,COMM,NC,NEXT)
    IENT=3
    LAST=NEXT+1
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
!...	   A HISTORY COMMAND. RETURN TO 3 IF A PREVIOUS COMMAND SHALL BE EXEC
       CALL NGHIST(SVAR,LAST)
       IF(LAST.EQ.0) GOTO 3
       NCOMP=0
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
       NCOMP=LIKA
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
       NCOMP2=LIKA
    ELSE
       NCOMP3=LIKA
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
  END FUNCTION NCOMP
!========
  SUBROUTINE NGHIST(LINE,LAST)
!...EXECUTES A HISTORY COMMAND
!       LAST IS SET TO 0 IF LINE IS SET TO A COMMAND FROM HISTORY LIST
!    CHARACTER HIST*80,LINE*(*),CH1*1
    CHARACTER LINE*(*),CH1*1
    LOGICAL IED
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
!          IF(KFLAGS(6).EQ.1) THEN
!500          WRITE(KOU,*)
!..CHECK IF WE CAN HANDLE THE LINE
!             NBUF=LEN_TRIM(LINE)
!             NEXT=NBUF
!             JLMARG=1
!             MBUF=MIN(LEN(LINE),79)
!             JAMARG=0
!             CALL EM1BR(JBAUD)
!             JECHO=0
!             JBREAK=0
!             JESC=0
!             JBLANK=1
!             JLF=0
!             JCP=MBUF
!             CALL EM1(JRTN,LINE,NBUF,JCP,MBUF,JLMARG,JAMARG, &
!                  JBAUD,JECHO,JBREAK,JESC,JBLANK,JLF,EM1NFH)
!             IF(NEXT.GT.NBUF)LINE(NBUF+1:)=' '
!             LINE(NBUF+1:)=' '
!..CTRL-N OR LF, CTRL-P EDIT NEXT/PREVIOUS LINE
!             IF(JRTN.EQ.14 .OR. JRTN.EQ.10 ) THEN
!                LOW=LOW+1
!                IF(LOW.GT.20) LOW=1
!                LINE=HIST(LOW)
!                GOTO 500
!             ELSEIF(JRTN.EQ.16) THEN
!                LOW=LOW-1
!                IF(LOW.LE.0) LOW=20
!                LINE=HIST(LOW)
!                GOTO 500
!             ENDIF
!          ELSE
!             CALL FOOLED(LINE)
!          ENDIF
       ENDIF em1
       LAST=0
    ENDIF
900 RETURN
  END SUBROUTINE NGHIST
! =================
  LOGICAL FUNCTION LABBR(NAMN1,NAMN2,EXAKT)
! SUBROUTINE LABBR
! ...TO COMPERE TWO CHARACTER STRINGS WHEN NAMN1 MAY BE AN ABBREVIATION
    CHARACTER NAMN1*(*),NAMN2*(*)
    CHARACTER*1 CH1,CHSEP1,CHSEP2
    PARAMETER (CHSEP1='_',CHSEP2='-')
!    LOGICAL EXAKT,EOLCH
    LOGICAL EXAKT
    L1=LEN_TRIM(NAMN1)
    L2=LEN_TRIM(NAMN2)
    LABBR=.FALSE.
    IF(L1.GT.L2) GOTO 900
    IF(L1.EQ.L2 .AND. NAMN1(1:L1).EQ.NAMN2(1:L1)) THEN
       EXAKT=.TRUE.
       GOTO 800
    ELSE
       EXAKT=.FALSE.
    ENDIF
    J=0
    COMPERE: DO I=1,L1
       CH1=NAMN1(I:I)
       J=J+1
       IF(CH1.EQ.CHSEP1 .OR. CH1.EQ.CHSEP2) THEN
90        IF(NAMN2(J:J).EQ.CHSEP1 .OR. NAMN2(J:J).EQ.CHSEP2) GOTO 100
          J=J+1
          IF(J.GT.L2) GOTO 900
          GOTO 90
       ENDIF
       IF(CH1.NE.NAMN2(J:J)) GOTO 900
100    CONTINUE
    enddo COMPERE
800 LABBR=.TRUE.
900 RETURN
  END FUNCTION LABBR
! ================================
  SUBROUTINE GPARQ(PROMPT,LINE,LAST)
!....ROUTINE TO GIVE ADDITIONAL PROMPTING BEFORE THE NORMAL QUESTION
    CHARACTER PROMPT*(*),LINE*(*)
!    LOGICAL EOLCH
    LAST=LAST+1
    IF(EOLCH(LINE,LAST)) WRITE(KOU,100)PROMPT
100 FORMAT(1X,A)
    LAST=LAST-1
900 RETURN
  END SUBROUTINE GPARQ
! =====================================
  SUBROUTINE GPARLD(PROMT,LINE,NEXT,ITERM,STRING,CDEF,FUNHLP)
! ..INTERACTIVE INPUT OF MULTILINE ANSWER
! ENTRY: PROMT   STRING WITH QUESTION
!        LINE    STRING READ FROM TERMINAL
!        NEXT    FIRST CHARACTER TO SCAN
!        HELP    HELP ROUTINE
! 	 ITERM   TERMINATOR (INTEGER VALUE OF CHARACTER)
! EXIT:  NEXT    UPDATED TO FUNCTION TERMINATOR
!        STRING  INPUT
! 
    EXTERNAL FUNHLP
    CHARACTER PROMT*(*),LINE*(*),STRING*(*),CDEF*(*)
    CHARACTER LINE2*80,TERM*1
!    LOGICAL EOLCH
    IENT=1
    GOTO 10
!...ENTRY POINT WITHOUT DISPLAY OF DEFAULT VALUE
    ENTRY GPARL(PROMT,LINE,NEXT,ITERM,STRING,CDEF,FUNHLP)
    IENT=2
10  IF(IENT.EQ.1) THEN
       CALL GPARCD(PROMT,LINE,NEXT,ITERM,LINE2,CDEF,FUNHLP)
    ELSE
       CALL GPARC(PROMT,LINE,NEXT,ITERM,LINE2,CDEF,FUNHLP)
    ENDIF
!...EMPTY INPUT LINE MEANS ACCEPT DEFAULT IF ANY
    JP=1
    IF(EOLCH(LINE2,JP)) THEN
       IF(CDEF.EQ.'NONE') GOTO 920
       STRING=CDEF
       GOTO 900
    ENDIF
    IP=1
    TERM=CHAR(ITERM)
100 IF(INDEX(LINE2,TERM).LE.0) THEN
       JP=LEN_TRIM(LINE2)
       STRING(IP:)=LINE2
       IP=IP+JP
       IF(IP.GT.LEN(STRING)) GOTO 910
       CALL GPARC('& ',LINE,NEXT,ITERM,LINE2,TERM,FUNHLP)
       GOTO 100
    ENDIF
    JP=LEN_TRIM(LINE2)
    STRING(IP:)=LINE2
    IF(IP+JP.GT.LEN(STRING)) GOTO 910
900 RETURN
!910 CALL ST2ERR(1057,'GPARL','TOO LONG INPUT')
910    buperr=1057
    GOTO 900
!920 CALL ST2ERR(1032,'GPAR','PARAMETER VALUE MISSING')
    920 buperr=1032
    GOTO 900
  END SUBROUTINE GPARLD
!==================================
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
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CH1*1,CDEF*(*),STRING*(*),SSD*30
    CHARACTER PPROMT*132,CH2*1
!    LOGICAL EOLCH,SG2ERR,WDEF,MATP
    LOGICAL WDEF,MATP
    EXTERNAL HELP
    ITYP=1
    WDEF=.TRUE.
    GOTO 9
!---------------
    ENTRY GQARI(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
    ITYP=1
    WDEF=.FALSE.
    GOTO 9
!-----------------------------------------------
!      TO READ A REAL NUMBER
    ENTRY GQARRD(PROMT,SVAR,LAST,VAL,RDEF,HELP)
!...SUBROUTINE GPARR
    ITYP=3
    WDEF=.TRUE.
    GOTO 9
!----------------
    ENTRY GQARR(PROMT,SVAR,LAST,VAL,RDEF,HELP)
    ITYP=3
    WDEF=.FALSE.
    GOTO 9
!-----------------------------------------------
!      TO READ A STRING VALUE
    ENTRY GQARCD(PROMT,SVAR,LAST,JTYP,STRING,CDEF,HELP)
!...SUBROUTINE GPARCD
!      JTYP DEFINES THE TERMINATION OF A STRING
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
    WDEF=.TRUE.
    GOTO 3
!--------------
    ENTRY GQARC(PROMT,SVAR,LAST,JTYP,STRING,CDEF,HELP)
    WDEF=.FALSE.
3   IF(JTYP.LE.0) THEN
!       CALL ST2ERR(1030,'GPARC','NO SUCH TYPE OPTION')
       buperr=1030
       GOTO 900
    ENDIF
    IF(JTYP.EQ.7) THEN
       MATP=.TRUE.
       ITYP=4
    ELSE
       MATP=.FALSE.
       IF(JTYP.LE.6) THEN
          ITYP=JTYP+3
       ELSE
          ITYP=10
          CH2=CHAR(JTYP)
       ENDIF
    ENDIF
! All routines converge here, update command level and save promt 
! for use by help routines
9   CONTINUE
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
             IF(IDEF.NE.NONE) THEN
                X=IDEF
                CALL WRINUM(PPROMT,JJP,10,0,X)
             ENDIF
          ELSEIF(ITYP.EQ.3) THEN
!             IF(RDEF.NE.RNONE) CALL WRINUM(PPROMT,JJP,10,0,RDEF)
! to avoid getting values as 0.0250000004  rather than 0.025 ...
             IF(RDEF.NE.RNONE) CALL WRINUM(PPROMT,JJP,8,0,RDEF)
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
             PPROMT(M+1:M+2)=': '
             JJP=M+2
          ENDIF
          CALL BOUTXT(KOU,PPROMT(1:JJP))
          IJP=JJP
       ELSE
          PPROMT=PROMT
          JJP=LEN_TRIM(PROMT)
          CALL BOUTXT(KOU,PROMT(1:MAX(1,JJP)))
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
!...THE PART HERE IS GLITCHY AS ICE
       L1=0
       L2=0
       GOTO(40,50,60,70,80,70,100),ITYP-3
40     L1=INDEX(SVAR(LAST:),',')
50     L2=INDEX(SVAR(LAST:),' ')
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
!...
60     L1=INDEX(SVAR(LAST:),'.')
70     L2=INDEX(SVAR(LAST:),';')
!...STRING INCLUDING THE ;
       IF(ITYP.EQ.9 .AND. L2.GT.0) L2=L2+1
       GOTO 400
!...WHOLE STRING
80     L1=LEN(SVAR)
       GOTO 400
100    L2=INDEX(SVAR(LAST:),CH2)
400    IF(L1.GT.0 .AND. L2.GT.0) THEN
          L1=LAST+MIN(L1,L2)-1
       ELSEIF(L1.LE.0 .AND. L2.LE.0) THEN
          L1=LEN(SVAR)+1
       ELSE
!...         BUG FOUND HERE: IF L1 IS LEN(SVAR) AND LAST>1 L1 SET >LEN(SVAR)+1
          L1=MIN(LAST+MAX(L1,L2)-1,LEN(SVAR)+1)
       ENDIF
500    IF(L1.GT.LAST) THEN
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
910 IF(ITYP.EQ.1) THEN
       IF(IDEF.EQ.NONE) GOTO 920
       IVAL=IDEF
    ELSEIF(ITYP.EQ.3) THEN
       IF(RDEF.EQ.RNONE) GOTO 920
       VAL=RDEF
    ELSEIF(CDEF.NE.CNONE) THEN
       STRING=CDEF
    ELSE
       GOTO 920
    ENDIF
!...SET POSITION IN STRING TO POSITION OF ,
    LAST=I
    GOTO 900
!...NO ANSWER AND NO DEFAULT VALUE, ERROR RETURN
!920 CALL ST2ERR(1032,'GPAR','PARAMETER VALUE MISSING')
920 buperr=1032
    GOTO 900
  END SUBROUTINE GQARID
! ===================
  subroutine set_echo(ion)
    integer ion
    jecho=ion
1000 continue
    return
  end subroutine set_echo
! ===================
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
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER SVAR*(*),CH1*1,CH2*1,CDEF*(*),STRING*(*)
    character*1, parameter ::  par(4)=['(','{','[','<']
    character*1, parameter :: ipar(4)=[')','}',']','>']
!    LOGICAL EOLCH,SG2ERR
3   IF(JTYP.LE.0) THEN
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
9   LAST=LAST+1
    IF(LAST.LT.1 .OR. LAST.GT.LEN(SVAR)) LAST=LEN(SVAR)+1
!...SKIP BLANKS STARTING FROM THE POSITION AFTER LAST
!      IF LAST OUTSIDE SVAR THEN ASK QUESTION
    I=LAST
10  CONTINUE
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
500 IF(L1.GT.LAST) THEN
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
! ===================
! deleted WRICE, MWRICE, AWRICE, BWRICE as better routine already
! ============
  SUBROUTINE GETINM(SVAR,LAST,IVAL)
! ...IDENTICAL TO GETINT EXCEPT THAT A TERMINATING COMMA ",", IS
!      SKIPPED
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER SVAR*(*)
!    LOGICAL SG2ERR
    CALL GETINT(SVAR,LAST,IVAL)
    IF(BUPERR.NE.0) RETURN
    IF(SVAR(LAST:LAST).EQ.',') LAST=LAST+1
    RETURN
  END SUBROUTINE GETINM
! ==================================================
  SUBROUTINE GETREM(SVAR,LAST,VAL)
! ...IDENTICAL TO GETREL EXCEPT THAT A TERMINATING COMMA "," IS SKIPPED
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER SVAR*(*)
!    LOGICAL SG2ERR
    CALL GETREL(SVAR,LAST,VAL)
    IF(BUPERR.NE.0) RETURN
    IF(SVAR(LAST:LAST).EQ.',') LAST=LAST+1
    RETURN
  END SUBROUTINE GETREM
! ==================
  SUBROUTINE GETINT(SVAR,LAST,IVAL)
!...DECODES AN INTEGER FROM A TEXT
!      IT MAY BE PRECCEDED BY SPACES AND A + OR -
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!    INTEGER GPS
!    LOGICAL EOLCH
    CHARACTER SVAR*(*)
    IF(EOLCH(SVAR,LAST)) THEN
!       CALL ST2ERR(1031,'GETINT','LINE EMPTY')
       buperr=1031
    ELSEIF(SVAR(LAST:MIN(LEN(SVAR),LAST+3)).EQ.'NONE') THEN
       IVAL=NONE
       IERR=0
    ELSE
       IERR=GPS(SVAR,LAST,VALUE)
       IF(IERR.EQ.0) THEN
          IF(VALUE.GT.FLOAT(MAXINT) .OR. VALUE.LT.FLOAT(MININT)) THEN
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
900 RETURN
  END SUBROUTINE GETINT
! ======================
  SUBROUTINE GETOCT(LINE,IP,IVAL)
!...DECODE AN OCTAL NUMBER
!    LOGICAL EOLCH
    CHARACTER LINE*(*)
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
900 RETURN
  END SUBROUTINE GETOCT
! ===========
  SUBROUTINE GETHEX(LINE,IP,IVAL)
!...DECODE A HEXADECIMAL NUMBER
!    LOGICAL EOLCH
    CHARACTER LINE*(*),CH1*1
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
    IF(ISIGN.EQ.1) CALL SETB(1,IVAL)
900 RETURN
  END SUBROUTINE GETHEX
! =============
  SUBROUTINE NOHELP(PROMT,LINE)
    CHARACTER PROMT*(*),LINE*(*)
    RETURN
  END SUBROUTINE NOHELP
! ==============
  SUBROUTINE TOPHLP(PROMPT,LINE)
    CHARACTER PROMPT*(*),LINE*(*)
! do not save the current command ...
!    helprec%level=helprec%level-1
!    write(*,11)helprec%level,(helprec%cpath(i)(1:8),i=1,helprec%level)
11  format('tophlp: ',i3,10(', ',a))
    LINE(2:2)='!'
    RETURN
  END SUBROUTINE TOPHLP
! ============================================
  SUBROUTINE FISEPA(STR,IP0,IP1)
!...FINDS A SEPARATOR AFTER POSITION IP0
!      A separator is:
!      Any character exept A-Z, 0-9 and _
    CHARACTER STR*(*),CH1*1
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
! ====================================
!       INTEGER FUNCTION LENS replaced by LEN_TRIM in F95
!...SUBROUTINE LENS
! ========================
  SUBROUTINE WRIHEX(STR,IVAL)
!...TO WRITE AN INTEGER AS HEXADECIMAL
!    LOGICAL TESTB
    CHARACTER STR*(*)
    J=IVAL
    IP=0
10  IP=IP+1
    K=0
    write(*,*)'calling testb from wrihex'
    IF(TESTB(4*IP-3,J)) K=8
    IF(TESTB(4*IP-2,J)) K=K+4
    IF(TESTB(4*IP-1,J)) K=K+2
    IF(TESTB(4*IP,J)) K=K+1
    IF(K.GT.9) THEN
       STR(IP:IP)=CHAR(K-10+ICHAR('A'))
    ELSE
       STR(IP:IP)=CHAR(K+ICHAR('0'))
    ENDIF
    IF(IP.LT.LEN(STR)) GOTO 10
    RETURN
  END SUBROUTINE WRIHEX
! ==============
  SUBROUTINE SSORT(CMD,NS,INDEX)
!...SORTING
    CHARACTER CMD(*)*(*),STR*40
    DIMENSION INDEX(*)
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
! ==================================================
  SUBROUTINE FDMTP(LINE1,IP,LINE2)
!...FINDS A MATCHING ) AFTER THAT AT IP. IP UPDATED TO POSITION AFTER )
    CHARACTER LINE1*(*),LINE2*(*)
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
! ====================
  LOGICAL FUNCTION YESCHK(CH1)
!...SUBROUTINE YESCHK
    CHARACTER CH1*1
    YESCHK=.FALSE.
    IF(CH1.EQ.'Y' .OR. CH1.EQ.'y') YESCHK=.TRUE.
    RETURN
  END FUNCTION YESCHK
! ==========================
  INTEGER FUNCTION KNDEX(LINE,IP,SS)
! SUBROUTINE KNDEX
!...SEARCHES FOR STRING SS IN LINE FROM IP
    CHARACTER LINE*(*),SS*(*)
    K=INDEX(LINE(IP:),SS)
    IF(K.GT.0) K=IP-1+K
    KNDEX=K
900 RETURN
  END FUNCTION KNDEX
! ==========================
  SUBROUTINE CPSSTR(STRING,LC)
!...THIS SUBROUINE COMPRESSES STRING BY REPLACING MULTIPLE SPACES
!	OR TABS WITH A SINGLE SPACE
    CHARACTER STRING*(*)
!    LOGICAL EOLCH
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
!========
  SUBROUTINE UNTAB(LINE)
!...REMOVES ALL TABS FROM LINE. INSERTS SPACES UP TO NEXT TAB STOP
!       TAB STOPS GIVEN IN ITABS. TABS AFTER POSITION 80 REPLACED
!       WITH A SPACE
    CHARACTER LINE*(*),CHTAB*1,XLINE*128
    DIMENSION ITABS(11)
    DATA ITABS/8,16,24,32,40,48,56,64,72,80,81/
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
900 RETURN
  END SUBROUTINE UNTAB
!========
  SUBROUTINE GPARID(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*)
    CHARACTER SLIN*512
    EXTERNAL HELP
!    CALL GQXENV(SVAR,ENVIR)
    CALL GQXENV(SVAR)
100 CALL GQARID(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
!    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARID
!
  SUBROUTINE GPARI(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*)
    CHARACTER SLIN*80
    EXTERNAL HELP
!    character envir(9)*60
!    COMMON /GENVIR/ENVIR
!    CALL GQXENV(SVAR,ENVIR)
    CALL GQXENV(SVAR)
100 CALL GQARI(PROMT,SVAR,LAST,IVAL,IDEF,HELP)
!    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARI
!
  SUBROUTINE GPARR(PROMT,SVAR,LAST,VAL,RDEF,HELP)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*)
    CHARACTER SLIN*80
    EXTERNAL HELP
!    CHARACTER ENVIR(9)*60
!    COMMON /GENVIR/ENVIR
!    CALL GQXENV(SVAR,ENVIR)
    CALL GQXENV(SVAR)
100 CALL GQARR(PROMT,SVAR,LAST,VAL,RDEF,HELP)
!    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARR
!
  SUBROUTINE GPARRD(PROMT,SVAR,LAST,VAL,RDEF,HELP)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*)
    CHARACTER SLIN*80
    EXTERNAL HELP
!    CHARACTER ENVIR(9)*60
!    COMMON /GENVIR/ENVIR
!    CALL GQXENV(SVAR,ENVIR)
    CALL GQXENV(SVAR)
100 CALL GQARRD(PROMT,SVAR,LAST,VAL,RDEF,HELP)
!    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM1(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARRD
!
  SUBROUTINE GPARC(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CDEF*(*),SVAL*(*)
    CHARACTER SLIN*80
    EXTERNAL HELP
!    CHARACTER ENVIR(9)*60
!    COMMON /GENVIR/ENVIR
!    CALL GQXENV(SVAR,ENVIR)
    CALL GQXENV(SVAR)
100 CALL GQARC(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
    IF(BUPERR.NE.0) GOTO 900
    SLIN=SVAL(1:max(1,LEN_TRIM(sval)))
!    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN)
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARC
!
  SUBROUTINE GPARCD(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER PROMT*(*),SVAR*(*),CDEF*(*),SVAL*(*)
    CHARACTER SLIN*80
    EXTERNAL HELP
!    CHARACTER ENVIR(9)*60
!    COMMON /GENVIR/ENVIR
!    CALL GQXENV(SVAR,ENVIR)
    CALL GQXENV(SVAR)
100 CALL GQARCD(PROMT,SVAR,LAST,JTYP,SVAL,CDEF,HELP)
    IF(BUPERR.NE.0) GOTO 900
    SLIN=SVAL
!    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN,ENVIR)
    CALL GPTCM2(IFLAG,SVAR,LAST,SLIN)
!...the next line was missing ...
    IF (IFLAG.NE.0) GOTO 100
900 RETURN
  END SUBROUTINE GPARCD
!========
  subroutine openlogfile(name,text,lun)
! opens a logfile for commands
    character name*(*),text*(*)
    if(lun.le.0) then
       if(logfil.gt.0) close(logfil)
       goto 1000
    endif
    kkp=index(name,'.')
    if(kkp.le.0) then
       kkp=len_trim(name)
       name(kkp+1:)='.LOG'
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
!========
!  SUBROUTINE GPTCM1(IFLAG,SVAR,LAST,SLIN,ENVIR)
  SUBROUTINE GPTCM1(IFLAG,SVAR,LAST,SLIN)
!...handling of MACRO directives like @& @? and @# etc
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER SVAR*(*),slin*(*)
    CHARACTER PP*30,CH1*1
!    CHARACTER ENVIR(9)*60,LABEL*8,LABLIN*60,SYMBOL*60
    CHARACTER LABEL*8,LABLIN*60,SYMBOL*60
    PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!    LOGICAL SG2ERR,TESTB,EOLCH
!
    IFLAG=0
!...IF NO ERROR RETURN
    IF(.NOT.BUPERR.NE.0) GOTO 900
    IF(LAST.LE.0 .OR. LAST.GE.LEN(SVAR)) GOTO 900
    SLIN=SVAR(LAST:)
!...IF FIRST CHARACTER NOT A @ RETURN
!
!    ENTRY GPTCM2(IFLAG,SVAR,LAST,SLIN,ENVIR)
    ENTRY GPTCM2(IFLAG,SVAR,LAST,SLIN)
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
!...       A PAUSE REQUESTED. CONTINUE AFTER A RETURN, skipp if
!          bit 1 of IOX(8) set
!       write(*,*)'calling testb from gptcm1'
!       IF(.NOT.TESTB(1,IOX(8))) THEN
          WRITE(KOUD,*)'Hit RETURN to continue'
          READ(KIUD,310)CH1
310       FORMAT(A1)
!       ENDIF
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
       call system(slin(2:))
    ENDIF
900 RETURN
!
910 WRITE(LER,911)LABEL
911 FORMAT('No label ',A,' terminating macro')
    SVAR='set interactive '
    GOTO 900
  END SUBROUTINE GPTCM1
!========
!  SUBROUTINE GQXENV(SVAR,ENVIR)
  SUBROUTINE GQXENV(SVAR)
!...EXCHANGES REFERENCES TO ENVIRONMENT MACRO VARIABLES TO ACTUAL VALUES
!       REFERENCES ARE FOR EXAMPLE ##4
!    CHARACTER SVAR*(*),ENVIR(*)*(*),CH1*1,LABLIN*60,LABEL*8
    CHARACTER SVAR*(*),CH1*1,LABLIN*60,LABEL*8
!    COMMON/TCMACRO/IUL,IUN(5),MACEXT
!    CHARACTER*3 MACEXT
    CHARACTER HOLDER*200
!    LOGICAL TESTB
!...BOS 00.10.12: add a possibility to pause if characters are @&
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
  
  SUBROUTINE MACBEG(LINE,LAST,OK)
!....subroutine to execute set-interactive allowing nesting of macros
!
! addera lablar i macro sa man kan ange MACRO fil LABEL
! och vid stop som @? eller @& man kan interaktivt ange GOTO label
! Ocksa en generisk subrutin som gor att man kan fa fram ett variabelvarde
! call macsymval(package,symbol,ival,rval,cval)
!
    CHARACTER LINE*(*),MACFIL*256,FIL*256,CH1*1
!    LOGICAL OK,FIRST,EOLCH,SG2ERR
    LOGICAL OK,FIRST
!    COMMON/TCMACRO/IUL,IUN(5),MACEXT
!    character*3 MACEXT,USEEXT
    character*3 USEEXT
!    common/bincdb/binfil
!    common/tercdb/terfil
!    character*128 binfil(3)
!    character*128 terfil(3)
!    EXTERNAL FILHLP
    SAVE FIRST,LUN
    DATA FIRST/.TRUE./
    IF(FIRST) THEN
       FIRST=.FALSE.
       IUL=0
       lun=50
    ENDIF
    MACFIL=' '
    useext=macext
!    ipos=index(line(max(1,last):),'.')
!    if (ipos.gt.0) then
!       if (LEN_TRIM(line(ipos:)).gt.1) then
!          useext=' '
!       endif
!    endif
!    CALL GPARFD('Macro filename: ',LINE,LAST,1,FIL,MACFIL,USEEXT,FILHLP)
    CALL GPARC('Macro filename: ',LINE,LAST,1,FIL,MACFIL,nohelp)
    CALL FXDFLT(FIL,MACEXT)
!    if (LEN_TRIM(fil).gt.0) call tcgffn(fil)
    IF(BUPERR.NE.0) GOTO 910
!    write(*,*)'open macro: ',iul
!    LUN=50
    OPEN(LUN,FILE=FIL,ACCESS='SEQUENTIAL',STATUS='OLD', &
         FORM='FORMATTED',IOSTAT=IERR,ERR=910)
    IF(IUL.LT.5) THEN
       IUL=IUL+1
       IUN(IUL)=KIU
    ELSE
!       CALL ST2ERR(1083,'MACBEG','TOO DEEPLY NESTED MACRO FILES')
       buperr=1083
       OK=.FALSE.
       GOTO 900
    ENDIF
    CALL GPARC(' ',LINE,LAST,1,CH1,'Y',FILHLP)
    IF(CH1.NE.'Y') THEN
       CALL SETB(1,IOX(8))
    ELSEIF(CH1.EQ.'&') THEN
       CALL CLRB(1,IOX(8))
    ENDIF
    OK=.TRUE.
    KIU=LUN
    LUN=LUN+1
    GOTO 900
!
    ENTRY MACEND(LINE,LAST,OK)
! set interactive gives back control to calling macro if any
    IF(KIU.NE.KIUD) THEN
       IF(KIU.NE.0) CLOSE(KIU)
!       write(*,*)'end of macro: ',kiu,kiud,iul
       IF(IUL.GT.0) THEN
!          write(*,*)'calling macro: ',iun(iul)
          KIU=IUN(IUL)
          IUL=IUL-1
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
900 RETURN
910 OK=.FALSE.
    buperr=1000+ierr
    GOTO 900
  END SUBROUTINE MACBEG
  
!============================================
! to be implemented
!============================================

  logical function testb(ibit,iword)
! obsolete, use library fundtion btest with inverse argument order
! and inverse bit numbering ...
    write(*,*)' **** warning using obsolete testb *** '
    stop 'use btest instead of testb'
    if(ibit.le.0 .or. ibit.gt.32) then
       buperr=1060
    else
       testb=btest(iword,ibit)
    endif
  end function testb

  subroutine setb(ibit,iword)
    stop 'setb'
!    call ibset(iword,ibit)
  end subroutine setb

  subroutine clrb(ibit,iword)
    stop 'clrb'
!    call ibclr(iword,ibit)
  end subroutine clrb

  subroutine boutxt(lut,line)
    character line*(*)
!    write(*,*)'boutxt; ',lut,line
    write(lut,10,advance='no')line
10  format(a)
    return
  end subroutine boutxt

  subroutine bintxt(lin,line)
    character line*(*)
    integer iostatus
    read(lin,10,iostat=iostatus)line
10  format(a)
    if(iostatus.lt.0) then
! reading beyond EOL
       line='fin '
    endif
    return
  end subroutine bintxt

  function ionoff(subr)
    character subr*(*)
    ionoff=0
  end function ionoff

  subroutine comnd(line)
    character line*(*)
    stop 'comnd'
  end subroutine comnd

  subroutine tcgffn(line)
    character line*(*)
    stop 'tcgffn'
  end subroutine tcgffn

  subroutine fxdflt(file,ext)
! add default file extention, no good as it thinks .. is an externtion
    character file*(*),ext*(*)
    kx=index(file,'.')
    if(kx.le.0) then
       kx=len_trim(file)
       file(kx+1:)='.'//ext
    endif
    return
  end subroutine fxdflt

  subroutine GPARFD(quest,LINE,LAST,kk,FIL,MACFIL,USEEXT,hlprut)
    character quest*(*),line*(*),fil*(*),macfil*(*),useext*(*)
    external hlprut
    stop 'gparfd'
  end subroutine GPARFD

  subroutine iniio
!    write(*,*)'iniio: ',koud,kiud,lerd
    kou=koud
    kiu=kiud
    ler=lerd
    iox=0
    return
  end subroutine iniio

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
! Rather final version of PUTFUN below     
!
! MODULE PUTFUNLIB
!
  SUBROUTINE PUTFUN(STRING,L,MAXS,SYMBOL,LOKV,LROT,ALLOWCH,NV)
!...READS AN EXPRESSION FROM STRING POSITION L AND CREATES AN BINARY TREE
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER STRING*(*),CH*1,SYMBOL(*)*(*),MESSG*40
    PARAMETER (ZERO=0.0D0)
    DIMENSION LOKV(*)
    integer allowch
!    type(putfun_symlink) :: symlist
    LOGICAL NOTPP
    TYPE(putfun_node), pointer :: LROT,NYNOD,NONOD,DUMMY,dummy2
! ENTRY:
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
10  continue
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

  SUBROUTINE NYBIN(kod,binnod,NOTPP)
!...INSERTS A NEW OPNODE IN THE TREE
    TYPE(putfun_node), pointer :: binnod,temp
    LOGICAL NOTPP
    integer, parameter :: zero=0.0d0
! ENTRY:
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

  SUBROUTINE NYUNI(KOD,negmark,uninod,IPN,NOTPP)
!   Creates a node with a unary operator
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    TYPE(putfun_node), pointer :: UNINOD
    LOGICAL NOTPP
    parameter (one=1.0D0)
!
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
900 RETURN
  END SUBROUTINE NYUNI

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE NYLP(uninod,IPN,NOTPP)
!...OPENING PARENTHESIS, push links on LEVEL. Also after unary operator
    TYPE(putfun_node), pointer :: uninod
    type(putfun_stack), pointer :: temp
    LOGICAL NOTPP
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

  subroutine NYRP(IPN,NOTPP)
!...CLOSING PARENTHESIS
    implicit double precision (a-h,o-z)
    TYPE(putfun_node), pointer :: uninod,subtree
    LOGICAL NOTPP
100 continue
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

  SUBROUTINE NYVAR(TEXT,L,IOPUNI,negmark,MAXS,SYMBOL,LOKV,allowch,dummy2)
!  SUBROUTINE NYVAR(TEXT,L,IOPUNI,negmark,MAXS,SYMBOL,LOKV,symlist)
! insert a symbol
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (NOPER=14)
    CHARACTER TEXT*(*),SYMBOL(*)*(*),CH*1,NAME*16
    PARAMETER (ZERO=0.0D0)
    LOGICAL DEL2
    DIMENSION LOKV(*)
    integer allowch
    type(putfun_node), pointer :: dummy2
    parameter (one=1.0d0)
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

  SUBROUTINE NYDAT(KOD,VAL,nynod,negmark)
! store a constant or symbol.  The address to the node is returned in lok
! which is used if the symbol is used several times.
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    TYPE(putfun_node), pointer :: nynod
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

  FUNCTION EVALF(LROT,VAR)
!      Calculates the value of an expression
!
! VAR is array with values of symbols that can be referenced
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION VAR(*),STACK(20)
    character ch1*1
    type(putfun_node), pointer :: lrot,llink,current,mlink
    TYPE PUTFUN_SAVE
       integer right
       type(putfun_node), pointer :: savecurrent
       type(putfun_save), pointer :: previous
    end TYPE PUTFUN_SAVE
    type(putfun_save), pointer :: topsave,temp
    PARAMETER (ZERO=0.0D0)
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
800 EVALF=STACK(1)
900 RETURN
  END FUNCTION EVALF

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE EUNARY(KOD,X)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (ONE=1.0D0)
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

  SUBROUTINE EBINRY(KOD,X,Y)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    ZERO=0.0D0
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

!...added by Zikui and also an updated ERF
  FUNCTION AIVAN(PECN)
!      SUBROUTINE FUNCTION AIVAN
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      CALCULATES THE DIMENSIONLESS SUPERCOOLING OF DIFFUSION BY
!      IVANTSOV'S SOLUTION
!      APPROXIMATIVE FORMULA FOR ERROR FUNCTION GIVEN BY:
!      ABRAMOWITZ AND STEGUN: HANDBOOK OF MATHEMATICAL FUNCTIONS,
!      NATIONAL BUREAU OF STANDARDS, 9TH EDITION, 1970
    PARAMETER (ONE=1.0D0,TWO=2.0D0,PI=3.141592654D0)
    IF(PECN.LE.8.5D0) THEN
       AIVAN=DSQRT(PI*PECN)*DEXP(PECN)*(ONE-ERF(DSQRT(PECN)))
    ELSE
       A=ONE
       C=ONE
       Q=ONE
       DO I=1,9
          A=A*(TWO*DFLOAT(I)-ONE)/TWO/PECN
          C=-C
          Q=Q+A*C
       enddo
       AIVAN=Q
    ENDIF
    RETURN
  END FUNCTION AIVAN

  FUNCTION PF_BSUM(FA)
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
!      SUBROUTINE FUNCTION BSUM
    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    PARAMETER(ZERO=0.0D+00,PI=3.14159D0)
    PARAMETER(PI3=PI*PI*PI)
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

  FUNCTION PF_HS(X)
!      SUBROUTINE FUNCTION HS
!      Calculates Heaviside function
    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    PARAMETER(ZERO=0.0D+00,ONE=1.0D+00)
    HS=ZERO
    IF (X.GE.ZERO) PF_HS=ONE
    RETURN
  END FUNCTION PF_HS

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  FUNCTION PF_ERF(X0)
!      SUBROUTINE FUNCTION ERF
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      CALCULATES ERRR-FUNCTION OF X, USING AN
!      APPROXIMATIVE FORMULA GIVEN BY:
!      ABRAMOWITZ AND STEGUN: HANDBOOK OF MATHEMATICAL FUNCTIONS,
!      NATIONAL BUREAU OF STANDARDS, 9TH EDITION, 1970
    PARAMETER (ONE=1.0D0,TWO=2.0D0)
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

  SUBROUTINE WRTFUN(STRING,IPOS,LROT,SYMBOL)
!      Writes an expression
!
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    type(putfun_node), pointer :: lrot,current,lnode,rnode,tnode
    TYPE PUTFUN_SAVE
       integer right
       type(putfun_node), pointer :: savecurrent
       type(putfun_save), pointer :: previous
    end TYPE PUTFUN_SAVE
    type(putfun_save), pointer :: topsave,temp,bug
!    type(putfun_node), dimension(20) :: link
!    integer doneboth(20)
   CHARACTER STRING*(*),SYMBOL(*)*(*)
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
76  format(a,2i4,1pe15.6)
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

  SUBROUTINE WRTLPQ(STRING,IPOS,LINK,KOD,LOD,negmark)
! write a left ( or unary operator followed by (
! the unary operator is in LOD
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER STRING*(*)
!    IF(LOD) 10,90,20
    IF(LOD.eq.0) goto 90
    if(LOD.gt.0) goto 20
!..unary operator
10  IF(LOD.EQ.-1) THEN
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

  SUBROUTINE WRTRPQ(STRING,IPOS,LINK,KOD,LOD)
!  write a right )  but if LOD<-1 do not write (
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER STRING*(*)
!    IF(LOD+1) 90,10,20
    IF(LOD+1.lt.0) goto 90
    if(LOD+1.gt.0) goto 20
!..negation need ) if KOD>0
10  IF(KOD.GT.0) CALL CONS(STRING,IPOS,')')
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

  SUBROUTINE WRTBIQ(STRING,IPOS,KOD)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER STRING*(*)
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

  SUBROUTINE WRTDAQ(STRING,IPOS,KOD,VAL,SYMBOL,negmark)
!     write a number, if KOD<0 a whole number
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER NAME*8,SYMBOL(*)*(*)
    CHARACTER STRING*(*)
    PARAMETER (ZERO=0.0D0)
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

  SUBROUTINE DELFUN(LROT,IWS)
!   delete a putfun expression :: not converted to structures
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION IWS(*),LINK(20)
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
!..right link visited, remove binary operator
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

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE CONS(STR1,IPOS,STR2)
!      CONS. TWO STRINGS, RESULT IN PARAMETER STR1
!      IPOS = POSITION IN STR1 WHERE STR2 SHOULD BE PUT
!      IPOS IS UPPDATED TO THE FIRST FREE POSITION AT THE END
!      OF STR1. TRAILING SPACES ARE STRIPPED OFF.
!      IF STR2 CONTAINES ONLY SPACES ONE SPACE IS WRITTEN
!      IN TO STR1.
    CHARACTER STR1*(*),STR2*(*),SPC*(1)
    PARAMETER (SPC=' ')
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

   SUBROUTINE GZRFUN(PROMT,LINE,NEXT,STRING,FUNHLP)
!x     GZRFUN(PROMT,LINE,NEXT,STRING,FUNHLP)
!x..INTERACTIVE INPUT OF A FUNCTION FROM THE TERMINAL
!xENTRY: PROMT   STRING WITH QUESTION
!x       LINE    STRING READ FROM TERMINAL
!x       NEXT    FIRST CHARACTER TO SCAN
!x       FUNHLP  HELP ROUTINE
!xEXIT:  NEXT    UPDATED TO FUNCTION TERMINATOR (;)
!x       STRING  FUNCTION
!x
    EXTERNAL FUNHLP
    CHARACTER PROMT*(*),LINE*(*),STRING*(*),LINE2*300
    IP=1
!    CALL GPARC(PROMT,LINE,NEXT,6,LINE2,';',FUNHLP)
100 IF(INDEX(LINE2,';').LE.0) THEN
       JP=LEN_trim(LINE2)
       STRING(IP:)=LINE2
       IP=IP+JP
       IF(IP.GT.LEN(STRING)) GOTO 910
!       CALL GPARC('& ',LINE,NEXT,6,LINE2,';',FUNHLP)
       GOTO 100
    ENDIF
    JP=LEN_trim(LINE2)
    STRING(IP:)=LINE2
    IF(IP+JP.GT.LEN(STRING)) GOTO 910
900 RETURN
910 continue
!910 CALL ST2ERR(1057,'GZRFUN','TOO LONG FUNCTION')
    GOTO 900
  END SUBROUTINE GZRFUN

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE EXPHLP(PROMPT,SVAR)
    CHARACTER PROMPT*(*),SVAR*(*)
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
900 RETURN
  END SUBROUTINE EXPHLP

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!

  SUBROUTINE PUTPRP(NAMN,MAXS,SYMBOL,PROMPT,ILEN)
!...CREATES A PROMPT
    CHARACTER NAMN*(*),PROMPT*(*),SYMBOL(*)*(*)
!...write a prompt with name of all variables
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

! end module putfunlib
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  SUBROUTINE HPCALC
!...EMULATES A HP CALCULATOR ON SCREEN
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION STK(4),REG(0:9)
    PARAMETER (NOP=40,MAXPRG=100)
    CHARACTER LINE*80,INPUT*80,OPER(NOP)*10,CH1*1,CH2*2
    CHARACTER PROG(MAXPRG+1)*20
!    LOGICAL SG2ERR,PROGT,OK,RUN,TRACE
    LOGICAL PROGT,OK,RUN,TRACE
    SAVE STK,REG,PROG
!    EXTERNAL HPHLP
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
         OPER(39)/'          '/,OPER(40)/'          '/
!
    WRITE(*,*)'REVERSE POLISH CALCULATOR'
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
          PROG(LPROG)=INPUT
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
       CALL GPARC('HPP>',LINE,LAST,1,INPUT,' ',HPHLP)
       CALL CAPSON(INPUT)
       IF(INPUT(1:1).EQ.' ') GOTO 100
       OK=.FALSE.
    ELSE
       WRITE(*,101)STK(1)
101    FORMAT(1PE15.7)
       CALL GPARC('HPC>',LINE,LAST,1,INPUT,' ',HPHLP)
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
!    IF(.NOT.SG2ERR(IERR)) THEN
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
    CALL NGHELP(INPUT,NEXT,OPER,NAXOP)
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
       CALL GPARI('REG#',LINE,LAST,I,-1,NOHELP)
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
       CALL GPARI('REG#',LINE,LAST,I,-1,NOHELP)
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
       CALL GPARI('STEP ',LINE,LAST,IP,0,NOHELP)
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
       CALL GPARC('CONDITION ( > = OR < )',LINE,LAST,1, &
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
421 KPROG=KPROG+1
    IF(KPROG.GT.LPROG) THEN
       WRITE(*,*)'PROGRAM ENDS AT ',LPROG
    ELSE
       WRITE(*,401)KPROG,PROG(KPROG)
    ENDIF
    OK=.FALSE.
    GOTO 100
!...EDIT
430 IF(JPROG.EQ.0) JPROG=LPROG-2
    CALL GPARI('STEP ',LINE,LAST,I,JPROG+1,NOHELP)
!    IF(I.GT.0 .AND. I.LT.LPROG) THEN
!       CALL FOOLED(PROG(I))
!    ELSE
       WRITE(*,*)'NO SUCH STEP'
!    ENDIF
    OK=.FALSE.
    GOTO 100
!...BACK
440 WRITE(*,*)'NOT IMPLEMENTED'
    GOTO 100
441 WRITE(*,401)KPROG,PROG(KPROG)
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
490 CONTINUE
500 CONTINUE
    GOTO 100
  END SUBROUTINE HPCALC

  SUBROUTINE HPHLP
    WRITE(*,10)
10  FORMAT(' This is a revese polish calculator'/&
         ' Input are numbers, + - * / and ^ and OPCODEs.',&
         ' Use HELP to list OPCODEs.'/' Several numbers an operations',&
         ' can be given on one line.'/' The content of the X register',&
         ' is displayed after each operation'//&
         ' Example input: 30000 8 1273 * / chs 1.5 3 ^ + exp 2 *'/&
         ' Computes 2*EXP(1.5**3-30000/(8*1273))'//)
    RETURN
  END SUBROUTINE HPHLP

!---------------------------------------------------------
! A new set of on-line help routines

  subroutine init_help(file)
! This routine is called from gtpini to inititate the on-line help system
    character*(*) file
    character*80 line
! test that file exists
    open(21,file=file,access='sequential',status='old',&
         err=900,iostat=jerr)
    read(21,10)line
10  format(a)
    if(line(2:14).eq.'documentclass') then
       helprec%type='latex   '
    else
       helprec%type='unknown '
    endif
    close(21)
    helprec%okinit=1
    helprec%filename=file
    goto 1000
900 continue
    write(*,910)file(1:len_trim(file))
910 format('Warning, cannot open ',a,' no on-line help')
    helprec%okinit=0
    helprec%filename=' '
1000 continue
    return
  end subroutine init_help

  subroutine helplevel1(line)
! This routine is called from the monitor for the top level command
! It initiates the path to find the correct help text
! In all gparx routines the help level in increased and the question saved
    character*(*) line
    helprec%level=1
    helprec%cpath(1)=line
1000 continue
    return
  end subroutine helplevel1

  subroutine q2help(prompt,line)
! This routine is called from submenu
! when the user types a ?
    character*(*) prompt,line
    character helpquest*16
    integer savedlevel
    savedlevel=helprec%level-1
! If the ? is followed by a text push that text on the helprec%cpath
!    write(*,*)'q2help: ',savedlevel,line(1:20)
    ip=2
    if(.not.eolch(line,ip)) then
!       write(*,*)helprec%level,helprec%cpath(helprec%level)
       helpquest=line(ip:)
       call capson(helpquest)
! use the saved helprec%level 
       helprec%level=savedlevel
       helprec%cpath(helprec%level)=helpquest
!       write(*,11)helprec%level,(helprec%cpath(i)(1:8),i=1,helprec%level)
11     format('q2help: ',i3,10(', ',a))
    endif
! write help text from help file and then return with ?! to get submenu
    call q1help(prompt,line)
    line='?!'
    return
  end subroutine q2help

  subroutine q1help(prompt,line)
! This routine is called from all gparx routines 
! when the user types a ?
! prompt is never used ...
    implicit none
    character*(*) prompt,line
    character hline*80,mtxt*80
    character subsec(4)*10,saved(4)*16
    integer nsaved(4)
    integer izz,jj,kk,kkk,level,nl,l2,np1,np2,nsub
!
!    write(*,*)'called on-line help',helprec%okinit
!    if(helprec%okinit.ne.0) then
!       write(*,*)'help file: ',helprec%filename(1:len_trim(helprec%filename))
!    endif
!    write(*,*)'we are on command path level: ',helprec%level
!    do i=1,helprec%level
!       write(*,10)i,helprec%cpath(i)
!    enddo
!10  format(i3,': ',a)
!
    subsec(1)='%\section{'
    subsec(2)='%\subsecti'
    subsec(3)='%\subsubse'
    subsec(4)='%\subsubsu'
    if(helprec%okinit.eq.0) then
       write(kou,*)'Sorry no help file'
       goto 1000
    endif
! for debugging list current search path:
!    do nl=1,helprec%level
!       write(*,17)'Search levels: ',nl,helprec%cpath(nl)(1:24)
!17     format(a,i3,2x,a)
!    enddo
!
    open(31,file=helprec%filename,status='old',access='sequential')
    nl=0
    level=3
    np1=0
    nsub=1
    if(helprec%type.eq.'latex   ') then
! plain LaTeX file, search for "%\section{" with command
100    continue
       if(level.gt.helprec%level) then
          level=level-1
       endif
       l2=len_trim(helprec%cpath(level))
!       write(*,107)level,l2,helprec%cpath(level)(1:l2)
107    format('Command/question: ',2i3,': ',a)
110    continue
       read(31,120,end=700)hline
120    format(a)
121    format('- ',a)
       nl=nl+1
       kk=index(hline,subsec(nsub))
       if(nsub.gt.1) then
          izz=index(hline,subsec(nsub-1))
       else
          izz=0
       endif
       if(kk.gt.0) then
!          write(*,*)'found section: ',hline(1:30),nl
          call capson(hline)
!          kk=kk+8
          kk=index(hline,'SECTION{')+8
          if(kk.le.8) then
             write(*,*)'Help file error, not correct LaTeX',kk
             np1=0
             goto 700
          endif
! remove everything after and including } in the help file
          jj=index(hline,'}')-1
          hline(jj+1:)=' '
          mtxt=hline(kk:jj)
!127       kkk=index(mtxt,'-')
!          if(kkk.gt.0) then
! NOTE that in LaTeX files the underscore _ is replaced by \_
! this means one will not find commands with underscore ....
!             mtxt(kkk:kkk)='_'
!             goto 127
!          endif
!          write(*,*)'found mtxt: ',mtxt
          do jj=1,4
             if(level.gt.jj) then
! remove previous command if any from the search text
                kkk=nsaved(jj)
                if(mtxt(1:kkk+1) .eq. saved(jj)(1:kkk+1)) then
!                   write(*,*)'Removing ',mtxt(1:kkk)
                   mtxt=mtxt(kkk+2:)
                endif
             endif
          enddo
!          write(*,130)level,mtxt(1:12),helprec%cpath(level)(1:l2),nl
!          write(*,130)level,hline(kk:jj),helprec%cpath(level)(1:l2),nl
130       format('comparing: ',i3,': ',a,' with ',a,' line ',i4)
!          if(hline(kk:kk+l2-1).eq.helprec%cpath(level)(1:l2)) then
          if(mtxt(1:12).eq.helprec%cpath(level)(1:l2)) then
! found one level, save line number for first line
             np1=nl
!             write(*,*)'We found one level at line: ',np1
             saved(level)=helprec%cpath(level)
             nsaved(level)=len_trim(saved(level))
             np2=0
             level=level+1
             if(level.gt.helprec%level) then
! we have no more input from user, list the text to next %\section or %\sub...
                level=level-1
                goto 100
             endif
! look for start line of next level
!             write(*,133)level,helprec%cpath(level)
133          format('Next cpath: ',i3,': ',a)
             if(index(helprec%cpath(level),'COMMAND: ').gt.0 .or. &
                  index(helprec%cpath(level),' what? ').gt.0) then
! if next "command/question" starts with COMMAND. or finishes with "what?" 
! skip that because it is a submenu question and if this is
! the last command line then just return as user interface will help
!                write(*,*)'level and helprec%level',level,helprec%level,nsub
                if(level.ge.helprec%level) goto 900
                level=level+1
                goto 100
             endif
             if(nsub.lt.3) nsub=nsub+1
             goto 100
          elseif(level.eq.helprec%level .and. np1.gt.0 .and. np2.eq.0) then
! save end of text
             np2=nl-1
             goto 700
          else
! continue searching
             goto 100
          endif
       elseif(izz.gt.0) then
! we have found the start of a higher order section, finish searching
          np2=nl-1
          goto 700
       endif
       goto 110
    else
! One should also implement HTML help files
       write(*,*)'Filetype not understood'
       goto 1000
    endif
! we should write lines from np1 to np2 from help file
700 continue
    if(np1.gt.0) then
       if(np2.lt.np1) then
          write(*,*)'Help text range error: ',np1,np2
          goto 900
       endif
! HERE ONE SHOULD OPEN A NEW WINDOW TO DISPLAY THE HELP TEXT BELOW, FOR EXAMPLE
!       call displayhelp(helprec%filename,np1,np2)
! write a blank line
       write(kou,*)
! OR IF A HELP WINDOW IS ALREADY OPEN ONE SHOULD FOCUS THE TEXT IN THAT WINDOW
! ON LINES NP1 to NP2.
! IN THE HELP WINDOW THE USER SHOULD BE ABLE TO SCROLL THE WHOLE HELP TEXT
! AND MAYBE MAKE NEW SEACHES.  THE USER CAN CLOSE THE WINDOW WHEN HE WANTS
       write(*,798)np1,np2
798    format(' >>> We should open a help window to display text: ',2i5)
       rewind(31)
       nl=0
800    continue
       read(31,120)hline
       nl=nl+1
       if(nl.ge.np2) then
          goto 900
       elseif(nl.ge.np1) then
! do not write lines starting with \ (backslash, ascii 92),
! they are LaTeX commands.  Ignore all except \item
          if(ichar(hline(1:1)).ne.92) then
             write(*,120)hline(1:len_trim(hline))
          elseif(hline(2:5).eq.'item') then
             write(*,121)hline(6:len_trim(hline))
          endif
       endif
       goto 800
    else
       write(*,*)'No help found'
    endif
900 continue
    close(31)
!
1000 continue
    return
  end subroutine q1help

  SUBROUTINE Q3HELP(LINE,LAST,COMM,NC)
!  SUBROUTINE Q3HELP(LINE,LAST,COMM,NC)
! used in submeny when user gives "? 'command' " taken as "help 'command'"
!...EXECUTES A HELP COMMAND
    CHARACTER LINE*(*),COMM(NC)*(*),CMD*40
    PARAMETER (MC=100)
    DIMENSION INDX(MC)
!    CALL GPARC('COMMAND: ',LINE,LAST,1,CMD,'*',FILHLP)
!    CALL GPARC('COMMAND: ',LINE,LAST,1,CMD,'*',q1help)
! To avoid storing "COMMAND" in the helprec%cpath
!    if(helprec%level.gt.2) helprec%level=helprec%level-1 !.. HELP HELP not OK
    if(helprec%level.gt.3) helprec%level=helprec%level-1
!    write(*,*)'q3help: asking for command for help:',line(1:2),last
    CALL GPARC('COMMAND: ',LINE,LAST,1,CMD,'*',tophlp)
!    write(*,*)'q3help: command:',cmd
    IF(CMD(1:1).EQ.'*') THEN
!...LIST ALL COMMANDS IN UNIX ALPHABETICAL ORDER
       NKPL=80/(LEN(COMM(1))+1)
       IF(NKPL*(LEN(COMM(1))+1).GE.80) NKPL=NKPL-1
100    IF(NC.LT.MC) THEN
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
          helprec%cpath(helprec%level)=CMD
!          write(*,11)helprec%level,(helprec%cpath(i)(1:8),i=1,helprec%level)
11        format('q3help: ',i3,10(', ',a))
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

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!
! WPACK
!
  SUBROUTINE WINIT(NWT,NWR,IWS)
!...INITIATES A WORKSPACE
! ENTRY: NWT IS THE DIMENSION OF THE WORKSPACE
!        NWR IS THE NUMBER OF WORDS TO BE EXCLUDED IN THE BEGINNING
!        IWS IS THE WORKSPACE
! EXIT:  THE FREE LIST IS INITIATED IN IWS
! ERRORS: NWR LESSER THAN ZERO
!         NWT LESSER THAN NWR+100
    DIMENSION IWS(*)
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
  
  SUBROUTINE WOLD(FIL,NW,IWS)
!...READS A FILE INTO A WORKSPACE. THE FILE MUST HAVE BEEN WRITTEN BY WSAVE
! ENTRY: FIL A CHARACTER WITH A LEGAL FILE NAME
!        NW THE DIMENSION OF IWS
!        IWS THE WORKSPACE
! CALLS: WRKCHK TO CHECK THE FREE LIST
! EXIT:  THE CONTENT OF THE FILE IS IN IWS. THE DIMENSION OF IWS IS SET TO
!            NW AND THE LAST FREE AREA IS CORRECTED
    CHARACTER FIL*(*)
    DIMENSION IWS(*)
    OPEN(UNIT=LUN,FILE=FIL,ACCESS='SEQUENTIAL',STATUS='OLD',&
         IOSTAT=IERR,ERR=910,FORM='UNFORMATTED')
    READ(LUN,END=100,ERR=100)J,(IWS(K),K=1,J)
!...CHECK THE WORKSPACE
    CALL WRKCHK(LAST,NW,IWS)
100 CLOSE(LUN)
900 RETURN
910 continue
    GOTO 100
  END SUBROUTINE WOLD

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE WSAVE(FIL,NW,IWS)
!...WRITES A WORKSPACE ON A FILE
! ENTRY: FIL IS A CHARACTER WITH A LEGAL FILE NAME
!        NW IS THE DIMENSION OF THE WORKSPACE
!        IWS IS THE WORKSPACE
! CALLS: WRKCHK TO CHECK THE WORKSPACE
! ERROR: IF THE WORKSPACE IS INCORRECT IT CANNOT BE SAVED
    DIMENSION IWS(*)
!    LOGICAL SG2ERR
    CHARACTER FIL*(*)
    I=IWS(2)
    CALL WRKCHK(LAST,I,IWS)
!    IF(SG2ERR(IERR)) GOTO 900
    if(buperr.ne.0) goto 900
    OPEN(UNIT=LUN,FILE=FIL,ACCESS='SEQUENTIAL',STATUS='UNKNOWN',&
         IOSTAT=IERR,ERR=910,FORM='UNFORMATTED')
    WRITE(LUN,ERR=910)LAST,(IWS(I),I=1,LAST)
800 CLOSE(LUN)
900 RETURN
910 continue
    GOTO 800
  END SUBROUTINE WSAVE

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE WPATCH(NW,IWS)
!...ROUTINE TO PATCH A WORKSPACE
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION IWS(*)
    CHARACTER LINE*80,CHX*(NBPW),CHHEX*(2*NBPW)
!    EXTERNAL WPHLP
!    LOGICAL SG1ERR,SG2ERR,EOLCH
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
!      READ(KIU,*)LINE
    CALL BINTXT(KIU,LINE)
    IP=1
    CALL GETINT(LINE,IP,IADR)
!    IF(SG2ERR(IERR)) THEN
    if(buperr.ne.0) then
!       CALL RESERR
       buperr=0
       IF(LINE(1:1).EQ.'?') CALL WPHLP(IP,LINE)
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
111 CALL GPARR('NEW VALUE: ',LINE,IP,Z,RNONE,WPHLP)
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

  SUBROUTINE WPHLP(ITYP,LINE)
!...HELP ROUTINE FOR WPATCH
    CHARACTER LINE*(*)
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

  SUBROUTINE WRKCHK(LAST,NW,IWS)
!...CHECKS THE FREE LIST IN A WORKSPACE
! ENTRY: NW IS THE DIMENSION
!        IWS IS THE WORKSPACE
! EXIT:  LAST IS PUT TO THE LAST WORD USED IN THE WORKSPACE
! ERRORS: ANY ERROR IN THE FREE LIST (POINTER OUTSIDE WORKSPACE ETC)
    DIMENSION IWS(*)
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

  SUBROUTINE WLIST(IWS)
!...LISTS THE FREE AREAS
    DIMENSION IWS(*)
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

  SUBROUTINE WTREST(NYB,NW,IWS)
!...RESERVES THE LAST PART OF THE WORKSPACE
! ENTRY: IWS IS A WORKSPACE
! EXIT:  NYB IS A POINTER TO THE RESERVED PART
!        NW IS THE NUMBER OF RESERVED WORDS
    DIMENSION IWS(*)
    LOK=1
100 LAST=LOK
    LOK=IWS(LAST)
    IF(LOK.GT.0) GOTO 100
    NW=IWS(LAST+1)
    CALL WTAKE(NYB,NW,IWS)
900 RETURN
  END SUBROUTINE WTREST

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE WTAKE(NYB,NW,IWS)
!......RESERVS NW WORDS IN THE WORKSPACE
! ENTRY: NW IS THE NUMBER OF WORDS TO BE RESERVED
!        IWS IS THE WORKSPACE
! EXIT:  NYB POINTS TO THE FIRST WORD THAT IS RESERVED
! ERROR: TOO SMALL OR TOO LARGE NUMBER OF WORDS TO BE RESERVED
    DIMENSION IWS(*)
!...THE FREE LIST START IN THE FIRST WORD
!      IN EACH FREE AREA THE FIRST WORD POINTS TO THE NEXT FREE AREA
!      AND THE SECOND GIVES THE NUMBER OF WORDS IN THIS AREA
!      THE FREE LIST ENDS WITH THE POINTER EQUAL TO ZERO
    IF(NW.LT.2) GOTO 910
    LOKB=1
4   LOKA=IWS(LOKB)
5   IF(LOKA.LE.0) GOTO 920
    IF(LOKA.GE.IWS(2)) GOTO 930
    IF(IWS(LOKA+1)-NW) 10,20,30
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

  SUBROUTINE WRELS(IDP,NW,IWS)
!......Returns NW words beginning from IDP to the free workspace list
!      The free workspace list is in increasing order
!      IWS(1) points to the first free space
!      IWS(2) gives the total number of words in the workspace
    DIMENSION IWS(*)
!......Check that the released space is at lest 2 words and that it is
!      inside the workspace (That is between 3 and IWS(2))
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
    IF(LOKA+IWS(LOKA+1)-LOKC) 120,110,930
!..The released space follows directly on LOKA => Merge LOKA and LOKC
110 LOKC=LOKA
    IWS(LOKC+1)=IWS(LOKC+1)+NW
    GOTO 130
!..Set the pointer from LOKC to LOKB and from LOKA to LOKC
120 IWS(LOKC)=LOKB
    IWS(LOKA)=LOKC
    IWS(LOKC+1)=NW
!..Check if LOKC now can be merged with LOKB!
130 IF(LOKC+IWS(LOKC+1)-LOKB) 900,140,940
!..Merge LOKC and LOKB
140 IWS(LOKC)=IWS(LOKB)
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

  INTEGER FUNCTION NWCH(NB)
! subroutine nwch
! number of words to store a character with nb bytes
! nbpw is the number of bytes in a word.  If not even multiple add 1 word
    i=nb/nbpw
    if(mod(nb,nbpw).gt.0) then
       i=i+1
    endif
    nwch=i
    return
  end FUNCTION NWCH

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE STORC(N,IWS,text)
! Stores a text character in an integer workspace at position N
! The length of the character to store is len(text)
    integer n,iws(*),llen
    character text*(*)
! maximal size of character, note used also to store functions and bibliography
    integer, parameter :: maxchar=2048,maxequiv=512
! NOTE BELOW DIMENSIONING BLEOW, maxchar=nbpw*maxequiv
!    character (len=:), allocatable :: localtxt
!    integer, allocatable, dimension(:) :: localint
    character*(maxchar) localtxt,localtxt2
! assumed 32 bit integres, 8 bits/character, 4 characters/word =nbpw
    integer localint(maxequiv),localint2(maxequiv)
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

  SUBROUTINE LOADC(N,IWS,text)
! copies a text from an integer workspace at position N into a character
! The number of characters to copy is len(text)
    integer n,iws(*),llen
    character text*(*)
!    character (len=:), allocatable :: localtxt
!    integer, allocatable, dimension(:) :: localint
! maximal size of character, note used also to store functions and bibliography
    integer, parameter :: maxchar=2048,maxequiv=512
! NOTE BELOW DIMENSIONING BELOW, maxchar=nbpw*maxequiv
    character*(maxchar) localtxt
! assumed 32 bit integer, 8 bits character, 4 characters/word
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

  SUBROUTINE STORR(N,IWS,VALUE)
!...STORES A REAL NUMBER IN A WORKSPACE
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION IWS(*)
    INTEGER JWS(2),int(2)
    DOUBLE PRECISION WS
    EQUIVALENCE (WS,JWS),(int,aws)
! move the exact bit pattern from real VALUE to integer IWS(N)
    WS=VALUE
    IWS(N)=JWS(1)
    IWS(N+1)=JWS(2)
!    int=jws
!    write(*,17)value,ws,aws
17  format('storr: ',3(1pe14.6),/10x,4i14)
    RETURN
  END SUBROUTINE STORR

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE LOADR(N,IWS,VALUE)
!...LOADS A REAL NUMBER FROM A WORKSPACE
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION IWS(*)
    INTEGER JWS(2)
    DOUBLE PRECISION WS
    EQUIVALENCE (WS,JWS)
! move the exact bit pattern from integer IWS(N) to real VALUE
    JWS(1)=IWS(N)
    JWS(2)=IWS(N+1)
    VALUE=WS
    RETURN
  END SUBROUTINE LOADR

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE STORRN(N,IWS,ARR)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION IWS(*),ARR(*)
    integer, parameter :: maxr=256
    double precision dlocal(maxr)
    integer ilocal(maxr*nwpr)
    equivalence (dlocal,ilocal)
    if(n.gt.256) then
       write(*,*)'STORRN cannot handle arrays larger than ',maxr
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

  SUBROUTINE LOADRN(N,IWS,ARR)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION IWS(*),ARR(*)
    integer, parameter :: maxr=256
    double precision dlocal(maxr)
    integer ilocal(maxr*nwpr)
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

  SUBROUTINE STORR1(ARR,VAL)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    ARR=VAL
    RETURN
  END SUBROUTINE STORR1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  SUBROUTINE LOADR1(ARR,VAL)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    VAL=ARR
    RETURN
  END SUBROUTINE LOADR1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\


END MODULE METLIB
