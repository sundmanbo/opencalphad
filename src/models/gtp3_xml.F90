!***************************************************************
! General Thermodynamic Package (GTP)
! for thermodynamic modelling and calculations
!
! MODULE GENERAL_THERMODYNAMIC_PACKAGE
!
! Copyright 2011-2022, Bo Sundman, France
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
! contact person: bo.sundman@gmail.com
!
!-----------------------------------------------------------------------
!
! for known unfinished/unchecked bugs and parallelization problems
! look for BEWARE
!
!-----------------------------------------------------------------------
!
! Description of the XTDB data structure
! definition of xml elements and attributes for XTDB files
!
! Versions
! 2024.11.19 Begun revise previous XTDB structure using TYPE
! 2025.01.13 Integrate with OC
!
! Some default values for the XTDB file, can be changed by user or read_xtdb
! These are also set in pmon6 when NEW Y command
  character (len=8), parameter :: XTDBversion='0.1.14   '
  character (len=8) :: lowtdef    ='298.15  '
  character (len=8) :: hightdef   ='6000    '
  character (len=64) :: bibrefdef  ='U.N. Known  '
  character (len=16) :: eldef     ='VA /-'
!  character (len=52) :: ModelAppendXTDB='C:\Users\bosun\Documents\OCHOME\ModelAppendXTDB.XTDB'
  character (len=20) :: ModelAppendXTDB='.\ModelOCAppend.XTDB'
  logical :: unary1991=.TRUE., includemodels=.FALSE.
  integer xtdberr
! this is set TRUE when MQMQA quads entered as species. If constituents .false.
  logical :: xtdbmqmqa=.true.
!
! Number of XDB tags 
  integer, parameter :: nxtdbtags=36
! contains all tag names in no particular order
! if tags extend 18 characters changes may be needed in gtp3EX.F90
  character (len=18), dimension(nxtdbtags), parameter :: xtdbtags=&
!        123456789.123456789.123456789.123456789.123456789.123456789.
       ['XTDB              ',&
        'Defaults          ',&
        'DatabaseInfo      ',&
        'AppendXTDB        ',& ! path to extra TDB file
        'Element           ',&
        'Species           ',& 
        'TPfun             ',& 
        'Trange            ',& 
! 8 above ------------------------------------------
        'Phase             ',& ! With several nested tags
        'Sublattices       ',& 
        'Constituents      ',&
        'CrystalStructure  ',&
        'AmendPhase        ',& ! can have models as magnetic etc
        'Appendix          ',& ! Surround tags in an AppendXTDB file 
        'DisorderedPart    ',& ! as TC DISORDERED_PART and/or NEVER model
        '                  ',& ! chanded Disordered_3Part to attribute
! 16 above end of phase tags------------------
        'Parameter         ',& ! if moved edit xmlpartag in gtp3EX.F90
        'Parameter2        ',& ! maybe never implemented in OC, need more tags
        'Bibliography      ',&
        'Bibitem           ',& ! inside Bibliography
! 20 above Models ----------------------------
        'Models            ',& ! With model tags following
        'Magnetic          ',& ! The modelss Model Parameter Identifiers MPID
        'Einstein          ',&
        'Liquid2state      ',&
        'Volume            ',&
        'EEC               ',& !
        'TernaryXpol       ',& ! Ternary extrapolation tag
! 27 above I think EBEF is not needed as a model, it is defined the parameters
        'UnarySystems      ',& ! These tags are optional for arrangeing data
        'BinarySystem      ',&
        'TernarySystem     ',&
! 30 above, ------------------- add new tags below
        '                  ',& ! Free
        '                  ',&
        '                  ',&
        '                  ',&
        '                  ',&
        '                  ']  !36
!------------------------------------
!
! All tag attributes are defined below to be easy to modify
! They are in the order of the tags
! XTDB tag attributes
!  character (len=18), dimension(nxmltags), parameter :: xmltags=&
  integer, parameter :: nxtdbatt=4
  character (len=9), dimension(nxtdbatt), parameter :: xtdbatt=&
          ['Version  ','Software ','Date     ','Signature']
!...........
!  Defaults 2
  integer, parameter :: ndefatt=9
  character (len=18), dimension(ndefatt), parameter :: defatt=&
       ['LowT             ','HighT            ','Bibref           ',&
        'Elements         ','DefaultModels    ','EEC              ',&
        '                 ','                 ','                 ']
!.............
! DatabaseInfo 3
  integer, parameter :: ninfoatt=3
  character (len=16), dimension(ninfoatt), parameter :: infoatt=&
       ['Software        ','Date            ','Signature       ']
!        123456789.123456...123456789.123456---123456789.123456
!...........
! AppendXTDB 4
  integer, parameter :: nappatt=5,lenappatt=16
  character (len=lenappatt), dimension(nappatt), parameter :: appatt=&
       ['Models          ','Parameters      ','TPfuns          ',&
        'Bibligraphy     ','Miscellaneous   '] 
!        123456789.123456...123456789.123456---123456789.123456
!............
! Element 5
  integer, parameter :: nelatt=5
  character (len=8), dimension(nelatt), parameter :: elatt=&
       ['Id      ','Refstate','Mass    ','H298    ','S298    ']
!        12345678...12345678---12345678---12345678---12345678
!................
! Species 6
  integer, parameter :: nspatt=4
  character (len=16), dimension(nspatt), parameter :: spatt=&
       ['Id              ','Stoichiometry   ','MQMQA           ',&
        'UNIQUAC         ']
!................
! Tpfun attributes 7
  integer, parameter :: ntpatt=4
  character (len=8),dimension(ntpatt), parameter :: tpatt=&
       ['Id      ','LowT    ','Expr    ','HighT   ']
!        12345678...12345678---12345678---12345678
! Trange attributes 8
  integer, parameter :: ntratt=2
  character (len=8),dimension(ntratt), parameter :: tratt=&
       ['Expr    ','HighT   ']
!.............
! Phase attributes 9
  integer, parameter :: nphatt=3
  character (len=16),dimension(nphatt), parameter :: phatt=&
       ['Id              ','Configuration   ','State           ']
!        123456789.123456...123456789.123456---123456789.123456
!.............
! Sublattice attributes (nested in the Phase element) 10
  integer, parameter :: nsubatt=3     
  character (len=16),dimension(nsubatt), parameter :: subatt=&
       ['NumberOf        ','Multiplicities  ','WyckoffPosition ']
!        123456789.123456...123456789.123456---123456789.123456
!.............
! Constituents attributes (used inside Phase element) maybe add NN index 11
  integer, parameter :: nconatt=2
  character (len=16),dimension(nconatt), parameter :: conatt=&
       ['Sublattice      ','List            ']
!...............
! CrystalStructure attributes (used inside Phase element) 12
  integer, parameter :: ncrystatt=4
  character (len=16),dimension(ncrystatt), parameter :: crystatt=&
       ['Prototype       ','PearsonSymbol   ','SpaceGroup      ',&
        'StructurBericht ']
!.............
! AmendPhase attributes 13
  integer, parameter :: namphatt=2
  character (len=12),dimension(namphatt), parameter :: amphatt=&
        ['Model       ','Permutation ']
!         123456789.12---123456789.12---123456789.12
!..............
! Appendix tag 14
! Begin and end tag in an appended XTDB file (except Bibliograpy)
!  integer, parameter :: npermatt=0
!  character (len=8),dimension(npermatt), parameter :: permatt=&
!...............
! DisorderedPart tag both NEVER model and Disordered_3Parts with Subtract 15
  integer, parameter :: ndis=3
  character (len=12),dimension(ndis), parameter :: disatt=&
       ['Disordered  ','Sum         ','Subtract    ']
!        123456789.12---123456789.12---123456789.12
!...
! Unused 16
!  integer, parameter :: ndis3=2
!  character (len=12),dimension(ndis3), parameter :: dis3att=&
!................
! Parameter attributes, Id is as in TDB files 17
  integer, parameter :: npar=5
  character (len=8),dimension(npar), parameter :: paratt=&
       ['Id      ','LowT    ','Expr    ','HighT   ','Bibref  ']
!..........
! Parameter2 attributes (not supported by OC) 18
  integer, parameter :: npar2=7
  character (len=8),dimension(npar2), parameter :: par2att=&
       ['Id      ','MPID    ','Phase   ','LowT    ','Expr    ',&
        'HighT   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
!................
! Bibliography has no attributes contains only Bibitem elements
!................
! Bibitem attributes.  They provide reference to parameters and models 19
  integer, parameter :: nbibatt=4
  character (len=8),dimension(nbibatt), parameter :: bibattt=&
       ['Id      ','Text    ','Date    ','Sign    ']
!        12345678---12345678...12345678---12345678---12345678
!...................
! UnarySystem 20
  integer, parameter :: nuniatt=2
  character (len=8), dimension(nuniatt), parameter :: uniatt=&
       ['Element ','Bibref  ']
!....................
! BinarySystem attributes.  The Species is two elements separated by a space 21
! The CalcDia attribute is a software depenednt command string
  integer, parameter :: nbinatt=3
  character (len=8),dimension(nbinatt), parameter :: binatt=&
       ['Species ','Bibref  ','CalcDia ']
!        12345678---12345678...12345678---12345678---12345678
!....................
! TernarySystem attributes.  The Species is 3 elements separated by a space 22
! The CalcDia attribute is a software depenednt command string
  integer, parameter :: nteratt=3
  character (len=8),dimension(nteratt), parameter :: teratt=&
       ['Species ','Bibref  ','CalcDia ']
!        12345678---12345678...12345678---12345678---12345678
!================================================================
! Attributes:
!================================================================
! The AmedPhase model attribute has these values
!......................
! Magnetic model attributes Id="IHJBCC" or IHJREST or IHJQX
  integer, parameter :: nmagatt=5
  character (len=8),dimension(nmagatt), parameter :: magatt=&
       ['Id      ','MPID1   ','MPID2   ','MPID3   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
!......................
! Einstein attributes Id="GEIN"
  integer, parameter :: ngeinatt=3
  character (len=8),dimension(ngeinatt), parameter :: geinatt=&
       ['Id      ','MPID    ','Bibref  ']
!        12345678---12345678...12345678
!.....................
! Liquid2state attributes, Id="LIQ2STATE"
  integer, parameter :: nliq2att=4
  character (len=8), dimension(nliq2att), parameter :: liq2att=&
       ['Id      ','MPID1   ','MPID2   ','Bibref  ']
!.......................
! Volume, ID="XGL05"          not implemented in OC
  integer, parameter :: nvolatt=5
  character (len=8),dimension(nvolatt), parameter :: volattt=&
       ['Id      ','MPID1   ','MPID2   ','MPID3   ','Bibref  ']
!        12345678...12345678...12345678...12345678...12345678
!...................  
! EEC attributes, Id="EEC"
  integer, parameter :: neccatt=2
  character (len=8), dimension(neccatt), parameter :: eecatt=&
       ['Id      ','Bibref  ']
!..................
! TernaryXpol attributes.  
  integer, parameter :: nterxpolatt=3
  character (len=8), dimension(nterxpolatt), parameter :: terxpolatt=&
       ['Phase   ','System  ','Xmode   ']
!        12345678---12345678...12345678---12345678---12345678
!....................
!=========================================================
!
! Current list of MPID in OC, related to the models
  integer, parameter :: noofmpid=9
  character (len=8), dimension(noofmpid), parameter :: mpidok=&
       ['G       ','TC      ','BMAG    ','CT      ','NT      ','IBM     ',&
        'LNTH    ','G2      ','L       ']
! 8      12345678---12345678...12345678---12345678---12345678---12345678
! The L is in princple allowed only for excess G parameters but treated as G
! IBM  was intended for element specific magneton number ....
! model    MPID index
! IHJBCC       2  3
! IHJREST      2  3
! IHJQX        3  4  5
! GEIN         7
! LIQUD2STATE  7  8
! 
! OLD list of MPID, some may have a constituent/element 
!  character (len=8), dimension(36), parameter :: mpidw=&
!       ['G       ','TC      ','BMAG    ','CT      ','NT      ','IBM     ',&
!        'LNTH    ','V0      ','VA      ','VB      ','VC      ','VS      ',&
!        'MQ      ','MF      ','MG      ','G2      ','THT2    ','DCP2    ',&
!        'LPX     ','LPY     ','LPZ     ','LPTH    ','EC11    ','EC12    ',&
!        'RHO     ','VISC    ','LAMB    ','HMVA    ','TSCH    ','CSCH    ',&
!        '        ','        ','        ','        ','        ','        ']
! 8      12345678---12345678...12345678---12345678---12345678---12345678
!
! The meaning of the model parameters is entered in init_gtp in gtp3A.F90
!
! An attempt to reconcile XTDB handling of models and additions with OC
! Some models/additions has no parameters.  Those which has are listed below.
! - DisorderdPart and EBEF has 2 separate sets of parameters (software)
! - Permutations usually use wildcard parameters (with *) to reduce the
!   number of duplicate model parameters (software)
! - EBEF is the same as DisorderedPart
!
! All parameters i OC has an MPID index, The parameters for the Gibbs energy
! G (or L) has index 1 (one)
!------------------------------------------------
!
! In OC each parameter has an MPID index stored which is summed
! independently and later used to calculate the addition.
!
! XTDB identifier and MPID        OC MPID index and name
! IHJBCC and IHJREST  Inden-Hillert-Jarl magnetic model, AFF=-1 and AFF=-3
!      TC                         2 TC     4          Curie/Neel T
!      BMAGN                      3 BMAG   3          Bohr magneton number
! IHJQX               Inden-Hillert-Jarl-Qing-Xiong magnetic model, AFF=0
!      CT                         4 CTA    4          Curie T
!      NT                         5 NTA    5          Neel T
!      BMAGN                      3 BMAG   3          Aver. Bohr magneton num
! GEIN                Einstein low T vibrational energy                     
!      LNTH                       7 LNTH   2          Einstein T
! LIQ2STATE           Merging amorphous low T phase and liguid
!      LNTH                       7 LNTH   2          Einstein T for amorph.
!      GD                        16 G2     6          Melting energy of amorph
! FCC4PERM
! BCC4PERM
! FCC4PERM
! FCC4PERM
!==================================================
! For the moment we have 9 models ...?
  INTEGER, parameter :: gtp_xtdbcompatibility_version=1
  type gtp_xtdbcompatibility
     character(len=:), allocatable :: modelid
! this character has the MPID used in the xtdb file
     character*8, dimension(:), allocatable :: mpid
     integer nmpid
! this character has the default MPID used by oc
     character*8, dimension(:), allocatable :: ocmpid
     integer, dimension(:), allocatable ::  ocix
  end type gtp_xtdbcompatibility
  type(gtp_xtdbcompatibility), dimension(:), allocatable :: xtdbmodel
  integer, parameter :: nxtdbmpids=9
!
!-------------------------- old below
! Models accepted by OC in the AmendPhase tag
  integer, parameter ::noofmodels=5
  character (len=16), dimension(noofmodels), parameter :: amphmodel=&
! 8      123456789.123456---123456789.123456---123456789.123456
       ['IHJBCC          ','IHJREST         ','IHJQX           ',&
        'GEIN            ','LIQ2STATE       ']
! Permutations accepted by OC in the AmendPhase tag
  integer, parameter :: noofpermut=2
  character (len=16), dimension(noofpermut), parameter :: amphpermut=&
        ['FCC4PERM        ','BCC4PERM        ']
!
! IHJBCC    1    Inden-Hillert-Jarl for BCC with Aff=-1
! IHJREST   1    Inden-Hillert-Jarl for other with Aff=-3
! IHJQX     2    Inden-Hillert-Jarl-Qing-Xiong with Aff=0
! GEIN      4    Einstein low T 
! LIQ2STATE 5    Liquid 2 state model
! VLOWP1    7    Low P volome model according to Lu?
!
! These have no parameters and are treated in another way
! DISORDEREDPART same as TDB file DISORDERED_PART and NEVER
! FCC4PERM  FCC symmetric tetrahedron permutations
! BCC4PERM  BCC asymmetric tetrahedron permutations
! EEC       Equi Entropy Criterion is set by Delfaults
! EBEF      Effective Bond Energy Formalism may use "species@sublattice"
!
!--------------- end of old

!=========================================================
! Predefined functions in TPfuns
  integer, parameter :: predeftpfun=5
  character*8, dimension(predeftpfun), parameter :: nottpfun= &
       ['LN      ','LOG     ','EXP     ','ERF     ','GEIN    ']
! LN and LOG is the same thing, LOG10 is not used. 
! TPfun have these hardcoded in xmlmake
!
!=========================================================   
! There is a need to handle the Model feature of XTDB with different
! software and data structure in applocation software.  The data structures
! here and below is for temporary use reading xtdb file
!========================================================
!
     type const
! list of constituents in each sublattice of a phase, used in phnest
       character (len=:), allocatable :: subx
       character (len=:), allocatable :: list
    end type const
!    type(const), allocatable, target :: constrec
    type phnest
! all data needed to create the phase record before entering parameters
       integer ncon
       character*1 state
       character (len=:), allocatable :: id
       character (len=:), allocatable :: confent
! when reading Noof here one allocates the dimension of clist !!
       character (len=:), allocatable :: Noof
! The mult remain character until the phase recond is allocated
       character (len=:), allocatable :: mult
! for each sublattice clist will be allocated with the constituents!       
! Each array element can have a differnt number of characters!!
       type(const), dimension(:), pointer :: clist 
       character (len=:), allocatable :: crystal
! The model Id is the amendph
       character (len=:), allocatable :: amendph
! this is the attributes from the XTDB file for disordered part
! disordered phase, sublattices to sum and if subtract ordered as disordered
       character (len=:), allocatable :: dispar
    end type phnest
    type(phnest), allocatable :: phrec
!
! Attributes for AppendXTDB files.  The *appy indicate if todo (-1) or done 1
  character*64 modelappx,parappx,tpfunappx,biblioappx,miscappx
  integer modelappy,parappy,tpfunappy,biblioappy,miscappy,allappy
!
  
!------------------------------------------------------------
! global parameter copied from modile xtdblib in xtdbread.F90 itself
!
!
! used positions in attpos ?? line number in current file
  integer attpos,fline
! maximum number of nested tags
  integer, parameter :: maxlevels=10,commenttag=999
! tagnest(level) is negative, -tagno, if more attributes for tagno to read 
  integer, dimension(maxlevels) :: tagnest
! tagnames have max length 18, the endoftag is set to '</tagname>'
  character(len=21), dimension(maxlevels) :: endoftag
!   
! an expression will be concatinated from TPfun/Trange and Parameter/Trange tags
  character(len=:), allocatable :: wholexpr

  character(len=:), allocatable :: cc

! set to .true. when reading subsets of the XTDB file
  logical ignorEOT

! TPfuns are used in parameters. All TPfun for entered parameters must be found
  integer, parameter :: maxtpfun=500
  integer ntp,missingtp,missingbib
! alltpfun are names of all tpfuns missing or entered

! the extracted data for software stored in these records
! these integers are the last entered element, species etc.
!  integer nselel,nselsp,nselph,nselpar,nseltp,nselbib
! nselph already defined in gtp3
  integer nselel,nselsp,nselpar,nseltp,nselbib

  type ocelement
     character*2 elname
     character(len=:), allocatable :: data
  end type ocelement
  type(ocelement), dimension(:), allocatable :: selel
    
  type ocspecies
     character*24 species
     character(len=:), allocatable :: data
     character*2, dimension(:), allocatable :: elnames
     double precision, dimension(:), allocatable :: stoicc
! electric charge
     double precision :: charge
! mqmqa or uniquac
     character(len=:), allocatable :: extra
  end type ocspecies
  type(ocspecies), dimension(:), allocatable :: selsp
! this array will have the selsp indices in alpahetical order
  integer, dimension(:), allocatable :: selspord
    
  type ocphases
     character*24 phasename
     integer nsublat
! in this array only selected constituents are entered
     character(len=:), allocatable :: mult
     character(len=:), allocatable :: const
     character(len=:), allocatable :: confent
     character(len=:), allocatable :: amendph
     character(len=:), allocatable :: dispar
     character(len=:), allocatable :: data
     type(octerxpol), pointer :: terxpol
  end type ocphases
  type(ocphases), dimension(:), allocatable :: selph

  type ocxparam
     character*64 parname
     character(len=:), allocatable :: data
  end type ocxparam
  type(ocxparam), dimension(:), allocatable :: selpar
  
  type octpfun
     character*16 tpfunname
     character(len=:), allocatable :: data
! seltp(*)%tatus is negative if missing
     integer status
  end type octpfun
  type(octpfun), dimension(:), allocatable :: seltpfun
    
  type ocbib
     character*8 bibitem
     character(len=:), allocatable :: data
     integer status
  end type ocbib
  type(ocbib), dimension(:), allocatable :: selbib

  character(len=:), allocatable :: defaultbib

  type octerxpol
! ternary extrapolation linked from the phase
! The ternaryXpol tags with selected constituents are stored in these records.
! If the phase is already selected they are added to the texpol list
! otherwise kept in the firstxpol list until the phase is entered
     character(len=:), allocatable :: phase
     character(len=:), allocatable :: sps
     character(len=:), allocatable :: xpol
     type(octerxpol), pointer :: next
  end type octerxpol
  logical debug

! this is the start of a linked list of ternary extrapolations
! waiting for the phase to be selected.  Typically it contains TernaryXpol
! what are defined inside the Phase tag itself
  type(octerxpol), pointer :: firstxpol,lastxpol,xpol

! When first phase entered we must not enter more elements/species
    logical nomorelements
! for dimensioning these arrays and amount used
    integer maxtdbel, maxtdbsp, maxtdbph, maxpar, maxtp, maxbib

  
