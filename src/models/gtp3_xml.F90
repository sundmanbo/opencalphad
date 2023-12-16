! definition of xml elements and attributes for XTDB files
!
! Some default values for the XTDB file, can be changed by user or read_xtdb
! These are also set in pmon6 when NEW Y command
  character (len=8), parameter :: XTDBversion='0.0.4   '
  character (len=8) :: lowtdef    ='298.15  '
  character (len=8) :: hightdef   ='6000    '
  character (len=64) :: bibrefdef  ='U.N. Known  '
  character (len=16) :: eldef     ='VA /-'
  logical :: unary1991=.TRUE., includemodels=.FALSE.
!
! Number of XML elements (tags)
  integer, parameter :: nxmlel=36
! Number of attributes per element, some has zero
  integer, parameter :: mxml01=1,mxml02=2,mxml03=4,mxml05=3,&
       mxml06=3,mxml07=3,mxml09=9,mxml10=2,mxml11=2,mxml12=2,&
       mxml13=5,mxml14=3,mxml15=2,mxml16=3,mxml17=4,mxml18=5,&
       mxml20=5,mxml21=2,mxml22=3,mxml23=3,mxml24=4,mxml25=2,&
       mxml26=3,mxml27=4,mxml28=2,mxml29=2,mxml30=5,mxml31=4
! All XML elements (not using this array when writing XML except 19 and 29)
  character (len=17), dimension(nxmlel), parameter :: xmlel=&
!        123456789.123456789.123456789.123456789.123456789.123456789.
       ['AmendPhase       ',&
        'AppendiXTDB      ',&
        'Bibitem          ',&
        'Bibliography     ',&
        'BinarySystem     ',&
        'Constituents     ',&
        'CrystalStructure ',&
        'DatabaseInfo     ',&
        'Defaults         ',&
        'Disordered_2Part ',& ! Same as TC DISORDERED_PART and NEVER model
! 10-----------------------------
        'Disordered_3Part ',& ! Same as TC DISORDERED_PART subracting ordered
        'EEC              ',&
        'Element          ',&
        'Einstein         ',&
        'HigherOrderSystem',&
        'KohlerModel      ',&
        'Liquid2state     ',&
        'Magnetic         ',&
        'Models           ',&
        'Parameter        ',& ! 20, index 19 is used in GTP3E for output of XTDB
! 20 ------------------------------------------
        'Permutations     ',&
        'Phase            ',&
        'Sublattices      ',& ! Using Sites not accepted
        'Species          ',& 
        'TernarySystem    ',&
        'ToopModel        ',& !25 was 28
        'TPfun            ',& !26 was 29 used in gtp3E.F90
        'Trange           ',& 
        'UnarySystem      ',&
        'Volume           ',&
! 30 ------------------------------------------
        'XTDB             ',&
        '                 ',&
        '                 ',&
        '                 ',&
        '                 ',&
        '                 ']  !36
!------------------------------------
! tags "Metadata", "writer" and "Database" not included ...
! Attributes in alphabetical order of the elements
! 1 AmendPhase attributes
! The model attribute can have several model Id separateb by spaces
  character (len=8),dimension(mxml01), parameter :: xmlapat=&
        ['Model  ']
! 2 AppendiXTDB attributes.  Extra files for the database
!  character (len=11),dimension(mxml02), parameter :: xmlasat=&
  character (len=11),dimension(mxml02), parameter :: xmlaxat=&
       ['File       ','Description']
! 3 Bibitem attributes.  They provide reference to parameters and models
  character (len=8),dimension(mxml03), parameter :: xmlbiat=&
       ['Id      ','Text    ','Date    ','Sign    ']
!        12345678---12345678...12345678---12345678---12345678
! 4 Bibliography has no attributes contains only Bibitem elements
! 5 BinarySystem attributes.  The system is two elements joined by a hyphen
! The Commands is a text with commands for the appropriate software. Optional
  character (len=8),dimension(mxml05), parameter :: xmlbsat=&
       ['System  ','Commands','Bibref  ']
! x6 ConstArray attributes.  Is has one or several mested SublConst elements
!  character (len=8),dimension(mxmlca), parameter :: xmlcaat=['Degree  ']
! 6 Constituents attributes (used inside Phase element) maybe add NN index
  character (len=16),dimension(mxml06), parameter :: xmlccat=&
       ['Sublattice      ','List            ','                ']
! 7 CrystalStructure attributes (used inside Phase element)
  character (len=16),dimension(mxml07), parameter :: xmlcsat=&
      ['Prototype       ','PearsonSymbol   ','SpaceGroup      ']
!       123456789.123456---123456789.123456...123456789.123456
! 8 DatabaseInfo has no attributes
! 9 Defaults attributes
  character (len=16), dimension(mxml09), parameter :: xmldfat=&
       ['LowT            ','HighT           ','Bibref          ',&
        'Elements        ','DefaultModels   ','                ',&
        '                ','                ','                ']
!        123456789.123456...123456789.123456---123456789.123456
! 10 Disordered_2Part tag (NEVER model)
  character (len=12),dimension(mxml10), parameter :: xmldp2at=&
       ['Disordered  ','Sum         ']
!        123456789.12...123456789.123456---123456789.123456
!================================ 10 above
! 11 Disordered_3Part attributes. 
  character (len=12),dimension(mxml11), parameter :: xmldp3at=&
       ['Disordered  ','Sum         ']
!        123456789.123456...123456789.123456---123456789.123456
! 12 EEC attributes, Id="EEC"
  character (len=8),dimension(mxml12), parameter :: xmleecat=&
       ['Id      ','Bibref  ']
!        12345678---12345678
! 13 Element attributes
  character (len=8), dimension(mxml13), parameter :: xmlelat=&
       ['Id      ','Refstate','Mass    ','H298    ','S298    ']
!        12345678...12345678---12345678---12345678---12345678
! 14 Einstein attributes Id="GLOWTEIN"
  character (len=8),dimension(mxml14), parameter :: xmlemat=&
       ['Id      ','MPID    ','Bibref  ']
!        12345678---12345678...12345678
! 15 HigherOrderSystem optional
  character (len=8),dimension(mxml15), parameter :: xmlhosat=&
       ['System  ','Bibref  ']
!        12345678---12345678...12345678
! 16 Kohler attributes
  character (len=8),dimension(mxml16), parameter :: xmlkmat=&
       ['Phase   ','System  ','Bibref  ']
!        12345678---12345678...12345678
! 17 Liquid2state attributes 
  character (len=8),dimension(mxml17), parameter :: xmll2at=&
       ['Id      ','MPID1   ','MPID2   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
! 18 Magnetic model attributes Id="IHJBCC" or IHJREST or IHJQX
  character (len=8),dimension(mxml18), parameter :: xmlmmat=&
       ['Id      ','MPID1   ','MPID2   ','MPID3   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
! 19 Models have no attributes
! 19 Parameter attributes, Id is TDB format
  character (len=8),dimension(mxml20), parameter :: xmlppat=&
       ['Id      ','LowT    ','Expr    ','HighT   ','Bibref  ']
! x20 Parameter2 attributes not included in XTDB?
!  character (len=8),dimension(mxmlp2), parameter :: xmlp2at=&
!       ['Id      ','MPID    ','Phase   ','LowT    ','Expr    ',&
!        'HighT   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
!        123456789.12---123456789.12
!=================================================== 20
! 21 Permutation model attributes 
! Id can be FCC4PERM or BCC4PERM in an AmendPhase element
  character (len=8),dimension(mxml21), parameter :: xmlpmat=&
       ['Id      ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
! 22 Phase attributes
  character (len=16),dimension(mxml22), parameter :: xmlphat=&
       ['Id              ','Configuration   ','State           ']
!        123456789.123456...123456789.123456---123456789.123456
! 23 Sublattice attributes (nested in the Phase element)
  character (len=16),dimension(mxml23), parameter :: xmlslat=&
       ['NumberOf        ','Multiplicities  ','WyckoffPosition ']
!        123456789.123456...123456789.123456---123456789.123456
! 24 Species attributes
  character (len=16), dimension(mxml24), parameter :: xmlspat=&
       ['Id              ','Stoichiometry   ','MQMQA           ',&
        'UNIQUAC         ']
! 25 TernarySystem optional
  character (len=12),dimension(mxml25), parameter :: xmltsat=&
       ['System  ','Bibref  ']
! 26 Toop attributes.  The Toop constituent should be the first one!
  character (len=8),dimension(mxml26), parameter :: xmltmat=&
       ['Phase   ','System  ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
! 27 TPfun attributes
  character (len=8),dimension(mxml27), parameter :: xmltpat=&
       ['Id      ','LowT    ','Expr    ','HighT   ']
!        12345678...12345678---12345678---12345678
! 28 Trange attributes
  character (len=8),dimension(mxml28), parameter :: xmltrat=&
       ['Expr    ','HighT   ']
!        12345678...12345678
! 29 UnarySystem 
  character (len=8), dimension(mxml29), parameter :: xmlusat=&
       ['Element ','Bibref  ']
! 30 Volume
  character (len=8),dimension(mxml30), parameter :: xmlvmat=&
       ['Id      ','MPID1   ','MPID2   ','MPID3   ','Bibref  ']
!        12345678...12345678
!=========================================== 30
! 31 XTDB tag
  character (len=8),dimension(mxml31), parameter :: xmlxxat=&
       ['Version ','Software','Date    ','Sign    ']
! 8      12345678---12345678...12345678---12345678---12345678---12345678
! 32 -
! 33 -
! 34 -
! 35 -
! 36 -
! Attributes for Volume model
! Attributes for XTDB
! 16     123456789.123456---123456789.123456
! 32     123456789.123456789.123456789.12---123456789.123456789.123456789.12
!=========================================================
!
! Model parameter id for OC
  character (len=8), dimension(36), parameter :: mpidw=&
       ['G       ','TC      ','BMAG    ','CT      ','NT      ','IBM     ',&
        'LNTH    ','V0      ','VA      ','VB      ','VC      ','VS      ',&
        'MQ      ','MF      ','MG      ','G2      ','THT2    ','DCP2    ',&
        'LPX     ','LPY     ','LPZ     ','LPTH    ','EC11    ','EC12    ',&
        'RHO     ','VISC    ','LAMB    ','HMVA    ','TSCH    ','CSCH    ',&
        '        ','        ','        ','        ','        ','        ']
! 8      12345678---12345678...12345678---12345678---12345678---12345678
!
! Addition  OCix identifiers
! IHJBCC    1    Inden-Hillert-Jarl for BCC with Aff=-1
! IHJREST   1    Inden-Hillert-Jarl for other with Aff=-3
! IHJQX     2    Inden-Hillert-Jarl-Qing-Xiong with Aff=0
! GLOWTEIN  4    Einsten low T 
! LIQ2STATE 5    Liquid 2 state model
! VOLOWP    7    Low P volome model according to Lu?
!
! These have no parameters and treated in another way
! DISORDEREDPART same as TDB file DISORDERED_PART and NEVER
! FCC4PERM  FCC symmetric tetrahedron permutations
! BCC4PERM  BCC asymmetric tetrahedron permutations
! EEC       Equi Entropy Criterion
! EBEF      Effective Bond Energy Formalism using "species@sublattice"
!
!   
