! definition of xml elements and attributes for XTDB files
!
! Some default values for the XTDB file, can be changed by user or read_xtdb
! These are also set in pmon6 when NEW Y command
  character (len=8), parameter :: XTDBversion='0.0.3   '
  character (len=8) :: lowtdef    ='298.15  '
  character (len=8) :: hightdef   ='6000    '
  character (len=64) :: bibrefdef  ='U.N. Known  '
  character (len=16) :: eldef     ='VA /-'
  logical :: unary1991=.TRUE., includemodels=.FALSE.
!
! Number of XML elements (tags)
  integer, parameter :: nxmlel=36
! Number of attributes per element
  integer, parameter :: mxmldf=9,mxmlel=5,mxmlsp=4,mxmltp=4,&
       mxmltr=2,mxmlph=3,mxmlsl=3,mxmlcc=3,mxmldp=6,mxmlap=1,&
       mxmlpp=5,mxmlp2=7,mxmlca=1,mxmlsc=2,mxmlbb=0,mxmlbi=4,&
       mxmlmd=0,mxmlmm=5,mxmlvm=5,mxmlpm=2,mxmlem=3,mxmll2=4,&
       mxmleec=1,mxmltm=3,mxmlkm=3,mxmlas=1,mxmlbs=2,mxmlts=1,&
       mxmlus=0,mxmlxx=4,mxmlcy=3,mxmlsys=2
! All XML elements (not using this array when writing XML except 19 and 29)
  character (len=17), dimension(nxmlel), parameter :: xmlel=&
!        123456789.123456789.123456789.123456789.123456789.123456789.
       ['AmendPhase       ',&
        'AssessedSystem   ',&
        'Bibitem          ',&
        'Bibliography     ',&
        'Binary           ',&
        'ConstArray       ',&
        'Constituents     ',&
        'CrystallStructure',&
        'DatabaseInfo     ',&
        'Defaults         ',& ! 10
        'DisorderedPart   ',& ! Same as TC DISORDERED_PART and NEVER model
        'EEC              ',&
        'Element          ',&
        'Einstein         ',&
        'KohlerTernary    ',&
        'Liquid2state     ',&
        'Magnetic         ',&
        'Models           ',&
        'Parameter        ',& ! 19 This indx used in GTP3E for output of XTDB
        'Parameter2       ',& ! 20
        'Permutations     ',&
        'Phase            ',&
        'Sites            ',&
        'Species          ',& ! replaces 'Sublattices      ',&
        'SublConst        ',&
        'System           ',&
        'Ternary          ',&
        'ToopTernary      ',& !28
        'TPfun            ',& !29 used in gtp3E.F90
        'Trange           ',& !30
        'UnAssSystem      ',&
        'Volume           ',&
        'XTDB             ',&
        '                 ',&
        '                 ',&
        '                 ']  !36
! Attributes in alphabetical order of the elements
! 1 AmendPhase attributes
! The model attribute can have several model Id separateb by spaces
  character (len=8),dimension(mxmlap), parameter :: xmlapat=&
        ['Model  ']
! 2 Assessed system attributes.  It has nested Binary and Ternary elements
  character (len=8),dimension(mxmlas), parameter :: xmlasat=&
       ['Software']
! 3 Bibitem attributes.  They provide reference to parameters and models
  character (len=8),dimension(mxmlbi), parameter :: xmlbiat=&
       ['Id      ','Text    ','Date    ','Sign    ']
!        12345678---12345678...12345678---12345678---12345678
! 4 Bibliography has no attributes contains only Bibitem elements
!  character (len=8),dimension(mxmlbb), parameter :: xmlbbat=&
! 5 Binary attributes.  The value of system is two elements joined by a hyphen
! The Commands is a text with commands for the appropriate software
  character (len=8),dimension(mxmlbs), parameter :: xmlbsat=&
       ['System  ','Commands']
! 6 ConstArray attributes.  Is has also one or several mested SublConst elements
  character (len=8),dimension(mxmlca), parameter :: xmlcaat=['Degree  ']
! 7 Constituents attributes (used inside Phase element)
  character (len=16),dimension(mxmlcc), parameter :: xmlccat=&
       ['Sublattice      ','List            ','                ']
! 8 Crystallography attributes (used inside Phase element)
  character (len=16),dimension(mxmlcy), parameter :: xmlcyat=&
       ['StructurBerict  ','WyckoffPositions','                ']
! 9 DatabaseInfo has no attributes
! 10 Defaults attributes
  character (len=16), dimension(mxmldf), parameter :: xmldfat=&
       ['LowT            ','HighT           ','Bibref          ',&
        'Elements        ','Model           ','                ',&
        '                ','                ','                ']
!        123456789.123456...123456789.123456---123456789.123456
! 11 DisorderedPart attributes. Subtract="Y" menas subtract ordered as disorderd
! If no subtract then the NEVER model
  character (len=16),dimension(mxmldp), parameter :: xmldpat=&
       ['Disordered      ','Subtract        ','Sum             ',&
        '                ','                ','                ']
!        123456789.123456...123456789.123456---123456789.123456
!================================ 10
! 12 EEC attributes, Id="EEC"
  character (len=8),dimension(mxmleec), parameter :: xmleecat=&
       ['Bibref  ']
!        12345678---12345678
! 13 Element attributes
  character (len=8), dimension(mxmlel), parameter :: xmlelat=&
       ['Id      ','Refstate','Mass    ','H298    ','S298    ']
!        12345678...12345678---12345678---12345678---12345678
! 14 Einstein attributes Id="GLOWTEIN"
  character (len=8),dimension(mxmlem), parameter :: xmlemat=&
       ['Id      ','MPID    ','Bibref  ']
!        12345678---12345678
! 15 Kohler attributes
  character (len=8),dimension(mxmlkm), parameter :: xmlkmat=&
       ['Phase   ','Constit ','Bibref  ']
! 16 Liquid2state attributes 
  character (len=8),dimension(mxmll2), parameter :: xmll2at=&
       ['Id      ','MPID1   ','MPID2   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
! 17 Magnetic model attributes Id="IHJBCC" or IHJREST or IHJQX
  character (len=8),dimension(mxmlmm), parameter :: xmlmmat=&
       ['Id      ','MPID1   ','MPID2   ','MPID3   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
! 18 Models have no attributes
!  character (len=8),dimension(mxmlmd), parameter :: xmlmdat=&
!       '                ','                ','                ',
! 19 Parameter attributes, Id is TDB format
  character (len=8),dimension(mxmlpp), parameter :: xmlppat=&
       ['Id      ','LowT    ','Expr    ','HighT   ','Bibref  ']
! 20 Parameter2 attributes
  character (len=8),dimension(mxmlp2), parameter :: xmlp2at=&
       ['Id      ','MPID    ','Phase   ','LowT    ','Expr    ',&
        'HighT   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
!        123456789.12---123456789.12
! 21 Permutation attributes 
! Id can be FCC4PERM or BCC4PERM in an AmendPhase element
  character (len=8),dimension(mxmlpm), parameter :: xmlpmat=&
       ['Id      ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
!=================================================== 20
! 22 Phase attributes
  character (len=16),dimension(mxmlph), parameter :: xmlphat=&
       ['Id              ','Configuration   ','State           ']
!        123456789.123456...123456789.123456---123456789.123456
! 23 Sites attributes (nested in the Phase element), replaces Sublattices
  character (len=16),dimension(mxmlsl), parameter :: xmlslat=&
       ['NumberOf        ','Multiplicities  ','                ']
! 24 Species attributes
  character (len=16), dimension(mxmlsp), parameter :: xmlspat=&
       ['Id              ','Stoichiometry   ','MQMQA           ',&
        'UNIQUAC         ']
! Sublattice/Sites attributes (nested in the Phase element) replaced by Sites
!  character (len=16),dimension(mxmlsl), parameter :: xmlslat=&
!       ['NumberOf        ','Multiplicities  ','                ']
! 25 SublConstituent attribute (nested in the parameter2 element)
  character (len=12),dimension(mxmlsc), parameter :: xmlscat=&
       ['Sublattice  ','Species     ']
! 26 System attributes to organize Parameters in unary, binary etc 
  character (len=8),dimension(mxmlsys), parameter :: xmlsysat=&
       ['Elements','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
! 27 Ternary attributes.
  character (len=8),dimension(mxmlts), parameter :: xmltsat=&
       ['System  ']
! 28 Toop attributes.  The Toop constituent should be the first one!
  character (len=8),dimension(mxmltm), parameter :: xmltmat=&
       ['Phase   ','Constit ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
! 29 TPfun attributes
  character (len=8),dimension(mxmltp), parameter :: xmltpat=&
       ['Id      ','LowT    ','Expr    ','HighT   ']
!        12345678...12345678---12345678---12345678
! 30 Trange attributes
  character (len=8),dimension(mxmltr), parameter :: xmltrat=&
       ['Expr    ','HighT   ']
!        12345678...12345678
! 31 UnAssSystem attributes are NONE
! 32 Volume
! 33 XTDB
! 34 -
! 35 -
! 36 -
!================================================ 30
! Attributes for Volume model
  character (len=8),dimension(mxmlvm), parameter :: xmlvmat=&
       ['Id      ','MPID1   ','MPID2   ','MPID3   ','Bibref  ']
! Attributes for XTDB
  character (len=8),dimension(mxmlxx), parameter :: xmlxxat=&
       ['Version ','Software','Date    ','Sign    ']
! 8      12345678---12345678...12345678---12345678---12345678---12345678
! 16     123456789.123456---123456789.123456
! 32     123456789.123456789.123456789.12---123456789.123456789.123456789.12
!=========================================================
!
! Model parameter id
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
