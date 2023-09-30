! definition of xml elements and attributes for XTDB files
!
! Some default values for the XTDB file, can be changed by user or read_xtdb
! These are also set in pmon6 when NEW Y command
  character (len=8), parameter :: XTDBversion='0.0.3   '
  character (len=8) :: lowtdef    ='298.15  '
  character (len=8) :: hightdef   ='6000    '
  character (len=64) :: refiddef  ='U.N. Known  '
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
       mxmlus=0,mxmlxx=4,mxmlcy=3
! All XML elements (not using this array when writing XML except 19 and 29)
  character (len=16), dimension(nxmlel), parameter :: xmlel=&
       ['AmendPhase      ',&
        'AssessedSystem  ',&
        'Bibitem         ',&
        'Bibliography    ',&
        'Binary          ',&
        'ConstArray      ',&
        'Constituents    ',&
        'Crystallograpy  ',&
        'DatabaseInfo    ',&
        'Defaults        ',& ! 10
        'DisorderdPart   ',& ! DisorderedPart is idetical to SplitPhase
        'EEC             ',&
        'Element         ',&
        'Einstein        ',&
        'KohlerTernary   ',&
        'Liquid2state    ',&
        'Magnetic        ',&
        'Models          ',&
        'Parameter       ',& ! 19
        'Parameter2      ',& ! 20
        'Permutations    ',&
        'Phase           ',&
        'Species         ',&
        'SplitPhase      ',& ! Splitphase is identical to DisorderedPart
        'Sublattices     ',&
        'SublConstituent ',&
        'Ternary         ',&
        'ToopTernary     ',&
        'TPfun           ',& !29
        'Trange          ',& !30
        'UnAssSystem     ',&
        'Volume          ',&
        'XTDB            ',&
        '                ',&
        '                ',&
        '                ']  !36
! Attributes in alphabetical order of the elements
! AmendPhase attributes
! The model attribute can have several model Id separateb by spaces
  character (len=8),dimension(mxmlap), parameter :: xmlapat=&
        ['Model  ']
! Assessed system attributes.  It has nested Binary and Ternary elements
  character (len=8),dimension(mxmlas), parameter :: xmlasat=&
       ['Software']
! Bibitem attributes.  They provide reference to parameters and models
  character (len=8),dimension(mxmlbi), parameter :: xmlbiat=&
       ['Id      ','Text    ','Date    ','Sign    ']
!        12345678---12345678...12345678---12345678---12345678
! Bibliography has no attributes contains only Bibitem elements
!  character (len=8),dimension(mxmlbb), parameter :: xmlbbat=&
! Binary attributes.  The system has two elements joined by hyphen as value.
! The Commands is a text with commands for the appropriate software
  character (len=8),dimension(mxmlbs), parameter :: xmlbsat=&
       ['System  ','Commands']
! ConstArray attributes.  Is has also one or several mested SublConst elements
  character (len=8),dimension(mxmlca), parameter :: xmlcaat=['Degree  ']
! Constituents attributes (used inside Phase element)
  character (len=16),dimension(mxmlcc), parameter :: xmlccat=&
       ['Sublattice      ','List            ','                ']
! Crystallography attributes (used inside Phase element)
  character (len=16),dimension(mxmlcy), parameter :: xmlcyat=&
       ['StructurBerict  ','WyckoffPositions','                ']
! DatabaseInfo has no attributes
! Defaults attributes
  character (len=16), dimension(mxmldf), parameter :: xmldfat=&
       ['LowT            ','HighT           ','Refid           ',&
        'Elements        ','Model           ','                ',&
        '                ','                ','                ']
!        123456789.123456...123456789.123456---123456789.123456
! DisorderedPart attributes. This element is identical to the SplitPhase element
  character (len=16),dimension(mxmldp), parameter :: xmldpat=&
       ['Disordered      ','Sum             ','                ',&
        '                ','                ','                ']
!================================ 10
! EEC attributes, identifier EEC
  character (len=8),dimension(mxmleec), parameter :: xmleecat=&
       ['Refid   ']
! Element attributes
  character (len=8), dimension(mxmlel), parameter :: xmlelat=&
       ['Id      ','Refstate','Mass    ','H298    ','S298    ']
!        12345678...12345678---12345678---12345678---12345678
! Einstein attributes
  character (len=8),dimension(mxmlem), parameter :: xmlemat=&
       ['Id      ','MPID    ','Refid   ']
!        12345678---12345678
! Kohler attributes
  character (len=8),dimension(mxmlkm), parameter :: xmlkmat=&
       ['Phase   ','Constit ','Refid   ']
! Liquid2state attributes 
  character (len=8),dimension(mxmll2), parameter :: xmll2at=&
       ['Id      ','MPID1   ','MPID2   ','Refid   ']
!        12345678---12345678...12345678---12345678---12345678
! Magnetic model attributes
  character (len=8),dimension(mxmlmm), parameter :: xmlmmat=&
       ['Id      ','MPID1   ','MPID2   ','MPID3   ','Refid   ']
!        12345678---12345678...12345678---12345678---12345678
! Models have no attributes
!  character (len=8),dimension(mxmlmd), parameter :: xmlmdat=&
!       '                ','                ','                ',
! Parameter attributes
  character (len=8),dimension(mxmlpp), parameter :: xmlppat=&
       ['Id      ','LowT    ','Expr    ','HighT   ','Bibref  ']
! Parameter2 attributes
  character (len=8),dimension(mxmlp2), parameter :: xmlp2at=&
       ['Id      ','MPID    ','Phase   ','LowT    ','Expr    ',&
        'HighT   ','Bibref  ']
!        12345678---12345678...12345678---12345678---12345678
!        123456789.12---123456789.12
! Permutation attributes NEW DESIGN
! Id can be FCC4Perm or BCC4Perm in an AmendPhase element
  character (len=8),dimension(mxmlpm), parameter :: xmlpmat=&
       ['Id      ','Refid   ']
!        12345678---12345678...12345678---12345678---12345678
!=================================================== 20
! Phase attributes
  character (len=16),dimension(mxmlph), parameter :: xmlphat=&
       ['Id              ','Configuration   ','State           ']
!        123456789.123456...123456789.123456---123456789.123456
! Species attributes
  character (len=16), dimension(mxmlsp), parameter :: xmlspat=&
       ['Id              ','Stoichiometry   ','MQMQA           ',&
        'UNIQUAC         ']
! Splitphase attributes are identical to DisorderedPart element
! Sublattice attributes (nested in the Phase element)
  character (len=16),dimension(mxmlsl), parameter :: xmlslat=&
       ['NumberOf        ','Multiplicities  ','                ']
! SublConstituent attribute (nested in the parameter2 element)
  character (len=12),dimension(mxmlsc), parameter :: xmlscat=&
       ['Sublattice  ','Species     ']
! Ternary attributes.
  character (len=8),dimension(mxmlts), parameter :: xmltsat=&
       ['System  ']
! Toop attributes.  The Toop constituent should be the first one!
  character (len=8),dimension(mxmltm), parameter :: xmltmat=&
       ['Phase   ','Constit ','Refid   ']
!        12345678---12345678...12345678---12345678---12345678
! TPfun attributes
  character (len=8),dimension(mxmltp), parameter :: xmltpat=&
       ['Id      ','LowT    ','Expr    ','HighT   ']
!        12345678...12345678---12345678---12345678
! Trange attributes
  character (len=8),dimension(mxmltr), parameter :: xmltrat=&
       ['Id      ','HighT   ']
!        12345678...12345678
! UnAssSystem attributes are NONE
!================================================ 30
! Attributes for Volume model
  character (len=8),dimension(mxmlvm), parameter :: xmlvmat=&
       ['Id      ','MPID1   ','MPID2   ','MPID3   ','Refid   ']
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
! IHJBCC    1    Inden-Hillert-Jarl BCC (-1) an REST (-3)
! IHJQX     2    Inden-Hillert-Jarl-Qing-Xiong 
! GLOWTEIN  4    Einsten low T 
! LIQ2STATE 5    Liquid 2 state model
! VOLOWP    7    Low P volome model
!
! These have no parameters and treated in another way
! SPLITPHASE or as I prefer DISORDEREDPART
! FCC4PERM  FCC symmetric tetrahedron permutations
! BCC4PERM  BCC asymmetric tetrahedron permutations
! EEC       Equi Entropy Criterion
! EBEF      Effective Bond Energy Formalism (including split phase)
!
!   
