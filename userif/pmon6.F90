!
MODULE cmon1oc
!
! Copyright 2012-2019, Bo Sundman, France
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
! Contact persion: Bo.Sundman@gmail.com
!
!------------------------------------------------------------------
!
!************************************
! command line monitor for OC version 3
!************************************
!
  use ocsmp
!  use liboceq
! 
! parallel processing, set in gtp3.F90
!$  use omp_lib
!
  implicit none
!
! option record
  TYPE ocoptions
! unit for listing, default is kou (screen)
     integer lut
  end TYPE ocoptions
  type(ocoptions) :: optionsset
!
contains
!
  subroutine oc_command_monitor(version,linkdate,narg,argline)
! command monitor
    implicit none
!
! linkdat is date when program was linked
! argline and narg are inline arguments
    character linkdate*(*),version*(*),argline(*)*64
    integer narg
! various symbols and texts
    character :: ocprompt*8='--->OC5:'
    character name1*24,name2*24,line*80,model*72,chshort
    integer, parameter :: ocmonversion=40
! for the on-line help, at present turn off by default, if a HTML file set TRUE
    character*128 browser,latexfile,htmlfile
    logical :: htmlhelp=.FALSE.
!    logical :: htmlhelp=.TRUE.
! element symbol and array of element symbols for database use
    character elsym*2,ellist(maxel)*2
! more texts for various purposes
    character text*72,string*256,ch1*1,chz*1,selection*27,funstring*1024
    character axplot(2)*24,axplotdef(2)*24,quest*20
    character longstring*2048,optres*40
! this is now declared in metlib package
!    character workingdir*128
! separate file names for remembering and providing a default
    character ocmfile*128,ocufile*128,tdbfile*128,ocdfile*128,filename*128
! home for OC and default directory for databases
    character ochome*64,ocbase*64
! prefix and suffix for composition sets
    character prefix*4,suffix*4
! element mass
    double precision mass
! constituent fractions of a phase
    double precision, dimension(maxconst) :: yarr
! stoichiometry of a specis and sublattice sites of a phase
    double precision, dimension(maxsubl) :: stoik
! calculated vaules of a function (G, G.T, G.P, G.T.T; G.T.P and G.P.P)
    double precision val(6)
! estimated chemical potentials after a grid minimization and TP for ref states
    double precision cmu(maxel),tpa(2)
! cpu time measurements
    double precision ending,starting
!>>>> has to be reorganized ------------------------------------
! axis variables and limits
! default values used for axis variables
    double precision dinc,dmin,dmax
! plot ranges, texts and defaults for graphics
    type(graphics_options) :: graphopt
    integer grunit
! path to start directory stored inside metlib!!
!    character macropath*128
! plot texts
!    type(graphics_textlabel), allocatable, target :: textlabel
    type(graphics_textlabel), pointer :: textlabel
    type(graphics_textlabel), pointer :: labelp
! axis data structures
    type(map_axis), dimension(5) :: axarr
! if more than one start equilibrium these are linked using the ceq%next index
    type(gtp_equilibrium_data), pointer :: starteq
! for map results
    type(map_node), pointer :: maptop,mapnode,maptopsave,maptopcheck
!    type(map_line) :: mapline
! seqxyz has initial values of seqx, seqy and seqz
    integer noofaxis,noofstarteq,seqxyz(3)
! this should be removed
!    TYPE(ssm_node), pointer :: resultlist
!<<<<<<<--------------------------------------------------------------
! used for element data and R*T
    double precision h298,s298,rgast
! temporary reals
    double precision xxx,xxy,xxz,totam,cpham
! input data for grid minimizer
    double precision, dimension(maxel) :: xknown,aphl
! arrays for grid minimization results
    integer, dimension(maxel) :: iphl,icsl,nyphl
! selected kommand and subcommands
    integer kom,kom2,kom3,kom4
! selected output mode for results and the default, list output unit lut
    integer listresopt,lrodef,lut,afo
! integers used for elements, phases, composition sets, equilibria, defaults
    integer iel,iph,ics,ieq,idef
! for gradients in MU and interdiffusivities
    integer nend
! dimension of mugrad for 16x16 matrix 
    double precision mugrad(300),mobilities(20)
!-------------------
! selection of minimizer and optimizer
    integer minimizer,optimizer
! plot unit for experimental data used in enter many_equilibria
    integer :: plotdataunit(9)=0,plotunit0=0
! temporary integer variables in loops etc
    integer i1,i2,j1,j2,iax,threads
! more temporary integers
    integer jp,kl,svss,language,last,leak,j3
! and more temporary integers
    integer ll,lokcs,lokph,lokres,loksp,lrot,maxax
! and more temporary integers
    integer mode,ndl,neqdef,noelx,nofc,nopl,nops,nv,nystat,times
! temporary matrix
!    double precision latpos(3,3)
! used to call init_gtp for the NEW command
    integer intv(10)
    double precision dblv(10)
!-------------------
! variables for lmdif
!    integer, parameter :: lwam=2500
    integer :: lwam=2500,nfev
    integer :: nopt1=100, mexp=0,nvcoeff=0,nopt,iflag,mexpdone=0,nvcoeffdone=0
    integer, allocatable, dimension(:) :: iwam
    double precision, allocatable, dimension(:) :: wam
! tccovar is the covariance matrix used to calculate RSD as in Thermo-Calc
    double precision, allocatable, dimension(:,:) :: fjac,cov1,cormat,tccovar
    double precision :: optacc=1.0D-3
    logical :: updatemexp=.true.
! saved parameters for analyze
    double precision, allocatable, dimension(:,:) :: savedcoeff
    double precision savesumerr,delta
    integer analyze,cormatix,nvcoeffsave,mexpsave,iz
! this is least square error from using LMDIF
! 1: previous value, 2 new value, 3 normalized error (divided by m-n)
    double precision err0(3)
! occational segmentation fault when deallocating www ....
!    double precision, dimension(maxw) :: www
!    double precision, dimension(:), allocatable :: www
    double precision, dimension(:), allocatable :: coefs
    double precision, dimension(:), allocatable :: errs
!    external new_assessment_calfun
!-------------------
! loop variable when entering constituents of a phase
    integer icon,flc
! array with constituents in sublattices when entering a phase
!    character, dimension(maxconst) :: const*24
! for macro and logfile and repeating questions
    logical logok,stop_on_error,once,wildcard,twice,startupmacro,temporary
    logical listzero
! unit for logfile input, 0 means no logfile
    integer logfil
! remember default for calculate phase
    integer defcp
! for state variables as conditions
    integer istv
    double precision coeffs(10),textfontscale
    TYPE(gtp_state_variable), target :: stvrvar
    TYPE(gtp_state_variable), pointer :: stvr
!    TYPE(gtp_state_variable), dimension(10) :: stvarr
    TYPE(gtp_condition), pointer :: pcond,firstc
! current equilibrium records
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq
    TYPE(gtp_phase_varres), pointer :: parres
! addition record used for listing calculated values
    type(gtp_phase_add), pointer :: addrec
!
    character actual_arg(2)*16
    character cline*128,option*80,aline*128,plotfile*128,eqname*24
! variable phase tuple
    type(gtp_phasetuple), pointer :: phtup
!----------------------------------------------------------------
! here are all commands and subcommands
!    character (len=64), dimension(6) :: oplist
    integer, parameter :: ncbas=30,nclist=21,ncalc=12,ncent=21,ncread=6
    integer, parameter :: ncam1=18,ncset=27,ncadv=15,ncstat=6,ncdebug=9
    integer, parameter :: nselect=6,nlform=6,noptopt=9,nsetbit=6
    integer, parameter :: ncamph=18,naddph=12,nclph=6,nccph=6,nrej=9,nsetph=6
    integer, parameter :: nsetphbits=15,ncsave=6,nplt=21,nstepop=6
    integer, parameter :: ninf=9
! basic commands
    character (len=16), dimension(ncbas), parameter :: cbas=&
       ['AMEND           ','CALCULATE       ','SET             ',&
        'ENTER           ','EXIT            ','LIST            ',&
        'QUIT            ','READ            ','SAVE            ',&
        'HELP            ','INFORMATION     ','BACK            ',&
        'NEW             ','MACRO           ','ABOUT           ',&
        'DEBUG           ','SELECT          ','DELETE          ',&
        'STEP            ','MAP             ','PLOT            ',&
        'HPCALC          ','FIN             ','OPTIMIZE        ',&
        'SHOW            ','ANALYZE_ASSESSMT','                ',&
        '                ','                ','                ']
! in French
!        'MODIFIEZ        ','CALCULEZ        ','REGLEZ          ',&
!        'ENTREZ          ','EXIT            ','AFFICHER        ',&
!        'QUIT            ','LIRE            ','SAUVGARDE       ',&
!        'AIDEZ           ','INFORMATION     ','RETURNEZ        ',&
!        'NOUVEAU         ','MACRO           ','ABOUT           ',&
!        'DEBUG           ','SELECTIONEZ     ','EFFACEZ         ',&
!        'STEP            ','MAP             ','DESSINEZ        ',&
!        'HPCALC          ','FIN             ','                ']
! NOTE a command line can contain options preceded by /
! for example "list /out=myfile.dat all_data" or
!-------------------
! subcommands to LIST
    character (len=16), dimension(nclist) :: clist=&
         ['DATA            ','SHORT           ','PHASE           ',&
         'STATE_VARIABLES ','BIBLIOGRAPHY    ','MODEL_PARAM_ID  ',&
         'AXIS            ','TPFUN_SYMBOLS   ','QUIT            ',&
         'PARAMETER       ','EQUILIBRIA      ','RESULTS         ',&
         'CONDITIONS      ','SYMBOLS         ','LINE_EQUILIBRIA ',&
         'OPTIMIZATION    ','MODEL_PARAM_VAL ','ERROR_MESSAGE   ',&
         'ACTIVE_EQUILIBR ','ELEMENTS        ','                ']
!-------------------
! subsubcommands to LIST DATA
    character (len=16), dimension(nlform) :: llform=&
        ['SCREEN          ','TDB             ','MACRO           ',&
         'LATEX           ','PDB             ','                ']
!-------------------
! subsubcommands to LIST PHASE
    character (len=16), dimension(nclph) :: clph=&
        ['DATA            ','CONSTITUTION    ','MODEL           ',&
         '                ','                ','                ']
!-------------------
! subsubcommands to LIST OPTIMIZE results
    character (len=16), dimension(noptopt) :: optopt=&
        ['SHORT           ','LONG            ','COEFFICIENTS    ',&
         'GRAPHICS        ','DEBUG           ','MACRO           ',&
         'EXPERIMENTS     ','CORRELATION_MTRX','TC_RSD          ']
!------------------- subcommands to CALCULATE
    character (len=16), dimension(ncalc) :: ccalc=&
         ['TPFUN_SYMBOLS   ','PHASE           ','NO_GLOBAL       ',&
         'TRANSITION      ','QUIT            ','GLOBAL_GRIDMIN  ',&
         'SYMBOL          ','EQUILIBRIUM     ','ALL_EQUILIBRIA  ',&
         'WITH_CHECK_AFTER','                ','                ']
!-------------------
! subcommands to CALCULATE PHASE
    character (len=16), dimension(nccph) :: ccph=&
         ['ONLY_G          ','G_AND_DGDY      ','ALL_DERIVATIVES ',&
          'CONSTITUTION_ADJ','DIFFUSION_COEFF ','QUIT            ']
!-------------------
! subcommands to ENTER
    character (len=16), dimension(ncent) :: center=&
         ['TPFUN_SYMBOL    ','ELEMENT         ','SPECIES         ',&
         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
         'CONSTITUTION    ','EXPERIMENT      ','QUIT            ',&
         'EQUILIBRIUM     ','SYMBOL          ','OPTIMIZE_COEFF  ',&
         'COPY_OF_EQUILIB ','COMMENT         ','MANY_EQUILIBRIA ',&
         'MATERIAL        ','PLOT_DATA       ','GNUPLOT_TERMINAL',&
         '                ','                ','                ']
!-------------------
! subcommands to READ
    character (len=16), dimension(ncread) :: cread=&
        ['UNFORMATTED     ','TDB             ','QUIT            ',&
         'DIRECT          ','PDB             ','                ']
!-------------------
! subcommands to SAVE
! note SAVE TDB, MACRO, LATEX part of LIST DATA !!
    character (len=16), dimension(ncsave) :: csave=&
        ['TDB             ','SOLGAS          ','QUIT            ',&
         'DIRECT          ','UNFORMATTED     ','PDB             ']
!-------------------
! subcommands to AMEND first level
! many of these should be subcommands to PHASE
    character (len=16), dimension(ncam1) :: cam1=&
         ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
         'TPFUN_SYMBOL    ','CONSTITUTION    ','QUIT            ',&
         'COMPONENTS      ','GENERAL         ','                ',&
         'OPTIMIZING_COEFS','EQUILIBRIUM     ','                ',&
         'LINE            ','                ','                ']
!-------------------
! subsubcommands to AMEND PHASE
    character (len=16), dimension(ncamph) :: camph=&
         ['ADDITION        ','COMPOSITION_SET ','DISORDERED_FRACS',&
         'UNIQUAC_MODEL   ','DIFFUSION       ','DEFAULT_CONSTIT ',&
         '                ','FCC_PERMUTATIONS','BCC_PERMUTATIONS',&
         '                ','GADDITION       ','AQUEUS_MODEL    ',&
         'QUASICHEM_MODEL ','FCC_CVM_TETRAHDR','FLORY_HUGG_MODEL',&
         '                ','                ','QUIT            ']
! original
!         ['MAGNETIC_CONTRIB','COMPOSITION_SET ','DISORDERED_FRACS',&
!         'TWOSTATE_LIQUID ','SCHOTTKY_ANOMATY','DEFAULT_CONSTIT ',&
!         'LOWT_CP_MODEL   ','FCC_PERMUTATIONS','BCC_PERMUTATIONS',&
!         'ELASTIC_MODEL_1 ','GADDITION       ','AQUEUS_MODEL    ',&
!         'QUASICHEM_MODEL ','FCC_CVM_TETRAHDR','FLORY_HUGG_MODEL',&
!         'CRYSTAL_BREAKDWN','SECOND_EINSTEIN ','QUIT            ']
!-------------------
! subsubsubcommands to PHASE ADDITION
    character (len=16), dimension(naddph) :: caddph=&
         ['MAGNETIC_CONTRIB','QUIT            ','                ',&
         'TWOSTATE_LIQUID ','SCHOTTKY_ANOMATY','                ',&
         'LOWT_CP_MODEL   ','                ','                ',&
         'ELASTIC_MODEL_1 ','                ','SMOOTH_CP_STEP  ']
!-------------------
! subcommands to SET
    character (len=16), dimension(ncset) :: cset=&
         ['CONDITION       ','STATUS          ','ADVANCED        ',&
         '                ','INTERACTIVE     ','REFERENCE_STATE ',&
         'QUIT            ','ECHO            ','PHASE           ',&
         'UNITS           ','LOG_FILE        ','WEIGHT          ',&
         'NUMERIC_OPTIONS ','AXIS            ','INPUT_AMOUNTS   ',&
         'VERBOSE         ','AS_START_EQUILIB','BIT             ',&
         'VARIABLE_COEFF  ','SCALED_COEFF    ','OPTIMIZING_COND ',&
         'RANGE_EXPER_EQU ','FIXED_COEFF     ','SYSTEM_VARIABLE ',&
         'INITIAL_T_AND_P ','                ','                ']
! subsubcommands to SET STATUS
    character (len=16), dimension(ncstat) :: cstatus=&
         ['ELEMENT         ','SPECIES         ','PHASE           ',&
         'CONSTITUENT     ','                ','                ']
!        123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET ADVANCED
    character (len=16), dimension(ncadv) :: cadv=&
         ['EQUILIB_TRANSFER','QUIT            ','                ',&
          'GRID_DENSITY    ','SMALL_GRID_ONOFF','MAP_SPECIAL     ',&
          'GLOBAL_MIN_ONOFF','OPEN_POPUP_OFF  ','WORKING_DIRECTRY',&
          'HELP_POPUP_OFF  ','EET_EXTRAPOL    ','LEVEL           ',&
          '                ','                ','                ']
!         123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET BITS
    character (len=16), dimension(nsetbit) :: csetbit=&
         ['EQUILIBRIUM     ','GLOBAL          ','PHASE           ',&
          '                ','                ','                ']
!          123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET PHASE
    character (len=16), dimension(nsetph) :: csetph=&
         ['QUIT            ','STATUS          ','DEFAULT_CONSTIT ',&
          'AMOUNT          ','BITS            ','CONSTITUTION    ']
!         123456789.123456---123456789.123456---123456789.123456
!-------------------
! subsubsubcommands to SET PHASE BITS
! Most of these now selected by AMEND PHASE ... 
!    character (len=16), dimension(nsetphbits) :: csetphbits=&
!         ['FCC_PERMUTATIONS','BCC_PERMUTATIONS','IONIC_LIQUID_MDL',&
!         'AQUEOUS_MODEL   ','QUASICHEMICAL   ','FCC_CVM_TETRADRN',&
!         'FACT_QUASICHEMCL','NO_AUTO_COMP_SET','QUIT            ',&
!         'EXTRA_DENSE_GRID','FLORY_HUGGINS   ','                ',&
!         '                ','                ','                ']
! The bits can still be set here by numbers but the text is no longer shwon
    character (len=16), dimension(nsetphbits) :: csetphbits=&
         ['                ','                ','                ',&
         '                ','                ','                ',&
         '                ','NO_AUTO_COMP_SET','QUIT            ',&
         'EXTRA_DENSE_GRID','                ','                ',&
         '                ','                ','                ']

!         123456789.123456---123456789.123456---123456789.123456
!-------------------
! subcommands to STEP
    character (len=16), dimension(nstepop) :: cstepop=&
         ['NORMAL          ','SEPARATE        ','QUIT            ',&
          'CONDITIONAL     ','                ','                ']
!         123456789.123456---123456789.123456---123456789.123456
!-------------------
! subcommands to DEBUG
    character (len=16), dimension(ncdebug) :: cdebug=&
         ['FREE_LISTS      ','STOP_ON_ERROR   ','ELASTICITY      ',&
          'TEST1           ','TPFUN           ','BROWSER         ',&
          'TRACE           ','                ','                ']
!-------------------
! subcommands to SELECT, maybe some should be CUSTOMMIZE ??
    character (len=16), dimension(nselect) :: cselect=&
         ['EQUILIBRIUM     ','MINIMIZER       ','GRAPHICS        ',&
         'LANGUAGE        ','OPTIMIZER       ','                ']
!-------------------
! subcommands to DELETE
    character (len=16), dimension(nrej) :: crej=&
         ['ELEMENTS        ','SPECIES         ','PHASE           ',&
          'QUIT            ','COMPOSITION_SET ','EQUILIBRIUM     ',&
          'STEP_MAP_RESULTS','                ','                ']
!-------------------
! subcommands to INFORMATION
    character (len=16), dimension(ninf) :: cinf=&
         ['ELEMENTS        ','SPECIES         ','PHASE           ',&
          'QUIT            ','COMPOSITION_SET ','EQUILIBRIUM     ',&
          'CHANGES         ','                ','                ']
!-------------------
! subcommands to PLOT OPTIONS/ GRAPHICS OPTIONS
    character (len=16), dimension(nplt) :: cplot=&
        ['RENDER          ','SCALE_RANGES    ','RATIOS_XY       ',&
         'AXIS_LABELS     ','MANIPULATE_LINES','TITLE           ',&
         'GRAPHICS_FORMAT ','OUTPUT_FILE     ','GIBBS_TRIANGLE  ',&
         'QUIT            ','POSITION_OF_KEYS','APPEND          ',&
         'TEXT            ','TIE_LINES       ','FONT_AND_COLOR  ',&
         'LOGSCALE        ','LINE_WITH_POINTS','PAUSE_OPTIONS   ',&
         'MISCELLANEOUS   ','                ','                ']
!-------------------
!        123456789.123456---123456789.123456---123456789.123456
! minimizers
    character (len=16), dimension(2) :: minimizers=&
         ['LUKAS_HILLERT   ','SUNDMAN_HILLERT ']
!------------------------------------------------------------------------
! optimizers
    character (len=16), dimension(2) :: optimizers=&
         ['LMDIF           ','VA05AD          ']
!------------------------------------------------------------------------
!
! before we come here gtp_init has been called in the main program
! some defaults
    language=1
    logfil=0
    defcp=1
    seqxyz=0
! save the working directory (where OC is started)
    call getcwd(workingdir)
!    write(*,*)'Working directory is: ',trim(workingdir)
! this is used to save the path to any directory where a macro is started
!    macropath=' '
! initiate command line history
    myhistory%hpos=0
! defaults for optimizer
    nvcoeff=0
! iexit(2)=1 means listing scaled coefficients (Va05AD)
!    iexit=0
!    iexit(2)=1
! present the software
    write(kou,10)version,trim(linkdate),ocmonversion,gtpversion,hmsversion,&
         smpversion
10  format(/'Open Calphad (OC) software version ',a,', linked ',a,/&
         'with command line monitor version ',i2//&
         'This program is available with a GNU General Public License.'/&
         'It includes the General Thermodynamic Package, version ',A,','/&
         "Hillert's equilibrium calculation algorithm version ",A,','/&
         'step/map/plot software version ',A,' using GNUPLOT 5.0 graphics.'/&
         'Numerical routines are extracted from LAPACK and BLAS and'/&
         'the assessment procedure uses LMDIF from ANL.'/)
!
! lines starting with !$ will be included when compiling with -fopenmp
!$    write(kou,11)
11  format('Linked with OpenMp for parallel execution')
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Default gnuterminals, edit these as they may not be same on your systems
! Screen
    graphopt%gnutermid=' '
    i1=1
    graphopt%gnutermid(i1)='SCREEN '
! MAX 80 characters to set terminal .... HERE FONT AND SIZE IS SET
#ifdef aqplt
! Aqua plot screen on some Mac systems
!    graphopt%gnuterminal(i1)='aqua size 900,600 font "arial,20"'
    graphopt%gnuterminal(i1)='aqua size 940,700 font "arial,16"'
! it should be #elif not #elseif .... suck
#elif qtplt
! Qt plot screen on some LINUX systems
    graphopt%gnuterminal(i1)='qt size 940,700 font "arial,16"'
!    graphopt%gnuterminal(i1)='qt size 900,600 font "arial,16"'
#else
! wxt default plot screen (used on most Window systems)
    graphopt%gnuterminal(i1)='wxt size 940,700 font "arial,16"'
!    graphopt%gnuterminal(i1)='wxt size 900,600 font "arial,16"'
#endif
    graphopt%filext(i1)='  '
! NOTE THAT IN SCREEN PLOT WINDOW ONE CAN SELECT FILE OUTPUT
! Postscript
    i1=2
    graphopt%gnutermid(i1)='PS  '
    graphopt%gnuterminal(i1)='postscript color solid fontscale 1.2'
    graphopt%filext(i1)='ps  '
! Adobe Portable Document Format (PDF)
    i1=3
    graphopt%gnutermid(i1)='PDF '
#ifdef qtplt
! On LINUX
    graphopt%gnuterminal(i1)='pdfcairo '
#else
! NOTE size is in inch
!   graphopt%gnuterminal(i1)='pdf color solid size 6,4 enhanced fontscale 0.7'
!   graphopt%gnuterminal(i1)='pdf color solid size 6,4 enhanced font "arial,20"'
    graphopt%gnuterminal(i1)='pdf color solid size 6,5 enhanced font "arial,16"'
#endif
    graphopt%filext(i1)='pdf  '
! Graphics Interchange Format (GIF)
    i1=4
    graphopt%gnutermid(i1)='GIF  '
    graphopt%gnuterminal(i1)='gif enhanced fontscale 0.7'
!    graphopt%gnuterminal(i1)='gif enhanced '
    graphopt%filext(i1)='gif  '
    graphopt%gnutermax=i1
! Portable graphics format (PNG)
    i1=5
    graphopt%gnutermid(i1)='PNG  '
    graphopt%gnuterminal(i1)='png enhanced fontscale 0.7'
    graphopt%filext(i1)='png  '
    graphopt%gnutermax=i1
! by default spawn plots
    graphopt%status=ibset(graphopt%status,GRKEEP)
! if winhlp set also GRWIN to allow GRKEEP
! Possible problem here ... GRWIN=0 if Windows, GRWIN=1 if not Windows ...
!    write(*,*)'Setting windows bit 1: ',GRWIN
#ifdef winhlp    
!    write(*,*)'UI: Setting windows bit 2: ',GRWIN
    graphopt%status=ibset(graphopt%status,GRWIN)
#endif
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! jump here after NEW to reinitiallize all local variables also
20  continue
! clear file names
    ocmfile=' '; ocufile=' '; tdbfile=' '
! initiallize ploted, it is not done in reset_plotoptions
    graphopt%plotend='pause mouse'
! reset plot ranges and their defaults
    call reset_plotoptions(graphopt,plotfile,textlabel)
    axplotdef=' '
! default list unit
    optionsset%lut=kou
    lut=kou
! default for list short
    chshort='A'
! set default minimizer, 2 is matsmin, 1 does not work ...
    minimizer=2
! set default optimimzer, 1 is LMDIF, 2 is VA05AD (no longer available)
    optimizer=1
! by default no stop on error and no logfile
    stop_on_error=.false.
    logfil=0
    buperr=0
!
! nopenpopup is declared in metlib3.F90 and dis/allow open popup windows
! it is initiated to FALSE, if user change it will not be reinitiated here
!    nopenpopup=.FALSE.
! in init_gtp the first equilibrium record is created and 
! firsteq has been set to that
!
!25  continue
! default values of T and P.  NOTE these are not set as conditions
    firsteq%tpval(1)=1.0D3
    firsteq%tpval(2)=1.0D5
!
! default list result option
    lrodef=1
! default axis limits set to be 0 and 1
    maxax=5
    noofaxis=0
! state variable for plot axis (only 2)
    do j1=1,2
       axplotdef(j1)=' '
    enddo
! remove any results from step and map
    nullify(maptop)
    nullify(mapnode)
    nullify(maptopsave)
! entered start equilibria
    nullify(starteq)
    noofstarteq=0
! set default fractions when entering composition
    xknown=one
! set default equilibrium to 1 and current equilibrium (CEQ) to firsteq
    neqdef=1
    ceq=>firsteq
! >>> we should remove all equilibria !! ??
! here one should read a user initialisation file as a macro
! file can be at current directory or at home directory
! initiate on-line help
! local environment: please create OCHOME as an environment variable
    ochome=' '
    call get_environment_variable('OCHOME ',ochome)
    startupmacro=.FALSE.
! if help is not set then set these filenames as blanks
    browser=' '
    latexfile=' '
    htmlfile=' '
    noochome: if(ochome(1:1).eq.' ') then
       inquire(file='ochelp.tex ',exist=logok)
       if(.not.logok) then
          write(*,*)'Warning, no help file'
       else
! help file on local directory ... no browser nor HMTL
          call init_help(' ','ochelp.tex ',' ')
       endif
    else ! there is a OCHOME environment variable
! both LINUX and WINDOWS accept / as separator between directory and file names
       write(*,*)'OC home directory (OCHOME): ',trim(ochome)
! when we are here OCHOME is not empty!
#ifdef winhlp
! THIS IS FOR WINDOWS
! first argument is browser
! second is LaTeX source file
! third argument is HTML file with hypertargets <a id="target" />
       browser='C:\PROGRA~1\INTERN~1\iexplore.exe '
! normal tex/html help files
       latexfile=trim(OCHOME)//'\'//'ochelp.tex'
       htmlfile=trim(OCHOME)//'\'//'ochelp.html'
! special for on-line HELP TESTING
!       write(*,*)' *** testing on-line help'
!       latexfile='c:\users\bosse\documents\oc\oc\src\manual\ochelp5.tex'
!       htmlfile='c:\users\bosse\documents\oc\oc\src\manual\ochelp5.html'
#elif lixhlp
! THIS IS FOR LINUX
!       browser='/usr/bin/firefox '
       browser='firefox '
       latexfile=trim(OCHOME)//'/'//'ochelp.tex'
       htmlfile=trim(OCHOME)//'/'//'ochelp.html'
#endif
       inquire(file=browser,exist=logok)
       if(.not.logok) then
          write(*,*)'No browser for on-line help'
          htmlhelp=.FALSE.
          call init_help(' ',latexfile,' ')
       else
!          write(*,*)'Using: ',trim(browser),' for on-line help'
! the *.tex is for search, *.html for browser
          inquire(file=htmlfile,exist=logok)
          if(.not.logok) then
! if we have no html file disable the help_popup
!             write(*,*)'No online popup helpfile: ',trim(htmlfile)
             htmlhelp=.FALSE.
             call init_help(' ',latexfile,' ')
          else
! htmlhelp is set TRUE inside init_help if there is a HMTL file in the call
!          htmlhelp=.TRUE.
             call init_help(browser,latexfile,htmlfile)
             write(*,*)'On-line help using ',trim(browser)
!          write(*,*)'Help LaTeX: ',trim(latexfile)
!          write(*,*)'Help popup: ',trim(htmlfile)
          endif
       endif
!       if(.not.htmlhelp) then
! no html file or user have explicitly set HELP_POPUP_OFF, respect that
!          call init_help(' ',latexfile,' ')
!       else
!          call init_help(browser,latexfile,htmlfile)
!       endif
! default directory for databases
       ocbase=trim(ochome)//'/databases'
       inquire(file=trim(ocbase),exist=logok)
!       if(logok) then
!          write(*,*)'There is a database directory: ',trim(ocbase)
!       else
!          write(*,*)'No database directory'
!       endif
! running a initial macro file
       cline=trim(ochome)//'/start.OCM '
       inquire(file=cline,exist=logok)
       if(logok) then
          write(*,*)'Reading your startup macro: ',trim(cline)
          last=0
! This just open the file and sets input unit to file
          call macbeg(cline,last,logok)
          startupmacro=.TRUE.
!          macropath=string
!       else
!          write(*,*)'No initiation file'
       endif
    endif noochome
    write(*,*)'Working directory is: ',trim(workingdir)
!
! finished initiallization
!
!
!============================================================
! return here for next command
100 continue
    if(gx%bmperr.ne.0) goto 990
    if(buperr.ne.0) goto 990
! turn off any options set
    call ocmon_reset_options(optionsset)
! initiate command level for help routines
    call helplevel1('Initiate help level for OC')
! handling of inline arguments ONCE
    if(.not.startupmacro .and. narg.gt.0) then
! at present accept only one argument assumed to be a macro file name
       narg=0
       cline=argline(1)
       last=0
!       if(cline(1:1).eq.'<') then
!          write(*,*)'OC reads can start a macro from the command line'
!       else
          call macbeg(cline,last,logok)
!          macropath=string
!       endif
    endif
! read the command line with gparc to have output on logfile
! NOTE read from macro file if set.
    last=len(aline)
    aline=' '
    cline=' '
!    call gparc(ocprompt,aline,last,5,cline,' ',tophlp)
    call gparc(ocprompt,aline,last,5,cline,' ',q2help)
    j1=1
    if(len_trim(cline).gt.80) then
       write(kou,101)
101    format(' *** Warning: long input lines may be truncated',&
            ' and cause errors')
    endif
! with empty line just prompt again
    if(eolch(cline,j1)) goto 100
! with macro command character just prompt again
    if(cline(j1:j1).eq.'@') goto 100
    kom=ncomp(cline,cbas,ncbas,last)
!  write(kou,*)'kom= ',kom,last
    if(kom.le.0) then
       if(kom.lt.0) then
          write(kou,*)'Ambiguous command, available commands are:'
       else
          write(kou,*)'No such command, available commands are:'
       endif
       last=1
       cline=' *'
       call q3help(cline,last,cbas,ncbas)
       write(*,*)'An OS command must be prefixed by @'
       goto 100
    else
! check for options .... some of these do not work yet
! one should check for options after each subcommand or value entered ??
!       call ocmon_set_options(cline,last,optionsset)
       nops=0
110    continue
       if(.not.eolch(cline,last)) then
          if(cline(last:last).eq.'/') then
! this is an option!
             call getext(cline,last,2,option,' ',nopl)
             if(buperr.ne.0) then
                write(kou,*)'Error reading option',buperr
                buperr=0; goto 100
             endif
             call ocmon_set_options(option,afo,optionsset)
             if(afo.ne.0) then
                write(kou,*)'Please give the command again'
                goto 100
             endif
             goto 110
          else
! set "last" back one character to prepare for next call of GPARx 
! as the first thing done by GPARx is to increment last by 1 to bypass a ,
             last=last-1   
          endif
       endif
    endif
! save command for help path
    if(helprec%level.lt.maxhelplevel) then
       helprec%level=helprec%level+1
       helprec%cpath(helprec%level)=cbas(kom)
    else
       write(*,*)'Warning, exceeded helprec%level limit 1'
    endif
! The IF loop is for handling of defaults in submenu. "l ,,,,," took all 
! defaults but "l,,,,," did not ....
! if last>1 and cline(last-1:last-1) is a space and cline(last:last) a comma,
! increment last
    if(last.eq.1) then
       last=last+1
    elseif(last.lt.len(cline)) then
       if(cline(last:last).ne.' ') then
          if(cline(last+1:last+1).eq.',') last=last+1
       endif
    endif
!
!================================================ separating main commands
!------------------------- separating subcommands
!......................... separating subsubcommands
! jump here if there is an inline argument
! 99  continue
    main: SELECT CASE(kom)
! command selection
!=================================================================
    CASE DEFAULT
       write(kou,*)'No such command'
       goto 100
!=================================================================
    CASE(1) ! AMEND
! amend subcommands
!       ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
!        'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
!        'TPFUN_SYMBOL    ','CONSTITUTION    ','QUIT            ',&
!        'COMPONENTS      ','GENERAL         ','                ',&
!        'OPTIMIZING_COEFS','EQUILIBRIUM     ','                ',&
!        'LINE            ','                ','                ']
! disable continue optimization
!       iexit=0
!       iexit(2)=1
       kom2=submenu(cbas(kom),cline,last,cam1,ncam1,4)
       amend: SELECT CASE(kom2)
       CASE DEFAULT
          write(kou,*)'No such answer'
          goto 100
!-------------------------
       case(1) ! amend symbol (of state variable function)
          call gparc('Symbol name: ',cline,last,1,name1,' ',q1help)
          call capson(name1)
          do svss=1,nosvf()
             if(name1(1:16).eq.svflista(svss)%name) goto 1020
          enddo
          write(kou,*)'No such symbol'
          goto 100
1020      continue
          if(btest(svflista(svss)%status,SVCONST)) then
! if symbol is just a numeric constant we can change its value
             actual_arg=' '
             xxx=evaluate_svfun_old(svss,actual_arg,1,ceq)
             call gparrd('Give new value: ',cline,last,xxy,xxx,q1help)
             call set_putfun_constant(svss,xxy)
             goto 100
          endif
! If the user answers anything but y or Y the bits are cleared
! symbol can be evaluated only explicitly (like variables in TC)
          call gparcd('Should the symbol only be evaluated explicitly?',&
               cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
             svflista(svss)%status=ibset(svflista(svss)%status,SVFVAL)
          else
             svflista(svss)%status=ibclr(svflista(svss)%status,SVFVAL)
          endif
! symbols can be evaluated in just a particular equilibrium
! like H298 for experimental data on H(T)-H298
! BEWARE: if equilibria are calculated in threads this must be calculated
! before the parallelization, testing bit EQNOTHREAD
          call gparcd(&
               'Should the symbol be evaluated for a specific equilibrium?',&
               cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
             ll=ceq%eqno
             call gparid('Specify equilibrium number:',cline,last,&
                  neqdef,ll,q1help)
             svflista(svss)%status=ibset(svflista(svss)%status,SVFEXT)
             svflista(svss)%eqnoval=neqdef
! set status bit that this equilibrium must be calculated before parallel calc
             ceq%status=ibset(ceq%status,EQNOTHREAD)
          else
             svflista(svss)%status=ibclr(svflista(svss)%status,SVFEXT)
             svflista(svss)%eqnoval=0
! this cannot be cleared, there may be other threadprotected symbols
!           ceq%status=ibclr(ceq%status,EQNOTHREAD)
          endif
!-------------------------
       case(2) ! amend element
          write(kou,*)'Not implemented yet'
!-------------------------
       case(3) ! amend species
          call gparc('Species symbol: ',cline,last,1,name1,' ',q1help)
          call find_species_record(name1,loksp)
          if(gx%bmperr.ne.0) goto 100
          write(*,'(a)')'You can only amend UNIQAC area and segments'
          call gparrd('UNIQAC surface area (q): ',cline,last,xxx,one,q1help)
          if(xxx.le.zero) then
             write(*,'(a)')'Area must be >0, set to default 1.00'
             xxx=one
          endif
          call gparrd('UNIQAC segments (r): ',cline,last,xxy,one,q1help)
          if(xxy.le.zero) then
             write(*,'(a)')'Segments must be >0, set to default 1.00'
             xxy=one
          endif
! mark UNIQUAC in species status word and allocate space for values
          call set_uniquac_species(loksp)
          if(gx%bmperr.ne.0) goto 100
          call enter_species_property(loksp,1,xxx)
          call enter_species_property(loksp,2,xxy)
          if(gx%bmperr.ne.0) goto 100
!-------------------------
       case(4) ! amend phase subcommands
          call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
          if(buperr.ne.0) goto 990
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
! sometimes lokph is used below
          call get_phase_record(iph,lokph)
          call get_phasetup_name(iph,name1)
!
          kom3=submenu('AMEND for phase '//trim(name1),&
               cline,last,camph,ncamph,2)
!          write(*,*)'Amend phase subcommand: ',kom3
          amendphase: SELECT CASE(kom3)
! subsubcommands to AMEND PHASE
!         ['ADDITION        ','COMPOSITION_SET ','DISORDERED_FRACS',&
!         '                ','                ','DEFAULT_CONSTIT ',&
!         '                ','FCC_PERMUTATIONS','BCC_PERMUTATIONS',&
!         '                ','GADDITION       ','AQUEUS_MODEL    ',&
!         'QUASICHEM_MODEL ','FCC_CVM_TETRAHDR','FLORY_HUGG_MODEL',&
!         '                ','                ','QUIT            ']
!....................................................
          CASE DEFAULT
             write(kou,*)'Amend phase subcommand error'
!....................................................
          case(1) ! amend phase addition
             kom4=submenu('Addition of',cline,last,caddph,naddph,1)
!          write(*,*)'Amend phase addition: ',kom4
!         ['MAGNETIC_CONTRIB','QUIT            ','                ',&
!         'TWOSTATE_LIQUID ','SCHOTTKY_ANOMATY','                ',&
!         'LOWT_CP_MODEL   ','                ','                ',&
!         'ELASTIC_MODEL_A ','QUASICHEM_MODEL ','FCC_CVM_TETRAHDR']
!
             amendphaseadd: SELECT CASE(kom4)
             case(1) ! amend phase <name> magnetic contribution
                idef=-3
! zero value of antiferromagnetic factor means Inden-Qing model
                call gparid('Antiferromagnetic factor: ',&
                     cline,last,j1,idef,q1help)
                if(buperr.ne.0) goto 990
                if(j1.eq.0) then
! Xiong modification of Inden-Hillert-Jarl magnetic model has AFF=0
                   call gparcd('BCC type phase: ',cline,last,1,ch1,'N',q1help)
                   call gparcd('Using individual Bohr magnetons: ',&
                        cline,last,1,ch1,'N',q1help)
                   if(ch1.ne.'N') then
                      call set_phase_status_bit(lokph,PHBMAV)
                   endif
                   j2=xiongmagnetic
                   call add_addrecord(lokph,ch1,xiongmagnetic)
                else
                   if(j1.eq.-1) then
! Inden magnetic for BCC
                      call add_addrecord(lokph,'Y',indenmagnetic)
                   else
! Inden magnetic for FCC
                      call add_addrecord(lokph,'N',indenmagnetic)
                   endif
                   j2=indenmagnetic
                endif
                call gparcd('Is the addition calculated for one mole? ',&
                     cline,last,1,ch1,'N',q1help)
! The magnetic model calculates a molar Gibbs energy, must be multiplied with
! the number of atoms in the phase. j2 set above to the addition type
                if(ch1.eq.'Y' .or. ch1.eq.'y') then
                   call setpermolebit(lokph,j2)
                endif
!....................................................
             case(2) ! QUIT
                goto 100
!....................................................
             case(3) ! not used
                continue
!....................................................
             case(4) ! amend phase <name> addition twostate_liquid model
                call add_addrecord(lokph,' ',twostatemodel1)
                write(kou,667)
667             format('This addition require THET parameters for the',&
                     ' Einstein T of the amorphous state'/'and G2 parameters',&
                     ' for the transition to the liquid state.')
!....................................................
             case(5) ! amend phase <name> addition Schottky anomality
                call add_addrecord(lokph,' ',schottkyanomality)
                write(*,668)
668             format('This addition requires the TSCH and CSCH parameters')
!....................................................
             case(6) ! not used
!....................................................
             case(7) ! amend phase <name> LowT_CP_model
                call add_addrecord(lokph,' ',einsteincp)
                write(*,*)'This addition requires the THET parameter'
                call gparcd('Is the addition calculated for one mole? ',&
                     cline,last,1,ch1,'Y',q1help)
! The CP model calculates a molar Gibbs energy, must be multiplied with
! the number of atoms in the phase. j2 set above to the addition type
                if(ch1.eq.'Y' .or. ch1.eq.'y') then
                   call setpermolebit(lokph,einsteincp)
                endif
!....................................................
             case(8) ! not used
!....................................................
             case(9) ! not used
!....................................................
             case(10) ! amend phase elastic model
!                call add_addrecord(lokph,' ',elasticmodel1)
                write(*,*)'This addition is not yet implemented'
                !....................................................
             case(11) ! amend phase ... unused
                continue
!....................................................
             case(12) ! amend phase ... smooth-Cp-step
                call add_addrecord(lokph,' ',secondeinstein)
                write(*,672)
672             format('This addition recures the THT2 and DCP2 parameters')
! The smooth CP model calculates a molar Gibbs energy, must be multiplied with
! the number of atoms in the phase. j2 set above to the addition type
                if(ch1.eq.'Y' .or. ch1.eq.'y') then
                   call setpermolebit(lokph,secondeinstein)
                endif
             end select amendphaseadd
!************************************ end of amend phase ... addition
!....................................................
          case(2) ! amend phase <name> composition set add/remove
             call gparcd('Add new set? ',cline,last,1,ch1,'Y ',q1help)
             if(buperr.ne.0) goto 990
             if(ch1.eq.'Y' .or. ch1.eq.'y') then
                call gparc('Prefix: ',cline,last,1,prefix,' ',q1help)
                call gparc('Suffix: ',cline,last,1,suffix,' ',q1help)
                call enter_composition_set(iph,prefix,suffix,ics)
                if(gx%bmperr.ne.0) goto 990
! list the number of new composition set
                write(kou,*)'New composition set number is ',ics
! ask for default constitution of new set
                call ask_default_constitution(cline,last,iph,ics,ceq)
             else
! remove the highest (max is 9).  Can be dangerous.  Can not be made if there
! are several equilibra unless second argument is changed to .TRUE.
                call remove_composition_set(iph,.FALSE.)
                if(gx%bmperr.ne.0) goto 990
             endif
!....................................................
          case(3) ! amend phase <name> disordered_fracset
             if(.not.allowenter(2)) then
                gx%bmperr=4125
                goto 990
             endif
! we should check the number of sublattices of the phase ...
             idef=4
             call gparid('Sum up to sublattice: ',cline,last,ndl,idef,q1help)
             if(buperr.ne.0) goto 990
! ch1 is parameter suffix, j1=0 means never completely disordred (sigma)
!             call gparcd('Can the phase be totally disordered? ',cline,last,&
!                  1,ch1,'N',q1help)
! Answer Y means the ordered parameters are calculated twice, once with
! the disordered composition.  Answer N for phase which are never disordered
!             call gparcd('Is the disordered part assessed indepently? ',&
!                  cline,last,1,ch1,'N',q1help)
             call gparcd('Should ordered part cancel when disordered? ',&
                  cline,last,1,ch1,'N',q1help)
             if(buperr.ne.0) goto 990
             if(ch1.eq.'N') then
! like sigma which is never completely disordered
                j1=0
             else
! like FCC ordering where the disordered state can be modeled independently
                j1=1
                write(kou,*)'This phase can be totally disordered'
             endif
             ch1='D'
             call add_fraction_set(iph,ch1,ndl,j1)
             if(gx%bmperr.ne.0) goto 990
!....................................................
          case(4) ! UNIQUAC model
             write(*,*)'Not implemented yet'
!....................................................
          case(5) ! DIFUSION properties
! copy the rest of the line to the subroutine
             text=cline(last:)
             call add_addrecord(lokph,text,DIFFCOEFS)
!....................................................
          case(6) ! amend phase <name> default_constitution
! to change default constitution of any composition set give #comp.set.
             call ask_default_constitution(cline,last,iph,ics,ceq)
!....................................................
          case(7) ! moved
!....................................................
          case(8) ! amend phase ... FCC_PERMUTATIONS
             if(check_minimal_ford(lokph)) goto 100
             call set_phase_status_bit(lokph,PHFORD)
!....................................................
          case(9) ! amend phase ... BCC_PERMUTATIONS
                if(check_minimal_ford(lokph)) goto 100
                call set_phase_status_bit(lokph,PHBORD)
!....................................................
          case(10) ! moved
!....................................................
          case(11) ! AMEND PHASE GADDITION
! add a constant term to G, value in J/FU
             lokcs=phasetuple(iph)%lokvares
             if(allocated(ceq%phase_varres(lokcs)%addg)) then
                xxy=ceq%phase_varres(lokcs)%addg(1)
             else
! maybe we will use more terms later ....
                allocate(ceq%phase_varres(lokcs)%addg(1))
             endif
             call gparrd('Addition to G in J/FU (formula units): ',&
                  cline,last,xxx,xxy,q1help)
             ceq%phase_varres(lokcs)%addg(1)=xxx
! set bit that this should be calculated
             ceq%phase_varres(lokcs)%status2=&
                  ibset(ceq%phase_varres(lokcs)%status2,CSADDG)
!....................................................
          case(12) ! amend phase ... aqueous model
             write(*,*)'Not implemented yet'
!             call set_phase_status_bit(lokph,PHAQ1)
!....................................................
          case(13) ! amend phase ... quasicemichal model (several)
             call gparid('Quasichemical type: ',cline,last,jp,3,q1help)
             if(jp.lt.0 .or. jp.gt.3) then
                write(*,*)'Value must be between 1 and 3'
             else
                qcmodel=jp
             endif
!             write(kou,*)'Not implemented yet'
! Future model bits
!                call set_phase_status_bit(lokph,PHQCE)
!                call set_phase_status_bit(lokph,PHFACTCE)
!....................................................
          case(14) ! amend phase ... FCC_CVM_TETRAHEDRON MODEL
             write(kou,*)'Not implemented yet'
!                call set_phase_status_bit(lokph,PHCVMCE)
!....................................................
          case(15) ! amend phase ... FLORY_HUGGINS_MODEL
             write(kou,*)'Not implemented yet'
!             call set_phase_status_bit(lokph,PHFHV)
!             call clear_phase_status_bit(lokph,PHID)
!....................................................
          case(16) ! moved
!....................................................
          case(17) ! moved
!....................................................
          case(18) ! amend phase ... quit
             goto 100
          END SELECT amendphase
!------------------------- end of amend phase
       case(5) ! amend parameter
          write(kou,*)'Not implemented yet, only ENTER PARAMETER'
!-------------------------
       case(6) ! amend bibliography
          call enter_bibliography_interactivly(cline,last,1,j1)
!-------------------------
       case(7) ! amend TPFUN symbol
          write(kou,*)' *** Dangerous if you have several equilibria!'
          call gparc('TP-fun symbol: ',cline,last,1,name1,' ',q1help)
          call find_tpsymbol(name1,idef,xxx)
          if(gx%bmperr.ne.0) then
             write(*,*)'Ambiguouos or unknown symbol'; goto 990
          endif
          if(idef.eq.0) then
! it is a function , this call just read the function starting with low T etc.
             call enter_tpfun_interactivly(cline,last,funstring,jp)
! this stores the tpfun, lrot<0 means the symbol already exists
             lrot=-1
             call enter_tpfun(name1,funstring,lrot,.FALSE.)
             if(gx%bmperr.ne.0) goto 990
! mark functions not calculated.  This should be done in all ceq ...
             ceq%eq_tpres(lrot)%tpused(1)=-one
             ceq%eq_tpres(lrot)%tpused(2)=-one
          elseif(idef.eq.2) then
             write(*,*)'You cannot change an optimizing coefficients'
             goto 100
          else
! it is a constant, you can change if
             call gparrd('Value: ',cline,last,xxy,xxx,q1help)
             call capson(name1)
             call enter_tpconstant(name1,xxy)
          endif
!-------------------------
       case(8) ! amend constitution (also as ENTER CONST and SET PHASE )
          call ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
          if(gx%bmperr.ne.0) goto 990
!-------------------------
       case(9) ! QUIT amend
          continue
!-------------------------
       case(10) ! amend components
          write(*,*)'WARNING: not fully implemented yet'
!          goto 100
          if(associated(ceq%lastcondition)) then
             write(kou,*)'Warning: All your conditions will be removed'
          endif
          i2=1
          line=' '
          do i1=1,noel()
             call get_component_name(i1,line(i2:),ceq)
             i2=len_trim(line)+2
          enddo
          aline=' '
          call gparcd('Give all new components: ',cline,last,&
               5,aline,line,q1help)
! option is a character with the new components ...
          call amend_components(aline,ceq)
          if(gx%bmperr.ne.0) goto 990
!-------------------------
       case(11) ! amend general
          call amend_global_data(cline,last)
!-------------------------
       case(12) ! amend unused
          write(*,*)'Not implemented yet'
!-------------------------
       case(13) ! amend OPTIMIZING_COEFF, (rescale or recover)
          call gparcd('Should the coefficients be rescaled?',&
               cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
! set start values to current values
!             firstash%coeffstart=firstash%coeffvalues*firstash%coeffscale
!             firstash%coeffscale=firstash%coeffstart
! Note the "current value" is "start value" times "scaling factor"
!             firstash%coeffvalues=one
             do j2=0,size(firstash%coeffstate)-1
                if(firstash%coeffstate(j2).ge.10) then
                   call get_value_of_constant_index(firstash%coeffindex(j2),xxx)
                   if(gx%bmperr.ne.0) then
                   write(*,*)'Error getting value of assessment coefficient',j2
                      goto 100
                   endif
!                   write(*,*)'Assessment coefficient value: ',xxx
! Set all values equal to the current value of the TP variable ...
                   firstash%coeffscale(j2)=xxx
                   firstash%coeffstart(j2)=xxx
                   firstash%coeffvalues(j2)=one
!                   call change_optcoeff(firstash%coeffindex(j2),xxx)
                endif
             enddo
             firstash%coeffrsd=zero
             call listoptcoeff(mexp,err0,.FALSE.,lut)
             if(allocated(cormat)) then
                deallocate(cormat)
                deallocate(tccovar)
             endif
          else
             call gparcd('Do you want to recover the coefficients values?',&
                  cline,last,1,ch1,'N',q1help)
             if(ch1.eq.'y' .or. ch1.eq.'Y') then
! set current optimizing values back to start values
!                firstash%coeffvalues=firstash%coeffstart*firstash%coeffscale
                do j2=0,size(firstash%coeffstate)-1
! This affects only current optimizing coefficients!!
                   if(firstash%coeffstate(j2).ge.10) then
                      xxx=firstash%coeffstart(j2)
                      firstash%coeffvalues(j2)=xxx/firstash%coeffscale(j2)
! we must also change the value of the associated TP variable ??
                      call change_optcoeff(firstash%coeffindex(j2),xxx)
                   endif
                enddo
! no change of start value or scaling factor but zero RSD and sum of squares
                firstash%coeffrsd=zero
                if(allocated(cormat)) then
                   deallocate(cormat)
                   deallocate(tccovar)
                endif
                err0(2)=zero
                call listoptcoeff(mexp,err0,.FALSE.,lut)
             else
                write(kou,557)
557             format('Nothing done as there are no other amend',&
                     ' optimizing option')
             endif
          endif
!-------------------------
       case(14) ! AMEND EQUILIBRIUM intended to add to experimental list
          write(*,*)'Not implemented yet'
!-------------------------
       case(15) ! Nothing defined
          write(*,*)'Not implemented yet'
!-------------------------
       case(16) ! AMEND LINEs of calculated equilibria
! possible amendment of all stored equilibria as ACTIVE or INACTIVE
          call amend_stored_equilibria(axarr,maptop)
!-------------------------
       case(17) ! Nothing defined
          write(*,*)'Not implemented yet'
!-------------------------
       case(18) ! Nothing defined
          write(*,*)'Not implemented yet'
       END SELECT amend
!=================================================================
! calculate subcommands
!         ['TPFUN_SYMBOLS   ','PHASE           ','NO_GLOBAL       ',&
!         'TRANSITION      ','QUIT            ','GLOBAL_GRIDMIN  ',&
!         'SYMBOL          ','EQUILIBRIUM     ','ALL_EQUILIBRIA  ',&
!         'WITH_CHECK_AFTER','                ','                ']
    CASE(2)
       kom2=submenu(cbas(kom),cline,last,ccalc,ncalc,8)
       calculate: SELECT CASE(kom2)
       CASE DEFAULT
          write(kou,*)'No such calculate command'
          goto 100
!-------------------------
       CASE(1) ! calculate TPFUN symbols , use current values of T and P
          call gparcd('name: ',cline,last,5,name1,'*',q1help)
          lrot=0
          iel=index(name1,'*')             
          if(iel.gt.1) name1(iel:)=' '
! as TP functions call each other force recalculation and calculate all
! even if just a single function is requested
          call change_optcoeff(-1,zero)
          do j1=1,notpf()
             call eval_tpfun(j1,ceq%tpval,val,ceq%eq_tpres)
             if(gx%bmperr.gt.0) goto 990
          enddo
          if(name1(1:1).ne.'*') then
             once=.TRUE.
2009         continue
             call find_tpfun_by_name(name1,lrot)
!             write(*,*)'cui: ',lrot,iel,gx%bmperr
             if(gx%bmperr.ne.0) then
                if(iel.eq.0) goto 990
                gx%bmperr=0
             else
! found function number from lrot ???
                j1=lrot
                call eval_tpfun(j1,ceq%tpval,val,ceq%eq_tpres)
                if(gx%bmperr.gt.0) goto 990
                if(once) then
                   once=.FALSE.
                   write(lut,2011)100000,ceq%tpval
                endif
                write(lut,2012)j1,val
!                call list_tpfun(lrot,0,longstring)
!                call wrice2(lut,0,12,78,1,longstring)
                if(iel.gt.1) goto 2009
             endif
          else
             write(lut,2011)notpf(),ceq%tpval
2011         format(/'Calculating ',i6,' functions for T,P=',F10.2,1PE15.7/&
                  3x,'No   F',11x,'F.T',9x,'F.P',9x,'F.T.T',&
                  7x,'F.T.P',7x,'F.P.P')
!             call cpu_time(starting)
             do j1=1,notpf()
                call eval_tpfun(j1,ceq%tpval,val,ceq%eq_tpres)
                if(gx%bmperr.gt.0) goto 990
                write(lut,2012)j1,val
2012            format(I5,1x,6(1PE12.4))
             enddo
!             call cpu_time(ending)
          endif
!          write(kou,2013)ending-starting
!2013      format('CPU time used: ',1pe15.6)
!---------------------------------------------------------------
       case(2) ! calculate phase, _all _only_g or _g_and_dgdy, etc
! asks for phase name and constitution.  DO NOT ALLOW * by setting iph=-1
! before calling!
          iph=-1
          call ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
          if(gx%bmperr.ne.0) goto 990
! if iph<0 then * has been given as phase name
          if(iph.lt.0) then
             write(kou,*)'Cannot loop for all phases'
             goto 100
          endif
! subcommands for calculate phase
!         ['ONLY_G          ','G_AND_DGDY      ','ALL_DERIVATIVES ',&
!          'CONSTITUTION_ADJ','DIFFUSION_COEFF ','QUIT            ']
!
          kom3=submenu('Calculate what for phase?',cline,last,ccph,nccph,defcp)
!        if(kom2.le.0) goto 100
!        ph-a ph-G ph-G+dg/dy
          defcp=kom3
          lut=optionsset%lut
! use current value of T and P
          if(kom3.ne.4) then
             write(*,2015)ceq%tpval
2015         format('Using T=',F9.2,' K and P=',1pe14.6,&
                  ' Pa, results in J/F.U.')
          endif
          rgast=globaldata%rgas*ceq%tpval(1)
! this is the number of moles formula units the user specified
          cpham=ceq%phase_varres(lokcs)%amfu
          calcphase: SELECT CASE(kom3)
!.......................................................
          CASE DEFAULT
             write(kou,*)'Calculate phase subcommand error'
!.......................................................
          case(1) ! calculate phase < > only G
             call calcg(iph,ics,0,lokres,ceq)
             if(gx%bmperr.ne.0) goto 990
             parres=>ceq%phase_varres(lokres)
             write(lut,2031)(cpham*rgast*parres%gval(j1,1),j1=1,4)
! G=H-T*S; H=G+T*S; S=-G.T; H = G + T*(-G.T) = G - T*G.T
             write(lut,2032)cpham*parres%gval(1,1)/parres%abnorm(1),&
                  cpham*(parres%gval(1,1)-ceq%tpval(1)*parres%gval(2,1))*rgast,&
                  parres%abnorm(1)
2031         format(/'G, dG/dT dG/dP d2G/dT2:',4(1PE14.6))
2032         format('G/RT, H, atoms/F.U:',3(1PE14.6))
! also list contributions from calculated additions ...!!!
             call list_addition_values(lut,parres)
!.......................................................
          case(2) ! calculate phase < >  G and dG/dy
             call calcg(iph,ics,1,lokres,ceq)
             if(gx%bmperr.ne.0) goto 990
             parres=>ceq%phase_varres(lokres)
             nofc=noconst(iph,ics,firsteq)
             write(lut,2031)(cpham*rgast*parres%gval(j1,1),j1=1,4)
             write(lut,2041)(rgast*parres%dgval(1,j1,1),j1=1,nofc)
2041         format('dG/dy:   ',4(1PE16.8),(/9x,4e16.8)/&
                  ' NOTE THAT dG/dy_i is NOT THE CHEMICAL POTENTIAL of i!')
!.......................................................
          case(3) ! calculate phase < > all derivatives
             call gparid('Number of times: ',cline,last,times,1,q1help)
             call tabder(iph,ics,times,ceq)
             if(gx%bmperr.ne.0) goto 990
             write(*,2042)
2042         format('Values are per mole formula unit'/&
                  ' NOTE THAT dG/dy_i is NOT THE CHEMICAL POTENTIAL of i!')
!             if(gx%bmperr.ne.0) goto 990
!.......................................................
          case(4,5) ! calculate phase with constitution_adjustment
! or derivatives of chemical potentials and mobility data
! convert to phase tuple here as that is used in the application call
             do jp=1,nooftup()
!                if(phasetuple(jp)%phaseix.eq.iph .and. &
                if(phasetuple(jp)%ixphase.eq.iph .and. &
                     phasetuple(jp)%compset.eq.ics) then
!                   write(*,*)'This is phase tuple ',jp
                   goto 2044
                endif
             enddo
             write(*,*)'No such tuple'
             goto 100
2044         continue
             phtup=>phasetuple(jp)
! Get current constitution of the phase
             call calc_phase_molmass(iph,ics,xknown,aphl,totam,xxy,xxx,ceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Error finding current composition'
                goto 990
             endif
! ask for overall composition
             totam=one
             quest='Mole fraction of XX:'
             do nv=1,noel()-1
                if(totam.gt.zero) then
! assume elements as components
                   call get_component_name(nv,elsym,ceq)
                   quest(18:19)=elsym
! prompt with current mole fraction:
                   call gparrd(quest,cline,last,xxy,xknown(nv),q1help)
                   if(buperr.ne.0) then
                      buperr=0; xxy=zero
                   endif
                   if(xxy.gt.totam) then
                      write(kou,*)'Fraction too large, set to ',totam
                      xxy=totam
                   endif
                else
                   xxy=zero
                endif
                xknown(nv)=xxy
                totam=totam-xxy
! yarr is used here to provide an array for the chemical potentials
                yarr(nv)=ceq%cmuval(nv)
             enddo
! after loop nv=noel()
             call get_component_name(nv,elsym,ceq)
             write(kou,2088)elsym,totam
2088         format('Mole fraction of ',a,' set to ',F8.5)
             xknown(nv)=totam
             yarr(nv)=ceq%cmuval(nv)
! use current T and P
             if(kom3.eq.4) then
! constituent adjustment, the FALSE means not quiet
                call equilph1b(phtup,ceq%tpval,xknown,xxx,yarr,.FALSE.,ceq)
                if(gx%bmperr.ne.0) goto 990
                write(kou,2087)xxx,(yarr(nv),nv=1,noel())
2087            format(/'Calculated Gibbs energy/FU/RT: ',1pe14.6,&
                     ' and the chemical potentials/RT:'/6(1pe12.4))
             else
!.............................................
! calculate phase diffusion: chem.pot derivatives and mobilities
! mugrad(I,J) are derivatives of the chemical potential of endmember I
!         with respect to endmember J
! mobilities(i) is mobility of component i
                mugrad=zero
                mobilities=zero
! derivatives of mu and mobilities, the FALSE means not quiet
                call equilph1d(phtup,ceq%tpval,xknown,yarr,.FALSE.,&
                     nend,mugrad,mobilities,ceq)
                if(gx%bmperr.ne.0) goto 990
                write(kou,2096)nend
2096            format(/'Chemical potential derivative matrix, dG_I/dn_J for ',&
                     i3,' endmembers')
                write(kou,2094)(nv,nv=1,nend)
2094            format(3x,6(6x,i6)/(3x,6i12))
                do nv=0,nend-1
! An extra LF is generated when just 6 components!! use ll, kp j1, i2
                   write(kou,2095)nv+1,(mugrad(nend*nv+jp),jp=1,nend)
!2095               format(i3,6(1pe12.4)/(3x,6e12.4))
2095               format(i3,6(1pe12.4)/(3x,6e12.4))
                enddo
                write(kou,2098)noel()
2098            format(/'Mobility values mols/m2/s ?? for',i3,' components')
                write(kou,2095)1,(mobilities(jp),jp=1,noel())
             endif
!.......................................................
          case(6) ! Quit
!             write(*,*)'Not implemeneted yet'
          END SELECT calcphase
! set bits to warn that listings may be inconsistent
          ceq%status=ibclr(ceq%status,EQNOEQCAL)
          ceq%status=ibset(ceq%status,EQINCON)
!---------------------------------- end of calculate phase
       case(3) ! calculate equilibrium without initial global minimization
          if(minimizer.eq.1) then
! Lukas minimizer, first argument=0 means do not use grid minimizer
!           call calceq1(0,ceq)
             write(kou,*)'Not implemented yet'
          else
             call calceq2(0,ceq)
! check that invmat allocated and stored
!           write(*,*)'inverted y: ',ceq%phase_varres(2)%cinvy(1,1)
          endif
          if(gx%bmperr.ne.0) then
             ceq%status=ibset(ceq%status,EQFAIL)
             goto 990
          endif
!----------------------------------
       case(4) ! calculate transition
          call calctrans(cline,last,ceq)
! clear this bit so there there will be no warning the listing is inconsistent
          if(gx%bmperr.ne.0) goto 990
          ceq%status=ibclr(ceq%status,EQINCON)
!----------------------------------
       case(5) ! quit
          goto 100
!-----------------------------------------------------------
       case(6) ! calculate global grid minimum
! extract values for mass balance calculation from conditions
          call extract_massbalcond(ceq%tpval,xknown,totam,ceq)
          if(gx%bmperr.ne.0) goto 990
! debug output
!          write(*,2101)totam,(xknown(j1),j1=1,noel())
!2101      format('UI N&x: ',F6.3,9F8.5)
! generate grid and find the phases and constitutions for the minimum.
! Note: global_gridmin calculates for total 1 mole of atoms, not totam
!          call global_gridmin(1,ceq%tpval,totam,xknown,nv,iphl,icsl,&
! iphl is dimensioned (1:maxel), maxel=100, it is destroyed inside merge_grid ..
!          call global_gridmin(1,ceq%tpval,xknown,nv,iphl,icsl,&
!               aphl,nyphl,yarr,cmu,iphl,ceq)
          call global_gridmin(1,ceq%tpval,xknown,nv,iphl,icsl,&
               aphl,nyphl,cmu,ceq)
          if(gx%bmperr.ne.0) goto 990
!          write(kou,2102)nv,(iphl(j1),icsl(j1),j1=1,nv)
! we should write phase tuples ... ?? 
          write(kou,2102)nv,(iphl(j1),icsl(j1),j1=1,nv)
2102      format('Number of stable phases ',i2/13(i4,i2))
! In some cases "c n" converges better if we scale with the total amount here
          do j1=1,nv
             call get_phase_compset(iphl(j1),icsl(j1),lokph,lokcs)
             ceq%phase_varres(lokcs)%amfu=totam*ceq%phase_varres(lokcs)%amfu
          enddo
! if set clear this bit so we can list the equilibrium
          if(btest(ceq%status,EQNOEQCAL)) ceq%status=ibclr(ceq%status,EQNOEQCAL)
!2103      format('Stable phase ',2i4,': ',a)
!---------------------------------------------------------------
       case(7) ! calculate symbol
!          call evaluate_all_svfun(kou,ceq)
! to calculate derivatives this must be in the minimizer module
          call gparcd('Name ',cline,last,1,name1,'*',q1help)
! always calculate all state variable functions as they may depend on eachother
!          write(*,*)'UI: calculating all functions'
          call meq_evaluate_all_svfun(-1,ceq)
! ignore error
          if(gx%bmperr.ne.0) then
             write(*,*)'UI: error calculating all functions', gx%bmperr
             gx%bmperr=0
          endif
          if(name1(1:1).eq.'*') then
! this calculate them again ... and lists the values
             call meq_evaluate_all_svfun(lut,ceq)
          else
! This code is also used in SHOW (command 25)
             call capson(name1)
!             call find_svfun(name1,istv,ceq)
             call find_svfun(name1,istv)
             if(gx%bmperr.ne.0) goto 990
             mode=1
             actual_arg=' '
!             write(*,*)'UI: calculating the requested function!'
             xxx=meq_evaluate_svfun(istv,actual_arg,mode,ceq)
             if(gx%bmperr.ne.0) goto 990
             write(*,2047)name1(1:len_trim(name1)),xxx
2047         format(a,'= ',1pe16.8)
          endif
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(8) ! calculate equilibrium for current equilibrium ceq
! using the grid minimizer
          if(minimizer.eq.1) then
! Lukas minimizer, first argument=1 means use grid minimizer
!           call calceq1(1,ceq)
             write(kou,*)'No longer available'
          else
             call calceq2(1,ceq)
             if(gx%bmperr.eq.4204) then
! if the error code is "too many iterations" try without grid minimizer
! it converges in many cases
                write(*,2048)gx%bmperr
2048            format('Error ',i5,', cleaning up and trying harder')
                gx%bmperr=0
                call calceq2(0,ceq)
             endif
          endif
! calceq2 set appropriate bits for listing
          if(gx%bmperr.ne.0) then
             ceq%status=ibset(ceq%status,EQFAIL)
             goto 990
          endif
!---------------------------------------------------------------
       case(9) ! calculate all equilibria
! rather complex to handle both parallel on non-parallel and with/without 
! griminimizer ...
          if(allocated(firstash%eqlista)) then
             call gparcd('With global minimizer? ',cline,last,1,ch1,'N',q1help)
! mode=0 is without grid minimizer ?? mode=-1 ??
             mode=1
             if(ch1.eq.'N' .or. ch1.eq.'n') mode=0
!             if(ch1.eq.'N' .or. ch1.eq.'n') mode=-1
! Seach for memory leaks
             call gparid('How many times? ',cline,last,leak,1,q1help)
! leak<0 means forever ... not allowed
             if(leak.lt.1) leak=1
! Minimize output
             listzero=.false.
!             call gparcd('List zero weight? ',cline,last,1,chz,'Y',q1help)
!             if(chz.eq.'N') listzero=.true.
! allow output file, if idef>1 no output
             idef=leak
             lut=optionsset%lut
             jp=0
             i2=0
! if compiled with parallel and gridminimizser set then calculate
! sequentially to create composition sets
! TEST THIS IN PARALLEL !!!
             call cpu_time(xxx)
             call system_clock(count=j1)
             threads=1
! OPENMP parallel start
!$             threads=omp_get_num_threads()
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! for parallelizing:
! YOU MUST UNCOMMENT USE OMP_LIB IN GTP3.F90 or PMON6.F90
! YOU MUST USE THE SWICH -fopenmp FOR COMPILATION AND WHEN LINKING
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!             gridmin: if(mode.eq.1) then
!
! return here until leak is zero, a negative leak will never stop
2060         continue
             gridmin: if(mode.eq.1 .or. &
                  btest(globaldata%status,GSNOPAR)) then
! if we use grid minimizer do not use parallel even if compiled with OpenMP
                do i1=1,size(firstash%eqlista)
                   neweq=>firstash%eqlista(i1)%p1
                   jp=jp+1
                   if(neweq%weight.eq.zero) then
                      if(listzero) write(kou,2050)neweq%eqno,neweq%eqname
2050                  format('Zero weight equilibrium ',i4,2x,a)
                   else
                      i2=i2+1
                      call calceq3(mode,.FALSE.,neweq)
                      if(gx%bmperr.ne.0) then
                         write(kou,2051)gx%bmperr,neweq%eqno,neweq%eqname
2051                     format(' *** Error code ',i5,' for equilibrium ',&
                              i5,': ',a,' reset ')
                         gx%bmperr=0
                      elseif(idef.eq.1) then
! extract names of stable phases
                         jp=1
                         line=' '
                         do j3=1,nooftup()
                            phtup=>phasetuple(j3)
         if(neweq%phase_varres(phtup%lokvares)%phstate.ge.PHENTSTAB) then
             call get_phasetup_name(j3,line(jp:))
             jp=len_trim(line)+2
          endif
                         enddo
                         write(lut,2052)neweq%eqno,&
                              neweq%eqname(1:len_trim(neweq%eqname)),&
                              neweq%tpval(1),trim(line)
2052                     format(i5,2x,a,', T=',F8.2,', ',a)
                      endif
                   endif
! extra symbol calculations ....
!                   write(*,*)'Listing extra'
                   if(idef.eq.1) then
                      call list_equilibrium_extra(lut,neweq,plotunit0)
                      if(gx%bmperr.ne.0) then
                         write(kou,*)'Error ',gx%bmperr,' reset'
                         gx%bmperr=0
                      endif
                   endif
                enddo
             else
! Here we calculate without grid minimizer, if parallel we must turn off
! creating/removing composition sets!! not safe to do that!!
!$             globaldata%status=ibset(globaldata%status,GSNOACS)
!$             globaldata%status=ibset(globaldata%status,GSNOREMCS)
!        !$OMP for an OMP directive
!        !$ as "sentinel", will be compiled if -fopenmp
! this statement must not be inside a parallel do ...
                svss=size(firstash%eqlista)
! NOTE: $OMP  threadprivate(gx) declared in TPFUN4.F90 ??
!$OMP parallel do private(neweq)
                paraloop: do i1=1,svss
!                do i1=1,size(firstash%eqlista)
! the error code must be set to zero for each thread ?? !!
                   jp=jp+1
                   gx%bmperr=0
                   neweq=>firstash%eqlista(i1)%p1
! it seems stupid to get this value each loop but outside it is unity
!$                   threads=omp_get_num_threads()
                   if(neweq%weight.eq.zero) then
                      if(listzero) write(kou,2050)neweq%eqno,neweq%eqname
                   else
! write output only for idef=1
!$                     if(.TRUE. .and. idef.eq.1) then
!$                      write(*,663)'Equil/loop/thread/maxth/error: ',&
!$                             neweq%eqname,i1,omp_get_thread_num(),&
!$                             threads,gx%bmperr
663                   format(a,a,5i5)
! calceq3 gives no output
!$                        call calceq3(mode,.FALSE.,neweq)
!$                     else
! note first argument zero means do not use grid minimizer
                          call calceq3(mode,.FALSE.,neweq)
!$                     endif
                      i2=i2+1
                      if(gx%bmperr.ne.0) then
                         write(kou,2051)gx%bmperr,neweq%eqno,neweq%eqname
!                         write(*,*)'Error: ',gx%bmperr
                         gx%bmperr=0
                      elseif(idef.eq.1) then
! extract names of stable phases
                         jp=1
                         line=' '
                         do j3=1,nooftup()
                            phtup=>phasetuple(j3)
         if(neweq%phase_varres(phtup%lokvares)%phstate.ge.PHENTSTAB) then
             call get_phasetup_name(j3,line(jp:))
             jp=len_trim(line)+2
          endif
                         enddo
                         write(lut,2052)neweq%eqno,&
                              neweq%eqname(1:len_trim(neweq%eqname)),&
                              neweq%tpval(1),trim(line)
                      endif
                   endif
! Listing extra'
                   if(idef.eq.1) then
                      call list_equilibrium_extra(lut,neweq,plotunit0)
                      if(gx%bmperr.ne.0) then
!                         write(kou,*)'*** Error ',gx%bmperr,' reset'
                         gx%bmperr=0
                      endif
                   endif
                enddo paraloop
!- $OMP end parallel do not needed???
! OPENMP parallel end loop
! allow composition sets to be created again
!$             globaldata%status=ibclr(globaldata%status,GSNOACS)
!$             globaldata%status=ibclr(globaldata%status,GSNOREMCS)
             endif gridmin
! repeat this until leak is zero, if leak negative never stop.
             leak=leak-1
             call cpu_time(xxz)
             if(leak.ne.0) then
                call system_clock(count=ll)
                xxy=ll-j1
! or should i2 be used ??
                write(*,669)i2,(xxz-xxx)/i2,xxy/i2
669        format(/'Calculated ',i8,' equlibria, average CPU and clock time',&
                F12.8,F10.6)
                goto 2060
             endif
!
             call system_clock(count=ll)
! ?? jp ??             write(kou,664)jp,xxz-xxx,ll-j1,threads
             write(kou,664)xxz-xxx,ll-j1,threads
!664          format('Calculated equilibria out of ',i5/&
664          format('Total CPU time: ',1pe12.4,' s and ',i7,' clockcycles',&
                  ' using ',i4,' thread(s)')
! this unit may have been used to extract calculated data for plotting
             if(plotunit0.gt.0) then
                write(kou,670)
670             format('Closing a GNUPLOT file oc_many0.plt'/&
                     'that may need some editing before plotting')
                write(plotunit0,665)graphopt%plotend
665             format('e'/a)
!665             format('e'/'pause mouse'/)
                close(plotunit0)
! UNFINISHED possibly we could reopen the file again and make oopies 
! of tha data to avoid manual editing
             endif
          else
             write(kou,*)'You must first SET RANGE of experimental equilibria'
          endif
!---------------------------------------------------------------
       case(10) ! calculate with_check_after
! this is same as "calculate no_grid" but with automatic grid check after
! we must set global bit 18 and then clear it
! If bit 20 is set we will clear it now and set it afterwards
          if(btest(globaldata%status,GSNORECALC)) then
             globaldata%status=ibclr(globaldata%status,GSNORECALC)
             temporary=.TRUE.
          else
             temporary=.FALSE.
          endif
          globaldata%status=ibset(globaldata%status,GSTGRID)
! calculate with no grid before but check after
          call calceq2(0,ceq)
          if(gx%bmperr.ne.0) then
             ceq%status=ibset(ceq%status,EQFAIL)
          endif
! reset bit GSTGRID and maybe GSNORECALC
          globaldata%status=ibclr(globaldata%status,GSTGRID)
          if(temporary) &
               globaldata%status=ibset(globaldata%status,GSNORECALC)
!-------------------------------------------------------
       case(11) ! not used (yet)
!-------------------------------------------------------
       case(12) ! not used (yet)
       END SELECT calculate
!=================================================================
! SET SUBCOMMANDS
!         ['CONDITION       ','STATUS          ','ADVANCED        ',&
!         'LEVEL           ','INTERACTIVE     ','REFERENCE_STATE ',&
!         'QUIT            ','ECHO            ','PHASE           ',&
!         'UNITS           ','LOG_FILE        ','WEIGHT          ',&
!         'NUMERIC_OPTIONS ','AXIS            ','INPUT_AMOUNTS   ',&
!         'VERBOSE         ','AS_START_EQUILIB','BIT             ',&
!         'VARIABLE_COEFF  ','SCALED_COEFF    ','OPTIMIZING_COND ',&
!         'RANGE_EXPER_EQU ','FIXED_COEFF     ','SYSTEM_VARIABLE ',&
!         'INITIAL_T_AND_P ','                ','                ']
    CASE(3) ! SET SUBCOMMANDS
! disable continue optimization
!       iexit=0
!       iexit(2)=1
       kom2=submenu(cbas(kom),cline,last,cset,ncset,1)
       if(kom2.le.0) goto 100
       set: SELECT CASE(kom2)
       CASE DEFAULT
          write(kou,*)'Set subcommand error'
!-----------------------------------------------------------------------
       CASE(1) ! set condition
          if(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'You have no data!'
             goto 100
          endif
          call set_condition(cline,last,ceq)
!------------------------------------------------------------------
       CASE(2) ! set status for elements, species, phases, constituents
          name1='STATUS of'
          kom3=submenu(name1,cline,last,cstatus,ncstat,3)
          setstatus: SELECT CASE(kom3)
!.................................................................
          CASE DEFAULT
             write(kou,*)'Set status subcommand error'
!.................................................................
          case(1) ! set status element suspend/restore
             call gparc('Element symbol: ',cline,last,1,name1,' ',q1help)
             call find_element_by_name(name1,iel)
             if(gx%bmperr.ne.0) goto 100
             call gparcd('New status: ',cline,last,1,ch1,'SUSPEND',q1help)
             call capson(ch1)
             if(ch1.eq.'S') then
                call change_element_status(name1,1,ceq)
             else
! restore element
                call change_element_status(name1,0,ceq)
             endif
!.................................................................
          case(2) ! set status species suspend/restore
             call gparc('Species symbol: ',cline,last,1,name1,' ',q1help)
             call find_species_record(name1,loksp)
             if(gx%bmperr.ne.0) goto 100
             call gparcd('New status: ',cline,last,1,ch1,'SUSPEND',q1help)
             call capson(ch1)
             if(ch1.eq.'S') then
                call change_species_status(name1,1,ceq)
             else
                call change_species_status(name1,0,ceq)
             endif
!.................................................................
          case(3) ! set status phase (ENTERED, FIX, DORMANT, SUSPEND or HIDDEN)
! Now allow multiple phase names and *S, *D and *E
!             call gparc('Phase name(s): ',cline,last,1,name1,' ',q1help)
! argument 5 means whole input line
             call gparc('Phase name(s): ',cline,last,5,line,'=',q1help)
             string=line
3017         continue
             ll=index(string,'=')
             if(ll.eq.0) then
                call gparc('More phase name(s): ',cline,last,5,line,'=',q1help)
                string(len_trim(string)+2:)=line
                goto 3017
             endif
!3018         continue
! exttract first letter after = (if any)
             j1=ll
             call getext(string,j1,1,name1,' ',iph)
             ch1=name1(1:1)
! if user has given "=e 0" then keep the amount is cline
             cline=string(j1:)
             string(ll:)=' '
!             write(*,*)'s1: ',j1,cline(1:len_trim(cline))
             if(ch1.eq.' ') then
! if ll==1 then input was finished by equal sign, ask for status
                call gparcd(&
                     'New status S(uspend), D(ormant), E(ntered) or F(ixed)?',&
                     cline,last,1,name1,'E',q1help)
                ch1=name1(1:1)
             else
                last=0
             endif
             nystat=99
             call capson(ch1)
! new values of status ??
             if(ch1.eq.'S') nystat=phsus
             if(ch1.eq.'D') nystat=phdorm
             if(ch1.eq.'E') nystat=phentered
             if(ch1.eq.'F') nystat=phfixed
!             if(ch1.eq.'H') nystat=phhidden
! no longer available if(ch1.eq.'N') nystat=5
             if(nystat.eq.99) then
                write(kou,*)'No such status'
                goto 100
             endif
             xxx=zero
!             write(*,*)'s2: ',last,cline(1:len_trim(cline))
             if(nystat.eq.phentered .or. nystat.eq.phfixed) then
                call gparrd('Amount: ',cline,last,xxx,zero,q1help)
             endif
             call change_many_phase_status(string,nystat,xxx,ceq)
             if(gx%bmperr.ne.0) goto 100
!.................................................................
          CASE(4) ! set status constituent 
             write(kou,*)'Not implemented yet'
!.................................................................
          case(5) ! set status subcommand status for ?
             write(kou,*)'Not implemented yet'
!.................................................................
          case(6) ! set status subcommand status for ?
             write(kou,*)'Not implemented yet'
          END SELECT setstatus
!-----------------------------------------------------------
       case(3) ! set ADVANCED
! default is DENSE_GRID
          name1='advanced command'
          kom3=submenu(name1,cline,last,cadv,ncadv,4)
          advanced: select case(kom3)
!.................................................................
          CASE DEFAULT
             write(kou,*)'Set advanced subcommand error'
!.................................................................
! SET ADVANCED EQUILIB_TRANSFER
! transfer a ceq record from map results%savedceq to eqlista
! so it can be used interactivly
          case(1)
             if(.not.associated(maptop)) then
                write(kou,*)'There are no results from map or step'
                goto 100
             else
                write(kou,3100)maptop%saveceq%free-1
3100            format('Saved ceq records from 1 to ',i3) 
             endif
             write(kou,*)'To transfer CEQ from result area to current'
             call gparid('Saved ceq number',cline,last,icon,1,q1help)
             if(icon.gt.0 .and. icon.lt.maptop%saveceq%free) then
                name1='COPIED_RESULTS_'
                i2=len_trim(name1)+1
                call wriint(name1,i2,icon)
                write(*,*)'Equilibrium name: ',i2,': ',name1
                call enter_equilibrium(name1,i1)
                if(gx%bmperr.ne.0) goto 990
                write(*,*)'Created equilibrium ',i1
! note... this will overwrite the name ...
                eqlista(i1)=maptop%saveceq%savedceq(icon)
! maybe not possible, eqlista is maybe protected ... no it is not
                write(*,*)'Trying to change name ...'
                eqlista(i1)%eqname=name1
                call selecteq(i1,ceq)
             else
                write(kou,*)'No such saved equilibrium'
             endif
! set bit that data may be inconsistent
             eqlista(i1)%status=ibset(eqlista(i1)%status,EQINCON)
!.................................................................
          case(2) ! quit
             continue
!.................................................................
          case(3) ! SET ADVANCED unused
             write(*,*)'Not implemented'
!.................................................................
          case(4) ! SET ADVANCED GRID_DENSITY
             call gparid('Level: ',cline,last,ll,1,q1help)
             if(ll.eq.0) then
! this set GSOGRID, small grid and clears GSXGRID
                globaldata%status=ibset(globaldata%status,GSOGRID)
                globaldata%status=ibclr(globaldata%status,GSXGRID)
                write(kou,3110)'Sparse','set'
             elseif(ll.eq.1) then
! this clears GSXGRID, bit 14, of global status word and GSOGRID
                globaldata%status=ibclr(globaldata%status,GSXGRID)
                globaldata%status=ibclr(globaldata%status,GSOGRID)
                write(kou,3110)'Normal','set'
3110            format(a,' grid ',a)
             elseif(ll.eq.2) then
! set GSXGRID (and clear GSOGRID)
                globaldata%status=ibclr(globaldata%status,GSOGRID)
                globaldata%status=ibset(globaldata%status,GSXGRID)
                write(kou,3110)'Dense','set'
             else
                write(*,*)'Only level 0, 1 and 2 implemented'
             endif
!.................................................................
          case(5) ! SET ADVANCED SMALL_GRID_ONOFF
! replaced by setting grid_density to 0
             write(*,*)'Please use SET ADVANCED GRID 0'
             continue
!             if(btest(globaldata%status,GSOGRID)) then
!                globaldata%status=ibclr(globaldata%status,GSOGRID)
!                write(kou,3110)'Small','reset'
!             else
! set GSOGRID and clear GSXGRID if set
!                globaldata%status=ibclr(globaldata%status,GSXGRID)
!                globaldata%status=ibset(globaldata%status,GSOGRID)
!                write(kou,3110)'Small','set'
!             endif
!.................................................................
          case(6) ! MAP_SPECIAL
!             if(nofixphfortip) then
!                write(*,*)'Always using fix phase when mapping'
!                nofixphfortip=.false.
!             else
!                write(*,*)'Map diagrams with tie-lines in phase ',&
!                     'without fix phase'
!                nofixphfortip=.true.
!             endif
             write(*,*)'Not implemented yet'
!.................................................................
          case(7) ! GLOBAL_MIN_ONOFF
             call gparc('Turn global minimization off?: ',cline,last,&
                  1,ch1,'N',q1help)
             if(ch1.eq.'Y' .or. ch1.eq.'y') then
                globaldata%status=ibset(globaldata%status,GSNOGLOB)
                write(*,*)'Global minimizer turned off'
             else
                globaldata%status=ibclr(globaldata%status,GSNOGLOB)
                write(*,*)'Global minimizer turned on'
             endif
!             if(btest(globaldata%status,GSNOGLOB)) then
!                globaldata%status=ibclr(globaldata%status,GSNOGLOB)
!                write(*,*)'Global minimizer turned on'
!             else
!                globaldata%status=ibset(globaldata%status,GSNOGLOB)
!                write(*,*)'Global minimizer turned off'
!             endif
!             write(*,*)'Not implemented yet'
!.................................................................
          case(8) ! OPEN_POPUP_OFF
             call gparcd('Turn off popup for open? ',cline,last,&
                  1,ch1,'Y',q1help)
             if(ch1.eq.'Y') then
! nopopup is declared in metlib3.F90 module
! nopenpopup is declared in metlib3.F90 module
                nopenpopup=.TRUE.
                write(kou,*)'Popup windows for open files turned off'
             else
                nopenpopup=.FALSE.
                write(kou,*)'Popup windows for open files enabled'
             endif
!.................................................................
          case(9) ! WORKING DIRECTORY
             write(kou,*)'Current working directory: ',trim(workingdir)
             write(kou,*)'To change please give full path'
             call gparc('New: ',cline,last,1,string,workingdir,q1help)
             inquire(file=string,exist=logok)
             if(.not.logok) then
                write(*,*)'No such directory'
             elseif(trim(workingdir).ne.trim(string)) then
                write(*,'(a,a)')'Working directory set to: ',trim(string)
                workingdir=string
             endif
!             write(*,*)'Cannot be changed'
!.................................................................
          case(10) ! HELP_POPUP_OFF
             call gparcd('Turn off popup help? ',cline,last,&
                  1,ch1,'Y',q1help)
             if(ch1.eq.'Y') then
                ochelp%htmlhelp=.FALSE.
                htmlhelp=.FALSE.
             else
                htmlhelp=.TRUE.
                string=browser
                call gparcd('Browser including full path ',&
                     cline,last,1,browser,string,q1help)
                string=latexfile
                call gparcd('LaTeX help file including full path ',&
                     cline,last,1,latexfile,string,q1help)
                string=htmlfile
                call gparcd('HTML help file including full path ',&
                     cline,last,1,htmlfile,string,q1help)
                call init_help(browser,latexfile,htmlfile)
             endif
!.................................................................
          case(11) ! SET ADVANCED EET_EXTRAPOL / HICKEL_EXTRAPOL for solids
             call gparcd('Turn on equi-entropy extrapolation (EET)?',&
                  cline,last,1,ch1,'Y',q1help)
             if(ch1.eq.'Y') then
! remove this bit for EET, use only sysparem(1) for T (as integer)
!                globaldata%status=ibset(globaldata%status,GSHICKEL)
                call gparrd('Low T limit?',cline,last,xxx,1.0D3,q1help)
!                if(thickel.lt.one) then
                if(xxx.gt.1.0D1) then
! set_hickel_check is in minimizer/matsmin.F90
                   call set_hickel_check(xxx)
!                   globaldata%sysreal(1)=1.0D1
                endif
             else
! remove the use of this bit for EET
!                globaldata%status=ibclr(globaldata%status,GSHICKEL)
                write(*,*)'EET extrapolation for solids turned off'
                call set_hickel_check(zero)
!                globaldata%sysreal(1)=zero
             endif
!.................................................................
          case(12) ! SET ADVANCED LEVEL
             call gparcd('I am an expert of OC: ',cline,last,1,ch1,'N',q1help)
             if(ch1.eq.'Y') then
                globaldata%status=ibset(globaldata%status,2)
                write(*,*)'Felicitations!'
             else
                write(*,*)'Sorry, not yet'
             endif
!.................................................................
          case(13) ! not used
!.................................................................
          case(14) ! not used
!.................................................................
          case(15) ! not used
          end select advanced
!-----------------------------------------------------------
       case(4) ! set LEVEL, MOVED TO SET ADBANCED
! 
          write(*,*)'Unused'
!-----------------------------------------------------------
! end of macro excution (can be nested)
       case(5) ! set INTERACTIVE
          call macend(cline,last,logok)  
! if this was the startupmacro set it false and possibly read an inline macro ..
! NOTE a startup macro can call other macros ...
          if(kiu.eq.kiud) startupmacro=.false.
!          macropath=' '
!          write(*,*)'Macro terminated'
!-----------------------------------------------------------
       case(6) ! set REFERENCE_STATE
          call gparc('Component name: ',cline,last,1,name1,' ',q1help)
          call find_component_by_name(name1,iel,ceq)
          if(gx%bmperr.ne.0) goto 100
          call gparc('Reference phase: ',cline,last,1,name1,'SER ',q1help)
          if(name1(1:4).eq.'SER ') then
             write(kou,*)'Reference state is stable phase at 298.15 K and 1 bar'
! this means no reference phase, SER is at 298.15K and 1 bar
             iph=-1
          else
             call find_phase_by_name(name1,iph,ics)
             if(gx%bmperr.ne.0) goto 100
! temperature * means always to use current temperature
             xxy=-one
             call gparr('Temperature: /*/: ',cline,last,xxx,xxy,q1help)
!             write(*,*)'problem: ',buperr,xxx,xxy,one
! when calling gparr the default was not "set" as default and rubbish returned
! now the default is always the default even if not shown
             if(buperr.ne.0) then
                buperr=0
                tpa(1)=-one
             elseif(xxx.le.zero) then
                tpa(1)=-one
             else
                tpa(1)=xxx
             endif
             xxy=1.0D5
             call gparrd('Pressure: ',cline,last,xxx,xxy,q1help)
             if(xxx.le.zero) then
                tpa(2)=xxy
             else
                tpa(2)=xxx
             endif
          endif
          call set_reference_state(iel,iph,tpa,ceq)
          if(gx%bmperr.eq.0) then
!             write(kou,3104)
3104         format(' You may have to make a new calculation before the',&
                  ' correct values'/&
                  ' of chemical potentials or other properties are shown.')
          endif
!-----------------------------------------------------------
       case(7) ! quit
          goto 100
!-----------------------------------------------------------
       case(8) ! set ECHO
          call gparcd('On?',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'Y' .or. ch1.eq.'y') then
             j1=1
          else
             j1=0
          endif
          call set_echo(j1)
!-----------------------------------------------------------
       case(9) ! set PHASE subcommands (constitution, status)
          call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) then
             if(name1(1:2).eq.'* ') then
                iph=-1
                gx%bmperr=0
             else
                goto 990
             endif
          endif
          kom3=submenu(cbas(kom),cline,last,csetph,nsetph,2)
          setphase: SELECT CASE(kom3)
          CASE DEFAULT
             write(kou,*)'Set phase status subcommand error'
             goto 100
!............................................................
          case(1) ! quit
             continue
!............................................................
! begin code copied from 3045
          case(2) ! SET PHASE STATUS <phase> <status>
             if(iph.gt.0) then
                j1=get_phase_status(iph,ics,text,i1,xxx,ceq)
                if(gx%bmperr.ne.0) goto 100
                if(xxx.ge.zero) then
                   write(kou,3046)text(1:i1),xxx
3046               format('Current status is ',a,' with ',1pe15.6,&
                        ' formula units.')
                else
                   write(kou,3047)text(1:i1)
3047               format('Current status is ',a)
                endif
             endif
             call gparcd(&
                  'Suspend, Dormant, Entered, Fixed, Hidden or Not hidden?',&
                  cline,last,1,ch1,'SUSPEND',q1help)
             nystat=99
             call capson(ch1)
! new values of status ??
             if(ch1.eq.'E') nystat=phentered
             if(ch1.eq.'S') nystat=phsus
             if(ch1.eq.'D') nystat=phdorm
             if(ch1.eq.'F') nystat=phfixed
             if(ch1.eq.'H') nystat=phhidden
! not avail if(ch1.eq.'N') nystat=5
             if(nystat.eq.99) then
                write(kou,*)'No such status'
                goto 100
             endif
             xxx=zero
             if(nystat.eq.phentered .or. nystat.eq.phfixed) then
                call gparrd('Amount: ',cline,last,xxx,zero,q1help)
             endif
             call change_phase_status(iph,ics,nystat,xxx,ceq)
             if(gx%bmperr.ne.0) goto 100
             if(iph.gt.0) then
                j1=get_phase_status(iph,ics,text,i1,xxy,ceq)
                if(gx%bmperr.ne.0) goto 100
                if(xxy.ge.zero) then
                   write(kou,3048)text(1:i1),xxy
3048               format('New status is ',a,' with ',1pe15.6,&
                        ' formula units.')
                else
                   write(kou,3049)text(1:i1)
3049               format('New status is ',a)
                endif
             else
                write(kou,*)'New status set for all phases'
             endif
! end code copied from 3045
!............................................................
          case(3:4) !set phase default_constitution wildcard allod, also AMOUNT
!             write(*,*)'SET PHASE AMOUNT or DEFAULT_CONST',kom3,iph,ics
             if(kom3.eq.3) then
! set default constituntion of phase
!                call set_default_constitution(iph,ics,ceq)
                call ask_default_constitution(cline,last,iph,ics,ceq)
             else
! set phase amount
                call gparrd('Amount: ',cline,last,xxx,zero,q1help)
                call set_phase_amounts(iph,ics,xxx,ceq)
             endif
!............................................................
! subsubsub command
          case(5) ! set phase bits
             if(iph.lt.0) then
                write(kou,*)'Wildcards not allowed in this case'
                goto 100
             endif
             call get_phase_record(iph,lokph)
             kom4=submenu('Set which bit?',cline,last,csetphbits,nsetphbits,9)
             phasebit: SELECT CASE(kom4)
             CASE DEFAULT
! allow any bit changes for experts ...
                if(btest(globaldata%status,GSADV)) then
                   call getint(cline,last,ll)
                   if(ll.ge.0 .and. ll.le.31) then
                      write(kou,*)'Ahh, you are an expert ... changing bit: ',ll
                      if(test_phase_status_bit(lokph,ll)) then
                         call clear_phase_status_bit(lokph,ll)
                         write(kou,*)'Clearing bit ',ll
                      else
                         call set_phase_status_bit(lokph,ll)
                      endif
                   else
                      write(kou,*)'Illegal bit number',ll
                   endif
                else
                   write(kou,*)'Set phase bit subcommand error'
                endif
!............................................................
             case(1) ! FCC_PERMUTATIONS FORD
! if check returns .true. it is not allowed to set FORD or BORD
                if(check_minimal_ford(lokph)) goto 100
                write(*,*)' *** WARNING: Depreceated command, use AMEND PHASE'
                call set_phase_status_bit(lokph,PHFORD)
             case(2) ! BCC_PERMUTATIONS BORD
                if(check_minimal_ford(lokph)) goto 100
                write(*,*)' *** WARNING: Depreceated command, use AMEND PHASE'
                call set_phase_status_bit(lokph,PHBORD)
             case(3) ! IONIC_LIQUID_MDL this may require tests and 
! other bits changed ..
                write(*,*)' *** WARNING: set by enter phase <name> I2SL'
!                write(kou,*)'Cannot be set interactivly yet, only from TDB'
!                call set_phase_status_bit(lokph,PHIONLIQ)
             case(4) ! AQUEOUS_MODEL   
                write(*,*)'Not implemented yet'
!                call set_phase_status_bit(lokph,PHAQ1)
             case(5) ! QUASICHEMICAL   
                write(*,*)'Not implemented yet'
!                call set_phase_status_bit(lokph,PHQCE)
             case(6) ! FCC_CVM_TETRADRN
                write(*,*)'Not implemented yet'
!                call set_phase_status_bit(lokph,PHCVMCE)
             case(7) ! FACT_QUASICHEMCL
                write(*,*)'Not implemented yet'
!                call set_phase_status_bit(lokph,PHFACTCE)
             case(8) ! NO_AUTO_COMP_SET, do not create compsets automatically
                call set_phase_status_bit(lokph,PHNOCS)
             case(9) ! QUIT
                write(kou,*)'No other bits changed'
             case(10) ! EXTRA_DENSE_GRID, this can be toggled ...
                if(test_phase_status_bit(lokph,PHXGRID)) then
                   write(kou,*)'Bit already set, is cleared'
                   call clear_phase_status_bit(lokph,PHXGRID)
                else
                   write(kou,*)'Extra gridpoints for this phase.'
                   call set_phase_status_bit(lokph,PHXGRID)
                endif
             case(11) ! Flory-Huggins polymer model
                call set_phase_status_bit(lokph,PHFHV)
                call clear_phase_status_bit(lokph,PHID)
             end SELECT phasebit
!............................................................
          case(6) ! SET PHASE ... CONSTITUTION iph and ics set above
             call ask_phase_new_constitution(cline,last,iph,ics,lokcs,ceq)
          END SELECT setphase
!-------------------------------------------------------------
       case(10) ! set UNIT (for state variables)
          write(kou,*)'Not implemented yet'
!-------------------------------------------------------------
       case(11) ! set LOG_FILE
          call gparfile('Log file name: ',cline,last,1,name1,'oclog',8,q1help)
          call capson(name1)
          if(name1(1:5).eq.'NONE ') then
             call openlogfile(' ',' ',-1)
             logfil=0
          else
             call gparc('Title: ',cline,last,5,model,' ',q1help)
             call openlogfile(name1,model,39)
             if(buperr.ne.0) then
                write(kou,*)'Error opening logfile: ',buperr
                logfil=0
             else
                logfil=39
             endif
          endif
!-------------------------------------------------------------
       case(12) ! set weight
          if(.not.allocated(firstash%eqlista)) then
             write(kou,*)'You must first set a range of experimental equilibria'
             goto 100
          endif
! NOTE mexp must be updated to the correct number of EXPERIMENTS
! that is done by OPTIMIZE
          updatemexp=.true.
          mexp=0
          call gparrd('Weight ',cline,last,xxx,one,q1help)
          if(buperr.ne.0) goto 100
! The weight must be 0 or positive
          xxx=abs(xxx)
          call gparcd('Equilibria (abbrev name) or range: ',cline,last,&
               1,name1,'CURRENT',q1help)
! THINK HOW TO UPDATE MEXP!!! <<<<<<<<<<<<<<<<<<
          if(name1(1:8).eq.'CURRENT ') then
             if(ceq%eqname(1:20).eq.'DEFAULT_EQUILIBRIUM ') then
                write(kou,*)'You cannot set weight for the default equilibrium'
             else
                ceq%weight=xxx
             endif
          elseif(name1(1:1).eq.'*') then
! set this weight for all
             i2=0
             do i1=1,size(firstash%eqlista)
                firstash%eqlista(i1)%p1%weight=xxx
                i2=i2+1
             enddo
             write(kou,3066)i2
          else
             ll=1
!             write(*,*)'trying to extract a number from: ',trim(name1)
             call getint(name1,ll,i1)
             bupp: if(buperr.eq.0) then
! user provide a singe number or a range, if range the negative number also
                call getint(name1,ll,i2)
                if(buperr.ne.0) then
! it was a single number                   
                   buperr=0
                   i2=-i1
                endif
                i2=-i2
                ll=0
!                setwei: do j1=i1,i2
                setwei: do j1=1,size(firstash%eqlista)
                   if(firstash%eqlista(j1)%p1%eqno.ge.i1 .and. &
                        firstash%eqlista(j1)%p1%eqno.le.i2) then
                      firstash%eqlista(j1)%p1%weight=xxx
!                      write(*,*)'Changing weight for equilibrium ',&
!                           firstash%eqlista(j1)%p1%eqno
                      ll=ll+1
                   endif
                enddo setwei
                write(kou,3066)ll
             else
! set this weight to all equilibria with name abbriviations fitting name1
                buperr=0
                call capson(name1)
                if(name1(1:1).ne.' ') then
                   write(*,*)'Equilibra with names matching: ',trim(name1)
                   i2=0
                   do i1=1,size(firstash%eqlista)
                      if(index(firstash%eqlista(i1)%p1%eqname,&
                           name1(1:len_trim(name1))).gt.0) then
                         firstash%eqlista(i1)%p1%weight=xxx
                         i2=i2+1
                      endif
                   enddo
                else
                   write(*,*)'No name given'
                endif
                write(kou,3066)i2
             endif bupp
3066         format('Changed weight for ',i5,' equilibria')
          endif
!-------------------------------------------------------------
! turn on/off global minimization, creating composition sets
! convergence limits, iterations, minimum constituent fraction, etc
       case(13) ! set NUMERIC_OPTIONS
          i2=ceq%maxiter
          call gparid('Max number of iterations: ',cline,last,i1,i2,q1help)
          if(i1.gt.0) then
             ceq%maxiter=i1
          endif
!------------
          xxx=ceq%xconv
          call gparrd('Max error in fraction: ',cline,last,xxy,xxx,q1help)
          if(xxy.gt.1.0D-30) then
             ceq%xconv=xxy
          else
             ceq%xconv=1.0D-30
          endif
!------------ what is this? not used in gtp3X.F90
          xxx=ceq%gdconv(1)
          call gparrd('Max cutoff driving force: ',cline,last,xxy,xxx,q1help)
          if(xxy.gt.1.0D-5) then
             ceq%gdconv(1)=xxy
          endif
!------------ if the point between two gridpoints in a phase is less then merge
          xxx=ceq%gmindif
          call gparrd('Min difference merging gridpoints: ',cline,last,&
               xxy,xxx,q1help)
          if(xxy.gt.-1.0D-5) then
             ceq%gmindif=xxy
          else
             ceq%gmindif=-1.0D-2
          endif
!-------------------------------------------------------------
       case(14) ! set axis
          if(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'You have no data!'
             goto 100
          endif
          i1=noofaxis+1
          call gparid('Axis number',cline,last,iax,i1,q1help)
          if(iax.lt.1 .or. iax.gt.maxax) then
             write(kou,3300)maxax
3300         format('Axis number must be between 1 and ',i1)
             goto 100
          endif
! by giving a value of iax lesser than noofaxis one can change an already
! defined axis, values larger than i1 (=noofaxis+1) not allowed.
          if(iax.gt.i1) then
             iax=i1
             write(kou,*)'Axis must be set in sequential order',&
                  ', axis number set to ',iax
          endif
! as condition one may give a condition number followed by :
! or a single state variable like T, x(o) etc.
          if(iax.lt.i1) then
! set the current condition as default answer
             jp=1
             call get_one_condition(jp,name1,axarr(iax)%seqz,ceq)
             if(gx%bmperr.ne.0) goto 990
             jp=index(name1,'=')
             name1(jp:)=' '
! set current axis limits as default
             dmin=axarr(iax)%axmin
             dmax=axarr(iax)%axmax
          else
! new axis, defaults 0 and 1
             name1=' '
             dmin=zero
             dmax=one
          endif
! reset default plot options
          call reset_plotoptions(graphopt,plotfile,textlabel)
          axplotdef=' '
          call gparcd('Condition varying along axis: ',cline,last,1,&
               text,name1,q1help)
          call capson(text)
!          if(text(1:1).eq.' ') goto 100
          removeaxis: if(text(1:1).eq.' ' .or. text(1:4).eq.'NONE') then
! this means remove an axis, shift any higher axis down
             if(iax.lt.noofaxis) write(kou,*)'Shifting axis down'
             do i2=iax,noofaxis
                axarr(i2)=axarr(i2+1)
             enddo
             if(noofaxis.gt.1) then
                noofaxis=noofaxis-1
                write(kou,*)'One axis removed'
! remove axplotdef for all axis!!! one may change from PD to step sep
                axplotdef=' '
             endif
             goto 100
          else ! add or change axis variable
             i1=len_trim(text)
             if(text(i1:i1).eq.':') then
! condition given as an index in the condition list terminated by : like "1:"
                i1=1
                call getrel(text,i1,xxx)
                if(buperr.ne.0) then
                   gx%bmperr=buperr; goto 990
                endif
                i2=int(xxx)
                firstc=>ceq%lastcondition
                if(associated(firstc)) then
                   firstc=>firstc%next
                   pcond=>firstc%next
                   i1=0
                   do while(.not.associated(pcond,firstc))
! increment i1 only for active conditions as listed by list_condition
                      if(pcond%active.eq.0) i1=i1+1
                      if(i1.eq.i2) goto 3310
                      pcond=>pcond%next
                   enddo
                   gx%bmperr=4131; goto 990
3310               continue
                else
                   gx%bmperr=4131; goto 990
                endif
! pcond points to condition record for axis, save in (map_axis) :: axarr
! check that it is not a fix phase condition (istv negative)
                if(pcond%statev.lt.0) then
                   write(*,*)'Cannot set fix phase as axis'
                   goto 100
                endif
! copy the state variable to the axis record
                allocate(axarr(iax)%axcond(1))
                axarr(iax)%axcond(1)=pcond%statvar(1)
! This is probably the only reference needed for the axis condition
                axarr(iax)%seqz=pcond%seqz
                axarr(iax)%more=0
! remove axplotdef for all axis!!! 
                axplotdef=' '
             else ! a condition given as text
! check if axis variable is a condition, maybe create it if allowed
!                write(*,*)'decoding axis condition: ',text(1:20)
                stvr=>stvrvar
! this call also accept state variable functions like t_c, cp (if entered)
! UNIFINISHED: but it also accepts unknown texts ... 
                call decode_state_variable(text,stvr,ceq)
                if(gx%bmperr.ne.0) goto 990
!                write(*,*)'check if this state variable is a condition'
                pcond=>ceq%lastcondition
                i1=1; coeffs(1)=one
!                call get_condition(i1,coeffs,istv,indices,iref,unit,pcond)
                call get_condition(i1,stvr,pcond)
                if(gx%bmperr.ne.0) then
! if new conditions are allowed then maybe enter this as condition
                   write(*,*)'You must set the variable as a condition',&
                        ' before setting it as axis'
                   goto 990
                endif
                axarr(iax)%nterm=pcond%noofterms
                axarr(iax)%istv=pcond%statev
                axarr(iax)%iref=pcond%iref
                axarr(iax)%iunit=pcond%iunit
! copy the state variable record to the axis record
                if(.not.allocated(axarr(iax)%axcond)) then
                   allocate(axarr(iax)%axcond(1))
                endif
                axarr(iax)%axcond(1)=pcond%statvar(1)
! These are now redundant ...
!                jp=size(pcond%condcoeff)
!                write(*,*)'axis indices size: ',jp
!                allocate(axarr(iax)%indices(4,jp))
!                axarr(iax)%indices=pcond%indices
!                allocate(axarr(iax)%coeffs(jp))
!                axarr(iax)%coeffs=pcond%condcoeff
!
                axarr(iax)%seqz=pcond%seqz
!                write(*,*)'Condition sequential index: ',axarr(iax)%seqz
                axarr(iax)%more=0
! remove axplotdef for all axis!!! 
                axplotdef=' '
             endif
          endif removeaxis
!          dmin=axvalold(1,iax)
!          dmin=zero
          once=.TRUE.
3570      continue
          call gparrd('Minimal value:',cline,last,xxx,dmin,q1help)
          if(buperr.ne.0) goto 100
          axarr(iax)%axmin=xxx
!          axval(1,iax)=xxx
!          dmax=axvalold(2,iax)
!          dmax=one
          call gparrd('Maximal value:',cline,last,xxx,dmax,q1help)
          if(buperr.ne.0) goto 100
          if(xxx.le.axarr(iax)%axmin) then
             write(kou,*)'Maximal value must be higher than minimal'
             if(once) then
                once=.FALSE.
                goto 3570
             else
                write(kou,*)'Return to command level'
                goto 100
             endif
          endif
          axarr(iax)%axmax=xxx
!          axval(2,iax)=xxx
! default step 1/100 of difference ?? several diagram failed ...
! default step 1/40 of difference
          dinc=0.025*(axarr(iax)%axmax-axarr(iax)%axmin)
          call gparrd('Increment:',cline,last,xxx,dinc,q1help)
          if(buperr.ne.0) goto 100
          axarr(iax)%axinc=xxx
! iax can be smaller than noofaxis if an existing axis has been changed
          if(iax.gt.noofaxis) noofaxis=iax
!  write(*,3602)(axval(i,iax),i=1,3)
!3602      format(/'axlimits: ',3(1pe12.4))
!-------------------------------------------------------------
       case(15) ! set input amounts
          call set_input_amounts(cline,last,ceq)
!-------------------------
       case(16) ! SET VERBOSE
! This toggles verbose for all commands.
! it is always turned of fwhen a command is finished ...
!          write(kou,3603)'on/off',globaldata%status,GSVERBOSE
          if(btest(globaldata%status,GSSILENT)) then
! turn off VERBOSE and turn on SILENT
!             globaldata%status=ibclr(globaldata%status,GSVERBOSE)
             globaldata%status=ibclr(globaldata%status,GSSILENT)
             write(kou,3603)'off',globaldata%status
          else
! turn on VERBOSE
!             globaldata%status=ibset(globaldata%status,GSVERBOSE)
             globaldata%status=ibset(globaldata%status,GSSILENT)
             write(kou,3603)'on',globaldata%status,GSSILENT
          endif
3603      format('Silent is turned ',a,2x,z8,i5)
!          if(ocv()) then
!             write(kou,*)'Verbose mode on'
!          else
!             write(kou,*)'Verbose mode off'
!          endif
!-------------------------
! the current set of condition sill be used as start equilibrium for map/step
! Calculate the equilibrium and ask for a direction.
       case(17) ! SET AS_START_EQUILIBRIUM
          if(noofaxis.lt.2) then
             write(kou,*)'You must set two axis first'
             goto 100
          endif
          call calceq2(1,ceq)
          if(gx%bmperr.ne.0) goto 990
          call gparid('Give an axis direction: ',cline,last,ndl,2,q1help)
          if(buperr.ne.0) goto 990
          if(abs(ndl).gt.noofaxis) then
             write(kou,*)'Direction must be +/- axis number'
             goto 100
          endif
! Store a copy of equilibrium and the direction in a equential list
! starting with starteq
          eqname='_START_EQUILIBRIUM_'
          jp=len_trim(eqname)+1
          noofstarteq=noofstarteq+1
          call wriint(eqname,jp,noofstarteq)
          call copy_equilibrium(neweq,eqname,ceq)
          if(gx%bmperr.ne.0) goto 990
          neweq%multiuse=ndl
          if(associated(starteq)) then
             starteq%nexteq=neweq%eqno
          else
             starteq=>neweq
             starteq%nexteq=0
             write(*,*)'Starteq next',starteq%nexteq
          endif
          write(*,*)'A copy of current equilibrium linked as start eqilibrium'
!-------------------------
       case(18) ! SET BIT (all kinds of bits) just global implemented
!          call gparid('Toggle bit (from 0-31, -1 quits):',&
!               cline,last,ll,-1,q1help)
!         ['EQUILIBRIUM     ','GLOBAL          ','PHASE           ',&
          kom3=submenu('Set which status word?',cline,last,csetbit,nsetbit,2)
          setbit: SELECT CASE(kom3)
          CASE DEFAULT
             write(kou,*)'SET BIT subcommand error'
!................................................................
          case(1) ! equilibrium status word
!        EQNOTHREAD=0, EQNOGLOB=1, EQNOEQCAL=2,  EQINCON=3, &
!        EQFAIL=4,     EQNOACS=5,  EQGRIDTEST=6, EQGRIDCAL=7
3610         continue
!             write(kou,*)'Current equlibrium status: ',ceq%status
             write(kou,3612)ceq%status
             call gparid('Which bit? ',cline,last,ll,-1,tophlp)
             if(cline(1:1).eq.'?') then
                write(kou,3612)ceq%status
3612            format('Set/reset bits of the equilibrium status word,',/&
                     'Bit If set means',/&
                     ' 0  No threads allowed (no parallel calculation)',/&
                     ' 1  No global minimization allowed',/&
                     ' 2  No equilibrium has been calculated',/&
                     ' 3  Conditions and results not consistent',/'-'/&
                     ' 4  Last equilibrium calculation failed',/&
                     ' 5  No automatic generation of composition sets',/&
                     ' 6  Equilibrim tested by global minimizer',/&
                     ' 7  Last results are from a grid minimization'/&
                     'Current value of status word: ',z8)
                goto 3610
             endif
             if(ll.lt.0 .or. ll.gt.7) then
                write(kou,*)'No such bit, no bit changed'
             else
                call gparcd('Do you want to set the bit?',cline,last,1,&
                     ch1,'Y',q1help)
                if(ch1.eq.'Y') then
                   ceq%status=ibset(ceq%status,ll)
                   write(kou,3614)'set',ceq%status
3614               format('Bit ',a,', new equilibrium status word: ',z8)
                else
                   ceq%status=ibclr(ceq%status,ll)
                   write(kou,3614)'cleared',ceq%status
                endif
             endif
!             write(*,*)'Not implemented yet'
!................................................................
! maybe change order of questions, maybe check name exits etc ....
          case(2) ! global status word
3708         continue
! subroutine TOPHLP forces return with ? in position cline(last:last)
             write(kou,3709)globaldata%status
3709         format('Current global status word (hexadecimal): ',z8)
             call gparid('Set/reset global status bit (from 0-31, -1 quits):',&
                  cline,last,ll,-1,tophlp)
             if(cline(1:1).eq.'?') then
                write(kou,3710)
3710            format('Set/reset bits of global status word ',&
                     ' (only experts should change these) '/&
                     'Bit If set means:'/&
                     ' 0  user is a beginner'/&
                     ' 1  user is experienced'/&
                     ' 2  user is an expert'/&
                     ' 3  global minimizer will not be used'/'-'/&
                     ' 4  global minimizer must not merge comp.sets.'/&
                     ' 5  there are no data'/&
                     ' 6  there are no phases'/&
                     ' 7  comp.sets must not be created automatically'/'-'/&
                     ' 8  comp.sets must not be deleted automatically'/&
                     ' 9  data has changed since last save'/&
                     '10  verbose is on'/&
                     '11  verbose is permanently on'/'-'/&
                     '12  supress warning messages'/&
                     '13  no cleanup after an equilibrium calculation'/&
                     '14  denser grid used in grid minimizer'/&
                     '15  calculations in parallel is not allowed'/'-'/&
                     '16  no global test at node points during STEP/MAP'/&
                     '17  the components are not the elements'/&
                     '18  test if equilibrium global AFTER calculation'/&
                     '19  use old grid minimizer'/'-'/&
                     '20  do not recalculate if global test AFTER fails'/&
                     '21  use old map algorithm'/&
                     '22  no automatic startpoints for MAP'/&
                     '23-31 unused')
                goto 3708
             endif
             if(ll.lt.0 .or. ll.gt.31) then
                write(kou,*)'No bit changed'
             elseif(btest(globaldata%status,GSADV) .or. ll.le.2) then
! user must have expert bit set to change any other bit than the user type bit
                call gparcd('Do you want to set the bit?',cline,last,1,&
                     ch1,'Y',q1help)
                if(ch1.eq.'Y') then
                   globaldata%status=ibset(globaldata%status,ll)
                   write(kou,3617)ll,' set',globaldata%status
3617               format('Bit ',i2,a,', new equilibrium status word: ',z8)
                else
                   globaldata%status=ibclr(globaldata%status,ll)
                   write(kou,3617)ll,' cleared',globaldata%status
                endif
! replaced by question above
!                if(btest(globaldata%status,ll)) then
!                   globaldata%status=ibclr(globaldata%status,ll)
!                   write(*,3711)'cleared',globaldata%status
!3711               format('Bit ',a,', new value of status word: ',z8)
!                else
!                   globaldata%status=ibset(globaldata%status,ll)
!                   write(*,3711)'set',globaldata%status
!                endif
                if(.not.btest(globaldata%status,GSADV)) then
! if expert/experienced bit is cleared ensure that occational user bit is set
                   globaldata%status=ibset(globaldata%status,GSOCC)
                endif
             else
                write(kou,*)'Cannot be changed unless you have expert status'
             endif
!....................................................
          case(3) ! set bit phase ...
             write(*,*)'Please use set phase ... bit '
          end select setbit
!-------------------------
       case(19) ! set variable_coefficient, 0 to 99
          if(.not.btest(firstash%status,AHCOEF)) then
             write(kou,*)'No optimizing coefficients'
             goto 100
          endif
! zero the relative standard deviation
          firstash%coeffrsd=zero
          call gpari('Coeffient index/range: ',cline,last,i1,-1,q1help)
          if(i1.lt.0 .or. i1.ge.size(firstash%coeffstate)) then
!             write(*,*)'Dimension ',size(firstash%coeffstate)
! coefficients have indices 0 to size(firstash%coeffstate)-1
             write(kou,*)'No such coefficient'
             goto 100
          endif
! upper limit must be negative and must follow directly on same line
!          write(*,*)'pmon: ',last,': ',cline(last:last)
          if(last.lt.len(cline) .and. cline(last:last).eq.'-') then
! pick up upper range limit as a negative value, 
! the question should thus never be asked ...
             last=last-1
             call gpari('Upper index (as negative): ',cline,last,i2,-i1,q1help)
             if(i2.lt.0) then
! a negative value, its positive value must be >=i1
                i2=-i2
                if(i2.lt.i1) then
                   i2=i1
                   write(kou,*)'Illegal range, setting variable just: ',i1
                endif
!             elseif(i1.ge.size(firstash%coeffstate)) then
! coefficients have indices 0 to size(firstash%coeffstate)-1
!                i2=size(firstash%coeffstate)-1
!                write(kou,*)'Setting all coefficients fixed after ',i1
             else
! any other value ignored
                i2=i1
                write(kou,*)'Not understood, setting variable just: ',i1
             endif
          else
             i2=i1
          endif
!          write(*,*)'pmon: ',i1,i2
! possible loop if i2>i1
          j1=i1
3740      continue
!          write(*,*)'pmon: ',i1,i2,j1
          xxy=firstash%coeffvalues(j1)*firstash%coeffscale(j1)
! this coefficient is not used, igore unless i1=i2
          if(i2.gt.i1 .and. firstash%coeffstate(j1).eq.0) goto 3745
          if(firstash%coeffstate(j1).lt.10) then
             nvcoeff=nvcoeff+1
          endif
          firstash%coeffstate(j1)=10
          if(i1.eq.i2) then
! when setting a single coefficient variable ask for value
             call gparrd('Start value: ',cline,last,xxx,xxy,q1help)
             if(buperr.ne.0) goto 100
! set new value
             call change_optcoeff(firstash%coeffindex(j1),xxx)
             if(gx%bmperr.ne.0) goto 100
             firstash%coeffvalues(j1)=one
             firstash%coeffscale(j1)=xxx
             firstash%coeffstart(j1)=xxx
! current value may be scaled  ??????????
!             xxy=firstash%coeffvalues(j1)*firstash%coeffscale(j1)
!             call gparrd('Start value: ',cline,last,xxx,xxy,q1help)
!             if(buperr.ne.0) goto 100
          else
! coefficient used, set it variable with current value
             xxx=xxy
          endif
3745      if(i2.gt.j1) then
             j1=j1+1
             goto 3740
          endif
          write(kou,*)'Number of variable coefficients are ',nvcoeff
!------------------------- 
       case(20) ! set scaled_coefficient
          write(*,*)'Not implemeneted yet'
!          if(firstash%coeffstate(i1).lt.10) then
!             nvcoeff=nvcoeff+1
!          endif
! zero the relative standard deviation
!          firstash%coeffrsd=zero
!-------------------------
       case(21) ! set optimizing_conditions
!          call gparrd('DSTEP (VA05AD): ',cline,last,xxx,dstep,q1help)
!          dstep=xxx
!          call gparrd('DMAX (VA05AD): ',cline,last,xxx,dmax2,q1help)
!          dmax2=xxx
          call gparrd('LMDIF accuracy: ',cline,last,xxx,optacc,q1help)
          optacc=xxx
!-------------------------
       case(22) ! set range_experimental_equilibria
          if(allocated(firstash%eqlista)) then
             write(kou,*)'Experimental equilibria already entered'
             goto 100
          endif
          call gparid('First equilibrium number: ',cline,last,i1,2,q1help)
          j1=noeq()
          call gparid('Last equilibrium number: ',cline,last,i2,j1,q1help)
          if(i2.lt.i1) then
             write(kou,*)'No equilibria?'
             goto 100
          endif
! allocate the firstash%eqlista array and store equilibrium numbers
          j1=i2-i1+1
          firstash%firstexpeq=i1
          write(*,*)'Allocating firstash%eqlista ',j1,i1
          allocate(firstash%eqlista(j1))
          do i2=1,j1
             firstash%eqlista(i2)%p1=>eqlista(i1)
             i1=i1+1
          enddo
! close the plotdataunits!
          do i1=1,9
             if(plotdataunit(i1).gt.0) then
                write(plotdataunit(i1),22)graphopt%plotend
22              format('e'/a)
!22              format('e'/'pause mouse'/)
                close(plotdataunit(i1))
                plotdataunit(i1)=0
             endif
          enddo
!          write(*,*)'Not implemeneted yet'
!-------------------------
       case(23) ! set fixed_coefficient
          if(.not.btest(firstash%status,AHCOEF)) then
             write(kou,*)'No optimizing coefficients'
             goto 100
          endif
! zero the relative standard deviation
          firstash%coeffrsd=zero
! lower limit or range
          call gpari('Coeffient index/range: ',cline,last,i1,-1,q1help)
          if(i1.lt.0 .or. i1.ge.size(firstash%coeffstate)) then
!             write(*,*)'Dimension ',size(firstash%coeffstate)
! coefficients have indices 0 to size(firstash%coeffstate)-1
             write(kou,*)'No such coefficient'
             goto 100
          endif
! allow writing range on same line as 5-7 but also as 5 -7 on separate lines
!          write(*,*)'pmon1: ',last,': ',cline(last:last)
          frange: if(last.lt.len(cline) .and. cline(last:last).eq.'-') then
             last=last-1
! upper limit must be negative
             call gpari('Upper index limit (as negative): ',&
                  cline,last,i2,-i1,q1help)
             if(i2.lt.0) then
! a negative value, its positive value must be >=i1
                i2=-i2
                if(i2.lt.i1) then
                   i2=i1
                   write(kou,*)'Illegal range, setting fixed just: ',i1
                elseif(i2.ge.size(firstash%coeffstate)) then
                   i2=size(firstash%coeffstate)-1
                endif
             elseif(i1.ge.size(firstash%coeffstate)) then
! coefficients have indices 0 to size(firstash%coeffstate)-1
                i2=size(firstash%coeffstate)-1
                write(kou,*)'Setting all coefficients fixed after ',i1
             else
! any other value ignored
                i2=i1
                write(kou,*)'Not understood, setting fixed just: ',i1
             endif
          else
             i2=i1
          endif frange
! possible loop if i2>i1
          j1=i1
!          write(*,*)'pmon2: ',i1,j1
3720      continue
          xxy=firstash%coeffvalues(j1)*firstash%coeffscale(j1)
          if(i1.eq.i2) then
! A single coefficient, when fixing a single coefficinet ask for value
             call gparrd('Start value: ',cline,last,xxx,xxy,q1help)
             if(buperr.ne.0) goto 100
! set new value
             call change_optcoeff(firstash%coeffindex(j1),xxx)
             if(gx%bmperr.ne.0) goto 100
             firstash%coeffvalues(j1)=one
             firstash%coeffscale(j1)=xxx
             firstash%coeffstart(j1)=xxx
          else
             call get_value_of_constant_index(firstash%coeffindex(j1),xxx)
          endif
! set as fixed without changing any min/max values (first time)
!          write(*,*)'pmon3: ',xxx,firstash%coeffstate(j1)
          if(firstash%coeffstate(j1).gt.13) then
             write(kou,*)'Coefficient state wrong, set to 1'
             firstash%coeffstate(j1)=1
             nvcoeff=nvcoeff-1
          elseif(firstash%coeffstate(j1).ge.10) then
             firstash%coeffstate(j1)=max(1,firstash%coeffstate(j1)-10)
             nvcoeff=nvcoeff-1
          elseif(xxx.ne.zero) then
! mark that the coefficient is fixed and nonzero 
             firstash%coeffstate(j1)=1
          else
             firstash%coeffstate(j1)=0
          endif
          if(i2.gt.j1) then
             j1=j1+1
             goto 3720
          endif
          write(kou,3730)nvcoeff
3730      format('Number of variable coefficients are now ',i3)
!------------------------- 
       case(24) ! SET SYSTEM_VARIABLE
          call gparid('System variable index: ',cline,last,ll,0,q1help)
          if(ll.gt.0 .and. ll.le.10) then
! sysparam(2) used during STEP/MAP often to check if equilibrium is stable
             call gparid('System variable value: ',cline,last,j1,0,q1help)
             globaldata%sysparam(ll)=j1
          else
             write(*,*)'Index must be between 1 and 10'
          endif
!------------------------- 
       case(25) ! SET INITIAL_T_AND_P start values?, NOT CONDITIONS!!
          write(kou,3750)ceq%tpval
3750      format(/'NOTE: these are only local values, not conditions',&
               2(1pe12.4)/)
          call gparrd('New value of T: ',cline,last,xxx,1.0D3,q1help)
          if(buperr.ne.0) goto 100
          ceq%tpval(1)=xxx
          call gparrd('New value of P: ',cline,last,xxx,1.0D5,q1help)
          if(buperr.ne.0) goto 100
          ceq%tpval(2)=xxx
!------------------------- 
       case(26) ! unused
          continue
!------------------------- 
       case(27) ! unused
          continue
       END SELECT set
!=================================================================
! ENTER with subcommand for element, species etc
!         ['TPFUN_SYMBOL    ','ELEMENT         ','SPECIES         ',&
!         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
!         'CONSTITUTION    ','EXPERIMENT      ','QUIT            ',&
!         'EQUILIBRIUM     ','SYMBOL          ','OPTIMIZE_COEFF  ',&
!         'COPY_OF_EQUILIB ','COMMENT         ','MANY_EQUILIBRIA ',&
!         'MATERIAL        ','PLOT_DATA       ','GNUPLOT_TERMINAL',&
!         '                ','                ','                ']
    CASE(4)
! disable continue assessment optimization (not reelevant)
!       iexit=0
!       iexit(2)=1
       kom2=submenu(cbas(kom),cline,last,center,ncent,11)
       enter: SELECT CASE(kom2)
       CASE DEFAULT
          write(kou,*)'Enter subcommand error'
!---------------------------------------------------------------
! maybe change order of questions, maybe check name exits etc ....
       CASE(1) ! enter TPFUN symbol (constants, functions, tables)
          call gparc('Symbol name: ',cline,last,1,name1,' ',q1help)
          if(buperr.ne.0) goto 990
!  if(badsymname(name1)) then
          if(.not.proper_symbol_name(name1,0)) then
             write(kou,*)'Bad symbol name'
             goto 990
          endif
! check if already entered, 
          call find_tpsymbol(name1,idef,xxx)
          if(gx%bmperr.ne.0) then
! new symbol, can be function, constant or table (??)
             gx%bmperr=0
             call gparcd('Function, constant or table? ',cline,last,1,name2,&
                  'FUNCTION ',q1help)
             if(buperr.ne.0) goto 990
             call capson(name2)
             if(compare_abbrev(name2,'FUNCTION ')) then
! this call just read the function
                call enter_tpfun_interactivly(cline,last,funstring,jp)
                if(gx%bmperr.ne.0) goto 990
! here the function is stored
                lrot=0
                call enter_tpfun(name1,funstring,lrot,.FALSE.)
                if(gx%bmperr.ne.0) goto 990
             elseif(compare_abbrev(name2,'CONSTANT ')) then
! Enter a numeric constant
                call gparrd('Value: ',cline,last,xxx,zero,q1help)
                call enter_tpconstant(name1,xxx)
             elseif(compare_abbrev(name2,'TABLE ')) then
                write(kou,*)'Tables are not implemented yet'
             else
                write(kou,*)'No such type of symbol'
             endif
          else
! symbol already exist, idef=0 function, =1 constant, =2 oprimizing coefficient
             if(idef.eq.0) then
                write(kou,*)'Use AMEND to change existing TP function'
             elseif(idef.eq.2) then
                write(kou,*)'You cannot change values of optimizing ',&
                     'coefficients this way'
             else
! Values of constants can be changed here
                call gparrd('Value: ',cline,last,xxy,xxx,q1help)
                if(buperr.ne.0) goto 990
                call capson(name1)
                call enter_tpconstant(name1,xxy)
             endif
          endif
!---------------------------------------------------------------
       case(2) ! enter element
          if(.not.allowenter(1)) then
             gx%bmperr=4125
             goto 990
          endif
          call gparc('Element symbol: ',cline,last,1,elsym,' ',q1help)
          call gparcd('Element full name: ',cline,last,1,name1,elsym,q1help)
          call gparc('Element reference phase: ',cline,last,1,&
               name2,'SER ',q1help)
          call gparrd('Element mass (g/mol): ',cline,last,mass,one,q1help)
          if(buperr.ne.0) goto 990
          call gparrd('Element H298-H0: ',cline,last,h298,zero,q1help)
          if(buperr.ne.0) goto 990
          call gparrd('Element S298: ',cline,last,s298,one,q1help)
          if(buperr.ne.0) goto 990
          call enter_element(elsym,name1,name2,mass,h298,s298)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(3) ! enter species
          if(.not.allowenter(1)) then
             gx%bmperr=4125
             goto 990
          endif
          call gparc('Species symbol: ',cline,last,1,name1,' ',q1help)
! NOTE: add check species name legal!
          call gparc('Species stoichiometry: ',cline,last,1,name2,' ',q1help)
          call decode_stoik(name2,noelx,ellist,stoik)
          if(gx%bmperr.ne.0) goto 990
          call enter_species(name1,noelx,ellist,stoik)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(4) ! enter phase
          if(.not.allowenter(2)) then
             gx%bmperr=4125
             goto 990
          endif
          call enterphase(cline,last)
!---------------------------------------------------------------
       case(5) ! enter parameter only if there are phases
          if(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'You must enter some phase before'
             goto 100
          endif
! the last 0 means enter
          call enter_parameter_interactivly(cline,last,0)
! Strange things may happen when entering parameters interactively 
! This was due to an error in tpfun package ... not yet fixed ...
! Recalculate all TP functions!!  TWICE!!
          call change_optcoeff(-1,zero)
          do j1=1,notpf()
             call eval_tpfun(j1,ceq%tpval,val,ceq%eq_tpres)
             if(gx%bmperr.gt.0) goto 990
          enddo
          call change_optcoeff(-1,zero)
!          write(*,*)'pmon: A second time',notpf()
!          do j1=1,notpf()
!             call eval_tpfun(j1,ceq%tpval,val,ceq%eq_tpres)
!             if(gx%bmperr.gt.0) goto 990
!          enddo
!          call force_recalculate_tpfuns
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(6) ! enter bibliography
          call enter_bibliography_interactivly(cline,last,0,j1)
          if(gx%bmperr.ne.0) goto 990
          write(kou,*)'Bibliography number is ',j1
!---------------------------------------------------------------
       case(7) ! enter constitution
          call ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(8) ! enter experiment
! almost the same as set_condition ...
          if(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'You have no data!'
             goto 100
          endif
! enter_experiments is in models/gtp3D ...
          call enter_experiment(cline,last,ceq)
!---------------------------------------------------------------
       case(9)  ! enter QUIT
          goto 100
!---------------------------------------------------------------
       case(10) ! enter equilibrium is always allowed if there are phases
          if(.not.allowenter(3)) then
             write(kou,*)'You must have entered your system first'
             goto 100
          endif
! generate a default names line EQ_x ehere x is eqfree
          call geneqname(quest)
          call gparcd('Name: ',cline,last,1,text,quest,q1help)
          if(buperr.ne.0) goto 100
          call enter_equilibrium(text,ieq)
          if(gx%bmperr.ne.0) goto 990
! by default also select this equilibrium
          write(kou,303)ieq
303       format('Equilibrium number is ',i3)
          call gparcd('Select this equilibrium: ',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'Y' .or. ch1.eq.'y') then
             call selecteq(ieq,ceq)
          endif
!---------------------------------------------------------------
       case(11) ! enter symbol (for state variables expressions)
          call enter_svfun(cline,last,ceq)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
! enter optimizing coefficients called A00 to A99 (or whatever set as max)
       case(12)
          if(.not.allocated(firstash%coeffstate)) then
             call gparid('Number of coefficients: ',cline,last,i1,100,q1help)
             if(buperr.ne.0) goto 100
             i1=i1-1
             if(i1.lt.1) then
                write(*,*)'You must have at least 1 coefficient'
                goto 100
             elseif(i1.gt.99) then
                write(*,*)'You cannot have more than 100 coefficient'
                goto 100
             endif
             allocate(firstash%coeffvalues(0:i1))
             allocate(firstash%coeffrsd(0:i1))
             allocate(firstash%coeffscale(0:i1))
             allocate(firstash%coeffstart(0:i1))
             allocate(firstash%coeffmin(0:i1))
             allocate(firstash%coeffmax(0:i1))
             allocate(firstash%coeffindex(0:i1))
             allocate(firstash%coeffstate(0:i1))
! coeffvalues should be of the order of one
             firstash%coeffvalues=one
             firstash%coeffrsd=zero
             firstash%coeffscale=zero
             firstash%coeffstart=zero
             firstash%coeffmin=zero
             firstash%coeffmax=zero
             firstash%coeffindex=0
             firstash%coeffstate=0
! create the corresponding TP constants for coeffvalues
             call enter_optvars(j1)
             call makeoptvname(name1,i1)
             write(kou,556)name1(1:3),i1
556          format(/'Coefficients entered with symbols A00 to ',a/&
                  'Note that indices are from 0 to ',i2)
             do i2=0,i1
                firstash%coeffindex(i2)=j1+i2
             enddo
             firstash%status=ibset(firstash%status,AHCOEF)
          else
             write(kou,553)size(firstash%coeffstate)
553          format('You have already ',i3,' optimizing coefficients entered')
          endif
          call gparid('Size of workspace: ',cline,last,lwam,2500,q1help)
          if(lwam.gt.2000) lwam=2000
          if(allocated(wam)) then
             deallocate(wam)
             deallocate(iwam)
          endif
!          write(*,551)firstash%status
!551       format('Assessment status word: ',z8)
!---------------------------------------------------------------
! enter copy_of equilibrium (for test!)
       case(13)
! Check if there is any phases, otherwise not allowed
          if(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'Not allowed unless you have data!'
             goto 100
          endif
          call gparc('Name of new equilibrium: ',cline,last,1,text,' ',q1help)
          if(buperr.ne.0) goto 100
          if(text(1:1).eq.' ') then
             write(*,*)'You must specify a unique name'
             goto 100
          endif
          call copy_equilibrium(neweq,text,ceq)
!          write(*,*)'Back from copy equilibrium'
          if(gx%bmperr.ne.0) goto 990
          write(kou,*)'New equilibrium no: ',neweq%eqno
!---------------------------------------------------------------
! enter COMMENT for current equilibrium
       case(14)
          write(*,*)'Current equilibrium name: ',ceq%eqname
          call gparc('One line text: ',cline,last,5,text,' ',q1help)
          ceq%comment=text
!---------------------------------------------------------------
! enter MANY_EQUILIBRIA
! The plotdataunit array should be zero at first call, then the unit is opened
! (if there are any plot_data commands).  It will remain open until
! a set range command is given
       case(15)
          call enter_many_equil(cline,last,plotdataunit)
!---------------------------------------------------------------
! enter MATERIAL
! ask for database, then major element, mass/mole fraction of elements
! read the database; jump possibly to SCHEIL/STEP calculation 
! or simply ask for T and calculate equilibrium; 
       case(16)
          call enter_material(cline,last,nv,xknown,ceq)
          if(gx%bmperr.ne.0) goto 990
          xxy=firsteq%tpval(1)
          call gparrd('Temperature ',cline,last,xxx,xxy,q1help)
! set T and P
          cline='P=1E5 T='
          i1=len_trim(cline)+1
          call wrinum(cline,i1,10,0,xxx)
          i1=0
          call set_condition(cline,i1,ceq)
! calculate the equilibrium
          call calceq2(1,ceq)
          if(gx%bmperr.ne.0) then
             ceq%status=ibset(ceq%status,EQFAIL)
             goto 990
          endif
!---------------------------------------------------------------
! enter PLOT DATA
! the file ocmanyi.plt with unit plotdataunit(i) must already be open!
! it is opened in the enter_many_equilibria if there is a plot_data command
       case(17)
          call gparid('Dataset number:',cline,last,i1,1,q1help)
! here only the normal plotdata units 1 to 9 are legal
          if(i1.gt.0 .and. i1.lt.10) then
             if(plotdataunit(i1).lt.10) then
                write(kou,*)'No plotdata file for this dataset'
                goto 100
             endif
             call gparrd('X coordinate:',cline,last,xxx,zero,q1help)
             call gparrd('Y coordinate:',cline,last,xxy,one,q1help)
             call gparid('Symbol:',cline,last,i2,1,q1help)
             write(plotdataunit(i1),171)i1,xxx,xxy,i2
171          format(i3,2(1pe14.6),i5,' have a nice day')
          else
             write(kou,*)'No plotdata file for dataset ',i1
          endif
!---------------------------------------------------------------
! ENTER GNUPLOT_TERMINAL
       case(18)
          write(kou,172)graphopt%gnutermax
172       format('GNUPLOT terminals are:',i2/&
               4x,'Name',5x,'> command',6x,'GNUPLOT options')
          write(kou,173)(i2,graphopt%gnutermid(i2),&
               trim(graphopt%gnuterminal(i2)),i2=1,graphopt%gnutermax)
173       format(i2,2x,a,' > set terminal ',a)
          write(kou,174)
174       format('Change (exact match required) or enter a new GNUPLOT termial')
          call gparc('Terminal id (8 chars):',cline,last,1,text,' ',q1help)
          call capson(text)
          if(text(1:1).eq.' ') goto 100
          do i1=1,graphopt%gnutermax
             if(text(1:8).eq.graphopt%gnutermid(i1)) then
                string=graphopt%gnuterminal(i1)
                write(*,*)'Modifying terminal ',graphopt%gnutermid(i1)
                goto 176
             endif
          enddo
! gnutermid not found, a new terminal
          call gparcd('You want to enter a new terminal "'//trim(text)//'"?',&
               cline,last,1,ch1,'Y',q1help)
          if(ch1.ne.'Y') then
             write(*,*)'Please try again'; goto 100
          endif
          if(graphopt%gnutermax.ge.8) then
             write(kou,*)'There can max be 8 terminals'
             goto 100
          endif
          i1=graphopt%gnutermax+1
          graphopt%gnutermax=i1
          string=' '
! enter a new set terminal id and definition
176       continue
          graphopt%gnutermid(i1)=text(1:8)
          call gparc('Text after set terminal (see GNUPLOT manual):',&
               cline,last,5,text,string,q1help)
          graphopt%gnuterminal(i1)=text
          if(i1.ne.1) then
! SCREEN has no file extention
             call gparc('File extention:',cline,last,1,text,' ',q1help)
             graphopt%filext(i1)=text(1:4)
          endif
          write(*,179)i1,graphopt%gnutermid(i1),trim(graphopt%gnuterminal(i1)),&
               trim(graphopt%filext(i1))
179       format('New terminal definition for plot '/&
               i2,2x,a,'set terminal ',a/4x,'with file extention: ',a)
!----------------------------------------------------------------
! enter not used
       case(19)
          write(*,*)'Not implemeneted yet'
! this is at amend components
!          i2=1
!          line=' '
!          do i1=1,noel()
!             call get_component_name(i1,line(i2:),ceq)
!             i2=len_trim(line)+2
!          enddo
!          call gparcd('Give all new components: ',cline,last,&
!               5,option,line,q1help)
!          call enter_components(option,ceq)
!          if(gx%bmperr.ne.0) goto 990
!----------------------------------------------------------------
! enter unused
       case(20)
!----------------------------------------------------------------
! enter unused
       case(21)
       END SELECT enter
!=================================================================
! exit
    CASE(5)
       call gparcd('Are you sure?',cline,last,1,ch1,'N',q1help)
       if(ch1.eq.'y' .or. ch1.eq.'Y') then
          if(logfil.gt.0) then
             write(logfil,*)'set interactive'
          endif
          call openlogfile(' ',' ',-1)
          stop 'Ha en bra dag'
       endif
!=================================================================
! list with subcommands
!        ['DATA            ','SHORT           ','PHASE           ',&
!         'STATE_VARIABLES ','BIBLIOGRAPHY    ','MODEL_PARAM_ID  ',&
!         'AXIS            ','TPFUN_SYMBOLS   ','QUIT            ',&
!         'PARAMETER       ','EQUILIBRIA      ','RESULTS         ',&
!         'CONDITIONS      ','SYMBOLS         ','LINE_EQUILIBRIA ',&
!         'OPTIMIZATION    ','MODEL_PARAM_VAL ','ERROR_MESSAGE   ',&
!         ,ACTIVE_EQUILIBR ','ELEMENTS        ','                ']
! SHOW is the same as LIST STATE_VARIABLES including also CALC SYMBOL !!
! SHOW is cammand 25
    CASE(6,25) ! LIST andSHOW
       if(kom.eq.25) then
          kom2=4
       else
          kom2=submenu(cbas(kom),cline,last,clist,nclist,12)
          if(kom2.le.0) goto 100
       endif
       lut=optionsset%lut
!       write(*,*)'PMON: show xliqni should come here ... YES ',kom,kom2
       list: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'LIST FORMAT subcommand error'
          goto 100
!-----------------------------------------------------------
       case(1) ! list data, not dependent on equilibrium!!
! NOTE output file for SCREEN can be set by /output=
! LIST DATA SCREEN/TDB/MACRO/LaTeX
! it is also possible to give SAVE TDB 
          kom3=submenu('Output format?',cline,last,llform,nlform,1)
          if(kom.gt.0) then
             call list_many_formats(cline,last,kom3,kou)
             if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
                write(kou,*)bmperrmess(gx%bmperr)
             elseif(gx%bmperr.ne.0) then
                write(kou,*)'Error code ',gx%bmperr
             endif
          else
             write(kou,*)'Unknown format'
          endif
!-----------------------------------------------------------
       case(2,20) ! list short with status bits, also LIST ELEMENT (kom2=20)
          if(kom2.eq.20) then
             ch1='C'
          else
             call gparcd('Option (A/C/M/P)',cline,last,1,ch1,chshort,q1help)
             call capson(ch1)
          endif
          write(lut,6022)ceq%eqname,globaldata%rgasuser,&
               globaldata%pnorm,globaldata%status,ceq%status
6022      format('Equilibrium name',9x,'Gas constant Pressure norm',&
               5x,'Status Global   Equilib'/&
               1x,a,1pe12.4,2x,1pe12.4,10x,z8,2x,z8)
!....................................................................
! options are A=all phases; P=some phases; C=components; M=phase models
          if(ch1.eq.'A') then
! A all
             chshort='A'
             call list_all_elements(lut)
             call list_all_species(lut)
             call list_all_phases(lut,ceq)
!....................................................................
          elseif(ch1.eq.'P') then
! just the phases
! P phases sorted: stable/ unstable in driving force order/ dormant the same
             chshort='P'
             call list_sorted_phases(lut,ceq)
             if(btest(ceq%status,EQFAIL)) write(lut,6305)'above'
          elseif(ch1.eq.'C') then
!....................................................................
! global values and the chemical potentials
             chshort='C'
             write(kou,*)
             call list_global_results(lut,ceq)
!             write(lut,6303)'Some component data ....................'
             write(lut,6303)'Some data for components ...............'
             j1=1
             if(listresopt.ge.4 .and. listresopt.le.7) then
                j1=2
             endif
             call list_components_result(lut,j1,ceq)
!....................................................................
          elseif(ch1.eq.'M') then
! list models for all phases
             do iph=1,noph()
                call list_phase_model(iph,1,lut,' ',ceq)
             enddo
!....................................................................
          else
             write(kou,*)'Only option A, C, M and P implemented'
          endif
!-----------------------------------------------------------
       case(3) ! list phase subcommands
          call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
          if(buperr.ne.0) goto 990
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
          kom3=submenu('List what for phase?',cline,last,clph,nclph,2)
          listphase: SELECT CASE(kom3)
!...............................................................
          CASE DEFAULT
             write(kou,*)'list phase subcommand error'
!...............................................................
          CASE(1) ! list phase data
             call list_phase_data(iph,' ',lut)
!...............................................................
! list phase constitution
          case(2) ! list phase constitution
!             idef=110
!             call gparid('Output mode: ',cline,last,mode,idef,q1help)
!             if(buperr.ne.0) goto 990
!  call list_phase_results(iph,ics,mode,kou,firsteq)
             write(lut,6051)ceq%eqno,ceq%eqname
6051         format(/'Output for equilibrium: ',i3,', ',a,5x,a4,'.',a2,'.',a2)
             mode=110
             once=.TRUE.
             call list_phase_results(iph,ics,mode,lut,once,ceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Last equilibrium calculation failed'
                goto 990
             endif
!...............................................................
          case(3) ! list phase model (including disordered fractions)
             write(kou,6070)'For ',ceq%eqno,ceq%eqname
6070      format(a,'equilibrium: ',i3,', ',a)
             call list_phase_model(iph,ics,lut,' ',ceq)
          END SELECT listphase
!------------------------------
! THIS IS ALSO THE SHOW command and list model-parameter-value case(17) of LIST
       case(4,17)  ! list state_variable or model_parameter_value, or SHOW
!6099      continue
          if(btest(ceq%status,EQNOEQCAL) .or. btest(ceq%status,EQFAIL)) then
             write(lut,6101)
6101         format(' *** Warning,',&
                'equilibrium not calculated, values are probably wrong')
          elseif(btest(ceq%status,EQINCON)) then
             write(lut,6102)
6102         format(' *** Warning, values can be inconsistent with',&
                ' current conditions')
          endif
          once=.TRUE.
! LOOP here for list state_variables or model_parameter_values or SHOW
6105      continue
!          write(*,*)'At label 6105: ',last,' "',trim(cline),'"',kom,kom2
! NOTE: 4th argument is 5 because otherwise a "," will terminate reading cline
! and state variables like x(fcc,cr) will not work.
          if(kom.eq.25) then
! the command is SHOW             
!             write(*,*)'PMON: show xliqni should come here ... YES '
             call gparc('property: ',cline,last,5,line,' ',q1help)
          else
! the command is LIST STATE_VARIABLES
             if(kom2.eq.4) then
                call gparc('State variable: ',cline,last,5,line,' ',q1help)
             else
! the command is LIST MODEL_PARAMETER_VALUE                
                if(once) then
                   write(kou,*)'Remember always to specify the phase!'
                   once=.FALSE.
                endif
                call gparc('Parameter ident: ',cline,last,5,line,' ',q1help)
             endif
          endif
! if line empty return to command level
          j1=1
          if(eolch(line,j1)) goto 100
          j1=index(line,',')
          if(j1.gt.0) then
! check if there is a , before a ( as that is not allowed.  There are
! state variables like x(fcc,cr) ... (this is not a strong test ...)
             ll=index(line,'(')
             if(ll.le.0 .or. ll.gt.j1) then
                write(*,*)'Please use a space as separator',&
                     ' except within ( ) as in x(liq,cr) !'
                goto 100
             endif
          endif
! model is just used to return texts
          model=' '
! we should extract the text from last up to first space and save rest in cline
          j1=index(line,' ')
          name1=line(1:j1)
          call capson(name1)
! note gparc etc increment last before looking for answer, keep space in cline
          cline=line(j1:)
          last=1
          if(index(name1,'*').gt.0) then
! generate many values
! i1 values are resturned in yarr with dimension maxconst. 
! longstring are the state variable symbols for the values ...
             call get_many_svar(name1,yarr,maxconst,i1,longstring,ceq)
             if(gx%bmperr.eq.0) then
! not a nice output, we should include the state variables FIX!!
                write(lut,6106)i1,longstring(1:len_trim(longstring))
6106            Format('Listing of ',i3,' state variables:'/a)
                write(lut,6107)(yarr(i2),i2=1,i1)
6107            format('Values: ',5(1pe14.6)/(8x,5(1pe14.6)))
                if(index(name1,'*,').gt.0) write(*,6121)trim(name1)
6121            format(' *** Note that for unstable phases ',a,&
                     ' is not shown or listed as zero')
             endif
          else
! the value of a state variable, symbol? or model parameter variable is returned
! STRANGE the symbol xliqni is accepted in get_state_var_value ???
!             write(*,*)'pmon show: ',name1
             call get_state_var_value(name1,xxx,model,ceq)
!             write(*,*)'PMON: show xliqni should come here 6 ... ',gx%bmperr
             if(gx%bmperr.eq.0) then
                write(lut,6108)model(1:len_trim(model)),xxx
6108            format(1x,a,'=',1PE15.7)
             else
                gx%bmperr=0
!                write(*,*)'PMON: show xliqni should come here ... NO!!!'
! If error then try to calculate a symbol ...
! below copied from calculate symbol, first calculate all symbols ignore errors
! calculate all symbols ignoring errors (note dot derivatives not calculated)
                call meq_evaluate_all_svfun(-1,ceq)
                if(gx%bmperr.ne.0) gx%bmperr=0
                call capson(line)
!                call find_svfun(name1,istv,ceq)
!                write(*,*)'PMON: calling find_svfun again ...'
                call find_svfun(name1,istv)
                if(gx%bmperr.ne.0) goto 990
                mode=1
                actual_arg=' '
                xxx=meq_evaluate_svfun(istv,actual_arg,mode,ceq)
!                write(*,*)'pmon error: calling meq_evaluate_svfun',istv,xxx
                if(gx%bmperr.ne.0) goto 990
                write(kou,2047)trim(name1),xxx
! this format statement elsewhere
!2047            format(a,'= ',1pe16.8)
             endif
          endif
          if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
             write(lut,*)bmperrmess(gx%bmperr)
          elseif(gx%bmperr.ne.0) then
             write(lut,*)'Error code ',gx%bmperr
          endif
          gx%bmperr=0
! try to pick up more properties etc from cline if not empty
          if(.not.eolch(cline,last)) then
! there are more symbols, state variables or model_parameters in cline
             last=last-1
             goto 6105
          elseif(kom.ne.25) then
! for list state_variables and list model_parameter_value ask for more
             goto 6105
          endif
!-----------------------------------------------------------
       case(5) ! list data bibliography
          call gparcd('Bibliographic id:',cline,last,1,name1,'ALL',q1help)
          if(name1.eq.'ALL ') name1=' '
          call list_bibliography(name1,lut)
!-----------------------------------------------------------
       case(6) ! list model_parameter_identifiers
          call list_defined_properties(lut)
!-----------------------------------------------------------
       case(7) ! list axis
          if(noofaxis.le.0) then
             write(kou,*)'No axis set'
             goto 100
          endif
          write(lut,6131)
6131      format(4x,'Axis variable',12x,'Min',9x,'Max',9x,'Max increment')
!6131      format(4x,'Axis variable',12x,'Start',7x,'Final',7x,'Increment')
          do iax=1,noofaxis
             jp=1
             call get_one_condition(jp,text,axarr(iax)%seqz,ceq)
             if(gx%bmperr.ne.0) then
                write(kou,*)'Condition sequential index: ',iax,axarr(iax)%seqz
                goto 990
             endif
! we just want the expression, remove the value including the = sign
             jp=index(text,'=')
             text(jp:)=' '
!             write(kou,6132)iax,axvar(iax),(axval(j1,iax),j1=1,3)
             write(lut,6132)iax,text(1:24),&
                  axarr(iax)%axmin,axarr(iax)%axmax,axarr(iax)%axinc
6132         format(i2,2x,a,3(1pe12.4))
          enddo
!-----------------------------------------------------------
       case(8) ! list tpfun symbol
          call gparcd('name: ',cline,last,5,name1,'*',q1help)
          lrot=0
          iel=index(name1,'*')             
          if(iel.gt.1) name1(iel:)=' '
          if(name1(1:1).ne.'*') then
6140         continue
             call find_tpfun_by_name(name1,lrot)
!             write(*,*)'cui: ',lrot,iel,gx%bmperr
             if(gx%bmperr.ne.0) then
                if(iel.eq.0) goto 990
                gx%bmperr=0
             else
                longstring=' '
                write(longstring,6142)lrot
6142            format(i5)
                jp=len_trim(longstring)+2
                call list_tpfun(lrot,0,longstring(jp:))
                call wrice2(lut,0,12,78,1,longstring)
                if(iel.gt.1) goto 6140
             endif
          else
             call list_all_funs(lut)
          endif
!------------------------------------------------------------
       case(9) ! list quit
!------------------------------------------------------------
       case(10) ! list parameter for a phase (just one). Last 1 means list
          call enter_parameter_interactivly(cline,last,1)
!-----------------------------------------------------------
       case(11,19) ! list equilibria and list ACTIVE_EQUILIBRIA (not result)
! if 19 then skip equilibria with zero weight
          nv=noeq()
! skip if there is just one equilibrium kom=6=LIST; kom2=19=ACTIVE-EQUIL
!          write(*,*)'PMON: ',kom,kom2,nv
          if(kom2.eq.19 .and. nv.eq.1) goto 100
          write(*,6212)
6212      format('Number  Name',25x,'T   Weight Comment')
          jp=0
          do iel=1,nv
             if(associated(ceq,eqlista(iel))) then
                name1='**'
             else
                name1=' '
             endif
!             write(*,*)'PMON: ',kom2,iel,eqlista(iel)%weight,jp
!             j1=len_trim(eqlista(iel)%comment)
             if(eqlista(iel)%weight.gt.zero) then
! always list equilibria with weight>0
                write(lut,6203)iel,name1(1:2),eqlista(iel)%eqname,&
                     eqlista(iel)%tpval(1),eqlista(iel)%weight,&
                     eqlista(iel)%comment(1:28)
!                     eqlista(iel)%comment(1:32)
!6203            format(i4,1x,a2,1x,a,' T=',F8.2,', weight=',F5.2,', ',a)
6203            format(i4,1x,a2,1x,a,1x,F8.2,1x,F5.2,1x,a)
             elseif(iel.eq.1 .or. kom2.eq.11) then
! for kom2=11 list all equilibria without including weight
! NOTE all equilibria outside "range" (default and step/map) has weight=-1.0
                write(lut,6202)iel,name1(1:2),eqlista(iel)%eqname,&
                     eqlista(iel)%tpval(1),eqlista(iel)%comment(1:28)
!                     eqlista(iel)%tpval(1),eqlista(iel)%comment(1:32)
!6202            format(i4,1x,a2,1x,a,' T=',F8.2,', ',a)
6202            format(i4,1x,a2,1x,a,1x,F8.2,7x,a)
             elseif(eqlista(iel)%weight.eq.zero) then
                jp=jp+1
             endif
!             if(j1.gt.1) then
!                write(lut,6204)eqlista(iel)%comment(1:j1)
!6204            format(12x,a)
!             endif
          enddo
          if(kom2.eq.19 .and. jp.gt.0) &
               write(*,'(/"Number of equilibria with zero weight: ",i4)')jp
!------------------------------
       case(12) ! list results
! skip if no calculation made
          if(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'No results as no data'
             goto 100
          elseif(btest(ceq%status,EQGRIDCAL)) then
             write(kou,*)' *** Last calculation was not a full equilibrium'
          endif
          call gparid('Results output mode: ',cline,last,&
               listresopt,lrodef,q1help)
          if(buperr.ne.0) then
             write(kou,*)'No such mode, using default'
             buperr=0
             listresopt=lrodef
          endif
! CCI extending the number of listing options
!          if(listresopt.gt.0 .and. listresopt.le.9) then
          if(listresopt.gt.0 .and. listresopt.le.11) then
             lrodef=listresopt
          endif
! CCI end          
          call date_and_time(optres,name1)
          write(lut,6051)ceq%eqno,ceq%eqname,optres(1:4),optres(5:6),optres(7:8)
! write comment line if any
          if(len_trim(ceq%comment).gt.0) then
             write(lut,6308)trim(ceq%comment)
6308         format(3x,a)
          endif
          if(btest(ceq%status,EQFAIL)) then
             write(lut,6305)'below'
6305         format(/' *** The results ',a,&
                  ' are not a valid equilibrium as last calculation failed'/)
!  elseif(btest(globaldata%status,GSNOEQCAL)) then
          elseif(btest(ceq%status,EQNOEQCAL)) then
             write(lut,6307)'below'
6307         format(/' *** The results listed ',a,'does not represent',&
                  ' a calculated equilibrium'/)
          elseif(btest(ceq%status,EQINCON)) then
             write(lut,6306)'below'
6306         format(/' *** The results listed ',a,' may be inconsistent',&
                  ' with the current conditions'/)
          endif
          write(lut,6302)'Conditions .............................'
6302      format(a,20('.'),':')
6303      format(/a,20('.'),':')
          call list_conditions(lut,ceq)
          write(lut,6303)'Some global data, reference state SER ..'
          call list_global_results(lut,ceq)
          if(btest(ceq%status,EQNOEQCAL)) then
             write(*,6277)ceq%status
6277         format(' *** No results as no equilibrium calculated! ',z8)
             goto 6363
          endif
!          write(lut,6303)'Some component data ....................'
          write(lut,6303)'Some data for components ...............'
          j1=1
          if(listresopt.ge.4 .and. listresopt.le.7) then
! j1=2 means mass fractions
             j1=2
          endif
          call list_components_result(lut,j1,ceq)
! Phase output starts with newline
!         write(lut,6304,advance='no')'Some Phase data ........................'
          write(lut,6304,advance='no')'Some data for phases ...................'
6304      format(/a,20('.'),':')
          if(listresopt.le.1) then
! 1: stable phases with mole fractions in value order 
             mode=1000
          elseif(listresopt.eq.2) then
! 2: stable phases with mole fractions and constitution in value order
             mode=1010
          elseif(listresopt.eq.3) then
! 3: stable phases with mole fractions and constitution in alphabetical order
             mode=1110
          elseif(listresopt.eq.4) then
! 4: stable phases with mass fractions in value order
             mode=1001
          elseif(listresopt.eq.5) then
! 5: stable phases with mass fractions in alphabetical order
             mode=1101
          elseif(listresopt.eq.6) then
! 6: stable phases with mass fractions and constitution in value order
             mode=1011
          elseif(listresopt.eq.7) then
! 7: all phases with mass fractions in value order
             mode=1
          elseif(listresopt.eq.8) then
! 9: all phases with mole fractions in alphabetical order
             mode=110
          elseif(listresopt.eq.9) then
! 9: all phases with mole fractions an constitutions in value order
             mode=10
          elseif(listresopt.eq.10) then
! CCI
! 10: stable phases with constituent fractions time FU of hase in value order
! SOLGASMIX type output
             mode=10000
          elseif(listresopt.eq.11) then
! 11: stable phases with constituent fractions time FU of hase in value order
             mode=10010
! CCI end             
          else
! all phase with with mole fractions
             mode=0
          endif
          ics=1
          once=.TRUE.
          do iph=1,noph()
             ics=0
6310         continue
             ics=ics+1
! moved to gtp3C
!             if(listresopt.ge.4 .and. listresopt.le.7) then
! use phase amount in mass
!                write(lut,6308)'Mass      '
!6308            format('Name                Status ',a,' Volume',&
!                 '    Form.U    At/FU     DGM    X/W:')
!                     '    Form.U    At/FU     DGM   Frac:')
!             else
! use phase amount in mole
!                write(lut,6308)'Moles     '
!             endif
             call list_phase_results(iph,ics,mode,lut,once,ceq)
             if(gx%bmperr.ne.0) then
! if error take next phase
                gx%bmperr=0
             else
! else take next composition set
                goto 6310
             endif
          enddo
! list experiments if any
6363      continue
          if(associated(ceq%lastexperiment)) then
             write(lut,491)ceq%weight
!491          format(/'Weight ',F6.2)
491          format('Weight ',F6.2)
! list all experiments ........................................
             call meq_list_experiments(lut,ceq)
             write(lut,*)
!          else
!             write(*,*)'No experiments found'
          endif
          if(btest(ceq%status,EQNOEQCAL)) goto 100
! list if anyting should be calculated or listed separately (not their values)
          if(allocated(ceq%eqextra)) then
             write(lut,492)ceq%eqextra(1)(1:len_trim(ceq%eqextra(1))),&
                  ceq%eqextra(2)(1:len_trim(ceq%eqextra(2)))
492          format('Calculate: ',a/'List: ',a)
!          else
!             write(*,*)'No extra lines'
          endif
! make sure phases with positive DGM listed
          call list_phases_with_positive_dgm(mode,lut,ceq)
          if(btest(ceq%status,EQFAIL)) then
             write(lut,6305)'above'
          elseif(btest(ceq%status,EQNOEQCAL)) then
             write(lut,6307)'above'
          elseif(btest(ceq%status,EQINCON)) then
             write(lut,6306)'above'
          endif
!------------------------------
       case(13) ! list conditions
          write(kou,6070)'Conditions for ',ceq%eqno,ceq%eqname
          call list_conditions(lut,ceq)
!------------------------------
       case(14) ! list symbols (state variable functions, not TP funs)
          call list_all_svfun(lut,ceq)
!------------------------------
! list lines, output of calculated and stored equilibria
       case(15)
! temporary listing of all stored equilibria as test
! IDEA: Add question for symbols and state variables to be listed!!
! Add a heading to make spece for more dara
! ceq #; Next;      T;  axis value; 0-n user symbols;           
!  9999  9999  20000.00 +1.2345E+00 1.2345E+00 1.2345E+00 1.2345E+00 1.2345E+00
          call list_stored_equilibria(lut,axarr,maptop)
!------------------------------
! list optimization, several suboptions
!    character (len=16), dimension(noptopt) :: optopt=&
!        ['SHORT           ','LONG            ','COEFFICIENTS    ',&
!         'GRAPHICS        ','DEBUG           ','MACRO           ',&
!         'EXPERIMENTS     ','CORRELATION_MTRX','                ']
       case(16)
          if(.not.allocated(firstash%coeffstate)) then
             write(kou,*)'No listing as no optimizing parameters'
             goto 100
          endif
          call date_and_time(optres,name1)
          kom2=submenu('List ',cline,last,optopt,noptopt,1)
! allow output file
          lut=optionsset%lut
          write(lut,600)optres(1:4),optres(5:6),optres(7:8),&
               name1(1:2),name1(3:4),err0(3)
600       format(/'Optimization results at ',a4,'.',a2,'.',a2,&
               ':',a2,'h',a2,', normalized sum of error: ',1pe12.4)
          listopt: SELECT CASE(kom2)
!..........................................................
             case DEFAULT
                write(kou,*)'No such option'
!...........................................................
! list optimization short
             case(1) ! short
!                if(updatemexp) then
!                   write(*,*)'You must OPTIMIZE first'
!                   goto 100
!                endif
                if(allocated(errs)) then
                   call listoptshort(lut,mexp,nvcoeff,errs)
                else
                   write(kou,*)'No current optimization'
                endif
                call listoptcoeff(mexp,err0,.FALSE.,lut)
!...........................................................
! list optimization long
             case(2) ! long
                write(*,*)'Not implemented yet'
!...........................................................
! list optimization coefficients
             case(3) ! coefficient values
                if(mexp.eq.mexpdone .and. nvcoeff.eq.nvcoeffdone) then
                   call listoptcoeff(mexp,err0,.TRUE.,lut)
                else
                   call listoptcoeff(mexp,err0,.FALSE.,lut)
                endif
!...........................................................
! list optimization graphics, plot calculated vs experiment values
             case(4) ! graphics
                write(*,*)'Not implemented yet'
!...........................................................
! list optimization debug ??
             case(5) ! debug
                if(nvcoeff.ne.nvcoeffdone .or. mexp.ne.mexpdone) then
                   write(*,*)'No optimization done with current set of ',&
                        'coefficients or experiments'
                   goto 100
                elseif(.not.allocated(fjac)) then
                   write(*,*)'No optimization done'
                   goto 100
                endif                
                write(*,*)'Listing the Jacobian: ',nvcoeff,mexp
!                iflag=2
!                call fdjac2(mexp,nvcoeff,coefs,errs,fjac,mexp,iflag,zero,wam)
!                write(*,*)'fjac: ',nvcoeff,mexp,iflag
                do i2=1,mexp
                   write(*,563)(fjac(i2,ll),ll=1,nvcoeff)
                enddo
                if(allocated(cov1)) then
                   write(*,*)'The covariance matrix Jac^T * Jac: '
                   do i2=1,nvcoeff
                      write(*,563)(cov1(i2,ll),ll=1,nvcoeff)
                   enddo
                endif
!...........................................................
! list optimization macro: create macro file with all experiments
             case(6) ! MACRO include experiments
                write(*,*)'Not implemented yet'
!...........................................................
! list optimization experiments
             case(7) ! experiments with weight>0
                if(allocated(errs)) then
                   call listoptshort(lut,mexp,nvcoeff,errs)
                else
                   write(kou,*)'No current optimization'
                endif
!...........................................................
! list optimization correlation matrix
             case(8) ! unused
                if(nvcoeff.eq.nvcoeffdone .and. allocated(cormat)) then
                   write(kou,*)'Correlation matrix is (symmetric):'
                   do i2=1,nvcoeff
                      write(kou,563)(cormat(i2,j2),j2=1,i2)
                   enddo
                   write(kou,*)'Covariance matrix is (symmetric): '
                   do i2=1,nvcoeff
                      write(kou,563)(cov1(i2,j2),j2=1,i2)
                   enddo
                else
                   write(*,*)'No correlation matrix calculated'
                endif
!...........................................................
! list optimization RSD (according to OC and TC)
             case(9) ! unused
                write(kou,3998)
3998            format(/'Relative Standard Deviation (RSD) values according',&
                     ' to OC and TC'/'Variable  OC          TC')
                i2=0
                do i1=0,size(firstash%coeffstate)-1
                   if(firstash%coeffstate(i1).ge.10) then
                      i2=i2+1
                      write(*,'(i7,2(1pe12.4))')i2,&
                           sqrt(abs(cov1(i2,i2))),&
                           sqrt(abs(tccovar(i2,i2)))
                   endif
                enddo
                write(kou,3999)sqrt(err0(3))
3999            format('The difference is the square root of the normalized',&
                     ' sum or errors: ',1pe12.4)
             end SELECT listopt
!------------------------------
! list model_parameter_values, part of case(4)
!       case(17)
!          write(*,*)'Not implemented yet'
!------------------------------
! list error message
       case(18)
          i2=4204
          call gparid('Error code: ',cline,last,i1,i2,q1help)
          if(i1.ge.4000 .and. i1.le.nooferm) then
             write(kou,4999)i1,bmperrmess(i1)
4999         format('The error code ',i4', means: '/a)
          else
             write(kou,*)'Not a standard OC error message'
          endif
!------------------------------
! list ?? nonzero_equilibria/active-equil merged with list equilibra, case 11
!       case(19)
!          write(*,*)'Not implemented yet'
!------------------------------
! list elements, same as LIST SHORT C (components) case 2
!       case(20)
!          write(*,*)'Not implemented yet'
!------------------------------
! list ??
       case(21)
          write(*,*)'Not implemented yet'
       end SELECT list
!=================================================================
! quit
    case(7)
       if(cline(1:1).eq.'q') then
          call gparcd('Are you sure?',cline,last,1,ch1,'N',q1help)
       else
! upper case Q will quit without question
          ch1='y'
       endif
       if(ch1.eq.'y' .or. ch1.eq.'Y') then
          if(logfil.gt.0) then
             write(logfil,*)'set interactive'
          endif
          call openlogfile(' ',' ',-1)
          stop 'Have a nice day'
       endif
!=================================================================
! READ subcommand
!        ['UNFORMATTED     ','TDB             ','QUIT            ',&
!         'DIRECT          ','                ','                ']
    case(8)
! disable continue optimization
!       iexit=0
!       iexit(2)=1
       if(noel().ne.0) then
          write(kou,*)'You already have data, read destroys your current data'
          write(kou,*)'You must give a NEW Y command to remove data first'
!          call gparcd('Do you want to continue? ',cline,last,1,ch1,'N',q1help)
!          if(.not.(ch1.eq.'y' .or. ch1.eq.'Y')) then
!             write(kou,*)'Command ignored, to read answer yes' 
          goto 100
!       else
! all records must be removed and init_gtp is called.  This is fragile ...
!             call new_gtp
!             if(gx%bmperr.ne.0) goto 990
!             write(kou,*)'All previous data deleted'
!          endif
       endif
       kom2=submenu(cbas(kom),cline,last,cread,ncread,2)
       read: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'Read subcommand error'
!-----------------------------------------------------------
       case(1) ! read unformatted file created by SAVE
          if(ocufile(1:1).ne.' ') then
             text=ocufile
             call gparcd('File name: ',cline,last,1,ocufile,text,q1help)
          else
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
             call gparfile('File name: ',cline,last,1,ocufile,' ',2,q1help)
!             call gparc('File name: ',cline,last,1,ocufile,' ',q1help)
          endif
          call gtpread(ocufile,text)
          if(gx%bmperr.ne.0) then
             ocufile=' '; goto 990
          endif
! This is written by the gtpread subroutine
!          kl=len_trim(text)
!          if(kl.gt.1) then
!             write(kou,8110)text(1:kl)
!          endif
!8110      format(/'Savefile text: ',a/)
! if there is an assessment record set nvcoeff ...
          if(allocated(firstash%eqlista)) then
!             write(*,*)'Reading the assessment record'
             if(allocated(firstash%coeffvalues)) then
                nvcoeff=0
                kl=size(firstash%coeffvalues)-1
                do j1=0,kl
                   if(firstash%coeffstate(j1).ge.10) then
                      nvcoeff=nvcoeff+1
                   endif
                enddo
                write(kou,3730)nvcoeff
             else
                write(*,*)'No coefficients allocated'
             endif
          endif
!---------------------------------------------------------
       case(2) ! read TDB
          if(tdbfile(1:1).ne.' ') then
! set previous tdbfil as default
             text=tdbfile
             call gparcd('File name: ',cline,last,1,tdbfile,text,q1help)
          else
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT, 8=LOG
! negative is for write, 0 read without filter, -100 write without filter
             call gparfile('File name: ',cline,last,1,tdbfile,' ',1,q1help)
!             call gparc('File name: ',cline,last,1,tdbfile,' ',q1help)
          endif
! if tdbfle starts with "ocbase/" replace that with content of ocbase!!
!          write(*,*)'PMON tdbfile: ',trim(tdbfile)
! check for replacement of OCBASE probably redundant now ...
          name1=tdbfile(1:7)
          call capson(name1)
          if(name1(1:7).eq.'OCBASE/' .or. name1(1:8).eq.'OCBASE\ ') then
             tdbfile=trim(ocbase)//tdbfile(7:)
             write(*,*)'database file: ',trim(tdbfile)
          endif
! this call checks the file exists and returns the elements
! it also lists the DATABASE_INFO text
          call checkdb2(tdbfile,'.tdb',jp,ellist)
          if(gx%bmperr.ne.0) then
             write(kou,*)'No database with this name'
             tdbfile=' '
             goto 990
          elseif(jp.eq.0) then
             write(kou,*)'No elements in the database'
             tdbfile=' '
             goto 100
          endif
          write(kou,8203)jp,(ellist(kl),kl=1,jp)
8203      format('Database has ',i2,' elements: ',18(a,1x)/(1x,28(1x,a)))
          ellist='  '
          write(kou,8205)
8205      format('Give the elements to select, finish with empty line')
          jp=1
          selection='Select elements /all/:'
8210      continue
          call gparc(selection,cline,last,1,ellist(jp),' ',q1help)
          if(jp.eq.1 .and. cline(1:4).eq.'all ') then
! this is if someone actually types "all".  If he types "ALL" that will be AL
             jp=0
          elseif(cline(1:1).eq.'q' .or. cline(1:1).eq.'Q' .or.&
               cline(1:4).eq.'NONE') then
! if user regets selection he can quit
             write(*,*)'Quitting, nothing selected'
             goto 100
          elseif(ellist(jp).ne.'  ') then
             call capson(ellist(jp))
             jp=jp+1
             if(jp.gt.size(ellist)) then
                write(kou,*)'Max number of elements selected: ',size(ellist)
             else
                ll=last
                if(eolch(cline,last)) then
! if empty line list current selection and prompt for more
                   write(*,8220)jp-1,(ellist(iel),iel=1,jp-1)
                else
! we must reset position in cline if there is more ...
                   last=ll
                endif
                selection='Select elements /no more/:'
                goto 8210
             endif
          else
             jp=jp-1
          endif
          if(jp.eq.0) then
             write(kou,*)'All elements selected'
          else
             write(*,8220)jp,(ellist(iel),iel=1,jp)
8220         format('Selected ',i2,' elements: ',20(a,1x))
          endif
          call readtdb(tdbfile,jp,ellist)
          if(gx%bmperr.ne.0) then
! inside readtdb any "buperr" will be set as gx%bmperr
             write(kou,*)'There were errors reading the TDB file', gx%bmperr
             if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
                write(kou,*)bmperrmess(gx%bmperr)
             endif
             write(kou,*)'Please correct these before continuing'
! ignore any type ahead
             last=len(cline)
             call gparcd('Do you want to continue anyway?',&
                  cline,last,1,ch1,'N',q1help)
             if(ch1.ne.'Y') then
                stop 'Good luck fixing the TDB file'
             endif
          endif
! also list the bibliography
          write(kou,*)
          call list_bibliography(' ',kou)
          write(kou,*)
          if(firsteq%multiuse.ne.0) then
             write(*,8221)
8221         format(/' *** There were warnings from reading the database'/&
                  ' *** If you run a macro file please scroll back and check!'/)
          endif
!-----------------------------------------------------------
!8300      continue
       case(3) ! read quit
          goto 100
!-----------------------------------------------------------
       case(4) ! read direct
          write(*,*)'Read direct not implemented yet'
!-----------------------------------------------------------
! this should be merged with read TDB
       case(5) ! read PDB 
          if(tdbfile(1:1).ne.' ') then
             text=tdbfile
             call gparcd('File name: ',cline,last,1,tdbfile,text,q1help)
          else
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
             call gparfile('File name: ',cline,last,1,tdbfile,' ',6,q1help)
!             call gparc('File name: ',cline,last,1,tdbfile,' ',q1help)
          endif
! this call checks the file exists and returns the elements
          call checkdb2(tdbfile,'.pdb',jp,ellist)
          if(gx%bmperr.ne.0) then
             write(kou,*)'No PDB database with this name'
             goto 990
          elseif(jp.eq.0) then
             write(Kou,*)'No elements in the database'
          endif
          write(kou,8203)jp,(ellist(kl),kl=1,jp)
          write(kou,8205)
          jp=1
          selection='Select elements /all/:'
8217      continue
          call gparc(selection,cline,last,1,ellist(jp),' ',q1help)
          if(ellist(jp).ne.'  ') then
             call capson(ellist(jp))
             jp=jp+1
             if(jp.gt.size(ellist)) then
                write(kou,*)'Max number of elements selected: ',size(ellist)
             else
                selection='Select elements /no more/:'
                goto 8217
             endif
          else
             jp=jp-1
          endif
          if(jp.eq.0) then
             write(kou,*)'All elements selected'
          else
             write(*,8220)jp,(ellist(iel),iel=1,jp)
          endif
! later we can add possible options
          name1=' '
          call readpdb(tdbfile,jp,ellist,name1)
! also list the bibliography
          call list_bibliography(' ',kou)
          write(kou,*)
!-----------------------------------------------------------
       case(6) ! read ??
          write(*,*)'Nothing yet'
       end SELECT read
!=================================================================
! SAVE in various formats (NOT MACRO and LATEX, use LIST DATA)
! It is a bit inconsistent as one READ TDB but not SAVE TDB ...
!        ['TDB             ','SOLGAS          ','QUIT            ',&
!         'DIRECT          ','UNFORMATTED     ','PDB             ']
    CASE(9)
! default is 3, unformatted
       kom2=submenu(cbas(kom),cline,last,csave,ncsave,5)
       if(kom2.le.0 .or. kom2.gt.ncsave) goto 100
!
!       if(kom2.gt.3) then
! removed this comment line altogether in version 5.019
! Do not ask this question for QUIT, TDB and SOLGASMIX files
!          call gparc('Comment line: ',cline,last,5,model,' ',q1help)
!       endif
       call date_and_time(optres,name1)
! optres(1:8) is year+month+day, name1(1:4) is hour and minutes
       model=' '//optres(1:4)//'.'//optres(5:6)//'.'//optres(7:8)//&
            ' '//name1(1:2)//'h'//name1(3:4)//' '
!       write(*,*)'comment text: ',trim(model)
       save: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'save subcommand error'
!-----------------------------------------------------------
       case(1) ! save TDB, same as list data TDB
! format 2 is TDB, see list data ...
          kom3=2
          call list_many_formats(cline,last,kom3,kou)
          if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
             write(kou,*)bmperrmess(gx%bmperr)
          elseif(gx%bmperr.ne.0) then
             write(kou,*)'Error code ',gx%bmperr
          endif
!-----------------------------------------------------------
       case(2) ! SOLGAS
          text=' '
! Give a warning that this must not be run on a LINUX computer
          write(*,'(/a/)')' WARNING: Do not run on LINUX/MAC'//&
               ' because END-OF-LINE different from Windows'
!' WARNING: Do not run on LINUX/MAC because END-OF-LINE different from Windows'
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
          call gparfile('File name: ',cline,last,1,filename,text,-7,q1help)
!          call gparc('File name: ',cline,last,1,filename,text,q1help)
          kl=max(index(filename,'.dat '),index(filename,'.DAT '))
          if(kl.le.0) then
             kl=len_trim(filename)+1
             if(kl.eq.1) then
                write(*,*)'Too short file name'
                goto 100
             endif
             filename(kl:)='.DAT '
          endif
          kl=1
          call save_datformat(filename,version,kl,ceq)   
!-----------------------------------------------------------
       case(3) ! save quit, do nothing
          continue
!-----------------------------------------------------------
       case(4) ! save DIRECT
          if(ocdfile(1:1).ne.' ') then
             text=ocdfile
             call gparcd('File name: ',cline,last,1,ocdfile,text,q1help)
          else
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
             call gparfile('File name: ',cline,last,1,ocdfile,' ',-4,q1help)
!             call gparc('File name: ',cline,last,1,ocdfile,' ',q1help)
          endif
          jp=0
          kl=index(ocdfile(2:),'.')+1
          if(kl.le.0) then
             jp=len_trim(ocdfile)
          elseif(ocdfile(kl+1:kl+1).eq.' ') then
! just ending a filename with . not accepted as extention
             jp=kl
          endif
          if(kl.le.1 .and. jp.le.0) then
             write(kou,*)'Missing file name, nothing saved'
             goto 100
          endif
          if(jp.gt.0) ocdfile(jp+1:)='.OCD '
          text='M '//model
          call gtpsave(ocdfile,text)
!-----------------------------------------------------------
       case(5) ! save unformatted
132       continue
          if(ocufile(1:1).ne.' ') then
             text=ocufile
             call gparcd('File name: ',cline,last,1,ocufile,text,q1help)
          else
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
             call gparfile('File name: ',cline,last,1,ocufile,' ',-2,q1help)
!             call gparc('File name: ',cline,last,1,ocufile,' ',q1help)
          endif
          jp=0
! ignore first letter as in macro files a file name may start with ./
          kl=index(ocufile(2:),'.')+1
          if(kl.le.1) then
             jp=len_trim(ocufile)
          elseif(ocufile(kl+1:kl+1).eq.' ') then
! just ending a filename with . not accepted as extention
             jp=kl
          endif
          if(kl.le.1 .and. jp.le.0) then
             write(kou,*)'Missing file name, nothing saved'
             goto 100
          endif
! I have no way to handle the extention to upper case ... inside C routine
!          if(jp.gt.0) ocufile(jp+1:)='.ocu '
          if(jp.gt.0) ocufile(jp+1:)='.OCU '
          inquire(file=ocufile,exist=logok)
          if(logok) then
             call gparcd('File exists, overwrite?',cline,last,1,ch1,'N',q1help)
             if(ch1.ne.'Y') then
                write(*,133)
133             format('Please use another file name')
                ocufile=' '
                goto 132
             endif
             write(*,134)trim(ocufile)
134          format(/'Overwriting previous results on ',a)
          endif
          text='U '//model
          call gtpsave(ocufile,text)
!-----------------------------------------------------------
       case(6) ! PDB
          write(kou,*)'PDB output not yet implemented'
          continue
       end SELECT save
!=================================================================
! help ... just list the commands
    case(10)
       call q3help(cline,last,cbas,ncbas)
       goto 100
!=================================================================
! subcommands to INFORMATION ... very little implemented
!        ['ELEMENTS         ','SPECIES         ','PHASE           ',&
!         'QUIT             ','COMPOSITION_SET ','EQUILIBRIUM     ',&
!         'CHANGES          ','                ','                ']
    case(11)
       kom2=submenu(cbas(kom),cline,last,cinf,ninf,7)
       information: select case(kom2)
!-------------------------------------------------------
          CASE DEFAULT
             write(*,*)'Information subcommand error'
!--------------------------------------------------------
          case(1)
             write(kou,210)
210          format('The elements are those from the periodic chart.'/&
                  'Normally the components are the same as the elemets but',&
                  ' the user',/'can define any rthogonal set of species as',&
                  ' components.')
!--------------------------------------------------------
          case(2)
             write(kou,211)
211          format('Species are molecular like aggregates of elements with',&
                  ' fixed stoichiometry.',/'The elements are the simplest',&
                  ' species.  The constituents of a phase',/' are a subset of',&
                  ' the species.')
!--------------------------------------------------------
! phases info
          case(3)
             write(*,*)'Not written yet'
!--------------------------------------------------------
! quit
          case(4)
!--------------------------------------------------------
! composition set
          case(5)
             write(*,*)'Not written yet'
!--------------------------------------------------------
! equilibrium
          case(6)
             write(*,*)'Not written yet'
!--------------------------------------------------------
! changes
          case(7)
             open(31,file='changes.txt ',access='sequential',err=990,&
                  iostat=buperr)
             do while(.TRUE.)
                do i1=1,40
                   read(31,17,end=244,err=990)line
                   write(kou,17)trim(line)
                enddo
                write(kou,*)'Press return to continue'
                read(kiu,17)ch1
             enddo
244          close(31)
!--------------------------------------------------------
          case(8) ! none
             write(*,*)'Not written yet'
!--------------------------------------------------------
          case(9) ! none
             write(*,*)'Not written yet'
          end select information
       goto 100
!=================================================================
! back / goto, return to calling (main) program
    case(12)
!       call gparcd('Are you sure?',cline,last,1,ch1,'N',q1help)
!       if(ch1.eq.'y' .or. ch1.eq.'Y') then
!          write(*,*)'Welcome back!'
!       endif
       return
!=================================================================
! NEW command, same as reinitiate
    case(13)
! one must deallocate everyting explicitly to use memory again
       call gparcd('All data will be removed, are you sure?',cline,last,&
            1,ch1,'N',q1help)
       if(ch1.ne.'Y') then
          write(kou,*)'*** NO CHANGE, upper case Y needed for NEW'
          goto 100
       endif
!------remove assessment data
!       write(*,*)'No segmentation fault 1'
       if(allocated(firstash%eqlista)) then
          write(*,*)'Assessment data removed, not deallocated: memory leak'
       endif
!       write(*,*)'No segmentation fault 2'
       if(allocated(firstash%eqlista)) deallocate(firstash%eqlista)
       deallocate(firstash)
!       write(*,*)'No segmentation fault 3, deallocate errs '
       if(mexp.gt.0) deallocate(errs)
       mexp=0
       updatemexp=.true.
       nvcoeff=0
!       iexit=0
!       iexit(2)=1
!       write(*,*)'No Segmentation fault 4'
!----- deleting map results ...
!       write(*,*)'Deleting map results'
       if(associated(maptopsave)) then
! this is necessary only if no plot of last step/map made ...
!          write(kou,*)'We link to maptopsave'
          maptop%plotlink=>maptopsave
          nullify(maptopsave)
       endif
       seqxyz=0
       call delete_mapresults(maptop)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error deleting map results! Report this error with macro!'
          stop
       endif
!       write(*,*)'Back from delete_mapresults'
! remove any results from step and map
!       if(associated(maptop)) then
!          write(*,*)'maptop nullified: ',maptop%next%seqx
!          maptop%next%seqx=0
!          maptop%next%seqy=0
!          maptop%seqx=0
!          maptop%seqy=0
!          nullify(maptop)
!       endif
!       write(*,*)'we are here'
       nullify(maptop)
       nullify(mapnode)
       nullify(maptopsave)
       seqxyz=0
!----- deallocate local axis records
       do jp=1,noofaxis
          if(allocated(axarr(jp)%axcond)) deallocate(axarr(jp)%axcond)
!          deallocate(axarr(jp)%indices)
!          deallocate(axarr(jp)%coeffs)
       enddo
       noofaxis=0
! remove some more defaults ...
       defcp=1
! deallocate does not work on pointers!!!
       nullify(starteq)
       noofstarteq=0
!
! this routine fragile, inside new_gtp init_gtp is called
!       write(*,*)'No segmentation fault 7'
       call new_gtp
       if(gx%bmperr.ne.0) then
          write(*,*)'Error deleting data! Report this error with macro!'
          stop
       endif
       write(kou,*)'All data removed'
       call init_gtp(intv,dblv)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error initiating! Report this error with macro!'
          stop
       endif
       write(kou,*)'Workspaces initiated'
!       ceq=>firsteq
       goto 20
!=================================================================
! macro begin
    case(14) ! file name is asked inside macbeg
       call macbeg(cline,last,logok)
       if(buperr.ne.0 .or. gx%bmperr.ne.0) goto 990
!=================================================================
! about
    case(15)
       write(kou,15010)version,linkdate
15010  format(/'This is OpenCalphad (OC), a free software for ',&
            'thermodynamic calculations as'/&
            'described in B Sundman, U R Kattner, M Palumbo and S G Fries, ',&
            'Integrating'/'Materials and Manuf. Innov. (2015) 4:1; ',&
            'B Sundman, X-G Lu and H Ohtani,'/'Comp Mat Sci, Vol 101 ',&
            '(2015) 127-137 and B Sundman et al., Comp Mat Sci, '/&
            'Vol 125 (2016) 188-196'//&
            'It is available for download at http://www.opencalphad.org or'/&
            'the sundmanbo/opencalphad repository at http://www.github.com'//&
            'This software is protected by the GNU General Public License'/&
            'You may freely distribute copies as long as you also provide ',&
            'the source code'/'and use the GNU GPL license also for your own',&
            ' additions and modifications.'//&
            'The software is provided "as is" without any warranty of any ',&
            'kind, either'/'expressed or implied.  ',&
            'The full license text is provided with the software'/&
            'or can be obtained from the Free Software Foundation ',&
            'http://www.fsf.org'//&
            'Copyright 2011-2018, Bo Sundman, Gif sur Yvette, France.'/&
            'Contact person Bo Sundman, bo.sundman@gmail.com'/&
            'This version ',a,' was compiled ',a/)
!=================================================================
! debug subcommands
    case(16)
       kom2=submenu(cbas(kom),cline,last,cdebug,ncdebug,1)
       debug: SELECT CASE (kom2)
!------------------------------
       CASE DEFAULT
          write(kou,*)'Debug subcommand error ',kom2
!------------------------------
! debug free lists
       CASE(1)
          write(*,*)'Check components masses'
          call compmassbug(ceq)
!          write(*,*)'Calculating equilibrium record size'
!          kom3=ceqsize(ceq)
!          write(kou,*)'Current equilibrium record memory use: ',kom3
! list all tuples
          write(kou,1617)
1617      format('Phase tuples content:'/&
               'Tuple lokph   compset ixphase lokvares nextcs phase name',&
               '       disfra vareslink')
          do jp=1,nooftup()
             call get_phasetup_name(jp,name1)
! this is a check that %ihaseix and lokvares are correct
!             if(phasetuple(jp)%compset.eq.1) then
!                call get_phase_compset(jp,1,lokph,lokcs)
!             else
!                call get_phase_compset(phasetuple(jp)%ixphase,&
!                     phasetuple(jp)%compset,lokph,lokcs)
!             endif
!             write(kou,16020)jp,phasetuple(jp),name1,lokph,lokcs
             write(kou,16020)jp,phasetuple(jp),name1,&
                  ceq%phase_varres(phasetuple(jp)%lokvares)%disfra%varreslink
!16020        format(i3,': ',2i7,2i9,3x,a/i12,18x,i7)
16020        format(i3,': ',2i7,2i9,i6,3x,a,2x,i6)
          enddo
          call list_free_lists(kou)
!------------------------------
! debug stop_on_error
       CASE(2)
          if(stop_on_error) then
             stop_on_error=.FALSE.
             write(kou,*)'No longer stop on error'
          else
             stop_on_error=.true.
          endif
!------------------------------
! debug elasticity (this is temporary)
       CASE(3)
          write(*,*)'Not implemented yet'
!          write(kou,*)'Input current lattice parameter values (3x3 matrix)',&
!               ' for phase 1'
!          iph=1
!          ics=1
!          xxx=7.1D-6
!          xxy=1.0D-12
!          call gparrd('lattice par (1,1):',cline,last,latpos(1,1),xxx,nohelp)
!          call gparrd('lattice par (1,2):',cline,last,latpos(1,2),xxy,nohelp)
!          call gparrd('lattice par (1,3):',cline,last,latpos(1,3),xxy,nohelp)
!          call gparrd('lattice par (2,1):',cline,last,latpos(2,1),xxy,nohelp)
!          call gparrd('lattice par (2,2):',cline,last,latpos(2,2),xxx,nohelp)
!          call gparrd('lattice par (2,3):',cline,last,latpos(2,3),xxy,nohelp)
!          call gparrd('lattice par (3,1):',cline,last,latpos(3,1),xxy,nohelp)
!          call gparrd('lattice par (3,2):',cline,last,latpos(3,2),xxy,nohelp)
!          call gparrd('lattice par (3,3):',cline,last,latpos(3,3),xxx,nohelp)
!          call set_lattice_parameters(iph,ics,latpos,ceq)
!----------------------------------
! debug test1 (whatever)
       case(4)
!          write(*,*)'Nothing here'
          do i1=1,nosp()
             call get_species_location(i1,loksp,name1)
             if(gx%bmperr.ne.0) goto 990
             call get_species_component_data(loksp,i2,iphl,stoik,xxx,&
                  xxy,ceq)
             if(gx%bmperr.ne.0) goto 990
             write(kou,1670)i1,loksp,name1,xxx,xxy,(iphl(j1),stoik(j1),j1=1,i2)
1670         format(2i4,1x,a12,1x,2F6.2,2x,10(i3,1x,F7.4))
          enddo
!---------------------------------
! debug tpfun (whatever)
       case(5)
          call gparid('Function index:',cline,last,ll,-1,nohelp)
          call list_tpfun_details(ll)
!---------------------------------
! debug browser
       case(6)
! testing using HTML helpfile with "anchors" like <a name="label" />
!   related to a question or command
! and the the help utility will search for a specific label as below
! NOTE in the LaTeX file \usepackage{hyperref}
! and in the text \hypertarget{selectname}
! using "path/browser" "file://path/helpfile#label" should position
! the html window at label!!
! the label "selectname" is in the html file ...
          call gparcd('File name: ',cline,last,5,model,&
               './manual\html\ochelp.html#selectelement ',q1help)
!          browser='"C:\Program Files\Mozilla firefox\firefox.exe" '
! this browser can be opened without ""
          browser='C:\PROGRA~1\INTERN~1\iexplore.exe '
! it works to open the ochelp.html on the first page
!          string=trim(browser)//&
!               ' -file ./manual/html/ochelp.html'
!          write(*,'(a)')trim(string)
! gnu fortran ...
!          call system(...)
!          call execute_command_line(string)
! now the complicated one ...
          string=trim(browser)//&
               '     "file://C:\Users\bosse\documents\oc\oc\src\'//&
               trim(model(3:))//'"'
          write(*,'(a)')trim(string)
! gnu fortran ...
!          call system(...)
          call execute_command_line(string)
! This command works in a Windows terminal window:
! "C:\program files\Mozilla firefox\firefox.exe" 
!  "file://c:\users\bosse\documents\oc\oc\src\manual\html\ochelp.html#selectelement"
! but problems using as command ...
! This works also:
!c:\Program Files\internet explorer\iexplore.exe" "file://c:\users\bosse\documents\oc\oc\src\manual\html\ochelp.html#selectelement"
! maybe possible to access by directory names with only 8 bytes ...
!linux!linux!linux!linux!linux!linux!linux!linux!linux!linux!linux!linux!linux
!
!>          call gparcd('File name: ',cline,last,5,model,&
!>            '/home/bosse/OC/OC5-20/manual/ochelp.html#selectelement ',q1help)
! this browser can be opened without ""
!>          browser='/usr/bin/firefox '
!>          string=trim(browser)//' "file:'//trim(model(1:))//'"'
!>          write(*,'(a)')trim(string)
! This command works in a Linux terminal window:
! /usr/bin/firefox -file /home/bosse/OC/OC5-20/manual/ochelp.html
! it does not work to add #selectelement at the end (no such file name)
! This work in Linux terminal window:
! /usr/bin/firefox "file:/home/bosse/.../manual/ochelp.html#selectelement"
!
!linux!linux!linux!linux!linux!linux!linux!linux!linux!linux!linux!linux!linux
!---------------------------------
! debug trace
       case(7)
             call gparcd('HTML help?',cline,last,1,ch1,'Y',q1help)
             if(ch1.eq.'Y') then
                helptrace=.TRUE.
             else
                helptrace=.FALSE.
             endif
             call gparcd('plotting?',cline,last,1,ch1,'N',q1help)
             if(ch1.eq.'Y') then
                plottrace=.TRUE.
             else
                plottrace=.FALSE.
             endif
!..................................
! not used
       case(8)
!..................................
! not used
       case(9)
       END SELECT debug
!=================================================================
! select command
    case(17)
       kom2=submenu(cbas(kom),cline,last,cselect,nselect,1)
       selct: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'Select subcommand error'
          goto 100
!-----------------------------------------------------------
       CASE(1) ! select equilibrium
          if(ceq%eqno.lt.noeq()) then
             name1='NEXT'
          else
             name1='DEFAULT'
          endif
          call gparcd('Give name or number?',cline,last,1,text,&
               name1,q1help)
          if(buperr.ne.0) goto 990
! if the user types "next" in lower case or an abbrev it does not work
          call capson(text)
          if(compare_abbrev(text,'NEXT ')) then
             i1=ceq%eqno+1
             call selecteq(i1,ceq)
             if(gx%bmperr.ne.0) goto 990
             neqdef=i1
          elseif(compare_abbrev(text,'PREVIOUS ')) then
             i1=max(ceq%eqno-1,1)
             call selecteq(i1,ceq)
             if(gx%bmperr.ne.0) goto 990
             neqdef=i1
          elseif(compare_abbrev(text,'LAST ')) then
             i1=noeq()
             call selecteq(i1,ceq)
             if(gx%bmperr.ne.0) goto 990
             neqdef=i1
          elseif(compare_abbrev(text,'FIRST ')) then
             i1=1
             call selecteq(i1,ceq)
             if(gx%bmperr.ne.0) goto 990
             neqdef=i1
          else
! check if number
             j1=1
             call getint(text,j1,i1)
             if(buperr.ne.0) then
                buperr=0
! findeq accepts PREVIOUS and FIRST (same as DEFAULT)
                ieq=ceq%eqno
                call findeq(text,ieq)
                if(gx%bmperr.ne.0) goto 990
                neqdef=ieq
                ceq=>eqlista(ieq)
             else
                call selecteq(i1,ceq)
                if(gx%bmperr.ne.0) goto 990
                neqdef=i1
             endif
          endif
          write(kou,*)'Current equilibrium ',ceq%eqname
!-----------------------------------------------------------
       CASE(2) ! select minimizer
          write(kou,*)'Sorry, only one available: ',minimizers(2)
!          call gparcd('Give name of minimizer?',&
!               cline,last,1,text,'LUKAS ',q1help)
!          call capson(text)
!          j1=len_trim(text)
!  write(*,*)'input: ',text(1:16),j1
!          if((text(1:j1).eq.minimizers(1)(1:j1))) then
!             minimizer=1
!          else
!             minimizer=2
!          endif
          write(kou,*)'Selected minimizer: ',minimizers(minimizer)
!-----------------------------------------------------------
       case(3) ! select graphics
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
       case(4) ! select language, at present only 1 English and 2 French
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
       case(5) ! select optimizer
          write(kou,845)optimizers(optimizer)
          write(kou,844)optimizers
844       format('Available optimizers: '/,(2x,a,2x,a,2x,a))
845       format('Current optimizer is: '/,2x,a)
          call gparcd('Do you want to use LMDIF?',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'Y') then
             optimizer=1
          else
             write(*,*)'Sorry VA05AD is no longer available'
!             call gparcd('Do you want to use VA05AD?',&
!                  cline,last,1,ch1,'Y',q1help)
!             if(ch1.eq.'Y') then
!                optimizer=2
!             endif
          endif
          write(kou,*)'You have selected ',optimizers(optimizer)
!-----------------------------------------------------------
       case(6)
          goto 100
       END SELECT selct
!=================================================================
! DELETE not much implemented ...
!         ['ELEMENTS        ','SPECIES         ','PHASE           ',&
!          'QUIT            ','COMPOSITION_SET ','EQUILIBRIUM     ',&
!          'STEP_MAP_RESULTS','                ','                ']
    CASE(18)
! disable continue optimization
!       iexit=0
!       iexit(2)=1
       kom2=submenu(cbas(kom),cline,last,crej,nrej,6)
       delete: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'Delete subcommand error'
          goto 100
!-----------------------------------------------------------
! delete element
       case(1)
          write(kou,18010)
18010     format(' *** Warning, this command will delete the data for the',&
            ' element, species or'/' phase specified and the data cannot',&
            ' be recovered unless read again from'/' file.  If you',&
            ' only want to temporarily remove some data use QUIT'/&
            ' from this command and then SET STATUS'/)
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! delete species
       case(2)
          write(kou,18010)
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! delete phase
       case(3)
          write(kou,18010)
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! quit
       case(4)
          goto 100
!-----------------------------------------------------------
! delete composition set, always that with higest number
       case(5)
          call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
          if(buperr.ne.0) goto 990
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
          call remove_composition_set(iph,.FALSE.)
          if(gx%bmperr.ne.0) goto 990
!-----------------------------------------------------------
! delete equilibria
       case(6)
          call gparc('Equilibrium name or abbr.:',cline,last,1,name1,' ',q1help)
          if(buperr.ne.0) goto 990
          call delete_equilibria(name1,ceq)
          if(gx%bmperr.ne.0) goto 990
!-----------------------------------------------------------
! delete step_map_results
       case(7)
          if(associated(maptopsave)) then
! this is necessary only if no plot of last step/map made ...
             maptop%plotlink=>maptopsave
             nullify(maptopsave)
             write(*,*)'maptopsave nullified'
          endif
          seqxyz=0
! this does not delete _mapnode and _mapline equilibria ???
          call delete_mapresults(maptop)
! remove any results from step and map
!          if(associated(maptop)) then
!             write(*,*)'maptop nullified: ',maptop%next%seqx
!             maptop%next%seqx=0
!             maptop%next%seqy=0
!             maptop%seqx=0
!             maptop%seqy=0
!             nullify(maptop)
!          endif
          nullify(maptop)
          nullify(mapnode)
          nullify(maptopsave)
!----- deallocate local axis records
          do jp=1,noofaxis
             if(allocated(axarr(jp)%axcond)) deallocate(axarr(jp)%axcond)
          enddo
          noofaxis=0
! remove some more defaults ...
          defcp=1
! deallocate does not work on pointers!!!
          nullify(starteq)
          noofstarteq=0
          call reset_plotoptions(graphopt,plotfile,textlabel)
          axplotdef=' '
!-----------------------------------------------------------
!
       case(8)
          continue
!-----------------------------------------------------------
!
       case(9)
          continue
       end SELECT delete
!=================================================================
! STEP, must be tested if compatible with assessments
!         ['NORMAL          ','SEPARATE        ','QUIT            ',&
!          'CONDITIONAL     ','                ','                ']
    case(19)
! disable continue optimization
!       iexit=0
!       iexit(2)=1
       if(noofaxis.ne.1) then
          write(kou,*)'You must set exactly one independent axis variable',&
               ' for a step calculation.'
          goto 100
       endif
       ll=degrees_of_freedom(ceq)
       if(ll.ne.0) then
          write(*,*)'Degrees of freedom not zero',ll
          goto 100
       endif
! IMPORTANT I have changed the order between option and reinitiate!!
       kom2=submenu('Options?',cline,last,cstepop,nstepop,1)
! check if adding results
       if(associated(maptop)) then
          write(kou,833)
833       format('There are previous results from step or map')
          call gparcd('Delete them?',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
! there should be a more careful deallocation to free memory
             call delete_mapresults(maptop)
!             deallocate(maptop%saveceq)
             nullify(maptop)
             nullify(maptopsave)
             write(kou,*)'Previous results removed'
! delete equilibria associated with STEP/MAP
             call delete_equilibria('_MAP*',ceq)
             seqxyz=0
! remove all graphopt settings
             call reset_plotoptions(graphopt,plotfile,textlabel)
             axplotdef=' '
          else
! for step separate it seems difficult to have correct seqx !!
!             seqxyz(1)=maptop%next%seqx
             seqxyz(1)=max(maptop%next%seqx,maptop%previous%seqx,maptop%seqx)
             seqxyz(2)=maptop%seqy
! list maptop for debugging
!             write(*,*)'PM maptop node: ',trim(maptop%nodeceq%eqname)
!             maptopcheck=>maptop%next
!             do while(.not.associated(maptopcheck,maptop))
!                write(*,*)'PM: maptop node: ',trim(maptopcheck%nodeceq%eqname)
!                maptopcheck=>maptopcheck%next
!                if(.not.associated(maptopcheck%previous%next,maptopcheck)) then
!                   write(*,*)'PM next and previous does not agree'
!                endif
!             enddo
!             if(associated(maptop%plotlink)) then
!                write(*,*)'PM plotlink: ',trim(maptop%plotlink%nodeceq%eqname)
!             endif
! end debugging
             maptopsave=>maptop
             nullify(maptop)
             write(*,'(a,2i4)')'Previous results kept',seqxyz
          endif
       endif
!       kom2=submenu('Options?',cline,last,cstepop,nstepop,1)
       step: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'No such step option'
!-----------------------------------------------------------
! STEP NORMAL
       case(1)
!       call gparcd('Option?',cline,last,1,text,'NORMAL ',q1help)
!       if(associated(resultlist)) then
! maptop is returned as main map/step record for results
! noofaxis is current number of axis, axarr is array with axis data
! starteq is start, equilibria, if empty set it to ceq
          if(.not.associated(starteq)) then
             starteq=>ceq
          endif
! can one have several STEP commands YES!
          if(associated(maptop)) then
             write(*,*)'Deleting previous step/map results missing'
             goto 100
          endif
! seqzyz are initial values for creating equilibria for lines and nodes
          call map_setup(maptop,noofaxis,axarr,seqxyz,starteq)
! mark that interactive listing of conditions and results may be inconsistent
          ceq%status=ibset(ceq%status,EQINCON)
          if(.not.associated(maptop)) then
! if one has errors in map_setup maptop may not be initiated, if one
! has saved previous calculations in maptopsave restore those
             if(associated(maptopsave)) then
                write(kou,*)'Restoring previous map results'
                maptop=>maptopsave
                nullify(maptopsave)
             endif
          elseif(associated(maptopsave)) then
! THIS ERROR WITH LOST CALCULATONS MAY BE THERE FOR STEP SEPERATE AND MAP
!             write(*,*)'PM linking previous results'
             write(kou,*)'Set link to previous step results'
             maptop%plotlink=>maptopsave
          endif
! debugging: last maptop/line used
!          write(*,'(a,2i4)')'PMON: sexy 1:',maptop%next%seqx,maptop%seqy
! remove start equilibria
          nullify(starteq)
          if(gx%bmperr.ne.0) goto 990
!-----------------------------------------------------------
! STEP SEPARATE
       case(2) ! calculate for each entered phase separately
          starteq=>ceq
! can one have several STEP commands??
          if(associated(maptop)) then
             write(*,*)'Deleting previous step/map results missing'
             goto 100
          endif
          call step_separate(maptop,noofaxis,axarr,seqxyz,starteq)
! mark that interactive listing of conditions and results may be inconsistent
          ceq%status=ibset(ceq%status,EQINCON)
          if(.not.associated(maptop)) then
! if one has errors in map_setup maptop may not be initiated, if one
! has saved previous calculations in maptopsave restore those
             if(associated(maptopsave)) then
                write(kou,*)'Restoring previous map results'
!                maptop=>maptopsave
                maptop%plotlink=>maptopsave
                nullify(maptopsave)
             endif
          elseif(associated(maptopsave)) then
             write(kou,*)'Set link to previous step results'
             maptop%plotlink=>maptopsave
          endif
! set default yaxis as GM(*)
          if(axplotdef(2)(1:1).eq.' ') then
             axplotdef(2)='GM(*)'
          endif
! update maptop%seqx to maptop%prvious%seqx+1 to allow more maptop records
          maptop%seqx=maptop%previous%seqx+1
!          write(*,'(a,4i4)')'PMON: separate seqx:',maptop%next%seqx,&
!               maptop%seqx,maptop%previous%seqx,maptop%seqy
! remove start equilibria
          nullify(starteq)
!-----------------------------------------------------------
! STEP QUIT
       case(3)
!-----------------------------------------------------------
! STEP CONDITIONAL
       case(4)
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! STEP unused
       case(5)
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! STEP unused
       case(6)
          write(kou,*)'Not implemented yet'
       end SELECT step
!=================================================================
! MAP, must be tested if compatible with assessments
    case(20)
! maybe disable continue optimization ??
       if(noofaxis.lt.2) then
          write(kou,*)'You must set two axis with independent variables'
          goto 100
       endif
       write(kou,20014)
20014   format('The map command is fragile, please send problematic diagrams',&
            ' to the',/'OC development team'/)
       if(associated(maptop)) then
          write(kou,833)
          call gparcd('Reinitiate?',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
             call delete_mapresults(maptop)
!             deallocate(maptop%saveceq)
             nullify(maptop)
             nullify(maptopsave)
! this removes all previous equilibria associated with STEP/MAP commands
! already done by delete_mapresults
!             call delete_equilibria('_MAP*',ceq)
!             if(gx%bmperr.ne.0) then
!                write(kou,*)'Error removing old MAP equilibria'
!                goto 990
!             endif
! initiate indexing nodes and lines
             seqxyz=0
! remove all graphopt settings
             call reset_plotoptions(graphopt,plotfile,textlabel)
             axplotdef=' '
          else
! start indexing new nodes/lines from previous 
!             write(*,*)'mapnode: ',maptop%seqx,maptop%previous%seqx,&
!                  maptop%next%seqx
             seqxyz(1)=maptop%next%seqx
             seqxyz(2)=maptop%seqy
!             seqxyz(3) can be used for something else ...
             maptopsave=>maptop
             nullify(maptop)
!             write(*,*)'seqxyz: ',seqxyz
          endif
! this should never be done ! It destroys the possibility to find old nodes
!          call delete_equilibria('_MAP*',ceq)
       endif
! maptop is returned as main map/step record for results
! noofaxis is current number of axis, axarr is array with axis data
! starteq is start equilibria, if empty set it to ceq
       if(.not.associated(starteq)) then
          starteq=>ceq
          starteq%nexteq=0
       endif
       ll=degrees_of_freedom(starteq)
       if(ll.ne.0) then
          write(*,*)'Degrees of freedom not zero ',ll
          goto 100
       endif
! maptop is first nullified inside map_setup, then alloctated to return result
       call map_setup(maptop,noofaxis,axarr,seqxyz,starteq)
       if(.not.associated(maptop)) then
! if one has errors in map_setup maptop may not be initiated, if one
! has saved previous calculations in maptopsave restore those
          if(associated(maptopsave)) then
             write(kou,*)'Restoring previous map results'
             maptop=>maptopsave
             nullify(maptopsave)
          endif
       elseif(associated(maptopsave)) then
          write(kou,*)'Set link to previous map results'
          maptop%plotlink=>maptopsave
          nullify(maptopsave)
       endif
! remove start equilibria
       nullify(starteq)
! mark that interactive listing of conditions and results may be inconsistent
       ceq%status=ibset(ceq%status,EQINCON)
       if(gx%bmperr.ne.0) goto 990
! end of MAP command
!=================================================================
! PLOT COMMAND with many options
! Always specify the axis when giving this command, default is previous!!
! Sunbommands comes after
    case(21)
       if(.not.associated(maptop)) then
          write(kou,*)'You must give a STEP or MAP command before PLOT'
          goto 100
       endif
       wildcard=.FALSE.
       do iax=1,2
! default scaling factor for the axis variable
          graphopt%scalefact(iax)=one
          plotdefault: if(axplotdef(iax)(1:1).eq.' ') then
! insert a default answer for plot axis
             if(iax.le.noofaxis) then
                jp=1
                call get_one_condition(jp,text,axarr(iax)%seqz,ceq)
                if(gx%bmperr.ne.0) then
                   write(*,*)'Error getting axis condition from index: ',&
                        iax,axarr(iax)%seqz
                   goto 990
                endif
! we just want the expression, remove the value including the = sign
                jp=index(text,'=')
                text(jp:)=' '
                if(maptop%tieline_inplane.eq.1) then
! if tie-lines in the plane is 1 and calculating axis was x(A)
! then plot axis should be x(*,cu) 
                   jp=index(text,'(')
                   if(jp.gt.0) then 
                      text=text(1:jp)//'*,'//text(jp+1:)
                   endif
                endif
             else
                text='NP(*) '
             endif
             axplotdef(iax)=text
          endif plotdefault
! the 4th argument to gparc means the following:
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
!------------------------------------------------------------------------
21000      continue
          if(iax.eq.1) then
             call gparcd('Horizontal axis variable',&
                  cline,last,7,axplot(iax),axplotdef(iax),q1help)
          else
             call gparcd('Vertical axis variable',&
                  cline,last,7,axplot(iax),axplotdef(iax),q1help)
          endif
          if(buperr.ne.0) goto 990
! extract a possible scaling factor like 0.001*GM(*)
          jp=1
          call getrel(axplot(iax),jp,xxx)
          if(buperr.eq.0) then
! there is a numerical factor
             graphopt%scalefact(iax)=xxx
! a number must be followed by a *
             if(axplot(iax)(jp:jp).ne.'*') then
                write(*,*)'Scaling factor must be followed by *'
                goto 990
             else
! Fortran allows overlapping strings in assignments
                axplot(iax)=axplot(iax)(jp+1:)
             endif
          else
! no scaling factor, graphopt%scalfactor(iax) already unity
             buperr=0
          endif
          if(index(axplot(iax),'*').gt.0) then
!             if(wildcard) then
!                write(*,*)'Wildcards allowed for one axis only'
!                goto 21000
!             else
                wildcard=.TRUE.
!             endif
          endif
          if(axplotdef(iax).ne.axplot(iax)) then
! if new axis variable then reset default plot options
! plot ranges and their defaults
             call reset_plotoptions(graphopt,plotfile,textlabel)
! check that axis variable is a correct state variable or symbol
! Code copied from show variable (case(4,17) around line 3273)
             if(index(axplot(iax),'*').gt.0) then
! generate many values
! the values are returned in yarr with dimension maxconst. 
! longstring are the state variable symbols for the values ...
                call get_many_svar(axplot(iax),yarr,maxconst,i1,longstring,ceq)
                if(gx%bmperr.ne.0) then
! if error go back to command level
                   write(kou,*)'Illegal axis variable!  Error code: ',gx%bmperr
                   goto 100
!                else
!                   write(*,*)'pmon test value: ',yarr(1)
                endif
             else
! the value of a state variable or model parameter variable is returned
! STRANGE the symbol xliqni is accepted in get_state_var_value ???
                call get_state_var_value(axplot(iax),xxx,model,ceq)
                if(gx%bmperr.ne.0) then
! if error check if it is a complicated symbol like CP=H.T
                   gx%bmperr=0
! If error then try to calculate a symbol ...
                   call capson(axplot(iax))
!                   call find_svfun(axplot(iax),istv,ceq)
                   call find_svfun(axplot(iax),istv)
                   if(gx%bmperr.ne.0) then
                      write(kou,*)'Illegal axis variable, error: ',gx%bmperr
                      goto 100
                   endif
                endif
             endif
          endif
! remember most recent axis as default (and to avoid reset)
          axplotdef(iax)=axplot(iax)
       enddo
!       call gparcd('GNUPLOT file name',cline,last,1,plotfile,'ocgnu',q1help)
!       form=' '
! first argument is the number of plot axis, always 2 at present
       jp=2
       if(associated(maptopsave)) then
          write(kou,*)'We link to maptopsave'
          maptop%plotlink=>maptopsave
!       else
!          write(*,*)'There is no maptopsave'
       endif
!-----------------------------------------------------------
! PLOT subcommands, default is PLOT, NONE does not work ...
! subcommands to PLOT OPTIONS/ GRAPHICS OPTIONS
!    character (len=16), dimension(nplt) :: cplot=&
!        ['RENDER          ','SCALE_RANGES    ','RATIOS_XY       ',&
!         'AXIS_LABELS     ','MANIPULATE_LINES','TITLE           ',&
!         'GRAPHICS_FORMAT ','OUTPUT_FILE     ','GIBBS_TRIANGLE  ',&
!         'QUIT            ','POSITION_OF_KEYS','APPEND          ',&
!         'TEXT            ','TIE_LINES       ','FONT_AND_COLOR  ',&
!         'LOGSCALE        ','LINE_WITH_POINTS','PAUSE_OPTION    ']
!         'MISCELLANEOUS   ','                ','                ']
!-------------------
! return here after each subcommand
21100   continue
       if(graphopt%gnutermsel.lt.1 .or. &
            graphopt%gnutermsel.gt.graphopt%gnutermax) then
          write(kou,*)'No such graphics terminal: ',graphopt%gnutermsel
       elseif(graphopt%gnutermsel.ne.1) then
          write(kou,2910)trim(graphopt%gnutermid(graphopt%gnutermsel)),&
               trim(plotfile),trim(graphopt%filext(graphopt%gnutermsel))
2910      format(/'Graphics output as ',a,' on file: ',a,'.',a)
       endif
       write(kou,21112)
21112  format(/'Note: give only one option per line!')
       kom2=submenu('Options?',cline,last,cplot,nplt,1)
       plotoption: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'No such plot option'
!-----------------------------------------------------------
! PLOT RENDER no more options to plot ...
       case(1)
!2190      continue
! use the graphics record to transfer data ...
          graphopt%pltax(1)=axplot(1)
          graphopt%pltax(2)=axplot(2)
          if(graphopt%gibbstriangle) then
! if gibbstriangle make sure min is 0
             graphopt%plotmin(1)=zero
             graphopt%dfltmin(1)=zero
             graphopt%plotmin(2)=zero
             graphopt%dfltmin(2)=zero
             if(graphopt%rangedefaults(1).ne.0 .or. &
                  graphopt%rangedefaults(2).ne.0) then
! if gibbstriangle and scaling make sure xmax and ymax are the same
                xxx=min(graphopt%plotmax(1),graphopt%plotmax(2))
                graphopt%plotmax(1)=xxx
                graphopt%dfltmax(1)=xxx
                graphopt%plotmax(2)=xxx
                graphopt%dfltmax(2)=xxx
             endif
          endif
          graphopt%filename=' '
!          write(*,*)' >>>>>>>>>>>>> ',trim(plotfile)
          graphopt%filename=plotfile
! added ceq in the call to make it possible to handle change of reference states
          call ocplot2(jp,maptop,axarr,graphopt,version,ceq)
          if(gx%bmperr.ne.0) goto 990
! always restore default plot file name and plot option to screem
          if(graphopt%gnutermsel.ne.1) &
               write(kou,*)'Restoring plot device to screen'
          graphopt%gnutermsel=1
          graphopt%filename='ocgnu '
          plotfile='ocgnu'
!          call gparc('Hardcopy (P for postscript)?',&
!               cline,last,1,form,'none',q1help)
!          if(form.ne.'none') then
!             call capson(form)
!             call ocplot2(jp,axplot,plotfile,maptop,axarr,graphopt,form,ceq)
!          endif
!-----------------------------------------------------------
! PLOT SCALE_RANGE of either X or Y
       case(2)
          call gparcd('For X or Y axis? ',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'X' .or. ch1.eq.'x') then
!             if(graphopt%axistype(1).eq.1) then
!                write(kou,*)'The x axis set to linear'
!                graphopt%axistype(1)=0
!             else
!                graphopt%axistype(1)=1
!             endif
             goto 21120
          elseif(ch1.eq.'Y' .or. ch1.eq.'y') then
!             if(graphopt%axistype(2).eq.1) then
!                write(kou,*)'The y axis set to linear'
!                graphopt%axistype(2)=0
!             else
!                graphopt%axistype(2)=1
!             endif
             goto 21130
          else
             write(kou,*)'Please answer X or Y'
          endif
          goto 21100
!............................................ user limits X axis (1)
21120     continue
          call gparcd('Default limits',cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'Y' .or. ch1.eq.'y') then
             graphopt%rangedefaults(1)=0
          else
             graphopt%rangedefaults(1)=1
             twice=.FALSE.
21104        continue
             call gparrd('Low limit',cline,last,xxx,graphopt%dfltmin(1),q1help)
             if(graphopt%gibbstriangle .and. xxx.ne.zero) then
                write(*,*)'Lower limit of a Gibbs triangle plot must be zero'
                goto 21100
             endif
             graphopt%plotmin(1)=xxx
             graphopt%dfltmin(1)=xxx
             once=.TRUE.
21105        continue
             call gparrd('High limit',cline,last,xxx,&
                  graphopt%dfltmax(1),q1help)
             if(xxx.le.graphopt%plotmin(1)) then
                if(once) then
                   write(kou,*)'Think before typing'
                   once=.FALSE.
                elseif(twice) then
                   write(kou,*)'Back to command level'
                   goto 100
                else
                   write(kou,*)'Please give the low limit again!'
                   twice=.TRUE.
                   goto 21104
                endif
                write(kou,21106)graphopt%plotmin(1)
21106           format('High limit must be higher than low: ',1pe14.6)
                goto 21105
             endif
             graphopt%plotmax(1)=xxx
             graphopt%dfltmax(1)=xxx
          endif
          goto 21100
!---------------------------------------------- user limits Y axis (2)
21130     continue
          call gparcd('Default limits',cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'Y' .or. ch1.eq.'y') then
             graphopt%rangedefaults(2)=0
          else
             graphopt%rangedefaults(2)=1
             twice=.FALSE.
21107        continue
             call gparrd('Low limit',cline,last,xxx,graphopt%dfltmin(2),q1help)
             if(graphopt%gibbstriangle .and. xxx.ne.zero) then
                write(*,*)'Lower limit of a Gibbs triangle plot must be zero'
                goto 21100
             endif
             graphopt%plotmin(2)=xxx
             graphopt%dfltmin(2)=xxx
             once=.TRUE.
21108        continue
             call gparrd('High limit',cline,last,xxx,&
                  graphopt%dfltmax(2),q1help)
             if(xxx.le.graphopt%plotmin(2)) then
                if(once) then
                   write(*,*)'Think before typing'
                   once=.FALSE.
                elseif(twice) then
                   write(kou,*)'Back to command level'
                   goto 100
                else
                   write(kou,*)'Please give the low limit again!'
                   twice=.TRUE.
                   goto 21107
                endif
                write(kou,21106)graphopt%plotmin(2)
                goto 21108
             endif
             graphopt%plotmax(2)=xxx
             graphopt%dfltmax(2)=xxx
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT RATIOS of axis, normal values 1,1 (probably quite useless....)
       case(3)
          call gparrd('X-axis plot ratio',cline,last,xxx,graphopt%xsize,q1help)
          if(xxx.le.0.1) then
             write(*,*)'Ratio set to 0.1'
             xxx=0.1D0
          endif
          graphopt%xsize=xxx
          call gparrd('Y-axis plot ratio',cline,last,xxx,graphopt%ysize,q1help)
          if(xxx.le.0.1) then
             write(*,*)'Ratio set to 0.1'
             xxx=0.1D0
          endif
          graphopt%ysize=xxx
!          write(*,*)'Not implemented yet'
          goto 21100
!-----------------------------------------------------------
! PLOT AXIS_LABELS
       case(4)
          call gparcd('For X or Y axis? ',cline,last,1,ch1,'X',q1help)
          if(ch1.eq.'X' .or. ch1.eq.'x') then
             call gparcd('Axis label: ',cline,last,5,&
                  graphopt%plotlabels(2),axplot(1),q1help)
! remember that plotlabel(1) is the title
             graphopt%labeldefaults(2)=len(graphopt%plotlabels(2))
          elseif(ch1.eq.'Y' .or. ch1.eq.'y') then
             call gparcd('Axis label: ',cline,last,5,&
                  graphopt%plotlabels(3),axplot(2),q1help)
! remember that plotlabel(1) is the title
             graphopt%labeldefaults(3)=len(graphopt%plotlabels(3))
          else
             write(kou,*)'Please answer X or Y'
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT MANIPULATE LINE COLORS
       case(5)
          write(kou,22400)
22400     format('OC uses GNUPLOT and it is possible to edit',&
               ' the file "ocgnu.plt" file'/&
               'generated by OC to use extensive facilities',&
               ' provided by GNUPLOT.'/&
               'Only a few of them is provided here.'/&
               'OC has 10 different colors to identify the lines plotted.',&
               ' Line 11 or'/' higher will repeat these colors.  With',&
               ' this command you can select'/' one of these 12 colors',&
               ' to be used for the first line plotted.')
          call gparid('The color index should be on the first line?',&
               cline,last,flc,1,q1help)
          if(flc.lt.1 .or. flc.gt.10) then
             write(*,*)'Number must be between 1 and 10'
          else
             graphopt%linett=flc
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT TITLE
       case(6)
          call gparcd('Plot title',cline,last,5,line,'DEFAULT',q1help)
          if(line(1:8).eq.'DEFAULT ') then
             graphopt%labeldefaults(1)=0
          else
             graphopt%plotlabels(1)=line
             graphopt%labeldefaults(1)=len_trim(graphopt%plotlabels(1))
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT GRAPHICS_FORMAT
! when setting graphics format always also ask for plot file
       case(7,8)
          if(kom2.eq.7) then
! subroutine TOPHLP forces return with ? in position cline(1:1)
29130        continue
             call gparid('Graphics format index:',cline,last,grunit,1,tophlp)
!             if(cline(1:1).eq.'?'
             if(cline(1:1).eq.'?' .or. &
                  grunit.lt.1.or.grunit.gt.graphopt%gnutermax) then
                write(kou,29133)
29133           format('Avalable graphics formats are:')
                write(kou,29135)(i1,graphopt%gnutermid(i1),&
                     i1=1,graphopt%gnutermax)
29135           format(i3,2x,a)
                goto 29130
             endif
             graphopt%gnutermsel=grunit
             write(kou,*)'Graphics format set to: ',graphopt%gnutermid(grunit)
          endif
!-----------------------------------------------------------
! PLOT OUTPUT_FILE, always asked when changing graphics terminal type
21140     continue
! default extension: 1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT, 8=LOG
! negative is for write, 0 read without filter, -100 write without filter
! DO NOT USE tinyfiledialog here ...
!          call gparfile('Plot file name: ',cline,last,1,plotfile,' ',-5,q1help)
          call gparcd('Plot file',cline,last,1,plotfile,'ocgnu',q1help)
          once=.false.
          if(plotfile(1:2).eq.'./') then
! save in macro directory if iumaclevl>0, else in current working directory
             if(iumaclevl.gt.0) then
                plotfile=trim(macropath(iumaclevl))//plotfile
             else
                plotfile=trim(workingdir)//plotfile
             endif
             write(*,*)'Saving on file: ',trim(plotfile)
             once=.true.
          endif
          if(plotfile(1:6).ne.'ocgnu ') then
             if(index(plotfile,'.').le.0) then
                if(graphopt%gnutermsel.ne.1) then
                   filename=trim(plotfile)//'.'//&
                        graphopt%filext(graphopt%gnutermsel)
                else
! just changing name of the GNUPLOT command file
                   filename=trim(plotfile)//'.plt '
                   plotfile=filename
                endif
             endif
!             filename=trim(plotfile)//'.plt '
             inquire(file=filename,exist=logok)
             if(logok) then
                call gparcd('File exists, overwrite?',&
                     cline,last,1,ch1,'N',q1help)
                if(.not.(ch1.eq.'Y' .or. ch1.eq.'y')) then
                   write(*,133)
                   plotfile=' '
                   goto 21140
                endif
                write(*,134)trim(filename)
                once=.true.
             endif
          endif
!          if(.not.once) write(*,'('P134)trim(filename)
! I am not sure how to inform user where the plot file is saved ....
          goto 21100
!-----------------------------------------------------------
! PLOT GIBBS_TRIANGLE
       case(9)
!          write(*,*)'Not implemented yet'
          chz='Y'
          if(graphopt%gibbstriangle) chz='N'
          call gparcd('A Gibbs triangle diagram?',cline,last,5,ch1,chz,q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
             graphopt%gibbstriangle=.TRUE.
             write(*,22500)
22500        format('The Gibbs triangle layout courtesy of',&
                  ' Catalina Pineda Heresi at RUB, Germany')
          else
             graphopt%gibbstriangle=.FALSE.
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT QUIT
       case(10)
! just return to command level
!-----------------------------------------------------------
! PLOT position of line labels (keys)
       case(11)
          write(kou,21200)
21200     format('Key to lines can be positioned: '/&
               'top/bottom left/center/right inside/outside on/off')
          call gparcd('Position?',cline,last,5,line,'top right',q1help)
          graphopt%labelkey=line
          call gparcd('Font,size: ',cline,last,5,line,'arial,12',q1help)
          graphopt%labelkey=trim(graphopt%labelkey)//' font "'//trim(line)//'"'
!          write(*,*)'pmon: ',trim(graphopt%labelkey)
          goto 21100
!-----------------------------------------------------------
! PLOT APPEND a gnuplot file
       case(12)
          write(kou,*)'Give a file name with graphics in GNUPLOT format'
! append plot file, specifying extension PLT
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
          call gparfile('File name',cline,last,1,filename,'  ',5,q1help)
!          call gparcd('File name',cline,last,1,text,'  ',q1help)
! check it is OK and add .plt if necessary ...
          jp=index(filename,'.plt ')
          if(jp.le.0) then
             jp=len_trim(filename)
             filename(jp+1:)='.plt'
          endif
          open(23,file=filename,status='old',access='sequential',err=21300)
          close(23)
          graphopt%appendfile=filename
          goto 21100
! error opening file, remove any previous appended file
21300     continue
          if(graphopt%appendfile(1:1).ne.' ') then
             write(*,21304)trim(graphopt%appendfile)
21304        format('Removing append file: ',a)
          else
             write(kou,*)'No such file name: ',trim(filename)
          endif
          graphopt%appendfile=' '
          goto 21100
!-----------------------------------------------------------
! PLOT TEXT anywhere on plot
       case(13)
          labelp=>graphopt%firsttextlabel
          if(associated(labelp)) then
             call gparcd('Modify existing text?',cline,last,1,ch1,'NO',q1help)
             if(ch1.eq.'y' .or. ch1.eq.'Y') then
                jp=0
                do while(associated(labelp))
                   jp=jp+1
                   write(kou,2310)jp,labelp%xpos,labelp%ypos,&
                        labelp%textfontscale,labelp%angle,&
                        trim(labelp%textline)
2310               format(i3,2(1pe12.4),2x,0pF5.2,2x,i4,5x,a)
                   labelp=>labelp%nexttextlabel
                enddo
                call gparid('Which text index?',cline,last,kl,1,q1help)
                if(kl.lt.1 .or. kl.gt.jp) then
                   write(*,*)'No such text label'
                   goto 21100
                endif
                labelp=>graphopt%firsttextlabel
                do jp=2,kl
                   labelp=>labelp%nexttextlabel
                enddo
                call gparcd('New text: ',cline,last,5,text,&
                     labelp%textline,q1help)
                labelp%textline=trim(text)
                call gparrd('New X position: ',cline,last,xxx,&
                     labelp%xpos,q1help)
                call gparrd('New Y position: ',cline,last,xxy,&
                     labelp%ypos,q1help)
                call gparrd('New Fontscale: ',cline,last,&
                     textfontscale,labelp%textfontscale,q1help)
                if(textfontscale.lt.0.2) textfontscale=0.2
                call gparid('New angle (degrees): ',cline,last,j1,&
                     labelp%angle,q1help)
                if(buperr.ne.0) then
                   write(*,*)'Error reading coordinates'; buperr=0; goto 21100
                endif
                labelp%xpos=xxx
                labelp%ypos=xxy
                labelp%textfontscale=textfontscale
                labelp%angle=j1
! ask for more options
                goto 21100
             endif
          endif
! input a new label
          call gparrd('X position: ',cline,last,xxx,zero,q1help)
          call gparrd('Y position: ',cline,last,xxy,zero,q1help)
          call gparrd('Fontscale: ',cline,last,textfontscale,0.8D0,q1help)
          if(textfontscale.le.0.2) textfontscale=0.2
          call gparid('Angle (degree): ',cline,last,j1,0,q1help)
          if(buperr.ne.0) then
             write(*,*)'Error reading coordinates'; buperr=0; goto 21100
          endif
          line=' '
          if(noofaxis.eq.2) then
! Calculate the equilibria at the specific point
             write(kou,22100)
22100        format(' *** Note: the positioning of the text will use the ',&
                  'axis variables for which',/11x,'the diagram was calculated',&
                  ' even if you plot with other variables!')
             call gparcd('Do you want to calculate the equilibrium? ',&
                  cline,last,1,ch1,'Y',q1help)
             if(ch1.eq.'y' .or. ch1.eq.'Y') then
! Check if plotted diagram (axplot) has same axis as calculated (axarr)??
! Or better, calculate using the plot axis ...
                line=' '
                call calc_diagram_point(axarr,axplot,xxx,xxy,line,ceq)
                if(gx%bmperr.ne.0) then
                   write(*,*)'Calculation failed ',gx%bmperr
                   gx%bmperr=0
                   line='Sorry calculation failed'
                endif
! when implemented add the stable phase names to "line" as default for text
             endif
          endif
! There is no gparcd which allows editing the existing text ... emacs!!
          text=' '
          call gparcd('Text: ',cline,last,5,text,line,q1help)
          if(text(1:1).eq.' ') then
             write(*,*)'Label ignored'
             goto 21100
          endif
! I know one should never allocate pointers but this is the only way ???
          allocate(textlabel)
          textlabel%xpos=xxx
          textlabel%ypos=xxy
          textlabel%textfontscale=textfontscale
          textlabel%angle=j1
          textlabel%textline=trim(text)
          if(associated(graphopt%firsttextlabel)) then
             textlabel%nexttextlabel=>graphopt%firsttextlabel
!             write(*,*)trim(graphopt%firsttextlabel%textline)
          else
             nullify(textlabel%nexttextlabel)
          endif
          graphopt%firsttextlabel=>textlabel
! the record is now linked from graphopt, nullify the pointer ...
! (some memory lost ...)
          nullify(textlabel)
! also clean the cline character otherwise labels may be overwritten
          cline=' '
          last=len(cline)
          goto 21100
!-----------------------------------------------------------
! PLOT TIE_LINES increment
       case(14)
          call gparid('Tie-line plot increment?',cline,last,kl,3,q1help)
          if(kl.lt.0) kl=0
          graphopt%tielines=kl
!          write(*,*)'No implemented yet'
          goto 21100
!-----------------------------------------------------------
! PLOT FONT_AND_COLOR ... and some more things ...
       case(15)
          call gparcd('Font ',cline,last,1,name1,'default',q1help)
          write(*,*)'Sorry this option not yet implemeted'
! monovariant and tielinecolor declared in smp2.F90
          call gparcd('Monovariant color ',cline,last,1,&
               name1,monovariant,q1help)
          call capson(name1)
          do kl=1,6
             if(name1(kl:kl).lt.'0' .or. name1(kl:kl).gt.'9') then
                if(name1(kl:kl).lt.'A' .or. name1(kl:kl).gt.'F') then
                   write(*,*)'The color must be a hexadecimal value',&
                        ' between 000000 (black) and FFFFFF (white)'
                   goto 21100
                endif
             endif
          enddo
          monovariant=name1(1:6)
          call gparcd('Tie-line color ',cline,last,1,&
               name1,tielinecolor,q1help)
          call capson(name1)
          do kl=1,6
             if(name1(kl:kl).lt.'0' .or. name1(kl:kl).gt.'9') then
                if(name1(kl:kl).lt.'A' .or. name1(kl:kl).gt.'F') then
                   write(*,*)'Wrong color, must be between 000000 and FFFFFF'
                   goto 21100
                endif
             endif
          enddo
          tielinecolor=name1(1:6)
          goto 21100
!-----------------------------------------------------------
! PLOT LOGSCALE
       case(16)
          call gparcd('For x or y axis? ',cline,last,1,ch1,'y',q1help)
          if(ch1.eq.'x') then
             if(graphopt%axistype(1).eq.1) then
                write(kou,*)'The x axis set to linear'
                graphopt%axistype(1)=0
             else
                graphopt%axistype(1)=1
! set range to defaults when changing to LOG 
                graphopt%rangedefaults(1)=0
             endif
          elseif(ch1.eq.'y') then
             if(graphopt%axistype(2).eq.1) then
                write(kou,*)'The y axis set to linear'
                graphopt%axistype(2)=0
             else
                graphopt%axistype(2)=1
! set range to defaults when changing to LOG 
                graphopt%rangedefaults(2)=0
             endif
          else
             write(kou,*)'Please answer x or y'
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT LINE_WITH_POINTS
       case(17)
          call gparcd('Plot a symbol at each calculated point?',&
               cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'Y' .or. ch1.eq.'y') then
             graphopt%linestyle=1
          else
             graphopt%linestyle=0
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT PAUSE_OPTIONS
       case(18)
          write(kou,*)'Specify option after pause !'
          call gparc('GNUPLOT pause option?',cline,last,5,text,' ',q1help)
          if(len_trim(text).eq.0) then
             write(kou,*)'Warning, plot will exit directly!'
!             text='-1'
          endif
          graphopt%plotend='pause '//text
          goto 21100
!-----------------------------------------------------------
! MISCELLANEOUS, for texts use option 1 in gparc to allow ,,,, as finish
       case(19)
          call gparc('Text in lower left corner?',cline,last,1,text,' ',q1help)
          graphopt%lowerleftcorner=text
          call gparcd('Spawn plot?',cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'Y') then
             graphopt%status=ibset(graphopt%status,GRKEEP)
          else
             graphopt%status=ibclr(graphopt%status,GRKEEP)
          endif
          call gparcd('Remove headings?',cline,last,1,ch1,'N',q1help)
          if(ch1.ne.'N') then
             write(*,*)'No title set!',ch1
             graphopt%status=ibset(graphopt%status,GRNOTITLE)
          else
             graphopt%status=ibclr(graphopt%status,GRNOTITLE)
          endif
          goto 21100
!-----------------------------------------------------------
! unused
       case(20)
          goto 21100
!-----------------------------------------------------------
! unused
       case(21)
          goto 21100
       end SELECT plotoption
!=================================================================
! HPCALC
    case(22)
       call hpcalc
       buperr=0
!=================================================================
! FIN, do not ask if sure, the French always know what they do ...
    case(23)
       if(logfil.gt.0) then
          write(logfil,*)'set interactive'
       endif
       call openlogfile(' ',' ',-1)
       stop 'Au revoir'
!=================================================================
! OPTIMIZE and CONTINUE.  Current optimizer is optimizers(optimizer)
    case(24)
       call gparid('Number of iterations: ',cline,last,i1,nopt1,q1help)
       if(buperr.ne.0) goto 100
       nopt1=i1
       nopt=i1
!       write(*,606)'dead 1',mexp,nvcoeff,iexit
606    format(a,10i4)
! some optimizers have no CONTINUE
!       if(optimizer.eq.1) iexit(4)=0
!       continue: if(mexp.gt.0 .and. iexit(4).eq.2) then
! iexit(4) from previous optimize allows continue with same Jacobian
!          call gparcd('Continue with same Jacobian? ',cline,last,1,&
!               ch1,'Y',q1help)
!          if(ch1.eq.'Y') then
!             ient=1
!             goto 987
!          endif
!       endif continue
! Initiate arrays when new optimization
!       ient=0
       if(.not.allocated(firstash%eqlista)) then
          write(kou,*)'There are no equilibria with experiments!'
          goto 100
       endif
!       write(*,*)'dead 2A',mexp,nvcoeff
!       if(allocated(www)) then
!          write(*,*)'Deallocating www: ',size(www),www(1)
!          deallocate(www)
!       endif
!       write(*,*)'dead 2B',mexp,nvcoeff
       if(allocated(coefs)) deallocate(coefs)
!       write(*,*)'dead 2C',mexp,nvcoeff
       if(allocated(errs))  deallocate(errs)
!       write(*,*)'dead 3',mexp,nvcoeff
! size of errors array, sum experiments for all equilibria
       mexp=0
       do i1=1,size(firstash%eqlista)
! skip equilibria with zero weight
          if(firstash%eqlista(i1)%p1%weight.eq.zero) cycle
          if(associated(firstash%eqlista(i1)%p1%lastexperiment)) then
             i2=firstash%eqlista(i1)%p1%lastexperiment%seqz
             mexp=mexp+i2
          else
             write(*,*)'No experiment in equilibrium ',i1
          endif
       enddo
!       write(*,*)'Number of experiments: ',mexp
       allocate(errs(mexp))
       updatemexp=.false.
! copy the variable coefficients to coefs
       if(nvcoeff.le.0) then
          write(*,*)'No coefficients to optimize'
          nvcoeff=0
       else
          i2=0
          allocate(coefs(nvcoeff))
          do i1=0,size(firstash%coeffstate)-1
             if(firstash%coeffstate(i1).ge.10) then
                i2=i2+1
                if(i2.gt.nvcoeff) then
                   write(kou,*)'More variable coefficients than expected',&
                        i2,nvcoeff
                   goto 100
                endif
                coefs(i2)=firstash%coeffvalues(i1)
!                coefs(i2)=firstash%coeffvalues(i1)*firstash%scale(i1)
! We do not have to bother about the associtated TP variable, it will
! be set by the calfun routine to coefs*firstashscale
!                call change_optcoeff(firstash%coeffindex(i1),&
!                     firstash%coeffvalues(i1))
!                     firstash%coeffvalues(i1))
                if(gx%bmperr.ne.0) then
                   write(*,*)'Error finding coefficient TP fun'
                   goto 100
                endif
             endif
          enddo
          if(i2.lt.nvcoeff) then
             write(kou,*)'Internal error for variable coefficients',&
                  i2,nvcoeff
             goto 100
          endif
       endif
! JUMP HERE IF CONTINUE optimization  ... NOT YET implemented
987    continue
! mexp    Number of experiments
! nvcoeff Number of coefficients to be optimized
! errs Array with differences with experiments and calculated values
! coefs Array with coefficients
       if(mexp.le.0 .or. nvcoeff.le.0) then
          write(kou,569)mexp,nvcoeff
569       format('Cannot optimize with zero experiments or coefficients',2i5)
          goto 100
       endif
       firstash%lwam=lwam
       write(*,558)mexp,nvcoeff,lwam
558    format(/'*************************************************************'/&
            '>>>   Start of optimization using LMDIF'/&
            '>>>   with ',i4,' experiments and ',i3,' coefficients ',/&
            '>>>   and allocated workspace ',i5/&
            '*************************************************************')
!
       j1=nopt
       if(.not.allocated(iwam)) then
! value of lwam set by user
          allocate(iwam(lwam))
          allocate(wam(lwam))
       endif
       if(allocated(fjac)) deallocate(fjac)
! fjac is used to calculate the Jacobian and other things
! err0(1) is set to the sum of errors squared for the initial values of coefs
573    format(a,6(1pe12.4))
       allocate(fjac(mexp,nvcoeff))
!       write(*,'(a,10(1pe12.4))')'lmdif1: ',(coefs(iz),iz=1,nvcoeff)
!->->->->->-> HERE THE OPTIMIZATION IS MADE <-<-<-<-<-<-
! nfev set to number of iterations
       call lmdif1(mexp,nvcoeff,coefs,errs,optacc,nopt,nfev,&
            iwam,wam,lwam,fjac,err0)
!->->->->->-> HERE THE OPTIMIZATION IS MADE <-<-<-<-<-<-
       mexpdone=mexp
       nvcoeffdone=nvcoeff
! on return nopt is set to a message but 
! first copy the coefs to coeffvalues ...
!       write(*,573)'Coeffs ut: ',(coefs(j2),j2=1,nvcoeff)
       i2=0
       do i1=0,size(firstash%coeffstate)-1
          if(firstash%coeffstate(i1).ge.10) then
             i2=i2+1
             firstash%coeffvalues(i1)=coefs(i2)
!             write(*,555)'final: ',i1,i2,&
!                  firstash%coeffvalues(i1)*firstash%coeffscale(i1),&
!                  coefs(i2),firstash%coeffscale(i1)
!555          format(a,2i3,3(1pe12.4))
          endif
       enddo
! then calculate final sum of errots
       xxx=zero
       do i2=1,mexp
          xxx=xxx+errs(i2)**2
       enddo
! this is the final sum of errors squared
       err0(2)=xxx
       if(mexp-nvcoeff.gt.0) then
! should I add or subract 1??
          err0(3)=xxx/real(mexp-nvcoeff)
       else
! when equal number of experiment and coefficients
          err0(3)=1.0D30
       endif
! The top nvcoeff*nvcoeff submatrix of fjac is R^T * R
!       write(*,*)'The unsymmetric R^T*R submatrix returned from lmfif1:'
!       do i2=1,nvcoeff
!          write(*,563)(fjac(j2,i2),j2=1,nvcoeff)
!       enddo
!       read(*,'(a)')ch1
! cormat will be the CORRELATION MATRIX if optimization successful
! otherwise it will not be allocated
       if(allocated(cormat)) then
          deallocate(cormat)
          deallocate(tccovar)
       endif
!--------------- begin calculate correlation matrix and RSD
! zero the relative standard deviations (RSD)
       firstash%coeffrsd=zero
       if(j1.gt.0 .and. nopt.gt.0 .and. nopt.le.6) then
! if there is a result calculate the Jacobian in fjac
! mexp,nvcoeff,coeffs,errs are same as in the call to lmdif1
! This will overwrite the fjac returned from the call to lmdif1
!          write(*,*)'Calculating the Jacobian: '
! allocate array to extract calculated values of experiments
          if(allocated(calcexp)) deallocate(calcexp)
          allocate(calcexp(mexp))
          iflag=2
! penulitima argument zero means use machine precision to calculate derivative
          call fdjac2(mexp,nvcoeff,coefs,errs,fjac,mexp,iflag,zero,wam)
! debug output ...
!          write(*,*)'pmon: fjac: ',nvcoeff,mexp,iflag
!          do i2=1,mexp
!             write(*,563)(fjac(i2,ll),ll=1,nvcoeff)
!          enddo
563       format(6(1pe12.4))
!          write(*,*)'End listing of Jacobian fjac calculated by fdjac2'
!          read(*,'(a)')ch1
! Next calculate M = (fjac)^T (fjac); ( ^T means transponat)
          if(allocated(cov1)) deallocate(cov1)
! the cov1 is symmetric and should have these dimensions:
          allocate(cov1(nvcoeff,nvcoeff))
          cov1=zero
          do i2=1,nvcoeff
             do j2=1,nvcoeff
                xxx=zero
                do ll=1,mexp
                   xxx=xxx+fjac(ll,i2)*fjac(ll,j2)
!                   write(*,564)'xxx: ',i2,j2,ll,xxx
564                format(a,3i5,1pe12.4)
                enddo
! this matrix is symmetric ... which index first ???
                cov1(j2,i2)=xxx
!                cov1(i2,j2)=xxx
             enddo
          enddo
!          write(*,*)'M = (Jac)^T (Jac); (^T means transponat)',nvcoeff
!          do i2=1,nvcoeff
!             write(*,563)(cov1(i2,ll),ll=1,nvcoeff)
!          enddo
! invert cov1 using LAPACK+BLAS via Lukas routine ...
          if(nvcoeff.gt.1) then
! cormat deallocated above, dimension is cormat(nvcoeff,nvcoeff) !!
             allocate(cormat(nvcoeff,nvcoeff))
             allocate(tccovar(nvcoeff,nvcoeff))
! symmetric?   call mdinv(nvcoeff,nvcoeff+1,cov1,cormat,nvcoeff,iflag)
! NOTE: cov1 and cormat should both have dimension cov1(nvcoeff,nvcoeff)
             call mdinv(nvcoeff,cov1,cormat,nvcoeff,iflag)
! invert unsymmetrical matrix
             if(iflag.eq.0) then
                write(*,*)'Failed invert matrix=Jac^T*Jac',iflag
             endif
! RSD depend on scaling factor of coefficient
!             write(*,*)'PMON norm.error and covariant matrix: ',err0(3)
!             do i1=1,nvcoeff
!                write(*,'(6(1pe12.4))')(cormat(i1,i2),i2=1,nvcoeff)
!             enddo
! all elements in the covariance matrix should be multiplied with err0(3)
             tccovar=cormat
             do i1=1,nvcoeff
                do i2=1,nvcoeff
! I get exactly the same RSD as TC if I ignore the normalized error !!
! but according to theory it should be multiplied with the normalized error
                  cormat(i1,i2)=err0(3)*cormat(i1,i2)
                enddo
             enddo
! divide all values with the square root of the  diagonal elements
! save covarance matrix n cov1
             cov1=cormat
             do i1=1,nvcoeff
                do i2=1,nvcoeff
                   xxx=sqrt(abs(cov1(i1,i1)*cov1(i2,i2)))
                   cormat(i1,i2)=cormat(i1,i2)/xxx
                enddo
             enddo
!             write(*,*)'Correlation after dividing with sqrt(abs(c_ii*c_jj))'
!             do i1=1,nvcoeff
!                write(*,'(6(1pe12.4))')(cormat(i1,i2),i2=1,nvcoeff)
!             enddo
          elseif(abs(cov1(1,1)).gt.1.0D-38) then
! cov1 is just a single value
             allocate(cormat(1,1))
             allocate(tccovar(1,1))
!             cormat(1,1)=one
! IF THERE IS A SINGLE VARIABLE ITS CORRELATION MATRIX MUST BE UNITY
             cormat(1,1)=one
             tccovar(1,1)=one
          else
             write(*,*)'Correlation matrix singular'
          endif
       endif
! write the correlation matrix  this is still very uncertain ,,,
!       if(allocated(cormat)) then
!          if(nvcoeff.gt.0) then
!             write(*,*)'Correlation matrix (symmetric):'
!             do i2=1,nvcoeff
!                write(kou,'(8(1pe10.2))')(cormat(i2,j2),j2=1,nvcoeff)
!             enddo
!          endif
!       endif
! zero all RSD values
       firstash%coeffrsd=zero
       if(allocated(cormat) .and. allocated(cov1)) then
! calculate the RSD (Relative Standard Deviation) for each parameter
! the last calculated values of the experiments in calcexp
!          write(*,*)'The sum of all calculated equilibria,',&
!               ' very different magnitudes ...'
          xxx=zero
          do i2=1,mexp
! the calculated value is stored in calcexp by fdjac if calcexp is allocated!
             xxx=xxx+calcexp(i2)
!             write(*,766)i2,calcexp(i2),xxx
766          format('pmon: Calculated value',i4,2(1pe12.4))
          enddo
!          ll=max(1,mexp-nvcoeff)
! This value may be negative!
!          xxy=xxx/real(ll)
! the difference between the calculated and experimental value is errs(1:mexp)
! err0(2) is sum of all errors squared          
!          xxx=err0(2)/real(ll)  ... this is err0(3)
! I am not sure about this ...
          i2=0
          do i1=0,size(firstash%coeffstate)-1
             if(firstash%coeffstate(i1).ge.10) then
! this is an optimized parameter, they are indexed starting from zero!!
                i2=i2+1
! But in cormat they are indexed from 1 .. nvcoeff
!                firstash%coeffrsd(i1)=sqrt(abs(cormat(i2,i2))*xxx)/xxy
!                write(*,'(a,3(1pe12.4))')'RSD: ',cov1(i2,i2),xxx,xxy
!                firstash%coeffrsd(i1)=abs(sqrt(abs(cov1(i2,i2))*xxx)/xxy)
! we have already multiplied all terms in covariance matrix with err0(3)
!                firstash%coeffrsd(i1)=abs(sqrt(abs(cov1(i2,i2))*err0(3)
                firstash%coeffrsd(i1)=abs(sqrt(abs(cov1(i2,i2))))
             endif
          enddo
       endif
! deallocate calcexp to avoid storing these values when running LMDIF
       if(allocated(calcexp)) deallocate(calcexp)
!--------------- end calculate correlation matrix and RSD
! some nice output .....
       write(kou,5020)
       if(j1.eq.0) then
          write(*,*)'Dry run with zero iterations'
       elseif(nopt.eq.0) then
          write(kou,5000)nopt
5000      format(/'*** No optimization due to improper input parameters',i3)
       elseif(nopt.eq.1) then
          write(kou,5001)nopt,optacc
5001      format(/'LMDIF return code ',i2/&
               'Relative error for sum of squares is within ',1pe10.2)
       elseif(nopt.eq.2) then
          write(kou,5002)nopt,optacc
5002      format(/'LMDIF return code ',i2/&
               'Relative error of parameters is within ',1pe10.2)
       elseif(nopt.eq.3) then
          write(kou,5003)nopt
5003      format(/'LMDIF return code ',i2,': successful optimization')
       elseif(nopt.eq.4) then
          write(kou,5004)nopt
5004      format(/'*** LMDIF return code ',i2/&
               'Sum of squares does not decrease')
       elseif(nopt.eq.5) then
          write(kou,5005)nopt,nfev
5005      format(/'*** LMDIF return code ',i2/&
               'Maximum calls of function ',i5,' exceeded')
       elseif(nopt.eq.6) then
          write(kou,5006)nopt,optacc
5006      format(/'*** LMDIF return code ',i2/&
               'Cannot reduce error, requested accuracy ',1pe10.2,' too small')
! '*** Cannot reduce error, requested accuracy 123456789. too small
       elseif(nopt.eq.6) then
          write(kou,5007)nopt,optacc
5007      format(/'*** LMDIF return code ',i2/&
              'Cannot improve result, requested accuracy ',1pe10.2,' too small')
       else
          write(kou,5008)nopt
5008      format('*** LMDIF return code ',i7/&
               'Unknown code, see LMDIF documentation.')
       endif
       write(kou,5010)nfev,err0
5010   format(/'Iterations ',i4,', sum of errors changed from ',&
            1pe14.6,' to ',1pe14.6/17x,'Normalized sum of errors:',20x,1pe14.6)
       write(kou,5020)
5020   format(/78('*'))
! finally list the coefficient values
       call listoptcoeff(mexp,err0,.FALSE.,lut)
! end of call to LMDIF
!=================================================================
! SHOW is immpemented as a special case of LIST STATE_VARIABLES
!    CASE(25)
!       write(kou,*)'Not implemented yet'
!=================================================================
! ANALYZE_ASSESSMT
    CASE(26)
       if(.not.allocated(firstash%eqlista)) then
          write(kou,*)'No assessment record'
          goto 100
       elseif(nvcoeff.le.0) then
          write(kou,*)'No variable optimizing coefficients'; goto 100
       elseif(nvcoeff.ne.nvcoeffdone) then
          write(kou,*)'No optimization made with these coefficients',&
               nvcoeff,nvcoeffdone
          goto 100 
       elseif(mexp.ne.mexpdone) then
          write(kou,*)'No optimization made with these experiments',&
               mexp,mexpdone
          goto 100 
       endif
       call gpari('Index of coefficent to change: ',cline,last,&
            analyze,NONE,q1help)
       if(buperr.ne.0) goto 990
       xxy=zero
       if(analyze.lt.0) then
! using give a negative coefficient, restore saved coefficients
! if nvcoefdone and mxexp same
          write(*,*)'Trying to restore saved coefficients'
          if(nvcoeffsave.eq.nvcoeff .and. mexpsave.eq.mexp) then
             if(allocated(savedcoeff)) then
! if analyze < 0 then restore sevedcoeff
                i2=0
                do j2=0,size(firstash%coeffstate)-1
                   if(firstash%coeffstate(j2).ge.10) then
! this a variable coefficient
                      i2=i2+1
                      firstash%coeffscale(j2)=savedcoeff(1,i2)
                      firstash%coeffstart(j2)=savedcoeff(2,i2)
! I am not sure if xxx should be savedcoeff or scale*start ... ???
                      xxx=savedcoeff(3,i2)
                      firstash%coeffvalues(j2)=xxx
                      firstash%coeffrsd(j2)=zero
! this should update all other places including TP function 
                      call change_optcoeff(firstash%coeffindex(j2),xxx)
                   endif
                enddo
                deallocate(savedcoeff)
                err0(2)=savesumerr
                write(*,*)'Restored saved coefficients'
             else
                write(*,*)'No coefficients saved'
             endif
          else
! giving a negative number makes it possible to use ANALYZE again
! for another set of coefficients and experiments
             write(kou,*)'Cannot restore as variable coefficients ',&
                  'or experiments changed'
             if(allocated(savedcoeff)) deallocate(savedcoeff)
          endif
          goto 100
       else
          if(.not.allocated(savedcoeff)) then
! when ANALYZE first time save the current set of variable coefficients
             allocate(savedcoeff(3,nvcoeff))
             mexpsave=0
! if already allocated mexpsave nonzero
          endif
          i2=0
          xxy=zero
          do j2=0,size(firstash%coeffstate)-1
! only active coefficients saved ... extract the one to be changed
             if(firstash%coeffstate(j2).ge.10) then
                i2=i2+1
                if(mexpsave.eq.0) then
                   savedcoeff(1,i2)=firstash%coeffscale(j2)
                   savedcoeff(2,i2)=firstash%coeffstart(j2)
                   savedcoeff(3,i2)=firstash%coeffvalues(j2)
!                 write(*,'(a,3(1pe14.6))')'saved: ',(savedcoeff(iz,i2),iz=1,3)
                   firstash%coeffrsd(j2)=zero
                endif
                if(analyze.eq.j2) then
                   cormatix=i2
                   xxy=savedcoeff(1,i2)*savedcoeff(3,i2)
!                   write(*,*)'Coefficient: ',cormatix,xxy
                endif
             endif
          enddo
          if(mexpsave.eq.0) then
             write(*,*)'Saved ',i2,'currently variable coefficients'
! save current sum of errors, nvcoeff and mexp
             savesumerr=err0(2)
             nvcoeffsave=nvcoeff; mexpsave=mexp
          endif
       endif
! if xxy is zero it is not an optimized coefficient
       if(xxy.eq.zero) then
          write(kou,*)'Specified coefficent not set as variable',analyze
          if(allocated(savedcoeff)) deallocate(savedcoeff); goto 100
       endif
! ask for new value with the current value as default
       call gparrd('New value: ',cline,last,xxx,xxy,q1help)
       delta=(xxx-xxy)/firstash%coeffscale(analyze)
!       write(*,*)'Delta: ',xxx-xxy,delta
! UNFINISHED
! Now all variable coefficients should be modified using the correlation matrix
       i2=0
       do j2=0,size(firstash%coeffstate)-1
! modify all other coefficient according to the correlation matrix       
! new_value_i =  old_value_i + correlation_matrix_ji * delta (where j=analyze)
          if(firstash%coeffstate(j2).ge.10) then
             i2=i2+1
             xxx=firstash%coeffvalues(j2)
             xxy=xxx+cormat(cormatix,i2)*delta
!             firstash%coeffvalues(j2)=xxy*firstash%coeffscale(j2)
! %coeffvalues should be of the order 1
! No change of %coeffstart and %coeffscale
             firstash%coeffvalues(j2)=xxy
             xxz=xxy*firstash%coeffscale(j2)
! optimizing coefficients are also TP functions, we must update the
! TP function value!! I do not understand this and "list tp" is wrong
! but it seems to work.  If I set the value *firstash%coeffscale it blows up!
             call change_optcoeff(firstash%coeffindex(j2),xxz)
!             call change_optcoeff(firstash%coeffindex(j2),xxy)
!        call change_optcoeff(firstash%coeffindex(j2),firstash%coeffvalues(j2))
!             write(*,'(a,2i4,4(1pe12.4))')'New value: ',i2,j2,&
!                  xxx,cormat(cormatix,i2),delta,firstash%coeffvalues(j2)
          endif
       enddo
       write(*,*)'To calculate a new set of errors use OPTIMIZE'
!       write(kou,*)'Not implemented yet'
!=================================================================
! unused
    CASE(27)
       write(kou,*)'Not implemented yet'
!=================================================================
! unused
    CASE(28)
       write(kou,*)'Not implemented yet'
!=================================================================
! unused
    CASE(29)
       write(kou,*)'Not implemented yet'
!=================================================================
! unused
    CASE(30)
       write(kou,*)'Not implemented yet'
!=================================================================
!
    END SELECT main
! command executed, prompt for another command unless error code
    if(gx%bmperr.eq.0) goto 100
!============================================================
! handling errors
990 continue
    write(kou,991)gx%bmperr,buperr
991 format(/'Error codes: ',2i6)
    if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
       write(kou,992)trim(bmperrmess(gx%bmperr))
992    format('Message: ',a/)
    else
       write(kou,*)'No defined error message, maybe I/O error'
    endif
    if(stop_on_error) then
! turn off macro but remain inside software
       call macend(cline,last,logok)  
       write(kou,*)'Stop_on_error set, press return to finish program'
       read(kiu,17)ch1
17     format(a)
       stop
    endif
    gx%bmperr=0; buperr=0
    goto 100
!
  end subroutine oc_command_monitor

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\bergin{verbatim}
  integer function submenu(query,cline,last,ccomnd,ncomnd,kdef)
! general subcommand decoder
! query is the prompt
! cline and last is user input and position
! ccomnd is the menu and ncomnd number of menu entries
! kdef is the default (to be added to query)
!  implicit double precision (a-h,o-z)
    implicit none
    character cline*(*),ccomnd(*)*(*),query*(*)
    integer last,kdef,ncomnd
!\end{verbatim}
!    external q2help
    character defansw*16,query1*64,text*256
    integer kom2,lend,lenq
    logical once
    lenq=len_trim(query)
    if(query(lenq:lenq).eq.'?') then
       query1=query(1:lenq)
    else
       query1=query(1:lenq)//' what?'
       lenq=lenq+6
    endif
    once=.true.
! this is to force loading of q2help on MacOS (did not help)
    submenu=0
! if cline(last:last) is "," skip one character
!  write(kou,10)'submenu 1: ',query(1:lenq),last,cline(last:last+5)
!10  format(a,a,i4,': ',a)
    if(last.lt.len(cline)) then
       if(cline(last:last).eq.',') last=last+1
    endif
100 continue
    if(kdef.lt.1 .or. kdef.gt.ncomnd) then
! no default answer
       if(eolch(cline,last)) then
! empty line, note fourth argument 5 copes whole of cline into text
!        call gparc(query1(1:lenq),cline,last,5,text,' ',tophlp)
          call gparc(query1(1:lenq),cline,last,5,text,' ',q2help)
          if(buperr.ne.0) goto 1000
          cline=text
       else
          cline=cline(last:)
       endif
    else
! there is a default answer
! this eolch skips spaces.  If only spaces it returns TRUE
       if(eolch(cline,last)) then
! there is no user input passed to this subroutine, write the question
          defansw=ccomnd(kdef)
          lend=len_trim(defansw)+1
! note fourth argument 5 copes whole of cline into text
!          write(*,102)'submenu 6: ',last,trim(cline)
          call gparcd(query1(1:lenq),cline,last,5,text,defansw,q2help)
          if(buperr.ne.0) goto 1000
          cline=text
       else
! if first character is , take default answer
!          write(*,102)'submenu 7: ',last,trim(cline)
102       format(a,i5,'"',a,'"')
          if(cline(last:last).eq.',') then
! a , means accept default answer
             submenu=kdef
             goto 1000
          else
             defansw=ccomnd(kdef)
             lend=len_trim(defansw)+1
! note fourth argument 5 copes whole of cline into text
! gparcd skips one character, backspace last, it does not matter if it is ,
             last=last-1
!             write(*,102)'submenu 8: ',last,trim(cline)
             call gparcd(query1(1:lenq),cline,last,5,text,defansw,q2help)
             if(buperr.ne.0) goto 1000
             cline=text
!             cline=cline(last:)
!             write(*,*)'sumbemu: ',trim(cline),last
! added 20190207 because "enter gamma ac(a)/x(a); gave segmentation fault
! but that was not the error, the error was missing =
!             once=.false.
          endif
       endif
    endif
!
!    write(*,102)'submenu 9: ',last,trim(cline)
    kom2=ncomp(cline,ccomnd,ncomnd,last)
    if(kom2.le.0) then
       if(once) then
          if(cline(1:1).ne.'?') once=.false.
          if(kom2.lt.0) write(kou,*)'Ambiguous answer, please try again'
          write(kou,*)'Possible answers are:'
          last=1
          cline=' *'
!        call nghelp(cline,last,ccomnd,ncomnd)
          call q3help(cline,last,ccomnd,ncomnd)
          last=len(cline)
          goto 100
       else
          write(kou,*)'Answer not understood, returning to upper level'
          goto 1000
       endif
    else
       submenu=kom2
       if(helprec%level.lt.maxhelplevel) then
          helprec%level=helprec%level+1
          helprec%cpath(helprec%level)=ccomnd(kom2)
       else
          write(*,*)'Warning, exceeded helprec%level limit 2'
       endif
    endif
!    write(*,102)'submenu last: ',last,trim(cline)
1000 continue
    return
  end function submenu

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ocmon_set_options(useroption,afo,optionsset)
    implicit none
    character*(*) useroption
    integer afo
    TYPE(ocoptions) :: optionsset
!\end{verbatim}
    integer next,kom,slen,errno,jj
    character option*64,string*64,dummy*128,date*8,time*10
    integer, parameter :: nopt=9
    character (len=16), dimension(nopt) :: copt=&
        ['OUTPUT          ','ALL             ','FORCE           ',&
         'VERBOSE         ','SILENT          ','APPEND          ',&
         '                ','                ','                ']
! copy "option" to a local string as it may be just a single character!!
    option=' '
    option=useroption
! /? will list options
    afo=0
    if(option(1:2).eq.'? ') then
       write(kou,10)
10     format('Available options (preceded by /) are:')
       next=1
       dummy=' * '
       call q3help(dummy,next,copt,nopt)
!       write(*,*)'Back from q3help'
       afo=1
       goto 1000
    endif
    kom=ncomp(option,copt,nopt,next)
    if(kom.le.0) then
       write(kou,*)'Unknown option ignored: ',option(1:len_trim(option))
       goto 1000
    else
       select case(kom)
       case default
          write(kou,*)'Option not implemented: ',option(1:len_trim(option))
          write(kou,10)
          next=1
          dummy=' * '
          call q3help(dummy,next,copt,nopt)
          afo=1
!-----------------------------------
       case(1) ! /output means open a file and ovewrite any previous content
!          write(*,*)'Option not implemented: ',option(1:len_trim(option))
! next argument after = must be a file name
! 6 means extension DAT
 !         jj=next+1
 !         if(eolch(option,jj)) then
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
             call gparfile('Output file',option,next,1,string,'  ',-7,q1help)
             if(string(1:1).eq.' ') then
                string='ocoutput.DAT'
                write(kou,*)' *** No file name given, will use: ',trim(string)
             endif
             slen=len_trim(string)
!          else
!             call getext(option,next,2,string,' ',slen)
!          endif
! add extention .dat if to extenstion provided
          if(index(string,'.').le.0) then
             string(slen+1:)='.DAT '
          endif
! close any previous output file          
          close(21)
          open(21,file=string,access='sequential',status='unknown',&
               err=900, iostat=errno)
          optionsset%lut=21
! write a header
          call date_and_time(date,time)
232       format(/'%%%%%%%%%% OC output ',a,a4,'-',a2,'-',a2,2x,a2,'h',a2)
          write(21,232)'written: ',date(1:4),date(5:6),date(7:8),&
               time(1:2),time(3:4)
          write(kou,231)'Output',trim(string)
!-----------------------------------
       case(2) ! /all ??
          write(*,*)'Option not implemented: ',trim(option)
!-----------------------------------
       case(3) ! /force
          write(*,*)'Option not implemented: ',trim(option)
!-----------------------------------
       case(4) ! /verbose
          globaldata%status=ibset(globaldata%status,GSVERBOSE)
          write(kou,*)'VERBOSE option set  ... but not really implemented'
!-----------------------------------
       case(5) ! /silent
          globaldata%status=ibclr(globaldata%status,GSVERBOSE)
          globaldata%status=ibset(globaldata%status,GSSILENT)
!-----------------------------------
       case(6) ! /APPEND, open file and write at end
!          write(*,*)'Option not implemented: ',option(1:len_trim(option))
! next argument after = must be a file name
!          jj=next
!          if(eolch(option,jj)) then
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=PDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
          call gparfile('Append to file:',option,next,&
               1,string,'  ',-7,q1help)
          if(string(1:1).eq.' ') then
             string='ocappend.DAT'
             write(kou,*)' *** No file name given, will use: ',trim(string)
          endif
!          else
!             call getext(option,next,2,string,' ',slen)
!          endif
! add extention .dat if to extension provided
          slen=len_trim(string)
          if(index(string,'.').le.0) then
             string(slen+1:)='.DAT '
          endif
! close any previous output file (should not be necessary)
          close(21)
          open(21,file=string,access='sequential',status='unknown',&
               err=900, iostat=errno)
          optionsset%lut=21
! read until end-of-file
200       continue
          read(21,210,end=220)dummy
210       format(a)
          goto 200
! write not allowed after finding EOF, we must backspace
220       continue
          backspace(21)
! write a header
          call date_and_time(date,time)
          write(21,232)'appended: ',date(1:4),date(5:6),date(7:8),&
               time(1:2),time(3:4)
          write(kou,231)'Append',trim(string)
231       format(a,' on file: ',a)
!-----------------------------------
       case(7) ! 
          continue
!-----------------------------------
       case(8) ! 
          continue
!-----------------------------------
       case(9) ! 
          continue
       end select
    endif
    goto 1000
! errors
900 continue
    write(kou,*)'Failed to open output file, error cofe=',errno
    goto 1000
1000 continue
    return
  end subroutine ocmon_set_options
    
!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ocmon_reset_options(optionsset)
    implicit none
    TYPE(ocoptions) :: optionsset
!    TYPE(ocoptions), pointer :: optionsset
!\end{verbatim}
    if(btest(globaldata%status,GSVERBOSE)) then
! reset verbose option
       if(.not.btest(globaldata%status,GSSETVERB)) then
! if user has SET VERBOSE do not resest VERBOSE
          globaldata%status=ibclr(globaldata%status,GSVERBOSE)
       endif
    endif
! reset output unit option
    if(optionsset%lut.ne.kou) then
       close(optionsset%lut)
       optionsset%lut=kou
       write(kou,"(a,i4)")'Output unit reset to screen: ',kou
    endif
!1000 continue
    return
  end subroutine ocmon_reset_options

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine listoptshort(lut,mexp,nvcoeff,errs)
! short listing of optimizing variables and result
    integer lut,mexp,nvcoeff
    double precision errs(*)
!    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: neweq
    integer i1,i2,j1,j2,j3,neq
    character name1*24,line*80
    double precision xxx,sum
    type(gtp_condition), pointer :: experiment
!
! list all experiments, only possible if there are experiments
    if(mexp.eq.0) then
       write(lut,666)
666    format(/'No experiments so no results'/)
       goto 1000
    endif
! list experiments, mexp is number of EXPERIMENTS, not equilibria!!
    write(lut,620)size(firstash%eqlista),mexp
620 format(/'List of ',i5,' equilibria with ',i5,&
         ' experimental data values'/&
         '  No Equil name      Weight Experiment $ calculated',19x,&
         'Error')
    j3=0
    allequil: do i1=1,size(firstash%eqlista)
! skip equilibria with zero weight
       neweq=>firstash%eqlista(i1)%p1
       if(neweq%weight.eq.zero) cycle allequil
       name1=neweq%eqname(1:12)
! LOOP for all experiments for this equilibrium (maybe none??)
       if(.not.associated(neweq%lastexperiment)) cycle allequil
       experiment=>neweq%lastexperiment%next
       if(.not.associated(experiment)) cycle allequil
700    continue
          i2=neweq%lastexperiment%seqz
!          write(*,*)'number of experiments: ',i2
          neq=neweq%eqno
          do j2=1,i2
! j1 is position in line to write experiment
             j1=1
             line=' '
! this subroutine returns experiment and calculated value: "H=1000:200 $ 5000"
             call meq_get_one_experiment(j1,line,j2,neweq)
             j3=j3+1
             if(neq.gt.0) then
!                write(lut,622)neq,name1(1:12),neweq%weight,line(1:40),errs(j3)
                write(lut,622)neq,name1(1:15),neweq%weight,line(1:40),errs(j3)
622             format(i4,1x,a,2x,F5.2,1x,a,1x,1pe12.4)
                neq=0
             else
                write(lut,623)line(1:40),errs(j3)
623             format(28x,a,1x,1pe12.4)
             endif
! list the equilibrium name just for the first (or only) experiment
!             name1=' '
          enddo
          experiment=>experiment%next
!590       if(.not.associated(experiment,neweq%lastexperiment)) then
!             experiment=>experiment%next
!             goto 700
          if(j2.lt.i2 .and. .not.associated(experiment)) then
             write(*,*)'Missing experiment in equilibrium ',neweq%eqno
             cycle allequil
          endif
    enddo allequil
! list sum of squares
    sum=zero
    do j1=1,mexp
       sum=sum+errs(j1)**2
    enddo
! same as PARROT
    j1=mexp-nvcoeff
    if(j1.gt.0) then
       write(lut,621)sum,mexp,nvcoeff,j1,sum/j1
    else
       write(lut,621)sum,mexp,nvcoeff,0,zero
    endif
621 format(/'Final sum of squared errors: ',1pe16.5,&
         ' using ',i4,' experiments and'/&
         i3,' coefficient(s).  Degrees of freedom: ',i4,&
         ', normalized error: ',1pe13.4/)
1000 continue
    return
  end subroutine listoptshort

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine calctrans(cline,last,ceq)
! calculate a phase transition
    character cline*(*)
    integer last
    type(gtp_equilibrium_data), pointer :: ceq
!\end
    character name1*30
    integer j1,iph,ics
    double precision xxx
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: stvr
!
    write(kou,2090)
2090 format('To calculate when a phase will appear/disappear',&
          ' by releasing a condition.')
    if(btest(ceq%status,EQNOEQCAL)) then
       write(kou,2095)
2095   format('You must make an equilibrium calculation before using',&
            ' this command.')
       goto 1000
    endif
    call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
    call find_phase_by_name(name1,iph,ics)
    if(gx%bmperr.ne.0) goto 1000
    j1=test_phase_status(iph,ics,xxx,ceq)
    if(j1.eq.PHFIXED) then
       write(kou,*)'Phase status already fixed'
       goto 1000
    endif
    call list_conditions(kou,ceq)
    write(kou,2097)
2097 format('You must release one condition, give its number')
    call gparid('Condition number',cline,last,j1,1,q1help)
    if(j1.le.0 .or. j1.gt.noel()+2) then
       write(kou,*)'No such condition'
       goto 1000
    endif
! this finds condition with given number
    call locate_condition(j1,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(pcond%active.eq.0) then
! the condition is active, deactivate it!
       pcond%active=1
    else
       write(kou,*)'This condition is not active!'
       goto 1000
    endif
! Condition released, now set the phase as fix with zero moles
    call change_phase_status(iph,ics,PHFIXED,xxx,ceq)
    if(gx%bmperr.ne.0) goto 1000
! Calculate equilibrium
    call calceq2(1,ceq)
    if(gx%bmperr.ne.0) goto 1000
! get the value of the released condition and set it to the new value
    stvr=>pcond%statvar(1)
    call state_variable_val(stvr,xxx,ceq)
    if(gx%bmperr.ne.0) goto 1000
    write(kou,2099)xxx
2099 format('The transition occurs at ',1pe16.8,', set as condition')
    pcond%prescribed=xxx
    pcond%active=0
! set phase back as entered and stable
!    write(*,*)'Set phase back as entered'
    call change_phase_status(iph,ics,PHENTSTAB,zero,ceq)
1000 continue
    return
  end subroutine calctrans

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

END MODULE cmon1oc

