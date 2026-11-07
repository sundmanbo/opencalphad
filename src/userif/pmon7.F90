MODULE cmon1oc
!
! Copyright 2012-2026, Bo Sundman, France
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
!*****************************
! command line monitor for OC 
!*****************************
!
  use ocsmp
  use liboceqplus
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
!------------------------------------------------------------------
! Data kept between commands, moved here from oc_command_monitor
! so they can be used also in subroutines called from the monitor
!
! variables for optimization (lmdif)
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
  integer analyze,cormatix,nvcoeffsave,mexpsave
! this is least square error from using LMDIF
! 1: previous value, 2 new value, 3 normalized error (divided by m-n)
  double precision err0(3)
  double precision, dimension(:), allocatable :: errs
  double precision, dimension(:), allocatable :: coefs
!
! variables for step, map and plot
! graphics record for plot ranges, texts and defaults
  type(graphics_options) :: graphopt
! plot texts
  type(graphics_textlabel), pointer :: textlabel
! axis data structures
  type(map_axis), dimension(5) :: axarr
  integer noofaxis
! for map results
  type(map_node), pointer :: maptop,maptopsave
! seqxyz has initial values of seqx, seqy and seqz
  integer seqxyz(3)
! default plot axis for some STEP command:
! 1 for SEPARATE, 2 SCHEIL, 3 TZERO, 4 PARAEQUIL, 5 LIQUID_EET
  logical stepspecial(5)
  character plotfile*256,axplot(2)*24,axplotdef(2)*24
! plot unit for experimental data used in enter many_equilibria
  integer :: plotdataunit(9)=0,plotunit0=0
!------------------------------------------------------------------
!
contains
!
!\begin{verbatim}
  subroutine oc_command_monitor(version,linkdate,narg,argline)
! command monitor
    implicit none
!
! linkdat is date when program was linked
! argline and narg are inline arguments
    character linkdate*(*),version*(*),argline(*)*(*)
    integer narg
!\end{verbatim}
! various symbols and texts, version 7
    character :: ocprompt*8='--->OC7:'
    character name1*24,name2*24,name3*24,dummy*24,line*80,model*72,chshort*1
    integer, parameter :: ocmonversion=75
! for the on-line help, at present turn off by default, if a HTML file set TRUE
    character*128 browser,latexfile,htmlfile,unformfile,xtdbdef
    logical :: htmlhelp=.FALSE.
! element symbol and array of element symbols for database use
    character elsym*2,ellist(maxel)*2,elbase(maxel)*2,parael*2
! more texts for various purposes
    character text*72,string*256,ch1*1,chz*1,selection*27,funstring*1024
    character quest*20
    character longstring*5000,optres*40
! measure calculate carefully
    double precision finish2,start2
    integer endoftime,startoftime
! separate file names for remembering and providing a default
    character ocmfile*128,ocufile*128,tdbfile*128,xtdbfile*128
    character ocdfile*128,filename*128
    character zext*8,mqmqass*60
! home for OC and default directory for databases
    character ochome*128,ocbase*128
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
! the beginning of a sequential list of all ternary methods
! cpu time measurements
    double precision ending,starting
!>>>> has to be reorganized ------------------------------------
! axis variables and limits
! default values used for axis variables
    double precision dinc,dmin,dmax
    integer grunit
! species for ternary extrapolation method
    character xspecies(3)*24,tkmode*6
! path to start directory declared inside metlib!!
! plot texts
    type(graphics_textlabel), pointer :: labelp
    type(map_node), pointer :: mapnode
! csv file converion
    integer ioc,ip
! this should be removed
! for paraequilibrium meqrec is needed also here
    TYPE(meq_setup), pointer :: meqrec
    TYPE(meq_setup), allocatable, target :: meqrec1
!<<<<<<<--------------------------------------------------------------
! used for element data and R*T
    double precision h298,s298,rgast
! temporary reals
    double precision xxx,xxy,xxz,totam,cpham,xpara(2),gms
! input data for grid minimizer
    double precision, dimension(maxel) :: xknown,aphl
! arrays for grid minimization results
    integer, dimension(maxel) :: iphl,icsl,nyphl
! selected kommand and subcommands
    integer kom,kom2,kom3,kom4
! selected output mode for results and the default, list output unit lut
    integer listresopt,lrodef,lut,afo
! integers used for elements, phases, composition sets, equilibria, defaults
    integer iel,iph,ics,ieq,idef,iph2,tupix(2),icond
! for gradients in MU and interdiffusivities
    integer nend
! for mqmqanend, a negative value needed att first call to mqmqa_species
! it is declared in gt3_dd2.F90
! dimension of mugrad for 16x16 matrix 
!CCI
    double precision, allocatable, dimension(:) :: mugrad,mobilities
    double precision, allocatable, dimension(:) :: nsites
    integer, allocatable, dimension(:) :: nkl
    integer nsub
!CCI
!-------------------
! selection of minimizer and optimizer
    integer minimizer,optimizer
! temporary integer variables in loops etc
    integer i1,i2,j4,j5,j2,iax,threads,modelx,jquad
! more temporary integers
    integer jp,kl,svss,language,last,leak,j3,tzcond,eetcond
! and more temporary integers
    integer ll,lokcs,lokph,lokres,loksp,lrot,maxax
! and more temporary integers
    integer mode,ndl,neqdef,noelx,nofc,nopl,nops,nv,nystat,times,fromeq
! temporary matrix
! used to call init_gtp for the NEW command
    integer intv(10)
    double precision dblv(10)
! debugging mqmqma_data%const lines ... lines 5236 ff
    integer ik,ij,kp,s1,thiscon
! ceq%phase_varres(lokcs)%mqmqaf%compvar(j4)%vk_ij,
    type(gtp_phase_varres), pointer :: mqmqavar
!-------------------
    integer iz,jz,ztyp
!    external new_assessment_calfun
!    external calfun
!-------------------
! loop variable when entering constituents of a phase
    integer icon,flc
! array with constituents in sublattices when entering a phase
! only used for interactive entering the mqmqa_constituent
    character, dimension(25) :: const*24
! This is for species in the mqmqa model which may contain commas ","
    character mqmqacon*24,ch3*3
! mqmqa quadbonds
    double precision quadbonds(4)
! for macro and logfile and repeating questions
    logical logok,stop_on_error,once,wildcard,twice,startupmacro,temporary
    logical listzero,maptopbug
! fast elements for Scheil, max 3
    character*2 fast(3)
! unit for logfile input, 0 means no logfile
    integer logfil
! remember default for calculate phase
    integer defcp
! set TRUE by NEW when all data removed
    logical reinit
! set by some extracted subroutines to 100 or 990 to jump to that label
    integer jump
! for state variables as conditions
    integer istv
    double precision coeffs(10),textfontscale
    TYPE(gtp_state_variable), target :: stvrvar
    TYPE(gtp_state_variable), pointer :: stvr
    TYPE(gtp_condition), pointer :: pcond,firstc
! current equilibrium records
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq
    TYPE(gtp_phase_varres), pointer :: parres
! addition record used for listing calculated values
    type(gtp_phase_add), pointer :: addrec
!
    character actual_arg(2)*16
    character cline*256,option*80,aline*128,eqname*24,aux*4
! variable phase tuple
    type(gtp_phasetuple), pointer :: phtup
!CCI
    integer :: indexPrecond, indexSplitsolver, typeOfChange
!CCI
!----------------------------------------------------------------
! here are all commands and subcommands
    integer, parameter :: ncbas=30,nclist=27,ncalc=18,ncent=21,ncread=9
    integer, parameter :: ncam1=18,ncset=27,ncadv=18,ncstat=6,ncdebug=12
    integer, parameter :: nselect=6,nlform=6,noptopt=9,mqmqacc=6,nsetbit=6
    integer, parameter :: ncamph=18,naddph=12,nclph=6,nccph=6,nrej=9,nsetph=6
    integer, parameter :: nsetphbits=15,ncsave=6,nplt=15,nstepop=9
    integer, parameter :: nplt2=18
    integer, parameter :: ninf=15
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
        'SHOW            ','                ','                ',&
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
! subcommands to AMEND first level
! many of these should be subcommands to PHASE
    character (len=16), dimension(ncam1) :: cam1=&
         ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
         'TPFUN_SYMBOL    ','CONSTITUTION    ','QUIT            ',&
         'COMPONENTS      ','GENERAL         ','ASSESSMENT_RESLT',&
         'OPTIMIZING_COEFS','EQUILIBRIUM     ','REDUNDANT_SETS  ',&
         'LINES           ','START_CONSTIT   ','                ']
!-------------------
! subsubcommands to AMEND PHASE
! the UNIQUAC model specified when entering the phase
    character (len=16), dimension(ncamph) :: camph=&
         ['ADDITION        ','COMPOSITION_SET ','DISORDERED_FRACS',&
         '                ','DIFFUSION       ','DEFAULT_CONSTIT ',&
         '                ','FCC_PERMUTATIONS','BCC_PERMUTATIONS',&
         'REMOVE_COMPSETS ','ASYMMETRIES     ','AQUEUS_MODEL    ',&
         'QUASICHEM_MODEL ','FCC_CVM_TETRAHDR','                ',&
         '                ','                ','QUIT            ']
!-------------------
! subsubsubcommands to AMEND PHASE ADDITION
    character (len=16), dimension(naddph) :: caddph=&
         ['MAGNETIC_CONTRIB','QUIT            ','GADDITION       ',&
         'TWOSTATE_LIQUID ','SCHOTTKY_ANOMALY','VOLUME_MODEL1   ',&
         'LOWT_CP_MODEL   ','LIQUID_2STATE   ','                ',&
         'ELASTIC_MODEL_1 ','                ','SMOOTH_CP_STEP  ']
!-------------------
! subcommands to LIST
    character (len=16), dimension(nclist) :: clist=&
         ['DATA            ','SHORT           ','PHASE           ',&
         'STATE_VARIABLES ','BIBLIOGRAPHY    ','MODEL_PARAM_ID  ',&
         'AXIS            ','TPFUN_SYMBOLS   ','                ',&
         'PARAMETER       ','EQUILIBRIA      ','RESULTS         ',&
         'CONDITIONS      ','SYMBOLS         ','LINE_EQUILIBRIA ',&
         'OPTIMIZATION    ','MODEL_PARAM_VAL ','ERROR_MESSAGE   ',&
         'ACTIVE_EQUILIBR ','ELEMENTS        ','EXCELL_CSV_FILE ',&
         'MQMQA_SPECIAL   ','ESTIMAT_ACCURACY','WORKING_DIRECTOR',&
         'QUADS           ','                ','                ']
!-------------------
! subsubcommands to LIST DATA
    character (len=16), dimension(nlform) :: llform=&
        ['SCREEN          ','                ','MACRO           ',&
         '                ','                ','                ']
!-------------------
! subsubcommands to LIST PHASE
! removed VARKAPPA and changed MQMQA_VARIABLES to ASYMMETRIES
    character (len=16), dimension(nclph) :: clph=&
        ['DATA            ','CONSTITUTION    ','MODEL           ',&
         'ASYMMETRIES     ','TERNARY_ASYM    ','                ']
!-------------------
! subsubcommands to LIST OPTIMIZE results
    character (len=16), dimension(noptopt) :: optopt=&
        ['SHORT           ','LONG            ','COEFFICIENTS    ',&
         'GRAPHICS        ','DEBUG           ','MACRO           ',&
         'EXPERIMENTS     ','CORRELATION_MTRX','TC_RSD          ']
!-------------------
! subsubcommands to LIST MQMQA_SPECIALS
    character (len=16), dimension(mqmqacc) :: mqmqalist=&
        ['QUADS           ','ASYMMETRIES     ','DEBUG           ',&
         'EXCESS          ','INTERNALS       ','                ']
!------------------- 
! subcommands to CALCULATE
    character (len=16), dimension(ncalc) :: ccalc=&
         ['TPFUN_SYMBOLS   ','PHASE           ','NO_GLOBAL       ',&
         'TRANSITION      ','QUIT            ','GLOBAL_GRIDMIN  ',&
         'SYMBOL          ','EQUILIBRIUM     ','ALL_EQUILIBRIA  ',&
         'WITH_CHECK_AFTER','TZERO_POINT     ','CAREFULLY       ',&
         'ONLY_GRIDMIN    ','BOSSES_METHOD   ','PARAEQUILIBRIUM ',&
         'LIQUID_EET      ','                ','                ']
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
         'DIRECT          ','XTDB            ','SELECTED_PHASES ',&
         'ENCRYPTED       ','                ','                ']
!-------------------
! subcommands to SAVE
! note SAVE TDB, MACRO, LATEX part of LIST DATA !!
    character (len=16), dimension(ncsave) :: csave=&
        ['TDB             ','                ','QUIT            ',&
         'DIRECT          ','UNFORMATTED     ','XTDB            ']
!-------------------
! subcommands to SET
    character (len=16), dimension(ncset) :: cset=&
         ['CONDITION       ','STATUS          ','ADVANCED        ',&
         '                ','INTERACTIVE     ','REFERENCE_STATE ',&
         'QUIT            ','ECHO            ','PHASE           ',&
         'UNITS           ','LOG_FILE        ','WEIGHT          ',&
         'NUMERIC_OPTIONS ','AXIS            ','INPUT_AMOUNTS   ',&
         'VERBOSE         ','AS_START_EQUILIB','BIT             ',&
         'OPTCOEFF_VARIABL','OPTCOEFF_SCALED ','LMDIF_ACCURACY  ',&
         'RANGE_EXPER_EQU ','OPTCOEFF_FIXED  ','SYSTEM_VARIABLE ',&
         'INITIAL_T_AND_P ','LINEAR_SYSTEM   ','GRID_GENERATOR  ']
! subsubcommands to SET STATUS
    character (len=16), dimension(ncstat) :: cstatus=&
         ['ELEMENT         ','SPECIES         ','PHASE           ',&
         'CONSTITUENT     ','                ','                ']
!        123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET ADVANCED
    character (len=16), dimension(ncadv) :: cadv=&
         ['EQUILIB_TRANSFER','QUIT            ','SYMBOL          ',&
          'GRID_DENSITY    ','SMALL_GRID_ONOFF','MAP_SPECIALS    ',&
          'GLOBAL_MIN_ONOFF','OPEN_POPUP_OFF  ','WORKING_DIRECTRY',&
          'HELP_POPUP_OFF  ','EEC_METHOD      ','LEVEL           ',&
          'NO_MACRO_STOP   ','PROTECTION      ','IGNORE_MACRO_ERR',&
          'XTDB_DEFAULTS   ','DEBUG           ','                ']
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
! Some bits can still be set here by numbers but the text is no longer shown
! most bits are set by AMEND PHASE command
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
          'CONDITIONAL     ','TZERO           ','LIQUID_EET      ',&
          'SCHEIL_GULLIVER ','PARAEQUILIBRIUM ','FAST            ']
!         123456789.123456---123456789.123456---123456789.123456
!-------------------
! subcommands to DEBUG
    character (len=16), dimension(ncdebug) :: cdebug=&
         ['FREE_LISTS      ','STOP_ON_ERROR   ','PARAMETER_STRUCT',&
          'SPECIES         ','TPFUN           ','BROWSER         ',&
          'TRACE           ','SYMBOL_VALUE    ','MAP_STARTPOINTS ',&
          'GRID            ','ASYMMETRIES     ','BOMBMATTA       ']
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
         ['ELEMENTS        ','SPECIES         ','PHASES          ',&
          'QUIT_INFO       ','COMPOSITION_SET ','EQUILIBRIUM     ',&
          'HELP_SYSTEM     ','CONDITIONS      ','DATABASES       ',&
          'CHANGES         ','PHASE_DIAGRAM   ','PROPERTY_DIAGRAM',&
          'STATE_VARIABLES ','                ','                ']
!-------------------
! subcommands to PLOT OPTIONS/ GRAPHICS OPTIONS
! Now there are two levels (using EXTRA) but still a mess
    character (len=16), dimension(nplt) :: cplot=&
        ['RENDER          ','SCALE_RANGES    ','FONT            ',&
         'AXIS_LABELS     ','                ','TITLE           ',&
         'GRAPHICS_FORMAT ','OUTPUT_FILE     ','                ',&
         'QUIT            ','POSITION_OF_KEYS','APPEND          ',&
         'TEXT_LABEL      ','                ','EXTRA           ']
! subsubcommands to PLOT EXTRA
    character (len=16), dimension(nplt2) :: cplot2=&
        ['COLOR           ','LOGSCALE        ','RATIOS_XY       ',&
         'LINE_TYPE       ','MANIPULATE_LINES','PAUSE_OPTION    ',&
         'LOWER_LEFT_TEXT ','TIE_LINES       ','GIBBS_TRIANGLE  ',&
         'QUIT            ','SPAWN           ','NO_HEADING      ',&
         'AXIS_FACTOR     ','GRID            ','                ',&
         '                ','                ','                ']
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
! ceq has no value here!!!  Moved this to gtp3A: initialize_global_parameters
! defaults for several step special
    stepspecial=.FALSE.
! save the working directory (where OC is started?)
    call getcwd(workingdir)
! this is used to save the path to any directory where a macro is started
! initiate command line history
    myhistory%hpos=0
! defaults for optimizer, number of variable coefficients
    nvcoeff=0
! present the software
    write(kou,10)version,trim(linkdate),ocmonversion,gtpversion,hmsversion,&
         smpversion
10  format(/'Open Calphad (OC) software version ',a,', linked ',a,/&
         'with command line monitor version ',i3//&
         'This program is available with a GNU General Public License.'/&
         'either version 2 of the License, or any later version.'/&
         'It includes the General Thermodynamic Package, version ',A,','/&
         "Hillert's equilibrium calculation algorithm version ",A,','/&
         'step/map/plot software version ',A,' using GNUPLOT 5.2 graphics.'/&
         'Numerical routines are extracted from LAPACK and BLAS and'/&
         'the assessment procedure uses LMDIF from ANL.'/&
         'Many users and AI Claude has been of great help to fix problems.'/)
!
! lines starting with !$ will be included when compiling with -fopenmp
!$    write(kou,11)
11  format('Linked with OpenMp for parallel execution')
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Default gnuterminals, edit these as they may not be same on your systems
    graphopt%gnutermid=' '
!    graphopt%status=0 initiated to zero
! Screen is terminal 1
    graphopt%gnutermid(1)='SCREEN '
! default font, not reinitiated if set explicitly
    graphopt%font='Arial '
! MAX 80 characters to set terminal .... HERE FONT AND SIZE IS SET
! test compilation ...
!
#ifdef aqplt
! Aqua plot screen on some Mac systems
    graphopt%gnuterminal(1)='aqua size 600,500 font "'//&
         trim(graphopt%font)//',16"'
! it should be #elif not #elseif .... suck
#elif qtplt
! Qt plot screen on some LINUX systems
    graphopt%gnuterminal(1)='qt size 600,500 font "'//&
         trim(graphopt%font)//',16"'
#elif x11
! x11 plot screen on other LINUX systems    
    graphopt%gnuterminal(1)='x11 size 840,700 font "'//&
         trim(graphopt%font)//',16"'
#else
! wxt default plot screen (used on most Window systems)
    graphopt%gnuterminal(1)='wxt size 840,700 font "'//&
         trim(graphopt%font)//',16"'
! This uses 'start /B ' in front of plot command to spawn plot windows
#endif
    graphopt%filext(1)='  '
! NOTE THAT THE SCREEN PLOT WINDOW ALLOWS YOU TO SELECT FILE OUTPUT
! Postscript
    i1=2
    graphopt%gnutermid(i1)='PS  '
    graphopt%gnuterminal(i1)='postscript color solid fontscale 1.2'
    graphopt%filext(i1)='ps  '
! Adobe Portable Document Format (PDF)
    i1=3
    graphopt%gnutermid(i1)='PDF '
!--------- #ifdef qtplt
! On LINUX ??
!----------#else
! NOTE size is in inch
    graphopt%gnuterminal(i1)='pdf color solid size 6,5 enhanced font "'//&
         trim(graphopt%font)//',16"'
!----------#endif
    graphopt%filext(i1)='pdf  '
! Graphics Interchange Format (GIF)
    i1=4
    graphopt%gnutermid(i1)='GIF  '
    graphopt%gnuterminal(i1)='gif enhanced fontscale 0.7'
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
! if winhlp set also GRKEEP
#ifdef winhlp    
! This uses 'start /B ' in front of plot command to spawn plot windows
    graphopt%status=ibset(graphopt%status,GRKEEP)
#endif
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! jump here after NEW to reinitiallize all local variables also
20  continue
! clear file names
    ocmfile=' '; ocufile=' '; tdbfile=' '; xtdbfile=' '
! clear some other variables
    dummy=' '; name1=' '; name2=' '; name3=' '
    tzcond=0
    eetcond=0
    parael=' '
! initiating the equence of mqmqaquads, declared in models/gtp3_dd2.F90
    mqmqanend=-100
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
! initiate the limit on number of equilibria saved during step/map
    totalsavedceq=0
!
! nopenpopup is declared in metlib3.F90 and dis/allow open popup windows
! it is initiated to FALSE, if user change it will not be reinitiated here
! in init_gtp the first equilibrium record is created and 
! firsteq has been set to that
!
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
    do j4=1,2
       axplotdef(j4)=' '
    enddo
! remove any results from step and map
    nullify(maptop)
    nullify(mapnode)
    nullify(maptopsave)
! entered start equilibria
    nullify(starteqs(1)%p1)
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
#ifdef winhlp
! BROWSER FOR WINDOWS
    browser='C:\"Program Files\Mozilla Firefox"\firefox.exe '
#elif lixhlp
! BROWSER FOR LINUX
    browser='/usr/bin/firefox '
#elif machlp
! BROWSER FOR MAC
    browser='/Applications/Firefox.app/Contents/MacOS/firefox '
#endif
    noochome: if(ochome(1:1).eq.' ') then
! there is no OCHOME environment variable, maybe a local ochelp.html?
       inquire(file='ochelp.html ',exist=logok)
       if(.not.logok) then
          write(*,*)'Warning, no environment variable OCHOME and no help file'
          htmlfile=' '
          htmlhelp=.FALSE.
       else
        write(*,*)'Warning, no environment variable OCHOME but local help file'
          htmlfile='ochelp.html'
          htmlhelp=.TRUE.
       endif
       call init_help(browser,htmlfile)
    else
! there is a OCHOME environment variable
! both LINUX and WINDOWS accept / as separator between directory and file names
       write(*,*)'Found OC home directory (OCHOME): ',trim(ochome)
#ifdef winhlp
! HTML FILE FOR WINDOWS
! normal tex/html help files
       htmlfile=trim(OCHOME)//'\'//'ochelp.html'
#elif lixhlp
! HTML FILE FOR LINUX
       htmlfile=trim(OCHOME)//'/'//'ochelp.html'
#elif machlp
! HTML FILE FOR MAC
       htmlfile=trim(OCHOME)//'/'//'ochelp.html'
#endif
       call init_help(browser,htmlfile)
       if(.not.ochelp%htmlhelp) then
          write(kou,*)'Warning, no file "ochelp.html" at OCHOME or no browser'
          write(kou,*)trim(browser)
          write(kou,*)trim(htmlfile)
          htmlhelp=.FALSE.
       else
          write(kou,*)'Online help by '//trim(browser)//&
               ' and ochelp.html'
       endif
! default directory for databases
       ocbase=trim(ochome)//'/databases'
       cline=trim(ochome)//'/start.OCM '
       inquire(file=cline,exist=logok)
       if(logok) then
          write(*,*)'Reading your startup macro: ',trim(cline)
          last=0
! This just open the file and sets input unit to file
          call macbeg(cline,last,logok)
          startupmacro=.TRUE.
       endif
    endif noochome
! initiate XML defaults
    lowTdef='298.15  '; hightdef='6000    '; bibrefdef='U.N. Known';
    eldef='VA /-'; unary1991=.TRUE.; includemodels=.FALSE.
! running a initial macro file
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
          call macbeg(cline,last,logok)
    endif
! read the command line with gparc to have output on logfile
! NOTE read from macro file if set.
    last=len(aline)
    aline=' '
    cline=' '
    call gparcx(ocprompt,aline,last,5,cline,' ','?TOPHLP')
    j4=0
    if(len_trim(cline).gt.80) then
       write(kou,101)
101    format(' *** Warning: long input lines may be truncated',&
            ' and cause errors')
    endif
! with empty line just prompt again, j4 incremented by eolch
    if(eolch(cline,j4)) goto 100
! with macro command prefix character just prompt again
    if(cline(j4:j4).eq.'@') goto 100
! with the new help facilities "tophlp" is difficult ...
    if(cline(j4:j4+1).eq.'? ') then
! just provide the menu as help
       j4=0
       call q3helpx(cline,j4,cbas,ncbas)
       goto 100
    endif
! Now finally detect the command
    kom=ncomp(cline,cbas,ncbas,last)
    if(kom.le.0) then
       if(kom.lt.0) then
          write(kou,*)'Ambiguous command, available commands are:'
       else
          write(kou,*)'No such command, available commands are:'
       endif
       last=1
       cline=' *'
       call q3helpx(cline,last,cbas,ncbas)
       write(*,*)'An OS command must be prefixed by @'
       goto 100
    else
! check for options .... some of these do not work yet
! one should check for options after each subcommand or value entered ??
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
! save command for help path MAYBE NOT NEEDED ANY LONGER ??
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
!- - - - - - - - - - - - - separating subsubsubcommands
! jump here if there is an inline argument
    main: SELECT CASE(kom)
! command selection
!=================================================================
    CASE DEFAULT
       write(kou,*)'No such command'
       goto 100
!=================================================================
    CASE(1) ! AMEND
!=================================================================
! amend subcommands
!       ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
!        'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
!        'TPFUN_SYMBOL    ','CONSTITUTION    ','QUIT            ',&
!        'COMPONENTS      ','GENERAL         ','ASSESSMENT_RESLT',&
!        'OPTIMIZING_COEFS','EQUILIBRIUM     ','REDUNDANT_SETS  ',&
!        'LINES           ','START_CONSTIT   ','DEBUG           ']
! disable continue optimization
       iph=0
! setting iph=0 means there is no phase in the subcommands below
       kom2=submenu(cbas(kom),cline,last,cam1,ncam1,4,'?TOPHLP')
       amend: SELECT CASE(kom2)
       CASE DEFAULT
          write(kou,*)'No such amendment',kom2
!-------------------------
       case(1) ! amend symbol (of state variable function)
          call gparcx('Symbol name: ',cline,last,1,name1,' ','?Amend symbol')
          call capson(name1)
          do svss=1,nosvf()
             if(name1(1:16).eq.svflista(svss)%name) goto 1020
          enddo
          write(kou,*)'No such symbol'; goto 100
1020      continue
          if(svflista(svss)%status.ne.0) then
! if any bit except SVCONST set we cannot amend it
             if(.not.btest(svflista(svss)%status,SVCONST)) then
                write(*,*)'Symbol is not amendable'; goto 100
             endif
          endif
! No bits or the SVCONST bit is set, it is amendable, get its value
          actual_arg=' '
          xxx=evaluate_svfun_old(svss,actual_arg,1,ceq)
          if(btest(svflista(svss)%status,SVCONST)) then
! symbol is a numeric constant or evaluated explicitly, we can change its value
! value must be set in all equilibria ??
             call gparrdx('Give new value: ',cline,last,xxy,xxx,'?Amend symbol')
             if(buperr.eq.0) then
                call set_putfun_constant(svss,xxy)
                goto 100
             else
! we want to amend something else
                buperr=0
             endif
          endif
! Now we can set one special bit! But first clear input line
          last=len(cline)
          write(kou,1021)
1021      format('You can specify:'/&
               ' V for a symbol evaluated only when referenced explicitly'/&
               ' X for a symbol to be evaluated at a particular equilibrium')
! with SET ADVANCED SYMBOL one can set EXPORT/IMPORT for assessments
          call gparcdx('Please specify V or X',&
               cline,last,1,ch1,'X','?Amend symbol')
          call capson(ch1)
          if(ch1.eq.'V') then
! If V then set bit to evaluate symbol only when explicitly referenced
             svflista(svss)%status=ibset(svflista(svss)%status,SVFVAL)
          elseif(ch1.eq.'X') then
! if X then evaluate symbol only at specific equilibrium?
! For example H298 for experimental data on H(T)-H298
! BEWARE: if equilibria are calculated in threads this must be calculated
! before the parallelization, testing bit EQNOTHREAD
             ll=ceq%eqno
             call gparidx('Specify equilibrium number:',cline,last,&
                  neqdef,ll,'?Amend symbol evaluated at equilib')
! UNFINISHED! Check equilibrium exist or only allow current?
             if(neqdef.le.1 .or. neqdef.gt.noeq()) then
                write(*,*)'No such equilibrium'; goto 100
             endif
             svflista(svss)%status=ibset(svflista(svss)%status,SVFEXT)
             svflista(svss)%eqnoval=neqdef
! set status bit that this equilibrium must be calculated before parallel calc
             ceq%status=ibset(ceq%status,EQNOTHREAD)
             write(*,*)'The value of this symbol calculated in equilibrium: ',&
                  neqdef
             goto 100
          else
             write(kou,*)'Illegal letter "',ch1,'"'
          endif
!-------------------------
       case(2) ! amend element
          call gparcx('Element symbol: ',cline,last,1,elsym,' ',&
               '?Amend element')
          call find_element_by_name(elsym,iel)
          if(gx%bmperr.ne.0) goto 100
          call get_element_data(iel,elsym,name1,dummy,mass,h298,s298)
          if(gx%bmperr.ne.0) goto 100
          write(*,'(a)')'You are only allowed to change the mass'
          call gparrdx('New mass: ',cline,last,xxx,mass,'?Amend element')
          call new_element_data(iel,elsym,name1,dummy,xxx,h298,s298)
!-------------------------
       case(3) ! amend species
          call gparcx('Species symbol: ',cline,last,1,name1,' ',&
               '?Amend species')
          call find_species_record(name1,loksp)
          if(gx%bmperr.ne.0) goto 100
          write(*,'(a)')'You can only amend UNIQAC area and segments'
          call gparrdx('UNIQAC surface area (q): ',cline,last,xxx,one,&
               '?Amend species')
          if(xxx.le.zero) then
             write(*,'(a)')'Area must be >0, set to default 1.00'
             xxx=one
          endif
          call gparrdx('UNIQAC segments (r): ',cline,last,xxy,one,&
               '?Amend species')
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
          call gparcx('Phase name: ',cline,last,1,name1,' ','?Amend for phase')
          if(buperr.ne.0) goto 990
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
! sometimes lokph is used below
          call get_phase_record(iph,lokph)
          call get_phasetup_name(iph,name1)
!
          kom3=submenu('Amend for phase '//trim(name1),&
               cline,last,camph,ncamph,2,'?TOPHLP')
          amendphase: SELECT CASE(kom3)
! subsubcommands to AMEND PHASE
!         ['ADDITION        ','COMPOSITION_SET ','DISORDERED_FRACS',&
!         '                 ','DIFFUSION       ','DEFAULT_CONSTIT ',&
!         'TERNARY_EXTRAPOL','FCC_PERMUTATIONS','BCC_PERMUTATIONS',&
!         'REMOVE_COMPSETS ','ASYMMETRIES     ','AQUEUS_MODEL    ',&
!         'QUASICHEM_MODEL ','FCC_CVM_TETRAHDR','                ',&
!         '                ','                ','QUIT            ']
!....................................................
          CASE DEFAULT
             write(kou,*)'Amend phase subcommand error'
!....................................................
          case(1) ! amend phase addition
             kom4=submenu('Addition of',cline,last,caddph,naddph,1,&
                  '?TOPHLP')
!         ['MAGNETIC_CONTRIB','QUIT            ','GADDITION      ',&
!         'TWOSTATE_LIQUID ','SCHOTTKY_ANOMALY','VOLUME_MODEL1   ',&
!         'LOWT_CP_MODEL   ','LIQUID_2STATE   ','                ',&
!         'ELASTIC_MODEL_A ','QUASICHEM_MODEL ','FCC_CVM_TETRAHDR']
!
             amendphaseadd: SELECT CASE(kom4)
             case default
                write(*,*)'No such addition'
! Inden magnetism
             case(1) ! amend phase <name> magnetic contribution
                call ui_amend_phase_magnetic(cline,last,lokph)
!....................................................
             case(2) ! QUIT
                goto 100
!....................................................
             case(3) ! amend phase ... addition gaddition
! different additions can be added for different composition sets
                call get_phase_compset(iph,ics,lokph,lokcs)
                if(gx%bmperr.ne.0) goto 100
                if(allocated(ceq%phase_varres(lokcs)%addg)) then
                   xxy=ceq%phase_varres(lokcs)%addg(1)
                else
! maybe we will use more terms later ....
                   xxy=zero
                   allocate(ceq%phase_varres(lokcs)%addg(1))
                endif
!\hypertarget{Amend phase Gaddition}{}
                call gparrdx('Addition to G in J/FU (formula units): ',&
                     cline,last,xxx,xxy,'?Amend Gaddition')
                ceq%phase_varres(lokcs)%addg(1)=xxx
! set bit that this should be calculated
                ceq%phase_varres(lokcs)%status2=&
                     ibset(ceq%phase_varres(lokcs)%status2,CSADDG)
!....................................................
             case(4,8) ! amend phase <name> add liquid_2state/twostate_liquid
                write(kou,667)
667             format('This addition require LNTH parameters for the',&
                     ' Einstein T of the amorphous state'/'and G2 parameters',&
                     ' for the transition to the liquid state.')
! WRONG IDEA to set bit to allow G2 to be composition independent
! ensure ch1 is a captial letter!
! if ch1 is N then the addition record will have the twostatemodel2(=12) value
!     and the PH2STATE in the phase record must be set also:
!     But as phlista is protected it is set inside add_addrecord
                modelx=twostatemodel1
! inside add_addrecord modelx can be changed to twostatemodel2 if G2 fixed
                ch1='Y'
                call add_addrecord(lokph,ch1,modelx)
                call gparcdx('Is the low T calculated per mole atoms?',&
                     cline,last,1,ch1,'Y','?Add per formula unit')
! The CP model calculates a molar Gibbs energy, must be multiplied with
! the number of atoms in the phase.
                if(ch1.eq.'Y' .or. ch1.eq.'y') then
                   call setpermolebit(lokph,modelx)
                endif
!....................................................
             case(5) ! amend phase <name> addition Schottky anomaly
                call add_addrecord(lokph,' ',schottkyanomaly)
                write(*,668)
668             format('This addition requires the TSCH and CSCH parameters')
!....................................................
! VOLUME MODEL1
             case(6) ! volume model1
                call add_addrecord(lokph,' ',volmod1)
                write(*,*)'Added volume model 1'
!....................................................
! Einstein low T model
             case(7) ! amend phase <name> LowT_CP_model
                call add_addrecord(lokph,' ',einsteincp)
                write(*,*)'This addition requires the LNTH parameter'
                call gparcdx('Is the addition calculated for one mole atoms? ',&
                     cline,last,1,ch1,'Y','?Add per formula unit')
! The CP model calculates a molar Gibbs energy, must be multiplied with
! the number of atoms in the phase. j2 set above to the addition type
                if(ch1.eq.'Y' .or. ch1.eq.'y') then
                   call setpermolebit(lokph,einsteincp)
                endif
!....................................................
             case(9,11) ! not used
!....................................................
             case(10) ! amend phase elastic model
                write(*,*)'This addition is not yet implemented'
!....................................................
             case(12) ! amend phase ... smooth-Cp-step
                call add_addrecord(lokph,' ',secondeinstein)
                call gparcdx('Is the addition calculated for one mole? ',&
                     cline,last,1,ch1,'Y','?Add per formula unit')
! The CP model calculates a molar Gibbs energy, must be multiplied with
! the number of atoms in the phase. j2 set above to the addition type
                if(ch1.eq.'Y' .or. ch1.eq.'y') then
                   call setpermolebit(lokph,secondeinstein)
                endif
                write(*,672)
672             format('This addition recures the THT2 and DCP2 parameters')
! The smooth CP model calculates a molar Gibbs energy, must be multiplied with
! the number of atoms in the phase. j2 set above to the addition type
             end select amendphaseadd
!-------------------------
          case(2) ! amend phase <name> composition set add/remove
             call gparcdx('Add new set? ',cline,last,1,ch1,'Y ','?Add new cs')
             if(buperr.ne.0) goto 990
             if(ch1.eq.'Y' .or. ch1.eq.'y') then
                call gparcx('Prefix: ',cline,last,1,prefix,' ',&
                     '?Add new cs')
                call gparcx('Suffix: ',cline,last,1,suffix,' ',&
                     '?Add new cs')
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
             call ui_amend_phase_disordered(cline,last,iph)
!....................................................
          case(4) ! Not used
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
          case(7) ! this subcommand is empty
!....................................................
!\hypertarget{Amend FCC-permutations}{}
          case(8) ! amend phase ... FCC_PERMUTATIONS
             if(check_minimal_ford(lokph)) goto 100
             call set_phase_status_bit(lokph,PHFORD)
!....................................................
!\hypertarget{Amend BCC-permutations}{}
          case(9) ! amend phase ... BCC_PERMUTATIONS
             if(check_minimal_ford(lokph)) goto 100
             call set_phase_status_bit(lokph,PHBORD)
!....................................................
          case(10) ! amend phase <...> remove_compsets
             write(*,*)'PMON: delete unstable composition sets'
             call delete_unstable_compsets(lokph,ceq)
!....................................................
          case(11) ! AMEND PHASE ... ASYMMETRIES for MQMQA phase
             call ui_amend_phase_asymmetries(cline,last,iph,ceq)
!....................................................
          case(12) ! amend phase ... aqueous model
             write(*,*)'Not implemented yet'
!....................................................
          case(13) ! amend phase ... quasicemichal model (2 versions)
             call gparidx('Quasichemical type: ',cline,last,jp,3,&
                  '?Amend quasichemical')
             if(jp.lt.0 .or. jp.gt.3) then
                write(*,*)'Value must be between 1 and 3'
             else
                qcmodel=jp
             endif
! Future model bits
!....................................................
          case(14) ! amend phase ... FCC_CVM_TETRAHEDRON MODEL
             write(kou,*)'Not implemented yet'
!....................................................
          case(15,16,17,18) ! amend phase ... unused
          END SELECT amendphase
!------------------------- end of amend phase
       case(5) ! amend parameter
! Claude 2026.10.15: as ENTER PARAMETER but an existing parameter is replaced
          if(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'You must enter a phase before'
             goto 100
          endif
          amendparameter=.true.
          call enter_parameter_interactivly(cline,last,0)
          amendparameter=.false.
!-------------------------
       case(6) ! amend bibliography (in gtp3D)
          call enter_bibliography_interactivly(cline,last,1,j4)
!-------------------------
       case(7) ! amend TPFUN symbol
          write(kou,*)' *** Dangerous if you have several equilibria!'
          call gparcx('TP-fun symbol: ',cline,last,1,name1,' ',&
               '?Amend TPfun')
          call find_tpsymbol(name1,idef,xxx)
          if(gx%bmperr.ne.0) then
             write(*,*)'Ambiguouos or unknown symbol'; goto 990
          endif
          if(idef.eq.0) then
! it is a function , this call just read the function starting with low T etc.
             call enter_tpfun_interactivly(cline,last,funstring,jp)
! this stores the tpfun, lrot<0 means the symbol already exists
             lrot=-1
! last argument -1 means not reading from TDB file
             call store_tpfun(name1,funstring,lrot,-1)
             if(gx%bmperr.ne.0) goto 990
! mark functions not calculated.  This should be done in all ceq ...
             ceq%eq_tpres(lrot)%tpused(1)=-one
             ceq%eq_tpres(lrot)%tpused(2)=-one
          elseif(idef.eq.2) then
             write(*,*)'You cannot change an optimizing coefficients'
             goto 100
          else
! it is a constant, you can change if
             call gparrdx('Value: ',cline,last,xxy,xxx,'?Amend TPfun')
             call capson(name1)
             call store_tpconstant(name1,xxy)
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
          call gparcdx('Give all new components: ',cline,last,&
               5,aline,line,'?Amend components')
! option is a character with the new components ...
          call amend_components(aline,ceq)
          if(gx%bmperr.ne.0) goto 990
!-------------------------
       case(11) ! amend general
! this asks for some general things ....
          call amend_global_data(cline,last)
!
!-------------------------
       case(12) ! amend assessment result
          call ui_amend_assessment_result(cline,last)
!-------------------------
       case(13) ! amend OPTIMIZING_COEFF, (rescale or recover)
          call ui_amend_optimizing_coeffs(cline,last,lut)
!-------------------------
       case(14) ! AMEND EQUILIBRIUM intended to add to experimental list
          write(*,*)'Not implemented yet'
!-------------------------
       case(15) ! AMEND REDUDANT composition sets
          write(*,*)'This will set all unstable additional composition sets',&
               ' as suspended'
          ll=0
          call suspend_unstable_sets(ll,ceq)
!-------------------------
       case(16) ! AMEND LINEs of calculated equilibria
! possible amendment of all stored equilibria as ACTIVE or INACTIVE
          call amend_stored_equilibria(axarr,maptop)
!-------------------------
       case(17) ! AMEND START_CONSTITUTION for assessments
! copy constitutions from one equilibrium to another, to handle miscibility gaps
! Default from "previous"
          if(.not.allocated(firstash%coeffstate)) then
             write(kou,*)'This is used during assessments'
! to copy start values from one experimental equilibria to another
             goto 100
          endif
          ll=max(1,ceq%eqno-1)
          call gparidx('From equilibrium number: ',cline,last,&
               fromeq,ll,'?Amend start_constitution')
! copy constitutions of non-suspended phases from fromeq to current equilibrium
          call copyfracs(fromeq,ceq)
!-------------------------
       case(18) !       removed AMEND DEBUG_MQMQA
       END SELECT amend
!=================================================================
! calculate subcommands
!=================================================================
!        ['TPFUN_SYMBOLS   ','PHASE           ','NO_GLOBAL       ',&
!         'TRANSITION      ','QUIT            ','GLOBAL_GRIDMIN  ',&
!         'SYMBOL          ','EQUILIBRIUM     ','ALL_EQUILIBRIA  ',&
!         'WITH_CHECK_AFTER','TZERO_POINT     ','CAREFULLY       ',&
!         'ONLY_GRIDMIN    ','BOSSES_METHOD   ','PARAEQUILIBRIUM ',&
!         'LIQUID_EET      ','                ','                ']
    CASE(2) ! CALCULATE
       kom2=submenu(cbas(kom),cline,last,ccalc,ncalc,8,'?TOPHLP')
       calculate: SELECT CASE(kom2)
       CASE DEFAULT
          write(kou,*)'No such calculate command'
          goto 100
!-------------------------
       CASE(1) ! calculate TPFUN symbols , use current values of T and P
          call ui_calculate_tpfun(cline,last,lut,ceq)
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
          defcp=1
          kom3=submenu('Calculate what for phase?',cline,last,ccph,nccph,defcp,&
               '?TOPHLP')
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
          case(1) ! calculate_phase .. calculate phase < > only G
             call calcg(iph,ics,0,lokres,ceq)
             if(gx%bmperr.ne.0) goto 990
             parres=>ceq%phase_varres(lokres)
             write(lut,2031)(cpham*rgast*parres%gval(j4,1),j4=1,4)
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
             write(lut,2031)(cpham*rgast*parres%gval(j4,1),j4=1,4)
             write(lut,2033)
2033         format('dG/dy: ... dG/dy_i is NOT THE CHEMICAL POTENTIAL of i!')
             write(lut,2041)(rgast*parres%dgval(1,j4,1),j4=1,nofc)
2041         format(9x,5e14.4)
!.......................................................
          case(3) ! calculate phase < > all derivatives
             call gparidx('Number of times: ',cline,last,times,1,&
                  '?Calculate phase ... loop')
! attempt to measure calcg_interal bottlenecks
             call tabder(iph,ics,times,ceq)
             if(gx%bmperr.ne.0) goto 990
! write 20 values
             write(*,2042)
2042         format('Values are per mole formula unit'/&
                  ' NOTE THAT dG/dy_i is NOT THE CHEMICAL POTENTIAL of i!')
!.......................................................
          case(4,5) ! calculate phase with constitution_adjustment
             call ui_calculate_phase_adjust(cline,last,iph,ics,kom3,lut,ceq)
!.......................................................
          case(6) ! Quit
          END SELECT calcphase
! set bits to warn that listings may be inconsistent
          ceq%status=ibclr(ceq%status,EQNOEQCAL)
          ceq%status=ibset(ceq%status,EQINCON)
!---------------------------------- end of calculate phase
       case(3) ! calculate equilibrium without initial global minimization
          if(minimizer.eq.1) then
! Lukas minimizer, first argument=0 means do not use grid minimizer
             write(kou,*)'Not implemented yet'
          else
             call calceq2(0,ceq)
! check that invmat allocated and stored
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
! This command is depreciated as it has been misunderstood
          write(*,2101)
2101      format('This command has been removed as it has been misunderstood.',&
               /'It calculated points on the Gibbs energy curves of the ',&
               ' different phases'/'to be used as start points for the',&
               ' actual minimization.'/'It is replaced by the command ',&
               '"CALCULATE ONLY_GRIDMIN".'//&
               'In order to calculate the equilibrium use "C E".'/)
          goto 100
!---------------------------------------------------------------
       case(7) ! calculate symbol
! to calculate derivatives this must be in the minimizer module
          call gparcdx('Name ',cline,last,1,name1,'*','?Calculate symbol')
! always calculate all state variable functions as they may depend on eachother
          call meq_evaluate_all_svfun(-1,ceq)
          if(gx%bmperr.ne.0) then
! ignore error unless inside macro
             write(*,*)'UI: error calculating all functions', gx%bmperr
             if(kiu.ne.kiud) goto 990
             gx%bmperr=0
          endif
          if(name1(1:1).eq.'*') then
! this calculate them again ... and lists the values
             call meq_evaluate_all_svfun(lut,ceq)
          else
! This code is also used in SHOW (command 25)
             call capson(name1)
             call find_svfun(name1,istv)
             if(gx%bmperr.ne.0) goto 990
             mode=1
             actual_arg=' ' 
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
             write(kou,*)'No longer available'
          else
             call calceq2(1,ceq)
             if(gx%bmperr.eq.4204) then
! if the error code is "too many iterations" try without grid minimizer
! it converges in some cases
                write(*,2048)gx%bmperr
2048            format('Error ',i5,', cleaning up and trying harder')
                gx%bmperr=0
                call calceq2(0,ceq)
             endif
          endif
! calceq2 set appropriate bits for listing
          if(gx%bmperr.ne.0) then
! Claude 2026.10.15: explain failures at exact salt composition
             j4=gx%bmperr; gx%bmperr=0
             if(mqmqa_exact_salt(ceq)) then
! generalized and modified Claude explanation
                write(kou,20501)
20501           format(' *** MM The overall composition is exactly that of a',&
                     ' single salt.'/&
                     ' *** or outside the range of the entered phases.',&
                     ' *** The chemical potentials are then undetermined and',&
                     ' the calculation'/' *** is meaningless.')
                gx%bmperr=j4
             else
                gx%bmperr=j4
                if(gx%bmperr.eq.4204) write(*,2049)
             endif
2049         format('If the conditions allow using the global minimizer, ',&
                  ' try "CALCULATE CAREFULLY"')
             ceq%status=ibset(ceq%status,EQFAIL)
             goto 990
          endif
!---------------------------------------------------------------
       case(9) ! calculate all equilibria
          call ui_calculate_all_equilibria(cline,last,lut)
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
       case(11) ! CALCUALTE TZERO
! The degrees of freedom must be zero
          ll=degrees_of_freedom(ceq)
          if(ll.ne.0) then
             write(*,*)'You must have zero degrees of freedoms for this'
             goto 100
          endif
          write(kou,*)'You should have calculated an equilibrium',&
               ' close to the T0 point'
! ask for 2 phases and which condition to vary
! try to remember the phases ... user may try several times
          if(dummy(1:1).ne.' ') dummy=name2
          call gparcdx('First phase ',cline,last,1,name2,dummy,'?Tzero')
          call find_phase_by_name(name2,iph,ics)
          if(gx%bmperr.ne.0) goto 990
          if(dummy(1:1).ne.' ') dummy=name3
          call gparcdx('Second phase ',cline,last,1,name3,dummy,'?Tzero')
          call find_phase_by_name(name3,iph2,ics)
          if(gx%bmperr.ne.0) goto 990
          dummy=name3
          call list_conditions(kou,ceq)
          if(tzcond.eq.0) then
             j2=1
          else
             j2=tzcond
          endif
          call gparidx('Release condition number',cline,last,tzcond,j2,'?Tzero')
          call tzero(iph,iph2,tzcond,xxx,ceq)
          if(gx%bmperr.ne.0) goto 990
          write(*,*)'Equal Gibbs energy at:'
          call list_conditions(kou,ceq)
! a warning when list equilibria
          ceq%status=ibset(ceq%status,EQINCON)
!-------------------------------------------------------
       case(13) ! Only gridmin no merge
! extract values for mass balance calculation from conditions
          call extract_massbalcond(ceq%tpval,xknown,totam,ceq)
          if(gx%bmperr.ne.0) goto 990
! debug output
! generate grid and find the phases and constitutions for the minimum.
! Note: global_gridmin calculates for total 1 mole of atoms, not totam
! iphl is dimensioned (1:maxel), maxel=100, it is destroyed inside merge_grid ..
          temporary=.false.
          if(.not.btest(globaldata%status,GSNOMERGE)) then
             temporary=.true.
             globaldata%status=ibset(globaldata%status,GSNOMERGE)
          endif
          call global_gridmin(1,ceq%tpval,xknown,nv,iphl,icsl,&
               aphl,nyphl,cmu,ceq)
          if(temporary) globaldata%status=ibclr(globaldata%status,GSNOMERGE)
          if(gx%bmperr.ne.0) goto 990
! In some cases "c n" converges better? if we scale with the total amount here??
          do j4=1,nv
             call get_phase_compset(iphl(j4),icsl(j4),lokph,lokcs)
             ceq%phase_varres(lokcs)%amfu=totam*ceq%phase_varres(lokcs)%amfu
          enddo
! if set clear this bit so we can list the equilibrium
          if(btest(ceq%status,EQNOEQCAL)) ceq%status=ibclr(ceq%status,EQNOEQCAL)
! we should write phase tuples ... ?? 
          write(kou,2102)nv,(iphl(j4),icsl(j4),j4=1,nv)
2102      format('Number of stable phases ',i2/13(i4,i2))
!-------------------------------------------------------
       case(12,14) ! Calculate carefully the equilibrium (bosses_method)
! extract values for mass balance calculation from conditions
          call system_clock(count=startoftime)
          call cpu_time(start2)
          call extract_massbalcond(ceq%tpval,xknown,totam,ceq)
          if(gx%bmperr.ne.0) goto 990
          call global_gridmin(1,ceq%tpval,xknown,nv,iphl,icsl,&
               aphl,nyphl,cmu,ceq)
          if(gx%bmperr.ne.0) goto 990
          j4=1
          if(kom.eq.14) j4=0
! first parameter 0 means bosses_method, 1 means carefully
          call calculate_carefully(j4,ceq)
          if(gx%bmperr.ne.0) goto 990
          call system_clock(count=endoftime)
          call cpu_time(finish2)
          call get_state_var_value('GS ',xxx,name1,ceq)
          if(gx%bmperr.ne.0) gx%bmperr=0
          write(*,654)finish2-start2,endoftime-startoftime,xxx
654       format('Final result: ',1pe12.4,' cpu seconds, ',&
               i7,' cc, G=',1pe15.7,' J/mol')
!-------------------------------------------------------
       case(15) ! CALCULATE PARAEQUILIBRIUM
          write(kou,876)
876       format('You should have calculated an equilibrium',&
               ' close to the paraequilibrium'/&
               'and suspended all but the two phases involved')
! ask for 2 phases and the fast diffusing element
! try to remember the phases ... user may use the command several times
          if(dummy(1:1).ne.' ') dummy=name2
          call gparcdx('Matrix phase ',cline,last,1,name2,dummy,&
               '?Calculate paraeq')
          call find_phasetuple_by_name(name2,tupix(1))
          if(gx%bmperr.ne.0) goto 990
          if(dummy(1:1).ne.' ') dummy=name3
          call gparcdx('Growing phase ',cline,last,1,name3,dummy,&
               '?Calculate paraeq')
          call find_phasetuple_by_name(name3,tupix(2))
          if(gx%bmperr.ne.0) goto 990
          dummy=name3
          call gparcdx('Fast diffusing element',cline,last,1,&
               elsym,parael,'?Calculate paraeq')
          call capson(elsym)
          call find_element_by_name(elsym,icond)
          parael=elsym
          call calc_paraeq(tupix,icond,xpara,meqrec,meqrec1,ceq)
! set a warning when list result
          ceq%status=ibset(ceq%status,EQINCON)
          if(gx%bmperr.ne.0) goto 990
          write(kou,877)trim(elsym),xpara
877       format(/'Paraequilibrium fractions of ',a,': ',2F10.6/&
               'Please note that the phase amounts are not adjusted,',&
               ' only the compositions'/)
!-------------------------------------------------------
       case(16) ! CALCULATE LIQUID_EET, check how TZERO is calculated
          ll=degrees_of_freedom(ceq)
          if(ll.ne.0) then
             write(*,*)'You must have zero degrees of freedoms for this'
             goto 100
          endif
          write(kou,878)
878       format(/'The Equi-Entropy T is the temperature where a solid phase ',&
               'has the same entropy'/'as the liquid phase.  This can be ',&
               'considered as the limit of stability of the'/'solid phase ',&
               'even if its Gibbs energy may become lower than the Gibbs ',&
               'energy'/'of the liquid at an even higher T.  This command ',&
               'varies T or a composition'/'of the phase to find EET. ',&
               'It may fail if there is no EET for the phase.'&
               //'You should have calculated an equilibrium',&
               ' close to the EET point already.')
! ask for liquid and another phase name and a condition to vary
! remember the phases ... a user may try several times
          dummy='LIQUID '
          call gparcdx('The liquid phase ',cline,last,1,name2,dummy,'?EET ')
          call find_phase_by_name(name2,iph,ics)
          if(gx%bmperr.ne.0) goto 990
          if(dummy(1:1).ne.' ') dummy=name3
          call gparcdx('The solid phase ',cline,last,1,name3,dummy,'?EET ')
          call find_phase_by_name(name3,iph2,ics)
          if(gx%bmperr.ne.0) goto 990
          dummy=name3
          call list_conditions(kou,ceq)
          if(eetcond.eq.0) then
             j2=1
          else
             j2=eetcond
          endif
         call gparidx('Release condition number',cline,last,eetcond,j2,'?EET')
          call liquid_eet(iph,iph2,eetcond,xxx,ceq)
          if(gx%bmperr.ne.0) then
! Failed calculation normally means T is negative, try to set it to 1
! does not change T ??? Why?
             goto 990
          endif
          write(*,'(/a)')'The two phases has equal entropy at:'
          call list_conditions(kou,ceq)
! a warning when list equilibria
          ceq%status=ibset(ceq%status,EQINCON)
!-------------------------------------------------------
       case(17) ! CALCULATE ??
          write(kou,*)'Not implemented yet'
!-------------------------------------------------------
       case(18) ! CALCULATE ??
          write(kou,*)'Not implemented yet'
       END SELECT calculate
!=================================================================
! SET SUBCOMMANDS
!=================================================================
!         ['CONDITION       ','STATUS          ','ADVANCED        ',&
!         'LEVEL           ','INTERACTIVE     ','REFERENCE_STATE ',&
!         'QUIT            ','ECHO            ','PHASE           ',&
!         'UNITS           ','LOG_FILE        ','WEIGHT          ',&
!         'NUMERIC_OPTIONS ','AXIS            ','INPUT_AMOUNTS   ',&
!         'VERBOSE         ','AS_START_EQUILIB','BIT             ',&
!         'OPTCOEFF_VARIABL','OPTCOEFF_SCALED ','LMDIF_ACCURACY  ',&
!         'RANGE_EXPER_EQU ','OPTCOEFF_FIXED  ','SYSTEM_VARIABLE ',&
!         'INITIAL_T_AND_P ','LINEAR_SYSTEM   ','GRID_GENERATOR  ']
    CASE(3) ! SET SUBCOMMANDS
! disable continue optimization
       kom2=submenu(cbas(kom),cline,last,cset,ncset,1,'?TOPHLP')
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
          kom3=submenu(name1,cline,last,cstatus,ncstat,3,'?TOPHLP')
          setstatus: SELECT CASE(kom3)
!.................................................................
          CASE DEFAULT
             write(kou,*)'Set status subcommand error'
!.................................................................
          case(1) ! set status element suspend/restore
             call gparcx('Element symbol: ',cline,last,1,name1,' ',&
                  '?Set status element')
             call find_element_by_name(name1,iel)
             if(gx%bmperr.ne.0) goto 100
             call gparcdx('New status: ',cline,last,1,ch1,'SUSPEND',&
                  '?Set status element')
             call capson(ch1)
             if(ch1.eq.'S') then
                call change_element_status(name1,1,ceq)
             else
! restore element
                call change_element_status(name1,0,ceq)
             endif
!.................................................................
          case(2) ! set status species suspend/restore
             call gparcx('Species symbol: ',cline,last,1,name1,' ',&
                  '?Set status species')
             call find_species_record(name1,loksp)
             if(gx%bmperr.ne.0) goto 100
             call gparcdx('New status: ',cline,last,1,ch1,'SUSPEND',&
                  '?Set status species')
             call capson(ch1)
             if(ch1.eq.'S') then
                call change_species_status(name1,1,ceq)
             else
                call change_species_status(name1,0,ceq)
             endif
!.................................................................
          case(3) ! set status phase (ENTERED, FIX, DORMANT, SUSPEND or HIDDEN)
             call ui_set_status_phases(cline,last,ceq)
!.................................................................
          CASE(4,5,6) ! set status constituent 
             write(kou,*)'Not implemented yet'
          END SELECT setstatus
!-----------------------------------------------------------
       case(3) ! set ADVANCED subcommands
! default is WORKING_DIRECT
! subsubcommands to SET ADVANCED
!    character (len=16), dimension(ncadv) :: cadv=&
!         ['EQUILIB_TRANSFER','QUIT            ','SYMBOL          ',&
!          'GRID_DENSITY    ','SMALL_GRID_ONOFF','MAP_SPECIALS    ',&
!          'GLOBAL_MIN_ONOFF','OPEN_POPUP_OFF  ','WORKING_DIRECTRY',&
!          'HELP_POPUP_OFF  ','EEC_METHOD      ','LEVEL           ',&
!          'NO_MACRO_STOP   ','PROTECTION      ','IGNORE_MACRO_ERR',&
!          'XTDB_DEFAULTS   ','MACRO_DIRECTORY ','DEBUG           ']
          name1='Advanced command'
! changed default to working_directory
          kom3=submenu(name1,cline,last,cadv,ncadv,9,'?TOPHLP')
          setadvanced: select case(kom3)
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
             call gparidx('Saved ceq number',cline,last,icon,1,&
                  'Set advanced transfer')
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
          case(3) ! SET ADVANCED SYMBOL to connect with TP constants
! Set the SVCONST bit and allow setting a new value at the same time
             if(.not.allocated(firstash%eqlista)) then
                write(kou,*)'Not allowed as no assessment coefficients'
                goto 100
             endif
             call gparcx('Symbol name: ',cline,last,1,name1,' ',&
                  '?Set adv symbol')
             call capson(name1)
             do svss=1,nosvf()
                if(name1(1:16).eq.svflista(svss)%name) exit
             enddo
             if(svss.gt.nosvf()) then
                write(kou,*)'No such symbol'; goto 100
             endif
             if(.not.btest(svflista(svss)%status,SVCONST)) then
                write(kou,*)'Can only be used for constants'; goto 100
             endif
! Here the symbols can be set to be EXPORTED or EXPORTED to assess coeff
             call gparix('Coefficient index, 0-99?',cline,last,jp,0,&
                  '?Export symbol')
             if(jp.le.0 .or. jp.gt.size(firstash%coeffvalues)) then
                write(*,*)'No such coefficent'; goto 100
             endif
! nv is index to TP function for coefficient
             nv=firstash%coeffindex(jp)
             call gparcdx('Export or Import?',cline,last,1,ch1,'E',&
                  'EXPORT SYMBOL')
             if(ch1.eq.'E') then
! UNFINISHED
                svflista(svss)%status=ibset(svflista(svss)%status,SVEXPORT)
! use firstash% record to handle value transfer                
! probably more firstash variables should be set
                firstash%coeffvalues(jp)=svflista(svss)%svfv
! trying to set bit and copy value to TPFUN
! impossible as tpfuns is private to general_thermodynamic_package !!
! save the index to coefficient in %eqnoval !!
                svflista(svss)%eqnoval=jp
             else
! UNFINISHED this must also set a bit in the TP function/assessment record
                svflista(svss)%status=ibset(svflista(svss)%status,SVIMPORT)
! trying to set bit and copy value to TPFUN
             endif
!.................................................................
          case(4) ! SET ADVANCED GRID_DENSITY
             call gparidx('Level: ',cline,last,ll,1,'?Set adv grid-density')
             if(ll.eq.0) then
! this set GSOGRID, small grid and clears GSXGRID
                globaldata%status=ibset(globaldata%status,GSOGRID)
                globaldata%status=ibclr(globaldata%status,GSXGRID)
                globaldata%status=ibclr(globaldata%status,GSYGRID)
                write(kou,3110)'Sparse','set'
             elseif(ll.eq.1) then
! DEFAULT, all gridbits are cleared
                globaldata%status=ibclr(globaldata%status,GSXGRID)
                globaldata%status=ibclr(globaldata%status,GSOGRID)
                globaldata%status=ibclr(globaldata%status,GSYGRID)
                write(kou,3110)'Normal','set'
3110            format(a,' grid ',a)
             elseif(ll.eq.2) then
! set GSXGRID (and clear GSOGRID and GSYGRID)
                globaldata%status=ibclr(globaldata%status,GSOGRID)
                globaldata%status=ibset(globaldata%status,GSXGRID)
                globaldata%status=ibclr(globaldata%status,GSYGRID)
                write(kou,3110)'Dense','set'
             elseif(ll.eq.3) then
! set GSYGRID (and clear GSXGRID and GSOGRID)
                globaldata%status=ibclr(globaldata%status,GSOGRID)
                globaldata%status=ibclr(globaldata%status,GSXGRID)
                globaldata%status=ibset(globaldata%status,GSYGRID)
                write(kou,3110)'Only level 0, 1 and 2 implemented'
             else
                write(*,*)'Only level 0, 1 and 2 implemented'
             endif
!.................................................................
          case(5) ! SET ADVANCED SMALL_GRID_ONOFF
! replaced by setting grid_density to 0
             write(*,*)'Please use SET ADVANCED GRID 0'
             continue
! set GSOGRID and clear GSXGRID if set
!.................................................................
          case(6) ! MAP_SPECIALS
             ll=mapglobalcheck
             if(ll.le.0) ll=10
             call gparidx('Global test interval during STEP/MAP?: ',&
                  cline,last,mapglobalcheck,ll,'?Set adv global onoff')
!.................................................................
          case(7) ! GLOBAL_MIN_ONOFF
             call gparcx('Turn global minimization off?: ',cline,last,&
                  1,ch1,'N','?Set adv global onoff')
             if(ch1.eq.'Y' .or. ch1.eq.'y') then
                globaldata%status=ibset(globaldata%status,GSNOGLOB)
                write(*,*)'Global minimizer turned off'
             else
                globaldata%status=ibclr(globaldata%status,GSNOGLOB)
                write(*,*)'Global minimizer turned on'
             endif
!.................................................................
          case(8) ! OPEN_POPUP_OFF
             call gparcdx('Turn off popup for open? ',cline,last,&
                  1,ch1,'Y','?Set adv open popup')
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
! case 17 used to be set advanced macro ...
          case(9) ! WORKING DIRECTORY and possibly MACRO
             call ui_set_working_directory(cline,last,kom3)
!.................................................................
          case(10) ! HELP_POPUP_OFF
             call gparcdx('Turn off popup help? ',cline,last,&
                  1,ch1,'Y','?Set adv help popup')
             if(ch1.eq.'Y') then
                ochelp%htmlhelp=.FALSE.
                htmlhelp=.FALSE.
             else
                htmlhelp=.TRUE.
                string=browser
                call gparcdx('Browser including full path ',&
                     cline,last,1,browser,string,'?Set adv help popup')
                string=htmlfile
                call gparcdx('HTML help file including full path ',&
                     cline,last,1,htmlfile,string,'?Set adv help popup')
                call init_help(browser,htmlfile)
                if(.not.ochelp%htmlhelp) write(kou,*)&
                     'Error initiating html help'
             endif
!.................................................................
          case(11) ! SET ADVANCED EEC_METHOD
             call gparcdx('Turn on equi-entropy criterion (EEC)?',&
                  cline,last,1,ch1,'Y','?Set adv EEC')
             if(ch1.eq.'Y' .or. ch1.eq.'y') then
!check if there is a phase with liquid but set!!
                anyliq: do j4=1,noph()
                   if(test_phase_status_bit(j4,PHLIQ)) exit anyliq
                enddo anyliq
! if we have not found any liquid j4>noph() here !!
                if(j4.gt.noph()) then
                   write(kou,*)'No liquid phase! Set bit 10 of liquid phase'
                   goto 100
                endif
                call gparrdx('Low T limit (min 10)?',cline,last,xxx,1.0D3,&
                     '?Set adv EEC')
                if(xxx.gt.1.0D1) then
! set_eec_check is in minimizer/matsmin.F90
                   globaldata%sysreal(1)=xxx
                endif
             else
                write(*,*)'EEC method for solids turned off as answer not Y'
                globaldata%sysreal(1)=zero
             endif
!.................................................................
          case(12) ! SET ADVANCED LEVEL
             call gparcdx('I am an beginner of OC: ',cline,last,1,ch1,'N',&
                  '?Set adv level')
             if(ch1.eq.'Y') then
                globaldata%status=ibset(globaldata%status,1)
                write(*,*)'Bon courage!'
             else
                call gparcdx('I am an expert of OC: ',cline,last,1,ch1,'N',&
                     '?Set adv level')
                if(ch1.eq.'Y') then
                   globaldata%status=ibset(globaldata%status,2)
                   write(*,*)'Felicitations!'
                else
                   write(*,*)'Sorry, not yet'
                endif
             endif
!.................................................................
          case(13) ! NO_MACRO_STOP on/off
             call gparcdx('Ignore macro @&: ',cline,last,1,ch1,'Y',&
                     '?Set adv no-macro-stop')
! iox(8) is declared in metlib4
             if(ch1.ne.'Y') then
                iox(8)=0
             else
                iox(8)=1
             endif
!.................................................................
          case(14) ! PROTECTION
             call gparrdx('Code',cline,last,proda,zero,'?Set adv protect')
             call gparrdx('Privilege',cline,last,privilege,zero,&
                  '?Set adv protect')
!.................................................................
          case(15) ! IGNORE_MACRE_ERRORS normally a macro error returns inter
             continue
!.................................................................
          case(16) ! SET ADVANCED XTDB_DEFAULTS
             dummy=lowtdef
             call gparcdx('Default low T limit',cline,last,1,lowtdef,dummy,&
                  '?Set adv xtdb')
             unary1991=.FALSE.
             if(lowtdef.eq.'298.15 ') unary1991=.TRUE.
             dummy=hightdef
             call gparcdx('Default high T limit',cline,last,1,hightdef,dummy,&
                  '?Set adv xtdb')
             ch1='N'
             call gparcdx('Include model descriptions (Y to include)',&
                  cline,last,1,chz,ch1,'?Set adv xtdb')
             if(chz.eq.'Y') then
                includemodels=.TRUE.
             else
                includemodels=.FALSE.
             endif
             dummy=eldef
             call gparcdx('Default elements (NONE to remove)',&
                  cline,last,1,eldef,dummy,'?Set adv xtdb')
             if(dummy(1:5).eq.'NONE ') eldef=' '
             dummy=bibrefdef
! argument 4 equal to 5 means the whole line is read
             call gparcdx('Default bibiographic reference (One line!)',&
                  cline,last,5,bibrefdef,dummy,'?Set adv xtdb')
!.................................................................
          case(17) ! SET ADVANCED_DEBUG
             call ui_set_mqmqa_debug(cline,last)
!.................................................................
          case(18) ! unused
             continue
          end select setadvanced
!-----------------------------------------------------------
! this is part of SET, 4 option removed
       case(4) 
          write(*,*)'Unused'
!-----------------------------------------------------------
! end of macro excution (can be nested)
       case(5) ! set INTERACTIVE
          call macend(cline,last,logok)  
! if this was the startupmacro set it false and possibly read an inline macro ..
! NOTE a startup macro can call other macros ...
          if(kiu.eq.kiud) startupmacro=.false.
!-----------------------------------------------------------
       case(6) ! set REFERENCE_STATE
          call ui_set_reference_state(cline,last,ceq)
!-----------------------------------------------------------
       case(7) ! quit set
          goto 100
!-----------------------------------------------------------
       case(8) ! set ECHO
          call gparcdx('On?',cline,last,1,ch1,'Y','?Set echo')
          if(ch1.eq.'Y' .or. ch1.eq.'y') then
             j4=1
          else
             j4=0
          endif
          call set_echo(j4)
!-----------------------------------------------------------
       case(9) ! set PHASE with subcommands (constitution, status)
          call gparcx('Phase name: ',cline,last,1,name1,' ',&
               '?Set phase')
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) then
             if(name1(1:2).eq.'* ') then
                iph=-1
                gx%bmperr=0
             else
                goto 990
             endif
          endif
! subsubbmenue of SET PHASE ...
          kom3=submenu('Set for phase what?',cline,last,csetph,nsetph,2,&
               '?TOPHLP')
          setphase: SELECT CASE(kom3)
!............................................................
          CASE DEFAULT
             write(kou,*)'Set phase status subcommand error'
             goto 100
!............................................................
          case(1) ! quit
             continue
!............................................................
          case(2) ! SET PHASE STATUS <phase> <status>
             call ui_set_phase_status(iph,ics,cline,last,ceq)
!............................................................
          case(3,4) !set phase default_constit wildcard allowed, also AMOUNT
             if(kom3.eq.3) then
! set default constituntion of phase
                call ask_default_constitution(cline,last,iph,ics,ceq)
             else
! set phase amount
                call gparrdx('Amount: ',cline,last,xxx,zero,&
                     '?Set phase constitution')
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
             kom4=submenu('Set which bit?',cline,last,csetphbits,nsetphbits,9,&
                  '?TOPHLP')
             phasebit: SELECT CASE(kom4)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(1) ! FCC_PERMUTATIONS FORD
! if check returns .true. it is not allowed to set FORD or BORD
                if(check_minimal_ford(lokph)) goto 100
                write(*,*)' *** WARNING: Depreceated command, use AMEND PHASE'
                call set_phase_status_bit(lokph,PHFORD)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(2) ! BCC_PERMUTATIONS BORD
                if(check_minimal_ford(lokph)) goto 100
                write(*,*)' *** WARNING: Depreceated command, use AMEND PHASE'
                call set_phase_status_bit(lokph,PHBORD)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(3) ! IONIC_LIQUID_MDL this may require tests and 
! other bits changed ..
                write(*,*)' *** WARNING: set by enter phase <name> I2SL'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(4) ! AQUEOUS_MODEL   
                write(*,*)'Not implemented yet'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(5) ! QUASICHEMICAL   
                write(*,*)'Not implemented yet'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(6) ! FCC_CVM_TETRADRN
                write(*,*)'Not implemented yet'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(7) ! FACT_QUASICHEMCL
                write(*,*)'Not implemented yet'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(8) ! NO_AUTO_COMP_SET, do not create compsets automatically
                call set_phase_status_bit(lokph,PHNOCS)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(9) ! QUIT
                write(kou,*)'No other bits changed'
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(10) ! EXTRA_DENSE_GRID, this can be toggled ...
                if(test_phase_status_bit(lokph,PHXGRID)) then
                   write(kou,*)'Bit already set, is cleared'
                   call clear_phase_status_bit(lokph,PHXGRID)
                else
                   write(kou,*)'Extra gridpoints for this phase.'
                   call set_phase_status_bit(lokph,PHXGRID)
                endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             case(11) ! PHEECLIQ bit for EEC phase
                write(*,*)'Bit for set for EEC liquid'
                call set_phase_status_bit(lokph,PHEECLIQ)
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
          call ui_set_logfile(cline,last,logfil)
!-------------------------------------------------------------
       case(12) ! set weight
          call ui_set_weight(cline,last,ceq)
!-------------------------------------------------------------
! turn on/off global minimization, creating composition sets
! convergence limits, iterations, minimum constituent fraction, etc
       case(13) ! set NUMERIC_OPTIONS
          call ui_set_numeric_options(cline,last,ceq)
!CCI
!-------------------------------------------------------------
       case(14) ! set axis
          call ui_set_axis(cline,last,maxax,ceq)
!-------------------------------------------------------------
       case(15) ! set input amounts
          call set_input_amounts(cline,last,ceq)
!-------------------------
       case(16) ! SET VERBOSE
! This toggles verbose for all commands.
! it is always turned of fwhen a command is finished ...
          if(btest(globaldata%status,GSSILENT)) then
! turn off VERBOSE and turn on SILENT
             globaldata%status=ibclr(globaldata%status,GSSILENT)
             write(kou,3603)'off',globaldata%status
          else
! turn on VERBOSE
             globaldata%status=ibset(globaldata%status,GSSILENT)
             write(kou,3603)'on',globaldata%status,GSSILENT
          endif
3603      format('Silent is turned ',a,2x,z8,i5)
!-------------------------
       case(17) ! SET AS_START_EQUILIBRIUM
! the current set of condition sill be used as start equilibrium for map/step
! Calculate the equilibrium and ask for a direction.
          if(noofaxis.lt.2) then
             write(kou,*)'You must set two axis first'
             goto 100
          endif
          call calceq2(1,ceq)
          if(gx%bmperr.ne.0) goto 990
          call gparidx('Give an axis direction: ',cline,last,ndl,2,&
               '?Set as start equil')
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
          starteqs(noofstarteq)%p1=>neweq
          write(*,*)'A copy of current equilibrium linked as start eqilibrium'
!---------------------------------------------------------------
       case(18) ! SET BIT (all kinds of bits) just global implemented
!         ['EQUILIBRIUM     ','GLOBAL          ','PHASE           ',
          kom3=submenu('Set which status word?',cline,last,csetbit,nsetbit,2,&
               '?TOPHLP')
          setbit: SELECT CASE(kom3)
          CASE DEFAULT
             write(kou,*)'SET BIT subcommand error'
!................................................................
          case(1) ! equilibrium status word
             call ui_set_bit_equilibrium(cline,last,ceq)
!................................................................
! maybe change order of questions, maybe check name exits etc ....
          case(2) ! global status word
             call ui_set_bit_global(cline,last)
!....................................................
          case(3) ! set bit phase ...
             write(*,*)'Please use set phase ... bit '
          end select setbit
!------------------------------------------------------------------
       case(19) ! set optcoeff_variabl, 0 to 99
          call ui_set_optcoeff_variable(cline,last)
!---------------------------------------------------------
       case(20) ! set optcoeff_scaled
          write(*,*)'Not implemeneted yet'
! zero the relative standard deviation
!----------------------------------------------------------
       case(21) ! set lmdif_accuracy, always propose the default!
          optacc=1.0D-3
          call gparrdx('LMDIF accuracy: ',cline,last,xxx,optacc,&
               '?Set optimizer conditions')
          write(kou,'("LMDIF accuracy set to ",1pe12.4)')xxx
          optacc=xxx
!----------------------------------------------------------
       case(22) ! set range_experimental_equilibria
          call ui_set_range_exper_equil(cline,last)
!---------------------------------------------------------
       case(23) ! set optcoeff_fixed
          call ui_set_optcoeff_fixed(cline,last)
!---------------------------------------------------------
       case(24) ! SET SYSTEM_VARIABLE
          write(kou,3733)
3733      format('Variable 2 is frequency of stability checks during step/map')
          call gparidx('System variable: ',cline,last,ll,0,&
               '?Set system variable')
          if(ll.eq.2) then
! sysparam(2) used during STEP/MAP often to check if equilibrium is stable
             call gparidx('System variable value: ',cline,last,j4,0,&
                  '?Set system variable')
             globaldata%sysparam(ll)=j4
          else
             write(*,*)'No other variable can be changed'
          endif
!-----------------------------------------------------------
       case(25) ! SET INITIAL_T_AND_P start values?, NOT CONDITIONS!!
          write(kou,3750)ceq%tpval
3750      format(/'NOTE: these are only local values, not conditions',&
               2(1pe12.4)/)
          call gparrdx('New value of T: ',cline,last,xxx,1.0D3,&
               'Set initial TP')
          if(buperr.ne.0) goto 100
          ceq%tpval(1)=xxx
          call gparrdx('New value of P: ',cline,last,xxx,1.0D5,&
               '?Set initial TP')
          if(buperr.ne.0) goto 100
          ceq%tpval(2)=xxx
!----------------------------------------------------------
!CCI
       case(26) ! SET LINEAR_SYSTEM
!  Splitsolver? is it really used ?
          indexSplitsolver = default_splitsolver
          call gparidx('Do you allow splitting of linear system? (1=Y; 0=N): ',&
               cline,last, indexSplitsolver,0,'?Set LINEAR_SYSTEM')
          if((indexSplitsolver.eq.0).or.(indexSplitsolver.eq.1)) then
            ceq%splitsolver = indexSplitsolver
          else
            ceq%splitsolver = default_splitsolver
          endif
!------------ Pre-conditioner?
          indexPrecond = default_precondsolver
          call gparidx('Enable the using Jacobi preconditioner? (1=Y; 0=N) : ',&
               cline,last, indexPrecond,0,'?Set LINEAR_SYSTEM')
          if((indexPrecond.eq.0).or.(indexPrecond.eq.1)) then
            ceq%precondsolver = indexPrecond
          else
            ceq%precondsolver = default_precondsolver
          endif
!------------ Scale the change of phase amount?
          typeOfChange = default_typechangephaseamount
          write(kou,555)
555       format('Scaleing changes in phase amount with total number',//&
               'of atoms?')
          call gparidx('Scaling (2=max, 1=sum; 0=one):',&
               cline,last, typeOfChange,0,'?Set LINEAR_SYSTEM')
          if((typeOfChange.eq.0).or.(typeOfChange.eq.1).or.&
               (typeOfChange.eq.2)) then
             ceq%type_change_phase_amount = typeOfChange
          else
             ceq%type_change_phase_amount = default_typechangephaseamount
          endif
!CCI
!------------------------- 
       case(27) ! SET GRID_GENERATOR, no choice available
          continue
       END SELECT set
!=================================================================
! ENTER with subcommand for element, species etc
!=================================================================
!         ['TPFUN_SYMBOL    ','ELEMENT         ','SPECIES         ',&
!         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
!         'CONSTITUTION    ','EXPERIMENT      ','QUIT            ',&
!         'EQUILIBRIUM     ','SYMBOL          ','OPTIMIZE_COEFF  ',&
!         'COPY_OF_EQUILIB ','COMMENT         ','MANY_EQUILIBRIA ',&
!         'MATERIAL        ','PLOT_DATA       ','GNUPLOT_TERMINAL',&
!         '                ','                ','                ']
    CASE(4)
! disable continue assessment optimization (not reelevant)
       kom2=submenu(cbas(kom),cline,last,center,ncent,11,'?TOPHLP')
       enter: SELECT CASE(kom2)
       CASE DEFAULT
          write(kou,*)'Enter subcommand error'
!---------------------------------------------------------------
! maybe change order of questions, maybe check name exits etc ....
       CASE(1) ! enter TPFUN symbol (constants, functions, tables)
          call ui_enter_tpfun(cline,last)
!---------------------------------------------------------------
       case(2) ! enter element
          if(.not.allowenter(1)) then
             gx%bmperr=4125
             goto 990
          endif
          call gparcx('Element symbol: ',cline,last,1,elsym,' ',&
               '?Enter element')
          if(buperr.ne.0) goto 990
          call capson(elsym)
          if(.not.(elsym(1:1).ge.'A' .and. elsym(1:1).le.'Z')) then
             write(*,*)'An element symbol must start with a letter A-Z'
             goto 100
          endif
          call gparcdx('Element full name: ',cline,last,1,name1,elsym,&
               '?Enter element')
          call gparcdx('Element reference phase: ',cline,last,1,&
               name2,'SER ','?Enter element')
          call gparrdx('Element mass (g/mol): ',cline,last,mass,one,&
               '?Enter element')
          if(buperr.ne.0) goto 990
          call gparrdx('Element H298-H0: ',cline,last,h298,zero,&
               '?Enter element')
          if(buperr.ne.0) goto 990
          call gparrdx('Element S298: ',cline,last,s298,one,'?Enter element')
          if(buperr.ne.0) goto 990
          call store_element(elsym,name1,name2,mass,h298,s298)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(3) ! enter species
! Allow entering species even if there are phases entered
! needed for the MQMQA model
!>>> There may be problems with MQMQA quad names such as Fe,Mn/SI1/4O,Al2/3O ?
          call gparcx('Species symbol: ',cline,last,1,name1,' ',&
               '?Enter species')
! check if it is an MQMQA quad
          name2=name1
          call capson(name2)
          iz=index(name2,'/')
          if(iz.gt.0 .and. &
               (name2(iz+1:iz+1).ge.'A' .and. name2(iz+1:iz+1).le.'Z')) then
! A MQMQA species has a letter after the /
! In MQMQA species, the two sublattices indicated by /
! but we musr also separate species in same, for example: Fe,Mn/Si1/4O,Al2/3O
             mqmqass=' '
             call gparcx('MQMQA specification: ',cline,last,5,&
                  mqmqass,' ','?Enter MQMQA species')
! typically mqmqass is Na/Cl 6.0 6.0 2.4 (last values the FNN/SNN ratio
! or Mg,Na/Cl  6.00 3.00 3.00 etc. as in Solgsmax DAT files
! check species used before quad numbers must already be entered!
! mqmqanend should be negative at first call
             write(*,575)trim(name1),trim(mqmqass),mqmqanend
575          format('Call mqmqa_species: "',a,'" "',a,'" ',i5)
             call mqmqa_species(name1,mqmqass,mqmqanend)
          else
! A species with no / or an ionic / followed by + - or number
             call gparcx('Species stoichiometry: ',cline,last,1,name2,' ',&
                  '?Enter species')
             call decode_stoik(name2,noelx,ellist,stoik)
             if(gx%bmperr.ne.0) goto 990
! all species must be entered
             call enter_species(name1,noelx,ellist,stoik)
          endif
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
             write(kou,*)'You must enter a phase before'
             goto 100
          endif
! the last 0 means enter
          call enter_parameter_interactivly(cline,last,0)
! Strange things may happen when entering parameters interactively 
! This was due to an error in tpfun package ... not yet fixed ... ??
          call change_optcoeff(-1,zero)
          do j4=1,notpf()
             call eval_tpfun(j4,ceq%tpval,val,ceq%eq_tpres)
             if(gx%bmperr.gt.0) goto 990
          enddo
          call change_optcoeff(-1,zero)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(6) ! enter bibliography
          call enter_bibliography_interactivly(cline,last,0,j4)
          if(gx%bmperr.ne.0) goto 990
          write(kou,*)'Bibliography number is ',j4
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
          call gparcdx('Name: ',cline,last,1,text,quest,'?Enter equilibrium')
          if(buperr.ne.0) goto 100
          call enter_equilibrium(text,ieq)
          if(gx%bmperr.ne.0) goto 990
! by default also select this equilibrium
          write(kou,303)ieq,trim(text)
303       format('Equilibrium number is ',i3,', name: ',a)
          call gparcdx('Select this equilibrium: ',cline,last,1,ch1,'Y',&
               '?Enter equilibrium')
          if(yeschk(ch1)) then
             call selecteq(ieq,ceq)
! COPY current values of entered symbols from first equilibrium
             ceq%svfunres=firsteq%svfunres
          endif
!---------------------------------------------------------------
       case(11) ! enter symbol (for state variables expressions)
! several questions asked inside this call
          call enter_svfun(cline,last,ceq)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
! enter optimizing coefficients called A00 to A99 (or whatever set as max)
       case(12)
          call ui_enter_optimize_coeff(cline,last)
!---------------------------------------------------------------
! enter copy_of equilibrium (for test!)
       case(13)
! Check if there is any phases, otherwise not allowed
          if(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'Not allowed unless you have data!'
             goto 100
          endif
          call gparcx('Name of new equilibrium: ',cline,last,1,text,' ',&
               '?Enter copyof')
          if(buperr.ne.0) goto 100
          if(text(1:1).eq.' ') then
             write(*,*)'You must specify a unique name'
             goto 100
          endif
          call copy_equilibrium(neweq,text,ceq)
          if(gx%bmperr.ne.0) goto 990
          write(kou,*)'New equilibrium no: ',neweq%eqno
!---------------------------------------------------------------
! enter COMMENT for current equilibrium
       case(14)
          write(*,*)'Current equilibrium name: ',ceq%eqname
          call gparcx('One line text: ',cline,last,5,text,' ',&
               '?Enter comment')
          ceq%comment=text
!---------------------------------------------------------------
! enter MANY_EQUILIBRIA
! The plotdataunit array should be zero at first call, then the unit is opened
! (if there are any plot_data commands).  It will remain open until
! a set range command is given
       case(15)
          call enter_many_equil(cline,last,plotdataunit)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
! enter MATERIAL
! ask for database, then major element, mass/mole fraction of elements
! read the database; jump possibly to SCHEIL/STEP calculation 
! or simply ask for T and calculate equilibrium; 
       case(16)
          write(*,*)'MM Not implemented yet'
          goto 100
! keep for future ambitions ...
          call enter_material(cline,last,nv,xknown,ceq)
          if(gx%bmperr.ne.0) goto 990
          xxy=firsteq%tpval(1)
          call gparrdx('Temperature ',cline,last,xxx,xxy,'?Enter material')
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
          call gparidx('Dataset number:',cline,last,i1,1,'?Enter plot data')
! here only the normal plotdata units 1 to 9 are legal
          if(i1.gt.0 .and. i1.lt.10) then
             if(plotdataunit(i1).lt.10) then
                write(kou,*)'No plotdata file for this dataset'
                goto 100
             endif
             call gparrdx('X coordinate:',cline,last,xxx,zero,&
                  '?Enter plot data')
             call gparrdx('Y coordinate:',cline,last,xxy,one,&
                  '?Enter plot data')
             call gparidx('Symbol:',cline,last,i2,1,&
                  '?Enter plot data')
             write(plotdataunit(i1),171)i1,xxx,xxy,i2
171          format(i3,2(1pe14.6),i5,' have a nice day')
          else
             write(kou,*)'No plotdata file for dataset ',i1
          endif
!---------------------------------------------------------------
! ENTER GNUPLOT_TERMINAL
       case(18)
          call ui_enter_gnuplot_terminal(cline,last)
!----------------------------------------------------------------
! enter unused
       case(19,20,21)
          write(*,*)'Not implemeneted yet'
       END SELECT enter
!=================================================================
! exit
    CASE(5)
!=================================================================
       call gparcdx('Are you sure?',cline,last,1,ch1,'N','?Exit')
       if(ch1.eq.'y' .or. ch1.eq.'Y') then
          if(logfil.gt.0) then
             write(logfil,*)'set interactive'
          endif
          call openlogfile(' ',' ',-1)
          stop 'Ha en bra dag'
       endif
!=================================================================
! LIST and SHOW with subcommands
!=================================================================
!         ['DATA            ','SHORT           ','PHASE           ',&
!         'STATE_VARIABLES ','BIBLIOGRAPHY    ','MODEL_PARAM_ID  ',&
!         'AXIS            ','TPFUN_SYMBOLS   ','QUIT            ',&
!         'PARAMETER       ','EQUILIBRIA      ','RESULTS         ',&
!         'CONDITIONS      ','SYMBOLS         ','LINE_EQUILIBRIA ',&
!         'OPTIMIZATION    ','MODEL_PARAM_VAL ','ERROR_MESSAGE   ',&
!         'ACTIVE_EQUILIBR ','ELEMENTS        ','EXCELL_CSV_FILE ',&
!         'MQMQA_SPECIAL   ','ESTIMAT_ACCURACY','WORKING_DIRECTOR',&
!         '                ','                ','                ']
! SHOW is main cammand 25
    CASE(6,25) ! LIST and SHOW
       if(kom.eq.25) then
! SHOW is the same as LIST STATE_VARIABLES including also CALC SYMBOL !!
          kom2=4
       else
! default for LIST is RESULT, number 12
          kom2=submenu(cbas(kom),cline,last,clist,nclist,12,'?TOPHLP')
          if(kom2.le.0) goto 100
       endif
       lut=optionsset%lut
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
!    character (len=16), dimension(nlform) :: llform=&
!        ['SCREEN          ','                ','MACRO           ',&
!         '                ','                ','                ']
          if(globaldata%encrypted.ne.0) then
             write(kou,*)'Illegal for encrypted databases'
             goto 100
          endif
          kom3=submenu('Output format for data?',cline,last,llform,nlform,1,&
               '?TOPHLP')
          if(kom3.eq.1) then
! list on screen
             call list_many_formats(cline,last,kom3,kou)
             if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
                write(kou,*)bmperrmess(gx%bmperr)
             elseif(gx%bmperr.ne.0) then
                write(kou,*)'Error code ',gx%bmperr
             endif
          elseif(kom3.eq.3) then
             write(*,*)'Output in MACRO format not yet implemented'
          else
! TDB format does not work here, use SAVE
             write(*,*)'Use SAVE to list with other formats than SCREEN'
          endif
!-----------------------------------------------------------
       case(2) ! list short with status bits
          if(kom2.eq.20) then
             ch1='C'
          else
! note D is a hidden option including the status bits
             call gparcdx('Option (A/C/M/P)',cline,last,1,ch1,chshort,&
                  '?List short')
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
          elseif(ch1.eq.'D') then
! just the phases
! P phases sorted: stable/ unstable in driving force order/ dormant the same
             chshort='P'
             call list_sorted_phases(lut,1,ceq)
             if(btest(ceq%status,EQFAIL)) write(lut,6305)'above'
!....................................................................
          elseif(ch1.eq.'P') then
! just the phases without status bits
! P phases sorted: stable/ unstable in driving force order/ dormant the same
             chshort='P'
             call list_sorted_phases(lut,0,ceq)
             if(btest(ceq%status,EQFAIL)) write(lut,6305)'above'
6305         format(/' *** The results ',a,&
                  ' are not a valid equilibrium as last calculation failed'/)
          elseif(ch1.eq.'C') then
!....................................................................
! global values and the chemical potentials
             chshort='C'
             write(kou,*)
             call list_global_results(lut,ceq)
             write(lut,6303)'Some data for components ...............'
6303         format(/a,20('.'),':')
             j4=1
             if(listresopt.ge.4 .and. listresopt.le.7) then
                j4=2
             endif
             call list_components_result(lut,j4,ceq)
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
       case(3) ! list phase subcommands SOME CHANGES HERE
!        ['DATA            ','CONSTITUTION    ','MODEL           ',&
!         'ASYMMETRIES     ','TERNARY_ASYM    ','                ']
!
          call gparcx('Phase name: ',cline,last,1,name1,' ','?List phase')
          if(buperr.ne.0) goto 990
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
          kom3=submenu('List what for phase?',cline,last,clph,nclph,2,&
               '?TOPHLP')
          listphase: SELECT CASE(kom3)
!...............................................................
          CASE DEFAULT
             write(kou,*)'list phase ... subcommand error'
!...............................................................
          CASE(1) ! list phase data
             call list_phase_data(iph,' ',lut)
!...............................................................
! list phase constitution
          case(2) ! list phase ... constitution
             write(lut,6051)version,ceq%eqno,ceq%eqname
6051         format(/'OC version',a,' equilibrium: ',i3,', ',&
                  a,3x,a4,'.',a2,'.',a2)
             mode=110
             once=.TRUE.
             call list_phase_results(iph,ics,mode,lut,once,ceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Last equilibrium calculation failed'
                goto 990
             endif
!...............................................................
          case(3) ! list phase ... model (including disordered fractions)
             write(kou,6070)'For ',ceq%eqno,ceq%eqname
6070      format(a,'equilibrium: ',i3,', ',a)
             call list_phase_model(iph,ics,lut,' ',ceq)
!...............................................................
          case(4) ! LIST PHASE ... ASYMMETRIES
! clumsy way check this is the MQMQA phase
             lokcs=phasetuple(iph)%lokvares
             parres=>ceq%phase_varres(lokcs)
! Claude 2026.11.07: separate listing for several anions
             if(mqmqa_data%ncon2.ge.2) then
                call list_mqmqa_var2(parres)
                goto 100
             endif
             if(.not.allocated(parres%mqmqaf%y_ik)) then
                write(*,*)'MM This is not the MQMQA liquid phase'
                goto 990
             endif
! parres is pointer to gtp_phase_varres for the MQMQA phase
             call list_mqmqa_variables(parres)
!...............................................................
          case(5) ! LIST PHASE TERNARY_ASYM
             call list_tersys
             goto 100
          END SELECT listphase
!------------------------------
! THIS IS ALSO THE SHOW command and list model-parameter-value case(17) of LIST
! SHOW STATE VARIABLE VALUE
       case(4,17)  ! list state_variable or model_parameter_value, or SHOW
          call ui_list_state_variables(cline,last,kom,kom2,lut,ceq)
!-----------------------------------------------------------
       case(5) ! list data bibliography
          call gparcdx('Bibliographic id:',cline,last,1,name1,'ALL',&
               '?List biblio')
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
          do iax=1,noofaxis
             jp=1
             call get_one_condition(jp,text,axarr(iax)%seqz,ceq)
             if(gx%bmperr.ne.0) then
                write(kou,*)'PMON: Condition sequential index: ',&
                     iax,axarr(iax)%seqz
                goto 990
             endif
! we just want the expression, remove the value including the = sign
             jp=index(text,'=')
             text(jp:)=' '
             write(lut,6132)iax,text(1:24),&
                  axarr(iax)%axmin,axarr(iax)%axmax,axarr(iax)%axinc
6132         format(i2,2x,a,3(1pe12.4))
          enddo
!-----------------------------------------------------------
       case(8) ! list tpfun symbol
          call gparcdx('name: ',cline,last,5,name1,'*','?List TPfun')
          lrot=0
          iel=index(name1,'*')             
          if(iel.gt.1) name1(iel:)=' '
          if(name1(1:1).ne.'*') then
6140         continue
             call find_tpfun_by_name(name1,lrot)
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
       case(9) ! list quit, removed 2026.11.07 (QUIT and QUADS both start with Q)
!------------------------------------------------------------
       case(10) ! list parameter for a phase (just one). Last 1 means list
          call enter_parameter_interactivly(cline,last,1)
!-----------------------------------------------------------
       case(11,19) ! list EQUILIBRIA and list ACTIVE_EQUILIBRIA (not result)
          call ui_list_equilibria(kom2,lut,ceq)
!------------------------------
       case(12) ! LIST RESULTS (maybe also LIST ESTIMATED_ACCURA?)
          call ui_list_results(cline,last,version,listresopt,lrodef,lut,ceq)
!------------------------------
       case(13) ! list conditions
          write(kou,6070)'Conditions for ',ceq%eqno,ceq%eqname
          call list_conditions(lut,ceq)
!------------------------------
       case(14) ! list symbols (state variable functions, not TP funs)
          call list_all_svfun(lut,ceq)
!------------------------------
! list line_equilibria, (line-equilibria) of calculated and stored equilibria
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
!         'EXPERIMENTS     ','CORRELATION_MTRX','MQMQA_QUAD      ']
       case(16)
          if(.not.allocated(firstash%coeffstate)) then
             write(kou,*)'No listing as no optimizing parameters'
             goto 100
          endif
          call date_and_time(optres,name1)
          kom2=submenu('List ',cline,last,optopt,noptopt,1,'?TOPHLP')
! allow output file
          lut=optionsset%lut
! if errs not allocated no optimization made
          if(allocated(errs)) then
! trying to avoid segmentation fault when errs destryed by PLOT with APPEND
             if(size(errs).ne.mexp) then
                write(*,*)'Allocation error of "errs"',size(errs),mexp
          write(*,*)' **** Warning, datastructure corrupted, save what you can'
                goto 100
             endif
             write(lut,600)optres(1:4),optres(5:6),optres(7:8),&
                  name1(1:2),name1(3:4),err0(3)
600          format(/'Optimization results at ',a4,'.',a2,'.',a2,&
                  ':',a2,'h',a2,', normalized sum of error: ',1pe12.4)
          else
             write(*,*)'No current optimization'
          endif
          listopt: SELECT CASE(kom2)
!..........................................................
             case DEFAULT
                write(kou,*)'No such option'
!...........................................................
! list optimization short
             case(1) ! short
                if(allocated(errs)) then
                   if(size(errs).eq.mexp) then
! in matsmin
                      call listoptshort(lut,mexp,nvcoeff,errs)
                   else
! After PLOT ... with APPEND of experimental data "errs" seems destroyed?? 
                      write(kou,*)'Allocation error: ',mexp,size(errs)
                      deallocate(errs)
                   endif
                endif
! in gtp3C
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
563                format(6(1pe12.4))
                else
                   write(*,*)'No correlation matrix calculated'
                endif
!...........................................................
! list optimization RSD (according to OC and TC)
             case(9)
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
!------------------------------
! list error message
          case(18)
             i2=4204
             call gparidx('Error code: ',cline,last,i1,i2,'?List error msg')
             if(i1.ge.4000 .and. i1.le.nooferm) then
                write(kou,4999)i1,bmperrmess(i1)
4999            format('The error code ',i4,', means: '/a)
             else
                write(kou,*)'Not a standard OC error message'
             endif
!------------------------------
! list ?? nonzero_equilibria/active-equil merged with list equilibra, case 11
!------------------------------
! list elements
          case(20)
             call list_all_elements(kou)
!------------------------------
! list Excell CSV file, code copied from PLOT
          case(21)
             call ui_list_excel_csv(cline,last,version,ceq)
!-------------------------------------------------------------------
! list MQMQA_SPECIALS phase ..
!    character (len=16), dimension(mqmqacc) :: mqmqalist=&
!        ['QUADS           ','ASYMMETRIES     ','DEBUG           ',&
!         'EXCESS          ','INTERNAL        ','                ']
          case(22)
! depreciated command ??
! trying to avoid segmentation fault when errs destryed by PLOT with APPEND
             call gparcx('Phase name: ',cline,last,1,name1,' ','?List phase')
             if(buperr.ne.0) goto 990
! special to debug read database
             if(name1.eq.'DEBUG') then
                mqmqtdb=.true.
                goto 100
             endif
             call find_phase_by_name(name1,iph,ics)
             if(gx%bmperr.ne.0) goto 990
! A stupid way to find lokvares ... (or lokcs, I have forgooten which)
             if(.not.allocated(mqmqa_data%contyp)) then
! if this not allocated there is no MQMQA data
                write(*,*)'No data for MQMQA phases'
                exit main
             endif
             lokcs=1
             do while(ceq%phase_varres(lokcs)%phlink.ne.iph)
                lokcs=lokcs+1
             enddo
!
             kom2=submenu('MQMQA special?',&
                  cline,last,mqmqalist,mqmqacc,2,'?TOPHLP')
!                                           ^ this 2 is the default answer
!..........................................................
             listphasemqmqa: select case(kom2)
             case DEFAULT
                write(kou,*)'No such option'
!...........................................................
! list mqmqa_speciai  <phase> ...
             case(1) ! Quads
                jquad=0
                write(*,1679)
1679            format('Specie Loc Name',16x,'Elements, stoichiometry, bonds')
                qlista: do i1=1,nosp()
                   call get_species_location(i1,loksp,name1)
                   if(gx%bmperr.ne.0) goto 990
                   if(index(name1,'-Q').le.0) cycle qlista
! this is a species that is a quad
                   jquad=jquad+1
                   call get_species_component_data(loksp,i2,iphl,stoik,&
                        xxx,xxy,ceq)
                   if(gx%bmperr.ne.0) goto 990
                   do j4=1,i2
! pick up element symbols
                      call get_element_data(iphl(j4),ellist(j4),name2,&
                           dummy,mass,h298,s298)
                   enddo
                   write(kou,1680)i1,loksp,name1,&
                        (ellist(j4),stoik(j4),j4=1,i2)
1680               format(3x,i3,i4,1x,a,1x,4(a2,F8.6,1x))
! this is in quad alphabetical order, not in alphabetical order of quad elements
! Claude 2026.11.07: mqmqa_quadbonds needs the constituent index, not jquad
!                   call mqmqa_quadbonds(jquad,quadbonds)
                   j4=mqmqa_quad_of_species(loksp)
                   if(j4.le.0) cycle qlista
                   call mqmqa_quadbonds(j4,quadbonds)
                   if(i2.eq.2) then
                      write(kou,1681)(quadbonds(j4),j4=1,3)
                   else
                      write(kou,1681)(quadbonds(j4),j4=1,i2)
                   endif
1681               format(26x,'bonds: ',4(F10.6,1x))
!    max length         8+25+4*11=33+44=77
                enddo qlista
!
                call quadprops
!
                write(kou,3001)size(mqmqa_data%cations)
3001            format(/'List of ',i3,' cations: *** to be fixed')
3002            format(15(a,', '))
                if(jquad.eq.0) write(kou,*)'No MQMQA quads found'
! maybe include listing of mqmqa_data%constoi(1..4,index)
!...........................................................
! list mqmqa_special <phase> ASYMMETRIES
! note: tersys, xquad, compvar are not linked from the phase!!!
             case(2)
! Claude 2026.11.07: several anions have their own listings
                if(mqmqa_data%ncon2.ge.2) then
                   call list_tersys
                   mqmqavar=>ceq%phase_varres(lokcs)
                   call list_mqmqa_var2(mqmqavar)
                   goto 4445
                endif
! list element names, numbers and quad indices, i1 set to number of quads
                call list_quads(i1)
! tersys is global data
                call list_tersys
! listing of fraction in alphbetical order
                write(kou,4123)mqmqa_data%nquad,&
                     (ceq%phase_varres(lokcs)%yfr(i1),i1=1,mqmqa_data%nquad)
! mqmqaf defined globally
                anoq: if(.not.&
                     allocated(ceq%phase_varres(lokcs)%mqmqaf%xquad)) then
                   write(*,*)'Quads not allocated'
! nothing more to list
                   goto 4445
                else
                   write(kou,4122)mqmqa_data%nquad,&
                        (ceq%phase_varres(lokcs)%mqmqaf%xquad(i1),&
                        i1=1,mqmqa_data%nquad)
4122               format('Fractions of ',i2,' in Quad order, x_ij:',/(12F6.3/))
4123               format('Fractions of ',i2,' in species OC alpha order:',/&
                        (12F6.3/))
! lokcs is MQMQA phase
                   mqmqavar=>ceq%phase_varres(lokcs)
! Now this incluse also y_i/k and x_ij equations and values
                   call list_mqmqa_variables(mqmqavar)
                   goto 4445
                endif anoq
4445            continue
!
!...........................................................
! list mqmqa_special <phase> DEBUG 
             case(3)
                write(*,*)'Use SET ADVANCED DEBUG to set debug level'
!...........................................................
! list mqmqa_special <phase> excess tree   I have no idea what this 
! list the constituents of the phase in the order they have in constitlink
! THIS IS AN EMERGY SUBROUTINE IN gtp3XQ NOT CONFORMING WITH THE STRUCTURE
             case(4)
                call get_phase_record(iph,lokph)
                call listpartree(lokph)

!...........................................................
! list mqmqa_special <phase> INTERNALS  AMEND_VARKAPPA
             case(5)
! CODE BELOW REPLACED by new code around line 1374 AMEND PHASE ... ASYMMETRY
! list quad and quad indices, i2 loops for all quads
                if(allocated(mqmqa_data%contyp)) then
                   write(*,'(/a)')'Listing mqmqa_data%contyp by quadprops'
                   write(*,*)'Missing call to quadprops'
1388               format('quad2compvar: ',20i3)
                   write(*,1388)mqmqa_data%quad2compvar
                else
                   write(*,*)'No phase with this property'
                endif
!...........................................................
! LIST mqmqa_special <phase> unused
             case(6)
                continue
             end SELECT listphasemqmqa
!------------------------------ end list mqmqa_specials
! LIST ESTIMATE_ACCURACY.  Additional calculations are made
! Eventually included in case(12) ???
          case(23) 
             if(btest(ceq%status,EQNOEQCAL) .or. btest(ceq%status,EQFAIL)) then
                write(kou,*)'You must calculate an equilibrium first'
                goto 100
             endif
             xxy=5.0
             call gparrdx('Estimated uncertainty in conditions (%): ',&
               cline,last,xxx,xxy,'?List confidence interval')
             i1=optionsset%lut
             if(i1.eq.0) i1=kou
             call calc_conf_interval(i1,xxx,ceq)
             if(gx%bmperr.ne.0) then
                ceq%status=ibset(ceq%status,EQFAIL)
                goto 990
             endif
!------------------------------
! list WORKING_DIR
          case(24)
             write(kou,1685)trim(workingdir)
1685         format('Current working directory is: ',a)
!------------------------------
! list QUADS, Claude 2026.11.07, AS_SPECIES: x_ABxy and species name,
! SHORT: the same 3 on each line with the quad index
          case(25)
             if(.not.allocated(mqmqa_data%contyp)) then
                write(kou,*)'No MQMQA phase'
                goto 100
             endif
             call gparcdx('How (AS_SPECIES or SHORT): ',cline,last,1,name1,&
                  'AS_SPECIES','?List quads')
             call capson(name1)
             if(name1(1:1).eq.'A') then
                call list_quads_as_species(lut,.FALSE.)
             elseif(name1(1:1).eq.'S') then
                call list_quads_as_species(lut,.TRUE.)
             else
                write(kou,*)'Only AS_SPECIES or SHORT implemented'
             endif
!------------------------------
! list ??
          case(26,27)
             write(kou,*)'Not implemented yet'
          end SELECT list
!=================================================================
! quit in the main list of select
       case(7)
!=================================================================
          if(cline(1:1).eq.'q') then
             call gparcdx('Are you sure?',cline,last,1,ch1,'N','?Quit')
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
!=================================================================
!        ['UNFORMATTED     ','TDB             ','QUIT            ',&
!         'DIRECT          ','XTDB            ','SELECTED_PHASES ']
       case(8)
! disable continue optimization
       if(noel().ne.0) then
          write(kou,*)'You already have data, read destroys your current data'
          write(kou,*)'You must give a NEW Y command to remove data first'
          goto 100
! all records must be removed and init_gtp is called.  This is fragile ...
       endif
       kom2=submenu(cbas(kom),cline,last,cread,ncread,2,'?TOPHLP')
       read: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          if(cline(len_trim(cline):len_trim(cline)).ne.'?') then
! This avoids error messages when ? is typed
             write(kou,*)'Read subcommand error: ',trim(cline)
          endif
!-----------------------------------------------------------
       case(1) ! read unformatted file created by SAVE
          if(ocufile(1:1).ne.' ') then
             text=ocufile
             call gparcdx('File name: ',cline,last,1,ocufile,text,&
                  '?Read unformatted')
          else
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
             ztyp=2
             call gparfilex('File name: ',cline,last,1,ocufile,' ',ztyp,&
                  '?Read unformatted')
          endif
          call gtpread(ocufile,text)
          if(gx%bmperr.ne.0) then
             ocufile=' '; goto 990
          endif
! This is written by the gtpread subroutine
! if there is an assessment record set nvcoeff ...
          if(allocated(firstash%coeffvalues)) then
             nvcoeff=0
             kl=size(firstash%coeffvalues)-1
             do j4=0,kl
                if(firstash%coeffstate(j4).ge.10) then
                   nvcoeff=nvcoeff+1
                endif
             enddo
             write(kou,3730)nvcoeff
3730         format('Number of variable coefficients are now ',i3)
          else
             write(*,*)'No coefficients allocated'
          endif
          if(allocated(firstash%eqlista)) then
             write(*,*)'There are experimental data'
          endif
!---------------------------------------------------------
       case(2,7) ! read TDB and read ENCRYPTED
          call ui_read_tdb(cline,last,kom2,tdbfile,ocbase)
!-----------------------------------------------------------
       case(3) ! read quit
          goto 100
!-----------------------------------------------------------
       case(4) ! read direct
          write(*,*)'Read direct not implemented yet'
!-----------------------------------------------------------
! read the new XTDB format for Calphad databases
       case(5) ! read XTDB 
          call ui_read_xtdb(cline,last,xtdbfile)
!-----------------------------------------------------------
       case(6) ! read SELECTED_PHASES
! part of read TDB          
          call ui_read_tdb(cline,last,kom2,tdbfile,ocbase)
!...............................................................
       case(8,9) ! read ?
       end SELECT read
!=================================================================
! SAVE in various formats (NOT MACRO and LATEX, use LIST DATA)
!=================================================================
! It is a bit inconsistent as one READ TDB but not SAVE TDB ...
!        ['TDB             ','                ','QUIT            ',&
!         'DIRECT          ','UNFORMATTED     ','XTDB            ']
    CASE(9)  ! SAVE
! default is 3, unformatted
       kom2=submenu(cbas(kom),cline,last,csave,ncsave,1,'?TOPHLP')
       if(kom2.le.0 .or. kom2.gt.ncsave) goto 100
!
       call date_and_time(optres,name1)
! optres(1:8) is year+month+day, name1(1:4) is hour and minutes
       model=' '//optres(1:4)//'.'//optres(5:6)//'.'//optres(7:8)//&
            ' '//name1(1:2)//'h'//name1(3:4)//' '
       save: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'save subcommand error'
!-----------------------------------------------------------
       case(1) ! save TDB, same as list data TDB
! format 1 is TDB, see list data ...
          if(globaldata%encrypted.ne.0) then
             write(kou,*)'Illegal for encrypted databases'
             goto 100
          endif
! gparfilex next to last argument is negative if 
          ztyp=-1
          call gparfilex('File name: ',cline,last,1,filename,text,ztyp,&
               '?Save TDB')
! Bosse do not understand ???
          kl=max(index(filename,'.dat '),index(filename,'.TDB '))
          if(kl.le.0) then
             kl=len_trim(filename)+1
             if(kl.eq.1) then
                write(*,*)'Too short file name'
                goto 100
             endif
! Bosse do not understand ???
          endif
! inside list_TDB_format
          call list_TDB_format(filename)
          if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
             write(kou,*)bmperrmess(gx%bmperr)
          elseif(gx%bmperr.ne.0) then
             write(kou,*)'Error code ',gx%bmperr
          endif
!-----------------------------------------------------------
       case(2) ! used to be SOLGAS no longer available
          continue
!-----------------------------------------------------------
       case(3) ! save quit, do nothing
          continue
!-----------------------------------------------------------
       case(4) ! save DIRECT
          write(*,*)'Not implemented'
!-----------------------------------------------------------
       case(5) ! save unformatted
          call ui_save_unformatted(cline,last,ocufile,model)
!-----------------------------------------------------------
       case(6) ! SAVE XTDB
          if(globaldata%encrypted.ne.0) then
             write(kou,*)'Illegal for encrypted databases'
             goto 100
          endif
! is there any data?
          if(noph().le.0) then
             write(kou,*)'There is no data to save'
             goto 100
          endif
! Ask for file name
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
          ztyp=-6
          call gparfilex('File name: ',cline,last,1,xtdbfile,' ',ztyp,&
               '?Save XTDB')
!
          zext='XTDB'
! this subrouine is in gtp3EX.F90
          call write_xtdbformat(xtdbfile,zext)
          if(gx%bmperr.ne.0) goto 990
       end SELECT save
!=================================================================
! help ... just list the commands
    case(10)
!=================================================================
       call q3helpx(cline,last,cbas,ncbas)
!=================================================================
! subcommands to INFORMATION ... very little implemented
!=================================================================
!        ['ELEMENTS         ','SPECIES         ','PHASES          ',&
!         'QUIT-INFO        ','COMPOSITION_SET ','EQUILIBRIUM     ',&
!         'HELP_SYSTEM      ','CONDITIONS      ','DATABASES       ',&
!         'CHANGES          ','PHASE_DIAGRAM   ','PROPERTY_DIAGRAM',&
!         'STATE_VARIABLES  ','                ','                ']
    case(11)
! initial default is CHANGES
       iz=10
! return here until quit
207    continue
       kom2=submenu('Topic?',cline,last,cinf,ninf,iz,'?TOPHLP')
! change default to quit
       iz=4
       information: select case(kom2)
!-------------------------------------------------------
       CASE DEFAULT
          write(*,*)'Information subcommand error'
!--------------------------------------------------------
! INFO elements
       case(1)
          write(kou,210)
210       format('The elements are those from the periodic chart.'/&
               'Normally the components are the same as the elemets but',&
               ' the user',/'can define any orthogonal set of species as',&
               ' components.')
          call q4help('Info elements',jp)
!--------------------------------------------------------
! info species
       case(2)
          write(kou,211)
211       format('Species are molecular like aggregates of elements with',&
               ' fixed stoichiometry.',/'The elements are the simplest',&
               ' species.'/'The constituents of a phase are a subset of',&
               ' the species.')
          call q4help('Info species',jp)
!--------------------------------------------------------
! info phases
       case(3)
          call q4help('Info phases',jp)
!--------------------------------------------------------
! quit, we exit to top level here as there is a loop
       case(4)
          goto 100
!--------------------------------------------------------
! info composition set
       case(5)
          call q4help('Info compset',jp)
!--------------------------------------------------------
! info equilibrium
       case(6)
          call q4help('Info equilibrium',jp)
!--------------------------------------------------------
! INFO help system
       case(7) ! none
          call q4help('Info helpsystem',jp)
!--------------------------------------------------------
! INFO conditions
       case(8) ! none
          call q4help('Info conditions',jp)
!--------------------------------------------------------
! INFO databases
       case(9) ! none
          call q4help('Info databases',jp)
!--------------------------------------------------------
! INFO changes
       case(10)
          write(kou,'(a/)')'Writing from "OCHOME/changes.txt"'
          open(31,file=trim(OCHOME)//'/changes.txt ',access='sequential',&
               err=990,iostat=buperr)
          changes: do while(.TRUE.)
             do i1=1,40
                read(31,17,end=244,err=990)line
                write(kou,17)trim(line)
17              format(a)
             enddo
             write(kou,*)'Press return to continue, q to quit'
             read(kiu,17)ch1
             if(ch1.eq.'q' .or. ch1.eq.'Q') exit changes
          enddo changes
244       close(31)
!--------------------------------------------------------
! INFO phase diagram
       case(11) 
          call q4help('Info phasediagram',jp)
!--------------------------------------------------------
! INFO property diagram
       case(12) 
          call q4help('Info propertydiagram',jp)
!--------------------------------------------------------
! INFO statevariables
       case(13) 
          call q4help('Info statevariables',jp)
!--------------------------------------------------------
! INFO missing
       case(14,15) ! none
!--------------------------------------------------------
       end select information
! loop unitil quit
       goto 207
!=================================================================
! back / goto, return to calling (main) program
    case(12)
!=================================================================
       write(*,*)'Welcome back!'
       return
!=================================================================
! NEW command, same as reinitiate
    case(13) ! NEW
       call ui_new(cline,last,defcp,reinit)
       if(reinit) goto 20
!=================================================================
! MACRO begin
    case(14) ! file name is asked inside macbeg
!=================================================================
       call macbeg(cline,last,logok)
       if(buperr.ne.0 .or. gx%bmperr.ne.0) goto 990
!=================================================================
! ABOUT
    case(15)
!=================================================================
       write(kou,15010)version,linkdate
15010  format(/'This is OpenCalphad (OC), a free software for ',&
            'thermodynamic calculations as described in:'/&
            'B Sundman, U R Kattner, M Palumbo and S G Fries, ',&
            'Int Mat and Manu Innov (2015) 4:1; '/&
            'B Sundman, X-G Lu and H Ohtani, Comp Mat Sci, Vol 101 ',&
            '(2015) 127-137'/'B Sundman, N Dupin and B Hallstedt, ',&
            'Calphad, Vol 75 (2021) 102330'//&
            'It is available for download at http://www.opencalphad.org or'/&
            'the sundmanbo/opencalphad repository at https://www.github.com'//&
            'This software is protected by the GNU General Public License'/&
            'either version 2 of the license, or any later version.'/&
            'You may use freely and distribute copies as long as you provide ',&
            'the source code'/'and use the GNU GPL license also for your own',&
            ' additions and modifications.'//&
            'The software is provided "as is" without any warranty of any ',&
            'kind, either'/'expressed or implied.  ',&
            'The full license text is provided with the software'/&
            'or can be obtained from the Free Software Foundation ',&
            'http://www.fsf.org'//&
            'Copyright 2011-2022, Bo Sundman, Gif sur Yvette, France.'/&
            'Contact person Bo Sundman, bo.sundman@gmail.com'/&
            'This version ',a,' was compiled ',a/)
!=================================================================
! debug and subcommands
    case(16)
!=================================================================
       kom2=submenu(cbas(kom),cline,last,cdebug,ncdebug,1,'?TOPHLP')
       debug: SELECT CASE (kom2)
!------------------------------
       CASE DEFAULT
          write(kou,*)'Debug subcommand error ',kom2
!------------------------------
! debug free lists
       CASE(1)
! testing quads to cats
          call gparidx('Quad index',cline,last,i1,i2,'?MQMQA debug')
          write(*,*)'MM i1, i2=',i1,i2
          call quad2cat(i1,i2,iphl)
          write(*,*)'MM ',i1, i2, iphl(1), iphl(2)
          goto 100
       CASE(2)
          if(stop_on_error) then
             stop_on_error=.FALSE.
             write(kou,*)'No longer stop on error'
          else
             write(kou,*)'Macro will stop if error'
             stop_on_error=.true.
          endif
!------------------------------
! debug parameter structure
       case(3) ! advanced listing of data structure
! open file
! negative is for write, 0 read without filter, -100 write without filter
          write(*,*)'PMON: DEBUG case 3'
          ztyp=-7
          call gparfilex('Output file',cline,last,1,string,'  ',ztyp,&
               '?Debug parameter')
          if(string(1:1).eq.' ') then
             string='parameter_output.DAT'
             write(kou,*)' *** No file name given, will use: ',trim(string)
          endif
! add extention .dat if to extenstion provided
! close any previous output file          
          close(21)
          open(21,file=string,access='sequential',status='unknown',&
               err=990, iostat=buperr)
          lut=21
! select phase
          phlistloop: do while(name1(1:1).ne.' ')
             call gparcx('Phase: ',cline,last,1,name1,' ',&
                  '?Set Advanced')
             if(name1(1:1).eq.' ') exit phlistloop
             call find_phase_by_name(name1,iph,ics)
             if(gx%bmperr.ne.0) then
                write(*,*)'No such phase; ',trim(name1)
                exit phlistloop
             endif
             call get_phase_record(iph,lokph)
! list parameter structure
             write(*,*)'MM parameter listing will be wrong or fail listing'
             call debug_phaseparameters(lokph,lut,ceq)
! another phase
          enddo phlistloop
          write(lut,*)'Closing file'
          write(kou,*)'Closing file ',trim(string)
          close(lut)
!----------------------------------
! debug species
       case(4)
          do i1=1,nosp()
             call get_species_location(i1,loksp,name1)
             if(gx%bmperr.ne.0) goto 990
             call get_species_component_data(loksp,i2,iphl,stoik,xxx,&
                  xxy,ceq)
             if(gx%bmperr.ne.0) goto 990
             write(kou,1670)i1,loksp,name1,xxx,xxy,(iphl(j4),stoik(j4),j4=1,i2)
1670         format(2i4,1x,a12,1x,2F6.2,2x,10(i3,1x,F7.4))
          enddo
!---------------------------------
! debug tpfun
       case(5)
          call gparidx('Function index:',cline,last,ll,-1,'?Debug TPfun')
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
          call gparcdx('File name: ',cline,last,5,model,&
               './manual\html\ochelp.html#selectelement ','?Debug browser')
! this browser can be opened without ""
          browser='C:\PROGRA~1\INTERN~1\iexplore.exe '
! it works to open the ochelp.html on the first page
! gnu fortran ...
! now the complicated one ...
          string=trim(browser)//&
               '     "file://C:\Users\bosse\documents\oc\oc\src\'//&
               trim(model(3:))//'"'
          write(*,'(a)')trim(string)
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
          call gparcdx('Read TDB debug?',cline,last,1,ch1,'Y','?Debug dbcheck')
          if(ch1.eq.'Y') then
             dbcheck=.TRUE.
          else
             dbcheck=.FALSE.
          endif
          call gparcdx('HTML help?',cline,last,1,ch1,'Y','?Debug trace')
          if(ch1.eq.'Y') then
             helptrace=.TRUE.
          else
             helptrace=.FALSE.
          endif
          call gparcdx('plotting?',cline,last,1,ch1,'N','?Debug plot')
          if(ch1.eq.'Y') then
             plottrace=.TRUE.
          else
             plottrace=.FALSE.
          endif
!..................................
! debug symbol value
       case(8)
          call ui_debug_symbol_value(cline,last,ceq)
!..................................
! debug map_startpoints commented away
       case(9)
!..................................
! debug grid.  This calculates grid for phases one by one and check
       case(10)
! the first argument is set inside check_all_phases, it cannot be a constant
          i1=0
          call check_all_phases(i1,ceq)
!..................................
! DEBUG ASYMMETRIES ... 
       case(11)
          call gparcx('Phase name: ',cline,last,1,name1,' ','?Debug asymmetry')
          if(buperr.ne.0) goto 990
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
!
          lokcs=phasetuple(iph)%lokvares
          parres=>ceq%phase_varres(lokcs)
          write(*,*)'MM not implemented'
! emergency debug as I forgotten the whole data structure
!........................
       case(12) ! test bombmatta for mapping
          nullify(starteqs(1)%P1)
          starteqs(1)%P1=>ceq
          call bombmatta(maptop,noofaxis,axarr,seqxyz,starteqs)
       END SELECT debug
!=================================================================
! SELECT command
    case(17)
!=================================================================
       kom2=submenu(cbas(kom),cline,last,cselect,nselect,1,'?TOPHLP')
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
          call gparcdx('Give name or number?',cline,last,1,text,&
               name1,'?Select equilibrium')
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
             j4=1
             call getint(text,j4,i1)
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
          write(kou,'(a,i4,", name: ",a)')'Current equilibrium no: ',&
               ceq%eqno,ceq%eqname
!-----------------------------------------------------------
       CASE(2) ! select minimizer
          write(kou,*)'Sorry, only one available: ',minimizers(2)
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
          call gparcdx('Do you want to use LMDIF?',cline,last,1,ch1,'Y',&
               '?Select optimizer')
          if(ch1.eq.'Y') then
             optimizer=1
          else
             write(*,*)'Sorry VA05AD is no longer available'
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
       kom2=submenu(cbas(kom),cline,last,crej,nrej,6,'?TOPHLP')
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
          call gparcx('Phase name: ',cline,last,1,name1,' ','?Delete phase')
          if(buperr.ne.0) goto 990
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
          call remove_composition_set(iph,.FALSE.)
          if(gx%bmperr.ne.0) goto 990
!-----------------------------------------------------------
! delete equilibria
       case(6)
          call gparcx('Equilibrium name or abbr.:',cline,last,1,name1,' ',&
               '?Delete equilibrium')
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
          nullify(starteqs(1)%p1)
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
!=================================================================
!         ['NORMAL          ','SEPARATE        ','QUIT            ',&
!          'CONDITIONAL     ','TZERO           ','LIQUID_EET      ',&
!          'SHEIL_GULLIVER  ','PARAEQUILIBRIUM ','FAST            ']
    case(19)
! disable continue optimization
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
! forget any previous step special
       stepspecial=.FALSE.
! IMPORTANT I have changed the order between option and reinitiate!!
       kom2=submenu('Step options?',cline,last,cstepop,nstepop,1,'?TOPHLP')
! skip here for step quit and other cases not implemented: case(3 and 4)
       if(kom2.eq.3 .or. kom2.eq.4) goto 100
! check if there are previous results
       call ui_delete_step_map(1,cline,last,ceq)
! indicate to graphics that we have a step calculation
       graphopt%noofcalcax=1
       step: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'No such step option'
!-----------------------------------------------------------
! STEP NORMAL
       case(1)
! maptop is returned as main map/step record for results
! noofaxis is current number of axis, axarr is array with axis data
! starteq is start, equilibria, if empty set it to ceq
          if(noofstarteq.eq.0) then
             noofstarteq=1
             starteqs(1)%p1=>ceq
          endif
          if(associated(maptop)) then
! can one have several STEP commands YES!
             write(*,*)'Deleting previous step/map results missing'
             goto 100
          endif
! seqzyz are initial values for creating equilibria for lines and nodes
          call map_setup(maptop,noofaxis,axarr,seqxyz,starteqs)
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
             write(kou,'(a)')'Link set to previous step/map results.'
             maptop%plotlink=>maptopsave
          endif
! debugging: last maptop/line used
! remove start equilibria
          nullify(starteqs(1)%p1)
          noofstarteq=0
          if(gx%bmperr.ne.0) goto 990
!-----------------------------------------------------------
! STEP SEPARATE
       case(2) ! calculate for each entered phase separately (one by one)
! it will always use the current equilibrium
! can one have several STEP commands??
          stepspecial(1)=.TRUE.
          if(associated(maptop)) then
             write(*,*)'Deleting previous step/map results missing'
             goto 100
          endif
          call step_separate(maptop,noofaxis,axarr,seqxyz,ceq)
! mark that interactive listing of conditions and results may be inconsistent
          ceq%status=ibset(ceq%status,EQINCON)
          if(.not.associated(maptop)) then
! if one has errors in map_setup maptop may not be initiated, if one
! has saved previous calculations in maptopsave restore those
             if(associated(maptopsave)) then
                write(kou,*)'Restoring previous map results'
                maptop%plotlink=>maptopsave
                nullify(maptopsave)
             endif
          elseif(associated(maptopsave)) then
             write(kou,'(a)')'Link set to previous map/step results'
             maptop%plotlink=>maptopsave
          endif
! set default yaxis as GM(*)
          if(axplotdef(2)(1:1).eq.' ') then
             axplotdef(2)='GM(*)'
          endif
! update maptop%seqx to maptop%prvious%seqx+1 to allow more maptop records
          maptop%seqx=maptop%previous%seqx+1
! remove start equilibria
          nullify(starteqs(1)%p1)
          noofstarteq=0
          stepspecial(1)=.TRUE.
!-----------------------------------------------------------
! STEP QUIT, note quitting already above
       case(3)
!-----------------------------------------------------------
! STEP CONDITIONAL (NOT for Scheil-Gulliver), note quitting already above
       case(4)
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! STEP TZERO plotlink
       case(5)
          call ui_step_tzero(cline,last,tzcond,name2,ceq)
!-----------------------------------------------------------
! STEP LIQUID_EET, not yet implemented
       case(6)
          write(kou,*)'Not implemented yet'
          goto 100
       case(7,9)
          call ui_step_scheil(cline,last,kom2,ceq)
!-----------------------------------------------------------
! STEP PARAEQUILIBRIUM
       case(8)
          write(kou,874)
874       format('Before this command you must have set an alloy composition',&
               ' and calculated',/&
               'and suspended all phases except the two involved and',&
               ' you should have'/&
               'calculated a paraequilibrium')
          call gparcx('Have you done all that?',cline,last,1,&
               name1,'NO','?Step paraeq')
          call capson(name1)
          if(name1(1:1).ne.'Y') goto 100
          if(dummy(1:1).ne.' ') dummy=name2
          call gparcdx('Matrix phase ',cline,last,1,name2,dummy,'?Step paraeq')
          call find_phasetuple_by_name(name2,tupix(1))
          if(gx%bmperr.ne.0) goto 990
          if(dummy(1:1).ne.' ') dummy=name3
          call gparcdx('Growing phase ',cline,last,1,name3,dummy,'?Step paraeq')
          call find_phasetuple_by_name(name3,tupix(2))
          if(gx%bmperr.ne.0) goto 990
          dummy=name3
          call gparcdx('Fast diffusing element',cline,last,1,&
               elsym,parael,'?Step paraeq')
          call find_element_by_name(elsym,icond)
          parael=elsym
          write(kou,875)trim(name1)
875       format('The simulation will vary the axis variable and calulate',&
               ' compositions'/'of the two phases which have the same',&
               ' chemical potential of ',a)
! Delete previous step/map results
          if(associated(maptop)) then
             write(kou,*)'Previous map/step results will be deleted'
             call delete_mapresults(maptop)
          endif
          nullify(maptop)
          nullify(maptopsave)
! This is to keep trace of the total number of saved equilibria
          totalsavedceq=0
! initiate indexing nodes and lines
          seqxyz=0
! remove all graphopt settings
          call reset_plotoptions(graphopt,plotfile,textlabel)
! set default plot axis
          axplotdef(1)='W(*,'//trim(parael)//') '
! one can calculate paraequilibria diagrams at constant T
! step paraequil
          stepspecial(4)=.TRUE.
          call step_paraequil(maptop,noofaxis,axarr,seqxyz,tupix,icond,ceq)
          if(gx%bmperr.ne.0) goto 990
! sum the points calculated
          write(kou,'(a,2i5,a)')'Paraequilibrium points: ',totalsavedceq
!-----------------------------------------------------------
! STEP FAST part of STEP SCHEIL
       end SELECT step
!=================================================================
! MAP, must be tested if compatible with assessments
!=================================================================
    case(20)
       call ui_map(cline,last,ceq)
! end of MAP command
!=================================================================
! PLOT COMMAND with many options and EXTRA
!=================================================================
!        ['RENDER          ','SCALE_RANGES    ','FONT            ',&
!         'AXIS_LABELS     ','                ','TITLE           ',&
!         'GRAPHICS_FORMAT ','OUTPUT_FILE     ','                ',&
!         'QUIT            ','POSITION_OF_KEYS','APPEND          ',&
!         'TEXT_LABEL      ','                ','EXTRA           ']
! Always specify the axis first when giving this command, default is previous!!
! loop with subcommands comes after
    case(21)
       if(.not.associated(maptop)) then
          write(kou,*)'You must give a STEP or MAP command before PLOT'
          goto 100
       endif
       call ui_plot_axes(cline,last,jump,ceq)
       if(jump.eq.100) goto 100
       if(jump.eq.990) goto 990
! first argument is the number of plot axis, always 2 at present
       jp=2
       if(associated(maptopsave)) then
          write(kou,'(a)')'Link set to maptopsave'
          maptop%plotlink=>maptopsave
       endif
! restore default graphopt%linetype
!-----------------------------------------------------------
! PLOT subcommands, default is PLOT, NONE does not work ...
! subcommands to PLOT OPTIONS/ GRAPHICS OPTIONS
! THIS IS A MESS, should be reorganized in levels
!    character (len=16), dimension(nplt) :: cplot=&
!        ['RENDER          ','SCALE_RANGES    ','FONT            ',&
!         'AXIS_LABELS     ','                ','TITLE           ',&
!         'GRAPHICS_FORMAT ','OUTPUT_FILE     ','                ',&
!         'QUIT            ','POSITION_OF_KEYS','APPEND          ',&
!         'TEXT_LABEL      ','                ','EXTRA           ']
!-------------------
! return here after each sub or subsub command
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
       kom2=submenu('Plot options?',cline,last,cplot,nplt,1,'?TOPHLP')
       plotoption: SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'No such plot option'
!-----------------------------------------------------------
! PLOT RENDER no more options to plot ...
       case(1)
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
          graphopt%filename=plotfile
          if(maptop%tieline_inplane.lt.0) then
! set the isopleth bit
             if(index(graphopt%pltax(1),'*').eq.0 .and. &
                  index(graphopt%pltax(1),'*').eq.0) then
                graphopt%status=ibset(graphopt%status,GRISOPLETH)
             else
! Probably meaningless to plot fractions ... but who knows?
                graphopt%status=ibclr(graphopt%status,GRISOPLETH)
             endif
          else
! for step and tie-lines in plane clear the bit
             graphopt%status=ibclr(graphopt%status,GRISOPLETH)
          endif
! added ceq in the call to make it possible to handle change of reference states
          call ocplot2(jp,maptop,axarr,graphopt,version,ceq)
          if(gx%bmperr.ne.0) goto 990
! always restore default plot file name and plot option to screem
          if(graphopt%gnutermsel.ne.1) &
               write(kou,*)'Restoring plot device to screen'
          graphopt%gnutermsel=1
          plotfile='ocgnu'
          graphopt%filename=plotfile
!-----------------------------------------------------------
! PLOT SCALE_RANGE of either X or Y
       case(2)
          call ui_plot_scale_ranges(cline,last,jump)
          if(jump.eq.100) goto 100
          goto 21100
!-----------------------------------------------------------
! PLOT unused select FONT
       case(3)
          call gparcdx('Font (check what your GNUPLOT has): ',&
               cline,last,1,name1,graphopt%font,'?Plot font')
          graphopt%font=name1
! font size ignored but it is better to have the question now ...
          call gparidx('Font size: ',cline,last,iz,16,'?Plot font')
          write(*,*)'Size is ignored at present ...'
! we have to change "font" in all terminals and key
          allgnu: do i1=1,graphopt%gnutermax
             iz=index(graphopt%gnuterminal(i1),'"')
             if(iz.le.0) cycle allgnu
             i2=index(graphopt%gnuterminal(i1)(iz:),',')
             name1=graphopt%gnuterminal(i1)(iz+i2-1:)
             graphopt%gnuterminal(i1)(iz+1:)=graphopt%font
             i2=len_trim(graphopt%gnuterminal(i1))
             graphopt%gnuterminal(i1)(i2+1:)=name1
          enddo allgnu
          iz=index(graphopt%labelkey,'"')
          if(iz.gt.0) then
             i2=index(graphopt%labelkey(iz:),',')
             name1=graphopt%labelkey(iz+i2-1:)
             graphopt%labelkey(iz+1:)=graphopt%font
             i2=len_trim(graphopt%labelkey)
             graphopt%labelkey(i2+1:)=name1
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT AXIS_LABELS
       case(4)
          call gparcdx('For X or Y axis? ',cline,last,1,ch1,'X',&
               '?Plot axis labels')
          if(ch1.eq.'X' .or. ch1.eq.'x') then
             call gparcdx('Axis label: ',cline,last,5,&
                  graphopt%plotlabels(2),axplot(1),'?Plot axis labels')
! remember that plotlabel(1) is the title
             graphopt%labeldefaults(2)=len(graphopt%plotlabels(2))
          elseif(ch1.eq.'Y' .or. ch1.eq.'y') then
             call gparcdx('Axis label: ',cline,last,5,&
                  graphopt%plotlabels(3),axplot(2),'?Plot axis labels')
! remember that plotlabel(1) is the title
             graphopt%labeldefaults(3)=len(graphopt%plotlabels(3))
          else
             write(kou,*)'Please answer X or Y'
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT unused
       case(5)
!-----------------------------------------------------------
! PLOT TITLE
       case(6)
          call gparcdx('Plot title',cline,last,5,line,'DEFAULT','?Plot title')
          if(line(1:8).eq.'DEFAULT ') then
             graphopt%labeldefaults(1)=0
          else
             graphopt%plotlabels(1)=line
             graphopt%labeldefaults(1)=len_trim(graphopt%plotlabels(1))
          endif
          goto 21100
!-----------------------------------------------------------
! PLOT GRAPHICS_FORMAT and PLOT OUTPUT_FILE
! when setting graphics format always also ask for plot file
       case(7,8)
          call ui_plot_format_file(cline,last,kom2)
          goto 21100
!-----------------------------------------------------------
! PLOT unused
       case(9)
!-----------------------------------------------------------
! PLOT QUIT
       case(10)
! just return to command level
!-----------------------------------------------------------
! PLOT position of line labels (position_of_keys)
       case(11)
          write(kou,21200)
21200     format('Key to lines can be positioned: '/&
               'top/bottom left/center/right inside/outside on/off')
          call gparcdx('Position?',cline,last,5,line,'top right','?Plot keys')
          graphopt%labelkey=trim(line)
          if(line(1:3).ne.'off') then
             call gparidx('Size: ',cline,last,iz,12,'?Plot keys')
             graphopt%labelkey=trim(graphopt%labelkey)//' font "'&
                  //trim(graphopt%font)//','
             ll=len_trim(graphopt%labelkey)+1
             write(graphopt%labelkey(ll:),'(i2,a)')iz,'"'
          endif
          write(*,*)'GNUPLOT will use: set key ',trim(graphopt%labelkey)
          goto 21100
!-----------------------------------------------------------
! PLOT APPEND a gnuplot file or csv file
       case(12)
          call ui_plot_append(cline,last,string)
          goto 21100
!-----------------------------------------------------------
! PLOT TEXT anywhere on plot
       case(13)
          call ui_plot_text(cline,last,ceq)
          goto 21100
!---------------------------------------------------------
! PLOT unused
       case(14)
!---------------------------------------------------------
! PLOT EXTRA, subsubcommand
       case(15)
! subsubcommands to PLOT (may not be updated!!)
!    character (len=16), dimension(nplt2) :: cplot2=&
!        ['COLOR           ','LOGSCALE        ','RATIOS_XY       ',&
!         'LINE_TYPE       ','MANIPULATE_LINES','PAUSE_OPTION    ',&
!         'LOWER_LEFT_TEXT ','TIE_LINES       ','GIBBS_TRIANGLE  ',&
!         'QUIT            ','SPAWN           ','NO_HEADING      ',&
!         'AXIS_FACTOR     ','GRID            ','                ',&
!         '                ','                ','                ']
!-------------------------------------------------------------------
! default set to GIBBS-TRIANGLE
          kom3=submenu('Extra options?',cline,last,cplot2,nplt2,9,'?TOPHLP')
          plotextra: SELECT CASE(kom3)
          case default
! this is typically when using a ? or ??
             write(*,*)'No such extra option'
             cline=' '
             last=len(cline)
             goto 21100
!...............................................
! PLOT EXTRA axis_factor for example to plot kJ or GPa instead of J and Pa
          case(13)
             call gparcdx('Which axis?',cline,last,1,ch1,'Y',&
                  '?Plot extra factor')
             call capson(ch1)
             if(ch1.eq.'Y' .or. ch1.eq.'X') then
                call gparrdx('Factor?',cline,last,xxx,1.0D-3,&
                     '?Plot extra factor')
                if(ch1.eq.'X') graphopt%scalefact(1)=abs(xxx)
                if(ch1.eq.'Y') graphopt%scalefact(2)=abs(xxx)
             else
                write(*,*)'No such axis'
             endif
             goto 21100
!...............................
! PLOT EXTRA COLOR ... and some more things ...
          case(1)
! monovariant and tielinecolor declared in smp2.F90
             call gparcdx('Monovariant color ',cline,last,1,&
                  name1,monovariant,'?Plot color')
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
             call gparcdx('Tie-line color ',cline,last,1,&
                  name1,tielinecolor,'?Plot font')
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
!...............................................
! PLOT EXTRA Gibbs triangle
          case(9)
             chz='Y'
             if(graphopt%gibbstriangle) chz='N'
             call gparcdx('A Gibbs triangle diagram?',cline,last,5,ch1,chz,&
                  'PLOT Gibbs triangle')
             if(ch1.eq.'y' .or. ch1.eq.'Y') then
                graphopt%gibbstriangle=.TRUE.
                write(*,22500)
22500           format('The Gibbs triangle layout courtesy of',&
                     ' Catalina Pineda Heresi at RUB, Germany')
             else
                graphopt%gibbstriangle=.FALSE.
             endif
             goto 21100
!...............................................
! PLOT EXTRA GRID
          case(14)
             call gparcdx('Plot grid?',cline,last,1,ch1,'Y',&
                  '?Plot extra factor')
             call capson(ch1)
             if(ch1.eq.'Y') then
                graphopt%setgrid=1
             else
                graphopt%setgrid=0
             endif
             goto 21100
!...............................................
! PLOT EXTRA LINE_TYPE
          case(4)
             j4=last
             if(eolch(cline,j4)) then
! write this only if the lime is empty
                write(*,22300)
22300           format('Default 1 restore normal line types:',/&
                     ' 0 means dashed lines,'/,' 1 means full line',/&
                     '>1 means symbol at each calculated point')
             endif
             call gparidx('Line type?',cline,last,iz,1,'?Plot line symbols')
             if(iz.eq.0) then
! this means dashed lines and possibly symbols if already set ..
                graphopt%linetype=0
             elseif(iz.gt.1) then
! this means symboles and possibly dashed lines if already set
                graphopt%linewp=iz
             else
! this means full lines and no symbols
                graphopt%linewp=1
                graphopt%linetype=1
             endif
             goto 21100
!...............................................
! PLOT EXTRA LOGSCALE
          case(2)
             call gparcdx('For x or y axis (or NONE)? ',cline,last,1,ch1,'y',&
                  '?Plot logax')
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
                write(kou,*)'Both axis set to be linear'
                graphopt%axistype(1)=0
                graphopt%axistype(2)=0
             endif
             goto 21100
!...............................................
! PLOT EXTRA text in lower left corner
          case(7)
             call gparcx('Text in lower left corner?',cline,last,1,text,' ',&
                  '?Extra lower-left-corner')
             graphopt%lowerleftcorner=text
             goto 21100
!...............................................
! PLOT EXTRA MANIPULATE LINE COLORS
          case(5)
             write(kou,22400)
22400        format('OC uses GNUPLOT and it is possible to edit',&
                  ' the file "ocgnu.plt" file'/&
                  'generated by OC to use extensive facilities',&
                  ' provided by GNUPLOT.'/&
                  'Only a few of them is provided here.'/&
                  'OC has 10 different colors to identify the lines plotted.',&
                  ' Line 11 or'/' higher will repeat these colors.  With',&
                  ' this command you can select'/' one of these 10 colors',&
                  ' to be used for the first line plotted.')
             call gparidx('The color index should be on the first line?',&
                  cline,last,flc,1,'?Plot manipulate colors')
             if(flc.lt.1 .or. flc.gt.10) then
                write(*,*)'Number must be between 1 and 10'
             else
                graphopt%linett=flc
             endif
             goto 21100
!...............................................
! PLOT EXTRA remove headings
          case(12)
             call gparcdx('Remove headings?',cline,last,1,ch1,'N',&
                  '?Plot no heading')
             if(ch1.ne.'N') then
                write(*,*)'No title set!',ch1
                graphopt%status=ibset(graphopt%status,GRNOTITLE)
             else
                graphopt%status=ibclr(graphopt%status,GRNOTITLE)
             endif
             goto 21100
!...............................................
! PLOT EXTRA PAUSE_OPTIONS uselss??
          case(6)
             write(kou,*)'Specify option after pause !'
             call gparcx('GNUPLOT pause option?',cline,last,5,text,' ',&
                  '?Plot pause')
             if(len_trim(text).eq.0) then
                write(kou,*)'Warning, plot will exit directly!'
             endif
             graphopt%plotend='pause '//text
             goto 21100
!...............................................
! PLOT EXTRA QUIT
          case(10)
             goto 21100
!...............................................
! PLOT EXTRA RATIOS
          case(3)
             call gparrdx('X-axis plot ratio',cline,last,xxx,graphopt%xsize,&
                  '?Plot ratios')
             if(xxx.le.0.1) then
                write(*,*)'Ratio set to 0.1'
                xxx=0.1D0
             endif
             graphopt%xsize=xxx
             call gparrdx('Y-axis plot ratio',cline,last,xxx,graphopt%ysize,&
                  'PLOT ratios')
             if(xxx.le.0.1) then
                write(*,*)'Ratio set to 0.1'
                xxx=0.1D0
             endif
             graphopt%ysize=xxx
             goto 21100
!...............................................
! PLOT EXTRA spawn plot
          case(11)
             call gparcdx('Spawn plot?',cline,last,1,ch1,'N','?Plot extra')
             if(ch1.eq.'Y') then
                graphopt%status=ibset(graphopt%status,GRKEEP)
             else
                graphopt%status=ibclr(graphopt%status,GRKEEP)
             endif
             goto 21100
!...............................................
! PLOT EXTRA Tie-line plot increment
          case(8)
             call gparidx('Tie-line plot increment?',cline,last,kl,3,&
                  '?Plot tieline')
             if(kl.lt.0) kl=0
             graphopt%tielines=kl
             goto 21100
!...............................................
! PLOT EXTRA unused
          case(15)
             goto 21100
!...............................................
! PLOT EXTRA unused
          case(16)
             goto 21100
!...............................................
! PLOT EXTRA unused
          case(17)
             goto 21100
!...............................................
! PLOT EXTRA unused
          case(18)
             goto 21100
!-----------------------------------------------------------
          end select plotextra
          goto 21100
!-----------------------------------------------------------
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
!=================================================================
    case(24)
       call ui_optimize(cline,last,lut)
! end of call to LMDIF
!=================================================================
! SHOW is immpemented as a special case of LIST STATE_VARIABLES
!    case(25) 
!=================================================================
! not used
    CASE(26,27,28,29,30)
       write(kou,*)'Not implemented yet'
!=================================================================
!
    END SELECT main
! command executed, prompt for another command unless error code
    if(gx%bmperr.eq.0) goto 100
!============================================================
! handling errors
990 continue
    write(kou,991)gx%bmperr,buperr,kiu
991 format(/'Error codes: ',3i6)
    if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
       write(kou,992)trim(bmperrmess(gx%bmperr))
992    format('Message: ',a/)
    else
       write(kou,*)'No defined error message, maybe I/O error'
    endif
    if(stop_on_error) then
! error running a macro, terminate macro and return interactive
       write(*,*)'Error running MACRO file, return to interactive mode?'
       if(iox(8).eq.0) then
! iox(8) is nonzero if one has set "no stop on @&"
! in such a case ignore any error
          read(*,993)ch1
993       format(a)
          if(ch1.eq.' ') then
             call macend(cline,last,logok)  
             kiu=kiud
          else
             write(*,*)'Continue on your own risk '
          endif
       endif
    endif
! turn off macro but remain inside software
    gx%bmperr=0; buperr=0
    goto 100
!
  end subroutine oc_command_monitor

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  integer function submenu(query,cline,last,ccomnd,ncomnd,kdef,menutarget)
! general subcommand decoder
! query is the prompt
! cline and last is user input and position
! ccomnd is the menu and ncomnd number of menu entries
! kdef is the default (to be added to query)
!  implicit double precision (a-h,o-z)
    implicit none
    character cline*(*),ccomnd(*)*(*),query*(*),menutarget*(*)
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
    if(last+2.lt.len(cline)) then
       if(cline(last:last).eq.',') last=last+1
    else
       last=len(cline)-2
    endif
    if(cline(last:last+2).eq.' ? ' .or. cline(last:last+1).eq.'? ') then
! This handles help for things like "set ?"
! if cline is just "?" just display menu but reset last to 1
! in order to provide help also specific for the command
       call q3helpx(cline,last,ccomnd,ncomnd)
       last=len(cline)
       goto 1000
    endif
100 continue
    if(kdef.lt.1 .or. kdef.gt.ncomnd) then
! no default answer
       if(eolch(cline,last)) then
! empty line, note fourth argument 5 copes whole of cline into text
! the hypertext is the submenu prompt, last argument set to mark TOPHLP!
          call gparcx(query1(1:lenq),cline,last,5,text,' ','?TOPHLP')
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
333       continue
! this is submenu command input
! note fourth argument 5 copes whole of cline into text
! the queary is the hypertarget, last argument to indicate TOPHLP
! if user answers a single ? retutn here, if ?? use webrowser
          call gparcdx(query1(1:lenq),cline,last,5,text,defansw,'?TOPHLP')
          if(buperr.ne.0) goto 1000
          if(text(1:1).eq.'?') then
             if(text(2:2).ne.'?') then
! if we have just a ? here we should display the menue
                call q3helpx(' * ',last,ccomnd,ncomnd)
                last=len(cline)
! note that two ?? should have been taken care of inside gparcdx
             endif
             goto 333
          endif
          cline=text
       else
! if first character is , take default answer
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
! in this case there is no user input in this call
             call gparcdx(query1(1:lenq),cline,last,5,text,defansw,'?TOPHLP')
             if(buperr.ne.0) goto 1000
             cline=text
! added 20190207 because "enter gamma ac(a)/x(a); gave segmentation fault
! but that was not the error, the error was missing =
          endif
       endif
    endif
!
    kom2=ncomp(cline,ccomnd,ncomnd,last)
    if(kom2.le.0) then
       if(once) then
          if(cline(1:1).ne.'?') once=.false.
          if(kom2.lt.0) write(kou,*)'Ambiguous answer, please try again'
          write(kou,*)'Possible answers are:'
          last=1
          cline=' *'
          call q3helpx(cline,last,ccomnd,ncomnd)
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
    integer next,kom,slen,errno,jj,ztyp
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
       call q3helpx(dummy,next,copt,nopt)
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
          call q3helpx(dummy,next,copt,nopt)
          afo=1
!-----------------------------------
       case(1) ! /output means open a file and ovewrite any previous content
! next argument after = must be a file name
! 6 means extension DAT
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
          ztyp=-7
          call gparfilex('Output file',option,next,1,string,'  ',ztyp,&
               '?Command options')
          if(string(1:1).eq.' ') then
             string='ocoutput.DAT'
             write(kou,*)' *** No file name given, will use: ',trim(string)
          endif
          slen=len_trim(string)
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
! next argument after = must be a file name
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
          ztyp=-7
          call gparfilex('Append to file:',option,next,&
               1,string,'  ',ztyp,'?Command options')
          if(string(1:1).eq.' ') then
             string='ocappend.DAT'
             write(kou,*)' *** No file name given, will use: ',trim(string)
          endif
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
       case(7,8,9) ! 
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
    return
  end subroutine ocmon_reset_options

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

! subroutines extracted from oc_command_monitor
#include "pmon7A.F90"

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

END MODULE cmon1oc
