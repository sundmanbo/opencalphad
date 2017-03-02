!
MODULE cmon1oc
!
! Copyright 2012-2015, Bo Sundman, France
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
  subroutine oc_command_monitor(version,linkdate)
! command monitor
    implicit none
!
! passed as date when program was linked
    character linkdate*(*),version*(*)
! various symbols and texts
    character :: ocprompt*4='OC4:'
    character name1*24,name2*24,line*80,model*72,chshort
    integer, parameter :: ocmonversion=30
! element symbol and array of element symbols for database use
    character elsym*2,ellist(maxel)*2
! more texts for various purposes
    character text*72,string*256,ch1*1,selection*27,funstring*1024
    character axplot(3)*24,axplotdef(3)*24,quest*20
    character plotform*32,longstring*2048,optres*40
! separate file names for remembering and providing a default
    character ocmfile*64,ocufile*64,tdbfile*64,ocdfile*64,filename*64
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
! plot ranges, texts and defaults
    type(graphics_options) :: graphopt
! plot texts
!    type(graphics_textlabel), allocatable, target :: textlabel
    type(graphics_textlabel), pointer :: textlabel
    type(graphics_textlabel), pointer :: labelp
! axis data structures
    type(map_axis), dimension(5) :: axarr
! if more than one start equilibrium these are linked using the ceq%next index
    type(gtp_equilibrium_data), pointer :: starteq
! for map results
    type(map_node), pointer :: maptop,mapnode,maptopsave
!    type(map_line) :: mapline
! seqxyz has initial values of seqx, seqy and seqz
    integer noofaxis,noofstarteq,seqxyz(3)
! this should be removed
!    TYPE(ssm_node), pointer :: resultlist
!<<<<<<<--------------------------------------------------------------
! used for element data and R*T
    double precision h298,s298,rgast
! temporary reals
    double precision xxx,xxy,totam
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
    integer i1,i2,j1,iax
! more temporary integers
    integer jp,kl,svss,language,last,leak
! and more temporary integers
    integer ll,lokcs,lokph,lokres,loksp,lrot,maxax
! and more temporary integers
    integer mode,ndl,neqdef,noelx,nofc,nopl,nops,nv,nystat
! temporary matrix
!    double precision latpos(3,3)
! used to call init_gtp for the NEW command
    integer intv(10)
    double precision dblv(10)
!-------------------
! variables for assessment using VA05AD
    integer :: nopt=100,iprint=1,ient,mexp=0,nvcoeff,nwc
    integer iexit(5)
    double precision :: dstep=1.0D-4,dmax2=1.0D2,acc=1.0D-3
    integer, parameter :: maxw=5000
! occational segmentation fault when deallocating www ....
    double precision, dimension(maxw) :: www
!    double precision, dimension(:), allocatable :: www
    double precision, dimension(:), allocatable :: coefs
    double precision, dimension(:), allocatable :: errs
!    external new_assessment_calfun
!-------------------
! loop variable when entering constituents of a phase
    integer icon
! array with constituents in sublattices when entering a phase
!    character, dimension(maxconst) :: const*24
! for macro and logfile and repeating questions
    logical logok,stop_on_error,once,wildcard,twice
! unit for logfile input, 0 means no logfile
    integer logfil
! remember default for calculate phase
    integer defcp
! for state variables as conditions
    integer istv
    double precision coeffs(10)
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
    character cline*128,option*80,aline*128,plotfile*64,eqname*24
! variable phase tuple
    type(gtp_phasetuple), pointer :: phtup
!----------------------------------------------------------------
! here are all commands and subcommands
!    character (len=64), dimension(6) :: oplist
    integer, parameter :: ncbas=30,nclist=21,ncalc=9,ncent=21,ncread=6
    integer, parameter :: ncam1=18,ncset=24,ncadv=6,ncstat=6,ncdebug=6
    integer, parameter :: nselect=6,nlform=6,noptopt=9,nsetbit=6
    integer, parameter :: ncamph=12,nclph=6,nccph=6,nrej=9,nsetph=6
    integer, parameter :: nsetphbits=15,ncsave=6,nplt=18,nstepop=6
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
        '                ','                ','                ',&
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
! NOTE a command line can contain options preceeded by /
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
         '                ','                ','                ']
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
         '                ','                ','                ']
!-------------------
! subcommands to CALCULATE
    character (len=16), dimension(ncalc) :: ccalc=&
         ['TPFUN_SYMBOLS   ','PHASE           ','NO_GLOBAL       ',&
         'TRANSITION      ','QUIT            ','GLOBAL_GRIDMIN  ',&
         'SYMBOL          ','EQUILIBRIUM     ','ALL_EQUILIBRIA  ']
!-------------------
! subcommands to CALCULATE PHASE
    character (len=16), dimension(nccph) :: ccph=&
         ['ONLY_G          ','G_AND_DGDY      ','ALL_DERIVATIVES ',&
          'CONSTITUTION_ADJ','DIFFUSION_COEFF ','                ']
!-------------------
! subcommands to ENTER
    character (len=16), dimension(ncent) :: center=&
         ['TPFUN_SYMBOL    ','ELEMENT         ','SPECIES         ',&
         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
         'CONSTITUTION    ','EXPERIMENT      ','QUIT            ',&
         'EQUILIBRIUM     ','SYMBOL          ','OPTIMIZE_COEFF  ',&
         'COPY_OF_EQUILIB ','COMMENT         ','MANY_EQUILIBRIA ',&
         'MATERIAL        ','PLOT_DATA       ','                ',&
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
        ['TDB             ','SOLGASMIX       ','UNFORMATTED     ',&
         'DIRECT          ','                ','QUIT            ']
!-------------------
! subcommands to AMEND first level
! many of these should be subcommands to PHASE
    character (len=16), dimension(ncam1) :: cam1=&
         ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
         'TPFUN_SYMBOL    ','CONSTITUTION    ','QUIT            ',&
         'COMPONENTS      ','GENERAL         ','CP_MODEL        ',&
         'ALL_OPTIM_COEFF ','EQUILIBRIUM     ','                ',&
         'LINE            ','                ','                ']
!-------------------
! subsubcommands to AMEND PHASE
    character (len=16), dimension(ncamph) :: camph=&
!         ['MAGNETIC_CONTRIB','COMPOSITION_SET ','DISORDERED_FRACS',&
         ['MAGNETIC_INDEN  ','COMPOSITION_SET ','DISORDERED_FRACS',&
         'TWOSTATE_MODEL_1','QUIT            ','DEFAULT_CONSTIT ',&
         'DEBYE_CP_MODEL  ','EINSTEIN_CP_MDL ','XIONG_INDEN_MAGN',&
         'ELASTIC_MODEL_1 ','GADDITION       ','                ']
!-------------------
! subcommands to SET
    character (len=16), dimension(ncset) :: cset=&
         ['CONDITION       ','STATUS          ','ADVANCED        ',&
         'LEVEL           ','INTERACTIVE     ','REFERENCE_STATE ',&
         'QUIT            ','ECHO            ','PHASE           ',&
         'UNITS           ','LOG_FILE        ','WEIGHT          ',&
         'NUMERIC_OPTIONS ','AXIS            ','INPUT_AMOUNTS   ',&
         'VERBOSE         ','AS_START_EQUILIB','BIT             ',&
         'VARIABLE_COEFF  ','SCALED_COEFF    ','OPTIMIZING_COND ',&
         'RANGE_EXPER_EQU ','FIXED_COEFF     ','                ']
! subsubcommands to SET STATUS
    character (len=16), dimension(ncstat) :: cstatus=&
         ['ELEMENT         ','SPECIES         ','PHASE           ',&
         'CONSTITUENT     ','                ','                ']
!        123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET ADVANCED
    character (len=16), dimension(ncadv) :: cadv=&
         ['EQUILIB_TRANSF  ','QUIT            ','EXTRA_PROPERTY  ',&
          'DENSE_GRID_ONOFF','                ','                ']
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
    character (len=16), dimension(nsetphbits) :: csetphbits=&
         ['FCC_PERMUTATIONS','BCC_PERMUTATIONS','IONIC_LIQUID_MDL',&
         'AQUEOUS_MODEL   ','QUASICHEMICAL   ','FCC_CVM_TETRADRN',&
         'FACT_QUASICHEMCL','NO_AUTO_COMP_SET','QUIT            ',&
         'EXTRA_DENSE_GRID','FLORY_HUGGINS   ','                ',&
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
          'TEST1           ','TEST2           ','                ']
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
! subcommands to PLOT OPTIONS/ GRAPHICS OPTIONS
    character (len=16), dimension(nplt) :: cplot=&
         ['RENDER          ','XRANGE          ','YRANGE          ',&
         'XTEXT           ','YTEXT           ','TITLE           ',&
         'GRAPHICS_FORMAT ','OUTPUT_FILE     ','GIBBS_TRIANGLE  ',&
         'QUIT            ','POSITION_OF_KEYS','APPEND          ',&
         'TEXT            ','TIE_LINES       ','KEEP            ',&
         'LOGSCALE        ','                ','                ']
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
! some defaults
    language=1
    logfil=0
    defcp=1
    seqxyz=0
! defaults for optimizer, iexit(2)=1 means listing scaled coefficients (Va05AD)
    nvcoeff=0
    iexit=0
    iexit(2)=1
!
    write(kou,10)version,trim(linkdate),ocmonversion,gtpversion,hmsversion,&
         smpversion
10  format(/'Open Calphad (OC) software version ',a,', linked ',a,/&
         'with command line monitor version ',i2//&
         'This program is available with a GNU General Public License.'/&
         'It includes the General Thermodynamic Package, version ',A,','/&
         "Hillert's equilibrium calculation algorithm version ",A,','/&
         'step/map/plot software version ',A,', ',&
         'LAPACK and BLAS numerical routines'/'and LMDIF from ANL ',&
         '(Argonne, USA) is used by the assessment procedure'/)
!aa
!$    write(kou,11)
11  format('Linked with OpenMp for parallel execution')
!
! jump here after NEW to reinitiallize all local variables also
20  continue
! clear file names
    ocmfile=' '; ocufile=' '; tdbfile=' '
! plot ranges and their defaults
    graphopt%gibbstriangle=.FALSE.
    graphopt%rangedefaults=0
    graphopt%labeldefaults=0
    graphopt%plotmin=zero
    graphopt%dfltmin=zero
    graphopt%plotmax=one
    graphopt%dfltmax=one
    graphopt%appendfile=' '
    graphopt%status=0
    graphopt%labelkey='top right'
    nullify(graphopt%firsttextlabel)
    nullify(textlabel)
    plotfile='ocgnu'
    plotform=' '
! default list unit
    optionsset%lut=kou
! default for list short
    chshort='A'
! set default minimizer, 2 is matsmin, 1 does not work ...
    minimizer=2
! set default optimimzer, 1 is LMDIF, 2 is VA05AD
    optimizer=1
! by default no stop on error and no logfile
    stop_on_error=.false.
    logfil=0
!
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
! state variable for plot axis
    do j1=1,3
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
! local environment
    ochome=' '
    call get_environment_variable('OCHOME ',ochome)
    if(ochome(1:1).eq.' ') then
       inquire(file='ochelp.hlp ',exist=logok)
       if(.not.logok) then
          write(*,*)'Warning, no help file'
       else
! help file on local directory?
          call init_help('ochelp.hlp ')
       endif
    else
! both LINUX and WINDOWS accept / as separator between directory and file names
!       write(*,*)'Help file: ',trim(ochome)//'\ochelp.hlp '
       call init_help(trim(ochome)//'/ochelp.hlp ')
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
!          write(*,*)'there is an initiation file!',trim(cline)
          last=0
          call macbeg(cline,last,logok)
!       else
!          write(*,*)'No initiation file'
       endif
    endif
!
! finished initiallization
88  continue
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
! read the command line with gparc to have output on logfile
    last=len(aline)
    aline=' '
    cline=' '
    call gparc(ocprompt,aline,last,5,cline,' ',tophlp)
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
       goto 100
    else
! check for options .... this does not work yet
! one should check for options after each subcommand or value entered ??
!       call ocmon_set_options(cline,last,optionsset)
       nops=0
110    continue
       if(.not.eolch(cline,last)) then
          if(cline(last:last).eq.'/') then
             call getext(cline,last,2,option,' ',nopl)
             if(buperr.ne.0) then
                write(kou,*)'Error reading option'
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
       write(*,*)'Warning, exceeded helprec%level limit'
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
    SELECT CASE(kom)
! command selection
!=================================================================
    CASE DEFAULT
       write(kou,*)'No such command'
       goto 100
!=================================================================
! amend subcommands
!       ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
!        'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
!        'TPFUN_SYMBOL    ','CONSTITUTION    ','QUIT            ',&
!        'COMPONENTS      ','GENERAL         ','                ',&
!         'ALL_OPTIM_COEFF','EQUILIBRIUM     ','                ',&
!         'LIST           ','                ','                ']
    CASE(1) ! AMEND
! disable continue optimization
       iexit=0
       iexit(2)=1
       kom2=submenu(cbas(kom),cline,last,cam1,ncam1,4)
       SELECT CASE(kom2)
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
             call gparrd('Give new value; ',cline,last,xxy,xxx,q1help)
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
! BEWARE: if equilibria are calculated in threads this must be calculated first
! even in parallel threads
          call gparcd(&
               'Should the symbol be evaluated for this equilibrium only?',&
               cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
             svflista(svss)%status=ibset(svflista(svss)%status,SVFEXT)
             svflista(svss)%eqnoval=neqdef
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
          write(kou,*)'Not implemented yet'
!-------------------------
! subsubcommands to AMEND PHASE
!       ['MAGNETIC_CONTRIB','COMPOSITION_SET ','DISORDERED_FRACS',&
!        'TWOSTATE_MODEL_1   ','QUIT            ','DEFAULT_CONSTIT ',&
!        'DEBYE_CP_MODEL  ','EINSTEIN_CP_MDL ','XIONG_INDEN_MAGM',&
!        'ELASTIC_MODEL_1 ','GADDITION       ','                ']
       case(4) ! amend phase subcommands
          call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
          if(buperr.ne.0) goto 990
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
!
          kom3=submenu(cbas(kom),cline,last,camph,ncamph,2)
          SELECT CASE(kom3)
!....................................................
          CASE DEFAULT
             write(kou,*)'Amend phase subcommand error'
!....................................................
          case(1) ! amend phase <name> magnetic contribution
             idef=-3
             call gparid('Antiferromagnetic factor: ',&
                  cline,last,j1,idef,q1help)
             if(buperr.ne.0) goto 990
             call get_phase_record(iph,lokph)
             if(j1.eq.-1) then
! Inden magnetic for BCC
                call add_addrecord(lokph,'Y',indenmagnetic)
             else
! Inden magnetic for FCC
                call add_addrecord(lokph,'N',indenmagnetic)
             endif
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
             idef=1
             call gpari('Sum up to sublattice: ',cline,last,ndl,idef,q1help)
             if(buperr.ne.0) goto 990
! ch1 is parameter suffix, j1=0 means never completely disordred (sigma)
             call gparcd('Can the phase be totally disordered? ',cline,last,&
                  1,ch1,'N',q1help)
             if(buperr.ne.0) goto 990
             if(ch1.eq.'N') then
! like sigma
                j1=0
             else
! like FCC ordering
                j1=1
             endif
             ch1='D'
             call add_fraction_set(iph,ch1,ndl,j1)
             if(gx%bmperr.ne.0) goto 990
!....................................................
          case(4) ! amend phase <name> twostate model 1
             call add_addrecord(iph,' ',twostatemodel1)
!....................................................
          case(5) ! amend phase quit
             goto 100
!....................................................
          case(6) ! amend phase <name> default_constitution
! to change default constitution of any composition set give #comp.set.
             call ask_default_constitution(cline,last,iph,ics,ceq)
!....................................................
          case(7) ! amend phase <name> Debye model
             call add_addrecord(iph,' ',debyecp)
!....................................................
          case(8) ! amend phase einstein cp model
             call add_addrecord(iph,' ',einsteincp)
!....................................................
          case(9) ! amend phase Inden-Xiong magnetic_model
             call gparcd('BCC type phase: ',cline,last,1,ch1,'N',q1help)
             call add_addrecord(iph,ch1,xiongmagnetic)
!....................................................
          case(10) ! amend phase elastic model
             call add_addrecord(iph,' ',elasticmodel1)
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
             call gparrd('Addition to G in J/FU: ',cline,last,xxx,xxy,q1help)
             ceq%phase_varres(lokcs)%addg(1)=xxx
! set bit that this should be calculated
             ceq%phase_varres(lokcs)%status2=&
                  ibset(ceq%phase_varres(lokcs)%status2,CSADDG)
!....................................................
          case(12) ! amend phase ??
             write(kou,*)'Not implemented yet'
          END SELECT
!-------------------------
       case(5) ! amend parameter
          write(kou,*)'Not implemented yet, only ENTER PARAMETER'
!-------------------------
       case(6) ! amend bibliography
          call enter_reference_interactivly(cline,last,1,j1)
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
             write(*,*)'You cannot change optimizing coefficients'
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
       case(12) ! amend Cp_model (Heat capacity)
          write(*,*)'Not implemented yet'
!-------------------------
       case(13) ! amend ALL_OPTIM_COEFF, (rescale or recover)
          call gparcd('Should the coefficients be rescaled?',&
               cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
! set start values to current values
             firstash%coeffstart=firstash%coeffvalues*firstash%coeffscale
             firstash%coeffscale=firstash%coeffstart
             firstash%coeffvalues=one
          else
             call gparcd('Do you want to recover the coefficients values?',&
                  cline,last,1,ch1,'N',q1help)
             if(ch1.eq.'y' .or. ch1.eq.'Y') then
! set current values back to start values
                firstash%coeffvalues=one
                firstash%coeffstart=firstash%coeffscale
             else
                write(kou,557)
557             format('Nothing done as there are no other amend option',&
                     'for coefficients')
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
       END SELECT
!=================================================================
! calculate subcommands
!         ['TPFUN_SYMBOLS   ','PHASE           ','NO_GLOBAL       ',&
!         'TRANSITION      ','QUIT            ','GLOBAL_GRIDMIN  ',&
!         'SYMBOL          ','EQUILIBRIUM     ','ALL_EQUILIBRIA  ']
    CASE(2)
       kom2=submenu(cbas(kom),cline,last,ccalc,ncalc,8)
       SELECT CASE(kom2)
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
2011         format(/'Calculating ',i4,' functions for T,P=',F10.2,1PE15.7/&
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
          SELECT CASE(kom3)
!.......................................................
          CASE DEFAULT
             write(kou,*)'Calculate phase subcommand error'
!.......................................................
          case(1) ! calculate case < > only G
             call calcg(iph,ics,0,lokres,ceq)
             if(gx%bmperr.ne.0) goto 990
             parres=>ceq%phase_varres(lokres)
             write(lut,2031)(rgast*parres%gval(j1,1),j1=1,4)
! G=H-T*S; H=G+T*S; S=-G.T; H = G + T*(-G.T) = G - T*G.T
             write(lut,2032)parres%gval(1,1)/parres%abnorm(1),&
                  (parres%gval(1,1)-ceq%tpval(1)*parres%gval(2,1))*rgast,&
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
             write(lut,2031)(rgast*parres%gval(j1,1),j1=1,4)
             write(lut,2041)(rgast*parres%dgval(1,j1,1),j1=1,nofc)
2041         format('dG/dy:   ',4(1PE16.8),(/9x,4e16.8))
!.......................................................
          case(3) ! calculate phase < > all derivatives
             call tabder(iph,ics,ceq)
             write(*,*)' NOTE THAT dG/dy_i is NOT THE CHEMICAL POTENTIAL of i!'
             if(gx%bmperr.ne.0) goto 990
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
! Get current composition of the phase
             call calc_phase_molmass(iph,ics,xknown,aphl,totam,xxy,xxx,ceq)
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
2087            format(/'Calculated Gibbs energy/RT: ',1pe14.6,&
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
          case(6) !
             write(*,*)'Not implemeneted yet'
          END SELECT
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
!2101      format('N&x: ',F6.3,9F8.5)
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
! we should write phase tuples ...
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
          call meq_evaluate_all_svfun(-1,ceq)
! ignore error
          if(gx%bmperr.ne.0) gx%bmperr=0
          if(name1(1:1).eq.'*') then
! this calculate them again ... and lists the values
             call meq_evaluate_all_svfun(lut,ceq)
          else
             call capson(name1)
             call find_svfun(name1,istv,ceq)
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
          if(minimizer.eq.1) then
! Lukas minimizer, first argument=1 means use grid minimizer
!           call calceq1(1,ceq)
             write(kou,*)'Not implemented yet'
          else
             call calceq2(1,ceq)
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
             call gparcd('With gridminimizer? ',cline,last,1,ch1,'N',q1help)
! mode=0 is without grid minimizer 
             mode=1
             if(ch1.eq.'N' .or. ch1.eq.'n') mode=0
! Seach for memory leaks
             call gparid('How many times? ',cline,last,leak,1,q1help)
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
! OPENMP parallel start
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
                      write(kou,2050)neweq%eqno,neweq%eqname
2050                  format('Zero weight equilibrium ',i4,2x,a)
                   else
                      i2=i2+1
                      call calceq3(mode,.FALSE.,neweq)
                      if(gx%bmperr.ne.0) then
                         write(kou,2051)gx%bmperr,neweq%eqname,mode
2051                     format(' *** Error code ',i5,' for equilibrium ',&
                              a,' reset',i2)
                         gx%bmperr=0
                      elseif(idef.eq.1) then
                         write(lut,2052)neweq%eqno,&
                              neweq%eqname(1:len_trim(neweq%eqname)),&
                              neweq%tpval(1)
2052                     format('Calculated equilibrium ',i5,2x,a,', T=',F8.2)
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
! We calculate without grid minimizer, if parallel we must turn off
! creating/removing composition sets!! not safe to do that!!
!$             globaldata%status=ibset(globaldata%status,GSNOACS)
!$             globaldata%status=ibset(globaldata%status,GSNOREMCS)
!        !$OMP for an OMP directive
!        !$ as sentinel
! NOTE: $OMP  threadprivate(gx) declared in TPFUN4.F90 ??
!$OMP parallel do private(neweq)
                do i1=1,size(firstash%eqlista)
! the error code must be set to zero for each thread ?? !!
                   jp=jp+1
                   gx%bmperr=0
                   neweq=>firstash%eqlista(i1)%p1
                   if(neweq%weight.eq.zero) then
                      write(kou,2050)neweq%eqno,neweq%eqname
                   else
! write output only for idef=1
!$                     if(.TRUE. .and. idef.eq.1) then
!$                      write(*,663)'Equil/loop/thread/maxth/error: ',&
!$                             neweq%eqname,i1,omp_get_thread_num(),&
!$                             omp_get_num_threads(),gx%bmperr
663                   format(a,a,5i5)
! calceq3 gives no output
!$                        call calceq3(mode,.FALSE.,neweq)
!$                     else
! note first argument zero means do not use grid minimizer
                          call calceq3(mode,.FALSE.,neweq)
!$                     endif
                      i2=i2+1
                      if(gx%bmperr.ne.0) then
                         write(kou,2051)gx%bmperr,neweq%eqname,mode
                         write(*,*)'Error: ',gx%bmperr
                         gx%bmperr=0
                      elseif(idef.eq.1) then
                         write(lut,2052)neweq%eqno,&
                              neweq%eqname(1:len_trim(neweq%eqname)),&
                              neweq%tpval(1)
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
                enddo
!- $OMP end parallel do not needed???
! OPENMP parallel end loop
! allow composition sets to be created again
!$             globaldata%status=ibclr(globaldata%status,GSNOACS)
!$             globaldata%status=ibclr(globaldata%status,GSNOREMCS)
             endif gridmin
! repeat this until leak is zero, if leak negative never stop.
             leak=leak-1
             if(leak.ne.0) then
                call system_clock(count=ll)
                xxy=ll-j1; xxy=xxy/i2
                write(*,669)i2,ll-j1,xxy
669             format(/' *** Number of equlibria calculated ',2i10,F8.3/)
                goto 2060
             endif
!
             call system_clock(count=ll)
             call cpu_time(xxy)
             write(kou,664)i2,jp,xxy-xxx,ll-j1
664          format('Calculated ',i5,' equilibria out of ',i5/&
                  'Total time: ',1pe12.4,' s and ',i7,' clockcycles')
! this unit may have been used to extract calculated data for plotting
             if(plotunit0.gt.0) then
                write(kou,670)
670             format('Closing a GNUPLOT file oc_many0.plt'/&
                     'that may need some editing before plotting')
                write(plotunit0,665)
665             format('e'/'pause mouse'/)
                close(plotunit0)
! UNFINISHED possibly we could reopen the file again and make oopies 
! of tha data to avoid manual editing
             endif
          else
             write(kou,*)'You must first SET RANGE of experimental equilibria'
          endif
       END SELECT
!=================================================================
! SET SUBCOMMANDS
!         ['CONDITION       ','STATUS          ','ADVANCED        ',&
!         'LEVEL           ','INTERACTIVE     ','REFERENCE_STATE ',&
!         'QUIT            ','ECHO            ','PHASE           ',&
!         'UNITS           ','LOG_FILE        ','WEIGHT          ',&
!         'NUMERIC_OPTIONS ','AXIS            ','INPUT_AMOUNTS   ',&
!         'VERBOSE         ','AS_START_EQUILIB','BIT             ',&
!         'VARIABLE_COEFF  ','SCALED_COEFF    ','OPTIMIZING_COND ',&
!         'EXPERIMENT_EQUIL','FIXED_COEFF     ','                ']
    CASE(3) ! SET SUBCOMMANDS
! disable continue optimization
       iexit=0
       iexit(2)=1
       kom2=submenu(cbas(kom),cline,last,cset,ncset,1)
       if(kom2.le.0) goto 100
       SELECT CASE(kom2)
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
          SELECT CASE(kom3)
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
          END SELECT
!-----------------------------------------------------------
       case(3) ! set ADVANCED
! default is DENSE_GRID
          name1='advanced command'
          kom3=submenu(name1,cline,last,cadv,ncadv,4)
          select case(kom3)
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
          case(3) ! SET ADVANCED EXTRTA_PROPERTY for a species
             call gparc('Species symbol: ',cline,last,1,name1,' ',q1help)
             call find_species_record(name1,loksp)
             if(gx%bmperr.ne.0) goto 100
             xxy=one
             call gparrd('Property value: ',cline,last,xxx,xxy,q1help)
             if(buperr.eq.0) then
                call enter_species_property(loksp,xxx)
             endif
!.................................................................
          case(4) ! SET ADVANCED DENSE_GRID_ONOFF
! this sets bit 14 of global status word, also if bit 2 (expert) not set
             if(btest(globaldata%status,GSXGRID)) then
                globaldata%status=ibclr(globaldata%status,GSXGRID)
                write(*,3110)'reset'
3110            format('Dense grid ',a)
             else
                globaldata%status=ibset(globaldata%status,GSXGRID)
                write(*,3110)'dense grid set'
             endif
!.................................................................
          case(5) ! nothing yet
             write(*,*)'Not implemented yet'
!.................................................................
          case(6) ! nothing yet
             write(*,*)'Not implemented yet'
          end select
!-----------------------------------------------------------
       case(4) ! set LEVEL, not sure what it will be used for ...
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! end of macro excution (can be nested)
       case(5) ! set INTERACTIVE
          call macend(cline,last,logok)  
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
          SELECT CASE(kom3)
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
             SELECT CASE(kom4)
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
                call set_phase_status_bit(lokph,PHFORD)
             case(2) ! BCC_PERMUTATIONS BORD
                if(check_minimal_ford(lokph)) goto 100
                call set_phase_status_bit(lokph,PHBORD)
             case(3) ! IONIC_LIQUID_MDL this may require tests and 
! other bits changed ..
                write(kou,*)'Cannot be set interactivly yet, only from TDB'
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
             end SELECT
!............................................................
          case(6) ! SET PHASE ... CONSTITUTION iph and ics set above
             call ask_phase_new_constitution(cline,last,iph,ics,lokcs,ceq)
       END SELECT
!-------------------------------------------------------------
       case(10) ! set UNIT (for state variables)
          write(kou,*)'Not implemented yet'
!-------------------------------------------------------------
       case(11) ! set LOG_FILE
          call gparcd('Log file name: ',cline,last,1,name1,'oclog',q1help)
          call capson(name1)
          if(name1(1:5).eq.'NONE ') then
             call openlogfile(' ',' ',-1)
             logfil=0
          else
             call gparc('Title: ',cline,last,5,model,' ',q1help)
             call openlogfile(name1,model,39)
             if(buperr.ne.0) then
                write(kou,*)'Error opening logifile: ',buperr
                logfil=0
             else
                logfil=39
             endif
          endif
!-------------------------------------------------------------
       case(12) ! set weight
          if(.not.allocated(firstash%eqlista)) then
             write(kou,*)'There are no experimental equilibria'
             goto 100
          endif
          call gparrd('Weight ',cline,last,xxx,one,q1help)
          if(buperr.ne.0) goto 100
          call gparcd('Equilibria (abbrev name) ',cline,last,&
               1,name1,'CURRENT',q1help)
          if(name1(1:8).eq.'CURRENT ') then
             if(ceq%eqname(1:20).eq.'DEFAULT_EQUILIBRIUM ') then
                write(kou,*)'You cannot change the weight for the default'
             else
                ceq%weight=abs(xxx)
! mark that list optimize may not work
                mexp=0
             endif
          elseif(name1(1:1).eq.'*') then
! set this weight for all
             i2=0
             do i1=1,size(firstash%eqlista)
                firstash%eqlista(i1)%p1%weight=xxx
                i2=i2+1
             enddo
             write(kou,*)'Changed weight for ',i2,' equilibria'
          else
! set this weight to all equilibria with name abbriviations fitting name1
             call capson(name1)
             i2=0
             do i1=1,size(firstash%eqlista)
                if(index(firstash%eqlista(i1)%p1%eqname,&
                     name1(1:len_trim(name1))).gt.0) then
                   firstash%eqlista(i1)%p1%weight=xxx
                   i2=i2+1
! mark that list optimize may not work
                   mexp=0
                endif
             enddo
             write(kou,*)'Changed weight for ',i2,' equilibria'
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
          xxx=ceq%xconv
          call gparrd('Max error in fraction: ',cline,last,xxy,xxx,q1help)
          if(xxy.gt.1.0D-30) then
             ceq%xconv=xxy
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
!          elseif(iax.lt.i1) then
! replacing an existing axis, reset defaults ... ???
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
             starteq%next=neweq%eqno
          else
             starteq=>neweq
             starteq%next=0
             write(*,*)'Starteq next',starteq%next
          endif
          write(*,*)'A copy of current equilibrium linked as start eqilibrium'
!-------------------------
       case(18) ! SET BIT (all kinds of bits) just global implemented
!          call gparid('Toggle bit (from 0-31, -1 quits):',&
!               cline,last,ll,-1,q1help)
!         ['EQUILIBRIUM     ','GLOBAL          ','PHASE           ',&
          kom3=submenu('Set which status word?',cline,last,csetbit,nsetbit,2)
          SELECT CASE(kom3)
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
3612            format('Toggles equilibrium status word, currently: ',z8,/&
                     'Bit If set means',/&
                     ' 0  No threads allowed (no parallel calculation)',/&
                     ' 1  No global minimization allowed',/&
                     ' 2  No equilibrium has been calculated',/&
                     ' 3  Conditions and results not consistent',/'-'/&
                     ' 4  Last equilibrium calculation failed',/&
                     ' 5  No automatic generation of composition sets',/&
                     ' 6  Equilibrim tested by gridminimizer',/&
                     ' 7  Current results are from a grid minimization'/)
                goto 3610
             endif
             if(ll.lt.0 .or. ll.gt.7) then
                write(kou,*)'No such bit, no bit changed'
             else
                if(btest(ceq%status,ll)) then
                   ceq%status=ibclr(ceq%status,ll)
                   write(kou,*)'Bit cleared'
                else
                   ceq%status=ibset(ceq%status,ll)
                   write(kou,*)'Bit set'
                endif
             endif
!             write(*,*)'Not implemented yet'
!................................................................
! maybe change order of questions, maybe check name exits etc ....
          case(2) ! global status word
3708         continue
! subroutine TOPHLP forces return with ? in position cline(last:last)
             write(kou,3709)globaldata%status
3709         format('Current global status word: ',z8)
             call gparid('Toggle global status bit (from 0-31, -1 quits):',&
                  cline,last,ll,-1,tophlp)
             if(cline(1:1).eq.'?') then
                write(kou,3710)globaldata%status
3710            format('Toggles global status word ',z8,&
                     ' (only experts should change these) '/&
                     'Bit If set means:'/&
                     ' 0  user is a beginner'/&
                     ' 1  user is experienced'/&
                     ' 2  user is an expert'/&
                     ' 3  gridminimizer will not be used'/'-'/&
                     ' 4  gridminimizer must not merge comp.sets.'/&
                     ' 5  there are no data'/&
                     ' 6  there are no phases'/&
                     ' 7  comp.sets must not be created automatically'/'-'/&
                     ' 8  comp.sets must not be deleted automatically'/&
                     ' 9  data has changed since last save'/&
                     '10  verbose is on'/&
                     '11  verbose is permanently on'/'-'/&
                     '12  no warnings'/&
                     '13  no cleanup after an equilibrium calculation'/&
                     '14  denser grid used in grid minimizer'/&
                     '15  calculations in parallel is not allowed'/'-'/&
                     '16  no global test at node points durung STEP/MAP'/&
                     '17  the components are not the elements'/&
                     '18  test calcúlated equilibrium with the grid minimizer')
                goto 3708
             endif
             if(ll.lt.0 .or. ll.gt.31) then
                write(kou,*)'No bit changed'
             elseif(btest(globaldata%status,GSADV) .or. ll.le.2) then
! user must have expert bit set to change any other bit than the user type bit
                if(btest(globaldata%status,ll)) then
                   globaldata%status=ibclr(globaldata%status,ll)
                   write(*,3711)'cleared',globaldata%status
3711               format('Bit ',a,', new value of status word: ',z8)
                else
                   globaldata%status=ibset(globaldata%status,ll)
                   write(*,3711)'set',globaldata%status
                endif
                if(.not.btest(globaldata%status,GSADV)) then
! if expert/experienced bit is cleared ensure that experienced bit is set
                   globaldata%status=ibset(globaldata%status,GSOCC)
                endif
             else
                write(kou,*)'You must have expert status to toggle this!'
             endif
!....................................................
          case(3) ! set bit phase ...
             write(*,*)'Please use set phase ... bit '
          end select
!-------------------------
       case(19) ! set variable_coefficent, 0 to 99
          if(.not.btest(firstash%status,AHCOEF)) then
             write(kou,*)'No optimizing coefficents'
             goto 100
          endif
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
          xxy=firstash%coeffvalues(j1)*firstash%coeffscale(j1)
          if(i1.eq.i2) then
! when setting a single coefficinet variable ask for value
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
! set coefficient variable with current value
             xxx=xxy
          endif
! set new value
          call change_optcoeff(firstash%coeffindex(j1),xxx)
          if(gx%bmperr.ne.0) goto 100
          firstash%coeffvalues(j1)=one
          firstash%coeffscale(j1)=xxx
          firstash%coeffstart(j1)=xxx
! mark an optimized coefficient without min/max
          if(firstash%coeffstate(j1).lt.10) then
             nvcoeff=nvcoeff+1
          endif
          firstash%coeffstate(j1)=10
          if(i2.gt.j1) then
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
!-------------------------
       case(21) ! set optimizing_conditions, more will be added
          call gparrd('DSTEP (VA05AD): ',cline,last,xxx,dstep,q1help)
          dstep=xxx
          call gparrd('DMAX (VA05AD): ',cline,last,xxx,dmax2,q1help)
          dmax2=xxx
          call gparrd('ACC (VA05AD): ',cline,last,xxx,acc,q1help)
          acc=xxx
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
                write(plotdataunit(i1),22)
22              format('e'/'pause mouse'/)
                close(plotdataunit(i1))
                plotdataunit(i1)=0
             endif
          enddo
!          write(*,*)'Not implemeneted yet'
!-------------------------
       case(23) ! set fixed_coefficient
          if(.not.btest(firstash%status,AHCOEF)) then
             write(kou,*)'No optimizing coefficents'
             goto 100
          endif
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
       case(24) ! unused
          write(*,*)'Not implemeneted yet'
       END SELECT
!=================================================================
! enter with subcommand for element, species etc
!         ['TPFUN_SYMBOL    ','ELEMENT         ','SPECIES         ',&
!         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
!         'CONSTITUTION    ','EXPERIMENT      ','QUIT            ',&
!         'EQUILIBRIUM     ','SYMBOL          ','OPTIMIZE_COEFF  ',&
!         'COPY_OF_EQUILIB ','COMMENT         ','MANY_EQUILIBRIA ',&
!         'MATERIAL        ','                ','                ']
    CASE(4)
! disable continue optimization
       iexit=0
       iexit(2)=1
       kom2=submenu(cbas(kom),cline,last,center,ncent,11)
       SELECT CASE(kom2)
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
          call gparc('Element full name: ',cline,last,1,name1,' ',q1help)
          call gparc('Element reference phase: ',cline,last,1,name2,' ',q1help)
          call gparrd('Element mass (g/mol): ',cline,last,mass,one,q1help)
          if(buperr.ne.0) goto 990
          call gparr('Element H298-H0: ',cline,last,h298,zero,q1help)
          if(buperr.ne.0) goto 990
          call gparr('Element S298: ',cline,last,s298,one,q1help)
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
       case(5) ! enter parameter is always allowed
          call enter_parameter_interactivly(cline,last,0)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(6) ! enter bibliography
          call enter_reference_interactivly(cline,last,0,j1)
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
          call gparc('Name: ',cline,last,1,text,' ',q1help)
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
             allocate(firstash%coeffscale(0:i1))
             allocate(firstash%coeffstart(0:i1))
             allocate(firstash%coeffmin(0:i1))
             allocate(firstash%coeffmax(0:i1))
             allocate(firstash%coeffindex(0:i1))
             allocate(firstash%coeffstate(0:i1))
! coeffvalues should be of the order of one
             firstash%coeffvalues=one
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
553          format('You have already ',i3,' optimizing coefficents entered')
          endif
!          write(*,551)firstash%status
!551       format('Assessment status word: ',z8)
!---------------------------------------------------------------
! enter copy of equilibrium (for test!)
       case(13)
          call gparc('Name: ',cline,last,1,text,' ',q1help)
          if(buperr.ne.0) goto 100
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
! enter not used
       case(18)
          write(*,*)'Not implemeneted yet'
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
       END SELECT
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
!         'OPTIMIZATION    ','MODEL_PARAM_VAL ','                ']
    CASE(6) ! LIST
       kom2=submenu(cbas(kom),cline,last,clist,nclist,12)
       if(kom2.le.0) goto 100
       lut=optionsset%lut
       SELECT CASE(kom2)
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
       case(2) ! list short with status bits
          call gparcd('Option (A/C/M/P)',cline,last,1,ch1,chshort,q1help)
          call capson(ch1)
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
          SELECT CASE(kom3)
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
             if(gx%bmperr.ne.0) goto 990
!...............................................................
          case(3) ! list phase model (including disordered fractions)
             write(kou,6070)'For ',ceq%eqno,ceq%eqname
6070      format(a,'equilibrium: ',i3,', ',a)
             call list_phase_model(iph,ics,lut,' ',ceq)
          END SELECT
!------------------------------
       case(4,17)  ! list state variable or parameter identifier value, loop.
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
6105      continue
! NOTE: 4th argument is 5 because otherwise , will terminate reading cline
! and state variables like x(fcc,cr) will not work.
          if(kom2.eq.4) then
             call gparc('State variable: ',cline,last,5,line,' ',q1help)
          else
             if(once) then
                write(kou,*)'Remember always to specify the phase!'
                once=.FALSE.
             endif
             call gparc('Parameter ident: ',cline,last,5,line,' ',q1help)
          endif
          j1=1
          if(eolch(line,j1)) goto 100
          model=' '
          if(index(line,'*').gt.0) then
! generate many values
! i1 values are resturned in yarr with dimension maxconst. 
! longstring are the state variable symbols for the values ...
             call get_many_svar(line,yarr,maxconst,i1,longstring,ceq)
             if(gx%bmperr.eq.0) then
! not a nice output FIX!!
                write(lut,*)longstring(1:len_trim(longstring))
                write(lut,6107)(yarr(i2),i2=1,i1)
6107            format('Values: ',5(1pe14.6)/(8x,5(1pe14.6)))
             endif
          else
! in model the state variable is returned as generated by the program
             call get_state_var_value(line,xxx,model,ceq)
             if(gx%bmperr.eq.0) then
                write(lut,6108)model(1:len_trim(model)),xxx
6108            format(1x,a,'=',1PE15.7)
             endif
          endif
          if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
             write(lut,*)bmperrmess(gx%bmperr)
          elseif(gx%bmperr.ne.0) then
             write(lut,*)'Error code ',gx%bmperr
          endif
          gx%bmperr=0
          goto 6105
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
       case(10) ! list parameter for a phase (just one)
          call enter_parameter_interactivly(cline,last,1)
!-----------------------------------------------------------
       case(11) ! list equilibria (not result)
          do iel=1,noeq()
             if(associated(ceq,eqlista(iel))) then
                name1='**'
             else
                name1=' '
             endif
!             j1=len_trim(eqlista(iel)%comment)
             if(eqlista(iel)%weight.le.zero) then
                write(lut,6202)iel,name1(1:2),eqlista(iel)%eqname,&
                     eqlista(iel)%tpval(1),eqlista(iel)%comment(1:33)
6202            format(i4,1x,a2,1x,a,' T=',F8.2,', ',a)
             else
                write(lut,6203)iel,name1(1:2),eqlista(iel)%eqname,&
                     eqlista(iel)%tpval(1),eqlista(iel)%weight,&
                     eqlista(iel)%comment(1:20)
6203            format(i4,1x,a2,1x,a,' T=',F8.2,', weight=',F5.2,', ',a)
             endif
!             if(j1.gt.1) then
!                write(lut,6204)eqlista(iel)%comment(1:j1)
!6204            format(12x,a)
!             endif
          enddo
!------------------------------
       case(12) ! list results
! if no calculation made skip
          if(btest(ceq%status,EQNOEQCAL)) then
             write(*,6277)ceq%status
6277         format(' *** No results as no equilibrium calculated! ',z8)
             goto 100
          elseif(btest(ceq%status,EQGRIDCAL)) then
             write(kou,*)' *** Last calculation was not a full equilibrium'
          endif
          call gparid('Output mode: ',cline,last,listresopt,lrodef,q1help)
          if(buperr.ne.0) then
             write(kou,*)'No such mode, using default'
             buperr=0
             listresopt=lrodef
          endif
          if(listresopt.gt.0 .and. listresopt.le.9) then
             lrodef=listresopt
          endif
          call date_and_time(optres,name1)
          write(lut,6051)ceq%eqno,ceq%eqname,optres(1:4),optres(5:6),optres(7:8)
! write comment line if any
          if(len_trim(ceq%comment).gt.0) then
             write(lut,6308)trim(ceq%comment)
6308         format(3x,a)
          endif
          if(btest(ceq%status,EQFAIL)) then
             write(lut,6305)'below'
6305         format(/' *** The listing ',a,&
                  ' are not a valid equilibrium as last calculation failed'/)
          elseif(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'No results as no data'
             goto 100
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
          call list_stored_equilibria(lut,axarr,maptop)
!------------------------------
! list optimization, several suboptions
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
               name1(1:2),name1(3:4)
600       format(/'Listing of optimization results: date ',a4,'.',a2,'.',a2,&
               ' : ',a2,'h',a2)
          SELECT CASE(kom2)
!..........................................................
             case DEFAULT
                write(kou,*)'No such option'
!...........................................................
             case(1) ! short
               call listoptshort(lut,mexp,errs)
!...........................................................
             case(2) ! long
                write(*,*)'Not implemented yet'
!...........................................................
             case(3) ! coefficent values
                call listoptcoeff(lut)
!...........................................................
             case(4) ! graphics
                write(*,*)'Not implemented yet'
!...........................................................
             case(5) ! debug
                write(*,*)'Not implemented yet'
!...........................................................
             case(6) ! MACRO include experiments
                write(*,*)'Not implemented yet'
!...........................................................
             case(7) ! unused
                write(*,*)'Not implemented yet'
!...........................................................
             case(8) ! unused
                write(*,*)'Not implemented yet'
!...........................................................
             case(9) ! unused
                write(*,*)'Not implemented yet'
             end SELECT
!------------------------------
! list model parameter values, part of case(4)
!       case(17)
!          write(*,*)'Not implemented yet'
!------------------------------
! list error message
       case(18)
          i2=4204
          call gparid('Error code: ',cline,last,i1,i2,q1help)
          if(i1.ge.4000 .and. i1.lt.4399) then
             write(kou,4999)i1,bmperrmess(i1)
4999         format('The error code ',i4', means: '/a)
          else
             write(kou,*)'No a standard OC error message'
          endif
!------------------------------
! list ??
       case(19)
          write(*,*)'Not implemented yet'
!------------------------------
! list ??
       case(20)
          write(*,*)'Not implemented yet'
!------------------------------
! list ??
       case(21)
          write(*,*)'Not implemented yet'
       end SELECT
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
! read subcommand
!        ['UNFORMATTED     ','TDB             ','QUIT            ',&
!         'DIRECT          ','                ','                ']
    case(8)
! disable continue optimization
       iexit=0
       iexit(2)=1
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
       SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'Read subcommand error'
!-----------------------------------------------------------
       case(1) ! read unformatted file created by SAVE
          if(ocufile(1:1).ne.' ') then
             text=ocufile
             call gparcd('File name: ',cline,last,1,ocufile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,ocufile,' ',q1help)
          endif
          call gtpread(ocufile,text)
          if(gx%bmperr.ne.0) goto 990
          kl=len_trim(text)
          if(kl.gt.1) then
             write(kou,8110)text(1:kl)
          endif
8110      format(/'Savefile text: ',a/)
! if there is an assessment record set nvcoeff ...
          if(allocated(firstash%eqlista)) then
             write(*,*)'There is an assessment record'
             nvcoeff=0
             kl=size(firstash%coeffvalues)-1
             do j1=0,kl
                if(firstash%coeffstate(j1).ge.10) then
                   nvcoeff=nvcoeff+1
                endif
             enddo
             write(kou,3730)nvcoeff
          endif
!---------------------------------------------------------
       case(2) ! read TDB
          if(tdbfile(1:1).ne.' ') then
! set tdbfil as default
             text=tdbfile
             call gparcd('File name: ',cline,last,1,tdbfile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,tdbfile,' ',q1help)
          endif
! if tdbfle starts with "ocbase/" replace that with content of ocbase!!
          name1=tdbfile(1:7)
          call capson(name1)
          if(name1(1:7).eq.'OCBASE/' .or. name1(1:7).eq.'OCBASE\') then
             tdbfile=trim(ocbase)//tdbfile(7:)
             write(*,*)'database file: ',trim(tdbfile)
          endif
! this call checks the file exists and returns the elements
          call checkdb(tdbfile,'.tdb',jp,ellist)
          if(gx%bmperr.ne.0) then
             write(kou,*)'No database with this name'
             goto 990
          elseif(jp.eq.0) then
             write(Kou,*)'No elements in the database'
          endif
          write(kou,8203)jp,(ellist(kl),kl=1,jp)
8203      format('Database has ',i2,' elements: ',18(a,1x)/(1x,28(1x,a)))
          write(kou,8205)
8205      format('Give the elements to select, finish with empty line')
          jp=1
          selection='Select elements /all/:'
8210      continue
!          call gparc('Select elements/all/: ',&
          call gparc(selection,cline,last,1,ellist(jp),' ',q1help)
          if(ellist(jp).ne.'  ') then
             call capson(ellist(jp))
             jp=jp+1
             if(jp.gt.size(ellist)) then
                write(kou,*)'Max number of elements selected: ',size(ellist)
             else
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
       case(5) ! read PDB 
          if(tdbfile(1:1).ne.' ') then
             text=tdbfile
             call gparcd('File name: ',cline,last,1,tdbfile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,tdbfile,' ',q1help)
          endif
! this call checks the file exists and returns the elements
          call checkdb(tdbfile,'.pdb',jp,ellist)
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
       end SELECT
!=================================================================
! save in various formats (NOT TDB, MACRO and LATEX, use LIST DATA)
! It is a bit inconsistent as one READ TDB but not SAVE TDB ...
!        ['TDB             ','SOLGASMIX       ','UNFORMATTED     ',&
!         'DIRECT          ','                ','QUIT            ']
    CASE(9)
! default i 3, unformatted
       kom2=submenu(cbas(kom),cline,last,csave,ncsave,3)
       if(kom2.le.0 .or. kom2.gt.ncsave) goto 100
!
       if(kom2.gt.2) then
! Do not ask this question for TDB and SOLGASMIX files
          call gparc('Comment line: ',cline,last,5,model,' ',q1help)
       endif
       SELECT CASE(kom2)
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
       case(2) ! SOLGASMIX
          text=' '
          call gparc('File name: ',cline,last,1,filename,text,q1help)
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
          call save_datformat(filename,kl,ceq)   
!-----------------------------------------------------------
       case(3) ! save unformatted
132       continue
          if(ocufile(1:1).ne.' ') then
             text=ocufile
             call gparcd('File name: ',cline,last,1,ocufile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,ocufile,' ',q1help)
          endif
          jp=0
          kl=index(ocufile,'.')
          if(kl.le.0) then
             jp=len_trim(ocufile)
          elseif(ocufile(kl+1:kl+1).eq.' ') then
! just ending a filename with . not accepted as extention
             jp=kl
          endif
          if(kl.le.0 .and. jp.le.0) then
             write(kou,*)'Missing file name, nothing saved'
             goto 100
          endif
          if(jp.gt.0) ocufile(jp+1:)='.ocu '
          inquire(file=ocufile,exist=logok)
          if(logok) then
             call gparcd('File exists, overwrite?',cline,last,1,ch1,'N',q1help)
             if(ch1.ne.'Y') then
                write(*,133)
133             format('You can use another file name')
                ocufile=' '
                goto 132
             endif
             write(*,134)trim(ocufile)
134          format(/'Overwriting previous results on ',a/)
          endif
          text='U '//model
          call gtpsave(ocufile,text)
!-----------------------------------------------------------
       case(4) ! save DIRECT
          if(ocdfile(1:1).ne.' ') then
             text=ocdfile
             call gparcd('File name: ',cline,last,1,ocdfile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,ocdfile,' ',q1help)
          endif
          jp=0
          kl=index(ocdfile,'.')
          if(kl.le.0) then
             jp=len_trim(ocdfile)
          elseif(ocdfile(kl+1:kl+1).eq.' ') then
! just ending a filename with . not accepted as extention
             jp=kl
          endif
          if(kl.le.0 .and. jp.le.0) then
             write(kou,*)'Missing file name, nothing saved'
             goto 100
          endif
          if(jp.gt.0) ocdfile(jp+1:)='.ocd '
          text='M '//model
          call gtpsave(ocdfile,text)
!-----------------------------------------------------------
       case(5) ! 
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
       case(6) ! save quit, do nothing
          continue
       end SELECT
!=================================================================
! help ... just list the commands
    case(10)
       call q3help(cline,last,cbas,ncbas)
       goto 100
!=================================================================
! information
    case(11)
       write(kou,*)'Not implemented yet'
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
          write(*,*)' *** Warning, assessment data not removed'
       endif
!       write(*,*)'No segmentation fault 2'
       if(allocated(firstash%eqlista)) deallocate(firstash%eqlista)
       deallocate(firstash)
!       write(*,*)'No segmentation fault 3'
       mexp=0
       nvcoeff=0
       iexit=0
       iexit(2)=1
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
! remove any results from step and map
       if(associated(maptop)) then
!          write(*,*)'maptop nullified: ',maptop%next%seqx
          maptop%next%seqx=0
          maptop%next%seqy=0
          nullify(maptop)
       endif
       nullify(mapnode)
       nullify(maptopsave)
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
       graphopt%rangedefaults=0
       graphopt%tielines=0
       graphopt%status=0
       graphopt%axistype=0
       graphopt%tielines=0
       graphopt%gibbstriangle=.FALSE.
       graphopt%labelkey='top right'
       graphopt%appendfile=' '
       nullify(graphopt%firsttextlabel)
!
! this routine fragile, inside new_gtp init_gtp is called
!       write(*,*)'No segmentation fault 7'
       call new_gtp
       write(kou,*)'All data removed'
       call init_gtp(intv,dblv)
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
       write(kou,15010)linkdate
15010  format(/'This is OpenCalphad (OC), a free software for ',&
            'thermodynamic calculations'/&
            'described by B Sundman, U R Kattner, M Palumbo and S G Fries, ',&
            'Integrating'/'Materials and Manu Innov (2015) 4:1 and ',&
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
            'Copyright 2011-2016, Bo Sundman, France.'/&
            'Contact person Bo Sundman, bo.sundman@gmail.com'/&
            'This version linked ',a/)
!=================================================================
! debug subcommands
    case(16)
       kom2=submenu(cbas(kom),cline,last,cdebug,ncdebug,1)
       SELECT CASE (kom2)
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
               'Tuple lokph   compset ixphase lokvares nextcs phase name')
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
!          call delete_all_conditions(0,ceq)
!---------------------------------
! debug test2 (whatever)
       case(5)
          write(*,*)'Nothing here'
!---------------------------------
! debug unused
       case(6)
          write(*,*)'Neither here'
       END SELECT
!=================================================================
! select command
    case(17)
       kom2=submenu(cbas(kom),cline,last,cselect,nselect,1)
       SELECT CASE(kom2)
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
             call gparcd('Do you want to use VA05AD?',&
                  cline,last,1,ch1,'Y',q1help)
             if(ch1.eq.'Y') then
                optimizer=2
             endif
          endif
          write(kou,*)'You have selected ',optimizers(optimizer)
!-----------------------------------------------------------
       case(6)
          goto 100
       END SELECT
!=================================================================
! delete not much implemented ...
!         ['ELEMENTS        ','SPECIES         ','PHASE           ',&
!          'QUIT            ','COMPOSITION_SET ','EQUILIBRIUM     ',&
!          'STEP_MAP_RESULTS','                ','                ']
    CASE(18)
! disable continue optimization
       iexit=0
       iexit(2)=1
       kom2=submenu(cbas(kom),cline,last,crej,nrej,3)
       SELECT CASE(kom2)
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
          if(associated(maptop)) then
             write(*,*)'maptop nullified: ',maptop%next%seqx
             maptop%next%seqx=0
             maptop%next%seqy=0
             nullify(maptop)
          endif
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
          graphopt%rangedefaults=0
          graphopt%tielines=0
          graphopt%status=0
          graphopt%axistype=0
          graphopt%tielines=0
          graphopt%gibbstriangle=.FALSE.
          graphopt%labelkey='upper right'
          graphopt%appendfile=' '
          nullify(graphopt%firsttextlabel)
!-----------------------------------------------------------
!
       case(8)
          continue
!-----------------------------------------------------------
!
       case(9)
          continue
       end SELECT
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
! check if adding results
       if(associated(maptop)) then
          write(kou,833)
833       format('There are previous results from step or map')
          call gparcd('Delete them?',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
! there should be a more careful deallocation to free memory
             deallocate(maptop%saveceq)
             nullify(maptop)
             nullify(maptopsave)
             write(kou,*)'Previous results removed'
! delete equilibria associated with STEP/MAP
             call delete_equilibria('_MAP*',ceq)
             seqxyz=0
          else
             seqxyz(1)=maptop%next%seqx
             seqxyz(2)=maptop%seqy
             maptopsave=>maptop
             nullify(maptop)
!             write(kou,*)'Previous results kept'
          endif
       endif
       kom2=submenu('Options?',cline,last,cstepop,nstepop,1)
       SELECT CASE(kom2)
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
! can one have several STEP commands??
          if(associated(maptop)) then
             write(*,*)'Deleting previous step/map results missing'
          endif
! seqzyz are initial values for creating equilibria for lines and nodes
! if previous results should be kept it should not be zeroed, see above
!          seqxyz=0
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
          endif
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
          endif
!          seqxyz=0
          call step_separate(maptop,noofaxis,axarr,seqxyz,starteq)
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
          endif
! set default yaxis as GM(*)
          if(axplotdef(2)(1:1).eq.' ') then
             axplotdef(2)='GM(*)'
          endif
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
       end SELECT
!=================================================================
! MAP, must be tested if compatible with assessments
    case(20)
! disable continue optimization
!       iexit=0
!       iexit(2)=1
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
             deallocate(maptop%saveceq)
             nullify(maptop)
             nullify(maptopsave)
! this removes all previous equilibria associated with STEP/MAP commands
             call delete_equilibria('_MAP*',ceq)
             if(gx%bmperr.ne.0) then
                write(kou,*)'Error removing old MAP equilibria'
                goto 990
             endif
! initiate indexing nodes and lines
             seqxyz=0
          else
! start indexing new noes/lines from previous 
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
          starteq%next=0
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
          write(Kou,*)'Set link to previous map results'
          maptop%plotlink=>maptopsave
          nullify(maptopsave)
       endif
! remove start equilibria
       nullify(starteq)
! mark that interactive listing of conditions and results may be inconsistent
       ceq%status=ibset(ceq%status,EQINCON)
       if(gx%bmperr.ne.0) goto 990
!=================================================================
! PLOT
    case(21)
       if(.not.associated(maptop)) then
          write(kou,*)'You must give a STEP or MAP command before PLOT'
          goto 100
       endif
       wildcard=.FALSE.
       do iax=1,2
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
! return here after each subcommand
21000      continue
          if(iax.eq.1) then
             call gparcd('Horizontal axis variable',&
                  cline,last,7,axplot(iax),axplotdef(iax),q1help)
          else
             call gparcd('Vertical axis variable',&
                  cline,last,7,axplot(iax),axplotdef(iax),q1help)
          endif
          if(buperr.ne.0) goto 990
          if(index(axplot(iax),'*').gt.0) then
!             if(wildcard) then
!                write(*,*)'Wildcards allowed for one axis only'
!                goto 21000
!             else
                wildcard=.TRUE.
!             endif
          endif
          if(axplotdef(iax).ne.axplot(iax)) then
! RESTORE DEFAULTS if not same axis variables !!! 
             graphopt%rangedefaults(iax)=0
             graphopt%appendfile=' '
             graphopt%gibbstriangle=.FALSE.
             graphopt%labelkey='top right'
! more options to restore ...
          endif
! remember axis as default
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
! plot options subcommand, default is PLOT, NONE does not work ...
21100   continue
       if(plotform(1:1).eq.'P') then
          write(kou,21110)plotfile(1:len_trim(plotfile))
21110     format(/' *** Graphics format is PS on: ',a,'.ps ')
       elseif(plotform(1:1).eq.'G') then
          write(kou,21111)plotfile(1:len_trim(plotfile))
21111     format(/' *** Graphics format is GIF on: ',a,'.gif ')
       elseif(plotform(1:1).eq.'A') then
          write(kou,21113)plotfile(1:len_trim(plotfile))
21113     format(/' *** Graphics format is PDF on: ',a,'.pdf ')
       endif
       write(kou,21112)
21112  format(/'Note: give only one option per line!')
       kom2=submenu('Options?',cline,last,cplot,nplt,1)
       SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'No such plot option'
!-----------------------------------------------------------
! RENDER no more options to plot ...
       case(1)
! added ceq in the call to make it possible to handle change of reference states
!2190      continue
          call ocplot2(jp,axplot,plotfile,maptop,axarr,graphopt,plotform,ceq)
          if(gx%bmperr.ne.0) goto 990
! always restore default plot file name!!
          plotfile='ocgnu'
!          call gparc('Hardcopy (P for postscript)?',&
!               cline,last,1,form,'none',q1help)
!          if(form.ne.'none') then
!             call capson(form)
!             call ocplot2(jp,axplot,plotfile,maptop,axarr,graphopt,form,ceq)
!          endif
!-----------------------------------------------------------
! XRANGE
       case(2)
          call gparcd('Default limits',cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'Y' .or. ch1.eq.'y') then
             graphopt%rangedefaults(1)=0
          else
             graphopt%rangedefaults(1)=1
             twice=.FALSE.
21104        continue
             call gparrd('Low limit',cline,last,xxx,graphopt%dfltmin(1),q1help)
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
!-----------------------------------------------------------
! YRANGE
       case(3)
          call gparcd('Default limits',cline,last,1,ch1,'N',q1help)
          if(ch1.eq.'Y' .or. ch1.eq.'y') then
             graphopt%rangedefaults(2)=0
          else
             graphopt%rangedefaults(2)=1
             twice=.FALSE.
21107        continue
             call gparrd('Low limit',cline,last,xxx,graphopt%dfltmin(2),q1help)
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
! XTEXT
       case(4)
          write(*,*)'Not implemented yet'
          goto 21100
!-----------------------------------------------------------
! YTEXT
       case(5)
          write(*,*)'Not implemented yet'
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
             call gparcd('Graphics format (ACROBAT(PDF)/PS/GIF/SCREEEN)',&
                  cline,last,1,ch1,'SCREEN',q1help)
             if(ch1.eq.'a' .or. ch1.eq.'A') then
                write(kou,*)'Graphics format set to ACROBAT PDF'
                plotform='A'
             elseif(ch1.eq.'p' .or. ch1.eq.'P') then
                write(kou,*)'Graphics format set to PS'
                plotform='P'
             elseif(ch1.eq.'g' .or. ch1.eq.'G') then
                write(kou,*)'Graphics format set to GIF'
                plotform='G'
             else
                write(kou,*)'Graphics format set to SCREEN'
                plotform=' '
             endif
          endif
!-----------------------------------------------------------
! PLOT OUTPUT_FILE, always asked when changing terminal
          call gparcd('Plot file',cline,last,1,plotfile,'ocgnu',q1help)
          goto 21100
!-----------------------------------------------------------
! PLOT GIBBS_TRIANGLE
       case(9)
          write(*,*)'Not implemented yet'
!          call gparcd('Triangular diagram?',cline,last,5,ch1,'NO',q1help)
!          if(ch1.eq.'y') then
!             graphopt%gibbstriangle=.TRUE.
!          else
             graphopt%gibbstriangle=.FALSE.
!          endif
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
          goto 21100
!-----------------------------------------------------------
! PLOT APPEND a gnuplot file
       case(12)
          write(kou,*)'Give a file name with graphics in GNUPLOT format'
          call gparcd('File name',cline,last,1,text,'  ',q1help)
! check it is OK and add .plt if necessary ...
          jp=index(text,'.plt ')
          if(jp.le.0) then
             jp=len_trim(text)
             text(jp+1:)='.plt'
          endif
          open(23,file=text,status='old',access='sequential',err=21300)
          close(23)
          graphopt%appendfile=text
          goto 21100
! error opening file, remove any previous appended file
21300     continue
          if(graphopt%appendfile(1:1).ne.' ') then
             write(*,21304)trim(graphopt%appendfile)
21304        format('Removing append file: ',a)
          else
             write(kou,*)'No such file name: ',trim(text)
          endif
          graphopt%appendfile=' '
          goto 21100
!-----------------------------------------------------------
! PLOT TEXT 
       case(13)
          labelp=>graphopt%firsttextlabel
          if(associated(labelp)) then
             call gparcd('Modify existing text?',cline,last,1,ch1,'NO',q1help)
             if(ch1.eq.'y' .or. ch1.eq.'Y') then
                jp=0
                do while(associated(labelp))
                   jp=jp+1
                   write(kou,2310)jp,labelp%xpos,labelp%ypos,labelp%textline
2310               format(i3,2(1pe12.4),5x,a)
                   labelp=>labelp%nexttextlabel
                enddo
                call gparid('Which text index?',cline,last,kl,1,q1help)
                if(kl.lt.1 .or. kl.gt.jp) then
                   write(*,*)'No such text label'
                   goto 100
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
                if(buperr.ne.0) then
                   write(*,*)'Error reading coordinates'; goto 100
                endif
                labelp%xpos=xxx
                labelp%ypos=xxy
! ask for more options
                goto 21100
             endif
          endif
! input a new label
          call gparrd('X position: ',cline,last,xxx,zero,q1help)
          call gparrd('Y position: ',cline,last,xxy,zero,q1help)
          if(buperr.ne.0) then
             write(*,*)'Error reading coordinates'; goto 100
          endif
          line=' '
          if(noofaxis.eq.2) then
! Calculate the equilibria at the specific point
             write(kou,*)'This is possible only when you plot with'//&
                  ' the same axis as you calculated!'
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
                   line=' '
                endif
! when implemented add the stable phase names to "line" as default for text
             endif
          endif
! There is no gparcd which allows editing the existing text ... emacs!!
          call gparcd('Text: ',cline,last,5,text,line,q1help)
          if(text(1:1).eq.' ') then
             write(*,*)'Label ignored'
             goto 21100
          endif
! I know one should never allocate pointers but this is the only way ???
          allocate(textlabel)
          textlabel%xpos=xxx
          textlabel%ypos=xxy
          textlabel%textline=trim(text)
          if(associated(graphopt%firsttextlabel)) then
             textlabel%nexttextlabel=>graphopt%firsttextlabel
             write(*,*)trim(graphopt%firsttextlabel%textline)
          else
             nullify(textlabel%nexttextlabel)
          endif
          graphopt%firsttextlabel=>textlabel
! the record is now linked from graphopt, nullify the pointer ...
          nullify(textlabel)
! also clean the cline character otherwise labels may be overwritten
          cline=' '
          last=len(cline)
          goto 21100
!-----------------------------------------------------------
! PLOT TIE_LINES increment
       case(14)
          call gparid('Tie-line increment?',cline,last,kl,0,q1help)
          if(kl.lt.0) kl=0
          graphopt%tielines=kl
!          write(*,*)'No implemented yet'
          goto 21100
!-----------------------------------------------------------
! PLOT KEEP does not work ...
       case(15)
          if(btest(graphopt%status,GRKEEP)) then
             graphopt%status=ibclr(graphopt%status,GRKEEP)
!             write(kou,*)'Graphics window closes with mouse'
          else
             graphopt%status=ibset(graphopt%status,GRKEEP)
!             write(kou,*)'Graphics window must be closed separately'
          endif
          write(*,*)'Not implemented yet'
          goto 21100
!-----------------------------------------------------------
! LOGSCALE
       case(16)
          call gparcd('For x or y axis? ',cline,last,1,ch1,'x',q1help)
          if(ch1.eq.'x') then
             if(graphopt%axistype(1).eq.1) then
                write(kou,*)'The x axis set to linear'
                graphopt%axistype(1)=0
             else
                graphopt%axistype(1)=1
             endif
          elseif(ch1.eq.'y') then
             if(graphopt%axistype(2).eq.1) then
                write(kou,*)'The y axis set to linear'
                graphopt%axistype(2)=0
             else
                graphopt%axistype(2)=1
             endif
          else
             write(kou,*)'Please answer x or y'
          endif
          goto 21100
!-----------------------------------------------------------
! unused
       case(17)
          goto 21100
!-----------------------------------------------------------
! unused
       case(18)
          goto 21100
       end SELECT
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
!     call gparcd('Are you sure?',cline,last,1,ch1,'N',q1help)
!     if(ch1.eq.'y' .or. ch1.eq.'Y') then
       stop 'Au revoir'
!=================================================================
! OPTIMIZE and CONTINUE.  Current optimizer is optimizers(optimizer)
    case(24)
       call gparid('Number of iterations: ',cline,last,i1,nopt,q1help)
       if(buperr.ne.0) goto 100
       nopt=i1
!       write(*,606)'dead 1',mexp,nvcoeff,iexit
606    format(a,10i4)
! some optimires have no CONTINUE
       if(optimizer.eq.1) iexit(4)=0
       continue: if(mexp.gt.0 .and. iexit(4).eq.2) then
! iexit(4) from previous optimize allows continue with same Jacobian
          call gparcd('Continue with same Jacobian? ',cline,last,1,&
               ch1,'Y',q1help)
          if(ch1.eq.'Y') then
             ient=1
             goto 987
          endif
       endif continue
! Initiate arrays when new optimization
       ient=0
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
!                   write(*,*)'Set scaled coefficient ',i2,' to ',&
!                        firstash%coeffvalues(i1)
             endif
          enddo
          if(i2.lt.nvcoeff) then
             write(kou,*)'Internal error for variable coefficients',&
                  i2,nvcoeff
             goto 100
          endif
       endif
! caculate how big www is needed, it depends on mexp and nvcoeff
! Size of www from VA05AD 
!             M*N         +(M+N)*N               +N      
       nwc=mexp*nvcoeff+(mexp+nvcoeff)*nvcoeff+nvcoeff
!               +M   +N      +N*N            +N+M+N
       j1=nwc+mexp+nvcoeff+nvcoeff*nvcoeff+2*nvcoeff+mexp
!       allocate(www(j1))
       if(maxw.lt.j1) then
          write(*,*)'Too big problem, increase maxw, current value',maxw
          goto 100
       endif
! JUMP HERE IF CONTINUE optimization
987    continue
! mexp    Number of experiments
! nvcoeff Number of coefficients to be optimized
! errs Array with differences with experiments and calculated values
! coefs Array with coefficients
! VA05AD variables: dstep, dmax2, acc, iterations, output unit, workspace
!                   entry mode, exit mode
       if(mexp.le.0 .or. nvcoeff.le.0) then
          write(kou,569)mexp,nvcoeff
569       format('Cannot optimize with zero experiments or coefficients',2i5)
          goto 100
       endif
       write(*,558)mexp,nvcoeff,maxw
558    format(/'>>>   Start of optimization   >>>'/&
            'Experiments, coefficients and workspace: ',3(1x,i5))
!
       iprint=1
! There is a va05ad emulator called lmdif ...
       call va05ad(mexp,nvcoeff,errs,coefs,dstep,dmax2,acc,nopt,iprint,www,&
            ient,iexit)
!       write(kou,559)iprint,iexit
!559    format(/'Back from optimization',10i3/)
! we must copy the current scaled coefficients back to firstash%coeffvalues
       i2=1
       do i1=0,size(firstash%coeffstate)-1
          if(firstash%coeffstate(i1).ge.10) then
             firstash%coeffvalues(i1)=coefs(i2)
!             write(*,558)i2,i1,firstash%coeffvalues(i1)
!558          format('Saving scaled coefficient ',i3,' to ',i3,1pe12.4)
             i2=i2+1
          endif
       enddo
! calculated errors 
!       write(*,*)'Errors: ',(errs(j1),j1=1,mexp)
!=================================================================
! unused
    CASE(25)
       write(kou,*)'Not implemented yet'
!=================================================================
! unused
    CASE(26)
       write(kou,*)'Not implemented yet'
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
    END SELECT
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
!  implicit double precision (a-h,o-z)
    implicit none
    character cline*(*),ccomnd(*)*(*),query*(*)
    character defansw*16,query1*64,text*256
    integer last,kdef,ncomnd,kom2,lend,lenq
!\end{verbatim}
    logical once
    lenq=len_trim(query)
    if(query(lenq:lenq).eq.'?') then
       query1=query(1:lenq)
    else
       query1=query(1:lenq)//' what?'
       lenq=lenq+6
    endif
    once=.true.
    submenu=0
!  write(kou,10)'submenu 1: ',query(1:lenq),last,cline(last:last+5)
!10  format(a,a,i4,': ',a)
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
       if(eolch(cline,last)) then
! there is no user input passed to this subroutine
          defansw=ccomnd(kdef)
          lend=len_trim(defansw)+1
! note fourth argument 5 copes whole of cline into text
          call gparcd(query1(1:lenq),cline,last,5,text,defansw,q2help)
          if(buperr.ne.0) goto 1000
          cline=text
       else
! skip one character, if next is , take default answer
! write(*,102)'submenu 7: ',last,cline(1:last+5)
!102       format(a,i5,'"',a,'"')
!        last=last+1
          if(cline(last:last).eq.',') then
             submenu=kdef
             goto 1000
          else
             cline=cline(last:)
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
          write(*,*)'Warning, exceeded helprec%level limit'
       endif
    endif
1000 continue
    return
  end function submenu

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ocmon_set_options(option,afo,optionsset)
    implicit none
    character*(*) option
    integer afo
    TYPE(ocoptions) :: optionsset
!\end{verbatim}
    integer next,kom,slen,errno
    character string*64,dummy*128,date*8,time*10
    integer, parameter :: nopt=9
    character (len=16), dimension(nopt) :: copt=&
        ['OUTPUT          ','ALL             ','FORCE           ',&
         'VERBOSE         ','SILENT          ','APPEND          ',&
         '                ','                ','                ']
!
! /? will list options
    afo=0
    if(option(1:2).eq.'? ') then
       write(kou,10)
10     format('Available options (preceeded by /) are:')
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
       case(1) ! /output means ovewrite any previous content
!          write(*,*)'Option not implemented: ',option(1:len_trim(option))
! next argument after = must be a file name
          call getext(option,next,2,string,' ',slen)
! add extention .dat if to extenstion provided
          if(index(string,'.').le.0) then
             string(slen+1:)='.dat '
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
          write(*,*)'output to file: ',string(1:len_trim(string)),optionsset%lut
!-----------------------------------
       case(2) ! /all ??
          write(*,*)'Option not implemented: ',option(1:len_trim(option))
!-----------------------------------
       case(3) ! /force
          write(*,*)'Option not implemented: ',option(1:len_trim(option))
!-----------------------------------
       case(4) ! /verbose
          globaldata%status=ibset(globaldata%status,GSVERBOSE)
          write(kou,*)'VERBOSE option set'
!-----------------------------------
       case(5) ! /silent
          globaldata%status=ibclr(globaldata%status,GSVERBOSE)
          globaldata%status=ibset(globaldata%status,GSSILENT)
!-----------------------------------
       case(6) ! /APPEND, open file and write at end
!          write(*,*)'Option not implemented: ',option(1:len_trim(option))
! next argument after = must be a file name
          call getext(option,next,2,string,' ',slen)
! add extention .dat if to extenstion provided
          if(index(string,'.').le.0) then
             string(slen+1:)='.dat '
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
! write not allowed after fininfg EOF, we must backspace
220       continue
          backspace(21)
! write a header
          call date_and_time(date,time)
          write(21,232)'appended: ',date(1:4),date(5:6),date(7:8),&
               time(1:2),time(3:4)
          write(kou,231)string(1:len_trim(string))
231       format('Output will be appended to file: ',a)
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
       write(*,*)'ouput unit reset to screen',optionsset%lut
    endif
!1000 continue
    return
  end subroutine ocmon_reset_options


!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine listoptshort(lut,mexp,errs)
! short listing of optimizing variables and result
    integer lut,mexp
    double precision errs(*)
!    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: neweq
    integer i1,i2,j1,j2,j3
    character name1*24,line*80
    double precision xxx
!
    write(lut,610)
610 format(/'List of coefficents with non-zero values'/&
         'Name  Current value   Start value    Scaling factor',&
         ' RSD')
    name1=' '
    do i1=0,size(firstash%coeffstate)-1
!       write(*,611)i1,firstash%coeffstate(i1)
!611                format('Coefficient ',i2,' state ',i2)
       coeffstate: if(firstash%coeffstate(i1).ge.10) then
! optimized variable, read from TP constant array
          call get_value_of_constant_index(firstash%coeffindex(i1),xxx)
          call makeoptvname(name1,i1)
          write(lut,615)name1(1:3),xxx,&
               firstash%coeffstart(i1),firstash%coeffscale(i1),zero
615       format(a,2x,4(1pe15.6))
          if(firstash%coeffstate(i1).eq.11) then
! there is a prescribed minimum
             write(lut,616)' minimum ',firstash%coeffmin(i1)
616          format(6x,'Prescribed ',a,': ',1pe12.4)
          elseif(firstash%coeffstate(i1).eq.12) then
! there is a prescribed minimum
             write(lut,616)' maximum ',firstash%coeffmax(i1)
          elseif(firstash%coeffstate(i1).eq.13) then
! there is a prescribed minimum
             write(lut,617)firstash%coeffmin(i1),firstash%coeffmax(i1)
617          format(6x,'Prescribed min and max: ',2(1pe12.4))
          elseif(firstash%coeffstate(i1).gt.13) then
             write(lut,*)'Wrong coefficent state, set to 10'
             firstash%coeffstate(i2)=10
          endif
       elseif(firstash%coeffstate(i1).gt.0) then
! fix variable status
          call get_value_of_constant_index(firstash%coeffindex(i1),xxx)
          call makeoptvname(name1,i1)
          write(lut,615)name1(1:3),xxx
       elseif(firstash%coeffscale(i1).ne.0) then
! coefficient with negative status, status set to 1
          call get_value_of_constant_index(firstash%coeffindex(i1),xxx)
          write(lut,619)i1,firstash%coeffscale(i1),xxx,zero
619       format('Wrong state for coefficient ',i3,4(1pe12.4))
          firstash%coeffstate(i1)=1
       endif coeffstate
    enddo
! list all experiments, only possible after first optimize
    if(mexp.eq.0) then
       write(lut,666)
666    format(/'No optimization so no results'/)
       goto 1000
    endif
    write(lut,620)size(firstash%eqlista),mexp
620 format(/'List of ',i5,' equilibria with ',i5,&
         ' experimental data values'/&
         'Equil name    Weight Experiment $ calculated',24x,&
         'Error')
    j3=0
    experim: do i1=1,size(firstash%eqlista)
! skip equilibria with zero weight
       neweq=>firstash%eqlista(i1)%p1
       if(neweq%weight.eq.zero) cycle experim
       name1=neweq%eqname(1:12)
       if(associated(neweq%lastexperiment)) then
          i2=neweq%lastexperiment%seqz
          do j2=1,i2
             j1=1
             line=' '
! this subroutine returns experiment and calculated value: "H=1000:200 $ 5000"
             call meq_get_one_experiment(j1,line,j2,neweq)
             j3=j3+1
             write(lut,622)name1(1:12),neweq%weight,line(1:45),&
                  errs(j3)
622          format(a,2x,F5.2,2x,a,1x,1pe12.4)
! list the equilibrium name just for the first (or only) experiment
             name1=' '
          enddo
       endif
    enddo experim
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

