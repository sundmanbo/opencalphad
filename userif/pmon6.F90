!
MODULE cmon1oc
!
! Copyright 2012-2014, Bo Sundman, France
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
! first version of a command monitor for OC
!************************************
!
! currently the step/map/plot is included in liboceq
  use ocsmp
!  use liboceq
!
  implicit none
!
! option record
  TYPE ocoptions
     integer lut
  end TYPE ocoptions
  type(ocoptions) :: optionsset
!
contains
!
  subroutine oc_command_monitor(linkdate)
! command monitor
    implicit none
!
! various symbols and texts
    character symbol*24,name1*24,name2*24,name3*24,line*80,model*72
    integer, parameter :: ocmonversion=20
! element symbol and array of element symbols
    character elsym*2,ellist(20)*2
! more texts for various purposes
    character text*72,string*256,chc*3,phtype*1,ch1*1,defansw*16
    character axplot(3)*24,axplotdef(3)*24
    character xquest*24,form*32,longstring*2048
! separate file names for remembering and providing a default
    character ocmfile*64,ocufile*64,tdbfile*64,ocdfile*64,texfile*64
! prefix and suffix for composition sets
    character prefix*4,suffix*4
! element mass
    double precision mass
! constituent indices in a phase
    integer, dimension(maxconst) :: knr
! constituent fractions of a phase
    double precision, dimension(maxconst) :: yarr
! stoichiometry of a specis and sublattice sites of a phase
    double precision, dimension(maxsubl) :: stoik,sites
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
! axis data structures
    type(map_axis), dimension(5) :: axarr
! if more than one start equilibrium these are linked using the ceq%next index
    type(gtp_equilibrium_data), pointer :: starteq
! for map results
    type(map_node), pointer :: maptop,mapnode,maptopsave
    type(map_line) :: mapline
    integer noofaxis,noofstarteq
! this should be removed
!    TYPE(ssm_node), pointer :: resultlist
!<<<<<<<--------------------------------------------------------------
! used for element data and R*T
    double precision h298,s298,rgast
! temporary reals
    double precision xxx,xxy
! array for constituent in endmember and interaction, parameter property type
    integer endm(10),lint(2,3),typty
! input data for grid minimizer
    double precision, dimension(maxel) :: xknown,aphl
! arrays for grid minimization results
    integer, dimension(maxel) :: iphl,icsl,nyphl
! selected kommand and subcommands
    integer kom,kom1,kom2,kom3,kom4
! selected output mode for results and the default
    integer listresopt,lrodef
! integers used for elements, phases, composition sets, equilibria, defaults
    integer iel,iph,ics,ieq,idef
!-------------------
! temporary integer variables in loops etc
    integer jl,i1,i2,iax
! more temporary integers
    integer jp,kl,svss,language,last
! and more temporary integers
    integer ll,lokcs,lokfil,lokph,lokres,loksp,lrot,maxax
! and more temporary integers
    integer minimizer,mode,ndl,neqdef,noelx,nofc,nopl,nops,nsl,nv,nystat
! temporary matrix
    double precision latpos(3,3)
!-------------------
! loop variable when entering constituents of a phase
    integer icon
! array with constituents in sublattices when entering a phase
    character, dimension(maxconst) :: const*24
! passed as date when program was linked
    character linkdate*(*)
! for macro and logfile and repeating questions
    logical logok,stop_on_error,once,wildcard,twice
! unit for logfile input, 0 means no logfile
    integer logfil
! remember default for calculate phase
    integer defcp,defnsl
! for state variables as conditions
    integer istv,indices(4,10),iref,unit
    double precision coeffs(10)
    TYPE(gtp_state_variable), pointer :: stvr
    TYPE(gtp_state_variable), dimension(10) :: stvarr
    TYPE(gtp_condition), pointer :: pcond,firstc
! current equilibrium records
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq
    TYPE(gtp_phase_varres), pointer :: parres
!
    character (len=34) :: quest1='Number of sites on sublattice xx: '
!
    character actual_arg(2)*16
    character cline*128,option*80,aline*128,plotfile*64,eqname*24
!----------------------------------------------------------------
! here are all commands and subcommands
!    character (len=64), dimension(6) :: oplist
    integer, parameter :: ncbas=27,nclist=15,ncalc=9,ncent=15,ncread=6
    integer, parameter :: ncam1=12,ncset=18,ncadv=3,ncstat=6,ncdebug=6
    integer, parameter :: nselect=6,nlform=6
    integer, parameter :: ncamph=12,nclph=6,nccph=3,nrej=6,nsetph=6
    integer, parameter :: nsetphbits=15,ncsave=6,nplt=6,nstepop=6
! basic commands
    character (len=16), dimension(ncbas), parameter :: cbas=&
       ['AMEND           ','CALCULATE       ','SET             ',&
        'ENTER           ','EXIT            ','LIST            ',&
        'QUIT            ','READ            ','SAVE            ',&
        'HELP            ','INFORMATION     ','BACK            ',&
        'NEW             ','MACRO           ','ABOUT           ',&
        'DEBUG           ','SELECT          ','DELETE          ',&
        'STEP            ','MAP             ','PLOT            ',&
        'HPCALC          ','FIN             ','                ',&
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
         'STATE_VARIABLES ','BIBLIOGRAPHY    ','PARAMETER_IDENTI',&
         'AXIS            ','TPFUN_SYMBOLS   ','QUIT            ',&
         '                ','EQUILIBRIA      ','RESULTS         ',&
         'CONDITIONS      ','SYMBOLS         ','LINE_EQUILIBRIA ']
!-------------------
! subsubcommands to LIST DATA
    character (len=16), dimension(nlform) :: llform=&
        ['SCREEN          ','TDB             ','MACRO           ',&
         'LATEX           ','                ','                ']
!-------------------
! subsubcommands to LIST PHASE
    character (len=16), dimension(nclph) :: clph=&
        ['DATA            ','CONSTITUTION    ','MODEL           ',&
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
         ['ONLY_G          ','G_AND_DGDY      ','ALL_DERIVATIVES ']
!-------------------
! subcommands to ENTER
    character (len=16), dimension(ncent) :: center=&
         ['TPFUN_SYMBOL    ','ELEMENT         ','SPECIES         ',&
         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
         'CONSTITUTION    ','EXPERIMENT      ','QUIT            ',&
         'EQUILIBRIUM     ','SYMBOL          ','OPTIMIZE_COEFF  ',&
         'COPY_OF_EQUILIB ','                ','                ']
!-------------------
! subcommands to READ
    character (len=16), dimension(ncread) :: cread=&
         ['UNFORMATTED     ','TDB             ','QUIT            ',&
         'DIRECT          ','                ','                ']
!-------------------
! subcommands to SAVE
    character (len=16), dimension(ncsave) :: csave=&
         ['UNFORMATTED     ','TDB             ','MACRO           ',&
         'DIRECT          ','LATEX           ','QUIT            ']
!-------------------
! subcommands to AMEND first level
! many of these should be subcommands to PHASE
    character (len=16), dimension(ncam1) :: cam1=&
         ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
         'PHASE           ','PARAMETER       ','BIBLIOGRAPHY    ',&
         'TPFUN_SYMBOL    ','CONSTITUTION    ','QUIT            ',&
         'COMPONENTS      ','GENERAL         ','DEBYE_MODEL     ']
!-------------------
! subsubcommands to AMEND PHASE
    character (len=16), dimension(ncamph) :: camph=&
         ['MAGNETIC_CONTRIB','COMPOSITION_SET ','DISORDERED_FRACS',&
         'GLAS_TRANSITION ','QUIT            ','DEFAULT_CONSTIT ',&
         'DEBYE_CP_MODEL  ','EINSTEIN_CP_MDL ','INDEN_WEI_MAGMOD',&
         'ELASTIC_MODEL_A ','                ','                ']
!-------------------
! subcommands to SET.  
    character (len=16), dimension(ncset) :: cset=&
         ['CONDITION       ','STATUS          ','ADVANCED        ',&
         'LEVEL           ','INTERACTIVE     ','REFERENCE_STATE ',&
         'QUIT            ','ECHO            ','PHASE           ',&
         'UNITS           ','LOG_FILE        ','WEIGHT          ',&
         'NUMERIC_OPTIONS ','AXIS            ','INPUT_AMOUNTS   ',&
         'VERBOSE         ','AS_START_EQUILIB','                ']
! subsubcommands to SET STATUS
    character (len=16), dimension(ncstat) :: cstatus=&
         ['ELEMENT         ','SPECIES         ','PHASE           ',&
         'CONSTITUENT     ','                ','                ']
!        123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET ADVANCED
    character (len=16), dimension(ncadv) :: cadv=&
         ['EQUILIB_TRANSF  ','QUIT            ','                ']
!         123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET PHASE
    character (len=16), dimension(nsetph) :: csetph=&
         ['QUIT            ','STATUS          ','DEFAULT_CONSTITU',&
          'AMOUNT          ','BITS            ','                ']
!         123456789.123456---123456789.123456---123456789.123456
!-------------------
! subsubsubcommands to SET PHASE BITS
    character (len=16), dimension(nsetphbits) :: csetphbits=&
         ['FCC_PERMUTATIONS','BCC_PERMUTATIONS','IONIC_LIQUID_MDL',&
         'AQUEOUS_MODEL   ','QUASICHEMICAL   ','FCC_CVM_TETRADRN',&
         'FACT_QUASICHEMCL','NO_AUTO_COMP_SET','                ',&
         '                ','                ','                ',&
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
         '                ','                ','                ']
!-------------------
! subcommands to SELECT, maybe some should be CUSTOMMIZE ??
    character (len=16), dimension(nselect) :: cselect=&
         ['EQUILIBRIUM     ','MINIMIZER       ','GRAPHICS        ',&
         'LANGUAGE        ','                ','                ']
!-------------------
! subcommands to DELETE
    character (len=16), dimension(nrej) :: crej=&
         ['ELEMENTS        ','SPECIES         ','PHASE           ',&
          'QUIT            ','COMPOSITION_SET ','EQUILIBRIUM     ']
!-------------------
! subcommands to PLOT OPTIONS
    character (len=16), dimension(nplt) :: cplot=&
         ['PLOT            ','XRANGE          ','YRANGE          ',&
         'XTEXT           ','YTEXT           ','TITLE           ']
!-------------------
!        123456789.123456---123456789.123456---123456789.123456
! minimizers
    character (len=16), dimension(2) :: minimizers=&
         ['LUKAS_HILLERT   ','SUNDMAN_HILLERT ']
!------------------------------------------------------------------------
!
! some defaults
    language=1
    logfil=0
    defcp=1
!
    write(kou,10)linkdate(1:len_trim(linkdate)),ocmonversion,&
         gtpversion,hmsversion,smpversion
10  format(/' Open Calphad (OC) software linked ',a,&
         ' with command line monitor ',i2//&
         ' This program is available with a GNU General Public License.'/&
         ' It includes the General Thermodynamic Package, version ',A/&
         " and Hillert's equilibrium calculation algorithm version ",A/&
         ' and step/map/plot software version ',A/)
!
! jump here after NEW to reinitiallize all local variables also
20  continue
! clear file names
    ocmfile=' '; ocufile=' '; tdbfile=' '
! plot ranges and their defaults
    graphopt%rangedefaults=0
    graphopt%labeldefaults=0
    graphopt%plotmin=zero
    graphopt%dfltmin=zero
    graphopt%plotmax=one
    graphopt%dfltmax=one
! initiate on-line help
    call init_help('ochelp.hlp ')
! set default minimizer
    minimizer=2
! by default no stop on error and no nogfile
    stop_on_error=.false.
    lokfil=0
!
! in init_gtp the first equilibrium record is created and 
! firsteq has been set to that
!
25  continue
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
    do jl=1,3
       axplotdef(jl)=' '
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
! >>> we should remove all equilibria !!
! here one should read a user initialisation file as a macro
! file can be at current directory or at home directory
!
!============================================================
! return here for next command
100 continue
    if(gx%bmperr.ne.0) goto 990
    if(buperr.ne.0) goto 990
! turn off options set
    call ocmon_reset_options(optionsset)
! initiate command level for help routines
    call helplevel1('Initiate help level for OC')
! read the command line with gparc to have output on logfile
    last=len(aline)
    aline=' '
    cline=' '
    call gparc('OC2: ',aline,last,5,cline,' ',tophlp)
    jl=1
! with empty line just prompt again
    if(eolch(cline,jl)) goto 100
! with macro command character just prompt again
    if(cline(jl:jl).eq.'@') goto 100
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
             call ocmon_set_options(option)
!             nops=nops+1
!             oplist(nops)=option
!             write(kou,*)'found option: ',oplist(nops)
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
!        'COMPONENTS      ','GENERAL         ','                ']
    CASE(1)
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
! this cannot be cleared here, there may be other
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
!        'GLAS_TRANSITION ','QUIT            ','DEFAULT_CONSTIT ',&
!        'DEBYE_CP_MODEL  ','EINSTEIN_CP_MDL ','INDEN_WEI_MAGNET',&
!        'ELASTIC_MODEL_A ','                ','                ']
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
                  cline,last,jl,idef,q1help)
             if(buperr.ne.0) goto 990
             call get_phase_record(iph,lokph)
             call add_magrec_inden(lokph,1,jl)
!....................................................
          case(2) ! amend phase <name> composition set add/remove
             call gparcd('Add new set? ',cline,last,1,ch1,'Y ',q1help)
             if(buperr.ne.0) goto 990
             if(ch1.eq.'Y' .or. ch1.eq.'y') then
                call gparc('Prefix: ',cline,last,1,prefix,' ',q1help)
                call gparc('Suffix: ',cline,last,1,suffix,' ',q1help)
!                call add_composition_set(iph,prefix,suffix,ics,ceq)
                call add_composition_set(iph,prefix,suffix,ics)
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
! ch1 is parameter suffix, jl=0 means never completely disordred (sigma)
             call gparcd('Can the phase be totally disordered? ',cline,last,&
                  1,ch1,'N',q1help)
             if(buperr.ne.0) goto 990
             if(ch1.eq.'N') then
! like sigma
                jl=0
             else
! like FCC ordering
                jl=1
             endif
             ch1='D'
             call add_fraction_set(iph,ch1,ndl,jl)
             if(gx%bmperr.ne.0) goto 990
!....................................................
          case(4) ! amend phase <name> glas_transition
             call add_addrecord(iph,glastransmodela)
!....................................................
          case(5) ! amend phase quit
             goto 100
!....................................................
          case(6) ! amend phase <name> default_constitution
             call ask_default_constitution(cline,last,iph,ics,ceq)
!....................................................
          case(7) ! amend phase <name> Debye model
             call add_addrecord(iph,debyecp)
!....................................................
          case(8) ! amend phase einstein cp model
             call add_addrecord(iph,einsteincp)
!....................................................
          case(9) ! amend phase wei_inden_magnetic_model
             call add_addrecord(iph,weimagnetic)
!....................................................
          case(10) ! amend phase elastic model
             call add_addrecord(iph,elasticmodela)
!....................................................
          case(11) ! amend phase ??
             write(kou,*)'Not implemented yet'
!....................................................
          case(12) ! amend phase ??
             write(kou,*)'Not implemented yet'
          END SELECT
!-------------------------
       case(5) ! amend parameter
          write(kou,*)'Not implemented yet'
!-------------------------
       case(6) ! amend bibliography
          call enter_reference_interactivly(cline,last,1,jl)
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
             call enter_tpfun_interactivly(cline,last,string,jp)
! this stores the tpfun, lrot<0 means the symbol already exists
             lrot=-1
             call enter_tpfun(name1,string,lrot,.FALSE.)
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
          if(associated(ceq%lastcondition)) then
             write(kou,*)'Warning, your conditions may not be valid after this'
          endif
          call amend_components(cline,last,ceq)
!-------------------------
       case(11) ! amend general
          call amend_global_data(cline,last)
!-------------------------
       case(12) ! Nothing defined
          continue
       END SELECT
!=================================================================
! calculate subcommands
    CASE(2)
       kom2=submenu(cbas(kom),cline,last,ccalc,ncalc,8)
       SELECT CASE(kom2)
       CASE DEFAULT
          write(kou,*)'calculate subcommand error'
          goto 100
!-------------------------
       CASE(1) ! calculate TPFUN symbols , use current values of T and P
          write(kou,2011)notpf(),ceq%tpval
2011      format(/'Calculating ',i3,' functions for T,P=',F10.2,1PE15.7/&
               3x,'No   F',11x,'F.T',9x,'F.P',9x,'F.T.T',7x,'F.T.P',7x,'F.P.P')
          call cpu_time(starting)
          do jl=1,notpf()
             call eval_tpfun(jl,ceq%tpval,val,ceq%eq_tpres)
             if(gx%bmperr.gt.0) goto 990
             write(kou,2012)jl,val
2012         format(I5,1x,6(1PE12.4))
          enddo
          call cpu_time(ending)
          write(kou,2013)ending-starting
2013      format('CPU time used: ',1pe15.6)
!---------------------------------------------------------------
       case(2) ! calculate phase, _all _only_g or _g_and_dgdy, separated later
          ! asks for phase name and constitution
          call ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
          if(gx%bmperr.ne.0) goto 990
! if iph<0 then * has been given as phase name
          if(iph.lt.0) then
             write(kou,*)'Cannot loop for all phases'
             goto 100
          endif
! use current value of T and P
          rgast=globaldata%rgas*ceq%tpval(1)
!
          kom3=submenu('Calculate what for phase?',cline,last,ccph,nccph,defcp)
!        if(kom2.le.0) goto 100
!       ph-a ph-G ph-G+dg/dy
          defcp=kom3
          SELECT CASE(kom3)
!.......................................................
          CASE DEFAULT
             write(kou,*)'Calculate phase subcommand error'
!.......................................................
          case(1) ! calculate phase < > only G
             call calcg(iph,ics,0,lokres,ceq)
             if(gx%bmperr.ne.0) goto 990
             parres=>ceq%phase_varres(lokres)
             write(kou,2031)(rgast*parres%gval(jl,1),jl=1,4)
             write(kou,2032)parres%gval(1,1)/parres%abnorm(1),parres%abnorm(1)
2031         format('G/N, dG/dT:',4(1PE16.8))
2032         format('G/N/RT, N:',2(1PE16.8))
!.......................................................
          case(2) ! calculate phase < >  G and dG/dy
             call calcg(iph,ics,1,lokres,ceq)
             if(gx%bmperr.ne.0) goto 990
             parres=>ceq%phase_varres(lokres)
             nofc=noconst(iph,ics,firsteq)
             write(kou,2031)(rgast*parres%gval(jl,1),jl=1,4)
             write(kou,2041)(rgast*parres%dgval(1,jl,1),jl=1,nofc)
2041         format('dG/dy:   ',4(1PE16.8))
!.......................................................
          case(3) ! calculate phase < > all
             call tabder(iph,ics,ceq)
             if(gx%bmperr.ne.0) goto 990
          END SELECT
!----------------------------------
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
!----------------------------------
       case(4) ! calculate transition
! set a phase fix and remove one condition.  One must have calculated an
! equilibrium
          write(kou,2090)
2090      format('To calculate when a phase will appear/dissapear',&
               ' by releasing a condition.')
          if(btest(ceq%status,EQNOEQCAL)) then
             write(kou,2095)
2095         format('You must make an equilibrium calculation before using',&
                  ' this command.')
             goto 100
          endif
          call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
          call find_phase_by_name(name1,iph,ics)
          if(gx%bmperr.ne.0) goto 990
          jl=test_phase_status(iph,ics,xxx,ceq)
          if(jl.eq.PHFIXED) then
             write(kou,*)'Phase status already fixed'
             goto 100
          endif
          call list_conditions(kou,ceq)
          write(kou,2097)
2097      format('You must release one condition, give its number')
          call gparid('Condition number',cline,last,jl,1,q1help)
          if(jl.le.0 .or. jl.gt.noel()+2) then
             write(kou,*)'No such condition'
             goto 100
          endif
! this finds condition with given number
          call locate_condition(jl,pcond,ceq)
          if(gx%bmperr.ne.0) goto 990
          if(pcond%active.eq.0) then
! the condition is active, deactivate it!
             pcond%active=1
          else
             write(kou,*)'This condition is not active!'
             goto 100
          endif
! Condition released, now set the phase as fix with zero moles
          call change_phase_status(iph,ics,PHFIXED,xxx,ceq)
          if(gx%bmperr.ne.0) goto 990
! Calculate equilibrium
          call calceq2(1,ceq)
          if(gx%bmperr.ne.0) goto 990
! get the value of the released condition and set it to the new value
          stvr=>pcond%statvar(1)
          call state_variable_val(stvr,xxx,ceq)
          if(gx%bmperr.ne.0) goto 990
          write(kou,2099)xxx
2099      format('The transition occurs at ',1pe16.8/'Now set as condition')
          pcond%prescribed=xxx
          pcond%active=0
! set phase back as entered and stable
          write(*,*)'Set phase back as entered'
          call change_phase_status(iph,ics,PHENTSTAB,zero,ceq)
          if(gx%bmperr.ne.0) goto 990
!
!          write(*,*)'Not implemented yet'
          goto 100
!----------------------------------
       case(5) ! quit
          goto 100
!-----------------------------------------------------------
       case(6) ! calculate global grid minimum
! extract values for mass balance calculation from conditions
          call extract_massbalcond(ceq%tpval,xknown,xxx,ceq)
          if(gx%bmperr.ne.0) goto 990
! debug output
!  write(*,2101)xxx,(xknown(i),i=1,noel())
2101      format('N&x: ',F6.3,9F8.5)
! generate grid and find the phases and constitutions for the minimum.
          call global_gridmin(1,ceq%tpval,xknown,nv,iphl,icsl,aphl,nyphl,&
               yarr,cmu,ceq)
          if(gx%bmperr.ne.0) goto 990
          write(kou,2102)'Phases: ',nv,(iphl(jl),icsl(jl),jl=1,nv)
2102      format(a,i2,11(i4,i2))
!---------------------------------------------------------------
       case(7) ! calculate symbol
!          call evaluate_all_svfun(kou,ceq)
! to calculate derivatives this must be in the minimizer module
          call gparcd('Name ',cline,last,1,name1,'*',q1help)
          if(name1(1:1).eq.'*') then
             call meq_evaluate_all_svfun(kou,ceq)
          else
             call capson(name1)
             call find_svfun(name1,istv,ceq)
             if(gx%bmperr.ne.0) goto 990
             mode=1
             actual_arg=' '
             xxx=meq_evaluate_svfun(istv,actual_arg,mode,ceq)
             write(*,2047)name1(1:len_trim(name1)),xxx
2047         format(a,'= ',1pe16.8)
          endif
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(8) ! calculate equilibrium for current equilibrium ceq
          if(minimizer.eq.1) then
! Lukas minimizer, first argiment=1 means use grid minimizer
!           call calceq1(1,ceq)
             write(kou,*)'Not implemented yet'
          else
             call calceq2(1,ceq)
          endif
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(9) ! calculate all equilibria
          write(kou,*)'Not implemented yet'
       END SELECT
!=================================================================
    CASE(3) ! set subcommands
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
             call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
             call find_phase_by_name(name1,iph,ics)
             if(gx%bmperr.ne.0) then
                if(name1(1:2).eq.'* ') then
                   iph=-1
                   gx%bmperr=0
                else
                   goto 990
                endif
             else
                jl=get_phase_status(iph,ics,text,i1,xxx,ceq)
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
! no longer available if(ch1.eq.'N') nystat=5
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
                jl=get_phase_status(iph,ics,text,i1,xxy,ceq)
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
          name1='advanced command'
          kom3=submenu(name1,cline,last,cadv,ncadv,1)
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
! set bit that data bay be inconsistent
             eqlista(i1)%status=ibset(eqlista(i1)%status,EQINCON)
!.................................................................
          case(2) ! quit
             continue
!.................................................................
          case(3) ! nothing yet
             write(*,*)'Not implemented yet'
          end select
!-----------------------------------------------------------
       case(4) ! set LEVEL
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! end of macro excution (can be nested)
       case(5) ! set INTERACTIVE
          call macend(cline,last,logok)  
!-----------------------------------------------------------
       case(6) ! set REFERENCE_STATE
!        write(kou,*)'Reference states not implemented yet'
          call gparc('Component name: ',cline,last,1,name1,' ',q1help)
          call find_component_by_name(name1,iel,ceq)
          if(gx%bmperr.ne.0) goto 100
          call gparc('Reference phase: ',cline,last,1,name1,' ',q1help)
          if(name1(1:4).eq.'SER ') then
! this means no reference phase, SER is at 298.15K and 1 bar
             iph=-1
          else
             call find_phase_by_name(name1,iph,ics)
             if(gx%bmperr.ne.0) goto 100
! temperature * means always to use current temperature
             call gparr('Temperature: /*/: ',cline,last,xxx,zero,q1help)
             if(xxx.le.zero) then
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
             write(kou,3104)
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
             jl=1
          else
             jl=0
          endif
          call set_echo(jl)
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
! copied from 3045
          case(2) ! SET PHASE STATUS <phase> <status>
             if(iph.gt.0) then
                jl=get_phase_status(iph,ics,text,i1,xxx,ceq)
                if(gx%bmperr.ne.0) goto 100
                if(xxx.ge.zero) then
                   write(kou,3046)text(1:i1),xxx
                else
                   write(kou,3047)text(1:i1)
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
                jl=get_phase_status(iph,ics,text,i1,xxy,ceq)
                if(gx%bmperr.ne.0) goto 100
                if(xxy.ge.zero) then
                   write(kou,3048)text(1:i1),xxy
                else
                   write(kou,3049)text(1:i1)
                endif
             else
                write(kou,*)'New status set for all phases'
             endif
! end copied from 3045
!............................................................
          case(3:4) !set phase default constitution wildcard allod, also AMOUNT
             if(kom3.eq.3) then
! set phase default constituntion
                call set_default_constitution(iph,ics,ceq)
             else
! set phase amount
                call gparrd('Amount: ',cline,last,xxx,zero,q1help)
                call set_phase_amounts(iph,ics,xxx,ceq)
             endif
!............................................................
! subsubsub command
          case(5) ! set phase bits
             if(iph.lt.0) then
                write(kou,*)'Wildcard not allowed'
                goto 100
             endif
             call get_phase_record(iph,lokph)
             kom4=submenu('Set which bit?',cline,last,csetphbits,nsetphbits,8)
             SELECT CASE(kom4)
             CASE DEFAULT
                write(kou,*)'Set phase bit subcommand error'
                goto 100
!............................................................
             case(1) ! FCC_PERMUTATIONS FORD
! if check returns .true. it is not allowed to set FORD or BORD
                if(check_minimal_ford(lokph)) goto 100
                call set_phase_status_bit(lokph,PHFORD)
             case(2) ! BCC_PERMUTATIONS BORD
                if(check_minimal_ford(lokph)) goto 100
                call set_phase_status_bit(lokph,PHBORD)
             case(3) ! IONIC_LIQUID_MDL 
                call set_phase_status_bit(lokph,PHIONLIQ)
             case(4) ! AQUEOUS_MODEL   
                call set_phase_status_bit(lokph,PHAQ1)
             case(5) ! QUASICHEMICAL   
                call set_phase_status_bit(lokph,PHQCE)
             case(6) ! FCC_CVM_TETRADRN
                call set_phase_status_bit(lokph,PHCVMCE)
             case(7) ! FACT_QUASICHEMCL
                call set_phase_status_bit(lokph,PHFACTCE)
             case(8) ! NO_AUTO_COMP_SET, not allod to create compsets automatic
                call set_phase_status_bit(lokph,PHNOCS)
             case(9) ! ELASTIC_MODEL_A
! set by amend??
                continue
             end SELECT
!............................................................
          case(6)
             write(kou,*)'Not implemented yet'
          END SELECT
!-------------------------------------------------------------
       case(10) ! set UNIT
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
          write(kou,*)'Not implemented yet'
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
! replacing an existing axis, get the current as defualt
!             call encode(...)
! deallocate axarr ?? redundant ??
!             deallocate(axarr(iax)%indices)
!             deallocate(axarr(iax)%coeffs)
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
             noofaxis=noofaxis-1
             write(kou,*)'One axis removed'
             goto 100
          else ! add or change axis variable
             i1=len_trim(text)
             if(text(i1:i1).eq.':') then
! condition given as an index in the condition list terminated by :
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
!                call decode_state_variable2(text,istv,indices,iref,unit,&
!                     stvr,ceq)
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
! default step 1/40
!          dinc=axvalold(3,iax)
!          dinc=0.025*(axval(2,iax)-axval(1,iax))
          dinc=0.025*(axarr(iax)%axmax-axarr(iax)%axmin)
          call gparrd('Increment:',cline,last,xxx,dinc,q1help)
          if(buperr.ne.0) goto 100
          axarr(iax)%axinc=xxx
! iax can be smaller than noofaxis if an existing axis has been changed
          if(iax.gt.noofaxis) noofaxis=iax
!  write(*,3602)(axval(i,iax),i=1,3)
3602      format(/'axlimits: ',3(1pe12.4))
!-------------------------------------------------------------
       case(15) ! set input amounts
          call set_input_amounts(cline,last,ceq)
!-------------------------
       case(16) ! SET VERBOSE
! This sets permanent verbose for all commands.  If on turn it off
          write(kou,*)'Not implemented yet'
!-------------------------
! the current set of condition sill be used as start equilibrium for map/step
! Calculate the equilibrium and ask for a direction.
       case(17) ! SET AS_START_EQUILIBRIUM
          if(noofaxis.lt.2) then
             write(kou,*)'You must set two axis forst'
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
       case(18)
          write(kou,*)'Not implemented yet'
       END SELECT
!=================================================================
! enter subcommand for element, species etc
    CASE(4)
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
                call enter_tpfun_interactivly(cline,last,string,jp)
                if(gx%bmperr.ne.0) goto 990
! here the function is stored
                lrot=0
                call enter_tpfun(name1,string,lrot,.FALSE.)
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
          call new_element(elsym,name1,name2,mass,h298,s298)
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
          call new_species(name1,noelx,ellist,stoik)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(4) ! enter phase
          if(.not.allowenter(2)) then
             gx%bmperr=4125
             goto 990
          endif
          call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
! ionic liquid require special sorting of constituents on anion sublattice
          call capson(name1)
          defnsl=1
          if(name1(1:4).eq.'GAS ') then
             phtype='G'
             model='CEF-RKM'
          elseif(name1(1:7).eq.'LIQUID ') then
             phtype='L'
             model='CEF-RKM'
          elseif(name1(1:9).eq.'IONIC_LIQ') then
             phtype='L'
             model='IONIC_LIQUID'
             defnsl=2
          else
             phtype='S'
             model='CEF-RKM'
          endif
          call gparid('Number of sublattices: ',cline,last,nsl,defnsl,q1help)
          if(buperr.ne.0) goto 990
          if(nsl.le.0) then
             write(kou,*)'At least one configurational space!!!'
             goto 100
          elseif(nsl.gt.10) then
             write(kou,*)'Maximum 10 sublattices'
             goto 100
          endif
          icon=0
          sloop: do ll=1,nsl
! 'Number of sites on sublattice xx: '
!  123456789.123456789.123456789.123
             once=.true.
4042         continue
             write(quest1(31:32),4043)ll
4043         format(i2)
             call gparrd(quest1,cline,last,sites(ll),one,q1help)
             if(buperr.ne.0) goto 990
             if(sites(ll).le.1.0D-2) then
                write(kou,*)'Number of sites must be larger than 0.01'
                if(once) then
                   once=.false.
                   goto 4042
                else
                   goto 100
                endif
             endif
! This should be extended to allow several lines of input
! 4 means up to ;
             once=.true.
4045         continue
             call gparc('All Constituents: ',cline,last,4,text,';',q1help)
             if(buperr.ne.0) goto 990
             knr(ll)=0
             jp=1
4047         continue
             if(eolch(text,jp)) goto 4049
             if(model(1:13).eq.'IONIC_LIQUID ' .and. ll.eq.1 &
                  .and. knr(1).eq.0) then
! a very special case: a single "*" is allowed on 1st sublattice for ionic liq
                if(text(jp:jp).eq.'*') then
                   icon=icon+1
                   const(icon)='*'
                   knr(1)=1
                   cycle sloop
                endif
             endif
             call getname(text,jp,name3,1,ch1)
             if(buperr.eq.0) then
                icon=icon+1
                const(icon)=name3
                knr(ll)=knr(ll)+1
! increment jp to bypass a separating , 
                jp=jp+1
                goto 4047
             elseif(once) then
                write(kou,*)'Input error ',buperr,', at ',jp,', please reenter'
                buperr=0; once=.false.; goto 4045
             else
                goto 100
             endif
             buperr=0
4049         continue
          enddo sloop
          call new_phase(name1,nsl,knr,const,sites,model,phtype)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(5) ! enter parameter is always allowed
          call enter_parameter_interactivly(cline,last)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(6) ! enter bibliography
          call enter_reference_interactivly(cline,last,0,jl)
          if(gx%bmperr.ne.0) goto 990
          write(kou,*)'Bibliography number is ',jl
!---------------------------------------------------------------
       case(7) ! enter constitution
          call ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(8) ! enter experiment
          write(kou,*)'Not implemented yet'
!---------------------------------------------------------------
       case(9)  ! enter QUIT
          goto 100
!---------------------------------------------------------------
       case(10) ! enter equilibrium is always allowed if there are phases
          if(.not.allowenter(3)) then
             write(kou,*)'There must be at least one phase'
             goto 100
          endif
          call gparc('Name: ',cline,last,1,text,' ',q1help)
          if(buperr.ne.0) goto 100
          call enter_equilibrium(text,ieq)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
       case(11) ! enter symbol (for state variables expressions)
          call enter_svfun(cline,last,ceq)
          if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
! enter optimizing coefficients
       case(12)
          call enter_optvars
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
! enter not used
       case(14)
!---------------------------------------------------------------
! enter not used
       case(15)
          goto 100
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
! list subcommands
    CASE(6)
       kom2=submenu(cbas(kom),cline,last,clist,nclist,12)
       if(kom2.le.0) goto 100
       SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'LIST FORMAT subcommand error'
          goto 100
!-----------------------------------------------------------
       case(1) ! list data for everything, not dependent on equilibrium!!
! NOTE output file set by /output=
          kom3=submenu('Output format?',cline,last,llform,nlform,1)
          if(kom.gt.0) then
             call list_many_formats(kom3,kou)
          else
             write(kou,*)'Unknown format'
          endif
!-----------------------------------------------------------
       case(2) ! list short with status bits
          write(kou,6022)ceq%eqname,globaldata%rgasuser,&
               globaldata%pnorm,globaldata%status
6022      format('Equilibrium name',9x,'Gas constant Pressure norm',&
               22x,'Status'/1x,a,1pe12.4,2x,1pe12.4,20x,z8)
          call list_all_elements(kou)
          call list_all_species(kou)
          call list_all_phases(kou,ceq)
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
             call list_phase_data(iph,kou)
!...............................................................
! list phase constitution
          case(2) ! list phase constitution
             idef=110
             call gparid('Output mode: ',cline,last,mode,idef,q1help)
             if(buperr.ne.0) goto 990
!  call list_phase_results(iph,ics,mode,kou,firsteq)
             write(kou,6051)ceq%eqno,ceq%eqname
6051         format('Output for equilibrium: ',i3,', ',a)
             call list_phase_results(iph,ics,mode,kou,ceq)
             if(gx%bmperr.ne.0) goto 990
!...............................................................
          case(3) ! list phase model (including disordered fractions)
             write(kou,6070)'For ',ceq%eqno,ceq%eqname
6070      format(a,'equilibrium: ',i3,', ',a)
             call list_phase_model(iph,ics,kou,ceq)
          END SELECT
!------------------------------
       case(4)  ! list state variable or parameter identifier value, loop.
6099      continue
          if(btest(ceq%status,EQNOEQCAL) .or. btest(ceq%status,EQFAIL)) then
             write(kou,6101)
6101         format(' *** Warning,',&
                'equilibrium not calculated, values are probably wrong')
          elseif(btest(ceq%status,EQINCON)) then
             write(kou,6102)
6102         format(' *** Warning, values can be inconsistent with',&
                ' current conditions')
          endif
6105      continue
! NOTE: 4th argument is 5 because otherwise , will terminate reading cline
! and state variables like x(fcc,cr) will not work.
          if(kom2.eq.4) then
             call gparc('State variable: ',cline,last,5,line,' ',q1help)
          else
             write(kou,*)'Remember always to give the phase!'
             call gparc('Parameter ident: ',cline,last,5,line,' ',q1help)
          endif
          jl=1
          if(eolch(line,jl)) goto 100
          model=' '
          if(index(line,'*').gt.0) then
! generate many values
! i1 values are resturned in yarr with dimension maxconst. "model" not used
             call get_many_svar(line,yarr,maxconst,i1,model,ceq)
             if(gx%bmperr.eq.0) then
                write(kou,6107)(yarr(i2),i2=1,i1)
6107            format('Values: ',5(1pe14.6)/(8x,5(1pe14.6)))
             endif
          else
! in model the state variable is returned as generated by the program
             call get_state_var_value(line,xxx,model,ceq)
             if(gx%bmperr.eq.0) then
                write(kou,6108)model(1:len_trim(model)),xxx
6108            format(1x,a,'=',1PE15.7)
             endif
          endif
          if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
             write(kou,*)bmperrmess(gx%bmperr)
          elseif(gx%bmperr.ne.0) then
             write(kou,*)'Error code ',gx%bmperr
          endif
          gx%bmperr=0
          goto 6105
!-----------------------------------------------------------
       case(5) ! list data bibliography
          call list_bibliography(kou)
!-----------------------------------------------------------
       case(6) ! list parameter symbols
          call list_defined_properties(kou)
!-----------------------------------------------------------
       case(7) ! list axis
          if(noofaxis.le.0) then
             write(kou,*)'No axis set'
             goto 100
          endif
          write(kou,6131)
6131      format(4x,'Axis variable',12x,'Min',9x,'Max',9x,'Max increment')
!6131      format(4x,'Axis variable',12x,'Start',7x,'Final',7x,'Increment')
          do iax=1,noofaxis
             jp=1
             call get_one_condition(jp,text,axarr(iax)%seqz,ceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Condition sequential index: ',iax,axarr(iax)%seqz
                goto 990
             endif
! we just want the expression, remove the value including the = sign
             jp=index(text,'=')
             text(jp:)=' '
!             write(kou,6132)iax,axvar(iax),(axval(jl,iax),jl=1,3)
             write(kou,6132)iax,text(1:24),&
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
                call list_tpfun(lrot,0,longstring)
                call wrice2(kou,0,12,78,1,longstring)
                if(iel.gt.1) goto 6140
             endif
          else
             call list_all_funs(kou)
          endif
!------------------------------------------------------------
       case(9) ! list quit
!------------------------------------------------------------
       case(10) ! unused list subcommand
        continue
!-----------------------------------------------------------
       case(11) ! list equilibria (not result)
          do iel=1,noeq()
             if(associated(ceq,eqlista(iel))) then
                write(kou,6202)iel,eqlista(iel)%eqname
6202            format(i5,2x,'**',2x,a)
             else
                write(kou,6203)iel,eqlista(iel)%eqname
6203            format(i5,6x,a)
             endif
          enddo
!------------------------------
       case(12) ! list results
          call gparid('Output mode: ',cline,last,listresopt,lrodef,q1help)
          lrodef=listresopt
          write(kou,6051)ceq%eqno,ceq%eqname
!  if(btest(globaldata%status,GSEQFAIL)) then
          if(btest(ceq%status,EQFAIL)) then
             write(kou,6305)
6305         format(/' *** The results listed are not a valid equilibrium',&
                  ' as last calculation failed'/)
          elseif(btest(globaldata%status,GSNOPHASE)) then
             write(kou,*)'No results as no data'
             goto 100
!  elseif(btest(globaldata%status,GSNOEQCAL)) then
          elseif(btest(ceq%status,EQNOEQCAL)) then
             write(kou,6307)
6307         format(/' *** The results listed does not represent',&
                  ' a calculated equilibrium'/)
          elseif(btest(ceq%status,EQINCON)) then
             write(kou,6306)
6306         format(/' *** The results listed may be inconsistent',&
                  ' with the current conditions'/)
          endif
          write(kou,6302)'Conditions .........'
6302      format(a,40('.'),':')
6303      format(/a,40('.'),':')
          call list_conditions(kou,ceq)
          write(kou,6303)'Some global data ...'
          call list_global_results(kou,ceq)
          write(kou,6303)'Some component data '
          jl=1
          if(listresopt.ge.4 .and. listresopt.le.7) then
             jl=2
          endif
          call list_components_result(kou,jl,ceq)
          write(kou,6303)'Some Phase data ....'
! mode >1000 lists stable phases only
          if(listresopt.le.1) then
! stable phases with mole fractions
             mode=1000
          elseif(listresopt.eq.2) then
! stable phases with mole fractions and constitution
             mode=1010
          elseif(listresopt.eq.3) then
! stable phases with mole fractions in alphabetical order
             mode=1100
          elseif(listresopt.eq.4) then
! stable phases with mass fractions
             mode=1001
          elseif(listresopt.eq.5) then
! stable phases with mass fractions in alphabetical order
             mode=1101
          elseif(listresopt.eq.6) then
! stable phases with mass fractions and constitution
             mode=1011
          elseif(listresopt.eq.7) then
! all phases with mass fractions
             mode=1
          elseif(listresopt.eq.8) then
! all phases with mole fractions and constitution
             mode=110
          elseif(listresopt.eq.9) then
! all phases with mole fractions and constitution in alphabetical order
             mode=10
          else
! all phase with with mole fractions
             mode=0
          endif
          ics=1
          do iph=1,noph()
             ics=0
6310         continue
             ics=ics+1
             call list_phase_results(iph,ics,mode,kou,ceq)
             if(gx%bmperr.ne.0) then
! if error take next phase
                gx%bmperr=0
             else
! else take next composition set
                goto 6310
             endif
          enddo
!  write(kou,*)
          call list_phases_with_positive_dgm(mode,kou,ceq)
!  if(btest(globaldata%status,GSEQFAIL)) then
          if(btest(ceq%status,EQFAIL)) then
             write(kou,6305)
!  elseif(btest(globaldata%status,GSNOEQCAL)) then
          elseif(btest(ceq%status,EQNOEQCAL)) then
             write(kou,6307)
!  elseif(btest(globaldata%status,GSINCON)) then
          elseif(btest(ceq%status,EQINCON)) then
             write(kou,6306)
          endif
!------------------------------
       case(13) ! list conditions
          write(kou,6070)'Conditions for ',ceq%eqno,ceq%eqname
          call list_conditions(kou,ceq)
!------------------------------
       case(14) ! list symbols (state variable functions, not TP funs)
          call list_all_svfun(kou,ceq)
!------------------------------
! list lines, output of calculated and stored equilibria
       case(15)
! temporary listing of all stored equilibria as test
          call list_stored_equilibria(kou,axarr,maptop)
       end SELECT
!=================================================================
! quit
    case(7)
       if(cline(1:1).eq.'q') then
          call gparcd('Are you sure?',cline,last,1,ch1,'N',q1help)
       else
! upper case Q will quit
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
    case(8)
       if(noel().ne.0) then
          write(kou,*)'You already have data, read destroys your current data'
          call gparcd('Do you want to continue? ',cline,last,1,ch1,'N',q1help)
          if(.not.(ch1.eq.'y' .or. ch1.eq.'Y')) then
             write(kou,*)'Command ignored, to read answer yes' 
             goto 100
          else
! all records must be removed and init_gtp is called.  This is fragile ...
             call new_gtp
             if(gx%bmperr.ne.0) goto 990
             write(kou,*)'All previous data deleted'
          endif
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
!---------------------------------------------------------
       case(2) ! read TDB
          if(tdbfile(1:1).ne.' ') then
             text=tdbfile
             call gparcd('File name: ',cline,last,1,tdbfile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,tdbfile,' ',q1help)
          endif
          write(kou,8202)
8202      format('Give the elements to select, finish with empty line')
          jp=1
8210      continue
          call gparc('Select elements/all/: ',&
               cline,last,1,ellist(jp),' ',q1help)
          if(ellist(jp).ne.'  ') then
             call capson(ellist(jp))
             jp=jp+1
             if(jp.gt.size(ellist)) then
                write(kou,*)'Max number of elements selected: ',size(ellist)
             else
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
!-----------------------------------------------------------
!8300      continue
       case(3) ! read quit
          goto 100
!-----------------------------------------------------------
       case(4) ! read direct
          write(*,*)'Not implemented yet'
!-----------------------------------------------------------
       case(5) ! read ??
          goto 100
!-----------------------------------------------------------
       case(6) ! read ??
          goto 100
       end SELECT
!=================================================================
! save unformatted
    case(9)
       kom2=submenu(cbas(kom),cline,last,csave,ncsave,1)
       if(kom2.le.0 .or. kom2.gt.ncsave) goto 100
!
       call gparc('Comment line: ',cline,last,5,model,' ',q1help)
       SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'save subcommand error'
!-----------------------------------------------------------
       case(1) ! save unformatted
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
             write(kou,*)'Missing file name'
             goto 100
          endif
          if(jp.gt.0) ocufile(jp+1:)='.ocu '
          text='U '//model
          call gtpsave(ocufile,text)
!-----------------------------------------------------------
       case(2) ! save TDB
          if(tdbfile(1:1).ne.' ') then
             text=tdbfile
             call gparcd('File name: ',cline,last,1,tdbfile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,tdbfile,' ',q1help)
          endif
          jp=0
          kl=index(tdbfile,'.')
          if(kl.le.0) then
             jp=len_trim(tdbfile)
          elseif(tdbfile(kl+1:kl+1).eq.' ') then
! just ending a filename with . not accepted as extention
             jp=kl
          endif
          if(kl.le.0 .and. jp.le.0) then
             write(kou,*)'Missing file name'
             goto 100
          endif
          if(jp.gt.0) tdbfile(jp+1:)='.TDB '
          text='T '//model
          call gtpsave(tdbfile,text)
!-----------------------------------------------------------
       case(3) ! save MACRO
          if(ocmfile(1:1).ne.' ') then
             text=ocmfile
             call gparcd('File name: ',cline,last,1,ocmfile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,ocmfile,' ',q1help)
          endif
          jp=0
          kl=index(ocmfile,'.')
          if(kl.le.0) then
             jp=len_trim(ocmfile)
          elseif(ocmfile(kl+1:kl+1).eq.' ') then
! just ending a filename with . not accepted as extention
             jp=kl
          endif
          if(kl.le.0 .and. jp.le.0) then
             write(kou,*)'Missing file name'
             goto 100
          endif
          if(jp.gt.0) ocmfile(jp+1:)='.OCM '
          text='M '//model
          call gtpsave(ocmfile,text)
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
             write(kou,*)'Missing file name'
             goto 100
          endif
          if(jp.gt.0) ocdfile(jp+1:)='.ocd '
          text='M '//model
          call gtpsave(ocdfile,text)
!-----------------------------------------------------------
       case(5) ! save LaTeX format for publishing
          if(texfile(1:1).ne.' ') then
             text=texfile
             call gparcd('File name: ',cline,last,1,texfile,text,q1help)
          else
             call gparc('File name: ',cline,last,1,texfile,' ',q1help)
          endif
          jp=0
          kl=index(texfile,'.')
          if(kl.le.0) then
             jp=len_trim(texfile)
          elseif(texfile(kl+1:kl+1).eq.' ') then
! just ending a filename with . not accepted as extention
             jp=kl
          endif
          if(kl.le.0 .and. jp.le.0) then
             write(kou,*)'Missing file name'
             goto 100
          endif
          if(jp.gt.0) texfile(jp+1:)='.tex '
          text='M '//model
          call gtpsave(texfile,text)
!-----------------------------------------------------------
       case(6) ! save quit
          continue
       end SELECT
!=================================================================
! help ... just list commands
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
!----- deallocate local axis records
       do jp=1,noofaxis
          deallocate(axarr(jp)%axcond)
          deallocate(axarr(jp)%indices)
          deallocate(axarr(jp)%coeffs)
       enddo
       noofaxis=0
       if(associated(maptop)) then
          deallocate(maptop)
       endif
       nullify(starteq)
       graphopt%rangedefaults=0
!
! this routine fragile, inside new_gtp init_gtp is called
       call new_gtp
       write(kou,*)'All data removed'
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
15010  format(/'This is Open Calphad (OC), a free software for ',&
            'thermodynamic calculations,'/&
            'available for download at http://www.opencalphad.com'//&
            'This software is protected by the GNU General Public License'/&
            'You may freely distribute copies as long as you also provide ',&
            'the source code.'/'The software is provided "as is" without ',&
            'any warranty of any kind, either'/'expressed or implied.'/&
            'The full license text is provided with the software or can be ',&
            'obtained from'/'the Free Software Foundation ',&
            'http://www.fsf.org'//&
            'Copyright 2010-2014, several persons.'/&
            'Contact person Bo Sundman, bo.sundman@gmail.com'/&
            'This version linked ',a/)
       goto 100
!=================================================================
! debug subcommands
    case(16)
!       write(*,*)'Calculating equilibrium record size'
       kom2=ceqsize(ceq)
       write(kou,*)'Equilibrium record size: ',kom2
       kom2=submenu(cbas(kom),cline,last,cdebug,ncdebug,0)
       SELECT CASE (kom2)
!------------------------------
       CASE DEFAULT
          write(kou,*)'Default case ',kom2
!------------------------------
! debug free lists
       CASE(1)
! list all tuples
          do jp=1,nooftuples
             write(kou,16020)jp,phasetuple(jp)%phase,phasetuple(jp)%compset
16020        format(i3,': ',2i4)
          enddo
          call list_free_lists(kou)
!------------------------------
! debug stop_on_error
       CASE(2)
          stop_on_error=.true.
!------------------------------
! debug elasticity
       CASE(3)
          write(kou,*)'Input current lattice parameter values (3x3 matrix)',&
               ' for phase 1'
          iph=1
          ics=1
          xxx=7.1D-6
          xxy=1.0D-12
          call gparrd('lattice par (1,1):',cline,last,latpos(1,1),xxx,nohelp)
          call gparrd('lattice par (1,2):',cline,last,latpos(1,2),xxy,nohelp)
          call gparrd('lattice par (1,3):',cline,last,latpos(1,3),xxy,nohelp)
          call gparrd('lattice par (2,1):',cline,last,latpos(2,1),xxy,nohelp)
          call gparrd('lattice par (2,2):',cline,last,latpos(2,2),xxx,nohelp)
          call gparrd('lattice par (2,3):',cline,last,latpos(2,3),xxy,nohelp)
          call gparrd('lattice par (3,1):',cline,last,latpos(3,1),xxy,nohelp)
          call gparrd('lattice par (3,2):',cline,last,latpos(3,2),xxy,nohelp)
          call gparrd('lattice par (3,3):',cline,last,latpos(3,3),xxx,nohelp)
          call set_lattice_parameters(iph,ics,latpos,ceq)
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
          else
! check if number
             jl=1
             call getint(text,jl,i1)
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
          call gparcd('Give name of minimizer?',&
               cline,last,1,text,'LUKAS ',q1help)
          call capson(text)
          jl=len_trim(text)
!  write(*,*)'input: ',text(1:16),jl
          if((text(1:jl).eq.minimizers(1)(1:jl))) then
             minimizer=1
          else
             minimizer=2
          endif
          write(kou,*)'Selected minimizer: ',minimizers(minimizer)
!-----------------------------------------------------------
       case(3) ! select graphics
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
       case(4) ! select languegae, at present only 1 English and 2 French
          write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
       case(5)
          goto 100
!-----------------------------------------------------------
       case(6)
          goto 100
       END SELECT
!=================================================================
! delete
    case(18)
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
! delete equilibrium
       case(6)
          call gparcd('Equilibrium name: ',cline,last,1,name1,'_MAP* ',q1help)
          if(buperr.ne.0) goto 990
          call delete_equilibria(name1,ceq)
          if(gx%bmperr.ne.0) goto 990
       end SELECT
!=================================================================
! step
    case(19)
       if(noofaxis.ne.1) then
          write(kou,*)'You must set exactly one independent axis variable',&
               ' for a step calculation.'
          goto 100
       endif
! check if adding results
       if(associated(maptop)) then
          write(kou,*)'There are some results already from step or map'
          call gparcd('Delete them?',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
! there should be a more careful deallocation to free memory
             deallocate(maptop%saveceq)
             nullify(maptop)
             nullify(maptopsave)
          else
             write(kou,*)'Results removed'
             maptopsave=>maptop
             nullify(maptop)
          endif
! this should preferably be done directly after map/step, but kept for debug
          call delete_equilibria('_MAP*',ceq)
       endif
       kom2=submenu('Options?',cline,last,cstepop,nstepop,1)
       SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'No such plot option'
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
          call map_setup(maptop,noofaxis,axarr,starteq)
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
          call step_separate(maptop,noofaxis,axarr,starteq)
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
! map
    case(20)
       if(noofaxis.lt.2) then
          write(kou,*)'You must set two axis with independent variables'
          goto 100
       endif
       write(kou,20014)
20014   format('The map command is fragile, please send problematic diagrams',&
            ' to the',/'OC development team'/)
       if(associated(maptop)) then
          write(kou,*)'There are some results already form step or map'
          call gparcd('Reinitiate?',cline,last,1,ch1,'Y',q1help)
          if(ch1.eq.'y' .or. ch1.eq.'Y') then
             deallocate(maptop%saveceq)
             nullify(maptop)
             nullify(maptopsave)
! this should preferably be done directly after map/step, but kept for debug
             call delete_equilibria('_MAP*',ceq)
             if(gx%bmperr.ne.0) then
                write(kou,*)'Error removing old MAP equilibria'
                goto 990
             endif
          else
             maptopsave=>maptop
             nullify(maptop)
          endif
! this should preferably be done directly after map/step, but kept for debug
          call delete_equilibria('_MAP*',ceq)
       endif
! maptop is returned as main map/step record for results
! noofaxis is current number of axis, axarr is array with axis data
! starteq is start equilibria, if empty set it to ceq
       if(.not.associated(starteq)) then
          starteq=>ceq
          starteq%next=0
       endif
       call map_setup(maptop,noofaxis,axarr,starteq)
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
! mark that interactive listing of conditions and results may be inconsistent
       ceq%status=ibset(ceq%status,EQINCON)
       if(gx%bmperr.ne.0) goto 990
!=================================================================
! plot
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
!2100      continue
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
             if(wildcard) then
                write(*,*)'Wildcards allowed for one axis only'
                goto 21000
             else
                wildcard=.TRUE.
             endif
          endif
          if(axplotdef(iax).ne.axplot(iax)) then
! if not same axis variable remove any ranges defined !!!
             graphopt%rangedefaults(iax)=0
          endif
! remember axis as default
          axplotdef(iax)=axplot(iax)
       enddo
       call gparcd('GNUPLOT file name',cline,last,1,plotfile,'ocgnu',q1help)
       form=' '
! first argument is the number of plot axis, always 2 at present
       jp=2
       if(associated(maptopsave)) then
          write(*,*)'We link to maptopsave'
          maptop%plotlink=>maptopsave
!       else
!          write(*,*)'There is no maptopsave'
       endif
!-----------------------------------------------------------
! plot options subcommand, default is PLOT, NONE does not work ...
21100   continue
! give an empty line as a sublte alert for plot options
       write(kou,*)
       kom2=submenu('Options?',cline,last,cplot,nplt,1)
       SELECT CASE(kom2)
!-----------------------------------------------------------
       CASE DEFAULT
          write(kou,*)'No such plot option'
!-----------------------------------------------------------
! no more options, just PLOT!
       case(1)
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
! PLOT LABEL
       case(6)
          call gparcd('Plot title',cline,last,5,line,'DEFAULT',q1help)
          if(line(1:8).eq.'DEFAULT ') then
             graphopt%labeldefaults(1)=0
          else
             graphopt%plotlabels(1)=line
             graphopt%labeldefaults(1)=len_trim(graphopt%plotlabels(1))
          endif
!-----------------------------------------------------------
       end SELECT
!-----------------------------------------------------------
2190   continue
       call ocplot2(jp,axplot,plotfile,maptop,axarr,graphopt,form)
       if(gx%bmperr.ne.0) goto 990
       call gparc('Hardcopy (P for postscript)?',&
            cline,last,1,form,'none',q1help)
       if(form.ne.'none') then
          call capson(form)
          call ocplot2(jp,axplot,plotfile,maptop,axarr,graphopt,form)
       endif
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
! no command defined yet
    case(24)
       write(kou,*)'Not implemented yet'
!=================================================================
! no command defined yet
    case(25)
       write(kou,*)'Not implemented yet'
!=================================================================
! unused
    case(26)
       write(kou,*)'Not implemented yet'
!=================================================================
! unused
    case(27)
       write(kou,*)'Not implemented yet'
!=================================================================
!
    END SELECT
! command executed, prompt for another command unless error code
    if(gx%bmperr.eq.0) goto 100
!============================================================
! handling errors
990 continue
    write(kou,*)
    write(kou,*)'Error codes: ',gx%bmperr,buperr
    if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
       write(kou,*)bmperrmess(gx%bmperr)
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
10  format(a,a,i4,': ',a)
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
102       format(a,i5,'"',a,'"')
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
  subroutine ocmon_set_options(option)
    implicit none
    character*(*) option
    TYPE(ocoptions) :: optionsset
!    TYPE(ocoptions), pointer :: optionsset
!\end{verbatim}
    integer next,kom,slen
    character string*64
    integer, parameter :: nopt=6
    character (len=16), dimension(nopt) :: copt=&
        ['OUTPUT          ','ALL             ','FORCE           ',&
         'VERBOSE         ','SILENT          ','                ']
!
    kom=ncomp(option,copt,nopt,next)
    if(kom.le.0) then
       write(kou,*)'Unknown option ignored: ',option(1:len_trim(option))
       goto 1000
    else
       select case(kom)
       case default
          write(*,*)'Option not implemented: ',option(1:len_trim(option))
!-----------------------------------
       case(1) ! /output means APPEND
!          write(*,*)'Option not implemented: ',option(1:len_trim(option))
! next argument after = must be a file name
          call getext(option,next,2,string,' ',slen)
          write(*,*)'should append to file: ',string(1:len_trim(string))
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
       case(6) ! /
          write(*,*)'Option not implemented: ',option(1:len_trim(option))
!-----------------------------------
       end select
    endif
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
! The only thing tested here is verbose
    if(btest(globaldata%status,GSVERBOSE)) then
       if(.not.btest(globaldata%status,GSSETVERB)) then
! if user has SET VERBOSE do not resest VERBOSE
          globaldata%status=ibclr(globaldata%status,GSVERBOSE)
       endif
    endif
1000 continue
    return
  end subroutine ocmon_reset_options

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

END MODULE cmon1oc

! ====================================================================
! Feature:
! Parameters that depend on a constituent like mobilities and the new
! Bohr magneton number.  Previously only average Bohr magneton numbers
! have been used.
! The Bohr magneton number of Fe in FCC can be a function of composition i.e.
! BMAGN&FE(BCC) = x_Fe*BMAGN&FE(BCC,FE)+x_CR*BMAGN&FE(BCC,CR)+
!                 X_FE*X_CR*BMAGN(BCC,CR,FE)+ ...
! A similar function BMAGN&CR(BCC)= ....
! The identifications "&FE" or "&CR" are stored in the property record.
!
! Each property has a unique number but usually only a few are needed
! When calculating for a phase the there is an internal propertytype counter
! incremented for each property found.  This keeps track of the actual
! property that has been calculated for each value of the property counter.
! Property type 1 is always the Gibbs energy.
!
! ===========================================================
! Parallellizing
! Fraction values may be different for the same phase if there are
! parallell calculation of equilibria.  
! Note that the number of sites for a sublattice are normally static 
! but for ionic liquids they are dynamic.
! Function values must also be stored separate in each thread
! as they are calculated for different values of T and P.
! The error code must be separate in each thread as it may
! not be fatal but can be handelled by the process
! 
