!
MODULE oc_cmon1
!
! Copyright 2012-13, Bo Sundman, France
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
!  use matsmin (used in smp)
  use smp
!  use lukasmin (with the calceq1 routine)
!
  implicit none
!
contains
!
  subroutine oc_command_monitor(linkdate)
! command monitor
    implicit none
!
! various symbols and texts
  character symbol*24,name1*24,name2*24,name3*24,line*80,model*72
! element symbol and array of element symbols
  character elsym*2,ellist(10)*2
! more texts for various purposes
  character text*72,string*256,chc*3,filename*32,phtype*1,ch1*1,defansw*16
  character xquest*24,form*32
! axis variables and limits
  character axvar(5)*24
  double precision axval(3,5),axvalold(3,5),tpa(2)
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
! estimated chemical potentials after a grid minimization
  double precision cmu(maxel)
! cpu time measurements
  double precision ending,starting
! used for axis variables
  double precision dinc,dmin,dmax
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
  logical logok,stop_on_error,once
! unit for logfile input, 0 means no logfile
  integer logfil
! remember default for calculate phase
  integer defcp
! current equilibrium records
  TYPE(gtp_equilibrium_data), pointer :: ceq
  TYPE(gtp_phase_varres), pointer :: parres
  TYPE(ssm_node), pointer :: resultlist
!
  character (len=34) :: quest1='Number of sites on sublattice xx: '
!
  character cline*128,option*64,aline*128,plotfile*64
  character (len=64), dimension(6) :: oplist
  integer, parameter :: ncbas=24,nclist=15,ncopt=6,ncalc=9,ncent=12,ncread=3
  integer, parameter :: ncam1=12,ncset=18,ncadv=3,ncstat=6,ncdebug=6,nselect=6
  integer, parameter :: ncamph=12,nclph=6,nccph=3,nrej=6,nsetph=6
  integer, parameter :: nsetphbits=15,ncsave=3
! basic commands
  character (len=16), dimension(ncbas) :: cbas=&
       ['AMEND           ','CALCULATE       ','SET             ',&
        'ENTER           ','EXIT            ','LIST            ',&
        'QUIT            ','READ            ','SAVE            ',&
        'HELP            ','INFORMATION     ','BACK            ',&
        'NEW             ','MACRO           ','ABOUT           ',&
        'DEBUG           ','SELECT          ','DELETE          ',&
        'STEP            ','MAP             ','PLOT            ',&
        'HPCALC          ','FIN             ','                ']
! in French
!        'MODIFIEZ        ','CALCULEZ        ','REGLEZ          ',&
!        'ENTREZ          ','EXIT            ','AFFICHER        ',&
!        'QUIT            ','LIRE            ','SAUVGARDE       ',&
!        'AIDEZ           ','INFORMATION     ','RETURNEZ        ',&
!        'NOUVEAU         ','MACRO           ','ABOUT           ',&
!        'DEBUG           ','SELECTIONEZ     ','EFFACEZ         ',&
!        'STEP            ','MAP             ','DESSINEZ        ',&
!        'HPCALC          ','FIN             ','                ']
! options preceeded by -  
! for example "list -out=myfile.dat all_data" or
! "list all_data -out=myfile.dat"
  character (len=16), dimension(ncopt) :: copt=&
       ['OUTPUT          ','ALL             ','FORCE           ',&
        'VERBOSE         ','SILENT          ','                ']
!-------------------
! subcommands to LIST
  character (len=16), dimension(nclist) :: clist=&
       ['DATA            ','SHORT           ','PHASE           ',&
        'STATE_VARIABLES ','REFERENCES      ','PARAMETER_IDENTI',&
        'AXIS            ','TPFUN_SYMBOLS   ','QUIT            ',&
        'VALUE_OF_PARA_ID','EQUILIBRIA      ','RESULTS         ',&
        'CONDITIONS      ','SYMBOLS         ','                ']
!-------------------
! subsubcommands to LIST PHASE
  character (len=16), dimension(nclph) :: clph=&
       ['DATA            ','CONSTITUTION    ','MODEL           ',&
        '                ','                ','                ']
!-------------------
! subcommands to CALCULATE
  character (len=16), dimension(ncalc) :: ccalc=&
       ['TPFUN_SYMBOLS   ','PHASE           ','NO_GLOBAL       ',&
        '                ','QUIT            ','GLOBAL_GRIDMIN  ',&
        'SYMBOL          ','EQUILIBRIUM     ','ALL_EQUILIBRIA  ']
!-------------------
! subcommands to CALCULATE PHASE
  character (len=16), dimension(nccph) :: ccph=&
       ['ONLY_G          ','G_AND_DGDY      ','ALL_DERIVATIVES ']
!-------------------
! subcommands to ENTER
  character (len=16), dimension(ncent) :: center=&
       ['TPFUN_SYMBOL    ','ELEMENT         ','SPECIES         ',&
        'PHASE           ','PARAMETER       ','REFERENCE       ',&
        'CONSTITUTION    ','EXPERIMENT      ','QUIT            ',&
        'EQUILIBRIUM     ','SYMBOL          ','                ']
!-------------------
! subcommands to READ
  character (len=16), dimension(ncread) :: cread=&
       ['UNFORMATTED     ','TDB             ','QUIT            ']
!-------------------
! subcommands to SAVE
  character (len=16), dimension(ncsave) :: csave=&
       ['UNFORMATTED     ','TDB             ','MACRO           ']
!-------------------
! subcommands to AMEND first level
! many of these should be subcommands to PHASE
  character (len=16), dimension(ncam1) :: cam1=&
       ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
        'PHASE           ','PARAMETER       ','REFERENCE       ',&
        'TPFUN_SYMBOL    ','CONSTITUTION    ','QUIT            ',&
        'COMPONENTS      ','GENERAL         ','DEBYE_MODEL     ']
!-------------------
! subsubcommands to AMEND PHASE
  character (len=16), dimension(ncamph) :: camph=&
       ['MAGNETIC_CONTRIB','COMPOSITION_SET ','DISORDERED_FRACS',&
        'GLAS_TRANSITION ','QUIT            ','DEFAULT_CONSTIT ',&
        'DEBYE_CP_MODEL  ','EINSTEIN_CP_MDL ','INDEN_WEI_MAGNET',&
        'ELASTIC_MODEL_A ','                ','                ']
!-------------------
! subcommands to SET.  
  character (len=16), dimension(ncset) :: cset=&
       ['CONDITION       ','STATUS          ','ADVANCED        ',&
        'LEVEL           ','INTERACTIVE     ','REFERENCE_STATE ',&
        'QUIT            ','ECHO            ','PHASE           ',&
        'UNITS           ','LOG_FILE        ','WEIGHT          ',&
        'NUMERIC_OPTIONS ','AXIS            ','INPUT_AMOUNTS   ',&
        '                ','                ','                ']
! subsubcommands to SET STATUS
  character (len=16), dimension(ncstat) :: cstatus=&
       ['ELEMENT         ','SPECIES         ','PHASE           ',&
        'CONSTITUENT     ','                ','                ']
!        123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET ADVANCED
  character (len=16), dimension(ncadv) :: cadv=&
       ['LEVEL           ','                ','                ']
!        123456789.123456---123456789.123456---123456789.123456
! subsubcommands to SET PHASE
  character (len=16), dimension(nsetph) :: csetph=&
       ['CONSTITUTION    ','STATUS          ','DEFAULT_CONSTITU',&
        'AMOUNT          ','BITS            ','                ']
!        123456789.123456---123456789.123456---123456789.123456
!-------------------
! subsubsubcommands to SET PHASE BITS
  character (len=16), dimension(nsetphbits) :: csetphbits=&
       ['FCC_PERMUTATIONS','BCC_PERMUTATIONS','IONIC_LIQUID_MDL',&
        'AQUEOUS_MODEL   ','QUASICHEMICAL   ','FCC_CVM_TETRADRN',&
        'FACT_QUASICHEMCL','NO_AUTO_COMP_SET','                ',&
        '                ','                ','                ',&
        '                ','                ','                ']
!        123456789.123456---123456789.123456---123456789.123456
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
        'QUIT            ','                ','                ']
!-------------------
!        123456789.123456---123456789.123456---123456789.123456
! minimizers
  character (len=16), dimension(2) :: minimizers=&
       ['LUKAS_HILLERT   ','SUNDMAN_HILLERT ']
!
! some defaults
  language=1
  logfil=0
  defcp=1
!
  write(kou,10)10,linkdate(1:len_trim(linkdate)),gtpversion,hmsversion
10 format(/' Command monitor ',i2,' of the Open Calphad (OC) software',&
        ' linked ',a/&
        ' This program is available with a GNU General Public License.'/&
        ' It includes the General Thermodynamic Package, version ',A/&
        ' and the minimizer ',A/)
20 continue
!  call init_gtp
  if(gx%bmperr.ne.0) then
     write(kou,*)'Initiation failed'
     if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
        write(kou,22)gx%bmperr,bmperrmess(gx%bmperr)
22      format(' *** error: ',i7/a)
     else
        write(kou,*)'Unknown error code: ',gx%bmperr
     endif
     stop 'Have a nice day anyway'
  endif
! initiate on-line help
  call init_help('ochelp.hlp ')
! set default minimizer
  minimizer=2
!
  stop_on_error=.false.
! no nogfile
  lokfil=0
!
! in init_gtp the first equilibrium record is created and 
! firsteq has been set to that
!
! jump here after NEW to reinitiallize all local variables also
25 continue
! default values of T and P.  NOTE these are not set as conditions
  firsteq%tpval(1)=1.0D3
  firsteq%tpval(2)=1.0D5
!
! default list result option
  lrodef=1
! default axis values
  maxax=0
  do jl=1,5
     axvalold(1,jl)=3.0D2
     axvalold(2,jl)=3.0D3
     axvalold(3,jl)=3.0D1
  enddo
! remove any results from step and map
  nullify(resultlist)
! set default fractions when entering composition
  xknown=one
! set default equilibrium to 1 and current equilibrium (CEQ) to firsteq
  neqdef=1
  ceq=>firsteq
! here one should read a user initialisation file as a macro
! file can be at current directory or at home directory
!
!============================================================
! return here for a new command
100 continue
  if(gx%bmperr.ne.0) goto 990
  if(buperr.ne.0) goto 990
! initiate command level for help routines
  call helplevel1('Initiate help level for OC')
! read the command line with gparc to have output on logfile
  last=len(aline)
  aline=' '
  cline=' '
  call gparc('OC1: ',aline,last,5,cline,' ',tophlp)
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
!     call nghelp(cline,last,cbas,ncbas)
     call q3help(cline,last,cbas,ncbas)
     goto 100
  else
! check for options .... this does not work yet
! and one should check for options after each subcommand or value entered
     nops=0
110  continue
     if(.not.eolch(cline,last)) then
        if(cline(last:last).eq.'-') then
           call getext(cline,last,2,option,' ',nopl)
           if(buperr.ne.0) then
              write(kou,*)'Error reading option'
              buperr=0; goto 100
           endif
           nops=nops+1
           oplist(nops)=option
           write(kou,*)'option: ',oplist(nops)
           goto 110
        else
! set back one character as GPARx always increment last by 1 to bypass a ,
           last=last-1   
        endif
     endif
  endif
! save command for help path
  helprec%level=helprec%level+1
  helprec%cpath(helprec%level)=cbas(kom)
! The if loop is for handleing of defaults in submenu. "l ,,,,," took all 
! defaults but "l,,,,," did not ....
! if last>1 and cline(last-1:last-1) is a space and cline(last:last) a comma,
! increment last
  if(last.eq.1) then
     last=last+1
  elseif(last.lt.len(cline)) then
!     write(*,*)'menu increment 1: ',last
     if(cline(last:last).ne.' ') then
        if(cline(last+1:last+1).eq.',') last=last+1
     endif
!     write(*,*)'menu increment 2: ',last
  endif
  SELECT CASE(kom)
! command selection
!=================================================================
  CASE DEFAULT
     write(kou,*)'No such command'
     goto 100
!=================================================================
! amend subcommands
!       ['SYMBOL          ','ELEMENT         ','SPECIES         ',&
!        'PHASE           ','PARAMETER       ','REFERENCE       ',&
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
1020    continue
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
           call gparid('Antiferromagnetic factor: ',cline,last,jl,idef,q1help)
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
              call add_composition_set(iph,prefix,suffix,ics,ceq)
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
     case(6) ! amend reference
        call enter_reference_interactivly(cline,last,1,jl)
!-------------------------
     case(7) ! amend TPFUN symbol
        write(kou,*)'Not implemented yet'
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
2011    format(/'Calculating ',i3,' functions for T,P=',F10.2,1PE15.7/&
             3x,'No   F',11x,'F.T',9x,'F.P',9x,'F.T.T',7x,'F.T.P',7x,'F.P.P')
        call cpu_time(starting)
        do jl=1,notpf()
           call eval_tpfun(jl,ceq%tpval,val,ceq%eq_tpres)
           if(gx%bmperr.gt.0) goto 990
           write(kou,2012)jl,val
2012       format(I5,1x,6(1PE12.4))
        enddo
        call cpu_time(ending)
        write(kou,2013)ending-starting
2013    format('CPU time used: ',1pe15.6)
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
2031       format('G/N, dG/dT:',4(1PE16.8))
2032       format('G/N/RT, N:',2(1PE16.8))
!.......................................................
        case(2) ! calculate phase < >  G and dG/dy
           call calcg(iph,ics,1,lokres,ceq)
           if(gx%bmperr.ne.0) goto 990
           parres=>ceq%phase_varres(lokres)
           nofc=noconst(iph,ics,firsteq)
           write(kou,2031)(rgast*parres%gval(jl,1),jl=1,4)
           write(kou,2041)(rgast*parres%dgval(1,jl,1),jl=1,nofc)
2041       format('dG/dy:   ',4(1PE16.8))
!.......................................................
        case(3) ! calculate phase < > all
           call tabder(iph,ics,ceq)
           if(gx%bmperr.ne.0) goto 990
        END SELECT
!----------------------------------
     case(3) ! calculate equilibrium without initial global minimization
        if(minimizer.eq.1) then
! Lukas minimizer, first argiment=0 means do not use grid minimizer
!           call calceq1(0,ceq)
           write(kou,*)'Not implemented yet'
        else
           call calceq2(0,ceq)
! check that invmat allocated and stored
!           write(*,*)'inverted y: ',ceq%phase_varres(2)%cinvy(1,1)
        endif
!----------------------------------
     case(4) ! calculate ?
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
2101    format('N&x: ',F6.3,9F8.5)
! generate grid and find the phases and constitutions for the minimum.
        call global_gridmin(1,ceq%tpval,xknown,nv,iphl,icsl,aphl,nyphl,&
             yarr,cmu,ceq)
        if(gx%bmperr.ne.0) goto 990
! multiply amount of phases with total number of moles
        call scale_phase_amounts(xxx,ceq)
        write(kou,2102)'Stable phases: ',nv,(iphl(jl),icsl(jl),jl=1,nv)
2102    format(a,i2,10(2x,i3,i2))
!---------------------------------------------------------------
     case(7) ! calculate symbol
        call evaluate_all_svfun(kou,ceq)
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
!3035    continue
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
3046             format('Current status is ',a,' with ',1pe15.6,&
                      ' formula units.')
              else
                 write(kou,3047)text(1:i1)
3047             format('Current status is ',a)
              endif
           endif
           call gparcd(&
                'Suspend, Dormant, Entered, Fixed, Hidden or Not hidden?',&
                cline,last,1,ch1,'SUSPEND',q1help)
           nystat=99
           call capson(ch1)
           if(ch1.eq.'E') nystat=0
           if(ch1.eq.'S') nystat=1
           if(ch1.eq.'D') nystat=2
           if(ch1.eq.'F') nystat=3
           if(ch1.eq.'H') nystat=4
           if(ch1.eq.'N') nystat=5
           if(nystat.eq.99) then
              write(kou,*)'No such status'
              goto 100
           endif
           xxx=zero
           if(nystat.eq.0 .or. nystat.eq.3) then
              call gparrd('Amount: ',cline,last,xxx,zero,q1help)
           endif
           call change_phase_status(iph,ics,nystat,xxx,ceq)
           if(gx%bmperr.ne.0) goto 100
           if(iph.gt.0) then
              jl=get_phase_status(iph,ics,text,i1,xxy,ceq)
              if(gx%bmperr.ne.0) goto 100
              if(xxy.ge.zero) then
                 write(kou,3048)text(1:i1),xxy
3048             format('New status is ',a,' with ',1pe15.6,' formula units.')
              else
                 write(kou,3049)text(1:i1)
3049             format('New status is ',a)
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
        write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
     case(4) ! set LEVEL
        write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
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
! this means no reference phase
           iph=-1
        else
           call find_phase_by_name(name1,iph,ics)
           if(gx%bmperr.ne.0) goto 100
        endif
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
        call set_reference_state(iel,iph,tpa,ceq)
        if(gx%bmperr.eq.0) then
           write(kou,3104)
3104     format(' You must make a new calculation before the correct values'/&
                ' of chemical potentials or energy properties are shown.')
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
        kom3=submenu(cbas(kom),cline,last,csetph,nsetph,2)
        SELECT CASE(kom3)
        CASE DEFAULT
           write(kou,*)'Set phase status subcommand error'
           goto 100
!............................................................
        case(1) ! SET PHASE CONSTITUTION <name> <comp.set> ....
           call ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
           if(gx%bmperr.ne.0) goto 990
!           goto 100
!............................................................
! copied from 3045
        case(2) ! SET PHASE STATUS <phase> <status>
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
              else
                 write(kou,3047)text(1:i1)
              endif
           endif
           call gparcd(&
                'Suspend, Dormant, Entered, Fixed, Hidden or Not hidden?',&
                cline,last,1,ch1,'SUSPEND',q1help)
           nystat=99
           call capson(ch1)
           if(ch1.eq.'E') nystat=0
           if(ch1.eq.'S') nystat=1
           if(ch1.eq.'D') nystat=2
           if(ch1.eq.'F') nystat=3
           if(ch1.eq.'H') nystat=4
           if(ch1.eq.'N') nystat=5
           if(nystat.eq.99) then
              write(kou,*)'No such status'
              goto 100
           endif
           xxx=zero
           if(nystat.eq.0 .or. nystat.eq.3) then
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
        case(3:4) !set phase default constitution wildcard allowed, also AMOUNT
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
           if(kom3.eq.3) then
! set phase default constituntion
              call set_default_constitution(iph,ics,0,ceq)
           else
! set phase amount
              call gparrd('Amount: ',cline,last,xxx,zero,q1help)
              call set_phase_amounts(iph,ics,xxx,ceq)
           endif
!............................................................
! subsubsub command
        case(5) ! set phase bits
           call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
           call find_phase_by_name(name1,iph,ics)
           if(gx%bmperr.ne.0) goto 990
           call get_phase_record(iph,lokph)
           kom4=submenu('Set which bit?',cline,last,csetphbits,nsetphbits,8)
           SELECT CASE(kom4)
           CASE DEFAULT
              write(kou,*)'Set phase bit subcommand error'
              goto 100
!............................................................
           case(1) ! FCC_PERMUTATIONS
! if check returns .true. it is not allowed to set FORD or BORD
              if(check_minimal_ford(lokph)) goto 100
              call set_phase_status_bit(lokph,PHFORD)
           case(2) ! BCC_PERMUTATIONS
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
           case(8) ! NO_AUTO_COMP_SET, not allowed to create compsets automatic
              call set_phase_status_bit(lokph,PHNOCS)
           case(9) ! ELASTIC_MODEL_A
! set by amend??
!              call set_phase_status_bit(lokph,PHNOCS)
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
!        goto 100
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
        call gparcd('Independent axis variable',cline,last,1,text,'T ',q1help)
! check there is a condition ...
        iax=1
        maxax=1
        axvar(iax)='T'
        dmin=axvalold(1,iax)
        call gparrd('Min value:',cline,last,xxx,dmin,q1help)
        if(buperr.ne.0) goto 100
        axval(1,iax)=xxx
        dmax=axvalold(2,iax)
        call gparrd('Max value:',cline,last,xxx,dmax,q1help)
        if(buperr.ne.0) goto 100
        axval(2,iax)=xxx
! default step 1/40
        dinc=axvalold(3,iax)
        call gparrd('Increment:',cline,last,xxx,dinc,q1help)
        if(buperr.ne.0) goto 100
        axval(3,iax)=xxx
        axvalold(1,iax)=axval(1,iax)
        axvalold(2,iax)=axval(2,iax)
        axvalold(3,iax)=axval(3,iax)
!  write(*,3602)(axval(i,iax),i=1,3)
3602    format(/'axlimits: ',3(1pe12.4))
!-------------------------------------------------------------
     case(15) ! set input amounts
        call set_input_amounts(cline,last,ceq)
!-------------------------
     case(16)
        write(kou,*)'Not implemented yet'
!-------------------------
     case(17)
        write(kou,*)'Not implemented yet'
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
        goto 100
!---------------------------------------------------------------
     CASE(1) ! enter TPFUN symbol (constants, functions, tables)
        call gparc('Symbol name: ',cline,last,1,name1,' ',q1help)
        if(buperr.ne.0) goto 990
!  if(badsymname(name1)) then
        if(.not.proper_symbol_name(name1,0)) then
           write(kou,*)'Bad symbol name'
           goto 990
        endif
        call gparcd('Function, value or table? ',cline,last,1,name2,&
             'FUNCTION ',q1help)
        if(buperr.ne.0) goto 990
        call capson(name2)
        if(compare_abbrev(name2,'FUNCTION ')) then
           call enter_tpfun_interactivly(cline,last,string,jp)
           if(gx%bmperr.ne.0) goto 990
           call enter_tpfun(name1,string,lrot)
           if(gx%bmperr.ne.0) goto 990
        elseif(compare_abbrev(name2,'VALUE ')) then
           write(kou,*)'Not implemented yet'
        elseif(compare_abbrev(name2,'TABLE ')) then
           write(kou,*)'Not implemented yet'
        else
           write(kou,*)'No such type of symbol'
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
        call gparr('Element H298-H0: ',cline,last,h298,zero,q1help)
        call gparr('Element S298: ',cline,last,s298,one,q1help)
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
        phtype='S'
        call gpari('Number of sublattices: ',cline,last,nsl,1,q1help)
        if(buperr.ne.0) goto 990
        if(nsl.gt.10) then
           write(kou,*)'Maximum 10 sublattices'
           goto 100
        endif
        icon=0
        sloop: do ll=1,nsl
! 'Number of sites on sublattice xx: '
!  123456789.123456789.123456789.123
           once=.true.
4042       continue
           write(quest1(31:32),4043)ll
4043       format(i2)
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
4045       continue
           call gparc('All Constituents: ',cline,last,4,text,';',q1help)
           if(buperr.ne.0) goto 990
           knr(ll)=0
           jp=1
4047       continue
           if(eolch(text,jp)) goto 4049
           call getname(text,jp,name3,1,ch1)
           if(buperr.eq.0) then
              icon=icon+1
              const(icon)=name3
              knr(ll)=knr(ll)+1
              goto 4047
           elseif(once) then
              write(kou,*)'Input error ',buperr,', please reenter'
              buperr=0; once=.false.; goto 4045
           else
              goto 100
           endif
           buperr=0
4049       continue
        enddo sloop
        model='CEF-RKM'
        call new_phase(name1,nsl,knr,const,sites,model,phtype)
        if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
     case(5) ! enter parameter is always allowed
        call enter_parameter_interactivly(cline,last)
        if(gx%bmperr.ne.0) goto 990
!---------------------------------------------------------------
     case(6) ! enter reference
        call enter_reference_interactivly(cline,last,0,jl)
        if(gx%bmperr.ne.0) goto 990
        write(kou,*)'Reference number is ',jl
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
! enter not used
     case(12)
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
        write(kou,*)'List subcommand error'
        goto 100
!-----------------------------------------------------------
     case(1) ! list data for everything
        call list_all_elements(kou)
        if(gx%bmperr.ne.0) goto 990
        call list_all_species(kou)
        if(gx%bmperr.ne.0) goto 990
        call list_all_funs(kou)
        if(gx%bmperr.ne.0) goto 990
        do iph=1,noph()
           call list_phase_data(iph,kou)
           if(gx%bmperr.ne.0) goto 990
        enddo
! list reference phase last
        iph=0
        call list_phase_data(0,kou)
! finally list data references
        write(kou,*)
        call list_references(kou)
!-----------------------------------------------------------
     case(2) ! list short with status bits
        write(kou,6022)globaldata%name,globaldata%rgasuser,&
             globaldata%pnorm,globaldata%status
6022    format('System',18x,'Gas constant Pressure norm',23x,'Status'/&
             a,1pe12.4,2x,1pe12.4,21x,z8)
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
6051       format('Output for equilibrium: ',i3,', ',a)
           call list_phase_results(iph,ics,mode,kou,ceq)
           if(gx%bmperr.ne.0) goto 990
!...............................................................
        case(3) ! list phase model (including disordered fractions)
           call list_phase_model(iph,ics,kou,ceq)
        END SELECT
!------------------------------
     case(4,10)  ! list state variable or parameter identifier value, loop.
6099 continue
        if(btest(ceq%status,EQNOEQCAL) .or. btest(ceq%status,EQFAIL)) then
           write(kou,6101)
6101       format(' *** Warning,',&
                'equilibrium not calculated, values are probably wrong')
        elseif(btest(ceq%status,EQINCON)) then
           write(kou,6102)
6102       format(' *** Warning, values can be inconsistent with',&
                ' current conditions')
        endif
6105    continue
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
! in model the state variable is returned as generated by the program
        call get_state_var_value(line,xxx,model,ceq)
        if(gx%bmperr.ne.0) then
           write(kou,*)'Error code ',gx%bmperr
           gx%bmperr=0; goto 6099
        endif
        write(kou,6108)model(1:len_trim(model)),xxx
6108    format(1x,a,'=',1PE15.7)
        goto 6105
!-----------------------------------------------------------
     case(5) ! list data references
        call list_references(kou)
!-----------------------------------------------------------
     case(6) ! list parameter symbols
        call list_defined_properties(kou)
!-----------------------------------------------------------
     case(7) ! list axis
        write(kou,6131)
6131    format(4x,'Axis variable',12x,'Minimum',5x,'Maximum',5x,'Increment')
        do iax=1,maxax
           write(kou,6132)iax,axvar(iax),(axval(jl,iax),jl=1,3)
6132       format(i2,2x,a,3(1pe12.4))
        enddo
!-----------------------------------------------------------
     case(8) ! list tpfun symbol
        call list_all_funs(kou)
!------------------------------------------------------------
     case(9) ! list quit
        goto 100
!------------------------------------------------------------
!     case(10) ! list parameter identifiers, same as state_variables
!        goto 100
!-----------------------------------------------------------
     case(11) ! list equilibria
        do iel=1,noeq()
           write(kou,6202)iel,eqlista(iel)%eqname
        enddo
6202    format(i5,2x,a)
!------------------------------
     case(12) ! list results
        call gparid('Output mode: ',cline,last,listresopt,lrodef,q1help)
        lrodef=listresopt
        write(kou,6051)ceq%eqno,ceq%eqname
!  if(btest(globaldata%status,GSEQFAIL)) then
        if(btest(ceq%status,EQFAIL)) then
           write(kou,6305)
6305       format(/' *** The results listed are not valid',&
                ' as last calculation failed'/)
        elseif(btest(globaldata%status,GSNOPHASE)) then
           write(kou,*)'No results as no data'
           goto 100
!  elseif(btest(globaldata%status,GSNOEQCAL)) then
        elseif(btest(ceq%status,EQNOEQCAL)) then
           write(kou,6307)
6307       format(/' *** The results listed does not represent',&
                ' a calculated equilibrium'/)
        elseif(btest(ceq%status,EQINCON)) then
           write(kou,6306)
6306       format(/' *** The results listed may be inconsistent',&
                ' with the current conditions'/)
        endif
        write(kou,6302)'Conditions .........'
6302    format(a,40('.'),':')
6303    format(/a,40('.'),':')
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
        if(listresopt.eq.1) then
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
           mode=10
        elseif(listresopt.eq.9) then
! all phases with mole fractions and constitution in alphabetical order
           mode=110
        else
! all phase with with mole fractions
           mode=0
        endif
        ics=1
        do iph=1,noph()
           ics=0
6310       continue
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
        call list_conditions(kou,ceq)
!------------------------------
     case(14) ! list symbols (state variable functions, not TP funs)
        call list_all_svfun(kou,ceq)
!------------------------------
! list ?
     case(15)
        goto 100
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
! all records must be removed and init_gtp called
           call new_gtp
           if(gx%bmperr.ne.0) goto 990
        endif
     endif
     kom2=submenu(cbas(kom),cline,last,cread,ncread,2)
     call gparc('File name: ',cline,last,1,filename,' ',q1help)
     SELECT CASE(kom2)
!-----------------------------------------------------------
     CASE DEFAULT
        write(kou,*)'Read subcommand error'
        goto 100
!-----------------------------------------------------------
     case(1) ! read unformatted file created by SAVE
        call gtpread(filename,text)
        if(gx%bmperr.ne.0) goto 990
        kl=len_trim(text)
        if(kl.gt.1) then
           write(kou,8110)text(1:kl)
        endif
8110    format(/'Savefile text: ',a/)
!---------------------------------------------------------
! read tdb
!8200 continue
     case(2) ! read TDB
        call readtdb(filename)
!-----------------------------------------------------------
8300 continue
     case(3) ! read ?
        goto 100
     end SELECT
!=================================================================
! save unformatted
  case(9)
     kom2=submenu(cbas(kom),cline,last,csave,ncsave,1)
     if(kom2.le.0 .or. kom2.gt.ncsave) goto 100
     call gparc('File name: ',cline,last,1,filename,' ',q1help)
     jp=0
     kl=index(filename,'.')
     if(kl.le.0) then
        jp=len_trim(filename)
     elseif(filename(kl+1:kl+1).eq.' ') then
! just ending a filename with . not accepted as extention
        jp=kl
     endif
     if(kl.le.0 .and. jp.le.0) then
        write(kou,*)'Missing file name'
        goto 100
     endif
     SELECT CASE(kom2)
!-----------------------------------------------------------
     CASE DEFAULT
        write(kou,*)'save subcommand error'
!-----------------------------------------------------------
     case(1) ! save unformatted
        if(jp.gt.0) filename(jp+1:)='.ocu '
        text='U '
        call gtpsave(filename,text)
!-----------------------------------------------------------
     case(2) ! save TDB
        if(jp.gt.0) filename(jp+1:)='.TDB '
        text='T '
        call gtpsave(filename,text)
!-----------------------------------------------------------
     case(3) ! save MACRO
        if(jp.gt.0) filename(jp+1:)='.BMM '
        text='M '
        call gtpsave(filename,text)
     end SELECT
!=================================================================
! help ... just list commands
  case(10)
!     call nghelp(cline,last,cbas,ncbas)
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
! new, same as reinitiate
  case(13)
! one must deallocate everyting and init_gtp is called inside new_gtp
     call gparcd('All data will be removed, are you sure?',cline,last,&
          1,ch1,'N',q1help)
     if(ch1.ne.'Y') then
        write(kou,*)'*** NO CHANGE, upper case Y needed for NEW'
        goto 100
     endif
! this routine is not updated and will return error code
     call new_gtp
!=================================================================
! macro begin
  case(14)
     call macbeg(cline,last,logok)
     if(buperr.ne.0 .or. gx%bmperr.ne.0) goto 990
!=================================================================
! about
  case(15)
     write(kou,15010)
15010 format(/'This is Open Calphad (OC), a free software for',&
           ' thermodynamic calculations'/&
           'This software is protected by the GNU General Public License'/&
           'You may freely distribute copies as long as you also provide ',&
           'the source code.'/'The software is provided "as is" without ',&
           'any warranty of any kind, either'/'expressed or implied.'//&
           'The full license text is provided with the software or can be ',&
           'obtained from'/'the Free Software Foundation http://www.fsf.org'//&
           'Copyright 2010-2012, several persons.'/&
           'Contact person Bo Sundman, bo.sundman@gmail.com'/)
     goto 100
!=================================================================
! debug subcommands
  case(16)
     kom2=submenu(cbas(kom),cline,last,cdebug,ncdebug,0)
     SELECT CASE (kom2)
!------------------------------
     CASE DEFAULT
        write(kou,*)'Default case ',kom2
!------------------------------
! debug free lists
     CASE(1)
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
        call gparc('Give name or number?',cline,last,1,text,' ',q1help)
        if(buperr.ne.0) goto 990
! if number select that
        jl=1
        call getint(text,jl,i1)
        if(buperr.ne.0) then
           buperr=0
           call findeq(text,ieq)
           if(gx%bmperr.ne.0) goto 990
           neqdef=ieq
!     write(kou,*)ieq,text,eqlista(ieq)%eqname
           ceq=>eqlista(ieq)
        else
           call selecteq(i1,ceq)
           if(gx%bmperr.ne.0) goto 990
           neqdef=i1
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
     write(kou,18010)
18010 format(' *** Warning, this command will delete the data for the',&
           ' element, species or'/' phase specified and the data cannot',&
           ' be recovered unless read again from'/' file.  If you',&
           ' only want to temporarily remove some data use QUIT'/&
           ' from this command and then SET STATUS'/)
     kom2=submenu(cbas(kom),cline,last,crej,nrej,3)
     SELECT CASE(kom2)
!-----------------------------------------------------------
     CASE DEFAULT
        write(kou,*)'Delete subcommand error'
        goto 100
!-----------------------------------------------------------
! delete element
     case(1)
        write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! delete species
     case(2)
        write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! delete phase
     case(3)
        write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
     case(4) !quit
        goto 100
!-----------------------------------------------------------
! delete ??
     case(5)
        write(kou,*)'Not implemented yet'
!-----------------------------------------------------------
! delete ??
     case(6)
        write(kou,*)'Not implemented yet'
     end SELECT
!=================================================================
! step
  case(19)
     call gparcd('Option?',cline,last,1,text,'NORMAL ',q1help)
     if(associated(resultlist)) then
        write(kou,*)'There are some results already form step or map'
        call gparcd('Reinitiate?',cline,last,1,ch1,'Y',q1help)
        if(ch1.eq.'y' .or. ch1.eq.'Y') then
           nullify(resultlist)
        endif
     endif
! step will use axval(1..3,1) as min, max and step
     axvar(1)='T'
     call sstep1(text,axvar,axval,resultlist,ceq)
! mark that conditions and any listed results may be inconsistent
     ceq%status=ibset(ceq%status,EQINCON)
!=================================================================
! map
!20000 continue
  case(20)
     write(kou,*)'Not implemented yet'
!=================================================================
! plot
  case(21)
! the 4th argument to gparc means the following:
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
     call gparcd('Dependent axis variable',cline,last,7,text,'NP(*) ',q1help)
     call gparcd('File name',cline,last,1,plotfile,'ocg',q1help)
     form=' '
     call splot1(text,plotfile,resultlist,form)
     if(gx%bmperr.ne.0) goto 100
     call gparc('Hardcopy (P for postscript)?',cline,last,1,form,'none',q1help)
     if(form.ne.'none') then
        call capson(form)
        call splot1(text,plotfile,resultlist,form)
     endif
!=================================================================
! HPCALC
  case(22)
     call hpcalc
     buperr=0
!=================================================================
! FIN, do not ask, the French always know what they do ...
  case(23)
     if(logfil.gt.0) then
        write(logfil,*)'set interactive'
     endif
     call openlogfile(' ',' ',-1)
!     call gparcd('Are you sure?',cline,last,1,ch1,'N',q1help)
!     if(ch1.eq.'y' .or. ch1.eq.'Y') then
        stop 'Au revoir'
!     endif
!=================================================================
! unused
  case(24)
     write(kou,*)'No such command'
!=================================================================
!
  END SELECT
! command executed, prompt for another command unless error code
  if(gx%bmperr.eq.0) goto 100
!============================================================
! handling errors
990  continue
  write(kou,*)
  write(kou,*)'Error codes: ',gx%bmperr,buperr
  if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
     write(kou,*)bmperrmess(gx%bmperr)
  endif
  if(stop_on_error) then
! turn off macro but remain inside software
     call macend(cline,last,logok)  
     write(kou,*)'Stop_on_error set, press return to finish program'
     read(kiu,17)ch1
17   format(a)
     stop
  endif
  gx%bmperr=0; buperr=0
  goto 100
!
end subroutine oc_command_monitor

integer function submenu(query,cline,last,ccomnd,ncomnd,kdef)
! general subcommand decoder
!  implicit double precision (a-h,o-z)
  implicit none
  character cline*(*),ccomnd(*)*(*),query*(*)
  character defansw*16,query1*64,text*256
  integer last,kdef,ncomnd,kom2,lend,lenq
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
10 format(a,a,i4,': ',a)
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
!        write(*,102)'submenu 7: ',last,cline(1:last+5)
102     format(a,i5,'"',a,'"')
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
     helprec%level=helprec%level+1
     helprec%cpath(helprec%level)=ccomnd(kom2)
  endif
1000 continue
  return
end function submenu

END MODULE oc_cmon1

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
