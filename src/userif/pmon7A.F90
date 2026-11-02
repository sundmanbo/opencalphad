!
! pmon6A.F90 is included in pmon6.F90 (module cmon1oc) after the
! subroutines oc_command_monitor and ocmon_set_options/reset_options
!
! It contains code extracted from oc_command_monitor to make the command
! monitor shorter.  All selection of commands and subcommands remain in
! oc_command_monitor, these subroutines are called from there.
! Data kept between commands are module variables in cmon1oc.
! A return from these subroutines means returning to the command level,
! if gx%bmperr or buperr is nonzero the error is handled there.
!
! 48 subroutines
!  subroutine ui_set_status_phases(cline,last,ceq)
!  subroutine ui_set_working_directory(cline,last,kom3)
!  subroutine ui_set_mqmqa_debug(cline,last)
!  subroutine ui_set_reference_state(cline,last,ceq)
!  subroutine ui_set_phase_status(iph,ics,cline,last,ceq)
!  subroutine ui_set_logfile(cline,last,logfil)
!  subroutine ui_set_weight(cline,last,ceq)
!  subroutine ui_set_numeric_options(cline,last,ceq)
!  subroutine ui_set_axis(cline,last,maxax,ceq)
!  subroutine ui_set_bit_equilibrium(cline,last,ceq)
!  subroutine ui_set_bit_global(cline,last)
!  subroutine ui_set_optcoeff_variable(cline,last)
!  subroutine ui_set_range_exper_equil(cline,last)
!  subroutine ui_set_optcoeff_fixed(cline,last)
!  subroutine ui_enter_tpfun(cline,last)
!  subroutine ui_enter_optimize_coeff(cline,last)
!  subroutine ui_enter_gnuplot_terminal(cline,last)
!  subroutine ui_read_tdb(cline,last,kom2,tdbfile,ocbase)
!  subroutine ui_read_xtdb(cline,last,xtdbfile)
!  subroutine ui_save_unformatted(cline,last,ocufile,model)
!  subroutine ui_list_state_variables(cline,last,kom,kom2,lut,ceq)
!  subroutine ui_list_equilibria(kom2,lut,ceq)
!  subroutine ui_list_results(cline,last,version,listresopt,lrodef,lut,ceq)
!  subroutine ui_debug_symbol_value(cline,last,ceq)
!  subroutine ui_amend_phase_magnetic(cline,last,lokph)
!  subroutine ui_amend_phase_disordered(cline,last,iph)
!  subroutine ui_calculate_tpfun(cline,last,lut,ceq)
!  subroutine ui_calculate_phase_adjust(cline,last,iph,ics,kom3,lut,ceq)
!  subroutine ui_amend_phase_asymmetries(cline,last,iph,ceq)
!  subroutine ui_calculate_all_equilibria(cline,last,lut)
!  subroutine ui_amend_assessment_result(cline,last)
!  subroutine ui_amend_optimizing_coeffs(cline,last,lut)
!  subroutine ui_optimize(cline,last,lut)
!  subroutine ui_new(cline,last,defcp,reinit)
!  subroutine ui_step_tzero(cline,last,tzcond,name2,ceq)
!  subroutine ui_step_scheil(cline,last,kom2,ceq)
!  subroutine ui_map(cline,last,ceq)
!  subroutine ui_plot_format_file(cline,last,kom2)
!  subroutine ui_plot_append(cline,last,string)
!  subroutine ui_plot_text(cline,last,ceq)
!  subroutine ui_list_excel_csv(cline,last,version,ceq)
!  subroutine ui_plot_scale_ranges(cline,last,jump)
!  subroutine ui_plot_axes(cline,last,jump,ceq)
!  subroutine ui_delete_step_map(kind,cline,last,ceq)
!

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_status_phases(cline,last,ceq)
! SET STATUS PHASE, several phases can be given
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,line*80,string*256,ch1*1
    integer iph,j4,ll,nystat
    double precision xxx
! Now allow multiple phase names and *S, *D and *E
! argument 5 means whole input line
    call gparcx('Phase name(s): ',cline,last,5,line,'=',&
         '?Set status phase')
    string=line
3017 continue
    ll=index(string,'=')
    if(ll.eq.0) then
       call gparcx('More phase name(s): ',cline,last,5,line,'=',&
            'Set status phase')
       string(len_trim(string)+2:)=line
       goto 3017
    endif
! exttract first letter after = (if any)
    j4=ll
    call getext(string,j4,1,name1,' ',iph)
    ch1=name1(1:1)
! if user has given "=e 0" then keep the amount is cline
    cline=string(j4:)
    string(ll:)=' '
    if(ch1.eq.' ') then
! if ll==1 then input was finished by equal sign, ask for status
       call gparcdx(&
            'New status S(uspend), D(ormant), E(ntered) or F(ixed)?',&
            cline,last,1,name1,'E','?Set status phase')
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
! no longer available if(ch1.eq.'N') nystat=5
    if(nystat.eq.99) then
       write(kou,*)'No such status'
       return
    endif
    xxx=zero
    if(nystat.eq.phentered .or. nystat.eq.phfixed) then
       call gparrdx('Amount: ',cline,last,xxx,zero,&
            'Set status phase amount')
    endif
    call change_many_phase_status(string,nystat,xxx,ceq)
    return
  end subroutine ui_set_status_phases

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_working_directory(cline,last,kom3)
! SET ADVANCED WORKING_DIRECTORY
    character cline*(*)
    integer last
    integer kom3
!\end{verbatim}
    character string*256,ch1*1
    integer i1,j4,ztyp
    logical logok
    write(kou,*)'Current working directory: ',trim(workingdir)
    write(kou,'(a)')'To change directory please select an OCM file'
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT
! IMPORTANT: extensions also in utilities/ftinydialog.F90/ftinyopen !!!!
! try to set current working directory as input to allow editing
! The promt here is never displayed ...
    ztyp=3
    call gparfilex('Select new working directory',&
         cline,last,1,string,' ',ztyp,'?Set adv workdir')
!
    inquire(file=string,exist=logok)
    if(.not.logok) then
       write(*,*)'No such directory, working directory not changed'
    elseif(trim(workingdir).ne.trim(string)) then
! strip away any file name (up to last / or \)
       j4=len_trim(string)
       ch1=string(j4:j4)
       do while(j4.gt.0 .and. .not.(ch1.eq.'/' .or. ch1.eq.'\'))
          j4=j4-1
          ch1=string(j4:j4)
       enddo
       string(j4:)=' '
! this is a gfortran special extension
       call chdir(string,i1)
       if(i1.ne.0) then
          write(*,*)'Failed change working directory',i1
       else
          write(*,'(a,a,i3)')'New working directory: ',&
               trim(string),kom3
          workingdir=string
       endif
    endif
    return
  end subroutine ui_set_working_directory

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_mqmqa_debug(cline,last)
! SET ADVANCED DEBUG, MQMQA debug flags
    character cline*(*)
    integer last
!\end{verbatim}
    integer i1,i2
! select macro and set default directory to the one with the macro for TDB
! check consistency of some data in mqma_data!!
! some also in gtp_mqmqa_var!!
! the command AMEND DEBUG now has this question
! temporarily it may remain part of LIST MQMQA <phase> DEBUG (and still is)
    write(*,1687)mqmqa_data%nconst,mqmqa_data%nquad,&
         mqmqa_data%ncon1,mqmqa_data%ncat,&
         mqmqa_data%ncon2,mqmqa_data%nan,&
         mqmqa_data%npair,mqmqa_data%lcat
1687 format(/'Values of some duplicate global data:',&
         'Number of quads: ',2i4/&
         'Number of cations: ',2i4/&
         'Number of anions: ',2i4/&
         'Number of pairs: ',i4/&
      'Value of cation*(cation+1)/2: ',i4/)
! the format labels currently also at the LIST MQMQA command
    i2=0
    mqmqdebug=.false.
    mqmqdebug2=.false.
    mqmqtdb=.false.
    mqmqxcess=.false.
    mqmqder=.false.
    mqverbose=.false.
    mqter=.false.
    mqenter=.false.
    write(*,1688)
1688 format('Debug settings:'/'0 no debug'/'1 debug asymmetry'/&
         '2 debug some more'/&
         '3 debug reading TDB'/&
         '4 debug parameter calculation'/&
         '5 debug partial derivative calculation'/&
         '6 verbose listing of asymmetries'/&
         '7 debug ternary'/&
         '8 debug enter ternary parameter')
    call gparidx('Set mqmqa debug?',cline,last,i1,i2,'?MQMQA debug')
! set several ....
1689 continue
    if(i1.eq.1) then
! asymmetry debug
       mqmqdebug=.true.
    elseif(i1.eq.2) then
       mqmqdebug2=.true.
    elseif(i1.eq.3) then
! tdb reading debug
       mqmqtdb=.true.
    elseif(i1.eq.4) then
! calculation debug
       mqmqxcess=.true.
    elseif(i1.eq.5) then
! partial derivative debug
       mqmqder=.true.
    elseif(i1.eq.6) then
! partial derivative debug more
       mqverbose=.true.
    elseif(i1.eq.7) then
! ternary debug
       mqter=.true.
    elseif(i1.eq.8) then
! bug entering ternary parameter
       mqenter=.true.
    endif
    call gparidx('Set more debug?',cline,last,i1,i2,'?MQMQA debug')
    if(i1.gt.0) goto 1689
    write(*,1676)mqmqdebug,mqmqdebug2,mqmqtdb,&
         mqmqxcess,mqmqder,mqverbose,mqter,mqenter
1676 format('MQ debug setting:',8l2)
    continue
    return
  end subroutine ui_set_mqmqa_debug

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_reference_state(cline,last,ceq)
! SET REFERENCE_STATE
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24
    integer iel,iph,ics
    double precision xxx,xxy,tpa(2)
    call gparcx('Component name: ',cline,last,1,name1,' ',&
         '?Set reference phase')
    call find_component_by_name(name1,iel,ceq)
    if(gx%bmperr.ne.0) return
    call gparcx('Reference phase: ',cline,last,1,name1,'SER ',&
         '?Set reference phase')
    if(name1(1:4).eq.'SER ') then
       write(kou,*)'Reference state is stable phase at 298.15 K and 1 bar'
! this means no reference phase, SER is at 298.15K and 1 bar
       iph=-1
    else
       call find_phase_by_name(name1,iph,ics)
       if(gx%bmperr.ne.0) return
! temperature * means always to use current temperature
       xxy=-one
       call gparrx('Temperature: /*/: ',cline,last,xxx,xxy,&
            '?Set reference phase')
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
       call gparrdx('Pressure: ',cline,last,xxx,xxy,&
            '?Set reference phase')
       if(xxx.le.zero) then
          tpa(2)=xxy
       else
          tpa(2)=xxx
       endif
    endif
    call set_reference_state(iel,iph,tpa,ceq)
    if(gx%bmperr.eq.0) then
3104   format(' You may have to make a new calculation before the',&
            ' correct values'/&
            ' of chemical potentials or other properties are shown.')
    endif
    return
  end subroutine ui_set_reference_state

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_phase_status(iph,ics,cline,last,ceq)
! SET PHASE <phase> STATUS, iph negative means all phases
    integer iph,ics
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character text*72,ch1*1
    integer i1,j4,nystat
    double precision xxx,xxy
    if(iph.gt.0) then
       j4=get_phase_status(iph,ics,text,i1,xxx,ceq)
       if(gx%bmperr.ne.0) return
       if(xxx.ge.zero) then
          write(kou,3046)text(1:i1),xxx
3046      format('Current status is ',a,' with ',1pe15.6,&
               ' formula units.')
       else
          write(kou,3047)text(1:i1)
3047      format('Current status is ',a)
       endif
    endif
    call gparcdx(&
         'Suspend, Dormant, Entered, Fixed, Hidden or Not hidden?',&
         cline,last,1,ch1,'SUSPEND','?Set phase status')
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
       return
    endif
    xxx=zero
    if(nystat.eq.phentered .or. nystat.eq.phfixed) then
       call gparrdx('Amount: ',cline,last,xxx,zero,&
            '?Set phase amount')
    endif
    call change_phase_status(iph,ics,nystat,xxx,ceq)
    if(gx%bmperr.ne.0) return
    if(iph.gt.0) then
       j4=get_phase_status(iph,ics,text,i1,xxy,ceq)
       if(gx%bmperr.ne.0) return
       if(xxy.ge.zero) then
          write(kou,3048)text(1:i1),xxy
3048      format('New status is ',a,' with ',1pe15.6,&
               ' formula units.')
       else
          write(kou,3049)text(1:i1)
3049      format('New status is ',a)
       endif
    else
       write(kou,*)'New status set for all phases'
    endif
    return
  end subroutine ui_set_phase_status

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_logfile(cline,last,logfil)
! SET LOG_FILE, logfil is the unit of the log file, 0 if none
    character cline*(*)
    integer last
    integer logfil
!\end{verbatim}
    character name1*24,line*80,model*72
    integer ztyp
    logical maptopbug
! tinyfiles_dialog has difficult returning a non-existant file name
! the argument "-8" means open a log file for output
    maptopbug=.true.
    if(associated(maptop)) then
       maptopbug=.false.
    endif
    ztyp=-8
    call gparfilex('Log file name: ',cline,last,1,model,'oclog',ztyp,&
         '?Set logfile')
    name1=model(1:5)
    call capson(name1)
    if(maptopbug .and. associated(maptop)) then
! for unkown reason maptop has become associated here but was not 8 lines above!
       nullify(maptop)
       write(*,*)'PMON clear link to maptop'
    endif
    if(name1(1:5).eq.'NONE ') then
! close log file
       call openlogfile(' ',' ',-1)
       logfil=0
       write(*,*)'Log file closed'
    else
       if(len_trim(model).eq.0) then
          model='OCLOG.LOG'
       elseif(index(model,'.LOG ').eq.0) then
          model=trim(model)//'.LOG'
       endif
       write(*,*)'Setting logfile to: "',trim(model),'"'
       call gparcx('Title: ',cline,last,5,line,' ','?Set logfile')
       call openlogfile(model,line,39)
       if(buperr.ne.0) then
          write(kou,*)'Error opening logfile: ',buperr
          logfil=0
       else
          write(*,'(a,a)')'Commands will be logged in file ',trim(model)
          logfil=39
       endif
    endif
    return
  end subroutine ui_set_logfile

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_weight(cline,last,ceq)
! SET WEIGHT of experimental equilibria
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24
    integer i1,i2,j4,ll
    double precision xxx
    if(.not.allocated(firstash%eqlista)) then
       write(kou,*)'You must first set a range of experimental equilibria'
       return
    endif
! NOTE mexp must be updated to the correct number of EXPERIMENTS
! that is done by OPTIMIZE
    updatemexp=.true.
    mexp=0
! the residuals from the last optimization are no longer valid
    if(allocated(errs)) deallocate(errs)
    call gparrdx('Weight ',cline,last,xxx,one,'?Set weight')
    if(buperr.ne.0) return
! The weight must be 0 or positive
    xxx=abs(xxx)
    call gparcdx('Equilibria (abbrev name) or range: ',cline,last,&
         1,name1,'CURRENT','?Set weight')
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
!                setwei: do j4=i1,i2
          setwei: do j4=1,size(firstash%eqlista)
             if(firstash%eqlista(j4)%p1%eqno.ge.i1 .and. &
                  firstash%eqlista(j4)%p1%eqno.le.i2) then
                firstash%eqlista(j4)%p1%weight=xxx
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
3066   format('Changed weight for ',i5,' equilibria')
    endif
    return
  end subroutine ui_set_weight

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_numeric_options(cline,last,ceq)
! SET NUMERIC_OPTIONS
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer i1,i2
    double precision xxx,xxy
    i2=ceq%maxiter
    call gparidx('Max number of iterations: ',cline,last,i1,i2,&
         '?Set numeric')
    if(i1.gt.0) then
       ceq%maxiter=i1
    endif
!------------
    xxx=ceq%xconv
    call gparrdx('Max error in fraction: ',cline,last,xxy,xxx,&
         '?Set numeric')
!CCI
    if(xxy.gt.default_minxconv) then
       ceq%xconv=xxy
    else
       ceq%xconv=default_minxconv
    endif
!CCI
!------------ what is this? not used in gtp3X.F90
    xxx=ceq%gdconv(1)
    call gparrdx('Max cutoff driving force: ',cline,last,xxy,xxx,&
         '?Set numeric')
    if(xxy.gt.default_mingdconv) then
       ceq%gdconv(1)=xxy
!CCI
    else
       ceq%gdconv(1)=default_mingdconv
!CCI
    endif
!------------ if the point between two gridpoints in a phase is less then merge
    xxx=ceq%gmindif
    call gparrdx('Min difference merging gridpoints: ',cline,last,&
         xxy,xxx,'?Set numeric')
!CCI
!strange old value was 1000 times lower -1e-5 in the test vs -1.e-2 as default
    if(xxy.gt.default_mingridmin) then
       ceq%gmindif=xxy
    else
       ceq%gmindif=default_mingridmin
    endif
    return
  end subroutine ui_set_numeric_options

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_axis(cline,last,maxax,ceq)
! SET AXIS, maxax is the maximum number of axis
    character cline*(*)
    integer last
    integer maxax
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,text*72
    integer i1,i2,iax,jp
    double precision xxx,dmin,dmax,dinc,coeffs(10)
    logical once
    TYPE(gtp_state_variable), target :: stvrvar
    TYPE(gtp_state_variable), pointer :: stvr
    TYPE(gtp_condition), pointer :: pcond,firstc
    if(btest(globaldata%status,GSNOPHASE)) then
       write(kou,*)'You have no data!'
       return
    endif
    i1=noofaxis+1
    call gparidx('Axis number',cline,last,iax,i1,'?Set axis')
    if(iax.lt.1 .or. iax.gt.maxax) then
       write(kou,3300)maxax
3300   format('Axis number must be between 1 and ',i1)
       return
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
       if(gx%bmperr.ne.0) return
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
    call gparcdx('Condition varying along axis: ',cline,last,1,&
         text,name1,'?Set axis')
    call capson(text)
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
       return
    elseif(trim(name1).eq.trim(text)) then
! check if same variable, quit this IF loop
       continue
    else ! changed axis variable, set default limits
       dmin=zero
       dmax=one
       i1=len_trim(text)
       if(text(i1:i1).eq.':') then
! condition given as an index in the condition list terminated by : like "1:"
          i1=1
          call getrel(text,i1,xxx)
          if(buperr.ne.0) then
             gx%bmperr=buperr; return
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
             gx%bmperr=4131; return
3310         continue
          else
             gx%bmperr=4131; return
          endif
! pcond points to condition record for axis, save in (map_axis) :: axarr
! check that it is not a fix phase condition (istv negative)
          if(pcond%statev.lt.0) then
             write(*,*)'Cannot set fix phase as axis'
             return
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
          stvr=>stvrvar
! this call also accept state variable functions like t_c, cp (if entered)
! UNIFINISHED: but it also accepts unknown texts ... 
          call decode_state_variable(text,stvr,ceq)
          if(gx%bmperr.ne.0) return
          pcond=>ceq%lastcondition
          i1=1; coeffs(1)=one
          call get_condition(i1,stvr,pcond)
          if(gx%bmperr.ne.0) then
! if new conditions are allowed then maybe enter this as condition
             write(*,*)'You must set the variable as a condition',&
                  ' before setting it as axis'
             return
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
          axarr(iax)%seqz=pcond%seqz
          axarr(iax)%more=0
! remove axplotdef for all axis!!! 
          axplotdef=' '
       endif
    endif removeaxis
    once=.TRUE.
3570 continue
    call gparrdx('Minimal value:',cline,last,xxx,dmin,'?Set axis')
    if(buperr.ne.0) return
    axarr(iax)%axmin=xxx
    call gparrdx('Maximal value:',cline,last,xxx,dmax,'?Set axis')
    if(buperr.ne.0) return
    if(xxx.le.axarr(iax)%axmin) then
       write(kou,*)'Maximal value must be higher than minimal'
       if(once) then
          once=.FALSE.
          goto 3570
       else
          write(kou,*)'Return to command level'
          return
       endif
    endif
    axarr(iax)%axmax=xxx
! default step 1/100 of difference ?? several diagram failed ...
! default step 1/40 of difference, same as TC ...
    dinc=0.025*(axarr(iax)%axmax-axarr(iax)%axmin)
! default step 1/50 of difference, somethimes better, sometimes worse ...
    call gparrdx('Increment:',cline,last,xxx,dinc,'?Set axis')
    if(buperr.ne.0) return
    if(xxx.lt.0.01*dinc) then
! someone (me) set xxx=0 and got a lot of trouble ...
       write(*,*)'Too small increment not allowed.'
       xxx=0.01*dinc
    endif
    axarr(iax)%axinc=xxx
! iax can be smaller than noofaxis if an existing axis has been changed
    if(iax.gt.noofaxis) noofaxis=iax
    return
  end subroutine ui_set_axis

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_bit_equilibrium(cline,last,ceq)
! SET BIT EQUILIBRIUM
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character ch1*1
    integer ll
!        EQNOTHREAD=0, EQNOGLOB=1, EQNOEQCAL=2,  EQINCON=3, &
3610 continue
    write(kou,3612)ceq%status
    call gparidx('Which bit? ',cline,last,ll,-1,'?Set status bit')
    if(cline(1:1).eq.'?') then
       write(kou,3612)ceq%status
3612   format('Set/reset bits of the equilibrium status word,',/&
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
       call gparcdx('Do you want to set the bit?',cline,last,1,&
            ch1,'Y','?Set status bit')
       if(ch1.eq.'Y') then
          ceq%status=ibset(ceq%status,ll)
          write(kou,3614)'set',ceq%status
3614      format('Bit ',a,', new equilibrium status word: ',z8)
       else
          ceq%status=ibclr(ceq%status,ll)
          write(kou,3614)'cleared',ceq%status
       endif
    endif
    return
  end subroutine ui_set_bit_equilibrium

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_bit_global(cline,last)
! SET BIT GLOBAL
    character cline*(*)
    integer last
!\end{verbatim}
    character ch1*1
    integer ll
3708 continue
! subroutine TOPHLP forces return with ? in position cline(last:last)
    write(kou,3709)globaldata%status
3709 format('Current global status word (hexadecimal): ',z8)
    call gparidx('Set/reset global status bit (from 0-31, -1 quits):',&
         cline,last,ll,-1,'?Global status bits')
    if(cline(1:1).eq.'?') then
       write(kou,3710)
3710   format('Set/reset bits of global status word ',&
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
       call gparcdx('Do you want to set the bit?',cline,last,1,&
            ch1,'Y','?Global status bits')
       if(ch1.eq.'Y') then
          globaldata%status=ibset(globaldata%status,ll)
          write(kou,3617)ll,' set',globaldata%status
3617      format('Bit ',i2,a,', new equilibrium status word: ',z8)
       else
          globaldata%status=ibclr(globaldata%status,ll)
          write(kou,3617)ll,' cleared',globaldata%status
       endif
! replaced by question above
       if(.not.btest(globaldata%status,GSADV)) then
! if expert/experienced bit is cleared ensure that occational user bit is set
          globaldata%status=ibset(globaldata%status,GSOCC)
       endif
    else
       write(kou,*)'Cannot be changed unless you have expert status'
    endif
    return
  end subroutine ui_set_bit_global

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_optcoeff_variable(cline,last)
! SET OPTCOEFF_VARIABLE, one coefficient or a range
    character cline*(*)
    integer last
!\end{verbatim}
    integer i1,i2,j4
    double precision xxx,xxy
    if(.not.btest(firstash%status,AHCOEF)) then
       write(kou,*)'No optimizing coefficients'
       return
    endif
! zero the relative standard deviation
    firstash%coeffrsd=zero
    call gparix('Coefficent index/range: ',cline,last,i1,-1,&
         '?Set variable coeff')
    if(i1.lt.0 .or. i1.ge.size(firstash%coeffstate)) then
! coefficients have indices 0 to size(firstash%coeffstate)-1
       write(kou,*)'No such coefficient'
       return
    endif
! upper limit must be negative and must follow directly on same line
    if(last.lt.len(cline) .and. cline(last:last).eq.'-') then
! pick up upper range limit as a negative value, 
! the question should thus never be asked ...
       last=last-1
       call gparix('Upper index (as negative): ',cline,last,i2,-i1,&
            '?Set variable coeff')
       if(i2.lt.0) then
! a negative value, its positive value must be >=i1
          i2=-i2
          if(i2.lt.i1) then
             i2=i1
             write(kou,*)'Illegal range, setting variable just: ',i1
          endif
! coefficients have indices 0 to size(firstash%coeffstate)-1
       else
! any other value ignored
          i2=i1
          write(kou,*)'Not understood, setting variable just: ',i1
       endif
    else
       i2=i1
    endif
! possible loop if i2>i1
    j4=i1
3740 continue
    xxy=firstash%coeffvalues(j4)*firstash%coeffscale(j4)
! this coefficient is not used, igore unless i1=i2
    if(i2.gt.i1 .and. firstash%coeffstate(j4).eq.0) goto 3745
    if(firstash%coeffstate(j4).lt.10) then
       nvcoeff=nvcoeff+1
    endif
    firstash%coeffstate(j4)=10
    if(i1.eq.i2) then
! when setting a single coefficient variable ask for value
       call gparrdx('Start value: ',cline,last,xxx,xxy,&
            '?Set variable coeff')
       if(buperr.ne.0) return
! set new value
       call change_optcoeff(firstash%coeffindex(j4),xxx)
       if(gx%bmperr.ne.0) return
       firstash%coeffvalues(j4)=one
       firstash%coeffscale(j4)=xxx
       firstash%coeffstart(j4)=xxx
    else
! coefficient used, set it variable with current value
       xxx=xxy
    endif
3745 if(i2.gt.j4) then
       j4=j4+1
       goto 3740
    endif
    write(kou,*)'Number of variable coefficients are ',nvcoeff
    return
  end subroutine ui_set_optcoeff_variable

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_range_exper_equil(cline,last)
! SET RANGE_EXPER_EQUIL
    character cline*(*)
    integer last
!\end{verbatim}
    integer i1,i2,j4
    if(allocated(firstash%eqlista)) then
       write(kou,*)'Experimental equilibria already entered'
       return
    endif
    call gparidx('First equilibrium number: ',cline,last,i1,2,&
         '?Set range')
    j4=noeq()
    call gparidx('Last equilibrium number: ',cline,last,i2,j4,&
         '?Set range')
    if(i2.lt.i1) then
       write(kou,*)'No equilibria?'
       return
    endif
! allocate the firstash%eqlista array and store equilibrium numbers
    j4=i2-i1+1
    firstash%firstexpeq=i1
    write(*,*)'Allocating firstash%eqlista ',j4,i1
    allocate(firstash%eqlista(j4))
    do i2=1,j4
       firstash%eqlista(i2)%p1=>eqlista(i1)
       i1=i1+1
    enddo
! close the plotdataunits!
    do i1=1,9
       if(plotdataunit(i1).gt.0) then
          write(plotdataunit(i1),22)graphopt%plotend
22        format('e'/a)
          close(plotdataunit(i1))
          plotdataunit(i1)=0
       endif
    enddo
    return
  end subroutine ui_set_range_exper_equil

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_set_optcoeff_fixed(cline,last)
! SET OPTCOEFF_FIXED, one coefficient or a range
    character cline*(*)
    integer last
!\end{verbatim}
    integer i1,i2,j4
    double precision xxx,xxy
! check not needed?
    if(.not.btest(firstash%status,AHCOEF)) then
       write(kou,*)'No optimizing coefficients'
       return
    endif
! zero the relative standard deviation
    firstash%coeffrsd=zero
! lower limit or range
    call gparix('Coeffient index/range: ',cline,last,i1,-1,&
         '?Set fix coeff')
    if(i1.lt.0 .or. i1.ge.size(firstash%coeffstate)) then
! coefficients have indices 0 to size(firstash%coeffstate)-1
       write(kou,*)'No such coefficient'
       return
    endif
! allow writing range on same line as 5-7 but also as 5 -7 on separate lines
    frange: if(last.lt.len(cline) .and. cline(last:last).eq.'-') then
       last=last-1
! upper limit must be negative
       call gparix('Upper index limit (as negative): ',&
            cline,last,i2,-i1,'?Set fix coeff')
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
    j4=i1
3720 continue
    xxy=firstash%coeffvalues(j4)*firstash%coeffscale(j4)
    if(i1.eq.i2) then
! A single coefficient, when fixing a single coefficinet ask for value
       call gparrdx('Fix value: ',cline,last,xxx,xxy,&
            '?Set fix coeff')
       if(buperr.ne.0) return
! set new value
       call change_optcoeff(firstash%coeffindex(j4),xxx)
       if(gx%bmperr.ne.0) return
       firstash%coeffvalues(j4)=one
       firstash%coeffscale(j4)=xxx
       firstash%coeffstart(j4)=xxx
    else
       call get_value_of_constant_index(firstash%coeffindex(j4),xxx)
    endif
! set as fixed without changing any min/max values (first time)
    if(firstash%coeffstate(j4).gt.13) then
       write(kou,*)'Coefficient state wrong, set to 1'
       firstash%coeffstate(j4)=1
       nvcoeff=nvcoeff-1
    elseif(firstash%coeffstate(j4).ge.10) then
       firstash%coeffstate(j4)=max(1,firstash%coeffstate(j4)-10)
       nvcoeff=nvcoeff-1
    elseif(xxx.ne.zero) then
! mark that the coefficient is fixed and nonzero 
       firstash%coeffstate(j4)=1
    else
! Fixed zero parameters are not listed
       firstash%coeffstate(j4)=-1
    endif
    if(i2.gt.j4) then
       j4=j4+1
       goto 3720
    endif
    write(kou,3730)nvcoeff
3730 format('Number of variable coefficients are now ',i3)
    return
  end subroutine ui_set_optcoeff_fixed

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_enter_tpfun(cline,last)
! ENTER TPFUN_SYMBOL, functions, constants (tables not implemented)
    character cline*(*)
    integer last
!\end{verbatim}
    character name1*24,name2*24,funstring*1024
    integer idef,jp,lrot
    double precision xxx,xxy
    call gparcx('TPfun name: ',cline,last,1,name1,' ','?Enter TPfun')
    if(buperr.ne.0) return
    if(.not.proper_symbol_name(name1,0)) then
       write(kou,*)'Bad symbol name'
       return
    endif
! check if already entered, 
    call find_tpsymbol(name1,idef,xxx)
    if(gx%bmperr.ne.0) then
! new symbol, can be function, constant or table (??)
       gx%bmperr=0
       call gparcdx('Function, constant or table? ',cline,last,1,name2,&
            'FUNCTION ','?Enter TPfun')
       if(buperr.ne.0) return
       call capson(name2)
       if(compare_abbrev(name2,'FUNCTION ')) then
! this call just read the function
          call enter_tpfun_interactivly(cline,last,funstring,jp)
          if(gx%bmperr.ne.0) return
! here the function is stored
          lrot=0
! last argument -1 means not reading from TDB file
          call store_tpfun(name1,funstring,lrot,-1)
          if(gx%bmperr.ne.0) return
       elseif(compare_abbrev(name2,'CONSTANT ')) then
! Enter a numeric constant
          call gparrdx('Value: ',cline,last,xxx,zero,'?Enter TPfun')
          call store_tpconstant(name1,xxx)
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
          call gparrdx('New value: ',cline,last,xxy,xxx,'?Enter TPfun')
          if(buperr.ne.0) return
          call capson(name1)
          call store_tpconstant(name1,xxy)
! we must evaluate all TPFUNS in all equilibria to be sure value propagates!
       endif
    endif
    return
  end subroutine ui_enter_tpfun

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_enter_optimize_coeff(cline,last)
! ENTER OPTIMIZE_COEFF, allocate optimizing coefficients A00 ...
    character cline*(*)
    integer last
!\end{verbatim}
    character name1*24
    integer i1,i2,j4
    if(.not.allocated(firstash%coeffstate)) then
       call gparidx('Number of coefficients: ',cline,last,i1,100,&
            '?Enter coeffs')
       if(buperr.ne.0) return
       i1=i1-1
       if(i1.lt.1) then
          write(*,*)'You must have at least 1 coefficient'
          return
       elseif(i1.gt.99) then
          write(*,*)'You cannot have more than 100 coefficient'
          return
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
       call enter_optvars(j4)
       call makeoptvname(name1,i1)
       write(kou,556)name1(1:3),i1
556    format(/'Coefficients entered with symbols A00 to ',a/&
            'Note that indices are from 0 to ',i2)
       do i2=0,i1
          firstash%coeffindex(i2)=j4+i2
       enddo
       firstash%status=ibset(firstash%status,AHCOEF)
    else
       write(kou,553)size(firstash%coeffstate)
553    format('You have already ',i3,' optimizing coefficients entered')
    endif
    call gparidx('Size of workspace: ',cline,last,lwam,2500,&
         '?Enter coeffs')
    if(allocated(wam)) then
       deallocate(wam)
       deallocate(iwam)
    endif
    return
  end subroutine ui_enter_optimize_coeff

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_enter_gnuplot_terminal(cline,last)
! ENTER GNUPLOT_TERMINAL, modify or add a GNUPLOT terminal
    character cline*(*)
    integer last
!\end{verbatim}
    character text*72,string*256,ch1*1
    integer i1,i2
    write(kou,172)graphopt%gnutermax
172 format('GNUPLOT terminals are:',i2/&
         4x,'Name',5x,'> command',6x,'GNUPLOT options')
    write(kou,173)(i2,graphopt%gnutermid(i2),&
         trim(graphopt%gnuterminal(i2)),i2=1,graphopt%gnutermax)
173 format(i2,2x,a,' > set terminal ',a)
    write(kou,174)
174 format('Change (exact match required) or enter a new GNUPLOT termial')
    call gparcx('Terminal id (8 chars):',cline,last,1,text,' ',&
         '?Enter GNUTERM')
    call capson(text)
    if(text(1:1).eq.' ') return
    do i1=1,graphopt%gnutermax
       if(text(1:8).eq.graphopt%gnutermid(i1)) then
          string=graphopt%gnuterminal(i1)
          write(*,*)'Modifying terminal ',graphopt%gnutermid(i1)
          goto 176
       endif
    enddo
! gnutermid not found, a new terminal
    call gparcdx('You want to enter a new terminal "'//trim(text)//'"?',&
         cline,last,1,ch1,'Y','?Enter GNUTERM')
    if(ch1.ne.'Y') then
       write(*,*)'Please try again'; return
    endif
    if(graphopt%gnutermax.ge.8) then
       write(kou,*)'There can max be 8 terminals'
       return
    endif
    i1=graphopt%gnutermax+1
    graphopt%gnutermax=i1
    string=' '
! enter a new set terminal id and definition
176 continue
    graphopt%gnutermid(i1)=text(1:8)
    call gparcx('Text after set terminal (see GNUPLOT manual):',&
         cline,last,5,text,string,'?Enter GNUTERM')
    graphopt%gnuterminal(i1)=text
    if(i1.ne.1) then
! SCREEN has no file extention
       call gparcx('File extention:',cline,last,1,text,' ',&
            '?Enter GNUTERM')
       graphopt%filext(i1)=text(1:4)
    endif
    write(*,179)i1,graphopt%gnutermid(i1),trim(graphopt%gnuterminal(i1)),&
         trim(graphopt%filext(i1))
179 format('New terminal definition for plot '/&
         i2,2x,a,'set terminal ',a/4x,'with file extention: ',a)
    return
  end subroutine ui_enter_gnuplot_terminal

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_read_tdb(cline,last,kom2,tdbfile,ocbase)
! READ TDB (kom2=2), READ ENCRYPTED (kom2=7) and READ SELECTED_PHASES (kom2=6)
! tdbfile is the previous TDB file (default), ocbase the database directory
    character cline*(*)
    integer last
    integer kom2
    character tdbfile*(*),ocbase*(*)
!\end{verbatim}
    character name1*24,text*72,line*80,selection*27,ch1*1
    character ellist(maxel)*2,elbase(maxel)*2
    integer i1,iel,j2,j4,j5,jp,kl,ll,ztyp
    logical selphases
! selphases is TRUE for READ SELECTED_PHASES
    selphases=kom2.eq.6
    if(selphases) then
! Ask if TDB or XTDB
       call gparcdx('Database format: ',&
            cline,last,1,name1,'TDB','?Read select phase')
       call capson(name1)
! here XTDB files are excluded temporarily
       if(name1(1:1).ne.'T') then
          write(*,*)'Selected phases implemented only for TDB files'
          gx%bmperr=4399; return
       endif
! indicate if the database is encrypted! (not changed for selected phases)
    elseif(kom2.eq.7) then
       globaldata%encrypted=1
    else
       globaldata%encrypted=0
    endif
    if(tdbfile(1:1).ne.' ') then
! set previous tdbfil as default
       text=tdbfile
       call gparcdx('File name: ',cline,last,1,tdbfile,text,'?Read TDB')
    else
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
       ztyp=1
       call gparfilex('File name: ',cline,last,1,tdbfile,' ',ztyp,&
            '?Read TDB')
    endif
! if tdbfle starts with "ocbase/" replace that with content of ocbase!!
! check for replacement of OCBASE probably redundant now ...
    name1=tdbfile(1:7)
    call capson(name1)
    if(name1(1:7).eq.'OCBASE/' .or. name1(1:8).eq.'OCBASE\ ') then
       tdbfile=trim(ocbase)//tdbfile(7:)
       write(*,*)'database file: ',trim(tdbfile)
    endif
! this call checks the file exists and returns the elements
! it also lists the DATABASE_INFO text
    call checkdb2(tdbfile,'.tdb',jp,elbase)
    if(gx%bmperr.ne.0) then
       write(kou,*)'No database with this name'
       tdbfile=' '
       return
    elseif(jp.eq.0) then
       write(kou,*)'No elements in the database'
       tdbfile=' '
       return
    elseif(jp.lt.0 .and. .not.selphases) then
! encrypted databases return jp=-1, we do not know number of elements ...
       write(kou,*)'Cannot list elements in encrypted databases'
       j4=20
       goto 8207
    endif
    j4=jp
    write(kou,8203)jp,(elbase(kl),kl=1,j4)
8203 format('Database has ',i2,' elements: ',18(a,1x)/(1x,28(1x,a)))
    ellist='  '
    write(kou,8205)
8205 format('Give the elements to select, finish with empty line')
8207 continue
    jp=1
    selection='Select elements /all/:'
8210 continue
    call gparcx(selection,cline,last,1,ellist(jp),' ','?Read TDB')
    if(jp.eq.1 .and. cline(1:4).eq.'all ') then
! this is if someone actually types "all".  If he types "ALL" that will be AL
       jp=0
    elseif(cline(1:1).eq.'q' .or. cline(1:1).eq.'Q' .or.&
         cline(1:4).eq.'NONE') then
! if user regets selection he can quit
       write(*,*)'Quitting, nothing selected'
       return
    elseif(ellist(jp).ne.'  ') then
       call capson(ellist(jp))
       jp=jp+1
       if(jp.gt.size(ellist)) then
          write(kou,*)'Max number of elements selected: ',size(ellist)
       else
          ll=last
! Check if element exist, unless encrypted ...
          if(selphases .or. globaldata%encrypted.eq.0) then
             elcheck: do j5=1,j4
                if(ellist(jp-1).eq.elbase(j5)) exit elcheck
             enddo elcheck
! if we come here with j4>j5 then ellist(jp) is not in elbase(1..j4)
             if(j5.gt.j4) then
                jp=jp-1
                if(selphases) then
                   write(kou,'(a,i3,1x,a)')'No such element: ',jp,ellist(jp)
                else
                   write(kou,'(a,i3,1x,a)')' *** WARNING: No such element:',&
                        jp,ellist(jp)
                endif
             endif
          endif
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
8220   format('Selected ',i2,' elements: ',20(a,1x))
    endif
    if(selphases) then
! SPECIAL SELECT_PHASES
! ask for phses to be selected, max 100, seltdbph global variable
       allocate(seltdbph(100))
       j4=0
       selection='Select phase(s) /all/:'
       selph: do while (.TRUE.)
          call gparcdx(selection,&
               cline,last,1,line,' ','?Read select phase')
          if(line(1:1).eq.' ') exit selph
          selection='Select more phase(s):'
! There is at least one phase name in line
          j2=1
          phname: do while(.not.eolch(line,j2))
             j4=j4+1
             if(j4.gt.100) then
                write(*,*)'Max 100 phases can be selected'
                exit selph
             endif
             j2=j2-1
! getext increments i2 by 1 at each call.  A space or , between each name
             call getext(line,j2,1,seltdbph(j4),' ',i1)
             call capson(seltdbph(j4))
          enddo phname
       enddo selph
! nselph is a global variabel
       nselph=j4
       if(nselph.gt.0) then
          write(*,*)nselph,' phase abbreviations specified '
       else
          write(*,*)'No phase specified, all will be read'
       endif
! Finally read the TDB file
! if seltdbph is allocated only those phase will be inlcuded
       call readtdb(tdbfile,jp,ellist)
       deallocate(seltdbph)
    else
       call readtdb(tdbfile,jp,ellist)
    endif
    if(gx%bmperr.ne.0) then
! inside readtdb any "buperr" will be set as gx%bmperr
       write(kou,*)'There were errors reading the TDB file', gx%bmperr
       if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
          write(kou,*)bmperrmess(gx%bmperr)
       endif
       write(kou,*)'Please correct these before continuing'
! ignore any type ahead
       last=len(cline)
       call gparcdx('Do you want to continue anyway?',&
            cline,last,1,ch1,'N','?Read TDB error')
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
8221   format(/' *** There were warnings from reading the database'/&
            ' *** If you run a macro file please scroll back and check!'/)
    endif
    return
  end subroutine ui_read_tdb

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_read_xtdb(cline,last,xtdbfile)
! READ XTDB, xtdbfile is the previous XTDB file (default)
    character cline*(*)
    integer last
    character xtdbfile*(*)
!\end{verbatim}
    character name1*24,text*72,selection*27
    character ellist(maxel)*2
    integer iel,jp,kl,ztyp
    if(xtdbfile(1:1).ne.' ') then
       text=xtdbfile
       write(*,*)'debug; ',trim(text)
       call gparcdx('File name: ',cline,last,1,xtdbfile,text,'?Read XTDB')
    else
! THESE TYPES ARE USED ALSO IN METLIB4
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT, 8=LOG
! negative is for write, 0 read without filter, -100 write without filter
       ztyp=6
       call gparfilex('File name: ',cline,last,1,xtdbfile,' ',ztyp,&
            '?Read XTDB')
    endif
    if(xtdbfile(1:1).eq.' ') return
! this call checks the file exists and returns the elements
! It is in gtp3EY and can handle the <Element keyword in XTDB files
    jp=0
    write(*,*)'Opening: ',trim(xtdbfile)
    call xtdbread(xtdbfile,jp,ellist)
    if(gx%bmperr.ne.0) then
       write(kou,*)'No XTDB database with this name'
       return
    elseif(jp.eq.0) then
       write(Kou,*)'No elements in the database'
       return
    endif
    name1=ellist(1)
    write(kou,8203)jp,(ellist(kl),kl=1,jp)
    write(kou,8205)
    kl=jp
    jp=1
    selection='Select elements /all/:'
8217 continue
    call gparcx(selection,cline,last,1,ellist(jp),' ','?Read XTDB')
    if(cline(1:4).eq.'NONE') then
! if user regets selection he can quit
       write(*,*)'Quitting, nothing selected'
       return
    endif
    if(ellist(jp).ne.'  ') then
       call capson(ellist(jp))
       jp=jp+1
       if(jp.gt.kl) then
       write(kou,*)'Max number of elements selected: ',jp,kl,size(ellist)
          selection='Select elements /all/:'
          goto 8217
       else
          selection='Select elements /no more/:'
          goto 8217
       endif
    elseif(index(selection,'all').gt.0) then
! user has selected all, restore ellist(1)
       ellist(1)=name1
       jp=jp-1
    endif
8212 continue
    if(jp.eq.0) then
       jp=kl
       write(kou,*)'All elements selected',jp
    else
       write(*,8220)jp,(ellist(iel),iel=1,jp)
    endif
! This should read the XTDB files in new XML format.  This is in gtp3EX/Y.F90
    call xtdbread(xtdbfile,jp,ellist)
! NOT YET WRITTEN transfer selected system from gtp3_xml.F90 to OC proper
! also list the bibliography
    write(kou,*)
8203 format('Database has ',i2,' elements: ',18(a,1x)/(1x,28(1x,a)))
8205 format('Give the elements to select, finish with empty line')
8220 format('Selected ',i2,' elements: ',20(a,1x))
    return
  end subroutine ui_read_xtdb

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_save_unformatted(cline,last,ocufile,model)
! SAVE UNFORMATTED, ocufile is the previous file (default), model has the date
    character cline*(*)
    integer last
    character ocufile*(*),model*(*)
!\end{verbatim}
    character text*72,ch1*1
    integer ieq,jp,kl,ztyp
    logical logok
132 continue
! save unformatted after step/map not recommended as equilibria
! unless equilibria with _MAPLINE and _MAPNODE not deleted
! Reading an unformatted file with these prevents any new new STEP/MAP
    call findeq('_MAPLINE_1 ',ieq)
    if(gx%bmperr.eq.0) then
       write(kou,*)'Please use DELETE STEP_MAP before unformatted save'
       return
    else
! there are no map/step equilibria, OK to save
       gx%bmperr=0
    endif
!
    if(ocufile(1:1).ne.' ') then
       text=ocufile
       call gparcdx('File name: ',cline,last,1,ocufile,text,&
            '?Save unformatted')
    else
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
       ztyp=-2
       call gparfilex('File name: ',cline,last,1,ocufile,' ',ztyp,&
            '?Save unformatted')
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
       return
    endif
! I have no way to handle the extention to upper case ... inside C routine
    if(jp.gt.0) ocufile(jp+1:)='.OCU '
    inquire(file=ocufile,exist=logok)
    if(logok) then
       call gparcdx('File exists, overwrite?',cline,last,1,ch1,'N',&
            '?Save overwite')
       if(ch1.ne.'Y') then
          write(*,133)
133       format('Please use another file name')
          ocufile=' '
          goto 132
       endif
       write(*,134)trim(ocufile)
134    format(/'Overwriting previous results on ',a)
    endif
    text='U '//model
    call gtpsave(ocufile,text)
    return
  end subroutine ui_save_unformatted

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_list_state_variables(cline,last,kom,kom2,lut,ceq)
! LIST STATE_VARIABLES (kom2=4), LIST MODEL_PARAM_VAL (kom2=17) and SHOW (kom=25)
! lut is the output unit
    character cline*(*)
    integer last
    integer kom,kom2,lut
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,line*80,model*72,longstring*5000,actual_arg(2)*16
    integer i1,i2,istv,j4,ll,mode
    double precision xxx
    double precision, dimension(maxconst) :: yarr
    logical once
    if(btest(ceq%status,EQNOEQCAL) .or. btest(ceq%status,EQFAIL)) then
       write(lut,6101)
6101   format(' *** Warning,',&
          'equilibrium not calculated, values are probably wrong')
    elseif(btest(ceq%status,EQINCON)) then
       write(lut,6102)
6102   format(' *** Warning, values can be inconsistent with',&
          ' current conditions')
    endif
    once=.TRUE.
! LOOP here for list state_variables or model_parameter_values or SHOW
6105 continue
! NOTE: 4th argument is 5 because otherwise a "," will terminate reading cline
! and state variables like x(fcc,cr) will not work.
    if(kom.eq.25) then
! SHOW: this execute the SHOW command
       call gparcx('Property: ',cline,last,5,line,' ','?Show property')
    else
! the command is LIST STATE_VARIABLES
       if(kom2.eq.4) then
          call gparcx('State variable: ',cline,last,5,line,' ',&
               '?List state variables')
       else
! the command is LIST MODEL_PARAMETER_VALUE                
          if(once) then
             write(kou,*)'Remember always to specify the phase!'
             once=.FALSE.
          endif
          call gparcx('Parameter ident: ',cline,last,5,line,' ',&
               '?List model parameter val')
       endif
    endif
! if line empty return to command level
    j4=1
    if(eolch(line,j4)) return
    j4=index(line,',')
    if(j4.gt.0) then
! check if there is a , before a ( as that is not allowed.  There are
! state variables like x(fcc,cr) ... (this is not a strong test ...)
       ll=index(line,'(')
       if(ll.le.0 .or. ll.gt.j4) then
          write(*,*)'Please use a space as separator',&
               ' except within ( ) as in x(liq,cr) !'
          return
       endif
    endif
! model is just used to return texts
    model=' '
! we should extract the text from last up to first space and save rest in cline
    j4=index(line,' ')
    name1=line(1:j4)
    call capson(name1)
! dot derivatives not allowed explicitly, must be entered as symbols
    if(index(name1,'.').gt.0) then
       write(kou,*)'You must enter dot derivatives as symbols!'
       return
    endif
! note gparc etc increment last before looking for answer, keep space in cline
    cline=line(j4:)
    last=1
! allow also DGM(#) to generate many values ...
    if(index(name1,'*').gt.0 .or. index(name1,'DGM(#)').gt.0) then
! generate many values
! i1 values are resturned in yarr with dimension maxconst. 
! longstring are the state variable symbols for the values ...
       call get_many_svar(name1,yarr,maxconst,i1,longstring,ceq)
       if(gx%bmperr.eq.0) then
! not a nice output, we should include the state variables FIX!!
          write(lut,6106)i1,longstring(1:len_trim(longstring))
6106      Format('Listing of ',i3,' state variables:'/a)
          write(lut,6107)(yarr(i2),i2=1,i1)
6107      format('Values: ',5(1pe14.6)/(8x,5(1pe14.6)))
          if(index(name1,'*,').gt.0) write(*,6121)trim(name1)
6121      format(' *** Note that for unstable phases ',a,&
               ' is not shown or listed as zero')
       endif
    else
! the value of a state variable, symbol? or model parameter variable is returned
! STRANGE the symbol xliqni is accepted in get_state_var_value ???
! get_state_var_value is in gtp3F.F90
       call get_state_var_value(name1,xxx,model,ceq)
       if(gx%bmperr.eq.0) then
          write(lut,6108)trim(model),xxx
6108      format(1x,a,'=',1PE15.7)
       else
          gx%bmperr=0
! If error then try to calculate a symbol ...
! below copied from calculate symbol, first calculate all symbols ignore errors
! calculate all symbols ignoring errors (note dot derivatives not calculated)
          call meq_evaluate_all_svfun(-1,ceq)
          if(gx%bmperr.ne.0) gx%bmperr=0
          call capson(line)
          call find_svfun(name1,istv)
          if(gx%bmperr.ne.0) return
          mode=1
          actual_arg=' '
          xxx=meq_evaluate_svfun(istv,actual_arg,mode,ceq)
          if(gx%bmperr.ne.0) return
          write(kou,2047)trim(name1),xxx
! this format statement elsewhere
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
2047 format(a,'= ',1pe16.8)
    return
  end subroutine ui_list_state_variables

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_list_equilibria(kom2,lut,ceq)
! LIST EQUILIBRIA (kom2=11) and LIST ACTIVE_EQUILIBRIA (kom2=19)
! lut is the output unit
    integer kom2,lut
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,text*72
    integer i2,iel,iz,jp,jz,nv
! if 19 then skip equilibria with zero weight
    nv=noeq()
! skip if there is just one equilibrium kom=6=LIST; kom2=19=ACTIVE-EQUIL
    if(kom2.eq.19 .and. nv.eq.1) return
    write(lut,6212)
6212 format('Number  Name',25x,'T   Weight Comment & phases')
    jp=0
    do iel=1,nv
       if(associated(ceq,eqlista(iel))) then
          name1='**'
       else
          name1=' '
       endif
       text=eqlista(iel)%comment
       jz=len_trim(text)
       if(jz.lt.20) then
! if there is space add names of stable phases
          if(jz.gt.0) then
             text(jz+1:)=' & '; jz=jz+4
          else
             jz=1
          endif
          do iz=1,nooftup()
             i2=phasetuple(iz)%lokvares
             if(eqlista(iel)%phase_varres(i2)%phstate.gt.0) then
                if(eqlista(iel)%phase_varres(i2)%phstate.eq.2) then
! prefix any FIX phase with *
                   text(jz:jz)='*'; jz=jz+1
                endif
                call get_phasetup_name(iz,text(jz:))
! text is limited to 72 characters and anyway only 32 are written
                jz=min(len_trim(text)+2,40)
             endif
          enddo
       endif
       if(eqlista(iel)%weight.gt.zero) then
! always list equilibria with weight>0
          write(lut,6203)iel,name1(1:2),eqlista(iel)%eqname,&
               eqlista(iel)%tpval(1),eqlista(iel)%weight,trim(text)
6203      format(i4,1x,a2,1x,a,1x,F8.2,1x,F5.2,1x,a)
       elseif(iel.eq.1 .or. kom2.eq.11) then
! for kom2=11 list all equilibria without including weight
! NOTE all equilibria outside "range" (default and step/map) has weight= -1.0
          write(lut,6202)iel,name1(1:2),eqlista(iel)%eqname,&
               eqlista(iel)%tpval(1),trim(text)
6202      format(i4,1x,a2,1x,a,1x,F8.2,7x,a)
       elseif(eqlista(iel)%weight.eq.zero) then
          jp=jp+1
       endif
    enddo
    if(kom2.eq.19 .and. jp.gt.0) &
         write(lut,'(/"Number of equilibria with zero weight: ",i4)')jp
    return
  end subroutine ui_list_equilibria

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_list_results(cline,last,version,listresopt,lrodef,lut,ceq)
! LIST RESULTS, listresopt is the output mode, lrodef its default
! lut is the output unit
    character cline*(*)
    integer last
    character version*(*)
    integer listresopt,lrodef,lut
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,optres*40
    integer ics,iph,j4,mode
    logical once
! skip if no calculation made
    if(btest(globaldata%status,GSNOPHASE)) then
       write(kou,*)'No results as no data'
       return
    elseif(btest(ceq%status,EQGRIDCAL)) then
       write(kou,*)' *** Last calculation was not a full equilibrium'
    endif
    call gparidx('Results output mode: ',cline,last,&
         listresopt,lrodef,'?List results')
    if(buperr.ne.0) then
       write(kou,*)'No such mode, using default'
       buperr=0
       listresopt=lrodef
    endif
! CCI extending the number of listing options
    if(listresopt.gt.0 .and. listresopt.le.12) then
       lrodef=listresopt
    endif
! CCI end          
    call date_and_time(optres,name1)
    write(lut,6051)version,ceq%eqno,ceq%eqname,&
         optres(1:4),optres(5:6),optres(7:8)
! write comment line if any
    if(len_trim(ceq%comment).gt.0) then
       write(lut,6308)trim(ceq%comment)
6308   format(3x,a)
    endif
    if(btest(ceq%status,EQFAIL)) then
       write(lut,6305)'below'
6305   format(/' *** The results ',a,&
            ' are not a valid equilibrium as last calculation failed'/)
    elseif(btest(ceq%status,EQNOEQCAL)) then
       write(lut,6307)'below'
6307   format(/' *** The results listed ',a,' does not represent',&
            ' a calculated equilibrium'/)
    elseif(btest(ceq%status,EQINCON)) then
       write(lut,6306)'below'
6306   format(/' *** The results listed ',a,' may be inconsistent',&
            ' with the current conditions'/)
    endif
    write(lut,6302)'Conditions .............................'
6302 format(a,20('.'),':')
6303 format(/a,20('.'),':')
    call list_conditions(lut,ceq)
    write(lut,6303)'Some global data, reference state SER ..'
    call list_global_results(lut,ceq)
    if(btest(ceq%status,EQNOEQCAL)) then
       write(*,6277)ceq%status
6277   format(' *** No results as no equilibrium calculated! ',z8)
       goto 6363
    endif
    write(lut,6303)'Some data for components ...............'
    j4=1
    if(listresopt.ge.4 .and. listresopt.le.7) then
! j4=2 means mass fractions
       j4=2
    endif
    call list_components_result(lut,j4,ceq)
! Phase output starts with newline
    write(lut,6304,advance='no')'Some data for phases ...................'
6304 format(/a,20('.'),':')
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
    elseif(listresopt.eq.12) then
! 12: just one line per phase, no compositions
       mode=10020
    else
! all phase with with mole fractions
       mode=0
    endif
    ics=1
    once=.TRUE.
    do iph=1,noph()
       ics=0
6310   continue
       ics=ics+1
! moved to gtp3C
! use phase amount in mass
!                     '    Form.U    At/FU     DGM   Frac:')
! use phase amount in mole
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
6363 continue
    if(associated(ceq%lastexperiment)) then
       write(lut,491)ceq%weight
491    format('Weight ',F6.2)
! list all experiments ........................................
       call meq_list_experiments(lut,ceq)
       write(lut,*)
    endif
    if(btest(ceq%status,EQNOEQCAL)) return
! list if anyting should be calculated or listed separately (not their values)
    if(allocated(ceq%eqextra)) then
       write(lut,492)ceq%eqextra(1)(1:len_trim(ceq%eqextra(1))),&
            ceq%eqextra(2)(1:len_trim(ceq%eqextra(2)))
492    format('Calculate: ',a/'List: ',a)
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
6051 format(/'OC version',a,' equilibrium: ',i3,', ',&
          a,3x,a4,'.',a2,'.',a2)
    return
  end subroutine ui_list_results

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_debug_symbol_value(cline,last,ceq)
! DEBUG SYMBOL_VALUE, test that a symbol has a value, used in test macros
! stops the program if not
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,line*80,model*72,actual_arg(2)*16
    integer istv,mode
    double precision xxx,xxy,xxz
! this allows a command "debug symbol value" which will test if symbol
! has the specified value (+/- 10^(-6).  
! Should be useful in the test macros ...
! 4th argument 2 means terminate at " ", ignore any ,
    call gparcx('Symbol: ',cline,last,2,name1,' ','?Debug symbol value')
    call gparrx('Value: ',cline,last,xxy,zero,'?Debug symbol value')
    call capson(name1)
! code below copied from SHOW command
    model=' '
    call get_state_var_value(name1,xxx,model,ceq)
    if(gx%bmperr.ne.0) then
       gx%bmperr=0
! If error then try to calculate a symbol ...
! below copied from calculate symbol, first calculate all symbols ignore errors
! calculate all symbols ignoring errors (note dot derivatives not calculated)
       call meq_evaluate_all_svfun(-1,ceq)
       if(gx%bmperr.ne.0) gx%bmperr=0
       call capson(line)
       call find_svfun(name1,istv)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error finding symbol'
          stop
       endif
       mode=1
       actual_arg=' '
       xxx=meq_evaluate_svfun(istv,actual_arg,mode,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error calculating symbol'
          stop
       endif
       write(kou,2047)trim(name1),xxx
    endif
! test for NaN, a NaN is never equal to itself
    if(xxx .ne. xxx) then
       write(kou,*)'Calculated value of ',trim(name1),' is NaN'
       stop
    endif
    xxz=1.0D-6
    if(abs(xxy).gt.1.0d0) xxz=xxz*abs(xxy)
    if(abs(xxx-xxy).gt.xxz) then
       write(kou,'(a,2(1pe12.4))')'Symbol value outside limit!',xxx,xxy
       stop
    else
       write(kou,*)'Testing symbol ',trim(name1),' value OK ++++++++'
    endif
2047 format(a,'= ',1pe16.8)
    return
  end subroutine ui_debug_symbol_value

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_amend_phase_magnetic(cline,last,lokph)
! AMEND PHASE <name> ADDITION MAGNETIC_CONTRIBUTION, lokph is the phase record
    character cline*(*)
    integer last
    integer lokph
!\end{verbatim}
    character aux*4,ch1*1,chz*1
    integer idef,j2,j4
    idef=-3
! zero value of antiferromagnetic factor means Inden-Qing model
    call gparidx('Antiferromagnetic factor: ',&
         cline,last,j4,idef,'?Amend magnetism')
    if(buperr.ne.0) return
    if(j4.eq.0) then
! Inden-Hillert-Qing-Xiong magnetic model has AFF=0
       call gparcdx('BCC type phase: ',cline,last,1,chz,'N',&
            '?Amend magnetism')
       call gparcdx('Using individual Bohr magnetons: ',&
            cline,last,1,ch1,'N','?Amend magnetism')
       if(.not.(ch1.eq.'Y' .or. ch1.eq.'y')) then
          call set_phase_status_bit(lokph,PHBMAV)
          aux=' '
       else
          write(*,*)'PMON mark use IBM parameter'
          aux(2:2)='I'
       endif
! xiongmagnetic is a predefined addition index, chz is Y or y for BCC
       j2=xiongmagnetic
       aux(1:1)=chz
       call add_addrecord(lokph,aux,xiongmagnetic)
    else
       if(j4.eq.-1) then
! Inden magnetic for BCC
          call add_addrecord(lokph,'Y',indenmagnetic)
       else
! Inden magnetic for FCC
          call add_addrecord(lokph,'N',indenmagnetic)
       endif
       j2=indenmagnetic
    endif
      call gparcdx('Is the addition calculated per mole of atoms?',&
         cline,last,1,ch1,'Y','?Add per formula unit')
! The magnetic model calculates a molar Gibbs energy, must be multiplied with
! the number of atoms in the phase. j2 set above to the addition type
    if(ch1.eq.'Y' .or. ch1.eq.'y') then
       call setpermolebit(lokph,j2)
    endif
    return
  end subroutine ui_amend_phase_magnetic

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_amend_phase_disordered(cline,last,iph)
! AMEND PHASE <name> DISORDERED_FRACS, iph is the phase index
    character cline*(*)
    integer last
    integer iph
!\end{verbatim}
    character ch1*1
    integer idef,j4,ll,lokcs,ndl
    double precision xxx
    if(.not.allowenter(2)) then
       gx%bmperr=4125
       return
    endif
! we should check the number of sublattices of the phase ...
    lokcs=phasetuple(iph)%lokvares
    idef=size(firsteq%phase_varres(lokcs)%sites)
    call gparidx('Sum up to sublattice: ',cline,last,ndl,idef,&
         '?Amend phase disordfrac')
    if(buperr.ne.0) return
    call gparcdx('Should the ordered part cancel when disordered? ',&
         cline,last,1,ch1,'N','?Amend phase disordfrac')
    if(buperr.ne.0) return
    if(ch1.eq.'N' .or. ch1.eq.'n') then
! like sigma which is never completely disordered
       j4=0
    else
! like FCC ordering where the disordered state can be modeled independently
       j4=1
       write(kou,*)'This phase can be totally disordered'
    endif
    ch1='D'
    call add_fraction_set(iph,ch1,ndl,j4)
! forgot to add the sites 
    lokcs=phasetuple(iph)%lokvares
    if(j4.eq.0) then
       xxx=zero
       do ll=1,ndl
          xxx=xxx+firsteq%phase_varres(lokcs)%sites(ll)
       enddo
       firsteq%phase_varres(lokcs)%disfra%fsites=xxx
    endif
    if(gx%bmperr.ne.0) return
    return
  end subroutine ui_amend_phase_disordered

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_calculate_tpfun(cline,last,lut,ceq)
! CALCULATE TPFUN_SYMBOLS at current T and P, lut is the output unit
    character cline*(*)
    integer last
    integer lut
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24
    integer iel,j4,lrot
    double precision val(6)
    logical once
    call gparcdx('name: ',cline,last,5,name1,'*','?Calculate TPfun')
    lrot=0
    iel=index(name1,'*')             
    if(iel.gt.1) name1(iel:)=' '
! as TP functions call each other force recalculation and calculate all
! even if just a single function is requested
    call change_optcoeff(-1,zero)
    do j4=1,notpf()
       call eval_tpfun(j4,ceq%tpval,val,ceq%eq_tpres)
       if(gx%bmperr.gt.0) return
    enddo
    if(name1(1:1).ne.'*') then
       once=.TRUE.
2009   continue
       call find_tpfun_by_name(name1,lrot)
       if(gx%bmperr.ne.0) then
          if(iel.eq.0) return
          gx%bmperr=0
       else
! found function number from lrot ???
          j4=lrot
          call eval_tpfun(j4,ceq%tpval,val,ceq%eq_tpres)
          if(gx%bmperr.gt.0) return
          if(once) then
             once=.FALSE.
             write(lut,2011)1,ceq%tpval
          endif
          write(lut,2012)j4,val
          if(iel.gt.1) goto 2009
       endif
    else
       write(lut,2011)notpf(),ceq%tpval
2011   format(/'Calculating ',i6,' functions for T,P=',F10.2,1PE15.7/&
            3x,'No   F',11x,'F.T',9x,'F.P',9x,'F.T.T',&
            7x,'F.T.P',7x,'F.P.P')
       do j4=1,notpf()
          call eval_tpfun(j4,ceq%tpval,val,ceq%eq_tpres)
          if(gx%bmperr.gt.0) return
          write(lut,2012)j4,val
2012      format(I5,1x,6(1PE12.4))
       enddo
    endif
    return
  end subroutine ui_calculate_tpfun

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_calculate_phase_adjust(cline,last,iph,ics,kom3,lut,ceq)
! CALCULATE PHASE <name> CONSTITUTION_ADJ (kom3=4) or DIFFUSION_COEFF (kom3=5)
! iph and ics are phase and composition set, lut is the output unit
    character cline*(*)
    integer last
    integer iph,ics,kom3,lut
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character elsym*2,quest*20
    integer jp,nend,nsub,nv
    double precision xxx,xxy,totam
    double precision, dimension(maxel) :: xknown,aphl
    double precision, dimension(maxconst) :: yarr
    double precision, allocatable, dimension(:) :: mugrad,mobilities,nsites
    integer, allocatable, dimension(:) :: nkl
    type(gtp_phasetuple), pointer :: phtup
! or derivatives of chemical potentials and mobility data
! convert to phase tuple here as that is used in the application call
    do jp=1,nooftup()
       if(phasetuple(jp)%ixphase.eq.iph .and. &
            phasetuple(jp)%compset.eq.ics) then
          goto 2044
       endif
    enddo
    write(*,*)'No such tuple'
    return
2044 continue
    phtup=>phasetuple(jp)
! Get current constitution of the phase
    call calc_phase_molmass(iph,ics,xknown,aphl,totam,xxy,xxx,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error finding current composition'
       return
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
          call gparrdx(quest,cline,last,xxy,xknown(nv),&
               '?Calculate phase adjust')
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
2088 format('Mole fraction of ',a,' set to ',F8.5)
    xknown(nv)=totam
    yarr(nv)=ceq%cmuval(nv)
! use current T and P
    if(kom3.eq.4) then
! constituent adjustment, the FALSE means not quiet
       call equilph1b(phtup,ceq%tpval,xknown,xxx,yarr,.FALSE.,ceq)
       if(gx%bmperr.ne.0) return
       write(kou,2087)xxx,(yarr(nv),nv=1,noel())
2087   format(/'Calculated Gibbs energy/FU/RT: ',1pe14.6,&
            ' and the chemical potentials/RT:'/6(1pe12.4))
    else
!.............................................
! calculate phase diffusion: chem.pot derivatives and mobilities
! mugrad(I,J) are derivatives of the chemical potential of endmember I
!         with respect to endmember J
! mobilities(i) is mobility of component i
! derivatives of mu and mobilities, the FALSE means not quiet
! CCI
       call get_sublattice_number(phtup%ixphase,nsub,ceq)
       allocate(nkl(nsub))
       allocate(nsites(nsub))
       call get_sublattice_structure(phtup%ixphase,phtup%compset,&
            nsub,nkl,nsites,ceq)
       nend=1
       do nv=1,nsub
         nend = nend * nkl(nv)
       enddo
       deallocate(nkl)
       deallocate(nsites)
       allocate(mugrad(nend*nend))
       allocate(mobilities(noel()))
       mugrad(:)=zero
       mobilities(:)=zero

       call equilph1d(phtup,ceq%tpval,xknown,yarr,.FALSE.,&
            nend,mugrad,mobilities,ceq)
       if(gx%bmperr.ne.0) return
       write(kou,2096)nend
2096   format(/'Chemical potential derivative matrix, dG_I/dn_J for ',&
            i3,' endmembers')
       write(kou,2094)(nv,nv=1,nend)
2094   format(3x,6(6x,i6)/(3x,6i12))
       do nv=0,nend-1
! An extra LF is generated when just 6 components!! use ll, kp j4, i2
          write(kou,2095)nv+1,(mugrad(nend*nv+jp),jp=1,nend)
2095      format(i3,6(1pe12.4)/(3x,6e12.4))
       enddo
       write(kou,2098)noel()
2098   format(/'Mobility values mols/m2/s ?? for',i3,' components')
       write(kou,2095)1,(mobilities(jp),jp=1,noel())
       !CCI
       deallocate(mugrad)
       deallocate(mobilities)
       call list_defined_properties(lut)
!CCI
    endif
    return
  end subroutine ui_calculate_phase_adjust

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_amend_phase_asymmetries(cline,last,iph,ceq)
! AMEND PHASE <name> ASYMMETRIES, change the asymmetry of an MQMQA ternary
! iph is the phase index
    character cline*(*)
    integer last
    integer iph
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character ch1*1
    integer asymter,new_toop,lokcs
    TYPE(gtp_phase_varres), pointer :: parres
! This is to change the asymmetry of a ternary
    if(.not.allocated(tersys)) then
! if there is an MQMQA phase then tersys should be allocated
       write(*,*)'No MQMQA phase asymmetries entered'
       return
    endif
    lokcs=phasetuple(iph)%lokvares
    write(*,*)'You can now try to amend MQMQA ternary'
! copied from gtp3XQ listconst
    parres=>ceq%phase_varres(lokcs)
    call list_quads_with_single_cation(parres)
! tersys is global data
    call list_tersys
! list all ternaries with quad names
    call list_ternary_cations(0,parres)
    write(kou,*)
    call gparidx('Index of ternary to modify (0=none)?',cline,last,&
         asymter,0,'?Asymmetry modify')
    if(asymter.le.0 .or. asymter.gt.size(tersys)) then
       write(*,3125)
3125   format(/'No such ternary, nothing done'/)
       return
    endif
    if(index(tersys(asymter)%asymm,'T').gt.0) then
       write(*,*)'ui_Amend: ternary is already asymmetric: '
       call gparcdx('Do you want to remove asymmetry?',cline,last,1,&
            ch1,'Y','?Asymmetry modify')
       if(ch1.eq.'Y' .or. ch1.eq.'y') then
! remove asymmetry
          tersys(asymter)%asymm='KKK'
          new_toop=0
          call new_ternary_asym(asymter,new_toop,parres,.false.)
          write(*,*)'Ternary set symmetric'
       else
          write(*,*)'Ternary asymmetry not changed, try again'
       endif
       return
    endif
! set isasym zero
    call list_ternary_cations(asymter,parres)
    write(*,414)
414 format('Use 1, 2 or 3 to specify the Toop cation in this ternary!')
    call gparidx('Use  0, 1, 2 or 3',cline,last,&
         new_toop,0,'?Asymmetry modify')
! change of asymmetry
    if(new_toop.le.0 .or. new_toop.gt.3) then
       write(kou,3127)
3127   format('No change of ternary asymmetry '/&
            'Use only numbers 1, 2 or 3 in the ternary cati order')
       asymter=0
       return
    endif
! expression for varkappa and other variables will be made automatically?
3122 continue
    newXupdate=newXupdate+1
! asymter is the index in the tersys array of the teranary with new asymmetry
    call gparcdx('Verbose?',cline,last,1,ch1,'N','?Asymmetry modify')
    if(ch1.eq.'Y') then
       call new_ternary_asym(asymter,new_toop,parres,.true.)
       write(*,*)'MM back from new_ternary_asym'
    else
       call new_ternary_asym(asymter,new_toop,parres,.false.)
    endif
! repeat short listing the asymmetries, format label 3101 below
    call list_tersys
    return
  end subroutine ui_amend_phase_asymmetries

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_calculate_all_equilibria(cline,last,lut)
! CALCULATE ALL_EQUILIBRIA, all experimental equilibria (in parallel if allowed)
! lut is set to the output unit
    character cline*(*)
    integer last
    integer lut
!\end{verbatim}
    character line*80,model*72,ch1*1
    integer i1,i2,idef,iz,j3,j4,jp,leak,ll,mode,svss,threads
    double precision gms,xxx,xxy,xxz
    logical listzero
    TYPE(gtp_equilibrium_data), pointer :: neweq
    type(gtp_phasetuple), pointer :: phtup
! rather complex to handle both parallel on non-parallel and with/without 
! griminimizer ...
     if(allocated(firstash%eqlista)) then
        call gparcdx('With global minimizer? ',cline,last,1,ch1,'N',&
             '?Calculate all')
! mode=0 is without grid minimizer ?? mode=-1 ??
        mode=1
        if(ch1.eq.'N' .or. ch1.eq.'n') mode=0
! Seach for memory leaks
        call gparidx('How many times? ',cline,last,leak,1,'Calculate all')
! leak<0 means forever ... not allowed but leak=-1 generates output
        iz=leak
        if(leak.lt.1) leak=1
! Minimize output
        listzero=.false.
! allow output file, if idef>1 no output
        idef=leak
        lut=optionsset%lut
        jp=0
        i2=0
! if compiled with parallel and gridminimizser set then calculate
! sequentially to create composition sets
! TEST THIS IN PARALLEL !!!
        call cpu_time(xxx)
        call system_clock(count=j4)
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
2060    continue
        gridmin: if(mode.eq.1 .or. &
             btest(globaldata%status,GSNOPAR)) then
! if we use grid minimizer do not use parallel even if compiled with OpenMP
           do i1=1,size(firstash%eqlista)
              neweq=>firstash%eqlista(i1)%p1
              jp=jp+1
              if(neweq%weight.eq.zero) then
                 if(listzero) write(kou,2050)neweq%eqno,neweq%eqname
2050             format('Zero weight equilibrium ',i4,2x,a)
              else
                 i2=i2+1
                 call calceq3(mode,.FALSE.,neweq)
                 if(gx%bmperr.ne.0) then
                    write(kou,2051)gx%bmperr,neweq%eqno,neweq%eqname
2051                format(' *** Error code ',i5,' for equilibrium ',&
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
                    call get_state_var_value('GMS ',gms,model,neweq)
                    write(lut,2052)neweq%eqno,&
                         trim(neweq%eqname),neweq%tpval(1),gms,trim(line)
2052                format(i4,2x,a,', T=',F8.2,', GMS= ',1pe12.4,', ',a)
                 endif
              endif
! extra symbol calculations ....
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
!----- $OMP parallel do private(neweq)
!$OMP parallel do private(neweq,gms)
           paraloop: do i1=1,svss
! the error code must be set to zero for each thread ?? !!
              jp=jp+1
              gx%bmperr=0
              neweq=>firstash%eqlista(i1)%p1
! it seems stupid to get this value each loop but outside it is unity
!$                   threads=omp_get_num_threads()
              if(neweq%weight.eq.zero) then
                 if(listzero) write(kou,2050)neweq%eqno,neweq%eqname
              else
! skip this output
!-!$                     if(.TRUE. .and. idef.eq.1) then
!$                     if(.TRUE. .and. iz.lt.0) then
! output only if "number of times" is negative above
!$                      write(*,663)'Equil/loop/thread/maxth/error: ',&
!$                             neweq%eqname,i1,omp_get_thread_num(),&
!$                             threads,gx%bmperr
663              format(a,a,5i5)
! calceq3 gives no output
!$                        call calceq3(mode,.FALSE.,neweq)
!$                     else
! note first argument zero means do not use grid minimizer
                 call calceq3(mode,.FALSE.,neweq)
!$                     endif
                 i2=i2+1
                 line=' '
                 if(gx%bmperr.ne.0) then
                    write(kou,2051)gx%bmperr,neweq%eqno,neweq%eqname
                    gx%bmperr=0
                 elseif(idef.eq.1) then
                    if(threads.eq.1) then
                       jp=1
                       do j3=1,nooftup()
                          phtup=>phasetuple(j3)
    if(neweq%phase_varres(phtup%lokvares)%phstate.ge.PHENTSTAB) then
        call get_phasetup_name(j3,line(jp:))
        jp=len_trim(line)+2
     endif
                       enddo
                    endif
                    call get_state_var_value('GMS ',gms,model,neweq)
                    write(lut,2052)neweq%eqno,&
                         trim(neweq%eqname),neweq%tpval(1),gms,trim(line)
! Listing extra'
                    call list_equilibrium_extra(lut,neweq,plotunit0)
                    if(gx%bmperr.ne.0) then
                       gx%bmperr=0
                    endif
                 endif
              endif
           enddo paraloop
!- $OMP end parallel do not needed???
! OPENMP parallel end loop
! allow composition sets to be created again
!$             globaldata%status=ibclr(globaldata%status,GSNOACS)
!$             globaldata%status=ibclr(globaldata%status,GSNOREMCS)
        endif gridmin
        call cpu_time(xxz)
        call system_clock(count=ll)
        xxy=ll-j4
! or should i2 be used ??
        write(*,669)i2,(xxz-xxx)/i2,xxy/i2
669   format(/'Calculated ',i8,' equlibria, average CPU and clock time',&
           F12.8,1x,F9.5)
! repeat this until leak is zero, if leak negative never stop.
        leak=leak-1
        if(leak.ne.0) then
           goto 2060
        endif
!
        call system_clock(count=ll)
! ?? jp ??             write(kou,664)jp,xxz-xxx,ll-j4,threads
        write(kou,664)xxz-xxx,ll-j4,threads
664     format('Total CPU time: ',1pe12.4,' s and ',i7,' clockcycles',&
             ' using ',i4,' thread(s)')
! this unit may have been used to extract calculated data for plotting
        if(plotunit0.gt.0) then
           write(kou,670)
670        format('Closing a GNUPLOT file oc_many0.plt'/&
                'that may need some editing before plotting')
           write(plotunit0,665)graphopt%plotend
665        format('e'/a)
           close(plotunit0)
! UNFINISHED possibly we could reopen the file again and make oopies 
! of tha data to avoid manual editing
        endif
     else
        write(kou,*)'You must first SET RANGE of experimental equilibria'
     endif
    return
  end subroutine ui_calculate_all_equilibria

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_amend_assessment_result(cline,last)
! AMEND ASSESSMENT_RESULT, analyze the effect of changing a coefficient
    character cline*(*)
    integer last
!\end{verbatim}
    integer i2,j2
    double precision xxx,xxy,xxz
    if(.not.allocated(firstash%eqlista)) then
       write(kou,*)'No assessment record'
       return
    elseif(nvcoeff.le.0) then
       write(kou,*)'No variable optimizing coefficients'; return
    elseif(nvcoeff.ne.nvcoeffdone) then
       write(kou,*)'No optimization made with these coefficients',&
            nvcoeff,nvcoeffdone
       return 
    elseif(mexp.ne.mexpdone) then
       write(kou,*)'No optimization made with these experiments',&
            mexp,mexpdone
       return 
    endif
    call gparix('Index of coefficent to change: ',cline,last,&
         analyze,NONE,'?Amend assess result')
    if(buperr.ne.0) return
    xxy=zero
    if(analyze.lt.0) then
! using a negative coefficient, restore saved coefficients
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
       return
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
                firstash%coeffrsd(j2)=zero
             endif
             if(analyze.eq.j2) then
                cormatix=i2
                xxy=savedcoeff(1,i2)*savedcoeff(3,i2)
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
       if(allocated(savedcoeff)) deallocate(savedcoeff); return
    endif
! ask for new value with the current value as default
    call gparrdx('New value: ',cline,last,xxx,xxy,'?Amend assess result')
    delta=(xxx-xxy)/firstash%coeffscale(analyze)
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
! %coeffvalues should be of the order 1
! No change of %coeffstart and %coeffscale
          firstash%coeffvalues(j2)=xxy
          xxz=xxy*firstash%coeffscale(j2)
! optimizing coefficients are also TP functions, we must update the
! TP function value!! I do not understand this and "list tp" is wrong
! but it seems to work.  If I set the value *firstash%coeffscale it blows up!
          call change_optcoeff(firstash%coeffindex(j2),xxz)
! set RSD to zero
          firstash%coeffrsd=zero
       endif
    enddo
    write(*,*)'To calculate a new set of errors use OPTIMIZE'
    return
  end subroutine ui_amend_assessment_result

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_amend_optimizing_coeffs(cline,last,lut)
! AMEND OPTIMIZING_COEFFS, rescale or recover coefficients, lut is the output unit
    character cline*(*)
    integer last
    integer lut
!\end{verbatim}
    character ch1*1
    integer j2
    double precision xxx
    if(.not.allocated(firstash%eqlista)) then
       write(*,*)'No assessment record allocated'; return
    endif
    call gparcdx('Should the coefficients be rescaled?',&
         cline,last,1,ch1,'N','?Amend optim coeff')
    if(ch1.eq.'y' .or. ch1.eq.'Y') then
! set start values to current values
! Note the "current value" is "start value" times "scaling factor"
       do j2=0,size(firstash%coeffstate)-1
          if(firstash%coeffstate(j2).ge.10) then
             call get_value_of_constant_index(firstash%coeffindex(j2),xxx)
             if(gx%bmperr.ne.0) then
             write(*,*)'Error getting value of assessment coefficient',j2
                return
             endif
! Set all values equal to the current value of the TP variable ...
             firstash%coeffscale(j2)=xxx
             firstash%coeffstart(j2)=xxx
             firstash%coeffvalues(j2)=one
          endif
       enddo
       firstash%coeffrsd=zero
       call listoptcoeff(mexp,err0,.FALSE.,lut)
       if(allocated(cormat)) then
          deallocate(cormat)
          deallocate(tccovar)
       endif
    else
       call gparcdx('Do you want to recover the coefficients values?',&
            cline,last,1,ch1,'N','?Amend optim coeffs')
       if(ch1.eq.'y' .or. ch1.eq.'Y') then
! set current optimizing values back to start values
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
557       format('Nothing done as there are no other amend',&
               ' optimizing option')
       endif
    endif
    return
  end subroutine ui_amend_optimizing_coeffs

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_optimize(cline,last,lut)
! OPTIMIZE (and continue) using LMDIF, lut is the output unit
    character cline*(*)
    integer last
    integer lut
!\end{verbatim}
    integer i1,i2,j2,j4,ll
    double precision xxx
    call gparidx('Number of iterations: ',cline,last,i1,nopt1,&
         '?Optimize')
    if(buperr.ne.0) return
    nopt1=i1
    nopt=i1
606 format(a,10i4)
!
    if(.not.allocated(firstash%eqlista)) then
       write(kou,*)'There are no equilibria with experiments!'
       return
    endif
    if(allocated(coefs)) deallocate(coefs)
    if(allocated(errs))  deallocate(errs)
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
                return
             endif
             coefs(i2)=firstash%coeffvalues(i1)
! We do not have to bother about the associtated TP variable, it will
! be set by the calfun routine to coefs*firstashscale
             if(gx%bmperr.ne.0) then
                write(*,*)'Error finding coefficient TP fun'
                return
             endif
          endif
       enddo
       if(i2.lt.nvcoeff) then
          write(kou,*)'Internal error for variable coefficients',&
               i2,nvcoeff
          return
       endif
    endif
! JUMP HERE IF CONTINUE optimization  ... NOT YET implemented
987 continue
! mexp    Number of experiments
! nvcoeff Number of coefficients to be optimized
! errs Array with differences with experiments and calculated values
! coefs Array with coefficients
    if(mexp.le.0 .or. nvcoeff.le.0) then
       write(kou,569)mexp,nvcoeff
569    format('Cannot optimize with zero experiments or coefficients',2i5)
       return
    endif
    firstash%lwam=lwam
    write(*,558)mexp,nvcoeff,mexp*nvcoeff+5*nvcoeff+mexp,lwam
558 format(/'*************************************************************'/&
         '>>>   Start of optimization using LMDIF'/&
         '>>>   with ',i4,' experiments and ',i3,' coefficients.',/&
         '>>>   Workspace needed ',i6,', out of allocated ',i6/&
         '*************************************************************')
!
    j4=nopt
    if(.not.allocated(iwam)) then
! value of lwam set by user
       allocate(iwam(lwam))
       allocate(wam(lwam))
    endif
    if(allocated(fjac)) deallocate(fjac)
! fjac is used to calculate the Jacobian and other things
! err0(1) is set to the sum of errors squared for the initial values of coefs
573 format(a,6(1pe12.4))
    allocate(fjac(mexp,nvcoeff))
!->->->->->-> HERE THE OPTIMIZATION IS MADE <-<-<-<-<-<-
! nfev set to number of iterations
    call lmdif1(calfun,mexp,nvcoeff,coefs,errs,optacc,nopt,nfev,&
         iwam,wam,lwam,fjac,err0)
!->->->->->-> HERE THE OPTIMIZATION IS MADE <-<-<-<-<-<-
    mexpdone=mexp
    nvcoeffdone=nvcoeff
! on return nopt is set to a message but 
! first copy the coefs to coeffvalues ...
    i2=0
    do i1=0,size(firstash%coeffstate)-1
       if(firstash%coeffstate(i1).ge.10) then
          i2=i2+1
          firstash%coeffvalues(i1)=coefs(i2)
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
! cormat will be the CORRELATION MATRIX if optimization successful
! otherwise it will not be allocated
    if(allocated(cormat)) then
       deallocate(cormat)
       deallocate(tccovar)
    endif
!--------------- begin calculate correlation matrix and RSD
! zero the relative standard deviations (RSD)
    firstash%coeffrsd=zero
    if(j4.gt.0 .and. nopt.gt.0 .and. nopt.le.6) then
! if there is a result calculate the Jacobian in fjac
! mexp,nvcoeff,coeffs,errs are same as in the call to lmdif1
! This will overwrite the fjac returned from the call to lmdif1
! allocate array to extract calculated values of experiments
       if(allocated(calcexp)) deallocate(calcexp)
       allocate(calcexp(mexp))
       iflag=2
! penulitima argument zero means use machine precision to calculate derivative
       call fdjac2(calfun,mexp,nvcoeff,coefs,errs,fjac,mexp,iflag,zero,wam)
! debug output ...
563    format(6(1pe12.4))
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
564             format(a,3i5,1pe12.4)
             enddo
! this matrix is symmetric ... which index first ???
             cov1(j2,i2)=xxx
          enddo
       enddo
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
       elseif(abs(cov1(1,1)).gt.1.0D-38) then
! cov1 is just a single value
          allocate(cormat(1,1))
          allocate(tccovar(1,1))
! IF THERE IS A SINGLE VARIABLE ITS CORRELATION MATRIX MUST BE UNITY
          cormat(1,1)=one
          tccovar(1,1)=one
       else
          write(*,*)'Correlation matrix singular'
       endif
    endif
! write the correlation matrix  this is still very uncertain ,,,
! zero all RSD values
    firstash%coeffrsd=zero
    if(allocated(cormat) .and. allocated(cov1)) then
! calculate the RSD (Relative Standard Deviation) for each parameter
! the last calculated values of the experiments in calcexp
       xxx=zero
       do i2=1,mexp
! the calculated value is stored in calcexp by fdjac if calcexp is allocated!
          xxx=xxx+calcexp(i2)
766       format('pmon: Calculated value',i4,2(1pe12.4))
       enddo
! This value may be negative!
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
! we have already multiplied all terms in covariance matrix with err0(3)
             firstash%coeffrsd(i1)=abs(sqrt(abs(cov1(i2,i2))))
          endif
       enddo
    endif
! deallocate calcexp to avoid storing these values when running LMDIF
    if(allocated(calcexp)) deallocate(calcexp)
!--------------- end calculate correlation matrix and RSD
! some nice output .....
    write(kou,5020)
    if(j4.eq.0) then
       write(*,*)'Dry run with zero iterations'
    elseif(nopt.eq.0) then
       write(kou,5000)nopt
5000   format(/'*** No optimization due to improper input parameters',i3)
    elseif(nopt.eq.1) then
       write(kou,5001)nopt,optacc
5001   format(/'LMDIF return code ',i2/&
            'Relative error for sum of squares is within ',1pe10.2)
    elseif(nopt.eq.2) then
       write(kou,5002)nopt,optacc
5002   format(/'LMDIF return code ',i2/&
            'Relative error of parameters is within ',1pe10.2)
    elseif(nopt.eq.3) then
       write(kou,5003)nopt
5003   format(/'LMDIF return code ',i2,': successful optimization')
    elseif(nopt.eq.4) then
       write(kou,5004)nopt
5004   format(/'*** LMDIF return code ',i2/&
            'Sum of squares does not decrease')
    elseif(nopt.eq.5) then
       write(kou,5005)nopt,nfev
5005   format(/'*** LMDIF return code ',i2/&
            'Maximum calls of function ',i5,' exceeded')
    elseif(nopt.eq.6) then
       write(kou,5006)nopt,optacc
5006   format(/'*** LMDIF return code ',i2/&
            'Cannot reduce error, requested accuracy ',1pe10.2,' too small')
! '*** Cannot reduce error, requested accuracy 123456789. too small
    elseif(nopt.eq.6) then
       write(kou,5007)nopt,optacc
5007   format(/'*** LMDIF return code ',i2/&
           'Cannot improve result, requested accuracy ',1pe10.2,' too small')
    else
       write(kou,5008)nopt
5008   format('*** LMDIF return code ',i7/&
            'Unknown code, see LMDIF documentation.')
    endif
    write(kou,5010)nfev,err0
5010 format(/'Iterations ',i4,', sum of errors changed from ',&
         1pe14.6,' to ',1pe14.6/17x,'Normalized sum of errors:',20x,1pe14.6)
    write(kou,5020)
5020 format(/78('*'))
! finally list the coefficient values
    call listoptcoeff(mexp,err0,.FALSE.,lut)
    return
  end subroutine ui_optimize

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_new(cline,last,defcp,reinit)
! NEW, remove all data.  reinit is TRUE if data removed and OC must be reinitiated
! defcp is the default for CALCULATE PHASE
    character cline*(*)
    integer last
    integer defcp
    logical reinit
!\end{verbatim}
    character ch1*1
    integer jp
    integer intv(10)
    double precision dblv(10)
    reinit=.false.
!=================================================================
! one must deallocate everyting explicitly to use memory again
    call gparcdx('All data will be removed, are you sure?',cline,last,&
         1,ch1,'N','?New')
    if(ch1.ne.'Y') then
       write(kou,*)'*** NO CHANGE, upper case Y needed for NEW'
       return
    endif
! remove global check during map
    mapglobalcheck=0
    stepspecial=.FALSE.
!------remove assessment data
    if(allocated(firstash%eqlista)) then
       write(*,*)'Assessment data removed, not deallocated: memory leak'
    endif
    if(allocated(firstash%eqlista)) deallocate(firstash%eqlista)
    deallocate(firstash)
    if(allocated(errs)) deallocate(errs)
    mexp=0
    updatemexp=.true.
    nvcoeff=0
! initiate the limit on number of equilibria saved during step/map
    totalsavedceq=0
!----- deleting map results ...
    if(associated(maptopsave)) then
! this is necessary only if no plot of last step/map made ...
       write(kou,*)'We link to maptopsave'
       maptop%plotlink=>maptopsave
       nullify(maptopsave)
    endif
    seqxyz=0
    call delete_mapresults(maptop)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error deleting map results! Report this error with macro!'
       stop
    endif
    nullify(maptop)
    nullify(maptopsave)
    seqxyz=0
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
!
! special MQMQA multivalent check
    mqmqa_multival=.false.
! this routine fragile, inside new_gtp init_gtp is called
    call new_gtp
    if(gx%bmperr.ne.0) then
       write(*,*)'Error deleting data! Report this error with macro!'
       stop
    endif
    write(kou,*)'All data removed, reinitiating'
    call init_gtp(intv,dblv)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error initiating! Report this error with macro!'
       stop
    endif
    reinit=.true.
    return
  end subroutine ui_new

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_step_tzero(cline,last,tzcond,name2,ceq)
! STEP TZERO, calculate a T0 line between two phases
! tzcond is the index of the condition to release
! name2 is returned as the second phase, used as default by CALCULATE TZERO
    character cline*(*)
    integer last
    integer tzcond
    character name2*(*)
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24
    integer ics,iph,iph2,j2,jp
    write(kou,871)
871 format('For this command you must already have used',&
         ' "calculate tzero"'/&
         'for the two phases you will specify below and you must',&
         ' have specified an axis'/&
         'with the composition of the fast diffusing element.')
    call gparcx('Have you done all that?',cline,last,1,&
         name1,'NO','?Step Tzero')
    call capson(name1)
    if(name1(1:1).ne.'Y') return
    call gparcx('First phase name: ',cline,last,1,name1,' ',&
         '?Step Tzero')
    if(buperr.ne.0) return
    call find_phase_by_name(name1,iph,ics)
    if(gx%bmperr.ne.0) return
    if(ics.ne.1) then
       write(*,*)'You must use first composition set'
       return
    endif
    call gparcx('Second phase name: ',cline,last,1,name2,' ',&
         '?Step Tzero')
    call find_phase_by_name(name2,iph2,ics)
    if(gx%bmperr.ne.0) return
    if(ics.ne.1) then
       write(*,*)'You must use first composition set'
       return
    endif
! normally T is the first condition
    j2=1
    call gparidx('Release condition number',cline,last,tzcond,j2,&
         '?Step Tzero')
! Delete previous step/map results
    call ui_delete_step_map(0,cline,last,ceq)
! step tzero
    stepspecial(3)=.TRUE.
    call step_tzero(maptop,noofaxis,axarr,seqxyz,iph,iph2,tzcond,ceq)
    if(gx%bmperr.ne.0) return
! sum the points calculated
    jp=maptop%linehead(1)%number_of_equilibria+&
         maptop%linehead(2)%number_of_equilibria
    write(kou,'(a,i5,a)')'Calculated ',jp,' points along the tzero line'
    return
  end subroutine ui_step_tzero

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_step_scheil(cline,last,kom2,ceq)
! STEP SCHEIL_GULLIVER (kom2=7) or STEP FAST (kom2=9, with fast diffusing elements)
    character cline*(*)
    integer last
    integer kom2
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,line*80,elsym*2
    character*2 fast(3)
    integer iel,nv
    write(kou,872)
872 format('Before this command you must have set an alloy composition',&
         ' and calculated'/&
         'an equilibrium in the liquid and have set an axis with T',&
         ' as variable.')
    call gparcx('Have you done all that?',cline,last,1,&
         name1,'NO','?Step Scheil')
    call capson(name1)
    if(name1(1:1).ne.'Y') return
! Delete previous step/map results
    call ui_delete_step_map(0,cline,last,ceq)
    write(kou,873)
873 format('The simulation will decrease T and change the liquid',&
         ' composition depending'/&
         'on the solids formed until there is no liquid stable.')
! step scheil
    stepspecial(2)=.TRUE.
    if(kom2.eq.9) then
! ask for fast diffusing elements in Scheil simulation
       nv=1
       fast=' '
       line='A fast diffusing element: '
       fastel: do while(.TRUE.)
          call gparcx(trim(line),cline,last,1,elsym,' ','?Step Scheil')
          if(elsym(1:1).ne.' ') then
             call capson(elsym)
             call find_element_by_name(elsym,iel)
             if(gx%bmperr.ne.0) return
             fast(nv)=elsym
             nv=nv+1
             if(nv.gt.3) exit fastel
             line='Another fast diffusing element: '
          else
             exit fastel
          endif
       enddo fastel
       call step_scheil2(maptop,noofaxis,axarr,seqxyz,fast,ceq)
    else
! step scheil with no fast diffusiing elements
       call step_scheil(maptop,noofaxis,axarr,seqxyz,ceq)
    endif
    if(gx%bmperr.ne.0) return
    return
  end subroutine ui_step_scheil

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_map(cline,last,ceq)
! MAP, calculate a phase diagram with 2 or 3 axis
! ceq may be changed to the start equilibrium
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character ch1*1
    integer ll
! maybe disable continue optimization ??
    if(noofaxis.lt.2) then
       write(kou,*)'You must set two axis with independent variables'
       return
    endif
    if(noofaxis.gt.3) then
       write(kou,*)'More than 3 axis not implemented yet'
       return
    endif
    stepspecial=.FALSE.
! indicate to graphics we have calculated with 2 axis
    graphopt%noofcalcax=noofaxis
    write(kou,20014)
20014 format('The map command is fragile, please send problematic diagrams',&
         ' to the',/'OC development team'/)
! when setting logfile the maptop became associated !! 
    call ui_delete_step_map(2,cline,last,ceq)
! maptop is returned as main map/step record for results
! noofaxis is current number of axis, axarr is array with axis data
! starteq is start equilibria, if empty set it to ceq
    if(noofstarteq.eq.0) then
       noofstarteq=1
       starteqs(1)%p1=>ceq
    endif
    ceq=>starteqs(1)%p1
    ll=degrees_of_freedom(ceq)
    if(ll.ne.0) then
       write(*,*)'Degrees of freedom not zero ',ll
       return
    endif
! maptop is first nullified inside map_setup, then alloctated to return result
    call map_setup(maptop,noofaxis,axarr,seqxyz,starteqs)
    if(gx%bmperr.ne.0) then
       write(kou,*)'Error return from MAP: ',gx%bmperr
       gx%bmperr=0
    endif
    if(.not.associated(maptop)) then
! if one has errors in map_setup maptop may not be initiated, if one
! has saved previous calculations in maptopsave restore those
       if(associated(maptopsave)) then
          write(kou,*)'Restoring previous map results'
          maptop=>maptopsave
          nullify(maptopsave)
       endif
    elseif(associated(maptopsave)) then
       write(kou,'(a)')'Link set to previous map results'
       maptop%plotlink=>maptopsave
       nullify(maptopsave)
    endif
! remove start equilibria
    noofstarteq=0
    nullify(starteqs(1)%p1)
! mark that interactive listing of conditions and results may be inconsistent
    ceq%status=ibset(ceq%status,EQINCON)
    if(gx%bmperr.ne.0) return
    return
  end subroutine ui_map

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_plot_format_file(cline,last,kom2)
! PLOT GRAPHICS_FORMAT (kom2=7) and/or OUTPUT_FILE (kom2=8)
    character cline*(*)
    integer last
    integer kom2
!\end{verbatim}
    character ch1*1,aline*128,filename*128
    integer grunit,i1,jp,ztyp
    logical logok,once
    if(kom2.eq.7) then
! subroutine TOPHLP forces return with ? in position cline(1:1)
29130  continue
       call gparidx('Graphics format index:',cline,last,grunit,1,&
            '?Plot formats')
       if(cline(1:1).eq.'?' .or. &
            grunit.lt.1 .or. grunit.gt.graphopt%gnutermax) then
          write(kou,29133)
29133     format('Avalable graphics formats are:')
          write(kou,29135)(i1,graphopt%gnutermid(i1),&
               i1=1,graphopt%gnutermax)
29135     format(i3,2x,a)
          goto 29130
       endif
       graphopt%gnutermsel=grunit
       write(kou,*)'Graphics format set to: ',graphopt%gnutermid(grunit)
    endif
!-----------------------------------------------------------
! PLOT OUTPUT_FILE, always asked when changing graphics terminal type
21140 continue
! default extension: 1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT, 8=LOG
! negative is for write, 0 read without filter, -100 write without filter
! DO NOT USE tinyfiledialog here ...
    write(*,*)'To use file the browser give just a <'
    call gparcdx('Plot file',cline,last,1,plotfile,'ocgnu','?Plot file')
! to avoid confusion abot use of > and <
    if(plotfile(1:1).eq.'<' .or. plotfile(1:1).eq.'>') then
! use the file browser
       ztyp=-5
       call gparfilex('File name: ',cline,last,1,plotfile,' ',ztyp,&
            '?Plot file')
! make sure there is a plt extention
       jp=index(plotfile,'.')
       if(jp.le.0) then
          jp=len_trim(plotfile)
          plotfile(jp+1:)='.plt'
       endif
       write(*,*)'Output will be on: ',trim(plotfile)
    endif
    once=.false.
    if(plotfile(1:2).eq.'./') then
! save in macro directory if iumaclevl>0, else in current working directory
       if(iumaclevl.gt.0) then
! we are executing a macro, skip the ./
          aline=plotfile(3:)
          plotfile=trim(macropath(iumaclevl))//aline
       else
! running interactivly prefix with working directory (default?)
          aline=plotfile(2:)
          plotfile=trim(workingdir)//aline
       endif
! trouble passing on ling file names ....
       write(*,*)'PMON working directory: ',trim(workingdir)
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
       inquire(file=filename,exist=logok)
       if(logok) then
          call gparcdx('File exists, overwrite?',&
               cline,last,1,ch1,'N','PLOT file')
          if(.not.(ch1.eq.'Y' .or. ch1.eq.'y')) then
             write(*,133)
133          format('Please use another file name')
             plotfile=' '
             goto 21140
          endif
          write(*,134)trim(filename)
134       format(/'Overwriting previous results on ',a)
          once=.true.
       endif
    endif
! I am not sure how to inform user where the plot file is saved ....
    return
    return
  end subroutine ui_plot_format_file

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_plot_append(cline,last,string)
! PLOT APPEND a gnuplot or csv file (a csv file is converted to csvappend.plt)
! string is only passed to keep the old (strange) 3rd header line of csvappend.plt
    character cline*(*)
    integer last
    character string*(*)
!\end{verbatim}
    character ch1*1,name1*24,filename*128,optres*40
    integer ioc,ip,jp,ztyp
    write(kou,*)'Give a file name with graphics in .plt or .csv format'
! append plot file, specifying extension PLT
! default extension (1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT, 6=XTDB, 7=DAT
! negative is for write, 0 read without filter, -100 write without filter
    ztyp=5
    call gparfilex('File name',cline,last,1,filename,'  ',ztyp,&
         '?Plot append')
! check file exits, convert .csv to .plt, and add .plt if necessary ...
    jp=max(index(filename,'.csv '),index(filename,'.CSV '))
    write(*,*)'Full file name: ',trim(filename)
    if(jp.gt.0) then
! A csv file must be converted to a plt file ------------- csv special begin
! default separator in FactSage is ";"
       ch1=';'
       call gparcdx('Separating character (, ; or ?)',cline,last,&
            1,ch1,';','?CSV separator')
       write(*,21205)trim(filename)
21205  format(/'Converting csv file: ',a,' to csvappend.plt')
       open(23,file=filename,status='old',access='sequential',err=21300)
! First column is x-axis, the remaining columns are y-axis
! create a new file for the GNUPLOT, overwrite any old
       open(31,file='csvappend.plt',status='unknown',access='sequential',&
            err=21305)
! write header
21210  format(a)
       call date_and_time(optres,name1)
       write(31,21220)trim(filename),optres(1:4),optres(5:6),&
            optres(7:8),name1(1:2),name1(3:4),trim(string)
21220  format('# Converted by OC from csv file: ',a/&
            '# ',a4,'-',a2,'-',a2,2x,a2,'h',a2/&
            '# ',a//&
            'set terminal wxt size 840,700 font "Arial,16"'/&
            'set title "OpenCalphad CSV" '/&
            'set origin 0.0, 0.0'/&
            'set size   1.0, 1.0'//&
            '$OCCSV2502000 << EOD')
! on the second line to the last one are values to be plotted as symbols
       read(23,21210,end=21250)string
       ioc=1; ip=1; jp=1
       commas: do while(.true.)
! ch1 should be the separating character
          ip=index(string(ip:),ch1)
          if(ip.le.0) exit commas
! replace ch1 by a space
          string(jp+ip-1:jp+ip-1)=' '
          ioc=ioc+1
          ip=jp+ip+1
          jp=ip
       enddo commas
! there are ioc values to be plotted
       do while(.true.)
! add a space at the end of line ... GNUPLOT need that
          write(31,'(a," ")')trim(string)
          ip=1; jp=1
! we have to replace ch1 with a space on every line ....
          read(23,21210,end=21250)string
          tty17: do while(ip.gt.0)
             jp=ip
             ip=index(string(ip:),ch1)
             if(ip.eq.0) exit tty17
             ip=jp+ip-1
             string(ip:ip)=' '
          enddo tty17
       enddo
! end of file
! we have read the whole csv file, close the PLT file
21250  continue
       close(23)
! add the GNUPLOT ending
       write(31,21260)ioc
21260  format('EOD'//&
            'plot for [i=2:',i3,&
            '] $OCCSV2502000 using 1:i with points pt (i-1) notitle'/)
       close(31)
! set the newly created PLT file as appendfile
       filename='csvappend.plt'
    endif
!------------------------ end csv special
    jp=index(filename,'.plt ')
    if(jp.le.0) then
       jp=len_trim(filename)
       filename(jp+1:)='.plt'
    endif
! test file exists by opening and closing it
    open(23,file=filename,status='old',access='sequential',err=21330)
    close(23)
    graphopt%appendfile=filename
    return
! error opening file, remove any previous appended file
21300 continue
    if(graphopt%appendfile(1:1).ne.' ') then
       write(*,21304)trim(graphopt%appendfile)
21304  format('Error, removing append file: ',a)
    else
       write(kou,21303)trim(filename)
21303  format('Error opening file: ',a/'Maybe path or "./" missing?: ')
    endif
    graphopt%appendfile=' '
    return
21305 continue
    write(*,21306)
21306 format('Error opening csvappend.plt file')
    return
21330 continue
    write(*,21331)trim(filename)
21331 format('Error opening file, maybe missing path or "./" ? ',a)
    return
    return
  end subroutine ui_plot_append

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_plot_text(cline,last,ceq)
! PLOT TEXT_LABEL, add or modify a text on the plot
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character ch1*1,line*80,text*72
    integer j4,jp,kl
    double precision textfontscale,xxx,xxy
    type(graphics_textlabel), pointer :: labelp
    labelp=>graphopt%firsttextlabel
    if(associated(labelp)) then
       call gparcdx('Modify existing text?',cline,last,1,ch1,'NO',&
            '?Plot texts')
       if(ch1.eq.'y' .or. ch1.eq.'Y') then
          jp=0
          do while(associated(labelp))
             jp=jp+1
             write(kou,2310)jp,labelp%xpos,labelp%ypos,&
                  labelp%textfontscale,labelp%angle,&
                  trim(labelp%textline)
2310         format(i3,2(1pe12.4),2x,0pF5.2,2x,i4,5x,a)
             labelp=>labelp%nexttextlabel
          enddo
          call gparidx('Which text index?',cline,last,kl,1,'?Plot texts')
          if(kl.lt.1 .or. kl.gt.jp) then
             write(*,*)'No such text label'
             return
          endif
          labelp=>graphopt%firsttextlabel
          do jp=2,kl
             labelp=>labelp%nexttextlabel
          enddo
          call gparcdx('New text: ',cline,last,5,text,&
               labelp%textline,'?Plot texts')
          labelp%textline=trim(text)
          call gparrdx('New X position: ',cline,last,xxx,&
               labelp%xpos,'?Plot texts')
          call gparrdx('New Y position: ',cline,last,xxy,&
               labelp%ypos,'?Plot texts')
          call gparrdx('New Fontscale: ',cline,last,&
               textfontscale,labelp%textfontscale,'?Plot texts')
          if(textfontscale.lt.0.2) textfontscale=0.2
          call gparidx('New angle (degrees): ',cline,last,j4,&
               labelp%angle,'?Plot texts')
          if(buperr.ne.0) then
             write(*,*)'Error reading coordinates'; buperr=0; return
          endif
          labelp%xpos=xxx
          labelp%ypos=xxy
          labelp%textfontscale=textfontscale
          labelp%angle=j4
! ask for more options
          return
       endif
    endif
! input a new label
    call gparrdx('X position: ',cline,last,xxx,zero,'?Plot texts')
    call gparrdx('Y position: ',cline,last,xxy,zero,'?Plot texts')
    call gparrdx('Fontscale: ',cline,last,textfontscale,0.8D0,&
         '?Plot texts')
    if(textfontscale.le.0.2) textfontscale=0.2
    call gparidx('Angle (degree): ',cline,last,j4,0,'?Plot texts')
    if(buperr.ne.0) then
       write(*,*)'Error reading coordinates'; buperr=0; return
    endif
    line=' '
    if(noofaxis.eq.2) then
! Calculate the equilibria at the specific point
       write(kou,22100)
22100  format(' *** Note: the positioning of the text will use the ',&
            'axis variables for which',/11x,'the diagram was calculated',&
            ' even if you plot with other variables!')
       call gparcdx('Do you want to calculate the equilibrium? ',&
            cline,last,1,ch1,'Y','?Plot texts')
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
! There is no gparcd which allows editing the existing text ... use emacs!!
    text=' '
    call gparcdx('Text: ',cline,last,5,text,line,'?Plot texts')
    if(text(1:1).eq.' ') then
       write(*,*)'Label ignored'
       return
    endif
! I know one should never allocate pointers but this is the only way ???
    allocate(textlabel)
    textlabel%xpos=xxx
    textlabel%ypos=xxy
    textlabel%textfontscale=textfontscale
    textlabel%angle=j4
    textlabel%textline=trim(text)
    if(associated(graphopt%firsttextlabel)) then
       textlabel%nexttextlabel=>graphopt%firsttextlabel
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
    return
    return
  end subroutine ui_plot_text

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_list_excel_csv(cline,last,version,ceq)
! LIST EXCEL_CSV, write STEP results as a CSV table
    character cline*(*)
    integer last
    character version*(*)
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character text*72,model*72,aline*128,longstring*5000
    integer i1,iax,istv,jp,ztyp
    double precision xxx
    double precision, dimension(maxconst) :: yarr
    logical wildcard
    if(noofaxis.gt.1 .or. .not.associated(maptop)) then
       write(kou,*)'You must give a STEP command before list excell_csv'
       return
    endif
    wildcard=.FALSE.
    iax=1
    jp=1
    call get_one_condition(jp,text,axarr(iax)%seqz,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'PMON Error getting axis condition from index: ',&
            iax,axarr(iax)%seqz
       return
    endif
! we just want the expression, remove the value including the = sign
    jp=index(text,'=')
    text(jp:)=' '
    axplotdef(1)=text
    if(maptop%tieline_inplane.eq.1) then
! if tie-lines in the plane is 1 (.e. YES) and calculating axis was x(A)
! then plot axis should be x(*,cu) 
       jp=index(text,'(')
       if(jp.gt.0) then 
          text=text(1:jp)//'*,'//text(jp+1:)
       endif
    endif
! default for second axis always NP(*)
    axplotdef(2)='NP(*)'
! the 4th argument to gparc means the following:
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!      7 TEXT TERMINATED BY SPACE OR "," BUT IGNORING SUCH INSIDE ( )
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
    iax=1
    call gparcdx('Independent variable',&
         cline,last,7,axplot(iax),axplotdef(iax),'?List excell CSV')
! dependent variables, can be wildcard
    iax=2
    call gparcdx('Dependent values',&
         cline,last,7,axplot(iax),axplotdef(iax),'?List excell CSV')
    if(buperr.ne.0) return
    if(index(axplot(iax),'*').gt.0 .or. index(axplot(iax),'#').gt.0) then
       wildcard=.TRUE.
    endif
    if(index(axplot(iax),'*').gt.0) then
! generate many values
! the values are returned in yarr with dimension maxconst. 
! longstring are the state variable symbols for the values ...
       call get_many_svar(axplot(iax),yarr,maxconst,i1,longstring,ceq)
       if(gx%bmperr.ne.0) then
! if error go back to command level
          write(kou,*)'Illegal axis variable!  Error code: ',gx%bmperr
          return
       endif
    elseif(index(axplot(iax),'#').gt.0) then
! generate many values including for metastable phases
! the values are returned in yarr with dimension maxconst. 
! longstring are the state variable symbols for the values ...
       call get_many_svar(axplot(iax),yarr,maxconst,i1,longstring,ceq)
       if(gx%bmperr.ne.0) then
! if error go back to command level
          write(kou,*)'Illegal axis variable!  Error code: ',gx%bmperr
          return
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
          call find_svfun(axplot(iax),istv)
          if(gx%bmperr.ne.0) then
             write(kou,*)'Illegal axis variable, error: ',gx%bmperr
             return
          endif
       endif
    endif
! output file, NOTE if /APPEND the file already open!
    if(optionsset%lut.ne.kou) then
       write(*,*)'Appended on CSV files not implemented yet',&
            optionsset%lut
       close(optionsset%lut); optionsset%lut=0
    endif
    if(buperr.ne.0) buperr=0
! What does -5 as argument mean?? Well, open for write!!
    ztyp=-5
    call gparfilex('Output file: ',cline,last,1,plotfile,' ',ztyp,&
         '?List excell CSV')
! make sure there is a file name
    if(len_trim(plotfile).le.0) then
       plotfile=' '
       if(buperr.ne.0) then
          buperr=0
       endif
       write(*,*)'Output on screen'
    else
       jp=index(plotfile,'.')
       if(jp.le.0) then
          jp=len_trim(plotfile)
          plotfile(jp+1:)='.csv'
       endif
       write(*,*)'Output will be on: ',trim(plotfile)
    endif
    if(plotfile(1:2).eq.'./') then
! save in macro directory if iumaclevl>0, else in current working directory
       if(iumaclevl.gt.0) then
! we are executing a macro, skip the ./
          aline=plotfile(3:)
          plotfile=trim(macropath(iumaclevl))//aline
       else
! running interactivly prefix with working directory (default?)
          aline=plotfile(2:)
          plotfile=trim(workingdir)//aline
       endif
! trouble passing on ling file names ....
    endif
1234 continue
! use the graphics record to transfer data ...
    graphopt%pltax(1)=axplot(1)
    graphopt%pltax(2)=axplot(2)
    graphopt%filename=plotfile
! this command only for tabulating STEP commands
    graphopt%status=ibset(graphopt%status,GRCSVTABLE)
    graphopt%status=ibclr(graphopt%status,GRISOPLETH)
! added ceq in the call to make it possible to handle change of reference states
    call ocplot2(jp,maptop,axarr,graphopt,version,ceq)
    graphopt%status=ibclr(graphopt%status,GRCSVTABLE)
    if(gx%bmperr.ne.0) return
    write(*,*)'CSV output saved on file: ',trim(plotfile)
    return
  end subroutine ui_list_excel_csv

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_plot_scale_ranges(cline,last,jump)
! PLOT SCALE_RANGES for X or Y axis
! jump is set to 100 if the monitor should go back to command level
    character cline*(*)
    integer last
    integer jump
!\end{verbatim}
    character ch1*1
    double precision xxx
    logical once,twice
    jump=0
    call gparcdx('For X or Y axis? ',cline,last,1,ch1,'Y','?Plot limits')
    if(ch1.eq.'X' .or. ch1.eq.'x') then
       goto 21120
    elseif(ch1.eq.'Y' .or. ch1.eq.'y') then
       goto 21130
    else
       write(kou,*)'Please answer X or Y'
    endif
    return
!............................................ user limits X axis (1)
21120 continue
    call gparcdx('Default limits',cline,last,1,ch1,'N','?Plot limits')
    if(ch1.eq.'Y' .or. ch1.eq.'y') then
       graphopt%rangedefaults(1)=0
    else
       graphopt%rangedefaults(1)=1
       twice=.FALSE.
21104  continue
       call gparrdx('Low limit',cline,last,xxx,graphopt%dfltmin(1),&
            '?Plot limits')
       if(graphopt%gibbstriangle .and. xxx.ne.zero) then
          write(*,*)'Lower limit of a Gibbs triangle plot must be zero'
          return
       endif
       graphopt%plotmin(1)=xxx
       graphopt%dfltmin(1)=xxx
       once=.TRUE.
21105  continue
       call gparrdx('High limit',cline,last,xxx,&
            graphopt%dfltmax(1),'?Plot limits')
       if(xxx.le.graphopt%plotmin(1)) then
          if(once) then
             write(kou,*)'Think before typing'
             once=.FALSE.
          elseif(twice) then
             write(kou,*)'Back to command level'
             jump=100; return
          else
             write(kou,*)'Please give the low limit again!'
             twice=.TRUE.
             goto 21104
          endif
          write(kou,21106)graphopt%plotmin(1)
21106     format('High limit must be higher than low: ',1pe14.6)
          goto 21105
       endif
       graphopt%plotmax(1)=xxx
       graphopt%dfltmax(1)=xxx
    endif
    return
!---------------------------------------------- user limits Y axis (2)
21130 continue
    call gparcdx('Default limits',cline,last,1,ch1,'N','?Plot limits')
    if(ch1.eq.'Y' .or. ch1.eq.'y') then
       graphopt%rangedefaults(2)=0
    else
       graphopt%rangedefaults(2)=1
       twice=.FALSE.
21107  continue
       call gparrdx('Low limit',cline,last,xxx,graphopt%dfltmin(2),&
            '?Plot limits')
       if(graphopt%gibbstriangle .and. xxx.ne.zero) then
          write(*,*)'Lower limit of a Gibbs triangle plot must be zero'
          return
       endif
       graphopt%plotmin(2)=xxx
       graphopt%dfltmin(2)=xxx
       once=.TRUE.
21108  continue
       call gparrdx('High limit',cline,last,xxx,&
            graphopt%dfltmax(2),'?Plot limits')
       if(xxx.le.graphopt%plotmin(2)) then
          if(once) then
             write(*,*)'Think before typing'
             once=.FALSE.
          elseif(twice) then
             write(kou,*)'Back to command level'
             jump=100; return
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
    return
    return
  end subroutine ui_plot_scale_ranges
!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_plot_axes(cline,last,jump,ceq)
! PLOT, ask for the horizontal and vertical axis variables
! jump is set to 100 or 990 if the monitor should go to that label
    character cline*(*)
    integer last
    integer jump
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,text*72,model*72,longstring*5000
    integer i1,iax,istv,jp
    double precision xxx
    double precision, dimension(maxconst) :: yarr
    logical wildcard
    jump=0
    wildcard=.FALSE.
! values of stepspecial ...
    pltaxdef: do iax=1,2
       plotdefault: if(axplotdef(iax)(1:1).eq.' ') then
! If there is no previous plot axis variable, propose one
          iaxval: if(iax.eq.1 .and. stepspecial(2)) then
! Scheil, PFL (Phase Fraction Liquid) is a special function
             if(iax.eq.1) text='PFL'
          elseif(iax.le.noofaxis) then
! extract the actual axis condition used for calculation
             jp=1
             call get_one_condition(jp,text,axarr(iax)%seqz,ceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Error getting axis condition from index: ',&
                     iax,axarr(iax)%seqz
                jump=990; return
             endif
             jp=index(text,'=')
             text(jp:)=' '
             if(.not.(text(1:2).eq.'MU' .or. text(1:2).eq.'AC' .or.&
                  text(1:4).eq.'LNAC')) then
                if(maptop%tieline_inplane.eq.1) then
! if tie-lines in the plane is 1 (.e. YES) and calculating axis was x(cu)
! then plot axis should be x(*,cu) 
                   jp=index(text,'(')
                   if(jp.gt.0) then 
                      text=text(1:jp)//'*,'//text(jp+1:)
                   endif
                endif
! do not modify axis variables MU(C), AC(C), LNAC(C) !!!
             endif
          else
! this the vertical axis of a STEP calculation, most often T as axis 1
! maybe change default for iax=1 also.  Most frequent vertical axis is NP(*)
             if(iax.eq.2) text='NP(*)'
             if(stepspecial(1)) then
! step separate, default vertical axis is GM, horizontal fraction
                if(iax.eq.2) text='GM(*)'
             elseif(stepspecial(2)) then
! Scheil, PFL (Phase Fraction Liquid) or PFS are special plot functions
                if(iax.eq.1) text='PFL'
                if(iax.eq.2) text='T'
             elseif(stepspecial(3)) then
! Tzero, fraction vs T
                if(iax.eq.2) text='w(c)'
             elseif(stepspecial(4)) then
! Paraequilibrium, fraction vs T
                if(iax.eq.2) text='T'
             elseif(stepspecial(5)) then
! step liquid_eet
                if(iax.eq.2) text='T'
             endif
             nullify(maptop%plotlink)
          endif iaxval
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
! Here the user specifies his axis for plotting
21000   continue
       if(iax.eq.1) then
          call gparcdx('Horizontal axis variable',&
               cline,last,7,axplot(iax),axplotdef(iax),'?Plot command')
! Note "7" means that a "," inside x(liq,fe) will not return just "x(liq"
       else
          call gparcdx('Vertical axis variable',&
               cline,last,7,axplot(iax),axplotdef(iax),'?Plot command')
       endif
       if(buperr.ne.0) then
          jump=990; return
       endif
! extract a possible scaling factor like 0.001*GM(*)
       jp=1
       call getrel(axplot(iax),jp,xxx)
       if(buperr.eq.0) then
! there is a numerical factor
          graphopt%scalefact(iax)=xxx
! a number must be followed by a *
          if(axplot(iax)(jp:jp).ne.'*') then
             write(*,*)'Scaling factor must be followed by *'
             jump=990; return
          else
! Fortran allows overlapping strings in assignments
             axplot(iax)=axplot(iax)(jp+1:)
          endif
       else
! no scaling factor, graphopt%scalfactor(iax) already unity
          buperr=0
       endif
       if(index(axplot(iax),'*').gt.0 .or. index(axplot(iax),'#').gt.0) then
             wildcard=.TRUE.
       endif
       if(axplotdef(iax).ne.axplot(iax)) then
! if new axis variable then reset default plot options
! plot ranges and their defaults
          call reset_plotoptions(graphopt,plotfile,textlabel)
! check that axis variable is a correct state variable or symbol
! Most code copied from show variable (case(4,17) around line 3273)
! Avoid capson of axplot(iax) for possible other problems later
          name1=axplot(iax)
          call capson(name1)
          if(name1(1:4).eq.'PFL ' .or. name1(1:4).eq.'PFS ') then
! this is a special function allowed in Scheil simulations for phase frac liq
             if(.not.stepspecial(2)) then
  write(*,*)'The PFL and PFS functions are allowed only for Scheil simulations'
                jump=100; return
             endif
          elseif(index(axplot(iax),'*').gt.0) then
! generate many values
! the values are returned in yarr with dimension maxconst. 
! longstring are the state variable symbols for the values ...
             call get_many_svar(axplot(iax),yarr,maxconst,i1,longstring,ceq)
             if(gx%bmperr.ne.0) then
! if error go back to command level
                write(kou,*)'Illegal axis variable!  Error code: ',gx%bmperr
                jump=100; return
             endif
          elseif(index(axplot(iax),'#').gt.0) then
! generate many values including for metastable phases
! the values are returned in yarr with dimension maxconst. 
! longstring are the state variable symbols for the values ...
             call get_many_svar(axplot(iax),yarr,maxconst,i1,longstring,ceq)
             if(gx%bmperr.ne.0) then
! if error go back to command level
                write(kou,*)'Illegal axis variable!  Error code: ',gx%bmperr
                jump=100; return
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
                call find_svfun(axplot(iax),istv)
                if(gx%bmperr.ne.0) then
                   write(kou,*)'Illegal axis variable, error: ',gx%bmperr
                   jump=100; return
                endif
             endif
          endif
       endif
! remember most recent axis as default (and to avoid reset)
       axplotdef(iax)=axplot(iax)
    enddo pltaxdef
    return
  end subroutine ui_plot_axes

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine ui_delete_step_map(kind,cline,last,ceq)
! remove previous STEP/MAP results before a new STEP/MAP calculation
! kind=0 delete without asking (STEP TZERO and SCHEIL_GULLIVER)
! kind=1 STEP, ask, previous results can be kept
! kind=2 MAP, ask, previous results can be kept
    integer kind
    character cline*(*)
    integer last
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character ch1*1
    if(kind.eq.0) then
       if(associated(maptop)) then
          write(kou,*)'Previous map/step results will be deleted'
          call delete_mapresults(maptop)
       endif
    else
       if(.not.associated(maptop)) return
       write(kou,833)
833    format('There are previous results from step or map')
       if(kind.eq.1) then
          call gparcdx('Delete them?',cline,last,1,ch1,'Y','?Step old data')
       else
          call gparcdx('Reinitiate?',cline,last,1,ch1,'Y','?Map old data')
       endif
       if(.not.(ch1.eq.'y' .or. ch1.eq.'Y')) then
! keep previous results, start indexing new nodes/lines from previous
          if(kind.eq.1) then
! for step separate it seems difficult to have correct seqx !!
             seqxyz(1)=max(maptop%next%seqx,maptop%previous%seqx,maptop%seqx)
          else
             seqxyz(1)=maptop%next%seqx
          endif
          seqxyz(2)=maptop%seqy
!             seqxyz(3) can be used for something else ...
          maptopsave=>maptop
          nullify(maptop)
          if(kind.eq.1) write(*,'(a,2i4)')'Previous results kept',seqxyz
          return
       endif
! there should be a more careful deallocation to free memory
       call delete_mapresults(maptop)
    endif
    nullify(maptop)
    nullify(maptopsave)
    if(kind.eq.1) write(kou,*)'Previous results removed'
! This is to keep trace of total number of saved equilibria
    totalsavedceq=0
! delete equilibria associated with STEP/MAP (for MAP done by delete_mapresults)
    if(kind.eq.1) call delete_equilibria('_MAP*',ceq)
! initiate indexing nodes and lines
    seqxyz=0
! remove all graphopt settings
    call reset_plotoptions(graphopt,plotfile,textlabel)
    axplotdef=' '
    return
  end subroutine ui_delete_step_map

