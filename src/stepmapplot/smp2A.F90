! These soubroutine calculate the diagram, smp2B plot it

!\addtotable subroutine map_setup
!\begin{verbatim}
  subroutine map_setup(maptop,nax,axarr,seqxyz,starteqs)
! main map/step routine
! THIS HAS BEEN SPLIT IN TWO PARTS
! This first part tranforms all user provided or automatic start points
! to start equilibria
! The second goes through the list of start equiliria until it is
! empty
!
! maptop is the main map_node record which will return all calculated lines.
! nax is the number of axis (can be just one for STEP)
! axarr is an array of records specifying the axis for the step/map
! seqxyz are intial values for number of nodes and lines
! starteqs is an array with equilibrium data record
! they are linked using the ceq%next index
    implicit none
    integer nax,seqxyz(*)
    type(map_axis), dimension(nax) :: axarr
!    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(starteqlista), dimension(*) :: starteqs
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq,starteq
    type(gtp_condition), pointer :: pcond
    TYPE(map_node), pointer :: tmp
    type(map_line), pointer :: mapline
! should this meqrec be a pointer or not??
    type(meq_setup), pointer :: meqrec
    type(map_fixph), allocatable :: mapfix
!    type(map_fixph), pointer :: mapfix
    double precision starting,finish2,axvalok,dgm,tsave,xxx,yyy,zzz
    integer starttid,endoftime,bytdir,seqz,nrestore,termerr,lastimethiserror
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(gtp_state_variable), target :: svrtarget

! save current conditions
    character savedconditions*1024
! for saving a copy of constitutions
    double precision, allocatable, dimension(:) :: copyofconst
! inactive are indices of axis conditions inactivated by phases set fixed
! inactive not used ...
    integer iadd,irem,isp,seqx,seqy,mode,halfstep,jj,ij,inactive(4),bytaxis
    integer ceqlista
! inmap=1 turns off converge control of T
    integer, parameter :: inmap=1
    character ch1*1
    logical firststep,onetime
!
!    write(*,*)'in map_setup'
! save all conditions 
!    call get_all_conditions(savedconditions,-1,starteqs(1)%p1)
    ij=1
    savedconditions=' '
    savecond: do jj=1,nax
!       write(*,*)'SMP2A get_one: ',ij,axarr(jj)%seqz
       call get_one_condition(ij,savedconditions,&
            axarr(jj)%seqz,starteqs(1)%p1)
       if(gx%bmperr.ne.0) then
          gx%bmperr=0; savedconditions=' '; exit savecond
       endif
       ij=len_trim(savedconditions)+2
    enddo savecond
!    write(kou,*)'SMP2A saved: ',trim(savedconditions)
    nrestore=0
    lastimethiserror=0
! first transform start points to start equilibria on zero phase lines
! All axis conditions except one are converted to fix phase conditions 
! (if there is just one axis skip this)
! One or more map_node records are created with mapline records each
    call cpu_time(starting)
    call system_clock(count=starttid)
    inactive=0
!
    if(ocv()) write(*,*)'Entering map_setup',nax
! if automatic statpoints requested they are generatet here
!    call auto_startpoints(maptop,nax,axarr,seqxyz,starteq)
!    ceq=>starteq
    ceq=>starteqs(1)%p1
    iadd=1
!    ceqlista=1
21  continue
!    write(*,'(a,a,3i4)')'SMP2A Start equilibrium: ',trim(ceq%eqname),&
!         ceq%eqno,ceq%nexteq,ceq%multiuse
!    if(ceq%nexteq.gt.0) then
!       ceq=>eqlista(ceq%nexteq)
!       iadd=iadd+1
!       goto 21
!    endif
! noofstarteq is a global variable in SMP, set by calling routine
    if(noofstarteq.gt.0) write(*,*)'There are ',noofstarteq,' start equilibria'
! loop to change all start equilibria to start points
! Store the start points in map_node records started from maptop
    do ceqlista=1,noofstarteq
       ceq=>starteqs(ceqlista)%p1
!       write(*,*)'SMP2A calling map_startpoint: ',trim(ceq%eqname),ceq%eqno
!       read(*,106)ch1
106    format(a)
! convert all axis conditions except one to fix phase
       call map_startpoint(maptop,nax,axarr,seqxyz,inactive,ceq)
       if(gx%bmperr.ne.0) then
          write(*,101)ceq%nexteq,gx%bmperr
101       format('Failed calculate a start point: ',i4,i7)
!             ceq=>eqlista(ceq%nexteq)
          gx%bmperr=0
          goto 900
       endif
! I have not really implemented several startpoint, I am not sure
! if each does each have separate maptop and savesec ....
! error if no startpoints 
       if(.not.associated(maptop)) then
          write(*,*)'Cound not find a single start equilibria for',ceqlista
!          gx%bmperr=4224; goto 1100
          goto 900
       endif
!       write(*,*)'There is a MAPTOP record ...'
! create array of equilibrium records for saving results
       seqy=maxsavedceq
       call create_saveceq(maptop%saveceq,seqy)
       if(gx%bmperr.ne.0) goto 1000
! initiate node counter done, line counter will be incremented
       if(maptop%seqx.gt.1) write(*,85)maptop%seqx,maptop%seqy+1
85     format('Previous step/map results saved'/&
            'New mapnode/line equilibria indices will start from: ',i3,i5)
!       maptop%seqy=0
!       write(*,*)'savesize: ',size(maptop%saveceq%savedceq)
! if there are more startpoints try to convert these to start equilibria
900    continue
!       write(*,*)'At label 900: ',gx%bmperr
    enddo
!    write(*,*)'SMP Finished loop',associated(maptop)
    if(associated(maptop)) then
       if(allocated(maptop%linehead)) then
! Clear any error code if we have linhead allocated
          if(gx%bmperr.ne.0) gx%bmperr=0
       else
          write(*,*)'Failed to find any lines to calculate'
          goto 1000
       endif
    else
! no maptop record
       write(*,*)'Failed finding startpoints for step/map'
       goto 1100
    endif
!-----------------------------------------------------
! now we should calculate all lines stored as start equilibria       
! but maybe there are no start equilibria??
! starteq is a ceq record, mapping will use maptop record ....
    write(*,*)'SMP2A call map_doallines'
    call map_doallines(maptop,nax,axarr,seqxyz,starteq)
!    write(*,*)'SMP2A back from map_doallines'
!-----------------------------------------------------
1000 continue
!--------------------------------------------------
! Here we have now finished the step/map.
! Set back inactive axis conditions How??
!    do ij=2,inactive(1)
!       call locate_condition(inactive(ij),pcond,ceq)
!       pcond%active=0
!    enddo
    call system_clock(count=endoftime)
    call cpu_time(finish2)
    if(gx%bmperr.ne.0) then
       write(*,1005)gx%bmperr
1005   format('STEP/MAP terminated with error code: ',i5)
       gx%bmperr=0
    else
       write(*,1010)maptop%saveceq%free-1,finish2-starting,endoftime-starttid
1010   format(/'Finished step/map with ',i5,' equilibria in ',&
            1pe12.4,' CPU s and ',i7,' cc')
    endif
    if(len_trim(savedconditions).gt.0) then
!       write(*,*)'SMP2A restore: ',trim(savedconditions)
!       if(index(savedconditions,'>=').gt.0) then
! conditions including a fix phase, do not try to restore 
!          write(*,*)'SMP2A cannot restore original conditions'
!          goto 1100
!       endif
!       write(*,*)'Restoring all initial conditions: '
!       write(*,*)trim(savedconditions)
! ij is incremented by 1 inside set_condition
       ij=0
! SUCK, I fixed that conditions with 2 terms was not entered again but
!       after other changes to handle condition with species such as O-2
!       the same problem!  Just remove all conditions and set those saved!!
!       write(*,*)'SMP2A conditions at end of step/map'
!       It may create loss of memory but ... what the heck ... buy more!
!       call list_conditions(kou,ceq)
!       write(*,*)'SMP2A remove all conditions'
!       if(nax.eq.1) then
!          write(*,*)'SMP2 Conditions can be changed by some STEP commands'
!       endif
!       goto 1100
!-----------------------------------------------------
! I am not sure it is critical to restore conditions ...
! it could be some cases when conditions are modified in STEP TZERO/SCHEIL/PARA
!-----------------------------------------------------
! this does not work because axis and maybe other things refer to
! conditions by index.  If I remove all condtions to restore them
! these indices become invalid    
!       call set_condition('*:=none ',ij,starteqs(1)%p1)
!       call list_conditions(kou,ceq)
       ij=0
!       write(*,*)'SMP2A restore axis cond: ',trim(savedconditions)
       call set_condition(savedconditions,ij,starteqs(1)%p1)
       if(gx%bmperr.ne.0) write(*,*)'Error restoring axis conditions',gx%bmperr
!       write(*,*)'SMP2A restored conditions:'
!       call list_conditions(kou,ceq)
    else
       write(*,*)'SMP2A axis conditions could not be restored'
    endif
1100 continue
    return
  end subroutine map_setup
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_doallines
!\begin{verbatim}
  subroutine map_doallines(maptop,nax,axarr,seqxyz,starteq)
! main map/step routine
! maptop is the main map_node record which will return all calculated lines.
! nax is the number of axis (can be just one for STEP)
! axarr is an array of records specifying the axis for the step/map
! seqxyz are intial values for number of nodes and lines
! starteq is an equilibrium data record, if there are more start equilibria
! they are linked using the ceq%next index
    implicit none
    integer nax,seqxyz(*)
    type(map_axis), dimension(nax) :: axarr
    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq
    type(gtp_condition), pointer :: pcond
    TYPE(map_node), pointer :: tmp
    type(map_line), pointer :: mapline
! should this meqrec be a pointer or not??
    type(meq_setup), pointer :: meqrec
    type(map_fixph), allocatable :: mapfix
!    type(map_fixph), pointer :: mapfix
    double precision starting,finish2,axvalok,dgm,tsave,xxx,yyy,zzz,axval
    integer starttid,endoftime,bytdir,seqz,nrestore,termerr,lastimethiserror
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(gtp_state_variable), target :: svrtarget
! save current conditions
    character savedconditions*1024
! for saving a copy of constitutions
    double precision, allocatable, dimension(:) :: copyofconst
! inactive are indices of axis conditions inactivated by phases set fixed
! inactive not used ...
    integer iadd,irem,isp,seqx,seqy,mode,halfstep,jj,ij,inactive(4),bytaxis
    integer ceqlista,phfix,haha,lastax,mapx,lokph,lokcs,bypass
! inmap=1 turns off converge control of T
    integer, parameter :: inmap=1
    character ch1*1
    logical firststep,onetime,noderrmess
!
!    write(*,*)'in map_doallines'
!-------------------------------
! return here for each new line to be calculated
! NOTE we can start a new thread for each line, when a node is found
! all threads stop.  
! If the node already exists the exit corresponing to the new line removed
! and the thread ends
! initiate phfix, looking for crash it seems to be used before set ...
    phfix=0
! If the node is new it is created and exits added and the thread ends.
    inactive=0
    nrestore=0
    lastimethiserror=0
300 continue
! this is to write a warning message once for each line
    onetime=.true.   
    bytaxis=0
    firststep=.TRUE.
! THREADPROTECTED CALL the map_findline will copy the ceq from mapnode
    if(ocv()) write(*,*)'Looking for a line to calculate'
    call map_findline(maptop,axarr,mapfix,mapline)
    if(gx%bmperr.ne.0) goto 1000
! if no line we are finished!
!   write(*,*)'Back from map_findline 1: ',associated(mapline),allocated(mapfix)
! segmentation fault crash later ...
    if(.not.associated(mapline)) goto 900
!    write(*,*)'We will start calculate line: ',mapline%lineid,mapline%axandir
    if(maptop%tieline_inplane.ne.0) then
! for mapping we need to check how all axis varies
       allocate(mapline%axvals(nax))
       allocate(mapline%axvalx(nax))
       if(maptop%tieline_inplane.gt.0) then
! with tie-lines in plane we must check axis variable for stable phase also
          allocate(mapline%axvals2(nax))
!       else
! any special  to do??          
       endif
    endif
! Each thread must have separate meqrec and ceq records
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    ceq=>mapline%lineceq
! ?? We may have incompatibility between ceq and meqrec if new compsets added
! maybe meqrec should not be a pointer?
    meqrec=>mapline%meqrec
    noderrmess=.true.
! No grid minimization and the phr is not deallocated with mode<0
! It is necessary to generate new meqrec for each line as there may be new
! composition sets created in other threads.  But we must also specify 
! phases set fix due to the mapping to replace axis conditions.  
! We must provide an array of phase tuples with fix phases. 
!    write(*,*)'Calling calceq7 for new line: ',mapline%lineid
    if(ocv()) write(*,*)'Calling calceq7 for new line: ',mapline%lineid
! mode=-1 means no gridminimization and do not deallocate phr
    mapline%problems=0
    mapline%lasterr=0
    mode=-1
    if(ocv()) write(*,*)'This call generates mapline%meqrec for this line'
    bytdir=0
! the save constitutions may be useful if problems ... ???
    if(allocated(copyofconst)) deallocate(copyofconst)
! segmentation fault in this subroutine ...
! because I checked only size(..) and not if it was allocated ...
    call save_constitutions(ceq,copyofconst)
! segmentation fault before this output ...
!    write(*,*)'called save_constitutions: ',size(copyofconst)
305 continue
! to be able to handle problems copy the constitutions!!
!    if(mapline%problems.gt.0) then
!       write(*,*)'problems',mapline%problems,ceq%tpval(1)
!    endif
!    write(*,*)'Calling calceq7 with T=',ceq%tpval(1),mapline%axandir
!    write(*,*)'Calling calceq7 with meqrec%status:',meqrec%status
    call calceq7(mode,meqrec,mapfix,ceq)
!    write(*,*)'SMP2A Back from calceq7 ',gx%bmperr,meqrec%status
    if(gx%bmperr.ne.0) then
! error 4187 is to set T or P to less than 0.1
       if(gx%bmperr.eq.4187) then
          goto 306
       endif
       if(mapline%number_of_equilibria.eq.0) then
! We can add/subtract a small amount of axis condition if error at first step
!          write(*,*)'Error at first equilibrium: ',gx%bmperr,mapline%axandir
          mapline%lasterr=gx%bmperr
          mapline%problems=mapline%problems+1
!          if(bytdir.eq.1) then
! we have tried adding a small step in axandir direction, now change direction
!             mapline%axandir=-mapline%axandir
!          elseif(bytdir.gt.1) then
! give up
!             goto 306
!          endif
! Extract the current value of the axis state variable items using pcond
          jj=abs(mapline%axandir)
!          write(*,*)'SMP: axandir: ',jj,gx%bmperr
          gx%bmperr=0
          if(jj.le.0 .or. jj.gt.2) then
             write(*,*)'SMP error: no axis direction! Set to 1'
             mapline%axandir=1
             jj=1
!             call list_conditions(kou,ceq)
          endif
          seqz=axarr(jj)%seqz
          call locate_condition(seqz,pcond,ceq)
          if(gx%bmperr.ne.0) goto 1000
          call condition_value(1,pcond,zzz,ceq)
          if(gx%bmperr.ne.0) goto 1000
          if(ocv())write(*,765)bytdir,jj,zzz,mapline%axvals(jj),axarr(jj)%axinc
765       format('Attempt to step 1: ',i2,i3,3(1pe16.8))
! first time bytdir=1, second time bytdir=2, compensate for first step ...
!          yyy=1.0D-2*bytdir*axarr(jj)%axinc
!          yyy=1.0D-3*bytdir*axarr(jj)%axinc
!          xxx=zzz+mapline%axandir*yyy
          xxx=zzz
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! restore constitutions, not a good idea ?? ...
!          write(*,*)'Restore constitutions 1'
          call restore_constitutions(ceq,copyofconst)
!
          call map_problems(maptop,mapline,axarr,xxx,1)
          if(gx%bmperr.ne.0) goto 306
          if(ocv()) write(*,737)'Error at first step: ',mapline%axandir,&
               mapline%nodfixph,zzz,xxx
737       format(a,2i3,6(1pe14.6))
!          read(*,738)ch1
738       format(a)
! set the condition value ... ???
          if(nax.gt.1) then
! run time error that axvals has dimension 0 ... when step
             mapline%axvals(abs(jj))=xxx
          endif
          call condition_value(0,pcond,xxx,ceq)
          if(gx%bmperr.ne.0) goto 1000
          if(ocv()) write(*,765)0,mapline%axandir,zzz,xxx,yyy
! call calceq7 again, we must deallocate meqrec%phr 
          deallocate(meqrec%phr)
          goto 305
       endif
306    continue
!       write(*,*)'SMP2 Generating mapline%meqrec failed 2: ',gx%bmperr
       call map_lineend(mapline,axarr(abs(mapline%axandir))%lastaxval,ceq)
! look for a new line to follow
       goto 300
    endif
!    write(*,*)'back from calceq7B'
! if all has gone well deallocate mapfix
    if(allocated(mapfix)) deallocate(mapfix)
!    write(*,*)'SMP successfully deallocated mapfix'
!--------------------------------
! limit the maximum change in T and P, should be small during step/map
    meqrec%tpmaxdelta(1)=2.0D1
    meqrec%tpmaxdelta(2)=1.0D1
    bypass=0
!--------------------------------
! return to label 310 after each new equilibrium calculated along the same line
! Follow the equilibria along a line.  For each equilibria calculated
! store the data.  If the phase set want to change (irem or iadd>0) calculate
! exactly the phase change, generate a node and terminate the line and then
! look for a new line to follow.
310 continue
    halfstep=0
! save current value of T if trouble later ...
    tsave=ceq%tpval(1)
! try saving constitutions ...
    if(allocated(copyofconst)) deallocate(copyofconst)
    call save_constitutions(ceq,copyofconst)
! emergency return when two phases want to change status
320 continue
    iadd=0
! Note setting iadd=-1 turn on verbose inside meq_sameset
321 continue
    irem=0
    mapline%meqrec%noofits=0
!    write(*,*)'Calling meq_sameset 7: ',mapline%number_of_equilibria,&
!         ceq%tpval(1),gx%bmperr
!
!    call list_conditions(kou,ceq)
!
!    write(*,*)'Calling meq_sameset ',mapline%more,mapline%number_of_equilibria
!    write(*,884)1,mapline%linefixph(1)%ixphase,&
!         mapline%linefixph(1)%compset,iadd,meqrec%nphase,abs(phfix)
!884 format('SMP fix phase ',i1,':',i3,i2,', new fix phase: ',i3,&
!            ', number of phases: ',i3,' abs(phfix): ',i3)
!--------------------------------------------------------------------------
! This is where most equilibrium calculations are made
!--------------------------------------------------------------------------
!
!    write(*,*)'smp2A calling meq_sameset from map_doallines',ceq%tpval(1)
    call meq_sameset(irem,iadd,mapx,mapline%meqrec,mapline%meqrec%phr,inmap,ceq)
!
!--------------------------------------------------------------------------
!    write(*,331)'SMP Back from meq_sameset ',mapline%number_of_equilibria,&
!         irem,iadd,gx%bmperr,phfix,ceq%tpval(1),ceq%phase_varres(4)%dgm
331 format(a,5i5,2(F10.2))
!    write(*,884)2,mapline%linefixph(1)%ixphase,&
!         mapline%linefixph(1)%compset,iadd,meqrec%nphase,abs(phfix)
!------------------------------------------------------------------
! we come back here if iadd was 0 but removed as 
3000 continue
! new global check for stable and metastable phases
!    write(*,*)'SMP error 6A:',mapline%number_of_equilibria,&
!         maptop%globalcheckinterval
    phasecheck: if(gx%bmperr.eq.0 .and. iadd.eq.0 .and. irem.eq.0) then
!       write(*,*)'SMP error 6B:',mapline%number_of_equilibria,&
!            maptop%globalcheckinterval
!       if(maptop%globalcheckinterval.le.0) then
!        write(*,*)'SMP maptop%globalcheckinterval:',maptop%globalcheckinterval
!          maptop%globalcheckinterval=10
!       endif
       checkinterval: if(maptop%globalcheckinterval.gt.0) then
          if(mod(mapline%number_of_equilibria,maptop%globalcheckinterval).eq.0)&
               then
! this may set error code if equilibrium should be recalculated
! and it may change constitutions of metastable phases
!          write(*,'(a,i5)')'SMP check_all_phases at equilibrium: ',&
!               mapline%number_of_equilibria
             jj=0
             call check_all_phases(jj,ceq)
             if(gx%bmperr.ne.0) then
!             if(associated(mapline%lineceq,ceq)) then
! This is true and dangerous but I will be careful programming ...
!                write(*,*)'SMP ceq is same as mapline%lineceq'
!             else
!                write(*,*)'SMP ceq is NOT same as mapline%lineceq'
!             endif
!             call get_phase_compset(iph,ics,lokph,lokres)
                if(gx%bmperr.eq.4366) then
! terminate line and call gridminimizer
!                   write(*,*)'SMP check_all_phases require gridminimizer',jj
                   gx%bmperr=0
                   call map_halfstep(halfstep,0,axvalok,mapline,axarr,ceq)
                   if(gx%bmperr.eq.0) goto 321
                elseif(gx%bmperr.eq.4365) then
!                write(*,*)'SMP check_all_phases error, call map_halfstep:',jj
                   gx%bmperr=0
! we have to convert jj=iph*10+ics to index in mapline%meqrec%phr
! Check if constitution is the one se in check_all_phases
!                write(*,95)(yarr(ii),ii=1,jj)
!95              format('3Y gridy: ',10(F7.4))
                   call map_halfstep(halfstep,0,axvalok,mapline,axarr,ceq)
                   if(gx%bmperr.eq.0) goto 321
                endif
! otherwise ignore any errors
                gx%bmperr=0
             endif
          endif
       endif checkinterval
    endif phasecheck
!------------------------------------------------------------------
!    write(*,*)'SMP looking for error 7:'
    sameseterror: if(gx%bmperr.ne.0) then
!       write(*,*)'Error in meq_sameset called from smp',gx%bmperr
! if error 4359 (slow convergence), 4204 (too many its) take smaller step ...
! error 4195 means negative phase amounts
491    continue
       if(gx%bmperr.eq.4195 .or. gx%bmperr.eq.4359 &
            .or. gx%bmperr.eq.4204) then
! I am not sure there is really any change for the equilibrium calculated ...
!          write(*,317)'Trying half step: ',halfstep,mapline%axandir,&
!               mapline%number_of_equilibria,lastimethiserror,ceq%tpval(1)
317       format(a,2i3,2i4,f9.2)
          if(mapline%number_of_equilibria-lastimethiserror.gt.10) then
             lastimethiserror=mapline%number_of_equilibria
!          if(mapline%meqrec%noofits-lastimethiserror.gt.10) then
!             lastimethiserror=mapline%meqrec%noofits
             gx%bmperr=0
             mapline%axfact=1.0D-2
             call map_halfstep(halfstep,0,axvalok,mapline,axarr,ceq)
!             write(*,*)'Back from halfstep 1',halfstep,gx%bmperr
             if(gx%bmperr.eq.0) goto 321
          endif
       elseif(gx%bmperr.eq.4364) then
! Two stoichiometric phases with same composition stable, we have
! to calculated an invariant equilibrium T in a different way.
! if tielines in plane create nodepoint otherwise I do not know what to do
          if(maptop%tieline_inplane.gt.0) then
! dummy values (for the moment)
             axval=ceq%tpval(1)
             haha=mapx
             phfix=iadd
             lastax=abs(mapline%axandir)
! maybe save last calculated equilibrium as endpoint of current line?
! Collecting values needed fot map_newnode
! irem is last fix phase, haha is entered phase, phfix is new stable phase
!             write(*,219)'SMP call map_newnode: ',lastax,axval,&
!                  ceq%tpval(1),meqrec%nstph,irem,haha,phfix
!219          format(a,i2,2F12.4,5i5)
! list current settings: is content of mapline%meqrec same as meqrec??
!             write(*,885)mapline%nfixphases,&
!                  mapline%linefixph(1)%ixphase,mapline%linefixph(1)%compset,&
!                  mapline%meqrec%nv,mapline%meqrec%iphl(1),&
!                  mapline%meqrec%iphl(2)
!885          format('SMP Fixed phase:',i2,': ',i3,i2,', entered: ',i2,': ',5i3)
             mapline%status=ibset(mapline%status,TWOSTOICH)
             call map_newnode(mapline,meqrec,maptop,axval,lastax,axarr,&
                  phfix,haha,ceq)
             if(gx%bmperr.ne.0) then
! give up on this line, map_lineend set error code to zero
                write(*,*)'Failed create node point, terminate and take next',&
                     gx%bmperr
                call map_lineend(mapline,axvalok,ceq)
                axvalok=zero
             endif
          endif
       endif
! give up this line, reset error code and check if there are more lines
       gx%bmperr=0
       goto 805
    endif sameseterror
!    write(*,323)'Calc line: ',gx%bmperr,irem,iadd,mapline%axandir,&
!         mapline%meqrec%noofits,mapline%meqrec%nstph,ceq%tpval(1)
    if(ocv())write(*,323)'Calc line: ',gx%bmperr,irem,iadd,mapline%axandir,&
         mapline%meqrec%noofits,mapline%meqrec%nstph,ceq%tpval(1)
323 format(a,i5,2i3,2i4,i3,f10.2)
    if(iadd.gt.0) then
! check if it is a closing miscibility gap or loss of ordering
! remove iadd if it is a phase with same composition as an already stable one
       if(same_composition(iadd,mapline%meqrec%phr,mapline%meqrec,ceq,dgm)) &
            iadd=0
    endif
!    write(*,*)'Check if same phase: ',iadd
330 continue
    if(gx%bmperr.eq.0 .and. irem.eq.0 .and. iadd.eq.0) then
! no error and no change of phase set, just store the calculated equilibrium.
! and calculate another point along the line
!       write(*,*)'hms: Storing equilibrium',&
!            mapline%number_of_equilibria,maptop%globalcheckinterval
       if(mapline%number_of_equilibria.gt.10 .and. mapline%nodfixph.gt.0) then
! we have managed 3 steps, set phase at start node as entered (if dormant)
          if(meqrec%phr(mapline%nodfixph)%phasestatus.eq.PHDORM) then
!             write(*,*)'Phase set entered ',mapline%nodfixph
             meqrec%phr(mapline%nodfixph)%phasestatus=PHENTUNST
          endif
       endif
!       mapline%problems=0
!       nrestore=0
       call map_store(mapline,axarr,nax,maptop%saveceq)
       if(gx%bmperr.ne.0 .or. mapline%more.eq.0) then
! Test if we are running out of memory 
          if(gx%bmperr.eq.4219) goto 1000
          if(gx%bmperr.eq.4360) then
! too big difference in some axis, take halfstep
!             write(*,*)'Take a half step',halfstep
             gx%bmperr=0; halfstep=halfstep+1
             call map_halfstep(halfstep,0,axvalok,mapline,axarr,ceq)
             if(gx%bmperr.eq.0) goto 321
          endif
! terminate line any error code will be cleared inside map_lineend. 
!          write(*,*)'Calling map_lineend 1'
          call map_lineend(mapline,axarr(abs(mapline%axandir))%lastaxval,ceq)
          goto 300
       endif
! stored last calculated equilibrium 
       mapline%problems=0
       nrestore=0
! check which axis variable changes most rapidly, maybe change step axis
! (for tie-lines in plane check axis values for all phases)
! and take a step in this axis variable making sure it inside the limits
! and continue, else terminate and take another start equilibrium
! Normally do not change the phase kept fix.
!       write(*,*)'hms: taking a step'
       call map_step(maptop,mapline,mapline%meqrec,mapline%meqrec%phr,&
            axvalok,nax,axarr,ceq)
!       write(*,*)'Back from map_step 1',mapline%more,&
!            mapline%number_of_equilibria,gx%bmperr
       if(gx%bmperr.ne.0) then
!          write(*,*)'SMP2A error return from map_step 1: ',gx%bmperr
          gx%bmperr=0
          if(meqrec%tpindep(1)) then
!             write(*,*)'SMP2A restore T 1: ',tsave,axvalok
             ceq%tpval(1)=tsave
          endif
          call map_halfstep(halfstep,0,axvalok,mapline,axarr,ceq)
          if(gx%bmperr.eq.0) then
! jump back without setting halfstep=0, setting iadd=-1 turn on debug output 
!          iadd=-1
             goto 321
          endif
       endif
! if mapline%more>0 continue, otherwise line has terminated at axis limit
! check if there are other nodes with lines to calculate
!       write(*,*)'Back from step:',gx%bmperr,mapline%more,ceq%tpval(1)
! if mapline%more>=0 there is no error and a new equilibrium to calculate
! if mapline%more<0 the line has ended at axis limit or there is an error
       if(mapline%more.ge.0) goto 310
       if(gx%bmperr.ne.0) then
!          write(*,*)'SMP2A Error stepping to next equilibria, ',gx%bmperr
       endif
! any error code will be cleared inside map_lineend.
!       write(*,*)'Calling map_lineend 1'
       call map_lineend(mapline,axarr(abs(mapline%axandir))%lastaxval,ceq)
! look for a new line to follow
       goto 300
! finish thread started at label 300 ??
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    elseif(gx%bmperr.ne.0) then
!       write(*,*)'Error return from meq_sameset: ',gx%bmperr,mapline%lasterr,&
!            ceq%tpval(1)
       termerr=gx%bmperr
       gx%bmperr=0
       if(meqrec%tpindep(1)) then
          if(ocv()) write(*,*)'Restoring T 2: ',tsave,axvalok
          ceq%tpval(1)=tsave
       endif
! also restore constitutions
       nrestore=nrestore+1
       if(nrestore.lt.3) then
!          write(*,*)'Restore constitutions 2',nrestore
          call restore_constitutions(ceq,copyofconst)
!
! take smaller steps!
          mapline%axfact=1.0D-2
!          write(*,552)'Call halfstep: ',bytaxis,nrestore,&
!               mapline%number_of_equilibria,axvalok
552       format(a,3i3,2(1pe12.4))
          call map_halfstep(halfstep,0,axvalok,mapline,axarr,ceq)
          if(gx%bmperr.eq.0) then
! jump back without setting halfstep=0, setting iadd=-1 turns on debug output 
!          iadd=-1
             goto 321
          endif
       elseif(nax.gt.1 .and. bytaxis.eq.0) then
!          write(*,*)'Restore last OK: ',mapline%number_of_equilibria,nrestore,&
!               axvalok
          call restore_constitutions(ceq,copyofconst)
          if(meqrec%tpindep(1)) then
             if(ocv()) write(*,*)'Restoring T 3: ',tsave,axvalok
             ceq%tpval(1)=tsave
          endif
          if(ocv()) write(*,555)'Repeated error 7, try to change axis',&
               gx%bmperr,ceq%tpval(1),axvalok,tsave
555       format(a,i5,3F8.2)
          gx%bmperr=0
          bytaxis=1
! Make sure that the current axis has the last successfully calculated value
! as prescribed value
          call locate_condition(axarr(abs(mapline%axandir))%seqz,pcond,ceq)
          if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to extract the value, 0 means to set the value
          call condition_value(1,pcond,xxx,ceq)
          if(gx%bmperr.ne.0) goto 1000
          call condition_value(0,pcond,axvalok,ceq)
          if(gx%bmperr.ne.0) goto 1000
          write(*,19)'Force changeaxis: ',mapline%axandir,gx%bmperr,axvalok,xxx
19        format(a,i3,i5,2(1pe14.6))
!
          call map_force_changeaxis(maptop,mapline,mapline%meqrec,&
               nax,axarr,axvalok,ceq)
!          write(*,*)'new changeaxis: ',mapline%axandir,gx%bmperr,axvalok
!          call list_conditions(kou,ceq)
          if(gx%bmperr.eq.0) goto 320
       endif
! Giv up, terminate the line and check if there are other lines to calculate
! macro map1 ends at composition axis end still with T axis as variable !!
!       write(*,*)'Calling map_lineend 3',nrestore,termerr
       gx%bmperr=termerr
       call map_lineend(mapline,axarr(abs(mapline%axandir))%lastaxval,ceq)
! find a new line
       goto 300
    endif
!------------------------------------------------------------
379 continue
    phasechange: if(irem.gt.0 .and. iadd.gt.0) then
! We can also have a stoichiometic phase with ALLOTROPIC transformation
! which will change form one to another at a fix T
       if(allotropes(irem,iadd,meqrec%noofits,ceq)) then
          irem=0
          goto 379
       endif
! if there is phase which wants to appear and another disappear then
! first check if they are the composition sets of the same phase
! calculate with half the step 5 times. If axvalok=0 no previous axis value
! BUG: Problems here for map5.OCM, when matsmin compiled with -O2
! two extra composition sets of BCC and LIQUID wanted to appear.
!  Will lok at that later ...
       if(onetime) then
!          write(*,22)'SMP: phases appear and disappear at same time: ',&
!               iadd,irem,phasetuple(iadd)%lokph,phasetuple(irem)%lokph
22        format(a,4i4)
          onetime=.false.
       endif
!       write(*,*)
! restore constitutions
!       write(*,*)'Restore constitutions 3',halfstep,axvalok,ceq%tpval(1)
       call restore_constitutions(ceq,copyofconst)
       call map_halfstep(halfstep,1,axvalok,mapline,axarr,ceq)
       if(gx%bmperr.eq.0) then
! jump back without setting halfstep=0
          goto 320
       elseif(nax.gt.1 .and. bytaxis.eq.0) then
! try to change axis with active condition.
          if(meqrec%tpindep(1)) then
             if(ocv()) write(*,*)'Restoring T 4: ',tsave,axvalok
             ceq%tpval(1)=tsave
          endif
          write(*,557)gx%bmperr,ceq%tpval(1),axvalok
557       format('Repeated error 8, try to change axis',i5,F8.2,1pe14.6)
          gx%bmperr=0
          bytaxis=1
          call map_force_changeaxis(maptop,mapline,mapline%meqrec,nax,axarr,&
               axvalok,ceq)
          if(gx%bmperr.eq.0) goto 320
          call map_lineend(mapline,axvalok,ceq)
       else
! there is an error, take another line
          call map_lineend(mapline,axvalok,ceq)
       endif
!-----------------------------------------------------
! phasechange elseif: a new phase stable or a stable wants to disappear
    elseif(irem.gt.0 .or. iadd.gt.0) then
!       write(*,*)'SMP2A new phase 2: ',iadd,irem,mapline%nodfixph,&
!            mapline%number_of_equilibria
       if(mapline%number_of_equilibria.lt.2 .and.&
            ((irem.gt.0 .and. irem.eq.mapline%nodfixph) .or. &
            (iadd.gt.0 .and. iadd.eq.mapline%nodfixph))) then
          mapline%axandir=-mapline%axandir
          write(*,*)'Ignore same phase as at startnode: ',1,mapline%nodfixph
          write(*,*)'Phase set dormant ',mapline%nodfixph
          meqrec%phr(mapline%nodfixph)%phasestatus=PHDORM
! if iadd or irem is equal to mapline%nodfixph change
! direction of the axis
          irem=0; iadd=0
          goto 320
       elseif(mapline%number_of_equilibria.le.5 .and.&
            ((irem.gt.0 .and. irem.eq.mapline%nodfixph) .or. &
            (iadd.gt.0 .and. iadd.eq.mapline%nodfixph))) then
!          write(*,*)'Startnode phase ignored: ',2,mapline%nodfixph,&
!               ceq%tpval(1)
          iadd=0; irem=0
! set the phase dormant and decrease step
!          write(*,559)mapline%nodfixph,axvalok
559       format('Phase set dormant ',i5,1pe14.6)
          meqrec%phr(mapline%nodfixph)%phasestatus=PHDORM
          call map_halfstep(halfstep,1,axvalok,mapline,axarr,ceq)
          if(gx%bmperr.eq.0) then
! jump back without setting halfstep=0
             goto 320
          elseif(nax.gt.1 .and. bytaxis.eq.0) then
! try to change axis with active condition.
             if(meqrec%tpindep(1)) then
                if(ocv()) write(*,*)'Restoring T 7: ',tsave,axvalok
                ceq%tpval(1)=tsave
             endif
! try to change fix phase ...
             write(*,*)'Trying to change fix phase'
             gx%bmperr=0
! if active axis condition is extensive we must change condition value!!
!
             bytaxis=abs(mapline%axandir)
             call locate_condition(axarr(bytaxis)%seqz,pcond,ceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Cannot locate condition: ',axarr(bytaxis)%seqz
                goto 1000
             endif
             svrrec=>pcond%statvar(1)
             call condition_value(1,pcond,zzz,ceq)
             if(gx%bmperr.ne.0) goto 1000
             if(svrrec%argtyp.eq.1 .and. svrrec%statevarid.ge.10) then
! 0 is not good check, it can be a component
! NOTE: If we extract value for currect fix phase we must change axvals/axvals2
!              i1=svr2%argtyp; i2=svr2%phase; i3=svr2%compset
                
                svrtarget=svrrec
                svrtarget%argtyp=3
                svrtarget%phase=mapline%stableph(1)%ixphase
                svrtarget%compset=mapline%stableph(1)%compset
! This extracts the composition of the entered phase for first new line
! we must use a pointer in state_variable_val
                svr2=>svrtarget
                call state_variable_val(svr2,xxx,ceq)
                if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to extract the value, 0 means to set the value
                call condition_value(0,pcond,xxx,ceq)
!                write(*,*)'Old/New axis condition: ',zzz,xxx,pcond%prescribed
!             else
!                write(*,*)'Axis is potential, same value',svrrec%statevarid
             endif
!-------------------------------------------------
             call map_bytfixphase(mapline,nax,mapline%meqrec,xxx,ceq)
             if(gx%bmperr.eq.0) then
                axvalok=zero; goto 320
             endif
!
             write(*,561)gx%bmperr,ceq%tpval(1),axvalok
561          format('Repeated error 9, try to change axis',i5,F8.2,1pe14.6)
             write(*,*)'Trying to change axis with acitive condition'
             gx%bmperr=0
             bytaxis=1
             call map_force_changeaxis(maptop,mapline,mapline%meqrec,&
                  nax,axarr,axvalok,ceq)
             if(gx%bmperr.eq.0) goto 320
             call map_lineend(mapline,axvalok,ceq)
          else
! there is a persistent error, take another line, set error code
             if(gx%bmperr.eq.0) then
                write(*,*)'SMP2A persistent error?'
                gx%bmperr=4399
             endif
             call map_lineend(mapline,axvalok,ceq)
          endif
       endif
       if(mapline%more.eq.0) then
! This is the last equilibrium at axis limit
          if(irem.gt.0) then
! terminate the line and check if there are other lines to calculate
             call map_lineend(mapline,axvalok,ceq)
             goto 300
          elseif(iadd.gt.0) then
             if(ocv()) write(*,*)'New phase at axis limit, IGNORE',iadd
             meqrec%phr(iadd)%dormlink=meqrec%dormlink
             meqrec%dormlink=iadd
             meqrec%phr(iadd)%phasestatus=PHDORM
!             meqrec%phr(iadd)%curd%status2=&
!                  ibset(meqrec%phr(iadd)%curd%status2,CSSUS)
!             meqrec%phr(iadd)%curd%status2=&
!                  ibset(meqrec%phr(iadd)%curd%status2,CSFIXDORM)
             goto 320
          endif
       endif
!       write(*,*)'New set of stable phases: ',iadd,irem,ceq%tpval(1)
! calculate the exact value of the variable axis for the phase change
! then check if we have already found this node point and if not
! generate new start points with and without the phase
! HERE WE CREATE A NODE WITH NEW EXIT LINES
       call map_calcnode(irem,iadd,maptop,mapline,mapline%meqrec,axarr,ceq)
! segmentation fault in map_calcnode 170518 !!
!       write(*,*)'Back from map_calcnode',gx%bmperr,irem,iadd
       if((gx%bmperr.ne.0 .or. irem.ne.0 .or. iadd.ne.0) .and. noderrmess) then
          write(*,777)gx%bmperr,irem,iadd,noderrmess,ceq%tpval(1)
777       format('SMP problem calculating node: ',3i5,l2,F8.2)
          noderrmess=.false.
       endif
       noderror: if(gx%bmperr.ne.0) then
! if error one can try to calculate using a shorter step or other things ...
!          write(*,*)'Error return from map_calcnode: ',gx%bmperr
          if(gx%bmperr.eq.4353) then
! this means node point not global, the line leading to this is set inactive
! and we should not generate any startpoint.             
             write(*,*)'Setting line inactive',mapline%lineid
             mapline%status=ibset(mapline%status,EXCLUDEDLINE)
             call map_lineend(mapline,axvalok,ceq)
             goto 805
          endif
          if(meqrec%tpindep(1)) then
! restore the original temperature, maybe also compositions ...
!             write(*,*)'Restored T 5: ',tsave,axvalok
             ceq%tpval(1)=tsave
          endif
! restore here creates an infinite loop with no axis increment in map2-crmo
!          write(*,*)'Restore constitutions 4'
          call restore_constitutions(ceq,copyofconst)
!          write(*,800)'map_calcnode error: ',gx%bmperr,mapline%problems,&
!               mapline%lasterr,axvalok
800       format(a,3i5,1pe12.4)
          gx%bmperr=0
          call map_halfstep(halfstep,0,axvalok,mapline,axarr,ceq)
!          write(*,*)'back from halfstep 2',halfstep,gx%bmperr
!          if(gx%bmperr.eq.0.and. halfstep.le.5) then
          if(gx%bmperr.eq.0.and. halfstep.le.4) then
             goto 320
          elseif(nax.gt.1 .and. bytaxis.eq.0) then
! try to change axis with active condition.
             if(ocv()) write(*,*)'Trying to change axis with active condition'
             gx%bmperr=0
             if(meqrec%tpindep(1)) then
                if(ocv()) write(*,*)'Restoring T 6: ',tsave,axvalok
                ceq%tpval(1)=tsave
             endif
             if(ocv()) write(*,803)'Repeated error 2, try to change axis',&
                  gx%bmperr,halfstep,ceq%tpval(1)
803          format(a,i5,i3,1pe12.4)
             bytaxis=1; gx%bmperr=0
             call map_force_changeaxis(maptop,mapline,mapline%meqrec,&
                  nax,axarr,axvalok,ceq)
             if(gx%bmperr.eq.0) goto 320
             call map_lineend(mapline,axvalok,ceq)
!          elseif(bypass.eq.0) then
! Problem with 18 component system a phase pops up and down
! Cleanup needed ...
!             bypass=1
!             write(*,*)'SMP2A problem calculate node, try bypass'
!             goto 310
          else
!             write(*,*)' *** Repeated errors calling map_calcnode,',&
!                  ' terminate line',gx%bmperr
! terminate line and follow another line, error reset inside map_lineend
             if(gx%bmperr.eq.0) gx%bmperr=4369
             call map_lineend(mapline,axvalok,ceq)
          endif
       endif noderror
! we come here if a new node has been calculated and stored
       axvalok=zero
    else
! phasechance: els: Here neither iadd or irem>0, we should never be here
! and nno error ... we should go back to label 3000
       write(*,*)'SMPA no phase change?',gx%bmperr,iadd,irem
       stop 'Report this error to the OC development team!'
    endif phasechange
! we have finished a line and look for another at label 300
805 continue
    write(kou,808)mapline%number_of_equilibria,ceq%tpval(1)
808 format('Finishing line with ',i5,' equilibria at T=',0pF8.2,' ')
    mapline%problems=0
    mapline%lasterr=0
    goto 300
!-----------------------------------------------------
! we come here when there are no more lines to calculate
900 continue
!-----------------------------------------------------
! jump here if errors above
1000 continue
!--------------------------------------------------
    return
  end subroutine map_doallines
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine bombmatta
!\begin{verbatim}
  subroutine bombmatta(maptop,nax,axarr,seqxyz,starteqs)
! calculate a number of equilibria inside the region of x and y
!
! nax is the number of axis (can be just one for STEP)
! axarr is an array of records specifying the axis for the step/map
! seqxyz are intial values for number of nodes and lines
! starteq is an equilibrium data record, if there are more start equilibria
! they are linked using the ceq%next index
    implicit none
    integer nax,seqxyz(*)
    type(map_axis), dimension(nax) :: axarr
!    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(starteqlista), dimension(*) :: starteqs
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq,starteq
    type(gtp_condition), pointer :: xcond,ycond
    type(gtp_phase_varres), pointer :: phres
    integer s1,s2,s3,n1,n2,lokcs,nel,globalstatus,iph,potax,touse,newset
    integer, allocatable, dimension(:,:) :: phstable,phused
    double precision xval,yval,xlen,ylen
    integer, parameter :: nss=5
! start in the middle, close to end points at the end
    double precision, dimension(nss), parameter :: axinc=&
         [0.49D0, 0.78D0, 0.22D0, 0.01D0, 0.99D0]
    character name*24
    double precision, dimension(nss*nss) ::  xuse,yuse
!
    starteq=>starteqs(1)%p1
    if(nax.ne.2) then
       write(*,*)'S2A only for map with 2 axis'
       goto 1000
    endif
    nel=noel()
    if(allocated(phstable)) then
       deallocate(phstable)
       deallocate(phused)
    endif
    newset=nooftup()
! there cannot be more than nel phases stable
    allocate(phstable(0:nel,nss*nss+5))
    allocate(phused(0:nel,2*nss))
    phstable=0
    write(*,*)'S2A allocate phstable: ',nel,50,size(phstable),newset
    ceq=>starteq
! supress messages from minimizer
    globalstatus=globaldata%status
    globaldata%status=ibset(globaldata%status,GSSILENT)
    potax=0
! extrahera axis variables and their min and max
! identify any potential axis, statevarid=1=T; 2=P; 3=MU, 4=AC; 5=LNAC
    call locate_condition(axarr(1)%seqz,xcond,ceq)
    if(xcond%statvar(1)%statevarid.le.5) potax=1
    call locate_condition(axarr(2)%seqz,ycond,ceq)
    if(ycond%statvar(1)%statevarid.le.5) potax=2
    if(gx%bmperr.ne.0) goto 1000
    if(potax.gt.0) write(*,*)'S2A potential axis: ',potax
    xlen=axarr(1)%axmax-axarr(1)%axmin
    ylen=axarr(2)%axmax-axarr(2)%axmin
    write(*,*)'S2A axis length: ',xlen,ylen
! start loop
    n1=0
    xloop: do s1=1,nss
! calculate at intervals 0.02 0.1 0.3 0.5 0.7 0.9 0.98 in x and y axis (49 eq)
! set condionon on x axis
       xval=axarr(1)%axmin+axinc(s1)*xlen
! first argument 0 is to set condition, 1 means extract value
       call condition_value(0,xcond,xval,ceq)
       if(gx%bmperr.ne.0) cycle xloop
       yloop: do s2=1,nss
          yval=axarr(2)%axmin+axinc(s2)*ylen
          write(*,'(a,i2,4(1pe12.4))')'S2A x,y: ',n1+1,xval,yval
! set condition on y axis
          call condition_value(0,ycond,yval,ceq)
          if(gx%bmperr.ne.0) cycle yloop
          call calceq2(1,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'S2A failed calculation',gx%bmperr
             gx%bmperr=0; cycle yloop
          endif
          n1=n1+1
          xuse(n1)=xval
          yuse(n1)=yval
! loop to extract stable phases, there can be new composition sets
          n2=0
! start from 2 as first phase_varres is the stable_el_refernce phase
          do lokcs=2,nooftup()
             phres=>ceq%phase_varres(lokcs)
             if(phres%phstate.ge.PHENTSTAB) then
                n2=n2+1
                iph=phres%phlink
                call get_phase_name(iph,1,name)
                if(gx%bmperr.ne.0) gx%bmperr=0
                write(*,'(a,2i2,i5,2x,a)')'S2A stable:',s1,n1,lokcs,trim(name)
! save lokcs as we can have several composition sets
                phstable(n2,n1)=lokcs
             endif
          enddo
! number of stable phases at this equilibrium
          phstable(0,n1)=n2
       enddo yloop
    enddo xloop
! we have calculate all 25 equilibria
    do s1=1,n1
       write(*,'(a,i3,2x,i2,2x,5i5)')'S2A equil: ',s1,(phstable(s2,s1),s2=0,nel)
    enddo
! now decide which points to use as start points, skip points with 
    phused=0
    touse=0
! skip points with phases already used
    all: do s1=1,n1
       if(phstable(0,s1).eq.0) cycle all
       write(*,'(a,5i5)')'S2A compare equil',s1,n1,phstable(0,s1),touse
       phases: do s2=1,phstable(0,s1)
          newset1: do s3=1,touse
! compare with saved equil, skip if an equilibrium has the same phases
             if(phstable(s2,s1).eq.phused(s2,s3)) cycle newset1
          enddo newset1
       enddo phases
! if s3 is less than touse we have an equil with a new set of phases
       write(*,*)'S2A skip as same: ',s3,touse
       if(touse.gt.0 .and. s3.gt.touse) cycle all
! this equilibrium has a new set of phases
       touse=touse+1
       do s3=1,phstable(0,s1)
          phused(s3,touse)=phstable(s3,s1)
       enddo
       write(*,'(a,i3,2x,2F12.5,i2,2x,5i5)')'S2A use: ',s1,xuse(s1),yuse(s1),&
            (phused(s2,s1),s2=0,nel)
    enddo all
    newset=nooftup()-newset
    if(newset.gt.0) write(*,*)'S2A created ',newset,' composition sets'
    write(*,*)'S2A equilibria to use: ',touse
1000 continue
! reset the globaldata%status
    globaldata%status=globalstatus
    return
  end subroutine bombmatta

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_startpoint
!\begin{verbatim}
  subroutine map_startpoint(maptop,nax,axarr,seqxyz,inactive,ceq)
! convert a start equilibrium to a start point replacing all but one axis
! conditions with fix phases.  The start equilibrium must be already
! calculated. ceq is a datastructure with all relevant data for the equilibrium
! A copy of ceq and the corresponing meqrec must be made and linked from maprec
! the axis conditions replaced by fix phases are inactive
! maptop is returned as a first nodepoint(although it is not a node)
! nax is number of axis, axarr records with axis information
! seqxyz is array with indices for numbering nodepoints and lines
! inactive is used for map to replaced axis by fix phase
!       and for step inactive(1) nonzero means create just one linehead
! ceq is equilibrium record
    implicit none
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(map_node), pointer :: maptop
    integer nax,seqxyz(*)
    integer inactive(*)
    type(map_axis), dimension(nax) :: axarr
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: neweq
    TYPE(gtp_condition), pointer :: condition,lastcond
    type(meq_setup), pointer :: meqrec
    TYPE(map_line), pointer :: mapline
    TYPE(map_line), dimension(3) :: tmpline
    TYPE(map_node), pointer :: mapnode,tmpnode
    type(gtp_phasetuple), dimension(3) :: forbidden
    type(map_fixph), allocatable :: mapfix
!    type(map_fixph), pointer :: mapfix
    type(gtp_phasetuple), dimension(:), allocatable :: mapfixph
    integer mode,axactive,iax,jp,ieq,naxvar,seqx,kp,zz,kpos,seqy
    character eqname*24
    double precision value
!
!    write(*,*)"Entering map_startpoint"
    nullify(tmpnode)
! replace all but one axis conditions with fix phases.  In ceq we have
! a calculated equilibrium with all conditions. make sure it works
! (without global minimization).  We will save the meq_setup record!
!    write(*,*)'meq_startpoint: allocating meqrec'
    allocate(meqrec)
    meqrec%status=0
! We must use mode=-1 for map_replaceaxis below has to calculate several equil
! and the phr array must not be deallocated.  mapfix will be used later to
! indicate fix and stable phases for different lines (maybe ...)
    mode=-1
    if(allocated(mapfix)) deallocate(mapfix)
!    nullify(mapfix)
!    write(*,*)'SMP2A meq_startpoint: after allocating meqrec 1'
    call calceq7(mode,meqrec,mapfix,ceq)
    if(gx%bmperr.ne.0) then
! try using grid minimizer
       gx%bmperr=0
! most data inside meqrec like meqrec%phr are deallocated inside calceq7
! but calling it with mode=-1 it is kept so it must be deallocated here 
! BUG here 2019.03.03 not allocated!
       if(allocated(meqrec%phr)) deallocate(meqrec%phr)
       call calceq7(1,meqrec,mapfix,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error calling calceq7 in map_startpoint A',gx%bmperr
          goto 1000
       endif
       call calceq7(mode,meqrec,mapfix,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error calling calceq7 in map_startpoint B',gx%bmperr
          goto 1000
       endif
    endif
! check if equilibrium inside axis limits ...
    do iax=1,nax
       call locate_condition(axarr(iax)%seqz,condition,ceq)
       if(gx%bmperr.ne.0) goto 1000
       call condition_value(1,condition,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
       if(value.lt.axarr(iax)%axmin .or. value.gt.axarr(iax)%axmax) then
          write(*,*)'Startpoint outside axis limits',iax,value
          gx%bmperr=4225; goto 1000
       endif
    enddo
!    write(*,1001)'After calceq7: ',(meqrec%phr(jp)%curd%amfu,&
!         jp=1,meqrec%nphase)
1001 format(a,6(1pe12.4))
200 continue
!---------------------------------- moved before creating first linehead
! create map_node normally with two exiting lines but in some cases more.
    if(associated(maptop)) then
! we have already a maptop record, add a new mapnode at the circular list end
! set appropriate next/previous/first links
       tmpnode=>maptop%previous
       allocate(maptop%previous)
! initiate all status bits to zero
       maptop%previous%status=0
       tmpnode%next=>maptop%previous
       mapnode=>maptop%previous
! initiate mapnode
       mapnode%noofstph=-1
       mapnode%previous=>tmpnode
       mapnode%next=>maptop
       mapnode%first=>maptop
       mapnode%seqx=tmpnode%seqx+1
       mapnode%nodefix%ixphase=0
       mapnode%status=0
       mapnode%artxe=0
       mapnode%globalcheckinterval=mapglobalcheck
!       write(*,*)'creating another mapnode record',mapnode%seqx
! nullify here indicates more than one node record
       nullify(tmpnode)
    else
! This is the first (and maybe only) mapnode record (later maptop)
!       write(*,*)'Creating first maptop'
! UNFINISHED: VALGRIND indicates loss of >24000 bytes in map_startpoint 
       allocate(maptop)
       mapnode=>maptop
! inititate status and links
       mapnode%status=0
       mapnode%noofstph=meqrec%nstph
       mapnode%savednodeceq=-1
!       mapnode%noofstph=-1
       mapnode%next=>mapnode
       mapnode%previous=>mapnode
       mapnode%first=>mapnode
       mapnode%number_ofaxis=nax
       mapnode%nodefix%ixphase=0
       mapnode%status=0
       mapnode%artxe=0
! type_of_node =1 step special; =2 step scheil; =3 step tzero;
!              =4 step paraequil; =5 step nple
! same indices used in stepspecial in pmon
       mapnode%type_of_node=0
       mapnode%globalcheckinterval=mapglobalcheck
! if there is a previous MAP/STEP then 
! seqx and seqy pass on the last used indices for _MAPNODE and _MAPLINE
!       write(*,*)'Seqxyz 1: ',seqxyz(1),seqxyz(2)
! seqx is set to 0 here, will be increemented by 1 at copy_equilibrium
       mapnode%seqx=seqxyz(1)
       mapnode%seqy=seqxyz(2)
       if(ocv()) write(*,*)'created maptop',maptop%seqx
! set the tieline_inplane or not
! For step calculation, tieline_inplane=0
! if there are more than one condition on an extensive_variable
! that is not an axis variable then no tielines in plane, tieline_inplane=-1
! If there are tie_lines in plane then tieline_inplane=1
       mapnode%tieline_inplane=tieline_inplane(nax,axarr,ceq)
       if(mapnode%tieline_inplane.lt.0) then
          write(*,*)'Mapping without tie-lines in the plane'
       endif
       tmpnode=>maptop
! forgetting to do this created a crash when plotting ...
       nullify(maptop%plotlink)
    endif
!
!-----------------------------------------------------------------
! if naxvar>1 find a phase to set fix to replace an axis variable
    naxvar=nax
    if(naxvar.gt.1) then
! in tmpline info on fix/stable phases to be stored in linehead records
       call map_replaceaxis(meqrec,axactive,ieq,nax,axarr,tmpline,inactive,&
            forbidden,ceq)
       if(ocv()) write(*,205)'Back from replaceaxis with: ',gx%bmperr,&
            axactive,ieq,&
            tmpline(1)%linefixph(1)%ixphase,tmpline(1)%linefixph(1)%compset,&
            tmpline(1)%stableph(1)%ixphase,tmpline(1)%stableph(1)%compset
205    format(a,3i5,5x,2(2i3))
       if(gx%bmperr.ne.0) goto 1000
       if(ieq.gt.2) then
          write(*,*)'Ignoring 3rd exit from invariant!'
          ieq=2
       endif
    else
! only one axis, i.e. a step command, create a map_node record with 2 lines
       axactive=1
       if(inactive(1).eq.0) then
          ieq=2
       else
          ieq=1
       endif
    endif
!    write(*,1001)'After replace: ',(meqrec%phr(jp)%curd%amfu,&
!         jp=1,meqrec%nphase)
!-----------------------------------------------------------------------
! finished converting a start equilibrium to a start point, 
!    mapnode%type_of_node=0
    mapnode%lines=ieq
! debug listing of links for maptop ...
!    write(*,*)'maptop: ',maptop%noofstph
!    write(*,*)'maptop next: ',maptop%next%noofstph
!    write(*,*)'maptop prev: ',maptop%previous%noofstph
!    write(*,*)'maptop next next: ',maptop%next%next%noofstph
!    write(*,*)'maptop prev prev: ',maptop%previous%previous%noofstph
!    write(*,*)'maptop next prev: ',maptop%next%previous%noofstph
!    if(associated(maptop,maptop%next)) then
!       write(*,*)'maptop and maptop%next is same record'
!    endif
!
! Save the T, P and chemical potentials
    allocate(mapnode%chempots(meqrec%nrel))
    do jp=1,meqrec%nrel
       mapnode%chempots(jp)=ceq%cmuval(jp)
    enddo
    mapnode%tpval=ceq%tpval
    mapnode%nodeceq=>ceq
!-----------------------------------------------------------------------
    if(ocv()) write(*,*)'allocating lineheads: ',ieq,maptop%seqy
! ensure mapnode%lines is correctly set
    allocate(mapnode%linehead(ieq))
    mapnode%lines=ieq
!    mapnode%type_of_node=0
! meqrec%status
    do jp=1,ieq
       mapnode%linehead(ieq)%meqrec%status=0
    enddo
! we can have 3 or more exits if starting inside a 3 phase triagle for isotherm
    if(ieq.lt.3) then
! STEP command: set one exit in each direction of the active axis axactive
! or we found a phase to set fix in a map command?
!       do jp=1,2
       do jp=1,ieq
!--------------------- code moved from map_findline
! make a copy of the equilibrium record
          if(ocv()) write(*,*)'We found a line from node: ',mapnode%seqx
          eqname='_MAPLINE_'
! kpos=10 means write number from position 10
          kpos=10
          seqy=maptop%seqy+1
          call wriint(eqname,kpos,seqy)
!          write(*,*)'Calling copy_equilibrium'
          call copy_equilibrium(neweq,eqname,mapnode%nodeceq)
!          write(*,*)'back from copy_equilibrium 6'
          if(gx%bmperr.ne.0) then
             write(*,*)'Error creating equilibrium: ',eqname
             goto 1000
          endif
          maptop%seqy=seqy
!------------------------------ end code copied
! one line has +axactive, the other -axactive
          if(ieq.eq.2) then
             mapnode%linehead(jp)%axandir=(3-2*jp)*axactive
          else
! this is used for Scheil-Gulliver step with just one axis
!            write(*,*)'SMP2A Scheil map_startpoint: ',inactive(1),jp,ieq
             mapnode%linehead(jp)%axandir=inactive(1)
          endif
          mapnode%linehead(jp)%number_of_equilibria=0
          mapnode%linehead(jp)%first=0
          mapnode%linehead(jp)%last=0
!          mapnode%linehead(jp)%axchange=0
! careful balance between map4 (U-O) and Fe-Mo (map5) macros
          mapnode%linehead(jp)%axchange=-1
!          mapnode%linehead(jp)%axchange=-2
! lineid is set when calculations along the line starts
!          mapnode%linehead(jp)%lineid=0
          mapnode%linehead(jp)%done=0
          mapnode%linehead(jp)%status=0
          mapnode%linehead(jp)%more=1
          mapnode%linehead(jp)%termerr=0
          mapnode%linehead(jp)%firstinc=zero
! saving equilibrium pointer in lineceq
          mapnode%linehead(jp)%lineceq=>neweq
!          mapnode%linehead(jp)%evenvalue=zero
! to ensure small initial steps
          mapnode%linehead(jp)%evenvalue=value+(3-2*jp)*axarr(1)%axinc
!          write(*,*)'evenvalue: ',mapnode%linehead(jp)%evenvalue,value
          mapnode%linehead(jp)%start=>mapnode
          mapnode%linehead(jp)%axfact=1.0D-2
! this is set to zero indicating the stable phases are saved in ceq record
          mapnode%linehead(jp)%nstabph=0
          mapnode%linehead(jp)%lineid=seqy
!          write(*,*)'mapline%lineid assigned',seqy
          mapnode%linehead(jp)%nodfixph=0
! %more is 1 while line is calculated, 0 means terminated at axis limit
! > 0 means error code <0 means exit removed ?? or is it %done ??
          mapnode%linehead(jp)%more=1
!-------------------------
          if(maptop%tieline_inplane.lt.0) then
! tie-lines not in plane, code just copied with some mods from tielines in plane
             kp=tmpline(1)%nfixphases
             mapnode%linehead(jp)%nfixphases=kp
             allocate(mapnode%linehead(jp)%linefixph(kp))
             allocate(mapnode%linehead(jp)%linefix_phr(kp))
!             write(*,454)jp,axactive,mapnode%linehead(jp)%axandir,kp
454          format('Axis direction etc: ',i2,2i4,2x,i3)
             do zz=1,kp
                mapnode%linehead(jp)%linefixph(zz)=tmpline(1)%linefixph(zz)
                mapnode%linehead(jp)%linefix_phr(zz)=tmpline(1)%linefix_phr(zz)
             enddo
! we can have many stable phases
             mapnode%linehead(jp)%nstabph=tmpline(1)%nstabph
             allocate(mapnode%linehead(jp)%stableph(tmpline(1)%nstabph))
             allocate(mapnode%linehead(jp)%stablepham(tmpline(1)%nstabph))
             allocate(mapnode%linehead(jp)%stable_phr(tmpline(1)%nstabph))
             do kp=1,mapnode%linehead(jp)%nstabph
                mapnode%linehead(jp)%stableph(kp)=tmpline(1)%stableph(kp)
                mapnode%linehead(jp)%stablepham(kp)=tmpline(1)%stablepham(kp)
                mapnode%linehead(jp)%stable_phr(kp)=tmpline(1)%stable_phr(kp)
             enddo
!             write(*,*)'allocated size of stableph 1: ',jp,&
!                  size(mapnode%linehead(jp)%stableph)
             if(ocv())write(*,27)'We have a startpoint for no tie-lines map:',&
!                  axactive,mapnode%linehead(jp)%linefixph(1)%phaseix,&
                  axactive,mapnode%linehead(jp)%linefixph(1)%ixphase,&
                  mapnode%linehead(jp)%linefixph(1)%compset,&
                  mapnode%linehead(jp)%nstabph,&
                  (mapnode%linehead(jp)%stableph(kp)%ixphase,&
                  mapnode%linehead(jp)%stableph(kp)%compset,&
                  kp=1,mapnode%linehead(jp)%nstabph)
27           format(a,i3,5x,2i3,5x,i3,2x,10(i5,i2))
!------------------------- below for tielines in plane
          elseif(maptop%tieline_inplane.gt.0) then
! if there are 2 axis there is one fix phase, if 3 axis there are two
! This is not really necessary here but for other nodes with branches it is
             kp=tmpline(1)%nfixphases
!             write(*,*)'tip: Number of fixed phases: ',jp,kp
             mapnode%linehead(jp)%nfixphases=kp
             allocate(mapnode%linehead(jp)%linefixph(kp))
             allocate(mapnode%linehead(jp)%linefix_phr(kp))
             do zz=1,kp
                mapnode%linehead(jp)%linefixph(zz)=tmpline(1)%linefixph(zz)
                mapnode%linehead(jp)%linefix_phr(zz)=tmpline(1)%linefix_phr(zz)
             enddo
! there is just one stable phase
             allocate(mapnode%linehead(jp)%stableph(1))
             allocate(mapnode%linehead(jp)%stable_phr(1))
             mapnode%linehead(jp)%nstabph=1
             mapnode%linehead(jp)%stableph=tmpline(1)%stableph
             mapnode%linehead(jp)%stable_phr=tmpline(1)%stable_phr
! WOW I forgot to allocate stablepham
             if(allocated(tmpline(1)%stableph)) then
                kp=size(tmpline(1)%stableph)
                allocate(mapnode%linehead(jp)%stablepham(kp))
             else
                write(*,*)'SMP: no stablepham array allocated'
                stop
             endif
             if(ocv()) write(*,25)'We have saved a startpoint for map:',&
!                  axactive,mapnode%linehead(jp)%linefixph(1)%phaseix,&
                  axactive,mapnode%linehead(jp)%linefixph(1)%ixphase,&
                  mapnode%linehead(jp)%linefixph(1)%compset,&
                  mapnode%linehead(jp)%nstabph,&
                  mapnode%linehead(jp)%stableph(1)%ixphase,&
                  mapnode%linehead(jp)%stableph(1)%compset
25           format(a,i3,5x,2i3,5x,i3,2x,2i3)
!------------------------- below for STEP
          else
! this is for STEP
             if(ocv()) write(*,*)'For STEP no need of fixed phases.'
!             write(*,*)'SMP2A Scheil step here?'
             mapnode%linehead(jp)%nfixphases=0
             allocate(mapnode%linehead(jp)%stableph(meqrec%nstph))
             allocate(mapnode%linehead(jp)%stable_phr(meqrec%nstph))
! UNFINISHED check why no allocation of stablepham ??
             mapnode%linehead(jp)%nstabph=meqrec%nstph
             do kp=1,mapnode%linehead(jp)%nstabph
                zz=meqrec%stphl(kp)
                mapnode%linehead(jp)%stableph(kp)%ixphase=meqrec%phr(zz)%iph
                mapnode%linehead(jp)%stableph(kp)%compset=meqrec%phr(zz)%ics
                mapnode%linehead(jp)%stable_phr(kp)=zz
             enddo
          endif
!-------------------------
          nullify(mapnode%linehead(jp)%end)
          mapnode%linehead(jp)%nodfixph=0
       enddo
    else
! when more than two exits the set of stable phases must be different for
! each line.  This can happen if we start in a three-phase region in an
! isothermal section with tie-lines in plane
       write(*,*)'Cannot handle more than two exits from start equilibrium'
       gx%bmperr=4226; goto 1000
    endif
! mapnode must have pointers to its own copies of ceq and meqrec
    eqname='_MAPNODE_'
    jp=10
! maptop%next is the most recent created mapnode ??
    seqx=maptop%next%seqx+1
!    write(*,*)'SMP2A New mapnode index: ',seqx,&
!         maptop%next%seqx,maptop%previous%seqx
    seqx=max(maptop%next%seqx,maptop%previous%seqx)+1
    maptop%next%seqx=seqx
!    write(*,666)seqx,maptop%seqx,maptop%next%seqx,maptop%previous%seqx
666 format('maptop seqx: ',10i3)
    call wriint(eqname,jp,seqx)
! make a copy of ceq in a new equilibrium record with the pointer neweq
! This copy is a record in the array "eqlista" of equilibrium record, thus
! it will be updated if new composition sets are created in other threads.
    call copy_equilibrium(neweq,eqname,ceq)
!    write(*,*)'Created MAPNODE ',seqx
    if(gx%bmperr.ne.0) then
       write(*,*)'Error in startpoint creating equilibrium: ',eqname
       goto 1000
    endif
    if(associated(mapnode,maptop)) maptop%seqx=seqx
    mapnode%nodeceq=>neweq
! If the new node has two stoichiometric phases then mapline%status
! Copy the current meqrec to mapnode, the mapline records
! will generate their own new meqrec records when they are activated
! if the phr array is allocated then deallocate it as it is no longer needed
    if(allocated(meqrec%phr)) then
       deallocate(meqrec%phr)
    endif
    mapnode%meqrec=meqrec
! trying to reduce memory loss
    deallocate(meqrec)
!    write(*,*)'We are here 15!'
! NOTE: The phr array has been deallocated, maybe it should be kept ...
! but then we must change mode to -1 in the call to calceq7 above
!---------------------
! The lines below must be done when creating the mapnod%linehead record
! we must have separate copies of meqrec and ceq for use in each thread
!    mapline%meqrec=mapnode%meqrec
!    mapline%ceq=mapnode%ceq
! finished what must be done when creating mapnode%linehead
!
    if(ocv()) write(*,*)'Exiting map_startpoint',gx%bmperr
1000 continue
    if(gx%bmperr.ne.0 .and. associated(tmpnode)) then
! we have created a maptop record but then had an error, nullify mapnode
       write(*,*)'Nullifying maptop: ',gx%bmperr
       nullify(maptop)
    endif
    return
  end subroutine map_startpoint

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_replaceaxis
!\begin{verbatim}
  subroutine map_replaceaxis(meqrec,axactive,ieq,nax,axarr,tmpline,&
       inactive,forbidden,ceq)
! replace an axis condition with a fix phase
! meqrec is equilibrium calculation record
! axactive is the axis with active condition, ieq is number of exiting lines
! ieq is the number of lines exiting from the startpoint
! nax is number of axis, axarr are description of the axis
! axarr is array with axis records
! tmpline is to transfer some line data to calling routine
! inactive is not really used.
! forbidden are phasetupes with forbidden phases
! ceq is equilibrium record
    implicit none
    type(meq_setup), pointer :: meqrec
    integer nax,axactive,ieq
    type(map_axis), dimension(nax) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
    type(map_line), dimension(3) :: tmpline
    integer inactive(*)
    type(gtp_phasetuple), dimension(*) :: forbidden
!\end{verbatim}
    integer iph,jph,naxvar,iax,tip,jj,jax,irem,iadd,kj,nrel,sj,kax
    integer ics,lokph,lokcs,kph,kcs,forbiddenix,sph,mapx
    double precision aval,avalm,xxx,yyy,savamfu(3)
! dummy phase tuple, maybe use nullify instead?
    type(gtp_phasetuple) :: zerotup
    type(gtp_condition), pointer :: pcond
    integer, dimension(:), allocatable :: axis_withnocond
! handle change of condition value
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(gtp_state_variable), target :: svrtarget
! turns off converge control for T
    integer, parameter :: inmap=1
!
    zerotup%lokph=0
    zerotup%compset=0
    zerotup%ixphase=0
    zerotup%lokvares=0
    zerotup%nextcs=0
!
    nrel=noel()
    tip=tieline_inplane(nax,axarr,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'In map_replaceaxis ',tip
!-----------------------------------------------------------------
! check if start point is an invariant equilibria, can easily happen in 
! ternary isotherms
    if(inveq(jj,ceq)) then
       if(tip.gt.0) then
! ignore this for less than 3 components
          if(nrel.eq.3) then
! we are in an isothermal triangle, 3 startlines
             write(*,*)'Start equilibrium is invariant',jj
             ieq=3
!             goto 1000
          elseif(nrel.gt.3) then
! I do not know what kind of equilibrium this is
             write(*,160)
160          format('Start equilibrium invariant with tie-lines in plane',&
                  ' but not 3 components'/'I do not know how to handle this')
             gx%bmperr=4399
             goto 1000
          endif
       else
! start equilibrium for a system without tie-lines in plane is invariant
! a rare case
          write(*,161)
161       format('Start equilibrium invariant without tie-lines in plane'/&
               'I do not know how to handle this')
          gx%bmperr=4399
          goto 1000
       endif
    endif
    naxvar=nax
! zero the number of fix phases and allocate data for the lines needed
    tmpline%nfixphases=0
    allocate(tmpline(1)%linefixph(naxvar-1))
    allocate(tmpline(1)%linefix_phr(naxvar-1))
!========================================================== tie-lines in plane
    tieline_in_plane: if(tip.eq.1) then
! We have tie-lines in the plane, only one stable phase in addition to fix
!       write(*,*)'map_replaceaxis: allocate: tmpline(1)%stableph(1)'
       allocate(tmpline(1)%stableph(1))
       allocate(tmpline(1)%stablepham(1))
       allocate(tmpline(1)%stable_phr(1))
       allocate(axis_withnocond(nax))
       axis_withnocond=0
       stablephases: if(meqrec%nstph.gt.1) then
! and we have two or more stable phase, we can directly generate startpoints
100       continue
          if(meqrec%nstph.eq.3) then
! this is a unique case when we must create 3 lines
!             write(*,*)'Startpoint inside invariant not yet implemented'
!             gx%bmperr=4399; goto 1000
!
! save some data ...??
!             do jph=1,3
!                jj=meqrec%stphl(jph)
!                savamfu(iph)=meqrec%phr(jj)%curd%amfu
!                stableph(1,jph)=meqrec%phr(jj)%iph
!                stableph(2,jph)=meqrec%phr(jj)%ics
!             enddo
! loop for the 3 stable phases setting one of them as fix in turn
             meqrec%nfixph=1
             meqrec%nstph=1
             forbiddenix=3
             fixphaseloop: do jph=1,3
! all phases are already set as stable
                kj=meqrec%stphl(jph)
!                write(*,*)'tmpline 1: ',jph,kj
                if(jph.gt.1) then
                   allocate(tmpline(jph)%linefixph(1))
                   allocate(tmpline(jph)%linefix_phr(1))
                   allocate(tmpline(jph)%stableph(1))
                   allocate(tmpline(jph)%stablepham(1))
                   allocate(tmpline(jph)%stable_phr(1))
                endif
! do we need to set values in meqrec??
!                meqrec%fixph(1,1)=meqrec%phr(kj)%iph
!                meqrec%fixph(2,1)=meqrec%phr(kj)%ics
!                meqrec%fixpham(1)=zero
                tmpline(jph)%nfixphases=1
                tmpline(jph)%linefixph(1)=zerotup
!                write(*,*)'tmpline 2A: ',jph,kj
!                write(*,*)'tmpline 2C: ',allocated(meqrec%phr)
!                write(*,*)'tmpline 2B: ',meqrec%phr(kj)%iph
!                write(*,*)'tmpline 2C: ',allocated(tmpline(jph)%linefixph)
                tmpline(jph)%linefix_phr(1)=kj
                tmpline(jph)%linefixph(1)%ixphase=meqrec%phr(kj)%iph
!                write(*,*)'SMP2A tmpline 3: ',jph,kj
                tmpline(jph)%linefixph(1)%compset=meqrec%phr(kj)%ics
!                write(*,*)'SMP2A tmpline 4: ',jph,kj
                tmpline(jph)%nstabph=1
                kph=jph+1
                if(kph.gt.3) kph=1
                sph=meqrec%stphl(kph)
                tmpline(jph)%axandir=1
                write(*,*)'tmpline',kph,sph,tmpline(jph)%axandir
                tmpline(jph)%stableph(1)=zerotup
                tmpline(jph)%stableph(1)%ixphase=meqrec%phr(sph)%iph
                tmpline(jph)%stableph(1)%compset=meqrec%phr(sph)%ics
                tmpline(jph)%stablepham(1)=one
                tmpline(jph)%stable_phr(1)=sph
! lines:  (fix,stable,forbidden) :  (1,2,3);   (2,3,1);   (3,1,2)
! we must mark the third phase as forbidden !!!
                jj=meqrec%stphl(forbiddenix)
                write(*,*)'tmpline 5',forbiddenix,jj
                forbidden(jph)=zerotup
                forbidden(jph)%ixphase=meqrec%phr(jj)%iph
                forbidden(jph)%compset=meqrec%phr(jj)%ics
                forbiddenix=forbiddenix+1
                if(forbiddenix.gt.3) forbiddenix=1
             enddo fixphaseloop
             write(*,65)'Lines: ',&
                  tmpline(1)%linefixph(1)%ixphase,&
                  tmpline(1)%stableph(1)%ixphase,forbidden(1)%ixphase,&
                  tmpline(2)%linefixph(1)%ixphase,&
                  tmpline(2)%stableph(1)%ixphase,forbidden(2)%ixphase,&
                  tmpline(3)%linefixph(1)%ixphase,&
                  tmpline(3)%stableph(1)%ixphase,forbidden(3)%ixphase
65           format(a,3i4,5x,3i4,5x,3i4)
! we should set the axis composition to the stable phase ...
! and we should test ...
             goto 1000
! this is end of generating startpoint from a ternary isothermal triangle
          endif
          write(*,*)'Tie-lines in the plane and start equilibrium with',&
               ' several stable phases'
          jax=0
!          call list_conditions(kou,ceq)
          do iax=1,nax
             call locate_condition(axarr(iax)%seqz,pcond,ceq)
! skip axis already removed
             if(pcond%active.eq.1) cycle
             if(pcond%statev.ge.10) then
! Best to replace an extensive variable with a fix phase
! But we cannot use a condition N=1 or B=1 for example.  It must depend on
! a component
!                write(*,*)'Condition 1: ',iax,pcond%seqz
!                if(pcond%indices(1,1).eq.0) cycle
                if(pcond%statvar(1)%argtyp.eq.0) cycle
!                write(*,*)'Condition 2: ',iax,pcond%indices(1,1)
                jax=iax; exit
             endif
          enddo
          if(jax.eq.0) then
! we must accept to replace a potential axis, use one depending on a component
! If we have a P-T diagram? This would not work
             do iax=1,nax
                call locate_condition(axarr(iax)%seqz,pcond,ceq)
                if(pcond%statev.gt.2) then
                   jax=iax; exit
                endif
             enddo
          endif
!----------------------------------------------
! determine a phase to set fix with zero amount
          avalm=1.0D5
          if(ocv()) write(*,*)'Removing axis ',jax,&
               ', looking for the phase to fix'
! select the phase with smallest amount ... phr has been deallocated ...
!                aval=ceq%phase_varres(lokcs)%amfu
          do jph=1,meqrec%nstph
             jj=meqrec%stphl(jph)
! amfu is amount formula units, abnorm(1) is atoms/formula units
             aval=meqrec%phr(jj)%curd%amfu*meqrec%phr(jj)%curd%abnorm(1)
             if(aval.lt.avalm) then
                kph=meqrec%phr(jj)%iph
                kcs=meqrec%phr(jj)%ics
                kj=jj
                avalm=aval
! we have 2 stable phases, jph is 1 or 2
                sj=3-jph
             endif
          enddo
! The phase meqrec%phr(kj)%iph/ics should be set fix 
          sj=meqrec%stphl(sj)
!          write(*,73)'Fix phase: ',kj,meqrec%phr(kj)%iph,meqrec%phr(kj)%ics,&
!               ' Stable phase: ',sj,meqrec%phr(sj)%iph,meqrec%phr(sj)%ics
73        format(a,3i4,a,3i4)
          meqrec%phr(kj)%curd%dgm=zero
          meqrec%phr(kj)%curd%amfu=zero
          meqrec%phr(kj)%stable=1
          meqrec%phr(kj)%phasestatus=PHFIXED
! The array fixph contains also phases with explicit condition to be fixed
          meqrec%nfixph=meqrec%nfixph+1
          meqrec%fixph(1,meqrec%nfixph)=kph
          meqrec%fixph(2,meqrec%nfixph)=kcs
          meqrec%fixpham(meqrec%nfixph)=zero
! and the axis condition pcond should be removed
          pcond%active=1
          inactive(1)=inactive(1)+1
          inactive(inactive(1))=pcond%seqz
!          meqrec%inactiveaxis(1)=pcond%seqz
!          write(*,77)jax,pcond%seqz,pcond%prescribed
77        format(' Removing condition: ',2i3,2(1pe12.4))
! We have tried not to replace T or P,
! but if this is done it must be indicated specially like this
          if(pcond%statev.eq.1) then
             meqrec%tpindep(1)=.TRUE.
             if(ocv()) write(*,*)'Marking that T is variable'
          elseif(pcond%statev.eq.2) then
             meqrec%tpindep(2)=.TRUE.
          endif
! set amount of stable phase
          meqrec%phr(sj)%curd%amfu=one
! if both axis are extensive (isothermal section) modify active axis condition
! to be the composition of the stable phase
          kax=3-jax
          call locate_condition(axarr(kax)%seqz,pcond,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Cannot locate condition: ',axarr(kax)%seqz
             goto 1000
          endif
! first argument 1 means extract value of condition
          call condition_value(1,pcond,xxx,ceq)
          if(pcond%statev.ge.10) then
!             write(*,*)'isothermal section'
             svrrec=>pcond%statvar(1)
! NOTE: If we change fix/entered phase we must change axvals/axvals2
             svrtarget=svrrec
             svrtarget%argtyp=3
             svrtarget%phase=meqrec%phr(sj)%iph
             svrtarget%compset=meqrec%phr(sj)%ics
! This extracts the composition of the entered phase for first new line
! we must use a pointer in state_variable_val
             svr2=>svrtarget
             call state_variable_val(svr2,yyy,ceq)
             if(gx%bmperr.ne.0) goto 1000
!             write(*,71)jax,xxx,yyy
71           format('Change ',i3,' axis condition from/to ',2F10.6)
! first argument 1 means to extract the value, 0 means to set the value
             call condition_value(0,pcond,yyy,ceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Cannot set axis condition'
                gx%bmperr=4399; goto 1000
             endif
          endif
!---------------------------------------------------
! calculate the equilibrium with the new set of conditions
          if(ocv()) write(*,*)'Calling meq_sameset inside  map_replaceaxis'
          irem=0; iadd=0;
!          write(*,*)'smp2A calling meq_sameset from map_replaceaxis'
          call meq_sameset(irem,iadd,mapx,meqrec,meqrec%phr,inmap,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error calling meq_sameset in startpoint: ',gx%bmperr
             goto 1000
          elseif(irem.gt.0 .or. iadd.gt.0) then
             write(*,*)'Change of phase set in startpoint...',irem,iadd
             gx%bmperr=4227; goto 1000
          endif
!------------------------------------------------------
          if(ocv()) write(*,*)'A successful calculation with one axis',&
               ' condition replaced by a fix phase.'
          if(ocv()) write(*,*)'Released axis: ',jax,' fix phase: ',kph,kcs
          axis_withnocond(jax)=1
          naxvar=naxvar-1
          if(naxvar.eq.1) then
! when we are here we have a start point and can determine the number of exits
! for the moment just assume 2nd axis is the remaining condition!!
             tmpline(1)%nfixphases=1
             tmpline(1)%linefixph=zerotup
! kj and kph set in loop above ... hope they have not changed
             tmpline(1)%linefixph%ixphase=kph
             tmpline(1)%linefixph%compset=kcs
             tmpline(1)%linefix_phr=kj
             tmpline(1)%nstabph=0
! Note meqrec%phr is a TYPE meq_phase with an link curd to phase_varres
! meqrec%phr is a more complex TYPE
             do jph=1,meqrec%nstph
                jj=meqrec%stphl(jph)
                if(meqrec%phr(jj)%iph.eq.kph .and.&
                     meqrec%phr(jj)%ics.eq.kcs) cycle
                tmpline(1)%stableph(1)=zerotup
                tmpline(1)%stableph(1)%ixphase=meqrec%phr(jj)%iph
                tmpline(1)%stableph(1)%compset=meqrec%phr(jj)%ics
                tmpline(1)%stable_phr(1)=jj
                tmpline(1)%nstabph=tmpline(1)%nstabph+1
!                tmpline(1)%nstabph=1
! why exit?? Maybe because there can only be a single phase!!
!                exit
             enddo
             if(tmpline(1)%nstabph.eq.0) then
                write(*,*)'No stable phase !!'
                stop
             endif
! This is the axis with active condition
             axactive=2
             ieq=2
          else
             write(*,*)'Not implemented more than 2 axis'
             gx%bmperr=4228; goto 1000
          endif
! ========================================== tie-lines in plane and one phase
       else ! we have just a single phase stable we must move in some direction
! ceq%multiuse is direction
!          write(*,*)'SMP2A Tie-line in plane and single phase,',&
!               ' This may not work ... '
          call map_startline(meqrec,axactive,ieq,nax,axarr,tmpline,ceq)
          if(gx%bmperr.ne.0) goto 1000
       endif stablephases
! ============================================= no tie-lines in plane
    else !tie-lines NOT in the plane
! I am not sure what stableph and axis_withnocond are used for ...
!   write(*,*)'SMP2A multiple startpoint without tie-lines in plane not allowed'
!       gx%bmperr=4399; goto 1000
       allocate(axis_withnocond(nax))
       axis_withnocond=0
       call map_startline(meqrec,axactive,ieq,nax,axarr,tmpline,ceq)
       if(gx%bmperr.ne.0) goto 1000
    endif tieline_in_plane
!       
! check if more axis must be released
900 continue
    if(nax.gt.2) then
       write(*,*)'Cannot handle more than 2 axis at present'
       gx%bmperr=4228
    endif
1000 continue
!    write(*,*)'Return from map_replaceaxis with conditions: '
!    call list_conditions(kou,ceq)
    return
  end subroutine map_replaceaxis

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_startline
!\begin{verbatim}
  subroutine map_startline(meqrec,axactive,ieq,nax,axarr,tmpline,ceq)
! find a phase to fix to replace an axis condition when we 
! do not have tie-lines in the plane or when we 
! have tie-lines in the plane but start in a single phase region
! meqrec is equilibrium record already initiated
! axactive is set to the axis with active condition
! ieq is the number of lines exiting from the startpoint
! nax is number of axis, axarr are description of the axis
! axarr are axis records
! tmpline is a line record ... not needed ... ??
    implicit none
    integer nax,axactive,ieq
    type(meq_setup), pointer :: meqrec
    type(map_line), dimension(2) :: tmpline
    type(map_axis), dimension(nax) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer jax,iax,idir,irem,iadd,iph,jj,jph,kph,ll,mapx
    integer :: maxtry=0
    integer, parameter :: nstabphdim=20
    double precision curval,startval
    type(gtp_condition), pointer :: pcond
! turns off converge control for T
    integer, parameter :: inmap=1
    save maxtry
!
!    write(*,*)'In map_startline, find a phase to set fix',ceq%multiuse
! start in negative direction unless direction given
    idir=-1
    if(ceq%multiuse.ne.0) then
       if(abs(ceq%multiuse).gt.nax) then
          write(*,*)'Error in direction, no such axis, ',ceq%eqno,ceq%multiuse
! this can happen for startpoints.  21 is lower left, 22 is lower right
! 23 is upper left, 24 is upper right and 30 is in the middle          
! Try to generate several directions for each, at present just one
          if(ceq%multiuse.eq.21) then
! directions +1 and +2
             call list_conditions(kou,ceq)
             jax=1; idir=1
          elseif(ceq%multiuse.eq.22) then
! directions +1 and -2
             call list_conditions(kou,ceq)
             jax=2; idir=-1
          elseif(ceq%multiuse.eq.23) then
! directions -1 and +2
             call list_conditions(kou,ceq)
             jax=1; idir=-1
          elseif(ceq%multiuse.eq.24) then
! directions -1 and -2
             call list_conditions(kou,ceq)
             jax=2; idir=-1
          elseif(ceq%multiuse.eq.30) then
! all 4 directions ...
             call list_conditions(kou,ceq)
          else
             write(*,*)'Error in direction, no such axis, ',ceq%multiuse
             gx%bmperr=4229; goto 1000
          endif
       else
! direction is +/-axis
!          write(*,*)'SMP2A direction: ',ceq%multiuse
          if(ceq%multiuse.gt.0) idir=1
          jax=abs(ceq%multiuse)
          call locate_condition(axarr(jax)%seqz,pcond,ceq)
!          write(*,*)'SMP2A axis condition: ',pcond%statev,gx%bmperr
          if(gx%bmperr.ne.0) goto 1000
       endif
    else
! no axis selected
!       write(*,*)'SMP2A no direction',ceq%multiuse,nax
       jax=0
       idir=-1
       findax: do iax=1,nax
          call locate_condition(axarr(iax)%seqz,pcond,ceq)
          if(gx%bmperr.ne.0) goto 1000
          if(pcond%statev.lt.10) then
! this means intensive variable (T,P chemical potential)
             idir=-1; jax=iax; exit findax
          endif
       enddo findax
! both axis are extensive, take the first axix
       if(jax.eq.0) jax=1
       call locate_condition(axarr(jax)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'Searching for phase to fix along axis: ',jax
    endif
    call condition_value(1,pcond,curval,ceq)
!    write(*,'(a,3i4,F10.2)')'SMP2A initial value: ',gx%bmperr,jax,idir,curval
    if(gx%bmperr.ne.0) goto 1000
! it seems OK until here ....
    startval=curval
! increment axis variable using axinc and calculate with meq_sameset
100 continue
       curval=curval+idir*axarr(jax)%axinc
       call condition_value(0,pcond,curval,ceq)
!       write(*,'(a,i5,F12.5)')'SMP2A current value: ',gx%bmperr,curval
       if(gx%bmperr.ne.0) goto 1000
       irem=0; iadd=0; meqrec%noofits=0
!       write(*,*)'SMP2A calling meq_sameset from map_startline 1'
       call meq_sameset(irem,iadd,mapx,meqrec,meqrec%phr,inmap,ceq)
!       if(ocv()) write(*,110)'Search for phase change: ',&
!       write(*,110)'Search for phase change: ',&
!            idir*jax,gx%bmperr,irem,iadd,ceq%tpval(1),curval,axarr(jax)%axinc
110    format(a,i2,3i5,2x,F8.2,2(1pe12.4))
       maxtry=maxtry+1
       if(maxtry.gt.1000) then
          write(*,*)'SMP2A eternal loop: ',maxtry
          stop
       endif
       if(gx%bmperr.ne.0) goto 1000
       nophasechange: if(irem.eq.0 .and. iadd.eq.0) then
          if(idir.lt.0) then
             if(curval.le.axarr(jax)%axmin) then
! change direction
                idir=1
                curval=startval
             endif
          elseif(idir.gt.0) then
             if(curval.ge.axarr(jax)%axmax) then
                write(*,*)'No phase change along this axis'
                goto 1010
             endif
          endif
          goto 100
       endif nophasechange
!----------------------------------------------------------
! we found a phase to set fix!
    meqrec%nfixph=meqrec%nfixph+1
! This is written to handle several axis i.e. several fix phases.
!    write(*,*)'SMP2A found a phase change: ',irem,iadd
    fixfas: if(irem.gt.0) then
       if(meqrec%nstph.eq.1) then
          write(*,*)'Attempt to set the only phase as fix!'
          gx%bmperr=4230; goto 1000
       endif
!       write(*,*)'Remove axis condition and set stable phase fix: ',irem
! phase already in lists, just mark it is no fixed with zero amount
       meqrec%phr(irem)%stable=1
       meqrec%phr(irem)%curd%amfu=zero
       meqrec%phr(irem)%curd%dgm=zero
! set that the phase has fixed amount
       meqrec%phr(irem)%phasestatus=PHFIXED
       meqrec%fixph(1,meqrec%nfixph)=meqrec%phr(irem)%iph
       meqrec%fixph(2,meqrec%nfixph)=meqrec%phr(irem)%ics
       kph=irem
!---------------------------------------------------------------
    else !fixfas iadd
!       write(*,*)'SMP2A set new phase fix: ',iadd
       if(meqrec%nstph.eq.meqrec%maxsph) then
          write(*,*)'Too many phases stable',meqrec%maxsph
          gx%bmperr=4231; goto 1000
       endif
! copied from meq_phaseset
! the phase must be added in sequential order of phase and composition set no
       findplace: do jph=1,meqrec%nstph
          jj=meqrec%stphl(jph)
          if(meqrec%phr(iadd)%iph.gt.meqrec%phr(jj)%iph) then
             cycle
          endif
          if(meqrec%phr(iadd)%iph.lt.meqrec%phr(jj)%iph) then
             exit
          endif
! if same phase number compare composition set numbers
          if(meqrec%phr(iadd)%iph.eq.meqrec%phr(jj)%iph) then
             if(meqrec%phr(iadd)%ics.gt.meqrec%phr(jj)%ics) then
                cycle
             else
                exit
             endif
          endif
       enddo findplace
       do kph=meqrec%nstph,jph,-1
          meqrec%stphl(kph+1)=meqrec%stphl(kph)
       enddo
!       write(*,*)'SMP2A still trying to fix conditions ...'
! phase added at jph, (note jph may be equal to nstph+1)
       meqrec%stphl(jph)=iadd
       meqrec%nstph=meqrec%nstph+1
       meqrec%phr(iadd)%itadd=meqrec%noofits
       meqrec%phr(iadd)%curd%dgm=zero
       meqrec%phr(iadd)%stable=1
! set that the phase has fixed amount
       meqrec%phr(iadd)%phasestatus=PHFIXED
       meqrec%fixph(1,meqrec%nfixph)=meqrec%phr(iadd)%iph
       meqrec%fixph(2,meqrec%nfixph)=meqrec%phr(iadd)%ics
       kph=iadd
    endif fixfas
! meqrec%nfixph is used to reduce the number of variables in the system
! matrix.  Fix phases have no variable amount.
    meqrec%fixpham(meqrec%nfixph)=zero
!
!    write(*,*)'Now release axis condition: ',kph,pcond%active
! Must not forget to set if T or P is variable!
    pcond%active=1
    if(pcond%statev.eq.1) then
       meqrec%tpindep(1)=.TRUE.
    elseif(pcond%statev.eq.2) then
       meqrec%tpindep(2)=.TRUE.
    endif
! calling meq_sameset with iadd=-1 turn on verbose
    irem=0; iadd=0
!    write(*,*)'SMP2A calling meq_sameset from map_startline 2'
    call meq_sameset(irem,iadd,mapx,meqrec,meqrec%phr,inmap,ceq)
!    if(ocv()) write(*,110)'meq_sameset calculated: ',&
    if(gx%bmperr.gt.0) then
       write(*,*)'Failed to calculate with fix phase',gx%bmperr
       goto 1000
    elseif(iadd.gt.0 .or. irem.gt.0) then
       write(*,*)'Another phase want to be stable: ',iadd,irem
       gx%bmperr=4232; goto 1000
    endif
!    write(*,110)'SMP2A start calculated: ',0,gx%bmperr,irem,iadd,ceq%tpval(1)
!    if(gx%bmperr.ne.0) goto 1000
!
! we must return some values
!    write(*,*)'SMP2A now create start node and line equilibria'
! two exits
    ieq=2
! active axis, the remaining one, if jax=1 then 2, if jax=2 then 1
    axactive=3-jax
! templine is map_line record, some data must be set
    tmpline(1)%nfixphases=1
!    tmpline(1)%linefixph%phaseix=meqrec%phr(kph)%iph
    tmpline(1)%linefixph%ixphase=meqrec%phr(kph)%iph
    tmpline(1)%linefixph%compset=meqrec%phr(kph)%ics
    tmpline(1)%linefix_phr=kph
! allocate space for all stable phases minus one as fix, may already be alloc
! The number of stable phases can vary for different MAP commands
    if(allocated(tmpline(1)%stableph)) then
       deallocate(tmpline(1)%stableph)
       deallocate(tmpline(1)%stablepham)
       deallocate(tmpline(1)%stable_phr)
    endif
!    write(*,*)'map_startline: allocate 2: ',nstabphdim
    allocate(tmpline(1)%stableph(nstabphdim))
    allocate(tmpline(1)%stablepham(nstabphdim))
    allocate(tmpline(1)%stable_phr(nstabphdim))
    ll=0
    tmpline(1)%nstabph=0
    do jph=1,meqrec%nstph
       jj=meqrec%stphl(jph)
!       write(*,*)'Stable phase: ',meqrec%nstph,kph,jj
       if(jj.eq.kph) cycle
!       if(meqrec%phr(jj)%iph.eq.kph .and.&
!            meqrec%phr(jj)%ics.eq.kcs) cycle
!       write(*,66)'smp3: upper bound: ',jph,jj,size(tmpline(1)%stableph),&
!            nstabphdim,meqrec%nstph
66     format(a,10i4)
       ll=ll+1
!       write(*,*)'Store stable phase: ',jj,ll
       tmpline(1)%stableph(ll)%ixphase=meqrec%phr(jj)%iph
       tmpline(1)%stableph(ll)%compset=meqrec%phr(jj)%ics
       tmpline(1)%stablepham(ll)=meqrec%phr(jj)%curd%amfu
       tmpline(1)%stable_phr(ll)=jj
       tmpline(1)%nstabph=tmpline(1)%nstabph+1
! why exit?
!       exit
    enddo
!    if(ocv()) write(*,300)axactive,kph,tmpline(1)%linefixph%phaseix,&
!    write(*,300)axactive,kph,tmpline(1)%linefixph%ixphase,&
!         tmpline(1)%linefixph%compset,tmpline(1)%nstabph,&
!         (tmpline(1)%stableph(jj)%ixphase,tmpline(1)%stableph(jj)%compset,&
!         jj=1,tmpline(1)%nstabph)
300 format('exit map_startline: ',i2,i3,2x,2i3,2x,i2,10(2x,i3,i2))
    if(tmpline(1)%nstabph.eq.0) then
       write(*,*)'No stable phase !!'
       stop
    endif
1000 continue
    return
1010 continue
    gx%bmperr=4233
    goto 1000
  end subroutine map_startline

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_checkstep
!\begin{verbatim}
  subroutine map_checkstep(mapline,value,jj,axarr,nax,saveceq)
! check if step too large
! mapline is line record
! axarr is array with axis records
! nax is number of axis
! saveceq is record for saved equilibria
    implicit none
    integer nax
    type(map_line), pointer :: mapline
    type(map_axis), dimension(nax) :: axarr
    type(map_ceqresults), pointer :: saveceq
!\end{verbatim}
    integer place,jph,jj
    type(meq_setup), pointer :: meqrec
    type(gtp_state_variable), target :: axstv1
    type(gtp_state_variable), pointer :: axstv
    double precision value
    character ch1*1
    logical saveonfile
! pointer to last calculated (can be zero) and last free
! store last calulated axis values in axarr(iax)%lastaxval
!    write(*,*)'In map_checkstep',mapline%start%number_ofaxis,nax
!    do jj=1,nax
!       axstv1=axarr(jj)%axcond(1)
!       axstv=>axstv1
!       call state_variable_val(axstv,value,mapline%lineceq)
!       if(gx%bmperr.gt.0) goto 1000
!       if(nax.gt.1) then
! when several axis check if any has a big change ...
!    if(mapline%number_of_equilibria.gt.3) then
    if(abs(axarr(jj)%lastaxval-value).gt.&
         1.0D-1*(axarr(jj)%axmax-axarr(jj)%axmin)) then
       write(*,17)jj,mapline%axandir,mapline%number_of_equilibria,&
            axarr(jj)%lastaxval,value
17     format(' *** Too large change in axis: ',2i3,' at ',i4,&
            2(1pe14.6))
       gx%bmperr=4360; goto 1000
    endif
1000 continue
    return
  end subroutine map_checkstep

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_store
!\begin{verbatim}
  subroutine map_store(mapline,axarr,nax,saveceq)
! store a calculated equilibrium
! mapline is line record
! axarr is array with axis records
! nax is number of axis
! saveceq is record for saved equilibria
    implicit none
    integer nax
    type(map_line), pointer :: mapline
    type(map_axis), dimension(nax) :: axarr
    type(map_ceqresults), pointer :: saveceq
!\end{verbatim}
    integer place,jph,jj,lokcs
    type(meq_setup), pointer :: meqrec
    type(gtp_state_variable), target :: axstv1
    type(gtp_state_variable), pointer :: axstv
    double precision value
    character ch1*1
    logical saveonfile,testforspinodal
! pointer to last calculated (can be zero) and last free
! store last calulated axis values in axarr(iax)%lastaxval ALLOCATE
!    write(*,*)'SMP in map_store',gx%bmperr,globaldata%sysparam(2)
! insert a test for spinodal at every iii equilibriia
    testforspinodal=.FALSE.
    if(globaldata%sysparam(2).gt.0) then
       if(mod(mapline%number_of_equilibria,globaldata%sysparam(2)).eq.0) &
            testforspinodal=.TRUE.
    endif
!
    do jj=1,nax
       axstv1=axarr(jj)%axcond(1)
       axstv=>axstv1
       call state_variable_val(axstv,value,mapline%lineceq)
       if(gx%bmperr.gt.0) goto 1000
!       write(*,*)'map_store: ',value
! this check could be moved before store to take halfstep??
       if(nax.gt.1 .and. mapline%number_of_equilibria.gt.3) &
            call map_checkstep(mapline,value,jj,axarr,nax,saveceq)
       if(gx%bmperr.ne.0) goto 1000
!       if(nax.gt.1) then
! when several axis check if any has a big change ...
!          if(mapline%number_of_equilibria.gt.3) then
!             if(abs(axarr(jj)%lastaxval-value).gt.&
!                  1.0D-1*(axarr(jj)%axmax-axarr(jj)%axmin)) then
!                  2.0D-2*(axarr(jj)%axmax-axarr(jj)%axmin)) then
!                write(*,17)' *** map_store large step in axis: ',&
!                     mapline%number_of_equilibria,jj,&
!                     mapline%axandir,axarr(jj)%lastaxval,value
!17              format(a,3i3,2(1pe14.6))
!                gx%bmperr=4360; goto 1000
!                read(*,'(a)')ch1
!             endif
!          endif
!        endif
       axarr(jj)%lastaxval=value
    enddo
    if(repeatederr.ge.2) then
! VERY STRANGE BEHAVIOUR HERE, repeatederr not reset ??
! maybe not store if repeatederr nonzero 
       jj=repeatederr; repeatederr=0
!       write(*,*)'SMP in map_store',jj,repeatederr,gx%bmperr
! Finnaly I will store the calculated equilibrium but skip it for plotting
! if lasterr nonzero.
!       gx%bmperr=4399
!       goto 1000
    endif
    repeatederr=0
!    write(*,18)'stored: ',mapline%number_of_equilibria,(axarr(jj)%lastaxval,&
!         jj=1,mapline%start%number_ofaxis)
!18  format(a,i3,5(1pe14.6))
!-----------------------
    saveonfile=.FALSE.
! >>>> begin treadprotected
!    write(*,*)'map_store: ',saveonfile
    call reserve_saveceq(place,saveceq)
    if(gx%bmperr.eq.4219) then
! the memory is full, save this equilibrium, clean up and empty all on file
       saveonfile=.TRUE.
       gx%bmperr=0
    elseif(gx%bmperr.ne.0) then
! some other fatal error
       goto 1000
    endif
    if(repeatederr.gt.0) then
! maybe not store if repeatederr nonzero
!       write(*,*)'SMP in map_store',repeatederr,gx%bmperr,place
       repeatederr=0
    endif
!    write(*,*)'map_store: ',place,allocated(mapline%meqrec%phr)
!    write(*,*)'map_store: ',place,assigned(mapline%meqrec)
! >>>> end threadprotected
!-----------------------
! when step_tzero and some other step procedures MEQREC is not used
    if(.not.allocated(mapline%meqrec%phr)) goto 600
! loop through all phases and if their status is entered set it as PHENTUNST
! then loop through all stable to set status PHENTSTAB
! That is important for extracting values later ...
    meqrec=>mapline%meqrec
    do jph=1,meqrec%nphase
!          write(*,*)'phase and status: ',jph,meqrec%phr(jph)%curd%phstate,&
!               PHENTSTAB
!       if(meqrec%phr(jph)%curd%phstate.ge.PHENTUNST .and. &
!            meqrec%phr(jph)%curd%phstate.le.PHENTSTAB) then
!          meqrec%phr(jph)%curd%phstate=PHENTUNST
       if(meqrec%phr(jph)%curd%phstate.ge.PHENTUNST .and. &
            meqrec%phr(jph)%curd%phstate.le.PHENTERED) then
          meqrec%phr(jph)%curd%phstate=PHENTUNST
!       else
!          write(*,*)'map_store found a phase with status: ',&
!                         meqrec%phr(jph)%curd%phstate
       endif
    enddo
!    write(*,*)'map_store, stable phases',meqrec%nstph,place
    do jph=1,meqrec%nstph
       jj=meqrec%stphl(jph)
       if(meqrec%phr(jj)%curd%phstate.lt.PHFIXED) then
          meqrec%phr(jj)%curd%phstate=PHENTSTAB
! check if phase is inside miscibility gap
          if(testforspinodal) then
             lokcs=phasetuple(meqrec%phr(jj)%curd%phtupx)%compset
             call calc_qf(lokcs,value,mapline%lineceq)
             write(*,'(a,i3,F8.2,4(1pe12.4))')'SMP qf: ',lokcs,&
                  mapline%lineceq%tpval(1),value
             if(gx%bmperr.ne.0) then
                write(*,*)'SMP error chacking for instability',lokcs
                gx%bmperr=0
             elseif(value.lt.zero) then
                write(*,*)'SMP detected phase inside spinodal: ',lokcs
                gx%bmperr=4399; goto 1000
             endif
          endif
!       else
!          write(*,*)'Fix phase 1: ',jj,meqrec%phr(jj)%iph,meqrec%phr(jj)%ics
       endif
    enddo
!    write(*,201)' map store, stable phases: ',&
!         (meqrec%phr(meqrec%stphl(jj))%iph,&
!         meqrec%phr(meqrec%stphl(jj))%ics,jj=1,meqrec%nstph)
201 format(a,10(2i3,2x))
!-----------------------------------------
600 continue
! this copies the whole data structure !!!
! LIKELY PLACE FOR SEGMENTATION FAULT !!!
!    write(*,*)'SMP storing equilibrium record: ',place
    saveceq%savedceq(place)=mapline%lineceq
    saveceq%savedceq(place)%nexteq=0
    if(mapline%last.gt.0) then
       saveceq%savedceq(mapline%last)%nexteq=place
    endif
    mapline%last=place
    mapline%number_of_equilibria=mapline%number_of_equilibria+1
    if(mapline%first.eq.0) mapline%first=place
! this counter is zeroed when starting a new map/step unless old saved kept
    totalsavedceq=totalsavedceq+1
    if(totalsavedceq.gt.maxsavedceq) then
       write(kou,202)totalsavedceq
202    format(78('*')/'SMP saved equilibria overflow ',i5&
            ' and save on file is not implemented.'/&
            'Use smaller increments or reinitiate before STEP or MAP'/78('*')/)
       gx%bmperr=4219
    endif
    if(saveonfile) then
! We have to wind up all unfinished lines to continue step/map
! but this is not yet implemented
       write(*,207)
207    format(/' *** Buffer full and save on file not yet implemented,',&
            ' terminating step/map'/)
       gx%bmperr=4219
    endif
1000 continue
! nothing allocated?
!    write(*,*)'SMP exit map_store',place
    return
  end subroutine map_store

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_lineend
!\begin{verbatim}
  subroutine map_lineend(mapline,value,ceq)
! terminates gracefully a line at an axis limit or an error.
! mapline probably not needed except for testing
! value is last calculated axis value
! ceq is equilibrium record
    implicit none
    integer mode
    type(map_line), pointer :: mapline
    type(gtp_equilibrium_data), pointer :: ceq
    double precision value
!\end{verbatim}
!    type(meq_setup), pointer :: meqrec
! this will be called when a line ends at an axis limit, nothing to do?
    if(gx%bmperr.ne.0) then
       write(kou,75)mapline%number_of_equilibria,value,gx%bmperr
75     format('Terminating line with ',i4,' equilibria at ',1pe12.4,&
            ' due to error',i5)
       mapline%termerr=gx%bmperr
       gx%bmperr=0
! maybe do some cleanup ??
    else
       write(kou,77)mapline%number_of_equilibria,value
77     format('Terminating line with ',i4,' equilibria at axis limit ',1pe12.4)
       mapline%termerr=0
    endif
! mark there is no node at the end
    nullify(mapline%end)
1000 continue
! This routine should clear any error code
    gx%bmperr=0
    return
  end subroutine map_lineend

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_changeaxis
!\begin{verbatim}
  subroutine map_changeaxis(mapline,nyax,oldax,nax,axarr,axval,bytax,ceq)
! changes the axis with active condition to nyax
! mapline is line record
! nyax is index of new active axis
! oldax is index of old active axis
! nax is number of axis (always 2?)
! axarr is array with axis records
! axval the value to set as condition on new axis
! bytax logical, if true ignore axval ?? also used to indicate change of fix ph
! ceq is equilibrium record
    type(map_line), pointer :: mapline
    type(gtp_equilibrium_data), pointer :: ceq
    type(map_axis), dimension(nax) :: axarr
    logical bytax
    integer nyax,nax,oldax
    double precision axval
!\end{verbatim} %+
    type(gtp_condition), pointer :: pcond,lastcond
    type(gtp_state_variable), pointer :: axcondrec
    integer jax,iadd,irem,ierr,mapx
    double precision value
! turns off converge control for T
    integer, parameter :: inmap=1
! look for the condition record for new axis
!    write(*,*)'In map_changeaxis: ',nyax,axval
    call locate_condition(axarr(nyax)%seqz,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
!-----------
120 continue
! calculate the value of the inactive axis (nyax) condition.  An inactive
! condition is not updated automatically. Set prescribed value and
! activate the condition.  
    if(pcond%active.eq.0) then
       write(*,*)'In map_changeaxis, new axis condition is already acive!'
       goto 1000
    endif
    if(ocv()) write(*,*)'Axis condition: ',axarr(nyax)%axcond(1)%oldstv
!    svrrec=>pcond%statvar(1)
    axcondrec=>pcond%statvar(1)
!    axcondrec=>axarr(nyax)%axcond(1)
127 format('Map_changeaxis: ',2i2,2x,i3,2x,4i3,2x,i5,2x,2i5)
    call state_variable_val(axcondrec,value,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(ocv()) write(*,130)'New axis, current and prescribed: ',nyax,&
         value,axval,mapline%axvalx(nyax)
130 format(a,i2,2(1pe12.4))
! when called from bytaxis we should ignore current value ...
    if(bytax) then
       pcond%prescribed=axval
    else
       pcond%prescribed=value
    endif
    pcond%active=0
! we must indicate if T or P are now fixed ...
    if(pcond%statev.eq.1) then
       mapline%meqrec%tpindep(1)=.FALSE.
       if(ocv()) write(*,*)'Marking that T is variable'
    elseif(pcond%statev.eq.2) then
       mapline%meqrec%tpindep(2)=.FALSE.
    endif
!-------------------------------------------
! this is the old axis with active condition, look for its condition
    call locate_condition(axarr(oldax)%seqz,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(pcond%active.ne.0) then
       if(ocv()) write(*,*)'Wow, old axis condition is still active'
    else
! deactivate condition
       if(ocv())write(*,*)'Current value of old axis cond: ',pcond%prescribed
       pcond%active=1
    endif
! we must indicate if T or P are not fixed ...
    if(pcond%statev.eq.1) then
! in one case the value ceq%tpval was zero whereas the condition was positive
! This was due to a failed calculation of a new invariant equilibrium.
       mapline%meqrec%tpindep(1)=.TRUE.
       if(ocv()) write(*,*)'Marking that T is variable',ceq%tpval(1)
       ceq%tpval(1)=pcond%prescribed
    elseif(pcond%statev.eq.2) then
       mapline%meqrec%tpindep(2)=.TRUE.
    endif
!----------------------------------------------------------
! now we calculate the same equilibrium but with different axis condition!
    irem=0
    iadd=0
! add=-1 turn on verbose in meq_sameset
!    iadd=-1
!    if(bytax) then
!       write(*,*)'Calling meq_sameset in map_bytaxis:'
!       call list_conditions(kou,ceq)
!       iadd=-1
!    endif
    if(ocv()) write(*,*)'Map_changeaxis call meq_sameset, T=',ceq%tpval(1)
!    write(*,*)'SMP2A calling meq_sameset from map_changeaxis'
    call meq_sameset(irem,iadd,mapx,mapline%meqrec,mapline%meqrec%phr,inmap,ceq)
    if(gx%bmperr.ne.0) then
!       write(*,*)'Error from meq_sameset when trying to change axis',gx%bmperr
    endif
!       ierr=gx%bmperr; gx%bmperr=0
!       write(*,*)'Error trying to change axis: ',ierr
!       call list_conditions(kou,ceq)
!       gx%bmperr=ierr
!       if(ocv()) write(*,*)'Something really wrong ...',gx%bmperr,ceq%tpval(1)
!    else
!       write(*,990)gx%bmperr,irem,iadd,ceq%tpval(1)
!990    format(//' *** sucess *** ',3i5,1pe15.7//)
!    endif
!
1000 continue
    return
  end subroutine map_changeaxis

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_force_changeaxis
!\begin{verbatim} %-
  subroutine map_force_changeaxis(maptop,mapline,meqrec,nax,axarr,axvalok,ceq)
! force change of axis with active condition.  Works only with 2 axis.
! (and for tie-line not in plane ??).  Calls map_changeaxis ...
! maptop is node record
! mapline is line record
! meqrec is equilibrium calculation record
! nax is number of axis, also in maptop record
! axarr is array with axis records
! axvalok is last successfully calculated axis value
! ceq is equilibrium record
    implicit none
    integer nax
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(meq_setup) :: meqrec
    type(gtp_equilibrium_data), pointer :: ceq
    type(map_axis), dimension(*) :: axarr
    double precision axvalok
!\end{verbatim}
    double precision axfact,slope,xxx,value,axval,zzz
    integer nyax,seqz,jaxwc,oldax
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: svrrec
! copied from map_step
    if(ocv()) write(*,*)'Force change of axis with active condition: ',&
         mapline%axandir
! We have to change the axis with active condition, assume 2 axis
    jaxwc=abs(mapline%axandir)
    nyax=3-jaxwc
    oldax=jaxwc
    if(ocv()) write(*,101)mapline%number_of_equilibria,jaxwc,&
         nyax,xxx,mapline%axvals(oldax),ceq%tpval(1)
101 format('Bytaxis slope ',i3,2x,2i2,6(1pe12.4))
    axfact=1.0D-2
!
! Extract the current value of the (old) axis state variable items using pcond
    seqz=axarr(nyax)%seqz
    call locate_condition(seqz,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
    zzz=pcond%prescribed
    svrrec=>pcond%statvar(1)
    call state_variable_val(svrrec,axval,ceq)
    if(gx%bmperr.ne.0) goto 1000
! find the direction
    if(mapline%axvals(nyax)-mapline%axvalx(nyax).lt.0) then
! set negative direction and a small step
!       write(*,*)'Force_changeaxis 1: ',mapline%axandir,-nyax
       mapline%axandir=-nyax
!       xxx=mapline%axvals(nyax)-1.0D-2*axarr(nyax)%axinc
    else
! set positive direction
!       write(*,*)'Force_changeaxis 2: ',mapline%axandir,nyax
       mapline%axandir=nyax
!       xxx=mapline%axvals(nyax)+1.0D-2*axarr(nyax)%axinc
    endif
    xxx=zzz
    if(ocv()) write(*,62)'New axis direction: ',mapline%axandir,&
         mapline%axvals(nyax),mapline%axvalx(nyax)
62  format(a,i3,2x,2(1pe14.6))
! set new axis value as prescribed ... otherwise problems in map_changeaxis
    pcond%prescribed=xxx
    if(ocv()) write(*,63)'Call map_changeaxis',nyax,mapline%axchange,&
         mapline%number_of_equilibria,axval,zzz,xxx,ceq%tpval(1)
63  format(a,i2,2i3,4(1pe12.4))
!    call list_conditions(kou,ceq)
    call map_changeaxis(mapline,nyax,oldax,nax,axarr,xxx,.TRUE.,ceq)
    if(gx%bmperr.ne.0) then
!       seqz=gx%bmperr; gx%bmperr=0
!       write(*,*)'Error back from map_changeaxis: ',seqz
!       call list_conditions(kou,ceq)
!       gx%bmperr=seqz
       goto 1000
    endif
! change pcond!!!
    seqz=axarr(nyax)%seqz
    call locate_condition(seqz,pcond,ceq)
!    write(*,*)'After map_change: ',nyax,pcond%seqz,pcond%statev
    jaxwc=nyax
    mapline%axchange=mapline%number_of_equilibria
! value below is assumed to be most recently calculated value
    value=mapline%axvals(jaxwc)
    if(ocv()) write(*,16)'Axis, old and new condition: ',&
         mapline%axandir,value,xxx,ceq%tpval(1)
16  format(a,i3,6(1pe12.4))
! take a step in the axis variable.  mapline%axandir is +/-jaxwc
! mark axis changed
    mapline%axchange=mapline%number_of_equilibria
    if(mapline%axandir.gt.0) then
       value=value+axfact*axarr(jaxwc)%axinc
    else
       value=value-axfact*axarr(jaxwc)%axinc
    endif
    if(ocv()) write(*,202)'In map_step new, step & T: ',jaxwc,&
         mapline%axandir,value,axfact*axarr(jaxwc)%axinc,ceq%tpval(1)
202 format(a,2i3,3(1pe14.6))
    mapline%more=1
! Make sure value is set for the active axis condition!!
    seqz=axarr(jaxwc)%seqz
    call locate_condition(seqz,pcond,ceq)
! this call sets value as condition on the axis!
    if(ocv()) write(*,207)'Axis condition: ',jaxwc,pcond%statev,value
207 format(a,i2,i4,1pe14.6)
    call condition_value(0,pcond,value,ceq)
    if(gx%bmperr.ne.0) goto 1000
!
1000 continue
    return
  end subroutine map_force_changeaxis

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_step
!\begin{verbatim}
  subroutine map_step(maptop,mapline,meqrec,phr,axvalok,nax,axarr,ceq)
! select old or new step method
    implicit none
    integer nax
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(meq_setup) :: meqrec
    type(meq_phase), dimension(*), target :: phr
    type(gtp_equilibrium_data), pointer :: ceq
    type(map_axis), dimension(*) :: axarr
    double precision axvalok
!\end{verbatim}
! User can set GSOLDMAP 
! When not tielines inplane select old map
!    if(btest(globaldata%status,GSOLDMAP) .or. maptop%tieline_inplane.lt.0) then
    if(btest(globaldata%status,GSOLDMAP)) then
       call map_step_old(maptop,mapline,meqrec,phr,axvalok,nax,axarr,ceq)
    else
       call map_step2(maptop,mapline,meqrec,phr,axvalok,nax,axarr,ceq)
    endif
  end subroutine map_step

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_step_old
!\begin{verbatim}
  subroutine map_step_old(maptop,mapline,meqrec,phr,axvalok,nax,axarr,ceq)
! This is the OLD map_step routine used until 2018.01.31
! used also for map as mapping is stepping in one axis with fix phase condition
! calculate the next equilibrium along a line.  New phases can appear.
! axis with active condition can change and the direction.
! maptop is map node record
! mapline is line record
! phr is new array phase status (just for debugging)
! axvalok is last successfully calculated axis value
! nax number of axis, redundant as also in maptop record
! axarr is array with axis records
! ceq is equilibrium record
    implicit none
    integer nax
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(meq_setup) :: meqrec
    type(meq_phase), dimension(*), target :: phr
    type(gtp_equilibrium_data), pointer :: ceq
    type(map_axis), dimension(*) :: axarr
    double precision axvalok
!\end{verbatim}
    type(gtp_condition), pointer :: pcond
    integer seqz,jaxwc,jax,cmode,cmix(10),nyax,oldax,maybecongruent
    integer istv,indices(4),iref,iunit,ip,i1,i2,i3
    double precision value,dax1(5),dax2(5),axval(5),axval2(5),axvalt(5)
    double precision laxfact,xxx,yyy,bigincfix
    double precision preval(5),curval(5),prefixval(5),curfixval(5)
    double precision, parameter :: endfact=1.0D-6
    character ch1*1,statevar*24,encoded*24 
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(gtp_state_variable), target :: svrtarget
    logical tnip
!
!    write(*,*)'In map_step1: ',mapline%number_of_equilibria
!================================================================== new step
! tnip emergency to stop mapping outside limit for non-active axis
    tnip=.FALSE.
    laxfact=one
    maybecongruent=0
    axis: if(mapline%more.eq.0) then
! this means the current equilibrium is the last, line is terminated
       mapline%more=-1
       goto 1000
!================================================================== new step
! this is for STEP with one axis
    elseif(nax.eq.1) then
       seqz=axarr(1)%seqz
!       write(*,*)'Condition index: ',seqx
       call locate_condition(axarr(1)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
       call condition_value(1,pcond,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
! save last sucessfully calculated value in axvalok and axarr(1)%lastaxval
       axvalok=value
       axarr(1)%lastaxval=value
! good check point
       if(ocv()) write(*,16)'In map_step: ',mapline%number_of_equilibria,&
            mapline%axandir,value
16     format(a,2i3,6(1pe14.6))
       if(mapline%evenvalue.ne.zero) then
! If there is a value in mapline%evenvalue this is the first step in a new
! region, take 3 very small steps before using that as next value on axis!
          if(mapline%number_of_equilibria.lt.3) then
             value=value+1.0D-3*(mapline%evenvalue-value)
          else
             value=mapline%evenvalue
             mapline%evenvalue=zero
          endif
       else
! just take a step in axis variable.  mapline%axandir is -1 or +1
          value=value+axarr(1)%axinc*mapline%axandir
       endif
       mapline%more=1
       if(value.le.axarr(1)%axmin) then
          value=axarr(1)%axmin
! mapline%more=0 means this is the last calculation
          mapline%more=0
       elseif(value.ge.axarr(1)%axmax) then
          value=axarr(1)%axmax
          mapline%more=0
       endif
       call condition_value(0,pcond,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
! check conditions
!       call list_conditions(kou,ceq)
!       write(*,*)'New axis value: ',value,ceq%tpval(1)
!=============================================================== new step
! This is for MAP with 2 or more axis, both tie-line in plane and not
    else 
! at regular intervals check that phases with 2 or more composition sets have
! not identical constitutions!!
       if(mod(mapline%number_of_equilibria,3).eq.0) then
          call separate_constitutions(ceq)
       endif
! this is the current axis with acitive condition
       jaxwc=abs(mapline%axandir)
       bigincfix=one
!       write(*,*)'map_step: Number of fix phases: ',mapline%meqrec%nfixph
!       write(*,*)'map_step: Fix phase: ',mapline%meqrec%fixph(1,1),&
!            mapline%meqrec%fixph(2,1)
! Here we must compare changes in all axis to determine the axis for
! next step and how long step.  Last axis values stored in mapline%axvals
! Save previous currently in mapline%axvals in axval2
       nyax=0
       loopaxis: do jax=1,nax
          seqz=axarr(jax)%seqz
!          write(*,*)'Locating axis condition: ',seqz
          call locate_condition(seqz,pcond,ceq)
          if(gx%bmperr.ne.0) goto 1000
!          write(*,*)'Found axis condition'
          svrrec=>pcond%statvar(1)
          call state_variable_val(svrrec,axval(jax),ceq)
          if(gx%bmperr.ne.0) goto 1000
!          write(*,53)'Axis value: ',svrrec%oldstv,svrrec%argtyp,svrrec%phase,&
!               svrrec%compset,svrrec%component,axval(jax),mapline%axvals(jax)
53        format(a,5i4,2(1pe12.4))
          if(mapline%number_of_equilibria.eq.1) then
! for first equilibria just save the axis value
             mapline%axvals(jax)=axval(jax)
             laxfact=1.0D-2
          else
! for later equilibria calculate the slope
             preval(jax)=mapline%axvals(jax)
             curval(jax)=axval(jax)
             dax1(jax)=(axval(jax)-mapline%axvals(jax))/axarr(jax)%axinc
!             write(*,*)'dax1: ',jax,dax1(jax)
             axval2(jax)=mapline%axvals(jax)
             mapline%axvalx(jax)=mapline%axvals(jax)
             mapline%axvals(jax)=axval(jax)
          endif
!----------------------------- below tie-line in/not in plane separate new step
          tip1: if(maptop%tieline_inplane.gt.0) then
! if we have tie-lines in plane we must find the value of the axis condition
! for the fix phase or if it is a phase or component dependent state variable
             svrtarget=svrrec
             istv=svrrec%oldstv
             if(istv.ge.10) then
! in svrrec we have the axis variable for an extensive phase variable.  
! The value in mapline%axvals is for the entered phase, extract the value
! for the fix phase.  
! NOTE: If we change fix/entered phase we must change axvals/axvals2
!                i1=svr2%argtyp; i2=svr2%phase; i3=svr2%compset
                svrtarget%argtyp=3
                svrtarget%phase=mapline%linefixph(1)%ixphase
                svrtarget%compset=mapline%linefixph(1)%compset
             endif
! we must use a pointer in state_variable_val
             svr2=>svrtarget
             call state_variable_val(svr2,xxx,ceq)
             if(gx%bmperr.ne.0) goto 1000
             if(mapline%number_of_equilibria.eq.1) then
! for first equilibria just save the axisvalue for the fix phase
                mapline%axvals2(jax)=xxx
             else
! for later equilibria calculate the slope and check if close to limit
                dax2(jax)=(xxx-mapline%axvals2(jax))/axarr(jax)%axinc
                axvalt(jax)=mapline%axvals2(jax)
                if(jax.ne.jaxwc .and. istv.ge.10) then
                   prefixval(jax)=xxx
                   curfixval(jax)=mapline%axvals2(jax)
                   if(abs(prefixval(jax)-curfixval(jax)).gt.&
                        0.5D0*axarr(jax)%axinc) then
                      bigincfix=5.0D-1
                   endif
! for axis with inactive condition check if next step would pass min/max limit
! If so reduce the step in the active axis but do not change active axis!!
! xxx is last axis value, mapline%axvals2(jax) is previous
                   if(mapline%number_of_equilibria-mapline%axchange.gt.3) then
                      if(2*xxx-mapline%axvals2(jax).lt.axarr(jax)%axmin) then
                         nyax=jax
                    elseif(2*xxx-mapline%axvals2(jax).gt.axarr(jax)%axmax) then
                         nyax=jax
                      endif
                   endif
                   if(nyax.gt.0) then
!                      write(*,12)'Change nyax: ',nyax,&
!                           mapline%number_of_equilibria,curfixval(nyax),&
!                           curval(nyax)
12                    format(a,2i3,6(1pe12.4))
! This restriction needed to calculate two-phase regions with almost 
! verical lines (in T) and with one composition close to the axis limit
! and the other quite far away (like U4O9-GAS in U-O system)
! it should perhaps be refined to check that the lines are vertical ...
                      if(abs(curfixval(jax)-curval(jax)).gt.&
                           axarr(jax)%axinc) then
!                         write(*,*)'Ignore axis chnage!! ',nyax
                         nyax=0
                      endif
                   endif
                else
                   prefixval(jax)=xxx
                   curfixval(jax)=mapline%axvals2(jax)
! This test is very sensitive and if maybecongruent is set nonzero
! it is too much to reduce the step by 1.0D-2 below.  If so the map5
! fails at low T and I calculate too many points.  I set the
! reduction to 1.0D-1 which seems OK.
                   if(istv.ge.10 .and. &
                        abs(curval(jax)-curfixval(jax)).lt.&
                        axarr(jax)%axinc) then
!                        0.1*axarr(jax)%axinc) then
! if phase compositions are close decrease step!!
!                      write(*,93)'Phase compositions close:',jax,&
!                           mapline%number_of_equilibria,curval(jax),&
!                           curfixval(jax)
93                    format(a,2i5,4(1pe12.4))
                      maybecongruent=jax
                   endif
                endif
                mapline%axvals2(jax)=xxx
! check which change is the largest
!                if(ocv()) write(*,99)mapline%number_of_equilibria,jax,jaxwc,&
                write(*,99)mapline%number_of_equilibria,jax,jaxwc,&
                     nyax,mapline%axvals(jax),dax1(jax),&
                     mapline%axvals2(jax),dax2(jax)
99              format('Slope: ',i3,2x,3i2,6(1pe12.4))
             endif
             if(nyax.gt.0) then
!                write(*,*)'axis change due to limits: ',nyax,jaxwc
                mapline%axchange=mapline%number_of_equilibria
             endif
! here we can calculate the extrapolated values of both phases
! last calculated value a
             if(istv.ge.10) then
!                write(*,32)'stp xextra: ',jax,jaxwc,nyax,mapline%axandir,&
!                     curval(jax),preval(jax),curfixval(jax),prefixval(jax)
32              format(a,3i2,i3,6(1pe12.4))
             endif
          else
!------------------------------------------------------------
! here we have not tie-lines in the plane, we may need to change active axis
! and length of the step.
!             write(*,98)jax,axval(jax),mapline%axvals(jax),dax1(jax)
98           format('Tie-line not in plane, slope: ',i3,2x,6(1pe12.4))
! action to check if axis outside limit or slope requires axis change
! is done at the tip2 statement below
          endif tip1
       enddo loopaxis
!-------------------------------------------------------------
!       write(*,73)'Saved: ',(jax,mapline%axvalx(jax),&
!            mapline%axvals(jax),jax=1,nax)
!73     format(a,2(i4,2(1pe14.6)))
! dax1(jaxwc) is for active axis, if dax2(jaxw) is larger
! we should decrease the step length accordingly
       value=axval(jaxwc)
       if(mapline%number_of_equilibria.eq.1) then
! for the first step no slopes to check but take a very small step
          laxfact=1.0D-3
       else
          tip2: if(maptop%tieline_inplane.gt.0) then
! We have tielines in plane
! check if we should reduce axis step or change axis with active condition
             xxx=abs(dax2(jaxwc))
             if(nyax.eq.0) then
                nyax=jaxwc
                do jax=1,nax
                   if(jax.ne.jaxwc) then
! good check point
                      if(ocv()) write(*,33)jax,jaxwc,nyax,0,dax2(jax),xxx
!                      write(*,33)jax,jaxwc,nyax,0,dax2(jax),xxx
                     if(mapline%number_of_equilibria-mapline%axchange.gt.3) then
                         if(abs(dax2(jax)).gt.2*xxx) then
!           write(*,*)'Change active axis due to slope to/from: ',jax,jaxwc
                            xxx=abs(dax2(jax)); nyax=jax
                         endif
                      endif
                   endif
                enddo
33              format('Checking slopes: ',4i2,6(1pe12.4))
             endif
             if(nyax.ne.jaxwc) then
! We have to change the axis with active condition
                write(*,101)mapline%number_of_equilibria,jaxwc,&
                     nyax,xxx,dax1(1),dax2(1),dax1(2),dax2(2)
                if(ocv()) write(*,101)mapline%number_of_equilibria,jaxwc,&
                     nyax,xxx,mapline%axvals(1),dax1(1),&
                     mapline%axvals2(2),dax2(2)
101             format('Slope 3: ',i3,2x,2i2,6(1pe12.4))
! decrease the axis step factor
                mapline%axfact=1.0D-3
                oldax=abs(mapline%axandir)
                if(dax1(nyax).lt.0) then
! set negative direction and a small step
                   mapline%axandir=-nyax
                   xxx=mapline%axvals(nyax)-1.0D-2*axarr(nyax)%axinc
                else
! set positive direction and small step
                   mapline%axandir=nyax
                   xxx=mapline%axvals(nyax)+1.0D-2*axarr(nyax)%axinc
                endif
!                write(*,*)'axandir: ',nyax,mapline%axandir,xxx
                if(ocv()) write(*,63)'Call map_changeaxis',nyax,&
                     mapline%axchange,&
                     mapline%number_of_equilibria,dax1(nyax),dax2(nyax),xxx
63              format(a,i2,2i3,4(1pe12.4))
                call map_changeaxis(mapline,nyax,oldax,nax,axarr,xxx,&
                     .FALSE.,ceq)
                if(gx%bmperr.ne.0) goto 1000
! change pcond!!!
                seqz=axarr(nyax)%seqz
                call locate_condition(seqz,pcond,ceq)
                if(ocv()) write(*,*)'After map_change: ',&
                     nyax,pcond%seqz,pcond%statev
                jaxwc=nyax
                mapline%axchange=mapline%number_of_equilibria
! value below is assumed to be most recently calculated value
                value=mapline%axvals(jaxwc)
                if(ocv()) write(*,16)'Axis, old and new condition: ',&
                     mapline%axandir,value,xxx,ceq%tpval(1)
             endif
! 
!-----------------------------------------------------------------
          elseif(maptop%tieline_inplane.lt.0) then
! Tie-lines not in the plane
             do jax=1,nax
                if(jax.eq.jaxwc) cycle
! check if outside axis limit of non-active condition
                if(axval(jax).le.axarr(jax)%axmin) then
                   tnip=.TRUE.
                   write(kou,310)'Below ',jax,axval(jax),axarr(jax)%axmin
310                format(a,' limit',i3,2(1pe14.6),' of non-active axis')
                elseif(axval(jax).ge.axarr(jax)%axmax) then
                   tnip=.TRUE.
                   write(kou,310)'Above ',jax,axval(jax),axarr(jax)%axmax
                endif
! check if bytaxis when tie-lines not in plane
                if(abs(dax1(jax)).gt.one) then
                   write(*,*)'map_step_old: Change active axis: ',jax
                   call map_force_changeaxis(maptop,mapline,mapline%meqrec,&
                        nax,axarr,axvalok,ceq)
                   if(gx%bmperr.eq.0) goto 1000
                endif
             enddo
          endif tip2
       endif
!----------------------------------------------------------------------
! Here we decide the step to take in the axis variable.  
! mapline%axandir is +/-jaxwc
! laxfact takes into account if the fix phase changes more rapidly
! if maybecongruent is jaxwc then take small step
       i3=mapline%number_of_equilibria - mapline%axchange
       if(nax.gt.1) then
          if(i3.lt.3) then
! take small steps when starting a line or after axis change
             laxfact=1.0D-2
          elseif(i3.lt.6) then
             laxfact=1.0D-1
!          else
! laxfact= 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.0
!             laxfact=min(1.0,2.0D0*laxfact)
          endif
!          write(*,*)'stepcheck: ',nax,maybecongruent,i3
          if(maybecongruent.gt.0 .and. i3.ge.3) then
!             mapline%axfact=1.0D-2
             mapline%axfact=1.0D-1
!            write(*,*)'Decrease step due to close compositions',mapline%axfact
          endif
       endif
       axvalok=value
! laxfact is not saved between calls
! bigincfix 0.5 if fix phase changes more than 0.5*axinc
       bigincfix=one
       if(mapline%axandir.gt.0) then
          value=value+bigincfix*laxfact*axarr(jaxwc)%axinc*mapline%axfact
       else
          value=value-bigincfix*laxfact*axarr(jaxwc)%axinc*mapline%axfact
       endif
!       write(*,313)'laxfact: ',jaxwc,laxfact,value,&
!            axarr(jaxwc)%axinc,mapline%axfact,axvalok
313    format(a,i3,6(1pe12.4))
! good point for checking
       if(ocv()) write(*,65)'map_step: ',mapline%number_of_equilibria,&
            mapline%axandir,laxfact,mapline%axfact,ceq%tpval(1),axvalok,value
65     format(a,2i3,2(1pe10.2),4(1pe14.6))
       if(ocv()) write(*,202)'In map_step new, step & T: ',jaxwc,&
            mapline%axandir,value,laxfact*axarr(jaxwc)%axinc,ceq%tpval(1)
202    format(a,2i3,3(1pe14.6))
       if(mapline%axfact.lt.one) then
! calculation OK and no problems, make sure mapline%axfact approaches unity
!                   write(*,*)'Incrementing mapline%axfact: ',mapline%axfact
!          mapline%axfact=min(one,1.2D0*mapline%axfact)
! Trying to make axfact decrease less (like line above) makes map worse
!          mapline%axfact=min(one,2.0D0*mapline%axfact)
! factor above works well but sometimes too big increase
          mapline%axfact=min(one,1.5D0*mapline%axfact)
       endif
!======================================================================
! if the new axis value exceeds the min or max limit calculate for the limit
       mapline%more=1
       if(value.le.axarr(jaxwc)%axmin) then
          value=axarr(jaxwc)%axmin
! if a condition is an extensive variable like mole fraction avoid calculate
! for x(a)=0 or x(a)=1
          call locate_condition(axarr(jaxwc)%seqz,pcond,ceq)
          if(pcond%statev.gt.10) then
             value=value+endfact*axarr(jaxwc)%axinc
          endif
! mapline%more=0 means this is the last calculation ... At axis low limit
          write(kou,23)'low',value
23        format('At axis ',a,' limit',1pe12.4)
          mapline%more=0
       elseif(value.ge.axarr(jaxwc)%axmax) then
          value=axarr(jaxwc)%axmax
! if a condition is an extensive variable like mole fraction avoid calculate
! for x(a)=0 or x(a)=1 ........ at axis high limit
          call locate_condition(axarr(jaxwc)%seqz,pcond,ceq)
          if(pcond%statev.gt.10) then
             value=value-endfact*axarr(jaxwc)%axinc
          endif
          write(*,23)'high',value
          mapline%more=0
       endif
       if(ocv()) write(*,205)'Axis limits: ',mapline%more,axarr(jaxwc)%axmin,&
            value,axarr(jaxwc)%axmax
205    format(a,i2,3(1pe12.4))
! Make sure value is set for the active axis condition!!
       seqz=axarr(jaxwc)%seqz
       call locate_condition(seqz,pcond,ceq)
! this call sets value as condition on the axis!
       if(ocv()) write(*,207)'Axis condition: ',jaxwc,pcond%statev,value
!       write(*,207)'Axis condition: ',jaxwc,pcond%statev,value
207    format(a,i2,i4,1pe14.6)
       call condition_value(0,pcond,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
    endif axis
!------------------------------------------
1000 continue
! tnip set TRUE above if inactive axis outside limits and tie-line not in plane
    if(maptop%tieline_inplane.lt.0 .and. tnip) mapline%more=0
! if error code set mapline%more<0
    if(gx%bmperr.ne.0) mapline%more=-1
!    if(associated(pcond)) then
!       write(*,*)'Exit map_step: ',nyax,pcond%seqz,ceq%tpval(1)
!    endif
! To know which phase has nonzero amount
!    write(*,1001)'step_am: ',(mapline%meqrec%phr(ip)%curd%amfu,&
!         ip=1,mapline%meqrec%nphase),ceq%tpval(1)
1001 format(a,6(1pe12.4))
    return
  end subroutine map_step_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_step2
!\begin{verbatim}
  subroutine map_step2(maptop,mapline,meqrec,phr,axvalok,nax,axarr,ceq)
! used for map and step, mapping is to step with all but one axis replaced
! by a fix phase condition.  Map with tie-lines in plane special
! For map check if we should change independent (active) step axis.
! For tie-lines in plance check if we should change fix phase
! Set condition for the next equilibrium along the axis.  New phases can appear.
! axis with active condition can change and the direction.
! maptop is map node record
! mapline is line record
! phr is new array phase status (just for debugging)
! axvalok is last successfully calculated axis value
! nax number of axis, redundant as also in maptop record
! axarr is array with axis records
! ceq is equilibrium record
    implicit none
    integer nax
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(meq_setup) :: meqrec
    type(meq_phase), dimension(*), target :: phr
    type(gtp_equilibrium_data), pointer :: ceq
    type(map_axis), dimension(*) :: axarr
    double precision axvalok
!\end{verbatim}
    type(gtp_condition), pointer :: pcond
    integer seqz,jaxwc,jax,cmode,cmix(10),nyax,oldax,maybecongruent,mapeqno
    integer istv,indices(4),iref,iunit,ip,i1,i2,i3,jxxx
    double precision value,dax1(5),dax2(5),axval(5),axval2(5)
    double precision laxfact,xxx,yyy,maxstep
    double precision preval(5),curval(5),prefixval(5),curfixval(5)
    double precision, parameter :: endfact=1.0D-6
! trying to change step axis for mapping with tie-lines in plane
    integer fixbyte(2),twoextensiveaxis
    double precision isoent(2,2),isofix(2,2),isoe,isof,isofact
    double precision lastaxisvalue,stepfact
    character ch1*1,statevar*24,encoded*24
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(gtp_state_variable), target :: svrtarget
    logical tnip,nyfixph,ignore,approach
! new check for large step when tie-lines in the plane
    double precision ysave
    save maxstep,isofix,isoent,fixbyte,ignore,ysave,approach
!
    mapeqno=mapline%number_of_equilibria
! the dgm variables are for Al3N2 in the Al-Ni system which is not found stable
!   write(*,'(a,i4,F8.2)')'In map_step2: ',mapeqno,ceq%tpval(1)
!   write(*,'(a,i5,3i3,5(F10.2))')'In map_step2: ',mapeqno,meqrec%nv,&
!         maptop%tieline_inplane,mapline%axandir,ceq%tpval(1),&
!         ceq%phase_varres(3)%dgm,ceq%phase_varres(4)%dgm,&
!         ceq%phase_varres(5)%dgm
!    call list_conditions(kou,ceq)
    if(mapline%more.eq.0) then
! this means the current equilibrium is the last, line is terminated
       mapline%more=-1
       goto 1000
    endif
! tnip emergency to stop mapping outside limit for non-active axis
    tnip=.FALSE.
    laxfact=one
    twoextensiveaxis=0
    maybecongruent=0
! new global check for stable and metastable phases
!    if(maptop%globalcheckinterval.gt.0 .and. &
!         mod(mapeqno,maptop%globalcheckinterval).eq.0) then
! this may set error code if equilibrium should be recalculated
! and it may change constitutions of metastable phases
!       call check_all_phases(0,ceq)
!       if(gx%bmperr.ne.0) then
! these errors mean a new stable phase detected, we should terminate line
!          if(gx%bmperr.eq.4364 .or. gx%bmperr.eq.4365) goto 1000
! otherwise ignore any errors
!          gx%bmperr=0
!       endif
!    endif
    if(nax.eq.1) then
!================================================================== new step
! this is for STEP with one axis
       seqz=axarr(1)%seqz
!       write(*,*)'Condition index: ',seqx
       call locate_condition(axarr(1)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
       call condition_value(1,pcond,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
! save last sucessfully calculated value in axvalok and axarr(1)%lastaxval
       axvalok=value
       axarr(1)%lastaxval=value
! good check point
       if(ocv()) write(*,16)'In map_step: ',mapeqno,mapline%axandir,value
16     format(a,2i3,6(1pe14.6))
       if(mapline%evenvalue.ne.zero) then
! If there is a value in mapline%evenvalue this is the first steps in a new
! region, take 3 very small steps before using that as next value on axis!
          if(mapeqno.lt.3) then
             if(mapeqno.eq.1) then
                maxstep=mapline%evenvalue-value
!                write(*,*)'SMP maxstep: ',mapeqno,maxstep
             endif
             value=value+1.0D-3*(mapline%evenvalue-value)
          elseif(mapline%evenvalue.ne.zero .and. mapeqno.lt.6) then
! take a few more small steps ...
             value=value+2.0D-1*maxstep
          else
             value=mapline%evenvalue
             mapline%evenvalue=zero
          endif
       else
! just take a step in axis variable.  mapline%axandir is -1 or +1
          value=value+axarr(1)%axinc*mapline%axandir
       endif
!       write(*,*)'Next axis value: ',value
       mapline%more=1
       if(value.le.axarr(1)%axmin) then
          value=axarr(1)%axmin
! mapline%more=0 means this is the last calculation
          mapline%more=0
       elseif(value.ge.axarr(1)%axmax) then
          value=axarr(1)%axmax
          mapline%more=0
       endif
       call condition_value(0,pcond,value,ceq)
       goto 1000
!       if(gx%bmperr.ne.0) goto 1000
    endif
!=============================================================== new map
! This is for MAP with 2 or more axis, both tie-line in plane and not
    if(mod(mapeqno,3).eq.0) then
! at regulaar intervals check that phases with 2 or more composition sets have
! not identical constitutions!!  Should fix Cr-Fe metastable extrapolation!!
! It does not change anything for the stable phases
       call separate_constitutions(ceq)
    endif
! this is the current axis with acitive condition
    jaxwc=abs(mapline%axandir)
!    bigincfix=one
!       write(*,*)'map_step: Number of fix phases: ',mapline%meqrec%nfixph
!       write(*,*)'map_step: Fix phase: ',mapline%meqrec%fixph(1,1),&
!            mapline%meqrec%fixph(2,1)
! Here we must compare changes in all axis to determine the axis for
! next step and how long step.  Last axis values stored in mapline%axvals
! Save previous currently in mapline%axvals in axval2
    nyax=0
! isofact is to keep check of changes in fix phase when tie-lines in plane
    isofact=one
    loopaxis: do jax=1,nax
       seqz=axarr(jax)%seqz
       call locate_condition(seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
       svrrec=>pcond%statvar(1)
       call state_variable_val(svrrec,axval(jax),ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,53)'Axis value 1: ',svrrec%oldstv,svrrec%argtyp,svrrec%phase,&
!            svrrec%compset,svrrec%component,axval(jax),mapline%axvals(jax)
53     format(a,5i4,2(1pe12.4))
       if(mapeqno.eq.1) then
! for first equilibria just save the axis value
          approach=.true.
          mapline%axvals(jax)=axval(jax)
          laxfact=1.0D-2
!          isoent(2,jax)=axval(jax)
       else
! for later equilibria calculate the slope
          preval(jax)=mapline%axvals(jax)
          curval(jax)=axval(jax)
! CHECK CHANGE OF AXIS AND FIX PHASE HERE FOR ENTERED PHASE 1 of 3
          if(ocv()) write(*,94)'New and old axis values: ',mapeqno,jax,jaxwc,&
!          write(*,94)'New and old axis values: ',mapeqno,jax,jaxwc,&
               curval(jax),preval(jax),curval(jax)-preval(jax),&
               (curval(jax)-preval(jax))/axarr(jax)%axinc
94           format(a,i2,2x,2i2,2F10.4,2(1pe12.4))
          dax1(jax)=(axval(jax)-mapline%axvals(jax))/axarr(jax)%axinc
          axval2(jax)=mapline%axvals(jax)
          mapline%axvalx(jax)=mapline%axvals(jax)
          mapline%axvals(jax)=axval(jax)
!          isoent(1,jax)=isoent(2,jax)
!          isoent(2,jax)=axval(jax)
       endif
!----------------------------- below tie-line in/not in plane separate new step
       tip1: if(maptop%tieline_inplane.gt.0) then
! if we have tie-lines in plane we must find the value of the axis condition
! for the fix phase or if it is a phase or component dependent state variable
          svrtarget=svrrec
          istv=svrrec%oldstv
! istv>=10 means extensive condition (not potential)
          extvar: if(istv.ge.10) then
! in svrrec we have the axis variable for an extensive phase variable.  
! The value in mapline%axvals is for the entered phase, extract the value
! for the fix phase.  
! NOTE: If we change fix/entered phase we must change axvals/axvals2
             twoextensiveaxis=twoextensiveaxis+1
             ignore=.false.
             jxxx=jax
             svrtarget%argtyp=3
             svr2=>svrtarget
! extract composition of entered phase
             svrtarget%phase=mapline%stableph(1)%ixphase
             svrtarget%compset=mapline%stableph(1)%compset
! we must use a pointer in state_variable_val
             call state_variable_val(svr2,yyy,ceq)
             if(gx%bmperr.ne.0) goto 1000
! extract composition of fix phase
             svrtarget%phase=mapline%linefixph(1)%ixphase
             svrtarget%compset=mapline%linefixph(1)%compset
! we must use a pointer in state_variable_val
             svr2=>svrtarget
             call state_variable_val(svr2,xxx,ceq)
             if(gx%bmperr.ne.0) goto 1000
!             write(*,99)'Axis value 2: ',1,jax,jaxwc,0,xxx
             if(mapeqno.eq.1) then
! for first equilibria just save the axisvalue for the fix phase
                mapline%axvals2(jax)=xxx
                isofix(1,jax)=zero
                isofix(2,jax)=xxx
                isoent(1,jax)=zero
                isoent(2,jax)=yyy
                fixbyte(1)=mapline%linefixph(1)%ixphase
                fixbyte(2)=mapline%linefixph(1)%compset
             else
! for later equilibria calculate the slope and check if close to limit
! dax2 is slope value for fix phase
                isofix(1,jax)=isofix(2,jax)
                isofix(2,jax)=xxx
                isoent(1,jax)=isoent(2,jax)
                isoent(2,jax)=yyy
                if(fixbyte(1).ne.mapline%linefixph(1)%ixphase .and.&
                     fixbyte(2).ne.mapline%linefixph(1)%compset) then
                   ignore=.true.
                   fixbyte(1)=mapline%linefixph(1)%ixphase
                   fixbyte(2)=mapline%linefixph(1)%compset
                endif
                dax2(jax)=(xxx-mapline%axvals2(jax))/axarr(jax)%axinc
! CHECK CHANGE OF AXIS AND FIX PHASE HERE FOR FIX PHASE 2/3
!                if(ocv()) write(*,94)'Fix phase values:        ',&
!                write(*,94)'Fix phase values:   ',&
!                     mapeqno,jax,jaxwc,&
!                     xxx,mapline%axvals2(jax),xxx-mapline%axvals2(jax),&
!                     (xxx-mapline%axvals2(jax))/axarr(jax)%axinc
                mapline%axvals2(jax)=xxx
                if(jax.ne.jaxwc .and. istv.ge.10) then
                   prefixval(jax)=xxx
                   curfixval(jax)=mapline%axvals2(jax)
! for axis with inactive condition check if next step would pass min/max limit
! If so reduce the step in the active axis but do not change active axis!!
! xxx is last axis value, mapline%axvals2(jax) is previous
                   if(mapeqno-mapline%axchange.gt.3) then
                      if(2*xxx-mapline%axvals2(jax).lt.axarr(jax)%axmin) then
                         nyax=jax
                     elseif(2*xxx-mapline%axvals2(jax).gt.axarr(jax)%axmax) then
                         nyax=jax
                      endif
!                      if(nyax.gt.0) write(*,*)'SMP: change axis 1',isofact
                   endif
! nothing happends here ...
                   if(nyax.gt.0) then
! This restriction needed to calculate two-phase regions with almost 
! verical lines (in T) and with one composition close to the axis limit
! and the other quite far away (like U4O9-GAS in U-O system)
! it should perhaps be refined to check that the lines are vertical ...
                      if(abs(curfixval(jax)-curval(jax)).gt.&
                           axarr(jax)%axinc) then
!                      write(*,*)'Ignore axis change!! ',nyax
                         nyax=0
                      endif
                   endif
                else
                   prefixval(jax)=xxx
                   curfixval(jax)=mapline%axvals2(jax)
! This test is very sensitive and if maybecongruent is set nonzero
! it is too much to reduce the step by 1.0D-2 below.  If so the map5
! fails at low T and I calculate too many points.  I set the
! reduction to 1.0D-1 which seems OK.  istv>=10 means extensive variable
                   if(istv.ge.10 .and. &
                        abs(curval(jax)-curfixval(jax)).lt.&
                        axarr(jax)%axinc) then
! if phase compositions are close decrease step!!
93                    format(a,2i5,4(1pe12.4))
                      maybecongruent=jax
                   endif
                endif
                mapline%axvals2(jax)=xxx
! check which change is the largest
!             if(ocv()) write(*,99)'Slope: ',mapeqno,jax,jaxwc,&
!                write(*,97)'Slope: ',mapeqno,jax,jaxwc,&
!                     mapline%axvals(jax),dax1(jax),&
!                     mapline%axvals2(jax),dax2(jax)
!97              format(a,11x,i4,2i2,2(F10.4),2x,2(F10.4))
!                write(*,99)'Slope: ',mapeqno,jax,jaxwc,nyax,&
!                  entph_dxy(1,jax),entph_dxy(2,jax),&
!                  fixph_dxy(1,jax),fixph_dxy(2,jax)
99              format(a,i3,3i2,4(F10.4),2x,2(F10.4))
             endif
             if(nyax.gt.0) then
                mapline%axchange=mapline%number_of_equilibria
             endif
          else
! this axis is not extensive variable, same value as dax1(jax)
             dax2(jax)=dax1(jax)
          endif extvar
! end special for tie-lines in plane
       endif tip1
    enddo loopaxis
!-------------------------------------------------------------
! trying to avoid too big steps when two extensive axis variables
    if(twoextensiveaxis.eq.2) then
! UNFINISHED: this assumes both axis are compositions (fractions) !!!!!!!!!
! what about a composition and an enthalpy axis ??
       isoe=sqrt((isoent(2,1)-isoent(1,1))**2+(isoent(2,2)-isoent(1,2))**2)
       isof=sqrt((isofix(2,1)-isofix(1,1))**2+(isofix(2,2)-isofix(1,2))**2)
       if(plottrace) write(*,888)'smp1: ',&
            mapeqno,isoent(2,1),isoent(1,1),isoent(2,2),&
            isoent(1,2),isofix(2,1),isofix(1,1),isofix(2,2),isofix(1,2)
       if(mapeqno.gt.1) then
          i3=abs(mapline%axandir)
          if(plottrace) write(*,888)'smp2: ',mapline%axandir,isoe,isof,&
               axarr(i3)%axinc
888       format(a,i3,4F8.5,2x,4F8.5)
!          if(.not.ignore .and. isof.gt.3.0D0*isoe) then
          if(.not.ignore) then
! isofact is set to unity above
! THE TESTS HERE ARE QUITE CRAY BUT THEY WORK REASONABLY FOR
! MAP-10, BEF-500Y, CRFEMO(1400K), BEF-1500 and BEF-2500
             if(isoe.gt.2.0D0*axarr(i3)%axinc) then
! change in entered phase larger than max step
                isofact=axarr(i3)%axinc/isoe
!             elseif(isof.gt.3.0D0*axarr(i3)%axinc) then
             elseif(isof.gt.3.0D0*isoe) then
!             if(isof.gt.3.0D0*isoe) then
!             if(isof.gt.3.0D0*axarr(i3)%axinc) then
                isofact=isoe/isof
!                isofact=axarr(i3)%axinc/isof
             endif
             if(plottrace) write(*,'(a,3(1pe12.4))')'smp3: ',isoe,isof,isofact
          endif
       endif
    endif
!-------------------------------------------------------------
! for understanding what is happening ....
!    if(maptop%tieline_inplane.gt.0) then
!       write(*,59)'tieline: ',mapeqno,jaxwc,jxxx,nyax,&
!            mapline%stableph(1)%ixphase,mapline%linefixph(1)%ixphase,&
!            mapline%axvals(jxxx),mapline%axvals2(jxxx),&
!            mapline%axvals(3-jxxx),preval(jxxx),prefixval(jxxx)
!59          format(a,i4,3i2,2i3,2F10.5,f10.2,2(f10.5))
!    endif
! list last calculated tie-line
! we should check for the step length accordingly
    value=axval(jaxwc)
    if(mapeqno.eq.1) then
! for the first step no slopes to check but take a very small step
       laxfact=1.0D-3
    else
       tip2: if(maptop%tieline_inplane.gt.0) then
! We have tielines in plane
! check if we should reduce axis step or change axis with active condition
!          xxx=abs(dax2(jaxwc))
! xxx is set to the slope for the current independent axis and fix phase
          xxx=abs(dax1(jaxwc))
          nyfixph=.false.
!          write(*,*)'Attention 1: ',mapeqno,nyax,jaxwc
          if(nyax.eq.0) then
             nyax=jaxwc
             do jax=1,nax
                if(jax.ne.jaxwc) then
! good check point ?? YES
!                   write(*,33)mapeqno,jaxwc,jax,nyax,mapline%axandir,&
!                        meqrec%nv,&
!                        dax2(jax),xxx,dax1(jax),ceq%tpval(1)
33                 format('Check 7: ',6i3,6(1pe12.4))
! MISSING check for changing of fix/stable phase but keep same axis!!
                   if(mapeqno.gt.3 .and. mapeqno-mapline%axchange.gt.3) then
                      isotest1: if(isofact.eq.one) then
! ignore changing axis if isofact not unity
                         if(abs(dax2(jax)).gt.2*xxx) then
! dependent axis changes more! change independent axis
                            xxx=abs(dax2(jax))
58                          format(a,2i3,2(1pe12.4))
                            nyfixph=.true.
                            nyax=jax
                         elseif(abs(dax1(jax)).gt.2*xxx) then
                            xxx=abs(dax1(jax))
                            nyax=jax
                         endif
                      endif isotest1
                   endif
                else
! if the independent axis is extensive check if we should change fix phase
                   seqz=axarr(jax)%seqz
                   call locate_condition(seqz,pcond,ceq)
                   if(gx%bmperr.ne.0) goto 1000
                   svrrec=>pcond%statvar(1)
!                   call state_variable_val(svrrec,axval(jax),ceq)
!                   if(gx%bmperr.ne.0) goto 1000
! If independent axis is an extensive variable check for fix phase change
! This does not seem to change anything!!!
                   if(svrrec%oldstv.ge.10) then
                      if(mapeqno-mapline%axchange.gt.3 .and. &
                           abs(dax2(jax)).gt.abs(dax1(jax))) then
! dependent axis for fix phase changes more, change axis and fix phase!
                         nyfixph=.true.
 !                        write(*,101)'Change fix phase?',mapeqno,jaxwc,&
 !                             nyax,mapline%linefixph(1)%ixphase,&
 !                             mapline%stableph(1)%ixphase,nyfixph,&
 !                             dax2(jax),dax1(jax)
                      endif
                   endif
                endif
             enddo
          endif
! This is all for tie-lines in the plane!!
!          if(nyax.ne.jaxwc) write(*,*)'Attention 2: ',mapeqno,nyax,jaxwc
!          write(*,152)'Attention 2: ',mapeqno,nyax,jaxwc,.FALSE.,&
!               dax1(nyax),dax2(nyax),dax1(3-nyax),dax2(3-nyax),ceq%tpval(1)
          limits: if(nyax.eq.jaxwc .and. jaxwc.ne.jxxx .and. &
               mapeqno-mapline%axchange.gt.3 .and. .not.nyfixph) then
! Problems in U-O system with gas and U3O8 when gas is almost pure O
! If the entered (not fixed) phase cannot vary its composition 
! that is bad but do nothing here
             if(fixedcomposition(mapline%stableph(1)%ixphase)) then
!                write(*,*)'Continue as entered phase has fixed composition!'
                exit limits
             endif
! check if phase compositions are close
             if(abs(mapline%axvals(jxxx)-mapline%axvals2(jxxx)).gt.&
                  axarr(jxxx)%axinc) then
!                write(*,69)'Continue as phase compositions not close',&
!                     mapline%axvals(jxxx),mapline%axvals2(jxxx)
69              format(a,2F10.6)
! They are not ... do nothing
                exit limits
             endif
! No changes, check if we are close to the end of the extensive variable axis
             if(2*mapline%axvals(jxxx)-preval(jxxx).gt.&
                  axarr(jxxx)%axmax) then
!                write(*,91)'high',jxxx,2*mapline%axvals(jxxx)-preval(jxxx)
                nyax=jxxx
!                write(*,*)'SMP nyax 4:',nyax
91              format('At ',a,' limit, change axis to: ',i2,F10.6)
             elseif(2*mapline%axvals(jxxx)-preval(jxxx).lt.&
                  axarr(jxxx)%axmin) then
!                write(*,91)'low',jxxx,2*mapline%axvals(jxxx)-preval(jxxx)
                nyax=jxxx
!                write(*,*)'SMP nyax 5:',nyax
             endif
          endif limits
!
!          write(*,152)'Attention 3: ',mapeqno,nyax,jaxwc,nyfixph,&
!               dax1(nyax),dax2(nyax),dax1(3-nyax),dax2(3-nyax),ceq%tpval(1)
          newaxis: if(nyax.ne.jaxwc) then
! We have to change the axis with active condition
!             write(*,101)'Slope 3: ',mapeqno,jaxwc,nyax,&
!                  mapline%linefixph(1)%ixphase,&
!                  mapline%stableph(1)%ixphase,nyfixph,&
!                  mapline%axvals(nyax),mapline%axvals2(nyax),&
!                  dax1(nyax),dax2(nyax)
101          format(a,5i3,l2,6(1pe12.4))
! decrease the axis step factor
             mapline%axfact=1.0D-3
             oldax=abs(mapline%axandir)
! emergency fix: if dax1(nyax) is zero we must change fix phase!
             if(dax1(nyax).eq.zero .and. .not.nyfixph) nyfixph=.TRUE.
!             write(*,152)'SMP: change active axis: ',nyax,mapline%axandir,&
!                  jaxwc,nyfixph,dax1(nyax),dax2(nyax),&
!                  dax1(3-nyax),dax2(3-nyax),ceq%tpval(1)
152          format(a,3i3,l2,5(1pe10.2))
             if(nyfixph) then
! We must set new fix phase, take the direction from dax2
                if(dax2(nyax).lt.0) then
! set negative direction and a small step
                   mapline%axandir=-nyax
                   xxx=mapline%axvals2(nyax)-1.0D-2*axarr(nyax)%axinc
                else
! set positive direction and small step
                   mapline%axandir=nyax
                   xxx=mapline%axvals2(nyax)+1.0D-2*axarr(nyax)%axinc
                endif
             else
                if(dax1(nyax).lt.0) then
! set negative direction and a small step
                   mapline%axandir=-nyax
                   xxx=mapline%axvals(nyax)-1.0D-2*axarr(nyax)%axinc
                else
! set positive direction and small step
                   mapline%axandir=nyax
                   xxx=mapline%axvals(nyax)+1.0D-2*axarr(nyax)%axinc
                endif
             endif
             if(ocv()) write(*,63)'Call map_changeaxis',nyax,&
                  mapline%axchange,&
                  mapeqno,dax1(nyax),dax2(nyax),xxx
63           format(a,i2,2i3,4(1pe12.4))
!  bytax is TRUE if axval is new axis condition
!             if(nyfixph) then
!                call list_conditions(kou,ceq)
!             endif
             if(nyfixph) then
! This routine switches the fix and entered phases
                if(plottrace) write(*,*)'new fix phase: ',nyfixph
                call map_bytfixphase(mapline,oldax,meqrec,xxx,ceq)
                if(gx%bmperr.ne.0) goto 1000
                ignore=.TRUE.
             endif
!             write(*,*)'New independent axis and value: ',nyax,xxx,nyfixph
             call map_changeaxis(mapline,nyax,oldax,nax,axarr,xxx,&
                  nyfixph,ceq)
!             call map_changeaxis(mapline,nyax,oldax,nax,axarr,xxx,&
!                     .FALSE.,ceq)
!             if(nyfixph) then
!                call list_conditions(kou,ceq)
!                write(*,*)'new fix phase ',mapline%axandir,ceq%tpval(1)
!                read(*,62)ch1
!62              format(a)
!             endif
             if(gx%bmperr.ne.0) goto 1000
! change pcond!!!
             seqz=axarr(nyax)%seqz
             call locate_condition(seqz,pcond,ceq)
             if(ocv()) write(*,*)'After map_change: ',&
                  nyax,pcond%seqz,pcond%statev
             jaxwc=nyax
             mapline%axchange=mapline%number_of_equilibria
! value below is assumed to be most recently calculated value
             value=mapline%axvals(jaxwc)
             if(ocv()) write(*,16)'Axis, old and new condition: ',&
                  mapline%axandir,value,xxx,ceq%tpval(1)
          endif newaxis
! 
!-----------------------------------------------------------------
       elseif(maptop%tieline_inplane.lt.0) then
! Tie-lines not in the plane
          do jax=1,nax
             if(jax.eq.jaxwc) cycle
! check if outside axis limit of non-active condition
             if(axval(jax).le.axarr(jax)%axmin) then
                tnip=.TRUE.
                write(kou,310)'Below ',jax,axval(jax),axarr(jax)%axmin
310             format(a,' limit',i3,2(1pe14.6),' of non-active axis')
             elseif(axval(jax).ge.axarr(jax)%axmax) then
                tnip=.TRUE.
                write(kou,310)'Above ',jax,axval(jax),axarr(jax)%axmax
             endif
! check if bytaxis when tie-lines not in plane
             if(abs(dax1(jax)).gt.one) then
!                write(*,*)'map_step: Change active axis: ',jax
                call map_force_changeaxis(maptop,mapline,mapline%meqrec,&
                     nax,axarr,axvalok,ceq)
                if(gx%bmperr.eq.0) goto 1000
             endif
          enddo
! end check for tie-lines in plane
       endif tip2
    endif
!----------------------------------------------------------------------
! Here we decide the step to take in the axis variable.  
! mapline%axandir is +/-jaxwc
! laxfact takes into account if the fix phase changes more rapidly
! if maybecongruent is jaxwc then take small step
    i3=mapline%number_of_equilibria - mapline%axchange
    if(nax.gt.1) then
       if(i3.lt.3) then
! take small steps when starting a line or after axis change
          laxfact=1.0D-2
       elseif(i3.lt.6) then
          laxfact=1.0D-1
       endif
       if(maybecongruent.gt.0 .and. i3.ge.3) then
          mapline%axfact=1.0D-1
       endif
    endif
    axvalok=value
! laxfact is not saved between calls
! bigincfix 0.5 if fix phase changes more than 0.5*axinc
!    bigincfix=one
    lastaxisvalue=value
    if(mapline%axandir.gt.0) then
       value=value+isofact*laxfact*axarr(jaxwc)%axinc*mapline%axfact
    else
       value=value-isofact*laxfact*axarr(jaxwc)%axinc*mapline%axfact
    endif
! good point for checking
    if(ocv()) write(*,65)'map_step: ',mapeqno,&
         mapline%axandir,laxfact,mapline%axfact,ceq%tpval(1),axvalok,value
65  format(a,2i3,2(1pe10.2),4(1pe14.6))
    if(ocv()) write(*,202)'In map_step new, step & T: ',jaxwc,&
         mapline%axandir,value,laxfact*axarr(jaxwc)%axinc,ceq%tpval(1)
202 format(a,2i3,3(1pe14.6))
    if(mapline%axfact.lt.one) then
! calculation OK and no problems, make sure mapline%axfact approaches unity
!                   write(*,*)'Incrementing mapline%axfact: ',mapline%axfact
!          mapline%axfact=min(one,1.2D0*mapline%axfact)
! Trying to make axfact decrease less (like line above) makes map worse
       mapline%axfact=min(one,2.0D0*mapline%axfact)
    endif
!======================================================================
! if the new axis value exceeds the min or max limit calculate for the limit
    mapline%more=1
    if(value.le.axarr(jaxwc)%axmin) then
       value=axarr(jaxwc)%axmin
! if a condition is an extensive variable like mole fraction avoid calculate
! for x(a)=0 or x(a)=1
       call locate_condition(axarr(jaxwc)%seqz,pcond,ceq)
       if(pcond%statev.gt.10) then
          value=value+endfact*axarr(jaxwc)%axinc
       endif
! mapline%more=0 means this is the last calculation
       write(kou,23)'low',value
23     format('At axis ',a,' limit',1pe12.4)
       mapline%more=0
    elseif(value.ge.axarr(jaxwc)%axmax) then
       value=axarr(jaxwc)%axmax
! if a condition is an extensive variable like mole fraction avoid calculate
! for x(a)=0 or x(a)=1
       call locate_condition(axarr(jaxwc)%seqz,pcond,ceq)
       if(pcond%statev.gt.10) then
          value=value-endfact*axarr(jaxwc)%axinc
       endif
       write(*,23)'high',value
       mapline%more=0
    endif
!....... special for axis limits of isothermal sections DOES NOT WORK
! check if we are close to an axis limit for isothermal sections
    if(mapeqno.gt.2 .and. twoextensiveaxis.eq.2) then
! The fraction of the third component of entered phase (where we step):
       call locate_condition(axarr(jaxwc)%seqz,pcond,ceq)
       if(pcond%statev.gt.10) then
          xxx=pcond%prescribed
       endif
       yyy=one-isoent(2,jaxwc)-isoent(2,3-jaxwc)
       if(yyy.le.0.5D0*axarr(jaxwc)%axinc) then
! changing the axis variable will make third fraction negative       
! we should decrease value ...
!          write(*,'(a,i3,F9.5,7F8.5)')'At boundary? ',mapeqno,yyy,&
!               isoent(2,jaxwc),isoent(2,3-jaxwc),&
!               isofix(2,jaxwc),isofix(2,3-jaxwc),xxx,value,value-xxx
          if(approach) then
             write(*,*)'Approaching limit of third component'
             approach=.false.
          endif
          if(yyy.gt.zero) then
             if(yyy.lt.ysave) then
                value=xxx+0.9D0*yyy
             endif
             ysave=yyy
          else
             write(*,*)'Impossible!',yyy
          endif
       endif
    endif
!......
    if(ocv()) write(*,205)'Axis limits: ',mapline%more,axarr(jaxwc)%axmin,&
         value,axarr(jaxwc)%axmax
205 format(a,i2,3(1pe12.4))
! Make sure value is set for the active axis condition!!
    seqz=axarr(jaxwc)%seqz
    call locate_condition(seqz,pcond,ceq)
! CHECK CHANGE OF AXIS AND FIX PHASE HERE 3/3
    if(ocv()) write(*,207)'New axis condition: ',jaxwc,pcond%statev,value,&
         value-lastaxisvalue
207 format(a,i2,i4,2(1pe14.6))
    call condition_value(0,pcond,value,ceq)
    if(gx%bmperr.ne.0) goto 1000
!------------------------------------------
1000 continue
! tnip set TRUE above if inactive axis outside limits and tie-line not in plane
    if(maptop%tieline_inplane.lt.0 .and. tnip) mapline%more=0
! if error code set mapline%more<0
    if(gx%bmperr.ne.0) mapline%more=-1
!    if(associated(pcond)) then
!       write(*,*)'Exit map_step: ',nyax,pcond%seqz,ceq%tpval(1)
!    endif
! To know which phase has nonzero amount
!    write(*,1001)'step_am: ',(mapline%meqrec%phr(ip)%curd%amfu,&
!         ip=1,mapline%meqrec%nphase),ceq%tpval(1)
1001 format(a,6(1pe12.4))
!    write(*,*)'Leaving map_step2 '
    return
! axis limit
  end subroutine map_step2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_bytfixphase
!\begin{verbatim}
  subroutine map_bytfixphase(mapline,axis,meqrec,xxx,ceq)
! Try to change the fix phase for axis
! the new axis value is in xxx (not needed??)
! mapline is map line record
! ceq is equilibrium record
    implicit none
    type(map_line), pointer :: mapline
    integer axis
    type(meq_setup) :: meqrec
    double precision xxx
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! REMEMBER: %stableph(..) and %linefixph are arrays of phase tuples !!!
! 5 integers: lokph,compset,ixphase,lokvares,nextcs, only ixphase/compset set
! we need meqrec!!!
    type(gtp_phasetuple) :: phtup1
    integer lokcs,phrix
    double precision x1,x2
! just as check
    x1=mapline%stablepham(1)
!    write(*,33)'Phase change 1:',meqrec%fixph(1,1),meqrec%fixph(2,1),&
!         meqrec%iphl(1),meqrec%icsl(1),meqrec%iphl(2),meqrec%icsl(2),&
!         meqrec%aphl(1),meqrec%aphl(2),xxx
33  format(a,3(i3,i2,2x),3F8.3)
! we must change in meqrec also!! This is for tie-lines in plane,
! one fix phase, one stable phase
    phtup1=mapline%linefixph(1)
    phrix=mapline%linefix_phr(1)
    mapline%linefixph(1)=mapline%stableph(1)
    mapline%linefix_phr(1)=mapline%stable_phr(1)
    if(meqrec%nfixph.ne.1) then
       write(*,*)'MAP wants to change ONE fix phase: ',meqrec%nfixph
       gx%bmperr=4399; goto 1000
    endif
    meqrec%fixph(1,1)=mapline%linefixph(1)%ixphase
    meqrec%fixph(2,1)=mapline%linefixph(1)%compset
    meqrec%fixpham(1)=zero
    meqrec%iphl(1)=mapline%linefixph(1)%ixphase
    meqrec%icsl(1)=mapline%linefixph(1)%compset
!------------- now the stable phase  ?? value of stable_phr=??
    mapline%stableph(1)=phtup1
    mapline%stable_phr(1)=phrix
!    write(*,*)'SMP2A phrix switching fix/stable phase: ',phrix
! nstabph is part of mapfix record ... saved in meqrec%nv
! we are not changing the number of fix or stable phases ...
    if(meqrec%nv.ne.2) then
       write(*,*)'MAP wants to change ONE stable phase: ',meqrec%nv
       gx%bmperr=4399; goto 1000
    endif
    meqrec%iphl(2)=phtup1%ixphase
    meqrec%icsl(2)=phtup1%compset
    meqrec%aphl(2)=x1
! we have changed the stable phase, set a positive amount
    mapline%stablepham(1)=x1
!    write(*,33)'Phase change 2:',meqrec%fixph(1,1),meqrec%fixph(2,1),&
!         meqrec%iphl(1),meqrec%icsl(1),meqrec%iphl(2),meqrec%icsl(2),&
!         meqrec%aphl(1),meqrec%aphl(2),xxx
1000 continue
    return
  end subroutine map_bytfixphase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_calcnode
!\begin{verbatim}
  subroutine map_calcnode(irem,iadd,maptop,mapline,meqrec,axarr,ceq)
! we have found a change in the set of stable phases.  check if this node
! already been found and if so eliminate a line record.  Otherwise 
! create a new node record with line records and continue mapping one
! of these.
! irem and iadd are indices (in phr?) of phase that will disappear/appear
! maptop is map node record
! mapline is map line record
! meqrec is equilibrium calculation record, ! Note changes in meqrec is local,
!      not copied to mapline%meqrec!!!
! axarr is array with axis records
! ceq is equilibrium record
    implicit none
    integer irem,iadd
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(meq_setup) :: meqrec
    type(map_axis), dimension(*) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    type(gtp_condition), pointer :: lastcond,pcond
    integer iremsave,iaddsave,iph,ics,jj,jph,kph,phfix,seqx,jax,haha
    integer what,type,cmix(10),maxstph,noplot,mode,addtupleindex,mapx,sameadd
    double precision, parameter :: addedphase_amount=1.0D-2
    double precision value,axval,axvalsave,tx,nodefixpham
    type(gtp_state_variable), pointer :: svrrec
    logical global
    double precision, allocatable, dimension(:) :: yfra
    type(gtp_equilibrium_data), target :: tceq
    type(gtp_equilibrium_data), pointer :: pceq
    character phname*32
! turns off converge control for T
    integer, parameter :: inmap=1
!
!    write(*,*)'In map_calcnode phase change add/remove: ',iadd,irem
! we have already called same_composition(iadd...)
    iremsave=irem
    iaddsave=iadd
    if(irem.gt.0) then
       if(iadd.gt.0) then
          write(*,*)'Confusion, both add and remove phases?'
! check if phase to be added is already stable
          if(same_composition(iadd,meqrec%phr,meqrec,ceq,zero)) then
             iadd=0
             phfix=-irem
          else
! go back and calculate with half the step ... 
             gx%bmperr=4220; goto 1000
          endif
       else
          phfix=-irem
          iadd=irem
       endif
    else
       phfix=iadd
    endif
!--------------------------------------------
! remove here the axis condition, abs(mapline%axandir) gives active axis
! axandir is the axis with active condition.  It can be negative
    jax=abs(mapline%axandir)
!    write(*,*)'Remove axis condition: ',jax,axarr(jax)%seqz
    lastcond=>ceq%lastcondition
    if(.not.associated(lastcond)) then
       write(*,*)'in map_calcnode, no conditions: ',jax
       gx%bmperr=4221; goto 1000
    endif
    pcond=>lastcond
60  continue
    pcond=>pcond%next
    if(pcond%seqz.eq.axarr(jax)%seqz) goto 70
    if(.not.associated(pcond,lastcond)) goto 60
    write(*,*)'in map_calcnode the axis condition not found: ',jax
    gx%bmperr=4221; goto 1000
!
70  continue
! this removes the condition, remember pcond as condition will be set again!!
    pcond%active=1
    axval=pcond%prescribed
    if(ocv()) write(*,77)pcond%seqz,pcond%prescribed,ceq%tpval(1),axval
77  format('Removing condition: ',i3,6(1pe12.4))
! if the condition is T or P this must be indicated specially
! if a potential condition released we can have one more stable phse
    maxstph=0
    if(pcond%statev.eq.1) then
       meqrec%tpindep(1)=.TRUE.
       if(ocv()) write(*,*)'Marking that T is variable'
       maxstph=1
    elseif(pcond%statev.eq.2) then
       meqrec%tpindep(2)=.TRUE.
       maxstph=1
    endif
!--------------------------------------------
! independently if iadd or irem >0 set this phase, phfix, fix with zero amount
! we may return here if there is problems calculate the node equilibrium
100 continue
! set phfix fix with amount nodefixpham
    nodefixpham=zero
! NOTE it must be added so meqrec%stphl in ascending order
    if(phfix.gt.0 .and. meqrec%nstph.eq.meqrec%maxsph+maxstph) then
! No more phases allowed, we must see if  some other phase may be removed
!       write(*,*)'Too many stable phases at nodepoint',meqrec%maxsph
! set back pcond as active, this saved top of miscibility gap in Cr-Mo !!!
       pcond%active=0
!       if(same_composition(iadd,meqrec%phr,meqrec,ceq,zero)) then
!          iadd=0; goto 201
!       endif
!       write(*,'(a,10i5)')'SMP node with too many stable phases: ',&
!            iremsave,iaddsave,phfix,meqrec%nstph,maxstph
       gx%bmperr=4223; goto 1000
    else
       if(ocv()) write(*,*)'Number of stable phases at nodepoint',&
            meqrec%nstph,maxstph
    endif
    if(phfix.gt.0) then
! the phase must be added in sequential order of phase and composition set no
       findplace: do jph=1,meqrec%nstph
          jj=meqrec%stphl(jph)
          if(meqrec%phr(phfix)%iph.gt.meqrec%phr(jj)%iph) then
             cycle
          endif
          if(meqrec%phr(phfix)%iph.lt.meqrec%phr(jj)%iph) then
             exit
          endif
! if same phase number compare composition set numbers
          if(meqrec%phr(phfix)%iph.eq.meqrec%phr(jj)%iph) then
             if(meqrec%phr(phfix)%ics.gt.meqrec%phr(jj)%ics) then
                cycle
             else
                exit
             endif
          endif
       enddo findplace
! one should come here at exit, iadd should be inserted before 
! meqrec%stphl(jph), jph can be nstph+1 if added phase should be the last
! otherwise shift previous phases one step up.
       do kph=meqrec%nstph,jph,-1
          meqrec%stphl(kph+1)=meqrec%stphl(kph)
       enddo
! phase added at jph, (note jph may be equal to nstph+1)
       meqrec%stphl(jph)=phfix
       meqrec%nstph=meqrec%nstph+1
       meqrec%phr(phfix)%itadd=meqrec%noofits
       meqrec%phr(phfix)%curd%dgm=zero
       meqrec%phr(phfix)%curd%amfu=nodefixpham
       meqrec%phr(phfix)%stable=1
! set that the phase has fixed amount
       meqrec%phr(phfix)%phasestatus=PHFIXED
    else
! we are removing a phase, abs(phfix) already in meqrec%phr
!       meqrec%stphl(jph)=phfix
!       meqrec%nstph=meqrec%nstph+1
!       write(*,*)'Removing a phase: ',phfix
       if(phfix.ge.0) then
          gx%bmperr=4234
          goto 1000
       endif
       meqrec%phr(-phfix)%itadd=meqrec%noofits
       meqrec%phr(-phfix)%curd%dgm=zero
       meqrec%phr(-phfix)%curd%amfu=nodefixpham
       meqrec%phr(-phfix)%stable=1
! set that the phase has fixed amount
       meqrec%phr(-phfix)%phasestatus=PHFIXED       
    endif
!--------------
! mark that the phase is fix, we have to be careful not to exceed size
! Sigh, the fixed phases must be in sequential order ??? ... not done here
! ... maybe not needed ??
!    write(*,*)'added fix phase: ',phfix
    meqrec%nfixph=meqrec%nfixph+1
    if(meqrec%nfixph.gt.size(meqrec%fixpham)) then
       write(*,*)'Too many phases set fixed during mapping',&
            meqrec%nfixph,size(meqrec%fixpham)
       gx%bmperr=4235; goto 1000
    endif
! meqrec%nfixph is used to reduce the number of variables in the system
! matrix.  Fix phases have no variable amount.
    meqrec%fixph(1,meqrec%nfixph)=meqrec%phr(abs(phfix))%iph
    meqrec%fixph(2,meqrec%nfixph)=meqrec%phr(abs(phfix))%ics
    meqrec%fixpham(meqrec%nfixph)=zero
!    write(*,*)'Set fixed phase: ',meqrec%nfixph,&
!         meqrec%phr(abs(phfix))%iph,meqrec%phr(abs(phfix))%ics,PHFIXED
! I am not sure what this make but error allocating svar inside meq_sameset
!       meqrec%nv=meqrec%nv+1
!---------------------------------------------------
! call meq_sameset with new set of phases and axis condition removed
! If there is a phase change (iadd or irem nonzeri) or error it exits 
    sameadd=0
200 continue
    iadd=0; irem=0
!    write(*,'(a,3i5)')'In map_calcnode calling sameset for new node: ',&
!         meqrec%nstph,phfix
!
!    write(*,*)'SMP2A Calling meq_sameset from map_calcnode'
    call meq_sameset(irem,iadd,mapx,meqrec,meqrec%phr,inmap,ceq)
!
!   write(*,202)'Calculated node with fix phase: ',gx%bmperr,irem,iadd,ceq%tpval
202 format(a,3i4,2(1pe12.4))
201 continue
!-------------------------------------------------
! trouble if error or another phase wants to be stable/dissapear
! We may have to calculate with the axis fix again, maybe even read up
! the previous calculated equilibrium
    if(gx%bmperr.ne.0) then
       if(ocv()) write(*,*)'Error trying to calculate a node point',gx%bmperr
       if(gx%bmperr.eq.4187) goto 1000
    elseif(irem.gt.0) then
       if(ocv()) write(*,*)'Another phase wants to disappear',irem
       gx%bmperr=4222
    elseif(iadd.gt.0) then
! another phase wants to be stable
!      write(*,*)'SMPNODE: Another phase wants to be stable ',iadd,sameadd,phfix
       if(same_composition(iadd,meqrec%phr,meqrec,ceq,zero)) then
          iadd=0; goto 201
       endif
!       write(*,'(a,3i5,F10.2)')'Error: new phase stable: ',&
!            iremsave,iaddsave,iadd,ceq%tpval(1)
       gx%bmperr=4223
    else
! It worked to calculate with a new fix phase releasing all axis condition!!!
! *************************************************************
! check that the node point is global using grid minimizer
! ceq is copied inside global_equil_check and not destroyed??.
! mode=0 means do not recalculate if gridpoint below is found
       mode=0
!       write(*,*)'NOT Calling global check'
!       global=.TRUE.
!       write(*,*)'Check if nodepoint global'
! make a copy of the whole equilibrium record and set a pointer to the copy
! Does this really make a copy of the conditions etc inside ceq?
!       tceq=ceq
!       pceq=>tceq
!       write(*,*)'SMP value of T: ',pceq%tpval(1)
! SEGMENTATION FAULT and other strange errors after this call
! very difficult to find ... puhhhhh
! --- BUT THERE is still a segmentation fault
!       global=global_equil_check1(mode,addtupleindex,yfra,pceq)
       global=global_equil_check1(mode,addtupleindex,yfra,ceq)
       if(.not.global) then
          write(*,*)'gridminimizer found node point not global'
! set this line as INACTIVE and do not generate any start points
          mapline%status=ibset(mapline%status,EXCLUDEDLINE)
          gx%bmperr=4353
          goto 1000
       endif
! *************************************************************
       goto 500
    endif
! Problems, the simplest is to go back and try a smaller step
! But we must first remove the fix phase and restore the axis condition
!    write(*,54)'Error calculating node point? ',gx%bmperr,mapline%lasterr,&
!         irem,iadd,phfix,pcond%statev,mapline%problems,axval
54  format(a,2i5,5i3,1pe12.4)
!    if(maptop%tieline_inplane.gt.0) then
! if <0 isopleth, 0 step, >0 tie-lines in plane
!       write(*,*)'Tie-lines in plane:'
! if T axis maybe change to extensive axis ...
!    endif
    if(ocv()) write(*,*)'Error calculating node point, take shorter step'
    pcond%active=0
    pcond%prescribed=axval
    if(pcond%statev.eq.1) then
       meqrec%tpindep(1)=.FALSE.
       if(ocv()) write(*,55)'Marking that T is a condition again',&
            axval,ceq%tpval(1)
55     format(a,6(1pe14.6))
    elseif(pcond%statev.eq.2) then
       meqrec%tpindep(2)=.FALSE.
!       ceq%tpval(2)=value
    endif
!    write(*,*)'error in map_calcnode, remove phfix: ',phfix
    if(phfix.gt.0) then
! we must remove phfix from the list of stable phases and shift down
       meqrec%nstph=meqrec%nstph-1
       do iph=1,meqrec%nstph
          jj=meqrec%stphl(iph)
          if(jj.ge.phfix) then
             meqrec%stphl(iph)=meqrec%stphl(iph+1)
          endif
       enddo
! we must zero the last stable phase !!
       meqrec%stphl(meqrec%nstph+1)=0
       meqrec%phr(phfix)%itrem=meqrec%noofits
       meqrec%phr(phfix)%prevam=zero
       meqrec%phr(phfix)%stable=0
       meqrec%phr(phfix)%curd%amfu=zero
! we do not need to do anyting if -phfix should have been removed, then it
! is should remain among the stable phases, just remove it as fixed
    endif
    meqrec%fixph(1,meqrec%nfixph)=meqrec%phr(abs(phfix))%iph
    meqrec%fixph(2,meqrec%nfixph)=meqrec%phr(abs(phfix))%ics
    meqrec%nfixph=meqrec%nfixph-1
    mapline%lasterr=gx%bmperr
!    write(*,*)'SMP lasterr: ',mapline%lasterr,&
!         gx%bmperr,phfix,meqrec%phr(phfix)%phasestatus
! I had forgotten this!!
    meqrec%phr(abs(phfix))%phasestatus=0
! exit as no node found
    goto 1000
!------------------------------------------------------
! When we are there we have successfully calculated an equilibrium with a
! new phase set create a node with this equilibrium and a new line records
500 continue
!    write(*,*)'SMP2 Successful calculation of a node point',phfix
! phfix is set negative if phase should be removed
! NOTE the phase set fix in the node may not be the same which
! wanted to disappear/appear when calling the map_calcnode!!
! If iremsave=phfix the fix phase is one to be removed.
!    write(*,*)'SM2A node with new fix phase: ',phfix,iremsave,iaddsave
! I do not understand the next IF statement/BoS 200222
    if(iremsave.eq.-phfix) then
!       write(*,*)'In SMP2A with strange assignment ...',iremsave,-phfix
       phfix=-abs(phfix)
    endif
! if the user wants to have global minimization during mapping this is
! time to test if the current equilibrium is the global one.  We can use
! a temporary ceq record and chech the set of phases and chemical potentials
!
! NOTE that after a global equilibrium new composition set can have been
! created ... that should not be allowed unless they are really stable ...
! and one may have the same phases but different composition sets ... it
! can be quite messy.
! We have to set back the axis condition, before or after creating the node?
! and the new value ...
    if(pcond%noofterms.gt.1) then
       write(*,*)'Cannot handle conditions with several terms'
       gx%bmperr=4236; goto 1000
    endif
! this sets the condition as active
    pcond%active=0
    svrrec=>pcond%statvar(1)
    call state_variable_val(svrrec,value,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(ocv()) write(*,510)'Checking condition value; ',lastcond%seqz,&
         value,pcond%prescribed,ceq%tpval(1)
510 format(a,i3,6(1pe12.4))
! set the new condition value on the axis
    pcond%prescribed=value
    if(pcond%statev.eq.1) then
       meqrec%tpindep(1)=.FALSE.
       ceq%tpval(1)=value
       if(ocv()) write(*,*)'Marking that T is a condition again',value
    elseif(pcond%statev.eq.2) then
       meqrec%tpindep(2)=.FALSE.
       ceq%tpval(2)=value
    endif
! Save this as the last equilibrium of the line
    if(maptop%tieline_inplane.gt.0) then
! remove phfix as fix, otherwise graphics will be strange!
517    format(a,2i3,5x,5(2i3))
! remove phfix as fix
       if(phfix.lt.0) then
          write(*,*)'SM2A negative phfix used as index?',phfix
       endif
       mapline%meqrec%phr(phfix)%curd%phstate=PHENTERED
! this is necessary not to have data from this phase interfering with the line
       if(ocv()) write(*,519)phfix,mapline%meqrec%phr(phfix)%iph,&
            mapline%meqrec%phr(phfix)%ics,phentunst
519    format('Removing ',i3,2x,2i3,' as stable as last line equil',i3)
!?????????????????????????????????????????
       mapline%meqrec%nstph=mapline%meqrec%nstph-1
    endif
!    write(*,*)'Storing last point on line',phfix,maptop%tieline_inplane
    call map_store(mapline,axarr,maptop%number_ofaxis,maptop%saveceq)
    if(gx%bmperr.ne.0) then
!       if(gx%bmperr.eq.4300) write(*,*)'Node point ignored'
       goto 1000
    endif
! If we have an error here it may be that the node axis has big jumps
! Do not save any node
! here we have stored the last equilibrium that lead to th enode
! now update all condition records related to axis
!--------------------
! now store all axis values as prescribed vaules in the condition records
! A rather clumsy way and cannot handle expressions ...
    lastcond=>ceq%lastcondition
    pcond=>lastcond
600 continue    
       pcond=>pcond%next
       do jax=1,maptop%number_ofaxis
          if(pcond%seqz.eq.axarr(jax)%seqz) then
!             write(*,*)'At node set axis ',jax,axarr(jax)%lastaxval
             pcond%prescribed=axarr(jax)%lastaxval
          endif
       enddo
       if(.not.associated(pcond,lastcond)) goto 600
!-------------------
    if(maptop%tieline_inplane.gt.0) then
! Now set phfix back again for storing at the node record!!
       iph=1
       do jj=mapline%meqrec%nstph,1,-1
          if(mapline%meqrec%stphl(jj).gt.phfix) then
             mapline%meqrec%stphl(jj+1)=mapline%meqrec%stphl(jj)
          else
             iph=jj+1; exit
          endif
       enddo
       mapline%meqrec%stphl(iph)=phfix
       mapline%meqrec%phr(meqrec%stphl(iph))%curd%phstate=PHENTSTAB
       mapline%meqrec%nstph=mapline%meqrec%nstph+1
       if(ocv())write(*,517)'In map_calcnode: ',phfix,meqrec%nstph,&
           (meqrec%phr(meqrec%stphl(jj))%iph,meqrec%phr(meqrec%stphl(jj))%ics,&
           jj=1,meqrec%nstph)
!           meqrec%phr(meqrec%stphl(jj))%phstate,jj=1,meqrec%nstph)
518    format(a,2i3,5x,5(2i3,i2,2x))
    endif
!--------------------------------------------------------
!
    if(mapline%evenvalue.ne.zero) then
! if we have taken halfsteps then use the original even step
       if(ocv()) write(*,*)'Using original even step: ',mapline%evenvalue
       axval=mapline%evenvalue
    endif
!
! Finally create the new node and with new exit lines
    haha=0
    if(maptop%tieline_inplane.lt.0) then
! test if invariant ...
       if(inveq(haha,ceq)) then
! haha is set to number of stable phases at invariant.
! the number of lines ending at an invariant isopleth is 2*haha
! current number of stable phases is meqrec%nstph. 
! sign(1,phfix) is 1 if phfix>0; -1 if phfix<0
!          write(*,21)meqrec%nstph,haha,phfix,meqrec%nstph-haha+sign(1,phfix)
21        format('SMP2A stable phases mm: ',3i5,i10)
       endif
    endif
!    write(*,*)'SMP2A Test for invariant equilibrium: ',haha
    call get_phase_name(meqrec%phr(abs(phfix))%iph,meqrec%phr(abs(phfix))%ics,&
         phname)
    if(gx%bmperr.ne.0) then
       write(*,*)'SMP2A illegal phase name: ',phfix
       goto 1000
    endif
    if(phfix.gt.0) then
       write(*,501)ceq%tpval(1),trim(phname)
501    format('Creating a node at ',F10.2,' where ',a,' appears')
    else
       write(*,502)ceq%tpval(1),trim(phname)
502    format('Creating a node at ',F10.2,' where ',a,' disappear')
    endif
!    write(*,*)'calling map_newnode: ',mapline%meqrec%nfixph,meqrec%nfixph,haha
!    if(haha.gt.1) &
!         write(*,*)'SMP2A invariant!! we should greate several exits ',haha
! inside map_newnode the approriate number of exits will be generated
    call map_newnode(mapline,meqrec,maptop,axval,jax,axarr,phfix,haha,ceq)
    if(gx%bmperr.ne.0) then
       if(ocv()) write(*,*)'Error return from map_newnode: ',gx%bmperr
    endif
!    write(*,*)'Back from map_newnode',phfix
! all done??
1000 continue
    return
  end subroutine map_calcnode

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_newnode
!\begin{verbatim}
  subroutine map_newnode(mapline,meqrec1,maptop,axval,lastax,axarr,&
       phfix,haha,ceq)
! must be partially THREADPROTECTED
! first check if a node with this equilibrium already exists
! if not add a new node with appropriate lineheads and arrange all links
! Take care if tie-lines in the plane all lines do not have to be calculated
! NOTE: meqrec1 not the same as mapline%meqrec !! ??
! mapline is line record for current line
! meqrec1 has information about last calculated equilibrium
! maptop is node record
! axval is the axis value attemped to calculate when phase set wanted to change
! lastax is index of last active axis
! axarr are axis records
! phfix is phase which is set fix at node point
! haha is larger than 1 if the calculated equilibrium is invariant
! ceq is equilibrium record
    implicit none
    type(map_node), pointer :: maptop
    type(meq_setup) :: meqrec1
    type(map_line), pointer :: mapline,nodexit
    type(map_axis), dimension(*) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
    integer phfix,lastax,haha
    double precision axval
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: newceq,tmpceq
    type(map_node), pointer :: mapnode,newnode
    type(map_line), pointer :: linenode,tmpline
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(map_fixph), allocatable :: mapfix
    type(meq_setup), pointer :: meqrec2
    type(gtp_state_variable), target :: svrtarget
    integer remph,addph,nel,iph,ics,jj,seqx,nrel,jphr,stabph,kph,kcs,kk,lfix
    integer zph,stepax,kpos,seqy,jp,nopotax,lokcs,lokph,haha2,linefphr
! there should be 8 significant digits, first step factor
!    double precision, parameter :: vz=1.0D-9,axinc1=1.0D-3
    double precision, parameter :: vz=1.0D-8,axinc1=1.0D-3
    character eqname*24,phases*60
    double precision stepaxval,middle,testv,xxx
! mark that line ended with two stoichometric phases and data for isopleth inv
    integer twostoichset,errall,jfix,jstab,jlast,kstab,zp,phrix,infix3,nexit
    integer onlyone,notone,jused,zz,tz,qy,savefix1,savefix2,nodein(2),nodeut(2)
! lifexix, nodefix and prevfix are used to fix pair of phases that have zero
! amount at the exit points of lines from an invariant equilibrium.
    integer linefix,nodefix,infix,infix2,doubline,twice,firstoutfix,outfix,qq
    integer, allocatable, dimension(:,:) :: invph,nodeout
!    double precision, allocatable, dimension(:) :: exitcomp,eqcopy
    logical stepinvariantnode
!
    lfix=0
! the phase kept fix with zero amount at the node is phfix  It can be
! negative at STEP if it is a phase that will dissapear.
    if(ocv()) write(*,87)'We are in map_newnode with a fix phase: ',&
         phfix,ceq%tpval(1)
87  format(a,i4,2x,1pe12.4)
!    write(*,*)'We have access to phr: ',meqrec1%phr(abs(phfix))%iph,&
!         meqrec1%phr(abs(phfix))%ics
! mapnode should be set to point at maptop
    twostoichset=0
    if(btest(mapline%status,TWOSTOICH)) then
       twostoichset=1
       write(*,'(a)')'SMP line ended with two stochiometric allotropes stable'
    endif
    mapnode=>maptop
    nrel=meqrec1%nrel
100 continue
!---------------------------------------------------------------------
! loop all mapnodes to check if any has the same chemical potentials
!---------------------------------------------------------------------
!       write(*,*)'Comparing with node: ',mapnode%seqx,nrel
!       write(*,105)'T diff: ',ceq%tpval(1),mapnode%tpval(1),&
!            abs(ceq%tpval(1)-mapnode%tpval(1)),abs(vz*mapnode%tpval(1))
!       write(*,105)'P diff: ',ceq%tpval(2),mapnode%tpval(2),&
!            abs(ceq%tpval(2)-mapnode%tpval(2)),abs(vz*mapnode%tpval(2))
       if(abs(ceq%tpval(1)-mapnode%tpval(1)).gt.abs(vz*mapnode%tpval(1)) .or.&
            abs(ceq%tpval(2)-mapnode%tpval(2)).gt.abs(vz*mapnode%tpval(2))) then
!          write(*,*)'Not same, compare with next'
          goto 110
       endif
       do nel=1,nrel
!          write(*,105)'Chempot: ',ceq%cmuval(nel),mapnode%chempots(nel),&
!               abs(ceq%cmuval(nel)-mapnode%chempots(nel)),&
!               abs(vz*mapnode%chempots(nel))
105       format(a,5(1pe16.8))
          if(abs(ceq%cmuval(nel)-mapnode%chempots(nel)).gt.&
               abs(2.0D1*vz*mapnode%chempots(nel))) then
!             write(*,'(a,3(1pe12.4))')'SMP not same chempots, at this node',&
!                  abs(ceq%cmuval(nel)-mapnode%chempots(nel)),&
!                  abs(2.0D1*vz*mapnode%chempots(nel))
             goto 110
          endif
       enddo
! We can come here with a STEP command without any fix phases
       if(maptop%tieline_inplane.eq.0) then
          write(*,*)'SMP2A map_calcnode: Step command'
          goto 800
       endif
! T, P and all chemical potentials the same, one should maybe check phases??
       iph=mapline%linefixph(1)%ixphase
       ics=mapline%linefixph(1)%compset
!       if(ocv()) write(*,107)'Node exist: ',&
!       write(*,107)'Node already exist: ',&
!            mapnode%seqx,size(mapnode%linehead),iph,ics
107    format(a,i5,i3,i5,i2)
! do not remove exits from invariant nodes ...
       if(btest(mapnode%status,MAPINVARIANT)) goto 800
       removexit: do jj=1,size(mapnode%linehead)
! loop for all exits
          nodexit=>mapnode%linehead(jj)
          if(ocv()) write(*,108)'Exit: ',jj,nodexit%done,&
               nodexit%linefixph(1)%ixphase,nodexit%linefixph(1)%compset
!               nodexit%linefixph(1)%phaseix,nodexit%linefixph(1)%compset
108       format(a,i4,i7,i5,i2)
          if(nodexit%done.le.0) cycle
          if(nodexit%linefixph(1)%ixphase.eq.iph .and. &
               nodexit%linefixph(1)%compset.eq.ics) then
!             write(*,*)'Number of stable phases: ',&
!                  nodexit%nstabph,mapline%nstabph
             if(nodexit%nstabph.eq.mapline%nstabph) then
! if we have same number of stable phases they must be checked (invariant)
!                write(*,*)'Can be an invariant equilibrium!',mapline%nstabph
             endif
             mapnode%linehead(jj)%done=-1
             write(*,106)mapnode%linehead(jj)%lineid,jj,mapnode%seqx
106          format('Removed line ',i2,', exit ',i3,' from node: ',i3)
             exit removexit
          endif
       enddo removexit
       goto 800
! take next mapnode
110    continue
! difficult error to detect, I had written mapnode=mapnde%next !!!
       mapnode=>mapnode%next
! the next links should form a circular list ...
       if(.not.associated(mapnode,maptop)) goto 100
!==================================================================
! 
120 continue
    mapnode=>maptop%next
    seqx=mapnode%seqx+1
! if maptop%next is maptop do not nullify this pointer !!
! Always add the new record as the next link to maptop
    if(associated(mapnode,maptop)) then
! a single maptop record
!       write(*,*)'allocate mapnone%next 1'
       allocate(mapnode%next)
       mapnode%next%status=0
    else
! there is more mapnode records ... allocation here means memory leak
! I do not know how to fix ... it seems one can deallocate pointers!! no leak
!       write(*,*)'allocate mapnone%next 2'
       allocate(maptop%next)
       maptop%next%status=0
    endif
    newnode=>maptop%next
    newnode%first=>maptop
    newnode%next=>mapnode
    newnode%previous=>maptop
    mapnode%previous=>newnode
    newnode%seqx=seqx
!    write(*,*)'Maptop and next: ',maptop%seqx,maptop%next%seqx,newnode%seqx
!
    eqname='_MAPNODE_'
    jj=10
!    write(*,*)'SMP2A map_newnode copy equilibrium: ',seqx,nrel
    call wriint(eqname,jj,seqx)
! This copy is a record in the array "eqlista" of equilibrium record, thus
! it will be updated if new composition sets are created in other threads.
!    write(*,*)'Check 1: ',mapline%meqrec%nfixph,meqrec%nfixph,mapline%lineid,&
!         mapnode%seqx
    call copy_equilibrium(newceq,eqname,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error creating equilibrium: ',eqname
       goto 1000
    endif
    newnode%nodeceq=>newceq
! if twostoichset is set then add a comment in the new equilibrium
    newnode%artxe=0
    if(twostoichset.eq.1) then
!       write(*,*)'SMP2A setting artxe'
       newnode%artxe=1
    endif
! save a copy of ceq also in result (reserve is threadprotected)
    if(ocv()) write(*,*)'Copies node ceq to saveceq'
    call reserve_saveceq(jj,maptop%saveceq)
    if(gx%bmperr.ne.0) goto 1000
    maptop%saveceq%savedceq(jj)=newceq
    newnode%savednodeceq=jj
!    write(*,*)'Copy successful'
!    write(*,*)'Before copying meqrec: ',mapline%meqrec%nfixph,meqrec%nfixph
! maybe it is not necessary to save meqrec and chemical potentials??
    newnode%meqrec=meqrec1
!    write(*,*)'New node index: ',newnode%seqx
    allocate(newnode%chempots(nrel))
    newnode%chempots=ceq%cmuval
    newnode%tpval=ceq%tpval
!    newnode%type_of_node=0
! correct value of lines will be set later
    newnode%lines=0
    newnode%tieline_inplane=maptop%tieline_inplane
! this seems to be wrong, maptop%number_ofaxis is zero when step separate
    newnode%number_ofaxis=maptop%number_ofaxis
! save index of the phase set fix at the node
!    write(*,*)'SMP Saving index of new fix phase: ',abs(phfix)
    if(phfix.lt.0) then
       newnode%nodefix%ixphase=-meqrec1%phr(abs(phfix))%iph
    else
       newnode%nodefix%ixphase=meqrec1%phr(abs(phfix))%iph
    endif
    newnode%nodefix%compset=meqrec1%phr(abs(phfix))%ics
!    write(*,*)'Saved node fix phase: ',newnode%nodefix%phase,&
!         newnode%nodefix%compset
! the set of stable phases
    newnode%noofstph=meqrec1%nstph
    allocate(newnode%stable_phases(newnode%noofstph))
    do jj=1,newnode%noofstph
!       newnode%stable_phases(jj)%phaseix=meqrec1%iphl(jj)
       newnode%stable_phases(jj)%ixphase=meqrec1%iphl(jj)
       newnode%stable_phases(jj)%compset=meqrec1%icsl(jj)
    enddo
!
! >>>>>>>>>>>>>>>>>>> add code here to generate 2*haha-1 exuts    
!    if(haha.gt.0) write(*,*)'SMP2A found invariant with exits: !!',2*haha-1
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
! Thats all in the newnode ... except the lineheads ....
! Hm, when taking the different exits we must know phase sets and axis
! directions, with some efforts one could check which axis variable will
! change most rapidly for each exit but that can wait.  But I must know
! which phase set to have stable in the different lines ... but not for step.
! For invariant equilibra with tie-lines not in the plane that can be quite
! messy if I remeber correctly
!-----------------------------------------------------
! create mapline records in newnode with different sets of stable phases
    if(ocv()) write(*,*)'now generate lineheads',maptop%tieline_inplane,&
         mapline%meqrec%nfixph,meqrec1%nfixph
    if(maptop%tieline_inplane.eq.0) then
! this is a step, just one line and one exit with the new set of stable phases
       newnode%lines=1
       if(noel().eq.1) then
! step with single phase: problems with phase change as old phase still stable
          call get_phase_compset(phfix,1,lokph,lokcs)
! this change will reomve the previously stable phase in newnode and 
! below also in meqrec1
          jj=newnode%stable_phases(1)%ixphase
          newnode%stable_phases(1)=phasetuple(phfix)
          phfix=-jj
          newnode%nodeceq%phase_varres(lokcs)%PHSTATE=PHENTSTAB
          newnode%nodeceq%phase_varres(lokcs)%amfu=one
          newnode%noofstph=1
!          write(*,*)'SMP phases 1A: ',phfix,newnode%stable_phases(1)%ixphase,&
!               newnode%stable_phases(2)%ixphase
! But I had to remove the previously stable phase also this way !!
          call get_phase_compset(-phfix,1,lokph,lokcs)
          newnode%nodeceq%phase_varres(lokcs)%PHSTATE=PHENTUNST
          newnode%nodeceq%phase_varres(lokcs)%amfu=zero
       endif
    elseif(maptop%tieline_inplane.gt.0) then
! mapping with tie-lines in plane. Always 3 lines meet ... 2 new exits ??
! the number of exits depends on number of axis,
! for 2 axis 3 lines meet, for 3 axis (one of which is a potential) 4 lines
       newnode%lines=2
    elseif(haha.gt.0) then
! for mapping without tie-lines in plane and haha is nonzero then we are at
! an invariant equlibrium with haha stable phases and 2*haha-1 exiting lines
!       newnode%lines=2*jj-1
       if(inveq(haha2,ceq)) then
!          newnode%lines=2*haha-1
! Only few of the exit lines will be in the plane f the diagram.  Assume 8
! i.e. there will be 7 exits          
          newnode%lines=7
       else
          newnode%lines=3
       endif
    else
! mapping without tie-lines in plane
! at other node points 4 lines meets, 3 exits
       write(*,*)'Unknown type of node create exit lines: ',newnode%lines
       newnode%lines=3
    endif
! set link to end node in mapline
    mapline%end=>newnode
!=============================================================================
! we must create sufficient linehead records and initiate their content
! differently depening on STEP (case 1), MAP with tie-lines in plane (case 2)
! and MAP without tie-lines in plane (case 3).  In the latter case special
! care must be taken for invariant nodes. (for case 2 all nodes are invariants)
! check if we have a potential axis and select that as axandir
    stepax=mapline%axandir
    nopotax=0
    if(maptop%number_ofaxis.gt.1) then
!       write(*,*)'Seach for step axis'
       kk=0
       do jj=1,maptop%number_ofaxis
          if(axarr(jj)%axcond(1)%statevarid.lt.5) then
! positive or negative direction is unknown
             stepax=jj
             nopotax=jj
! the value of this condition is hopefully in the axarr(jj)%lastaxval ??
! It was stored there after calculating the node
!             write(*,*)'Found axis and value: ',axarr(jj)%lastaxval
             stepaxval=axarr(jj)%lastaxval
          endif
! save the axis with the value closest to the "middle" of the axis
          if(kk.eq.0) then
             kk=1
             middle=abs(5.0D-1-axarr(jj)%lastaxval/&
                  (axarr(jj)%axmax-axarr(jj)%axmin))
!             write(*,*)'middle: ',kk,middle
          else
             testv=abs(5.0D-1-axarr(jj)%lastaxval/&
                  (axarr(jj)%axmax-axarr(jj)%axmin))
             if(testv.lt.middle) then
                middle=testv
                kk=jj
             endif
!             write(*,*)'middle: ',kk,middle,testv
          endif
       enddo
       if(nopotax.eq.0) then
          stepax=kk
          stepaxval=axarr(kk)%lastaxval
       endif
!       write(*,*)'Set step axis to: ',stepax,&
!            axarr(stepax)%axcond(1)%statevarid
    endif
!
!
!    write(*,*)'SMP2A creating lineheads: ',haha,newnode%lines
!    if(newnode%lines.gt.3) write(*,*)'SMP: generate exit lines: ',newnode%lines
    allocate(newnode%linehead(newnode%lines),stat=errall)
    if(errall.ne.0) then
       write(*,*)'SMP2A Allocation error 1: ',errall
       gx%bmperr=4370; goto 1000
    endif
!    newnode%type_of_node=0
!
    do jp=1,newnode%lines
!--------------------- code moved from map_findline
! COPY of the equilibrium record from newnode to newnode%linehead(jp)%lineceq
       if(ocv()) write(*,*)'We found a line from node: ',mapnode%seqx
       newnode%linehead(jp)%meqrec%status=0
       eqname='_MAPLINE_'
       kpos=10
       seqy=maptop%seqy+1
       call wriint(eqname,kpos,seqy)
       call copy_equilibrium(newnode%linehead(jp)%lineceq,eqname,&
            newnode%nodeceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error creating equilibrium: ',eqname
          goto 1000
       endif
!       write(*,*)'SMP phases 2: ',seqy,phfix,newnode%stable_phases(1)%ixphase,&
!               newnode%nodeceq%phase_varres(lokcs)%phstate,newnode%noofstph
       maptop%seqy=seqy
       newnode%linehead(jp)%lineid=seqy
       newnode%linehead(jp)%nodfixph=0
! mapline%more is positive for line to be calculated, 0 means end at axis limit
       newnode%linehead(jp)%more=1
    enddo
!------------------------------ end code copied
!
!   write(*,*)'*** Trying to create node with # exit lines: ',haha,newnode%lines
! STEP has just 1 exit; 
! MAP tie-line in plane 2; isopleth non-invariant 3; isopleth invariant >3
    kpos=min(newnode%lines,4)
!    select case(newnode%lines)
    exits: select case(kpos)
!==========================================================================
    case default
       write(*,*)'SMP2A node error: exit lines= ',newnode%lines
       gx%bmperr=4237; goto 1000
!==========================================================================
    case(1)! step node with just one exit
! If phfix negative the fix phase wants to dissapear
       if(inveq(jj,ceq)) write(*,'(a)')'This is an invariant node'
       changephaseset: if(phfix.lt.0) then
! remove a phase ---------------------------
          remph=-phfix
!          write(kou,88)remph,' disappears,',meqrec1%nstph
88        format('SMP a node created where phase ',i3,a,' stable phases:',i3)
          if(meqrec1%nstph.eq.1) then
             write(*,*)'Attempt to remove the only stable phase!!!'
             gx%bmperr=4238; goto 1000
          endif
! shift phases after remph up?? in meqrec1%stphlnewnode%lines)
! irem is index to meqrec1%phr(), meqrec1%stphl(jph) is index to meqrec1%phr
          meqrec1%nstph=meqrec1%nstph-1
          do iph=1,meqrec1%nstph
             jj=meqrec1%stphl(iph)
             if(jj.ge.remph) then
                meqrec1%stphl(iph)=meqrec1%stphl(iph+1)
             endif
          enddo
! we must zero the last phase, hm itrem is not really relevant ...
          meqrec1%stphl(meqrec1%nstph+1)=0
! occational error because "remph" has illegal index value for meqrec1%phr 
          if(remph.le.0) then
             write(*,*)'Occational error around line 4487',remph
             remph=-remph
          endif
          meqrec1%phr(remph)%itrem=meqrec1%noofits
          meqrec1%phr(remph)%prevam=zero
          meqrec1%phr(remph)%stable=0
          meqrec1%phr(remph)%curd%amfu=zero
!          write(*,*)'SMP lineeq 3: ',meqrec1%nstph,meqrec1%stphl(1)
       elseif(phfix.gt.0) then
! we have to add phase phfix to the stable phase set, that is no problem
! as it is already in all lists, just remove that it should be fix
!          write(kou,88)phfix,' appears,   ',meqrec1%nstph
! meqrec1%nfixph seems not to be used .... ??
          if(meqrec1%nfixph.gt.0) then
             meqrec1%fixph(1,meqrec1%nfixph)=0
             meqrec1%fixph(2,meqrec1%nfixph)=0
             meqrec1%phr(phfix)%phasestatus=PHENTSTAB
             meqrec1%nfixph=meqrec1%nfixph-1
          endif
          stepinvariantnode=.FALSE.
          if(inveq(jj,ceq)) then
!
! if node is invariant we must remove one phase, which? NOT phfix
             stepinvariantnode=.TRUE.
             newnode%status=ibset(newnode%status,STEPINVARIANT)
!             write(*,*)'SMP2A invariant node at step',phfix,newnode%status,&
!                  meqrec1%nstph
! set the invariant bit in the node and calculate en equilibrium at
! a very small step above the invariant to find the new set of phases
             newnode%linehead(1)%meqrec=meqrec1
             tmpline=>newnode%linehead(1)
!             do kk=1,meqrec1%nstph
!                jj=meqrec1%stphl(kk)
!                write(*,294)'SMP initial set of phases: ',kk,jj,&
!                     meqrec1%phr(jj)%iph,meqrec1%phr(jj)%curd%amfu
!             enddo
             meqrec2=>tmpline%meqrec
!             do kk=1,meqrec2%nstph
!                jj=meqrec2%stphl(kk)
!                write(*,294)'SMP same initial set of phases: ',kk,jj,&
!                     meqrec2%phr(jj)%iph,meqrec2%phr(jj)%ics,&
!                     meqrec2%phr(jj)%curd%amfu
!             enddo
294          format(a,3i5,i2,1pe14.6)
             call locate_condition(axarr(1)%seqz,pcond,tmpline%lineceq)
             if(gx%bmperr.ne.0) goto 100
!             call list_conditions(kou,tmpline%lineceq)
!             call list_sorted_phases(kou,tmpline%lineceq)
!             if(gx%bmperr.ne.0) goto 100
             pcond%prescribed=pcond%prescribed+&
                  1.0D-3*stepax*axarr(1)%axinc
!             call list_conditions(kou,tmpline%lineceq)
!             write(*,*)'SMP small step invariant to find phase which disappear'
             call calceq3(0,.FALSE.,tmpline%lineceq)
! first argument -1 to keep the datastructure in meqrec2
!             call calceq7(-1,meqrec2,mapfix,tmpline%lineceq)
!             write(*,*)'Back from calceqx',gx%bmperr,meqrec1%nstph
!             call list_sorted_phases(kou,tmpline%lineceq)
! NOTE the content of meqrec2 has not been updated as calceq3 creates a new
! independent meqrec structure.  We must copy the values of phase amounts
! using a pointer directly to the phase_varres record
! list amount of phases after this small step.  However, the layout of
! the meqrec records are the same, we can use phase indices and other things
             do kk=1,meqrec2%nstph-1
                jj=meqrec2%stphl(kk)
                call get_phase_compset(meqrec2%phr(jj)%iph,meqrec2%phr(jj)%ics,&
                     lokph,lokcs)
                xxx=tmpline%lineceq%phase_varres(lokcs)%amfu
!                write(*,294)'SMP new set of phases: ',kk,jj,&
!                     meqrec2%phr(jj)%iph,meqrec2%phr(jj)%ics,xxx
                if(xxx.gt.zero) then
                   meqrec2%phr(jj)%curd%amfu=xxx
                else
                   do zz=kk,meqrec2%nstph-1
                      meqrec2%stphl(zz)=meqrec2%stphl(zz+1)
                   enddo
                endif
             enddo
             meqrec2%nstph=meqrec2%nstph-1
!             do kk=1,meqrec2%nstph
!                jj=meqrec2%stphl(kk)
!                write(*,294)'SMP final initial set of phases: ',kk,jj,&
!                     meqrec2%phr(jj)%iph,meqrec2%phr(jj)%ics,&
!                     meqrec2%phr(jj)%curd%amfu
!             enddo
! finally copy to meqrec1 ...
             meqrec1%nstph=meqrec2%nstph
             do kk=1,meqrec2%nstph
                meqrec1%stphl(kk)=meqrec2%stphl(kk)
! I assume the amounts are not needed, they should already be in lineceq ...??
             enddo
! rearrange the array of stable phases, one should be removed
!             stop 'all OK?'
          endif
       else
          write(*,*)'This is another never never error',phfix
          gx%bmperr=4239; goto 1000
       endif changephaseset
! set values in linhead record
       if(ocv()) write(*,*)'Creating linehead node record in: ',newnode%seqx
       newnode%linehead(1)%number_of_equilibria=0
       newnode%linehead(1)%first=0
       newnode%linehead(1)%last=0
!       newnode%linehead(1)%lineid=0
!       newnode%linehead(1)%axchange=1
       newnode%linehead(1)%axchange=-1
       newnode%linehead(1)%done=1
       newnode%linehead(1)%status=0
       newnode%linehead(1)%more=1
       newnode%linehead(1)%termerr=0
       newnode%linehead(1)%axfact=1.0D-2
       newnode%linehead(1)%nfixphases=0
! try to get a nice output of stable phases below
!       if(stepinvariantnode) then
!          allocate(newnode%linehead(1)%stableph(meqrec2%nstph))
!          allocate(newnode%linehead(1)%stable_phr(meqrec2%nstph))
!          newnode%linehead(1)%nstabph=0
!          do iph=1,meqrec2%nstph
!             newnode%linehead(1)%nstabph=newnode%linehead(1)%nstabph+1
!             jj=meqrec2%stphl(iph)
!             newnode%linehead(1)%stableph(iph)%ixphase=meqrec2%phr(jj)%iph
!             newnode%linehead(1)%stableph(iph)%compset=meqrec2%phr(jj)%ics
!             newnode%linehead(1)%stable_phr(iph)=jj
!          enddo
!       else
          allocate(newnode%linehead(1)%stableph(meqrec1%nstph))
          allocate(newnode%linehead(1)%stable_phr(meqrec1%nstph))
          newnode%linehead(1)%nstabph=0
          do iph=1,meqrec1%nstph
             newnode%linehead(1)%nstabph=newnode%linehead(1)%nstabph+1
             jj=meqrec1%stphl(iph)
             newnode%linehead(1)%stableph(iph)%ixphase=meqrec1%phr(jj)%iph
             newnode%linehead(1)%stableph(iph)%compset=meqrec1%phr(jj)%ics
             newnode%linehead(1)%stable_phr(iph)=jj
          enddo
!       endif
! end attempt
       newnode%linehead(1)%firstinc=1.0D-2*axarr(1)%axinc*mapline%axandir
!       newnode%linehead(1)%firstinc=1.0D-3*axarr(1)%axinc*mapline%axandir
       newnode%linehead(1)%evenvalue=axval
       newnode%linehead(1)%start=>newnode
       nullify(newnode%linehead(1)%end)
       if(ocv()) write(*,333)mapline%axandir,newnode%linehead(1)%firstinc,&
            newnode%linehead(1)%evenvalue
333    format('linehead: ',i3,2(1pe15.6))
       newnode%linehead(1)%axandir=mapline%axandir
!============================================================================
    case(2) ! Step node with two exits: Tie-lines in plane node, 3 lines meet,
!              2 new exits
!       write(*,*)'Trying to implement "tie-lines in plane" nodes'
       if(ocv()) write(*,*)'Creating linehead node record in: ',newnode%seqx
!       write(*,*)'Creating linehead node record in: ',newnode%seqx
! this is probably redundant, fixph already reset
       if(meqrec1%nfixph.gt.0) then
          meqrec1%fixph(1,meqrec1%nfixph)=0
          meqrec1%fixph(2,meqrec1%nfixph)=0
          meqrec1%phr(phfix)%phasestatus=PHENTSTAB
          meqrec1%nfixph=meqrec1%nfixph-1
       endif
!-------------- 
! no need for loop here I guess ... but I am oldfashioned
! begin doublecheck
       if(newnode%lines.ne.size(newnode%linehead)) then
          write(*,*)'SMP2A Trouble ahead!!'
          stop
       endif
! end doublecheck
       do jj=1,2
! initiate data in map_line record
          newnode%linehead(jj)%number_of_equilibria=0
          newnode%linehead(jj)%first=0
          newnode%linehead(jj)%last=0
!          newnode%linehead(jj)%lineid=0
!          newnode%linehead(jj)%axchange=1
          newnode%linehead(jj)%axchange=-1
          newnode%linehead(jj)%done=1
          newnode%linehead(jj)%status=0
          newnode%linehead(jj)%more=1
          newnode%linehead(jj)%termerr=0
          newnode%linehead(jj)%axfact=1.0D-2
!          newnode%linehead(jj)%axandir=mapline%axandir
! ??????????? can stepax be negative ???????????????
          newnode%linehead(jj)%axandir=stepax
          newnode%linehead(jj)%nfixphases=1
! this dimensioning is OK for two axis, if 3 it should be 2 etc.
          allocate(newnode%linehead(jj)%linefixph(1))
          allocate(newnode%linehead(jj)%linefix_phr(1))
! with tie-lines in the plane there is always just one stable phase
          allocate(newnode%linehead(jj)%stableph(1))
          allocate(newnode%linehead(jj)%stablepham(1))
          allocate(newnode%linehead(jj)%stable_phr(1))
! a small first step in same axis as used to find the node 
! We may have to change direction, in particular if the nodephase reappears
!          newnode%linehead(jj)%firstinc=1.0D-2*axinc*mapline%axandir
          newnode%linehead(jj)%firstinc=axinc1*axarr(abs(stepax))%axinc
          newnode%linehead(jj)%evenvalue=zero
! node records at start and end
          newnode%linehead(jj)%start=>newnode
          nullify(newnode%linehead(jj)%end)
       enddo
! This node represent a point where 3 lines meet 4 if 3 axis), each with a 
! different phase fix with zero amount.  One line is the one we followed
! to find the node, no need to generate an exit for that.
! It seems we do not have to bother so much with nfixph and fixph ...
! In meqrec1%phr there are currently two fixed phases, one which was fixed
! along the line (LFIX in mapline%linefixph), the other fixed for the node
! point, given by PHFIX which is an index to meqrec1%phr.  The third phase
! was stable with positive amount along the line LENT)
! The three lines are: FIX    STABLE    UNSTABLE
! already done         LFIX   LENT      PHFIX
! exit 1               PHFIX  LFIX      LENT
! exit 2               LENT   PHFIX     LFIX
       jphr=0
       if(allocated(mapline%linefixph)) then
          if(size(mapline%linefixph).gt.1) then
! If there are 3 axis this would be OK
             write(*,*)'SMP2B too many fix phases ...',size(mapline%linefixph)
             gx%bmperr=4240; goto 290
          endif
       endif
!       write(*,888)mapline%linefixph(1)%ixphase,mapline%linefixph(1)%compset,&
!            phfix,meqrec1%nphase,abs(phfix)
888    format('Old fix phase 2: ',i3,i2,', new fix phase: ',i3,&
            ', number of phases: ',i3,' abs(phfix): ',i3)
       do jj=1,mapline%meqrec%nphase
! loop through whole phr array to be sure nothing is wrong
          if(mapline%meqrec%phr(jj)%stable.eq.1) then
             if(jj.eq.abs(phfix) .or.&
                  (meqrec1%phr(jj)%iph.eq.mapline%linefixph(1)%ixphase .and.&
                  meqrec1%phr(jj)%ics.eq.mapline%linefixph(1)%compset)) cycle
             if(jphr.gt.0) then
                write(*,*)'Problems, two entered phases: ',jj,jphr
                gx%bmperr=4241; goto 290
             else
                jphr=jj
!                write(*,*)'Found entered phase: ',jphr
             endif
          endif
       enddo
! jphr is the phase that was stable along the line
       if(jphr.eq.0) then
          write(*,*)'Problems, not a single entered phase!'
          gx%bmperr=4242; goto 290
       endif
       zph=0
       do jj=1,meqrec1%nphase
          if(meqrec1%phr(jj)%iph.eq.mapline%linefixph(1)%ixphase .and. &
               meqrec1%phr(jj)%ics.eq.mapline%linefixph(1)%compset) then
! this is the index in phr for the phase that was fix along the line
             zph=jj
          endif
       enddo
       if(zph.eq.0) then
          write(*,203)' *** warning: cannot find the fix phase: ',zph,&
               mapline%linefixph(1)%ixphase,mapline%linefixph(1)%compset
203       format(a,10i4)
       endif
! For isothermal sections with no potential axis we must change the axis
! condition when following a new line
!       if(nopotax.eq.0) then
!          write(*,*)'Changing the axis variable for the new entered phase'
!       endif
! In mapnode there is a nfixph and array linefixph
!       if(ocv()) write(*,207)mapline%linefixph(1)%phaseix,&
!       if(ocv()) write(*,207)mapline%linefixph(1)%ixphase,&
!       write(*,207)mapline%linefixph(1)%ixphase,&
!            mapline%linefixph(1)%compset,&
!            meqrec%phr(phfix)%iph,meqrec%phr(phfix)%ics,&
!            meqrec%phr(jphr)%iph,meqrec%phr(jphr)%ics
207    format('LFIX: ',2i3,5x,' PHFIX: ',2i3,5x,' LENT: ',2i3)
! The two exits are:   FIX    STABLE    UNSTABLE
! exit 1               PHFIX  LFIX      LENT
! exit 2               LENT   PHFIX     LFIX
! for alcrni           
       iph=mapline%linefixph(1)%ixphase
       ics=mapline%linefixph(1)%compset
       linefphr=mapline%linefix_phr(1)
       phrix=mapline%linefix_phr(1)
       newnode%linehead(1)%linefixph%ixphase=meqrec1%phr(phfix)%iph
       newnode%linehead(1)%linefixph%compset=meqrec1%phr(phfix)%ics
       newnode%linehead(1)%linefix_phr=phfix
       newnode%linehead(1)%nstabph=1
! the previously fix phase is set as entered with stablepham as initial amount
       newnode%linehead(1)%stableph(1)%ixphase=iph
       newnode%linehead(1)%stableph(1)%compset=ics
! value of %stable_phr=??
       newnode%linehead(1)%stable_phr(1)=phrix
       newnode%linehead(1)%stablepham(1)=one
! store the phase number that must not become stable in nodfixph
       newnode%linehead(1)%nodfixph=jphr
!       newnode%linehead(1)%nodfixtup=meqrec1%phr(jphr)%phtupix
!       write(*,*)'SM2A nodfix:',phasetuple(meqrec1%phr(jphr)%phtupix)%ixphase,&
!            meqrec1%phr(jphr)%ics
!-----------
       newnode%linehead(2)%linefixph%ixphase=meqrec1%phr(jphr)%iph
       newnode%linehead(2)%linefixph%compset=meqrec1%phr(jphr)%ics
       newnode%linehead(2)%linefix_phr=jphr
       newnode%linehead(2)%nstabph=1
       newnode%linehead(2)%stableph(1)%ixphase=meqrec1%phr(phfix)%iph
       newnode%linehead(2)%stableph(1)%compset=meqrec1%phr(phfix)%ics
       newnode%linehead(2)%stable_phr(1)=phfix
       newnode%linehead(2)%stablepham(1)=one
! store the phase number that must not become stable in nodfixph
       newnode%linehead(2)%nodfixph=zph
!       newnode%linehead(1)%nodfixtup=meqrec1%phr(zph)%phtupix
       if(nopotax.eq.0) then
! If we have no potential axis we MUST change the axis condition
! to represent the axis composition of the new stable phase
          write(*,712)stepax,axarr(abs(stepax))%axcond(1)%statevarid,stepaxval
712       format('Creating nodepoint with no potential axis: ',2i4,1pe12.4)
! we have to change the axis condition to be the current composition of the
! new stable phase.  
!          write(*,*)'Conditions at node point'
!          call list_conditions(kou,newnode%nodeceq)
! change condition value for the lines exiting this node point
          tmpceq=>newnode%linehead(1)%lineceq
          call locate_condition(axarr(stepax)%seqz,pcond,tmpceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Cannot locate condition: ',axarr(stepax)%seqz
             goto 1000
          endif
          svrrec=>pcond%statvar(1)
!          call state_variable_val(svrrec,xxx,tmpceq)
!          if(gx%bmperr.ne.0) goto 1000
!          write(*,*)'Condition/State variable value: ',xxx
! NOTE: If we change fix/entered phase we must change axvals/axvals2
!              i1=svr2%argtyp; i2=svr2%phase; i3=svr2%compset
          svrtarget=svrrec
          svrtarget%argtyp=3
          svrtarget%phase=newnode%linehead(1)%stableph(1)%ixphase
          svrtarget%compset=newnode%linehead(1)%stableph(1)%compset
! This extracts the composition of the entered phase for first new line
! we must use a pointer in state_variable_val
          svr2=>svrtarget
          call state_variable_val(svr2,xxx,tmpceq)
          if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to extract the value, 0 means to set the value
          call condition_value(0,pcond,xxx,tmpceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error setting new conditions in lineceq'
             goto 1000
          endif
!          write(*,*)'Setting condition for line 1 to ',xxx
!          write(*,211)'New conditions at line: ',newnode%linehead(1)%lineid,&
!               trim(newnode%linehead(1)%lineceq%eqname),&
!               svr2%phase,svr2%compset,xxx
211       format(a,i3,a,2x,' phase/set: ',2i3,2x,1pe12.4)
!          call list_conditions(kou,newnode%linehead(1)%lineceq)
!--------- second exit
          tmpceq=>newnode%linehead(2)%lineceq
          call locate_condition(axarr(stepax)%seqz,pcond,tmpceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Cannot locate condition: ',axarr(stepax)%seqz
             goto 1000
          endif
          svrrec=>pcond%statvar(1)
          svrtarget=svrrec
          svrtarget%argtyp=3
          svrtarget%phase=newnode%linehead(2)%stableph(1)%ixphase
          svrtarget%compset=newnode%linehead(2)%stableph(1)%compset
! This extracts the composition of the entered phase for second new line
! ONLY CHANGE ... has no influence on the problem ...
          svr2=>svrtarget
! the line above should be there I think but was missing ... xxx below wrong??
          call state_variable_val(svr2,xxx,tmpceq)
          if(gx%bmperr.ne.0) goto 1000
          call condition_value(0,pcond,xxx,tmpceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error setting new conditions in lineceq'
             goto 1000
          endif
!          write(*,*)'Setting condition for line 2 to ',xxx
!          write(*,211)'New conditions at line: ',newnode%linehead(2)%lineid,&
!               trim(newnode%linehead(2)%lineceq%eqname),&
!               svr2%phase,svr2%compset,xxx
!          call list_conditions(kou,newnode%linehead(2)%lineceq)
       endif
! list the exits:
       if(ocv()) write(*,56)'Created linehead 1 for node: ',newnode%seqx,&
!       write(*,56)'Created linehead 1 for node: ',newnode%seqx,&
            newnode%linehead(1)%linefixph%ixphase,&
            newnode%linehead(1)%linefixph%compset,&
            newnode%linehead(1)%stableph(1)%ixphase,&
            newnode%linehead(1)%stableph(1)%compset,&
            newnode%linehead(1)%nodfixph
       if(ocv()) write(*,56)'Created linehead 2 for node: ',newnode%seqx,&
!       write(*,56)'Created linehead 2 for node: ',newnode%seqx,&
            newnode%linehead(2)%linefixph%ixphase,&
            newnode%linehead(2)%linefixph%compset,&
            newnode%linehead(2)%stableph(1)%ixphase,&
            newnode%linehead(2)%stableph(1)%compset,&
            newnode%linehead(2)%nodfixph
56     format(a,i3,5x,2i3,5x,2i3,5x,2i3)
       if(newnode%lines.ne.2) then
          write(*,*)'SMP2A setting newnode%lines'
          newnode%lines=2
       endif
! the fix and stable phases must be copied to meqrec1 when line is started
!       write(*,*)'Created node with 2 exits: ',newnode%seqx,ceq%tpval(1)
! prevent the lines from being used as that makes the program crash
290    continue
!======================================================================
    case(3) ! Normal node in a phase diagram without tie-lines in plane
! Two crossing lines, one in and 3 exits
! THERE IS NO CASE WHEN FINDING AN INVARIANT
! this is probably redundant, fixph already reset
       if(meqrec1%nfixph.gt.0) then
!          write(*,*)'Not redundant ...'
          meqrec1%fixph(1,meqrec1%nfixph)=0
          meqrec1%fixph(2,meqrec1%nfixph)=0
          meqrec1%phr(abs(phfix))%phasestatus=PHENTSTAB
          meqrec1%nfixph=meqrec1%nfixph-1
       endif
!-------------- 
! no need for loop here I guess ... but I am oldfashioned
       do jj=1,3
! initiate data in map_line record
          newnode%linehead(jj)%number_of_equilibria=0
          newnode%linehead(jj)%first=0
          newnode%linehead(jj)%last=0
!          newnode%linehead(jj)%lineid=0
!          newnode%linehead(jj)%axchange=1
          newnode%linehead(jj)%axchange=-1
          newnode%linehead(jj)%done=1
          newnode%linehead(jj)%status=0
          newnode%linehead(jj)%more=1
          newnode%linehead(jj)%termerr=0
          newnode%linehead(jj)%axfact=1.0D-2
!          newnode%linehead(jj)%axandir=mapline%axandir
          newnode%linehead(jj)%axandir=stepax
          newnode%linehead(jj)%nfixphases=1
! this dimensioning is OK for two axis, if 3 axis it should be 2 etc.
          allocate(newnode%linehead(jj)%linefixph(1))
          allocate(newnode%linehead(jj)%linefix_phr(1))
! There will be different number of stable phases in the lines 
! if new phase appear, 2 lines with jphr+1, one with jphr
! if old phase dissapear, 2 lines with jphr-1, one with jphr
          jphr=mapline%nstabph
          if(jphr.eq.1 .and. phfix.lt.0) then
             write(*,*)'Trying to remove the only entered phase !'
             gx%bmperr=4238; goto 1000
          endif
          if(jj.eq.3) then
!             write(*,*)'Allocating stableph: ',jj,jphr
             allocate(newnode%linehead(jj)%stableph(jphr))
             allocate(newnode%linehead(jj)%stablepham(jphr))
             allocate(newnode%linehead(jj)%stable_phr(jphr))
          else
             if(phfix.lt.0) then
!                write(*,*)'Allocating stableph: ',jj,jphr-1
                allocate(newnode%linehead(jj)%stableph(jphr-1))
                allocate(newnode%linehead(jj)%stablepham(jphr-1))
                allocate(newnode%linehead(jj)%stable_phr(jphr-1))
             else
!                write(*,*)'Allocating stableph: ',jj,jphr+1
                allocate(newnode%linehead(jj)%stableph(jphr+1))
                allocate(newnode%linehead(jj)%stablepham(jphr+1))
                allocate(newnode%linehead(jj)%stable_phr(jphr+1))
             endif
          endif
! a small first step in same axis as used to find the node 
! We may have to change direction, in particular if the nodephase reappears
! evenvalue important only for STEP with one axis
!          newnode%linehead(jj)%firstinc=1.0D-2*axinc*mapline%axandir
!          newnode%linehead(jj)%firstinc=axinc1*axinc*mapline%axandir
          newnode%linehead(jj)%firstinc=axinc1*axarr(abs(stepax))%axinc
          newnode%linehead(jj)%evenvalue=zero
! links to node records at start and end of line
          newnode%linehead(jj)%start=>newnode
          nullify(newnode%linehead(jj)%end)
       enddo
!
       if(allocated(mapline%linefixph)) then
          if(size(mapline%linefixph).gt.1) then
! error if 2 axis but would be OK if 3 axis
             write(*,*)'Problem, many fix phases!',size(mapline%linefixph)
             gx%bmperr=4240; goto 390
          endif
       endif
       stabph=0
       do jj=1,meqrec1%nphase
! loop through whole phr array to be sure nothing is wrong
          if(meqrec1%phr(jj)%stable.eq.1) then
! there should be 2 fixed phases, one along the line and one at the node
! If 3 or more axis there will be more fixed phases, phfix can be negative
!             if(jj.eq.abs(phfix) .or.&
!                  (meqrec1%phr(jj)%iph.eq.mapline%linefixph(1)%phase .and.&
!                   meqrec1%phr(jj)%ics.eq.mapline%linefixph(1)%compset)) cycle
! we should include phfix in stabph!!
!             if(meqrec1%phr(jj)%iph.eq.mapline%linefixph(1)%phaseix .and.&
             if(meqrec1%phr(jj)%iph.eq.mapline%linefixph(1)%ixphase .and.&
                  meqrec1%phr(jj)%ics.eq.mapline%linefixph(1)%compset) cycle
             stabph=stabph+1
          endif
       enddo
! Hm, stabph calculated this way is wrong, use mapline%nstabph
!       write(*,312)meqrec1%nphase,stabph,mapline%nstabph,phfix,&
!            mapline%linefixph(1)%phase,mapline%linefixph(1)%compset
312    format('In map_newnode 312: ',i3,i5,i3,i5,i3,i2,10i5)
       stabph=mapline%nstabph
       if(stabph.eq.0) then
          write(*,*)'Problems, no entered phase!'
          gx%bmperr=4242; goto 390
       endif
! 4 lines meet in all nodes except invariants
! We have 1 fix phase and f=stabph enterend phases and 1 new/old (+/-PHFIX)
! LFIX is phase fix along line up to node
! +PHFIX is a new phase becommong stable, -PHFIX is stable phase dissapearing
! if not invariant generate 3 exits, with PHFIX>0 these are
! The three exits are:   FIX    STABLE phases             nodfixph
! exit 1                 LFIX   f+1 (+PHFIX)              PHFIX
! exit 2                 PHFIX  f+1 (+LFIX)               LFIX
! exit 3                 PHFIX  f   (-PHFIX and add LFIX) LFIX
! If PHFIX<0 is an old phase becommong unstable
! The three exits are:   FIX    STABLE    not allwed appear/disappear
! exit 1                 LFIX   f-1 (-PHFIX)              PHFIX
! exit 2                 PHFIX  f-1 (-PHFIX not LFIX)     LFIX
! exit 3                 PHFIX  f   (-PHFIX and add LFIX) LFIX
       iph=mapline%linefixph(1)%ixphase
       ics=mapline%linefixph(1)%compset
       linefphr=mapline%linefix_phr(1)
! for use below I need to know the position of iph+ics in meqrec1%phr ...
       flfix: do jj=1,meqrec1%nstph
          if(meqrec1%phr(jj)%iph.eq.iph .and. meqrec1%phr(jj)%ics.eq.ics) then
             lfix=jj; exit flfix
          endif
       enddo flfix
       kph=mapline%meqrec%phr(abs(phfix))%iph
       kcs=mapline%meqrec%phr(abs(phfix))%ics
! exit 1 has same linefix as incomming line ------------------------
       newnode%linehead(1)%linefixph%ixphase=iph
       newnode%linehead(1)%linefixph%compset=ics
       newnode%linehead(1)%linefix_phr=linefphr
       if(phfix.gt.0) then
!          write(*,*)'allocated size of stableph 2: ',size(mapline%stableph)
          do jj=1,stabph
             newnode%linehead(1)%stableph(jj)%ixphase=&
                  mapline%stableph(jj)%ixphase
             newnode%linehead(1)%stableph(jj)%compset=&
                  mapline%stableph(jj)%compset
             newnode%linehead(1)%stablepham(jj)=mapline%stablepham(jj)
             newnode%linehead(1)%stable_phr(jj)=mapline%stable_phr(jj)
          enddo
! add phfix as stable phase
          jj=stabph+1
          newnode%linehead(1)%stableph(jj)%ixphase=kph
          newnode%linehead(1)%stableph(jj)%compset=kcs
          newnode%linehead(1)%stablepham(jj)=zero
          newnode%linehead(1)%stable_phr(jj)=abs(phfix)
! UNFINISHED check why stable_phr and nodefxph same??
          newnode%linehead(1)%nodfixph=abs(phfix)
!          newnode%linehead(1)%nodfixtup=meqrec1%phr(abs(phfix))%phtupix
          newnode%linehead(1)%nstabph=jj
       else
! phfix is negative, a phase disappear
          kk=0
          do jj=1,stabph-1
! remove -phfix as stable phase
             if(mapline%stableph(jj)%ixphase.eq.kph .and.&
                  mapline%stableph(jj)%compset.eq.kcs) then
                kk=jj+1
             else
                kk=kk+1
             endif
             newnode%linehead(1)%stableph(jj)%ixphase=&
                  mapline%stableph(kk)%ixphase
             newnode%linehead(1)%stableph(jj)%compset=&
                  mapline%stableph(kk)%compset
             newnode%linehead(1)%stablepham(jj)=mapline%stablepham(kk)
             newnode%linehead(1)%stable_phr(jj)=mapline%stable_phr(kk)
          enddo
          newnode%linehead(1)%nodfixph=abs(phfix)
          newnode%linehead(1)%nstabph=stabph-1
       endif
!
! exit 2 has PHFIX as linefix ----------------------------------
       newnode%linehead(2)%linefixph%ixphase=kph
       newnode%linehead(2)%linefixph%compset=kcs
       newnode%linehead(2)%linefix_phr=abs(phfix)
       if(phfix.gt.0) then
          do jj=1,stabph
             newnode%linehead(2)%stableph(jj)%ixphase=&
                  mapline%stableph(jj)%ixphase
             newnode%linehead(2)%stableph(jj)%compset=&
                  mapline%stableph(jj)%compset
             newnode%linehead(2)%stablepham(jj)=mapline%stablepham(jj)
             newnode%linehead(2)%stable_phr(jj)=mapline%stable_phr(jj)
          enddo
! add LFIX as stable phase
          jj=stabph+1
          newnode%linehead(2)%stableph(jj)%ixphase=iph
          newnode%linehead(2)%stableph(jj)%compset=ics
          newnode%linehead(2)%stablepham(jj)=zero
          newnode%linehead(2)%stable_phr(jj)=lfix
          newnode%linehead(2)%nodfixph=lfix
          newnode%linehead(2)%nstabph=jj
       else
          kk=0
          do jj=1,stabph-1
! remove -phfix as stable phase
             if(mapline%stableph(jj)%ixphase.eq.kph .and.&
                  mapline%stableph(jj)%compset.eq.kcs) then
                kk=jj+1
             else
                kk=kk+1
             endif
             newnode%linehead(2)%stableph(jj)%ixphase=&
                  mapline%stableph(kk)%ixphase
             newnode%linehead(2)%stableph(jj)%compset=&
                  mapline%stableph(kk)%compset
             newnode%linehead(2)%stable_phr(jj)=mapline%stable_phr(kk)
          enddo
          newnode%linehead(2)%nodfixph=lfix
!          newnode%linehead(1)%nodfixtup=meqrec1%phr(lfix)%phtupix
          newnode%linehead(2)%nstabph=stabph-1
       endif
!
! exit 3 has PHFIX as linefix ----------------------------------
       newnode%linehead(3)%linefixph%ixphase=kph
       newnode%linehead(3)%linefixph%compset=kcs
       newnode%linehead(3)%linefix_phr=abs(phfix)
       do jj=1,stabph
          if(mapline%stableph(jj)%ixphase.eq.kph .and. &
               mapline%stableph(jj)%compset.eq.kcs) then
! exchange PHFIX for LFIX as stable phase
             newnode%linehead(3)%stableph(jj)%ixphase=iph
             newnode%linehead(3)%stableph(jj)%compset=ics
             newnode%linehead(3)%stablepham(jj)=zero
             newnode%linehead(3)%stable_phr(jj)=abs(phfix)
          else
             newnode%linehead(3)%stableph(jj)%ixphase=&
                  mapline%stableph(jj)%ixphase
             newnode%linehead(3)%stableph(jj)%compset=&
                  mapline%stableph(jj)%compset
             newnode%linehead(3)%stablepham(jj)=mapline%stablepham(jj)
             newnode%linehead(3)%stable_phr(jj)=mapline%stable_phr(jj)
          endif
       enddo
       newnode%linehead(3)%nodfixph=lfix
       newnode%linehead(3)%nstabph=stabph
!
       if(ocv()) then
          do jj=1,3
             write(*,356)jj,newnode%seqx,&
!                  newnode%linehead(jj)%linefixph%phaseix,&
                  newnode%linehead(jj)%linefixph%ixphase,&
                  newnode%linehead(jj)%linefixph%compset,&
                  newnode%linehead(jj)%nodfixph,&
                  newnode%linehead(jj)%nstabph,&
                  (newnode%linehead(jj)%stableph(kk)%ixphase,&
                  newnode%linehead(jj)%stableph(kk)%compset,&
                  kk=1,newnode%linehead(jj)%nstabph)
          enddo
356       format('Tie-line NOT in plane node exits: ',&
               i2,i3,i4,i2,i5,i3,10(i4,i2))
       endif
! the fix and stable phases must be copied to meqrec1 when line is started
!       write(*,*)'Created node with 2 exits: ',newnode%seqx,ceq%tpval(1)
390    continue
!---------------------------------------------------------------
! invariant isopleth, more than 3 exits
    case(4) ! isopleth invariants for isopleths, inveq
! number of stable phases equal to components+1
! number of adjacent regions with "components" stable phases is "components+1"
! number of exit lines are 2*(components+1) ?? limit to 8 (minus 1 for entry)
! each line has a fix phase and one of the phases is stable at the invariant
! (set as not "nodefix").  The remaining phases  are entered.
! Each phase is fix for two lines and "nodefix" for two others
! This is the way to generate the exit lines:
! - loop for all phases to set a phase fix (for two lines)
! - loop for the next two phases to set one phase not stable
! the remaining phases are set entered (amount?) generate a line startpoint
! take care of remobing line into the invariant
!
! How to know if the node is invariant? Gibbs phase rule, Degrees of freedom
! f = n + 2 - p
! where n is number of components, 2 if T and P variable, 1 if T or P variable,
! 0 if both T and P fixed, p is number of stable phases.
!       write(*,*)'SMP2A Generating exits from isopleth invariant',newnode%lines
! Two crossing lines, one in and 3 exits
! this is probably redundant, fixph already reset
! phfix is the new stable phase! Must be positive
! mapline is the just finished line
       if(meqrec1%nfixph.gt.0) then
!          write(*,*)'Invariant isopleth:',meqrec1%nfixph,phfix,mapline%nstabph
          meqrec1%fixph(1,meqrec1%nfixph)=0
          meqrec1%fixph(2,meqrec1%nfixph)=0
          meqrec1%phr(abs(phfix))%phasestatus=PHENTUNST
          meqrec1%nfixph=meqrec1%nfixph-1
       endif
! determine LFIX, the phase which was fix along incomming line
       iph=mapline%linefixph(1)%ixphase
       ics=mapline%linefixph(1)%compset
! for use below I need to know the position of iph+ics in meqrec1%phr ...
       lfix=mapline%linefix_phr(1)
       flfix2: do jj=1,meqrec1%nstph
! this loop is only for stable phase it does not include the fix
          if(meqrec1%phr(jj)%iph.eq.iph .and. meqrec1%phr(jj)%ics.eq.ics) then
             lfix=jj
             meqrec1%phr(lfix)%phasestatus=PHENTUNST
             exit flfix2
          endif
       enddo flfix2
       if(lfix.eq.0) stop 'ERROR'
! this is total number of phases at each the invariant
! 1 fix and stabph-2 should be stable at each exit
       stabph=mapline%nstabph
       if(stabph.eq.0) then
          write(*,*)'Problems, no entered phase!'
          gx%bmperr=4242; goto 490
       endif
! Collect all stable phases to be used as different exits.
! invph(1,jj) is iph,, invph(2,jj) is ics; invph(3,jj) is index in meqrec1%phr
! invph(4,jj) is to count number of times jj has been linefix
! invph(5,jj) is to count number of times jj has been nodefix
! invph(6,jj) is index to phase_varres
       allocate(invph(6,stabph+2))
       invph=0
       do jj=1,stabph
! stableph is a phase_tuple
          invph(1,jj)=mapline%stableph(jj)%ixphase
          invph(2,jj)=mapline%stableph(jj)%compset
          invph(3,jj)=mapline%stable_phr(jj)
! stable_phr is used to find the index in phr and index to phase_varres
! I DO NOT TRUST THE VALUE, "stable_phr"
! I SHOULD REORGANIZE PHE TO BE IN PHASE TUPLE ORDER.
! THERE ARE ALWAYS PROBLEM IS IF NEW COMPOSIION SETS ARE CREATED DURING MAPPING
          do zz=1,meqrec1%nphase
             if(meqrec1%phr(zz)%iph.eq.invph(1,jj) .and.&
                  meqrec1%phr(zz)%ics.eq.invph(2,jj)) then
                invph(3,jj)=zz
!                if(zz.ne.mapline%stable_phr(jj)) &
!                     write(*,*)'SMP correction: ',jj,zz,mapline%stable_phr(jj)
             endif
          enddo
          call get_phase_compset(invph(1,jj),invph(2,jj),lokph,lokcs)
          if(gx%bmperr.ne.0) goto 1000
          invph(6,jj)=lokcs
       enddo
! at the end of loop jj=stabph+1; store phfix, the new fix phase, 
       invph(1,jj)=meqrec1%phr(phfix)%iph
       invph(2,jj)=meqrec1%phr(phfix)%ics
       invph(3,jj)=phfix
       call get_phase_compset(invph(1,jj),invph(2,jj),lokph,lokcs)
       if(gx%bmperr.ne.0) goto 1000
       invph(6,jj)=lokcs
! this is the phase fix at incomming line, only one exit line with this fix
       invph(1,jj+1)=iph
       invph(2,jj+1)=ics
       invph(3,jj+1)=lfix
       call get_phase_compset(invph(1,jj+1),invph(2,jj+1),lokph,lokcs)
       if(gx%bmperr.ne.0) goto 1000
       invph(6,jj+1)=lokcs
       jlast=stabph+2
! STABLE PHASES HAS TO BE IN PHR ORDER!! SORT invph
!       do kk=1,stabph+2
!          write(*,'(a,i3,2x,i3,i2,4i4)')'SMP invph:',kk,(invph(zz,kk),zz=1,6)
!       enddo
! second argument is first dimenstion of invph!!
       call sort_invph(jlast,6,invph)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'SMP sorted invph: '
!       do kk=1,stabph+2
!          write(*,'(a,i3,2x,i3,i2,4i4)')'SMP invph:',kk,(invph(zz,kk),zz=1,6)
!       enddo
       do jj=1,jlast
          phases=' '
          call get_phase_name(invph(1,jj),invph(2,jj),phases)
! keep track of the phase found at the invariant and the linefix phase
! nodein(1) is linefix
          if(invph(3,jj).eq.lfix) then
             nodein(1)=jj
             linefix=jj
             invph(4,nodein(1))=1
! nodein(2) is nodefix
          elseif(invph(3,jj).eq.phfix) then
             nodein(2)=jj
             nodefix=jj
             invph(5,nodein(2))=1
          endif
       enddo
! the entering line found at node, nodefix, mark it is used
! the entering line had this phase fix with zero amount
!       write(*,'(a,5i4)')'SMP2: linefix and nodefix: ',&
!            onlyone,lfix,notone,phfix
! all the others should be fixed one two exits
!-------------- 
       tmpceq=>newnode%nodeceq
! max 20 exit lines ....
       allocate(nodeout(2,10))
!       write(*,*)'SMP call to find all exits',tmpceq%tpval(1)
       call find_inv_exits(nexit,nodeout,nodein,stabph,invph,6,axarr,tmpceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,'(a,i3,5(i5,i3))')'SMP back from find_inv_exit: ',nexit,&
!            nodein(1),nodein(2),(nodeout(1,jj),nodeout(2,jj),jj=1,nexit)
!       stop 'SMP does it work? YES!'
       if(nexit.gt.2*mapline%nstabph) then
          write(*,*)'SMP too many exit lines: ',nexit,newnode%lines
       elseif(nexit.le.0) then
          write(*,*)'SMP no exits found?, just continue with one line'
          newnode%lines=1
       else
          newnode%lines=2*nexit+1
       endif
! There are nexit pairs of phases in nodeout for all exits, total number
! of exits are 2*nexit+1  (one exit eliminated because that was entering)
!
!-------------- 
! We have to generate newnode%lines exits!!
       jphr=mapline%nstabph
       jfix=1
! set bit in mapnode!
       if(newnode%status.ne.0) write(*,*)'SMP2 nodestatus: ',newnode%status
       newnode%status=ibset(newnode%status,MAPINVARIANT)
!       write(*,*)'SMP2 number of exit lines: ',newnode%lines,jphr
       allexit: do jj=1,newnode%lines
! initiate common data in map_line record in all exit lines
          newnode%linehead(jj)%number_of_equilibria=0
          newnode%linehead(jj)%first=0
          newnode%linehead(jj)%last=0
          newnode%linehead(jj)%axchange=-1
          newnode%linehead(jj)%done=1
          newnode%linehead(jj)%status=0
          newnode%linehead(jj)%more=1
          newnode%linehead(jj)%termerr=0
          newnode%linehead(jj)%axfact=1.0D-2
          newnode%linehead(jj)%axandir=stepax
! this dimensioning is OK for two axis, if 3 axis it should be 2 etc.
          newnode%linehead(jj)%nfixphases=1
          if(allocated(newnode%linehead(jj)%linefixph)) then
             write(*,*)'SMP2A line 5537: Strange allocated error in map17',&
                  jj,jphr
             deallocate(newnode%linehead(jj)%linefixph)
             deallocate(newnode%linehead(jj)%linefix_phr)
             if(allocated(newnode%linehead(jj)%stableph)) then
                write(*,*)'SMP2A line 5537: skipping!'
             endif
          endif
          allocate(newnode%linehead(jj)%linefixph(1))
          allocate(newnode%linehead(jj)%linefix_phr(1))
! There will be the same number of stable phases in all lines 
          allocate(newnode%linehead(jj)%stableph(jphr))
          allocate(newnode%linehead(jj)%stablepham(jphr))
          allocate(newnode%linehead(jj)%stable_phr(jphr))
! a small first step in same axis as used to find the node 
! We may have to change direction, in particular if the nodephase reappears
! evenvalue important only for STEP with one axis
          newnode%linehead(jj)%firstinc=axinc1*axarr(abs(stepax))%axinc
          newnode%linehead(jj)%evenvalue=zero
! links to node records at start and end of line
          newnode%linehead(jj)%start=>newnode
          nullify(newnode%linehead(jj)%end)
! number of stable phases along all lines.  Additionally a fix and a forbidden
          newnode%linehead(jj)%nstabph=stabph
! possible problem with meqrec%status
          if(newnode%linehead(jj)%meqrec%status.ne.0) then
             write(*,*)'SMP zero meqrec%status for newnode%linehead',&
                  newnode%linehead(jj)%meqrec%status
             newnode%linehead(jj)%meqrec%status=0
          endif
       enddo allexit
!
! ------------------------ ISOPLETHAL INVARIANTS EXITS -----------
! we have set LINEFIX and NODEFIX phases for each line
! We know the LINEFIX and NODEFIX phases for the line INTO THE INVARIAT
! For the first exit line we just swich these as they are at the same point
! with zero amount of both phases
!
!  (C is ?) Cfix         Afix  Dfix       Bfix (B is ?)
!            \   ABCE..   \BCE./  BCDE.. /
!             \            \  /         /
!              \____________\/_________/  
!      ABDE..  !_________ ABCDE..______!   CDE...
!              /            /\         \
!             /  ABDE..    /  \  ACDE.. \
!            /            /ADE.\         \
!           Dfix         Bfix   Cfix      Afix
!
! For the other exits FIND_INV_EXITS above have found all points along
! the invariant line that have two phases with zero amount.
! If next<0 just one exit will be generated with linefix/nodefix changed
!
       jj=1
       phases=' '
!       write(*,717)newnode%nodeceq%tpval(1)
717    format(/' *************** invariant node at ',F10.2)
! now code to create correct combination of linefix and nodefix
! first a line with nodefix as linefix and vice versa
! and old linefix as nodefix phase (sorry very confusing for me too)
       newnode%linehead(jj)%linefixph%ixphase=invph(1,nodefix)
       newnode%linehead(jj)%linefixph%compset=invph(2,nodefix)
       newnode%linehead(jj)%linefix_phr=invph(3,nodefix)
! this is just to understand what is happening
       call get_phase_name(invph(1,nodefix),invph(2,nodefix),phases)
       zp=len_trim(phases)+3
       phases(zp-1:zp-1)='('
!       write(*,*)'smp2: fix: ',trim(phases),onlyone,invph(3,onlyone)
! set linefix=LFIX as phase forbidden to become stable when line starts
! enclose the nodefix phase with ( .... )
       newnode%linehead(jj)%nodfixph=invph(3,linefix)
       call get_phase_name(invph(1,linefix),invph(2,linefix),phases(zp:))
       zp=len_trim(phases)+3
       phases(zp-2:zp-2)=')'
! at this line we have switched nodefix/linefix; nodefix is fix along the line
! the incomming line had the same fix phases so set both 4 and 5 to 1
       invph(4,nodefix)=1; invph(5,nodefix)=1
       invph(4,linefix)=1; invph(5,linefix)=1
! NOW add the stable phases excluding linefix and nodefix
       kk=0
       names: do zz=1,stabph+2
          if(zz.eq.linefix .or. zz.eq.nodefix) cycle names
          kk=kk+1
          newnode%linehead(jj)%stableph(kk)%ixphase=invph(1,zz)
          newnode%linehead(jj)%stableph(kk)%compset=invph(2,zz)
          newnode%linehead(jj)%stablepham(kk)=1.0D-2
          newnode%linehead(jj)%stable_phr(kk)=invph(3,zz)
          call get_phase_name(invph(1,zz),invph(2,zz),phases(zp:))
          zp=len_trim(phases)+2
       enddo names
! note again: nodefix here is fix along the line, linefix is stable at invarant
!       write(*,430)jj,nodefix,linefix,trim(phases)
       write(*,430)jj,trim(phases)
430    format('SMP2A invexit ',i3,' >>> ',a)
! The code above is for the FIRST exit line
! Here we will create 2 exit lines using nodeout(1,jj) and nodeout(2,jj)
! with jj=1..nexit, repeat with switched linefix/nodefix
       do qq=1,nexit
          linefix=nodeout(1,qq)
          nodefix=nodeout(2,qq)
          doubline=0
392       continue
          jj=jj+1
          newnode%linehead(jj)%linefixph%ixphase=invph(1,linefix)
          newnode%linehead(jj)%linefixph%compset=invph(2,linefix)
          newnode%linehead(jj)%linefix_phr=invph(3,linefix)
          invph(4,linefix)=invph(4,linefix)+1
          call get_phase_name(invph(1,linefix),invph(2,linefix),phases)
          zp=len_trim(phases)+3
          phases(zp-1:zp-1)='('
          newnode%linehead(jj)%nodfixph=invph(3,nodefix)
          invph(5,nodefix)=invph(5,nodefix)+1
          call get_phase_name(invph(1,nodefix),invph(2,nodefix),phases(zp:))
          zp=len_trim(phases)+3
          phases(zp-2:zp-2)=')'
          kk=0
          names2: do zz=1,stabph+2
             if(zz.eq.linefix .or. zz.eq.nodefix) cycle names2
             kk=kk+1
             if(kk.le.stabph) then
                newnode%linehead(jj)%stableph(kk)%ixphase=invph(1,zz)
                newnode%linehead(jj)%stableph(kk)%compset=invph(2,zz)
                newnode%linehead(jj)%stablepham(kk)=1.0D-2
                newnode%linehead(jj)%stable_phr(kk)=invph(3,zz)
             else
                write(*,'(a,10i5)')'SMP2 too many stable phases: ',jj,kk,zz,&
                     invph(1,zz),invph(2,zz),linefix,nodefix
             endif
             call get_phase_name(invph(1,zz),invph(2,zz),phases(zp:))
             zp=len_trim(phases)+2
          enddo names2
          write(*,430)jj,trim(phases)
          if(doubline.eq.0) then
! we have switch linefix and nodefis to create one more exit line
             doubline=linefix
             linefix=nodefix; nodefix=doubline
             phases=' '
             goto 392
          endif
       enddo
490    continue
!       stop ' *** Unfinished invariant isopleth node exits *** '
    end select exits
!=========================================================================
    goto 1000
!------------------------------------------- 
! we have found a node with same chemical potentials
! we should perhaps also check the set of phases ... ???
800 continue
    if(ocv()) write(*,*)'This node already found',mapnode%seqx
! we set a link in the mapline record to this node and has finished!
    mapline%end=>mapnode
    if(ocv()) write(*,*)'Line: ',mapline%lineid,' ends in node: ',mapnode%seqx
! >>> We must also mark the "%done=-1" in the linehead record corresponding to
! the line we just followed. 
!    
1000 continue
    return
  end subroutine map_newnode ! redefined argument mecreq to mecreq1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine sort_invph
!\begin{verbatim}
  subroutine sort_invph(nitems,ndim,array)
! primitive sorting of array
    implicit none
    integer nitems,ndim,array(ndim,*)
!\end{verbatim}
! sort array in acending order of value in array(3,*)    
    integer ia,ib,ic,more
    more=nitems
    ia=1
    do while(more.gt.0)
       more=0
       do ia=2,nitems
          if(array(3,ia-1).gt.array(3,ia)) then
             more=more+1
             do ic=1,ndim
                ib=array(ic,ia-1); array(ic,ia-1)=array(ic,ia); array(ic,ia)=ib
             enddo
          endif
       enddo
    enddo
1000 continue
  end subroutine sort_invph

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_inv_exits
!\begin{verbatim}
  subroutine find_inv_exits(nexit,phut,phin,stabph,invph,dim1,axarr,thisceq)
! Find the phases to be linefix when nodefix with zero amount at the invariant
! NEW IDEA (N is numbe of elements)
! 1. extract the composition of all stable phases at invariant (N+1)
! 2. set up a system of linear equations M_j x_ij = c_i
!    where x_ij is composition of component i in phase j,  c_i is the condition
!    for component i, (N-1 conditions) and M_j amount of phase j
! 3. The cases when his system has a solution represent exits!
! This also solves the problem with the number of exits!
! The conditions must be simple, such as x(cr)=0.05 ...
! Normal conditions are T, P, N and x, with N-1 x conditions.
! For 5 components, 6 phases at invariant, 4 mass balance conditions
! input/output:
! nexit is number of pairs that form two exit lines
! phut(1,*) and (2,*) are two phases with zero amount at the invariant
! phin is on enter the two phase with zero amount when the invariant is found
! stabph is number of stable phases along the lines, at node point 2 more stable
! invph is matrix with all phases at the node, dim1 is its first dim
! dim1 is first dimension of invph
! axarr has axis information (needed to maniplulate conditions)
! eqcopy1 is array with phase anounts and constitutions at new point (not used)
! linerec is a line record with all necessary data to calculate en equil
    implicit none
    integer nexit,phut(2,*),phin(2),stabph,dim1,invph(dim1,*),par(2)
    type(map_axis), dimension(*) :: axarr
    double precision, dimension(:), allocatable :: eqcopy1
    type(gtp_equilibrium_data), pointer :: thisceq
!    type(map_line), pointer :: linerec
!\end{verbatim}
    double precision, allocatable, dimension(:) :: condval,rhs
    double precision, allocatable, dimension(:,:) :: phaseval,test
    type(gtp_condition), pointer :: lastcond,pcond,axcond
    type(gtp_state_variable), pointer :: statevar
! assume less than 20 components ...
    integer, parameter :: mcomp=20
    character text*32,ch1*1
    integer, allocatable, dimension(:) :: ipiv,jphase
    type(gtp_state_variable), dimension(mcomp), target :: stvarray
!
    integer ii,jj,kk,mm,ip,seqz,ncomp,ncomp1,info,ldb,zz
! dim is number of phases at the invariant, dim-1 is number of components
! For 5 components, 6 phases at invariant, select 4 phases to find exit
! x_ij is fraction of j in phase i
!    x_11 + x_21 + x_31     N_1       C_1    
!  ( x_12 + x_22 + x_32 ) ( N_2 ) = ( C_2 )    
!    x_13 + x_23 + x_33     N_3       C_3    
!    1      1      1        N_4       1
! to find N_i (which must be >0).  The two excluded phases represet the exit
!
! 1. extract all condition values, skip the axis condition
    allocate(condval(stabph-1))
    lastcond=>thisceq%lastcondition
    pcond=>lastcond%next
    ncomp=0
    cloop: do while(.true.)
       cskip: if(pcond%active.eq.0) then
! this is an active condition, extract the state variable record
          statevar=>pcond%statvar(1)
          if(statevar%statevarid.ge.10 .and. statevar%argtyp.eq.1) then
! statevarid>=10 is extensive condition on a component
             seqz=pcond%seqz
! skip axis conditions
             if(seqz.eq.axarr(1)%seqz .or. seqz.eq.axarr(1)%seqz) exit cskip
! There must not be any terms
             if(pcond%noofterms.gt.1) cycle cloop
             ip=1
             ncomp=ncomp+1
! remember the state variable to be used extracting values from phases
             stvarray(ncomp)=statevar
             if(ncomp.gt.mcomp) stop 'SMP Too many components'
             condval(ncomp)=pcond%prescribed
!             call get_one_condition(ip,text,seqz,thisceq)
!             write(*,*)'SMP condition: ',trim(text),ncomp,condval(ncomp)
          endif
       endif cskip
       if(associated(pcond,lastcond)) exit cloop
       pcond=>pcond%next
    enddo cloop
! we must have extracted stabph-1 extensive conditions ...
    if(ncomp.ne.stabph-1) then
       write(*,*)'SMP too few conditions for invariants',ncomp,stabph-2
       gx%bmperr=4399; goto 1000
    endif
! 2. extract all phase compositions, for stabph+2 phases and ncomp compositions
    allocate(phaseval(ncomp,stabph+2))
    do ii=1,stabph+2
       do jj=1,ncomp
! insert phase index to get phase composition
          stvarray(jj)%argtyp=3
          stvarray(jj)%phase=invph(1,ii)
          stvarray(jj)%compset=invph(2,ii)
! this subroutine uses character as argument
!          call get_state_var_value(stvarray(jj),phaseval(jj,ii),text,thisceq)
! this subroutine uses state variable as argument
          statevar=>stvarray(jj)
          call state_variable_val(statevar,phaseval(jj,ii),thisceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'SMP error extracting phase value',trim(text),ii,jj
             goto 1000
          endif
       enddo
    enddo
!----------------------------------
! This is the invariant used to debug this 
!Output for equilibrium:   1, DEFAULT_EQUILIBRIUM          2020.04.28
!Conditions .................................................:
!  1:P=100000, 2:N=1, 3:W%(CR)=5, 4:W%(MO)=8, 5:W%(V)=1, 6:<M23C6>=0,
!    7:<KSI_CARBIDE>=0
! Degrees of freedom are   0
!
!Some global data, reference state SER ......................:
!T=   1232.69 K (   959.54 C), P=  1.0000E+05 Pa, V=  6.4311E-06 m3
!N=   1.0000E+00 moles, B=   5.4289E+01 g, RT=   1.0249E+04 J/mol
!G= -6.09563E+04 J, G/N=-6.0956E+04 J/mol, H= 3.4334E+04 J, S= 7.730E+01 J/K
!
!Some data for components ...................................:
!Component name    Moles      Mass-fr  Chem.pot/RT  Activities  Ref.state
!C                 7.1177E-02  0.01575 -3.2298E+00  3.9566E-02  SER (default)
!CR                5.2205E-02  0.05000 -7.7539E+00  4.2908E-04  SER (default)
!FE                8.2069E-01  0.84425 -5.8469E+00  2.8888E-03  SER (default)
!MO                4.5269E-02  0.08000 -8.0031E+00  3.3442E-04  SER (default)
!V                 1.0657E-02  0.01000 -1.4256E+01  6.4376E-07  SER (default)
!
!Some data for phases .......................................:
!Name                Status Mass       Volume    Form.Units Cmp/FU dGm/RT  Comp:
!FCC_A1#1................ E  4.802E-02  6.36E-06  8.48E-01    1.03  0.00E+00  W:
! FE     9.34742E-01  MO     2.09038E-02  C      7.38762E-03  V      1.04403E-03
! CR     3.59223E-02
!
!MC_FCC_A1#2............. E  7.219E-04  6.12E-09  9.35E-03    1.86  0.00E+00  W:
! MO     4.38724E-01  C      1.33163E-01  CR     4.68908E-02  FE     6.00860E-03
! V      3.75214E-01
!
!HCP_A3.................. E  3.954E-03  6.74E-08  4.54E-02    1.50  0.00E+00  W:
! MO     7.20128E-01  FE     7.62966E-02  C      6.87569E-02  V      4.49350E-02
! CR     8.98830E-02
!
!KSI_CARBIDE............. F  0.000E+00  0.00E+00  0.00E+00    4.00  0.00E+00  W:
! FE     4.58671E-01  CR     9.59266E-02  C      5.64977E-02  V      0.00000E+00
! MO     3.88905E-01
!
!M23C6................... F  0.000E+00  0.00E+00  0.00E+00   29.00  0.00E+00  W:
! FE     5.94954E-01  MO     1.04967E-01  C      5.17714E-02  V      1.20517E-04
! CR     2.48187E-01
!
!M7C3.................... E  1.590E-03  0.00E+00  3.67E-03   10.00  0.00E+00  W:
! FE     4.01696E-01  MO     1.09953E-01  C      8.31013E-02  V      2.77992E-02
! CR     3.77451E-01
!
! In alphabetical order of phases and components (mass percent)
! Phase   Cr        Mo        V          Mass of phase ??
! FCC      3.59223   2.09038   0.1044    0.04802
! MC_FCC   4.68908  43.8724   37.5214    0.0007219
! HCP      8.98830  72.0128    4.4925    0.003954
! KSI      9.59266  38.8905    0         0
! M23     24.8187   10.4967    0.00012   0
! M7      37.7451   10.9953    2.77992   0.001590    
!-----------------------------------
! debug:
!    do ii=1,stabph+2
!       write(*,88)'SMP all: ',ii,(phaseval(jj,ii),jj=1,ncomp)
!    enddo
88  format(a,i3,6(1pe12.4))
! The code gives: (in mass percent)  
! 3.5822 2.0904 0.1044 etc ...
! ALL CORRECT!! WOW
! (missing in list result is mass %! I have only mass fraction)
!
! 3. select submatrix with dim-2 phases and solve for phase fractions.
!    For solutions with phase phase fractions >0 the 2 excluded phases
!    are exits.
    nexit=0
    ncomp1=ncomp+1
    allocate(test(ncomp1,ncomp1))
    allocate(rhs(ncomp1))
    allocate(ipiv(ncomp1))
! I am not certain of this dimensioning ...
    allocate(jphase(ncomp1*(ncomp1+1)))
!    allocate(lukas(ncomp1+1,ncomp1))
! All possible ncomp x ncomp marices from phaseval are solved for phase amounts
! for the correct content of the components.  One should find phin!
! This means we actually have 3 phases with zero amount at the lines??
! The matrix phaseval has stabph+2 rows and columns
! We must copy this to test eliminating 3 rows
! THIS WAY OF GENERATING ALL COMBINATIONS OF ncomp x nacomp MATRICES IS
! involved but seems to work .... WoW
    do ii=1,ncomp1
       jphase(ii)=ii
    enddo
    zz=0
    kloop: do while(.true.)
       zz=zz+1
!       write(*,'(a,10i3)')'SMP subset ----------- ',zz,jphase
! test is destroyed when solving the system of linear equations
! and must be regenerated totally each time
       test=zero
! copy fractions from phaseval(P,*) represent fractions in phase P
! test(*,K) is fractions in all phases for component K
       do jj=1,ncomp1
          do kk=1,ncomp
             test(kk,jj)=phaseval(kk,jphase(jj))
          enddo
       enddo
! last line should be a row of 1.0
       do kk=1,ncomp1
          test(ncomp1,kk)=one
       enddo
!       do jj=1,ncomp1
!          write(*,88)'SMP sub: ',zz,(test(jj,kk),kk=1,ncomp1)
!       enddo
!------------------------
! this should be the solving ...
! LAPACK routine to L*U factorize A, the original A is destroyed
!    call dgetrf(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
! ipiv is array with N pivot 
       ldb=ncomp1
       call dgetrf(ncomp1,ncomp1,test,ldb,ipiv,info)
       if(info.ne.0) then
          write(*,*)'SMP error from dgetrf',info
! some combination of phases may not work, just skip
          goto 100
!          gx%bmperr=4399; goto 1000
       endif
! solve the system of linear equations, X is overwritten by solution
!    call dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
       do kk=1,ncomp
          rhs(kk)=condval(kk)
       enddo
       rhs(ncomp1)=one
!       write(*,'(a,10(1pe10.2))')'SMP rhs: ',rhs
       call dgetrs('N',ncomp1,1,test,ldb,ipiv,rhs,ldb,info)
       if(info.ne.0) then
! some combination of phases may not work, just skip
          goto 100
!          write(*,*)'SMP error from dgetrs',info
!          gx%bmperr=4399; goto 1000
       endif
600    continue
!       write(*,'(a,10(1pe10.2))')'SMP phase amounts: ',rhs
! check if all amounts greater than zero
       do kk=1,ncomp1
          if(rhs(kk).le.zero) goto 100
       enddo
!       write(*,'(a,10(1pe11.3))')'SMP phase amounts: ',rhs
! Wow, now it works, but I must find which phases are excluded
! A very clumsy set to find which two phases that are excluded ...
       zz=0
       ex1: do jj=1,stabph+2
          do kk=1,stabph+2
             if(jphase(kk).eq.jj) cycle ex1
          enddo
          zz=zz+1; par(zz)=jj
       enddo ex1
!       write(*,'(a,20i3)')'SMP solution ',zz,jphase,0,phin,par
!       write(*,'(a,20i3)')'SMP solution ',zz,phin,par
! check if solution equal to phin
       if((par(1).eq.phin(1) .or. par(1).eq.phin(2)) .and. &
            (par(2).eq.phin(1) .or. par(2).eq.phin(2))) then
          continue
!          write(*,'(a,2i3,3x,2i2)')'SMP same as phin: ',par,phin
       else
          nexit=nexit+1
          phut(1,nexit)=par(1)
          phut(2,nexit)=par(2)
       endif
!       read(*,'(a)')ch1
! here we have solved the system of linear equations
!------------------------
100    continue
! exclude a different phase in jphase ...
       jj=ncomp1
       kk=0
       jloop: do while(.true.)
          jphase(jj)=jphase(jj)+1
          if(jphase(jj).gt.stabph+2-kk) then
             jj=jj-1
             kk=kk+1
             if(jj.ge.1) cycle jloop
             exit kloop
          else
             exit jloop
          endif
! increment all values in jphase after jj
       enddo jloop
       do kk=jj+1,ncomp1
          jphase(kk)=jphase(kk-1)+1
       enddo
    enddo kloop
! now we have found the all exits ...
!    write(*,'(a,i3,5(i5,i3))')'SMP exits: ',nexit,&
!         (phut(1,jj),phut(2,jj),jj=1,nexit)
!
1000 continue
    return
  end subroutine find_inv_exits

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_inv_nodephase
!\begin{verbatim}
  subroutine find_inv_nodephase(phut,phin,stabph,invph,dim1,axarr,&
       eqcopy1,linerec)
! Find the phase to set as linefix when nodefix has zero amount
! phut will on exit be two phases with zero amount at the invariant
! phin is on enter the two phase with zero amount when the invariant is found
! stabph is number of stable phases along the lines, at node point 2 more stable
! invph is matrix with all phases at the node, dim1 is its first dim
! axarr has axis information (needed to manipulate conditions)
! eqcopy1 is array with phase anounts and constitutions at new point (not used)
! linerec is a line record with all necessary data to calculate en equil
    implicit none
    integer phut(2),phin(2),stabph,dim1,invph(dim1,*)
    type(map_axis), dimension(*) :: axarr
    double precision, dimension(:), allocatable :: eqcopy1
    type(map_line), pointer :: linerec
!\end{verbatim}
! linefix and nodefix is 2nd index in invph(5,*)
    type(gtp_condition), pointer :: lastcond,pcond,pcondx,pcondt
!    type(gtp_equilibrium_data), target :: thisceq1
    type(gtp_equilibrium_data), pointer :: thisceq
    integer linefix,nodefix,reset,ph1,ph2,lokcs1,lokcs2,onemoretry
    integer ii,jj,kk,nodeph,mapx,iadd,irem,iz,jax,mode,okcond,lokph,zeroam
    character*24 phname1,phname2,phname3,phname4
    double precision, parameter :: mamfu=1.0D-6,mdgm=1.0D-4
    double precision, dimension(:), allocatable :: eqcopy
!
! we must set the axis condition on T and
! remove the axis condition on the composition on the other axis
!-------------------- copied from somewhere ...................
!    write(*,10)phin,stabph,linerec%lineceq%tpval(1)
10 format(/'SMP find_inv exits from isopleth invariant: ',2i4,2x,i2,F10.2)
    reset=globaldata%status
! supress messages from calceq3 done inside find_inv
    globaldata%status=ibset(globaldata%status,GSSILENT)
    globaldata%status=ibclr(globaldata%status,GSSILENT)
! I wonder how much is copied when I use = ...?
! lineceq is an equilibrium record in eqlista ...
    thisceq=>linerec%lineceq
    write(*,*)'SMP exit equilibrum: ',thisceq%eqname
! constitution is OK here
!    call list_conditions(kou,linerec%lineceq)
!    call list_sorted_phases(kou,thisceq)
!    write(*,*)'SMP constitution entering find_inv'
! suspend all phases not involved in the invariant, 1 means suspend
    call suspend_somephases(1,invph,6,stabph+2,thisceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'SMP error calling suspend_somephases'; goto 1000
    endif
    if(allocated(eqcopy)) deallocate(eqcopy)
    call save_phase_constitutions(0,thisceq,eqcopy)
    if(gx%bmperr.ne.0) goto 1000
! loop both axis to extract condition pointer
    okcond=0
    do jax=1,2
! Set the condition on T and remember the condition on composition
       lastcond=>thisceq%lastcondition
       if(.not.associated(lastcond)) then
          write(*,*)'in find_inv, no conditions: ',jax
          gx%bmperr=4221; goto 1000
       endif
       pcond=>lastcond
60     continue
       pcond=>pcond%next
       if(pcond%seqz.eq.axarr(jax)%seqz) goto 70
       if(.not.associated(pcond,lastcond)) goto 60
       write(*,*)'in find_inv the axis condition not found: ',jax
       gx%bmperr=4221; goto 1000
!
70     continue
       if(pcond%statev.ge.10) then
! save pointer to extensive condition and remove it
          pcond%active=1
          okcond=okcond+1
          pcondx=>pcond
       elseif(pcond%statev.eq.1) then
! set condition on T as active
          okcond=okcond+1
          pcondt=>pcond
! try setting two fix phases .... remove condition on T
!          pcond%active=1
       endif
    enddo
    if(okcond.ne.2) then
       write(*,*)'Conditions not T and X, quitting'
       gx%bmperr=4399
       goto 1000
    endif
!----------------- end copy from somewhere.................
! most of these variables are just for debugging
    ph1=phin(1)
    ph2=phin(2)
! extract the name of the linefix and nodefix phases
    call get_phase_name(invph(1,ph1),invph(2,ph1),phname1)
    call get_phase_name(invph(1,ph2),invph(2,ph2),phname2)
! set small amounts of ph1 ??
!    call change_many_phase_status(phname1,PHENTSTAB,mamfu,thisceq)
!    call change_many_phase_status(phname1,PHENTSTAB,zero,thisceq)
!    if(gx%bmperr.ne.0) then
!       write(*,*)'SMP error setting zero of line/nodefix'
!       goto 1000
!    endif
! below we check the amounts and driving forces of both linefix and nodefix
! Note the phase_varres indices are the same in all equilibria!
    call get_phase_compset(invph(1,ph1),invph(2,ph1),lokph,lokcs1)
    call get_phase_compset(invph(1,ph2),invph(2,ph2),lokph,lokcs2)
    if(gx%bmperr.ne.0) then
       write(*,*)'SMP find_inv failed get phase_varres index'
       goto 1000
    endif
! remove the condition on compostion and set fix T
    pcondx%active=1
    pcondt%active=0
!-------------------------------------------------------
!
! the loop below tries set fix one phase at a time to discover an exit
! where one other phase is stable with zero amount
! Initial amounts of phases and constitutions are in eqcopy
! all calculations below are made at fixed T with different fix phases
! afterwards we check that the amount of nodefix is still zero (or very small)
! total number of phases at nodepoint is stabph+2
    call list_sorted_phases(kou,0,thisceq)
    write(*,*)'SMP above listing for initial set of stable phases',gx%bmperr
    loop: do ii=1,stabph+2
       unused: if(invph(4,ii).eq.0 .and. invph(5,ii).eq.0) then
! only test phases which has not already been used, extract its name
          call get_phase_name(invph(1,ii),invph(2,ii),phname3)
          onemoretry=0
! jump back here to try with a small amount of phname1 or phname2 (no good)
50        continue
! restore inital amounts and constitutions 1 copies from eqcopy to thisceq
          call save_phase_constitutions(1,thisceq,eqcopy)
          call list_sorted_phases(kou,0,thisceq)
          write(*,76)trim(phname3),ii,thisceq%tpval(1),gx%bmperr
76        format(/'SMP find_inv ******** testing as fixed: ',a,i4,F10.2,i7)
! set ii fix with zero amount
          call change_many_phase_status(phname3,PHFIXED,zero,thisceq)
! debug listing
          call list_conditions(kou,thisceq)
          write(*,*)'SMP find_inv call calceq3:',gx%bmperr,onemoretry
          mode=0
          call calceq3(mode,.FALSE.,thisceq)
          if(gx%bmperr.ne.0) then
! calculation error, remove phase ii as fix and try another
             write(*,*)'SMP find_inv error finding exit with fix: ',&
                  trim(phname3),gx%bmperr
             gx%bmperr=0
             call list_conditions(kou,thisceq)
             call list_sorted_phases(kou,0,thisceq)
! remove this phase as fix and continue loop
             goto 120
          endif
! debug listing
          call list_conditions(kou,thisceq)
          call list_sorted_phases(kou,0,thisceq)
          write(*,*)'SMP find_inv: phases AFTER calculations',gx%bmperr
          jj=0
          zeroam=0
          zeroloop: do kk=1,stabph+2
! At any exit all phases invph should have almost zero dgm
             if(abs(thisceq%phase_varres(invph(6,kk))%dgm).gt.mdgm) then
                write(*,*)'SMP too negative dgm: ',kk,&
                     thisceq%phase_varres(invph(6,kk))%dgm
!                exit unused
             endif
             if(thisceq%phase_varres(invph(6,kk))%amfu.lt.mamfu) then
! we have a phase with zero amount
                if(kk.ne.ii) then
                   if(zeroam.eq.0) then
                      zeroam=kk
                   else
                      write(*,*)'SMP two or more zero amount phases'
                   endif
                endif
             else
                jj=jj+1
             endif
          enddo zeroloop
!          write(*,*)'SMP number of stable phases',jj,stabph
!          call list_sorted_phases(kou,thisceq)
          if(jj.ne.stabph) then
             write(*,'(a,2i3,F10.2)')'SMP Skip wrong number of stable phases',&
                  jj,stabph+2,thisceq%tpval(1)
             goto 120
          endif
          write(*,'(a,3(i5,i3))')'SMP phin: ',jj,stabph,invph(4,phin(1)),&
               invph(5,phin(1)),invph(4,phin(2)),invph(5,phin(2))
! select either phin(1) or phin(2) if free and zero amount with this phase
          if(invph(4,phin(1)).eq.1 .and.&
               thisceq%phase_varres(lokcs1)%amfu.lt.mamfu) then
! phase invph(*,ii) and invph(*,phin(1)) have zero amount, use as exit?
!  ADD CHECK IF amfu for BOTH phin(1) and ph(2) are zero we must check dgm ...
             phut(1)=ii; phut(2)=phin(1)
             write(*,112)trim(phname3)//'+'//trim(phname1),phut(1),phut(2),&
                  thisceq%tpval(1)
112          format('SMP find_inv **** success: ',a,2i4,F10.2)
             goto 200
          elseif(invph(4,phin(2)).eq.1 .and. &
               thisceq%phase_varres(lokcs2)%amfu.lt.mamfu) then
! phase invph(*,ii) and invph(*,phin(1)) have zero amount, use this as exit
             phut(1)=ii; phut(2)=phin(2)
             write(*,112)trim(phname3)//'+'//trim(phname2),phut(1),phut(2),&
                  thisceq%tpval(1)
             goto 200
          elseif(invph(4,phin(1)).eq.1 .and. invph(4,phin(2)).eq.1) then
! This is first call to find_inv and we have found two new phases with
! zero amount.  We do not need to find any more!
             call get_phase_name(invph(1,zeroam),invph(2,zeroam),phname4)
! Indicae this by setting it negative!
             phut(1)=ii; phut(2)=-zeroam
             write(*,114)trim(phname3)//'+'//trim(phname4),phut(1),-phut(2),&
                  thisceq%tpval(1)
114          format('SMP find_inv **** success BUT IGNORED: ',a,2i4,F10.2)
          else
             call list_sorted_phases(kou,0,thisceq)
             write(*,113)trim(phname3),trim(phname1)//' nor '//trim(phname2)
113          format('SMP Skipping ',a,' as neither ',a,' has zero amount')
          endif
120       continue
! Failed but make two more tries.  Constitutions restored above
!          if(onemoretry.eq.0) then
! try once more ...
!             write(*,*)'SMP try one more time with fix ',trim(phname3)
!             call change_many_phase_status(phname1,PHENTSTAB,1.0D-2,thisceq)
!             onemoretry=1
!             goto 50
!          elseif(onemoretry.eq.1) then
! try once more ...
!             onemoretry=2
!             write(*,*)'SMP try one more time with fix ',trim(phname3)
!             call change_many_phase_status(phname2,PHENTSTAB,1.0D-2,thisceq)
!             goto 50
!          endif
! giv up on this fix phase
       endif unused
! try another fix phase ...
       call change_many_phase_status(phname3,PHENTERED,zero,thisceq)
       if(gx%bmperr.ne.0) goto 1000
    enddo loop
! we have not found any set of phases for an exit line
! if we arrive here we should maybe try 2 fix phases and release T?
    write(*,*)'SMP find_inv failed to find two phases with zero amount'
    gx%bmperr=4399; goto 1000
200 continue
! we have a pair of phases in phut, reset phname3 as entered
    call change_many_phase_status(phname3,PHENTERED,zero,thisceq)
! copy constitution from thisceq to eqcopy, then copy to linerec%lineceq
!   call save_phase_constitutions(1,thisceq,eqcopy)
    if(allocated(eqcopy)) deallocate(eqcopy)
    call save_phase_constitutions(0,thisceq,eqcopy)
    call save_phase_constitutions(1,linerec%lineceq,eqcopy)
    if(gx%bmperr.ne.0) write(*,*)'Problem to copy constitutions'
!----------------------------- exit
1000 continue
    ii=gx%bmperr; gx%bmperr=0
! restore axis conditions, set x condition and remove T condition
    pcondx%active=0
    pcondt%active=1
! If we have found an exit the phase set and constitution are in thisceq
! Restore phases earlier suspended, 0 menas restore
    call suspend_somephases(0,invph,6,stabph+2,thisceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'SMP error calling suspend_somephases'; goto 1000
    endif
    if(allocated(eqcopy)) deallocate(eqcopy)
    gx%bmperr=ii
! reset globaldata%status
    globaldata%status=reset
    return
  end subroutine find_inv_nodephase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine reserve_saveceq
!\begin{verbatim}
  subroutine reserve_saveceq(location,saveceq)
! must be THREADPROTECTED
! location index of reserved ceq record in saveceq
    implicit none
    integer location
    type(map_ceqresults), pointer :: saveceq
!\end{verbatim}
    location=saveceq%free
!    write(*,*)'SMP reserve record: ',location,saveceq%size
!    write(*,*)'Reserve place for equilibrium: ',location,saveceq%size
    if(location.eq.saveceq%size-10) then
! indicate overflow with 5 places left if some emergency saving needed
       write(*,*)'Close to overflow in saveceq: ',saveceq%free
       gx%bmperr=4219; goto 1000
    endif
    saveceq%free=location+1
! end THREADPROTECT
1000 continue
    return
  end subroutine reserve_saveceq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_findline
!\begin{verbatim}
  subroutine map_findline(maptop,axarr,mapfix,mapline)
! must be THREADPROTECTED
! Searches all node records from maptop for a map_line record to be calculated
! ?? already been found and if so eliminate a line record ??
! maptop map node record
! axarr array with axis records
! mapfix returned fixph record with phases to be ste as fixed for this line
! mapline returned mapline record for line to be calculated
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(map_axis), dimension(*) :: axarr
    type(map_fixph), allocatable :: mapfix
! memory leak as mapfix is allocated below ...
!    type(map_fixph), pointer :: mapfix
!\end{verbatim}
    type(map_node), pointer :: mapnode
    type(gtp_condition), pointer :: pcond
    type(meq_setup), pointer :: meqrec
    type(gtp_equilibrium_data), pointer :: ceq2
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(gtp_state_variable), target :: svrtarget
    integer nyline,jp,seqy,iph,ics,lokph,lokcs,ip,mapx
    integer mode,nofecond,jax,activax,irem1,iadd1,irem,iadd
    integer, parameter :: inmap=1
    double precision finc,fixpham,xstab,xfix,xcorr,natpermol
    character eqname*24
! sometimes there are many phases with long names ...
    character phaseset*512
    logical hullerombuller
!
    mapnode=>maptop
! for the moment skip this for tie-lines not in the plane
!    write(*,*)'In map_findline: ',mapnode%tieline_inplane
100 continue
!       write(*,*)'Looking a lines exiting from all nodes:'
!       write(*,*)'mapnode index: ',mapnode%seqx
       if(.not.allocated(mapnode%linehead)) then
          write(*,*)'ERROR found mapnode without exits'
          nullify(mapline)
          goto 1000
       endif
!       write(*,*)'map_findline: ',mapnode%lines
       do nyline=1,mapnode%lines
! if done is >=0 then this is a line to be calculated
          if(mapnode%linehead(nyline)%done.ge.0) then
             mapline=>mapnode%linehead(nyline)
             mapline%done=-1
             goto 500
          endif
       enddo
!       write(*,*)'findline 1: ',mapnode%seqx
       mapnode=>mapnode%next
!       write(*,*)'findline 2: ',mapnode%seqx
       if(.not.associated(mapnode,maptop)) goto 100
! no more lines to calculate
    if(ocv()) write(*,*)'nullifying mapline as no more lines to calculate'
    nullify(mapline); goto 1000
! jump here if we found a nyline
!--------------------------------------------------------------------
500 continue
! we must copy the equilibrium record ceq to the line record
    goto 503
!----------------------------------------------------------------------
! code deleted
!---------------------------------------------------------------------
503 continue
!    write(*,*)'At label 503',mapline%firstinc
    if(mapline%firstinc.ne.zero) then
! update the axis condition if mapline%firstinc is nonzero
       jp=abs(mapline%axandir) 
       call locate_condition(axarr(jp)%seqz,pcond,mapline%lineceq)
       if(gx%bmperr.ne.0) goto 1000
! Wow, here is the problem !!!
       pcond%prescribed=pcond%prescribed+mapline%firstinc
       finc=mapline%firstinc
       if(ocv()) write(*,501)'Selecting line and condition: ',&
            seqy,pcond%prescribed,finc
501    format(a,i3,2(1pe14.6))
!       if(pcond%active.ne.0) then
!          write(*,*)'Error: axis condition not active!'
!          pcond%active=0
!       endif
    endif
! check that correct axis condition is active, maybe I have not made sure
! that the map_line records are independent ...
!    write(*,*)'axies ',maptop%number_ofaxis
    do jp=1,maptop%number_ofaxis
       call locate_condition(axarr(jp)%seqz,pcond,mapline%lineceq)
!       write(*,*)'condition located',jp,axarr(jp)%seqz
       if(gx%bmperr.ne.0) goto 1000
       if(pcond%active.ne.0) then
          if(jp.eq.abs(mapline%axandir)) then
!             write(*,*)'Setting axis condition active!',jp,mapline%axandir
             pcond%active=0
!             write(*,*)'map_findline: ',jp,pcond%prescribed
          endif
       else
          if(jp.ne.abs(mapline%axandir)) then
!             write(*,*)'Setting axis condition NOT active!',jp,mapline%axandir
             pcond%active=1
          endif
       endif
    enddo
! we may have a set of stable phases in mapnode%stable_phases, maybe they
! should be set, at least when mapping. 
! a meqrec record will be created by calceq7 at the first calculation
! for mapping set values in mapfix about which phases that should be fix
! or stable when calling calceq7, at present ingnore that
!-------------------------------------------------------------
!    write(*,*)'tielines: ',maptop%tieline_inplane
    if(maptop%tieline_inplane.lt.0) then
! ISOPLETH
       if(allocated(mapfix)) then
          deallocate(mapfix)
       endif
       allocate(mapfix)
! with only 2 axis we have just 1 fix phase for mapping, fixph is a tuple
       allocate(mapfix%fixph(1))
       mapfix%status=0
       mapfix%nfixph=1
       mapfix%fixph=mapline%linefixph(1)
! we can have several stable phases when no tie-lines in plane
       ip=mapline%nstabph
       allocate(mapfix%stableph(ip))
       allocate(mapfix%stablepham(ip))
       allocate(mapfix%stable_phr(ip))
!       write(*,*)'Findline: Tie-lines not in plane: ',nyline,ip
! create a heading text for the line
       phaseset=' '
       call get_phasetuple_name(mapfix%fixph(1),phaseset)
       if(gx%bmperr.ne.0) goto 1000
       ip=len_trim(phaseset)
       phaseset(ip+1:ip+10)=', stable: '
       ip=len_trim(phaseset)+2
       if(mapnode%linehead(nyline)%nstabph.le.0) then
          write(*,*)'No stable phases for a line'
          write(*,*)'Error 14:  ',nyline,mapnode%linehead(nyline)%nstabph,&
               mapnode%linehead(nyline)%stableph(1)%ixphase,&
               mapnode%linehead(nyline)%stableph(1)%compset
          mapfix%nstabph=0
          gx%bmperr=4242; goto 1000
       endif
       mapfix%nstabph=mapnode%linehead(nyline)%nstabph
!       write(*,*)'Findline: stable phases: ',mapfix%nstabph
       do jp=1,mapfix%nstabph
! this is stored only for "real" nodes
          mapfix%stableph(jp)=mapnode%linehead(nyline)%stableph(jp)
          mapfix%stable_phr(jp)=mapnode%linehead(nyline)%stable_phr(jp)
          call get_phasetuple_name(mapfix%stableph(jp),phaseset(ip:))
          if(gx%bmperr.ne.0) goto 1000
! this values hould perhaps be in linehead??
!          mapfix%stablepham(jp)=mapnode%linehead(nyline)%stablepham(jp)
          mapfix%stablepham(jp)=one
!          call get_phase_compset(mapfix%stableph(1)%phase,&
!               mapfix%stableph(1)%compset,lokph,lokcs)
!          mapline%lineceq%phase_varres(lokcs)%amfu=one
          ip=len_trim(phaseset)+2
       enddo
       write(kou,520)mapline%lineid,mapline%lineceq%tpval(1),phaseset(1:ip)
520    format(/'Line ',i3,' T=',F8.2,' fix: ',a)
!-------------------------------------------------------------
    elseif(maptop%tieline_inplane.gt.0) then
! TIE-LINES IN PLANE, NOTE: meqrec not allocated!!
!       if(mapnode%nodefix%phase.gt.0) then
!          write(*,*)'We have to set the fix phase along the line: ',&
!               mapnode%nodefix%phase,mapnode%nodefix%compset
! phr no longer allocated ...
!          iph=mapnode%meqrec%phr(mapnode%fixph)%iph
!          ics=mapnode%meqrec%phr(mapnode%fixph)%ics
!          call get_phase_compset(iph,ics,lokph,lokcs)
!          write(*,*)'Setting phase status as fixed, phase_varres: ',lokcs
!          mapline%lineceq%phase_varres(lokcs)%status2=&
!               ibset(mapline%lineceq%phase_varres(lokcs)%status2,PHFIXED)
!       endif
! mapline here should be identical to mapnode%linehead(nyline)
!       if(ocv()) write(*,505)'In findline: add phase set for',&
!       write(*,505)'In findline: add phase set for',&
!            ' tie-lines in plane, node:',&
!            mapnode%seqx,nyline,mapnode%linehead(nyline)%nstabph,nofecond
!505    format(/a,a,10i4)
       if(allocated(mapfix)) then
!          write(*,*)'Deallocating mapfix'
          deallocate(mapfix)
       endif
       allocate(mapfix)
       allocate(mapfix%fixph(1))
       allocate(mapfix%stableph(1))
       allocate(mapfix%stablepham(1))
       allocate(mapfix%stable_phr(1))
       mapfix%nfixph=1
       mapfix%status=0
!.........................................................
! trying to impove mapping with two extensive axis variables
       nofecond=0
       do jax=1,maptop%number_ofaxis
          call locate_condition(axarr(jax)%seqz,pcond,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
! active condition means pcond%active=0 !!
          if(pcond%active.eq.0) activax=jax
          if(pcond%statev.gt.10) nofecond=nofecond+1
       enddo
! default value
       fixpham=zero
! fix nonzero phase select!!
!       if(nofecond.eq.2) then    ! trying to have nonzero fix phase ...
       if(nofecond.eq.17) then  ! skip trying to have nonzero fix phase
! ISOTHERMAL, two extensive axis variables
! test if we can have non-zero fix phase amount.  Calculate the equilibrium
! set positive amount both in mapfix and in phase_varres ...??
          write(*,47)activax,mapline%linefixph(1)%ixphase
47        format(/'*** Find_line nonzero fix phase',2i3,2(1pe12.4))
          mapfix%fixph=mapline%linefixph(1)
          mapfix%stablepham(1)=one
          mapfix%stableph(1)=mapnode%linehead(nyline)%stableph(1)
          mapfix%stable_phr(1)=mapnode%linehead(nyline)%stable_phr(1)
!          ceq2=>mapline%lineceq
          meqrec=>mapline%meqrec
          mapfix%nstabph=1
          write(*,54)mapfix%fixph%ixphase,mapfix%fixph%compset
54        format('SMP fix phase: ',2i5,1pe12.4)
          write(*,55)mapfix%nstabph,mapfix%stableph%ixphase,&
               mapfix%stableph%compset,mapfix%stablepham
55        format('SMP stable phase: ',i2,2i5,1pe12.4)
! meqrec is allocated inside calceq7. mode=1 use gridminimizer
! mode=0 no gridminimizer; mode=-1 no grid and no deallocation of phr
          mode=-1
          call calceq7(mode,meqrec,mapfix,mapline%lineceq)
          write(*,7)gx%bmperr,iadd,irem
7         format('Calculated equilibrium in map_findline',3i5)
          if(gx%bmperr.ne.0) then
             write(*,*)'SMP Failed calculate equilibrium try to adjust amounts'
             gx%bmperr=0
             iadd1=0; irem1=0
!             write(*,*)'SMP2A calling meq_sameset from map_findline 1'
             call meq_sameset(irem1,iadd1,mapx,mapline%meqrec,&
                  mapline%meqrec%phr,inmap,mapline%lineceq)
             write(*,*)'Check with meq_sameset: ',gx%bmperr,irem1,iadd1
             if(gx%bmperr.ne.0) then
                goto 1000
             elseif(iadd1.ne. 0 .or. irem1.ne.0) then
                write(*,*)'ignore nozero iadd or irem'
             endif
          endif
! try to change the amount of the fix phase by selecting a composition
! along the tieline with 30% of the fix phase
! Now we must change a condition ...
          call locate_condition(axarr(activax)%seqz,pcond,mapline%lineceq)
!          write(*,*)'SMP2A Located condition',activax
          svrrec=>pcond%statvar(1)
! NOTE: If we change fix/entered phase we must change axvals/axvals2
          svrtarget=svrrec
          svrtarget%argtyp=3
! calculate composition of entered phase
!          svrtarget%phase=meqrec%phr(sj)%iph
!          svrtarget%compset=meqrec%phr(sj)%ics
          svrtarget%phase=mapfix%stableph(1)%ixphase
          svrtarget%compset=mapfix%stableph(1)%compset
! This extracts the composition of the entered phase for first new line
! we must use a pointer in state_variable_val
          svr2=>svrtarget
          call state_variable_val(svr2,xstab,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
          svrtarget%phase=mapfix%fixph(1)%ixphase
          svrtarget%compset=mapfix%fixph(1)%compset
! This extracts the composition of the entered phase for first new line
! we must use a pointer in state_variable_val
          svr2=>svrtarget
          call state_variable_val(svr2,xfix,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
! set fix phase amount to 0.3 as we may find a third phase along the line ..
! but we must take into account how many moles of atoms in fix phase
!          natpermol=meqrec%phr(??fixphase)%curd%abnorm(1)
          iadd=mapfix%fixph(1)%ixphase
          write(*,*)'SMP Natpermol: ',iadd,meqrec%phr(iadd)%curd%abnorm(1)
          natpermol=one
          fixpham=0.3D0/natpermol
          xcorr=(one-fixpham)*xstab+fixpham*xfix
          write(*,71)fixpham,xstab,xfix,xcorr
71        format('Change: ',4(1pe16.8))
! first argument 1 means to extract the value, 0 means to set the value
          call condition_value(0,pcond,xcorr,mapline%lineceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Cannot set axis condition'
             gx%bmperr=4399; goto 1000
          endif
! Then call meq_sameset ignoring any new phases that tries to be stable
          iadd=0; irem=0
!          write(*,*)'SMP2A Calling meq_sameset from map_findline 2'
          call meq_sameset(irem,iadd,mapx,mapline%meqrec,&
               mapline%meqrec%phr,inmap,mapline%lineceq)
          if(gx%bmperr.ne.0) then
             gx%bmperr=0; goto 1000
          elseif(irem.gt.0 .or. irem.gt.0) then
             write(*,*)'ignoring new phases: ',irem,iadd
          endif
! change the amount of the fix phase
          allocate(mapfix%fixphamap(1))
          mapfix%fixphamap(1)=fixpham
! if hullerombuller true below then it will change fix and stable phase
          hullerombuller=.FALSE.
          mapfix%stablepham(1)=one-fixpham
          write(*,*)'find mapline conditions: '
          call list_conditions(kou,mapline%lineceq)
!          goto 1000
!..................................
       else
!-----------------------------------------------------------------------
! with a potential axis ?
! we should check that the not-fixed phase can vary composition ...
!          write(*,*)'SMP2A using code with nofecond.ne.17 !!'
          ip=mapnode%linehead(nyline)%stableph(1)%ixphase
! fixedcomposition is a logical funtion in gtp3F.F90
          if(fixedcomposition(ip)) then
             mapfix%fixph=mapnode%linehead(nyline)%stableph(1)
             hullerombuller=.TRUE.
!          write(*,*)'Selecting other phase as fix',mapfix%fixph%ixphase,&
!               mapfix%fixph%compset
          else
!          write(*,*)'Changing fix phase: ',mapline%linefixph(1)%ixphase,&
!               mapline%linefixph(1)%compset
             mapfix%fixph=mapline%linefixph(1)
             hullerombuller=.FALSE.
          endif
       endif
! create a heading text for the line
       phaseset=' '
       call get_phasetuple_name(mapfix%fixph(1),phaseset)
       if(gx%bmperr.ne.0) goto 1000
       ip=len_trim(phaseset)+4
       phaseset(ip-2:ip-2)='+'
! It seems to be diffcult to reset tjis variable ....
       repeatederr=0
!       write(*,*)'Fixed phase: ',mapfix%nfixph,&
!            mapfix%fixph%ixphase,mapfix%fixph%compset
       if(mapnode%linehead(nyline)%nstabph.gt.0) then
! this is stored only for "real" nodes
          mapfix%nstabph=1
          if(hullerombuller) then
             mapfix%stableph(1)=mapline%linefixph(1)
             mapfix%stable_phr(1)=mapline%linefix_phr(1)
          else
             mapfix%stableph(1)=mapnode%linehead(nyline)%stableph(1)
             mapfix%stable_phr(1)=mapnode%linehead(nyline)%stable_phr(1)
          endif
          call get_phasetuple_name(mapfix%stableph(1),phaseset(ip:))
          if(gx%bmperr.ne.0) goto 1000
! set positive amount both in mapfix and in phase_varres ...??
          mapfix%stablepham(1)=one-fixpham
          ip=len_trim(phaseset)
          if(ip.gt.1) then
             write(kou,516)mapline%lineid,&
                  mapline%lineceq%tpval(1),phaseset(1:ip)
516          format(/'New line: ',i3,' T=',F8.2,' with: ',a)
!             write(*,507)' *** Phase fix: ',mapfix%fixph(1)%ixphase,&
!                  mapfix%fixph(1)%compset,', entered: ',&
!                  mapfix%stableph(1)%ixphase,&
!                  mapfix%stableph(1)%compset,', old node: ',mapline%nodfixph
507          format(a,2i3,a,2i3,a,2i3)
          else
             write(kou,521)
521          format(/'Line with unknown phases, wow')
          endif
       else
          write(*,*)'No stable phase!! why??'
          write(*,*)'stable 4:  ',nyline,mapnode%linehead(nyline)%nstabph,&
               mapnode%linehead(nyline)%stableph(1)%ixphase,&
               mapnode%linehead(nyline)%stableph(1)%compset
          mapfix%nstabph=0
       endif
!       write(*,*)'SMP looking for segmentation fault'
!-------------------------------------------------------------
    else
! For STEP we should set a small positive amount of a new stable phase
!       if(mapnode%nodefix%phaseix.gt.0) then
       if(mapnode%nodefix%ixphase.gt.0) then
! If the fix phase at the node was disappearing the phase index is negative
!          write(*,*)'Add a small amount to the new stable phase: ',&
!               mapnode%nodefix%phase,mapnode%nodefix%compset
!          call get_phase_compset(abs(mapnode%nodefix%phaseix),&
          call get_phase_compset(abs(mapnode%nodefix%ixphase),&
               mapnode%nodefix%compset,lokph,lokcs)
          mapline%lineceq%phase_varres(lokcs)%amfu=1.0D-2
       endif
!
       phaseset=' '
       ip=1
       do jp=1,mapnode%linehead(1)%nstabph
          call get_phasetuple_name(mapnode%linehead(1)%stableph(jp),&
               phaseset(ip:))
          if(gx%bmperr.ne.0) goto 1000
          ip=len_trim(phaseset)+2
       enddo
       if(ip.gt.1) then
! just to get current value of axis condition
          call locate_condition(axarr(1)%seqz,pcond,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
          call condition_value(1,pcond,finc,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
          write(kou,522)mapline%lineid,finc,phaseset(1:ip)
522       format(/'Line ',i3,' from ',1pe14.6,' with: ',a)
       else
          write(*,*)'Line with unkonwn stable phases: ',&
               mapnode%linehead(1)%nstabph
       endif
!       write(*,*)'SMP is mapfix allocated? ',allocated(mapfix)
!       if(.not.allocated(mapfix)) then
! for STEP calculations mapfix was normally not allocated but I need the status
! but instead of adding this set a bit in the meqrec record after first
! call to calceq7
!          allocate(mapfix)
!          mapfix%nfixph=0
!          mapfix%status=0
!          if(btest(mapnode%status,STEPINVARIANT)) then
!             write(*,*)'SMP invarant step node',mapnode%status
!             mapfix%status=ibset(mapfix%status,STEPINVARIANT)
!          endif
!       endif
    endif
1000 continue
    return
  end subroutine map_findline

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_findline_old
!\begin{verbatim}
  subroutine map_findline_old(maptop,axarr,mapfix,mapline)
! must be THREADPROTECTED
! Searches all node records from maptop for a map_line record to be calculated
! ?? already been found and if so eliminate a line record ??
! maptop map node record
! axarr array with axis records
! mapfix returned fixph record with phases to be ste as fixed for this line
! mapline returned mapline record for line to be calculated
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(map_axis), dimension(*) :: axarr
    type(map_fixph), allocatable :: mapfix
! memory leak as mapfix is allocated below ...
!    type(map_fixph), pointer :: mapfix
!\end{verbatim}
    type(map_node), pointer :: mapnode
    type(gtp_condition), pointer :: pcond
    type(meq_setup), pointer :: meqrec
    type(gtp_equilibrium_data), pointer :: ceq2
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(gtp_state_variable), target :: svrtarget
    integer nyline,jp,seqy,iph,ics,lokph,lokcs,ip,mapx
    integer mode,nofecond,jax,activax,irem1,iadd1,irem,iadd
    integer, parameter :: inmap=1
    double precision finc,fixpham,xstab,xfix,xcorr,natpermol
    character eqname*24
! sometimes there are many phases with long names ...
    character phaseset*512
    logical hullerombuller
!
    mapnode=>maptop
! for the moment skip this for tie-lines not in the plane
!    write(*,*)'In map_findline: ',mapnode%tieline_inplane
100 continue
!       write(*,*)'Looking a lines exiting from all nodes:'
!       write(*,*)'mapnode index: ',mapnode%seqx
       if(.not.allocated(mapnode%linehead)) then
          write(*,*)'ERROR found mapnode without exits'
          nullify(mapline)
          goto 1000
       endif
       do nyline=1,mapnode%lines
! if done is >=0 then this is a line to be calculated
          if(mapnode%linehead(nyline)%done.ge.0) then
             mapline=>mapnode%linehead(nyline)
             mapline%done=-1
             goto 500
          endif
       enddo
!       write(*,*)'findline 1: ',mapnode%seqx
       mapnode=>mapnode%next
!       write(*,*)'findline 2: ',mapnode%seqx
       if(.not.associated(mapnode,maptop)) goto 100
! no more lines to calculate
    if(ocv()) write(*,*)'nullifying mapline as no more lines to calculate'
    nullify(mapline); goto 1000
! jump here if we found a nyline
!--------------------------------------------------------------------
500 continue
! we must copy the equilibrium record ceq to the line record
    goto 503
!----------------------------------------------------------------------
! code below moved to ....
    if(ocv()) write(*,*)'We found a line from node: ',mapnode%seqx
    eqname='_MAPLINE_'
    jp=10
    seqy=maptop%seqy+1
    call wriint(eqname,jp,seqy)
    call copy_equilibrium(mapline%lineceq,eqname,mapnode%nodeceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error creating equilibrium: ',eqname
       goto 1000
    endif
    maptop%seqy=seqy
    mapline%lineid=seqy
    mapline%nodfixph=0
! mapline%more is positive while line is calculated, 0 means at axis limit
    mapline%more=1
! end code moved
!---------------------------------------------------------------------
503 continue
!    write(*,*)'At label 503',mapline%firstinc
    if(mapline%firstinc.ne.zero) then
! update the axis condition if mapline%firstinc is nonzero
       jp=abs(mapline%axandir) 
       call locate_condition(axarr(jp)%seqz,pcond,mapline%lineceq)
       if(gx%bmperr.ne.0) goto 1000
! Wow, here is the problem !!!
       pcond%prescribed=pcond%prescribed+mapline%firstinc
       finc=mapline%firstinc
       if(ocv()) write(*,501)'Selecting line and condition: ',&
            seqy,pcond%prescribed,finc
501    format(a,i3,2(1pe14.6))
!       if(pcond%active.ne.0) then
!          write(*,*)'Error: axis condition not active!'
!          pcond%active=0
!       endif
    endif
! check that correct axis condition is active, maybe I have not made sure
! that the map_line records are independent ...
!    write(*,*)'axies ',maptop%number_ofaxis
    do jp=1,maptop%number_ofaxis
       call locate_condition(axarr(jp)%seqz,pcond,mapline%lineceq)
!       write(*,*)'condition located',jp,axarr(jp)%seqz
       if(gx%bmperr.ne.0) goto 1000
       if(pcond%active.ne.0) then
          if(jp.eq.abs(mapline%axandir)) then
!             write(*,*)'Setting axis condition active!',jp,mapline%axandir
             pcond%active=0
!             write(*,*)'map_findline: ',jp,pcond%prescribed
          endif
       else
          if(jp.ne.abs(mapline%axandir)) then
!             write(*,*)'Setting axis condition NOT active!',jp,mapline%axandir
             pcond%active=1
          endif
       endif
    enddo
! we may have a set of stable phases in mapnode%stable_phases, maybe they
! should be set, at least when mapping. 
! a meqrec record will be created by calceq7 at the first calculation
! for mapping set values in mapfix about which phases that should be fix
! or stable when calling calceq7, at present ingnore that
!-------------------------------------------------------------
!    write(*,*)'tielines: ',maptop%tieline_inplane
    if(maptop%tieline_inplane.lt.0) then
! ISOPLETH
       if(allocated(mapfix)) then
          deallocate(mapfix)
       endif
       allocate(mapfix)
       mapfix%status=0
! with only 2 axis we have just 1 fix phase for mapping
       allocate(mapfix%fixph(1))
       mapfix%nfixph=1
       mapfix%fixph=mapline%linefixph(1)
! we can have several stable phases when no tie-lines in plane
       ip=mapline%nstabph
       allocate(mapfix%stableph(ip))
       allocate(mapfix%stablepham(ip))
       allocate(mapfix%stable_phr(ip))
!       write(*,*)'Findline: Tie-lines not in plane: ',nyline,ip
! create a heading text for the line
       phaseset=' '
       call get_phasetuple_name(mapfix%fixph(1),phaseset)
       if(gx%bmperr.ne.0) goto 1000
       ip=len_trim(phaseset)+4
       phaseset(ip-1:ip-2)='+'
       if(mapnode%linehead(nyline)%nstabph.le.0) then
          write(*,*)'No stable phases for a line'
          write(*,*)'Error 14:  ',nyline,mapnode%linehead(nyline)%nstabph,&
               mapnode%linehead(nyline)%stableph(1)%ixphase,&
               mapnode%linehead(nyline)%stableph(1)%compset
          mapfix%nstabph=0
          gx%bmperr=4242; goto 1000
       endif
       mapfix%nstabph=mapnode%linehead(nyline)%nstabph
!       write(*,*)'Findline: stable phases: ',mapfix%nstabph
       do jp=1,mapfix%nstabph
! this is stored only for "real" nodes
          mapfix%stableph(jp)=mapnode%linehead(nyline)%stableph(jp)
          mapfix%stable_phr(jp)=mapnode%linehead(nyline)%stable_phr(jp)
          call get_phasetuple_name(mapfix%stableph(jp),phaseset(ip:))
          if(gx%bmperr.ne.0) goto 1000
! this values hould perhaps be in linehead??
          mapfix%stablepham(jp)=one
!          call get_phase_compset(mapfix%stableph(1)%phase,&
!               mapfix%stableph(1)%compset,lokph,lokcs)
!          mapline%lineceq%phase_varres(lokcs)%amfu=one
          ip=len_trim(phaseset)+2
       enddo
       write(kou,520)mapline%lineid,mapline%lineceq%tpval(1),phaseset(1:ip)
520    format(/'Line ',i3,' T=',F8.2,' with: ',a)
!-------------------------------------------------------------
    elseif(maptop%tieline_inplane.gt.0) then
! TIE-LINES IN PLANE, NOTE: meqrec not allocated!!
!       if(mapnode%nodefix%phase.gt.0) then
!          write(*,*)'We have to set the fix phase along the line: ',&
!               mapnode%nodefix%phase,mapnode%nodefix%compset
! phr no longer allocated ...
!          iph=mapnode%meqrec%phr(mapnode%fixph)%iph
!          ics=mapnode%meqrec%phr(mapnode%fixph)%ics
!          call get_phase_compset(iph,ics,lokph,lokcs)
!          write(*,*)'Setting phase status as fixed, phase_varres: ',lokcs
!          mapline%lineceq%phase_varres(lokcs)%status2=&
!               ibset(mapline%lineceq%phase_varres(lokcs)%status2,PHFIXED)
!       endif
! mapline here should be identical to mapnode%linehead(nyline)
!       if(ocv()) write(*,505)'In findline: add phase set for',&
!       write(*,505)'In findline: add phase set for',&
!            ' tie-lines in plane, node:',&
!            mapnode%seqx,nyline,mapnode%linehead(nyline)%nstabph
!505    format(a,a,10i4)
       if(allocated(mapfix)) then
!          write(*,*)'Deallocating mapfix'
          deallocate(mapfix)
       endif
       allocate(mapfix)
       allocate(mapfix%fixph(1))
       allocate(mapfix%stableph(1))
       allocate(mapfix%stablepham(1))
       allocate(mapfix%stable_phr(1))
       mapfix%nfixph=1
       mapfix%status=0
!.........................................................
! trying to impove mapping with two extensive axis variables
       nofecond=0
       do jax=1,maptop%number_ofaxis
          call locate_condition(axarr(jax)%seqz,pcond,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
! active condition means pcond%active=0 !!
          if(pcond%active.eq.0) activax=jax
          if(pcond%statev.gt.10) nofecond=nofecond+1
       enddo
! default value
       fixpham=zero
! fix nonzero phase select!!
!       if(nofecond.eq.2) then    ! trying to have nonzero fix phase ...
       if(nofecond.eq.17) then  ! skip trying to have nonzero fix phase
! ISOTHERMAL, two extensive axis variables
! test if we can have non-zero fix phase amount.  Calculate the equilibrium
! set positive amount both in mapfix and in phase_varres ...??
          write(*,47)activax,mapline%linefixph(1)%ixphase
47        format(/'*** Find_line nonzero fix phase',2i3,2(1pe12.4))
          mapfix%fixph=mapline%linefixph(1)
          mapfix%stablepham(1)=one
          mapfix%stableph(1)=mapnode%linehead(nyline)%stableph(1)
          mapfix%stable_phr(1)=mapnode%linehead(nyline)%stable_phr(1)
!          ceq2=>mapline%lineceq
          meqrec=>mapline%meqrec
          mapfix%nstabph=1
          write(*,54)mapfix%fixph%ixphase,mapfix%fixph%compset
54        format('SMP fix phase: ',2i5,1pe12.4)
          write(*,55)mapfix%nstabph,mapfix%stableph%ixphase,&
               mapfix%stableph%compset,mapfix%stablepham
55        format('SMP stable phase: ',i2,2i5,1pe12.4)
! meqrec is allocated inside calceq7. mode=1 use gridminimizer
! mode=0 no gridminimizer; mode=-1 no grid and no deallocation of phr
          mode=-1
          call calceq7(mode,meqrec,mapfix,mapline%lineceq)
          write(*,7)gx%bmperr,iadd,irem
7         format('Calculated equilibrium in map_findline',3i5)
          if(gx%bmperr.ne.0) then
             write(*,*)'SMP Failed calculate equilibrium try to adjust amounts'
             gx%bmperr=0
             iadd1=0; irem1=0
!             write(*,*)'SMP2A Calling meq_sameset from map_findline_old 1'
             call meq_sameset(irem1,iadd1,mapx,mapline%meqrec,&
                  mapline%meqrec%phr,inmap,mapline%lineceq)
             write(*,*)'Check with meq_sameset: ',gx%bmperr,irem1,iadd1
             if(gx%bmperr.ne.0) then
                goto 1000
             elseif(iadd1.ne. 0 .or. irem1.ne.0) then
                write(*,*)'ignore nozero iadd or irem'
             endif
          endif
! try to change the amount of the fix phase by selecting a composition
! along the tieline with 30% of the fix phase
! Now we must change a condition ...
          call locate_condition(axarr(activax)%seqz,pcond,mapline%lineceq)
          write(*,*)'Located condition',activax
          svrrec=>pcond%statvar(1)
! NOTE: If we change fix/entered phase we must change axvals/axvals2
          svrtarget=svrrec
          svrtarget%argtyp=3
! calculate composition of entered phase
!          svrtarget%phase=meqrec%phr(sj)%iph
!          svrtarget%compset=meqrec%phr(sj)%ics
          svrtarget%phase=mapfix%stableph(1)%ixphase
          svrtarget%compset=mapfix%stableph(1)%compset
! This extracts the composition of the entered phase for first new line
! we must use a pointer in state_variable_val
          svr2=>svrtarget
          call state_variable_val(svr2,xstab,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
          svrtarget%phase=mapfix%fixph(1)%ixphase
          svrtarget%compset=mapfix%fixph(1)%compset
! This extracts the composition of the entered phase for first new line
! we must use a pointer in state_variable_val
          svr2=>svrtarget
          call state_variable_val(svr2,xfix,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
! set fix phase amount to 0.3 as we may find a third phase along the line ..
! but we must take into account how many moles of atoms in fix phase
!          natpermol=meqrec%phr(??fixphase)%curd%abnorm(1)
          iadd=mapfix%fixph(1)%ixphase
          write(*,*)'SMP Natpermol: ',iadd,meqrec%phr(iadd)%curd%abnorm(1)
          natpermol=one
          fixpham=0.3D0/natpermol
          xcorr=(one-fixpham)*xstab+fixpham*xfix
          write(*,71)fixpham,xstab,xfix,xcorr
71        format('Change: ',4(1pe16.8))
! first argument 1 means to extract the value, 0 means to set the value
          call condition_value(0,pcond,xcorr,mapline%lineceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Cannot set axis condition'
             gx%bmperr=4399; goto 1000
          endif
! Then call meq_sameset ignoring any new phases that tries to be stable
          iadd=0; irem=0
!          write(*,*)'SMP2A Calling meq_sameset from map_findline_old 2'
          call meq_sameset(irem,iadd,mapx,mapline%meqrec,&
               mapline%meqrec%phr,inmap,mapline%lineceq)
          if(gx%bmperr.ne.0) then
             gx%bmperr=0; goto 1000
          elseif(irem.gt.0 .or. irem.gt.0) then
             write(*,*)'ignorin new phases: ',irem,iadd
          endif
! change the amount of the fix phase
          allocate(mapfix%fixphamap(1))
          mapfix%fixphamap(1)=fixpham
! if hullerombuller true below then it will change fix and stable phase
          hullerombuller=.FALSE.
          mapfix%stablepham(1)=one-fixpham
          write(*,*)'find mapline conditions: '
          call list_conditions(kou,mapline%lineceq)
!          goto 1000
!..................................
       else
!-----------------------------------------------------------------------
! with potential axis ?
! we should check that the not-fixed phase can vary composition ...
!          write(*,*)'SMP2A using old code 2 for nofecond.ne.17'
          ip=mapnode%linehead(nyline)%stableph(1)%ixphase
! fixedcomposition is a logical funtion in gtp3F.F90
          if(fixedcomposition(ip)) then
             mapfix%fixph=mapnode%linehead(nyline)%stableph(1)
             hullerombuller=.TRUE.
!          write(*,*)'Selecting other phase as fix',mapfix%fixph%ixphase,&
!               mapfix%fixph%compset
          else
!          write(*,*)'Changing fix phase: ',mapline%linefixph(1)%ixphase,&
!               mapline%linefixph(1)%compset
             mapfix%fixph=mapline%linefixph(1)
             hullerombuller=.FALSE.
          endif
       endif
! create a heading text for the line
       phaseset=' '
       call get_phasetuple_name(mapfix%fixph(1),phaseset)
       if(gx%bmperr.ne.0) goto 1000
       ip=len_trim(phaseset)+4
       phaseset(ip-2:ip-2)='+'
!       write(*,*)'Fixed phase: ',mapfix%nfixph,&
!            mapfix%fixph%ixphase,mapfix%fixph%compset
       if(mapnode%linehead(nyline)%nstabph.gt.0) then
! this is stored only for "real" nodes
          mapfix%nstabph=1
          if(hullerombuller) then
             mapfix%stableph(1)=mapline%linefixph(1)
             mapfix%stable_phr(1)=mapline%linefix_phr(1)
          else
             mapfix%stableph(1)=mapnode%linehead(nyline)%stableph(1)
             mapfix%stable_phr(1)=mapnode%linehead(nyline)%stable_phr(1)
          endif
          call get_phasetuple_name(mapfix%stableph(1),phaseset(ip:))
          if(gx%bmperr.ne.0) goto 1000
! set positive amount both in mapfix and in phase_varres ...??
          mapfix%stablepham(1)=one-fixpham
          ip=len_trim(phaseset)
          if(ip.gt.1) then
             write(kou,516)mapline%lineid,&
                  mapline%lineceq%tpval(1),phaseset(1:ip)
516          format(/'Line: ',i3,' T=',F8.2,' with: ',a)
!             write(*,507)' *** Phase fix: ',mapfix%fixph(1)%ixphase,&
!                  mapfix%fixph(1)%compset,', entered: ',&
!                  mapfix%stableph(1)%ixphase,&
!                  mapfix%stableph(1)%compset,', old node: ',mapline%nodfixph
507          format(a,2i3,a,2i3,a,2i3)
          else
             write(kou,521)
521          format(/'Line with unknown phases, wow')
          endif
       else
          write(*,*)'No stable phase!! why??'
          write(*,*)'stable 4:  ',nyline,mapnode%linehead(nyline)%nstabph,&
               mapnode%linehead(nyline)%stableph(1)%ixphase,&
               mapnode%linehead(nyline)%stableph(1)%compset
          mapfix%nstabph=0
       endif
!-------------------------------------------------------------
    else
! For STEP we should set a small positive amount of a new stable phase
!       if(mapnode%nodefix%phaseix.gt.0) then
       if(mapnode%nodefix%ixphase.gt.0) then
! If the fix phase at the node was disappearing the phase index is negative
!          write(*,*)'Add a small amount to the new stable phase: ',&
!               mapnode%nodefix%phase,mapnode%nodefix%compset
!          call get_phase_compset(abs(mapnode%nodefix%phaseix),&
          call get_phase_compset(abs(mapnode%nodefix%ixphase),&
               mapnode%nodefix%compset,lokph,lokcs)
          mapline%lineceq%phase_varres(lokcs)%amfu=1.0D-2
       endif
!
       phaseset=' '
       ip=1
       do jp=1,mapnode%linehead(1)%nstabph
          call get_phasetuple_name(mapnode%linehead(1)%stableph(jp),&
               phaseset(ip:))
          if(gx%bmperr.ne.0) goto 1000
          ip=len_trim(phaseset)+2
       enddo
       if(ip.gt.1) then
! just to get current value of axis condition
          call locate_condition(axarr(1)%seqz,pcond,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
          call condition_value(1,pcond,finc,mapline%lineceq)
          if(gx%bmperr.ne.0) goto 1000
          write(kou,522)mapline%lineid,finc,phaseset(1:ip)
522       format(/'Line ',i3,' from ',1pe14.6,' with: ',a)
       else
          write(*,*)'Line with unkonwn stable phases: ',&
               mapnode%linehead(1)%nstabph
       endif
! for step calculation mapfix is not needed
!       nullify(mapfix)
    endif
1000 continue
    return
  end subroutine map_findline_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_saveceq
!\begin{verbatim}
  subroutine create_saveceq(ceqres,size)
! creates an array of equilibrium records to save calculated lines for step
! and map.  This can be very big
    type(map_ceqresults), pointer :: ceqres
    integer size
!\end{verbatim}
!    write(*,*)'In create saveceq',size
    integer errall
    allocate(ceqres)
    ceqres%size=size
    ceqres%free=1
    allocate(ceqres%savedceq(size),stat=errall)
    if(errall.ne.0) then
       write(*,*)'SMP2A Allocation error 1: ',errall
       gx%bmperr=4370; goto 1000
    endif
1000 continue
    return
  end subroutine create_saveceq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine delete_mapresults
!\begin{verbatim}
  subroutine delete_mapresults(maptop)
! delete all saved results created by step or map
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    type(map_ceqresults), pointer :: saveceq
    TYPE(map_node), pointer :: current,nexttop,mapnode,delnode
    TYPE(map_line), pointer :: linehead
    TYPE(gtp_equilibrium_data), pointer :: ceq
    integer ieq,jj
!    integer place,lastused
!
    if(.not.associated(maptop)) then
       write(*,*)'No step or map results to delete'
       goto 1000
    endif
!    write(*,*)'smp2A in delete_mapresults'
    current=>maptop
!    deloop: do while(associated(current))
!       write(*,*)'smp2A maybe no saveceq?',associated(current%saveceq)
!       if(associated(current%saveceq)) &
!            write(*,*)'Saved equilibria:',current%saveceq%free-1
!       current=>current%plotlink
!    enddo deloop
!    write(*,*)'All equilibria saved in mapnodes listed'
! all mapnodes has a pointer to first where the saveceq is allocated
    current=>maptop
    do while(associated(current))
!       write(*,*)'smp2a current associated'
       if(associated(current%saveceq)) then
          if(allocated(current%saveceq%savedceq)) then
             write(*,*)'SMP: deleting saved step/map line equilibria: ',&
                  current%saveceq%free-1
             deallocate(current%saveceq%savedceq)
          endif
       endif
! adding this write avoided a segmentation fault ... no longer ...
!       write(*,*)'SMP: are there more mapnode records?',&
!            associated(current%plotlink),associated(current%next)
       nexttop=>current%plotlink
       mapnode=>current%next
       do while(.not.associated(mapnode,current))
!          write(*,*)'SMP: cleaning up more',mapnode%lines
          if(allocated(mapnode%linehead)) then
!             write(*,*)'SMP: cleaning maplines: ',size(mapnode%linehead)
             do jj=1,mapnode%lines
! should these be deallocated explicitly??
                linehead=>mapnode%linehead(jj)
                if(allocated(linehead%axvals)) deallocate(linehead%axvals)
                if(allocated(linehead%axvals2)) deallocate(linehead%axvals2)
                if(allocated(linehead%axvalx)) deallocate(linehead%axvalx)
             enddo
             deallocate(mapnode%linehead)
          endif
          delnode=>mapnode
          mapnode=>mapnode%next
          deallocate(delnode)
       enddo
       delnode=>current
       current=>nexttop
! deallocate the last mapnode
       if(associated(current)) deallocate(delnode)
    enddo
    write(*,*)'Deleting _MAPx equilibria'
    ceq=>firsteq
    call delete_equilibria('_MAP*',ceq)
1000 continue
    return
  end subroutine delete_mapresults

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function tieline_inplane
!\begin{verbatim}
  integer function tieline_inplane(nax,axarr,ceq)
! returns -1 if tielines are not in the plane (isopleth)
!          0 for step calculations (nax=1)
!          1 if tielines in the plane (binary T-X, ternary isopleths
!          set if more than one extensive variable is not axis variables
! nax number of axis
! axarr array with axis records
    integer nax
    type(map_axis), dimension(nax) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    type(gtp_condition), pointer :: lastcond,pcond
    integer noc,inplane,nexv,iax
!
    inplane=0
    if(nax.eq.1) goto 1000
    lastcond=>ceq%lastcondition
    if(.not.associated(lastcond)) then
       write(*,*)'Whops, mapping with no conditions?'
       gx%bmperr=4243; goto 1000
    endif
    nexv=0
    pcond=>lastcond
100 continue
       pcond=>pcond%next
       if(pcond%statev.gt.9) then
! statev>10 means extensive variable, maximum one not axis variable 
! For example binary T-X has extra conditions, P,N; ternary X-X isoterm T,P,N
! A fix chemical potential OK, a fix phase is the same as activity condition
          if(pcond%active.eq.0) then
! active=0 means it is an active condition
             do iax=1,nax
                if(axarr(iax)%seqz.eq.pcond%seqz) goto 200
             enddo
! we have a condition on an extensive variable that is not an axis
             nexv=nexv+1
200          continue
          endif
       endif
       if(.not.associated(pcond,lastcond)) goto 100
       inplane=-1
       if(nexv.le.1) inplane=1
1000 continue
    tieline_inplane=inplane
!    if(ocv()) write(*,*)'tie-line in plane return: ',inplane
    return
  end function tieline_inplane

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_map_equilibrium
!\begin{verbatim}
  subroutine list_map_equilibrium(maptop,mapline,axarr,xxx,typ)
! output of all relevant infor for a failed equilibrium calculation
! maptop map node record
! mapline current line record
! axarr array with axis records
! xxx current active axis value that caused problems to calculate
! typ indicates the type of problem
    integer typ
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(map_axis), dimension(*) :: axarr
    double precision xxx
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    integer jj,nph,lokph,lokcs,fixph,fixcs
    character name*24
    double precision yyy
! list current conditions (indicate active axis variable)
! list stable phases
!    write(*,*)'SMP map problems: ',typ,gx%bmperr,mapline%nodfixph
    ceq=>mapline%lineceq
    call list_conditions(kou,ceq)
! There is only one fix phase at all mapping at present !!
    jj=1
    fixph=mapline%meqrec%fixph(1,jj) 
    fixcs=mapline%meqrec%fixph(2,jj)
!
!    nph=noofphasetuples()
    nph=nooftup()
    write(*,66,advance='no')
66  format('Phases: ')
    do jj=1,nph
       lokcs=phasetuple(jj)%lokvares
       if(ceq%phase_varres(lokcs)%phstate.eq.PHENTSTAB) then
          yyy=ceq%phase_varres(lokcs)%amfu
          call get_phase_name(phasetuple(jj)%ixphase,phasetuple(jj)%compset,&
               name)
          if(phasetuple(jj)%ixphase.eq.fixph .and. &
               phasetuple(jj)%compset.eq.fixcs) then
             write(*,67,advance='no')'*'//trim(name)//'=',yyy
          else
             write(*,67,advance='no')trim(name)//'=',yyy
          endif
67        format(a,F4.1,1x)
       elseif(ceq%phase_varres(lokcs)%phstate.eq.PHFIXED) then
! Ahhh, the fix phase is not set as condition in ceq!!
          call get_phase_name(phasetuple(jj)%ixphase,phasetuple(jj)%compset,&
               name)
          write(*,67,advance='no')'*'//trim(name)//' '
       endif
    enddo
    write(*,77)'SMP: ',fixph,fixcs,mapline%axandir,xxx
77  format(/,a,3i3,1pe14.6)
! try for the AL-Cr-Ni case ... tuple 16, FCC_L12#2, should not be stable ...
!    lokcs=phasetuple(16)%lokvares
!    ceq%phase_varres(lokcs)%phstate=PHENTERED
! 15 is FCC_L12 is fix with 1 mole! try changing amounts
!    ceq%phase_varres(14)%amfu=one
!    ceq%phase_varres(15)%amfu=zero
! 15 is FCC_L12 is fix with 1 mole! try changing fix phase to 14 BCC
! I think meqrec is deallocated after this we have to change somewhere else
!    jj=1
!    mapline%meqrec%fixph(1,jj)=14
!    mapline%meqrec%fixph(2,jj)=1
!
1000 continue
    return
  end subroutine list_map_equilibrium

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_problems
!\begin{verbatim}
  subroutine map_problems(maptop,mapline,axarr,xxx,typ)
! jump here for different problems
! maptop map node record
! mapline current line record
! axarr array with axis records
! xxx current active axis value that caused problems to calculate
! typ indicates the type of problem
    integer typ
    type(map_node), pointer :: maptop
    type(map_line), pointer :: mapline
    type(map_axis), dimension(*) :: axarr
    double precision xxx
!\end{verbatim}
    character ch1*1
    integer oldaxis
    double precision yyy
! skip debug output
!    write(*,7)'In map_problem: ',typ,mapline%problems,mapline%lasterr,&
!         mapline%axandir,mapline%nodfixph,maptop%number_ofaxis,xxx
7   format(a,6i4,6(1pe14.6))
! we can list the current conditions here, 
! note fix phases for mapping not included as condition here!!
!    write(*,*)'Map problems 1'
!    call list_conditions(kou,mapline%lineceq)
!    call list_short_results(kou,mapline%lineceq)
!    read(*,10)ch1
10  format(a)
! for debugging:
!    call list_map_equilibrium(maptop,mapline,axarr,xxx,typ)
!    if(mapline%problems.gt.5) then
    if(mapline%problems.gt.2) then
       if(mapline%nodfixph.gt.0) then
!          call list_conditions(kou,mapline%lineceq)
!          if(gx%bmperr.ne.0) then
!             write(*,*)'Error listing conditions'
!             gx%bmperr=0
!          endif
          write(*,11)mapline%lineid,trim(mapline%lineceq%eqname)
11        format('SMP2A giving up on this line',i3,': ',a)
!          write(*,11)mapline%lineid,trim(mapline%lineceq%eqname),&
!               mapline%meqrec%phr(mapline%nodfixph)%iph,&
!               mapline%meqrec%phr(mapline%nodfixph)%ics
!11        format('I give up on this line',i3,2x,a,' with fix phase ',2i4)
!       else
!          write(*,11)mapline%nodfixph,mapline%lineceq%eqname,0,0
       endif
       gx%bmperr=4244; goto 1000
    endif
!    write(*,*)'Map problems 2'
!---------------------------------------------
! list current conditions
!    call list_conditions(kou,mapline%lineceq)
    if(maptop%number_ofaxis.eq.1) then
! for step only take smaller steps or calculate with grid minimizer
       if(typ.eq.1) then
! take a smaller step
! current axis condition value is xxx, mapline%firstinc is the step taken
          xxx=xxx-0.999*mapline%firstinc
       else
          write(*,*)'Unknown problem ',typ
          gx%bmperr=4245
       endif
       goto 1000
    endif
!=======================================================
! two or more axis
    select case(typ)
    case default
       write(*,*)'Unknown problem ',typ
       gx%bmperr=4245
!------------------------------------------------------
    case(1) ! error at first step, for map opposite direction
! current axis condition value is xxx, mapline%firstinc is the step taken
       yyy=xxx
!       write(*,*)'First increment: ',mapline%axandir,mapline%firstinc
       if(mapline%problems.eq.1) then
! first time here take the step in opposite direction
!          xxx=yyy-0.99D0*mapline%firstinc
!>>        xxx=yyy-1.01D0*mapline%firstinc        best tested value
          xxx=yyy-1.01D0*mapline%firstinc
          mapline%axandir=-mapline%axandir
       elseif(mapline%problems.eq.2) then
! second time take a small step in previous direction
!          xxx=yyy-0.02D0*mapline%firstinc
!>>        xxx=yyy+0.02D0*mapline%firstinc         best tested value
          xxx=yyy+0.02D0*mapline%firstinc
          mapline%axandir=-mapline%axandir
       elseif(mapline%problems.eq.3) then
! third time take small step in other axis
!          write(*,*)'Changing active axis'
! we must extract axis value, change condition etc. assume only 2 axis
!          oldaxis=mapline%axandir
!          mapline%axandir=3-mapline%axandir
!          call list_conditions(kou,mapline%lineceq)
       elseif(mapline%problems.eq.4) then
! fourth time take small step in opposite direction (of axis set with 3)
!          xxx=yyy-0.02D0*mapline%firstinc
!          mapline%axandir=-mapline%axandir
       endif
       mapline%axfact=1.0D-2
! the returned value xxx will be set as condition
!       call condition_value(0,pcond,xxx,ceq)
!       if(gx%bmperr.ne.0) goto 1000
!------------------------------------------------------
    end select
1000 continue
    return
  end subroutine map_problems

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_halfstep_bad
!\begin{verbatim}
  subroutine map_halfstep_bad(halfstep,type,axvalok,mapline,axarr,ceq)
! THIS MADE MANY MAP macro FAIL: 3,6,7,8,12,13 and finally crash ...
! Used when an error calculating a normal step or a node point
! take back the last sucessfully calculated axis value and take smaller step
! possibly one should also restore the ceq record.
! halfstep number of times halfstep has been called for this problem
! axvalok last cucessfully calculated value of active axis
! mapline line record
! axarr array with axis records
! ceq equilibrium record
    implicit none
    integer halfstep,type
    double precision axvalok
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(map_line), pointer :: mapline
    type(map_axis), dimension(*) :: axarr
!\end{verbatim}
    type(gtp_condition), pointer :: pcond
    double precision value
    double precision :: sfact=1.0D-2
    integer jax
!    write(*,*)'In map_halfstep_bad',halfstep
    if(halfstep.eq.1) then
       sfact=0.5D0
    else
       sfact=sfact*sfact
    endif
    halfstep=halfstep+1
    if(type.eq.1 .and. (axvalok.eq.zero .or. halfstep.ge.3)) then
!       write(*,*)'Two phases competing to appear/disappear',axvalok,halfstep
       gx%bmperr=4246
    else
! Previous axis value should be axvalok, find current
       jax=abs(mapline%axandir)
       call locate_condition(axarr(jax)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to get the value
       call condition_value(1,pcond,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'Current active axis value: ',value
! at first call remember the original axis value
       if(halfstep.eq.1) then
          if(ocv()) write(*,67)'First call to map_half, value:',value,axvalok
67        format(a,2(1pe14.6))
          mapline%evenvalue=value
       elseif(halfstep.gt.3) then
!          write(*,*)'SMP2A Tried halfstep 3 times, giving up'
          gx%bmperr=4368
       endif
       if(mapline%axfact.le.1.0D-6) then
! error initiallizing axfact ???
          write(*,*)'Too small value of mapline%axfact: ',mapline%axfact
          mapline%axfact=1.0D-3
       endif
! take a small step
       if(mapline%axandir.gt.0) then
          value=axvalok+sfact*mapline%axfact*axarr(jax)%axinc
       else
          value=axvalok-sfact*mapline%axfact*axarr(jax)%axinc
       endif
       write(*,97)'Halfstep axis value: ',mapline%axandir,value,axvalok,&
            mapline%axfact,axarr(jax)%axinc
97     format(a,i2,5(1pe14.6))
! first argument 0 means to set the value
       call condition_value(0,pcond,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
       if(ocv()) write(*,*)'Taking a small step, new axis value: ',jax,value
    endif
1000 continue
    return
  end subroutine map_halfstep_bad

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine map_halfstep
!\begin{verbatim}
  subroutine map_halfstep(halfstep,type,axvalok,mapline,axarr,ceq)
! Used when an error calculating a normal step or a node point
! take back the last sucessfully calculated axis value and take smaller step
! possibly one should also restore the ceq record.
! halfstep number of times halfstep has been called for this problem
! axvalok last cucessfully calculated value of active axis
! mapline line record
! axarr array with axis records
! ceq equilibrium record
    implicit none
    integer halfstep,type
    double precision axvalok
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(map_line), pointer :: mapline
    type(map_axis), dimension(*) :: axarr
!\end{verbatim}
    type(gtp_condition), pointer :: pcond
    double precision value
    double precision, parameter :: sfact=1.0D-2
    integer jax
    repeatederr=repeatederr+1
!    write(*,*)'In map_halfstep',halfstep,repeatederr
    halfstep=halfstep+1
    if(type.eq.1 .and. (axvalok.eq.zero .or. halfstep.ge.3)) then
!       write(*,*)'Two phases competing to appear/disappear',axvalok,halfstep
       gx%bmperr=4246
    else
! Previous axis value should be axvalok, find current
       jax=abs(mapline%axandir)
       call locate_condition(axarr(jax)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to get the value
       call condition_value(1,pcond,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'Current active axis value: ',value
! at first call remember the original axis value
       if(halfstep.eq.1) then
          if(ocv()) write(*,67)'First call to map_half, value:',value,axvalok
67        format(a,2(1pe14.6))
          mapline%evenvalue=value
       elseif(halfstep.gt.3) then
!          write(*,*)'SMP2A Tried halfstep 3 times, giving up'
          gx%bmperr=4368
       endif
       if(mapline%axfact.le.1.0D-6) then
! error initiallizing axfact ???
          write(*,*)'Too small value of mapline%axfact: ',mapline%axfact
          mapline%axfact=1.0D-3
       endif
! take a small step
       if(mapline%axandir.gt.0) then
          value=axvalok+1.0D-1*mapline%axfact*axarr(jax)%axinc
       else
          value=axvalok-1.0D-1*mapline%axfact*axarr(jax)%axinc
       endif
!       write(*,97)'Halfstep axis value: ',mapline%axandir,value,axvalok,&
!            mapline%axfact,axarr(jax)%axinc
97     format(a,i2,5(1pe14.6))
! first argument 0 means to set the value
       call condition_value(0,pcond,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
       if(ocv()) write(*,*)'Taking a small step, new axis value: ',jax,value
    endif
1000 continue
    return
  end subroutine map_halfstep

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine step_separate
!\begin{verbatim}
  subroutine step_separate(maptop,noofaxis,axarr,seqxyz,starteq)
! calculates for each phase separately along an axis (like G curves)
! There can not be any changes of the stable phase ...
! maptop map node record
! noofaxis must be 1
! axarr array of axis records
! seqxyz indices for map and line records
! starteq equilibrium record for starting
    implicit none
    integer noofaxis,seqxyz(*)
    type(map_axis), dimension(noofaxis) :: axarr
    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq
    integer ntup,itup,iph,ics,nystat,inactive(4),notop,seqy,mode,jk
    integer jj,seqz,iadd,irem,nv,saveq,lokcs,mapx
    type(gtp_phasetuple), dimension(:), allocatable :: entphcs
    integer, dimension(:), allocatable :: stsphcs
    type(map_line), pointer :: mapline
    type(map_fixph), allocatable :: mapfix
!    type(map_fixph), pointer :: mapfix
!    TYPE(map_node), pointer :: curtop
    type(meq_setup), pointer :: meqrec
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: svr
    type(meq_phase), pointer :: phr
    double precision val,xxx,yyy,axvalok
    logical firstline
!    integer, parameter :: maxsavedceq=2000
! decreased to 1800 as I sometimes run out of memeory
    integer, parameter :: maxsavedceq=1800
    character name*24
! turns off convergence control for T
    integer, parameter :: inmap=1
!
!    write(*,*)'In step_separate'
    if(noofaxis.ne.1) then
       write(kou,*)'Step separate only with one axis variable'
       goto 1000
    endif
! this subroutine returnes the total number of phase and composition sets
!    call sumofphcs(ntup,ceq)
!    ntup=totalphcs(starteq)
!    ntup=nonsusphcs(starteq)
    ntup=nooftup()
    allocate(entphcs(ntup))
    allocate(stsphcs(ntup))
    itup=0
    ceq=>starteq
! collect all current phase status and set all phases suspended
    nystat=-3
    val=zero
    do iph=1,noph()
       do ics=1,noofcs(iph)
          itup=itup+1
!          entphcs(itup)%phaseix=iph
          entphcs(itup)%ixphase=iph
          entphcs(itup)%compset=ics
          stsphcs(itup)=test_phase_status(iph,ics,val,ceq)
!          write(*,*)'step-sep ',iph,noofcs(iph),ics,itup,stsphcs(itup)
          if(gx%bmperr.ne.0) goto 1000
! phase status -1, 0 and 1 are all saved as 0
          if(stsphcs(itup).ge.-1 .and. stsphcs(itup).le.1) stsphcs(itup)=0
! do not change status of dormant phases ...
          if(stsphcs(itup).ne.-2) then
             call change_phase_status(iph,ics,nystat,val,ceq)
             if(gx%bmperr.ne.0) goto 1000
          endif
       enddo
    enddo
!    write(*,'(a,10i3)')'Suspended all phases',stsphcs
! indicator if maptop allocated
!    nullify(curtop)
    notop=0
! loop through all phases with stsphcs less than 3
!    nystat=0
!============================================================
    phaseloop: do itup=1,ntup
! nystat:-4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fix
!       write(*,*)'SS Phase ',itup,' status ',stsphcs(itup),ntup
       if(stsphcs(itup).gt.-2) then
! set default constitution, if none specified in the middle
!          write(*,*)'loop for phase phase ',itup,' stable'
          iph=entphcs(itup)%ixphase
          call set_default_constitution(entphcs(itup)%ixphase,&
               entphcs(itup)%compset,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Failed setting default constitution'
             goto 500
          endif
! set phase as entered
!          write(*,*)'Set phase stable ',itup,entphcs(itup)%phase
          call change_phase_status(entphcs(itup)%ixphase,&
               entphcs(itup)%compset,1,one,ceq)
!               entphcs(itup)%compset,0,one,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Failed setting phase entered',gx%bmperr
             goto 500
          endif
! debug listing of phase constitution to check
!          call list_phase_model(entphcs(itup)%ixphase,entphcs(itup)%compset,&
!               kou,ceq)
! here we should set the condition for overall composition to that of the phase
! Extract the current value of the axis state variable items using pcond
!          write(*,*)'Extracting axis condition value '
          seqz=axarr(1)%seqz
!          write(*,*)'Locating condition ',seqz
          call locate_condition(seqz,pcond,ceq)
          if(gx%bmperr.ne.0) goto 500
! if condition is a composition set it to be the current value with the
! default composition of the phase, 17 is mole fraction
          svr=>pcond%statvar(1)
          if(svr%statevarid.eq.17) then
! this call calculates the value of the axis condition with default composition
             call state_variable_val(svr,val,ceq)
             if(gx%bmperr.ne.0) goto 500
             call get_phasetuple_name(entphcs(itup),name)
! axis variable is composition, skip hases with no variance
             call get_phase_variance(entphcs(itup)%ixphase,nv)
             if(nv.eq.0) then
                write(*,71)name(1:len_trim(name)),val
71              format(/'Ignoring phase with fixed composition: ',a,F10.6)
!----------------
                lokcs=phasetuple(iph)%lokvares
                write(*,*)'indices: ',iph,phasetuple(iph)%ixphase,lokcs
                goto 500
! handle stoichiometric phases in step_separate ....
! we need to initiate a line with just one point
! special call to map_startpoint/map_findline for just one point
!                inactive=0
!                call map_startpoint(maptop,noofaxis,axarr,seqxyz,inactive,ceq)
!                if(gx%bmperr.ne.0) goto 500
!                call map_findline(maptop,axarr,mapfix,mapline)
!                if(gx%bmperr.ne.0) goto 500
!                ceq=>mapline%lineceq
!                meqrec=>mapline%meqrec
! this call gives error meqrec allready allocated
!                ceq=>??
!                call calceq7(mode,meqrec,mapfix,ceq)
!                if(gx%bmperr.ne.0) then
!                   write(*,*)'Error calculating stoichiometric phase',gx%bmperr
!                endif
! store the value of G
!                call map_store(mapline,axarr,maptop%number_ofaxis,&
!                     maptop%saveceq)
!                if(gx%bmperr.ne.0) then
!                   write(*,*)'Error storing equilibrium',gx%bmperr
!                   goto 900
!                endif
! change the calculated value of G by adding 1.0D4 and store
!                mapline%lineceq%phase_varres(lokcs)%gval(1,1)=&
!                     mapline%lineceq%phase_varres(lokcs)%gval(1,1)+1.0D3
!                call map_lineend(mapline,val,ceq)
!                goto 500
!----------------
             endif
!             if(ocv()) write(*,73)name(1:len_trim(name)),val
! check if val is within axis limits
             if(val.lt.axarr(1)%axmin .or. val.gt.axarr(1)%axmax) then
! write adjusting startpoint to be inside limits
                val=axarr(1)%axmin+0.1D0*(axarr(1)%axmax-axarr(1)%axmin)
             endif
             write(*,73)name(1:len_trim(name)),val
73           format(/'Setting start condition for ',a,f10.5)
! first argument 1 means to extract the value, 0 means to set the value
             call condition_value(0,pcond,val,ceq)
             if(gx%bmperr.ne.0) goto 500
          endif
          mode=-1
!
          if(notop.eq.0) then
             notop=-1
! create maptop and things for storing results
! map_startpoint calculates the equilibrium and generates two start points
!             write(*,*)'Creating start point',itup
             inactive=0
             call map_startpoint(maptop,noofaxis,axarr,seqxyz,inactive,ceq)
             if(gx%bmperr.ne.0) goto 500
! create array of equilibrium records for saving results
! if larger than 500 I get segmentation fault ,,,,
             saveq=maxsavedceq
             call create_saveceq(maptop%saveceq,saveq)
             if(gx%bmperr.ne.0) goto 900
! initiate line counter (redundant) ... maybe if several step separate?
!             if(seqxyz(2).ne.0) then
!                write(*,*)'step_separate seqy: ',seqxyz(2)
!             endif
          else
! we generate a second or later startpoint for another phase
! note that maptop is allocated a new map_node linked from this
!             write(*,*)'Creating next start point',itup
             inactive=0
             call map_startpoint(maptop,noofaxis,axarr,seqxyz,inactive,ceq)
             if(gx%bmperr.ne.0) then
                goto 500
             endif
          endif
          firstline=.TRUE.
! find a stored line to calculate
! in this subroutine we have only one axis variable
200       continue
!          write(*,*)'Calling findline:'
          call map_findline(maptop,axarr,mapfix,mapline)
          if(gx%bmperr.ne.0) goto 500
!          write(*,*)'Back from map_findline in STEP',associated(mapline)
          ceq=>mapline%lineceq
          meqrec=>mapline%meqrec
! this is the first equilibrium along the line, create meqrec in step_separate
305       continue
!          do jk=1,ntup
!             if(stsphcs(jk).eq.-2) write(*,*)'SS phase ',jk,' dormant B'
!             if(stsphcs(jk).ge.0) write(*,*)'SS phase ',jk,' stable B'
!          enddo
!          write(*,*)'smp2A calling calceq7 for first point'
          call calceq7(mode,meqrec,mapfix,ceq)
!          write(*,*)'smp2A back from calceq7',gx%bmperr
          if(gx%bmperr.ne.0) then
! error 4187 is to set T or P less than 0.1
             if(gx%bmperr.eq.4187) then
                write(*,*)'We jump to 333'
                goto 333
             endif
             if(mapline%number_of_equilibria.eq.0) then
! We can add/subtract a small amount of axis condition if error at first step
                write(*,*)'Error at first equilibrium: ',gx%bmperr,&
                     mapline%axandir
             endif
!             write(*,*)'SMP error: ',gx%bmperr
! if step turn on grid minimizer
             write(*,*)'Turn on grid minimizer'
             if(maptop%number_ofaxis.eq.1) then
                call calceq7(mode,meqrec,mapfix,ceq)
                if(gx%bmperr.ne.0) then
                   write(kou,*)'Failed calling grid minimizer',gx%bmperr
                   gx%bmperr=0
                endif
             endif
! reset error code and take another line
!             write(*,*)'SMP2 Generating mapline%meqrec failed 1: ',gx%bmperr
             gx%bmperr=0; goto 333
          else
! calculation OK, do it again (why?) without creating meqrec, save and
! return here after taking a step using the same meqrec
380          continue
             iadd=0
             irem=0
!             write(*,*)'SMP2A calling meq_sameset from step_separate'
             call meq_sameset(irem,iadd,mapx,mapline%meqrec,&
                  mapline%meqrec%phr,inmap,ceq)
             if(gx%bmperr.ne.0) then
!                write(*,*)'SMP2A Error calling meq_sameset',gx%bmperr
                goto 333
             elseif(iadd.ne.0 .or. irem.ne.0) then
                write(*,*)'Change of phases not allowed! ',iadd,irem
                goto 333
             endif
! store the result
             call map_store(mapline,axarr,maptop%number_ofaxis,maptop%saveceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Error storing equilibrium',gx%bmperr
                goto 900
             endif
!             do jk=1,ntup
!                if(stsphcs(jk).eq.-2) write(*,*)'SS phase ',jk,' dormant C'
!                if(stsphcs(jk).ge.0) write(*,*)'SS phase ',jk,' stable C'
!             enddo
             call map_step(maptop,mapline,mapline%meqrec,mapline%meqrec%phr,&
                  axvalok,noofaxis,axarr,ceq)
!             write(*,*)'Back from map_step 2 ',mapline%more,&
!                  mapline%number_of_equilibria
             if(gx%bmperr.ne.0) then
! if error just terminate line
                write(*,*)'Error return from map_step 2: ',gx%bmperr
                mapline%more=-1
                gx%bmperr=0; goto 333
             endif
             if(mapline%more.ge.0) goto 380
          endif
333       continue
!          write(*,*)'Calling map_linend 2'
          call map_lineend(mapline,axarr(abs(mapline%axandir))%lastaxval,ceq)
          if(firstline) then
! follow the axis in the other direction
             if(gx%bmperr.ne.0) then
                write(*,*)'Removing error code',gx%bmperr
             endif
             firstline=.FALSE.
             goto 200
          endif
! finished step in both directions
500       continue
! remove any error before calculating next phase
          if(gx%bmperr.ne.0) then
             write(*,*)'Removing error code to calculate next phase',gx%bmperr
             gx%bmperr=0
          endif
       endif
! set current phase as suspended and calculate for next phase
!       call change_phase_status(entphcs(itup)%phaseix,entphcs(itup)%compset,&
       call change_phase_status(entphcs(itup)%ixphase,entphcs(itup)%compset,&
            -3,zero,ceq)
!       write(*,*)'At end of phase loop itup=',itup
    enddo phaseloop
!============================================================
! Terminate but restore all phase status, even if error above
900 continue
    val=zero
!    write(*,*)'SMP Trying to restoring original phase status',ntup
! reset ceq to be starteq !! otherwise nothing is changed
    ceq=>starteq
    do itup=1,ntup
!       write(*,910)itup,entphcs(itup)%ixphase,entphcs(itup)%compset,&
!            stsphcs(itup)
910    format('Restoring all phase status: ',4i5)
!       call change_phtup_status(itup,stsphcs(itup),val,ceq)
       call change_phase_status(entphcs(itup)%ixphase,&
            entphcs(itup)%compset,stsphcs(itup),val,ceq)
       if(gx%bmperr.ne.0) goto 1000
! trying to set status entered ...
!       write(*,911)'step_sep: restored? ',itup,entphcs(itup)%ixphase,&
!            entphcs(itup)%compset,stsphcs(itup),val
911    format(a,3i4,i6,1pe12.4)
    enddo
1000 continue
    return
  end subroutine step_separate

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine step_special_setup
!\begin{verbatim}
  subroutine step_special_setup(maptop,seqxyz,exits,starteq)
! create mapnode and tzero line and other special step commands
! maptop map node record
! seqxyz indices for map and line records
! exits is 1 or more depending on type of step
    implicit none
    integer seqxyz(*),exits
    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq
    integer jj,jp,seqz,iadd,irem,nv,saveq,lokcs,mapx,idir,seqx,seqy,kpos
    type(map_node), pointer :: mapnode
    type(map_line), pointer :: mapline
    type(map_fixph), allocatable :: mapfix
    type(meq_setup), pointer :: meqrec
    type(gtp_state_variable), pointer :: svr
    type(meq_phase), pointer :: phr
    type(gtp_condition), pointer :: pcond
    double precision xxx,yyy,zzz,fact
!    logical firstline
    character eqname*24
    integer, parameter :: maxsavedceq=1800
!
!    write(*,*)'In step_special_setup',exits
!======================================================
! create maptop, maplines and things for storing results
! we cannot use map_startpoint as we are not calculating equilibria ...
! we must allocate a maptop and its next and previous to point at itself
    allocate(maptop)
    mapnode=>maptop
! inititate status and links, maybe some of these change for other applications
    mapnode%status=0
    mapnode%noofstph=2
    mapnode%savednodeceq=-1
    mapnode%next=>mapnode
    mapnode%previous=>mapnode
    mapnode%first=>mapnode
    mapnode%number_ofaxis=1
    mapnode%nodefix%ixphase=0
    mapnode%status=0
! mapnone%lines incremented when created ??
    mapnode%lines=0
! %artxe nonzero if node with two stoichiometric phases with same composition
    mapnode%artxe=0
    mapnode%globalcheckinterval=0
    mapnode%seqx=seqxyz(1)
    mapnode%seqy=seqxyz(2)
!
! skip saving chemical potentials?
    mapnode%tpval=starteq%tpval
    mapnode%nodeceq=>starteq
    eqname='_MAPNODE_'
    jp=10
! maptop%next is the the same mapnode !!!
    seqx=maptop%next%seqx+1
!    seqy=maptop%next%seq+1
    maptop%next%seqx=seqx
    call wriint(eqname,jp,seqx)
! make a copy of ceq in a new equilibrium record with the pointer neweq
! This copy is a record in the array "eqlista" of equilibrium record, thus
! it will be updated if new composition sets are created in other threads.
    call copy_equilibrium(neweq,eqname,starteq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'Created MAPNODE ',seqx
! set the tieline_inplane or not
! For step calculation, tieline_inplane=0
! if there are more than one condition on an extensive_variable
! that is not an axis variable then no tielines in plane, tieline_inplane=-1
! If there are tie_lines in plane then tieline_inplane=1
    mapnode%tieline_inplane=0
! forgetting to do this created a crash when plotting ...
    nullify(maptop%plotlink)
! we must store 1 or 2 (=exits) lineceq using starteq
    mapnode%lines=exits
    allocate(mapnode%linehead(mapnode%lines))
!    write(*,*)'step_special_setup',maptop%seqx,exits
!    mapnode%type_of_node=0
    idir=1
    do jp=1,exits
       mapnode%linehead(jp)%axandir=idir
       idir=-1
       mapnode%linehead(jp)%number_of_equilibria=0
       mapnode%linehead(jp)%first=0
       mapnode%linehead(jp)%last=0
       mapnode%linehead(jp)%axchange=-1
       mapnode%linehead(jp)%done=0
       mapnode%linehead(jp)%status=0
       mapnode%linehead(jp)%more=1
       mapnode%linehead(jp)%termerr=0
       mapnode%linehead(jp)%firstinc=zero
! saving equilibrium pointer in lineceq
       mapnode%linehead(jp)%lineceq=>starteq
       mapnode%linehead(jp)%start=>mapnode
       mapnode%linehead(jp)%axfact=1.0D-2
! this is set to zero indicating the stable phases are saved in lineceq record
       mapnode%linehead(jp)%nstabph=0
       mapnode%linehead(jp)%lineid=seqy
       mapnode%seqy=seqy+1
       mapnode%linehead(jp)%nodfixph=0
! %more is 1 while line is calculated, 0 means terminated at axis limit
! > 0 means error code <0 means exit removed ?? or is it %done ??
       mapnode%linehead(jp)%more=1
       mapnode%lines=exits
       nullify(mapnode%linehead(jp)%end)
    enddo
!
! create array of equilibrium records for saving results
!    write(*,*)'step_special_setup create saveceq:',maxsavedceq
    saveq=maxsavedceq
    call create_saveceq(maptop%saveceq,saveq)
    if(gx%bmperr.ne.0) goto 1000
! in this subroutine we have only one axis variable
1000 continue
    return
  end subroutine step_special_setup

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine step_tzero
!\begin{verbatim}
  subroutine step_tzero(maptop,noofaxis,axarr,seqxyz,iph1,iph2,tzcond,starteq)
! calculates t for two phases where they have same Gibbs energy
! second version using step_special_setup
! There can not be any other phases
! maptop map node record
! noofaxis must be 1
! axarr array of axis records
! seqxyz indices for map and line records
! iph1 and iph2 should be phase index (compset 1 in both)
! tzcond should be condition number for T
    implicit none
    integer noofaxis,seqxyz(*),iph1,iph2,tzcond
    type(map_axis), dimension(noofaxis) :: axarr
    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq
    integer jj,jp,seqz,iadd,irem,nv,saveq,lokcs,mapx,idir,seqx,seqy,kpos
    type(map_node), pointer :: mapnode
    type(map_line), pointer :: mapline
    type(map_fixph), allocatable :: mapfix
    type(meq_setup), pointer :: meqrec
    type(gtp_state_variable), pointer :: svr
    type(meq_phase), pointer :: phr
    type(gtp_condition), pointer :: pcond
    double precision xxx,yyy,zzz,fact
!    logical firstline
    character eqname*24
    integer, parameter :: maxsavedceq=1800
! turns off convergence control for T
    integer, parameter :: inmap=1
!
!    write(*,*)'In step_tzero',iph1,iph2
    if(noofaxis.ne.1) then
       write(kou,*)'Step tzero only with one axis variable'
       goto 1000
    endif
! check that we have a tzero point
    ceq=>starteq
!
    call tzero(iph1,iph2,tzcond,yyy,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Start point is not on a tzero line'
       gx%bmperr=4399; goto 1000
    endif
! extract axis condition value
    call locate_condition(axarr(1)%seqz,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to get the value
    call condition_value(1,pcond,xxx,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,88)xxx,yyy
88  format('At x=',F10.6,' Tzero=',F10.2,10x,1pe12.4)
!======================================================
    call step_special_setup(maptop,seqxyz,2,starteq)
    if(gx%bmperr.ne.0) goto 1000
!
    mapnode=>maptop
!    write(*,*)'step_tzero creating maplines'
    tzstep: do jp=1,2
       mapline=>mapnode%linehead(jp)
       eqname='_MAPLINE_'
       kpos=10
       seqy=maptop%seqy+1
       call wriint(eqname,kpos,seqy)
       call copy_equilibrium(mapnode%linehead(jp)%lineceq,eqname,&
            mapnode%nodeceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error creating equilibrium: ',eqname
          goto 1000
       endif
!       write(*,*)'step_tzero created mapline ',seqy
       maptop%seqy=seqy
       mapnode%linehead(jp)%lineid=seqy
       mapnode%linehead(jp)%nodfixph=0
! mapline%more is positive for line to be calculated, 0 means end at axis limit
       mapnode%linehead(jp)%more=1
       ceq=>mapline%lineceq
! A very small first axis increment, extract axis condition value
       call locate_condition(axarr(1)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to get the value
       call condition_value(1,pcond,xxx,ceq)
       if(gx%bmperr.ne.0) goto 1000
       fact=1.0D-2
       idir=mapline%axandir
!       write(*,*)'axis direction: ',idir,xxx
       tzlimits: do while(.TRUE.)
          xxx=xxx+fact*idir*axarr(1)%axinc
          if(xxx.lt.axarr(1)%axmin .or. xxx.gt.axarr(1)%axmax) exit tzlimits
          call condition_value(0,pcond,xxx,ceq)
          call tzero(iph1,iph2,tzcond,yyy,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'TZERO step ',jp,' ended with error ',gx%bmperr
             gx%bmperr=0; cycle tzstep
!          else
!             write(*,88)xxx,yyy,fact
          endif
          call map_store(mapline,axarr,maptop%number_ofaxis,maptop%saveceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error storing equilibrium',gx%bmperr
             gx%bmperr=0; cycle tzstep
          endif
! save missing .........
          fact=min(2.0d0*fact,1.0d0)
       enddo tzlimits
       if(xxx.lt.axarr(1)%axmin) then
          xxx=max(axarr(1)%axmin,1.0D-6)
          call condition_value(0,pcond,xxx,ceq)
          call tzero(iph1,iph2,tzcond,yyy,ceq)
       elseif(xxx.gt.axarr(1)%axmax) then
          xxx=min(axarr(1)%axmax,0.999999D0)
          call condition_value(0,pcond,xxx,ceq)
          call tzero(iph1,iph2,tzcond,yyy,ceq)
       endif
!       write(*,88)xxx,yyy
       call map_store(mapline,axarr,maptop%number_ofaxis,maptop%saveceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error storing equilibrium',gx%bmperr
          gx%bmperr=0; cycle tzstep
       endif
    enddo tzstep
!
1000 continue
    return
  end subroutine step_tzero

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine step_tzero1
!\begin{verbatim}
  subroutine step_tzero1(maptop,noofaxis,axarr,seqxyz,iph1,iph2,tzcond,starteq)
! ORIGINAL tzero step NO LONGER USED
! calculates t for two phases where they have same Gibbs energy
! There can not be any other phases
! maptop map node record
! noofaxis must be 1
! axarr array of axis records
! seqxyz indices for map and line records
! iph1 and iph2 should be phase index (compset 1 in both)
! tzcond should be condition number for T
    implicit none
    integer noofaxis,seqxyz(*),iph1,iph2,tzcond
    type(map_axis), dimension(noofaxis) :: axarr
    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq
    integer jj,jp,seqz,iadd,irem,nv,saveq,lokcs,mapx,idir,seqx,seqy,kpos
    type(map_node), pointer :: mapnode
    type(map_line), pointer :: mapline
    type(map_fixph), allocatable :: mapfix
    type(meq_setup), pointer :: meqrec
    type(gtp_state_variable), pointer :: svr
    type(meq_phase), pointer :: phr
    type(gtp_condition), pointer :: pcond
    double precision xxx,yyy,zzz,fact
!    logical firstline
    character eqname*24
    integer, parameter :: maxsavedceq=1800
! turns off convergence control for T
    integer, parameter :: inmap=1
!
!    write(*,*)'In step_tzero',iph1,iph2
    if(noofaxis.ne.1) then
       write(kou,*)'Step tzero only with one axis variable'
       goto 1000
    endif
! check that we have a tzero point
    ceq=>starteq
!
    call tzero(iph1,iph2,tzcond,yyy,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Start point is not on a tzero line'
       gx%bmperr=4399; goto 1000
    endif
! extract axis condition value
    call locate_condition(axarr(1)%seqz,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to get the value
    call condition_value(1,pcond,xxx,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,88)xxx,yyy
88  format('At x=',F10.6,' Tzero=',F10.2,10x,1pe12.4)
!======================================================
! create maptop, maplines and things for storing results
!    write(*,*)'Creating start point'
! we cannot use map_startpoint as we are not calculating equilibria ...
!    call map_startpoint(maptop,noofaxis,axarr,seqxyz,inactive,ceq)
!    if(gx%bmperr.ne.0) goto 500
! we must allocate a maptop and its next and previous to point at itself
    allocate(maptop)
    mapnode=>maptop
! inititate status and links
    mapnode%status=0
    mapnode%noofstph=2
    mapnode%savednodeceq=-1
    mapnode%next=>mapnode
    mapnode%previous=>mapnode
    mapnode%first=>mapnode
    mapnode%number_ofaxis=noofaxis
    mapnode%nodefix%ixphase=0
    mapnode%status=0
! mapnone%lines incremented when created ??
    mapnode%lines=0
! %artxe nonzero if node with two stoichiometric phases with same composition
    mapnode%artxe=0
    mapnode%globalcheckinterval=0
    mapnode%seqx=seqxyz(1)
    mapnode%seqy=seqxyz(2)
!
! skip saving chemical potentials?
    mapnode%tpval=ceq%tpval
    mapnode%nodeceq=>ceq
    eqname='_MAPNODE_'
    jp=10
! maptop%next is the most recent created mapnode ??
    seqx=maptop%next%seqx+1
    maptop%next%seqx=seqx
    call wriint(eqname,jp,seqx)
! make a copy of ceq in a new equilibrium record with the pointer neweq
! This copy is a record in the array "eqlista" of equilibrium record, thus
! it will be updated if new composition sets are created in other threads.
    call copy_equilibrium(neweq,eqname,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'Created MAPNODE ',seqx
! set the tieline_inplane or not
! For step calculation, tieline_inplane=0
! if there are more than one condition on an extensive_variable
! that is not an axis variable then no tielines in plane, tieline_inplane=-1
! If there are tie_lines in plane then tieline_inplane=1
    mapnode%tieline_inplane=0
! forgetting to do this created a crash when plotting ...
    nullify(maptop%plotlink)
! we must store 2 lineceq using ceq, one in each direction
    mapnode%lines=2
    allocate(mapnode%linehead(2))
!    write(*,*)'step_tzero created maptop',maptop%seqx
    idir=-1
    do jp=1,2
       mapnode%linehead(jp)%axandir=idir
       idir=1
       mapnode%linehead(jp)%number_of_equilibria=0
       mapnode%linehead(jp)%first=0
       mapnode%linehead(jp)%last=0
       mapnode%linehead(jp)%axchange=-1
       mapnode%linehead(jp)%done=0
       mapnode%linehead(jp)%status=0
       mapnode%linehead(jp)%more=1
       mapnode%linehead(jp)%termerr=0
       mapnode%linehead(jp)%firstinc=zero
! saving equilibrium pointer in lineceq
       mapnode%linehead(jp)%lineceq=>ceq
       mapnode%linehead(jp)%start=>mapnode
       mapnode%linehead(jp)%axfact=1.0D-2
! this is set to zero indicating the stable phases are saved in ceq record
       mapnode%linehead(jp)%nstabph=0
       seqy=mapnode%seqy
       mapnode%linehead(jp)%lineid=seqy
       mapnode%seqy=seqy+1
       mapnode%linehead(jp)%nodfixph=0
! %more is 1 while line is calculated, 0 means terminated at axis limit
! > 0 means error code <0 means exit removed ?? or is it %done ??
       mapnode%linehead(jp)%more=1
       nullify(mapnode%linehead(jp)%end)
    enddo
!
! suck
!
! create array of equilibrium records for saving results
    saveq=maxsavedceq
    call create_saveceq(maptop%saveceq,saveq)
    if(gx%bmperr.ne.0) goto 1000
! in this subroutine we have only one axis variable
200 continue
!
    tzstep: do jp=1,2
       mapline=>mapnode%linehead(jp)
       eqname='_MAPLINE_'
       kpos=10
       seqy=maptop%seqy+1
       call wriint(eqname,kpos,seqy)
       call copy_equilibrium(mapnode%linehead(jp)%lineceq,eqname,&
            mapnode%nodeceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error creating equilibrium: ',eqname
          goto 1000
       endif
!       write(*,*)'step_tzero created mapline ',seqy
       maptop%seqy=seqy
       mapnode%linehead(jp)%lineid=seqy
       mapnode%linehead(jp)%nodfixph=0
! mapline%more is positive for line to be calculated, 0 means end at axis limit
       mapnode%linehead(jp)%more=1
       ceq=>mapline%lineceq
! A very small first axis increment, extract axis condition value
       call locate_condition(axarr(1)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to get the value
       call condition_value(1,pcond,xxx,ceq)
       if(gx%bmperr.ne.0) goto 1000
       fact=1.0D-2
       idir=mapline%axandir
!       write(*,*)'axis direction: ',idir,xxx
       tzlimits: do while(.TRUE.)
          xxx=xxx+fact*idir*axarr(1)%axinc
          if(xxx.lt.axarr(1)%axmin .or. xxx.gt.axarr(1)%axmax) exit tzlimits
          call condition_value(0,pcond,xxx,ceq)
          call tzero(iph1,iph2,tzcond,yyy,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'TZERO step ',jp,' ended with error ',gx%bmperr
             gx%bmperr=0; cycle tzstep
!          else
!             write(*,88)xxx,yyy,fact
          endif
          call map_store(mapline,axarr,maptop%number_ofaxis,maptop%saveceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error storing equilibrium',gx%bmperr
             gx%bmperr=0; cycle tzstep
          endif
! save missing .........
          fact=min(2.0d0*fact,1.0d0)
       enddo tzlimits
       if(xxx.lt.axarr(1)%axmin) then
          xxx=max(axarr(1)%axmin,1.0D-6)
          call condition_value(0,pcond,xxx,ceq)
          call tzero(iph1,iph2,tzcond,yyy,ceq)
       elseif(xxx.gt.axarr(1)%axmax) then
          xxx=min(axarr(1)%axmax,0.999999D0)
          call condition_value(0,pcond,xxx,ceq)
          call tzero(iph1,iph2,tzcond,yyy,ceq)
       endif
!       write(*,88)xxx,yyy
       call map_store(mapline,axarr,maptop%number_ofaxis,maptop%saveceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error storing equilibrium',gx%bmperr
          gx%bmperr=0; cycle tzstep
       endif
    enddo tzstep
!
1000 continue
! ORIGINAL tzero step NO LONGER USED
    return
  end subroutine step_tzero1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine step_scheil
!\begin{verbatim}
  subroutine step_scheil(maptop,noofaxis,axarr,seqxyz,starteq)
! calculates a Scheil-Gulliver solidification simulation
! maptop map node record
! noofaxis must be 1
! axarr array of axis records
! seqxyz indices for map and line records
! starteq is an equilibrium with just liquid stable
    implicit none
    integer noofaxis,seqxyz(*)
    type(map_axis), dimension(noofaxis) :: axarr
    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq
    integer jj,jp,seqz,iadd,irem,nv,saveq,lokcs,mapx,idir,seqx,seqy,kpos
    integer inactive(4),mode,nc,nsch,liquid
    type(map_node), pointer :: mapnode
    type(map_line), pointer :: mapline
    type(map_fixph), allocatable :: mapfix
    type(meq_setup), pointer :: meqrec
    type(gtp_state_variable), pointer :: svr
    type(meq_phase), pointer :: phr
    type(gtp_condition), pointer :: pcond,firstcond,axcond
    double precision xxx,yyy,zzz,fact,fact1,axvalok,npliqval,liqfrac(20)
    character eqname*24,phname*24,npliq*24,encoded*72
    integer, parameter :: maxsavedceq=1800
! turns off convergence control for T
    integer, parameter :: inmap=1
    logical solids
! needed to store links to condition values
    TYPE smp_scheil_condval
! these pointers must be updated for each new line (equilibrium)
       type(gtp_condition), pointer :: p1
    end type smp_scheil_condval
! These two arrays keep track of conditions and liquid compositis
! the first is pointers to the condition record, the second is statevariable id
    type(smp_scheil_condval), dimension(20) :: scheilval
    TYPE(gtp_state_variable), target, dimension(20) :: scheilsvr
!
    write(*,*)'In step_scheil'
    if(noofaxis.ne.1) then
       write(kou,*)'Scheil simulations use one axis variable'
       goto 1000
    endif
! axis condition must be T, extract its value
    call locate_condition(axarr(1)%seqz,pcond,starteq)
    if(gx%bmperr.ne.0) goto 1000
    if(pcond%statev.ne.1) then
! pcond%statev=1 means T
       write(*,*)'Axis condition must be T'
       gx%bmperr=4399; goto 1000
    endif
! first argument 1 means to get the value
    axcond=>pcond
    call condition_value(1,pcond,xxx,starteq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,'(a,F10.2)')'Scheil start T=',xxx
!
    inactive=0
! inactive(1)=-1 means only one exit point with direcition -1
    inactive(1)=-1
! generate step/map datastructure needed for plotting and phase set changes.
    call map_startpoint(maptop,noofaxis,axarr,seqxyz,inactive,starteq)
    if(gx%bmperr.ne.0) goto 1000
! There should be two maplines generated, the stable phase should be the liquid
! but do not be fuzzy, one may quech a two-phase mixture
!    write(*,*)'Scheil step 1',allocated(maptop%linehead)
!    write(*,*)'Scheil lineheads: ',size(maptop%linehead),&
!         maptop%linehead(1)%axandir
! create array of equilibrium records for saving results
    seqy=maxsavedceq
    call create_saveceq(maptop%saveceq,seqy)
    if(gx%bmperr.ne.0) goto 1000
! Mark this as a Scheil step
    maptop%type_of_node=3
! ensure plotlink is nullified!!
    nullify(maptop%plotlink)
! initiate node counter done, line counter will be incremented
    if(maptop%seqx.gt.1) write(*,85)maptop%seqx,maptop%seqy+1
85  format('Previous step/map results saved'/&
         'New mapnode/line equilibria indices will start from: ',i3,i5)
! take the first (only) line created by map_startpoint
!    write(*,*)'Calling map_findline'
    call map_findline(maptop,axarr,mapfix,mapline)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'Back from map_findline in Scheil'
    ceq=>mapline%lineceq
    meqrec=>mapline%meqrec
    mode=-1
    call calceq7(mode,meqrec,mapfix,ceq)
    if(gx%bmperr.ne.0) goto 1000
    xxx=ceq%tpval(1)
    if(meqrec%nstph.gt.1) then
       write(*,*)'More than one phase stable at startpoint'
       gx%bmperr=4399; goto 1000
    endif
! check stable phase is liquid
    call get_phasetup_name(meqrec%phr(meqrec%stphl(1))%phtupix,phname)
    if(gx%bmperr.ne.0) goto 1000
    liquid=meqrec%phr(meqrec%stphl(1))%iph
    write(*,*)'Stable phase at start: ',trim(phname),liquid
    npliq='NP('//trim(phname)//') '
!=======================================================
! create special result array to save current fraction of liquid
!    allocate(mapline%stepresultid(1))
!    mapline%stepresultid(1)=npliq
! extract relevant conditions and store in scheilval and scheilsvr
    firstcond=>ceq%lastcondition%next
    pcond=>firstcond
    nc=0
    nsch=0
    ploop: do while(.TRUE.)
! if %active nonzero the condition is not active
       if(pcond%active.ne.0) cycle ploop
       nc=nc+1
! to prevent eternal loop
       if(nc.gt.20) exit ploop
!       write(*,'(a,i3,a,i5)')'Condition ',nc,' type ',pcond%statev
       if(pcond%statev.lt.0) then
          write(*,*)'Fix phases not allowed as conditions'
          gx%bmperr=4399; goto 1000
       endif
       svr=>pcond%statvar(1)
!       write(*,*)'State variable id: ',svr%statevarid,svr%argtyp
! statvarid<10 means potential, allow and ignore
       if(svr%statevarid.le.10) goto 100
! 11 <= statvarid <=15 are G, H etc, not allowed.  Neither is Y
       if(svr%statevarid.le.15 .or. svr%statevarid.ge.20) then
          write(*,*)'Illegal condition for Scheil simulation',svr%statevarid
          gx%bmperr=4399; goto 1000
       endif
! Allowed state variables are N, X, B and W without phase specification
! argtyp=0 means total such as N=1
       if(svr%argtyp.eq.0) goto 100
! argtyp=1 means component, >1 other means phase or compset specification
       if(svr%argtyp.gt.1) then
          write(*,*)'Condition has wrong type of arguments: ',svr%argtyp
          gx%bmperr=4399; goto 1000
       endif
       if(pcond%symlink1.gt.0) then
! value must not be a symbol
          write(*,*)'Condition value must not be a symbol'
          gx%bmperr=4399; goto 1000
       endif
       nsch=nsch+1
       scheilval(nsch)%p1=>pcond
! save state variable but change it to include liquid phase index
       scheilsvr(nsch)=svr
! replace argtyp and add phase and compset
       scheilsvr(nsch)%argtyp=3
       scheilsvr(nsch)%phase=liquid
       scheilsvr(nsch)%compset=1
!       write(*,'(a,i3,F10.6)')'Condition value: ',nsch,pcond%prescribed
! Puuuuuh, condition allowed, link to its current value
100    continue
       pcond=>pcond%next
! current value
       if(associated(pcond,firstcond)) exit ploop
    enddo ploop
!    write(*,'(a,i3,a,i3)')'Found ',nc,' active conditions and saved ',nsch
! test that we can extract (and set) liquid conditions and state variable
    do nc=1,nsch
       svr=>scheilsvr(nc)
       call state_variable_val(svr,xxx,ceq)
!       write(*,'(a,i3,2F10.6)')'Liquid initial conditions: ',&
!            nc,scheilval(nc)%p1%prescribed,xxx
    enddo
! initial
    npliqval=one
    solids=.FALSE.
! Now find T when first solid phase will appear
! mapx does not seem to be used, inmap=1 turn off T convergence test(?)
! all data in meqrec was set calling calceq7 above
! axis conditio
! first argument 1 means to get the value
    call condition_value(1,axcond,xxx,ceq)
    if(gx%bmperr.ne.0) goto 1000
    irem=0; iadd=0; nc=0
! iadd=-1 turn on verbose in meq_sameset
!    iadd=-1
! large step before first solid appears
    fact1=1.0D1
    axarr(1)%axinc=fact1*axarr(1)%axinc
    axvalok=xxx
!==================================================== big loop
    node: do while(.TRUE.)
!   follow axis including nodepoints with phase changes
!   start with small steps
!       fact=1.0D-2
       line: do while(iadd.le.0 .and. irem.eq.0)
!         follow line until a nodepoint
!          axarr(1)%axval=axarr(1)%axval-axarr(1)%axinc
          if(solids) then
! update the liquid composition
! We have located the pcond records for each new line below
!             write(*,*)'Update liquid composition at T=',ceq%tpval(1)
             do nc=1,nsch
! this call extract the liquid composition
                svr=>scheilsvr(nc)
                call state_variable_val(svr,liqfrac(nc),ceq)
                if(gx%bmperr.ne.0) then
                   write(*,*)'Error extracting liquid composition'
                   goto 1000
                endif
! and this sets it as the overall composition
                call condition_value(0,scheilval(nc)%p1,liqfrac(nc),ceq)
                if(gx%bmperr.ne.0) then
                   write(*,*)'Error setting new liquid composition'
                   goto 1000
                endif
             enddo
             call get_state_var_value(npliq,yyy,encoded,ceq)
             if(gx%bmperr.ne.0) gx%bmperr=0
             npliqval=npliqval*yyy
             write(*,'(a,F7.2,"% ",F7.2,": ",10(1x,F8.4))')'Liquid:',&
                  1.0D2*npliqval,ceq%tpval(1),(liqfrac(nc),nc=1,nsch)
! turn on debug info in meq_sameset
!             iadd=-1
          endif
! take a step in the axis variable T
          call map_step2(maptop,mapline,meqrec,meqrec%phr,axvalok,1,axarr,ceq)
          if(gx%bmperr.ne.0) goto 1000
          if(ceq%tpval(1).lt.axarr(1)%axmin) then
             write(*,*)'At low T limit ',axarr(1)%axmin
             goto 900
          endif
! calculate until a phase change
!          write(*,*)'Calling meq_sameset',ceq%tpval(1),npliqval
          call meq_sameset(irem,iadd,mapx,mapline%meqrec,mapline%meqrec%phr,&
               inmap,ceq)
!          write(*,*)'Back from meq_sameset',ceq%tpval(1),gx%bmperr
          if(iadd.eq.0 .and. irem.eq.0) then
! Store the equilibrium along the line
             call map_store(mapline,axarr,1,maptop%saveceq)
!             write(*,*)'Stored calculated equilibrium'
             if(gx%bmperr.ne.0) then
                write(*,*)'Error storing equilibria',gx%bmperr
                goto 1000
             endif
          endif
       enddo line
! exit line loop when iadd or irem nonzero, i.e. new set of phases
       if(.not.solids) then
! if solids FALSE set it TRUE
          solids=.TRUE.
          axarr(1)%axinc=axarr(1)%axinc/fact1
          fact1=1.0D0
       endif
! Maybe not store here because the T is not correct
!       call map_store(mapline,axarr,1,maptop%saveceq)
!       write(*,*)'Stored calculated equilibrium'
       if(gx%bmperr.ne.0) then
          write(*,*)'Error storing equilibria',gx%bmperr
          goto 1000
       endif
! use map_calcnode to create new mapnode and mapline.  
! We should not set any fix phases, just continue along the axis
! as with a step command with different sets of stable phases
       call map_calcnode(irem,iadd,maptop,mapline,meqrec,axarr,ceq)
! in map_calcnode a new _MAPNODE and _MAPLINE is created with the new set
! of phases.  Store the end point of the line
       nullify(maptop%plotlink)
! Terminate the current line, must be after calcnode ...
       call map_lineend(mapline,axarr(1)%lastaxval,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Rest error ',gx%bmperr
          gx%bmperr=0
       endif
       write(*,*)'Per cent liquid and T',1.0D2*npliqval,ceq%tpval(1)
       if(.not.(npliqval.gt.0.01)) then
! terminate if npliqval<0.01 BUT IT DOES NOT WORK ???
          write(*,*)'Terminating as liquid fraction less than 1%'
          goto 900
!       else
!          if(npliqval.gt.0.01) then
!             write(*,*)'Terminating as liquid fraction less than 1%'
!             goto 900
!          endif
       endif
! The Scheil simulation continue along the same axis with new set of phases.
       call map_findline(maptop,axarr,mapfix,mapline)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error return from map_findline, terminating'
          goto 1000
       endif
       ceq=>mapline%lineceq
!       write(*,*)'SMP2A calling calceq7 after findline,',allocated(mapfix),mode
! Evidently we have to call calceq7 to initiate meqrec ??
       meqrec=>mapline%meqrec
       call calceq7(mode,meqrec,mapfix,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Failed calling calceq7',gx%bmperr
          goto 1000
       endif
! check if zero fraction of liquid here
       call get_state_var_value(npliq,yyy,encoded,ceq)
       write(*,*)'SMP2A Scheil liquid fraction: ',yyy
       if(yyy.lt.0.03) then
! Terminate the current line
          call map_lineend(mapline,axarr(1)%lastaxval,ceq)
          goto 900
       endif
! we have to locate the condition records for the liquid comp in the new ceq
       firstcond=>ceq%lastcondition%next
       pcond=>firstcond
       ploop2: do while(.TRUE.)
          if(pcond%active.ne.0) cycle ploop2
          svr=>pcond%statvar(1)
          do nc=1,nsch
             if(svr%statevarid.eq.scheilsvr(nc)%statevarid .and. &
                  svr%argtyp.eq.1 .and.&
                  svr%component.eq.scheilsvr(nc)%component) then
                scheilval(nc)%p1=>pcond
!                write(*,*)'Found scheil condition in new ceq: ',nc 
             endif
          enddo
          pcond=>pcond%next
          if(associated(pcond,firstcond)) exit ploop2
!          write(*,*)'Looping conditions in new ceq'
       enddo ploop2
!       write(*,*)'Node T=',ceq%tpval(1)
    enddo node
    write(*,*)'Never here!'
!
!===========================================
! exit here if no liquid left of at low T limit
900 continue
! maybe clean up?
1000 continue
    return
  end subroutine step_scheil

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine step_paraequil
!\begin{verbatim}
  subroutine step_paraequil(maptop,noofaxis,axarr,seqxyz,tupix,fastelno,starteq)
! calculates a paraequilibrium diagram
! maptop map node record
! noofaxis must be 1
! axarr array of axis records
! seqxyz indices for map and line records
! starteq is an equilibrium with just two phases stable
! tupix are phasetuple indices of two phases
! fastelno fast diffusing component index
!
! TO BE MODIFIED
!
! we will use the same overall conditions except for the carbon
    implicit none
    integer noofaxis,seqxyz(*),tupix(*),fastelno
    type(map_axis), dimension(noofaxis) :: axarr
    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq
    integer jj,jp,seqz,iadd,irem,nv,saveq,lokcs,mapx,idir,seqx,seqy,kpos
    integer inactive(4),mode,nc,npara,liquid
    type(map_node), pointer :: mapnode
    type(map_line), pointer :: mapline
    type(map_fixph), allocatable :: mapfix
    type(meq_setup), pointer :: meqrec
    type(gtp_state_variable), target :: fastxsvr,fastmusvr
    type(gtp_state_variable), target :: matrixsvr,growxsvr
    type(gtp_state_variable), pointer :: svr
    type(meq_phase), pointer :: phr
    type(gtp_condition), pointer :: pcond,axcond
    double precision xxx,yyy,zzz,fact,fact1,axvalok
    character eqname*24,phname*24,npliq*24,encoded*72,setmucond*64
    integer, parameter :: maxsavedceq=1800
! temporary storage of results
    double precision xpara(2)
! turns off convergence control for T
    integer, parameter :: inmap=1
! needed to store links to condition values
    TYPE smp_paraequil_condval
! these pointers must be updated for each new line (equilibrium)
       type(gtp_condition), pointer :: p1
    end type smp_paraequil_condval
! These two arrays keep track of conditions and liquid compositis
! the first is pointers to the condition record, the second is statevariable id
!    type(smp_paraequil_condval), dimension(20) :: paraval
!    TYPE(gtp_state_variable), target, dimension(20) :: parasvr
!
!    write(*,*)'In step_paraequil',tupix(1),tupix(2),fastelno
    if(noofaxis.ne.1) then
       write(kou,*)'Paraequilibrium simulations one axis variable'
       goto 1000
    endif
!    ceq=>starteq
    jp=1
    findxcond: do while(.true.)
! find the condition on the amount of the fast diffusing element
! ?? does this loop through all conditions number 1..n? YES
       call locate_condition(jp,pcond,starteq)
       if(gx%bmperr.eq.4295) then
! this error code means no more conditions
          gx%bmperr=0; exit findxcond
       endif
       if(gx%bmperr.ne.0) goto 1000
! skip conditions not active
       if(pcond%active.ne.0) cycle findxcond
       svr=>pcond%statvar(1)
!       write(*,*)'findcond: ',jp,svr%statevarid,svr%argtyp,svr%component
       if(svr%component.eq.fastelno) then
          if(svr%argtyp.ne.1) then
             write(*,*)'Problem, condition not on overall fraction'
             stop 'paraeq 3'
          endif
!          fastxcondno=jp
! here it should be assigned, not a pointer
          fastxsvr=svr
       endif
! avoid eternal loop?
       jp=jp+1
       if(jp.gt.20) stop 'eternal loop in step_paraeq'
    enddo findxcond
!    write(*,*)'Calling calc_paraeq first time',tupix(1),tupix(2),fastelno
! check we can calculate a paraequilibrium
    call calc_paraeq(tupix,fastelno,xpara,starteq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Sorry, cannot calculate an initial paraequilibrium',gx%bmperr
       goto 1000
    endif
!    write(*,'(a,2F10.6)')'first paraeq:',xpara(1),xpara(2)
!
!    gx%bmperr=4399; goto 1000
! =================================================================    
    inactive=0
! inactive(1)=-1 is used when only one exit point with direcition -1
! generate step/map datastructure needed for plotting and phase set changes.
! in map_startpoint an equilibrium will be calculated and maplines created
    call map_startpoint(maptop,noofaxis,axarr,seqxyz,inactive,starteq)
    if(gx%bmperr.ne.0) goto 1000
! create storage area for results
!    write(*,*)'Back from map_startpoint'
    call create_saveceq(maptop%saveceq,maxsavedceq)
    if(gx%bmperr.ne.0) goto 1000
! Mark this as a paraequil step
    maptop%type_of_node=4
! ensure plotlink is nullified!!
    nullify(maptop%plotlink)
!    write(*,*)'Taking the first line'
! take the first line created by map_startpoint
    call map_findline(maptop,axarr,mapfix,mapline)
    if(gx%bmperr.ne.0) goto 1000
    ceq=>mapline%lineceq
! meqrec contain information from the calculated equilibrium
    meqrec=>mapline%meqrec
    mode=-1
    call locate_condition(axarr(1)%seqz,axcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
!----------------------------------------------- line loop
    jp=0
    lineloop: do while(.TRUE.)
! there will be no phase changes during the STEP command, no new nodes
       jp=jp+1
!       write(*,*)'Calculating paraequilibrium',jp
       call calc_paraeq(tupix,fastelno,xpara,ceq)
       if(gx%bmperr.ne.0) then
! terminate the line and check if more lines
          goto 500
       endif
! first argument 1 means to get the value
       call condition_value(1,axcond,xxx,ceq)
       if(gx%bmperr.ne.0) goto 1000
       write(*,'(a,F12.6,": ",2F10.6)')'paraeq:',xxx,xpara(1),xpara(2)
! calculation OK, save it
       call map_store(mapline,axarr,1,maptop%saveceq)
!       write(*,*)'Stored calculated equilibrium'
       if(gx%bmperr.ne.0) then
          write(*,*)'Error storing equilibria',gx%bmperr
          goto 1000
       endif
! take a step, at second line the step is zero ... why??
       call map_step2(maptop,mapline,meqrec,meqrec%phr,axvalok,1,axarr,ceq)
       if(gx%bmperr.ne.0) goto 500
! when outside limits aapline%more is negative
       if(mapline%more.lt.0) then
! this indicate outside axis limits, call map_findline or finish
          call map_lineend(mapline,axarr(abs(mapline%axandir))%lastaxval,ceq)
          goto 510
       endif
       cycle lineloop
! treating problems 
500    continue
       if(gx%bmperr.ne.0) then
          write(*,*)'SMP2A error in step_paraequil',gx%bmperr
! terminate the line, error code cleared
          call map_lineend(mapline,axarr(mapline%axandir)%lastaxval,ceq)
! some errors maybe fatal 
       endif
510    continue
! take another line created by map_startpoint
       
       call map_findline(maptop,axarr,mapfix,mapline)
       if(gx%bmperr.ne.0) goto 1000
       if(.not.associated(mapline)) then
!          write(*,*)'SMP2A no more lines'
!          call list_conditions(kou,ceq)
          exit lineloop
       endif
       ceq=>mapline%lineceq
! axcond changed because ceq changed!!
!       write(*,*)'New line, change axis condition record'
!       call list_conditions(kou,ceq)
       call locate_condition(axarr(1)%seqz,axcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
! Wow, forgot > 
       svr=>axcond%statvar(1)
       call state_variable_val(svr,xxx,ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'Next line start at: ',xxx
!       call list_conditions(kou,ceq)
! first argument 0 means to set the value NOT ALWAYS T BEWARE
!       call condition_value(0,axcond,xxx,ceq)
!       if(gx%bmperr.ne.0) goto 1000
!       call list_conditions(kou,ceq)
! meqrec contain information from the calculated equilibrium
       meqrec=>mapline%meqrec
    enddo lineloop
!===========================================
! exit here when followed the line in both directions  remove all axcond
900 continue
! maybe clean up?
! Allow plotting tie-lines
    maptop%tieline_inplane=1
1000 continue
!    write(*,*)'Finished step_paraequil, list condition?'
!    call list_conditions(kou,ceq)
    return
  end subroutine step_paraequil

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine auto_startpoints
!\begin{verbatim}
  subroutine auto_startpoints(maptop,noofaxis,axarr,seqxyz,starteqs)
! Calculates 5 equilibria and store them as start points for mapping
! maptop map node record
! noofaxis must be 2
! axarr array of axis records
! seqxyz indices for map and line records
! starteq equilibrium record for starting
    implicit none
    integer noofaxis,seqxyz(*)
    type(map_axis), dimension(noofaxis) :: axarr
!    TYPE(gtp_equilibrium_data), pointer :: starteq
    TYPE(starteqlista), dimension(*) :: starteqs
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
! genrate one startpoint in each corner and one in the center
! For the corner add one directions along each axis
! For the center add 4 directions, totally 12 lines
! For isothermal sections one corner startpoint will be lost     
! startpoint 0.02x, 0.02y; direction +x and +y
! startpoint 0.94x, 0.02y; direction -x and +y
! startpoint 0.94x, 0.94y; direction -x and -y (lost in isothermal section)
! startpoint 0.02x, 0.94y; direction +x and -y
! startpoint 0.3x, 0.3y; all 4 directions (should work also in isothermal)
    integer seqz1,seqz2,j1,j2,mode,nss
    double precision xx1,xx2
    TYPE(gtp_equilibrium_data), pointer :: ceq,neweq,starteq
    type(gtp_condition), pointer :: pcond1,pcond2
    double precision, dimension(2), parameter :: x1=[0.02,0.92]
    double precision, dimension(2), parameter :: x2=[0.02,0.92]
    character*24 eqname
!
    if(noofaxis.ne.2 .or. &
         btest(globaldata%status,GSNOAUTOSP)) goto 1000
!    goto 1000
! the rest here works but not converting the startpoint to lines.
    write(*,*)'SMP *** in auto_startpoints'
    ceq=>starteqs(1)%p1
    mode=1
    eqname='_STARTEQ_00'
    nss=0
! loop for corners
100 continue
    cycle1: do j1=1,2
!-----------
       xx2=axarr(2)%axmin+x2(j2)*(axarr(2)%axmax-axarr(2)%axmin)
       seqz2=axarr(2)%seqz
       call locate_condition(seqz2,pcond2,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'SMP failed find 2nd condition ',j1,j2
          gx%bmperr=0
          cycle cycle1
       endif
! first argument 1 means get value, 0 means set value
       call condition_value(0,pcond2,xx2,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error setting start point condition',gx%bmperr
          gx%bmperr=0; cycle cycle1
       endif
       cycle2: do j2=1,2
          write(*,*)'SMP auto ',j1,j2
          nss=nss+1
          xx1=axarr(1)%axmin+x1(j1)*(axarr(1)%axmax-axarr(1)%axmin)
          seqz1=axarr(1)%seqz
          call locate_condition(seqz1,pcond1,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'SMP failed find first condition ',j1,j2
             gx%bmperr=0
             cycle cycle2
          endif
! first argument 1 means get value, 0 means set value
          call condition_value(0,pcond1,xx1,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error setting start point condition',gx%bmperr
             gx%bmperr=0; cycle cycle2
          endif
! calculate equilibrium
!          write(*,130)'SMP startpoint: ',nss,xx1,xx2
130       format(a,i3,2(1pe14.4))
!          call list_conditions(kou,ceq)
          call calceq2(mode,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'SMP failed calculate startpoint'
             gx%bmperr=0
          else
! enter a start equilibrum with two directions
             write(*,*)'SMP eqname: ',eqname
             call incunique(eqname(10:11))
             write(*,*)'SMP eqname: ',eqname
             call copy_equilibrium(neweq,eqname,ceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Failed to store starteq: ',trim(eqname),gx%bmperr
                gx%bmperr=0; cycle cycle2
             endif
             write(*,*)'SMP Created equilibrium: ',trim(eqname),neweq%eqno
             neweq%multiuse=20+nss
! create the list, ceq is always same equilibrium as stareq
             neweq%nexteq=ceq%nexteq
             starteq%nexteq=neweq%eqno
          endif
       enddo cycle2
    enddo cycle1
! a start point in the middle
500 continue
    xx1=0.7*axarr(1)%axmin+0.3*axarr(1)%axmax
    seqz1=axarr(1)%seqz
    call locate_condition(seqz1,pcond1,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'SMP failed find first condition',3,3
       gx%bmperr=0
       goto 1000
    endif
! first argument 1 means get value, 0 means set value
    call condition_value(0,pcond1,xx1,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error setting first central point condition',gx%bmperr
       gx%bmperr=0; goto 1000
    endif
    xx2=0.6*axarr(2)%axmin+0.4*axarr(2)%axmax
    seqz2=axarr(2)%seqz
    call locate_condition(seqz2,pcond2,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'SMP failed find 2nd condition ',j1,j2
       gx%bmperr=0
       goto 1000
    endif
! first argument 1 means get value, 0 means set value
    call condition_value(0,pcond2,xx2,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error setting second central point condition',gx%bmperr
       gx%bmperr=0; goto 1000
    endif
! calculate equilibrium
!    write(*,130)'SMP startpoint: ',5,xx1,xx2
!    call list_conditions(kou,ceq)
    call calceq2(mode,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'SMP failed calculate startpoint'
       gx%bmperr=0
    else
! enter a start equilibrum with two directions
       call incunique(eqname(10:11))
       call copy_equilibrium(neweq,eqname,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Failed to store start equilibrium',gx%bmperr
          gx%bmperr=0; goto 1000
       endif
       neweq%multiuse=30 
       neweq%nexteq=starteq%nexteq
       starteq%nexteq=neweq%eqno
    endif
1000 continue
    write(*,*)'SMP *** leaving auto_startpoint'
    return
  end subroutine auto_startpoints

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine reset_plotoptions
!\begin{verbatim}
  subroutine reset_plotoptions(graphopt,plotfile,textlabel)
! if new axis then reset default plot options
! plot ranges and their defaults
    character plotfile*(*)
    type(graphics_options) :: graphopt
    type(graphics_textlabel), pointer :: textlabel
!\end{verbatim}
    integer savebit
    graphopt%gibbstriangle=.FALSE.
    graphopt%rangedefaults=0
! axistype 0 is linear, 1 is logarithmic
    graphopt%axistype=0
! labeldefaults(1) is the title!!!
    graphopt%labeldefaults=0
    graphopt%tielines=0
    graphopt%plotmin=zero
    graphopt%dfltmin=zero
    graphopt%plotmax=one
    graphopt%scalefact=one
    graphopt%dfltmax=one
    graphopt%appendfile=' '
! do not reset font!
!    graphopt%font='Arial'
! This is confused ... GRWIN=0 if WIndows, GRWIN=1 if not windows ... SUCK
!    if(btest(graphopt%status,GRWIN)) savebit=1
! if the bit GRKEEP is set it should remain set
!    savebit=0
!    if(btest(graphopt%status,GRKEEP)) savebit=1
!    if(savebit.ne.0) graphopt%status=ibset(graphopt%status,GRKEEP)
! remove all texts ... loosing some memory ...
    nullify(graphopt%firsttextlabel)
    graphopt%labelkey='top right font "'//trim(graphopt%font)//',12" '
    nullify(graphopt%firsttextlabel)
    nullify(textlabel)
    plotfile='ocgnu'
! reset status but by default spawn plots
    graphopt%status=0
    graphopt%status=ibset(graphopt%status,GRKEEP)
! lowerleftcorner
    graphopt%lowerleftcorner=' '
! default plot terminal
    graphopt%gnutermsel=1
! plot linetype default 1
    graphopt%linetype=1
! no plot symbols
    graphopt%linewp=0
! axis tics size etc
    graphopt%textonaxis=0
! setgrid
    graphopt%setgrid=0
! do not reset plotend if set
!    plotend=plotenddefault
!    write(*,*)'Plot options reset'
    return
  end subroutine reset_plotoptions

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

