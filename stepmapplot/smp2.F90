! Data structures and routines for step/map/plot (using gnuplot)
!
! Modifications:
! For step there are 3 initial short steps after phase change
! ?? Conditions saved in node points (not fix phase but axis values)
! ?? Conditions saved in line equilibria
!
MODULE ocsmp
!
! Copyright 2012-2016, Bo Sundman, France
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
!------------------------------
!
  use liboceqplus
!
  implicit none
  character*8, parameter :: smpversion='SMP-2.10'
!
! note the type map_fixph declared in matsmin.F90 (in liboceq)
!
! MAP_NODE records are created whenever the set of stable phases changes.
! It can have links to two or more MAP_LINE records with calculated equilibra.
! Initially at least one of these map_line records are empty with just 
! information of the axis to vary and direction.
! There is a FIRST map_node record and subsequent are linked by a double linked
! list, next and previous.  All map_node records have a pointer to the first.
! There are a special use of map_node records when the node has as many stable
! phases as the line and it has only two lines connected:
! - when starting following a line from a start point
! - when the array of calculated equilibria must be saved on a file a map_node
! record is created for each unifinished line.  The map_node and map_line
! records must be saved on the random file also and initiated for the next.
! All map_node records must be kept for the next block.  It is possible that
! a line is being calculated leading to a node with no exits.  If that node 
! is removed it will be created again and all already calculated lines exiting
! will be calculated again ...
!
! MAP_LINE records are created for each line followed during the step/map.  It 
! contains links to stored gtp_equilibrium_data records and some
! additional info.
! The stored gtp_equilibrium_data records belong to an array that can be saved
! on a random file and the links are the index to reords in this array.  The
! gtp_equilibrium_data records are linked internallw with these indices also.
! The map_line record has links to two map_node records representing the
! two ends of the line.  A map_line that terminates at the end of an axis will
! have a zero link for that end of the line.
!
! These are bits in the map_line status word
! if set the whole line is inactive
  integer, parameter :: EXCLUDEDLINE=0
!
!\begin{verbatim}
  TYPE map_line
! This is record contains a list of calculated equilibria along a line
! These are pointers to map_node records at start and end of line.
     type(map_node), pointer :: start,end
! For threading this record must be separate from other threads
! This record is created by calling calceq7 the first equilibrium of line
     type(meq_setup) :: meqrec
! the active ceq record must be stored together with all others in order to be
! copied from the ceq record saved with the node record when line starts
     type(gtp_equilibrium_data), pointer :: lineceq
! this is the number of calculated equilibria for this line and the index
! of the first and last one stored.  
! The stored equilibria has an internal next link.
! lineid is a sequential index of the lines. done is negative if done
! nfixphases are the number of fixed phases replacing axis conditions
! status kan be used to delete a line
     integer number_of_equilibria,first,last,lineid,done,nfixphases,status
! This is used during mapping to identify lines that have the same fixed phases
! if we have 3 or more axis there can be 2 or more fix phases along the line??
     type(gtp_phasetuple), dimension(:), allocatable :: linefixph
! Save the phase tuplet representing the phase fix at start node here
! If it wants to be stable at first step along a line change axis direction
!     type(gtp_phasetuple) :: nodfixph
! This is the phase index in the phr array
     integer nodfixph
! We must also save the number and set of stable phases and theit amounts
! as we will have different stable phases for different lines
     integer nstabph
     type(gtp_phasetuple), dimension(:), allocatable :: stableph
     double precision,     dimension(:), allocatable :: stablepham
! axandir is set when linenode is created to the axis and direction for first
! step from the node.  It can be changed to another axis and direction
! during map and indicate the current axis with active condition
! axchange remember the equilibrum number for an axis change
     integer axandir,axchange
! more is 1 while following the line, 0 for last equilbrium, -1 when finished
! termerr is zero unless line terminated with error, -1 means exit not used
! problem is nonzero if map_problems has been called
! lasterr is the last error occured calculating this line
     integer more,termerr,problems,lasterr
! firstinc is a value to add the the axis variable for the first equilibrium
! to avoid finding the node point again.  Evenvalue is the next value
! to calculate during a step calculation.  Both set when creating the node.
! At start nodes they are zero
     double precision firstinc,evenvalue
! During map the last axis values for ALL axis are stored here
     double precision, dimension(:), allocatable :: axvals
! If tie-lines in the plane we must also check the axis values for
! the other line as we may have to change the fix phase
     double precision, dimension(:), allocatable :: axvals2
! save previous values of axvals to handle axis changes ...
     double precision, dimension(:), allocatable :: axvalx
! factor to control length of step in axis with axtive condition
     double precision :: axfact
  end TYPE map_line
!\end{verbatim}
!
!-------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE map_node
! this record organizes the step/map results.  Whenever there is a 
! change of the set of stable phases a node record is created and it
! can have links to several map_line records.  The map node record has a
! link to a gtp_equilibrium_data record (ceq) for the equilibrium at the node.
! This is copied to the map_line record when this is activated.  
! In the map_line record an axis and direction to start is stored.
! NOTE all gtp_equilibrium_data (ceq) records are pointers to the global
! array as new composition sets may be created along any line.
! The node record is identified by the set of stable phases and the
! chemical potentials of the components.  One must be able to identify the
! node as one may fins the same node following different lines.
! locally stored linerecords for lines exiting the node
     type(map_line), dimension(:), allocatable :: linehead
! links to other nodes
! plotlink is used to overlay two or more map or step commands
     type(map_node), pointer :: first,next,previous,plotlink
! saved copy of the meqrec record used to calculate the node 
     type(meq_setup) :: meqrec
! link to saved copy of the equilibrium record
     type(gtp_equilibrium_data), pointer :: nodeceq
! link to array of saved equilibrium record.   (only maptop?)
     type(map_ceqresults), pointer :: saveceq
! copy of nodeceq in saveceq (composition sets not updated but needed for plot)
     integer savednodeceq
! type_of_node not used?? should identify invariants when tie-line not in plane
! lines are number of line records
! noofstph is number of stable phases (copied from meqrec)
! tieline_inplane is 1 if so, 0 if step, -1 if no tie-lines (only maptop)
! number_ofaxis is the number of axis, 1=step;  (only maptop)
     integer type_of_node,lines,noofstph,tieline_inplane,number_ofaxis
! seqx is unique identifier for a map node
! seqy unique identifier for maplines, incremented for each line (only maptop)
     integer seqx,seqy
! nodefix is the phase held fix when calculating node (negative if removed)
     type(gtp_phasetuple) :: nodefix
! Value of T and P, copied from meqrec
     double precision, dimension(2) :: tpval
! chemical potentials, copied from meqrec
     double precision, dimension(:), allocatable :: chempots
! stable phase+compset, copied from meqrec (not used?)
     type(gtp_phasetuple), dimension(:), allocatable :: stable_phases
  end TYPE map_node
!\end{verbatim}
!
!-------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE map_axis
! description of the axis variables used for step/map
! The axis condition in bits and pieces
     integer nterm,istv,iref,iunit
     integer, dimension(:,:), allocatable :: indices
     type(gtp_state_variable), dimension(:), allocatable :: axcond
     double precision, dimension(:), allocatable :: coeffs
! the min, max and increment along the axis
     double precision axmin,axmax,axinc
! more must be initiated to 0, if nonzero replaced by a fixed phase
! seqz is the sequential index of the condition in the list (this is not
! changed if conditions are added (at the end) or deleted (active=1)
! we cannot use a pointer as that depend on the current equilibrium.
     integer more,seqz
! This is the last succesfully calculated axis value
     double precision lastaxval
  end TYPE map_axis
!\end{verbatim}
! decrlared as an array with each axis as one element of the array
!
!-------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE map_ceqresults
! stores calculated equilibrium records
     integer size,free
     TYPE(gtp_equilibrium_data), dimension(:), allocatable :: savedceq
  end TYPE map_ceqresults
!\end{verbatim}
!
!--------------------------------------------------------------
!
!\begin{verbatim}
  TYPE graphics_textlabel
! To put labels on a graph we must store these in a list
     TYPE(graphics_textlabel), pointer :: nexttextlabel
     double precision xpos,ypos
     character*40 textline
  end type graphics_textlabel
!
!\end{verbatim}
!
!------------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE graphics_options
! setting options for the plotting, this replaces most arguments in the call
! to ocplot2(ndx,pltax,filename,maptop,axarr,form)
! ndx is mumber of plot axis, pltax is text with plotaxis variables
! filename is intermediary file (maybe not needed)
! maptop is map_node record with all results
! form is type of output (screen or postscript or gif)
     integer :: status=0,rangedefaults(3)=0,axistype(2)=0
     double precision, dimension(3) :: plotmin,plotmax
     double precision, dimension(3) :: dfltmin,dfltmax
! labeldefaults 0 if default, 1 if text provided in plotlabels
! linetype 0 is black full line, >100 is symbols
     integer :: labeldefaults(3),linetype,tielines=0
! label 1 is heading, 2 is x-axis text, 3 is y-axis text 
     character*64, dimension(3) :: plotlabels
     logical gibbstriangle
! the set key command in GNUPLOT specifies where the line id is written
! it can be on/off, placed inside/outside, left/right/center, top/bottom/center,
! and some more options that may be implemented later ...
     character labelkey*24,appendfile*72
! text label to be written at a given position
     TYPE(graphics_textlabel), pointer :: firsttextlabel
! many more options can easily be added when desired, linetypes etc
  end TYPE graphics_options
!\end{verbatim}
!
! fix status during mapping, normally 2 means fix (not used)
  integer, parameter :: MAPPHASEFIX=3
! maximum number of saved equilibria during step/map
  integer, parameter :: maxsavedceq=2000
! for graphopt status word
! GRKEEP is set if graphics windows kept
  integer, parameter :: GRKEEP=0
!
!-------------------------------------------------
!
CONTAINS

!\begin{verbatim}
  subroutine map_setup(maptop,nax,axarr,seqxyz,starteq)
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
    type(map_fixph), pointer :: mapfix
    double precision starting,finish2,axvalok,dgm,tsave,xxx,yyy,zzz
    integer starttid,endoftime,bytdir,seqz,nrestore,termerr
! save current conditions
    character savedconditions*1024
! for copy of constituions
    double precision, allocatable, dimension(:) :: copyofconst
! inactive are indices of axis conditions inactivated by phases set fixed
! inactive not used ...
    integer iadd,irem,isp,seqx,seqy,mode,halfstep,jj,ij,inactive(4),bytaxis
! inmap=1 turns off converge control of T
    integer, parameter :: inmap=1
    character ch1*1
    logical firststep
!
!    write(*,*)'in map_setup'
    call get_all_conditions(savedconditions,-1,starteq)
    if(gx%bmperr.ne.0) then
       write(kou,*)'Cannot save current conditions'
       savedconditions=' '
!    else
!       write(*,*)'Saved: ',trim(savedconditions)
    endif
    nrestore=0
! first transform start points to start equilibria on zero phase lines
! All axis conditions except one are converted to fix phase conditions 
! (if there is just one axis skip this)
! One or more map_node records are created with mapline records each
    call cpu_time(starting)
    call system_clock(count=starttid)
    inactive=0
!
    if(ocv()) write(*,*)'Entering map_setup',nax
! loop to change all start equilibria to start points
! Store the start points in map_node records started from maptop
100 continue
       ceq=>starteq
!       write(*,*)'calling map_startpoint'
       call map_startpoint(maptop,nax,axarr,seqxyz,inactive,ceq)
!       write(*,*)'back from map_startpoint'
       if(gx%bmperr.ne.0) then
          if(ceq%next.gt.0) then
             write(*,101)ceq%next,gx%bmperr
101          format('Failed calculate a start point: ',i4,i7)
             ceq=>eqlista(ceq%next)
             gx%bmperr=0; goto 100
          endif
       endif
! error if no startpoints 
       if(.not.associated(maptop)) then
          write(*,*)'Cound not find a single start equilibria'
          gx%bmperr=4224; goto 1000
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
!    write(*,*)'savesize: ',size(maptop%saveceq%savedceq)
!-------------------------------
! return here for each new line to be calculated
! NOTE we can start a new thread for each line, when a node is found
! all threads stop.  
! If the node already exists the exit corresponing to the new line removed
! and the thread ends
! If the node is new it is created and exits added and the thread ends.
300 continue
    bytaxis=0
    firststep=.TRUE.
! THREADPROTECTED CALL the map_findline will copy the ceq from mapnode
    if(ocv()) write(*,*)'Looking for a line to calculate'
    call map_findline(maptop,axarr,mapfix,mapline)
    if(gx%bmperr.ne.0) goto 1000
! if no line we are finished!
!    write(*,*)'Back from map_findline 1'
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
       endif
    endif
! Each thread must have separate meqrec and ceq records
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    ceq=>mapline%lineceq
! ?? We may have incompatibility between ceq and meqrec if new compsets added
! maybe meqrec should not be a pointer?
    meqrec=>mapline%meqrec
! No grid minimization and the phr is not deallocated with mode<0
! It is necessary to generate new meqrec for each line as there may be new
! composition sets created in other threads.  But we must also specify 
! phases set fix due to the mapping to replace axis conditions.  
! We must provide an array of phase tuples with fix phases. 
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
!    write(*,*)'Stored: ',size(copyofconst)
305 continue
! to be able to handle problems copy the constitutions!!
!    if(mapline%problems.gt.0) then
!       write(*,*)'problems',mapline%problems,ceq%tpval(1)
!    endif
!    write(*,*)'Calling calceq7 with T=',ceq%tpval(1),mapline%axandir
    call calceq7(mode,meqrec,mapfix,ceq)
!    write(*,*)'Back from calceq7 ',gx%bmperr
    if(gx%bmperr.ne.0) then
! error 4187 is to set T or P to less than 0.1
       if(gx%bmperr.eq.4187) then
          write(*,*)'We are here 1'
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
          gx%bmperr=0
!          if(mapline%firstinc.ne.zero) then
! mapline%firstinc if nonzero the active axis conditioon may already
! been incremented this, go back double the distance in the axis
!             mapline%axandir=-mapline%axandir
!             bytdir=1
!          endif
! first time here bytdir=0,  second time bytdir=1,
! axandir has changed sign and take double negative step
!          bytdir=bytdir+1
! Extract the current value of the axis state variable items using pcond
          jj=abs(mapline%axandir)
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
       write(*,*)'Generating mapline%meqrec failed: ',gx%bmperr
       call map_lineend(mapline,axarr(abs(mapline%axandir))%lastaxval,ceq)
! look for a new line to follow
       goto 300
    endif
!    write(*,*)'back from calceq7'
!--------------------------------
! limit the maximum change in T and P, should be small during step/map
    meqrec%tpmaxdelta(1)=2.0D1
    meqrec%tpmaxdelta(2)=1.0D1
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
    call meq_sameset(irem,iadd,mapline%meqrec,mapline%meqrec%phr,inmap,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error in meq_sameset called from smp',gx%bmperr
!       goto 1000
! check if there are more lines
       gx%bmperr=0
       goto 805
    endif
!    write(*,323)'Calc line: ',gx%bmperr,irem,iadd,mapline%axandir,&
!         mapline%meqrec%noofits,mapline%meqrec%nstph,ceq%tpval(1)
    if(ocv())write(*,323)'Calc line: ',gx%bmperr,irem,iadd,mapline%axandir,&
         mapline%meqrec%noofits,mapline%meqrec%nstph,ceq%tpval(1)
323 format(a,i5,2i3,2i4,i3,f10.2)
! check if it is a closing miscibility gap or loss of ordering
    if(iadd.gt.0) then
       if(same_composition(iadd,mapline%meqrec%phr,mapline%meqrec,ceq,dgm)) &
            iadd=0
    endif
!    write(*,*)'Check if same phase: ',iadd
330 continue
    if(gx%bmperr.eq.0 .and. irem.eq.0 .and. iadd.eq.0) then
! no error and no change of phase set, just store the calculated equilibrium.
!       write(*,*)'hms: Storing equilibrium'
       if(mapline%number_of_equilibria.gt.10 .and. mapline%nodfixph.gt.0) then
! we have managed 3 steps, set phase at start node as entered (if dormant)
          if(meqrec%phr(mapline%nodfixph)%phasestatus.eq.PHDORM) then
!             write(*,*)'Phase set entered ',mapline%nodfixph
             meqrec%phr(mapline%nodfixph)%phasestatus=PHENTUNST
          endif
       endif
       mapline%problems=0
       nrestore=0
!       mapline%lasterr=0
       call map_store(mapline,axarr,nax,maptop%saveceq)
       if(gx%bmperr.ne.0) then
! Test if we are running out of memory 
          if(gx%bmperr.eq.4219) goto 1000
! or that a line has a too big jump
! terminate line any error code will be cleared inside map_lineend. 
!          write(*,*)'Calling map_lineend 1'
          call map_lineend(mapline,axarr(abs(mapline%axandir))%lastaxval,ceq)
          goto 300
       endif
! check which axis variable changes most rapidly, maybe change step axis
! (for tie-lines in plane check axis values for all phases)
! and take a step in this axis variable making sure it inside the limits
! and continue, else terminate and take another start equilibrium
! Normally do not change the phase kept fix.
!       write(*,*)'hms: taking a step'
       call map_step(maptop,mapline,mapline%meqrec,mapline%meqrec%phr,&
            axvalok,nax,axarr,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error return from map_step: ',gx%bmperr
          gx%bmperr=0
          if(meqrec%tpindep(1)) then
             write(*,*)'Restore T 1: ',tsave,axvalok
             ceq%tpval(1)=tsave
          endif
          call map_halfstep(halfstep,axvalok,mapline,axarr,ceq)
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
          write(*,*)'Error stepping to next equilibria, ',gx%bmperr
       endif
! any error code will be cleared inside map_lineend.
!       write(*,*)'Calling map_lineend 2'
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
          call map_halfstep(halfstep,axvalok,mapline,axarr,ceq)
          if(gx%bmperr.eq.0) then
! jump back without setting halfstep=0, setting iadd=-1 turn on debug output 
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
!          write(*,19)'Force changeaxis: ',mapline%axandir,gx%bmperr,axvalok,xxx
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
    if(irem.gt.0 .and. iadd.gt.0) then
! if there is phase which wants to appear and another disappear then
! calculate with half the step 5 times. If axvalok=0 no previous axis value
       write(*,*)'Phases want to appear/disappear: ',iadd,irem
! restore constitutions
!       write(*,*)'Restore constitutions 3',axvalok,ceq%tpval(1)
       call restore_constitutions(ceq,copyofconst)
       call map_halfstep(halfstep,axvalok,mapline,axarr,ceq)
       if(gx%bmperr.eq.0) then
! jump back without setting halfstep=0
          goto 320
       elseif(nax.gt.1 .and. bytaxis.eq.0) then
! try to change axis with active condition.
          if(meqrec%tpindep(1)) then
             if(ocv()) write(*,*)'Restoring T 4: ',tsave,axvalok
             ceq%tpval(1)=tsave
          endif
          if(ocv()) write(*,557)gx%bmperr,ceq%tpval(1),axvalok
557       format('Repeated error 8, try to change axis',i5,F8.2,1pe14.6)
          gx%bmperr=0
          bytaxis=1
          call map_force_changeaxis(maptop,mapline,mapline%meqrec,nax,axarr,&
               axvalok,ceq)
          if(gx%bmperr.eq.0) goto 320
          call map_lineend(mapline,zero,ceq)
       else
! there is an error, take another line
          call map_lineend(mapline,zero,ceq)
       endif
!-----------------------------------------------------
! a new phase stable or a stable wants to disappear
    elseif(irem.gt.0 .or. iadd.gt.0) then
!       write(*,*)'New phase 2: ',iadd,irem,mapline%nodfixph,&
!            mapline%number_of_equilibria
!       if(mapline%number_of_equilibria.lt.2 .and.&
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
          write(*,*)'Startnode phase ignored: ',2,mapline%nodfixph,&
               ceq%tpval(1)
          iadd=0; irem=0
! set the phase dormant and decrease step
          write(*,559)mapline%nodfixph,axvalok
559       format('Phase set dormant ',i5,1pe14.6)
          meqrec%phr(mapline%nodfixph)%phasestatus=PHDORM
          call map_halfstep(halfstep,axvalok,mapline,axarr,ceq)
          if(gx%bmperr.eq.0) then
! jump back without setting halfstep=0
             goto 320
          elseif(nax.gt.1 .and. bytaxis.eq.0) then
! try to change axis with active condition.
             if(meqrec%tpindep(1)) then
                if(ocv()) write(*,*)'Restoring T 7: ',tsave,axvalok
                ceq%tpval(1)=tsave
             endif
             write(*,557)gx%bmperr,ceq%tpval(1),axvalok
             gx%bmperr=0
             bytaxis=1
             call map_force_changeaxis(maptop,mapline,mapline%meqrec,&
                  nax,axarr,axvalok,ceq)
             if(gx%bmperr.eq.0) goto 320
             call map_lineend(mapline,zero,ceq)
          else
! there is a persistent error, take another line
             call map_lineend(mapline,zero,ceq)
          endif
       endif
       if(mapline%more.eq.0) then
! This is the last equilibrium at axis limit
          if(irem.gt.0) then
! terminate the line and check if there are other lines to calculate
             call map_lineend(mapline,zero,ceq)
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
!       write(*,*)'New phase 3: ',iadd,irem,ceq%tpval(1)
! calculate the exact value of the variable axis for the phase change
! then check if we have already found this node point and if not
! generate new start points with and without the phase
       call map_calcnode(irem,iadd,maptop,mapline,mapline%meqrec,axarr,ceq)
       if(gx%bmperr.ne.0) then
! if error one can try to calculate using a shorter step or other things ...
!          write(*,*)'Error return from map_calcnode: ',gx%bmperr
          if(gx%bmperr.eq.4353) then
! this means node point not global, line leading to this is set inactive
! and we should not generate any startpoint.             
             write(*,*)'Setting line inactive',mapline%lineid
             mapline%status=ibset(mapline%status,EXCLUDEDLINE)
             call map_lineend(mapline,zero,ceq)
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
          call map_halfstep(halfstep,axvalok,mapline,axarr,ceq)
          if(gx%bmperr.eq.0) then
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
             call map_lineend(mapline,zero,ceq)
          else
             if(ocv()) write(*,*)' *** Repeated errors calling map_calcnode,',&
                  ' terminate line',gx%bmperr
! terminate line and follow another line, error reset inside map_lineend
             call map_lineend(mapline,zero,ceq)
          endif
       endif
       axvalok=zero
    endif
! we have finished a line
805 continue
    write(*,808)mapline%number_of_equilibria,ceq%tpval(1)
808 format('Finished line with ',i5,' equilibria at T=',F8.2)
    mapline%problems=0
    mapline%lasterr=0
    goto 300
!-----------------------------------------------------
! we come here when there are no more lines to calculate, maybe save on file?
900 continue
!-----------------------------------------------------
1000 continue
!--------------------------------------------------
! Here we have now finished the step/map.
! Set back inactive axis conditions
!    do ij=2,inactive(1)
!       call locate_condition(inactive(ij),pcond,ceq)
!       pcond%active=0
!    enddo
    call system_clock(count=endoftime)
    call cpu_time(finish2)
    if(gx%bmperr.ne.0) then
       write(*,1005)gx%bmperr
1005   format('STEP/MAP terminated with error code: ',i5)
    else
       write(*,1010)finish2-starting,endoftime-starttid
1010   format('Finished step/map ',1pe12.4,' s and ',i7,' clockcycles')
    endif
    return
  end subroutine map_setup

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
! inactive is not really used (conditions replaced by fix phase)
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
    TYPE(map_line), dimension(2) :: tmpline
    TYPE(map_node), pointer :: mapnode,tmpnode
!    type(gtp_phasetuple) :: phfix
    type(map_fixph), pointer :: mapfix
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
    allocate(meqrec)
! We must use mode=-1 for map_replaceaxis below has to calculate several equil
! and the phr array must not be deallocated.  mapfix will be used later to
! indicate fix and stable phases for different lines (maybe ...)
    mode=-1
    nullify(mapfix)
    call calceq7(mode,meqrec,mapfix,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error calling calceq7 in map_startpoint',gx%bmperr
       goto 1000
    endif
! check if equilibrium inside axis limits ...
    do iax=1,nax
       call locate_condition(axarr(iax)%seqz,condition,ceq)
       if(gx%bmperr.ne.0) goto 1000
       call condition_value(1,condition,value,ceq)
       if(gx%bmperr.ne.0) goto 1000
       if(value.lt.axarr(iax)%axmin .or. value.gt.axarr(iax)%axmax) then
          write(*,*)'Startpoint outside axis limits'
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
       tmpnode%next=>maptop%previous
       mapnode=>maptop%previous
       mapnode%noofstph=-1
       mapnode%previous=>tmpnode
       mapnode%next=>maptop
       mapnode%first=>maptop
       mapnode%seqx=tmpnode%seqx+1
!       mapnode%nodefix%phaseix=0
       mapnode%nodefix%ixphase=0
!       write(*,*)'creating another mapnode record',mapnode%seqx
! nullify here indicates more than one node record
       nullify(tmpnode)
    else
! This is the first (and maybe only) mapnode record (later maptop)
!       write(*,*)'Creating first maptop'
       allocate(maptop)
       mapnode=>maptop
       mapnode%noofstph=meqrec%nstph
       mapnode%savednodeceq=-1
!       mapnode%noofstph=-1
       mapnode%next=>mapnode
       mapnode%previous=>mapnode
       mapnode%first=>mapnode
       mapnode%number_ofaxis=nax
!       mapnode%nodefix%phaseix=0
       mapnode%nodefix%ixphase=0
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
       call map_replaceaxis(meqrec,axactive,ieq,nax,axarr,tmpline,inactive,ceq)
       if(ocv()) write(*,205)'Back from replaceaxis with: ',gx%bmperr,&
            axactive,ieq,&
            tmpline(1)%linefixph(1)%ixphase,tmpline(1)%linefixph(1)%compset,&
            tmpline(1)%stableph(1)%ixphase,tmpline(1)%stableph(1)%compset
!            tmpline(1)%linefixph(1)%phaseix,tmpline(1)%linefixph(1)%compset,&
!            tmpline(1)%stableph(1)%phaseix,tmpline(1)%stableph(1)%compset
205    format(a,3i5,5x,2(2i3))
       if(gx%bmperr.ne.0) goto 1000
    else
! only one axis, i.e. a step command, create a map_node record with 2 lines
       axactive=1
       ieq=2
       continue
    endif
!    write(*,1001)'After replace: ',(meqrec%phr(jp)%curd%amfu,&
!         jp=1,meqrec%nphase)
!-----------------------------------------------------------------------
! finished converting a start equilibrium to a start point, 
    mapnode%type_of_node=0
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
    allocate(mapnode%linehead(ieq))
! we can have 3 or more exits if starting inside a 3 phase triagle for isotherm
    if(ieq.eq.2) then
! STEP command: set one exit in each direction of the active axis axactive
       do jp=1,2
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
          mapnode%linehead(jp)%axandir=(3-2*jp)*axactive
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
! tie-lines not in plane, just copied from tielines in plane
             kp=tmpline(1)%nfixphases
             mapnode%linehead(jp)%nfixphases=kp
             allocate(mapnode%linehead(jp)%linefixph(kp))
!             write(*,454)jp,axactive,mapnode%linehead(jp)%axandir,kp
454          format('Axis direction etc: ',i2,2i4,2x,i3)
             do zz=1,kp
                mapnode%linehead(jp)%linefixph(zz)=tmpline(1)%linefixph(zz)
             enddo
! we can have many stable phases
             mapnode%linehead(jp)%nstabph=tmpline(1)%nstabph
             allocate(mapnode%linehead(jp)%stableph(tmpline(1)%nstabph))
             allocate(mapnode%linehead(jp)%stablepham(tmpline(1)%nstabph))
             do kp=1,mapnode%linehead(jp)%nstabph
                mapnode%linehead(jp)%stableph(kp)=tmpline(1)%stableph(kp)
                mapnode%linehead(jp)%stablepham(kp)=tmpline(1)%stablepham(kp)
             enddo
!             write(*,*)'allocated size of stableph 1: ',jp,&
!                  size(mapnode%linehead(jp)%stableph)
             if(ocv())write(*,27)'We have a startpoint for no tie-lines map:',&
!                  axactive,mapnode%linehead(jp)%linefixph(1)%phaseix,&
                  axactive,mapnode%linehead(jp)%linefixph(1)%ixphase,&
                  mapnode%linehead(jp)%linefixph(1)%compset,&
                  mapnode%linehead(jp)%nstabph,&
!                  (mapnode%linehead(jp)%stableph(kp)%phaseix,&
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
             do zz=1,kp
                mapnode%linehead(jp)%linefixph(zz)=tmpline(1)%linefixph(zz)
             enddo
! there is just one stable phase
             allocate(mapnode%linehead(jp)%stableph(1))
             mapnode%linehead(jp)%nstabph=1
             mapnode%linehead(jp)%stableph=tmpline(1)%stableph
             if(ocv()) write(*,25)'We have saved a startpoint for map:',&
!                  axactive,mapnode%linehead(jp)%linefixph(1)%phaseix,&
                  axactive,mapnode%linehead(jp)%linefixph(1)%ixphase,&
                  mapnode%linehead(jp)%linefixph(1)%compset,&
                  mapnode%linehead(jp)%nstabph,&
!                  mapnode%linehead(jp)%stableph(1)%phaseix,&
                  mapnode%linehead(jp)%stableph(1)%ixphase,&
                  mapnode%linehead(jp)%stableph(1)%compset
25           format(a,i3,5x,2i3,5x,i3,2x,2i3)
!------------------------- below for STEP
          else
! this is for STEP
             if(ocv()) write(*,*)'For STEP no need of fixed phases.'
             mapnode%linehead(jp)%nfixphases=0
             allocate(mapnode%linehead(jp)%stableph(meqrec%nstph))
             mapnode%linehead(jp)%nstabph=meqrec%nstph
             do kp=1,mapnode%linehead(jp)%nstabph
                zz=meqrec%stphl(kp)
!                mapnode%linehead(jp)%stableph(kp)%phaseix=meqrec%phr(zz)%iph
                mapnode%linehead(jp)%stableph(kp)%ixphase=meqrec%phr(zz)%iph
                mapnode%linehead(jp)%stableph(kp)%compset=meqrec%phr(zz)%ics
             enddo
!             mapnode%linehead(jp)%nstabph
!             mapnode%linehead(jp)%stableph(1)%phase
!             mapnode%linehead(jp)%stableph(1)%compset
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
!    seqx=maptop%next%seqx+1
!    write(*,*)'Maptop index: ',maptop%next%seqx,maptop%previous%seqx
    seqx=max(maptop%next%seqx,maptop%previous%seqx)+1
!    write(*,666)seqx,maptop%seqx,maptop%next%seqx,maptop%previous%seqx
666 format('maptop seqx: ',10i3)
    call wriint(eqname,jp,seqx)
! NOT A GOOD IDEA TO SET PHASE FIX HERE ...
!    write(*,*)'stp MAP Copying equilibrium record: ',ceq%eqname,ceq%eqno
!    write(*,*)'stp to new  equilibrium record: ',eqname,jp
!    write(*,*)'stp fixph: ',mapnode%linehead(1)%linefixph(1)%phaseix,&
!         mapnode%linehead(1)%linefixph(1)%compset
! list conditions, fix phases due to mapping missing
!    call list_conditions(kou,ceq)
! add the fix phase, try to use status 3 for a mapfix phase
!    call change_phase_status(mapnode%linehead(1)%linefixph(1)%phaseix,&
!         mapnode%linehead(1)%linefixph(1)%compset,MAPPHASEFIX,zero,ceq)
! list conditions, fix phases due to mapping missing
!    write(*,*)'Maybe now the phase is fix?'
!    call list_conditions(kou,ceq)
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
! Copy the current meqrec to mapnode, the mapline records
! will generate their own new meqrec records when they are activated
! if the phr array is allocated then deallocate it as it is no longer needed
    if(allocated(meqrec%phr)) then
       deallocate(meqrec%phr)
    endif
    mapnode%meqrec=meqrec
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

!\begin{verbatim}
  subroutine map_replaceaxis(meqrec,axactive,ieq,nax,axarr,tmpline,&
       inactive,ceq)
! replace an axis condition with a fix phase
! meqrec is equilibrium calculation record
! axactive is the axis with active condition, ieq is number of exiting lines
! ieq is the number of lines exiting from the startpoint
! nax is number of axis, axarr are description of the axis
! axarr is array with axis records
! tmpline is to transfer some line data to calling routine
! inactive is not really used.
! ceq is equilibrium record
    implicit none
    type(meq_setup), pointer :: meqrec
    integer nax,axactive,ieq
    type(map_axis), dimension(nax) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
    type(map_line), dimension(2) :: tmpline
    integer inactive(*)
!\end{verbatim}
    integer iph,jph,naxvar,iax,tip,jj,jax,irem,iadd,kj,nrel
    integer ics,lokph,lokcs,kph,kcs
    double precision aval,avalm
    type(gtp_condition), pointer :: pcond
    integer, dimension(:), allocatable :: axis_withnocond
! turns off converge control for T
    integer, parameter :: inmap=1
!
    nrel=noel()
    tip=tieline_inplane(nax,axarr,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(ocv()) write(*,*)'In map_replaceaxis ',tip
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
! zero the number of fix phases and allocate data for those needed
    tmpline%nfixphases=0
    allocate(tmpline(1)%linefixph(naxvar-1))
!========================================================== tie-lines in plane
    tieline_in_plane: if(tip.eq.1) then
! We have tie-lines in the plane, only one stable phase in addition to fix
       allocate(tmpline(1)%stableph(1))
       allocate(tmpline(1)%stablepham(1))
       allocate(axis_withnocond(nax))
       axis_withnocond=0
       stablephases: if(meqrec%nstph.gt.1) then
! and we have two or more stable phase, we can directly generate startpoints
100       continue
          write(*,*)'Tie-lines in the plane and start equilibrium with',&
               ' several stable phases'
          jax=0
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
             aval=meqrec%phr(jj)%curd%amfu
             if(aval.lt.avalm) then
                kph=meqrec%phr(jj)%iph
                kcs=meqrec%phr(jj)%ics
                kj=jj
                avalm=aval
             endif
          enddo
! The phase meqrec%phr(kj)%iph/ics should be set fix 
          if(ocv()) write(*,*)'Fixing phase: ',kj,&
               meqrec%phr(kj)%iph,meqrec%phr(kj)%ics
          meqrec%phr(kj)%curd%dgm=zero
! 
!          if(tieline_inplane(nax,axarr,ceq).eq.1) then
! if tie-lines in the plane set amount 0.1
!             write(*,*)'Setting fix phase nonzero amount'
             meqrec%phr(kj)%curd%amfu=1.0D-1
!          else
             meqrec%phr(kj)%curd%amfu=zero
!          endif
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
          if(ocv()) write(*,77)pcond%seqz,pcond%prescribed,ceq%tpval(1)
77        format(' Removing condition: ',i3,2(1pe12.4))
! We have tried not to replace T or P,
! if this is done it must be indicated specially like this
          if(pcond%statev.eq.1) then
             meqrec%tpindep(1)=.TRUE.
             if(ocv()) write(*,*)'Marking that T is variable'
          elseif(pcond%statev.eq.2) then
             meqrec%tpindep(2)=.TRUE.
          endif
!---------------------------------------------------
! calculate the equilibrium with the new set of conditions
          if(ocv()) write(*,*)'Calling meq_sameset inside  map_replaceaxis'
          irem=0; iadd=0;
          call meq_sameset(irem,iadd,meqrec,meqrec%phr,inmap,ceq)
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
!             tmpline(1)%linefixph%phaseix=kph
             tmpline(1)%linefixph%ixphase=kph
             tmpline(1)%linefixph%compset=kcs
             tmpline(1)%nstabph=0
             do jph=1,meqrec%nstph
                jj=meqrec%stphl(jph)
                if(meqrec%phr(jj)%iph.eq.kph .and.&
                     meqrec%phr(jj)%ics.eq.kcs) cycle
!                tmpline(1)%stableph(1)%phaseix=meqrec%phr(jj)%iph
                tmpline(1)%stableph(1)%ixphase=meqrec%phr(jj)%iph
                tmpline(1)%stableph(1)%compset=meqrec%phr(jj)%ics
                tmpline(1)%nstabph=tmpline(1)%nstabph+1
!                tmpline(1)%nstabph=1
! why exit??
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
          call map_startline(meqrec,axactive,ieq,nax,axarr,tmpline,ceq)
!          write(*,*)'Tie-line in plane and single phase,',&
!               ' no check yet if this works ... '
       endif stablephases
! ============================================= no tie-lines in plane
    else !tie-lines NOT in the plane
! I am not sure what stableph and axis_withnocond are used for ...
!       allocate(tmpline(1)%stableph(1)) this is allocated in map_startline
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
    return
  end subroutine map_replaceaxis

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
    integer jax,iax,idir,irem,iadd,iph,jj,jph,kph,ll
    integer, parameter :: nstabphdim=20
    double precision curval,startval
    type(gtp_condition), pointer :: pcond
! turns off converge control for T
    integer, parameter :: inmap=1
!
!    write(*,*)'In map_startline, find a phase to set fix'
! start in negative direction unless direction given
    idir=-1
    if(ceq%multiuse.ne.0) then
       if(abs(ceq%multiuse).gt.nax) then
          write(*,*)'Error in direction, no such axis, ',ceq%multiuse
          gx%bmperr=4229; goto 1000
       endif
       if(jax.gt.0) idir=1
       jax=abs(ceq%multiuse)
       call locate_condition(axarr(jax)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
    else
       jax=0
       do iax=1,nax
          call locate_condition(axarr(iax)%seqz,pcond,ceq)
          if(gx%bmperr.ne.0) goto 1000
          if(pcond%statev.lt.10) then
! this means intensive variable (T,P chemical potential)
             jax=iax; exit
          endif
       enddo
       if(jax.eq.0) jax=1
       call locate_condition(axarr(jax)%seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'Searching for phase to fix along axis: ',jax
    endif
    call condition_value(1,pcond,curval,ceq)
    if(gx%bmperr.ne.0) goto 1000
    startval=curval
! increment axis variable using axinc and calculate with meq_sameset
100 continue
       curval=curval+idir*axarr(jax)%axinc
       call condition_value(0,pcond,curval,ceq)
       if(gx%bmperr.ne.0) goto 1000
       irem=0; iadd=0; meqrec%noofits=0
       call meq_sameset(irem,iadd,meqrec,meqrec%phr,inmap,ceq)
!       if(ocv()) write(*,110)'Search for phase change: ',&
!       write(*,110)'Search for phase change: ',&
!            idir*jax,gx%bmperr,irem,iadd,ceq%tpval(1),curval
110    format(a,i2,3i5,2x,F8.2,1pe14.6)
       if(gx%bmperr.ne.0) goto 1000
       nophasechange: if(irem.eq.0 .and. iadd.eq.0) then
          if(idir.lt.0) then
             if(curval.le.axarr(jax)%axmin) then
! change direction
                idir=+1
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
    fixfas: if(irem.gt.0) then
       if(meqrec%nstph.eq.1) then
          write(*,*)'Attempt to set the only phase as fix!'
          gx%bmperr=4230; goto 1000
       endif
       write(*,*)'Remove axis condition and set stable phase fix: ',irem
! phase already in lists, just mark it is no fixed with zero amount
!       meqrec%phr(irem)%itrem=meqrec%noofits
!       meqrec%phr(irem)%prevam=zero
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
!       write(*,*)'Remove axis condition and set new phase fix: ',iadd
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
!    write(*,*)'No axis condition and call meq_sameset with fix phase: ',kph
! Must not forget to set if T or P is variable!
    pcond%active=1
    if(pcond%statev.eq.1) then
       meqrec%tpindep(1)=.TRUE.
    elseif(pcond%statev.eq.2) then
       meqrec%tpindep(2)=.TRUE.
    endif
! calling meq_sameset with iadd=-1 turn on verbose
    irem=0; iadd=0
    call meq_sameset(irem,iadd,meqrec,meqrec%phr,inmap,ceq)
!    if(ocv()) write(*,110)'meq_sameset calculated: ',&
    if(gx%bmperr.gt.0) then
       write(*,*)'Filed to calculate with fix phase',gx%bmperr
       goto 1000
    elseif(iadd.gt.0 .or. irem.gt.0) then
       write(*,*)'Another phase want to be stable: ',iadd,irem
       gx%bmperr=4232; goto 1000
    endif
!    write(*,110)'meq_sameset calculated: ',0,gx%bmperr,irem,iadd,ceq%tpval(1)
!    if(gx%bmperr.ne.0) goto 1000
!
!
! we must return some values
! two exits
    ieq=2
! active axis, the remaining one, if jax=1 then 2, if jax=2 then 1
    axactive=3-jax
! templine is map_line record, some data must be set
    tmpline(1)%nfixphases=1
!    tmpline(1)%linefixph%phaseix=meqrec%phr(kph)%iph
    tmpline(1)%linefixph%ixphase=meqrec%phr(kph)%iph
    tmpline(1)%linefixph%compset=meqrec%phr(kph)%ics
! allocate space for all stable phases minus one as fix, may already be alloc
!    allocate(tmpline(1)%stableph(meqrec%nstph-1))
!    allocate(tmpline(1)%stablepham(meqrec%nstph-1))
! The number of stable phases can vary for different MAP commands
    if(allocated(tmpline(1)%stableph)) then
       deallocate(tmpline(1)%stableph)
       deallocate(tmpline(1)%stablepham)
    endif
    allocate(tmpline(1)%stableph(nstabphdim))
    allocate(tmpline(1)%stablepham(nstabphdim))
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
!       tmpline(1)%stableph(ll)%phaseix=meqrec%phr(jj)%iph
       tmpline(1)%stableph(ll)%ixphase=meqrec%phr(jj)%iph
       tmpline(1)%stableph(ll)%compset=meqrec%phr(jj)%ics
       tmpline(1)%stablepham(ll)=meqrec%phr(jj)%curd%amfu
       tmpline(1)%nstabph=tmpline(1)%nstabph+1
! why exit?
!       exit
    enddo
!    if(ocv()) write(*,300)axactive,kph,tmpline(1)%linefixph%phaseix,&
    if(ocv()) write(*,300)axactive,kph,tmpline(1)%linefixph%ixphase,&
         tmpline(1)%linefixph%compset,tmpline(1)%nstabph,&
!         (tmpline(1)%stableph(jj)%phaseix,tmpline(1)%stableph(jj)%compset,&
         (tmpline(1)%stableph(jj)%ixphase,tmpline(1)%stableph(jj)%compset,&
         jj=1,tmpline(1)%nstabph)
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
    integer place,jph,jj
    type(meq_setup), pointer :: meqrec
    type(gtp_state_variable), target :: axstv1
    type(gtp_state_variable), pointer :: axstv
    double precision value
    logical saveonfile
! pointer to last calculated (can be zero) and last free
! store last calulated axis values in axarr(iax)%lastaxval
!    write(*,*)'In map_store',mapline%start%number_ofaxis,nax
    do jj=1,nax
       axstv1=axarr(jj)%axcond(1)
       axstv=>axstv1
       call state_variable_val(axstv,value,mapline%lineceq)
       if(gx%bmperr.gt.0) goto 1000
       if(nax.gt.1) then
! when several axis check if any has a big change ...
          if(mapline%number_of_equilibria.gt.3) then
             if(abs(axarr(jj)%lastaxval-value).gt.&
                  1.0D-1*(axarr(jj)%axmax-axarr(jj)%axmin)) then
!                  2.0D-2*(axarr(jj)%axmax-axarr(jj)%axmin)) then
                write(*,17)'map_step found too large step in axis: ',jj,&
                     mapline%number_of_equilibria,axarr(jj)%lastaxval,value
17              format(a,i2,i3,2(1pe14.6))
                gx%bmperr=4300; goto 1000
             endif
          endif
       endif
       axarr(jj)%lastaxval=value
    enddo
!    write(*,18)'stored: ',mapline%number_of_equilibria,(axarr(jj)%lastaxval,&
!         jj=1,mapline%start%number_ofaxis)
!18  format(a,i3,5(1pe14.6))
!-----------------------
    saveonfile=.FALSE.
! >>>> begin treadprotected
    call reserve_saveceq(place,saveceq)
    if(gx%bmperr.eq.4219) then
! the memory is full, save this equilibrium, clean up and empty all on file
       saveonfile=.TRUE.
       gx%bmperr=0
    elseif(gx%bmperr.ne.0) then
! some other fatal error
       goto 1000
    endif
! >>>> end threadprotected
!-----------------------
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
!       else
!          write(*,*)'Fix phase 1: ',jj,meqrec%phr(jj)%iph,meqrec%phr(jj)%ics
       endif
    enddo
!    write(*,201)' map store, stable phases: ',&
!         (meqrec%phr(meqrec%stphl(jj))%iph,&
!         meqrec%phr(meqrec%stphl(jj))%ics,jj=1,meqrec%nstph)
201 format(a,10(2i3,2x))
!-----------------------------------------
! this copies the whole data structure !!!
! LIKELY PLACE FOR SEGMENTATION FAULT !!!
!    write(*,*)'SMP storing equilibrium record: ',place
    saveceq%savedceq(place)=mapline%lineceq
    saveceq%savedceq(place)%next=0
    if(mapline%last.gt.0) then
       saveceq%savedceq(mapline%last)%next=place
    endif
    mapline%last=place
    mapline%number_of_equilibria=mapline%number_of_equilibria+1
    if(mapline%first.eq.0) mapline%first=place
    if(saveonfile) then
! We have to wind up all unfinished lines to continue step/map
! but this is not yet implemented
       write(*,207)
207    format(/' *** Buffer full and save on file not yet implemented,',&
            ' terminating step/map'/)
       gx%bmperr=4219
    endif
1000 continue
    return
  end subroutine map_store

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
! maybe do some cleanup
    else
       write(kou,77)mapline%number_of_equilibria,value
77     format('Terminating line with ',i4,' equilibria at axis limit ',1pe12.4)
       mapline%termerr=0
    endif
! mark there is no node at the end
    nullify(mapline%end)
! Not much more to do I guess??
!    write(*,*)'Terminated line: ',mapline%lineid,mapline%number_of_equilibria
1000 continue
! This reoutine should never return an error code
    gx%bmperr=0
    return
  end subroutine map_lineend

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine map_changeaxis(mapline,nyax,oldax,nax,axarr,axval,bytax,ceq)
! changes the axis with active condition to nyax
! mapline is line record
! nyax is index of new active axis
! oldax is index of old active axis
! nax is number of axis (always 2?)
! axarr is array with axis records
! axval the value to set as condition on new axis
! bytax logical, if true ignore axval
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
    integer jax,iadd,irem,ierr
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
       mapline%meqrec%tpindep(1)=.FALSE.
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
       mapline%meqrec%tpindep(1)=.TRUE.
    endif
!----------------------------------------------------------
! now we calculate the same equilibrium but with different axis condition!
    irem=0
! add=-1 turn on verbose in meq_sameset
!    iadd=-1
    if(ocv()) write(*,*)'Map_changeaxis call meq_sameset, T=',ceq%tpval(1)
    call meq_sameset(irem,iadd,mapline%meqrec,mapline%meqrec%phr,inmap,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error from meq_sameset when trying to change axis',gx%bmperr
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
       mapline%axandir=-nyax
!       xxx=mapline%axvals(nyax)-1.0D-2*axarr(nyax)%axinc
    else
! set positive direction
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

!\begin{verbatim}
  subroutine map_step(maptop,mapline,meqrec,phr,axvalok,nax,axarr,ceq)
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
!    write(*,*)'In map_step: ',mapline%number_of_equilibria
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
!                svrtarget%phase=mapline%linefixph(1)%phaseix
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
                if(ocv()) write(*,99)mapline%number_of_equilibria,jax,jaxwc,&
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
!                write(*,101)mapline%number_of_equilibria,jaxwc,&
!                     nyax,xxx,dax1(1),dax2(1),dax1(2),dax2(2)
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
!                   write(*,*)'Change active axis: ',jax
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
! Trying to male axfact decrease less (like line above) makes map worse
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
23        format('At axis ',a,' limit',1pe12.4)
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
  end subroutine map_step

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
    integer what,type,cmix(10),maxstph,noplot,mode
    double precision, parameter :: addedphase_amount=1.0D-2
    double precision value,axval,axvalsave
    type(gtp_state_variable), pointer :: svrrec
    logical global
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
! axandir is the avis with active condition.  It can be negative
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
       meqrec%tpindep(1)=.TRUE.
       maxstph=1
    endif
!--------------------------------------------
! independently if iadd or irem >0 set this phase, phfix, fix with zero amount
! we may return here if there is problems calculate the node equilibrium
100 continue
! set phfix fix with zero amount.
! NOTE it must be added so meqrec%stphl in ascending order
    if(phfix.gt.0 .and. meqrec%nstph.eq.meqrec%maxsph+maxstph) then
! No more phases allowed, we must see if  some other phase may be removed
!       write(*,*)'Too many stable phases at nodepoint',meqrec%maxsph
! set back pcond as active, this saved top of miscibility gap in Cr-Mo !!!
       pcond%active=0
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
       meqrec%phr(phfix)%curd%amfu=zero
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
       meqrec%phr(-phfix)%curd%amfu=zero
       meqrec%phr(-phfix)%stable=1
! set that the phase has fixed amount
       meqrec%phr(-phfix)%phasestatus=PHFIXED       
    endif
!--------------
! mark that the phase is fix, we have to be careful not to exceed size
! Sigh, the fixed phases must be in sequential order ??? ... not done here
! ... maybe not needed ??
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
200 continue
    iadd=0; irem=0
!    write(*,*)'In map_calcnode calling sameset for new node: ',&
!         meqrec%nstph,meqrec%nfixph
!
    call meq_sameset(irem,iadd,meqrec,meqrec%phr,inmap,ceq)
!
!    write(*,202)'Calculated with fix phase: ',gx%bmperr,irem,iadd,ceq%tpval
202 format(a,3i4,2(1pe12.4))
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
       if(ocv()) write(*,*)'Another phase wants to be stable',iadd
       gx%bmperr=4223
    else
! It worked to calculate with a new fix phase releasing all axis condition!!!
! *************************************************************
! check that the node point is global using grid minimizer
! ceq is not destroyed
       mode=0
       global=global_equil_check1(mode,ceq)
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
       meqrec%tpindep(1)=.FALSE.
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
!    write(*,*)'lasterr: ',mapline%lasterr,gx%bmperr
    goto 1000
!------------------------------------------------------
! When we are there we have successfully calculated an equilibrium with a
! new phase set create a node with this equilibrium and necessary line records
500 continue
! phfix is set negative if phase should be removed
! NOTE the phase set fix in the node may not be the same which
! wanted to disappear/appear when calling the map_calcnode!!
! If iremsave=phfix the fix phase is one to be removed.
    if(iremsave.eq.-phfix) then
       phfix=-abs(phfix)
    endif
! if the user wants to have global minimization during mapping this is
! time to test if the current equilibrium is the global one.  We can use
! a temporary ceq record and chech the set of phases and chemical potentials
!
! >>> add test of the global equilibrium here ... it has been added elsewhere
!
    haha=0
    if(maptop%tieline_inplane.lt.0) then
! test if invariant ...
       if(inveq(haha,ceq)) then
          write(*,*)'Invariant equilibrium ignored',haha
! haha is set to number of stable phases at invariant.
! the number of lines ending at this is 2*haha
       endif
    endif
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
       meqrec%tpindep(1)=.FALSE.
       ceq%tpval(2)=value
    endif
! Save this as the last equilibrium of the line
    if(maptop%tieline_inplane.gt.0) then
! remove phfix as fix, otherwise graphics will be strange!
!       write(*,517)'In map_calcnode, mapline%meqrec%nstph 1: ',phfix,&
!            mapline%meqrec%nstph,&
!            (mapline%meqrec%phr(mapline%meqrec%stphl(jj))%iph,&
!            mapline%meqrec%phr(mapline%meqrec%stphl(jj))%ics,&
!            jj=1,mapline%meqrec%nstph)
517    format(a,2i3,5x,5(2i3))
! remove phfix
!       do jj=1,mapline%meqrec%nstph-1
!       do jj=1,mapline%meqrec%nstph
!          if(mapline%meqrec%stphl(jj).ge.phfix) then
! ?? what is done here
!             mapline%meqrec%stphl(jj)=mapline%meqrec%stphl(jj+1)
! this is necessary not to have data from this phase interfering with the line
!             write(*,519)mapline%meqrec%phr(jj)%iph,&
!                  mapline%meqrec%phr(jj)%ics,phentunst
!519          format('Removing ',2i3,' as stable as last line node',i3)
!             mapline%meqrec%phr(jj)%curd%phstate=PHENTUNST
!          endif
!       enddo
!????????????????????????????????????????
! remove phfix
!       mapline%meqrec%phr(phfix)%curd%phstate=PHENTUNST
       mapline%meqrec%phr(phfix)%curd%phstate=PHENTERED
! this is necessary not to have data from this phase interfering with the line
       if(ocv()) write(*,519)phfix,mapline%meqrec%phr(phfix)%iph,&
            mapline%meqrec%phr(phfix)%ics,phentunst
519    format('Removing ',i3,2x,2i3,' as stable as last line equil',i3)
!?????????????????????????????????????????
       mapline%meqrec%nstph=mapline%meqrec%nstph-1
!       write(*,517)'In map_calcnode, mapline%meqrec%nstph 2: ',phfix,&
!            mapline%meqrec%nstph,&
!            (mapline%meqrec%phr(mapline%meqrec%stphl(jj))%iph,&
!            mapline%meqrec%phr(mapline%meqrec%stphl(jj))%ics,&
!            jj=1,mapline%meqrec%nstph)
    endif
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
! Finally create the new node and with new lines
!    write(*,*)'calling map_newnode: ',mapline%meqrec%nfixph,meqrec%nfixph
    call map_newnode(mapline,meqrec,maptop,axval,jax,axarr,phfix,haha,ceq)
    if(gx%bmperr.ne.0) then
       if(ocv()) write(*,*)'Error return from map_newnode: ',gx%bmperr
    endif
! all done??
1000 continue
    return
  end subroutine map_calcnode

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine map_newnode(mapline,meqrec,maptop,axval,lastax,axarr,&
       phfix,haha,ceq)
! must be partially THREADPROTECTED
! first check if a node with this equilibrium already exists
! if not add a new node with appropriate lineheads and arrange all links
! Take care if tie-lines in the plane all lines do not have to be calculated
! NOTE: meqrec not the same as mapline%meqrec !! ??
! mapline is line record for current line
! meqrec has information about last calculated equilibrium
! maptop is node record
! axval is the axis value attemped to calculate when phase set wanted to change
! lastax is index of last active axis
! axarr are axis records
! phfix is phase which is set fix at node point
! haha is nonzero if the calculated equilibrium is invariant
! ceq is equilibrium record
    implicit none
    type(map_node), pointer :: maptop
    type(meq_setup) :: meqrec
    type(map_line), pointer :: mapline,nodexit
    type(map_axis), dimension(*) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
    integer phfix,lastax,haha
    double precision axval
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: newceq,tmpceq
    type(map_node), pointer :: mapnode,newnode
    type(map_line), pointer :: linenode
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: svrrec,svr2
    type(gtp_state_variable), target :: svrtarget
    integer remph,addph,nel,iph,ics,jj,seqx,nrel,jphr,stabph,kph,kcs,kk,lfix
    integer zph,stepax,kpos,seqy,jp,nopotax,lokcs,lokph
! there should be 8 significant digits, first step factor
    double precision, parameter :: vz=1.0D-9,axinc1=1.0D-3
    character eqname*24
    double precision stepaxval,middle,testv,xxx
    lfix=0
! the phase kept fix with zero amount at the node is phfix  It can be
! negative at STEP if it is a phase that will dissapear.
!    write(*,*)'Found a node point with phase change',phfix,ceq%tpval(1)
    if(ocv()) write(*,87)'We are in map_newnode with a fix phase: ',&
         phfix,ceq%tpval(1)
87  format(a,i4,2x,1pe12.4)
!    write(*,*)'We have access to phr: ',meqrec%phr(abs(phfix))%iph,&
!         meqrec%phr(abs(phfix))%ics
! mapnode should be set to point at maptop
    mapnode=>maptop
    nrel=meqrec%nrel
100 continue
! loop all mapnodes to check if any has the same chemical potentials
!       write(*,*)'Comparing chemical potentials with node: ',mapnode%seqx,nrel
!       write(*,105)'T diff: ',ceq%tpval(1),mapnode%tpval(1),&
!            abs(ceq%tpval(1)-mapnode%tpval(1)),abs(vz*mapnode%tpval(1))
!       write(*,105)'P diff: ',ceq%tpval(2),mapnode%tpval(2),&
!            abs(ceq%tpval(2)-mapnode%tpval(2)),abs(vz*mapnode%tpval(2))
       if(abs(ceq%tpval(1)-mapnode%tpval(1)).gt.abs(vz*mapnode%tpval(1)) .or.&
            abs(ceq%tpval(2)-mapnode%tpval(2)).gt.abs(vz*mapnode%tpval(2))) &
            goto 110
       do nel=1,nrel
!          write(*,105)'Chempot: ',ceq%cmuval(nel),mapnode%chempots(nel),&
!               abs(ceq%cmuval(nel)-mapnode%chempots(nel)),&
!               abs(vz*mapnode%chempots(nel))
105       format(a,5(1pe16.8))
          if(abs(ceq%cmuval(nel)-mapnode%chempots(nel)).gt.&
               abs(vz*mapnode%chempots(nel))) goto 110
       enddo
! We can come here with a STEP command without any fix phases
       if(maptop%tieline_inplane.eq.0) then
!          write(*,*)'Step command'
          goto 800
       endif
! T, P and all chemical potentials the same, one should maybe check phases??
!       iph=mapline%linefixph(1)%phaseix
       iph=mapline%linefixph(1)%ixphase
       ics=mapline%linefixph(1)%compset
       if(ocv()) write(*,107)'Node exist: ',&
            mapnode%seqx,size(mapnode%linehead),iph,ics
107    format(a,i5,i3,i5,i2)
       removexit: do jj=1,size(mapnode%linehead)
! loop for all exits
          nodexit=>mapnode%linehead(jj)
          if(ocv()) write(*,108)'Exit: ',jj,nodexit%done,&
               nodexit%linefixph(1)%ixphase,nodexit%linefixph(1)%compset
!               nodexit%linefixph(1)%phaseix,nodexit%linefixph(1)%compset
108       format(a,i4,i7,i5,i2)
          if(nodexit%done.le.0) cycle
!          if(nodexit%linefixph(1)%phaseix.eq.iph .and. &
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
!---------------------------------------------------
    if(maptop%number_ofaxis.gt.1) then
       if(ocv()) write(*,*)' *** map_newnode not finished handling mapping ...'
    endif
!---------------------------------------------------
! We have to create a new mapnode, add as next to maptop
! as I have understood it it does not work to write ... but maybe ...
!    allocate(newnode)
!    newnode%next=>maptop%next
!    maptop%next=>newnode
!    newnode%previous=>maptop
!    newnode%next%previous=>newnode
! so I have written it as ...
120 continue
    mapnode=>maptop%next
    seqx=mapnode%seqx+1
! if maptop%next is maptop do not nullify this pointer !!
! Always add the new record as the next link to maptop
    if(associated(mapnode,maptop)) then
! a single maptop record
       allocate(mapnode%next)
    else
! there is more mapnode records ... 
       allocate(maptop%next)
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
!    write(*,*)'copy equilibrium: ',seqx,nrel
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
! save a copy of ceq also in result (reserve is threadprotected)
    if(ocv()) write(*,*)'Copies node ceq to saveceq'
    call reserve_saveceq(jj,maptop%saveceq)
    if(gx%bmperr.ne.0) goto 1000
    maptop%saveceq%savedceq(jj)=newceq
    newnode%savednodeceq=jj
!    write(*,*)'Copy successful'
!    write(*,*)'Before copying meqrec: ',mapline%meqrec%nfixph,meqrec%nfixph
! maybe it is not necessary to save meqrec and chemical potentials??
    newnode%meqrec=meqrec
!    write(*,*)'New node index: ',newnode%seqx
    allocate(newnode%chempots(nrel))
    newnode%chempots=ceq%cmuval
    newnode%tpval=ceq%tpval
    newnode%type_of_node=0
! correct value of lines will be set later
    newnode%lines=0
    newnode%tieline_inplane=maptop%tieline_inplane
! this seems to be wrong, maptop%number_ofaxis is zero when step separate
    newnode%number_ofaxis=maptop%number_ofaxis
! save index of the phase set fix at the node
!    write(*,*)'Saving index of fix phase: ',phfix
    if(phfix.lt.0) then
!       newnode%nodefix%phaseix=-meqrec%phr(abs(phfix))%iph
       newnode%nodefix%ixphase=-meqrec%phr(abs(phfix))%iph
    else
!       newnode%nodefix%phaseix=meqrec%phr(abs(phfix))%iph
       newnode%nodefix%ixphase=meqrec%phr(abs(phfix))%iph
    endif
    newnode%nodefix%compset=meqrec%phr(abs(phfix))%ics
!    write(*,*)'Saved node fix phase: ',newnode%nodefix%phase,&
!         newnode%nodefix%compset
! the set of stable phases
    newnode%noofstph=meqrec%nstph
    allocate(newnode%stable_phases(newnode%noofstph))
    do jj=1,newnode%noofstph
!       newnode%stable_phases(jj)%phaseix=meqrec%iphl(jj)
       newnode%stable_phases(jj)%ixphase=meqrec%iphl(jj)
       newnode%stable_phases(jj)%compset=meqrec%icsl(jj)
    enddo
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
         mapline%meqrec%nfixph,meqrec%nfixph
    if(maptop%tieline_inplane.eq.0) then
! this is a step, just one line and one exit with the new set of stable phases
       newnode%lines=1
       if(noel().eq.1) then
! step with single phase: problems with phase change as old phase still stable
          call get_phase_compset(phfix,1,lokph,lokcs)
! this change will reomve the previously stable phase in newnode and 
! below also in meqrec
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
! mapping with tie-lines in plane. Always 3 lines ... 2 new exits ??
! depends on number of axis, for 2 axis OK, for 3 axis (one pot) 4 lines meet
       newnode%lines=2
    elseif(haha.gt.0) then
! for mapping without tie-lines in plane haha is nonzero if we are at
! an invariant equlibrium with haha stable phases.
!       write(*,*)'For the moment ignore this'
!       newnode%lines=2*jj-1
       newnode%lines=3
    else
! mapping without tie-lines in plane
! at other node points 4 lines meets, 3 exits
!       write(*,*)'Invariant equilibrium found, exit lines: ',newnode%lines
       newnode%lines=3
    endif
! set link to end node in maplie
    mapline%end=>newnode
!=============================================================================
! we have created sufficient linehead records, now initiate their content
! different depening on STEP (case 1), MAP with tie-lines in plane (case 2)
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
    allocate(newnode%linehead(newnode%lines))
    do jp=1,newnode%lines
!--------------------- code moved from map_findline
! COPY of the equilibrium record from newnode to newnode%linehead(jp)%lineceq
       if(ocv()) write(*,*)'We found a line from node: ',mapnode%seqx
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
    select case(newnode%lines)
!==========================================================================
    case default
       write(*,*)'*** Trying to create node with lines: ',newnode%lines
       gx%bmperr=4237; goto 1000
!==========================================================================
    case(1)! step node with just one exit
! If phfix negative the fix phase wants to dissapear
       changephaseset: if(phfix.lt.0) then
! remove a phase ---------------------------
          remph=-phfix
!          write(kou,88)remph,' disappears,',meqrec%nstph
88        format('SMP a node created where phase ',i3,a,' stable phases:',i3)
          if(meqrec%nstph.eq.1) then
             write(*,*)'Attempt to remove the only stable phase!!!'
             gx%bmperr=4238; goto 1000
          endif
! shift phases after remph up?? in meqrec%stphlnewnode%lines)
! irem is index to meqrec%phr(), meqrec%stphl(jph) is index to meqrec%phr
          meqrec%nstph=meqrec%nstph-1
          do iph=1,meqrec%nstph
             jj=meqrec%stphl(iph)
             if(jj.ge.remph) then
                meqrec%stphl(iph)=meqrec%stphl(iph+1)
             endif
          enddo
! we must zero the last phase, hm itrem is not really relevant ...
          meqrec%stphl(meqrec%nstph+1)=0
          meqrec%phr(remph)%itrem=meqrec%noofits
          meqrec%phr(remph)%prevam=zero
          meqrec%phr(remph)%stable=0
          meqrec%phr(remph)%curd%amfu=zero
!          write(*,*)'SMP lineeq 3: ',meqrec%nstph,meqrec%stphl(1)
       elseif(phfix.gt.0) then
! we have to add phase phfix to the stable phase set, that is no problem
! as it is already in all lists, just remove that it should be fix
!          write(kou,88)phfix,' appears,   ',meqrec%nstph
! meqrec%nfixph seems not to be used ....
          if(meqrec%nfixph.gt.0) then
             meqrec%fixph(1,meqrec%nfixph)=0
             meqrec%fixph(2,meqrec%nfixph)=0
             meqrec%phr(phfix)%phasestatus=PHENTSTAB
             meqrec%nfixph=meqrec%nfixph-1
          endif
       else
          write(*,*)'This is another never never error'
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
       allocate(newnode%linehead(1)%stableph(meqrec%nstph))
       newnode%linehead(1)%nstabph=0
       do iph=1,meqrec%nstph
          newnode%linehead(1)%nstabph=newnode%linehead(1)%nstabph+1
          jj=meqrec%stphl(iph)
!          newnode%linehead(1)%stableph(iph)%phaseix=meqrec%phr(jj)%iph
          newnode%linehead(1)%stableph(iph)%ixphase=meqrec%phr(jj)%iph
          newnode%linehead(1)%stableph(iph)%compset=meqrec%phr(jj)%ics
       enddo
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
! this is probably redundant, fixph already reset
       if(meqrec%nfixph.gt.0) then
          meqrec%fixph(1,meqrec%nfixph)=0
          meqrec%fixph(2,meqrec%nfixph)=0
          meqrec%phr(phfix)%phasestatus=PHENTSTAB
          meqrec%nfixph=meqrec%nfixph-1
       endif
!-------------- 
! no need for loop here I guess ... but I am oldfashioned
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
! with tie-lines in the plane there is always just one stable phase
          allocate(newnode%linehead(jj)%stableph(1))
          allocate(newnode%linehead(jj)%stablepham(1))
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
! In meqrec%phr there are currently two fixed phases, one which was fixed
! along the line (LFIX in mapline%linefixph), the other fixed for the node
! point, given by PHFIX which is an index to meqrec%phr.  The third phase
! was stable with positive amount along the line LENT)
! The three lines are: FIX    STABLE    UNSTABLE
! already done         LFIX   LENT      PHFIX
! exit 1               PHFIX  LFIX      LENT
! exit 2               LENT   PHFIX     LFIX
       jphr=0
       if(allocated(mapline%linefixph)) then
          if(size(mapline%linefixph).gt.1) then
! If there are 3 axis this would be OK
             write(*,*)'Problem, too many fix phases ...'
             gx%bmperr=4240; goto 290
          endif
       endif
       do jj=1,meqrec%nphase
! loop through whole phr array to be sure nothing is wrong
          if(meqrec%phr(jj)%stable.eq.1) then
             if(jj.eq.abs(phfix) .or.&
!                  (meqrec%phr(jj)%iph.eq.mapline%linefixph(1)%phaseix .and.&
                  (meqrec%phr(jj)%iph.eq.mapline%linefixph(1)%ixphase .and.&
                  meqrec%phr(jj)%ics.eq.mapline%linefixph(1)%compset)) cycle
             if(jphr.gt.0) then
                write(*,*)'Problems, two entered phases: ',jj,jphr
                gx%bmperr=4241; goto 290
             else
                jphr=jj
             endif
          endif
       enddo
       zph=0
       do jj=1,meqrec%nphase
!          if(meqrec%phr(jj)%iph.eq.mapline%linefixph(1)%phaseix .and. &
          if(meqrec%phr(jj)%iph.eq.mapline%linefixph(1)%ixphase .and. &
               meqrec%phr(jj)%ics.eq.mapline%linefixph(1)%compset) then
! this is the index in phr for the phase that was fix along the line
             zph=jj
          endif
       enddo
       if(zph.eq.0) then
          write(*,203)' *** warning: cannot find the fix phase: ',zph,&
               mapline%linefixph(1)%ixphase,mapline%linefixph(1)%compset
!               mapline%linefixph(1)%phaseix,mapline%linefixph(1)%compset
203       format(a,10i4)
!       else
!          write(*,203)' Found the fix phase: ',zph,&
!               mapline%linefixph(1)%phase,mapline%linefixph(1)%compset
       endif
! jphr is the phase that was stable along the line
       if(jphr.eq.0) then
          write(*,*)'Problems, not a single entered phase!'
          gx%bmperr=4242; goto 290
       endif
! For isothermal sections with no potential axis we must change the axis
! condition when following a new line
!       if(nopotax.eq.0) then
!          write(*,*)'Changing the axis variable for the new entered phase'
!       endif
! In mapnode there is a nfixph and array linefixph
!       if(ocv()) write(*,207)mapline%linefixph(1)%phaseix,&
       if(ocv()) write(*,207)mapline%linefixph(1)%ixphase,&
            mapline%linefixph(1)%compset,&
            meqrec%phr(phfix)%iph,meqrec%phr(phfix)%ics,&
            meqrec%phr(jphr)%iph,meqrec%phr(jphr)%ics
207    format('LFIX: ',2i3,5x,' PHFIX: ',2i3,5x,' LENT: ',2i3)
! The two exits are:   FIX    STABLE    UNSTABLE
! exit 1               PHFIX  LFIX      LENT
! exit 2               LENT   PHFIX     LFIX
!       iph=mapline%linefixph(1)%phaseix
       iph=mapline%linefixph(1)%ixphase
       ics=mapline%linefixph(1)%compset
!       newnode%linehead(1)%linefixph%phaseix=meqrec%phr(phfix)%iph
       newnode%linehead(1)%linefixph%ixphase=meqrec%phr(phfix)%iph
       newnode%linehead(1)%linefixph%compset=meqrec%phr(phfix)%ics
       newnode%linehead(1)%nstabph=1
! the previously fix phase is set as entered with stablepham as initial amount
!       newnode%linehead(1)%stableph(1)%phaseix=iph
       newnode%linehead(1)%stableph(1)%ixphase=iph
       newnode%linehead(1)%stableph(1)%compset=ics
       newnode%linehead(1)%stablepham(1)=one
! store the phase number that must not become stable in nodfixph
       newnode%linehead(1)%nodfixph=jphr
!-----------
!       newnode%linehead(2)%linefixph%phaseix=meqrec%phr(jphr)%iph
       newnode%linehead(2)%linefixph%ixphase=meqrec%phr(jphr)%iph
       newnode%linehead(2)%linefixph%compset=meqrec%phr(jphr)%ics
       newnode%linehead(2)%nstabph=1
!       newnode%linehead(2)%stableph(1)%phaseix=meqrec%phr(phfix)%iph
       newnode%linehead(2)%stableph(1)%ixphase=meqrec%phr(phfix)%iph
       newnode%linehead(2)%stableph(1)%compset=meqrec%phr(phfix)%ics
       newnode%linehead(2)%stablepham(1)=one
       newnode%linehead(2)%nodfixph=zph
       if(nopotax.eq.0) then
! If we have no potential axis we MUST change the axis condition
! to represent the axis composition of the new stable phase
          write(*,712)stepax,axarr(abs(stepax))%axcond(1)%statevarid,stepaxval
712       format('MAPPING WITH NO POTENTIAL AXIS: ',2i4,1pe12.4)
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
!          svrtarget%phase=newnode%linehead(1)%stableph(1)%phaseix
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
!          write(*,*)'Setting condition for new line to ',xxx
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
!          svrtarget%phase=newnode%linehead(2)%stableph(1)%phaseix
          svrtarget%phase=newnode%linehead(2)%stableph(1)%ixphase
          svrtarget%compset=newnode%linehead(2)%stableph(1)%compset
! This extracts the composition of the entered phase for second new line
          call state_variable_val(svr2,xxx,tmpceq)
          if(gx%bmperr.ne.0) goto 1000
          call condition_value(0,pcond,xxx,tmpceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error setting new conditions in lineceq'
             goto 1000
          endif
!          write(*,211)'New conditions at line: ',newnode%linehead(2)%lineid,&
!               trim(newnode%linehead(2)%lineceq%eqname),&
!               svr2%phase,svr2%compset,xxx
!          call list_conditions(kou,newnode%linehead(2)%lineceq)
       endif
!
       if(ocv()) write(*,56)'Created linehead 1 for node: ',newnode%seqx,&
!            newnode%linehead(1)%linefixph%phaseix,&
            newnode%linehead(1)%linefixph%ixphase,&
            newnode%linehead(1)%linefixph%compset,&
!            newnode%linehead(1)%stableph(1)%phaseix,&
            newnode%linehead(1)%stableph(1)%ixphase,&
            newnode%linehead(1)%stableph(1)%compset,&
            newnode%linehead(1)%nodfixph
       if(ocv()) write(*,56)'Created linehead 2 for node: ',newnode%seqx,&
!            newnode%linehead(2)%linefixph%phaseix,&
            newnode%linehead(2)%linefixph%ixphase,&
            newnode%linehead(2)%linefixph%compset,&
!            newnode%linehead(2)%stableph(1)%phaseix,&
            newnode%linehead(2)%stableph(1)%ixphase,&
            newnode%linehead(2)%stableph(1)%compset,&
            newnode%linehead(2)%nodfixph
56     format(a,i3,5x,2i3,5x,2i3,5x,2i3)
! the fix and stable phases must be copied to meqrec when line is started
!       write(*,*)'Created node with 2 exits: ',newnode%seqx,ceq%tpval(1)
! prevent the lines from being used as that makes the program crash
290    continue
!======================================================================
    case(3) ! Normal node in a phase diagram without tie-lines in plane
! Two crossing lines, one in and 3 exits
! THERE IS NO CASE WHEN FINDING AN INVARIANT
! this is probably redundant, fixph already reset
       if(meqrec%nfixph.gt.0) then
!          write(*,*)'Not redundant ...'
          meqrec%fixph(1,meqrec%nfixph)=0
          meqrec%fixph(2,meqrec%nfixph)=0
          meqrec%phr(abs(phfix))%phasestatus=PHENTSTAB
          meqrec%nfixph=meqrec%nfixph-1
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
          else
             if(phfix.lt.0) then
!                write(*,*)'Allocating stableph: ',jj,jphr-1
                allocate(newnode%linehead(jj)%stableph(jphr-1))
                allocate(newnode%linehead(jj)%stablepham(jphr-1))
             else
!                write(*,*)'Allocating stableph: ',jj,jphr+1
                allocate(newnode%linehead(jj)%stableph(jphr+1))
                allocate(newnode%linehead(jj)%stablepham(jphr+1))
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
             write(*,*)'Problem, too many fix phases ...'
             gx%bmperr=4240; goto 390
          endif
       endif
       stabph=0
       do jj=1,meqrec%nphase
! loop through whole phr array to be sure nothing is wrong
          if(meqrec%phr(jj)%stable.eq.1) then
! there should be 2 fixed phases, one along the line and one at the node
! If 3 or more axis there will be more fixed phases, phfix can be negative
!             if(jj.eq.abs(phfix) .or.&
!                  (meqrec%phr(jj)%iph.eq.mapline%linefixph(1)%phase .and.&
!                   meqrec%phr(jj)%ics.eq.mapline%linefixph(1)%compset)) cycle
! we should include phfix in stabph!!
!             if(meqrec%phr(jj)%iph.eq.mapline%linefixph(1)%phaseix .and.&
             if(meqrec%phr(jj)%iph.eq.mapline%linefixph(1)%ixphase .and.&
                  meqrec%phr(jj)%ics.eq.mapline%linefixph(1)%compset) cycle
             stabph=stabph+1
          endif
       enddo
! Hm, stabph calculated this way is wrong, use mapline%nstabph
!       write(*,312)meqrec%nphase,stabph,mapline%nstabph,phfix,&
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
!       iph=mapline%linefixph(1)%phaseix
       iph=mapline%linefixph(1)%ixphase
       ics=mapline%linefixph(1)%compset
! for use below I need to know the position of iph+ics in meqrec%phr ...
       flfix: do jj=1,meqrec%nstph
          if(meqrec%phr(jj)%iph.eq.iph .and. meqrec%phr(jj)%ics.eq.ics) then
             lfix=jj; exit flfix
          endif
       enddo flfix
       kph=mapline%meqrec%phr(abs(phfix))%iph
       kcs=mapline%meqrec%phr(abs(phfix))%ics
! exit 1 has same linefix as incomming line ------------------------
!       newnode%linehead(1)%linefixph%phaseix=iph
       newnode%linehead(1)%linefixph%ixphase=iph
       newnode%linehead(1)%linefixph%compset=ics
       if(phfix.gt.0) then
!          write(*,*)'allocated size of stableph 2: ',size(mapline%stableph)
          do jj=1,stabph
!             newnode%linehead(1)%stableph(jj)%phaseix=&
!                  mapline%stableph(jj)%phaseix
             newnode%linehead(1)%stableph(jj)%ixphase=&
                  mapline%stableph(jj)%ixphase
             newnode%linehead(1)%stableph(jj)%compset=&
                  mapline%stableph(jj)%compset
             newnode%linehead(1)%stablepham(jj)=mapline%stablepham(jj)
          enddo
! add phfix as stable phase
          jj=stabph+1
!          newnode%linehead(1)%stableph(jj)%phaseix=kph
          newnode%linehead(1)%stableph(jj)%ixphase=kph
          newnode%linehead(1)%stableph(jj)%compset=kcs
          newnode%linehead(1)%stablepham(jj)=zero
!          newnode%linehead(1)%nodfixph=meqrec%phr(abs(phfix))%iph
          newnode%linehead(1)%nodfixph=abs(phfix)
          newnode%linehead(1)%nstabph=jj
       else
          kk=0
          do jj=1,stabph-1
! remove -phfix as stable phase
!             if(mapline%stableph(jj)%phaseix.eq.kph .and.&
             if(mapline%stableph(jj)%ixphase.eq.kph .and.&
                  mapline%stableph(jj)%compset.eq.kcs) then
                kk=jj+1
             else
                kk=kk+1
             endif
!             newnode%linehead(1)%stableph(jj)%phaseix=&
!                  mapline%stableph(kk)%phaseix
             newnode%linehead(1)%stableph(jj)%ixphase=&
                  mapline%stableph(kk)%ixphase
             newnode%linehead(1)%stableph(jj)%compset=&
                  mapline%stableph(kk)%compset
             newnode%linehead(1)%stablepham(jj)=mapline%stablepham(kk)
          enddo
!          newnode%linehead(1)%nodfixph=meqrec%phr(abs(phfix))%iph
          newnode%linehead(1)%nodfixph=abs(phfix)
          newnode%linehead(1)%nstabph=stabph-1
       endif
!
! exit 2 has PHFIX as linefix ----------------------------------
!       newnode%linehead(2)%linefixph%phaseix=kph
       newnode%linehead(2)%linefixph%ixphase=kph
       newnode%linehead(2)%linefixph%compset=kcs
       if(phfix.gt.0) then
          do jj=1,stabph
!             newnode%linehead(2)%stableph(jj)%phaseix=&
!                  mapline%stableph(jj)%phaseix
             newnode%linehead(2)%stableph(jj)%ixphase=&
                  mapline%stableph(jj)%ixphase
             newnode%linehead(2)%stableph(jj)%compset=&
                  mapline%stableph(jj)%compset
             newnode%linehead(2)%stablepham(jj)=mapline%stablepham(jj)
          enddo
! add LFIX as stable phase
          jj=stabph+1
!          newnode%linehead(2)%stableph(jj)%phaseix=iph
          newnode%linehead(2)%stableph(jj)%ixphase=iph
          newnode%linehead(2)%stableph(jj)%compset=ics
          newnode%linehead(2)%stablepham(jj)=zero
!          newnode%linehead(2)%nodfixph=meqrec%phr(lfix)%iph
          newnode%linehead(2)%nodfixph=lfix
          newnode%linehead(2)%nstabph=jj
       else
          kk=0
          do jj=1,stabph-1
! remove -phfix as stable phase
!             if(mapline%stableph(jj)%phaseix.eq.kph .and.&
             if(mapline%stableph(jj)%ixphase.eq.kph .and.&
                  mapline%stableph(jj)%compset.eq.kcs) then
                kk=jj+1
             else
                kk=kk+1
             endif
!             newnode%linehead(2)%stableph(jj)%phaseix=&
!                  mapline%stableph(kk)%phaseix
             newnode%linehead(2)%stableph(jj)%ixphase=&
                  mapline%stableph(kk)%ixphase
             newnode%linehead(2)%stableph(jj)%compset=&
                  mapline%stableph(kk)%compset
             newnode%linehead(2)%stablepham(jj)=mapline%stablepham(kk)
          enddo
!          newnode%linehead(2)%nodfixph=meqrec%phr(lfix)%iph
          newnode%linehead(2)%nodfixph=lfix
          newnode%linehead(2)%nstabph=stabph-1
       endif
!
! exit 3 has PHFIX as linefix ----------------------------------
!       newnode%linehead(3)%linefixph%phaseix=kph
       newnode%linehead(3)%linefixph%ixphase=kph
       newnode%linehead(3)%linefixph%compset=kcs
       do jj=1,stabph
!          if(mapline%stableph(jj)%phaseix.eq.kph .and. &
          if(mapline%stableph(jj)%ixphase.eq.kph .and. &
               mapline%stableph(jj)%compset.eq.kcs) then
! exchange PHFIX for LFIX as stable phase
!             newnode%linehead(3)%stableph(jj)%phaseix=iph
             newnode%linehead(3)%stableph(jj)%ixphase=iph
             newnode%linehead(3)%stableph(jj)%compset=ics
             newnode%linehead(3)%stablepham(jj)=zero
          else
!             newnode%linehead(3)%stableph(jj)%phaseix=&
!                  mapline%stableph(jj)%phaseix
             newnode%linehead(3)%stableph(jj)%ixphase=&
                  mapline%stableph(jj)%ixphase
             newnode%linehead(3)%stableph(jj)%compset=&
                  mapline%stableph(jj)%compset
             newnode%linehead(3)%stablepham(jj)=mapline%stablepham(jj)
          endif
       enddo
!       newnode%linehead(3)%nodfixph=meqrec%phr(lfix)%iph
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
!                  (newnode%linehead(jj)%stableph(kk)%phaseix,&
                  (newnode%linehead(jj)%stableph(kk)%ixphase,&
                  newnode%linehead(jj)%stableph(kk)%compset,&
                  kk=1,newnode%linehead(jj)%nstabph)
          enddo
356       format('Tie-line NOT in plane node exits: ',&
               i2,i3,i4,i2,i5,i3,10(i4,i2))
       endif
! the fix and stable phases must be copied to meqrec when line is started
!       write(*,*)'Created node with 2 exits: ',newnode%seqx,ceq%tpval(1)
390    continue
! invariants for isopleths, number of stable phases equal to components+1
! number of adjacent regions with "components" stable phases is "components+1"
! number of exit lines are 2*(components+1)
! each line has a fix phase and one of the phases stable at the invariant
! as not stable.  The remaining phases  are entered.
! Each phase is fix for two lines and not stable for two others
! This is the way to generate the exit lines:
! - loop for all phases to set a phase fix (for two lines)
! - loop for the next two phases to set one phase not stable
! the remaining phases are set entered (amount?) generate a line startpoint
! take care of remobing line into the invariant

! How to know if the node is invariant? Gibbs phase rule, Degrees of freedom
! f = n + 2 - p
! where n is number of components, 2 if T and P variable, 1 if T or P variable,
! 0 if both T and P fixed, p is number of stable phases.
    end select
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
  end subroutine map_newnode

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
    type(map_fixph), pointer :: mapfix
!\end{verbatim}
    type(map_node), pointer :: mapnode
    type(gtp_condition), pointer :: pcond
    integer nyline,jp,seqy,iph,ics,lokph,lokcs,ip
    double precision finc
    character eqname*24
! sometimes there are many phases with long names ...
    character phaseset*512
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
       allocate(mapfix)
! with only 2 axis we have just 1 fix phase for mapping
       allocate(mapfix%fixph(1))
       mapfix%nfixph=1
       mapfix%fixph=mapline%linefixph(1)
! we can have several stable phases when no tie-lines in plane
       ip=mapline%nstabph
       allocate(mapfix%stableph(ip))
       allocate(mapfix%stablepham(ip))
!       write(*,*)'Findline: Tie-lines not in plane: ',nyline,ip
! create a heading text for the line
       phaseset=' '
       call get_phasetup_name_old(mapfix%fixph(1),phaseset)
       if(gx%bmperr.ne.0) goto 1000
       ip=len_trim(phaseset)+4
       phaseset(ip-1:ip-2)='+'
       if(mapnode%linehead(nyline)%nstabph.le.0) then
          write(*,*)'No stable phases for a line'
          write(*,*)'Error 14:  ',nyline,mapnode%linehead(nyline)%nstabph,&
!               mapnode%linehead(nyline)%stableph(1)%phaseix,&
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
          call get_phasetup_name_old(mapfix%stableph(jp),phaseset(ip:))
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
520    format(/'Line ',i3,' T=',F8.2,' with: ',a)
!       nullify(mapfix)
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
       allocate(mapfix)
       allocate(mapfix%fixph(1))
       allocate(mapfix%stableph(1))
       allocate(mapfix%stablepham(1))
       mapfix%nfixph=1
       mapfix%fixph=mapline%linefixph(1)
! create a heading text for the line
       phaseset=' '
       call get_phasetup_name_old(mapfix%fixph(1),phaseset)
       if(gx%bmperr.ne.0) goto 1000
       ip=len_trim(phaseset)+4
       phaseset(ip-2:ip-2)='+'
!       write(*,*)'Fixed phase: ',mapfix%nfixph,&
!            mapfix%fixph%ixphase,mapfix%fixph%compset
       if(mapnode%linehead(nyline)%nstabph.gt.0) then
! this is stored only for "real" nodes
          mapfix%nstabph=1
          mapfix%stableph(1)=mapnode%linehead(nyline)%stableph(1)
          call get_phasetup_name_old(mapfix%stableph(1),phaseset(ip:))
          if(gx%bmperr.ne.0) goto 1000
! set positive amount both in mapfix and in phase_varres ...
!          mapfix%stablepham(1)=mapnode%linehead(nyline)%stablepham(1)
          mapfix%stablepham(1)=one
!          call get_phase_compset(mapfix%stableph(1)%phase,&
!               mapfix%stableph(1)%compset,lokph,lokcs)
!          mapline%lineceq%phase_varres(lokcs)%amfu=one
          ip=len_trim(phaseset)
          if(ip.gt.1) then
             write(kou,516)mapline%lineid,&
                  mapline%lineceq%tpval(1),phaseset(1:ip)
516          format(/'Line: ',i3,' T=',F8.2,' with: ',a)
!             write(*,507)' *** Phase numbers: ',mapfix%fixph(1),&
!                  mapfix%stableph(1),mapline%nodfixph
507          format(a,2i3,2x,2i3,2x,2i3)
          else
             write(kou,521)
521          format(/'Line with unknown phases, wow')
          endif
       else
          write(*,*)'No stable phase!! why??'
          write(*,*)'stable 4:  ',nyline,mapnode%linehead(nyline)%nstabph,&
!               mapnode%linehead(nyline)%stableph(1)%phaseix,&
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
          call get_phasetup_name_old(mapnode%linehead(1)%stableph(jp),&
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
       nullify(mapfix)
    endif
1000 continue
    return
  end subroutine map_findline

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine create_saveceq(ceqres,size)
! creates an array of equilibrium records to save calculated lines for step
! and map
    type(map_ceqresults), pointer :: ceqres
    integer size
!\end{verbatim}
    if(ocv()) write(*,*)'In create saveceq'
    allocate(ceqres)
    ceqres%size=size
    ceqres%free=1
    allocate(ceqres%savedceq(size))
1000 continue
    return
  end subroutine create_saveceq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine delete_mapresults(maptop)
! delete all saved results created by step or map
    TYPE(map_node), pointer :: maptop
!\end{verbatim}
    type(map_ceqresults), pointer :: saveceq
    TYPE(map_node), pointer :: current
    TYPE(gtp_equilibrium_data), pointer :: ceq
    integer ieq
!    integer place,lastused
!
    if(.not.associated(maptop)) then
       write(*,*)'No step or map results to delete'
       goto 1000
    endif
    current=>maptop
!    deloop: do while(associated(current))
!       write(*,*)'Number of stored equilibria: ',current%saveceq%free-1
!       current=>current%plotlink
!    enddo deloop
!    write(*,*)'All mapnodes listed'
! all mapnodes has a pointer to first where the saveceq is allocated
    current=>maptop
    do while(associated(current))
!       write(*,*)'deleting map results: ',current%saveceq%free-1
       deallocate(current%saveceq%savedceq)
       current=>current%plotlink
    enddo
!    write(*,*)'SMP: deleting _MAPx equilibria'
    ceq=>firsteq
    call delete_equilibria('_MAP*',ceq)
!    write(*,*)'Done delete_mapresults'
1000 continue
    return
  end subroutine delete_mapresults

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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

!\begin{verbatim}
  logical function inveq(phases,ceq)
! Only called for tie-lines not in plane.  If tie-lines in plane then all
! nodes are invariants.
! UNFINISHED
    integer phases
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer nrel,ii,nostph,tpvar,degf
    type(gtp_condition), pointer :: pcond,lastcond
    type(gtp_state_variable), pointer :: stvr
! How to know if the ceq is invariant? Gibbs phase rule, Degrees of freedom
! f = n + 2 - p
! where n is number of components, 2 if T and P variable, 1 if T or P variable,
! 0 if both T and P fixed, p is number of stable phases.
!    write(*,*)'in inveq'
    nrel=noel()
! sum up nubler of stable phases and check if T and P are fixed
    nostph=0
    do ii=1,noofphasetuples()
       if(ceq%phase_varres(phasetuple(ii)%lokvares)%phstate.gt.0) &
            nostph=nostph+1
    enddo
! loop all conditions
    lastcond=>ceq%lastcondition
    pcond=>lastcond
    tpvar=2
100 continue
       if(pcond%active.eq.0) then
! condtion is active
          stvr=>pcond%statvar(1)
! statevarid 1 is T and 2 is P
          if(stvr%statevarid.eq.1 .or. stvr%statevarid.eq.2) then
! Hm, ceq is not the equilibrium record for the node point ...
             tpvar=tpvar-1
          endif
       endif
       pcond=>pcond%next
       if(.not.associated(pcond,lastcond)) goto 100
!
    degf=nrel+tpvar-nostph
!    write(*,*)'in inveq 2',degf,nrel,tpvar,nostph
    if(degf.eq.0) then
! number of exit phases is equal to the number of stable phases ???
       phases=nostph
       inveq=.true.
!       write(*,200)'We have an invariant equilibrium!',nrel,tpvar,nostph,phases
200    format(a,5i7)
    else
!       write(*,210)degef,nrel,tpvar,phases
210     format('Not inveq, elements, stable phases: ',4i4)
!       if not invariant there are 3 exits (2 lines crossing)
        phases=degf
       inveq=.false.
    endif
1000 continue
    return
  end function inveq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
    if(mapline%problems.gt.5) then
       if(mapline%nodfixph.gt.0) then
!          call list_conditions(kou,mapline%lineceq)
!          if(gx%bmperr.ne.0) then
!             write(*,*)'Error listing conditions'
!             gx%bmperr=0
!          endif
          write(*,11)mapline%lineid,trim(mapline%lineceq%eqname),&
               mapline%meqrec%phr(mapline%nodfixph)%iph,&
               mapline%meqrec%phr(mapline%nodfixph)%ics
11        format('I give up on this line',i3,2x,a,' with fix phase ',2i4)
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

!\begin{verbatim}
  subroutine map_halfstep(halfstep,axvalok,mapline,axarr,ceq)
! Used when an error calculating a normal step or a node point
! take back the last sucessfully calculated axis value and take smaller step
! possibly one should also restore the ceq record.
! halfstep number of times halfstep has been called for this problem
! axvalok last cucessfully calculated value of active axis
! mapline line record
! axarr array with axis records
! ceq equilibrium record
    implicit none
    integer halfstep
    double precision axvalok
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(map_line), pointer :: mapline
    type(map_axis), dimension(*) :: axarr
!\end{verbatim}
    type(gtp_condition), pointer :: pcond
    double precision value
    double precision, parameter :: sfact=1.0D-2
    integer jax
!    write(*,*)'In map_halfstep'
    halfstep=halfstep+1
    if(axvalok.eq.zero .or. halfstep.ge.3) then
       if(ocv()) write(*,*)'Problem resolving two phases competing to',&
            ' appear/disappear'
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
!       write(*,97)'Setting axis value: ',mapline%axandir,value,axvalok,&
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

!\begin{verbatim}
  subroutine list_stored_equilibria(kou,axarr,maptop)
! list all nodes and lines from step/map
! use amed to exclude/include lines
! kou output unit
! axarr array with axis records
! maptop map node record
    integer kou
    type(map_node), pointer :: maptop
    type(map_axis), dimension(*) :: axarr
!\end{verbatim} %+
    type(map_ceqresults), pointer :: results
    type(map_node), pointer :: mapnode,localtop
    type(gtp_equilibrium_data), pointer :: thisceq
    type(gtp_condition), pointer :: pcond
    integer kl,ll,jax,nax
    double precision, dimension(:), allocatable :: axxx
    type(gtp_state_variable), pointer :: svrrec
!
    if(.not.associated(maptop)) then
       write(kou,*)'No stored equilibria'
       goto 1000
    endif
    if(associated(maptop%plotlink)) then
       write(*,*)'The plotlink is set !!! '
    endif
    results=>maptop%saveceq
    if(.not.associated(results)) then
       write(kou,*)'No stored equilibria'
       goto 900
    endif
    nax=maptop%number_ofaxis
    allocate(axxx(nax))
    write(kou,90)
90  format('List of all stored equilibria')
! if there has been several STEP/MAP there can be several localtop
    localtop=>maptop
!
! return here if there has been several step/map commands
99  continue
    mapnode=>localtop
! list all mapnodes for this step/map command
100 continue
!    mapnode=>localtop
    write(kou,101)mapnode%seqx,mapnode%noofstph,mapnode%savednodeceq
101 format(' Mapnode: ',i5,' with ',i2,' stable phases, ceq saved in ',i5)
    write(*,*)'Number of exit lines: ',mapnode%lines
    do kl=1,mapnode%lines
       if(.not.associated(mapnode%linehead(kl)%end)) then
          if(mapnode%linehead(kl)%termerr.gt.0) then
             write(kou,105)mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria,&
                  mapnode%linehead(kl)%termerr
105          format('  Line ',i3,' with ',i3,&
                  ' equilibria ended with error: ',i6)
          else
             write(kou,110)mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria
110          format('  Line ',i3,' with ',i3,&
                  ' equilibria ending at axis limit')
          endif
       else
          ll=mapnode%linehead(kl)%end%seqx
          write(kou,120)mapnode%linehead(kl)%lineid,&
               mapnode%linehead(kl)%number_of_equilibria,ll
120       format('  Line ',i3,' with ',i3,' equilibria ending at node ',i3)
       endif
       if(btest(mapnode%linehead(kl)%status,EXCLUDEDLINE)) then
          write(*,*)'Line excluded'
          cycle
       endif
       ll=mapnode%linehead(kl)%first
!       write(*,*)'list first equilibrium ',ll
!       write(*,*)'axis: ',mapnode%number_ofaxis
       if(.not.allocated(results%savedceq)) then
          write(*,*)'Cannot find link to saved equilibria! '
       else
          ceqloop: do while(ll.gt.0)
             thisceq=>results%savedceq(ll)
             do jax=1,nax
                call locate_condition(axarr(jax)%seqz,pcond,thisceq)
                if(gx%bmperr.ne.0) goto 1000
                svrrec=>pcond%statvar(1)
                call state_variable_val(svrrec,axxx(jax),thisceq)
                if(gx%bmperr.ne.0) goto 1000
             enddo
             write(kou,150)ll,thisceq%next,thisceq%tpval(1),axxx
150          format('     Saved ceq ',2i5,f8.2,5(1pe14.6))
             ll=thisceq%next
          enddo ceqloop
       endif
    enddo
!    write(kou,160)mapnode%seqx,mapnode%previous%seqx
160 format('Current node: ',i2,' followed by: ',i2)
    mapnode=>mapnode%previous
!    localtop=>localtop%previous
    if(.not.associated(mapnode,localtop)) goto 100
! plotlink needed if there has been several step/map commands
900 continue
    if(associated(localtop%plotlink)) then
       write(lut,910)
910    format(/'Results from a previous step/map command,',&
            ' equilibrium numbers will overlap')
       localtop=>localtop%plotlink
       write(*,*)'Setting result link'
       results=>localtop%saveceq
       if(.not.associated(results)) then
          write(kou,*)'No stored equilibria'
          goto 900
       endif
       goto 99
    endif
!
    write(kou,*)'That is all'
1000 continue
    return
  end subroutine list_stored_equilibria

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
  subroutine amend_stored_equilibria(axarr,maptop)
! allows amending inactive/acive status of all lines from step/map
    type(map_node), pointer :: maptop
    type(map_axis), dimension(*) :: axarr
!\end{verbatim}
    character cline*12,status*8,ch1*1,phline*78
    type(map_ceqresults), pointer :: results
    type(map_node), pointer :: mapnode,localtop,testnode
    type(gtp_equilibrium_data), pointer :: thisceq
    type(gtp_condition), pointer :: pcond
    integer kl,ll,jax,last,nax
    double precision, dimension(:), allocatable :: axxx
    type(gtp_state_variable), pointer :: svrrec
!
    if(.not.associated(maptop)) then
       write(kou,*)'No stored equilibria'
       goto 1000
    endif
    if(associated(maptop%plotlink)) then
       write(*,*)'There is more than one maptop record'
    endif
    localtop=>maptop
!    if(associated(localtop,maptop)) write(*,*)'maptop same as localtop'       
    nax=localtop%number_ofaxis
    allocate(axxx(nax))
    write(kou,90)
    nullify(results)
90  format('Amend all stored equilibria ... from several maptops')
!
! return here if %plotlink is not empty
! each plotlink has its own results link to saveceq
99  continue
! there can be lines associated with several maptops
! but I have trouble finding them.  They should be linked by plotlink
    mapnode=>localtop
    results=>mapnode%saveceq
100 continue
!    mapnode=>maptop
    if(associated(localtop,mapnode)) write(*,*)'mapnode same as localtop'       
    if(.not.associated(results)) then
       write(kou,*)'No stored equilibria'
       goto 900
    endif
    status=' '
    write(kou,101)mapnode%seqx,mapnode%noofstph,mapnode%savednodeceq
101 format(' Mapnode: ',i5,' with ',i2,' stable phases, ceq saved in ',i5)
    write(*,*)'Number of exit lines: ',mapnode%lines
    lineloop: do kl=1,mapnode%lines
       if(mapnode%linehead(kl)%number_of_equilibria.eq.0) then
          write(*,*)'Skipping empty line >>>'
          cycle lineloop
       endif
       if(btest(mapnode%linehead(kl)%status,EXCLUDEDLINE)) then
          status='EXCLUDED'
       else
          status='INCLUDED'
       endif
       call line_with_phases_withdgm0(phline,mapnode%linehead(kl)%lineceq)
       if(gx%bmperr.ne.0) then
          phline='Sorry cannot list stable phases'
          gx%bmperr=0
       endif
       if(.not.associated(mapnode%linehead(kl)%end)) then
          if(mapnode%linehead(kl)%termerr.gt.0) then
             write(kou,105)mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria,&
                  mapnode%linehead(kl)%termerr,status,trim(phline)
105          format('  Line ',i3,' with ',i3,&
                  ' equilibria ended with error: ',i6,'.  ',a,' with phases:'/a)
          else
             write(kou,110)mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria,status,trim(phline)
110          format('  Line ',i3,' with ',i3,&
                  ' equilibria ending at axis limit.  ',a,' with phases: ',/a)
          endif
       else
          ll=mapnode%linehead(kl)%end%seqx
          write(kou,120)mapnode%linehead(kl)%lineid,&
               mapnode%linehead(kl)%number_of_equilibria,ll,status,trim(phline)
120       format('  Line ',i3,' with ',i3,' equilibria ending at node ',i3,&
               '.  ',a,' with phases: ',/a)
       endif
!       write(*,*)'mapnode%linehead%first: '
       ll=mapnode%linehead(kl)%first
!       write(*,*)'axis conditions: ',ll,axarr(1)%seqz,axarr(2)%seqz
       if(.not.btest(mapnode%linehead(kl)%status,EXCLUDEDLINE)) then
          if(.not.allocated(results%savedceq)) then
             write(*,*)'Cannot find link to saved equilibria! '
          else
!             write(*,*)'Loop for all equlibria'
             do while(ll.gt.0)
                thisceq=>results%savedceq(ll)
!                write(*,*)'thisceq pointer set',ll
                do jax=1,nax
                   call locate_condition(axarr(jax)%seqz,pcond,thisceq)
                   if(gx%bmperr.ne.0) goto 1000
!                   write(*,*)'local_condition NOT OK'
                   svrrec=>pcond%statvar(1)
                   call state_variable_val(svrrec,axxx(jax),thisceq)
                   if(gx%bmperr.ne.0) goto 1000
!                   write(*,*)'state variable found OK'
                enddo
                write(kou,150)ll,thisceq%next,thisceq%tpval(1),axxx
150             format('     Saved ceq ',2i5,f8.2,5(1pe14.6))
                ll=thisceq%next
             enddo
          endif
       endif
! if deleted ask for Restore, else ask for Keep or Delete
       last=len(cline)
       if(btest(mapnode%linehead(kl)%status,EXCLUDEDLINE)) then
          call gparcd('Include this line? ',cline,last,1,ch1,'N',q1help)
       else
          call gparcd('Exclude this line? ',cline,last,1,ch1,'N',q1help)
       endif
       if(biglet(ch1).eq.'Y') then
          if(btest(mapnode%linehead(kl)%status,EXCLUDEDLINE)) then
             mapnode%linehead(kl)%status=&
                  ibclr(mapnode%linehead(kl)%status,EXCLUDEDLINE)
             write(kou,*)'Line activated'
          else
             mapnode%linehead(kl)%status=&
                  ibset(mapnode%linehead(kl)%status,EXCLUDEDLINE)
             write(kou,*)'Line inactivated'
          endif
       elseif(biglet(ch1).eq.'Q') then
          goto 1000
       endif
    enddo lineloop
!    write(kou,160)mapnode%seqx,mapnode%previous%seqx
!160 format('Current node: ',i2,' followed by: ',i2)
    mapnode=>mapnode%previous
!    if(.not.associated(mapnode,maptop)) goto 100
    if(.not.associated(mapnode,localtop)) then
!       if(associated(mapnode,maptop)) write(*,*)'mapnode same as maptop'
       goto 100
    endif
900 continue
! plotlink needed if there has been several step/map commands
    if(associated(localtop%plotlink)) then
       write(lut,910)
910    format(/'Results from a previous step/map command,',&
            ' equilibrium numbers will overlap')
       localtop=>localtop%plotlink
       goto 99
    endif
    write(kou,*)'That is all'
1000 continue
    return
  end subroutine amend_stored_equilibria

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine ocplot2(ndx,pltax,filename,maptop,axarr,graphopt,pform,ceq)
! Main plotting routine, generates a GNUPLOT data file for a step/map calc
! NOTE for isothermal section ocplot3 is used (when 2 axis with wildcards)
! ndx is mumber of plot axis, 
! pltax is text with plotaxis variables
! filename is the name of the GNUPLOT file
! maptop is map_node record with all results
! axarr is array of axis records
! graphopt is graphical option record
! pform is type of output (screen/acrobat/postscript/gif)
! ceq is equilibrium record
    implicit none
    integer ndx
    character pltax(*)*(*),filename*(*),pform*(*)
    type(map_axis), dimension(*) :: axarr
    type(map_node), pointer :: maptop
    type(graphics_options) :: graphopt
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
    type(map_ceqresults), pointer :: results
    TYPE(gtp_equilibrium_data), pointer :: curceq
    type(map_node), pointer :: mapnode,invar,localtop
    type(map_line), pointer :: mapline
    logical wildcard
    character ch1*1,gnuplotline*256,pfd*64,pfc*64
    character pfh*64,dummy*24,applines(10)*128,appline*128
    double precision, dimension(:,:), allocatable :: anp
    double precision, dimension(:), allocatable :: xax,yyy
    integer, parameter :: maxval=10000
    integer, dimension(:), allocatable :: nonzero,linzero,linesep
!    integer, dimension(:), allocatable :: linesep
! encoded2 stores returned text from get_many ... 2048 is too short ...
! selphase used when plotting data just for a selected phase like y(fcc,*)
    character statevar*64,encoded1*1024,encoded2*4096,selphase*24
    character*64, dimension(:), allocatable :: phaseline
    integer i,ic,jj,k3,kk,kkk,lokcs,nnp,np,nrv,nv,nzp,ip,nstep,nnv
    integer nr,line,next,seqx,nlinesep,ksep,iax,anpax,notanp,appfil
    double precision xmax,xmin,ymax,ymin,value,anpmin,anpmax
! lhpos is last used position in lineheader
    integer giveup,nax,ikol,maxanp,lcolor,lhpos,repeat,anpdim,qp
    character date*8,mdate*12,title*128,backslash*2,lineheader*1024
    character deftitle*128,labelkey*64
    logical overflow,first,last,novalues,selectph,varofun,moretops
! textlabels
    type(graphics_textlabel), pointer :: textlabel
! line identification (title)
    character*16, dimension(:), allocatable :: lid
!    logical nolid
!
!    write(*,*)'In ocplot2, looking for segmentation fault 1'
    if(index(pltax(1),'*').gt.0 .and. index(pltax(2),'*').gt.0) then
       call ocplot3(ndx,pltax,filename,maptop,axarr,graphopt,pform,ceq)
       goto 1000
    endif
    moretops=.FALSE.
    seqx=0
    call date_and_time(date)
    mdate=" "//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//" "
    deftitle='Open Calphad 4.0 prerelease '//mdate//': with GNUPLOT'
    if(graphopt%labeldefaults(1).eq.0) then
       title=deftitle
    else
! alwas inlcude open calphad and date, add user title at the end
!       123456789.123456789.123456789
!      'Open Calphad 3.0 2015-03-16 : with GNUPLOT'
       jj=len_trim(deftitle)
       title=deftitle(1:jj+1)//graphopt%plotlabels(1)
    endif
!
    if(.not.associated(maptop)) then
       write(kou,*)'In ocplot2 but nothing to plot'
       gx%bmperr=4247; goto 1000
    endif
!    write(*,*)'Entering OC plot version 2: ',&
!         pltax(1)(1:len_trim(pltax(1))),', ',pltax(2)(1:len_trim(pltax(2))),&
!         maptop%next%seqx
    if(graphopt%rangedefaults(1).ne.0) then
       write(*,11)'x',graphopt%plotmin(1),graphopt%plotmax(1)
11     format('Limits set by user for ',a,': ',2(1pe14.6))
    endif
    if(graphopt%rangedefaults(2).ne.0) then
       write(*,11)'y',graphopt%plotmin(2),graphopt%plotmax(2)
    endif
! allocate as many items in linesep as there are mapnodes.
! Hm, if merging plots the number of separators needed can be any value
    jj=100+10*maptop%next%seqx+1
    allocate(linesep(jj))
    nax=maptop%number_ofaxis
    linesep=0
! allocate texts to identify the lines on the gnuplot file
    allocate(phaseline(jj))
!    if(maptop%number_ofaxis.gt.1) then
!       write(*,*)'Warning: may not not handle map graphics correctly',jj
!    endif
!
    giveup=0
    nrv=maxval
    allocate(xax(nrv))
! to insert MOVE at axis terminations
    nlinesep=1
    phaseline(1)=' '
    phaseline(2)=' '
    nv=0
! min and max not used by gnuplot but may be useful if plotpackage change
! or for manual scaling.
    xmin=1.0D20
    xmax=-1.0D20
    ymin=1.0D20
    ymax=-1.0D20
    maxanp=1000
    np=maxanp
    qp=1
    wildcard=.FALSE.
    selectph=.FALSE.
    do iax=1,2
!       write(*,*)'Allocating for axis: ',iax
       call capson(pltax(iax))
       if(index(pltax(iax),'*').gt.0) then
          if(wildcard) then
             write(*,*)'in OCPLOT2 one axis variable with wildcard allowed'
             goto 1000
          endif
! wildcards allowed only on one axis, we do not know how many columns needed
! allocate as many array elements as columns
          anpdim=np
          allocate(anp(np,nrv))
          allocate(nonzero(np))
          allocate(linzero(np))
          allocate(yyy(np))
! nzp should be dimension of yyy, np returns the number of values in yyy
          nzp=np
          nonzero=0
          wildcard=.TRUE.
          anpax=iax
! when we plot things like y(fcc,*) we should only select equilibria
! with fcc stable. Such state variables contain a ","
          ikol=index(pltax(iax),',')
          if(ikol.gt.0) then
! if the * is after the , then extract the phase name before             
             if(pltax(iax)(ikol+1:ikol+1).eq.'*') then
                nrv=index(pltax(iax),'(')
                if(nrv.lt.ikol) then
                   selphase=pltax(iax)(nrv+1:ikol-1)
!                   write(*,*)'Selected phase: ',trim(selphase)
                   selectph=.TRUE.
                endif
             endif
          endif
       endif
    enddo
    if(.not.wildcard) then
       anpdim=1
       allocate(anp(1,nrv))
       wildcard=.FALSE.
       nnp=1
       anpax=2
    endif
! zero anp, maybe waste of time but ...
    anp=zero
! if anpax is 1, notanp is 2, if anpax is 2, notanp is 1
    notanp=3-anpax
    localtop=>maptop
!    write(*,*)'In ocplot2, looking for segmentation fault 2'
!-------------
! come back here if there is another localtop in plotlink!
77  continue
! change all "done" marks in mapnodes to zero
    ikol=0
    do nrv=1,localtop%lines
       localtop%linehead(nrv)%done=0
    enddo
    mapnode=>localtop%next
    do while(.not.associated(mapnode,localtop))
       do nrv=1,mapnode%lines
          mapnode%linehead(nrv)%done=0
       enddo
       mapnode=>mapnode%next
    enddo
!-----------
    results=>localtop%saveceq
    mapnode=>localtop
    line=1
! extract the names of stable phases for this lone

!    write(*,*)'mapnode index: ',mapnode%seqx
!    write(*,*)'Before label 100: ',results%free
!------------------------------------------- begin loop 100
100    continue
       mapline=>localtop%linehead(line)
! initiate novalues to false for each line
       novalues=.false.
!       write(*,*)'In ocplot2, looking for segmentation fault 3'
! skip line if EXCLUDEDLINE set
       if(btest(mapline%status,EXCLUDEDLINE)) then
!          write(*,*)'Skipping a line 3'
          if(line.lt.mapnode%lines) then
             line=line+1
             goto 100
          else
             goto 500
          endif
       endif
110    continue
! mark line is plotted
!       write(*,*)'Values from mapline ',mapline%lineid
! loop for all calculated equilibra, phases and composition sets
!       write(*,*)'Before label 150: ',mapline%lineid
!150    continue
       nr=mapline%first
       if(mapline%done.ne.0) goto 220
       mapline%done=-1
       if(ocv()) write(*,*)'Plotting line: ',&
            mapline%lineid,mapline%number_of_equilibria
       first=.TRUE.
       last=.FALSE.
! we may have empty lines due to bugs ...
!       write(*,*)'Axis with wildcard and not: ',anpax,notanp
!200    continue
! this is the loop for all equilibria in the line
       plot1: do while(nr.gt.0)
          nv=nv+1
          if(nv.ge.maxval) then
             write(*,*)'Too many points to plot',maxval
             goto 600
          endif
          curceq=>results%savedceq(nr)
          if(ocv()) write(*,201)'Current equilibrium: ',nr,nv,curceq%tpval(1)
201       format(a,2i5,F8.2,1pe14.6)
! extract the names of stable phases to phaseline from the first equilibrium
! Note that the information of the fix phases are not saved
          if(first) then
!             write(*,*)'In ocplot2, looking for segmentation fault 4A'
             first=.false.
             kk=1
             do jj=1,noph()
                do ic=1,noofcs(jj)
                   k3=test_phase_status(jj,ic,value,curceq)
                   if(gx%bmperr.ne.0) goto 1000
                   if(k3.gt.0) then
! this phase is stable or fix
                      call get_phase_name(jj,ic,dummy)
                      if(gx%bmperr.ne.0) goto 1000
                      if(selectph) then
! this is an attempt to remove lines from irrelevant equilibria when plotting
! data for a specific phase like y(fcc#4,*)
                         if(abbr_phname_same(dummy,trim(selphase))) then
                            novalues=.FALSE.
!                            write(*,*)'Setting novalues FALSE',mapline%lineid
                         else
                            novalues=.TRUE.
!                            write(*,*)'Setting novalues TRUE',mapline%lineid
                         endif
                      endif
! I think this phaseline is no longer used ??
                      phaseline(nlinesep)(kk:)=dummy
                      kk=len_trim(phaseline(nlinesep))+2
                   endif
                enddo
             enddo
!             write(*,117)nlinesep,phaseline(nlinesep)(1:kk-1)
117          format('Stable phases on line ',i5/a)
! the segmentation fault was that linzero not always allocated ....
             if(allocated(linzero)) linzero=0
          endif
! no wildcards allowed on this axis
          statevar=pltax(notanp)
          call meq_get_state_varorfun_value(statevar,value,encoded1,curceq)
!          write(*,*)'SMP axis variable 1: ',trim(encoded1),value
          if(gx%bmperr.ne.0) then
! this error should not prevent plotting the other points FIRST SKIPPING
             write(*,212)'SMP skipping a point, error evaluating: ',&
                  statevar(1:10),curceq%tpval(1),nv,nr
212          format(a,a,f10.2,2i5)
             gx%bmperr=0
             nv=nv-1; goto 215
          endif
          xax(nv)=value
!          write(*,201)'at 202: ',nr,nv,curceq%tpval(1),value
!          xax(nv)=curceq%tpval(1)
!          write(*,*)'After label 200: ',mapline%lineid,nr,nv
          if(xax(nv).lt.xmin) xmin=xax(nv)
          if(xax(nv).gt.xmax) xmax=xax(nv)
! second axis
          statevar=pltax(anpax)
!          varofun=.FALSE.
          if(wildcard) then
! NEW ignore data for this equilibrium if NOVALUES is TRUE
! because selphase not equal to the stable phase found above
             if(novalues) then
!                write(*,*)'Ignoring equilibria without ',trim(selphase)
                yyy=zero
                np=1
             else
!                write(*,*)'Getting a wildcard value 1: ',nr,trim(statevar)
                call get_many_svar(statevar,yyy,nzp,np,encoded2,curceq)
! compiling without -finit-local-zero gives a segmentation fault here
! running the MAP11 macro
                qp=np
!             write(*,*)'wildcard value 1: ',nr,trim(statevar)
!                write(*,223)'Values: ',np,(yyy(i),i=1,np)
223             format(a,i3,8F8.4)
                if(gx%bmperr.ne.0) then
                   write(*,*)'yaxis error: "',trim(statevar),'"'
                   goto 1000
                endif
             endif
!             write(*,213)trim(encoded2),np,(yyy(ic),ic=1,np)
213          format('windcard: ',a,i3,6(1pe12.4))
!             write(*,16)'val: ',kp,nr,gx%bmperr,(yyy(i),i=1,np)
16           format(a,2i3,i5,/6(1pe11.3/))
             anpmin=1.0D20
             anpmax=-1.0D20
             lcolor=0
             do jj=1,np
                if(last) then
                   if(linzero(jj).ne.0) then
! in last equilibria we may have a value from the new phase at the node
                      anp(jj,nv)=yyy(jj)
!                   elseif(yyy(jj).ne.zero) then
!                      write(*,*)'SMP skipping a value for line ',nlinesep
                   endif
                else
                   anp(jj,nv)=yyy(jj)
                   if(abs(yyy(jj)).gt.zero) then
                      nonzero(jj)=1
                      linzero(jj)=1
! save the first column with nonzero for use with invariants
                      if(ikol.eq.0) ikol=jj
                      if(anp(jj,nv).gt.anpmax) anpmax=anp(jj,nv)
                      if(anp(jj,nv).lt.anpmin) anpmin=anp(jj,nv)
! extract state variable jj
                      if(.not.allocated(lid)) then
                         allocate(lid(np+5))
!                         write(*,*)'Allocate lid: ',np
                      endif
! getext( , ,2, , , ) returns next text item up to a space
                      call getext(encoded2,lcolor,2,encoded1,'x',lhpos)
                      lid(jj)=encoded1
                      kk=len_trim(encoded1)
                      if(kk.gt.len(lid(jj))) then
                         lid(jj)(7:)='..'//encoded1(kk-6:kk)
                      else
                      endif
                   else
! skip state variable
                      call getext(encoded2,lcolor,2,encoded1,'x ',lhpos)
                   endif
                endif
             enddo
!             write(*,*)'OK Point: ',nr,nv,xax(nv)
          else
! More than one state variable or function value like CP
! UNFINISHED PROBLEM WITH NEGATIVE CP HERE 
! try skipping this value (below) if last equilibrium on the line 
!             varofun=.TRUE.
             call meq_get_state_varorfun_value(statevar,value,encoded1,curceq)
!             write(*,*)'SMP axis variable 2: ',statevar(1:3),value
             if(gx%bmperr.ne.0) then
! SECOND Skipping
                write(*,212)'SMP Skipping a point, error evaluating: ',&
                     statevar(1:10),curceq%tpval(1),nv,nr
                nv=nv-1; goto 215
             endif
             if(results%savedceq(nr)%next.eq.0) then
! THIRD ?? Skipping
!                write(*,212)'SMP skip last evaluated symbol: ',&
!                     trim(statevar),curceq%tpval(1),nv,nr
                nv=nv-1; goto 215
             endif
!             if(gx%bmperr.ne.0) goto 1000
             anp(1,nv)=value
!             write(*,201)'at 19: ',nr,nv,curceq%tpval(1),value
!             write(*,19)'Bug: ',nr,nv,seqx,xax(nv),anp(1,nv)
19           format(a,3i4,2(1pe12.4))
             anpmin=anp(1,nv)
             anpmax=anp(1,nv)
          endif
          if(anpmin.lt.ymin) ymin=anpmin
          if(anpmax.gt.ymax) ymax=anpmax
215       continue
! reset any previous error code
          if(gx%bmperr.ne.0) then
!             write(*,*)'SMP reset error code ',gx%bmperr
             gx%bmperr=0
          endif
          nr=curceq%next
          if(nr.gt.0) then
             if(results%savedceq(nr)%next.eq.0) then
!                write(*,*)'We have found last equilibria along the line: ',last
                last=.TRUE.
             endif
          endif
!>>>>>>>>>>>>>>>>>>>>>>>>>> starting a line
!          write(*,*)'Next equilibrium: ',nr,nv,xax(nv)
!          read(*,17)ch1
17           format(a)
       enddo plot1
220    continue
! finished one line
       if(nax.gt.1) then
!------------------ special for invariant lines
! for phase diagram always move to the new line 
          map1: if(nlinesep.ge.1) then
             if(linesep(nlinesep).lt.nv) then
! we should never have several linesep for the same value of nv!
                nlinesep=nlinesep+1
                linesep(nlinesep)=nv
! phaseline(nlinesep) is already filled with spaces
                phaseline(nlinesep+1)=' '
!                write(*,*)'adding empty line 1',nlinesep,linesep(nlinesep)
                inv: if(localtop%tieline_inplane.gt.0 .and. &
                     associated(mapline%end)) then
!                   write(*,*)'Tie-lines in plane, an invariant equil here'
! extract values for invariant equilibrium
                   invar=>mapline%end
                   if(ocv()) write(*,*)'Invariant eq: ',&
                        invar%seqx,invar%savednodeceq
                   if(invar%savednodeceq.lt.0) then
                      write(*,*)'SMP equilibrium not saved, skipping'
                      goto 222
                   endif
                   curceq=>results%savedceq(invar%savednodeceq)
!-------------------
! get the names of stable phases
                   kk=1
                   do jj=1,noph()
                      do ic=1,noofcs(jj)
                         k3=test_phase_status(jj,ic,value,curceq)
                         if(k3.gt.0) then
! this phase is stable or fix
                            call get_phase_name(jj,ic,dummy)
                            phaseline(nlinesep)(kk:)=dummy
                            kk=len_trim(phaseline(nlinesep))+2
                         endif
                      enddo
                   enddo
!                   write(*,117)nlinesep,phaseline(nlinesep)(1:kk-1)
!-------------------
! axis without wildcard
                   statevar=pltax(notanp)
                   call meq_get_state_varorfun_value(statevar,value,&
                        encoded1,curceq)
!                   write(*,*)'SMP axis variable 3: ',encoded1(1:3),value
                   if(gx%bmperr.ne.0) then
! THIRD skipping
                      write(*,212)'SMP skipping a point, error evaluating ',&
                           statevar,curceq%tpval(1),nv,0
                      goto 222
                   endif
                   nv=nv+3
                   if(nv.ge.maxval) then
                      write(*,*)'Too many points to plot 2',maxval
                      goto 600
                   endif
                   xax(nv-2)=value
                   xax(nv-1)=value
                   xax(nv)=value
!                   write(*,335)'New line: ',nlinesep,statevar(1:5),value
335                format(a,i3,' <',a,'> ',3(1pe14.6))
! axis with possible wildcard
                   statevar=pltax(anpax)
                   if(wildcard) then
! this cannot be a state variable derivative
!                    write(*,*)'Getting a wildcard value 2: ',nr,statevar(1:20)
                      call get_many_svar(statevar,yyy,nzp,np,encoded2,curceq)
                      if(gx%bmperr.ne.0) goto 1000
! save one non-zero value per line, 3 lines
                      ic=0
                      do jj=1,np
! we have put anp to zero above
!                         anp(jj,nv-2)=zero
!                         anp(jj,nv-1)=zero
!                         anp(jj,nv)=zero
                         if(yyy(jj).ne.zero) then
                            anp(ikol,nv-ic)=yyy(jj)
!                            write(*,7)nv-ic,jj,ikol,yyy(jj),anp(ikol,nv-ic)
7                           format('Nonzero value line ',i3,' column ',i3,&
                                 ' placed in column ',i3,2(1pe12.4))
                            ic=ic+1
                         endif
                      enddo
                   else
! if no wildcard there is no invariant line
!                      write(*,*)'It can be an invariant point here!!'
                      nv=nv-3
                      goto 225
!----------------------------------
! code below until label 222 redundant
                      np=1
                      call meq_get_state_varorfun_value(statevar,&
                           value,encoded1,curceq)
                      if(gx%bmperr.ne.0) then
! FOURTH Skipping never executed
                      write(*,212)'SMP skipping a point, error evaluating: ',&
                              statevar(1:10),curceq%tpval(1),nv,0
                         nv=nv-1; goto 222
                      endif
!                      if(gx%bmperr.ne.0) goto 1000
                      anp(1,nv)=value
!                      write(*,201)'at 225: ',nr,nv,curceq%tpval(1),value
                   endif
! this must be written on 3 lines terminated with an empty line
                   nlinesep=nlinesep+1
                   linesep(nlinesep)=nv
!                   write(*,*)'adding empty line 1',nlinesep,linesep(nlinesep)
                endif inv
! code from "goto 225" at 17 lines above until here is never used
!-------------------------------
! exit here if error and skip point
222             continue
             endif
          endif map1
! jump here if no wildcard
225       continue
       endif
!---- take next node along the same line
230    continue
!       write(*,*)'In ocplot2, looking for segmentation fault 4L'
       kk=seqx
       if(associated(mapline%end)) then
          seqx=mapline%end%seqx
       else
          seqx=0
       endif
240    continue
!       write(*,*)'Next node: ',seqx
       if(seqx.eq.0) then
          if(nlinesep.gt.0) then
             if(linesep(nlinesep).lt.nv) then
! we should never have several linesep for the same value of nv!
                nlinesep=nlinesep+1
                linesep(nlinesep)=nv
!                write(*,*)'adding empty line 2',nlinesep,linesep(nlinesep)
             endif
          endif
          if(line.eq.2) goto 500
          line=2
          goto 100
       else
          if(kk.eq.seqx) then
             if(giveup.gt.3) then
                write(*,*)'infinite loop ?'
                seqx=0
                goto 240
             endif
             giveup=giveup+1
          endif
          mapnode=>localtop%next
!          write(*,*)'In ocplot2, looking for segmentation fault 4M'
! loop through all mapnodes
250       continue
          if(mapnode%seqx.eq.seqx) then
! >>> this is just for step, for map one must find line connected
             mapline=>mapnode%linehead(1)
! skip line if EXCLUDEDLINE set
             if(.not.btest(mapline%status,EXCLUDEDLINE)) then
                if(mapline%done.eq.0) goto 110
             else
!                write(*,*)'Skipping a line 2'
             endif
!             if(mapline%done.eq.0) goto 110
          endif
          if(.not.associated(mapnode,localtop)) then
             mapnode=>mapnode%next
             goto 250
          else
! we have gone through all mapnodes without finding one with index seqx!!
!             write(*,*)'Cannot find node: ',seqx
! this seems not to be a problem .... probably already found.
             seqx=0; goto 240
          endif
       endif
!------------------------------------------- end loop 100
! check if we can find any lines not starting from localtop to be plotted
500    continue
!       write(*,*)'Checking for unplotted lines'
!       write(*,*)'In ocplot2, looking for segmentation fault 5'
       mapnode=>localtop%next
       do while(.not.associated(mapnode,localtop))
          jjline: do jj=1,mapnode%lines
             if(mapnode%linehead(jj)%done.eq.0) then
                if(ocv()) write(*,*)'Found a line in node: ',mapnode%seqx,jj
                line=jj
                if(nlinesep.gt.0) then
                   if(linesep(nlinesep).lt.nv) then
! we should never have several linesep for the same value of nv!
                      nlinesep=nlinesep+1
                      linesep(nlinesep)=nv
!                     write(*,*)'adding empty line:',nlinesep,linesep(nlinesep)
                   endif
                endif
                mapline=>mapnode%linehead(line)
!                goto 110
! skip line if EXCLUDEDLINE set
                if(btest(mapline%status,EXCLUDEDLINE)) then
!                   write(*,*)'Skipping a line 1'
!                   mapnode%linehead(line)%done=0
                   cycle jjline
                else
                   goto 110
                endif
             endif
          enddo jjline
          mapnode=>mapnode%next
!          write(*,*)'Looking at node: ',mapnode%seqx
       enddo
!--------------------------------------------
! end extracting data
600    continue
       overflow=.FALSE.
! but we may have another maptop !!
       if(associated(localtop%plotlink)) then
          if(.not.moretops) then
             write(*,*)'More than one maptop record'
             moretops=.true.
          endif
          localtop=>localtop%plotlink
          goto 77
       endif
!       write(*,*)'Number of points: ',nv
       if(.not.wildcard) then
          np=1; nrv=nv
!          write(*,*)'Extracted values: ',nrv
          goto 800
       endif
!============================================================
! remove columns that are only zeroes
!       write(*,*)'Now remove colums with just zeros',nv,nrv
!       read(*,17)ch1
       ic=0
! if a selected phase has been plotten np and qp may be different
! select the largest!
       nnp=max(np,qp)
!       write(*,*)'wildcard 3: ',wildcard,np,qp,nnp
!------------------------------------------ begin loop 650
650    ic=ic+1
660       continue
          if(ic.gt.nnp) goto 690
          if(nonzero(ic).eq.0) then
! shift all values from ic+1 to np
             if(nnp.ge.maxanp) then
                write(kou,*)'Too many points in anp array 1',maxanp,nv,nnp
                overflow=.TRUE.
                nnp=maxanp
             endif
             if(nv.ge.maxval) then
                write(kou,*)'Too many points in anp array 2',maxval,nv
                overflow=.TRUE.
                nv=maxval
             endif
             do jj=ic,nnp-1
                do nnv=1,nv
                   anp(jj,nnv)=anp(jj+1,nnv)
                enddo
                nonzero(jj)=nonzero(jj+1)
! also shift label
                lid(jj)=lid(jj+1)
             enddo
             nnp=nnp-1
             goto 660
          endif
! there is no more space in arrays to plot
          if(overflow) then
             write(*,*)'overflow',nv,nnp
             goto 690
          endif
          goto 650
!------------------------------------------ end loop 650
! nnp is the number of columns to plot
! nv is the number of separate lines
690 continue
       nrv=nv
       np=nnp
!       goto 800
!============================================ generate gnuplot file
800 continue
    write(*,808)np,nv,maxanp,maxval
808 format('plot data used: ',2i7,' out of ',2i7)
    if(np.eq.0) then
       write(kou,*)'No data to plot'
       gx%bmperr=4248
       goto 1000
    endif
!------------------------------------------------------------
! copy from lineheader to lid, each item separated by a space
!    allocate(lid(np))
! at present I have no idea what is in the lines ... just set numbers
!    kk=1
!    do i=1,np
!       call wriint(lineheader,kk,i)
!       kk=kk+1
!    enddo
!    kk=0
!    do i=1,np
!       call getext(lineheader,kk,1,lid(i),'x ',lhpos)
!    enddo
    if(.not.allocated(lid)) then
!       if(np.ge.1) then
! lid should always be allocated if np>1, but ... one never knows 
       allocate(lid(np))
       do i=1,np
          lid(i)='calculated '
       enddo
    endif
!------------------------------------------------------------
!    write(*,*)'We are here removing _ ',np
! replace _ by - in lid
    do i=1,np
! lid may contain phase names with _
! replace _ by - in lid because _ is interpreted as subscript (as LaTeX)
798    continue
!       write(*,*)'lid: ',i,': ',trim(lid(i))
       nv=index(lid(i),'_')
       if(nv.gt.0) then
          lid(i)(nv:nv)='-'
          goto 798
       endif
    enddo
! move data output to the end of PLT file ...
!    write(*,*)'We jump to 2000'    
    goto 2000
2000 continue
!    write(*,*)'We are at 2000 '
!----------------------------------------------------------------------
!
    call get_plot_conditions(encoded1,maptop%number_ofaxis,axarr,ceq)
!
    call ocplot2B(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,anp,lid,&
         title,filename,graphopt,pform,encoded1)
!    goto 900
! deallocate, not really needed for local arrays ??
    deallocate(anp)
    deallocate(xax)
    deallocate(linesep)
    if(allocated(yyy)) then
       deallocate(yyy)
       deallocate(nonzero)
    endif
1000 continue
    return
  end subroutine ocplot2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
  subroutine ocplot2B(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,anp,lid,&
       title,filename,graphopt,pform,conditions)
! called from icplot2 to generate the GNUPLOT file after extracting data
! np is number of columns (separate lines), if 1 no labelkey
! nrv is number of values to plot?
! nlinesep is the number of separate lines (index to linesep)
! linesep is the row when the line to plot finishes
! pltax
! xax array of values for single valued axis (T or mu etc)
! anpax=2 if axis with single value is column 2 and (multiple) values in
!         columns 3 and higher
! anp array of values for axis with multiple values (can be single values also)
! lid array with GNUPLOT line types for the different lines
! title Title of the plot
! filename GNUPLOT file name, (also used for pdf/ps/gif file)
! graphopt is graphical option record
! pform is output form (scree/acrobat/postscript/git)
! conditions is a character with the conditions for the diagram
    implicit none
    integer np,anpax,nlinesep
    integer ndx,nrv,linesep(*),anpdim
    character pltax(*)*(*),filename*(*),pform*(*),lid(*)*(*),title*(*)
    character conditions*(*)
    type(graphics_options) :: graphopt
    double precision xax(*),anp(anpdim,*)
    type(graphics_textlabel), pointer :: textlabel
!\end{verbatim}
!----------------------------------------------------------------------
! internal
    integer ii,jj,kk,lcolor,appfil,nnv,ic,repeat,ksep,nv,k3,kkk
    character pfc*64,pfh*64,backslash*2,appline*128
    character applines(20)*128,gnuplotline*80,labelkey*64
! write the gnuplot command file with data appended
!
!    write(*,10)'in ocplot2B: ',np,anpax,nrv,pform(1:1),trim(title),&
!         nlinesep,(linesep(kk),kk=1,nlinesep)
10  format(a,3i5,' "',a,'" '/a/i3,2x,15i4)
    if(index(filename,'.plt ').le.0) then 
       kk=len_trim(filename)
       pfc=filename(1:kk)//'.'//'plt '
       kkk=kk+4
    endif
    open(21,file=pfc,access='sequential',status='unknown')
    write(21,1600)
1600 format('# GNUPLOT file generated by Open Calphad'/)
! if there is just one curve do not write any key.  May be overriiden later ..
    if(pform(1:1).eq.'P') then
       pfh=filename(1:kk)//'.'//'ps '
       write(21,850)pfh(1:len_trim(pfh))
850    format('set terminal postscript color'/'set output "',a,'"')
    elseif(pform(1:1).eq.'A') then
       pfh=filename(1:kk)//'.'//'pdf '
       write(21,851)pfh(1:len_trim(pfh))
851    format('set terminal pdf color'/'set output "',a,'"')
    elseif(pform(1:1).eq.'G') then
       pfh=filename(1:kk)//'.'//'gif '
       write(21,852)pfh(1:len_trim(pfh))
852    format('set terminal gif'/'set output "',a,'"')
    endif
! this part is independent of which axis is a single value
!------------------ some GNUPLOT colors:
! colors are black: #000000, red: #ff000, web-green: #00C000, web-blue: #0080FF
! dark-yellow: #C8C800, royal-blue: #4169E1, steel-blue #306080,
! gray: #C0C0C0, cyan: #00FFFF, orchid4: #804080, chartreuse: 7CFF40
! if just one line set key off for that line.
    if(np.eq.1 .and. graphopt%appendfile(1:1).eq.' ') then
       labelkey=' off'
    else
       labelkey=graphopt%labelkey
    endif
    write(21,860)trim(title),trim(conditions),&
         trim(pltax(1)),trim(pltax(2)),labelkey
860 format('set title "',a,' \n ',a,'"'/&
            'set xlabel "',a,'"'/'set ylabel "',a,'"'/&
            'set key ',a/&
            'set style line 1 lt 2 lc rgb "#000000" lw 2'/&
            'set style line 2 lt 2 lc rgb "#FF0000" lw 2'/&
            'set style line 3 lt 2 lc rgb "#00C000" lw 2'/&
            'set style line 4 lt 2 lc rgb "#0080FF" lw 2'/&
            'set style line 5 lt 2 lc rgb "#C8C800" lw 2'/&
            'set style line 6 lt 2 lc rgb "#4169E1" lw 2'/&
            'set style line 7 lt 2 lc rgb "#C0C0C0" lw 2'/&
            'set style line 8 lt 2 lc rgb "#00FFFF" lw 2'/&
            'set style line 9 lt 2 lc rgb "#804080" lw 2'/&
            'set style line 10 lt 2 lc rgb "#7CFF40" lw 2')
!
    if(graphopt%rangedefaults(1).ne.0) then
! user defined ranges for x axis
       write(21,870)'x',graphopt%plotmin(1),graphopt%plotmax(1)
870    format('set ',a1,'range [',1pe12.4,':',1pe12.4,'] ')
    endif
    if(graphopt%rangedefaults(2).ne.0) then
! user defined ranges for y axis
       write(21,870)'y',graphopt%plotmin(2),graphopt%plotmax(2)
    endif
!----------------------
! logarithmic axis
    if(graphopt%axistype(1).eq.1) then
       write(21,151)'x'
151    format('set logscale ',a)
    endif
    if(graphopt%axistype(2).eq.1) then
       write(21,151)'y'
    endif
!------------------------------------------------------------
! set labels
    textlabel=>graphopt%firsttextlabel
    do while(associated(textlabel))
       jj=len_trim(textlabel%textline)
       write(21,1505)textlabel%textline(1:jj),textlabel%xpos,textlabel%ypos
1505   format('set label "',a,'" at ',1pe14.6,', ',1pe14.6)
!       write(*,1505)textlabel%textline(1:jj),textlabel%xpos,textlabel%ypos
       textlabel=>textlabel%nexttextlabel
    enddo
!---------------------------------------------------------------
! handle appended files here ....
!
    if(graphopt%appendfile(1:1).eq.' ') then
       appfil=0
    else
       appfil=23
       write(kou,*)'Appending data from: ',trim(graphopt%appendfile)
       open(appfil,file=graphopt%appendfile,status='old',&
            access='sequential',err=1750)
!
       write(21,1720)'# APPENDED from '//trim(graphopt%appendfile)
! copy all lines up to "plot" to new graphics file
       nnv=0
1710   continue
       read(appfil,1720,end=1750)appline
1720   format(a)
!------------------------------------------------------------------
! ignore some lines with "set" in the append file
! set title
! set xlabel
! set ylabel
! set xrange
! set yrange
! set key
       if(appline(1:10).eq.'set title ' .or.&
            appline(1:11).eq.'set xlabel ' .or.&
            appline(1:11).eq.'set ylabel ' .or.&
            appline(1:11).eq.'set xrange ' .or.&
            appline(1:11).eq.'set yrange ' .or.&
            appline(1:8).eq. 'set key ') then
!          write(*,*)'ignoring append line ',trim(appline)
          goto 1710
       endif
!------------------------------------------------------------------
       if(index(appline,'plot "-"').gt.0) then
          applines(1)=appline
          ic=1
1730      continue
! if line ends with \ then read more
          ii=len_trim(appline)
!          write(*,*)'There are more? ',appline(i:ii),ii,ic
          if(appline(ii:ii).eq.'\') then
! continuation lines
             read(appfil,1720,end=1750)appline
             ic=ic+1
             applines(ic)=appline
             goto 1730
          endif
!          write(*,*)'appline: ',trim(appline),ic
!          close(appfil)
!          appfil=0
          goto 1770
       endif
       write(21,1720)trim(appline)
       nnv=nnv+1
       goto 1710
! error oppening append file
1750   continue
       write(kou,*)'Error opening or reading the append file, skipping it'
       close(appfil)
       appfil=0
1770   continue
    endif
!---------------------------------------------------------------
    backslash=',\'
    lcolor=1
!    write(*,*)'backslash "',backslash,'" '
    if(anpax.eq.2) then
! anpax=2 means the single valued axis is colum 2 and possibly multiple
! values in column 3 and higher
! this is the multiple plot, file name only given once!!
! last line tuple on a separate format statement, if np>2
! np is number of columns
       if(np.eq.1 .and. appfil.eq.0) then
          write(21,880)lcolor,' ',' '
       else
          write(21,880)lcolor,trim(lid(1)),backslash
880       format('plot "-" using 2:3 with lines ls ',i2,' title "'a,'"',a)
       endif
       do ii=2,np-1
          lcolor=lcolor+1
          if(lcolor.gt.10) lcolor=1
          write(21,882)ii+2,lcolor,trim(lid(ii)),backslash
882       format('"" using 2:',i3,' with lines ls ',i2,' title "'a,'"',a)
       enddo
       lcolor=lcolor+1
       if(lcolor.gt.10) lcolor=1
       if(appfil.eq.0) then
          if(np.ge.2) write(21,882)np+2,lcolor,trim(lid(np)),' '
       else
! write the last calculated curve if np>1
          if(np.ge.2) write(21,882)ii+2,lcolor,trim(lid(ii)),backslash
! we should append data, change plot "-" to just "" in appline(1)
          ii=index(applines(1),'plot "-"')
          applines(1)(1:ii+7)='""'
          do ii=1,ic
             write(21,884)trim(applines(ii))
884          format(a)
          enddo
       endif
    else
! anpax=2 means the single valued axis is colum 2
! this writes the file name only once and last line separate if np>2
       if(np.eq.1 .and. appfil.eq.0) then
          write(21,890)lcolor,' ',' '
       else
          write(21,890)lcolor,trim(lid(1)),backslash
890       format('plot "-" using 3:2 with lines ls ',i2,' title "',a,'"',a)
       endif
       lcolor=2
       do ii=2,np-1
          write(21,892)ii+2,lcolor,trim(lid(ii)),backslash
892       format('"" using ',i3,':2 with lines ls ',i2,' title "'a,'"',a)
          lcolor=lcolor+1
          if(lcolor.gt.10) lcolor=1
       enddo
       if(appfil.eq.0) then
          if(np.ge.2) write(21,892)np+2,lcolor,trim(lid(np)),' '
       else
! write the last calculated curve if np>1
          if(np.ge.2) write(21,892)ii+2,lcolor,trim(lid(ii)),backslash
! we should append data, change plot "-" to just "" in appline(1)
          ii=index(applines(1),'plot "-"')
          applines(1)(1:ii+7)='""'
          do ii=1,ic
             write(21,884)trim(applines(ii))
          enddo
       endif
    endif
!------------------------------------------------------------------
! Add output of data as many times as there are lines above
! if anpax=1 then we must put the first colum after the colon in gnuplot
    repeat=0
1800 continue
    ksep=2
    if(nlinesep.lt.2) linesep(2)=nrv
    repeat=repeat+1
    jj=0
!    write(*,*)'Writing repeat, rows, columns ',repeat,nrv,np+2
    do nv=1,nrv
       write(21,1820)nv,xax(nv),(anp(jj,nv),jj=1,np)
1820    format(i4,1000(1pe14.6))
       if(nv.eq.linesep(ksep)) then
! an empty line in the dat file means a MOVE to the next point.
          if(nv.lt.nrv) then
             if(jj.gt.1) then
! phaseline is never used in the plot, just included in the file
                write(21,1819)ksep
1819            format('# end of line '//'# Line ',i5,&
                     ', with these stable phases:'/'# ',a)
             else
                write(21,1822)ksep
1822            format('# end of line '//'# Line ',i5,&
                     ', with unknown phases')
             endif
          else
! try to avoid rubbish
             write(21,1821)
1821         format('# end of line '/)
          endif
! test of uninitiallized variable, ksep must not exceed nlinesep
          ksep=min(ksep+1,nlinesep)
       endif
    enddo
    if(repeat.lt.np) then
! we must have an e and repeat the data output for each line to plot
! eventually we may select just the columns with data used !!
       write(21,1833)repeat+1
1833   format('e '//'#--- repeat lines ',i3,' ------------------------------')
       goto 1800
    else
! we must have a simple "e" at the end of last
       write(21,1834)
1834   format('e ')
    endif
!-----------------------------------------------------
! finally copy the data from the append file, it should be correctly formatted
    if(appfil.gt.0) then
       ic=0
1900   continue
       read(appfil,884,end=1910)appline
       ic=ic+1
       if(appline(1:11).eq.'pause mouse') goto 1900
       write(21,884)trim(appline)
       goto 1900
1910   continue
!       write(*,*)'Appended ',ic,' data lines'
       close(appfil)
       appfil=0
    endif
!------------------------------------------------------
    if(pform(1:1).eq.' ') then
! if not hardcopy pause gnuplot.  Mouse means clicking in the graphics window
! will close it. I would like to have an option to keep the graphics window...
       write(21,990)
990    format('pause mouse')
    endif
    close(21)
    if(appfil.ne.0) close(appfil)
    appfil=0
!-------------------------------------------------------------------
! execute the GNUPLOT command file
!
    gnuplotline='gnuplot '//pfc(1:kkk)//' & '
! if gnuplot cannot be started with gnuplot give normal path ...
!    gnuplotline='"c:\program files\gnuplot\bin\wgnuplot.exe" '//pfc(1:kkk)//' '
    k3=len_trim(gnuplotline)+1
    write(kou,*)'Gnuplot command file: ',pfc(1:kk+4)
    if(pform(1:1).ne.' ') then
       write(kou,*)'Graphics output file: ',pfh(1:kk+4)
    endif
! call system without initial "gnuplot " keeps the window !!!
    if(btest(graphopt%status,GRKEEP)) then
!       write(*,*)'plot command: "',gnuplotline(9:k3),'"'
       write(*,*)'Trying to spawn'
       call system(gnuplotline(9:k3))
    else
!       write(*,*)'plot command: "',gnuplotline(1:k3),'"'
       call system(gnuplotline)
    endif
1000 continue
    return
  end subroutine ocplot2B

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine ocplot3(ndx,pltax,filename,maptop,axarr,graphopt,pform,ceq)
! special to plot isothermal sections (two columns like x(*,cr) x(*,ni))
! ndx is mumber of plot axis, 
! pltax is text with plotaxis variables
! filename is intermediary file (maybe not needed)
! maptop is map_node record with all results
! axarr are axis records
! graphopt is graphics record (should be extended to include more)
! pform is graphics form
! pform is type of output (screen or postscript or gif)
    implicit none
    integer ndx
    character pltax(*)*(*),filename*(*),pform*(*)
    type(map_axis), dimension(*) :: axarr
    type(map_node), pointer :: maptop
    type(graphics_options) :: graphopt
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
    integer, parameter :: maxsame=200
    type(map_node), pointer :: plottop,curtop,endnode
    type(map_line), pointer :: curline
    type(gtp_equilibrium_data), pointer :: curceq
    type(map_ceqresults), pointer :: results
    integer ii,jj,point,plotp,lines,eqindex,lasteq,nooftup,lokcs,jph,same,kk
    integer, parameter :: maxval=2000,mazval=100
    double precision, allocatable :: xval(:,:),yval(:,:),zval(:,:),tieline(:,:)
    integer offset,nofeq,sumpp,last,nofinv,ntieline,mtieline
    double precision xxx,yyy
    integer, allocatable :: plotkod(:),lineends(:)
    character xax1*8,xax2*24,yax1*8,yax2*24,axis1*32,axisx(2)*32,axisy(2)*32
    character phname*32,encoded*1024,axis*32
    character lid(2,maxsame)*24
!
! xval and yval and ccordinates to plot, 
! points on one line is       xval(1,jj),yval(1,jj)
! points on the other line is xval(2,jj),yval(2,jj)
! zval is an occational point zval(1,kk),zval(2,kk) at line ends (invariants)
! plotkod is a specific code for each point (not used?)
! lineends is the index in xval/yval for an end of line
    allocate(xval(2,maxval))
    allocate(yval(2,maxval))
    allocate(plotkod(maxval))
    allocate(zval(2,maxsame))
    allocate(lineends(maxsame))
    if(.not.associated(maptop)) then
       write(*,*)'No data to plot'
       goto 1000
    endif
    write(*,17)
17  format(//'Using ocplot3')
! extract the axis variables
    jph=index(pltax(1),'*')
    xax1=pltax(1)(1:jph-1)
    xax2=pltax(1)(jph+1:)
    jph=index(pltax(2),'*')
    yax1=pltax(2)(1:jph-1)
    yax2=pltax(2)(jph+1:)
!
! initiate loop to extract values
    point=0
    plotp=0
    nofinv=0
! if graphopt%tielines not zero check if the tielines are in plane ...
! %tieline_tieline_inplane <0 means step, 0 means isopleth
    if(graphopt%tielines.gt.0) then
       if(maptop%tieline_inplane.le.0) then
          write(kou,*)'No tie-lines can be plotted'
          graphopt%tielines=0
       endif
    endif
! same is incremented for each line
    same=0
    sumpp=0
    plottop=>maptop
    curtop=>plottop
    nooftup=noofphasetuples()
100 continue
    if(.not.allocated(curtop%linehead)) goto 500
    lines=size(curtop%linehead)
    results=>plottop%saveceq
    node: do ii=1,lines
       curline=>curtop%linehead(ii)
!       write(*,*)'Data from line: ',curline%lineid,lines
       eqindex=curline%first
       lasteq=curline%last
       if(eqindex.le.0 .or. eqindex.gt.lasteq) cycle node
       last=same
       same=same+1
       if(same.gt.maxsame) then
          write(*,*)'Too many lines to plot ',maxsame
          cycle node
       endif
       nofeq=lasteq+1-eqindex
       axisx=' '
       axisy=' '
       ntieline=0
!       write(*,*)'Plotting tie-lines: ',graphopt%tielines
       if(graphopt%tielines.gt.0) then
! the number of tie-lines to extract
          mtieline=nofeq/graphopt%tielines
!          write(*,*)'Number of tie-lines: ',mtieline
          allocate(tieline(4,mtieline+5))
!          write(*,*)'Allocating for tielines: ',mtieline+1
! UNFINISHED: try to have equal number of equilibria at the beginning and end 
       endif
       line: do eqindex=eqindex,lasteq
! extract for each stable phase the state variable in pltax       
          point=point+1
          curceq=>results%savedceq(eqindex)
! find the stable phases (max 3)
          plotp=plotp+1
          if(plotp.gt.maxval) then
             write(*,*)'Too many points to plot ',maxval
             cycle node
          endif
          jj=1
          equil: do jph=1,nooftup
             lokcs=phasetuple(jph)%lokvares
             if(curceq%phase_varres(lokcs)%phstate.ge.PHENTSTAB) then
                if(jj.ge.3) then
                   if(eqindex.eq.lasteq) then
! skip this point if it is the last
                      plotp=plotp-1
!                      write(*,*)'ocplot3 indexing error, skipping last point'
                   endif
                   cycle equil
                endif
! the generation of the axis state variable is needed just once for all points 
! if we save the axis text ...
                plotkod(plotp)=same
                call get_phasetup_name(jph,phname)
                if(same.gt.last) then
                   lid(jj,same)=phname
                endif
                if(axisx(1)(1:1).eq.' ') then
                   axisx(1)=trim(xax1)//trim(phname)//trim(xax2)
                   axisy(1)=trim(yax1)//trim(phname)//trim(yax2)
                elseif(axisx(2)(1:1).eq.' ') then
                   axisx(2)=trim(xax1)//trim(phname)//trim(xax2)
                   axisy(2)=trim(yax1)//trim(phname)//trim(yax2)
                endif
                call meq_get_state_varorfun_value(axisx(jj),xxx,encoded,curceq)
                xval(jj,plotp)=xxx
                call meq_get_state_varorfun_value(axisy(jj),xxx,encoded,curceq)
                yval(jj,plotp)=xxx
!                write(*,19)'X/Y axis variable: ',plotp,trim(axis),&
!                     xval(jj,plotp),xxx,point
!19              format(a,i5,2x,a,2F10.6,2i5)
                jj=jj+1
             endif
          enddo equil
          if(graphopt%tielines.gt.0) then
! exact coordinates for tielines each ntieline equilibria
!             write(*,*)'saving tieline?',eqindex,&
!                  mod(eqindex,graphopt%tielines),ntieline,plotp
             if(mod(eqindex,graphopt%tielines).eq.0) then
                ntieline=ntieline+1
                tieline(1,ntieline)=xval(1,plotp)
                tieline(2,ntieline)=yval(1,plotp)
                tieline(3,ntieline)=xval(2,plotp)
                tieline(4,ntieline)=yval(2,plotp)
             endif
          endif
!          write(*,23)'phase line: ',same,last,trim(lid(1,same)),&
!               trim(lid(2,same))
          last=same
!          write(*,21)xval(1,plotp),yval(1,plotp),&
!               xval(2,plotp),yval(2,plotp),plotp
!21        format('phase 1: ',2F10.7,10x,'phase 2: ',2F10.7,i7)
          lineends(same)=plotp
       enddo line
! check if line ends in a node and add its coordinates for the same phases
       endnode=>curline%end
       if(associated(endnode)) then
! there is a nod at the end, extracts its ceq record
          curceq=>endnode%nodeceq
          plotp=plotp+1
          do jj=1,2
             call meq_get_state_varorfun_value(axisx(jj),xxx,encoded,curceq)
             if(gx%bmperr.ne.0) then
                write(*,*)'Error extracting end points'
                goto 1000
             endif
             xval(jj,plotp)=xxx
             call meq_get_state_varorfun_value(axisy(jj),yyy,encoded,curceq)
             if(gx%bmperr.ne.0) goto 1000
             yval(jj,plotp)=yyy
          enddo
!          write(*,*)'Endnode x,y: ',plotp,xxx,yyy
! correct lineends!!
          lineends(same)=plotp
       endif
!       write(*,23)'phase line: ',same,plotp,trim(lid(1,same)),trim(lid(2,same))
23     format(a,2i5,3x,a,' and ',a)
       if(ntieline.gt.0) then
! All tielines on the same line with a space in between
          same=same+1
          do eqindex=1,ntieline
             plotp=plotp+1
             xval(1,plotp)=tieline(1,eqindex)
             yval(1,plotp)=tieline(2,eqindex)
! this means the tie-lines will be plotted twice ... but mhy not??
             xval(2,plotp)=tieline(1,eqindex)
             yval(2,plotp)=tieline(2,eqindex)
             plotkod(plotp)=-100
             plotp=plotp+1
             xval(1,plotp)=tieline(3,eqindex)
             yval(1,plotp)=tieline(4,eqindex)
             xval(2,plotp)=tieline(3,eqindex)
             yval(2,plotp)=tieline(4,eqindex)
             lineends(same)=plotp
             plotkod(plotp)=-101
          enddo
          lid(1,same)='tieline'
          lid(2,same)='tieline'
       endif
! no longer any use of tieline
       if(allocated(tieline)) deallocate(tieline)
    enddo node
!----------------------------------------------------------------------
! finished all lines in this curtop, take next
! but first generate the invariant
    curtop=>curtop%next
    if(.not.associated(curtop,maptop)) then
!       write(*,*)'Extracting data from node'
       curceq=>curtop%nodeceq
       if(associated(curceq)) then
!          write(*,*)'Extracting data from node equilibrium'
          plotp=plotp+1
          jj=1
          same=same+1
          nodeequil: do jph=1,nooftup
             lokcs=phasetuple(jph)%lokvares
             if(curceq%phase_varres(lokcs)%phstate.ge.PHENTSTAB) then
! plotkod set to -1 to indicate invariant
                plotkod(plotp)=-1
                call get_phasetup_name(jph,phname)
!                write(*,*)'Stable phase ',trim(phname),jj
                if(jj.lt.3) lid(jj,same)='invariant'
                axis=trim(xax1)//trim(phname)//trim(xax2)
                call meq_get_state_varorfun_value(axis,xxx,encoded,curceq)
                if(jj.ge.3) then
                   nofinv=nofinv+1
                   zval(1,nofinv)=xxx
                else
                   xval(jj,plotp)=xxx
                endif
                axis=trim(yax1)//trim(phname)//trim(yax2)
                call meq_get_state_varorfun_value(axis,xxx,encoded,curceq)
                if(jj.ge.3) then
                   zval(2,nofinv)=xxx
                else
                   yval(jj,plotp)=xxx
                endif
                jj=jj+1
             endif
          enddo nodeequil
!          write(*,333)'invariant:',nofinv,plotp,xval(1,plotp),yval(1,plotp),&
!               xval(2,plotp),yval(2,plotp),zval(1,nofinv),zval(2,nofinv)
!333       format(a,i3,i4,3(2F8.4,2x))
! create a line for just this point
!          same=same+1
          lineends(same)=plotp
       endif
       goto 100
    endif
!    do jj=1,same
!       write(*,23)'phases: ',same,jj,trim(lid(1,jj)),trim(lid(2,jj))
!    enddo
!------------------------------------------------
! there can be more maptops linked via plotlink
    if(associated(plottop%plotlink)) then
       plottop=>plottop%plotlink
       goto 100
    endif
!========================================================
    call get_plot_conditions(encoded,maptop%number_ofaxis,axarr,ceq)
! now we should have all data to plot in xval and yval arrays
500 continue
!    write(*,*)'found lines/points to plot: ',same,plotp
!    write(*,502)(lineends(ii),ii=1,same)
502 format(10i5)
    call ocplot3B(same,nofinv,lineends,2,xval,2,yval,2,zval,plotkod,pltax,&
         lid,filename,graphopt,pform,encoded)
    deallocate(xval)
    deallocate(yval)
    deallocate(plotkod)
1000 continue
    return
  end subroutine ocplot3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
  subroutine ocplot3B(same,nofinv,lineends,nx1,xval,ny1,yval,nz1,zval,plotkod,&
       pltax,lid,filename,graphopt,pform,conditions)
! called by ocplot3 to write the GNUPLOT file for two wildcard columns
! same is the number of lines to plot
! nofinv number of invariants
! lineends array with row numbers where ecah line ends
! nx1 first dimension of xval
! xval 2D matrix with values to plot on x axis
! ny1 first dimension of yval
! yval 2D matrix with values to plot on y axis
! nz1 first dimension of zval
! zval 2D matrix with third point of invariant triangles
! plotkod integer array indicating the type of line (-1 skip line)
! pltax text for axis
! lid array with GNUPLOT line types
! filename is intermediary file (maybe not needed)
! graphopt is graphics option record
! maptop is map_node record with all results
! pform is type of output (screen/acrobat/postscript/gir)
! conditions is a text with conditions for the calculation
    implicit none
    character pltax(*)*(*),filename*(*),pform*(*),lid(nx1,*)*(*),conditions*(*)
    type(graphics_options) :: graphopt
    integer same,plotkod(*),nx1,ny1,nz1,nofinv
    integer lineends(*)
    double precision xval(nx1,*),yval(ny1,*),zval(nz1,*)
!\end{verbatim}
    integer, parameter :: maxcolor=200
    integer ii,jj,kk,jph,offset,n1
    type(graphics_textlabel), pointer :: textlabel
    character gnuplotline*64,date*12,mdate*12,title*128,deftitle*64,backslash*2
    character labelkey*24,applines(10)*128,appline*128,pfc*80,pfh*80
    integer sumpp,np,appfil,ic,nnv,kkk,lcolor(maxcolor),iz,again
    integer done(maxcolor),foundinv,fcolor,k3
    character color(maxcolor)*24
!  
    write(*,*)'Using the rudimentary graphics in ocplot3B!'
!
    call date_and_time(date)
    mdate=" "//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//" "
    deftitle='Open Calphad 4.0 prerelease '//mdate//': with GNUPLOT'
    if(graphopt%labeldefaults(1).eq.0) then
       title=deftitle
    else
! alwas inlcude open calphad and date, add user title at the end
       title=trim(deftitle)//' '//graphopt%plotlabels(1)
    endif
! np should be the number of different lines to be plotted
! if there is just one line do not write any key.  May be overriiden later ..
    np=same
    if(np.eq.1 .and. graphopt%appendfile(1:1).eq.' ') then
       labelkey=' off'
    else
       labelkey=graphopt%labelkey
    endif
!
    if(index(filename,'.plt ').le.0) then 
       kk=len_trim(filename)
       pfc=filename(1:kk)//'.'//'plt '
       kkk=kk+4
    else
       pfc=filename
    endif
    write(*,*)'filename: ',trim(pfc),' ',pform(1:1)
    open(21,file=pfc,access='sequential',status='unknown')
    write(21,110)
110 format('# GNUPLOT file generated by Open Calphad'/)
    if(pform(1:1).eq.'P') then
       pfh=filename(1:kk)//'.'//'ps '
       write(21,120)trim(pfh)
120    format('set terminal postscript color'/'set output "',a,'"')
    elseif(pform(1:1).eq.'A') then
       pfh=filename(1:kk)//'.'//'pdf '
       write(*,121)trim(pfh)
       write(21,121)trim(pfh)
121    format('set terminal pdf color'/'set output "',a,'"')
    elseif(pform(1:1).eq.'G') then
       pfh=filename(1:kk)//'.'//'gif '
       write(21,122)trim(pfh)
122    format('set terminal gif'/'set output "',a,'"')
    endif
!    open(21,file='ocgnu.plt ',access='sequential ',status='unknown ')
!
    write(21,130)trim(title),trim(conditions),&
         trim(pltax(1)),trim(pltax(2)),labelkey
130 format('set title "',a,' \n ',a,'"'/&
            'set xlabel "',a,'"'/'set ylabel "',a,'"'/&
            'set key ',a/&
            'set style line 1 lt 2 lc rgb "#000000" lw 2'/&
            'set style line 2 lt 2 lc rgb "#804080" lw 2'/&
            'set style line 3 lt 2 lc rgb "#00C000" lw 2'/&
            'set style line 4 lt 2 lc rgb "#FF0000" lw 2'/&
            'set style line 5 lt 2 lc rgb "#0080FF" lw 2'/&
            'set style line 6 lt 2 lc rgb "#C8C800" lw 2'/&
            'set style line 7 lt 2 lc rgb "#4169E1" lw 2'/&
            'set style line 8 lt 2 lc rgb "#7CFF40" lw 2'/&
            'set style line 9 lt 2 lc rgb "#C0C0C0" lw 2'/&
            'set style line 10 lt 2 lc rgb "#00FFFF" lw 2'/&
            'set style line 11 lt 2 lc rgb "#8F8F8F" lw 3'/&
            'set style line 12 lt 2 lc rgb "#8F8F8F" lw 1')
! The last two styles (11 and 12) are for invariants and tielines
!
! ranges for x and y
    if(graphopt%rangedefaults(1).ne.0) then
! user defined ranges for x axis
       write(21,150)'x',graphopt%plotmin(1),graphopt%plotmax(1)
150    format('set ',a1,'range [',1pe12.4,':',1pe12.4,'] ')
    endif
    if(graphopt%rangedefaults(2).ne.0) then
! user defined ranges for y axis
       write(21,150)'y',graphopt%plotmin(2),graphopt%plotmax(2)
    endif
!----------------------
! logarithmic axis
    if(graphopt%axistype(1).eq.1) then
       write(21,151)'x'
151    format('set logscale ',a)
    endif
    if(graphopt%axistype(2).eq.1) then
       write(21,151)'y'
    endif
!----------------------
! line labels
! set labels
    textlabel=>graphopt%firsttextlabel
    do while(associated(textlabel))
       jj=len_trim(textlabel%textline)
       write(21,160)textlabel%textline(1:jj),textlabel%xpos,textlabel%ypos
160    format('set label "',a,'" at ',1pe14.6,', ',1pe14.6)
!       write(*,1505)textlabel%textline(1:jj),textlabel%xpos,textlabel%ypos
       textlabel=>textlabel%nexttextlabel
    enddo
!---------------------------------------------------------------
! handle heading of appended files here ....
!
    if(graphopt%appendfile(1:1).eq.' ') then
       appfil=0
    else
       appfil=23
       write(kou,*)'Appending data from: ',trim(graphopt%appendfile)
       open(appfil,file=graphopt%appendfile,status='old',&
            access='sequential',err=280)
!
       write(21,210)'# APPENDED from '//trim(graphopt%appendfile)
! copy all lines up to "plot" to new graphics file
       nnv=0
200   continue
       read(appfil,210,end=290)appline
210   format(a)
!------------------------------------------------------------------
! ignore some lines with "set" in the append file
! set title
! set xlabel
! set ylabel
! set xrange
! set yrange
! set key
       if(appline(1:10).eq.'set title ' .or.&
            appline(1:11).eq.'set xlabel ' .or.&
            appline(1:11).eq.'set ylabel ' .or.&
            appline(1:11).eq.'set xrange ' .or.&
            appline(1:11).eq.'set yrange ' .or.&
            appline(1:8).eq. 'set key ') then
!          write(*,*)'ignoring append line ',trim(appline)
          goto 200
       endif
!------------------------------------------------------------------
       if(index(appline,'plot "-"').gt.0) then
          applines(1)=appline
          ic=1
230       continue
! if line ends with \ then read more
          ii=len_trim(appline)
!          write(*,*)'There are more? ',appline(i:ii),ii,ic
          if(appline(ii:ii).eq.'\') then
! continuation lines
             read(appfil,210,end=290)appline
             ic=ic+1
             applines(ic)=appline
             goto 230
          endif
          goto 290
       endif
       write(21,210)trim(appline)
       nnv=nnv+1
       goto 200
! error oppening append file
280    continue
       write(kou,*)'Error opening or reading the append file, skipping it'
       close(appfil)
       appfil=0
290    continue
    endif
! coordinate the content of lid with the colors
!    do ii=1,same
!       write(*,*)'phases: ',ii,trim(lid(1,ii)),' ',trim(lid(2,ii))
!    enddo
    nnv=0
    iz=0
    color(1)=' '
!    if(2*same.gt.maxcolor) then
!       write(*,*)'Number of lines: ',2*same
!    endif
    pair: do jj=1,2
       point: do ii=1,same
          iz=iz+1
! plotkod -1 negative means ignore
! plotkod -100 and -101 used for tie-lines
          if(jj.eq.2 .and. plotkod(iz).eq.-1) cycle pair
             do ic=1,nnv
                if(trim(lid(jj,ii)).eq.trim(color(ic))) then
                   if(iz.gt.maxcolor) then
                      write(kou,*)'lcolor dimension overflow',iz
                   else
                      lcolor(iz)=ic
                   endif
                   goto 295
                endif
             enddo
! no match, increment nnv and assign that color to lcolor
! skip colors 11 and 12, reserved for invariants and tielines
          nnv=nnv+1
          lcolor(iz)=nnv
          color(nnv)=lid(jj,ii)
!          write(*,293)'Assigned: ',nnv,jj,ii,iz,trim(color(nnv))
!293       format(a,4i5,' "',a,'"')
295       continue
       enddo point
    enddo pair
! replace _ by - (in phase names)
    do kk=1,nnv
297    continue
       jj=index(color(kk),'_')
       if(jj.gt.0) then
          color(kk)(jj:jj)='-'
          goto 297
       endif
    enddo
!---------------------------------------------------------------
! check for invariant and tieline and replace color!
!    do ii=1,2*same
!       write(*,*)'original: ',ii,lcolor(ii),trim(color(lcolor(ii)))
!    enddo
    lcolor1: do ii=1,2*same
       jj=lcolor(ii)
       if(trim(color(jj)).eq.'invariant') then
!          write(*,*)'Changing ',ii,jj,' to 11'
          lcolor(ii)=11
          color(11)='invariant'
          do kk=ii+1,2*same
             if(lcolor(kk).eq.jj) then
!                write(*,*)'subsequent: ',kk,lcolor(kk),jj,11
                lcolor(kk)=11
             endif
          enddo
          exit lcolor1
       endif
    enddo lcolor1
    lcolor2: do ii=1,2*same
       jj=lcolor(ii)
       if(trim(color(jj)).eq.'tieline') then
          lcolor(ii)=12
          color(12)='tie-line'
          do kk=ii+1,2*same
             if(lcolor(kk).eq.jj) lcolor(kk)=12
          enddo
          exit lcolor2
       endif
    enddo lcolor2
!    do ii=1,2*same
!       write(*,*)'Final: ',ii,lcolor(ii),trim(color(lcolor(ii)))
!    enddo
!----------------------------------------------------------------
! Here we generate the datafile with coordinates to plot
! if nx1 or ny1 is 1 plot all on other axis versus single axis coordinate
! if nx1=ny1 plot the pairs xval(1..nx1,jj) yval(1..ny1,jj)
! else error
!----------------------
    backslash=',\'
! empty line before the plot command
!    write(*,*)'lines: ',same,nofinv
    write(21,*)
    ii=0
    done=0
    do kk=1,2
       do jj=1,same
          ii=ii+1
          if(ii.eq.1) then
             write(21,309)lcolor(ii),trim(color(lcolor(ii))),backslash
309          format('plot "-" using 1:2 with lines ls ',i2,' title "',a,'"',a)
             done(lcolor(1))=1
          else
! the invariant "lines" are just a point and occur only once
! the last line for the plot command has no backslash
!             if(ii.eq.2*(same-nofinv)+nofinv) backslash=' '
             if(ii.eq.2*same) backslash=' '
! we can only use linestyles 1 to 10 except for invariants and tie-lines
             fcolor=lcolor(ii)
             if(fcolor.gt.12) then
                fcolor=mod(fcolor,10)
                if(fcolor.eq.0) fcolor=10
             endif
             if(done(lcolor(ii)).eq.1) then
                write(21,320)fcolor,backslash
320             format('"" using 1:2 with lines ls ',i2,' notitle ',a)
             else
                write(21,310)fcolor,trim(color(lcolor(ii))),backslash
310             format('"" using 1:2 with lines ls ',i2,' title "',a,'"',a)
                done(lcolor(ii))=1
             endif
          endif
       enddo
    enddo
!    goto 500
! loop for all line coordinates
!    sumpp=0
    do kk=1,2
!       again=sumpp
       foundinv=0
       sumpp=0
       do jj=1,same
! handle invariants, just one point, just once
! the invariants will be plotted twice ...
          if(sumpp+1.eq.lineends(jj)) then
             foundinv=foundinv+1
             sumpp=sumpp+1
!             write(*,*)'Assumed to be an invariant!',foundinv,kk
             write(21,549)xval(1,sumpp),yval(1,sumpp)
             write(21,549)xval(2,sumpp),yval(2,sumpp)
             write(21,549)zval(1,foundinv),zval(2,foundinv)
             write(21,549)xval(1,sumpp),yval(1,sumpp)
! we are at the end of a line, write a blank line
             write(21,548)jj
548          format('e '//'# end of invariant',2i5)
          else
             do while(sumpp.lt.lineends(jj))
                sumpp=sumpp+1
                write(21,549)xval(kk,sumpp),yval(kk,sumpp)
549             format(2f12.6,4i7)
! plotkod -101 means tieline
                if(plotkod(sumpp).eq.-101) write(21,552)
552             format(1x)
             enddo
! we are at the end of a line, write a blank line
             write(21,551)jj
551          format('e '//'# end of line',2i5)
          endif
       enddo
    enddo
!------------------------------------------------------------------------
! finally copy the data from the append file, it should be correctly formatted
    if(appfil.gt.0) then
       ic=0
1900   continue
       read(appfil,884,end=1910)appline
884    format(a)
       ic=ic+1
       if(appline(1:11).eq.'pause mouse') goto 1900
       write(21,884)trim(appline)
       goto 1900
1910   continue
!       write(*,*)'Appended ',ic,' data lines'
       close(appfil)
       appfil=0
    endif
!------------------------------------------------------
    if(pform(1:1).eq.' ') then
! if not hardcopy pause gnuplot.  Mouse means clicking in the graphics window
! will close it. I would like to have an option to spawn the graphics window...
! so it is kept while continuing the program.
       write(21,990)
990    format('pause mouse')
!990    format('e'//'pause mouse')
    endif
    close(21)
    if(appfil.ne.0) close(appfil)
    appfil=0

!    write(21,565)
!565 format('e'//'pause mouse'/)
!    close(21)
!
!    gnuplotline='gnuplot ocgnu.plt '
    gnuplotline='gnuplot '//pfc
! if gnuplot cannot be started with gnuplot give normal path ...
!    gnuplotline='"c:\program files\gnuplot\bin\wgnuplot.exe" '//pfc(1:kkk)//' '
!    k3=len_trim(gnuplotline)+1
    write(*,*)'Gnuplot command line: ',trim(gnuplotline)
    if(pform(1:1).ne.' ') then
       write(*,*)'Graphics output file: ',trim(pfh)
    endif
    if(btest(graphopt%status,GRKEEP)) then
       write(*,*)'Trying to spawn'
       call system(gnuplotline(9:))
    else
       call system(gnuplotline)
    endif
!900 continue
1000 continue
    return
  end subroutine ocplot3B

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine calc_diagram_point(axarr,pltax,xxx,xxy,line,ceq)
! calculates the equilibrium for axis coordinates xxx,xxy
! to obtain the set of stable phases
! axarr specifies calculation axis, 
! pltax plot axis
! xxx and xxy are axis coordinates for calculating a point
! line is a character where the stable phases at the point is returned
! ceq is the current equilibrium, should be the default with axis conditions
! ONLY COORDINATES FOR CALCULATION AXIS ALLOWED
    implicit none
    type(map_axis), dimension(*) :: axarr
    double precision xxx,xxy
    character line*(*),pltax(*)*(*)
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer jax,nax,kk,jj,ic,k3
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: svrrec
    double precision value
    character dummy*24
!
!    write(*,*)'Not implemented yet'
!    goto 1000
!
! We should check if plotaxis are the same as those used for calculation!!! 
! x-axis
    dummy=' '
    if(xxx.lt.axarr(1)%axmin .or. xxx.gt.axarr(1)%axmax) then
       write(*,11)'X value outside axis limits',xxx,&
            axarr(1)%axmin,axarr(1)%axmax
11     format(a,3(1pe14.5))
       gx%bmperr=4399; goto 1000
    endif
    call locate_condition(axarr(1)%seqz,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
! first argument 1 means to extract the value, 0 means to set the value
    call condition_value(0,pcond,xxx,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(pcond%active.ne.0) then
! active=0 means condition is active
       pcond%active=0
! we must indicate if T or P are now fixed ...
!       if(pcond%statev.eq.1) then
!          mapline%meqrec%tpindep(1)=.FALSE.
!       elseif(pcond%statev.eq.2) then
!          mapline%meqrec%tpindep(1)=.FALSE.
!       endif
    endif
! y-axis
    if(xxy.lt.axarr(2)%axmin .or. xxy.gt.axarr(2)%axmax) then
       write(*,11)'Y value outside axis limits',xxy,&
            axarr(2)%axmin,axarr(2)%axmax
       gx%bmperr=4399; goto 1000
    endif
    call locate_condition(axarr(2)%seqz,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
    call condition_value(0,pcond,xxy,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(pcond%active.ne.0) then
! active=0 means condition is active
       pcond%active=0
! we must indicate if T or P are now fixed ... ??
!       if(pcond%statev.eq.1) then
!          mapline%meqrec%tpindep(1)=.FALSE.
!       elseif(pcond%statev.eq.2) then
!          mapline%meqrec%tpindep(1)=.FALSE.
!       endif
    endif
! calculate the equilibrium without global minimization
!    call list_conditions(kou,ceq)
    call calceq3(0,.FALSE.,ceq)!
    if(gx%bmperr.ne.0) then
       write(*,*)'Calculation failed'
       goto 1000
    endif
! extract the names of the stable phases
    kk=1
    do jj=1,noph()
       do ic=1,noofcs(jj)
          k3=test_phase_status(jj,ic,value,ceq)
          if(k3.gt.0) then
! this phase is stable or fix
             call get_phase_name(jj,ic,dummy)
             line(kk:)=dummy
             kk=len_trim(line)+2
!             write(*,*)'Stable phases ',trim(line)
          endif
       enddo
    enddo
! replace _ with - in phase names
100 continue
    kk=index(line,'_')
    if(kk.gt.0) then
       line(kk:kk)='-'
       goto 100
    endif
1000 continue
    return
  end subroutine calc_diagram_point

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
    integer ntup,itup,iph,ics,nystat,inactive(4),notop,seqy,mode
    integer jj,seqz,iadd,irem,nv,saveq,lokcs
    type(gtp_phasetuple), dimension(:), allocatable :: entphcs
    integer, dimension(:), allocatable :: stsphcs
    type(map_line), pointer :: mapline
    type(map_fixph), pointer :: mapfix
!    TYPE(map_node), pointer :: curtop
    type(meq_setup), pointer :: meqrec
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: svr
    type(meq_phase), pointer :: phr
    double precision val,xxx,yyy,axvalok
    logical firstline
    integer, parameter :: maxsavedceq=2000
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
    do iph=1,noph()
       do ics=1,noofcs(iph)
          itup=itup+1
!          write(*,*)'SMP ',iph,noofcs(iph),ics,itup
!          entphcs(itup)%phaseix=iph
          entphcs(itup)%ixphase=iph
          entphcs(itup)%compset=ics
          stsphcs(itup)=test_phase_status(iph,ics,val,ceq)
          if(gx%bmperr.ne.0) goto 1000
          val=zero
          call change_phase_status(iph,ics,nystat,val,ceq)
          if(gx%bmperr.ne.0) goto 1000
       enddo
    enddo
!    write(*,*)'Suspended all phases'
! indicator if maptop allocated
!    nullify(curtop)
    notop=0
! loop through all phases with stsphcs less than 3
!    nystat=0
!============================================================
    phaseloop: do itup=1,ntup
! nystat:-4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fix
       if(stsphcs(itup).gt.-2) then
! set default constitution, if none specified in the middle
!          iph=entphcs(itup)%phaseix
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
             call get_phasetup_name_old(entphcs(itup),name)
! axis variable is composition, skip hases with no variance
!             call get_phase_variance(entphcs(itup)%phaseix,nv)
             call get_phase_variance(entphcs(itup)%ixphase,nv)
             if(nv.eq.0) then
                write(*,71)name(1:len_trim(name)),val
71              format(/'Ignoring phase with fixed composition: ',a,F10.6)
!----------------
                lokcs=phasetuple(iph)%lokvares
!                write(*,*)'indices: ',iph,phasetuple(iph)%phaseix,lokcs
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
          call map_findline(maptop,axarr,mapfix,mapline)
          if(gx%bmperr.ne.0) goto 500
!          write(*,*)'Back from map_findline 2'
!          write(*,*)'We have a line'
          ceq=>mapline%lineceq
          meqrec=>mapline%meqrec
! this is the first equilibrium along the line, create meqrec in step_separate
305       continue
          call calceq7(mode,meqrec,mapfix,ceq)
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
             write(*,*)'Generating mapline%meqrec failed: ',gx%bmperr
             gx%bmperr=0; goto 333
          else
! calculation OK, do it again (why?) without creating meqrec, save and
! return here after taking a step using the same meqrec
380          continue
             iadd=0
             irem=0
             call meq_sameset(irem,iadd,mapline%meqrec,&
                  mapline%meqrec%phr,inmap,ceq)
             if(gx%bmperr.ne.0) then
!                write(*,*)'Error calling meq_sameset',gx%bmperr
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
             call map_step(maptop,mapline,mapline%meqrec,mapline%meqrec%phr,&
                  axvalok,noofaxis,axarr,ceq)
             if(gx%bmperr.ne.0) then
! if error just terminate line
                write(*,*)'Error return from map_step: ',gx%bmperr
                gx%bmperr=0; goto 333
             endif
             if(mapline%more.ge.0) goto 380
          endif
333       continue
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
!    write(*,*)'Restoring previous phase status'
    do itup=1,ntup
!       write(*,910)itup,entphcs(itup)%phase,entphcs(itup)%compset,stsphcs(itup)
910    format('Restoring all phase status: ',4i5)
!       call change_phase_status(entphcs(itup)%phaseix,&
       call change_phase_status(entphcs(itup)%ixphase,&
            entphcs(itup)%compset,stsphcs(itup),val,ceq)
       if(gx%bmperr.ne.0) goto 1000
    enddo
1000 continue
    return
  end subroutine step_separate

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  logical function abbr_phname_same(full,short)
! return TRUE if short is a correct abbreviation of full
! This is used in macro step4 to plot fractions in different composition sets
    implicit none
    character*(*) full,short
!\end{verbatim}
    logical same
    integer k1,k2
    character*1 ch1,ch2
!    write(*,*)'Comparing ',trim(full)//' : '//trim(short)
    same=.false.
! unequal if full has no # or different index after
    k1=index(short,'#')
    if(k1.gt.0) then
       k2=index(full,'#')
! full has no compset
       if(k2.le.0) then
! if short has #1 then the full phase without # should be accepted
!          write(*,*)'full has no compset:  ',short(k1+1:k1+1),k2
          if(short(k1+1:k1+1).eq.'1') then
             same=.true.
             goto 1000
          endif
       else
! the character after # must be the same
          if(short(k1+1:k1+1).eq.full(k2+1:k2+1)) then
             same=.true.
             goto 1000
          endif
       endif
    endif
! if short is without # then all compsets match
    if(compare_abbrev(short,full)) then
       same=.true.
    endif
1000 continue
!    write(*,*)'Comparing ',trim(short)//' with '//trim(full),' is ',same
    abbr_phname_same=same
    return
  end function abbr_phname_same

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine get_plot_conditions(text,ndx,axarr,ceq)
! extacts the conditions from ceq and removes those that are axis variables
    implicit none
    character text*(*)
    integer ndx
    type(map_axis), dimension(ndx) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer jj,seqz,ip,jp
    character symbol*24
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: svrrec,svr2
!
! get all conditions
    call get_all_conditions(text,-1,ceq)
!    write(*,*)'PC1: ',trim(text),ndx
    do jj=1,ndx
! replace the values of those that are axis with X or Y
       seqz=axarr(jj)%seqz
       call locate_condition(seqz,pcond,ceq)
       if(gx%bmperr.ne.0) goto 1000
       svrrec=>pcond%statvar(1)
       symbol=' '
       ip=1
       call encode_state_variable(symbol,ip,svrrec,ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'PC2: ',symbol(1:ip-1),jj
       jp=index(text,symbol(1:ip-1))
       if(jp.gt.0) then
          seqz=jp+index(text(jp:),'=')-1
          ip=jp+index(text(jp:),' ')-1
          if(jj.eq.1) then
             text(seqz:)='=X, '//text(ip:)
          else
             text(seqz:)='=Y, '//text(ip:)
          endif
!       else
!          write(*,*)'Cannot find: ',symbol(1:ip-1),' in ',trim(text)
       endif
    enddo
1000 continue
    return
  end subroutine get_plot_conditions

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

END MODULE ocsmp

! map/step subroutines:
! map_setup should be OK
! map_startpoint to convert from start equilibria to point on a line
! map_replaceaxis ok
! map_startline ?? stupid name > map_sputility ??
! map_store should be OK
! map_lineend should be OK
! map_changeaxis change of active axis, new axis specified
! map_force_changeaxis works only with 2 axis and isopleths(?)
! map_step OK
! map_calcnode should be OK
! map_newnode cannot handle isopleths with 3 exits or invariants
! map_reserve_saveceq
! map_findline should be ok
! create_saveceq
! delete_mapresults
! integer function tieline_inplane
! logical function inveq (invariant_equilibrium)
! map_problems called when problems
! map_halfstep when convergence trouble
! list_stored_equilibria
! amend_stored_equilibria
! ocplot3 extract data for diagrams with two wildcard axis
! calc_diagram_point calculates an equilibrium for specified coordinates
! ocplot3B generate GNUPLOT file with data from ocplot3
! ocplot2 extracting data for all diagrams except isothermal sections
! ocplot2B generate GNUPLOT file with data from ocplot2
! step_separate calculates each phase separately
! abbr_phname_same checks if a phase name is an abbreviation including compset

! There should be a reorganizing of mapnodes at the end
!

