! included in smp2.F90
!\begin{verbatim}
  subroutine ocplot2(ndx,maptop,axarr,graphopt,version,ceq)
!  subroutine ocplot2(ndx,pltax,filename,maptop,axarr,graphopt,pform,&
!       version,ceq)
! Main plotting routine, generates a GNUPLOT data file for a step/map calc
! NOTE for isothermal section ocplot3 is used (when 2 axis with wildcards)
! ndx is mumber of plot axis, 
! pltax is text with plotaxis variables
! filename is the name of the GNUPLOT file
! maptop is map_node record with all results
! axarr is array of axis records
! graphopt is graphical option record
! NOT USED: pform is type of output (screen/acrobat/postscript/gif)
! ceq is equilibrium record
    implicit none
    integer ndx
!    character pltax(*)*(*),filename*(*),pform*(*)
    character version*(*)
    type(map_axis), dimension(*) :: axarr
    type(map_node), pointer :: maptop
    type(graphics_options) :: graphopt
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
! local variables set from graphopt
!    character pltax(2)*64,filename*64,pform*32
    character pltax(2)*64,filename*64
!
    type(map_ceqresults), pointer :: results
    TYPE(gtp_equilibrium_data), pointer :: curceq
    type(map_node), pointer :: mapnode,invar,localtop
    type(map_line), pointer :: mapline
    logical wildcard
!    integer, parameter :: mofapl=100
    character ch1*1,gnuplotline*256,pfd*64,pfc*64
!    character pfh*64,dummy*24,applines(mofapl)*128,appline*128
    character pfh*64,dummy*24
    double precision, dimension(:,:), allocatable :: anp
    double precision, dimension(:), allocatable :: xax,yyy
! Too big??
!    integer, parameter :: maxval=10000
    integer, parameter :: maxval=2000
    integer, dimension(:), allocatable :: nonzero,linzero,linesep
!    integer, dimension(:), allocatable :: linesep
! encoded2 stores returned text from get_many ... 2048 is too short ...
! selphase used when plotting data just for a selected phase like y(fcc,*)
    character statevar*64,encoded1*1024,encoded2*4096,selphase*24
    character*64, dimension(:), allocatable :: phaseline
    integer i,ic,jj,k3,kk,kkk,lokcs,nnp,np,nrv,nv,nzp,ip,nstep,nnv,nofapl
    integer nr,line,next,seqx,nlinesep,ksep,iax,anpax,notanp,appfil
    double precision xmax,xmin,ymax,ymin,value,anpmin,anpmax
! lhpos is last used position in lineheader
    integer giveup,nax,ikol,maxanp,lcolor,lhpos,repeat,anpdim,qp
    character date*8,mdate*12,title*128,backslash*2,lineheader*1024
    character deftitle*128,labelkey*64
    logical overflow,first,last,novalues,selectph,varofun,moretops
    logical, allocatable, dimension(:) ::  nevernone
! textlabels
    type(graphics_textlabel), pointer :: textlabel
! line identification (title)
    character*16, dimension(:), allocatable :: lid
!
!    write(*,*)'In ocplot2, looking for segmentation fault 1'
! transfer from graphics record to local variables
    pltax(1)=graphopt%pltax(1)
    pltax(2)=graphopt%pltax(2)
    filename=graphopt%filename
!    write(*,*)' ocplot2 >> plot file: ',trim(filename)
!    pform=graphopt%pform
! continue as before ...
    if(index(pltax(1),'*').gt.0 .and. index(pltax(2),'*').gt.0) then
!       write(*,*)'Using ocplot3'
       call ocplot3(ndx,pltax,filename,maptop,axarr,graphopt,&
            version,ceq)
       goto 1000
    endif
    moretops=.FALSE.
    seqx=0
    call date_and_time(date)
    mdate=" "//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//" "
    deftitle='OpenCalphad '//version//': '//mdate//': with GNUPLOT'
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
!    write(*,*)'SMP: Allocating linesep: ',jj
    allocate(linesep(jj))
    nax=maptop%number_ofaxis
    linesep=0
! allocate texts to identify the lines on the gnuplot file
!    write(*,*)'SMP: Allocating phaseline: ',jj
    allocate(phaseline(jj))
!    if(maptop%number_ofaxis.gt.1) then
!       write(*,*)'Warning: may not not handle map graphics correctly',jj
!    endif
!
    giveup=0
    nrv=maxval
!    write(*,*)'SMP: allocating xax: ',nrv
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
    selphase=' '
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
!          write(*,*)'SMP: allocating anp1: ',np*nrv
          allocate(anp(np,nrv))
!          write(*,*)'SMP: allocating anp2: ',np
! nonzero indicates for each column if there is any nonzero value
! columns with only zero values will be eliminated before plotting
          allocate(nonzero(np))
!          write(*,*)'SMP: allocating nonzero: ',np
! linzero indicate for the present block of equilibria for each column
! if this column contain nonzero values
          allocate(linzero(np))
!          write(*,*)'SMP: allocating yyy: ',np
          allocate(yyy(np))
          nzp=np
! nzp should be dimension of yyy, np returns the number of values in yyy
! yyy is to extract state variable values for the column with wildcard
! NOTE binary phase diagrams are plotted with wildcard axis like x(*,cr) vs T
! nevernone is an attempt to remove columns that are zero by the value NaN
          allocate(nevernone(np))
          nevernone=.FALSE.
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
                   write(*,*)'SMP2B wildcard selected phase: ',trim(selphase)
                   selectph=.TRUE.
                endif
             endif
          else
             write(*,*)'Wildcard without ,!'
          endif
       endif
    enddo
    if(.not.wildcard) then
       anpdim=1
!       write(*,*)'SMP: allocating anp2: ',np
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
!    write(*,*)'SMP at label 77A: ',localtop%lines
    ikol=0
    do nrv=1,localtop%lines
       localtop%linehead(nrv)%done=0
    enddo
! we sometimes have a segmentation fault when several maptops ...
    if(associated(localtop%next)) then
       mapnode=>localtop%next
    else
       write(*,*)'Mapnode next link missing 1'
       goto 79
    endif
!    write(*,*)'SMP at label 77B: ',localtop%lines
    thisloop: do while(.not.associated(mapnode,localtop))
       do nrv=1,mapnode%lines
          mapnode%linehead(nrv)%done=0
       enddo
       if(.not.associated(mapnode%next)) then
          write(*,*)'Mapnode next link missing 2'
          exit thisloop
       endif
       mapnode=>mapnode%next
    enddo thisloop
!-----------
79  continue
    results=>localtop%saveceq
    mapnode=>localtop
    line=1
! looking for segmentation fault running map11.OCM as part of all.OCM in oc4P
! This error may be due to having created (or not created) composition sets ...
!    write(*,*)'ocplot2 after label 79'
! extract the names of stable phases for this lone

!    write(*,*)'mapnode index: ',mapnode%seqx
!    write(*,*)'Before label 100: ',results%free
!------------------------------------------- begin loop 100
100    continue
       mapline=>localtop%linehead(line)
! initiate novalues to false for each line
       novalues=.false.
! We have a segmentation fault sfter this in oc4P when running map11.OCM
! at the end of running all macros.  
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
! segmentation fault after here
!             write(*,*)'In ocplot2, looking for segmentation fault 4A'
             first=.false.
             kk=1
!             if(selectph) novalues=.TRUE.
             phloop: do jj=1,noph()
!                write(*,*)'In ocplot2, segmentation fault 4Ax',jj,noofcs(jj)
                do ic=1,noofcs(jj)
                   k3=test_phase_status(jj,ic,value,curceq)
                   if(gx%bmperr.ne.0) goto 1000
                   if(k3.gt.0) then
! this phase is stable or fix
                      call get_phase_name(jj,ic,dummy)
                      if(gx%bmperr.ne.0) goto 1000
                      if(selectph) then
! this is an attempt to remove lines from irrelevant equilibria when plotting
! data for a specific phase like y(fcc#4,*) which is not stable??
                         if(abbr_phname_same(dummy,trim(selphase))) then
                            novalues=.FALSE.
!                            write(*,*)'SMP2B novalues FALSE',mapline%lineid
                            exit phloop
                         else
                            novalues=.TRUE.
!                            write(*,*)'SMP2B novalues TRUE',mapline%lineid
                         endif
                      endif
                      if(.not.novalues) then
! I think this phaseline is no longer used ?? YES it is
!                         write(*,*)'SMP addto phaseline 1: ',trim(dummy),&
!                              nlinesep
                         phaseline(nlinesep)(kk:)=dummy
                         kk=len_trim(phaseline(nlinesep))+2
                      endif
                   endif
                enddo
             enddo phloop
!             write(*,117)nlinesep,phaseline(nlinesep)(1:kk-1)
117          format('Stable phases on line ',i5/a)
! the segmentation fault was that linzero not always allocated ....
             if(allocated(linzero)) linzero=0
          endif
! no wildcards allowed on this axis
!          write(*,*)'In ocplot2, segmentation fault 4Ay1',nr,nv
          statevar=pltax(notanp)
          call meq_get_state_varorfun_value(statevar,value,encoded1,curceq)
!          write(*,*)'SMP axis variable 1: ',trim(encoded1),value
          if(gx%bmperr.ne.0) then
! this error should not prevent plotting the other points FIRST SKIPPING
             write(*,212)'SMP skipping a point, error evaluating: ',&
                  statevar(1:10),curceq%tpval(1),nv,nr
212          format(a,a,f10.2,2i5)
! buperr resets putfun error 
             gx%bmperr=0; buperr=0
             nv=nv-1; goto 215
          endif
          xax(nv)=value
!          write(*,201)'at 202: ',nr,nv,curceq%tpval(1),value
!          xax(nv)=curceq%tpval(1)
! macro step1 run in parallel plotting cp has segm fault after this line
!          write(*,*)'After label 200: ',mapline%lineid,nr,nv
          if(xax(nv).lt.xmin) xmin=xax(nv)
          if(xax(nv).gt.xmax) xmax=xax(nv)
! second axis
          statevar=pltax(anpax)
!          varofun=.FALSE.
!          write(*,*)'In ocplot2, segmentation fault after 4B ',wildcard
          if(wildcard) then
! NEW ignore data for this equilibrium if NOVALUES is TRUE
! because "selphase" not equal to the stable phase found above
             if(novalues) then
!                write(*,*)'Ignoring equilibria without ',trim(selphase)
                yyy=zero
!                np=1
! skip this equilibrium, nv=nv-1, and take next equilibrium, increement nr!!
                nv=nv-1
                goto 199
!                cycle plot1
             else
!                write(*,*)'SMP2B wildcard value 1: ',nr,trim(statevar)
!                write(*,*)'In ocplot2, segmentation fault after 4C1: ',&
!                     trim(statevar),nooftup()
! segmentation fault is inside this call for map11.OCM
! probably because new composition set created
                call get_many_svar(statevar,yyy,nzp,np,encoded2,curceq)
!                write(*,*)'In ocplot2, segmentation fault search: '
! compiling without -finit-local-zero gives a segmentation fault here
! running the MAP11 macro
                qp=np
!                write(*,*)'SMP2B wildcard value 2: ',nr,trim(statevar)
!                write(*,223)'SMP2B Values: ',np,(yyy(i),i=1,np)
!                if(selectph) then
!                   write(*,*)'SMP2B: number of values: ',trim(selphase),np,nv
!223                format(a,i3,8F8.4)
!                endif
                if(gx%bmperr.ne.0) then
                   write(*,*)'yaxis error: "',trim(statevar),'"'
                   goto 1000
                endif
             endif
!             write(*,*)'On ocplot2, segmentation fault 4D1'
!             write(*,213)trim(encoded2),np,(yyy(ic),ic=1,np)
213          format('windcard: ',a,i3,6(1pe12.4))
!             write(*,16)'val: ',kp,nr,gx%bmperr,(yyy(i),i=1,np)
16           format(a,2i3,i5,/6(1pe11.3/))
             anpmin=1.0D20
             anpmax=-1.0D20
             lcolor=0
!             write(*,*)'On ocplot2, segmentation fault before 4D2',np
! this is a loop for all values for this equilibria
! Here we may try to replace zero values by RNONE ???
!             write(*,*)'SMP2B RNONE: ',RNONE
             do jj=1,np
                if(last) then
                   if(linzero(jj).ne.0) then
! in last equilibria we may have a value from the new phase at the node
                      anp(jj,nv)=yyy(jj)
!                   elseif(yyy(jj).ne.zero) then
!                      write(*,*)'SMP skipping a value for line ',nlinesep
                   endif
                else
! trying to avoid plotting a line at zero for unstable/unused state variables 
                   if(yyy(jj).eq.zero) then
!                      if(nax.eq.1) then
! FOR STEP calculations try to make the ending of a property at zero ...
! This is for STEP 1 figure 3
! BUT it did not work, many strange lines appeared in other plots ...
!                         if(nv.gt.1 .and. abs(anp(jj,nv-1)).ne.zero) then
! Hm, we cannot set anp(jj,nv) to zero because then all will be zero ...
!                            anp(jj,nv-1)=zero
!                         else
!                            anp(jj,nv)=rnone
!                         endif
!                      else
! for MAP calculations just ignore the point ... also for STEP ...
                         anp(jj,nv)=rnone
!                      endif
                   else
! Hm, jumps from zero to finite values in step1, fig 3 plotting w(phase,cr) ..
!                      if(nv.gt.1 .and. anp(jj,nv-1).eq.rnone) then
!                         anp(jj,nv-1)=zero
!                      endif
                      anp(jj,nv)=yyy(jj)
                   endif
                   if(abs(yyy(jj)).gt.zero) then
                      nonzero(jj)=1
                      linzero(jj)=1
! save the first column with nonzero for use with invariants
                      if(ikol.eq.0) ikol=jj
                      if(anp(jj,nv).gt.anpmax) anpmax=anp(jj,nv)
!                      if(anp(jj,nv).lt.anpmin) anpmin=anp(jj,nv)
                      if(anp(jj,nv).ne.rnone .and. &
                           anp(jj,nv).lt.anpmin) anpmin=anp(jj,nv)
! extract state variable jj
                      if(.not.allocated(lid)) then
!                         write(*,*)'Allocating lid: ',np+5
                         allocate(lid(np+5))
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
! A single state variable or function value like CP
! I HAVE HAD PROBLEMS WITH NEGATIVE CP HERE 
! try skipping this value (below) if last equilibrium on the line 
!             varofun=.TRUE.
!             write(*,*)'SMP: calling meq_get_state_varofun ',trim(statevar)
! there is a segmentation fault in this call
             call meq_get_state_varorfun_value(statevar,value,encoded1,curceq)
!             write(*,*)'SMP axis variable 2: ',statevar(1:3),value
             if(gx%bmperr.ne.0) then
! SECOND Skipping
                write(*,212)'SMP Skipping a point, error evaluating: ',&
                     statevar(1:10),curceq%tpval(1),nv,nr
                nv=nv-1; goto 215
             endif
             if(results%savedceq(nr)%nexteq.eq.0) then
! THIRD ?? Skipping
!                write(*,212)'SMP skip last evaluated symbol: ',&
!                     trim(statevar),curceq%tpval(1),nv,nr
                if(trim(statevar).ne.trim(encoded1)) then
! If "statevar" not equal to "encoded1" skip last point
! This is a clumsy way to avoid negative CP=H.T values at end of lines ...
                   nv=nv-1; goto 215
                endif
             endif
!             if(gx%bmperr.ne.0) goto 1000
             anp(1,nv)=value
! macro test step1 run in parallel has segme fault plotting cp before this line
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
199       continue
          nr=curceq%nexteq
          if(nr.gt.0) then
             if(results%savedceq(nr)%nexteq.eq.0) then
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
!       write(*,*)'In ocplot2, segmentation fault 4F'
       if(nax.gt.1) then
!---------------------------------------------------------------
!------------------ special for invariant lines
!---------------------------------------------------------------
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
!                            write(*,*)'SMP addto phaseline 2: ',trim(dummy),&
!                                 nlinesep
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
!                   write(*,335)'New line: ',nlinesep,nv,linesep(nlinesep),&
!                        statevar(1:5),value
335                format(a,3i4,' <',a,'> ',3(1pe14.6))
! axis with possible wildcard
                   statevar=pltax(anpax)
!                   write(*,*)'In ocplot2, segmentation fault 4H'
                   if(wildcard) then
! this cannot be a state variable derivative
!                    write(*,*)'Getting a wildcard value 2: ',nr,statevar(1:20)
                      call get_many_svar(statevar,yyy,nzp,np,encoded2,curceq)
                      if(gx%bmperr.ne.0) goto 1000
! save one non-zero value per line, 3 lines
                      ic=0
                      do jj=1,np
! we have put anp to zero above ??
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
! for RNONE = NaN add empty line after invariant ...
                      if(linesep(nlinesep).lt.nv) then
! we should never have several linesep for the same value of nv!
                         nlinesep=nlinesep+1
                         linesep(nlinesep)=nv
!                         write(*,*)'Empty line after invariant: ',nlinesep,nv
                      endif
                   else
! if no wildcard there is no invariant line
!                      write(*,*)'It can be an invariant point here!!'
                      nv=nv-3
                      goto 225
                   endif
                endif inv
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
             if(nnp.gt.maxanp) then
                write(kou,*)'Too many points in anp array 1',maxanp,nv,nnp
                overflow=.TRUE.
                nnp=maxanp
             endif
             if(nv.gt.maxval) then
                write(kou,*)'Too many points in anp array 2',maxval,nv
                overflow=.TRUE.
                nv=maxval
             endif
             if(.not.allocated(lid)) then
!                write(*,*)'SMP allocating lid 3: ',np
                allocate(lid(np+5))
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
             write(*,*)'plot data overflow',nv,nnp
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
!       write(*,*)'SMP: allocate lid 4: ',np
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
! NOW pltax should be the the axis labels if set manually
    if(graphopt%labeldefaults(2).ne.0) pltax(1)=graphopt%plotlabels(2)
    if(graphopt%labeldefaults(3).ne.0) pltax(2)=graphopt%plotlabels(3)
!    write(*,*)' >>>>>>>>**>>> plot file: ',trim(filename)
    call ocplot2B(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,anp,lid,&
         phaseline,title,filename,graphopt,version,encoded1)
!         title,filename,graphopt,pform,version,encoded1)
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
  subroutine ocplot2B_old(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,&
       anp,lid,phaseline,title,filename,graphopt,version,conditions)
!       title,filename,graphopt,pform,version,conditions
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
! NOT USED: pform is output form (scree/acrobat/postscript/git)
! conditions is a character with the conditions for the diagram
    implicit none
    integer np,anpax,nlinesep
    integer ndx,nrv,linesep(*),anpdim
!    character pltax(*)*(*),filename*(*),pform*(*),lid(*)*(*),title*(*)
    character pltax(*)*(*),filename*(*),lid(*)*(*),title*(*)
    character conditions*(*),version*(*),phaseline(*)*(*)
    type(graphics_options) :: graphopt
    double precision xax(*),anp(anpdim,*)
    double precision scale1,scalem
    type(graphics_textlabel), pointer :: textlabel
!\end{verbatim}
!----------------------------------------------------------------------
! internal
    integer ii,jj,kk,lcolor,appfil,nnv,ic,repeat,ksep,nv,k3,kkk,nofapl
    integer, parameter :: mofapl=100
    character pfc*64,pfh*64,backslash*2,appline*128
    character applines(mofapl)*128,gnuplotline*80,labelkey*64,rotate*16
    character labelfont*16,linespoints*12
! write the gnuplot command file with data appended
!
!    write(*,10)'in ocplot2B: ',np,anpax,nrv,pform(1:1),trim(title),&
!         nlinesep,(linesep(kk),kk=1,nlinesep)
10  format(a,3i5,' "',a,'" '/a/i3,2x,15i4)
!    write(*,*)'In ocplot2B filename: ',trim(filename)
    if(index(filename,'.plt ').le.0) then 
       kk=len_trim(filename)
       pfc=filename(1:kk)//'.'//'plt '
       kkk=kk+4
    else
       pfc=filename
       kk=index(pfc,'.')-1
       kkk=len_trim(filename)+1
    endif
!    write(*,*)'In OCPLOT2B opening ',trim(pfc)
    open(21,file=pfc,access='sequential',status='unknown')
!    write(*,*)'In OCPLOT2B after opening ',trim(pfc)
    write(21,1600)trim(title)
1600 format('# GNUPLOT file generated by OpenCalphad'/'# ',a)
! if there is just one curve do not write any key.  May be overriiden later ..
    if(graphopt%gnutermsel.lt.1 .or. &
         graphopt%gnutermsel.gt.graphopt%gnutermax) then
       write(*,*)'Unknown graphics terminal: ',graphopt%gnutermsel
       goto 1000
    elseif(graphopt%gnutermsel.gt.1) then
! terminal 1 is screen without any output file
       pfh=filename(1:kk)//'.'//graphopt%filext(graphopt%gnutermsel)
! set the screen as a comment ...
       write(21,840)trim(graphopt%gnuterminal(1)),&
            trim(graphopt%gnuterminal(graphopt%gnutermsel)),trim(pfh)
840    format('#set terminal ',a/'set terminal ',a/'set output "',a,'"')
    else
! terminal 1 is screen without any output file
       write(21,841)trim(graphopt%gnuterminal(graphopt%gnutermsel))
841    format('set terminal ',a)
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
! OC logo oclogo added by Catalina Pineda
!    write(*,*)'Plot heading 1? ',btest(graphopt%status,GRNOTITLE)
    if(btest(graphopt%status,GRNOTITLE)) then
       write(21,858)trim(title),trim(conditions)
    else
       write(21,859)trim(title),trim(conditions)
    endif
    write(21,860)graphopt%xsize,graphopt%ysize,&
         trim(pltax(1)),trim(pltax(2)),trim(labelkey)
858 format('# set title "',a,' \n #',a,'" font "arial,10" ')
859 format('set title "',a,' \n ',a,'" font "arial,10" ')
860 format('set origin 0.0, 0.0 '/&
         'set size ',F8.4', ',F8.4/&
         'set xlabel "',a,'"'/'set ylabel "',a,'"'/&
         'set label "O" at graph -0.090, -0.100 font "Garamond bold,20"'/&
         'set label "C" at graph -0.080, -0.100 font "Garamond bold,20"'/&
         'set key ',a/&
         'set style line 1 lt 2 lc rgb "#000000" lw 2 pt 10'/&
         'set style line 2 lt 2 lc rgb "#4169E1" lw 2 pt 6'/&
         'set style line 3 lt 2 lc rgb "#00C000" lw 2 pt 3'/&
         'set style line 4 lt 2 lc rgb "#FF0000" lw 2 pt 2'/&
         'set style line 5 lt 2 lc rgb "#0080FF" lw 2 pt 4'/&
         'set style line 6 lt 2 lc rgb "#C8C800" lw 2 pt 5'/&
         'set style line 7 lt 2 lc rgb "#C0C0C0" lw 2 pt 7'/&
         'set style line 8 lt 2 lc rgb "#00FFFF" lw 2 pt 8'/&
         'set style line 9 lt 2 lc rgb "#804080" lw 2 pt 9'/&
         'set style line 10 lt 2 lc rgb "#7CFF40" lw 2 pt 1')
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
! set labels (user added text in diagram)
    textlabel=>graphopt%firsttextlabel
    do while(associated(textlabel))
       rotate=' '
       if(textlabel%angle.ne.0) write(rotate,177)textlabel%angle
177    format(' rotate by ',i5)
       labelfont=' '
!       write(*,*)'textfontscale: ',textlabel%textfontscale
       if(textlabel%textfontscale.ne.one) then
          write(labelfont,178)int(10*textlabel%textfontscale)
178       format(' font "Sans,',i2,'" ')
       endif
!       if(textlabel%angle.eq.0) then
       write(21,1505)trim(textlabel%textline),textlabel%xpos,textlabel%ypos,&
            trim(labelfont),trim(rotate)
1505   format('set label "',a,'" at ',1pe12.4,', ',1pe12.4,a,a)
!       else
!         write(21,1506)trim(textlabel%textline),textlabel%xpos,textlabel%ypos,&
!               textlabel%angle
!1506      format('set label "',a,'" at ',1pe12.4,', ',1pe12.4,&
!               ' rotate by ',i5)
!       endif
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

! note if append file is GIBBSTRIANGLE
       if(appline(1:10).eq.'# GIBBSTRI') then
          write(*,*)'Warning: appended file is in Gibbstriangle format,',&
               ' plot will be strange!'
          goto 1710
       endif
!------------------------------------------------------------------
! ignore some lines with "set" in the append file
! set title
! set xlabel
! set ylabel
! set xrange
! set yrange
! set terminal
! set origin
! set size
! set key
       if(appline(1:10).eq.'set title ' .or.&
            appline(1:11).eq.'set xlabel ' .or.&
            appline(1:11).eq.'set ylabel ' .or.&
            appline(1:11).eq.'set xrange ' .or.&
            appline(1:11).eq.'set yrange ' .or.&
            appline(1:11).eq.'set output ' .or.&
            appline(1:13).eq.'set terminal ' .or.&
            appline(1:11).eq.'set origin ' .or.&
            appline(1: 9).eq.'set size ' .or.&
            appline(1: 8).eq.'set key ') then
!          write(*,*)'ignoring append line ',trim(appline)
          goto 1710
       endif
!------------------------------------------------------------------
       if(index(appline,'plot "-"').gt.0 .or. &
            index(appline,"plot '-'").gt.0 .or.&
            index(appline,"plot for ").gt.0) then
! here we save the actual plot commands from the appendfile!!
          applines(1)=appline
          ic=1
1730      continue
! if line ends with \ then read more
          ii=len_trim(appline)
!          write(*,*)'There are more? ',appline(i:ii),ii,ic
          if(appline(ii:ii).eq.'\') then
! continuation lines !! NOTE EACH plot expected at the beginning of the line
             read(appfil,1720,end=1750)appline
             ic=ic+1
             if(ic.ge.mofapl) then
                write(*,*)'Too many header lines in append file',ic
             else
                applines(ic)=appline
             endif
             goto 1730
          endif
! debug output of saved plot command
          nofapl=ic
!          write(*,*)(trim(applines(jj)),jj=1,nofapl)
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
!-----------------------------------------------
! text in lower left corner
    ii=len_trim(graphopt%lowerleftcorner)
    if(ii.gt.0) then
! in square diagram below figure
       write(21,209)trim(graphopt%lowerleftcorner)
209    format('set label "',a,'" at graph -0.10, -0.08 ')
    endif
! om lowerleftcorner is empty ignore it
!---------------------------------------------------------------
! this is subroutine ocplot2B_old
!---------------------------------------------------------------
    backslash=',\'
! here we should start from the value in graphopt%linett
!    lcolor=1
    lcolor=graphopt%linett
    if(lcolor.lt.1 .or. lcolor.gt.10) then
       write(*,*)'Illegal line color',lcolor
       lcolor=1
    endif
! if graphopt%linestyle=0 use lines, otherwise linespoints
    if(graphopt%linestyle.eq.0) then
       linespoints='lines'
    else
       linespoints='linespoints'
    endif
!    write(*,*)'ocplot2B linett: ',lcolor
!    write(*,*)'backslash "',backslash,'" '
    if(anpax.eq.2) then
! anpax=2 means the single valued axis is colum 2 and possibly multiple
! values in column 3 and higher
! this is the multiple plot, file name only given once!!
! last line tuple on a separate format statement, if np>2
! np is number of columns
       if(np.eq.1 .and. appfil.eq.0) then
          write(21,880)trim(linespoints),lcolor,' ',' '
       else
!          write(21,880)lcolor,trim(lid(1)),backslash
!880       format('plot "-" using 2:3 with lines ls ',i2,' title "'a,'"',a)
          write(21,880)trim(linespoints),lcolor,trim(lid(1)),backslash
880       format('plot "-" using 2:3 with ',a,' ls ',i2,' title "'a,'"',a)
       endif
       do ii=2,np-1
          lcolor=lcolor+1
          if(lcolor.gt.10) lcolor=1
!          write(21,882)ii+2,lcolor,trim(lid(ii)),backslash
!882       format('"" using 2:',i3,' with lines ls ',i2,' title "'a,'"',a)
          write(21,882)ii+2,trim(linespoints),lcolor,trim(lid(ii)),backslash
882       format('"" using 2:',i3,' with ',a,' ls ',i2,' title "'a,'"',a)
       enddo
       lcolor=lcolor+1
       if(lcolor.gt.10) lcolor=1
       if(appfil.eq.0) then
!          if(np.ge.2) write(21,882)np+2,lcolor,trim(lid(np)),' '
       if(np.ge.2) write(21,882)np+2,trim(linespoints),lcolor,trim(lid(np)),' '
       else
! write the last calculated curve if np>1
!          if(np.ge.2) write(21,882)ii+2,lcolor,trim(lid(ii)),backslash
          if(np.ge.2) write(21,882)ii+2,trim(linespoints),&
               lcolor,trim(lid(ii)),backslash
! we should append data, change plot "-" to just "" in appline(1)
          ii=index(applines(1),'plot "-"')
          applines(1)(1:ii+7)='""'
          write(*,*)'Inserting the plot commands from append file',nofapl
          do ii=1,nofapl
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
       lcolor=graphopt%linett+1
       if(lcolor.lt.2 .or. lcolor.gt.10) then
          write(*,*)'Illegal line color',lcolor
          lcolor=2
       endif
!       lcolor=2
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
!          write(*,*)'Adding the plot commands from append file',nofapl,ic
          do ii=1,nofapl
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
! ANPAX is axis with multiple values
    if(anpax.ne.0) then
       scalem=graphopt%scalefact(anpax)
       scale1=graphopt%scalefact(3-anpax)
    else
       scalem=graphopt%scalefact(1)
       scale1=graphopt%scalefact(2)
    endif
!
    do nv=1,nrv
!---------------------------------------------------------------
! trying to handle RNONE
! values written multiplied with graphopt%scalefact, 
! first value is single valued axis (can be X or Y axis) multiplied with scale1
! remaining values multiplied svalem
       write(21,2820,advance='no')nv,scale1*xax(nv)
       do jj=1,np-1
! second and later columns represent Y axis
          if(anp(jj,nv).ne.rnone) then
             write(21,2821,advance='no')scalem*anp(jj,nv)
          else
             write(21,2822,advance='no')
          endif
       enddo
       if(anp(jj,nv).ne.rnone) then
          write(21,2821)scalem*anp(jj,nv)
       else
          write(21,2822)
       endif
2820   format(i4,1pe16.6)
2821   format(1pe16.6)
2822   format(' NaN ')
!---------------------------------------------------------------
       if(nv.eq.linesep(ksep)) then
! an empty line in the dat file means a MOVE to the next point.
          if(nv.lt.nrv) then
             if(jj.gt.1) then
! phaseline is never used in the plot, just included in the file
                write(21,1819)ksep,trim(phaseline(ksep))
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
!    write(*,*)'In OCPLOT2B finished 1 "',pform(1:1),'"'
!    if(pform(1:1).eq.' ') then
    if(graphopt%gnutermsel.eq.1) then
! if not hardcopy pause gnuplot.  Mouse means clicking in the graphics window
! will close it. I would like to have an option to keep the graphics window...
       write(21,990)trim(graphopt%plotend)
!990    format('pause mouse')
990    format(a)
    endif
    close(21)
!    write(*,*)'In OCPLOT2B closed ',trim(pfc),kkk
    if(appfil.ne.0) close(appfil)
    appfil=0
!-------------------------------------------------------------------
! execute the GNUPLOT command file
    gnuplotline='gnuplot '//pfc(1:kkk)//' & '
! Uncomment the following line for OS having gnuplot5 and 
! comment the above line with gnuplotline='gnuplot '//pfc(1:kkk)//' & '
!  gnuplotline='gnuplot5 '//pfc(1:kkk)//' & '
! Reason - Choose the gnuplot command as per the Operating System's 
! existing command "gnuplot", "gnuplot5"etc..
! Or, give the path of the gnuplot file as described below:
! if gnuplot cannot be started with gnuplot give normal path ...
!    gnuplotline='"c:\program files\gnuplot\bin\wgnuplot.exe" '//pfc(1:kkk)//' '
    k3=len_trim(gnuplotline)+1
    write(kou,*)'Gnuplot command file: ',pfc(1:kk+4)
!    if(pform(1:1).ne.' ') then
    if(graphopt%gnutermsel.ne.1) then
       write(kou,*)'Graphics output file: ',pfh(1:kk+4)
    endif
! GRWIN set by compiler option, if 1 we are running on Microsoft windows
    if(grwin.eq.1) then
! call system without initial "gnuplot " keeps the window !!!
       if(btest(graphopt%status,GRKEEP)) then
!          write(*,*)'plot command: "',gnuplotline(9:k3),'"'
!          write(*,*)'Trying to spawn: ',trim(gnuplotline)
!          call system(gnuplotline(9:k3))
! spawn plot on Windows ?? NOT ISO-TERMAL DIAGRAM
!          write(*,*)'executing command: "start /B '//trim(gnuplotline)
!          call execute_command_line('start /B '//gnuplotline(9:k3))
          write(*,*)'executing command: "start /B '//trim(gnuplotline)//'"'
          call execute_command_line('start /B '//trim(gnuplotline))
! WORKS WITH OCPLOT3B
!          call execute_command_line('start /B '//trim(gnuplotline))
       else
!          write(*,*)'plot command: "',gnuplotline(1:k3),'"'
!          call system(gnuplotline)
          write(*,*)'executing command: "'//trim(gnuplotline)//'"'
          call execute_command_line(gnuplotline)
       endif
    else
! plot on non-windows system, do not use "start /B"
! how to implement GRKEEP?
       write(*,*)'executing command: '//trim(gnuplotline)
       call execute_command_line(gnuplotline)
    endif
1000 continue
    return
  end subroutine ocplot2B_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
  subroutine ocplot2B(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,anp,lid,&
       phaseline,title,filename,graphopt,version,conditions)
!       title,filename,graphopt,pform,version,conditions)
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
! NOT USED: pform is output form (scree/acrobat/postscript/git)
! conditions is a character with the conditions for the diagram
    implicit none
    integer np,anpax,nlinesep
    integer ndx,nrv,linesep(*),anpdim
    character pltax(*)*(*),filename*(*),lid(*)*(*),title*(*)
    character conditions*(*),version*(*),phaseline(*)*(*)
    type(graphics_options) :: graphopt
    double precision xax(*),anp(anpdim,*)
    double precision scale1,scalem
    type(graphics_textlabel), pointer :: textlabel
!\end{verbatim}
!----------------------------------------------------------------------
! internal
    integer ii,jj,kk,lcolor,appfil,nnv,ic,repeat,ksep,nv,k3,kkk,nofapl
    integer, parameter :: mofapl=100
    integer appfiletyp
    character pfc*64,pfh*64,backslash*2,appline*128
    character applines(mofapl)*128,gnuplotline*80,labelkey*64,rotate*16
    character labelfont*16,linespoints*12,tablename*16,year*16,hour*16
! write the gnuplot command file with data appended
!
!    write(*,10)'in ocplot2B: ',np,anpax,nrv,pform(1:1),trim(title),&
!         nlinesep,(linesep(kk),kk=1,nlinesep)
10  format(a,3i5,' "',a,'" '/a/i3,2x,15i4)
!    write(*,*)'In ocplot2B filename: ',trim(filename)
    if(index(filename,'.plt ').le.0) then 
       kk=len_trim(filename)
       pfc=filename(1:kk)//'.'//'plt '
       kkk=kk+4
    else
       pfc=filename
       kk=index(pfc,'.')-1
       kkk=len_trim(filename)+1
    endif
!    write(*,*)'In OCPLOT2B opening ',trim(pfc)
    open(21,file=pfc,access='sequential',status='unknown')
!    write(*,*)'In OCPLOT2B after opening ',trim(pfc)
    write(21,1600)trim(title)
1600 format('# GNUPLOT file generated by OpenCalphad'/'# ',a)
! if there is just one curve do not write any key.  May be overriiden later ..
    if(graphopt%gnutermsel.lt.1 .or. &
         graphopt%gnutermsel.gt.graphopt%gnutermax) then
       write(*,*)'Unknown graphics terminal: ',graphopt%gnutermsel
       goto 1000
    elseif(graphopt%gnutermsel.gt.1) then
! terminal 1 is screen without any output file
       pfh=filename(1:kk)//'.'//graphopt%filext(graphopt%gnutermsel)
! set the screen as a comment ...
       write(21,840)trim(graphopt%gnuterminal(1)),&
            trim(graphopt%gnuterminal(graphopt%gnutermsel)),trim(pfh)
840    format('#set terminal ',a/'set terminal ',a/'set output "',a,'"')
    else
! terminal 1 is screen without any output file
       write(21,841)trim(graphopt%gnuterminal(graphopt%gnutermsel))
841    format('set terminal ',a)
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
    call date_and_time(year,hour)
!    write(*,*)'"',year,'"  "',hour,'"'
    tablename='OCT'//year(3:8)//hour(1:6)
! OC logo oclogo added by Catalina Pineda
!    write(*,*)'Plot heading 2? ',btest(graphopt%status,GRNOTITLE)
    if(btest(graphopt%status,GRNOTITLE)) then
       write(21,858)trim(title),trim(conditions)
    else
       write(21,859)trim(title),trim(conditions)
    endif
    write(21,860)graphopt%xsize,graphopt%ysize,&
         trim(pltax(1)),trim(pltax(2)),trim(labelkey)
858 format('#set title "',a,' \n #',a,'" font "arial,10" ')
859 format('set title "',a,' \n ',a,'" font "arial,10" ')
860 format('set origin 0.0, 0.0 '/&
         'set size ',F8.4', ',F8.4/&
         'set xlabel "',a,'"'/'set ylabel "',a,'"'/&
         'set label "O" at graph -0.090, -0.100 font "Garamond bold,20"'/&
         'set label "C" at graph -0.080, -0.100 font "Garamond bold,20"'/&
         'set key ',a/&
         'set style line 1 lt 2 lc rgb "#000000" lw 2 pt 10'/&
         'set style line 2 lt 2 lc rgb "#4169E1" lw 2 pt 6'/&
         'set style line 3 lt 2 lc rgb "#00C000" lw 2 pt 3'/&
         'set style line 4 lt 2 lc rgb "#FF0000" lw 2 pt 2'/&
         'set style line 5 lt 2 lc rgb "#0080FF" lw 2 pt 4'/&
         'set style line 6 lt 2 lc rgb "#C8C800" lw 2 pt 5'/&
         'set style line 7 lt 2 lc rgb "#C0C0C0" lw 2 pt 7'/&
         'set style line 8 lt 2 lc rgb "#00FFFF" lw 2 pt 8'/&
         'set style line 9 lt 2 lc rgb "#804080" lw 2 pt 9'/&
         'set style line 10 lt 2 lc rgb "#7CFF40" lw 2 pt 1')
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
! set labels (user added text in diagram)
    textlabel=>graphopt%firsttextlabel
    do while(associated(textlabel))
       rotate=' '
       if(textlabel%angle.ne.0) write(rotate,177)textlabel%angle
177    format(' rotate by ',i5)
       labelfont=' '
!       write(*,*)'textfontscale: ',textlabel%textfontscale
       if(textlabel%textfontscale.ne.one) then
          write(labelfont,178)int(10*textlabel%textfontscale)
178       format(' font "Sans,',i2,'" ')
       endif
!       if(textlabel%angle.eq.0) then
       write(21,1505)trim(textlabel%textline),textlabel%xpos,textlabel%ypos,&
            trim(labelfont),trim(rotate)
1505   format('set label "',a,'" at ',1pe12.4,', ',1pe12.4,a,a)
!       else
!         write(21,1506)trim(textlabel%textline),textlabel%xpos,textlabel%ypos,&
!               textlabel%angle
!1506      format('set label "',a,'" at ',1pe12.4,', ',1pe12.4,&
!               ' rotate by ',i5)
!       endif
       textlabel=>textlabel%nexttextlabel
    enddo
!---------------------------------------------------------------
! handle appended files here ....
!
    appfil1: if(graphopt%appendfile(1:1).eq.' ') then
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
! note if append file is GIBBSTRIANGLE
       if(appline(1:10).eq.'# GIBBSTRI') then
          write(*,*)'Warning: appended file is in Gibbstriangle format,',&
               ' plot will be strange!'
          goto 1710
       endif
!------------------------------------------------------------------
! ignore some lines with "set" in the append file
! set title
! set xlabel
! set ylabel
! set xrange
! set yrange
! set terminal
! set origin
! set size
! set key
       if(appline(1:10).eq.'set title ' .or.&
            appline(1:11).eq.'set xlabel ' .or.&
            appline(1:11).eq.'set ylabel ' .or.&
            appline(1:11).eq.'set xrange ' .or.&
            appline(1:11).eq.'set yrange ' .or.&
            appline(1:11).eq.'set output ' .or.&
            appline(1:13).eq.'set terminal ' .or.&
            appline(1:11).eq.'set origin ' .or.&
            appline(1: 9).eq.'set size ' .or.&
            appline(1: 8).eq.'set key ') then
!          write(*,*)'ignoring append line ',trim(appline)
          goto 1710
       endif
!------------------------------------------------------------------
! changes here in the new subroutine ocplot2B
!       write(*,*)'SMP: appline1: ',trim(appline)
       ii=index(appline,'plot "-"')
       if(ii.gt.0) then
          applines(1)=appline
          appfiletyp=1
       else
          ii=index(appline,"plot '-'")
          if(ii.gt.0) then
             applines(1)=appline
             appfiletyp=1
          else
! we can also have a "plot for ... " do not change applines(1)
             ii=index(appline,"plot for ")
             if(ii.le.0) then
! just copy the file to ocgnu.plt
                write(21,'(a)')trim(appline)
                goto 1710
             else
! this is a "plot for" appfile using a table that has already been copied
                applines(1)=appline
                appfiletyp=2
             endif
          endif
       endif
! we have now found the plot command in the append file. There can be more
!       write(*,*)'SMP appfiletyp: ',appfiletyp
! here we save the actual plot commands from the appendfile!!
       applines(1)=appline
!       write(*,*)'SMP appline1: ',trim(applines(1))
       ic=1
1730   continue
! if line ends with \ then read more
       ii=len_trim(appline)
! write(*,*)'There are more? ',appline(i:ii),ii,ic
!       if(appline(ii:ii).eq.'\') then
       if(appline(ii:ii).eq.' ') then
! continuation lines !! NOTE EACH plot expected at the beginning of the line
          read(appfil,1720,end=1750)appline
          ic=ic+1
          if(ic.ge.mofapl) then
             write(*,*)'Too many header lines in append file',ic
          else
             applines(ic)=appline
          endif
          goto 1730
       endif
! debug output of saved plot command
       nofapl=ic
!          write(*,*)(trim(applines(jj)),jj=1,nofapl)
!          write(*,*)'appline: ',trim(appline),ic
!          close(appfil)
!          appfil=0
       goto 1770
!       endif
!       write(*,*)'SMP appline1B: ',trim(appline)
       write(21,1720)trim(appline)
       nnv=nnv+1
       goto 1710
! error oppening append file
1750   continue
       write(kou,*)'Error opening or reading the append file, skipping it'
       close(appfil)
       appfil=0
1770   continue
    endif appfil1
! end of appendfile special
!-----------------------------------------------
! text in lower left corner
    ii=len_trim(graphopt%lowerleftcorner)
    if(ii.gt.0) then
! in square diagram below figure
       write(21,209)trim(graphopt%lowerleftcorner)
209    format('set label "',a,'" at graph -0.10, -0.08 ')
    endif
! if lowerleftcorner is empty ignore it
!---------------------------------------------------------------
! this is new subroutine ocplot2B
!---------------------------------------------------------------
!========================================= begin using plot for
! ANPAX is axis with multiple values
    if(anpax.ne.0) then
       scalem=graphopt%scalefact(anpax)
       scale1=graphopt%scalefact(3-anpax)
    else
       scalem=graphopt%scalefact(1)
       scale1=graphopt%scalefact(2)
    endif
! now write all data once as a table ended with EOD
    write(21,3000)trim(tablename)
3000 format('$',a,' << EOD')
!
! columnheaders used as keys
    write(21,3100)'KEYS: ',trim(pltax(3-anpax)),(trim(lid(jj)),jj=1,np)
3100 format(a,a,' ',100(a,' '))
    ksep=2
    do nv=1,nrv
!---------------------------------------------------------------
! values written multiplied with graphopt%scalefact, 
! first value is single valued axis (can be X or Y axis) multiplied with scale1
! remaining values multiplied svalem
       write(21,'(i4,1pe16.6)',advance='no')nv,scale1*xax(nv)
       do jj=1,np-1
! second and later columns represent Y axis
          if(anp(jj,nv).ne.rnone) then
             write(21,2821,advance='no')scalem*anp(jj,nv)
          else
             write(21,'(a)',advance='no')' NaN '
          endif
       enddo
       if(anp(jj,nv).ne.rnone) then
          write(21,2821)scalem*anp(jj,nv)
       else
          write(21,'(a)')' NaN '
       endif
!2820   format(i4,1pe16.6)
2821   format(1pe16.6)
!2822   format(' NaN ')
!---------------------------------------------------------------
       if(nv.eq.linesep(ksep)) then
! an empty line in the dat file means a MOVE to the next point.
          if(nv.lt.nrv) then
             write(21,3819)ksep-1
3819         format('# end of line ',i2//)
          else
! try to avoid rubbish
             write(21,3821)
3821         format('# end of line '//)
          endif
! test of uninitiallized variable, ksep must not exceed nlinesep
          ksep=min(ksep+1,nlinesep)
       endif
    enddo
    write(21,'(a)')'EOD'
    if(appfil.gt.0) then
! if there is an appendfile add set multiplot
! The "writeback" is important for uniform scaling of multiplots
       write(21,3828)
3828   format('set multiplot'/&
            'set xrange [] writeback'/'set yrange [] writeback')
    endif
! if anpax is axis with single value (1=x, 2=y)
    if(anpax.eq.1) then
       write(21,3900)np+2,trim(tablename)
3900   format('plot for [i=3:',i2,'] $',a,' using i:2',&
            ' with lines ls (i-2) title columnheader(i)') 
    else
       write(21,3910)np+2,trim(tablename)
3910   format('plot for [i=3:',i2,'] $',a,' using 2:i',&
            ' with lines ls (i-2) title columnheader(i)') 
    endif
! plot command from appfil
    if(appfil.gt.0) then
! try to avoid overlapping keys ...
! The "restore" for x/yrange means the scaling from the "plot for"
! will be used also for the appended data
       write(21,3912)
3912   format('set key bottom right font "arial,12"'/&
            'set xrange restore'/'set yrange restore')
       if(appfiletyp.eq.2) then
! just one line with plot for ... 
! the data to append is already copied as a table
          write(21,'(a)')trim(applines(1))
          close(appfil)
          appfil=0
          write(21,'(a)')'unset multiplot'
       else
          do jj=1,nofapl
             write(21,'(a)')trim(applines(jj))
          enddo
       endif
    endif
! if the plot command is "plot '-' ... then
! copy the data from the append file, it should be correctly formatted
    if(appfil.gt.0) then
       ic=0
1900   continue
! this is copying the actual data to plot from the append file.
       read(appfil,884,end=1910)appline
884    format(a)
       ic=ic+1
       if(appline(1:11).eq.'pause mouse') goto 1900
       write(21,884)trim(appline)
       goto 1900
1910   continue
!       write(*,*)'Appended ',ic,' data lines'
       write(21,'(a)')'unset multiplot'
       close(appfil)
       appfil=0
    endif
!------------------------------------------------------
!    write(*,*)'In OCPLOT2B finished 1 "',pform(1:1),'"'
    if(graphopt%gnutermsel.eq.1) then
! if not hardcopy pause gnuplot.  Mouse means clicking in the graphics window
! will close it. I would like to have an option to keep the graphics window...
       write(21,990)trim(graphopt%plotend)
!990    format('pause mouse')
990    format(a)
    endif
    close(21)
!    write(*,*)'In OCPLOT2B closed ',trim(pfc),kkk
    if(appfil.ne.0) close(appfil)
    appfil=0
!-------------------------------------------------------------------
! execute the GNUPLOT command file
    gnuplotline='gnuplot '//pfc(1:kkk)//' & '
! Uncomment the following line for OS having gnuplot5 and 
! comment the above line with gnuplotline='gnuplot '//pfc(1:kkk)//' & '
!  gnuplotline='gnuplot5 '//pfc(1:kkk)//' & '
! Reason - Choose the gnuplot command as per the Operating System's 
! existing command "gnuplot", "gnuplot5"etc..
! Or, give the path of the gnuplot file as described below:
! if gnuplot cannot be started with gnuplot give normal path ...
!    gnuplotline='"c:\program files\gnuplot\bin\wgnuplot.exe" '//pfc(1:kkk)//' '
    k3=len_trim(gnuplotline)+1
    write(kou,*)'Gnuplot command file: ',pfc(1:kk+4)
    if(graphopt%gnutermsel.ne.1) then
       write(kou,*)'Graphics output file: ',pfh(1:kk+4)
    endif
! GRWIN is set by compiler option, if 1 we are running microspft windows
    if(grwin.eq.1) then
! call system without initial "gnuplot " keeps the window !!!
       if(btest(graphopt%status,GRKEEP)) then
!          write(*,*)'plot command: "',gnuplotline(9:k3),'"'
!          write(*,*)'Trying to spawn: ',trim(gnuplotline)
!          call system(gnuplotline(9:k3))
! spawn plot on Windows ?? NOT ISO-TERMAL DIAGRAM
!          write(*,*)'executing command: "start /B '//trim(gnuplotline)
!          call execute_command_line('start /B '//gnuplotline(9:k3))
          write(*,*)'executing command: "start /B '//trim(gnuplotline)//'"'
          call execute_command_line('start /B '//trim(gnuplotline))
! WORKS WITH OCPLOT3B
!          call execute_command_line('start /B '//trim(gnuplotline))
       else
!          write(*,*)'plot command: "',gnuplotline(1:k3),'"'
!          call system(gnuplotline)
          write(*,*)'executing command: "'//trim(gnuplotline)//'"'
          call execute_command_line(gnuplotline)
       endif
    else
! plot on non-windows system without "start /B ...
! how to implement GRKEEP?
       write(*,*)'executing command: '//trim(gnuplotline)
       call execute_command_line(gnuplotline)
    endif
1000 continue
    return
  end subroutine ocplot2B

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
!  subroutine ocplot3(ndx,pltax,filename,maptop,axarr,graphopt,pform,&
  subroutine ocplot3(ndx,pltax,filename,mastertop,axarr,graphopt,&
       version,ceq)
! special to plot isothermal sections (two columns like x(*,cr) x(*,ni))
! ndx is mumber of plot axis, 
! pltax is text with plotaxis variables
! filename is intermediary file (maybe not needed)
! mastertop is the map_node record with all results
! axarr are axis records
! graphopt is graphics record (should be extended to include more)
! pform is graphics form
! NOT USED: pform is type of output (screen or postscript or gif)
    implicit none
    integer ndx
    character pltax(*)*(*),filename*(*),version*(*)
    type(map_axis), dimension(*) :: axarr
    type(map_node), pointer :: mastertop
    type(graphics_options) :: graphopt
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
    integer, parameter :: maxsame=200
    type(map_node), pointer :: plottop,curtop,endnode,maptop
    type(map_line), pointer :: curline
    type(gtp_equilibrium_data), pointer :: curceq
    type(map_ceqresults), pointer :: results
    integer ii,jj,point,plotp,lines,eqindex,lasteq,nooftup,lokcs,jph,same,kk
    integer, parameter :: maxval=2000,mazval=100
    double precision, allocatable :: xval(:,:),yval(:,:),zval(:,:),tieline(:,:)
    integer offset,nofeq,sumpp,last,nofinv,ntieline,mtieline,noftielineblocks
    double precision xxx,yyy
! TO BE REPLACED BY GNUTERMSEL: pform
!    character pform*8
    integer, allocatable :: plotkod(:),lineends(:)
    character xax1*8,xax2*24,yax1*8,yax2*24,axis1*32,axisx(2)*32,axisy(2)*32
    character phname*32,encoded*1024,axis*32
    character lid(2,maxsame)*24
    integer nooflineends
!
! do not change mastertop!
    maptop=>mastertop
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
    nooflineends=0
!    write(*,17)
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
    noftielineblocks=0
!    write(*,*)'SMP Number of lines: ',lines
    node: do ii=1,lines
! up to version 5.014: (june 2018)
!       curline=>curtop%linehead(ii)
! because crash plotting a ternary with 2 start points I changed
!       curline=>plottop%linehead(ii)
! BUT the line above does not work for map10 and map3 (last plot of H)
! so until further investigations I keep:
       curline=>curtop%linehead(ii)
       if(btest(curline%status,EXCLUDEDLINE)) then
          write(*,*)'Excluded line: ',ii,curline%lineid,lines
          cycle node
       endif
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
! estimate the number of tie-lines to extract
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
          if(.not.associated(curceq)) then
             write(*,*)'SMP error, equilibrium missing?: ',&
                  curtop%seqx,curline%lineid,eqindex
             cycle line
          endif
! find the stable phases (max 3)
          plotp=plotp+1
          if(plotp.gt.maxval) then
             write(*,*)'Too many points to plot ',maxval
             cycle node
          endif
          jj=1
          equil: do jph=1,nooftup
             lokcs=phasetuple(jph)%lokvares
             call get_phasetup_name(jph,phname)
! crash as lokcs not valid ... the 3rd of 4th time plotted ...
!             write(*,321)'SMP bug? ',curtop%seqx,curline%lineid,eqindex,&
!                  jph,lokcs,trim(phname)
321          format(a,2i4,i5,2i4,' : ',a)
! crash next line in alcrni-1200 mapping ...
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
!
!>>>>>>> here we should allow a wildcard axis like ac(*)
!>>>>>>> without any phase label!!                
!                
                if(axisx(1)(1:1).eq.' ') then
                   axisx(1)=trim(xax1)//trim(phname)//trim(xax2)
                   axisy(1)=trim(yax1)//trim(phname)//trim(yax2)
                elseif(axisx(2)(1:1).eq.' ') then
                   axisx(2)=trim(xax1)//trim(phname)//trim(xax2)
                   axisy(2)=trim(yax1)//trim(phname)//trim(yax2)
                endif
!                write(*,*)'We are here 1X: ',trim(axisx(jj))
                call meq_get_state_varorfun_value(axisx(jj),xxx,encoded,curceq)
                xval(jj,plotp)=xxx
!                write(*,*)'We are here 1Y: ',trim(axisy(jj))
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
! there is a node at the end, extracts its ceq record
          curceq=>endnode%nodeceq
          plotp=plotp+1
          if(plotp.gt.maxval) then
             write(*,*)'Too many points to plot 1:',maxval
             gx%bmperr=4399; goto 1000
          endif
          do jj=1,2
!             write(*,*)'We are here 2:',trim(axisx(jj))
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
       noftielineblocks=noftielineblocks+1
       if(ntieline.gt.0) then
! All tielines on the same line with a space in between
          same=same+1
          do eqindex=1,ntieline
             plotp=plotp+1
             if(plotp.gt.maxval) then
                write(*,*)'Too many points to plot 2:',maxval
                gx%bmperr=4399; goto 1000
             endif
             xval(1,plotp)=tieline(1,eqindex)
             yval(1,plotp)=tieline(2,eqindex)
! this means the tie-lines will be plotted twice ... but why not??
             xval(2,plotp)=tieline(1,eqindex)
             yval(2,plotp)=tieline(2,eqindex)
             plotkod(plotp)=-100
             plotp=plotp+1
             if(plotp.gt.maxval) then
                write(*,*)'Too many points to plot 3:',maxval
                gx%bmperr=4399; goto 1000
             endif
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
! but first generate the monovariant (not invariant)
    curtop=>curtop%next
    if(.not.associated(curtop,maptop)) then
!       write(*,*)'Extracting data from node'
       curceq=>curtop%nodeceq
       if(associated(curceq)) then
!          write(*,*)'Extracting data from node equilibrium'
          plotp=plotp+1
          if(plotp.gt.maxval) then
             write(*,*)'Too many points to plot 4:',maxval
             gx%bmperr=4399; goto 1000
          endif
          jj=1
          same=same+1
          nodeequil: do jph=1,nooftup
             lokcs=phasetuple(jph)%lokvares
             if(curceq%phase_varres(lokcs)%phstate.ge.PHENTSTAB) then
! plotkod set to -1 to indicate monovariant (not invariant)
                plotkod(plotp)=-1
                call get_phasetup_name(jph,phname)
!                write(*,*)'Stable phase ',trim(phname),jj
!                if(jj.lt.3) lid(jj,same)='invariant'
                if(jj.lt.3) lid(jj,same)='monovariant'
                axis=trim(xax1)//trim(phname)//trim(xax2)
!                write(*,*)'We are here 3: ',trim(axis)
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
          lineends(same)=plotp
       endif
!       write(*,*)'jump back to 100, same and plotp',same,plotp
       goto 100
    endif
!    do jj=1,same
!       write(*,23)'phases: ',same,jj,trim(lid(1,jj)),trim(lid(2,jj))
!    enddo
!------------------------------------------------
! there can be more maptops linked via plotlink
    if(associated(plottop%plotlink)) then
       jj=plottop%seqx
       plottop=>plottop%plotlink
       write(*,*)'ocplot3B current and next maptop: ',jj,plottop%seqx
! this added 180918 to plot results from several MAP commands
       curtop=>plottop
       maptop=>plottop
!       write(*,*)'Current number of lines: ',same,plotp,&
!            allocated(curtop%linehead)
       goto 100
    endif
!========================================================
    call get_plot_conditions(encoded,maptop%number_ofaxis,axarr,ceq)
! now we should have all data to plot in xval and yval arrays
500 continue
    write(*,808)same,plotp,maxsame,maxval
808 format('plot data used: ',2i7,' out of ',2i7)
!    write(*,*)'found lines/points to plot: ',same,plotp,nofinv
!    write(*,502)(lineends(ii),ii=1,same)
502 format(10i5)
! NOW pltax should be the the axis labels if set manually
    if(graphopt%labeldefaults(2).ne.0) pltax(1)=graphopt%plotlabels(2)
    if(graphopt%labeldefaults(3).ne.0) pltax(2)=graphopt%plotlabels(3)
    call ocplot3B(same,nofinv,lineends,2,xval,2,yval,2,zval,plotkod,pltax,&
         lid,filename,graphopt,version,encoded)
    deallocate(xval)
    deallocate(yval)
    deallocate(plotkod)
1000 continue
    return
  end subroutine ocplot3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
  subroutine ocplot3B(same,nofinv,lineends,nx1,xval,ny1,yval,nz1,zval,plotkod,&
       pltax,lid,filename,graphopt,version,conditions)
! called by ocplot3 to write the GNUPLOT file for two wildcard columns
! same is the number of lines to plot
! nofinv number of monovariants (not invariants)
! lineends array with row numbers where each line ends
! nx1 first dimension of xval
! xval 2D matrix with values to plot on x axis
! ny1 first dimension of yval
! yval 2D matrix with values to plot on y axis
! nz1 first dimension of zval
! zval 2D matrix with third point of monovariant (not invariant) triangles
! plotkod integer array indicating the type of line (-1 skip line)
! pltax text for axis
! lid array with GNUPLOT line types
! filename is intermediary file (maybe not needed)
! graphopt is graphics option record
! maptop is map_node record with all results
! REPLACED BY gnutermsel pform is type of output (screen/acrobat/postscript/gir)
! conditions is a text with conditions for the calculation
    implicit none
!    character pltax(*)*(*),filename*(*),pform*(*),lid(nx1,*)*(*),conditions*(*)
    character pltax(*)*(*),filename*(*),lid(nx1,*)*(*),conditions*(*)
    character version*(*)
    type(graphics_options) :: graphopt
    integer same,plotkod(*),nx1,ny1,nz1,nofinv
    integer lineends(*)
    double precision xval(nx1,*),yval(ny1,*),zval(nz1,*)
!\end{verbatim}
    integer, parameter :: maxcolor=200,mofapl=100
    integer ii,jj,kk,jph,offset,n1,nofapl
    type(graphics_textlabel), pointer :: textlabel
    character gnuplotline*64,date*12,mdate*12,title*128,deftitle*64,backslash*2
    character labelkey*64,applines(mofapl)*128,appline*128,pfc*80,pfh*80
    integer sumpp,np,appfil,ic,nnv,kkk,lcolor(maxcolor),iz,again
    integer done(maxcolor),foundinv,fcolor,k3
    character color(maxcolor)*24,rotate*16,labelfont*16,linespoints*12
    integer naptitle,apptitles(maxcolor)
! we plot monovariant twice, once with border once with filledcurves!!
!    integer noofmono,jjj,monovariant2(100)
    integer, parameter :: monovariantborder=11
    integer xtieline,xmonovariant
! now a global variable
!    character monovariant*6
! Gibbs triangle variables
    logical plotgt,appgt
    double precision sqrt3,xxx,yyy,xmax,ltic
!  
    write(*,*)'Using the rudimentary graphics in ocplot3B!',graphopt%linett
! light green            'fc "#B0FFB0" notitle ',a)
! very faint green       'fc "#F0FFF0" notitle ',a)
! yellow                 'fc "#FFFF00" notitle ',a)
! goldenrod              'fc "#DAA520" notitle ',a)
! dark green    monovariant='50FF50'
!    monovariant='00FFFF'
!    write(*,*)'Filename: ',trim(ocgnu)
!
    call date_and_time(date)
    mdate=" "//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//" "
    deftitle='OpenCalphad '//version//': '//mdate//': with GNUPLOT'
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
! problems with identifying monovariants (not invariants) and tie-lines lines
    lcolor=0
!
    if(index(filename,'.plt ').le.0) then 
       kk=len_trim(filename)
       pfc=filename(1:kk)//'.'//'plt '
       kkk=kk+4
    else
       pfc=filename
    endif
!    write(kou,*)'Filename: ',trim(pfc),', format: ',pform(1:1)
!    write(kou,*)'Graphics format: ',graphopt%gnutermsel
    open(21,file=pfc,access='sequential',status='unknown')
    write(21,110)trim(title)
110 format('# GNUPLOT file generated by OpenCalphad'/'# ',a)
    if(graphopt%gnutermsel.lt.1 .or. &
         graphopt%gnutermsel.gt.graphopt%gnutermax) then
       write(*,*)'Unknown graphics terminal: ',graphopt%gnutermsel
       goto 1000
    elseif(graphopt%gnutermsel.gt.1) then
! terminal 1 is screen without any output file
       pfh=filename(1:kk)//'.'//graphopt%filext(graphopt%gnutermsel)
! set the screen as a comment ...
       write(21,840)trim(graphopt%gnuterminal(1)),&
            trim(graphopt%gnuterminal(graphopt%gnutermsel)),trim(pfh)
840    format('#set terminal ',a/'set terminal ',a/'set output "',a,'"')
! 840    format('set terminal ',a/'set output "',a,'"')
    else
! terminal 1 is screen without any output file
       write(21,841)trim(graphopt%gnuterminal(graphopt%gnutermsel))
841    format('set terminal ',a)
    endif
    if(graphopt%gibbstriangle) then
!       write(*,*)'Gibbs triangle diagrams are under development!'
! xmax should be the largest scale value of x and y axis (same length!!)
       xmax=one
       plotgt=.true.
       if(graphopt%rangedefaults(1).ne.0) then
          xmax=min(xmax,graphopt%plotmax(1))
       endif
       if(graphopt%rangedefaults(2).ne.0) then
          xmax=min(xmax,graphopt%plotmax(2))
       endif
       ltic=0.01*xmax
       sqrt3=0.5D0*sqrt(3.0D0)
       write(21,844)sqrt3*xmax,xmax
!       write(21,844)sqrt3*xmax,xmax,sqrt3*xmax+0.02
! These maybe not necessary ... 0.866 is 0.5sqrt(3)
844    format('# GIBBSTRIANGLE '/&
            'set bmargin 3'/'set lmargin 3'/'set rmargin 3'/'set tmargin 3'/&
            'set origin 0.0, 0.0 '/&
            'set size ratio 0.866'/&
            'set yrange [0:',F10.6,']'/'set xrange [0:',F10.6,']'/&
            'set noborder'/'set noxtics'/'set noytics')
!            "set label 'Z' at 0, -0.03 center"/&
!            "set label 'X' at 1, -0.03 center"/&
!            "set label 'Y' at 0.5,",F10.6," center")
! This replaces axis without tics, only a tic in the middle
       write(21,845)xmax, xmax, 0.5*xmax, sqrt3*xmax, 0.5*xmax, sqrt3*xmax,&
! next 3 values are value and position for max values of Y axis
!            xmax, 0.343*xmax, sqrt3*(xmax+0.15d0), &
            xmax, 0.343*one, sqrt3*(one+0.15d0), &
! next 3 values are value and positions of max values for X axis
!            xmax, 0.92*xmax, -5.0*ltic, &
            xmax, 0.92*one, -5.0*ltic, &
! these are rudimentary ticmarks
            0.25*xmax-2*ltic, 0.5*sqrt3*xmax, 0.25*xmax, 0.5*sqrt3*xmax, &
            0.25*xmax-ltic,0.5*sqrt3*xmax+1.5*ltic,&
            0.25*xmax,0.5*sqrt3*xmax,0.75*xmax+2*ltic,&
            0.5*sqrt3*xmax,0.75*xmax,0.5*sqrt3*xmax,&
            0.75*xmax+ltic,0.5*sqrt3*xmax+1.5*ltic,0.75*xmax,&
            0.5*sqrt3*xmax,0.5*xmax,-ltic,0.5*xmax,0.0
845    format('set style line 90 lt 1 lw 3 pt -1 ps 1'/&
            'set style line 91 lt 1 lw 2 pt -1 ps 1'/&
            'set arrow 1 from 0,0 to ',F8.4,', 0.0 nohead linestyle 90'/&
            'set arrow 2 from ',F8.4,', 0 to ',F8.4,',',F8.4,&
            ' nohead linestyle 90'/&
            'set arrow 3 from ',F8.4,',',F8.4,' to 0,0 nohead linestyle 90'/&
            '# axis max values ...'/&
            'set label "',F6.2,'" at graph ',F8.4,',',F8.4/& 
            'set label "',F6.2,'" at graph ',F8.4,',',F8.4/& 
            '# tickmarks ...'/&
            'set arrow 4 from ',F8.4,',',F8.4,' to ',F8.4,',',F8.4,&
            ' nohead linestyle 91'/&
            'set arrow 5 from ',F8.4,',',F8.4,' to ',F8.4,',',F8.4,&
            ' nohead linestyle 91'/&
            'set arrow 6 from ',F8.4,',',F8.4,' to ',F8.4,',',F8.4,&
            ' nohead linestyle 91'/&
            'set arrow 7 from ',F8.4,',',F8.4,' to ',F8.4,',',F8.4,&
            ' nohead linestyle 91'/&
            'set arrow 8 from ',F8.4,',',F8.4,' to ',F8.4,',',F8.4,&
            ' nohead linestyle 91'/&
            '# end most of special Gibbs triangle')
    else
       plotgt=.false.
    endif
!
    write(*,*)'Plot heading 3? ',btest(graphopt%status,GRNOTITLE)
    if(btest(graphopt%status,GRNOTITLE)) then
       write(21,128)trim(title),trim(conditions)
    else
       write(21,129)trim(title),trim(conditions)
    endif
    write(21,130)graphopt%xsize,graphopt%ysize,&
         trim(pltax(1)),trim(labelkey)
128 format('#set title "',a,' \n #',a,'" font "arial,10"')
129 format('set title "',a,' \n ',a,'" font "arial,10"')
130 format('set origin 0.0, 0.0 '/&
         'set size ',F8.4', ',F8.4/&
         'set xlabel "',a,'"'/&
         'set key ',a)
    if(plotgt) then
! OC logo added by Catalina Pineda
! when Gibbs triangle the ylabel and logo must be placed carefully
! THIS IS THE Y-AXIS WITH 60 degrees angle
       write(21,131)trim(pltax(2)), 0.15*xmax, 0.37*xmax
131    format('set label "',a,'" at ',F8.4,',',F8.4,' rotate by 60 '/&
!            'set label "O" at screen 0.130, 0.027 font "Garamond bold,20"'/&
!            'set label "C" at screen 0.139, 0.027 font "Garamond bold,20"')
            'set label "O" at graph -0.090, -0.100 font "Garamond bold,20"'/&
            'set label "C" at graph -0.077, -0.100 font "Garamond bold,20"')
! we should also enforce same length of X and Y axis !!!
    else
! SQUARE DIAGRAM
       write(21,132)trim(pltax(2))
132    format('set ylabel "',a,'"'/&
            'set label "O" at graph -0.090, -0.100 font "Garamond bold,20"'/&
            'set label "C" at graph -0.080, -0.100 font "Garamond bold,20"')
    endif
    write(21,133)
133 format('# if the value after solid is 0 the monovariants are transparent'/&
         'set style fill transparent solid 1'/&
         'set style line 1 lt 2 lc rgb "#000000" lw 2 pt 10'/&
         'set style line 2 lt 2 lc rgb "#00C000" lw 2 pt 2'/&
         'set style line 3 lt 2 lc rgb "#4169E1" lw 2 pt 7'/&
         'set style line 4 lt 2 lc rgb "#FF0000" lw 2 pt 3'/&
!         'set style line 5 lt 2 lc rgb "#8F8F8F" lw 2 pt 4'/&
!        'set style line 5 lt 2 lc rgb "#FF4500" lw 2 pt 4'/&
          'set style line 5 lt 2 lc rgb "#00FFFF" lw 2 pt 10'/&
         'set style line 6 lt 2 lc rgb "#0080FF" lw 2 pt 5'/&
         'set style line 7 lt 2 lc rgb "#804080" lw 2 pt 6'/&
!         'set style line 7 lt 2 lc rgb "#FF4500" lw 2 pt 6'/&
         'set style line 8 lt 2 lc rgb "#00C000" lw 2 pt 8'/&
         'set style line 9 lt 2 lc rgb "#C0C0C0" lw 2 pt 1'/&
        'set style line 10 lt 2 lc rgb "#DAA520" lw 2 pt 4')
!         'set style line 10 lt 2 lc rgb "#00FFFF" lw 2 pt 10'/&
! orange is #FF4500
! goldenrod hex: DAA520"
! line style 11 is monovariant (not invariant), 12 tieline
!         'set style line 11 lt 2 lc rgb "goldenrod" lw 3'/&
!         'set style line 11 lt 2 lc rgb "#DAA520" lw 3'/&
!         'set style line 12 lt 2 lc rgb "goldenrod" lw 1')
!         'set style line 11 lt 2 lc rgb "#804080" lw 3'/&
!         'set style line 12 lt 2 lc rgb "#804080" lw 1')
!         'set style line 11 lt 2 lc rgb "#7CFF40" lw 3'/&
! for monovariants use filledcurves fc "#xxxxxx" AND linestyle 11
    write(21,134)tielinecolor,tielinecolor
! line style 11 is for the limits of the monovariants, 12 for tie-lines
134 format('set style line 11 lt 2 lc rgb "#',a,'" lw 3'/&
         'set style line 12 lt 2 lc rgb "#',a,'" lw 1')
!         'set style line 12 lt 2 lc rgb "#7CFF40" lw 1')
! The last two styles (11 and 12) are for monovariants (not invariants)
!   and tielines
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
       if(plotgt) then
          xxx=textlabel%xpos+0.5D0*textlabel%ypos
!          sqrt3=0.5D0*sqrt(3.0D0)
          yyy=sqrt3*textlabel%ypos
       else
          xxx=textlabel%xpos
          yyy=textlabel%ypos
       endif
!       write(*,*)'SMP text: ',textlabel%textline,textlabel%xpos,xxx,plotgt
!----------------------- new
       rotate=' '
       if(textlabel%angle.ne.0) write(rotate,177)textlabel%angle
177    format(' rotate by ',i5)
       labelfont=' '
!       write(*,*)'textfontscale: ',textlabel%textfontscale
       if(textlabel%textfontscale.ne.one) then
          write(labelfont,178)int(10*textlabel%textfontscale)
178       format(' font "Sans,',i2,'" ')
       endif
!       if(textlabel%angle.eq.0) then
!       write(21,1505)trim(textlabel%textline),textlabel%xpos,textlabel%ypos,&
       write(21,1505)trim(textlabel%textline),xxx,yyy,&
            trim(labelfont),trim(rotate)
1505   format('set label "',a,'" at ',1pe12.4,', ',1pe12.4,a,a)
!       else
!         write(21,1506)trim(textlabel%textline),textlabel%xpos,textlabel%ypos,&
!               textlabel%angle
!1506      format('set label "',a,'" at ',1pe12.4,', ',1pe12.4,&
!               ' rotate by ',i5)
!       endif
!       textlabel=>textlabel%nexttextlabel
! OLD below
!       if(textlabel%angle.eq.0) then
!          write(21,160)trim(textlabel%textline),xxx,yyy,&
!160       format('set label "',a,'" at ',1pe14.6,', ',1pe14.6)
!       else
!          write(21,161)trim(textlabel%textline),xxx,yyy,textlabel%angle
!161      format('set label "',a,'" at ',1pe14.6,', ',1pe14.6,&
!               ' rotate by ',i5)
!       endif
       jj=len_trim(textlabel%textline)
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
       appgt=.false.
! copy all lines up to "plot" to new graphics file
       nnv=0
200   continue
       read(appfil,210,end=290)appline
210   format(a)
!       if(appline(1:1).eq.'#') then
!          write(*,*)'input: ',trim(appline),appgt
!       endif
       if(appline(1:15).eq.'# GIBBSTRIANGLE') then
          appgt=.true.
          if(.not.plotgt) then
             write(*,*)'Append file is in Gibbstriangle format'
             goto 280
          endif
       endif
!------------------------------------------------------------------
! ignore some lines with "set" in the append file
! set title
! set xlabel
! set ylabel
! set xrange
! set yrange
! set output
! set terminal
! set size
! set key
       if(appline(1:10).eq.'set title ' .or.&
            appline(1:11).eq.'set xlabel ' .or.&
            appline(1:11).eq.'set ylabel ' .or.&
            appline(1:11).eq.'set xrange ' .or.&
            appline(1:11).eq.'set yrange ' .or.&
            appline(1:11).eq.'set output ' .or.&
            appline(1:13).eq.'set terminal ' .or.&
            appline(1:11).eq.'set origin ' .or.&
            appline(1: 9).eq.'set size ' .or.&
            appline(1: 8).eq.'set key ') then
!          write(*,*)'ignoring append line ',trim(appline)
          goto 200
       endif
!------------------------------------------------------------------
       if(index(appline,'plot "-"').gt.0 .or. &
            index(appline,"plot '-'").gt.0 .or. &
            index(appline,"plot for ").gt.0) then
! this is ocplot3B, reading plot command lines
          if(plotgt) then
! check if append file has square or triangular coordinates ...
             if(.not.appgt) then
! If GIBBSTRIANGLE they must be converted unless already a triangle ....
                write(*,*)'Please use append file with Gibbs triangle'
                goto 280
             endif
          elseif(appgt) then
! a triangular append file must be transformed to square ...
             write(*,*)'Please use append file with square coordinates'
             goto 280
          endif
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
             if(ic.ge.mofapl) then
                write(*,*)'Too many head lines in append file',ic
             else
                applines(ic)=appline
             endif
             goto 230
          endif
          nofapl=ic
          goto 290
       else
! ignore all lines until "plot "-" ...
          goto 200
       endif
! These are coordinate lines
       write(21,210)trim(appline)
       nnv=nnv+1
       goto 200
! error oppening append file
280    continue
       write(kou,*)' *** Cannot open or read the append file, skipping it'
       close(appfil)
       appfil=0
290    continue
! do not close the append file, we have to read the data also!
!       write(*,*)'Finished reading appendfile head: ',nofapl,ic
    endif
! coordinate the content of lid with the colors
!    do ii=1,same
!       write(*,*)'phases: ',ii,trim(lid(1,ii)),' ',trim(lid(2,ii))
!    enddo
    nnv=0
    iz=0
!    color(1)=' '
    color=' '
!    if(2*same.gt.maxcolor) then
!       write(*,*)'Number of lines: ',2*same,maxcolor
!    endif
    pair: do jj=1,2
! maybe same must be incremented with the number of tieline blocks??
! and ivariants??  They have separate plot commands ...
       point: do ii=1,same
          iz=iz+1
! plotkod -1 negative means ignore
! plotkod -100 and -101 used for tie-lines
          if(jj.eq.2 .and. plotkod(iz).eq.-1) then
             write(*,*)'Ignoring this line ',jj,iz,plotkod(iz)
!             cycle pair
             cycle point
          endif
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
! skip colors 11 and 12, reserved for monovariants and tielines
          nnv=nnv+1
          lcolor(iz)=nnv
          color(nnv)=lid(jj,ii)
!          write(*,293)'color select: ',nnv,jj,ii,iz,trim(color(nnv))
293       format(a,4i5,' "',a,'"')
295       continue
       enddo point
    enddo pair
!    write(*,*)'Finished assigning colors',iz
! replace _ by - (in phase names).  nnv is number of titles 
    do kk=1,nnv
297    continue
       jj=index(color(kk),'_')
       if(jj.gt.0) then
          color(kk)(jj:jj)='-'
          goto 297
       endif
    enddo
!---------------------------------------------------------------
! check for monovariants (not invariant) and tieline and replace color!
    do ii=1,2*same
       if(lcolor(ii).le.0) then
          write(*,*)'missing color in ',ii,' out of ',2*same
!       else
!          write(*,*)'original: ',ii,lcolor(ii),trim(color(lcolor(ii)))
       endif
    enddo
    lcolor1: do ii=1,2*same
       jj=lcolor(ii)
       if(trim(color(jj)).eq.'invariant') then
          write(*,*)' *** wrong: invariant color ',jj,trim(color(jj))
       endif
       if(trim(color(jj)).eq.'monovariant') then
!          write(*,*)'found monovariant ',ii,jj,2*same
          lcolor(ii)=11
!          color(11)='invariant'
          color(11)='monovariant'
          do kk=ii+1,2*same
             if(lcolor(kk).eq.jj) then
!                write(*,*)'subsequent: ',kk,lcolor(kk),jj,11
                lcolor(kk)=11
             endif
          enddo
! why exit?
!          exit lcolor1
       endif
    enddo lcolor1
    lcolor2: do ii=1,2*same
       jj=lcolor(ii)
       if(jj.le.0 .or. jj.gt.maxcolor) then
! this is a line that should not be plotted ...
          write(*,*)'smp2B: problem: ',ii,jj
          lcolor(ii)=11
          cycle lcolor2
       endif
       if(trim(color(jj)).eq.'tieline') then
!          write(*,*)'found tie-line ',ii,jj
          lcolor(ii)=12
          color(12)='tie-line'
          do kk=ii+1,2*same
             if(lcolor(kk).eq.jj) lcolor(kk)=12
          enddo
! why exit?
!          exit lcolor2
       endif
    enddo lcolor2
!    do ii=1,2*same
!       write(*,*)'Final: ',ii,lcolor(ii),trim(color(lcolor(ii)))
!    enddo
!----------------------------------------------------------------
    if(plotgt) then
! convert all coordinates to a Gibbs trangle, ax and ay are the square coordin
! x = ax + 0.5*ay
! y = 0.5*sqrt(3)*ay
!
! As the 3rd point of the monovariant is connected to the xval/yval I must
! make the same loop as below when plotting ...
!       write(*,*)'Converting coordinates to Gibbs Triangle',same,lineends(1)
       foundinv=0
!       sqrt3=0.5D0*sqrt(3.0D0)
       sumpp=0
       do jj=1,same
          if(sumpp+1.eq.lineends(jj)) then
             foundinv=foundinv+1
!             write(*,*)'Monovariant at ',foundinv,sumpp
             sumpp=sumpp+1
             xval(1,sumpp)=xval(1,sumpp)+5.0D-1*yval(1,sumpp)
             xval(2,sumpp)=xval(2,sumpp)+5.0D-1*yval(2,sumpp)
             zval(1,foundinv)=zval(1,foundinv)+5.0D-1*zval(2,foundinv)
             zval(2,foundinv)=sqrt3*zval(2,foundinv)
             yval(1,sumpp)=sqrt3*yval(1,sumpp)
             yval(2,sumpp)=sqrt3*yval(2,sumpp)
          else
             do while(sumpp.lt.lineends(jj))
                sumpp=sumpp+1
                xval(1,sumpp)=xval(1,sumpp)+5.0D-1*yval(1,sumpp)
                xval(2,sumpp)=xval(2,sumpp)+5.0D-1*yval(2,sumpp)
                yval(1,sumpp)=sqrt3*yval(1,sumpp)
                yval(2,sumpp)=sqrt3*yval(2,sumpp)
             enddo
          endif
       enddo
    endif
!----------------------------------------------------------------
! text in lower left corner
    ii=len_trim(graphopt%lowerleftcorner)
    if(graphopt%gibbstriangle) then
       if(ii.gt.3) then
          write(21,208)trim(graphopt%lowerleftcorner),-0.14
208       format('set label "',a,'" at graph ',F10.4,', -0.05 ')
       elseif(ii.gt.0) then
          write(21,208)trim(graphopt%lowerleftcorner),-0.08
       endif
    elseif(ii.gt.0) then
! in square diagram below figure
       write(21,209)trim(graphopt%lowerleftcorner)
209    format('set label "',a,'" at graph -0.10, -0.08 ')
    endif
! if lowerleftcorner is empty ignore it
!----------------------------------------------------------------
!----------------------------------------------------------------
! Finished all options, now deal with the data to plot!
! this is subroutine ocplot3B for two extensive axis
!----------------------------------------------------------------
! Here we generate the datafile with coordinates to plot
! if nx1 or ny1 is 1 plot all on other axis versus single axis coordinate
! if nx1=ny1 plot the pairs xval(1..nx1,jj) yval(1..ny1,jj)
!----------------------
    backslash=',\'
! empty line before the plot command
!    write(*,*)'lines: ',same,nofinv
    write(21,*)
! here we should start from the value in graphopt%linett
    ii=0
    kk=graphopt%linett-1
    if(kk.ne.0) then
       write(*,*)'Ignoring manipulation of line colors'
    endif
! if graphopt%linestyle=0 use lines, otherwise linespoints
    if(graphopt%linestyle.eq.0) then
       linespoints='lines'
    else
       linespoints='linespoints'
    endif
    done=0
!    noofmono=0
! Here we write all plot "-" using ... and subsequent "" using ...
    naptitle=0
    xtieline=0
    xmonovariant=0
    kkloop: do kk=1,2
       jjloop: do jj=1,same
          ii=ii+1
          if(ii.eq.1) then
             if(lcolor(ii).eq.11) then
! this is monovariant!
! this is the first plot command, ii=1 so kk must be 1!!
                if(kk.eq.1) then
                   write(21,306)monovariant,backslash
306                format('plot "-" using 1:2 with filledcurves ',&
                        'fc "#',a,'" notitle ',a)
                   xmonovariant=jj
!                   write(*,*)'SMP monovariant 1: ',xmonovariant
                else
! this else branch is impossible, when ii=1 then kk=1 !!! but ...
                   write(21,307)monovariant,trim(color(lcolor(ii))),backslash
307                format('plot "-" using 1:2 with filledcurves ',&
                        'fc "#',a,'" title "',a,'"',a)
! light green                     'fc "#B0FFB0" title "',a,'"',a)
! very faint green                'fc "#F0FFF0" title "',a,'"',a)
! faint yellow                     'fc "#EEFFCC" title "',a,'"',a)
                endif
             elseif(lcolor(ii).eq.12) then
! tie-line
                write(21,308)'lines',lcolor(ii),&
                        trim(color(lcolor(ii))),backslash
                xtieline=jj
!                write(*,*)'SMP xtieline 1: ',xtieline
308             format('plot "-" using 1:2 with ',a,' ls ',i2,' notitle ',a)
             else
! normal line with label
                write(21,309)trim(linespoints),lcolor(ii),&
                     trim(color(lcolor(ii))),backslash
             endif
             naptitle=naptitle+1
             apptitles(naptitle)=lcolor(ii)
309          format('plot "-" using 1:2 with ',a,' ls ',i2,' title "',a,'"',a)
             done(lcolor(1))=1
          else
! all lines except the first plotted here
! the last line for the plot command has no backslash --- except if append file
             if(ii.eq.2*same .and. appfil.eq.0) backslash=' '
! we can only use linestyles 1 to 10 except for monovariants and tie-lines
             fcolor=lcolor(ii)
             if(fcolor.gt.12) then
                fcolor=mod(fcolor,10)
                if(fcolor.eq.0) fcolor=10
! fixed Nath MoNiRe isotherm at 1500 K had some lines with no lcolor assignment!
             elseif(fcolor.le.0) then
                lcolor(ii)=1
                fcolor=1
             endif
             cone: if(done(lcolor(ii)).eq.1) then
! we have already a title for this line ... except tie-lines and monovariant
                if(lcolor(ii).eq.11) then
                   if(kk.eq.1) then
! first time plotting an invariant use thick lines
                      write(21,320)'lines',monovariantborder,backslash
! save the index of the last monovariant to add title!
                      xmonovariant=jj
!                      write(*,*)'SMP monovariant 2: ',xmonovariant
                   else
! if kk=2 check if this is last monovariant, if so add title
                      if(jj.eq.xmonovariant) then
                         write(21,318)monovariant,trim(color(11)),backslash
318                      format('"" using 1:2 with filledcurves ',&
                              'fc "#',a,'" title "',a,'" ',a)
!                         write(*,*)'SMP monovariant 5: ',xmonovariant,jj
                      else
                         write(21,319)monovariant,backslash
319                      format('"" using 1:2 with filledcurves ',&
                              'fc "#',a,'" notitle ',a)
                      endif
                   endif
                elseif(lcolor(ii).eq.12) then
! tie-line, if kk==2 and xtieline==jj add label
                   if(kk.eq.1) then
                      write(21,320)'lines',fcolor,backslash
                      xtieline=jj
!                      write(*,*)'SMP xtieline 2: ',xtieline
                   elseif(xtieline.ne.jj) then
                      write(21,320)'lines',fcolor,backslash
                   else
                      write(21,299)'lines',fcolor,trim(color(12)),backslash
299                   format('"" using 1:2 with ',a,' ls ',i2,&
                           ' title "',a,'" ',a)
!                      write(*,*)'SMP xtieline 5:',jj,xtieline
                   endif
                else
! normal line with no title
                   write(21,320)trim(linespoints),fcolor,backslash
                endif
320             format('"" using 1:2 with ',a,' ls ',i2,' notitle ',a)
             else 
! we have a new line withou title
                if(fcolor.eq.11) then
                   if(kk.eq.1) then
! first time plotting a monovariant use thick lines
                      write(21,320)'lines',monovariantborder,backslash
                      xmonovariant=jj
!                      write(*,*)'SMP monovariant 3: ',xmonovariant
                   else
                      if(jj.eq.xmonovariant) then
! this is the last monovariant, add title
                         write(21,321)monovariant,&
                              trim(color(lcolor(ii))),backslash
321                      format('"" using 1:2 with filledcurves ',&
                              'fc "#',a,'" title "',a,'"',a)
!                         write(*,*)'SMP monovariant 4: ',xmonovariant,jj
                      else
! not the last monovariant, notitle
                         write(21,325)monovariant,backslash
325                      format('"" using 1:2 with filledcurves ',&
                              'fc "#',a,'" notitle ',a)
                      endif
                   endif
                elseif(fcolor.eq.12) then
! this is a tie-line without title
                   if(kk.eq.1) then
! if kk=1 no not add title, just keep track of last tie-line
                      write(21,320)'lines',fcolor,backslash
                      xtieline=jj
!                      write(*,*)'SMP xtieline 3: ',xtieline
                   else
! if kk=2 add title if xtieline=jj
                      if(jj.eq.xtieline) then
                         write(21,331)'lines',fcolor,&
                              trim(color(lcolor(ii))),backslash
!                         write(*,*)'SMP xtieline 4:',jj,xtieline
                      else
                         write(21,320)'lines',fcolor,backslash
                      endif
                   endif
                else
! any normal line add title
                   write(21,331)trim(linespoints),fcolor,&
                        trim(color(lcolor(ii))),backslash
                endif
                naptitle=naptitle+1
                apptitles(naptitle)=lcolor(ii)
331             format('"" using 1:2 with ',a,' ls ',i2,' title "',a,'"',a)
                done(lcolor(ii))=1
             endif cone
          endif
       enddo jjloop
    enddo kkloop
! if we have an append file we must add the plotcommands in applines(1:nofapl)
! we made sure a few lines above that there is a backslash at last line above
    if(appfil.gt.0) then
! match the titles of all applines with those used above in lcolor
! and make sure the line with the same title has the same color
! note there is only one color for monovariants ... problem with lcolor ...
!       call ocappfixlabels(nofapl,applines,same,color,lcolor,nnv)
       call ocappfixlabels(nofapl,applines,same,color,apptitles,nnv)
! replace the 'plot "-" ' by just '"" ' 
       applines(1)='"" '//applines(1)(9:)
!       write(*,*)'from append file',nofapl,trim(applines(1))
       do ii=1,nofapl
          write(21,'(a)')trim(applines(ii))
       enddo
    endif
!    goto 500
! loop for all line coordinates
!    sumpp=0
    do kk=1,2
!       again=sumpp
       foundinv=0
       sumpp=0
       do jj=1,same
! handle monovariants, just one point, just once
! the monovariants (not invariants) will be plotted twice ...
          if(sumpp+1.eq.lineends(jj)) then
             foundinv=foundinv+1
             sumpp=sumpp+1
!             write(*,*)'Assumed to be an monovariant!',foundinv,kk
             write(21,600)jj,lcolor(jj)
600          format('# Line ',2i5,' representing a monovariant')
             write(21,549)xval(1,sumpp),yval(1,sumpp)
             write(21,549)xval(2,sumpp),yval(2,sumpp)
             write(21,549)zval(1,foundinv),zval(2,foundinv)
             write(21,549)xval(1,sumpp),yval(1,sumpp)
! we are at the end of a line, write a blank line
             write(21,548)jj
548          format('e '//'# end of monovariant',2i5)
!548          format('e '//'# end of invariant',2i5)
          else
! this is the beginning of a line to be plotted
             if(lcolor(jj).eq.12) then
                write(21,605)jj,lcolor(jj)
605             format('# Line ',2i5,' representing tielines')
             else
                write(21,610)jj,lcolor(jj),trim(color(lcolor(jj)))
610             format('# Line ',2i5,' representing phase: ',a)
             endif
             do while(sumpp.lt.lineends(jj))
                sumpp=sumpp+1
                write(21,549)xval(kk,sumpp),yval(kk,sumpp)
!549             format(2e15.6,4i7)
549             format(2(1pe16.6),4i7)
! plotkod -101 means tieline
! UNFINISHED: VALGRIND indicates plotkod(sumpp) is uninitiallized??
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
       if(appline(1:11).eq.'pause mouse') then
          write(*,*)'reading appendfile ends at "puase mouse"'
          goto 1910
       else
          write(21,884)trim(appline)
          goto 1900
       endif
1910   continue
!       write(*,*)'Appended ',ic,' data lines'
       close(appfil)
       appfil=0
    endif
!------------------------------------------------------
!    if(pform(1:1).eq.' ') then
    if(graphopt%gnutermsel.eq.1) then
! if not hardcopy pause gnuplot.  Mouse means clicking in the graphics window
! will close it. I would like to have an option to spawn the graphics window...
! so it is kept while continuing the program.
       write(21,990)trim(graphopt%plotend)
990    format(a)
!990    format('pause mouse')
!990    format('e'//'pause mouse')
    else
! add pause mouse as comment
       write(21,991)
991    format('# pause mouse')
    endif
    close(21)
    if(appfil.ne.0) close(appfil)
    appfil=0

!    write(21,565)
!565 format('e'//'pause mouse'/)
!    close(21)
!
!    gnuplotline='gnuplot ocgnu.plt '
    gnuplotline='gnuplot '//trim(pfc)//' & '
! if gnuplot cannot be started with gnuplot give normal path ...
!    gnuplotline='"c:\program files\gnuplot\bin\wgnuplot.exe" '//pfc(1:kkk)//' '
    k3=len_trim(gnuplotline)+1
    write(*,*)'Gnuplot command line: ',trim(gnuplotline)
!    if(pform(1:1).ne.' ') then
    if(graphopt%gnutermsel.ne.1) then
       write(*,*)'Graphics output file: ',trim(pfh)
    endif
! grwin set by compiler option, 1 means windows 
    if(grwin.eq.1) then
       if(btest(graphopt%status,GRKEEP)) then
! this is a TERNARY PLOT with 2 extenive axis
!          write(*,*)'executing command '//trim(gnuplotline(9:))
!          call system(gnuplotline(9:))
          write(*,*)'Executing Command: "start /B '//trim(gnuplotline)//'"'
! WORKS WITH OCPLOT3B
          call execute_command_line('start /B '//trim(gnuplotline))
       else
!          write(*,*)'executing command '//trim(gnuplotline)
          call execute_command_line(gnuplotline)
       endif
    else
! plot on non-windows system
! how to implement GRKEEP?
       write(*,*)'executing command '//trim(gnuplotline)
       call execute_command_line(gnuplotline)
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
       if(k2.le.0) then
! full has no compset
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
    logical once
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
105          format('  Line ',i3,' with ',i5,&
                  ' equilibria ended with error: ',i6)
          else
             write(kou,110)mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria
110          format('  Line ',i3,' with ',i5,&
                  ' equilibria ending at axis limit')
          endif
       else
          ll=mapnode%linehead(kl)%end%seqx
          write(kou,120)mapnode%linehead(kl)%lineid,&
               mapnode%linehead(kl)%number_of_equilibria,ll
120       format('  Line ',i3,' with ',i5,' equilibria ending at node ',i3)
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
          once=.true.
          write(kou,140)
140       format('Saved ceq     link       T            X')
          ceqloop: do while(ll.gt.0)
             thisceq=>results%savedceq(ll)
             do jax=1,nax
                call locate_condition(axarr(jax)%seqz,pcond,thisceq)
                if(gx%bmperr.ne.0) goto 300
                svrrec=>pcond%statvar(1)
                call state_variable_val(svrrec,axxx(jax),thisceq)
                if(gx%bmperr.ne.0)then
                   if(once) then
                   write(*,*)' *** Error ',gx%bmperr,&
                        ' reset, data may be missing'
                      once=.false.
                   endif
                   gx%bmperr=0
                endif
             enddo
             write(kou,150)ll,thisceq%nexteq,thisceq%tpval(1),axxx
150          format(2i9,f9.2,5(1pe13.5))
             ll=thisceq%nexteq
          enddo ceqloop
       endif
300    continue
       if(gx%bmperr.ne.0) then
          write(*,*)' *** Error ',gx%bmperr,' reset, data maybe missing'
          gx%bmperr=0
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
105          format('  Line ',i3,' with ',i5,&
                  ' equilibria ended with error: ',i6/2x,a,' with phases: ',a)
          else
             write(kou,110)mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria,status,trim(phline)
110          format('  Line ',i3,' with ',i5,&
                  ' equilibria ending at axis limit.'/2x,a,' with phases: ',a)
          endif
       else
          ll=mapnode%linehead(kl)%end%seqx
          write(kou,120)mapnode%linehead(kl)%lineid,&
               mapnode%linehead(kl)%number_of_equilibria,ll,status,trim(phline)
120       format('  Line ',i3,' with ',i5,' equilibria ending at node ',i3/&
               2x,a,' with phases: ',a)
       endif
!       write(*,*)'mapnode%linehead%first: '
       ll=mapnode%linehead(kl)%first
!       write(*,*)'axis conditions: ',ll,axarr(1)%seqz,axarr(2)%seqz
!       if(.not.btest(mapnode%linehead(kl)%status,EXCLUDEDLINE)) then
!          if(.not.allocated(results%savedceq)) then
!             write(*,*)'Cannot find link to saved equilibria! '
!          else
!             write(*,*)'Loop for all equlibria'
!             do while(ll.gt.0)
!                thisceq=>results%savedceq(ll)
!                write(*,*)'thisceq pointer set',ll
!                do jax=1,nax
!                   call locate_condition(axarr(jax)%seqz,pcond,thisceq)
!                   if(gx%bmperr.ne.0) goto 1000
!                   write(*,*)'local_condition NOT OK'
!                   svrrec=>pcond%statvar(1)
!                   call state_variable_val(svrrec,axxx(jax),thisceq)
!                   if(gx%bmperr.ne.0) goto 1000
!                   write(*,*)'state variable found OK'
!                enddo
!                write(kou,150)ll,thisceq%next,thisceq%tpval(1),axxx
!150             format('     Saved ceq ',2i5,f8.2,5(1pe14.6))
!                ll=thisceq%next
!             enddo
!          endif
!       endif
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

!\begin{berbatim}
  subroutine ocappfixlabels(nofapl,applines,same,color,appcol,nnv)
! check if there are the same labels in applines.  They are normally phase
! names and should be the same !! use the same color for the same phase!!
! if we find a title matching we must change all lines with the same ls x
! to the value of x in the new file !!
! nofapl:   number of lines to plot in appfile 
! applines: the plot command lines from appfile
! same:     number of lines to plot in the current calculation
! color:    the titles in the current calculation
! appcol:   the linestyle for color in the current calculation
! clast:    the last used linestyle in the current calculation
    implicit none
    integer nofapl,same,nnv
    character*(*) applines(nofapl),color(*)
    integer appcol(*)
!\end{berbatim}
    integer i1,j1,k1,nols,ip,jp,found,oldls,newls
    integer, parameter :: mofapl=100
    integer nyttls(nofapl)
    character endofline*(2),title*24
    integer changels(2,nofapl)
! the line label is color(lcolor(jj)) ... suck.  Use the index of color as nyls
!    write(*,*)'ocappfixlabels: ',nofapl,same,nnv
!    do j1=1,same
!       if(color(j1)(1:1).ne.' ') write(*,'(i3,2x,a,i4)')j1,trim(color(j1)),&
!            appcol(j1)
!    enddo
!    do i1=1,nofapl
!       write(*,'(i3,2x,a)')i1,trim(applines(i1))
!    enddo
! we have to check if there are new phases in the appfile.  They should
! be given color/labels that start from nnv
! unique for the appfile.  They should be given colors nnv and higher!
!
! still problem for CHI with the bef-500-gibbs testfile
!
! The applines look like:
! "" using 1:2 with filledcurves fc "#EEFFCC" title "monovariant",\
! "" using 1:2 with lines ls  4 title "BCC-A2",\
! "" using 1:2 with lines ls 12 title "tie-line",\
! "" using 1:2 with lines ls  4 notitle ,\
! each title (in color) is associated with a ls value (in lcolor)
! if we find a matching title in applines change the ls value to lcolor
! for all lines with the ls value
! for titles in applines with no title in color use last free ls
    endofline=',\'
    found=0
    newls=nnv
    loop1: do j1=1,nofapl
       if(j1.eq.nofapl) endofline=' '
! check if the line has a title
       ip=index(applines(j1),' title ')
       if(ip.gt.0) then
!          write(*,*)' *** Found title: ',trim(applines(j1)),j1
          title=applines(j1)(ip+8:)
          jp=index(title,'"')
          if(jp.eq.0) then
             write(*,*)'Missing ": ',trim(applines(j1)),j1
             stop
          else
             title(jp:)=' '
          endif
          if(trim(title).eq.'tie-line' .or.&
               trim(title).eq.'monovariant') then
! titles "tie-line" and "monovariant" are just replaced by notitle
             applines(j1)(ip:)=' notitle'//endofline
!             write(*,*)'removed tie/mono: ',trim(applines(j1)),j1
             cycle loop1
          endif
! compare with color(1..same)
!          write(*,*)'Comparing with old labels'
          loop2: do i1=1,same
             if(color(i1)(1:1).eq.' ') cycle loop2
             if(trim(color(i1)).eq.trim(title)) then
! we have found the applines title in current labels
!                write(*,*)'Found title: ',trim(title)
                applines(j1)(ip:)=' notitle'//endofline
!                write(*,*)'removed: ',trim(applines(j1)),j1
! we must also save the value after ' line ls ' to change other lines
! "" using 1:2 with lines ls  4 title "BCC-A2",\
                ip=index(applines(j1),' lines ls ')
                if(ip.gt.0) then
                   jp=ip+9
                   ip=ip+10
                   call getint(applines(j1),jp,oldls)
                   if(buperr.ne.0) then
                      write(*,*)'No ls number: ',trim(applines(j1))
                      stop
                   endif
                   found=found+1
                   changels(1,found)=oldls
                   changels(2,found)=i1
                   write(applines(j1)(ip:ip+1),'(i2)')i1
!                   write(*,69)'Changing ls: ',trim(applines(j1)),oldls,i1,found
69                 format(a,a,5i4)
                else
                   write(*,*)'missing "lines ls" in: ',trim(applines(j1))
                   stop
                endif
                cycle loop1
             endif
!             write(*,*)'No match: ',trim(title),' ',trim(color(i1)),i1,same
          enddo loop2
! if we come here we have a new title in the appfiles,
! that should be assigned newls
          write(*,*)'New label in appfiles: ',trim(applines(j1))
          stop
       else
! this is a line without title but we may have to change number after "line ls"
          ip=index(applines(j1),' lines ls ')
          if(ip.gt.0) then
!             write(*,'(a,10(i4,i3))')'changels: ',&
!                  (changels(1,k1),changels(2,k1),k1=1,found)
             jp=ip+9
             ip=ip+10
             call getint(applines(j1),jp,oldls)
             if(buperr.ne.0) then
                write(*,*)'Cannot find ls number: ',trim(applines(j1)),found
                stop
             endif
!             write(*,*)'Found ls number: ',oldls
! ignore 11 and 12 as they are tie-lines or monovariant
             if(oldls.eq.11 .or. oldls.eq.12) cycle loop1
! seach for replacement
             getls: do k1=1,found
!                if(changels(1,k1).ne.oldls) cycle getls
                if(changels(1,k1).eq.oldls) goto 100
             enddo getls
! if k1>found then we have not found oldls
             if(k1.gt.found) then
                write(*,79)'Cannot find oldls: ',oldls,k1,found,j1,&
                     trim(applines(j1))
79              format(a,4i3,' in ',a)
                write(*,'(10(i2,i3))')(changels(1,k1),changels(2,k1),k1=1,found)
                stop
             endif
! write the new ls number in applines(j1)
100          continue
             write(applines(j1)(ip:ip+1),'(i2)')changels(2,k1)
!             write(*,*)'Changed ls: ',trim(applines(j1))
!          else
!             write(*,*)'skipping: ',trim(applines(j1))
          endif
       endif
    enddo loop1
!    write(*,*)'New applines'
!    do j1=1,nofapl
!       write(*,*)j1,trim(applines(j1))
!    enddo
1000 continue
    return
  end subroutine ocappfixlabels
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

