! included in smp2.F90.  Generating graphics using GNUPLOT

!\addtotable subroutine ocplot2
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
    character pltax(2)*64,filename*128
!
    type(map_ceqresults), pointer :: results
    TYPE(gtp_equilibrium_data), pointer :: curceq
    type(map_node), pointer :: mapnode,invar,localtop
    type(map_line), pointer :: mapline
    logical wildcard,hashtag
    character ch1*1,gnuplotline*256,pfd*128,pfc*256
    character pfh*128,dummy*24
    double precision, dimension(:,:), allocatable :: anp
    double precision, dimension(:), allocatable :: xax,yyy
! save isopleth invariants in special array, max 50 invariants
    double precision xyinv(4,50)
    integer ninv
! Too big??
!    integer, parameter :: maxval=10000
! plotting isothermal section Cr-Fe-Mo required more than 2000
    integer, parameter :: maxval=4000
    integer, dimension(:), allocatable :: nonzero,linzero,linesep
!    integer, dimension(:), allocatable :: linesep
! encoded2 stores returned text from get_many ... 2048 is too short ...
! selphase used when plotting data just for a selected phase like y(fcc,*)
    character statevar*64,encoded1*1024,encoded2*4096,selphase*24,funsym*24
    character*128, dimension(:), allocatable :: phaseline
    integer i,ic,jj,k3,kk,kkk,lokcs,nnp,np,nrv,nv,nzp,ip,nstep,nnv,nofapl
    integer nr,line,next,seqx,nlinesep,ksep,iax,anpax,notanp,appfil,errall
    double precision xmax,xmin,ymax,ymin,value,anpmin,anpmax
! used for Scheil
    double precision npflval
    logical scheilorder
! lhpos is last used position in lineheader
    integer giveup,nax,ikol,maxanp,lcolor,lhpos,repeat,anpdim,qp
    integer nix,stoichfix,invlines,invnode,nrett,mfix
    integer, allocatable, dimension(:) :: ixpos
! setting color on isopleth lines?  Dimension is max different fix phases
    integer, allocatable, dimension(:,:) :: phamfu
    integer fixphasecolor
! trying to understand
    integer ttunodeid,ttuheads,ttutoplines,ttuline,ttuplotline,haha
    character date*8,mdate*12,title*128,backslash*2,lineheader*1024
    character deftitle*128,labelkey*64
    logical overflow,first,last,novalues,selectph,varofun,moretops,isopleth
    logical, allocatable, dimension(:) ::  nevernone
! dot derivatives should not be calcuated at first point of a range
    logical skipdotder
! textlabels
    type(graphics_textlabel), pointer :: textlabel
! line identification (title)
    character*16, dimension(:), allocatable :: lid
!    character*32, dimension(:), allocatable :: lid
!
!    write(*,*)'In ocplot2 graphopt%status: ',maptop%status,MAPINVARIANT
! transfer from graphics record to local variables
! initiate lines_excluded
    lines_excluded=0
    scheilorder=.FALSE.
! create the terminal plot_line record
    allocate(lastplotline)
    nullify(lastplotline%nextline)
    lastplotline%type=-1
    plotline1=>lastplotline
! when creating a new plotline: ??
! 1: allocate(plotline%nextline)
! 2: plotline%nextline%nextline=>plotline1
! 3: ploline1=>plotline1%nextline
! transfer from graphics record to local variables
    pltax(1)=graphopt%pltax(1)
    pltax(2)=graphopt%pltax(2)
    isopleth=btest(graphopt%status,GRISOPLETH)
!    write(*,*)'ocplot2 wildcard: ',trim(pltax(1)),' & ',trim(pltax(2))
!    if(index(pltax(1),'*').gt.0 .or. index(pltax(2),'*').gt.0) then
! fixed in PMON6
! allow plotting phase compositions also for isopleths ...
!       isopleth=.FALSE.
!    endif
    if(isopleth) write(*,*)'smp2b plotting isopleth'
    filename=graphopt%filename
    funsym=' '
! for isopleths this value determine the line color
    fixphasecolor=1
! If wildcard on two axis use ocplot3 to extract data (tie-lines in plane)
    if(index(pltax(1),'*').gt.0 .and. index(pltax(2),'*').gt.0) then
!       write(*,*)'Using ocplot3'
       call ocplot3(ndx,pltax,filename,maptop,axarr,graphopt,&
            version,ceq)
       goto 1000
    endif
! for tzero lines there is no meqrec record, meqrec%phr not allocated
!    write(*,*)'In ocplot2: ',maptop%lines,allocated(maptop%linehead)
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
11     format('SMP limits set by user for ',a,': ',2(1pe14.6))
    endif
    if(graphopt%rangedefaults(2).ne.0) then
       write(*,11)'y',graphopt%plotmin(2),graphopt%plotmax(2)
    endif
! allocate as many items in linesep as there are mapnodes.
! Hm, if merging plots the number of separators needed can be any value
    jj=100+10*maptop%next%seqx+1
!    write(*,*)'SMP: Allocating linesep: ',jj
    allocate(linesep(jj),stat=errall)
    if(errall.ne.0) then
       write(*,*)'SMP2B Allocation error 1: ',errall
       gx%bmperr=4370; goto 1000
    endif
    nax=maptop%number_ofaxis
    linesep=0
! allocate texts to identify the lines on the gnuplot file
!    write(*,*)'SMP: Allocating phaseline: ',jj
    allocate(phaseline(jj),stat=errall)
    allocate(phamfu(2,jj),stat=errall)
    if(errall.ne.0) then
       write(*,*)'SMP2B Allocation error 2: ',errall
       gx%bmperr=4370; goto 1000
    endif
! zero array fr isopleth fix phase
    phamfu=0
!    if(maptop%number_ofaxis.gt.1) then
!       write(*,*)'Warning: may not not handle map graphics correctly',jj
!    endif
!
    giveup=0
    nrv=maxval
!    write(*,*)'SMP: allocating xax: ',nrv
    allocate(xax(nrv),stat=errall)
    if(errall.ne.0) then
       write(*,*)'SMP2B Allocation error 3: ',errall
       gx%bmperr=4370; goto 1000
    endif
! to insert MOVE at axis terminations
    nlinesep=1
    phaseline(1)=' '
    phaseline(2)=' '
! nv is number of values to plot
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
    ninv=0
    wildcard=.FALSE.
    selectph=.FALSE.
    hashtag=.FALSE.
    selphase=' '
    graphopt%specialdiagram=0
    if(maptop%type_of_node.eq.3) then
! this change the order of plotting the lines, maybe needed only for PFL/PFS ??
       scheilorder=.TRUE.
    endif
    do iax=1,2
!       write(*,*)'Allocating for axis: ',iax
       call capson(pltax(iax))
       if(pltax(iax)(1:4).eq.'PFL ' .or. pltax(iax)(1:4).eq.'PFS ') then
          if(maptop%type_of_node.eq.3) then
! this is a function only used for plotting phase fraction liquid or solids
! in Scheil simulations.
             npflval=one
! this indicates to ocplot2B that one must use plot "-" ...etc
! to have different colors and labels on different lines
!             write(*,*)'ocplot2: Setting graphopt%specialdiagram=2'
             graphopt%specialdiagram=2
          else
             write(*,*)'Plot axis PFL/PFS are reserved for Scheil simulations'
             gx%bmperr=4399; goto 1000
          endif
       endif
!       wildcard1: if(index(pltax(iax),'*').gt.0) then
       wildcard1: if(index(pltax(iax),'*').gt.0 .or. &
            index(pltax(iax),'(#)').gt.0) then
! searching for (#) avoids problem when # is used for comp. sets or sublattatice
          i=index(pltax(iax),'#')
          if(i.gt.0 .and.&
               (pltax(iax)(i+1:i+1).eq.')'.or.pltax(iax)(i+1:i+1).eq.',')) then
! this means the phase name is #, indicating all phases including dormant
! Note that # is used to indicate composition sets, thus ignore #2 etc
             hashtag=.TRUE.
!             write(*,*)'SMP2B hastag set true',trim(pltax(iax)),i
          endif
          if(wildcard) then
             write(*,*)'in OCPLOT2 one axis variable with wildcard allowed'
             goto 1000
          endif
! wildcards allowed only on one axis, we do not know how many columns needed
! allocate as many array elements as columns
          anpdim=np
!          write(*,*)'SMP: allocating anp1: ',np*nrv
          allocate(anp(np,nrv),stat=errall)
!          write(*,*)'SMP: allocating anp2: ',np
! nonzero indicates for each column if there is any nonzero value
! columns with only zero values will be eliminated before plotting
          allocate(nonzero(np),stat=errall)
!          write(*,*)'SMP: allocating nonzero: ',np
! linzero indicate for the present block of equilibria for each column
! if this column contain nonzero values
          allocate(linzero(np))
!          write(*,*)'SMP: allocating yyy: ',np
          allocate(yyy(np),stat=errall)
          if(errall.ne.0) then
             write(*,*)'SMP2B Allocation error 5: ',errall
             gx%bmperr=4370; goto 1000
          endif
          nzp=np
! nzp should be dimension of yyy, np returns the number of values in yyy
! yyy is to extract state variable values for the column with wildcard
! NOTE binary phase diagrams are plotted with wildcard axis like x(*,cr) vs T
! nevernone is an attempt to remove columns that are zero by the value NaN
          allocate(nevernone(np),stat=errall)
          if(errall.ne.0) then
             write(*,*)'SMP2B Allocation error 6: ',errall
             gx%bmperr=4370; goto 1000
          endif
          nevernone=.FALSE.
          nonzero=0
          wildcard=.TRUE.
          anpax=iax
! we can have wildcards as np(*), w(fcc,*) or w(*,cr) or y(gas,*)
! IT IS NOT ALLOWED TO HAVE y(*,*) ... only one wildcard
! when we plot things like y(fcc,*) we should only select equilibria
! with fcc stable. Check if * is before or after a ,
! NOTE that a single * without , can be phase or component:
! MU(*) is for a component
          ikol=index(pltax(iax),',')
!          write(*,*)'smp2b selectph: ',trim(pltax(iax)),ikol
          if(ikol.gt.0) then
! if the * is after the , then extract the phase name before             
! and set selecrph to TRUE
             if(pltax(iax)(ikol+1:ikol+1).eq.'*') then
                nrv=index(pltax(iax),'(')
                if(nrv.lt.ikol) then
                   selphase=pltax(iax)(nrv+1:ikol-1)
!                   write(*,*)'SMP2B wildcard selected phase: ',trim(selphase)
                   selectph=.TRUE.
                endif
             endif
!          else
! this is perfectly possible, for example NP(*)
!             write(*,*)'SMP2B: Wildcard without ,!'
          endif
       endif wildcard1
    enddo
    if(.not.wildcard) then
       anpdim=1
       allocate(anp(anpdim,nrv),stat=errall)
       if(errall.ne.0) then
          write(*,*)'SMP2B Allocation error 7: ',errall
          gx%bmperr=4370; goto 1000
       endif
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
!    write(*,*)'SMP2B ocplot2 at label 77A: ',localtop%lines
    ikol=0
    do nrv=1,localtop%lines
       if(allocated(localtop%linehead)) localtop%linehead(nrv)%done=0
    enddo
! we sometimes have a segmentation fault when several maptops ...
    if(associated(localtop%next)) then
       mapnode=>localtop%next
       invnode=0
       if(btest(mapnode%status,MAPINVARIANT)) then
          invnode=size(mapnode%linehead)
!          write(*,*)'ocplot2 invariant node 1',invnode
       endif
    else
       write(*,*)'Mapnode next link missing 1'
       goto 79
    endif
!    write(*,*)'SMP2B ocplot2 at label 77B: ',localtop%lines
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
    if(.not.associated(localtop%saveceq)) then
       write(*,*)'Plot data structure has no results to plot'
       gx%bmperr=4399; goto 1000
    endif
    results=>localtop%saveceq
    mapnode=>localtop
    line=1
! looking for segmentation fault running map11.OCM as part of all.OCM in oc4P
! This error may be due to having created (or not created) composition sets ...
!    write(*,*)'SMP2B ocplot2 after label 79'
! extract the names of stable phases for this lone

!    write(*,*)'mapnode index: ',mapnode%seqx
!    write(*,*)'Before label 100: ',results%free
!------------------------------------------- begin loop 100
! loop back here from ??
100    continue
!       write(*,*)'SMP2B ocplot2 at label 100',localtop%seqx,line,&
!            size(localtop%linehead)
       mapline=>localtop%linehead(line)
! TRYING TO UNDERSTAND WHAT IS HAPPENING HERE ....
       ttunodeid=localtop%seqx
       ttuheads=size(localtop%linehead)
       ttutoplines=localtop%lines
       ttuline=line
       ttuplotline=nv
!       write(*,101)'At 100',ttunodeid,ttuheads,ttutoplines,ttuline,ttuplotline
101    format('Line selected: ',a,': nodeid ',i3,', heads/lines: ',2i2,&
            ' index: ',i2,', plotline: ',2i4)
! initiate novalues to false for each line
       novalues=.false.
! skip first point if dot derivative
       skipdotder=.TRUE.
! We have a segmentation fault sfter this in oc4P when running map11.OCM
! at the end of running all macros.  
!       write(*,*)'In ocplot2, looking for segmentation fault 3'
! skip line if EXCLUDEDLINE set
       if(btest(mapline%status,EXCLUDEDLINE)) then
!          write(*,*)'Skipping a line 3'
          lines_excluded=lines_excluded+1
          if(line.lt.mapnode%lines) then
             line=line+1
             goto 100
          else
             goto 500
          endif
       endif
! We jump here from where?? ... several places
110    continue
! mark line is plotted
!       write(*,*)'Values from mapline ',mapline%lineid
! loop for all calculated equilibria, phases and composition sets
!       write(*,*)'Before label 150: ',mapline%lineid
!150    continue
       nr=mapline%first
       if(mapline%done.ne.0) goto 220
       mapline%done=-1
!       if(ocv()) write(*,*)'Plotting line: ',&
!       write(*,*)'Plotting line: ',&
!            mapline%lineid,mapline%number_of_equilibria,mapline%termerr
!--------------
       ttunodeid=localtop%seqx
       ttuheads=size(localtop%linehead)
       ttutoplines=localtop%lines
       ttuline=line
       ttuplotline=nv
!       write(*,101)'at 110',ttunodeid,ttuheads,ttutoplines,ttuline,&
!            ttuplotline
!--------------
       if(mapline%lineid.le.0) then
          write(*,*)'Skipping line with id less or equal to zero'
          goto 500
       elseif(mapline%number_of_equilibria.le.0) then
!          write(*,*)'Skipping line with no equilibria.'
          goto 500
       endif
       first=.TRUE.
! last set true when we reach the last equilibrium on the line
       last=.FALSE.
! we may have empty lines due to bugs ...
!       write(*,*)'Axis with wildcard and not: ',anpax,notanp
!200    continue
! this is the loop for all equilibria in the line
!       write(*,*)'SMP2: nr and nv: ',nr,nv
! Possibly skip last if mapline%termerr not zero
!       if(mapline%termerr.ne.0) write(*,*)'SMP2B termerr:',mapline%termerr
       nrett=nv+1
       plot1: do while(nr.gt.0)
! nr is index to stored equilibrium
          if(last.and. mapline%termerr.ne.0) then
! skip this equilibrium!!
             nr=0
             write(*,*)'Skipping last point of a line in the plot'
! same as cycle plot1 but maybe safer??
             goto 220
!             cycle plot1
          endif
          nv=nv+1
          if(nv.ge.maxval) then
             write(*,*)'Too many points to plot',maxval
             goto 600
          endif
          curceq=>results%savedceq(nr)
          if(ocv()) write(*,201)'Current equilibrium: ',nr,nv,curceq%tpval(1)
201       format(a,2i5,F8.2,1pe14.6)
! extract the names of stable phases to phaseline from the first equilibrium
! Note that the information of the fix phases have not been saved
          if(first) then
! segmentation fault after here
!             write(*,*)'In ocplot2, looking for segmentation fault 4A'
             first=.false.
             kk=1
! leave space to write fixphasecolor index first
             if(isopleth) kk=5
!             if(selectph) novalues=.TRUE.
             phloop: do jj=1,noph()
!                write(*,*)'In ocplot2, segmentation fault 4Ax',jj,noofcs(jj)
                do ic=1,noofcs(jj)
                   k3=test_phase_status(jj,ic,value,curceq)
                   if(gx%bmperr.ne.0) goto 1000
                   stableph1: if(k3.gt.0) then
! this phase is stable or fix
                      call get_phase_name(jj,ic,dummy)
                      if(gx%bmperr.ne.0) goto 1000
                      if(selectph) then
! this is an attempt to remove lines from irrelevant equilibria when plotting
! data for a specific phase like y(fcc#4,*) which is not stable??
                         if(abbr_phname_same(dummy,trim(selphase))) then
                            novalues=.FALSE.
!                            write(*,217)'SMP2B novalues set FALSE',&
!                                 mapline%lineid,trim(dummy),trim(selphase)
217                         format(a,i5,2x,a,2x,a)
                            exit phloop
                         else
                            novalues=.TRUE.
!                            write(*,*)'SMP2B novalues set TRUE',mapline%lineid
                         endif
                      endif
                      if(.not.novalues) then
! I think this phaseline is no longer used ?? YES it is
!                         write(*,*)'SMP addto phaseline 1: ',trim(dummy),&
!                              nlinesep
                         phaseline(nlinesep)(kk:)=dummy
                         kk=len_trim(phaseline(nlinesep))+2
                      endif
                      linecolor: if(isopleth .and. value.eq.zero) then
! attempt to have the same color for all lines with same fix phase
                         kkk=10*jj+ic
                         do mfix=1,nlinesep
                            if(kkk.eq.phamfu(1,mfix)) then
                               phamfu(1,nlinesep)=kkk
                               phamfu(2,nlinesep)=phamfu(2,mfix)
!                               write(*,*)'smp2b same: ',kkk,phamfu(1,mfix)
                               exit linecolor
                            endif
                         enddo
                         if(mfix.gt.nlinesep) then
                            phamfu(1,nlinesep)=kkk
                            phamfu(2,nlinesep)=fixphasecolor
                            fixphasecolor=fixphasecolor+1
! save phase names in lid for KET, how many phases?
                            if(.not.allocated(lid)) then
                               allocate(lid(20))
                               lid=' '
                            endif
                            lid(fixphasecolor-1)=dummy
                         endif
!                         write(*,'(a,7i5)')'smp2b phasecolor: ',nlinesep,mfix,&
!                              kkk,jj,fixphasecolor
                      endif linecolor
                   endif stableph1
                enddo
             enddo phloop
! finding place to change color of line in Scheil simulations
!             write(*,'(a,i3,2x,a)')'ocplot2 extract stable phases ',&
!                  nlinesep,trim(phaseline(nlinesep))
!             do kkk=1,nlinesep
!                write(*,'(a,i3,3i5)')'smp2b color: ',nlinesep,&
!                     kkk,phamfu(1,kkk),phamfu(2,kkk)
!             enddo
! This destroyed phaseline for map with tie-lines in plasne ...
             if(isopleth) &
                  write(phaseline(nlinesep)(1:3),'(i3)')phamfu(2,nlinesep)
!             write(*,117)'smp2b phaseline: ',nlinesep,phamfu(2,nlinesep),&
!                  trim(phaseline(nlinesep))
117          format(a,2i5,2x,a)
! the segmentation fault was that linzero not always allocated ....
             if(allocated(linzero)) linzero=0
          endif
! no wildcards allowed on this axis
          statevar=pltax(notanp)
!          write(*,'(a,a,3i4,1pe12.4)')'In ocplot2, segmentation fault 4Ay1: ',&
!               trim(statevar),nr,nv,notanp,curceq%tpval(1)
          if(skipdotder) then
! If skipdotder and this is a dot derivative return error code to skip value
!             write(*,*)'smp2B skipping point A: ',trim(statevar),nv
!             skipdotder=.FALSE.; special_circumstances=1
             special_circumstances=1
          endif
          if(statevar(1:4).eq.'PFL ' .or. statevar(1:4).eq.'PFS ') then
             call meq_get_state_varorfun_value('NPM(LIQUID) ',&
                  value,encoded1,curceq)
             value=npflval*value
             npflval=value
             if(statevar(1:4).eq.'PFS ') value=one-value
          else
             call meq_get_state_varorfun_value(statevar,value,encoded1,curceq)
!          write(*,*)'SMP axis variable 1: ',trim(encoded1),value
             if(gx%bmperr.ne.0) then
! this error should not prevent plotting the other points FIRST SKIPPING
                write(*,212)'SMP skipping a point 1, error evaluating: ',&
                     statevar(1:10),curceq%tpval(1),nv,nr
212             format(a,a,f10.2,2i5)
! buperr resets putfun error 
                gx%bmperr=0; buperr=0
                nv=nv-1; goto 215
             endif
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
!          write(*,*)'In ocplot2: wildcard, selectph and novalues:',&
!               wildcard,selectph,novalues
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
             else
!                write(*,*)'SMP2B wildcard value 1: ',nr,trim(statevar)
!                write(*,*)'In ocplot2, segmentation fault after 4C1: ',&
!                     trim(statevar),nooftup()
! segmentation fault is inside this call for map11.OCM
! probably because new composition set created
!                write(*,*)'SMP2B get_many ',trim(statevar),nzp,selectph,hashtag
! nzp is dimentsion of yyy, np is number of values
!                write(*,*)'SMP2B calling get_many_svar: ',trim(statevar)
                call get_many_svar(statevar,yyy,nzp,np,encoded2,curceq)
!                write(*,*)'In ocplot2, segmentation fault search: '
                if(gx%bmperr.ne.0) then
                   write(*,*)'Error return from "get_many_svar',gx%bmperr
                   goto 1000
                endif
! problem that part of encoded2 desctroyed in late calls, it is OK here
!                write(*,737)len_trim(encoded2),trim(encoded2)
737             format('smp2b: debug encoded2 ',i5/a/)
!                write(*,738)(yyy(qp),qp=1,np)
738             format('SMP mm: ',(10F7.3))
! compiling without -finit-local-zero gives a segmentation fault here
! running the MAP11 macro
                qp=np
!                write(*,*)'SMP2B wildcard value 2: ',nr,trim(statevar)
!                write(*,223)'SMP2B Values: ',np,(yyy(i),i=1,np)
!                if(selectph) then
!                   write(*,*)'SMP2B: number of values: ',trim(selphase),np,nv
!223                format(a,i3,8F8.4)
!                endif
                nix=np
! we must allocate the array to indicate which values that should e plotted
                if(allocated(ixpos)) deallocate(ixpos)
                allocate(ixpos(nix),stat=errall)
                if(errall.ne.0) then
                   write(*,*)'SMP2B Allocation error 8: ',errall
                   gx%bmperr=4370; goto 1000
                endif
! This quite complicated IF is to handle the case when the wildcard is
! is a phase or component/constituent
!                write(*,*)'SMP stvarix: ',trim(statevar),selectph
                if(statevar(1:2).EQ.'MU' .or. &
                     statevar(1:2).EQ.'AC' .or. statevar(1:4).EQ.'LNAC' .or. &
                     selectph.and.(&
                     statevar(1:2).EQ.'N(' .or. statevar(1:2).EQ.'B(' .or. &
                     statevar(1:2).EQ.'X(' .or. statevar(1:2).EQ.'X%' .or. &
                     statevar(1:2).EQ.'W(' .or. statevar(1:2).EQ.'W%')) then
! for the state variables MU(*), AC(*), LNAC(*), N(*), X(*), W(*) the *
! means component, not phase, set selectph=.FALSE.
!                   write(*,*)'SMP Wildcard means component! ',trim(statevar),np
! not calling stvarix, all values should be included
! Allocation should be number of components
                   ixpos=1
                elseif(.not.selectph .and. .not.hashtag) then
!                elseif(selectph) then
! this routine supress values for phases that are not relevant ...
! IT SHOULD NOT BE USED FOR CASES LIKE Y(GAS,*)
!                   write(*,737)len_trim(encoded2),trim(encoded2)
! allocation should be number of phases
! if 
!                   write(*,'(a,a/a,2i5)')'SMP2B stvarix1: ',trim(statevar),&
!                        trim(encoded2),nlinesep,nix
                   call stvarix(statevar,phaseline(nlinesep),&
                        encoded2,nix,ixpos)
                   if(gx%bmperr.ne.0) then
                      write(*,*)'SMP2B yaxis error: "',trim(statevar),'"'
                      goto 1000
                   endif
                else
! here we may have hashtag TRUE
!                   write(*,*)'SMP2B Hashtag: ',hashtag
! We should supress values for suspended phases !!!
                   if(hashtag) then
                      call hashtag_susphix(statevar,phaseline(nlinesep),&
                           encoded2,nix,ixpos,curceq)
                      if(gx%bmperr.ne.0) then
                         write(*,*)'SMP2B failed handle hashtag'
                         goto 1000
                      endif
                   else
                      ixpos=1
                   endif
                endif
             endif
!             if(hashtag) then
! all values are used
!                ixpos=1
!             endif
!             write(*,*)'On ocplot2, segmentation fault 4D1'
!             write(*,213)trim(encoded2),np,(yyy(ic),ic=1,np)
213          format('WILDCARD: ',a,i3/6(1pe12.4))
!             write(*,214)np,(ixpos(ic),ic=1,np)
214          format('ixpos: ',12i3)
!             write(*,16)'val: ',kp,nr,gx%bmperr,(yyy(i),i=1,np)
16           format(a,2i3,i5,/6(1pe11.3/))
             anpmin=1.0D20
             anpmax=-1.0D20
             lcolor=0
!             write(*,'(a,90i3)')'SMP ixpos: ',(ixpos(jj),jj=1,np)
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
! now we have used stvarix to identify the relevant phases
                   if(yyy(jj).eq.zero) then
                      if(ixpos(jj).eq.0) then
! for STEP calculations try to make the ending a property at zero
! for MAP calculations just ignore the point ... also for STEP ...
                         anp(jj,nv)=rnone
                      else
                         anp(jj,nv)=zero
                      endif
                   else
! Hm, jumps from zero to finite values in step1, fig 3 plotting w(phase,cr) ..
!                      if(nv.gt.1 .and. anp(jj,nv-1).eq.rnone) then
!                         anp(jj,nv-1)=zero
!                      endif
                      anp(jj,nv)=yyy(jj)
                   endif
! difficult ...
                   if(yyy(jj).ne.rnone .and. ixpos(jj).ne.0) then
! ths is the trick to supress lines for phases that are never stable
                      nonzero(jj)=1
                      linzero(jj)=1
! save the first column with nonzero for use with invariants
                      if(ikol.eq.0) ikol=jj
                      if(anp(jj,nv).gt.anpmax) anpmax=anp(jj,nv)
                      if(anp(jj,nv).ne.rnone .and. &
                           anp(jj,nv).lt.anpmin) anpmin=anp(jj,nv)
! extract state variable jj used for table headings and key
                      if(.not.allocated(lid)) then
                         allocate(lid(np+5),stat=errall)
                         if(errall.ne.0) then
                            write(*,*)'SMP2B Allocation error 9: ',errall
                            gx%bmperr=4370; goto 1000
                         endif
                      endif
                      if(.not.isopleth) then
! if not isopleth save variable symbol in lid for headings (KEY)
! getext( , ,2, , , ) returns next text item up to a space
                         call getext(encoded2,lcolor,2,encoded1,'x',lhpos)
                         lid(jj)=encoded1
                         kk=len_trim(encoded1)
                         if(kk.gt.len(lid(jj))) then
                            lid(jj)(7:)='..'//encoded1(kk-6:kk)
                         endif
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
! segmentation fault if I plot Cp after shifting  to a NEW MAPTOP record
!             write(*,'(a,F8.2,a,a)')'SMP meq_get_state_varofun 7: T=',&
!                  curceq%tpval(1),' axis: ',trim(statevar)
! encoded1 not set correctly for dot derivative !!
!             encoded1='dummy'
! there is a segmentation fault in this call
! segmentation fault if maptop changed (several map/step commands)
             if(skipdotder) then
! return error code if calculating a derivative at first point of line
!                write(*,*)'smp2B skipping point B: ',trim(statevar),nv
                skipdotder=.FALSE.; special_circumstances=1
             endif
! special Scheil simulation 
             if(statevar(1:4).eq.'PFL ' .or. statevar(1:4).eq.'PFS ') then
                call meq_get_state_varorfun_value('NPM(LIQUID) ',&
                     value,encoded1,curceq)
                value=npflval*value
                npflval=value
                if(statevar(1:4).eq.'PFS ') value=one-value
                encoded1=statevar
             else
! end special Scheil, evaluate function normally
                call meq_get_state_varorfun_value(statevar,value,&
                     encoded1,curceq)
! encoded1 here is wrong?? not Cp when it should be, also when no error
!             write(*,*)'SMP axis value 7: ',trim(encoded1),value
                if(gx%bmperr.ne.0) then
! SECOND Skipping
                   if(gx%bmperr.ne.4373) &
                       write(*,212)'SMP Skipping a point 2,error evaluating: ',&
                        statevar(1:10),curceq%tpval(1),nv,nr
                   nv=nv-1; goto 215
                endif
             endif
! save to use in lid if not allocated
             funsym=encoded1
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
!                write(*,*)'We have found last equilibria along the line: ',nr
                last=.TRUE.
! skipdotder=TRUE means skip dot derivatives at first equilibrium in next line
                skipdotder=.TRUE.
             endif
          endif
!>>>>>>>>>>>>>>>>>>>>>>>>>> starting a line
!          write(*,*)'Next equilibrium: ',nr,nv,xax(nv)
!          read(*,17)ch1
17           format(a)
       enddo plot1
220    continue
! finished one line
!       write(*,*)'SMP2B at 220: nr and nv: ',nr,nv
       invariant_lines: if(nax.gt.1) then
!---------------------------------------------------------------
!------------------ special for invariant lines ?? and others
!---------------------------------------------------------------
! for phase diagram always move to the new line 
          map1: if(nlinesep.ge.1) then
             newsep: if(linesep(nlinesep).lt.nv) then
! we should never have several linesep for the same value of nv!
                nlinesep=nlinesep+1
                linesep(nlinesep)=nv
! names of phases on new line
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
! This is a check for node with 2 stoichiometric phases, if so skip first line
!                   if(invar%artxe.eq.1) then
!                      write(*,*)'Found node with 2 stoichiomeric phases'
!                   endif
                   curceq=>results%savedceq(invar%savednodeceq)
!-------------------
! get the names of stable phases from the node equilibrium record
                   kk=1
                   stoichfix=0
                   extractphnames: do jj=1,noph()
                      do ic=1,noofcs(jj)
! value is amount of phase?
                         k3=test_phase_status(jj,ic,value,curceq)
                         if(k3.gt.0) then
! this phase is stable or fix
!                            write(*,*)'SMP addto phaseline 2: ',trim(dummy),&
!                                 nlinesep
                            call get_phase_name(jj,ic,dummy)
                            if(kk.lt.100) phaseline(nlinesep)(kk:)=dummy
                            kk=len_trim(phaseline(nlinesep))+2
                            stoichfix=stoichfix+1
                         endif
                      enddo
                   enddo extractphnames
! invariant         write(*,*)'SMP2 phases: ',trim(phaseline(nlinesep)),nlinesep
                   if(stoichfix.gt.3) then
                      write(*,*)'SMP2B too many stable phases at invariant',&
                           stoichfix
                      stoichfix=3
                   endif
!-------------------
! axis without wildcard
                   statevar=pltax(notanp)
                   if(skipdotder) then
! skip calculating a derivative if this is first point of a region
!                      write(*,*)'smp2B skipping point C: ',trim(statevar),nv
                      skipdotder=.FALSE.; special_circumstances=1
                   endif
                   call meq_get_state_varorfun_value(statevar,value,&
                        encoded1,curceq)
!                   write(*,*)'SMP axis variable 3: ',encoded1(1:3),value
                   if(gx%bmperr.ne.0) then
! THIRD skipping
                      write(*,212)'SMP skipping a point 3, error evaluating ',&
                           statevar,curceq%tpval(1),nv,0
                      goto 222
                   endif
! save symbol name if lid not allocated
                   funsym=encoded1
! This is tielines inplane, normally 3 lines to generate
! but when 2 stoichiometric phass with same composition one is not set stable
!                   nv=nv+3
! NOTE if not wildcard nv is decremented in next "else" statement
                   nv=nv+stoichfix
                   if(nv.ge.maxval) then
                      write(*,*)'Too many points to plot 2',maxval
                      goto 600
                   endif
                   do invlines=0,stoichfix-1
                      xax(nv-invlines)=value
                   enddo
!                   write(*,335)'New line: ',nlinesep,nv,linesep(nlinesep),&
!                        statevar(1:5),value
335                format(a,3i4,' <',a,'> ',3(1pe14.6))
! axis with possible wildcard
                   statevar=pltax(anpax)
!                   write(*,*)'In ocplot2, segmentation fault 4H'
                   wild2: if(wildcard) then
! this cannot be a state variable derivative
!                    write(*,*)'Getting a wildcard value 2: ',nr,statevar(1:20)
                      call get_many_svar(statevar,yyy,nzp,np,encoded2,curceq)
                      if(gx%bmperr.ne.0) goto 1000
! we have to handle axis values that are zero what is np here???
                      nix=np
! to supress suspended phases
!                      write(*,'(a,a/a,2i5)')'SMP2B stvarix2: ',trim(statevar),&
!                           trim(encoded2),nlinesep,nix
                      call stvarix(statevar,phaseline(nlinesep),&
                           encoded2,nix,ixpos)
                      if(gx%bmperr.ne.0) goto 1000
! save one non-zero value per line, 3 lines
                      ic=0
                      do jj=1,np
! np is the number of values retrieved by get_many_svar
! only those with nonzero values in ixpos should be used, one per line.
                         if(ixpos(jj).ne.0) then
                            anp(ikol,nv-ic)=yyy(jj)
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
! if no wildcard extract the phase with zero amount

                      nv=nv-stoichfix
                      goto 225
                   endif wild2
222                continue
!                else
!                   write(*,*)'SMP no else link ...'
                endif inv
             endif newsep
          endif map1
! jump here if no wildcard
225       continue
       endif invariant_lines
!       
       if(invnode.ne.0) then
!          write(*,'(a,3i4,2(1pe12.4))')'ocplot2 Invariant isopleth node: ',&
!               invnode,ninv,nrett,xax(nrett),anp(1,nrett)
          if(ninv.eq.0) then
! the first invariant isopleth 
             ninv=ninv+1
             xyinv(1,1)=xax(nrett); xyinv(2,1)=anp(1,nrett)
             xyinv(3,1)=xax(nrett); xyinv(4,1)=anp(1,nrett)
          else
! check if anp value same as already saved invariant
             cinv: do kk=1,ninv
                if(abs(anp(1,nrett)-xyinv(2,kk)).lt.1.0D-4) then
! same ... check if x values lesser than xyinv(1,kk) or greater than xyinv(3,kk)
                   if(xax(nrett).lt.xyinv(1,kk)) xyinv(1,kk)=xax(nrett)
                   if(xax(nrett).gt.xyinv(3,kk)) xyinv(3,kk)=xax(nrett)
                   goto 227
                endif
             enddo cinv
! this is a new invariant, do not accept zero values
             if(abs(xax(nrett)).gt.1.0D-6.and.abs(anp(1,nrett)).gt.1.0D-6) then
                ninv=ninv+1
                xyinv(1,ninv)=xax(nrett); xyinv(2,ninv)=anp(1,nrett)
                xyinv(3,ninv)=xax(nrett); xyinv(4,ninv)=anp(1,nrett)
!                write(*,*)'New invariant:',ninv,xyinv(1,ninv),xyinv(2,ninv)
!             else
!                write(*,*)'Invariant with zero values ignored'
             endif
227          continue
          endif
       endif
!---- take next node along the same line
! Then jump back to label 100 and plot other lines ... a bit stupid ...
230    continue
!       write(*,*)'SMP2B at 230: nr and nv: ',nr,nv
       kk=seqx
       if(associated(mapline%end)) then
          seqx=mapline%end%seqx
       else
          seqx=0
       endif
240    continue
!       write(*,'(a,5i5,l2)')'ocplot2 next node: ',seqx,nlinesep,&
!            linesep(nlinesep),nv,line,scheilorder
       if(seqx.eq.0) then
          if(nlinesep.gt.0) then
             if(linesep(nlinesep).lt.nv) then
! we should never have several linesep for the same value of nv!
                nlinesep=nlinesep+1
                linesep(nlinesep)=nv
!                write(*,*)'adding empty line 2',nlinesep,linesep(nlinesep)
             endif
          endif
! Hm, this was not designed for multicomponent isopleths ....
          if(line.eq.2) then
!             write(*,*)'ocplot2 jump to label 500',line
             goto 500
          endif
          if(scheilorder) then
! The mapnodes must be followed in numeric order
! ane they have just one line each.
              line=1
!             write(*,*)'SMP2B scheilorder',localtop%seqx,&
!                  localtop%next%seqx,localtop%previous%seqx
             localtop=>localtop%previous
             mapline=>localtop%linehead(1)
!             write(*,*)'ocplot2 newline ',nlinesep,linesep(nlinesep)
             goto 110
          endif
          line=2
! jump back to label 100 for next line
!          write(*,*)'ocplot2 jump to label 100 for line 2'
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
          invnode=0
          if(btest(mapnode%status,MAPINVARIANT)) then
             invnode=size(mapnode%linehead)
!             write(*,*)'ocplot2 invariant node 2',invnode
          endif
! loop through all mapnodes
250       continue
!          write(*,*)'ocplot2, at label 250: ',mapnode%seqx,seqx
!          if(mapnode%seqx.eq.seqx) then
! If just for STEP then check number of axis for calculation
          if(graphopt%noofcalcax.eq.1 .and. mapnode%seqx.eq.seqx) then
! >>> this is just for step, for map one must find line connected
             do haha=1,size(mapnode%linehead)
                mapline=>mapnode%linehead(haha)
                if(mapline%done.eq.0) then
                   if(.not.btest(mapline%status,EXCLUDEDLINE)) then
!                      write(*,*)'ocplot2 jump to label 110',&
!                           seqx,mapline%number_of_equilibria
                      goto 110
                   else
                      lines_excluded=lines_excluded+1
                   endif
                endif
             enddo
          endif
!          write(*,*)'ocplot2 associated?  ',associated(mapnode,localtop),&
!               mapnode%seqx,seqx
          if(.not.associated(mapnode,localtop)) then
             mapnode=>mapnode%next
             invnode=0
             if(btest(mapnode%status,MAPINVARIANT)) then
                invnode=size(mapnode%linehead)
! Does all invariant nodes have the same number of stable phases YES!
! But they can have different number of line exits
!                write(*,*)'ocplot2 invariant node 3',invnode,mapnode%seqx
             endif
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
!       write(*,*)'ocplot2 at 500 ',seqx
!       write(*,*)'In ocplot2, looking for segmentation fault 5'
!       mapnode=>localtop%next
! when we have several plots localtop is the important one!!
       mapnode=>localtop
       invnode=0
       if(btest(mapnode%status,MAPINVARIANT)) then
          invnode=size(mapnode%linehead)
!          if(invnode.ne.mapnode%lines) write(*,*)'SMP2B check invnodes 1!'
!          write(*,*)'ocplot2 invariant node 4',invnode
       endif
! Check for unplotted lines
       anymoretoplot: do while(.TRUE.)
!          write(*,*)'>>>>>Checking unplotted lines at node: ',mapnode%seqx
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
! skip line if EXCLUDEDLINE set
                if(btest(mapline%status,EXCLUDEDLINE)) then
                   cycle jjline
                endif
!                write(*,*)'SMP2B jump to 110 for line ',jj,&
!                     ' in mapnode ',mapnode%seqx
! %done=-1 means already plotted ...
!                mapnode%linehead(jj)%done=-1
                goto 110
             endif
          enddo jjline
          mapnode=>mapnode%next
          invnode=0
          if(btest(mapnode%status,MAPINVARIANT)) then
             invnode=size(mapnode%linehead)
!             if(invnode.ne.mapnode%lines) write(*,*)'SMP2B check invnodes 2!'
!             write(*,*)'ocplot2 invariant node 5',invnode
          endif
          if(associated(mapnode,localtop)) exit anymoretoplot
       enddo anymoretoplot
!--------------------------------------------
! end extracting data
600    continue
       overflow=.FALSE.
! but we may have another maptop !!
       if(associated(localtop%plotlink) .and. .not.scheilorder) then
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
!       write(*,651)'SMP wildcard, np=, qp, nnp, ic: ',wildcard,np,qp,nnp,ic
651    format(a,l2,10i5)
!------------------------------------------ begin loop 650
650    ic=ic+1
660       continue
          if(ic.gt.nnp) goto 690
! nonzero(ic) is a column with nonzero value to plot ... redundant??
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
                allocate(lid(np+5),stat=errall)
                if(errall.ne.0) then
                   write(*,*)'SMP2B Allocation error 10: ',errall
                   gx%bmperr=4370; goto 1000
                endif
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
!       write(*,651)'SMP at 690:                     ',wildcard,np,qp,nnp,ic
       nrv=nv
       np=nnp
!       goto 800
!============================================ generate gnuplot file
800 continue
! are there any isopleth invariants?
       if(ninv.gt.0) then
          do kk=1,ninv
! nlinesep is last line with data, linesep(nlinsesp) is index of last data line
!             write(*,'("Isoinv: ",3i4,4(1pe12.4))')nlinesep,linesep(nlinesep),&
!                  kk,(xyinv(jj,kk),jj=1,4)
! add these to lines to be plotted
             kkk=nlinesep
             if(len_trim(phaseline(kkk)).gt.0) then
                write(*,*)'smp2b phaseline: "',trim(phaseline(kkk)),'"',kkk
             endif
             phaseline(kkk)='100 invariant equilibrium'
             phaseline(kkk+1)=' '
!             write(*,*)'smp2b isoinv: "',trim(phaseline(kkk)),'"',kkk
             jj=linesep(kkk)
             nlinesep=nlinesep+1
! the line nlinesep contain  2 points, beginning and end of invariant line
             linesep(nlinesep)=jj+2
             xax(jj+1)=xyinv(1,kk); anp(1,jj+1)=xyinv(2,kk)
             xax(jj+2)=xyinv(3,kk); anp(1,jj+2)=xyinv(4,kk)
             nrv=nrv+2
          enddo
       endif
! add the invariant lines to be plotted
! two points for each invariants (X and Y)
!
!    write(*,808)np,nv,nlinesep,maxanp,maxval
!    write(*,'(a,(16i4))')'SMP pp: ',(linesep(kk),kk=1,nlinesep)
808 format('plot data used: ',3i7,' out of ',2i7)
    if(np.eq.0) then
       write(kou,*)'No data to plot',np
       gx%bmperr=4248
       goto 1000
    endif
!------------------------------------------------------------
    if(.not.allocated(lid)) then
! lid is "LineIdenifier and changes color for each line with different meaning
!       if(np.ge.1) then
! lid should always be allocated if np>1, but ... one never knows 
!       write(*,*)'ocplot2 allocate lid 4: ',np,nlinesep
!       if(scheilorder) then
! for Scheil simulations phaseline(1..nlinesep) are phase names
!          allocate(lid(nlinesep-1),stat=errall)
!          do i=1,nlinesep-1
!             jj=len_trim(phaseline(i))
!             if(jj.le.32) then
!                lid(i)=phaseline(i)
!             else
! Too long list, replace the middle by ...
!                lid(i)=phaseline(i)(1:22)//'...'//phaseline(i)(jj-6:jj)
!             endif
!             write(*,*)'ocplot2 phaseline: ',trim(lid(i))
!          enddo
! Wow, set np=nlinesep-1 to have separate colors of Scheil lines
!          write(*,*)'ocplot2 change np to nphaseline: '
!          np=nlinesep-1
!       else
! normally np=1 if we come here, plotting a single value
!          write(*,*)'Allocating lid: ',trim(encoded1),':',trim(funsym),np
       allocate(lid(np),stat=errall)
       if(errall.ne.0) then
          write(*,*)'SMP2B Allocation error 11: ',errall
          gx%bmperr=4370; goto 1000
       endif
       do i=1,np
          lid(i)=funsym
       enddo
    endif
!    endif
!------------------------------------------------------------
2000 continue
!    write(*,*)'We are at 2000 '
!----------------------------------------------------------------------
!
    call get_plot_conditions(encoded1,maptop%number_ofaxis,axarr,ceq)
!
! option to create a CSV table
    if(btest(graphopt%status,GRCSVTABLE)) then
       call list_csv(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,anp,lid,&
            phaseline,title,filename,version,encoded1)
    else
! NOW pltax should be the the axis labels if set manually
       if(graphopt%labeldefaults(2).ne.0) pltax(1)=graphopt%plotlabels(2)
       if(graphopt%labeldefaults(3).ne.0) pltax(2)=graphopt%plotlabels(3)
! write(*,*)' >>>>>>>>**>>> plot file: ',trim(filename)
       call ocplot2B(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,anp,lid,&
            phaseline,title,filename,graphopt,version,encoded1,fixphasecolor)
!         title,filename,graphopt,pform,version,encoded1)
!    goto 900
    endif
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

!\addtotable subroutine ocplot2B
!\begin{verbatim} %-
  subroutine ocplot2B(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,anp,lid,&
       phaseline,title,filename,graphopt,version,conditions,fixphasecolor)
! called from icplot2 to generate the GNUPLOT file after extracting data
! np is number of columns (separate lines), if 1 no labelkey
! nrv is number of values to plot
! nlinesep is the number of separate lines (index to linesep)
! linesep is the row index when a line to be plotted finishes (1..nlinesep)
! pltax are axis labels
! xax array of values for single valued axis (T or mu etc)
! anpax=2 if axis with single value is column 2 and (multiple) values in
!         columns 3 and higher
! anp array of values for axis with multiple values (can be single values also)
! lid array with a label for the different lines
! title Title of the plot
! filename GNUPLOT file name, (also used for pdf/ps/gif file)
! graphopt is graphical option record
! conditions is a character with the conditions for the diagram
! fixphasecolor is used for isopleths to know how many columns are needed
    implicit none
    integer np,anpax,nlinesep,fixphasecolor
    integer ndx,nrv,linesep(*),anpdim,npx
    character pltax(*)*(*),filename*(*),lid(*)*(*),title*(*)
    character conditions*(*),version*(*),phaseline(*)*(*)
    type(graphics_options) :: graphopt
    double precision xax(*),anp(anpdim,*)
    double precision scale1,scalem
    type(graphics_textlabel), pointer :: textlabel
!\end{verbatim}
!----------------------------------------------------------------------
! internal
    integer ii,jj,kk,lcolor,appfil,nnv,ic,repeat,ksep,nv,k3,kkk,nofapl,iz
    integer, parameter :: mofapl=100
    integer, parameter :: maxmultiplotlines=100
! ltf1 is a LineTypeoFfset for current plot when appending a plot, 0 default
! linewp is plotting with points along the line
    integer appfiletyp,lz,ltf1,linewp
! appending an multiplot in ocplot2B
    integer multibuffline
    character multibuffer(maxmultiplotlines)*128
    logical appendmultiplot
! other things ...
    character pfc*128,pfh*128,backslash*2,appline*128,inline*8,colord*5
    character applines(mofapl)*128,gnuplotline*256,labelkey*64,rotate*16
    character labelfont*32,linespoints*12,tablename*16,year*16,hour*16
    logical isoplethplot
! write the gnuplot command file with data appended
!
!    write(*,10)'in ocplot2B: ',np,anpax,nrv,nlinesep,trim(title),&
!         (linesep(kk),kk=1,nlinesep)
!    write(*,*)'smp2b isoplethplot 2: ',btest(graphopt%status,GRISOPLETH)
10  format(a,4i5,a/(15i4))
! graphopt%specialdiagram=2 is a Scheil plot phase amount/vs T
!    if(graphopt%specialdiagram.eq.2) then
!       write(*,*)'In ocplot2B scheil: ',graphopt%specialdiagram
! just check all is OK
!       do kk=1,nlinesep-1
!          write(*,*)'phaseline: ',kk,': ',trim(phaseline(kk))
!       enddo
!    endif
    multibuffline=0
    nofapl=0
    appendmultiplot=.FALSE.
    ltf1=0
    if(graphopt%appendfile(1:1).ne.' ') ltf1=10
    if(index(filename,'.plt ').le.0) then 
       kk=len_trim(filename)
       pfc=filename(1:kk)//'.'//'plt '
       kkk=kk+4
    else
       pfc=filename
       kk=index(pfc,'.')-1
       kkk=len_trim(filename)+1
    endif
    open(21,file=pfc,access='sequential',status='unknown')
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
! terminal 1 is screen without any output file, add PDF as comment
       write(21,841)trim(graphopt%gnuterminal(graphopt%gnutermsel)),&
            trim(graphopt%gnuterminal(3))
841    format('set terminal ',a/'#set terminal ',a/'#set output "ocgnu.pdf"')
!841    format('set terminal ',a)
    endif
! for isopleths we will generate one column for each phase even if there
! are just a single value in anpax.  To handle this we need different
! values of "np" for these columns
    isoplethplot=btest(graphopt%status,GRISOPLETH)
    npx=np
    if(isoplethplot) then
! This is the number of columns with phases in the isopleth incl. invariant
       npx=fixphasecolor
    endif
! this part is independent of which axis is a single value
!------------------ some GNUPLOT colors:
! colors are black: #000000, red: #ff000, web-green: #00C000, web-blue: #0080FF
! dark-yellow: #C8C800, royal-blue: #4169E1, steel-blue #306080,
! gray: #C0C0C0, cyan: #00FFFF, orchid4: #804080, chartreuse: 7CFF40
! if just one line set key off for that line.
    if(graphopt%specialdiagram.eq.2) then
! for Scheil
       labelkey=' on font "Arial,12" '
    elseif(npx.eq.1 .and. graphopt%appendfile(1:1).eq.' ') then
       labelkey=' off'
    else
       labelkey=graphopt%labelkey
    endif
    call date_and_time(year,hour)
!    write(*,*)'"',year,'"  "',hour,'"'
    tablename='OCT'//year(3:8)//hour(1:6)
!    write(*,*)'Plot heading 2? ',btest(graphopt%status,GRNOTITLE)
    call replace_uwh(conditions)
    if(btest(graphopt%status,GRNOTITLE)) then
       write(21,858)trim(title),trim(conditions),trim(graphopt%font)
    else
       write(21,859)trim(title),trim(conditions),trim(graphopt%font)
    endif
858 format('#set title "',a,' \n #',a,'" font "',a,',10" ')
859 format('set title "',a,' \n ',a,'" font "',a,',10" ')
    lz=graphopt%linetype
! replace _ and & in axis texts by "\_" and "\&"
!    write(*,*)'Replacing _ and &: ',trim(pltax(2))
    call replace_uwh(pltax(1))
    call replace_uwh(pltax(2))
    if(isoplethplot) then
       write(21,8601)graphopt%xsize,graphopt%ysize,&
            trim(pltax(1)),trim(pltax(2)),trim(labelkey),&
            ltf1+1,ltf1+2,ltf1+3,ltf1+4,ltf1+5,&
            ltf1+6,ltf1+7,ltf1+8,ltf1+9,ltf1+10
8601   format('set origin 0.0, 0.0 '/&
            'set size ',F8.4', ',F8.4/&
            'set xlabel "',a,'"'/'set ylabel "',a,'"'/&
! Help with stackoverflow to fix nice logo independent of plot size!
!        'set label "~O{.0  C}" at graph -0.1, -0.1 font "Garamond Bold,20"'/&
        'set label "~O{.0  C}" at screen 0.02, 0.03 font "Garamond Bold,20"'/&
            'set key ',a/&
            'set linetype ',i2,' lc rgb "#000000" lw 2 pt 10'/&
            'set linetype ',i2,' lc rgb "#4169E1" lw 2 pt 6'/&
            'set linetype ',i2,' lc rgb "#00C000" lw 2 pt 3'/&
            'set linetype ',i2,' lc rgb "#FF0000" lw 2 pt 2'/&
            'set linetype ',i2,' lc rgb "#FF00FF" lw 2 pt 4'/&
            'set linetype ',i2,' lc rgb "#C8C800" lw 2 pt 5'/&
            'set linetype ',i2,' lc rgb "#C0C0C0" lw 2 pt 7'/&
            'set linetype ',i2,' lc rgb "#00FFFF" lw 2 pt 8'/&
            'set linetype ',i2,' lc rgb "#804080" lw 2 pt 9'/&
            'set linetype ',i2,' lc rgb "#7CFF40" lw 2 pt 1'/&
            '# for invariants, orange'/&
            'set linetype 100 lc rgb "#FFC000" lw 3 pt 1')
       else
          write(21,860)graphopt%xsize,graphopt%ysize,&
               trim(pltax(1)),trim(pltax(2)),trim(labelkey),&
               ltf1+1,lz,ltf1+2,lz,ltf1+3,lz,ltf1+4,lz,ltf1+5,lz,&
               ltf1+6,lz,ltf1+7,lz,ltf1+8,lz,ltf1+9,lz,ltf1+10,lz
860       format('set origin 0.0, 0.0 '/&
               'set size ',F8.4', ',F8.4/&
               'set xlabel "',a,'"'/'set ylabel "',a,'"'/&
! Help with stackoverflow to fix nice logo independent of plot size!
!          'set label "~O{.0  C}" at graph -0.1, -0.1 font "Garamond Bold,20"'/&
        'set label "~O{.0  C}" at screen 0.02, 0.03 font "Garamond Bold,20"'/&
               'set key ',a/&
               'set style line ',i2,' lt ',i2,' lc rgb "#000000" lw 2 pt 10'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#4169E1" lw 2 pt 6'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#00C000" lw 2 pt 3'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#FF0000" lw 2 pt 2'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#FF00FF" lw 2 pt 4'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#C8C800" lw 2 pt 5'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#C0C0C0" lw 2 pt 7'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#00FFFF" lw 2 pt 8'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#804080" lw 2 pt 9'/&
               'set style line ',i2,' lt ',i2,' lc rgb "#7CFF40" lw 2 pt 1'/&
               '# for invariants, faded read'/&
               'set style line 100 lt 1 lc rgb "#FF3333" lw 3 pt 1')
       endif
! add some useful things for maniplulation of graph
    write(21,8000)
8000 format(/'# Some useful GNUPLOT commands for editing the figure'/&
          '# This is a dashed line (on pdf/wxt):'/&
          '# set style line 15 lt 0 lc rgb "#C8C800" lw 2 pt 2'//&
          '# set pointsize 0.6'/&
          '# set label "text" at 0.5, 0.5 rotate by 60 font "arial,12"'/&
          '# set xrange [0.5 : 0.7] '/&
          '# Adding manually a line and keep scaling:'/&
          '# set arrow x0, y0 to x1,y1 nohead linestyle 1'/&
          '# set multiplot'/&
          '# set xrange [] writeback'/&
          '#  ... plot someting'/&
          '# set xrange restore'/&
          '#  ... plot more using same axis scaling '/&
          '# unset multiplot'/)
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
!          write(labelfont,178)int(10*textlabel%textfontscale)
!178       format(' font "Sans,',i2,'" ')
!          write(*,1178)trim(graphopt%font),int(10*textlabel%textfontscale)
!1178       format(' SMP2B font "',a,',',i2,'" ')
          write(labelfont,178)trim(graphopt%font),&
               int(10*textlabel%textfontscale)
178       format(' font "',a,',',i2,'" ')
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
    appfildata: if(graphopt%appendfile(1:1).eq.' ') then
       appfil=0
    else
       appfil=23
       write(kou,*)'Appending data from: ',trim(graphopt%appendfile)
       open(appfil,file=graphopt%appendfile,status='old',&
            access='sequential',err=1750)
!
       write(21,1719)trim(graphopt%appendfile)
1719   format(//78('#')/'# APPENDED from ',a)
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
! save lines between "set multiplot" and "unset multiplot" to a buffer
! appending a file which contains an already appended part
! copy all lines to a buffer
       if(appline(1:14).eq.'set multiplot ') then
!          write(*,*)'ocplot2B Found multiplot in appended file'
          if(multibuffline.gt.0) then
             write(*,*)'smp2B appending twice, will probably fail',multibuffline
          else
             multibuffline=1
!             write(*,*)'ocplot2B appending a "set multiplot"'
          endif
          appendmultiplot=.TRUE.
          goto 1710
       endif
       if(appline(1:16).eq.'unset multiplot ') then
!          write(*,*)'ocplot2B found "unset multiplot", saved ',&
!               multibuffline,' lines'
          appendmultiplot=.FALSE.
!          do iz=1,multibuffline-1
!             write(*,*)'ocplot2B: 'trim(multibuffer(iz))
!          enddo
! pause mouse ??
          exit appfildata
!          goto 1710
       endif
       if(appendmultiplot) then
          if(multibuffline.gt.maxmultiplotlines) then
             write(*,*)'Too many appendbuffer lines',multibuffline
          else
             multibuffer(multibuffline)=appline
             multibuffline=multibuffline+1
          endif
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
       write(*,*)'SMP appfiletyp: ',appfiletyp
! here we save the actual plot commands from the appendfile!!
       applines(1)=appline
!       write(*,*)'SMP appline1: ',trim(applines(1))
       ic=1
1730   continue
! if line ends with \ then read more
       ii=len_trim(appline)
! write(*,*)'There are more? ',appline(i:ii),ii,ic
       if(appline(ii:ii).eq.'\') then
!       if(appline(ii:ii).eq.' ') then
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
       nofapl=ic
       write(*,*)'Read applines header lines: ',nofapl
! debug output of saved plot command
!       nofapl=ic
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
       write(21,1719)'# end of append',multibuffline
!       write(*,1719)' end of append',multibuffline
       write(kou,*)'Error opening or reading the append file, skipping it'
       close(appfil)
       appfil=0
1770   continue
    endif appfildata
! end of appendfile special
!-----------------------------------------------
! text in lower left corner
    ii=len_trim(graphopt%lowerleftcorner)
    if(ii.gt.0) then
! in square diagram below figure
!       write(21,209)trim(graphopt%lowerleftcorner)
!209    format('set label "',a,'" at graph -0.10, -0.08 ')
       write(21,209)trim(graphopt%lowerleftcorner)
209    format('set label "',a,'" at graph -0.05, -0.05 ')
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
    if(graphopt%specialdiagram.eq.2) then
!============================================== start Scheil
! Handle Scheil diagram with color changes along a single line
! ONLY if one axis is PFL or PFS!!! NOT if one plot composition of a phase
!       write(*,*)'ocplot2B special for Scheil with PFL or PFS',nlinesep
! if there is an appended file we must set multiplot here ...
       if(appfil.gt.0) then
          write(*,*)'ocplot2B adding set multiplot'
          write(21,3828)
       endif
       backslash=',\'
       inline='plot "-"'
! All lines to plot, colord can be 2:3 or 3:2 depending what is on the axis
          colord=' 2:3 '
! no need to change, evidently xax and anp are already shifted ... suck
!       write(*,'(a,a,2F12.6)')'ocplot2B colord: ',colord,xax(1),anp(1,1)
       do kk=1,nlinesep-1
          call replace_uwh(phaseline(kk))
          write(21,2000)inline,colord,kk,trim(phaseline(kk)),backslash
2000      format(a,' using 'a,' with lines ls ',i2,' title "',a,'"',a)
          inline='""'
          if(kk.eq.nlinesep-2) backslash=' '
       enddo
! Then all data with an empty line and a line with a single  "e"       
! between each line.
! nlinesep is number of separate lines to plot (with different sets of phases)
! linesep(1..nlinesep) is number of lines of data for each line to plot
! nrv is total lines with lines with data to write
!       write(*,*)'ocplot2B linesep: ',(linesep(jj),jj=1,nlinesep)
       write(21,'(a)')'# Line   1, phases: '//trim(phaseline(1))
       jj=2
       ltw: do nv=1,nrv
!          write(*,'(3i4,1pe12.4)')'ocgnu2B data: ',jj,nv,linesep(jj),xax(nv)
          if(nv.eq.linesep(jj)) then
! a new line start
             if(jj.eq.nlinesep) exit ltw
             write(21,'(i4,2F12.6)')nv,xax(nv),anp(1,nv)
             write(21,2100)jj,trim(phaseline(jj))
2100         format(/'e'/'# Line ',i3,' phases ',a)
             jj=jj+1
          endif
          write(21,'(i4,2F12.6)')nv,xax(nv),anp(1,nv)
       enddo ltw
       write(21,2110)
2110   format(/'e'/)
!2110   format(/'e'/'pause mouse')
!       close(21)
! Finished writing plot file
!       stop 'test'
       goto 4000
!============================================== end Scheil
    endif
! now write all data once as a table ended with EOD
    write(21,3000)nrv,trim(tablename)
3000 format(//'# begin of data with lines',i7/'$',a,' << EOD')
!
! A digit before the first phase gives number of columns to plot
!    write(*,*)'smp2b: isopleth? ',isoplethplot,np,npx
    if(isoplethplot) read(phaseline(1),'(i3)')fixphasecolor
!    write(*,*)'SMP2B replace _ in keys: ',npx
    do jj=1,npx
! remove _ in keys
       call replace_uwh(lid(jj))
    enddo
    if(graphopt%tielines.gt.0) then
       write(*,*)'ocplot2 does not plot tielines,',&
            ' they are perpendicular to the potential axis'
    endif
! Plot grid?
    if(graphopt%setgrid.eq.1) write(21,777)
777 format('set grid')
! columnheaders used as keys
    if(isoplethplot) then
! This column headin is not set before
       lid(npx)='Invariant'
!       write(*,3100)'KEYS: ',trim(pltax(3-anpax)),(trim(lid(jj)),jj=1,npx)
       write(21,3100)'KEYS: ',trim(pltax(3-anpax)),(trim(lid(jj)),jj=1,npx)
!2900 format(a,i3,2x,10(a,2x))
    else
       write(21,3100)'KEYS: ',trim(pltax(3-anpax)),(trim(lid(jj)),jj=1,np)
3100   format(a,' ',100(a,' '))
    endif
    ksep=2
    write(21,*)'# First line: ',trim(phaseline(1))
!    write(*,*)'smp2b isoplethplot 3: ',isoplethplot,&
!         btest(graphopt%status,GRISOPLETH),fixphasecolor
    do nv=1,nrv
!---------------------------------------------------------------
! values written multiplied with graphopt%scalefact, 
! first value is single valued axis (can be X or Y axis) multiplied with scale1
! remaining values multiplied svalem
       write(21,'(i4,1pe16.6)',advance='no')nv,scale1*xax(nv)
! note that isopletpots have np=1
       do jj=1,np-1
! second and later columns represent Y axis
          if(anp(jj,nv).ne.rnone) then
             write(21,2821,advance='no')scalem*anp(jj,nv)
          else
             write(21,'(a)',advance='no')' NaN '
          endif
       enddo
       if(isoplethplot) then
! fixphasecolor 100 means invariant
          if(fixphasecolor.lt.100) then
! dummy column values up to fixphasecolor
             do jj=1,fixphasecolor-1
                write(21,'(" NaN ")',advance='no')
             enddo
! This column has real value, then maybe additional columns with dummy values
! note only anp(1,nv) has any value!!
             write(21,2821,advance='no')scalem*anp(1,nv)
             do jj=1,npx-fixphasecolor-1
                write(21,'(" NaN ")',advance='no')
             enddo
! The last colum is for invariants which are written seperatey below
             write(21,'(" NaN ")')
          else
! this is an invariant, last column has values.  There may be no invariant!
             do jj=1,npx-1
                write(21,'(" NaN ")',advance='no')
             enddo
             write(21,2821)scalem*anp(1,nv)
          endif
       elseif(anp(jj,nv).ne.rnone) then
          write(21,2821)scalem*anp(jj,nv)
       else
          write(21,'(a)')' NaN '
       endif
!2820   format(i4,1pe16.6)
2821   format(1pe16.6)
!2822   format(' NaN ')
!---------------------------------------------------------------
! here we shift to another line and color
       if(nv.eq.linesep(ksep)) then
! an empty line in the plot file means a MOVE to the next point.
          if(nv.lt.nrv) then
             write(21,3819)ksep-1,trim(phaseline(ksep-1)),trim(phaseline(ksep))
3819         format('# end of line ',i3,2x,a//'# new line ',a)
!             write(*,*)'SMP2B readfixcolor: ',trim(phaseline(ksep)),ksep
             if(isoplethplot) read(phaseline(ksep),'(i3)')fixphasecolor
          else
! try to avoid rubbish
             write(21,3821)ksep-1,trim(phaseline(ksep-1))
3821         format('# end of line ',i3,2x,a//)
          endif
! test of uninitiallized variable, ksep must not exceed nlinesep
          ksep=min(ksep+1,nlinesep)
       endif
    enddo
    write(21,3823)
3823 format('EOD'//)
    if(appfil.gt.0) then
! if there is an appendfile add set multiplot
       write(*,*)'ocpolt2B trying to include appfile ...'
! The "writeback" is important for uniform scaling of multiplots
! NOTE this is also used for Scheil above
       write(21,3828)
3828   format('set multiplot'/&
            'set xrange [] writeback'/'set yrange [] writeback')
    endif
! If no file appended the line types are (i-2)
! if a file appended then line types are (i-2+ltf1(=10))
! if anpax is axis with single value (1=x, 2=y)
!    if(isoplethplot) then
!       write(21,3800)trim(tablename)
!3800   format('plot $',a,' using 2:3:4 with lines lc variable notitle')
!    elseif(anpax.eq.1) then
! here we use npx for both isopleths and others!
    if(graphopt%linewp.le.1) then
       if(anpax.eq.1) then
! linewp=0 is dashed and =1 is line without point
          write(21,3900)npx+2,trim(tablename),ltf1
3900      format('plot for [i=3:',i2,'] $',a,' using i:2',&
               ' with lines ls (i-2+',i2,') title columnheader(i)') 
       else
          write(21,3910)npx+2,trim(tablename),ltf1
3910      format('plot for [i=3:',i2,'] $',a,' using 2:i',&
               ' with lines ls (i-2+',i2,') title columnheader(i)') 
       endif
    else
! plot line with points at every linewp-1 calculated point
       if(anpax.eq.1) then
          write(21,3600)npx+2,trim(tablename),ltf1,graphopt%linewp-1
3600      format('plot for [i=3:',i2,'] $',a,' using i:2',&
               ' with lp ls (i-2+',i2,') pi ',i3,' title columnheader(i)') 
       else
          write(21,3610)npx+2,trim(tablename),ltf1,graphopt%linewp-1
3610      format('plot for [i=3:',i2,'] $',a,' using 2:i',&
               ' with lp ls (i-2+',i2,') pi ',i3,' title columnheader(i)') 
       endif
    endif
!    write(*,*)'SMP2 linespoint increment 1:',graphopt%linewp-1
!=================================================================
! we come here after plotting a Scheil diagram above, try including appfiles
4000 continue
!=================================================================
! plot command from appfil
    if(appfil.gt.0) then
!       write(*,*)'ocplot2B appending a file at label 3912'
! try to avoid overlapping keys ...
! The "restore" for x/yrange means the scaling from the "plot for"
! will be used also for the appended data
       write(21,3912)trim(graphopt%font)
3912   format('set key bottom right font "',a,',12"'/&
            'set xrange restore'/'set yrange restore')
       if(appfiletyp.eq.2) then
! just one line with plot for ... 
! the data to append is already copied as a table
          write(21,'(a)')trim(applines(1))
! if applines>0 write those lines before "unset"
          if(multibuffline.gt.0) then
             do iz=1,multibuffline-1
                write(21,'(a)')trim(multibuffer(iz))
             enddo
          endif
          multibuffline=0
          write(21,'(a)')'unset multiplot #appfiletype 2'
!          write(21,'(a)')'unset multiplot'
          close(appfil)
          appfil=0
       elseif(multibuffline.gt.0) then
! these are "plot "-" ... lines, not connected to "set multiplot"
          write(21,'(a,2i7)')'# not appfiletype 2, multibufline, nofapl',&
               multibuffline,nofapl
          do iz=1,multibuffline-1
             write(21,'(a)')trim(multibuffer(iz))
          enddo
          multibuffline=0
          do jj=1,nofapl
             write(21,'(a)')trim(applines(jj))
          enddo
          write(21,'(a)')'unset multiplot #appfiletype 1'
          close(appfil)
          appfil=0
       endif
    endif
! if the plot command is "plot '-' ... then
! copy the data from the append file, it should be correctly formatted
! as I understand appfil muste always be 0 here ... no
    if(appfil.gt.0) then
       if(multibuffline.gt.0) then
          write(*,*)' *** ocplot2B appending multiplot',multibuffline
          do iz=1,multibuffline-1
             write(21,'(a)')trim(multibuffer(iz))
          enddo
          multibuffline=0
!          write(21,'(a)')'unset multiplot # closing appfil'
       endif
! appfile header lines
       write(*,*)'Appfile header lines: ',nofapl
       ic=0
       write(21,'(a)')'# Copying appfile data'
! these line contain the 'plot "-" ..."
       do jj=1,nofapl
          write(21,'(a)')trim(applines(jj))
       enddo
1900   continue
! this is copying the actual data to plot from the append file.
!       write(*,*)'Copying appfile data'
       read(appfil,884,end=1910)appline
884    format(a)
       ic=ic+1
! skip pause mouse
       if(appline(1:12).eq.'pause mouse ') goto 1900
       write(21,884)trim(appline)
       goto 1900
! end of copying appfile data
1910   continue
       write(*,*)'Appended ',ic,' data lines'
!       if(multibuffline.gt.0) then
! ocplot2B  add multiple plot "multplot" commands ....
!          write(*,*)'ocplot2B adding prevous multplot ...',multibuffline
!          do iz=1,multibuffline-1
!             write(21,884)trim(multibuffer(iz))
!          enddo
!       endif
!       write(21,'(a)')'unset multiplot'
       write(21,'(a)')'unset multiplot # closing appfil 3'
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
!    write(kou,*)'Gnuplot command file: ',pfc(1:kk+4)
    if(graphopt%gnutermsel.ne.1) then
       write(kou,*)'Graphics output file: ',pfh(1:kk+4)
    endif
! plotonwin is set by compiler option, if 1 we are running microspft windows
    if(lines_excluded.gt.0) write(kou,11)lines_excluded
11  format('SMP Some calculated lines are excluded from the plot',i5)
    if(plotonwin.eq.1) then
! call system without initial "gnuplot " keeps the window !!!
       if(btest(graphopt%status,GRKEEP)) then
!          write(*,*)'plot command: "',gnuplotline(9:k3),'"'
!          write(*,*)'Trying to spawn: ',trim(gnuplotline)
!          call system(gnuplotline(9:k3))
! spawn plot on Windows ?? NOT ISO-TERMAL DIAGRAM
!          write(*,*)'executing command: "start /B '//trim(gnuplotline)
!          call execute_command_line('start /B '//gnuplotline(9:k3))
    write(*,*)'ocplot2B executing command: "start /B '//trim(gnuplotline)//'"'
          call execute_command_line('start /B '//trim(gnuplotline))
! WORKS WITH OCPLOT3B
!          call execute_command_line('start /B '//trim(gnuplotline))
       else
!          write(*,*)'plot command: "',gnuplotline(1:k3),'"'
!          call system(gnuplotline)
          write(*,*)'ocplot2B executing command: "'//trim(gnuplotline)//'"'
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

!\addtotable subroutine ocplot3
!\begin{verbatim}
!  subroutine ocplot3(ndx,pltax,filename,maptop,axarr,graphopt,pform,&
  subroutine ocplot3(ndx,pltax,filename,mastertop,axarr,graphopt,&
       version,ceq)
! special to plot isothermal sections (two columns like x(*,cr) x(*,ni))
!   or other diagrams with two extensive variable on the axis
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
    integer ii,jj,point,plotp,lines,eqindex,lasteq,nooftups,lokcs,jph,same,kk
    integer, parameter :: maxval=4000,mazval=100
    double precision, allocatable :: xval(:,:),yval(:,:),zval(:,:),tieline(:,:)
    integer offset,nofeq,sumpp,last,nofinv,ntieline,mtieline,noftielineblocks
    double precision xxx,yyy
    character*64, dimension(:), allocatable :: phaseline
! TO BE REPLACED BY GNUTERMSEL: pform
!    character pform*8
    integer, allocatable :: plotkod(:),lineends(:)
    character xax1*8,xax2*24,yax1*8,yax2*24,axis1*32,axisx(2)*32,axisy(2)*32
    character phname*32,encoded*1024,axis*32
    character lid(2,maxsame)*24
    integer nooflineends,ephl,invnode
!
! do not change mastertop!
    maptop=>mastertop
! xval and yval and ccordinates to plot, 
! points on one line is       xval(1,jj),yval(1,jj)
! points on the other line is xval(2,jj),yval(2,jj)
! zval is an occational point zval(1,kk),zval(2,kk) at line ends (invariants)
! plotkod is a specific code for each point
!    normal point just the index, -1 means point should be suppressed
!    -100 or -101 is a tieline; -1000 an invariant
! lineends is the index in xval/yval for an end of line
    if(.not.associated(maptop)) then
       write(*,*)'No data to plot'
       goto 1000
    endif
    allocate(xval(2,maxval))
    allocate(yval(2,maxval))
    allocate(plotkod(maxval))
    allocate(zval(2,maxsame))
    allocate(lineends(maxsame))
! phaseline is the name of the phases stable along a line
    allocate(phaseline(maxsame))
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
          write(kou,*)' Warning, tie-lines may be wrong'
!          graphopt%tielines=0
       endif
    endif
! same is incremented for each line
    same=0
    sumpp=0
    plottop=>maptop
    curtop=>plottop
    nooftups=nooftup()
!    nooftup=noofphasetuples()
100 continue
    if(.not.allocated(curtop%linehead)) goto 500
    lines=size(curtop%linehead)
    results=>plottop%saveceq
    invnode=0
    if(btest(plottop%status,MAPINVARIANT)) then
       invnode=size(plottop%linehead)
       write(*,*)'ocplot3 invariant node 1',invnode
    endif
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
          ephl=1
          equil: do jph=1,nooftups
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
                phaseline(same)(ephl:)=phname
                ephl=min(len_trim(phaseline(same))+2,62)
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
          phaseline(same)='tie-line'
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
          ephl=1
          nodeequil: do jph=1,nooftups
             lokcs=phasetuple(jph)%lokvares
             if(curceq%phase_varres(lokcs)%phstate.ge.PHENTSTAB) then
! plotkod set to -1 to indicate monovariant (not invariant)
                plotkod(plotp)=-1
                call get_phasetup_name(jph,phname)
                phaseline(same)(ephl:)=phname
                ephl=min(len_trim(phaseline(same))+2,62)
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
! jump back to label 100
       goto 100
    endif
!    do jj=1,same
!       write(*,23)'phases: ',same,jj,trim(lid(1,jj)),trim(lid(2,jj))
!    enddo
!------------------------------------------------
! Jump here if there is a line with illegal lineid
772 continue
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
         lid,phaseline,filename,graphopt,version,encoded)
    deallocate(xval)
    deallocate(yval)
    deallocate(plotkod)
1000 continue
    return
  end subroutine ocplot3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
  subroutine ocplot3B(same,nofinv,lineends,nx1,xval,ny1,yval,nz1,zval,plotkod,&
       pltax,lid,phaseline,filename,graphopt,version,conditions)
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
! phaseline is phase names of phase stable along the line
! filename is intermediary file (maybe not needed)
! graphopt is graphics option record
! maptop is map_node record with all results
! REPLACED BY gnutermsel pform is type of output (screen/acrobat/postscript/gir)
! conditions is a text with conditions for the calculation
    implicit none
!    character pltax(*)*(*),filename*(*),pform*(*),lid(nx1,*)*(*),conditions*(*)
    character pltax(*)*(*),filename*(*),lid(nx1,*)*(*),conditions*(*)
    character version*(*)
    character phaseline(*)*(*)
    type(graphics_options) :: graphopt
    integer same,plotkod(*),nx1,ny1,nz1,nofinv
    integer lineends(*)
    double precision xval(nx1,*),yval(ny1,*),zval(nz1,*)
!\end{verbatim}
    integer, parameter :: maxcolor=200,mofapl=100
    integer ii,jj,kk,jph,offset,n1,nofapl,ltf2
    type(graphics_textlabel), pointer :: textlabel
    character gnuplotline*256,date*12,mdate*12,title*128
    character deftitle*64,backslash*2
! pointincremet emergy fix
    character piincrement*8
    character labelkey*64,applines(mofapl)*128,appline*128,pfc*128,pfh*128
    integer sumpp,np,appfil,ic,nnv,kkk,lcolor(maxcolor),iz,again
    integer done(maxcolor),foundinv,fcolor,k3
    character color(maxcolor)*24,rotate*16,labelfont*32,linespoints*12
    integer naptitle,apptitles(maxcolor)
! we plot monovariant twice, once with border once with filledcurves!!
!    integer noofmono,jjj,monovariant2(100)
    integer, parameter :: monovariantborder=11
    integer xtieline,xmonovariant,lz
! now a global variable
!    character monovariant*6
! Gibbs triangle variables
    logical plotgt,appgt
    double precision sqrt3,xxx,yyy,xmax,ltic,xf,yf
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
    ltf2=0
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
! terminal 1 is screen without any output file, add PDF as comment
       write(21,841)trim(graphopt%gnuterminal(graphopt%gnutermsel)),&
            trim(graphopt%gnuterminal(3))
841    format('set terminal ',a/'#set terminal ',a/'#set output "ocgnu.pdf"')
!841    format('set terminal ',a)
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
!    write(*,*)'Plot heading 3? ',btest(graphopt%status,GRNOTITLE)
    if(btest(graphopt%status,GRNOTITLE)) then
       write(21,128)trim(title),trim(conditions),trim(graphopt%font)
    else
       call replace_uwh(conditions)
       write(21,129)trim(title),trim(conditions),trim(graphopt%font)
    endif
    call replace_uwh(pltax(1))
! Plot grid?
    if(graphopt%setgrid.eq.1) write(21,777)
777 format('set grid')
    write(21,130)graphopt%xsize,graphopt%ysize,&
         trim(pltax(1)),trim(labelkey)
128 format('#set title "',a,' \n #',a,'" font "',a,',10"')
129 format('set title "',a,' \n ',a,'" font "',a,',10"')
130 format('set origin 0.0, 0.0 '/&
         'set size ',F8.4', ',F8.4/&
         'set xlabel "',a,'"'/&
         'set key ',a)
    call replace_uwh(pltax(2))
    if(plotgt) then
! OC logo added by Catalina Pineda
! when Gibbs triangle the ylabel and logo must be placed carefully
! THIS IS THE Y-AXIS WITH 60 degrees angle
       write(21,131)trim(pltax(2)), 0.15*xmax, 0.37*xmax
131    format('set label "',a,'" at ',F8.4,',',F8.4,' rotate by 60 '/&
! Help with stackoverflow to fix nice logo independent of plot size!
         'set label "~O{.0  C}" at screen 0.02, 0.03 font "Garamond Bold,20"')
!         'set label "~O{.0  C}" at graph -0.1, -0.1 font "Garamond Bold,20"')
! we should also enforce same length of X and Y axis !!!
    else
! SQUARE DIAGRAM
       write(21,132)trim(pltax(2))
132    format('set ylabel "',a,'"'/&
! Help with stackoverflow to fix nice logo independent of plot size!
         'set label "~O{.0  C}" at screen 0.02, 0.03 font "Garamond Bold,20"')
!        'set label "~O{.0  C}" at graph -0.1, -0.1 font "Garamond Bold,20"')
    endif
    lz=graphopt%linetype
    write(21,133)lz,lz,lz,lz,lz,lz,lz,lz,lz,lz
133 format('# if the value after solid is 0 the monovariants are transparent'/&
         'set style fill transparent solid 1'/&
         'set style line 1 lt ',i2,' lc rgb "#000000" lw 2 pt 10'/&
         'set style line 2 lt ',i2,' lc rgb "#00C000" lw 2 pt 2'/&
         'set style line 3 lt ',i2,' lc rgb "#4169E1" lw 2 pt 7'/&
         'set style line 4 lt ',i2,' lc rgb "#FF0000" lw 2 pt 3'/&
         'set style line 5 lt ',i2,' lc rgb "#00FFFF" lw 2 pt 10'/&
         'set style line 6 lt ',i2,' lc rgb "#FF00FF" lw 2 pt 5'/&
         'set style line 7 lt ',i2,' lc rgb "#804080" lw 2 pt 6'/&
         'set style line 8 lt ',i2,' lc rgb "#00C000" lw 2 pt 8'/&
         'set style line 9 lt ',i2,' lc rgb "#C0C0C0" lw 2 pt 1'/&
        'set style line 10 lt ',i2,' lc rgb "#DAA520" lw 2 pt 4')
! add some useful things for manual maniplulation of graph
    write(21,8000)
8000 format(/'# Some useful GNUPLOT commands for editing the figure'/&
          '# This is a dashed line (on pdf/wxt):'/&
          '# set style line 15 lt 0 lc rgb "#C8C800" lw 2 pt 2'//&
          '# set pointsize 0.6'/&
          '# set label "text" at 0.5, 0.5 rotate by 60 font "arial,12"'/&
          '# set xrange [0.5 : 0.7] '/&
          '# Adding manually a line and keep scaling:'/&
          '# set arrow x0, y0 to x1,y1 nohead linestyle 1'/&
          '# set multiplot'/&
          '# set xrange [] writeback'/&
          '#  ... plot someting'/&
          '# set xrange restore'/&
          '#  ... plot more using same axis scaling '/&
          '# set nomultiplot'/)
!
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
!          write(labelfont,178)int(10*textlabel%textfontscale)
!178       format(' font "Sans,',i2,'" ')
          write(labelfont,178)trim(graphopt%font),&
               int(10*textlabel%textfontscale)
178       format(' font "',a,',',i2,'" ')
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
!          write(*,*)'ignoring append file line ',trim(appline)
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
          write(*,*)'missing color in ',ii,' out of ',2*same,' set to 9'
          lcolor(ii)=9
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
!          write(21,208)trim(graphopt%lowerleftcorner),-0.14
!208       format('set label "',a,'" at graph ',F10.4,', -0.05 ')
          write(21,208)trim(graphopt%lowerleftcorner),-0.05
208       format('set label "',a,'" at graph ',F10.4,', -0.03 ')
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
    write(21,*)
! here we should start from the value in graphopt%linett
    ii=0
    kk=graphopt%linett-1
    if(kk.ne.0) then
       write(*,*)'Ignoring manipulation of line colors'
    endif
! if graphopt%linestyle=0 use lines, otherwise linespoints
! this never worked ...
!    if(graphopt%linestyle.eq.0) then
    linespoints='lines ls'
    piincrement=' '
    if(graphopt%linewp.gt.1) then
! this should add a symbol at each calculated line but it does not work (yet)
!       linespoints='lp lt '
       linespoints='lp ls '
       write(piincrement,'(" pi ",i3,1x)')graphopt%linewp-1
    endif
!    write(*,*)'smp2B linesplot increment 2: "',piincrement,'"',graphopt%linewp
!    write(*,*)'SMP2B plotting lines: ',trim(linespoints),graphopt%linewp
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
! SUDDENLY lcolor(ii) is not set null !!
!                write(*,*)'SMP2B ocplot3B 1: ',ii,'"',linespoints,'"'
!                write(*,*)'SMP2B ocplot3B 2: ',lcolor(ii),color(lcolor(ii))
                write(21,309)trim(linespoints),lcolor(ii),&
                     trim(piincrement),trim(color(lcolor(ii))),backslash
             endif
             naptitle=naptitle+1
             apptitles(naptitle)=lcolor(ii)
!309          format('plot "-" using 1:2 with ',a,' ls ',i2,' title "',a,'"',a)
309          format('plot "-" using 1:2 with ',a,1x,i2,1x,a,' title "',a,'"',a)
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
                      write(21,320)'lines ls',monovariantborder,' ',backslash
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
                      write(21,320)'lines ls',fcolor,' ',backslash
                      xtieline=jj
!                      write(*,*)'SMP xtieline 2: ',xtieline
                   elseif(xtieline.ne.jj) then
                      write(21,320)'lines ls',fcolor,' ',backslash
                   else
                      write(21,299)'lines',fcolor,trim(color(12)),backslash
299                   format('"" using 1:2 with ',a,' ls ',i2,&
                           ' title "',a,'" ',a)
!                      write(*,*)'SMP xtieline 5:',jj,xtieline
                   endif
                else
! normal line with no title
!                   write(*,320)trim(linespoints),fcolor,backslash
                   write(21,320)trim(linespoints),fcolor,&
                        trim(piincrement),backslash
                endif
!320             format('"" using 1:2 with ',a,' ls ',i2,' notitle ',a)
320             format('"" using 1:2 with ',a,1x,i2,1x,a,' notitle ',a)
             else 
! we have a new line withou title
                if(fcolor.eq.11) then
                   if(kk.eq.1) then
! first time plotting a monovariant use thick lines
                      write(21,320)'lines ls',monovariantborder,' ',backslash
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
                      write(21,320)'lines ls',fcolor,' ',backslash
                      xtieline=jj
!                      write(*,*)'SMP xtieline 3: ',xtieline
                   else
! if kk=2 add title if xtieline=jj
                      if(jj.eq.xtieline) then
                         write(21,331)'lines ls',fcolor,&
                              trim(color(lcolor(ii))),backslash
!                         write(*,*)'SMP xtieline 4:',jj,xtieline
                      else
                         write(21,320)'lines ls',fcolor,' ',backslash
                      endif
                   endif
                else
! any normal line add title
                   write(21,331)trim(linespoints),fcolor,&
                        trim(piincrement),trim(color(lcolor(ii))),backslash
                endif
                naptitle=naptitle+1
                apptitles(naptitle)=lcolor(ii)
!331             format('"" using 1:2 with ',a,' ls ',i2,' title "',a,'"',a)
331             format('"" using 1:2 with ',a,1x,i2,1x,a,' title "',a,'"',a)
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
! we have to include the scaling factors graphopt%scalefact(1,2)
    xf=one; yf=one
    if(graphopt%scalefact(1).ne.one) xf=graphopt%scalefact(1)
    if(graphopt%scalefact(2).ne.one) yf=graphopt%scalefact(2)
!    write(*,*)'SMP2B x/y factors: ',xf,yf
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
             write(21,549)xf*xval(1,sumpp),yf*yval(1,sumpp)
             write(21,549)xf*xval(2,sumpp),yf*yval(2,sumpp)
             write(21,549)xf*zval(1,foundinv),yf*zval(2,foundinv)
             write(21,549)xf*xval(1,sumpp),yf*yval(1,sumpp)
! we are at the end of a line, write a blank line
             write(21,548)jj,trim(phaseline(jj))
!             write(21,548)jj
548          format('e '//'# end of monovariant',i5,2x,a/)
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
                write(21,549)xf*xval(kk,sumpp),yf*yval(kk,sumpp)
!549             format(2e15.6,4i7)
549             format(2(1pe16.6),4i7)
! plotkod -101 means tieline
! UNFINISHED: VALGRIND indicates plotkod(sumpp) is uninitiallized??
                if(plotkod(sumpp).eq.-101) write(21,552)
552             format(1x)
             enddo
! we are at the end of a line, write a blank line
             write(21,551)jj,trim(phaseline(jj))
!             write(21,551)jj
551          format('e '//'# end of line',i5,2x,a/)
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
       if(appline(1:12).eq.'pause mouse ') then
          write(*,*)'reading appendfile ends at "puase mouse"'
          goto 1910
       else
          write(21,884)trim(appline)
          goto 1900
       endif
1910   continue
       write(*,*)'Appended ',ic,' data lines'
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
    if(lines_excluded.gt.0) write(kou,11)lines_excluded
11  format('SMP Some calculated lines excluded from the plot',i5)
! plotonwin set by compiler option, 1 means windows 
    if(plotonwin.eq.1) then
       if(btest(graphopt%status,GRKEEP)) then
! this is a TERNARY PLOT with 2 extenive axis
!          write(*,*)'executing command '//trim(gnuplotline(9:))
!          call system(gnuplotline(9:))
   write(*,*)'ocplot3B executing Command: "start /B '//trim(gnuplotline)//'"'
! WORKS WITH OCPLOT3B
          call execute_command_line('start /B '//trim(gnuplotline))
       else
          write(*,*)'ocplot3B executing command '//trim(gnuplotline)
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

!\addtotable subroutine calc_diagram_point
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

!\addtotable logical function abbr_phname_same
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

!\addtotable subroutine get_plot_conditions
!\begin{verbatim}
  subroutine get_plot_conditions(text,ndx,axarr,ceq)
! extracts the conditions from ceq and replaces those that are axis variables
    implicit none
    character text*(*)
    integer ndx
    type(map_axis), dimension(*) :: axarr
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
! if line too long (>200) divide in middle
    jj=len_trim(text)
    if(jj.gt.250) then
! text 1:ip
       ip=jj/3
       call find_space_in_text(text,ip,10)
       text(ip+4:)=text(ip:)
       text(ip:ip+3)=' \n '
! maybe some character get lost ... no one will notice
! text ip+4:2*ip+8
       jp=2*ip+4
       call find_space_in_text(text,jp,10)
       text(jp+4:)=text(jp:)
       text(jp:jp+3)=' \n '
!       write(*,*)'Dividing condition text 3 parts',ip,jp,jj
    elseif(jj.gt.100) then
!       write(*,*)'Dividing condition text in the middle',jj
       jj=jj/2
       call find_space_in_text(text,jj,10)
       text(jj+4:)=text(jj:)
       text(jj:jj+3)=' \n '
    endif
1000 continue
    return
  end subroutine get_plot_conditions

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_space_in_text
!\begin{verbatim}
  subroutine find_space_in_text(text,jp,maxc)
! moves jp max +/-maxc charactres to find a space (? or , or : or )
! If none found do not change jp
    implicit none
    character text*(*)
    integer jp,maxc
!\end{verbatim}
    integer ap
    ap=jp
    add: do ap=jp,jp+maxc
       if(text(ap:ap).eq.' ') goto 900
    enddo add
    sub: do ap=jp-1,jp-maxc,-1
       if(text(ap:ap).eq.' ') goto 900
    enddo sub
! no space found
    ap=jp
900 continue
!    write(*,*)'SMP2B: ',text(jp-maxc:jp+maxc),jp,ap
    jp=ap
    return
  end subroutine find_space_in_text

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable double precision function get_axis_phase_value
!\begin{verbatim}
  double precision function get_axis_phase_value(phase,axis,axarr,ceq)
! extacts the condition for one axis and if is something like X(A)
! it changes to x(phase,A) and extracts and returns that value !!
! if it is a potential like T it just returns its value
    implicit none
    character phase*(*)
    integer axis
    type(map_axis), dimension(*) :: axarr
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer seqz,ip
    character symbol*64,dummy*64
    double precision value
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: svrrec,svr2
!
    value=zero
! find the condition for axis "axis"
    seqz=axarr(axis)%seqz
    call locate_condition(seqz,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
    svrrec=>pcond%statvar(1)
!    write(*,*)'Value of axis',axis,' for phase ',trim(phase),svrrec%statevarid
    if(svrrec%statevarid.le.9) then
! it is a potential, extract its current value
       symbol=' '
       ip=1
       call encode_state_variable(symbol,ip,svrrec,ceq)
       if(gx%bmperr.ne.0) goto 1000
       call get_state_var_value(symbol,value,dummy,ceq)
    elseif(svrrec%argtyp.eq.1) then
! it is an extensive variable (%statevarid>=10) for a component such as X(A)
! A smarter way is to modify svrrec to insert phase/compset index ...
! edit the phase into the symbol
       symbol=' '
       ip=1
       call encode_state_variable(symbol,ip,svrrec,ceq)
       if(gx%bmperr.ne.0) goto 1000
       ip=index(symbol,'(')
       dummy=symbol(ip+1:)
       symbol(ip+1:)=trim(phase)//','
       ip=len_trim(symbol)
       symbol=symbol(1:ip)//dummy
!       write(*,*)'SMP2B value of: ',trim(symbol)
!       call get_stable_state_var_value(symbol,value,dummy,ceq)
       call get_state_var_value(symbol,value,dummy,ceq)
       if(gx%bmperr.ne.0) goto 1000
    else
       write(*,*)'Illegal axis type: ',svrrec%statevarid,svrrec%argtyp
       gx%bmperr=4399
    endif
1000 continue
!    write(*,*)'SMP2B value of: ',trim(symbol),value
    get_axis_phase_value=value
    return
  end function get_axis_phase_value
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_stored_equilibria
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
    type(meq_setup), pointer :: lineeq
    integer kl,ll,jax,nax,jj,jk,phtupix
    double precision, dimension(:), allocatable :: axxx
    type(gtp_state_variable), pointer :: svrrec
    logical once
    character*100 phases
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
    write(kou,101)mapnode%seqx,mapnode%nodeceq%tpval(1),mapnode%noofstph,&
         mapnode%savednodeceq,mapnode%lines
!         mapnode%savednodeceq,mapnode%status,mapnode%lines
101 format(' Mapnode: ',i3,' at T=',F10.2,', ',i2,' phases, ceq saved ',&
         i5,', exting lines: ',i2)
    do kl=1,mapnode%lines
       if(.not.associated(mapnode%linehead(kl)%end)) then
          if(mapnode%linehead(kl)%termerr.gt.0) then
             write(kou,105)kl,mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria,&
                  mapnode%linehead(kl)%termerr
105          format('  Line ',i3,', id: ',i3,' with ',i5,&
                  ' equilibria ended with error: ',i6)
          else
             write(kou,110)kl,mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria
110          format('  Line ',i3,', id: ',i3,' with ',i5,&
                  ' equilibria ending at axis limit')
          endif
       else
          ll=mapnode%linehead(kl)%end%seqx
          write(kou,120)kl,mapnode%linehead(kl)%lineid,&
               mapnode%linehead(kl)%number_of_equilibria,ll
120       format('  Line ',i3,', id: ',i3,' with ',i5,&
               ' equilibria ending at node ',i3)
       endif
       if(btest(mapnode%linehead(kl)%status,EXCLUDEDLINE)) then
          write(*,*)'Line excluded'
          cycle
       endif
       ll=mapnode%linehead(kl)%first
!       write(*,*)'SMP2B allcrach 1: ',ll
! BOS 191224 add phase names
! tzero lines has no meqrec%phr allocated NOTE kl is K-EL not K-ETT
! for tzero lines it is listed line 3 and 4 although there are only 2 ????
       if(ll.gt.0 .and. allocated(mapnode%linehead(kl)%meqrec%phr)) then
! only if there is an link to a linehead
          lineeq=>mapnode%linehead(kl)%meqrec
          phases=' '
          jk=1
          do jj=1,lineeq%nstph
!             write(*,*)'SMP2B allcrach 2: ',jj,lineeq%stphl(jj),&
!                  lineeq%phr(lineeq%stphl(jj))%phtupix
             phtupix=lineeq%phr(lineeq%stphl(jj))%phtupix
             if(jk.lt.72) then
                call get_phasetup_name(phtupix,phases(jk:))
                jk=len_trim(phases)+2
             else
                phases(jk:)=' ... more'
             endif
          enddo
          write(*,*)'Phases: ',trim(phases)
! BOS 191224 end add phase names
!       write(*,*)'list first equilibrium ',ll
!       write(*,*)'axis: ',mapnode%number_ofaxis
          if(.not.allocated(results%savedceq)) then
             write(*,*)'Cannot find link to saved equilibria! '
          else
             once=.true.
             write(kou,140)
140          format('Saved ceq     link       T            X')
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
150             format(2i9,f9.2,5(1pe13.5))
                ll=thisceq%nexteq
             enddo ceqloop
          endif
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

!\addtotable subroutine amend_stored_equilibria
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
    logical all
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
    last=len(cline)
    call gparcdx('Only excluded? ',cline,last,1,ch1,'Y','?PLOT options')
    if(ch1.eq.'Y' .or. ch1.eq.'y') then
       all=.FALSE.
    else
       all=.TRUE.
    endif
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
    write(kou,102)mapnode%seqx,mapnode%nodeceq%tpval(1),&
         mapnode%noofstph,mapnode%savednodeceq
102 format(' Mapnode: ',i5,' at T=',F10.2,' with ',i2,&
         ' stable phases, ceq saved in ',i5)
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
             write(kou,105)kl,mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria,&
                  mapnode%linehead(kl)%termerr,status,trim(phline)
105          format('  Line ',i3,', id: ',i3,' with ',i5,&
                  ' equilibria ended with error: ',i6,2x,a/'  with phases: ',a)
          else
             write(kou,110)kl,mapnode%linehead(kl)%lineid,&
                  mapnode%linehead(kl)%number_of_equilibria,status,trim(phline)
110          format('  Line ',i3,', id: ',i3,' with ',i5,&
                  ' equilibria ending at axis limit.',2x,a/'  with phases: ',a)
          endif
       else
          ll=mapnode%linehead(kl)%end%seqx
          write(kou,120)kl,mapnode%linehead(kl)%lineid,&
               mapnode%linehead(kl)%number_of_equilibria,ll,status,trim(phline)
120       format('  Line ',i3,', id: ',i3,' with ',&
               i5,' equilibria ending at node ',i3,2x,a/'  with phases: ',a)
       endif
       ll=mapnode%linehead(kl)%first
! if deleted ask for Restore, else ask for Keep or Delete
       last=len(cline)
       if(btest(mapnode%linehead(kl)%status,EXCLUDEDLINE)) then
          call gparcdx(' *** Include this line? ',cline,last,1,ch1,'N',&
               '?PLOT options')
       elseif(all) then
          call gparcdx('Exclude this line? ',cline,last,1,ch1,'N',&
               '?PLOT options')
       else
          cycle lineloop
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
    mapnode=>mapnode%previous
    if(.not.associated(mapnode,localtop)) then
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

!\addtotable subroutine list_csv_table
!\begin{verbatim} %
  subroutine list_csv(np,nrv,nlinesep,linesep,pltax,xax,anpax,anpdim,anp,lid,&
       phaseline,title,filename,version,encoded1)
! list results from a STEP as an CSV (Comma Separated Values) table
! called from ocplot2 when all values extraced
! np number of regions
! nrv total number of lines with calculated values
! nlinesep number of lines for each region
! linesep(j..nlinesep) is index of last line for region j
! pltax are heading of table columns, each region separate headings
! xax values on single value axis
! anpax not used
! anpdim first dimension of anp
! anp values of axis with (possibly) multiple values
! lid column headings for each region ?
! phaseline ??
! title of plot
! filename output file
! graphopt record
! version of OC
! encoded1 all conditions
!     
! use ov
    integer np,nrv,nlinesep,anpax,anpdim
    integer linesep(*)
    character filename*(*),title*(*),version*(*),encoded1*(*)
    character pltax(*)*(*),phaseline(*)*(*),lid(*)*(*)
    double precision xax(*),anp(anpdim,*)
! not needed?
!    character lid(*)*(*)
!\end{verbatim}
    integer jj,kk
!
!    write(*,'(a,10i5)')'SMP2B: in list_csv',np,nrv,nlinesep,anpax,anpdim
!    write(*,'(a,a,a,i3,1pe12.4)')'SMP2B: plotfile: "',trim(filename),'"',&
!         len_trim(filename),rnone
    if(filename(1:1).ne.' ') then
       open(22,file=filename,access='sequential',status='unknown',err=1100)
       lut=22
    else
       lut=kou
    endif
! header
    if(np.eq.1) then
       write(lut,100)trim(pltax(1)),trim(lid(np))
    elseif(np.le.100) then
       write(lut,101,advance='no')trim(pltax(1)),(trim(lid(jj)),jj=1,np-1)
       write(lut,102)trim(lid(np))
    else
       write(*,*)'Cannot tablulate when more than 100 variables ...'
    endif
100 format('"',a,'", "',a,'"')
101 format('"',a,'"',100(',"',a,'"'))
102 format(a,'"')
! loop for lines
    do jj=1,nrv
       write(lut,200,advance='no')xax(jj)
200    format(1PE13.5)
201    format(',',1PE13.5)
202    format(',')
! loop for columns except first and last value       
       do kk=1,np-1
          if(anp(kk,jj).ne.rnone) then
             write(lut,201,advance='no')anp(kk,jj)
          else
             write(lut,202,advance='no')
          endif
       enddo
       if(anp(np,jj).ne.rnone) then
          write(lut,201)anp(np,jj)
       else
          write(lut,202)
       endif
    enddo
    if(lut.ne.kou) close(lut)
1000 continue
    return
! failed open output file
1100 continue
    write(*,*)'Cannot open file: ',trim(filename)
    goto 1000
  end subroutine list_csv

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine ocappfixlabels
!\begin{verbatim}
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
!\end{verbatim}
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
          ip=index(applines(j1),' lines ls ')+10
          if(ip.le.0) then
             write(*,*)'Old color in appfile: ',trim(applines(j1))
             stop
          else
! reuse the same color for something else
! changing like this created problems below ... try reuse old color
!             write(applines(j1)(ip:ip+1),'(i2)')same+1
             applines(j1)(ip:ip+1)='01'
             write(*,*)'New color in appfile: ',trim(applines(j1))
          endif
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
                write(*,79)'Cannot find old ls: ',oldls,k1,found,j1,&
                     trim(applines(j1))
79              format(a,4i3,' in ',a)
                write(*,'(10(i2,i3))')(changels(1,k1),changels(2,k1),k1=1,found)
! replace colow with 01
                ip=index(applines(j1),' lines ls ')
                applines(j1)(ip+10:ip+11)='01'
!                stop
             endif
! write the new ls number in applines(j1)
100          continue
!             write(applines(j1)(ip:ip+1),'(i2)')changels(2,k1)
! line above must be wrong, changed to that below 2021.03.08/BoS, then removed 
!             write(applines(j1)(ip+10:ip+11),'(i2)')changels(2,k1)
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

!\addtotable subroutine stvarix 
!\begin{verbatim}
  subroutine stvarix(statevar,phaseline,encoded,nix,ixpos)
! extract the indices of corresponding wildcard state variables
! For example X(*,O) replace "*" with phase names in phline
! and search which index correspond to X(C1_MO2#2,O) in encoded ... suck
! statevar: with ONE wildcard (for a phase)
! phaseline: phase names separated by a space
! encoded: state variables returned by get_many(...) separated by a space
! nix: number of state variables in encoded
! ixpos: integer array with corresponding index for values
! NOTE is there is a # in the statevar then ALL values should be included
    implicit none
    integer nix,ixpos(*)
    character*(*) statevar, phaseline, encoded
!\end{verbatim}
    character sstring*48,phase*32,prefix*24,suffix*24,cha*1
    integer ip,jp,kp,lp,pix,vix
!    write(*,*)'stvarix: ',trim(statevar)
!    write(*,*)'phases:  "',trim(phaseline),'"'
!    write(*,*)'encoded: "',trim(encoded),'"'
! initiate ipos to zero
!    write(*,8)1,trim(statevar),len_trim(encoded),trim(encoded)
!8   format('smp2B: ',i1,' searching for "',a,'" in encoded with length ',i5/a/)
!
    ixpos(1:nix)=0
    ip=index(statevar,'*')
! if no wildcard skip
    if(ip.le.0) goto 1000
! debug check search string:
!    write(*,8)2,trim(statevar),len_trim(encoded),trim(encoded)
    prefix=statevar(1:ip-1)
    suffix=statevar(ip+1:)
!
    ip=1
    jp=index(phaseline,' ')
    outside: do while(jp.gt.ip)
! create search string with phase name replacing wildcard
       sstring=trim(prefix)//phaseline(ip:jp-1)//trim(suffix)
!       write(*,*)'Seach for: "',trim(sstring),'"',jp
       pix=0
       kp=1
       lp=index(encoded,' ')
!       write(*,*)'encoded item: ',trim(encoded(kp:lp-1)),kp,lp-1
       inside: do while(lp.gt.kp) 
          pix=pix+1
!          write(*,*)'Same? ',trim(sstring),' and ',trim(encoded(kp:lp-1)),pix
!          read(*,10)cha
10        format(a)
          if(trim(sstring).eq.trim(encoded(kp:lp-1))) then
!             vix=vix+1
! it seems simpler to indicate for each possible value that it is relevent
             ixpos(pix)=1
! than to return an array with the relevant values ...
!             ixpos(vix)=pix
!             write(*,*)'Found: ',trim(sstring),vix,ixpos(vix)
! select next phase name
             ip=jp+1
             jp=jp+index(phaseline(jp+1:),' ')
             cycle outside
          endif
! compare with next item in encoded
          kp=lp+1
          lp=lp+index(encoded(lp+1:),' ')
!          write(*,*)'Next encoded item: ',trim(encoded(kp:lp-1)),kp,lp-1
!          read(*,10)cha
       enddo inside
       write(*,*)'SMP2B: Cannot find: ',trim(sstring)
       gx%bmperr=4399; goto 1000
    enddo outside
!    write(*,70)nix,(ixpos(vix),vix=1,nix)
70  format('SMP nix: ',i3,30i2)
!    write(*,*)'vix mm: ',vix,(ixpos(vix),vix=1,nix)
1000 continue
    return
  end subroutine stvarix

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine hashtag_susphix 
!\begin{verbatim}
  subroutine hashtag_susphix(statevar,phaseline,encoded,nix,ixpos,ceq)
! replace # with non-suspended phase names in DGM(#)
! statevar: DGM(#) wildcard (for a phase)
! phaseline: phase names separated by a space
! encoded: state variables returned by get_many(...) separated by a space
! nix: number of state variables in encoded
! ixpos: integer array with corresponding index for values
    implicit none
    integer nix,ixpos(*)
    character*(*) statevar, phaseline, encoded
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character sstring*48,phase*32,prefix*24,suffix*24,cha*1
    integer ip,jp,kp,lp,pix,lpre,lenc,iph,ics
    double precision amfu
!    write(*,*)'SMP2B hastag_suspix: ',trim(statevar),nix
!    write(*,*)'phases:  "',trim(phaseline),'"'
!    write(*,*)'encoded: "',trim(encoded),'"'
! initiate ipos to zero
!    write(*,8)1,trim(statevar),len_trim(encoded),trim(encoded)
!8   format('smp2B: ',i1,' searching for "',a,'" in encoded with length ',i5/a/)
!
    lenc=len_trim(encoded)
!    write(*,*)'SMP2B in hashtag_susphix: ',trim(statevar)
!    write(*,'(a,a)')'SMP2B phaseline: ',trim(phaseline)
!    write(*,'(a,a)')'SMP2B encoded: ',lenc
!
    ixpos(1:nix)=0
    ip=index(statevar,'#')
! if no hashtag skip
    if(ip.le.0) goto 1000
! if hashtag not followed ) or , it indicate a composition set, skip
    if(.not.(statevar(ip+1:ip+1).eq.')' .or. statevar(ip+1:ip+1).eq.',')) &
         goto 1000
! debug check search string:
!    write(*,8)2,trim(statevar),len_trim(encoded),trim(encoded)
    prefix=statevar(1:ip-1)
    lpre=len_trim(prefix)
!    write(*,*)'SMP2B hastag1: "',prefix(1:lpre),'"',len(encoded)
!
    ip=1
    pix=0
    ixpos(1:nix)=0
!  we can ignore phaseline and we skip all phases in encoded that are suspended
    do while(ip.lt.lenc)
       jp=index(encoded(ip:),prefix(1:lpre))
       kp=index(encoded(ip:),' ')
       phase=encoded(ip+jp+lpre-1:ip+kp-3)
!       write(*,'(a,a,a,3i5)')'SMP2B hashtag 2: "',trim(phase),'"',ip,jp,kp
! update ip for next phase
       ip=ip+kp
!       write(*,*)'SMP2B rest: ',trim(encoded(ip:))
! find phase number/set       
!       call find_phase_by_name_exact(phase,iph,ics)
       call find_phase_by_name(phase,iph,ics)
       if(gx%bmperr.ne.0) then
          write(*,*)'SMP2B hashtag found nonexisting phase name: ',phase
       endif
       pix=pix+1
       if(test_phase_status(iph,ics,amfu,ceq).ge.PHDORM) then
! this phase is not suspended, it should be included
! it seems simpler to indicate for each possible value that it is relevent
          ixpos(pix)=1
! than to return an array with the relevant values ...
!          write(*,*)'not suspended: ',trim(phase),pix,ixpos(pix)
!       else
!          write(*,*)'suspended: ',trim(phase)
       endif
    enddo
1000 continue
    return
  end subroutine hashtag_susphix

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine replace_UWH
!\begin{verbatim}
  subroutine replace_uwh(text)
! replaces underscore by a hyphen for texts used in GNUPLOT
! replaces ampersand, &, by @
    implicit none
    character*(*) text
!\end{verbatim}
    integer jj
! replace _ by - in lid
    jj=-1
    do while(jj.ne.0)
! replace _ by - in lid because _ is interpreted as subscript (as LaTeX)
       if(jj.gt.0) text(jj:jj)='-'
       jj=index(text,'_')
    enddo
    jj=-1
    do while(jj.ne.0)
! replace & by z in lid because & is treated strangely by GNUPLOT
       if(jj.gt.0) text(jj:jj)='%'
       jj=index(text,'&')
    enddo
!    write(*,*)'SMP2B text without "_": ',trim(text)
    return
  end subroutine replace_uwh
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

