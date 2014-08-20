! Sundman's STEP, MAP and PLOT package (SMP)
!
MODULE smp
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
! Done: a first simple step and plot routine implemented
!
!-----------------------------
!
! To be done:
! - design of a storage structure of calculated equilibria
! - step (one variable axis)
! - map (two or more variable axis, follow lines with zero phasefraction)
! - graphical postprocessing using the stored equilibria
! 
!------------------------------
!
  use matsmin
!
  implicit double precision (a-h,o-z)
!
  TYPE ssm_node
! these records form a circular list with next pointing the the next record
! created.  The root record is the last created and root%next is the first
! ceqlink points to a record with all data for the equilibrium
! seqno and tval are redundant
! status can be used to suspend an equilinrium in plot or listings
     integer seqno,status
     type(ssm_node), pointer :: next
     type(gtp_equilibrium_data), pointer :: ceqlink
!     type(gtp_equilibrium_data), allocatable :: ceqlink
     double precision tval
  end TYPE ssm_node
!  
CONTAINS
  
  subroutine sstep1(stepopt,axvar,axval,resultlist,ceq)
! steps in variable axvar from axval(1) to axval(2) with increment axval(1)
! ceq is a datastructure with all relevant thermodynamic and related data
    implicit double precision (a-h,o-z)
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(ssm_node), pointer :: resultlist
    TYPE(ssm_node), pointer :: lnod
    dimension axval(3,*)
    character axvar(*)*(*),stepopt*(*),ch1*1,condtext*60
    logical form
!
    write(kou,7)'Entering step version 1: ',stepopt(1:len_trim(stepopt)),&
         axvar(1)(1:len_trim(axvar(1))),(axval(i,1),i=1,3)
7   format(/50('-')/a,1x,a,1x,a,3(1pe12.4))
! do the very simplest, assume T as axis, start from axval(1,1) to axval(2,1)
! with increments axval(3,1) and store calculated equilibria (complete ceq
! datastructures)
! check one can calculate an equilibrium at lower T
    ceq%tpval(1)=axval(1,1)
    write(kou,77)ceq%tpval(1)
77  format('Starting at T=',F10.2)
!    call calceq2(1,ceq)
!    if(gx%bmperr.ne.0) then
!       write(*,*)'Error calculating first equilibrium for step',gx%bmperr
!       goto 1000
!    endif
!    write(*,*)'Calculated with ceq no: ',ceq%eqno
! start the curcular link
!    nullify(resultlist)
    nstep=0
    ceq%tpval(1)=axval(1,1)
! turn off post gridtest after each calculation
!    if(btest(ceq%status,EQGRIDTEST)) then
!       write(*,*)'gridtest turned on by default',ceq%status
!    else
!       write(*,*)'gridtest turned off by default',ceq%status
!    endif
    ceq%status=ibclr(ceq%status,EQGRIDTEST)
!    if(btest(ceq%status,EQGRIDTEST)) then
!       write(*,*)'gridtest turned on',ceq%status
!    else
!       write(*,*)'gridtest turned off',ceq%status
!    endif
!-----------------------------------------------------
! return here for each new calculation (except last)
100 continue
    write(kou,*)' ----------- Calculate at T=',ceq%tpval(1)
! always reset unstable phases to their default constitution
!    call set_default_constitution(-1,0,0,ceq)
!    if(gx%bmperr.ne.0) goto 1000
! call without gridminimizer
    call calceq2(0,ceq)
    if(gx%bmperr.ne.0) then
! if this fails then try ahain without the grid minmizer.  
       gx%bmperr=0
       call calceq2(0,ceq)
       if(gx%bmperr.ne.0) then
! if this fails then try with the grid minmizer.  
! NOTE: composition sets may change
          write(kou,105)gx%bmperr,ceq%tpval(1),bmperrmess(gx%bmperr)
105       format(' *** Error : ',i5,' at T=',f10.2/a/'Trying harder')
          gx%bmperr=0
          call calceq2(1,ceq)
          if(gx%bmperr.ne.0) then
             gx%bmperr=0
! Here reset all phases to their default constitution
! Do not change their amounts
!    call set_default_constitution(-1,0,-1,ceq)
!    if(gx%bmperr.ne.0) goto 1000
             call calceq2(0,ceq)
             if(gx%bmperr.ne.0) then
                goto 1000
             endif
          endif
       endif
    endif
    call add_to_steplista(resultlist,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'ssm lista: ',resultlist%seqno,resultlist%next%seqno
! take a step
    nstep=nstep+1
    ceq%tpval(1)=ceq%tpval(1)+axval(3,1)
!    write(*,*)'New T=',ceq%tpval(1),axval(3,1),nstep
!    read(*,807)ch1
807 format(a)
! >>>>>
! insert the axis value as condition
    condtext='T='
    ip=3
    call wrinum(condtext,ip,15,0,ceq%tpval(1))
    ip=0
    call set_condition(condtext,ip,ceq)
    if(gx%bmperr.ne.0) goto 1000
! <<<<<
    if(ceq%tpval(1).lt.axval(2,1)) goto 100
! last point
    ceq%tpval(1)=axval(2,1)
    write(kou,*)' ----------- Last calculation at T=',ceq%tpval(1)
!    call calceq2(1,ceq)
! with first argument 0 there is no global minimization
    call calceq2(0,ceq)
    if(gx%bmperr.ne.0) then
! if this fails then try with the grid minmizer
!       write(*,*)'Trying harder'
       gx%bmperr=0
       call calceq2(1,ceq)
       if(gx%bmperr.ne.0) then
          gx%bmperr=0
          call calceq2(0,ceq)
          if(gx%bmperr.ne.0) then
             goto 1000
          endif
       endif
    endif
    call add_to_steplista(resultlist,ceq)
    nstep=nstep+1
! verify list
!    write(*,*)'Number of records: ',resultlist%seqno
!    lnod=>resultlist%next
!800 continue
!    write(*,*)'lnod loop: ',lnod%seqno,lnod%tval,lnod%ceqlink%tpval(1)
!    if(lnod%seqno.ne.resultlist%seqno) then
!       lnod=>lnod%next
!       goto 800
!    endif
!
! save on file
    write(*,*)'Saving not implemented yet'
    goto 1000
!    
    write(*,*)'Saving unformatted on ocstep.gtp'
    open(21,file='ocstep.gtp ',access='sequential',status='unknown',&
         form='unformatted')
    form=.FALSE.
    lnod=>resultlist%next
810 continue
! saving on file not implememted
!    write(*,*)'saving: ',lnod%seqno,lnod%tval,lnod%ceqlink%tpval(1)
!    if(lnod%seqno.ne.resultlist%seqno) then
!       call saveequil(21,form,lnod%ceqlink)
!       lnod=>lnod%next
!       goto 810
!    endif
    close(21)

1000 continue
    return
  end subroutine sstep1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine add_to_steplista(node,ceq)
! This subrotine creates a linked list.  node is always the last created node
! It should be nullified att first call and node%next is then set to itself
! At second and later calls node%next is allocated and the previous node
! will point at the new node.  The new node will have its next pointer set
! to point at the previous node%next which is the first node created.
! Each node will have a pointer to the ceq record representing the results
! the node numbering is just a way to check all is OK
! IT IS NOT EASY TO HANDLE POINTER AND SUCH TO CREATE THIS LIST ...
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(ssm_node), pointer :: node,slink,nynod
!    TYPE(gtp_equilibrium_data), allocatable, target :: ceqsave
!    write(*,*)'in ssm_lista: ',ceq%tpval(1)
! create a copy of the whole ceq gtp_equilibrium_data for this equilibrium
! possibly only some parts should be copied
! allocate an new ssm_node to point at ceqsave
    if(.not.associated(node)) then
       allocate(node)
       node%seqno=1
       node%tval=ceq%tpval(1)
       node%next=>node
    else
       slink=>node%next
       allocate(node%next)
       nynod=>node%next
       nynod%seqno=node%seqno+1
       nynod%tval=ceq%tpval(1)
       nynod%next=>slink
       node=>node%next
    endif
    allocate(node%ceqlink)
! does this copy ceq?
    node%ceqlink=ceq
! cannot find size of structure this way ...
!    k=sizeof(ceq)
!    write(*,*)'size: ',k
!    write(*,*)'node: ',node%seqno,node%tval,node%next%seqno
1000 continue
    return
  end subroutine add_to_steplista

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine splot1(depax,filename,resultlist,form)
! depax is text with dependent axis
! write plot state variable depax using gnuplot
! file name is intermediary file (maybe not needed)
! resultlist is circular list with calculated equilibria
! form is type of output (screen or postscript)
! ceq is current equilibria
    implicit double precision (a-h,o-z)
    TYPE(gtp_equilibrium_data), pointer :: ceq,curceq
    TYPE(ssm_node), pointer :: resultlist, lnod
    character depax*(*),ch1*1,gnuplotline*256,filename*(*),pfd*64,pfc*64
    character form*(*),pfh*64
    double precision, dimension(:,:), allocatable :: anp
    double precision, dimension(:), allocatable :: xax,yyy
    integer, dimension(:), allocatable :: nonzero
    character statevar*64,encoded*64
!\end{verbatim}
!
    write(*,*)'Entering OC plot version 1: ',depax(1:len_trim(depax)),gx%bmperr
    nrv=resultlist%seqno
! debug check of results present
!    write(*,*)'Number of records: ',nrv
!    lnod=>resultlist%next
!10  continue
!    write(*,*)'loop: ',lnod%seqno,lnod%tval,lnod%ceqlink%tpval(1)
!    if(lnod%seqno.ne.resultlist%seqno) then
!       lnod=>lnod%next
!       goto 10
!    endif
!
    allocate(xax(nrv))
    lnod=>resultlist%next
    nv=0
    xmin=1.0D20
    xmax=-1.0D20
    ymin=1.0D20
    ymax=-1.0D20
    wildcard: if(index(depax,'*').gt.0) then
! we do not know how many columns needed ...
! allocate as many array elements as columns
       np=1000
       allocate(anp(np,nrv))
       allocate(nonzero(np))
       allocate(yyy(np))
       nonzero=0
       write(*,*)'allocated: ',np,nrv
!------------------------------------------- begin loop 100
100    continue
! loop for all calculated equilibra, phases and composition sets
          nv=nv+1
          curceq=>lnod%ceqlink
          xax(nv)=curceq%tpval(1)
!          write(*,*)'T=',xax(nv)
          jj=0
          if(xax(nv).lt.xmin) xmin=xax(nv)
          if(xax(nv).gt.xmax) xmax=xax(nv)
          statevar=depax
! nzp should be dimension of yyy, np returns the number of values in yyy
          nzp=1000
          call get_many_svar(statevar,yyy,nzp,np,encoded,curceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'yaxis error: "',statevar(1:20),'"'
             goto 1000
!          else
!             write(*,*)'get_many: ',np,yyy(1)
          endif
!          write(*,16)'values: ',np,gx%bmperr,(yyy(i),i=1,np)
16        format(a,i3,i5,6(1pe12.4))
          if(gx%bmperr.ne.0) goto 1000
          do jj=1,np
             anp(jj,nv)=yyy(jj)
             if(abs(yyy(jj)).gt.zero) nonzero(jj)=1
          enddo
          if(lnod%seqno.ne.resultlist%seqno) then
             lnod=>lnod%next
             goto 100
          endif
!------------------------------------------- end loop 100
! end extracting data
! remove columns that are only zeroes
       ic=0
       nnp=np
!------------------------------------------ begin loop 150
150    ic=ic+1
160       continue
          if(ic.gt.nnp) goto 190
          if(nonzero(ic).eq.0) then
! shift all values from ic+1 to np
             do jj=ic,nnp-1
                do nnv=1,nv
                   anp(jj,nnv)=anp(jj+1,nnv)
                enddo
                nonzero(jj)=nonzero(jj+1)
             enddo
             nnp=nnp-1
             goto 160
          endif
          goto 150
!------------------------------------------ end loop 150
190 continue
       write(*,*)'reduced columns: ',np,nnp
       np=nnp
!======================================================================
    elseif(depax(1:3).eq.'cp ') then ! special single phase cp
       statevar=depax
       np=1000
       allocate(anp(1,nrv))
       curceq=>lnod%ceqlink
! find a stable phase
       do lokcs=1,noph()
!          if(curceq%phase_varres(lokcs)%amount(1).gt.zero) exit
          if(curceq%phase_varres(lokcs)%amfu.gt.zero) exit
       enddo
       write(*,*)'plot: ',lokcs
       nv=0
200    continue
          nv=nv+1
          xax(nv)=curceq%tpval(1)
          jj=0
          if(xax(nv).lt.xmin) xmin=xax(nv)
          if(xax(nv).gt.xmax) xmax=xax(nv)
!          call get_state_var_value(statevar,value,encoded,curceq)
!          if(gx%bmperr.ne.0) goto 1000
!          value=curceq%phase_varres(lokcs)%gval(1,1)
!          anp(1,nv)=value*globaldata%rgas*xax(nv)
          value=curceq%phase_varres(lokcs)%gval(4,1)
!          anp(1,nv)=xax(nv)*value*globaldata%rgas*xax(nv)
          anp(1,nv)=-xax(nv)*value*globaldata%rgas*xax(nv)
          write(*,211)'splot: ',nv,depax(1:10),nv,curceq%tpval(1),anp(1,nv)
211       format(a,i4,2x,a,i4,2(1pe12.4))
          if(lnod%seqno.ne.resultlist%seqno) then
             lnod=>lnod%next
             curceq=>lnod%ceqlink
             goto 200
          endif
! finished extracting all values for statevar
       np=1
!======================================================================
    else ! state variable without wildcard
       statevar=depax
       allocate(anp(1,nrv))
280    continue
! loop for all calculated equilibra, phases and composition sets
          nv=nv+1
          curceq=>lnod%ceqlink
          xax(nv)=curceq%tpval(1)
          jj=0
          if(xax(nv).lt.xmin) xmin=xax(nv)
          if(xax(nv).gt.xmax) xmax=xax(nv)
!          write(*,*)'splot ',statevar(1:20),curceq%eqno,curceq%tpval(1)
          call get_state_var_value(statevar,value,encoded,curceq)
          if(gx%bmperr.ne.0) goto 1000
          anp(1,nv)=value
          write(*,290)'splot: ',depax(1:10),nv,curceq%tpval(1),value
290       format(a,a,i4,2(1pe12.4))
          if(lnod%seqno.ne.resultlist%seqno) then
             lnod=>lnod%next
             goto 280
          endif
! finished extracting all values for statevar
       np=1
    endif wildcard
!--------------------------------------------------------------
    if(np.eq.0) then
       write(kou,*)'No data to plot'
       gx%bmperr=7777
       goto 1000
    endif
! gnuplot output
300 continue
    if(np+1.gt.1000) then
       write(*,*)'More than 1000 data points ...'
       goto 1000
    endif
    kk=len_trim(filename)
    pfd=filename(1:kk)//'.'//'dat '
!    write(*,*)'Gnuplot output on ',pdf(1:kk+4)
    open(21,file=pfd,access='sequential',status='unknown')
    write(21,310)(i,i=1,np+1)
310 format('# gnuplot output'/'#',i7,(1000i14))
    do nv=1,nrv
       write(21,320)xax(nv),(anp(jj,nv),jj=1,np)
320    format(f8.2,1000(1pe14.6))
    enddo
    write(21,330)
330 format('# plot "ocg.dat" using 1:2 with lines,',&
         ' "ocg.dat" using 1:3 with lines, ...'/'# set term postscript'/&
         '# set output "ocg.ps"'/'# plot ... '/'# ps2pdf ocg.ps')
    close(21)
    write(*,*)'Gnuplot data file   : ',pfd(1:kk+4)
! write the gnuplot command file
    kk=len_trim(filename)
    pfc=filename(1:kk)//'.'//'plt '
    kkk=kk+4
    open(21,file=pfc,access='sequential',status='unknown')
    if(form(1:1).eq.'P') then
       pfh=filename(1:kk)//'.'//'ps '
       write(21,366)pfh(1:len_trim(pfh))
366    format('set terminal postscript'/'set output "',a,'"')
    endif
    write(21,370)depax(1:len_trim(depax)),(pfd(1:kkk),i,i=2,np+1)
370 format('set title "Open Calphad with Gnuplot"'/&
         'set ylabel "',a,'"'/&
         'set xlabel "Temperature (K)"'/&
         'plot "',a,'" using 1:',i3,' with lines,',&
         999('"',a,'" using 1:',i3,' with lines,'))
    if(form(1:1).eq.' ') then
! if not hardcopy pause gnuplot
       write(21,371)
371    format('pause -1')
    endif
    close(21)
    gnuplotline='gnuplot '//pfc(1:kkk)//' '
    k3=len_trim(gnuplotline)+1
    write(*,*)'Gnuplot command file: ',pfc(1:kk+4)
    if(form(1:1).eq.'P') then
       write(*,*)'Gnuplot postscript file: ',pfh(1:kk+4)
!       write(*,*)'Gnuplot hardcopy file  : ',pfh(1:kk+4)
    endif
    call system(gnuplotline(1:k3))
1000 continue
    return
  end subroutine splot1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

END MODULE smp


