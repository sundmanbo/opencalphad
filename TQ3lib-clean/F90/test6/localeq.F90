!
! Test program to calculate a single equilibrium and extract some values
! Can be useful for a reactor simulator
!
program simple
!
  implicit none
  integer, parameter :: maxel=100
  integer nsel,ip,ierr,localequil,ipos,jpos,ie
  character selel(maxel)*2,ch1*1
  double precision tp(2),nmz(maxel)
  character todo*60,outspec*60
  double precision, dimension(:), allocatable :: outval,nrest
! elements
  nsel=7
  selel(1)='AL' 
  selel(2)='CR'
  selel(3)='CU'
  selel(4)='FE'
  selel(5)='MG'
  selel(6)='SI'
  selel(7)='ZN' 
! set values of temperature and pressure
  tp(1)=850
  tp(2)=1.0D5
! set amount of elements in moles.  Here it is important to have the
! amounts in the same order as the elements above.
! using numbers they must refer to teh elements in alphabetical order
  nmz(1)=0.87d0
  nmz(2)=0.003d0
  nmz(3)=0.03d0
  nmz(4)=0.002d0
  nmz(5)=0.04d0
  nmz(6)=0.05d0
  nmz(7)=0.005d0
! allocate space for results N(liquid,*), N(*), MU(*) and H
  allocate(outval(2*nsel+1))
  allocate(nrest(nsel))
!-------------------------------------
  todo='verbose open cost507r.tdb '
  outspec='n(liquid,*) n(*) mu(*) h'
!
  IERR=localequil(todo,nsel,selel,tp,nmz,outspec,outval)
  if(ierr.ne.0) then
     write(*,*)'Abnormal return from localequil ',ierr
     stop 1
  endif
!--------------------------------------------------------------------
! list results
  write(*,100)tp(1),tp(2),(nmz(ip),ip=1,nsel)
! amounts in liquid
  write(*,101)'N(liq,*) ',(outval(ip),ip=1,nsel)
! you have to calculate the solids yourself as N(*) - N(liquid,*)
  jpos=nsel
  do ie=1,nsel
     nrest(ie)=outval(jpos+ie)-outval(ie)
  enddo
  write(*,101)'N(sol,*) ',(nrest(ip),ip=1,nsel)
! chemical poterntails
  ipos=jpos+nsel
  write(*,102)'MU: ',(outval(ipos+ip),ip=1,nsel)
  ipos=ipos+nsel
  write(*,103)outval(ipos+1)
100 format('Input: ',F8.2,1pe14.6(/9x,0pF8.4,8F8.4))
101 format(a,(9F8.4))
102 format(a,(6e12.4))
103 format('Enthalpy: ',1pe14.6)
!
!===============================================================
! calculating again using same system at a higher T
  write(*,201)
201 format(/'>>>>>>>>>>>>>> Press RETURN to perform the second calculation')
  read(*,202)ch1
202 format(a)
  todo='NOVERBOSE '
  tp(1)=800.0D0
  IERR=localequil(todo,nsel,selel,tp,nmz,outspec,outval)
  if(ierr.ne.0) then
     write(*,*)'Abnormal return from localequil ',ierr
     stop 2
  endif
!---------------------------------------------------------------------
! list results
  write(*,100)tp(1),tp(2),(nmz(ip),ip=1,nsel)
! amounts in liquid
  write(*,101)'N(liq,*) ',(outval(ip),ip=1,nsel)
! you have to calculate the solids yourself as N(*) - N(liquid,*)
  jpos=nsel
  do ie=1,nsel
     nrest(ie)=outval(jpos+ie)-outval(ie)
  enddo
  write(*,101)'N(sol,*) ',(nrest(ip),ip=1,nsel)
  ipos=jpos+nsel
  write(*,102)'MU: ',(outval(ipos+ip),ip=1,nsel)
  ipos=ipos+nsel+1
  write(*,103)outval(ipos)
!------------------------------------
  write(*,*)'Finished normally'
end program simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function localequil(todo,nsel,selel,tp,nmz,outspec,outval)
!
! This makes all data structures and routines in OC available
  use liboctq
!
  implicit none
! maxel and maxph defined in pmod package
! zero and one already defined in OC
! double precision, parameter :: zero=0.0,one=1.0
  integer n,n1,n2,n3,nrph,nsel,mode,iph,ics,lent,loop,ip,jp,ie,ierr
  integer cnum(maxel+2)
  character filename*60,todo*(*),outspec*(*)
  character condition*60,statevar*60,line*80
  character target*60,phcsname*36,selel(nsel)*2
  double precision value,tp(2),nmz(maxel),outval(*)
  double precision nval(maxel)
  logical once
! inital value of verbose
  logical, save :: verbose=.FALSE.
! this is a pointer to the data structure for the equilibrium.
  type(gtp_equilibrium_data), save, pointer :: ceq
!
!  write(*,*)'todo: ',trim(todo)
! -----------------------------------------------------------
! execute commands in todo
  ip=0
  loop=0
10 continue
! extract command in todo, eolch skips until next nonspace, if empty go to 100
!  if(eolch(todo,ip)) goto 100
! 3rd argument 2 means getext will extract text up to next space character
  call getext(todo,ip,2,statevar,' ',lent)
! if empty jump to 100
  if(lent.le.1) goto 100
! set UPPER CASE
  call capson(statevar)
!
  select case(statevar)
  case default
     write(*,*)'Unknown command: ',statevar,lent
     ierr=1; goto 1000
  case('NOCALCULATE ')
! this can maybe be used to extract more data from last calculation
     goto 500
  case('VERBOSE ')
! turn on verbose mode
     verbose=.true.
  case('NOVERBOSE ')
! turn on verbose mode
     verbose=.false.
  case('OPEN ')
! initiate OC and open database and read in data for elements
     call tqini(n,ceq)
! ALWAYS CHECK FOR ERROR CODE AFTER CALL TO OC ROUTINE
     if(gx%bmperr.ne.0) goto 1000
     call getext(todo,ip,2,filename,' ',lent)
     if(lent.le.1) then
        write(*,*)'No file name after OPEN'
        ierr=2; goto 1000
     endif
     if(verbose) then
        write(*,*)'Reading the database: ',trim(filename)
     else
! turn off warnings from reading the TDB file
        call readtdbsilent
     endif
     call tqrpfil(filename,nsel,selel,ceq)
     if(gx%bmperr.ne.0) then
        write(*,*)'Error reading database: ',trim(filename)
        ierr=3; goto 1000
     endif
     if(verbose) then
! List all phases read from the database
        call tqgnp(nrph,ceq)
        write(*,20)nrph
20      format(/'Number of phases read from database: ',i3)
        line=' '
        n1=1
        n3=0
        do n=1,nrph
           call tqgpn(n,phcsname,ceq)
           if(n1+len_trim(phcsname)+1.gt.len(line)) then
              write(*,21)line(1:len_trim(line))
21            format(a)
              line=phcsname(1:len_trim(phcsname))//','
              n1=len_trim(phcsname)+2
           else
              line(n1:)=phcsname(1:len_trim(phcsname))//','
              n1=len_trim(line)+2
           endif
           if(phcsname(1:7).eq.'STABLE ') then
              n3=n
           endif
        enddo
        write(*,21)line(1:len_trim(line)-2)
     endif
!------------------------------------------------------------------
! suspend the phase STABLE, its phase tuple index is n3
! status values: -3 suspend; -2 dormant; 
!                -1, 0, 1 entered and not stable/unknown/stable; 2 fixed
! the value relevant only for fix status
     value=zero
     if(n3.gt.0) call tqphtupsts(n3,-2,value,ceq)
     if(gx%bmperr.ne.0) goto 1000
! one can add more cases ... for example different types of conditions
  end select
  goto 10
!------------------------------------------------------------------
! always set conditions
100 continue
  if(verbose) then
     write(*,*)
     write(*,*)'T and P: ',tp(1),tp(2)
  endif
  n1=0
  n2=0
  condition='T'
  call tqsetc(condition,n1,n2,tp(1),cnum(1),ceq)
  if(gx%bmperr.ne.0) goto 1000
  condition='P'
  call tqsetc(condition,n1,n2,tp(2),cnum(2),ceq)
  if(gx%bmperr.ne.0) goto 1000
  condition='N'
  n2=0
  do n1=1,nsel
     call tqsetc(condition,n1,n2,nmz(n1),cnum(2+n1),ceq)
     if(gx%bmperr.ne.0) goto 1000
  enddo
! ----------------------------------------------------------
! if verbose use an internal routine to list the conditions ....
! lut or * is stadard output unit in Fortran
  if(verbose) then
     lut=6
     write(*,*)
     write(*,*)'Current set of conditions:'
     call list_conditions(6,ceq)
  endif
! ----------------------------------------------------------
! now calculate the equilibria
! n1=0 means call grid minimizer, n1=-1 do not call grid minimizer
  target=' '
  n1=0
  n2=0
  write(*,*)
  write(*,*)'Calculating the equilibrium'
  call tqce(target,n1,n2,value,ceq)
  if(gx%bmperr.ne.0) then
     write(*,*)'Calculation failed, error code: ',gx%bmperr
     ierr=gx%bmperr; goto 1000
  elseif(verbose) then
     write(*,*)'Successful calculation'
  endif
  if(verbose) then
! extra output
     write(lut,80)
80   format(/'Some global data: ')
     call list_global_results(lut,ceq)
! n1=1 means mole fractions
     n1=1
     write(lut,81)
81   format(/'Some component data: ')
     call list_components_result(lut,n1,ceq)
! ALWAYS CHECK FOR ERROR CODE AFTER CALL TO OC ROUTINE
     if(gx%bmperr.ne.0) then
        write(*,*)'Reset error 1: ',gx%bmperr
        gx%bmperr=0
     endif
! The phase results can be formatted in different ways with mode
! mode 1000 stable phases, only composition, mole fractions in value order
!      1100 stable phases, only composition, mole fractions in alpabetical ord
!      1110 stable phases, constitution, mole fractions in alpabetical ord
!      1011 stable phases, only composition, mole fractions in value order
! etc
     mode=1000
! once is used to write the heading once
     once=.TRUE.
     write(lut,82)
82   format(/'Some phase data: ')
! loop for all phases
     do iph=1,noph()
! loop for all composition sets of the phase
        do ics=1,noofcs(iph)
! this list stable phases only  
! once TRUE means heading written.  Automatically set to FALSE
           call list_phase_results(iph,ics,mode,lut,once,ceq)
        enddo
     enddo
     if(gx%bmperr.ne.0) then
        write(*,*)'Reset error 2: ',gx%bmperr
        gx%bmperr=0
     endif
  endif
! ----------------------------------------------------
! now handle requested output of results
500 continue
  if(verbose) write(*,*)'Output: ',trim(outspec)
! extract result according to outspec
  ip=0
  loop=0
! return here for more output, jump to 900 when finished
510 continue
!  if(eolch(outspec,ip)) goto 900
  call getext(outspec,ip,2,statevar,' ',lent)
  if(lent.le.1 .and. statevar(1:1).eq.' ') goto 900
  call capson(statevar)
  if(verbose) write(*,*)'Output variable: ',trim(statevar),ip
  select case(statevar)
!--------------------------------------------
  case default
     write(*,*)'Cannot extract unknown state variable: ',trim(statevar)
     ierr=10; goto 1000
! one may merge these N(*) cases to simplify the code ...
!--------------------------------------------
  case('N(*) ')
     n3=maxel
     call tqgetv(statevar,-1,0,n3,nval,ceq)
     if(gx%bmperr.ne.0) then
        write(*,*)'Error extracting ',trim(statevar)
        ierr=gx%bmperr; goto 1000
     endif
     if(verbose) then
        write(*,*)'Result ',trim(statevar)
        write(*,520)(nval(ie),ie=1,n3)
520     format(7(1pe11.3))
     endif
     do ie=1,nsel
        outval(loop+ie)=nval(ie)
     enddo
     loop=loop+nsel
!--------------------------------------------
  case('N(GAS,*) ')
     n3=maxel
     call tqgetv(statevar,-1,0,n3,nval,ceq)
     if(gx%bmperr.ne.0) then
        write(*,*)'Error extracting ',trim(statevar)
        ierr=gx%bmperr; goto 1000
     endif
     if(verbose) then
        write(*,*)'Result ',trim(statevar)
        write(*,520)(nval(ie),ie=1,n3)
     endif
     do ie=1,nsel
        outval(loop+ie)=nval(ie)
     enddo
     loop=loop+nsel
!--------------------------------------------
  case('N(LIQUID,*) ')
     n3=maxel
     call tqgetv(statevar,-1,0,n3,nval,ceq)
     if(gx%bmperr.ne.0) then
        write(*,*)'Error extracting ',trim(statevar)
        ierr=gx%bmperr; goto 1000
     endif
     if(verbose) then
        write(*,*)'Result ',trim(statevar)
        write(*,520)(nval(ie),ie=1,n3)
     endif
     do ie=1,nsel
        outval(loop+ie)=nval(ie)
     enddo
!     loop=loop+nsel
! we may have several liquids ...
     jp=index(statevar,',')
     statevar(jp:)='#2,*) '
!     write(*,*)'A second liquid? ',trim(statevar)
     n3=maxel
     call tqgetv(statevar,-1,0,n3,nval,ceq)
     if(gx%bmperr.ne.0) then
! this error OK if there is no liquid#2, just ignore
        gx%bmperr=0
     else
! otherwise add together
        do ie=1,nsel
           outval(loop+ie)=outval(loop+ie)+nval(ie)
        enddo
     endif
     loop=loop+nsel
!--------------------------------------------
  case('MU(*) ')
     n3=maxel
     call tqgetv(statevar,-1,0,n3,nval,ceq)
     if(gx%bmperr.ne.0) then
        write(*,*)'Error extracting ',trim(statevar)
        ierr=gx%bmperr; goto 1000
     endif
     if(verbose) then
        write(*,*)'Result ',trim(statevar)
        write(*,520)(nval(ie),ie=1,n3)
     endif
     do ie=1,nsel
        outval(loop+ie)=nval(ie)
     enddo
     loop=loop+nsel
!--------------------------------------------
  case('H ')
     n3=maxel
     call tqgetv(statevar,0,0,n3,nval,ceq)
     if(verbose) then
        write(*,*)'Result ',trim(statevar),loop+1
        write(*,520)(nval(ie),ie=1,n3)
     endif
     if(gx%bmperr.ne.0) then
        write(*,*)'Error extracting ',trim(statevar)
        ierr=gx%bmperr; goto 1000
     endif
     outval(loop+1)=nval(1)
     loop=loop+1
  end select
! look for another output variable
  goto 510
!
900 continue
  localequil=0
  return
! error
1000 continue
  localequil=ierr
  return
end function localequil

