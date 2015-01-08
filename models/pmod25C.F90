!
! included in pmod25.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     7. State variable manipulations
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine get_state_var_value(statevar,value,encoded,ceq)
! called with a state variable character
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   character statevar*(*),encoded*(*)
   double precision value
!\end{verbatim}
!   integer indices(4)
   integer iunit,ip,lrot,mode
   type(gtp_state_variable), pointer :: svr
   character actual_arg(2)*16,name*16
!
   iunit=0
   call decode_state_variable(statevar,svr,ceq)
!   write(*,20)statevar(1:len_trim(statevar)),svr%oldstv,svr%norm,&
!        svr%argtyp,svr%component
20  format('25C gsvv 1: ',a,' : ',4i3)
   if(gx%bmperr.ne.0) then
!      goto 1000
! it can be a state variable symbol ...
!
! Possible problem ... this can cause nesting as a state variable will
! normally evaluate some state variables or other state variable functions
!
      gx%bmperr=0
      name=statevar
      call capson(name)
      call find_svfun(name,lrot,ceq)
      if(gx%bmperr.ne.0) then
         write(*,*)'Neither state variable or symbol'
         gx%bmperr=8888; goto 1000
      else
! get the value of the symbol, may involve other symbols and state variablse
! The actual_arg is a facility not yet implemented and not allowed here
! if mode=0 the stored value may be used, mode=1 always evaluate
         actual_arg=' '
         mode=1
! this is OK as no derivative
         value=evaluate_svfun_old(lrot,actual_arg,mode,ceq)
         encoded=name
      endif
   else
! it is a real state variable
      call state_variable_val(svr,value,ceq)
      if(gx%bmperr.ne.0) goto 1000
      ip=1
      call encode_state_variable(encoded,ip,svr,ceq)
      if(gx%bmperr.ne.0) then
         write(*,*)'encode error: ',gx%bmperr
         gx%bmperr=0; encoded='dummy'
      endif
   endif
1000 continue
   return
 end subroutine get_state_var_value

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine get_many_svar(statevar,values,mjj,kjj,encoded,ceq)
! called with a state variable name with woldcards allowed like NP(*), X(*,CR)
! mjj is dimension of values, kjj is number of values returned
! encoded used to specify if phase data in phasetuple order ('Z')
! >>>> BIG problem: How to do with phases that are note stable?
! If I ask for w(*,Cr) I only want the fraction in stable phases
! but whenthis is used for GNUPLOT the values are written in a matix
! and the same column in that phase must be the same phase ...
! so I have to have the same number of phases from each equilibria.
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   character statevar*(*),encoded*(*)
   double precision values(*)
   integer mjj,kjj
!\end{verbatim}
   integer indices(4),modind(4)
   double precision xnan,xxx
   integer jj,lokph,lokcs,k1,k2,k3,iref,jl,iunit,istv
   type(gtp_state_variable), pointer :: svr
!   logical phtupord
! calculate the NaN bit pattern
   xnan=0.0d0
!   xnan=0.0d0/xnan
   if(gx%bmperr.ne.0) then
      write(*,*)'Error entering get_many_svar ',gx%bmperr,xnan
   endif
!------------------------
   iunit=0
   modind=0
!   phtupord=.FALSE.
!   if(encoded(1:1).eq.'Z') then
! when called from TQ interface the phase order should be as for phase tuples
!      phtupord=.TRUE.
!   endif
! called from minimizer for testing
!   write(*,*)'gmv 1: ',statevar(1:20)
!   call decode_state_variable(statevar,istv,indices,iref,iunit,svr,ceq)
   call decode_state_variable(statevar,svr,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'Failed decode statevar in get_many_svar',gx%bmperr
      goto 1000
   endif
! translate svr data to old indices etc
   istv=svr%oldstv
   iref=svr%phref
   iunit=svr%unit
! svr%argtyp specifies values in indices:
! svr%argtyp: 0=no arguments; 1=comp; 2=ph+cs; 3=ph+cs+comp; 4=ph+cs+const
   indices=0
   if(svr%argtyp.eq.1) then
      indices(1)=svr%component
   elseif(svr%argtyp.eq.2) then
      indices(1)=svr%phase
      indices(2)=svr%compset
   elseif(svr%argtyp.eq.3) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%component
   elseif(svr%argtyp.eq.4) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%constituent
!   else
!      write(*,*)'state variable has illegal argtyp: ',svr%argtyp
!      gx%bmperr=7775; goto 1000
   endif
!
!   write(*,20)istv,indices,iref,gx%bmperr
20 format('gmsvar 1: ',i5,4i4,3i7)
! -----------------------------------------
! Indices 1: one or all components (-1)
!          Indices 2+3: 0 or phase+set 
! Indices 1+2: phase+set
!          Indices 3: 0 or component (-1) or constituent (-2)
! indices 4 never used
! -----------------------------------------
! -1 means element or component
! -2 species or constituent
! -3 phase
! -4 composition set
   jj=0
   if(indices(1).ge.0) then
      if(indices(2).ge.0) then
         if(indices(3).ge.0) then
! all indices given, a single value
            jj=jj+1
            if(jj.gt.mjj) goto 1100
            call state_variable_val3(istv,indices,iref,&
                 iunit,values(jj),ceq)
            if(gx%bmperr.ne.0) goto 1000
         elseif(indices(3).eq.-1) then
! loop for components, indices 1+2 must be phase+compset
            do k3=1,noofel
               indices(3)=k3
               jj=jj+1
               if(jj.gt.mjj) goto 1100
               call state_variable_val3(istv,indices,iref,&
                    iunit,values(jj),ceq)
               if(gx%bmperr.ne.0) goto 1000
            enddo
         elseif(indices(3).eq.-2) then
! loop for constituents, indices 1+2 must be phase+compset
            call get_phase_record(indices(1),lokph)
            do k3=1,phlista(lokph)%tnooffr
               indices(3)=k3
               jj=jj+1
               if(jj.gt.mjj) goto 1100
               call state_variable_val3(istv,indices,iref,&
                    iunit,values(jj),ceq)
               if(gx%bmperr.ne.0) goto 1000
            enddo
         else
! indices(3) must be -2, -1 or >=0 so if we are here there is an error
            write(*,17)'Illegal set of indices 1',(indices(jl),jl=1,4)
17          format(a,4i4)
            gx%bmperr=7777; goto 1000
         endif
      elseif(indices(2).eq.-3) then
! if indices(1)>=0 then indices(2)<0 must means a loop for all phase+compset
         do k2=1,noofph
            indices(2)=k2
            call get_phase_record(indices(2),lokph)
            do k3=1,phlista(lokph)%noofcs
               indices(3)=k3
               jj=jj+1
               if(jj.gt.mjj) goto 1100
               call get_phase_compset(indices(2),indices(3),lokph,lokcs)
! if composition set not stable so return NaN
!               if(test_phase_status(indices(2),indices(3),xxx,ceq).ge.2) then
               if(test_phase_status(indices(2),indices(3),xxx,ceq).le. &
                    PHENTUNST) then
                  values(jj)=xnan
!               elseif(.not.btest(ceq%phase_varres(lokcs)%status2,&
!                    CSSTABLE)) then
!                  values(jj)=xnan
               elseif(ceq%phase_varres(lokcs)%dgm.lt.zero) then
! the phase must not have negative driving force
                  values(jj)=xnan
               else
! problem that get_many returns values for unstable phases
                  call state_variable_val3(istv,indices,iref,&
                       iunit,values(jj),ceq)
                  if(gx%bmperr.ne.0) goto 1000
!                  write(*,23)'25C many 1: ',indices,values(jj),&
!                       ceq%phase_varres(lokcs)%dgm
23                format(a,2i3,2(1pe14.6))
               endif
            enddo
         enddo
      else
! if indices(1)>=0 then indices(2) must be -3 or >=0, so if here it is error
         write(*,17)'Illegal set of indices 2',(indices(jl),jl=1,4)
         gx%bmperr=7777; goto 1000
      endif
   elseif(indices(1).eq.-1) then
! loop for component as first indices, 2+3 can be fix phase+compset
      if(indices(2).ge.0) then
         do k1=1,noofel
            indices(1)=k1
            jj=jj+1
            if(jj.gt.mjj) goto 1100
            call state_variable_val3(istv,indices,iref,&
                 iunit,values(jj),ceq)
            if(gx%bmperr.ne.0) goto 1000
         enddo
      elseif(indices(2).eq.-3) then
! loop for components and phase+compset
         do k1=1,noofel
            indices(1)=k1
            do k2=1,noofph
               indices(2)=k2
               call get_phase_record(indices(2),lokph)
               do k3=1,phlista(lokph)%noofcs
                  indices(3)=k3
                  jj=jj+1
                  if(jj.gt.mjj) goto 1100
                  call get_phase_compset(indices(2),indices(3),lokph,lokcs)
! if composition not stable so return NaN
                  if(test_phase_status(indices(2),indices(3),xxx,ceq).le. &
                       PHENTSTAB) then
                     values(jj)=xnan
!                  elseif(.not.btest(ceq%phase_varres(lokcs)%status2,&
!                       CSSTABLE)) then
!                     values(jj)=xnan
                  elseif(ceq%phase_varres(lokcs)%dgm.lt.zero) then
! the phase must not have negative driving force
                     values(jj)=xnan
                  else
                     call state_variable_val3(istv,indices,iref,&
                          iunit,values(jj),ceq)
                     if(gx%bmperr.ne.0) goto 1000
!                     write(*,23)'25C many 2: ',indices(1),indices(2),&
!                          values(jj),ceq%phase_varres(lokcs)%dgm
                  endif
               enddo
            enddo
         enddo
      else
! if we come here it must be an error
         write(*,17)'Illegal set of indices 3',(indices(jl),jl=1,4)
         gx%bmperr=7777; goto 1000
      endif
   elseif(indices(1).eq.-3) then
! loop for phase+compset as indices(1+2)
! here we must be careful not to destroy original indices, use modind
!      write(*,*)'get_many NP(*) 1: ',gx%bmperr,indices(3)
!      write(*,*)'Loop for many phases',indices(1)
      do k1=1,noofph
         modind(1)=k1
         modind(2)=0
         call get_phase_record(modind(1),lokph)
!         write(*,19)'test 17',modind,gx%bmperr,xnan
         if(gx%bmperr.ne.0) goto 1000
         do k2=1,phlista(lokph)%noofcs
            modind(2)=k2
            jj=jj+1
            if(jj.gt.mjj) goto 1100
            call get_phase_compset(modind(1),modind(2),lokph,lokcs)
!            write(*,19)'test 2: ',modind,gx%bmperr,xnan
19          format(a,4i3,i7,1pe12.4)
            if(gx%bmperr.ne.0) goto 1000
            call get_phase_compset(modind(1),modind(2),lokph,lokcs)
            if(gx%bmperr.ne.0) then
!               write(*,19)'error 2',modind,gx%bmperr
               goto 1000
            endif
            if(test_phase_status(modind(1),modind(2),xxx,ceq).le. &
                 PHENTUNST) then
! if phase not entered or fix return NaN
               values(jj)=xnan
            else
               if(indices(3).eq.0) then
! This is typically listing of NP(*) for all phases
                  modind(3)=indices(3)
!                  write(*,16)16,(modind(i),i=1,4),gx%bmperr
!16                format('bug: ',i3,4i5,i7)
                  call state_variable_val3(istv,modind,iref,&
                       iunit,values(jj),ceq)
!                  write(*,*)'get_many NP(*): 2',gx%bmperr
                  if(gx%bmperr.ne.0) goto 1000
               elseif(ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) then
!               elseif(.not.btest(ceq%phase_varres(lokcs)%status2,&
!                    CSSTABLE)) then
! if wildcard for index 3 the phase must be stable
!                  write(*,*)'Not stable: ',modind(1),modind(2),&
!                       ceq%phase_varres(lokcs)%phstate
                  values(jj)=xnan
               elseif(indices(3).gt.0) then
! This is typically listing of w(*,cr), only in stable range of phases
                  modind(3)=indices(3)
!                  write(*,16)16,(modind(i),i=1,4),gx%bmperr
!16                format('bug: ',i3,4i5,i7)
                  call get_phase_compset(modind(1),modind(2),lokph,lokcs)
                  if(gx%bmperr.ne.0) goto 1000
!                  write(*,23)'25C many 3: ',modind(1),modind(2),values(jj),&
!                       ceq%phase_varres(lokcs)%dgm
                  if(ceq%phase_varres(lokcs)%dgm.lt.zero) then
! the phase must not have negative driving force
                     values(jj)=xnan
                  else
                     call state_variable_val3(istv,modind,iref,&
                          iunit,values(jj),ceq)
                  endif
                  if(gx%bmperr.ne.0) goto 1000
               elseif(indices(3).eq.-1) then
! loop for components of all phases
                  do k3=1,noofel
                     modind(3)=k3
                     call state_variable_val3(istv,modind,iref,&
                          iunit,values(jj),ceq)
                     if(gx%bmperr.ne.0) goto 1000
                  enddo
               elseif(indices(3).eq.-2) then
! loop for constituents of all phases
                  do k3=1,phlista(lokph)%tnooffr
                     modind(3)=k3
                     call state_variable_val3(istv,modind,iref,&
                          iunit,values(jj),ceq)
                     if(gx%bmperr.ne.0) goto 1000
                  enddo
               else
! error if here
                  write(*,17)'Illegal set of indices 4',(indices(jl),jl=1,4)
                  gx%bmperr=7777; goto 1000
               endif
               if(gx%bmperr.ne.0) then
                  write(*,19)'error 3',modind,gx%bmperr
                  goto 1000
               endif
            endif
         enddo
      enddo
   else
! error if here
      write(*,17)'Illegal set of indices 5',(indices(jl),jl=1,4)
      gx%bmperr=7777; goto 1000
   endif
!   ip=1
!   call encode_state_variable(encoded,ip,istv,indices,iunit,iref,ceq)
!   if(gx%bmperr.ne.0) then
!      write(*,*)'encode error: ',gx%bmperr
!      gx%bmperr=0; encoded='dummy'
!   endif
1000 continue
   kjj=jj
   return
1100 continue
   write(*,*)'Overflow in array to get_state_variables'
   gx%bmperr=7777; goto 1000
 end subroutine get_many_svar

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine decode_state_variable(statevar,svr,ceq)
! converts a state variable character to state variable record
   character statevar*(*)
   type(gtp_state_variable), pointer :: svr
   type(gtp_equilibrium_data), pointer :: ceq
! this subroutine using state variable records is a front end of the next:
!\end{verbatim} %+
!   type(gtp_state_variable) :: svrec   
   integer istv,indices(4),iref,iunit
   call decode_state_variable3(statevar,istv,indices,iref,iunit,svr,ceq)
1000 continue
   return
 end subroutine decode_state_variable

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine decode_state_variable3(statevar,istv,indices,iref,iunit,svr,ceq)
! converts an old state variable character to indices 
! Typically: T, x(fe), x(fcc,fe), np(fcc), y(fcc,c#2), ac(h2,bcc), ac(fe)
! NOTE! model properties like TC(FCC),MQ&FE(FCC,CR) must be detected
! NOTE: added storing information in a gtp_state_variable record svrec !!
!
! this routine became as messy as I tried to avoid
! but I leave it to someone else to clean it up ...
!
! state variable and indices
! Symbol  no   index1 index2 index3 index4  
! T       1    -
! P       2    -
! MU      3    component or phase,constituent
! AC      4    component or phase,constituent
! LNAC    5    component or phase,constituent
!                                          index (in svid array)
! U       10   (phase#set)                    6     Internal energy (J)
! UM      11    "                             6     per mole components
! UW      12    "                             6     per kg
! UV      13    "                             6     per m3
! UF      14    "                             6     per formula unit
! S       2x    "                             7     entropy
! V       3x    "                             8     volume
! H       4x    "                             9     enthalpy
! A       5x    "                            10     Helmholtz energy
! G       6x    "                            11     Gibbs energy
! NP      7x    "                            12     moles of phase
! BP      8x    "                            13     mass of moles
! DG      9x    "                            15     Driving force
! Q       10x   "                            14     Internal stability
! N       11x (component/phase#set,component) 16  moles of components
! X       111   "                            17     mole fraction of components
! B       12x   "                            18     mass of components
! W       122   "                            19     mass fraction of components
! Y       13    phase#set,constituent#subl   20     constituent fraction
!----- model variables <<<< these now treated differently
! TC      -     phase#set                    -      Magnetic ordering T
! BMAG    -     phase#set                    -      Aver. Bohr magneton number
! MQ&     -     element, phase#set           -      Mobility
! THET    -     phase#set                    -      Debye temperature
!
   implicit none
   integer, parameter :: noos=20
   character*4, dimension(noos), parameter :: svid = &
       ['T   ','P   ','MU  ','AC  ','LNAC','U   ','S   ','V   ',&
        'H   ','A   ','G   ','NP  ','BP  ','DG  ','Q   ','N   ',&
        'X   ','B   ','W   ','Y   ']
!        1      2      3      4      4      6      7      8         
   character statevar*(*)
   integer istv,iref,iunit
   integer, dimension(4) :: indices
   type(gtp_equilibrium_data), pointer :: ceq
! I shall try to use this record type instead of separate arguments: !!
!   type(gtp_state_variable), pointer :: svrec
   type(gtp_state_variable), pointer :: svr
!\end{verbatim}
!   type(gtp_state_variable), allocatable, target :: svr
   integer is,jp,kp,iph,ics,icon,icomp,norm,narg,icc
   double precision cmass,asum
!
   character argument*60,arg1*24,arg2*24,ch1*1,lstate*60,propsym*60
   integer typty
   logical deblist
! initiate svr internal variables
   deblist=.FALSE.
!   deblist=.TRUE.
   if(ocv()) deblist=.TRUE.
   if(deblist) write(*,*)'25C entering decode_statevariable: ',&
        statevar(1:len_trim(statevar))
!   write(*,*)'25C svr allocated'
   allocate(svr)
!   write(*,*)'25C svr assignment start'
   svr%oldstv=0
   svr%norm=0
   svr%unit=0
! svr%argtyp: 0=no arguments; 1=comp; 2=ph+cs; 3=ph+cs+comp; 4=ph+cs+const
   svr%argtyp=0
   svr%phref=0
   svr%phase=0
   svr%compset=0
   svr%component=0
   svr%constituent=0
!   write(*,*)'25C svr assignment end'
!
! For wildcard argument "*" return:
! -1 for element or component
! -2 for species or constituent
! -3 for phase
! -4 for composition set
   istv=-1
   indices=0
   iref=0
   iunit=0
   iph=0
   ics=0
   norm=0
! local character for state variable
   lstate=statevar
   call capson(lstate)
   if(deblist) write(*,*)'25C decode_state_var 1: ',lstate(1:20)
! compare first character
   ch1=lstate(1:1)
   do is=1,noos
      if(ch1.eq.svid(is)(1:1)) goto 50
   enddo
! it may be a property, parameter identifier
   goto 600
!------------------------------------------------------------
50 continue
   if(deblist) write(*,*)'25C dsv 1: ',is,lstate(1:30)
   if(is.eq.1) then
      if(lstate(2:2).ne.' ') then
! it must be a property like TC or THET
         goto 600
      endif
! T
      istv=1; svr%oldstv=1; svr%statevarid=1; goto 1000
   elseif(is.eq.2) then
! P
      if(lstate(2:2).ne.' ') goto 600
      istv=2; svr%oldstv=2; svr%statevarid=2; goto 1000
   elseif(is.gt.5) then
      goto 100
   endif
!------------------------------------------------------------
! MU      3    component, possible suffix S for SER reference
   chemp: if(is.eq.3) then
      if(lstate(1:2).ne.'MU') then
         goto 600
      endif
      istv=3
      jp=3
   elseif(is.eq.4) then
! AC is 4 but just A or AM, AV etc can mean Helmholtz Energy or a property
      if(lstate(1:2).ne.'AC') then
         is=8; goto 100
      endif
      istv=4
      jp=3
 elseif(is.eq.5) then
! LNAC    5    component
      if(lstate(1:4).ne.'LNAC') goto 600
      istv=5
      jp=5
   endif chemp
! MU, AC and LNAC can have a suffix 'S', reference state, iref=0 is default
   if(lstate(jp:jp).eq.'S') then
! This iref has not been treated correctly so far.  The idea is now that
! iref=0 means user defined reference state, if the user has not defined any
! reference state it means SER.  If the user specifies a suffix S it means
! always SER even if the user has defined another reference state.
! Maybe iref>0 will have some other meaing in the future ...
      iref=-1
      jp=jp+1
   endif
! extract the argument, can be one or two indices
   svr%oldstv=istv; svr%statevarid=istv
   if(lstate(jp:jp).ne.'(') goto 1130
   kp=index(lstate,')')
   if(kp.lt.jp) goto 1140
   argument=lstate(jp+1:kp-1)
   kp=index(argument,',')
   if(kp.gt.0) then
! >>> if two arguments first is phase ??? different from TC
      arg1=argument(1:kp-1)
      arg2=argument(kp+1:)
      if(arg1(1:2).eq.'* ') then
         iph=-3
      else
         call find_phase_by_name(arg1,iph,ics)
         if(gx%bmperr.ne.0) goto 1150
      endif
      if(arg2(1:2).eq.'* ') then
         icon=-2
      else
         call find_constituent(iph,arg2,cmass,icon)
         if(gx%bmperr.ne.0) goto 1160
         call set_constituent_reference_state(iph,icon,asum)
         if(gx%bmperr.ne.0) then
            gx%bmperr=4112; goto 1000
         endif
      endif
! composition set irrelevant as chempot depend only on species stoichiometry
      indices(1)=iph
      indices(2)=icon
      svr%phase=iph
      svr%compset=1
      svr%constituent=icon
      svr%argtyp=4
   else
      if(argument(1:2).eq.'* ') then
         icomp=-1
      else
         call find_component_by_name(argument,icomp,ceq)
         if(gx%bmperr.ne.0) goto 1170
      endif
      indices(1)=icomp
      svr%component=icomp
      svr%argtyp=1
   endif
   goto 1000
!=================================================================
! extensive variable, is=6..20 or a model property
100 continue
   jp=2
! check second letter for some state variables
   if(deblist) write(*,105)is,norm,jp
105 format('25C dsv 4: ',3i4)
   letter2: if(is.eq.12 .and. lstate(jp:jp).ne.'P') then
! This is for Nx or a property
      is=16
   elseif(is.eq.13) then
! this can be Bx for component, BP for phase or BMAG for Bohr magnetons
      if(lstate(jp-1:jp).eq.'BP') then
         jp=jp+1
      else
! this is Bx or a property
         is=18
      endif
   elseif(is.eq.14 .and. lstate(jp-1:jp).ne.'DG') then
! this is for Dx, can be a property
!      gx%bmperr=4107; goto 1000
      goto 600
   elseif(is.eq.12 .or. is.eq.14) then
! This is NP or DG, increment jp to check the second character
      jp=jp+1
   elseif(is.eq.17 .or. is.eq.19) then
! X and W can have a suffix % to indicate percentage
      if(lstate(jp:jp).eq.'%') then
         iunit=100
         jp=jp+1
         svr%unit=iunit
      endif
   endif letter2
!---------------------------------------------------------------------
! If we come here the first (and sometimes second) letter must have been:
!               A,  B, BP,  D,  G, H,  N, NP,  Q, S, U,  W,  X,  Y
! and "is" is  10, 18, 13, 14, 11, 9, 16, 12, 15, 7, 6, 19, 17, 20
! NOTE: for N and B the second character has been checked and jp incremented
!       if equal to P.  The third (for NP and BP forth) character must 
!       be normallizing (MWVF), a space or a (, otherwise it is a property
   if(deblist) write(*,*)'25C lstate: ',lstate(1:20)
! these have no normalizing: Q, X, W, Y
   nomalize: if(is.le.14 .or. is.eq.16 .or. is.eq.18) then
! ZM      x1   (phase)                             per mole components
! ZW      x2   (phase)                             per kg
! ZV      x3   (phase)                             per m3
! ZF      x4   phase must be specified             per formula unit
      ch1=lstate(jp:jp)
      jp=jp+1
      if(ch1.eq.'M') then
         norm=1
      elseif(ch1.eq.'W') then
         norm=2
      elseif(ch1.eq.'V') then
         norm=3
      elseif(ch1.eq.'F') then
         norm=4
      else
! no or default normalization, backspace
         jp=jp-1
      endif
      svr%norm=norm
      if(deblist) write(*,*)'25C Normalize 1: ',is,jp,ch1,norm
   endif nomalize
!---------------------------------------------------------------------
! extract arguments if any. If arguments then lstate(jp:jp) should be (
! Typically G(fcc#2), N(Cr), BP(fcc), Y(sigma#2,cr#3), TC(BCC#2)
!300 continue
   if(deblist) write(*,*)'25C args: ',jp,lstate(1:jp+10)
   narg=0
   args: if(lstate(jp:jp).eq.'(') then
      kp=index(lstate,')')
      if(kp.le.0) then
         if(deblist) write(*,110)'25C dsv 5: ',is,jp,kp,lstate(1:20)
110      format(a,3i3,a)
         gx%bmperr=4103; goto 1000
      endif
      argument=lstate(jp+1:kp-1)
      kp=index(argument,',')
      arg: if(kp.gt.0) then
         arg1=argument(1:kp-1)
         arg2=argument(kp+1:)
         narg=2
         kp=index(arg2,',')
         if(kp.gt.0) then
! too many arguments to a state variable
            gx%bmperr=4097; goto 1000
         endif
      else !no arg
         narg=1
         arg1=argument
      endif arg
   elseif(lstate(jp:jp).ne.' ') then
! if additional character then it must be a property
      goto 600
   endif args
!------------------
! transform arguments to indices, different arguments for 6-
! Handle arguments: U, S, V, H, A, G, NP, BP,DG, Q, N, X, B, W, Y
!                   6, 7, 8, 9,10, 11,12, 13,14,15,16,17,18,19,20
   if(narg.eq.1) then
      if(is.le.15 .or. is.ge.21) then
! single argument is phase+composition set
         if(arg1(1:2).eq.'* ') then
            iph=-3
            ics=-4
         else
            call find_phase_by_name(arg1,iph,ics)
            if(gx%bmperr.ne.0) goto 1000
         endif
         indices(1)=iph
         indices(2)=ics
         svr%phase=iph
         svr%compset=ics
         svr%argtyp=2
      elseif(is.eq.20) then
! state variable Y must have 2 arguments
         gx%bmperr=4098; goto 1000
      else
! single argument is component for is=16-19
         if(arg1(1:2).eq.'* ') then
            icomp=-1
         else
            call find_component_by_name(arg1,icomp,ceq)
            if(gx%bmperr.ne.0) goto 1000
         endif
         indices(1)=icomp
         svr%component=icomp
         svr%argtyp=1
      endif
   elseif(narg.eq.2) then
! two arguments only for is=16-20, first phase, second component or constit
      if(is.le.15 .or. is.ge.21) then
         gx%bmperr=4110; goto 1000
      endif
      if(arg1(1:2).eq.'* ') then
         iph=-3
         ics=-4
      else
         call find_phase_by_name(arg1,iph,ics)
         if(gx%bmperr.ne.0) goto 1000
      endif
      indices(1)=iph
      indices(2)=ics
      svr%phase=iph
      svr%compset=ics
      if(is.eq.20) then
         if(arg2(1:2).eq.'* ') then
            icc=-2
         else
            call find_constituent(iph,arg2,cmass,icc)
            if(gx%bmperr.ne.0) goto 1000
         endif
         svr%constituent=icc
         svr%argtyp=4
      else
         if(arg2(1:2).eq.'* ') then
            icc=-1
         else
            call find_component_by_name(arg2,icc,ceq)
            if(gx%bmperr.ne.0) goto 1000
         endif
         svr%component=icc
         svr%argtyp=3
      endif
! note indices(4) never used as icc is constituent index, arg2 must have
! a #sublattice to find the correct, otherwise always the first occurence
! In a sigma (Fe)(Cr)(Cr,Fe) y(sigma,cr)=1 but y(sigma,cr#3) gives Cr in third
      indices(3)=icc
   elseif((is.ge.12 .and. is.le.15) .or. is.eq.17 .or. is.ge.19) then
! There must be an argument for NP, BP, DG, Q, X, W, Y, TC and BMAG
      gx%bmperr=4111; goto 1000
   elseif(norm.eq.4) then
! there must be a phase specification for a quantity per formula unit
      gx%bmperr=4115; goto 1000
   endif
!   if(is.eq.17 .or. is.eq.19) then
!      is=is-1
!      svr%norm=1
   if(is.eq.16) svr%norm=1
   if(is.eq.18) svr%norm=2
!   endif
!-----------------------
500 continue
!-----------------------------------------------------------------------
! U       1x   (phase,composition set)             Internal energy (J)
! S       2x                                       entropy
! V       3x                                       volume
! H       4x                                       enthalpy
! A       5x                                       Helmholtz energy
! G       6x                                       Gibbs energy
! NP      7x   phase                               moles of phase
! BP      8x   phase                               mass of phase
! N       9x   (component/phase,component)         moles           >>14
! X       9x   component/phase,component           mole fraction   >>15
! B       10x  (component/phase,component)         mass            >>16
! W       10x                                      mass fraction   >>17
! Y       11    phase,constituent#sublattice       constituent fraction >>18
! Q       12                                       Internal stability   >>19
! DG      13x                                      Driving force
! TC, BM, MQ& etc (model variables)
   svr%statevarid=is
   extensive: if(is.eq.6) then
! U       1x   (phase)                             Internal energy (J)
      istv=10+norm
   elseif(is.eq.7) then
! S       2x                                       entropy
      istv=20+norm
   elseif(is.eq.8) then
! V       3x                                       volume
      istv=30+norm
   elseif(is.eq.9) then
! H       4x                                       enthalpy
      istv=40+norm
   elseif(is.eq.10) then
! A       5x                                       Helmholtz energy
      istv=50+norm
   elseif(is.eq.11) then
! G       6x                                       Gibbs energy
      istv=60+norm
   elseif(is.eq.12) then
! NP      7x    phase                              moles of phase
      istv=70+norm
   elseif(is.eq.13) then
! BP      8x   phase                               mass of phase
      istv=80+norm
   elseif(is.eq.14) then
! DG      9x                                       Driving force
      istv=90+norm
   elseif(is.eq.15) then
! Q      10x                                       Internal stability
      istv=100+norm
   elseif(is.eq.16 .or. is.eq.17) then
! N       11x    (component/phase,component)       moles
! X=NM    111                                      mole fraction
! X%      111, iunit=100                           mole percent
      if(is.eq.16) then
         istv=110+norm
      else
         istv=111
      endif
   elseif(is.eq.18 .or. is.eq.19) then
! B       12x    (component/phase,component)        mass
! W=BW    122                                       mass fraction
! W%      122, iunit=100                            mass percent
      if(is.eq.18) then
         istv=120+norm
      else
         istv=122
      endif
   elseif(is.eq.20) then
! Y       130    phase#comp.set,constituent#sublat  constituent fraction
      istv=130
   else
! the symbol may be a property
      if(deblist) write(*,*)'maybe a property ',is
      goto 600
   endif extensive
   goto 1000
!------------------------------------------------
! handling of properties like TC, BMAGN, MQ etc
600 continue
! the symbol may be a property symbol
   propsym=statevar
! second argument 0 means a symbol
   call find_defined_property(propsym,0,typty,iph,ics)
   if(deblist) write(*,*)'25C at 600: ',propsym(1:len_trim(propsym)),typty
   if(gx%bmperr.ne.0) then
      svr%oldstv=-1; goto 1000
   endif
   indices(1)=iph
   indices(2)=ics
   svr%phase=iph
   svr%compset=ics
!----------------------------- unfinished ?????
   if(typty.gt.100) then
! typty: third argument is constituent (or component??)
      istv=-typty/100
      indices(3)=typty+100*istv
      svr%argtyp=4
   elseif(typty.gt.1) then
      istv=-typty
      svr%argtyp=3
      svr%argtyp=2
   else
! unknown propery
      write(*,*)'Unknown state variable or property',typty
      gx%bmperr=7777; goto 1000
   endif
   svr%oldstv=istv
   svr%statevarid=istv
   svr%constituent=indices(3)
   if(deblist) write(*,611)'Property: ',is,istv,typty,indices
611 format(a,10i4)
!------------------------------------------------
1000 continue
! accept the current istv as svr%oldstv, store a suffix S on MU as phref<0
   svr%oldstv=istv
   svr%phref=iref
   if(deblist) write(*,1001)'25C exit decode: ',istv,(indices(is),is=1,4),&
        norm,iref,iunit,svr%oldstv,svr%phase,svr%compset,svr%component,&
        svr%constituent,svr%norm,svr%phref,svr%unit,svr%argtyp,&
        svr%statevarid,gx%bmperr
1001 format(a,i5,4i3,2x,3i5/17x,i5,4i3,2x,6i5)
   return
!---------------- errors -------------------------------
! Wrong first character of state variable
1100 continue
   gx%bmperr=4099; goto 1000
! M not followed by U
!1110 continue
!   gx%bmperr=4100; goto 1000
! L not followed by NAC
!1120 continue
!   gx%bmperr=4101; goto 1000
! No opening ( for arguments
1130 continue
   gx%bmperr=4102; goto 1000
! No closing ) for arguments
1140 continue
   gx%bmperr=4103; goto 1000
! Unknown phase used as argument in state variable
1150 continue
   gx%bmperr=4104; goto 1000
! No such constituent
1160 continue
   gx%bmperr=4105; goto 1000
! No such component
1170 continue
   gx%bmperr=4106; goto 1000
 end subroutine decode_state_variable3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_phase_molmass(iph,ics,xmol,wmass,totmol,totmass,amount,ceq)
! calculates mole fractions and mass fractions for a phase#set
! xmol and wmass are fractions of components in mol or mass
! totmol is total number of moles and totmass total mass of components.
! amount is number of moles of components per formula unit.
   implicit none
   TYPE(gtp_equilibrium_data) :: ceq
   integer iph,ics
   double precision, dimension(*) :: xmol,wmass
   double precision amount,totmol,totmass
!\end{verbatim}
   integer ic,jc,lokph,lokcs,ll,iel,lokel,ie,kk,loksp
   double precision as,yz,xsum,wsum
   double precision, dimension(maxel) :: x2mol,w2mass
!
   do ic=1,noofel
      xmol(ic)=zero
      wmass(ic)=zero
   enddo
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   ic=0
!
! bug here when calculating Cr-Fe because we create new composition set ...
   if(ocv()) write(*,14)'25c cpm: ',iph,ics,lokph,lokcs
14 format(a,10i5)
   allsubl: do ll=1,phlista(lokph)%noofsubl
      as=ceq%phase_varres(lokcs)%sites(ll)
      allcons: do kk=1,phlista(lokph)%nooffr(ll)
         ic=ic+1
         if(.not.btest(ceq%phase_varres(lokcs)%constat(ic),CONSUS)) then
            yz=ceq%phase_varres(lokcs)%yfr(ic)
            loksp=phlista(lokph)%constitlist(ic)
! isq just for debug output
!            isq=splista(loksp)%alphaindex
!             write(*,11)'cpm 3: ',lokph,lokcs,loksp,splista(loksp)%noofel
!11           format(a,5i3)
            do iel=1,splista(loksp)%noofel
               lokel=splista(loksp)%ellinks(iel)
               ie=ellista(lokel)%alphaindex
               if(ie.ne.0) then
                  xmol(ie)=xmol(ie)+&
                       as*yz*splista(loksp)%stoichiometry(iel)
               endif
            enddo
!            if(ie.gt.0) then
!               write(*,711)ic,loksp,isq,lokel,ie,yz,xmol(ie)
!            else
!               write(*,711)ic,loksp,isq,lokel,ie,yz
!            endif
!711         format('cpmm: ',5i9,2F7.4)
         endif
      enddo allcons
   enddo allsubl
! normallize, All ok here
!   write(*,713)'A',noofel,(xmol(iq),iq=1,noofel)
713 format('25c x:',a,i2,10f7.4)
!800 continue
   xsum=zero
   wsum=zero
! here xmol(i) is equal to the number of moles of element i per formula unit
! set wmass(i) to the mass of of element i per mole formula unit and sum
   do ic=1,noofel
      wmass(ic)=xmol(ic)*ellista(elements(ic))%mass
      xsum=xsum+xmol(ic)
      wsum=wsum+wmass(ic)
   enddo
!   write(*,713)'F',noofel,xsum,(xmol(iq),iq=1,noofel)
   do ic=1,noofel
      xmol(ic)=xmol(ic)/xsum
      wmass(ic)=wmass(ic)/wsum
   enddo
! This is the current number of formula unit of the phase, zero if not stable
   amount=ceq%phase_varres(lokcs)%amfu
! ceq%phase_varres(lokcs)%abnorm(1) is moles atoms for one formula unit
! ceq%phase_varres(lokcs)%abnorm(2) is mass for one formula unit
   totmol=amount*xsum
   totmass=amount*wsum
!   write(*,713)'G',noofel,totmol,totmass,amount
! all seems OK here
!   write(*,811)xsum,ceq%phase_varres(lokcs)%abnorm(1),&
!        wsum,ceq%phase_varres(lokcs)%abnorm(2),amount,totmass
!   write(*,811)xsum,ceq%phase_varres(lokcs)%abnorm(1),&
!        wsum,ceq%phase_varres(lokcs)%abnorm(2),amount,totmass
811 format('cphmm: ',6(1pe12.4))
!   write(*,*)'cpmm: ',totmol,totmass
! all calculation so far in elements, convert to current components
! NOTE: sum of mole fractions can be zero or negative with other 
! components than elements
76 format(a,10F7.4)
78 format(a,2i3,3(1PE12.4))
!   do ic=1,noofel
!      write(*,298)(ceq%invcompstoi(jc,ic),jc=1,noofel)
!   enddo
!298 format('25C: ',6(1pe12.4))
   goto 1000
! what is this ... converting to user defined components ... (not implemented)
   x2mol=zero
   w2mass=zero
   do ic=1,noofel
      do jc=1,noofel
         x2mol(ic)=x2mol(ic)+ceq%invcompstoi(jc,ic)*xmol(jc)
!         write(*,78)'addon: ',ic,jc,x2mol(ic),ceq%invcompstoi(jc,ic),xmol(jc)
         w2mass(ic)=w2mass(ic)+ceq%invcompstoi(ic,jc)*wmass(jc)
      enddo
   enddo
!   do ic=1,noofel
!      write(*,99)'ci: ',(ceq%invcompstoi(jc,ic),jc=1,noofel)
!   enddo
99 format(a,7e11.3)
!   write(*,76)'cmm2: ',(x2mol(ic),ic=1,noofel)
   do ic=1,noofel
      xmol(ic)=x2mol(ic)
      wmass(ic)=w2mass(ic)
   enddo
! something wrong between writing label 713 above and here !!!!!!!!!!!!!
!   write(*,713)'B',noofel,(xmol(iq),iq=1,noofel)
1000 continue
   return
 end subroutine calc_phase_molmass

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_phase_mol(iph,xmol,ceq)
! calculates mole fractions for phase iph, compset 1 in equilibrium ceq
! used for grid generation and some other things
! returns current constitution in xmol equal to mole fractions of components
   implicit none
   integer iph
   double precision xmol(*)
   TYPE(gtp_equilibrium_data) :: ceq
!\end{verbatim}
   integer ic,lokph,lokcs,ll,kk,loksp,lokel,iel,ie
   double precision as,yz,xsum
   do ic=1,noofel
      xmol(ic)=zero
   enddo
   lokph=phases(iph)
   lokcs=phlista(lokph)%linktocs(1)
   ic=0
   allsubl: do ll=1,phlista(lokph)%noofsubl
      as=ceq%phase_varres(lokcs)%sites(ll)
      allcons: do kk=1,phlista(lokph)%nooffr(ll)
         ic=ic+1
         if(.not.btest(ceq%phase_varres(lokcs)%constat(ic),CONSUS)) then
            yz=ceq%phase_varres(lokcs)%yfr(ic)
            loksp=phlista(lokph)%constitlist(ic)
            do iel=1,splista(loksp)%noofel
               lokel=splista(loksp)%ellinks(iel)
               ie=ellista(lokel)%alphaindex
               if(ie.ne.0) then
                  xmol(ie)=xmol(ie)+&
                       as*yz*splista(loksp)%stoichiometry(iel)
               endif
            enddo
         endif
      enddo allcons
   enddo allsubl
! normallize
   xsum=zero
   do ic=1,noofel
      xsum=xsum+xmol(ic)
   enddo
   do ic=1,noofel
      xmol(ic)=xmol(ic)/xsum
   enddo
1000 continue
   return
 end subroutine calc_phase_mol

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_molmass(xmol,wmass,totmol,totmass,ceq)
! summing up N and B for each component over all phases with positive amount
! Check that totmol and totmass are correct ....
   implicit none
   double precision, dimension(*) :: xmol,wmass
   double precision totmol,totmass
   TYPE(gtp_equilibrium_data) :: ceq
!\end{verbatim}
   double precision am,amult,tmol,tmass
   double precision, dimension(maxel) :: xph,wph
   integer ic,iph,lokph,ics,lokcs
   do ic=1,noofel
      xmol(ic)=zero
      wmass(ic)=zero
   enddo
   totmol=zero
   totmass=zero
   allph: do iph=1,noofph
      lokph=phases(iph)
      if(.not.btest(phlista(lokph)%status1,phhid)) then
         allcs: do ics=1,phlista(lokph)%noofcs
            lokcs=phlista(lokph)%linktocs(ics)
! ceq%phase_varres(lokcs)%amfu is current number of formula units
! ceq%phase_varres(lokcs)%abnorm(1) is number of real atoms in a formula unit
            am=ceq%phase_varres(lokcs)%amfu*&
                 ceq%phase_varres(lokcs)%abnorm(1)
            if(am.gt.zero) then
               call calc_phase_molmass(iph,ics,xph,wph,tmol,tmass,amult,ceq)
               if(gx%bmperr.ne.0) goto 1000
!               write(*,17)'25c amult:',iph,ics,am,amult,tmol,tmass
!               write(*,18)'25c x0: ',(xph(ic),ic=1,noofel)
!               write(*,18)'25c w0: ',(wph(ic),ic=1,noofel)
17             format(a,2i4,6(1pe14.6))
18             format(a,8(F9.5))
               do ic=1,noofel
                  xmol(ic)=xmol(ic)+am*xph(ic)
                  wmass(ic)=wmass(ic)+tmass*wph(ic)
               enddo
               totmass=totmass+tmass
               totmol=totmol+tmol
            endif
         enddo allcs
      endif
   enddo allph
! we have summed the number of moles and mass of all elements in all phases
!   xsum=zero
!   wsum=zero
!   do ic=1,noofel
!      xsum=xsum+xmol(ic)
!      wsum=wsum+wmass(ic)
!   enddo
!   write(*,21)'25C x1: ',xsum,totmol,(xmol(ic),ic=1,noofel)
!   write(*,21)'25C w2: ',wsum,totmass,(wmass(ic),ic=1,noofel)
!21 format(a,2(1pe12.4),10(0pF9.4))
   if(totmass.gt.zero) then
      do ic=1,noofel
         xmol(ic)=xmol(ic)/totmol
         wmass(ic)=wmass(ic)/totmass
      enddo
!   else
!      write(*,*)'There is no mass at all in the system!'
!      gx%bmperr=4185; goto 1000
   endif
!   write(*,21)'25C x1: ',totmol,(xmol(ic),ic=1,noofel)
!   write(*,21)'25C w1: ',totmass,(wmass(ic),ic=1,noofel)
21 format(a,1pe12.4,8(0pF9.5))
!   else
! this is not an error if no calculation has been made
!      write(*,28)'25C: calc_molmass: No mole fractions',totmol,totmass,xsum,&
!           (xmol(ic),ic=1,noofel)
28    format(a,3(1pe12.4)/'25C. ',10f7.4)
!      gx%bmperr=4185; goto 1000
!   endif
!   wsum=zero
!   do ic=1,noofel
!      wmass(ic)=xmol(ic)*ellista(elements(ic))%mass
!      wsum=wsum+wmass(ic)
!      write(*,44)'cmm4: ',ic,xmol(ic),wmass(ic),&
!           ellista(elements(ic))%mass,wsum,totmass
44    format(a,i3,6(1pe12.4))
!   enddo
!   if(wsum.gt.zero) then
!      do ic=1,noofel
!         wmass(ic)=wmass(ic)/wsum
!      enddo
!   endif
1000 continue
   return
 end subroutine calc_molmass

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine sumprops(props,ceq)
! summing up G, S, V, N and B for all phases with positive amount
! Check if this is correct
   implicit none
   TYPE(gtp_equilibrium_data) :: ceq
   double precision props(5)
!\end{verbatim}
   integer lokph,lokcs,ics
   double precision am
   if(gx%bmperr.ne.0) write(*,*)'error entering sumprops ',gx%bmperr
   props=zero
   allph: do lokph=1,noofph
!      write(*,*)'25c sumprops: ',lokph
      if(.not.btest(phlista(lokph)%status1,phhid)) then
!         lokcs=phlista(lokph)%cslink
         allcs: do ics=1,phlista(lokph)%noofcs
! phase_varres(lokcs)%amfu is the amount formula units of the phase
! phase_varres(lokcs)%abnorm(1) is the moles of real atoms/formula unit
! am is the number of moles of real atoms of the phase
            lokcs=phlista(lokph)%linktocs(ics)
            am=ceq%phase_varres(lokcs)%amfu*&
                 ceq%phase_varres(lokcs)%abnorm(1)
!            write(*,*)'25c sumprops: ',lokph,am
            if(am.gt.zero) then
! properties are G, G.T=-S, G.P=V and moles and mass of real atoms
! Note gval(*,1) is per mole formula unit and ceq%phase_varres(lokcs)%abnorm(1)
! is the number of real atoms per formula unit
               props(1)=props(1)+am*ceq%phase_varres(lokcs)%gval(1,1)/&
                    ceq%phase_varres(lokcs)%abnorm(1)
               props(2)=props(2)+am*ceq%phase_varres(lokcs)%gval(2,1)/&
                    ceq%phase_varres(lokcs)%abnorm(1)
               props(3)=props(3)+am*ceq%phase_varres(lokcs)%gval(3,1)/&
                    ceq%phase_varres(lokcs)%abnorm(1)
               props(4)=props(4)+am
! ceq%phase_varres(lokcs)%abnorm(2) should be the current mass
! but I am not sure it is updated for the current composition
!               props(5)=props(5)+am*ceq%phase_varres(lokcs)%abnorm(2)/&
!                    ceq%phase_varres(lokcs)%abnorm(1)
! I think abnorm(2) is actual mass
!               props(5)=props(5)+ceq%phase_varres(lokcs)%abnorm(2)
               props(5)=props(5)+am*ceq%phase_varres(lokcs)%abnorm(2)
!               write(*,11)'25C sumprops: ',lokcs,props(1),props(4),props(5),&
!                    ceq%phase_varres(lokcs)%abnorm(2)
!               write(*,11)'sumprops ',lokcs,am,props(4),&
!                    ceq%phase_varres(lokcs)%abnorm(1)
11             format(a,i4,6(1pe12.4))
            endif
         enddo allcs
      endif
   enddo allph
1000 continue
   if(gx%bmperr.ne.0) write(*,*)'error exiting sumprops ',gx%bmperr
   return
 end subroutine sumprops

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine encode_state_variable(text,ip,svr,ceq)
! writes a state variable in text form position ip.  ip is updated
   character text*(*)
   integer ip
   type(gtp_state_variable), pointer :: svr
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer istv,indices(4),iunit,iref
   iref=svr%phref
   iunit=svr%unit
! if svr%oldstv>=10 then istv should be 10*(svr%oldstv-5)+svr%norm
!   if(svr%oldstv.ge.10) then
!      istv=10*(svr%oldstv-5)+svr%norm
!      write(*,*)'25C encode: ',svr%oldstv,svr%norm,istv
!   else
      istv=svr%oldstv
!   endif
! svr%argtyp specifies values in indices:
! svr%argtyp: 0=no arguments; 1=comp; 2=ph+cs; 3=ph+cs+comp; 4=ph+cs+const
   indices=0
   if(svr%argtyp.eq.1) then
      indices(1)=svr%component
   elseif(svr%argtyp.eq.2) then
      indices(1)=svr%phase
      indices(2)=svr%compset
   elseif(svr%argtyp.eq.3) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%component
   elseif(svr%argtyp.eq.4) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%constituent
!   else
!      write(*,*)'state variable has illegal argtyp: ',svr%argtyp
!      gx%bmperr=7775; goto 1000
   endif
   call encode_state_variable3(text,ip,istv,indices,iunit,iref,ceq)
1000 continue
   return
 end subroutine encode_state_variable

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine encode_state_variable3(text,ip,istv,indices,iunit,iref,ceq)
! writes a state variable in text form position ip.  ip is updated
! the internal coding provides in istv, indices, iunit and iref
! ceq is needed as compopnents can be different in different equilibria ??
! >>>> unfinished as iunit and iref not really cared for
   implicit none
   integer, parameter :: noos=20
   character*4, dimension(noos), parameter :: svid = &
       ['T   ','P   ','MU  ','AC  ','LNAC','U   ','S   ','V   ',&
        'H   ','A   ','G   ','NP  ','BP  ','DG  ','Q   ','N   ',&
        'X   ','B   ','W   ','Y   ']
   character*(*) text
   integer, dimension(4) :: indices
   integer istv,ip,iunit,iref
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer jp,ics,kstv,iph,norm
!
   character stsymb*60
   character*1, dimension(4), parameter :: cnorm=['M','W','V','F']
!
   if(istv.le.0) then
! this is a parameter property symbol: TC (-2), BM (-3), MQ&FE(FCC) (-4) etc
! translate to 21, 22, 23 ...
      kstv=19-istv
      goto 200
!      gx%bmperr=4116; goto 1000
   endif
! T or P
   if(istv.le.2) then
      text(ip:ip)=svid(istv)
      ip=ip+1
      goto 1000
   endif
   stsymb=' '
!   potential: if(istv.le.6) then
   potential: if(istv.le.5) then
! Potential, MU, AC or LNAC, possible suffix 'S' for SER
      stsymb=svid(istv)
      jp=len_trim(stsymb)+1
!      if(iref.ne.0) then
      if(iref.lt.0) then
! New use of svr%phref and iref, <0 means use SER as reference state
         stsymb(jp:jp)='S'
         jp=jp+1
      endif
      stsymb(jp:jp)='('
      jp=jp+1
      if(indices(2).eq.0) then
! problem ... component names different in different equilibria ....
         call get_component_name(indices(1),stsymb(jp:),ceq)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
      else
! always use composition set 1
         ics=1
         call get_phase_name(indices(1),ics,stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=','; jp=jp+1
         call get_phase_constituent_name(indices(1),indices(2),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
      endif
      stsymb(jp:jp)=')'
      goto 800
   endif potential
   if(istv.lt.10) then
!      write(*,*)'unknown potential'
      gx%bmperr=4158; goto 1000
   endif
! Extensive property has istv>=10
   norm=mod(istv,10)
   kstv=(istv+1)/10+5
!    write(*,3)'encode 3: ',kstv,indices
!3   format(a,5i5)
   if(kstv.eq.16 .and. norm.eq.1) then
! NM should be X
      if(indices(1).ne.0) kstv=17
   elseif(kstv.eq.17) then
! BW should be W
      if(norm.eq.2 .and. indices(1).ne.0) then
         kstv=19
      else
         kstv=18
      endif
   elseif(kstv.ge.18) then
! Y
!      kstv=kstv+2
      kstv=20
   endif
!   write(*,11)'esv 7: ',istv,kstv,indices
11 format(a,10i4)
   stsymb=svid(kstv)
   jp=len_trim(stsymb)+1
!   write(*,*)'25C norm 1A: ',kstv,norm
   if(kstv.le.16 .or. kstv.eq.18) then
      if(norm.gt.0 .and. norm.le.4) then
!         write(*,*)'25C norm 1B: ',kstv,norm
         stsymb(jp:jp)=cnorm(norm)
         jp=jp+1
      elseif(norm.ne.0) then
!         write(*,*)'25C norm 1C: ',kstv,norm
         gx%bmperr=4118; goto 1000
      endif
   endif
   goto 500
!----------------- 
! parameter property symbols
200 continue
   iph=indices(1)
   ics=indices(2)
   if(indices(3).ne.0) then
      kstv=-100*istv+indices(3)
   else
      kstv=-istv
   endif
! this call creates the symbol or gives an error
   call find_defined_property(stsymb,1,kstv,iph,ics)
   if(gx%bmperr.ne.0) goto 1000
   jp=len_trim(stsymb)+1
   goto 800
!------------------
! handle indices
500 continue
   noind: if(indices(3).gt.0) then
! 3 indices, phase, comp.set and constituent allowed for Y
! or phase, comp.set and component, allowed for N, X, B and W
! or phase, comp.set and constituent allowed for MQ&
      if(kstv.eq.20) then
! this is Y
         stsymb(jp:jp)='('
         jp=jp+1
         call get_phase_name(indices(1),indices(2),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=','
         jp=jp+1
         call get_phase_constituent_name(indices(1),indices(3),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=')'
         jp=jp+1
      elseif(kstv.ge.16 .and. kstv.le.19) then
! allow for percent or %
         if(iunit.eq.100) then
            stsymb(jp:jp+1)='%('
            jp=jp+2
         else
            stsymb(jp:jp)='('
            jp=jp+1
         endif
         call get_phase_name(indices(1),indices(2),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=','
         jp=jp+1
         call get_component_name(indices(3),stsymb(jp:),ceq)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=')'
         jp=jp+1
      else
         gx%bmperr=4117; goto 1000
      endif
   elseif(indices(2).gt.0) then
! 2 indices, can only be phase and comp.set
         stsymb(jp:jp)='('
         jp=jp+1
         call get_phase_name(indices(1),indices(2),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=')'
         jp=jp+1
   elseif(indices(1).gt.0) then
! 1 index, can only be component
! allow for percent or %
         if(iunit.eq.100) then
            stsymb(jp:jp+1)='%('
            jp=jp+2
         else
            stsymb(jp:jp)='('
            jp=jp+1
         endif
         call get_component_name(indices(1),stsymb(jp:),ceq)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=')'
         jp=jp+1
! >>>> unfinished
   endif noind
!
800 continue
   text(ip:ip+jp-1)=stsymb
   ip=ip+jp
   if(text(ip:ip).eq.' ') then
! remove a trailing space occuring in some cases
      ip=ip-1
   endif
1000 continue
   return
 end subroutine encode_state_variable3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine encode_state_variable_record(text,ip,svr,ceq)
! writes a state variable in text form position ip.  ip is updated
! the svr record provide istv, indices, iunit and iref
! ceq is needed as compopnents can be different in different equilibria ??
! >>>> unfinished as iunit and iref not really cared for
   implicit none
   integer, parameter :: noos=20
   character*4, dimension(noos), parameter :: svid = &
       ['T   ','P   ','MU  ','AC  ','LNAC','U   ','S   ','V   ',&
        'H   ','A   ','G   ','NP  ','BP  ','DG  ','Q   ','N   ',&
        'X   ','B   ','W   ','Y   ']
   character*(*) text
   type(gtp_state_variable) :: svr
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer jp,ics,kstv,iph,norm
   integer, dimension(4) :: indices
   integer istv,ip,iunit,iref
!
   character stsymb*60
   character*1, dimension(4), parameter :: cnorm=['M','W','V','F']
!
   istv=svr%oldstv
   norm=svr%norm
   iunit=svr%unit
   indices=0
   if(svr%argtyp.eq.1) then
      indices(1)=svr%component
   elseif(svr%argtyp.eq.2) then
      indices(1)=svr%phase
      indices(2)=svr%compset
   elseif(svr%argtyp.eq.3) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%component
   elseif(svr%argtyp.eq.4) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%constituent
   endif
! there is some cloudy thinking here.  If the user has defined his own
! reference state that should be used.  The information is stored in the
! component record (ceq%complist(i)%phlink
! But if the user specifies MUS(i) one should use SER ... how to transfer that
! information to the calculating routines?
! By default svr%phref=0, then use user defined.  If phref<0 use SER ??
   iref=svr%phref
!
   if(istv.le.0) then
! this is a parameter property symbol: TC (-2), BM (-3), MQ&FE(FCC) (-4) etc
! translate to 21, 22, 23 ...
      kstv=19-istv
      goto 200
!      gx%bmperr=4116; goto 1000
   endif
! T or P
   if(istv.le.2) then
      text(ip:ip)=svid(istv)
      ip=ip+1
      goto 1000
   endif
   stsymb=' '
!   potential: if(istv.le.6) then
   potential: if(istv.le.5) then
! Potential, MU, AC or LNAC, possible suffix 'S' for SER
      stsymb=svid(istv)
      jp=len_trim(stsymb)+1
!      if(iref.ne.0) then
      if(iref.lt.0) then
! new use of phref and iref, <0 means use SER and suffix S
         stsymb(jp:jp)='S'
         jp=jp+1
      endif
      stsymb(jp:jp)='('
      jp=jp+1
      if(indices(2).eq.0) then
! problem ... component names can be different in different equilibria ....
         call get_component_name(indices(1),stsymb(jp:),ceq)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
      else
! always use composition set 1
         ics=1
         call get_phase_name(indices(1),ics,stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=','; jp=jp+1
         call get_phase_constituent_name(indices(1),indices(2),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
      endif
      stsymb(jp:jp)=')'
      goto 800
   endif potential
   if(istv.lt.10) then
!      write(*,*)'unknown potential'
      gx%bmperr=4158; goto 1000
   endif
! Extensive property has istv>=10
   norm=mod(istv,10)
   kstv=(istv+1)/10+5
!    write(*,3)'encode 3: ',kstv,indices
!3   format(a,5i5)
   if(kstv.eq.16 .and. norm.eq.1) then
! NM should be X
      if(indices(1).ne.0) kstv=17
   elseif(kstv.eq.17) then
! BW should be W
      if(norm.eq.2 .and. indices(1).ne.0) then
         kstv=19
      else
         kstv=18
      endif
   elseif(kstv.ge.18) then
! Y
!      kstv=kstv+2
      kstv=20
   endif
!    write(*,11)'esv 7: ',istv,kstv,indices
11  format(a,10i4)
   stsymb=svid(kstv)
   jp=len_trim(stsymb)+1
   write(*,*)'25C norm 2: ',kstv,norm
   if(kstv.le.16 .or. kstv.eq.18) then
      if(norm.gt.0 .and. norm.le.4) then
         stsymb(jp:jp)=cnorm(norm)
         jp=jp+1
      elseif(norm.ne.0) then
         gx%bmperr=4118; goto 1000
      endif
   endif
   goto 500
!----------------- 
! parameter property symbols
200 continue
   iph=indices(1)
   ics=indices(2)
   if(indices(3).ne.0) then
      kstv=-100*istv+indices(3)
   else
      kstv=-istv
   endif
! this call creates the symbol or gives an error
   call find_defined_property(stsymb,1,kstv,iph,ics)
   if(gx%bmperr.ne.0) goto 1000
   jp=len_trim(stsymb)+1
   goto 800
!------------------
! handle indices
500 continue
   noind: if(indices(3).gt.0) then
! 3 indices, phase, comp.set and constituent allowed for Y
! or phase, comp.set and component, allowed for N, X, B and W
! or phase, comp.set and constituent allowed for MQ&
      if(kstv.eq.20) then
! this is Y
         stsymb(jp:jp)='('
         jp=jp+1
         call get_phase_name(indices(1),indices(2),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=','
         jp=jp+1
         call get_phase_constituent_name(indices(1),indices(3),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=')'
         jp=jp+1
      elseif(kstv.ge.16 .and. kstv.le.19) then
! allow for percent or %
         if(iunit.eq.100) then
            stsymb(jp:jp+1)='%('
            jp=jp+2
         else
            stsymb(jp:jp)='('
            jp=jp+1
         endif
         call get_phase_name(indices(1),indices(2),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=','
         jp=jp+1
         call get_component_name(indices(3),stsymb(jp:),ceq)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=')'
         jp=jp+1
      else
         gx%bmperr=4117; goto 1000
      endif
   elseif(indices(2).gt.0) then
! 2 indices, can only be phase and comp.set
         stsymb(jp:jp)='('
         jp=jp+1
         call get_phase_name(indices(1),indices(2),stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=')'
         jp=jp+1
   elseif(indices(1).gt.0) then
! 1 index, can only be component
! allow for percent or %
         if(iunit.eq.100) then
            stsymb(jp:jp+1)='%('
            jp=jp+2
         else
            stsymb(jp:jp)='('
            jp=jp+1
         endif
         call get_component_name(indices(1),stsymb(jp:),ceq)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=')'
         jp=jp+1
! >>>> unfinished
   endif noind
!
800 continue
   text(ip:ip+jp-1)=stsymb
   ip=ip+jp
   if(text(ip:ip).eq.' ') then
! remove a trailing space occuring in some cases
      ip=ip-1
   endif
1000 continue
   return
 end subroutine encode_state_variable_record

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine state_variable_val(svr,value,ceq)
! calculate the value of a state variable in equilibrium record ceq
! It transforms svr data to old format and calls state_variable_val3
   type(gtp_state_variable), pointer :: svr
   double precision value
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer istv, indices(4),iref,iunit
!
   iref=svr%phref
   iunit=svr%unit
!   if(svr%oldstv.gt.10) then
!      istv=10*(svr%oldstv-5)+svr%norm
!   else
      istv=svr%oldstv
!   endif
! svr%argtyp specifies values in indices:
! svr%argtyp: 0=no arguments; 1=comp; 2=ph+cs; 3=ph+cs+comp; 4=ph+cs+const
   indices=0
   if(svr%argtyp.eq.1) then
      indices(1)=svr%component
   elseif(svr%argtyp.eq.2) then
      indices(1)=svr%phase
      indices(2)=svr%compset
   elseif(svr%argtyp.eq.3) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%component
   elseif(svr%argtyp.eq.4) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%constituent
   elseif(svr%argtyp.ne.0) then
      write(*,*)'25C state variable has illegal argtyp: ',svr%argtyp
      gx%bmperr=7775; goto 1000
   endif
!   write(*,910)'25C svv: ',istv,indices,iref,iunit,value
910 format(a,i3,2x,4i3,2i3,1pe14.6)
   call state_variable_val3(istv,indices,iref,iunit,value,ceq)
   if(gx%bmperr.ne.0) then
      write(*,920)'25C error: ',gx%bmperr,istv,svr%oldstv,svr%argtyp
920   format(a,i5,2x,2i4,i2)
!   else
!      write(*,*)'25C value: ',value
   endif
1000 continue
   return
 end subroutine state_variable_val

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine state_variable_val3(istv,indices,iref,iunit,value,ceq)
! calculate the value of a state variable in equilibrium record ceq
! istv is state variable type (integer)
! indices are possible specifiers
! iref indicates use of possible reference state, 0 current, -1 SER
! iunit is unit, (K, oC, J, cal etc). For % it is 100
! value is the calculated values. for state variables with wildcards use
! get_many_svar
   implicit none
   integer, dimension(4) :: indices
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer istv,iref,iunit
   double precision value
!\end{verbatim}
   double precision props(5),xmol(maxel),wmass(maxel),stoi(10),cmpstoi(10)
   double precision vt,vp,amult,vg,vs,vv,div,aref,vn,bmult,tmass,tmol
   double precision qsp,gref,spmass,rmult
   integer kstv,norm,lokph,lokcs,icx,jp,ncmp,ic,iprop,loksp,nspel
   integer endmember(maxsubl),ielno(maxspel)
   value=zero
   ceq%rtn=globaldata%rgas*ceq%tpval(1)
!   write(*,10)'25C svval: ',istv,indices,iref,iunit,gx%bmperr,value
10 format(a,i4,4i4,3i5,1PE17.6)
   potentials: if(istv.lt.0) then
! negative istv indicate parameter property symbols
      kstv=-istv
      goto 200
!      gx%bmperr=4097; goto 1000
   elseif(istv.ge.10) then
      goto 50
   elseif(istv.eq.1) then
      value=ceq%tpval(1)
   elseif(istv.eq.2) then
      value=ceq%tpval(2)
   else
      if(istv.eq.3) then
! MUx(component) or MU(phase,constituent), x can be S for SER
         goto 500
      elseif(istv.eq.4) then
! ACx(component) or AC(phase,constituent)
         goto 500
      elseif(istv.eq.5) then
! LNACx(component) or LNAC(phase,constituent)
         goto 500
      endif
! wrong or state variable not implemented
      write(*,10)'25C not impl: ',istv,indices,iref,iunit,gx%bmperr,value
      goto 1100
   endif potentials
! normal return
   goto 1000
!----------------------------------------------------------
! extensive variable (N, X, G ...) or model variable (TC, BMAG)
50  continue
   norm=mod(istv,10)
   kstv=istv/10
! this may not be necessary in all cases but do it anyway:
! sum over all stable phases, props(1..3) are G, G.T and G.P,
! props(4) is amount of moles of components, props(5) is mass of components
   call sumprops(props,ceq)
   if(gx%bmperr.ne.0) goto 1000
! if verbose on
   if(ocv()) write(*,51)'25C stv A: ',props
51 format(a,5(1PE12.3))
! kstv can be 1 to 15 for different properties
! norm can be 1, 2, 3 or 4 for normalizing. 0 for not normallizing
!             M  W  V    F
! OLD: iref can be 0 or 1 for reference state
! iref can be 0 for using current referennce state
! iref <0 for default reference state (SER)
   le10: if(kstv.le.10) then
! kstv=       1  2  3  4  5  6  7   8   9     10
! state var;  U, S, V, H, A, G, NP, BP, DG and Q
      vt=ceq%tpval(1)
      vp=ceq%tpval(2)
!      ceq%rtn=globaldata%rgas*ceq%tpval(1)
      amult=ceq%rtn
!      write(*,*)'stv B: ',vt,vp,amult
      if(indices(1).eq.0) then
         vg=props(1)
         vs=-props(2)
         vv=props(3)
! normalizing
         if(norm.eq.1) then
            div=props(4)
         elseif(norm.eq.2) then
            div=props(5)
         elseif(norm.eq.3) then
            div=props(3)
            if(div.eq.zero) then
               gx%bmperr=4114; goto 1000
            endif
         elseif(norm.eq.4) then
            gx%bmperr=4115; goto 1000
         else
            div=one
         endif
! for phase specific the aref should be independent of amult and div ??
! for system wide these are unity
         rmult=one
      else
! phase specific, indices are phase and composition set
         call get_phase_compset(indices(1),indices(2),lokph,lokcs)
         if(gx%bmperr.ne.0) goto 1000
         vg=ceq%phase_varres(lokcs)%gval(1,1)
         vs=-ceq%phase_varres(lokcs)%gval(2,1)
         vv=ceq%phase_varres(lokcs)%gval(3,1)
         if(norm.eq.1) then
            div=ceq%phase_varres(lokcs)%abnorm(1)
            rmult=div
         elseif(norm.eq.2) then
! abnorm(2) should be the mass per formulat unit
            div=ceq%phase_varres(lokcs)%abnorm(2)
            rmult=div
         elseif(norm.eq.3) then
            div=ceq%phase_varres(lokcs)%gval(3,1)
            if(div.eq.zero) then
               gx%bmperr=4114; goto 1000
            endif
            rmult=div
         elseif(norm.eq.4) then
! per formula unit
            div=one
            rmult=div
         else
! no normalizing for a specific phase, value for current amount
! NOTE amult is alreadt RT
            amult=amult*ceq%phase_varres(lokcs)%amfu
            rmult=ceq%phase_varres(lokcs)%amfu
            div=one
!             div=ceq%phase_varres(lokcs)%abnorm(1)
         endif
! for phase specific the aref is for one mole of atoms and should 
! be independent of amult and div ??
!         if(amult.eq.zero) then
!            rmult=zero
!         else
!            rmult=div/amult
!         endif
      endif
! here the reference state should be considered
!      aref=zero
      if(iref.eq.0) then
! >>>> unfinished
         call calculate_reference_state(kstv,indices(1),indices(2),aref,ceq)
         if(gx%bmperr.ne.0) goto 1000
!         write(*,53)'25 C Reference state:',iref,aref,rmult
      elseif(iref.lt.0) then
         aref=zero
      else
         write(*,*)'25C Reference state undefined',iref
         aref=zero
      endif
! if phase specific the scaling for phase specific must be compensated
      aref=rmult*aref
!      write(*,53)'at kstv1: ',kstv,props,aref,div
53    format(a,i3,5(1PE12.3))
      kstv1: if(kstv.eq.1) then
! 1: U = G + TS - PV = G - T*G.T - P*G.P
         value=amult*(vg+vt*vs-vp*vv-aref)/div
      elseif(kstv.eq.2) then
! 2: S = -G.T
         value=amult*(vs-aref)/div
!         write(*,54)value,amult,vs,aref,div
54       format('25C svv: ',5(1pe12.4))
      elseif(kstv.eq.3) then
! 3: V = G.P
         value=amult*(vv-aref)/div
      elseif(kstv.eq.4) then
! 4: H = G + TS = G - T*G.T
         if(ocv()) write(*,177)'25C H:',vg+vt*vs,aref,amult,div,rmult
177      format(a,6(1pe12.4))
         value=amult*(vg+vt*vs-aref)/div
      elseif(kstv.eq.5) then
! 5: A = G - PV = G - P*G.P
         value=amult*(vg-vp*vv-aref)/div
      elseif(kstv.eq.6) then
! 6: G
!         write(*,177)'25C G:',vg,aref,amult,div
         value=amult*(vg-aref)/div
      elseif(kstv.eq.7) then
! 7: NP
         value=ceq%phase_varres(lokcs)%abnorm(1)* &
              ceq%phase_varres(lokcs)%amfu/div
      elseif(kstv.eq.8) then
! 8: BP
! abnorm(2) should be the mass per formula unit
         value=ceq%phase_varres(lokcs)%abnorm(2)* &
              ceq%phase_varres(lokcs)%amfu/div
      elseif(kstv.eq.9) then
! 9: DG (driving force)
!         write(*,202)'svval DG:',lokcs,ceq%phase_varres(lokcs)%dgm,div
202      format(a,i5,2(1pe12.4))
         value=ceq%phase_varres(lokcs)%dgm/div
      elseif(kstv.eq.10) then
! 10: Q (stability, thermodynamic factor), not implemented
         gx%bmperr=4081; goto 1000
!      else
!         write(*,*)'svval after 10:',kstv
      endif kstv1
      goto 1000
   endif le10
!----------------------------------------------------------------------
! here with kstv>10
! kstv=       11  12  13 
! state var:   N   B   Y   
   le12: if(kstv.le.12) then
! normallizing for N (kstv=11) and B (kstv=12)
!      write(*,88)'25c svv 12: ',indices(1),norm,props(5)
88    format(a,2i3,6(1pe12.4))
      if(indices(1).eq.0) then
! no first index means the sum over all phases
! props(4) is amount of moles of components, props(5) is mass of components
         if(kstv.eq.11) then
            vn=props(4)
         else
            vn=props(5)
         endif
! normalizing
         if(norm.eq.1) then
            div=props(4)
         elseif(norm.eq.2) then
            div=props(5)
         elseif(norm.eq.3) then
! we may not have any volume data ...
            div=props(3)
            if(div.eq.zero) then
               gx%bmperr=4114; goto 1000
            endif
         elseif(norm.eq.4) then
            gx%bmperr=4115; goto 1000
         else
            div=one
         endif
! This is N or B without index but possibly normallized
!         write(*,89)'25C svv, N or B: ',vn,div
89       format(a,5(1pe12.4))
         value=vn/div
      else
! one or two indices, overall of phase specific component amount
         if(indices(2).eq.0) then
! one index is component specific, N(comp.), B(CR) etc. Sum over all phases
! props(4) is amount of moles of components, props(5) is mass of components
            call calc_molmass(xmol,wmass,tmol,tmass,ceq)
            if(gx%bmperr.ne.0) goto 1000
!            write(*,89)'25c mm: ',tmol,tmass
!            write(*,93)'25c x: ',(xmol(icx),icx=1,noofel)
!            write(*,93)'25c w: ',(wmass(icx),icx=1,noofel)
93          format(a,9F7.4)
            icx=1
            if(kstv.eq.11) then
               bmult=props(4)
            else
               bmult=props(5)
            endif
         else
! two indices is phase and component specific. bmult is amount of phase
            call calc_phase_molmass(indices(1),indices(2),&
                 xmol,wmass,tmol,tmass,bmult,ceq)
            icx=3
         endif
         if(gx%bmperr.ne.0) goto 1000
!         write(*,13)'gsvv 19: ',norm,(xmol(iq),iq=1,noofel)
777      format('gsvv 77: ',10(f7.4))
         if(kstv.eq.11) then
! total moles of component
            vn=xmol(indices(icx))
            amult=tmol
!            write(*,777)kstv,icx,indices(icx),norm,vn,amult,bmult
!777         format('N(i): ',4i4,3(1pe12.4))
         else
! total mass of component
            vn=wmass(indices(icx))
            amult=tmass
         endif
!         write(*,13)'gsvv 8: ',norm,vn,amult,bmult,tmol,tmass
13       format(a,i3,7(1PE10.2))
         norm3: if(norm.eq.1) then
! NM or X
            if(tmol.ne.zero) then
               value=amult*vn/tmol
            else
! problem at x(phase,component) was zero when phase fix with zero amount
!               value=zero
               value=vn
            endif
! percent %
!            write(*,*)'x%: ',iunit,value
            if(iunit.eq.100) value=1.0D2*value
         elseif(norm.eq.2) then
! NW or W
            if(tmass.gt.zero) then
               value=amult*vn/tmass
            else
               value=zero
            endif
! percent %
            if(iunit.eq.100) value=1.0D2*value
         elseif(norm.eq.3) then
! NV
            if(props(3).gt.zero) then
               value=amult*vn/props(3)
            else
               gx%bmperr=4114
            endif
         elseif(norm.eq.4) then
! NF or BF with one or two indices
            if(indices(2).eq.0) then
               gx%bmperr=4115; goto 1000
            else
               value=vn
            endif
         else
! N(comp), N(phase,comp), B(comp) or B(phase,comp)
            value=bmult*vn
         endif norm3
      endif
      goto 1000
   endif le12
!-----------------------------------------------------------------
! special for Y
   if(kstv.eq.13) then
! 13: Y
      call get_phase_compset(indices(1),indices(2),lokph,lokcs)
      if(gx%bmperr.ne.0) goto 1000
      value=ceq%phase_varres(lokcs)%yfr(indices(3))
   else
! wrong state variable specification
      value=zero
      gx%bmperr=4113
   endif
   goto 1000
!-----------------------------------------------------------------
! values of parameter property symbols
! >>> this can easily be generallized ... next time around ...
! here with state variable <0, syetm and user defined properties
200   continue
   select case(kstv)
   case default
      write(kou,*)'Unknown parameter identifier: ',kstv
!.......................................
   case(2:5,7,9:19) 
! 2: TC (Curie/Neel Temperature)
! 3: BM (Average Bohr magneton number)
! 4: CTA just Curie Temperature
! 5: NTA just Neel temperature
! 7: THET Debye or Einstein temperature
! 9: RHO electrical resistivity
! 10: MAGS Magnetic suseptibility
! 11: GTT Glas transition temperature
! 12: VISC viscosity
! 13: LPX Lattice parameter in X direction
! 14: LPY Lattice parameter in Y direction
! 15: LPZ Lattice parameter in Z direction
! 16: LPTH Lattice angle
! 17: EC11 Elastic constant C11
! 18: EC12 Elastic constant C12
! 19: EC44 Elastic constant C44
      call get_phase_compset(indices(1),indices(2),lokph,lokcs)
      if(gx%bmperr.ne.0) goto 1000
! nprop is number of properties calculated.  Property 1 is always G
      find1: do jp=2,ceq%phase_varres(lokcs)%nprop
! the listprop array contain identification of the property stored there
         if(ceq%phase_varres(lokcs)%listprop(jp).eq.kstv) then
            value=ceq%phase_varres(lokcs)%gval(1,jp)
            goto 1000
         endif
      enddo find1
!.......................................
   case(6,8) 
! 6: IBM& Individual Bohr magneton number
! 8: MQ& mobility value
      call get_phase_compset(indices(1),indices(2),lokph,lokcs)
      if(gx%bmperr.ne.0) goto 1000
! property is kstv*100+indices(3) (constituent identifier)
      iprop=100*kstv+indices(3)
      find2: do jp=2,ceq%phase_varres(lokcs)%nprop
         if(ceq%phase_varres(lokcs)%listprop(jp).eq.iprop) then
            value=ceq%phase_varres(lokcs)%gval(1,jp)
            goto 1000
         endif
      enddo find2
   end select
!.......................................
   gx%bmperr=4113; goto 1000
!-----------------------------------------------------------------
! chemical potentials, activites etc, istv is 3, 4 or 5 for MU, AC and LNAC
! there can be a reference state
500 continue
!   ceq%rtn=globaldata%rgas*ceq%tpval(1)
! if one argument that is a component, if two these are phase and constituent
   if(indices(2).ne.0) then
      lokph=phases(indices(1))
      loksp=phlista(lokph)%constitlist(indices(2))
! split the species in elements, convert to components, add chemical potentials
      call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp)
      if(gx%bmperr.ne.0) goto 1000
      if(qsp.ne.zero) then
!         write(*,*)'Cannot calculate potential of charged species'
         gx%bmperr=4159; goto 1000
      endif
!      write(*,*)'converting to component',nspel,ielno(1),stoi(1)
      call elements2components(nspel,stoi,ncmp,cmpstoi,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,*)'converting to component',ncmp,cmpstoi(1)
      value=zero
      do ic=1,ncmp
         value=value+cmpstoi(ic)*ceq%complist(ic)%chempot(1)
      enddo
! >>>> subtract reference state: i.e. calculate G for the phase with 
! just this constituent
      endmember(1)=indices(2)
      call calcg_endmember(indices(1),endmember,gref,ceq)
      if(gx%bmperr.ne.0) goto 1000
      value=value-gref*ceq%rtn
!      write(*,511)'25C refstate: ',endmember(1),indices(1),gref,value
511   format(a,2i3,6(1pe14.6))
! possibly convert to AC or LNAC
      goto 700
   else
! MU(i) should be in position i, not indexed by splink ??
! loop through components, different components in each equilibrium
!      do ic=1,noofel
!         write(*,*)'state var value: ',indices(1),ceq%complist(ic)%splink,&
!              ceq%complist(ic)%chempot(1)
!         if(indices(1).eq.ceq%complist(ic)%splink) then
      if(indices(1).le.0 .or. indices(1).gt.noofel) then
!         write(*,*)'Asking for nonexisting chemical potential'
         gx%bmperr=4171; goto 1000
      endif
! iref=0 is default i.e. mean MU and iref<0 should mean MUS
! If a component has a defined reference state that is in complist(indices(1))
!      write(*,*)'25C Reference state: ',iref,ceq%complist(indices(1))%phlink
      if(iref.eq.0 .and. ceq%complist(indices(1))%phlink.ne.0) then
!      if(iref.eq.1) then
! phlink is phase, endmember is enmember, tpref<0 means current T
! we should also have a stoichiometry factor ??
         endmember(1)=indices(2)
         aref=ceq%tpval(1)
         if(ceq%complist(indices(1))%tpref(1).gt.zero) then
! reference state is at a fixed T, negative tpref(1) means current T
            ceq%tpval(1)=ceq%complist(indices(1))%tpref(1)
         endif
!         write(*,*)'25C calling calcg_endmember: ',&
!              ceq%complist(indices(1))%phlink,&
!              ceq%complist(indices(1))%endmember
         call calcg_endmember(ceq%complist(indices(1))%phlink,&
              ceq%complist(indices(1))%endmember,gref,ceq)
         if(gx%bmperr.ne.0) goto 1000
         ceq%tpval(1)=aref
         aref=ceq%complist(indices(1))%chempot(1)
         value=ceq%complist(indices(1))%chempot(1)-gref*ceq%rtn
!         write(*,513)'25C gref: ',indices(1),value,aref,gref*ceq%rtn
513      format(a,i3,5(1pe14.6))
      else
! this value should always be referenced to SER
! the value in chempot(2) is probably redundant after this change
         value=ceq%complist(indices(1))%chempot(1)
      endif
!      write(*,*)'25C chempot: ',indices(1),&
!           ceq%complist(indices(1))%chempot(1),&
!           ceq%complist(indices(1))%chempot(2)
      goto 700
   endif
! convert from MU to AC or LNAC if necessary
700 continue
!   ceq%rtn=globaldata%rgas*ceq%tpval(1)
   if(istv.eq.4) then
! AC = exp(mu/RT)
      value=exp(value/ceq%rtn)
   elseif(istv.eq.5) then
! LNAC = mu/RT
      value=value/ceq%rtn
   endif
!-----------------------------------------------------------------
1000 continue
   return
1100 continue
   gx%bmperr=4078
!   write(*,*)'State variable value not implemented yet'
   goto 1000
 end subroutine state_variable_val3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!-\begin{verbatim}
 subroutine state_var_value_derivative_old(svr1,svr2,value,ceq)
! THIS SUBROUTINE MOVED TO MINIMIZER
! subroutine state_var_value_derivative(istv,indices,iref,iunit,&
!      istv2,indices2,iref2,iunit2,value,ceq)
! calculates a state variable value derivative NOT IMPLEMENTED YET
! istv and istv2 are state variable type (integer)
! indices and indices2 are possible specifiers
! iref and iref2 are possible reference state
! iunit and iunit2 are units, (K, oC, J, cal etc)
! value is calculated value
! ceq is current equilibrium
   implicit none
   TYPE(gtp_state_variable), pointer :: svr1,svr2
   TYPE(gtp_equilibrium_data) :: ceq
!   integer :: istv,iref,iunit,istv2,iref2,iunit2
!   integer, dimension(4) :: indices,indices2
   double precision value
!-\end{verbatim}
!
   value=zero
   write(*,17)svr1%statevarid,svr1%argtyp,svr2%statevarid,svr2%argtyp
17 format('25C: state_var_value_derivative: ',10i4)
! this must be calculated in the minimizer
!   call meq_state_var_value_derivative(svr1,svr2,value,ceq)
   gx%bmperr=4078
1000 continue
   return
 end subroutine state_var_value_derivative_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calculate_reference_state(kstv,iph,ics,aref,ceq)
! Calculate the user defined reference state for extensive properties
! kstv is the typde of property: 1 U, 2 S, 3 V, 4 H, 5 A, 6 G
! It can be phase specific (iph.ne.0) or global (iph=0)
   implicit none
   integer kstv,iph,ics
   double precision aref
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! kstv=1  2  3  4  5  6 other values cared for elsewhere
!      U  S  V  H  A  G
   integer iel,phref
   double precision gref(6),bref(6),xmol(maxel),wmass(maxel),xxx(6)
   double precision tmol,tmass,bmult
!
!   write(*,*)'Reference states not implemented yet'; goto 1000
!   write(*,*)'25C reference state:',kstv,iph,ics
   if(kstv.lt.1 .or. kstv.gt.6) then
!      write(*,*)'No reference state for kstv: ',kstv
      goto 1000
   endif
   aref=zero
   bref=zero
   gref=zero
   xxx=zero
! loop for all components to extract the value of their reference states
! Multiply that with the overall composition (iph=0) or the phase composition
   xmol=zero
   do iel=1,noofel
! this is the reference phase for component iel
      phref=ceq%complist(iel)%phlink
      if(phref.gt.0) then
! special endmember call that returns G, G.T, G.P, G.T.T, G.T.P and G.P.P
!         write(*,73)'25C R state: ',iel,phref,ceq%complist(iel)%endmember
73       format(a,2i3,2x,10i4)
         call calcg_endmember6(phref,ceq%complist(iel)%endmember,gref,ceq)
         if(gx%bmperr.ne.0) goto 1000
         if(iph.gt.0) then
! multiply with mole fractions of phase iph,ics
            call calc_phase_molmass(iph,ics,xmol,wmass,tmol,tmass,bmult,ceq)
         else
! multiply with overall mole fractions
            call calc_molmass(xmol,wmass,tmol,tmass,ceq)
         endif
! note xxx, bref and gref are arrays
         xxx=bref+xmol(iel)*gref
!         write(*,70)'25Crs: ',bref,gref,xxx,(xmol(ij),ij=1,noofel)
70       format(a,6(1pe12.4)/,2(7x,6e12.4/),8(0pF8.4))
         bref=xxx
      else
! this is not really needed, it is bref that is used below
         gref=zero
      endif
   enddo
! calculate the correct correction depending on kstv
   if(kstv.eq.1) then
! U = G - T*G.T - P*G.P
      aref=bref(1)-ceq%tpval(1)*bref(2)-ceq%tpval(2)*bref(3)
   elseif(kstv.eq.2) then
! S = - G.T
      aref=-bref(2)
      
   elseif(kstv.eq.3) then
! V
      aref=bref(3)
      
   elseif(kstv.eq.4) then
! H = G - T*G.T
      aref=bref(1)-ceq%tpval(1)*bref(2)
      
   elseif(kstv.eq.5) then
! A = G - P*G.P
      aref=bref(1)-ceq%tpval(2)*bref(3)
      
   elseif(kstv.eq.6) then
! G
      aref=bref(1)
   endif
!   write(*,75)kstv,aref
75 format('25C ref:',i3,6(1pe12.4))
1000 continue
   return
 end subroutine calculate_reference_state

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 subroutine sortinphtup(n,m,xx)
! subroutine to sort the values in xx which are in phase and compset order
! in phase tuple order.  This is needed by the TQ interface
! The number of values belonging to the phase is m (for example composition)
   integer n,m
   double precision xx(n*m)
!
   integer iz,jz,kz,lz,lokph,aha
   double precision, dimension(:), allocatable :: dum
! I assume the values are NP(*), maybe there are other cases ...
   allocate(dum(n*m))
   kz=0
   do iz=1,noofph
      lokph=phases(iz)
      do jz=1,noofcs(lokph)
! in xx the values are sequentially for all composition sets for this phase
! But they should be stored in tuple order and compset 2 etc comes at the end
! the index to the tuple is in %phtups
! phlista(lokph)%linktocs(jz) is index of phase_varres record for compset
! firsteq%phase_varres(..)%phtupx is index of phase tuple for compset
! There can be m values (for example compositions) for each phase
         aha=(firsteq%phase_varres(phlista(lokph)%linktocs(jz))%phtupx-1)*m
         do lz=1,m
            dum(aha+lz)=xx(kz+lz)
         enddo
         kz=kz+m
      enddo
   enddo
   xx=dum
   deallocate(dum)
1000 continue
   return
 end subroutine sortinphtup

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

