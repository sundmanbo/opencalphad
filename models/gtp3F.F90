!
! gtp3E included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     10. State variable manipulations
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine get_stable_state_var_value(statevar,value,encoded,ceq)
! called with a state variable character
! If the state variable includes a phase it checks if the phase is stable ...
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   character statevar*(*),encoded*(*)
   double precision value
!\end{verbatim} %+
   integer lokcs,ics,ip
!   type(gtp_state_variable), pointer :: svr
   type(gtp_state_variable), target :: svrvar
   type(gtp_state_variable), pointer :: svr
   character modstatevar*28
!
! memory leak fix
   svr=>svrvar
   call decode_state_variable(statevar,svr,ceq)
   if(gx%bmperr.ne.0) goto 1000
! check if state variable inclused a phase
! argtyp=0: no arguments
! argtyp=1: component
! argtyp=2: phase+compset
! argtyp=3: phase+compset+component
! argtyp=4: phase+compset+constituent
   modstatevar=statevar
!   write(*,*)'3F stable: ',modstatevar,svr%argtyp,phlista(svr%phase)%noofcs
   if(svr%argtyp.eq.2) then
! if compset > 1 specified do nothing
      if(svr%compset.ne.1) goto 1000
!      svr%phase,svr%compset
      lokcs=phlista(svr%phase)%linktocs(svr%compset)
!      write(*,*)'3F phase: ',svr%compset,phlista(svr%phase)%noofcs,&
!           ceq%phase_varres(lokcs)%phstate,PHENTSTAB
      if(ceq%phase_varres(lokcs)%phstate.ne.PHENTSTAB) then
! phase+compset is not stable, chek if there is other stable compset
         loop: do ics=1,phlista(svr%phase)%noofcs
            lokcs=phlista(svr%phase)%linktocs(ics)
!            write(*,*)'3F looping: ',ics,lokcs,ceq%phase_varres(lokcs)%phstate
            if(ceq%phase_varres(lokcs)%phstate.eq.PHENTSTAB) then
! add a composition set index after phase name
               ip=index(modstatevar,',')
               if(ip.eq.0) ip=index(modstatevar,')')
! maybe there is a #1 ??
               if(modstatevar(ip-2:ip-2).eq.'#') then
                  modstatevar(ip-1:ip-1)=char(ics+ichar('0'))
               else
                  modstatevar(ip:)='#'//char(ics+ichar('0'))
                  modstatevar(ip+2:)=statevar(ip:)
               endif
!               write(*,*)'3F Modfied statevar: ',modstatevar
               exit loop
            endif
         enddo loop
      endif
   endif
   call get_state_var_value(modstatevar,value,encoded,ceq)
1000 continue
! possible memory leak
   nullify(svr)
   return
 end subroutine get_stable_state_var_value

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine get_state_var_value(statevar,value,encoded,ceq)
! called with a state variable character
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   character statevar*(*),encoded*(*)
   double precision value
!\end{verbatim}
!   integer indices(4)
   integer iunit,ip,lrot,mode
! memory leak
   type(gtp_state_variable), target :: svrvar
   type(gtp_state_variable), pointer :: svr
   character actual_arg(2)*16,name*16
!
!   write(*,*)'3F In get_state_variable_value: ',statevar
   iunit=0
   svr=>svrvar
   call decode_state_variable(statevar,svr,ceq)
!   write(*,20)statevar(1:len_trim(statevar)),svr%oldstv,svr%norm,&
!        svr%argtyp,svr%component
20  format('3F gsvv 1: ',a,' : ',4i3)
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
!         write(*,*)'3F Neither state variable or symbol'
         goto 1000
      else
! get the value of the symbol, may involve other symbols and state variablse
! The actual_arg is a facility not yet implemented and not allowed here
! if mode=0 the stored value may be used, mode=1 always evaluate
!         write(*,*)'3F Found function: ',lrot
         actual_arg=' '
         mode=1
! this is OK if it is not a derivative
         value=evaluate_svfun_old(lrot,actual_arg,mode,ceq)
         if(gx%bmperr.eq.4217) goto 1000
         encoded=name
      endif
   else
! it is a real state variable
!      write(*,*)'3F calling state_variable_val'
      call state_variable_val(svr,value,ceq)
      if(gx%bmperr.ne.0) goto 1000
      ip=1
      encoded=' '
      call encode_state_variable(encoded,ip,svr,ceq)
      if(gx%bmperr.ne.0) then
         write(*,*)'3F encode error: ',trim(encoded),gx%bmperr
         gx%bmperr=0; encoded='dummy'
      endif
   endif
1000 continue
! possible memory leak
   nullify(svr)
   return
 end subroutine get_state_var_value

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine get_many_svar(statevar,values,mjj,kjj,encoded,ceq)
! called with a state variable name with woldcards allowed like NP(*), X(*,CR)
! mjj is dimension of values, kjj is number of values returned
! encoded used to specify if phase data in phasetuple order ('Z')
! >>>> BIG question: How to do with phases that are note stable?
! If I ask for w(*,Cr) I only want the fraction in stable phases
! but whenthis is used for GNUPLOT the values are written in a matix
! and the same column in that phase must be the same phase ...
! so I have to have the same number of phases from each equilibria.
!
! CURRENTLY if x(*,*) and x(*,A) mole fractions only in stable phases
! Proposal: use * for all phases, use *S o $ or something else for all stable??
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   character statevar*(*),encoded*(*)
   double precision values(*)
   integer mjj,kjj
!\end{verbatim}
   integer indices(4),modind(4)
   double precision xnan,xxx
   integer jj,lokph,lokcs,k1,k2,k3,iref,jl,iunit,istv,enpos
! memory leak
   type(gtp_state_variable), target :: svrvar
   type(gtp_state_variable), pointer :: svr
!   logical phtupord
! calculate the NaN bit pattern
   xnan=0.0d0
!   xnan=0.0d0/xnan
   encoded=' '
   enpos=1
   if(gx%bmperr.ne.0) then
      write(*,*)'3F Error entering get_many_svar ',gx%bmperr,xnan
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
!   write(*,*)'3Y gmv 1: ',trim(statevar)
   svr=>svrvar
   call decode_state_variable(statevar,svr,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'3F Failed decode statevar in get_many_svar',gx%bmperr
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
   endif
!
!   write(*,20)istv,indices,iref,iunit,gx%bmperr
20 format('3F many 1: ',i5,4i4,3i7)
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
            call encode_state_variable3(encoded,enpos,istv,indices,&
                 iunit,iref,ceq)
            if(gx%bmperr.ne.0) goto 1000
            enpos=enpos+1
            call state_variable_val3(istv,indices,iref,&
                 iunit,values(jj),ceq)
            if(gx%bmperr.ne.0) goto 1000
         elseif(indices(3).eq.-1) then
! loop for components, indices 1+2 must be phase+compset
            do k3=1,noofel
               indices(3)=k3
               jj=jj+1
               if(jj.gt.mjj) goto 1100
               call encode_state_variable3(encoded,enpos,istv,indices,&
                    iunit,iref,ceq)
               if(gx%bmperr.ne.0) goto 1000
               enpos=enpos+1
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
               call encode_state_variable3(encoded,enpos,istv,indices,&
                    iunit,iref,ceq)
               if(gx%bmperr.ne.0) goto 1000
               enpos=enpos+1
               call state_variable_val3(istv,indices,iref,&
                    iunit,values(jj),ceq)
               if(gx%bmperr.ne.0) goto 1000
            enddo
         else
! indices(3) must be -2, -1 or >=0 so if we are here there is an error
            write(*,17)'3F Illegal set of indices 1',(indices(jl),jl=1,4)
17          format(a,4i4)
            gx%bmperr=4317; goto 1000
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
               call encode_state_variable3(encoded,enpos,istv,indices,&
                    iunit,iref,ceq)
               if(gx%bmperr.ne.0) goto 1000
               enpos=enpos+1
! if composition set not stable so return NaN (in xnan)
               if(test_phase_status(indices(2),indices(3),xxx,ceq).le. &
                    PHENTUNST) then
                  values(jj)=xnan
               elseif(ceq%phase_varres(lokcs)%dgm.lt.zero) then
! the phase must not have negative driving force
                  values(jj)=xnan
               else
! problem that get_many returns values for unstable phases
                  call state_variable_val3(istv,indices,iref,&
                       iunit,values(jj),ceq)
                  if(gx%bmperr.ne.0) goto 1000
23                format(a,2i3,2(1pe14.6))
               endif
            enddo
         enddo
      else
! if indices(1)>=0 then indices(2) must be -3 or >=0, so if here it is error
         write(*,17)'3F Illegal set of indices 2',(indices(jl),jl=1,4)
         gx%bmperr=4317; goto 1000
      endif
   elseif(indices(1).eq.-1) then
! loop for component as first indices, 2+3 can be fix phase+compset
! NOTE: loop for x(*,*) is below with indices(1).eq.-3
!      write(*,*)'3F indices: ',indices
      if(indices(2).ge.0) then
         do k1=1,noofel
            indices(1)=k1
            jj=jj+1
            if(jj.gt.mjj) goto 1100
            call encode_state_variable3(encoded,enpos,istv,indices,&
                 iunit,iref,ceq)
            if(gx%bmperr.ne.0) goto 1000
            enpos=enpos+1
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
                  call encode_state_variable3(encoded,enpos,istv,indices,&
                       iunit,iref,ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  enpos=enpos+1
                  if(test_phase_status(indices(2),indices(3),xxx,ceq).le. &
                       PHENTSTAB) then
! xnan means "no value"
                     values(jj)=xnan
                  elseif(ceq%phase_varres(lokcs)%dgm.lt.zero) then
! the phase must not have negative driving force
                     values(jj)=xnan
                  else
! here the value is extracted
                     call state_variable_val3(istv,indices,iref,&
                          iunit,values(jj),ceq)
                     if(gx%bmperr.ne.0) goto 1000
                  endif
               enddo
            enddo
         enddo
      else
! if we come here it must be an error
         write(*,17)'3F Illegal set of indices 3',(indices(jl),jl=1,4)
         gx%bmperr=4317; goto 1000
      endif
   elseif(indices(1).eq.-3) then
! loop for phase+compset as indices(1+2)
! here we must be careful not to destroy original indices, use modind
!      write(*,*)'3F get_many NP(*) etc 1: ',gx%bmperr,indices(3)
!      write(*,*)'Loop for many phases',indices(1)
      phloop: do k1=1,noofph
         modind(1)=k1
         modind(2)=0
         call get_phase_record(modind(1),lokph)
!         write(*,19)'3F test 17',modind,gx%bmperr,xnan
         if(gx%bmperr.ne.0) goto 1000
         csloop: do k2=1,phlista(lokph)%noofcs
            modind(2)=k2
            call get_phase_compset(modind(1),modind(2),lokph,lokcs)
            if(gx%bmperr.ne.0) goto 1000
            if(indices(3).eq.0) then
! This is typically listing of NP(*) for all phases
               modind(3)=0
               call encode_state_variable3(encoded,enpos,istv,modind,&
                    iunit,iref,ceq)
               if(gx%bmperr.ne.0) goto 1000
               enpos=enpos+1
               jj=jj+1
               if(jj.gt.mjj) goto 1100
               if(ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) then
                  values(jj)=xnan
               else
                  call state_variable_val3(istv,modind,iref,&
                       iunit,values(jj),ceq)
                  if(gx%bmperr.ne.0) goto 1000
               endif
!            elseif(ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) then
!               call encode_state_variable3(encoded,enpos,istv,modind,&
!                    iunit,iref,ceq)
!               if(gx%bmperr.ne.0) goto 1000
!               enpos=enpos+1
!               values(jj)=xnan
            elseif(indices(3).gt.0) then
! This is typically listing of w(*,cr), only in stable range of phases
               modind(3)=indices(3)
               call get_phase_compset(modind(1),modind(2),lokph,lokcs)
               if(gx%bmperr.ne.0) goto 1000
               call encode_state_variable3(encoded,enpos,istv,modind,&
                    iunit,iref,ceq)
               if(gx%bmperr.ne.0) goto 1000
               enpos=enpos+1
               jj=jj+1
               if(jj.gt.mjj) goto 1100
               if(ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) then
! if phase is unstable set dummy value
                  values(jj)=xnan
               else
                  call state_variable_val3(istv,modind,iref,&
                       iunit,values(jj),ceq)
               endif
!               write(*,73)'3F listing w(*,A): ',istv,modind,iref,iunit,&
!                    ceq%phase_varres(lokcs)%phstate,jj,values(jj)
73             format(a,i5,2x,4i3,2x,2i4,2i5,1pe12.4)
               if(gx%bmperr.ne.0) goto 1000
            elseif(ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) then
! loop for all components or constitunets of stable phases
! Maybe it should be included to have same number of values in all ranges?
               cycle csloop
            elseif(indices(3).eq.-1) then
! loop for all components of all phases, skip unstable phases
               elloop: do k3=1,noofel
                  modind(3)=k3
                  call encode_state_variable3(encoded,enpos,istv,modind,&
                       iunit,iref,ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  enpos=enpos+1
                  jj=jj+1
                  if(jj.gt.mjj) goto 1100
                  call state_variable_val3(istv,modind,iref,&
                       iunit,values(jj),ceq)
                  if(gx%bmperr.ne.0) goto 1000
               enddo elloop
            elseif(indices(3).eq.-2) then
! loop for constituents of all phases
               conloop: do k3=1,phlista(lokph)%tnooffr
                  modind(3)=k3
                  call encode_state_variable3(encoded,enpos,istv,modind,&
                       iunit,iref,ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  enpos=enpos+1
                  jj=jj+1
                  if(jj.gt.mjj) goto 1100
                  call state_variable_val3(istv,modind,iref,&
                       iunit,values(jj),ceq)
                  if(gx%bmperr.ne.0) goto 1000
               enddo conloop
            else
! error if here
               write(*,17)'3F Illegal set of indices 4',(indices(jl),jl=1,4)
               gx%bmperr=4317; goto 1000
            endif
            if(gx%bmperr.ne.0) then
               write(*,19)'3F error 3',modind,gx%bmperr
19             format(a,4i4,i7)
               goto 1000
            endif
         enddo csloop
      enddo phloop
   else
! error if here
      write(*,17)'3F Illegal set of indices 5',(indices(jl),jl=1,4)
      gx%bmperr=4317; goto 1000
   endif
1000 continue
! possible memory leak, nullify does not release memory
   nullify(svr)
   kjj=jj
   return
1100 continue
   write(*,*)'3F Overflow in array to get_state_variables'
   gx%bmperr=4317; goto 1000
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
!   write(*,*)'3F in decode ',trim(statevar),istv
! initiate svr internal variables
   deblist=.FALSE.
!   deblist=.TRUE.
   if(ocv()) deblist=.TRUE.
   if(deblist) write(*,*)'3F entering decode_statevariable: ',&
        statevar(1:len_trim(statevar))
!   write(*,*)'3F svr allocated'
! memory leak
!   allocate(svr)
!   write(*,*)'3F svr assignment start'
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
!   write(*,*)'3F svr assignment end'
!
! For wildcard argument "*" return:
! -1 for element or component
! -2 for species or constituent
! -3 for phase
! -4 for composition set
   istv=-1
   indices=0
! iref=0 means user defined reference state
   iref=0
! unit is not implemented (can apply to T, P, V, mass, etc)
   iunit=0
   iph=0
   ics=0
   norm=0
! local character for state variable
   lstate=statevar
   call capson(lstate)
   if(deblist) write(*,*)'3F decode_state_var 1: ',lstate(1:20)
! compare first character
   ch1=lstate(1:1)
!   write(*,*)'3F decoding: ',trim(lstate),is,' ',ch1
   do is=1,noos
      if(ch1.eq.svid(is)(1:1)) goto 50
   enddo
! it may be a property, parameter identifier
   if(deblist) write(*,*)'3F jump to 600!'
   goto 600
!------------------------------------------------------------
50 continue
   if(deblist) write(*,*)'3F dsv 1: ',is,ch1
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
105 format('3F dsv 4: ',3i4)
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
   if(deblist) write(*,*)'3F lstate: ',lstate(1:20),jp
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
      if(deblist) write(*,*)'3F Normalize 1: ',is,jp,ch1,norm
   endif nomalize
!---------------------------------------------------------------------
! reference state can be specified by an S for SER
! If no S the user specified reference states applies
   if(lstate(jp:jp).eq.'S') then
      jp=jp+1
      iref=-1
   elseif(lstate(jp:jp).eq.'R') then
      write(*,*)'Ignoring suffix "R" on ',trim(statevar),&
           ', user reference is default'
      jp=jp+1
   endif
!---------------------------------------------------------------------
! extract arguments if any. If arguments then lstate(jp:jp) should be (
! Typically H(fcc#2), N(Cr), BP(fcc), Y(sigma#2,cr#3), TC(BCC#2)
!300 continue
   if(deblist) write(*,*)'3F args: ',jp,lstate(1:jp+10)
   narg=0
   args: if(lstate(jp:jp).eq.'(') then
      kp=index(lstate,')')
      if(kp.le.0) then
         if(deblist) write(*,110)'3F dsv 5: ',is,jp,kp,lstate(1:20)
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
      if(deblist) write(*,*)'3F maybe a property ',is
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
   if(deblist) write(*,*)'3F at 600: ',propsym(1:len_trim(propsym)),typty
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
      write(*,*)'3F Unknown state variable or property',typty
      gx%bmperr=4318; goto 1000
   endif
   svr%oldstv=istv
   svr%statevarid=istv
   svr%constituent=indices(3)
   if(deblist) write(*,611)'3F Property: ',is,istv,typty,indices
611 format(a,10i4)
!------------------------------------------------
1000 continue
! accept the current istv as svr%oldstv, store a suffix S on MU as phref<0
   svr%oldstv=istv
   svr%phref=iref
   if(deblist) write(*,1001)'3F exit decode: ',istv,(indices(is),is=1,4),&
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
!
 end subroutine decode_state_variable3 !allocate

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_phase_molmass(iph,ics,xmol,wmass,totmol,totmass,amount,ceq)
! calculates mole fractions and mass fractions for a phase#set
! xmol and wmass are fractions of components in mol or mass
! totmol is total number of moles and totmass total mass of components.
! amount is number of moles of components per formula unit.
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer iph,ics
   double precision, dimension(*) :: xmol,wmass
   double precision amount,totmol,totmass
!\end{verbatim}
   integer ic,jc,lokph,lokcs,ll,iel,lokel,ie,kk,loksp,nspel
   integer compnos(maxspel)
   double precision as,yz,xsum,wsum,stoi(maxspel),smass,qsp
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
   if(ocv()) write(*,14)'3F cpm: ',iph,ics,lokph,lokcs
14 format(a,10i5)
   allsubl: do ll=1,phlista(lokph)%noofsubl
      as=ceq%phase_varres(lokcs)%sites(ll)
      allcons: do kk=1,phlista(lokph)%nooffr(ll)
         ic=ic+1
         if(.not.btest(ceq%phase_varres(lokcs)%constat(ic),CONSUS)) then
            yz=ceq%phase_varres(lokcs)%yfr(ic)
            loksp=phlista(lokph)%constitlist(ic)
            if(.not.btest(globaldata%status,GSNOTELCOMP)) then
! the elements are the components
               do iel=1,splista(loksp)%noofel
                  lokel=splista(loksp)%ellinks(iel)
                  ie=ellista(lokel)%alphaindex
                  if(ie.ne.0) then
                     xmol(ie)=xmol(ie)+&
                          as*yz*splista(loksp)%stoichiometry(iel)
                  endif
               enddo
            else
! when we have other components than the elements
! we must convert the element stoichiometry to component stoichiometry
!               write(*,*)'3F other components than elements'
               call get_species_component_data(loksp,nspel,compnos,stoi,&
                    smass,qsp,ceq)
               do iel=1,nspel
                  xmol(compnos(iel))=xmol(compnos(iel))+as*yz*stoi(iel)
               enddo
            endif
         endif
      enddo allcons
   enddo allsubl
! normallize, All ok here
!   write(*,713)'A',noofel,(xmol(iq),iq=1,noofel)
713 format('3F x:',a,i2,10f8.4)
!800 continue
   xsum=zero
   wsum=zero
! here xmol(i) is equal to the number of moles of element i per formula unit
! set wmass(i) to the mass of of element i per mole formula unit and sum
   do ic=1,noofel
!      wmass(ic)=xmol(ic)*ellista(elements(ic))%mass
      wmass(ic)=xmol(ic)*ceq%complist(ic)%mass
      xsum=xsum+xmol(ic)
      wsum=wsum+wmass(ic)
   enddo
!   write(*,713)'3F cpmm: ',noofel,xsum,(xmol(ic),ic=1,noofel)
   do ic=1,noofel
      xmol(ic)=xmol(ic)/xsum
      wmass(ic)=wmass(ic)/wsum
   enddo
!   write(*,713)'3F cpmm: ',noofel,xsum,(xmol(ic),ic=1,noofel)
! This is the current number of formula unit of the phase, zero if not stable
   amount=ceq%phase_varres(lokcs)%amfu
! ceq%phase_varres(lokcs)%abnorm(1) is moles atoms for one formula unit
! ceq%phase_varres(lokcs)%abnorm(2) is mass for one formula unit
   totmol=amount*xsum
   totmass=amount*wsum
!   write(*,717)'3F z:',noofel,lokcs,totmol,totmass,amount,&
!        wsum,ceq%phase_varres(lokcs)%abnorm(2)
717 format(a,i3,i6,6(1pe12.4))
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
!298 format('3F: ',6(1pe12.4))
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
   TYPE(gtp_equilibrium_data),pointer :: ceq
!\end{verbatim}
   integer ic,lokph,lokcs,ll,kk,loksp,lokel,iel,ie,compnos(maxspel),nspel
   double precision as,yz,xsum,smass,qsp,stoi(maxspel)
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
            if(.not.btest(globaldata%status,GSNOTELCOMP)) then
! the elements are the components
               do iel=1,splista(loksp)%noofel
                  lokel=splista(loksp)%ellinks(iel)
                  ie=ellista(lokel)%alphaindex
                  if(ie.ne.0) then
                     xmol(ie)=xmol(ie)+&
                          as*yz*splista(loksp)%stoichiometry(iel)
                  endif
               enddo
            else
! when we have other components than the elements
! we must convert the element stoichiometry to component stoichiometry
!               write(*,*)'3F other components than elements'
               call get_species_component_data(loksp,nspel,compnos,stoi,&
                    smass,qsp,ceq)
               do iel=1,nspel
                  xmol(compnos(iel))=xmol(compnos(iel))+as*yz*stoi(iel)
               enddo
            endif

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
   TYPE(gtp_equilibrium_data), pointer :: ceq
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
!               write(*,17)'3F amult:',iph,ics,am,amult,tmol,tmass
!               write(*,18)'3F x0: ',(xph(ic),ic=1,noofel)
!               write(*,18)'3F w0: ',(wph(ic),ic=1,noofel)
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
   if(totmass.gt.zero) then
      do ic=1,noofel
         xmol(ic)=xmol(ic)/totmol
         wmass(ic)=wmass(ic)/totmass
      enddo
!      write(*,21)'3F x1: ',totmol,(xmol(ic),ic=1,noofel)
!      write(*,21)'3F w2: ',totmass,(wmass(ic),ic=1,noofel)
21    format(a,1pe12.4,8(0pF9.5))
!   else
!      write(*,*)'There is no mass at all in the system!'
!      gx%bmperr=4185; goto 1000
   endif
!   write(*,21)'3F x1: ',totmol,(xmol(ic),ic=1,noofel)
!   write(*,21)'3F w1: ',totmass,(wmass(ic),ic=1,noofel)
!   else
! this is not an error if no calculation has been made
!      write(*,28)'3F: calc_molmass: No mole fractions',totmol,totmass,xsum,&
!           (xmol(ic),ic=1,noofel)
28    format(a,3(1pe12.4)/'3F. ',10f7.4)
!      gx%bmperr=4185; goto 1000
!   endif
!   wsum=zero
!   do ic=1,noofel
!      wmass(ic)=xmol(ic)*ellista(elements(ic))%mass
!      wsum=wsum+wmass(ic)
!      write(*,44)'3F cmm4: ',ic,xmol(ic),wmass(ic),&
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
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision props(5)
!\end{verbatim}
   integer lokph,lokcs,ics
   double precision am
   if(gx%bmperr.ne.0) write(*,*)'3F error entering sumprops ',gx%bmperr
   props=zero
   allph: do lokph=1,noofph
!      write(*,*)'3F sumprops: ',lokph
      if(.not.btest(phlista(lokph)%status1,phhid)) then
!         lokcs=phlista(lokph)%cslink
         allcs: do ics=1,phlista(lokph)%noofcs
! phase_varres(lokcs)%amfu is the amount formula units of the phase
! phase_varres(lokcs)%abnorm(1) is the moles of real atoms/formula unit
! am is the number of moles of real atoms of the phase
            lokcs=phlista(lokph)%linktocs(ics)
! skip phases that are not entered
            if(ceq%phase_varres(lokcs)%phstate.eq.phdorm) cycle allcs
            am=ceq%phase_varres(lokcs)%amfu*&
                 ceq%phase_varres(lokcs)%abnorm(1)
!            write(*,7)'3F sumprops 2: ',lokph,lokcs,am,&
!                 ceq%phase_varres(lokcs)%abnorm(1),&
!                 ceq%phase_varres(lokcs)%abnorm(2),props(5)
7           format(a,2i5,6(1pe12.4))
! valgrind complains this jump if for an uninitiallized valiable ??
            if(am.gt.zero) then
! properties are G, G.T=-S, G.P=V and moles and mass of real atoms
! Note gval(*,1) is per mole formula unit and ceq%phase_varres(lokcs)%abnorm(1)
! is the number of real atoms per formula unit
!               write(*,13)'3F props1:',lokcs,props(1),props(2),props(3),am,&
!                    ceq%phase_varres(lokcs)%abnorm(1)
               props(1)=props(1)+am*ceq%phase_varres(lokcs)%gval(1,1)/&
                    ceq%phase_varres(lokcs)%abnorm(1)
!               write(*,10)props(1),am,ceq%phase_varres(lokcs)%gval(1,1),&
!                    ceq%phase_varres(lokcs)%abnorm(1)
10             format('3F props: ',6(1pe12.4))
               props(2)=props(2)+am*ceq%phase_varres(lokcs)%gval(2,1)/&
                    ceq%phase_varres(lokcs)%abnorm(1)
               props(3)=props(3)+am*ceq%phase_varres(lokcs)%gval(3,1)/&
                    ceq%phase_varres(lokcs)%abnorm(1)
               props(4)=props(4)+am
!               write(*,13)'3F props2:',lokcs,props(1),props(2),props(3),&
!                    props(4),ceq%phase_varres(lokcs)%gval(3,1)
13             format(a,i3,6(1pe12.4))
! ceq%phase_varres(lokcs)%abnorm(2) should be the current mass
! %abnorm(2) is actual mass, its should be multiplied with %amfu, not am!!
! This value is calculated in set_constitution ... check there if problems
               props(5)=props(5)+ceq%phase_varres(lokcs)%amfu*&
                    ceq%phase_varres(lokcs)%abnorm(2)
!               write(*,75)'3F sumprops: ',lokcs,am,&
!                    ceq%phase_varres(lokcs)%abnorm(2),props(5)
75             format(a,i4,6(1pe12.4))
!               write(*,11)'3F sumprops: ',lokcs,props(1),props(4),props(5),&
!                    ceq%phase_varres(lokcs)%abnorm(2)
!               write(*,11)'3F sumprops ',lokcs,am,props(4),&
!                    ceq%phase_varres(lokcs)%abnorm(1)
11             format(a,i4,6(1pe12.4))
            endif
         enddo allcs
      endif
   enddo allph
1000 continue
   if(gx%bmperr.ne.0) write(*,*)'3F error exiting sumprops ',gx%bmperr
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
!   write(*,*)'3F ************* encode: '
   iref=svr%phref
   iunit=svr%unit
! if svr%oldstv>=10 then istv should be 10*(svr%oldstv-5)+svr%norm
!   if(svr%oldstv.ge.10) then
!      istv=10*(svr%oldstv-5)+svr%norm
!      write(*,*)'3F encode: ',svr%oldstv,svr%norm,istv
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
   integer jp,ics,kstv,iph,norm,sublat
   double precision mass
!
   character stsymb*60
   character*1, dimension(4), parameter :: cnorm=['M','W','V','F']
!
   sublat=0
!   write(*,*)'3F encode: ',istv
   if(istv.le.0) then
! this is a parameter property symbol: TC (-2), BM (-3), MQ&FE(FCC) (-4) etc
! translate to 21, 22, 23 ...
      kstv=19-istv
      goto 200
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
         stsymb(jp:jp)=')'
         jp=jp+1
      else
! always use composition set 1
         ics=1
         call get_phase_name(indices(1),ics,stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=','; jp=jp+1
!         call get_phase_constituent_name(indices(1),indices(2),&
!              stsymb(jp:),sublat)
! I am not sure if indices(2) is constituent numbered for each sublattice
! or numbered from the beginning, assume the latter !!
         call get_constituent_name(indices(1),indices(2),&
              stsymb(jp:),mass)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         if(sublat.gt.1) then
            stsymb(jp:)='#'//char(ichar('0')+sublat)//')'
            jp=jp+3
         else
            stsymb(jp:jp)=')'
            jp=jp+1
         endif
      endif
      goto 800
   endif potential
   if(istv.lt.10) then
!      write(*,*)'3F unknown potential'
      gx%bmperr=4158; goto 1000
   endif
! Extensive property has istv>=10
   norm=mod(istv,10)
   kstv=(istv+1)/10+5
!    write(*,3)'3F encode 3: ',kstv,indices
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
!   write(*,11)'3F esv 7: ',istv,kstv,indices
11 format(a,10i4)
   stsymb=svid(kstv)
   jp=len_trim(stsymb)+1
!   write(*,*)'3F norm 1A: ',kstv,norm
   if(kstv.le.16 .or. kstv.eq.18) then
      if(norm.gt.0 .and. norm.le.4) then
!         write(*,*)'3F norm 1B: ',kstv,norm
         stsymb(jp:jp)=cnorm(norm)
         jp=jp+1
      elseif(norm.ne.0) then
!         write(*,*)'3F norm 1C: ',kstv,norm
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
!   write(*,*)'3F parameter property symbol: ',kstv,iph,ics
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
!         call get_phase_constituent_name(indices(1),indices(3),&
!              stsymb(jp:),sublat)
         call get_constituent_name(indices(1),indices(3),&
              stsymb(jp:),mass)
         if(gx%bmperr.ne.0) goto 1000
! sublattice is the last argument
         jp=len_trim(stsymb)+1
         if(sublat.gt.1) then
            stsymb(jp:)='#'//char(ichar('0')+sublat)//')'
            jp=jp+3
         else
            stsymb(jp:jp)=')'
            jp=jp+1
         endif
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
   if(ip+jp.gt.len(text)) then
      write(*,*)'State variable value output exceed character variable length'
      gx%bmperr=4319; goto 1000
   endif
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
   integer jp,ics,kstv,iph,norm,sublat
   integer, dimension(4) :: indices
   integer istv,ip,iunit,iref
   double precision mass
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
!         call get_phase_constituent_name(indices(1),indices(2),&
!              stsymb(jp:),sublat)
         call get_constituent_name(indices(1),indices(2),&
              stsymb(jp:),mass)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         if(sublat.gt.1) then
            stsymb(jp:)='#'//char(ichar('0')+sublat)//')'
            jp=jp+3
         else
            stsymb(jp:jp)=')'
            jp=jp+1
         endif
      endif
      stsymb(jp:jp)=')'
      goto 800
   endif potential
   if(istv.lt.10) then
!      write(*,*)'3F unknown potential'
      gx%bmperr=4158; goto 1000
   endif
! Extensive property has istv>=10
   norm=mod(istv,10)
   kstv=(istv+1)/10+5
!    write(*,3)'3F encode 3: ',kstv,indices
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
!    write(*,11)'3F esv 7: ',istv,kstv,indices
11  format(a,10i4)
   stsymb=svid(kstv)
   jp=len_trim(stsymb)+1
!   write(*,*)'3F norm 2: ',kstv,norm
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
!         call get_phase_constituent_name(indices(1),indices(3),&
!              stsymb(jp:),sublat)
         call get_constituent_name(indices(1),indices(3),&
              stsymb(jp:),mass)
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         if(sublat.gt.1) then
            stsymb(jp:)='#'//char(ichar('0')+sublat)//')'
            jp=jp+3
         else
            stsymb(jp:jp)=')'
            jp=jp+1
         endif
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
! searching for experimental bug
!   write(*,*)'3F state_variable_val: ',svr%statevarid,iref,iunit
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
      write(*,*)'3F state variable has illegal argtyp: ',svr%argtyp
      gx%bmperr=4320; goto 1000
   endif
!   write(*,910)'3F svv: ',istv,indices,iref,iunit,value
910 format(a,i3,2x,4i3,2i3,1pe14.6)
   call state_variable_val3(istv,indices,iref,iunit,value,ceq)
   if(gx%bmperr.ne.0) then
!      write(*,920)'3F error 7: ',gx%bmperr,istv,svr%oldstv,svr%argtyp
!920   format(a,i5,2x,2i15,i2)
!   else
!      write(*,*)'3F value: ',value
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
   double precision qsp,gref,spmass,rmult,tsave,rtn,spextra
   integer kstv,norm,lokph,lokcs,icx,jp,ncmp,ic,iprop,loksp,nspel,iq
   integer endmember(maxsubl),ielno(maxspel)
   value=zero
   ceq%rtn=globaldata%rgas*ceq%tpval(1)
!   write(*,10)'3F svval3: ',istv,indices,iref,iunit,gx%bmperr,value
10 format(a,i4,4i4,3i5,1PE17.6)
   potentials: if(istv.lt.0) then
! negative istv indicate parameter property symbols
      kstv=-istv
      goto 200
!      gx%bmperr=4097; goto 1000
   elseif(istv.ge.10) then
      goto 50
   elseif(istv.eq.1) then
! this is T
      value=ceq%tpval(1)
   elseif(istv.eq.2) then
! this is P
      value=ceq%tpval(2)
   elseif(istv.le.5) then
! the check of reference state state is made at label 500 
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
      write(*,10)'3F not impl: ',istv,indices,iref,iunit,gx%bmperr,value
      goto 1100
   else
! wrong or state variable not implemented
      write(*,10)'3F not impl: ',istv,indices,iref,iunit,gx%bmperr,value
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
   if(ocv()) write(*,51)'3F stv A: ',props
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
!      write(*,*)'3F stv B: ',vt,vp,amult
      if(indices(1).eq.0) then
! global value for the whole system
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
!         write(*,555)'3F pp: ',vg,props
555      format(a,6(1pe12.4))
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
! iref=0 means user defined reference state >>>> unfinished
!??????????????????????????????????
! UNFINISHED
! If O has reference state but no other elements then ignore refstate
! for integral quantities UNLESS all components has the same reference
! state ....
!         write(*,52)'3F Ref state:',iref,kstv,indices(1),indices(2),rmult
52       format(a,4i4,1pe12.4)
! IMPORTANT !!! calculate reference state may destroy valies in %gval
         call calculate_reference_state(kstv,indices(1),indices(2),aref,ceq)
         if(gx%bmperr.ne.0) goto 1000
! value here seems OK
!         write(*,53)'3F Reference state:',iref,aref,rmult
      elseif(iref.lt.0) then
         aref=zero
      else
! positive value of iref is undefined
         write(*,*)'3F Reference state undefined: ',iref
         aref=zero
      endif
! if phase specific the scaling for phase specific must be compensated
      aref=rmult*aref
!      write(*,53)'3F at kstv1: ',kstv,props,aref,div
53    format(a,i3,5(1PE12.3))
      kstv1: if(kstv.eq.1) then
! 1: U = G + TS - PV = G - T*G.T - P*G.P
         value=amult*(vg+vt*vs-vp*vv-aref)/div
      elseif(kstv.eq.2) then
! 2: S = -G.T
         value=amult*(vs-aref)/div
!         write(*,54)value,amult,vs,aref,div
54       format('3F svv: ',5(1pe12.4))
      elseif(kstv.eq.3) then
! 3: V = G.P
         value=amult*(vv-aref)/div
!         write(*,54)amult,vv,aref,div,value
      elseif(kstv.eq.4) then
! 4: H = G + TS = G - T*G.T
! Problem with vg here when reference state is set
         if(ocv()) write(*,177)'3F H:',vg+vt*vs,aref,amult,div,rmult
177      format(a,6(1pe12.4))
         value=amult*(vg+vt*vs-aref)/div
      elseif(kstv.eq.5) then
! 5: A = G - PV = G - P*G.P
         value=amult*(vg-vp*vv-aref)/div
      elseif(kstv.eq.6) then
! 6: G
!         write(*,177)'3F G:',vg,aref,amult,div
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
!         write(*,202)'3F svval DG:',lokcs,ceq%phase_varres(lokcs)%dgm,div
202      format(a,i5,2(1pe12.4))
         value=ceq%phase_varres(lokcs)%dgm/div
      elseif(kstv.eq.10) then
! 10: Q (stability, thermodynamic factor), not implemented
!         gx%bmperr=4081; goto 1000
         call calc_qf(lokcs,value,ceq)
!      else
!         write(*,*)'3F svval after 10:',kstv
      endif kstv1
      goto 1000
   endif le10
!----------------------------------------------------------------------
! here with kstv>10
! kstv=       11  12  13 
! state var:   N   B   Y   
   le12: if(kstv.le.12) then
! normallizing for N (kstv=11) and B (kstv=12)
!      write(*,88)'3F svv 12: ',indices(1),norm,props(4),props(5)
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
!         write(*,89)'3F svv, N or B: ',vn,div
89       format(a,5(1pe12.4))
         value=vn/div
      else
! one or two indices, overall of phase specific component amount
         if(indices(2).eq.0) then
! one index is component specific, N(comp.), B(CR) etc. Sum over all phases
! props(4) is amount of moles of components, props(5) is mass of components
            call calc_molmass(xmol,wmass,tmol,tmass,ceq)
            if(gx%bmperr.ne.0) goto 1000
!            write(*,89)'3F mm: ',tmol,tmass
!            write(*,93)'3F x: ',(xmol(icx),icx=1,noofel)
!            write(*,93)'3F w: ',(wmass(icx),icx=1,noofel)
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
!            write(*,92)'3F cpmm: ',indices(icx),tmol,tmass,bmult,&
!                 wmass(1),wmass(2),wmass(3)
92          format(a,i3,3(1pe12.4),3(0pF8.5))
         endif
         if(gx%bmperr.ne.0) goto 1000
!         write(*,13)'3F gsvv 19: ',norm,(xmol(iq),iq=1,noofel)
777      format('gsvv 77: ',10(f7.4))
         if(kstv.eq.11) then
! total moles of component
            vn=xmol(indices(icx))
            amult=tmol
! added next line 2015-08-20 to get correct N(sigma,mo)
            bmult=tmol
!            write(*,777)kstv,icx,indices(icx),norm,vn,amult,bmult
!777         format('3F N(i): ',4i4,3(1pe12.4))
         else
! total mass of component
            vn=wmass(indices(icx))
            amult=tmass
! added next line 2015-08-20 to get correct N(sigma,mo)
            bmult=tmass
         endif
!         write(*,13)'3F gsvv 8: ',norm,vn,amult,bmult,tmol,tmass
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
!            write(*,*)'3F x%: ',iunit,value
            if(iunit.eq.100) value=1.0D2*value
         elseif(norm.eq.2) then
! NW or W
            if(tmass.gt.zero) then
               value=amult*vn/tmass
            else
               value=zero
            endif
! problem when plotting w(*,C) for phase fix with 0 amount
!            value=wmass(indices(icx))
            value=vn
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
!   write(*,*)'3F svv3 at 200:',kstv
   select case(kstv)
   case default
      write(kou,*)'Unknown parameter identifier: ',kstv
!.......................................
   case(2:5,7,9:19) 
!-------------------------------------------------------------------
! 2: TC (Curie/Neel Temperature)
! 3: BM (Average Bohr magneton number)
! 4: CTA just Curie Temperature
! 5: NTA just Neel temperature
! 7: THET Debye or Einstein temperature
! 8: MQ& mobility
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
! 20: Flory-Huggins model parameter
!-------------------------------------------------------------------
!
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
   case(6,8,20) 
! 6: IBM& Individual Bohr magneton number
! 8: MQ& mobility value
! 20: FHV  Flory Huggins volume
!      write(*,*)'3F svv3 mob1: ',indices(1),indices(2),iprop
      call get_phase_compset(indices(1),indices(2),lokph,lokcs)
      if(gx%bmperr.ne.0) goto 1000
! property is kstv*100+indices(3) (constituent identifier)
      iprop=100*kstv+indices(3)
!      write(*,485)'3F svv3 mob2: ',indices(1),indices(2),iprop,&
!           ceq%phase_varres(lokcs)%nprop
485   format(a,2i3,10i5)
      find2: do jp=2,ceq%phase_varres(lokcs)%nprop
!         write(*,485)'3F calcprop: ',ceq%phase_varres(lokcs)%listprop(jp)
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
! here indices(2) is considered to specify a reference state ...???
!   write(*,502)'3F refstate 500: ',iref,indices(1),indices(2)
502 format(a,10i4)
   if(indices(2).ne.0) then
! This has nothing to do with reference state ... ??? see else link for that
! I wonder if this code is ever used ...
      lokph=phases(indices(1))
      loksp=phlista(lokph)%constitlist(indices(2))
! split the species in elements, convert to components, add chemical potentials
      call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp,spextra)
      if(gx%bmperr.ne.0) goto 1000
      if(qsp.ne.zero) then
!         write(*,*)'3F Cannot calculate potential of charged species'
         gx%bmperr=4159; goto 1000
      endif
! other components than elements not implemented
!      write(*,*)'3F converting to component',nspel,ielno(1),stoi(1)
      call elements2components(nspel,stoi,ncmp,cmpstoi,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,*)'3F converting to component',ncmp,cmpstoi(1)
      value=zero
      do ic=1,ncmp
         value=value+cmpstoi(ic)*ceq%complist(ic)%chempot(1)
      enddo
! >>>> subtract reference state: i.e. calculate G for the phase with 
! just this constituent.  Note indices(1) is phase record, change to index
!      write(*,*)'3F refphase: ',indices(1),phlista(indices(1))%alphaindex,value
      ic=phlista(indices(1))%alphaindex
      endmember(1)=indices(2)
!      write(*,*)'3F callcg_endmember 1: ',indices(1)
      call calcg_endmember(indices(1),endmember,gref,ceq)
      if(gx%bmperr.ne.0) goto 1000
      value=value-gref*ceq%rtn
!      write(*,511)'3F refstate: ',endmember(1),indices(1),gref,value
511   format(a,2i3,6(1pe14.6))
! possibly convert to AC or LNAC
      goto 700
   else
!      write(*,*)'3F elselink: ',indices
      if(indices(1).le.0 .or. indices(1).gt.noofel) then
!         write(*,*)'3F Asking for nonexisting chemical potential'
         gx%bmperr=4171; goto 1000
      endif
! iref=0 is default meaning user defined reference state,
! if iref<0 use SER as reference state, ignoring user defined reference state
!
! If a component has a defined reference state that is in complist(indices(1))
      if(iref.eq.0 .and. ceq%complist(indices(1))%phlink.ne.0) then
!         write(*,*)'3F Reference state: ',indices(1),indices(2),&
!              ceq%complist(indices(1))%phlink
! phlink is phase, endmember is enmember, tpref<0 means current T
! we should also have a stoichiometry factor ??
!         endmember(1)=indices(2)
         tsave=ceq%tpval(1)
         if(ceq%complist(indices(1))%tpref(1).gt.zero) then
! reference state is at a fixed T, negative tpref(1) means current T
            ceq%tpval(1)=ceq%complist(indices(1))%tpref(1)
         endif
!         write(*,*)'3F calling calcg_endmember: ',&
!              ceq%complist(indices(1))%phlink,&
!              ceq%complist(indices(1))%endmember
         ic=ceq%complist(indices(1))%phlink
!         ic=phlista(ic)%alphaindex
!         write(*,*)'3F refphase: ',indices(1),ic,phlista(indices(1))%alphaindex
!         ic=phlista(indices(1))%alphaindex
! the first index should be phase index, not location
! We may have to restore gval after this!!!
!         write(*,*)'3F callcg_endmember 2: ',-ic
         call calcg_endmember(-ic,ceq%complist(indices(1))%endmember,gref,ceq)
         if(gx%bmperr.ne.0) then
            write(*,*)'3F Error calculating refstate for chemical pot'
            goto 1000
         endif
! RT for current T
         rtn=globaldata%rgas*ceq%tpval(1)
         ceq%tpval(1)=tsave
         aref=ceq%complist(indices(1))%chempot(1)
!         value=ceq%complist(indices(1))%chempot(1)-gref*rtn
         value=aref-gref*rtn
!         write(*,513)'3F gref: ',indices(1),value/rtn,aref/rtn,gref,rtn
513      format(a,i3,5(1pe14.6))
         ceq%complist(indices(1))%chempot(2)=value
      else
! the value in chempot(1) should always be referenced to SER
! the value in chempot(2) should always be for the user reference
         value=ceq%complist(indices(1))%chempot(1)
      endif
!      write(*,*)'3F chempot: ',indices(1),&
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
!   write(*,*)'3F State variable value not implemented yet'
   goto 1000
 end subroutine state_variable_val3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_qf(lokcs,value,ceq)
! calculates eigenvalues of the second derivative matrix, stability function
! lokcs is index of phase_varres
! value calculated value returned
! ceq is current equilibrium
! For ionic liquid and charged crystalline phases one should
! calculate eigenvectors to find neutral directions.
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer :: lokcs
   double precision value
!\end{verbatim}
   integer info,nofc,nofc2,nmax,ii,jj,lokph,nsl,cc,rr,zz
   double precision, allocatable :: d2g(:),eigenv(:),work(:)
   double precision dummy(1,1)
   integer, allocatable :: skip(:)
! number of constituents
! ignore sublattices with single constituents ....
   nofc=size(ceq%phase_varres(lokcs)%yfr)
   lokph=ceq%phase_varres(lokcs)%phlink
   nsl=phlista(lokph)%noofsubl
   allocate(skip(nsl+1))
   info=1
   nmax=0
   do ii=1,nsl
      if(phlista(lokph)%nooffr(ii).eq.1) then
         nmax=nmax+1
         skip(nmax)=info
!         write(*,*)'QF skipping column/row ',info
      endif
      info=info+phlista(lokph)%nooffr(ii)
   enddo
   if(nmax.eq.nsl) then
! phase has no variable composition
      value=one
      goto 1000
   endif
   skip(nmax+1)=phlista(lokph)%tnooffr+1
!   write(*,*)'QF dimension: ',nofc,nmax
   nofc=nofc-nmax
   nofc2=nofc*(nofc+1)/2
   allocate(d2g(nofc2))
   allocate(eigenv(nofc))
   allocate(work(3*nofc))
   nsl=1
   cc=0
   rr=0
   row: do ii=1,nofc+nmax
      if(ii.eq.skip(nsl)) then
! skip this column
         nsl=nsl+1
         cycle row
      endif
      cc=cc+1
      rr=cc-1
      column: do jj=ii,nofc+nmax
         do zz=1,nmax
            if(jj.eq.skip(zz)) cycle column
         enddo
         rr=rr+1
!         write(*,17)'QF assigning ',ii,jj,' to ',cc,rr
!17       format(a,2i4,a,2i4)
         d2g(ixsym(cc,rr))=ceq%phase_varres(lokcs)%d2gval(ixsym(ii,jj),1)
      enddo column
   enddo row
!
!-------------------------------------------------------------------
! uncomment the call to dspev in order to make Q work
! AND link to LAPACK
!-------------------------------------------------------------------
! use LAPACK routine, note d2g is destroyed inside dspev
!   write(*,*)'LAPACK routine DSPEV not implemented'
   call dspev('N','U',nofc,d2g,eigenv,dummy,1,work,info)
!   info=-1000
   if(info.eq.0) then
!      write(*,120)(eigenv(ii),ii=1,nofc)
!120   format('Eigenvalues: ',6(1pe10.2))
! return the most negative value
      value=eigenv(1)
   else
      write(*,*)'Error calculating eigenvalues of phase matrix',info
      gx%bmperr=4321
   endif
1000 continue
   return
 end subroutine calc_qf

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calculate_reference_state(kstv,iph,ics,aref,ceq)
! Calculate the user defined reference state for extensive properties
! kstv is the typde of property: 1 U, 2 S, 3 V, 4 H, 5 A, 6 G
! It can be phase specific (iph.ne.0) or global (iph=0)
! IMPORTANT
! For integral quantitites (like calculated here) the reference state
! is ignored unless all components have the same phase as reference (like Hmix)
   implicit none
   integer kstv,iph,ics
   double precision aref
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! BIG BUG, the values of %gval is not restored!!
! kstv=1  2  3  4  5  6 other values cared for elsewhere
!      U  S  V  H  A  G
   integer iel,phref,allcomp
   double precision gref(6),bref(6),xmol(maxel),wmass(maxel),xxx(6)
   double precision tmol,tmass,bmult
!
!   write(*,*)'Reference states not implemented yet'; goto 1000
!   write(*,*)'3F reference state:',kstv,iph,ics
   if(kstv.lt.1 .or. kstv.gt.6) then
!      write(*,*)'3F No reference state for kstv: ',kstv
      goto 1000
   endif
   aref=zero
   bref=zero
   gref=zero
   xxx=zero
   allcomp=0
! loop for all components to extract the value of their reference states
! Multiply that with the overall composition (iph=0) or the phase composition
   xmol=zero
   do iel=1,noofel
! this is the reference phase for component iel
      phref=ceq%complist(iel)%phlink
!      write(*,*)'3F Reference for: ',iel,phref
! added when starting to handle P as variable.  V should not depend
! on a reference state unless all have the same phase as reference
      if(allcomp.eq.0) then
         if(phref.gt.0) then
            allcomp=phref
!            write(*,*)'3F Setting allcomp: ',allcomp
         else
! at least one component has no reference phase, ignore all refernce states
            aref=zero
            goto 900
         endif
      elseif(phref.ne.allcomp) then
! different reference phases for the components, ignore the reference state
!         writing(*,*)'3F Ignoring reference state as not same for all'
         aref=zero
         goto 900
!      else
! phref is same, continue the loop
! ignore any user defined reference state for the other components
      endif
! UNFINISHED ?? For integral properties, kstv=1..
      if(phref.gt.0) then
! we should use the phase index, not location in call below
!         write(*,*)'3F ref.ph: ',phref,phlista(phref)%alphaindex
         phref=phlista(phref)%alphaindex
! special endmember call that returns G, G.T, G.P, G.T.T, G.T.P and G.P.P
!         write(*,73)'3F R state: ',iel,phref,ceq%complist(iel)%endmember
73       format(a,2i3,2x,10i4)
!         write(*,*)'3F callcg_endmember 3: ',phref
         call calcg_endmember6(phref,ceq%complist(iel)%endmember,gref,ceq)
         if(gx%bmperr.ne.0) then
            write(*,*)'3F Error return: ',gx%bmperr
            goto 1000
         endif
         if(iph.gt.0) then
! multiply with mole fractions of phase iph,ics
            call calc_phase_molmass(iph,ics,xmol,wmass,tmol,tmass,bmult,ceq)
         else
! multiply with overall mole fractions
            call calc_molmass(xmol,wmass,tmol,tmass,ceq)
         endif
! note xxx, bref and gref are arrays
         xxx=bref+xmol(iel)*gref
!         write(*,70)'3F rs: ',bref,gref,xxx,(xmol(ij),ij=1,noofel)
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
900 continue
!   write(*,75)kstv,aref
75 format('3F ref:',i3,6(1pe12.4))
1000 continue
   return
 end subroutine calculate_reference_state

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine sortinphtup(n,m,xx)
! subroutine sortinphtup(n,m,xx,ceq)
! subroutine to sort the values in xx which are in phase and compset order
! in phase tuple order.  This is needed by the TQ interface
! The number of values belonging to the phase is m (for example composition)
! argument ceq added as new composition sets can be created ...
   integer n,m
!   double precision xx(n*m)
   double precision xx(*)
!   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
!
   integer iz,jz,kz,lz,lokph,aha
   double precision, dimension(:), allocatable :: dum
! I assume the values are NP(*), maybe there are other cases ...
! Karl had overflow error in dum ... no problem to make it a little larger
! but then I cannot set xx=dum below ...
   allocate(dum(n*m+10))
!   write(*,*)'3F corrected sortinphtup',n,m
!   write(*,10)'3F in: ',(xx(iz),iz=1,n*m)
10 format(a,10(f7.4))
   kz=0
   do iz=1,noofph
      lokph=phases(iz)
      do jz=1,phlista(lokph)%noofcs
!         if(jz.gt.1) then
! in xx the values are sequentially for all composition sets for this phase
! But they should be stored in tuple order and compset 2 etc comes at the end
! the index to the tuple is in %phtups
! phlista(lokph)%linktocs(jz) is index of phase_varres record for compset
! firsteq%phase_varres(..)%phtupx is index of phase tuple for compset
! There can be m values (for example compositions) for each phase
! BUG FIXED: Sigli example gives hard error here
! index '0' of array 'firsteq' below lower boundary of 1
         aha=(firsteq%phase_varres(phlista(lokph)%linktocs(jz))%phtupx-1)*m
!         if(aha.ne.kz) then
!            write(*,*)'3F shifting from, to, values: ',kz,aha,m
!         endif
         do lz=1,m
            dum(aha+lz)=xx(kz+lz)
         enddo
         kz=kz+m
      enddo
   enddo
!   xx=dum
   do iz=1,n*m
      xx(iz)=dum(iz)
   enddo
   deallocate(dum)
!   write(*,10)'3F ut: ',(xx(iz),iz=1,n*m)
1000 continue
   return
 end subroutine sortinphtup

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_svfun(cline,last,ceq)
! enter a state variable function
   implicit none
   integer last
   character cline*(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer, parameter :: npfs=20
   integer ks,maxsym,ipos,jt,js,kdot,nsymb,allowch
   character name2*16,pfsym(npfs)*60,string*128,pfsymdenom*60
!   integer istv(npfs),indstv(4,npfs),iref(npfs),iunit(npfs),lokv(npfs)
   integer iarr(10,npfs),lokv(npfs)
! memory leak
   type(gtp_state_variable), target :: svrvar
   type(gtp_state_variable), pointer :: svr
   type(putfun_node), pointer :: lrot,datanod
!    
! maxsym is negative to allow the user to enter abs(maxs) symbols
! pfsym are the entered symbols
! lokv is only internal strage in putfun
! lrot is the root node of expression
! nsymb is the number of user entered symbols
!    write(kou,17)'enter svgun ',last,cline(1:20),nsvfun
17 format(a,i3,2x,a,i3)
   call gparc('Symbol name: ',cline,last,ichar('='),name2,' ',q1help)
   call capson(name2)
!   write(*,*)'3F enter_svfun: ',last,name2,':',cline(1:10)
   if(.not.proper_symbol_name(name2,0)) goto 1000
   do ks=1,nsvfun
      if(name2.eq.svflista(ks)%name) then
         gx%bmperr=4136; goto 1000
      endif
   enddo
! added allowch to handle symbols including & and #
   allowch=1
! TO BE IMPLEMENTED: enter symbols with dummy arguments like CP(@P1)=HM(@P1).T
! where @Pi is a phase, @Ci is a component and @Si is a species
! these dummy variables must be defined in symbol name ?? why ?? maybe not
   call gparc('Expression, end with ";" :',cline,last,6,string,';',q1help)
   maxsym=-npfs
   ipos=1
   call putfun(string,ipos,maxsym,pfsym,lokv,lrot,allowch,nsymb)
   if(pfnerr.ne.0) then
      pfnerr=0; gx%bmperr=4134; goto 1000
   endif
! on return nsymb is the number of external symbols used in the function
! these can be other functions or state variables or used defined identifiers
! like Curie temperature etc.  The symbols are in pfsym(1..nsymb)
!
!   write(*,11)nsymb,(pfsym(js)(1:len_trim(pfsym(js))),js=1,nsymb)
11 format('3F args: ',i2,': ',10(1x,a,','))
! identify symbols as state variables, if derivative there is a dot
   iarr=0
   jt=0
   svr=>svrvar
   do js=1,nsymb
      kdot=index(pfsym(js),'.')
      if(kdot.gt.0) then
! derivatives must be stored as two state variables
!         write(*,*)'3F Found dot derivative: ',kdot,pfsym(js)
! Only allow a single symbol in this case!!!
         if(nsymb.gt.1) then
!            write(*,*)'3F Only a single symbol allowed!'
            gx%bmperr=4320; goto 1000
         endif
         jt=1
! denominator, variable after . for with the derivative is taken
         pfsymdenom=pfsym(js)(kdot+1:)
         pfsym(js)(kdot:)=' '
         call decode_state_variable(pfsym(js),svr,ceq)
         if(gx%bmperr.ne.0) goto 1000
! store in the old way in iarr for two state variables
         iarr(1,js)=svr%oldstv
         iarr(2,js)=svr%norm
         iarr(3,js)=svr%unit
         iarr(4,js)=svr%phref
         iarr(5,js)=svr%argtyp
         iarr(6,js)=svr%phase
         iarr(7,js)=svr%compset
         iarr(8,js)=svr%component
         iarr(9,js)=svr%constituent
         iarr(10,js)=jt
         call decode_state_variable(pfsymdenom,svr,ceq)
         if(gx%bmperr.ne.0) goto 1000
! store in the old way in iarr for two state variables
         iarr(1,js+1)=svr%oldstv
         iarr(2,js+1)=svr%norm
         iarr(3,js+1)=svr%unit
         iarr(4,js+1)=svr%phref
         iarr(5,js+1)=svr%argtyp
         iarr(6,js+1)=svr%phase
         iarr(7,js+1)=svr%compset
         iarr(8,js+1)=svr%component
         iarr(9,js+1)=svr%constituent
      else
         call decode_state_variable(pfsym(js),svr,ceq)
         if(gx%bmperr.ne.0) then
! symbol not a state variable, may be another function
!            write(*,*)'3F not state variable: ',gx%bmperr,' "',&
!                 pfsym(js)(1:len_trim(pfsym(js))),'"'
            gx%bmperr=0
            do ks=1,nsvfun
               if(pfsym(js).eq.svflista(ks)%name) then
!                  write(*,*)'3F found another function: ',trim(pfsym(js))
                  iarr(1,js)=-ks
                  goto 390
               endif
            enddo
            write(*,*)'3F not a function: "',&
                 pfsym(js)(1:len_trim(pfsym(js))),'"'
            gx%bmperr=4135; goto 1000
         else
! This can be a parameter identifier like mobility: mq&fe(fcc)
!            write(*,*)'3F decoded 1: ',trim(pfsym(js)),svr%oldstv
!            write(*,*)'3F decoded 2: ',svr%statevarid
! to avoid confusing this with another function index subtract 1000
! Store in the old way in iarr
!            iarr(1,js)=svr%oldstv-1000
            if(svr%oldstv.lt.0) then
               iarr(1,js)=svr%oldstv-1000
            else
               iarr(1,js)=svr%oldstv
            endif
            iarr(2,js)=svr%norm
            iarr(3,js)=svr%unit
            iarr(4,js)=svr%phref
            iarr(5,js)=svr%argtyp
            iarr(6,js)=svr%phase
            iarr(7,js)=svr%compset
            iarr(8,js)=svr%component
            iarr(9,js)=svr%constituent
         endif
      endif
390   continue
   enddo
! for derivatives two iarr arrays
! Found bug in store_putfun if just a variable entered, coefficient set to 0.0
   call store_putfun(name2,lrot,nsymb+jt,iarr)
   if(nsymb.eq.0) then
! this is just a constant numeric value ... store it locally
      if(.not.associated(lrot%left) .and. .not.associated(lrot%left)) then
!         write(*,*)'3F just a constant!!'
! set bit to allow change the value but do not allow R to be changed
         if(nsvfun.gt.3) then
            svflista(nsvfun)%status=ibset(svflista(nsvfun)%status,SVCONST)
         endif
      endif
   endif
1000 continue
! possible memory leak
   nullify(svr)
   return
 end subroutine enter_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine set_putfun_constant(svfix,value)
! changes the value of a putfun constant
! svfix is index, value is new value
   implicit none
   integer svfix
   double precision value
!\end{verbatim} %+
   type(putfun_node), pointer :: lrot
   if(.not.btest(svflista(svfix)%status,SVCONST)) then
      write(*,*)'Symbol is not a constant'
      gx%bmperr=4323
   else
      lrot=>svflista(svfix)%linkpnode
      write(*,*)'3F current and new: ',lrot%value,value
      lrot%value=value
   endif
1000 continue
   return
 end subroutine set_putfun_constant

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine store_putfun(name,lrot,nsymb,iarr)
! enter an expression of state variables with name name with address lrot
! nsymb is number of formal arguments
! iarr identifies these
! idot if derivative
   implicit none
   character name*(*)
   type(putfun_node), pointer :: lrot
   integer nsymb
   integer iarr(10,*)
!\end{verbatim} %+
   integer jf,jg,idot
!   write(*,*)'3F: store_putfun ',nsvfun
   nsvfun=nsvfun+1
   if(nsymb.gt.0) then
      allocate(svflista(nsvfun)%formal_arguments(10,nsymb))
      idot=10
! dot derivatives have two consequtive symbols for the variable before/after
      do jf=1,nsymb
! the order is: 1: state variable (negative means index to another symbol)
! 2-5: norm, unit, phref, argtyp, 
! 6-10: phase, compset, component, constituent, derivative
         do jg=1,idot
            svflista(nsvfun)%formal_arguments(jg,jf)=iarr(jg,jf)
         enddo
!         write(*,77)(iarr(jg,jf),jg=1,idot)
77       format('3F: store_putfun: ',20i3)
      enddo
   endif
   svflista(nsvfun)%name=name
   svflista(nsvfun)%linkpnode=>lrot
   svflista(nsvfun)%status=0
   svflista(nsvfun)%narg=nsymb
! this is the number of actual argument needed (like @P, @C and @S)
   svflista(nsvfun)%nactarg=0
! eqnoval indicate which equilibrium to use to get its value.
! default is 0 meaning current equilibria, can be changed by AMEND SYMBOL
   svflista(nsvfun)%eqnoval=0
1000 continue
   return
 end subroutine store_putfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine store_putfun_old(name,lrot,nsymb,&
       istv,indstv,iref,iunit,idot)
! enter an expression of state variables
! name: character, name of state variable function
! lrot: pointer, to a putfun_node that is the root of the stored expression
! nsymb: integer, number of formal arguments
! istv: integer array, formal argument state variables typ
! indstv: 2D integer array, indices for the formal state variables
! iref: integer array, reference for the formal state variables
! iunit: integer array, unit of the formal state variables
   implicit none
   type(putfun_node), pointer :: lrot
   integer nsymb
   integer, dimension(*) :: istv,iref,iunit,idot
   integer, dimension(4,*) :: indstv
   character name*(*)
!\end{verbatim}
   integer jf
!    write(*,*)'3F store_putfun ',nsvfun
   nsvfun=nsvfun+1
   if(nsymb.gt.0) then
      allocate(svflista(nsvfun)%formal_arguments(10,nsymb))
      do jf=1,nsymb
         svflista(nsvfun)%formal_arguments(1,jf)=istv(jf)
         svflista(nsvfun)%formal_arguments(2,jf)=indstv(1,jf)
         svflista(nsvfun)%formal_arguments(3,jf)=indstv(2,jf)
         svflista(nsvfun)%formal_arguments(4,jf)=indstv(3,jf)
         svflista(nsvfun)%formal_arguments(5,jf)=indstv(4,jf)
         svflista(nsvfun)%formal_arguments(6,jf)=iref(jf)
         svflista(nsvfun)%formal_arguments(7,jf)=iunit(jf)
         svflista(nsvfun)%formal_arguments(8,jf)=idot(jf)
      enddo
   endif
   svflista(nsvfun)%name=name
   svflista(nsvfun)%linkpnode=>lrot
   svflista(nsvfun)%status=0
   svflista(nsvfun)%narg=nsymb
! this is the number of actual argument needed (like @P, @C and @S)
   svflista(nsvfun)%nactarg=0
! eqnoval indicate which equilibrium to use to get its value.
! default is 0 meaning current equilibria, can be changed by AMEND SYMBOL
   svflista(nsvfun)%eqnoval=0
1000 continue
   return
 end subroutine store_putfun_old

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine find_svfun(name,lrot,ceq)
! finds a state variable function called name (no abbreviations)
   implicit none
   character name*(*)
   integer lrot
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
! name must be in UPPER CASE and exact match required
   do lrot=1,nsvfun
!      write(*,*)'3F find_svfun: ',name,svflista(lrot)%name,lrot
      if(name.eq.svflista(lrot)%name) goto 500
   enddo
   write(*,*)'3F No such state variable: ',name
   gx%bmperr=4188; goto 1000
!
500 continue
! nothing more to do!
1000 continue
   return
 end subroutine find_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine list_svfun(text,ipos,lrot,ceq)
! list a state variable function
   implicit none
   character text*(*)
   integer ipos,lrot
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! copied svflista(lrot)%formal_arguments(2..5,jt) to indices as gfortran error
!   integer indstv(4)
   type(gtp_state_variable), target :: svr2
   type(gtp_state_variable), pointer :: svr
   character symbols(20)*32,afterdot*32
   integer js,jt,ip,istv,kl,mm
!    write(*,*)'3F list_svfun 1:',svflista(lrot)%narg
   if(lrot.le.0 .or. lrot.gt.nsvfun) then
      gx%bmperr=4140; goto 1000
   endif
   if(svflista(lrot)%narg.eq.0) goto 500
   js=0
   jt=0
100 continue
      jt=jt+1
      js=js+1
      ip=1
      symbols(js)=' '
      istv=svflista(lrot)%formal_arguments(1,jt)
!      write(*,*)'3F list_svfun: ',istv,js,jt
      if(istv.gt.-1000 .and. istv.lt.0) then
! istv<-1000 means this is a parameter property identifier
! function refer to another function (assuming never to have 1000 symbols ...
         symbols(js)=svflista(-istv)%name
      else
! the 1:10 was a new bug discovered in GNU fortran 4.7 and later
         svr=>svr2
         call make_stvrec(svr,svflista(lrot)%formal_arguments(1:10,jt))
         if(gx%bmperr.ne.0) then
            write(*,*)'3F failed creating state variable record'
            goto 1000
         endif
         call encode_state_variable(symbols(js),ip,svr,ceq)
         if(gx%bmperr.ne.0) then
            write(*,*)'3F failed encode state variable'
            goto 1000
         endif
         if(svflista(lrot)%formal_arguments(10,jt).ne.0) then
! a derivative!!!
!            write(*,111)'3F A dot derivative of ',js,jt,symbols(js)
111         format(a,2i3,': ',a)
            jt=jt+1
            afterdot=' '
            ip=1
            svr=>svr2
            call make_stvrec(svr,svflista(lrot)%formal_arguments(1:10,jt))
            call encode_state_variable(afterdot,ip,svr,ceq)
!            write(*,111)'3F wrt state variable  ',js,jt,afterdot
!            symbols(js)=symbols(js)(1:len_trim(symbols(js)))//'.'//afterdot
            symbols(js)=trim(symbols(js))//'.'//afterdot
!            write(*,111)'3F alltogether ',js,jt,symbols(js)
         endif
      endif
      if(jt.lt.svflista(lrot)%narg) goto 100
500 continue
   kl=len_trim(svflista(lrot)%name)
   text(ipos:ipos+kl+1)=svflista(lrot)%name(1:kl)//'= '
   ipos=ipos+kl+2
!   write(*,502)'3F wrtfun: ',jt,(trim(symbols(mm)),mm=1,jt)
502 format(a,i3,10(' "',a,'", '))
   call wrtfun(text,ipos,svflista(lrot)%linkpnode,symbols)
! where is pfnerr defined??
   if(pfnerr.ne.0) then
      write(kou,*)'Putfun error listing funtion ',pfnerr
      gx%bmperr=4142; goto 1000
   endif
!   text(ipos:ipos)=';'
!   ipos=ipos+1
1000 continue
   return
 end subroutine list_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine make_stvrec(svr,iarr)
! stores appropriate values from a formal argument list to a state variable
! function in a state variable record
   implicit none
   type(gtp_state_variable), pointer :: svr
   integer iarr(10)
!\end{verbatim}
   integer jt,norm,currid
!
! memory leak
!   allocate(svr)
   if(iarr(1).lt.-1000) then
! Handling of parameter property symbols like TC, BMAGN etc
! NOTE inside symbols  -1000 used to separate from other symbols 
      currid=iarr(1)+1000
      svr%statevarid=currid
   elseif(iarr(1).le.0) then
      write(*,*)'3F illegal argument to make_stvrec: ',iarr(1)
   elseif(iarr(1).lt.10) then
! This is T, P, MU, AC, LNAC
!         1  2  3   4   5
      svr%statevarid=iarr(1)
      currid=iarr(1)
   else
! This is U, S, V, H, A,  G,  NP, BP, DG, Q,   N,  X,  B,  W,  Y   symbol
!         6  7  8  9  10, 11, 12, 13, 14, 15,  16, 17, 18, 19. 20   new code
!         10 20 30 40 50  60  70  80  90  100  110 111 120 122 130  old code
! dvs iarr()=10 means U etc.
      jt=iarr(1)/10+5
      norm=mod(iarr(1),10)
! special for x and w, note norm is set to normallizing
      if(jt.eq.16 .and. norm.eq.1) jt=17
      if(jt.eq.18 .and. norm.eq.2) jt=19
      svr%statevarid=jt
      currid=iarr(1)
!      write(*,*)'3F make: ',iarr(1),jt
   endif
!   write(*,11)iarr
11 format('3F Arguments: ',10i5)
!   svr%oldstv=iarr(1)
   svr%oldstv=currid
   svr%norm=iarr(2)
   svr%unit=iarr(3)
   svr%phref=iarr(4)
   svr%argtyp=iarr(5)
   svr%phase=iarr(6)
   svr%compset=iarr(7)
   svr%component=iarr(8)
   svr%constituent=iarr(9)
1000 continue
   return
 end subroutine make_stvrec

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine list_all_svfun(kou,ceq)
! list all state variable funtions
   implicit none
   integer kou
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character text*256
   integer ks,ipos
   write(kou,17)
17 format('List of all state variable symbols'/' No Name = expression ;')
   do ks=1,nsvfun
      ipos=1
      call list_svfun(text,ipos,ks,ceq)
      if(pfnerr.ne.0) then
         gx%bmperr=4142; pfnerr=0; goto 1000
      endif
      write(kou,76)ks,text(1:ipos-1)
76    format(i3,2x,a)
   enddo
1000 continue
   return
 end subroutine list_all_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine evaluate_all_svfun_old(kou,ceq)
! THIS SUBROUTINE MOVED TO MINIMIZER but kept for initiallizing
! cannot be used for state variable functions that are derivatives ...
! evaluate and list values of all functions but it is still used somewhere
   implicit none
   integer kou
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   character actual_arg(10)*24
   integer kf
   double precision val
   if(kou.gt.0) write(kou,75)
75 format('No  Name ',12x,'Value')
   do kf=1,nsvfun
! actual arguments needed if svflista(kf)%nactarg>0
      val=evaluate_svfun_old(kf,actual_arg,0,ceq)
      if(gx%bmperr.ne.0) goto 1000
      if(kou.gt.0) write(kou,77)kf,svflista(kf)%name,val
77    format(i3,1x,a,1x,1PE15.8)
   enddo
1000 continue
   return
 end subroutine evaluate_all_svfun_old

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 double precision function evaluate_svfun_old(lrot,actual_arg,mode,ceq)
! THIS SUBROUTINE MOVED TO MINIMIZER
! but needed in some cases in this module ... ???
! envaluate all funtions as they may depend on each other
! actual_arg are names of phases, components or species as @Pi, @Ci and @Si
! needed in some deferred formal parameters  (NOT IMPLEMENTED YET)
   implicit none
   integer lrot,mode
   character actual_arg(*)*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   double precision argval(20)
   type(gtp_state_variable), target :: svr2
   type(gtp_state_variable), pointer :: svr
   integer jv,jt,istv,ieq
   double precision value
   argval=zero
   value=zero
!    write(*,*)'3F evaluate_svfun ',lrot,svflista(lrot)%narg,svflista(lrot)%name
! locate function
   if(lrot.le.0 .or. lrot.gt.nsvfun) then
      gx%bmperr=4140; goto 1000
   endif
   if(svflista(lrot)%narg.eq.0) goto 300
! get values of arguments
   jv=0
   jt=0
100 continue
      jt=jt+1
      istv=svflista(lrot)%formal_arguments(1,jt)
      if(istv.lt.0) then
! if eqnoval nonzero it indicates from which equilibrium to get its value
         ieq=svflista(lrot)%eqnoval
         if(ieq.eq.0) then
            value=ceq%svfunres(-istv)
         else
            value=eqlista(ieq)%svfunres(-istv)
         endif
!          write(*,*)'3F evaluate_svfun symbol',ieq,value
      else
! the 1:10 was a new bug discovered in GNU fortran 4.7 and later
         svr=>svr2
         call make_stvrec(svr,svflista(lrot)%formal_arguments(1:10,jt))
         if(gx%bmperr.ne.0) goto 1000
         if(svflista(lrot)%formal_arguments(10,jt).eq.0) then
! get state variable value
            call state_variable_val(svr,value,ceq)
         else
! state variable derivative, error code here should be handelled by calling
! routine and use meq_evaluate_evaluate
!            write(*,*)'3F In evaluate_svfun_old!!!'
!            write(*,*)'Use "calculate symbol" for state variable derivatives!'
            gx%bmperr=4217
!            call make_stvrec(svr2,svflista(lrot)%formal_arguments(1:10,jt))
!            call state_var_value_derivative(svr,svr2,value,ceq)
!            call meq_state_var_value_derivative(svr,svr2,value,ceq)
         endif
         if(gx%bmperr.ne.0) goto 1000
      endif
      jv=jv+1
      argval(jv)=value
      if(jt.lt.svflista(lrot)%narg) goto 100
! all arguments evaluated (or no arguments needed)
300 continue
!    write(*,333)'evaluate_svfun ',svflista(lrot)%name,argval(1),argval(2)
!333 format(a,a,2(1PE15.6))
!   write(*,340)'evaluate svfun 1: ',mode,lrot
340 format(a,5i4)
   modeval: if(mode.eq.0 .and. btest(svflista(lrot)%status,SVFVAL)) then
! If mode=0 and SVFVAL set return the stored value
      value=ceq%svfunres(lrot)
!      write(*,350)'evaluate svfun 2: ',0,lrot,value
   elseif(mode.eq.0 .and. btest(svflista(lrot)%status,SVFEXT)) then
! if mode=0 and SVFEXT set use value from equilibrium eqno
      ieq=svflista(lrot)%eqnoval
      if(ceq%eqno.eq.ieq) then
         value=evalf(svflista(lrot)%linkpnode,argval)
         if(pfnerr.ne.0) then
            write(*,*)'evaluate_svfun putfunerror ',pfnerr
            gx%bmperr=4141; goto 1000
         endif
         ceq%svfunres(lrot)=value
!         write(*,350)'evaluate svfun 3: ',ieq,lrot,value
      else
         value=eqlista(ieq)%svfunres(lrot)
      endif
!      write(*,350)'evaluate svfun 4: ',ieq,lrot,value
350 format(a,2i3,1pe12.4)
   else
! if mode=1 always evaluate
      value=evalf(svflista(lrot)%linkpnode,argval)
      if(pfnerr.ne.0) then
         write(*,*)'evaluate_svfun putfunerror ',pfnerr
         gx%bmperr=4141; goto 1000
      endif
   endif modeval
! save value in current equilibrium
   if(lrot.gt.0) ceq%svfunres(lrot)=value
1000 continue
!   write(*,*)'3F eval_svfun: ',lrot,value,size(ceq%svfunres)
   evaluate_svfun_old=value
   return
 end function evaluate_svfun_old

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
