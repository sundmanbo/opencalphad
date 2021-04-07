!
! GTP3F included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     10. State variable manipulations
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_stable_state_var_value
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
! looking for bug ... not here
!   write(*,*)'3F calling get_state_var_value'
   call get_state_var_value(modstatevar,value,encoded,ceq)
!   write(*,*)'3F back from get_state_var_value',value,' ',trim(encoded)
1000 continue
! possible memory leak
!   write(*,*)'3F exit get_stable_state_var_value'
   nullify(svr)
   return
 end subroutine get_stable_state_var_value

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_state_var_value
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
!check if there is a "." (dot) neaing it is a dot derivative
   if(index(statevar,'.').gt.0) then
      write(*,*)'3F dot derivatives must be entered as symbols: ',trim(statevar)
      gx%bmperr=4399; goto 1000
   endif
   call decode_state_variable(statevar,svr,ceq)
!   write(*,20)statevar(1:len_trim(statevar)),svr%oldstv,svr%norm,&
!        svr%argtyp,svr%component,gx%bmperr
20  format('3F gsvv 1: ',a,' : ',5i3)
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
!      call find_svfun(name,lrot,ceq)
      call find_svfun(name,lrot)
      if(gx%bmperr.ne.0) then
         write(*,*)'3F Neither state variable or symbol, maybe model-param-id?'
         gx%bmperr=4399; goto 1000
      else
! get the value of the symbol, may involve other symbols and state variablse
! The actual_arg is a facility not yet implemented and not allowed here
! if mode=0 the stored value may be used, mode=1 always evaluate
!         write(*,*)'3F Found function: ',lrot
         actual_arg=' '
         mode=1
         if(btest(svflista(lrot)%status,SVFDOT)) then
            gx%bmperr=4399; goto 1000
         endif
! this is OK if it is not a derivative
! BUT be careful!! it can be a value that must be calculated explicitly!!
         if(btest(svflista(lrot)%status,SVFVAL)) then
            value=ceq%svfunres(lrot)
!            write(*,*)'3F Extracting saved value for: ',trim(name),value
         else
!            write(*,*)'3F call svaluate_svfun_old 1'
            value=evaluate_svfun_old(lrot,actual_arg,mode,ceq)
            if(gx%bmperr.eq.4217) goto 1000
         endif
         encoded=name
      endif
   else
! it is a real state variable
!      write(*,*)'3F calling state_variable_val from get_state_var_value'
      call state_variable_val(svr,value,ceq)
!      write(*,*)'3F back from state_variable_val',value
      if(gx%bmperr.ne.0) goto 1000
      ip=1
      encoded=' '
      call encode_state_variable(encoded,ip,svr,ceq)
      if(gx%bmperr.ne.0) then
         write(*,*)'3F encode error: ',trim(encoded),gx%bmperr
         gx%bmperr=0; encoded='dummy'
      endif
!      write(*,*)'3F get_state_var_value encoded: ',trim(encoded)
   endif
1000 continue
! possible memory leak
   nullify(svr)
   return
 end subroutine get_state_var_value

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_many_svar
!\begin{verbatim}
 subroutine get_many_svar(statevar,values,mjj,kjj,encoded,ceq)
! called with a state variable name with wildcards allowed like NP(*), X(*,CR)
! mjj is dimension of values, kjj is number of values returned
! encoded used to specify if phase data in phasetuple order ('Z')
! >>>> BIG question: How to do with phases that are note stable?
! If I ask for w(*,Cr) I only want the fraction in stable phases
! but when this is used for GNUPLOT the values are written in a matix
! and the same column in that phase must be the same phase ...
! so I have to have the same number of phases from each equilibria.
! tentative added feature: # instead of * means also metastable phases
! BEWHERE # is used also for composition set and sublattice index!
!
! CURRENTLY if x(*,*) and x(*,A) mole fractions only in stable phases
!
! >>>>>>>>>>>>>>>> there is a segmentation fault in this subroutine when
! called from ocplot2 in the map11.OCM
! for the second plot as part of all.OCM
! but not when called by itself.  SUCK
! probably caused by the fact that the number of composition sets are different
! >>>>>>>>>>>>>>>>
! A new segmentation fault for map2 when plotting with 2 maptops and the
! first does not have a new composition set LIQUID_AUTO#2 created in the 
! second map.  I do not understand how that has ever worked??
! >>>>>>>>>>>>>>>>
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   character statevar*(*),encoded*(*)
   double precision values(*)
   integer mjj,kjj
!\end{verbatim}
   integer indices(4),modind(4)
   double precision xnan,xxx
   integer jj,lokph,lokcs,k1,k2,k3,iref,jl,iunit,istv,enpos,maxen
   logical onlystable
! memory leak
   type(gtp_state_variable), target :: svrvar
   type(gtp_state_variable), pointer :: svr
!   logical phtupord
! check for character overflow, leave at least 100 at end
   maxen=len(encoded)-30
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
   svr=>svrvar
   call decode_state_variable(statevar,svr,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'3F Failed decode statevar in get_many_svar',gx%bmperr
      goto 1000
   endif
!   write(*,*)'3F get_many_svar 1: ',trim(statevar),svr%argtyp,svr%phase
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
!   write(*,*)'3F get_many_svar 2: ',trim(statevar),svr%argtyp,indices
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
! indices(1)=-10 phase & compset means all phases also metastable "#"
   jj=0
   onlystable=.true.
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
! check for overflow in encoded
            if(enpos.gt.maxen) goto 1100
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
! check for overflow in encoded
               if(enpos.gt.maxen) goto 1100
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
! check for overflow in encoded
               if(enpos.gt.maxen) goto 1100
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
      elseif(indices(2).eq.-3 .or. indices(2).eq.-10) then
! if indices(1)>=0 then indices(2)<0 must means a loop for all phase+compset
!         write(*,*)'3F seg.fault ',noofph
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
! check for overflow in encoded
               if(enpos.gt.maxen) goto 1100
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
! check for overflow in encoded
            if(enpos.gt.maxen) goto 1100
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
! check for overflow in encoded
                  if(enpos.gt.maxen) goto 1100
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
!   elseif(indices(1).eq.-3) then
   elseif(indices(1).eq.-3 .or. indices(1).eq.-10) then
! loop for phase+compset as indices(1+2)
! here we must be careful not to destroy original indices, use modind
      if(indices(1).eq.-10) onlystable=.FALSE.
!      write(*,*)'3F get_many NP(*) etc 1: ',gx%bmperr,indices(1),indices(3),&
!           onlystable,noofph
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
! skip if phase+compset suspended
!            if(ceq%phase_varres(lokcs)%phstate.le.PHSUS)
!            if(indices(3).eq.0) then
            if(indices(3).le.0) then
! This is typically listing of NP(*) for all phases
               modind(3)=0
               call encode_state_variable3(encoded,enpos,istv,modind,&
                    iunit,iref,ceq)
               if(gx%bmperr.ne.0) then
                  write(*,*)'3F error encoding state variable'; goto 1000
               endif
               enpos=enpos+1
! check for overflow in encoded
               if(enpos.gt.maxen) goto 1100
               jj=jj+1
               if(jj.gt.mjj) goto 1100
! if the wildcard is # include also metastable
               if(onlystable .and. &
                    ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) then
                  values(jj)=xnan
               else
                  call state_variable_val3(istv,modind,iref,&
                       iunit,values(jj),ceq)
                  if(gx%bmperr.ne.0) then
                     write(*,*)'3F error calling __val3'; goto 1000
                  endif
               endif
!            elseif(ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) then
!               call encode_state_variable3(encoded,enpos,istv,modind,&
!                    iunit,iref,ceq)
!               if(gx%bmperr.ne.0) goto 1000
!               enpos=enpos+1
! check for overflow in encoded
               if(enpos.gt.maxen) goto 1100
!               values(jj)=xnan
            elseif(indices(3).gt.0) then
! This is typically listing of w(*,cr), only in stable range of phases
               modind(3)=indices(3)
!               write(*,*)'3F statevar 1A: ',modind(1),modind(2),modind(3)
               call get_phase_compset(modind(1),modind(2),lokph,lokcs)
               if(gx%bmperr.ne.0) goto 1000
               call encode_state_variable3(encoded,enpos,istv,modind,&
                    iunit,iref,ceq)
               if(gx%bmperr.ne.0) goto 1000
               enpos=enpos+1
! check for overflow in encoded
               if(enpos.gt.maxen) goto 1100
               jj=jj+1
!               write(*,*)'3F statevar 1B: ',trim(encoded),jj,&
!                    lokph,ceq%phase_varres(lokcs)%phlink
!               write(*,*)'3F statevar 1B: ',jj,&
!                    lokph,ceq%phase_varres(lokcs)%phlink
               if(jj.gt.mjj) goto 1100
!--------------------------------------------------------------
! Beware of segmentation faults at next call to state_variable_val3 !!!
! This code should take care of the problem when new composition sets
! have been created in different step/map commands and we then try
! to extract wildcard state variables from these to plot.
! ceq here can be a previous local ceq used for the step/map 
! without a compset created later.
!               write(*,333)'3F this composition set may not exist',&
!                    lokcs,ceq%phase_varres(lokcs)%phstate,PHENTSTAB,&
!                    ceq%phase_varres(lokcs)%phlink,lokph
333            format(a,i4,2x,2i3,2x,2i4)
               if(.not.allocated(ceq%phase_varres(lokcs)%yfr)) then
! if yfr is not allocated the composition set does not exist, skip this phase
!                  write(*,*)'3F this composition set does not exist',&
!                       ceq%phase_varres(lokcs)%phstate,PHENTSTAB
                  values(jj)=xnan; goto 600
               endif
!--------------------------------------------------------------
               if(ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) then
! if phase is not stable (phstate= -2, -1 or 0)set dummy value
                  values(jj)=xnan
               else
                  call state_variable_val3(istv,modind,iref,&
                       iunit,values(jj),ceq)
               endif
!               write(*,*)'3F statevar 1C: ',jj,values(jj)
!               write(*,*)'3F statevar 1C: ',trim(encoded),jj,values(jj)
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
! check for overflow in encoded
                  if(enpos.gt.maxen) goto 1100
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
! check for overflow in encoded
                  if(enpos.gt.maxen) goto 1100
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
600         continue
            if(gx%bmperr.ne.0) then
               write(*,19)'3F error 3',modind,gx%bmperr
19             format(a,4i4,i7)
               goto 1000
            endif
         enddo csloop
      enddo phloop
!      write(*,*)'3F jj: ',jj
   else
! error if here
      write(*,17)'3F Illegal set of indices 5',(indices(jl),jl=1,4)
      gx%bmperr=4317; goto 1000
   endif
1000 continue
! possible memory leak, BUT nullify does not release memory
   nullify(svr)
   kjj=jj
   return
1100 continue
   write(*,1102)enpos,maxen,jj,mjj
1102 format('3F Overflow using get_many_svar: ',2i6,5x,2i6)
   gx%bmperr=4317; goto 1000
 end subroutine get_many_svar

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine decode_state_variable
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

!\addtotable subroutine decode_state_variable3
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
!        9      10     11     12  ...
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
   istv=0
!   write(*,*)'3F in decode3 "',trim(statevar),'" ',istv
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
! -10 for also metastable phases and composition sets, using hash #
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
!   write(*,*)'3F decode_state_var 1: "',trim(lstate),'"'
! compare first character
   ch1=lstate(1:1)
!   write(*,*)'3F decoding: ',trim(statevar),is,' ',ch1
   do is=1,noos
      if(ch1.eq.svid(is)(1:1)) goto 50
   enddo
! it may be a property, parameter identifier
   if(deblist) write(*,*)'3F jump to 600!'
   goto 600
!------------------------------------------------------------
50 continue
! There is an ambiguous case with first letter: "AC" and "A"
! If is=4 it means we found an A, check if second letter is "C"
   if(is.eq.4) then
      if(lstate(2:2).ne.'C') then
         is=10
      endif
   endif
!   write(*,*)'3F we are here statevar ',trim(statevar),', ch1: ',ch1
   if(deblist) write(*,*)'3F dsv 1: ',ch1,is
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
   if(kp.lt.jp) then
!      write(*,*)'3F cannot find )',trim(lstate),jp,kp
      goto 1140
   endif
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
!         write(*,*)'3F findconst for ref: ',arg2,cmass,icc
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
! All this is for extensive properties 
! If we come here the first (and sometimes second) letter must have been:
!               A,  B, BP,  D,  G, H,  N, NP,  Q, S, U,  W,  X,  Y
! and "is" is  10, 18, 13, 14, 11, 9, 16, 12, 15, 7, 6, 19, 17, 20
! NOTE: for N and B the second character has been checked and jp incremented
!       if equal to P.  The third (for NP and BP forth) character must 
!       be normallizing (MWVF), a space or a (, otherwise it is a property
   if(deblist) write(*,*)'3F lstate: ',lstate(1:20),jp,is
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
! UNLESS MIXED REFERENCE STATES FOR THE ELEMENTS
   if(lstate(jp:jp).eq.'S') then
      jp=jp+1
      iref=-1
   elseif(lstate(jp:jp).eq.'R') then
      write(*,*)'Ignoring suffix "R" on ',trim(statevar),&
           ', user reference is default'
      jp=jp+1
   endif
   if(btest(ceq%status,EQMIXED)) then
! user has different phases as reference state for the elements use SER
!      write(*,*)'3F Mixed reference state for the elements, SER used'
! This is only set for integral U(6), S(7), H(9), A(10), G(11) (not V(8))
      if(is.ge.6 .and. is.le.11) iref=-1
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
         write(*,110)'3F dsv 5: ',is,jp,kp,lstate(1:20)
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
         elseif(arg1(1:2).eq.'# ') then
! include also values for metastable phases, only for single argument variables
            iph=-10
            ics=-10
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
!            write(*,*)'3F findconst 2: ',trim(arg2),cmass,icc
            if(gx%bmperr.ne.0) goto 1000
         endif
         svr%component=0
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
!      write(*,'(a,5i5)')'3F Y: ',istv,icc,indices(1),indices(2),svr%constituent
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
! typty>100 means a model-parameter-id with associated component such as MQ&FE
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

!\addtotable subroutine calc_phase_molmass
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
! bug here when calculating MAP11 because we create new composition sets
! when we map the different lines
   if(ocv()) write(*,14)'3F cpm: ',iph,ics,lokph,lokcs
14 format(a,10i5)
   allsubl: do ll=1,phlista(lokph)%noofsubl
! an error here in MAP11.OCM as the number of composition sets varied
! in the different map commands
      if(.not.allocated(ceq%phase_varres(lokcs)%sites)) then
         if(ics.eq.2) then
! fix to allow list and amend lines from mapping when a phase
! may have added composition sets ...
            call get_phase_compset(iph,1,lokph,lokcs)
         else
            write(*,*)'phase ',trim(phlista(lokph)%name),&
                 ' has no composition set ',ics
         endif
         amount=zero; gx%bmperr=4072; goto 700
      endif
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
! This is the current number of formula unit of the phase, 
! It is zero if not stable
   amount=ceq%phase_varres(lokcs)%amfu
! ceq%phase_varres(lokcs)%abnorm(1) is moles atoms for one formula unit
! ceq%phase_varres(lokcs)%abnorm(2) is mass for one formula unit
700 continue
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
!---------------------------------------------
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

!\addtotable subroutine calc_phase_mol
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

!\addtotable subroutine calc_molmass
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

!\addtotable subroutine sumprops
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
! segmentation fault here ?? during plotting using several STEP/MAP
! when new comp.sets may be allocated
! skip composition sets with no allocated yfr ....
            if(.not.allocated(ceq%phase_varres(lokcs)%yfr)) then
!               write(*,*)'3F skipping unallocated comp.set.',lokcs
               cycle allcs
            endif
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

!\addtotable subroutine encode_state_variable
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

!\addtotable subroutine encode_state_variable3
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
! problem ... component names can be different in different equilibria ....
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
         call findsublattice(indices(1),indices(3),sublat)
         if(gx%bmperr.ne.0) goto 1000
!         call get_phase_constituent_name(indices(1),indices(2),&
!              stsymb(jp:),sublat)
! I am not sure if indices(2) is constituent numbered for each sublattice
! or numbered from the beginning, assume the latter !!
!         call get_constituent_name(indices(1),indices(2),&
!              stsymb(jp:),mass)
! modified 190710/BoS as constituent index is in 3
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
   if(iref.lt.0) then
! we can have reference states for G H etc.
      stsymb(jp:jp)='S'
      jp=jp+1
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
         call findsublattice(indices(1),indices(3),sublat)
         if(gx%bmperr.ne.0) goto 1000
!         call get_phase_constituent_name(indices(1),indices(3),&
!              stsymb(jp:),sublat)
         call get_constituent_name(indices(1),indices(3),&
              stsymb(jp:),mass)
         if(gx%bmperr.ne.0) goto 1000
!         write(*,*)'3F encode y:',indices,sublat
! sublattice is the last argument
         jp=len_trim(stsymb)+1
         if(sublat.gt.1) then
            stsymb(jp:)='#'//char(ichar('0')+sublat)//')'
            jp=jp+3
         else
            stsymb(jp:jp)=')'
            jp=jp+1
         endif
!         write(*,*)'3F encode: "',trim(stsymb),'"'
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

!\addtotable subroutine findsublattice
!\begin{verbatim}
 subroutine findsublattice(iph,constix,sublat)
! find sublattice of constituent constix in phase lokph
! is lokph index to gtp_phaserecord or gtp_phase_varres??
! constix is constituent index
   implicit none
   integer iph,constix,sublat
!\end{verbatim}
   integer ix,lokph,nn
   if(iph.gt.0 .and. iph.le.noofph) then
      lokph=phases(iph)
   else
      gx%bmperr=4050; goto 1000
   endif
!   write(*,*)'3F args: ',iph,lokph,constix
!   write(*,*)'3F phase: ',phlista(lokph)%name
   if(constix.le.0) then
      write(*,*)'3F no such constituent in this phase',constix
      gx%bmperr=4399; goto 1000
   endif
!   nn=1
! BUG!! found 21.03.18 after 5 years !!
   nn=0
   loop: do sublat=1,phlista(lokph)%noofsubl
      nn=nn+phlista(lokph)%nooffr(sublat)
      if(constix.le.nn) exit loop
   enddo loop
   if(constix.gt.nn) then
      write(*,*)'3F no such constituent in this phase',constix
      gx%bmperr=4399
   endif
1000 continue
   return
 end subroutine findsublattice

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine encode_state_variable_record
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
! always use composition set 1 and assume sublattice 1 ??
         ics=1
         sublat=1
         call get_phase_name(indices(1),ics,stsymb(jp:))
         if(gx%bmperr.ne.0) goto 1000
         jp=len_trim(stsymb)+1
         stsymb(jp:jp)=','; jp=jp+1
         call findsublattice(indices(1),indices(3),sublat)
         if(gx%bmperr.ne.0) goto 1000
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
         call findsublattice(indices(1),indices(3),sublat)
         if(gx%bmperr.ne.0) goto 1000
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

!\addtotable subroutine state_variable_val
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
!   write(*,*)'3F calling state_variable_val3: '
   call state_variable_val3(istv,indices,iref,iunit,value,ceq)
   if(gx%bmperr.ne.0) then
!      write(*,920)'3F error 7: ',gx%bmperr,istv,svr%oldstv,svr%argtyp
!920   format(a,i5,2x,2i15,i2)
!   else
!      write(*,*)'3F value: ',value
   endif
!   write(*,*)'3F back from state_variable_val3: ',value
1000 continue
   return
 end subroutine state_variable_val

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine state_variable_val3
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
   double precision props(5),xmol(maxel),wmass(maxel),stoi(10)
   double precision, allocatable :: cmpstoi(:)
   double precision vt,vp,amult,vg,vs,vv,div,aref,vn,bmult,tmass,tmol
   double precision qsp,gref,spmass,rmult,tsave,rtn,spextra(10)
   integer kstv,norm,lokph,lokcs,icx,jp,ncmp,ic,iprop,loksp,nspel,iq,nspx
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
      write(*,10)'3F not impl 1: ',istv,indices,iref,iunit,gx%bmperr,value
      goto 1100
   else
! wrong or state variable not implemented
!      write(*,10)'3F not impl 2: ',istv,indices,iref,iunit,gx%bmperr,value
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
!      write(*,*)'3F stv B: ',norm,kstv,indices(1),vt,vp,amult
      if(indices(1).eq.0) then
! global value for the whole system
         vg=props(1)
         vs=-props(2)
         vv=props(3)
! normalizing: 0 none,1=M (moles), 2=W (mass), 3=W(volume), 4=F(formula unit)
!         write(*,*)'3F norm: ',norm,props(1)
         if(norm.eq.1) then
! props(4) is total number of moldes
            div=props(4)
         elseif(norm.eq.2) then
! props(5) is total mass
            div=props(5)
         elseif(norm.eq.3) then
! Normalizing per volume, there frequently no volume data
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
! normalizing: 0 none,1=M (moles), 2=W (mass), 3=W(volume), 4=F(formula unit)
! I have to think more here should normalizing be per phase or for total?
! GM(BCC) is for one mole BCC, M refer per mole of phase even if not stable
! NPM(BCC) is more fraction of BCC relative to total amount in system
!         write(*,*)'3F norm 2: ',norm,props(1),div
         if(norm.eq.1) then
! phase property normalized per phase like HM or SM
            div=ceq%phase_varres(lokcs)%abnorm(1)
! phase property normalized for whole system as NPM
            if(kstv.eq.7) div=props(4)
            rmult=div
         elseif(norm.eq.2) then
! abnorm(2) should be the mass per formulat unit
            div=ceq%phase_varres(lokcs)%abnorm(2)
! phase property normalized for whole system as NPM
            if(kstv.eq.8) div=props(5)
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
!      write(*,53)'3F more:',0,amult,vg,vp,vv,amult*(vg-vp*vv-aref)/div
53    format(a,i3,7(1PE11.3))
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
!         write(*,53)'3F more:',0,amult,vg,vp,vv,value
      elseif(kstv.eq.6) then
! 6: G
!         write(*,177)'3F G:',vg,aref,amult,div
         value=amult*(vg-aref)/div
      elseif(kstv.eq.7) then
! 7: NP
!         write(*,*)'3F npm:',norm,div
! div is normalizing, can be 1.0 or total volume
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
! 10: Q (stability, thermodynamic factor)
!         gx%bmperr=4081; goto 1000
         call calc_qf(lokcs,value,ceq)
!      else
!         write(*,*)'3F svval after 10:',kstv
      endif kstv1
!      write(*,53)'3F more:',-1,amult,vg,vp,vv,value
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
! here with state variable <0, system and user defined properties
200   continue
!   write(*,*)'3F svv3 at 200:',kstv,ndefprop
!   if(ndefprop.ne.33) then
! THIS IS A VERY CRUDE CHECK! Please check also the SELECT below !!!
! it may need to be modified !!!
!   if(ndefprop.ne.31) then modified to 32 to include VS
   if(ndefprop.ne.33) then
      write(*,*)'3F The model parameter identifiers has been changed!',33
      write(*,*)'3F You must correct state_variable_val3 in GTP3F.F90!'
! you may also have to change the case indices!!
      stop
   endif
!   write(*,*)'3F kstv: ',kstv
   select case(kstv)
   case default
      write(kou,*)'Unknown parameter identifier: ',kstv
! I need to separate out mpi's that have constituent index ...
! updated 2019.12.14
!-------------------------------------------------------------------
! These are model_parameter_ident in June 2018:
!   1 G     T P                                   0 Energy
!   2 TC    - P                                   2 Combined Curie/Neel T
!   3 BMAG  - -                                   1 Average Bohr magneton numb
!   4 CTA   - P                                   2 Curie temperature
!   5 NTA   - P                                   2 Neel temperature
!   6 IBM   - P &<constituent#sublattice>;       12 Individual Bohr magneton num
!   7 THET  - P                                   2 Debye or Einstein temp
!   8 V0    - -                                   1 Volume at T0, P0
!   9 VA    T -                                   4 Thermal expansion
!  10 VB    T P                                   0 Bulk modulus
!  11 VC    T P                                   0 Extra volume parameter
!  12 VS    - -                                   1 Diffusion volume 
!  13 MQ    T P &<constituent#sublattice>;       10 Mobility activation energy
!  14 MF    T P &<constituent#sublattice>;       10 RT*ln(mobility freq.fact.)
!  15 MG    T P &<constituent#sublattice>;       10 Magnetic mobility factor
!  16 G2    T P                                   0 Liquid two state parameter
!  17 THT2  - P                                   2 Smooth slope function T
!  18 DCP2  - P                                   2 Smooth slope funtion step
!  19 LPX   T P                                   0 Lattice param X axis
!  20 LPY   T P                                   0 Lattice param Y axis
!  21 LPZ   T P                                   0 Lattice param Z axis
!  22 LPTH  T P                                   0 Lattice angle TH
!  23 EC11  T P                                   0 Elastic const C11
!  24 EC12  T P                                   0 Elastic const C12
!  25 EC44  T P                                   0 Elastic const C44
!  26 UQT   T P &<constituent#sublattice>;        0 UNIQUAC residual parameter
!  27 RHO   T P                                   0 Electric resistivity
!  28 VISC  T P                                   0 Viscosity
!  29 LAMB  T P                                   0 Thermal conductivity
!  30 HMVA  T P                                   0 Enthalpy of vacancy form.
!  31 TSCH  - P                                   2 Schottky anomality T
!  32 CSCH  - P                                   2 Schottky anomality Cp/R.
!  33 QCZ   - -                                   1 MQMQA coordination factor
! I am not sure how to handle changes ...
!-------------------------------------------------------------------
!...................................... without constituent index
!   case(1:5,7:12,16:25,27:31)  OLD
!   case(1:5,7:12,16:25,27:32) 
   case(1:5,7:12,16:25,27:33) 
! not with constituent index: 6: individual Bohr magneton number
! not with constituent index: 13:15: Mobilities
! not with constituent index: 26: UNIQUAC 
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
!....................................... with constituent index
! These have a constituent index
   case(6,13:15,26)
! 6: IBM& Individual Bohr magneton number
! 13-15: MQ& etc mobility values
! 26: UNIQUAC parameter tau
!      write(*,*)'3F svv3 mob1: ',indices(1),indices(2)
      call get_phase_compset(indices(1),indices(2),lokph,lokcs)
      if(gx%bmperr.ne.0) goto 1000
! property is kstv*100+indices(3) (constituent index)
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
      write(*,*)'3F model parameter value has not been calculated'
      gx%bmperr=4361
   end select
!.......................................
! all legal case values goto somewhere else
!   gx%bmperr=4113; goto 1000
   goto 1000
!-----------------------------------------------------------------
! chemical potentials, activites etc, istv is 3, 4 or 5 for MU, AC and LNAC
! there can be a reference state
500 continue
!   ceq%rtn=globaldata%rgas*ceq%tpval(1)
! if one argument that is a component, if two these are phase and constituent
! here indices(2) is considered to specify a reference state ...???
!   write(*,502)'3F refstate 500: ',iref,indices(1),indices(3)
502 format(a,10i4)
!   if(indices(2).ne.0) then
! species index is in indices(3) !!!!
   if(indices(3).ne.0) then
! This has nothing to do with reference state ... ??? see else link for that
! I wonder if this code is ever used ...
!      write(*,502)'3F species: ',iref,indices(1),indices(3)
      lokph=phases(indices(1))
      loksp=phlista(lokph)%constitlist(indices(3))
! split the species in elements, convert to components, add chemical potentials
      call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp,nspx,spextra)
      if(gx%bmperr.ne.0) goto 1000
      if(qsp.ne.zero) then
!         write(*,*)'3F Cannot calculate potential of charged species'
         gx%bmperr=4159; goto 1000
      endif
      allocate(cmpstoi(noofel))
      cmpstoi=zero
! get_species_data gives only elements with non-zero stoiciometry
      do ic=1,nspel
         cmpstoi(ielno(ic))=stoi(ic)
      enddo
!      write(*,507)'3F elstoi:',loksp,nspel,indices(3),(cmpstoi(ic),ic=1,noofel)
507   format(a,3i3,10(1pe12.4))
! elements2components1 is in gtp3G
! ncmp returned as number of elements, cmpstoi is stoichiometry of ALL elements
! stoi is no longer used ...
      call elements2components1(nspel,stoi,ncmp,cmpstoi,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,508)'3F el2comp:',loksp,nspel,(cmpstoi(ic),ic=1,noofel)
508   format(a,2i3,10(1pe12.4))
      value=zero
      do ic=1,ncmp
         value=value+cmpstoi(ic)*ceq%complist(ic)%chempot(1)
      enddo
! >>>> subtract reference state: i.e. calculate G for the phase with 
! just this constituent.  Note indices(1) is phase record, change to index
!      write(*,*)'3F refphase: ',indices(1),value
      ic=phlista(indices(1))%alphaindex
! set endmember=0 to allow vacancies ...
! HM, here I think it can only be a single species .... 190710/BoS
      endmember=0
! changed from indices(2) which is composition set number
      endmember(1)=indices(3)
! This routine returns G for current number of atoms
      call calcg_endmemberx(indices(1),endmember,gref,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,'(a,i3,2x,10i3)')'3F callcg_endmember 1: ',indices(1),endmember
      value=value-gref*ceq%rtn
!      write(*,511)'3F refstate: ',indices(1),indices(3),gref*ceq%rtn,value
511   format(a,2i3,6(1pe12.4))
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
!   write(*,53)'3F more:',-1,amult,vg,vp,vv,value
   if(allocated(cmpstoi)) deallocate(cmpstoi)
   return
1100 continue
   gx%bmperr=4078
!   write(*,*)'3F State variable value not implemented yet'
   goto 1000
 end subroutine state_variable_val3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_qf
!\begin{verbatim}
 subroutine calc_qf(lokcs,value,ceq)
! calculates eigenvalues of the second derivative matrix, stability function
! using the Darken matrix with second derivatives: OK FOR SUBSTITUTIONAL
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
!   integer lokph,nsl
!   lokph=ceq%phase_varres(lokcs)%phlink
!   nsl=phlista(lokph)%noofsubl
!   write(*,*)'3F calc_qf: ',lokph,lokcs,nsl
!   if(nsl.eq.1) then
! For substitutional solutions
!      write(*,*)'3F nsl 1: ',nsl
!      call calc_qf_old(lokcs,value,ceq)
!   else
! For any onther model      
!      write(*,*)'3F nsl 2: ',nsl
      call calc_qf_romain(lokcs,value,ceq)
!   endif
 end subroutine calc_qf

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 subroutine calc_qf_romain(lokcs,value,ceq)
! calculates eigenvalues of the second derivative matrix, stability function
! using Romain Le Tellier proposal eliminating one dependent fraction in
! each sublattice and also one ion if charge balance.  Also ignore sublattices
! with a single constituent
! lokcs is index of phase_varres
! value calculated value returned
! ceq is current equilibrium
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer :: lokcs
   double precision value
!\end{verbatim}
! Algorithm: use the Hessian modified by eliminating one constituent "k_n"
! in each sublattice and one ionic constituent (if any)
! The terms of this "redused" Hessaian will be
! \sum_i\ne k_n \sum_j\ne k_n d2G/dy_idy_j - 
!        \sum_k_n( d2G/dy_idy_kn+ d2G/dy_jdy_kn - \sum_k_m d2G/dy_kndy_km)
! 
   integer lokph,ii,jj,nsl,ncol,jrow,lrow,subrow,jcol,lcol,subcol,tnfr,nn,info
! skip the last constituent in each sublattice, 
! skip sublattices with a single constituent
! skip one charged constituent ... NOT IMPLEMENTED
! The reduced Hessian (symmetric)
   double precision, allocatable, dimension(:) :: hessian
! work and result arrays. mm is needed when external charge balance
   double precision, allocatable, dimension(:) :: work,eigenval,ionfact,mm
! this is needed as argument but never used
   double precision eigenvect(1)
   type(gtp_phase_varres), pointer :: varres
   logical debug,excb
!
!   debug=.true.
   debug=.false.
   varres=>ceq%phase_varres(lokcs)
   lokph=varres%phlink
   nsl=phlista(lokph)%noofsubl
   tnfr=phlista(lokph)%tnooffr
   ncol=tnfr-nsl
   if(ncol.eq.0) then
! fixed composition
      value=one; goto 1000
   endif
   if(btest(phlista(lokph)%status1,PHEXCB)) then
! external charge balance, we need mm
!      allocate(mm(phlista(lokph)%tnooffr)
!      excb=.TRUE.
      write(*,*)'Stability check of charged phases not implemented'
      value=1.0D2
      goto 1000
   else
      excb=.FALSE.
   endif
!   write(*,*)'3F allocate: ',ncol*(ncol+1)/2,tnfr,nsl
   allocate(hessian(ncol*(ncol+1)/2))
! loop for all rows
   ii=0
   subrow=1
   lrow=phlista(lokph)%nooffr(subrow)
   jrow=0
   nn=0
   row: do while(ii.lt.tnfr)
      ii=ii+1
      if(ii.eq.lrow) then
         subrow=subrow+1
         if(subrow.gt.nsl) exit row
         lrow=lrow+phlista(lokph)%nooffr(subrow)
         if(phlista(lokph)%nooffr(subrow).eq.1) then
            ii=ii+1; subrow=subrow+1
            if(subrow.gt.nsl) exit row
            lrow=lrow+phlista(lokph)%nooffr(subrow)
         endif
         cycle row
      endif
! loop for all columns
      jcol=jrow
      jrow=jrow+1
      jj=ii-1
      subcol=subrow
      lcol=lrow
      col: do while(jj.lt.tnfr)
         jj=jj+1
         if(jj.eq.lcol) then
            subcol=subcol+1
            if(subcol.gt.nsl) exit col
            lcol=lcol+phlista(lokph)%nooffr(subcol)
            if(phlista(lokph)%nooffr(subcol).eq.1) then
               jj=jj+1
               subcol=subcol+1
               if(subcol.gt.nsl) exit col
               lcol=lcol+phlista(lokph)%nooffr(subcol)
            endif
            cycle col
         endif
         nn=nn+1
         jcol=jcol+1
!         write(*,'(a,3(2i4,2x))')'3F hessian: ',jrow,jcol,ii,jj,lrow,lcol
         hessian(ixsym(jrow,jcol))=varres%d2gval(ixsym(ii,jj),1)-&
              varres%d2gval(ixsym(ii,lrow),1)-&
              varres%d2gval(ixsym(jj,lcol),1)+&
              varres%d2gval(ixsym(lrow,lcol),1)
      enddo col
   enddo row
   if(debug) then
      do ii=1,ncol
         write(*,'("3F H: ",5(1pe12.4))')(hessian(ixsym(jj,ii)),jj=1,ncol)
      enddo
   endif
   if(ncol.eq.1) then
      value=hessian(1)
      goto 1000
   endif
! use LAPACK routine, note Hessian is destroyed inside dspev
   allocate(eigenval(ncol))
! work is work array at least 2*ncol, info is return code
   allocate(work(2*ncol))
   info=0
! 'N' means only eigenvalues, 'U' means Hessian is upper triangle
! ncol is dimension of Hessian, eigenval is calculated, 
! dummy values for eigenvect, 1
   call dspev('N','U',ncol,hessian,eigenval,eigenvect,1,work,info)
   if(info.eq.0) then
      if(debug) write(*,120)(eigenval(ii),ii=1,ncol)
120   format('Eigenvalues: ',6(1pe10.2))
! return the first value, negative if unstable
      value=eigenval(1)
   else
! <0 the "info" argument has illegal value
! >0 "info" off-diagonal elements if intermediate tridiagonal did not converge 
      value=zero
      write(*,*)'Error calculating eigenvalues of phase matrix',info
      gx%bmperr=4321
   endif
!   
1000 continue
! basically we are only interested if value is <0 or >0 ...
   if(value.gt.1.0D2) then
      value=1.0D2
   elseif(value.lt.-1.0D2) then
      value=-1.0D2
   endif
   return
 end subroutine calc_qf_romain

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 subroutine calc_qf_romain_old(lokcs,value,ceq)
! calculates eigenvalues of the second derivative matrix, stability function
! using Romain Le Tellier proposal eliminating one dependent fraction in
! each sublattice and also one ion if charge balance.  Also ignore sublattices
! with a single constituent
! lokcs is index of phase_varres
! value calculated value returned
! ceq is current equilibrium
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer :: lokcs
   double precision value
!\end{verbatim}
! Algorithm: use the Hessian modified by eliminating one constituent "k_n"
! in each sublattice and one ionic constituent (if any)
! The terms of this "redused" Hessaian will be
! \sum_i\ne k_n \sum_j\ne k_n d2G/dy_idy_j - 
!        \sum_k_n( d2G/dy_idy_kn+ d2G/dy_jdy_kn - \sum_k_m d2G/dy_kndy_km)
! 
   integer ncol,nsl,nion,ii,jj,ll,lokph,iy,nf,kk,info,iskip,jskip,tnfr
! kskip needed when external charge balance as one may have two skipped
! constituents in same sublattice
   integer icol,jrow,kskip
! constituents that should be skipped, one per sublattice and charge, max 10
   integer skipped(10)
! The reduced Hessian (symmetric)
   double precision, allocatable, dimension(:) :: hessian
! work and result arrays. mm is needed when external charge balance
   double precision, allocatable, dimension(:) :: work,eigenval,ionfact,mm
! this is needed as argument but never used
   double precision eigenvect(1)
   type(gtp_phase_varres), pointer :: varres
   double precision xxx,backstop,yyy
   logical excb,debug
   character phname*24
!
   write(*,*)'3F entering calc_gf_romain',lokcs
   value=zero
   varres=>ceq%phase_varres(lokcs)
   lokph=phasetuple(varres%phtupx)%ixphase
   if(btest(phlista(lokph)%status1,PHEXCB)) then
! external charge balance, we need mm
!      allocate(mm(phlista(lokph)%tnooffr)
!      excb=.TRUE.
      write(*,*)'Stability check of charged phases not implemented'
      value=1.0D2
      goto 1000
   else
      excb=.FALSE.
   endif
! Step 1: eliminate one constituent per sublattice (highest fraction)
! if external charge balance we need to eliminate one ion (the first)
! NOTE ionic liquid does not have external charge balance and d2G/dyidyj
! include the variation of sites(?)
   skipped=0
   nsl=phlista(lokph)%noofsubl
   iy=1
   kskip=1
   do ll=1,nsl
!      write(*,'(a,3i3,2x,10i3)')'3F debug ',lokph,nsl,ll,skipped
      if(phlista(lokph)%nooffr(ll).gt.1) then
         xxx=zero
         do nf=1,phlista(lokph)%nooffr(ll)
            if(varres%yfr(iy).gt.xxx) then
! skip constituent with largest fraction in each sublattice
! this gives sometimes a cusp in the middle for binaries ... not nice
               xxx=varres%yfr(iy)
               if(excb) then
                  skipped(kskip)=iy
               else
                  skipped(ll)=iy
               endif
            endif
            iy=iy+1
         enddo
      else
! totally skip sublattices with a sngle constituent
! negative value means second derivative wrt this fraction totally skipped
         skipped(ll)=-iy
         iy=iy+1
      endif
   enddo
! dimension of the Hessian. One more should be added if ionic (not now!)
   tnfr=phlista(lokph)%tnooffr
   ncol=tnfr-nsl
   if(ncol.gt.1) then
      call get_phase_name(lokph,1,phname)
      write(*,*)'3F Q for ',trim(phname),' with ',ncol,' Hessian'
      debug=.true.
   else
      debug=.false.
   endif
   allocate(hessian((ncol*(ncol+1))/2))
!   write(*,*)'3F allocated Hessian',lokcs,ncol*(ncol+1)/2
   nf=nsl
!   write(*,'(a,10i4)')'3F skipped: ',skipped
! in all terms of the Hessian we have to add the sum of all pairs of
! fractions that are skipped
! PROBABLY WRONG, WE MUST HAVE SEVERAL BACKSTOPS DEPENDING ON SUBLATTICES ...
   backstop=zero
! if we have ionic constituents nf is larger than nsl one more ...
   back1: do ii=1,nf
      if(skipped(ii).lt.0) cycle back1
      back2: do jj=ii,nf
         if(skipped(jj).lt.0) cycle back2
! this is sum_k1 \sum_k2 d2G/dy_k1dy_k2
         backstop=backstop+varres%d2gval(ixsym(skipped(ii),skipped(jj)),1)
!         write(*,*)'3F backstop 1: ',skipped(ii),skipped(jj),backstop
      enddo back2
   enddo back1
!   write(*,*)'3F backstop 2: ',0,0,backstop
!   do ii=1,tnfr
!      write(*,'(a,5(1pe12.4))')'3F d2G/dyidyj: ',&
!           (varres%d2gval(ixsym(ii,jj),1),jj=1,phlista(lokph)%tnooffr)
!   enddo
! The terms of this "reduced" Hessaian will be
! \sum_i\ne k_n \sum_j\ne k_n d2G/dy_idy_j - 
!        \sum_k_n( d2G/dy_idy_kn+ d2G/dy_jdy_kn - \sum_k_m d2G/dy_kndy_km)
   iskip=1
   jskip=1
   icol=1
   jrow=1
   loop1: do ii=1,tnfr
      if(ii.eq.abs(skipped(iskip))) then
!         write(*,*)'3F Skipping ii: ',ii,skipped(ii)
         iskip=iskip+1
         jskip=iskip
         cycle loop1
      endif
      xxx=zero
!      write(*,'(a,2i4,2x,10i3)')'3F bug 1: ',ii,nf,skipped
      do kk=1,nf
! PROBABLY WRONG, INCLUDE ONLY THOSE IN CURRENT SUBLATTICE ...
         if(skipped(kk).gt.0) then
            xxx=xxx+varres%d2gval(ixsym(ii,skipped(kk)),1)
         endif
      enddo
!      write(*,*)'3F Calculate xxx for ii: ',ii,xxx
! here the big loop to calculate the Hessian
      loop2: do jj=ii,tnfr
         if(jj.eq.abs(skipped(jskip))) then
! PROBABLY WRONG, INCLUDE ONLY THOSE IN CURRENT SUBLATTICE ...
!            write(*,*)'3F Skipping jj: ',jj,skipped(jj)
            jskip=jskip+1
            cycle loop2
         endif
         yyy=xxx
!         write(*,'(a,2i4,2x,10i3)')'3F bug 2: ',ii,nf,skipped
         do kk=1,nf
            if(skipped(kk).gt.0) then
               yyy=yyy+varres%d2gval(ixsym(jj,skipped(kk)),1)
            endif
         enddo
!         write(*,*)'3F calculate term',ii,jj
!         write(*,*)'3F Calculate yyy for ii and jj: ',ii,jj,yyy
! kxysm is indexing a symmetrix matrix when jj >= ii
! ixsym is indexing a symmetrix matrix whatever values of ii and jj
! multiply with the fractions to avoid extrapolation to infinity 
! at low fractions  ... meaningless
         xxx=varres%d2gval(ixsym(ii,jj),1)-yyy+backstop
!         if(xxx.gt.1.0d2) xxx=1.0d2
!         write(*,*)'3F index Hessian: ',icol,jrow,ixsym(icol,jrow)
         hessian(ixsym(icol,jrow))=xxx
!              varres%yfr(ii)*varres%yfr(jj)*(varres%d2gval(ixsym(ii,jj),1)-&
!              yyy+backstop)
!         write(*,'(a,4i3,4(1pe12.4))')'3F hessian ',icol,jrow,ii,jj,&
!              hessian(ixsym(icol,jrow)),varres%d2gval(ixsym(ii,jj),1),&
!              -yyy,backstop
! limit the terms in the Hession ...
         icol=icol+1
      enddo loop2
      jrow=jrow+1
      icol=jrow
   enddo loop1
! if ncol > 1 calculate eigenvalues
   if(ncol.eq.1) then
      value=hessian(1)
      goto 1000
   endif
   if(debug) then
      do ii=1,ncol
         write(*,'("3F Hessian: ",5(1pe12.4))')(hessian(ixsym(jj,ii)),jj=1,ncol)
      enddo
   endif
! use LAPACK routine, note Hessian is destroyed inside dspev
   allocate(eigenval(ncol))
! work is work array at least 2*ncol, info is return code
   allocate(work(2*ncol))
   info=0
! 'N' means only eigenvalues, 'U' means Hessian is upper triangle
! ncol is dimension of Hessian, eigenval is calculated, 
! dummy values for eigenvect, 1
   call dspev('N','U',ncol,hessian,eigenval,eigenvect,1,work,info)
   if(info.eq.0) then
      if(debug) write(*,120)(eigenval(ii),ii=1,ncol)
120   format('Eigenvalues: ',6(1pe10.2))
! return the first value, negative if unstable
      value=eigenval(1)
   else
      value=zero
      write(*,*)'Error calculating eigenvalues of phase matrix',info
      gx%bmperr=4321
   endif
!   
1000 continue
! basically we are only interested if value is <0 or >0 ...
   if(value.gt.1.0D2) then
      value=1.0D2
   elseif(value.lt.-1.0D2) then
      value=-1.0D2
   endif
   return
 end subroutine calc_qf_romain_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_qf_otis
!\begin{verbatim}
 subroutine calc_qf_otis(lokcs,value,ceq)
! NOT USED -----------------------------
! calculates eigenvalues of the second derivative matrix, stability function
! using Otis reduced Hessian method.  Should work indpent of model!!
! lokcs is index of phase_varres
! value calculated value returned
! ceq is current equilibrium
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer :: lokcs
   double precision value
!\end{verbatim}
! Algorithm:
! 1. Create a Jacobian matrix with massbalance and constraint n columns, m rows
! 2. Use an QR factorization of this, select the first n-m columns as Z
! 3. Calculate F = Z^T H Z where H is all second derivatives of G wrt fractions
! 4. Calculate eigenvalues of F.  If all positive no problem!
! 
   integer ncol,mrow,lda,info,ll,nsl,lokph,ii,jj
   double precision, allocatable, dimension(:,:) :: jac,zeta,fff
   double precision, allocatable, dimension(:) :: tau,work,eigenv
   double precision dummy(1,1)
   type(gtp_phase_varres), pointer :: varres
!
   varres=>ceq%phase_varres(lokcs)
! Step 1: Jacobian: ncol columns, mrow rows
   value=zero
   allocate(jac(ncol,mrow))
   jac=zero
!
! Step 2: QR factorisation, n>m
   allocate(tau(ncol))
   allocate(work(ncol))
   lda=ncol
   info=0
   call dgeqr2(mrow,ncol,jac,lda,tau,work,info)
   if(info.ne.0) then
      write(*,*)'Error return from DGEGR2: ',info,lokph
      goto 1000
   endif
! How to extract Q?? Documentation of DGEQR2:
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!>  and tau in TAU(i).
! zeta should be a matrix with the first ncol x ncol-mrow part of Q
   allocate(zeta(ncol,ncol-mrow))
   zeta=zero
!
! Step 3: multiply Z^T * H * Z
   allocate(fff(ncol,ncol))
   fff=zero
   do ii=1,ncol
      do jj=1,ncol
         fff(ii,jj)=fff(ii,jj)+zeta(ii,jj)*varres%d2gval(ixsym(ii,jj),1)
      enddo
   enddo
   jac=zero
   do ii=1,ncol
      do jj=1,ncol
         jac(ii,jj)=jac(ii,jj)+fff(ii,jj)*zeta(ii,jj)
      enddo
   enddo
!
! Step 4: calculate eigenvalues
! use LAPACK routine, note d2g is destroyed inside dspev
!   write(*,*)'LAPACK routine DSPEV not implemented'
   allocate(eigenv(ncol))
   info=0
   call dspev('N','U',ncol,fff,eigenv,dummy,1,work,info)
   if(info.eq.0) then
!      write(*,120)(eigenv(ii),ii=1,ncol)
120   format('Eigenvalues: ',6(1pe10.2))
! return the most negative value
      value=eigenv(1)
   else
      write(*,*)'Error calculating eigenvalues of phase matrix',info
      gx%bmperr=4321
   endif
!   
1000 continue
   return
 end subroutine calc_qf_otis

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_qf_sub
!\begin{verbatim}
 subroutine calc_qf_sub(lokcs,value,ceq)
! NOT USED -----------------------------
! calculates eigenvalues of the second derivative matrix, stability function
! using the Darken matrix with second derivatives: OK FOR SUBSTITUTIONAL
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
   integer info,nofc,nofc2,nmax,ii,jj,lokph,nsl,cc,rr,zz,pp,qq
   double precision, allocatable :: d2g(:),eigenv(:),work(:)
   double precision dummy(1,1),dmuidyii
   integer, allocatable :: skip(:)
   type(gtp_phase_varres), pointer :: varres
! number of constituents
! ignore sublattices with single constituents ....
   varres=>ceq%phase_varres(lokcs)
   nofc=size(varres%yfr)
   lokph=varres%phlink
!   nofc=size(ceq%phase_varres(lokcs)%yfr)
!   lokph=ceq%phase_varres(lokcs)%phlink
   nsl=phlista(lokph)%noofsubl
   allocate(skip(nsl+1))
   info=1
   nmax=0
   do ii=1,nsl
      if(phlista(lokph)%nooffr(ii).eq.1) then
         nmax=nmax+1
         skip(nmax)=info
!         write(*,*)'3F QF skipping column/row ',info
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
! Calculation of matrix elements modeified 2019.11.01/BoS
! See documentation for Darken stability in minimizer documentation
! dGA/dxB = d2GM/dxAdxB -
!           \sum_C x_C (d2GM/dxAdxC+d2GM/dxBdxC + \sum_D x_D d2GM/dxCdxD)
! This may not work when sublattices ?... but difficult to calculate dG_I/dn_I
! PROBLMES WHEN USER DEFINED REFERENCE STATES !! ??
! For Cr-Fe when plot of gm(bcc) before q(bcc) then q(bcc) is rubbish
!
   row: do ii=1,nofc+nmax
      if(ii.eq.skip(nsl)) then
! skip this row.  A clumsy way to skip sublattices with a single constituent
         nsl=nsl+1
         cycle row
      endif
      cc=cc+1
      rr=cc-1
      column: do jj=ii,nofc+nmax
         do zz=1,nmax
! skip this column.  A clumsy way to skip sublattices with a single constituent
            if(jj.eq.skip(zz)) cycle column
         enddo
         rr=rr+1
!         write(*,17)'QF calc: ',ii,jj,' to ',cc,rr
17       format(a,2i4,a,2i4)
         dmuidyii=varres%d2gval(ixsym(ii,jj),1)
!         write(*,33)'3F start: ',ii,jj,0,dmuidyii
33       format(a,3i3,3(1pe12.4))
!         d2g(ixsym(cc,rr))=ceq%phase_varres(lokcs)%d2gval(ixsym(ii,jj),1)
! extra summation over all constituents (except those alone in a sublattice)
         extra1: do qq=1,nofc+nmax
            do zz=1,nmax
! skip this term.  A clumsy way to skip sublattices with a single constituent
               if(qq.eq.skip(zz)) cycle extra1
            enddo
            dmuidyii=dmuidyii-varres%yfr(qq)*(varres%d2gval(ixsym(ii,qq),1)+&
                 varres%d2gval(ixsym(jj,qq),1))
!            write(*,33)'3F minus: ',ii,jj,qq,dmuidyii,&
!                 varres%yfr(qq)*varres%d2gval(ixsym(ii,qq),1),&
!                 varres%yfr(qq)*varres%d2gval(ixsym(jj,qq),1)
            extra2: do pp=1,nofc+nmax
               do zz=1,nmax
! skip this term.  A clumsy way to skip sublattices with a single constituent
                  if(pp.eq.skip(zz)) cycle extra2
               enddo
               dmuidyii=dmuidyii+varres%yfr(qq)*varres%yfr(pp)*&
                    varres%d2gval(ixsym(pp,qq),1)
!               write(*,33)'3F adding: ',0,pp,qq,dmuidyii,&
!                    varres%yfr(qq)*varres%yfr(pp)*varres%d2gval(ixsym(pp,qq),1)
            enddo extra2
         enddo extra1
!         write(*,33)'3F result: ',cc,rr,0,dmuidyii
         d2g(ixsym(cc,rr))=dmuidyii
      enddo column
   enddo row
!   do ii=1,nofc
!      write(*,21)'3F d2Gdy2: ',(varres%d2gval(ixsym(ii,jj),1),jj=1,nofc)
!   enddo
!   do ii=1,nofc
!      write(*,21)'3F dmudy: ',(d2g(ixsym(ii,jj)),jj=1,nofc)
!   enddo
21 format(a,6(1pe12.4))
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
120   format('3F Eigenvalues: ',6(1pe10.2))
! return the most negative value
      value=eigenv(1)
   else
      write(*,*)'Error calculating eigenvalues of phase matrix',info
      gx%bmperr=4321
   endif
1000 continue
   return
 end subroutine calc_qf_sub

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_qf_old
!\begin{verbatim}
 subroutine calc_qf_old(lokcs,value,ceq)
! NOT USED -----------------------------
! calculates eigenvalues of the second derivative matrix, stability function
! this old version that seems to work for Ag-Cu ...
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
!   write(*,*)'3F in calc_qf_old: '
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
120   format('Eigenvalues: ',6(1pe10.2))
! return the most negative value
      value=eigenv(1)
   else
      write(*,*)'Error calculating eigenvalues of phase matrix',info
      gx%bmperr=4321
   endif
1000 continue
   return
 end subroutine calc_qf_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calculate_reference_state
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
!\end{verbatim} %+
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
      if(kstv.eq.3) then
! added when starting to handle P as variable.  V should not depend
! on a reference state unless all have the same phase as reference
!         write(*,*)'3F Reference for: ',iel,phref
! removed as we should allow different reference stated for G and H
         if(allcomp.eq.0) then
            if(phref.gt.0) then
               allcomp=phref
            else
! at least one component has no reference phase, ignore all refernce states
               aref=zero
               goto 900
            endif
         elseif(phref.ne.allcomp) then
! different reference phases for the components, ignore the reference state
!            write(*,*)'3F Ignoring reference state as not same for all'
            aref=zero
            goto 900
! phref is same, continue the loop
! ignore any user defined reference state for the other components
         endif
      endif
! UNFINISHED ?? For integral properties, kstv=1..
100   continue
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

!\addtotable subroutine calculate_reference_state_old
!\begin{verbatim}
 subroutine calculate_reference_state_old(kstv,iph,ics,aref,ceq)
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
 end subroutine calculate_reference_state_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine enter_svfun
!\begin{verbatim}
 subroutine enter_svfun(cline,last,ceq)
! enter a state variable function
   implicit none
   integer last
   character cline*(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer, parameter :: npfs=20
   integer ks,maxsym,ipos,jt,js,kdot,nsymb,allowch,lbuf
   character name2*16,pfsym(npfs)*60,string*128,pfsymdenom*60,fbuff*256
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
   call gparcx('Symbol name: ',cline,last,ichar('='),name2,' ','?Enter symbol')
   call capson(name2)
   if(name2(1:1).eq.' ') then
      gx%bmperr=4137; goto 1000
   endif
!   write(*,*)'3F enter_svfun: ',last,name2,':',trim(cline)
   if(.not.proper_symbol_name(name2,0)) goto 1000
! nsvfun is a global variable giving current number of state variable functions
   do ks=1,nsvfun
      if(name2.eq.svflista(ks)%name) then
         gx%bmperr=4136; goto 1000
      endif
   enddo
   kdot=0
   lbuf=0
   fbuff=' '
! added allowch to handle symbols including & and #
   allowch=1
! TO BE IMPLEMENTED: enter symbols with dummy arguments like CP(@P1)=HM(@P1).T
! where @Pi is a phase, @Ci is a component and @Si is a species
! these dummy variables must be defined in symbol name ?? why ?? maybe not
!   write(*,*)'3F symbol: "',trim(cline),'"',last
77 continue
   call gparcx('Expression, end with ";" :',cline,last,6,string,';',&
        '?Enter symbol')
! there can be multiple lines, last end by ; or empty line
   if(index(string,';').le.0) then
      fbuff(lbuf+1:)=string
      lbuf=len_trim(fbuff)
      string=' '
      write(*,*)'3F Continue: '
      goto 77
   elseif(lbuf.gt.0) then
      string=fbuff
   endif
   if(index(string,';').eq.1) then
      write(*,*)'3F empty expression, maybe forgotten =?'
      gx%bmperr=4134; goto 1000
   endif
!   write(*,*)'3F expression: ',trim(string)
   maxsym=-npfs
   ipos=1
   call putfun(string,ipos,maxsym,pfsym,lokv,lrot,allowch,nsymb)
   if(pfnerr.ne.0 .or. .not.associated(lrot)) then
      write(*,*)'3F error in putfun: ',pfnerr,associated(lrot)
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
! NOT a derivative
         call decode_state_variable(pfsym(js),svr,ceq)
         if(gx%bmperr.ne.0) then
! symbol not a state variable, may be another function
!            write(*,*)'3F not state variable: ',gx%bmperr,' "',&
!                 pfsym(js)(1:len_trim(pfsym(js))),'"'
            do ks=1,nsvfun
               if(pfsym(js).eq.svflista(ks)%name) then
!                  write(*,*)'3F found another function: ',trim(pfsym(js))
                  iarr(1,js)=-ks
                  gx%bmperr=0
                  goto 390
               endif
            enddo
! here it can be a model parameter id such as THET(BCC) or MQ&FE(BCC)
            write(*,*)'3F argument not understood: "',&
                 pfsym(js)(1:len_trim(pfsym(js))),'"',gx%bmperr
            gx%bmperr=4135; goto 1000
         else
! It is a state variable or a model parameter identifier
! to avoid confusing this with another function index subtract 1000
!            write(*,*)'3F state variable or model parameter id: "',&
!                 pfsym(js)(1:len_trim(pfsym(js))),'"',gx%bmperr
!            write(*,'(a,10i5)')'3F svr: ',svr%oldstv,svr%norm,svr%unit,&
!                 svr%phref,svr%argtyp,svr%phase,svr%compset,svr%component,&
!                 svr%constituent
! Store in the old way in iarr
!            iarr(1,js)=svr%oldstv-1000
!            write(*,*)'3F state variable or model parameter id',svr%oldstv
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
! The call above updates the global value of nsvfun so it means the new symbol
   if(nsymb.eq.0) then
! this is just a constant numeric value ... store it locally.  Why .and. ??
      if(.not.associated(lrot%left) .and. .not.associated(lrot%left)) then
!         write(*,*)'3F just a constant!!'
! set bit to allow change the value but do not allow R to be changed
         if(nsvfun.gt.3) then
            svflista(nsvfun)%status=ibset(svflista(nsvfun)%status,SVCONST)
         endif
      endif
   endif
   if(kdot.gt.0) then
! this is a dot derivative, set bits
      svflista(nsvfun)%status=ibset(svflista(nsvfun)%status,SVFVAL)
      svflista(nsvfun)%status=ibset(svflista(nsvfun)%status,SVFDOT)
!      write(*,*)'3F setting explicit bit: ',SVFDOT
!   endif
   else
! this created a crash when entering a dot derivative, only notmal functions
! there seems to be a problem that already existing state variable functions
! are not evaluated so they give a correct value
      call evaluate_all_svfun_old(-1,ceq)
      if(gx%bmperr.ne.0) then
! ignore any errors
!         write(*,*)' Error calculating the state variable functions!',gx%bmperr
         gx%bmperr=0
      endif
   endif
! If a function is entered that cannot be calculated we get values such as NaN 
!   ceq%svfunres(nsvfun)=zero
!   write(*,*)'3F store zero in svfunres',nsvfun,ceq%svfunres(nsvfun)
1000 continue
! NOTE eqnoval should be zeroed
! NOTE svfval should be set if only calculated when explicitly referenced
! possible memory leak
   nullify(svr)
!   write(*,*)'3F exit enter_svfun'
   return
 end subroutine enter_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine set_putfun_constant
!\begin{verbatim} %-
 subroutine set_putfun_constant(svfix,value)
! changes the value of a putfun constant
! svfix is index, value is new value
! THIS CAN BE A FUNCTION WITH SVFVAL bit set, 
! in that case change it to a constant
   implicit none
   integer svfix
   double precision value
!\end{verbatim} %+
   type(putfun_node), pointer :: lrot
!   write(*,*)'We are in set_putfun_constant 1'
   if(btest(svflista(svfix)%status,SVFVAL)) then
! converting a symbol from an expression to a constant
! this means some loss of memory used for the expression
      svflista(svfix)%status=ibset(svflista(svfix)%status,SVCONST)
      svflista(svfix)%status=ibclr(svflista(svfix)%status,SVFVAL)
! set number of arguments to zero ... this will make a mess ...
!      svflista(svfix)%narg=0
! remove link to expression in in linkpnode ??
!      lrot=>svflista(svfix)%linkpnode
! do we have to delete the expression?  memory loss negligable ...
   endif
   if(.not.btest(svflista(svfix)%status,SVCONST)) then
      write(*,*)'Symbol is not a constant'
      gx%bmperr=4323
   else
      lrot=>svflista(svfix)%linkpnode
!      write(*,*)'3F constant: ',lrot%value,value
      svflista(svfix)%svfv=value
! duplicate value, I am not sure where ...
      lrot%value=value
   endif
1000 continue
   return
 end subroutine set_putfun_constant

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine store_putfun
!\begin{verbatim} %-
 subroutine store_putfun(name,lrot,nsymb,iarr)
! enter an expression of state variables with name name with address lrot
! nsymb is number of formal arguments
! iarr identifies these
   implicit none
   character name*(*)
   type(putfun_node), pointer :: lrot
   integer nsymb
   integer iarr(10,*)
!\end{verbatim} %+
! idot set if if derivative
   integer jf,jg,idot
!   write(*,*)'3F: store_putfun ',nsvfun
   nsvfun=nsvfun+1
   svflista(nsvfun)%status=0
   svflista(nsvfun)%tplink=0
   svflista(nsvfun)%eqnoval=0
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
   svflista(nsvfun)%narg=nsymb
! this is the number of actual argument needed (like @P, @C and @S)
   svflista(nsvfun)%nactarg=0
! eqnoval indicate which equilibrium to use to get its value.
! default is 0 meaning any equilibria, can be changed by AMEND SYMBOL
   svflista(nsvfun)%eqnoval=0
1000 continue
   return
 end subroutine store_putfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine store_putfun_old
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

!\addtotable subroutine find_svfun
!\begin{verbatim}
 subroutine find_svfun(name,lrot)
! finds a state variable function called name (no abbreviations)
! ceq not needed!!??
   implicit none
   character name*(*)
   integer lrot
!   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
! name must be in UPPER CASE and exact match required
   do lrot=1,nsvfun
!      write(*,*)'3F find_svfun: ',name,svflista(lrot)%name,lrot
      if(name.eq.svflista(lrot)%name) goto 500
   enddo
   write(*,*)'3F No such state variable function: ',name
   gx%bmperr=4188; goto 1000
!
500 continue
! nothing more to do!
1000 continue
   return
 end subroutine find_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine find_symbol_with_equilno
!\begin{verbatim}
 subroutine find_symbol_with_equilno(lrot,eqno)
! finds a state variable function with equilibrium index
   implicit none
   integer lrot,eqno
!\end{verbatim} %+
! skip the first 3 functions, R, RT and TC
   if(lrot.lt.0) lrot=3
   eqno=0
!   write(*,*)'3F find_sweq 1: ',lrot,nsvfun
   allfun: do while(lrot.lt.nsvfun)
      lrot=lrot+1
      if(svflista(lrot)%eqnoval.gt.0) then
         eqno=svflista(lrot)%eqnoval
! for debugging
!         write(*,*)'3F symbol ',svflista(lrot)%name,' at equilibrium ',eqno
         goto 1000
      endif
   enddo allfun
   lrot=0
1000 continue
   return
 end subroutine find_symbol_with_equilno

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine list_svfun
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
! istv<-1000 means this is a model_parameter_identifier
! function refer to another function (assuming never to have 1000 symbols ...
         symbols(js)=svflista(-istv)%name
      else
! the 1:10 was a new bug discovered in GNU fortran 4.7 and later
! PROBABLE MY BUG 2020-08-31/BOS, not declared allocatable ... SUCK
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
!         write(*,*)'3F list_svfun: ',trim(symbols(js)),js,jt
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
! add special information, first fill with blanks
   text(ipos:)=' '
   if(svflista(lrot)%eqnoval.gt.0) then
! symbol should only be evaluated in equilibrium EQNOVAL
      write(text(ipos:ipos+3),470)svflista(lrot)%eqnoval
470   format(i4)
   elseif(svflista(lrot)%tplink.gt.0) then
! symbol is imported from or exported to TP function
      write(text(ipos:ipos+3),470)svflista(lrot)%tplink
   endif
   js=ipos+4
   ipos=ipos+7
! Mark with a letter in position 5!
   if(btest(svflista(lrot)%status,SVNOAM)) then
! symbol must not be amended (for R, RT and T_C)
      text(js:js)='N'
   elseif(btest(svflista(lrot)%status,SVCONST)) then
! symbol is a constant (can be amended)
      text(js:js)='C'
   elseif(btest(svflista(lrot)%status,SVFDOT)) then
! symbol is a dot derivative calculated only when explicitly referenced
      text(js:js)='D'
   elseif(btest(svflista(lrot)%status,SVFVAL)) then
! symbol calculated only when explicitly referenced
      text(js:js)='V'
   elseif(btest(svflista(lrot)%status,SVFEXT)) then
! symbol evaluated for an equilibrium, the equilibrium number already written
      text(js:js)='X'
   elseif(btest(svflista(lrot)%status,SVEXPORT)) then
! symbol imported from TP function, function index is specified
      text(js:js)='E'
   elseif(btest(svflista(lrot)%status,SVIMPORT)) then
! symbol exported to a TP function (assess coeff), function index is specified
      text(js:js)='I'
   endif
! name and expression
!   kl=len_trim(svflista(lrot)%name)
!   text(ipos:ipos+kl+1)=svflista(lrot)%name(1:kl)//'= '
   text(ipos:)=trim(svflista(lrot)%name)//'='
   ipos=len_trim(text)+2
! svflista(lrot)%linkpode is a pointer to a pufun_node record
   if(.not.associated(svflista(lrot)%linkpnode)) then
      text(ipos:)=' = no expression; '
   else
!      write(*,502)'3F wrtfun: ',jt,(trim(symbols(mm)),mm=1,jt)
502   format(a,i3,10(' "',a,'", '))
      call wrtfun(text,ipos,svflista(lrot)%linkpnode,symbols)
! where is pfnerr defined??
      if(pfnerr.ne.0) then
         write(kou,*)'Putfun error listing funtion ',pfnerr
         gx%bmperr=4142; goto 1000
      endif
   endif
1000 continue
   return
 end subroutine list_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine make_stvrec
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
   currid=0
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

!\addtotable subroutine list_all_svfun
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
17 format('List of all state variable symbols'/' No Special Name= expression ;')
!17 format('List of all state variable symbols'/' No        Name= expression ;')
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

!\addtotable subroutine evaluate_all_svfun_old
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
!      write(*,*)'3F call svaluate_svfun_old 2'
      val=evaluate_svfun_old(kf,actual_arg,0,ceq)
      if(gx%bmperr.ne.0) goto 1000
      if(kou.gt.0) write(kou,77)kf,svflista(kf)%name,val
77    format(i3,1x,a,1x,1PE15.8)
   enddo
1000 continue
   return
 end subroutine evaluate_all_svfun_old

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable double precision function evaluate_svfun_old
!\begin{verbatim} %-
 double precision function evaluate_svfun_old(lrot,actual_arg,mode,ceq)
! THIS SUBROUTINE MOVED TO MINIMIZER
! but needed in some cases in this module ... ???
! envaluate all funtions as they may depend on each other
! actual_arg are names of phases, components or species as @Pi, @Ci and @Si
! needed in some deferred formal parameters  (NOT IMPLEMENTED YET)
! if mode=1 always evaluate
   implicit none
   integer lrot,mode
   character actual_arg(*)*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   double precision argval(20)
   type(gtp_state_variable), target :: svr2
   type(gtp_state_variable), pointer :: svr
   integer jv,jt,istv,ieq
! added to handle symbols that are model parameter id
   integer indices(4),iref,iunit
   double precision value
   argval=zero
   value=zero
! calculate symbol * does not come here for H298 ...
!   write(*,*)'3F evaluate_svfun ',lrot,svflista(lrot)%narg,svflista(lrot)%name
! locate function
   if(lrot.le.0 .or. lrot.gt.nsvfun) then
      gx%bmperr=4140; goto 1000
   endif
   if(btest(svflista(lrot)%status,SVFDOT)) then
!      write(*,*)'3F Warning has SVFDOT set, return error ',lrot
      gx%bmperr=4399; goto 1000
   elseif(btest(svflista(lrot)%status,SVFVAL) .and. mode.ne.1) then
! this symbol is keeps it value unless evaluated explicitly (mode=1)
!      write(*,*)'3F Warning has SVFVAL set: ',lrot,svflista(lrot)%name,value
      value=ceq%svfunres(lrot)
      goto 1000
   elseif(btest(svflista(lrot)%status,SVFEXT)) then
! the symbol is associated with a specific equilibrium we must fetch
! its value from that equilibrium unless that is ceq!!
      ieq=svflista(lrot)%eqnoval
!      write(*,*)'3F SVFEXT set: ',lrot,ieq,svflista(lrot)%name
      if(ieq.gt.0 .and. ieq.ne.ceq%eqno) then
         value=eqlista(ieq)%svfunres(lrot)
! save its value also in this equilibrium
         goto 900
      endif
   endif
   if(svflista(lrot)%narg.eq.0) goto 300
!--------------------------------------------------------------------
! get values of arguments ... THIS IS NOT IMPLEMENTED ... I think ??
   jv=0
   jt=0
100 continue
      jt=jt+1
      istv=svflista(lrot)%formal_arguments(1,jt)
!      write(*,*)'3F get argument: ',istv,lrot,svflista(lrot)%eqnoval
      if(istv.lt.0) then
! evidently istv<1 can also mean this is a model parameter identifier
! how to know?  Here only when entering the symbol?
         if(istv.lt.-1000) then
            ieq=-istv-1000
!            write(*,*)'3F model parameter identifier: ',ieq
!            write(*,*)'3F allocated: ',size(svflista(lrot)%formal_arguments)
!            write(*,'(a,10i5)')'3F svflista: ',&
!                 svflista(lrot)%formal_arguments(5,jt),&
!                 svflista(lrot)%formal_arguments(6,jt),&
!                 svflista(lrot)%formal_arguments(7,jt)
! VERY VERY CLUMSY, must be changes to use svr state variable record
! indices are PHASE, COMPSET, COMPONENT, 
            indices(1)=svflista(lrot)%formal_arguments(5,jt)
            indices(2)=svflista(lrot)%formal_arguments(6,jt)
            indices(3)=svflista(lrot)%formal_arguments(7,jt)
            indices(4)=0
            iref=0
            iunit=0
            call state_variable_val3(-ieq,indices,iref,iunit,value,ceq)
!            value=zero
         else
! if eqnoval nonzero it indicates from which equilibrium to get its value
            ieq=svflista(lrot)%eqnoval
            if(ieq.eq.0) then
               value=ceq%svfunres(-istv)
            else
               value=eqlista(ieq)%svfunres(-istv)
            endif
         endif
      else
! the 1:10 was a new bug discovered in GNU fortran 4.7 and later
! FOUND PROBABLE BUG 2020-08-31/BOS %formal_arguments never allocated ???
         svr=>svr2
         call make_stvrec(svr,svflista(lrot)%formal_arguments(1:10,jt))
         if(gx%bmperr.ne.0) goto 1000
         if(svflista(lrot)%formal_arguments(10,jt).eq.0) then
! get state variable value
            call state_variable_val(svr,value,ceq)
         else
! state variable derivative, error code set above, it must be handelled
!  by calling other routine and use meq_evaluate_svfun
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
! all (if any) arguments evaluated (or no arguments needed)
!--------------------------------------------------------------------
300 continue
!    write(*,333)'evaluate_svfun ',svflista(lrot)%name,argval(1),argval(2)
!333 format(a,a,2(1PE15.6))
!   write(*,340)'3F evaluate svfun 1: ',mode,lrot
340 format(a,5i4)
   modeval: if(mode.eq.0 .and. btest(svflista(lrot)%status,SVFVAL)) then
! If mode=0 and SVFVAL set return the stored value
      value=ceq%svfunres(lrot)
!      write(*,350)'3F evaluate svfun 2: ',0,lrot,value
   elseif(mode.eq.0 .and. btest(svflista(lrot)%status,SVFEXT)) then
! if mode=0 and SVFEXT set use value from equilibrium eqno
      ieq=svflista(lrot)%eqnoval
      if(ceq%eqno.eq.ieq) then
         value=evalf(svflista(lrot)%linkpnode,argval)
         if(pfnerr.ne.0) then
            write(*,*)'3F evaluate_svfun putfunerror 1',pfnerr
            gx%bmperr=4141; pfnerr=0; buperr=0; goto 1000
         endif
         ceq%svfunres(lrot)=value
!         write(*,350)'3F evaluate svfun 3: ',ieq,lrot,value
      else
! Hm, we already did this earlier ... redundant?
         value=eqlista(ieq)%svfunres(lrot)
      endif
!      write(*,350)'3F evaluate svfun 4: ',ieq,lrot,value
350 format(a,2i3,1pe12.4)
   else
! if mode=1 always evaluate unless another equilibrium, we jumped to 900 above
      value=evalf(svflista(lrot)%linkpnode,argval)
      if(pfnerr.ne.0) then
         write(*,*)'3F evaluate_svfun putfunerror 2',pfnerr
         gx%bmperr=4141; pfnerr=0; buperr=0; goto 1000
      endif
   endif modeval
!   if(btest(svflista(lrot)%status,SVFVAL)) then
!    if(lrot.gt.4) write(*,*)'3F evaluated symbol: ',lrot,value
!   endif
! save value in current equilibrium
900 continue
   if(lrot.gt.0) then
      ceq%svfunres(lrot)=value
!      if(lrot.gt.4) write(*,*)'Saved symbol ',lrot,' in equil ',ceq%eqno,value
   endif
1000 continue
!   write(*,*)'3F eval_svfun: ',lrot,value,size(ceq%svfunres)
   evaluate_svfun_old=value
   return
 end function evaluate_svfun_old

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
