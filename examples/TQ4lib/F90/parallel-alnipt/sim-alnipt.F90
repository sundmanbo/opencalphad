!
! Small program to simulate diffusion in 1D using OC
! Ternary system Al-Ni-Pt coating of superallys
!
program sim2
!
  use liboctqa
!$ use omp_lib
!
  implicit none
  integer, parameter :: gpmax=100,cmax=10,sfile=31
  integer notused,cnum(10),nsel,gpix(gpmax),gp,gp1,gp2,nt,cc1,cc2
  integer howmany,plut,modmod,maxloop,half,phtup,jj,jp,gpcur,nrow,ioerr
  character database*60,selel(cmax)*2,gpname*24,line*60
  character setupfile*60,profiles*60,phasename*24
  double precision tpval(2),nval,xval(cmax),dx(cmax),dxmax
  double precision xxx,xxy,dmu,dff,sum,gpxval(cmax)
  double precision sumpos,sumneg,sdxp,sdxn,mobi(cmax),xin(cmax)
  double precision eqcpu1,eqcpu2,eqcpusum
  double precision xleft(cmax),xright(cmax),xsum
  integer eqcc1, eqcc2,eqccsum,lrest,rrest,iel,rest
  double precision, parameter :: mumin=1.0d-4
  logical parallel
  type(gtp_equilibrium_data), pointer :: ceq,gridceq
!
  type grid
! This has a pointer to one equilibrium per gridpoint
! and has an array of mole fractions in that gridpoint for diffusion
     type(gtp_equilibrium_data), pointer :: eqp
     double precision xval(cmax)
  end type grid
  type(grid), dimension(gpmax), target :: gpp
  type(grid), pointer :: gridpoint
!
  write(*,1)
1 format('Simulation of uphill diffusion Al-Ni-Pt, OC example 2021'/&
        'Surface coating of Pt-Al on an Al-Ni turbine blade (see PDF)'//&
        'In the input file is specified on separate lines:'/&
        'database file, symbol of each element'/&
        'name of stable phase,'/&
        'number of gridpoints'/&
        'Al and Pt mole fractions in the blade (Ni rest),'/&
        'Al and Ni mole fractions of the surface layer (Pt rest),'/&
        'Mobilities (mol/(J s)) of Al, Ni and Pt,'/&
        'Maximum number of timesteps'/&
        'Minimal change of fractions'/&
        'file name for intermediate output (empty line means no output),'/&
        'empty line means run in parallel,'//&
        'Any lines starting with a hash caracter "#" are ignored.'//)
!
  write(*,2,advance='no')
2 format('Setup file name:')
  read(*,3)setupfile
3 format(a)
! open input file
  open(sfile,file=setupfile,access='sequential',status='old',iostat=ioerr)
  if(ioerr.ne.0) stop 'Error opening input file'
!
!-------------------------------------
! initiate
  call tqini(notused,ceq)
  if(gx%bmperr.ne.0) stop 'Cannot initiate TQ interface'
! read database and elements
  call readline(sfile,database)
  nsel=1
  eloop: do while(.true.)
     call readline(sfile,selel(nsel))
     if(selel(nsel)(1:1).eq.' ') exit eloop
     call capson(selel(nsel))
     nsel=nsel+1
  enddo eloop
  nsel=nsel-1
  rest=nsel
! read database with selected elements, rest is zet to zero!
  call tqrpfil(database,rest,selel,ceq)
  if(gx%bmperr.ne.0) stop 'Cannot read selected elememts from database'
!-----------------
! suspend all phases except the specified phase
  call tqphsts(-1,-2,zero,ceq)
  if(gx%bmperr.ne.0) stop 'TQ error suspending phases'
! set specified phase as ENTERED
  call readline(sfile,phasename)
  call tqgpi(phtup,phasename,ceq)
  if(gx%bmperr.ne.0) stop 'TQ error finding stable phase'
  call tqphsts(phtup,1,one,ceq)
  if(gx%bmperr.ne.0) stop 'TQ error entering stable phase'
!-----------------
! set phasename as reference state for all elements ?
!-----------------
! number of gridpoints
  call readline(sfile,line)
  jp=1
  call getint(line,jp,gpcur)
  if(buperr.ne.0) stop 'Input error 1'
  if(gpcur.lt.10 .or. gpcur.gt.100) then
     write(*,*)'Gridpoints must be between 10 and 100'
     stop 'Error in input file'
  endif
!-----------------
! Left side (interior) initial composition in alphabetical order
  call readline(sfile,line)
  jp=1
  xsum=zero
  xleft=zero
  lrest=0
  do iel=1,nsel
     call getrel(line,jp,xleft(iel))
! the metlib routines have a global error code buperr
     if(buperr.ne.0) then
        if(line(jp:jp).eq.'*') then
           lrest=iel
           buperr=0
           jp=jp+1
        else
           write(*,*)'Error reading left side fractions for element ',iel
           stop
        endif
     elseif(xleft(iel).le.zero .or. xleft(iel).ge.one) then
        write(*,'(a,i3)')'Mole fraction must be between 0 and 1 for element',iel
        stop 'Error in setup file'
     endif
     xsum=xsum+xleft(iel)
  enddo
  if(lrest.eq.0) stop 'No element defined as "rest" on left hand side'
  xleft(lrest)=one-xsum
  if(xleft(lrest).ge.one) then
     write(*,*)'Fractions should add up to unity'
     stop 'Error in setup file'
  endif
!-----------------
! Right side (surface) initial composition in alphabetical order
  call readline(sfile,line)
  jp=1
  xsum=zero
  xright=zero
  rrest=0
  do iel=1,nsel
     call getrel(line,jp,xright(iel))
     if(buperr.ne.0) then
        if(line(jp:jp).eq.'*') then
           rrest=iel
           buperr=0
           jp=jp+1
        else
           write(*,*)'Error reading right side fractions for element ',iel
           stop
        endif
     elseif(xright(iel).le.zero .or. xright(iel).ge.one) then
        write(*,'(a,i3)')'Mole fraction must be between 0 and 1 for element',iel
        stop 'Error in setup file'
     endif
     xsum=xsum+xright(iel)
  enddo
  if(rrest.eq.0) stop 'No element defines as "rest" on right hand side'
  xright(rrest)=one-xsum
  if(xright(rrest).le.zero) then
     write(*,*)'Fractions on right hand side should add up to unity'
     stop 'Error in setup file'
  endif
!----------------------------
! Temperature
  call readline(sfile,line)
  jp=1
  call getrel(line,jp,tpval(1))
  if(buperr.ne.0) stop 'Input errot 7'
  if(tpval(1).lt.2.0D2 .or. tpval(1).gt.3.0D3) then
     write(*,*)'T must be between 200 and 3000 K'
     stop 'Error in setup file'
  endif
! Pressure
  tpval(2)=1.0D5
!----------------------------
! Mobilities in alphabetical order
  call readline(sfile,line)
  jp=1
  do iel=1,nsel
     call getrel(line,jp,mobi(iel))
     if(buperr.ne.0) stop 'Input error 8'
     if(mobi(iel).gt.0.1 .or. mobi(iel).lt.1.0D-6) then
        write(*,*)'Mobility data out of range for element ',iel
        stop 'Error in setup file'
     endif
  enddo
!----------------------------
! Max timesteps
  call readline(sfile,line)
  jp=1
  call getint(line,jp,maxloop)
  if(buperr.ne.0) stop 'Input error 8'
  if(maxloop.lt.10) stop 'Too few timesteps'
! Minimal change in fractions
  call readline(sfile,line)
  jp=1
  call getrel(line,jp,dff)
  if(dff.le.zero .or. dff.gt.0.1) stop 'Fraction change out of range'
!----------------------------
! intermediate output
  call readline(sfile,profiles)
  if(profiles(1:1).ne.' ') then
     modmod=1
     plut=30
     open(plut,file=profiles,access='sequential',status='unknown',iostat=ioerr)
     if(ioerr.ne.0) stop 'Cannot open profile output file'
! start GNUPLOT graphics
     write(30,10)
10   format('# GNUPLOT file for plotting profiles'/&
          'set terminal wxt size 840,700 font "Arial,16"'/&
          'set title "OpenCalphad   simulator"'/&
          'set origin 0.0, 0.0 '/'set size   1.0,   1.0'/&
          'set xlabel "Gridpoint"'/'set ylabel "Fractions"'/&
          'set key top right font "Arial,12"')
  else
     plut=0
     modmod=maxloop
  endif
!----------------------------
! run sequentially/parallel
  call readline(sfile,line)
  if(line(1:1).ne.' ') then
     parallel=.FALSE.
  else
     parallel=.TRUE.
  endif
!------------------------------------- end of input
! Echo input on screen and output file (if any)
  write(*,17)trim(database),(selel(iel),iel=1,nsel)
  write(*,18)trim(phasename),gpcur,(xleft(iel),iel=1,nsel)
  write(*,19)(xright(iel),iel=1,nsel)
  write(*,20)(mobi(iel),iel=1,nsel)
  write(*,21)tpval(1),maxloop,dff
  if(plut.gt.0) write(*,22)trim(profiles)
17 format(/'# Input values: '/'# ',a,1x,10(1x,a))
18 format('# ',a,i5/'# ',10(F7.4))
19 format('# ',10(F7.4))
20 format('# ',10(1PD12.4))
21 format('# ',F8.2,i8,D12.4)
22 format('# ',a//)
  if(plut.gt.0) then
     write(plut,17)trim(database),(selel(iel),iel=1,nsel)
     write(plut,18)trim(phasename),gpcur,(xleft(iel),iel=1,nsel)
     write(plut,19)(xright(iel),iel=1,nsel)
     write(plut,20)(mobi(iel),iel=1,nsel)
     write(plut,21)tpval(1),maxloop,dff
! Start data for graphics 
     if(plut.gt.0) write(plut,23)
23   format(/'$profile << EOD')
  endif
  eqcpusum=zero
  eqccsum=0
!---------------------------------------
! set conditions
  nval=1.0D0
! cnum not used at all!!
  call tqsetc('T ',0,0,tpval(1),cnum(1),ceq)
  if(gx%bmperr.ne.0) goto 1000
  call tqsetc('P ',0,0,tpval(2),cnum(2),ceq)
  if(gx%bmperr.ne.0) goto 1000
  call tqsetc('N ',0,0,nval,cnum(3),ceq)
  if(gx%bmperr.ne.0) goto 1000
!
  do iel=1,nsel
     if(iel.ne.lrest) then
        call tqsetc('X ',iel,0,xleft(iel),cnum(4),ceq)
     endif
  enddo
  if(gx%bmperr.ne.0) goto 1000
!---------------------------------------
! calculate an equilibrium just to test
  write(*,*)'Calculating equilibrium with gridminimizer'
  call tqlc(kou,ceq)
  call tqce(' ',0,0,zero,ceq)
  if(gx%bmperr.ne.0) goto 1000
!---------------------------------------
! list the equilibrium results as a check
  call tqlr(kou,ceq)
  if(gx%bmperr.ne.0) goto 1000
!---------------------------------------
! create grid with gpcur equilibria with initial composition
! half the gridpoints has the left hand composition, the other the right hand
  gpname='GP_001'
  half=gpcur/2
  write(*,43)
43 format(/'Creating grid'/'Gridpoint  name   Fractions in alphabetical order')
  xval=xleft
  rest=lrest
  do gp=1,gpcur
! copy the "ceq" equilibrium to a new record in gpp(gp)%eqp
! the list of conditions are not copied but linked to original equilibrium
     call tqcceq(gpname,gpix(gp),gpp(gp)%eqp,ceq)
     if(gx%bmperr.ne.0) goto 1000
! To ensure the conditions are unique in each gridpoint nullify the list
! (it creates som lost memory ...)
     nullify(gpp(gp)%eqp%lastcondition)
! Set new conditions in equilibrium for this gridpoint
     call tqsetc('T ',0,0,tpval(1),cnum(1),gpp(gp)%eqp)
     if(gx%bmperr.ne.0) goto 1000
     call tqsetc('P ',0,0,tpval(2),cnum(2),gpp(gp)%eqp)
     if(gx%bmperr.ne.0) goto 1000
     call tqsetc('N ',0,0,nval,cnum(3),gpp(gp)%eqp)
     if(gx%bmperr.ne.0) goto 1000
! halfway change compositions
     if(gp.gt.half) then
        xval=xright
        rest=rrest
     endif
     do iel=1,nsel
        if(iel.ne.rest) then
           call tqsetc('X ',iel,0,xval(iel),cnum(4),gpp(gp)%eqp)
        endif
     enddo
     if(gx%bmperr.ne.0) goto 1000
! this is a way to check the gridpoint conditions
!     write(*,*)'Conditions for equilibrium ',gp
!     call tqlc(kou,gpp(gp)%eqp)
!     if(gx%bmperr.ne.0) goto 1000
! calculate without gridmin
! because start values for the phases were copied above
     call tqce(' ',-1,0,zero,gpp(gp)%eqp)
     if(gx%bmperr.ne.0) goto 1000
!
     write(*,44)gp,trim(gpname),(xval(iel),iel=1,nsel)
44   format(i5,6x,a,10F7.4)
     do iel=1,nsel
        gpp(gp)%xval(iel)=xval(iel)
     enddo
! increment the equilibrium name index
     call incname(gpname,6)
  enddo
! write intitial profile for all 3 components for all gridpoints
  nt=0
  nrow=1
  write(*,*)'Initial composition profile:',nsel,gpcur
! 
  do iel=1,nsel
     write(*,51)(gpp(gp)%xval(iel),gp=1,gpcur)
  enddo
  if(plut.gt.0) then
     write(plut,'("# ",2i10)')nt,nrow
     do iel=1,nsel
        write(plut,51)(gpp(gp)%xval(iel),gp=1,gpcur)
     enddo
  endif
51 format(100(0pF6.3))
52 format(i10,': ',100(0pF6.3))
! 
  if(parallel) then
     write(*,*)'Calculations will be made in parallel'
  else
     write(*,*)'Calculations will be made sequentially'
  endif
!--------------------------------------- simulation starts
  call cpu_time(xxx)
  call system_clock(count=cc1)
! Take a time step, calculate diffusion and modify compositions
! Use the difference in chemical potential between two adjacent gridpoints
! to calculate the flow of elements
  simulate: do while(.TRUE.)
     nt=nt+1
     dxmax=zero
     diff: do gp1=1,gpcur-1
        gp2=gp1+1
        sumneg=zero
        sumpos=zero
        do jj=1,nsel
           gpxval(jj)=gpp(gp1)%xval(jj)
           dmu=gpp(gp2)%eqp%cmuval(jj)-gpp(gp1)%eqp%cmuval(jj)
           dx(jj)=mobi(jj)*dmu
           if(abs(dx(jj)).gt.1.0d-2) then
              write(*,*)'Very strong diffusion!',nt,jj,dx(jj)
           endif
! dxmax is used to check convergence, if max dxmax small then terminate
           if(abs(dx(jj)).gt.dxmax) dxmax=abs(dx(jj))
           if(dx(jj).gt.zero) then
              sumpos=sumpos+dx(jj)
           else
              sumneg=sumneg-dx(jj)
           endif
        enddo
! The sum of the fractions should always be unity, make the sum of all
        if(sumpos.le.1.0D-12 .or. sumneg.le.1.0D-12) then
! There is no diffusion
           sdxp=zero; sdxn=zero
        elseif(sumpos.gt.sumneg) then
! scale the maximal flow to be the same as the minimal
           sdxp=sumneg/sumpos
           sdxn=one
        else
           sdxp=one
           sdxn=sumpos/sumneg
        endif
! move the atoms!!
        do jj=1,nsel
           if(dx(jj).ge.zero) then
              gpp(gp1)%xval(jj)=gpp(gp1)%xval(jj)+dx(jj)*sdxp
              gpp(gp2)%xval(jj)=gpp(gp2)%xval(jj)-dx(jj)*sdxp
           else
              gpp(gp1)%xval(jj)=gpp(gp1)%xval(jj)+dx(jj)*sdxn
              gpp(gp2)%xval(jj)=gpp(gp2)%xval(jj)-dx(jj)*sdxn
           endif
        enddo
! Check fractions are in range and sum is unity
        sum=zero
        do jj=1,nsel
           if(gpp(gp1)%xval(jj).ge.one .or.gpp(gp1)%xval(jj).le.zero) then
              write(*,69)gp1,jj,gpp(gp1)%xval(jj)
69            format('Fraction outside limits at gridpoint ',2i3,F10.4)
!              stop 'katastrof 1!'
              if(gpp(gp1)%xval(jj).ge.one) gpp(gp1)%xval(jj)=1.0-1.0D-8
              if(gpp(gp1)%xval(jj).le.zero) gpp(gp1)%xval(jj)=1.0D-8
           endif
           sum=sum+gpp(gp1)%xval(jj)
        enddo
        if(abs(sum-one).gt.1.0D-7) then
           write(*,*)'Sum of fractions not unity at gridpoint ',gp1,sum
           write(*,'(a,3F10.5)')'Fractions: ',(gpxval(jj),jj=1,3)
           write(*,'(a,3F10.5)')'Fractions: ',(gpp(gp1)%xval(jj),jj=1,3)
           stop 'katastrof 2!'
        endif
     enddo diff
! modmod controls output
!  initially write each profile, then every 10, then every 100, then 1000
     if(nt.eq.10) then
        modmod=10
     elseif(nt.eq.100) then
        modmod=100
     elseif(nt.eq.1000) then
        modmod=1000
     elseif(nt.eq.10000) then
        modmod=10000
     endif
     if(modmod.gt.maxloop) modmod=maxloop
     if(mod(nt,modmod).eq.0) then
        if(plut.gt.0) then
           modmod=2*modmod
           nrow=nrow+1
           write(plut,'("# ",2i10)')nt,nrow
           do iel=1,nsel
              write(plut,51)(gpp(gp)%xval(iel),gp=1,gpcur)
           enddo
        endif
        write(*,*)'Done ',nt,' timesteps'
! debug output
!        write(*,50)nt,dxmax,(gpp(gp)%xval(1),gp=1,gpcur)
!        do iel=2,nsel
!           write(*,51)(gpp(gp)%xval(iel),gp=1,gpcur)
!        enddo
     endif
!---------------------------------------
! the equilibrium with the new composition in all gridpoints (parallel?)
     rest=lrest
     newx: do gp1=1,gpcur
        gridpoint=>gpp(gp1)
        if(gp1.gt.half) then
           rest=rrest
        endif
        do iel=1,nsel
           if(iel.ne.rest) then
              call tqsetc('X ',iel,0,gridpoint%xval(iel),cnum(4),gridpoint%eqp)
              if(gx%bmperr.ne.0) then
                 write(*,*)'Error setting condition: ',gridpoint%eqp%eqname
                 gx%bmperr=0
              endif
           endif
        enddo
! possible debug output of conditions ....
!        write(*,*)'Conditions for equilibrium ',gp1
!        call tqlc(kou,gridpoint%eqp)
!        if(gx%bmperr.ne.0) goto 1000
     enddo newx
! alternative running without parallel
     call cpu_time(eqcpu1)
     call system_clock(count=eqcc1)
     pos: if(parallel) then
!$OMP parallel do private(gridceq)
        neweq: do gp1=1,gpcur
           gridceq=>gpp(gp1)%eqp
! This give new chemical potentials
! the number of threads must be obtained inside the loop otherwize 1 or 0
!$        howmany=omp_get_num_threads()
! calculate without gridmin
           call tqce(' ',-1,0,zero,gridceq)
           if(gx%bmperr.ne.0) then
              write(*,*)'Error calculating equil: ',gridceq%eqname
              gx%bmperr=0
           endif
        enddo neweq
!$omp end parallel do
     else
        howmany=1
        seqloop: do gp1=1,gpcur
           gridceq=>gpp(gp1)%eqp
! Calculate equilibrium for new chemical potentials for next diffusion calc
           call tqce(' ',-1,0,zero,gridceq)
           if(gx%bmperr.ne.0) then
              write(*,*)'Error calculating equil: ',gridceq%eqname,gx%bmperr
              gx%bmperr=0
           endif
        enddo seqloop
     endif pos
     call cpu_time(eqcpu2)
     call system_clock(count=eqcc2)
     eqcpusum=eqcpusum+eqcpu2-eqcpu1
     eqccsum=eqccsum+eqcc2-eqcc1
! loop back until simulation timestep exceeded or no change in composition
     if(abs(dxmax).lt.1.0D-5 .or. nt.gt.maxloop) exit simulate
  enddo simulate
!---------------------------------------
! Output of results and repeat input (if already forgotten)
  call system_clock(count=cc2)
  call cpu_time(xxy)
  write(*,190)nt,dxmax,xxy-xxx,cc2-cc1
  write(*,193)eqcpusum,eqccsum,howmany
  write(*,17)trim(database),(selel(iel),iel=1,nsel)
  write(*,18)trim(phasename),gpcur,(xleft(iel),iel=1,nsel)
  write(*,19)(xright(iel),iel=1,nsel)
  write(*,20)(mobi(iel),iel=1,nsel)
  write(*,21)tpval(1),maxloop,dff
  if(plut.gt.0) then
! Results saved on file, finish the graphics output and CPU times
     write(plut,195)nsel
195  format('EOD'//&
          'set style line 1 linetype 1 linecolor rgb "#000000" linewidth 1',&
          ' pointtype 10'/&
          'set style line 2 lt 1 lc rgb "#4169E1" lw 1 pt 7'/&
          'set style line 3 lt 1 lc rgb "#00C000" lw 1 pt 6'//&
          'plot for [myRow=0:',i3,'] $profile  matrix using 1:3 ',&
          ' every :::myRow::myRow with linespoints linestyle 1+myRow ',&
          ' title sprintf("Row number %d",myRow)'//'pause mouse')
     write(plut,190)nt,dxmax,xxy-xxx,cc2-cc1
190  format(/'# Timesteps:',i6,', Dmax=',1E12.4/&
          '# CPU time:',F12.4,'s, clockcycles: ',i8)
     write(plut,193)eqcpusum,eqccsum,howmany
193  format(/'# For equilibrium calculation: CPU time ',F12.4,' s and cc: ',i8/&
          '# Number of thread(s): ',i3)
     close(plut)
  endif
!
! Remind the file name for graphics
  write(*,'(/"Graphics output on ",a)')trim(profiles)
  write(*,991)
991 format(/'All well that ends well'/)
  stop
!-----------------------------------------
! OC error
1000 continue
  write(*,1001)gx%bmperr,trim(bmperrmess(gx%bmperr))
1001 format('Error code ',i5/'Message: ',a)
  stop
!
end program sim2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

subroutine incname(name,pos)
! increment sequential index in equilibrium name
  implicit none
  character name*(*)
  integer pos,jpos,iv
  jpos=pos
  loop: do while(jpos.gt.1)
     iv=ichar(name(jpos:jpos))-ichar('0')
     if(iv.lt.9) then
        iv=iv+1
        name(jpos:jpos)=char(iv+ichar('0')); exit loop
     else
        name(jpos:jpos)='0'
        jpos=jpos-1
     endif
  enddo loop
  return
end subroutine incname

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

subroutine readline(setup,line)
! read seup file ingnoring comment lines
  implicit none
  integer setup
  integer :: nl=0
  save nl
  character line*(*)
  do while(.TRUE.)
     read(setup,10,end=1000)line
10   format(a)
     nl=nl+1
!     write(*,*)'Echo input line: "',trim(line),'"'
     if(line(1:1).ne.'#') return
  enddo
1000 continue
  write(*,*)'Found EOF of setup file after line ',nl
  return
end subroutine readline

