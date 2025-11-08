!----------------------------
! Reading the XTDB file with AppendXTDB: Models, parameters, TPfun and Biblio
! 1. All data on original file.  It may have to be rewinded for missing TPfuns
! 2. One can select a subset of elements.  All species which can form from
!    this subset are entered.
! 3. An external Model file, this is read initially and never again
! 4. The original file is then read to the end and phases, parameters entered
!    If TPfun missing the file may be rewinded and reread severel times
!    until the number of missing TPfuns/bibitems are zero or constant.
! 5. One can have an external parameter file which will be read once.
!    Needed parameters are entered and file closed.  Preferably no TPfuns
! 6. One can have an external TPfuns file which is read and rewinded
!    to extract missing TPfuns until no TPfun entered from this file.
!    Warning of duplicates of same TPfuns.
! 7. One can have a file with biblitems which is opende if there are
!    missing bibitems at the end.
! 8. If there ate missing TPfun or bibitem at the end that is reported.
!----------------------------
!
! The values for modelappy, parappy, tpfunappy, biblioappy initiates to 0
! if there are no associated file:
! modelappy=1 if there is an AppendXTB Model file.  Read once and set to zero
! parappy=1 if there is an AppendXTB parameter file.  Read once and set to zero
! tpfunappy=1 if there is an AppendXTB TPfun file.  Read once and 
!     incremented by one for each rewind.  Stop rewinding and set to 0
!     when number of missing TPfuns is zero or do not decrease.
! biblioappy=1 if the is a biblogapy file.  Read once and set to zero.
! allappy is 0 initially.  It is set to 1 while reading parappy and tpfunappy
!     files.  If the masterfile is rewinded allappy is set to 2 o 3
!     to prevent reading anythng except TPfun/Trange and bibitem tags.
!
!----------------------------
! To do:
! 1. Clean up a bit, remove duplicate variables
! 2. Reread masterfile for TPfuns if no AppendXDB file
! 3. Store complete attribute for each tag in selxx temporary storage

!------------------------------
!
! Module xtdblib
! program to read an xml file, in particular XTDB
!  implicit none
! this contais tags, attributes and some global variables
! #include "gtp3_xml.F90"
!
!  integer nselph
!
!==================================================================
!  contains
!==================================================================
!    
! subroutines and functions
!
! xtdbread(masterfile) 
! read a whole database.  This may open and read other AppendXTDB files
!
! xtdbtag(unit,fline,tagname,matt,pretag,attributes)
! extract all attributes of tag from the file in sequential order.
! fline is the line number.  The tag can depend om other nested tags.
!
! logical function getatt(line,ip,attname,value)
! extract sequentially the attributes and their values within "" from line
! The value should be all spaces and for the first attribute ip should be 1
!
! - The important tags on a single XTDB file are:
!    Model, Element, Species, Phase (and nested tags), Parameter, Tpfun, Trange 
! but as a TPfun may use other TPfun the file may be rewinded to find all
! - A first sequential read will pick up all phases parameters selected but
! as TPfuns may call other TPfuns and it may be necessary to rewind the file.
! - The program keeps track of all needed TPfuns and may rewind to find them.
! - All bibliographic references are stored until reading the bibliography.
! - An AppendXTDB file for models is recommended if the database may be
! used by different software.
!
! For large databases the the primary XTDB file should contain all
!       AppendXTDB, Element, Species, Phase (and nested tags)
! On separate AppendXTDB file one can have:
! a file sith model information 
! a file with Parameter and Trange tags (which is read only once)
! a file with TPfun and Trange tags which is rewinded until all TPfuns extracted
! a file with bibliograpgic refenences (read once at the end)
! The Parameter file can contain TPfuns but it will not be rewinded
! The TPfun file may contain Parameters but it should be as short as possible
! as it may be rewinded several times.
! An AppendXTDB file should start with <Appendix> and end with </Appendix>
!
! The XTDB tags, attributes and some global data and variables 
! are in gtp3_xml.F90
!
! Some routines in this file:
! subroutine xtdbread(filename,nel,elnames)
! open filename and read an XTDB file including ApendXTDB files and
! extracts tag and data for nel elemsnte for storing in thermodynamic softaware
! if nel=0 the elements in the satabase is returned
!
! subroutine xtdbtag(unit,fline,tagname,matt,pretag,attributes)
! reads lines from unit until the end of attributes of a tag
! If matt=1 a nested tag found with tagname and attributes
!    matt=0 a complete tag found with tagname and attributes
!    matt=-1 the end of a nested tag found and some action needed
! Maybe pretag in not needed inside xtdbtag?
!
! logical function getatt(attributes,ip,attname,values)
! extracts sequentially from attributes the attname and its values using ip
! True if it finds and attribute
!
! logical function xeolch(line,ip)
! Skipping spaces and TAB characters in "line" from "ip".  True if no data
!
!* logical function check_mpid(mpid,phase)
! Check that phase has a model corresponding to the MPID of a parameter
!
!--------------------------------------------------------
  
  subroutine xtdbread(masterfile,nel,el)
! To read xtdb files and extract tags and attributes including nesting
    implicit none
    character*(*) masterfile
    character*2 el(*)
    integer nel
! if nel=0 extract all elements from database
!----------------------------
! Reading the XTDB file with AppendXTDB: Models, Parameters, TPfun and Biblio
! 1. All data on original file.  It may be rewinded for missing TPfuns
! 2. If there is an external Model file it is read and closed before any data.
! 3. The original file is then read to EOF and phases, parameters and 
!    missing TPfuns and bibliographic refs are read. File will not be closed!
!    DUPLICATE TPfuns will always be reported as errors, keeping first case.
! 4. If there is an external parameter file (WITHOUT ANY TPFUNS) this is
!    opened and parameters needed are entered and the file closed at EOF.
!    Any TPfuns found will be create a warning but entered if needed.
!    This file will not be rewound. 
! 5. If there is an external TPfuns file this is opened and read to extract
!    missing TPfuns, it may be rewinded several times until no 
!    missing TPfun found and then it will be closed.
! 6. If there are TPfuns or biblitems missing the original XTDB file will
!    be rewinded and read again until no missing TPfuns or bibitems found
!    and the file closed
! 7. If TPfun or bibitem missing now that will be reported.
!----------------------------
! line is line read from file
    character line*256,tagname*64
! tag attributes
!    character(len=:), allocatable :: attributes ! this does not work
!    character attributes*(1024)
    character(len=:), allocatable :: attributes
! the global character  wholexpr is concatinated Expr for TPfun or Parameter
! The different files used for reading the XTDB files
! The database can be a single (primary) file or split on several files.
!    All elements and phases and extra files must be in the primary file
! The extra files can define models, parameters, tpfuns and bibliography
! The model file is not needed but defines the MPIDs
! The parameter file should only have Parameter and Trange tags
! The tpfun file should only have TPfun and Trange tags
! The bibliography file should only have bibliograpgy and bibitem tags
    integer unit   ! the current file used by xtdbget, any of those below
    integer, parameter :: unit1=21  ! main XTDB file, 
                   ! with everyting or just element and phases
                   ! open the whole time, may be rewinded at the end
    integer, parameter :: unit2=22 ! parameter file read once and then closed
    integer, parameter :: unit3=23 ! tpfun file may be rewinded several times
!
! for nested tags level is current level,  matt=0 if tagend missing (nesting)
    character pretag*24,values*256,attname*18
    integer tagno,matt
! for storing data during nested tags, each at a higher level ?? needed
    integer, parameter :: maxatt=5
    character*256, dimension(maxatt) ::  saveatt
! tags with nested tagas are TPfun 1: Trange
!    character phaseid*24
! Phase sublattices, 10xconstituents, crystalstructure, amendphase, disord ...)
    character phid*24,phconfig*24,phmodels*128,phstate*1
! check of missing tpfun and missingbib
    integer prevmisstp,prevmissbib
! extending an allocated variable ??
!
!  extending an allocatable and already allocated array (data(m))
!     data = [ data, ( 0, kk=1,n ) ]     see ternary extrapolation!
!      

!    curent line in a database file (mater or AppendXTDB
    integer fline
! for tempraty storage of phase data
    integer phnsub,phisub
    character (len=:), allocatable :: lowt, default_lowt,tpfun
    character (len=:), allocatable :: hight,default_hight,expr
    character (len=:), allocatable :: parid,bibref
    character*16 usedtpfun   ! is tpfun found in Expr of Parameter or TPfun

    logical listphases,add,ternaryxwarning

! handling of AppendXTDB
! The Models must be read first and only once
! The Parameters after reading the phases and only once
! The TPfun after all paramaters and maybe several times
! The Bibiliography read last and and only once
! Bibilographic references
    character(len=:), allocatable :: currentfile
    character ch1*1,spel*24,spname*24,stoichiometry*80,mqmqa*40,uniquac*12
    character stoisp*40,refstate*40,phname*24
    character(len=:), allocatable :: addon
    character(len=:), allocatable :: xpoldata
! for disorderedpart
    character (len=:), allocatable :: sum,disph,sub
    character*2 element
!
! for various purposes
    double precision mass, h298, s298
    integer ip,iq,ir,is,it,level,prevlev,lk,kk,semicolon,lph,jp
    integer missbib,unknowntp,jel,skip,rewinds
    type(octerxpol), allocatable, target :: terxpol
!
! initiate position in arrays for storing selected information
! We automatically introdce Va
! We may have to order elements and species alphabetically!!
    nselel=0; nselsp=1; nselph=0; nselpar=0; nseltp=0; nselbib=0

! software defaults
    default_lowt='298.15'
    default_hight='6000'
    debug=.false.
    rewinds=0
    ternaryxwarning=.true.
!
!    write(*,*)'3EY xtdbread: ',trim(masterfile),nel
! data for missing and found AppendXTDB, initiated to 0 below
! when specified set -1, while reading +1, when read 0
! allappy initially set to 0, set to 1 when any other AppendXTDB file read
    if(nel.eq.0) then
!       write(*,*)'Just to know which elements are in the database'
       ignorEOT=.true.
    else
       ignorEOT=.false.
       do jel=1,nel
          call capson(el(jel))
       enddo
       write(*,8)(el(jel),jel=1,nel)
8      format(/'Extract data for: ',20(a,x))
! array for bibliographic references, selbib, for selected parameter below
! allocate arrays for storing data, related to number of elements
       listphases=.true.
       maxtdbel=nel+2; maxtdbsp=10*nel*(nel-1); maxtdbph=5*nel*(nel+1); 
       maxpar=20*nel*(nel+1)*(nel+2)
       maxtp=30*nel*(nel+1); maxbib=5*nel*(nel+1)*(nel+2)
       if(allocated(selel)) then
! deallocate all
          deallocate(selel)
          deallocate(selsp)
          deallocate(selph)
          deallocate(selpar)
          deallocate(seltpfun)
          deallocate(selbib)
       endif
       allocate(selel(-1:maxtdbel))
       if(debug) write(*,9)'Dimensioning:',&
            maxtdbel,maxtdbsp,maxtdbph,maxpar,maxtp,maxbib
9      format(a,6i7)
! Always include Va and /- as element -1 and 0
       selel(-1)%elname='/-'; selel(-1)%data=' '
       selel(0)%elname='VA'; selel(0)%data=' '
       allocate(selsp(maxtdbsp))
!       allocate(selspord(maxtdbsp)) allocated later
       selsp(1)%species='VA'; selsp(1)%data='VA'; selsp(1)%charge=0.0D0
! needed for Va as species
       allocate(selsp(1)%elnames(1));
       allocate(selsp(1)%stoicc(1));
       selsp(1)%elnames='Va';  selsp(1)%stoicc(1)=1.0D0
!
       allocate(selph(maxtdbph))
       allocate(selpar(maxpar))
       allocate(seltpfun(maxtp))
       allocate(selbib(maxbib))
    endif
! When just returning all elements
    jel=0
! AppendXTB files
    allappy=0
    modelappy=0; parappy=0; tpfunappy=0; biblioappy=0
    prevmisstp=0; missingtp=0; prevmissbib=0; missingbib=0
!
    nselbib=0
!
!    write(*,*)'Opening XTDB file: ',trim(masterfile)
    unit=unit1
    fline=0
    open(unit,file=masterfile,&
         access='sequential',form='formatted',status='old')
! zero main file line number
    currentfile=masterfile
! zero position in characters containing attributes
    attpos=0
    attributes=' '
    tpfun=' '
! initiate tag nesting
    level=0
    matt=0
    pretag=' '
    xtdberr=0
    nomorelements=.FALSE.
!
! big loop for reading everything
!==============================================================
    readall: do while(.true.)
! only one tag per line, attributes can be on following lines
       if(line(1:2).eq.'<?') then
          write(*,*)'Ignoring tag <?'
          fline=fline+1
          cycle readall
       endif
!--------------------------
! the call below will extract lines until the end of arrributes of a tag
17     continue
       if(fline.lt.0) then
          write(*,*)'3EY fline',fline
          exit readall
       endif
!       write(*,*)'3EY len(attributes): ',len(attributes),len(tpfun)
       call xtdbtag(unit,fline,tagname,matt,pretag,attributes)
!       write(*,18)unit,matt,xtdberr,trim(tagname),trim(pretag),&
!            len(attributes),len(tpfun)
18     format('Read xtdbtag: ',3i5,' "',a,'" "',a,'" ',i5)
       if(xtdberr.ne.0) then
! we may have end of file (4700) reading an AppendXTDB files here ...
          if(nel.eq.0) then
             nel=jel
          elseif(xtdberr.ne.4700) then
             write(*,*)'3EY Back with xtdberr: ',xtdberr,parappy,tpfunappy
             goto 999
          endif
          xtdberr=0
       endif
! if nel=0 just extract elements, do not bother about anything else
! but read to endoffile
       if(nel.eq.0) then
          if(tagname.eq.'Element ') then
             jel=jel+1
!             write(*,*)'nel: ',nel,jel,' tagname: ',trim(tagname)
             ip=1; values=' '
             do while(getatt(attributes,ip,attname,values))
                if(attname(1:2).eq.'Id') el(jel)=trim(values)
             enddo
          endif
          cycle readall
! end reading file just to obtain the elements
       else
! we must read all elements and species before the phases
          nomorelements=.FALSE.
          if(tagname.eq.'Element ') then
             refstate='SER'; 
             ip=1; values=' '
             do while(getatt(attributes,ip,attname,values))
                ignore: if(attname(1:2).eq.'Id') then
! Va and /- are preselected
                   if(values(1:2).eq.'VA'.or.values(1:2).eq.'/-') cycle readall
                    do jel=1,nel
                      if(values.eq.el(jel)) exit ignore
                   enddo
! this element not in the selection, ignore
!                   write(*,*)'Ignoring: ',trim(values),len_trim(values)
                   cycle readall
                endif ignore
                element=values
! this element is selected, extract its data, 
! what about the same element several times in the xtdbfile?
                if(attname(1:8).eq.'Refstate') refstate=trim(values)
                if(attname(1:4).eq.'Mass') then
                   jp=1
                   call getrel(values,jp,mass)
                   if(gx%bmperr.ne.0) then
                      write(*,*)'3EY Illegal mass for ',trim(values)
                      gx%bmperr=0; mass=1.0D0
                   endif
                   if(mass.le.zero) mass=1.0D0
                endif
                if(attname(1:4).eq.'H298') then
                   jp=1
                   call getrel(values,jp,h298)
                   if(gx%bmperr.ne.0) then
                      write(*,*)'3EY Illegal h298 for ',trim(values)
                      gx%bmperr=0; h298=0.0D0
                   endif
                endif
                if(attname(1:4).eq.'S298') then
                   jp=1
                   call getrel(values,jp,mass)
                   if(gx%bmperr.ne.0) then
                      write(*,*)'3EY Illegal S298 for ',trim(values)
                      gx%bmperr=0; s298=0.0D0
                   endif
                endif
! all attributes saved here if element selected
!                mass=1.0D0; h298=0.0D0; s298=0.0D0
                call OCenterel(element,tpfun,refstate,mass,h298,s298)
             cycle readall
             enddo
          endif
       endif
!---------------------------
! matt= +1 new tag nest, 0 tag finished, -1 tag nest finished
!       write(*,44)trim(tagname),trim(pretag),matt,fline
44     format('44 Back from xtdbtag: "',a,'"  pretag: "',a,'"',i3,i7)
!
       if(matt.eq.1) then
! if matt=1 means new tag with attributes but no endoftag, increase level
          level=level+1
          if(level.gt.maxatt) then
             write(*,46)fline,(ip,endoftag(ip),ip=1,maxatt)
46           format('3EY Line ',i7,' Tags: ',5(i3,2x,a))
          endif
          saveatt(level)=attributes
          endoftag(level)=tagname
          pretag='</'//trim(tagname)//'>'
!========================= handling end of nested tags
       elseif(matt.eq.-1) then
!========================= nested TPfun =========
! if matt=-1 the line is end of nested pretag, decrease level !!??
!          write(*,*)' >>>Nested tag end: "',&
!               trim(tagname),'" pretag: "',trim(pretag),'" ',level,fline
          if(trim(pretag).eq.'</XTDB>' .or. trim(pretag).eq.'</Appendix>') then
!================== end of file: </XTDB> or </Appendix>
! Here we are at the end of the primary file or an AppendXTDB file
! The ModelXTDB is read when found, the other will be opened for read here
             if(debug) write(*,*)'End of file for: ',trim(currentfile)
             prevmisstp=missingtp
             missingtp=0
             do ntp=1,nseltp
                if(seltpfun(ntp)%status.lt.0) missingtp=missingtp+1
             enddo
             prevmissbib=missingbib
             missingbib=0
             do ntp=1,nselbib
                if(selbib(ntp)%status.lt.0) missingbib=missingbib+1
             enddo
             if(debug) then
                write(*,*)'Number of TPfuns found and missing:',nseltp,missingtp
                write(*,*)'Number of bibitem found and missing:',&
                     nselbib,missingbib
! The principal XTDB file is not closed, it may be rewinded in the end
                write(*,47)parappy,tpfunappy,biblioappy,allappy
47              format('Any AppendXTDB files to read?',4i3)
                write(*,*)'Press return to continue 1',tpfunappy
                read(*,'(a)')ch1
             endif
! if there is an AppendXTDB for Parameters read that now!
! The Parameter file is done first (well actually the ModelAPpendXTDB)
             if(parappy.eq.1) then
                if(debug) write(*,*)'Appending Parameter file: ',trim(parappx)
!----------------------------------------------------------
! allappy set to unity indicate only parameter, tpfun and trange tags allowed
                allappy=1
! parappy set to 0 as it it read only once
                parappy=0
                unit=unit2
                fline=0
                open(unit,file=parappx,access='sequential',&
                     status='old',form='formatted')
                if(debug) write(*,*)' *** Opened ',trim(parappx),' line',fline
                currentfile=parappx
! Now read parameters from appendixfile
                cycle readall
             elseif(tpfunappy.eq.1) then
!----------------------------------------------- tpfuns
! if there is an AppendXTDB for Tpfuns read that, maybe including rewinds
! This file may have to be rewinded if there are unknown TPfuns
! missingtp and missingbib set to zero initially
                ignorEOT=.true.
                allappy=2
                tpfunappy=tpfunappy+1
                if(unit.eq.unit2) then
! We may not have opend any parameterappendix but of so, close it
                   if(debug) write(*,*)'Closing append file: ',trim(parappx)
                   close(unit2)
                endif
                if(debug) write(*,*)'Appending TPfun file: ',trim(tpfunappx)
                unit=unit3
                fline=0
                open(unit,file=tpfunappx,access='sequential',&
                     status='old',form='formatted')
                if(debug) write(*,*)' *** Opened ',trim(tpfunappx)
                currentfile=tpfunappx
! Now read TPfuns from appendixfile
                cycle readall
             elseif(tpfunappy.ge.2) then
!----------------------------------------- maybe rewind the TPfun file?
! we may have to rewind the TPfunappx file if there are unknown TPfuns
                if(debug) write(*,401)nseltp,prevmisstp,missingtp
401             format('Of ',i4,' TPfun missing changed from ',i4,' to ',i4)
                if(missingtp.gt.0 .and. missingtp-prevmisstp.gt.0) then
                   prevmisstp=missingtp
                   if(debug) then
                      write(*,*)'Press return to Rewind ',trim(tpfunappx)
                      read(*,'(a)')ch1
                   endif
! Rewind appendixfile and read again (maybe several times)
! as TPfuns may use other TPfuns
                   fline=0
                   rewind(unit)
                   cycle readall
                else
! The number of missing TPfuns does not decrease after rewinding
! Either all TPfuns found or the missing are not on tpfunappx
                   if(debug) write(*,*)'Closing tpfunfile: ',trim(tpfunappx)
                   tpfunappy=0
                   close(unit)
                endif
! if there is a biblioappy read this here 
! It is also read in the following elseif but that is outside this scope
                if(biblioappy.gt.0) then
                   fline=0
                   if(debug) write(*,410)trim(biblioappx)
410                format(/'Appending bibliography 1: ',a)
! the file is opended and read separatelym only bibitems allowed
                   biblioappy=0
                   call xtdbbiblio(biblioappx)
                   if(debug) write(*,*)'Back from reading bibliography'
                   missingbib=0
                   do ntp=1,nselbib
                      if(selbib(ntp)%status.lt.0) missingbib=missingbib+1
                   enddo
! if there are missing tpfuns and/or biblographics rewind and read the main file
                   if(missingtp.gt.0 .or. missingbib.gt.0) then
                      if(debug) write(*,491)missingtp,missingbib,masterfile
491                   format('Missing data, ',2i4,' rewinding: ',a)
! allappy=3 means only Bibitem and TPfun tags read
                      ignorEOT=.true.
                      allappy=3
                      fline=0
                      unit=unit1 
                      currentfile=masterfile
                      rewind(unit)
                      cycle readall
                   endif
                elseif(missingtp.gt.0 .or. missingbib.gt.0) then
! a final try rewinding and rewind and read again from primary file
                   rewinds=rewinds+1
!                   write(*,*)'Rewinding masterfile A',rewinds
                   unit=unit1
                   currentfile=masterfile
                   rewind(unit)
                   cycle readall
                endif
             elseif(biblioappy.eq.1) then
! We can be here if there is no AppendXTDB for TPfuns
! if there is an AppendXTDB for bibliography read that now!
! There may also be a <Bibliography> tag in the main file? 
                if(debug) write(*,*)'Appending bibliography: ',trim(biblioappx)
                biblioappy=0
                call xtdbbiblio(biblioappx)
                if(debug) write(*,*)'Back from reading bibliography'
                prevmisstp=missingtp
                missingtp=0
                do ntp=1,nseltp
                   if(seltpfun(ntp)%status.lt.0) missingtp=missingtp+1
                enddo
                prevmissbib=missingbib
                missingbib=0
                do ntp=1,nselbib
                   if(selbib(ntp)%status.lt.0) missingbib=missingbib+1
                enddo
                if(missingtp.gt.0 .and. missingtp-prevmisstp.gt.0) then
! to avoid rewind masterfile forever
                   prevmisstp=0
                   fline=0
                   if(debug) write(*,*)'Found missing data, rewinding: ',&
                        masterfile,fline
! we must ignore all tags except TPfun/Trange and Bibitem
                   allappy=3
                   rewinds=rewinds+1
!                   write(*,*)'Rewinding masterfile B',rewinds
                   unit=unit1
                   currentfile=masterfile
                   rewind(unit)
                   cycle readall
                else
! time to close when TPfuns as well as masterfile rewinded
                   fline=0
                   if(debug) write(*,321)trim(masterfile)
321                format('3EY Finished reading XTDB file: ',a)
                   unit=unit1
                   goto 990
                endif
             else
! reading bibliography AppendXTB if there is no AppendXTDB for TPfuns
                if(debug) write(*,*)'Applied all AppendXTDB (or none)'
                missingtp=0
                do ntp=1,nseltp
                   if(seltpfun(ntp)%status.lt.0) missingtp=missingtp+1
                   enddo
                prevmissbib=missingbib
                missingbib=0
                do ntp=1,nselbib
                   if(selbib(ntp)%status.lt.0) missingbib=missingbib+1
                enddo
                if(rewinds.gt.0.and.missingtp.eq.0) write(*,324)trim(masterfile)
324             format('All TPfuns found, closing ',a)
!                if(missingtp.gt.0 .and. missingtp-prevmisstp.gt.0) then
                if(missingtp.gt.0 .and. rewinds.lt.3) then
                   allappy=3
                   fline=0
                   rewinds=rewinds+1
                   write(*,322)missingtp-prevmisstp,trim(masterfile),rewinds
322                format('3EY Missing ',i3,' TPfuns, rewinding: ',a,i3)
! try reading the masterfile once again, 
!                   write(*,*)'Rewinding masterfile C',rewinds
                   currentfile=masterfile
                   unit=unit1
                   rewind(unit)
                   cycle readall
                else
! time to close when there is no TPfuns file and masterfile has been rewinded
                   if(debug) then
                      write(*,*)'Finished 2, closing: ',trim(masterfile)
                      write(*,491)missingtp,missingbib,masterfile
                      write(*,321)trim(masterfile)
                   endif
                   unit=unit1
                   write(*,*)'Closing file 1'
                   goto 990
                endif
             endif
!========================= Below handling end of nested tags
          elseif(trim(pretag).eq.'</TPfun>') then
! a nested tag has the end of tag on a single line
!======================== end of nested </Tpfun> and Trange
! simple TPfun without Trange are entered below
! check if the TPfun should be add(ed) (used in a Parameter or other TPfun
             call xtdbentertpfun(tpfun,add)
             if(add) then
                call OCentertpfun(tpfun)
! the expression is stored in the global wholexpr
             endif
          elseif(trim(pretag).eq.'</Phase>') then
!======================== end of Phase, always nested
! The end of Phase tag (always nested)             
             wholexpr=phrec%confent//' '//phrec%state//' '//phrec%clist(1)%list
! inside OCenterphase the phase may be added to selph if constituents entered
             call OCenterphase(phrec)
!             write(*,201)fline,phrec%id,phrec%confent,phrec%state,&
!                  phnsub,phrec%mult
201          format(/'Collected all phase data: ',i7/&
                  'Phase name: ',a,5x,a,5x,a,5x,i2,5x,a)
!             do ip=1,phnsub
!                write(*,202)phrec%clist(ip)%subx,phrec%clist(ip)%list
202             format('Sublattice: ',a,'  Constituents: ',a)
!             enddo
!>>>>>>> create record for phase and store data in OC <<<<<<<<<<<<<
! clear for next phase by deallocate the whole phrec
             deallocate(phrec)
          elseif(trim(pretag).eq.'</Sublattices>') then
!======================== end of Sublattices tag, part pf Phase
! Update level here for next tag, done separately here as code is messy
! Note that the <Constituent tag is inside sublattices without nested tags
! Data will be taken care of at the end of the <Phase tag 
!             write(*,*)'  End of nested Sublattice tag: "',level,matt,phnsub
             continue
          elseif(trim(pretag).eq.'</Parameter>') then
!========================= end of nested Parameter, may include Trange
! The end of a nested Parameter tag, non-nested parameters treated below 
! A nested parameter contains Trange records
!             write(*,*)'Found end of nested </Parameter> add ',bibref
! the ' N ' to be compatible with TDB files ... add bibref inside xtdbOCfun
! NOTE wholexpr is an allocated character
             wholexpr=wholexpr//' N'
             if(.not.allocated(bibref)) bibref=' '
             call OCenterpar(parid,skip,bibref)
! if skip negative the parameter is not entered
          elseif(trim(pretag).eq.'</Bibliography>') then
!========================= end of nested Bibliography, nothing to do
             continue
!             write(*,*)'  End of Bibliography'
          else
!  still with matt=-1, nested tags without special action ...
!             write(*,667)trim(tagname),trim(pretag),matt,level,allappy, fline
! Probably end of file ...
             if(debug) write(*,*)'End of file at ',fline
             continue
667          format('   Unknown end of nested tag: "',a,'" "',a,'"',5i7)
          endif
!============================ prepare to read next tag, matt=-1
! We have to decrease the nesting level and prepare for reading a new tag
          level=level-1
          if(level.gt.0) then
! it must be wrong to set pretag same as tagname ??
             pretag='</'//trim(endoftag(level))//'>'
          else
! Level 0 means we have read the whole xtdb file, we should never be here
!             write(*,*)'We have read past end of </XTDB>'
             continue
          endif
          cycle readall
! Done all for end of nested tag, matt=-1
!============================ end of actions for all nested tags, matt=-1
       endif
! Here matt=0 or 1, maybe we need to save attributes
!       write(*,660)trim(tagname),trim(pretag),matt,fline
660    format('660 Tag: "',a,'" pretag "',a,'" ',i3,i7)
!
       if(matt.ge.0) then
          if(tagname(1:1).eq.' ') then
! This should not appear
! for rewind of masterfile it is OK
             if(allappy.eq.0) write(*,*)'xtdbtag found no tag on line ',&
                  fline,matt,allappy
             cycle readall
          endif
       else
! we have found the end of a nested tag
          write(*,*)'Handle end of nested tag: ',trim(tagname),matt,level,fline
       endif
! check if we have reached endoffile, maybe rewind and read again?
       if(fline.lt.0) then
! unless there are errors one should never come here
          write(*,*)'Closing file 3'
          goto 990
       endif
!======================================= select case for a new tag
! The tags and attribues are defined in gtp3_xml.F90
       lk=len_trim(tagname)+1
       tagno=0
       findtag: do kk=1,nxtdbtags
! compare ignoring trailing spaces
          if(tagname(1:lk).eq.xtdbtags(kk)(1:lk)) then
             tagno=kk;
! Abbreviated tags not allowed
             if(xtdbtags(kk).eq.' ') exit findtag
          endif
       enddo findtag
!       write(*,99)tagname(1:lk),tagno,fline,level,matt
99     format('99 Tag: ',a,', number',i3', line: ',i7,', level: ',2i3)
!
!
! detect tag
!------------------------------------------
       select case(tagno)
!------------------------------------------
       case default
! process the XTD tag, sometimes nested tags depend on previous attributes
          write(*,*)' *** Unknown tag ignored: <',trim(tagname),'>'
!------------------------------------------
       case(1) ! XTDB
          ip=1; values=' '
!          write(*,*)'XTDB: '
! the getatt function extracts attribute name and value sequentially from ip
          do while(getatt(attributes,ip,attname,values))
             if(debug) write(*,38)trim(attname),trim(values)
38           format(3x,'Att: ',a,' = ',a,i5)
          enddo
!------------------------------------------
       case(2) ! Defaults
          ip=1; values=' '
!          write(*,*)'Defaults: '
          do while(getatt(attributes,ip,attname,values))
! these should replace the software defaults such as LowT and HighT
!             write(*,38)trim(attname),trim(values)
             if(attname(1:8).eq.'Bibref') defaultbib=values
          enddo
! lowT and highT should be used for parameters and TPfuns
!------------------------------------------
       case(3) ! DatabaseInfo
          ip=1; values=' '
!          write(*,*)'DatabaseInfo: '
!          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
!          enddo
!------------------------------------------
       case(4) ! AppendXTDB
          ip=1; values=' '
!          write(*,*)'AppendTDB files: '
          do while(getatt(attributes,ip,attname,values))
! this will extract the attributes in any order
!             write(*,38)trim(attname),trim(values)
! save appfiles in appropriate variables (declared in gtp3_xml
! negative *appy means it is has to be done
             if(attname(1:6).eq.'Models') then
                modelappx=trim(values)
                modelappy=1
             endif
             if(attname(1:10).eq.'Parameters') then
                parappx=trim(values)
                parappy=1
             endif
             if(attname(1:6).eq.'TPfuns') then
                tpfunappx=trim(values)
                tpfunappy=1
             endif
             if(attname(1:12).eq.'Bibliography') then
                biblioappx=trim(values)
                biblioappy=1
             endif
! Not implemented.  Maybe useful for non-thermodynamic data?
!             if(attname(1:13).eq.'Miscellaneous') then
!                miscappx=trim(values)
!                miscappy=-1
!             endif
          enddo
! list AppendXTDB file set
!          write(*,*)'AppendXTDB: ',modelappy,parappy,tpfunappy,biblioappy
! if there is a Models appendix read that now!
          if(modelappy.eq.1) then
             call xtdbmodels(modelappx)
! setting modelappy=0 means we have read it
! The other AppendXTDB files opened after reading the primary file once
             modelappy=0
          endif
!------------------------------------------
       case(5) ! Element, 5 attributes. some may be missing ...
!          ip=1; values=' '; mass=1.0D0; h298=0.0D0; s298=0.0D0
!          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
!             if(attname(1:2).eq.'Id') tpfun=trim(values)
!------------------------------------------
       case(6) ! Species
          ip=1; values=' '
! we need to know the elements of the species
          mqmqa=' '; uniquac=' '
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
             if(attname(1:2).eq.'Id') spname=trim(values)
             if(attname(1:13).eq.'Stoichiometry') stoisp=trim(values)
             if(attname(1:5).eq.'MQMQA') mqmqa=trim(values)
             if(attname(1:7).eq.'UNIQUAC') uniquac=trim(values)
          enddo
!          write(*,*)'Calling OCenterspecies: ',trim(spname),'" "',&
!               trim(stoisp),'"'
! in OCenterspecies there is a check if species can be formed from elements
          call OCenterspecies(spname,stoisp,mqmqa,uniquac)
!------------------------------------------
       case(7) ! TPfun
          ip=1; values=' '
          if(matt.eq.0) then
! matt=0 means no nested tags, matt=-1 is taken care of earler
! We may have to read all TPfuns several times until all used has been found
! matt=0  means there are no nested tranges
             lowt=default_lowt
             hight=default_hight
             do while(getatt(attributes,ip,attname,values))
!                   write(*,38)trim(attname),trim(values)
                if(attname(1:2).eq.'Id') tpfun=trim(values)
                if(attname(1:4).eq.'LowT') lowt=trim(values)
                if(attname(1:4).eq.'Expr') expr=trim(values)
                if(attname(1:5).eq.'HighT') hight=trim(values)
             enddo
! semicolon is accepted and added if not present
             semicolon=len_trim(expr)
             if(expr(semicolon:semicolon).eq.';') then
                expr(semicolon:semicolon)=' '
                semicolon=semicolon-1
             endif
! terminate with N as in TDB files
             wholexpr=lowt//' '//expr(1:semicolon)//'; '//hight//' N'
! check if the TPfun is needed, tpfun is Id of <TPfun tag, TPfuns has no bibref
!             write(*,*)'3EY xtdbread TPfun: ',trim(tpfun)
             call xtdbentertpfun(tpfun,add)
             if(add) then
                call OCentertpfun(tpfun)
             endif
          else
! this TPfun has nested Trange tags, extract the current data and read more
             lowt=default_lowt
             hight=default_hight
             do while(getatt(attributes,ip,attname,values))
!                write(*,38)trim(attname),trim(values)
                if(attname(1:2).eq.'Id') tpfun=trim(values)
                if(attname(1:4).eq.'LowT') then
                   lowt=trim(values)
!                   write(*,*)' *** Explicit lowt: ',trim(tpfun),': ',lowt
                endif
                if(attname(1:4).eq.'Expr') expr=trim(values)
                if(attname(1:5).eq.'HighT') hight=trim(values)
             enddo
             semicolon=len_trim(expr)
             if(expr(semicolon:semicolon).eq.';') then
                expr(semicolon:semicolon)=' '
                semicolon=semicolon-1
             endif
!             if(attname(1:5).eq.'HighT') hight=trim(values)
! this wholexpr will be extended by Trange records
!             wholexpr=expr(1:semicolon)//'; '//trim(hight)//
! in the TDB format there is a Y after highT if there is expression after
! This Y is added when reading the Trange expression
             wholexpr=trim(lowt)//' '//expr(1:semicolon)//'; '//trim(hight)
!             write(*,*)'TPfun: ',trim(tpfun),' ',trim(wholexpr)
          endif
!------------------------------------------
       case(8) ! Trange  
          if(matt.gt.0) then
             write(*,*)'Trange tags cannot have nested tags',fline
             xtdberr=4900; goto 1000
          endif
          ip=len_trim(wholexpr)
!          write(*,88)fline,ip
88        format(/'In Trange',2i7)
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
             hight=default_hight
             if(attname(1:4).eq.'Expr') expr=trim(values)
             if(attname(1:5).eq.'HighT') hight=trim(values)
          enddo
          semicolon=len_trim(expr)
          if(expr(semicolon:semicolon).eq.';') then
             expr(semicolon:semicolon)=' '
             semicolon=semicolon-1
          endif
          if(attname(1:5).eq.'HighT') hight=trim(values)
! there is already an expression with a highT limit in wholexpr, join with Y
          wholexpr=trim(wholexpr)//' Y '//expr(1:semicolon)//'; '//trim(hight)
! len_trim of parid which is not assigned can be infinit!!
!          write(*,19)trim(wholexpr)
19        format('In Trange wholexpr: ',a)
!------------------------------------------
       case(9) ! Phase
          ip=1; values=' '
          if(allocated(phrec)) then
! remove data for any previous phase, this should already have been done
             write(*,*)'3Y Failed deallocate previous phase data!'
             deallocate(phrec)
          endif
          allocate(phrec)
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
             phrec%state='S'
             if(attname(1:2).eq.'Id') then
                phrec%id=trim(values)
             endif
             if(attname(1:13).eq.'Configuration') phrec%confent=trim(values)
             if(attname(1:5).eq.'State') phrec%state=trim(values)
! more data in phrec in nested tags <Sublattices <Constutents <AmendPhase ...
          enddo
!------------------------------------------
       case(10) ! Sublattices
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
! we need to know the number of sublattices to allocate phrec%clist
             if(attname(1:8).eq.'NumberOf') phrec%noof=trim(values)
             if(attname(1:14).eq.'Multiplicities') phrec%mult=trim(values)
          enddo
! use unformatted read to allocate pointer array for sublattice constituents
          read(phrec%noof,*)phnsub
          allocate(phrec%clist(phnsub))
 ! initiate phisub, incremented for each <Constituent tag
          phisub=0
!------------------------------------------
       case(11) ! Constituents
! there is one of these tag for each sublattice
          ip=1; values=' '
! phisub=0 set by sublattice tag, incremented for each Constituent record
          phisub=phisub+1
          if(phisub.gt.phnsub) then
             write(*,*)'Too many sublattices ',phisub,phnsub,', line: ',fline
             stop
          endif
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
! this is just sublattice number, will be ignored
             if(attname(1:10).eq.'Sublattice') &
                  phrec%clist(phisub)%subx=trim(values)
! This is list of constituents in the sublattice
! The character will be allocated automatically
             if(attname(1:4).eq.'List') phrec%clist(phisub)%list=trim(values)
          enddo
!          write(*,991)phisub,len_trim(phrec%clist(phisub)%list),&
!               phrec%clist(phisub)%list
991       format('Subl: ',i2,i5,', const: ',a)
!------------------------------------------
       case(12) ! CrystalStructure
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
             if(attname(1:9).eq.'Prototype') write(*,*)'Prottype: ',trim(values)
             if(attname(1:13).eq.'PearsonSymbol') &
                  write(*,*)'Pearson symbol: ',trim(values)
             if(attname(1:10).eq.'SpaceGroup') &
                  write(*,*)'Spacegroup: ',trim(values)
          enddo
          write(*,*)'Ignored in OC'
!------------------------------------------
       case(13) ! AmendPhase
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
             if(attname(1:6).eq.'Models') phrec%amendph=trim(values)
          enddo
!------------------------------------------
       case(14) ! Appendix
! Initial tag on AppendXTDB files 
          continue
!          write(*,*)'Appendix has no attributes'
!------------------------------------------
       case(15) ! DisorderedPart 15
!          write(*,*)'Found DisorderedPart line ',fline
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
! this is the name of the disordered phase if it is a separate phase
             if(attname(1:10).eq.'Disordered') phrec%dispar=trim(values)
             if(attname(1:3).eq.'Sum') sum=trim(values)
             if(attname(1:8).eq.'Subtract') sub=trim(values)
          enddo
!          write(*,*)'Finished extracting attributes "',sum,'"',allocated(sum)
          if(allocated(sum)) then
! informatted read of the number of sublattices to sum for ordered phase
             read(sum,*)kk
!             write(*,*)'Ordered sublattices: ',kk,phnsub
             if(kk.gt.phnsub .or. kk.lt.2) then
                write(*,*)kk,phnsub,fline
356             format('DisorderedPart sum is too large or small:',2i3,&
                     ' line ',i7)
                xtdberr=4444; goto 1000
             endif
          endif
          phrec%dispar=sum
! Fortran automatically extend allocated characters ... wow!
          if(allocated(disph)) then
             phrec%dispar=phrec%dispar//' '//disph
          endif
          if(allocated(sub)) then
             phrec%dispar=phrec%dispar//' '//sub
          endif
!------------------------------------------
       case(16) ! This was Disordered_3Part, now merged inside DisorderedPart
          write(*,*)'Not implemented tag  16 unused '
!------------------------------------------
       case(17) ! Parameter
! just for check, when first parameter entered list all entered phases
!          listphases=.false.
!          if(listphases) then
!             listphases=.false.
!             write(*,68)nselph
68         format(/'In xtdbread list of ',i5,' selected phases:')
! we have not created %const with the entered constituents!!!
!             do ip=1,nselph
!                write(*,69)trim(selph(ip)%phasename),trim(selph(ip)%const)
69              format('Phase: "',a,'" constituents: "',a,'"')
!             enddo
!          endif
          ip=1; values=' '
          lowt=default_lowt
          hight=default_hight
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
             if(attname(1:2).eq.'Id') parid=trim(values)
             if(attname(1:4).eq.'LowT') lowt=trim(values)
             if(attname(1:4).eq.'Expr') expr=trim(values)
             if(attname(1:5).eq.'HighT') hight=trim(values)
             if(attname(1:6).eq.'Bibref') bibref=trim(values)
          enddo
          semicolon=len_trim(expr)
          if(expr(semicolon:semicolon).eq.';') then
             expr(semicolon:semicolon)=' '
             semicolon=semicolon-1
          endif
          if(matt.eq.0) then
! save the parameter in OC, the ' N ' to be compatible with TDB files
!             wholexpr=lowt//' '//expr(1:semicolon)//'; '//hight//' N '//bibref
             wholexpr=lowt//' '//expr(1:semicolon)//'; '//hight//' N'
!             write(*,1700)parid,trim(wholexpr)
1700         format('Parameter ',a/a)
! save the parameter if phase and constituents enetered, skip not used
! If no biref set use the default
             if(.not.allocated(bibref)) bibref=' '
             call OCenterpar(parid,skip,bibref)
          else
! There is one or more Trange records
! IMPORTANT, bibref must be added when matt=-1 for this parameter!!!
             wholexpr=lowt//' '//expr(1:semicolon)//'; '//hight
!             write(*,*)'   <<<Nested parameter: ',bibref,matt,fline
          endif
!------------------------------------------
       case(18) ! Parameter2
          write(*,*)'Not implemented tag  18: Parameter2'
!------------------------------------------
       case(19) ! Bibliography
!          there are no attributes
!          ip=1; values=' '
!          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
!          enddo
          continue
!------------------------------------------
       case(20) ! Bibitem from masterfile.
! From an AppendXTDB the subroutine xtdbbiblio is used
! assume Id is always before text
          ip=1; values=' '
          if(getatt(attributes,ip,attname,values)) then
             if(attname(1:2).eq.'Id') bibref=trim(values)
             lk=len_trim(bibref)
!             write(*,*)'Bibitem: "',bibref(1:lk),'" line:',fline
! next attribute should be the reference text
             values=' '
             if(getatt(attributes,ip,attname,values)) then
! check if reference missing
                do kk=1,nselbib
                   if(bibref(1:lk).eq.selbib(kk)%bibitem(1:lk)) then
                      if(selbib(kk)%status.lt.0) then
                         selbib(kk)%status=1
! ignore attname, just take the next attribute value
                         selbib(kk)%data=trim(values)
                      endif
                   endif
                enddo
             else
! there is no text for this bibitem
                write(*,375)bibref(1:lk),nselbib,kk
375             format('*** Warning bibref "',a,'" has no text')
                selbib(kk)%data='No text'
             endif
          endif
!------------------------------------------
       case(21) ! Model maybe ModelInfo 
          write(*,*)'Not implemented tag  21 ModelInfo'
!------------------------------------------
! All model attributes read by subroutine xtdbmodels on a separate file
! Except the TernaryXpol
       case(22) ! Magnetism
          write(*,*)'Modelinfo magnetism  22'
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
             write(*,38)trim(attname),trim(values)
          enddo
!------------------------------------------
       case(23) ! Einstein 
          write(*,*)'Modelinfo Einstein  23'
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
             write(*,38)trim(attname),trim(values)
          enddo
!------------------------------------------
       case(24) ! Liquid2state
          write(*,*)'Modelinfo Liquid2State  24'
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
             write(*,38)trim(attname),trim(values)
          enddo
!------------------------------------------
       case(25) ! Volume
          write(*,*)'Modelinfo Volume  25'
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
             write(*,38)trim(attname),trim(values)
          enddo
!------------------------------------------
       case(26) ! EEC
          write(*,*)'Modelinfo EEC 26'
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
             write(*,38)trim(attname),trim(values)
          enddo
!------------------------------------------
       case(27) ! TernaryXpol
! This can be inside a phase tag ??
!          write(*,345)trim(attributes)
345       Format('Modelinfo TernaryXpol 27: '/a)
          ip=1; values=' '
          lph=0; xpoldata=' '
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values),ip
!---------------------------------------------------
             if(attname(1:5).eq.'Phase') then
                phname=trim(values)
                call capson(phname)
!                write(*,*)'3EY read all phases?',nselph
                if(nselph.gt.0) then
!                   write(*,*)'Calling findabbrphname 3: "',phname,'"'
                   call findabbrphname(phname,lph)
! if TernaryXpol is inside a <Phase tag it is not found ... it should be saved
! the phase may be entered later, check the constituents
!                   if(lph.le.0) cycle readall
! If inside a Phase tag the phases are not yet in selph ...
!                   write(*,*)'Found phase "',trim(values),'" "',&
!                        selph(lph)%phasename,'" ',lph
! we found the phase, lph is the phase index in selph if present in selph
!---------------------------------------------------
                endif
             elseif(attname(1:12).eq.'Constituents') then
! Check the species and if they are selected
!                write(*,*)'  Const: ',trim(values)
                call capson(values)
                lk=1
                is=0
                addon=' '
                all3: do while(is.lt.3)
! assuming a single space between species names
                   ir=index(values(lk:),' ')
                   if(ir.le.1) then
!                      write(*,*)'Ignoring TernaryXpol at line ',ir,fline
                      cycle readall
                   endif
                   spname=values(lk:lk+ir-1)
                   lk=lk+ir
!                   write(*,*)'Species "',trim(spname),'"',nselsp,is,lk
                   do kk=1,nselsp
                      spel=selsp(kk)%species
!                      write(*,*)'"',trim(spname),'" ? "',trim(spel),'"',kk,is
!                      write(*,*)'"',spname,'" ? "',spel,'"',kk
                      if(trim(spel).eq.trim(spname)) then
                         if(trim(selsp(kk)%species).eq.trim(spname)) then
                            is=is+1
                            addon=addon//' '//trim(spname)
!                            write(*,*)'Found: ',addon(3:),is,lk
                            cycle all3
!                         else
!                            write(*,*)'not same'
                         endif
                      endif
                   enddo
! this species it not selected, ignore the ternary extrapolation
!                   write(*,*)'Skip TernaryXpol as species "',&
!                        trim(spname),'" not selected',is
                   cycle readall
                enddo all3
! we arrive here if all 3 constituents are selected otherwise skipped
!---------------------------------------------------
             elseif(attname(1:4).eq.'Xpol') then
!                write(*,*)'Extracting Xpol: "',trim(values),'"'
                xpoldata=trim(values)
             else
                write(*,*)'Unknown attribute: "',attname,'" line',fline
             endif
          enddo ! end of all attributes
          if(xpoldata(1:1).eq.' ') then
             write(*,*)'3Y TernaryXpol error, no Xpol value, "',xpoldata,&
                  '" line',fline
             cycle readall
          endif
!---------------------------------------------------
! if we come here we can create a TernaryXpol record
! save it at the phase if lph is nonzero
          if(ternaryxwarning) then
             write(*,*)'Warning: no check on conflicting ternaryXpol tags'
             ternaryxwarning=.false.
          endif
          if(lph.gt.0) then
! we can directly add it to a selected phase
!             write(*,*)'3EY add TernaryXpol record to selected phase',lph
             if(associated(selph(lph)%terxpol)) then
                xpol=>selph(lph)%terxpol
                allocate(selph(lph)%terxpol)
                selph(lph)%terxpol%next=>xpol
                xpol=>selph(lph)%terxpol
             else
                allocate(selph(lph)%terxpol)
                nullify(selph(lph)%terxpol%next)
                xpol=>selph(lph)%terxpol
             endif
! no need to save the phase name, the data is linked from the phase
             xpol%phase=phname
             xpol%sps=addon(3:)
             xpol%xpol=xpoldata
          elseif(allocated(phrec)) then
! the phase is not (yet) selected (maybe the xpol tag is inside the phase tag)
             if(.not.associated(firstxpol)) then
! this is the first TernaryXpol without a selected phase
                allocate(firstxpol)
                nullify(firstxpol%next)
                xpol=>firstxpol
             else
! add the new Xpol record first
                allocate(xpol)
                xpol%next=>firstxpol
                firstxpol=>xpol
             endif
! enter data in xpol, we must use the phase name from current phase tag
             xpol%phase=trim(phrec%id)
!             write(*,*)'TernaryXpol inside phasetag: "',&
!                  xpol%phase,'" "',addon(3:),'" "',xpoldata,'"'
             xpol%sps=addon(3:)
             xpol%xpol=xpoldata
          else
             write(*,*)'TernaryXpol for ',trim(phname),' ignored, line',fline
             cycle readall
          endif
! if we arrive here the constituents are OK but we may not have a phase (lph=0)
!          write(*,*)'3EY TernaryXpol "',xpol%phase,'" ',lph
!------------------------------------------
       case(28) ! UnarySystem has no interest for the software
          write(*,*)'Found UnarySystems 28'
!------------------------------------------
       case(29) ! BinarySystem has no interest for the software
          write(*,*)'Found BinarySystem  29'
!------------------------------------------
       case(30) ! TernarySystem has no interest for the software
          write(*,*)'Found TernarySystem 30'  
       end select
!------------------------------------------
! end of file but maybe rewind and read for some special tag? such as TPfun
900    continue
!       write(*,*)'At label 900, select case: ',tagno
       cycle readall
!
    enddo readall
! errors and end of file-  ignoreEOT is true if just extracting elements
990 continue
!    write(*,*)'At label 990'
    if(.not. ignorEOT) then
! There can be TernaryXpol records not linked to a phase
       xpol=>firstxpol
!       nullify(lastxpol)
       fixpol: do while(associated(xpol))
! Either xpol will be linked to a selph(lph)%terxpol or igored
!          write(*,*)'3EY There is an unlinked TernaryXpol for: ',xpol%phase
!          write(*,*)'Calling findabbrphname 2: "',xpol%phase,'"'
          call findabbrphname(xpol%phase,lph)
          if(lph.le.0) then
             write(*,*)'Ignore TernaryXpol as ',xpol%phase,' not selected'
             if(.not.associated(lastxpol)) then
                xpol=>xpol%next
             endif
          else
! link the ternaryxpol to the phase
             write(*,*)'TernaryXpol included in phase ',&
                  trim(selph(lph)%phasename)
             lastxpol=>xpol%next
             if(associated(selph(lph)%terxpol)) then
                xpol%next=>selph(lph)%terxpol
             else
                nullify(xpol%next)
             endif
             selph(lph)%terxpol=>xpol
             xpol=>lastxpol
          endif
       enddo fixpol
!       if(debug) then
!          write(*,*)'3Y End of file',unit
       write(*,993)nselel,nselsp,nselph,nselpar,nseltp,nselbib
993       format(/'3EY nselel ',i4,' nselsp ',i4,' nselph ',i4,&
               ' nselpar ',i4,' nseltp ',i5,' nselbib ',i4/)
!       endif
       missingtp=0
       do ntp=1,nseltp
          if(seltpfun(ntp)%status.lt.0) missingtp=missingtp+1
       enddo
       prevmissbib=missingbib
       missingbib=0
       do ntp=1,nselbib
          if(selbib(ntp)%status.lt.0) missingbib=missingbib+1
       enddo
    endif
    if(debug .or. missingtp.gt.0 .or. missingbib.gt.0) then
       write(*,*)'Number of TPfuns and missing ones:',nseltp,missingtp
       write(*,*)'Number of bibref and missing ones:',nselbib,missingbib
    endif
!
999 continue
    close(unit)
1000 continue
! list or save all data in OC data structures
    if(.not.ignorEOT) then
       write(*,312)trim(masterfile),xtdberr
312 format(/'Listing of selected data from XTDB file: ',a,', error: ',i5)
       call list_selected_xtdbdata
    endif
    return
1100 continue
    write(*,*)'Error reading xtdbfile ',xtdberr
    goto 990
!
  end subroutine xtdbread

! \/!\/!\/!\/!\/!/!\/!\/!\/!\/!/!\/!\/!\/!\/!/!\/!\/!\/!\/!/!\/!\/!\/!\/!/!

  subroutine list_selected_xtdbdata
! temporary listing from xtdb local arrays
!
    implicit none
    character line*128, charge*24, spline*24,ch1*1
    integer jj,kk,ip,lk,nn
    type(octerxpol), pointer :: terxpol
!    
!    if(missingtp.gt.0) then
!       write(*,491)missingtp,missingbib
491    format('Missing TPfuns, bibitems: ',2i4)
!    endif
!
    write(*,12)
12  format(//'Finished reading the XTDB file, temporary listing')
! list all elements, species, phases with constit, parameters, tpfuns and biblio
! list all elements -------------------------------------------------
    write(*,311)
311 format(//'Elements:')
    write(*,10)(selel(nn)%elname,nn=1,nselel)
10  format('Elements entered: '20(a,1x))
!
! list all species --------------------------------------------------
    write(*,299)
299 format(/'Species:')
    if(allocated(selsp)) then
       do ntp=1,nselsp
          jj=selspord(ntp)
          charge=' '
          if(jj.le.0) write(*,*)'3EY selspord: ',ntp,jj,selsp(ntp)%species
          if(selsp(jj)%charge.ne.0.0d0) &
               write(charge,'(F10.6)')selsp(jj)%charge
          spline=selsp(jj)%species
          if(allocated(selsp(jj)%elnames)) then
             nn=size(selsp(jj)%elnames)
             lk=27
             do kk=1,nn
                write(line(lk:),23)selsp(jj)%elnames(kk),selsp(jj)%stoicc(kk)
23              format(a,2x,F10.7)
                lk=lk+15
             enddo
          if(selsp(jj)%charge.ne.0.0D0)spline=trim(spline)//' charge '//charge
          else
             lk=len_trim(line)
             line(lk+3:)='MQMQA quad'
          endif
          write(*,24)ntp,trim(spline)
24        format('Species: ',i3,2x,a)
       enddo
    endif
!
! list all selected phases --------------------------------------------
    write(*,*)'Press return to list selected phases'
    read(*,'(a)')ch1
    write(*,298)
298 format(/'Phases:')
    do ip=1,nselph
       write(*,69)trim(selph(ip)%phasename),selph(ip)%confent,&
            selph(ip)%nsublat,trim(selph(ip)%mult),&
            trim(selph(ip)%const)
69     format('Phase: ',a,' cfg: ',a,'  Subl: ',i2,'  Mult: ',a/&
            '   Constituents: "',a,'"')
       if(allocated(selph(ip)%amendph)) write(*,76)selph(ip)%amendph
76     format('   AmendPhase: ',a)
       if(allocated(selph(ip)%dispar)) &
            write(*,'("   DisPart ",a)')selph(ip)%dispar
       terxpol=>selph(ip)%terxpol
!       if(.not.associated(terxpol)) write(*,*)'No ternary methods'
       do while(associated(terxpol))
          write(*,78)trim(selph(ip)%phasename),terxpol%sps,terxpol%xpol
78        format('   TernaryXpol Phase="',a,&
               '" Constituents="',a,'" Xpol="',a,'"')
          terxpol=>terxpol%next
       enddo
    enddo
!
! list all parameters selected ---------------------------------------
    write(*,*)'Press return to list selected parameters'
    read(*,'(a)')ch1
    if(nselpar.gt.0) then
       write(*,63)nselpar
63     format(/'All ',i5,' parameters entered')
       do ip=1,nselpar
          if(len_trim(selpar(ip)%parname)+len_trim(selpar(ip)%data).gt.60) then
             write(*,72)ip,trim(selpar(ip)%parname),trim(selpar(ip)%data)
72           format(i4,2x,a,' ='/10x,a)
          else
             write(*,73)ip,trim(selpar(ip)%parname),trim(selpar(ip)%data)
73           format(i4,2x,a,' = ',a)
          endif
       enddo
    endif
!
! list asll TPfuns ---------------------------------------------------
    write(*,*)'Press return to list selected TPfuns'
    read(*,'(a)')ch1
    if(nseltp.gt.0) then
       write(*,82)nseltp
82     format(/'All ',i5,' TPfun entered or missing (-1)')
       do ip=1,nseltp
          if(len_trim(seltpfun(ip)%data).le.50) then
             write(*,77)ip,seltpfun(ip)%tpfunname,seltpfun(ip)%status,&
                  trim(seltpfun(ip)%data)
77           format(i3,2x,a,2x,i2,3x,a)
          else
             write(*,87)ip,seltpfun(ip)%tpfunname,seltpfun(ip)%status,&
                  trim(seltpfun(ip)%data)
87           format(i3,2x,a,2x,i2/5x,a)
          endif
       enddo
    endif
!
! list all biteams ---------------------------------------------------
    write(*,*)'Press return to list selected bibilography'
    read(*,'(a)')ch1
    if(nselbib.gt.0) then
       write(*,772)nselbib
772    format(/'All ',i4,' bibitems entered or missing (-1)')
       do ip=1,nselbib
          write(*,79)ip,selbib(ip)%bibitem,selbib(ip)%status,selbib(ip)%data
79        format('Bibitem ',i4,' "',a,'" ',i3/5x,a)
       enddo
    endif
!
    return
  end subroutine list_selected_xtdbdata
    
! \/!\/!\/!\/!\/!/!\/!\/!\/!\/!/!\/!\/!\/!\/!/!\/!\/!\/!\/!/!\/!\/!\/!\/!/!

  subroutine xtdbtag(unit,fline,tagname,matt,pretag,attributes)
! this subroutine extract a tag and its attributes from lines read from file.  
! Tag begins with "<tagname>" if no attributes, otherwise "<tagname "
! ONLY ONE TAG PER LINE but the attributes can be on following lines.
! End of tag attributes is ">" or "/>", the latter also means end of tag
! End of tag may be on a separate line as "/>" or if nested </tagname>
! USING XEOLCH
    
! Attributes are one or more identifier="values" decoded in calling routine
!-------------------
!    character tagname*(*),pretag*(*),attributes*(*)
    character(len=:), allocatable :: attributes
!    character attributes*(1024)
    character tagname*(*),pretag*(*)
    integer unit,fline,matt
!-------------
    character line*512
    integer ep,ip,jp,kp,tp,taglen,tagend,lines,lentagname
    logical comment,rewindtp
! problem initiating rewindtp? It should be false unless rewinding
! and the toggle rewind within TPfun tags
    save rewindtp
!-------
    ! matt  tagend meaning
!  0    -2     looking for <tag
!  0    -1     end of attributes not found 
!  1     0     end of attribues found but not end of tag, nesting
! -1     0     end of nested tag found, decrease level
    tagend=-2
    comment=.false.
    attpos=-1            ! attpost incremeneted by 2 when used
    tagname=' '
    attributes=' '
    lines=0
! allappy is 0 when reading primary XTDB file first time.
! After rewinding it is toggled false/true when reading nested TPfun
    if(allappy.le.1 .and. rewindtp) then
       rewindtp=.false.
!       write(*,*)'Setting rewindtp .FALSE.',allappy,rewindtp,fline
    endif
!    if(pretag(1:1).ne.' ') write(*,*)'xtdbtag with pretag: ',trim(pretag)
!================================================
! we may have to read the file several times to pick up all tags ....
    readtag: do while(.true.)
! maybe read several lines until all attributes extracted
       read(unit,100,end=1100)line
100    format(a)
       fline=fline+1
       lines=lines+1
!       write(*,102)trim(pretag),fline,trim(line)
102    format(/'xtdbtag mode "',a,'" line ',i6,' "',a,'"')
! skip continuation of comment lines
       if(len_trim(line).gt.500) then
! warning of line longer than 500 characters
          write(*,*)' *** Warning, XML tag longer than 500 characters.',fline
       endif
       if(comment) then
          ip=index(line,'--')
          if(ip.gt.0) then
             if(line(ip+2:ip+2).ne.'>') then
                write(*,*)' *** Error, use of -- inside comment, line',fline
                xtdberr=4518
                exit readtag
             endif
! we found end of multiline comment -->
!             write(*,*)'End of multiline comment',fline
             attributes=' '
             attpos=-1
             comment=.false.
! skip rest of line
             if(len_trim(line).gt.ip+3) then
                write(*,103)fline
103             format('Skipping text trailing comment on line ',i7)
             endif
! reset line count for next tag
             tagend=-2
             lines=0
          endif
          cycle readtag !------- skip rest of line and cycle
       endif
!================================================
! use xeolch to find first character on line
       ip=1
       if(xeolch(line,ip)) then
! if line empty read next line
          cycle readtag
       endif
! first character of a tag must be <
       if(line(ip:ip).eq.'<') then
          if(tagend.ne.-2) then
             write(*,*)' *** ERROR XTDB has only one tag per line,',fline
             xtdberr=4515; exit readtag
          endif
          if(line(ip:ip+3).eq.'<!--') then
! we have found start of a comment, if not finish on same line set comment
             if(.not.index(line(ip:),'-->').gt.0) comment=.TRUE.
             cycle readtag
          endif
          if(line(ip:ip+1).eq.'</') then
! we have </ it is the end of a nested tag and should be equal to pretag
             jp=len_trim(pretag)
             if(line(ip:ip+jp-1).eq.pretag(1:jp)) then
! end of nested tag, skipping any text after </tagname>, negative matt
!                write(*,666)line(ip:ip+jp-1),trim(pretag),fline
666             format('PRETAG: "',a,'" "',a,'" ',i7)
                matt=-1
                if(allappy.gt.1 .and. line(ip:ip+7).eq.'</TPfun>') then
! A TPfun with nested Trange, return as normal after clearing rewindtp
                   rewindtp=.FALSE.
                endif
             elseif(allappy.gt.1) then
! we are rewinding file loocking only for <TPfun and </TPfun>
!                write(*,*)'Ignore all endoftags except TPfun'
                if(line(ip:ip+7).ne.'</TPfun>') then
                   cycle readtag
                endif
             elseif(allappy.le.1) then
                write(*,*)' *** Error, illegal end of tag, line',fline
                xtdberr=4514
             endif
! just end of reading, maybe indicate an error?
             if(.not.(ignorEOT .or. (trim(line(ip:)).eq.trim(pretag)))) then
                if(allappy.eq.0) write(*,345)trim(line(ip:)),trim(pretag),fline
345             format('Found unexpected end of tag: "',a,'" "',a,'" ',i3,i7)
             endif
             exit readtag
          endif
! We have found a new tag.
! we can be after rewinding the master file. If so allappy=3
! and only TPfuns or bibitem are accepted,  ignore all other tags 
!          write(*,*)'Check rewinding: 1: ',allappy,rewindtp,trim(line)
          if(allappy.gt.1) then
!             write(*,777)allappy,rewindtp,line(ip:ip+25),ip,matt,fline
777          format('Allappy ',i3,l3,' line "',a,'" ',3i5)
             if(line(ip:ip+7).eq.'<Bibitem') then
                tagname='Bibitem'
                attributes=trim(line(ip+8:))
!                attpos=len_trim(attributes)
                matt=0
! no nested tags
                exit readtag
             elseif(line(ip:ip+6).eq.'<TPfun ') then
! Save this to allow for a Trange tag
! extract attributes and return as normal except set rewindtp
                if(matt.eq.0) then
                   rewindtp=.TRUE.
!                   write(*,*)'Setting rewindtp .TRUE.',fline
                endif
             else
                if(line(ip:ip+7).eq.'<Trange ') then
!                   if(allappy.gt.0) write(*,*)'Preatag: ',trim(pretag)
                   continue
                else
! ignore this tag
!                   write(*,39)line(ip:ip+8),rewindtp,allappy,fline
39                 format('Rewind, ignoring this tags: "',a,'" ',l3,3i7)
                   cycle readtag
                endif
             endif
          endif
!=========================== below is treatment at first read
! new tag found, check if it has attributes, find first space
!          write(*,*)'Check rewinding 2: ',allappy,rewindtp,trim(line)
          jp=index(line(ip:),' ')
!          write(*,*)'Found a space: ',line(ip+1:ip+jp-1),jp !-------------
          if(line(ip+jp-2:ip+jp-2).eq.'>') then
! Tag name end with > this is a tag without attributes but with nested tags
! skip any trailing text
             tagname=line(ip+1:ip+jp-3)
!             write(*,*)'TAGNAME1: "',trim(tagname),'"'
             matt=1
             exit readtag
          else
             tagname=line(ip+1:ip+jp-2)
!             write(*,*)'TAGNAME2: "',trim(tagname),'"'
          endif
! we have found the tagname, tagend=-1 indicate attributes can be on next linws
          lentagname=len_trim(tagname)
          tagend=-1
          ip=ip+jp
       endif
! we are looking for attributes in line after ip
       if(xeolch(line,ip)) then
! line empty and we have not found end of attributes
          cycle readtag
       else
          jp=index(line(ip:),'>')
! this indicate end of attributes
          if(jp.le.0) then
! no end of attributes, copy all to attributes and read next line
             if(ip.gt.0) then
!                attributes(attpos+2:)=trim(line(ip:))
                attributes=attributes//' '//trim(line(ip:))
                attpos=len_trim(attributes)
             endif
          else
! Found > as end of attributes.  maybe end of tag, matt=0 means no nested tags
             if(line(ip+jp-2:ip+jp-1).eq.'/>') then
!                attributes(attpos+2:)=line(ip:ip+jp-3)
                attributes=attributes//' '//line(ip:ip+jp-3)
                matt=0
             else
! If no end of tag nested tags may follow on next lines
!                attributes(attpos+2:)=line(ip:ip+jp-2)
                attributes=attributes//' '//line(ip:ip+jp-2)
                tp=ip+jp+1
                matt=1
! check for rubbish and maybe full endoftag after rubbish
                if(.not.xeolch(line,tp)) then
                   jp=index(line(tp:),'<')
                   if(jp.gt.0) then
! The text has a <, if part of </tagname> means end of tag
             write(*,69)line(tp+jp+1:tp+jp+lentagname),trim(tagname),tp
69                    format('3EY Are "',a,'" and "',a,'" equal?',i7)
                      if(line(tp+jp:tp+jp).eq.'/' .and. &
             line(tp+jp+1:tp+jp+lentagname).eq.trim(tagname)) then
! end of tag found after >, no nesting for this tag
                write(*,*)'3EY Trash found between > and tagend on line ',fline
                         matt=0
                      else
                       write(*,*)'3EY *** Error, new tag after > on line ',fline
                         xtdberr=4766
                         exit readtag
                      endif
                   else
! ignore rubbish after >
                      write(*,88)fline,trim(line)
88                    format('3EY *** Rubbish after > line ',i7,&
                           '. Quit reading.'/a/)
                      xtdberr=4766
                      exit readtag
                   endif
                endif
             endif
! we have found > or /> or </tagname>
             exit readtag
          endif
       endif
! maybe add
       if(lines.gt.2) then
! check for tags with multiple lines
          if(tagname(1:1).eq.' ') write(*,104)fline-1
104       format('3EY *** Warning, line without tag, line ',i7)
          if(lines.gt.3 .and. matt.lt.0) write(*,106)fline-3
106       format('3EY *** Warning, very long list of attributes, line ',i7)
       endif
    enddo readtag
! puh.................
1000 continue
!    write(*,1002)matt,fline,trim(tagname)
1002 format('Exiting xtdbtag: ',i2,i7,' ',a)
    return
1100 continue
!    write(*,*)'End of file in xtdbtag'
    xtdberr=4700
    fline=-1
    goto 1000
  end subroutine xtdbtag

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  logical function xeolch(line,ip)
! identical to eolch in METLIB, except it can handle endofline
!...End Of Line CHeck, TO SKIP SPACES/TAB FROM IP. RETURNS .TRUE. IF another
    character line*(*)
    integer ip
!
    integer, parameter :: itab=9
    xeolch=.true.
    if(ip.le.0) ip=1
    loop: do while(ip.le.len(line))
       if(line(ip:ip).eq.' ' .or. ichar(line(ip:ip)).eq.itab) then
          ip=ip+1
       else
          exit loop
       endif
    enddo loop
! with allocated characters there is maybe no space before EOL
! only when ip is larger than len(line) it returns true
    if(ip.le.len(line)) xeolch=.false.
900 RETURN
  END function xeolch

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  logical function getatt(text,ip,attname,value)
! extract "values" of any XML "attname" in text from position ip
    character*(*) text,attname,value
    integer ip
!
    integer jp,kp,attlen
!
    getatt=.false.
!    write(*,*)'3EY getatt 1',len(text)
    find: if(.not.xeolch(text,ip)) then
       if(text(ip:ip).eq.'=') exit find
       jp=index(text(ip+1:),'=')
       if(jp.gt.0) then
! the attname is terminated by a = (possibly preceeded by spaces)
          attname=text(ip:ip+jp-1)
!          write(*,*)'3EY getatt 2 ',trim(attname)
          ip=ip+jp+1
          if(.not.xeolch(text,ip)) then
! the values are preceeded by a " (possibly prceeded by spaces)
             if(text(ip:ip).eq.'"') then
                jp=index(text(ip+1:),'"')
                if(jp.gt.0) then
                   value=text(ip+1:ip+jp-1)
                   ip=ip+jp+1
                   getatt=.true.
                else
                   xtdberr=4601
                endif
             else
                xtdberr=4602
             endif
          else
             xtdberr=4603
          endif
       else
          xtdberr=4607
       endif
! no error if line empty
    endif find
    return
  end function getatt

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  logical function check_mpid(mpid,phase)
! check that phase has a model corresponding to the mpid of a parameter
    character*(*) mpid,phase
!
    integer np,ip
! 1. loop to find the phase
! 2. loop models for the phase to find one with the MPID
! 3. return TRUE if found, FALSE if not    
    check_mpid=.FALSE.
1000 continue
    return
  end function check_mpid


!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine getxtdbatt(attname,text,ip,values)
! extract "values" of attribute "att" of "tag" from text at position ip
    character*(*) text,attname,values
    integer att,ip
!
    integer jp,kp,attlen
! seach for attribute att is the array with attributes
    attlen=len_trim(attname)
!    write(*,*)'3EY getxtdbatt 1 >',text(ip:ip+25),'<',ip,attlen
    jp=index(text(ip:),attname(1:attlen))
    if(jp.le.0) goto 1100
! set jp to position after attname, the attribute must finish with ="
    jp=ip+jp+attlen-1
!    write(*,*)'3EY getxtdbatt 2 >',text(jp:jp+5),'<',jp
    if(text(jp:jp+1).ne.'="') goto 1110
    kp=index(text(jp+2:),'"')
!    write(*,*)'3EY getxtdbatt 3 >',text(jp:kp+5),'<',jp,kp
    if(kp.le.0) goto 1120
    values=text(jp+2:jp+kp)
! update ip to position after "
    ip=jp+kp+1
!
    1000 continue
!    write(*,*)'3EY exit getxtdbatt ',text(ip:ip+5),ip,xtdberr
    return
! cannot find attribute
1100 xtdberr=4501
    goto 1000
! attribute has no trailing ="
1110 xtdberr=4502
    goto 1000
! attribute value has no final "
1120 xtdberr=4503
    goto 1000
  end subroutine getxtdbatt

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!
! xtdbinitmpid removed
!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine xtdbmodels(appfile)
! reads an AppendXTDB file with models and MPID
! may change a default MPID
    character*(*) appfile
!
!    character(len=:), allocatable :: attributes
!    character attributes*(1024)
    character(len=:), allocatable :: attributes
    character tagname*18,pretag*24,values*24,attname*24
    integer unit,matt,lc,lt,ip,sm,mpid,level
    logical modeltag
!    
    unit=26
    write(*,5)trim(appfile)
5   format('In xtdbmodels reading: ',a)
!
!    xtdbmpid is allocated in gtp3EX with 9 models
!
    fline=0
    pretag=' '
    open(unit,file=appfile,access='sequential',form='formatted',&
         err=1200,status='old')
!
! initiating xtdbmpid
    if(.not.allocated(xtdbmodel)) then
! maybe deallocate and allocate new?
       write(*,*)'Allocating model MPID'
       call xtdbinitmpid(nxtdbmpids)
    else
       write(*,*)'3EY using current model MPID'
    endif
!
!    write(*,*)'Model file opened'
    level=0
    modeltag=.false.
    pretag='</Models>'
!
    models: do while(.true.)
       if(fline.lt.0) exit models
!       write(*,*)'Call with pretag "',trim(pretag),'"'
       call xtdbtag(unit,fline,tagname,matt,pretag,attributes)
       if(xtdberr.ne.0) goto 1100
!       write(*,10)trim(tagname),trim(attributes),fline,matt,modeltag
10     format('Model tag: "',a,'" att "',a,'"',i7,i3,l2)
       lt=len_trim(tagname)+1
! make sure we prepare an endoftag
! This should use levels but here we have only Models and Bibliography
       if(matt.eq.1) then
          pretag='</'//tagname(1:lt-1)//'>'
          level=level+1
       elseif(matt.eq.-1) then
          level=level-1
!          write(*,*)'Leveling down: ',trim(tagname),' ',trim(pretag),level
          pretag='</Models>'
! at level 0 we close the file
          if(level.eq.0) goto 1000
       endif
       if(.not.modeltag) then
          if(tagname(1:lt).eq.'Models ') then
             modeltag=.true.
          else
!             write(*,*)'Expecting only nested modeltags'
             xtdberr=4700; goto 1100
          endif
          cycle models
       endif
       if(matt.lt.0) then
! just skip the end of a model tag
          cycle models
       endif
!
!       write(*,*)'Model tag "',tagname(1:lt),'" pretag "',trim(pretag),'"'
! tags expected are
       if(tagname(1:lt).eq.'Magnetic ') then
!------------------------------------------!       case(22) ! Magnetism
!          write(*,40)trim(tagname),trim(attributes)
40        format('*** Addition: ',a,' "',a,'"')
          ip=1; values=' '; sm=0; mpid=0
          do while(getatt(attributes,ip,attname,values))
             lc=len_trim(values)
!             write(*,38)trim(attname),trim(values),lc,ip
38           format(3x,'Att: "',a,'" = "',a,'" ',i3,i7)
             if(sm.eq.0) then
! first attribute must be modelid
                if(attname(1:2).eq.'Id') then
                   do sm=1,3
                      if(values(1:lc).eq.xtdbmodel(sm)%modelid) exit
                   enddo
                endif
                if(sm.le.0 .or. sm.gt.3) then
!                   write(*,*)'3EY Unknown magnetic model "',trim(values),'"'
                   xtdberr=7777; goto 1000
                endif
             elseif(sm.eq.1) then
! Hmmmm in this way the database can use different symbols for each phase ???
! IHJBCC: second and third attribiures are MPD1 and MPD2
                mpid=mpid+1
!                write(*,*)'Storing MPID in xtdbmodel',sm,mpid
! accept the name of the mpid in the XTDB file
                xtdbmodel(sm)%mpid(mpid)=trim(values)
!                write(*,63)sm,mpid,xtdbmodel(sm)%mpid(mpid)
63              format('Model ',i1,', mpid',i1,' in xtdb file is: "',a,'"')
             elseif(sm.eq.2) then
! IHJREST: second and third attribiures are MPD1 and MPD2
                mpid=mpid+1
                xtdbmodel(sm)%mpid(mpid)=trim(values)
!                write(*,63)sm,mpid,xtdbmodel(sm)%mpid(mpid)
             elseif(sm.eq.3) then
! IHJQX: second, third and forth attribiure, MPD1, MPD2 and MPD3
                mpid=mpid+1
                xtdbmodel(sm)%mpid(mpid)=trim(values)
!                write(*,63)sm,mpid,trim(xtdbmodel(sm)%mpid(mpid))
             endif
! skip bibitem
             if((mpid.eq.2 .and. sm.lt.3) .or. mpid.eq.3) exit
          enddo
          if(sm.eq.3) then
!             write(*,67)tagname(1:lt),sm,trim(xtdbmodel(sm)%modelid),&
!                  trim(xtdbmodel(sm)%mpid(1)),trim(xtdbmodel(sm)%mpid(2)),&
!                  trim(xtdbmodel(sm)%mpid(3))
67           format('Model tag: ',a,i3,' Id: "',a,'" MPIDs: ',3('"',a,'" '))
          else
!             write(*,66)tagname(1:lt),sm,trim(xtdbmodel(sm)%modelid),&
!                  trim(xtdbmodel(sm)%mpid(1)),trim(xtdbmodel(sm)%mpid(2))
66           format('Model tag: ',a,i3,' Id: "',a,'" MPIDs: "',a,'" "',a,'"')
          endif
       elseif(tagname(1:lt).eq.'Einstein ') then
!------------------------------------------
!       case(23) ! Einstein 
!          write(*,*)'Modelinfo Einstein  23'
          ip=1; values=' '; sm=4
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
             if(attname(1:4).eq.'MPID') xtdbmodel(sm)%mpid(1)=trim(values)
          enddo
!          write(*,65)tagname(1:lt),sm,trim(xtdbmodel(sm)%modelid),&
!               trim(xtdbmodel(sm)%mpid(1))
65           format('Model tag: ',a,i3,' Id: "',a,'" MPIDs: "',a,'"')
       elseif(tagname(1:lt).eq.'Liq2State ') then
!------------------------------------------
!       case(24) ! Liq2state
!          write(*,*)'Modelinfo Liq2State  24'
          ip=1; values=' '; sm=5
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
             if(attname(1:5).eq.'MPID1') xtdbmodel(sm)%mpid(1)=trim(values)
             if(attname(1:5).eq.'MPID2') xtdbmodel(sm)%mpid(2)=trim(values)
          enddo
!          write(*,66)tagname(1:lt),sm,trim(xtdbmodel(sm)%modelid),&
!               trim(xtdbmodel(sm)%mpid(1)),trim(xtdbmodel(sm)%mpid(2))
       elseif(tagname(1:lt).eq.'Volume ') then
!------------------------------------------
!       case(25) ! Volume, not implemented
!          write(*,*)'Modelinfo Volume  25'
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
          enddo
       elseif(tagname(1:lt).eq.'EEC ') then
!------------------------------------------
!       case(26) ! EEC
!          write(*,*)'Modelinfo EEC 26'
          ip=1; values=' '
          do while(getatt(attributes,ip,attname,values))
!             write(*,38)trim(attname),trim(values)
          enddo
       elseif(tagname(1:lt-1).eq.'Bibliography') then
!------------------------------------------
!          write(*,*)'  Ignoring bibliography of models, closing file'
          goto 1000
       else
!------------------------------------------
! Models which has no MPID or is otherwise without model tag
          if(tagname(1:lt).eq.'DisorderedPart ') then
             ip=1; values=' '
             do while(getatt(attributes,ip,attname,values))
!                write(*,38)trim(attname),trim(values)
             enddo
!             write(*,*)'  The DisorderedPart tag is used in the Phase tag'
          elseif(tagname(1:lt).eq.'Permutations ') then
             ip=1; values=' '
             do while(getatt(attributes,ip,attname,values))
!                write(*,38)trim(attname),trim(values)
             enddo
!            write(*,*)'  The permutations Id is in AppendPhase Model attribute'
          elseif(tagname(1:lt).eq.'TernaryXpol ') then
             ip=1; values=' '
             do while(getatt(attributes,ip,attname,values))
!                write(*,38)trim(attname),trim(values)
             enddo
!             write(*,*)'  TernaryXpol tags are specified for each ternary'
          elseif(tagname(1:lt).eq.'EBEF ') then
!             write(*,*)'  EBEF is indicated by the parameters of the phase'
             continue
          else
!             write(*,*)'Model "',trim(tagname),'" not kown by software'
             continue
          endif
       endif
!
    enddo models
!-----------------
1000 continue
    write(*,1005)trim(appfile),fline
1005 format('Closing model file: ',a,i7/)
    close(unit)
1010 continue
    return
1100 continue
    write(*,*)'Error ',xtdberr,' reset after reading ',trim(appfile)
    xtdberr=0
    goto 1000
1200 continue
    write(*,*)'Error opening ',trim(appfile),' default MPID used'
    goto 1010
!
  end subroutine xtdbmodels

 !\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine xtdbbiblio(appfile)
! Listing bibitems from appfile with bibitems inside a bibliography tag
!    integer nselbib,maxbib
    character*(*) appfile
    character tagname*18,pretag*24
    character(len=:), allocatable :: attributes
!    character*1024 attributes
    character attname*24,values*128
!    character(len=:), allocatable :: attname,values
    integer unit,lc,ip,matt,nbib,lbib,found
!    integer, dimension(:), allocatable :: notfound
    logical bibtag
!    
    unit=25
    fline=0
    ignorEOT=.true.
    if(debug) write(*,5)nselbib,trim(appfile)
5   format(/'Trying to extract ',i5,' bibliographic references from: ',a,i7)
    open(unit,file=appfile,access='sequential',form='formatted',status='old')
!
!    bibtag=.false.
    bibtag=.true.
    matt=0
    found=0
    bibitems: do while(.true.)
       call xtdbtag(unit,fline,tagname,matt,pretag,attributes)
       if(xtdberr.ne.0) goto 1100
!*       write(*,10)trim(tagname),fline,matt,bibtag
10     format('Tag: "',a,'" ',i7,i3,l2)
       lc=len_trim(tagname)+1
       if(.not.bibtag) then
! with allappy.gt.1 the bibliography will never be returned from xtdbtag!!
          if(tagname(1:lc).eq.'Bibliography ') then
             pretag='</Bibliography>'
             bibtag=.true.
          else
             cycle bibitems
          endif
       else
          if(matt.eq.0) then
             ip=1; values=' '
             if(getatt(attributes,ip,attname,values)) then
! extract the bibitem Id, no trailing spaces allowed
!*                write(*,544)trim(attname),trim(values)
544             format('Read from file: "',a,'" "',a,'"')
! Assume first attribute is "Id"
!                if(attname(1:2).eq.'Id') then
!                write(*,*)' *** bibitem: ',trim(values)
                find: do nbib=1,nselbib
!                   lbib=len_trim(selbib(nbib)%bibitem)
!                   write(*,555)values(1:lbib),trim(selbib(nbib)%bibitem),nbib
555                format('Compare "',a,'" "',a,'"',i7)
!                   if(values(1:lbib).eq.selbib(nbib)%bibitem) then
                   if(values.eq.selbib(nbib)%bibitem) then
! This Id is used, extract the reference 
!                      write(*,556)values(1:lbib),trim(selbib(nbib)%bibitem)
556                   format('Found: "',a,'" "',a,'"')
                      selbib(nbib)%status=1
                      found=found+1
! We found the bibitem, extract the reference as next value, ignore attname
                      if(getatt(attributes,ip,attname,values)) then
!                         write(*,50)selbib(nbib)%bibitem,values
50                       format('Reference: ',a,' is ',a) 
                         selbib(nbib)%data=trim(values)
                         exit find
                      else
                         write(*,60)trim(appfile),fline
60                       format('3EY Formatting error file: "',a,'" line ',i7)
                      endif
                   endif
                enddo find
             endif
          endif
       endif
! make sure we prepare an endoftag
    enddo bibitems
1000 continue
    if(debug) write(*,*)'Found ',found,' relevant bibliographic references'
    close(unit)
    return
! errors model
1100 continue
    goto 1000
  end subroutine xtdbbiblio

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine addmissingtp(origin)
! extract TPfuns called inside wholexpr and add those missing to seltpfun
! with seltpfun(*)%status -1
    character*(*) origin
! all variables global
!
    integer ip,jp,kp,ntp
    character symbol*16
    ip=1
!    write(*,*)'Checking for new tpfuns'
    big: do while(ip.lt.len_trim(wholexpr))
! istpfun extracts symbols inside wholexpr ip is updated inside istpfun
       call istpfun(wholexpr,ip,symbol)
       if(symbol(1:1).ne.' ') then
          do ntp=1,nseltp
!             if(symbol.eq.alltpfun(ntp)) then
             if(symbol.eq.seltpfun(ntp)%tpfunname) then
                symbol=' ';cycle big
             endif
!             write(*,55)ntp,ip,symbol,seltpfun(ntp)%tpfunname
55           format('seltpfun ',2i5,' "',a,'" and "',a,'"')
          enddo
! this symbol is missing
!          write(*,*)' >>> addmissingtp tpfun: ',trim(symbol)
          nseltp=nseltp+1
!          write(*,60)trim(origin),trim(symbol),nseltp
60        format(' >>> TPfun used by "',a,'" added "',a,'" ',i4)
          seltpfun(nseltp)%tpfunname=symbol
          seltpfun(nseltp)%status=-1
       endif
    enddo big
1000 continue
    return
  end subroutine addmissingtp

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine istpfun(line,ip,symbol)
! extract unknown symbols (TPfuns) from an expression after position ip
    character line*(*),symbol*(*)
    integer ip
!------
! Strange error here when everyting was accepted as TPfuns ....
    integer jp,kp,mp
    character ch1*1
    call capson(line)
    kp=0
!    write(*,*)'3EY Looking for tpfun in: "',trim(line),'"',ip,kp
    symbol=' '
!    write(*,*)'In istpfun',ip,trim(line),len(line)
    issymbol: do while(ip.lt.len_trim(line))
! symbols must start with letter A-Z and can contain letters, digits and "_"
       ch1=line(ip:ip)
       ip=ip+1
!       write(*,*)'Testing "',ch1,'": symbol: "',trim(symbol),'" ',ip,kp
       if(ch1.ge.'A' .and. ch1.le.'Z') then
! first character must be a letter
          kp=kp+1
          symbol(kp:kp)=ch1
!          write(*,*)'Acceping "',ch1,'": symbol: "',trim(symbol),'" ',ip,kp
          cycle issymbol
       elseif(kp.ge.1 .and. &
            ((ch1.ge.'0' .and. ch1.le.'9') .or. ch1.eq.'_')) then
! 2nd and later character can be number or "_"
          kp=kp+1
          if(kp.gt.len(symbol)) then
             write(*,77)symbol,ip,kp,trim(line)
77           format('In istpfun: Too long symbol: "',a,'" ',2i4,' on line '/a)
             stop
          endif
!          write(*,*)'Accept "',ch1,'": symbol: "',trim(symbol),'" ',ip,kp
          symbol(kp:kp)=ch1
!          write(*,*)'istpfun symbol: ',trim(symbol),kp,ip
          cycle issymbol
       elseif(kp.ge.2) then
! We found a character illegal in a symbols, check if in nottpfun
          if(kp.ge.2 .and. kp.le.4) then
             do mp=1,5
! Skip symbols: LN LOG EXP ERF GEIN
                if(symbol(1:4).eq.nottpfun(mp)) goto 300
             enddo
! symbol is not predefined function: LN LOG EXP ERF GEIN
!             write(*,*)'reference to tpfun: ',trim(symbol),ip
          endif
          exit issymbol
       endif
! 
300    continue
       symbol=' '; kp=0
    enddo issymbol
!    write(*,*)'istpfun found tpfun: "',trim(symbol),'" ',ip
1000 continue
    return
  end subroutine istpfun

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

!\addtotable subroutine capson & Convert character to UPPER case
!\begin{verbatim}
!  SUBROUTINE capson(text)
! converts lower case ASCII a-z to upper case A-Z, no other changes
!    implicit none
!    character text*(*)
!\end{verbatim}
!    integer, parameter :: lowa=ichar('a'),lowz=ichar('z'),&
!         iup=ICHAR('A')-ICHAR('a')
!    integer i,ich1
!    DO i=1,len(text)
!       ich1=ichar(text(i:i))
!       IF(ich1.ge.lowa .and. ich1.le.lowz) THEN
!          text(i:i)=char(ich1+iup)
!       ENDIF
!    ENDDO
!  END SUBROUTINE capson
!
!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine xtdbentertpfun(tpfun,add)
! check if TPfun is needed and of so enter it in seltpfun
    character*(*) tpfun
    logical add
!
    integer lentp,jj
    lentp=len_trim(tpfun)
    add=.false.
    do ntp=1,nseltp
!       write(*,*)'3EY in xtdbentertpfun: ',&
!            tpfun(1:lentp),' ? ',seltpfun(ntp)%tpfunname(1:lentp),' ',ntp
       if(tpfun.eq.seltpfun(ntp)%tpfunname(1:lentp)) then
          if(seltpfun(ntp)%status.eq.-1) then
             add=.true.
             seltpfun(ntp)%status=1
! sometimes a final N has to be added ....
!             write(*,10)1,trim(wholexpr)
10           format('In xtdbentertpfun: ',i2,' "',a,'"')
             jj=len_trim(wholexpr)
             if(wholexpr(jj:jj).ne.'N') then
! Wow, wholexpr is allocated an it must ne extended like this ...
                wholexpr=wholexpr(1:jj)//' N'
!                write(*,10)2,trim(wholexpr)
             endif
             seltpfun(ntp)%data=trim(wholexpr)
! check if this TPfun need other TPfuns wholexpr is a global variable
             call addmissingtp(tpfun)
             goto 1000
          endif
       endif
    enddo
1000 continue
!    write(*,*)'exit xtdbentertpfun: ',add
    return
  end subroutine xtdbentertpfun

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine xtdbaddbibref(bibref,selbib,nselbib,maxbib)
! adding new bibrefs to selbib if not already there
    integer nselbib,maxbib
    type(ocbib) :: selbib(*)
    character*(*) bibref
    integer nbib,lbib,kk,same
    lbib=len_trim(bibref)
    do nbib=1,nselbib
       if(bibref.eq.selbib(nbib)%bibitem(1:lbib)) goto 1000
    enddo
! it is a new reference
!    if(nselbib.gt.0) write(*,*)'Adding bibref: "',bibref,'" "',&
!         trim(selbib(nselbib)),'"'
    if(nselbib.lt.maxbib) then
       nselbib=nselbib+1
       selbib(nselbib)%bibitem=bibref
       selbib(nselbib)%status=-1
! this is never written?
       write(*,*)' >>> xtdbaddbibref added: "',bibref,'"',nselbib
    else
! check if the same appears several times ...
       write(*,*)'Too many bibliographic references',nselbib
       same=0
       do nbib=1,nselbib
          write(*,*)'Bibref: ',selbib(nbib)%bibitem
          do kk=1,nselbib
             if(kk.ne.nbib) then
                if(selbib(nbib)%bibitem.eq.selbib(kk)%bibitem) then
                   write(*,*)'Duplicate reference',kk,nbib,selbib(kk)%bibitem
                   same=same+1
                endif
             endif
          enddo
       enddo
       write(*,*)'Duplicate references: ',same
    endif
1000 continue
    return
  end subroutine xtdbaddbibref

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

! Application (OC) dependent subroutines called by xtdbread

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine OCenterel(spel,name,refstate,mass,h298,s298)
!  subroutine OCenterel(spel,data)
! enter Element/species from the xtdb file.  Name is full name in OC, ignored
    character spel*2,name*(*),refstate*(*)
    double precision mass,h298,s298
! data is a text: Mass="5.199600E+01" H298="4.050000E+03" S298="2.354290E+01"
! /- and Va introduced automatically
    if(spel.eq.'/-' .or. spel.eq.'VA') goto 1000
    if(nomorelements) then
       write(*,*)'No more elements allowed'; goto 1000
    endif
    nselel=nselel+1
    selel(nselel)%elname=spel
    selel(nselel)%data=' '
! Hm alphabetical order!!  later ....
!    write(*,10)nselel,spel
10  format('Selected element; ',i2,2x,a)
! elements are also entered as species, except /-
    if(spel.ne.'/-') then
       nselsp=nselsp+1
       selsp(nselsp)%species=spel
       selsp(nselsp)%data=spel
       allocate(selsp(nselsp)%elnames(1))
       allocate(selsp(nselsp)%stoicc(1))
       selsp(nselsp)%elnames(1)=spel
       selsp(nselsp)%stoicc(1)=1.0D0
       selsp(nselsp)%charge=0.0D0
!       write(*,*)'Entering species ',spel,nselsp
       selsp(nselsp)%extra=' '
    endif
! enter element also in OC, needed to check species and mqmqa quads
! element symbol, name, reference state, mass, h298, s298
    call store_element(spel,name,refstate,mass,h298,s298)
1000 continue
    return
  end subroutine OCenterel

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine OCenterspecies(spname,stoi,mqmqa,uniquac)
! enter Element/species in xtdb file
    implicit none
    character*(*) spname,stoi,mqmqa,uniquac
! this requires extracting the stoichiometry ... 
    integer ia,ib,nel,kk,mel,kp
    character*2 el(10)
    double precision coef(10),qq
! this nend should be reinitiated when the NEW command
    integer, save :: nend=-1
!
    qq=0.0D0
    if(nomorelements) then
       write(*,*)'No more elements or species allowed'; goto 1000
    endif
    if(mqmqa(1:1).eq.' ') then
! This is NOT MQMQA quad
!    write(*,10)trim(spname),trim(stoi),nselsp
10  format('In OCenterspecies: "',a,'" stoi: "',a,'"',i4)
! check species is not a duplicate
       do ia=1,nselsp
          if(spname.eq.selsp(ia)%species) then
!          write(*,*)'Species "',trim(spname),'" already entered'
             goto 1000
          endif
       enddo
! extract the elements, the electron is included as /+/- here but removed below
       call extractstoi(stoi,nel,el,coef)
       if(xtdberr.ne.0) goto 1000
! xtdberr=5000 means element not entered, this species should be ignored
       mel=nel
! The elements must already be entered!
       thissp: do ia=1,nel
          entered: do ib=1,nselel
             if(el(ia).eq.selel(ib)%elname) cycle thissp
          enddo entered
          if(el(ia).eq.'/+') then
             qq=coef(ia)
             mel=nel-1
          elseif(el(ia).eq.'/-') then
             qq=-coef(ia)
             mel=nel-1
          else
! this species has an unknown element, ignore it
             goto 1000
          endif
       enddo thissp
! enter the species in OC
       call enter_species(spname, nel, el, coef)
       if(gx%bmperr.ne.0) goto 1000
! enter the species in temporary here
       nselsp=nselsp+1
       selsp(nselsp)%species=spname
       write(*,*)'3EY New species "',trim(selsp(nselsp)%species),'" ',nselsp
! Hm here I allocate a place for the electon, maybe use mel?
       allocate(selsp(nselsp)%elnames(mel))
       allocate(selsp(nselsp)%stoicc(mel))
       do ia=1,mel
          selsp(nselsp)%elnames(ia)=el(ia)
          selsp(nselsp)%stoicc(ia)=coef(ia)
       enddo
! charge
       if(qq.ne.0.0d0) then
          selsp(nselsp)%charge=qq
! Here I add charge as last constituent
!       selsp(nselsp)%elnames(nel)='/-'
!       selsp(nselsp)%stoicc(nel)=qq
       else
          selsp(nselsp)%charge=0.0D0
       endif
    else
!----------------------------------------------------
! this IS AN MQMQA QUAD, code copied from gtp3E (for TDB) lines 3868-3888 
       call capson(spname)
       kp=index(spname,'/')
       if(kp.gt.0 .and. &
            spname(kp+1:kp+1).ge.'A' .and. spname(kp+1:kp+1).le.'Z') then
! this is an MQMQA quad, an ion has /+ or /- or /digit           
!          write(*,572)trim(spname),trim(mqmqa)
572       format('3EY Call mqmqa_species: "',a,'" "',a,'" ')
! mqmqa_species in gtp3B.  It will check everything and enter the species in OC
! provided the elements are entered!!!
          call mqmqa_species(spname,mqmqa,nend)
          if(gx%bmperr.ne.0) then
             write(*,*)'3E error creating MQMQA quad',gx%bmperr
             goto 1000
          endif
       endif
!       write(*,*)'3EY spname modified? ',spname,nselsp
! modify the name
       nselsp=nselsp+1
       selsp(nselsp)%species=spname
!       if(eolch(longline,ip)) then
!         if(.not.silent) write(kou,*)'WARNING No stoichiometry for species: ',&
!               trim(name1)
!          tdbwarning=.TRUE.
!         write(*,*)'3E tdbwarning set true 3'
       selsp(nend)%extra=mqmqa
    endif
!    if(uniquac(1:1).ne.' ') selsp(nselsp)%extra=mqmqa
1000 continue
    return
  end subroutine OCenterspecies

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine extractstoi(spstoi,nel,el,coef)
! decode a species stoiciometry as element/number/ without spaces
! The element can be one or two letters and the number an integer, real
! or a quotent such as a/b:  AL2O3 or ALO1.5 or ALO3/2 are the same species
! No parenthesis are allowed Al2(SO4)3 must be written Al2S2O12
! The element names are case insensitve
! A single letter element must have a followng number (or a final space)
! A two letter element without a following number is assumbed to have unity
! A charge is /+ or /- followed by a digit
!
    character*(*) spstoi
    character*2 el(*)
    integer nel
    double precision coef(*)
!
    double precision stf,nom,kvot
    double precision, parameter :: one=1.0D0, zero=0.0D0
    integer lens,kk
    logical ddot,slash,ion
!
    integer ip
    character ch1*1
!    
    debug=.false.
!    debug=.true.
    call capson(spstoi)
    lens=len_trim(spstoi)
!    if(lens.le.1) then
!       write(*,*)'The stoichiometry must be at least 1 character!'
!       xtdberr=5000; goto 1000
! This can be a single element species with default stochiometry 1       
!       coef(1)=one
!    endif
    kvot=zero
    nom=one
    ion=.false.
! nel is incremeneter for each element, ip for each position in spstoi
    nel=0
    ip=1
    ch1=spstoi(ip:ip)
! lge and lle use ASCII character set
    if(.not.(LGE(CH1,'A') .AND. LLE(CH1,'Z'))) then
       write(*,10)trim(spstoi),ip,ch1
10     format(' *** Error in "',a,'" expected element at ',i3,' found "',a,'"')
       xtdberr=5001; goto 1000
    endif
!======================================= big extract loop
    extract: do while(.true.)
! The first letter of an element is in ch1
! It can be first letter of a second element or ion /
       nel=nel+1
! second character of element a space
       el(nel)=ch1
! Set default stoichiometric factor to 1.0
       stf=one
! default stochiometric factor
       coef(nel)=stf
       ip=ip+1
       if(ip.gt.lens) exit extract
! this ch1 is the first character of el(nel), ch1 is updated below
! for ip=1 ch1 is always a letter but
       if(ch1.eq.'/') then
          if(ion) then
             write(*,14)trim(spstoi),ip
14           format('*** Error in "',a,'" species already charged.',i3)
             xtdberr=5017; goto 1000
          endif
! This cannot be the first slement or a / in a factor as 1/3 
! But it can be the first letter of second element or ion "/"
! An ion, "/"  must be followed by a + or - to represent an ion
! ip already incremented above
          if(ip.gt.lens) then
             write(*,11)trim(spstoi),ip
11           format('*** Error ending "',a,'"  with /',i3)
             xtdberr=5011; goto 1000
          elseif(spstoi(ip:ip).eq.'+' .or. spstoi(ip:ip).eq.'-') then
             ion=.true.; ddot=.true.
! This is an ion, set a default charge 1
!             write(*,300)trim(spstoi),ip,nel,(coef(kk),kk=1,nel)
300          format('Ion 1: "',a,'" '2i3,10F8.3)
             el(nel)=spstoi(ip-1:ip)
             coef(nel)=one
! jump to extract a digit, the valence must be the end of an ion
             ip=ip+1
             if(ip.gt.lens) exit extract
             ch1=spstoi(ip:ip)
! there cannot be any new elements after a /+ or /-
             goto 500
          else
             write(*,13)trim(spstoi),ip
13           format('*** Error charge of "',a,'" must be "/+" or "/-" ',i3)
             xtdberr=5009; goto 1000
          endif
       endif
! ip was incremented above, second letter of an element, a factor or charge
       ch1=spstoi(ip:ip)
       if(debug) write(*,*)'second letter 1: ',ch1,ip
! If second letter indicate an ion jump back to extract the whole symbol
       if(ch1.eq.'/') cycle extract
       if(lge(ch1,'A') .and. lle(ch1,'Z')) then
! this must be the second letter of the element
! Save full name, the charge taken care of above
          el(nel)(2:2)=ch1
          ip=ip+1
          if(ip.gt.lens) exit extract
          ch1=spstoi(ip:ip)
! if the next character is a letter cycle, otherwise coefficient
          if((lge(ch1,'A') .and. lle(ch1,'Z')) .or. ch1.eq.'/') cycle extract
       endif
!-------------------------------------
! We arrive here to extract the stoichiometry or valency
! It can start with a digit or a decimal point
! reset the default stochiometry
       if(debug) write(*,*)'third letter 2: ',ch1,ip
       stf=zero; nom=1.0D0
       ddot=.false.; slash=.false.
! jump here if ion
500    continue
       stf=zero
       coefficient: do while(.true.)
          if(debug) write(*,*)'third letter 3: ',ch1,ip
          if(ch1.eq.'.') then
! Handle a decimal point inside a real number
             if(ddot) then
! not allowed for ions either
                if(debug) write(*,20)trim(spstoi),ip
20              format('*** Error in "',a,'" two decimal dots ',i3)
                xtdberr=5002; goto 1000
! the following numbers will be decimals
             endif
             ddot=.true.
             nom=1.0D-1
             ip=ip+1
             if(ip.gt.lens .and. stf.eq.zero) then
! stoichiometry ends with a ., if stf=0.0 no previous digit
                if(debug) write(*,30)trim(spstoi),ip
30              format('*** Error in "',a,'" missing stoichiometry ',i3)
                xtdberr=5003; goto 1000
             endif
             ch1=spstoi(ip:ip)
          endif
! now there must be a number!!!
          if(debug) write(*,*)'third character 4: ',ch1,ip
          stoik: do while(lge(ch1,'0') .and. lle(ch1,'9'))
! extract stoiciometic factor digit by digit
! stf is the previous numbers in the stochiometric factor, initially 0.0
             if(debug) write(*,*)'third character 5: ',ch1,ip
             if(ddot) then
                stf=stf+nom*(ichar(ch1)-ichar('0'))
                nom=1.0D-1*nom
             else
                stf=stf*nom+ichar(ch1)-ichar('0')
                nom=1.0D1*nom
             endif
             if(debug) write(*,35)'third character 6: ',ch1,ip,nel,stf,nom,kvot
35           format(1x,a,a,2i3,5F10.4)
             ip=ip+1
! we have an element with a stochiometric factor, OK if no more
! the stoichometry of last element set after extract
             if(ip.gt.lens) exit extract
             ch1=spstoi(ip:ip)
             if(debug) write(*,36)1,ch1,ip,nel,stf,nom,kvot
36           format('fourth character ',i2,': "',a,'"',2i3,5F10.4)
! there can be a / in stoichiometry, for example in AlO3/2,  This is not an ion
             if(ch1.eq.'/') then
                if(spstoi(ip+1:ip+1).eq.'+' .or. &
                     spstoi(ip+1:ip+1).eq.'-') then
! this is not a division, it is the electone!
                   goto 600
                endif
                if(slash) then
                   write(*,40)trim(spstoi),ip
40                 format('**** Error in "',a,'" two slashes',i3)
                   xtdberr=5004; goto 1000
                endif
                slash=.true.
                if(ddot) then
                   write(*,50)trim(spstoi),ip
50                 format('**** Error in "',a,'" both slash and dot!',i3)
                   xtdberr=5005; goto 1000
                endif
                if(stf.eq.zero) then
                   write(*,60)trim(spstoi),ip
60                 format('**** Error in "',a,'" no digits before slash!',i3)
                   xtdberr=5006; goto 1000
                endif
! kvot is set to current value of stf, stf will be value to divide with
                kvot=stf
                stf=zero
                ip=ip+1
                if(ip.gt.lens) then
                   write(*,70)trim(spstoi),ip
70                 format('**** Error in "',a,'" no digits after slash',i3)
                endif
                if(debug) write(*,16)2,ch1,ip,nel,stf,nom,kvot
16              format('Position ',i2,' letter "',a,'" ',2i3,10F8.4)
                ch1=spstoi(ip:ip)
! after a / there must be a digit, a "." or letter not allowed
                if(debug) write(*,16)3,ch1,ip,nel,stf,nom,kvot
                cycle stoik
             endif  ! we have taken care of a /
             if(debug) write(*,16)4,ch1,ip,nel,stf,nom,kvot
! if ch1 is not a digit exit here.
             if(stf.eq.zero) then
                write(*,80)trim(spstoi),ip,ch1
80              format('*** Error in "',a,'" digit error at',i3,' "',a,'"')
                xtdberr=5009; goto 1000
             endif
          enddo stoik
          if(debug) write(*,16)5,ch1,ip,nel,stf,nom,kvot
! there can be a . inside a stoichiometric factor followed by digits
          if(ch1.eq.'.') then
             if(slash) then
                write(*,50)trim(spstoi),ip
                xtdberr=5007; goto 1000
             endif
             cycle coefficient
          endif
600       continue
          if(debug) write(*,16)6,ch1,ip,nel,stf,nom,kvot
! we have to calculate the stoichiometric factor
! If kvot=0.0 then coef is stf
          if(kvot.eq.zero) then
             coef(nel)=stf
          else
             coef(nel)=kvot/stf
          endif
          if(debug) write(*,16)7,ch1,ip,nel,stf,nom,kvot,coef(nel)
          stf=zero
          kvot=zero
          nom=one
          if(debug) write(*,16)8,ch1,ip,nel,stf,nom,kvot,coef(nel)
          cycle extract
       enddo coefficient
! We come here if ch1 is not a digit, it can be a letter or / or ???
       if(debug) write(*,16)9,ch1,ip,nel,stf,nom,kvot,coef(nel)
       if(ch1.eq."/" .or. (ch1.ge.'A' .and. ch1.le.'Z')) cycle extract
!
       write(*,90)trim(spstoi),ip
90     format('*** Error in "',a,'" illegal character at ',i3)
       xtdberr=5023; goto 1000
    enddo extract
    if(ip.gt.lens) then
! we have to calculate the stoichiometric factor
! If kvot=0.0 then coef is stf
       if(kvot.eq.zero) then
          coef(nel)=stf
       else
          coef(nel)=kvot/stf
       endif
       if(debug) write(*,16)8,ch1,ip,nel,stf,nom,kvot,coef(nel)
    endif
!
!    write(*,*)'Leaving extractstoi: ',trim(spstoi),nel
!    do ip=1,nel
!       write(*,100)el(ip),coef(ip)
100    format('Element: "',a,'" Stochiomentry ',1pe16.6)
!    enddo
1000 continue
    return
  end subroutine extractstoi
  
!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine OCenterphase(phrec)
! listing phase data to be entered in OC
!    type(phnest), pointer :: phrec
    type(phnest) :: phrec
    integer ms,ns,ip,jp,js,lk
    character*24 constituent,ch1*1
    character*1024 ccxx
    logical none,giveup
!
    if(.not.nomorelements) then
! when first phase enetered we cannot entoer more elements/species
! we should arrange elements and species in alphanetical order!
       nomorelements=.TRUE.
       call alphabetical_order
    endif
!
    ms=size(phrec%clist)
!    write(*,7)trim(phrec%Id),ms,(phrec%clist(ns)%list,ns=1,ms)
7   format(/' in OCenterphase: ',a,'" ',i3/20(' "',a,'"'))
! we must check if the phase can exist, i.e. is there at least one
! entered species in each sublattice
    nselph=nselph+1
! maybe not needed?
!    allocate(cc(10))
!    write(*,*)'Try assigning ccxx'
    ccxx=' : '
!    write(*,*)'Success assigning ccxx'
    selph(nselph)%nsublat=ms
    selph(nselph)%mult=phrec%mult
! not any ternary extrapolation methods yet
    nullify(selph(nselph)%terxpol)
    sublatt: do ns=1,ms
! loop for all constituent in the phase in the database
       giveup=.false.
       none=.true.
!       selph(nselph)%const(ns)=' '
       ip=1
! extract database constituents in sublattice ns from first position
! with allocated characters there are not always a space at the end!!
       lk=len_trim(phrec%clist(ns)%list)
       any: do while(.not.xeolch(phrec%clist(ns)%list,ip))
          jp=index(phrec%clist(ns)%list(ip:),' ')
!          write(*,*)'Loop sublattice ',ns,ip,lk,jp
! ip is start of constituent, find terminating space, if endofline
!
! big problem here because %list does not terminate with spaces
! one has to be careful with jp=0 as there are no trailing spaces ....
          if(jp.eq.0 .and. ip.le.lk) then
             constituent=phrec%clist(ns)%list(ip:lk)
!             write(*,*)'Constituent ip:lk :',constituent,ip,lk
! giveup set to true here as jp=0 means it is tha last species on the line
! it would probably be safer to add a trailing space to %list
             giveup=.TRUE.
             jp=lk-ip+1
! infinite number of problems here with C and V and character strings
! that does not terminate with spaces !!! SUCK
          else
             constituent=phrec%clist(ns)%list(ip:ip+jp-1)
          endif
!          write(*,71)'Constituent: "',trim(constituent),'"',ip,jp,ip+jp,lk
71        format(a,a,a,5i3)
          ip=min(ip+jp,lk)
          do js=1,nselsp
! if the constituent is entered phase can exitst
             if(jp.gt.0) then
!                write(*,*)ns,'compare: "',constituent(1:jp),&
!                     '" "',selsp(js)%species(1:jp),'"',ip,js
! wow V was accepted twice as it matched VA also ....
                if(constituent(1:jp).eq.trim(selsp(js)%species)) then
! this constituent is entered, append it to selph(nselph)%const
                   ccxx=trim(ccxx)//' '//trim(constituent)
!                   write(*,*)'Saving constituent in ccxx',trim(ccxx)
                   none=.false.
                endif
             else
                exit any
             endif
          enddo
!          write(*,*)'Loop index ',ip,giveup
          if(giveup) exit any
       enddo any
       if(none) goto 1100
       ccxx=trim(ccxx)//' : '
!       write(*,*)'Next sublattice: "',trim(ccxx),'" ',ns
    enddo sublatt
! there is at least one entered constituent in each sublattice in selph%const
    selph(nselph)%phasename=phrec%Id
!    write(*,10)selph(nselph)%phasename,ms,trim(ccxx)
10  format('OK enter phase: ',a,i3,a)
!    write(*,*)'The remaining problem is to transfer cc(1:mn) to selph%const!'
    selph(nselph)%const=trim(ccxx)
    selph(nselph)%confent=trim(phrec%confent)
    if(allocated(phrec%amendph)) then
!       write(*,*)'Phase ',trim(phrec%id),' has models ',trim(phrec%amendph)
       selph(nselph)%amendph=phrec%amendph
!    else
!       write(*,*)'No amend phase allocated for: ',trim(phrec%Id)
    endif
    if(allocated(phrec%dispar)) then
!       write(*,*)'Phase ',trim(phrec%id),' has dispart ',trim(phrec%dispar)
       selph(nselph)%dispar=phrec%dispar
    endif
! this is not needed, at the end the records in firstxpol will be searched
!    xpol=>firstxpol
!    nullify(lastxpol)
! look for any ternary extrapolations for this selph(nselph)
!    if(associated(xpol)) then
!       write(*,*)'3EY xpol: "',xpol%phase,'" and "',selph(nselph)%phasename.'"'
!       if(xpol%phase.eq.selph(nselph)%phasename) then
!          write(*,*)'3EY adding ternary method',xpol
!          if(xpol.eq.firstxpol) then
!             firstxpol=>xpol%next
!          else
!             xpol2=>xpol
!          xpol2%next=>selph(nselph)%terxpol
!          selph(nselph)%terxpol=>xpol2
!       endif
!       xpol=xpol%next
!    endif
1000 continue
!    write(*,*)'Press return to continue'
!    read(*,'(a)')ch1
!
    return
1100 continue
! We arrive here there is a sublattice with no constituents entered in a subl
! The phase cannot be entered
!    write(*,5)phrec%Id
!    deallocate(selph(nselph)%const)
5   format('The phase ',a,' cannot exist in this system?')
    nselph=nselph-1
    goto 1000
  end subroutine OCenterphase

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine alphabetical_order
! this should arrange elements and species in alphabetical order as in TDB
    integer nn,ia,ib,ic,mel,kk,jj
    character bytel*2,ch1*1,byt2*24
    character*24, dimension(:), allocatable :: ss
    logical, save :: notdone=.true.
    integer, dimension(:), allocatable :: orig
!    write(*,*)'In alphabetical order',notdone
    if(notdone) then
!       write(*,*)'In alphabetical order to arrange elements and species'&
!            ' in alphabetical order'
       notdone=.false.
       ic=1

       ord1: do while(ic.gt.0)
          ic=0
! element /- is -1 and Va is 0
          do ia=1,nselel-1
             if(selel(ia)%elname.gt.selel(ia+1)%elname) then
                ic=1
! shift data in ia and ia+1, very clumsy
                selel(maxtdbel)%elname=selel(ia)%elname
                selel(maxtdbel)%data=selel(ia)%data
                selel(ia)%elname=selel(ia+1)%elname
                selel(ia)%data=selel(ia+1)%data
                selel(ia+1)%elname=selel(maxtdbel)%elname
                selel(ia+1)%data=selel(maxtdbel)%data
             endif
          enddo
          if(ic.eq.0) exit ord1
       enddo ord1
!       write(*,10)(selel(nn)%elname,nn=1,nselel)
10     format('Elements entered: '20(a,1x))
! this is very clumsy
!       write(*,*)'Order species alphabetically in selspord'
!       allocate(selspord(nselsp+1))
!       allocate(ss(nselsp+1))
!       allocate(orig(nselsp+1))
       allocate(selspord(maxtdbsp))
       allocate(ss(maxtdbsp))
       allocate(orig(maxtdbsp))
       do ia=1,nselsp
          selspord(ia)=ia
          orig(ia)=ia
          ss(ia)=selsp(ia)%species
       enddo
!
       kk=0
       ord2: do while(.true.)
! there is no /- species and Va is in alphabetical order
          ic=0
!          kk=kk+1
          do ia=1,nselsp-1
!             write(*,*)'Comparing ',ia,' "',trim(ss(ia)),&
!                     '" and "',trim(ss(ia+1)),'"'
             if(ss(ia).gt.ss(ia+1)) then
                ic=1
                byt2=ss(ia)
                ss(ia)=ss(ia+1)
                ss(ia+1)=byt2
                jj=orig(ia)
                orig(ia)=orig(ia+1)
                orig(ia+1)=jj
                selspord(ia+1)=jj
             endif
          enddo
          if(ic.eq.0) exit ord2
       enddo ord2
       do ia=1,nselsp
!          write(*,90)ia,trim(ss(ia)),orig(ia)
! selspord is in alphabetical order and has index to species
          selspord(ia)=orig(ia)
       enddo
    endif
    return
  end subroutine alphabetical_order

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine OCenterpar(parname,skip,bibref)
! enter Parameter in OC if phase and constituents entered, FUNCTION IN wholexpr
    character*(*) parname,bibref
! bibfref can be just a space if not used
    integer skip
! return skip -1 or -2 if not needed
!
    character*24 phase, constituent
    character MPID*8
! the call to addmissingtp should only be used if the parameter is entered
    character*24 phasename
    integer nn,lk,lm,ip,lph,kk

! Check that the parameter needed, i.e. phase and constituents entered
!    write(*,7)trim(parname)
7   format(/'In OCenterpar: ',a)
    call checkifparisneeded(parname,lph)
    
    if(lph.lt.0) then
! lph=-1 means phase not entered
!       if(lph.eq.-1) write(*,*)'Parameter not entered as phase not present'
!       if(lph.eq.-2) write(*,*)'Parameter not entered as  constituent missing'
       skip=lph
       goto 1000
    endif
!-----------------------
! Enter this parameter, the parameter may need new TPfuns, check wholexpr
    call addmissingtp(parname)
!
! enter expr with bibref and add bibref to be added
    wholexpr=wholexpr//' '//bibref
    nselpar=nselpar+1
    if(nselpar.gt.maxpar) then
       write(*,*)'Too many parameters ',maxpar
    else
       selpar(nselpar)%parname=trim(parname)
       selpar(nselpar)%data=trim(wholexpr)
    endif
! if entered in OC then add bibref to selbibref
    lk=len_trim(bibref)
    if(lk.gt.0) then
       do nn=1,nselbib
! exit if already in selbib
          if(bibref.eq.selbib(nn)%bibitem(1:lk)) goto 300
       enddo
       if(nselbib.ge.maxbib) then
          write(*,*)'Too many bibliographic references 2 ',maxbib
! check if the same appears several times ...
          do nn=1,nselbib
             write(*,*)'Bibref: ',selbib(nn)%bibitem
             do kk=1,nselbib
                if(kk.ne.nn) then
                   if(selbib(nn)%bibitem.eq.selbib(kk)%bibitem) then
                      write(*,*)'Duplicate reference',kk,nn,selbib(kk)%bibitem
                   endif
                endif
             enddo
          enddo
       else
          nselbib=nselbib+1
          selbib(nselbib)%bibitem=bibref
! mark this as missing
          selbib(nselbib)%status=-1
       endif
    endif
300 continue
! exit here if parameter not needed
1000 continue
    return
  end subroutine OCenterpar

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine checkifparisneeded(parname,lph)
! check if phase and constituents present in this parameter
! The phase name may be abbreviated between "_" characters
    implicit none
    character*(*) parname
    character*24 phasename,constituent
    character(len=:), allocatable :: mpid
    integer lk,lm,ln,lph,lcomma,lcolon,lsemic,lright,lend
    integer js,jmod,jmpid
!
    lph=0
    lk=index(parname,'(')
    if(lk.le.0) then
       write(*,*)'3EY Missing ( in parameter: ',trim(parname)
       xtdberr=4544; goto 1000
    endif
    mpid=parname(1:lk-1)
!    write(*,*)'MPID: "',mpid,'"',lk
! check MPID is OK (must not be abbreviated)
!    do ln=1,noofmpid
!       if(mpid.eq.mpidok(ln)) goto 110
!    enddo
!---------------------------------------------------------
!  VERY SPECIAL, always accept a single L as G.  Even for endmembers ...
    if(mpid.eq.'L') mpid='G'
    if(mpid.eq.'G') goto 110
! otherwise we have to find which model it belongs to ...
!---------------------------------------------------------

! There are 5 models with MPID, 3 magnetic, Einstein, liq2state
    do jmod=1,5
!       write(*,*)'Testing ',mpid,'. number of MPIDs for jmod ',&
!            jmod,xtdbmodel(jmod)%nmpid
       do jmpid=1,xtdbmodel(jmod)%nmpid
          if(mpid.eq.xtdbmodel(jmod)%mpid(jmpid)(1:lk-1)) then
!             write(*,*)'3EY Found model of MPID: "',mpid,'" is ',jmod,jmpid
             goto 110
          endif
!          write(*,108)mpid,xtdbmodel(jmod)%mpid(jmpid)(1:lk-1),jmod,jmpid
108       format('3EY MPIDs: "',a,'" and "',a,'" indices: ',2i4) 
       enddo
    enddo
! now we use gtp_xtdbcompatibility: xtdbmodel
! we have an unknown MPID
    write(*,*)'3EY Unknown Model Parameter IDentification (MPID): "',mpid,'"'
    write(*,*)'Parameter: "',trim(parname),'"'
    xtdberr=4546
    goto 1000
!---------------------------------------------------------
110 continue
! there must be a , after the phase name lk+1 is the first letter
    lm=index(parname(lk+1:),',')
    if(lm.le.0) then
       write(*,*)'3EY Missing "," after phase name in parameter: ',trim(parname)
       write(*,*)'Parameter: "',trim(parname),'"'
       xtdberr=4545; goto 1000
    endif
!
    phasename=parname(lk+1:lk+lm-1)
    lk=lk+lm+1
! check this is a selected phase, allowing abbreviations
!    write(*,*)'Calling findabbrphname 1: "',phasename,'"'
    call findabbrphname(phasename,lph)
    if(lph.le.0) then
       write(*,*)'Phase not selected "',phasename,'" parameter ignored'
       goto 1000
    endif
! if lph>0 check also that the constituents are present in correct sublattices
! a comma deparates constituents in a sublattice, a colon in different
! the semicolon and ) is the end of constituents
    lcomma=index(parname(lk:),',')
    lcolon=index(parname(lk:),':')
    lsemic=index(parname(lk:),';')
    lright=index(parname(lk:),')')
!    write(*,195)trim(parname(lk:)),lcolon,lcomma,lsemic,lright,lk
195 format('Constituent array "',a,'" ',5i5)
    if(lcomma.eq.0 .and. lcolon.eq.0) then
! the single constituent is terminated by ; or )
       if(lsemic.eq.0) then
! there is always )
          constituent=parname(lk:lk+lright-2)
       else
! if there is a ; that comes before )
          constituent=parname(lk:lk+lsemic-2)
       endif
!       write(*,*)'Constituent 1: ',trim(constituent)
! A wildcard * in a single sublattice phase (rare!) is accepted
       if(constituent(1:1).eq.'*') goto 500
       do js=1,nselsp
          if(constituent.eq.selsp(js)%species) goto 500
       enddo
       goto 1100
    else
! we have to handle several sublattices and interacting constituents
! I do not understand the -4, I thought -2 would be correct
       if(lsemic.eq.0) then
          lend=lk+lsemic-4
       else
          lend=lk+lright-4
       endif
!       write(*,*)'More: "',trim(parname(lk:)),'"',lcolon,lcomma,lend
! loop to extract all constituents terminated by , or : before lend
       more: do while(.true.)
          if(lcolon.gt.0 .and. lcomma.gt.0) then
             if(lcolon.lt.lcomma) then
                constituent=parname(lk:lk+lcolon-2)
!                write(*,*)'Constituent 2: ',trim(constituent)
                lk=lk+lcolon
                lcomma=lcomma-lcolon
! lk is updated to position after the constutent found
                lcolon=index(parname(lk:),':')
!                if(lcolon.gt.0) lcolon=lk+lcolon
!               write(*,*)'Remaining: "',trim(parname(lk:)),'"',lk,lcolon,lcomma
             else
                constituent=parname(lk:lk+lcomma-2)
!                write(*,*)'Constituent 3: ',trim(constituent)
                lk=lk+lcomma
                lcolon=lcolon-lcomma
                lcomma=index(parname(lk:),',')
!                if(lcomma.gt.0) lcomma=lk+lcomma
!               write(*,*)'Remaining: "',trim(parname(lk:)),'"',lk,lcolon,lcomma
             endif
! accept a wildcard "*" only if all other constituents selected
             if(constituent(1:1).ne.'*') then
                do js=1,nselsp
                   if(constituent.eq.selsp(js)%species) cycle more
                enddo
             endif
! constituent is not not selected
             goto 1100
          endif
!          write(*,*)'Only comma or colon',lcolon,lcomma,lend
          if(lcolon.gt.0) then
             constituent=parname(lk:lk+lcolon-2)
!             write(*,*)'Constituent 4: ',trim(constituent)
             lk=lk+lcolon
             lcolon=index(parname(lk:),':')
!             if(lcolon.gt.0) lcolon=lk+lcolon
!             write(*,*)'Remaining: "',trim(parname(lk:)),'"',lk,lcolon,lcomma
          elseif(lcomma.gt.0) then
             constituent=parname(lk:lk+lcomma-2)
!             write(*,*)'Constituent 5: ',trim(constituent)
             lk=lk+lcomma
             lcomma=index(parname(lk:),',')
!             if(lcomma.gt.0) lcomma=lk+lcomma
!             write(*,*)'Remaining: "',trim(parname(lk:)),'"',lk,lcolon,lcomma
          else
             constituent=parname(lk:lend)
!             write(*,*)'Constituent 6: ',trim(constituent),lk,lend,lsemic
             lk=lend
          endif
! always accept wildcard *
          if(constituent(1:1).eq.'*') goto 500
          do js=1,nselsp
             if(constituent.eq.selsp(js)%species) then
                if(lk.eq.lend) then
                   goto 500
                else
                   cycle more
                endif
             endif
          enddo
! constituent is not entered
          goto 1100
! maybe there are more constituents, unless lk is lk is ge lend
200       continue
          if(lk.ge.lend) exit more
          cycle more
       enddo more
    endif
! we have the phase and species, enter parameter
500 continue
!    write(*,*)'Phase and constituents OK'
1000 continue
    return
! parameter contain species not selected
1100 continue
!    write(*,*)' *** Parameter skipped as constituent not selected ***'
    lph=-2
    goto 1000
    return
  end subroutine checkifparisneeded

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine findabbrphname(phasename,lp)
    character phasename*(*)
    integer lp ! is position after phase name ?? or phase index??
!----------------
    integer lpp,lpx,ok1,lenph
    logical debug
    character*1 chp,chx
!    write(*,9)phasename
9   format('In findabbrphname: "',a,'"')
    debug=.false.
    lenph=len(phasename)
    ok1=0  ! index of any previous matching phase
!
! phase names may be abbreviated between each "_". ! It must start with A-Z
! An array with selcted phases in selph(array)%phasename
! It can be a bit complicated, the phase names are in arbitrary order
    lpp=0 ! position in phasename
    lp=1  ! array index in selph(array)
    lpx=0 ! position in selph(array)%phasename(lpx:lpx)
    bigloop: do while(.true.)
       lpp=lpp+1
       if(lpp.gt.lenph) then
! this is lenght of provide phase and we have match up to this position, accept
! normally a phase name ends with a space but with allocated characters ...
!          write(*,*)'Is "',phasename,'" same as "',selph(lp)%phasename,'"?'
          lph=lp; goto 1000
       endif
       chp=phasename(lpp:lpp)
       lpx=lpx+1
       chx=selph(lp)%phasename(lpx:lpx)
! if same character, compare next characters in lpp and lpx
!       write(*,*)'Letters "',chx,'" "',chp,'" position ',lpp,lpx
       if(chp.eq.chx .and. chp.ne.' ') cycle bigloop
! The first character MUST be the same
       if(lpp.gt.1) then
          if(chp.eq.' ') then
! trailing characters in selph(lp)%phasename irrelevant but check ambiguous
! up to now we have found the same characters but it can be ambigous ...
             if(ok1.gt.0) then
                write(*,30)phasename,selph(ok1)%phasename,selph(lp)%phasename
30              format('Ambiguous phase name: "',a,'"'/'"',a,'" and "',a,'"')
             else
! save this then check the remaining phases
                ok1=lp;
                lp=lp+1; lpp=0; lpx=0
             endif
          elseif(chp.eq.'_') then
! chx is not "_", skip in selph(lp)%ohasename up to "_".  If no _ skip
! if we find a "_" continue compare the characters following this
             do while(chx.ne.'_')
                if(lpx.ge.len(selph(lp)%phasename)) then
! we have reached the end of selph(lp), skip this phase
                   lp=lp+1; lpp=0; lpx=0
                   cycle bigloop
                endif
                lpx=lpx+1; chx=selph(lp)%phasename(lpx:lpx)
                if(chx.eq.' ') then
! we do not find any "_" but a space, skip this phase
                   lp=lp+1; lpp=0; lpx=0
                   cycle bigloop
                endif
             enddo
! we found a "_" in selph(lp)%phasename. Backtrack both lpp and lpx
             lpx=lpx-1; lpp=lpp-1
             cycle bigloop
          endif
       endif
! not the same character, compare with next phase in selph(array)
       lp=lp+1
       if(lp.le.nselph) then
          lpp=0; lpx=0
          cycle bigloop
       else
! we have compared with all phases in selph(1..nselph)
          if(ok1.gt.0) then
!             write(*,*)
50           format('We found phase: "',a,'" and "',a,'"')
             lp=ok1
          else
             lp=-1
          endif
          exit bigloop
       endif
    enddo bigloop
1000 continue
    return
  end subroutine findabbrphname

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine OCentertpfun(tpfuname)
! enter TPfun data in OC.  MISSING check if needed. Function in wholexpr
    character*(*) tpfuname
! the call to addmissingtp should only be used if the parameter is entered
    integer nn
! check if tpfun uses other tpfun
    call addmissingtp(tpfuname)  ! the wholexpr is a global variable
! add bibref here otherwise it may be entered as missing TPfun       
!    write(*,10)tpfuname,trim(wholexpr)
10  format('In xtdbOCfun: ',a,2x,a,2x,a)
    return
  end subroutine OCentertpfun

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine OCenterbibitem(bibitem,text)
! enter referenced bibitem in OC
    character*(*) bibitem,text
! should only be used if the parameter is entered
    integer nbib,ll
    write(*,*)'In OCenterbibitem'
    ll=len_trim(bibitem)
    do nbib=1,nselbib
       if(bibitem.eq.selbib(nbib)%bibitem(1:ll)) goto 200
    enddo
! bibitem not found
    goto 1000
! found bibitem
200 continue
    nselbib=nselbib+1
    selbib(nselbib)%bibitem=bibitem
    selbib(nselbib)%data=text
    selbib(nselbib)%status=1
1000 continue
    return
  end subroutine OCenterbibitem

! end module xtdblib

