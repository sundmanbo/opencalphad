! program to extract dokumentation of a Fortran source file
! written by Bo Sundman 2016-2019
!
Program makedok
  integer, parameter :: maxtab=500
  character pfil*64,dfil*64,sfil*64,line*80,dline*80,curfil*64
  character beginxverb*18,endxverb*16,ch1*1
  character verbbuff(500)*80,lastverb*64,sameverb*64,texverb*64
  character verbanew(500)*80
  character nounderscore1*80,nounderscore2*40
  character tablentries(maxtab)*80,tablefile(maxtab)*40
  character sectext*80
  integer, dimension(maxtab) :: tabord
  integer, dimension(3,500) :: nsecverb
  logical dokfil,once,lend,dend,EOF,merge,newverb,percentplus
  integer foundverb,includedverb
  beginxverb='!\begin{verbatim} '
  endxverb='!\end{verbatim} '
  write(*,*)'ASCII value of backslash: ',ichar('\')
  write(*,10)
  texverb=' '
  ntablentries=0
  includedverb=0
  percentplus=.false.
!
! Below the idea is to extract the text between 
! 1. the closest \section BEFORE current \verbatim and
! 2. the closest BEFORE THE NEXT \verbatim.
! maybe better change to the text
! 1. the closest \section AFTER THE PREVIOUS \verbatim to
! 2. to closest \section AFTER THE CURRENT \verbatim
! The important thing is to keep text belonging together together
! and not to miss any documentation text
!
! Added that a source file can have (one level) of include files
!
! The extraction stops if the first line after \verbatim in the documentation
! does not fit the first line found in the source code and ask for manual
! editing of this.  It turned out to be very complicated to handle when new
! subroutines had been added or shifted place.
!
10 format(///'   This is a program to generate and update documentation',&
        ' of software.'//&
        '     written by Bo Sundman 2016-2019 for OpenCalphad, version 2',//&
        'It expects two files, one with the software code',/&
        'and one with existing LaTeX documentation (which can be empty)',//&
        'The idea is that "critical" parts of the code should be included ',&
        'in the '/&
        'documentation and that whenever such a section in the source ',&
        'code has been '/&
        'updated such changes will be detected by this program and replaced.'/&
        'A critical part is typically global data declarations and ',&
        'subroutines and '/&
        'functions with their arguments.'//&
        'Updating the documentation of a developing software is a complex ',&
        'task '/&
        'and this software tries to help with this.  It searches the LaTeX ',&
        'and code files '/&
        'sequentially comparing the critical parts in the two files:'/&
        '- If a critical part has changed in the source code ',&
        'has changed'/'  it is simply replaced in the LaTeX file.'/&
        '- If an critical part existing in the LaTeX file is missing in ',&
        'the source code'/&
        '  the program stops and demands a manual update of the LaTeX file.'/&
        '- If a critical part existing in the source code is missing ',&
        'in the LaTeX file'/&
        '  the program stops and demands a manual update of the LaTeX file.'//&
        'The critical parts of the software that is included in the',&
        ' documentation'/'must be enclosed by lines "!\begin{verbatim}"',&
        ' and "!\end{verbatim}".'/&
        'There must be an exact match of the text on the line after the '/&
        '"\begin{verbatim} ..." line in the LaTeX and code files.'//&
        'There can also be a line !\addtotable <text> in the code ',&
        'for things '/&
        '(as subroutines) that should be added to a LaTeX table'//&
        'Each critical part must be a separate "\(sub)section" ',&
        'in the LaTeX file.'//&
        '"!\end{verbatim} %+" means the next critical part is ',&
        'merged with the current.'//&
        'All text in the LaTeX file between ',&
        'the preceeding \{sub}..section up to the next '/&
        '\{sub}..section preceeding the next \begin{verbatim}, '/&
        'will be copied to the new LaTeX documentation'/'file.'//&
        'If there are differences inside a critical part the old version ',&
        'be included as a LaTeX'/&
        'coment together with the new in the documentation file.'//&
        'The input files will never be changed.'//&
        'Program file name (.F90):')
  lastokverb=0
  read(*,20)pfil
  k=index(pfil,'.')
  if(k.eq.0) then
     k=len_trim(pfil)
     pfil(k+1:)='.F90 '
  endif
  curfil=pfil
20 format(a)
  write(*,30)
30 format('Current documentation file name (.tex, if none give return):')
  read(*,20)dfil
  if(dfil(1:1).eq.' ') then
! no previous documentation file
     sfil=pfil(1:k)//'.tex'
     dokfil=.false.
     write(*,*)'no previous dokfile, writing: ',sfil(1:len_trim(sfil))
  else
     dokfil=.true.
     k=index(dfil,'.')
     if(k.eq.0) then
        sfil=dfil(1:len_trim(dfil))//'_new.tex '
        dfil(len_trim(dfil)+1:)='.tex '
     else
        sfil=dfil(1:k-1)//'_new.tex '
     endif
  endif
! source code file
  lunlevel=1
  lunf90=31
  open(lunf90,file=pfil,access='sequential',status='old')
! new documentation file
  write(*,*)'Writing on new documentation file: ',trim(sfil)
  open(23,file=sfil,access='sequential',status='unknown')
  nverbsec=0
!
37 continue
  if(dokfil) then
! This is the LaTeX file
     open(22,file=dfil,access='sequential',status='old')
! search for first \begin{verbatim}
     linoldok=0
     dend=.false.
     idum=-1
     call searchverb(22,linoldok,idum,sectext,dend)
     rewind(22)
     if(dend) then
! no \begin{verbatim} in current dokumentation ... just ignore old dokfile
        write(*,*)'No verbatim, old documentation file ignored 1'
        dokfil=.false.
        close(22)
        goto 37
     endif
!     write(*,*)'Found first \begin{verbatim} at line: ',linoldok
! search for last \{sub}section before first \begin{verbatim}
     linsec=0
40   continue
     lastlinsec=linsec-1
     call searchsection(22,linsec,.FALSE.,dend)
     if(dend) then
! no \section{ in current dokumentation ... just ignore old dokfile
        write(*,*)'Old documentation file ignored 2'
        dokfil=.false.
        close(22)
        goto 37
     endif
     if(linsec.lt.linoldok) goto 40
! we have bypassed the first \begin{verbatim}, last before is line lastlinsec
     rewind(22)
     nw=1
     i=1
!     write(*,*)'calling copytonewdoc',lastlinsec
     nld=lastlinsec+1
     call copytonewdoc(22,23,i,lastlinsec,nw,.false.)
! now scan the whole document 
!     write(*,*)'Calling scandoc at ',nld
     call scandoc(22,nld,nverb,nsecverb)
     write(*,44)nld,nverb
44   format('Old documentation has ',i5,' lines with ',i4,' verbatim sections'/)
!     do i=1,nverb
!        write(*,*)'Lines: ',(nsecverb(j,i),j=1,3)
!     enddo
  else
! no old docfile, write default preamble on new
     write(23,90)
90 format('\documentclass[12pt]{article}'/'\usepackage[latin1]{inputenc}'/&
        '\topmargin -1mm'/'\oddsidemargin -1mm'/'\evensidemargin -1mm'/&
        '\textwidth 155mm'/'\textheight 220mm'/'\parskip 2mm'/&
        '\parindent 3mm'/'%\pagestyle{empty}'//'\begin{document}'//)
  endif
!
!  stop 'at present'
!--------------------------------------------------------------------
! Now scan the code file to generate a new docfile
  nw=11
  nl=0
  merge=.false.
100 continue
! we are looking for a \begin{verbatim} in source file
  read(lunf90,20,end=900)line
  nl=nl+1
  kink=index(line,'include "')
  if(kink.gt.0) then
     lunf90=lunf90+1
     kend=index(line,'" ')
     curfil=line(kink+9:kend-1)
     write(*,101)line(kink+9:kend-1)
101  format(/' >>> opening include file: ',a/)
     open(lunf90,file=line(kink+9:kend-1),access='sequential',status='old')
     lunlevel=lunlevel+1
  endif
  if(line(1:13).eq.'!\addtotable ') then
! an entry to the table of functions and subroutines
     ntablentries=ntablentries+1
     write(*,*)'Add to table: ',trim(line(14:)),ntablentries
     if(ntablentries.gt.500) then
        write(*,*)'Too many table entries',ntablentries
        goto 990
     endif
     tablentries(ntablentries)=line(14:)
     tablefile(ntablentries)=curfil
     goto 100
  elseif(line(1:18).ne.beginxverb) then
     goto 100
  endif
! we have found a !\begin{verbatim} in source code, remove the "!"
!  write(*,*)'Found a critical part',nl,jl,merge
  newverb=.true.
  if(.not.merge) then
     line=line(2:)
     nbeg=nl
! store text in verbbuff
     jl=0
  else
!     write(*,*)'We have a critical part to merge with previous',jl
     goto 130
  endif
  percentplus=.false.
120 jl=jl+1
     if(jl.gt.500) then
        write(*,*)'Too large verbatim section, source lines: ',nl
!        stop
        goto 990
     endif
! save verbatim text in buffer
     verbbuff(jl)=line
! Here we are reading a critical part until !\end{verbatim}
130  continue
     read(lunf90,20,end=800)line
     if(newverb) then
! save the first line after verbatim 
        lastverb=line
        newverb=.false.
     endif
     nl=nl+1
     if(line(1:16).ne.endxverb) goto 120
     if(line(17:18).eq.'%+') then
        percentplus=.true.
! if there is a %+ after !\end{verbatim} merge with next verbatim
! insert a line ------------ between
        jl=jl+1
        verbbuff(jl)='---------------'
        merge=.true.
        includedverb=includedverb+1
!        write(*,*)' +++ We have to seach for a critical part to add',jl
        goto 100
     else
!        write(*,*)'Total lines in critical part: ',jl
        merge=.false.
     endif
! we should now search for section in LaTeX file where this fits
! The second line must be idential in the LaTeX file (excluding !)
  line=line(2:)
  jl=jl+1
  verbbuff(jl)=line
  if(.not.dokfil) then
! no dokfil just write verbatim text on new dokfile
     do il=1,jl
        write(23,20)verbbuff(il)(1:len_trim(verbbuff(il)))
        nw=nw+1
     enddo
     goto 100
  endif
! search in old dokfil for line matching verbbuff(2)
  rewind(22)
  nd=0
  iskip=0
200 continue
! This is searching in the LaTeX file, nverbsec is last found section
! The new verbatim section must be directly after this
     foundverb=nverbsec
     call searchverb(22,nd,foundverb,sectext,dend)
     if(dend) then
! no matching verbatim text stop and ask user to edit documentation file
        write(*,201)lastokverb,texverb,&
        verbbuff(2)(1:len_trim(verbbuff(2))),curfil
201     format(72('*')/&
             'New verbatim missing in old documentation after line ',i7,/&
             'with verbatim: ',a/&
             'New verbatim text: ',a/'Source file: ',a/72('*')/&
             'Please edit documentation or source and restart'/&
             'Please note OUTPUT FILE IS NOT COMPLETE!!')
!        stop
        goto 990
     endif
     nbegverb=nd
     read(22,20,end=810)dline
     nd=nd+1
     if(dline.ne.verbbuff(2)) then
        write(*,207)trim(verbbuff(2)),trim(dline),trim(texverb),trim(sectext)
207     format(/' *** New verbatim in source code with first line: '/a/&
             10x,'not equal to next verbatim in the documentation: '/a//&
             'Modify documentation after verbatim:'/a/&
             'before the section:'/a/&
             'or reorganize the source code'/&
             'Please note OUTPUT FILE IS NOT COMPLETE!!'/)
!        stop
        goto 990
     endif
! remember the last verbatim section with correct first line
     lastokverb=nd
     texverb=lastverb
     nverbsec=nverbsec+1
! !\begin{verbatim} 
     write(*,217)nverbsec,nd,trim(verbbuff(2)(2:40))
217  format('Matched verbatim ',i4,' at line: ',i6,': ',a,'...')
  once=.true.
  il=2
  nbeg=nd
220 continue
     il=il+1
! dline is read from the OLD LaTeX file
     read(22,20,end=810)dline
     nd=nd+1
! verbbuff is from new source code, should be jl lines
     if(il.le.jl .and. dline.eq.verbbuff(il)) then
        continue
     elseif(il.lt.jl) then
        if(once) then
!           write(*,*)'  *** Some verbatim lines different from old *** '
           once=.false.
        endif
        nw=nw+1
     endif
!     write(*,'(a,2i4,a/a)')'vvvv: ',il,jl,dline(1:50),verbbuff(il)(1:50)
!     verbbuff(il)=dline
     if(il.lt.jl) goto 220
!     if(dline(1:15).ne.'\end{verbatim} ') goto 220
250 continue
     do iverb=1,nverb
        if(nsecverb(3,iverb).eq.nbegverb) goto 270
     enddo
     write(*,266)nbegverb,nverb,(nsecverb(3,i),i=1,nverb)
266  format('Error searching for preceeding section:'/2i5,2x,10i5)
!     stop 'error in scanning old docfile'
     goto 990
270  continue
! mark that section written on new docfile
     nsecverb(3,iverb)=-nsecverb(3,iverb)
! this call copies the section INCLUDING THE OLD VERBATIM
! if once is false the old verbatim is commented away with %
!     write(*,*)'verbatim?',iverb,nsecverb(1,iverb),nsecverb(2,iverb)
     call copytonewdoc(22,23,nsecverb(1,iverb),nsecverb(2,iverb),nw,once)
! if once is .false. then add the new verbatim
     if(.not.once) then
        write(*,*)' *** Replacing with new verbatim, skipping first line '
        do jjj=2,jl
           write(23,'(a)')trim(verbbuff(jjj))
        enddo
! add a } to match {\small and an empty line
        write(23,280)
280     format('}'/)
     endif
!  endif
  goto 100
! all done??
!--------------------------------------
800 continue
  write(*,*)'EOF while searching for !\end{verbatim} in source, line: ',nbeg
  if(lunlevel.gt.1) then
     close(lunf90)
     lunf90=lunf90-1
     lunlevel=lunlevel-1
     goto 130
  else
!     stop
     goto 990
  endif
810 continue
  write(*,*)'EOF while searching for \end{verbatim} in old dokfile, line: ',&
       nbeg
!  stop
  goto 990
!--------------------------------------
! end of file in source file
900 continue
  close(lunf90)
  if(lunlevel.gt.1) then
     lunf90=lunf90-1
     lunlevel=lunlevel-1
! continue reading from source file
     goto 100
  endif
! check if any old docfile parts not written on new
  do iverb=1,nverb
     if(nsecverb(3,iverb).gt.0) then
!        write(*,911)iverb,(nsecverb(j,iverb),j=1,3),trim(verbbuff(iverb))
        write(*,911)iverb,(nsecverb(j,iverb),j=1,3)
911     format('%!%!%!%!%!%!%!%!%! Section with verbatim missing in code',4i5)
        write(23,912)(nsecverb(j,iverb),j=1,3)
912     format('%!%!%!%!%!%!%!%!%! Section with verbatim missing in code',3i5)
        call copytonewdoc(22,23,nsecverb(1,iverb),nsecverb(2,iverb),nw,.true.)
     endif
  enddo
990 continue
! finished or terminating due to error
  if(dokfil) then
     close(22)
  endif
  if(ntablentries.gt.0) then
! write all table entries in alphabetical order
     tabord=0
!     write(*,*)'Sorting table ',ntablentries
     call ssort(tablentries,ntablentries,tabord)
!     write(*,'(a,(15i4))')'Sorted table ',(tabord(jj),jj=1,ntablentries)
     write(23,991)ntablentries
991  format('\newpage'/'Tables with ',i5,' functions and subroutines'//&
          '\begin{tabular}{ll}'/'Name & File \\\hline')
     kk=0
     do jj=1,ntablentries
! table the table entries in alphabetical order
        nounderscore1=tablentries(tabord(jj))
!        write(*,'(a,2i4,2x,a)')'table: ',jj,tabord(jj),trim(nounderscore1)
! replace any _ by \_
        k1=1
980     continue
        k2=index(nounderscore1(k1:),'_')
        if(k2.gt.0) then
           nounderscore1(k1+k2:)=nounderscore1(k1+k2-1:)
           nounderscore1(k1+k2-1:k1+k2-1)='\'
           k1=k1+k2+2
!           write(*,'(a,i3,1x,a)')'underscore',k1,trim(nounderscore1)
           goto 980
        endif
        write(23,992)trim(nounderscore1),trim(tablefile(tabord(jj)))
992     format(a,' & ',a,'\\')
        kk=kk+1
        if(kk.gt.40) then
           write(23,995)
995        format('\end{tabular}'//'\begin{tabular}{ll}'/&
                'Name & File \\\hline')
           kk=0
        endif
     enddo
     write(23,993)
993  format('\end{tabular}')
     write(*,*)'Total number of table entries: ',ntablentries
  endif
  write(23,20)'\end{document}'
  close(23)
  write(*,500)nl,nw,sfil(1:len_trim(sfil))
500 format('Read ',i7' lines, written ',i7,' lines on ',a)
end Program makedok

!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!!/!!\!

subroutine skiptoline(lin,linno)
! rewind lin and then skip lino lines
  character line*80
  rewind(lin)
  nl=0
10 continue
  if(nl+1.eq.linno) goto 1000
  read(lin,20)line
20 format(a)
  nl=nl+1
  goto 10
1000 continue
  return
end subroutine skiptoline

subroutine searchverb(lin,linno,lastfound,sectext,EOF)
! searches for a line with \begin{verbatim} from current line
! return the line number in linno
! skp lines less than lastfound
  character line*80,sectext*(*)
  logical EOF
  iseqbeg=0
  EOF=.false.
10 continue  
  read(lin,20,end=1100)line
20 format(a)
  linno=linno+1
! save last \(sub)section
  if(line(1:8).eq.'\section' .or. line(1:11).eq.'\subsection' .or.&
       line(1:14).eq.'\subsubsection') then
     sectext=line
  endif
  if(line(1:17).ne.'\begin{verbatim} ') goto 10
! skip all lines with verbatim already found
  if(lastfound.ge.0) then
     iseqbeg=iseqbeg+1
     if(iseqbeg.le.lastfound) goto 10
  endif
1000 continue
  return
! no more \begin{verbatim}, just return last line and set EOF true
1100 continue
  EOF=.true.
  goto 1000
end subroutine searchverb

subroutine searchsection(lin,linno,NOBEG,EOF)
! searches for a line with \section{... or \subsection{... or \subsub...
! from current line
! return the line number in linno
  character line*80
  logical EOF,NOBEG
!  write(*,*)'Enter searchsection: ',linno
  lin1=linno
  EOF=.false.
10 continue  
  read(lin,20,end=1100)line
20 format(a)
  linno=linno+1
! ASCII value of \
  if(ichar(line(1:1)).ne.92) goto 10
! when searching for first \begin{verbatim} lin1=0
  if(NOBEG .and. line(1:17).eq.'\begin{verbatim} ') then
     write(*,33)lin1,linno
33   format(' *** WARNING, two verbatim sections with no \{sub}section',&
          ' in between'/' between lines ',i5,' and ',i5)
     stop 'fix LaTeX file and rerun'
  endif
!  write(*,*)'Al line ',linno,line(1:20)
  if(line(1:9).eq.'\section{' .or. line(1:12).eq.'\subsection{' .or.&
       line(1:15).eq.'\subsubsection{' .or. &
       line(1:18).eq.'\subsubsubsection{') goto 1000
  goto 10
1000 continue
  return
! no \section or \subs... set EOF true
1100 continue
  EOF=.true.
  goto 1000
end subroutine searchsection

subroutine copytonewdoc(lin,lut,firstline,lastline,nw,sameverb)
! copies old LaTeX file from firstline to lastline in lin to lut
  integer firstline,lastline
! This line must be long enough to hold a whole paragraph
  character line*256
  logical sameverb,verbic,mark
! sameverb is .false. if old verbatim text should be commented
  verbic=.true.
!  write(*,*)'Enter copytonew',firstline,lastline,nw
  call skiptoline(lin,firstline)
  iz=0
  lc=firstline-1
  mark=.true.
10 continue
  if(lc.ge.lastline) goto 1000
     read(lin,20,end=1100)line
20   format(a)
     k=len_trim(line)
     if(k.gt.250) then
        write(*,21)lc,kc
21      format('Beware, line ,'i4,' longer than 250 characters',i4/&
             'output may be truncated')
     endif
     if(verbic) then
! start 
        if(.not.sameverb) then
           if(index(line,'begin{verbatim}').gt.0) then
              verbic=.false.
           endif
        endif
! do not write a line with "\end{document}"
        if(index(line,'\end{document}').gt.0) then
           write(*,*)'Skipping origial end of document'
        else
           write(lut,20)line(1:k)
        endif
     else
! old verbatim as comments, DO NOT WRITE THE \end{verbatim}, it creates trouble
        if(mark) then
           write(lut,98)
98         format('%! THE LINES BELOW ARE TO BE DELETED WHEN TEXT UPDATED')
99         format('%! THE LINES ABOVE ARE TO BE DELETED WHEN TEXT UPDATED'/&
                '%! THE LINES BELOW ARE FROM THE NEW SOURCE CODE')
           mark=.false.
        endif
        if(line(1:14).ne.'\end{verbatim}') then
           write(lut,23)line(1:k)
23         format('%',a)
!        else
!           write(*,*)'Skipping commented \end{verbatim}'
        endif
     endif
     lc=lc+1
     nw=nw+1
     goto 10
1000 continue
     if(.not.mark) write(lut,99)
  return
! this should not happen
1100 continue
  write(*,1110)lc,lastline
1110 format('EOF when copying old docfile to new ',2i7)
  stop 
end subroutine copytonewdoc
  
subroutine scandoc(lind,nl,nverb,linsec)
! subroutine to scan old doc file and exrtact "verb" sections and
! the appropriate surrounding {sub}section limits.
  character line*80
  dimension linsec(3,*)
  logical EOF
!  write(*,*)'Enter scandoc at line ',nl
  nverb=0
  nsec=nl
10 continue
  read(lind,20,end=900)line
20 format(a)
  nl=nl+1
! ASCII value of \ is 92. EMACS did not like '\'
  if(ichar(line(1:1)).ne.92) goto 10
100 continue
  if(line(1:17).eq.'\begin{verbatim} ') then
! a \begin found, store the line of the section start
     nverb=nverb+1
     nv=nl
!     write(*,*)'Found: ',nl,nverb,line(1:len_trim(line))
     linsec(1,nverb)=nsec
     linsec(3,nverb)=nl-1
!     write(*,*)'Storing start ',nverb,nsec
110  continue
     read(lind,20,end=800)line
     nl=nl+1
     if(line(1:15).eq.'\end{verbatim} ') then
        lend=nl
!        write(*,*)'Calling searchsection ',lend
        call searchsection(lind,lend,.TRUE.,EOF)
        linsec(2,nverb)=lend-2
!        write(*,*)'Line for next section after \end{verbatim} ',lend
        nsec=lend-1
        nl=lend
        if(EOF) goto 1000
     else
        goto 110
     endif
!  else
!     write(*,*)'skipping'
  endif
  goto 10
!------------------------
! EOF here means missing \end{verbatim}
800 continue
  write(*,*)'Missing \end{verbatim} in old doc file, begin at line: ',nv
900 continue
! this makes last part of file will be written as "unmatched"
  nverb=nverb+1
  linsec(1,nverb)=nsec
  linsec(2,nverb)=nl-1
  linsec(3,nverb)=nl-1
1000 continue
  return
end subroutine scandoc

subroutine ssort(CMD,NS,INDEX)
!...SORTING characters max 80 characters
! index returns the alphabetical order of CMD, no change in CMD
  CHARACTER CMD(*)*(*),STR*80
  DIMENSION INDEX(*)
  L=LEN(CMD(1))
  ITOP=1
  INDEX(ITOP)=1
!  write(*,*)'ssort: ',L,NS
100 ITOP=ITOP+1
  IF(ITOP.GT.NS) GOTO 900
  STR=CMD(ITOP)
  IF(STR(1:L).GE.CMD(INDEX(ITOP-1))) THEN
     INDEX(ITOP)=ITOP
     GOTO 100
  ENDIF
!  write(*,*)'find place ',itop
  J1=1
  J2=ITOP
  J=(J1+J2)/2
200 IF(STR(1:L).LT.CMD(INDEX(J))) THEN
     J2=J
  ELSEIF(J.GT.J1) THEN
     J1=J
  ELSE
     J=J2
     GOTO 300
  ENDIF
  IF(J1.NE.J2) THEN
     K=J
     J=(J1+J2)/2
     IF(K.NE.J) GOTO 200
     J=J2
  ENDIF
!...PLACE FOUND
300 CONTINUE
  MOVE: DO K=ITOP-1,J,-1
     INDEX(K+1)=INDEX(K)
  enddo MOVE
  INDEX(J)=ITOP
  GOTO 100
900 RETURN
END SUBROUTINE SSORT

