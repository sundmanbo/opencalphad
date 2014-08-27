! program to extract dokumentation of a Fortran source file
!
Program makedok
  character pfil*64,dfil*64,sfil*64,line*80,dline*80
  character beginxverb*18,endxverb*16,ch1*1
  character verbbuff(500)*80
  integer, dimension(3,500) :: nsecverb
  logical dokfil,once,lend,dend,EOF,merge
  beginxverb='!\begin{verbatim} '
  endxverb='!\end{verbatim} '
  write(*,*)ichar('\')
  write(*,10)
!
! Below the idea is to extract the text between 
! 1. the closest \section BEFORE current \verbation and
! 2. the closest BEFORE THE NEXT \verbatim.
! maybe better change to the text
! 1. the closest \section AFTER THE PREVIOUS \verbatim to
! 2. to closest \section AFTER THE CURRENT \verbatim
! The important thing is to keep belong together together
! and not to miss any documentation text
!
! Added that a source file can have (one level) of include files
!
10 format(/'   This is a program to generate and update documentation',&
        ' of software.'//&
        'It expects two files, one with the software code',&
        ' and one with existing'/'LaTeX documentation (this can be empty)',//&
        'In the software the parts that should be included in the',&
        ' documentation'/'must be enclosed by lines !\begin{verbatim}',&
        ' and !\end{verbatim}.'//&
        'If there is "!\end{verbatiom} %+" this verbatim will be merged with',&
        ' the next one.'//&
        'Such a section will be searched for',&
        ' in the documentation file and a match'/&
        'is found with the line after the \begin{verbatim}',&
        ' in the old documentation,'/&
        'all text between',&
        ' the preceeding \{sub}..section up to the',&
        ' \{sub}..section'/&
        'preceeding the next \begin{verbatim},'&
        'will be copied to the new',&
        ' documentation'/'file.'//&
        'If there are differences a message will be written and both old',&
        ' and new'/'verbatim texts written on new documentation file.'//&
        'If no match with second line a message will be generated and the',&
        ' message'/'and the verbatiom part will be written to the new',&
        ' documentation file.'//&
        'The input files will not be changed.'//&
        'Program file name:')
  read(*,20)pfil
  k=index(pfil,'.')
  if(k.eq.0) then
     k=len_trim(pfil)
     pfil(k+1:)='.F90 '
  endif
20 format(a)
  write(*,30)
30 format('Current documentation file name (if none give return):')
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
! source file
  lunlevel=1
  lunf90=31
  open(lunf90,file=pfil,access='sequential',status='old')
! new documentation file
  open(23,file=sfil,access='sequential',status='unknown')
!
37 continue
  if(dokfil) then
     open(22,file=dfil,access='sequential',status='old')
! search for first \begin{verbatim}
     linoldok=0
     dend=.false.
     call searchverb(22,linoldok,dend)
     rewind(22)
     if(dend) then
! no \begin{verbatim} in current dokumentation ... just ignore old dokfile
        write(*,*)'Old documentation file ignored 1'
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
     call copytonewdoc(22,23,i,lastlinsec,nw)
! now scan the whole document 
!     write(*,*)'Calling scandoc at ',nld
     call scandoc(22,nld,nverb,nsecverb)
     write(*,*)nverb,nld
44   format('Found ',i4,' verbatim sections in old documentation, lines',i7)
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
! Now actually write new docfile
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
     write(*,*)'Opening include file: ',line(kink+9:kend-1)
     open(lunf90,file=line(kink+9:kend-1),access='sequential',status='old')
     lunlevel=lunlevel+1
  endif
  if(line(1:18).ne.beginxverb) goto 100
! we have found a !\begin{verbatim} in source, remove the !
  if(.not.merge) then
     line=line(2:)
     nbeg=nl
! store text in verbbuff
     jl=0
  else
     goto 130
  endif
120 jl=jl+1
     if(jl.gt.500) then
        write(*,*)'Too large verbatim section, source line: ',nl
        stop
     endif
! save verbatim text in buffer
     verbbuff(jl)=line
130  continue
     read(lunf90,20,end=800)line
     nl=nl+1
     if(line(1:16).ne.endxverb) goto 120
     if(line(17:18).eq.'%+') then
! if there is a %+ after !\end{verbatim} merge with next verbation
        merge=.true.
        goto 100
     else
        merge=.false.
     endif
  line=line(2:)
  jl=jl+1
  verbbuff(jl)=line
  if(.not.dokfil) then
! no dokfil write verbatim text on new dokfile
     do il=1,jl
        write(23,20)verbbuff(il)(1:len_trim(verbbuff(il)))
        nw=nw+1
     enddo
     goto 100
  endif
! search in old dokfil for matching line after \begin{verbatim
  rewind(22)
  nd=0
200 continue
     call searchverb(22,nd,dend)
     if(dend) then
! no matching verbatim text
        write(*,*)'New verbatim text missing in old, line: ',nw
        write(23,202)
202     format('%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!'/&
             '%! New verbatim text missing in old')
        nw=nw+2
        do il=1,jl
           write(23,222)'new: ',verbbuff(il)(1:len_trim(verbbuff(il)))
           nw=nw+1
        enddo
        write(23,202)
        nw=nw+1
        goto 100
     endif
     nbegverb=nd
     read(22,20,end=810)dline
     nd=nd+1
     if(dline.ne.verbbuff(2)) goto 200
! found matching line after \begin{verbatim} in existing documentation file
!  write(23,20)'\begin{verbatim}'
!  write(23,20)dline(1:len_trim(dline))
!  nw=nw+2
! accept this as match but compare until \end{verbatim}
  once=.true.
  il=2
  nbeg=nd
220 continue
     il=il+1
     if(il.gt.jl) then
        write(*,*)'Old verbatim text longer than new, line: ',nw
        write(23,231)nw
231     format('%! Old verbatim text longer than new line: ',i5)
     endif
!     write(23,20)dline(1:len_trim(dline))
!     nw=nw+1
     read(22,20,end=810)dline
     nd=nd+1
!     write(*,20)dline(1:len_trim(dline))
!     write(*,20)verbbuff(il)(1:len_trim(verbbuff(il)))
     if(il.le.jl .and. dline.eq.verbbuff(il)) then
!        write(23,20)dline(1:len_trim(dline))
!        nw=nw+1
        continue
     elseif(il.lt.jl) then
        if(once) then
           write(*,*)'New verbatim text different from old, line: ',nw
           write(23,236)nw
236        format('%! New verbatim text different from old, line: ',i5)
        endif
        write(*,*)'Old: ',dline(1:len_trim(dline))
        write(*,*)'New: ',verbbuff(il)(1:len_trim(verbbuff(il)))
        write(23,222)'old: ',dline(1:len_trim(dline))
        write(23,222)'new: ',verbbuff(il)(1:len_trim(verbbuff(il)))
222     format('%! ',a,a)
        nw=nw+1
     endif
     if(dline(1:15).ne.'\end{verbatim} ') goto 220
  if(il.lt.jl) then
! new verbatim text longer than old
     write(*,*)'New verbatim text longer than old, line: ',nw
238  format('%! New verbatim text longer than old, line: ',i5)
     write(23,238)nw
     do il=il+1,jl
        write(23,222)verbbuff(il)(1:len_trim(verbbuff(il)))
        nw=nw+1
     enddo
  endif
! continue write from old docfile until next section containing \begin{verb..
!  nextverb=nd
!  call searchverb(22,nextverb,dend)
250 continue
!  write(*,*)'At label 250',nbegverb
!  if(dend) then
! copy to endoffile
!     call copytonewdoc(22,23,nd,nextverb,nw)
!  else
! search for preceeding \section etc
     do iverb=1,nverb
        if(nsecverb(3,iverb).eq.nbegverb) goto 270
     enddo
     write(*,266)nbegverb,nverb,(nsecverb(3,i),i=1,nverb)
266  format('Error searching for preceeding section:'/2i5,2x,10i5)
     stop 'error in scanning old docfile'
270  continue
! mark that section written on new docfile
     nsecverb(3,iverb)=-nsecverb(3,iverb)
     call copytonewdoc(22,23,nsecverb(1,iverb),nsecverb(2,iverb),nw)
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
     stop
  endif
810 continue
  write(*,*)'EOF while searching for \end{verbatim} in old dokfile, line: ',&
       nbeg
  stop
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
! missing check if any old docfile parts not written on new
  do iverb=1,nverb
     if(nsecverb(3,iverb).gt.0) then
911     format('%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!'/&
             '%! Section with unreferenced verbatim ',5i5)
        write(23,911)(nsecverb(j,iverb),j=1,3)
        write(*,911)(nsecverb(j,iverb),j=1,3)
        call copytonewdoc(22,23,nsecverb(1,iverb),nsecverb(2,iverb),nw)
     endif
  enddo
  if(dokfil) then
     close(22)
  endif
  write(23,20)'\end{document}'
  close(23)
  write(*,500)nl,nw,sfil(1:len_trim(sfil))
500 format('Read ',i7' lines, written ',i7,' lines on ',a)
end Program makedok

subroutine copydoc(lin,lut,linno)
! copies text from current position in file lin to file lut
! up to (but not including) the line containing the last \{sub}...section
! before the first \begin{verbatiom}

end subroutine copydoc

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

subroutine searchverb(lin,linno,EOF)
! searches for a line with \begin{verbatim} from current line
! return the line number in linno
  character line*80
  logical EOF
  EOF=.false.
10 continue  
  read(lin,20,end=1100)line
20 format(a)
  linno=linno+1
  if(line(1:17).ne.'\begin{verbatim} ') goto 10
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
     stop 'fix and rerun'
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

subroutine copytonewdoc(lin,lut,firstline,lastline,nw)
! copies from firstline to lastline in lin to lut
  integer firstline,lastline
  character line*80
!  write(*,*)'Enter copytonew',firstline,lastline,nw
  call skiptoline(lin,firstline)
  lc=firstline-1
10 continue
  if(lc.ge.lastline) goto 1000
     read(lin,20,end=1100)line
20   format(a)
     write(lut,20)line(1:len_trim(line))
     lc=lc+1
     nw=nw+1
!     write(*,*)'copy line: ',nw,line(1:40)
     goto 10
1000 continue
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
  write(*,*)'Missing }end{verbatim} in old doc file, begin at line: ',nv
900 continue
! this makes last part of file will be written as "unmatched"
  nverb=nverb+1
  linsec(1,nverb)=nsec
  linsec(2,nverb)=nl-1
  linsec(3,nverb)=nl-1
1000 continue
  return
end subroutine scandoc
