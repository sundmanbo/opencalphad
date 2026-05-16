program linkocdate
! extract current date and inserts it in the source code of the main program
  character date*8,time*10,mdate*18,line*80
  call date_and_time(date=date,time=time)
  write(*,*)'Stored linking date+time: ',date,' ',time(1:4)
  mdate="'"//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '//&
        time(1:2)//':'//time(3:4)//"'"
  open(21,file='src/pmain1-save.F90',access='sequential',status='old')
  open(22,file='src/pmain1.F90',access='sequential',status='unknown')
100 continue
  read(21,110,end=200)line
  k=index(line,'linkdate=')
  if(k.gt.0) then
     line(k+9:)=mdate
  endif
  write(22,110)line(1:len_trim(line))
110 format(a)
  goto 100
200 continue
  close(21)
  close(22)
end program linkocdate
