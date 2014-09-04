!
! liboctqc.F90
!
! Copyright Bo Sundman, 2013
!
! This is a single Fotrtan subroutine that can be called from C/C++ program
! to access the TQlib subroutines.  In the future one may implement
! several additional subroutines here, this is just a test.
! The function value returned is the error code (0 no error)
! other values can be returned in the arguments.
!
! 140214 BOS modified to new names of libraries and some other things
!
!\begin{verbatim}
integer function liboctqc1(call,iv,tv,rv,iq)
!
! this is the name of the Fortran TQ library
  use liboctq
!
! This is written as a single subroutine because all need acces to the
! ceq variable with the Fortran 90 structure the an equilibrium.  
! If the C program could declare a module this variable could be global
! inside this module.
! It is intersting to have several ceq variables as they are independent
! and can be executed in parallell by separate threads.
! Currently the integer variable ceq identifies the which ceq variable to
! access bu there is no provision yet (from C) to have several ceq.
! In the Fortran main program it works.
!
  implicit none
! call is which service
! iv is an integer argument
! tv is a text argument
! rv is a double argument (can it be an array?)
! iq idetifies the equilibrium (ceq)
  integer call,iv,iq
  double precision rv
  byte tv(*)
!
! contains calls to
! 1 tqini    initiate OC workspapces
! 2 tqrfil   read a TDB file
! 3 tqgcom   get number and names of components
! 4 tqgcom   get number of phases
! 5 tqgcom   get phase name
! 6 tqsetc   set a condition
! 7 tqce     calculate equilibrium
! 8 tqgetv   extract values
! 9          reset error code
! 
!\end{verbatim}
!
  integer i,j,k,ierr,idum,nc,ij,cnum(20)
  character*256 filename,cnames(40)*24
  type(gtp_equilibrium_data), pointer :: ceq
  double precision localv
  double precision, dimension(20) :: npf
! This is the Fortran strucure with all information about phase, conditions etc
! one can have several in the Fortran program but I do not know how to manage
! that i C
  save ceq
!
!  write(*,*)'entering liboctqc1: ',call
  ierr=0
  ij=call
!-------------------------------------------
  select case(ij)
!-------------------------------------------
  case default
! this if call is outside the defined range
     ierr=7777
!-------------------------------------------
  case(1) ! initiate
! idum has no meaning
     idum=20
     call tqini(idum,ceq)
     if(gx%bmperr.ne.0) goto 900
!     write(*,*)'back from tqini',ierr
!-------------------------------------------
  case(2) ! read tdbfile
     filename=' '
     tochar: do i=1,len(filename)
        if(tv(i).eq.0) exit tochar
        filename(i:i)=char(tv(i))
     enddo tochar
     i=len_trim(filename)
!     write(*,*)'Filename: ',i,filename(1:i+10)
!     filename(i+1:)=' '
!     filename='steel1.TDB '
     call tqrfil(filename,ceq)
     if(gx%bmperr.ne.0) goto 900
!-------------------------------------------
  case(3) ! number of components and their names
     do i=1,20
        cnames(i)=' '
     enddo
     call tqgcom(nc,cnames,ceq)
     if(gx%bmperr.ne.0) goto 900
     iv=nc
! tv dimensioned (20,24) in C meaning 24 characters for each element
     do i=1,nc
        j=len_trim(cnames(i))
        do k=1,j
           tv(24*(i-1)+k)=ichar(cnames(i)(k:k))
        enddo
! terminate each name with a zero to please the C program
        tv(24*(i-1)+j+1)=0
     enddo
!-------------------------------------------
  case(4) ! number of phases
     call tqgnp(nc,ceq)
     if(gx%bmperr.ne.0) goto 900
     iv=nc
!-------------------------------------------
  case(5) ! list names of phanse
     cnames(1)=' '
     nc=iv+1
     call tqgpn(nc,cnames(1),ceq)
     if(gx%bmperr.ne.0) goto 900
     j=len_trim(cnames(1))
!     write(*,*)'Name of phase ',nc,': ',cnames(1),j
     do k=1,j
        tv(k)=ichar(cnames(1)(k:k))
     enddo
!     write(*,*)'First letter: ',char(tv(1))
! terminate each name with a zero to please the C program
     tv(j+1)=0
!-------------------------------------------
  case(6) ! set condition
     filename=' '
     tochar2: do i=1,len(filename)
        if(tv(i).eq.0) exit tochar2
        filename(i:i)=char(tv(i))
     enddo tochar2
     i=iv
     j=0
     localv=rv
!     write(*,*)'CASE 6; ',filename(1:len_trim(filename)),i,localv
     call tqsetc(filename,i,j,localv,cnum(1),ceq)
     if(gx%bmperr.ne.0) goto 900
!-------------------------------------------
  case(7) ! calculate equilibrium
     i=0
     j=0
     filename=' '
     localv=0
     call tqce(filename,i,j,localv,ceq)
     if(gx%bmperr.ne.0) goto 900
!-------------------------------------------
  case(8) ! get value
     filename=' '
     tochar3: do i=1,len(filename)
        if(tv(i).eq.0) exit tochar3
        filename(i:i)=char(tv(i))
     enddo tochar3
     i=iv
     j=0
     k=size(npf)
     call tqgetv(filename,i,j,k,npf,ceq)
     if(gx%bmperr.ne.0) goto 900
     rv=npf(1)
!--------------------------------------------
  case(9) ! reset error vcode
     gx%bmperr=0
  end select
  goto 1000
! If there is an error I terminate here as I do not know how to do it in C
900 continue
  write(*,*)'Error code ',gx%bmperr,' for call ',call
  ierr=gx%bmperr
  stop
!
1000 continue
  liboctqc1=ierr
  return
end function liboctqc1
!
