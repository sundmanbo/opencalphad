!
! gtp3XQ for for MQMQA
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
!>     15B. Section: calculate G and other things for MQMQA and Toop/Kohler
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

! removed debug output

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine config_entropy_mqmqa
!\begin{verbatim}
 subroutine config_entropy_mqmqa1(phvar,moded,lokph,tval)
! dummy
   implicit none
   type(gtp_phase_varres), pointer :: phvar
   integer moded,lokph
   double precision tval
! modified arguments for call
   integer ncon
   type(gtp_species) :: sprec
! fq is max number of quads
! fz max number of constituents on a sublattice
! f1s dimension for other arrays ceff1 etc
!   integer, parameter :: fq=99, fz=20, f1=50
   integer, parameter :: fq=20, fz=10, f1s=50
! max allower error in sum ceqf1=1 and ceqf2=1
   double precision, parameter :: ceqferr=1.0D-7
! number of pairs and sublattice fractions
   integer noofpair,ncons1,ncons2
! not needed ....
!   integer loksp,nspel,ielno(10),nextra,ncation
! these are used to as index of species on sublatte 1 (ee,ff) and 2 (gg,hh)
   integer ee,ff,gg,hh
! loop variables
   integer s1,s2,s3,s4,em,c1
! pointer? to mqmqaf record with all fraction records
   type(gtp_mqmqa_var), pointer :: mqf
!   type(gtp_mqmqa_var) :: mqf
! site fractions and amounts
!   double precision yy1(fz),yy2(fz),nn1(fz),nn2(fz)
   double precision yy1(fz),yy2(fz)
! fractions in sublattices
   double precision sum1,sum2,sum3,sum4,half
! contyp(1-4,i) specify sublattice +/- of element and if alone or mixing 2/1
! contyp(5,i) is is the pair index for a quadrupols that is a pair
! contyp(6-7,i) for a pair are species index
! contyp(8-9,i) for a pair are ZERO
! contyp(6-9,i) for other quadrupols are pair indices (2 or 4 indicies)
! contyp(10,i)  should be i ... just as a check
! contyp(11,i)  for a pair is constituent index in sublattice 1
! contyp(12,i)  for a pair is constituent index in sublattice 1
! contyp(11-12,i) for other quadrupoles are zero
   integer em1,em2,em3,em4,mpj
! %pinq(pair) is index in %contyp for a pair
! cridx(pair_index) is the index of corresponding quad %contyp(5,q) is pair
!   integer cridx(f1) REDUNDANT
! Index to the 2-4 sublattice fractions associated with a quad
!  integer fyqix(2,fq),fyqix2(2,fq)
! pair and coord.equiv fractions for pairs in a  quad
   double precision pair(fq),ceqf1(fq),ceqf2(fq)
! test correct way to calculate pair fraction
   double precision cpair(fq),dcpair(fq,fq),cpairsum,dcpairsum(fq),dp(fq,fq)
   double precision spair(fq)
! various factors
   double precision sm1,term,fffy,fff1,fff2
   double precision ffem,ffceq1,ffceq2,once1,once2
! indicate which species that are involved in a quadrupole
!   integer involved(noofsp),stoix1,stoix2
! species in sublatice 1 and 2
   integer nspin(2),eesub,ggsub
! first and second derivatives wrt constituents ...
!   double precision dma1 is coefficent of site fraction in subl 1 for quad 
!   double precision dms1 is sum of coefficents in subl 1 for a quad i
! sum each part separately
   double precision tsub,dvvv(fq,fq),lsub(fq),tend
   double precision ssub,dssub(fq),send,dsend(fq),squad,dsquad(fq)
   double precision d2ssub(fq*(fq+1)/2)
! first index is sublattice constituent, second is quad index
   double precision b1iA(fz,fq),b2iX(fz,fq),b1iAB(fq),b2iXY(fq),sum1AB,sum2XY
! this should give stoichiometry of (species,quad) on the two sublattices
   double precision dmy1(fz,fq),dmy2(fz,fq)
! second derivatives d2xx of site fractions ...
   double precision d2yy1(fz,fq*(fq+1)/2),d2yy2(fz,fq*(fq+1)/2)
   double precision dpair(f1s,fq),dceqf(f1s,fq),yfrac,dummy1,dq1,dq2,dq3
   double precision dyy1(fz,fq),dyy2(fz,fq),dceqf1(f1s,fq),dceqf2(f1s,fq)
   double precision dsm1(fq),d2sm1(fq*(fq-1)/2),dterm(fq),ojoj,alone1,alone2
   character conname*24,endname*24,spname1*24,spname2*24,connames(fq)*24
   double precision endkvot(fq),dendkvot(fq,fq),d2endkvot(fq,fq*(fq+1)/2)
   double precision mulceq(f1s),dmulceq(f1s,fq*(fq+1)/2),divisor
! this is a scaling with total amount of atoms
   double precision invnorm,fqq,pairceq
! quad entropy rewritten ...
   integer pair1,pair2,pair3,pair4,e2,f2,g2,h2
! save here the indices of constituents in sublattices of pairs
! needed for the charge equivalent fractions, ceqf1 and ceqf2
! MAYBE NOT NEEDED when %contyp(11..12,quad) have constituent indices?
   integer eij(2,fq),nomix,all2,q1,s7,s8
! modfied AB/XY loop requires, pq is pair indices, subcon is sublattices indices
! pqq is pair index in %contyp ...
! fq is index in corresponding %constoi
   integer pq(4),pqq(4),sq1(2),sq2(2),fq1(2),fq2(2)
! to avoid adding quadrupols twice
   logical done,ddebug
!
! This is a maybe a reasonable place to initiate csumx for excess parameters?
!   if(allocated(mqmqa_data%csumx)) then
! this is used in calc_mqmqa to skip excess terms with very small fractions.
!      write(*,*)'3XQ Maybe initiating csumx to FALSE'
! initiating here leads to failed convergece, initiated now in calcg_internal
!      mqmqa_data%csumx=.FALSE.
!   endif
   ncon=phlista(lokph)%tnooffr
   ddebug=.FALSE.
!   ddebug=.TRUE.
   if(ddebug) write(*,*)'3X in config_entropy_mqmqa1',lokph,moded,ncon
!   phrec=phlista(lokph)
   invnorm=phvar%abnorm(1)
!   invnorm=one
!   phvar%abnorm(1)=one
!   phvar%abnorm(1)=one
! We should probably update abnorm(2) and (3) also ...
!   phvar%abnorm(2)=invnorm*phvar%abnorm(2)
!   phvar%abnorm(3)=invnorm*phvar%abnorm(3)
!
!   write(*,'(a,i3,1pe12.3)')'3X in MQMQA, version 5: ',ncon,one/invnorm
!
   if(.not.allocated(mqmqa_data%contyp)) then
      write(*,*)'3X MQMQA missing constituent information'
      gx%bmperr=4399; goto 1000
   endif
   if(ncon.ne.mqmqa_data%nconst) then
      write(*,*)'3Xncon, %nconst: ',ncon,mqmqa_data%nconst
      stop '3X constituent problems in mqmqa ...'
   endif
11    format(a,4(F5.2,2x))
!   write(*,*)'3X error return as unfinished'
!   gx%bmperr=4399; goto 1000
!   if(.not.allocated(phvar%mqmqaf%yy1)) then
! THIS MOVED BELOW BUT SHOULD EVENTUALLY BE HERE
! allocating fraction arrayes for use in entropy an excess calculations
!      write(*,*)'3X allocating phase_varres%mqmqaf arrays'
!      allocate(phvar%mqmqaf%yy1(20))
!      allocate(phvar%mqmqaf%yy2(20))
!... add more ...
!   endif
! to avoid typing too much (should mqmqaf be a target? no compiler error)
! problem allocating arrays to this pointer !!!
!   mqf=>phvar%mqmqaf
   do s1=1,ncon
! wow, using phase constituent order to find quad name !! Keep it at present
!      conname=splista(phrec%constitlist(mqmqa_data%contyp(10,s1)))%symbol
   conname=splista(phlista(lokph)%constitlist(mqmqa_data%contyp(10,s1)))%symbol
      connames(s1)=conname
      if(ddebug) write(*,3)s1,(mqmqa_data%contyp(s2,s1),s2=1,14),&
           (mqmqa_data%constoi(s2,s1),s2=1,4),phvar%yfr(s1),trim(conname)
3     format('3X mq:',i2,1x,4i2,1x,i3,1x,4i2,1x,i2,4i3,4F5.1,F5.2,1x,a)
   enddo
   if(ddebug) then
      do s1=1,ncon
         write(*,4)s1,(mqmqa_data%pp(s2,s1),s2=1,4),trim(connames(s1))
4        format('3X pp:',i2,4(F8.5),2x,a)
      enddo
   endif
!   write(*,'(a,20i3)')'3X pinq: ',(mqmqa_data%pinq(s1),s1=1,mqmqa_data%npair)
!   write(*,6)phvar%yfr
6  format('3X y: ',9F7.4)
! maybe use mqf variables?  Need allocation
! local variables can be replaced by those stored in phvar
! local fraction variables and derivatives
   yy1=zero; yy2=zero; pair=zero; ceqf1=zero; ceqf2=zero;
   b1iA=zero; b2iX=zero; dpair=zero; dceqf1=zero; dceqf2=zero
   b1iAB=zero; b2iXY=zero; dmy1=zero; dmy2=zero
   cpair=zero; dcpair=zero
! write(*,431)'3X d2S/Rx:',((phvar%d2gval(ixsym(s2,s3),1),s3=s2,ncon),s2=1,ncon)
! any species used below is indicated by a 1 or 2 depending on sublattice
!   do s1=1,ncon
!      if(mqmqa_data%contyp(5,s1).ne.s1) then
!         write(*,*)'3X *** Warning %contyp index 10 not correct'
!      endif
!   enddo
! these count the sum of element and pair stoichiometries for a quad
!   fyp=zero
!----------------------------------------------------
! the array species in each sublattice will have missing values
!----------------------------------------------------
! we must calculate a number of auxilliary fraction variables from the
! site fractions using mqmqa_data%contyp
!   do s1=1,ncon
!      write(*,14)'3X %%contyp: ',s1,(mqmqa_data%contyp(s7,s1),s7=1,14),&
!           trim(connames(s1))
14    format(a,i3,1x,4i2,1x,i3,1x,4i3,1x,i3,1x,4i3,1x,a)
!   enddo
   mpj=mqmqa_data%npair
   if(ddebug) write(*,15)'3X pinq: ',mpj,(mqmqa_data%pinq(s1),s1=1,mpj)
15 format(a,i3,2x,20i4)
   nspin(1)=mqmqa_data%ncon1
   nspin(2)=mqmqa_data%ncon2
!   noofpair=mqmqa_data%npair
!   write(*,'(a,2i3,2x,i3)')'3X subl const and pairs: ',nspin,mpj
   noofpair=0
! BIG LOOP OVER ALL QUADS, calculating fracions of pairs, sublattices etc
   sumfrac: do s1=1,ncon
      conname=connames(s1)
      if(mqmqa_data%contyp(10,s1).ne.s1) then
         write(*,*)'3X mqmqa_data%contyp(10,s1) not redundant'
      endif
      s3=mqmqa_data%contyp(5,s1)
      typ: if(s3.gt.0) then
! AN PAIR quadrupol AA:XX, increment the pair counter
! the index of the quadrupole fraction is in ALSO in contyp(10,s1) 
!         yfrac=phvar%yfr(mqmqa_data%contyp(10,s1))
         yfrac=phvar%yfr(s1)
! Pair fractions has to be normallized later, here multiply with %pp
         noofpair=noofpair+1
         pair(s3)=pair(s3)+yfrac
         dpair(s3,s1)=one
         cpair(s3)=cpair(s3)+yfrac*mqmqa_data%pp(1,s1)
! dcpair( pairindex, quadindex )
         dcpair(s3,s1)=mqmqa_data%pp(1,s1)
! Calculating pairs
!         write(*,'(a,2i3,7F8.4)')'3X pair1: ',s1,s3,one,yfrac,pair(s3),&
!              mqmqa_data%pp(1,s1),cpair(s3),dcpair(s3,s1),dpair(s3,s1)
! all second derivatives of pair is zero
! the index of constituent in first sublattice is in %contyp(11,s1)
! ee is species index, eesub is index of species in sublattice
         ee=mqmqa_data%contyp(6,s1)
         eesub=mqmqa_data%contyp(11,s1)
         eij(1,noofpair)=eesub
! the index of constituent in second sublattice is in %contyp(12,s1)
         gg=mqmqa_data%contyp(7,s1)
         ggsub=mqmqa_data%contyp(12,s1)
         eij(2,noofpair)=ggsub
! ee and gg are pair indices, eesub and ggsub sublatice const. indices
!         write(*,50)'3X decode1: ',s1,ee,gg,eesub,ggsub
50       format(a,i3,5x,2i4,5x,2i4)
         spname1=splista(ee)%symbol
         spname2=splista(gg)%symbol
! remember which species that are used by marking them (only needed for pairs)
! this is the stoichiometric factors of the species in the pair
         fff1=2.0d0/mqmqa_data%constoi(1,s1)
         fff2=2.0d0/mqmqa_data%constoi(2,s1)
!         else
!            write(*,*)'3X contyp error 1: ',mqmqa_data%contyp(1,s1)
!            gx%bmperr=4399; goto 1000
!         endif
! SAVE the location in sublattice array of species eesub in quad s1
! eesub and ggsub are sublattice indices
! >>>>>>>>>>>>>>>>>>>>>>> .............. EQUATION B15 part 1
         
         yy1(eesub)=yy1(eesub)+fff1*yfrac
         b1iA(eesub,s1)=fff1
!         write(*,12)'3X yy1 add1:',s1,eesub,1,yy1(eesub),fff1,yfrac
12       format(a,3i3,5F10.6)
! there is a single contribution from this quad to the site fractions
         b1iAB(s1)=fff1
         yy2(-ggsub)=yy2(-ggsub)+fff2*yfrac
         b2iX(-ggsub,s1)=fff2
         b2iXY(s1)=fff2
!         write(*,12)'3X yy2 add1:',s1,ggsub,2,yy2(-ggsub),fff2,yfrac
! equivalent sublattice fraction for the sublattice constituents
! >>>>>>>>>>>>>>>>>>>>>>>> .............. EQUATION B19 part 1
         ceqf1(eesub)=ceqf1(eesub)+yfrac
         ceqf2(-ggsub)=ceqf2(-ggsub)+yfrac
         dceqf1(eesub,s1)=one
         dceqf2(-ggsub,s1)=one
! Calculating ceq
!         write(*,333)'3X ceqf1e:',s1,0,0,1,eesub,ceqf1(eesub),&
!              yfrac,one,1,trim(spname1),trim(connames(s1))
!         write(*,333)'3X ceqf2e:',s1,0,0,2,ggsub,ceqf2(-ggsub),&
!              yfrac,one,2,trim(spname2),trim(connames(s1))
333      format(a,1x,5i3,3F10.6,' ceq',i1,'(',a,')  ',a)
! end of pair summations
      else
!--------------------------------------------------------------
! this is a quadrupol AB:XY consisting of 2 or 4 pairs typ A:X and B:Y
! the pair indices in %contyp are indicated in contyp(6..9,s1)
! IT IS A BIT INVOLVED AND CAN (certainly) BE SIMPLIFIED ....
         ffem=0.5D0
         fffy=one
         yfrac=phvar%yfr(s1)
         if(mqmqa_data%contyp(9,s1).gt.0) then
! contyp(9,s1) nonzero for quadrupoles with 4 pairs A:X, A:Y, B:X, B:Y
! set ffem=0.25 if 4 pairs
            ffem=0.25D0
! set fffy=0.5 to avoid adding same fraction twice
            fffy=0.5D0
         endif
! these refer to constituent species, ff, gg in first; gg hh in second
!         ee=0; ff=0; gg=0; hh=0
! s2 loops over the species involved in the quadrupol, it can be 3 or 4
! in %contyp(1..4,s1) is indicated if same species twice (2) or not (1)
! in %constoi(1..4,s1) is the coordination number
! s2 loops positions 1..4 in contyp and constoi
! these are used to find correct stoichimetry index
! which constoi to use? AA:XY should have (1,2) and (1,3) for AA:XX and AA:YY
! which constoi to use? AB:XX should have (1,3) and (2,3) for AA:XX and BB:XX
! which constoi to use? AB:XY should have (1,3), (1,4), (2,3) and (2,4) for ...
! position 6, 7, 8, 9 are indices to pairs, s2 incremented at loop end
! in the pairs %contyp(11,pairindex) and %contyp(12,pairindex) and subl index
         once1=one; once2=one; alone1=2.0d0; alone2=2.0d0
         ffceq1=0.5D0; ffceq2=0.5D0
         pq=0; sq1=0; sq2=0
! pq are the pair indices, 2 or 4
! but below we use pq as indices to mqmqa_data ... we need pinq(pq(j))
         pq(1)=mqmqa_data%contyp(6,s1)
         pq(2)=mqmqa_data%contyp(7,s1)
         pqq(1)=mqmqa_data%pinq(pq(2))
         pqq(2)=mqmqa_data%pinq(pq(2))
! here we saved A and X assuming mixing in first sublattice
! we must also save the stoichiometric factors of the sublattice species
         sq1(1)=mqmqa_data%contyp(11,pqq(1))
! fq1 this is index to %constoi for this sublattice constituent
         fq1(1)=1
         sq2(1)=mqmqa_data%contyp(12,pqq(1))
         fq2(1)=3
         if((mqmqa_data%contyp(1,s1).eq.2)) then
! quadruplet AA:XY, pairs AA:XX and AA:YY
! Same constituents in first sublattice, indices in %contyp(11, %contyp(6,s1))
!                                               and %contyp(12, %contyp(7,s1))
! mixing in second sublattice, same constituent twice in first
            sq1(2)=sq1(1)
            fq1(2)=fq1(1)
! replace stoichiometric factor
            fq2(1)=2
            sq2(2)=mqmqa_data%contyp(12,pqq(2))
            fq2(2)=3
            alone2=one
            nomix=1
!            write(*,'(a,2i3,2x,2i3)')'3X mixing in 2: ',sq1,sq2
         elseif(abs(mqmqa_data%contyp(3,s1)).eq.2) then
! quadrupole AB:XX, first pair AA:XX, second BB:XX
! Same constituents in second sublattice, indices in %contyp(11, %contyp(6,s1))
!                                                and %contyp(12, %contyp(7,s1))
! mixing in first sublattice, same constituent twice in second
            sq2(2)=sq2(1)
            fq2(2)=fq2(1)
! add second sublattice constituent twice
            sq1(2)=mqmqa_data%contyp(11,pqq(2))
            fq1(2)=2
            alone1=one
            nomix=2
!            write(*,'(a,2i3,2x,2i3)')'3X mixing in 1: ',sq1,sq2
         else
! quadupole AB:XY, 4 pairs used, AA:XX; AA:YY; BB:XX BB:YY
! 4 pairs, we have to add 2 more
            pq(3)=mqmqa_data%contyp(8,s1)
            pq(4)=mqmqa_data%contyp(9,s1)
            pqq(3)=mqmqa_data%pinq(pq(3))
            pqq(4)=mqmqa_data%pinq(pq(4))
            sq1(2)=mqmqa_data%contyp(11,pqq(3))
! 
            fq1(2)=2
            fq2(2)=4
! I am not sure how the pairs are arranged 
! but testing 3 pairs the sublattice constituent must be different
            if(sq1(2).eq.sq1(1)) sq1(2)=mqmqa_data%contyp(11,pqq(2))
            sq2(2)=mqmqa_data%contyp(12,pqq(2))
            if(sq2(2).eq.sq2(1)) sq2(2)=mqmqa_data%contyp(12,pqq(2))
            alone1=one; alone2=one
            nomix=4
!            write(*,*)'3X reciprocal cluster',mqmqa_data%contyp(2,s1)
         endif
!         write(*,313)'3X pq mm: ',s1,pq,pqq,sq1,sq2,fq1,fq2
313      format(a,i3,2x,4i2,2x,4i2,4x,2i2,2x,2i2,4x,2i2,2x,2i2)
! contribution from all pairs included in this quadruple, nonzero pq
         pqloop: do s2=1,4
!            write(*,'(a,2i3)')'3x pqloop: ',s2,pq(s2)
            if(pq(s2).eq.0) exit pqloop
            pair(pq(s2))=pair(pq(s2))+ffem*yfrac
            dpair(pq(s2),s1)=ffem
! EMERGENCY, how to know which %pp to use for each pair???
! modified in gtp3B to ensure that pairs are correlated with %constoi ??
! s2 is assumed to be %pp index, pq(s2) constittuent index ...
            cpair(pq(s2))=cpair(pq(s2))+yfrac*mqmqa_data%pp(s2,s1)
! dcpair( pairindex, quadindex )
            dcpair(pq(s2),s1)=mqmqa_data%pp(s2,s1)
!            write(*,'(a,3i3,2F10.7)')'3X dcpair2: ',pq(s2),s1,s2,&
!                 yfrac,dcpair(pq(s2),s1)
! Calculating pairs in SNN
!            write(*,'(a,3i3,6F10.6)')'3X pair2: ',s1,s2,pq(s2),ffem,yfrac,&
!                 pair(pq(s2)),mqmqa_data%pp(s2,s1),cpair(pq(s2)),&
!                 dcpair(pq(s2),s1)
         enddo pqloop
!         write(*,'(a,i3,2x,2i3,2x,2i3)')'3X sqi: ',s1,sq1,sq2
         s7=0
         subloop: do s2=1,2
! For the site fractions and equivalent fraction ceqfi we have to
! extract all constituent species of the quadrupol s1 using the pair s3
! divided by with the coordination factor in s2 for quadrupol s1
! the species in first sublattice of the pair
            if(sq1(s2).le.0) then
               write(*,*)'3X no constituent in first sublattice!!!',s1,s2,sq1
               stop
            else
! We have to use the correct sublattice index and coordination factor !!
! eesub should be in mqmqa_data%contyp(10+s2,s1) ??  What is sq1(s2)?
               eesub=sq1(s2)
               eesub=mqmqa_data%contyp(10+s2,s1)
!                    write(*,'(a,3i3,F8.3)')'3X sublattice index: ',&
!                    eesub,sq1(s2),fq1(s2),mqmqa_data%constoi(s2,s1)
! SAVE the sublattice location of species eesub for quad s1
               fff1=fffy*alone1/mqmqa_data%constoi(fq1(s2),s1)
               yy1(eesub)=yy1(eesub)+fff1*yfrac
               b1iA(eesub,s1)=fff1
!               write(*,13)'3X yy1 add2:',s1,s2,eesub,yy1(eesub),fff1*yfrac,&
!                    fff1,yfrac,fffy,alone1,mqmqa_data%constoi(fq1(s2),s1)
13             format(a,3i3,3F10.6,5(F6.3))
! there can be more than one contribution to site fraction from this quad
! nomix=1 if single in 1
               if(nomix.ne.1) then
                  b1iAB(s1)=b1iAB(s1)+fff1
               else
                  b1iAB(s1)=fff1
               endif
               ceqf1(eesub)=ceqf1(eesub)+fffy*ffceq1*yfrac
               dceqf1(eesub,s1)=fffy*ffceq1
            endif
!---------- second sublattice
            if(sq2(s2).gt.0) then
! constituent index is negative in second sublattice!!
               write(*,*)'3X no constituent in second sublattice!!!',s1,s2,sq2
               write(*,14)'3X %contyp: ',s1,(mqmqa_data%contyp(s7,s1),s7=1,14)
               gx%bmperr=4399; goto 1000
            else
! NOW the species in second sublattice of the pair NOTE negative
               ggsub=sq2(s2)
! SAVE the sublattice location of species eesub  and ggsub for quad s1
! fq1(s2) specify stoichiometry index of const. in 1st sublattice in AB/XY
! the species indexing in %contyp(11..14) is the same as for %constoi(1..4)
! fq2(s2) specify stoichiometry index of const. in 2nd sublattice in AB/XY
               fff2=fffy*alone2/mqmqa_data%constoi(fq2(s2),s1)
               yy2(-ggsub)=yy2(-ggsub)+fff2*yfrac
!               write(*,13)'3X yy2 add2:',s1,ggsub,s2,yy2(-ggsub),fff2,yfrac,&
!                    fffy,alone2,mqmqa_data%constoi(fq2(s2),s1)
52             format(a,i3,5F8.5)
               b2iX(-ggsub,s1)=fff2
! nomix=2 if single in sublattice 2
               if(nomix.ne.2) then
                  b2iXY(s1)=b2iXY(s1)+fff2
               else
                  b2iXY(s1)=fff2
               endif
331         format('3Xq n(',a2,'): ',3i3,2i3,4F7.4,2x,a)
! equivalent site fraction, each mixing element will be counted twice
! for quadrupole with 4 pairs fffy=0.25; otherwice 0.5
! >>>>>>>>>>>>>>>>>>> ................ EQUATION B17 part 2
               ceqf2(-ggsub)=ceqf2(-ggsub)+fffy*ffceq2*yfrac
               dceqf2(-ggsub,s1)=fffy*ffceq2
!               write(*,333)'3X ceqf1q:',s1,0,s3,1,eesub,ceqf1(eesub),&
!                    yfrac,fffy*ffceq1,1,trim(spname1),trim(connames(s1))
!               write(*,333)'3X ceqf2q:',s1,0,s3,2,ggsub,ceqf2(-ggsub),&
!                    yfrac,fffy*ffceq2,2,trim(spname2),trim(connames(s1))
! increment s2 for next pair in quadrupole s1
               endif
            enddo subloop
         endif typ
! problem with pair fractions ...
!         do s3=1,mpj
!            write(*,'(a,2i3,5F10.7)')'3X loop:',s1,s3,(dpair(s3,s2),s2=1,ncon)
!            write(*,'(a,2i3,5F10.7)')'3X loop:',s1,s3,(dcpair(s3,s2),s2=1,ncon)
!         enddo
      enddo sumfrac
!      write(*,*)'3X sumfrac done'
!
!------------------------------ end BIG LOOP over all quads
!      do s3=1,nspin(1)
!         write(*,342)'3X b1iA(m,n):',s3,s1,(b1iA(s3,s4),s4=1,ncon)
!      enddo
!      write(*,341)'3X b1iAB(n)    :',s1,(b1iAB(s4),s4=1,ncon)
!      do s3=1,nspin(2)
!         write(*,342)'3X b2iX(m,n):',s3,s1,(b2iX(s3,s4),s4=1,ncon)
!      enddo
!      write(*,341)'3X b2iXY(n)    :',s1,(b2iXY(s4),s4=1,ncon)
!      write(*,340)'3X yy1: ',(yy1(s4),s4=1,3)
!      write(*,340)'3X yy2: ',(yy2(s4),s4=1,3)
340   format(a,7F10.7)
342   format(a,2i2,7F7.4)
720   format(a,i3,4(4I3,2x))

! debug listings:
!   write(*,*)'3X summed all amounts, next normallize'
!   write(*,720)'3X contyp:  ',0,((mqmqa_data%contyp(s2,s1),s2=11,14),s1=1,ncon)
!      write(*,200)'3X p_AB/XY:',(phvar%yfr(s1),s1=1,ncon)
!      write(*,200)'3X n1     :',(yy1(s1),s1=1,nspin(1))
!      write(*,200)'3X n2     :',(yy2(s1),s1=1,nspin(2))
!      write(*,200)'3X pairs  :',(pair(s1),s1=1,noofpair)
!      write(*,200)'3X cpairs :',(cpair(s1),s1=1,noofpair)
!      do s1=1,noofpair
!         write(*,200)'3X dcpairs:',(dcpair(s1,s2),s2=1,ncon)
!      enddo
!      write(*,200)'3X ceqf1  :',(ceqf1(s1),s1=1,nspin(1))
!      write(*,200)'3X ceqf2  :',(ceqf2(s1),s1=1,nspin(2))
!   stop
!   do s3=1,nspin(1)
!      write(*,342)'3X b1iA(m,n):',s3,s1,(b1iA(s3,s4),s4=1,ncon)
!   enddo
!   write(*,341)'3X b1iAB(n)    :',s1,(b1iAB(s4),s4=1,ncon)
!   do s3=1,nspin(2)
!      write(*,342)'3X b2iX(m,n):',s3,s1,(b2iX(s3,s4),s4=1,ncon)
!   enddo
!   write(*,341)'3X b2iXY(n)    :',s1,(b2iXY(s4),s4=1,ncon)
341 format(a,i2,7F10.7)
!-------------- we have extracted all comp.variables and their deriv wrt quads
! Now sum amounts and normallize
!
! NOTE in b1iA and b1iA the first index is subl.const, second is quad 
!    sometimes I mix them up ...
!
!   write(*,*)'Sublattice fractions and detivatives:
! first sublattice
   sum1AB=zero
!   write(*,*)'3X nspin: ',nspin
   do s1=1,nspin(1)
      sum1AB=sum1AB+yy1(s1)
!      write(*,88)'3X subl: ',s1,yy1(s1),(b1iA(s1,s2),s2=1,ncon)
   enddo
88 format(a,i2,F7.3,2x,9(F8.4))
!   write(*,'(a,F7.3,a)')'3X sum1AB: ',sum1AB
   do s1=1,nspin(1)
      yy1(s1)=yy1(s1)/sum1AB
   enddo
! second sublattice
   sum2XY=zero
   do s1=1,nspin(2)
      sum2XY=sum2XY+yy2(s1)
!      write(*,88)'3X sub2: ',s1,yy2(s1),(b2iX(s1,s2),s2=1,ncon)
   enddo
   do s1=1,nspin(2)
      yy2(s1)=yy2(s1)/sum2XY
   enddo
!   write(*,*)'3X nspin2: ',nspin
! derivatives of sublattice fractions wrt quads
   all2=ncon*(ncon+1)/2
   d2yy1=zero
   dummy1=one/sum1AB**2
!
!   write(*,*)'3X d2yy1 size: ',fz,fq*(fq+1)/2,fz*fq*(fq+1)/2,all2
      yder1: do s1=1,nspin(1)
      do s2=1,ncon
! b1iAB may contain contributions from two constituents in same quad
         dyy1(s1,s2)=(b1iA(s1,s2)-yy1(s1)*b1iAB(s2))/sum1AB
!         cycle yder1
! this gives phase matrix singuler
         do s3=1,ncon
            d2yy1(s1,ixsym(s2,s3))=&
                 (-b1iA(s1,s2)*b1iAB(s3)-b1iA(s1,s3)*b1iAB(s2)+&
                 2.0D0*yy1(s1)*b1iAB(s2)*b1iAB(s3))*dummy1
!            write(*,19)'3X dyy: ',s1,s2,s3,b1iA(s1,s2),b1iAB(s3),&
!                 b1iA(s1,s3),b1iAB(s2),2.0D0*yy1(s1),d2yy1(s1,ixsym(s2,s3))
19          format(a,3i2,6(1pe10.2))
         enddo
! try ... gives also phase matrix singular ...
!         d2yy1(s1,s1)=one/yy1(s1)
      enddo
   enddo yder1
! debug
!   do s1=1,nspin(1)
!      do s3=s1,all2
!         s8=ixsym(s3,s1)
!        write(*,'(a,4i4,1pe12.4)')'3X mqmqa d2yy1: ',s1,s3,s8,all2,d2yy1(s1,s8)
!         write(*,87)'3X d2yyj: ',1,s1,(d2yy1(s1,s2),s2=1,all2)
!      enddo
!   enddo
87 format(a,2i3,6(1pe10.2))
   d2yy2=zero
   dummy1=one/sum1AB**2
   yder2: do s1=1,nspin(2)
! the line below works when there are no SRO quads (species)
!      dyy2(s1,s1)=one; cycle yder2
! below needed when yy2 calculated from quads
      do s2=1,ncon
! b2iXY may contain contributions from two constituents in same quads
         dyy2(s1,s2)=(b2iX(s1,s2)-yy2(s1)*b2iXY(s2))/sum2XY
         cycle yder2
         do s3=1,ncon
            if(nspin(2).eq.1) then
! single sublattice fractions should not have any second derivaties ??
               d2yy2(s1,ixsym(s2,s3))=zero
            else
! appoximate ...
               d2yy2(s1,ixsym(s2,s3))=&
                    (-b2iX(s1,s2)*b2iXY(s3)-b2iX(s1,s3)*b2iXY(s2)+&
                    2.0D0*yy2(s1)*b2iXY(s2)*b2iXY(s3))*dummy1
            endif
         enddo
      enddo
   enddo yder2
!   do s1=1,nspin(2)
!      write(*,87)'3X d2yyj: ',2,s1,(d2yy2(s1,s2),s2=1,all2)
!   enddo
! ------------------------------------------
! calculate sublattice sites related to formula units
!   dummy1=invnorm/(sum1AB+sum2XY)
!   sum1AB=sum1AB*dummy1
!   sum2XY=sum2XY*dummy1
!   sum1AB=invnorm*sum1AB
!   sum2XY=invnorm*sum2XY
! We have to sum and normalize cpair
   cpairsum=zero
   dcpairsum=zero
   dp=zero
   do s1=1,noofpair
      spair(s1)=cpair(s1)
      cpairsum=cpairsum+cpair(s1)
      do s2=1,ncon
         dcpairsum(s2)=dcpairsum(s2)+dcpair(s1,s2)
         dp(s1,s2)=dcpair(s1,s2)
      enddo
   enddo
!   write(*,'(a,F10.6,2x,10(F8.4))')'3X cpsum:',cpairsum,&
!        (dcpairsum(s2),s2=1,ncon)
   do s1=1,noofpair
      cpair(s1)=cpair(s1)/cpairsum
      do s2=1,ncon
         dcpair(s1,s2)=(cpairsum*dp(s1,s2)-spair(s1)*dcpairsum(s2))/cpairsum**2
      enddo
! replacing pair here creates problems .... do it later
!      pair(s1)=cpair(s1)
! Calculate derivatives of pairs wrt quads, NEEDED FOR REFERENCE STATE
   enddo
!   do s1=1,noofpair
!      write(*,119)'3X cpair: ',s1,cpair(s1),(dcpair(s1,s2),s2=1,ncon)
!   enddo
119 format(a,i2,F10.7,2x,8F10.6)
!
!   check pairs are unity ... this pair fraction is wrong anyway ...
!   write(*,*)'3X pair fractions and derivatives:'
   dummy1=zero
! loop over all pairs
   do s1=1,noofpair
! Check sum is unity
      dummy1=dummy1+pair(s1)
!      write(*,120)s1,pair(s1),(dpair(s1,s2),s2=1,ncon)
   enddo
120 format('3X pairs:',i3,F7.4,1x,10F6.3)
   if(abs(dummy1-one).gt.1.0D-12) then
      write(*,*)'3X pair fractions does not add up to unity',dummy1
      write(*,'(a,10F7.4)')'3X pf: ',(pair(s1),s1=1,noofpair)
      gx%bmperr=4399; goto 1000
   endif
!
! NOW list the Charge Equivalent Fractions, related to sublattices
!   write(*,*)'3X Charge Equivalent fractions and derivatives:'
   dummy1=zero
   do s1=1,nspin(1)
! Check sum is unity
      dummy1=dummy1+ceqf1(s1)
!      write(*,81)'3X ceqf:',1,ceqf1(s1),(dceqf1(s1,s2),s2=1,ncon)
   enddo
   if(abs(dummy1-one).gt.ceqferr) then
      write(*,*)'3X Sum of charge equivalent fractions on subl 1 not 1:',dummy1
      write(*,'(a,7(F10.7))')'3X ceqf1: ',(ceqf1(s2),s2=1,nspin(1))
! assume this will be the fixed when converged ....
!      gx%bmperr=4399; goto 1000
   endif
   dummy1=zero
   do s1=1,nspin(2)
! Check sum is unity
      dummy1=dummy1+ceqf2(s1)
!      write(*,81)'3X ceqf:',2,ceqf2(s1),(dceqf2(s1,s2),s2=1,ncon)
   enddo
   if(abs(dummy1-one).gt.ceqferr) then
      write(*,*)'3X Sum of charge equivalent fractions on subl 2 not 1',dummy1
      write(*,'(a,7(F10.7))')'3X ceqf2: ',(ceqf2(s2),s2=1,nspin(2))
! assume this will be the fixed when converged ....
!      gx%bmperr=4399; goto 1000
   endif
81 format(a,i2,F7.4,1x,(10F7.4))
!   write(*,*)'3X all normallized fractions calculated'
!   write(*,*)'3X error return as unfinished'
!   gx%bmperr=4399
!   goto 1000
!---------------------------------------------------------------------------
! 2021.08.24 derivatives of site fractions wrt quadrupoles??
!---------------------------------------------------------------------------
!   write(*,*)'3X quitting as not finished below'
!   gx%bmperr=4399
!   goto 1000
! fraction listings
!   write(*,200)'3X p_AB/XY:',(phvar%yfr(s1),s1=1,ncon)
!   write(*,200)'3X sites/FU  :',sum1AB,sum2XY
!   write(*,200)'3X y1     :',(yy1(s1),s1=1,nspin(1))
!   write(*,200)'3X y2     :',(yy2(s1),s1=1,nspin(2))
!   do s1=1,nspin(1)
!      write(*,202)'3X dy1/dpi:',s1,(dyy1(s1,s2),s2=1,ncon)
!   enddo
!   do s1=1,nspin(2)
!      write(*,202)'3X dy2/dpi:',s1,(dyy2(s1,s2),s2=1,ncon)
!   enddo
! same as above
!   write(*,200)'3X x_A/B  :',(pair(s1),s1=1,noofpair)
!   write(*,200)'3X ceqf1  :',(ceqf1(s1),s1=1,nspin(1))
!   write(*,200)'3X ceqf2  :',(ceqf2(s1),s1=1,nspin(2))
200 format(a,(10F7.4))
202 format(a,i2,(10F7.4))
!   write(*,*)'3X now the entropy: >>>>>>>>>>>>>'
!--------------------------------------------------------------------------
   ! Problems here!!
! COPY ALL FRACTIONS VARIABLES AND DERIVATIVES TO MQMQAF for use in parameters
! allocate all arrays
   if(.not.allocated(phvar%mqmqaf%yy1)) then
! allocate first time only!!
! mqf is phvar%mqmqaf
      phvar%mqmqaf%nquad=ncon; phvar%mqmqaf%npair=noofpair; 
      phvar%mqmqaf%ns1=nspin(1); phvar%mqmqaf%ns2=nspin(2)
!      write(*,*)'3X allocating phvar%mqmqaf arrays',nspin,ncon,noofpair
      allocate(phvar%mqmqaf%yy1(nspin(1)))
      allocate(phvar%mqmqaf%yy2(nspin(2)))
      allocate(phvar%mqmqaf%dyy1(nspin(1),ncon))
      allocate(phvar%mqmqaf%dyy2(nspin(2),ncon))
      allocate(phvar%mqmqaf%d2yy1(nspin(1),ncon*(ncon+1)/2))
      allocate(phvar%mqmqaf%d2yy2(nspin(2),ncon*(ncon+1)/2))
      allocate(phvar%mqmqaf%ceqf1(nspin(1)))
      allocate(phvar%mqmqaf%ceqf2(nspin(2)))
      allocate(phvar%mqmqaf%dceqf1(nspin(1),ncon))
      allocate(phvar%mqmqaf%dceqf2(nspin(2),ncon))
      allocate(phvar%mqmqaf%pair(noofpair))
      allocate(phvar%mqmqaf%dpair(noofpair,ncon))
!   else
!      write(*,*)'3X copying data to phvar%mqmqaf arrays'
   endif
!   write(*,*)'3X mqf arrays allocated'
!   mqf=>phvar%mqmqaf
!
!   write(*,*)'3X d2yy1: ',nspin(1),all2,size(phvar%mqmqaf%d2yy1)
!   write(*,*)'3X d2yy1: ',nspin(1),all2,nspin(1)*all2
   phvar%mqmqaf%yy1(1)=yy1(1)
   do s1=1,nspin(1)
     phvar%mqmqaf%yy1(s1)=yy1(s1)
     phvar%mqmqaf%ceqf1(s1)=ceqf1(s1)
      do s2=1,ncon
         phvar%mqmqaf%dyy1(s1,s2)=dyy1(s1,s2)
         phvar%mqmqaf%dceqf1(s1,s2)=dceqf1(s1,s2)
      enddo
      do s3=s1,ncon
         s8=ixsym(s3,s1)
!         write(*,'(a,2i3,3i4)')'3X mqmqa: ',s1,s3,s8,ncon*(ncon+1)/2,all2
!         phvar%mqmqaf%d2yy1(s1,ixsym(s3,s2))=d2yy1(s1,ixsym(s3,s2))
! This statement kills whole subroutine
         phvar%mqmqaf%d2yy1(s1,s8)=d2yy1(s1,s8)
      enddo
  enddo
!
   do s1=1,nspin(2)
      phvar%mqmqaf%yy2(s1)=yy2(s1)
      phvar%mqmqaf%ceqf2(s1)=ceqf2(s1)
      do s2=1,ncon
         phvar%mqmqaf%dyy2(s1,s2)=dyy2(s1,s2)
         phvar%mqmqaf%dceqf2(s1,s2)=dceqf2(s1,s2)
      enddo
      do s3=s1,ncon
         s8=ixsym(s3,s1)
         phvar%mqmqaf%d2yy2(s1,s8)=d2yy2(s1,ixsym(s1,s8))
      enddo
   enddo
!   do s1=1,noofpair
! this will later be replaced by cpair!! for entropy the old pair works better
!      phvar%mqmqaf%pair(s1)=pair(s1)
!      do s2=1,ncon
!         phvar%mqmqaf%dpair(s1,s2)=dcpair(s1,s2)
! try using dpair ....
!         phvar%mqmqaf%dpair=dpair(s1,s2)
!      enddo
!   enddo
!   write(*,777)'3X mqf sub1 1 copied:',(phvar%mqmqaf%yy1(s1),s1=1,nspin(1))
!   write(*,777)'3X mqf sub1 2 copied:',(phvar%mqmqaf%yy2(s1),s1=1,nspin(2))
!   write(*,777)'3X mqf pair copied:',(phvar%mqmqaf%pair(s1),s1=1,noofpair)
!   do s1=1,noofpair
!   write(*,777)'3X mqf dpair:',phvar%mqmqaf%pair(s1),&
!        (phvar%mqmqaf%dpair(s1,s2),s2=1,ncon)
!   enddo
777 format(a,F10.7,2x,5(F10.6),(/5x,6F10.6))
!---------------------------------------------------------------------------
! ENTROPY CALCULATION
!---------------------------------------------------------------------------
! separate documentation, i,j in first subl, k,l in second subl
! p_ijkl is cluster fraction; x_i site fraction; v_ik pair fraction
! w_i coordination equivalent site fraction;
! \sum_i y'_i ln(y'_i) + \sum_j y"_j ln(y"_j)+        subattice fractions
!
! \sum_i\sum_k v_ik ln(v_ik/(w_i w_k))+                 pair fractions
!
! \sum_i\sum_k p_iikk ln(p_iikk/((v^4_ik/(w^2_i w^2_k)))+           
! \sum_i\sum_j\sum_k p_ijkk ln(p_ijkk/(2(v^2_ik v^2_jk)/(w_i w_j w^2_k)))+  
! \sum_i\sum_k\sum_l p_iikl ln(p_iikl/(2(v^2_ik v^2_il)/(w^2_i w_k w_l)))+
! \sum_i\sum_j\sum_k\sum_l p_ijkl ln(
!                         p_ijkl/(4(v_ik v_il v_jk v_jl)/(w_i w_j w_k w_l)))
!---------------------------------------------------------------------------
! Discovered 21/10/20 with help by Mac Poschmann:
! The entropy is distributed on the quads, dS/dquad is the sum of
! the entropy contribution from sublattices, pairs and the quads
! is related to each separate quad!  Use the dyy1(*,quadindex) etc
!-----------------------------------------------------------------
! Here we calculte for one formula unit (FU) of the phase
! at the end we multiply with current number of atomes/FU
!-----------------------------------------------------------------
!
   ssub=zero; dssub=zero
   dvvv=zero
!   write(*,'(a,6(1pe12.4))')'3X quads: ',(phvar%yfr(q1),q1=1,ncon)
! NEW CODE, loop over all quads
   qsub: do q1=1,ncon
! Entropy from sublattices
      tsub=zero
! replace dsub with dvvv
      s7=0
      quady: do s1=1,4
! Entropy contribution from sublattice constituents for the quad
         s7=s7+1
         s2=mqmqa_data%contyp(10+s1,q1)
         fqq=one
         if(s2.gt.0) then
! Specie in first sublattice >0, if a single species fqq=2
            if(mqmqa_data%contyp(1,q1).eq.2) fqq=2.0d0
            tsub=tsub+fqq*log(yy1(s2))/mqmqa_data%constoi(s7,q1)
!            write(*,700)'3X ssub1: ',q1,s1,s2,s7,tsub,&
!                 fqq*log(yy1(s2))/mqmqa_data%constoi(s7,q1),fqq,yy1(s2),&
!                 mqmqa_data%constoi(s7,q1)
700         format(a,4i3,2(1pe12.4),4(0PF10.6))
! the derivative of fqq*log(yy1(s2))/mqmqa_data%constoi wrt all quads!
            do s3=1,ncon
               dvvv(s3,q1)=dvvv(s3,q1)+&
                    fqq*dyy1(s2,s3)/(yy1(s2)*mqmqa_data%constoi(s7,q1))
!               write(*,706)'3X dvvv1: ',q1,s2,s3,dvvv(s3,q1)
706            format(a,3i3,4(1pe12.4))
            enddo
         elseif(s2.lt.0) then
! if a single species in second sublattice fqq=2
            if(mqmqa_data%contyp(s1,q1).eq.2) fqq=2.0d0
            tsub=tsub+log(yy2(-s2))/mqmqa_data%constoi(s7,q1)
! the derivative of fqq*log(yy2(s2))/mqmqa_data%constoi wrt all quads!
            do s3=1,ncon
               dvvv(s3,q1)=dvvv(s3,q1)+&
                    fqq*dyy2(-s2,s3)/(yy2(-s2)*mqmqa_data%constoi(s7,q1))
!               write(*,706)'3X dvvv2: ',q1,s2,s3,dvvv(s3,q1)
            enddo
         else
! no more sublattice constituents
            exit quady
         endif
! exit if this is a pair
         if(mqmqa_data%contyp(5,q1).gt.0) exit quady
      enddo quady
! first derivatives, dSsub/dquad
      lsub(q1)=tsub
      ssub=ssub+phvar%yfr(q1)*tsub
!      write(*,702)'3X ssub2: ',q1,ssub,phvar%yfr(q1),tsub
702   format(a,i3,5(1pe12.4))
   enddo qsub
! correct first derivatives with respect to quads using dvvv
!   do q1=1,ncon
!      write(*,701)'3X dvvv: ',(dvvv(s1,q1),s1=1,ncon)
!   enddo
   do q1=1,ncon
      dssub(q1)=lsub(q1)
! add on all derivatives wrt q1 from other entropy terms 
      do s1=1,ncon
         dssub(q1)=dssub(q1)+phvar%yfr(s1)*dvvv(s1,q1)
      enddo
   enddo
! OK here
!   write(*,701)'3X dssub: ',(dssub(q1),q1=1,ncon)
!   write(*,701)'3X SSUB:',ssub,ssub*phvar%amfu,phvar%amfu,phvar%abnorm(1),&
!        phvar%amfu*phvar%abnorm(1)
701 format(a,6(1pe12.4))
!   stop
600 format(a,1pe12.4,2x,6(1pe10.2))
!===============
! skip the pair and quad contributions
!   write(*,*)'3X Done sublattice entropy, skipping rest',squad
!   goto 900
!
!-------------------------------------------------------
! pair entropy
   send=zero; dsend=zero
   quadcef: do q1=1,ncon
      tend=zero
! loop of all pairs of this quad
      s1=5
!      allpairs: do while(.TRUE. .and. s1.lt.10)
      allpairs: do while(.TRUE. .and. s1.lt.9)
! mqmqa_data%contyp(5,q1) is nonzero if the quad is a pair
         s2=mqmqa_data%contyp(s1,q1)
         if(s1.eq.5 .and. s2.ne.0) then
! the quad q1 is a pair with index s2, only one calculation with s2=q1
            fqq=4.0D0
            s1=10
         else
            s1=s1+1
            s2=mqmqa_data%contyp(s1,q1)
! s2 is now the index a pair in this SNN quad is in %contyp(6..9,q1) 
! exit here ifthere is no pair
            if(s2.eq.0) exit allpairs
! fqq depends on q1
            fqq=1.0D0
            if(mqmqa_data%contyp(1,q1).eq.2) then
               fqq=2.0D0
            elseif(mqmqa_data%contyp(3,q1).eq.-2) then
               fqq=2.0D0
            endif
         endif
! Here s2 is a pair of the quadrupole q1.  The pair fraction is pair(s2)
! which should be divided by ceqf1(1,s2)*ceqf2(2,s2)
! The logarithm should be multiplied by qfnnsnn for the pair.  no more??
! Entropy: quadfrac*\sum_s2 fqq*ln( pair(s2)/v_s2k/(w_i w_k))/%qfnnsnn(s2)
! MAYBE save values of "pair/(ceqf1*ceqf2)" and derivaties for later use??
! REMEMBER ceqf1 is equivalent sublattice fraction ... what is eij(1,s2)??
! eij(1..2,s2) are species in first and second sublattice of the pair
! BUT they are now in %contype(11,s2) and %contyp(12,s2) ???
! KEEP eij as it is used as link from pair to sublattice constituents
!         write(*,'(a,i3,2x,2i3,2x,2i3)')'3X keep eij?: ',s2,eij(1,s2),&
!              eij(2,s2),mqmqa_data%contyp(11,s2),mqmqa_data%contyp(12,s2)
         ee=eij(1,s2); gg=-eij(2,s2)
         dq1=ceqf1(ee)*ceqf2(gg)
         mulceq(s2)=dq1
         endkvot(s2)=pair(s2)/dq1
         fqq=fqq/mqmqa_data%qfnnsnn(s2)
! >>>>>>>>>>>>>>>>  ............. EQUATION B21 2nd line first half
! This is the entropy contribution from a pair of this quad
! %qfnnsnn is read from database
! %dfnnsnn can be different for different pairs, composition dependence???
! But it should be a sum? or is that taken care of by the sum over p_AB/XY ??
         tend=tend+fqq*log(endkvot(s2))
!         write(*,421)'3X pairs: ',q1,s1,s2,tend,endkvot(s2),&
!              fqq/mqmqa_data%qfnnsnn(s2),fqq,mqmqa_data%qfnnsnn(s2)
421      format(a,3i3,5(1pe11.3))
! first derivatives, note multiplied by p_AB/XY ....
         do s3=1,ncon
            if(s3.eq.q1) dsend(s3)=dsend(s3)+fqq*log(endkvot(s2))
            dsend(s3)=dsend(s3)+fqq/endkvot(s2)*(&
                 dpair(s2,q1)/(mulceq(s2))**2-&
                 2.0d0*pair(s2)/mulceq(s2)**4*(&
                 ceqf1(ee)*dceqf2(gg,q1)+dceqf1(ee,q1)*ceqf2(gg)))
! skip 2nd derivatives ...
         enddo
      enddo allpairs
! Finally we must multiply the tend with the quad fraction
      send=send+phvar%yfr(q1)*tend
! derivatives of send wrt quad
   enddo quadcef
!
! ternary error before this
!   write(*,600)'3X SEND: ',send,(dsend(s1),s1=1,ncon)
!========================================================================
! skip quad entropies
!   write(*,*)'3X done pair entropies'
!   write(*,*)'3X skipping quad entropies'
!   goto 900
!========================== begin loop for all quads
!   write(*,*)'3X quadropole entropies:'
!   do s1=1,noofpair
!      write(*,440)'3X dpair/dq: ',q1,(dpair(s1,s2),s2=1,ncon)
!   enddo
440 format(a,i2,6(1pe10.2),(/20x,6e10.2))
   squad=zero; dsquad=zero
! replaced s1 by q1
   quadloop: do q1=1,ncon
      if(q1.ne.mqmqa_data%contyp(10,q1)) then
! TEST: the value in contyp(10,q1) should be q1 ...
         write(*,*)'3X problems in %contyp with quad indexing 7'
         gx%bmperr=4399; goto 1000
      endif
      lsub=zero
! New code for the general case
!                  p_i
! p_i * log( ------------------------------------)
!               xi_A/X*xi_B/X*xi_B/X*xi_B/Y
!               ---------------------------
!                  w_A * w_B * w_X * w_Y
!
      s1=mqmqa_data%contyp(5,q1)
      if(s1.gt.0) then
! this is a pair
         pair1=s1
         pair2=pair1
         pair3=pair1
         pair4=pair1
         ee=eij(1,pair1)
         ff=ee
         gg=-eij(2,pair1)
         hh=gg
! before adding this write statement hh was sometines not same as gg
! as it should be SUCK
!         write(*,*)'3X gg hh: ',gg,hh,ceqf2(gg),ceqf2(hh)
         fqq=one
!         write(*,'(a,10i3)')'3X quad1: ',q1,pair1,pair2,pair3,pair4,ee,ff,gg,hh
      elseif(mqmqa_data%contyp(9,q1).eq.0) then
! here either ee=ff or gg=hh
         pair1=mqmqa_data%contyp(6,q1)
         pair2=pair1
         ee=eij(1,pair1)
         gg=-eij(2,pair1)
         pair3=mqmqa_data%contyp(7,q1)
         pair4=pair3
         ff=eij(1,pair3)
         hh=-eij(2,pair3)
         fqq=2.0d0
!         write(*,'(a,10i3)')'3X quad2: ',q1,pair1,pair2,pair3,pair4,ee,ff,gg,hh
      else
! all ee, ff, gg, hh should be different, not certain if they are
         pair1=mqmqa_data%contyp(6,q1)
         ee=eij(1,pair1)
         gg=-eij(2,pair1)
         pair2=mqmqa_data%contyp(7,q1)
         ff=eij(1,pair2)
         hh=-eij(2,pair2)
         pair3=mqmqa_data%contyp(8,q1)
         if(ee.eq.ff) ff=eij(1,pair3)
         if(gg.eq.hh) hh=-eij(2,pair3)
         pair4=mqmqa_data%contyp(9,q1)
         fqq=4.0D0
!         write(*,'(a,10i3)')'3X quad4: ',q1,pair1,pair2,pair3,pair4,ee,ff,gg,hh
      endif
!
!      write(*,'(a,8F8.4)')'3X quadx: ',pair(pair1),ceqf1(ee),&
!           pair(pair2),ceqf1(ff),pair(pair3),ceqf2(gg),pair(pair4),ceqf2(hh)
      pairceq=fqq*pair(pair1)/ceqf1(ee)*pair(pair2)/ceqf1(ff)*&
           pair(pair3)/ceqf2(gg)*pair(pair4)/ceqf2(hh)
!      write(*,'(a,9i3,1pe12.4)')'3X quadx: ',q1,pair1,pair2,pair3,pair4,&
!           ee,ff,gg,hh,pairceq
!
      squad=squad+phvar%yfr(q1)*log(phvar%yfr(q1)/pairceq)
!      write(*,440)'3X squad: ',q1,squad,phvar%yfr(q1),pairceq
!
! New code for the general case
!                  p_i
! p_i * log( ------------------------------------)
!               xi_A/X*xi_B/X*xi_B/X*xi_B/Y
!               ---------------------------
!                  w_A * w_B * w_X * w_Y
!
! loop for derivatives
      do s1=1,ncon
         if(s1.eq.q1) lsub(s1)=log(phvar%yfr(q1)/pairceq)+one
         if(s1.eq.q1) dsquad(s1)=dsquad(s1)+log(phvar%yfr(q1)/pairceq)+one
! derivative for just q1 is OK
         lsub(s1)=lsub(s1)-phvar%yfr(q1)*&
              (dpair(pair1,s1)/pair(pair1)+dpair(pair2,s1)/pair(pair2)+&
              dpair(pair3,s1)/pair(pair3)+dpair(pair4,s1)/pair(pair4)-&
              dceqf1(ee,s1)/ceqf1(ee)-dceqf1(ff,s1)/ceqf1(ff)-&
              dceqf2(gg,s1)/ceqf2(gg)-dceqf2(hh,s1)/ceqf2(hh))
! Skipping this means I ignore effect of variable fracrion on pair and ceqf
!         dsquad(s1)=dsquad(s1)-phvar%yfr(q1)*&
!              (dpair(pair1,s1)/pair(pair1)+dpair(pair2,s1)/pair(pair2)+&
!              dpair(pair3,s1)/pair(pair3)+dpair(pair4,s1)/pair(pair4)-&
!              dceqf1(ee,s1)/ceqf1(ee)-dceqf1(ff,s1)/ceqf1(ff)-&
!              dceqf2(gg,s1)/ceqf2(gg)-dceqf2(hh,s1)/ceqf2(hh))
! skip 2nd derivatives
!         write(*,440)'3X lsub: ',s1,(lsub(s2),s2=1,ncon)
      enddo
!      write(*,440)'3X SQUAD: ',q1,squad,(dsquad(s1),s1=1,ncon)
   enddo quadloop
!
!   write(*,600)'3X SQUAD: ',squad,(dsquad(s1),s1=1,ncon)
! first derivatives are wrong ....
!   dsquad=zero
!   write(*,*)'3X done quad derivatives'
!   goto 900
!
!***********************************************************************
900 continue
! we have multiplied with amounts above, (?) set invnorm=one
!   write(*,*)'3X second derivatives are approximate.  Atoms/FU: ',invnorm
! Values should be per formula unit!
   invnorm=one
! store results in appropriate places, values divided by RT
! This is G/RT
   phvar%gval(1,1)=phvar%gval(1,1)+invnorm*(ssub+send+squad)
! derivative of G wrt T, i.e. -S/R
   phvar%gval(2,1)=phvar%gval(2,1)+invnorm*(ssub+send+squad)/tval
   if(moded.gt.0) then
! This is if first derivatives are requested (must be exact)
!      write(*,*)'3X start quad loop'
      do s1=1,ncon
         phvar%dgval(1,s1,1)=phvar%dgval(1,s1,1)+&
              invnorm*(dssub(s1)+dsend(s1)+dsquad(s1))
         phvar%dgval(2,s1,1)=phvar%dgval(2,s1,1)+&
              invnorm*(dssub(s1)+dsend(s1)+dsquad(s1))/tval
         if(moded.gt.1) then
! this is if second derivatives are requested
!            do s2=s1,ncon
!               phvar%d2gval(ixsym(s1,s2),1)=phvar%d2gval(ixsym(s1,s2),1)+&
!                    invnorm*d2sm1(ixsym(s1,s2))
!            enddo
! We just set 1/quad
            dummy1=phvar%yfr(s1)
            if(dummy1.lt.1.0D-12) dummy1=1.0D-12
            phvar%d2gval(ixsym(s1,s1),1)=one/dummy1
         endif
      enddo
!      write(*,*)'3X done quad loop'
!      write(*,431)'3X dS/Rq  :',(phvar%dgval(1,s1,1),s1=1,ncon)
!      write(*,431)'3X d2S/Rq2:',(phvar%d2gval(s1,1),s1=1,all2)
431   format(a,6(1pe12.4),(/6x,6e12.4))
   endif
!   mqf=>phvar%mqmqaf ??
!   write(*,*)'3X pair do loop npair: ',phvar%mqmqaf%npair
!   write(*,*)'3X pair do loop mqf%pair: ',allocated(phvar%mqmqaf%pair)
!   write(*,'(a,3(1pe14.6))')'3X MQMQA:',phvar%gval(1,1),&
!        phvar%gval(1,1)*8.31451,phvar%gval(1,1)*8.31451*phvar%amfu
! replace pair by cpair to handle endmembers
! Creates problems calculating the entropy in this routine ... SUCK
!   write(*,*)'3X pair do loop mqf%dpair: ',allocated(phvar%mqmqaf%dpair)
   do s1=1,phvar%mqmqaf%npair
      phvar%mqmqaf%pair(s1)=cpair(s1)
      do s2=1,ncon
         phvar%mqmqaf%dpair(s1,s2)=dcpair(s1,s2)
! converge problems, maybe use dp?
!         mqf%dpair(s1,s2)=dp(s1,s2)
      enddo
   enddo
   if(ddebug) write(*,*)'3X Done MQMQA configurational entropy'
! TEST temporary fix
!   do s1=1,mqf%npair
!      write(*,'(a,F9.6,2x,10F10.6)')'3X cpair: ',mqf%pair(s1),&
!           (mqf%dpair(s1,s2),s2=1,mqf%nquad)
!   enddo
!
1000 continue
   return
 end subroutine config_entropy_mqmqa1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
 
!\addtotable subroutine calc_mqmqa
!\begin{verbatim}
 subroutine calc_mqmqa(lokph,phres,ceq)
! Called from calcg_internal to calculate nonconfig G for the mqmqa phase
! another subroutine calculates the entropy using all data in phres%mqf
   implicit none
   integer lokph
   type(gtp_phase_varres), pointer :: phres
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! Most variables here are the same as in calcg_internal ...
   integer, parameter :: f1=50
   integer mqmqj,kend,s1,s2,s3,id,nofc2,ipy,lokfun,typty,itp,zp,nrealem,mqendx
   double precision vals(6),pyq,rtg,aff
   double precision, dimension(:), allocatable :: dpyq(:),d2pyq(:),d2vals(:)
   double precision, dimension(:,:), allocatable :: dvals(:,:),affarr(:)
! for saving FNN reference energies
   double precision refg(f1,f1)
   double precision dummy1,dummy2
! for MQMQA minimal fractions
   double precision, parameter :: MINMQMQA=1.0D-5
   TYPE(gtp_parcalc) :: gz
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_endmember), pointer :: endmemrec
   TYPE(gtp_interaction), pointer :: intrec
   TYPE(gtp_pystack), pointer :: pystack
   TYPE(gtp_phase_add), pointer :: addrec
   TYPE(gtp_mqmqa_var), pointer :: mqf
   TYPE(gtp_tooprec), pointer :: tooprec
! for handling excess parameters, just binary, use no mqmqa_data ksi arrays
   integer ij,jd,jq,qq1,qq2,ass,mpow,isumx,tsize,tch
   double precision ksi,sumx,dsumx
   double precision dksi(3),d2ksi(3)
!   logical ddebug
!------------------------------------- 
! tch is level of debug output, 0=none, 3=max
   tch=0
!   ddebug=.FALSE.
!   ddebug=.TRUE.
   if(tch.ge.1) write(*,*)'3XQ in calc_mqmqa nonconfig G'
   gz%nofc=phlista(lokph)%tnooffr
   nofc2=gz%nofc*(gz%nofc+1)/2
!   write(*,*)'3X allocating:',gz%nofc,nofc2
   allocate(dpyq(gz%nofc))
   allocate(d2pyq(nofc2))
   allocate(dvals(3,gz%nofc))
   allocate(d2vals(nofc2))
! this shortcut may be bad - but it works ---------------------------------
!   write(*,*)'3XQ assigning mqf pointer'
   mqf=>phres%mqmqaf
!   write(*,*)'3XQ assigning mqf pointer OK'
!-------------------
   allocate(affarr(mqf%npair))
   affarr=zero
   nullify(pystack)
   rtg=globaldata%rgas*ceq%tpval(1)
!   do s1=1,mqmqa_data%nconst
!      write(*,599)s1,(mqmqa_data%contyp(s2,s1),s2=1,14)
!599   format('3XQ contyp 7: ',i2,1x,4i2,1x,i3,1x,4i2,1x,i2,4i3)
!   enddo
   nrealem=0
!   refg=zero
   dummy2=zero
! list %pp
! %pp( quad , FNN index )
!   do mqmqj=1,mqmqa_data%nconst
!      write(*,17)'3XQ %pp: ',mqmqj,(mqmqa_data%pp(s1,mqmqj),s1=1,4)
!   enddo
17 format(a,i3,4(1pe12.4))
!--------------------------------------
! first loop over ALL endmembers
   mqmqj=0
   endmemrec=>phlista(lokph)%ordered
! This should be number of atoms for scaling G
!   dummy1=phres%abnorm(1)/rtg       this was OK before ...
   dummy1=one/rtg
! %amfu * %abnorm(1) is number of moles in the liquid
! in the test case we have 6 atoms in the liquid phase
!   dummy1=6.0D0/rtg
!   dummy1=one/(phres%abnorm(1)*rtg)
!   write(*,'(a,3(1pe14.6))')'3XQ mqmqa scaling: ',dummy1,&
!        phres%amfu,phres%abnorm(1)
! This first loop: all endmember parameters
! this can give SRO contribution and excess from SNN parameters
! or it makes it possible to calculate the G for the FNN parameters
   endmemloop1: do while(associated(endmemrec))
      mqmqj=mqmqj+1
      if(mqmqj.gt.mqmqa_data%nconst) exit endmemloop1
! We do not know if mqmqj is associated with this endmember!!
! there can be gaps in the endmember list?? 
! we must take kend from the endmember record, it is sored in %antalem
      mqendx=endmemrec%antalem
      kend=mqmqa_data%contyp(5,mqendx)
!      write(*,*)'3XQ endmemloop1A: ',mqmqj,mqendx,kend,nrealem
      if(kend.le.0) then
! This is an SNN parameter we calculate and add SNN energy and interactions ...
!         write(*,*)'3XQ SNN endmember record found',mqmqj
         proprec=>endmemrec%propointer
         mqsnn: do while(associated(proprec))
! This loop is not really necessay, in mqmqa the only property is G at present
            typty=proprec%proptype
            if(typty.ne.1) stop '3XQ illegal typty in mqmqa model'
            ipy=1
            lokfun=proprec%degreelink(0)
            call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
            if(gx%bmperr.ne.0) goto 1000
!            write(*,'(a,2i3,2(1pe12.4))')'3XQ SNN endmember',mqmqj,kend,&
!                 pyq,vals(1)
!            write(*,'(a,6(1Pe12.4))')'3XQ vals1:',vals
            if(ipy.eq.1) then
               vals=vals*dummy1
! This is an SNN ordering parameter, reference state addel in second loop
            endif
            pyq=phres%yfr(mqmqj)
! Should I use any factor??
!         aff=mqmqa_data%pp(1,mqmqj)
            aff=one
! NOTE the reference state contribution to this SNN added in next loop
! for all quads!!
            do itp=1,3
               phres%dgval(itp,mqmqj,ipy)=phres%dgval(itp,mqmqj,ipy)+vals(itp)
            enddo
! Initially ignore 2nd derivatives, d2G/dy2=1/y set by entropy calculation
! ipy is property, ipy=1 means G, ipy=2 means Curie T etc.
! %gval(1,1) is total G, %gval(2,1) is total dG/dT  etc.
            do itp=1,6
               phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*vals(itp)
            enddo
!            write(*,210)'3XQ SRO G, dG/dqi: ',mqmqj,mqmqj,pyq,aff,&
!                 phres%gval(1,1),(phres%dgval(1,s1,1),s1=1,gz%nofc)
            proprec=>proprec%nextpr
         enddo mqsnn
!600      continue
!         write(*,*)'3XQ any excess parameters will be handled in 3rd loop'
         endmemrec=>endmemrec%nextem
         cycle endmemloop1
      endif
! This is an FNN parameter, we calculate and save the value for later use
      nrealem=nrealem+1
!      write(*,*)'3XQ endmemloop1B: ',mqmqj,kend,nrealem
      proprec=>endmemrec%propointer
      aff=one/mqmqa_data%pp(1,mqmqj)
      mq1: do while(associated(proprec))
         typty=proprec%proptype
         if(typty.ne.1) stop 'illegal typty in mqmqa model'
         ipy=1
         lokfun=proprec%degreelink(0)
         call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
         if(gx%bmperr.ne.0) goto 1000
!         write(*,'(a,i3,F7.4,3(1Pe10.2))')'3XQ refg:',mqmqj,aff,vals(1),vals(2)
! we should divide this by the aff of this pair and we will multiply this
! FNN same aff but SNN fractions linking to this pair use another aff
!         write(*,'(a,2i3,2(1pe12.4))')'3XQ FNN endmember',mqmqj,kend,&
!                 pyq,vals(1)
         if(ipy.eq.1) then
            vals=vals*dummy1*aff
! save values of reference state for use with SNN parameters ??
! kend is FNN (pair) index 
            do s1=1,6
               refg(kend,s1)=vals(s1)
            enddo
         endif
! next property record (should not be any ...)
         proprec=>proprec%nextpr
         if(associated(proprec)) then
            write(*,*)'3XQ Warning: ignoring second mqmqa property recotd!'
         endif
!         write(*,200)'3XQ FNN G, dG/dqi: ',phres%gval(1,1),&
!              (phres%dgval(1,s1,1),s1=1,gz%nofc)
200      format(a,1pe12.4,2x,6(1pe12.4))
      enddo mq1
      endmemrec=>endmemrec%nextem
   enddo endmemloop1
!   write(*,*)'3XQ finished endmemloop1'
!--------------------------------------------------- end first endmember loop
! All endmembers with a single element in each sublattice must have a parameter
! these are counted above in endmemloop1
!   write(*,'(a,3i3)')'3XQ number of sublattice constituents and FNN: ',&
!        mqf%ns1,mqf%ns2,nrealem
   if(nrealem.ne.mqf%ns1*mqf%ns2) then
! This test is not foolproof one can enter an interaction parameter
! which creates an empty endmember record but that seems crazy
      write(*,216)mqf%ns1*mqf%ns2,nrealem
216   format('Some FNN constituents (A/X) have no parameter!, should be',&
           i3,' found only ',i3)
      gx%bmperr=4399; goto 1000
   endif
! second loop over all constutents (quads), ignore FNN endmember records
! but add reference state parameters to all SNN and reciprocal constituents
   ipy=1
   if(tch.ge.3) write(*,*)'3XQ adding reference to SNN endmembers'
   qloop: do mqmqj=1,gz%nofc
! this is quad fraction, multiply with all FNN reference energies
      pyq=phres%yfr(mqmqj)
      zp=mqmqa_data%contyp(5,mqmqj)
      pair: if(zp.gt.0) then
! this is an FNN  pair, reference energy in refg(zp,1..6), only one y derivative
! %pp(1..4,mqmqj) is stoichiometric factor for the pair
         aff=mqmqa_data%pp(1,mqmqj)
         do itp=1,3
            phres%dgval(itp,mqmqj,ipy)=phres%dgval(itp,mqmqj,ipy)+&
                 aff*refg(zp,itp)
         enddo
! Initially ignore 2nd derivatives, d2G/dy2=1/y set by entropy calculation
         do itp=1,6
            phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*aff*refg(zp,itp)
         enddo
!         write(*,205)'3XQ FNN: qix, FNN, aff, pyq, fun, DG: ',mqmqj,zp,aff,&
!              pyq,refg(zp,1),pyq*aff*refg(zp,1)
205      format(a,2i3,F8.5,2x,3(1pe12.4))
         if(tch.ge.3) &
          write(*,210)'3XQ FNN G:     ',mqmqj,mqmqj,pyq,aff,pyq,phres%gval(1,1)
210      format(a,2i3,2F8.5,1pe12.4,2x,6(1pe10.2))
      else
! this is an SNN with two or more pairs
! For each SNN pair add the contribution to the FNN reference state
! %contyp(1..4,mqmqj) is index of FNN reference energy
         if(tch.ge.3) write(*,'(a,i3,1x,4i3,4F8.5)')'3XQ pp2: ',mqmqj,&
              (mqmqa_data%contyp(s1,mqmqj),s1=6,9),&
              (mqmqa_data%pp(s1,mqmqj),s1=1,4)
         snnloop: do s1=6,9
! zp is index to an FNN record, there can be 2 or 4 FNN records
            zp=mqmqa_data%contyp(s1,mqmqj)
            if(zp.eq.0) exit snnloop
! %pp(1..4,mqmqj) is stoichiometric factor for the pair
            aff=mqmqa_data%pp(s1-5,mqmqj)
!            write(*,211)1,mqmqj,ipy,phres%dgval(1,mqmqj,ipy)
211         format('3XQ SNN dG/dy:',3i3,1(1pe12.4))
            do itp=1,3
               phres%dgval(itp,mqmqj,ipy)=phres%dgval(itp,mqmqj,ipy)+&
                    aff*refg(zp,itp)
            enddo
!           write(*,212)zp,mqmqj,ipy,phres%dgval(1,mqmqj,ipy),aff,aff*refg(zp,1)
212         format('3XQ SNN dG/dy reference added:',3i3,3(1pe12.4))
! Initially ignore 2nd derivatives, d2G/dy2=1/y set by entropy calculation
            if(tch.ge.3) write(*,213)s1,zp,mqmqj,ipy,pyq,aff,phres%gval(1,ipy)
213         format('3XQ SNN G ref:',4i3,2F10.5,1pe12.4)
            do itp=1,6
               phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*aff*refg(zp,itp)
            enddo
!            write(*,214)zp,mqmqj,ipy,phres%gval(1,ipy),pyq,aff,aff*refg(zp,1)
!            write(*,214)zp,mqmqj,ipy,phres%gval(1,ipy),pyq,aff,refg(zp,1)*rtg
214         format('3XQ SNN G ref added:',3i3,4(1pe12.4))
!            write(*,205)'3XQ SNN: qix, FNN, aff, pyq, fun, DG: ',mqmqj,zp,aff,&
!                 pyq,refg(zp,1),pyq*aff*refg(zp,1)
!                 (phres%dgval(1,s2,1),s2=1,gz%nofc)
         enddo snnloop
      endif pair
   enddo qloop
   if(tch.ge.3) write(*,*)'3XQ finished loop for endmembers'
! if this goto then excess is ignored and result correct
!   goto 800
!---------------------------------------------------------------------
! code below needed for excess parameters ONLY, all SNN FNN endmembers done
! NOTE some of them may not have a parameter
! This is to allocate csumx for handling quads with small fractions.
   isumx=0
!
   mqmqj=0
   endmemrec=>phlista(lokph)%ordered
   endmemloop2: do while(associated(endmemrec))
      if(mqmqj.gt.0) endmemrec=>endmemrec%nextem
      mqmqj=mqmqj+1
      if(tch.ge.3) write(*,*)'3XQ endmemloop2:',&
           mqmqj,mqmqa_data%nconst,associated(endmemrec)
      if(mqmqj.gt.mqmqa_data%nconst .or. .not.associated(endmemrec)) then
         exit endmemloop2
      endif
      kend=mqmqa_data%contyp(5,mqmqj)
      if(tch.ge.3) write(*,311)associated(endmemrec),&
           mqmqj,mqmqa_data%nconst,kend
311   format('3XQ second loop for endmember records: ',l2,3i5)
      intrec=>endmemrec%intpointer
! interaction parameters are NOT linked from SNN endmembers ?? really?
! They are stored in alphabetical order of the constituents
!      write(*,*)'3XQ Check interaction parameters 1',associated(endmemrec),&
!           associated(intrec),mqmqj,kend
! if we cycle here the results are the same as without excess parameters
!      cycle endmemloop2
!
      if(.not.associated(intrec)) cycle endmemloop2
! There are excess parameters, any Tooprecords?
      if(associated(intrec%tooprec)) then
! the allocatable arrays Toop1, Toop2 and Kohler have all same size
! equal to the number of binary combination of constituents
         tooprec=>intrec%tooprec
         if(tch.ge.3) then
            write(*,'(a,2i3,l2)')'3XQ A Toop/Kohler record, id:',&
                 tooprec%toopid,tooprec%endmemel,associated(tooprec%binint)
            if(allocated(tooprec%toop1)) then
               tsize=size(tooprec%toop1)
               write(*,320)'Toop1 ',(tooprec%toop1(jd),jd=1,tsize)
               write(*,320)'Toop2 ',(tooprec%toop1(jd),jd=1,tsize)
               write(*,320)'Kohler ',(tooprec%kohler(jd),jd=1,tsize)
320            format('3XQ ',a,': ',10i3)
            endif
         endif
! this is an excess parameter with possible excess parameters
!      write(*,'(a,2i3)')'3XQ endmember with excess parameter:',mqmqj
! just excess parameters, we must calculate product of fractions
! BRANCH for intrec%highlink and intrec%nexlink
!      write(*,'(a,i2,F10.6,6(1pe12.4))')'3XQ SNN df/dy: ',id,pyq,&
!           (dpyq(itp),itp=1,gz%nofc)
!-------------------------------------
! content of %contyp and %pinq
!      do jd=1,mqmqa_data%nconst
!         write(*,599)jd,(mqmqa_data%contyp(id,jd),id=1,14)
!599      format('3XQ contyp: ',i2,1x,4i2,1x,i3,1x,4i2,1x,i2,4i3)
!      enddo
!      write(*,*)'3XQ pinq: ',mqmqa_data%pinq
! extract fractions from the endmember and check if AB/X or A/XY or A/X 
      end if
!-------------------------------------- code below ignore Toop/Kohler
      id=endmemrec%fraclinks(1,1)
! jump back here for next interaction record (if any)
600   continue
! Note it is arbitrary if the cluster is endmember or interaction
      jd=intrec%fraclink(1)
! We must keep track of which endmember is separate!!!
      if(mqmqa_data%contyp(5,id).eq.0) then
! id is a cluster, jd is separate fraction, jq is additional salt OK
         ass=id
! %contyp(6,..9) are index of FNN, pairs, FNN pairs index in cintyp in PINQ
         jq=mqmqa_data%pinq(mqmqa_data%contyp(6,ass))
         if(jq.eq.jd) jq=mqmqa_data%pinq(mqmqa_data%contyp(7,ass))
         qq1=jd
         qq2=jq
!         write(*,'(a,6i3)')'3XQ ass, sep, sum 1:',ass,qq1,qq2
      elseif(mqmqa_data%contyp(5,jd).eq.0) then
! jd is the cluster, id is interaction endmember WRONG
         ass=jd
         jq=mqmqa_data%pinq(mqmqa_data%contyp(6,ass))
         if(jq.eq.id) jq=mqmqa_data%pinq(mqmqa_data%contyp(7,ass))
         qq1=id
         qq2=jq
!         write(*,'(a,6i3)')'3XQ ass, sep, sum 2:',ass,qq1,qq2
      else
! Interactions are only between clusters AB/X and endmembers A/X or B/X
         write(*,*)'3XQ interaction between two endmembers illegal'
         gx%bmperr=4399; goto 1000
      endif
!      write(*,428)phres%yfr
428   format('3XQ all yfr: ',20(1x,F8.6))
      if(tch.ge.3) write(*,430)id,jd,jq,qq1,qq2,ass,&
           phres%yfr(id),phres%yfr(jd),phres%yfr(jq)
430   format('3XQ interaction: ',3i3,3x,3i3,3x,3(1x,F8.6))
!------------------------------------- extract parameter value
      proprec=>intrec%propointer
      typty=proprec%proptype
      if(typty.ne.1) stop 'illegal typty in mqmqa model'
      ipy=1
! several powers  we must loop here -------------- not yet done
      if(proprec%degree.gt.0) write(*,*)'3XQ degree: ',proprec%degree
      mpow=0
700   continue
! first power is in link 0
      lokfun=proprec%degreelink(mpow)
      mpow=mpow+1
      if(mpow.gt.9) then
         write(*,*)'3XQ too high interaction power'
         gx%bmperr=4399; goto 1000
      endif
! some powers may not have a parameter, max 9.  If no function loop
      if(lokfun.le.0) goto 700
      call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
      if(gx%bmperr.ne.0) goto 1000
      if(tch.ge.3) write(*,'(a,3i4/4x,6(1Pe12.4))')'3XQ excess1:',&
           lokfun,mqmqj,mpow,vals(1)
      if(ipy.eq.1) then
! Nath has implemented this half in the converter
!         vals=0.5D0*vals/rtg
         vals=vals/rtg
      endif
! skip excess 1
!      cycle endmemloop2
!----------------------- multiply with fractions
! the parameter should be multiplied with cluster fractions and
! the separate endmember qq1 fraction normalized 
      isumx=isumx+1
      sumx=phres%yfr(qq1)+phres%yfr(qq2)+phres%yfr(ass)
      ksi=phres%yfr(qq1)/sumx
      if(mpow.eq.1) then
         pyq=phres%yfr(ass)*ksi
! most of the derivatives of pyq is zero
         dpyq=zero
         dsumx=-sumx**(-2)
! only those involving id, jd and jq are nonzero.
! the species qq1, qq2 and ass has one more term, qq2 is only in the sumx
!         dpyq(qq1)=pyq*dsumx+phres%yfr(ass)/sumx
!         dpyq(qq2)=pyq*dsumx
!         dpyq(ass)=pyq*dsumx+ksi
! corrected derivatives ...
         dpyq(qq1)=(phres%yfr(ass)-pyq)/sumx
         dpyq(ass)=(phres%yfr(qq1)-pyq)/sumx
         dpyq(qq2)=-pyq/sumx
      else
! NOT CORRECTED THESE ... suck
         pyq=phres%yfr(ass)*(ksi**mpow)
         dpyq=zero
         dsumx=-mpow*sumx**(-mpow-1)
         dpyq(qq1)=pyq*dsumx+mpow*phres%yfr(ass)*ksi**(mpow-1)
         dpyq(qq2)=pyq*dsumx
         dpyq(ass)=pyq*dsumx+ksi*mpow
      endif
! here the fraction product is calculated
!      write(*,650)ass,qq1,qq2,mpow,ksi,phres%yfr(ass),pyq,sumx,vals(1)*rtg
650   format('3XQ excess: ',4i3,4(1x,F8.6),1pe12.4)
!      write(*,'(a,2(1pe14.6))')'3XQ excess G:',pyq,pyq*vals(1)
! skip excess 2
!      cycle endmemloop2
!
! ---------------------------------
! add to G and first derivatives of G, ipy is property, ipy=1 is G
! 2nd derivatives ignored
! ---------------------------------
      do s1=1,gz%nofc
         do itp=1,3
            phres%dgval(itp,s1,ipy)=phres%dgval(itp,s1,ipy)+&
                 dpyq(s1)*vals(itp)
         enddo
      enddo
      do itp=1,6
         phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*vals(itp)
      enddo
! maybe several fraction powers of this property
!      write(*,*)'3XQ several powers? ',mpow,proprec%degree
      if(mpow.lt.proprec%degree) goto 700
!------------- next property for same interaction,
! each property can have different number of powers ... not implemented
      proprec=>proprec%nextpr
      if(associated(proprec)) then
! more than one property ... not implemented
         write(*,*)'3XQ MQMQA parameter with several propertyes!',mqmqj
      endif
      if(associated(intrec%highlink)) then
! a higher interaction ... not allowed
         write(*,*)'3XQ ternary MQMQA parameters not implemented',mqmqj
      endif
! there can be more interactions on this level
      intrec=>intrec%nextlink
      if(associated(intrec)) then
! There can be more than one interaction linked from an endmember
         if(tch.ge.3) write(*,*)'3XQ more interaction for an endmember',mqmqj
         goto 600
      endif
!      write(*,*)'3XQ done excess for endmember',mqmqj
! next endmember .... is set at the beginning

   enddo endmemloop2
!----------------------------------------------------- end SNN loop
800 continue
!   write(*,990)'3XQ exit calc_mqmqa G:',phres%gval(1,1),&
!        (phres%dgval(1,s1,1),s1=1,gz%nofc)
990 format(a,5(1pe14.6))
1000 return
 end subroutine calc_mqmqa

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine calc_toop
! called from cgint(lokph,lokpty,moded,vals,dvals,d2vals,gz,ceq)
!\begin{verbatim}
 subroutine calc_toop(lokph,lokpty,moded,vals,dvals,d2vals,gz,TOOPX,ceq)
! NOT USED FOR MQMQA liquid model ... done in calc_mqmqa
! This routine replaces all calculations inside cgint for Toop/Kohler excess
! binary interaction parameter with Toop or Kohler extrapolation
! toopx is the pointer to the kohler-Toop record
! toopx%binint is pointer back to calling subroutine
! A single composition dependent binary parameter is calculated
! But in the Toop/Kohler we can have additional fraction variables
   implicit none
   integer moded,lokph
   TYPE(gtp_property), pointer :: lokpty
   TYPE(gtp_parcalc) :: gz
! all fraction variable can be involved in derivatives of vals ...
   double precision vals(6),dvals(3,gz%nofc)
   double precision d2vals(gz%nofc*(gz%nofc+1)/2)
   TYPE(gtp_tooprec), pointer :: toopx
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! we use this to save the pointer from toopx
   TYPE(gtp_phase_varres), pointer :: phres
! fraction values to be used in RK series
   double precision x12,x21,sigma,dxrk,dxrk0
   double precision, allocatable, dimension(:) :: dsigma, dx12, dx21
! ternary fraction index
   integer jj(3),j1,j2,j3,link,count,toopconst,limit,jdeg,lfun,nyfr,tkdeb
! loop veriables 
   integer qz,ic,cc
! to avoid calculating derivatives if no constituents in toop1, toop2 or kohler
   logical not1,not2,nok
! for the RK calculation with Toop/Kohler fractions!
   double precision valtp(6)
   double precision dx,dx0,dx1,dx2,dxi,dxj,fff,rtg
! The first part here is to modify the fractions to be used in the RK series
! the gz record has information which elements involved
! gz%iq(1) and gz%iq(2) are index of the binary constituents
! We must also handle first and second derivatives wrt all fractions.
!    as the binary fractions are modified by adding or subtractions
! we come here from a binary interaction record will only deal with this
!
! These are UNUSED arrays with additional fractions to calculate derivatives
   integer dtoop1(5),dtoop2(5),dkohler(10),ntp1,ntp2,nkh
! These are arrays to eliminate cases with duplicate fractions in Toop1/2/Kohler
   integer, allocatable, dimension(:) :: ctoop1,ctoop2,ckohler
   integer nz
!
! Use the phres passed on via toopx%phres if there are more toopx records
! this link to phres is copied to toopx%phres before the call.
! In gtp3X, subroutine calcg_internal around line 858.
! This makes it possible to have several composition sets (I hope)
   if(associated(toopx%phres)) then
      phres=>toopx%phres
   else
      write(*,*)'3QX phres pointer is not assigned entering calc_toop'
      gx%bmperr=4399; goto 1000
   endif
! debug level 0 nothing, 1 minimum, 2 Toop, 5 all
   tkdeb=2
   if(tkdeb.ge.2) write(*,*)'3XQ in calc_toop ',lokpty%degreelink(0)
! NOTE vals, dvals and d2vals set to zero in calcg before calling this routine
   rtg=gz%rgast
   if(lokpty%degree.eq.0) then
! quick exit if no composition dependence
      lfun=lokpty%degreelink(0)
      call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
      if(gx%bmperr.ne.0) goto 1000
      if(lokpty%proptype.eq.1) then
         valtp=valtp/rtg
      endif
! this is multiplied with y_i y_j (and their derivatives) at the return
      vals=vals+valtp
      goto 1000
   endif
! we come here if there are RK terms >0
   if(tkdeb.gt.0) then
      write(*,10)gz%iq(1),gz%iq(2),lokpty%degree
10    format(/'3XQ in calc_toop & Kohler with binary;',2i3,' degrees: ',i2)
   endif
!   do nz=1,3
! it seems that dvals are not properly initiatiad to zero?
!      write(*,7)nz,(dvals(nz,ic),ic=1,gz%nofc)
!7     format('3XQ initial dvals: ',i2,10(1pe12.4))
!   enddo
! We have to calculate the reduced fractions, it can involve many fractions
   nyfr=size(phres%yfr)
   allocate(dsigma(nyfr))
   allocate(dx12(nyfr))
   allocate(dx21(nyfr))
! default value of sigma is unity
   sigma=one
! these are default zero, i.e. derivatives with respect to no extra fractions
   dx12=zero
   dx21=zero
   dsigma=zero
! constituents are ordered alphabetically, x12 is the first in the endmember
   x12=gz%yfrem(gz%intlat(1))
   x21=gz%yfrint(1)
   if(tkdeb.ge.2) write(*,15)x12,x21,gz%iq(1),gz%iq(2)
15 format('3XQ initial fractions: ',2f8.4,2i5)
! We have a binary excess parameter which depend on x_A and x_B
! and a Redlich-Kister polynom (x_A -x_B)/sigma
! When the data for the system was entered some ternaries were
! specified as Toop or Kohler and the toopx record created with the
! information needed for the calculations below
! For all ternaries A-B-K where the composition of B is constant (Toop)
! the fraction of K should be added to A, i.e. x12
   if(phlista(lokph)%toopfirst%endmemel.ne.0) then
! CHECK FOR DUPLICATE FRACTION INDICES, an add ternary may add same fraction!!
! phlista(lokph)%firsttoop%free=-1 in add_ternary... (in gtp3H.F90)
! if phlista(lokph)%firsttoop%free=-1 check and remove redundant fractions!!
! This phlista(lokph)%firsttoop%free=0 at the end of gcalc (in gtp3X.F90)
      if(tkdeb.ge.1) write(*,16)phlista(lokph)%toopfirst%endmemel
16    format('3XQ Checking duplicates as phlista(lokph)%toopfirst%endmemel:',i2)
! The check made only once, this value is zeroed at end of calcg subroutine
      allocate(ctoop1(phlista(lokph)%toopfirst%free))
      allocate(ctoop2(phlista(lokph)%toopfirst%free))
      allocate(ckohler(phlista(lokph)%toopfirst%free))
   endif
   not1=.TRUE.
   if(tkdeb.ge.2) then
      write(*,8)toopx%free,nyfr
8     format('3XQ Number of Toop/Kohler ternaries: ',i3,&
           ' Total number of fractions: ',i3)
      write(*,12)phres%yfr
12    format('3XQ All yfr: ',20F7.4)
   endif
! toopx%free is last used index in %Toop1, %Toop2 and %Kohler
   not1=.TRUE.;    not2=.TRUE.;    nok=.TRUE.
   allcorr: do ic=1,toopx%free
      if(tkdeb.ge.2) &
           write(*,33)ic,toopx%toop1(ic),toopx%toop2(ic),toopx%kohler(ic)
33    format('3XQ List of Toop/Kohler constituents: ',i2,2x,3i3)
!------------ Toop1
      cc=toopx%toop1(ic)
      if(allocated(ctoop1)) then
! if ctoop1 allocated then check to eliminate duplicates
         if(tkdeb.ge.1) write(*,'(a)')'3XQ Check for duplicated fractions'
         do nz=1,ic-1
            if(cc.gt.0 .and. cc.eq.toopx%toop1(nz)) then
               jdeg=toopx%toop1(ic); toopx%toop1(ic)=0; cc=0
               write(*,69)'Toop1',nz,jdeg
69             format('3XQ eliminated duplicate ',a,' fraction',2i4)
            endif
         enddo
      endif
      if(cc.gt.0) then
! In this binary i-j with ternary k where i (endmember) is constant (Toop)
! Add the fraction of x_k to x_i
         x12=x12+phres%yfr(cc); dx12(cc)=one; not1=.FALSE.
         if(tkdeb.ge.2) &
              write(*,34)'x12   ',ic,cc,toopx%toop1(ic),phres%yfr(cc),x12
34       format('3XQ Added fraction to ',a,3i3,2E15.7)
      endif
!------------ Toop2
      cc=toopx%toop2(ic)
      if(allocated(ctoop2)) then
! if ctoop2 allocated check to eliminate duplicates
         do nz=1,ic-1
            if(cc.gt.0 .and. cc.eq.toopx%toop2(nz)) then
               jdeg=toopx%toop2(ic); toopx%toop2(ic)=0; cc=0
               write(*,69)'Toop2',nz,jdeg
            endif
         enddo
      endif
      if(cc.gt.0) then
! In this binary i-j with ternary k where i (interaction) is constant (Toop)
! Add the fraction of x_k to x_j
         x21=x21+phres%yfr(cc); dx21(cc)=one; not2=.FALSE.
         if(tkdeb.ge.2) &
              write(*,34)'x21   ',ic,cc,toopx%toop2(ic),phres%yfr(cc),x21
      endif
!------------ Kohler
      cc=toopx%Kohler(ic)
      if(allocated(ckohler)) then
! if ckohler allocated check to eliminate duplicates
         do nz=1,ic-1
            if(cc.lt.0 .and. cc.eq.toopx%kohler(nz)) then
               jdeg=toopx%kohler(ic); toopx%kohler(ic)=0; cc=0
               write(*,69)'3Kohler',nz,jdeg
            endif
         enddo
      endif
      if(cc.lt.0) then
! In this ternary i-j-k the i-j extrapolates as Kohler
! the composition of k should be subtracted from sigma (initiated to 1.0 above)
         sigma=sigma-phres%yfr(-cc); dsigma(-cc)=-one; nok=.FALSE.
         if(tkdeb.ge.2) write(*,35)ic,cc,toopx%kohler(ic),phres%yfr(-cc),sigma
35       format('3XQ subtracted fraction for sigma ',3i3,2E15.7)
      endif
   enddo allcorr
   if(x21.ge.one) then
      write(*,*)'3XQ Error: x21 larger than 1.0 in Toop/Kohler extrapolation!'
      gx%bmperr=4399; goto 1000
   endif
   if(x21.ge.one) then
      write(*,*)'3XQ Error: x21 larger than 1.0 in Toop/Kohler extrapolation!'
      gx%bmperr=4399; goto 1000
   endif
   if(sigma.le.zero) then
      write(*,*)'3XQ Error: negative sigma in Toop/Kohler extrapolation!'
      gx%bmperr=4399; goto 1000
   endif
! This is the RK fraction difference, sigma is the Kohler divisor
   dxrk0=(x12-x21)/sigma
! dxrk is the Tredlich-Kister term, it is raised to powers jdeg=0...n
! The derivative of dxrk**n is:
!   n*dxrk**((m-1)*[ (dx12-dx21)/sigma - (x12-x21)*dsigma/sigma**2 ]
! where dx12, dx21 be 0 or 1 and dsigma 0 or -1 for several fraction variables
! were set above.  
!
! dxrk=1.0 for jdeg=0
   dxrk=one
   if(tkdeb.ge.2) then
      write(*,17)'3XQ x12:   ', x12,',   dx12:   ',dx12,dxrk0
      write(*,17)'3XQ x21:   ', x21,',   dx21:   ',dx21
      write(*,17)'3XQ sigma: ', sigma,', dsigma:   ',dsigma
17    format(a,F8.6,a,10F7.3)!
   endif
!-----------------------------------------------------------------
! No documentation of code below (at present), see paper by Pelton 2001
!-----------------------------------------------------------------
! in toopx there are 3 arrays
! toop1 with toop constitunents to be added to iq(1)
! toop2 with toop constitunents to be added to iq(2)
! Kohler with constitunents to be subtracted from sigma
! Calculate the corrected the binary fractions x12 and x21 and sigma
   if(tkdeb.gt.0) write(*,20)x12,x21,sigma,dxrk,moded
20 format('3XQ fractions: ',2F8.4,' sigma,dxrk: ',2F8.4,' moded: ',i1)
! gz%iq(1) is first constitution, gz%iq(2) in interaction
   dx12(gz%iq(1))=one/sigma
   dx21(gz%iq(2))=one/sigma
   RK: do jdeg=0,lokpty%degree
      lfun=lokpty%degreelink(jdeg)
      call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
      if(gx%bmperr.ne.0) goto 1000
      if(lokpty%proptype.eq.1) then
         valtp=valtp/rtg
      endif
      vals=vals+dxrk*valtp
      if(tkdeb.ge.2) write(*,9)'3XQ vals1: ',jdeg,dxrk,valtp(1),rtg*valtp(1),&
           vals(1),rtg*vals(1)
9     format(a,i2,5(1PE13.5))
      noder5: if(moded.gt.0) then
! moded=0 no derivative, =1 first, =2 second; gz%iq(1) is endmember
! derivatives with respect to original x12 and x12
! qz=1 is parameter value, qz=2 is parameter derivative wrt T, qz(3) wrt P
         do qz=1,3
            dvals(qz,gz%iq(1))=dvals(qz,gz%iq(1))+dx12(gz%iq(1))*valtp(qz)
            dvals(qz,gz%iq(2))=dvals(qz,gz%iq(2))-dx21(gz%iq(2))*valtp(qz)
         enddo
! derivatives wrt Toop1 constintuents, use dx12, dx21 and dsigma
! all approximate ....... negative sign of dx21 taken care of when "added"
! NOTE dx12, dx21 and sigma are arrays as any constituent can be involved
         dx12(gz%iq(1))=(jdeg+1)*dxrk
         dx21(gz%iq(2))=(jdeg+1)*dxrk
! This part takes care of derivatives wrt fractions "k" in x12, x21 and sigma
! They have dx12(k)=dx21(k)=1 and dsigma)k)=-1
! dxrk**n * valtp is the Redlich-Kister term, valtp(1,2,3) is the parameter
! The derivative of dxrk**n * valtp is:
! n*dxrk**((n-1)*valtp*[ dx12/sigma -dx21/sigma -(x12-x21)*dsigma/sigma**2 ]
! valtp(1) is parameter value, valtp(2,3) is derivative wrt T and P respectivly
! where dx12, dx21 are 0 or 1 and dsigma is 0 or -1 for the fraction variables
! The fractions "k" involved have nonzero %toop1(ic), %toop2 or %kohler indices
         extraderivatives: do ic=1,toopx%free
! the arrays %toop1, %toop2 and %kohler have the same dimensions
! they have fraction indices in toop1, toop2 or kohler (most of which is 0)
! -------------------- derivatives for toop1
            cc=toopx%toop1(ic)
            ltoop1: if(.not.not1) then
! there is a fraction added to x12, fraction index in toopx%toop1(ic)
               if(cc.gt.0) then
                  do qz=1,3
! this fraction is added to x12, dx12=1 but we have to divide with sigma
                     dvals(qz,cc)=dvals(qz,cc)+(jdeg+1)*dxrk*valtp(qz)/sigma
                  enddo
                  if(tkdeb.ge.2) write(*,44)'Toop1 ',cc,dvals(1,cc)
44                format('3XQ ',a,' derivative: ',i2,1pe14.6)
! Any second derivatives is ignored (it may slow down convergence)
               endif
            endif ltoop1
!--------------------- derivatives for Toop2
            cc=toopx%toop2(ic)
            ltoop2: if(.not.not2) then
! there is a fraction added to x21, fraction index in toopx%toop2(ic)
               if(cc.gt.0) then
                  do qz=1,3
! dx21(ic) is unity here but divide with sigma.  OBS negative sign
                     dvals(qz,cc)=dvals(qz,cc)-(jdeg+1)*dxrk*valtp(qz)/sigma
                  enddo
                  if(tkdeb.ge.2) write(*,44)'Toop2 ',cc,dvals(1,cc)
! Any second derivatives ignored (it may slow down convergence)
               endif
            endif ltoop2
!---------------------- derivatives for Kohler, negative index of fraction!!!
            cc=toopx%kohler(ic)
            lkohler: if(.not.nok) then
! there is a fraction subtracted from sigma, fraction -index in toopx%kohler(ic)
               if(cc.lt.0) then
                  if(tkdeb.ge.2) write(*,54)cc,jdeg,dvals(1,-cc),&
                       (jdeg+1)*dxrk*valtp(1)*(x12-x21)/sigma**2,&
                       dxrk,valtp(1),(x12-x21),sigma
54  format('3XQ Kohler derivative: ',2i2,2(1pe12.4)/4x,4(1pe12.4))
                  do qz=1,3
! dxrk**n * valtp is the Redlich-Kister term, valtp is the parameter
! n*dxrk**((n-1)*valtp*[ dx12/sigma -dx21/sigma -(x12-x21)*dsigma/sigma**2 ]
! dsigma is unity here but divide with sigma**2
                     dvals(qz,-cc)=dvals(qz,-cc)-&
                          (jdeg+1)*dxrk*valtp(qz)*(x12-x21)/sigma**2
                  enddo
! Any second derivatives ignored (it may slow down convergence)
               endif
            endif lkohler
         enddo extraderivatives
! dxrk has one more power for next term
         dxrk=dxrk*dxrk0
      endif noder5
   enddo RK
!--------- maybe almost finished ???
!
!   if(tkdeb.ge.1) write(*,30)'3XQ vals2: ',vals(1),rtg*vals(1),&
!        gz%iq(1),gz%iq(2),rtg*dvals(1,gz%iq(1)),rtg*dvals(1,gz%iq(2))
30 format(a,2F12.4,2i2,2F12.4)
1000 continue
!------------------------------------------------------------------
! this calculates the whole  \sum_i (\xi_A - \xi_B)/sigma_AB)^i iL_AB
! and derivatives ....
!------------------------------------------------------------------
! The result is multiplied with the fractions x_A'x_B in the calling routine
   if(tkdeb.gt.0) write(*,'(a,2i3,F12.4)')'3XQ RT*vals: ',&
        gz%iq(1),gz%iq(2),rtg*vals(1)
!  if(tkdeb.gt.0) write(*,'(a,i3,2x,5F8.5)')'3XQ dxrk mm:',&
!       jdeg,rtg*vals(1),dxrk0,dxrk
   return
 end subroutine calc_toop

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

