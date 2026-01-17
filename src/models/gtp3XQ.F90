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
! test of indexing problem
   integer line757
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
         write(*,212)s1,mqmqa_data%contyp(10,s1)
212      format('3XQ Warning: mqmqa_data%contyp(10,s1) =/= s1:',2i4)
! emergecy fix 17/12 2025 does not work 
!         mqmqa_data%contyp(10,s1)=s1
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
! ixsym finds the sequential storage place of (i,j) in a symmetrical array
!   write(*,538)ncon,ixsym(ncon,ncon),ixsym(5,3),nspin
538 format('3XQ entropy: ',3i5,' nspin: ',20i3)
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
!      write(*,207)nspin(1),nspin(2),ncon,noofpair
207   format('3XQ allocating phvar%mqmqaf arrays',2i3,4i5)
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
!      write(*,*)'3XQ allocation of d2yy2:',size(phvar%mqmqaf%d2yy2)
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
!  write(*,771)nspin(1),nspin(2),ncon
   do s1=1,nspin(2)
      phvar%mqmqaf%yy2(s1)=yy2(s1)
      phvar%mqmqaf%ceqf2(s1)=ceqf2(s1)
      do s2=1,ncon
         phvar%mqmqaf%dyy2(s1,s2)=dyy2(s1,s2)
         phvar%mqmqaf%dceqf2(s1,s2)=dceqf2(s1,s2)
      enddo
!****************************************************************
!      write(*,*)'3XQ line 757 skipping a 2nd derivative'
!****************************************************************
      do s3=s1,ncon
         s8=ixsym(s3,s1)
! large dimension problem here ixsym is a function to access a symetric array
!         write(*,671)s1,s3,s8,ixsym(s1,s8),size(d2yy2)
671      format('3XQ accessing d2yy2: ',3i4,2i7)
!         phvar%mqmqaf%d2yy2(s1,s8)=d2yy2(s1,ixsym(s1,s8))
         line757=max(line757,s1*ixsym(s1,s8))
      enddo
   enddo
!   write(*,*)'3XQ line 771: ',line757,s1*ixsym(s1,s8)
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
! Here we calculate for one formula unit (FU) of the phase
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
! TEST: the value in contyp(10,q1) should be q1 ...  260111/BoS WHY??
         write(*,441)q1,mqmqa_data%contyp(10,q1),mqmqa_data%contyp(14,q1)
441      format('3X problems in %contyp with quad indexing:',3i5)
!         gx%bmperr=4399; goto 1000
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
   integer ij,jd,jq,qq1,qq2,ass,mpow,isumx,tsize,tch,iiz,mqmqcon,mqmqjy
   integer noofex,nqx,ncv,icv
   double precision ksi,sumx,dsumx
   double precision dksi(3),d2ksi(3)
!   logical ddebug
!------------------------------------- 
! tch is level of debug output, 0=none, 3=max
   tch=0
   noofex=0
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
!17 format(a,i3,4(1pe12.4))
!--------------------------------------
! Trying to understand the data structure.  List all species and some data
!   do mqmqj=1,noofsp
!      write(*,13)mqmqj,splista(mqmqj)%symbol,splista(mqmqj)%alphaindex,&
!           splista(mqmqj)%quadindex
!13    format('3XQ specie: ',i3,2x,a,2x,5i5)
!   enddo
!--------------------------------------
! debug output of varkappa mm moved to beginning of calc_mqmqa
   if(mqmqxcess .and. btest(phlista(lokph)%status1,PHMQMQX)) then
      write(*,*)'3XQ Debug output of quads, \varkappa_ij, \xi_ij and y_i/k'
!
! these variables are in the TYPE GTP_MQMQA_VAR
      nqx=mqmqa_data%nquad
      write(*,82)nqx
82    format('3XQ Quad fractions:',i3)
      write(*,84)(mqf%xquad(icv),icv=1,nqx)
84    format((8F8.5))
      ncv=size(mqf%compvar)
      write(*,78)ncv
78 format('3XQ varkappa_ij     varkappa_ji          xi_ij           xi_ji',i12)
!          123456789.123456123456789.123456.....123456789.123456123456789.123456
      do icv=1,ncv
         write(*,80)mqf%compvar(icv)%vk_ij,mqf%compvar(icv)%vk_ji,&
              mqf%compvar(icv)%xi_ij,mqf%compvar(icv)%xi_ji
80       format(2x,2(1pe16.8),5x,2(1pe16.8))
      enddo
      write(*,86)mqmqa_data%ncat,(mqf%y_ik(icv),icv=1,mqmqa_data%ncat)
86    format(/'3XQ y_i/k ',i2,': ',(7F9.6))
   endif
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
!   write(kou,299)
299 format('3QX endmember energy and entropy calculated, excess to be done')
!   goto 800
!---------------------------------------------------------------------
! code below needed for excess parameters ONLY, all SNN FNN endmembers done
! NOTE some of them may not have a reference energy parameter
! This is to allocate csumx for handling quads with small fractions.
   isumx=0
! debug output of G for check of excess
   if(mqmqxcess) then
      write(*,288)(phres%gval(itp,1),itp=1,4)
288   format('3XQ before excess:'/'G, dG/dT dG/dP d2G/dT2:',4(1pe14.6))
   endif
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
      if(tch.ge.3) write(*,311)mqmqj,mqmqa_data%nconst,kend,&
           associated(endmemrec%intpointer)
311   format(/'3XQ in loop for excess parameters: ',3i5,l2)
      intrec=>endmemrec%intpointer
! interaction parameters are NOT linked from SNN endmembers ?? really?
! They are stored in alphabetical order of the constituents
!      write(*,*)'3XQ Check interaction parameters 1',associated(endmemrec),&
!           associated(intrec),mqmqj,kend
! if we cycle here the results are the same as without excess parameters
!      cycle endmemloop2
!
      if(.not.associated(intrec)) then
         cycle endmemloop2
      endif 
      if(.not.btest(phlista(lokph)%status1,PHMQMQX)) then
! this is the first MQMQA implementation with correct reference energy and
! configurational entropy but very messy Toop/Kohler implementation
         goto 499
      endif
!
!      if(mqmqxcess) then
! if parameter errors in interactions below these are the endmemberquads A/X
!         write(*,313)(mqmqa_data%emquad(iiz),iiz=1,mqmqa_data%ncat)
313      format('3XQ endmember quads: ',15i3)
!      endif
!  
! mqmqj is NOT the mqmqa constituent index, it is just an endmember counter
! look for the constituent in fraction record, sublattice 1, constituent 1
! WOW !!!! it does not crash
      mqmqjy=endmemrec%fraclinks(1,1)
!
! we must find its position in the quad list      
!      write(*,314)mqmqj,mqmqjy,size(phlista(lokph)%constitlist)
!           associated(endmemrec%oendmemarr),associated(endmemrec%dendmemarr)
314   format(/'3XQ endmember data: ',3i3)
!      write(*,315)phlista(lokph)%constitlist
315   format('3XQ constituents: ',20i3)
!      
      if(mqmqxcess) write(*,318)phres%gval(1,ipy),&
           (phres%dgval(1,jq,ipy),jq=1,gz%nofc)
!
      noofex=noofex+1
      if(mqmqxcess) write(*,*)'3XQ excess with endmember constituent: ',&
           mqmqjy,ipy
      call new_mqmqa_excess(lokph,intrec,mqmqjy,vals,dvals,d2vals,gz,ceq)
      if(gx%bmperr.ne.0) goto 1000
! intrec is nullified inside new_mqmqa_excess
      if(mqmqxcess) write(*,*)'3XQ back with excess from endmember ',mqmqjy,&
           associated(endmemrec)
!      if(mqmqxcess) write(*,317)gz%nofc,vals(1),vals(2),&
!      write(*,317)gz%nofc,rtg*vals(1),rtg*vals(2),&
!           (rtg*dvals(1,jq),jq=1,gz%nofc)
317   format('3XQ Back from new_mqmqa:  ',i3,2(1pe12.4)/6(1pe12.4))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
! BIG STEP ... add vals, dvals to G and dG/dy
! what about aff?
!
! what is ipy? It is the property 1 is G, 2 is BMAGN or something else
      ipy=1
      if(mqmqxcess) write(*,318)phres%gval(1,ipy),&
           (phres%dgval(1,jq,ipy),jq=1,gz%nofc)
318   format('3XQ G & G.y: ',1pe12.4/6(1pe12.4))
      do itp=1,6
! loop for G, G.T, G.P, G.T.T, G.T.P, G.P.P to add excess contribution
         phres%gval(itp,ipy)=phres%gval(itp,ipy)+vals(itp)
      enddo
! TEMPORARILY REMOVED SOME LOOPS
      do jq=1,gz%nofc
! skip loop for dG/dy, d2G/dydT, d2G/dydP, only for constituents
         do itp=1,3
            phres%dgval(itp,jq,ipy)=phres%dgval(itp,jq,ipy)+dvals(itp,jq)
         enddo
      enddo
!      
! ignore 2nd derivatives as not calculated for excess
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!
      if(.not.associated(intrec)) cycle endmemloop2
!
!*************** remove all code below when excess code above OK *******
!******************* new code above should replace code below *******
!
! There are excess parameters, any Tooprecords?
499      continue
!      write(*,319)
319   format(//'3XQ *** this is the old mqmqa excess model **'//)
      if(associated(intrec%tooprec)) then
! the allocatable arrays Toop1, Toop2 and Kohler have all same size
! equal to the number of binary combination of constituents
!
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
         write(*,*)'3XQ MQMQA parameter with several properties!',mqmqj
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
!   write(*,990)'3XQ exit calc_mqmqa G:',rtg*phres%gval(1,1),rtg*vals(1)
!        (phres%dgval(1,s1,1),s1=1,gz%nofc)
990 format(a,5(1pe14.6))
1000 continue
   return
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
!
! new MQMQA excess subroutines below
! using a separate data structury for asymmetries
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine new_mqmqa_excess
! called from calc_mqmqa line 1429.  CALCULATES MQMQA excess
 !\begin{verbatim}
 ! subroutine new_mqmqa_excess(lokph,intrecin,mqmqj,vals,dvals,d2vals,gz,ceq)
 subroutine new_mqmqa_excess(lokph,intrecin,mqmqj,vals,dvals,d2vals,gz,ceq)
! vals(1..6) are G, dG.T, dG.P, d2G.T.T, d2G.T.P and d2G.P.P for parameter
! dvals(1,i) are first derivatives wrt fracton and 2nd wrt fraction, T or P
!          dvals(1,i) is dG.yi, dval2(2,i) is d2G.yi.T, dvals(3,i) is d2G.yi.P
! d2vals(i,j) are second derivatives to 2 fractions, IGNORED HERE
! gz%nofc is number of fraction variables multiplied with this parameter(?)
!
! To be written using the gtp_allinone data structure for asymmetric excess
   implicit none
! mqmqj is index of first constituent in endmemberrecord
   integer lokph,mqmqj
!   type(gtp_property), pointer :: lokpty
   type(gtp_parcalc) :: gz
   type(gtp_phase_varres), pointer :: phres
   TYPE(gtp_mqmqa_var), pointer :: mqf
   double precision vals(6),dvals(3,gz%nofc)
   double precision d2vals(gz%nofc*(gz%nofc+1)/2)
! pointer to first interaction record from an endmember
   TYPE(gtp_interaction), pointer :: intrecin,intrec
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! needed locally?
   TYPE(gtp_intstack), dimension(:), allocatable :: savedint
   TYPE(gtp_pystack), pointer :: pystack
   TYPE(gtp_phase_add), pointer :: addrec
   TYPE(gtp_terdata), pointer :: ternaries
   TYPE(gtp_property), pointer :: proprec
!   type(gtp_allinone), pointer :: compvar
!
   character*120 text
   double precision :: rtg
   logical :: once=.true.
   integer ppow,qpow,rpow,intlev,iiz,jj,pairquad,jp
   integer, save :: proprecno=0
   integer parquad(4),nprr,nfr
   integer :: nex=0
   integer, dimension(:), allocatable :: ylinks,qlinks
   integer ncv,icv,nqx,lokcs,lokfun,xq,cxq,mm
   double precision mqmqx_deltag,compprod,nomin,divisor,ternary
!
   logical, save :: ternaryonce=.true.
!
! composition derivatives are only relative to quads  !!!!!!!
! the composition variables for a parameters are expressions using asymmetric
! y_ik, \xi or \varkappa which depend on quads
! we have to sort out how this affects the derivatives
! dy_ik are factors for y_ik relative to quads, can be 1 or less
! dxi_ij and dx_ji and  dvk_ij and dvk_ji are 1 or less
! a parameter P multiplied with vk_ij(ij) has several contributions to the
! derivatives dP(zz), dvk_ij(ij,zz),zz=1,nquad
  integer idyix(5,mqmqa_data%nquad)
  integer zkij
  double precision haha,one1,dnomin,ddivisor,dternary
  double precision dyix(5,mqmqa_data%nquad)
!
   character*1 ptyp1
! The previous MQMQA excess implementation arrive here
! If mqmqa_data%exlevel is zero we should return and old code will still work.
   if(mqmqa_data%exlevel.eq.0) then
!      if(once) write(*,6)mqmqa_data%exlevel
6     format('3XQ *** this system use the old excess model ***',i5)
      goto 1000
   endif
! we are here because this endmember has an intercation link
!   if(mqmqxcess) write(*,5)mqmqj
!   write(*,5)mqmqj
5  format(/'3XQ in new_mqmqa_excess with endmember: ',i3)
! initiate ylinks for this tree with the endmember fraction
   intrec=>intrecin
! this is needed to move to next endmemeber
   nullify(intrecin)
!
! divide values with rtg?
   rtg=globaldata%rgas*ceq%tpval(1)
!---------------------------------
! not more than 10 interactions ....
   allocate(savedint(10))
   allocate(ylinks(10))
   allocate(qlinks(10))
! this is the endmember constituent
   nfr=1
   ylinks(1)=mqmqj
! there can only a one quad with two cations (pair) in an interaction
   pairquad=0
!   ifem: do jj=1,mqmqa_data%ncat
!      if(ylinks(1).eq.mqmqa_data%emquad(jj)) goto 17
!   enddo ifem
! this quad is evidently a pair AB/X
!   pairquad=ylinks(1)
17 continue
!
! THIS IS THE CALCULATION ROUTINE WITH extensive DEBUG LISTING ADDED
!
! loop here until all excess records from this endmember calculated
! The new excess model implementation using allinone etc below
! The quad fractions and related composition variables such as
! quadfractions and asymmetrical variables  have been set by set_constitution
!
! ceq%phase_varres(lokcs)%mqmqaf%compvar(icv)%vi_ij etc
! access to composition variables
!
   lokcs=phlista(lokph)%linktocs(1)
   mqf=>ceq%phase_varres(lokcs)%mqmqaf
!
!   if(btest(phlista(lokph)%status1,PHMQMQX)) then
   if(.false.) then
! these variables are in the TYPE GTP_MQMQA_VAR
!      write(*,*)'3XQ this is a check we have access to xquad, compvar etc?'
!
      nqx=mqmqa_data%nquad
      write(*,82)nqx
82    format('3XQ Quad fractions:',i3)
      write(*,84)(mqf%xquad(icv),icv=1,nqx)
84    format((8F8.5))
      ncv=size(mqf%compvar)
      write(*,78)ncv
78 format('3XQ varkappa_ij     varkappa_ji          xi_ij           xi_ji',i12)
!          123456789.123456123456789.123456.....123456789.123456123456789.123456
      do icv=1,ncv
         write(*,80)mqf%compvar(icv)%vk_ij,mqf%compvar(icv)%vk_ji,&
              mqf%compvar(icv)%xi_ij,mqf%compvar(icv)%xi_ji
80       format(2x,2(1pe16.8),5x,2(1pe16.8))
      enddo
      write(*,86)mqmqa_data%ncat,(mqf%y_ik(icv),icv=1,mqmqa_data%ncat)
86    format(/'3XQ y_i/k:',i2,3x,7F9.6)
!
      write(*,*)'Press return to continue'
      read(*,*)
   endif
!
! when we are here intrec must be associated
   intlev=0
   intloop: do while(associated(intrec))
!
! This has to be developed, intrec must be associated here
!
! there is a single set of sites, save constituent first index
! We may come back here for another interaction with same endmember
! Set ylinks to be indices of the OC fractions
!      write(*,*)'3XQ Starting intloop with component ',intlev,intrec%fraclink
! save name of interacting constituent even if no property
      nfr=nfr+1
      ylinks(nfr)=intrec%fraclink(1)
      proprec=>intrec%propointer
! loop for all property records
      listprop: do while(associated(proprec))
! we have found an excess parameter !!!
! there can be several property record for the same set of constituents
! there is a propointer here! error in gtp3B or empty interaction record
         proprecno=proprecno+1
         if(proprec%proptype.eq.34) ptyp1='G'
         if(proprec%proptype.eq.35) ptyp1='Q'
         if(proprec%proptype.eq.36) ptyp1='B'
!
         ternary=1.0D0
         ppow=proprec%asymdata%ppow
         qpow=proprec%asymdata%qpow
         rpow=proprec%asymdata%rpow
         if(mqmqxcess) then
! helps to understand what the parameter it is ....
            jp=1
            text=' '
            call mqmqa_excesspar_name(lokph,intlev,nfr,ylinks,text,jp)
            text(jp-1:)=';'//ptyp1//','//char(ichar('0')+ppow)//&
                 ','//char(ichar('0')+qpow)//','//char(ichar('0')+rpow)//')'
            write(*,115)trim(text),ppow,qpow,rpow
115         format(/'3XQ param: ',a,', pqr:',3i2)
! extract the quad pointers
         endif
!
         lokfun=proprec%degreelink(0)
         xq=proprec%asymdata%quad
! cxq transforms the quad index to an index in compvar (which as no diagonal)
         cxq=mqmqa_data%quad2compvar(xq)
         par3: if(nfr.gt.3) then
            if(ternaryonce) write(*,116)
116         format(/'3XQ *** ternary parameters not yet implemented ***'/)
            ternaryonce=.false.
            mqmqx_deltag=0.0d0
            vals=0.0d0
            dvals=0.0d0
            goto 800
! code below skipped            
! note ylinks(nfr) is not the correct ternary fraction, it should be the quad
! which is neither equal to xq, or the two quads in vk_ij or x_ij
            call ternary_factor(xq,mqf%compvar(cxq)%cat1,mqf%compvar(cxq)%cat2,&
                 ylinks,mm,ternary,proprec)
            if(gx%bmperr.ne.0) goto 1000
! mm is the index of the cation in the C/X quad, the factor is
!  y_ik/(\xi_ij) if mm is \gamma in i-j-\gamma or
!  y_ik/(\xi_ji) if mm is \nu in    i-j-\nu or just
!  y_ik          if neither
! we can check this in compvar(cxq) .... later
            write(*,104)mm
104         format('3XQ *** no check if the ternary i-j-',i2,' is asymmetric')
            ternary=mqf%y_ik(mm)
!            write(*,106)mm,ternary
106         format('3XQ *** ternary factor "y_m/(\xi_ij)", m=',i2,5x,1pe13.6)
         endif par3
! we cannot use quad index in compvar.  We should use the index
! related to the two elemnts i and j represented by the quad index
! cxq=quad2comvar(xq) converts the xq index for AB/X to the two quads A/X B/X 
!         if(mqmqxcess) write(*,117)xq,cxq,nfr,intlev,lokfun,&
         if(mqmqxcess) then
            if(ptyp1.eq.'G') then
               write(*,117)xq,cxq,nfr,lokfun,mqf%xquad(xq),&
                    mqf%compvar(cxq)%vk_ij,mqf%compvar(cxq)%vk_ji,ternary
117            format('3XQ quad, compvar:',2i3,', nfr ',i3,&
                    ', parameter link: ',i5/&
                    '3XQ Quad ',1pe11.4,' \vk_ij ',1pe11.4,' \vk_ji',1pe11.4,&
                    ' tern ',1pe11.4)
            elseif(ptyp1.eq.'Q') then
               write(*,112)xq,cxq,nfr,lokfun,mqf%xquad(xq),&
                    mqf%compvar(cxq)%xi_ij,mqf%compvar(cxq)%xi_ji,ternary
112            format('3XQ quad, compvar:',2i3,', nfr ',i3,&
                    ', parameter link: ',i5/&
                    '3XQ Quad ',1pe11.4,' \xi_ij ',1pe11.4,' \xi_ji',1pe11.4,&
                    ' tern ',1pe11.4)
            else
               write(*,*)'3XQ Unknown equation'
            endif
         endif
! this return 6 values in vals: G, G.T, G.P,  G.T.T, G.T.P, G.P.P
         call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
         if(gx%bmperr.ne.0) goto 1000
!--------------------------------------------------------------------
! divide all parameter values with rtg!!
         vals=vals/rtg
!--------------------------------------------------------------------
! calculate the parameter with composition variables
         if(ptyp1.eq.'G') then
! ppow is for ij or ji ? USING varkappa
!            write(*,*)'3XQ ppow: ',ppow,qpow
            nomin=(mqf%compvar(cxq)%vk_ij)**ppow*&
                 (mqf%compvar(cxq)%vk_ji)**qpow
!
            compprod=mqf%xquad(xq)*nomin*ternary
            mqmqx_deltag=compprod*vals(1)
!            write(*,120)xq,cxq,ppow,qpow,mqf%xquad(xq),&
!                 mqf%compvar(cxq)%vk_ij,mqf%compvar(cxq)%vk_ji,vals(1),&
!                 mqmqx_deltag
120         format('3XQ Quad, compvar :',2i3,', powers ',2i3,' quad ',1pe11.4/&
                 ' \vk_ij ',1pe11.4,' \vk_ji',1pe11.4,' para ',1pe11.4,&
                 ' Delta G ',1pe12.5)
            if(mqmqxcess) write(*,118)mqf%xquad(xq),nomin,1.0d0,compprod
118         format('3XQ comprod: ',1pe15.8,'*',1pe15.8,'/',1pe15.8,' =',1pe15.8)
! partial derivative begin, we have to loop for all quads <<<<<<<<<<<<<<<<<<<<
            do zkij=1,mqmqa_data%nquad
! loop for all quads, one added if zkij=xq; cxq is index of varkappa
! All mqf%compvar(cxq)%dvk_ij are already calculated ??!!  one is global 1.0D0
               one1=0.0d0
               if(zkij.eq.xq) then
!                  one1=one/xq
                  one1=one/mqf%xquad(xq)
               endif
! the partial derivative of quad*(\vk_ij**ppow)*(\vk_ji)**qpow*L
! (dq/dz + pdvk_ij/dz + qdvk_ji/dz )*(vk_ij**ppow)*(vk_ji**qpow)*L
!
!               haha=one1/xq+&
               haha=one1+&
                    ppow*mqf%compvar(cxq)%dvk_ij(xq)/mqf%compvar(cxq)%vk_ij+&
                    qpow*mqf%compvar(cxq)%dvk_ji(xq)/mqf%compvar(cxq)%vk_ji
!
!               write(*,308)haha,one1,ppow,mqf%compvar(cxq)%dvk_ij(xq),&
!                    mqf%compvar(cxq)%vk_ij,qpow,mqf%compvar(cxq)%dvk_ji(xq),&
!                    mqf%compvar(cxq)%vk_ji
308            format('3xQ haha: ',2(1pe11.3),i2,2(e11.3),i2,2(e11.3))
!
               if(mqmqxcess) then
                  write(*,310)zkij,xq,mqf%compvar(cxq)%dvk_ij(xq),&
                       mqf%compvar(cxq)%dvk_ji(xq),haha,vals(1),vals(2)
310               format('3XQ dvk_ij/ji: ',2i3,5(F8.4))
               endif
! dvals(1,xq) is derivative wrt zkij, dvals(2,xq) is derivative wrt zkij and T
               dvals(1,zkij)=dvals(1,zkij)+(haha-one1)*vals(1)
               dvals(2,zkij)=dvals(2,zkij)+(haha-one1)*vals(2)
            enddo
! partial derivative end >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         elseif(ptyp1.eq.'Q') then
            nomin=(mqf%compvar(cxq)%xi_ij)**ppow*(mqf%compvar(cxq)%xi_ji)**qpow
            divisor=(mqf%compvar(cxq)%xi_ij+mqf%compvar(cxq)%xi_ji)**(ppow+qpow)
            compprod=mqf%xquad(xq)*nomin*ternary/divisor
            mqmqx_deltag=compprod*vals(1)
            do zkij=1,mqmqa_data%nquad
! loop for all quads has to be added
               one1=1.0d0
               if(zkij.eq.xq) then
!                  one1=one/xq
                  one1=one/mqf%xquad(xq)
               endif
! excess is quad(xq)*\xi_ij**ppow\xi_ji**qpow/(\xi_ij+\xi_ji)**(ppow+qpow)
! next line is the derivative of compprod, nomin and divisor above wrt quad(xq) 
         dnomin=(ppow*mqf%compvar(cxq)%dxi_ij(xq)/mqf%compvar(cxq)%xi_ij+&
                 qpow*mqf%compvar(cxq)%dxi_ij(xq)/mqf%compvar(cxq)%xi_ji)*nomin
               ddivisor=(ppow+qpow)*(mqf%compvar(cxq)%dxi_ij(xq)+&
                                     mqf%compvar(cxq)%dxi_ji(xq))*divisor
               haha=one1+(dnomin/divisor-nomin*ddivisor/divisor**2)
               dvals(1,zkij)=dvals(1,zkij)+(haha-one1)*vals(1)
               dvals(2,zkij)=dvals(2,zkij)+(haha-one1)*vals(2)
            enddo
!            write(*,*)'3XQ \xi excess derivatives missing'
         else
            write(*,*)'3XQ proptyp B is not implemented: ',ptyp1
            stop
         endif
! vals is the returned value
         vals(1)=compprod*vals(1)
         vals(2)=compprod*vals(2)
         if(mqmqxcess) write(*,118)mqf%xquad(xq),nomin,divisor,compprod
! we have to add values to the derivatives as well as function
800      continue
         if(mqmqxcess) write(*,119)proprec%antalprop,(ylinks(jj),jj=1,nfr)
         if(mqmqxcess) write(*,121)vals(1),compprod,mqmqx_deltag
!         write(*,121)vals(1),compprod,mqmqx_deltag
119      format('3XQ Delta G to be added to G, propid:',i4,', yix:',10i4)
121      format('3XQ param ',1pe15.8,', comprod ',1pe15.8,', Delta G ',1pe15.8)
!---------------------------------------------------------------------
! there can be several parameters for the same set of constituents
! They differ by the powers but not the constituents
!         
         proprec=>proprec%nextpr
         if(mqmqxcess) then
            if(.not.associated(proprec)) then
               write(*,*)'3XQ no more properties'
            else
               write(*,*)'3XQ one more property!'
            endif
         endif
         nex=nex+1
      enddo listprop
!      write(*,*)'Press return to continue',intrec%fraclink(1)
!      read(*,*)
!
! we have listed all property records for these constituents
      push_ornext: if(associated(intrec%highlink)) then
! first try a higher level of interaction but save link to next
         intlev=intlev+1
         if(associated(intrec%nextlink)) then
!            write(*,97)intlev,intrec%nextlink%fraclink
97          format('3XQ saved nextlink at intlev: ',2i3)
            if(intlev.gt.9) then
               write(*,*)'Interaction level record overflow',intlev
               gx%bmperr=4399; goto 1000
            endif
            savedint(intlev)%saved=>intrec%nextlink
         else
            nullify(savedint(intlev)%saved)
         endif
         intrec=>intrec%highlink
      else
         intrec=>intrec%nextlink
! too many constituents ...
         nfr=nfr-1
         pop: do while(.not.associated(intrec))
!            write(*,*)'3XQ pop stack',intlev,nfr
            if(intlev.gt.0) then
               intrec=>savedint(intlev)%saved
               intlev=intlev-1
               nfr=nfr-1
            else
               exit intloop
            endif
!            if(associated(intrec)) &
!                 write(*,*)'3XQ take nextlink ',intrec%fraclink(1)
         enddo pop
         cycle intloop
      endif push_ornext
   enddo intloop
!
!------------ return to next endmember
   
1000  continue
   if(mqmqxcess) write(*,1001)proprecno
1001 format('3XQ Back to endmemeber record',i5)
   return
 end subroutine new_mqmqa_excess

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine convert_y2quadx
!\begin{verbatim}
 subroutine mqmqa_excesspar_name(lokph,intlev,nfr,ylinks,text,jp)
!
   integer lokph,intlev,nfr
   integer ylinks(*),jp
   character text*(*)
! write the set of constituents, complex as ylinks are quadindex
! and ylinks are constiuent order (which may be the same but not always)
!\end{verbatim}
   integer ii,jj,kk
   character*24, dimension(10) :: const
! strange intlev is 1 here, it is 0 in calling routine ... only endmember quad
!   write(*,10)nfr,(ylinks(ii),ii=1,nfr)
!10 format('3XQ *** no of const: ',i2,', ylinks: ',10i3)
!   write(*,20)phlista(lokph)%constitlist
!20 format('3XQ *** phase const: ',25i3)
! phlista(lokph)%constitlist is index in splista <<<<<<<<<<<<<<<
!   write(*,30)trim(splista(phlista(lokph)%constitlist(1))%symbol)
!30 format('3XQ *** phase const names: ',a)
!
! A useful excersize to remember how data in OC are stored!!!
!
   text='G(MSCL,'; jp=8
   do jj=1,nfr
    text(jp:)=trim(splista(phlista(lokph)%constitlist(ylinks(jj)))%symbol)//','
    jp=len_trim(text)+1
   enddo
!   write(*,40)(trim(splista(phlista(lokph)%constitlist(ylinks(jj)))%symbol),&
!        jj=1,intlev)
!        trim(splista(phlista(lokph)%constitlist(ylinks(3)))%symbol)
!40 format('3XQ *** phase const name: ',10(a,','))
!
   return
 end subroutine mqmqa_excesspar_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine ternary_factor
!\begin{verbatim}
 subroutine ternary_factor(xq,cat1,cat2,ylinks,mm,hejhopp,proprec)
! calculates the ternary factor of a parameter
   integer xq,cat1,cat2,mm,ylinks(*)
   double precision hejhopp
   type(gtp_property), pointer :: proprec
!\end{verbatim}
! xq is the AB/X quad index
! cxq in the index in compvar which gives 2 quad indices for A/X and B/X
! ylinks are the OC fraction indices
! mm is the unknown 4th quad
! hejhopp is the value to return, possibly 1.0D0
   integer ii
!   write(*,'(a,3i3,2x,10i3)')'3XQ trying to find mm',xq,cat1,cat2,&
!        (ylinks(ii),ii=1,4)
! but ylinks are OC fraction indices, not necessarily same as quad indices
! BUT at present, check which one of the last 2 in ylinks that is an A/X quad
   do ii=1,size(mqmqa_data%emquad)
      if(ylinks(3).eq.mqmqa_data%emquad(ii)) goto 100
   enddo
   do ii=1,size(mqmqa_data%emquad)
      if(ylinks(4).eq.mqmqa_data%emquad(ii)) goto 100
   enddo
   write(*,*)'3XQ cannot find the ternary C/X quad'
   gx%bmperr=4399; goto 1000
! return the index of the cation in the C/X quad
100 mm=ii
1000 continue
   return
 end subroutine ternary_factor

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine convert_y2quadx
!\begin{verbatim}
 subroutine convert_y2quadx(sem,nint,jord,pquad)
! This is to fix the constitution variables for an MQMQX excess parameter.
! It has one AB/X quad index and 2 A/X and B/X quads and possibly a C/X one
! convert y fraction indexed in sem and jord to quad indices in parquad
!
! I am really really fedup with this model
!
   implicit none
   integer sem,nint,jord(2,*),pquad(*)
!\end{verbatim}
   integer ii,jj,nq3,kk,pair,qorder(4),lowa,highb,temp(4),nbx
!
! input data from database
!   write(*,*)'3XQ *** fixing MQMQA parameter composition variables',&
!        size(mqmqa_data%emquad),size(mqmqa_data%con2quad)
!   write(*,5)sem,jord(2,1:nint)
5  format('3XQ fixing MQMQA parameter composition variables',10i3)
!   write(*,20)'emquad',(mqmqa_data%emquad(jj),jj=1,mqmqa_data%ncat)
!   write(*,20)'con2quad',(mqmqa_data%con2quad(jj),jj=1,mqmqa_data%nquad)
20 format('3XQ ',a,21i3)
! pquad(1) should be the pair quad AB/X among sem, jord(2,1..nint)
   temp(1)=mqmqa_data%con2quad(sem)
! pquad(2) should be the alphabetically first  in quad AB/X, i.e A/X
! pquad(3) should be the alphabetically second in quad AB/X, i.e B/X
! pquad(4) should be the 4th quand, not including A or B 
   temp(2)=mqmqa_data%con2quad(jord(2,1))
! one may have jord(2,2)=0 here if vacancies
   if(jord(2,2).eq.0) then
      write(*,*)'3XQ Vacancy not allowed in MQMQA quad'
      gx%bmperr=4399; goto 1000
   endif
   temp(3)=mqmqa_data%con2quad(jord(2,2))
! the quad index related to temp(2) and temp(3) should be temp(1) ???
! check:
! emquad have the quad indices of all A/X quads, there are ncat of them.
! The index of a quad (i,j) where j>i is emquad(i)+j-i   
!   write(*,*)'3XQ values temp: ',temp(1),temp(2),temp(3)
   if(temp(2).gt.temp(3)) then
      write(*,*)'3XQ parameter has wrong order of A/X and B/X quads'
      stop 76
!   else
!      ii=mqmqa_data%emquad(temp(2))+temp(3)-temp(2)
!      write(*,*)'3XQ values of mixed quad index: ',temp(1),temp(2),temp(3),ii
!      if(temp(1).ne.ii) then
!         write(*,*)'3XQ problems with quad indices'
!         stop 77
!      endif
   endif
   nq3=3
   if(nint.eq.3) then
      nq3=4; temp(4)=mqmqa_data%con2quad(jord(2,3))
   endif
!   write(*,10)'first',(temp(ii),pquad(ii),ii=1,nq3)
10 format('3XQ ',a,' quads ',2i3,', first ',2i3,', second ',2i3,', maybe ',2i3)
! find the AB/X quad and the arrange the others
! all but one of the quads in temp(1..nq3) should be A/X quads
! and temp(2) should have the lowest index of the AB/x quad and temp(3)
! the highest.  Any temp(4) quad should not be A/X or B/X
! this code is horrible
!
   pair=0; lowa=0
   loop4: do ii=1,nq3
      pquad(ii)=temp(ii)
      qorder(ii)=ii
!
! qx is quad
! Calculate: quad(temp(1)*\xi(temp(3),temp(2)**ppow
!      write(*,*)'3XQ is this the pair?',pquad(ii),qorder(ii)
      loopax: do jj=1,mqmqa_data%ncat
!
! cycle loop4 if temp(ii) is an A/X quad
!
         if(temp(ii).eq.mqmqa_data%emquad(jj)) cycle loop4
      enddo loopax
! if we arrive here temp(ii) is a AB/X quad
      if(pair.eq.0) then
! do not exit as we want to check there is not a second pair
!         pair=ii; lowa=jj-1
! ERROR: we have to loop mequad again to find jj! Or program smarter
         pair=ii
!         write(*,*)'3XQ loop to find the A/X quad index'
         notneeded: do jj=1,mqmqa_data%ncat
            if(temp(ii).lt.mqmqa_data%emquad(jj)) exit notneeded
!            if(temp(ii).gt.mqmqa_data%emquad(jj)) then
!               lowa=jj-1
!               exit notneeded
!            endif
         enddo notneeded
         lowa=jj-1
! lowa saves the quad index of the A/X quad for the AB/X quad
!         write(*,*)'3XQ the pair is quad ',ii,lowa
      else
         write(*,*)'3XQ convert_y2quads found two pair fractions in a parameter'
         gx%bmperr=4399; goto 1000
      endif
   enddo loop4
! if lowa=0 we have not found the AB/X quad
   if(lowa.eq.0) then
      write(*,*)'3XQ cannot find the AB/X quad',(temp(ii),ii=1,nq3),&
      ', among ',(mqmqa_data%emquad(ii),ii=1,mqmqa_data%ncat)
      stop
   endif
! set the pair as first quad in pquad; maybe change qorder
!   write(*,30)pair,lowa,(qorder(ii),ii=1,nq3)
30 format('3XQ we found the pair: 'i3,', lowa:',i3,', qorder:',15i3)
!   write(*,40)'3XQ pquad  before:',(pquad(ii),ii=1,nq3)
!   write(*,40)'3XQ qorder before',(qorder(ii),ii=1,nq3)
   if(pair.ne.1) then
! shift positions
      jj=pquad(1); kk=qorder(1)
      pquad(1)=pquad(pair); qorder(1)=qorder(pair)
      pquad(pair)=jj; qorder(pair)=kk
   endif
!   write(*,40)'3XQ pquad after:',(pquad(ii),ii=1,nq3)
!   write(*,40)'3XQ qorder after',(qorder(ii),ii=1,nq3)
40 format(a,4i3)
! it seems OK here ..................
! now pquad(1) is the pair AB/X. make pquad(2) to be A/X and pquad(3) as B/X
! Probably there is a smart way but I am just fed up with this
! All other constituents must be single cations: A/X, B/X or C/X
!   write(*,*)'3XQ value of nq3',nq3
   if(nq3.eq.3) then
! It should be sufficient that temp(2) < temp(3)
! But if there is a 4th quad one has to eliminate the quad without A and B
      if(pquad(2).gt.pquad(3)) then
         if(pquad(3).ne.lowa) then
            write(*,*)'3XQ problems finding A/X quad',lowa,pquad(2)
            jj=pquad(2); pquad(2)=pquad(3); pquad(3)=jj
         endif
      endif
!      write(*,*)'3XQ order of pquad:',(pquad(kk),kk=1,nq3)
   else
! lowa must be the A/X quad because AB/X must be after A/X
! the difference between quad AB/X and A/X must be related to the B/X
! pquad(1) is the index of AB/X quad, the A/X quad is lowa
!      write(*,20)'emquad again',(mqmqa_data%emquad(jj),jj=1,mqmqa_data%ncat)
      highb=pquad(1)-mqmqa_data%emquad(lowa)
!      write(*,*)'3XQ value of highb',pquad(1),mqmqa_data%emquad(lowa),highb
! the B/X quad should be highb indices in emquad higher than lowa
      nbx=mqmqa_data%emquad(lowa+highb)
!      write(*,*)'3XQ tables are turning:',pquad(3),pquad(4),nbx
      if(pquad(3).ne.nbx) then
         if(pquad(4).ne.nbx) then
            write(*,*)'3XQ circles are square'
            stop
         endif
         jj=pquad(4); pquad(4)=jj; pquad(3)=jj
      endif
!      write(*,*)'3XQ order of pquad:',(pquad(kk),kk=1,nq3)
   endif
! list everything
!   write(*,20)'emquad again',(mqmqa_data%emquad(jj),jj=1,mqmqa_data%ncat)
!   write(*,10)'final',(temp(ii),pquad(ii),ii=1,nq3)
!   write(*,666)(pquad(ii),ii=1,nq3)
666 format('3XQ fixed MQMQA parameter, quad is ',i3,', asymmetrical: ',10i3)
!   write(*,*)'3XQ hit return to handle next parameter'
!   read(*,*)
!
1000 continue
   return
 end subroutine convert_y2quadx

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine convert_y2quadx_old
!\begin{verbatim}
 subroutine convert_y2quadx_old(sem,nint,jord,pquad)
! This is to fix the constitution variables for an MQMQX excess parameter.
! It has one AB/X quad index and 2 A/X and B/X quads and possibly a C/X one
! convert y fraction indexed in sem and jord to quad indices in parquad
!
! I am really really fedup with this model
!
   implicit none
   integer sem,nint,jord(2,*),pquad(*)
!\end{verbatim}
   integer ii,jj,nq3,kk,pair,qorder(4),lowa,highb,temp(4),nbx
!
! input data from database
!   write(*,*)'3XQ *** fixing MQMQA parameter composition variables',&
!        size(mqmqa_data%emquad),size(mqmqa_data%con2quad)
!   write(*,5)sem,jord(2,1:nint)
5  format('3XQ fixing MQMQA parameter composition variables',10i3)
!   write(*,20)'emquad',(mqmqa_data%emquad(jj),jj=1,mqmqa_data%ncat)
!   write(*,20)'con2quad',(mqmqa_data%con2quad(jj),jj=1,mqmqa_data%nquad)
20 format('3XQ ',a,21i3)
! pquad(1) should be the pair quad AB/X among sem, jord(2,1..nint)
   temp(1)=mqmqa_data%con2quad(sem)
! pquad(2) should be the alphabetically first  in quad AB/X, i.e A/X
! pquad(3) should be the alphabetically second in quad AB/X, i.e B/X
! pquad(4) should be the 4th quand, not including A or B 
   temp(2)=mqmqa_data%con2quad(jord(2,1))
   temp(3)=mqmqa_data%con2quad(jord(2,2))
! the quad index related to temp(2) and temp(3) should be temp(1) ???
! check:
! emquad have the quad indices of all A/X quads, there are ncat of them.
! The index of a quad (i,j) where j>i is emquad(i)+j-i   
   write(*,*)'3XQ values temp: ',temp(1),temp(2),temp(3)
   if(temp(2).gt.temp(3)) then
      write(*,*)'3XQ parameter has wrong order of A/X and B/X quads'
      stop 76
   else
      ii=mqmqa_data%emquad(temp(2))+temp(3)-temp(2)
      write(*,*)'3XQ values of mixed quad index: ',temp(1),temp(2),temp(3),ii
      if(temp(1).ne.ii) then
         write(*,*)'3XQ problems with quad indices'
         stop 77
      endif
   endif
   nq3=3
   if(nint.eq.3) then
      nq3=4; temp(4)=mqmqa_data%con2quad(jord(2,3))
   endif
!   write(*,10)'first',(temp(ii),pquad(ii),ii=1,nq3)
10 format('3XQ ',a,' quads ',2i3,', first ',2i3,', second ',2i3,', maybe ',2i3)
! find the AB/X quad and the arrange the others
! all but one of the quads in temp(1..nq3) should be A/X quads
! and temp(2) should have the lowest index of the AB/x quad and temp(3)
! the highest.  Any temp(4) quad should not be A/X or B/X
! this code is horrible
!
   pair=0; lowa=0
   loop4: do ii=1,nq3
      pquad(ii)=temp(ii)
      qorder(ii)=ii
!
! qx is quad
! Calculate: quad(temp(1)*\xi(temp(3),temp(2)**ppow
!      write(*,*)'3XQ is this the pair?',pquad(ii),qorder(ii)
      loopax: do jj=1,mqmqa_data%ncat
!
! cycle loop4 if temp(ii) is an A/X quad
!
         if(temp(ii).eq.mqmqa_data%emquad(jj)) cycle loop4
      enddo loopax
! if we arrive here temp(ii) is a AB/X quad
      if(pair.eq.0) then
! do not exit as we want to check there is not a second pair
!         pair=ii; lowa=jj-1
! ERROR: we have to loop mequad again to find jj! Or program smarter
         pair=ii
!         write(*,*)'3XQ loop to find the A/X quad index'
         notneeded: do jj=1,mqmqa_data%ncat
            if(temp(ii).lt.mqmqa_data%emquad(jj)) exit notneeded
!            if(temp(ii).gt.mqmqa_data%emquad(jj)) then
!               lowa=jj-1
!               exit notneeded
!            endif
         enddo notneeded
         lowa=jj-1
! lowa saves the quad index of the A/X quad for the AB/X quad
!         write(*,*)'3XQ the pair is quad ',ii,lowa
      else
         write(*,*)'3XQ convert_y2quads found two pair fractions in a parameter'
         gx%bmperr=4399; goto 1000
      endif
   enddo loop4
! if lowa=0 we have not found the AB/X quad
   if(lowa.eq.0) then
      write(*,*)'3XQ cannot find the AB/X quad',(temp(ii),ii=1,nq3),&
      ', among ',(mqmqa_data%emquad(ii),ii=1,mqmqa_data%ncat)
      stop
   endif
! set the pair as first quad in pquad; maybe change qorder
!   write(*,30)pair,lowa,(qorder(ii),ii=1,nq3)
30 format('3XQ we found the pair: 'i3,', lowa:',i3,', qorder:',15i3)
!   write(*,40)'3XQ pquad  before:',(pquad(ii),ii=1,nq3)
!   write(*,40)'3XQ qorder before',(qorder(ii),ii=1,nq3)
   if(pair.ne.1) then
! shift positions
      jj=pquad(1); kk=qorder(1)
      pquad(1)=pquad(pair); qorder(1)=qorder(pair)
      pquad(pair)=jj; qorder(pair)=kk
   endif
!   write(*,40)'3XQ pquad after:',(pquad(ii),ii=1,nq3)
!   write(*,40)'3XQ qorder after',(qorder(ii),ii=1,nq3)
40 format(a,4i3)
! it seems OK here ..................
! now pquad(1) is the pair AB/X. make pquad(2) to be A/X and pquad(3) as B/X
! Probably there is a smart way but I am just fed up with this
! All other constituents must be single cations: A/X, B/X or C/X
!   write(*,*)'3XQ value of nq3',nq3
   if(nq3.eq.3) then
! It should be sufficient that temp(2) < temp(3)
! But if there is a 4th quad one has to eliminate the quad without A and B
      if(pquad(2).gt.pquad(3)) then
         if(pquad(3).ne.lowa) then
            write(*,*)'3XQ problems finding A/X quad',lowa,pquad(2)
            jj=pquad(2); pquad(2)=pquad(3); pquad(3)=jj
         endif
      endif
!      write(*,*)'3XQ order of pquad:',(pquad(kk),kk=1,nq3)
   else
! lowa must be the A/X quad because AB/X must be after A/X
! the difference between quad AB/X and A/X must be related to the B/X
! pquad(1) is the index of AB/X quad, the A/X quad is lowa
!      write(*,20)'emquad again',(mqmqa_data%emquad(jj),jj=1,mqmqa_data%ncat)
      highb=pquad(1)-mqmqa_data%emquad(lowa)
      write(*,*)'3XQ value of highb',pquad(1),mqmqa_data%emquad(lowa),highb
! the B/X quad should be highb indices in emquad higher than lowa
      nbx=mqmqa_data%emquad(lowa+highb)
!      write(*,*)'3XQ tables are turning:',pquad(3),pquad(4),nbx
      if(pquad(3).ne.nbx) then
         if(pquad(4).ne.nbx) then
            write(*,*)'3XQ circles are square'
            stop
         endif
         jj=pquad(4); pquad(4)=jj; pquad(3)=jj
      endif
!      write(*,*)'3XQ order of pquad:',(pquad(kk),kk=1,nq3)
   endif
! list everything
!   write(*,20)'emquad again',(mqmqa_data%emquad(jj),jj=1,mqmqa_data%ncat)
!   write(*,10)'final',(temp(ii),pquad(ii),ii=1,nq3)
!   write(*,666)(pquad(ii),ii=1,nq3)
666 format('3XQ fixed MQMQA parameter, quad is ',i3,', asymmetrical: ',10i3)
!   write(*,*)'3XQ hit return to handle next parameter'
!   read(*,*)
!
1000 continue
   return
 end subroutine convert_y2quadx_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine set_quadfractions(phres,verbose,yfr)
!\begin{verbatim}
 subroutine set_quadfractions(phres,verbose,yfra)
! copy values from yfr to xquad, y_ik etc using con2quad
! mqmqa_data%initaties phase variables for new mqmqa excess model
! the normal fractions, used by the config entropy, already set
   implicit none
   type(gtp_phase_varres), pointer :: phres
   type(gtp_mqmqa_var), pointer :: mqmqaf
   double precision yfra(*)
   logical verbose
!   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ia,iq
!   write(*,10)
10 format('3XQ in set_quadfractions, use con2quad for yfr to xquad'/&
        'then call calcasymvar to set \varkappa, \xi and Y_ik.',&
        'Latt som en platt')
   mqmqaf=>phres%mqmqaf
   if(.not.associated(mqmqaf)) then
      write(*,*)'3XQ there is no mqmqaf record for this phase'
      stop
   end if
   if(verbose) write(*,20)size(mqmqa_data%con2quad),&
        (mqmqa_data%con2quad(ia),ia=1,mqmqa_data%nquad)
20 format('3XQ mqmqaf%con2quad: ',i3,2x,20i3)
   do ia=1,mqmqa_data%nquad
      iq=mqmqa_data%con2quad(ia)
! I am not sure how to copy from yfr to mqmqaf%xquad
      mqmqaf%xquad(ia)=phres%yfr(iq)
      if(verbose) write(*,26)ia,phres%yfr(ia),iq,mqmqaf%xquad(ia)
26    format('3XQ the OC fraction: ',i3,1pe14.6,&
           ' is set to MQMQA quad: ',i3,1pe14.6)
   enddo
!  if(verbose) write(*,*)'3XQ calling calcasymvar for \varkappa_ij, \xi_ij etc.'
   call calcasymvar(phres)
!   if(verbose) write(*,*)'3XQ back from calcasymvar'
!
1000 continue
   return
 end subroutine set_quadfractions
 
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine init_excess_asymm
!\begin{verbatim}
 subroutine init_excess_asymm(lokph)
! subroutine init_excess_asymm(lokph,ic,ia)
! initaties phase variables for new mqmqa excess model
! called from gtp3B create_asymmetry
! number of independent quads, ic cations, ia anions (max 1)
   implicit none
! ic is number of cations, ia number of anions, there are also avalable globally
   integer ic,ia,lokph
!   type(gtp_phase_record), pointer :: phase
   type(gtp_ternary_asymmetry), pointer :: asym3rec
! there is a global mqmqa_data record to use!! <<<<<<<<<<<<<<<<<<,
!\end{verbatim}
   integer i,j,k,nseq,mm,apos,nbinsys,ntercat
!   integer i,j,k,nseq,mm,apos,lcat,lnan,nbinsys,ntercat
! how to create xquad mm when we need a pointer to gtp_phase_varres?
   type(gtp_equilibrium_data), pointer :: ceq
   type(gtp_phase_varres), pointer :: phres
   type(gtp_mqmqa_var), pointer :: mqf
   character*6 defasym
! Many properties are symmetric, for example xquad which has a single index
! and is indexed by ijkl(i,j,k,l) where ijkl(i,j,k,l)=ijkl(j,i,k,l)
! but other are unsymmetric such as varkappa and xi
!
!   write(*,*)'3QX In init_excess_asymm <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
!
   ceq=>firsteq
! I have forgotten how OC works.  When entering phases one can creat
! data structures in gtp_equilibrium_data (record pointer ceq)
! and these will be copied when new equilibrium records created
! (for example parallel calculations).  When a second gtp_equilibium_data
! has been created one is not allowed to change this data_structure
! the line below creates a pointer to the gtp_mqmqa_var data inside ceq
! maybe problem with the array here ...
   i=1
5  continue
      i=i+1
! this is very clumsy, but I have no better way
      phres=>ceq%phase_varres(i)
!      write(*,*)'loop: ',i,lokph,phres%phlink
      if(phres%phlink.ne.lokph) goto 5
!
!   write(*,*)'Found phase_varres!',i
   mqf=>phres%mqmqaf
! xquad is declared globally in OC BUT maybe better if
! xquad is declared within the gtp3_phase_varres mqmqaf record ??????/
! The mqf below is part of gtp_phase_varres 
! initiate with equal amount
! The xquads in a sequental array and used ONLY to calculate excess parameters 
! number of binary cation combination, in a binary 11-12-22
! in xquad the order is sequential in the cation order
!    1   2   3   4  ..  n   ! n+1 n+2 .. 2n-1 ! 2n  2n+1 .. ! ... ! n(n+1)/2
!    1/1 1/2 1/3 1/4    1/n ! 2/2 2/3 .. 2/n  ! 3/3 3/4  .. ! ... ! n/n   
! the function ijklx(i,j,k,l) calculates the sequential index
! we have to inititate con2quad below with the corresponding cation indices
! now we can create the xquad array and other things in mqf 
   if(.not.allocated(mqf%xquad)) then
!      write(*,*)'3XQ allocating xquad',mqmqa_data%nquad,mqmqa_data%nconst
      allocate(mqf%xquad(mqmqa_data%nquad))
      mqf%xquad=1.0d0/mqmqa_data%nquad
   endif
!
!   write(*,*)'3XQ *** Creation of cross indices for fractions and quads below'
!
!   ncat=ic
!   nan=ia
!   write(*,10)trim(phlista(lokph)%name),mqmqa_data%ncat,mqmqa_data%nan
10 format(/'3XQ Initiating mqmqa model for ',a,' with ',i3,' cations and ',&
        i2,' anion')
! initiate also values in mqmqa_data
!   write(*,*)'3XQ init_excess check:',mqmqa_data%ncon1,mqmqa_data%ncat,&
!        mqmqa_data%ncon2,mqmqa_data%nan,mqmqa_data%lcat
!   mqmqa_data%ncon1=ncat
!   mqmqa_data%ncon2=nan
!   mqmqa_data%lcat=ncat*(ncat+1)/2
! FNN/SNN ratio same for all pairs ...   in first work: qfnnsnn
!   allocate(mqmqa_data%etafs(ncat*nan))
!   mqmqa_data%etafs=2.4D0
! same as qfnnsnn
! the molefration xquad(1,2) is the same as xquad(2,1) and xquad
!   lcat=ncat*(ncat+1)/2
! lnan=1 if only one anion
!   write(*,*)'3xq value of lnan: ',mqmqa_data%lnan
!   lnan=nan*(nan+1)/2
! total number of quads, 
!>>>>>>> nquad, ncat, nan, lcat and lnan are global variables !!!!!!!!!!
! CHANGE TO USE VALUES IN MQMQA_DATA!!!
!   write(*,11)mqmqa_data%ncat,mqmqa_data%nan,mqmqa_data%nquad,&
!        mqmqa_data%ncon1,mqmqa_data%ncon2
11 format('3XQ mqmqa_data: ',10i4)
!   if(mqmqa_data%ncat.gt.1 .and. mqmqa%data%nan.gt.1) then
! cations 1 and 2 form quads 1/1 1/2 2/2 but xquad(2,1) same as xquad(1,2)
! 11, 12, 22 are separate quad fractions
!      nquad=ncat*(ncat+1)/2*nan*(nan+1)/2
!   elseif(nan.eq.1) then
! frequantly there will be a single anion
!      nquad=ncat*(ncat+1)/2
!   endif
!!
!-------------------------------------------
! now initate record with asymmetries
!   write(*,*)'Allocating asymmetries',mqmqa_data%ncat
   nseq=0
   if(mqmqa_data%nan.eq.1) then
      if(mqmqa_data%ncat.gt.1) then
         nbinsys=mqmqa_data%ncat*(mqmqa_data%ncat-1)/2
!...allocate ternary structure with asymmetry data
         if(mqmqa_data%ncat.gt.2) then
            ntercat=mqmqa_data%ncat*(mqmqa_data%ncat-1)*(mqmqa_data%ncat-2)/6
            allocate(tersys(ntercat))
! insert element indices
            mm=0
            do i=1,mqmqa_data%ncat-2
               do j=i+1,mqmqa_data%ncat-1
                  do k=j+1,mqmqa_data%ncat
! initiate all ternaries as symmetrical el(1) < el(2) < el(3)
                     mm=mm+1
                     if(mm.gt.ntercat) then
                        write(*,*)'wrong allocation of ntercat',mm,ntercat
                        stop
                     endif
                     tersys(mm)%seq=mm
                     tersys(mm)%el(1)=i
                     tersys(mm)%el(2)=j
                     tersys(mm)%el(3)=k
                     tersys(mm)%asymm='KKK'
                     tersys(mm)%isasym=0
                  enddo
               enddo
            enddo
!            write(*,17)mm
17          format('init_excess_ asymm allocated ternary structures ',i3)
            if(mm.ne.ntercat) then
               stop 'ternary allocation error'
            endif
!         else
!            write(*,*)'3XQ No ternary data structures needed'
         endif
      else
         write(*,*)'A liquid with a single cation and anion not implemented'
         stop
      endif
   else
      write(*,*)'Systems with multiple anions not implemented'
      stop
   endif
! varkappa and xi_ijis now part of allinone
!
!   if(mqmqa_data%ncat.eq.2) goto 80
!   write(*,67)mqmqa_data%ncat*(mqmqa_data%ncat-1)*mqmqa_data%nan/2
67 format('3XQ init_excess_asymm allocating asymmetrical compvar array: ',i5)
! we have to intitiate several variables in each compvar
!   allocate(compvar(ncat*(ncat-1)*nan/2))
!
   allocate(mqf%compvar(mqmqa_data%ncat*(mqmqa_data%ncat-1)/2*mqmqa_data%nan))
!   write(*,*)'3QX, initiating compvar for excess model variables',&
!        mqmqa_data%ncat,size(mqf%compvar)
   if(allocated(mqmqa_data%el2ancat)) then
!      write(*,69)
69    format('Heureca! el2ancat allocated')
!      write(*,70)size(mqmqa_data%el2ancat),mqmqa_data%ncat,mqmqa_data%el2ancat
70    format('3XQ el2ancat: ',2i3,5x,20i3)
   else 
      write(*,*)'3XQ line 2168: The array mqmqa_data%el2ancat not allocated!'
      write(*,*)'3XQ should have been done in correlate_const_and_quads'
      gx%bmperr=4399; goto 1000
   endif
!
   nseq=0
   mm=0
! it would have been better allocate compvar as this ...
   allocate(mqmqa_data%quad2compvar(mqmqa_data%ncat*(mqmqa_data%ncat+1)/2))
   dum1: do i=1,mqmqa_data%ncat
      dum2: do j=i,mqmqa_data%ncat
         nseq=nseq+1
         if(i.ne.j) then
            mm=mm+1
            mqmqa_data%quad2compvar(nseq)=mm
         else
            mqmqa_data%quad2compvar(nseq)=10000
         endif
      enddo dum2
   enddo dum1
!   write(*,71)mqmqa_data%quad2compvar
71 format('3XQ check quad2compvar',50i3)
!
   nseq=0
   first: do i=1,mqmqa_data%ncat-1
      second: do j=i+1,mqmqa_data%ncat
! initiallize allinone record, allocated as compvar array
         nseq=nseq+1
         mqf%compvar(nseq)%seq=nseq
! these indices are from 1 to n-1 ignoring anions
         mqf%compvar(nseq)%cat1=i
         mqf%compvar(nseq)%cat2=j
! these are the element indices in OC
         if(i.gt.mqmqa_data%xanionalpha) mqf%compvar(nseq)%elcat1=i+1
         if(j.gt.mqmqa_data%xanionalpha) mqf%compvar(nseq)%elcat2=j+1
! note it is negative of element alphabetical index
         mqf%compvar(nseq)%elan=-mqmqa_data%xanionalpha
         mqf%compvar(nseq)%anion=1
         mqf%compvar(nseq)%lastupdate=-1
! ivk_ij, jvi_ji, kvk_ijk, xi_ij etc allocated at each calculation
! NOTE vk_ij, xi_ij are single variables in each box, no need to allocate
         mqf%compvar(nseq)%vk_ij=0.0d0
         mqf%compvar(nseq)%vk_ji=0.0d0
         mqf%compvar(nseq)%xi_ij=0.0d0
         mqf%compvar(nseq)%xi_ji=0.0d0
! For identifying m used in eq.25 or 26 in Max paper for ternary excess
! in varkappa1 allocate arrays for which quad fractions vk and xi depend
! they can be different for each compvar
! %dvk_ij and %vdk_ji are single variables, arrays for derivatives
! %dvkx_ij and %vdkx_ji are type(zquad) ??, alternative arrays for derivatives
! allocated at first calculation
         allocate(mqf%compvar(nseq)%dxi_ij(mqmqa_data%nquad)) ! dxi_ij/dquad_k
         allocate(mqf%compvar(nseq)%dxi_ji(mqmqa_data%nquad)) ! dxi_ji/dquad_k
! The arrays dxi_ij are allocated here but xi_ij are single values in compvar
         mqf%compvar(nseq)%dxi_ij=0.0d0
         mqf%compvar(nseq)%dxi_ji=0.0d0
!         write(*,77)nseq,i,j
77       format(i4,2i5)
      enddo second
   enddo first
80 continue

! initiate newXupdate, there is a newupdate I do not know where it is declared
!   write(*,*)'3XQ newXupdate for varkappa and xi ',newXupdate
!   write(*,*)' *** Where is newupdate declared? ',newupdate
   newXupdate=0
! allocate quadz with zi_ijkl, for a single anion
! There are some data for Zv_ij/kl declared in mqmqa_data, use that!!
!   write(*,79)ncat*(ncat+1)/2
79 format('Allocating ',i3,' quadz array, for zv_ijkl data')
!   allocate(quadz(ncat*(ncat+1)/2))
! create crossreferences beween OC datastructure and MQMQX asymmetric
!   write(*,*)'3qx **** init_excess_asymm calls correlate_const_and_quads'
!
! THIS ROUTINE CALLS THIS ONE   call correlate_const_and_quads(lokph)
!   if(gx%bmperr.ne.0) goto 1000
!
!   write(*,90)mqmqa_data%ncat*mqmqa_data%nan
90 format('Allocating pair fraction array y_i/k: ',i4)
! y_ik varies with the current constitution
   allocate(mqf%y_ik(mqmqa_data%ncat*mqmqa_data%nan))
! with multiple anion derivatives add dimension nan also
! its content is set in varkappa1
!   write(*,*)'3XQ allocating mqf%dy_ik: ',ncat,nan,nquad, assume nan=1
! dy_ik is a structure information, independent of current constitution
   allocate(mqmqa_data%dy_ik(mqmqa_data%ncat,mqmqa_data%nquad))
   call pairfracs(.false.,mqf)
!
1000 continue
!
! REMOVE ncat from global data structure
  write(*,99)mqmqa_data%ncat,mqmqa_data%nquad,size(phres%yfr),size(mqf%compvar)
99 format(/'3QX **** leaving init_excess_asym : ',10i4/)
!
   return
 end subroutine init_excess_asymm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine pairfracs
!\begin{verbatim}
 subroutine pairfracs(list,mqf)
! calculate all pair fractions from a set of quad fractions
! pair fractions are y_v/x = \sum_ik/kl x_ij/kl*(delta_iv+delta_jv)/etafs
! if there is a single anion
   implicit none
   logical list
   type(gtp_mqmqa_var), pointer :: mqf
!\end{verbatim}
   integer i,j,v,dd,seq
!   double precision sum,nofs(ncat),yfs(ncat),sum1,sum2,zz
   double precision sum,sum1,sum2,zz
   double precision nofs(mqmqa_data%ncat)
   double precision yfs(mqmqa_data%ncat)
! how to move variables in the mqmqa_var record ??
! mqf is a pointer!!
!
   seq=0
   if(.not.allocated(mqf%xquad)) then
      write(*,*)'xquad not allocated'
      stop
   endif
!   write(*,*)'3XQ In pairfracts ncat and nan: ',ncat,nan
   if(list) write(*,6)
6  format(/'Calculating pair fractions'/&
        6x,'seq    i  j sum   sum2     1     2     3     4      5     6')
   sum1=0.0d0
   sum2=0.0d0
   do i=1,mqmqa_data%ncat
! loop will count each quad once including 11, 22 etc.
      do j=i,mqmqa_data%ncat
         if(mqmqa_data%nan.ne.1) then
            write(*,*)'Cannot calculate pair fractions with 2 or more anions'
            stop
         endif
         seq=seq+1
         zz=0.5d0*mqf%xquad(ijklx(i,j,1,1))
         if(seq.ne.ijklx(i,j,1,1)) then
! test for bugs ...
            write(*,*)'In pairfracs, ijklx and seq does not agree',seq
            stop
         endif
! if i=j they are added here
         nofs(i)=nofs(i)+zz
         nofs(j)=nofs(j)+zz
! y_ik(i) is the sum of all quads fractions with element i divided by /etafs
! dy_ik(i,z) is 0.5/etafs(i) for quad z
         yfs(i)=yfs(i)+zz/mqmqa_data%qfnnsnn(i)
         yfs(j)=yfs(j)+zz/mqmqa_data%qfnnsnn(j)
! These are constants, only calculate once, seq is the quad index
!        dy_ik(i,seq)=0.5d0/etafs(i)
!        dy_ik(j,seq)=0.5d0/etafs(j)
! ignore etafs ... but we must take stoichiometry Zv_ijkl into account!
         if(i.eq.j) then
            mqmqa_data%dy_ik(i,seq)=1.0d0
         else
            mqmqa_data%dy_ik(i,seq)=0.5d0
            mqmqa_data%dy_ik(j,seq)=0.5d0
         endif
!
         sum1=sum1+2*zz
         sum2=sum2+zz/mqmqa_data%qfnnsnn(i)+zz/mqmqa_data%qfnnsnn(j)
         if(list) then
            write(*,7)seq,i,j,sum1,sum2,nofs
7           format('y_ik: ',i3,2x,2i3,2F6.3,2x,(10F6.3/))
         endif
      enddo
   enddo
! a lot of trouble but it seems to work now ....
! These y_ik are for symmetrical systems ....unsure if this fnnsnn used ???
   do i=1,mqmqa_data%ncat
      mqf%y_ik(i)=yfs(i)/sum2
   enddo
!   write(*,*)'3XQ line 3179 Calculated y_i/k from quad fractions',mqf%y_ik(1)
   if(list) then
      write(*,10)'\etafs   ',mqmqa_data%qfnnsnn,sum2
      write(*,10)'y_i/k:   ',mqf%y_ik,sum1
10    format(a,10F7.4)
      do i=1,mqmqa_data%ncat
         write(*,12)i,(mqmqa_data%dy_ik(i,dd),dd=1,mqmqa_data%nquad)
12       format('dy_ik/dq: ',i2,12F6.3)
      enddo
   endif
 end subroutine pairfracs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable function ijklx(i,j,k,l)
!\begin{verbatim}
 integer function ijklx(i,j,k,l)
! i and j are cations, the order irrelevant
! k and l are anions, the order irrelevant
! always use the lowest value of i and j as first index below
! always use the lowest value of k and lj as first index below
   implicit none
   integer i,j,k,l
!\end{verbatim}
   integer iquad,kquad,a,b,x,y
!
   iquad=0
! Important order!!!
! Note i,j same as j,i and k,l same as l,k, lowest index always first
! Index order: 
!          1          2              ncat     ncat+1    ncat+2
!   (1,1,1,1), (1,2,1,1), ... (1,ncat,1,1), (2,2,1,1) (2,3,1,1) ... 
!   (2,ncat,1,1), (3,3,1,1), ... (3,ncat,1,1), (4,4,1,1) ... (ncat,ncat,1,1)
!   (1,1,1,2), (1,2,1,2), ... (1,ncat,1,2), ... (ncat,ncat,1,2), 
!   (1,1,2,2), (1,2,2,2), ... (ncat,ncat,2,2), (1,1,3,3), (1,2,3,3),  ... 
!   (ncat,ncat,nan,nan)
! indices (2,1,1,1) is same as (1,2,1,1) !
!------------------------------
! confusion where ncat,nan etc are stored
!   write(*,*)'Calling ijklx with: ',i,j,k,l
!   write(*,5)mqmqa_data%ncon1,mqmqa_data%ncon2,mqmqa_data%nconst,&
!        mqmqa_data%npair,mqmqa_data%lcat
!5  format('ijklx fixed values',2i4,2x,5i4)
! The cation index i,j order i<=j to find (i-1)*ncat-i*(i-1)/2+j
! The anion index  k,l order k<=l to find (k-1)*nan-k*(k-1)/2+l
! For each set of anion indices there are lcat=ncat*(ncat+1)/2 cation fractions
   if(i.le.0 .or. i.gt.mqmqa_data%ncon1 .or. &
        j.le.0 .or. j.gt.mqmqa_data%ncon1) goto 2000
   if(k.le.0 .or. k.gt.mqmqa_data%ncon2 .or. &
        l.le.0 .or. l.gt.mqmqa_data%ncon2) goto 2000
!
   if(l.lt.k) then
      kquad=(l-1)*mqmqa_data%ncon2-l*(l-1)/2+k-1
!      write(*,10)l,k,mqmqa_data%ncon2,kquad
   else
      kquad=(k-1)*mqmqa_data%ncon2-k*(k-1)/2+l-1
!      write(*,10)k,l,mqmqa_data%ncon2,kquad
   endif
10 format('Anion index in ijklx:  ',2i3,2i10)
!        
   if(j.lt.i) then
      iquad=(j-1)*mqmqa_data%ncon1-j*(j-1)/2+i
!      write(*,20)j,i,mqmqa_data%ncon1,kquad
   else
      iquad=(i-1)*mqmqa_data%ncon1-i*(i-1)/2+j
!      write(*,20)i,j,mqmqa_data%ncon1,kquad
   endif
20 format('Cation index in ijklx: ',2i3,2i10)
   iquad=kquad*mqmqa_data%lcat+iquad
!   write(*,30)iquad,kquad,lcat,i,j,k,l
30 format('Index in xquad: ',i5,5x,3i5,5x,2i5)
   if(iquad.gt.mqmqa_data%nconst) goto 1000
   ijklx=iquad
!   write(*,*)'Return from ijklx with:',iquad
   return
! error
1000 write(*,1010)i,j,k,l,mqmqa_data%ncon1,mqmqa_data%ncon2,mqmqa_data%lcat,&
          kquad,iquad
1010 format(' *** Indexing error in ijklx',4i4,2x,7i5,/'Stop!!!!')
   stop
!   read(*,*)iquad
!   ijklx=iquad
   !   return
2000 continue
   write(*,*)'Values outside limits',i,j,k,l
   stop
 end function ijklx

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine binsym
! calculates sequential index of a binary system
!\begin{verbatim}
 integer function binsym(i,j)
! SEPARATE FOR CATIONS AND ANION BINARIES, maybe merge with gtp_allinone?
! The binary systems form a symmetric matrix where (i,j) is the same as (j,i)
! and data for this system is stored as a linear array where where i > j always
! This function return the sequantial index for the binary (i,j)
! it is essentially the same as ijklx but only one set of indices
! The dimension of the binary cation matrix is the global variable ncat ...
!
! i/j    1   2   3   4   5   6   end of previous row  ncat*(ncat-1)/2 = 6*5/2
!  1     0   1   2   3   4   5    0       (ncat-j)*(ncat-j-1)/2 10  4*5/2 = 10
!  2     1   -   6   7   8   9    5  15 - (6-2)*(6-1)/2 = 15-4*5/2 = 5
!  3     2   6   -  10  11  12    9  15 - (6-3)*(6-2)/2 = 15-3*4/2 = 9
!  4     3   7  10   -  13  14   12  15 - (6-4)*(6-3)/2 = 15-2*3/2 = 12
!  5     4   8  11  13   -  15   14  15 - (6-5)*(6-4)/2 = 15-1     =14
!  6     5   9  12  14  15   -   note (6,6) is not a binary!!!
   implicit none
   integer i,j
!\end{verbatim}
!
   integer ix,iy
   if(i.le.0 .or. i.gt.mqmqa_data%ncat) goto 1100
   if(j.le.0 .or. j.gt.mqmqa_data%ncat) goto 1100
   ix=0
   if(j.lt.i) then
      if(j.gt.1) then
         ix=mqmqa_data%ncat*(mqmqa_data%ncat-1)/2 -&
              (mqmqa_data%ncat-j)*(mqmqa_data%ncat-j+1)/2
      endif
      iy=ix+i-j
   else
! j > i
      if(i.gt.1) then
         ix=mqmqa_data%ncat*(mqmqa_data%ncat-1)/2 - &
              (mqmqa_data%ncat-i+1)*(mqmqa_data%ncat-i)/2
      endif
      iy=ix+j-i
   endif
   binsym=iy
1000 continue
   return
1100 write(*,*)'Indexing error in binsym ',i,j,iy
   iy=-1
   goto 1000
 end function binsym

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine order3
!\begin{verbatim}
 subroutine order3(i,j,v,a,b,c)
! subroutine to rearrange i, j, v in increasing order in a, b, c
   implicit none
   integer i,j,k,a,b,c,v
!\end{verbatim}
! Return i, j, k ordered in a<b<c, do no change i, j, v
   if(i.lt.j) then
      if(j.lt.v) then
! i < j < v
         a=i; b=j; c=v                      ! i j v
      elseif(v.lt.j) then
         if(i.lt.v) then
! i < v < j
            a=i; b=v; c=j                   ! i v j
         elseif(i.gt.v) then
! v < i < j
            a=v; b=i; c=j                   ! v i j
         else
! i=v
            write(*,10)'1: i=v', i,j,v
10          format('order3 error, two indices same ',a,2x,3i4)
            goto 1100
         endif
      else
! j=v
         write(*,10)'2: j=v',i,j,v
         goto 1100
      endif
   elseif(j.lt.v) then
! here when i >= j and v > j thus j is smallest
      a=j
      if(i.lt.v) then            
         b=i; c=v
      elseif(v.lt.i) then
         b=v; c=i
      else
         write(*,10)'3: i=v',i,j,v
         goto 1100
      endif
   elseif(v.lt.j) then
! here when i>j and j>v
      a=v; b=j; c=i
   else !
! two or more numbers equal
      goto 1100
   endif
   return
!    
1100 continue
   write(*,*)' *** Error in call to order3: ',i,j,v
   a=-1
 end subroutine order3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine order3KKK
!\begin{verbatim}
 subroutine order3KKK(i,j,v,a,b,c,kkk)
! subroutine to rearrange i, j, v in increasing order in a, b, c
   implicit none
   integer i,j,k,a,b,c,v,jj
   character kkk*3,kopia*3,ch1*1
!\end{verbatim}
! Return i, j, v ordered in a<b<c, do no change i, j, v
! i,j,k are 1 2 or 3, kkk is rearranged to correspond to the re-arrangement
   kopia=kkk
!
!   write(*,11)i,j,v,kkk
11 format('3XQ entering order3KKK ',3i4,5x,a)
!
   if(i.lt.j) then
      if(j.lt.v) then
! i < j < v
         a=i; b=j; c=v                      ! i j v
      elseif(v.lt.j) then
         if(i.lt.v) then
! i < v < j
            a=i; b=v; c=j                   ! i v j
         elseif(i.gt.v) then
! v < i < j
            a=v; b=i; c=j                   ! v i j
         else
! i=v
            write(*,10)'1: i=v', i,j,v
10          format('order3 error, two indices same ',a,2x,3i4)
            goto 1100
         endif
      else
! j=v
         write(*,10)'2: j=v',i,j,v
         goto 1100
      endif
   elseif(j.lt.v) then
! here when i >= j and v > j thus j is smallest
      a=j
      if(i.lt.v) then            
         b=i; c=v
      elseif(v.lt.i) then
         b=v; c=i
      else
         write(*,10)'3: i=v',i,j,v
         goto 1100
      endif
   elseif(v.lt.j) then
! here when i>j and j>v
      a=v; b=j; c=i
   else 
! two or more numbers equal
      goto 1100
   endif
1000 continue
! rearrange kkk to the new order of cations.  
! KTK means the Toop element should be the second, TKK third and KKT first.
! programming this makes me sick  Just for a single Toop element
   fix: do jj=1,3
      ch1=kopia(jj:jj)
      if(ch1.ne.'T') cycle fix
      if(jj.eq.1) then
! Txx: the Toop element was originally third
         if(c.eq.v) then
! and still is, no change
            exit fix
         elseif(a.eq.v) then
! the first element is now the Toop element, change to xxT
            kkk='KKT'; exit fix
         else
! the Toop element must now be the second element
            kkk='KTK'; exit fix
         endif
      elseif(jj.eq.2) then
! xTx: the Toop element was the second            
         if(a.eq.j) then
! the first element is now the Toop element, change to xxT
            kkk='KKT'; exit fix
         elseif(b.eq.j) then
! no change
            exit fix
         else
! it must be the third element
            kkk='TKK'; exit fix
         endif
      else
! xxT: the Toop element was the first ... exit if it still is
         if(a.eq.i) exit fix
         if(a.eq.j) then
! it is now the second
            kkk='KTK'
         else
! or finally it is now the third
            kkk='TKK'
         endif
      endif
   enddo fix
!   
!   write(*,3)kkk,kopia,a,b,c
3  format('3XQ rearranged? "',a,'" original "',a,'"  ',3i3)
   return
!    
1100 continue
   write(*,*)' *** Error in call to order3: ',i,j,v
   a=-1
   goto 1000
 end subroutine order3KKK

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable function terind
!\begin{verbatim}
 integer function terind(i,j,v)
! integer function terind(i,j,v,ncat)
! find sequential index of ternary system i, j, k
! simplified version, SEPARATE FOR CATIONS AND ANION mixing
! The ternary systems form a symmetric matrix where (i,j,k) is same as (j,k,i)
! and data for the ternary is stored as a linear array where where i<j<k
! The terind function return the sequantial index for the ternary
!
! Example of the sequantial ordering of ternary system for 6 elements
!  i  j/v  3   4   5   6
!  1   2   1   2   3   4 
!  1   2   -   5   6   7
!  1   2   -   -   8   9
!  1   2   -   -   -  10 ---- end of first index 1, first 10 sequential values
!  2   3   -  11  12  13
!  2   3   -   -  14  15
!  2   3   -   -   -  16 ---- end of first index 2, first 16 sequential values
!  3   4   -   -  17  18
!  3   4   -   -   -  19 ---- end of first index 3, first 15 sequential values
!  4   5   -   -   -  20 (4,5,6) is the last ternary, 6*5*4/6=20
!
   implicit none
! ncat is the global variable for the number of cations ... suck
!   integer i,j,v,ncat
   integer i,j,v
!\end{verbatim}
   integer ix,iy,iz,a,b,c,bin,bp,cp
!   write(*,*)'Enter terind ',i,j,v,mqmqa_data%ncat
   if(i.le.0 .or. i.gt.mqmqa_data%ncat .or. &
        j.le.0 .or. j.gt.mqmqa_data%ncat .or. &
        v.le.0 .or. v.gt.mqmqa_data%ncat) goto 1100
!
   if(mqmqa_data%ncat.eq.3) then
      iz=1; goto 1000
   endif
! rearrange i, j k to indices a b c in increasing order
   call order3(i,j,v,a,b,c)
   if(a.lt.0) goto 1000
!
! the lowest index is a >=1, ix is number of skipped ternary systems
   ix=mqmqa_data%ncat*(mqmqa_data%ncat-1)*(mqmqa_data%ncat-2)/6 - &
        (mqmqa_data%ncat-a+1)*(mqmqa_data%ncat-a)*(mqmqa_data%ncat-a-1)/6
! we now have a binary matrix for i,v with dimension bin, indexed by (bp,cp)
   bin=mqmqa_data%ncat-a
   bp=b-a
   cp=c-a
   iy=bin*(bin-1)/2-(bin-bp+1)*(bin-bp)/2+cp-bp
   iz=ix+iy
!   write(*,10)a,b,c,mqmqa_data%ncat,ix,bin,bp,cp,iy,iz
10 format('terind: ',4i4,8i6)
1000 continue
   terind=iz
   return
1100 write(*,*)'Indexing error in terind ',i,j,v
   iz=-1
   goto 1000
 end function terind

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine asymm
!\begin{verbatim}
 integer function asymm(t,i,j,v)
! The ternary specified by t has 3 elements i-j-v.  v is redundant ... ?
! return 0 if neither element i nor j are asymmetric elements in this ternary
! return 1 if element i is an asymmetric element
! return 2 if element j is an asymmetric element
! return 3 if element both i and j are asymmetric elements
   implicit none
   integer t,i,j,v
!\end{verbatim}
   integer a,selectij
!    testasym=.true.
! the TDB format KKK or T3KK has been converted to element index
!   write(*,10)t,tersys(t)%isasym,i,j
10 format('3XQ In asymm: ternary ',i3,' asymmetry: ',3i3,' binary ',2i3,' OK')
   selectij=0
   do a=1,3
! The i-j is the binary.  
! if i is asymmetric but not j return 1
! if j is asymmetric but not i return 2
! if both i and j are asymmetric return 3
!       if(tersys(t)%isasym(a).eq.i) then
! according to Nathalie message 2025/10.20 I mixed up i and j, so I changed
      if(tersys(t)%isasym(a).eq.j) then
         selectij=selectij+2
      elseif(tersys(t)%isasym(a).eq.i) then
         selectij=selectij+1
      endif
20    format('in asymm: ',i1,i5,5x,2i3,' no more')
   enddo
100 continue
!   write(*,110)selectij,tersys(t)%asymm,tersys(t)%isasym,i,j
   asymm=selectij
110 format('3XQ exit asymm: ',i2,5x,a,2x,3i3,5x,2i3)
   return
 end function asymm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine calcasymvar
!\begin{verbatim}
 subroutine calcasymvar(phres)
! subroutine calcasymvar(mqmqavar)
! This must be called whenever the quad fractions has changed
! It updates varkappaij, xiij etc for the whole system
! and stores them in compvar(bin) datastructure
! Currently programmed ONLY for a single anion
   implicit none
   type(gtp_phase_varres), pointer :: phres
!   type(gtp_mqmqa_var), pointer :: mqmqavar
!   type(gtp_mqmqa_var), pointer :: mqmqavar
!\end{verbatim}
   integer i,j,ia,seq,k,l,m,ny,abrakadabra
   type(gtp_mqmqa_var), pointer :: mqf
! how to create xquad mm when we need a pointer to gtp_phase_varres?
!   type(gtp_equilibrium_data), pointer :: ceq
!   type(gtp_phase_varres), pointer :: phres
!   type(gtp_mqmqa_var), pointer :: mqf
   type(gtp_allinone), pointer :: box
! attempt to move mqmqa variables into the mqmqa_var record
!   phres=>ceq%phase_varres(lokvar)
!   mqf=>phres%mqmqaf
!
!   if(.not.allocated(mqmqavar%xquad)) then
!      write(*,*)'3XQ No xquad array allocated'
!      goto 1000
!   endif
!   if(.not.allocated(mqmqavar%compvar)) then
!      write(*,*)'3XQ No compvar array allocated'
!      goto 1000
!   endif
!
   ia=1
!   if(allocated(phres%mqmqaf%compvar)) then
!      write(*,*)'3XQ in calcasym: compvar: ',size(phres%mqmqaf%compvar)
!   else
!      write(*,*)'3XQ in calasym: phres%mqmqaf%compvar not allocated'
!   endif
! the separate array of binaries redundant?
! when a change of ternary asymmetries is made the newXupdate is incremented
   seq=0
   do i=1,mqmqa_data%ncat-1
      do j=i+1,mqmqa_data%ncat
! seq specifies a binary set of elements
! results are stored in compvar(seq) for use in Gibbs energy calculations
         seq=seq+1
!         write(*,*)'Calling varkappa1 ',i,j,seq
!         call varkappa1(mqmqavar%compvar(seq))
!         call varkappa1(seq,mqmqavar)
!         call varkappa1(seq,mqf)
!         write(*,*)'3XQ calcasymvar call varkappa1'
         call varkappa1(seq,phres)
      enddo
   enddo
! inside varkappa1 one adds quads to vk_ij and vk_ji and 
! if one has ijklx(vz1,vz1,ia,ia) in vk_ij and ijklx(vz2,vz2,ia,ia) in vk_ij
! then the %kvk_ij needs an additional ijkl(vz1,vz2,ia,ia)
! Check that here .... (this is due to bad initial programming)
!   write(*,790)1
790 format('3XQ **** DOUBLE CHECK KVK_IJK',i3)
!   mqf=>phres%mqmqaf
!   box%lastupdate=-1
!   write(*,*)'3XQ box%lastupdate: ',box%lastupdate
!   write(*,790)
!   if(box%lastupdate.ne.newXupdate) then
!      box%lastupdate=newXupdate
!      write(*,1001)box%seq,box%lastupdate,newXupdate
!1001  format('3XQ allinone record ',i3,' updated to new asymmetries ',i5)
!   else
!      write(*,*)'3XQ line 3707: mixed asymmetries added'
!   endif
!
!   write(*,*)'3XQ code below skipped as moved to varkappa1'
   goto 1000
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Code below moved to varkappa1 ... but some numerical problems persist ...
   write(*,790)2
   write(*,791)mqmqa_data%emquad
791 format('3XQ Attempt to add mixed quads, em2quad: ',25i3)
   do i=1,size(mqf%compvar)
!  if in vk_ij one has added (vz1,vz1,ia,ia)
!  and in vk_ji added        (vz2,vz2,ia,ia)
! one must add (vz1,vz2,ia,ia) to the kvk_ij (now done in calling routine)
! BUT this quad may already be present  !!!!!!!!!!!
      box=>mqf%compvar(i)
      write(*,792)box%seq,box%lastupdate,newXupdate
792   format('3XQ newXupdate: ',i3,2i5)
!      write(*,800)i,box%cat1,box%cat2
!      write(*,805)'ivk_ij  ',box%ivk_ij
!      write(*,805)'jvk_ij  ',box%jvk_ji
!      write(*,805)'kvk_ijk ',box%kvk_ijk
      do j=2,size(box%ivk_ij)
         do k=1,size(mqmqa_data%emquad)
            if(box%ivk_ij(j).eq.mqmqa_data%emquad(k)) then
! we have an endmember quad in ivk_ij (in addition to the first)
! Check if we have another endmember quad in jvk_ji
               do l=1,size(box%jvk_ji)
                  neverending: do m=1,size(mqmqa_data%emquad)
                     if(box%jvk_ji(l).eq.mqmqa_data%emquad(m)) then
                        if(k.ne.m) then
! we have 2 different endmember quads in ivk_ij and jvk_ji, 
! if the mixed quad is not alreay present add it
                           ny=ijklx(k,m,ia,ia)
                           do abrakadabra=1,size(box%kvk_ijk)
! check if this quad not already in box_kvk_ijk

                           enddo
! add this quad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           box%kvk_ijk=[box%kvk_ijk, ijklx(k,m,ia,ia)]
                           write(*,806)i,k,m,ijklx(k,m,ia,ia)
                           write(*,805)'kvk_ijk ',box%kvk_ijk
                        endif
                     endif
                  enddo neverending
               enddo
            endif
         enddo
      enddo
! a quad representing a vz,vz,ia,ia quad is part of emquad
      box%lastupdate=newXupdate
   enddo
800 format('3XQ compvar: ',i3,2x,2i3)
805 format(a,20i3)
806 format('3XQ adding mixed quad to kvk_ijk',i3,2x,2i3,2x,i3)
!
1000 continue
   return
 end subroutine calcasymvar

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine varkappa1
!\begin{verbatim}
! subroutine varkappa1(seq,mqf)
 subroutine varkappa1(seq,phres)
! box is a record of the type(gtp_allinone)
! this routine may initiate, calculate and store varkappaij, varkappaji, xiij
! and xiji for symmetric and asymmetric systems with Kohler/Toop
! It is programmed for a single anion and just for the MQMQX phase!
!
! It will inintiate all data in box if box%lastupdate ne newXupdate
!
! should it be a pointer?  Does it matter?
   implicit none
   integer seq
   type(gtp_phase_varres), pointer :: phres
!   type(gtp_mqmqa_var), pointer :: 
!\end{verbatim}
!
! replaced original i and j by icat and jcat below!!    integer i,j,ia,bin
!
! these are quad indices of i,i, i,j abd j,j
   integer mii,mij,mjj,ia
   type(gtp_allinone), pointer :: box
   type(gtp_mqmqa_var), pointer :: mqf
!
! ia represent the single anion
! varkappaij and varkappaji are the 2 composition variables to be multiplied
! with a binary i-j parameter or ternary.  
! It is modified depending on the types of
! extrapolations for each ternary it is involved: Kohler or Toop.
! initially varkappaij=x_ii and varkappaji=x_jj and sum=x_ij 
!           for the Kohler (Muggianu not implemented)
! and nugamma is set to zero
! - if element i is Toop in j-i-v the x_iv is added to nugamma
! - if element j is Toop in i-j-v the x_jv is added to nugamma
! for those involved in asymmetric ternaries the divison must include \nu\gamma
! At the end both varkappaij and varkappaji are divided by sum where
!                                   sum+varkappaij+varkappaji+\nu\gamma
! Maybe check if nugamma already included in sum ...
!
   integer i,ii,vz,v,w,vv,ternary,ll,lasthope,di,icat,jcat,nnn
   double precision varkappaij,varkappaji,sum,initialij,initialji,nugamma
   double precision xi_ij,xi_ji,sum1,sum2
   logical asymmetric
! added nov 3/2025.  See this date below
! in mixnugamma all vz that have asymmetric ternary with icat or jcat are saved
! because their mixed quad fractions should be added to kvk_ijk
   integer, dimension(:), allocatable :: mixnugamma
   integer selectij,qz1,qz2
! mixed update
   integer j,k,l,m,ny,abrakadabra
! If a binary i-j is part of 2 or more asymmetric ternaries i-j-\nu, i-j-\gamma
! the quad fraction x_\nu\gamma should be added to kvk_ijk (the denomonator)
! of kvk_ijk
! saving separate asymmetrical cations for a binary
!   integer, dimension(:), allocatable :: savevz
! debug output
   integer nn1,nn2,nn3,nn4,nn5,nn6,nn7
   logical nysym
! local variables used for updating quad indices for iasymm, jasymm, etc
!    integer, dimension(:), allocatable :: vk_ij,vk_ji,vk_ijk,xi_ij,xi_ji
! all asymmetric quad indices needed are stored in each separate gtp_allinone
!    integer nvk_ij,nvk_ji,nvk_ijk,nxi_ij,nxi_ji
!
! how to create xquad mm when we need a pointer to gtp_phase_varres?
!   type(gtp_equilibrium_data), pointer :: ceq
!   type(gtp_phase_varres), pointer :: phres
!   type(gtp_mqmqa_var), pointer :: mqf
! attempt to move mqmqa variables into the mqmqa_var record
!   ceq=>firsteq
!
! I have forgotten to set y_ik !!!!
! maybe problem with the array here ...
   mqf=>phres%mqmqaf
!   write(*,*)'3XQ line 3767 updating vk_ij, xi_ij and y_ik with new quad fracs'
!   write(*,10)'3XQ old',(mqf%y_ik(v),v=1,mqmqa_data%ncat)
10 format(a,15(f8.5))
   do v=1,mqmqa_data%ncat
      mqf%y_ik(v)=0.0d0
!      write(*,20)'3XQ dy_ik',(mqmqa_data%dy_ik(v,w),w=1,mqmqa_data%nquad)
20    format(a,(20F5.2))
      do w=1,mqmqa_data%nquad
         mqf%y_ik(v)=mqf%y_ik(v)+mqmqa_data%dy_ik(v,w)*mqf%xquad(w)
      enddo
   enddo
!   write(*,10)'3XQ line 3731 y_ik:',(mqf%y_ik(v),v=1,mqmqa_data%ncat)
!
!   write(*,*)'3XQ in varkappa1',seq
!   write(*,*)'3XQ in varkappa1',mqf%nquad
   if(.not.allocated(mqf%compvar)) then
      write(*,*)'3XQ line 3076 in varkappa: compvar not allocated, problems'
      gx%bmperr=4399; goto 1000
!   else
!      write(*,*)'3XQ varkappa allocated OK'
   endif
!   phres=>ceq%phase_varres(1)
!   mqf=>phres%mqmqaf
!   box=>mqf%compvar(seq)
   box=>mqf%compvar(seq)
! icat and jcat represent cations ... duplicated here (and many other places)
   icat=box%cat1
   jcat=box%cat2
   ia=box%anion
!   write(*,*)'In varkappa1: ',icat,jcat,ia
! the xquad values i,j and j,i are the same but for varkappa they are different
   if(icat.gt.jcat) then
      write(*,3)icat,jcat
3     format(/'In varkappa1: wrong order of elements ',2i4)
      stop
   endif
   mii=ijklx(icat,icat,ia,ia)
   mij=ijklx(icat,jcat,ia,ia)
   mjj=ijklx(jcat,jcat,ia,ia)
! how to deallocate box%asymm_nu and box%asymm_gamma?
!   deallocate(box%asymm_nu)
!   deallocate(box%asymm_gamma)
!
   nysym=.false.
   if(newXupdate.gt.box%lastupdate) then
      write(*,5)box%seq,box%lastupdate,newXupdate
5     format(/'Updating allinone record ',i5,' from ',i5,' to ',i5)
! *** this code part needed only once when all asymmetries are defined
! the arrays below are allocated, the initial 0 is overwritten if used
! This makes use of the new Fortran 2003 facility using [ ]
! Setting an allocatable array to single value means previous values deleted
!      box%ivk_ij=[0]; box%jvk_ji=[0]; box%kvk_ijk=[0]
!
! vk derivatives are quad indices, also denominator (same vk_ij and vk_ji)
! the statements below allocate and assign initial quad index
      box%ivk_ij=[mii]; box%jvk_ji=[mjj]; box%kvk_ijk=[mij]
! to handle derivatives the denominator is summed separately
      box%all_ijk=[mii, mjj, mij]
! xi are the Y_i/k fractions, for derivatives save quad indices in dxi_ij
!
      do di=1,mqmqa_data%nquad
! the derivatives of xi_ relative to quad index di
! The derivatives involves several quads, given by dy_ik(icat)
         box%dxi_ij(di)=mqmqa_data%dy_ik(icat,di)
         box%dxi_ji(di)=mqmqa_data%dy_ik(jcat,di)
      enddo
! calculate xi_  ... why??
      box%xi_ij=0.0d0; box%xi_ji=0.0d0
      do di=1,mqmqa_data%nquad
         box%xi_ij=box%xi_ij+box%dxi_ij(di)*mqf%xquad(di)
         box%xi_ji=box%xi_ji+box%dxi_ji(di)*mqf%xquad(di)
      enddo
! *** end of symmetric initialization of vk_ij, vk_ji, xi_ij and xi_ji
!
!  if in vk_ij one has added (vz1,vz1,ia,ia)
!  and in vk_ji added        (vz2,vz2,ia,ia)
! one must add (vz1,vz2,ia,ia) to the kvk_ij (now done in calling routine)
!         allocate(savevz(mqmqa_data%ncat,mqmqa_data%ncat))
! Now take care of asymmetries and update for later use
! Asymmetric vk and xi are updated in the vz loop AND at the end of the loop
! Nathalie correction initialized <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! The term above the division line for varkappa is same but in the denominator
! (below the division line) one quad fraction is added.  
! Either x(vz,i) or x(vz,j) where i or j is the element missing in Bosses
! expression
!       nugamma=0.0d0
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! this is not a debug variable, just for test output!!
      asymmetric=.false.
!
!      write(*,*)'3QX in varkappa1'
      vzloop: do vz=1,mqmqa_data%ncat
! loop for all ternary systems to find with asymmetric i-j-vz and j-i-vz
!         write(*,403)'3XQ in vzloop 1: ',vz,icat,jcat
403      format(a,i3,2x,2i3,2x,i3,2x,5i3)
         if(vz.eq.icat .or. vz.eq.jcat) cycle vzloop
! check the ternary icat-jcat-vz exists
!         ternary=terind(icat,jcat,vz,ncat)
         ternary=terind(icat,jcat,vz)
!         write(*,403)'3XQ in vzloop 2: ',vz,icat,jcat,ternary
         if(ternary.le.0) goto 1100
!
         selectij=asymm(ternary,icat,jcat,vz)
!         write(*,403)'3XQ in vzloop 3: ',vz,icat,jcat,ternary,selectij
! asymm returns 1 if icat is an asymmetric element in icat-jcat-vz  (gamma)
! asymm returns 2 if jcat is an asymmetric element in icat-jcat-vz  (nu)
! asymm returns 3 if both icat and jcat are asymmetric in icat-jcat-vz
! to be considered:  asymmetric i-j-nu and i-j-gamma requires x_\nu\gamma
!                    in the denominator.  For this the savevz is used
!
! ********* selectij=0 means no asymmetry in this ternary***************
         if(selectij.eq.0) cycle vzloop
         if(.not.asymmetric) then
! the asymmetric logical is to just for debug output of initial varkappa values
            asymmetric=.true.
         endif
!
!******************** asymmetric ternary *****************************
!         write(*,420)selectij,icat,jcat,vz
420      format('3XQ Asymmetric ternary ',i2,3i3)
         asymmetry: select case(selectij)
!
         case default
            write(*,*)'Illegal asymmetry ',selectij
            stop
!-------------------------------------------------------------------
         case(1) ! *************************************************
! icat is asymmetric, save in ivk_ij 
! an elegant Fortran assignment of an additional items in an allocatable
            box%jvk_ji=[box%jvk_ji, ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! All quad fractions added to ivk_ij must also be added to denominator
            box%kvk_ijk=[box%kvk_ijk, ijklx(icat,vz,ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(jcat,vz,ia,ia), &
                 ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! The asymmetric xi is depend on y_ik update dxi_ij and dxi_ji
            do nnn=1,mqmqa_data%nquad
!                box%dxi_ij(nnn)=box%dxi_ij(nnn)+dy_ik(icat,nnn)
               box%dxi_ji(nnn)=box%dxi_ji(nnn)+mqmqa_data%dy_ik(vz,nnn)
            enddo
!
!            exit asymmetry
!---------------------------------------------------------------------
         case(2) ! ***************************************************
! jcat is asymmetric, same as for icat just change icat to jcat!!!!
! and save in jvk_ji ...
            box%ivk_ij=[box%ivk_ij, ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! why not            box%kvk_ijk=[box%kvk_ijk, ijklx(jcat,vz,ia,ia)]
!            box%kvk_ijk=[box%kvk_ijk, ijklx(jcat,vz,ia,ia)]
! Nath noted missing  ijklx(vz1,vz2,ia,ia) if icat and jcat are asymmetrical
            box%kvk_ijk=[box%kvk_ijk, ijklx(jcat,vz,ia,ia)]  ! missing
            box%all_ijk=[box%all_ijk, ijklx(icat,vz,ia,ia), &
                 ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! The asymmetric xi is depend on y_ik update dxi_ij and dxi_ji
            do nnn=1,mqmqa_data%nquad
!                box%dxi_ij(nnn)=box%dxi_ij(nnn)+dy_ik(jcat,nnn)
               box%dxi_ij(nnn)=box%dxi_ij(nnn)+mqmqa_data%dy_ik(vz,nnn)
            enddo
!
!            exit asymmetry
!---------------------------------------------------------------------
         case(3) ! **************************************************
! Both icat and jcat are asymmetric
            box%ivk_ij=[box%ivk_ij, ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
            box%jvk_ji=[box%jvk_ji, ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! This is complicated, do not add ijklx(icat,vz,ia,ia), ijklx(jcat,vz,ia,ia)
! and only once ijkl(vz,vz,ia,ia) .....
! maybe not at all ?????????????
!            box%kvk_ijk=[box%kvk_ijk, ijklx(vz,vz,ia,ia)]
!            box%kvk_ijk=[box%kvk_ijk, ijklx(icat,vz,ia,ia), &
!                 ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! BUT x_(vz,vz,ia,ia) appears twice in the denominator ....(and twice on top)
            box%all_ijk=[box%all_ijk, ijklx(icat,vz,ia,ia), &
                 ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! The asymmetric xi is depend on y_ik update dxi_ij and dxu_ji
            do nnn=1,mqmqa_data%nquad
               box%dxi_ij(nnn)=box%dxi_ij(nnn)+mqmqa_data%dy_ik(icat,nnn)
               box%dxi_ji(nnn)=box%dxi_ji(nnn)+mqmqa_data%dy_ik(jcat,nnn)
            enddo
!
!            exit asymmetry
         end select asymmetry
! if selectij is 1 the asymmetry fixed for this ternary, if 2 or 3 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! code handling kvk_ijk terms due to extra x_ii and x_jj in ivk_ij and jvk_ji
! copied from end of calcasymvar to avoid it is repeted at all calculations
! skip first ivk_ij
         addkvkterm: do j=2,size(box%ivk_ij)
            do k=1,size(mqmqa_data%emquad)
               if(box%ivk_ij(j).eq.mqmqa_data%emquad(k)) then
! we have an endmember quad in ivk_ij (in addition to the first)
! Check if we have another endmember quad in jvk_ji, skip first jvk_ji
!                  do l=1,size(box%jvk_ji)
                  do l=2,size(box%jvk_ji)
                     neverending: do m=1,size(mqmqa_data%emquad)
                        if(box%jvk_ji(l).eq.mqmqa_data%emquad(m)) then
                           if(k.ne.m) then
! we have 2 different endmember quads in ivk_ij and jvk_ji, 
! if the mixed quad is not alreay present add it
                              ny=ijklx(k,m,ia,ia)
                              do abrakadabra=1,size(box%kvk_ijk)
! check if this quad not already in box_kvk_ijk
                                 write(*,*)'3XQ check duplicate line 4041 !!'
                              enddo
! add this quad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              box%kvk_ijk=[box%kvk_ijk, ijklx(k,m,ia,ia)]
                              write(*,806)i,k,m,ijklx(k,m,ia,ia)
                              write(*,805)'kvk_ijk ',box%kvk_ijk
                           endif
                        endif
                     enddo neverending
                  enddo
               endif
            enddo
         enddo addkvkterm
805 format(a,20i3)
806      format('3XQ adding mixed quad to kvk_ijk',i3,2x,2i3,2x,i3)
! end copied code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo vzloop
! the loop above should only be done whenever the asymmetry changes
      write(*,*)'3XQ line 4063, asymmetry updated for: ',box%cat1,box%cat2
!--------------------------------------------------------------------
! end of asymmetry detection loop
!--------------------------------------------------------------------
!
   endif ! end of update loop

!--------------------------------------------------------------------
! Below arrays box%ivk_ij, box%jvk_ji, box%dxi_ij are used to
! calculate \varkappa and \xi and their derivatives
!--------------------------------------------------------------------
!
! Now use the structures ivk_ij, jvk_ji, kvk_ijk and dxi_ij, dxi_ji
!   write(*,*)'3QX in varkappa1 line 3900',allocated(box%ivk_ij),&
!        allocated(box%dvk_ij)
   varkappaij=0.0d0; varkappaji=0.0d0; sum=0.0d0; nugamma=0.0d0
   do ii=1,size(box%ivk_ij)
      varkappaij=varkappaij+mqf%xquad(box%ivk_ij(ii))
!       write(*,697)'ivk_ij',ii,box%ivk_ij(ii),varkappaij,xquad(box%ivk_ij(ii))
697   format('Summing ',a,': ',2i3,2(1pe14.6))
   enddo
600 format('Partial sum: ',i3,a,1pe12.4,' quad: ',5i3)
   do ii=1,size(box%jvk_ji)
      varkappaji=varkappaji+mqf%xquad(box%jvk_ji(ii))
!       write(*,697)'jvk_ji',ii,box%jvk_ji(ii),varkappaji,xquad(box%jvk_ji(ii))
   enddo
   do ii=1,size(box%kvk_ijk)
      sum=sum+mqf%xquad(box%kvk_ijk(ii))
!       write(*,697)'sum',ii,box%kvk_ijk(ii),sum,xquad(box%kvk_ijk(ii))
   enddo
! all quad indices
!    write(*,696)' all_ijk: ',box%all_ijk
696 format('Quad indices in',a,': ',20i4)
   sum=sum+varkappaij+varkappaji+nugamma
!    write(*,601)sum,nugamma
601 format('Total value      Denominator: ',1pe12.4,' nugamma: ',1pe12.4)
! save normalized values and save also sum for use with derivatives
! at initiation sum=0.0, fix that
   if(sum.eq.0.0d0) sum=1.0d0
   box%vk_ij=varkappaij/sum
   box%vk_ji=varkappaji/sum
! the denominantor needed for derivatives
   box%denominator=sum
!    write(*,605)' vk_ij and vk_ji: ',box%vk_ij,box%vk_ji
605 format(' ** Normalized values of ',a,2(1pe12.4))
! and the derivatives ....
!
! Calculation of xi_ij using dxi
   sum1=0.0d0; sum2=0.0d0
   do di=1,mqmqa_data%nquad
      sum1=sum1+box%dxi_ij(di)*mqf%xquad(di)
      sum2=sum2+box%dxi_ji(di)*mqf%xquad(di)
   enddo
   box%xi_ij=sum1
   box%xi_ji=sum2
!
! debug output, ivk_ij, jvk_ji, kvk_ijk, dxi_ij, dxi_ji ---------------------
!    
   if(mqmqdebug .or. mqmqxcess) then
      nn1=size(box%ivk_ij); nn2=size(box%jvk_ji); nn3=size(box%kvk_ijk)
      nn4=mqmqa_data%nquad; nn5=mqmqa_data%nquad;
      if(allocated(box%asymm_nu)) then
         nn6=size(box%asymm_nu)
      else
         nn6=0
      endif
      if(allocated(box%asymm_gamma)) then
         nn7=size(box%asymm_gamma)
      else
         nn7=0
      endif
      write(*,700)2,nn1,nn2,nn3,nn4,nn5,nn6,nn7,nugamma
700   format('3XQ Sizes: ',i1,': ',7i3,1pe12.4)
      write(*,710)'ivk_ij  ',(box%ivk_ij(i),i=1,nn1)
      write(*,710)'jvk_ji  ',(box%jvk_ji(i),i=1,nn2)
      write(*,710)'kvk_ijk ',(box%kvk_ijk(i),i=1,nn3)
      write(*,709)'dxi_ij  ',(box%dxi_ij(i),i=1,nn4)
      write(*,709)'dxi_ji  ',(box%dxi_ji(i),i=1,nn5)
      if(nn6.gt.0) write(*,708)'nu      ',(box%asymm_nu(i),i=1,nn6)
      if(nn7.gt.0) write(*,708)'gamma   ',(box%asymm_gamma(i),i=1,nn7)
709   format('Factors ',a,': ',10f6.3)
708   format('Ternary quad asymmetry ',a,': ',5i4)
710   format('Quad in ',a,': ',5i4)
!
      write(*,607)3,box%vk_ij,box%vk_ji
607   format('Current values of vk_ij, vk_ji ',i2,2x,2(1pe15.5))
   endif
! end debug output ----------------------------------------------------------
! The asymmetric information collected as saved as quad index in local
! ivk_ij, jvk_ji, kvk_ijk for the \varkappa variables
! dxi_ij and dxi_ji for the \xi variables
!---------also for xi 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! this code use the updated data structure to calculate quickly
! This should be called by set constitution!!
!   write(*,*)'Calling dexcess_dq to allocate and set %dvk_ij etc'
   call dexcess_dq(box)
!   write(*,800)allocated(box%dvk_ij)
800 format(' *** Back from dexcess_dq to allocate %dvk_ij etc',l2)
   goto 900
!
500 continue    
!!!!!!!!!! here we use the asymmetry saved in box%asym1 and %asym2
!------------------------------------------------
! Here we calculate the derivatives using %asym1 and %asym2 ???
! ??
900 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   if(box%lastupdate.ne.newXupdate) then
      box%lastupdate=newXupdate
!      write(*,1001)box%seq,box%lastupdate
1001  format('3XQ allinone record ',i3,' updated to new asymmetries ',i5)
   endif
!
1000 continue
   return
1100 continue
   write(*,1105)icat,jcat,v
1105 format('Error return from tersym for elements: ',3i4)
   goto 1000
 end subroutine varkappa1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine dexcess_dq(box)
!\begin{verbatim}
 subroutine dexcess_dq(box)
! calculate the partial derivatives of a \varkappa or \xi variable
   implicit none
   type(gtp_allinone) :: box
! in ivk_ij, jvk_ji etc specify the indices of quad fractions involved for vk_ij
! A derivative wrt to a quad fractions included means it is 1, otherwise 0
! vk_ij has a numerator and a denominator, both are sums of quad fractions
! dvk_ij/dq_k is the result of both
! use type(zquad) for storing derivatives of ivk_ij, jvk_ji .. ???
! if_ij, if_ji, if_ijk are 0/1 depeding on the quad indices in each term
!    integer k,v,dtij,dtji,dtdiv
!\end{verbatim}
   integer k,v,d_ij,d_ji,d_ijk
   double precision numerator, denominator
   if(.not.allocated(box%dvk_ij)) then
! first time inititate arrays
      allocate(box%dvk_ij(mqmqa_data%nquad))
      allocate(box%dvk_ji(mqmqa_data%nquad))
   endif
! 
! calculate the derivatives of all vk_ij, vk_ji with respect to quads
! The quad indices are stored in ivk_ij, jvk_ji and kvk_ijk
!
!          \sum x_i    numerator           ivk_ij
! f=vk_i = --------- = ----------   = -------------------   \delta_mv=1 if m=k
!          \sum x_k    denominator    ivkij+jvkji+kvk_ijk
!
!           denominator*\delta_iv - numerator*\delta_ijkv
! df/dx_v = ---------------------------------------------   \delta_mv=1 if m=k
!                        denominator**2
! note value of numerator stored in vk_ij etc is already divided by denominator,
!                  \delta_iv       (numerator/denominator)*\delta_ijkv
! thus   df/dx_v = ------------  - -----------------------------------
!                   denominator             denominator
!
! many df/dx_v are zero ... trying to be smart? save only non-zero df/dx_v
!----------------------------------------------------------
! the arrays ivk_ij have only indices for the quads q they depend on
! vk_ij is the sum of those quads.  Many dvk_ij should be zero
! the denominator always depend on the same fractions as the numerator
   box%dvk_ij=0.0d0
!   write(*,10)box%seq,box%all_ijk
10 format('In dexcess_dq: allinone ',i3,' depend on quads: ',2x,20i3)
   kloop: do k=1,mqmqa_data%nquad
! we have to check all_ijk if vk depend on quad k
!       dvk_ij(k)=(if_ij/denominator_ij - if_ijk*numerator_ij)/denominator_ij
!       denominator_ijk and numerator_ij are sum of quad fractions
!
      d_ijk=0; d_ij=0; d_ji=0
!       write(*,15)box%all_ijk
15    format('kvk%ijk',20i3)
      tdloop: do v=1,size(box%all_ijk)
         if(k.eq.box%all_ijk(v)) then
! k is part of v_ij, this assignment actually redundant
            d_ijk=1; goto 17
         endif
      enddo tdloop
! varkappa independent of quad k
      box%dvk_ij(k)=0.0d0
      box%dvk_ji(k)=0.0d0
      cycle kloop
!
!       nonzero: if(d_ijk.eq.1) then
17    continue
! this varkappa depend on quad fraction k, calculate derivative
!      write(*,20)'vk_ij ',v,box%denominator
20    format('Denominator of ',a,' wrt quad ',i3,2x,1pe12.4) 
      t1loop: do v=1,size(box%ivk_ij)
         if(k.eq.box%ivk_ij(v)) then
            d_ij=1; exit t1loop
         endif
      enddo t1loop
!      if(d_ij.eq.1) write(*,30)'ivk_ij loop ',v,box%vk_ij
30    format('Numerator ',a,' wrt quad ',i3,1pe12.4)
      t2loop: do v=1,size(box%jvk_ji)
         if(k.eq.box%jvk_ji(v)) then
            d_ji=1; exit t2loop
         endif
      enddo t2loop
!      write(*,35)v,d_ijk,d_ij,d_ji
35    format('All d_xyz: ',i3,4i4)
!      if(d_ji.eq.1) write(*,30)'jvk_ji loop ',v,box%vk_ji
! Note that vk_ij and vk_ji are already divided by denominator
      box%dvk_ij(k)=(d_ij - box%vk_ij)/box%denominator
      box%dvk_ji(k)=(d_ji - box%vk_ji)/box%denominator
   enddo kloop
! debug output of the derivatives
   if(mqmqdebug) then
      do k=1,mqmqa_data%nquad
         write(*,100)k,box%dvk_ij(k),box%dvk_ji(k)
      enddo
100   format('3XQ In dexcess_dq: dvk_ij, dvk_ji wrt quad: ',i3,2(1pe14.6))
   endif
! now derivatives of xi with respect to quads      NOT DONE ???????
!
   return
 end subroutine dexcess_dq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine set_ternary_asymmetry(line)
!\begin{verbatim}
 subroutine set_ternary_asymmetry(line)
   implicit none
   character*(*) line
! to set asymmetries in a text
!\end{verbatim}
   integer i,j,ip,iq,ia,ib,ic,mm,icc(3),nc,kk,vz,toop(3)
   integer iph,ics,icon,ipm
   double precision mass
   character missingcon*60
   character phase*24,const(3)*24,asymcode*6,asymoc*3
   type(gtp_phaserecord), pointer :: phrec
!
   phase=' '
! called from gtp3E around line 5493
   write(*,10)trim(line)
10 format('3XQ set_ternary_asymmetry called from gtp3E: "',a,'"')
!   write(*,*)'3E set_ternary_asymmetry to be fixed'
! extract constituent indices and call setsym'
! text is extacted from frist nonblank position ip to first space
! first the phase, then 3 constituents finally the asymcode
   ip=0
   call getext(line,ip,2,phase,' ',iq)
   write(*,20)trim(phase),iq
20 format('3XQ Phase name: ',a,5x,i3)
   if(phase(1:1).ne.' ') then
      call find_phase_by_name(phase,iph,ics)
      if(gx%bmperr.ne.0) then
         write(*,21)trim(phase)
21       format(/'3XQ Ternary asymmetries for phase "',a,&
              '" ignored as phase not selected')
         gx%bmperr=0
         goto 1000
      endif
   else
      goto 1100
   endif
   nc=size(phlista(iph)%constitlist)
!   write(*,*)'3XQ in set_ternary_asymmetry, found phase ',iph,nc
!
! debug listing of mqmqa_data%contyp
!   do vz=1,nc
!      write(*,33)(mqmqa_data%contyp(i,vz),i=1,14)
!              4    5=FNN   6-7    8-9  same   11-12  13-14
33    format('3XQ: ',4i3,2x,i3, 2x,2i3,2x,2i3,2x,i3,2x,2i3,2x,2i3)
!   enddo
!
   extract_asymmetries: do while(ip.lt.len(line))
! save constituent indices in icc
      icc=0
      asymcode=' '; const=' '; missingcon=' '
      ipm=1
      find3: do i=1,3
         call getext(line,ip,2,const(i),' ',iq)
         if(gx%bmperr.ne.0) goto 1000
! The terminator is usually "!" or "/"
         if(const(i).eq.'!' .or. const(i).eq.'/') goto 1000
! we have to go through the list of constituents of the phase to find the
! constituents as we need their sequental index .... SUCK
!         mm=len_trim(const(i))+1 as constituents in MQMQA the name has suffix
! A small risk that it is an abbreviation ...
         mm=len_trim(const(i))
         if(mm.le.1) goto 1000
         compare: do j=1,nc
! the names of the constituents are in the species structure!
            kk=phlista(iph)%constitlist(j)
!            write(*,50)const(i)(1:mm),splista(kk)%symbol(1:mm)
50          format('Comparing "',a,'" and "',a,'" ',i3)
            if(const(i)(1:mm).eq.splista(kk)%symbol(1:mm)) then
! Hmmmmm, it is not species index we want, we want the number of the
! this species is phase constituent i, use mqmqa_data%contyp(5,i) !!
! in mqmqa_data%contyp(5,j) is the pair index !!??
               vz=mqmqa_data%contyp(5,j)
!               write(*,53)'3XQ Found asymmetric: ',i,splista(kk)%symbol,j,vz
53             format(a,i3,2x,a,5i4)
               icc(i)=vz
!               write(*,60)i,const(i),vz
60             format('3XQ cation index: ',i3,2x,a,2x,i4)
               cycle find3
            endif
         enddo compare
!         write(*,*)'3XQ Asymmetric constituent not found',i,const(i)(1:mm)
         missingcon(ipm:)=const(i)
         ipm=len_trim(missingcon)+2
! we have to read all constituents becaise some asymmetries may involve
! constitutents not selected
      enddo find3
! error if end of line
      if(ip.ge.len_trim(line)) goto 1100
      call getext(line,ip,2,asymcode,' ',iq)
! TKK means the third quad has the Toop element
! KTK means the second quad has the Toop element
! KKT means the first quad has the Toop element
! convert for example T3KT3 to TKK related to the 3 elements in the quads
! the value to be saved is the element index of the first, second or third quad
!      write(*,*)'3XQ Asymmetric code: ',asymcode
! if any icc is 0 skip
      do j=1,3
         if(icc(j).eq.0) then
            write(*,*)'3XQ missing asymmetry constituents: ',trim(missingcon)
            cycle extract_asymmetries
         endif
      enddo
!      write(*,100)trim(phase),icc(1),icc(2),icc(3),trim(asymcode)
100   format('3XQ asymmetric ternary in ',a,' elements ',3i4,5x,a)
! convert full asymmetry to OC
      call convert_asymm(asymcode,asymoc,icc,toop)
! output from convert_asymm seems OK but with some redundat data
! there can be several asymmetric ternaries
!      write(*,*)'3XQ Arrange actual order of cations in setasym'
! do we need toop?
!      call setasym(iph,icc,toop,nquad,asymoc)
      call setasym(iph,icc,toop,asymoc)
      if(gx%bmperr.ne.0) goto 1000
!      stop 'debug'
   enddo extract_asymmetries
!
   newXupdate=newXupdate+1
! tersym is declared globally, it should be within a phase record
! as each phase can have ternary symmetries
!
1000 continue
   return
1100 write(*,1110)line(min(1,ip-10):ip+10)
1110 format('Problem extracting ternary asymmetry: ',a)
   gx%bmperr=4499
   goto 1000
 end subroutine set_ternary_asymmetry

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine convert_asymm
!\begin{verbatim}
 subroutine convert_asymm(code1,code2,icc,toop)
! convert from 6 to 3 letters
   implicit none
   character code1*6,code2*3
! icc are the 3 cations in the ternary ... toop is ?
   integer icc(3),toop(3)
!\end{verbatim}
   character cha*1,chb*1
   integer ia,ib,iv,iw,ntoop
   code2='KKK'
   iw=0
   ntoop=0
   toop=0
   do iv=1,3
      iw=iw+1
      if(code1(iw:iw).eq.'T') then
         ntoop=ntoop+1
         ia=ichar(code1(iw+1:iw+1))-ichar('0')
         if(ia.gt.0 .and. ia.le.3) then
! The T is followed by a digit indicating the constituent, 1, 2 or 3
            toop(ntoop)=iv
            code2(4-ia:4-ia)='T'
! skip one position in code1
            iw=iw+1
         else
! toop elementet is indicated by the 3 cation positions
            toop(ntoop)=iv
            code2(4-iv:4-iv)='T'
         endif
!         write(*,*)'3XQ Toop cation is: ',icc(ia),' position ',toop(ntoop)
      endif
!      write(*,10)code1,code2,iv,icc
   enddo
!   write(*,10)code1,code2,icc,toop
10 format('3XQ asymm: "',a,'" to "',a,'" cations: ',3i2,' Toop: ',3i3)
   return
 end subroutine convert_asymm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine setasym
!\begin{verbatim}
 subroutine setasym(iph,icc,toop,kkk)
! set asymmetry of a ternary
! the cation indices can be in any order, must be ordered.
   implicit none
   integer iph,icc(3), toop(3)
   character*3 kkk
! 
! REDUNDANT AS INTEGRATED IN PMON6.F90
!
!\end{verbatim}
   integer i,j,k,dim3,ntercat
   integer a,b,c,mm,v
! default is 'KKK' which is symmetrical for the 3 binaries 1-2, 1-3 and 2-3
!            'TKK' means element 3 is asymmetrical for 1-2
!            'TKT' means element 3 is asymmetrical for 1-2 and element 1 for 2-3
!
! icc are the cation indices, toop is zero unless one or more toop cations
!   write(*,60)icc,toop,kkk
60 format('3XQ ENTERING SETASYM: icc: ',3i3,' toop: ',3i3,' kkk: ',a)
!  format(a,3i3,3x,3i3,2x,a)
!
   i=icc(1); j=icc(2); k=icc(3)
!
!   write(*,*)'3XQ calling order3KKK '
!
   call order3KKK(i,j,k,a,b,c,kkk)
!
   if(a.lt.0) then
      write(*,*)'Problems 10 in order3 ',i,j,k,a,b,c
      stop
   endif
! rearranged i, j, k
!   write(*,70)a,b,c,kkk
70 format('3XQ rearranged order in setasym: ',3i4,5x,a)
!
! if order changed, change KKK, assume only one T
!
! any phase may have asymmetric ternaries but at present only MQMQA
!   if(.not.allocated(phlista(iph)%tersys)) then
!      stop
!   endif
!
! emergency ... should be checked, a system with 3 constituent has 1 ternary
!   dim3=size(tersys)
!   write(*,333)a,b,c,mqmqa_data%ncat,mm,size(tersys)
333 format('3XQ In setasym: ',8i4)
   mm=terind(a,b,c)
   if(mm.le.0) then
      write(*,*)'3XQ terind cannot find this system',a,b,c
      stop
   end if
!   write(*,333)a,b,c,mqmqa_data%ncat,mm,size(tersys)
!
   newXupdate=newXupdate+1
! tersym is declared globally, it should be within a phase record
! as each phase can have ternary symmetries
!   write(*,511)mm,' old ',tersys(mm)%asymm,tersys(mm)%isasym,a,b,c
511 format('3XQ in setasym ternary: ',i4,a,' asymmetry <',a,'>   ',3i3,2x,3i3)
   tersys(mm)%asymm=kkk
   tersys(mm)%isasym=0
! or should one use i, j, k ???
! the indices in tersys(mm)%el are the 3 element indices of the ternary
!    write(*,300)mm,tersys(mm)%el
300 format('Element numbers in ternary ',i3,' are ',3i3)
   if(kkk(1:1).eq.'T') tersys(mm)%isasym(1)=tersys(mm)%el(3)
   if(kkk(2:2).eq.'T') tersys(mm)%isasym(2)=tersys(mm)%el(2)
   if(kkk(3:3).eq.'T') tersys(mm)%isasym(3)=tersys(mm)%el(1)
   write(*,511)mm,' new ',tersys(mm)%asymm,tersys(mm)%isasym,a,b,c
!
! for debugging list whole array
!   write(*,310)dim3
310 format(/'Listing of the ',i3,' ternary systems and their asymmetry',&
         /'  i  seq   cat1 cat2 cat3       T/0 T/0 T/0    asymmetry code')
!   ntercat=mqmqa_data%ncat*(mqmqa_data%ncat-1)*(mqmqma_data%ncat-2)/6
!   do i=1,ntercat
!      write(*,320)i,tersys(i)%seq,(tersys(i)%el(j),j=1,3),&
!           tersys(i)%isasym,tersys(i)%asymm
320   format(i3,i5,2x,3(1x,i4),5x,3i4,5x,a)
!   enddo
!   write(*,330)
330 format(/'Number in T/0 column is actual asymmetric element')
!   
   return
 end subroutine setasym

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine correlate_const_and_quads
!\begin{verbatim}
 subroutine correlate_const_and_quads(lokph)
! this subroutine should for each mqmqa constituent create their
! quad index element order for handling asymmetric variables in compvar
! el1  1 1 1 .. 1   ! 2   2  .. 2    ! 3 .. ! n-1
! el2  1 2 3 .. n-1 ! 2   3  .. n-1  ! 3 .. ! n-1
! quad 1 2 3    n   ! n+1 n+2   2n-1 ! 2n   ! (n-1)n/2
! With n elements and  one anion there are n-1 cations
! The anion element index can be any value from 1 to n
!
! i=el2ancat(j) is cation index of element j, a negative value mean anion
! i=con2quad(j) is index in quad fraction array of constituent j
!    it is populated using ijklx(cat1,cat2,1,1) where the 1 is the anion
! OC saves fractions in phase_varres%yfr(1..n) for a single sublattice
! there is no need to transfer fractions from quad to phase_varres%yfr
! A quad has 1 or 2 cations ALWAYS indexed from 1 .. n-1 (no anion fraction)
! i=el2ancat(j) is the cation index element j. If j is anion a negative value
! i=quadofel(j) is the cation index of an element i
! The anion element index is not used 
!  ... but its elllink is saved in xanione and element index in xanionalpha
!   
   implicit none
   integer iph,lokph,loksp,lokcs,nfr,isp,iel,jp,el1,el2,icon,endmem,mm
   integer cat1,cat2
   integer missing,ll,nocon
   logical noanion
   integer, allocatable, dimension(:) :: invert,inverse
   integer, allocatable, dimension(:) :: findan
   character*24 quadname
!
! called from create_asymmetry in gtp3B
!   write(*,7)lokph,nfr,noofel
7  format(/'3XQ In correlate_const_and_quad',3i5//)
!
   nfr=phlista(lokph)%nooffr(1)
   allocate(findan(noofel))
!   lokcs=phlista(lokph)%linktocs(1) composition set?
! note element numbers are not in order, the anion may be anywhere
!
! first step, find anion, it is present in all constituents
! Stupid to do this here, it has already been found but lost
   findan=0
   do isp=1,nfr
      loksp=phlista(lokph)%constitlist(isp)
      iel=size(splista(loksp)%ellinks)
!      if(iel.eq.3) then
!         write(*,3)isp,loksp,iel,splista(loksp)%symbol,&
!              (ellista(splista(loksp)%ellinks(jp))%symbol,jp=1,iel),&
!              (ellista(splista(loksp)%ellinks(jp))%alphaindex,jp=1,iel),&
!              (splista(loksp)%ellinks(jp),jp=1,iel)
!      else
!         write(*,2)isp,loksp,iel,splista(loksp)%symbol,&
!              (ellista(splista(loksp)%ellinks(jp))%symbol,jp=1,iel),&
!              (ellista(splista(loksp)%ellinks(jp))%alphaindex,jp=1,iel),&
!              (splista(loksp)%ellinks(jp),jp=1,iel)
!      endif
!2     format('3XQ const: ',3i3,2x,a12,2x,2(a,2x),4x,2(i3),5x,2(i3))
!3     format('3XQ const: ',3i3,2x,a12,2x,3(a,2x),3(i3),2x,3(i3))
      do jp=1,iel
         el1=splista(loksp)%ellinks(jp)
         findan(el1)=findan(el1)+1
      enddo
   enddo
!   write(*,4)'3XQ elements: ',findan
4  format(a,20i3)
! count the number of times an element occurs
!   write(*,22)(jp,elements(jp),ellista(jp)%alphaindex,ellista(jp)%symbol,&
!        jp=1,noofel)
22  format(/'3XQ elements :',10(3i2,1x,a,';')/)
   el1=0
! The anion should be present in all quads!
   do jp=1,noofel
      if(findan(jp).gt.el1) then
         el1=findan(jp); el2=jp;
      endif
   enddo
! el2 is the element index is in ellista, el1 is the alphabetical order
!   write(*,*)'3QX anion is element: ',el1,el2,findan(el2)
!      
   mqmqa_data%xanione=el2
   mqmqa_data%xanionalpha=ellista(el2)%alphaindex
!   write(*,6)ellista(mqmqa_data%xanione)%symbol,&
!        mqmqa_data%xanione,mqmqa_data%xanionalpha
6  format(/'3XQ line 3383 anion: ',a,' ellink: ',i3,' alphabetically: ',i3/)
!
! set up translation table for cations from ellink to 1..ncat
! the anion has a negative value in el2ancat, the cations index 1..ncat
!   write(*,*)'3XQ in correlate_const_and_quads ... allocating el2ancat'
! mqmqa_data is not allocated ... suck
!   if(allocated(mqmqa_data%el2ancat(noofel))) &
!        deallocate(mqmqa_data%el2ancat(noofel))
   allocate(mqmqa_data%el2ancat(noofel))
!   write(*,*)'Size of mqmqa_data%el2ancat ',size(mqmqa_data%el2ancat)
   do jp=1,noofel
      if(jp.lt.mqmqa_data%xanionalpha) then
!      if(jp.lt.mqmqa_data%xanione) then
         mqmqa_data%el2ancat(jp)=jp
      elseif(jp.gt.mqmqa_data%xanionalpha) then
!      elseif(jp.gt.mqmqa_data%xanione) then
         mqmqa_data%el2ancat(jp)=jp-1
      else
         mqmqa_data%el2ancat(jp)=-jp
      endif
!      write(*,*)'3xq mqmqa_data%el2cat: ',jp,mqmqa_data%el2ancat(jp)
   enddo
!   write(*,16)'3XQ Elements alphabetically:  ',&
!        ((ellista(elements(jp))%symbol),jp=1,noofel)
!   write(*,17)'3XQ Elements in ellista order:',(elements(jp),jp=1,noofel)
!   write(*,17)'3XQ Element alpha indices:    ',(jp,jp=1,noofel)
!   write(*,17)'3XQ Cation  alpha indices:    ',&
!        (mqmqa_data%el2ancat(jp),jp=1,noofel)
15  format(a,20(i2,a2))
16  format(a,20(1x,a2))
17  format(a,20i3)
!
! We need to know how to transfer compositions from phase_varres%yfr to xquad
   allocate(mqmqa_data%con2quad(nfr))
! loop though all constituents in the %constitlist, extract cations and
! calculate its index in the xquad.  Only done once!
   con2quad: do isp=1,nfr
      loksp=phlista(lokph)%constitlist(isp)
! there are 2 or 3 element links, one of which is an anion
      el1=ellista(splista(loksp)%ellinks(1))%alphaindex
      cat1=mqmqa_data%el2ancat(el1)
!      cat1=mqmqa_data%el2ancat(splista(loksp)%ellinks(1))
!      write(*,18)'First:   ',el1,cat1,cat1,splista(loksp)%symbol
18    format(a,3i4,5x,a)
      first: if(cat1.lt.0) then
! first link was to the anion, next must be a cation
         el1=ellista(splista(loksp)%ellinks(2))%alphaindex
         cat1=mqmqa_data%el2ancat(el1)
!         write(*,18)'Second:  ',el1,&
!              mqmqa_data%el2ancat(splista(loksp)%ellinks(2)),cat1
         more1: if(size(splista(loksp)%ellinks).gt.2) then
! there can be 1 or 2 cations, the first ellink was to an anion
            el1=ellista(splista(loksp)%ellinks(3))%alphaindex
            cat2=mqmqa_data%el2ancat(el1)
!            write(*,18)'Third:   ',splista(loksp)%ellinks(3),&
!                 mqmqa_data%el2ancat(splista(loksp)%ellinks(3)),cat2
         else
! if there is no third element the single cation is doubled
            cat2=cat1
         endif more1
      else
! we found one cation, the next ellink can be an anion or cation        
         el1=ellista(splista(loksp)%ellinks(2))%alphaindex
         cat2=mqmqa_data%el2ancat(el1)
!         write(*,18)'Fourth:  ',el1,&
!                 mqmqa_data%el2ancat(el1),cat2
         second: if(cat2.lt.0) then
            more2:if(size(splista(loksp)%ellinks).gt.2) then
! there can be 1 or 2 cations, the second ellink can be to the anion
               el1=ellista(splista(loksp)%ellinks(3))%alphaindex
               cat2=mqmqa_data%el2ancat(el1)
!               write(*,18)'Fifth:   ',el1,&
!                    mqmqa_data%el2ancat(el1),cat2
            else
! the single cation is doubled
               cat2=cat1
            endif more2
         endif second
      endif first
! when we come here we hav one or two cations
      mqmqa_data%con2quad(isp)=ijklx(cat1,cat2,1,1)
!      write(*,55)isp,cat1,cat2,mqmqa_data%con2quad(isp)
55    format('3xq loop: ',i3,2x,2i3,2x,i5)
   enddo con2quad
! allocate also array with all A/X quads
   allocate(mqmqa_data%emquad(mqmqa_data%ncat))
! enter data in emquad
   cat1=1
   cat2=mqmqa_data%ncat
   do isp=1,mqmqa_data%ncat
      mqmqa_data%emquad(isp)=cat1; cat1=cat1+cat2; cat2=cat2-1
   enddo
! list quads (why?)
   write(*,68)(mm,mm=1,mqmqa_data%nquad)
68 format('3XQ quads: ',21i3)
   write(*,57)'3XQ emquads:',(mqmqa_data%emquad(isp),isp=1,mqmqa_data%ncat)
57 format(a,25i4)
!
! loop for all constituents of the mqmqa phase
! we should populate all structures of the %alphaindex of the element
! skipping the alphaindex of the anion
!   write(*,60)'3XQ constituents   :',(jp,jp=1,nfr),&
!              'mqmqa_data%con2quad:',(mqmqa_data%con2quad(jp),jp=1,nfr)
60 format(/a,10(i3,1x)/a,10(i3,1x))
! icon is index of constituent in phase 1..n
! splista(icon)%symbol is species symbol
!   write(*,65)
65 format(/'3XQ Constituents in alphabetical order:')
!   write(*,70)(trim(splista(phlista(lokph)%constitlist(jp))%symbol),jp=1,nfr)
!
!   write(*,*)'3XQ Constituents in quad order:'
!   write(*,70)(trim(splista(phlista(lokph)%constitlist(mqmqa_data%con2quad(jp)))%symbol),jp=1,nfr)
!
!70 format('3XQ: ',10(a,', '))
!71 format('3XQ: ',2i3,3x,a)
!
   allocate(inverse(nfr))
!   write(*,87)
87 format(/'3XQ     OC fraction order   MQMQA quad order')
  do jp=1,nfr
      qqq: do el1=1,nfr
         cat1=mqmqa_data%con2quad(el1)
         if(cat1.eq.jp) then
            quadname=splista(phlista(lokph)%constitlist(el1))%symbol
            inverse(jp)=el1
            exit qqq
         endif
      enddo qqq
!      write(*,88)jp,trim(splista(phlista(lokph)%constitlist(jp))%symbol),&
!           el1,trim(quadname)
88    format('Order ',i3,3x,a12,i5,2x,a)
   enddo
!
! this is if we need to convert from xquad array to yfr
!   write(*,89)(inverse(jp),jp=1,nfr)
89 format('3XQ Quad2con: ',20i3)
1000 continue
   return
 end subroutine correlate_const_and_quads

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine list_quads(kk)
!\begin{verbatim}
 subroutine list_quads(kk)
! emergency subroutine because phlista protected in pmon6
   implicit none
   integer kk
!\end{verbatim}
   integer nel
!
   kk=0
! negative is anion
   write(*,2)(ellista(elements(nel))%symbol,nel=1,noofel)
2  format(/'3XQ Element names:      ',20(a2,1x))
   do nel=1,noofel
      if(mqmqa_data%el2ancat(nel).lt.0) kk=nel
   enddo
! elements as quad numbers
   write(*,3)size(mqmqa_data%el2ancat),mqmqa_data%el2ancat
3  format('3XQ Cation indices:',i2,2x,20i3)
!3  format('3XQ el2ancat:     ',i3,2x,20i3)
   if(kk.eq.0) then
      write(*,*)'You have a strange MQMQA system without any anion'
   else
      write(*,4)ellista(elements(kk))%symbol,mqmqa_data%xanionalpha,&
           mqmqa_data%xanione
4     format('3XQ The anion element name, index and link: ',a,2i3)
   endif
!
   return
 end subroutine list_quads

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine listconst
!\begin{verbatim}
 subroutine listconst(iph)
! emergency subroutine because phlista protected in pmon6
   implicit none
   integer iph
!\end{verbatim}
   type(gtp_phase_varres), pointer :: phres
   integer lokph,lokcs,isp,iel,elx(4),elxx(4),jp,j1,j4,nel,cations(2),jj,kk
   character :: elsym(noofel)*2

   elsym=' '
   kk=0
   do nel=1,noofel
      elsym(nel)=ellista(elements(nel))%symbol
      if(mqmqa_data%el2ancat(nel).lt.0) kk=nel
   enddo
   write(*,2)(elsym(jj),jj=1,noofel)
2  format(/'Element names:      ',20(a2,1x))
   write(*,3)size(mqmqa_data%el2ancat),mqmqa_data%el2ancat
3  format('3XQ el2ancat: ',i3,2x,20i3)
   if(kk.eq.0) then
      write(*,*)'You have a strange MQMQA system without any anion'
   endif
   write(*,4)elsym(kk),mqmqa_data%xanionalpha,&
        mqmqa_data%xanione
4     format('3XQ The anion element name, index and link: ',a,2i3)
!
   lokph=phases(iph)
   lokcs=phlista(lokph)%linktocs(1)
   isp=0
!   mqmqa_data%xanione=splista(j4)%ellinks(nel)   
!   write(kou,5)lokcs,mqmqa_data%xanione,mqmqa_data%xanionalpha
   write(kou,5)
5  format(/'Con  Quad Nel Elements      Elem index',2x,'Species name',&
        15x,'Cations')
   specie: do jp=1,phlista(lokph)%nooffr(1)
      isp=isp+1
      j4=phlista(lokph)%constitlist(jp)
      nel=size(splista(j4)%ellinks)
      elsym='  '
      elxx=1000
      jj=0
      element: do iel=1,nel
         elx(iel)=splista(j4)%ellinks(iel)
         elsym(iel)=ellista(elx(iel))%symbol
         elxx(iel)=ellista(splista(j4)%ellinks(iel))%alphaindex
         j1=mqmqa_data%el2ancat(elxx(iel))
         if(j1.gt.0) then
! con2cat(i) is the cation index of i, negative if anion
            jj=jj+1
            cations(jj)=j1
         endif
      enddo element
      if(jj.eq.1) cations(2)=cations(1)
      if(noofel.le.3) then
         write(kou,19)isp,mqmqa_data%con2quad(isp),nel,(elsym(kk),kk=1,3),&
              elxx,splista(j4)%symbol,cations
19       format(i3,i4,2x,i4,1x,3(a,2x),4(1x,i2),2x,a,2x,2i3)
      else
         write(kou,20)isp,mqmqa_data%con2quad(isp),nel,(elsym(kk),kk=1,4),&
              elxx,splista(j4)%symbol,cations
20       format(i3,i4,2x,i4,1x,4(a,2x),4(1x,i2),2x,a,2x,2i3)
      endif
   enddo specie
   write(*,*)'The quads are in the alphabetical order of the quad elements'
   return
 end subroutine listconst

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine listpartree
!\begin{verbatim}
 subroutine listpartree(lokph)
! list all endmembers and excess parameter records for a phase
! in order to understand the MQMQX phase
   implicit none
   integer lokph
!\end{verbatim}
   !\end{verbatim}
   type (gtp_endmember), pointer :: endmemrec,em
   type (gtp_interaction), pointer :: intrec
   type (gtp_property), pointer :: proprec
   type (gtp_asymprop), pointer :: asymdata
   integer intlevel,nofr,fracs(10),npr,intsave,ii,nint,powers(3)
!   double precision vals(6)
   character*3 tab1
   character*6 tab2
   character*9 tab3
!
   type stack
      type(gtp_interaction), pointer :: current
   end type stack
   type(stack), dimension(:), allocatable :: intstack
!   
   intlevel=0; fracs=0
   tab1='---'
   tab2='------'
   tab3='---------'
   allocate(intstack(5))
   write(*,5)
5  format('3XQ list of the excess parameter tree')
!   
   endmemrec=>phlista(lokph)%ordered
!   if(associated(endmemrec)) write(*,*)'3XQ there is an endmember'
   emloop: do while(associated(endmemrec))
      nofr=1
      fracs(nofr)=endmemrec%fraclinks(1,1)
      intrec=>endmemrec%intpointer
!
      if(associated(intrec)) write(*,10)fracs(1)
10    format('3XQ interactions from endmember ',i3)
      intsave=0
      nofr=nofr+1
      intloop:do while(associated(intrec))
         fracs(nofr)=intrec%fraclink(1)
         proprec=>intrec%propointer
         if(.not.associated(proprec)) then
            write(*,20)intsave+1,fracs(nofr)
20          format('3XQ interaction record level',i3,', constituent',i3)
         else
            proploop: do while(associated(proprec))
               if(.not.associated(proprec%asymdata)) then
                  powers=0
               else
                  powers(1)=proprec%asymdata%ppow
                  powers(2)=proprec%asymdata%qpow
                  powers(3)=proprec%asymdata%rpow
               endif
               npr=proprec%antalprop
               if(intsave.eq.0) then
                  write(*,100)' ',intsave+1,npr,powers,(fracs(ii),ii=1,nofr)
               elseif(intsave.eq.1) then
                  write(*,100)tab1,intsave+1,npr,powers,(fracs(ii),ii=1,nofr)
               elseif(intsave.eq.2) then
                  write(*,100)tab2,intsave+1,npr,powers,(fracs(ii),ii=1,nofr)
               elseif(intsave.eq.3) then
                  write(*,100)tab3,intsave+1,npr,powers,(fracs(ii),ii=1,nofr)
               else
                  write(*,100)'---',intsave+1,npr,powers,(fracs(ii),ii=1,nofr)
               endif
100            format('3XQ ',a,' at level ',i1,', func: ',i2,&
                    ', powers: ',3i2,', constituents ',9i3)
               proprec=>proprec%nextpr
            enddo proploop
         endif
         if(associated(intrec%highlink)) then
! save intrec%nextlink and jump to higher level
            intsave=intsave+1
            intstack(intsave)%current=>intrec%nextlink
            intrec=>intrec%highlink
            nofr=nofr+1
            fracs(nofr)=intrec%fraclink(1)
         else
! check the nextlink, pop saved if empty
            intrec=>intrec%nextlink
            pop: do while(.not.associated(intrec))
               write(*,*)'3XQ pop stack'
               if(intsave.gt.0) then
                  intrec=>intstack(intsave)%current
                  intsave=intsave-1
                  nofr=nofr-1
               else
                  exit intloop
               endif
            enddo pop
            cycle intloop
         endif
! if we come here there are no more interaction records for this endmember
      enddo intloop
!
!      write(*,*)'3XQ next endmember'
      endmemrec=>endmemrec%nextem
!
   enddo emloop
   write(*,*)'No more parameters'
1000 continue
      return
 end subroutine listpartree

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!addtotable subroutine varkappa
!\begin(verbatim} subroutine varkappa
 subroutine varkappadefs(phres)
! subroutine to list \varakappa, \xi and y_ik definitions
   type(gtp_phase_varres), pointer :: phres
!   type(gtp_phase_varres), pointer :: mqmqavar
!\end{verbatim}
   integer nv,cat1,cat2,i,j,k,ip,i0,jp
   type(gtp_mqmqa_var), pointer :: mqf
   type(gtp_allinone), pointer :: box
   character*80 line1,line2,qline
   character*2, dimension(:), allocatable :: quadcat
!
   mqf=>phres%mqmqaf
! copied from pmon
!   write(kou,4124)mqmqa_data%nquad,mqmqa_data%ncat
4124 format(/'The ',i3,' quads for ',i2,' cations are arranged ',&
          'in order of the n cations:'/&
          'Quad  ',9x,'1   2  ...  n | n+1 n+2 ... 2n-1 | 2n .. | n(n+1)/2'/&
          'Cation',9x,'1   1  ...  1 | 2   2   ...  2   | 3  .. | n'/&
          'Cation',9x,'1   2  ...  n | 2   3   ...  n   | 3  .. | n')
!
! identify the actual cations in all quads as above
!   write(*,50)(i,i=1,mqmqa_data%nquad)
! create the quadcat indeices used for the vk_ij quad dependences 
   allocate(quadcat(mqmqa_data%nquad))
   line1='Cat1:'
   ip=6
   i0=ichar('0')
   k=1
! To fix problems here see around line 4100 about box%ivk_ij, %jvk_ji %kvk_ijk
   do i=1,mqmqa_data%ncat
      do j=i,mqmqa_data%ncat
         line1(ip:ip+2)='  '//char(i0+i)
         quadcat(k)(1:1)=char(i0+i)
         k=k+1
         ip=ip+3
      enddo
      ip=ip+1
   enddo
51 format(a)
   line2='Cat2:'
   qline='Quad:'
   ip=6
   k=1
   i0=ichar('0')
   do i=1,mqmqa_data%ncat
      do j=i,mqmqa_data%ncat
         line2(ip:ip+2)='  '//char(i0+j)
         quadcat(k)(2:2)=char(i0+j)
         if(k.lt.10) then
            qline(ip:ip+2)='  '//char(i0+k)
         else
            qline(ip:ip+2)=' 1'//char(i0+k-10)
         endif
         ip=ip+3
         k=k+1
      enddo
      ip=ip+1
   enddo
! nice output of quads and cation dependencies
   write(*,51)trim(qline)
   write(*,51)trim(line1)
   write(*,51)trim(line2)
!
! quadcat(k)(1:2) are the 2 cation indices (as characters) in quad k
! ivk_ij, ivk_ji, kvk_ijk arrays of quad indices indices
!   vkloop: do nv=1,size(mqf%compvar)
!      box=>mqf%compvar(nv)
! box%ivk_ij(1..n) are indices of quads to be added 
!      write(*,100)'vk_ij',(box%ivk_ij(cat1),cat1=1,size(box%ivk_ij))
!      write(*,100)'vk_ji',(box%jvk_ji(cat1),cat1=1,size(box%jvk_ji))
!      write(*,100)'denom',(box%kvk_ijk(cat1),cat1=1,size(box%kvk_ijk))
!100   format(a,10i3)
!   enddo vkloop
   write(*,99)
99 format('3XQ some ternary asymmetries may still be wrong')
   vkloop2: do nv=1,size(mqf%compvar)
! _ij
      box=>mqf%compvar(nv)
      line1='x_'//quadcat(box%ivk_ij(1))
      ip=len_trim(line1)+1
      k=2
      do while(k.le.size(box%ivk_ij))
         line1(ip:)='+x_'//quadcat(box%ivk_ij(k))
         k=k+1
         ip=ip+5
      enddo
! To fix problems here see around line 4100 about box%ivk_ij, %jvk_ji %kvk_ijk
      write(*,105)'vk_'//char(i0+box%cat1)//char(i0+box%cat2)//' = '//&
           trim(line1)
! _ji
      line2='x_'//quadcat(box%jvk_ji(1))
      ip=len_trim(line2)+1
      k=2
      do while(k.le.size(box%jvk_ji))
         line2(ip:)='+x_'//quadcat(box%jvk_ji(k))
         k=k+1
         ip=ip+5
      enddo
      write(*,105)'vk_'//char(i0+box%cat2)//char(i0+box%cat1)//' = '//&
           trim(line2)
! _denom
! NOTE some quad fractions appear twice!! should be removed
      qline=trim(line1)//'+'//trim(line2)//' +x_'//quadcat(box%kvk_ijk(1))
      ip=len_trim(qline)+1
      k=2
      do while(k.le.size(box%kvk_ijk))
         qline(ip:)='+x_'//quadcat(box%kvk_ijk(k))
         k=k+1
         ip=ip+5
      enddo
      write(*,105)'denom = '//trim(qline)
105   format(a)
!
!      write(*,110)'vk_',box%cat1,box%cat2,quadcat(box%ivk_ij(1)),
!           ((quadcat(box%jvk_ji(cat1)),cat1=2,size(box%jvk_ji))
!
!      write(*,110)'vk_',box%cat1,box%cat2,&
!           (quadcat(box%jvk_ji(cat1)),cat1=1,size(box%jvk_ji))
!      write(*,120)'denom = vk_ij+vk_ji + ',&
!           (quadcat(box%kvk_ijk(cat1)),cat1=1,size(box%kvk_ijk))
!110   format(a,2i1,' = x_',a,'+'))
!120   format(a,20('x_',a,'+'))
   enddo vkloop2
1000 continue
   return
 end subroutine varkappadefs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

