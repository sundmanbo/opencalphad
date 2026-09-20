!
! gtp3XQ for for MQMQA
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
!>     15B. Section: calculate G and other things for MQMQA and Toop/Kohler
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
!
! Remaining things to do:
! *. fix the lost control of ternary asymmetries
! *. When TDB file generates ternaries make sure ternary quads are species
! *. Make sure TDB file can set ternary symmetries
! *. Check expression for the asymmetrik ivk_ij etc ... savenu, savegamma!!!
! 1. Write the code for the ternary_factor for ternary parameters
! 2. Test calculations of several systems and compare with FS
!
! I have added a new routine to set a single asymmetric ternary secifying
! one corner as the Toop quad and generate the vk_ij for the two binary sides.
! I believe the previous routine was wrong because it added asymmetry
! to all binaries attached to a Toop quad.  But I have not finished this.
!
!-----------------------------------------------------
! I lost control of these subroutines, a list with some info:
! config_entropy_mqmqa written 2020, no recent changes
! calc_mqmqa           2026 excess
! calc_toop            old, still used by one example, keep
! new_mqmqa_excess     new excess
! dvkij_dzijk          calculates partial derivatives of binary excess old?
! mqmqa_excesspar_name writes a excess parameter name with all constituents
! ternary_factor1      calculates the ternary parameter composition dependence
! ternary_factor2      calculates the ternary parameter TP dependance
! convert_y2quadx      ??
! set_quadfractions    copy OC yfr (CEQ) to quad fractions
! init_excess_asym     calculates composition variable values incl asymmetries
!                      or maybe setup the logics to calculate them?
! pairfracs            updates values of y_i/k and \xi_ij ??
! ibin function        find binaries associated with a ternary
! ijklx function       finds quad index of quadfraction i,j,k,l
! binsym function      sequential index of binary system ??
! order3               rearranges mqmqa constituent indices
! terind function      finds ternary asymmetry record
! new_ternary_asym     new code to set asymmeric constists of a single ternary
! calcasymvar          calculates vk_ij, x_ij and y_i/k (and derivatives?)
!                      and allocates ivk_ij etc from quads
! varkappa1            set values of vk_ij, xi_ij and y_k from quadfractions
! set_ternary_asymmetry sets a ternary asymmetry ??
! correlate_const_and_quads called from gtp3B when reading a TDB
! list_quads           handle listing from pmon6 (user i/f)
! list_quads_with_single_cation handle listing from pmon6 (user i/f)
! listconst            list constituents
! listpartree          list the parameters of an mqmqa phase
! quadprops            called from pmon6
! list_ternary_cations listing of vk_ij, xi_ij and y_k expressions
! list_mqmqa_variables list MQMQA variables expressions and values
! -------------------- no more
!
! AI Calude note:
! 2026-08-31 dead-code cleanup: removed calc_ternarymq, calc_newdvkij_values,
! add_ternary_asym, idonotunderstand, varkappa7, extract_asym, convert_asymm,
! setasym, order3KKK, list_quads_short, toop_ternary -- all had zero call
! sites anywhere in the tree (checked incl. comments/doc markers, and the
! order3KKK/setasym pair which only called each other). Full pre-cleanup
! copy kept at Claude2/gtp3XQ.F90.before-cleanup.
!

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
! Problem with fq=20 ...
   integer, parameter :: fq=50, fz=50, f1s=50
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
! contyp is initiated in gtp3B.F90 when reading the database
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
! 2026-05-11 diag: dchain accumulates the would-be chain-rule contribution
!                  to dsquad without applying it, so we can see its size
   double precision dchain(fq)
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
! ensure ternary paifractions is updated
   mqmqa_data%mqmqa_terasym1=mqmqa_data%mqmqa_terasym1+1
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
! 2026-05-10 fix: pqq(1) was pinq(pq(2)) (typo, both used pq(2)).
! That mis-mapped cation 1 of every cross quad to cation 2's pair record,
! corrupting yy1, ceqf1, and the pair-entropy term.  Only had a visible
! effect when two cations share an element (NOTOK).
         pqq(1)=mqmqa_data%pinq(pq(1))
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
   dchain=zero  
! 2026-05-11 diag: track would-be chain-rule contribution to dsquad
! replaced s1 by q1
   quadloop: do q1=1,ncon
      if(q1.ne.mqmqa_data%contyp(10,q1)) then
! TEST: the value in contyp(10,q1) should be q1 ...  260111/BoS WHY??
         write(*,441)q1,mqmqa_data%contyp(10,q1),mqmqa_data%contyp(14,q1)
441      format('3XQ problems in %contyp with quad indexing:',3i5)
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
! 2026-05-11 diag: compact summary of entropy components and gradients
   if(ddebug) then
      write(*,'(a,9(1pe12.4))')'3XD y    :',(phvar%yfr(s1),s1=1,ncon)
      write(*,'(a,9(1pe12.4))')'3XD ssub :',ssub,(dssub(s1),s1=1,ncon)
      write(*,'(a,9(1pe12.4))')'3XD send :',send,(dsend(s1),s1=1,ncon)
      write(*,'(a,9(1pe12.4))')'3XD squad:',squad,(dsquad(s1),s1=1,ncon)
      write(*,'(a,9(1pe12.4))')'3XD chain:',zero,(dchain(s1),s1=1,ncon)
      write(*,'(a,9(1pe12.4))')'3XD total:',ssub+send+squad,&
           (dssub(s1)+dsend(s1)+dsquad(s1),s1=1,ncon)
      write(*,'(a,9(1pe12.4))')'3XD pair :',(pair(s1),s1=1,noofpair)
      write(*,'(a,9(1pe12.4))')'3XD ceqf1:',(ceqf1(s1),s1=1,nspin(1))
      write(*,'(a,9(1pe12.4))')'3XD ceqf2:',(ceqf2(s1),s1=1,nspin(2))
   endif
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
!   double precision, dimension(:,:), allocatable :: dvals(:,:),affarr(:)
   double precision, dimension(:,:), allocatable :: dvals(:,:)
   double precision, dimension(:), allocatable :: affarr(:)
! for saving FNN reference energies
   double precision refg(f1,f1)
   double precision dummy1,dummy2
! for MQMQA minimal fractions
   double precision, parameter :: MINMQMQA=1.0D-5
   TYPE(gtp_parcalc) :: gz
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_endmember), pointer :: endmemrec
   TYPE(gtp_interaction), pointer :: intrec,terrec
   TYPE(gtp_pystack), pointer :: pystack
   TYPE(gtp_phase_add), pointer :: addrec
   TYPE(gtp_mqmqa_var), pointer :: mqf
   TYPE(gtp_tooprec), pointer :: tooprec
! for handling excess parameters, just binary, use no mqmqa_data ksi arrays
   integer ij,jd,jq,qq1,qq2,ass,mpow,isumx,tsize,tch,iiz,mqmqcon,mqmqjy
   integer noofex,nqx,ncv,icv,dd
   double precision ksi,sumx,dsumx
   double precision dksi(3),d2ksi(3)
!   logical calc_alldvkij
!   logical ddebug
   logical :: oldmqmqa_model = .true.
!   save oldmqmqa_model
!------------------------------------- 
! tch is level of debug output, 0=none, 3=max
   tch=0
   noofex=0
!   calc_alldvkij=.TRUE.
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
!      write(*,86)mqmqa_data%ncat,(mqf%y_ik(icv),icv=1,mqmqa_data%ncat)
!86    format(/'3XQ y_i/k ',i2,': ',(7F9.6))
!      do icv=1,ncv
!         write(*,80)mqf%compvar(icv)%vk_ij,mqf%compvar(icv)%vk_ji,&
!              mqf%compvar(icv)%xi_ij,mqf%compvar(icv)%xi_ji
!80       format(2x,2(1pe16.8),5x,2(1pe16.8))
!      enddo
!      
      if(mqmqder) then
         write(*,*)'3XQ 2215 Derivatives of vk_ij relative to quads',ncv
         do icv=1,ncv
            write(*,79)'ij',(mqf%compvar(icv)%dvk_ij(dd),dd=1,nqx)
            write(*,79)'ji',(mqf%compvar(icv)%dvk_ji(dd),dd=1,nqx)
79          format('3XQ dvk',a,10(1pe12.4))
         enddo
!         write(*,86)mqmqa_data%ncat,(mqf%y_ik(icv),icv=1,mqmqa_data%ncat)
86       format(/'3XQ y_i/k:',i2,3x,7F9.6)
      endif
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
      proprec=>endmemrec%propointer
      aff=one/mqmqa_data%pp(1,mqmqj)
!     write(*,'(a,3i3,1pe12.4)')'3XQ endmemloop1B: ',mqmqj,kend,nrealem,aff
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
!        write(*,'(a,2i3,2(1pe12.4))')'3XQ FNN endmember',mqmqj,kend,&
!                pyq,vals(1)
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
   if(tch.ge.3) write(*,203)phres%yfr
203 format('3XQ adding reference to SNN endmembers',20(F8.5))
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
         if(tch.ge.3) write(*,205)'3XQ FNN: qix, FNN, aff, pyq, fun, DG: ',&
              mqmqj,zp,aff,pyq,refg(zp,1),pyq*aff*refg(zp,1)
205      format(a,2i3,F8.5,2x,3(1pe12.4))
         if(tch.ge.3) &
          write(*,210)'3XQ FNN G:     ',mqmqj,mqmqj,pyq,aff,pyq,phres%gval(1,1)
210      format(a,2i3,2F8.5,1pe12.4,2x,6(1pe10.2))
      else
! this is an SNN with two or more AA/XX records
! For each SNN pair add the contribution to the FNN reference state
! %contyp(1..4,mqmqj) is index of FNN reference energy
         if(tch.ge.3) write(*,'(a,i3,1x,4i3,4F8.5)')'3XQ pp2: ',mqmqj,&
              (mqmqa_data%contyp(s1,mqmqj),s1=6,9),&
              (mqmqa_data%pp(s1,mqmqj),s1=1,4)
         snnloop: do s1=6,9
! zp is index to an FNN record, there can be 2 or 4 FNN records
            zp=mqmqa_data%contyp(s1,mqmqj)
            if(zp.eq.0) exit snnloop
! %pp(1..4,mqmqj) is stoichiometric factors for the pair
            aff=mqmqa_data%pp(s1-5,mqmqj)
            if(tch.ge.3) write(*,211)1,mqmqj,ipy,phres%dgval(1,mqmqj,ipy)
211         format('3XQ SNN dG/dy:',3i3,(1pe12.4))
            do itp=1,3
               phres%dgval(itp,mqmqj,ipy)=phres%dgval(itp,mqmqj,ipy)+&
                    aff*refg(zp,itp)
            enddo
            if(tch.ge.3) write(*,212)zp,mqmqj,ipy,phres%dgval(1,mqmqj,ipy),&
                 aff,aff*refg(zp,1)
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
299 format('3XQ endmember energy and entropy calculated, excess to be done')
!   goto 800
!---------------------------------------------------------------------
! code below needed for excess parameters ONLY, 
! all SNN and FNN endmembers already done
! NOTE some of them may not have a reference energy parameter
! This is to allocate csumx for handling quads with small fractions.
   isumx=0
! debug output of G for check of excess
   if(mqmqxcess) then
      write(*,288)(phres%gval(itp,1),itp=1,4)
288   format(/'3XQ line 1439 before excess:'/'G, dG/dT dG/dP d2G/dT2:',&
           4(1pe14.6))
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
! THIS IS NEW EXCESS MQMQA CODE
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
!      if(calc_alldvkij) then
! calculate partial derivatives of all vk_ij etc
!         call calc_newdvkij_values(phres,ceq)
!         calc_alldvkij=.FALSE.
!      endif
!
      noofex=noofex+1
      if(mqmqxcess) write(*,*)'3XQ excess with endmember constituent: ',&
           mqmqjy,ipy
!      write(*,316)mqmqjy,ipy
      call new_mqmqa_excess(lokph,intrec,mqmqjy,vals,dvals,d2vals,gz,ceq)
!      write(*,*)'3XQ back from new_mqmqa_excess',gx%bmperr
      if(gx%bmperr.ne.0) goto 1000
!------------- important --------------------
! vals, dvals and d2vals is the SUM OF ALL EXCESS parameters for this endmember
! gz is pointer to gtp_parcalc .... for phases with parameter permutations
!------------- important --------------------
      if(mqmqxcess) write(*,316)mqmqjy,ipy,vals(1)
316   format('3XQ endmember excess: ', 2i5,' calculated: ',1pe12.4)
! intrec is nullified inside new_mqmqa_excess
!      if(mqmqxcess) write(*,*)'3XQ back with excess from endmember ',mqmqjy,&
!           associated(endmemrec)
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
! dvals(1,jq=1,gz%nofc) set here
!      write(*,3181)gz%nofc,(dvals(1,jq),jq=1,gz%nofc)
3181   format('3XQ dex: ',i3,25(1pe12.4))
      do itp=1,6
! loop for G, G.T, G.P, G.T.T, G.T.P, G.P.P to add excess contribution
         phres%gval(itp,ipy)=phres%gval(itp,ipy)+vals(itp)
      enddo
! TEMPORARILY REMOVED SOME LOOPS
      if(mqmqder) write(*,3183)gz%nofc,(dvals(1,jq),jq=1,gz%nofc)
3183  format('3XQ dvals: ',i3,20(1pe12.4))
      do jq=1,gz%nofc
! skip loop for dG/dy, d2G/dydT, d2G/dydP, only for constituents
!         do itp=1,3   this skips 2nd derivative wrt T ??  ipy=1 is G
         if(mqmqder .and. abs(dvals(1,jq)).gt.1.0D-3) then
! this line is strange, supress it temporarily
            write(*,3182)jq,dvals(1,jq),phres%dgval(1,jq,1)
         endif
         do itp=1,3
            phres%dgval(itp,jq,ipy)=phres%dgval(itp,jq,ipy)+dvals(itp,jq)
         enddo
!         write(*,3182)jq,dvals(1,jq),phres%dgval(1,jq,1)
      enddo
3182  format('3XQ line 1567 addexcess: ',i3,2(1pe12.4))
!      
! ignore 2nd derivatives as not calculated for excess
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!
!      write(*,289)(phres%gval(itp,1),itp=1,4)
289   format('3XQ line 1532 after excess:'/'G, dG/dT dG/dP d2G/dT2:',4(1pe14.6))
      if(.not.associated(intrec)) cycle endmemloop2
!
!*************** remove all code below when excess code above OK *******
!******************* new code above should replace code below *******
!
! There are excess parameters, any Tooprecords?
!
! THIS IS OLD CODE WHISH SHOULD NO LONGER BE USED
!
!******************* new code above should replace code below *******
499      continue
      if(oldmqmqa_model) then
         write(*,319)
319      format(/'3XQ *** this is the old mqmqa excess model **'/)
!      stop "gtp3XQ line 1542"
      endif
!
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
         stop '3XQ aborting 1'
      endif
      if(associated(intrec%highlink)) then
! a higher interaction TO BE ADDED AND TESTED
         write(*,*)'3XQ ternary parameters not yet implemented?  Line 1841'
         stop '3XQ aborting 2'
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
      write(*,*)'3XQ phres pointer is not assigned entering calc_toop'
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
   if(x12.ge.one) then
      write(*,*)'3XQ Error: x12 larger than 1.0 in Toop/Kohler extrapolation!'
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
 subroutine new_mqmqa_excess(lokph,intrecin,mqmqj,vals,dvals,d2vals,gz,ceq)
! vals(1..6) are G, dG.T, dG.P, d2G.T.T, d2G.T.P and d2G.P.P for parameter
! dvals(1,i) are first derivatives wrt fracton and 2nd wrt fraction, T or P
!          dvals(1,i) is dG.yi, dval2(2,i) is d2G.yi.T, dvals(3,i) is d2G.yi.P
! d2vals(i,j) are second derivatives to 2 fractions, IGNORED HERE
! gz%nofc is number of fraction variables multiplied with this parameter(?)
!
! written using the gtp_allinone data structure for asymmetric excess
   implicit none
! mqmqj is index of first constituent in endmemberrecord
   integer lokph,mqmqj
!   type(gtp_property), pointer :: lokpty
   type(gtp_parcalc) :: gz
   type(gtp_phase_varres), pointer :: phres
   TYPE(gtp_mqmqa_var), pointer :: mqf
! dvals(1,x) is derivative wrt constituent x, dvals(2,x) is d2G/dTdx .....
! dvals(3,x) is derivative wrt d2G/dPdx
   double precision vals(6),dvals(3,gz%nofc)
! **** d2vals NOT USED
   double precision d2vals(gz%nofc*(gz%nofc+1)/2)
! pointer to first interaction record from an endmember
! intrecin is copied to intrec and then nullified. intrec may be updated below,
   TYPE(gtp_interaction), pointer :: intrecin,intrec,intrecfirst
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! needed locally?
   TYPE(gtp_intstack), dimension(:), allocatable :: savedint
   TYPE(gtp_pystack), pointer :: pystack
   TYPE(gtp_phase_add), pointer :: addrec
   TYPE(gtp_terdata), pointer :: ternaries
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_interaction), pointer :: ternaryexcess
!   type(gtp_allinone), pointer :: compvar
!
   character*120 text
   double precision :: rtg
   logical :: once=.true.
   integer ppow,qpow,rpow,intlev,iiz,jj,pairquad,jp
   integer, save :: proprecno=0
   integer parquad(4),nprr,nfr,dd
   integer :: nex=0
   integer, dimension(:), allocatable :: ylinks,qlinks
   integer ncv,icv,nqx,lokcs,lokfun,xq,cxq,mm
   integer termm
! include ternary and d(ternary)/dT
   double precision compprod,nomin,ternary(2),tpfun(6)
!
   logical, save :: ternaryonce=.true.
!
! composition derivatives are only relative to quads  !!!!!!!
! the composition variables for a parameters are expressions using asymmetric
! y_ik, \xi or \varkappa which depend on quads
! we have to sort out how this affects the derivatives
! dy_ik are factors for y_ik relative to quads, can be 1 or less
! dxi_ij and dx_ji and  dvk_ij and dvk_ji are 1 or less
! A parameter P multiplied with vk_ij(ij) has several contributions to the
! derivatives dP(zz), dvk_ij(ij,zz),zz=1,nquad
   integer idyix(5,mqmqa_data%nquad)
   integer zkij,nvkappa,ijx,nexrec,nooftps
   double precision term1,term2,dterm1,dterm2,dsum
   double precision haha,one1,dnomin,ddivisor,dternary
   double precision dyix(5,mqmqa_data%nquad),values(6)
   double precision dvalxq(mqmqa_data%nquad)
! pder(nquad) is the derivative of parameter wrt to each quad ij
!   double precision pder(mqmqa_data%nquad*(mqmqa_data%nquad-1)/2
   double precision debugder(mqmqa_data%nquad)  ! for debug only
! dvkijz(1,*) are dvk_ij/dxquad and dvkijz(2,*) are dvk_ji/dxquad
   double precision dvkijz(2,mqmqa_data%nquad)
! These are short for the values of vk_ij and vk_ji
   double precision vk_ij,vk_ji
! partial derivative for one parameter contribution
   double precision d1vals(mqmqa_data%nquad)
   double precision dtvals(mqmqa_data%nquad)
! FactSage Factor
!   double precision :: FSF=1.0d0
!
   integer jix,ii
   character*1 ptyp1
! The previous MQMQA excess implementation arrive here
! If mqmqa_data%exlevel is zero we should return and old code will still work.
!   write(*,*)'3XQ in new_mqmqa_excess',mqmqa_data%exlevel
   if(mqmqder) write(*,*)'3XQ in new_mqmqa_excess',mqmqa_data%exlevel
   if(mqmqa_data%exlevel.eq.0) then
!      if(once) write(*,6)mqmqa_data%exlevel
6     format('3XQ *** this system use the old excess model ***',i5)
      goto 1000
   endif
!--------------------------------------------------------------
! we are here because this endmember has an intercation link
!   if(mqmqxcess .and. associated(intrecin)) write(*,5)mqmqj
   if(mqmqxcess) write(*,5)mqmqj
!   write(*,5)mqmqj
5  format(/'3XQ in new_mqmqa_excess ',i3,' with intreraction record')
! initiate ylinks for this tree with the endmember fraction
   intrecfirst=>intrecin
   intrec=>intrecin
! this is needed to move to next endmemeber
   nullify(intrecin)
!
! below divide values with rtg?
   rtg=globaldata%rgas*ceq%tpval(1)
!---------------------------------
! not more than 10 interactions ....
   allocate(savedint(10))
   allocate(ylinks(10))
   allocate(qlinks(10))
! fraction indices
   ylinks=0
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
! THIS IS THE EXCESS CALCULATION ROUTINE WITH extensive DEBUG LISTING ADDED
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
! when we are here intrec must be associated
! All interaction records from a single endmember record will be calculated
   nexrec=0
   intlev=0
   nooftps=0
! if no parameters return zero
! otherwise vals and dvals sums up all excess parameters for this endmember
!   write(*,*)'3XQ line 2218: zero excess contribution'
!   mqmqx_deltag=0.0d0
! These values return sum of all excess parameters for this endmember
   vals=zero
! dvals has dimension dvals(3,mqmqa_data%nquad) dG/dy_i, d2G/dTdy_i, d2G/dPdy_i
   dvals=zero
! summing partial derivatives for each separate parameter, no P derivative
!   d1vals=zero
!   dtvals=zero
!
   intloop: do while(associated(intrec))
!
! intrec must be associated here, nexrec just counts interaction records
      nexrec=nexrec+1
!
! there is a single set of sites, save constituent first index
! We may come back here for another interaction with same endmember
! Set ylinks to be indices of the OC fractions x_11, x_12, ... x_nn
!      write(*,88)intlev,associated(intrec%propointer),intrec%fraclink
88    format('3XQ Starting intloop with component ',i3,l2,5i4)
! save name of interacting constituent even if no property
      nfr=nfr+1
      ylinks(nfr)=intrec%fraclink(1)
      proprec=>intrec%propointer
! termm is the index of a possible ternary cation, initiate to zero
      termm=0
! loop for all property records for same set of constituents
      proplist: do while(associated(proprec))
! we have found an excess parameter !!!
! This is used also to calculate the ternay parameter!!!
         lokfun=proprec%degreelink(0)
         if(lokfun.gt.0) nooftps=nooftps+1
         call eval_tpfun(lokfun,ceq%tpval,tpfun,ceq%eq_tpres)
         if(gx%bmperr.ne.0) goto 1000
         if(mqmqxcess) then
            write(*,113)ptyp1,lokfun,rtg,tpfun(1),tpfun(2)
113         format('3XQ tpfun: ',a,i4,6(1pe12.4))
         endif
         if(tpfun(1).eq.0.0d0) then
! skip if there is no TP function (no TPFUN used during testing)
            proprec=>proprec%nextpr
            nex=nex+1
            cycle proplist
         endif
! divide all parameter values with rtg!!
         tpfun=tpfun/rtg
! calculate all d(varkappa_ij)/dx_kl
         nvkappa=size(mqf%compvar)
! there can be several property record for the same set of constituents
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
!         if(.true.) then
! emergency debug LIST PARAMETER to understand what the parameter it is ....
            jp=1
            text=' '
            call mqmqa_excesspar_name(lokph,intlev,nfr,ylinks,text,jp)
            text(jp-1:)=';'//ptyp1//','//char(ichar('0')+ppow)//&
                 ','//char(ichar('0')+qpow)//','//char(ichar('0')+rpow)//')'
            write(*,115)trim(text),ppow,qpow,rpow
115      format(/'3XQ param: ',a,', pqr:',3i2)
! extract the quad pointers
         endif
!
!         lokfun=proprec%degreelink(0)
! can xq be zero here ??????
         xq=proprec%asymdata%quad
! cxq transforms the quad index to an index in compvar (which as not diagonal)
!
! here vk_ij are used to multiply with the parameter
! where are they calculated?
         cxq=mqmqa_data%quad2compvar(xq)
         vk_ij=mqf%compvar(cxq)%vk_ij
         vk_ji=mqf%compvar(cxq)%vk_ji
!         write(*,1160)xq,mqf%xquad(xq),vk_ij,vk_ji,mqf%compvar(cxq)%denominator
1160     format('3XQ xq etc: ',i3,4(1pe14.6))
!------------------------------------------------------------- ternary
         ternary(1)=one
         ternary(2)=zero
         terparam: if(nfr.gt.3) then
!            if(ternaryonce) write(*,116)
!            write(*,116)
116         format('3XQ line 2416 found ternary parameter')
!            ternaryonce=.false.
!            goto 1000
! locating argument for the ternary factor, is there any asymmetry?
! the mqf pointer is:   mqf=>ceq%phase_varres(lokcs)%mqmqaf, lokcs is the mqmqa
!            mqmqa_data%ncat
!            write(*,117)(mqf%y_ik(jix),jix=1,mqmqa_data%ncat)
117         format('3XQ line 2487 y_i/k',10(1pe12.4))
!            write(*,118)cxq,mqf%compvar(cxq)%xi_ij,mqf%compvar(cxq)%xi_ji
118         format('3XQ line 2489: ',i3,' xi_ij: ',1pe12.4,' xi_ji: ',1pe12.4)
! this is a call which currently does not calculate anything
!            write(*,*)'3XQ calling ternary_factor',size(ternary),termm
!            write(*,77)lokfun,tpfun(1)*rtg
77          format('3XQ tpfun before ternary_factor1 ',i3,1pf15.6)
            call ternary_factor1(lokph,mqf,xq,cxq,&
                 ylinks,termm,ternary,proprec,ppow,qpow,rpow,ceq)
            if(gx%bmperr.ne.0) goto 1000
! calculated ternary factor
!            write(*,130)ternary
130         format('3XQ back from ternary: :',2(1pe14.6))
! list parameter values here are divided by R*T, multiply to see parameter
            write(*,114)nfr,ptyp1,lokfun,rtg,rtg*tpfun(1),rtg*tpfun(2)
114         format('3XQ line 2466 ternary: ',i2,2x,a,i4,3(1pe12.4))
         else
!
! This should be a binary parameter (with 3 composition variables)
!         write(*,20)(ylinks(ii),&
!              trim(splista(phlista(lokph)%constitlist(ylinks(ii)))%symbol),&
!              ii=1,3),ppow,qpow,lokfun
20       format('3XQ3 L(PH',3(',',i1,':',a),') pows:',2i2,' fun: ',i3)

!
         endif terparam
!------------------------------------------------------------- end ternary
! Maybe a scaling difference with FactSage, multiply tpfun by FSF
!         FSF=1.5D0
!         write(*,*)'3XQ Scaling with ',FSF
!         tpfun(1)=FSF*tpfun(1)
!
!--------------------------------------------------------------------
! multiply the parameter with the composition variables
         ptyp: if(ptyp1.eq.'G') then
! ppow is for varkappa_ij, qpow is for varkappa_ji, term1 and term2 used below
! vk_ij and vk_ji are (sum of quands)/(sum of quads)
            term1=1.0d0
            term2=1.0d0
! if ppow or qpow is zero the term is unity
            if(ppow.gt.0) term1=vk_ij**ppow
            if(qpow.gt.0) term2=vk_ji**qpow
            nomin=term1*term2
! default ternary = 1.00, rtg=R*T
            compprod=mqf%xquad(xq)*nomin*ternary(1)
! vals(1) is the sum of all excess parameters linked from this endmember
            vals(1)=vals(1)+compprod*tpfun(1)
! What about derivatives wrt T?
            vals(2)=vals(2)+mqf%xquad(xq)*nomin*ternary(2)
! list 2 indices, 2 powers, 3 constitutions, tpfun, constituents*tpfun, vals
            if(xq.gt.0) then
! list value of excess parameter
               if(mqmqxcess) then
                  write(*,991)xq,nexrec,ppow,qpow,&
!                    mqf%xquad(xq),vk_ij,vk_ji,&
                       mqf%xquad(xq),term1,term2,&
                       rtg*tpfun(1),compprod*tpfun(1),vals(1)
991               format('3XQ line 2357:',i3,3i2,3F7.4,3(1pE12.4))
               endif
            else
               write(*,*)'3XQ no quad index!'
               stop
            endif
!--------------------------------------------------------------------
! BEGIN calculate partial derivatives ...........
! any quad can be involved in compvar(cxq)%vk_ij
! For ternary parameters some additional derivatives may be needed
            if(mqmqder) then
               write(*,992)xq,cxq,mqf%compvar(cxq)%cat1,mqf%compvar(cxq)%cat2
992            format('3XQ derivatives of quad: ',i2,', and vk_ij and vk_ji: ',&
                    i3,2x,2i3)
               write(*,*)'3XQ calling dvkij_dzijkl for varkappa: ',cxq
            endif
            nqx=mqmqa_data%nquad
            ncv=size(mqf%compvar)
! 
! cxq is varkappa involved with this parameter
! calculate all partial derivatives of this wrt nqx quad fractions
! The vk_ij/vk_ji are used for several parameters and their derivatives
! should calculated only once
! dvkijz(1, 1..nqx) are derivatives of vk_ij: dvk_ij/dxz
! dvkijz(2, 1..nqx) are derivatives of vk_ji: dvk_ji/dxz
            call dvkij_dzijk(mqf,cxq,dvkijz)
            if(gx%bmperr.ne.0) goto 1000
!
! loop for derivatives of parameter for all quads
            zkijloop: do zkij=1,nqx
!
! EG = xq * (vk_ij**pp) * (vk_ji**qq) * param
!
! dEG/dxz = xq * pp*(vk_ij**(pp-1))*dvk_ij/dxz * (vk_ji**qq) * param +
!           xq * (vk_ij**pp) * qq*(vk_ji**(qq-1))*dvk_ji/dxz * param +
!           dxq/dxz * (vk_ij**pp) * (vk_ji**qq) * param 
!
               if(ppow.eq.0) then
                  dterm1=term2
               elseif(ppow.eq.1) then
                  dterm1=dvkijz(1,zkij) * term2
               else
                  dterm1=ppow*vk_ij**(ppow-1)*dvkijz(1,zkij)*term2
               endif
! dvkijz(1,zkij) is dvk_ij/dxz and 
! dvkijz(2,zkij) is dvk_ji/dxz
               if(qpow.eq.0) then
                  dterm2=term1
               elseif(qpow.eq.1) then
                  dterm2=term1*dvkijz(2,zkij)
               else
                  dterm2=term1*qpow*vk_ji**(qpow-1)*dvkijz(2,zkij)
               endif
               if(zkij.eq.xq) then
                  dsum=(dterm1+dterm2)*mqf%xquad(xq)+term1*term2
               else
                  dsum=(dterm1+dterm2)*mqf%xquad(xq)
               endif
! dvkijz are the derivative of EG with respect to xqz
! dvals(1,...) is dG/dy, 
               dvals(1,zkij)=dvals(1,zkij)+dsum*tpfun(1)
! dvals(2,...) is d2G/dydT, dvals(3,...) is d2G/dydP
               dvals(2,zkij)=dvals(2,zkij)+dsum*tpfun(2)
! this is just for debug output below
               debugder(zkij)=dsum
            enddo zkijloop
! debug output of all fraction product derivatives
! rtg&tpfun(1) and vals(1) listed at line 2357
!            write(*,997)rtg*tpfun(1),vals(1),(debugder(zkij),zkij=1,nqx)
!            write(*,997)(debugder(zkij),zkij=1,nqx)
997         format('3XQ df/dx:',20(1pe11.3))
! partial derivative end >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if(mqmqder) write(*,122)(dvals(1,zkij),zkij=1,nqx)
122         format('3XQ dEG/dxz ',20(1pe11.3))
         elseif(ptyp1.eq.'Q') then
!
            write(*,*)'3XQ the Q parameter not implemented yet',ptyp1
            stop
         else
            write(*,*)'3XQ the B parameter not implemented: ',ptyp1
            stop
         endif ptyp
         if(mqmqxcess) write(*,*)'3XQ end of a parameter: ',ptyp1,lokfun,nexrec
!
800      continue
         if(mqmqxcess .and. associated(proprec)) write(*,96)' more ',vals(1)
96       format('3XQ ',a,' Current Excess G: ',1pe12.4)
         proprec=>proprec%nextpr
         nex=nex+1
      enddo proplist
!
! we have calculated all property records for one excess parameter
! there can be a ternary or more binary parameters
!
      ternaryexcess=>intrec%highlink
! All mqmqa parameters are "ternary" or higher
!      write(*,811)associated(ternaryexcess)
!811   format('3XQ is there a link to higher excess?',l2)
!
!         terrec=>intrec%highlink
!         call calc_ternarymq(lokph,phres,terrec,vals,dvals,d2vals,gz,ceq)
!      
!      if(associated(ternaryexcess)) then
! The ternary excess parameter handled above terparam: if .... endif terparam
!         write(*,*)'3XQ this must be an error, not implemeneted'
!         nullify(ternaryexcess)
!
!      endif
      push_ornext: if(associated(intrec%highlink)) then
! go to  higher level of interaction but save link to next for other parameters
         intlev=intlev+1
         if(associated(intrec%nextlink)) then
            if(mqmqxcess) write(*,97)intlev,intrec%nextlink%fraclink
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
!         write(*,98)associated(intrec)
!98       format('3XQ wow, a ternary parameter? ',l2)
!         stop '3XQ line 2622 not implemented yet'
! for ternary parameters we just cycle intloop with one more constituent
      else
         if(mqmqxcess) write(*,*)'3XQ any more excess on level?',intlev,nexrec
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
         if(.not.associated(intrec)) exit intloop
! why cycle?
!         cycle intloop
      endif push_ornext
   enddo intloop
!
!------------ return to next endmember
1000  continue
   if(mqmqxcess) then
      if(associated(intrecfirst)) then
         proprec=>intrecfirst%propointer
         write(*,1001)nexrec,vals(1) ! ,(dvals(1,mm),mm=1,mqmqa_data%nquad)
1001     format('3XQ exit new_mqmqa_excess, excess records: ',i5,2x,1pe12.4)
      endif
   endif
!   write(*,1099)vals(1),nexrec,nooftps
!   if(mqmqxcess) write(*,1099)vals(1),nexrec,nooftps
1099 format('3XQ exit new_mqmqa_excess with G=',1pe12.4,&
          ', ',i3,' parameters and ',i3,' TPFUNs')
   return
 end subroutine new_mqmqa_excess

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine dvkij_dzijk
!\begin{verbatim}
 subroutine dvkij_dzijk(mqf,cxq,dvkijk)
! calculates all partial derivatives of a parameter multiplied with
!          xquad * varkappa**ppow * varkappa**qpow
   implicit none
!   type(gtp_phase_varres), pointer :: phres
   type(gtp_mqmqa_var), pointer :: mqf
   type(gtp_allinone), pointer :: box
   integer cxq
! there are mqmqa_data%nquad variable and derivatives 
   double precision dvkijk(2,mqmqa_data%nquad)
! cxq is the varkappa index
! dvkijk is the the 2D array with derivatives of vk_ij and vk_ji with respect
! to all quad fractions.  Many of them will be zero
!\end{verbatim}
!
! Looking for errors in ternaries, suspect missing derivative wrt
! derivatives of quad fraction in denominator _kvk not included !!!???
! df/dx = -nominator/denominator**2 = -g/h**2
!
   integer ijkl,vkix,vkdenom,kk,mxq,dgij,dgji,dijk,cat1,cat2
   double precision sumi, sumj, sumk, dvkij, dvkji, dvkdenom
   logical skip,dksum
!   integer, allocatable :: qdone(:)
!   integer, allocatable :: indenom
   integer denomx,qq,ii
!
! derivative of a quotient  d(g/h) = 1/h*dg/dx - (g/h**2)*dh/dx =
!                                     (h*dg/dx - g*dh/dx)/h**2
! dgij=0 if no derivative of vkij
! dgji=0 if no derivative of vkji
! dijk=0 if no derivative of kvkall ?
! NOTE both vkij and vkji has the same denominator, specified by kvkijk !!!
!
! This routine calculates both d(vk_ij)/dx and d(vk_ji)/dx   
! with respect to all quadrupole fractions
! 
   dvkijk=0.0d0
!
   mxq=mqmqa_data%nquad
   box=>mqf%compvar(cxq)
   denomx=size(box%all_ijk)
   cat1=box%cat1
   cat2=box%cat2
! initiate done with all quad indices
!   allocate(qdone(denomx))
!   qdone=mqf%compvar(cxq)%all_ijk
!   write(*,7)size(dvkijk),qdone
7  format('3XQ *** enter dvkij_dzijk, qdone: ',i5,2x,20i3)
!   write(*,8)dvkijk
8  format('3XQ dvkijk:',6(1pe10.2))
!
! set all partical derivatives to zero as default return
   if(mqmqder) then
      write(*,*)'3XQ *** entering dvkij_dzijkl',cxq,mxq
!      write(*,*)'3XQ in dvkij_dzijkl',cxq,mxq
!      write(*,10)mqf%compvar(cxq)%ivk_ij
!      write(*,20)mqf%compvar(cxq)%jvk_ji
!      write(*,30)mqf%compvar(cxq)%kvk_ijk
!      write(*,30)mqf%compvar(cxq)%all_ijk
      write(*,10)box%ivk_ij
      write(*,20)box%jvk_ji
      write(*,30)box%kvk_ijk
      write(*,40)box%all_ijk
10    format('3XQ ivk_ij: ',10i3)
20    format('3XQ jvk_ji: ',10i3)
30    format('3XQ kvk_ijk:',10i3)
40    format('3XQ kvk_all:',10i3)
   endif
! initiate all partial derivaties to zero
   dksum=.true.
   sumi=zero
   sumj=zero
!
!   vk_ij = \sum_i xquad(ivk_ij) / (\sum_k xquad(kvk_ijk)+\sum_k xquad(ivk_ij))
!   vk_ji = \sum_i xquad(jvk_ji) / (\sum_k xquad(kvk_ijk)+\sum_k xquad(jvk_ji))
!
! This routine calculates derivatives of variables vk_ij does not depend on!
!
!   quadloop: do ijkl=1,mxq
! The loop for ksum needed only once ....................
!      do_dksum: if(dksum) then
! if derivative of sum_k xquad(kvk_ijk), if zero all derivatives zero
! the sumk include sum of all fractions in \sum_i and \sum_j
!         sumk=zero
! derivative of denominator is zero or one
!         dgij=0
!         dgji=0
!         dijk=0
!         skip=.true.
!
   sumk=0.0d0
   denominator: do kk=1,size(box%all_ijk)
! a quad fraction term can only apper once
      vkix=box%all_ijk(kk)
      sumk=sumk+mqf%xquad(vkix)
! listing of derivative calculations
! 3XQ kloop for mqf%compvar( 1)%all_ijk( 1)    1 sum  3.5933E-02 1 1 1
! 3XQ kloop for mqf%compvar( 1)%all_ijk( 2)    1 sum  6.7187E-01 1 3 1
! 3XQ kloop for mqf%compvar( 1)%all_ijk( 3)    1 sum  1.0000E+00 1 2 1
! 3XQ iloop for mqf%compvar( 1)%ivk_ij( 1)     1 sum  3.5933E-02 1 1 1
! 3XQ jloop for mqf%compvar( 1)%jvk_ji( 1)     1 sum  6.3593E-01 0 3 1
! 3XQ dvk:   1 1 0  3.5933E-02  6.3593E-01  1.0000E+00    9.6407E-01 -6.3593E-01
! 3XQ kloop for mqf%compvar( 1)%all_ijk( 1)    2 sum  3.5933E-02 1 1 1
!                                   1   2      3          4  
!      if(mqmqder) write(*,50)'k',cxq, ')%all_ijk(', kk,') ',&
!           ijkl,'k',sumk,1,vkix,ijkl
!                 5    6   7   8  9    10
!
!                         1                          2  3   4 
!50          format('3XQ ',a,'loop for mqf%compvar(',i2, a, i2,a,&
!                 i4,' sum',a,': ',1pe12.4,2x,3i2)
!                 5        6       7      8-10
!  mqf%compvar(cxq)%all_ijk(kk),ijkl,  sumk         
!
!               skip=.false.
!               write(*,*)'3XQ ********* ',kk,vkix,ijkl
!               dijk=1
!            endif
   enddo denominator
!   write(*,51)size(box%all_ijk),sumk,box%all_ijk
51 format('3XQ Summed ',i2,' quads in kvk_ijk ',1pe12.4,10i3)
   
! if sumk is zero there are no derivatives with respect to this quad
! this is not an error, just a message
!         write(*,54)cxq,ijkl
!54       format('3XQ vk(',i2,') does not depend on quad ',i2)
!
! below only if the vk_ij and vk_ji do not depend on ijkl
!
! The loop above summed all fraction variables of denominator, needed below
! We have to take care of the nominators if varkappa_ij and varkappa_ji      
!      sumi=zero
!
   sumi=0.0d0
   nominator1: do kk=1,size(box%ivk_ij)
!      dgij=0
      vkix=box%ivk_ij(kk)
      sumi=sumi+mqf%xquad(vkix)
!         if(mqf%compvar(cxq)%ivk_ij(kk).eq.ijkl) then 
!         if(box%ivk_ij(kk).eq.ijkl) then
!         if(vkix.eq.ijkl) then
!            dgij=1
!         endif
!                                1   2      3          4  
!      if(mqmqder) write(*,50)'i',cxq, ')%vk_ij(', kk,')   ',&
!           ijkl,'i',sumi,1,vkix,ijkl
!              5    6   7   8  9    10
   enddo nominator1
!   write(*,55)size(box%ivk_ij),sumi,box%ivk_ij
55 format('3XQ Summed ',i2,' quads in ivk_ij',1pe12.4,10i3)
!
! note sumi and/or sumj can be zero
!
   sumj=0.0d0
   nominator2: do kk=1,size(box%jvk_ji)
!      dgji=0
      vkix=box%jvk_ji(kk)
      sumj=sumj+mqf%xquad(vkix)
!         if(mqf%compvar(cxq)%jvk_ji(kk).eq.ijkl) then
!         if(box%jvk_ji(kk).eq.ijkl) then
!         if(vkix.eq.ijkl) then
!            dgji=1
!         endif
!                                1   2      3          4  
!         if(mqmqder) write(*,50)'j',cxq, ')%vk_ji(', kk,')   ',&
!              ijkl,'j',sumj,1,vkix,ijkl
!              5    6   7   8  9    10
   enddo nominator2
!   write(*,56)size(box%jvk_ji),sumj,box%jvk_ji
56 format('3XQ Summed ',i2,' quads in jvk_ji',1pe12.4,10i3)
!
! derivative of a quotient  d(g/h) = (1/h)*dg/dx - (g/h**2)*dh/dx = 
!              = (h*dg/dx - g*dh/dx)/h**2
!
! the derivarive value g=sumi/sumk; h=sumj/sumk;  di and dj can be zero  
! the derivarive value dj*sumi-di*sumj
! 
   if(sumk.eq.zero) then
      write(*,*)'3XQ line 2691, division by zero, check source code!!!'
      sumk=1.0d0
   endif
! Attempt 2026.03.29 to fix problem derivatives wrt fractions in denomonator
!
! the derivatives are calculated here, dg/dx and dh/dx is 0 or 1
! derivative of a quotient  d(g/h) = 1/h*(dg/dx) - (g/h**2)*dh/dx
!
! the denominatorof a vk_ij contains all quads
   dijk=1.0d0
   loopdenom: do kk=1,size(box%all_ijk)
! all quads in the vk are present in the denominator
      ijkl=box%all_ijk(kk)
      checknom1: do ii=1,size(box%ivk_ij)
         if(box%ivk_ij(ii).eq.ijkl) then
            dvkijk(1,ijkl)=(sumk - sumi*dijk)/sumk**2
         else
            dvkijk(1,ijkl)= -sumi*dijk/sumk**2
         endif
      enddo checknom1
      checknom2: do ii=1,size(box%jvk_ji)
! do not use dgij and dgji are 0 or 1 depending on quads in vk_ij or vk_ji
         if(box%jvk_ji(ii).eq.ijkl) then
            dvkijk(2,ijkl)=(sumk - sumj*dijk)/sumk**2
         else
            dvkijk(2,ijkl)= - sumj*dijk/sumk**2
         endif
      enddo checknom2
!
!      write(*,69)cxq,ijkl,dvkijk(1,ijkl),dvkijk(2,ijkl)
69    format('3XQ dvk(',i2,')_ij&_ji/dq(',i2,') = ',2(1pe12.4))
   enddo loopdenom
   if(mqmqder) then
      write(*,70)ijkl,sumi,sumj,sumk,&
           dvkijk(1,ijkl),dvkijk(2,ijkl)
70    format('3XQ dvk2: ',i2,3(1pe12.4),1x,2(1pe12.4))
   endif
!---------------------------------------------------------
! return derivatives of dvarkappa(ceq)/dxquad for all quads
!
1000 continue
   if(mqmqder) then
      write(*,1090)cat1,cat2,(dvkijk(1,ijkl),ijkl=1,mxq)
      write(*,1090)cat2,cat1,(dvkijk(2,ijkl),ijkl=1,mxq)
1090    format('3XQ dvk(',i1,',',i1,')/dq: ',20(1pe10.2))
      write(*,*)'3XQ *** exit dvkij_dzijk'
   endif
   return
 end subroutine dvkij_dzijk

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine mqmqa_excesspar_name
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

!\addtotable subroutine ternary_factor1
!\begin{verbatim}
! subroutine ternary_factor1(lokph,mqf,xq,cat1,cat2,ylinks,termm,hejhopp,&
!      proprec,ppow,qpow,rpow)
 subroutine ternary_factor1(lokph,mqf,xq,cxq,ylinks,termm,hejhopp,proprec,&
      ppow,qpow,rpow,ceq)
! calculates the ternary factor of a parameter
   implicit none
   integer lokph,xq,cxq,termm,ylinks(*),ppow,qpow,rpow
! return parameter and d(parameter)/dT  NOT d(parameter)/d(fracs) ??
! xq is the index of the quad2compvar record
! cxq is the index of the compvar record with vk_ij and vk_ji
! ylinks is array with quadfractions (OC fraction variables)
! termm is the index if the y_ik constituent variable (set below)
! hejhopp returns the calculated parameter and its T derivative   
! proprec is the property record with the TPFUN expression and other data
! ppow, qpow and tpow are the powers of the constituent variables
! ceq is the global data pointer ... 
   double precision hejhopp(2)
   type(gtp_property), pointer :: proprec
   TYPE(gtp_mqmqa_var), pointer :: mqf
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
!
! This deals with the composition dependence of the ternary factor
! it calls ternary_factor2 for the actual parameter
!
! xq is the AB/X quad index
! cxq in varkappa index which gives 2 quad indices for A/X and B/X
! ylinks are the OC fraction indices, not really used in this routine
! termm is index of the unknown 4th quad where is it?
! hejhopp is the value to return, possibly 1.0D0
! This subroutine calculates the factor in Max Poschmann eq.25 and eq.26
! If the ternary constituent "termm" is asymmetric for i-j-\gamma
!
!   Y_m/k             Y_j/k
! (---------) ( 1 -  --------- )**(r-1)       if i is asymmetric in i-j-m ?
!   xi_ji/k           xi_ji/k
!
!
!   Y_m/k             Y_i/k
! (---------) ( 1 -  --------- )**(r-1)       if j is asymmetric in i-j-m ?
!   xi_ij/k           xi_ij/k
!
!
!   Y_m/k ( 1 - xi_ij/k - xi_ji/k )**(r-1)    if neither case above
!
! This factor (and its first derivatives) must be multiplied with the binary
!    parameter expression already calculated in the subroutine calling this
!
! In 26.09.17 I come back to this code and I have learned a bit more,  Y_mk
! are in the TYPE gtp_mqmqa_var usually accessed by the pointer mqf%
! and they are called y_ik and their values are in mqf%y_ik(1..n)
! but they are not indexed by ylinks which is the constituent fractions
! The xi_ij variables are part of the compvar%xi_ij data structure
! related to the varkappa vk_ij and vk_ji (with asymmetries)
! Probably a new subroutine is needed above ternary_factor to handle fractions
!
!---------------------------------------------------------------------------
! Test case Ce-Cl-Li-Mg   THIS RESULT IS WRONG
!
!OC version  6.129  equilibrium:   1, DEFAULT_EQUILIBRIUM        2026.09.20
!Conditions .................................................:
!  1:T=1000, 2:P=100000, 3:N(LI)=1, 4:N(CE)-N(LI)=0, 5:N(MG)=1, 6:AC(CL)=1
! Degrees of freedom are   0
!
!Some global data, reference state SER ......................:
!T=   1000.00 K (   726.85 C), P=  1.0000E+05 Pa, V=  0.0000E+00 m3
!N=   9.0000E+00 moles, B=   3.8408E+02 g, RT=   8.3145E+03 J/mol
!G= -2.56821E+06 J, G/N=-2.8536E+05 J/mol, H=-1.8260E+06 J, S= 7.422E+02 J/K
!
!Some data for components ...................................:
!Component name    Moles      Mole-fr  Chem.pot/RT  Activities  Ref.state
!CE                1.0000E+00  0.11111 -1.5333E+02  2.5785E-67  SER (default)   
!CL                6.0000E+00  0.66667  0.0000E+00  1.0000E+00  SER (default)   
!LI                1.0000E+00  0.11111 -6.1282E+01  2.4288E-27  SER (default)   
!MG                1.0000E+00  0.11111 -9.4274E+01  1.1409E-41  SER (default)   
!
!Some data for phases .......................................:
!Name                Status Moles      Volume    Form.Units Cmp/FU dGm/RT  Comp:
!MSCL.................... E  9.000E+00  0.00E+00  6.70E+00    1.34  0.00E+00  X:
! CL     6.66667E-01  MG     1.11111E-01  CE     1.11111E-01  LI     1.11111E-01
!Constitution: There are     6 constituents:
! CEMG/CL-Q03  3.54502E-01  CE/CL-Q01    1.78078E-01  MG/CL-Q06    1.12093E-01
! CELI/CL-Q02  1.84988E-01  LIMG/CL-Q05  1.39707E-01  LI/CL-Q04    3.06330E-02
!
! no asymmetry and ternary parameter included a binary
!-------------------------------------------------------------------------
! reult below not correct
!   
!OC version  6.130  equilibrium:   1, DEFAULT_EQUILIBRIUM        2026.09.20
!Conditions .................................................:
!  1:T=1000, 2:P=100000, 3:N(LI)=1, 4:N(CE)-N(LI)=0, 5:N(MG)=1, 6:AC(CL)=1
! Degrees of freedom are   0
!
!Some global data, reference state SER ......................:
!T=   1000.00 K (   726.85 C), P=  1.0000E+05 Pa, V=  0.0000E+00 m3
!N=   9.0000E+00 moles, B=   3.8408E+02 g, RT=   8.3145E+03 J/mol
!G= -2.56049E+06 J, G/N=-2.8450E+05 J/mol, H=-4.3861E+05 J, S= 2.122E+03 J/K
!
!Some data for components ...................................:
!Component name    Moles      Mole-fr  Chem.pot/RT  Activities  Ref.state
!CE                1.0000E+00  0.11111 -1.5299E+02  3.6234E-67  SER (default)   
!CL                6.0000E+00  0.66667  0.0000E+00  1.0000E+00  SER (default)   
!LI                1.0000E+00  0.11111 -6.1066E+01  3.0148E-27  SER (default)   
!MG                1.0000E+00  0.11111 -9.3903E+01  1.6539E-41  SER (default)   
!
!Some data for phases .......................................:
!Name                Status Moles      Volume    Form.Units Cmp/FU dGm/RT  Comp:
!MSCL.................... E  9.000E+00  0.00E+00  6.68E+00    1.35  0.00E+00  X:
! CL     6.66667E-01  CE     1.11111E-01  LI     1.11111E-01  MG     1.11111E-01
!Constitution: There are     6 constituents:
! CEMG/CL-Q03  3.59555E-01  CE/CL-Q01    1.75071E-01  MG/CL-Q06    1.10525E-01
! CELI/CL-Q02  1.88671E-01  LIMG/CL-Q05  1.37986E-01  LI/CL-Q04    2.81917E-02
!
! with modifications but WITHOUT TERNARY PARAMETER  ........... suck
! BUT with the  TERNARY PARAMETER  included as binay........... suck
!--------------------------------------------------------------------------
!
! Now the ternary perameters hould be toally ignored
!
!OC version  6.130  equilibrium:   1, DEFAULT_EQUILIBRIUM        2026.09.20
!Conditions .................................................:
!  1:T=1000, 2:P=100000, 3:N(LI)=1, 4:N(CE)-N(LI)=0, 5:N(MG)=1, 6:AC(CL)=1
! Degrees of freedom are   0
!
!Some global data, reference state SER ......................:
!T=   1000.00 K (   726.85 C), P=  1.0000E+05 Pa, V=  0.0000E+00 m3
!N=   9.0000E+00 moles, B=   3.8408E+02 g, RT=   8.3145E+03 J/mol
!G= -2.56005E+06 J, G/N=-2.8445E+05 J/mol, H= 7.7218E+05 J, S= 3.332E+03 J/K
!
!Some data for components ...................................:
!Component name    Moles      Mole-fr  Chem.pot/RT  Activities  Ref.state
!CE                1.0000E+00  0.11111 -1.5297E+02  3.6950E-67  SER (default)   
!CL                6.0000E+00  0.66667  0.0000E+00  1.0000E+00  SER (default)   
!LI                1.0000E+00  0.11111 -6.1054E+01  3.0522E-27  SER (default)   
!MG                1.0000E+00  0.11111 -9.3881E+01  1.6896E-41  SER (default)   
!
!Some data for phases .......................................:
!Name                Status Moles      Volume    Form.Units Cmp/FU dGm/RT  Comp:
!MSCL.................... E  9.000E+00  0.00E+00  6.68E+00    1.35  0.00E+00  X:
! CL     6.66667E-01  MG     1.11111E-01  CE     1.11111E-01  LI     1.11111E-01
!Constitution: There are     6 constituents:
! CEMG/CL-Q03  3.59845E-01  CE/CL-Q01    1.74899E-01  MG/CL-Q06    1.10434E-01
! CELI/CL-Q02  1.88879E-01  LIMG/CL-Q05  1.37885E-01  LI/CL-Q04    2.80575E-02
!
! Hopefully correct result when totally ignoring the ternary parameter
!
!--------------------------------------------------------------------------
! Looking at values of G is almost the same but H and S varies a lot STRANGE
!--------------------------------------------------------------------------
!
! It seems frequently that r=1, i.e. The terma (...)**(r-1) can be ignored
! which means the ternary parameter depend only the Y_m/k term
! (multiplied with the x_ij (NOT xi_ij) i.e. almost a binary parameter ...
! Although Y_m/k may depend on several x_ij fractions.
!
! The Y_m/k values are in gtp_mqmqa_var, by pointer mqf
! 
   integer ii,jj,cat1,cat2,cat3,t1,toop
   integer,save :: noter=0
   integer, save :: lastupdate=0
!
   double precision tpfun(6),rtg,asymp
   integer lokfun,nooftps
   integer low,middle,high
!
! we must update xi_ij, xi_ji and v_ik ONCE for each iteration
! It updates all xi_ij etc so called only once per iteration
! mqmqa_pairupdate is in gtp3_dd2
!   write(*,2)lastupdate,mqmqa_data%mqmqa_terasym1
2  format('3XQ in ternary_factor',2i4)
!
   if(lastupdate.eq.mqmqa_data%mqmqa_terasym1) then
! this should be fixed some time ....
      write(*,*)'3XQ skipping ternary_factor 1'
      goto 1010
   endif
!   
   cat1=mqf%compvar(cxq)%cat1
   cat2=mqf%compvar(cxq)%cat2
! ylinks are OC phase constituent indices, not necessarily same as quad indices
!
! mqmqa_data%emquad(1..n) are indices of quads with single caton: A/X, B/X etc.
! The quadruplet index for the ternary quad is in ylinks(4)
!   do ii=1,size(mqmqa_data%emquad)
!      if(ylinks(4).eq.mqmqa_data%emquad(ii)) goto 12
!   enddo
!   write(*,*)'3XQ cannot find the ternary C/X quad'
!   gx%bmperr=4399; goto 1000
!
! return the index of the ternary cation in the parameter ... needed for y_ik
!12 termm=ii
!
! At the first call here the variables proprec%tersysix and proprec%dat3
! are initiated and used for future calls
!
   if(proprec%asymdata%tersysix.eq.0) then
!      write(*,7)
7     format('3XQ Ternary parameter asymmetry initiatiated, only once')
! check if this ternary is asymmetric, default is -1 meaning no Toop
!      proprec%asymdata%tersysix=-1
! we can have ylinks(1..4), one of them is a mixed quad, 
! two is equal to cat1 and cat2, we need the third one!
      yloop: do jj=1,4
         if(ylinks(jj).eq.cat1) cycle yloop
         emloop: do ii=1,size(mqmqa_data%emquad)
            if(ylinks(jj).eq.mqmqa_data%emquad(ii)) goto 12
         enddo emloop
      enddo yloop
      write(*,*)'3XQ cannot find the ternary C/X quad'
      gx%bmperr=4399; goto 1000
!
! save the index of the third cation in the C/X quad
12    continue
      cat3=ii
! save this for later calls
      proprec%asymdata%cat3=ii
! termm is the variable used initially below, do not mess things up!!!
      termm=ii
! then we must check if the ternary has cat1 or cat2 as Toop
!      write(*,15)cat1,cat2,cat3
15    format('3XQ ternary cations: ',3i3)
! n*(n-1)*(n-3)/6  5*4*3/6=10
! 1 2 3; 1 2 4; 1 2 5; 1 3 4; 1 3 5; 1 4 5; 2 3 4; 2 3 5; 2 4 5; 3 4 5
! certainly cat1 < cat2 but cat3 can be lower, higher or in between 
      if(cat3.lt.cat1) then
         low=cat3
         middle=cat1
         high=cat2
      elseif(cat3.lt.cat2) then
         low=cat1
         middle=cat3
         high=cat2
      else
         low=cat1
         middle=cat2
         high=cat3
      endif
      allter: do t1=1,size(tersys)
         if(tersys(t1)%el(1).ne.low) cycle allter
         if(tersys(t1)%el(2).ne.middle) cycle allter
         if(tersys(t1)%el(3).ne.high) cycle allter
! when we arrive here tersys(t1) has the cations cat1, cat2 and cat3
! But we need to save this only if cat1 or cat2 is Toop
         exit allter
      enddo allter
! here we have the t1 array with cat1, cat2 and cat3
      if(t1.le.size(tersys)) then
!         write(*,*)'3XQ found ternary',t1
         toop=index(tersys(t1)%asymm,'T')
         if(toop.le.0) then
            proprec%asymdata%tersysix=-1
            proprec%asymdata%toop=0
         elseif(tersys(t1)%el(toop).eq.cat1) then
! In cat1-cat2-\gamma we have cat1 as Toop, use xi_ji (still to be added)
            write(*,*)'3XQ found asymmetric ternary',t1
            proprec%asymdata%tersysix=t1
            proprec%asymdata%toop=cat1
         elseif(tersys(t1)%el(toop).eq.cat2) then
! In cat1-cat2-\nu  we have cat2 as Toop, use xi_ij (still to be added)
            write(*,*)'3XQ found asymmetric ternary',t1
            proprec%asymdata%tersysix=t1
            proprec%asymdata%toop=cat2
         endif
      end if
!      write(*,17)'3XQ initiated ',proprec%asymdata%cat3,&
!           proprec%asymdata%tersysix,proprec%asymdata%toop
17    format(a,' system with m=',i2,' and ternary',i3,' and Toop ',i2)
!   else
! The second and later calls just use the saved values
!      write(*,17)'3XQ using ',proprec%asymdata%cat3,proprec%asymdata%tersysix,&
!           proprec%asymdata%toop
   endif
! we must have the ternary cation index and if cat1 or cat2 is toop
   termm=proprec%asymdata%cat3
   toop=proprec%asymdata%toop
!------------------------------------------------------
!   write(*,13)termm,ylinks(4),&
!        (mqmqa_data%emquad(ii),ii=1,size(mqmqa_data%emquad))
13 format('3XQ line 3043 ternary constituent: ',2i2,5x,10i3)
!
!   if(noter.eq.0) then
!      write(*,*)'3XQ ternary parameters not implemented'
!      noter=1
!   else
      noter=noter+1
      if(noter.eq.1) then
! debug just the indices
!         write(*,10)ii,xq,cat1,cat2, termm,size(mqmqa_data%emquad),&
!              (ylinks(jj),jj=1,3),ylinks(4),associated(proprec)
10       format('3XQ line 2987 ternary: ',i3,5x,3i3,3x,2i3,3x,4i3,3x,L)
! debug list the constituents also
! ylinks is phase constituent index, 
         lokfun=proprec%degreelink(0)
!         write(*,20)(ylinks(ii),&
!              trim(splista(phlista(lokph)%constitlist(ylinks(ii)))%symbol),&
!              ii=1,4),termm,ppow,qpow,rpow,lokfun
!20       format('3XQ L(',4(i2,': ',a,1x),') m:',i1,', pows: ',3i1)
20       format('3XQ4 L(PH',4(',',i1,':',a),') m:',i1,' pows:',3i2,i4)
         noter=0
!         write(*,30)(phlista(lokph)%constitlist(ylinks(ii)),ii=1,4)
30       format('3XQ Redundant indices to x_ii fractions: ',4i3)
! here are y_i/k and xi_ij ??
!         write(*,40)mqf%y_ik
!         write(*,50)mqf%compvar(cxq)%xi_ij,mqf%compvar(cxq)%xi_ji
40       format('3XQ y_ik: ',10f10.6)
50       format('3XQ x_ij: ',f10.6,', x_ji: ',f10.6)
! current values of all fractions involved:
! parameter value, there is just a single one .... I hope
         nooftps=0
!         lokfun=proprec%degreelink(0)
         call ternary_factor2(lokfun,mqf,xq,cxq,ylinks,termm,hejhopp,proprec)
!         write(*,55)lokfun
55       format('3XQ lokfun: ',i5)
         if(lokfun.gt.0) nooftps=nooftps+1
!         write(*,*)'3XQ ceq%tpfun: ',ceq%tpval
!         write(*,*)'3XQ tpfun: ',tpfun
!         write(*,*)'3XQ ceq%eq_tpres: ',ceq%eq_tpres
! 6 values returned, L, dL/dT, dL/dP, d2L/dT2, d2L/dP2, d2L/dTdP
! maybe eval_tpfun is used in calling routine and here only fractions ???
         if(lokfun.le.0) then
            tpfun=0.0d0
         else
!            write(*,*)'3XQ calling eval_tpfun',lokfun
!                           int     TP       result  pointer
! wow, in ceq%eq_tpres save all calculated values with eval_tpfun ....!!!
            call eval_tpfun(lokfun,ceq%tpval,tpfun,ceq%eq_tpres)
            if(gx%bmperr.ne.0) goto 1000
         endif
!         write(*,77)lokfun,tpfun(1)
77       format('3XQ tpfun in ternary_factor1 ',i3,1pf15.6)
!
!         write(*,*)'3XQ back from tpfun'
!         write(*,100)termm,mqf%y_ik(termm),&
!              mqf%compvar(cxq)%xi_ij,mqf%compvar(cxq)%xi_ji,&
!              tpfun(1)
100      format('3XQ m and y_m: ',i2,f10.6,', xi_ij, xi_ji: ',2f10.6,1pe15.4)
         rtg=globaldata%rgas*ceq%tpval(1)
         tpfun=tpfun/rtg
!         write(*,110)tpfun(1),tpfun(2)
110      format('3XQ L/RT and (dL/dT)/RT ',2(1pe14.6))
!   Y_m/k ( 1 - xi_ij/k - xi_ji/k )**(r-1)    if neither case above
         if(rpow.gt.1) then
            asymp=(1.0d0-mqf%compvar(cxq)%xi_ij-mqf%compvar(cxq)%xi_ji)**rpow
         else
            asymp=1.0d0
         endif
         hejhopp=mqf%y_ik(termm)*asymp*tpfun(1)
      elseif(noter.gt.10) then
         noter=1
      endif
!   endif
!
1000 continue
!      write(*,*)'3XQ No forced update of ternary factor ...?'
!     write(*,*)'3XQ Forced ternary_factor updates: ',mqmqa_data%mqmqa_terasym1
!     lastupdate=mqmqa_data%mqmqa_terasym1
1010 continue
!      write(*,*)'3XQ leaving ternary_factor',termm
!      write(*,*)'3XQ ternary code not updated: mqmqa_data%mqmqa_terasym1',&
!           mqmqa_data%mqmqa_terasym1
   return
 end subroutine ternary_factor1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine ternary_factor2
!\begin{verbatim}
 subroutine ternary_factor2(lokph,mqf,xq,cxq,ylinks,termm,hejhopp,proprec)
! calculates the ternary factor of a parameter
   implicit none
   integer lokph,xq,cxq,termm,ylinks(*),ppow,qpow,rpow
! return parameter and d(parameter)/dT
   double precision hejhopp(2)
   type(gtp_property), pointer :: proprec
   TYPE(gtp_mqmqa_var), pointer :: mqf
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
!
! This deals with the parameter value, not with the composition dependance
!
! xq is the AB/X quad index
! cxq in varkappa index which gives 2 quad indices for A/X and B/X
! ylinks are the OC fraction indices, not really used in this routine
! termm is index of the unknown 4th quad
! hejhopp is the value to return, possibly 1.0D0
! This subroutine calculates the factor in Max Poschmann eq.25 and eq.26
! If the ternary constituent "termm" is asymmetric for i-j-\gamma
!
!   Y_m/k             Y_j/k
! (---------) ( 1 -  --------- )**(r-1)       if i is asymmetric in i-j-m ?
!   xi_ji/k           xi_ji/k
!
!
!   Y_m/k             Y_i/k
! (---------) ( 1 -  --------- )**(r-1)       if j is asymmetric in i-j-m ?
!   xi_ij/k           xi_ij/k
!
!
!   Y_m/k ( 1 - xi_ij/k - xi_ji/k )**(r-1)    if neither case above
!
! This factor (and its first derivatives) must be multiplied with the binary
!    parameter expression already calculated in the subroutine calling this
!
! In 26.09.17 I come back to this code and I have learned a bit more,  Y_mk
! are in the TYPE gtp_mqmqa_var usually accessed by the pointer mqf%
! and they are called y_ik and their values are in mqf%y_ik(1..n)
! but they are not indexed by ylinks which is the constituent fractions
! The xi_ij variables are part of the compvar%xi_ij data structure
! related to the varkappa vk_ij and vk_ji (with asymmetries)
! Probably a new subroutine is needed above ternary_factor to handle fractions
!
! It seems frequently that r=1, i.e. The terma (...)**(r-1) can be ignored
! which means the ternary parameter depend only the Y_m/k term
! (multiplied with the x_ij (NOT xi_ij) i.e. almost a binary parameter ...
! Although Y_m/k may depend on several x_ij fractions.
!
! The Y_m/k values are in gtp_mqmqa_var, by pointer mqf
! 
   integer ii,jj,cat1,cat2
   integer,save :: noter=0
   integer, save :: lastupdate=0
!
   double precision tpfun(6),rtg,asymp
   integer lokfun,nooftps
! we must update xi_ij, xi_ji and v_ik ONCE for each iteration
! It updates all xi_ij etc so called only once per iteration
! mqmqa_pairupdate is in gtp3_dd2
!   write(*,*)'3XQ Unfinished ternary_factor2'
   goto 1000
!
! code below not finished, has to be rewritten completely
!
! ********* BE CAREFUL ternary_factor2 NOT YET IMPLEMENTED *********
!
   if(lastupdate.eq.mqmqa_data%mqmqa_terasym1) then
      goto 1010
   endif
!   
   cat1=mqf%compvar(cxq)%cat1
   cat2=mqf%compvar(cxq)%cat2
! ylinks are OC phase constituent indices, not necessarily same as quad indices
! ylinks(4) is the ternary quad, transform to y_i/k index
! 
!   do ii=1,size(mqmqa_data%emquad)
!      if(ylinks(3).eq.mqmqa_data%emquad(ii)) goto 12
!   enddo
!   write(*,*)'3XQ line 2980 ternary constituent: ',termm
!   if(noter.eq.0) then
!      write(*,*)'3XQ ternary parameters not implemented'
!      noter=1
!   else
      noter=noter+1
      if(noter.eq.1) then
! debug just the indices
!         write(*,10)ii,xq,cat1,cat2, termm,size(mqmqa_data%emquad),&
!              (ylinks(jj),jj=1,3),ylinks(4),associated(proprec)
10       format('3XQ line 2987 ternary: ',i3,5x,3i3,3x,2i3,3x,4i3,3x,L)
! debug list the constituents also
! ylinks is phase constituent index, 
!         write(*,20)(ylinks(ii),&
!              trim(splista(phlista(lokph)%constitlist(ylinks(ii)))%symbol),&
!              ii=1,4),termm,ppow,qpow,rpow
!20       format('3XQ L(',4(i2,': ',a,1x),') m:',i1,' pows:',3i2)
         noter=0
!         write(*,30)(phlista(lokph)%constitlist(ylinks(ii)),ii=1,4)
30       format('3XQ Redundant indices to x_ii fractions: ',4i3)
! here are y_i/k and xi_ij ??
!         write(*,40)mqf%y_ik
!         write(*,50)mqf%compvar(cxq)%xi_ij,mqf%compvar(cxq)%xi_ji
40       format('3XQ y_ik: ',10f10.6)
50       format('3XQ x_ij: ',f10.6,', x_ji: ',f10.6)
! current values of all fractions involved:
! parameter value, there is just a single one .... I hope
         nooftps=0
         lokfun=proprec%degreelink(0)
!         write(*,55)lokfun
55       format('3XQ lokfun: ',i5)
         if(lokfun.gt.0) nooftps=nooftps+1
!         write(*,*)'3XQ ceq%tpfun: ',ceq%tpval
!         write(*,*)'3XQ tpfun: ',tpfun
!         write(*,*)'3XQ ceq%eq_tpres: ',ceq%eq_tpres
! 6 values returned, L, dL/dT, dL/dP, d2L/dT2, d2L/dP2, d2L/dTdP
! maybe eval_tpfun is used in calling routine and here only fractions ???
         if(lokfun.le.0) then
            tpfun=0.0d0
         else
!            write(*,*)'3XQ calling eval_tpfun',lokfun
!                           int     TP       result  pointer
! wow, in ceq%eq_tpres save all calculated values with eval_tpfun ....!!!
            call eval_tpfun(lokfun,ceq%tpval,tpfun,ceq%eq_tpres)
            if(gx%bmperr.ne.0) goto 1000
         endif
!
!         write(*,*)'3XQ back from tpfun'
!         write(*,100)termm,mqf%y_ik(termm),&
!              mqf%compvar(cxq)%xi_ij,mqf%compvar(cxq)%xi_ji,&
!              tpfun(1)
100      format('3XQ m and y_m: ',i2,f10.6,', xi_ij, xi_ji: ',2f10.6,1pe15.4)
         rtg=globaldata%rgas*ceq%tpval(1)
         tpfun=tpfun/rtg
!         write(*,110)tpfun(1),tpfun(2)
110      format('3XQ L/RT and (dL/dT)/RT ',2(1pe14.6))
!   Y_m/k ( 1 - xi_ij/k - xi_ji/k )**(r-1)    if neither case above
         if(rpow.gt.1) then
            asymp=(1.0d0-mqf%compvar(cxq)%xi_ij-mqf%compvar(cxq)%xi_ji)**rpow
         else
            asymp=1.0d0
         endif
         hejhopp=mqf%y_ik(termm)*asymp*tpfun(1)
      elseif(noter.gt.10) then
         noter=1
      endif
!   endif
!
1000 continue
!      write(*,*)'3XQ No forced update of ternary factor ...?'
!     write(*,*)'3XQ Forced ternary_factor updates: ',mqmqa_data%mqmqa_terasym1
!     lastupdate=mqmqa_data%mqmqa_terasym1
1010 continue
!      write(*,*)'3XQ leaving ternary_factor',termm
!      write(*,*)'3XQ ternary code not updated: mqmqa_data%mqmqa_terasym1',&
!           mqmqa_data%mqmqa_terasym1
   return
 end subroutine ternary_factor2

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
30 format('3XQ we found the pair: ',i3,', lowa:',i3,', qorder:',15i3)
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
            write(*,*)'3XQ line 3261 circles are square'
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
        ' Latt som en platt')
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
!      write(*,*)'3XQ line 3262 U species: ',ia,iq,size(mqmqaf%xquad)
      mqmqaf%xquad(ia)=phres%yfr(iq)
      if(verbose) write(*,26)ia,phres%yfr(ia),iq,mqmqaf%xquad(iq)
26    format('3XQ the OC fraction: ',i3,1pe14.6,&
           ' is set to MQMQA quad: ',i3,1pe14.6)
   enddo
!  if(verbose) write(*,*)'3XQ calling calcasymvar for \varkappa_ij, \xi_ij etc.'
!   call calcasymvar(phres,0)
! I think second argument nonzero means NOT INITIATE
   call calcasymvar(phres,1)
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
   integer ic,ia,lokph,jix,jiy
!   type(gtp_phase_record), pointer :: phase
   type(gtp_ternary_asymmetry), pointer :: asym3rec
! there is a global mqmqa_data record to use!! <<<<<<<<<<<<<<<<<<,
!\end{verbatim}
   integer i,j,k,nseq,mm,apos,nbinsys,ntercat,nn,binx
!   integer i,j,k,nseq,mm,apos,lcat,lnan,nbinsys,ntercat
! how to create xquad mm when we need a pointer to gtp_phase_varres?
   type(gtp_equilibrium_data), pointer :: ceq
   type(gtp_phase_varres), pointer :: phres
   type(gtp_mqmqa_var), pointer :: mqf
   character*6 defasym,ch1*1
! Many properties are symmetric, for example xquad which has a single index
! and is indexed by ijkl(i,j,k,l) where ijkl(i,j,k,l)=ijkl(j,i,k,l)
! but other are unsymmetric such as varkappa and xi
!
   if(mqverbose) write(*,*)'3XQ In init_excess_asymm',lokph
!
   ceq=>firsteq
! I have forgotten how OC works.  When entering phases one can creat
! data structures in gtp_equilibrium_data (record pointer ceq)
! and these will be copied when new equilibrium records created
! (for example parallel calculations).  When a second gtp_equilibium_data
! has been created one iMUST NOT change this data_structure
! the line below creates a pointer to the gtp_mqmqa_var data inside ceq
! maybe problem with the array here ...
   i=1
! we loop back here from a few lines below
5  continue
      i=i+1
! this is very clumsy, but I have no better way
      phres=>ceq%phase_varres(i)
!      write(*,*)'loop: ',i,lokph,phres%phlink
      if(phres%phlink.ne.lokph) goto 5
!----------------------------------------------------------------
! we have found the mscl phase
!      write(*,*)'Found phase_varres for the mqmqa liquid phase',i
! set the pointer to the mqmqaf record
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
! we have to initiate con2quad below with the corresponding cation indices
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
! frequently there will be a single anion
!      nquad=ncat*(ncat+1)/2
!   endif
!!
! varkappa and xi_ijis now part of allinone
66 continue
!
!   if(mqmqa_data%ncat.eq.2) goto 80
!   write(*,67)mqmqa_data%ncat*(mqmqa_data%ncat-1)*mqmqa_data%nan/2
67 format('3XQ init_excess_asymm allocating asymmetrical compvar array: ',i5)
! we have to intitiate several variables in each compvar
!   allocate(compvar(ncat*(ncat-1)*nan/2))
!
! HERE THE COMPVAR ARRAY IS CREATED line 3439, later called BOX
   allocate(mqf%compvar(mqmqa_data%ncat*(mqmqa_data%ncat-1)/2*mqmqa_data%nan))
!   write(*,*)'3XQ, initiating compvar for excess model variables',&
!        mqmqa_data%ncat,size(mqf%compvar)
   if(allocated(mqmqa_data%el2ancat)) then
!      write(*,69)
69    format('3XQ Heureca! el2ancat allocated')
!      write(*,70)size(mqmqa_data%el2ancat),mqmqa_data%ncat,mqmqa_data%el2ancat
70    format('3XQ el2ancat: ',2i3,5x,20i3)
   else 
      write(*,4)
4     format('3XQ line 3663: The array mqmqa_data%el2ancat not allocated!'/&
           'Should have been done in correlate_const_and_quads')
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
!-------------------------- moved initiation of y_ik before compvar
!   write(*,90)mqmqa_data%ncat*mqmqa_data%nan
90 format('Allocating pair fraction array y_i/k: ',i4)
   if(mqmqa_data%nan.ne.1) then
      write(*,91)
91    format(/'3XQ Sorry can only handle systems with 1 anion'/)
      stop mqmqa_data%nan
   endif
! y_ik varies with the current constitution
   allocate(mqf%y_ik(mqmqa_data%ncat*mqmqa_data%nan))
! for understanding the code later on, also save names of cations
!   write(*,*)'3XQ adding some information to the mqf data structure'
! data structures in gtp_mqmqa
!   write(*,'(a,10i3)')'3XQ line 3667 gtp_mqmqa%cat2species: ',&
!        mqmqa_data%cat2species
! emquad is constituent index of endmember species
!   write(*,'(a,10i3)')'3XQ line 3669 allocating names_y_ik'
   allocate(mqf%names_y_ik(mqmqa_data%ncat))
   allocate(mqf%spix_y_ik(mqmqa_data%ncat))
   allocate(mqf%spqx_y_ik(mqmqa_data%ncat))
! needed for asymmetries <<<<<<<<<<<<<<< this is where it is allocated
   allocate(mqmqa_data%el2quad(mqmqa_data%ncat))
!   write(*,*)'3XQ **** ALLOCATED EL2QUAD **** ',size(mqmqa_data%el2quad)
!   allocate(mqf%names_y_ik(mqmqa_data%ncat*mqmqa_data%nan))
! we want the name of constituent emquads for the y_ik
!   write(*,*)'3XQ line 3748 looking for y_ik species names',&
!        lokph,size(mqf%names_y_ik)
! list of constituent species are in phlista(lokph)%constitlist
!   write(*,*)'3XQ number of quad fractions: ',size(phlista(lokph)%constitlist)
!
! the name of y_ik should be the element names of emquads with a single cation
!      mqf%names_y_ik(jix)='cat_0'//char(ichar('0')+jix)
! this cannot be correct, we must fetch them from the phase record
   jiy=0
   nn=0
   do i=1,mqmqa_data%ncat
      do j=i,mqmqa_data%ncat
         nn=nn+1
         mm=phlista(lokph)%constitlist(nn)
         if(j.eq.i) then
! only save names when i=j         
            jiy=jiy+1
            mqf%names_y_ik(jiy)=splista(mm)%symbol
! remove numbers trailing -Q
            k=index(mqf%names_y_ik(jiy),'-Q')
            if(k.le.0) then
               write(*,*)'3XQ illegal quadname: ',mqf%names_y_ik(jiy)
            else
               mqf%names_y_ik(jiy)(k+2:)=' '
            endif
! Hm, mm is not correct, probably the order they were added
! maybe splista(mm)%alphaorder ?
!            write(*,82)nn,mm,splista(mm)%alphaindex,splista(mm)%symbol
82          format('3XQ line 3572 basic cation: ',3i3,' ',a)
            mqf%spix_y_ik(jiy)=splista(mm)%alphaindex
!            write(*,*)'3XQ y_ik index',mm,splista(mm)%alphaindex,nn
!            mqf%spix_y_ik(jiy)=splista(mm)
            mqf%spix_y_ik(jiy)=splista(mm)%alphaindex
!  mm is species index, nn is constituent index
            mqf%spqx_y_ik(jiy)=nn
! needed for asymmetry !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            mqmqa_data%el2quad(jiy)=nn
         endif
      enddo
   enddo
! check
!-------------------------------------------
! now initate record with asymmetries
!   write(*,*)'Allocating asymmetries',mqmqa_data%ncat
   if(allocated(tersys)) then
      write(*,*)'3XQ tersys already allocated'
      goto 66
   endif
   nseq=0
   binx=1
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
! tersys(mm)%emquads are the indices in mqf%spqx_y_ik(i) ....?
                     if(allocated(mqf%spqx_y_ik)) then
                        tersys(mm)%emquad(1)=mqf%spqx_y_ik(i)
                        tersys(mm)%emquad(2)=mqf%spqx_y_ik(j)
                        tersys(mm)%emquad(3)=mqf%spqx_y_ik(k)
                     else
                        write(*,*)'3XQ sorry mqf%spqx_y_ik)) not allocated'
                     endif
! set symmetric
                     tersys(mm)%asymm='KKK'
                     tersys(mm)%isasym=0
                     tersys(mm)%noasym=0
! index of the 3 binary systems associated with this ternary
                     tersys(mm)%binsys(1)=ibin(i,j,mqmqa_data%ncat)
                     tersys(mm)%binsys(2)=ibin(i,k,mqmqa_data%ncat)
                     tersys(mm)%binsys(3)=ibin(j,k,mqmqa_data%ncat)
!                     write(*,177)mm,tersys(mm)%el,tersys(mm)%binsys
177                  format('3XQ tersys : ',i3,' %el: ',3i3,', %bin: ',3i3)
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
! ------------- copy the species indices to compvar(
!   do nn=1,size(mqf%names_y_ik)
!   write(*,95)size(mqf%names_y_ik),&
!        (trim(mqf%names_y_ik(nn)),nn=1,size(mqf%names_y_ik))
95    format('3XQ ',i2,' quadruplets ii: ',10(a,', '))
!   enddo
!
! with multiple anion derivatives add dimension nan also
! its content is set in varkappa1
!
!   write(*,*)'3XQ line 3743 calling pairfracs'
!   call pairfracs(.true.,phres)
! the .true. forces some listing inside pairfracs
   call pairfracs(.false.,phres)
!
!-------------------------- moved initiation of y_ik above
   nseq=0
!   write(*,*)'3XQ initiating compvar in line 3524:',size(mqf%compvar)
   first: do i=1,mqmqa_data%ncat-1
      second: do j=i+1,mqmqa_data%ncat
! initiallize allinone record, allocated as compvar array, called box later
         nseq=nseq+1
         mqf%compvar(nseq)%seq=nseq
! these indices are from 1 to n-1 ignoring anions and actual element/species
         mqf%compvar(nseq)%cat1=i
         mqf%compvar(nseq)%cat2=j
! these are the element indices in OC
! These should be cation species as set above!!
!         if(i.gt.mqmqa_data%xanionalpha) mqf%compvar(nseq)%elcat1=i+1
!         if(j.gt.mqmqa_data%xanionalpha) mqf%compvar(nseq)%elcat2=j+1
! these are species indices of cations
!         mqf%compvar(nseq)%elcat1=mqf%spix_y_ik(i)
!         mqf%compvar(nseq)%elcat2=mqf%spix_y_ik(j)
         mqf%compvar(nseq)%quadicat1=mqf%spix_y_ik(i)
         mqf%compvar(nseq)%quadicat2=mqf%spix_y_ik(j)
! these are used when determining asymmetries
         mqf%compvar(nseq)%elcat1=mqf%spqx_y_ik(i)
         mqf%compvar(nseq)%elcat2=mqf%spqx_y_ik(j)
! note it is negative of element alphabetical index
         mqf%compvar(nseq)%elan=-mqmqa_data%xanionalpha
         mqf%compvar(nseq)%anion=1
         mqf%compvar(nseq)%boxlastupdate=-1
! should be cation/anon species names
!         write(*,65)i,j,mqf%names_y_ik(i),mqf%names_y_ik(j)
65       format('3XQ line 3616 cation: ',2i3,2x,a,2x,a)
         mqf%compvar(nseq)%qcat1=mqf%names_y_ik(i)
         mqf%compvar(nseq)%qcat2=mqf%names_y_ik(j)
         mqf%compvar(nseq)%qan='X  '
! ivk_ij, jvi_ji, kvk_ijk, xi_ij etc allocated in init_varkappa_symmetry
!         or when asymmetric in loop_tersys
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
! The arrays dxi_ij are allocated here and assigned below
         mqf%compvar(nseq)%dxi_ij=0.0d0
         mqf%compvar(nseq)%dxi_ji=0.0d0
!         write(*,77)nseq,i,j
77       format(i4,2i5)
      enddo second
   enddo first
80 continue

! initiate newXupdate, there is a newupdate I do not know where it is declared
   newXupdate=0
!
!888 continue
! Add code to handle asymmetries for the MQMQA phase if there is any
! characters in the mqmqa_data%tdbasymmetries
! this is called in set_ternary_asymmetry_tdb
!   if(allocated(mqmqa_data%tdbasymmetries)) then
!      write(*,*)'3XQ line 3674 there are asymmetries!',&
!           len_trim(mqmqa_data%tdbasymmetries)
!      call set_tdbasymmetries(phres)
!   endif
!
1000 continue
!
! REMOVE ncat from global data structure
!  write(*,99)mqmqa_data%ncat,mqmqa_data%nquad,size(phres%yfr),size(mqf%compvar)
99 format('3XQ leaving init_excess_asym: ',4i4/)
!
   return
 end subroutine init_excess_asymm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine pairfracs
!\begin{verbatim}
 subroutine pairfracs(list,phres)
! calculate all pair fractions from a set of quad fractions
! if there is a single anion
   implicit none
   logical list
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   integer i,j,v,dd,seq
   type(gtp_mqmqa_var), pointer :: mqf
! mqf is a pointer!!
!
   seq=0
   mqf=>phres%mqmqaf
!
   if(.not.allocated(mqf%xquad)) then
      write(*,*)'xquad not allocated'
      stop
   endif
   if(list) write(*,5)mqmqa_data%ncat,1,mqmqa_data%nquad
5  format('3XQ In pairfracts ncat=',i2,' and nan=',i2,' and nquad=',i3)
!
! For Kohler model:
! y_ik(i)  is 0.5*x_ii + 0.5\sum_j 0.5* x_ij
! When asymmetric the y_i/k can have mode terms, adding a second y_\nu/k
! Following Max paper eq.11 (if 2 or more anions eq.20 should be used)
!
! Calculate y_ik based on dy_ik which are arrays of constants  ncat*(ncat-1)/2
!         quad pairs:   11  12  13  14  22  23  24  33  34  44
!       quad indices:   1   2   3   4   5   6   7   8   9   10  ijklx(i,j,1,1)
! dy_ik(1:1..nquad) is  1.0 0.5 0.5 0.5 -   -   -   -   -   -
! dy_ik(2:1..nquad) is  -   0.5 -   -   1.0 0.5 0.5 0   -   -
! dy_ik(3:1..nquad) is  -   -   0.5 0   -   0.5 -   1   0.5 -
! dy_ik(4:1..nquad) is  -   -   -   0.5 -   -   0.5 -   0.5 1
!
! With asymmetries y_ik = \sum y_ik + y_i\nu  
!
   if(mqmqa_data%nan.ne.1) then
      write(*,*)'Cannot calculate pair fractions more than a single anions'
      stop
   endif
! write(*,*)'3XQ allocating mqmqa_data%%dy_ik: ',ncat,nan,nquad, assume nan=1
! dy_ik is a structure information, independent of current constitution
! This failes when mqmqa_data%nan > 1, check above that nan=1
   if(.not.allocated(mqmqa_data%dy_ik)) then
!    write(*,*)'3XQ line 3859 already allocated dy_ik',size(mqmqa_data%dy_ik),&
!           mqmqa_data%ncat*mqmqa_data%nquad
!   else
!      write(*,*)'3XQ line 3872 allocating dy_ik'
      allocate(mqmqa_data%dy_ik(mqmqa_data%ncat,mqmqa_data%nquad))
   endif
   mqmqa_data%dy_ik=0.0d0
!
   catloop1: do i=1,mqmqa_data%ncat
! loop will count each quad once including 11, 22 etc. Loops work for i=1
      v=0
      catloop2: do j=i,mqmqa_data%ncat
! the function ijklx(i,j,1,1) returns index of quad 
         v=v+1
         seq=ijklx(i,j,1,1)
         if(i.eq.j) then
            mqmqa_data%dy_ik(i,seq)=1.0D0
         elseif(j.gt.i) then
! this works !!!
            mqmqa_data%dy_ik(i,seq)=0.5D0
            mqmqa_data%dy_ik(j,seq)=0.5D0
         endif
      enddo catloop2
   enddo catloop1
!
! a lot of trouble finding this
   if(list) then
      write(*,6)
6     format(/'3XQ line 3843 Listing pair fractions'/&
           'dy_ik  quad: 1   2   3   4   5   6')
      do i=1,mqmqa_data%ncat
         write(*,12)i,(mqmqa_data%dy_ik(i,dd),dd=1,mqmqa_data%nquad)
12       format('dy_ik(',i1,',*): ',20F4.1)
      enddo
   endif
! now fix the xi_ij expressions as functions of x_ij assuming Kohler
!   write(*,*)'3XQ line 3858 initiating xi_ij variables.',&
!        size(mqf%compvar),size(mqf%y_ik)
!
! TYPE gtp_mqmqa_var has declations for y_ik
! gtp_allinone has declarations for vk_ij and xi_ij 
! mqf    is pointer to gtp_phase_var
! parres is pointer to gtp_phase_varres
!   phres=>mqf%phresq

   if(list) write(*,*)'3XQ exiting pairfracs'
1000 continue
   return
!
 end subroutine pairfracs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable function ibin(i,j,n)
!\begin{verbatim}
 integer function ibin(i,j,n)
! i and j are cations, there are n cations
! to find the indices of binaries a-b, a-c and c-b for ternary a-b-c
! when there are n elements, a<b<c
! iab=ibin(a,b,n)
! iac=iab+(c-b)
! ibc=ibin(b,c,n)
! provided by AI Claude ....
   implicit none
   integer i,j,n
!\end{verbatim}
   integer lo,hi
   lo=min(i,j); hi=max(i,j)
   ibin=(lo-1)*(2*n-lo)/2+(hi-lo)
 end function ibin
   
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
        j.le.0 .or. j.gt.mqmqa_data%ncon1) then
      write(*,7)'cation',i,j,mqmqa_data%ncon1
7     format('3XQ in ijklx wrong ',a,' indices: ',3i4)
      goto 2000
   endif
   if(k.le.0 .or. k.gt.mqmqa_data%ncon2 .or. &
        l.le.0 .or. l.gt.mqmqa_data%ncon2) then
      write(*,7)'anion',l,k,mqmqa_data%ncon1
      goto 2000
   endif
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
20 format('3XQ Cation index in ijklx: ',2i3,2i10)
   iquad=kquad*mqmqa_data%lcat+iquad
!   write(*,30)iquad,kquad,lcat,i,j,k,l
30 format('3XQ Index in xquad: ',i5,5x,3i5,5x,2i5)
   if(iquad.gt.mqmqa_data%nconst) goto 1000
   ijklx=iquad
!   write(*,*)'Return from ijklx with:',iquad
!
77 continue
   return
! errors 
1000 write(*,1010)i,j,k,l,mqmqa_data%ncon1,mqmqa_data%ncon2,mqmqa_data%lcat,&
          kquad,iquad
1010 format(' *** Indexing error in ijklx',4i4,2x,7i5,/'Stop!!!!')
   gx%bmperr=4399
   goto 77
!  
2000 continue
   write(*,2010)i,j,k,l,mqmqa_data%ncon1,mqmqa_data%ncon2
2010 format('3XQ In ijklx quad indices outside limits',4i3,5x,2i3)
   gx%bmperr=4399
   goto 77
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
!  4   5   -   -   -  20 (4,5,6) is the last ternary, 6*5*4/(2*3)=20
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
1100 write(*,*)'3XQ Indexing error in terind ',i,j,v
   iz=-1
   goto 1000
 end function terind

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine new_ternary_asym(asymter,new_toop,phres)
!\begin{verbatim}
 subroutine new_ternary_asym(asymter,new_toop,phres,xverbose)
   implicit none
   integer asymter,new_toop
   type(gtp_phase_varres), pointer :: phres
   logical xverbose
!\end{verbatim}
! This symbroutine modifies the vk_ij and xi_ij for a ternary with a new Toop.
! asymter is the index of the ternary a-b-c in tersys and new_toop
! is the cation which is Toop constituent.  The 2 binaries with the Toop
! cation must have their vk_ij modified.  
! IMPORTANT: existing Toop elements for these binaries are taken into account
! We have to find the 2 compvar for these binaries with existing aymmetries
! THERE IS NO WAY TO HANDLE REMOVING AN ASYMMETRY
   integer acat1,acat2,acat3,bin12, bin13, bin23,t1
   type(gtp_mqmqa_var), pointer :: mqf
! tersys is a global variable
!   type(gtp_terdata), pointer :: tersys
!
! Declartion of variables from varkappa1
!   integer mii,mij,mjj,ia
! this is the compvar record ... box should not be used
   type(gtp_allinone), pointer :: boxij,boxji,box
   logical verbose
!
   integer i,ii,vz,v,w,vv,ternary,ll,lasthope,di,icat,jcat,nnn,seq
   double precision varkappaij,varkappaji,sum,initialij,initialji,nugamma
   double precision xi_ij,xi_ji,sum1,sum2
   character*3 asymmetric
   integer ia
!
   integer, dimension(:), allocatable :: mixnugamma
   integer binij,binji,cat1,cat2
   character*3 asymmetry
! mixed update
   integer j,k,l,m,ny,abrakadabra
! If a binary i-j is part of 2 or more asymmetric ternaries i-j-\nu, i-j-\gamma
! the quad fraction x_\nu\gamma should be added to kvk_ijk (the denomonator)
! of kvk_ijk
! debug output
   integer nn1,nn2,nn3,nn4,nn5,nn6,nn7,gg,thisasym
   integer q1,q2,q3,qtoop,haha,nn
   logical nysym
   integer, allocatable, dimension(:) :: el2quadx
!
!
! debug
   verbose=xverbose
   if(verbose) write(*,*)'3XQ in new_ternary_asymmetry',asymter,new_toop
   asymmetric='KKK'
   if(new_toop.lt.1 .or. new_toop.gt.3) then
      write(*,*)'3XQ new_ternary_asymmetry illegal Toop',new_toop
      goto 1000
   endif
   asymmetric(new_toop:new_toop)='T'
   ia=1
   mqf=>phres%mqmqaf
   t1=asymter
   tersys(t1)%asymm=asymmetric
   tersys(t1)%isasym(new_toop)=1
!
   if(verbose) write(*,17)t1,tersys(t1)%asymm,tersys(t1)%isasym
17 format('3XQ line 4264 in_new_ternary_asym ',i3,'  "',a,'" ',3i3)
!
! A new asymetric ternary, reinitiate varkappa
!      call nathalie_asym(phres,verbose)
!
   do nn=1,size(mqf%compvar)
      box=>mqf%compvar(nn)
      cat1=box%cat1; cat2=box%cat2
      if(verbose) write(*,18)box%seq,cat1,cat2
18    format('3XQ set asymmetries for varkappa ',i2,' with cations: ',2i3)
! inititiate as symmetric, use the [ ... ] facility to store integers
! ijklx generate the index of the quadruple to the cation
!   this_ij=[ijklx(cat1,cat1,ia,ia)]
!   this_ji=[ijklx(cat2,cat2,ia,ia)]
!   denom_ij=[ijklx(cat1,cat2,ia,ia)]
! remove any previous values in box%
   if(allocated(box%ivk_ij)) deallocate(box%ivk_ij)
   if(allocated(box%jvk_ji)) deallocate(box%jvk_ji)
   if(allocated(box%kvk_ijk)) deallocate(box%kvk_ijk)
! Here this_ij and this_ji are inititaed for a symmetric ternary
   box%ivk_ij=[ijklx(cat1,cat1,ia,ia)]
   box%jvk_ji=[ijklx(cat2,cat2,ia,ia)]
! All quadruplets for ivk_ij and jvk_ji are added automatically to divisor
! but we may have to add explicitly mixed quadruplets below
!   box%kvk_ijk=[ijklx(cat1,cat2,ia,ia)]
! inside this loop set asymmetries for all ternaries
      call loop_all_ternaries(phres,box,verbose)
   enddo
! Adding mixed terms ....how?
!
!-------------------------------------------------------------------
1000 continue
   if(verbose) then
      call list_compvar(phres)
      write(*,*)'3XQ exit new_ternary_asym'
   endif
!
   return
 end subroutine new_ternary_asym

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
!
!\addtotable subroutine loop_all_ternaries
!\begin{verbatim}
 subroutine loop_all_ternaries(phres,box,inverbose)
!
! This is at least the 5th time I rewrite this routine
! It is called for a binary in "box" and tests if this binary i-j 
! has any asymmetries by looping all ternaries in tersys to check
! Any ternary with an asymmetry can modify the expressions for 
! the box%ivk_ij and box%jvk_ji and their denominator %kvk_ijk
! The %ivk_ij, %jvk_ji and %kvk_ijk arrays are sums of x_ij fraction variables
!
! Inside this routine three local allocatable arrayes keeps track of the
! asymmetries and at the end they are copied to the global box$ivk_ij etc
!
   implicit none
   type(gtp_phase_varres), pointer :: phres
   type(gtp_allinone), pointer :: box
   logical inverbose
!\end{verbatim}
!
   integer, allocatable, dimension(:) :: this_ij, this_ji, this_denom
   type(gtp_mqmqa_var), pointer :: mqf
   integer ia,cat1,cat2,t1,toopel,toop,elk,elm,nn,third(3),mix,mjx,myx,mm
   logical verbose
!   
   mqf=>phres%mqmqaf
   verbose=inverbose
!
   if(verbose) write(*,10)box%seq,box%cat1,box%cat2
10 format(/'3XQ loop ternaries for varkappa ',i3,5x,2i3)
   ia=1
   cat1=box%cat1; cat2=box%cat2
! inititiate as symmetric, use the [ ... ] facility to store integers
! ijklx generate the index of the quadruple to the cation
!   this_ij=[ijklx(cat1,cat1,ia,ia)]
!   this_ji=[ijklx(cat2,cat2,ia,ia)]
!   denom_ij=[ijklx(cat1,cat2,ia,ia)]
   this_ij=[cat1]
   this_ji=[cat2]
! MOVED TO CALLING ROUTINE remove any previous values in box%
!   if(allocated(box%ivk_ij)) deallocate(box%ivk_ij)
!   if(allocated(box%jvk_ji)) deallocate(box%jvk_ji)
!   if(allocated(box%kvk_ijk)) deallocate(box%kvk_ijk)
! Here this_ij and this_ji are inititaed for a symmetric ternary
!   box%ivk_ij=[ijklx(cat1,cat1,ia,ia)]
!   box%jvk_ji=[ijklx(cat2,cat2,ia,ia)]
! All quadruplets for ivk_ij and jvk_ji are added automatically to divisor
! but we may have to add explicitly mixed quadruplets below
!   box%kvk_ijk=[ijklx(cat1,cat2,ia,ia)]
!
!
   terloop:do t1=1,size(tersys)
! first check if this varkappa is part of this ternary!!
      
      toop=index(tersys(t1)%asymm,'T')
!      if(tersys(t1)%asymm.eq.'KKK') then
      if(toop.le.0 .or. toop.gt.3) then
         cycle terloop
      endif
      toopel=tersys(t1)%el(toop)
      if(verbose) write(*,7)t1,cat1,cat2,toop,toopel
7     format('3XQ ternary ',i2,' is asymmetric ',2i3,5x,2i3)
! check if this box is part of this ternary
! note cat1<cat2 and tersys(t1)%el(1..3) ordered increasingly
      do nn=1,3
         if(cat1.eq.tersys(t1)%el(nn)) then
            do mm=nn+1,3
               if(cat2.eq.tersys(t1)%el(mm)) then
                  goto 12
               endif
            enddo
         endif
      enddo
! this can be detected earlier ....
      if(verbose) write(*,11)cat1,cat2,tersys(t1)%el
11    format('3XQ the binary ',2i2,' is not part of the ternary ',3i2)
      cycle terloop
!
! this box is part of this ternary !!
12    continue
! toop is 1, 2 or 3.  toopel is actual index of Toop cation
! we have a Toop ternary %el(1..3) and a binary cat1-cat2
! if neither cat1 nor cat2 is a Toop element this box is unchanged
      if(toopel.ne.cat1 .and. toopel.ne.cat2) then
         cycle terloop
      endif
! This seems clumsy but we need to know the third cation
!  tersys(t1)%el is ordered increasingly and cat2 > cat1
      find3rd: do nn=1,3
         if(tersys(t1)%el(nn).eq.cat1 .or. &
              tersys(t1)%el(nn).eq.cat2) cycle find3rd
         goto 50
      enddo find3rd
      cycle terloop
!
! found the third cation
50    continue
      elk=tersys(t1)%el(nn)
! If toopel is neither box%cat1 not box%cat2 loop
!---------------------------------------------
! Nathalie algorithm, similar to the one I already tried but messed up
! Loop on the binaries ij
!    vk_ij = x_ii
!    vk_ji= x_jj
!    denominator = x_ii + x_jj
!    Loop on the ijk
!        if i Toop
!            vk_ji += x_kk 
!            denominator += x_kk
!        if j Toop    
!            vk_ij += x_kk
!            denominator += x_kk
!            Add the mixed terms
!----------------------------------
      if(toopel.eq.cat1) then
! cat1 is the Toop, then add elk to this_ji
! addquad is .TRUE. if the cation elk is not alreay in this_ij
!         if(addquad(box%jvk_ji,elk,elk) then
         if(addcat(this_ji,elk)) then
! addcat return .TRUE. if elk is not already in this_ji
            this_ji = [this_ji, elk]
         endif
      elseif(toopel.eq.cat2) then
         if(addcat(this_ij, elk)) then
! addcat return .TRUE. if elk is not already in this_ij
            this_ij = [this_ij, elk]
         endif
      endif
   enddo terloop
!----------------------------------------------------------
!
! When we are here we have looped all ternaries for one varkappa box
! and collected all asymmetries involving this box.
! Generate the arrays ivk_ij, jvk_ji and kvk_ijk with quadruplet indices
! if no asymmetries size(this_ij) and size(this_ji) is unity 
   if(verbose) then
      write(*,300)'this_ij: ',size(this_ij),this_ij
      write(*,300)'this_ji: ',size(this_ji),this_ji
300   format('3XQ After terloop: ',a,' with ',i2,' cations: ',10i3)
   endif
   mixed1a: do mix=2,size(this_ij)
! Add any asymmetric cations quadruplets and their mixed ones
      elk=this_ij(mix)
      box%ivk_ij=[box%ivk_ij, ijklx(elk,elk,ia,ia), ijklx(cat1,elk,ia,ia)]
      mixed1b: do myx=mix+1,size(this_ij)
! if there are 2 or more ternaries with toop ... add their mixed quadruplets
         elm=this_ij(myx)
         if(addcat(this_ij, elm)) then
            box%ivk_ij=[box%ivk_ij, ijklx(elk,elm,ia,ia)]
         endif
      enddo mixed1b
   enddo mixed1a
!
   mixed2a: do mjx=2,size(this_ji)
      elk=this_ji(mjx)
      box%jvk_ji=[box%jvk_ji, ijklx(elk,elk,ia,ia), ijklx(cat2,elk,ia,ia)]
      mixed2b: do myx=mjx+1,size(this_ji)
! if there are 2 or more ternaries with toop ... add their mixed quadruplets
         elm=this_ji(myx)
         if(addcat(this_ji, elm)) then
            box%jvk_ji=[box%jvk_ji, ijklx(elk,elm,ia,ia)]
         endif
      enddo mixed2b
   enddo mixed2a
!
! The denominator kvk_ijk is the sum of ivk_ij and jvk_ji + some mixed terns
   if(verbose) write(*,180)size(this_ij),size(this_ji)
180 format('3XQ now the denominator',2i3)
   box%kvk_ijk=[ijklx(cat1,cat2,ia,ia)]
!
   mixed1c: do mix=1,size(this_ij)
      elk=this_ij(mix)
! problem with testcase for vk_34/vk_43 with x_11 in ij and just x_44 in ji
      mixed2c: do mjx=1,size(this_ji)
! loop from 1 as first item in this_ij x_ij is default
! BUT this may generate many duplicate ij
         elm=this_ji(mjx)
         if(elk.eq.elm) cycle mixed2c
         if(verbose) write(*,190)mix,mjx,elk,elm,size(this_ji)
190      format('3XQ Indices: ',2i3,' representing pair ',2i3,' size: ',i3)
         if(addcat(this_ij, elm) .and. &
              .not.(elk.eq.cat1 .and. elm.eq.cat2)) then
! the .not. above to avoid duplicate ...
            if(verbose) write(*,200)elk,elm,box%seq
200         format('3XQ cation pair ',2i3,' added to denominator ',i3)
            box%kvk_ijk=[box%kvk_ijk, ijklx(elk,elm,ia,ia)]
         else
            if(verbose) write(*,210)elk,elm,box%seq
210         format('3XQ cation pair ',2i3,' already in denominator: ',i3)
         endif
      enddo mixed2c
   enddo mixed1c
!
   if(verbose) write(*,*)'Exit loop_all_ternaries',box%seq,size(box%kvk_ijk)
1000 continue
   return
 end subroutine loop_all_ternaries

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable function addcat
!\begin{verbatim}
 logical function addcat(catlist, icat)
! This return TRUE if icat is NOT included in icat
   implicit none
   integer, dimension(:) :: catlist
   integer icat
!\end{verbatim}
   integer jj,kk
!   
   kk=size(catlist)
!   write(*,100)icat,kk,(catlist(jj),jj=1,kk)
100 format('3XQ *** In addcat: ',i2,' size: ',i2,' cations: ',10i2)
   if(kk.eq.0) goto 900
   do jj=1,kk
      if(icat.eq.catlist(jj)) goto 1100
   enddo
900 continue
   addcat=.true.
!
1000 continue
!   write(*,1010)
1010 format('3XQ leaving addcat')
   return
!
1100 continue
   write(*,1020)icat,jj,(catlist(jj),jj=1,kk)
1020 format('3XQ addcat supressed duplicate ',2i3,' array ',10i3)
   addcat=.false.
   goto 1000
 end function addcat

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
!
!\addtotable subroutine loop_all_ternaries_failed
!\begin{verbatim}
   subroutine loop_all_ternaries_failed(phres,box,inverbose)
!
! This is at least the 5th time I rewrite this routine
! It is called for a binary in "box" and tests if this binary i-j 
! and loops all ternaries in tersys to check if it is part of!
! any ternary with asymmetry which would modify the expressions for 
! the box%ivk_ij and box%jvk_ji and their denominator %kvk_ijk for the binary
! The %ivk_ij, %jvk_ji and %kvk_ijk arrays has sums of x_ij fraction variables

   implicit none
   type(gtp_phase_varres), pointer :: phres
   type(gtp_allinone), pointer :: box
   logical inverbose
!\end{verbatim}
!
!-------------------------- my problem --------------------------------
! ternaries: 1             2             3            4              =4
! cations    1-2-3         1-2-4         1-3-4        2-3-4
! binaries   1-2 1-3 2-3   1-2 1-4 2-4   1-3 1-4 3-4  2-3 2-4 3-4
!            1   2   4     1   3   5     2   3   6    4   5   6
!
! 11-22-33           11-22-44           11-33-44           22-33-44
! 11-22 11-33 22-33  11-22 11-44 22-44  11-33 11-44 33-44  22-33 22-44 33-44
! 1     2     4      1     3     5      2     3     6      4     5     6
!
! binaries: 1     2     3     4     5     6                         =6
!           1-2   1-3   1-4   2-3   2-4   3-4
!
! quads, the independent fraction variabled, include constituents between each
! index     11  12  13  14  22  23   24  33  34  44
! sequence  1   2   4   4   5   6    7   8   9   10                =10
! cations:  1               2            3       4                 =4
!---------------------------------------------------------------------
!
! Find the binary compvar and cations involved
!   
! we should loop all ternaries with all compvar and set asymmetries
   type(gtp_mqmqa_var), pointer :: mqf
   type(gtp_allinone), pointer :: box1,box2
! savenu and savegamma has x_ij crossterms for multiple asymmetric ternaries
   integer, dimension(:), allocatable :: savenu
   integer, dimension(:), allocatable :: savegamma
! vz is the ternary cation, t1 is the ternary index
   integer haha,vz,toopel,toopem
   integer icat,jcat,ia,gg,emcat1,emcat2
   integer mii,mij,mjj,new_toop,nnn
   logical verbose
!
   integer abrakadabra,i,j,k,l,m,ny,iv,jv,t1,tcat1,tcat2,tcat3,vzem
!
! first maybe remove all ivk_ij etc? from box?
! that is ONLY necessary if one can remove an asymmetry. Assume not!
! set the appropriate indices in ivk_ij, jvk_ij, kvk_ijk
! 
! while debugging
   if(inverbose) verbose=.true.
!   verbose=.true.
   icat=box%cat1
   jcat=box%cat2
! wow, I did not know i saved these also ...'
! these are the cation quad fraction indices in x_ij
   emcat1=mqmqa_data%emquad(icat)
   emcat2=mqmqa_data%emquad(jcat)
! maybe useful sometimes ....
   if(verbose) write(*,1)icat,jcat,emcat1,emcat2
1  format(/'3XQ in loop_all_ternaries for binary: ',2i2,' or emquads: ',2i2)
!
! These are the cations of the ternary, their values are 1..n
   if(allocated(savenu)) then
      deallocate(savenu)
   endif
   if(allocated(savegamma)) then
      deallocate(savegamma)
   endif
!   write(*,777)box%elcat1,box%elcat2,box%quadicat1,box%quadicat2
!777 format('3XQ box%elcat: ',2i3,' box%quadicat: ',2i3)
! this is the single anion
   ia=1
!   if(verbose) write(*,3)icat,jcat,emcat1,emcat2
3  format('3XQ loop_all_ternaries to set asymmetries for box: ',2i3,2x,2i3)
!   write(*,*)'3XQ in loop_all_ternaries 2'
! is tersys initiated?
!   do t1=1,size(tersys)
!      if(verbose) write(*,555)t1,tersys(t1)%el,tersys(t1)%emquad,&
!           tersys(t1)%asymm
!555   format('3XQ tersys: ',i2,2x,3i3,2x,3i3,3x,a)
!   enddo
!   
   bigloop: do t1=1,size(tersys)
! we assume KKK asymmetries are default
      if(tersys(t1)%asymm.eq.'KKK') then
! no Toop element
         if(verbose) write(*,7)t1
7        format('Ternary ',i2,' is symmetrical, skipped')
         cycle bigloop
      endif
!
      new_toop=index(tersys(t1)%asymm,'T')
      if(new_toop.lt.1 .or. new_toop.gt.3) then
         write(*,*)'3XQ Asymmetric but no T: ',tersys(t1)%asymm,t1
         stop 'fatal error'
      endif
! the value of new_toop is always 1,2,3.  In toopel set cation
      toopel=tersys(t1)%el(new_toop)
      toopem=tersys(t1)%emquad(new_toop)
      if(verbose) write(*,4)t1,tersys(t1)%el,tersys(t1)%emquad,&
           tersys(t1)%asymm,emcat1,emcat2,toopem
4     format('3XQ line 4215 ternary ',i3,5x,3i3,5x,3i3,' "',a,'" ',2i3,3x,i3)
!
! This routine is NOT time critcal, a lot of sloppy coding ...
! because my head start to turn each time I try to program this ....
! Check if the binary i-j is part of this ternary
      iin3: do iv=1,3
!         if(icat.eq.tersys(t1)%emquad(iv)) goto 70
         if(emcat1.eq.tersys(t1)%emquad(iv)) goto 70
      enddo iin3
      if(verbose) write(*,69)t1,icat,emcat1
69    format('3XQ Ternary ',i2,' skipped as ',2i3,' not part of it')
      cycle bigloop
70    continue
      jin3: do jv=1,3
!         if(jcat.eq.tersys(t1)%emquad(jv)) goto 80
         if(emcat2.eq.tersys(t1)%emquad(jv)) goto 80
      enddo jin3
      if(verbose) write(*,69)t1,jcat,emcat2
      cycle bigloop
!
!
! the binary i-j is part of this ternary, fix asymmetries *******
!
80    continue
      if(verbose) write(*,81)iv,jv,icat,jcat,emcat1,emcat2,toopel,toopem
81    format('3XQ in loop_all_ternaries 4: ',2i3,5x,2i3,5x,2i3,5x,2i3)
! iv and jv indicate the binary, find the 3rd cation ....  very clumsy ...
      findvz1: do vz=1,3
         if(vz.eq.iv) exit findvz1
      enddo findvz1
      if(verbose) write(*,*)'3XQ looking for third cation 1: ',vz,iv
      if(vz.gt.3) cycle bigloop
      findvz2: do vz=1,3
         if(vz.eq.iv) cycle findvz2
         if(vz.ne.jv) exit findvz2
      enddo findvz2
      if(verbose) write(*,*)'3XQ looking for third cation 2: ',vz,jv
      if(vz.gt.3) then
         if(verbose) &
         write(*,*)'3XQ skip this ternary as it does not cotain this binary'
         cycle bigloop
      endif
!
! The binary icat-jcat is included in this ternary
      vzem=tersys(t1)%emquad(vz)
      if(verbose) then
         write(*,82)vz,vzem
82       format('3XQ Found third cation ',2i3)
         write(*,83)emcat1,emcat2,vzem,emcat1,emcat2
83       format('3XQ All 3 cations ',3i3,' for binary ',i1,'-',i1)
         write(*,85)icat,jcat,t1,vz
85       format('3XQ the binary ',i1,'-',i1,' in ternary ',i2,', with ',i2)
      endif
! These are the cations of the ternary
! asymmetries to this box can be added from several ternaries
! toopel is 1, 2 or 3.  toopem is the actual cation index, 1 to n
!
      toopel=tersys(t1)%el(new_toop)
      toopem=tersys(t1)%emquad(new_toop)
      tcat1=tersys(t1)%emquad(1)
      tcat2=tersys(t1)%emquad(2)
      tcat3=tersys(t1)%emquad(3)
      if(verbose) write(*,90)tcat1,tcat2,tcat3,emcat1,emcat2,toopel,toopem
!      write(*,90)tcat1,tcat2,tcat3,emcat1,emcat2,vzem,toopel,toopem
90    format('3XQ Ternary ',3i3,', binary: ',3i3,', Toop: ',2i3)
! select 
      if(toopem.eq.emcat1) goto 200
      if(toopem.eq.emcat2) goto 150
      if(verbose) write(*,*)'3XQ Toop cation not in this binary'
      cycle bigloop
!
150   continue
! here toopel is jcat or emcat2, add x_jcat,vz and \nu
! icat is asymmetric in this ternary, add terms in jvk_ij and in savegamma
      if(verbose) write(*,155)emcat1,vzem,'jvk_ji'
155   format('3XQ add ',2i3,' to ',a,' ********* unfinished')
!      cycle bigloop
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(addquad(box%jvk_ji, jcat,vz)) then
! addquad is FALSE if jcat,vz already in box%jvk_ji
         if(verbose) write(*,17)4239,jcat,vz, vz,vz, jcat,icat
17       format('3XQ line ',i4,' x_',2i1,' and x_',2i1,' added to jvk_',2i1)
! an elegant Fortran assignment of an additional items in an allocatable
         box%jvk_ji=[box%jvk_ji, ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
!-------------added to jvk_ji
! quad fractions added to jvk_ij added also to denominator, add ijklx(icat,vz)
!         if(verbose) write(*,18)4245,icat,vz
18       format('3XQ line ',i5,' added x_',2i1,' to box%kvk_ijk')
!         box%kvk_ijk=[box%kvk_ijk, ijklx(icat,vz,ia,ia)]
!         if(verbose) write(*,18)4247,icat,vz
!         box%all_ijk=[box%all_ijk, ijklx(jcat,vz,ia,ia), &
!              ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
      else
! skip adding to savenu and savegamma also
         if(verbose) &
              write(*,99)jcat,vz,ijklx(jcat,vz,ia,ia),' jvk_ji:',box%jvk_ji
99       format('3XQ quad x_',2i1,' or ',i3,a,10i3)
         goto 500
      endif
!      write(*,*)'3XQ line 4615 added quads to box%all_ijk'
!
! savenu is cross terns related to ij, savegamma to ji  ---------------------
!      write(*,*)'3XQ line 4617 check savenu'
      if(allocated(savegamma)) then
         if(verbose) write(*,373)'jcat use \nu',size(savegamma),savegamma
373      format('3XQ ',a,' mixed asymmetry terms',i3,': ',10i3)
         do gg=1,size(savenu)
! the mixed terms with \gamma should should be added to jvk_ji
            box%all_ijk=[box%all_ijk, ijklx(vz,savenu(gg),ia,ia)]
!                  write(*,374)'3XQ added ji',savenu(gg),vz
!                  write(*,375)'jvk_ji ',box%jvk_ji
374         format(a,' x_',2i1)
375         format('3XQ ',a,'=',10i4)
         enddo
         savenu=[savenu, vz ]
      else
! otherwize just add vz to savenu
         savenu=[vz]
!        write(*,373)'3XQ line 4377 savednu i ',size(savenu),savenu
      endif
!------------------------------------------------------ now savegamma
! savegamma is related to ji, maybe add denominator terms
!      write(*,*)'3XQ line 4632 savegamma'
      if(allocated(savegamma)) then
         if(verbose) write(*,373)'case 1 use \gamma',size(savegamma),savegamma
         do gg=1,size(savegamma)
! the mixed terms with \gamma should should be added to kvk_ijk
            box%kvk_ijk=[box%kvk_ijk, ijklx(vz,savegamma(gg),ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(vz,savegamma(gg),ia,ia)]
!                  write(*,374)'3XQ added kvk_ijk',savegamma(gg),vz
!                  write(*,375)'kvk_ji ',box%kvk_ijk
         enddo
! do not save vz as it does no relates to ij
!               savegamma=[savegamma, vz ]
      else
! and we must add vz to savegamma
         savegamma=[vz]
!         write(*,373)'saved i ',size(savevz),savevz
      endif
! The asymmetric xi_ji is depend on y_ik update dxi_ji and dxi_ji
! Hm, dxi_ji are just sums of y_jk?
! We ignore dxi_i ......
      do nnn=1,mqmqa_data%nquad
!                box%dxi_ij(nnn)=box%dxi_ij(nnn)+dy_ik(icat,nnn)
         box%dxi_ji(nnn)=box%dxi_ji(nnn)+mqmqa_data%dy_ik(vz,nnn)
      enddo
      if(gx%bmperr.ne.0) then
         write(*,*)'3XQ ijklx index error line 4533'
         stop
      endif
      goto 500
!
!---------------------------------------------------------------------
! >>>>>>>>>>>>  THE CODE AROUND HERE IS VERY UNCERTAIN  <<<<<<<<<<<<<<
!---------------------------------------------------------------------
!     case(2) ! ***************************************************
200   continue
      if(verbose) write(*,155)emcat2,vzem,'ivk_ij'
!      cycle bigloop
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! jcat is asymmetric in box, same as for icat just change icat to jcat!!!!
! and save in jvk_ji ...
      if(addquad(box%ivk_ij, icat,vz)) then
         if(verbose) write(*,27)4318, icat,vz, vz,vz, icat,jcat
27       format('3XQ line ',i5,' adding x_',2i1,' and x_',2i1,' to ivk_',2i1)
         box%ivk_ij=[box%ivk_ij, ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! Nath noted missing  ijklx(vz1,vz2,ia,ia) if icat and jcat are asymmetrical
!         if(verbose) write(*,28)4322,jcat,vz
!         box%kvk_ijk=[box%kvk_ijk, ijklx(jcat,vz,ia,ia)] ??????????
28       format('3XQ line ',i5,' added x_',2i1,' to box%kvk_ijk')
! ???        if(verbose) write(*,18)4320,jcat,vz
! ********** I think this is wrong ??????????????????????????????
!         if(verbose) write(*,18)4320,icat,vz
!         box%all_ijk=[box%all_ijk, ijklx(icat,vz,ia,ia), &
!              ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
!
      else
!         write(*,99)jcat,vz,ijklz(jcat,vz,ia,ia),box%jvk_ji
         if(verbose) &
              write(*,99)icat,vz,ijklx(icat,vz,ia,ia),' in ivk_ij: ',box%ivk_ij
         goto 500
      endif
! if savegamma allocated we must add terms to ivk_ijk
      if(allocated(savegamma)) then
         if(verbose) write(*,373)'case 2 use \gamma',size(savegamma),savegamma
         do gg=1,size(savegamma)
! mixed terms with \gamma should be added to ivk_ij
            box%ivk_ij=[box%ivk_ij, ijklx(vz,savegamma(gg),ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(vz,savegamma(gg),ia,ia)]
!                  write(*,374)'3XQ added ij',savegamma(gg),vz
!                  write(*,375)'ivk_ij ',box%ivk_ij
         enddo
         savegamma=[savegamma, vz ]
      else
! and we must add vz to savevz
         savegamma=[ vz ]
!        write(*,373)'savedgamma j ',size(savegamma),savegamma
      endif
!------------------------------------------------------ now savenu
! savegamma is related to ij, maybe add denominator terms
      if(allocated(savenu)) then
         if(verbose) write(*,373)'case 2 use \nu',size(savenu),savenu
         do gg=1,size(savenu)
! the mixed terms with \nu should should be added to kvk_ijk
            box%kvk_ijk=[box%kvk_ijk, ijklx(vz,savenu(gg),ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(vz,savenu(gg),ia,ia)]
!                  write(*,374)'3XQ added kvk_ijk',savenu(gg),vz
!                  write(*,375)'jvk_ji ',box%kvk_ijk
         enddo
      else
! add vz to savenu
         savenu=[ vz ]
      endif
! The asymmetric xi_ij is depend on y_ik update dxi_ij and dxi_ji
! Hm, dxi_ji are just sums of y_ik?
      do nnn=1,mqmqa_data%nquad
         box%dxi_ij(nnn)=box%dxi_ij(nnn)+mqmqa_data%dy_ik(vz,nnn)
      enddo
      if(gx%bmperr.ne.0) then
         write(*,*)'3XQ ijklx index error line 4578'
         stop
      endif
!
! endcases
!-----------------------------------------------------------
! *****************************************
500 continue
! maybe add code below ??
! savenu and savegamma are needed for multiple asymmetries
!      if(allocated(savenl cu)) then
!         write(*,*)'3XQ line 4725 add mixed quades with savenu'
!      endif
!      if(allocated(savegamma)) then
!         write(*,*)'3XQ line 4728 add mixed quades with savegamma'
!      endif
      if(verbose) then
         write(*,505)'3XQ savenu:    ',savenu
         write(*,505)'3XQ savegamma: ',savegamma
505      format('3XQ maybe ',a,' cross terms:',10i3)
      endif
      goto 700
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! code below skipped but needed to integrate savenu and savegamma in vk_ij etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loops below now redundant when we added savevz loops above ..... ????
! code handling kvk_ijk terms due to extra x_ii and x_jj in ivk_ij and jvk_ji
! copied from end of calcasymvar to avoid it is repeted at all calculations
! skip first ivk_ij
      addkvkterm: do j=2,size(box%ivk_ij)
         do k=1,size(mqmqa_data%emquad)
            if(box%ivk_ij(j).eq.mqmqa_data%emquad(k)) then
! we have an endmember quad in ivk_ij (in addition to the first)
! Check if we have another endmember quad in jvk_ji, skip first jvk_ji
               do l=2,size(box%jvk_ji)
                  neverending: do m=1,size(mqmqa_data%emquad)
                     if(box%jvk_ji(l).eq.mqmqa_data%emquad(m)) then
                        if(k.ne.m) then
! we have 2 different endmember quads in ivk_ij and jvk_ji, 
! if the mixed quad is not alreay present add it
                           ny=ijklx(k,m,ia,ia)
                           do abrakadabra=1,size(box%kvk_ijk)
! check if this quad not already in box_kvk_ijk
                              write(*,*)'3XQ check duplicate line 4764 !!'
                           enddo
! add this quad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           box%kvk_ijk=[box%kvk_ijk, ijklx(k,m,ia,ia)]
!                              write(*,806)i,k,m,ijklx(k,m,ia,ia)
!                              write(*,805)'kvk_ijk ',box%kvk_ijk
                        endif
                     endif
                  enddo neverending
                  if(gx%bmperr.ne.0) then
                     write(*,*)'3XQ ijklx index error line 4645'
                     stop
                  endif
               enddo
            endif
         enddo
      enddo addkvkterm
805   format(a,20i3)
806   format('3XQ adding mixed quad to kvk_ijk',i3,2x,2i3,2x,i3)
! end copied code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! currently skipping code above
700      continue
!
   enddo bigloop
!***********************************************************************
!
1000  continue
   if(verbose) write(*,1010)emcat1,emcat2
!   write(*,*)emcat1,emcat2
1010 format('3XQ The binary ',2i2,' has made loop_all_ternaries')
   return
 end subroutine loop_all_ternaries_failed
 
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\begin{verbatim}
 subroutine loop_all_ternaries_old(phres,box,verbose)
! Check asymmetries in all ternaries for cations in box 
! and initiate vk_ij if necessary
! save information in the binary box%ivk_ij and box%ix_ij
! the code extracted from the old varkappa1
   implicit none
   type(gtp_phase_varres), pointer :: phres
   logical verbose
!\end{verbatim}
   type(gtp_mqmqa_var), pointer :: mqf
   type(gtp_allinone), pointer :: box,box1,box2
   integer, dimension(:), allocatable :: savenu
   integer, dimension(:), allocatable :: savegamma
   integer haha,vz,toopel,toopem
   integer icat,jcat,ia,gg,t1
   integer mii,mij,mjj,new_toop,nnn
!
   integer abrakadabra,i,j,k,l,m,ny
!
! first maybe remove all ivk_ij etc? from box?
! that is ONLY necessary if one can remove an asymmetry. Assume not!
! set the appropriate indices in ivk_ij, jvk_ij, kvk_ijk
! 
! savenu and savegamma are cross asymmetries which will be picked up at the end
   if(allocated(savenu)) then
      deallocate(savenu)
   endif
   if(allocated(savegamma)) then
      deallocate(savegamma)
   endif
!
   ia=1
   if(verbose) write(*,3)box%cat1,box%cat2,&
        allocated(savenu),allocated(savegamma)
3  format(/'3XQ loop_all_ternaries to set asymmetries for box: ',2i3,2x,2l2)

! is tersys initiated?
!   do t1=1,size(tersys)
!      if(verbose) write(*,555)t1,tersys(t1)%el,tersys(t1)%emquad,&
!           tersys(t1)%asymm
555   format('3XQ tersys: ',i2,2x,3i3,2x,3i3,3x,a)
!   enddo
   
   bigloop: do t1=1,size(tersys)
! we assume KKK asymmetries are default
      if(verbose) write(*,4)t1,tersys(t1)%el,tersys(t1)%asymm
!      write(*,4)t1,tersys(t1)%el,tersys(t1)%asymm
4     format('3XQ line 4215 ternary ',i3,5x,3i3,' "',a,'"')
      if(tersys(t1)%asymm.eq.'KKK') then
! no Toop element
         cycle bigloop
      endif
!
      new_toop=index(tersys(t1)%asymm,'T')
! I assume there is just one T
      if(new_toop.le.0 .or. new_toop.gt.3) then
         write(*,*)'3XQ Asymmetric but no T: ',tersys(t1)%asymm,t1
         stop 'fatal error'
      endif
! the varkappa involved are those belonging to this ternary
! asymmetries to the same box will be aded at several loops ...
      toopel=tersys(t1)%el(new_toop)
      toopem=tersys(t1)%emquad(new_toop)
! toopel is 1, 2 or 3
! toopem is actual cation index 1..n
      vz=tersys(t1)%el(3)
      icat=box%cat1
      jcat=box%cat2
! check if this box is involved in this ternary, is icat or jcat Toop?
! toopem is the cation index ...
      if(verbose) write(*,5)icat,jcat,vz,new_toop,toopel,toopem
      if(icat.ne.toopem .and. jcat.ne.toopem) then
!         if(verbose) &
              write(*,*)'3XQ this ternary do not involve these varkappa'
         cycle bigloop
      endif
! if toopel is part of this box change vk_ij or vk_ji
5     format('3XQ line 4232 icat,jcat,vz: ',3i2,' new_toop,toopel, toopem ',3i3)
      if(toopel.eq.icat) goto 200
      if(toopel.ne.jcat) cycle bigloop
!
! here toopel is jcat, add x_jcat,vz and \nu
! icat is asymmetric in this ternary, add terms in jvk_ij and in savegamma
! addquad is FALSE if jcat,vz already in box%jvk_ji
      if(addquad(box%jvk_ji, jcat,vz)) then
         if(verbose) write(*,17)4239,jcat,vz, vz,vz, jcat,icat
17       format('3XQ line ',i4,' x_',2i1,' and x_',2i1,' added to jvk_',2i1)
! an elegant Fortran assignment of an additional items in an allocatable
         box%jvk_ji=[box%jvk_ji, ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
!-------------added to jvk_ji
! quad fractions added to jvk_ij added also to denominator, add ijklx(icat,vz)
!         if(verbose) write(*,18)4245,icat,vz
18       format('3XQ line ',i5,' added x_',2i1,' to box%kvk_ijk')
!         box%kvk_ijk=[box%kvk_ijk, ijklx(icat,vz,ia,ia)]
!         if(verbose) write(*,18)4247,icat,vz
!         box%all_ijk=[box%all_ijk, ijklx(jcat,vz,ia,ia), &
!              ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
      else
! skip adding to savenu and savegamma also
         write(*,99)jcat,vz,ijklx(jcat,vz,ia,ia),' jvk_ji:',box%jvk_ji
99       format('3XQ quad x_',2i1,' or ',i3,a,10i3)
         goto 500
      endif
!      write(*,*)'3XQ line 4615 added quads to box%all_ijk'
!
! savenu is cross terns related to ij, savegamma to ji  ---------------------
!      write(*,*)'3XQ line 4617 check savenu'
      if(allocated(savegamma)) then
         if(verbose) write(*,373)'jcat use \nu',size(savegamma),savegamma
373      format('3XQ ',a,' mixed asymmetry terms',i3,': ',10i3)
         do gg=1,size(savenu)
! the mixed terms with \gamma should should be added to jvk_ji
            box%all_ijk=[box%all_ijk, ijklx(vz,savenu(gg),ia,ia)]
!                  write(*,374)'3XQ added ji',savenu(gg),vz
!                  write(*,375)'jvk_ji ',box%jvk_ji
374         format(a,' x_',2i1)
375         format('3XQ ',a,'=',10i4)
         enddo
         savenu=[savenu, vz ]
      else
! otherwize just add vz to savenu
         savenu=[vz]
!        write(*,373)'3XQ line 4377 savednu i ',size(savenu),savenu
      endif
!------------------------------------------------------ now savegamma
! savegamma is related to ji, maybe add denominator terms
!      write(*,*)'3XQ line 4632 savegamma'
      if(allocated(savegamma)) then
         if(verbose) write(*,373)'case 1 use \gamma',size(savegamma),savegamma
         do gg=1,size(savegamma)
! the mixed terms with \gamma should should be added to kvk_ijk
            box%kvk_ijk=[box%kvk_ijk, ijklx(vz,savegamma(gg),ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(vz,savegamma(gg),ia,ia)]
!                  write(*,374)'3XQ added kvk_ijk',savegamma(gg),vz
!                  write(*,375)'kvk_ji ',box%kvk_ijk
         enddo
! do not save vz as it does no relates to ij
!               savegamma=[savegamma, vz ]
      else
! and we must add vz to savegamma
         savegamma=[vz]
!         write(*,373)'saved i ',size(savevz),savevz
      endif
! The asymmetric xi_ji is depend on y_ik update dxi_ji and dxi_ji
! Hm, dxi_ji are just sums of y_jk?
! We ignore dxi_i ......
      do nnn=1,mqmqa_data%nquad
!                box%dxi_ij(nnn)=box%dxi_ij(nnn)+dy_ik(icat,nnn)
         box%dxi_ji(nnn)=box%dxi_ji(nnn)+mqmqa_data%dy_ik(vz,nnn)
      enddo
      if(gx%bmperr.ne.0) then
         write(*,*)'3XQ ijklx index error line 4533'
         stop
      endif
      goto 500
!
!---------------------------------------------------------------------
! >>>>>>>>>>>>  THE CODE AROUND HERE IS VERY UNCERTAIN  <<<<<<<<<<<<<<
!---------------------------------------------------------------------
!     case(2) ! ***************************************************
200   continue
! jcat is asymmetric in box, same as for icat just change icat to jcat!!!!
! and save in jvk_ji ...
      if(addquad(box%ivk_ij, icat,vz)) then
         if(verbose) write(*,27)4318, icat,vz, vz,vz, icat,jcat
27       format('3XQ line ',i5,' adding x_',2i1,' and x_',2i1,' to ivk_',2i1)
         box%ivk_ij=[box%ivk_ij, ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! Nath noted missing  ijklx(vz1,vz2,ia,ia) if icat and jcat are asymmetrical
!         if(verbose) write(*,28)4322,jcat,vz
!         box%kvk_ijk=[box%kvk_ijk, ijklx(jcat,vz,ia,ia)] ??????????
28       format('3XQ line ',i5,' added x_',2i1,' to box%kvk_ijk')
! ???        if(verbose) write(*,18)4320,jcat,vz
! ********** I think this is wrong ??????????????????????????????
!         if(verbose) write(*,18)4320,icat,vz
!         box%all_ijk=[box%all_ijk, ijklx(icat,vz,ia,ia), &
!              ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
!
      else
!         write(*,99)jcat,vz,ijklz(jcat,vz,ia,ia),box%jvk_ji
         write(*,99)icat,vz,ijklx(icat,vz,ia,ia),' in ivk_ij: ',box%ivk_ij
         goto 500
      endif
! if savegamma allocated we must add terms to ivk_ijk
      if(allocated(savegamma)) then
         if(verbose) write(*,373)'case 2 use \gamma',size(savegamma),savegamma
         do gg=1,size(savegamma)
! mixed terms with \gamma should be added to ivk_ij
            box%ivk_ij=[box%ivk_ij, ijklx(vz,savegamma(gg),ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(vz,savegamma(gg),ia,ia)]
!                  write(*,374)'3XQ added ij',savegamma(gg),vz
!                  write(*,375)'ivk_ij ',box%ivk_ij
         enddo
         savegamma=[savegamma, vz ]
      else
! and we must add vz to savevz
         savegamma=[ vz ]
!        write(*,373)'savedgamma j ',size(savegamma),savegamma
      endif
!------------------------------------------------------ now savenu
! savegamma is related to ij, maybe add denominator terms
      if(allocated(savenu)) then
         if(verbose) write(*,373)'case 2 use \nu',size(savenu),savenu
         do gg=1,size(savenu)
! the mixed terms with \nu should should be added to kvk_ijk
            box%kvk_ijk=[box%kvk_ijk, ijklx(vz,savenu(gg),ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(vz,savenu(gg),ia,ia)]
!                  write(*,374)'3XQ added kvk_ijk',savenu(gg),vz
!                  write(*,375)'jvk_ji ',box%kvk_ijk
         enddo
      else
! add vz to savenu
         savenu=[ vz ]
      endif
! The asymmetric xi_ij is depend on y_ik update dxi_ij and dxi_ji
! Hm, dxi_ji are just sums of y_ik?
      do nnn=1,mqmqa_data%nquad
         box%dxi_ij(nnn)=box%dxi_ij(nnn)+mqmqa_data%dy_ik(vz,nnn)
      enddo
      if(gx%bmperr.ne.0) then
         write(*,*)'3XQ ijklx index error line 4578'
         stop
      endif
!
! endcases
!-----------------------------------------------------------
! *****************************************
500 continue
! maybe add code below ??
! savenu and savegamma are needed for multiple asymmetries
!      if(allocated(savenl cu)) then
!         write(*,*)'3XQ line 4725 add mixed quades with savenu'
!      endif
!      if(allocated(savegamma)) then
!         write(*,*)'3XQ line 4728 add mixed quades with savegamma'
!      endif
      if(verbose) then
         write(*,505)'3XQ savenu:    ',savenu
         write(*,505)'3XQ savegamma: ',savegamma
505      format('3XQ maybe ',a,' cross terms:',10i3)
      endif
      goto 700
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! code below skipped but needed to integrate savenu and savegamma in vk_ij etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loops below now redundant when we added savevz loops above ..... ????
! code handling kvk_ijk terms due to extra x_ii and x_jj in ivk_ij and jvk_ji
! copied from end of calcasymvar to avoid it is repeted at all calculations
! skip first ivk_ij
      addkvkterm: do j=2,size(box%ivk_ij)
         do k=1,size(mqmqa_data%emquad)
            if(box%ivk_ij(j).eq.mqmqa_data%emquad(k)) then
! we have an endmember quad in ivk_ij (in addition to the first)
! Check if we have another endmember quad in jvk_ji, skip first jvk_ji
               do l=2,size(box%jvk_ji)
                  neverending: do m=1,size(mqmqa_data%emquad)
                     if(box%jvk_ji(l).eq.mqmqa_data%emquad(m)) then
                        if(k.ne.m) then
! we have 2 different endmember quads in ivk_ij and jvk_ji, 
! if the mixed quad is not alreay present add it
                           ny=ijklx(k,m,ia,ia)
                           do abrakadabra=1,size(box%kvk_ijk)
! check if this quad not already in box_kvk_ijk
                              write(*,*)'3XQ check duplicate line 4764 !!'
                           enddo
! add this quad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           box%kvk_ijk=[box%kvk_ijk, ijklx(k,m,ia,ia)]
!                              write(*,806)i,k,m,ijklx(k,m,ia,ia)
!                              write(*,805)'kvk_ijk ',box%kvk_ijk
                        endif
                     endif
                  enddo neverending
                  if(gx%bmperr.ne.0) then
                     write(*,*)'3XQ ijklx index error line 4645'
                     stop
                  endif
               enddo
            endif
         enddo
      enddo addkvkterm
805   format(a,20i3)
806   format('3XQ adding mixed quad to kvk_ijk',i3,2x,2i3,2x,i3)
! end copied code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! currently skipping code above
700      continue
!      write(*,*)'3XQ line 4763 end bigloop in loop_all_ternaries'
!      
   enddo bigloop
!
1000  continue
!   write(*,*)'3XQ leaving loop_all_ternaries'
   return
 end subroutine loop_all_ternaries_old
            
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable function addquad
!\begin{verbatim}
 logical function addquad(quadlist, jcat,vz)
! This return TRUE if the quad with cations jcat, vz not in quadlist
   implicit none
   integer, allocatable, dimension(:) :: quadlist
   integer jcat,vz,jj
!\end{verbatim}
   integer quad,ia
   ia=1
   addquad=.false.
   quad=ijklx(jcat,vz,ia,ia)
   do jj=1,size(quadlist)
      if(quad.eq.quadlist(jj)) goto 1000
   enddo
   addquad=.true.
1000 continue
   return
 end function addquad

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine calcasymvar
!\begin{verbatim}
 subroutine calcasymvar(phres,init)
! subroutine calcasymvar(phres)
! subroutine calcasymvar(mqmqavar)
! This must be called whenever the quad fractions has changed
! If init not 0, update varkappaij, xiij expressions etc for the whole system
!                and stores them in compvar(bin) datastructure
! if init=0 just calculate new values of vk_ij, xi_ij etc
! Currently programmed ONLY for a single anion
   implicit none
   integer init
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   integer i,j,ia,seq,k,l,m,ny,abrakadabra,initvk
   type(gtp_mqmqa_var), pointer :: mqf
   character*3 asymmetry
! how to create xquad mm when we need a pointer to gtp_phase_varres?
!   type(gtp_equilibrium_data), pointer :: ceq
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
!   write(*,*)'3XQ in calcasymvar',init
   ia=1
   seq=0
! if init=0 do not initita
   if(init.ne.0) initvk=0
! initiate all asymmetry 0, earlier in init_excess_asymm, line 3134
! Database may have set some asymmetries ....
!   do i=1,size(tersys)
!      write(*,17)'before',i,tersys(i)%el,tersys(i)%isasym,tersys(i)%asymm
!17    format('3XQ ',a,' varkappa1: ',i3,2x,3i3,2x,3i2,' "',a,'"')
!   enddo
!
! this loop should not always inititate all varkappa ....
!   write(*,*)'3XQ line 4441 replace call to varkapp1 in calcasymvar'
   do i=1,mqmqa_data%ncat-1
      do j=i+1,mqmqa_data%ncat
! seq specifies a binary set of elements
! results are stored in compvar(seq) for use in Gibbs energy calculations
         seq=seq+1
! asymmetry is KKK or Tx where x=1, 2 or 3
         if(mqmqder) write(*,*)'3XQ calcasymvar call varkappa1'
!         write(*,*)'3XQ call varkappa1',seq
! the 3rd argument 0 below means no asymmetry change or set all symmetrical
!         write(*,13)seq,initvk,tersys(1)%asymm
13       format('3XQ line 4446 call varkappa1 from calcasymvar',2i4,' "',a,'"')
! call varkappa7 to to set vk, xi and y using quad fractions
! taking account of asymmetries
! varkappa1 no longer has code to set asymmetries, moved to new_asymm
         call varkappa1(seq,phres,1)
      enddo
   enddo
!   do i=1,size(tersys)
!      write(*,17)'after',i,tersys(i)%el,tersys(i)%isasym,tersys(i)%asymm
!   enddo
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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
1000 continue
   return
 end subroutine calcasymvar

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine test_asymm
!\begin{verbatim}
 integer function test_asymm(t,i,j,v)
! This original version was tried to be streamlined but restored
! The ternary specified by t has 3 elements i-j-v.  v is redundant ... ?
! return 0 if neither element i nor j are asymmetric elements in this ternary
! return 1 if element i is an asymmetric element
! return 2 if element j is an asymmetric element
! return 3 if element both i and j are asymmetric elements
   implicit none
   integer t,i,j,v
!\end{verbatim}
! ONLY ONE TOOP ELEMENT PER TERNARY
!---------------------------------------------------------------
! This is previous OC version when I had correct asymmetries, asymm=TKK
! The asymmetric element here is 1
! Varkappa index:   1, summing quads: 
!   nomin: vk_12 =(x_11)/denom
!   nomin: vk_21 =(x_22+x_23+x_33)/denom
!   denom: = x_11+x_22+x_23+x_33 +x_12+x_13
!Varkappa index:   2, summing quads: 
!   nomin: vk_13 =(x_11)/denom
!   nomin: vk_31 =(x_33+x_23+x_22)/denom
!   denom: = x_11+x_33+x_23+x_22 +x_13+x_12
!Varkappa index:   3, summing quads: <<<<<<<<<<<<<<<<< no change
!   nomin: vk_23 =(x_22)/denom
!   nomin: vk_32 =(x_33)/denom
!   denom: = x_22+x_33 +x_23
!
!---------------------------------------------------------------
!
! This is previous OC version when I had correct asymmetries, asymm=KTK
! The asymmetric element here is 2
!   Varkappa index:   1, summing quads: 
!   nomin: vk_12 =(x_11+x_13+x_33)/denom
!   nomin: vk_21 =(x_22)/denom
!   denom: = x_11+x_13+x_33+x_22 +x_12+x_23
!Varkappa index:   2, summing quads:  <<<<<<<<<<<<<<<<< no change
!   nomin: vk_13 =(x_11)/denom
!   nomin: vk_31 =(x_33)/denom
!   denom: = x_11+x_33 +x_13
!Varkappa index:   3, summing quads: 
!   nomin: vk_23 =(x_22)/denom
!   nomin: vk_32 =(x_33+x_13+x_11)/denom
!   denom: = x_22+x_33+x_13+x_11 +x_23+x_12
!
!-------------------------------------------------------------   
!
! This is previous OC version when I had correct asymmetries, asymm=KKT
!Varkappa index:   1, summing quads:  <<<<<<<<<<<<<<<<< no change
!   nomin: vk_12 =(x_11)/denom
!   nomin: vk_21 =(x_22)/denom
!   denom: = x_11+x_22 +x_12
!Varkappa index:   2, summing quads: 
!   nomin: vk_13 =(x_11+x_12+x_22)/denom
!   nomin: vk_31 =(x_33)/denom
!   denom: = x_11+x_12+x_22+x_33 +x_13+x_23
!Varkappa index:   3, summing quads: 
!   nomin: vk_23 =(x_22+x_12+x_11)/denom
!   nomin: vk_32 =(x_33)/denom
!   denom: = x_22+x_12+x_11+x_33 +x_23+x_13
!   
!-----------------------------------------------------------
   integer asymmetric1,asymmetric2,hejhopp,tpos,a
   integer, save :: third=1
! we have to check %asymm, if %asymm(i:i) is not 'T' return 0
! this routine is called 3 times for each ternary, with i, j and v permuted
! EXAMPLE: a single ternary 1-2-3 with 2 as Toop
! Call i   j   v
!   1  1   2   3...n   
!   2  1   3   2 4..n
!   3  1   4   2,3 5..n
!   4  1   5   2,3,4 6..n        this is impossible ... must be reorganized
! n-1  1   n   2..n-1            how will KKT be treated ....
!   n  2   3   1 4..n
! n+1  2   4   1 3 5..n     
!  ..  2   5   1 3,4 6..n
! old   if(tersys(t)%asymm(1:1).eq.'T') then assuming T1, T2 or T3
!
! TEST_ASYMM totally removed, asymmetry expressions set separately
! we should call test_asymm 3 times for each ternary
   write(*,10)i,j,v
10 format('3XQ REDUNDANT test_asymm i,j,v: ',3i3)
   stop 'test_asymm redundant'
! 
! i,j in call is binary, v can vary from 1..n excluding i and j
! the %asymm character refers to the first two quadruplets
! MODIFIED RESTORED VERSION assuming %asymm is TKK, KTK or KKT   
   tpos=0
! Problem indexing, i and j can be larger than 3 ...
   if(i.gt.3 .or. j.gt.3) then
      write(*,17)t,i,j,v,tersys(t)%asymm,tpos,i,j
17    format('3XQ TEST_ASYMM WITH ',i3,3x,3i3,' "',a,'" gives ',i2,10x,2i3)
   endif
! icat and jcat represent two binary quads for this ternary, they represent 
! vk_ij%cat1 and vk_ij%cat2 and should be just 1 or 2?
   if(tersys(t)%asymm(1:1).eq.'T') then
      tpos=1
   elseif(tersys(t)%asymm(2:2).eq.'T') then
      tpos=2
   endif
! This means KKT is meaningless <<<<<<<<<<<<<<<<<<<<<<<<,
!
! tersys(t)%seq, %el(3), %asymm, %isasym(3)
!   tersys(t)%isasym=0
!   tersys(t)%isasym(tpos)=1
!   test_asymm=tpos
!                1          3              3           1
!                             1     3      3          1
   test_asymm=tpos
100 continue
   return
 end function test_asymm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine init_symmetric_varkappa(phres)
!\begin{verbatim}
 subroutine init_symmetric_varkappa(phres,verbose)
! New routine initiate all varkappa as symmetric
!
   implicit none
   type(gtp_phase_varres), pointer :: phres
   logical verbose
!\end{verbatim}
   type(gtp_allinone), pointer :: box
   type(gtp_mqmqa_var), pointer :: mqf
   integer seq,icat,jcat,mii,mij,mjj,ia
   
! anion
   ia=1
!
   mqf=>phres%mqmqaf
!   write(*,*)'3XQ init_symmetric_varkappa list varkappa:'
! loop for all vk_ij, where is box initiated?
! the varkappa records represent binaries but also include ternary asymmetries
   do seq=1,size(mqf%compvar)
! What is in box?
      box=>mqf%compvar(seq)
      if(verbose) write(*,10)seq,box%qcat1,box%qcat2,box%qan, &
           box%cat1,box%cat2,-box%anion,&
           box%elcat1,box%elcat2,box%elan,&
           box%quadicat1,box%quadicat2,box%quadian,&
           box%vk_ij,box%vk_ji,box%denominator,box%xi_ij,box%xi_ji
10    format('3XQ box:',i3,1x,3(a,', '),'" '/&
           3i3,3x,3i3,3x,3i3,&
           'vk: ',2(f10.6),',  denom: ',f10.6,', xi: ',2(f10.6))
! box%icat1, %icat2 are just sequential indices staring from 1 ???
! initiate symetric compvar(bin)%ivk etc.
      icat=box%cat1; jcat=box%cat2
      mii=ijklx(icat,icat,ia,ia)
      mij=ijklx(icat,jcat,ia,ia)
      mjj=ijklx(jcat,jcat,ia,ia)
      box%ivk_ij=[mii]; box%jvk_ji=[mjj]; box%kvk_ijk=[mij]
! to simplify handling derivatives the denominator is summed separately
      box%all_ijk=[mii, mjj, mij]
      box%ivk_ij=[mii]
      if(verbose) then
         if(allocated(box%ivk_ij)) then
            write(*,20)-size(box%ivk_ij),box%ivk_ij,-size(box%jvk_ji),&
                 box%jvk_ji,-size(box%kvk_ijk),box%kvk_ijk
20          format('3XQ ivk_ij mm:',20i3)
         else
            write(*,*)'3XQ line 5030 box%ivk_ij not allocated'
         endif
         if(allocated(box%asymm_nu)) then
            write(*,30)'3XQ asymm_nu:    ',box%asymm_nu
         endif
         if(allocated(box%asymm_gamma)) then
            write(*,30)'3XQ asymm_gamma: ',box%asymm_gamma
         endif
30       format(a,20i3)
! Are the dxi_ij, dy_ik etc already initiated from database ??
!      write(*,*)'3XQ line 4698 dy_ik:',mqmqa_data%dy_ik
!      write(*,*)'3XQ line 4698 dxi_ij:',mqf%compvar(1)%dxi_ij
!......
      endif
   enddo
1000 continue
!   write(*,*)'3XQ initiated compvar'
   return
 end subroutine init_symmetric_varkappa

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine loop_tersys
!\begin{verbatim}
 subroutine loop_tersys(phres)
! New routine to add a ternary asymmetry replacing varkappa1
!
! Called when the asymmetry changes, initiates ALL vk_ij, xi_i etc expressions
!
! Algorithm loops tersys and adds cations in box%iasym, %jasym etc
! These are in a second loop used to create ivk_ij, jvk_ij, kvk_ijk etc
!
   implicit none
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   type(gtp_allinone), pointer :: box
   type(gtp_mqmqa_var), pointer :: mqf
   integer ia,ni,nj,toop,toopel,nter,ix,cat1,cat2j,cat3,catT,seq,vz,binsys
   integer cati,catj,icat,jcat,mii,mij,mjj,mji,ncat
   integer emcat1,emcat2,ter
! save indices of x_ii involved in asymmetries for i-j-v ternaries
! in box there are arrays for asymm_nu and asymm_gamma
! maybe these arrays are not needed.  Varkappa1 had also savenu, savegamma
!   integer, allocatable, dimension(:) :: icat,jcat
!   integer, allocatable, dimension(:) :: jnu, igamma
   logical first
!
   ncat=mqmqa_data%ncat
   write(*,*)'3XQ testing loop_tersys',ncat
   mqf=>phres%mqmqaf
!
! add default initiation of all compvar (or make that before calling ...
   write(*,*)'3XQ loop_tersys assumes all compvar and vk_ij initiated'
!
! We have quad indices 1....n representng all quads in mqmqa phase: cat1, cat2
! We have "endmember quad" indices with a single cation, emcat1,emcat2
! We have binaries indices with all quads
! Whe have ternary indices with just the endmember quads
!
! List all ternaries
!   quadloop: do ter=1,size(tersys)
!      write(*,700)ter,tersys(ter)%el,tersys(ter)%isasym,tersys(ter)%asymm
!   enddo quadloop
700 format('3XQ tersys: ',i3,3x,3i3,3x,3i3,3x,a)
! there is a single anion with index 1
   ia=1
! save all v for i-j-v where j is Toop in jnu and
!          v for i-j-v where i is Toop in igamma   
!   if(allocated(jnu)) then
!      these are local arrays
!      deallocate(jnu); deallocate(igamma)
!   endif
   bigloop: do nter=1,size(tersys)
! loop throug all ternaries and extract asymmetric information
      write(*,30)nter,tersys(nter)%el,tersys(nter)%isasym,tersys(nter)%asymm
30    format('3XQ tersys ',i2,5x,3i3,5x,3i2,' "',a,'"')
      if(tersys(nter)%asymm.eq.'KKK') then
! this ternary is symmetrical, maybe initiate all box here?
! index to the box? is
         emcat1=tersys(nter)%el(1)
         emcat2=tersys(nter)%el(2)
         mii=ijklx(emcat1,emcat2,ia,ia)
         mij=ijklx(emcat1,emcat2,ia,ia)
         mjj=ijklx(emcat1,emcat2,ia,ia)
         cycle bigloop
      endif
! there is an asymmetry
      toop=index(tersys(nter)%asymm,'T')
      toopel=tersys(nter)%el(toop)
      if(toop.lt.1 .or. toop.gt.3) then
! UNFINISHED BELOW
! impossible error ...
         write(*,*)'3XQ illegal %asym in tersys ',toop,nter,tersys(nter)%asymm
         stop
      endif
! This is the Toop cation in this ternary
!      write(*,*)'3XQ Found asymmetric ternary',nter,toop
      catT=tersys(nter)%el(toop)
      first=.true.
! extract the actual cation index from tersys
! i-j-vz with j=toop or i-j-gamma with j=toop afects the binaries i-j
      findallcats: do ix=1,3
         if(ix.eq.toop) then
            vz=tersys(nter)%el(ix)
         else
! this is the first quadruplet in the asymmetrical ternary
            if(first) then
               first=.false.
               cati=tersys(nter)%el(ix)
! Below add x_toop/cati to %ivk_ij and %kvk_ijk to the appropriate box cati,catj
            else
               catj=tersys(nter)%el(ix)
! Below add x_catj/toop to %jvk_ji and %kvk_ijk to the appropriate box cati,catj
            endif
         endif
      enddo findallcats
! This documentation because I will soon forget ...
! Ternary has quadruplets (cati, catj, toop) in sequential order
! ternary 1  2 ...... ! ? .....     ! ... ! last ternary n*(n-1)*(n-2)/6
! quad    1  1 .. 1   ! 2  2 .. 2   ! ... ! n-2
! quad    2  2 .. n-1 ! 3  3 .. n-1 ! ... ! n-1
! quad    3  4 .. n   ! 4  5 .. n   ! ... ! n
!
! there more binary compvar vk_ij records (which are called box below)
!  1 2    n-1 ! n+1 .......................! last n*(n-1)/2 
!  1 1 .. 1   ! 2   2 .. 2 ! 3 ... 3 ! ... ! n-1    
!  2 3    n   ! 3   4    n ! 4 ... n ! ... ! n
!
! and the number of quads include also those with two identical cations
!  1 2    n ! n+1 ..................................! last n*(n+1)/2 
!  1 1 .. 1 ! 2   2 .. 2 ! 3 3 .. 3 ! ... ! n-1 n-1 ! n   
!  1 2    n ! 2   3    n ! 3 4    n ! ... ! n-1 n   ! n
! we can use the function ijklx to find the quadruplet index of 2 cations
!
! When we are here 2 binaries with an asymmetric cation icat-vz and vz-jcat
! vk_cati,catj should have x_catT,catj added
! vk_catj,cati should have x_cati,catT added
! we must also think of asymm_nu, asymm_gamma for cross-asymmetries ....
! locate the mqf%box, they are arranged 
      seq=ijklx(cati,catj,ia,ia)
      write(*,*)'3XQ cati, catj, toop: ',cati,catj,seq
! box=>mqf%compvar(seq)
      binsys=ibin(cati,catj,ncat)

      box=>mqf%compvar(seq)
   enddo bigloop
! Now implement the asymmetries

! All asymmetries installed above
!   stop 'unfinished loop_tersys'
!
! maybe not needed   call varkappa7(phres,icat,jcat,jnu,igamma)
!
1000 continue
   write(*,1010)
1010 format('3XQ leaving loop_tersys')
   return
 end subroutine loop_tersys

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine varkappa1
!\begin{verbatim}
 subroutine varkappa1(seq,phres,asymter)
! seq is an index of varkappa in mqf%compvar array of all varkappa
! asymter 
!
! This is the original varkappa1 which handles asymmetric ternaries
! It is restored because the bad version failed
! OLD: when asymter=0 then seq is ignored and looped over all boxes
! now asymter is ignored
!
! phres is pointer to gtp_phase_varres for the mqmqa phase
! should phres it be a pointer?  Does it matter?  It seems to work
!
! This routine should only be called after a change of asymmetries
! But that is just to speed up calculations and will be made later 
!
! ---------->>> removed function
! *** phres is called parres in calling routine
! asymter is the index of integer array, if zero set all symmetrical
!         if nonzero the asymmetric constituent ins already set in %asymm
!         %asymm is 'TKK' where positon of T is the assymmetric constituent
! in varkappa the cations are ordered (1,1) (1,2) ... (2,2) ... (n,n)
! in tersys the cations are ordered (1,2,3) (1,2,4) ... (2,3,4) ... (n-2,n-1,n)
! box is a record of the type(gtp_allinone)
! tersym is a structure with all combination of 3 cations for the asymmetries
! tersym(tt)%el(1) %el(2) and %el(3) are cation indices in the ternary
! tersym(tt)%isasym(1) %isasym(2) %isasym(3) is 0 or asymmetric cation index
! tersym(tt)%asymm is a 3 character variable of asymmetry DO NOT USE
! this routine may initiate, calculate and store varkappa_ij, varkappa_ji, and
!    xi_ij and xi_ji for symmetric and asymmetric systems with Kohler/Toop
! It is programmed for a single anion and just for the MQMQX phase!
!
! It will inintiate all data in box if box%lastupdate neq newXupdate
!
! some initial thinking
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
! calculate the derivatives of all vk_ij, vk_ji with respect to quads
!
   implicit none
   integer seq,asymter
!   integer seq,asymter,new_toop!
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}

! replaced original i and j by icat and jcat below!!    integer i,j,ia,bin
!
! these are quad indices of i,i, i,j abd j,j
   integer mii,mij,mjj,ia
   type(gtp_allinone), pointer :: box
   type(gtp_mqmqa_var), pointer :: mqf
!
! ia represent the single anion
! varkappa_ij and varkappa_ji are the 2 composition variables to be multiplied
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
!                                   sum = sum+varkappaij+varkappaji+\nu\gamma
! CHECK if \nu\gamma already included in sum ...
!
   integer i,ii,vz,v,w,vv,ternary,ll,lasthope,di,icat,jcat,nnn
   double precision varkappaij,varkappaji,sum,initialij,initialji,nugamma
   double precision xi_ij,xi_ji,sum1,sum2
   logical asymmetric
! added nov 3/2025.  See this date below
! in mixnugamma all vz that have asymmetric ternary with icat or jcat are saved
! because their mixed quad fractions should be added to kvk_ijk
   integer, dimension(:), allocatable :: mixnugamma
!   integer selectij,qz1,qz2,cloop
   integer selectij,qz1,qz2
! mixed update
   integer j,k,l,m,ny,abrakadabra
! If a binary i-j is part of 2 or more asymmetric ternaries i-j-\nu, i-j-\gamma
! the quad fraction x_\nu\gamma should be added to kvk_ijk (the denomonator)
! of kvk_ijk
! saving multiple asymmetrical cations for a binary
   integer, dimension(:), allocatable :: savenu
   integer, dimension(:), allocatable :: savegamma
! debug output
   integer nn1,nn2,nn3,nn4,nn5,nn6,nn7,gg,thisasym
   logical nysym
!   character*3 nyasym
! local variables used for updating quad indices for iasymm, jasymm, etc
!    integer, dimension(:), allocatable :: vk_ij,vk_ji,vk_ijk,xi_ij,xi_ji
! all asymmetric quad indices needed are stored in each separate gtp_allinone
!    integer nvk_ij,nvk_ji,nvk_ijk,nxi_ij,nxi_ji
!
! how to create xquad mm when we need a pointer to gtp_phase_varres?
!   type(gtp_equilibrium_data), pointer :: ceq
! attempt to move mqmqa variables into the mqmqa_var record
!   ceq=>firsteq
!
!   write(*,*)'3XQ line 4970 calls to varkappa1 should be removed'
! varkappa1 is now modified and can be kept
!   write(*,407)1,vz,selectij,tersys(1)%el
!
! Check if y_ik set ...!!!!
   mqf=>phres%mqmqaf
!   write(*,7)'3XQ in varkappa1: ',seq,asymter,size(mqf%y_ik)
7  format(a,3i3)
   if(asymter.ne.0) then
      if(mqmqder) write(*,2)asymter
2     format(/'3XQ in varkappa1, updating asymmetries: ',2i5)
!   else
!      write(*,1)
1     format('3XQ initiating varkappa 1')
   endif
   if(seq.eq.1) then
! this is called for all compvar, for the first remove old savenu and savegamma
      if(allocated(savenu)) deallocate(savenu)
      if(allocated(savegamma)) deallocate(savegamma)
   endif
!   if(mqmqder) &
!        write(*,*)'3XQ line 4462 vk_ij, xi_ij and y_ik with new quad fracs'
!   write(*,*)'3XQ 4558',size(mqf%y_ik),mqmqa_data%ncat
!10 format(a,15(f8.5))
!
!   if(asymter.eq.0) then
! if asymter is 0 there is no need to update expressions for y_ik, vk_ij etc
!      if(seq.eq.0) then
!         write(*,*)'3XQ illegal seq and asymter in varkappa1:',seq,asymter
!         gx%bmperr=4900; goto 1000
!      else
!         write(*,*)'3XB line 4811 create datastructure for varkappa ',seq
!         goto 13
!      endif
! Here all ivk_ij etc have been initiated, we need just calculate vk_ij etc
! That is made from label 600 so we inita
!      cloop=1
!      box=>mqf%compvar(cloop)
!     write(*,*)'3XQ calculate using saved box(i)%ivk_ij data structures',cloop
! This is when box% and box%ivk_ij(... are allocated and fractios calculated
! at label 600 we can calculate values of vk_ij using stored ivk_ij, jvk_ij etc
!      goto 600
!   goto 1000
!----------------------------------------------------------------------   
!
!   write(*,*)'3XQ initiate the data structures for vk_ij etc.',asymter
!
! 2026.04.08: When a ternary asymmetry is changed, all varkappa must be updated
! A ternary asymmetri can be KKK, TKK, KTK or KKT where the asymmetric
! constituent is the first, second or third constituent.
! I do not remember how this is indicated in the loop below   
! But obviously there is some error as KTK and KKT is not registered correctly
! I do not remember how one identifies the asymmetric constituent below
!
!   write(*,10)'3XQ line 3731 y_ik:',(mqf%y_ik(v),v=1,mqmqa_data%ncat)
!
13 continue
  if(.not.allocated(mqf%compvar)) then
     write(*,*)'3XQ line 3076 in varkappa: compvar not allocated, problems'
     gx%bmperr=4399; goto 1000
!   else
!      write(*,*)'3XQ varkappa allocated OK'
   endif
!   write(*,407)2,vz,selectij,tersys(1)%el
!
! seq is a compvar structure for vk_ij and xi_ij
! it has 2 cations and the loop below goes through att ternaries where
! there are ternary asymmetries.  This should be rearranged 
! to loop though all ternaries and set appropriate vk_ij equations.
!
! But that is for later ... on step at a time
!
! Set cloop=0 to prevent looping back for next varkappa
!   cloop=0
! varkappa seems to initiat ivk_ij etc at each calculation! waste of CPU
! content of box will be allocated below using the [ ... ] notation
!   write(*,*)
   box=>mqf%compvar(seq)
! icat and jcat represent cations ... duplicated here (and many other places)
   icat=box%cat1
   jcat=box%cat2
   ia=box%anion
!   write(*,8)seq,size(box%ivk_ij),size(box%jvk_ji),size(box%kvk_ijk)
8  format('3XQ are we here 5496?',i3,2x,3i3)
!   write(*,'(a,4i4,5x,a)')'3XQ varkappa1 line 4267: ',seq,icat,jcat,ia
! the xquad values i,j and j,i are the same but for varkappa they are different
   if(icat.gt.jcat) then
      write(*,3)icat,jcat
3     format(/'In varkappa1: wrong order of elements ',2i4)
      stop
   endif
!--------------------- 
!
   if(allocated(box%ivk_ij)) goto 600
!--------------------- code below set symmetic
!   write(*,30)seq,icat,jcat,ia
30 format('3XQ initiating mii using ijklx, line 4327',i3,2x,3i3)
   mii=ijklx(icat,icat,ia,ia)
!   write(*,*)'3XQ back from ijklx',mii
   mij=ijklx(icat,jcat,ia,ia)
   mjj=ijklx(jcat,jcat,ia,ia)
   if(gx%bmperr.ne.0) then
      write(*,*)'3XQ ijklx index error line 4331'
      stop
   endif
!
   nysym=.false.
! deafult is 0, to update set box%lastupdate to -1
!
!   write(*,407)3,vz,selectij,tersys(1)%el
! below is code to update asymmetry
! and after that the code to calculate varkappa for current molefractions
!   if(box%lastupdate.ne.newXupdate) then
! ignore box%lastupdate as it is not used (yet)
!      if(asymter.eq.0) then
!         write(*,4)seq,box%lastupdate,newXupdate
!4        format('3XQ initiating varkappa newXupdate: ',i3,2x,2i4)
!      else
!         write(*,5)box%lastupdate,newXupdate
!         if(mqmqder) write(*,5)box%lastupdate,newXupdate
!5   format('3XQ *** Updating allinone record from ',i5,' to ',i5,' new: ',2i5)
!      endif
!   endif
!   write(*,*)'3XQ in varkappa1 updates: ',newXupdate,box%lastupdate
!
!   write(*,*)'3XQ in varkappa1 selectij ',newXupdate,box%boxlastupdate
   vzloopupdate:if(newXupdate.gt.box%boxlastupdate) then
!
! Always update! <<<<<<<<<<< box%ivk_ij etc updated separately !!!
!
! *** this if ... endif code part needed only when new asymmetries defined
! Below the arrays below are allocated, the initial 0 is overwritten if used
! This makes use of the new Fortran 2003 facility using [ ]
! Setting an allocatable array to single value means previous values deleted
!      box%ivk_ij=[0]; box%jvk_ji=[0]; box%kvk_ijk=[0]
!
! new asymmetry defined  IGNORE
! this will be handled by testing %asymm variable below (or above)
!      write(*,*)'In varkappa1 new asymmetry in ternary: ',asymter
!      if(asymter.gt.0) then
!         write(*,12)tersys(asymter)%asymm,len(tersys(asymter)%asymm)
12       format('3XQ tersys(%asymter)%asymm "',a,'"',i3)
!         tersys(asymter)%asymm
77       format('3XQ line 5728 varkappa1: ',i2,3x,3i3,3x,3i3,' "',a,'"')
!      endif
! repeating Max equations for vakappa_AB in ternary A-B-C
! -----------if A is asymmetric, \gamma in documentation
! v_AB    x_AA   
! v_BA    x_BB+x_BC+x_CC
! denom=  x_AA+x_BB+x_AB+x_BC+x_AC+x_CC             (=1 if only one ternary)
! ----------- if B is asymmetric, \nu in documentation  
! v_AB    x_AA+x_AC+x_CC
! v_BA    x_BB
! denom=  x_AA+x_BB+x_AB+x_BC+x_AC+x_CC             (=1 if only one ternary)
!------------ if C is asymmetric .........ignore
! if A and B are asymmetric in several ternaries the v_AB and v_BA 
! can include more quadruplets.  One has to update all ternary asymmetries
! at the same time because it is complicated to remove things in the [ ... ]
!
!      if(allocated(savevz)) deallocate(savevz)
      box%boxlastupdate=newXupdate
!      write(*,380)seq,size(box%ivk_ij),size(box%jvk_ji),size(box%kvk_ijk)
380   format('3XQ are we here 5496?',i3,2x,3i3)
! default nyasym is KKK, no asymmetry
!      if(asymter.gt.0) write(*,381)asymter,new_toop
!381   format('3XQ in varkappa1 new asymmetry: ',i2)
! allocate(box%ivk_ij ... done by the elegant [ ... ] statement
! vk derivatives are quad indices, also denominator (same vk_ij and vk_ji)
! the statements below allocate and assign initial quad index
      box%ivk_ij=[mii]; box%jvk_ji=[mjj]; box%kvk_ijk=[mij]
! to simplify handling derivatives the denominator is summed separately
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
! Now take care of asymmetries and update for later use
! Asymmetric vk and xi are updated in the vz loop AND at the end of the loop
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!      if(mqmqder) write(*,*)'3XQ in varkappa1',icat,jcat
!      write(*,*)'3XQ in varkappa1',icat,jcat,size(tersys)
!
! below vz loops through all ternaries ...
! and below that
! 
!      write(*,*)'3XQ with asymmetric cation ',thisasym
! this subroutine is called with the sequentially ordered box%icat,box%jcat
! it must create the basic Kohler model and possibly Toop asymmetries
!
!      write(*,*)'3XQ unfinished varkappa1 code around line 4384'
!      
! The loop below is for all pairs of varkappa records identifying Toop cations
! in ternaries i-j-vz
! and adjusting the expression to calculate varkappa_ij and varkappa_ji
!
! This loop should use the array ternary, vz=1, 2 or 3 testing asymmetries
      vzloop: do vz=1,mqmqa_data%ncat
! loop for all ternary systems to find those with asymmetric i-j-vz and j-i-vz
!         write(*,403)icat,jcat,vz,thisasym
403      format('3XQ in vzloop A: ',2i3,2x,i3,2x,i3,2x,5i3)
! if vz is icat or jcat it is not a ternary
         if(vz.eq.icat .or. vz.eq.jcat) cycle vzloop
! find the sequential order of the ternary icat-jcat-vz 
         ternary=terind(icat,jcat,vz)
! error if ternary not >0
         if(ternary.le.0) goto 1100
! if icat is Toop in this ternary add quadfractions of x_ivz to varkappa_ij
! if jcat is Toop in this ternary add quadfractions of x_jvz to varkappa_ji
! ********* selectij=0 means no asymmetry in this ternary***************
!
!         write(*,119)icat,jcat,vz
119      format('3XQ varkappa1 call test_asymm: ',3i4)
!         selectij=test_asymm(ternary,icat,jcat,vz)
! what we should test is the first 2 letters in tersys(ternary)%asymm
! my head starts to rotate again ....
!         write(*,*)'3XQ line 5554 calling test_asymm'
!
! TOTALLY REMOVED TEST_ASYMM.  Asymmetry set in other ways, not here
         selectij=0
!         selectij=test_asymm(ternary,icat,jcat,vz)
!         write(*,121)'3XQ back from test_asymm:',selectij,ternary,icat,jcat,vz
121       format(a,i2,5x,i3,5x,3i3)
!         if(vz.eq.thisasym) then
!
!         write(*,404)icat,jcat,vz,thisasym,ternary,selectij
404      format('3XQ in vzloop C: ',2i3,2x,i3,2x,i3,2x,5i3)
! ********* selectij=0 means no asymmetry in this ternary***************
! asymm returns 1 if icat is an asymmetric element in icat-jcat-vz  (gamma)
! asymm returns 2 if jcat is an asymmetric element in icat-jcat-vz  (nu)
! asymm returns 3 if both icat and jcat are asymmetric in icat-jcat-vz
! to be considered:  asymmetric i-j-nu and i-j-gamma requires x_\nu\gamma
!                    in the denominator.  For this the savenu/gamma is used
!
! looking for bug that tersys(1)%el(3) is destroyed somewhere
!         write(*,407)5,vz,selectij,tersys(1)%el
407      format('3XQ in varkappa1: ',3i3,', tersys(1)%el: ',3i12)
!         
         if(selectij.eq.0) cycle vzloop
!
!******************** asymmetric ternary *****************************
         write(*,420)selectij,icat,jcat,vz
420      format('3XQ set varkappa ternary asymmetry typ:',i2,' cations: ',3i3)
         asymmetry: select case(selectij)
!
         case default
            write(*,*)'Illegal asymmetry ',selectij
            stop
!-------------------------------------------------------------------
         case(1) ! *************************************************
! icat is asymmetric, save in jvk_ij and in savenu
! an elegant Fortran assignment of an additional items in an allocatable
            write(*,*)'3XQ line 5640 calling ijklx',jcat,vz
            box%jvk_ji=[box%jvk_ji, ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! Below quad fractions added to jvk_ij added to denominator, add ijklx(icat,vz
            box%kvk_ijk=[box%kvk_ijk, ijklx(icat,vz,ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(jcat,vz,ia,ia), &
                 ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! savenu is related to ij, savegamma to ji
            if(allocated(savenu)) then
!               write(*,373)'case 1 use \nu',size(savenu),savenu
373            format('3XQ ',a,' mixed asymmetry terms',i3,': ',10i3)
374            format(a,' x_',2i1)
               do gg=1,size(savenu)
! the mixed terms with \nu should should be added to jvk_ji
                  box%all_ijk=[box%all_ijk, ijklx(vz,savenu(gg),ia,ia)]
!                  write(*,374)'3XQ added ji',savenu(gg),vz
!                  write(*,375)'jvk_ji ',box%jvk_ji
375               format('3XQ ',a,'=',10i4)
               enddo
               savenu=[savenu, vz ]
            else
! otherwize just add vz to savenu
               savenu=[vz]
!               write(*,373)'3XQ line 4377 savednu i ',size(savenu),savenu
            endif
! savegamma is related to ji, maybe add denominator terms
            if(allocated(savegamma)) then
!               write(*,373)'case 1 use \gamma',size(savegamma),savegamma
               do gg=1,size(savegamma)
! the mixed terms with \gamma should should be added to kvk_ijk
                  box%kvk_ijk=[box%kvk_ijk, ijklx(vz,savegamma(gg),ia,ia)]
                  box%all_ijk=[box%all_ijk, ijklx(vz,savegamma(gg),ia,ia)]
!                  write(*,374)'3XQ added kvk_ijk',savegamma(gg),vz
!                  write(*,375)'kvk_ji ',box%kvk_ijk
               enddo
! do not save vz as it does no relates to ij
!               savegamma=[savegamma, vz ]
!            else
! and we must add vz to savegamma
!               savegamma=[vz]
!               write(*,373)'saved i ',size(savevz),savevz
            endif
! The asymmetric xi is depend on y_ik update dxi_ij and dxi_ji
            do nnn=1,mqmqa_data%nquad
!                box%dxi_ij(nnn)=box%dxi_ij(nnn)+dy_ik(icat,nnn)
               box%dxi_ji(nnn)=box%dxi_ji(nnn)+mqmqa_data%dy_ik(vz,nnn)
            enddo
            if(gx%bmperr.ne.0) then
               write(*,*)'3XQ ijklx index error line 4533'
               stop
            endif
!
!---------------------------------------------------------------------
         case(2) ! ***************************************************
! jcat is asymmetric, same as for icat just change icat to jcat!!!!
! and save in jvk_ji ...
            box%ivk_ij=[box%ivk_ij, ijklx(icat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! Nath noted missing  ijklx(vz1,vz2,ia,ia) if icat and jcat are asymmetrical
            box%kvk_ijk=[box%kvk_ijk, ijklx(jcat,vz,ia,ia)]
            box%all_ijk=[box%all_ijk, ijklx(icat,vz,ia,ia), &
                 ijklx(jcat,vz,ia,ia), ijklx(vz,vz,ia,ia)]
! if savegamma allocated we must add terms to jvk_ijk
            if(allocated(savegamma)) then
!               write(*,373)'case 2 use \gamma',size(savegamma),savegamma
               do gg=1,size(savegamma)
                  box%ivk_ij=[box%ivk_ij, ijklx(vz,savegamma(gg),ia,ia)]
                  box%all_ijk=[box%all_ijk, ijklx(vz,savegamma(gg),ia,ia)]
!                  write(*,374)'3XQ added ij',savegamma(gg),vz
!                  write(*,375)'ivk_ij ',box%ivk_ij
               enddo
               savegamma=[savegamma, vz ]
            else
! and we must add vz to savevz
               savegamma=[ vz ]
!               write(*,373)'savedgamma j ',size(savegamma),savegamma
            endif
! savenu is related to ij, maybe add denominator terms
            if(allocated(savenu)) then
!               write(*,373)'case 2 use \nu',size(savenu),savenu
               do gg=1,size(savegamma)
! the mixed terms with \nu should should be added to kvk_ijk
                  box%kvk_ijk=[box%kvk_ijk, ijklx(vz,savenu(gg),ia,ia)]
                  box%all_ijk=[box%all_ijk, ijklx(vz,savenu(gg),ia,ia)]
!                  write(*,374)'3XQ added kvk_ijk',savenu(gg),vz
!                  write(*,375)'jvk_ji ',box%kvk_ijk
               enddo
            endif
! The asymmetric xi is depend on y_ik update dxi_ij and dxi_ji
            do nnn=1,mqmqa_data%nquad
!                box%dxi_ij(nnn)=box%dxi_ij(nnn)+dy_ik(jcat,nnn)
               box%dxi_ij(nnn)=box%dxi_ij(nnn)+mqmqa_data%dy_ik(vz,nnn)
            enddo
            if(gx%bmperr.ne.0) then
               write(*,*)'3XQ ijklx index error line 4578'
               stop
            endif
!
!---------------------------------------------------------------------
         case(3) ! **************************************************
! Both icat and jcat are asymmetric NOT IMPLEMENTED
            write(*,788)icat,jcat,vz
788         format('3XQ *** Illegal with 2 asymmetric cations ',2i3,' with ',i3)
            gx%bmperr=4399; goto 1000
! tentative code below
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
            if(gx%bmperr.ne.0) then
               write(*,*)'3XQ ijklx index error line 4606'
               stop
            endif
!
         end select asymmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         write(*,778)icat,jcat,vz
778      format('3XQ asymmetry set ',3i3,' box%all: ',10i3)
         cycle vzloop
!         goto 747
!
! loops below now redundant when we added savevz loops above ..... ????
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
!                              write(*,806)i,k,m,ijklx(k,m,ia,ia)
!                              write(*,805)'kvk_ijk ',box%kvk_ijk
                           endif
                        endif
                     enddo neverending
                     if(gx%bmperr.ne.0) then
                        write(*,*)'3XQ ijklx index error line 4645'
                        stop
                     endif
                  enddo
               endif
            enddo
         enddo addkvkterm
805 format(a,20i3)
806      format('3XQ adding mixed quad to kvk_ijk',i3,2x,2i3,2x,i3)
! end copied code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!747      continue
      enddo vzloop
! the vzloop above should be done whenever the asymmetry changes
!
!      write(*,748)box%cat1,box%cat2
748   format('3XQ Asymmetry updated for varkappa_ij: ',2i3)
!--------------------------------------------------------------------
! end of asymmetry detection loop
!--------------------------------------------------------------------
!
!   else
!      write(*,*)'3XQ using current asymmetry'
!
   endif vzloopupdate
!
! Code below is to set values in vk_ij, xi_ij and y_j/k from x_ij
600 continue
!
!--------------------------------------------------------------------
! Below arrays box%ivk_ij, box%jvk_ji, box%dxi_ij are used to
! calculate \varkappa and \xi and their derivatives
!--------------------------------------------------------------------
!
! Now use the structures ivk_ij, jvk_ji, kvk_ijk and dxi_ij, dxi_ji
!   write(*,*)'3XQ in varkappa1 line 3900',allocated(box%ivk_ij),&
!        allocated(box%dvk_ij)
! if seq nonzero we have just set indices for asymmetry etc in %ivk_ij etc
! if seq=0 we must loop for cloop box using saved %ivk_ij indices
! incrementing cloop and return to label 600 until all done
   varkappaij=0.0d0; varkappaji=0.0d0; sum=0.0d0; nugamma=0.0d0
   do ii=1,size(box%ivk_ij)
      varkappaij=varkappaij+mqf%xquad(box%ivk_ij(ii))
!       write(*,697)'ivk_ij',ii,box%ivk_ij(ii),varkappaij,xquad(box%ivk_ij(ii))
697   format('Summing ',a,': ',2i3,2(1pe14.6))
   enddo
603 format('Partial sum: ',i3,a,1pe12.4,' quad: ',5i3)
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
696   format('Quad indices in',a,': ',20i4)
   sum=sum+varkappaij+varkappaji+nugamma
!    write(*,601)sum,nugamma
601   format('Total value      Denominator: ',1pe12.4,' nugamma: ',1pe12.4)
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
! Here the asymmetric contributions to vk should be added
!   write(*,666)seq,size(box%ivk_ij),size(box%jvk_ji),size(box%kvk_ijk),&
!        box%vk_ij,box%vk_ji
666 format('3XQ line 5834 vk_ij:',i2,3i3,2f10.6)
!
! DERIVATIVES MISSING?
!
! Values of xi_ij and y_ik can be calculated without looping
888 continue
!
! Calculating all mqf%y_ik values from xquad
   do v=1,mqmqa_data%ncat
      mqf%y_ik(v)=0.0d0
!      write(*,20)'3XQ dy_ik',(mqmqa_data%dy_ik(v,w),w=1,mqmqa_data%nquad)
!      write(*,20)'3XQ dy_ik',(mqf%dy_ik(v,w),w=1,mqmqa_data%nquad)
20    format(a,(20F5.2))
      do w=1,mqmqa_data%nquad
         mqf%y_ik(v)=mqf%y_ik(v)+mqmqa_data%dy_ik(v,w)*mqf%xquad(w)
      enddo
   enddo
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! this code use the updated data structure to calculate quickly
! This should be called by set constitution!!
!   write(*,*)'In varkappa1 calling dexcess_dq to allocate and set %dvk_ij?'
   goto 900
!
500 continue    
!
900 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   if(box%boxlastupdate.ne.newXupdate) then
      box%boxlastupdate=newXupdate
!      write(*,1001)box%seq,box%boxlastupdate
1001  format('3XQ allinone record ',i3,' updated to new asymmetries ',i5)
   endif
!
1000 continue
   if(mqmqder) write(*,*)'3XQ Leaving varkappa1'
!   write(*,*)'3XQ Leaving varkappa1'
   return
!
1100 continue
   write(*,1105)icat,jcat,v
1105 format('Error return from tersym for elements: ',3i4)
   goto 1000
 end subroutine varkappa1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine set_ternary_asymmetry_tdb
!\begin{verbatim}
 subroutine set_ternary_asymmetry_tdb(line)
   implicit none
   character*(*) line
! The phase name is extracted as first text on line
! This is used for TYPE_DEF when asymmetries are read from a TDB file 
! to set asymmetries in a text
!\end{verbatim}
!
   type(gtp_equilibrium_data), pointer :: ceq
   type(gtp_phase_varres), pointer :: phres
!   type(gtp_mqmqa_var), pointer :: mqf
!
   integer i,j,ip,iq,ia,ib,ic,mm,icc(3),nc,kk,vz,toop(3),elimq,spix,i3
   integer missasym,asymter
   integer iph,ics,icon,ipm,length,pcon,lokcs,iv1,qq,jj,selected
   double precision mass
   character phase*24
!
   ip=0
   call getext(line,ip,2,phase,' ',iq)
!   write(*,5)trim(phase)
5  format(/'3XQ In set_ternary_asymmetry_tdb called from gtp3E for ',a)
!   if(mqmqdebug) write(*,10)trim(line)
!   write(*,*)'3E set_ternary_asymmetry to be fixed'
! extract constituent indices and call setsym'
! text is extacted from frist nonblank position ip to first space
! first the phase, then 3 constituents finally the asymcode
!
!   write(*,20)trim(phase),iq,ip
20 format('3XQ Phase name: ',a,5x,i3)
   if(phase(1:1).ne.' ') then
      call find_phase_by_name(phase,iph,ics)
      if(gx%bmperr.ne.0) then
         write(*,21)trim(phase)
21       format('3XQ Ternary asymmetries for phase "',a,&
              '" ignored as phase not selected')
         gx%bmperr=0
         goto 1000
      endif
!      write(*,*)'3XQ Found phase ',iph
      if(phase(1:4).ne.'MSCL') then
         write(*,*)'3XQ skipping asymmetries of phase ',phase
         goto 1000
      endif
   else
      write(*,*)'3XQ No phase name in set_ternary_asymmetry, skipping'
      goto 1000
   endif
!
!10 format('3XQ line 5934 gtp3E calls set_ternary_asymmetry for phase: ',a)
! from iph we must find phres ... iph same as lokph?
!
!   length=len_trim(line)
!   write(*,22)trim(phase),iph,length
22 format('3XQ Ternary asymmetries for ',a,' must be set manually.',2i5)
! we have to extract the 3 constituents for each asymmetry and check
! if they are entered
!   write(*,31)line(ip:length)
31 format('3XQ from TDB: ',a)
!
! double sigh
! phres not allocated yet ............................... ????
!   write(*,*)'3XQ looking for phase link',iph
!   lokcs=phasetuple(iph)%lokvares
   ceq=>firsteq
! safer way to find phres ....
   i=1
7  continue
      i=i+1
!      write(*,*)'3XQ looking for phres',i,iph
      phres=>ceq%phase_varres(i)
!      write(*,*)'3XQ looking for phres',i,iph,phres%phlink
      if(phres%phlink.ne.iph) goto 7
!      write(*,*)'3XQ found phase ',iph
!----------------------------------------------------      
!>>>>>>>>>>>>>>>>>> crash here
!   phres=>ceq%phase_varres(lokcs)
!   write(*,*)'3XQ found phres '
!   mqf=>phres%mqmqaf
!   write(*,*)'3XQ found mqf '
!
! what we can do is to save the line with asymmetries to be processed later
!   write(*,*)'3XQ store the asymmetry line in the mqmqa_data%tdbasymmetries'
! the phase name is first in the line and it is terminated by "!"
! maybe several phases will be appended?
! SKIP THIS
!   if(allocated(mqmqa_data%tdbasymmetries)) then
!      write(*,*)'3XQ tdbasymmetries already allocated, IGNORE!'
!   else
!      mqmqa_data%tdbasymmetries=trim(line)//' !'
!      write(*,41)mqmqa_data%tdbasymmetries
41    format('3XQ saved: "',a,'" in mqmqa_data%tdbasymmetries')
!   endif
!
  i=len_trim(line)  
!  write(*,*)'3XQ call set_tdbasymmetries'
! skip phase name in list
   call set_tdbasymmetries(phres,line(ip:i))
!
1000 continue
   return
 end subroutine set_ternary_asymmetry_tdb

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine set_tdbasymmetries
!\begin{verbatim}
 subroutine set_tdbasymmetries(phres,line)
! hande tdb asymmetries stored in mqmqa_data%tdbasymmetries
   implicit none
   type(gtp_phase_varres), pointer :: phres
   character line*(*)
!\end{verbatim}
   character asymcode*6,ocasym*3,const*24,missingcon*24,pconid*24,ch1*1
   type(gtp_mqmqa_var), pointer :: mqf
   integer ip,jj,qq,asp,nc,pcon,spord(3),selected,toop,length,kkp,asymter
   logical verbose
!
!   write(*,*)'3XQ in set_tdbasymmetries extracting TDB asymmetries'
!
   mqf=>phres%mqmqaf
   if(.not.allocated(mqf%names_y_ik)) then
      write(*,*)'3XQ names_y_ik not allocated'
      stop
   endif
!------------------------------------------------------------

! there must be something which must be initiated before this routine ???
!   verbose=.false.
   verbose=mqverbose
!
! it seems this routine should inititae symmetric compvar
   call init_symmetric_varkappa(phres,verbose) 
!
!   write(*,*)'3XQ listing of vk_ij, vk_ji'
!   call list_compvar(phres)

!------------------------------------------------------------
   length=len_trim(line)
!   write(*,10)length,trim(line)
10 format(/'3XQ line length ',i4,' from tdbfile with asymmetries: "',a,'"')
!
!...extract from line sequences of 3 constituents followed by an asymmetry code
! the names of the "real" x_ii constituents are in mqf%names_y_ik !!!
!   write(*,88)(trim(mqf%names_y_ik(jj)),jj=1,size(mqf%names_y_ik))
88 format('3XQ const: ',10(a,2x))
   ip=1
   extractasym: do while(ip.lt.length)
! typically  MG/CL NA/CL CE/CL T1T1K ..... loop 3 times for species
      kkp=ip
      spord=0
      selected=0
!      write(*,90)kkp
90    format(/'3XQ loop to extract 3 quads with asymmetry ',i4)
      searchfor3: do asp=1,3
         if(eolch(line,ip)) exit extractasym
         if(line(ip:ip).eq.'!') exit extractasym
! extract name of species from line
         call skip_2_name(line,ip,const,1,ch1)
         if(const(1:1).eq.'!') then
            write(*,*)'3XQ found terminating character'
            exit extractasym
         endif
!         write(*,91)ip,const
91       format('3XQ position ',i5,' extracted quadname: "',a,'"')
         nameok: do jj=1,size(mqf%names_y_ik)
            qq=index(mqf%names_y_ik(jj),'-Q')-1
!            write(*,93)const(1:qq),mqf%names_y_ik(jj)(1:qq)
93          format('3XQ names: "',a,'" and "',a,'"')
            if(mqf%names_y_ik(jj)(1:qq).eq.const(1:qq)) then
               selected=selected+1
               spord(selected)=jj
!               write(*,96)asp,spord
96             format('3XQ found quad',i3,5x,3i3)
               cycle searchfor3
            endif
         enddo nameok
      enddo searchfor3
! always extract asymmetry code, can be T1T1K or TKK, take it letter by letter
      call skip_2_name(line,ip,asymcode,1,ch1)
!      write(*,*)'3XQ asymcode: "',asymcode,'"',ip
      if(selected.ne.3) then
!         write(*,222)line(kkp:ip),ip,spord
222      format('3XQ unselected: "',a,'" ',i4,3x,3i3)
      else
! decode asymmetry only if all 3 constituents selected
         if(len_trim(asymcode).eq.3) then
            ocasym=asymcode(1:3)
            toop=index(ocasym,'T')
         else
            jj=index(asymcode,'T')
            if(asymcode(jj+1:jj+1).eq.'1') then
               ocasym='TKK'; toop=1
            elseif(asymcode(jj+1:jj+1).eq.'2') then
               ocasym='KTK'; toop=2
            elseif(asymcode(jj+1:jj+1).eq.'3') then
               ocasym='KKT'; toop=3
            endif
         endif
! if selected=3 then all 3 asymmetric quads are constituents
! we have found 3 species names and checked if they are constituents
         write(*,223)line(kkp:ip),ip,spord,toop,ocasym
223      format('3XQ Asymmetri by TDB: "',a,'" ',i4,3x,3i3,3x,i2,2x,a)
!         write(*,38)selected,spord,asymcode
38       format('3XQ skipping selected ternary ',i1,5x,3i3,5x,a)
! make sure spord(1..3) are in increasing order (bubbelsort)
!
         if(spord(1).lt.spord(2) .and. &
              spord(2).lt.spord(3)) then
            if(mqverbose) write(*,*)'3XQ asymmetric quads are in order'
         else
            write(*,*)'3XQ problem: the asymmetrical quads in not in order'
            stop
         endif
! we have to find this ternary in tersys
! tersys has the ternaries ordered 
!  1 2 3;  1 2 4;  1 3 4; 2 3 4;  with 1, 2, 3, 4 in the order of names_y_ik
! loop tersys until we find the one with spord(1..3)
         tloop: do jj=1,size(tersys)
!            write(*,300)jj,tersys(jj)%el
300         format('3XQ line 6172 tersys ',i3,', cations: ',3i3)
            if(tersys(jj)%el(1).eq.spord(1) .and. &
                 tersys(jj)%el(2).eq.spord(2) .and. &
                 tersys(jj)%el(3).eq.spord(3)) then
               asymter=jj
            endif
         enddo tloop
! this is the asymmetrical ternary
!         write(*,*)'3XQ set ternary asymmetric: ',asymter,toop
! Something must be initiated before we can add a ternary asymmetry, what?
! the call to new_ternary_asym has asymter as index of ternary
         if(verbose) write(*,*)'3XQ calling new_ternary_asym'
         call new_ternary_asym(asymter,toop,phres,verbose)
         if(verbose) write(*,*)'3XQ back from new_ternary_asym'
      endif
! continue extract ternary asymmetries
!      write(*,*)'3XQ at position ',ip,' in line'
!      read(*,101)ch1
101   format(a)
   enddo extractasym
!
!
   goto 1000
!   
1000 continue
   if(mqverbose) write(*,*)'3XQ Exit from set_tdbasymmetries'
   return
 end subroutine set_tdbasymmetries

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine skip_2_name
!\begin{verbatim}
 subroutine skip_2_name(line,ip,name,mode,ch1)
! extract name after skipping spaces
   implicit none
   character line*(*),name*(*),ch1*1
   integer ip,mode
!   logical eolch
!\end{verbatim}
! skip to finst nonspace 
   if(eolch(line,ip)) then
      name=' '
   else
      call getname(line,ip,name,mode,ch1)
   endif
900 return
 end subroutine skip_2_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine correlate_const_and_quads
!\begin{verbatim}
! subroutine correlate_const_and_quads(loksp)
 subroutine correlate_const_and_quads(lokph)
! this subroutine should for each mqmqa constituent create their data
! called from gtp3B.F90.  NOT NEEDED FOR CALCULATIONS
! lokph is index of mqmqa phase record
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
   integer cat1,cat2,kk,val2
   integer missing,ll,nocon,s1,iv1,iv2,iv3
   logical noanion
   integer, allocatable, dimension(:) :: invert,inverse
   integer, allocatable, dimension(:) :: findan
   character quadname*24,elname*4,elval*4
   character*4, allocatable, dimension(:) :: catnames
   integer, allocatable, dimension(:) :: catindex
   integer, allocatable, dimension(:) :: multipleval
   integer multival,noofcations
!
! called from create_asymmetry in gtp3B
   nfr=phlista(lokph)%nooffr(1)
!   write(*,7)lokph,nfr,noofel
7  format(/'3XQ In correlate_const_and_quad',3i5/)
!
   allocate(findan(noofel))
!   lokcs=phlista(lokph)%linktocs(1) composition set?
! note element numbers are not in order, the anion may be anywhere
!
! first step, find the anion, it is present in all constituents
! Stupid to do this here, it has already been found but lost
   findan=0
   do isp=1,nfr
      loksp=phlista(lokph)%constitlist(isp)
      iel=size(splista(loksp)%ellinks)
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
! 2026-05-10 redesign: identify the anion via the first AA/XX quad's
! contyp(14,k) (which after the gtp3B reorder pass holds the anion species'
! splista index).  This is robust when several elements tie on findan
! (e.g. U-Cl with multiple U valencies has findan(U)==findan(CL)).
   el1=0; el2=0
   do isp=1,nfr
      if(mqmqa_data%contyp(5,isp).gt.0) then
         iv1=mqmqa_data%contyp(14,isp)
         if(iv1.ge.1 .and. iv1.le.noofsp) then
            el2=splista(iv1)%ellinks(1)
            exit
         endif
      endif
   enddo
   if(el2.le.0) then
! fallback to the count-based detection if no AA/XX found
      do jp=1,noofel
         if(findan(jp).gt.el1) then
            el1=findan(jp); el2=jp;
         endif
      enddo
   endif
! Maybe there are elements not dissolved in MQMQA as He or Ar?
!   write(*,*)'3XQ anion species index: ',mqmqa_data%anionspix
!
!   write(*,*)'3XQ multivalent? ',noofel,(noofel-1)*(noofel-2),mqmqa_data%nconst
   mqmqa_data%xanione=el2
   mqmqa_data%xanionalpha=ellista(el2)%alphaindex
!   write(*,6)ellista(mqmqa_data%xanione)%symbol,&
!        mqmqa_data%xanione,mqmqa_data%xanionalpha
6  format(/'3XQ line 3383 anion: ',a,' ellink: ',i3,' alphabetically: ',i3/)
!
! 2026-05-10 redesign: cation indexing now by species, not element.
! Build cat2species (cation_idx -> splista_idx), cat2el (cation_idx -> element
! alpha-idx), sp2cat (splista_idx -> cation_idx, neg if anion species),
! el2ancat (element alpha-idx -> first cation_idx for that element, neg if
! anion element, 0 if neither).  All paths populate the same arrays so
! downstream code (entropy, ternary asymmetry, list_quads) does not have to
! know whether the system has multivalent elements.
!
! The cation alphabetical order (1..ncat) is the order encoded in contyp(5,isp)
! for AA/XX quads and in contyp(6,isp), contyp(7,isp) for cross quads.  After
! the gtp3B alphabetical reorder pass, the splista index of the cation lives
! in contyp(13,isp) for AA/XX and contyp(11,isp), contyp(12,isp) for cross.
!
! Step 1: walk constituents and pick the splista index of each cation from
! its AA/XX (pure) quad.  Every cation in MQMQA must have a pure quad so
! this is sufficient.  For AA/XX, contyp(5,isp) is the cation alpha-index
! (1..ncat) and contyp(13,isp) is the cation splista index.  Cross quads
! are skipped here because gtp3B's "replace" pass has overwritten
! contyp(11,12) with sublattice indices, not splista indices.
   allocate(catindex(2*nfr))
   catindex=0
   do isp=1,nfr
      if(mqmqa_data%contyp(5,isp).gt.0) then
         iv1=mqmqa_data%contyp(5,isp)
         if(iv1.ge.1 .and. iv1.le.size(catindex)) &
              catindex(iv1)=mqmqa_data%contyp(13,isp)
      endif
   enddo
! Step 2: count contiguous cations.
   noofcations=0
   do iv1=1,size(catindex)
      if(catindex(iv1).gt.0) noofcations=iv1
   enddo
   if(noofcations.le.0) then
      write(*,*)'3XQ ERROR: no cations found from contyp'
      gx%bmperr=4399
      goto 1000
   endif
! Step 3: allocate and fill cat2species, cat2el.
   if(allocated(mqmqa_data%cat2species)) deallocate(mqmqa_data%cat2species)
   if(allocated(mqmqa_data%cat2el)) deallocate(mqmqa_data%cat2el)
   allocate(mqmqa_data%cat2species(noofcations))
   allocate(mqmqa_data%cat2el(noofcations))
   do iv1=1,noofcations
      mqmqa_data%cat2species(iv1)=catindex(iv1)
! cation species has a single non-anion element
   mqmqa_data%cat2el(iv1)=ellista(splista(catindex(iv1))%ellinks(1))%alphaindex
   enddo
   deallocate(catindex)
! Step 4: el2ancat sized noofel (element alpha-idx -> first cat or anion mark).
   if(allocated(mqmqa_data%el2ancat)) deallocate(mqmqa_data%el2ancat)
   allocate(mqmqa_data%el2ancat(noofel))
   mqmqa_data%el2ancat=0
   do jp=1,noofel
      if(jp.eq.mqmqa_data%xanionalpha) then
         mqmqa_data%el2ancat(jp)=-jp
      else
         do iv1=1,noofcations
            if(mqmqa_data%cat2el(iv1).eq.jp) then
               mqmqa_data%el2ancat(jp)=iv1
               exit
            endif
         enddo
      endif
   enddo
! Step 5: sp2cat sized noofsp (species-> cation idx, -1 for anion species).
   if(allocated(mqmqa_data%sp2cat)) deallocate(mqmqa_data%sp2cat)
   allocate(mqmqa_data%sp2cat(noofsp))
   mqmqa_data%sp2cat=0
   do iv1=1,noofcations
      if(mqmqa_data%cat2species(iv1).ge.1 .and. &
           mqmqa_data%cat2species(iv1).le.noofsp) then
         mqmqa_data%sp2cat(mqmqa_data%cat2species(iv1))=iv1
!         write(*,604)iv1,mqmqa_data%cat2species(iv1)
      endif
   enddo
   if(mqmqa_data%anionspix.ge.1 .and. mqmqa_data%anionspix.le.noofsp) then
      mqmqa_data%sp2cat(mqmqa_data%anionspix)=-1
   endif
!
   if(mqmqtdb) then
      write(*,604)noofcations,(mqmqa_data%cat2species(iv1),iv1=1,noofcations)
604   format('3XQ cat2species: ',i3,5x,20(i3,1x))
      write(*,605)noofcations,(mqmqa_data%cat2el(iv1),iv1=1,noofcations)
605   format('3XQ cat2el     : ',i3,5x,20(i3,1x))
      write(*,606)noofel,(mqmqa_data%el2ancat(jp),jp=1,noofel)
606   format('3XQ el2ancat   : ',i3,5x,20(i3,1x))
!
      write(*,16)'3XQ Elements alphabetically:  ',&
           ((ellista(elements(jp))%symbol),jp=1,noofel)
16    format(a,20(1x,a2))
   endif
!
! Step 6: build con2quad using cation indices already in contyp(5..7).
! contyp(5,isp) > 0 means AA/XX (single cation, doubled); contyp(5)=0 means
! cross quad with cation alphas in (6,7).  This logic is independent of how
! many distinct cations any one element has, so it is correct for both
! mono- and multi-valent systems.
   if(allocated(mqmqa_data%con2quad)) deallocate(mqmqa_data%con2quad)
   allocate(mqmqa_data%con2quad(nfr))
   con2quadvalency: do isp=1,nfr
      try1: if(mqmqa_data%contyp(5,isp).gt.0) then
         cat1=mqmqa_data%contyp(5,isp)
         cat2=cat1
      else
         cat1=mqmqa_data%contyp(6,isp)
         cat2=mqmqa_data%contyp(7,isp)
      endif try1
      mqmqa_data%con2quad(isp)=ijklx(cat1,cat2,1,1)
      if(gx%bmperr.ne.0) then
         write(*,47)cat1,cat2
47       format('3XQ element index wrong',2i3)
         write(*,*)'3XQ ijklx index error',cat1,cat2
      endif
   enddo con2quadvalency
!   write(*,711)mqmqa_data%con2quad
711 format('3XQ con2quad: ',20(i3))
!
! common code with or without valencies
!
888 continue
! allocate also array with all A/X quads
   allocate(mqmqa_data%emquad(mqmqa_data%ncat))
! enter data in emquad
   cat1=1
   cat2=mqmqa_data%ncat
   do isp=1,mqmqa_data%ncat
      mqmqa_data%emquad(isp)=cat1; cat1=cat1+cat2; cat2=cat2-1
   enddo
! list quads (why?)
!   write(*,68)(mm,mm=1,mqmqa_data%nquad)
68 format('3XQ quads:    ',21i3)
!   write(*,57)'3XQ emquads: ',(mqmqa_data%emquad(isp),isp=1,mqmqa_data%ncat)
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
!   write(*,70)(trim(splista(phlista(lokph)%&
!        constitlist(mqmqa_data%con2quad(jp)))%symbol),jp=1,nfr)
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
!   write(*,*)'3XQ leaving correlate_const_and_quads',lokph
   return
 end subroutine correlate_const_and_quads

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine list_quads(kk)
!\begin{verbatim}
 subroutine list_quads(kk)
! emergency subroutine because phlista protected in pmon6
! 2026-05-10 redesign: now multivalence-aware; iterates cat2species so an
! element with several cation species (e.g. U+ and U2 dimer) lists each.
   implicit none
   integer kk
!\end{verbatim}
   integer nel,iv1,iv2
!
   kk=0
   write(*,2)(ellista(elements(nel))%symbol,nel=1,noofel)   
2  format(/'3XQ In list_quads: Element names:      ',20(a2,1x))
! anion alpha-index has el2ancat < 0
   do nel=1,noofel
      if(mqmqa_data%el2ancat(nel).lt.0) kk=nel
   enddo
   if(kk.eq.0) then
      write(*,*)'You have a strange MQMQA system without any anion'
      goto 1000
   else
      write(*,4)ellista(elements(kk))%symbol,mqmqa_data%xanionalpha,&
           mqmqa_data%xanione
4     format('The anion element name, index and link: ',a,2i3)
   endif
   write(*,3)size(mqmqa_data%el2ancat),mqmqa_data%el2ancat
3  format('3XQ el2ancat (element-alpha -> first cation idx):',i3,2x,20i3)
   if(allocated(mqmqa_data%cat2species)) then
      write(*,5)size(mqmqa_data%cat2species), &
           (splista(mqmqa_data%cat2species(iv1))%symbol, &
            mqmqa_data%cat2species(iv1), &
            mqmqa_data%cat2el(iv1), iv1=1,size(mqmqa_data%cat2species))
5     format('3XQ cation species cat -> symbol(splista/element-alpha)',i3,':'/&
           20(a4,'(',i3,'/',i2,') '))
   endif
! one line per cation
   if(allocated(mqmqa_data%cat2species)) then
      write(*,10)
10    format(/'Cati Name',10x,'Species index')
      do iv1=1,size(mqmqa_data%cat2species)
! ONLY ENDMEMBER CATIONS, NO MIXED CATIONS LISTED HERE
! is cat2species really correct for multivalent cations as U
         write(*,20)iv1,splista(mqmqa_data%cat2species(iv1))%symbol,&
              mqmqa_data%cat2species(iv1)
20       format(i3,2x,a,2x,a)
      enddo
   endif
!
1000 continue
   return
 end subroutine list_quads

 !/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine list_quads_with_single_cation
!\begin{verbatim}
 subroutine list_quads_with_single_cation(phres)
! emergency subroutine because phlista protected in pmon6
! 2026-05-10 redesign: now multivalence-aware; iterates cat2species so an
! element with several cation species (e.g. U+ and U2 dimer) lists each.
   implicit none
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   type(gtp_mqmqa_var), pointer :: mqf   
!   type(gtp_allinone), pointer :: box
   integer nel,iv1,iv2,ip,emq,kk,i0

   character line1*300
   mqf=>phres%mqmqaf
   i0=ichar('0')
!
   write(*,10)
10 format('List of all cations also as species'/&
        'Cati Quad "no -Qij"    Quad index  Quad sequential  Species index')
   emq=0
   kk=1
! I may have to resort to list the phase constituents ...........
   loop1: do iv1=1,size(mqf%names_y_ik)
      write(*,20)iv1,mqf%names_y_ik(iv1),'x_'//char(i0+iv1)//char(i0+iv1),&
           mqf%spqx_y_ik(iv1),mqf%spix_y_ik(iv1)
20    format(i3,2x,a,a,2i16)
   enddo loop1
1000 continue
   return
 end subroutine list_quads_with_single_cation

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine listconst
!\begin{verbatim}
 subroutine listconst(iph)
! emergency subroutine because phlista protected in pmon6
! 2026-05-10 redesign: cation indices for each quad come from contyp(5..7)
! (already encoded as 1..ncat after the alphabetical reorder pass), not via
! ellinks->el2ancat which collapses two cations sharing one element.
   implicit none
   integer iph
!\end{verbatim}
   type(gtp_phase_varres), pointer :: phres
   integer lokph,lokcs,isp,iel,elx(4),elxx(4),jp,j4,nel,cations(2),jj,kk
! 2026-05-10: size elsym to at least 4 so the noofel<=3 print loop (which
! always indexes elsym(1..3)) doesn't go out of bounds when noofel<3.
   character :: elsym(max(4,noofel))*2

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
   write(kou,5)
5  format(/'Con  Quad Nel Elements      Elem index',2x,'Species name',&
        15x,'Cations')
   specie: do jp=1,phlista(lokph)%nooffr(1)
      isp=isp+1
      j4=phlista(lokph)%constitlist(jp)
      nel=size(splista(j4)%ellinks)
      elsym='  '
      elxx=1000
      element: do iel=1,nel
         elx(iel)=splista(j4)%ellinks(iel)
         elsym(iel)=ellista(elx(iel))%symbol
         elxx(iel)=ellista(splista(j4)%ellinks(iel))%alphaindex
      enddo element
! Cation indices in cation-alpha space (1..ncat) come straight from contyp:
!   col 5 > 0: AA/XX, both cation slots are this index
!   col 5 = 0: cross, the two cation alphas are in cols 6 and 7
      cations=0
      if(mqmqa_data%contyp(5,isp).gt.0) then
         cations(1)=mqmqa_data%contyp(5,isp)
         cations(2)=cations(1)
      else
         cations(1)=mqmqa_data%contyp(6,isp)
         cations(2)=mqmqa_data%contyp(7,isp)
      endif
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

!addtotable subroutine quadprops
!\begin(verbatim} subroutine quadprops
! subroutine quadprops(lokcs,spix)
 subroutine quadprops(phvar,spix)
! subroutine to list mqmqa_data%contyp used in pmon6
!   integer lokcs
   integer, dimension(*) :: spix
   type(gtp_phase_varres), pointer :: phvar
!\end{verbatim}
!   type(gtp_allinone), pointer :: box
!   character*300 line1,line2,qline
!   character*2, dimension(:), allocatable :: quadcat
   integer i1,i2,loksp,loksp2,loksp3,loksp4
! lokph ?
! splista(phlista(lokph)%constitlist(nk))%symbol
! unfinished here
   write(*,*)'3XQ quadprops unfinished'
!   phvar=>ceq%phase_varres(lokcs)
   do i2=1,mqmqa_data%nconst
!      loksp=phlista(lokph)%constitlist(i2)
!      loksp2=mqmqa_data%con2quad(i2)
!      loksp3=mqmqa_data%emquad(i2)
!      loksp4=mqmqa_data%quad2compvar(i2)
!
      write(kou,3433)i2,(mqmqa_data%contyp(i1,i2),i1=1,14),&
           trim(mqmqa_data%quadlist(i2))
! to be added mqmqa_data%constoi(i1,i2),i1=1,4 with format 4F5.1
!             phvar%yfr(i2)                    with format F5.2
!      loksp,loksp2,loksp3,loksp4,
3433  format('Quad ',i3,': ',4i3,1x,i4,1x,4i3,1x,i3,1x,4i3,1x,a)
!                            4      1     4      1     4        spix
   enddo
   write(*,*)mqmqa_data%emquad
   return
 end subroutine quadprops

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!addtotable subroutine list_ternary_cations
!\begin(verbatim} subroutine list_ternary_cations
 subroutine list_ternary_cations(asymter,new_toop,phres)
! subroutine to list \varkappa, \xi and y_ik definitions
   implicit none
   integer asymter,new_toop
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   integer i1,i2,i3
   type(gtp_mqmqa_var), pointer :: mqf
!
   mqf=>phres%mqmqaf
!   write(*,*)'3XQ line 6934 asymter: ',asymter
!   write(*,10)(trim(mqf%names_y_ik(i1)),i1=1,size(mqf%names_y_ik))
10 format('3XQ line 6906 Basic quads: ',10(a,' '))
   if(asymter.eq.0) then
      write(*,11)
11    format(/'3XQ Ternary Quadruplets',36x,'Quadruplet indices')
      do i1=1,size(tersys)
         write(*,20)i1,trim(mqf%names_y_ik(tersys(i1)%el(1))),&
              trim(mqf%names_y_ik(tersys(i1)%el(2))),&
              trim(mqf%names_y_ik(tersys(i1)%el(3))),&
              mqf%spqx_y_ik(tersys(i1)%el(1)),&
              mqf%spqx_y_ik(tersys(i1)%el(2)),&
              mqf%spqx_y_ik(tersys(i1)%el(3))
!              mqmqa_data%el2quad(tersys(i1)%el(1)),&
!              mqmqa_data%el2quad(tersys(i1)%el(2)),&
!              mqmqa_data%el2quad(tersys(i1)%el(3))
20       format(i7,',   1:',a,'   2: ',a,'   3: ',a,20x,3i3)
      enddo
   else
! list just one ternary with new asymmetry
      i1=asymter
         write(*,20)i1,trim(mqf%names_y_ik(tersys(i1)%el(1))),&
              trim(mqf%names_y_ik(tersys(i1)%el(2))),&
              trim(mqf%names_y_ik(tersys(i1)%el(3)))
   endif
1000 continue
   return
 end subroutine list_ternary_cations

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!addtotable subroutine list_mqmqa_variables
!\begin(verbatim} subroutine varkappadefs
 subroutine list_mqmqa_variables(phres)
! subroutine to list \varkappa, \xi and y_ik definitions
   implicit none
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   integer nv,cat1,cat2,i,j,k,ip,i0,jp,kp,nv2,ix,jx,iz,j4,lenhead
   type(gtp_allinone), pointer :: box
   character*300 line1,line2,qline,xijheader
   character*2, dimension(:), allocatable :: quadcat
   character*3 fact(0:2),coeff
   character*4 dxij_ij
   character*10 vk_val
   type(gtp_mqmqa_var), pointer :: mqf   
!
   fact(0)='0.0'
   fact(1)='1.0'
   fact(2)='0.5'
   mqf=>phres%mqmqaf
!
!23456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
!        '23456789.123456789.123456789.123456789.123456789.123456789.123456789.
! 72-80  '123456789.------    
   write(*,100)
100 format(/'3XQ Listing the composition variables in the MQMQA model and ',&
         'their relations.'/&
         'The MQMQA model for salts has 2 sublattices but describes its ',&
         'composition using ',/&
       'a single set of sites with quadruplet fractions with a single anion.',/&
         'The quadruples are denoted x_ij here and their sum is unity and ',&
         'they provide ',/&
         'the correct mass balance. The subscript _ij denotes two cations ',&
         'and the indexing',/&
         'is normally x_(i=1..n,j=i..n) for n cations, but x_ij and x_ji ',&
         'is the same.'/&
       'For the other constitution variables the order of the index matters!',/&
         'The variables vk_ij, xi_ij and y_i/k do not fullfil ',&
         'massbalance relations.',/&
         'vk stands for the greek varkappa letter and xi for the ksi letter.'//&
         'The reference state of the MQMQA phase has one endmember parameter ',&
         'for each',/&
         'quadruplet fraction x_ij because this implementation assumes a ',&
         'single anion.',/&
         'The y_i/k variables has a redundant anion index k and is related ',&
         'to the amount ',/&
         'of cation i but it is not the correct fraction of cation i.',//&
         'The xi_ij and y_i/k are the same unless j is a Toop element in ',&
         'the ternary',/&
         'i-j-m, in such cases y_m/k is added to xi_ij.'//&
         'The quasichemical configurational entropy have additional ',&
         'constituent variables ',/&
         'which are not used for the excess Gibbs energy.',//&
         'Is this enough? ')

   write(kou,4124)mqmqa_data%nquad,mqmqa_data%ncat
4124 format(/'3XQ Listing of quads and asymmetries using varkappadefs:'/&
        'The ',i3,' quads for ',i2,' cations are arranged ',&
        'in order of the n cations:'/&
        'Quad  ',9x,'1   2  ...  n | n+1 n+2 ... 2n-1 | 2n .. | n(n+1)/2'/&
        'Cation',9x,'1   1  ...  1 | 2   2   ...  2   | 3  .. | n'/&
        'Cation',9x,'1   2  ...  n | 2   3   ...  n   | 3  .. | n')
!
! identify the actual cations in all quads as above
! create the local quadcat indices array used for the vk_ij quad dependences 
   if(.not.allocated(quadcat)) then
      allocate(quadcat(mqmqa_data%nquad))
!      write(*,*)'3XQ size of quadcat: ',size(quadcat)
   endif
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
      line1(ip:ip)='| '
     ip=ip+2
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
      line2(ip:ip)='| '
      qline(ip:ip)='| '
      ip=ip+2
   enddo
! nice output of quads and cation dependencies
   write(*,51)trim(qline)
   write(*,51)trim(line1)
   write(*,51)trim(line2)
!
! List the quades with a single cation
   call list_quads_with_single_cation(phres)
!
! This call initiates values in ivk_ij jvk_ji and denom used below
   write(*,*)'3XQ in list_mqmqa_variables call calcasymvar force reinitiate'
   call calcasymvar(phres,1)
   if(gx%bmperr.ne.0) goto 1000
!
   write(*,11)mqmqa_data%nquad
11 format(/'The ',i2,' x_ij fractions ordered: x_11, x_12, .. x_1n, ',&
        'x_22, .. x2n, ... x_nn')
   write(*,12)(mqf%xquad(ip),ip=1,mqmqa_data%nquad)
12 format(6f10.6)
   write(*,13)
13 format('x_ij or x_ji is the same value which is not true for vk_ij etc..'/)
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

! Problem that box%ivk_ij etc are not allocated ... we must call calcasymvar
!
! allocate(box%ivk_ij( ... done by the elegant [ ... ] statement
! allocate(box%jvk_ij( ... done by the elegant [ ... ] statement
!
   write(*,8778)
8778 format('3XQ All other constituent variables are functions of these ',&
          'quadruplets.',/&
          'IMPORTANT: x_ij=x_ji, the order of indices for x_ij is irrelevant.'/)

!   write(*,101)mqmqa_data%ncat*(mqmqa_data%ncat-1)/2,size(mqf%compvar)
   write(*,101)size(mqf%compvar)
101 format('3XQ The expressions for the',i3,' varkappa constituent variables.'/&
         'IMPORTANT: vk_ij is not equal to vk_ji.')
!
   vkloop2: do nv=1,size(mqf%compvar)
! _ij
      box=>mqf%compvar(nv)
!      write(*,103)nv,size(box%ivk_ij),size(box%jvk_ji),&
!           size(box%all_ijk),size(box%kvk_ijk)
!103   format('Varkappa record: ',i3,', function of quads: ',4i4)
      write(*,103)nv
103   format('Varkappa variable: ',i3,', summing quads: ')
      line1='x_'//quadcat(box%ivk_ij(1))
      ip=len_trim(line1)+1
      k=2
      do while(k.le.size(box%ivk_ij))
         line1(ip:)='+x_'//quadcat(box%ivk_ij(k))
         k=k+1
         ip=ip+5
      enddo
! To fix problems here see around line 4100 about box%ivk_ij, %jvk_ji %kvk_ijk
      write(vk_val,104)box%vk_ij
104   format(f10.6)
      write(*,105)'vk_'//char(i0+box%cat1)//char(i0+box%cat2)//&
           ' = '//vk_val//' ('//trim(line1)//')/denom'
! _ji
105   format(a)
      line2='x_'//quadcat(box%jvk_ji(1))
      ip=len_trim(line2)+1
      k=2
      do while(k.le.size(box%jvk_ji))
         line2(ip:)='+x_'//quadcat(box%jvk_ji(k))
         k=k+1
         ip=ip+5
      enddo
      write(vk_val,104)box%vk_ji
      write(*,105)'vk_'//char(i0+box%cat2)//char(i0+box%cat1)//&
           ' = '//vk_val//' ('//trim(line2)//')/denom'
! _denom
! NOTE some quad fractions may appear twice!! should be removed
      qline=trim(line1)//'+'//trim(line2)//' +x_'//quadcat(box%kvk_ijk(1))
      ip=len_trim(qline)+1
      k=2
      do while(k.le.size(box%kvk_ijk))
         qline(ip:)='+x_'//quadcat(box%kvk_ijk(k))
         k=k+1
         ip=ip+5
      enddo
      write(*,106)'denom: = '//trim(qline)
106   format(19x,a)
   enddo vkloop2
!
!=====================================================================
! Now list expression and value of for y_ik
! We have x_11, x_12, .... x_1n, x_22, x_23, ... x_2n, x_33, ... x_nn
!  y_i/k = 0.5 * x_ii + 0.5 * x_ij  
! y_1/k = x_11 + 0.5x_12 + 0.5*x_13 + ... + 0.5*x_15 + ... + 0*x_22 + ...
! y_5/k = 
!
! x_11, x_12, x_13 / x_12, x_22, x_23 / x_13, x_23, x_33;  NOTE x_ij=x_ji
!
   write(*,109)mqmqa_data%ncat
109 format(/'3XQ The',i3,' y_i/k constituent variables as functions of quads.',&
         /'Remember x_ij and x_ji represent the same quadruplet fraction!')
   xijheader='y_i/k:  value       '
   ip=18
   header1: do nv=1,mqmqa_data%ncat
      header2: do nv2=nv,mqmqa_data%ncat
         xijheader(ip:ip+4)='x_'//char(i0+nv)//char(i0+nv2)
         ip=ip+5
      enddo header2
   enddo header1
   lenhead=ip
!
   write(*,110)xijheader(1:lenhead)
   ally_ik: do ix=1,mqmqa_data%ncat
      line1='y_'//char(i0+ix)//'/k:  '
      write(line1(7:16),'(F10.6)')mqf%y_ik(ix)
      ip=17
      coeffs1: do nv=1,mqmqa_data%nquad
         if(mqmqa_data%dy_ik(ix,nv).ne.0.0D0) then
            write(line1(ip:ip+4),808)mqmqa_data%dy_ik(ix,nv)
808         format(F4.1)
         else
            line1(ip:ip+4)='  - '
         endif
         ip=ip+5
      enddo coeffs1
      write(*,110)line1(1:ip)
110   format(a)
   enddo ally_ik
   write(*,810)
810 format('VERY IMPORTANT: The y_i/k variable is not the correct fraction ',&
         'of the element i'/&
         'in the phase because it ignores the stoichiometry'/&
         'of the quadruplet species.')
!
!=====================================================================
! Now list the xi_ij and their asymmetries
!   write(*,*)'3XQ ******************************************************'
!   write(*,*)'3XQ line 6309 size dxij_ij: ',size(mqf%compvar(1)%dxi_ij)
!   write(*,*)'3XQ ******************************************************'
!
   write(*,200)2*size(mqf%compvar)
200 format(/'3XQ line 6301 The ',i3,' xi_ij/k constituent variables.'/&
         'For symmetric systems xi_ij/k = y_i/k. ',&
         'But if "j" is a Toop element in',/&
         'the ternary i-j-\nu then y_\nu/k is added to xi_ij/k.'/&
         'REMEMBER xi_ij and xi_ji are NOT THE SAME !!!',/&
         'Below the xi_ij, xi_ji are expressed as sums of the symmetric ',&
         'x_ij fractions.')
!         'It is important to remember that xi_ij IS NOT EQUAL to xi_ji! ')
!
! Set xi_ij coefficient for a symmetric system xi_ij = y_i/k; xi_ji = y_j/k
!
   do ix=1,size(mqf%compvar)
      do jx=1,mqmqa_data%nquad
        mqf%compvar(ix)%dxi_ij(jx)=mqmqa_data%dy_ik(mqf%compvar(ix)%cat1,jx)    
        mqf%compvar(ix)%dxi_ji(jx)=mqmqa_data%dy_ik(mqf%compvar(ix)%cat2,jx)    
      enddo
   enddo
!
! Below is just listing
   xijheader(1:8)='xi_ij:  '
   write(*,110)xijheader(1:lenhead)
   ix=1
   jx=2
   xiloop: do nv=1,size(mqf%compvar)
      line1='xi_'//char(i0+ix)//char(i0+jx)//':'
      write(line1(7:16),'(F10.6)')mqf%compvar(nv)%xi_ij
      ip=17
! xij_ij should be the same as y_ik
      do nv2=1,mqmqa_data%nquad
         if(mqf%compvar(nv)%dxi_ij(nv2).eq.0.0D0) then
            line1(ip:ip+4)='   - '
         else
            write(line1(ip:ip+4),777)mqf%compvar(nv)%dxi_ij(nv2)
777         format(f5.1)
         endif
         ip=ip+5
      enddo
      write(*,110)line1(1:ip)
      line1='xi_'//char(i0+jx)//char(i0+ix)//':'
      write(line1(7:16),'(F10.6)')mqf%compvar(nv)%xi_ji
      ip=17
! xji_ij should be the same as y_jk
      do nv2=1,mqmqa_data%nquad
         if(mqf%compvar(nv)%dxi_ji(nv2).eq.0.0d0) then
            line1(ip:ip+4)='   - '
         else
            write(line1(ip:ip+4),777)mqf%compvar(nv)%dxi_ji(nv2)
         endif
         ip=ip+5
      enddo
      write(*,110)line1(1:ip)
      jx=jx+1
      if(jx.gt.mqmqa_data%ncat) then
         ix=ix+1
         jx=ix+1
      endif
   enddo xiloop
!
   write(*,333)
333 format(/'3XQ The symmetries of the ternary systems')
   call list_tersys
   goto 1000
!
! code below in copied in subroutine above
   ts: if(allocated(tersys)) then
! this popular listing is repearad 2 in gtp3XQ.F90 and 3 times in pmon6.F90
      write(*,3101)size(tersys)
3101  format(/'3XQ Listing of the',i3,' ternary systems and their asymmetries',&
  /'  i  seq   cat1 cat2 cat3',6x,'emquads',6x,'T/0 T/0 T/0    asymmetry code')
!           /'  i tern   cat1 cat2 cat3       T/0 T/0 T/0    asymmetry code')
      do iz=1,size(tersys)
         write(*,3201)iz,tersys(iz)%seq,(tersys(iz)%el(j4),j4=1,3),&
              tersys(iz)%emquad,tersys(iz)%isasym,tersys(iz)%asymm
3201     format(i3,i5,2x,3(1x,i4),2x,3(1x,i3),3x,3i4,5x,a)
      enddo
      write(*,3301)
3301  format('Number in cat1/2/3 columns is cation index,'/&
           'Number 1, 2 or 3 in T/0 columns refers to ',&
           'the cat1/2/3 COLUMN, NOT CATION INDEX!'/&
           'Asymmetry code is KKK for symmetric, KKT if cat3 is Toop etc.'//&
           'Change the asymmetry with the command AMEND PHASE ... ASYM')
   else
      write(kou,*)'No ternary asymmetry data allocated'
   endif ts
!
1000 continue
!   write(*,1100)
1100 format(/'3XQ leaving list_mqmqa_variables'/)
   return
 end subroutine list_mqmqa_variables

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine list_tersys
!\begin{verbatim}
 subroutine list_tersys
! list expressions for vk, xi and y and current values for all compvar
   implicit none
!   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
!   type(gtp_mqmqa_var), pointer :: mqf
!   type(gtp_allinone), pointer :: box
   integer iz,j4
!
   ts: if(allocated(tersys)) then
! this popular listing is repeated twice in gtp3XQ.F90, 3 times in pmon6.F90
      write(*,3101)size(tersys)
3101  format(/'3XQ Listing of the',i3,' ternary systems and their asymmetries',&
  /'  i  seq   cat1 cat2 cat3',6x,'emquads',6x,'T/0 T/0 T/0    asymmetry code')
!           /'  i tern   cat1 cat2 cat3       T/0 T/0 T/0    asymmetry code')
      do iz=1,size(tersys)
         write(*,3201)iz,tersys(iz)%seq,(tersys(iz)%el(j4),j4=1,3),&
              tersys(iz)%emquad,tersys(iz)%isasym,tersys(iz)%asymm
3201     format(i3,i5,2x,3(1x,i4),2x,3(1x,i3),3x,3i4,5x,a)
      enddo
      write(*,3301)
3301  format('Number in cat1/2/3 columns is cation.'/&
           ' A 1 in T/0 columns mean Toop')
   else
      write(kou,*)'No ternary asymmetry data allocated'
   endif ts
!
   return
 end subroutine list_tersys

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine nathalie_asym
!\begin{verbatim}
 subroutine nathalie_asym(phres,xverbose)
! called from new_ternary_asym around line 4108
! to generate varkappa_ij and varkappa_ji from exiting asymmetries
   implicit none
   type(gtp_phase_varres), pointer :: phres
   logical xverbose
!\end{verbatim}
! Nathalies algorith,=m
! Initialize all the vk_nm=x_nn
! For each ternary ijk
!    If i Toop, add to numerator
!        of vk_ji : x_kk
!        of vk_ki : x_jj
!    If j Toop, add to numerator
!        of vk_jj : x_kk
!        of vk_kj : x_ii
!    If k Toop, add to numerator
!        of vk_jk : x_ii
!        of vk_ik : x_jj
! After the loop on each ternary, 
! For each binary ab, add all the x_qq in the numerator of vk_ab and vk_ba 
!               to get the denominator
!
! For all the numerators and denominators, 
! add the mixed terms x_pq for each p and q present 
!   
   integer t1,toopix,toopem,ia,jb,vk1,vk2,icat,jcat,kcat
   type(gtp_mqmqa_var), pointer :: mqf   
   type(gtp_allinone), pointer :: box1,box2
   logical verbose
!
! run time error if I try to set xverbose, ignore it and just set it true here
   verbose=.true.
   if(verbose) write(*,10)
10 format('3XQ using Nathalies algorithm to generate asymmetric varkappa!')
   mqf=>phres%mqmqaf
   terloop1: do t1=1,size(tersys)
      if(tersys(t1)%asymm.eq.'KKK') cycle terloop1
!
      toopix=index(tersys(t1)%asymm,'T')
! toopem is cation index for Toop 
      toopem=tersys(t1)%el(toopix)
      if(tersys(t1)%noasym.ne.0) then
         write(*,*)'3XQ ternary asymmetry already set for ',t1
         cycle terloop1
      endif
      tersys(t1)%noasym=1
!      toopem=tersys(t1)%emquad(toopix),
! %el is cations  1-2-3, 1-2-4 etc   %emquad is the actual cation
!      icat=tersys(t1)%emquad(1);
!      jcat=tersys(t1)%emquad(2);
!      kcat=tersys(t1)%emquad(3);
      icat=tersys(t1)%el(1)
      jcat=tersys(t1)%el(2);
      kcat=tersys(t1)%el(3);
! single anion
      ia=1
      if(verbose) write(*,100)t1,toopix,&
           tersys(t1)%el,icat,jcat,kcat,tersys(t1)%binsys
100   format('3XQ ternary:',i2,' Toop: ',i2,&
           ' cations' ,3i3,3x,3i3,' binaries: ',3i3)
! ternary cations are in increasing order associated with binaries (boxes)
! 1-2-3         1-2-4,       1-3-4,       2-3-4           tersys
! 1-2 1-3 2-3   1-2 1-4 2-4  1-3 1-4 3-4  2-3 2-4 3-4     compvar()%icat %jcat
! 1   2   4     1   3   5    2   3   6    4   5   6       compvar(
!
! test 4 cation case is TKK, KKT, KKT, KTK
! this find the binaries involved  .... elaborated below
      if(toopix.eq.1) then
         vk1=tersys(t1)%binsys(1); vk2=tersys(t1)%binsys(2)   ! 1-2 and 1-3
      elseif(toopix.eq.2) then
         vk1=tersys(t1)%binsys(1); vk2=tersys(t1)%binsys(3)   ! 1-2 and 2-3
      else
         vk1=tersys(t1)%binsys(2); vk2=tersys(t1)%binsys(3)   ! 1-3 and 2-3
      endif
!
!   stop 'we are here 1'
!
      write(*,*)'3XQ binaries with asymmetries: ',vk1,vk2
      box1=>mqf%compvar(vk1); box2=>mqf%compvar(vk2)
! what information in the boxes?
      write(*,70)vk1,box1%cat1,box1%cat2,box1%elcat1,box1%elcat2,&
           box1%quadicat1,box1%quadicat2
      write(*,70)vk2,box2%cat1,box2%cat2,box2%elcat1,box2%elcat2,&
           box2%quadicat1,box2%quadicat2
70    format('3XQ box ',i2,': ',2i3,5x,2i3,5x,2i3)
! evidently %cat1 and %cat2 are the cation indices, %elcat the quad indices
!
! these are the inital varkappa asymmetries
      write(*,110)vk1,box1%ivk_ij
      write(*,110)vk1,box1%jvk_ji
      write(*,110)vk2,box2%ivk_ij
      write(*,110)vk2,box2%jvk_ji
110   format('3XQ box ',i3,' already sums quads:',10i3)
!
! now try to do something sensible inside the if-statement
      write(*,*)'3XQ toopix: ',toopix,toopem
      if(toopix.eq.1) then
         vk1=tersys(t1)%binsys(1); vk2=tersys(t1)%binsys(2)   ! 1-2 and 1-3
         box1=>mqf%compvar(vk1); box2=>mqf%compvar(vk2)
! For each ternary ijk  %cat2 
!    If i Toop, add to numerator
!        of vk_ji : x_kk
!        of vk_ki : x_jj
! REMEMBER order of icat,jcat in ijklx is irrelevant, x_ij is symmetric
! all quads of numerators included in the denominator
! ijklx(icat,jcat,ia,ia) gives index of the quadruplet with icat+jcat
         box1%ivk_ij=[box1%ivk_ij, ijklx(kcat,kcat,ia,ia)]
         box2%ivk_ij=[box2%ivk_ij, ijklx(jcat,jcat,ia,ia)]
!         box1%ivk_ij=[box1%ivk_ij, ijklx(jcat,jcat,ia,ia)]
!         box2%ivk_ij=[box2%ivk_ij, ijklx(kcat,kcat,ia,ia)]
!         box1%ivk_ij=[box1%ivk_ij, ijklx(jcat,jcat,ia,ia), &
!              ijklx(icat,jcat,ia,ia)]
!         box2%ivk_ij=[box2%ivk_ij, ijklx(kcat,kcat,ia,ia), &
!              ijklx(icat,kcat,ia,ia)]
      elseif(toopix.eq.2) then ! --------------------------------------
!
!    If j Toop, add to numerator  NOTE NOT SAME icat, kcat as previous!!!
!        of vk_ij : x_kk
!        of vk_kj : x_ii
         vk1=tersys(t1)%binsys(1); vk2=tersys(t1)%binsys(3)   ! 1-2 and 2-3
         box1=>mqf%compvar(vk1); box2=>mqf%compvar(vk2)
         box1%jvk_ji=[box1%jvk_ji, ijklx(icat,icat,ia,ia)]
         box2%ivk_ij=[box2%ivk_ij, ijklx(kcat,kcat,ia,ia)]
!         box1%jvk_ji=[box1%jvk_ji, ijklx(icat,icat,ia,ia)
!              ijklx(icat,jcat,ia,ia)]
!         box2%ivk_ij=[box2%ivk_ij, ijklx(kcat,kcat,ia,ia), &
!              ijklx(jcat,kcat,ia,ia)]
      else                    ! --------------------------------------
!    If k Toop, add to numerator
!        of vk_jk : x_ii
!        of vk_ik : x_jj
         write(*,190)icat,jcat,kcat,box2%ivk_ij
190      format('3XQ k is Toop: ',3i3,5x,10i3)
         vk1=tersys(t1)%binsys(2); vk2=tersys(t1)%binsys(3)   ! 1-3 and 2-3
         box1=>mqf%compvar(vk1); box2=>mqf%compvar(vk2)
         box1%jvk_ji=[box1%jvk_ji, ijklx(icat,icat,ia,ia)]
         box2%jvk_ji=[box2%jvk_ji, ijklx(jcat,jcat,ia,ia)]
!         box1%jvk_ji=[box1%jvk_ji, ijklx(icat,icat,ia,ia), &
!              ijklx(icat,kcat,ia,ia)]
!         box2%ivk_ij=[box2%ivk_ij, ijklx(kcat,kcat,ia,ia), &
!              ijklx(jcat,kcat,ia,ia)]
      endif
      write(*,210)vk1,mqf%compvar(vk1)%ivk_ij
      write(*,210)vk1,mqf%compvar(vk1)%jvk_ji
      write(*,210)vk2,mqf%compvar(vk2)%ivk_ij
      write(*,210)vk2,mqf%compvar(vk2)%jvk_ji
210   format('3XQ box: ',i3,' now sums quads: ',10i3)
!
   enddo terloop1
!

!   box1=>mqf%compvar(vk1); box2=>mqf%compvar(vk2)
!   write(*,200)toopem,vk1,box1%cat1,box1%cat2,vk2,box2%cat1,box2%cat2
200 format('3XQ Toop is ',i2,' binaries with Toop ',i3,': ',2i3,&
                                            ' and ',i3,': ',2i3)
!
! we have toopem in vk1 we have toop either for vi_ij
!   box1%vk_ij=[box1%vk_ij , x_ii]

!   stop 'we are here 2'
!
! For each ternary ijk
!    If i Toop, add to numerator
!        of vk_ji : x_kk
!        of vk_ki : x_jj
!    If j Toop, add to numerator
!        of vk_jj : x_kk            should be        of vk_ij : x_kk
!        of vk_kj : x_ii
!    If k Toop, add to numerator
!        of vk_jk : x_ii
!        of vk_ik : x_jj
!
!
! After the loop of all ternaries, 
! For each binary ab, add all the x_qq in the numerator of vk_ab and vk_ba 
!               to get the denominator
!
! For all the numerators and denominators, 
! add the mixed terms x_pq for each p and q present 
!
   terloop2: do t1=1,size(tersys)
! add mixed x_ia, x_jb
      write(*,*)'3XQ should now loop all tersys for crossterms ... '
   enddo terloop2
! maybe some final data ...

1000 continue
   write(*,*)'3XQ all done!'
   return
 end subroutine nathalie_asym

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine list_compvar
!\begin{verbatim}
 subroutine list_compvar(phres)
! list expressions for vk, xi and y and current values for all compvar
   implicit none
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   type(gtp_mqmqa_var), pointer :: mqf
   type(gtp_allinone), pointer :: box   
   character*300 line1,line2,qline
   character*10 vk_val
   character*2, dimension(:), allocatable :: quadcat
   integer nv,i0,ip,i,j,k,q1,q2,q3
!
   write(*,*)'3XQ in list_compvar'
   mqf=>phres%mqmqaf
!
   allocate(quadcat(mqmqa_data%nquad))
   i0=ichar('0')
   k=1
! initiate quadcat to be 2 indices appropriate to quadindices
   iniquadcat: do i=1,mqmqa_data%ncat
      do j=i,mqmqa_data%ncat
         quadcat(k)=char(i0+i)//char(i0+j)
         k=k+1
      enddo
   enddo iniquadcat
   write(*,10)quadcat
10 format('3XQ x_ij   : ',20(a,' '))
   write(*,11)(i,i=1,size(quadcat))
11 format('3XQ indices: ',20(i2,1x))
!
!
   vkloop2: do nv=1,size(mqf%compvar)
! _ij
      box=>mqf%compvar(nv)
      write(*,103)nv,size(box%ivk_ij),size(box%jvk_ji),size(box%kvk_ijk)
103   format('Varkappa variable: ',i3,', sizes: ',3i3,', summing quads: ')
      line1='x_'//quadcat(box%ivk_ij(1))
      ip=len_trim(line1)+1
      k=2
      do while(k.le.size(box%ivk_ij))
         line1(ip:)='+x_'//quadcat(box%ivk_ij(k))
         k=k+1
         ip=ip+5
      enddo
! To fix problems here see around line 4100 about box%ivk_ij, %jvk_ji %kvk_ijk
      write(vk_val,104)box%vk_ij
104   format(f10.6)
      write(*,105)'vk_'//char(i0+box%cat1)//char(i0+box%cat2)//&
           ' = '//vk_val//' ('//trim(line1)//')/denom'
! _ji
105   format(a)
      line2='x_'//quadcat(box%jvk_ji(1))
      ip=len_trim(line2)+1
      k=2
      do while(k.le.size(box%jvk_ji))
         line2(ip:)='+x_'//quadcat(box%jvk_ji(k))
         k=k+1
         ip=ip+5
      enddo
      write(vk_val,104)box%vk_ji
      write(*,105)'vk_'//char(i0+box%cat2)//char(i0+box%cat1)//&
           ' = '//vk_val//' ('//trim(line2)//')/denom'
! _denom
! NOTE some quad fractions may appear twice!! should be removed
      qline=trim(line1)//'+'//trim(line2)//' +x_'//quadcat(box%kvk_ijk(1))
      ip=len_trim(qline)+1
      k=2
      do while(k.le.size(box%kvk_ijk))
         qline(ip:)='+x_'//quadcat(box%kvk_ijk(k))
         k=k+1
         ip=ip+5
      enddo
      write(*,106)'denom: = '//trim(qline)
106   format(11x,a)
   enddo vkloop2
!------------------------------ indices
   write(*,200)
200 format('3XQ same information using quad indices')
   write(*,10)quadcat
!10 format('3XQ line 7339: ',20(a,' '))
   write(*,11)(i,i=1,size(quadcat))
!11 format('3XQ indices  : ',20(i2,1x))
!
! repeat using quad fraction indices
   vkloop3: do nv=1,size(mqf%compvar)
      box=>mqf%compvar(nv)
! ij
      line1=' '
      q2=1
      do q1=1,size(box%ivk_ij)
         write(line1(q2:q2+3),205)box%ivk_ij(q1)
205      format(' +',i2)
         q2=q2+4
      enddo
      write(*,105)'vk_'//char(i0+box%cat1)//char(i0+box%cat2)//&
           ' = ('//trim(line1)//')/denom'
! ji
      line2=' '
      q2=1
      do q1=1,size(box%jvk_ji)
         write(line2(q2:q2+3),205)box%jvk_ji(q1)
         q2=q2+4
      enddo
      write(*,105)'vk_'//char(i0+box%cat2)//char(i0+box%cat1)//&
           ' = ('//trim(line2)//')/denom'
! denom
      qline=trim(line1)//trim(line2)
      q2=len_trim(qline)+1
      do q1=1,size(box%kvk_ijk)
         write(qline(q2:q2+3),205)box%kvk_ijk(q1)
         q2=q2+4
      enddo
      write(*,105)'      demom = '//trim(qline)
   enddo vkloop3
!
   return
 end subroutine list_compvar
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine list_asymmetries
!\begin{verbatim}
 subroutine list_asymmetries(phres)
! list expressions for vk, xi and y and current values for all compvar
   implicit none
   type(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   write(*,*)'3XQ list_asymmetries_dummy'
   return
 end subroutine list_asymmetries

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

! asymmetry code
! j is the Toop element in i-j-\nu
! i is the Toop element in i-j-\gamma
!
!                 \sum_a=(i,\nu) \sum_b=(i,\nu) x_ab/kk                ivk_ij
! vk_ij/kk = ------------------------------------------------------- = -------
!            \sum_a=(i,j,\nu,\gamma) \sum_b=(i,j,\nu,\gamma) x_ab/kk   denom_ij
!
!                 \sum_a=(j,\gamma) \sum_b=(j,\gamma) x_ab/kk          jvk_ji
! vk_ji/kk = ------------------------------------------------------- = ------
!            \sum_a=(i,j,\nu,\gamma) \sum_b=(i,j,\nu,\gamma) x_ab/kk   denom_ij
!
! NOTE x_ij = x_ji and occures only once in sums !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ivk_ij = x_i,i + x_i,\nu + x_\nu,\nu
! jvk_ji = x_j,j + x_j,\gamma + x_\gamma,\gamma
! denom  = x_i,j + x_i,\nu+x_j,\gamma+x_\nu,\nu+x_\nu,\gamma+x_gamma,gamma
!
! initiate: ivk_ij=[x_ii]; jvk_ji=[x_jj]; denom=[x_ij]
!
! extradenom=[ ]
! binary loop vk: do i-j
!   ternary loop: do g=1,n   ------------------   g can be \nu, \gamma or both
!     if(g=i or g=j) cycle ternary loop
!     if(i is Toop in i-j-g) then  ...............g is \gamma
!       jvk_ij=[ jvk_ij , x_gg, x_jg ]
! denom will at the end have jvk_ji and ivk_ij added.  Add only x_ig
!       denom_ij = [ denom_ij, x_ig]
!       if(j is Toop in i-j-g) then ..............g is both \nu and \gamma
!         ivj_ji=[ ivk_ij, x_gg, x_ig, x_jg ]
!       endif
! there can have been previous \gamma or \nu, add extra x_\gamma,\nu
!       do h=1,size(extradenom)
!         denom_ij = [ denom_ij, x_gh ]
!       enddo
!       extradenom = [extradenom, g ]
!-----------
!     elseif(j is Toop in i-j-g) then ...........g is \nu
!       ivj_ji=      [ ivk_ij, x_gg, x_jg ]
!       denom_ij = [ denom_ij, x_gg, x_jg, x_ig ]
!     endif
!   enddo ternary loop
! enddo binary loop
!
!--------------------------- correct asymmetry:

! ONLY ONE TOOP ELEMENT PER TERNARY
!---------------------------------------------------------------
! This is previous OC version when I had correct asymmetries, asymm=TKK
! The asymmetric element here is 1
! Varkappa index:   1, summing quads: 
!   nomin: vk_12 =(x_11)/denom
!   nomin: vk_21 =(x_22+x_23+x_33)/denom
!         denom: = x_11+x_22+x_23+x_33 +x_12+x_13
!Varkappa index:   2, summing quads: 
!   nomin: vk_13 =(x_11)/denom
!   nomin: vk_31 =(x_33+x_23+x_22)/denom
!         denom: = x_11+x_33+x_23+x_22 +x_13+x_12
!Varkappa index:   3, summing quads: <<<<<<<<<<<<<<<<< no change
!   nomin: vk_23 =(x_22)/denom
!   nomin: vk_32 =(x_33)/denom
!         denom: = x_22+x_33 +x_23
!
!---------------------------------------------------------------
!
! This is previous OC version when I had correct asymmetries, asymm=KTK
! The asymmetric element here is 2
!   Varkappa index:   1, summing quads: 
!   nomin: vk_12 =(x_11+x_13+x_33)/denom
!   nomin: vk_21 =(x_22)/denom
!         denom: = x_11+x_13+x_33+x_22 +x_12+x_23
!Varkappa index:   2, summing quads:  <<<<<<<<<<<<<<<<< no change
!   nomin: vk_13 =(x_11)/denom
!   nomin: vk_31 =(x_33)/denom
!         denom: = x_11+x_33 +x_13
!Varkappa index:   3, summing quads: 
!   nomin: vk_23 =(x_22)/denom
!   nomin: vk_32 =(x_33+x_13+x_11)/denom
!         denom: = x_22+x_33+x_13+x_11 +x_23+x_12
!
!-------------------------------------------------------------   
!
! This is previous OC version when I had correct asymmetries, asymm=KKT
!Varkappa index:   1, summing quads:  <<<<<<<<<<<<<<<<< no change
!   nomin: vk_12 =(x_11)/denom
!   nomin: vk_21 =(x_22)/denom
!         denom: = x_11+x_22 +x_12
!Varkappa index:   2, summing quads: 
!   nomin: vk_13 =(x_11+x_12+x_22)/denom
!   nomin: vk_31 =(x_33)/denom
!         denom: = x_11+x_12+x_22+x_33 +x_13+x_23
!Varkappa index:   3, summing quads: 
!   nomin: vk_23 =(x_22+x_12+x_11)/denom
!   nomin: vk_32 =(x_33)/denom
!         denom: = x_22+x_12+x_11+x_33 +x_23+x_13
!   
!-----------------------------------------------------------

