! Data structures and routines for step/map/plot (using gnuplot)
!
MODULE ocsmp
!
! Copyright 2012-2021, Bo Sundman, France
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!------------------------------
!
  use liboceqplus
!
  implicit none
  character*8, parameter :: smpversion='SMP-2.30'
!
! note the type map_fixph declared in matsmin.F90 (in liboceq)
!
! Thought for smp3.  A new method to calculate diagrams with tie-line in plane
! No fix phase but instead a fixed compostion in the middle of the tie-line
! That means no problems changing fix phase!
!
! MAP_NODE records are created whenever the set of stable phases changes.
! It can have links to two or more MAP_LINE records with calculated equilibra.
! Initially at least one of these map_line records are empty with just 
! information of the axis to vary and direction.
! There is a FIRST map_node record and subsequent are linked by a double linked
! list, next and previous.  All map_node records have a pointer to the first.
! There are a special use of map_node records when the node has as many stable
! phases as the line and it has only two lines connected:
! - when starting following a line from a start point
! - when the array of calculated equilibria must be saved on a file a map_node
! record is created for each unifinished line.  The map_node and map_line
! records must be saved on the random file also and initiated for the next.
! All map_node records must be kept for the next block.  It is possible that
! a line is being calculated leading to a node with no exits.  If that node 
! is removed it will be created again and all already calculated lines exiting
! will be calculated again ...
!
! MAP_LINE records are created for each line followed during the step/map.  It 
! contains links to stored gtp_equilibrium_data records and some
! additional info.
! The stored gtp_equilibrium_data records belong to an array that can be saved
! on a random file and the links are the index to reords in this array.  The
! gtp_equilibrium_data records are linked internallw with these indices also.
! The map_line record has links to two map_node records representing the
! two ends of the line.  A map_line that terminates at the end of an axis will
! have a zero link for that end of the line.
!
!\begin{verbatim}
! These are bits in the map_line status word
! if EXCLUDEDLINE set the whole line is inactived
! TWOSTOICH set if there are tie-lines inplane and line ends with stoich phases
! with same composition (see for example U-O)
  integer, parameter :: EXCLUDEDLINE=0, TWOSTOICH=1
! Bit in the MAP_NODE record
  integer, parameter :: MAPINVARIANT=0,STEPINVARIANT=1
!\end{verbatim}
!
!\begin{verbatim}
  TYPE map_line
! This is record contains a list of calculated equilibria along a line
! These are pointers to map_node records at start and end of line.
     type(map_node), pointer :: start,end
! For threading this record must be separate from other threads
! This record is created by calling calceq7 the first equilibrium of line
     type(meq_setup) :: meqrec
! the active ceq record must be stored together with all others in order to be
! copied from the ceq record saved with the node record when line starts
     type(gtp_equilibrium_data), pointer :: lineceq
! this is the number of calculated equilibria for this line and the index
! of the first and last one stored.  
! The stored equilibria has an internal next link.
! lineid is a sequential index of the lines. done is negative if done
! nfixphases are the number of fixed phases replacing axis conditions
! status kan be used to delete a line
     integer number_of_equilibria,first,last,lineid,done,nfixphases,status
! This is used during mapping to identify lines that have the same fixed phases
! if we have 3 or more axis there can be 2 or more fix phases along the line??
! With 2 axis there can only be one fix phase!
     type(gtp_phasetuple), dimension(:), allocatable :: linefixph
! also save index to phr!! (do not trust ...)
     integer, dimension(:), allocatable :: linefix_phr
! Save the phase tuplet representing the phase fix at start node here
! If it wants to be stable at first step along a line change axis direction
!     type(gtp_phasetuple) :: nodfixph <<<< not used
! This is the phase index in the phr array (phr has both phase and compset)
     integer nodfixph
! We must also save the number and set of stable phases and theit amounts
! as we will have different stable phases for different lines
     integer nstabph
     type(gtp_phasetuple), dimension(:), allocatable :: stableph
     double precision,     dimension(:), allocatable :: stablepham
! added also index to phr as that seems useful
     integer,              dimension(:), allocatable :: stable_phr
! axandir is set when linenode is created to the axis and direction for first
! step from the node.  It can be changed to another axis and direction
! during map and indicate the current axis with active condition
! axchange remember the equilibrum number for an axis change
     integer axandir,axchange
! more is 1 while following the line, 0 for last equilbrium, -1 when finished
! termerr is zero unless line terminated with error, -1 means exit not used
! problem is nonzero if map_problems has been called
! lasterr is the last error occured calculating this line
     integer more,termerr,problems,lasterr
! firstinc is a value to add the the axis variable for the first equilibrium
! to avoid finding the node point again.  Evenvalue is the next value
! to calculate during a step calculation.  Both set when creating the node.
! At start nodes they are zero
     double precision firstinc,evenvalue
! During map the last axis values for ALL axis are stored here
     double precision, dimension(:), allocatable :: axvals
! If tie-lines in the plane we must also check the axis values for
! the other line as we may have to change the fix phase
     double precision, dimension(:), allocatable :: axvals2
! save previous values of axvals to handle axis changes ...
     double precision, dimension(:), allocatable :: axvalx
! save previous changes in axis values, for tie-line in plane
! dxval(phase,axis)
     double precision, dimension(:,:), allocatable :: dxval
! factor to control length of step in axis with axtive condition
     double precision :: axfact
! data particular to a step calculation, for example scheil
!     character*24, allocatable, dimension(:) :: stepresultid
  end TYPE map_line
!\end{verbatim}
!
!-------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE map_node
! this record organizes the step/map results.  Whenever there is a 
! change of the set of stable phases a node record is created and it
! can have links to several map_line records.  The map node record has a
! link to a gtp_equilibrium_data record (ceq) for the equilibrium at the node.
! This is copied to the map_line record when this is activated.  
! In the map_line record an axis and direction to start is stored.
! NOTE all gtp_equilibrium_data (ceq) records are pointers to the global
! array as new composition sets may be created along any line.
! The node record is identified by the set of stable phases and the
! chemical potentials of the components.  One must be able to identify the
! node as one may fins the same node following different lines.
! locally stored linerecords for lines exiting the node
     type(map_line), dimension(:), allocatable :: linehead
! links to other nodes
! plotlink is used to overlay two or more map or step commands
     type(map_node), pointer :: first,next,previous,plotlink
! saved copy of the meqrec record used to calculate the node 
     type(meq_setup) :: meqrec
! link to saved copy of the equilibrium record
     type(gtp_equilibrium_data), pointer :: nodeceq
! link to array of saved equilibrium record.   (only maptop?)
     type(map_ceqresults), pointer :: saveceq
! copy of nodeceq in saveceq (composition sets not updated but needed for plot)
     integer savednodeceq
! type_of_node not used?? Proposal
! =1 step normal; =2 step_separate; =3 Scheil; =4 Tzero; =5 Paraeq; =6 NPLE
! =10 map_tieline_inplane; =11 map_isotherm; ! =20 map_isopleth
! lines are number of line records
! noofstph is number of stable phases (copied from meqrec)
! tieline_inplane is 1 if so, 0 if step, -1 if no tie-lines (only maptop)
! number_ofaxis is the number of axis, 1=step;  (only maptop)
! artxe (extra) to indicate that the node has two stoichiom phases
! status for some bits maybe
! globalcheckinterval set when created from integer mapglobalcheck
     integer type_of_node,lines,noofstph,tieline_inplane,number_ofaxis,artxe
     integer status,globalcheckinterval
! seqx is unique identifier for a map node
! seqy unique identifier for maplines, incremented for each line (only maptop)
     integer seqx,seqy
! nodefix is the phase held fix when calculating node
     type(gtp_phasetuple) :: nodefix
! Value of T and P, copied from meqrec
     double precision, dimension(2) :: tpval
! chemical potentials, copied from meqrec
     double precision, dimension(:), allocatable :: chempots
! stable phase+compset, copied from meqrec (not used?)
     type(gtp_phasetuple), dimension(:), allocatable :: stable_phases
  end TYPE map_node
!\end{verbatim}
!
!-------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE map_axis
! description of the axis variables used for step/map
! The axis condition in bits and pieces
     integer nterm,istv,iref,iunit
     integer, dimension(:,:), allocatable :: indices
     type(gtp_state_variable), dimension(:), allocatable :: axcond
     double precision, dimension(:), allocatable :: coeffs
! the min, max and increment along the axis
     double precision axmin,axmax,axinc
! more must be initiated to 0, if nonzero replaced by a fixed phase
! seqz is the sequential index of the condition in the list (this is not
! changed if conditions are added (at the end) or deleted (active=1)
! we cannot use a pointer as that depend on the current equilibrium.
     integer more,seqz
! This is the last succesfully calculated axis value
     double precision lastaxval
  end TYPE map_axis
!\end{verbatim}
! decrlared as an array with each axis as one element of the array
!
!-------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE map_ceqresults
! stores calculated equilibrium records
     integer size,free,index
     TYPE(gtp_equilibrium_data), dimension(:), allocatable :: savedceq
  end TYPE map_ceqresults
!\end{verbatim}
!
!--------------------------------------------------------------
!
!\begin{verbatim}
  TYPE graphics_textlabel
! To put labels on a graph we must store these in a list
     TYPE(graphics_textlabel), pointer :: nexttextlabel
     double precision xpos,ypos,textfontscale
     integer angle
     character*80 textline
  end type graphics_textlabel
!
!\end{verbatim}
!
!------------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE plot_line
! here various information about a line to be plotted should be stored 
! for the moment it is under construction in parallel with old structures
! the "plot_line" records form a linked list starting at plotline1
     type(plot_line), pointer :: nextline
! linetype 1=normal; 2=binary invariant; 3=ternary monovariant; 4=tieline
! linetype -1=end of plotlines
     integer type
     integer active
  end type plot_line
!\end{verbatim}
!
!------------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE starteqlista
! links to equilibria that are used as start points for step or map
     type(gtp_equilibrium_data), pointer :: p1
  end type starteqlista
!\end{verbatim}
  type (starteqlista), dimension(20) :: starteqs
  integer noofstarteq
!
!------------------------------------------------------------------------
!
!\begin{verbatim}
  TYPE graphics_options
! setting options for the plotting, this replaces most arguments in the call
! to ocplot2(ndx,pltax,filename,maptop,axarr,form)
! ndx is mumber of plot axis, pltax is text with plotaxis variables
! filename is intermediary file (maybe not needed)
! maptop is map_node record with all results
! form is type of output (screen/postscript/pdf(acrobat)/gif)
!------------------------------------------------------------------
! status contain bits, BITS defined below (GRKEEP etc)
! rangedefaults(i) nonzero if min/max for axis i set by user
! axistype(i) is 1 if axis i is logscale
! plotmin/max are user definied min/max
! defltmin/max are default min/max (generated by the ploting software)
     integer :: status=0,rangedefaults(3)=0,axistype(2)=0,setgrid=0
     double precision, dimension(3) :: plotmin,plotmax
     double precision, dimension(3) :: dfltmin,dfltmax
! number of axis used for calculation (STEP=1 MAP=2 or more)
     integer noofcalcax
! scalefact is by defailt 1.0 and can be used to scale ais value, fore
! example to plot kJ rather than J for reasonable axis
     double precision, dimension(3) :: scalefact=one
! these define realative plot size for X and Y, normally 1.0 or less
     double precision :: xsize=1.0D0,ysize=1.0D0
! labeldefaults(i) for axis i 0 means default text, 1 text in plotlabels
! tielines>0 means plot a tieline every tielines calculated equilibrium
     integer :: labeldefaults(3),linett=1,tielines=0
! plotlabel(1) is heading, 2 is x-axis text, 3 is y-axis text 
     character*64, dimension(3) :: plotlabels
! linetype is 0 for dashed lines, 1 for full lines
!     integer linestyle
     integer linetype
! if linepoints >0 plot a symbol at each linewp point
     integer :: linewp=0
! if true plot a triangular diagram (isothermal section)
     logical gibbstriangle
! the set key command in GNUPLOT specifies where the line id is written
! it can be on/off, placed inside/outside, left/right/center, top/bottom/center,
! and some more options that may be implemented later ...
     character labelkey*48, font*32
! filename is file to write the GNUPLOT command and data file
! appendfile is a file name that will be appended unless empty
     character filename*256,appendfile*256
! gnuplot terminals and keys, gnuselterm is selected terminal type (1..8)
     integer gnutermsel,gnutermax
     character*80 gnuterminal(8)
     character*8 filext(8)
     character*8 gnutermid(8)
! firstextlabel is a pointer to a list of text label(s) 
! to be written at a given position
     TYPE(graphics_textlabel), pointer :: firsttextlabel
! pltax are the state variables to be plotted
! NOT USED: pform is the graphics format replaced by gnutermsel
     character pltax(2)*24,pform*32
! The default and current ending of a plot
     character*12 :: plotenddefault='pause mouse '
     character plotend*36
! added 18.09.24 text at lower left corner
     character (len=6) :: lowerleftcorner='      '
! added to have larger axis texts and line titles
     integer:: textonaxis=0
! nothing special 0, other as stepspecial: 1=separate; 2=Scheil; 3=Tzero;
!                                          4=paraequil; 5=NPLE
! but stepseparate, Tzero and paraequil works without using this.
! For Schiel I am trying to change line color for different parts of the line
     integer :: specialdiagram=0
! many more options can easily be added when desired, linetypes etc
  end TYPE graphics_options
!\end{verbatim}
!
! fix status during mapping, normally 2 means fix (not used)
  integer, parameter :: MAPPHASEFIX=3
! OS dependent values NOT BITS
#ifdef notwin
  integer, parameter :: PLOTONWIN=0
#else
  integer, parameter :: PLOTONWIN=1
#endif
!-------------------------------------------------
! BITS for graphopt status word, do not use bit 0 and 1 ...
! these bits are very confused ...
! GRKEEP is set if graphics windows kept (does not really matter)
! GRNOTITLE is set if no title plotted 
! GRISOPLETH if plot is an isopleth (no tielines)
! GRTABLE list results in a CSV table
  integer, parameter :: GRKEEP=2,   GRNOTITLE=3,  GRISOPLETH=4, GRCSVTABLE=5
!--------------------------------------------------
! default for some colors
  character (len=6) :: monovariant='7CFF40'       ! this is light green
  character (len=6) :: tielinecolor='7CFF40'
! for trace
  logical :: plottrace=.FALSE.
! Using memory for stored equilibria to avoid memory crash
! Totalsaved includes all equilibria saved during multiple map 
  integer totalsavedceq
  integer, parameter :: maxsavedceq=1999
! To warn that some calculated lines are excluded from plot
  integer :: lines_excluded=0
!
!-------------------------------------------------
!
  type(plot_line), pointer :: lastplotline,plotline1
! set by user for globalcheck during STEP/MAP
  integer :: mapglobalcheck=0
! repeated errors
  integer :: repeatederr=0
!
!-------------------------------------------------
!
CONTAINS

  ! routines to calculate the diagrams
  include "smp2A.F90"

  ! routines to plot the diagrams
  include "smp2B.F90"

END MODULE ocsmp

