!
!***************************************************************
! General Thermodynamic Package (GTP)
! for thermodynamic modelling and calculations
!
MODULE GENERAL_THERMODYNAMIC_PACKAGE
!
! Copyright 2011-2021, Bo Sundman, France
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
! contact person: bo.sundman@gmail.com
!
!-----------------------------------------------------------------------
!
!
! for known unfinished/unchecked bugs and parallelization problems
! look for BEWARE
!
! Using open MP parallelization (also added to metlib4.F90 for error code)
!$ use OMP_LIB
!
  use ocnum
  use metlib
  use ocparam
!
!! overall OC version number
  character (len=8), parameter :: version='  6.098 '
!
!
! data structure for non-encrypted TP functions
!
! use #include rather than include to have preprocessor options
!  include "gtp3_dd1.F90" for TP functions without decrypted databases
#include "gtp3_dd1.F90"
!
! most global data structure definitions
!
!  include "gtp3_dd2.F90" all other data structures
#include "gtp3_dd2.F90"
!
! XML elements and attributes
#include "gtp3_xml.F90"
!
CONTAINS

! 1-5: initialization, how many, find things, get things, set things, 
!include "gtp3A.F90"
#include "gtp3A.F90"

! 12: enter data
!include "gtp3B.F90"
#include "gtp3B.F90"

! 10: list data
!include "gtp3C.F90"
#include "gtp3C.F90"

! 11: save and read from files
!include "gtp3D.F90"
#include "gtp3D.F90"

! 9A: Read/write TDB/UNFORMATTED
!include "gtp3E.F90"
#include "gtp3E.F90"

! 9B: Read/write XML
!include "gtp3EX.F90"
#include "gtp3EX.F90"
#include "gtp3EY.F90"

! 7-8: state variable functions, interactive things
!include "gtp3F.F90"
#include "gtp3F.F90"

! 13-15: status for things, unfinished things, internal stuff
!include "gtp3G.F90"
#include "gtp3G.F90"

! 16: Additions (magnetic and others)
!include "gtp3H.F90"
#include "gtp3H.F90"

! 6: calculate things, gtp3XQ for MQMQA
!include "gtp3X.F90"
#include "gtp3X.F90"
#include "gtp3XQ.F90"

! 17-18: Grid minimizer and miscellaneous
!include "gtp3Y.F90"
#include "gtp3Y.F90"

! 19: TPFUN routines for non-encrypted databases
!include "gtp3Z.F90"
#include "gtp3Z.F90"

END MODULE GENERAL_THERMODYNAMIC_PACKAGE

