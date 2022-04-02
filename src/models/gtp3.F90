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
  use ocnum
  use metlib
  use ocparam
!
! for known unfinished/unchecked bugs and parallelization problems
! look for BEWARE
!
! Using open MP parallelization
!$ use OMP_LIB
!
! overall version number
  character (len=8), parameter :: version='  6.047 '
!
!
! data structure nor non-encrypted TP functions
!
  include "gtp3_dd1.F90"
!
! most global data structure definitions
!
  include "gtp3_dd2.F90"
!  
CONTAINS

! 1-5: initialization, how many, find things, get things, set things, 
include "gtp3A.F90"

! 12: enter data
include "gtp3B.F90"

! 10: list data
include "gtp3C.F90"

! 11: save and read from files
include "gtp3D.F90"

! 7: state variable manipulations
include "gtp3E.F90"

! 8-9: state variable functions, interactive things
include "gtp3F.F90"

! 13-15: status for things, unfinished things, internal stuff
include "gtp3G.F90"

! 16: Additions (magnetic and others)
include "gtp3H.F90"

! 6: calculate things
include "gtp3X.F90"

! 17-18: Grid minimizer and miscellaneous
include "gtp3Y.F90"

! 19: Assessment subroutine 
include "gtp3Z.F90"


END MODULE GENERAL_THERMODYNAMIC_PACKAGE

