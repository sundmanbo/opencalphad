!***************************************************************
! General Thermodynamic Package (GTP)
! for thermodynamic modelling and calculations
!
MODULE GENERAL_THERMODYNAMIC_PACKAGE
!
! Copyright 2011-2013, Bo Sundman, France
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
! for known unfinished/unchecked bugs and parallelization problems
! look for BEWARE
!
!-----------------------------------------------------------------------
!
! Description of data structure
!
! For all elements, species and phases there are two arrays defined.
! The first (data) array contains the elements species etc and all their data
! in the order they were entered and the data are never moved.
! The second (index) array contain the elements, species etc in alphabetical 
! (or whatever) order and is updated whenever a new element, species etc 
! is added. This array is an integer array with the index of the data array.
! Most links inside the different records to elements, species etc
! are indices to the data array which is never changed.  
! TPFUNS used in parameters are also stored in an array and the index
! to this array is stored in the property record to specify the function.
!
! For parameters inside each phase record there is one (or 2 if disordered set)
! lists with endmember parameters.  Each endmember record can be the root
! of a binary tree with interaction parameters.  Each of these records
! can have a property list with various data like G, TC, MQ etc.  These
! records are created dynamically and can only be found by following the links.
!
! Each phase has one or more composition sets.  These are part of the
! equilibrium data structure which also contains conditions, calculated
! values of TP functions and other symbols.
!
! An equilibrium record has been introduced.  One such record is created
! in init_gtp and it is called FIRSTEQ which is a global variable.
! There is also an array EQLISTA which should contain all allocated
! equilibrium records, FIRSTEQ is a pointer to the first element in this array.
! More equilibrium records can be initiated by enter_equilibrium subroutine.
! This copies the relevant data from FIRSTEQ.
! After a second equilibrium is created it is forbidden to enter elements,
! species? or phases and create additional fraction sets, i.e. one must
! not change the data structure except to add/remove composition sets
! (but not implemented yet for multiple equilibria).  Composition sets must be
! created in all equilibrium records at the same time (if done in a thread
! then all threads must stop while this is done).  During step/map calculation
! each calculated equilibria is saved for later use in plotting och other
! postprocessing.  These saved equilibria may have different number of
! composition sets so great care must be taken using them.
!
! The equilibrium data record is "stand alone" and contains all necessary
! data to describe the equilibrium (except the model parameters and other
! static data).  In parallel processing each thread will have its own
! equilibrium data record.
!
! The intention is that several equilibra can be created both to store 
! individual experimental data in assessments and for each thread in parallel.
! In the equilibrium record there are conditions, components (with chemical
! potentials) and an error code and most important, the phase_varres record
! array with one or more record for each phase.  This array must be identical
! in all equilibria recods.  Each composition set has a phase_varres record
! and they are linked from the phase record by the LIKTOCS array.  As the
! phase_varres records are in an array the link is simply an integer index of
! this array.  There is a free list inside the phase_varres array to be used
! when adding or removing a composition set.
! >>> NOT DONE: when composition sets are created they must be
!               created in all equilibrium records at the same time.
! The EQ_TPRES array is declared inside the equilibrium record.
! The index to a function in EQ_TPRES is the same as the index to the TPFUN
! array declared globally in TPFUN.  The TPFUN array has the actual expression,
! and EQ_TPRES has the last calculated results, which can be different in each
! equilibrium.  The TPFUN index is used in property records to specify
! the function of a parameter.
! 
! In many subroutines the equilibrium record called CEQ (Current EQilibrium)
! is an argument which means it operates on the data in that 
! equilibrium record only.
!
! In the record array PHASE_VARRES (including disordered) each phase 
! and composition set has a record.  If no parallel calculation and no
! experiments the equilibrium record FIRSTEQ is enough.
!
! In programming for parallel processing THREADPRIVATE
! should be avoided as it usually has a very slow implementation.
!
! Some routines exist both with and without the CEQ argument.  A programmer
! can create his own array of equilibrium data records and use any of
! them in such calls. ???? Maybe not, then how to update when a new
! composition set is needed???
!
! Thread specific data are needed for conditions, phase status, constitution,
! function values and calc results like G and derivatives for each phase,
! amounts of phases etc.  When calling a subroutine to get mole frations etc 
! the equilibrium record CEQ must be supplied.
!
! The global error code is defined in tpfunlib, that is not very good.  There
! must be an error code specific to each equilibrium.  Or can one declare
! the error code as "local" to the thread?
!
! routines for inverting matrix etc
  use lukasnum
! TP functions, tpfunlib makes USE of metlib
  use tpfunlib
! for parallel processing, not used yet
!  use OMP_LIB
!
!--------------------------------------------------------------------------
!
! To be added or fixed (in no special order, * means priority):
!  1 reference states for G, H, MU etc
!  2* ionic liquid model (incl normal quasichemical)
!  3 IRSID slag model (not by me)
!  4 volume model
!  5* parameter permutations for ordering (option B, done for F)
!  6 consistent units conversion between user/software (C/F/K,cal/BTU/eV/J etc)
!  7 reciprocal composition dependent parameters
!  8* HSS corrected quasichemical liquid model
!  9 FACT modified quasichemical liquid model
! 10 multicomponent CVM tetrahedron fcc model
! 11 dynamic components? components for each phase? The 2H2+O2=2H2O case
! 12* wildcards in state variables (like show x(*), partially done)
! 13 Implement SER phase data.  When calculating with the SER phase one must
!    have different magnetic models for different elements !!
! 14 amend_data for elements, species, phases, symbols
! 15 more symbols (variables, value-functions, PUTFUN functions, coefficients)
! 16 Find the reason for the difference with TC for the magnetic enthalpy
! 17 New magnetic model (Wei)
! 18 New heat capacity models for low T (from Mauro)
! 19 Enable that state variables functions can call TP functions
! 20 Improve the grid for gridminimization
! 21* fix the error code for parallel processing
! 22 save/read of data
! 23 parallellisation
! 24 Modify gcalc so one can calcule just a single property like the mobility
!
!-------------------------------------------------------------------------
!
! Outside this package
! 1 Decide output format for plotting (gnuplot/excell)
! 2 Step procedure
! 3 TQ interface
! 4 Map procedure
! 5 Assessment procedure
! 6 Retire
!
!---------------------------------------------------------------------------
! Currently working on
! TQ interface
! 
!---------------------------------------------------------------------------
! To be fixed
! - Phase status number and type, not same when changing and getting
! - Composition sets with same name (but different from the original) do
!   not require a #digit to be unique, that is wrong
!
!---------------------------------------------------------------------------
!
! Done so far
!  1 enter elements, species, phases, TP-functions, parameters
!  2 Calculate G for CEF phase with binary RK parameters and tabulate data
!  3 calculate with ternary composition dependent parameters
!  4 separation of data for parallell processing
!  5 read a TDB file (except a few things and with some restrictions)
!  6 Save/read data on file with new format, formatted and unformatted
!  7 TC and BM parameters calculated separately
!  8 magnetic model implemented.
!  9 primitive user i/f
! 10 more types of parameters (incl mobilities) possible
! 11 composition sets
! 12 some status bits (suspended/fix etc)
! 13 some status changes implemented
! 14 entering data in PMLFKT arrays while reading from a TDB file
! 15 routines to extract a value of a state variable including summed over all
!    phases for use in conditions, show and assessments
! 17 A special reference phase (SER) with index 0 that is the default reference
!    for all elements. (It must be a phase and not just a tpfun as magnetism)
!    All entered elements can have a set of parameters there.
! 18 wildcard parameters like G(fcc,*), L(fcc,a,b:*), L(fcc,*,a:b)
! 19 enter parameters for disordered fraction sets
! 20 calculations with several fraction sets in each phase (ordering)
!    both for sigma and fcc
! 21 references for parameters fixed
! 22 expanded primitive user i/f
! 23 Gridminimization part 1 done
! 24 Found a bug in TC for wildcard parameters (second derivatives)
! 25 Found a bug in TC for partitioned phases (second derivatives)
! 26 Added an array for max/min constituent fractions in phase_varres record
!    this will be used to select which composition set sould be used for a
!    given constitution.  A negative value means a min value, positive a max
! 27 Global gridminimization done (MISSING: do not find and store the 
!    best constitution of metastable phases)
! 28 Managed parallel gridminimization calculating each phase 
!    in separate threads for fixed T and P and one equilibrium record.
! 29 Created equilibrium record and made the PHASE_VARRES record part of
!    this.  A pointer CEQ (Current EQuilibrium) is passed to those subroutines
!    that need to know which equilibrium record data to operate on.
! 30 Moved TPFUN to a separate module so it is easier to rewrite
! 31 tpfun_parres array is made part of the equilibrium record. 
! 32 Entering, listing and finding conditions including expressions.
!    Using the condition for gridminimization
! 33 Added state variable functions (PUTEXP) like in TC POLY.
!    Some of these functions can be associated to a specific equilibrium
!    so one can transfer values between equilibria (using AMEND command).
! 34 Modified PD_LEVER/_INIT so it can understand some conditions set by GTP
! 35 Transform some conditions to PMLFKT format
! 36 Creating new equilibria seems to work.  Must handle disordered
!    fraction sets also with several equilibrium records and in parallel ...
! 37 Can handle new composition sets when two or more equilibria
! 38 FIX phase status now added to conditions
! 39 Extended the use of parameter references and corrected some errors
! 40 TPFUN and PUTFUN allocate dynamically their internal structures.
! 41 the components list moved to the equilibrium record and
!    user defined components possible (for each equilibrium)
! 42 replaced arrays and free lists for property, interaction and endmember
!    records with dynamic lists.  Keept arrays for elements, species, phases,
!    sublattices, phase_varres, tpfuns and statevar symbols.
! 43 Allow using condition number when changing or removing it, like 
!    set cond 1:=none. Allow set cond *=none (does not change phase 
!    conditions of fix phases)
! 44 Save/read unformatted implemented again (with some bugs)
! 45 Equilibrium results copied to gtp results arrays.  
! 46 Metastable phases from a gridcalculation have DGM and their 
!    best composition calculated and stored in GTP result arrays.
! 47 Made easier to extend the number of property parameter (G, TC, BM etc)
!    and for additions which properties they need (magnetism need TC, BM)
! 48 PD_LEVER now tests all phases to find the stable phase set (buggy)
! 49 Fixed constituent related property parameter types like mobilities
! 50 Some parameters like TC, BMAGN etc can be constants of just depend on
!    P or just on T.  This can be specified as bits.
! 51 Composition sets can be created when necesssary after a gridgridmin
!    Other threads must not use EQCALC when a composition set is created
! 52 eqcalc2 can now calculate equilibria with massbalance conditions but
!    need gridminimizer for good start values.  It needs fine-tuning for many
!    cases and fixing default start constitutions.
! 53 Changed most arrays declared as pointers to be declared allocatable.
! 54 Step with T axis works.
! 55 Plot with gnuplot works.
! 56 Default start constitutions, upper/lower limits of fractions
! 57 Wildcards for state variables for plotting
! 58 included sublattice structure inside phase structure, no SUBLISTA array
! 59 Added set-input-amounts
! 60 Changed grid generation in grid minimmizer, still not very good
! 61 Chnaged to IMPLICIT NONE
! 62 Fixed bug in state variable for amounts like N(FE)
! 63 Moved some status bits from globaldata to gtp_equilibrium_data
! 64 Added an array LINKTOCS in gtp_phase with links to all compsets. The NEXT
!    link in phase_varres removed.  This should simplify use of compsets.
! 65 Fixed so all parameter property symbols can be accessed as state variables
! 66 Most permutations for FCC done for entering and calculations (puh)
! 67 Added a lot ot parameter property symbols
! 68 Added a rudimentary on-line help facility
! 69 The phase status/bit setting/checking improved
! 70 Added a sequential number in all endmember and interaction records.
! 71 Split the error message data not to excees max 256 continuation
!
!=================================================================
!
!  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  IMPLICIT NONE
!
!-----------------------------------------------------------------
! error messages
! numbers 4000 to 4220 defined.  gx%bmperr is set to message index
  integer, parameter :: nooferm=4220
  character (len=64), dimension(4000:nooferm) :: bmperrmess
  data bmperrmess(4000:4199)&
      /'Too many coefficients in a TP function.                         ',&
       'Illegal character in a TP function, digit expected.             ',&
       'Unkown symbol in TP function                                    ',&
       'Expected ( after unary function                                 ',&
       'Too many ) in a TP function                                     ',&
       'Illegal character in a TP function                              ',&
       'Too few ) in a TP function                                      ',&
       'Too many ( in exponent                                          ',&
       'Illegally placed ( in the exponent of a TP function             ',&
       'No digits after ( in the exponent of a TP function              ',&
! 4010:
       'Illegally placed ) in exponent in TP function                   ',&
       'Too high power in a TP function, max 99, min -99                ',&
       'Missing ) in the exponent of a TP function                      ',&
       'Illegal termination of a TP function reading a TDB file         ',&
       'No more free TP root records                                    ',&
       'No more free TP expression records                              ',&
       'Illegal expression inside unary argument of a TP function       ',&
       'Illegal code found when evaluating a TP function                ',&
       'Found a coefficent zero in a term of a TP function              ',&
       'Illegal code in a TP function                                   ',&
! 4020:
       'Negative argument to logarithm in a TP function                 ',&
       'Unknown unary function in evaluation for a TP function          ',&
       'Too many symbols in a TP function term                          ',&
       'Two unary functions in a TP function term                       ',&
       'Too complicated TP function term                                ',&
       'Too many temperature ranges in a TP function                    ',&
       'TP function with same name already entered                      ',&
       'Symbol referenced in a parameter does not exist                 ',&
       'Missing separator between phase and constituent array in paramet',&
       'Cannot enter disordered fraction set when several composition se',&
! 4030
       'Cannot enter disordered fraction set when suspended constituents',&
       'Wildcards in interaction parameters not yet implemented         ',&
       'Interaction between 2 wildcards are illegal                     ',&
       'Illegal character in element symbol                             ',&
       'Element with this symbol already entered                        ',&
       'Element symbol and name must start with letter A-Z              ',&
       'Reference state must start with letter A-Z                      ',&
       'Element mass must not be negative                               ',&
       'Enthalpy difference H298-H0 must be positive                    ',&
       'Entropy at 298.15 must be positive                              ',&
! 4040
       'Too many elements                                               ',&
       'Too many species                                                ',&
       'No such element                                                 ',&
       'Text position outside text                                      ',&
       'Species symbol must start with letter A-Z                       ',&
       'No elements or too many elements in species formula             ',&
       'Unknown element in species formula                              ',&
       'Negative stoichiometric factor in species                       ',&
       'The charge must be the final "element"                          ',&
       'Species already entered                                         ',&
! 4050
       'No such phase                                                   ',&
       'No such species                                                 ',&
       'No such component                                               ',&
       'Phase name must start with letter A-Z                           ',&
       'Phase already entered                                           ',&
       'Model not implemented yet                                       ',&
       'Too few or too many sublattices                                 ',&
       'Sites on a sublattice must be positive                          ',&
       'Too few or too many constituents in a sublattice                ',&
       'Too many constituents                                           ',&
! 4060
       'No such TP function                                             ',&
       'Expected constituent array, found nothing                       ',&
       'Illegal character in constituent array                          ',&
       'Illegal degree of parameter, must be 0-9                        ',&
       'No free interaction records                                     ',&
       'Wrong number of sublattices                                     ',&
       'No such constituent in a sublattice                             ',&
       'No such interacting constituent                                 ',&
       'This phase has no disordered fraction set                       ',&
       'Wrong number of sublattices in disordered fraction set          ',&
! 4070
       'No free endmember records                                       ',&
       'No free property records                                        ',&
       'No such composition set                                         ',&
       'Inconsistent composition set specifications                     ',&
       'Overflow in push                                                ',&
       'Undeflow in pop                                                 ',&
       'Sublattice out of range for entering disordered fraction set    ',&
       'Disordered fraction set already entered                         ',&
       'Not implemented yet                                             ',&
       'Ionic liquid Not implemented yet                                ',&
! 4080
       'Suspended constituents not implemented yet                      ',&
       'Stability factor not implemented yet                            ',&
       'No such composition dependent property parameter                ',&
       'Empty line, expected species stoichiometry                      ',&
       'No element in species stoichiometry                             ',&
       'Species cannot be entered as it is implicitly suspended         ',&
       'Excess model not implemented yet                                ',&
       'Bad name for a symbol                                           ',&
       'Too deeply nested TP functions                                  ',&
       'Reading unkown addition type from file                          ',&
! 4090
       'Addition already entered                                        ',&
       'No more addition records                                        ',&
       'Maximum 9 composition sets                                      ',&
       'Illegal composition set number                                  ',&
       'No more records for phases or composition sets.                 ',&
       'Hidden phase cannot be ENTERED, SUSPENDED, DORMANT or FIXED     ',&
       'No such constituents                                            ',&
       'Too many argument to a state variable                           ',&
       'This state variable must have two arguments                     ',&
       'First character of a state variable is wrong                    ',&
! 4100
       'State variable starting with M not followed by U                ',&
       'State variable starting with L not followed by NAC              ',&
       'Missing ( for arguments for state variable                      ',&
       'Missing ) after arguments of state variable                     ',&
       'Unknown phase used as state varible argument                    ',&
       'Unknown constituent used as state variable argument             ',&
       'Unknown component used as state variable argument               ',&
       'State variable starting with D not followed by G                ',&
       'State variable starting with T follwed by other character than C',&
       'State variable starting with B missing P, MAG, M, V, W or F     ',&
! 4110
       'This state variable cannot not have two arguments               ',&
       'This state variable must have an argument                       ',&
       'Impossible reference state for this constituent                 ',&
       'No such property for this phase                                 ',&
       'Cannot calculate property value per volume as no volume data    ',&
       'Property per formula unit only for a single phase               ',&
       'State variable number must be larger than zero                  ',&
       'Only state variable Y can have 3 indices                        ',&
       'Illegal normalization of state variable                         ',&
       'Phase is hidden                                                 ',&
! 4120
       'Wrong syntax for mobility variable                              ',&
       'Ambiguous phase name                                            ',&
       'Illegal name for an equilibrium                                 ',&
       'Equilibrium with this name already entered                      ',&
       'No such equilibrium                                             ',&
       'Not allowed to enter more model data                            ',&
       'No state variable supplied                                      ',&
       'Illegal state variable for conditions                           ',&
       'Only one kind of state variable in expressions                  ',&
       'Illegal value for a state variable                              ',&
! 4130 line below
       'Factor in front of a condition must be followed by *            ',&
       'No such condition                                               ',&
       'Function name must start with a letter A-Z                      ',&
       'Function name and expression must be separated by "="           ',&
       'Error in function expression (putfun)                           ',&
       'Unknown symbol used in function                                 ',&
       'Symbol with this name already entered                           ',&
       'Symbol name must start with letter A-Z                          ',&
       'Illegal character in symbol name                                ',&
       'Cannot check name of unknown kind of symbol                     ',&
! 4140
       'No such symbol                                                  ',&
       'Error evaluating symbol value                                   ',&
       'Error listing symbol expression                                 ',&
       'No conditions at all                                            ',&
       'Degrees of freedom not zero                                     ',&
       'Unknown type of addition                                        ',&
       'Quitting due to repeated input error                            ',&
       'Gridminimizer found gridpoint outside range                     ',&
       'Gridminimizer error when generating endmember values            ',&
       'Gridminimizer found an element without gridpoint                ',&
! 4150 next line
       'Gridminimizer have no gridpoint for a pure element              ',&
       'Conditions not only T, P and massbalance                        ',&
       'Illegal to set all phases as fix                                ',&
       'Cannot enter a new equilibrium if there are no phases           ',&
       'Trying to enter an illegal reference                            ',&
       'A reference must have an identifier                             ',&
       'Reference identifier already exists                             ',&
       'Error in TDB file, species terminator error                     ',&
       'Unknown potential                                               ',&
       'Cannot calculate potentials for charged constituents            ',&
! 4160 next line
       'Illegal endmember for reference state                           ',&
       'End member without atoms                                        ',&
       'Same species twice in component list                            ',&
       'Component stoichiometry matrix singular                         ',&
       'Too many interaction levels                                     ',&
       'Error reading save file                                         ',&
       'Error reading save file at EOF                                  ',&
       'Composition set prefix must start with a letter                 ',&
       'This property has no specifier                                  ',&
       'Parameter specifier missing                                     ',&
! 4170
       'Properties needed for Inden magnetic model not defined          ',&
       'Request for non-existing chemical potential                     ',&
       'Removing current data not implemented                           ',&
       'Grid minimization not allowed                                   ',&
       'Illegal composition in call to grid minimization                ',&
       'Too many gridpoints                                             ',&
       'No phases and no gridpoints in call to grid minimization        ',&
       'Grid minimizer want to create composition set but is not allowed',&
       'Non-existing fix phase                                          ',&
       'State variable N or B cannot have two indices for grid minimizer',&
! 4180
       'Condition on B is not allowed for grid minimizer                ',&
       'Element has no composition in grid minimizer                    ',&
       'Too complicated mass balance conditions                         ',&
       'Two mass balance conditions for same element                    ',&
       'Cannot handle conditions on both N and B                        ',&
       'No mole fractions when summing composition                      ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
! 4190
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                '/
! 4200 errors in minimizer
  data bmperrmess(4200:4220)&
      /'No phase that can be set stable                                 ',&
       'Attempt to set too many phaes as stable                         ',&
       'Total amount is negative                                        ',&
       'Error solving system matrix                                     ',&
       'Too many iterations                                             ',&
       'Phase matrix singular                                           ',&
       'Cannot handle models without analythical second derivativatives ',&
       'Ionic liquid model not implemented yet                          ',&
       '                                                                ',&
       '                                                                ',&
! 4210
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
! 4220
       '                                                                '/
! last used error codes above
!
!=================================================================
!
!\begin{verbatim}
!-Bits in global status word (GS) in globaldata record
! level of user: beginner, occational, advanced; NOGLOB: no global gridmin calc
! NOMERGE: no merge of gridmin result, NODATA: not any data, 
! NOPHASE: no phase in system, NOACS: no automatic creation of composition set
! NOREMCS: do not remove any redundant unstable composition sets
! NOSAVE: data changed after last save command
! >>>> som of these should be moved to the gtp_equilibrium_data record
  integer, parameter :: &
       GSBEG=0,     GSOCC=1,     GSADV=2,     GSNOGLOB=3, &
       GSNOMERGE=4, GSNODATA=5,  GSNOPHASE=6, GSNOACS=7, &
       GSNOREMCS=8, GSNOSAVE=9
!-Bits in element record
  integer, parameter :: &
       ELSUS=0
!-Bits in species record
! Suspended, implicitly suspended, species is element, species is vacancy
! species have charge, species is component
  integer, parameter :: &
       SPSUS=0, SPIMSUS=1, SPEL=2, SPVA=3, &
       SPION=4, SPSYS=5
!\end{verbatim}
!\begin{verbatim}
!-Bits in phase record
! hidden, implictly hidden, ideal, no concentration variation (NOCV),
! Phase has parameters entered (PHHASP), 
! F option (FORD), B option (BORD), Sigma ordering (SORD),
! multiple/disordered fraction sets (MFS), gas, liquid, ionic liquid, 
! aqueous, dilute config. entropy (DILCE), quasichemical (QCE), CVM,
! FACT,  not create comp. sets (NOCS), Helmholz energy model (HELM),
! Model without 2nd derivatives (PHNODGDY2), Elastic model A,
! explicit charge balance needed (XCB),
  integer, parameter :: &
       PHHID=0,     PHIMHID=1,  PHID=2,    PHNOCV=3, &     ! 1 2 4 8
       PHHASP=4,    PHFORD=5,   PHBORD=6,  PHSORD=7, &
       PHMFS=8,     PHGAS=9,    PHLIQ=10,  PHIONLIQ=11, &   
       PHAQ1=12,    PHDILCE=13, PHQCE=14,  PHCVMCE=15,&
       PHFACTCE=16, PHNOCS=17,  PHHELM=18, PHNODGDY2=19,&
       PHELMA=20,   PHEXCB=21
! 
!-Bits in constituent fraction (phase_varres) record STATUS2
! CSDFS is set if record is for disordred fraction set, then one must use
!     sublattices from fraction_set record
! CSDLNK: a disordred fraction set in this phase_varres record
! CSSUS: set if comp. set if must not be stable, 
! CSFIXDORM: set if fix or dormant, 
! CSCONSUS set if one or more constituents suspended (status array constat
!     specify constituent status)
! CSORDER: set if fractions are ordered (only used for BCC/FCC ordering
!     with a disordered fraction set).
! CSSTABLE: set if phase is stable after an equilibrium calculation
! CSAUTO set if composition set created during calculations
! CSDEFCON set if there is a default constitution
! NOTE phase_status ENTERED means both CSSUS and CSFIXDORM are not set
   integer, parameter :: &
        CSDFS=0,    CSDLNK=1,  CSSUS=2,    CSFIXDORM=3, &
        CSCONSUS=4, CSORDER=5, CSSTABLE=6, CSAUTO=7, &
        CSDEFCON=8
!\end{verbatim}
!\begin{verbatim}
!-Bits in constat array for each constituent
! For each constituent: is suspended, is implicitly suspended, is vacancy
   integer, parameter :: &
        CONSUS=0,CONIMSUS=1,CONVA=2
!-Bits in state variable functions (svflista)
! SVFVAL symbol evaluated only explicitly (mode=1 in call)
   integer, parameter :: &
        SVFVAL=0,SVFEXT=1
!-Bits in gtp_equilibrium_data record
! EQNOTHREAD set if equilibrium must be calculated before threading 
! (in assessment) for example if a symbol must be evaluated in this 
! equilibrium before used in another like H(T)-H298
   integer, parameter :: &
        EQNOTHREAD=0, EQNOGLOB=1, EQNOEQCAL=2, EQINCON=3, &
        EQFAIL=4,     EQNOACS=5,  EQGRIDTEST=6
!-Bits in parameter property type record (gtp_propid)
! constant (no T or P dependence), only P, property has an element suffix
! (like mobility), property has a constituent suffix
   integer, parameter :: &
        IDNOTP=0, IDONLYP=1, IDELSUFFIX=2, IDCONSUFFIX=3
!- Bits in condition status word (some set in onther ways??)
! singlevar means T=, x(el)= etc, singlevalue means value is a number
! phase means the condition is a fix phase
  integer, parameter :: &
       ACTIVE=0,SINGLEVAR=1,SINGLEVALUE=2,PHASE=3
!
! >>> Bits for symbols and TP functions missing ???
!\end{verbatim}
!
!----------------------------------------------------------------------
!
!\begin{verbatim}
! Parameters defining the size of arrays etc.
! max elements, species, phases, sublattices, constituents (ideal phase)
 integer, parameter :: maxel=20,maxsp=1000,maxph=200,maxsubl=10,maxconst=1000
! maximum number of consitutents in non-ideal phase
 integer, parameter :: maxcons2=100
! maximum number of elsements in a species
 integer, parameter :: maxspel=10
! maximum number of references
 integer, private, parameter :: maxrefs=1000
! maximum number of equilibria
 integer, private, parameter :: maxeq=100
! some dp values, default precision of Y and default minimum value of Y
! zero and one set in tpfun
 double precision, private, parameter :: YPRECD=1.0D-6,YMIND=1.0D-30
! dimension for push/pop in calcg, max composition dependent interaction
 integer, private, parameter :: maxpp=1000,maxinter=3
! max number of symbols
 integer, private, parameter :: maxtpf=10*maxph
! max number of properties (G, TC, BMAG MQ%(...) etc)
 integer, private, parameter :: maxprop=50
! max number of state variable functions
 integer, private, parameter :: maxsvfun=500
! version number
! changes in last 2 digits means no change in SAVE/READ format
 character*8, parameter :: gtpversion='GTP-1.00'
!\end{verbatim}
!=================================================================
!\begin{verbatim}
! The number of additions to the Gibbs energy of a phase is increasing
! This is a way to try to organize them.  Each addtion has a unique
! number identifying it when created, listed or calculated.  These
! numbers are defined here
 integer, public, parameter :: indenmagnetic=1
 integer, public, parameter :: debyecp=2
 integer, public, parameter :: weimagnetic=3
 integer, public, parameter :: einsteincp=4
 integer, public, parameter :: elasticmodela=5
 integer, public, parameter :: glastransmodela=6
! Note that additions often use parameters like Curie or Debye temperatures
! defined by parameter identifiers stored in gtp_propid
!\end{verbatim}
!=================================================================
!
! below here are data structures and global data in this module
!
!=================================================================
!\begin{verbatim}
 TYPE gtp_global_data
! status should contain bits how advanced the user is and other defaults
! it also contain bits if new data can be entered (if more than one equilib)
    integer status
    character name*24
    double precision rgas,rgasuser,pnorm
 END TYPE gtp_global_data
 TYPE(gtp_global_data) :: globaldata
!\end{verbatim}
!=================================================================
!
! below here are thermodynamic model data structures
!
!=================================================================
!
!\begin{verbatim}
 TYPE gtp_element
! data for each element: symbol, name, reference state, mass, h298-h0, s298
    character :: symbol*2,name*12,ref_state*24
    double precision :: mass,h298_h0,s298
! splink: index of corresponing species in array splink
! Status bits are stored in the integer status
! alphaindex: the alphabetical order of this elements
! refstatesymbol: indicates H0 (1), H298 (0, default) or G (2) for endmembers
    integer :: splink,status,alphaindex,refstatesymbol
 END TYPE gtp_element
! allocated in init_gtp
 TYPE(gtp_element), private, allocatable :: ellista(:)
 INTEGER, private, allocatable :: ELEMENTS(:)
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_species
! data for each species: symnol, mass, charge, status
! mass is in principle redundant as calculated from element mass
    character :: symbol*24
    double precision :: mass,charge
! alphaindex: the alphabetical order of this species
! noofel: number of elements
    integer :: noofel,status,alphaindex
! Use an integer array ellinks to indicate the elements in the species
! The corresponing stoichiometry is in the array stochiometry
! ???? these should not be pointers, changed to allocatable ????
    integer, dimension(:), allocatable :: ellinks
    double precision, dimension(:), allocatable :: stoichiometry
 END TYPE gtp_species
! allocated in init_gtp
 TYPE(gtp_species), private, allocatable :: splista(:)
 INTEGER, private, allocatable :: SPECIES(:)
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_components
! The components are simply an array of indices to species records
! the components must be "orthogonal".  There is always a "systems components"
! that by default is the elements.
! Later one may implement that the user can define a different "system set"
! and also specific sets for each phase.
! The reference state is set as a phase and value of T and P.
! The name of the phase and its link and the link to the constituent is stored
! the endmember array is for the reference phase to calculate GREF
! The last calculated values of the chemical potentials (for user defined
! and default reference states) should be stored here.
    integer :: splink,phlink,status
    character*16 :: refstate
    integer, dimension(:), allocatable :: endmember
    double precision, dimension(2) :: tpref
    double precision, dimension(2) :: chempot
    double precision mass
 END TYPE gtp_components
! allocated in gtp_equilibrium_data
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_endmember
! end member parameter record, note ordered phases can have
! several permutations of fraction pointers like for B2: (Al:Fe) and (Fe:Al).
! There are links (i.e. indices) to next end member and to the interactio tree
! and to a list of property record
! The phase link is needed for SAVE/READ as one cannot know the number of
! sublattices otherwise.  One could just store nsl but a link back to the
! phase record might be useful in other cases.
! noofpermut: number of permutations (for ordered phases: (Al:Fe) and (Fe:Al)
! phaselink: index of phase record
! antalem: sequenial order of creation, useful to keep track of structure
! propointer: link to properties for this endmember
! nextem: link to next endmember
! intponter: root of interaction tree of parameters
! fraclinks: indices of fractions to be multiplied with the parameter
    integer :: noofpermut,phaselink,antalem
    TYPE(gtp_property), pointer :: propointer
    TYPE(gtp_endmember), pointer :: nextem
    TYPE(gtp_interaction), pointer :: intpointer
! there is at least one fraclinks per sublattice
! the second index of fraclinks is the permutation (normally only one)
! the first indec of fraclinks points to a fraction for each sublattice.
! The fractions are numbered sequentially independent of sublattices, a
! sigma phase with (FE:CR,MO:CR,FE,MO) has 6 fractions (incl one for FE in
! first sublattice) and the end member (FE:MO:CR) has the fraclinks 1,3,4
! This means these values can be used as index to the array with fractions.
! The actual species can be found via the sublattice record
!    integer, dimension(:,:), pointer :: fraclinks
    integer, dimension(:,:), allocatable :: fraclinks
 END TYPE gtp_endmember
! dynamically allocated when entering a parameter
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_interaction
! this record constitutes the parameter tree. There are links to NEXT
! interaction on the same level (i.e. replace current fraction) and
! to HIGHER interactions (i.e. includes current interaction)
! There can be several permutations of the interactions (both sublattice
! and fraction permuted, like interaction in B2 (Al:Al,Fe) and (Al,Fe:Al))
! The number of permutations of interactions can be the same, more or fewer
! comparaed to the lower order parameter (endmember or other interaction).
! The necessary information is stored in noofip.  It is not easy to keep
! track of permutations during calculations, the smart way to store the last
! permutation calculated is in this record ... but that will not work for
! parallell calsulations ...
! status: may be useful eventually
! antalint: sequential number of interaction record, to follow the structure
! order: for permutations one must have a sequential number in each node
! propointer: link to properties for this parameter
! nextlink: link to interaction on same level (replace interaction)
! highlink: link to interaction on higher level (include this interaction)
! sublattice: (array of) sublattices with interaction fraction
! fraclink: (array of) index of fraction to be multiplied with this parameter
! noofip: (array of) number of permutations, see above.
    integer status,antalint,order
    TYPE(gtp_property), pointer :: propointer
    TYPE(gtp_interaction), pointer :: nextlink,highlink
    integer, dimension(:), allocatable :: sublattice,fraclink,noofip
 END TYPE gtp_interaction
! allocated dynamically and linked from endmember records
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_property
! This is the property record.  The end member and interaction records
! have pointer to this.  Severall different properties can be linked
! from a parameter record like G, TC, BMAGN, VA, MQ etc.
! Some properties are connected to a constituent (or component?) like the
! mobility and also the Bohr mangneton number.
! Allocated as linked from endmembers and interaction records
! reference: can be used to indicate the source of the data
! refix: can be used to indicate the source of the data
! nextpr: link to next property record
! extra: TOOP and KOHLER can be implemented inside the property record
! proptype: type of propery, 1 is G, other see parameter property
! degree: if parameter has Redlich-Kister or similar degrees (powers)
! degreelink: indices of TP functions for different degrees (0-9)
! protect: can be used to prevent listing of the parameter
! antalprop: probably redundant (from the time of arrays of propery records)
    character*16 reference
    TYPE(gtp_property), pointer :: nextpr
    integer proptype,degree,extra,protect,refix,antalprop
    integer, dimension(:), allocatable :: degreelink
 END TYPE gtp_property
! property records, linked from endmember and interaction records, allocated
! when needed.  Each propery like G, TC, has a property record linking
! a TPFUN record (by index to tpfun_parres)
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_datareference
! store data references
! reference: can be used for search of reference
! refspec: free text
    character*16 reference
    character*64, dimension(:), allocatable :: refspec
 END TYPE gtp_datareference
! allocated in init_gtp
 TYPE(gtp_datareference), private, allocatable :: reflista(:)
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_propid
! this identifies different properties that can depend on composition
! Property 1 is the Gibbs energy and the others are usually used in
! some function to contribute to the Gibbs energy like TC or BMAGN
! But one can also have properties used for other things like mobilities
! with additional especification like MQ&FE
! symbol: property identifier like G for Gibbs energy
! note: short description for listings
! prop_elsymb: additional for element dependent properties like mobilities
    character symbol*4,note*16,prop_elsymb*2
! Each property has a unique value of idprop.  Status can state if a property
! has a constituent specifier or if it can depend on T or P
    integer status
! this can be a constituent specification for Bohr mangetons or mobilities
! such specification is stored in the property record, not here
!    integer prop_spec,listid
! >>> added "listid" as a conection to the "state variable" listing here.
! This replaces TC, BMAG, MQ etc included as "state variables" in order to
! list their values.  In this way all propids become available
 end TYPE gtp_propid
! the value TYPTY stored in property records is "idprop" or
! if IDELSUFFIX set then 100*"idprop"+ellista index of element
! if IDCONSUFFIX set then 100*"idprop"+constituent index
! When the parameter is read the suffix symbol is translated to the
! current element or constituent index
 TYPE(gtp_propid), dimension(:), private, allocatable :: propid
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_phase_add
! record for additions to the Gibbs energy for a phase like magnetism
! addrecno: ?
! aff: antiferomagnetic factor (Inden model)
! need_property: depend on these properties (like Curie T)
! explink: function to calculate with the properties it need
! nextadd: link to another addition
    integer type,addrecno,aff
    integer, dimension(:), allocatable :: need_property
    TYPE(tpfun_expression), dimension(:), pointer :: explink
    TYPE(gtp_phase_add), pointer :: nextadd
    type(gtp_elastic_modela), pointer :: elastica
 END TYPE gtp_phase_add
! allocated when needed and linked from phase record
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! addition record to calculate the elastic energy contribution
! declared as allocatable in gtp_phase_add
 TYPE gtp_elastic_modela
! lattice parameters (configuration) in 3 dimensions
    double precision, dimension(3,3) :: latticepar
! epsilon in Voigt notation
    double precision, dimension(6) :: epsa
! elastic constant matrix in Voigt notation
    double precision, dimension(6,6) :: cmat
! calculated elastic energy addition (with derivative to T and P?)
    double precision, dimension(6) :: eeadd
! maybe more
 end TYPE gtp_elastic_modela
!\end{verbatim}
!-----------------------------------------------------------------
! NOTE: if one wants to model bond energies beteween sites in a phase
! like in a 3 sublattice sigma one can enter parameters like G(sigma,A:B:*)
! which will mean the bond energy between an A atom in first sublattice and
! B in the second.  The parameter G(sigma,B:A:*) will be different.  Such
! parameters are added to the Gibbs energy even if there are G(sigma,A:B:C)
!-----------------------------------------------------------------
!\begin{verbatim}
! a smart way to have an array of pointers used in gtp_phase
 TYPE endmemrecarray 
    type(gtp_endmember), pointer :: p1
 end TYPE endmemrecarray
!
 TYPE gtp_phase
! this is the record for phase model data. It points to many other records.
! Phases are stored in order of creation in phlista(i) and can be found
! in alphabetical order through the array phases(i)
! For TC freaks: Treat index to phlista(i) as LOKPH in iws
! sublista is now removed and all data included in phlista
! sublattice and constituent data (they should be merged)
! The constitent link is the index to the splista(i), same function
! as LOKSP in iws.  Species in alphabetcal order is in species(i)
! One can allocate a dynamic array for the constituent list, done
! by subroutine create_constitlist.
! Note that the phase has a dynamic status word status2 in gtp_phase_varres
! which can be differnt in different parallell calculations.
! This status word has the FIX/ENT/SUS/DORM status bits for example
! name: phase name, note composition sets can have pre and suffixes
! model: free text
! phletter: G for gas, L for liquid
! alphaindex: the alphabetcal order of the phase (excluding gas and liquids)
    character name*24,models*72,phletter*1
    integer status1,alphaindex
! noofcs: number of composition sets, 
! nooffs: number of fraction sets (replaces partitioned phases in TC)
    integer noofcs,nooffs
! additions: link to addition record list
! ordered: link to endmember record list
! disordered: link to endmember list for disordered fractions (if any)
    TYPE(gtp_phase_add), pointer :: additions
    TYPE(gtp_endmember), pointer :: ordered,disordered
! To allow parallel processing of endmembers, store a pointer to each here
    integer noemr,ndemr
    TYPE(endmemrecarray), dimension(:), allocatable :: oendmemarr,dendmemarr
!-----------------------------------------------------------------
! this used to be sublista but is now incorporated in gtp_phase !!!
! static data, contains pointers to constituent record and sites
! noofsubl: number if sublattices
! cslink: is index to first composition set (deleted)
! linktocs: array with indices to phase_varres records (to replace clink)
! tnooffr: total number of fractions (constituents)
! nooffr: array with number of constituents in each sublattice
! sites: array with site rations (? dynamic for ionic liquid)
! constitlist: indices of species that are constituents (in all soblattices)
    integer noofsubl,tnooffr
    integer, dimension(9) :: linktocs
    integer, dimension(:), allocatable :: nooffr
    double precision, dimension(:), allocatable :: sites
    integer, dimension(:), allocatable :: constitlist
! allocated in init_gtp.
 END TYPE gtp_phase
! NOTE phase with index 0 is the reference phase for the elements
! The array sublista is now merged into phlista
! allocated in init_gtp
 TYPE(gtp_phase), private, allocatable :: phlista(:)
 INTEGER, private, allocatable :: PHASES(:)
!\end{verbatim}
!-----------------------------------------------------------------
!
!===================================================================
!
! below here are data structures for equilibrium description
!
!===================================================================
!
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_condition
! these records form a circular list linked from gtp_equilibrium_data records
! each record contains a condition to be used for calculation
! it is a state variable equation or a phase to be fixed
! The state variable is stored as an integer with indices
! NOTE: some state variables cannot be used as conditions: Q=18, DG=19, 25, 26
! There can be several terms in a condition (like x(liq,c)-x(fcc,c)=0)
! noofterms: number of terms in condition expression
! statev: the type of state variable (must be the same in all terms)
!           negative value of statev means phase index for fix phase
! active: zero if condition is active, nonzero for other cases
! unit: is 100 if value in percent, can also be used for temperature unit etc.
! nid: identification sequential number (in order of creation), redundant
! iref: part of the state variable (iref can be comp.set number)
! iunit: ? confused with unit?
! symlink: index of symbol for prescribed value (1) and uncertainity (2)
! condcoeff: there is a coefficient and set of indices for each term
! prescribed: the prescribed value
! NOTE: if there is a symlink value that is the prescribed value
! current: the current value (not used?)
! uncertainity: the uncertainity (for experiments)
    integer :: noofterms,statev,active,unit,nid,iref,iunit
!    TYPE(putfun_node), pointer :: symlink1,symlink2
! better to let condition symbol be index in svflista array
    integer symlink1,symlink2
    integer, dimension(:,:), allocatable :: indices
    double precision, dimension(:), allocatable :: condcoeff
    double precision prescribed, current, uncertainity
    TYPE(gtp_condition), pointer :: next, previous
 end TYPE gtp_condition
! declared inside the gtp_equilibrium_data record
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_state_variable
! this is to specify a formal or real argument to a function of state variables
! istv: state variable index
! indices: additional specifiers like phase, component, etc.
! iref: if a specified reference state (for chemical potentials)
! iunit: 100 for percent
    integer istv,indices(4),iref,iunit
 end TYPE gtp_state_variable
! declared inside gtp_function
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
  TYPE gtp_putfun_lista
! these are records for state variable functions.  The function itself
! is handelled by the putfun package.
! linkpnode: pointer to start node of putfun expression
! narg: number of symbols in the function
! nactarg: number of actual parameter specifications needed in call
!   (like @P, @C and @S
! status: can be used for various things
! status bit SVFVAL=0 means value evaluated only when called with mode=1
! eqnoval: used to specify the equilibrium the value should be taken from
!    (for handling what is called "variables" in TC)
! name: name of symbol
     integer narg,nactarg,status,eqnoval
     type(putfun_node), pointer :: linkpnode
     character name*16
! this array has identification of state variable (and other function) symbols 
     integer, dimension(:,:), pointer :: formal_arguments
  end TYPE gtp_putfun_lista
! this is the global array with state variable functions
  TYPE(gtp_putfun_lista), dimension(:), allocatable :: svflista
! NOTE the value of a function is stored locally in each equilibrium record
! in array svfunres.
! The number of entered state variable functions. Used to store a new one
  integer, private :: nsvfun
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_fraction_set
! info about disordred fractions for some phases like ordered fcc, sigma etc
! latd: the number of sublattices added to first disordred sublattice
! ndd: sublattices for this fraction set, 
! tnoofxfr: number of disordered fractions
! tnoofyfr: same for ordered fractions (=same as in phlista).
! varreslink: index of disordered phase_varres, 
! phdapointer: pointer to the same phase_varres record as varreslink
!    (Note that there is a bit set indicating that the sublattices should 
!    be taken from this record)
! totdis: 0 indicates no total disorder (sigma), 1=fcc, bcc or hcp
! id: parameter suffix, D for disordered
! dsites: number of sites in sublattices, disordred fractions stored in
!    another phase_varres record linked from phdapointer
! splink: pointers to species record for the constituents
! nooffr: the number of fractions in each sublattice
! y2x: the conversion from sublattice constituents to disordered and
! dxidyj: are the the coeff to multiply the y fractions to get the disordered
!        xfra(y2x(i))=xfra(y2x(i))+dxidyj(i)*yfra(i)
! disordered fractions stored in the phase_varres record with index varreslink
!    (also pointed to by phdapointer).  Maybe phdapointer is redundant??
! arrays originally declared as pointers now changed to allocatable
    integer latd,ndd,tnoofxfr,tnoofyfr,varreslink,totdis
    character*1 id
    double precision, dimension(:), allocatable :: dsites
    integer, dimension(:), allocatable :: nooffr
    integer, dimension(:), allocatable :: splink
    integer, dimension(:), allocatable :: y2x
    double precision, dimension(:), allocatable :: dxidyj
! in parallel processing the disordered phase_varres record is linked
! by this pointer,  used in parcalcg and calcg_internal
    TYPE(gtp_phase_varres), pointer :: phdapointer
 END TYPE gtp_fraction_set
! these records are declared in the phase_varres record as DISFRA for 
! each composition set and linked from the phase_varres record
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_phase_varres
! Data here must be different in equilibria representing different experiments
! or calculated in parallel or results saved from step or map.
! nextfree: In unused phase_varres record it is the index to next free record
!    The global integer csfree is the index of the first free record
! phlink: is index of phase record for this phase_varres record
! status2: has phase status bits like ENT/FIX/SUS/DORM
! constat: array with status word for each constituent, any can be suspended
! yfr: the site fraction array
! mmyfr: min/max fractions
! abnorm(1): amount moles of atoms for a formula unit of the composition set
! abnorm(2): amount in mass of the formula unit (both set by set_constitution)
! sites: site ratios (which can vary for ionic liquids)
! prefix and suffix are added to the name for composition sets 2 and higher
! disfra: a structure describing the disordered fraction set (if any)
    integer nextfree,phlink,status2
    double precision, dimension(2) :: abnorm
    character*4 prefix,suffix
! for ionic liquid derivatives of sites wrt fractions
!    double precision, dimension(:,:), pointer :: dsitesdy
!    double precision, dimension(:,:), pointer :: d2sitesdy2
! changed to allocatable
    integer, dimension(:), allocatable :: constat
    double precision, dimension(:), allocatable :: yfr
    real, dimension(:), allocatable :: mmyfr
    double precision, dimension(:), allocatable :: sites
! for ionic liquid derivatives of sites wrt fractions
    double precision, dimension(:,:), allocatable :: dsitesdy
    double precision, dimension(:,:), allocatable :: d2sitesdy2
! for extra fraction sets, better to go via phase record index above
! this TYPE(gtp_fraction_set) variable is a bit messy.  Declaring it in this
! way means the record is stored inside this record.
    type(gtp_fraction_set) :: disfra
! It seems difficult to get the phdapointer in disfra record to work
! ---
! arrays for storing calculated results for each phase (composition set)
! old: amount(1): is amount formula units of the composition set
! old: amount(2): is net charge of phase
! amfu: is amount formula units of the composition set (calculated result)
! netcharge: is net charge of phase
! dgm: driving force (calculated result)
! amcom: not used
! damount: set to last change of phase amount in equilibrium calculations
! qqsave: values of qq calculated in set_constitution
!    double precision amount(2),dgm,amcom,damount,qqsave(3)
    double precision amfu,netcharge,dgm,amcom,damount,qqsave(3)
! Other properties may be that: gval(*,2) is TC, (*,3) is BMAG, see listprop
! nprop: the number of different properties (set in allocate)
! ncc: total number of site fractions (redundant but used in some subroutines)
! listprop(1): is number of calculated properties
! listprop(2:listprop(1)): identifies the property stored in gval(1,ipy) etc
!   2=TC, 3=BMAG. Properties defined in the gtp_propid record
    integer nprop,ncc
    integer, dimension(:), allocatable :: listprop
! gval etc are for all composition dependent properties, gval(*,1) for G
! gval(*,1): is G, G.T, G.P, G.T.T, G.T.P and G.P.P
! dgval(1,j,1): is first derivatives of G wrt fractions j
! dgval(2,j,1): is second derivatives of G wrt fractions j and T
! dgval(3,j,1): is second derivatives of G wrt fractions j and P
! d2gval(ixsym(i,j),1): is second derivatives of G wrt fractions i and j
    double precision, dimension(:,:), allocatable :: gval
    double precision, dimension(:,:,:), allocatable :: dgval
    double precision, dimension(:,:), allocatable :: d2gval
! added for strain/stress, current values of lattice parameters
    double precision, dimension(3,3) :: curlat
! saved values from last equilibrium calculation
!    double precision, dimension(:), allocatable :: dsf
    double precision, dimension(:,:), allocatable :: cinvy
    double precision, dimension(:), allocatable :: cxmol
    double precision, dimension(:,:), allocatable :: cdxmol
 END TYPE gtp_phase_varres
! this record is created inside the gtp_equilibrium record
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
 TYPE gtp_equilibrium_data
! this contains all data specific to an equilibrium like conditions,
! status, constitution and calculated values of all phases etc
! Several equilibria may be calculated simultaneously in parallell threads
! so each equilibrium must be independent 
! NOTE: the error code must be local to each equilibria!!!!
! During step and map thses records with results are saved
! values of T and P, conditions etc.
! Values here are normally set by external conditions or calculated from model
! local list of components, phase_varres with amounts and constitution
! lists of element, species, phases and thermodynamic parameters are global
! tpval(1) is T, tpval(2) is P, rgas is R, rtn is R*T
! status: not used yet?
! gtperr: local errr code (not used yet)
! eqno: sequential number assigned when created
! next: ?
! eqname: name of equilibrium
! tpval: value of T and P
! rtn: value of R*T
    integer status,gtperr,eqno,next
    character eqname*24
    double precision tpval(2),rtn
! this array has the local results of the state variable functions
! svfunres: the values of state variable functions for this equilibrium
    double precision, dimension(:), allocatable :: svfunres
! the experiments are used in assessments and stored like conditions 
! lastcondition: link to condition list
! lastexperiment: link to experiment list
    TYPE(gtp_condition), pointer :: lastcondition,lastexperiment
! components and conversion matrix from components to elements
! complist: array with components
! compstoi: stoichiometric matrix of compoents relative to elements
! invcompstoi: inverted stoichiometric matrix
    TYPE(gtp_components), dimension(:), allocatable :: complist
    double precision, dimension(:,:), allocatable :: compstoi
    double precision, dimension(:,:), allocatable :: invcompstoi
! one record for each phase that can be calculated
! index in this array is the same as in sublat record
! phase_varres: here all data for the phase is stored
    TYPE(gtp_phase_varres), dimension(:), allocatable :: phase_varres
! index to the tpfun_parres array is the same as in the global array tpres 
! eq_tpres: here calculated values of TP functions are stored
    TYPE(tpfun_parres), dimension(:), pointer :: eq_tpres
! current values of chemical potentials stored in component record
! xconc: convergence criteria for constituent fractions and other things
    double precision xconv
! delta-G value for merging gridpoints in grid minimizer
! smaller value creates problem for test step3.BMM, MC and austenite merged
    double precision :: gmindif=-5.0D-2
! maxiter: maximum mnumber of iterations allowed
    integer maxiter
 END TYPE gtp_equilibrium_data
! The primary copy of this structures is declared globally as FIRSTEQ here
! Others may be created when needed for storing experimental data or
! for parallel processing. A global array of these are
! TYPE(gtp_equilibrium_data), dimension(:), allocatable, target :: eqlista
 TYPE(gtp_equilibrium_data), dimension(:), allocatable, target :: eqlista
 TYPE(gtp_equilibrium_data), pointer :: firsteq
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! for each permutation in the binary interaction tree of an endmember one must
! keep track of the permutation and the permutation limit.
! It is not possible to push the value on pystack as one must remember
! them when changing the endmember permutation
! integer, parameter :: permstacklimit=150
 TYPE gtp_parcalc
! This record contains temporary data that must be separate in different
! parallell processes when calculating G and derivatives for any phase.
! There is nothing here that need to be saved after the calculation is finished
! global variables used when calculating G and derivaties
! sublattice with interaction, interacting constituent, endmember constituents
! PRIVATE inside this structure not liked by some compilers....
! endcon must have maxsubl dimension as it is used for all phases
    integer :: intlat(maxinter),intcon(maxinter),endcon(maxsubl)
! interaction level and number of fraction variables
    integer :: intlevel,nofc
! explained above, used for FCC and BCC permutations
!    integer, dimension(permstacklimit) :: lastperm,permlimit
! interacting constituents (max 4) for composition dependent interaction
! iq(j) indicate interacting constituents
! for binary RK+Muggianu iq(3)=iq(4)=iq(5)=0
! for ternary Muggianu in same sublattice iq(4)=iq(5)=0
! for reciprocal composition dependent iq(5)=0
! for Toop, Kohler and simular iq(5) non-zero (not implemented)
    integer :: iq(5)
! fraction variables in endmember (why +2?) and interaction
    double precision :: yfrem(maxsubl+2),yfrint(maxinter)
! local copy of T, P and RT for this equilibrium
    double precision :: tpv(2),rgast
!    double precision :: ymin=1.0D-30
 end TYPE gtp_parcalc
! this record is declared locally in subroutine calcg_nocheck
!\end{verbatim}
!-------------------------------------------------------------------
!\begin{verbatim}
   TYPE gtp_pystack
! records created inside the subroutine push/pop_pystack
! data stored during calculations when entering an interaction record
! previous: link to previous record in stack
! ipermutsave: permutation must be saved
! intrecsave: link to interaction record
! pysave: saved value of product of all constituent fractions
! dpysave: saved value of product of all derivatives of constituent fractions
! d2pysave: saved value of product of all 2nd derivatives of constit fractions
      TYPE(gtp_pystack), pointer :: previous
      integer :: pmqsave
      TYPE(gtp_interaction), pointer :: intrecsave
      double precision :: pysave
      double precision, dimension(:), allocatable :: dpysave
      double precision, dimension(:), allocatable :: d2pysave
   end TYPE gtp_pystack
! declared inside the calcg_internal subroutine
!\end{verbatim}
!-----------------------------------------------------------------
!
! a global array to provide information about composition sets
! phcs(nph) is the composition set counter for phase nph
  integer, dimension(maxph) :: phcs
!
!===================================================================
!
! Below are private global variables like free lists etc.
!
!===================================================================

! Several arrays with lists have a free list: csfree,addrecs,eqfree,reffree
!
!\begin{verbatim}
! counters for elements, species and phases initiated to zero
 integer, private :: noofel=0,noofsp=0,noofph=0
! counters for property and interaction records, just for fun
 integer, private :: noofprop,noofint,noofem
! free lists in phase_varres records and addition records
 integer, private :: csfree,addrecs
! free list of references and equilibria
 integer, private :: reffree,eqfree
! maximum number of properties calculated for a phase
 integer, private :: maxcalcprop=20
! highest used phase_varres record (for saving on file)
 integer, private :: highcs
! Trace for debugging (not used)
 logical, private :: ttrace
! minimum constituent fraction
 double precision :: bmpymin
! number of defined property types like TC, BMAG etc
 integer, private :: ndefprop
!\end{verbatim}

CONTAINS

! 1-5: initialization, number of things, find things, get things, set things, 
include "pmod25A.F90"

! 6: calculate things
include "pmod25B.F90"

! 7: state variable manipulations
include "pmod25C.F90"

! 8-9: state variable functions, interactive things
include "pmod25D.F90"

! 10: list things
include "pmod25E.F90"

! 11: save and read
include "pmod25F.F90"

! 12: enter data
include "pmod25G.F90"

! 13-15: status for things, unfinished things, internal stuff
include "pmod25H.F90"

! 16: Additions (magnetic and others)
include "pmod25I.F90"

! 17-18: Grid minimizer and miscellaneous
include "pmod25J.F90"


END MODULE GENERAL_THERMODYNAMIC_PACKAGE

