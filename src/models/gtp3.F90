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
! equilibrium data structure which also contains conditions,
! calculated values of TP functions and other symbols.  To identify a
! phase+comp.set a phasetuple has been intoduced.  This contains two
! integers, th first is the phase number, the second the composition
! set number.  The second or higher composition set of a phase will
! have a tupel index higher that the number of phases.
!
! One equilibrium record is created in init_gtp and it is called
! FIRSTEQ which is a global variable.  There is also an array EQLISTA
! which should contain all allocated equilibrium records, FIRSTEQ is a
! pointer to the first element in this array.  More equilibrium
! records can be initiated by the enter_equilibrium subroutine.  This
! copies the relevant data from FIRSTEQ.  After a second equilibrium
! is created it is forbidden to enter elements, species and phases and
! create additional fraction sets, i.e. one must not change the data
! structure except to add/remove composition sets (but this should
! anyway be avoided).  Composition sets must be created in all
! equilibrium records at the same time (if done in a thread then all
! threads must stop while this is done).  During step/map calculation
! each calculated equilibria is saved for later use in plotting och
! other postprocessing.  These saved equilibria may have different
! number of composition sets so great care must be taken using them.
!
! The equilibrium data record is "stand alone" and contains all necessary
! data to describe the equilibrium (except the model parameters and other
! static data).  In parallel processing each thread will have its own
! equilibrium data record.
!
! The intention is that several equilibra can be created both to store
! individual experimental data in assessments and for each thread in
! parallel.  In the equilibrium record there are conditions,
! components (with chemical potentials) and an error code and most
! important, the phase_varres record array with one or more record for
! each phase.  This array must be identical in all equilibria recods.
! Each composition set has a phase_varres record and they are linked
! from the phase record by the LIKTOCS array.  As the phase_varres
! records are in an array the link is simply an integer index of this
! array.  There is a free list (in FIRSTEQ) in the phase_varres array
! to be used when adding or removing a composition set.  The EQ_TPRES
! array is declared inside the equilibrium record for calculated
! results of the TPFUNS as these can be different in each equilibria.
! The index to a function in EQ_TPRES is the same as the index to the
! TPFUN array declared globally in TPFUN.  The TPFUN array has the
! actual expression, and EQ_TPRES has the last calculated results,
! which can be different in each equilibrium.  The TPFUN index is used
! in property records to specify the function of a parameter.
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
!--------------------------------------------------------------------------
!
! EXTERNAL MODULES
! metlib package
  use metlib
!
! routines for inverting matrix, solving system of eqs, eigenvalues etc
!  use lukasnum  ! for Lukas solver
  use ocnum      ! for LAPACK and BLAS
!
! for parallel processing
!$  use OMP_LIB
!
!--------------------------------------------------------------------------
!
! Versions
! date       item
! 2013.03.01 Release version 1
! 2015.01.07 Release version 2
! 2016.02.14 Release version 3
! 2017.02.10 Release version 4
! 2018.03.02 Release version 5
! 2020.03.12 Release version 6
! after version on github numbered 4.011, incremented for each update
!
!=================================================================
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------
! error messages
! numbers 4000 to 4220 defined.  gx%bmperr is set to message index
! A lot of error flags set have no messages ....
  integer, parameter :: nooferm=4399
  character (len=64), dimension(4000:nooferm) :: bmperrmess
! The first 30 error messages mainly for TP functions
  data bmperrmess(4000:4199)&
      /'Too many coefficients in a TP function.                         ',&
       'Illegal character in a TP function, digit expected.             ',&
       'Unknown symbol in TP function                                   ',&
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
! These error mainly in GTP
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
       'Species symbol contain illegal letter or not letter A-Z as first',&
       'No elements or too many elements in species formula             ',&
       'Unknown element in species formula                              ',&
       'Negative stoichiometric factor in species                       ',&
       'The charge must be the final "element"                          ',&
       'Species already entered                                         ',&
! 4050
       'No such phase                                                   ',&
       'Unknown or ambiguous species name                               ',&
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
       'Reading unknown addition type from file                         ',&
! 4090
       'Addition already entered                                        ',&
       'No more addition records                                        ',&
       'Maximum 9 composition sets                                      ',&
       'Illegal composition set number                                  ',&
       'No more records for phases or composition sets.                 ',&
       'Hidden phase cannot be ENTERED, SUSPENDED, DORMANT or FIXED     ',&
       'Ambiguous or unknown constituent                                ',&
       'Too many argument to a state variable                           ',&
       'This state variable must have two arguments                     ',&
       'First character of a state variable is wrong                    ',&
! 4100
       'State variable starting with M not followed by U                ',&
       'State variable starting with L not followed by NAC              ',&
       'Missing ( for arguments of state variable                       ',&
       'Missing ) after arguments of state variable                     ',&
       'Unknown phase used as state variable argument                   ',&
       'Unknown constituent used as state variable argument             ',&
       'Unknown component used as state variable argument               ',&
       'State variable starting with D not followed by G                ',&
       'State variable starting with T follwed by other character than C',&
       'State variable starting with B missing P, MAG, M, V, W or F     ',&
! 4110
       'This state variable cannot not have two arguments               ',&
       'This state variable must have an argument                       ',&
       'Impossible reference state for this component                   ',&
       'No such property calculated for this phase                      ',&
       'Property normallized by volume impossible as no volume data     ',&
       'Property per formula unit is phase specific                     ',&
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
       'No such condition or experiment                                 ',&
       'Function name must start with a letter A-Z                      ',&
       'Function name and expression must be separated by "="           ',&
       'Error in function expression (putfun)                           ',&
       'Unknown symbol used in function                                 ',&
       'Symbol with this name already entered                           ',&
       'Symbol name must start with letter A-Z and not be reserved      ',&
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
       'Grid minimizer cannot be used with the current set of conditions',&
       'Too many gridpoints                                             ',&
       'No phases and no gridpoints for grid minimization               ',&
       'Grid minimizer wants but must not create composition sets       ',&
       'Non-existing fix phase                                          ',&
       'N, X, B or W cannot have two indices for use of grid minimizer  ',&
! 4180
       'Condition on B is not allowed for grid minimizer                ',&
       'An element has no composition in grid minimizer                 ',&
       'Too complicated mass balance conditions                         ',&
       'Two mass balance conditions for same element                    ',&
       'Cannot handle conditions on both N and B                        ',&
       'No mole fractions when summing composition                      ',&
       'Error in TDB file, missing function                             ',&
       'Temperature (K) or pressure (Pa) values must be larger than 0.1 ',&
       'No such state variable                                          ',&
       'Too many conditions on potentials                               ',&
! 4190
       'File already exist, overwriting not allowed                     ',&
       'Activity conditions must be larger than zero                    ',&
       'Cannot handle two fix phases                                    ',&
       'Too many stable phases                                          ',&
       'This phase must not be stable                                   ',&
       'Attempt to remove the only stable phase                         ',&
       'Enthalpy condition on unstable phase                            ',&
       'Illegal wildcard constituent in ionic liquid model              ',&
       'No equilibrium calculated, cannot calculate dot derivative      ',&
       'Error calculating equilibrium matrix for dot derivative         '/
! 4200 mainly errors in minimizer
  data bmperrmess(4200:4399)&
      /'No phase that can be set stable                                 ',&
       'Attempt to set too many phases as stable                        ',&
       'Total amount is negative                                        ',&
       'Error solving equilibrium matrix                                ',&
       'Too many iterations                                             ',&
       'Phase matrix singular                                           ',&
       'Cannot handle models without analytical second derivativatives  ',&
       'This type of condition not yet implemented                      ',&
       'This type of condition is not allowed                           ',&
       'Error setting up system matrix, too many equations              ',&
! 4210
       'Phase change not allowed                                        ',&
       'Attempt to delete composition sets when many equilibria         ',&
       'Too many equation in equilibrium matrix                         ',&
       'Derivatives with respect to T and P only are allowed            ',&
       'Error creating system matrix in initiate meqrec subroutine      ',&
       'This dot derivative not yet implemented                         ',&
       'Wildcard not allowed in dot derivative                          ',&
       'Use "calculate symbol" for state variable symbols               ',&
       'This experiment is not acivated                                 ',&
       'Too many equilibria in STEP/MAP, save on file not implemented   ',&
! mainly errors in STEP/MAP
! 4220 step/map
       'STEP/MAP error calculating node point, trying to decrease step  ',&
       'STEP/MAP error calculating node point, axis condition not found ',&
       'STEP/MAP error calculating node point, another phase stable     ',&
       'STEP/MAP error calculating node point, too many stable phases   ',&
       'Cannot find start equilibrium for step/map                      ',&
       'Startpoint for step/map outside axis limits                     ',&
       'Cannot yet handle nodepoints with more than 2 exits             ',&
       'Phase set changed in start point                                ',&
       'Only two axis implemented currently                             ',&
       'Axis direction error, no such axis                              ',&
! 4230
       'STEP/MAP tries to set the only stable phase as fix              ',&
       'Too many stable phases during mapping                           ',&
       'Another phase wants to be stable at node point                  ',&
       'No phase change searching along an axis for a start point       ',&
       'Internal error handling fix phases at node point                ',&
       'Too many phases set fix during mapping                          ',&
       'Mapping cannot handle expressions as conditions                 ',&
       'Node with no exit lines                                         ',&
       'Attempt to remove the only stable phase                         ',&
       'Yet another never never error                                   ',&
! 4240
       'Too many fix phases during mapping                              ',&
       'More than one entered phase                                     ',&
       'Not a single entered phase                                      ',&
       'Whops, mapping without conditions ...                           ',&
       'I give up on this line                                          ',&
       'Unknown problem                                                 ',&
       'Two phases compete to be stable                                 ',&
       'Nothing to plot in ocplot                                       ',&
       'No data so no plot                                              ',&
       'No experiments                                                  ',&
! more error messages for GTP and other modules
! 4250
       'Too many parameter identifiers, increase maxprop                ',&
       'Calling mass_of with illegal component number                   ',&
       'No such phase tuple index                                       ',&
       'Internal error, not a single lattice for a phase                ',&
       'Illegal phase index                                             ',&
       'The partially ionic liquid model must have two sublattices      ',&
       'This phase cannot be reference phase for this component         ',&
       'Internal error, constituent index outside range                 ',&
       'Same constituent twice in one sublattice                        ',&
       'Too many phases, increase dimension of phlista                  ',&
! 4260
       'The partially ionic liquid model has only cations in first subl.',&
       'Illegal parameter with wildcards mixed with cations             ',&
       'The partially ionic liquid model not only wildcard on 2nd subl. ',&
       'The partially ionic liquid model has no catioons on 2nd subl.   ',&
       'Only neutrals on 2nd sublattice of I2SL if wildcard on first    ',&
       'Illegal interaction parameter                                   ',&
       'Same constituent twice in interaction parameter                 ',&
       'There must be at least 4 sublattices for a phase with F/B option',&
       'Maximum two interaction levels using the F option               ',&
       'Internal error, unknown case for endmember permutation          ',&
! 4270
       'Interaction must be on first sublattice using option F or B     ',&
       'Cannot find endmember element for permutation                   ',&
       'Internal error, unknown case for permutations                   ',&
       'Internal error, too complicated                                 ',&
       'Internal error generating fcc permutations                      ',&
       'This excess parameter not yet implemented in option F or B      ',&
       'Internal error generating permutations for option F             ',&
       'BCC permutations (option B) not yet implemented                 ',&
       'Subcommand error when enter many_equilibria                     ',&
       'Too many columns when entering many_equilibria row              ',&
! 4280
       'Table row missing in colum when entering many_equilbria         ',&
       'Number expected after specifying fix phase                      ',&
       'Phase name expected after status command                        ',&
       'Too many equilibra, increase dimension of eqlista               ',&
       'Equilibrium name must start with a letter A-Z                   ',&
       'Cannot overwrite the default equilibrium                        ',&
       'Illegal use of wildcard                                         ',&
       'Error in constituent dependence for parameter idenifier         ',&
       'Yet another never never error                                   ',&
       'Charge must be given as /+ or /-                                ',&
! 4290
       'Error in parameter identifier                                   ',&
       'Phase missing in parameter                                      ',&
       'No such property name or index                                  ',&
       'Illegal to have a symbol as value of T or P                     ',&
       'Illegal to set a fix phase as experiment                        ',&
       'Calling locate_condition with illegal index                     ',&
       'Calling apply_condition with illegal option                     ',&
       'Species names must be surrounded by ( ) for set input_amounts   ',&
       'Illegal to enter property to a species that is an element       ',&
       'Saved file not same version as program                          ',&
! 4300
       'Data record format on save file not the same as in program      ',&
       'Bibliographic record too long on save file                      ',&
       'Error reading records for a phase from save file                ',&
       'Failed entering function from save file                         ',&
       'Too long line on save file                                      ',&
       'No element symbol after ELEMENT keyword in TDB file             ',&
       'No information after SPECIES keyword on TDB file                ',&
       'No terminator after FUNCTION keyword on TDB file                ',&
       'The CONSTITUENT keyword must follow directly after PHASE keyword',&
       'Error extracting constituents for a phase                       ',&
! 4310
       'Error that final : for constituents missing                     ',&
       'Empty line after FUNCTION keyword                               ',&
       'Line with PARAMETER keyword does not finish with !              ',&
       'Empty reference line on TDB file                                ',&
       'References must be surrounded by citation marks                 ',&
       'Function name must be on same line as FUNCTION keyword          ',&
       'End of file while searching for end of keyword in TDB file      ',&
       'Indices error in old state variable format                      ',&
       'Unknown state variable or property                              ',&
       'Character variable length insufficient for output of values     ',&
! 4320
       'State variable has illegal argument type                        ',&
       'Error calculating eigenvalues of phase matrix                   ',&
       'Only a single symbol allowed                                    ',&
       'Symbol must be a constant                                       ',&
       'Value of PHSTATE not correct                                    ',&
       'Illegal bit number for phase status                             ',&
       'Illegal phase for setting status bit                            ',&
       'Illegal selection of old phase status                           ',&
       'Condition specified by number must be followed by :=            ',&
       'Calling create_interaction with too many permutations           ',&
! 4330
       'No such addition type                                           ',&
       'Cp model not yet implemented                                    ',&
       'Magnetic model with separate Curie and Neel T not yet implement ',&
       'Addition model not yet implemented                              ',&
       'Not implemented this way                                        ',&
       'Parameter identifier not found                                  ',&
       'Value for model parameter identifier not found                  ',&
       'Flory-Huggins model must have one lattice and site              ',&
       'Too many parameter properties for this phase                    ',&
       'Internal error, listprop not allocated                          ',&
! 4340
       'Max two levels of interactions allowed                          ',&
       'Wildcard parameters not allowed in 2nd sublattice of I2SL model ',&
       'Illegal interaction parameter                                   ',&
       'Ternary cation interactions not yet implemented in I2SL         ',&
       'Too many phases for the global gridminimizer                    ',&
       'Global minimization with a fix phase not possible               ',&
       'Internal problems in grid minimizer                             ',&
       'Interaction levels more than 5 levels deep                      ',&
       'A TP function with this name already entered                    ',&
       'Illegal value of TP function index                              ',&
! 4350
       'A never never error evaluating a TP function                    ',&
       'Cannot find this TP function                                    ',&
       'Current equilibrium not global, gridmin found gridpoint below   ',&
       'Nodepoint not global, line ignored                              ',&
       'Illegal numerical value in equilibrium matrix                   ',&
       'Wrong version of data on unformatted file                       ',&
       'Error reserving space for unformatted save                      ',&
       'Error saving unformatted data file                              ',&
       'Recalculate as gridpoint below current equilibrium              ',&
       'Slow convergence with same set of stable phases                 ',&
! 4360
       'Too large change on axis, terminating mapping                   ',&
       'Model parameter value not calculated                            ',&
       'New set of components are not independent                       ',&
       'No equilibrium, a restored phase should be stable               ',&
       'Two phases with same composition stable at nodepoint            ',&
       'Gridtest indicate global minimization needed                    ',&
       'Gridtest request recalculation without gridminimizer            ',&
       'Missing property for calculating addition                       ',&
       'Tried halfstep 3 times, giving up on this line                  ',&
       'Repeated error calling map_calcnode, line terminated            ',&
! 4370
       'Error allocating data, no free memory                           ',&
       'Nonlinear equation solver HYBRD1 error                          ',&
       'Error from DGETRS/F generating isopleth invariant exits         ',&
       'Supressed value due to special circumstances                    ',&
       'Mobility parameters must not have wildcard constituents         ',&
       '5                                                               ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
! 4380
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '5                                                               ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
! 4390
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       '5                                                               ',&
       '                                                                ',&
       '                                                                ',&
       '                                                                ',&
       'No message assigned                                             '/
! last used error codes above
!
!=================================================================
!
! STATUS BITS are numbered 0-31
!\begin{verbatim}
!-Bits in GLOBAL status word (GS) in globaldata record
! level of user: beginner, occational, advanced; NOGLOB: no global gridmin calc
! NOMERGE: no merge of gridmin result, 
! NODATA: not any data, 
! NOPHASE: no phase in system, 
! NOACS: no automatic creation of composition set for any phase
! NOREMCS: do not remove any redundant unstable composition sets
! NOSAVE: data changed after last save command
! VERBOSE: maximum of listing
! SETVERB: permanent setting of verbose
! SILENT: as little output as possible
! NOAFTEREQ: no manipulations of results after equilibrium calculation
! XGRID: extra dense grid for all phases
! NOPAR: do not run in parallel
! NOSMGLOB do not test global equilibrium at node points
! NOTELCOMP the elements are not the components
! TGRID use grid minimizer to test if global after calculating equilibrium
! OGRID use old grid generator
! NORECALC do not recalculate equilibria even if global test after fails
! OLDMAP use old map algorithm
! NOAUTOSP do not generate automatic start points for mapping
! GSYGRID extra dense grid
! GSVIRTUAL (CCI) enables calculations with a virtual element
! >>>> some of these should be moved to the gtp_equilibrium_data record
  integer, parameter :: &
       GSBEG=0,       GSOCC=1,        GSADV=2,      GSNOGLOB=3,  &
       GSNOMERGE=4,   GSNODATA=5,     GSNOPHASE=6,  GSNOACS=7,   &
       GSNOREMCS=8,   GSNOSAVE=9,     GSVERBOSE=10, GSSETVERB=11,&
       GSSILENT=12,   GSNOAFTEREQ=13, GSXGRID=14,   GSNOPAR=15,  &
       GSNOSMGLOB=16, GSNOTELCOMP=17, GSTGRID=18,   GSOGRID=19,  &
       GSNORECALC=20, GSOLDMAP=21,    GSNOAUTOSP=22,GSYGRID=23,  &
       GSVIRTUAL=24
!----------------------------------------------------------------
!-Bits in ELEMENT record
  integer, parameter :: &
       ELSUS=0,       ELDEL=1
!----------------------------------------------------------------
!-Bits in SPECIES record
! SUS   Suspended,
! IMSUS implicitly suspended (when element suspended)
! EL    species is element, 
! VA    species is the vacancy
! ION   species have charge, 
! SYS   species is (system) component
! UQAC  species used in uniquac model (2 extra reals for area and volume)
  integer, parameter :: &
       SPSUS=0, SPIMSUS=1, SPEL=2, SPVA=3, &
       SPION=4, SPSYS=5,   SPUQC=6
!\end{verbatim}
!----------------------------------------------------------------
! Many not implemented
!\begin{verbatim}
!-Bits in PHASE record STATUS1 there are also bits in each phase_varres record!
! HID phase is hidden (not implemented)
! IMHID phase is implictly hidden (not implemented)
! ID phase is ideal, substitutional and no interaction
! NOCV phase has no concentration variation (I am not sure it is set or used)
! HASP phase has at least one parameter entered
! FORD phase has 4 sublattice FCC ordering with parameter permutations
! BORD phase has 4 sublattice BCC ordering with parameter permutations
! SORD phase has TCP type ordering (do not subract ordered as disordered, NEVER)
! MFS phase has a disordered fraction set
! GAS this is the gas phase (first in phase list) 
! LIQ phase is liquid (can be several but listed directly after gas)
! IONLIQ phase has ionic liquid model (I2SL)
! AQ1 phase has aqueous model (not implemented)
! STATE elemental liquid twostate (2-state) model parameter
! QCE phase has corrected quasichemical entropy (Hillerst-Selleby-Sundman)
! CVMCE phase has some CVM ordering entropy (not implemented)
! EXCB phase need explicit charge balance (has ions)
! XGRID use extra dense grid for this phase
! FACTCE phase has FACT quasichemical SRO model (not implemented)
! NOCS not allowed to create composition sets for this phase
! HELM parameters are for a Helmholz energy model (not implemented),
! PHNODGDY2 phase has model with no analytical 2nd derivatives
! not implemented ELMA phase has elastic model A (not implemented)
! EECLIQ this is the condensed phase (liquid) that should have highest entropy
! PHSUBO special use testing models DO NOT USE
! PALM interaction records numbered by PALMTREE NEEDED FOR PERMUTATIONS !!!
! MULTI may be used with care
! BMAV Xion magnetic model with average Bohr magneton number
! UNIQUAC The UNIQUAC fluid model
! TISR phase has the TSIR entropy model (E Kremer)
  integer, parameter :: &
       PHHID=0,     PHIMHID=1,    PHID=2,      PHNOCV=3, &     ! 1 2 4 8 : 0/F
       PHHASP=4,    PHFORD=5,     PHBORD=6,    PHSORD=7, &     ! 
       PHMFS=8,     PHGAS=9,      PHLIQ=10,    PHIONLIQ=11, &  ! 
       PHAQ1=12,    PH2STATE=13,  PHQCE=14,    PHCVMCE=15,&    ! 
       PHEXCB=16,   PHXGRID=17,   PHFACTCE=18, PHNOCS=19,&     !
       PHHELM=20,   PHNODGDY2=21, PHEECLIQ=22, PHSUBO=23,&     ! 
       PHPALM=24,   PHMULTI=25,   PHBMAV=26,   PHUNIQUAC=27, & !
       PHTISR=28                                          !
!
!----------------------------------------------------------------
!-Bits in PHASE_VARRES (constituent fraction) record STATUS2
! CSDFS is set if record is for disordred fraction set, then one must use
!     sublattices from fraction_set record
! CSDLNK: a disordred fraction set in this phase_varres record
! CSDUM2 and CSDUM3 not used
! CSCONSUS set if one or more constituents suspended (status array constat
!     specify constituent status)
! CSORDER: set if fractions are ordered (only used for BCC/FCC ordering
!     with a disordered fraction set).
! CSABLE: set if phase is stable after an equilibrium calculation ?? needed
! CSAUTO set if composition set created during calculations
! CSDEFCON set if there is a default constitution
! CSTEMPAR set if created by grid minimizer and can be suspended afterwards
!       when running parallel
! CSDEL set if record is not used but has been and then deleted (by gridmin)
! CSADDG means there are terms to be added to G 
! CSTEMPDOR means this compset was temporarily set dormant at an 
!       equilibrium calculation
   integer, parameter :: &
        CSDFS=0,    CSDLNK=1,  CSDUM2=2,    CSDUM3=3, &
        CSCONSUS=4, CSORDER=5, CSABLE=6,    CSAUTO=7, &
        CSDEFCON=8, CSTEMPAR=9,CSDEL=10,    CSADDG=11,&
        CSTEMPDOR=12
!\end{verbatim}
!----------------------------------------------------------------
!\begin{verbatim}
!-Bits in CONSTAT array for each constituent
! For each constituent: 
! SUS constituent is suspended (not implemented)
! IMSUS is implicitly suspended, 
! VA is vacancy
! QCBOND the constituent is a binary quasichemical cluster
   integer, parameter :: &
        CONSUS=0,   CONIMSUS=1,  CONVA=2,    CONQCBOND=3
!----------------------------------------------------------------
!-Bits in STATE VARIABLE FUNCTIONS (svflista)
! SVFVAL V symbol evaluated only when explicitly referenced (mode=1 in call)
! SVFEXT X symbol value taken from equilibrium %eqnoval
! SVCONST C symbol is a constant (can be changed with AMEND)
! SVFTPF - bit not used, replaced by export/import
! SVFDOT D symbol is a DOT function, like cp=h.t (also SVFVAL bit)
! SVFNOAM N symbol cannot be amended (only R, RT and T_C)
! SVEXPORT E symbol value exported to assessment coeff (TP constant)
! SVIMPORT I symbol value imported from TP-function (incl assessment coeff)
! ONLY ONE BIT CAN BE SET except for D and C+I and C+E,
! OTHER COMBINATIONS ARE NOT ALLOWED!!
!
   integer, parameter :: &
        SVFVAL=0,     SVFEXT=1,     SVCONST=2,     SVFTPF=3,&
        SVFDOT=4,     SVNOAM=5,     SVEXPORT=6,    SVIMPORT=7
!----------------------------------------------------------------
!-Bits in CEQ record (gtp_equilibrium_data)
! EQNOTHREAD set if equilibrium must be calculated before threading 
! (in assessment) for example if a symbol must be evaluated in this 
! equilibrium before used in another like H(T)-H298
! EQNOGLOB set if no global minimization
! EQNOEQCAL set if no successful equilibrium calculation made
! EQINCON set if current conditions inconsistent with last calculation
! EQFAIL set if last calculation failed
! EQNOACS set if no automatic composition sets ?? not used !! see GSNOACS
! EQGRIDTEST set if grid minimizer should be used after equilibrium
! EQGRIDCAL set if last calculation was using only gridminimizer
! EQMIXED set if mixed reference state for the elements
   integer, parameter :: &
        EQNOTHREAD=0, EQNOGLOB=1, EQNOEQCAL=2,  EQINCON=3, &
        EQFAIL=4,     EQNOACS=5,  EQGRIDTEST=6, EQGRIDCAL=7, &
        EQMIXED=8
!----------------------------------------------------------------
!-Bits in parameter property type record (gtp_propid)
! no T or P dependence (constant)
! only P dependence
! only T dependence
! there is an element suffix (like mobility),
! there is a constituent suffix
! Property has no addition (used when entering and listing data)
   integer, parameter :: &
        IDNOTP=0, IDONLYP=1, IDONLYT=2, IDELSUFFIX=3, IDCONSUFFIX=4,&
        IDNOADD=5
!----------------------------------------------------------------
!- Bits in condition status word (some set in onther ways??)
! singlevar means T=, x(el)= etc, singlevalue means value is a number
! phase means the condition is a fix phase
  integer, parameter :: &
       ACTIVE=0, SINGLEVAR=1, SINGLEVALUE=2, PHASE=3
!----------------------------------------------------------------
!- Bits in assessment head record status
! ahcoef set means coefficients are entered
  integer, parameter :: &
       AHCOEF=0
!
!----------------------------------------------------------------
!- Bits in addition record status word gtp_phase_add
! havepar set if the phase has parameters for this addition
! if not set the addition is not listed
! permol set if addition should be muliplied with number of atoms
  integer, parameter :: &
       ADDHAVEPAR=0, ADDPERMOL=1,ADDBCCMAG=2
!
! >>> Bits for symbols and TP functions missing ???
!\end{verbatim}
!
!----------------------------------------------------------------------
!
! Defining the phase status is very important, maybe a status for MAPFIX
! should be added.  Added EECDORM for solids with higher entropy than liquid
!\begin{verbatim}
! some constants, phase status
  integer, parameter :: EECDORM=-5
  integer, parameter :: PHHIDDEN=-4
  integer, parameter :: PHSUS=-3
  integer, parameter :: PHDORM=-2
  integer, parameter :: PHENTUNST=-1
  integer, parameter :: PHENTERED=0
  integer, parameter :: PHENTSTAB=1
  integer, parameter :: PHFIXED=2
  character (len=12), dimension(-5:2), parameter :: phstate=&
       (/'EEC_DORMANT ','HIDDEN      ','SUSPENDED   ','DORMANT     ',&
         'ENTERED UNST','ENTERED     ','ENTERED STBL','FIXED       '/)
!\end{verbatim}
!
!----------------------------------------------------------------------
!
!\begin{verbatim}
! Parameters defining the size of arrays etc.
! max elements, species, phases, sublattices, constituents (ideal phase)
!  integer, parameter :: maxel=100,maxsp=1000,maxph=600,maxsubl=10,maxconst=1000
! NOTE increasing maxph to 600 and maxtpf to 80*maxph made the equilibrium
! record very big and created problems storing equilibria at STEP/MAP!!!
  integer, parameter :: maxel=100,maxsp=1000,maxph=600,maxsubl=10,maxconst=1000
! maximum number of constituents in non-ideal phase
  integer, parameter :: maxcons2=300
! maximum number of elements in a species
  integer, parameter :: maxspel=10
! maximum number of references
  integer, private, parameter :: maxrefs=1000
! maximum number of equilibria
  integer, private, parameter :: maxeq=900
! some dp values, default precision of Y and default minimum value of Y
! zero and one set in tpfun
  double precision, private, parameter :: YPRECD=1.0D-6,YMIND=1.0D-30
! dimension for push/pop in calcg, max composition dependent interaction
  integer, private, parameter :: maxpp=1000,maxinter=3
! max number of TP symbols, TOO BIG VALUE MAKES SAVE AT STEP/MAP DIFFICULT
  integer, private, parameter :: maxtpf=20*maxph
!  integer, private, parameter :: maxtpf=80*maxph
! max number of properties (G, TC, BMAG MQ%(...) etc)
  integer, private, parameter :: maxprop=50
! max number of state variable functions
  integer, private, parameter :: maxsvfun=500
! version number of GTP (not OC)
  character*8, parameter :: gtpversion='GTP-3.30'
! THIS MUST BE CHANGED WHENEVER THE UNFORMATTED FILE FORMAT CHANGES!!!
  character*8, parameter :: savefile='OCF-3.20'
!\end{verbatim}
!=================================================================
!\begin{verbatim}
! The number of additions to the Gibbs energy of a phase is increasing
! This is a way to try to organize them.  Each addtion has a unique
! number identifying it when created, listed or calculated.  These
! numbers are defined here
  integer, public, parameter :: INDENMAGNETIC=1
  integer, public, parameter :: XIONGMAGNETIC=2
  integer, public, parameter :: DEBYECP=3
  integer, public, parameter :: EINSTEINCP=4
  integer, public, parameter :: TWOSTATEMODEL1=5
  integer, public, parameter :: ELASTICMODEL1=6
  integer, public, parameter :: VOLMOD1=7
  integer, public, parameter :: UNUSED_CRYSTALBREAKDOWNMOD=8
  integer, public, parameter :: SECONDEINSTEIN=9
  integer, public, parameter :: SCHOTTKYANOMALY=10
  integer, public, parameter :: DIFFCOEFS=11
! with composition independent G2 parameter
  integer, public, parameter :: TWOSTATEMODEL2=12
! name of additions:
  character(len=24) , public, dimension(12), parameter :: additioname=&
       ['Inden-Hillert magn model','Inden-Xiong magn model  ',&
       'Debye CP model          ','Einstein Cp model       ',&
       'Liquid 2-state model    ','Elastic model A         ',&
       'Volume model A          ','Unused CBT model        ',&
       'Smooth CP step          ','Schottky Anomaly        ',&
       'Diffusion coefficients  ','                        ']
!       123456789.123456789.1234   123456789.123456789.1234
! Note that additions often use extra parameters like Curie or Debye
! temperatures defined by model parameter identifiers stored in gtp_propid
!\end{verbatim}
!\begin{verbatim}  
! Model parameter identifiers entered in gtp3A.F90 and used mainly in gtp3H
! to calculate additions.  Used also in gtp3B when entering parameters
! Index Name Used in addition/model
!  1     G    Gibbs energy for endmembers or interactions
!  2     TC   Curie T in Inden-Hillert magnetic model
!  3 BMAG  - -                                   1 Average Bohr magneton numb
!  4 CTA   - P                                   2 Curie temperature
!  5 NTA   - P                                   2 Neel temperature
!  6 IBM   - P &<constituent#sublattice>;       12 Individual Bohr magneton num
!  7 THET  - P                                   2 Debye or Einstein temp
!  8 V0    - -                                   1 Volume at T0, P0
!  9 VA    T -                                   4 Thermal expansion
! 10 VB    T P                                   0 Bulk modulus
! 11 VC    T P                                   0 Alternative volume parameter
! 12 VS    T P                                   0 Diffusion volume parameter
! 13 MQ    T P &<constituent#sublattice>;       10 Mobility activation energy
! 14 MF    T P &<constituent#sublattice>;       10 RT*ln(mobility freq.fact.)
! 15 MG    T P &<constituent#sublattice>;       10 Magnetic mobility factor
! 16 G2    T P                                   0 Liquid two state parameter
! 17 THT2  - P                                   2 Smooth step function Tcrit
! 18 DCP2  - P                                   2 Smooth step function increm.
! 19 LPX   T P                                   0 Lattice param X axis
! 20 LPY   T P                                   0 Lattice param Y axis
! 21 LPZ   T P                                   0 Lattice param Z axis
! 22 LPTH  T P                                   0 Lattice angle TH
! 23 EC11  T P                                   0 Elastic const C11
! 24 EC12  T P                                   0 Elastic const C12
! 25 EC44  T P                                   0 Elastic const C44
! 26 UQT   T P &<constituent#sublattice>;       10 UNIQUAC residual parameter
! 27 RHO   T P                                   0 Electric resistivity
! 28 VISC  T P                                   0 Viscosity
! 29 LAMB  T P                                   0 Thermal conductivity
! 30 HMVA  T P                                   0 Enthalpy of vacancy form.
! 31 TSCH  - P                                   2 Schottky anomaly T
! 32 CSCH  - P                                   2 Schottky anomaly Cp/R.
! 33 QCZ   - -                                   1 MQMQA cluster coord factor
! DO NOT CHANGE THE ORDER in gtp3A, that would require changes elsewhere too
!\end{verbatim}  
! =================================================================
!
! below here are data structures and global data in this module
!
! First those belonging to the TPFUN package
!
!=================================================================
! VARIABLES and STRUCTURES originally in TPFUN
! length of a function symbol
  integer, parameter :: lenfnsym=16
  integer, private :: freetpfun
  double precision, parameter :: zero=0.0D0,one=1.0D0
!
!\begin{verbatim}
  TYPE gtp_parerr
! This record contains the global error code.  In parallel processing each
! parallel processes has its own error code copied to this if nonzero
! it should be replaced by gtperr for separate errors in treads
     INTEGER :: bmperr
  END TYPE gtp_parerr
  TYPE(gtp_parerr) :: gx
! needed to have error code as private in threads
!$OMP  threadprivate(gx)
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
  integer, parameter :: tpfun_expression_version=1
  TYPE tpfun_expression
! Coefficients, T and P powers, unary functions and links to other functions
     integer noofcoeffs,nextfrex
     double precision, dimension(:), pointer :: coeffs
! each coefficient kan have powers of T and P/V and links to other TPFUNS
! and be multiplied with a following LOG or EXP term. 
! wpow USED FOR MULTIPLYING WITH ANOTHER FUNCTION!!
     integer, dimension(:), pointer :: tpow
     integer, dimension(:), pointer :: ppow
     integer, dimension(:), pointer :: wpow
     integer, dimension(:), pointer :: plevel
     integer, dimension(:), pointer :: link
  END TYPE tpfun_expression
! These records are allocated when needed, not stored in arrays
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! BITS in TPFUN
! TPCONST     set if a constant value
! TPOPTCON    set if optimizing value
! TPNOTENT    set if referenced but not entered (when reading TDB files)
! TPVALUE     set if evaluated only explicitly (keeping its value)
! TPEXPORT    set if value should be exported to symbol
! TPIMPORT    set if value should be imported from symbol (only for constants)
! TPINTEIN    set if value should always be calculated
  integer, parameter :: &
       TPCONST=0,    TPOPTCON=1,   TPNOTENT=2,    TPVALUE=3, &
       TPEXPORT=4,   TPIMPORT=5
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
  integer, parameter :: tpfun_root_version=1
  TYPE tpfun_root
! Root of a TP function including name with links to coefficients and codes
! and results.  Note that during calculations which can be parallelized
! the results can be different for each parallel process
     character*(lenfnsym) symbol
! Why are limits declared as pointers?? They cannot be properly deallocated
! limits are the low temperature limit for each range
! funlinks links to expression records for each range
! each range can have its own function, status indicate if T and P or T and V
! nextorsymbol is initiated to next index, then possible symbol link!
! forcenewcalc force new calculation when optimizing variable changed
! rewind is used to check for duplicates reading from TDB file
! not saved on unformatted files
! If bit TPIMPORT set the function must be a constant
!    and nextorsymbol is index of symbol
! If bit TPEXPORT set then the value of the function (not the derivatives)
!    and nextorsymbol is index of symbol
!     integer noofranges,nextfree,status,forcenewcalc
     integer noofranges,nextorsymbol,status,forcenewcalc,rewind
     double precision, dimension(:), pointer :: limits
     TYPE(tpfun_expression), dimension(:), pointer :: funlinks
     double precision hightlimit
  END TYPE tpfun_root
! These records are stored in arrays as the actual function is global but each
! equilibrium has its own result array (tpfun_parres) depending on the local
! values of T and P/V.  The same indiex is used in the global and local arrays.
! allocated in init_gtp
  TYPE(tpfun_root), private, dimension(:), pointer :: tpfuns
!\end{verbatim}

!-----------------------------------------------------------------
!\begin{verbatim}
  integer, parameter :: tpfun_parres_version=1
  TYPE tpfun_parres
! Contains TP results, 6 double for results and 2 doubles for T and P 
! values used to calculate the results
! Note that during calculations which can be parallelized the
! results can be different for each tread
     integer forcenewcalc
     double precision, dimension(2) :: tpused
     double precision, dimension(6) :: results
  END TYPE tpfun_parres
! This array is local to the gtp_equilibrium_data record
! index is the same as the function
!\end{verbatim}
!
! =============================== end of TPFUN data structures
!
!\begin{verbatim}
  TYPE gtp_global_data
! status should contain bits how advanced the user is and other defaults
! it also contain bits if new data can be entered (if more than one equilib)
! sysparam are variables for different things
! sysparam(1) unused
! sysparam(2) number of equilibria between each check of spinodal at STEP/MAP??
! sysparem(3-10) unused ...
! sysreal(1) is the minimum T for EET check (equi-entopy T, Hickel)
!            if zero no EET c
     integer status
     character name*24
     double precision rgas,rgasuser,pnorm
! these are explicitly set to zero in new_gtp
     double precision, dimension(10) :: sysreal=zero
     integer :: sysparam(10)=0
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
! this constant must be incremented whenever a change is made in gtp_element
  INTEGER, parameter :: gtp_element_version=1
  TYPE gtp_element
! data for each element: symbol, name, reference state, mass, h298-h0, s298
     character :: symbol*2,name*12,ref_state*24
     double precision :: mass,h298_h0,s298
! splink: index of corresponding species in array splink
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
! this constant must be incremented whenever a change is made in gtp_species
  INTEGER, parameter :: gtp_species_version=2
  TYPE gtp_species
! data for each species: symbol, mass, charge, extra, status
! mass is in principle redundant as calculated from element mass
     character :: symbol*24
     double precision :: mass,charge
! alphaindex: the alphabetical order of this species
! noofel: number of elements
! nextra: number of extra properties (size of spextra)
     integer :: noofel,status,alphaindex
! Use an integer array ellinks to indicate the elements in the species
! The corresponing stoichiometry is in the array stochiometry
     integer, dimension(:), allocatable :: ellinks
     double precision, dimension(:), allocatable :: stoichiometry
! Can be used for extra species properties as in UNIQUAC models (area, volume)
     double precision, dimension(:), allocatable :: spextra
  END TYPE gtp_species
! allocated in init_gtp
  TYPE(gtp_species), private, allocatable :: splista(:)
  INTEGER, private, allocatable :: SPECIES(:)
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented whenever a change is made in gtp_component
  INTEGER, parameter :: gtp_component_version=1
  TYPE gtp_components
! The components are simply an array of indices to species records
! the components must be "orthogonal".  There is always a set of "systems
! components" that by default is the elements.
! Later one may implement that the user can define a different "user set"
! and maybe also specific sets for each phase.
! The reference state is set as a phase and value of T and P.
! The name of the phase and its link and the link to the constituent is stored
! the endmember array is for the reference phase to calculate GREF
! The last calculated values of the chemical potentials (for user defined
! and default reference states) should be stored here.
! molat is the number of moles of components in the defined reference state
     integer :: splink,phlink,status
     character*16 :: refstate
     integer, dimension(:), allocatable :: endmember
     double precision, dimension(2) :: tpref
     double precision, dimension(2) :: chempot
     double precision mass,molat
  END TYPE gtp_components
! allocated in gtp_equilibrium_data
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented whenever a change is made in gtp_endmember
  INTEGER, parameter :: gtp_endmember_version=1
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
! antalem: sequenial order of creation, useful to keep track of structure ??
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
! this constant must be incremented when a change is made in gtp_interaction
  INTEGER, parameter :: gtp_interaction_version=1
  TYPE gtp_interaction
! this record constitutes the parameter tree. There are links to NEXT
! interaction on the same level (i.e. replace current fraction) and
! to HIGHER interactions (i.e. includes current fraction)
! There can be several permutations of the interactions (both sublattice
! and fraction permuted, like interaction in B2 (Al:Al,Fe) and (Al,Fe:Al))
! The number of permutations of interactions can be the same, more or fewer
! comparared to the lower order parameter (endmember or other interaction).
! The necessary information is stored in noofip.  It is not easy to keep
! track of permutations during calculations, the smart way to store the last
! permutation calculated is in this record ... but that will not work for
! parallel calculations as this record is static ...
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
     TYPE(gtp_tooprec), pointer :: tooprec
     integer, dimension(:), allocatable :: sublattice,fraclink,noofip
  END TYPE gtp_interaction
! allocated dynamically and linked from endmember records and other
! interaction records (in a binary tree)
!\end{verbatim}
!--------------------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_phasetuple
  INTEGER, parameter :: gtp_tooprec_version=1
  TYPE gtp_tooprec
! this indicates that an binary interaction parameter has Kohler/Toop model
! and which constituents involved.  A binary interaction can have a list of
! several gtp_tooprec records for each ternary system this record is involed.
! The binary fractions multiplied with the parameters will be calculated
! using this list of ternaries indicated by this list, see the documentation.
! const1, const2 and const3 are the constituents in alphabetical order
! (which is also the numerical order). 
! toop is 0 if Kohler extrapolation, if 1, 2 or 3 it indicates the Toop element 
! extra is used when listing data
! uniqid is a unique identification of the record, used for debugging
     integer toop,const1,const2,const3,extra,uniqid
! Each gtp_tooprec is part of 3 lists for binary interactions
! indicated by the 3 constituents in the gtp_tooprec record.
! next12 is the next link for the 2 (alphaetically) first constituents,
! next13 is the next link for first and third constituents
! next23 is the next link for the second and third constituents
! seq is a sequential link in the order the records created (try to fix bug)
     type(gtp_tooprec), pointer :: next12,next13,next23,seq
! this is very useful to obtain information in the calc_toop subroutine
     type(gtp_phase_varres), pointer :: phres
  end type gtp_tooprec
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_property
  INTEGER, parameter :: gtp_property_version=2
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
! this added to avoid problems if model param id has changed between saving
! and reading an unformatted file
     character*4 modelparamid
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
! this constant must be incremented when a change is made in gtp_biblioref
! old name gtp_datareference
  INTEGER, parameter :: gtp_biblioref_version=1
  TYPE gtp_biblioref
! store data references
! reference: can be used for search of reference
! refspec: free text
     character*16 reference
!     character*64, dimension(:), allocatable :: refspec
! this is Fortran 2003/2008 standard, not available in GNU 4.8
!     character(len=:), allocatable :: nyrefspec
! Use wpack routines!!!
     integer, dimension(:), allocatable :: wprefspec
  END TYPE gtp_biblioref
! allocated in init_gtp
  TYPE(gtp_biblioref), private, allocatable :: bibrefs(:)
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_propid
  INTEGER, parameter :: gtp_propid_version=1
  TYPE gtp_propid
! this identifies different properties that can depend on composition
! Property 1 is the Gibbs energy and the others are usually used in
! some function to contribute to the Gibbs energy like TC or BMAGN
! But one can also have properties used for other things like mobilities
! with additional especification like MQ&FE
! symbol: property identifier like G for Gibbs energy
! note: short description for listings
! prop_elsymb: additional for element dependent properties like mobilities
     character symbol*4,note*28,prop_elsymb*2
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
! These are the properties defined 2020-11-27/BoS defined in init_gtp
!   1 G     T P                                   0 Energy
!   2 TC    - P                                   2 Combined Curie/Neel T
!   3 BMAG  - -                                   1 Average Bohr magneton numb
!   4 CTA   - P                                   2 Curie temperature
!   5 NTA   - P                                   2 Neel temperature
!   6 IBM   - P &<constituent#sublattice>;       12 Individual Bohr magneton num
!   7 THET  - P                                   2 Debye or Einstein temp
!   8 V0    - -                                   1 Volume at T0, P0
!   9 VA    T -                                   4 Thermal expansion
!  10 VB    T P                                   0 Bulk modulus
!  11 VC    T P                                   0 Alternative volume parameter
!  12 VS    T P                                   0 Diffusion volume parameter
!  13 MQ    T P &<constituent#sublattice>;       10 Mobility activation energy
!  14 MF    T P &<constituent#sublattice>;       10 RT*ln(mobility freq.fact.)
!  15 MG    T P &<constituent#sublattice>;       10 Magnetic mobility factor
!  16 G2    T P                                   0 Liquid two state parameter
!  17 THT2  - P                                   2 Smooth step function T
!  18 DCP2  - P                                   2 Smooth step function value
!  19 LPX   T P                                   0 Lattice param X axis
!  20 LPY   T P                                   0 Lattice param Y axis
!  21 LPZ   T P                                   0 Lattice param Z axis
!  22 LPTH  T P                                   0 Lattice angle TH
!  23 EC11  T P                                   0 Elastic const C11
!  24 EC12  T P                                   0 Elastic const C12
!  25 EC44  T P                                   0 Elastic const C44
!  26 UQT   T P &<constituent#sublattice>;       10 UNIQUAC residual parameter
!  27 RHO   T P                                   0 Electric resistivity
!  28 VISC  T P                                   0 Viscosity
!  29 LAMB  T P                                   0 Thermal conductivity
!  30 HMVA  T P                                   0 Enthalpy of vacancy form.
!  31 TSCH  - P                                   2 Schottky anomaly T
!  32 CSCH  - P                                   2 Schottky anomaly Cp/R.
!  33 QCM   - -                                   1 Modif Quasichem model ratio
!  
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_phase_add
  INTEGER, parameter :: gtp_phase_add_version=2
  TYPE gtp_phase_add
! record for additions to the Gibbs energy for a phase like magnetism
! addrecno: ?
! aff: antiferomagnetic factor (Inden model)
! constants: for some constants needed ?? NEW
! status: BIT 0 set if there are parameters
!         BIT 1 set if magnetic model is for BCC
! need_property: depend on these properties (like Curie T)
! explink: function to calculate with the properties it need (not allocatable?)
! nextadd: link to another addition
     integer type,addrecno,aff,status
     integer, dimension(:), allocatable :: need_property
     double precision, dimension(:), allocatable :: constants
     TYPE(tpfun_expression), dimension(:), pointer :: explink
! The following declaration is illegal ... but above OK and I can allocate
!     TYPE(tpfun_expression), dimension(:), allocatable, pointer :: explink
     TYPE(gtp_phase_add), pointer :: nextadd
     type(gtp_elastic_modela), pointer :: elastica
     type(gtp_diffusion_model), pointer :: diffcoefs
! calculated contribution to G, G.T, G.P, G.T.T, G.T.P and G.P.P
     double precision, dimension(6) :: propval
  END TYPE gtp_phase_add
! allocated when needed and linked from phase record
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! addition record to calculate the elastic energy contribution
! declared as allocatable in gtp_phase_add
! this constant must be incremented when a change is made in gtp_elastic_modela
  INTEGER, parameter :: gtp_elastic_modela_version=1
  TYPE gtp_elastic_modela
! lattice parameters (configuration) in 3 dimensions
     double precision, dimension(3,3) :: latticepar
! epsilon in Voigt notation
     double precision, dimension(6) :: epsa
! elastic constant matrix in Voigt notation
     double precision, dimension(6,6) :: cmat
! calculated elastic energy addition (with derivative to T and P?)
     double precision, dimension(6) :: eeadd
! maybe more ...
  end TYPE gtp_elastic_modela
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! addition record to calculate diffusion coefficients
! declared as allocatable in gtp_phase_add
! this constant must be incremented when a change is made in gtp_elastic_modela
  INTEGER, parameter :: gtp_diffusion_model_version=1
  TYPE gtp_diffusion_model
! status bit 0 set means no calculation of this record
! dilute, simple or magnetic
     integer difftypemodel,status
!  alpha values for magnetic diffusion (for interstitials in constituent order)
     double precision, allocatable, dimension(:) :: alpha
! indices of dependent constituent in each sublattices
     integer, allocatable, dimension(:) :: depcon
! indices of constituents with zerovolume
     integer, allocatable, dimension(:) :: zvcon
! calculated diffusion matrix
     double precision, allocatable, dimension(:,:) :: dcoef
! Maybe we need one for each composition set?? at least to save the matrix
     type(gtp_diffusion_model), pointer :: nextcompset
! maybe more ...
  end TYPE gtp_diffusion_model
!\end{verbatim}
!------------------------------------------------------------------
!\begin{verbatim}
  TYPE gtp_tpfun_as_coeff
! this is a TPFUN converted to coefficents without any references to other
! functions.  Each function can have several T ranges and coefficents for T**n
! USED FOR SOLGASMIX
     double precision, dimension(:), allocatable :: tbreaks
     double precision, dimension(:,:), allocatable :: coefs
     integer, dimension(:,:), allocatable :: tpows
! this is used only during conversion
!     type(gtp_tpfun_as_coeff), pointer :: nextcrec
  end type gtp_tpfun_as_coeff
!
!--------------------------------------------------------------------------
  INTEGER, parameter :: gtp_tpfun2dat_version=1
  TYPE gtp_tpfun2dat
! this is a temporary storage of TP functions converted to arrays of
! coefficients.  Allocated as an array when necessary and the index in
! this array is the same index as for the TPfun
! USED FOR SOLGASMIX
     integer nranges
!     type(gtp_tpfun_as_coeff) :: tpfuncoef
     type(gtp_tpfun_as_coeff) :: cfun
  end type gtp_tpfun2dat
!\end{verbatim}
!--------------------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_phasetuple
  INTEGER, parameter :: gtp_phasetuple_version=1
  TYPE gtp_phasetuple
! for handling a single array with phases and composition sets
! ixphase is phase index (often lokph), compset is composition set index
! ADDED also index in phlista (lokph) and phase_varres (lokvares) and
! nextcs which is nonzero if there is a higher composition set of the phase
! A tuplet index always refer to the same phase+compset.  New tuples with
! the same phase and other compsets are added at the end.
! BUT if a compset>1 is deleted tuples with higher index will be shifted down!
! CONFUSING ixphase is usually iph, phases in alphabetical order in phases
!           lokph is usually lokph, location in phlista
     integer lokph,compset,ixphase,lokvares,nextcs
! >>>>>>>>>>>> old     integer phaseix,compset,ixphase,lokvares,nextcs
  end TYPE gtp_phasetuple
!\end{verbatim}
  TYPE(gtp_phasetuple), target, allocatable :: PHASETUPLE(:)
! -----------------------------------------------------------------
! NOTE: if one wants to model bond energies beteween sites in a phase
! like in a 3 sublattice sigma one can enter parameters like
! G(sigma,A:B:*) which will mean the bond energy between an A atom in
! first sublattice and B in the second.  The parameter G(sigma,B:A:*)
! will be different.  Such parameters, multiplied with the fractions of
! the constutuents, are added to the Gibbs energy even if there are 
! also endmember parameters like G(sigma,A:B:C)
! -----------------------------------------------------------------
!\begin{verbatim}
! a smart way to have an array of pointers used in gtp_phase
  TYPE endmemrecarray 
     type(gtp_endmember), pointer :: p1
  end TYPE endmemrecarray
!-----------------------------------------------------------------
! this constant must be incremented when a change is made in gtp_phase
  INTEGER, parameter :: gtp_phase_version=1
  TYPE gtp_phaserecord
! this is the record for phase model data. It points to many other records.
! Phases are stored in order of creation in phlista(i) and can be found
! in alphabetical order through the array phases(i)
! sublista is now removed and all data included in phlista
! sublattice and constituent data (they should be merged)
! The constituent link is the index to the splista(i), same function
! as LOKSP in iws.  Species in alphabetcal order is in species(i)
! One can allocate a dynamic array for the constituent list, done
! by subroutine create_constitlist.
! Note that the phase has a dynamic status word status2 in gtp_phase_varres
! which can be differnt in different parallel calculations.
! This status word has the FIX/ENT/SUS/DORM status bits for example
! name: phase name, note composition sets can have pre and suffixes
! model: free text
! phletter: G for gas, L for liquid
! alphaindex: the alphabetcal order of the phase (gas and liquids first)
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
! noofsubl: number if sublattices
! tnooffr: total number of fractions (constituents)
! linktocs: array with indices to phase_varres records
! nooffr: array with number of constituents in each sublattice 
! Note that sites are stored in phase_varres as they may vary with the
! constitution for ionic liquid)
     integer noofsubl,tnooffr
     integer, dimension(9) :: linktocs
     integer, dimension(:), allocatable :: nooffr
! number of sites in phase_varres record as it can vary with composition
! constitlist: indices of species that are constituents (in all sublattices)
     integer, dimension(:), allocatable :: constitlist
! used in ionic liquid:
! i2slx(1) is index of Va, i2slx(2) is index if last anion (both can be zero)
     integer, dimension(2) :: i2slx
! allocated in init_gtp.
  END TYPE gtp_phaserecord
! NOTE phase with index 0 is the reference phase for the elements
! allocated in init_gtp
  TYPE(gtp_phaserecord), private, allocatable :: phlista(:)
  INTEGER, private, allocatable :: PHASES(:)
!\end{verbatim}
!-----------------------------------------------------------------
! data for liquid phase with mqmqa model (only one!)
  TYPE gtp_mqmqa
! contains special information for liquid models with MQMQA
! nconst is umber of constituents, ncon1 on 1st subl, ncon2 on 2nd subl
! contyp(1..4,const) 1,2 species in first sublattibe, -1,-2 in second sublattice
! contyp(5,const) non-zero for endmember
! contyp(6,7,const) for endmembers specifies 2 species
! contyp(6..9,const) specifies 2 or 4 endmembers of the quadrupole
! contyp(10,const) index of fraction (species in alphabetical order)
! constoi(4,const) real with stoichiometry of species in quadrupole
     integer nconst,ncon1,ncon2
     integer, allocatable, dimension(:,:) :: contyp
     double precision, allocatable, dimension(:,:) :: constoi
  end TYPE gtp_mqmqa
!  TYPE(gtp_mqmqa), private :: mqmqa_data
! it should be private when everything work and can be removed from pmon6
  TYPE(gtp_mqmqa) :: mqmqa_data
!
!===================================================================
!
! below here are data structures for equilibrium description
!
!===================================================================
!
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_state_variable
  INTEGER, parameter :: gtp_state_variable_version=1
  TYPE gtp_state_variable
! this is to specify a formal or real argument to a function of state variables
! statevarid/istv: state variable index >=9 is extensive
! phref/iref: if a specified reference state (for chemical potential
! unit/iunit: 100 for percent, no other defined at present
! argtyp together with the next 4 integers represent the indices(4), only 0-4
! argtyp=0: no indices (T or P)
! argtyp=1: component
! argtyp=2: phase and compset
! argtyp=3: phase and compset and component
! argtyp=4: phase and compset and constituent
! ?? what is norm ?? normalizing like M in HM ?
     integer statevarid,norm,unit,phref,argtyp
! these integers represent the previous indices(4)
     integer phase,compset,component,constituent
! a state variable can be part of an expression with coefficients
! the coefficient can be stored here.  Default value is unity.
! In many cases it is ignored
     double precision coeff
! NOTE this is also used to store a condition of a fix phase
! In such a case statev is negative and the absolute value of statev
! is the phase index.  The phase and compset indices are also stored in
! "phase" and "compset" ??
! This is a temporary storage of the old state variable identifier
     integer oldstv
  end TYPE gtp_state_variable
! used for state variables/properties in various subroutines
!\end{verbatim}
! statevarid: defined in decode_state_variable3 in gtp3F.F90
! potentials: 1=T;   2=P;   3=MU;  4=AC;  5=LNAC
! extensive:  6=U;   7=S;   8=V;   9=H;  10=A;    11=G;
! phase:     12=NP; 13=BP; 14=Q ; 15=DG
! amounts:   16=N;  17=X;  18=B;  19=W;  20=Y
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_condition
! NOTE on unformatted SAVE files the conditions are written as texts
  INTEGER, parameter :: gtp_condition_version=1
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
! seqz is a sequential index of conditions, used for axis variables
! experimettype: inequality (< 0 or > 0) and/or percentage (-101, 100 or 101)
! symlink: index of symbol for prescribed value (1) and uncertainty (2)
! condcoeff: there is a coefficient and set of indices for each term
! prescribed: the prescribed value
! NOTE: if there is a symlink value that is the prescribed value
! current: the current value (not used?)
! uncertainty: the uncertainty (for experiments)
     integer :: noofterms,statev,active,iunit,nid,iref,seqz,experimenttype
!    TYPE(putfun_node), pointer :: symlink1,symlink2
! better to let condition symbol be index in svflista array
     integer symlink1,symlink2
     integer, dimension(:,:), allocatable :: indices
     double precision, dimension(:), allocatable :: condcoeff
     double precision prescribed, current, uncertainty
! confusing with record statevar and integer statev
     TYPE(gtp_state_variable), dimension(:), allocatable :: statvar
     TYPE(gtp_condition), pointer :: next, previous
  end TYPE gtp_condition
! used inside the gtp_equilibrium_data record and elsewhere
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_putfun_lista
  INTEGER, parameter :: gtp_putfun_lista_version=2
  TYPE gtp_putfun_lista
! these are records for STATE VARIABLE FUNCTIONS.  The function itself
! is handelled by the putfun package.
! linkpnode: pointer to start node of putfun expression
! narg: number of symbols in the function
! nactarg: number of actual parameter specifications needed in call
!   (like @P, @C and @S
! status: can be used for various things
! status bit SVFVAL set means value evaluated only when called with mode=1
! SVCONST bit set if symbol is just a constant value (linknode is zero)
! eqnoval: used to specify the equilibrium the value should be taken from
!    (for handling what is called "variables" in TC, SVFEXT set also)
! SVFTPF set if symbol is a TP function, eqnoval is TPFUN index
! if SVIMPORT set then the symbol is set equal to a TP function (only value
!     no derivatives).  TP function index is in TPLINK
! if SVEXPORT set the the value of the symbol is copied to a TP function
!     (must be a constant).   TP function index is in TPLINK
! name: name of symbol
     integer narg,nactarg,status,eqnoval,tplink
     type(putfun_node), pointer :: linkpnode
     character name*16
! THIS IS OLY USED FOR CONSTANTS, VALUES ARE ALSO STORED IN CEQ%SVFUNRES
     double precision svfv
! this array has identification of state variable (and other function) symbols 
! It is allocated in various subroutines, maybe be allocatable? 2020-08-31/BoS
     integer, dimension(:,:), pointer :: formal_arguments
  end TYPE gtp_putfun_lista
! this is the global array with state variable functions, "symbols"
  TYPE(gtp_putfun_lista), dimension(:), allocatable :: svflista
! NOTE the value of a function is stored locally in each equilibrium record
! in array svfunres.
! The number of entered state variable functions. Used when a new one stored
  integer, private :: nsvfun
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_fraction_set
  INTEGER, parameter :: gtp_fraction_set_version=1
  TYPE gtp_fraction_set
! info about disordered fractions for some phases like ordered fcc, sigma etc
! latd: the number of sublattices added to first disordred sublattice
! ndd: sublattices for this fraction set, 
! tnoofxfr: number of disordered fractions
! tnoofyfr: same for ordered fractions (=same as in phlista).
! varreslink: index of disordered phase_varres, 
! totdis: 0 indicates no total disorder (sigma), 1=fcc, bcc or hcp
! id: parameter suffix, D for disordered
! dsites: number of sites in sublattices, disordred fractions stored in
!    another phase_varres record with index varreslink (above)
! splink: pointers to species record for the constituents
! nooffr: the number of fractions in each sublattice
! y2x: the conversion from sublattice constituents to disordered and
! dxidyj: are the the coeff to multiply the y fractions to get the disordered
!        xfra(y2x(i))=xfra(y2x(i))+dxidyj(i)*yfra(i)
! disordered fractions stored in the phase_varres record with index varreslink
! arrays originally declared as pointers now changed to allocatable
     integer latd,ndd,tnoofxfr,tnoofyfr,varreslink,totdis
     character*1 id
     double precision, dimension(:), allocatable :: dsites
     integer, dimension(:), allocatable :: nooffr
     integer, dimension(:), allocatable :: splink
     integer, dimension(:), allocatable :: y2x
     double precision, dimension(:), allocatable :: dxidyj
! formula unit factor needed when calculating G for disordered sigma etc
     double precision fsites
  END TYPE gtp_fraction_set
! these records are declared in the phase_varres record as DISFRA for 
! each composition set and linked from the phase_varres record
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this constant must be incremented when a change is made in gtp_phase_varres
! added quasichemical bonds
  INTEGER, parameter :: gtp_phase_varres_version=2
  TYPE gtp_phase_varres
! Data here must be different in equilibria representing different experiments
! or calculated in parallel or results saved from step or map.
! nextfree: In unused phase_varres record it is the index to next free record
!    The global integer csfree is the index of the first free record
!    The global integer highcs is the higest varres index used
! phlink: is index of phase record for this phase_varres record
! status2: has composition set status bits CSxyz
! phstate: indicate state: fix/stable/entered/unknown/dormant/suspended/hidden
!                           2   1      0        -1      -2      -3       -4
! phtupx: phase tuple index
     integer nextfree,phlink,status2,phstate,phtupx
! abnorm(1): moles of components per formula unit of the phase/composition set
! abnorm(2): mass of components per formula unit
! abnorm(3): moles atoms per formula unit (all abnorm set by SET_CONSTITUTION)
! prefix and suffix are added to the name for composition sets 2 and higher
     double precision, dimension(3) :: abnorm
     character*4 prefix,suffix
! constat: array with status word for each constituent, any can be suspended
! yfr: the site fraction array
! mmyfr: min/max fractions, negative is a minumum
! sites: site ratios (which can vary for ionic liquids)
     integer, dimension(:), allocatable :: constat
     double precision, dimension(:), allocatable :: yfr
     real, dimension(:), allocatable :: mmyfr
     double precision, dimension(:), allocatable :: sites
! for ionic liquid derivatives of sites wrt fractions (it is the charge), 
! 2nd derivates only when one constituent is vacancy
! 1st sublattice P=\sum_j (-v_j)*y_j + Qy_Va
! 2nd sublattice Q=\sum_i v_i*y_i
! dpqdy is the abs(valency) of the species, set in set_constitution
! for the vacancy it is the same as the number of sites on second subl.
! used in the minimizer and maybe elsewhere
     double precision, dimension(:), allocatable :: dpqdy
     double precision, dimension(:), allocatable :: d2pqdvay
! disfra: a structure describing the disordered fraction set (if any)
! for extra fraction sets, better to go via phase record index above
! this TYPE(gtp_fraction_set) variable is a bit messy.  Declaring it in this
! way means the record is stored inside this record.
     type(gtp_fraction_set) :: disfra
! ---
! stored calculated results for each phase (composition set)
! amfu: is amount formula units of the composition set (calculated result)
! netcharge: is net charge of phase
! dgm: driving force
! qcbonds: quasichemical bonds (NOT SAVED ON UNFORMATTED)
     double precision amfu,netcharge,dgm,qcbonds
! qcsro: current value of SRO (for quasichemical model)
     double precision, allocatable, dimension(:) :: qcsro
! Other properties may be that: gval(*,2) is TC, (*,3) is BMAG, see listprop
! nprop: the number of different properties (set in allocate)
! listprop(1): is number of calculated properties
! listprop(2:listprop(1)): identifies the property stored in gval(1,ipy) etc
!   2=TC, 3=BMAG. Properties defined in the gtp_propid record
     integer nprop
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
! saved values from last equilibrium for dot derivative calculations
     double precision, dimension(:,:), allocatable :: cinvy
     double precision, dimension(:), allocatable :: cxmol
     double precision, dimension(:,:), allocatable :: cdxmol
! terms added to G if bit CSADDG nonzero
     double precision, dimension(:), allocatable :: addg
! integer containing the iteration when invsaved updated
     integer invsavediter
! arrays to save time in calc_dgdyterms, do not need to be saved on unformatted
     double precision, dimension(:,:), allocatable ::invsaved
  END TYPE gtp_phase_varres
! this record is created inside the gtp_equilibrium_data record
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! this must be incremented when a change is made in gtp_equilibrium_data
  INTEGER, parameter :: gtp_equilibrium_data_version=1
  TYPE gtp_equilibrium_data
! this contains all data specific to an equilibrium like conditions,
! status, constitution and calculated values of all phases etc
! Several equilibria may be calculated simultaneously in parallel threads
! so each equilibrium must be independent 
! NOTE: the error code must be local to each equilibria!!!!
! During step and map each equilibrium record with results is saved
! values of T and P, conditions etc.
! Values here are normally set by external conditions or calculated from model
! local list of components, phase_varres with amounts and constitution
! lists of element, species, phases and thermodynamic parameters are global
! status: not used yet?
! multiuse: used for various things like direction in start equilibria
! eqno: sequential number assigned when created
! next: index of next free equilibrium record
!       also index of next equilibrium in a list during step/map calculation.
! eqname: name of equilibrium
! comment: a free text, for example reference for experimental data.
! tpval(1) is T, tpval(2) is P, rgas is R, rtn is R*T
! rtn: value of R*T
! weight: weight value for this experiment, default unity
!     integer status,multiuse,eqno,next
     integer status,multiuse,eqno,nexteq
     character eqname*24,comment*72
     double precision tpval(2),rtn
     double precision :: weight=one
! svfunres: the values of state variable functions valid for this equilibrium
     double precision, dimension(:), allocatable :: svfunres
! the experiments are used in assessments and stored like conditions 
! lastcondition: link to condition list
! lastexperiment: link to experiment list
     TYPE(gtp_condition), pointer :: lastcondition,lastexperiment
! components and conversion matrix from components to elements
! complist: array with components (species index or location)??
! compstoi: stoichiometric matrix of compoents relative to elements
! invcompstoi: inverted stoichiometric matrix
     TYPE(gtp_components), dimension(:), allocatable :: complist
     double precision, dimension(:,:), allocatable :: compstoi
     double precision, dimension(:,:), allocatable :: invcompstoi
! one record for each phase+composition set that can be calculated
! phase_varres: here all calculated data for the phases are stored
     TYPE(gtp_phase_varres), dimension(:), allocatable :: phase_varres
! index to the tpfun_parres array is the same as in the global array tpres 
! eq_tpres: here local calculated values of TP functions are stored
! should be allocatable, not a pointer
     TYPE(tpfun_parres), dimension(:), allocatable :: eq_tpres
! current values of chemical potentials stored in component record but
! duplicated here for easy acces by application software
     double precision, dimension(:), allocatable :: cmuval
! xconv: convergence criteria for constituent fractions and other things
! dgconv(1) is controlling decrease of DGM for unstable phases
! dgconv(2) not used (yet)
     double precision xconv,gdconv(2)
! delta-G value for merging gridpoints in grid minimizer
! smaller value creates problem for test step3.OCM, MC and austenite merged
!     double precision :: gmindif=-5.0D-2
! testing merging again 190604/BoS
     double precision :: gmindif=-1.0D-2
! maxiter: maximum number of iterations allowed
     integer maxiter
! CCI number of iterations needed for the equilibrium calculation
     integer conv_iter
! This is to store additional things not really invented yet ...
! It may be used in ENTER MANY_EQUIL for things to calculate and list
     character (len=80), dimension(:), allocatable :: eqextra
! this is to save a copy of the last calculated system matrix, needed ??
! to calculate dot derivatives, initiate to zero
     integer :: sysmatdim=0,nfixmu=0,nfixph=0
     integer, allocatable :: fixmu(:)
     integer, allocatable :: fixph(:,:)
     double precision, allocatable :: savesysmat(:,:)
! This is temporary data for EEC but must be separate for parallelization
! index of phase_varres for liquid
     integer eecliq
     double precision eecliqs
! temporary array to handle converge problems with change of stable phase set
     integer, dimension(:,:), allocatable :: phaseremoved
  END TYPE gtp_equilibrium_data
! The primary copy of this structures is declared globally as FIRSTEQ here
! Others may be created when needed for storing experimental data or
! for parallel processing. A global array of these are
  TYPE(gtp_equilibrium_data), dimension(:), allocatable, target :: eqlista
  TYPE(gtp_equilibrium_data), pointer :: firsteq
! This array of equilibrium records are used for storing results during
! STEP and MAP calculations.
  TYPE(gtp_equilibrium_data), dimension(:), allocatable :: eqlines
!\end{verbatim}
!-----------------------------------------------------------------
!\begin{verbatim}
! for each permutation in the binary interaction tree of an endmember one must
! keep track of the permutation and the permutation limit.
! It is not possible to push the value on pystack as one must remember
! them when changing the endmember permutation
! integer, parameter :: permstacklimit=150
! this constant must be incremented when a change is made in gtp_parcalc
  INTEGER, parameter :: gtp_parcalc_version=1
  TYPE gtp_parcalc
! This record contains temporary data that must be separate in different
! parallel processes when calculating G and derivatives for any phase.
! There is nothing here that need to be saved after the calculation is finished
! global variables used when calculating G and derivaties
! sublattice with interaction, interacting constituent, endmember constituents
! PRIVATE inside this structure not liked by some compilers....
! endcon must have maxsubl dimension as it is used for all phases
     integer :: intlat(maxinter),intcon(maxinter),endcon(maxsubl)
! interaction level and number of fraction variables
     integer :: intlevel,nofc
! interacting constituents (max 4) for composition dependent interaction
! iq(j) indicate interacting constituents
! for binary RK+Muggianu iq(3)=iq(4)=iq(5)=0
! for ternary Muggianu in same sublattice iq(4)=iq(5)=0
! for reciprocal composition dependent iq(5)=0
! 2020/BoS not used: Toop, Kohler and simular iq(5) non-zero (not implemented) 
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
! this constant must be incremented when a change is made in gtp_pystack
  INTEGER, parameter :: gtp_pystack_version=1
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
!===================================================================
!
! below here are data structures for various applications
! They indicate data that may need to be saved together with
! the thermodynamic data.  Exactly how this will be handelled
! will have to be solved later
!
!===================================================================
!
!-----------------------------------------------------------------
!\begin{verbatim}
  INTEGER, parameter :: gtp_eqnode_version=1
  TYPE gtp_eqnode
! This record is to arrange calculated equilibria, for example results
! from a STEP or MAP calculation, in an ordered way.  The equilibrium records
! linked from an eqnode record should normally represent one or more lines
! in a diagram but may be used for other purposes.
! ident is to be able to find a specific node
! nodedtype is to specify invariant, middle, end etc.
! status can be used to supress a line
! color can be used to sepecify color or linetypes (dotted, thick ... etc)
! exits are the number of lines that should exit from the node
! done are the number of calculated lines currently exiting from the node
     integer ident,nodetype,status,color,exits,done
! this node can be in a multilayerd list of eqnodes
     type(gtp_eqnode), pointer :: top,up,down,next,prev
! nodeq is a pointer to the equilibrium record at the node
     type(gtp_equilibrium_data), pointer :: nodeq
! eqlista are pointers to line of equilibria starting or ending at the node
! The equilibrium records are linked with a pointer inside themselves
     type(gtp_equilibrium_data), dimension(:), pointer :: eqlista
! axis is the independent axis variable for the line, negative means decrement
! noeqs gives the number of equilibria in each eqlista, a negative value
! indicates that the node is an endpoint (each line normally has a
! start point and an end point)
     integer, dimension(:), allocatable :: axis,noeqs
! This is a possibility to specify a status for each equilibria in each line
!    integer, dimension(:,:), allocatable :: eqstatus
  end TYPE gtp_eqnode
! can be allocated in a gtp_applicationhead record
!\end{verbatim}
!------------------------------------------------------------------
!\begin{verbatim}
! a smart way to have an array of pointers used in gtp_assessmenthead
  TYPE equilibrium_array
     type(gtp_equilibrium_data), pointer :: p1
  end TYPE equilibrium_array
  INTEGER, parameter :: gtp_assessment_version=1
  TYPE gtp_assessmenthead
! This record should summarize the essential information about assessment data
! using GTP.  How it should link to other information is not clear.  
! status is status word, AHCOEF is used
! varcoef is the number of variable coefficients
! firstexpeq is the first equilibrium with experimental data
! lwam is allocated workspace at last call to lmdif1
     integer status,varcoef,firstexpeq,lwam
     character*64 general,special
     type(gtp_assessmenthead), pointer :: nextash,prevash
! This is list of pointers to equilibria to be used in the assessnent
! size(eqlista) is the number of equilibria with experimental data
     type(equilibrium_array), dimension(:), allocatable :: eqlista
! These are the coefficients values that are optimized,
! current values, scaling, start values, RSD and optionally min and max
     double precision, dimension(:), allocatable :: coeffvalues
     double precision, dimension(:), allocatable :: coeffscale
     double precision, dimension(:), allocatable :: coeffstart
     double precision, dimension(:), allocatable :: coeffrsd
     double precision, dimension(:), allocatable :: coeffmin
     double precision, dimension(:), allocatable :: coeffmax
! These are the corresponding TP-function constants indices
     integer, dimension(:), allocatable :: coeffindex
! This array indicate currently optimized variables:
!  -1=unused, 0=fix, 1=fix with min, 2=fix with max, 3=fix with min and max
!  10=optimized, 11=opt with min, 12=opt with max, 13=opt with min and max
     integer, dimension(:), allocatable :: coeffstate
! Work arrays ...
     double precision, dimension(:), allocatable :: wopt
  end TYPE gtp_assessmenthead
! this record should be allocated for assessments when necessary
  type(gtp_assessmenthead), allocatable :: ashrecord
!  type(gtp_assessmenthead), pointer :: firstash,lastash
! but this is later allocated, to avoid memory loss ashrecord should be used
! and then this pointer should be set to that record
  type(gtp_assessmenthead), pointer :: firstash
!\end{verbatim}
!------------------------------------------------------------------
!\begin{verbatim}
  INTEGER, parameter :: gtp_applicationhead_version=1
  TYPE gtp_applicationhead
! This record should summarize the essential information about an application
! using GTP.  How it should link to other information is not clear.  
! The character variables should be used to indicate that.
     integer apptyp,status
     character*64 general,special
! These can be used to define axis and other things
     integer, dimension(:), allocatable :: ivals
     double precision, dimension(:), allocatable :: rvals
     character*64, dimension(:), allocatable :: cvals
     type(gtp_applicationhead), pointer :: nextapp,prevapp
! The headnode can be the start of a structure of eqnodes with lines
     type(gtp_eqnode) :: headnode
! this is the start of a list of nodes with calculated lines or
! single equilibria that belong to the application.
     type(gtp_eqnode), dimension(:), allocatable :: nodlista
  end TYPE gtp_applicationhead
! this record is allocated when necessary
  type(gtp_applicationhead), pointer :: firstapp,lastapp
!\end{verbatim}
!===================================================================
!
! Below are private global variables like free lists etc.
!
!===================================================================
!
! Several arrays with lists have a free list: csfree,addrecs,eqfree,reffree
! it is not really consistent how to handle deleted equilibria etc
! as the eqlista or phase_varres arrays  may have "holes" with deleted data
!
!\begin{verbatim}
! counters for elements, species and phases initiated to zero
  integer, private :: noofel=0,noofsp=0,noofph=0
! counter for phase tuples (combination of phase+compset)
  integer, private :: nooftuples=0
! counters for property and interaction records, just for fun
  integer, private :: noofprop,noofint,noofem
! free lists in phase_varres records and addition records
  integer, private :: csfree,addrecs
! free list of references and equilibria
  integer, private :: reffree,eqfree
! maximum number of properties calculated for a phase
  integer, private :: maxcalcprop=20
! highcs is highest used phase_varres record (for copy equil etc)
  integer, private :: highcs
! Trace for debugging (not used)
  logical, private :: ttrace
! Output for debugging gridmin
  integer, private :: lutbug=0
  double precision :: proda=zero,privilege=zero
! minimum constituent fraction
  double precision :: bmpymin
! number of defined property types like TC, BMAG etc
  integer, private :: ndefprop
! this is the index of mobility data, set in init_gtp in subroutine gtp3A
  integer, private :: mqindex
! quasichemical model type, 1=classic, 2=corrceted type 1, 3=corrected type 2
  integer :: qcmodel=1
! this is to remember how manytimes find_gridmeen needs to search all gridp
  integer :: ngridseek
! this is to handle EEC in the grid minimizer NOT GOOD FOR PARALLELIZATION
!  integer :: neecgrid
  double precision :: sliqmax,sliqmin,gliqeec,sliqeec
! this is for warnings about using unkown model parameter identifiers
  integer, parameter :: mundefmpi=10
  integer nundefmpi
  character undefmpi(mundefmpi)*4
! this is set zero by new_gtp and incremented each time a Toop record
! is created in any phase
   integer uniqid
! this is to allow select_phases from database files
   integer nselph
   character (len=24), allocatable, dimension(:) :: seltdbph
! This is to indicate mobility parameters, no wildcared fractions allowed
   integer nowildcard(3)
!\end{verbatim}

! undocumented debug indicator
   integer :: gtpdebug=0

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

