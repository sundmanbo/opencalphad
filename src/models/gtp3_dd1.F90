!
! Data structures for the TPFUN package
!
!=================================================================
! VARIABLES and STRUCTURES originally in TPFUN
! length of a function symbol
  integer, parameter :: lenfnsym=16
  integer, private :: freetpfun
!
! \begin{verbatim}
! ************* this declaration moved to metlib4
!  TYPE gtp_parerr
! This record contains the global error code.  In parallel processing each
! parallel processes has its own error code copied to this if nonzero
! it should be replaced by gtperr for separate errors in treads
!     INTEGER :: bmperr
!  END TYPE gtp_parerr
!  TYPE(gtp_parerr) :: gx
! needed to have error code as private in threads, also moved to metlib4
!--- $OMP  threadprivate(gx)
! \end{verbatim}
!-----------------------------------------------------------------
!
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
