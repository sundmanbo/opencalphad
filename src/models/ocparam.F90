
MODULE OCPARAM
	
    IMPLICIT NONE
    !----------------------------------------------------------------------
    ! Version numbers
    !----------------------------------------------------------------------
    ! version number of GTP (not OC)
    character*8, parameter :: gtpversion='GTP-3.30'
    ! THIS MUST BE CHANGED WHENEVER THE UNFORMATTED FILE FORMAT CHANGES!!!
    character*8, parameter :: savefile='OCF-3.20'
    !
    !----------------------------------------------------------------------
    !
    ! Parameters defining the size of arrays etc.
    ! max elements, species, phases, sublattices, constituents (ideal phase)
    ! NOTE increasing maxph to 600 and maxtpf to 80*maxph made the equilibrium
    ! record very big and created problems storing equilibria at STEP/MAP!!!
    integer, parameter :: maxel=100,maxsp=1000,maxph=600,maxsubl=10,maxconst=1000
    ! maximum number of constituents in non-ideal phase
    integer, parameter :: maxcons2=300
    ! maximum number of elements in a species
    integer, parameter :: maxspel=10
    ! maximum number of references
    integer, parameter :: maxrefs=1000
    ! maximum number of equilibria
    integer, parameter :: maxeq=900
    ! some dp values, default precision of Y and default minimum value of Y
    ! zero and one set in tpfun
    double precision, parameter :: YPRECD=1.0D-6,YMIND=1.0D-30
    ! dimension for push/pop in calcg, max composition dependent interaction
    integer, parameter :: maxpp=1000,maxinter=3
    ! max number of TP symbols, TOO BIG VALUE MAKES SAVE AT STEP/MAP DIFFICULT
    integer, parameter :: maxtpf=20*maxph
    !  integer, private, parameter :: maxtpf=80*maxph
    ! max number of properties (G, TC, BMAG MQ%(...) etc)
    integer, parameter :: maxprop=50
    ! max number of state variable functions
    integer, parameter :: maxsvfun=500

    double precision, parameter :: zero=0.0D0,one=1.0D0,two=2.0D0,ten=1.D1


    !----------------------------------------------------------------------
    ! Numerical parameters
    !----------------------------------------------------------------------
    integer, parameter :: default_splitsolver = 0
    ! 1 to allow to split the linear system when conditions leads to square matrix
    integer, parameter :: default_precondsolver = 0
    ! 1 to allow to use a Jacobi preconditionner for solving the linear system
    !
    !----------------------------------------------------------------------
    ! convergence criteria used in matsmin.F90
    !----------------------------------------------------------------------
    !
    !!!!!!!
    !! meq_phaseset subroutine
    !!!!!!!
    integer, parameter :: default_nochange = 4
    ! minimum number of iterations between a change of the set of stable phases
    ! should not be smaller than default_minadd/default_minrem
    integer, parameter :: default_minadd=4
    integer, parameter :: default_minrem=4
    !
    double precision, parameter :: default_addchargedphase = 1.D-2
    !Used to verify: charge(phase) > 1.0D-2
    !That is, checking that the phase to be added does not have a net charge
    !
    !!!!!!!
    !! meq_sameset subroutine
    !!!!!!!
    integer, parameter :: default_typechangephaseamount = 0
    ! By default, 0 leads to default_scalechangephaseamount=1.0
    ! 1 leads to default_scalechangephaseamount=sum of prescribed conditions -/+ 1
    ! 2 leads to default_scalechangephaseamount=max of (1, max of prescribed conditions)

    double precision, parameter :: default_scalechangephaseamount = 1.0
    ! scale all changes in phase amount with total number of atoms. By default,
    ! assume this is unity.


    double precision, parameter :: default_ylow = 1.D-3
    ! parameter added to avoid too drastic jumps in small site fractions
    ! normalizing factor, if y < ylow ....
    !
    double precision, parameter :: default_ymin = 1.D-12
    ! parameter added to avoid too drastic jumps in small site fractions
    !
    double precision, parameter :: default_ymingas = 1.D-30
    ! parameter added for gases, since the phase one must allow smaller constituent fractions
    ! normalizing factor, if y < critymingas then y = cirtymingas
    !
    double precision, parameter :: default_ionliqyfact = 3.D-1
    ! this is an emergecy fix to improve convergence for ionic liquid
    ! correction to site fractions in ionic liquids
    !
    double precision, parameter :: default_deltaTycond=2.5D1
    ! this is set each time the set of phases changes, controls change in T
    ! when there is a condition on y
    !
    integer, parameter :: default_nophasechange = 100
    !Criterion on the maximum number of iterations that should go by with no change in the set of phases
    ! That is, the system should have at least one phase change every default_nophasechange iterations
    !
    double precision, parameter :: default_maxphaseamountchange = 1.0E-10
    !Criterion on the minimum of amount of phase change (DeltaN) vis-a-vis slow convergence
    ! That is, if the set of stable phases doesn't change, and the change in stable phases is lower than
    ! default_maxphaseamountchange, then this is considered a 'slow convergence case'
    !
    double precision, parameter :: default_bigvalues = 1.0D+50
    !Criterion on the maximum element of smat matrix
    ! Most probably, if something in smat is bigger than default_bigvalues, a calculation error has occurred
    !
    double precision, parameter :: default_minimalchangesT = 1.0D-2
    ! minimal change in Temperature allowed when Temperature is variable
    !
    double precision, parameter :: default_limitchangesT = 0.2D0
    double precision, parameter :: default_deltaT = 1.0D1
    !modified xconv criterion CHECK
    !Used to verify:  DeltaT > default_delatT*%xconv (converged = 8)
    !Case where T is variable
    !
    double precision, parameter :: default_limitchangesP = 0.2D0
    double precision, parameter :: default_deltaP = 1.0D4
    !modified xconv criterion CHECK
    !Used to verify:   DeltaP > default_deltaP*%xconv (converged = 8)
    !Case where P is variable
    !
    double precision, parameter :: default_minimalchangesP = 1.0D-2
    ! minimal change in Pressure allowed when Pressure is variable
    !
    double precision, parameter :: default_chargefact = 1.0
    ! term added to the correction in site fraction due to electric charge
    !
    integer, parameter :: default_noremove=3
    !Criterion on the minimum number of iterations with
    ! N-DeltaN<0  Before removing the phase in question
    ! That is, the phase must have a negative quantity for default_noremove iterations
    ! before being removed
    !
    double precision, parameter :: default_yvar1 = 1.0D-4
    ! first limitation to change in site fraction
    !
    double precision, parameter :: default_yvar2 = 1.0D-13
    ! second limitation to change in site fraction
    !
    double precision, parameter :: default_upperyvar1 = 1.0D-3
    ! limitation to change in site fraction
    ! normalizing factor, if yvar1 > default_upperyvar1 ....
    !                     then yvar1 =  default_upperyvar1
    !
    double precision, parameter :: default_upperyvar2 = 1.0D-13
    ! limitation to change in site fraction
    ! normalizing factor, if yvar2 > default_upperyvar2 ....
    !                     then yvar2 =  default_upperyvar2
    !
    double precision, parameter :: default_correctionfactorYS = 1.0D1
    !multiplier in a criterion
    !Used to verify:
    !|Delta y(phase, constituent)(recursion =k)| > default_correctionfactorYS*|Delta y(phase, constituent)(recursion =k-1)|(converged = 3)
    !
    double precision, parameter :: default_correctionfactorXCONV = 1.0D2
    !multiplier in a criterion CHECK
    !Used to verify: In an unstable phase:
    !   Delta y(phase, constituent) > default_correctionfactorXCONV*%xconv(converged = 4)
    !
    double precision, parameter :: default_correctionfactorDGM = 1.0
    !default_correctionfactorDGM criterion
    !Used to verify:
    !   dgm(recursion=k) - dgm(recursion=k-1) > default_correctionfactorDGM *gdconv(1) (converged = 4)
    !Case where more than 10 constituents in the phases are present, (apparently) warranting a bigger gdconv(1)
    !
    double precision, parameter :: default_upperycormax2 = 1.0D-4
    ! check on max stepsize, determining whether or not it is too small
    !
    integer, parameter :: default_minimaliterations = 4
    !Criterion on the minimum number of iterations for the code as a whole
    !
    !!!!!!!
    !! userif/pmon6.F90, gtp3A.F90 Fortran files
    !!!!!!!
    double precision, parameter :: default_maxiter = 500
    ! default maximum number of iteration
    !
    double precision, parameter :: default_xconv = 1.D-6
    double precision, parameter :: default_minxconv = 1.D-30
    ! default and minimal values for ceq%xconv criterion
    !
    double precision, parameter :: default_mingdconv = 1.D-5
    double precision, parameter :: default_gdconv1 = 4.D-3
    double precision, parameter :: default_gdconv2 = 0.D0
    ! default and  minimal value for ceq%gdconv(1) criterion
    !
    double precision, parameter :: default_mingridmin = -1.D-2
    ! minimal value for ceq%gmindif criterion
    !

    !----------------------------------------------------------------------
    ! Physical parameters
    !----------------------------------------------------------------------
    double precision, parameter :: PI = 3.141592653589793D0
END MODULE OCPARAM

