This test example shows how to extract mobility data 

in the TDB file crfe+mob.data we have mobility data for BCC

 TYPE_DEFINITION & GES A_P_D BCC_A2 MAGNETIC  -1  0.4000!
 PHASE BCC_A2 %&  2   1.000   3.000 !
     CONSTITUENT BCC_A2 :CR FE: VA:!
     PARAMETER G(BCC_A2,CR:VA;0)  298.15 +GHSERCR+GPCRBCC; 6000 N REF283 !
     PARAMETER TC(BCC_A2,CR:VA;0)  298.15 -311.5; 6000 N REF281 !
     PARAMETER BMAG(BCC_A2,CR:VA;0)  298.15 -.01; 6000 N REF281 !
     PARAMETER MQ&FE#1(BCC_A2,CR:VA;0)  298.15 +1E-13+1E-17*T; 6000 N BOSSE !
     PARAMETER MQ&CR#1(BCC_A2,CR:VA;0)  298.15 +1E-09+1E-10*T; 6000 N BOSSE !
     PARAMETER G(BCC_A2,FE:VA;0)  298.15 +GHSERFE+GPFEBCC; 6000 N REF283 !
     PARAMETER TC(BCC_A2,FE:VA;0)  298.15 +1043; 6000 N REF281 !
     PARAMETER BMAG(BCC_A2,FE:VA;0)  298.15 +2.22; 6000 N REF281 !
     PARAMETER MQ&FE#1(BCC_A2,FE:VA;0)  298.15 +1E-12+1E-15*T; 6000 N BOSSE !
     PARAMETER MQ&CR#1(BCC_A2,FE:VA;0)  298.15 +1E-08+1E-10*T; 6000 N BOSSE !
    PARAMETER G(BCC_A2,CR,FE:VA;0)  298.15 +20500-9.68*T; 6000 N REF107 !
    PARAMETER TC(BCC_A2,CR,FE:VA;0)  298.15 +1650; 6000 N REF107 !
    PARAMETER TC(BCC_A2,CR,FE:VA;1)  298.15 +550; 6000 N REF107 !
    PARAMETER BMAG(BCC_A2,CR,FE:VA;0)  298.15 -0.85; 6000 N REF107 !

where 
MQ&FE(BCC,FE:VA) is the tracer mobility(?) of Fe in pure Fe 
MQ&FE(BCC,CE:VA) is the infinite dilute mobility(?) of Fe in pure Cr
MQ&CR(BCC,FE:VA) is the infinite dilute mobility(?) of Cr in pure Fe
MQ&CR(BCC,CR:VA) is the tracer mobility(?) of Cr in pure Fe 

We calculate an equilibrium at x(cr)=.3 and t?1200 (pure bcc)
and then extract the mobilites using tqcph2

Note we must use the array listprop to detect which property that is
mobilities, we have also two other properties, the Curie T and Bohr
mangeton numbers

