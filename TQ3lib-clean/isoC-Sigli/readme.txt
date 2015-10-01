Application software example in C++ for OCTQ version 3

C Sigli and B Sundman

The example shows calculation of phases precipitating during cooling
of an Al alloy.  A free database cost507R.TDB and subset thereof,
507ss.TDB is provided.  

The COST507 database is a collection of assessments made in a COST
project some 20 year ago.  Some of the binary and ternary systems are
very good but this database is not recommended for multicomponent
alloy calculations.  If you want a high quality database you must buy
them from a commercial vendor.

To link the program using OpenMP for parallel eqxcution convert the
linkpm.txt to a Windows batch file (cmd) and execute it.  On UNIX
systems create a Makefile using the same commands.

Before you link it you must compile the main OC program in order to
have the library files liboceq.a and liboceqplus.mod.

The program will ask for a database file, provide cost507ss or the
larger cost507r or some database file you have yourself in a TDB
format.  It will ask for the major element and the alloying
composition.

Then it will calculate the liquidus and solidus temperatures and all
phases that will form during the cooling to room temperature.

This shows a typical calculation for an alloy 2Cu-2Mg-6Zn-Al:

 Enter the name of the TDB file (*.tdb) : cost507ss

 the following elements are in the database:
AL / CU / MG / ZN /
 give the name of the reference element (solvant for example) : Al

W[CU]=0 change it? y/n y
 value in weight %= 2

W[MG]=0 change it? y/n y
 value in weight %= 2

W[ZN]=0 change it? y/n y
 value in weight %= 6

 tqini created: DEFAULT_EQUILIBRIUM
  1:P=100000, 2:N=1, 3:W(AL)=0.9, 4:W(CU)=.02, 5:W(MG)=.02, 6:T=1500
 Degrees of freedom are   0

 =====================================
    New Equilibrium at : 1226.85 C
  1:P=100000, 2:N=1, 3:W(AL)=0.9, 4:W(CU)=.02, 5:W(MG)=.02, 6:T=1500
 Degrees of freedom are   0
 ---------------------------------------
 ---------------------------------------
         LIQUID fat%= 100
 ---------------------------------------
        AL = 94.1958 (at%)
        CU = 0.888802 (at%)
        MG = 2.32379 (at%)
        ZN = 2.59161 (at%)
number of loops8
CEQ_0 2
CEQ_1 3
CEQ_2 4
CEQ_3 5
CEQ_4 6
====================================================
         loop n:1 increment of T=-250
====================================================
treating transition : 1500
====================================================
         loop n:2 increment of T=-62.5
====================================================
treating transition : 1000
treating transition : 750.001
====================================================
         loop n:3 increment of T=-15.625
====================================================
treating transition : 937.502
treating transition : 875.002
treating transition : 687.502
treating transition : 625.002
====================================================
         loop n:4 increment of T=-3.90625
====================================================
treating transition : 921.878
treating transition : 828.128
treating transition : 640.628
treating transition : 609.378
====================================================
         loop n:5 increment of T=-0.976563
====================================================
treating transition : 910.16
treating transition : 828.129
treating transition : 636.723
treating transition : 609.379
====================================================
         loop n:6 increment of T=-0.244141
====================================================
treating transition : 908.208
treating transition : 825.2
treating transition : 634.771
treating transition : 608.403
====================================================
         loop n:7 increment of T=-0.0610352
====================================================
treating transition : 908.209
treating transition : 824.469
treating transition : 634.039
treating transition : 607.916
====================================================
         loop n:8 increment of T=-0.0152588
====================================================
treating transition : 908.027
treating transition : 824.348
treating transition : 633.918
treating transition : 607.917
====================================================
 Parallel: 1 number of threads: 5
====================================================
                     composition
 --------------------------------------------------
AL (w%): 0.9
CU (w%): 0.02
MG (w%): 0.02
ZN (w%): 0.06
 --------------------------------------------------
0 634.847 LIQUID +
1 634.831 LIQUID + FCC_A1 +
2 551.167 FCC_A1 +
3 360.738 ALCU_THETA + FCC_A1 +
4 334.752 ALCU_THETA + FCC_A1 + MGZN2 +
elapsed time (s)= 0.520501

