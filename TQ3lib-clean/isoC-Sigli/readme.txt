Application software example for OCTQ version 3 in C++

C Sigli and B Sundman

The example shows calculation of phases precipitating during cooling
of an Al alloy.  A free database cost507R.TDB and subset thereof,
507ss.TDB is provided.  

The COST507 database is a collection of assessments made in a COST
project some 20 year ago.  Some of the binary and ternary systems are
very good but this database is not recommended for multicomponent
alloy calculations.  If you want high quality databases you must buy
them from a commercial vendor.

To link the program convert the linksigli.txt to a Windows batch file
(cmd) and execute it.  On UNIX systems create a Makefile using the
same commands.  If you do not have Open MP installed remove the
directive -fopenmp.

Before you link it you must compile the main OC program in order to
have the library files liboceq.a and liboceqplus.mod.

The program will ask for a database file, provide cost507ss or the
larger cost507r or some database file you have yourself.  It will ask
for the major element and the alloying composition.

Then it will calculate the liquidus and solidus temperatures and all
phases that will form during the cooling to room temperature.

This shows a typical calculation for an alloy 2Cu-2Mg-6Zn-Al:

 Enter the name of the TDB file (*.tdb) : cost507ss

 the following elements are in the database:
AL / CU / MG / ZN /
 give the name of the reference element (solvant for example) : al

W[CU]=0 change it? y/n y
 value in weight %= 2

W[MG]=0 change it? y/n y
 value in weight %= 2

W[ZN]=0 change it? y/n y
 value in weight %= 6

 tqini created: DEFAULT_EQUILIBRIUM
  1:P=100000, 2:N=1, 3:W(AL)=0.9, 4:W(CU)=.02, 5:W(MG)=.02, 6:T=1500
 Degrees of freedom are   0
 Composition set(s) created:            3
 Grid minimization:     684 gridpoints   3.1200E-02 s and      20 clockcycles
Phase change: its/add/remove:     5    0    1
Phase change: its/add/remove:    20    0   37
Phase change: its/add/remove:    25    0   38

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
number of loops5
====================================================
         loop n:1 increment of T=-100
====================================================
treating transition : 1500.15
Phase change: its/add/remove:     5   24    0
Phase change: its/add/remove:     5    0    1
Phase change: its/add/remove:     5    8    0
Phase change: its/add/remove:    10   35    0
====================================================
         loop n:2 increment of T=-10
====================================================
treating transition : 1000.15
Phase change: its/add/remove:     5    0    8
Phase change: its/add/remove:    10    0   35
Phase change: its/add/remove:    15    1    0
Phase change: its/add/remove:    20    0   24
Phase change: its/add/remove:     5   24    0
treating transition : 900.15
Phase change: its/add/remove:     5    0    1
treating transition : 700.15
Phase change: its/add/remove:     5    8    0
Phase change: its/add/remove:     5   35    0
====================================================
         loop n:3 increment of T=-1
====================================================
treating transition : 910.15
Phase change: its/add/remove:     5    0    8
Phase change: its/add/remove:    10    0   35
Phase change: its/add/remove:    15    1    0
Phase change: its/add/remove:    20    0   24
Phase change: its/add/remove:     5   24    0
treating transition : 830.15
Phase change: its/add/remove:     5    0    1
treating transition : 640.15
Phase change: its/add/remove:     5    8    0
treating transition : 610.15
Phase change: its/add/remove:     5   35    0
====================================================
         loop n:4 increment of T=-0.1
====================================================
treating transition : 908.15
Phase change: its/add/remove:     5    0    8
Phase change: its/add/remove:    10    0   35
Phase change: its/add/remove:    15    1    0
Phase change: its/add/remove:    20    0   24
Phase change: its/add/remove:     5   24    0
treating transition : 825.15
Phase change: its/add/remove:     5    0    1
treating transition : 634.15
Phase change: its/add/remove:     5    8    0
treating transition : 608.15
Phase change: its/add/remove:     5   35    0
====================================================
         loop n:5 increment of T=-0.01
====================================================
treating transition : 908.05
Phase change: its/add/remove:     5    0    8
Phase change: its/add/remove:    10    0   35
Phase change: its/add/remove:    15    1    0
Phase change: its/add/remove:    20    0   24
Phase change: its/add/remove:     5   24    0
treating transition : 824.35
Phase change: its/add/remove:     5    0    1
treating transition : 633.95
Phase change: its/add/remove:     5    8    0
treating transition : 607.95
Phase change: its/add/remove:     5   35    0
====================================================
                     composition
 --------------------------------------------------
AL (w%): 0.9
CU (w%): 0.02
MG (w%): 0.02
ZN (w%): 0.06
 --------------------------------------------------
Phase change: its/add/remove:     5    0    8
Phase change: its/add/remove:    10    0   35
Phase change: its/add/remove:    15    1    0
Phase change: its/add/remove:    22    0   24
Phase change: its/add/remove:     5   24    0
Phase change: its/add/remove:     6    0    1
Phase change: its/add/remove:     5    8    0
Phase change: its/add/remove:     5   35    0
0 634.85 LIQUID +
1 634.84 LIQUID + FCC_A1 +
2 551.17 FCC_A1 +
3 360.73 ALCU_THETA + FCC_A1 +
4 334.75 ALCU_THETA + FCC_A1 + MGZN2 +
elapsed time (s)= 0.930501

C:\Users\bosse\arbete\OC\src\TQ3lib-test\isoC-Sigli>

