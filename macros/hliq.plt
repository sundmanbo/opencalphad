#
# One can use expressions to convert values in tables: 
# (1-$3) above means the x-value willl be "1-value in column 3"
#
# pt pointtype 1 +, 2 x, 3 star, 4 square, 5 fill square, 6 circle
#              7 filled circle, 8 triangle up, 9 filled triangle up
#             10 triangle down, 11 filled triangle down, 12 romb
#             13 filled romb, 14 pentad, 15 filled pentad, 16 same as 1
# ps pointsize
# color ???
#
#
# Mixing enthalpy in Cu-Mg liquid
#
set xlabel "x(Mg)"
set ylabel "Mixing enthalpy (J.mol-1)"
set title "Mixing enthalpy in Cu-Mg liquid"
set key bottom right
#
# How to change pointtype in the middle?
#
plot "-" using 2:3 with points pt 2 ps 1.5 title "Calorimetry - Batalin 1987",\
"" using 2:3  with points pt 8 ps 1.5 title "Caolorimetry - Hultgren",\
"" using 2:3  with points pt 3 ps 1.5 title "Calorimetry - Sommer 1983"
  5    1.022200E-01 -3.724000E+03  2
  5    2.080700E-01 -5.963900E+03  2
  5    3.013400E-01 -7.027000E+03  2
  5    4.011800E-01 -7.596000E+03  2
  5    5.011100E-01 -7.855800E+03  2
  5    6.011900E-01 -7.558900E+03  2
  5    7.014600E-01 -6.581600E+03  2
  5    8.019000E-01 -4.923900E+03  2
  5    9.024600E-01 -2.833300E+03  2
e
  5    1.948100E-01 -7.323100E+03  8
  5    3.039800E-01 -9.192200E+03  8
  5    3.972900E-01 -1.007000E+04  8
  5    4.908400E-01 -1.008100E+04  8
  5    5.975000E-01 -9.290400E+03  8
  5    7.075300E-01 -8.005100E+03  8
  5    8.047800E-01 -6.223300E+03  8
  5    8.957000E-01 -4.007700E+03  8
e
  5    9.250000E-01 -1.900000E+03  3
  5    8.700000E-01 -3.200000E+03  3
  5    8.500000E-01 -3.500000E+03  3
  5    7.900000E-01 -4.800000E+03  3
  5    7.550000E-01 -5.600000E+03  3
  5    7.300000E-01 -5.800000E+03  3
  5    6.700000E-01 -7.000000E+03  3
  5    6.700000E-01 -6.600000E+03  3
  5    6.200000E-01 -7.500000E+03  3
  5    5.750000E-01 -7.900000E+03  3
  5    5.700000E-01 -8.100000E+03  3
  5    5.300000E-01 -8.400000E+03  3
  5    5.250000E-01 -8.300000E+03  3
  5    4.950000E-01 -8.700000E+03  3
  5    4.850000E-01 -8.600000E+03  3
  5    4.800000E-01 -8.500000E+03  3
  5    4.600000E-01 -8.900000E+03  3
  5    4.600000E-01 -8.900000E+03  3
  5    4.350000E-01 -8.500000E+03  3
  5    4.100000E-01 -9.000000E+03  3
  5    3.900000E-01 -8.900000E+03  3
  5    3.650000E-01 -8.950000E+03  3
  5    3.500000E-01 -8.850000E+03  3
  5    3.300000E-01 -8.800000E+03  3
  5    3.150000E-01 -8.650000E+03  3
  5    3.000000E-01 -8.550000E+03  3
e

pause mouse

