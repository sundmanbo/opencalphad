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
set xrange [ 0: 1 ]
set yrange [ -10000: 0 ]
set title "Mixing enthalpy in Cu-Mg liquid"
# set key outside right
set key top left
#
# How to change pointtype in the middle?
#
plot "-" using 2:1 with points pt 5 ps 2 title "Calorimetry - Batalin 1987",\
"" using (1-$2):1 with points pt 7 ps 2 title "Calorimetry - Sommer 1983"
  -3700  .1 
  -5610  .2 
  -6800  .3 
  -7600  .4 
  -7950  .5 
  -7700  .6 
  -6700  .7 
  -5100  .8 
  -2600  .9 
e
  -1900  .075 
  -3200  .13  
  -3500  .15  
  -4800  .21 
  -5600  .245 
  -5800  .27 
  -7000  .33 
  -6600  .33 
  -7500  .38 
  -7900  .425 
  -8100  .43 
  -8400  .47 
  -8300  .475 
  -8700  .505 
  -8600  .515 
  -8500  .52 
  -8900  .54 
  -8900  .54 
  -8500  .565 
  -9000  .59 
  -8900  .61 
  -8950  .635 
  -8850  .65 
  -8800  .67 
  -8650  .685 
  -8550  .7 
e
pause mouse

