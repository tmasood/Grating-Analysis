set cntrparam levels 10
set xrange [0:1]
set yrange [-0.2:0.8]
set view 60,80
#set auto
set zrange [0.0:1.0]
set data style lines
set title "3D gnuplot demo - contour of data grid plotting"
set xlabel "Lateral Position x"
set ylabel "Propagation Direction z"
set parametric
set pm3d; set palette
#splot "fields.dat"
#splot "fields15.dat"
splot "fields15-25.dat"
pause -1 "Hit return to continue (19)"
set pm3d solid
set noparametric
replot
pause -1 "Hit return to continue (20)"
set pm3d at b
replot
pause -1 "Hit return to continue (21)"
set pm3d at st solid
replot
pause -1 "Hit return to continue (22)"
unset surface
set pm3d at st solid
replot
pause -1 "Hit return to continue (23)"
set nocblabel
set view 50,50
set pm3d at bstbst
replot
pause -1 "Press Enter. (24)"
set pm3d transparent

set title "gray map"
set pm3d at b
set palette gray
set view 180,0,1.2; set yrange [*:*] reverse
set samples 100; set isosamples 100
replot
pause -1 "Press Enter."
set title "gray map, negative"
set pm3d at b
set palette gray negative
set view 180,0,1.2; set yrange [*:*] reverse
replot
pause -1 "Press Enter."

set title "colour map, using default rgbformulae 7,5,15 ... traditional pm3d (black-blue-red-yellow)"
set palette color positive
set view 180,0,1.2; set yrange [*:*] reverse
set samples 50; set isosamples 50
replot
pause -1 "Press Enter."

set title "colour, rgbformulae 3,11,6 ... green-red-violet"
set palette rgbformulae 3,11,6
replot
pause -1 "Press Enter."
set title "colour, rgbformulae 23,28,3  ... ocean (green-blue-white); OK are also all other permutations"
set palette rgbformulae 23,28,3
replot
pause -1 "Press Enter."

set title "colour, rgbformulae 30,31,32 ... color printable on gray (black-blue-violet-yellow-white)"
set palette rgbformulae 30,31,32
replot
pause -1 "Press Enter."

set title "rgbformulae 31,-11,32: negative formula number=inverted color"
set palette rgbformulae 31,-11,32
replot
pause -1 "Press Enter."
set yrange [*:*] noreverse



#
reset
