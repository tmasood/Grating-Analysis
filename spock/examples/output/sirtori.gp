set term postscript eps color solid "Times-Roman" 20
set nokey
set time
set notime
set output "sirtori.eps"
set xlabel "X (um)"
set xrange [-1.0:6]
set ylabel "mode intensity (a.u)"
set y2label "refractive index"
set ytics border nomirror norotate autofreq 
set y2tics border nomirror norotate autofreq 
set title "Plot of power outcoupled per pass (%) .vs. grating depth and outcoupler size"
set notitle
plot "sirtori.nf"  using 1:2 axes x1y1 title "OC size = 10um" with lines lw 4, \
"sirtori.dat" using 1:2 axes x1y2 title "OC size = 10um" with lines lw 2

