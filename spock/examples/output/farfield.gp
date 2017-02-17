set term postscript eps color solid "Times-Roman" 20
set nokey
set time
set notime
set output "farfield.eps"
set xlabel "X (um)"
set ylabel "Far field "
set title "Plot of power outcoupled per pass (%) .vs. grating depth and outcoupler size"
set notitle
plot "spencer.ff"  using 1:2 title "OC size = 10um" with lines lw 4

