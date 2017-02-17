set term postscript eps color solid "Times-Roman" 20
set nokey
set time
set notime
set output "nearfield.eps"
set xlabel "X (um)"
set ylabel "Near field "
set title "Plot of power outcoupled per pass (%) .vs. grating depth and outcoupler size"
set notitle
plot "sirtori.nf"  using 1:2 title "OC size = 10um" with lines lw 4

