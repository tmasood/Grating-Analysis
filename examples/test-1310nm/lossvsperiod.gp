set term postscript eps color solid "Times-Roman" 20
set nokey
set time
set notime
set yrange [-300:800]
set xrange [1750:4300]
set output "lossvsperiod.eps"
set xlabel "Grating Period (A)"
set ylabel "Alpha [Arb. Units]"
set title "Attenuation .vs. Grating Period for 600A Tooth Depth \
at wavelength = 1310nm"
set notitle
plot "lossvsprd.dat"  using 1:2 title "TD = 600A, lambda=1310nm" with lines lw 4
