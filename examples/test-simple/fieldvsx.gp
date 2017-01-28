set term postscript eps color solid "Times-Roman" 20
set nokey
set time
set notime
set xrange [-2.0:2.0]
set output "fieldvsx.eps"
set xlabel "X-axis (um)"
set ylabel "Field [Arb. Units]"
set title "Attenuation .vs. Grating Period for 600A Tooth Depth \
at wavelength = 1um"
set notitle
plot "sol"  using 1:2 title "lambda=1um" with lines lw 4
