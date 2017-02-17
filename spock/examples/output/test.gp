set term postscript landscape color "Times-Roman-bold" 
set time 
set output "test.ps" 
set xlabel "X (um)" 
set ylabel "Refractive index" 
set title "Plot of Refractive index .vs. X (um)" 
plot "test.dat" using 1:2 title "Refractive Index " with lines lw 4 
