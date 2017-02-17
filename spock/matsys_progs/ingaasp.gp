set term postscript landscape color
set output "ingaasp.ps"
set xlabel "Wavelength (um)"
set ylabel "Refractive Index"
set title "Refractive index vs Wavelength"
plot "index.data" with lines
