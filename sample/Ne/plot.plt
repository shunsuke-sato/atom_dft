set term postscript eps enhanced color font "Arial, 24" dl 5
set encoding iso_8859_1
Ba = 0.529
eV = 27.2
set sample 600

set xrange [0:2]


set xlabel "Radius ({\305})"
set ylabel "Radial wavefunctions"

set output "Ne_wavefunctions.eps"
set title "Wavefunction [Ne]"
set yrange [*:*]
p "results.dat" u ($1*Ba):6 smooth cspline title "1s" w l lt 3 lc 1 lw 2\
,"" u ($1*Ba):(-$7) smooth cspline title "2s" w l lt 1 lc 1 lw 2\
,"" u ($1*Ba):($8) smooth cspline title "2p" w l lt 1 lc rgb "forest-green" lw 2\
,0.0 title "" w l lt 1 lc 7


set output "Ne_density.eps"
set title "Density [Ne]"
set ylabel "r^2{/Symbol r}(r) (a.u.)"
set yrange [*:*]
p "results.dat" u ($1*Ba):($2*$1**2) smooth cspline title "" w l lt 1 lc 1 lw 2\
,0.0 title "" w l lt 1 lc 7


set output "Ne_potential.eps"
set title "Potential [Ne]"
set ylabel "Energy (eV)"
set yrange [-200:200]
p "results.dat" u ($1*Ba):(($3+$4+$5)*eV) smooth cspline title "V_{KS}" w l lt 1 lc 1 lw 2\
,"" u ($1*Ba):($3*eV) smooth cspline title "-Z/r" w l lt 1 lc rgb "forest-green" lw 2\
,"" u ($1*Ba):($4*eV) smooth cspline title "Hartree" w l lt 1 lc 3 lw 2\
,"" u ($1*Ba):($5*eV) smooth cspline title "Exchange-correlation" w l lt 1 lc 4 lw 2\
,0.0 title "" w l lt 1 lc 7



unset output
