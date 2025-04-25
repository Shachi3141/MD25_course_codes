# plot_fes_phi.gnu
set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'fes_psi_projection.png'

set title "1D Projected Free Energy: F(ψ)"
set xlabel "Psi (ψ) [degrees]"
set ylabel "Free Energy (kJ/mol)"

plot 'fes_psi.dat' using 1:2 with lines lw 2 lc rgb "blue" title "F(ψ)"
