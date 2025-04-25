# plot_fes_2d.gnu
set terminal pngcairo size 800,700 enhanced font 'Verdana,12'
set output 'fes_phi_psi.png'

set title "Free Energy Surface (FES): φ vs ψ"
set xlabel "Phi (ϕ) [degrees]"
set ylabel "Psi (ψ) [degrees]"

set pm3d map
set palette rgbformulae 22,13,-31  # nice color gradient
set cblabel "Free Energy (kJ/mol)"


set key off
splot 'fes.dat10.dat' using 1:2:3 with pm3d

