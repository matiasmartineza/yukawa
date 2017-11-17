set term pdf
set output 'image.pdf'
set title 'Tiempo M2P para distintos números de partículas para 4 threads' font ",12"

set xlabel "Número de partículas" font ",8" offset 0, 2
set ylabel "Tiempo M2P [s]" font ",8" offset 2, 0

set xtics font ",6"
set ytics font ",6"

plot "M2P_times" u 1:2 title 'Paralelizada' with lines lc rgb "blue", "M2P_times" u 1:3 title 'Serial' with lines lc rgb "red"
