reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\ExDiag\\data\\'
set pointsize 2.5
set bar 2.5
set style line 1  lt 1 lc rgb '#000000' lw 4 pt 11
set style line 2  lt 1 lc rgb '#FF0000' lw 4 pt 12
set style line 3  lt 1 lc rgb '#00FF00' lw 4 pt 13
set style line 4  lt 1 lc rgb '#0000FF' lw 4 pt 14
set style line 5  lt 1 lc rgb '#FF00FF' lw 4 pt 15
set style line 6  lt 1 lc rgb '#FFA500' lw 4 pt 16
set style line 7  lt 1 lc rgb '#666666' lw 4 pt 3
set style line 11 lt 2 lc rgb '#000000' lw 4 pt 5
set style line 12 lt 2 lc rgb '#FF0000' lw 4 pt 7
set style line 13 lt 2 lc rgb '#00FF00' lw 4 pt 9
set style line 14 lt 2 lc rgb '#0000FF' lw 4 pt 11
set style line 15 lt 2 lc rgb '#FF00FF' lw 4 pt 13
set style line 16 lt 2 lc rgb '#00FFFF' lw 4 pt 15
set style line 17 lt 2 lc rgb '#888888' lw 4 pt 1

set xlabel 'h'
set ylabel 'R'
set out    'G:\\LAT\\ExDiag\\R_ratio.eps'
 
plot \
'R_L8.dat'  using ($1):($2):($3) title 'L=8'  with yerrorbars ls 2,\
'R_L10.dat' using ($1):($2):($3) title 'L=10' with yerrorbars ls 3,\
'R_L12.dat' using ($1):($2):($3) title 'L=12' with yerrorbars ls 4