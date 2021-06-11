reset
#set term postscript color enhanced landscape "Helvetica" 24
#set term png truecolor enhanced font "Times,15"
set term pdf enhanced font "Helvetica,10"

cd 'G:\\LAT\\ExDiag\\data\\OTOC\\'
pwd

set pointsize 0.5
set bar 0.5
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

L = 12
dt = 0.01
hs = "1.0000 6.0000"
bs = "10.0 20.0"

suffix  = sprintf("L%i_dt%2.2E", L, dt)

set xlabel 't'
set ylabel '-<[s_3(0,t),s_3(0,0)]^2>'

set out    'G:\\LAT\\ExDiag\\doc\\plots\\otoc_'.suffix.'.pdf'

plot \
for [ih = 1:words(hs)] \
for [ib = 1:words(bs)] \
's3s3_L'.L.sprintf("_beta%2.2E_dt%2.2E", word(bs, ib)+0.0, dt).'_h'.word(hs,ih).'.dat' using ($1):($2-$3):($2+$3) title 'h='.substr(word(hs, ih),1,3).', {/Symbol b}='.word(bs, ib) with filledcurves ls (2*(ib-1)+(ih-1)+2) fs transparent solid 0.2,\
for [ih = 1:words(hs)] \
for [ib = 1:words(bs)] \
's3s3_L'.L.sprintf("_beta%2.2E_dt%2.2E", word(bs, ib)+0.0, dt).'_h'.word(hs,ih).'.dat' using ($1):($2)            notitle                                                           with lines        ls (2*(ib-1)+(ih-1)+2) lw 3 



#for [ih = 1:words(hs)] \
#'s3s3_'.suffix.'_h'.word(hs,ih).'.dat' 

set out 'C:\\Temp\\1.pdf'
plot sin(x) notitle 
