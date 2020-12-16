reset
#set term postscript color enhanced landscape "Helvetica" 24
#set term png truecolor enhanced font "Times,15"
set term pdf enhanced font "Helvetica,10"

cd 'G:\\LAT\\ExDiag\\data\\trotter\\'
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

LS = "10 12 14 16"
dts = "5.00E-02 2.50E-02"
hs = "1.0000 6.0000"

set xlabel 't'

do for [iscale = 1:2] {
	if(iscale==1){
		set logscale y
		suffix = "_log"
		set key bottom right
	}else{
		unset logscale y
		suffix = ""
		set key top left
	}
	do for [L in LS] {
 
		set title "L = ".L

		set ylabel '{/Symbol e}_T/({/Symbol d}t)^2'
		set out    'G:\\LAT\\ExDiag\\doc\\plots\\trotter_err_L'.L.suffix.'.pdf'
		 
		plot \
		for [ih = 1:words(hs)] \
		for [idt = 1:words(dts)] \
		'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'.dat' using ($1):(($2-$3)/word(dts, idt)**2):(($2+$3)/word(dts, idt)**2) \
		title 'h='.substr(word(hs, ih),1,3).', {/Symbol d}t='.word(dts, idt) with filledcurves ls (2*ih+idt+1) fs transparent solid 0.2 ,\
		for [ih = 1:words(hs)] \
		for [idt = 1:words(dts)] \
		'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'.dat' using ($1):($2/word(dts, idt)**2) \
		notitle with lines ls (2*ih+idt+1) lw 3

		if(iscale==1){
			set ylabel '{/Symbol e}_U'
			set out    'G:\\LAT\\ExDiag\\doc\\plots\\norm_err_L'.L.suffix.'.pdf'
		 
			plot \
			for [ih = 1:words(hs)] \
			for [idt = 1:words(dts)] \
			'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'.dat' using ($1):($4-$5):($4+$5) \
			title 'h='.substr(word(hs, ih),1,3).', {/Symbol d}t='.word(dts, idt) with filledcurves ls (2*ih+idt+1) fs transparent solid 0.2 ,\
			for [ih = 1:words(hs)] \
			for [idt = 1:words(dts)] \
			'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'.dat' using ($1):($4) \
			notitle with lines ls (2*ih+idt+1) lw 3
		}
	}
}


unset title
set ylabel '{/Symbol e}_T/({/Symbol d}t)^2'
set out    'G:\\LAT\\ExDiag\\doc\\plots\\trotter_err_vs_L.pdf'

set key top left maxrow 4
 
plot \
for [ih = 1:words(hs)] \
for [il = 1:words(LS)] \
'trotter_errs_L'.word(LS, il).'_dt2.50E-02_h'.word(hs, ih).'.dat' using ($1):(($2-$3)/0.025**2):(($2+$3)/0.025**2) \
title 'h='.substr(word(hs, ih),1,1).',L='.word(LS, il) with filledcurves ls (4*ih+il+1) fs transparent solid 0.1 ,\
for [il = 1:words(LS)] \
for [ih = 1:words(hs)] \
'trotter_errs_L'.word(LS, il).'_dt2.50E-02_h'.word(hs, ih).'.dat' using ($1):($2/0.025**2) \
notitle with lines ls (4*ih+il+1) lw 3

L=12
set ylabel '{/Symbol e}_T/({/Symbol d}t)^2'
set out    'G:\\LAT\\ExDiag\\doc\\plots\\trotter_err_L12_long.pdf'
		 
plot \
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'_long.dat' using ($1):(($2-$3)/word(dts, idt)**2):(($2+$3)/word(dts, idt)**2) \
title 'h='.substr(word(hs, ih),1,3).', {/Symbol d}t='.word(dts, idt) with filledcurves ls (2*ih+idt+1) fs transparent solid 0.1,\
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'_long.dat' using ($1):($2/word(dts, idt)**2) \
notitle with lines ls (2*ih+idt+1) lw 3

##########################################################################

set ylabel '{/Symbol e}_T/{/Symbol d}t^2'
set out    'G:\\LAT\\ExDiag\\doc\\plots\\trotter_err_L12_frobenius.pdf'

af = sqrt(924)
		 
plot \
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'_frobenius.dat' using ($1):(af*($2-$3)/word(dts, idt)**2):(af*($2+$3)/word(dts, idt)**2) \
title 'h='.substr(word(hs, ih),1,3).',{/Symbol d}t='.word(dts, idt) with filledcurves ls (2*ih+idt+1) fs transparent solid 0.1,\
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'_frobenius.dat' using ($1):(af*$2/word(dts, idt)**2) \
notitle with lines ls (2*ih+idt+1) lw 3,\
\
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'.dat' using ($1):($2/word(dts, idt)**2) \
notitle with lines ls (2*ih+idt+1) lw 3 dt (1,1)

############################################################################

set xlabel 't/10^3'

set ylabel '{/Symbol e}_T'
set out    'G:\\LAT\\ExDiag\\doc\\plots\\trotter_err_L12_superlong.pdf'
		 
plot \
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'_superlong.dat' using ($1/1000):(($2-$3)):(($2+$3)) \
title 'h='.substr(word(hs, ih),1,3).', {/Symbol d}t='.word(dts, idt) with filledcurves ls (2*ih+idt+1) fs transparent solid 0.1,\
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'_superlong.dat' using ($1/1000):($2) \
notitle with lines ls (2*ih+idt+1) lw 3

set ylabel '{/Symbol e}_T'
set out    'G:\\LAT\\ExDiag\\doc\\plots\\norm_err_L12_superlong.pdf'
		 
plot \
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'_superlong.dat' using ($1/1000):(($4-$5)):(($4+$5)) \
title 'h='.substr(word(hs, ih),1,3).', {/Symbol d}t='.word(dts, idt) with filledcurves ls (2*ih+idt+1) fs transparent solid 0.1,\
for [ih = 1:words(hs)] \
for [idt = 1:words(dts)] \
'trotter_errs_L'.L.'_dt'.word(dts, idt).'_h'.word(hs, ih).'_superlong.dat' using ($1/1000):($4) \
notitle with lines ls (2*ih+idt+1) lw 3

#/({/Symbol d}t)^2
#/word(dts, idt)**2

set out 'C:\\Temp\\1.pdf'
plot sin(x) notitle 