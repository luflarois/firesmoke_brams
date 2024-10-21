set datafile separator ','
set xdata time
set timefmt "%Y-%m-%dT%H:%M:%S"
set key autotitle columnhead
set ylabel "Temp [C]" 
set xlabel 'Tempo'
set y2tics
set ytics nomirror
set style line 100 lt 1 lc rgb "grey" lw 0.5 # linestyle for the grid
set grid ls 100 # enable grid with specific linestyle
set ytics 0.5 # smaller ytics
set xtics 1   # smaller xtics
set title "Temperatura - 01Dez2020 00h - 24h de simulacao - Delta T=1seg"
set key title "Niveis em [m]"
plot "Guarani_TEMPC_2020120100.csv" using 1:2 with lines, '' using 1:3 with lines, '' using 1:4 with lines, '' using 1:5 with lines, '' using 1:6 with lines, '' using 1:7 with lines, '' using 1:8 with lines, '' using 1:9 with lines