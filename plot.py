#/usr/bin/python2

import subprocess
#import equations

gnuplot_settings='''set ticslevel 0
set style line 100 lt -1 lw 1
set pm3d hidden3d 100
set pm3d depthorder
unset hidden3d
unset colorbox
set palette model HSV functions 0.250+0.500*gray, 0.700, 0.800
set term pop
set xlabel "x"
set ylabel "y"
set zlabel "z"
set datafile missing "NaN"
splot "-" title "Function" with pm3d'''


