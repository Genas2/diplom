#/usr/bin/python2

# gp in var names means GnuPlot

import os
import subprocess
#import equations

gp_settings = '''set ticslevel 0
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
splot "-" title "Function" with pm3d
'''
gp_input = open('test.gnuplot', 'r').read()
gp_filename = 'temp.gnuplot'
gp_fd = open(gp_filename, 'w')
gp_fd.write(gp_settings)
gp_fd.write(gp_input)
gp_fd.close()

p = subprocess.Popen(['gnuplot', '-p', gp_filename])
p.wait()

try:
    os.remove(gp_filename)
except OSError as (errno, strerror):
    print(gp_filename + ': ' + strerror)
