#/usr/bin/python2

# gp in var names means GnuPlot

import os
import subprocess

import sympy

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
def plot3d(func, var1_range, var2_range, mode='cartesian', scale_x=10, scale_y=10):
    if type(func).__name__ in ('int', 'float'):
        func = sympy.Number(func)

    var1 = sympy.Symbol(var1_range[0])
    var2 = sympy.Symbol(var2_range[0])

    if mode == 'cartesian':
        step_x = (var1_range[-1] - var1_range[1])/scale_x
        step_y = (var2_range[-1] - var2_range[1])/scale_y
    elif mode == 'spherical':
        step_x = ((sympy.sin(var1_range[-1])*sympy.cos(var2_range[-1]) 
                   - sympy.sin(var1_range[1])*sympy.cos(var2_range[1])) / scale_x)
        step_y = ((sympy.sin(var1_range[-1])*sympy.sin(var2_range[-1]) 
                   - sympy.sin(var1_range[1])*sympy.sin(var2_range[1])) / scale_y)

    output = ''
    for i in range(scale_x):
        var1_val = float(i*step_x)
        for j in range(scale_y):
            var2_val = float(j*step_y)
            try:
                output += str(var1_val) + ' ' + str(var2_val) + ' ' + str(func.subs(var2, var2_val).subs(var1, var1_val)) + '\n'
            except AttributeError:
                return 'undef'
        output += '\n'
    
    return output

x = sympy.Symbol('x')
y = sympy.Symbol('y')

#gp_input = open('test.gnuplot', 'r').read()
gp_input = plot3d(1, ['x',0,10], ['y',0,10], 'spherical', 5, 5)

if gp_input != 'undef':
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
