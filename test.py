#/usr/bin/env python2

import equations
import sympy

e = equations.Equations()

wf = e.Wave_Function(n=2, l=1, m=1)
#wf = sympy.trigsimp(wf)
#wf = sympy.radsimp(wf)
#wf = sympy.ratsimp(wf)
#wf = sympy.powsimp(wf)

#sympy.pretty_print(wf)
#rp = e.Radial_Part(2, 1)
#sympy.preview(rp)

wf = wf.subs(e.Z, 1)
#wf = wf.subs(e.r, sympy.sqrt(e.x**2 + e.y**2 + e.z**2))
#wf = wf.subs(e.theta, sympy.atan(sympy.sqrt(e.x**2 + e.y**2)/e.z))
wf = wf.subs(e.r, sympy.sqrt(e.x**2 + e.y**2))
wf = wf.subs(e.theta, sympy.atan(sympy.sqrt(e.x**2 + e.y**2)*2))
wf = wf.subs(e.phi, sympy.atan(e.y/e.x))
wf = wf.subs(e.a0, sympy.Rational(1,2))
#wf = sympy.separate(wf)
#wf = wf.subs(e.z, 0.0000000000000000000000001)
wf = wf.subs(e.y, 1)

#sympy.preview(wf)
#sympy.preview(rp)
p = sympy.Plot()
p.axes._label_axes = True
p.axes._label_ticks = True
p.axes._stride = [1,0.5,0]

rp = e.Radial_Part(n=2, l=1)
p[1] = (e.r**2 * rp**2, [e.r,0,5,100])
#rp1 = e.Radial_Part(n=2, l=1)
#p[2] = (e.r**2 * rp1**2, [e.r,0,5,100])
#rp2 = e.Radial_Part(n=3, l=2)
#p[3] = (e.r**2 * rp2**2, [e.r,0,5,100])

#sympy.pretty_print(wf)

#sympy.Plot(wf)
#f = sympy.exp( - sympy.sqrt(e.x**2 )) / sympy.sqrt( 1 - e.y**2/e.x**2)
#sympy.Plot(wf, [e.x,-1,1,30], [e.y,-1,1,30])

#sympy.Plot((1/(4*sympy.sqrt(2*sympy.pi))) * sympy.exp(-e.Z * sympy.sqrt(e.x**2)) * e.x * (e.Z/0.5)**(3.0/2.0)  , [e.x,-2,2])
