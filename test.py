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
rp = e.Radial_Part(n=2, l=0)
#sympy.preview(rp)
p = sympy.Plot()
p.axes._label_axes = True
p.axes._stride[0] = 1
p[1] = (rp, [e.r,0,10,35])
rp1 = e.Radial_Part(n=2, l=1)
p[2] = (rp1, [e.r,0,10,35])
#sympy.pretty_print(wf)

#sympy.Plot(wf)
#f = sympy.exp( - sympy.sqrt(e.x**2 )) / sympy.sqrt( 1 - e.y**2/e.x**2)
#sympy.Plot(wf, [e.x,-1,1,30], [e.y,-1,1,30])

#sympy.Plot((1/(4*sympy.sqrt(2*sympy.pi))) * sympy.exp(-e.Z * sympy.sqrt(e.x**2)) * e.x * (e.Z/0.5)**(3.0/2.0)  , [e.x,-2,2])
