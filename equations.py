#!/usr/bin/python2

import sympy
from sympy import Symbol

# Declaring variables
x = Symbol('x')
n = Symbol('n')
l = Symbol('l')
m = Symbol('m')
r = Symbol('r')
phi = Symbol('phi')
theta = Symbol('theta')

# PHI-equation for real numbers
# PHI = lambda m, phi: sympy.sqrt(1.0/(2 * sympy.pi)) * (sympy.cos(sympy.abs(m) * phi) if m>=0 else sympy.sin(sympy.abs(m) * phi))

#TODO: determine mode by phi type
def PHI (m=0, phi=Symbol('phi'), mode='numer'): 
    if mode == 'numer':
        lp = sympy.sqrt(1./(2. * sympy.pi)).evalf()
        rp = (sympy.cos(sympy.abs(m) * phi) if m>=0 else sympy.sin(sympy.abs(m) * phi)).evalf()
    elif mode == 'analit':
        lp = 1/sympy.sqrt(2 * sympy.pi)
        rp = sympy.cos(sympy.abs(m) * phi) if m>=0 else sympy.sin(sympy.abs(m) * phi)
    else:
        return False

    return lp * rp
    #return sympy.sqrt(1./(2. * sympy.pi)) * (sympy.cos(sympy.abs(m) * phi) if m>=0 else sympy.sin(sympy.abs(m) * phi))

# Legendre polinomials
#P = lambda n,x: (1./(2.**n * sympy.factorial(n))) * sympy.diff((x**2 - 1)**n, x, n)

def Legendre(n=0, x=Symbol('x')):
    t = type(x).__name__  
    if t == 'float' or t = 'int':
        rat_part = 1./(2.**n * sympy.factorial(n)).evalf()
        diff_part = sympy.diff((x**2 - 1)**n, x, n).evalf()
    elif x.is_Symbol or x.is_Add or x.is_Mul or x.is_Function:
        rat_part = 1/(2**n * sympy.factorial(n))
        diff_part = sympy.diff((x**2 - 1)**n, x, n)
    else:
        return False
    
    return (rat_part * diff_part).expand()

for i in range(0,10):
    print(Legendre(n=i, mode='analit').__str__() + ' | ' + sympy.legendre(i,x).__str__())

# generalized Legendre polinomials
gP = lambda n,m,x: (1 - (x**2))**(sympy.abs(m)/2) * sympy.diff(P(n,x), x, abs(m))

# THETA-equation solution
THETA = lambda l,m,theta: sympy.sqrt(sympy.Rational((2 * l + 1),2.0) \
                       * sympy.Rational(sympy.factorial(l - sympy.abs(m)),sympy.factorial(l + sympy.abs(m)))) \
                       * gP(l,m,theta).subs(theta, sympy.cos(theta))

#/*
# Angular part
#*/
#Y(l,m,theta,phi) := THETA(l, m, theta) * PHI(m, phi)$
#



#/*
# Angular part plot
#*/
#/*
#plot3d(Y(0,0,theta,phi), [theta,0,2*%pi], [phi,0,2*%pi], [transform_xy, make_transform([theta,phi,r],r*sin(phi)*sin(theta), r*cos(phi)*sin(theta),r*cos(theta))])$
#


## Generalized
#plot3d(Y(2,1,theta,phi), [theta,0,2*%pi], [phi,0,2*%pi], [transform_xy, make_transform([theta,phi,r],r*sin(phi)*sin(theta), r*cos(phi)*sin(theta),r*cos(theta))], [grid, 50, 50], [plot_format, gnuplot])$
#*/
#


#/*
# Laguerre polynomials
#*/
#L(n,r) := exp(r) * diff(r^n * exp(-r), r, n)$
#


#/*
# generalized Legendre polinomials
#*/
#aL(n,k,r) := diff(L(n, r), r, k)$
#


#/*
# Radial part
#*/
#R(n,l,r) := -(((2 * Z)/(n * a0))^3 * (n - l - 1)!/(2 * n * ((n+l)!)^3)) ^ (1/2) * exp(-(Z * r)/(a0 * n)) * at(aL(n+l, 2*l+1, t), t=(2*Z*r)/(a0*n))$
#


#/*
# wavefunction
#*/
#PSI(n,l,m,r,theta,phi) := R(n, l, r) * Y(l, m, theta, phi)$
#


########################################################################################################################################
# Temp
#R(n,l,r) := -(((2 * Z) / (n * a0))^3 * (factorial((n - l - 1)) / (2 * n * (factorial(n + l))^3)))^(1/2) * exp(-(Z * r)/(a0 * n)) * aL(n+l, 2*l+1, t)$
#kill(t)$ Rez:R(2,0,r)$ t:(2*Z*r)/(a0*n)$ ratsimp(ev(Rez));
#
