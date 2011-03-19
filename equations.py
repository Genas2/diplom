#!/usr/bin/python2

import sympy
from sympy import Symbol, factorial, Rational

# Declaring variables
x = Symbol('x')
n = Symbol('n')
l = Symbol('l')
m = Symbol('m')
r = Symbol('r')
phi = Symbol('phi')
theta = Symbol('theta')

def select_exec_mode(x):
    t = type(x).__name__

    if t == 'float' or t == 'int':
        return 'numer'
    elif t == 'str':
        return 'custom_var'

    try:
        if x.is_Symbol or x.is_Add or x.is_Mul or x.is_Function:
            return 'analit'
    except AttributeError:
        return 'undef'

# PHI-equation for real numbers

def Phi_Equation (m=0, phi=Symbol('phi')): 
    m = float(int(m))
    mode = select_exec_mode(phi)

    if mode == 'numer': 
        phi_val = phi
        phi = Symbol('phi')
    elif mode == 'custom_var':
        phi = Symbol(phi)
    elif mode == 'undef':
        return False

    left_part  = 1/sympy.sqrt(2 * sympy.pi)
    right_part = sympy.cos(sympy.abs(m) * phi) if m>=0 else sympy.sin(sympy.abs(m) * phi)

    if mode == 'numer':
        return (left_part * right_part).evalf()
    elif mode == 'analit' or mode == 'custom_var':
        return (left_part * right_part).expand()

# Legendre polinomials

def Legendre(n=0, x=Symbol('x')):
    t = type(x).__name__  
    mode = select_exec_mode(x)

    if mode == 'numer': 
        x_val = x
        x = Symbol('x')
    elif mode == 'custom_var':
        x = Symbol(x)

    if mode != 'undef':
        rat_part = 1/(2**n * sympy.factorial(n))
        diff_part = sympy.diff((x**2 - 1)**n, x, n)
    
    if mode == 'numer':
        return (rat_part.subs(x, x_val) * diff_part.subs(x, x_val)).evalf()
    elif mode == 'analit' or mode == 'custom_var':
        return (rat_part * diff_part).expand()
    else:
        raise TypeError

# generalized Legendre polinomials
# gP = lambda n,m,x: (1 - (x**2))**(sympy.abs(m)/2) * sympy.diff(P(n,x), x, abs(m))
def Generalized_Legendre(n=0, m=0, x=Symbol('x')):
    mode = select_exec_mode(x)

    if mode != 'undef':
        legendre_polinomial = Legendre(n=n, x=x)
        rat_part = (1 - (x**2))**(sympy.abs(m)/2)
        m=int(m)
        diff_part = sympy.diff(legendre_polinomial, x, abs(m))

        return (rat_part * diff_part).expand()
    else:
        return False

# THETA-equation solution
# THETA = lambda l,m,theta: sympy.sqrt(sympy.Rational((2 * l + 1),2.0) \
#                             * sympy.Rational(sympy.factorial(l - sympy.abs(m)),sympy.factorial(l + sympy.abs(m)))) \
#                        * gP(l,m,theta).subs(theta, sympy.cos(theta))

def Theta_Equation(l, m, theta):
    ''' l - orbital quantum number
        m - magnetic quantum number '''
    
    # Prevents integer division and float l,m
    l = int(l)
    m = float(int(m))

    mode = select_exec_mode(theta)

    gL = Generalized_Legendre(l, m, theta).subs(theta, sympy.cos(theta))
    
    if gL:
        rat_part = sympy.sqrt(Rational((2 * l + 1),2) * Rational(factorial(l - sympy.abs(m)),factorial(l + sympy.abs(m))))
        if mode == 'numer':
            return (rat_part * gL).evalf()
        elif mode == 'analit':
            return (rat_part * gL).expand()

    return False

#/*
# Angular part
#*/
#Y(l,m,theta,phi) := THETA(l, m, theta) * PHI(m, phi)$
#

def Angular_Part(l, m, theta, phi):
    ''' Executes Angular part of Shregenger equation in spherical coordinates '''

    THETA = Theta_Equation(l, m, theta)
    PHI = Phi_Equation(m, phi)

    return THETA * PHI


print(Angular_Part(1,1,theta,phi))
#/*
# Angular part plot
#*/
#/*
#plot3d(Y(0,0,theta,phi), [theta,0,2*%pi], [phi,0,2*%pi], [transform_xy, make_transform([theta,phi,r],r*sin(phi)*sin(theta), r*cos(phi)*sin(theta),r*cos(theta))])$
#
#Plot(sqrt(6)/(2*sqrt(4*pi)) * cos(phi) * sin(theta)**2, [phi,0,2*pi,35], [theta,-pi,pi,35], 'mode=spherical; color=zfade4')
#Plot((1/(4*sqrt(2*pi))) * exp(-Z * sqrt(x**2)) * x * (Z/0.5)**(3.0/2.0) * , [x,-2,2])


## Generalized
#plot3d(Y(2,1,theta,phi), [theta,0,2*%pi], [phi,0,2*%pi], [transform_xy, make_transform([theta,phi,r],r*sin(phi)*sin(theta), r*cos(phi)*sin(theta),r*cos(theta))], [grid, 50, 50], [plot_format, gnuplot])$
#*/



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
