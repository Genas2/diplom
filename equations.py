#!/usr/bin/python2

import sympy
from sympy import Symbol, factorial, Rational

class Equations:
    # Declaring variables
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    n = Symbol('n')
    l = Symbol('l')
    m = Symbol('m')
    r = Symbol('r')
    R = Symbol('R')
    phi = Symbol('phi')
    theta = Symbol('theta')
    #Z = Symbol('Z')
    #a0 = Symbol('a0')
    Z = 1
    a0 = 0.52917706
    R = 2.4
    
    types = {}
    __equations__ = []

    def __init__(self):
        # Setting initial quantum numbers values
        self.n_val = 1
        self.l_val = 0
        self.m_val = 0

        #List of available equations
        self.types['Angular part'] = self.Angular_Part
        self.types['Phi equation'] = self.Phi_Equation
        self.types['Ion wavefunction positive'] = self.Wave_Function_Ion_Positive
        self.types['Ion wavefunction negative'] = self.Wave_Function_Ion_Negative

        self.__equations__.append(('Angular part', self.Angular_Part, 'spherical'))
        self.__equations__.append(('Phi Equation', self.Phi_Equation, 'spherical'))
        self.__equations__.append(('Radial part', self.Radial_Part, 'spherical'))
        self.__equations__.append(('Ion wavefunction positive', self.Wave_Function_Ion_Positive, 'cartesian'))
        self.__equations__.append(('Ion wavefunction negative', self.Wave_Function_Ion_Negative, 'cartesian'))
        self.__equations__.append(('Square ion wavefunction positive', self.Square_Wave_Function_Ion_Positive, 'cartesian'))
        self.__equations__.append(('Square ion wavefunction negative', self.Square_Wave_Function_Ion_Negative, 'cartesian'))

    def select_exec_mode(self, x=x):
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
    
    def Phi_Equation(self, m='undef', phi=phi): 
        if m == 'undef':
            m = self.m_val

        m = int(m)
        mode = self.select_exec_mode(phi)
    
        if mode == 'numer': 
            phi_val = self.phi
            self.phi = Symbol('phi')
        elif mode == 'custom_var':
            self.phi = Symbol(phi)
        elif mode == 'undef':
            return False
    
        left_part  = 1/sympy.sqrt(2 * sympy.pi)
        right_part = sympy.cos(sympy.abs(m) * self.phi) if m>=0 else sympy.sin(sympy.abs(m) * self.phi)
    
        if mode == 'numer':
            return (left_part * right_part).evalf()
        elif mode == 'analit' or mode == 'custom_var':
            return (left_part * right_part).expand()
    
    # Legendre polinomials
    
    def Legendre(self, n=0, x=x):
        t = type(x).__name__  
        mode = self.select_exec_mode(x)
    
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
    def Generalized_Legendre(self, n=0, m=0, x=x):
        mode = self.select_exec_mode(x)
    
        if mode != 'undef':
            legendre_polinomial = self.Legendre(n=n, x=x)
            rat_part = (1 - (x**2))**(sympy.abs(m)/2.0)
            diff_part = sympy.diff(legendre_polinomial, x, abs(m))
    
            return (rat_part * diff_part).expand()
        else:
            return False
    
    # THETA-equation solution
    # THETA = lambda l,m,theta: sympy.sqrt(sympy.Rational((2 * l + 1),2.0) \
    #                             * sympy.Rational(sympy.factorial(l - sympy.abs(m)),sympy.factorial(l + sympy.abs(m)))) \
    #                        * gP(l,m,theta).subs(theta, sympy.cos(theta))
    
    def Theta_Equation(self, l='undef', m='undef', theta=theta):
        ''' l - orbital quantum number
            m - magnetic quantum number '''
        
        if l == 'undef':
            l = self.l_val
        if m == 'undef':
            m = self.m_val

        # Prevents integer division and float l,m
        l = int(l)
        m = int(m)
    
        mode = self.select_exec_mode(theta)
    
        gL = self.Generalized_Legendre(l, m, theta).subs(theta, sympy.cos(theta))
        
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
    
    def Angular_Part(self, l='undef', m='undef', theta=theta, phi=phi):
        ''' Generates Angular part of Shregenger equation in spherical coordinates '''
    
        if l == 'undef':
            l = self.l_val
        if m == 'undef':
            m = self.m_val

        THETA = self.Theta_Equation(l, m, theta)
        PHI = self.Phi_Equation(m, phi)
    
        return THETA * PHI
    
    #/*
    # Laguerre polynomials
    #*/
    #L(n,r) := exp(r) * diff(r^n * exp(-r), r, n)$
    def Laguerre(self, n=1, r=Symbol('r')):
        return sympy.radsimp(sympy.exp(r) * sympy.diff(r**n * sympy.exp(-r), r, n))
    
    #/*
    # generalized Legendre polinomials
    #*/
    #aL(n,k,r) := diff(L(n, r), r, k)$
    def Generalized_Laguerre(self, n=1, k=0, r=Symbol('r')):
        return sympy.diff(self.Laguerre(n, r), r, k)
    
    #/*
    # Radial part
    #*/
    #R(n,l,r) := -(((2 * Z)/(n * a0))^3 * (n - l - 1)!/(2 * n * ((n+l)!)^3)) ^ (1/2) * exp(-(Z * r)/(a0 * n)) * at(aL(n+l, 2*l+1, t), t=(2*Z*r)/(a0*n))$
    def Radial_Part(self, n='undef', l='undef', r=Symbol('r')):
        ''' Generates radial part of Shregenger equation '''

        if l == 'undef':
            l = self.l_val
        if n == 'undef':
            n = self.n_val

        #self.a0 = Symbol('a0')

        laguerre_part = self.Generalized_Laguerre(n + l, 2*l + 1, r).subs(r, (2.0 * self.Z * r)/(self.a0 * n))
        left_part = -(((2.0 * self.Z)/(n * self.a0))**3 * Rational(factorial(n - l -1), 2*n * factorial(n+l)**3))**Rational(1,2) 
        exp_part = sympy.exp(-(self.Z * r)/(self.a0 * n)) * ((2 * r)/ (n * self.a0)) ** l
    
        return left_part * exp_part * laguerre_part
    
    #/*
    # wavefunction
    #*/
    #PSI(n,l,m,r,theta,phi) := R(n, l, r) * Y(l, m, theta, phi)$
    #
    def Wave_Function(self, n='undef', l='undef', m='undef', r=Symbol('r'), theta=Symbol('theta'), phi=Symbol('phi')):
        ''' Generates wave function for Hydrogen atom '''
        
        if n == 'undef':
            n = self.n_val
        if l == 'undef':
            l = self.l_val
        if m == 'undef':
            m = self.m_val

        rad_part = self.Radial_Part(n=n, l=l)
        ang_part = self.Angular_Part(l=l, m=m)

        return rad_part * ang_part
    
#######################################################################################################################################
# MOLECULAR ION

# Draft
#S = exp(-R) * (1+ R + R**2/3)
#J = exp(-2*R) * (1 + 1/R)
#K = exp(-R) * (1/R - 2 * R/3)

#Epl = Rational(-1,2) + (J+K)/(1+S)

    def Wave_Function_Ion_Positive(self, R='undef'):
        if R == 'undef':
            R = self.R 
        return sympy.exp(-sympy.sqrt(self.x**2 + self.y**2)/2) + sympy.exp(-sympy.sqrt((self.x-R)**2 + self.y**2)/2)

    def Wave_Function_Ion_Negative(self, R='undef'):
        if R == 'undef':
            R = self.R 
        return sympy.exp(-sympy.sqrt(self.x**2 + self.y**2)/2) - sympy.exp(-sympy.sqrt((self.x-R)**2 + self.y**2)/2)

    def Square_Wave_Function_Ion_Positive(self, R='undef'):
        if R == 'undef':
            R = self.R 
        return sympy.abs(self.Wave_Function_Ion_Positive(R))**2

    def Square_Wave_Function_Ion_Negative(self, R='undef'):
        if R == 'undef':
            R = self.R 
        return sympy.abs(self.Wave_Function_Ion_Negative(R))**2

        #Rgl = 2.4
        #x,y = symbols('xy')
        #psiP = exp(-sqrt(x**2 + y**2)/2) + exp(-sqrt((x-Rgl)**2 + y**2)/2)
        #psiM = exp(-sqrt(x**2 + y**2)/2) - exp(-sqrt((x-Rgl)**2 + y**2)/2)
        #x,y = symbols('xy')
        #Plot(psiP, [x,-4,6], [y,-4,4])
        #Plot(psiP**2, [x,-4,6], [y,-4,4])
        #Plot(psiM, [x,-4,6], [y,-4,4])
        #Plot(psiM**2, [x,-4,6], [y,-4,4])

#######################################################################################################################################
# Temp
#R(n,l,r) := -(((2 * Z) / (n * a0))^3 * (factorial((n - l - 1)) / (2 * n * (factorial(n + l))^3)))^(1/2) * exp(-(Z * r)/(a0 * n)) * aL(n+l, 2*l+1, t)$
#kill(t)$ Rez:R(2,0,r)$ t:(2*Z*r)/(a0*n)$ ratsimp(ev(Rez));
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
    
#
