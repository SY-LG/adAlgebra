from sympy import Symbol,poly
from tools import euclid

x=Symbol('x')
f=poly(x**4+3*x**3-x**2-4*x+3)
g=poly(3*x**3+10*x**2+2*x-3)

euclid(f,g)