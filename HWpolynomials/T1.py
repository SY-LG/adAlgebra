from sympy import Symbol,poly,div
from tools import eq_print

x = Symbol('x')
f=poly(2*x**4-3*x**3+4*x**2)
g=poly(x**2-3*x+1)
q,r=div(f,g)
eq_print('f=g*q+r',f=f,g=g,q=q,r=r)