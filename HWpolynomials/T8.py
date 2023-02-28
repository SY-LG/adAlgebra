from sympy import Symbol,poly
from tools import rational_root

x=Symbol('x')
f=poly(4*x**4-7*x**2-5*x-1)
rational_root(f)