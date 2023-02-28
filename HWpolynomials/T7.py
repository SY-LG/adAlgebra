from sympy import Symbol,poly
from tools import sturm

x=Symbol('x')
f=poly(x**5-4*x+2)
print(sturm(f))