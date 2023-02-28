from sympy import Symbol,poly,div

x = Symbol('x')
f=poly(2*x**4-3*x**3+4*x**2)
g=poly(x**2-3*x+1)
q,r=div(f,g)
print(f'q={q.as_expr()}, r={r.as_expr()}')
