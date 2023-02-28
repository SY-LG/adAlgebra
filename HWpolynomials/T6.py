from sympy import Symbol,poly,div

x=Symbol('x')
f=poly(x**5+7*x**4+19*x**3+26*x**2+20*x+8)
g=poly(x+2)

root_cnt=0
while True:
	q,r=div(f,g)
	if r!=0:
		break
	f=q
	root_cnt+=1

print(f'f=({g.as_expr()})**{root_cnt}*({f.as_expr()})')