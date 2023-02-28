from sympy import Symbol,poly,div
from tools import polys_format

x=Symbol('x')
f=poly(x**5+7*x**4+19*x**3+26*x**2+20*x+8)
f0=f.copy()
g=poly(x+2)

root_cnt=0
while True:
	q,r=div(f,g)
	if r!=0:
		break
	f=q
	root_cnt+=1

print(polys_format('{}={}**{}*{}',f0,g,root_cnt,q))