from sympy import Symbol,symbols,poly,div,divisors
import math

# 辗转相除
# r_sig_positive=True: f=qg+r
# else: f=qg-r
def division_algorithm(f,g,r_sig_positive=True,verbose=False):
	sig='+' if r_sig_positive else '-'
	q=[]
	r=[]
	while True:
		q_temp,r_temp=div(f,g)
		if verbose:
			print(f'{str(f.as_expr()):20}=({str(q_temp.as_expr()):40})({str(g.as_expr()):20}){sig}({str(r_temp.as_expr()):20})')
		if r_temp==0:
			break
		if not r_sig_positive:
			r_temp=-r_temp
		q.append(q_temp)
		r.append(r_temp)
		f=g
		g=r_temp
	return q,r


# 两个多项式f,g（无论数域）的最大公因式（首一）为d，存在u,v使d=uf+vg
def euclid(f,g):
	q,r=division_algorithm(f,g)
	q.reverse()
	u_last=0
	v_last=1
	for qi in q:
		u=v_last
		v=u_last-qi*v_last
		u_last=u
		v_last=v
	d=r[-1]
	lc=d.LC()
	return d/lc,u/lc,v/lc

# 实系数多项式的实根范围
def boundary(f):
	def ub(f):
		coeffs=f.all_coeffs()
		if coeffs[0]<0:
			coeffs=[-coeff for coeff in coeffs]
		for k in range(len(coeffs)):
			if coeffs[k]<0:
				break
		else:
			return 0
		b=-min(coeffs)
		return 1+math.pow(b/coeffs[0],1/k)
	# g=f(-x)
	x=list(f.free_symbols)[0]
	y=Symbol('y',exclude=x)
	g=poly(f.subs(x,-y).subs(y,x))
	return -ub(g),ub(f)

# 实系数多项式的实根个数
def sturm(f,bound=None):
	if not bound:
		bound=boundary(f)
	def V(f,c):
		x=list(f.free_symbols)[0]
		df=f.diff()
		d,_,_=euclid(f,df)
		f,_=div(f,d)
		g=f.diff()
		_,r=division_algorithm(f,g,False)
		r=[g]+r
		temp_sig=1 if f.subs(x,c)>0 else -1
		v=0
		for g in r:
			if temp_sig*g.subs(x,c)<0:
				v+=1
				temp_sig*=-1
		return v
	lower_bound,upper_bound=bound
	return V(f,lower_bound)-V(f,upper_bound)

# 整系数多项式的有理根
def rational_root(f):
	roots=[]
	x=list(f.free_symbols)[0]
	a_n=f.all_coeffs()[0]
	a_0=f.all_coeffs()[-1]
	for p in divisors(a_n):
		for q in divisors(a_0):
			for sig in [-1,1]:
				try_root=q/p*sig
				if f.subs(x,try_root)==0:
					roots.append(try_root)
	return roots