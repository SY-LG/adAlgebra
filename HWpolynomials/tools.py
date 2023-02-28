from sympy import Symbol,symbols,poly,div,divisors
import math,numpy,re

def new_print(verbose,indent):
	return lambda x:print('  '*indent+x) if verbose else lambda x:None

def polys_format(string,*args,L=None):
	new_args=[]
	for arg in args:
		if str(type(arg))=="<class 'sympy.polys.polytools.Poly'>":
			arg='('+str(arg.as_expr())+')'
		else:
			arg=str(arg)
		new_args.append(arg)
	return string.format(*new_args)

def eq_print(eq,head='in equation: ',print=print,**kwargs):
	print(head+eq)
	for poly in re.split(r'[+-\\*/=]',eq):
		print(f'{poly}={kwargs[poly].as_expr()}')

def series_print(name,series,print=print):
	output=name
	for item in series:
		output+=' '+str(item.as_expr())+','
	print(output[:-1])

# 辗转相除
# r_sig_positive=True: f=qg+r
# else: f=qg-r
def division_algorithm(f,g,r_sig_positive=True,verbose=True,indent=0):
	print=new_print(verbose,indent)
	print(polys_format('running division algorithm for f={} and g={}',f,g))

	q=[]
	r=[]
	while True:
		q_temp,r_temp=div(f,g)
		if not r_sig_positive:
			r_temp=-r_temp
		if r_sig_positive:
			print(polys_format('{}={}*{}+{}',f,g,q_temp,r_temp))
		else:
			print(polys_format('{}={}*{}-{}',f,g,q_temp,r_temp))
		if r_temp==0:
			break
		q.append(q_temp)
		r.append(r_temp)
		f=g
		g=r_temp
	return q,r


# 两个多项式f,g（无论数域）的最大公因式（首一）为d，存在u,v使d=uf+vg
def euclid(f,g,verbose=True,indent=0):
	print=new_print(verbose,indent)
	print(polys_format('running euclid for f={} and g={}',f,g))

	q,r=division_algorithm(f,g,verbose=verbose,indent=indent+1)
	q.reverse()
	u_last=0
	v_last=1
	for qi in q:
		u=v_last
		v=u_last-qi*v_last
		u_last=u
		v_last=v
	d0=r[-1]
	lc=d0.LC()
	d,u,v=d0/lc,u/lc,v/lc

	eq_print('d=u*f+v*g',head='in euclid: ',print=print,f=f,g=g,u=u,v=v,d=d)
	return d,u,v

# 实系数多项式实根的简单上界
def upper_bound(f,verbose=True,indent=0):
	print=new_print(verbose,indent)
	print(polys_format('calculating upper bound of root of {}',f))

	coeffs=f.all_coeffs()
	if coeffs[0]<0:
		print('a_n<0, reversing f')
		coeffs=[-coeff for coeff in coeffs]
	for k in range(len(coeffs)):
		if coeffs[k]<0:
			break
	else:
		print('forall i, a_i>=0, upper bound is 0')
		return 0
	b=-min(coeffs)
	a_n=coeffs[0]

	print(f'find k={k}, b={b}, a_n={a_n}')
	return 1+math.pow(b/a_n,1/k)

# 实系数多项式的实根范围
def boundary(f,verbose=True,indent=0):
	print=new_print(verbose,indent)
	print(polys_format('calculating root boundary of {}',f))

	# g=f(-x)
	x=list(f.free_symbols)[0]
	y=Symbol('y',exclude=x)
	g=poly(f.subs(x,-y).subs(y,x))
	print(polys_format('f(-x)={}',g))

	lb=-upper_bound(g,verbose=verbose,indent=indent+1)
	ub=upper_bound(f,verbose=verbose,indent=indent+1)

	print(polys_format('root boundary of {} is [{},{}]',f,lb,ub))
	return lb,ub

# sturm序列的变号数
def V(sturm_series,c):
	v=temp_sig=0
	for g in sturm_series:
		symbols=g.free_symbols
		x=list(symbols)[0] if len(symbols) else Symbol('x')
		g_value=g.subs(x,c)
		if temp_sig==0:# start only
			temp_sig=numpy.sign(g_value)
		elif temp_sig*g_value<0:
			v+=1
			temp_sig*=-1
	return v

# 实系数多项式的实根个数
def sturm(f,bound=None,verbose=True,indent=0):
	print=new_print(verbose,indent)
	if not bound:
		print(polys_format('running sturm for f={}, boudary unspecified',f))
		bound=boundary(f,verbose=verbose,indent=indent+1)
	else:
		print(polys_format('running sturm for f={}, boudary specified as {}',f,bound))
	lower_bound,upper_bound=bound
	x=list(f.free_symbols)[0]

	df=f.diff()
	print(polys_format('derivative of f is {}',df))
	d,_,_=euclid(f,df,verbose=verbose,indent=indent+1)
	f,_=div(f,d)
	print(polys_format('new f for sturm calculation: {}',f))

	g=f.diff()
	print(polys_format('derivative of new f is {}',g))
	_,r=division_algorithm(f,g,False,verbose=verbose,indent=indent+1)
	sturm_series=[f,g]+r
	series_print("sturm series:",sturm_series,print=print)

	v_low=V(sturm_series,lower_bound)
	v_up=V(sturm_series,upper_bound)
	root_num=v_low-v_up

	print(f'V({lower_bound})={v_low}, V({upper_bound})={v_up}, root num={root_num}')
	return root_num

# 整系数多项式的有理根
def rational_root(f,verbose=True,indent=0):
	print=new_print(verbose,indent)
	print(polys_format('calculating rational roots of {}',f))

	roots=[]
	x=list(f.free_symbols)[0]
	a_n=f.all_coeffs()[0]
	a_0=f.all_coeffs()[-1]
	divisors_a_n=divisors(a_n)
	divisors_a_0=divisors(a_0)
	print(f'a_n={a_n}, divisors: {divisors_a_n}')
	print(f'a_0={a_0}, divisors: {divisors_a_0}')
	for p in divisors_a_n:
		for q in divisors_a_0:
			for sig in [-1,1]:
				try_root=q/p*sig
				value=f.subs(x,try_root)
				print(f'f({try_root})={value}')
				if value==0:
					roots.append(try_root)

	print(f'in conclusion, rational roots are {roots}')
	return roots