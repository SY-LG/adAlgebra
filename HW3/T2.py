from sympy import symbols,Matrix
import numpy as np

a,b,c,d,r,s,t=symbols('a b c d r s t')

A=[	[a,r,r,r],
	[a,b,s,s],
	[a,b,c,t],
	[a,b,c,d] ]

A=np.array(A)

def lu(A):
	n=A.shape[0]
	for k in range(n-1):
		if A[k,k]==0:
			raise RuntimeError
		A[k+1:,k]=A[k+1:,k]/A[k,k]
		A[k+1:,k+1:]=A[k+1:,k+1:]-A[k+1:,k]*A[k,k+1:]
		print(A)
	return A

print(lu(A))