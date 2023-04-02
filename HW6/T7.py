import numpy as np

A=[np.array([[1,-1,3],[-2,2,2],[-1,-1,5]]),np.identity(3)]
B=[np.identity(4),np.array([[1,-1,3,0],[-2,2,2,0],[-1,-1,5,0],[0,1,-1,1]])]
C=np.array([[2,3,1,0],[2,0,1,1],[0,-1,-2,3]])

m,n=A[0].shape
p,q=B[0].shape

G=np.zeros((q*m,p*n))
for Ai,Bi in zip(A,B):
	G+=np.kron(Bi.T,Ai)
c=np.reshape(C.T,(m*q,1))

print("solving: Gx=c")
print("G:",G,sep='\n')
print("c:",c,sep='\n')
x=np.linalg.solve(G,c)

print("solve finish.")
X=np.reshape(x,(q,m)).T
print("X:",X,sep='\n')

print("validating")
result=np.zeros(C.shape)
for Ai,Bi in zip(A,B):
	result+=np.dot(np.dot(Ai,X),Bi)
result-=C
print("result(delta):",result,sep='\n')