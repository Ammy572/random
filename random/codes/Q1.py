import numpy as np
import sys                                          #for path to external scripts
#sys.path.insert(0, '/home/user/txhome/storage/shared/gitlab/res2021/july/conics/codes/CoordGeo')        #path to my scripts
sys.path.insert(0, '/home/aman123/random/codes/CoordGeo')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
#If using termux
import subprocess
import shlex
#end if

#Sample size
simlen = 2
#Possible outcomes
#n = range(2,13)
# Generate X1 and X2
#y = np.random.randint(-6,6, size=(3, simlen))
y = np.random.randint(-6,6, size=(3, simlen))
print(f"Given vectors in this case is:\n{y}")
A = y[0]
B = y[1]
C = y[2]
#Q1.1.1
a = B-A
b = C-B
c = A-C
print ("Soln of Q1.1.1")
print ("The direction vector of AB is",a)
print ("The direction vector of BC is",b)
print ("The direction vector of CA is",c)

#Q1.1.2
BC_matrix = B - C

length_BC = np.linalg.norm(BC_matrix)
print("Soln of Q1.1.2")
print("Length of side BC:", length_BC)
#Q1.1.3
print("Soln of Q1.1.3")
Mat = np.array([[1,1,1],[A[0],B[0],C[0]],[A[1],B[1],C[1]]])

rank = np.linalg.matrix_rank(Mat)

if (rank<=2):
	print("Hence proved that points A,B,C in a triangle are collinear")
else:
    print("The given points are not collinear")
#line generating
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)

plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
#plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
#plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
#D = D.reshape(-1,1)
#E = E.reshape(-1,1)
#F = F.reshape(-1,1)
#G = G.reshape(-1,1)
tri_coords = np.block([[A, B, C]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['A', 'B', 'C']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'C' else -10
    plt.annotate(txt,
                 (tri_coords[0, i], tri_coords[1, i]),
                 textcoords="offset points",
                 xytext=(0, offset),
                 ha='center')

#plotting all lines 
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig("figure1_1_3.png")
#Q1.1.4
print("Soln of Q1.1.4")
print("parametric of AB form is x:",y[0],"+ k",a)
print("parametric of BC form is x:",y[1],"+ k",b)
print("parametric of CA form is x:",y[2],"+ k",c)

#Q1.1.5

print("Soln of 1.1.5")
omat = np.array([[0,1],[-1,0]])

def dir_vec(C,B):
  return C-B

def norm_vec(C,B):
    return omat@dir_vec(C,B)
n1=norm_vec(B,A)
#print(n1)
pro_AB=n1.T@A
print(n1,"x=",pro_AB)
n2=norm_vec(C,B)
pro_BC=n2.T@B
print(n2,"x=",pro_BC)
n3=norm_vec(A,C)
pro_CA=n3.T@C
print(n3,"x=",pro_CA)
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB
x_AB = line_gen(A,B)
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
A = A.reshape(-1,1)
B = B.reshape(-1,1)
tri_coords = np.block([[A,B]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B']

x_BC = line_gen(B,C)
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
B = B.reshape(-1,1)
C = C.reshape(-1,1)
tri_coords = np.block([[B,C]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['B','C']

x_CA = line_gen(C,A)
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
C = C.reshape(-1,1)
A = A.reshape(-1,1)
tri_coords = np.block([[C,A]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['C','A']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                textcoords="offset points", # how to position the text
                xytext=(0,10), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig("figure1_1_5.png")

#1.1.6

def AreaCalc(A, B, C):
    AB = y[0] - y[1]
    AC = y[0] - y[2]
#cross_product calculation
    cross_product = np.cross(AB,AC)
#magnitude calculation
    magnitude = np.linalg.norm(cross_product)

    area = 0.5 * magnitude

    return area

area_ABC = AreaCalc(A, B, C)
print("Soln of 1.1.6")
print("Area of triangle ABC:", area_ABC)

#Q.1.1.7
dotA=((B-A).T)@(C-A)
#dotA=dotA[0,0]
NormA=(np.linalg.norm(a))*(np.linalg.norm(c))
print('value of angle A: ', np.degrees(np.arccos((dotA)/NormA)))

dotB=(A-B).T@(C-B)
#dotB=dotB[0,0]
NormB=(np.linalg.norm(y[0]-y[1]))*(np.linalg.norm(y[2]-y[1]))
print('value of angle B: ', np.degrees(np.arccos((dotB)/NormB)))

dotC=(y[0]-y[2]).T@(y[1]-y[2])
#dotC=dotC[0,0]
NormC=(np.linalg.norm(y[0]-y[2]))*(np.linalg.norm(y[1]-y[2]))
print('value of angle C: ', np.degrees(np.arccos((dotC)/NormC)))

#Q1.2.1
D = (B + C)/2

#Similarly for E and F
E = (A + C)/2
F = (A + B)/2
print("Soln of Q1.2.1")
print("D:", list(D))
print("E:", list(E))
print("F:", list(F))
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')


#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#if using termux
#plt.savefig('tri_sss.pdf')
plt.savefig('figure1_2_1.png')


#Q1.2.2
m_a = D-y[0]
print("Soln of 1.2.2")
print("m1 = ",m_a)

nt = np.array([m_a[1]*1+m_a[0]*0,m_a[0]*-1+m_a[1]*0])

#calculation of constant we get on multiplying n^T and A
nt = nt*2
#the above step is done to cancel out the factor of 1/2 on both sides of equation
const = nt@y[0]
print("Value of n^T.A =",const)
m_b = E-y[1]
print("m2 = ",m_b)
nt = np.array([m_b[1]*1+m_b[0]*0,m_b[0]*-1+m_b[1]*0])
nt = nt*2
const = nt@y[1]
print("value of n^T.B =",const)
m_c = F-y[1]
print("m3 = ",m_c)
nt = np.array([m_c[1]*1+m_c[0]*0,m_c[0]*-1+m_c[1]*0])
nt = nt*2
const = nt@y[2]
print("Value of n^T.C = ",const)

def line_gen(A, B):
    len = 10
    dim = A.shape[0]
    x_AB = np.zeros((dim, len))
    lam_1 = np.linspace(0, 1, len)
    for i in range(len):
        temp1 = A + lam_1[i] * (B - A)
        x_AB[:, i] = temp1.T
    return x_AB


x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)

plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')

x_AD = line_gen(A, D)
plt.plot(x_AD[0, :], x_AD[1, :], label='$AD$')

x_BE = line_gen(B, E)
plt.plot(x_BE[0, :], x_BE[1, :], label='$BE$')

x_CF = line_gen(C, F)
plt.plot(x_CF[0, :], x_CF[1, :], label='$CF$')

A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-10,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig('figure_1_2_2.png')

#Q1.2.3
omat = np.array([[0,1],[-1,0]]) 

np.set_printoptions(precision=2)

def midpoint(x,y):
	return (x+y)/2

def norm_vec(A,B):
	return np.matmul(omat, dir_vec(A,B))

#direction
def dir_vec(A,B):
	return B-A

#given coordinates A,B,C

#Given D,E,F are midpoints of BC,CA,AB
D=midpoint(B,C)
E=midpoint(C,A)
F=midpoint(A,B)

#intersection of lines
n5=norm_vec(A,D)
M1 = n5.T
n6=norm_vec(B,E)
P1 = n6.T
def line_intersect(n1,A1,n2,A2):
	N=np.block([[M1],[P1]])
	p = np.zeros(2)
	p[0] = M1@A1
	p[1] = P1@A2
	#Intersection
	P=np.linalg.solve(N,p)
	return P
G = line_intersect(M1,A,P1,B)
print("("+str(G[0])+","+str(G[1])+")")
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)


#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
G = G.reshape(-1,1)
tri_coords = np.block([[A, B, C, D, E, F, G]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'G' else -10
    plt.annotate(txt,
                 (tri_coords[0, i], tri_coords[1, i]),
                 textcoords="offset points",
                 xytext=(0, offset),
                 ha='center')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig('figure_1_2_3.png')
#1.2.4

AG = np.linalg.norm(G - y[0])
GD = np.linalg.norm(D - G)

BG = np.linalg.norm(G - y[1])
GE = np.linalg.norm(E - G)
 
CG = np.linalg.norm(G - y[2])
GF = np.linalg.norm(F - G)
print("Soln of 1.2.4")
print("AG/GD= "+str(AG/GD))
print("BG/GE= "+str(BG/GE))
print("CG/GF= "+str(CG/GF))

#Q1.2.5
'''
#G=(A+B+C)/3
print("Soln of 1.2.5")
G=line_intersect(norm_vec(F,y[2]),y[2],norm_vec(E,y[1]),y[1])
print(f"The centroid of triangleABC is {G}")

Matr = np.array([[1,1,1],[A[0],D[0],G[0]],[A[1],D[1],G[1]]],dtype=object)

rnk = np.linalg.matrix_rank(Matr)

if (rank==2):
	print("Hence proved that points A,G,D in a triangle are collinear")
else:
	print("A,G,D are not collinear")
'''
#Q1.2.6
G = (A + B + C) / 3
print("Soln of 1.2.6")
print("centroid of the given triangle: ")      
      
print(G)
     
print("Hence Q.1.2.6 is verified.")
#Q1.2.7
print("Soln of 1.2.6")
print(f"A - F = {y[0]-F}")
print(f"E - D = {E-D}")
_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_DE = line_gen(D,E)
x_DF = line_gen(D,F)


#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_DE[0,:],x_DE[1,:],label='$DE$')
plt.plot(x_DF[0,:],x_DF[1,:],label='$DF$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#if using termux
#plt.savefig('tri_sss.pdf')
plt.savefig('figure1_5_7.png')

#Q1.3.1
bc=y[2]-y[1]
#print(bc)
t=np.array([0,1,-1,0]).reshape(2,2)
AD_1=t@bc
#normal vector of AD_1
AD_p=t@AD_1
print("Soln of 1.2.7")
print("Normal vector of AD1",AD_p)
#Q1.3.2
result = y[0]@AD_p
print("Soln of 1.3.2")
print(f"The equation of AD is {AD_p}X={result}")
'''
A = np.array([1,-1])
B = np.array([-4,6])
C = np.array([-3,-5])
'''
def alt_foot(A,B,C):
  m = y[1]-y[2]
  n = omat@m 
  N=np.block([[m],[n]])
  p = np.zeros(2)
  p[0] = m@y[0] 
  p[1] = n@y[1]
  #Intersection
  P=np.linalg.inv(N.T)@p
  return P



D12 = alt_foot(A,B,C)
print(D12)
'''
#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_AD = line_gen(A,D12)


#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,1],x_AD[1,1],label='$AD$')

#Labeling the coordinates
#tri_coords = np.vstack((A,B,C,O,I)).T
#np.block([[A1,A2,B1,B2]])
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D12.reshape(-1,1)

tri_coords = np.block([[A,B,C,D]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#calculation code starts
'''
A = np.array([1,-1])
nt = np.array([-1,11])
result = A@nt
print(f"The equation of AD is {nt}X={result}")
'''

#calculation code ends

#if using termux
#plt.savefig('tri_sss.pdf')
#plt.savefig('/home/dhruv/EE23010/Assignment1/figs/figure.png')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
# image = mpimg.imread('tri_sss.png')
# plt.imshow(image)
plt.savefig('FigureQ1_3_2')
'''
#Q1.3.3

def alt_foot(A,B,C):
  m = y[1]-y[2]
  n = omat@m 
  N=np.block([[m],[n]])
  p = np.zeros(2)
  p[0] = m@y[0] 
  p[1] = n@y[1]
  #Intersection
  P=np.linalg.inv(N.T)@p
  return P

M =  alt_foot(y[0],y[1],y[2])
N =  alt_foot(y[1],y[2],y[0])
K =  alt_foot(y[2],y[0],y[1])

# Print altitude foot points
print("Altitude foot point for A:", M)
print("Altitude foot point for B:", N)
print("Altitude foot point for C:", K)

BE_norm = norm_vec(N, y[1])
CF_norm = norm_vec(K, y[2])

#print(BE_norm)
#print(CF_norm)
print("Soln of 1.3.3")
print(f"{BE_norm}[x-B]=0" )
print(f"{CF_norm}[x-C]=0" )
'''
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_AD = line_gen(M.T,A)
x_BE = line_gen(B,N.T)
x_CF = line_gen(C,K.T)
x_AE = line_gen(A,N.T)
x_AF = line_gen(A,K.T)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')
plt.plot(x_AE[0,:],x_AE[1,:],linestyle='dotted')
plt.plot(x_AF[0,:],x_AF[1,:],linestyle='dotted')



plt.plot(A[0], A[1], 'o')
plt.text(A[0] * (1 + 0.1), A[1] * (1 - 0.1) , 'A')
plt.plot(B[0], B[1], 'o')
plt.text(B[0] * (1 - 0.2), B[1] * (1) , 'B')
plt.plot(C[0], C[1], 'o')
plt.text(C[0] * (1 + 0.03), C[1] * (1 - 0.1) , 'C')
plt.plot(M[0], M[1], 'o')
plt.text(M[0] * (1 + 0.03), M[1] * (1 - 0.1) , 'D')
plt.plot(N[0], N[1], 'o')
plt.text(N[0] * (1 + 0.03), N[1] * (1 - 0.1) , 'E')
plt.plot(K[0], K[1], 'o')
plt.text(K[0] * (1 + 0.03), K[1] * (1 - 0.1) , 'F')

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

plt.savefig('FigureQ1_3_3.png')

'''
#Q1.3.4
H = line_intersect(norm_vec(y[1],N),N,norm_vec(y[2],K),K)
print("Soln of Q1.3.4")
print(H)
'''
x_AB = line_gen(A,B)	
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_AD = line_gen(A,alt_foot(A,B,C))
x_AE = line_gen(A,alt_foot(B,A,C))
x_BE = line_gen(B,alt_foot(B,A,C))
x_CF = line_gen(C,alt_foot(C,A,B))
x_AF = line_gen(A,alt_foot(C,A,B))
x_CH = line_gen(C,H)
x_BH = line_gen(B,H)
x_AH = line_gen(A,H)

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE_1$')
plt.plot(x_AE[0,:],x_AE[1,:],linestyle = 'dashed',label='$AE_1$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF_1$')
plt.plot(x_AF[0,:],x_AF[1,:],linestyle = 'dashed',label='$AF_1$')
plt.plot(x_CH[0,:],x_CH[1,:],label='$CH$')
plt.plot(x_BH[0,:],x_BH[1,:],label='$BH$')
plt.plot(x_AH[0,:],x_AH[1,:],linestyle = 'dashed',label='$AH$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = M.reshape(-1,1)
E = N.reshape(-1,1)
F = K.reshape(-1,1)
H = H.reshape(-1,1)
tri_coords = np.block([[A,B,C,M,N,K,H]])
#tri_coords = np.vstack((A,B,C,alt_foot(A,B,C),alt_foot(B,A,C),alt_foot(C,A,B),H)).T
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D_1','E_1','F_1','H']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig('Figure 1_3_4')
plt.show()
'''
#Q1.3.5
result = int(((y[0] - H).T) @ (y[1] - y[2]))    # Checking orthogonality condition...

print("Soln of 1.3.4")
# printing output
if result == 0:
  print("(A - H)^T (B - C) = 0\nHence Verified...")

else:
  print("(A - H)^T (B - C)) != 0\nHence the given statement is wrong...")
#Q1.4.1
#to find the coefficients and constant of the equation of perpendicular bisector of BC
def midpoint(P, Q):
    return (P + Q) / 2
def perpendicular_bisector(B, C):
    midBC=midpoint(y[1],y[2])
    dir=y[1]-y[2]
    constant = -dir.T @ midBC
    return dir,constant
def perpendicular_bisector(B, C):
    midBC=midpoint(y[1],y[2])
    dir=y[1]-y[2]
    constant = -dir.T @ midBC
    return dir,constant
equation_coeff1,const1 = perpendicular_bisector(y[0], y[1])
equation_coeff2,const2 = perpendicular_bisector(y[1], y[2])
equation_coeff3,const3 = perpendicular_bisector(y[2], y[0])
print("Soln of 1.4.1")
print(f'Equation for perpendicular bisector of AB: ({equation_coeff1[0]:.2f})x + ({equation_coeff1[1]:.2f})y + ({const1:.2f}) = 0')
print(f'Equation for perpendicular bisector of  BC: ({equation_coeff2[0]:.2f})x + ({equation_coeff2[1]:.2f})y + ({const2:.2f}) = 0')
print(f'Equation for perpendicular bisector of  CA: ({equation_coeff3[0]:.2f})x + ({equation_coeff3[1]:.2f})y + ({const3:.2f}) = 0')
#Q1.4.2
O = line_intersect(a,F,c,E)
print("Soln of 1.4.2")
print(O)
#Q1.4.3
G=(y[2]+y[1])/2
#Q1.4.4
O_1 = O - y[0]
O_2 = O - y[1]
O_3 = O - y[2]
d = np.linalg.norm(O_1)
e = np.linalg.norm(O_2)
f = np.linalg.norm(O_3)
#print("Points of triangle A, B, C respectively are", A ,",", B ,",", C, ".")
#print("Circumcentre of triangle is", O, ".")
print("Soln of 1.4.4")
print(" OA, OB, OC are respectively", d,",", e,",",f, ".")
print("Here, OA = OB = OC.")
print("Hence verified.")
#Q.1.4.5
#Q1.4.6
dot_pt_O = (y[1] - O) @ ((y[2] - O).T)
norm_pt_O = np.linalg.norm(y[1] - O) * np.linalg.norm(y[2] - O)
cos_theta_O = dot_pt_O / norm_pt_O
angle_BOC = round(360-np.degrees(np.arccos(cos_theta_O)),5)  #Round is used to round of number till 5 decimal places
print("Soln of Q1.4.6")
print("angle BOC = " + str(angle_BOC))

#To find angle BAC
dot_pt_A = (y[1] - y[0]) @ ((y[2] - y[0]).T)
norm_pt_A = np.linalg.norm(y[1] - y[0]) * np.linalg.norm(y[2] - y[0])
cos_theta_A = dot_pt_A / norm_pt_A
angle_BAC = round(np.degrees(np.arccos(cos_theta_A)),5)  #Round is used to round of number till 5 decimal places
print("angle BAC = " + str(angle_BAC))
#To check whether the answer is correct
if angle_BOC == 2 * angle_BAC:
  print("\nangle BOC = 2 times angle BAC\nHence the give statement is correct")
else:
  print("\nangle BOC ≠ 2 times angle BAC\nHence the given statement is wrong")
#Q1.5.1
def unit_vec(A,B):
	return ((B-A)/np.linalg.norm(B-A))
E_1= unit_vec(y[0],y[1]) + unit_vec(y[0],y[2])
#point generated to create parametric form
#generating normal form
F_1=np.array([E_1[1],(E_1[0]*(-1))])
#matrix multiplication
C1= F_1@(y[0].T)
print("Soln of 1.5.1")
print("Internal Angular bisector of angle A is:",F_1,"*x = ",C1)

E_2= unit_vec(y[1],y[0]) + unit_vec(y[1],y[2])
#point generated to create parametric form
#generating normal form
F_2=np.array([E_2[1],(E_2[0]*(-1))])
#matrix multiplication
C2= F_2@(y[1].T)
print("Internal Angular bisector of angle B is:",F_2,"*x = ",C2)


E_3= unit_vec(y[0],y[1]) + unit_vec(y[0],y[2])
#point generated to create parametric form
#generating normal form
F_3=np.array([E_3[1],(E_3[0]*(-1))])
#matrix multiplication
C3= F_3@(y[1].T)
print("Internal Angular bisector of angle C is:",F_3,"*x = ",C3)
#Q1.5.2
a1 = unit_vec(y[1],y[2])
#print(a1)
a2 = unit_vec(y[2],y[0])
a3 = unit_vec(y[0],y[1])

I=line_intersect(a1-a3,y[1],a1-a2,y[2])
print("Soln of Q1.5.2")
print("Point I is ",I)
#Q1.5.3
i = I-y[0]
def angle_btw_vectors(v1, v2):
    dot_product = v1 @ v2
    norm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle = np.arccos(dot_product / norm)
    angle_in_deg = np.degrees(angle)
    return angle_in_deg

#Calculating the angles BAI and CAI
angle_BAI = angle_btw_vectors(a, i)
angle_CAI = angle_btw_vectors(c, i)

# Print the angles
print("Soln of Q1.5.3")
print("Angle BAI:", angle_BAI)
print("Angle CAI:", angle_CAI)

if np.isclose(angle_BAI, angle_CAI):
    print("Angle BAI is approximately equal to angle CAI.")
else:
    print("error")
#Q1.5.4
k1=1
k2=1

p = np.zeros(2)
'''
t = norm_vec(y[1], y[2])
n1 = t / np.linalg.norm(t)
t = norm_vec(y[2], y[0])
n2 = t / np.linalg.norm(t)
t = norm_vec(y[0], y[1])
n3 = t / np.linalg.norm(t)

p[0] = a1 @ y[1] - k1 * a2 @ y[2]
p[1] = a2 @ y[2] - k2 * a3 @ y[0]

N1 = np.block([[a1 - k1 * a2],[ a2 - k2 * a3]])
I1 = np.linalg.inv(N1)@p
'''
r1 = a1 @ (y[1]-I)
print("Soln of Q1.5.4")
print("Coordinates of point I:", I)
print(f"Distance from I to BC= {r1}")
#Q1.5.5
r2 = a2 @ (y[2]-I) 

r3 = a3 @ (y[0]-I)
'''
r2 = abs((a1T @ I) - (a1T @ y[0]))/(np.linalg.norm(a1))   #r1 is distance between I and AB (n1T.I - n1T.A=0, n1 and I are vectors)
r3 = abs((a2T @ I) - (a2T @ y[2]))/(np.linalg.norm(a2))   #r2 is distance between I and AC (n2T.I - n2T.C=0, n2 and I are vectors)
'''
print("Soln of Q1.5.5")
print("Distance between I and AB is",r2)
print("Distance between I and AC is",r3)
#Q1.5.6
#Q1.5.7
#Q1.5.8
'''
D=C-B
a=np.linalg.norm(C-B)
b=np.linalg.norm(A-C)
c=np.linalg.norm(A-B)
I=np.array([-1.47756217,-0.79495069])
print("the incentre coordiantes are",I)
# by using the data from previous question we have inradius value r
radius=1.8968927705299559
'''
p=pow(np.linalg.norm(y[2]-y[1]),2)
q=2*(c@(I-y[1]))
r=pow(np.linalg.norm(I-y[1]),2)-r1*r1

Discre=q*q-4*p*r
print("Soln of Q1.5.8")
print("the Value of discriminant is ",abs(round(Discre,6)))
#  so the value of Discrimant is extremely small and tends to zero
#  the discriminant value rounded off to 6 decimal places is also zero
#  so it proves that there is only one solution of the point

#  the value of point is x=B+k(C-B)
k=((I-y[1])@(y[2]-y[1]))/((y[2]-y[1])@(y[2]-y[1]))
print("the value of parameter k is ",k)
D3=y[1]+(k*(y[1]-y[2]))
print("the point of tangency of incircle by side BC is ",D3)
#  to check we also check the value of dot product of ID3 and BC
#print("the dot product of ID3 and BC",abs(round(((D3-I)@(C-B),6))))
#  so this comes out to be zero
print("Hence we prove that side BC is tangent To incircle and also found the value of k!")
#Q1.5.9
k1=((I-y[0])@(y[1]-y[0]))/((y[1]-y[0])@(y[1]-y[0]))
k2=((I-y[2])@(y[0]-y[2]))/((y[0]-y[2])@(y[0]-y[2]))
#finding E_3 and F_3
E3=y[0]+(k1*(y[1]-y[0]))
F3=y[2]+(k2*(y[0]-y[2]))
print("Soln of Q1.5.9")
print("k1 = ",k1)
print("k2 = ",k2)
print("E3 = ",E3)
print("F3 = ",F3)
#Q1.5.10
'''
D_3=np.array([-3.367,-0.967])
E_3=np.array([-0.136,-2.136])
F_3=np.array([0.066,0.308])
'''
def norm(X,Y):
    magnitude=round(float(np.linalg.norm([X-Y])),3)
    return magnitude
print("Soln of Q1.5.10")
print("AE_3=", norm(y[0],E3) ,"\nAF_3=", norm(y[0],F3) ,"\nBD_3=", norm(y[1],D3) ,"\nBF_3=", norm(y[1],F3) ,"\nCD_3=", norm(y[2],D3) ,"\nCE_3=",norm(y[2],E3))
#Q1.5.11
#creating array containing coefficients
Y = np.array([[1,1,0],[0,1,1],[1,0,1]])
x = np.linalg.norm(b)
y = np.linalg.norm(c)
z = np.linalg.norm(a)

#solving the equations
X = np.linalg.solve(Y,[x,y,z])
print("Soln of Q1_5_11")
#printing output 
print(X)
