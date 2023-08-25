# #Code by GVV Sharma
# #December 7, 2019
# #released under GNU GPL
# #Drawing a triangle given 3 sides

# import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
# image = mpimg.imread('exit-ramp.jpg')
# plt.imshow(image)
# plt.show()

import sys                                          #for path to external scripts
#sys.path.insert(0, '/home/user/txhome/storage/shared/gitlab/res2021/july/conics/codes/CoordGeo')        #path to my scripts
sys.path.insert(0, '/home/aman123/random/code1/CoordGeo')        #path to my scripts
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
from params import omat

#sys.path.insert(0, '/home/user/txhome/storage/shared/gitlab/res2021/july/conics/codes/CoordGeo')        #path to my scripts

#if using termux
#import subprocess
#import shlex
#end if

import numpy as np
import math as mt
import matplotlib.pyplot as plt
#If using termux
#import subprocess
#import shlex
#end if

def norm_(A):
	return(np.linalg.norm(A))
	
#Sample size
simlen = 2
#Possible outcomes
#n = range(0,7)
# Generate X1 and X2
A = np.array([0,0])
B = np.array([0, 0])
C = np.array([0,0])
y = np.random.randint(-6,7, size=(3, simlen))
print(y)
print(omat)
#Given that, 
A = y[0]
B = y[1]
C = y[2]


#A = np.array([1,-1])
#B = np.array([-4, 6])
#C = np.array([-3,-5])

A2X1 = A.reshape(2,1)
B2X1 = B.reshape(2,1)
C2X1 = C.reshape(2,1)

#Problem 1.1.1
m_AB = dir_vec(A,B)
m_BC = dir_vec(B,C)
m_CA = dir_vec(C,A)
print("The direction vector of AB is ",m_AB)
print("The direction vector of BC is ",m_BC)
print("The direction vector of CA is ",m_CA)

#Problem 1.1.2
BC_matrix = m_BC
length_BC = norm_(BC_matrix)
print("Length of side BC:", length_BC)

#Problem 1.1.3
one_mat=np.ones(3)
Mat = np.block([[one_mat],[A2X1,B2X1,C2X1]])
rank = np.linalg.matrix_rank(Mat)
if (rank<=2):
	print("Hence proved that points A,B,C in a triangle are collinear")
else:
	print("The given points are not collinear")

#Problem 1.1.4
print(f"parametric of AB form is x: {A} + k {m_AB}")
print(f"parametric of BC form is x: {B} + k {m_BC}")
print(f"parametric of CA form is x: {C} + k {m_CA}")

#Problem 1.1.5
n_AB = norm_vec(A,B)
n_BC = norm_vec(B,C)
n_CA = norm_vec(C,A)
print(f"The normal form equation of AB is {n_AB.T}x={n_AB.T@B}") 
print(f"The normal form equation of BC is {n_BC.T}x={n_BC.T@C}") 
print(f"The normal form equation of CA is {n_CA.T}x={n_CA.T@A}") 

#Problem 1.1.6
cross_product = np.cross(m_AB,dir_vec(A,C))
magnitude = norm_(cross_product)
Area_ABC = 0.5 * magnitude
print(f"Area of triangle ABC: {Area_ABC}")

#Problem 1.1.7
cosA = (m_AB.T@(dir_vec(A,C)))/(norm_(m_AB)*norm_(dir_vec(A,C)))
cosB = (m_BC.T@(dir_vec(B,A)))/(np.linalg.norm(m_BC)*np.linalg.norm(dir_vec(B,A)))
cosC = (m_CA.T@(dir_vec(C,B)))/(np.linalg.norm(m_CA)*np.linalg.norm(dir_vec(C,B)))
print(f"The cosine of angle A is {cosA}")
print(f"The cosine of angle B is {cosB}")
print(f"The cosine of angle C is {cosC}")

#Problem 1.2.1
D = (B + C)/2
E = (A + C)/2 #Mid-points
F = (A + B)/2
print("D:", list(D))
print("E:", list(E))
print("F:", list(F))

#Problem 1.2.2
n_AD = norm_vec(A,D)
n_BE = norm_vec(B,E)
n_CF = norm_vec(C,F)
print(f"The normal form equation of AD is {list(n_AD.T)}x={n_AD.T@A}") 
print(f"The normal form equation of BE is {list(n_BE.T)}x={n_BE.T@B}") 
print(f"The normal form equation of CF is {list(n_CF.T)}x={n_CF.T@C}") 

#Problem 1.2.3
P = np.block([[n_BE.T],[n_CF.T]])
R = np.block([[n_BE.T@B],[n_CF.T@C]])
G = np.linalg.solve(P,R)
G = G.T
print(G.tolist())

#Problem 1.2.4
print(f"The value of BG/GE is {norm_(dir_vec(B,G))/norm_(dir_vec(G,E))}")
print(f"The value of CG/GF is {norm_(dir_vec(C,G))/norm_(dir_vec(G,F))}")
print(f"The value of AG/GD is {norm_(dir_vec(A,G))/norm_(dir_vec(G,D))}")

#Problem 1.2.5
D2X1=D.reshape(2,1)
G2X1=G.T
Mat2 = np.block([[one_mat],[A2X1,D2X1,G2X1]])
rank = np.linalg.matrix_rank(Mat2)
if (rank==2):
	print("Hence proved that points A,G,D in a triangle are collinear")
else:
	print("Error")
	
#Problem 1.2.6
print(f"The value of G is {G.tolist()}")
print(f"The value of (A+B+C)/3 is {((A+B+C)/3).tolist()}") #Hence proved 

#Problem 1.2.7
print(f"The values of A - F = {list(dir_vec(F,A))} and E - D = {list(dir_vec(D,E))}")

#Problem 1.3.1
D_1 = alt_foot(A,B,C)
print("D1",D_1)
E_1 = alt_foot(B,C,A)
print("E1",E_1)
F_1 = alt_foot(C,A,B)
print("F1",F_1)
n_AD_1 = m_BC #As AD_1 is perp to BC
print(f"The normal vector of AD_1 is {n_AD_1}") 
#We can also cross verify using omat@AD_1

#Problem 1.3.2
print(f"The equation of line AD_1 is {list(n_AD_1.T)}x={n_AD_1.T@A}")

#Problem 1.3.3
n_BE_1=m_CA
n_CF_1=m_AB
print(f"The equation of line BE_1 is {list(n_BE_1.T)}x={n_BE_1.T@B}")
print(f"The equation of line CF_1 is {list(n_CF_1.T)}x={n_CF_1.T@C}")

#Problem 1.3.4
H = line_intersect(norm_vec(B,E_1),E_1,norm_vec(C,F_1),F_1) 
#P=np.block([[n_BE_1],[n_CF_1]])
#R=np.block([[n_BE_1.T@B],[n_CF_1.T@C]])
#H = np.linalg.solve(P,R) 
#H = H.T
#print(H.tolist()) #using solve linalg cross verify
print(H)

#Problem 1.3.5
product = (A - H)@((B - C).T)
print(product)

#Problem 1.4.1
print(f"The perpendicular bisector of AB is {dir_vec(B,A)}x = {dir_vec(B,A).T@(A+B)/2}")
print(f"The perpendicular bisector of BC is {dir_vec(C,B)}x = {dir_vec(C,B).T@(C+B)/2}")
print(f"The perpendicular bisector of CA is {dir_vec(A,C)}x = {dir_vec(A,C).T@(A+C)/2}")

#Problem 1.4.2
O = line_intersect(m_AB,(A+B)/2,dir_vec(A,C),(A+C)/2)
#P=np.block([[dir_vec(B,A)],[dir_vec(A,C)]])
#R=np.block([[dir_vec(B,A).T@(A+B)/2],[dir_vec(A,C).T@(A+C)/2]])
#O = np.linalg.solve(P,R) 
#O = O.T
#print(O.tolist()) #using solve linalg cross verify
print("O",O)
OA1 = dir_vec(O,A)
print("nOA",OA1 )
OB1 = dir_vec(O,B)
print("nOB",OB1)
OC1 = dir_vec(O,C)
print("nOC",OC1)
#Problem 1.4.3
print(f"{(dir_vec(B,A)@O)-dir_vec(B,A).T@(A+B)/2}")
#Hence verified 

#Problem 1.4.4
R = norm_(dir_vec(O,A)) #Circumradius
print(f"OA={norm_(dir_vec(O,A))}")
print(f"OB={norm_(dir_vec(O,B))}")
print(f"OC={norm_(dir_vec(O,C))}")

#Problem 1.4.5 
	#Figures 

#Problem 1.4.6
cosBAC = (m_AB.T@(dir_vec(A,C)))/(norm_(m_AB)*norm_(dir_vec(A,C)))
cosBOC = (dir_vec(O,B).T@(dir_vec(O,C)))/(norm_(dir_vec(O,B))*norm_(dir_vec(O,C)))
Angle_BAC = mt.degrees(mt.acos(cosBAC))
Angle_BOC = mt.degrees(mt.acos(cosBOC))
print(f"The angle BAC is {Angle_BAC}")
print(f"The angle BOC is {Angle_BOC}")

#Problem 1.4.7
const_mat=dir_vec(O,C) #L.H.S
mat_AO = dir_vec(O,A)  #R.H.S
var_mat = [[mat_AO[0],-mat_AO[1]],[mat_AO[1], mat_AO[0]]] #modifying equation such that X Y is a 2X1 rather than 4X4
angles = np.linalg.solve(var_mat,const_mat)
Angle_cos=abs(mt.degrees(mt.acos(angles[0]))) #Both should be same 
Angle_sin=abs(mt.degrees(mt.asin(angles[1])))
print(Angle_cos) #Both should be either equal or sum upto 180.. As sin(180-x) = sin(x)
print(Angle_sin)

#Problem 1.5.1
u_norm_AB = n_AB/norm_(n_AB)
u_norm_BC = n_BC/norm_(n_BC)
u_norm_CA = n_CA/norm_(n_CA)
c_1=n_AB.T@B
c_2=n_BC.T@C
c_3=n_CA.T@A
print(f"The angle bisector of angle A is {(u_norm_AB - u_norm_CA)}x = {(c_1/norm_(n_AB))-(c_3/norm_(n_CA))}")
print(f"The angle bisector of angle B is {(u_norm_AB - u_norm_BC)}x = {(c_1/norm_(n_AB))-(c_2/norm_(n_BC))}")
print(f"The angle bisector of angle C is {(u_norm_BC - u_norm_CA)}x = {(c_2/norm_(n_BC))-(c_3/norm_(n_CA))}")

#Problem 1.5.2
I = line_intersect((u_norm_AB - u_norm_BC),B,(u_norm_BC - u_norm_CA),C)
#P=np.block([[(u_norm_AB - u_norm_BC)],[(u_norm_BC - u_norm_CA)]])
#R=np.block([[(c_1/norm_(n_AB))-(c_2/norm_(n_BC))],[(c_2/norm_(n_BC))-(c_3/norm_(n_CA))]])
#I = np.linalg.solve(P,R) 
#I = I.T
#print(I.tolist()) #using solve linalg cross verify
print(I)

#Problem 1.5.3
cosBAI= (m_AB.T@(dir_vec(A,I)))/(norm_(m_AB)*norm_(dir_vec(A,I)))
cosCAI = (dir_vec(A,C).T@(dir_vec(A,I)))/(norm_(dir_vec(A,C))*norm_(dir_vec(A,I)))
Angle_BAI = mt.degrees(mt.acos(cosBAI))
Angle_CAI = mt.degrees(mt.acos(cosCAI))
print(f"The angle BAI is {Angle_BAI}")
print(f"The angle CAI is {Angle_CAI}")

#Problem 1.5.4
foot_perp_IBC = alt_foot(I,B,C)
pdist_IBC = norm_(dir_vec(foot_perp_IBC,I)) #Perpendiular distance from I to BC 
print(f"Perpendiular distance from I to BC is {pdist_IBC}")

#Problem 1.5.5
foot_perp_IAB = alt_foot(I,A,B)
foot_perp_ICA = alt_foot(I,C,A)
pdist_IAB = norm_(dir_vec(foot_perp_IAB,I)) #Perpendiular distance from I to AB 
pdist_ICA = norm_(dir_vec(foot_perp_ICA,I)) #Perpendiular distance from I to CA 
print(f"Perpendiular distance from I to AB is {pdist_IAB}")
print(f"Perpendiular distance from I to CA is {pdist_ICA}")

#Problem 1.5.6
I,r = icircle(A,B,C)
print(f"The inradius of the circle ABC is {abs(r)}")

#Problem 1.5.7
	#Figures
	
#Problem 1.5.8
print(f"parametric of BC form is x: {B} + k {m_BC}")
#Circle eq is ||x-I||^2 = r^2 
#Substituting x in the circle eq we get the quadratic eq ak^2 + bk + c where,
m = m_BC
a = norm_(m)**2
b = 2*(m@(B-I))
c = norm_(B-I)**2 - r**2 
Det = b**2 - 4*a*c
print(Det) #D=0 so same roots 
k = -b/(2*a)
D_3 = B + k*m_BC
print(f"point D_3 is {D_3}")

#Problem 1.5.9
#E_3
print(f"parametric of CA form is x: {C} + k {m_CA}")
#Circle eq is ||x-I||^2 = r^2 
#Substituting x in the circle eq we get the quadratic eq ak^2 + bk + c where,
m = m_CA
a = norm_(m)**2
b = 2*(m@(C-I))
c = norm_(C-I)**2 - r**2 
Det = b**2 - 4*a*c
print(Det) #D=0 so same roots 
k = -b/(2*a)
E_3 = C + k*m_CA
print(f"point E_3 is {E_3}")

#F_3
print(f"parametric of AB form is x: {A} + k {m_AB}")
#Circle eq is ||x-I||^2 = r^2 
#Substituting x in the circle eq we get the quadratic eq ak^2 + bk + c where,
m = m_AB
a = norm_(m)**2
b = 2*(m@(A-I))
c = norm_(A-I)**2 - r**2 
Det = b**2 - 4*a*c
print(Det) #D=0 so same roots 
k = -b/(2*a)
F_3 = A + k*m_AB
print(f"point F_3 is {F_3}")

#Problem 1.5.10
AE3 = norm_(dir_vec(A,E_3))
AF3 = norm_(dir_vec(A,F_3))
m = AE3
print(m) #They both are equal 
BD3 = norm_(dir_vec(B,D_3))
BF3 = norm_(dir_vec(B,F_3))
n = BD3
print(n) #They both are equal 
CD3 = norm_(dir_vec(C,D_3))
CE3 = norm_(dir_vec(C,E_3))
p = CD3
print(p) #They both are equal 

#Problem 1.5.11
a = norm_(m_BC)
b = norm_(dir_vec(A,C))
c = norm_(m_AB)
#From F3 being a point on AB, AF3 + BF3 = AB similarly 
# m + n = c
# n + p = a 
# p + m = b
mtrx_LHS = [[1,1,0],[0,1,1],[1,0,1]]
mtrx_RHS = [c ,a ,b]
Matrix_mnp = np.linalg.solve(mtrx_LHS,mtrx_RHS)
print(Matrix_mnp)
#m = (c+b-a)/2
#n = (c+a-b)/2
#p = (b+a-c)/2


#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)
x_AD = line_gen(A,D)
x_AD_1 = line_gen(A,D_1)
x_AE_1 = line_gen(A,E_1)
x_BE_1 = line_gen(B,E_1)
x_CF_1 = line_gen(C,F_1)
x_AF_1 = line_gen(A,F_1)
x_CH = line_gen(C,H)
x_BH = line_gen(B,H)
x_AH = line_gen(A,H)
x_OA = line_gen(O,A)
x_OB = line_gen(O,B)
x_OC = line_gen(O,C)
x_ccirc= circ_gen(O,R)
x_IA = line_gen(I,A)
x_IB = line_gen(I,B)
x_IC = line_gen(I,C)
IA_direction = I - A
IB_direction = I - B
IC_direction = I - C
# Normalize the direction vectors
IA_direction /= np.linalg.norm(IA_direction)
IB_direction /= np.linalg.norm(IB_direction)
IC_direction /= np.linalg.norm(IC_direction)
# Extend lines IA, IB, and IC from point I 
x_IA_extended = A + 8 * IA_direction
x_IB_extended = B + 6 * IB_direction #If u want to extend IA,IB,IC we can use these
x_IC_extended = C + 8 * IC_direction
x_icirc= circ_gen(I,r)
x_ID_3 = line_gen(I,D_3)
x_IE_3 = line_gen(I,E_3)
x_IF_3 = line_gen(I,F_3)

#FIgure 1
#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
#Labeling the coordinates
A1 = A.reshape(-1,1)
B1 = B.reshape(-1,1)
C1 = C.reshape(-1,1)
tri_coords = np.block([[A1, B1, C1]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['A', 'B', 'C']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'C' else -10
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
#if using termux
#plt.savefig('tri_sss.pdf')
plt.savefig('/home/aman123/random/code1/figs/figure_1.png')

#Figure 2
plt.figure() #Start a new figure 
#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
#Labeling the coordinates
D1 = D.reshape(-1,1)
E1 = E.reshape(-1,1)
F1 = F.reshape(-1,1)
G = G.reshape(-1,1)
tri_coords = np.block([[A1, B1, C1, D1, E1, F1, G]])
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
plt.savefig('/home/aman123/random/code1/figs/figure_2.png')

#Figure3
plt.figure() #Start a new figure 
#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD_1[0,:],x_AD_1[1,:],label='$AD_1$')
plt.plot(x_BE_1[0,:],x_BE_1[1,:],label='$BE_1$')
plt.plot(x_AE_1[0,:],x_AE_1[1,:],linestyle = 'dashed',label='$AE_1$')
plt.plot(x_CF_1[0,:],x_CF_1[1,:],label='$CF_1$')
plt.plot(x_AF_1[0,:],x_AF_1[1,:],linestyle = 'dashed',label='$AF_1$')
plt.plot(x_CH[0,:],x_CH[1,:],label='$CH$')
plt.plot(x_BH[0,:],x_BH[1,:],label='$BH$')
plt.plot(x_AH[0,:],x_AH[1,:],linestyle = 'dashed',label='$AH$')
#Labeling the coordinates
D_1 = D_1.reshape(-1,1)
E_1 = E_1.reshape(-1,1)
F_1 = F_1.reshape(-1,1)
H = H.reshape(-1,1)
tri_coords = np.block([[A1,B1,C1,D_1,E_1,F_1,H]])
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
plt.savefig('/home/aman123/random/code1/figs/figure_3.png')

#Figure4
plt.figure() #Start a new figure 
#Generating all lines 
plt.plot(x_AB[0, :], x_AB[1, :], label='$AB$')
plt.plot(x_BC[0, :], x_BC[1, :], label='$BC$')
plt.plot(x_CA[0, :], x_CA[1, :], label='$CA$')
plt.plot(x_OA[0,:],x_OA[1,:],label='$OA$')
plt.plot(x_OB[0,:],x_OB[1,:],label='$OB$')
plt.plot(x_OC[0,:],x_OC[1,:],label='$OC$')
#Plotting the circumcircle
plt.plot(x_ccirc[0,:],x_ccirc[1,:],label='$circumcircle$')
# Function to create an arrow with extended length
def create_arrow(ax, start, end, color, label, length_factor=1.0):
    ext_end = start + (end - start) * length_factor
    arrow = plt.arrow(start[0], start[1], ext_end[0] - start[0], ext_end[1] - start[1],
                      head_width=0.1, head_length=0.2, fc=color, ec=color, label=label)
    return arrow

# Arrow for OD, starting at O and extending over D
arrow_OD = create_arrow(plt, O, D, 'blue', '$OD$', length_factor=1.5)

# Arrow for OE, starting at O and extending over E
arrow_OE = create_arrow(plt, O, E, 'red', '$OE$', length_factor=1.5)

# Arrow for OF, starting at O and extending over F
arrow_OF = create_arrow(plt, O, F, 'green', '$OF$', length_factor=1.5)

O = O.reshape(-1,1)

tri_coords = np.block([[A1,B1,C1,O,D1,E1,F1]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','O','D','E','F']
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
plt.savefig('/home/aman123/random/code1/figs/figure_4.png')

#Figure1.4.5
plt.figure()
#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_OA[0,:],x_OA[1,:],label='$OA$')
plt.plot(x_OB[0,:],x_OB[1,:],label='$OB$')
plt.plot(x_OC[0,:],x_OC[1,:],label='$OC$')
#Plotting the circumcircle
plt.plot(x_ccirc[0,:],x_ccirc[1,:],label='$circumcircle$')
#Labeling the coordinates
tri_coords = np.block([[A1,B1,C1,O]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','O']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() 
plt.axis('equal')
plt.savefig('/home/aman123/random/code1/figs/figQ1_4_5.png')

#Figure5
plt.figure()
# Generating all lines
plt.plot(x_AB[0, :], x_AB[1, :], label='$AB$')
plt.plot(x_BC[0, :], x_BC[1, :], label='$BC$')
plt.plot(x_CA[0, :], x_CA[1, :], label='$CA$')
plt.plot(x_IA[0, :], x_IA[1, :], label='$IA$')
plt.plot(x_IB[0, :], x_IB[1, :], label='$IB$')
plt.plot(x_IC[0, :], x_IC[1, :], label='$IC$')
plt.plot(x_ID_3[0, :], x_ID_3[1, :], label='$ID_3$')
plt.plot(x_IE_3[0, :], x_IE_3[1, :], label='$IE_3$')
plt.plot(x_IF_3[0, :], x_IF_3[1, :], label='$IF_3$')


# Plot extended lines IA, IB, and IC
#plt.plot([A[0], x_IA_extended[0]], [A[1], x_IA_extended[1]], label='$IA Extended$')
#plt.plot([B[0], x_IB_extended[0]], [B[1], x_IB_extended[1]], label='$IB Extended$') #Extended 
#plt.plot([C[0], x_IC_extended[0]], [C[1], x_IC_extended[1]], label='$IC Extended$')
#Plotting the circumcircle
plt.plot(x_icirc[0,:],x_icirc[1,:],label='$incircle$')
D_3 = D_3.reshape(-1,1)
E_3 = E_3.reshape(-1,1)
F_3 = F_3.reshape(-1,1)
I = I.reshape(-1,1)
tri_coords = np.block([[A1,B1,C1,D_3,E_3,F_3,I]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D_3','E_3','F_3','I']
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
plt.savefig('/home/aman123/random/code1/figs/figure_5.png')
