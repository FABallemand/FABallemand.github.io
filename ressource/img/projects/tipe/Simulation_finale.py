## Simulation Principale:

## Initialisation:

import matplotlib.pyplot as plt
from math import *
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import copy as c

R=2.40 #Grand rayon du tore
rho=0.7 #Rayon du tore
I=1400 #Intensité du courant dans les bobines
N=18*2028 #Nombre de spires
q=1.6*(10)**-19 #Charge élémentaire
µ=4*pi*10**-7 #Perméabilité du vide
m=9.1*10**-31 #Masse e-
B=(q*µ*N*I)/(2*pi*m)

#v0=1.3*10**4 #Conservation de la quantité de mouvement
#v0=10**11
#v0=1.2*10**5
#v0=8*10**7 #D'après Thermodynamique (T=150*10**6K)
v0=10**6

t0=0
r0=R
theta0=0
z0=0

p_dicho=10**-10

T_calcul=5.4*10**-11
nbr_points=10000
t=np.linspace(0,T_calcul,nbr_points)
#maxstep=50000000
maxstep=10
#maxstep=1000000 # Maximum number of (internally defined) steps allowed for each integration point in t.

Racine = "F:/TIPE/Info/"
FileName = str(v0)+"_"+str(T_calcul)+"_"+str(nbr_points)+"_"+str(maxstep)

def FileWriteCsv(file,array):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            file.write(str(array[i,j]))
            if (j != array.shape[1]-1):
                file.write(",")
        file.write("\n")

def linsplt(l):
    return list(map(float,l.split(",")))

def FileReadCsv(file):
    lines = file.readlines()
    return np.array(list(map(linsplt,lines)))

def f(r):
    #fonction f issue du syst.
    return 2*(v0**2)-((R*v0)**2/(r**2))-(B*log(r/R))**2

def dichotomie(f,r1,r2,p):
    #Prend en argument une fonction f, un intervalle [r1,r2] et une précision et détermine le point où f s'annule avec la précision p
    while r2-r1>p:
        m=(r1+r2)*0.5
        if f(m)*f(r1)<0:
            r2=m
        else:
            r1=m
    return((r1+r2)*0.5)

def valeurs_negatives(f,R):
    #Prend en argument une fonction f et R et détermine deux valeurs à gauche et à droite de R pour lesquelles la fonction f est négative
    x1,x2=R,R
    while f(x1)>0:
        x1=x1/2
    while f(x2)>0:
        x2=2*x2
    x1=x1/2
    x2=2*x2
    return(x1,x2)

def zone_tore(f,R,p):
    #prend en argument une fonction f, le grand rayon du tore R et une "précision" p et détermine l'intervalle de r que peut atteindre la particule
    x1i,x2i=valeurs_negatives(f,R)
    x1=dichotomie(f,x1i,R,p)
    x2=dichotomie(f,R,x2i,p)
    return(x1,x2)

def instant_collision(dans_tore):
    i=0
    intervalle_temp=T_calcul/nbr_points
    for i in range(len(dans_tore)):
        if dans_tore[i]<=0:
            return(i,i*intervalle_temp)
    return("pas de sortie de tore","pas de sortie de tore")

## Enregistrement:

def f_s(r):
    return (((R*v0)**2)/r**3)-(B**2)*(log(r/R)/r)

def f_theta(r):
    return (R*v0)/(r**2)

def f_z(r):
    return B*log(r/R)

def dY(Y,t):
    r,s,theta,z=Y
    dYdt=[s,f_s(r),f_theta(r),f_z(r)]
    return dYdt

Y0=[r0,v0,theta0,z0] #CI

from scipy.integrate import odeint
sol = odeint(dY,Y0,t,mxstep=maxstep)
#sol = odeint(dY,Y0,t)

dans_tore=[rho**2 -(sol[i,0]-R)**2-sol[i,3]**2 for i in range(len(sol))]
Indice_collision,T_collision=instant_collision(dans_tore)
print(Indice_collision,T_collision)

file = open(Racine+FileName+".csv","w")
FileWriteCsv(file,sol)
file.close()

## Lecture:

#file = open(Racine+FileName+".csv","r")
file = open(Racine+"10000000000_9.000000000000001e-09_10000_1000000.csv","r")
sol = FileReadCsv(file)
file.close()

dans_tore=[rho**2 -(sol[i,0]-R)**2-sol[i,3]**2 for i in range(len(sol))]
Indice_collision,T_collision=instant_collision(dans_tore)
print(Indice_collision,T_collision)

## Visualisation:

## Générale:
solg=plt.figure("Méthode scipy",figsize = (16, 9))
plt.gcf().subplots_adjust(left = 0.13, bottom = 0.1,right = 0.98, top = 0.95, wspace = 0.2, hspace = 0.4)

ax = solg.add_subplot(3, 3, 1) #Table
data=[[R,"m"],
      [rho,"m"],
      [I,"A"],
      [N,""],
      [q,"C"],
      [m,"kg"],
      [B,""],
      [v0,"m/s"],
      [(r0,theta0,z0),"(m, rad, m)"],
      [T_calcul,"s"],
      [T_collision,"s"]
      ]
column_labels=["Valeur Numérique","Unité"]
row_labels=["Grand rayon R","Petit rayon p","Intensité courant bobines I","Nombre de Bobines","Charge particule q", "Masse particule m",r"Constante $\beta$","Vitesse initiale particule","Coordonnée initiales particules","Intervalle de temps de calcul","Sortie de tore"]
plt.axis('tight')
plt.axis('off')
plt.table(cellText=data,colLabels=column_labels,rowLabels=row_labels,loc="center")

ax = solg.add_subplot(3, 3, 2)
plt.plot(t,sol[:,0],label="r en fonction de t")
plt.legend()

ax = solg.add_subplot(3, 3, 3)
plt.plot(t,sol[:,1],label="s en fonction de t")
plt.legend()

ax = solg.add_subplot(3, 3, 4) #Portrait de phase
plt.plot(sol[:,0],sol[:,1],label="Portrait de phase",color="lightblue")
plt.legend()

ax = solg.add_subplot(3, 3, 5)
plt.plot(t,sol[:,2],label=r"$\theta$ en fonction de t")
plt.legend()

ax=solg.add_subplot(3,3,6)
plt.plot(t,sol[:,3],label="z en fonction de t")
plt.legend()

ax=solg.add_subplot(3,3,8)
plt.plot(t,dans_tore,label="positif<=>particule dans le tore")
plt.yticks([0])
plt.grid()

plt.show()

## r(t):
plt.figure()
plt.plot(t, sol[:, 0], label='r(t)')
plt.legend(loc='best')
plt.xlabel('t',fontsize=30)
plt.ylabel('r(t)',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.xlim(0,10**-6)
plt.tight_layout()
plt.show()

## r(t) (avec horizontales):
x1p,x2p=zone_tore(f,R,p_dicho)
print(x1p,x2p)

plt.figure()
Y1=[x1p for i in range(len(sol[:, 0]))]
Y2=[x2p for i in range(len(sol[:, 0]))]
plt.plot(t, sol[:, 0], label='r(t)')
plt.plot(t,Y1,"red")
plt.plot(t,Y2,"red")
plt.legend(loc='best')
plt.xlabel('t',fontsize=30)
plt.ylabel('r(t)',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
plt.show()

## z(t):
plt.figure()
plt.plot(t, sol[:, 3], label='z(t)')
plt.legend(loc='best')
plt.xlabel('t',fontsize=30)
plt.ylabel('z(t)',fontsize=30)
plt.xticks([0*10**-4,0.5*10**-4,10**-4],fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
plt.show()

## z(t) (réduit):
plt.figure()
plt.plot(t, sol[:, 3], label='z(t)')
plt.legend(loc='best')
plt.xlabel('t',fontsize=30)
plt.ylabel('z(t)',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0,10**-7)
plt.ylim(0,2.5*10**-5)
plt.tight_layout()
plt.show()

## Négatif <=> Sortie du tore:

plt.figure()
tore=[rho**2 -(sol[i,0]-R)**2-sol[i,3]**2 for i in range(len(sol))]

plt.plot(t,tore)
plt.xticks([0,5*10**-5,10**-4],fontsize=20)
plt.yticks(fontsize=20)
plt.grid()
plt.tight_layout()
plt.show()

##Négatif <=> Sortie du tore (réduit):
plt.figure()
plt.plot(t,dans_tore,label="positif<=>particule dans le tore")
plt.plot(t,[0 for i in range(len(t))],color="red")
plt.plot([T_collision for i in range(10)],[i-1 for i in range(10)],color="red")
plt.yticks([0],fontsize=20)
plt.xticks([0,5.4*10**-11,6*10**-11],fontsize=20)
#plt.xticks(fontsize=20)
plt.xlim(0,6*10**-11)
plt.ylim(-0.1,0.5)
plt.grid()
plt.tight_layout()
plt.show()

## Tore:

def pol_cart(list):
    r,theta,z = list[0],list[2],list[3]
    x=r*cos(theta)
    y=r*sin(theta)
    return(x,y,z)

def l_pol_cart(sol):
    X,Y,Z=[],[],[]
    for i in range(len(sol)):
        x,y,z=pol_cart(sol[i])
        X.append(x)
        Y.append(y)
        Z.append(z)
    return(X,Y,Z)

X,Y,Z=l_pol_cart(sol)
tore = plt.figure("Tore",clear=True)
ax = tore.add_subplot(1, 1, 1, projection='3d')

(theta, phi) = np.meshgrid(np.linspace(0, 2 * np.pi, 41),
                           np.linspace(0, 2 * np.pi, 41))

x = (1.9 + np.cos(phi)) * np.cos(theta)
y = (1.9 + np.cos(phi)) * np.sin(theta)
z = 0.7*np.sin(phi)

dplot = ax.plot_surface(x,y,z,color='whitesmoke' ,alpha=0.5)
ax.set(xlabel='x',
       ylabel='y',
       zlabel='z',
       xlim = [-3.2, 3.2],
       ylim = [-3.2, 3.2],
       zlim = [-1, 1],
       xticks = [-3.1,-2.4,-1.7,0,1.7,2.4,3.1],
       yticks = [-3.1,-2.4,-1.7,0,1.7,2.4,3.1],
       zticks = [-0.7, 0, 0.7],
       title='Tore')
ax.scatter(X, Y, Z, marker='p')
tore.tight_layout()
#ax.view_init(10, 45)
#ax.view_init(90,90)
ax.view_init(0,90)
plt.show()
##
for angle in range(0, 360):
    ax.view_init(angle, angle)
    plt.draw()
    plt.pause(.001)