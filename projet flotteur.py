import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math

fichier = open("liste_rho.txt")
lignes = fichier.readlines()
lignes.pop(0)
fichier.close()

L = []
for i in lignes:
    L.append(i.split(" "))

profondeur = []
rho = []

for i in L:
    i[1] = i[1][:-1]
    profondeur.append(float(i[0]))
    rho.append(float(i[1]))


f = interpolate.interp1d(profondeur, rho, kind = 'quadratic')
x = np.linspace(min(profondeur), max(profondeur), 50)
y = f(x)

def tracer_rho():
    plt.plot(np.array(profondeur), np.array(rho), '+', label = 'expérimentale')
    plt.plot(x, y, label = 'interpolation')
    plt.title('rho = f(z)')
    plt.xlabel("profondeur en m")
    plt.ylabel("rho en kg/m3")
    plt.grid()
    plt.legend()
    plt.show()

def f1(depth):
    if depth < min(profondeur) or depth > max(profondeur):
        if depth < min(profondeur):
            depth = min(profondeur)
        if depth > max(profondeur) :
            depth = max(profondeur)
    return f(depth)


mf = 60 #en kg
V = mf/f(10) #en m3
D = 11.3e-2 #en m
S = np.pi*(D/2)**2 #en m3
C = 1
N = 10000
n = 0.00108
t, dt = np.linspace(0, 10000, N+1, retstep=True)
h=dt
g = 9.8 #m/s-2



def g(X, V):
    depth = X[0]
    rho_interp = f1(depth)
    ma = rho_interp*4*np.pi*(D/2)**3 /3
    m = ma+mf
    return np.array([X[1], (-rho_interp*V/m + mf/m)*9.8- (rho_interp*S*C*X[1]*abs(X[1])+3*np.pi*n*D*X[1])/(2*m)])


y = np.zeros(N+1)
X = np.zeros((N+1, 2))
v = np.zeros(N+1)

X[0] = [x[0], 0]
y[0] = X[0][0]
v[0] = X[0][1]
def euler():
    for i in range(N):
        X[i+1] = X[i] + h*g(X[i])
        y[i+1] = X[i+1][0]
        v[i+1] = X[i+1][1]
    plt.plot(t, np.array(y), '--', label='Profondeur')
    plt.plot(t, v, label='vitesse')
    plt.xlabel("temps")
    plt.title("z = f(t)")
    plt.grid()
    plt.legend()
    plt.show()



Patm =1.01325e5 #Pa
K = 2.2e9 #Pa

def Pression(z):
    return  z/10

def f2(z):
    P = Pression(z)
    return 1025.22*(1+(P-Patm)/K)
dVmax = (np.pi*(42.5e-3)**2) *0.3/120
print(dVmax)

nvm = 0.9
pas = 5e-3  #en m
diam= 25e-3 #en m
tani = pas/(np.pi*diam)
phi = math.atan(tani/nvm)-math.atan(tani)
nem = math.tan(math.atan(tani)-phi)/tani
print('nem', nem,"1/nem", 1/nem, 'nvm', nvm)




def correc(kp, ki, kd, z_cible):
    V = mf/1025.22
    V1 = np.zeros(N+1)
    V1[0] = V
    z_point = np.zeros(N+1)
    z_2point = np.zeros(N+1)
    dVmax = 2.2578e-6*dt
    dV1 = np.zeros(N+1)
    E = np.zeros(N+1)
    for i in range(N):
        k1 = h*g(X[i], V)
        k2 = h*g(X[i] + k1/2, V)
        k3 = h*g(X[i]+ k2/2, V)
        k4 = h*g(X[i]+k3, V)
        X[i+1] = X[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        y[i+1] = X[i+1][0]
        v[i+1] = X[i+1][1]
        z_point[i+1] = (y[i+1]-y[i])/dt
        z_2point[i+1] = (z_point[i+1]-z_point[i])/dt
        dV = (kp * z_point[i+1] - ki*(z_cible - y[i + 1]) + kd * z_2point[i+1])
        dV = min(max(dV, -dVmax), dVmax)
        dV1[i+1] = dV
        V += dV*dt
        V1[i+1] = V
        if dV>0:
            E[i+1] = E[i] + dV*dt*Pression(y[i+1])*10**5/nem
        if dV<0:
            E[i+1] = E[i] + nvm*dV*dt*Pression(y[i+1])*10**5
    plt.plot(t, np.array(y), '--', label='Profondeur')
    plt.plot(t, v, label='vitesse')
    plt.xlabel("temps")
    plt.title("z = f(t)")
    plt.grid()
    plt.legend()
    plt.show()
    plt.plot(t, dV1, label='Debit')
    plt.xlabel("temps")
    plt.title("dV = f(t)")
    plt.grid()
    plt.legend()
    plt.show()
    plt.plot(t, V1, label='Volume')
    plt.xlabel("temps")
    plt.title("dV = f(t)")
    plt.grid()
    plt.legend()
    plt.show()
    plt.plot(t, E, label='Energie')
    plt.xlabel("temps")
    plt.title("E = f(t)")
    plt.grid()
    plt.legend()
    plt.show()

def correc_traj(kp, ki, kd, z_cible_tab, N_tab):
    V = mf/1025.22
    for i in range(N):
        z_cible = trajectoire(i,z_cible_tab, N_tab)
        k1 = h*g(X[i], V)
        k2 = h*g(X[i] + k1/2, V)
        k3 = h*g(X[i]+ k2/2, V)
        k4 = h*g(X[i]+k3, V)
        X[i+1] = X[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        y[i+1] = X[i+1][0]
        v[i+1] = X[i+1][1]
        V1 = mf / f2(Pression(y[i + 1]))
        V_prec = mf / f1(y[i])
        dV = kp * (mf/f2(Pression(z_cible))-V1) + ki*(mf/f2(Pression(z_cible))-V1)*dt + kd * ((mf/f2(Pression(z_cible))-V1) - (mf/f2(Pression(z_cible))-V_prec))/dt
        V += dV
    plt.plot(t, np.array(y), '--', label='Profondeur')
    plt.plot(t, v, label='vitesse')
    plt.xlabel("temps")
    plt.title("z = f(t)")
    plt.grid()
    plt.legend()
    plt.show()

def trajectoire(i, z_cible_tab, N_tab):
    if i <= N_tab[0]:
        z_cible = z_cible_tab[0]
    elif i > N_tab[1]:
        z_cible = z_cible_tab[2]
    else:
        z_cible = z_cible_tab[1]
    return z_cible



correc(1e-4, 1e-7, 1e-2, 20)