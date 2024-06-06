import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

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


mf = 40 #en kg
V = mf/f(10) #en m3
D = 11.3e-2 #en m
S = np.pi*(D/2)**2 #en m3
C = 1
N = 30000
h = 0.05
n = 0.00108
t = np.linspace(0, N, N+1)
g = 9.8 #m/s-2
dt = 100/N
Vmin = mf/1032.09
Vmax = mf/1025.22


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

#euler()

def correc(kp, ki, kd, z_cible_tab, N_tab):
    V = mf/f1(0)
    for i in range(N):
        z_cible = trajectoire(i,z_cible_tab, N_tab)
        k1 = h*g(X[i], V)
        k2 = h*g(X[i] + k1/2, V)
        k3 = h*g(X[i]+ k2/2, V)
        k4 = h*g(X[i]+k3, V)
        X[i+1] = X[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        y[i+1] = X[i+1][0]
        v[i+1] = X[i+1][1]
        V = mf / f1(y[i + 1])
        V_prec = mf / f1(y[i])
        dV = kp * (mf/f1(z_cible)-V) + ki*(mf/f1(z_cible)-V)*dt + kd * ((mf/f1(z_cible)-V) - (mf/f1(z_cible)-V_prec))/dt
        V += dV
    plt.plot(t, -np.array(y), '--', label='Profondeur')
    plt.plot(t, v, label='vitesse')
    plt.xlabel("temps")
    plt.title("z = f(t)")
    plt.grid()
    plt.legend()
    plt.show()

def P2(z_cible):
    P = (z_cible / 10 + 1)
    return P


def correc2(kp, ki, kd, z_cible_tab, N_tab, Vmin, Vmax):
    V = mf/f1(0)
    V_tab = [V]
    for i in range(N):
        z_cible = trajectoire(i,z_cible_tab, N_tab)
        k1 = h*g(X[i], V)
        k2 = h*g(X[i] + k1/2, V)
        k3 = h*g(X[i] + k2/2, V)
        k4 = h*g(X[i]+k3, V)
        X[i+1] = X[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        y[i+1] = X[i+1][0]
        v[i+1] = X[i+1][1]
        P = P2(y[i+1])
        P_prec = P2(y[i])
        dV = kp * (P-P2(z_cible)) + ki*(P-P2(z_cible))*dt+ kd * ((P-P2(z_cible)) - (P_prec-P2(z_cible)))/dt
        dV = regulation_dV(dV)
        V += dV
        V = regulation_V(V, Vmin, Vmax)
        V_tab.append(V)
    plt.plot(t, -np.array(y), '--', label='Profondeur')
    plt.plot(t, v, label='vitesse')
    plt.xlabel("temps")
    plt.title("z = f(t)"+', kp = '+str(kp)+', ki = '+str(ki)+', kd = '+str(kd))
    plt.grid()
    plt.legend()
    plt.show()

    # plt.plot(t, V_tab, label='Volume')
    # plt.xlabel("temps")
    # plt.title("Volume du flotteur")
    # plt.grid()
    # plt.legend()
    # plt.show()

def trajectoire(i, z_cible_tab, N_tab):
    if i <= N_tab[0]:
        z_cible = z_cible_tab[0]
    elif i > N_tab[1]:
        z_cible = z_cible_tab[2]
    else:
        z_cible = z_cible_tab[1]
    # z_cible  = z_cible_tab[0]
    return z_cible

def regulation_V(V, Vmin, Vmax):
    if V < Vmin:
        V = Vmin
    if V > Vmax:
        V = Vmax
    return V

def regulation_dV(dV):
    facteur = 2/7
    limite = 1.25*42.5**2*3.14*0.000000001*facteur
    if abs(dV) > limite:
        if dV < 0 :
            dV = -limite
        else :
            dV = limite
    return dV


correc2(50, 10, 120, [10,25,2], [5*N/11,7*N/11], Vmin, Vmax)