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
N = 20000
h = 0.05
n = 0.00108
t = np.linspace(0, 100, N+1)
g = 9.8 #m/s-2
dt = 100/N


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

Patm =1.01325e5 #Pa
K = 2.2e9 #Pa

def Pression(z):
    return 10* z + Patm

def f2(z):
    P = Pression(z)
    return f1(min(profondeur))*(1+(P-Patm)/K)

def correc(kp, ki, kd, z_cible):
    V = mf/1025.22
    for i in range(N):
        k1 = h*g(X[i], V)
        k2 = h*g(X[i] + k1/2, V)
        k3 = h*g(X[i]+ k2/2, V)
        k4 = h*g(X[i]+k3, V)
        X[i+1] = X[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        y[i+1] = X[i+1][0]
        v[i+1] = X[i+1][1]
        V = mf / f1(y[i + 1])
        V_prec = mf / f2(y[i])
        dV = kp * (mf/f2(z_cible)-V) + ki*(mf/f2(z_cible)-V)*dt + kd * ((mf/f2(z_cible)-V) - (mf/f2(z_cible)-V_prec))/dt
        V += dV
    plt.plot(t, np.array(y), '--', label='Profondeur')
    plt.plot(t, v, label='vitesse')
    plt.xlabel("temps")
    plt.title("z = f(t)")
    plt.grid()
    plt.legend()
    plt.show()

correc(1, 1, 1, 10)