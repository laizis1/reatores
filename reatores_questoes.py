#!/usr/bin/env python
# coding: utf-8

# In[17]:


import numpy as np
from scipy.integrate import odeint 
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import scipy.optimize as opt


# ## EXERCICIO 2
#     
# 
# 

# In[25]:


def integral (X, t, k, cS0, tB, vb):
    cB = cS0*(tB - vb*X)
    cS = cS0*(1 - X)
    dXdt = (k*cS*cB)/cS0
    return dXdt


# In[68]:


# S = ESTIRENO; B = BUTADIENO

# DADOS DA QUESTÃO
mms = 104
mmb = 54
v = 27 
ws = 2200
wb = 5000
k = 0.036
t = np.linspace(0,10,11)

# CONCENTRAÇÕES INICIAIS
nS0 = ws/mms
nB0 = wb/mmb
cS0 = nS0/v
cB0 = nB0/v #estireno limitante

# coeficientes estequiometricos

s = 1
b = 3.2
vb = b/s
tB = cB0/cS0
X = odeint(integral,0,t, args = (k, cS0, tB, vb))



np.set_printoptions(precision = 3)
print("t =", t[10], "X = ", X[10]) 

cS = np.zeros(len(t))
cB = np.zeros(len(t))
cP = np.zeros(len(t))
for i in range(len(t)):
    cB[i] = cS0*(tB - vb*X[i])
    cS[i] = cS0*(1 - X[i])
    cP[i] = cS0*(X[i])
    
print("cS = %.2f mol/L" % cS[10])
print("cB = %.2f mol/L" % cB[10])
print("cP = %.2f mol/L" % cP[10])

#PLOT
plt.figure(0, figsize = (8,6))
x, = plt.plot(t,X, label = "X")
plt.title("Conversão", fontsize = 20)
plt.xlabel("t (h)", fontsize = 14)
plt.ylabel("X", fontsize = 14)
plt.legend(handles=[x], fontsize = 12)

plt.figure(1, figsize = (8,6))
a, = plt.plot(t, cS, label = "cS")
b, = plt.plot(t, cB, label = "cB")
c, = plt.plot(t, cP, label = "cP")
plt.legend(handles=[a,b,c], fontsize = 12)
plt.title("Concentração", fontsize = 20)
plt.xlabel("t (h)", fontsize = 14)
plt.ylabel("Concentração (mol/L)", fontsize = 14)


# ## EXERCÍCIO 3

# In[32]:


# A = C4H8O; B = H2O; C = C4H10O2

def integral2 (X, t, cA0, cB0, k):
    cA = cA0*(1-X)
    dXdt = (k*cA*cB0)/cA0
    return dXdt

#DADOS DA QUESTÃO
cA0 = 0.25
k = 8.3e-4
t = np.linspace(1e-5,45,16)
mmb = 18
v = 1000
cB0 = v/mmb #como agua está em em excesso cB0 = cB
X = odeint(integral2,1e-5,t, args = (cA0, cB0, k))

cA = np.zeros(len(t))

cC = np.zeros(len(t))
for i in range(len(t)):
    cA[i] = cA0*(1 - X[i])
    cC[i] = cA0*(X[i])

print("cA = %.3f mol/L" % cA[15])
print("cB = %.1f mol/L" % cB0)
print("cC = %.2f mol/L" % cC[15])


#  ##  EXERCÍCIO 4

# In[33]:



#DADOS DA QUESTAO

nbase = 2.5 # mol/min
E = 85000
R = 8.314
T0 = 323
k0 = 1e-4
X = 0.9
T = 400
P = 101e4 # Pa

#CALCULO DO k e t reacao
k = k0*np.exp((E/R)*(1/T0 - 1/T))
tr = np.log(1/(1-X))*(1/k)
#calculo do Volume
nA = nbase*tr
V = ((nA*R*T)/P)*1000
nciclos = 3600/nA
nhoras = (tr*nciclos)/60
print("tempo de reação = %.0f min" % tr)
print("ciclos = %.0f, horas = %.2f" % (nciclos, nhoras))
print("V = %.2f L" %V)


# ## EXERCÍCIO 5

# In[34]:


def edo1(cA, t, k, cA0, fA0, V):
    Vt = V + fA0*t
    dcAdt = (-k*cA*Vt + cA0*fA0 + fA0*cA)/Vt
    return dcAdt

def edo2(cA, t, k, V):
    dcAdt =  -k*cA
    return dcAdt


# In[69]:


Vs = 50
V = 120
cA0 = 1
fA0 = 10
te = (V-Vs)/fA0
t = np.linspace(0,te,16)
k = 1.5
cA2 = np.zeros(len(t))
# LETRA A - CONSIDERAMOS COMO UM SEMI BATELADA
cA1 = odeint(edo1,0,t, args = (k, cA0, fA0, V))


#LETRA B - CONSIDERAMOS COMO BATELADA
cA2 = odeint(edo2,Vs*cA0/V,t, args = (k, V))

#PLOT
plt.figure(0, figsize = (8,6))
a, = plt.plot(t, cA1, label = "cA1")
plt.legend(handles=[a])
plt.title("Concentração", fontsize =20)
plt.xlabel("t (min)", fontsize =14)
plt.ylabel("Concentração (mol/L)", fontsize = 14)

plt.figure(1, figsize = (8,6))
a, = plt.plot(t, cA2, label = "cA2", color = "g")
plt.legend(handles=[a])
plt.title("Concentração", fontsize = 20)
plt.xlabel("t (min)", fontsize = 14)
plt.ylabel("Concentração (mol/L)", fontsize = 14)


# ## EXERCICIO 6

# In[38]:



#DADOS DA QUESTÃO

def f(X, K, vb, cA0):
    cB = cA0*(vb*X)
    cA = cA0*(1-X)
    y = K - (cB)**2/cA
    return y

T = 340
R = 0.082
P = 2
K = 0.1
vb = 2
cA0 = P/(R*T)
k = np.linspace(1, 5, 5)

X = opt.newton(f, 0.5, args =(K, vb, cA0))
print("X = %.3f" % X)
x = np.linspace(0, X-0.01, 11)
t = [[0]*len(x)]*len(k)


##for j in range (len(k)):
    ##for i in range(len(x)):
    ## r = k*(cA0*(1-x[i]) - (cA0*vb*x[i])**2/K)
    ##t[i] = (cA0/r)*x[i]
    ##t[i] = np.log(1/(1-x[i]))/(k*(1 - 0.1*K))

t = [[np.log(1/(1-x[i]))/(k[j]*(1 - 0.1*K)) for i in range(len(x))] for j in range(len(k))] 
print("Valores de k = ", k)
ff = "{0:.2f}"
print("Valores de t em X +/- = Xeq, para os valores de k usados:")
for i in range(len(k)):
    print("k =", k[i], " ", "t =", ff.format(t[i][10]))

fig = plt.figure(0, figsize = (8,6))
for i in range(len(t)):
    plt.plot(t[i],x)

plt.xlabel("t(?)", fontsize = 18)
plt.ylabel("Conversão", fontsize = 18)
plt.legend(["k = 1", "k = 2", "k = 3", "k = 4", "k = 5"])

    


# ## EXERCÍCIO 7

# In[66]:


#DADOS DA QUESTÃO

def g(X, V, k, cA0, v0):
    dXdV = k*(cA0*(1-X))**2/(cA0*v0)
    return dXdV

k = 1
V = np.linspace(0,2,11)
v0 = 0.01
cA0 = 0.02

# LETRA A

X1 = odeint(g,0,V, args = (k, cA0, v0))

print("Conversão Tubular = ", X1[10])

# LETRA B
tau = V[10]/v0 
Da = tau*k*cA0 #segunda ordem
X2 = ((1 + 2*Da)-np.sqrt(1+4*Da))/(2*Da)

print("Conversão CSTR = %.2f" % X2)

