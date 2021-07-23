#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import scipy.integrate as ing
from scipy.integrate import odeint 
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import scipy.optimize as opt


# In[2]:


#DADOS FORNECIDOS

R1 = 0.082
R2 = 8.314
FA0 = 2.5 #mol/min
k0 = 5e-3 #min^-1 
T0 = 150 + 273.15 #K
T = 227 + 273.15 #K
E = 85000 #J/mol
P = 10 #atm
V = 830 #L
vA = 1
vB = 1
vC = 2

#DADOS ADICIONAIS CALCULADOS

yA0 = 1
delta = vB + vC - vA
e = yA0*delta
cA0 = P/(R1*T)
k = k0*np.exp((E/R2)*(1/T0 - 1/T))


# ## QUESTÃO 1 ##
# 

# In[3]:


#CONCENTRAÇÃO DE A E TAXA DE REAÇÃO
X = np.linspace(0, 0.9, 10)
cA = np.zeros(len(X))
Vcstr = np.zeros(len(X))
for i in range(len(X)):
    cA[i] = cA0*(1-X[i])/(1 + e*X[i])
    rA = k*cA[i]
    Vcstr[i] = FA0*X[i]/rA

print("Vcstr = %.0f L \nV = %.0f L" % (Vcstr[9] , V))
print("O reator existente não é adequado")

#PLOT
plt.figure(0, figsize=(10, 8), dpi=80)
plt.plot(Vcstr, X,)
plt.xlabel("Volume (L)")
plt.ylabel("X")
plt.savefig("q1a")


# ## QUESTÃO 2 ##

# In[4]:


#Foram considerados dois arranjos, dois cstrs em série e um cstr seguido por um pfr

#CONSIDERANDO DOIS CSTRS EM SÉRIE

def f(X, V, cA0, e, k, FA0):
    cA = cA0*(1-X)/(1 + e*X)
    y = V - FA0*X/(k*cA)
    return y

def g(V, X, FA0, cA0, k, e):
    dVdX = (FA0*(1+e*X))/(k*cA0*(1-X))
    return dVdX 
    
V1 = np.linspace(0, V, 10)
X1 = np.zeros(len(V1))
cA1 = np.zeros(len(V1))

for i in range(len(V1)):
    X1[i] = opt.bisect(f, 0, 0.99, args = (V1[i], cA0, e, k, FA0)) #método da bisseção usado para zerar a função f, encontrando para cada uma dos valores de V qual é a conversão
    cA1[i] = cA0*(1-X1[i])/(1 + e*X1[i]) #cálculo da concentração

X2 = np.linspace(X1[9], 0.9, 10) #pro segundo reator a conversão vai de 0,86 até 0,9
V2 = np.zeros(len(X2))
cA2 = np.zeros(len(X2))

for i in range(len(X2)):
    cA2[i] = cA0*(1-X2[i])/(1 + e*X2[i]) #cálculo da concentração
    V2[i] = FA0*(X2[i] - X2[0])/(k*cA2[i]) #cálculo do V


#CONSIDERANDO CSTR + PFR
V21 = odeint(g,0, X2, args = (cA0, e, k, FA0)) #cálculo usando a biblioteca do odeint para equações diferenciais presente no python
cA21 = np.zeros(len(X2))

for i in range(len(X2)):
    cA21[i] = cA0*(1-X2[i])/(1 + e*X2[i])

print("(V1, X1) = (%.2f, %.2f)" % (V1[9], X1[9]))
print("(V2, X2) CSTR = (%.2f, %.2f)" % (V2[9], X2[9])) 
print("(V2, X2) Tubular = (%.2f, %.2f )" % (V21[9], X2[9]))

#PLOT
fig, ax1 = plt.subplots(figsize=(10,8))
ax1.set_title("CSTR 1")
ax1.set_xlabel('Volume(L)')
ax1.set_ylabel('X', color='tab:red')
ax1.plot(V1, X1, color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax2 = ax1.twinx() 
ax2.set_ylabel('Concentração', color="tab:blue")  
ax2.plot(V1, cA1, color="tab:blue")
ax2.tick_params(axis='y', labelcolor="tab:blue")
fig.tight_layout()
plt.savefig("q2a1")

fig1, ax3 = plt.subplots(figsize=(10,8))
ax3.set_title("CSTR 2")
ax3.set_xlabel('Volume(L)')
ax3.set_ylabel('X', color='tab:red')
ax3.plot(V2, X2, color='tab:red')
ax3.tick_params(axis='y', labelcolor='tab:red')
ax4 = ax3.twinx() 
ax4.set_ylabel('Concentração', color="tab:blue")  
ax4.plot(V2, cA2, color="tab:blue")
ax4.tick_params(axis='y', labelcolor="tab:blue")
plt.savefig("q2a2")

fig1, ax3 = plt.subplots(figsize=(10,8))
ax3.set_title("PFR 2")
ax3.set_xlabel('Volume(L)')
ax3.set_ylabel('X', color='tab:red')
ax3.plot(V21, X2, color='tab:red')
ax3.tick_params(axis='y', labelcolor='tab:red')
ax4 = ax3.twinx() 
ax4.set_ylabel('Concentração', color="tab:blue")  
ax4.plot(V21, cA21, color="tab:blue")
ax4.tick_params(axis='y', labelcolor="tab:blue")
plt.savefig("q2a3")


# ## QUESTÃO 3

# In[5]:


# SOLUÇÃO ANALÍTICA

V3 = odeint(g, 0, X, args = (FA0, cA0, k, e))
cA3 = np.zeros(len(X))
for i in range(len(X)):
    cA3[i] = cA0*(1-X[i])/(1 +e*X[i])

# SOLUÇÃO NUMÉRICA
foo = np.zeros(len(X))
for i in range(len(X)):
    rA = k*cA0*(1-X[i])/(1+e*X[i]) 
    foo[i] = FA0/rA #cálculo para o gráfico de levenspiel
    
fig, ax1 = plt.subplots(figsize=(10,8))
ax1.set_title("Levenspiel")
ax1.set_xlabel('X')
ax1.set_ylabel('FA0/rA')
ax1.scatter(X, foo, color='tab:red')
fig.tight_layout()  
plt.savefig("q3a1")


# In[6]:


#SOLUÇÃO NUMÉRICA

#Calculamos usando a curva de FA0/rA, dividindo o cálculo da área da curva em duas partes: primeiro otimizamos para achar uma função para o cálculo da área até 0,5 de conversão e depois outra função de 0,5 até 0,9
X31 = X[0:6] # função 1
X32 = X[5:10] # função 2
lev1 = foo[0:6] # função 1
lev2 = foo[5:10] # função 2
cA31 = np.zeros(len(X31)) 
cA32 = np.zeros(len(X32))
deg1 = 2 #otimizamos para achar um polinomio de grau 2 no intervalo de conversao 0-0,5
deg2 = 3 #otimizamos para achar um polinomio de grau 3 no intervalo de conversao 0,5-0,9
c1 = np.polynomial.polynomial.polyfit(X31, lev1, deg1) # polyfit acha os valores das constantes do polinômio
c2 = np.polynomial.polynomial.polyfit(X32, lev2, deg2) # polyfit acha os valores das constantes do polinômio

def g(X,a,b,c,d,e,f):
    return a + b*X + c*X**2 + d*X**3 + e*X**4 + f*X**5 #funcao polinomial para achar os valores das constantes

#CÁLCULO DO R-SQUARED
yfit1 = np.zeros(len(X31))
yfit2 = np.zeros(len(X32))

for i in range(len(X31)):
    yfit1[i] = g(X31[i] ,c1[0], c1[1], c1[2] , 0,0, 0)
for i in range(len(X32)):
    yfit2[i] = g(X32[i] , c2[0], c2[1], c2[2], c2[3] , 0,0)

yr1 = lev1 - yfit1
yr2 = lev2 - yfit2
ssr1 = sum(pow(yr1,2))
ssr2 = sum(pow(yr2,2))
sst1 = len(lev1)*np.var(lev1)
sst2 =  len(lev2)*np.var(lev2)
rsq1 = 1 - ssr1/sst1
rsq2 = 1 - ssr2/sst2
print("R1 = %.4f, R2 = %.4f" % (rsq1,rsq2))


c11 = np.zeros(len(c1))
c21 = np.zeros(len(c2))

for i in range(len(c1)):
    c11[i] = c1[i]/(i+1) # cálculo das constantes depois da integração do polinômio

for i in range(len(c2)):
    c21[i] = c2[i]/(i+1) # cálculo das constantes depois da integração do polinômio

V4 = np.zeros(len(X31)) 
V5 = np.zeros(len(X32))
for i in range(len(V4)):
    V4[i] = c11[0]*X31[i] + c11[1]*X31[i]**2 + c11[2]*X31[i]**3 
    cA31[i] = cA0*(1 - X31[i])/(1 +e*X31[i])
for i in range(len(V5)):
    V5[i] = c21[0]*(X32[i] - X31[5]) + c21[1]*(X32[i]**2 - X31[5]**2) + c21[2]*(X32[i]**3 - X31[5]**3) + c21[3]*(X32[i]**4 - X31[5]**4)
    V5[i] = V5[i] + V4[5]
    cA32[i] = cA0*(1 - X32[i])/(1 +e*X32[i])

V6 = np.concatenate((V4, V5), axis = None)
V6 = np.delete(V6, 5)
cA6 = np.concatenate((cA31, cA32), axis = None)
cA6 = np.delete(cA6, 5)


# In[ ]:


# GRAFICOS DE CONVERSÃO POR VOLUME DE TUBULAR X CSTR
fig, ax1 = plt.subplots(figsize=(10,8))
ax1.set_title("Conversão")
ax1.set_xlabel('Volume (L)')
ax1.set_ylabel('X')
ax1.plot(V6, X, color='tab:red')
ax1.plot(Vcstr, X, color='tab:blue')
ax1.legend(["tub", "cstr"], loc = 0)
fig.tight_layout()
plt.savefig("q3a2")

# GRAFICOS DE CONCENTRAÇÃO POR VOLUME DE TUBULAR X CSTR
fig, ax1 = plt.subplots(figsize=(10,8))
ax1.set_title("Concentração")
ax1.set_xlabel('Volume (L)')
ax1.set_ylabel('mol/L')
ax1.plot(V6, cA6, color='tab:red')
ax1.plot(Vcstr, cA, color='tab:blue')
ax1.legend(["tub", "cstr"], loc = 1)
fig.tight_layout()
plt.savefig("q3a3")


# ## QUESTÃO 4

# In[ ]:


def f(t, X, cA0, k):
    dtdX= 1/(k*(1-X))
    return dtdX

P0 = 2.7 # atm
cA0 = P0/(R1*T)
P = np.zeros(len(X))
cA = np.zeros(len(X))
cB = np.zeros(len(X))
cC = np.zeros(len(X))
for i in range(len(X)):
    cA[i] = (cA0*(1-X[i]))
    cB[i] = (cA0*vB*X[i])
    cC[i] = (cA0*vC*X[i])
    P[i] = (cA[i] + cB[i] + cC[i])*R1*T


t = odeint(f, 0, X, args= (cA0, k)) # usamos odeint para resolver a equação diferencial em "f", para os valores de X entre 0-0,9

print("Pressão em X(0.9) = %.2f atm" % P[9])
plt.figure(0, figsize = (10,8))
plt.title("Batch", fontsize = 16)
a, = plt.plot(P, cA, label = "cA", color = "b")
b, = plt.plot(P, cB, label = "cB", color = "r")
c, = plt.plot(P, cC, label = "cC", color = "g")
plt.legend(handles = [a,b,c])
plt.xlabel("P (atm)", fontsize = 12)
plt.ylabel("Concentração (mol/L)", fontsize =12)
plt.savefig("q4a1")

plt.figure(1, figsize = (10,8))
plt.title("Batch", fontsize = 16)
a, = plt.plot(t, cA, label = "cA", color = "b")
plt.legend(handles = [a])
plt.xlabel("t (min)", fontsize = 12)
plt.ylabel("mol/L", fontsize =12)
plt.savefig("q4a2")

plt.figure(2, figsize = (10,8))
plt.title("Batch", fontsize = 16)
b, = plt.plot(t, X, label = "X", color = "b")
plt.legend(handles = [b])
plt.xlabel("t (min)", fontsize = 12)
plt.ylabel("X", fontsize =12)
plt.savefig("q4a3")


# In[ ]:


nbase = 2.5
tr = np.log(1/(1-X))*(1/k)
t = tr + 0.75*60
#calculo do Volume
nA0 = nbase*t
V7 =((nA0*R1*T)/P0)
cA0 = nA0/V7
nA = nbase*t
cA = np.zeros(len(X))
nciclos = 3600/nA
print("tempo de reação = %.0f min" % tr[9])
print("ciclos = %.1f" % nciclos[9])
print("V = %.2f L" % V7[9])
for i in range(len(nA0)):
    cA[i] = cA0[i]*(1-X[i])
    
#PLOT 
fig, ax1 = plt.subplots(figsize=(10,8))
ax1.set_title("Batch")
ax1.set_xlabel('t (min)')
ax1.set_ylabel('X', color='tab:red')
ax1.plot(t, X, color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax2 = ax1.twinx() 
ax2.set_ylabel('Concentração(mol/L)', color="tab:blue")  
ax2.plot(t, cA, color="tab:blue")
ax2.tick_params(axis='y', labelcolor="tab:blue")
fig.tight_layout()
plt.savefig("q4a5")


# ## QUESTÃO 5

# In[ ]:


V8 = V + V2
V8 = np.concatenate((V1, V8), axis = None)
V8 = np.delete(V8, 9)
X3 = np.concatenate((X1,X2), axis = None)
X3 = np.delete(X3, 9)
V9 = V + V21
V9 = np.concatenate((V1, V9), axis = None)
V9 = np.delete(V9, 9)

plt.figure(0, figsize = (10,8))
plt.title("Comparação", fontsize = 16)
a, = plt.plot(Vcstr, X, label = "cstr unico", color = "b")
b, = plt.plot(V6, X, label = "tubular", color = "r")
c, = plt.plot(V7, X, label = "batelada", color = "g")
d, = plt.plot(V8, X3, label = "2 cstr série acumulado", color = "black")
i, = plt.plot(V9, X3, label = "cstr + pfr série acumulado", color = "purple")
plt.legend(handles = [a,b,c,d,i])
plt.xlabel("Volume (L)", fontsize = 12)
plt.ylabel("X", fontsize =12)
plt.savefig("q5a6")

