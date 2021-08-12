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


#S - SANGUE; E - ESTOMAGO

def dif (c, t):
    k1 = 0.15
    k2 = 0.6
    k4 = 0.2
    cE = c[0]
    cS = c[1]
    if cS < 0:
        k3 = k1*cE - k4*cS
    else:
        k3 = 0.1
    dcedt = (-k1-k2)*cE
    dcsdt = k1*cE - k3 - k4*cS
    return [dcedt, dcsdt]
        

cE0 = 250/40
x0 = [cE0, 0]
t = np.linspace(0, 24, 25)
c = odeint(dif, x0,t)
cE = c[:,0]
cS = c[:,1]

plt.figure(0, figsize = (10,8))
plt.title("Tarzlon no sangue", fontsize = 16)
a, = plt.plot(t, cS, label = "cS", color = "b")
plt.legend(handles = [a])
plt.xlabel("t (hora)", fontsize = 12)
plt.ylabel("Concentração (mol/L)", fontsize =12)
plt.savefig("q31")


# In[3]:


n = 2
cE0 = n*250/40
x0 = [cE0, 0]
c = odeint(dif, x0,t)
cE = c[:,0]
cS = c[:,1]
print(cE[4])
plt.figure(0, figsize = (10,8))
plt.title("Tarzlon no sangue", fontsize = 16)
a, = plt.plot(t, cS, label = "cS", color = "b")
plt.legend(handles = [a])
plt.xlabel("t (hora)", fontsize = 12)
plt.ylabel("Concentração (mol/L)", fontsize =12)
plt.savefig("q32")


# In[4]:


print(cE[4], cS[4])
cE0 = 250/40
x0 = [cE[4] + cE0, cS[4]]
t = np.linspace(4, 24, 21)
c = odeint(dif, x0,t)
cE = c[:,0]
cS = c[:,1]
print(cS)
plt.figure(0, figsize = (10,8))
plt.title("Tarzlon no sangue", fontsize = 16)
a, = plt.plot(t, cS, label = "cS", color = "b")
plt.legend(handles = [a])
plt.xlabel("t (hora)", fontsize = 12)
plt.ylabel("Concentração (mol/L)", fontsize =12)
plt.savefig("q33")


# In[5]:


print(cE[4], cS[4])
cE0 = 250/40
x0 = [cE[4] + cE0, cS[4]]
t = np.linspace(8, 24, 21)
c = odeint(dif, x0,t)
cE = c[:,0]
cS = c[:,1]
print(cS)
plt.figure(0, figsize = (10,8))
plt.title("Tarzlon no sangue", fontsize = 16)
a, = plt.plot(t, cS, label = "cS", color = "b")
plt.legend(handles = [a])
plt.xlabel("t (hora)", fontsize = 12)
plt.ylabel("Concentração (mol/L)", fontsize =12)
plt.savefig("q34")


# In[6]:


def dif1 (x, t):
    k1 = 0.15
    k2 = 0.6
    k4 = 0.2
    cE = x[0]
    cS = x[1]
    if cS < 0:
        k3 = k1*cE - k4*cS
    else:
        k3 = 0.1
    dcedt = (-k1-k2)*cE
    dcsdt = k1*cE - k3 - k4*cS
    return [dcedt, dcsdt]
        

cE0 = 250/40
x0 = [2*cE0, 0]
t = np.linspace(0, 4, 10)
c = odeint(dif1, x0,t)
cE = c[:,0]
cS = c[:,1]
foo1 = []
foo2 = []
foo2.insert(0, cS)
foo1.insert(0,cE)
var1 = cE[9]
for i in range(1, 21):
    t = np.linspace(4*i, 4*(i+1), 10)
    x0 = [var1 + cE0, cS[9]]
    c  = odeint(dif1, x0,t)
    cE = c[:,0]
    cS = c[:,1]
    foo1.insert(i, cE)
    foo2.insert(i, cS)
    var1 = cE[9]

cS = np.array(foo2)
cE = np.array(foo1)
cS = cS.flatten()
cE = cE.flatten()
for i in range(1, 21):
    cS = np.delete(cS, i*10 -i+1)
    cE = np.delete(cE, i*10 - i + 1)

t = np.linspace(0, 48, 190) 
plt.figure(0, figsize = (10,8))
plt.title("Tarzlon no sangue", fontsize = 16)
a, = plt.plot(t, cS, label = "cS", color = "b")
plt.legend(handles = [a])
plt.xlabel("t (hora)", fontsize = 12)
plt.ylabel("Concentração (mol/L)", fontsize =12)
plt.savefig("q35") 

plt.figure(1, figsize = (10,8))
plt.title("Tarzlon no estômago", fontsize = 16)
a, = plt.plot(t, cE, label = "cE", color = "r")
plt.legend(handles = [a])
plt.xlabel("t (hora)", fontsize = 12)
plt.ylabel("Concentração (mol/L)", fontsize =12)
plt.savefig("q36") 


# In[ ]:




