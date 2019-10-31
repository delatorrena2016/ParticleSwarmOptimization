# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal.
"""
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
class Particula:
    def __init__(self,pid,X,V,Xb,fun):
        self.pid, self.X, self.V, self.Xb = pid, X, V, Xb
        self.fitness = fun(self.X)
    def mostrarDatos(self):
        print(f"id:{self.pid} fit:{self.fitness}")
#%%  function sphere      
def sphere(X):
    Z = []
    for i in range(X.shape[0]):
        Z.append(sum(val**2 for val in X[i,:]))
    return np.array(Z)
#%% Parameters
param = {'N':10000, 'D':2, 'phi1':2, 'phi2':2, 'omega':0.7, 'it':15, 'fun':sphere}
print('PSO parameters:', param)
#%% Test object creation
X = np.random.rand(param.get('N'),param.get('D')) *1000
X = np.random.normal(5,0.01,size=(10000,2))

X = 20 * np.random.rand(200,2) - 10
#[0,1] -> [-10,10] y = mx + b
Z = sphere(X)
plt.scatter(X[:,0], X[:,1], c='b',s=20)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X[:,0], X[:,1], Z, color='b')
plt.show()
#X = np.array([5,4])
p1 = Particula(1,X,X,X,param.get('fun'))
p1.mostrarDatos()
