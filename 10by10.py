from __future__ import division
import numpy as np
import random as rd

# initialization

N = 10 # number of nodes 

A = np.zeros((N,N)) # connectivity matrix

for i in range(N):
    for j in range(N):
        if j == i:
            A[i,j] = 0
        else:
            A[i,j] = 1

dt = 0.01 # time step

sig = 0.25 # level of noise

sd = sig / np.sqrt(dt) # SD of the discrete noise

T = 100000 # number of time steps

x = np.zeros((N,T)) # x_i(t)
v = np.zeros((N,T)) # v_i(t)

for i in range(N):
    x[i,0] = rd.gauss(0,sig) # initial setup

# simulating

y = np.zeros((N, T))

for i in range(N):
    for j in range(N):
        y[i,0] = y[i,0] + A[i,j] * (x[j,0] - x[i,0])
    v[i,0] = y[i,0] + sd * rd.gauss(0, sig)

for t in range(T-1):
    for i in range(N):
        x[i,t+1] = x[i,t] + dt * v[i,t]
    for i in range(N):
        for j in range(N):
            y[i,t+1] = y[i,t+1] + A[i,j] * (x[j,t+1] - x[i,t+1])
        v[i,t+1] = y[i,t+1] + sd * rd.gauss(0, sig)

# reconstructing

M = np.zeros((N-1,N-1))
u = np.zeros(N-1)

for i in range(N-1):
    for j in range(N-1):
        for t in range(T):
            M[i,j] = M[i,j] + (x[i+1,t] - x[0,t]) * (x[j+1,t] - x[0,t])
            u[i] = u[i] + v[0,t] * (x[i+1,t] - x[0,t])

# det = np.linalg.det(xx)
# print det

B = np.linalg.solve(M, u)
print B
