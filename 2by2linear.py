from __future__ import division
import numpy as np
import random as rd
import pylab as plt


# initialization

alpha = 1 # connectivity matrix

beta = 2 # connectivity

dt = 0.01 # time step

sig = 0.5 # level of noise

sd = sig/np.sqrt(dt) # SD of the discrete noise

T = 10 # total time step

M = 1000 # number of time steps

TR = 1000 # number of trials

N = 2 # number of nodes

x1 = rd.gauss(0, sig) # initial setup
x2 = rd.gauss(0, sig)
v1 = 0
v2 = 0

# simulating the dynamics to obtain data and performing least square fit

xxsum = 0
v1xsum = 0
v2xsum = 0
a = np.zeros(TR)
b = np.zeros(TR)

for m in range(TR):
    for n in range(M):
        v1 = alpha * (x2 - x1) + sd * rd.gauss(0, sig)
        v2 = beta  * (x1 - x2) + sd * rd.gauss(0, sig)
        xxsum = xxsum + (x1 - x2) * (x1 - x2)
        v1xsum = v1xsum + v1 * (x2 - x1)
        v2xsum = v2xsum + v2 * (x1 - x2)
        x1 = x1 + v1 * dt
        x2 = x2 + v2 * dt

    a[m] = v1xsum / xxsum # estimate of alpha
    b[m] = v2xsum / xxsum # estimate of beta

print np.average(a)
print np.average(b)
print np.std(a)
print np.std(b)

# plot

bins1 = np.arange(0, 2, 0.01)
plt.xlim([0, 3])
plt.hist(a, bins=bins1)
plt.ylabel('count')


bins2 = np.arange(1, 3, 0.01)
plt.hist(b, bins=bins2)
plt.title('Distribution of reconstructed connectivity')
plt.ylabel('count')

plt.show()
