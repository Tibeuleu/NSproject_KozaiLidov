#!/usr/bin/python
# -*- coding:utf-8 -*-
from sys import exit as sysexit
import numpy as np
import matplotlib.pyplot as plt
from lib.integrator import frogleap
from lib.objects import Body, System

globals()['G'] = 6.67e-11 #Gravitational constant in SI units
globals()['Ms'] = 2e30 #Solar mass in kg
globals()['au'] = 1.5e11 #Astronomical unit in m

def main():
    #initialisation
    m = np.array([1, 1, 0.1])*Ms  # Masses in Solar mass
    mu = m[0]*m[1]/(m[0]+m[1])
    a = np.array([1., 1., 5.])*au   # Semi-major axis in astronomical units
    psi = np.array([0., 0., 80.])*np.pi/180.    # Inclination of the orbital plane in degrees

    x1 = np.array([-1., 0., 0.])*a[0]
    x2 = np.array([1., 0., 0.])*a[1]
    x3 = np.array([np.cos(psi[2]), 0., np.sin(psi[2])])*a[2]
    q = np.array([x1, x2, x3])

    v1 = np.array([0., -np.sqrt(G*mu/np.sqrt(np.sum(x1**2))), 0])
    v2 = np.array([0., np.sqrt(G*mu/np.sqrt(np.sum(x2**2))), 0.])
    v3 = np.array([0., np.sqrt(G*(m[0]+m[1])*(2./np.sqrt(np.sum(x3**2))-1./a[2])), 0.])
    v = np.array([v1, v2, v3])

    bodylist = []
    for i in range(m.shape[0]):
        bodylist.append(Body(m[i], q[i], v[i]))
    dyn_syst = System(bodylist)
    dyn_syst.COMShift()

    duration, step = 0.5*3e7, 1e1
    E, L = frogleap(duration, step, dyn_syst, recover_param=True)#, display=True)
    
    fig1 = plt.figure(figsize=(30,15))
    ax1 = fig1.add_subplot(111)
    ax1.plot(np.arange(E.shape[0])/duration, E, label=r"$E_m$")
    ax1.legend()
    fig1.savefig("plots/Em.png",bbox_inches="tight")
    fig2 = plt.figure(figsize=(30,15))
    ax2 = fig2.add_subplot(111)
    ax2.plot(np.arange(L.shape[0])/duration, np.sum(L**2,axis=1), label=r"$L^2$")
    ax2.legend()
    fig2.savefig("plots/L2.png",bbox_inches="tight")
    plt.show(block=True)
    
    return 0

if __name__ == '__main__':
    sysexit(main())
