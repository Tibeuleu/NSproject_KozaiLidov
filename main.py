#!/usr/bin/python
# -*- coding:utf-8 -*-
from sys import exit as sysexit
import numpy as np
import matplotlib.pyplot as plt
from lib.objects import Body, System
from lib.LeapFrog import leapfrog
from lib.hermite import hermite
from lib.plots import display_parameters
from lib.units import *

def main():
    #initialisation
    m = np.array([1., 1., 1e-1],dtype=np.longdouble)*Ms/Ms  # Masses in Solar mass
    a = np.array([1., 1., 5.],dtype=np.longdouble)*au/au   # Semi-major axis in astronomical units
    e = np.array([0., 0., 0.],dtype=np.longdouble)   # Eccentricity
    psi = np.array([0., 0., 0.],dtype=np.longdouble)*np.pi/180.    # Inclination of the orbital plane in degrees

    x1 = np.array([0., -1., 0.],dtype=np.longdouble)*a[0]*(1.+e[0])
    x2 = np.array([0., 1., 0.],dtype=np.longdouble)*a[1]*(1.+e[1])
    x3 = np.array([np.cos(psi[2]), 0., np.sin(psi[2])],dtype=np.longdouble)*a[2]*(1.+e[2])
    q = np.array([x1, x2, x3],dtype=np.longdouble)

    v1 = np.array([np.sqrt(Ga*m[0]*m[1]/((m[0]+m[1])*np.sqrt(np.sum((q[0]-q[1])**2)))), 0., 0.],dtype=np.longdouble)
    v2 = np.array([-np.sqrt(Ga*m[0]*m[1]/((m[0]+m[1])*np.sqrt(np.sum((q[0]-q[1])**2)))), 0., 0.],dtype=np.longdouble)
    v3 = np.array([0., np.sqrt(Ga*(m[0]+m[1])*(2./np.sqrt(np.sum(q[2]**2))-1./a[2])), 0.],dtype=np.longdouble)
    v = np.array([v1, v2, v3],dtype=np.longdouble)

    #integration parameters
    duration, step = 100*yr, np.array([60.],dtype=np.longdouble) #integration time and step in seconds
    step = np.sort(step)[::-1]
    integrator = "leapfrog"
    n_bodies = 3
    display = False
    gif = False
    savename = "{0:d}bodies_{1:s}".format(n_bodies, integrator)

    #simulation start
    E, L = [], []
    for i,step0 in enumerate(step):
        bodylist = []
        for j in range(n_bodies):
            bodylist.append(Body(m[j], q[j], v[j]))
        bin_syst = System(bodylist[0:2])
        dyn_syst = System(bodylist, main=True)

        if i != 0:
            display = False
        if integrator.lower() in ['leapfrog', 'frogleap', 'frog']:
            E0, L0, sma, ecc = leapfrog(dyn_syst, bin_syst, duration, step0, recover_param=True, display=display, savename=savename, gif=gif)
        elif integrator.lower() in ['hermite','herm']:
            E0, L0, sma, ecc = hermite(dyn_syst, bin_syst, duration, step0, recover_param=True, display=display, savename=savename, gif=gif)
        E.append(E0)
        L.append(L0)

    parameters = [duration, step, dyn_syst, integrator]
    display_parameters(E, L, sma, ecc, parameters=parameters, savename=savename)
    return 0

if __name__ == '__main__':
    sysexit(main())
