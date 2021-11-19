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
    m = np.array([1., 1., 1e-5])*Ms  # Masses in Solar mass
    a = np.array([1., 1., 5.])*au   # Semi-major axis in astronomical units
    e = np.array([0., 0., 1./4.])   # Eccentricity
    psi = np.array([0., 0., 0.])*np.pi/180.    # Inclination of the orbital plane in degrees

    x1 = np.array([0., -1., 0.])*a[0]
    x2 = np.array([0., 1., 0.])*a[1]
    x3 = np.array([np.cos(psi[2]), 0., np.sin(psi[2])])*a[2]
    q = np.array([x1, x2, x3])

    v1 = np.array([np.sqrt(G*m[1]**2/((m[0]+m[1])*np.sqrt(np.sum((q[0]-q[1])**2)))), 0., 0.])
    v2 = np.array([-np.sqrt(G*m[0]**2/((m[0]+m[1])*np.sqrt(np.sum((q[0]-q[1])**2)))), 0., 0.])
    v3 = np.array([0., np.sqrt(G*(m[0]+m[1])*(2./np.sqrt(np.sum(q[2]**2))-1./a[2])), 0.])
    v = np.array([v1, v2, v3])

    #integration parameters
    duration, step = 100*yr, np.array([1./(365.25*2.), 1./(365.25*2.), 1./365.25])*yr #integration time and step in years

    integrator = "leapfrog"
    n_bodies = 2
    display = True
    savename = "{0:d}bodies_{1:s}".format(n_bodies, integrator)

    #simulation start
    bodylist = []
    for i in range(n_bodies):
        bodylist.append(Body(m[i], q[i], v[i]))
    bin_syst = System(bodylist[0:2])
    dyn_syst = System(bodylist, main=True)

    E, L = [], []
    for step0 in step:
        if integrator.lower() in ['leapfrog', 'frogleap', 'frog']:
            E0, L0, sma, ecc = leapfrog(dyn_syst, bin_syst, duration, step0, recover_param=True, display=display, savename=savename)
        elif integrator.lower() in ['hermite','herm']:
            E0, L0, sma, ecc = hermite(dyn_syst, bin_syst, duration, step0, recover_param=True, display=display, savename=savename)
        E.append(E0)
        L.append(L0)

    parameters = [duration, step, dyn_syst, integrator]
    display_parameters(E, L, sma, ecc, parameters=parameters, savename=savename)

    return 0

if __name__ == '__main__':
    sysexit(main())
