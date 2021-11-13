#!/usr/bin/python
# -*- coding:utf-8 -*-
from sys import exit as sysexit
import numpy as np
import matplotlib.pyplot as plt
from lib.objects import Body, System
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

    bodylist = []
    for i in range(3):
        bodylist.append(Body(m[i], q[i], v[i]))
    dyn_syst = System(bodylist)
    dyn_syst.COMShift()

    duration, step1, step2 = 100*yr, 1e4, 1e5
    E1, L1 = dyn_syst.leapfrog(duration, step1, recover_param=True, display=True)
    E2, L2 = dyn_syst.leapfrog(duration, step2, recover_param=True)#, display=True)
    #E1, L1 = dyn_syst.hermite(duration, step1, recover_param=True)#, display=True)
    #E2, L2 = dyn_syst.hermite(duration, step2, recover_param=True)#, display=True)
    parameters = [duration, [step1, step2], dyn_syst, "leapfrog"]
    display_parameters([E1, E2], [L1, L2], parameters=parameters, savename="3bodies_mass_leapfrog")

    return 0

if __name__ == '__main__':
    sysexit(main())
