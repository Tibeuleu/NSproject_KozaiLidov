#!/usr/bin/python
# -*- coding:utf-8 -*-
from sys import exit as sysexit
import numpy as np
import matplotlib.pyplot as plt
from lib.integrator import frogleap
from lib.objects import Body, System

def main():
    #initialisation
    m = np.array([1, 1, 1e-5])

    x1 = np.array([-1, 0, 0])
    x2 = np.array([1, 0, 0])
    x3 = np.array([100, 0, 0])
    q = np.array([x1, x2, x3])

    v1 = np.array([0, -0.35, 0])
    v2 = np.array([0, 0.35, 0])
    v3 = np.array([0, 0, 0])
    v = np.array([v1, v2, v3])

    bodylist = []
    for i in range(2):
        bodylist.append(Body(m[i], q[i], v[i]))
    dyn_syst = System(bodylist)
    dyn_syst.COMShift()

    duration, step = 100, 0.01
    E, L = frogleap(duration, step, dyn_syst, recover_param=True, display=True)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(np.arange(E.shape[0])/duration, E, label=r"$E_m$")
    ax1.legend()
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(np.arange(L.shape[0])/duration, np.sum(L**2,axis=1), label=r"$L^2$")
    ax2.legend()
    plt.show(block=True)
    return 0

if __name__ == '__main__':
    sysexit(main())
