#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Hermite integrator
"""


from os import system
import numpy as np
from lib.plots import DynamicUpdate
from lib.units import *


def Drift(dyn_syst, dt):
    for body in dyn_syst.bodylist:
        body.q = body.q + dt * body.v


def Kick(dyn_syst, dt):
    for body in dyn_syst.bodylist:
        body.a = np.zeros(3,dtype=np.longdouble)
        for otherbody in dyn_syst.bodylist:
            if body != otherbody:
                rij = np.linalg.norm(body.q - otherbody.q)
                body.a = body.a - (body.q - otherbody.q) * Ga * otherbody.m / (rij ** 3)
        body.v = body.v + dt * body.a


def LP(dyn_syst, dt):
    dyn_syst.COMShift()
    Drift(dyn_syst, dt / 2)
    Kick(dyn_syst, dt)
    Drift(dyn_syst, dt / 2)
    dyn_syst.time = dyn_syst.time + dt

def leapfrog(dyn_syst, bin_syst, duration, dt, recover_param=False, display=False, savename=None, gif=False):
    if display:
        try:
            system("mkdir tmp")
        except IOError:
            system("rm tmp/*")
        d = DynamicUpdate(dyn_syst)
        d.launch(dyn_syst.blackstyle)

    N = np.ceil(duration / dt).astype(int)
    E = np.zeros(N,dtype=np.longdouble)
    L = np.zeros((N, 3),dtype=np.longdouble)
    sma = np.zeros(N,dtype=np.longdouble)
    ecc = np.zeros(N,dtype=np.longdouble)
    for j in range(N):
        LP(dyn_syst,dt)

        E[j] = dyn_syst.E
        L[j] = dyn_syst.L
        sma[j] = bin_syst.sma
        ecc[j] = bin_syst.ecc

        if display and j % 10 == 0:
            if gif:
                step = j
            else:
                step = None
            # display progression
            if len(dyn_syst.bodylist) == 1:
                d.on_running(dyn_syst, step=step, label="{0:.2f} years".format(j*dt/yr))
            else:
                d.on_running(dyn_syst, step=step, label="{0:.2f} years".format(j*dt/yr))
    if display:
        d.close()
        if gif:
            system("convert -delay 5 -loop 0 tmp/??????.png tmp/temp.gif && rm tmp/??????.png")
            system("convert tmp/temp.gif -fuzz 10% -layers Optimize plots/{0:s}_dynsyst.gif".format(savename))

    if recover_param:
        return E, L, sma, ecc