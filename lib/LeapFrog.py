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
        body.a = np.zeros(3)
        for otherbody in dyn_syst.bodylist:
            if body != otherbody:
                rij = np.linalg.norm(body.q - otherbody.q)
                body.a = body.a - (body.q - otherbody.q) * G * otherbody.m / (rij ** 3)
        body.v = body.v + dt * body.a


def LP(dyn_syst, dt):
    dyn_syst.COMShift()
    Drift(dyn_syst, dt / 2)
    Kick(dyn_syst, dt)
    Drift(dyn_syst, dt / 2)
    dyn_syst.time = dyn_syst.time + dt
    for body in dyn_syst.bodylist:
        body.p = body.v * body.m


def leapfrog(dyn_syst, bin_syst, duration, dt, recover_param=False, display=False, savename=None):
    if display:
        try:
            system("mkdir tmp")
        except IOError:
            system("rm tmp/*")
        d = DynamicUpdate(dyn_syst)
        d.launch(dyn_syst.blackstyle)

    N = np.ceil(duration / dt).astype(int)
    E = np.zeros(N)
    L = np.zeros((N, 3))
    sma = np.zeros(N)
    ecc = np.zeros(N)
    for j in range(N):
        LP(dyn_syst,dt)

        E[j] = dyn_syst.E
        L[j] = dyn_syst.L
        sma[j] = bin_syst.sma
        ecc[j] = bin_syst.ecc

        if display and j % 5 == 0:
            # display progression
            if len(dyn_syst.bodylist) == 1:
                d.on_running(dyn_syst, step=j, label="step {0:d}/{1:d}".format(j, N))
            else:
                d.on_running(dyn_syst, step=j, label="step {0:d}/{1:d}".format(j, N))
    if display:
        d.close()
        if not savename is None:
            system("convert -delay 5 -loop 0 tmp/??????.png tmp/temp.gif && rm tmp/??????.png")
            system("convert tmp/temp.gif -fuzz 10% -layers Optimize plots/{0:s}_dynsyst.gif".format(savename))

    if recover_param:
        return E, L, sma, ecc