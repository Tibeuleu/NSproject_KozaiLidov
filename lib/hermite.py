#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Hermite integrator
"""


from os import system
import numpy as np
from lib.plots import DynamicUpdate
from lib.units import *


def Update_a(dyn_syst):  # update acceleration of bodies in system
    for body in dyn_syst.bodylist:
        body.a = np.zeros(3)
        for otherbody in dyn_syst.bodylist:
            if body != otherbody:
                rij = np.linalg.norm(body.q - otherbody.q)
                body.a = body.a - (body.q - otherbody.q) * G * otherbody.m / (rij ** 3)


def Update_j(dyn_syst):  # update jerk of bodies in system
    for body in dyn_syst.bodylist:
        body.j = np.zeros(3)
        for otherbody in dyn_syst.bodylist:
            if body != otherbody:
                rij = np.linalg.norm(body.q - otherbody.q)
                deltav = (body.v - otherbody.v)
                deltar = (body.q - otherbody.q)
                vr = deltav + 3. * deltar * np.inner(deltav, deltar) / (rij ** 2)
                body.j = body.j - G * otherbody.m / (rij ** 3) * vr


def Predict(dyn_syst, dt):  # update predicted position and velocities of bodies in system
    for body in dyn_syst.bodylist:
        body.qp = body.q + dt * body.v + ((dt ** 2) * body.a / 2.) + ((dt ** 3) * body.j / 6.)
        body.vp = body.v + dt * body.a + ((dt ** 2) * body.j / 2.)


def Update_ap(dyn_syst):  # update acceleration of bodies in system
    for body in dyn_syst.bodylist:
        body.ap = np.zeros(3)
        for otherbody in dyn_syst.bodylist:
            if body != otherbody:
                rij = np.linalg.norm(body.qp - otherbody.qp)
                body.ap = body.ap - (body.qp - otherbody.qp) * G * otherbody.m / (rij ** 3)


def Update_jp(dyn_syst):  # update jerk of bodies in system
    for body in dyn_syst.bodylist:
        body.jp = np.zeros(3)
        for otherbody in dyn_syst.bodylist:
            if body != otherbody:
                rij = np.linalg.norm(body.qp - otherbody.qp)
                deltav = (body.vp - otherbody.vp)
                deltar = (body.qp - otherbody.qp)
                vr = deltav + 3. * deltar * np.inner(deltav, deltar) / (rij ** 2)
                body.jp = body.jp - G * otherbody.m / (rij ** 3) * vr


def Correct(dyn_syst, dt):  # correct position and velocities of bodies in system
    for body in dyn_syst.bodylist:
        a2 = (6. * (body.a - body.ap) + dt * (4 * body.j + 2 * body.jp)) / (dt ** 2)
        a3 = (12. * (body.a - body.ap) + dt * 6. * (body.j + body.jp)) / (dt ** 3)

        body.q = body.qp + ((dt ** 4) * a2 / 24.) + ((dt ** 5) * a3 / 120.)
        body.v = body.vp + ((dt ** 3) * a2 / 6.) + ((dt ** 4) * a3 / 24.)


def HPC(dyn_syst, dt):  # update position and velocities of bodies in system with hermite predictor corrector
    COMShift(dyn_syst)
    Update_a(dyn_syst)
    Update_j(dyn_syst)
    Predict(dyn_syst, dt)
    Update_ap(dyn_syst)
    Update_jp(dyn_syst)
    Correct(dyn_syst, dt)
    dyn_syst.time = dyn_syst.time + dt


def hermite(dyn_syst, duration, dt, recover_param=False, display=False, savename=None):
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
        HPC(dyn_syst, dt)

        E[j] = dyn_syst.E
        L[j] = dyn_syst.L
        sma[j] = bin_syst.sma
        ecc[j] = bin_syst.ecc

        if display and j % 100 == 0:
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