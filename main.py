#!/usr/bin/python
# -*- coding:utf-8 -*-
from sys import exit as sysexit
import numpy as np
from lib.integrator import frogleap
from lib.objects import Body, System

def main():
    #initialisation
    m = np.array([1e10, 1, 0])

    x1 = np.array([0, 0, 0])
    x2 = np.array([1, 0, 0])
    x3 = np.array([10, 0, 0])
    q = np.array([x1, x2, x3])

    v1 = np.array([0, 0, 0])
    v2 = np.array([0, 0, 0])
    v3 = np.array([0, 0, 0])
    v = np.array([v1, v2, v3])

    bodylist = []
    for i in range(3):
        bodylist.append(Body(m[i], q[i], v[i]))
    dyn_syst = System(bodylist)

    new_dyn_syst = frogleap(10, 0.01, dyn_syst, display=True)
    return 0

if __name__ == '__main__':
    sysexit(main())
