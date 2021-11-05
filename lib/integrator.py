#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Implementation of the various integrators for numerical integration.

Comes from the assumption that the problem is analytically defined in position-momentum (q-p) space for a given hamiltonian H.
"""
from os import system
import time
import numpy as np
from lib.plots import DynamicUpdate

globals()['G'] = 6.67e-11 #Gravitational constant in SI units
globals()['Ms'] = 2e30 #Solar mass in kg
globals()['au'] = 1.5e11 #Astronomical unit in m

def dv_dt(m_array, q_array):
    """
    Time derivative of the velocity, given by the position derivative of the Hamiltonian.
    dv/dt = -1/m*dH/dq
    """
    dv_array = np.zeros(q_array.shape)
    for i in range(q_array.shape[0]):
        q_j = np.delete(q_array, i, 0)
        m_j = np.delete(m_array, i, 0)
        dv_array[i] = -G*np.sum((m_j*(q_j-q_array[i])).T/np.sqrt(np.sum((q_j-q_array[i])**2, axis=1))**3, axis=1).T
    dv_array[np.isnan(dv_array)] = 0.
    return dv_array

def frogleap(duration, step, dyn_syst, recover_param=False, display=False):
    """
    Leapfrog integrator for first order partial differential equations.
    iteration : half-step drift -> full-step kick -> half-step drift
    """
    N = np.ceil(duration/step).astype(int)
    q_array = dyn_syst.get_positions()
    v_array = dyn_syst.get_velocities()
    masses = dyn_syst.get_masses()
    m_array = np.ones(q_array.shape)
    for i in range(q_array.shape[0]):
        m_array[i,:] = masses[i]
    
    E = np.zeros(N)
    L = np.zeros((N,3))
    
    if display:
        try:
            system("mkdir tmp")
        except IOError:
            system("rm tmp/*")
        d = DynamicUpdate()
        d.on_launch()
    for j in range(N):
        # half-step drift
        q_array, v_array = q_array + step/2*v_array , v_array
        # full-step kick
        q_array, v_array = q_array , v_array - step*dv_dt(m_array, q_array)
        # half-step drift
        q_array, v_array = q_array + step/2*v_array , v_array
        
        for i, body in enumerate(dyn_syst.bodylist):
            body.q = q_array[i]
            body.v = v_array[i]
            body.p = body.v*body.m
        dyn_syst.COMShift()
        
        E[j] = dyn_syst.Eval()
        L[j] = dyn_syst.Lval()

        if display:
            # display progression
            if len(dyn_syst.bodylist) == 1:
                d.on_running(q_array[0], q_array[1], q_array[2], step=j, label="step {0:d}/{1:d}".format(j,N))
            else:
                d.on_running(q_array[:,0], q_array[:,1], q_array[:,2], step=j, label="step {0:d}/{1:d}".format(j,N))
            time.sleep(1e-5)
    if display:
        system("convert -delay 5 -loop 0 tmp/?????.png tmp/temp.gif && rm tmp/?????.png")
        system("convert tmp/temp.gif -fuzz 30% -layers Optimize plots/dynsyst.gif && rm tmp/temp.gif")
        
    if recover_param:
        return E, L
