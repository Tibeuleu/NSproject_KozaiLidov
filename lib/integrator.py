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

globals()["G"] = 1. #Gravitationnal constant

def dp_dt(m_array, q_array):
    """
    Time derivative of the momentum, given by the position derivative of the Hamiltonian.
    dp/dt = -dH/dq
    """
    dp_array = np.zeros(q_array.shape)
    for i in range(q_array.shape[0]):
        q_j = np.delete(q_array, i, 0)
        m_j = np.delete(m_array, i, 0)#.reshape((q_j.shape[0],1))
        dp_array[i] = -G*m_array[i]*np.sum(m_j/np.sum(np.sqrt(np.sum((q_j-q_array[i])**2, axis=0)))**3*(q_j-q_array[i]), axis=0)
    dp_array[np.isnan(dp_array)] = 0.
    return dp_array

def frogleap(duration, step, dyn_syst, recover_param=False, display=False):
    """
    Leapfrog integrator for first order partial differential equations.
    iteration : half-step drift -> full-step kick -> half-step drift
    """
    N = np.ceil(duration/step).astype(int)
    q_array = dyn_syst.get_positions()
    p_array = dyn_syst.get_momenta()
    masses = dyn_syst.get_masses()
    m_array = np.ones(p_array.shape)
    for i in range(p_array.shape[0]):
        m_array[i,:] = masses[i]
    
    E = np.zeros(N)
    L = np.zeros((N,3))
    
    if display:
        try:
            system("mkdir tmp")
        except IOError:
            system("rm tmp/*")
        d = DynamicUpdate()
        d.min_x, d.max_x = -1.5*np.abs(q_array).max(), +1.5*np.abs(q_array).max()
        d.on_launch()
    for j in range(N):
        # half-step drift
        q_array, p_array = q_array + step/2*p_array/m_array , p_array
        # full-step kick
        q_array, p_array = q_array , p_array - step*dp_dt(m_array, q_array)
        # half-step drift
        q_array, p_array = q_array + step/2*p_array/m_array , p_array
        
        for i, body in enumerate(dyn_syst.bodylist):
            body.q = q_array[i]
            body.p = p_array[i]
            if body.m != 0.:
                body.v = body.p/body.m
        dyn_syst.COMShift()
        
        E[j] = dyn_syst.Eval()
        L[j] = dyn_syst.Lval()

        if display:
            # In center of mass frame
            q_cm = np.array([0,0,0])#np.sum(m_array*q_array, axis=0)/masses.sum()
            # display progression
            d.on_running(q_array[:,0]-q_cm[0], q_array[:,1]-q_cm[1], q_array[:,2]-q_cm[2], step=j, label="step {0:d}/{1:d}".format(j,N))
            time.sleep(1e-4)
    if display:
        system("convert -delay 5 -loop 0 tmp/?????.png tmp/temp.gif && rm tmp/?????.png")
        system("convert tmp/temp.gif -fuzz 30% -layers Optimize plots/dynsyst.gif && rm tmp/temp.gif")
        
    if recover_param:
        return E, L
